// SPDX-License-Identifier: GPL-3.0-or-later
// Persistent native sessions (T19): one worker runs the native pipeline; the host polls
// progress, answers typed solve requests and collects the result. Workers never hold
// Python objects or call back into the host: every exchange goes through the queues.
#include "cleanseams/session.hpp"
#include "pipeline_services.hpp"
#include "cleanseams/layout_postprocess.hpp"
#include "cleanseams/cache.hpp"
#include "cleanseams/pipeline.hpp"
#include "cleanseams/errors.hpp"
#include "cleanseams/report.hpp"
#include <atomic>
#include <cmath>
#include <condition_variable>
#include <deque>
#include <map>
#include <set>
#include <mutex>
#include <thread>
namespace cleanseams {
namespace {
std::atomic<std::uint64_t> next_session_id{1};
using Json=nlohmann::json;
}
namespace {
struct SessionState;
// The host solver service seen from the worker: builds the source BlenderBackend's
// isolated chart mesh, posts one typed request and waits for its answer. The planar
// n-gon maps and the conformal retry policy stay native (source adapter.BlenderBackend).
class HostUnwrapBackend final:public UnwrapBackend {
public:
    HostUnwrapBackend(SessionState& session,std::shared_ptr<const MeshSnapshot> mesh,std::string method,std::uint64_t iterations):
        session_(session),mesh_(std::move(mesh)),method_(std::move(method)),iterations_(iterations) {}
    std::vector<Vec2> unwrap(std::span<const std::uint8_t> cut,std::optional<std::span<const MeshIndex>> charts,const PinMap& pins,
                             std::optional<std::uint64_t> iterations,std::optional<std::span<const Vec2>>) override;
    std::uint64_t iterations() const override {return iterations_;}
    // Blender stores UVs in float32: the map it answers is the float32 layer read as float64.
    bool float32_storage() const override {return true;}
    // Source BlenderBackend counters, reported as unwrap_calls and solver_paths.
    std::uint64_t calls=0,planar_maps=0,conformal_fallbacks=0;
    std::string primary; // method the host reported for the primary solves
private:
    SessionState& session_;
    std::shared_ptr<const MeshSnapshot> mesh_;
    std::string method_;
    std::uint64_t iterations_;
    // Scratch-mesh solve; returns the map and the method the host actually used.
    std::pair<std::vector<Vec2>,std::string> solve(const std::shared_ptr<const CutSurface>&,std::span<const MeshIndex> charts,const PinMap& pins,
                                                   const std::string& method,std::uint64_t iterations);
};
struct SessionState {
    SessionConfig config;
    std::uint64_t id=next_session_id++;
    mutable std::mutex mutex;
    std::mutex join_mutex; // close(): one joiner at a time
    std::condition_variable changed;
    enum class State {idle,running,done,closed} state=State::idle;
    std::deque<ProgressEvent> progress;
    std::optional<HostRequest> pending;
    bool pending_taken=false;
    std::optional<HostResponse> response;
    RequestId next_request=1;
    std::optional<RunResult> result;
    Control control;
    std::thread worker;
    Revision geometry_revision=0;
    std::shared_ptr<const MeshSnapshot> mesh;
    std::chrono::steady_clock::time_point started;
    // Worker side: publish a request and wait for the matching response or cancellation.
    HostResponse exchange(HostRequest request,const std::function<bool()>& stopped={}) {
        std::unique_lock lock(mutex);
        control.check();
        auto expired=[&] {
            if(!stopped || !stopped())return false;
            pending.reset();response.reset();changed.notify_all();return true;
        };
        auto timeout=[] {HostResponse r;r.error=ErrorInfo{"host_failure","host_request","optional host phase expired"};return r;};
        if(expired())return timeout();
        request.identity.session_id=id;request.identity.request_id=next_request++;request.identity.geometry_revision=geometry_revision;
        pending=std::move(request);pending_taken=false;response.reset();changed.notify_all();
        auto answered=[&] {return response.has_value() || control.cancel.requested() || state==State::closed;};
        if(stopped) {
            while(!answered()) {
                if(expired())return timeout();
                changed.wait_for(lock,std::chrono::milliseconds(10),answered);
            }
        } else changed.wait(lock,answered);
        // Cancellation wins even when a response was queued before this worker woke.
        if(control.cancel.requested() || state==State::closed) {
            pending.reset();response.reset();throw ContractError(ErrorCode::cancelled,"session","cancelled while waiting for the host");
        }
        if(expired())return timeout();
        auto out=std::move(*response);response.reset();pending.reset();return out;
    }
    void publish(const Json& event) {
        control.check(); // cancellation checkpoint (source Control.check in emit)
        ProgressEvent e;e.stage=event.value("stage","");
        if(event.contains("done") && event.at("done").is_number())e.done=event.at("done").get<std::uint64_t>();
        if(event.contains("total") && event.at("total").is_number())e.total=event.at("total").get<std::uint64_t>();
        e.elapsed_s=event.value("elapsed",0.0);if(event.contains("eta"))e.eta_s=event.at("eta").get<double>();
        e.counters=event;e.counters.erase("elapsed");e.counters.erase("eta");
        std::lock_guard lock(mutex);progress.push_back(std::move(e));
    }
};
std::pair<std::vector<Vec2>,std::string> HostUnwrapBackend::solve(const std::shared_ptr<const CutSurface>& cs,std::span<const MeshIndex> charts,
                                                                  const PinMap& pins,const std::string& method,std::uint64_t iterations) {
    // Host requests leave from the pipeline thread only, never from a native worker (T24).
    if(in_worker_task())throw ContractError(ErrorCode::host_failure,"host.solve","host solves are requested from the pipeline thread only");
    // The source's scratch mesh: one vertex per cut sector of the requested charts, their
    // faces in face order, planar n-gon maps as the start and pins on scratch corners.
    const auto& m=*mesh_;const auto& ct=cs->tables();const auto& in=m.input();
    std::set<MeshIndex> wanted(charts.begin(),charts.end());std::vector<MeshIndex> faces;
    for(MeshIndex f=0;f+1<in.face_offsets.size();++f)if(wanted.contains(ct.chart_of_face[f]))faces.push_back(f);
    std::vector<Vec2> out(ct.sector_of_corner.size(),Vec2{0,0});
    if(faces.empty())return {out,method};
    ++calls;
    HostRequest request;request.kind="solve";HostSolveRequest solve;solve.method=method;solve.iterations=iterations;
    std::map<MeshIndex,TransferIndex> local;
    for(auto f:faces)for(auto c:m.face_corners(f))local.emplace(ct.sector_of_corner[c],0);
    TransferIndex next=0;for(auto& [sector,index]:local) {index=next++;solve.mesh.positions.push_back(in.positions[ct.sector_vert[sector]]);}
    for(auto f:faces) {
        for(auto c:m.face_corners(f)) {solve.mesh.corner_vertices.push_back(local.at(ct.sector_of_corner[c]));solve.scratch_to_original_corner.push_back(c);}
        solve.mesh.face_offsets.push_back(solve.mesh.corner_vertices.size());
    }
    const auto geometry=quality_geometry(cs);
    std::vector<Vec2> initial(solve.scratch_to_original_corner.size(),Vec2{0,0});
    for(auto chart:charts) {
        ChartSolveRequest planar;planar.geometry=geometry;planar.chart=chart;planar.validate=false;
        if(auto map=planar_ngon_map(planar)) {
            std::map<MeshIndex,Vec2> by_corner;for(std::size_t i=0;i<map->corner_ids.size();++i)by_corner[map->corner_ids[i]]=map->corner_uv[i];
            for(std::size_t i=0;i<initial.size();++i)
                if(auto hit=by_corner.find(static_cast<MeshIndex>(solve.scratch_to_original_corner[i]));hit!=by_corner.end())initial[i]=hit->second;
        }
    }
    solve.initial_uv=std::move(initial);
    for(std::size_t i=0;i<solve.scratch_to_original_corner.size();++i)
        if(auto hit=pins.find(static_cast<MeshIndex>(static_cast<std::uint64_t>(solve.scratch_to_original_corner[i])));hit!=pins.end())solve.pins.push_back({static_cast<TransferIndex>(i),hit->second});
    const auto corners=solve.scratch_to_original_corner.size();
    request.solve=std::move(solve);
    auto answer=session_.exchange(std::move(request));
    if(!answer.success)throw ContractError(ErrorCode::host_failure,"host.solve",answer.error?answer.error->message:"host solve failed");
    if(!answer.uv || answer.uv->size()!=corners)throw ContractError(ErrorCode::host_failure,"host.solve","host returned a map of the wrong size");
    // Map scratch corners back to the original corners.
    std::size_t i=0;for(auto f:faces)for(auto c:m.face_corners(f))out[c]=(*answer.uv)[i++];
    return {out,answer.method.empty()?method:answer.method};
}
std::vector<Vec2> HostUnwrapBackend::unwrap(std::span<const std::uint8_t> cut,std::optional<std::span<const MeshIndex>> charts_in,const PinMap& pins,
                                            std::optional<std::uint64_t> iterations,std::optional<std::span<const Vec2>>) {
    const auto cs=cut_surface(mesh_,cut);const auto& ct=cs->tables();
    std::vector<MeshIndex> charts;
    if(charts_in)charts.assign(charts_in->begin(),charts_in->end());else for(MeshIndex c=0;c<ct.n_charts;++c)charts.push_back(c);
    const auto geometry=quality_geometry(cs);
    std::vector<Vec2> out(ct.sector_of_corner.size(),Vec2{0,0});std::vector<MeshIndex> native;
    for(auto chart:charts) {
        bool has_pins=false;for(const auto& [c,p]:pins)has_pins=has_pins || ct.chart_of_corner[c]==chart;
        std::optional<ChartSolveResult> planar;
        if(!has_pins) {ChartSolveRequest r;r.geometry=geometry;r.chart=chart;r.validate=false;planar=planar_ngon_map(r);}
        if(planar) {for(std::size_t i=0;i<planar->corner_ids.size();++i)out[planar->corner_ids[i]]=planar->corner_uv[i];++planar_maps;}
        else native.push_back(chart);
    }
    if(native.empty())return out;
    const auto steps=iterations.value_or(iterations_);
    auto [values,used]=solve(cs,native,pins,method_,steps);primary=used;
    const std::set<MeshIndex> native_set(native.begin(),native.end());
    for(std::size_t c=0;c<out.size();++c)if(native_set.contains(ct.chart_of_corner[c]))out[c]=values[c];
    // The retry applies when the host solved with minimum stretch (RNA-discovered there).
    if(used!="MINIMUM_STRETCH")return out;
    // SLIM can badly deform coarse n-gons and tapered strips: compare a conformal solve.
    auto measure=[&](const std::vector<Vec2>& uv,MeshIndex chart) {
        const std::array<MeshIndex,1> one{chart};
        auto rep=QualitySnapshot(geometry,CornerUvView{uv},{},std::span<const MeshIndex>(one)).report();
        const auto crossings=boundary_crossings(*cs,CornerUvView{uv},{},std::span<const MeshIndex>(one)).counts[chart];
        return std::pair(rep,crossings);
    };
    struct Rank {bool invalid;std::array<double,3> values;};
    auto rank=[](const std::pair<QualityReport,std::uint64_t>& m) {
        const auto& r=m.first;const bool invalid=r.flipped || m.second || r.collapsed_uv || r.nonfinite;
        return Rank{invalid,{r.anisotropy_p95,r.area_p95,r.anisotropy_max}};
    };
    // Python tuple `<`: the first position whose items differ (`!=`, so NaN always differs)
    // decides with that item's `<`. A chart with no measurable triangle (a map that is
    // non-finite everywhere) has NaN percentiles; std::tuple would skip past them.
    auto less=[](const Rank& a,const Rank& b) {
        if(a.invalid!=b.invalid)return a.invalid<b.invalid;
        for(std::size_t i=0;i<a.values.size();++i)if(a.values[i]!=b.values[i] || std::isnan(a.values[i]))return a.values[i]<b.values[i];
        return false;
    };
    std::vector<MeshIndex> retry;std::map<MeshIndex,std::pair<QualityReport,std::uint64_t>> measured;
    for(auto chart:native) {
        measured[chart]=measure(out,chart);const auto& r=measured[chart].first;
        if(rank(measured[chart]).invalid || r.anisotropy_p95>1.3 || r.area_p95>1.4 || r.anisotropy_max>4)retry.push_back(chart);
    }
    if(!retry.empty()) {
        const auto alternative=solve(cs,retry,pins,"CONFORMAL",steps).first;
        for(auto chart:retry)
            if(less(rank(measure(alternative,chart)),rank(measured[chart]))) {
                for(std::size_t c=0;c<out.size();++c)if(ct.chart_of_corner[c]==chart)out[c]=alternative[c];
                ++conformal_fallbacks;
            }
    }
    return out;
}
RunResult to_run_result(const PipelineResult& res,const RunInput& input,const std::shared_ptr<const MeshSnapshot>& mesh,const std::string& backend) {
    auto out=run_result(res,input.options,input.pins,backend);
    if(out.geometry_fingerprint.empty())out.geometry_fingerprint=mesh->input().geometry_fingerprint.empty()?mesh->geometry_key():mesh->input().geometry_fingerprint;
    return out;
}
// Every index, length and value the worker will use is checked here, synchronously in
// start(), against the owned snapshot's validated counts (T25): the bindings build a
// RunInput without the wire schema, so nothing reaches the worker unchecked.
PipelineInput to_pipeline_input(const RunInput& input,const std::shared_ptr<const MeshSnapshot>& mesh) {
    PipelineInput in;in.mesh=mesh;const auto ne=mesh->tables().edge_vertices.size();
    const auto nv=mesh->input().positions.size(),nc=mesh->input().corner_vertices.size();
    auto bad=[](const char* field,const char* message) {throw ContractError(ErrorCode::invalid_input,field,message);};
    auto mask=[&](const std::vector<std::uint8_t>& given,const std::vector<EdgeEndpoints>& pairs,const char* field) {
        if(!given.empty()) {
            if(!pairs.empty())bad(field,"give endpoint pairs or an edge mask, not both");
            if(given.size()!=ne)throw ContractError(ErrorCode::invalid_input,field,"edge mask does not match the mesh edges");
            if(std::any_of(given.begin(),given.end(),[](auto x) {return x>1;}))bad(field,"edge mask values must be zero or one");
            return given;
        }
        EdgeMask out(ne,0);
        for(const auto& [a,b]:pairs) {
            const auto e=mesh->edge_index(a,b);
            if(e<0)throw ContractError(ErrorCode::invalid_input,field,"edge is not a mesh edge");
            out[static_cast<std::size_t>(e)]=1;
        }
        return out;
    };
    in.seams=mask(input.seam_mask,input.seam_edges,"seam_edges");
    const auto locked=mask(input.lock_mask,input.locked_edges,"locked_edges");
    if(std::any_of(locked.begin(),locked.end(),[](auto x) {return x!=0;}))in.locked=locked;
    // Non-finite incoming UVs are kept: the pipeline then measures no baseline, as the
    // source's `np.isfinite(uv_in).all()` does (Gate A F6). Pins stay finite (below).
    if(input.uv && input.uv->size()!=nc)bad("uv","UV count must equal corner count");
    in.uv=input.uv;for(const auto& p:input.pins) {
        if(p.corner<0)throw ContractError(ErrorCode::invalid_input,"pins","negative corner");
        if(static_cast<std::uint64_t>(p.corner)>=nc)bad("pins","pinned corner out of range");
        if(!std::isfinite(p.uv[0]) || !std::isfinite(p.uv[1]))bad("pins","pin coordinates must be finite");
        // A repeated corner would silently overwrite the first pin.
        if(!in.pins.emplace(static_cast<MeshIndex>(p.corner),p.uv).second)bad("pins","duplicate corner pin");
    }
    for(auto a:input.anchors)if(a<0 || static_cast<std::uint64_t>(a)>=nv)bad("anchors","anchor vertex out of range");
    pipeline_options_json(input.options); // rejects unknown enums and nonfinite or out-of-range numbers
    in.anchors=input.anchors;in.options=input.options;
    in.operation=input.operation;in.selected_faces=input.selected_faces;
    check_layout_operation(in);
    // Painted constraints (T29): typed masks and weights in core order, raw host attributes.
    // Values must be finite here, so a NaN weight is refused at start(); domains, types,
    // lengths, clamping and conflicts are the native resolver's. GB M2: only enhanced_v1 reads
    // them; under the other policies they are neither checked nor passed on, as at Gate A.
    if(input.options.policy=="enhanced_v1") {
        check_constraint_input(input.constraints,ne,nv);
        in.constraints=input.constraints;
    }
    return in;
}
}
void check_constraint_input(const ConstraintInput& k,std::uint64_t ne,std::uint64_t nv) {
    auto bad=[](const char* field,const char* message) {throw ContractError(ErrorCode::invalid_input,field,message);};
    for(const auto* m:{&k.required_edges,&k.forbidden_edges}) {
        if(!m->empty() && m->size()!=ne)bad("constraints","painted edge mask does not match the mesh edges");
        if(std::any_of(m->begin(),m->end(),[](auto x) {return x>1;}))bad("constraints","painted edge mask values must be zero or one");
    }
    for(const auto* w:{&k.allow_points,&k.avoid_points}) {
        if(!w->empty() && w->size()!=nv)bad("constraints","painted point weights do not match the mesh vertices");
        if(std::any_of(w->begin(),w->end(),[](double x) {return !std::isfinite(x);}))bad("constraints","painted weights must be finite");
    }
    for(const auto& a:k.attributes)
        if(std::any_of(a.values.begin(),a.values.end(),[](double x) {return !std::isfinite(x);}))
            throw ContractError(ErrorCode::invalid_input,"constraints.attributes."+a.name,"painted values must be finite");
}
struct Session::Impl:SessionState {};
Session::Session(SessionConfig config):impl_(std::make_unique<Impl>()) {impl_->config=resolve_session_config(std::move(config));}
Session::~Session() {close();}
void Session::start(RunInput input) {
    auto& s=*impl_;
    {
        std::lock_guard lock(s.mutex);
        if(s.state!=SessionState::State::idle)throw ContractError(ErrorCode::stale_state,"session","a session runs one input; create a new session");
        s.state=SessionState::State::running;
    }
    try {
        // Owned copies are made here, synchronously: no caller buffer is read after start().
        s.mesh=std::make_shared<const MeshSnapshot>(input.mesh,s.config.memory_budget_bytes);
        // Gate A F3: a mesh whose run tables cannot fit the session's budget is refused
        // here (memory_limit), synchronously, before any worker or stage starts.
        check_run_memory(table_counts(*s.mesh),s.config.memory_budget_bytes);
        s.geometry_revision=input.geometry_revision;s.started=std::chrono::steady_clock::now();
        auto pipeline=to_pipeline_input(input,s.mesh);pipeline.memory_budget=s.config.memory_budget_bytes;
        pipeline.scheduling.threads=s.config.worker_count;pipeline.scheduling.batch_max_charts=s.config.batch_max_charts;
        pipeline.scheduling.batch_max_triangles=s.config.batch_max_triangles;
        const auto backend=s.config.backend;
        s.worker=std::thread([this,in=std::move(pipeline),input=std::move(input),backend] {
            auto& s=*impl_;RunResult result;
            try {
                std::unique_ptr<HostUnwrapBackend> host;
                // The caller selects the primary method; the host reports the
                // method actually used. Native retry requests remain explicit.
                if(backend=="host")
                    host=std::make_unique<HostUnwrapBackend>(s,s.mesh,s.config.host_unwrap_method,static_cast<std::uint64_t>(in.options.solver.iterations));
                const ProgressSink sink=[&s](const Json& e) {s.publish(e);};
                LayoutPostprocessServices services;
                Revision fragment_revision=0;std::vector<std::string> fragment_methods;
                const bool fragment_host=std::find(s.config.host_capabilities.begin(),s.config.host_capabilities.end(),"region_proposal")!=s.config.host_capabilities.end();
                if(fragment_host) {
                    services.fragment_unwrap_method=s.config.fragment_unwrap_method;
                    services.region_proposal=[&](const MeshInput& mesh,std::span<const TransferIndex> corners,const std::function<bool()>& stopped) {
                        HostRequest request;request.kind="region_proposal";request.identity.state_revision=++fragment_revision;
                        request.region_proposal=HostRegionProposalRequest{mesh,{corners.begin(),corners.end()}};
                        auto response=s.exchange(std::move(request),stopped);
                        return response.success?std::move(response.uv):std::nullopt;
                    };
                    services.fragment_unwrap=[&](const HostSolveRequest& solve,const std::function<bool()>& stopped) {
                        HostRequest request;request.kind="solve";request.identity.state_revision=fragment_revision;request.solve=solve;
                        auto response=s.exchange(std::move(request),stopped);
                        if(response.success)fragment_methods.push_back(response.method);
                        return response.success?std::move(response.uv):std::nullopt;
                    };
                }
                const auto res=run_pipeline_with_services(in,host.get(),sink,s.control,services);
                result=to_run_result(res,input,s.mesh,backend);
                if(host) {
                    result.report["unwrap_calls"]=host->calls;
                    result.report["solver_paths"]={{"primary",host->primary.empty()?s.config.host_unwrap_method:host->primary},{"planar_maps",host->planar_maps},
                                                   {"conformal_fallbacks",host->conformal_fallbacks},{"isolated_scratch_mesh",true}};
                    // T41: the enhanced lines name the host's solver as well.
                    if(result.report.contains("messages"))result.report["messages"]=run_diagnostics(result.report);
                }
                if(in.operation=="postprocess_layout") {
                    auto& fragments=result.report["recover_fragments"];
                    fragments["host_service_available"]=fragment_host;
                    fragments["solver_backend"]=fragments.value("unwrap_calls",0)>0?"host":"not_run";
                    fragments["unwrap_methods"]=fragment_methods;
                    result.report["host_proposal_calls"]=fragments.value("proposal_calls",0);
                    result.report["unwrap_calls"]=result.report.value("unwrap_calls",std::uint64_t(0))+fragments.value("unwrap_calls",std::uint64_t(0));
                }
            } catch(const ContractError& error) {
                result.accepted=false;result.reason=std::string(error_code_name(error.code));
                result.cancelled=error.code==ErrorCode::cancelled;result.backend=backend;
                result.diagnostics.push_back({std::string(error_code_name(error.code)),error.field,error.what()});
            } catch(const std::exception& error) {
                result.accepted=false;result.reason="numerical_failure";result.backend=backend;
                result.diagnostics.push_back({"numerical_failure","",error.what()});
            }
            std::lock_guard lock(s.mutex);
            s.result=std::move(result);if(s.state!=SessionState::State::closed)s.state=SessionState::State::done;s.pending.reset();
            s.changed.notify_all();
        });
    } catch(...) {
        std::lock_guard lock(s.mutex);s.state=SessionState::State::idle;throw;
    }
}
std::vector<ProgressEvent> Session::poll_progress() {
    std::lock_guard lock(impl_->mutex);std::vector<ProgressEvent> out(impl_->progress.begin(),impl_->progress.end());impl_->progress.clear();return out;
}
std::optional<HostRequest> Session::next_host_request() {
    std::lock_guard lock(impl_->mutex);
    if(impl_->state!=SessionState::State::running || impl_->control.cancel.requested() || !impl_->pending || impl_->pending_taken || impl_->response)return std::nullopt;
    impl_->pending_taken=true;return impl_->pending;
}
void Session::submit_host_response(HostResponse response) {
    auto& s=*impl_;std::lock_guard lock(s.mutex);
    // Stale, duplicate, foreign or malformed answers are rejected; the request stays open.
    if(s.state!=SessionState::State::running || s.control.cancel.requested() || !s.pending || !s.pending_taken || s.response)throw ContractError(ErrorCode::stale_state,"host_response","no host request is waiting for this response");
    const auto& id=s.pending->identity;
    // A host built for another protocol (the F30 ABI cases) is refused outright.
    if(response.identity.version!=id.version)
        throw ContractError(ErrorCode::unsupported_version,"host_response.version","host protocol version "+std::to_string(response.identity.version)+" is not supported");
    if(response.identity.session_id!=id.session_id || response.identity.request_id!=id.request_id)
        throw ContractError(ErrorCode::stale_state,"host_response.identity","response does not answer the pending request");
    if(response.identity.geometry_revision!=id.geometry_revision)
        throw ContractError(ErrorCode::stale_state,"host_response.geometry_revision","response was computed for another geometry revision");
    if(response.identity.state_revision!=id.state_revision)
        throw ContractError(ErrorCode::stale_state,"host_response.state_revision","response was computed for another state revision");
    if(response.identity.chart_keys!=id.chart_keys)
        throw ContractError(ErrorCode::stale_state,"host_response.chart_keys","response was computed for different charts");
    if(response.kind!=s.pending->kind)throw ContractError(ErrorCode::invalid_input,"host_response.kind","response kind does not match the request");
    // A failure carries its structured error and a success none (as the wire contract).
    if(response.success==response.error.has_value())
        throw ContractError(ErrorCode::invalid_input,"host_response.error","failure requires a structured error; success forbids it");
    // The method the host reports must be one of Blender's unwrap methods (or omitted).
    if(response.kind=="region_proposal") {
        if(response.success && response.method!="SMART_PROJECT")
            throw ContractError(ErrorCode::invalid_input,"host_response.method","region proposal must report SMART_PROJECT");
    } else if(!response.method.empty() && response.method!="ANGLE_BASED" && response.method!="CONFORMAL" && response.method!="MINIMUM_STRETCH")
        throw ContractError(ErrorCode::invalid_input,"host_response.method","unknown host unwrap method");
    // A non-finite map is an answer, not a protocol error (Gate A F9): the source's
    // BlenderBackend ranks it invalid, retries CONFORMAL and the repair ladder follows.
    if(response.success && s.pending->solve) {
        if(!response.uv || response.uv->size()!=s.pending->solve->scratch_to_original_corner.size())
            throw ContractError(ErrorCode::invalid_input,"host_response.uv","solve response needs one UV per scratch corner");
    }
    if(response.success && s.pending->region_proposal) {
        if(!response.uv || response.uv->size()!=s.pending->region_proposal->scratch_to_original_corner.size())
            throw ContractError(ErrorCode::invalid_input,"host_response.uv","proposal response needs one UV per scratch corner");
    }
    s.response=std::move(response);s.changed.notify_all();
}
bool Session::ready() const {std::lock_guard lock(impl_->mutex);return impl_->result.has_value();}
RunResult Session::result() const {
    std::lock_guard lock(impl_->mutex);
    if(!impl_->result)throw ContractError(ErrorCode::stale_state,"session","the run has not finished");
    return *impl_->result;
}
RunResult Session::result_for_apply(std::string_view current_geometry_fingerprint) const {
    auto result=this->result();
    // The host recomputes the native fingerprint of its live data right before writing.
    // A run without a host fingerprint, or data changed since start(), is never applied.
    const auto& started=impl_->mesh->input().geometry_fingerprint;
    if(started.empty())
        throw ContractError(ErrorCode::stale_state,"apply.geometry_fingerprint","the run has no host geometry fingerprint to verify");
    if(current_geometry_fingerprint!=started)
        throw ContractError(ErrorCode::stale_state,"apply.geometry_fingerprint","the geometry changed after the run started; nothing applied");
    // Pins are compared exactly as given (the host's stored values), never after aligning
    // maps (Gate A F28, formerly host glue). HOST-01: the 0.1.0 add-on accepted a map whose
    // pins are within its pipeline tolerance (np.allclose) and wrote the solver's values, so
    // a host unwrapper that returns a pin a few float32 steps off (Blender CONFORMAL) was
    // applied there and refused here. Such a pin is now put back exactly on its stored value
    // before writing, provided no face of the map flips or collapses; a pin beyond that
    // tolerance is still refused and nothing is applied.
    if(result.accepted && result.uv) {
        auto& uv=*result.uv;
        std::vector<std::pair<std::size_t,Vec2>> snaps;double deviation=0;
        for(const auto& pin:result.pins) {
            const auto c=static_cast<std::size_t>(pin.corner);
            if(pin.corner<0 || c>=uv.size())
                throw ContractError(ErrorCode::constraint_conflict,"apply.pins","the result moved a pinned UV; nothing applied");
            if(uv[c]==pin.uv)continue;
            for(std::size_t k=0;k<2;++k) {
                const double d=std::abs(uv[c][k]-pin.uv[k]);
                if(!(d<=1e-8+1e-5*std::abs(pin.uv[k])))
                    throw ContractError(ErrorCode::constraint_conflict,"apply.pins","the result moved a pinned UV; nothing applied");
                deviation=std::max(deviation,d);
            }
            snaps.emplace_back(c,pin.uv);
        }
        if(!snaps.empty()) {
            const auto& mesh=impl_->mesh->input();
            const auto& off=mesh.face_offsets;
            auto area=[&](std::size_t f,const std::vector<Vec2>& map) {
                double a=0;const auto b=static_cast<std::size_t>(off[f]),e=static_cast<std::size_t>(off[f+1]);
                for(auto i=b;i<e;++i) {const auto& p=map[i];const auto& q=map[i+1<e?i+1:b];a+=p[0]*q[1]-q[0]*p[1];}
                return a;
            };
            std::vector<std::size_t> corner_face(uv.size());
            for(std::size_t f=0;f+1<off.size();++f)
                for(auto i=static_cast<std::size_t>(off[f]);i<static_cast<std::size_t>(off[f+1]) && i<uv.size();++i)corner_face[i]=f;
            auto snapped=uv;for(const auto& [c,p]:snaps)snapped[c]=p;
            std::set<std::size_t> faces;for(const auto& [c,p]:snaps)faces.insert(corner_face[c]);
            for(const auto f:faces) {
                const double before=area(f,uv),after=area(f,snapped);
                if(after==0 || (before>0)!=(after>0))
                    throw ContractError(ErrorCode::constraint_conflict,"apply.pins","putting a pinned UV back on its pin would flip or collapse a face; nothing applied");
            }
            uv=std::move(snapped);
            result.report["host_apply"]={{"pins_restored",snaps.size()},{"max_pin_deviation",deviation},
                                         {"rule","pins within the pipeline's tolerance (|d| <= 1e-8 + 1e-5 |pin|) are written exactly at their stored value (HOST-01)"}};
        }
    }
    return result;
}
void check_stored_uv(std::span<const Vec2> uv,std::span<const TransferIndex> host_loop_indices,std::span<const float> stored) {
    if(host_loop_indices.size()!=uv.size())
        throw ContractError(ErrorCode::invalid_input,"apply.host_loop_indices","one host loop per map corner required");
    if(stored.size()%2!=0)throw ContractError(ErrorCode::invalid_input,"apply.stored_uv","two stored values per host loop required");
    const auto loops=stored.size()/2;
    for(std::size_t c=0;c<uv.size();++c) {
        const auto loop=host_loop_indices[c];
        if(loop<0 || static_cast<std::uint64_t>(loop)>=loops)
            throw ContractError(ErrorCode::invalid_input,"apply.host_loop_indices","host loop outside the stored layer");
        const auto l=static_cast<std::size_t>(loop);
        // Round-to-nearest float32 conversion: exactly what storing the float64 map produces.
        if(stored[2*l]!=static_cast<float>(uv[c][0]) || stored[2*l+1]!=static_cast<float>(uv[c][1]))
            throw ContractError(ErrorCode::host_failure,"apply.uv_storage","UV transfer to the host changed values beyond float32 rounding");
    }
}
void Session::request_cancel() {
    impl_->control.cancel.request();std::lock_guard lock(impl_->mutex);impl_->changed.notify_all();
}
void Session::close() {
    if(!impl_)return;
    // Concurrent close() calls (a host thread and a module shutdown hook) serialize on
    // join_mutex: exactly one joins, the others wait until the worker has ended.
    std::lock_guard join(impl_->join_mutex);
    {
        std::lock_guard lock(impl_->mutex);
        if(impl_->state==SessionState::State::closed && !impl_->worker.joinable())return;
        impl_->control.cancel.request();impl_->state=SessionState::State::closed;impl_->changed.notify_all();
    }
    if(impl_->worker.joinable() && impl_->worker.get_id()!=std::this_thread::get_id())impl_->worker.join();
}
}
