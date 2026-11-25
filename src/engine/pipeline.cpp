// SPDX-License-Identifier: GPL-3.0-or-later
// Source core/pipeline.py: CoreBackend, backend_evaluator, backend_interior_solver,
// gate, passthrough_faces, _metrics, _contour_validator, _run and run.
#include "cleanseams/pipeline.hpp"
#include "cleanseams/layout_postprocess.hpp"
#include "pipeline_services.hpp"
#include "cleanseams/cache.hpp"
#include "cleanseams/errors.hpp"
#include "cleanseams/detail/numpy_compat.hpp"
#include "cleanseams/classify.hpp"
#include "cleanseams/necessity.hpp"
#include "cleanseams/symmetric_moves.hpp"
#include "pipeline_hooks.hpp"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <ctime>
#include <set>
namespace cleanseams {
namespace {
using Json=nlohmann::json;
using Clock=MonotonicClock;
std::vector<TransferIndex> face_list(std::span<const std::uint8_t> mask) {
    std::vector<TransferIndex> out;for(std::size_t f=0;f<mask.size();++f)if(mask[f])out.push_back(static_cast<TransferIndex>(f));return out;
}
QualityPolicy exempt_policy(std::span<const std::uint8_t> exempt) {QualityPolicy p;p.exempt_faces=face_list(exempt);return p;}
std::size_t count(std::span<const std::uint8_t> mask) {return static_cast<std::size_t>(std::count_if(mask.begin(),mask.end(),[](auto x) {return x!=0;}));}
double cpu_seconds() {return static_cast<double>(std::clock())/CLOCKS_PER_SEC;}
double seconds_since(Clock::time_point start) {return std::chrono::duration<double>(Clock::now()-start).count();}
std::vector<MeshIndex> chart_corners_in_face_order(const CutSurface& cs,MeshIndex chart) {
    std::vector<MeshIndex> out;for(auto f:cs.tables().chart_faces[chart])for(auto c:cs.mesh().face_corners(f))out.push_back(c);return out;
}
}
PinCompletion pin_policy(const PipelineOptions& o) {return o.policy=="legacy_010"?PinCompletion::legacy_010:PinCompletion::preserve;}
UnwrapBackend::~UnwrapBackend()=default;
CoreUnwrapBackend::CoreUnwrapBackend(std::shared_ptr<const MeshSnapshot> mesh,PipelineOptions options,CancellationToken cancel):
    mesh_(mesh),options_(std::move(options)),cancel_(cancel),backend_(std::move(mesh),[&] {
        NativeArapBackend::Settings s;s.iterations=10;s.pin_completion=pin_policy(options_);s.cancel=cancel_;return s;}()) {}
std::vector<Vec2> CoreUnwrapBackend::unwrap(std::span<const std::uint8_t> cut,std::optional<std::span<const MeshIndex>> charts,const PinMap& pins,
                                            std::optional<std::uint64_t> iterations,std::optional<std::span<const Vec2>> init) {
    std::vector<CornerPin> list;for(const auto& [c,p]:pins)list.push_back({static_cast<TransferIndex>(c),p});
    return backend_.unwrap(cut,charts,list,iterations,init).uv;
}
ChartEvaluator CoreUnwrapBackend::evaluator() const {
    EvaluatorSettings s;s.pin_completion=pin_policy(options_);s.cancel=cancel_;return default_evaluator(mesh_,options_.layout,s);
}
ChartEvaluator CoreUnwrapBackend::guard_evaluator() const {
    EvaluatorSettings s;s.pin_completion=pin_policy(options_);s.iterations=iterations();s.cancel=cancel_;return default_evaluator(mesh_,options_.layout,s);
}
InteriorSolver CoreUnwrapBackend::interior_solver() const {
    // Gate A F30: the run's pin completion, cancellation token and memory allowance.
    InteriorSolverSettings s;s.pin_completion=pin_policy(options_);s.cancel=cancel_;s.memory_budget=active_memory_budget();
    return default_interior_solver(options_.boundary,s);
}
ChartEvaluator backend_evaluator(std::shared_ptr<const MeshSnapshot> mesh,UnwrapBackend& backend,const PipelineOptions& options) {
    return [mesh,&backend,layout=options.layout](const std::shared_ptr<const CutSurface>& cs,std::span<const MeshIndex> charts,
                                                  const std::map<MeshIndex,double>& relax) {
        const std::vector<MeshIndex> list(charts.begin(),charts.end());
        // The solve itself stays on this thread (a host backend answers on the host's
        // main thread); only the native measurements of the returned map run in parallel.
        const auto uv=backend.unwrap(cs->tables().cut,std::span<const MeshIndex>(list),PinMap{});
        const auto geometry=quality_geometry(cs);ChartMeasurements out;
        std::vector<std::optional<ChartMeasurement>> slots(list.size());
        parallel_for(list.size(),[&](std::size_t i) {
            const auto ch=list[i];const std::array<MeshIndex,1> one{ch};ChartMeasurement m;
            m.report=QualitySnapshot(geometry,CornerUvView{uv},{},std::span<const MeshIndex>(one)).report();
            m.crossings=boundary_crossings(*cs,CornerUvView{uv},{},std::span<const MeshIndex>(one)).counts[ch];
            m.corner_ids=chart_corners_in_face_order(*cs,ch);for(auto c:m.corner_ids)m.corner_uv.push_back(uv.at(c));
            const auto found=relax.find(ch);m.within_budget=within_budget(m.report,m.crossings,layout,found==relax.end()?1.0:found->second);
            m.method="backend";slots[i]=std::move(m);
        });
        for(std::size_t i=0;i<list.size();++i)out[list[i]]=std::move(*slots[i]);
        return out;
    };
}
InteriorSolver backend_interior_solver(std::shared_ptr<const MeshSnapshot> mesh,UnwrapBackend& backend) {
    return [mesh,&backend](const std::shared_ptr<const QualityGeometry>& geometry,MeshIndex chart,std::span<const SectorPin> pins,std::span<const Vec2> uv) {
        // pins are given per sector; the backend takes them per corner
        const auto& cs=geometry->cut();std::map<MeshIndex,Vec2> by_sector;for(const auto& p:pins)by_sector[p.sector]=p.uv;
        const auto ids=chart_corners_in_face_order(cs,chart);PinMap corner_pins;bool free=false;
        for(auto c:ids) {
            const auto found=by_sector.find(cs.tables().sector_of_corner[c]);
            if(found!=by_sector.end())corner_pins[c]=found->second;else free=true;
        }
        std::vector<Vec2> values;
        if(!free) {
            // nothing to solve (a single polygon): the contour is the whole map
            std::vector<Vec2> out(uv.begin(),uv.end());for(const auto& [c,p]:corner_pins)out[c]=p;
            for(auto c:ids)values.push_back(out[c]);
            return std::pair(ids,values);
        }
        const std::array<MeshIndex,1> one{chart};
        const auto solved=backend.unwrap(cs.tables().cut,std::span<const MeshIndex>(one),corner_pins);
        for(auto c:ids)values.push_back(solved.at(c));
        return std::pair(ids,values);
    };
}
GateResult gate(const MeshSnapshot& m,const PipelineOptions& o) {
    const auto d=m.diagnose();const auto nf=static_cast<std::int64_t>(m.input().face_offsets.size()-1);
    GateResult out;
    out.info={{"faces",nf},{"components",d.at("components")},
              {"largest_component_fraction",detail::python_round(d.at("largest_component_fraction").get<double>(),4)},
              {"degenerate",d.at("degenerate")},{"degenerate_by_kind",d.at("by_kind")},
              {"degenerate_fraction",detail::python_round(d.at("degenerate_fraction").get<double>(),4)},{"first_degenerate",d.at("first_degenerate")}};
    if(nf>o.gate_max_faces) {out.ok=false;out.reason="too_many_faces";}
    else if(d.at("components").get<std::int64_t>()>o.gate_max_components ||
            d.at("largest_component_fraction").get<double>()<o.gate_min_main_fraction) {out.ok=false;out.reason="fragment_cloud";}
    else if(d.at("degenerate_fraction").get<double>()>o.gate_max_degenerate) {out.ok=false;out.reason="scan_debris";}
    return out;
}
std::vector<std::uint8_t> passthrough_faces(const MeshSnapshot& m,const PipelineOptions& o) {
    const auto& g=m.geometry();const auto nf=g.face_area.size();std::vector<std::uint8_t> out(nf,0);
    if(nf==0 || g.n_components<=1)return out;
    std::vector<double> area(g.n_components,0.0);for(std::size_t f=0;f<nf;++f)area[g.face_component[f]]+=g.face_area[f];
    const double floor=o.passthrough_fraction*detail::numpy_sum(g.face_area);
    for(std::size_t f=0;f<nf;++f)out[f]=area[g.face_component[f]]<floor;
    return out;
}
Json pipeline_metrics(const std::shared_ptr<const CutSurface>& cs,std::span<const Vec2> uv,std::span<const std::uint8_t> exempt) {
    return QualityEvaluation(quality_geometry(cs),CornerUvView{uv},exempt_policy(exempt)).metrics();
}
namespace {
// Source _contour_validator: the contour cleaner's validator with the run's exempt faces.
ContourValidator contour_validator(const PipelineOptions& o,std::vector<std::uint8_t> exempt) {
    return [o,exempt=std::move(exempt)](const std::shared_ptr<const QualityGeometry>& geometry,MeshIndex chart,std::span<const Vec2> uv,
                                        const ContourMeasurement& baseline) {
        const auto& bo=o.boundary;const std::array<MeshIndex,1> one{chart};
        ContourMeasurement m;const auto policy=exempt_policy(exempt);
        m.report=QualitySnapshot(geometry,CornerUvView{uv},policy,std::span<const MeshIndex>(one)).report();
        BoundaryCrossingOptions ignore;ignore.ignored_faces=policy.exempt_faces;
        m.crossings=boundary_crossings(geometry->cut(),CornerUvView{uv},ignore,std::span<const MeshIndex>(one)).counts[chart];
        const auto& r=m.report;const auto& b=baseline.report;
        if(r.nonfinite || r.flipped || m.crossings || r.collapsed_uv)return std::pair(false,m);
        if(o.target_islands && !within_budget(r,m.crossings,o.layout))return std::pair(false,m);
        if(r.anisotropy_p95>std::min(b.anisotropy_p95*bo.p95_growth,bo.abs_anisotropy_p95))return std::pair(false,m);
        if(r.area_p95>std::min(b.area_p95*bo.p95_growth,bo.abs_area_p95))return std::pair(false,m);
        if(r.anisotropy_max>b.anisotropy_max*bo.max_growth)return std::pair(false,m);
        return std::pair(true,m);
    };
}
}
// GB N3 (supersedes the M4 rows): the chart_over_budget rows of an adopted enhanced state by
// native_baseline's own listing rule, so that `budget_satisfied` changes only when a listed chart
// changes. `native_rows` are native_baseline's rows of candidate zero (`origin`), `source` the
// rule that produced them: "chart_gate" (the final chart gate: P95 anisotropy over
// max_anisotropy_p95 * final_gate_tolerance on charts of at least split_min_faces faces),
// "repair_absorption" (the recount after absorbed repair fragments: the same P95 rule on every
// chart) or "not_judged" (no gate ran: no rows). Both add FIX-04's area P95 rows (area budget,
// same face floor) for charts not listed by anisotropy; the P99 and maximum caps of the enhanced
// budget are not part of native's rule. A chart the pass left exactly as candidate zero had it
// (same chart key and map) keeps native_baseline's row or absence; every other chart of the kept
// map is judged by that rule on its measurement. Rows are [chart, P95 anisotropy to 3 places
// (null when not finite)] in the kept map's chart order. report.enhanced.retained.budget stays the
// enhanced budget (enhanced_chart_budget).
nlohmann::json native_rule_over_budget_rows(const RunState& kept,const RunState& origin,const nlohmann::json& native_rows,const std::string& source,
                                            const PipelineOptions& o,nlohmann::json* record_out) {
    if(!kept.surface || !origin.surface)throw ContractError(ErrorCode::invalid_input,"enhanced.budget_rows","measured states are required");
    if(source!="chart_gate" && source!="repair_absorption" && source!="not_judged")
        throw ContractError(ErrorCode::invalid_input,"enhanced.budget_rows","unknown listing rule "+source);
    struct {Json rows=Json::array();Json record;} out;
    const double budget=o.layout.max_anisotropy_p95*o.final_gate_tolerance,area_budget=o.layout.max_area_p95*o.final_gate_tolerance;
    std::map<MeshIndex,Json> listed;
    for(const auto& row:native_rows)if(row.is_array() && !row.empty() && row[0].is_number_unsigned())listed[row[0].get<MeshIndex>()]=row;
    std::map<std::string,MeshIndex> by_key;for(MeshIndex ch=0;ch<origin.charts.size();++ch)by_key.emplace(origin.charts[ch].key,ch);
    // Charts of the kept map that are candidate zero's own: same key, every corner's UV unchanged.
    const auto& kt=kept.surface->tables();std::vector<std::uint8_t> moved(kt.n_charts,0);
    const bool same_size=kept.uv.size()==origin.uv.size();
    for(std::size_t c=0;c<kept.uv.size();++c)if(!same_size || kept.uv[c]!=origin.uv[c])moved[kt.chart_of_corner[c]]=1;
    Json inherited=Json::array(),measured=Json::array();
    for(MeshIndex ch=0;ch<kept.charts.size();++ch) {
        const auto& c=kept.charts[ch];
        const auto found=moved[ch]?by_key.end():by_key.find(c.key);
        if(found!=by_key.end()) {
            if(const auto row=listed.find(found->second);row!=listed.end()) {Json r=row->second;r[0]=ch;out.rows.push_back(std::move(r));inherited.push_back(ch);}
            continue;
        }
        if(source=="not_judged" || c.exempt_only || !c.measured)continue;
        const bool floor=static_cast<std::int64_t>(c.faces)>=o.layout.split_min_faces;
        const bool anisotropy=std::isfinite(c.anisotropy_p95) && c.anisotropy_p95>budget && (source=="repair_absorption" || floor);
        const bool area=!anisotropy && o.policy!="legacy_010" && floor && std::isfinite(c.area_p95) && c.area_p95>area_budget; // FIX-04
        if(!anisotropy && !area)continue;
        out.rows.push_back(Json::array({ch,std::isfinite(c.anisotropy_p95)?Json(detail::python_round(c.anisotropy_p95,3)):Json(nullptr)}));
        measured.push_back(ch);
    }
    out.record={{"rule","native_baseline"},{"source",source},{"anisotropy_budget",budget},{"area_budget",area_budget},
                {"split_min_faces",o.layout.split_min_faces},{"inherited_charts",inherited},{"measured_charts",measured},
                {"note","rows of charts unchanged since candidate zero are native_baseline's; other charts are judged by the rule of `source` "
                        "(chart_gate: P95 on charts of at least split_min_faces faces; repair_absorption: P95 on every chart; FIX-04 area rows); "
                        "report.enhanced.retained.budget is the enhanced budget"}};
    if(record_out)*record_out=std::move(out.record);
    return out.rows;
}
namespace {
// The charts holding pass-through faces (repairs_summary.passthrough_charts) on a surface.
Json passthrough_charts_of(const CutSurface& cs,std::span<const std::uint8_t> passthrough) {
    std::set<MeshIndex> charts;for(std::size_t f=0;f<passthrough.size();++f)if(passthrough[f])charts.insert(cs.tables().chart_of_face[f]);
    return std::vector<MeshIndex>(charts.begin(),charts.end());
}
// One run of the source's _run with the report it builds.
class Run {
public:
    Run(const PipelineInput& input,UnwrapBackend* backend,const ProgressSink& progress,const Control& control):
        in_(input),m_(*input.mesh),o_(input.options),progress_(progress),control_(control),backend_(backend) {}
    PipelineResult execute();
private:
    const PipelineInput& in_;
    const MeshSnapshot& m_;
    PipelineOptions o_;
    const ProgressSink& progress_;
    const Control& control_;
    UnwrapBackend* backend_;
    std::unique_ptr<CoreUnwrapBackend> core_;
    std::shared_ptr<BatchLog> batches_=std::make_shared<BatchLog>();
    std::shared_ptr<ResolvedPreferences> prefs_; // enhanced_v1 painted constraints and preferences (T29)
    std::vector<SplitRecord> splits_;             // enhanced_v1: the ordinary gate's over-budget splits (T32, read-only record)
    Json retained_over_budget_=Json::array();     // enhanced_v1: the kept state's over-budget rows (GB M4, N3: native's rule)
    Clock::time_point started_=Clock::now();
    Deadline deadline_;
    std::string active_stage_="gate";
    Json rep_=Json::object();
    Deadline limit(double fraction) const {
        if(!deadline_.at())return Deadline{};
        const auto span=std::chrono::duration<double>(*deadline_.at()-started_)*fraction;
        return Deadline(started_+std::chrono::duration_cast<Clock::duration>(span));
    }
    void emit(const Json& event) {
        control_.check();
        Json message=event;message["elapsed"]=seconds_since(started_);
        const auto truthy=[&](const char* k) {return message.contains(k) && message.at(k).is_number() && message.at(k).get<double>()!=0;};
        if(truthy("done") && truthy("total")) {
            const double elapsed=message.at("elapsed").get<double>(),done=message.at("done").get<double>(),total=message.at("total").get<double>();
            message["eta"]=std::max(0.0,elapsed*(total-done)/done);
        }
        if(progress_)progress_(message);
    }
    void stage(Json record) {rep_["stages"].push_back(record);emit(record);}
    // valid: every check including the strict overlap check; admissible: valid, or feasible with
    // only candidate zero's unrepaired strict-only defects left (reported, never a refusal
    // native_baseline would not make: Section 11.3 policy).
    struct EnhancedOutcome {bool ran=false,adopted=false,feasible=true,valid=true,admissible=true;};
    // enhanced_v1: the enhancement pass on the complete candidate zero `state` (T28).
    EnhancedOutcome enhance(RepairState& state,const PinMap& pins,const EdgeMask& locked,bool baseline_accepted,const std::string& baseline_reason);
    // The pass itself; GB m4: enhance() falls back to native_baseline's result when it fails.
    EnhancedOutcome enhance_pass(RepairState& state,const PinMap& pins,const EdgeMask& locked,Clock::time_point t0,double cpu0);
    void say(const std::string& msg) {
        for(const auto& [word,name]:std::initializer_list<std::pair<const char*,const char*>>{{"generating","generate"},{"refining","refine"},
                {"unwrapping","unwrap"},{"cleaning","contours"},{"repairing","repair"},{"over budget","chart_gate"},{"enhancing","enhanced"}})
            if(msg.find(word)!=std::string::npos) {active_stage_=name;break;}
        emit({{"stage",active_stage_},{"status","running"},{"msg",msg}});
    }
};
PipelineResult Run::execute() {
    const auto t_all=Clock::now();const auto& mesh=in_.mesh;const auto nf=m_.input().face_offsets.size()-1,ne=m_.tables().edge_vertices.size();
    const auto band=pipeline_hooks::island_band(o_.target_islands);
    if(o_.profile!="custom" && o_.profile!="auto")
        throw ContractError(ErrorCode::invalid_option,"profile","profile must be custom or auto; explicit settings belong in layout and boundary");
    deadline_=Deadline::after(o_.deadline_s,started_);
    if(control_.deadline.at() && (!deadline_.at() || *control_.deadline.at()<*deadline_.at()))deadline_=control_.deadline;
    PipelineResult res;auto& rep=rep_;
    rep["timings"]=Json::object();rep["cpu_timings"]=Json::object();rep["stages"]=Json::array();rep["repairs"]=Json::array();
    auto& timings=rep["timings"];auto& cpu_timings=rep["cpu_timings"];
    {
        auto resolved=resolve_auto_profile(*mesh,o_,deadline_);
        if(!resolved.report.is_null()) {o_=std::move(resolved.options);rep["profile"]=std::move(resolved.report);}
    }
    if(in_.seams.size()!=ne)throw ContractError(ErrorCode::invalid_input,"seams","edge mask does not match the mesh edges");
    auto seams=in_.seams;EdgeMask locked=in_.locked?*in_.locked:EdgeMask(ne,0);
    if(locked.size()!=ne)throw ContractError(ErrorCode::invalid_input,"locked","edge mask does not match the mesh edges");
    const auto& pins=in_.pins;
    for(const auto& [c,p]:pins)
        if(c>=m_.tables().corner_face.size() || !std::isfinite(p[0]) || !std::isfinite(p[1]))throw ContractError(ErrorCode::invalid_input,"pins","pins need a valid corner and finite coordinates");
    // ------------------------------------------ enhanced_v1 constraints (T29)
    // Painted constraints and preferences are resolved first (Section 21). A hard conflict
    // refuses the run with the exact IDs, before any stage. Required edges join the input
    // seams as locked seams and forbidden edges leave them, so the baseline starts from an
    // admissible cut; the enhancement pass enforces them on every candidate.
    if(o_.policy=="enhanced_v1") {
        prefs_=std::make_shared<ResolvedPreferences>(resolve_preferences(m_,in_.constraints,o_.enhanced.preferences,in_.locked,in_.anchors,
                                                                          o_.layout.feature_angle,o_.layout.hide_dir));
        if(!prefs_->conflicts.empty()) {
            rep["enhanced"]={{"policy",enhanced_policy_version},{"status","refused"},{"reason","constraint_conflict"},{"preferences",prefs_->summary()}};
            stage({{"stage","constraints"},{"status","refused"},{"reason","constraint_conflict"},{"conflicts",prefs_->conflicts}});
            res.accepted=false;res.reason="constraint_conflict";res.seams=seams;
            timings["total"]=seconds_since(t_all);res.report=std::move(rep_);res.scheduling=batches_->json();return res;
        }
        for(std::size_t e=0;e<ne;++e) {
            if(prefs_->required[e]) {seams[e]=1;locked[e]=1;}
            else if(prefs_->forbidden[e])seams[e]=0;
        }
    }
    const bool core=backend_==nullptr;
    ChartEvaluator evaluate;InteriorSolver interior;
    if(core) {core_=std::make_unique<CoreUnwrapBackend>(mesh,o_,control_.cancel);backend_=core_.get();evaluate=core_->evaluator();interior=core_->interior_solver();}
    else {evaluate=backend_evaluator(mesh,*backend_,o_);interior=backend_interior_solver(mesh,*backend_);}
    const auto iterations=backend_->iterations()?backend_->iterations():10;
    const auto policy=pin_policy(o_);
    // ------------------------------------------------------------------ gate
    auto t0=Clock::now();auto cpu0=cpu_seconds();
    const auto gated=gate(m_,o_);rep["gate"]=gated.info;
    std::vector<std::uint8_t> exempt(nf,0);for(std::size_t f=0;f<nf;++f)exempt[f]=m_.geometry().face_degenerate_kind[f]!=0;
    const auto degenerate=count(exempt);
    const auto passthrough=passthrough_faces(m_,o_);for(std::size_t f=0;f<nf;++f)exempt[f]=exempt[f] || passthrough[f];
    rep["exempt_faces"]={{"degenerate",degenerate},{"passthrough",count(passthrough)}};
    timings["gate"]=seconds_since(t0);cpu_timings["gate"]=cpu_seconds()-cpu0;
    if(o_.gate && !gated.ok) {
        res.accepted=false;res.reason=gated.reason;res.exempt=exempt;
        stage({{"stage","gate"},{"status","refused"},{"reason",gated.reason}});
        timings["total"]=seconds_since(t_all);res.report=std::move(rep_);res.scheduling=batches_->json();return res;
    }
    // -------------------------------------------------------------- baseline
    t0=Clock::now();cpu0=cpu_seconds();
    auto cs=cut_surface(mesh,seams);
    Json baseline=nullptr;
    if(in_.uv && cs->tables().n_charts>0 &&
       std::all_of(in_.uv->begin(),in_.uv->end(),[](const Vec2& p) {return std::isfinite(p[0]) && std::isfinite(p[1]);}))
        baseline=pipeline_metrics(cs,*in_.uv,exempt);
    rep["baseline"]=baseline;
    rep["input"]={{"seams",count(seams)},{"charts",cs->tables().n_charts},{"disks",count(cs->tables().chart_is_disk)},{"mesh",m_.summary()}};
    timings["baseline"]=seconds_since(t0);cpu_timings["baseline"]=cpu_seconds()-cpu0;
    // ------------------------------------------------------------- seam cost
    t0=Clock::now();cpu0=cpu_seconds();
    std::optional<std::vector<double>> cost;
    if(o_.layout.use_exposure && !limit(0.55).expired()) {
        try {cost=seam_cost_field(m_,o_.layout);}
        catch(const std::exception& error) {rep["exposure_error"]=std::string("ContractError: ")+error.what();}
    }
    timings["exposure"]=seconds_since(t0);cpu_timings["exposure"]=cpu_seconds()-cpu0;
    // -------------------------------------------------------------- generate
    std::map<std::string,CachedChartMap> uv_cache;bool generation_timeout=false;std::optional<std::vector<Vec2>> uv_generated;
    const auto emit_sink=[this](const Json& event) {emit(event);};
    if(o_.generate) {
        say("generating seam layout");
        t0=Clock::now();cpu0=cpu_seconds();
        EdgeMask initial(ne);for(std::size_t e=0;e<ne;++e)initial[e]=seams[e] || locked[e];
        SeedingInput request;request.initial_cut=initial;request.locked=locked;request.cost=cost;request.evaluate=&evaluate;
        request.deadline=limit(0.55);request.pins=pins;request.pin_policy=policy;request.progress=emit_sink;
        request.scheduling=in_.scheduling;request.cancel=control_.cancel;request.batches=batches_;
        GeneratedLayout generated;
        if(!band)generated=generate_seams(mesh,o_.layout,request);
        else {
            auto target=pipeline_hooks::generate_for_target(mesh,o_,*band,core,backend_,request);
            generated=std::move(target.generated);rep["island_target"]=std::move(target.info);o_.layout=target.layout;
            if(core_)core_->set_options(o_);
            evaluate=core?core_->evaluator():backend_evaluator(mesh,*backend_,o_);
        }
        auto& lay=generated.report;
        generation_timeout=lay.timed_out;
        rep["layout"]=lay.as_dict();
        for(auto f:lay.exempt_faces)exempt[f]=1;
        stage({{"stage","generate"},{"status",lay.timed_out?"timed_out":"done"},{"charts",generated.cs->tables().n_charts},
               {"seam_edges",count(generated.seams)}});
        seams=generated.seams;
        if(lay.protected_seams)for(std::size_t e=0;e<ne;++e)locked[e]=locked[e] || (*lay.protected_seams)[e];
        uv_cache=lay.uv_cache;
        timings["generate"]=seconds_since(t0);cpu_timings["generate"]=cpu_seconds()-cpu0;
        uv_generated=lay.uv;
    }
    // ---------------------------------------------------------------- refine
    const ChartEvaluator guard=o_.refine_guard?(core?core_->guard_evaluator():evaluate):ChartEvaluator{};
    const ChartEvaluator* guard_ptr=guard?&guard:nullptr;
    if(o_.refine_paths && !generation_timeout) {
        say("refining seam chains");
        t0=Clock::now();cpu0=cpu_seconds();
        cs=cut_surface(mesh,seams);
        Json chain_results=Json::array();std::int64_t n_opt=0;FlankingCache previous;RoutingWorkspace workspace;
        auto paths=o_.paths;
        if(band)paths.require_budget=true;
        std::size_t quads=0;for(std::size_t f=0;f<nf;++f)quads+=m_.input().face_offsets[f+1]-m_.input().face_offsets[f]==4;
        if(static_cast<double>(quads)/static_cast<double>(nf)<.6) {
            paths.w_flow=0.;paths.w_turn=2*paths.w_turn;
            // The fixed variants exist only under legacy_010; FIX-01 derives them from the weights above.
            if(o_.policy=="legacy_010")for(auto& v:paths.variants)v={2*v[0],0.,v[2]};
        }
        std::int64_t passes=0,deferred=0;
        const std::vector<TransferIndex> anchors(in_.anchors.begin(),in_.anchors.end());
        for(std::int64_t sweep=0;sweep<std::max<std::int64_t>(1,std::min<std::int64_t>(3,o_.refine_passes));++sweep) {
            if(limit(0.8).expired())break;
            ++passes;
            const SeamGraph sg(mesh,seams,locked,anchors,cs);
            std::vector<const Chain*> chains;for(const auto& chain:sg.chains())chains.push_back(&chain);
            std::stable_sort(chains.begin(),chains.end(),[&](const Chain* a,const Chain* b) {
                return std::pair(-a->edge_length(m_),a->id)<std::pair(-b->edge_length(m_),b->id);
            });
            std::set<MeshIndex> dirty;const auto optimized_before=n_opt;
            for(std::size_t index=0;index<chains.size();++index) {
                const auto& chain=*chains[index];
                control_.check();
                if(limit(0.8).expired())break;
                if(std::any_of(chain.verts.begin(),chain.verts.end(),[&](auto v) {return dirty.contains(v);})) {++deferred;continue;}
                Json r;EdgeMask seams_new=seams;RefineResult result;
                if(chain.edges.size()<3)r={{"chain",chain.id},{"status","kept"},{"reason","too_short"}};
                else {
                    RefineRequest request;request.cut=cs;request.seams=seams;request.options=paths;request.protected_edges=std::span<const std::uint8_t>(locked);
                    if(cost)request.exposure=*cost;
                    request.evaluate=guard_ptr;request.previous=&previous;request.workspace=&workspace;
                    request.fixes.weight_variants=o_.policy!="legacy_010"; // FIX-01 (native_baseline)
                    request.fixes.corridor_bounds=o_.policy!="legacy_010"; // FIX-02
                    request.fixes.full_schedule=o_.policy!="legacy_010";   // FIX-03
                    request.speculative=core; // GB performance: the native guard measures charts independently
                    auto outcome=chain.cycle && o_.refine_cycles?refine_cycle(request,chain):refine_chain(request,chain);
                    result=std::move(outcome.result);seams_new=std::move(outcome.seams);r=result.to_json();r.erase("tried");
                }
                r["pass"]=sweep+1;
                chain_results.push_back(r);
                if(r.at("status")=="optimized") {
                    seams=std::move(seams_new);cs=result.cut;dirty.insert(result.dirty.begin(),result.dirty.end());++n_opt;
                }
                emit({{"stage","refine"},{"status","running"},{"done",index+1},{"total",chains.size()},{"pass",sweep+1}});
            }
            if(n_opt==optimized_before)break;
        }
        std::uint64_t spurs=0;
        if(o_.layout.spur_removal && !limit(0.8).expired()) {
            const SeamGraph sg(mesh,seams,locked,anchors,cs);
            if(std::any_of(sg.chains().begin(),sg.chains().end(),[](const Chain& chain) {return chain.free_end;})) {
                LayoutReport cleanup;
                auto swept=final_spur_sweep(mesh,cs->tables().cut,cs,o_.layout,cleanup,guard_ptr,std::span<const std::uint8_t>(locked),nullptr,
                                            policy==PinCompletion::preserve?pins:PinMap{});
                cs=swept.cs;for(std::size_t e=0;e<ne;++e)seams[e]=swept.cut[e] && !m_.tables().edge_boundary[e];
                spurs=cleanup.spurs_removed;
            }
        }
        rep["chains"]=chain_results;
        std::size_t kept=0;for(const auto& c:chain_results)kept+=c.at("status")=="kept";
        stage({{"stage","refine"},{"status",limit(0.8).expired()?"timed_out":"done"},{"chains",chain_results.size()},{"passes",passes},
               {"deferred",deferred},{"spurs_removed",spurs},{"optimized",n_opt},{"kept",kept}});
        timings["refine"]=seconds_since(t0);cpu_timings["refine"]=cpu_seconds()-cpu0;
        if(n_opt || spurs)uv_generated.reset(); // the layout changed: re-solve
    }
    // ---------------------------------------------------------------- unwrap
    say("unwrapping");
    t0=Clock::now();cpu0=cpu_seconds();
    RepairState state;state.seams=seams;state.exempt=exempt;state.cs=cut_surface(mesh,seams);
    // Gate A F7: user pins with different coordinates on one sector of the surface to
    // unwrap. legacy_010 solves them as the source does (the last pin written wins) and
    // the final pin check refuses with pins_changed; preserve never overwrites a pin, so
    // the run is refused here, before any pinned solve, with that reason and the report.
    if(policy==PinCompletion::preserve) {
        std::map<MeshIndex,Vec2> first;std::set<MeshIndex> conflicts;
        for(const auto& [c,p]:pins) {
            const auto s=state.cs->tables().sector_of_corner[c];const auto [at,inserted]=first.emplace(s,p);
            if(!inserted && at->second!=p)conflicts.insert(s);
        }
        if(!conflicts.empty()) {
            rep["pins_kept"]=false;
            stage({{"stage","unwrap"},{"status","refused"},{"reason","pin_conflict"},{"sectors",std::vector<MeshIndex>(conflicts.begin(),conflicts.end())}});
            timings["unwrap"]=seconds_since(t0);cpu_timings["unwrap"]=cpu_seconds()-cpu0;
            res.accepted=false;res.reason="pins_changed";res.seams=state.seams;res.cut=state.cs;res.exempt=state.exempt;
            timings["total"]=seconds_since(t_all);res.report=std::move(rep_);res.scheduling=batches_->json();
            return res;
        }
    }
    RepairContext ctx;ctx.mesh=mesh;ctx.backend=backend_;ctx.options=&o_;ctx.uv_cache=&uv_cache;ctx.pins=&pins;ctx.iterations=iterations;ctx.pin_completion=policy;
    ctx.scheduling=in_.scheduling;ctx.cancel=control_.cancel;ctx.batches=batches_.get();ctx.recovery_locked=&locked;
    std::uint64_t projected=0;Json deadline_fallback=Json::object();
    if(deadline_.at() || (o_.policy=="repair_v1" && o_.repair.region_recovery_seconds>0)) {
        projected=budget_unwrap(ctx,state,uv_generated,deadline_,deadline_fallback);
        rep["deadline_projected_faces"]=projected;rep["deadline_chart_fallback"]=deadline_fallback;
    } else if(uv_generated && pins.empty() &&
              std::all_of(uv_generated->begin(),uv_generated->end(),[](const Vec2& p) {return std::isfinite(p[0]) && std::isfinite(p[1]);}))
        state.uv=*uv_generated;
    else state.uv=backend_->unwrap(state.cs->tables().cut,std::nullopt,pins); // pins stay where the user put them
    auto m_unwrap=pipeline_metrics(state.cs,state.uv,state.exempt);
    rep["after_unwrap"]=m_unwrap;
    timings["unwrap"]=seconds_since(t0);cpu_timings["unwrap"]=cpu_seconds()-cpu0;
    auto invalid_of=[&]() {
        const auto v=QualityEvaluation(state.geometry(),CornerUvView{state.uv},exempt_policy(state.exempt)).validity();
        std::vector<MeshIndex> bad;
        for(MeshIndex ch=0;ch<state.cs->tables().n_charts;++ch)if(!v.valid[ch] && !state.cs->tables().chart_faces[ch].empty())bad.push_back(ch);
        return bad;
    };
    auto finish=[&](PipelineResult& out) {
        out.seams=state.seams;out.uv=state.uv;out.cut=state.cs;out.exempt=state.exempt;
        timings["total"]=seconds_since(t_all);out.report=std::move(rep_);out.scheduling=batches_->json();
    };
    if(!m_unwrap.at("valid").get<bool>() && !o_.accept_partial) {
        res.accepted=false;res.reason="unwrap_invalid";
        const auto v=QualityEvaluation(state.geometry(),CornerUvView{state.uv},exempt_policy(state.exempt)).validity();
        for(MeshIndex ch=0;ch<state.cs->tables().n_charts;++ch)if(!v.valid[ch])res.invalid_charts.push_back(ch);
        stage({{"stage","unwrap"},{"status","failed"},{"metrics",m_unwrap}});
        finish(res);return res;
    }
    const bool fallback_charts=deadline_fallback.value("conformal_charts",0)!=0 || deadline_fallback.value("harmonic_charts",0)!=0;
    const bool any_fallback=projected || fallback_charts || deadline_fallback.value("local_charts",0)!=0;
    const char* unwrap_status=any_fallback?"timed_out":"done";
    if(o_.policy=="repair_v1" && any_fallback)unwrap_status=deadline_.expired()?"timed_out":"fallback";
    stage({{"stage","unwrap"},{"status",unwrap_status},{"projected_faces",projected},{"metrics",m_unwrap}});
    const auto say_fn=[this](const std::string& msg) {say(msg);};
    auto repairs_summary=[](const Json& records,const char* status) {std::size_t n=0;for(const auto& r:records)n+=r.at("status")==status;return n;};
    // ---------------------------------------------------------------- repair
    if(!m_unwrap.at("valid").get<bool>()) {
        say("repairing invalid islands");
        t0=Clock::now();cpu0=cpu_seconds();
        auto repaired=repair_all(ctx,state,say_fn);
        for(auto& r:repaired.records)rep["repairs"].push_back(r);
        const auto m_rep=pipeline_metrics(state.cs,state.uv,state.exempt);
        stage({{"stage","repair"},{"status",repaired.invalid.empty()?"done":"failed"},{"repaired",repairs_summary(repaired.records,"repaired")},
               {"exempted",repairs_summary(repaired.records,"exempted")},{"unresolved",repairs_summary(repaired.records,"unresolved")},
               {"invalid_charts",repaired.invalid},{"metrics",m_rep}});
        timings["repair"]=seconds_since(t0);cpu_timings["repair"]=cpu_seconds()-cpu0;
        if(!repaired.invalid.empty()) {
            res.accepted=false;res.reason="chart_unrepairable";res.invalid_charts=repaired.invalid;
            finish(res);return res;
        }
    }
    // ------------------------------------------------------- per-chart gate
    if(o_.final_gate && !deadline_.expired()) {
        t0=Clock::now();cpu0=cpu_seconds();
        // enhanced_v1 (T32): the state before the ordinary over-budget split is recorded so the
        // enhancement pass can evaluate relief slits before that split. Read-only: the gate runs
        // exactly as native_baseline and candidate zero is its result.
        std::shared_ptr<const CutSurface> gate_before;std::vector<Vec2> gate_uv;
        if(o_.policy=="enhanced_v1") {gate_before=state.cs;gate_uv=state.uv;}
        auto ginfo=final_chart_gate(ctx,state,say_fn);
        if(gate_before && state.cs!=gate_before) {
            std::vector<MeshIndex> over;std::vector<double> p95;
            for(const auto& row:ginfo.at("over")) {over.push_back(row[0].get<MeshIndex>());p95.push_back(row[1].is_number()?row[1].get<double>():quality_nan);}
            splits_=record_split("final_chart_gate",gate_before,gate_uv,*state.cs,over,p95);
        }
        rep["chart_gate"]=ginfo;rep["chart_over_budget"]=ginfo.at("chart_over_budget");
        const bool area_listed=ginfo.contains("area_over") && !ginfo.at("area_over").empty(); // FIX-04
        if(!ginfo.at("over").empty() || area_listed) {
            // FIX-04: the area-only rows follow the anisotropy rows; they never make a split.
            const auto& listed=ginfo.at("chart_over_budget");Json head=Json::array();
            for(std::size_t i=0;i<std::min(listed.size(),ginfo.at("over").size());++i)head.push_back(listed[i]);
            const bool changed=(o_.policy=="legacy_010"?listed:head)!=ginfo.at("over");
            const bool split=!ginfo.at("split").empty() && ginfo.at("split")[0].at("valid").get<bool>() && changed;
            stage({{"stage","chart_gate"},{"status",split?"split":"listed"},{"over",ginfo.at("over")},{"still_over",ginfo.at("chart_over_budget")}});
        }
        timings["chart_gate"]=seconds_since(t0);cpu_timings["chart_gate"]=cpu_seconds()-cpu0;
    } else rep["chart_over_budget"]=Json::array();
    if(!rep["repairs"].empty() && !deadline_.expired()) {
        t0=Clock::now();cpu0=cpu_seconds();
        auto absorption=absorb_repair_fragments(ctx,state,&evaluate,locked,deadline_);
        for(std::size_t e=0;e<ne;++e)state.seams[e]=state.cs->tables().cut[e] && !m_.tables().edge_boundary[e];
        rep["repair_absorption"]=absorption;
        if(absorption.at("absorbed").get<std::uint64_t>()) {
            const auto p95=QualitySnapshot(state.geometry(),CornerUvView{state.uv},exempt_policy(state.exempt)).chart_anisotropy_p95();
            Json over=Json::array();const double budget=o_.layout.max_anisotropy_p95*o_.final_gate_tolerance;
            for(MeshIndex ch=0;ch<state.cs->tables().n_charts;++ch)if(p95[ch]>budget)over.push_back(Json::array({ch,detail::python_round(p95[ch],3)}));
            if(o_.policy!="legacy_010") {
                // FIX-04: area-only violators of the absorbed surface are listed too.
                std::set<MeshIndex> listed;for(const auto& row:over)listed.insert(row[0].get<MeshIndex>());
                std::vector<MeshIndex> area_only;
                for(const auto& row:area_over_budget(state,o_))if(!listed.contains(row[0].get<MeshIndex>()))area_only.push_back(row[0].get<MeshIndex>());
                for(auto& row:area_only_rows(area_only,p95))over.push_back(std::move(row));
            }
            rep["chart_over_budget"]=over;
        }
        timings["repair_absorption"]=seconds_since(t0);cpu_timings["repair_absorption"]=cpu_seconds()-cpu0;
    }
    m_unwrap=pipeline_metrics(state.cs,state.uv,state.exempt);
    // ------------------------------------------------------- clean contours
    if(o_.clean_contours && !generation_timeout) {
        say("cleaning UV contours");
        t0=Clock::now();cpu0=cpu_seconds();
        ContourRequest request;request.options=o_.boundary;for(const auto& [c,p]:pins)request.pinned_corners.insert(c);
        for(auto a:in_.anchors)request.manual_anchors.insert(static_cast<MeshIndex>(a));
        const auto validator=contour_validator(o_,state.exempt);request.solve=&interior;request.validate=&validator;request.pin_policy=policy;
        request.parallel=core; // GB performance: the native interior solver runs on workers
        auto cleaned=clean_all(state.geometry(),state.uv,request,std::nullopt,limit(0.95),emit_sink);
        const auto m2=pipeline_metrics(state.cs,cleaned.uv,state.exempt);
        auto d=[](const Json& j,const char* k) {return j.at(k).is_null()?std::numeric_limits<double>::quiet_NaN():j.at(k).get<double>();};
        const bool ok=m2.at("valid").get<bool>() &&
            d(m2,"anisotropy_p95")<=std::min(d(m_unwrap,"anisotropy_p95")*o_.p95_growth,std::max(o_.abs_anisotropy_p95,d(m_unwrap,"anisotropy_p95"))) &&
            d(m2,"area_p95")<=std::min(d(m_unwrap,"area_p95")*o_.p95_growth,std::max(o_.abs_area_p95,d(m_unwrap,"area_p95"))) &&
            d(m2,"anisotropy_max")<=d(m_unwrap,"anisotropy_max")*o_.max_growth+1e-9;
        Json st=Json::object();
        for(const auto& r:cleaned.charts)st[r.status]=st.value(r.status,0)+1;
        Json contours=Json::array();
        for(const auto& r:cleaned.charts) {
            Json arcs=Json::array();
            for(const auto& a:r.arcs)arcs.push_back(Json::array({a.mode,a.n,a.jag_before,std::isfinite(a.jag_after)?Json(a.jag_after):Json(nullptr)}));
            contours.push_back({{"chart",r.chart},{"status",r.status},{"reason",r.reason},{"intensity",r.intensity},{"arcs",arcs}});
        }
        rep["contours"]=contours;
        if(ok) {
            state.uv=std::move(cleaned.uv);
            stage({{"stage","contours"},{"status",limit(0.95).expired()?"timed_out":"done"},{"charts",st},{"metrics",m2}});
        } else stage({{"stage","contours"},{"status","reverted"},{"charts",st},{"metrics",m2},{"reason","deformation_budget"}});
        timings["contours"]=seconds_since(t0);cpu_timings["contours"]=cpu_seconds()-cpu0;
    }
    // --------------------------------------------------------------- validate
    auto final_metrics=pipeline_metrics(state.cs,state.uv,state.exempt);
    auto invalid=invalid_of();
    if(!invalid.empty() && o_.accept_partial) {
        // a rung after the repair stage broke a chart: one more pass of the ladder
        auto repaired=repair_all(ctx,state,say_fn);
        for(auto& r:repaired.records)rep["repairs"].push_back(r);
        invalid=repaired.invalid;
        final_metrics=pipeline_metrics(state.cs,state.uv,state.exempt);
    }
    if(projected && deadline_.expired())rep["output_hygiene"]={{"status","skipped"},{"reason","emergency_face_projection"}};
    else if(final_metrics.at("valid").get<bool>() && invalid.empty() && (o_.uniform_density || o_.orient_islands)) {
        t0=Clock::now();cpu0=cpu_seconds();
        auto output_uv=state.uv;Json hygiene=Json::object();
        if(o_.uniform_density) {auto r=pipeline_hooks::normalise_density(state.geometry(),output_uv,pins,state.exempt);output_uv=std::move(r.first);hygiene["density"]=r.second;}
        if(o_.orient_islands) {auto r=pipeline_hooks::orient_islands(state.geometry(),output_uv,pins);output_uv=std::move(r.first);hygiene["orientation"]=r.second;}
        const auto measured=pipeline_metrics(state.cs,output_uv,state.exempt);
        if(measured.at("valid").get<bool>()) {state.uv=std::move(output_uv);final_metrics=measured;hygiene["status"]="done";}
        else hygiene["status"]="reverted";
        rep["output_hygiene"]=hygiene;
        timings["output_hygiene"]=seconds_since(t0);cpu_timings["output_hygiene"]=cpu_seconds()-cpu0;
    }
    rep["final"]=final_metrics;
    // GB M4: enhanced_v1 finishes the island target again on an adopted state, from the generation's record.
    const Json generation_target=o_.policy=="enhanced_v1" && rep.contains("island_target")?rep["island_target"]:Json();
    if(rep.contains("island_target"))pipeline_hooks::finish_target(rep["island_target"],*band,state,o_,final_metrics);
    auto pins_held=[&] {
        bool kept=true;
        for(const auto& [c,p]:pins) {const std::array<double,2> got{state.uv[c][0],state.uv[c][1]},want{p[0],p[1]};kept=kept && detail::numpy_allclose(got,want);}
        return kept;
    };
    bool pins_kept=pins_held();
    rep["pins_kept"]=pins_kept;
    Json passthrough_charts=Json::array();
    if(count(passthrough)) {
        std::set<MeshIndex> charts;for(std::size_t f=0;f<nf;++f)if(passthrough[f])charts.insert(state.cs->tables().chart_of_face[f]);
        passthrough_charts=std::vector<MeshIndex>(charts.begin(),charts.end());
    }
    rep["repairs_summary"]={{"repaired",repairs_summary(rep["repairs"],"repaired")},{"exempted",repairs_summary(rep["repairs"],"exempted")},
        {"unresolved",repairs_summary(rep["repairs"],"unresolved")},{"unrepairable",repairs_summary(rep["repairs"],"unrepairable")},
        {"exempt_faces",count(state.exempt)},{"passthrough_charts",passthrough_charts}};
    bool accepted=false;std::string reason;
    auto decide=[&] {
        accepted=final_metrics.at("valid").get<bool>() && invalid.empty() && pins_kept;reason.clear();
        if(!accepted)reason=!pins_kept?"pins_changed":!invalid.empty() && o_.accept_partial?"chart_unrepairable":"final_invalid";
        else if(!baseline.is_null() && baseline.at("valid").get<bool>() && !o_.generate) {
            // editing existing seams: never degrade beyond the budget
            auto d=[](const Json& j,const char* k) {return j.at(k).is_null()?std::numeric_limits<double>::quiet_NaN():j.at(k).get<double>();};
            if(d(final_metrics,"anisotropy_p95")>std::max(d(baseline,"anisotropy_p95")*o_.p95_growth,o_.abs_anisotropy_p95) ||
               d(final_metrics,"area_p95")>std::max(d(baseline,"area_p95")*o_.p95_growth,o_.abs_area_p95)) {accepted=false;reason="degraded_beyond_budget";}
        }
    };
    decide();
    // ------------------------------------------------ enhanced_v1 (Phase B, T28)
    // Candidate zero is the complete native_baseline state above; no other policy gets here.
    if(o_.policy=="enhanced_v1") {
        const auto outcome=enhance(state,pins,locked,accepted,reason);
        if(outcome.adopted) {
            final_metrics=pipeline_metrics(state.cs,state.uv,state.exempt);rep["final"]=final_metrics;
            invalid=invalid_of();pins_kept=pins_held();rep["pins_kept"]=pins_kept;decide();
            // GB M4: every other field that describes the result's map is measured again on the
            // kept state (budget rows, island target, pass-through charts); candidate zero's
            // values are kept in report.enhanced.baseline_result.
            Json baseline_result={{"chart_over_budget",rep["chart_over_budget"]},{"passthrough_charts",rep["repairs_summary"]["passthrough_charts"]}};
            rep["chart_over_budget"]=retained_over_budget_;
            rep["repairs_summary"]["passthrough_charts"]=passthrough_charts_of(*state.cs,passthrough);
            if(!generation_target.is_null()) {
                baseline_result["island_target"]=rep["island_target"];rep["island_target"]=generation_target;
                pipeline_hooks::finish_target(rep["island_target"],*band,state,o_,final_metrics);
            }
            rep["enhanced"]["baseline_result"]=std::move(baseline_result);
        }
        if(outcome.ran && accepted && !outcome.feasible) {accepted=false;reason="constraint_conflict";}
        else if(outcome.ran && accepted && !outcome.admissible) {accepted=false;reason="enhanced_invalid";}
        rep["enhanced"]["accepted"]=accepted;
    }
    res.accepted=accepted;res.reason=reason;res.invalid_charts=invalid;
    finish(res);
    return res;
}
}
Run::EnhancedOutcome Run::enhance(RepairState& state,const PinMap& pins,const EdgeMask& locked,bool baseline_accepted,const std::string& baseline_reason) {
    EnhancedOutcome outcome;
    if(!baseline_accepted) {
        // Only a complete, accepted candidate zero enters the pass; the refusal stands.
        rep_["enhanced"]={{"policy",enhanced_policy_version},{"status","not_run"},{"reason","baseline_refused"},{"baseline_reason",baseline_reason}};
        return outcome;
    }
    const auto t0=Clock::now();const auto cpu0=cpu_seconds();
    // GB m4: the pass never fails a run native_baseline accepts. A failure anywhere in it (measuring
    // candidate zero, a report, the post-pass necessity or envelope domain, an error a stage did not
    // isolate) keeps native_baseline's decision and result, with the failure in report.enhanced;
    // only a cancellation propagates (a cancelled run applies nothing).
    const RepairState original=state;const auto stages=rep_["stages"].size();
    try {return enhance_pass(state,pins,locked,t0,cpu0);}
    catch(const std::exception& error) {
        if(const auto* contract=dynamic_cast<const ContractError*>(&error);contract && contract->code==ErrorCode::cancelled)throw;
        throw_if_cancelled(control_.cancel);
        state=original;
        while(rep_["stages"].size()>stages)rep_["stages"].erase(rep_["stages"].size()-1);
        retained_over_budget_=Json::array();
        rep_["enhanced"]={{"policy",enhanced_policy_version},{"status","failed"},{"reason","pass_failure"},{"failure",failure_record(error)},
                          {"baseline_retained",true}};
        stage({{"stage","enhanced"},{"status","failed"},{"failure",rep_["enhanced"]["failure"]}});
        rep_["timings"]["enhanced"]=seconds_since(t0);rep_["cpu_timings"]["enhanced"]=cpu_seconds()-cpu0;
        return EnhancedOutcome{};
    }
}
Run::EnhancedOutcome Run::enhance_pass(RepairState& state,const PinMap& pins,const EdgeMask& locked,Clock::time_point t0,double cpu0) {
    EnhancedOutcome outcome;
    say("enhancing seam layout");
    const auto& t=m_.tables();const auto ne=t.edge_vertices.size();
    ProposalContext ctx;ctx.mesh=in_.mesh;ctx.options=o_;ctx.exempt=state.exempt;ctx.cancel=control_.cancel;ctx.deadline=deadline_;
    ctx.float32_storage=backend_->float32_storage(); // GB: the host's float32 UV storage (strict pin check)
    if(!prefs_)prefs_=std::make_shared<ResolvedPreferences>(resolve_preferences(m_,in_.constraints,o_.enhanced.preferences,in_.locked,in_.anchors,
                                                                                    o_.layout.feature_angle,o_.layout.hide_dir));
    ctx.constraints=hard_constraints(m_,locked,in_.anchors,pins,prefs_->required,prefs_->forbidden);
    // Paint changes invalidate proposals (a new revision), never geometry or topology caches.
    ctx.constraints.revision=std::stoull(prefs_->digest.substr(0,15),nullptr,16);
    Json exposure={{"measured",false}};
    if(o_.layout.exposure_K>0) {
        // Raw measured exposure for the objective's V term and the routing heuristic; paint
        // and discounts never modify it.
        try {
            const auto field=exposure_field(m_,static_cast<std::uint64_t>(o_.layout.exposure_K));
            ctx.raw_exposure=field.edge;
            // GB (M3, seam_objective.v2): the sampling uncertainty of the same measured field.
            ctx.exposure_uncertainty=std::make_shared<const ExposureUncertainty>(in_.mesh,field);
            exposure={{"measured",true},{"directions",o_.layout.exposure_K},{"digest",values_digest(ctx.raw_exposure)},
                      {"uncertainty",ctx.exposure_uncertainty->to_json()}};
        } catch(const ContractError& error) {if(error.code==ErrorCode::cancelled)throw;exposure={{"measured",false},{"error",error.what()}};}
    }
    prefs_->raw_exposure=ctx.raw_exposure;
    ctx.preferences=prefs_; // generators route and price seams with the run's paint and viewer (T31, T32, T34+)
    const auto iterations=backend_->iterations()?backend_->iterations():10;
    // GB (m3): the backend's unwrap has one numerical contract, its default solve; a proposal
    // naming another initialization is refused before any work (supports_initialization).
    ctx.supports_initialization=[](const std::string& initialization) {return initialization.empty();};
    ctx.solve=[this,iterations](const std::shared_ptr<const CutSurface>& cs,std::span<const MeshIndex> charts,const PinMap& p,const std::string& initialization) {
        if(!initialization.empty())
            throw ContractError(ErrorCode::invalid_input,"enhanced.initialization","the backend solver honours only its default initialization, not '"+initialization+"'");
        return backend_->unwrap(cs->tables().cut,charts,p,iterations);
    };
    const auto origin=begin_enhancement(ctx,state.cs,state.uv,{{"source","native_baseline"}});
    // T32: relief is scheduled before the ordinary over-budget split (Section 21), so for the
    // relief proposal that welds a recorded split the run origin of the faces that split cut is
    // the chart before it; the absolute budgets and the origin's chart keys and map are
    // unchanged. GB m10: staged, not applied: every other generator's trial is measured against
    // candidate zero's split pieces until a state adopts that reunion.
    Json split_report=Json::array();
    for(const auto& record:splits_) {
        const std::array<MeshIndex,1> one{record.chart};
        const QualitySnapshot before(quality_geometry(record.surface),CornerUvView{record.uv},exempt_policy(ctx.exempt),std::span<const MeshIndex>(one));
        ctx.envelope.stage_rebase(record.faces,distortion_samples(before)); // GB m10: relief trials only
        split_report.push_back(record.summary());
    }
    ctx.baseline_splits=splits_;
    // T33: one evidence cache for the pass's spur-removal rounds and the final classification.
    ctx.necessity=std::make_shared<SeamNecessityCache>();
    // Section 21: constraint repair first (only when painted hard constraints exist), then
    // the built-in operations, then any caller's generators. T40: one bounded search over the
    // whole schedule (estimate-ordered, deduplicated across the pass, per-chart and per-feature
    // bounds, proven bounds before solves); the baseline pathways all ran in candidate zero.
    std::vector<ProposalGenerator> generators;
    // GB M5: enhanced.constraint_repair.enabled switches the repair off for the painted ablations.
    if(prefs_->hard_constraints() && o_.enhanced.constraint_repair.enabled)generators.push_back(constraint_repair_generator(prefs_,o_.enhanced.preferences));
    for(auto& g:enhanced_generators(o_))generators.push_back(std::move(g));
    for(const auto& g:in_.proposal_generators)generators.push_back(g);
    ctx.search=std::make_shared<SearchLedger>();
    auto result=run_enhancement(origin,generators,ctx);
    // A cancellation stops the pass with its incumbent intact; the run's own contract then
    // holds (cancelled runs apply nothing). A deadline keeps the incumbent and the run goes on.
    throw_if_cancelled(control_.cancel);
    auto kept=result.state;auto report=result.report;report["status"]="done";report["raw_exposure"]=exposure;
    report["schedule"]=enhanced_schedule(o_,prefs_->hard_constraints());
    report["baseline_splits"]=std::move(split_report);
    // Why edges are forbidden, required or attractive: every term of the routing cost of the
    // painted hard edges and of the retained seams the pass added (bounded).
    report["preferences"]=prefs_->summary();
    {
        Json explained=Json::array();std::set<MeshIndex> listed;
        auto add=[&](MeshIndex e) {if(explained.size()<256 && listed.insert(e).second)explained.push_back(prefs_->explain(m_,e,o_.enhanced.preferences));};
        for(std::size_t e=0;e<ne;++e)if(prefs_->required[e] || prefs_->forbidden[e])add(static_cast<MeshIndex>(e));
        for(std::size_t e=0;e<ne;++e)if(kept.cuts[e] && !origin.cuts[e])add(static_cast<MeshIndex>(e));
        report["preferences"]["edges"]=std::move(explained);
    }
    outcome.ran=true;outcome.adopted=kept.revision!=origin.revision;
    // GB (m4): a map writer after the pass that fails keeps the state it was given; the failure is
    // recorded with the pass's (cancellation and the caller's contract errors propagate).
    auto writer_failure=[&](const char* stage,const std::exception& error,Json& record) {
        if(!isolated_failure(error))throw;
        auto entry=failure_record(error);record["failed"]=entry;entry["stage"]="writer";entry["generator"]=stage;
        if(!ctx.diagnostics.contains("failures"))ctx.diagnostics["failures"]=Json::array();
        ctx.diagnostics["failures"].push_back(std::move(entry));
    };
    // T40, Section 21 "contour cleanup and ordinary repair/gate": the charts the pass changed
    // get the baseline's UV contour cleanup (T15, same request and validator as the baseline
    // stage) as a map writer through the transaction, kept only when admissible and ranking no
    // worse than the uncleaned state. The repair ladder and the chart gate need no second run:
    // a trial is never repaired or exempted (it is refused), and every trial already passed the
    // gate's P95 budget rule with the P99/maximum caps.
    // GB N5: a deadline that expired before it is a time-budget outcome of the pass (recorded).
    auto deadline_cut=[&] {if(ctx.search->stopped.empty())ctx.search->stopped="deadline";};
    if(outcome.adopted && o_.clean_contours && ctx.deadline.expired()) {
        report["contour_cleanup"]={{"charts",ctx.envelope.domain_charts(*kept.surface,kept.uv).size()},{"skipped","deadline"},{"committed",false}};
        deadline_cut();
    } else if(outcome.adopted && o_.clean_contours) {
        const auto charts=ctx.envelope.domain_charts(*kept.surface,kept.uv);
        Json record={{"charts",charts.size()}};
        if(!charts.empty())try {
            ContourRequest request;request.options=o_.boundary;for(const auto& [c,p]:pins)request.pinned_corners.insert(c);
            for(auto a:in_.anchors)request.manual_anchors.insert(static_cast<MeshIndex>(a));
            const auto interior=core_?core_->interior_solver():backend_interior_solver(in_.mesh,*backend_);
            const auto validator=contour_validator(o_,state.exempt);request.solve=&interior;request.validate=&validator;request.pin_policy=pin_policy(o_);
            auto cleaned=clean_all(quality_geometry(kept.surface),kept.uv,request,std::span<const MeshIndex>(charts),ctx.deadline);
            Json statuses=Json::object();for(const auto& r:cleaned.charts)statuses[r.status]=statuses.value(r.status,0)+1;
            record["statuses"]=statuses;
            if(cleaned.uv!=kept.uv) {
                Json trial;const auto before=kept.revision;
                kept=rewrite_state_uv(kept,std::move(cleaned.uv),"contour_cleanup",ctx,&trial,true);
                record["committed"]=kept.revision!=before;record["trial"]={{"admissible",trial.value("admissible",false)},{"rejection",trial.value("rejection","")},
                                                                    {"budget",trial.value("budget",Json::object())},{"changed_charts",trial.value("changed_charts",Json::array())}};
            } else record["committed"]=false;
        } catch(const std::exception& error) {writer_failure("contour_cleanup",error,record);}
        report["contour_cleanup"]=std::move(record);report["retained"]=run_state_summary(kept);report["provenance"]=kept.provenance;
    }
    if(outcome.adopted && (o_.uniform_density || o_.orient_islands)) {
        // Output hygiene of the committed state is a map writer: measured again through the
        // transaction and kept only when the rewritten state is still admissible.
        Json record;
        try {
        const auto geometry=quality_geometry(kept.surface);auto uv=kept.uv;
        if(o_.uniform_density)uv=pipeline_hooks::normalise_density(geometry,uv,pins,state.exempt).first;
        if(o_.orient_islands)uv=pipeline_hooks::orient_islands(geometry,uv,pins).first;
        // GB (Section 11.3 policy): a chart with an unrepaired strict defect is never written
        // over except by a validated trial; the writer keeps its map (candidate zero's, which the
        // baseline's own hygiene already normalised and oriented).
        Json held=Json::array();
        if(!kept.strict_defects.empty()) {
            const auto& kt=kept.surface->tables();std::set<MeshIndex> hold(kept.strict_defects.begin(),kept.strict_defects.end());
            for(std::size_t c=0;c<uv.size();++c)if(hold.contains(kt.chart_of_corner[c]))uv[c]=kept.uv[c];
            held=kept.strict_defects;
        }
        kept=rewrite_state_uv(kept,std::move(uv),"output_hygiene",ctx,&record);
        if(!held.empty())record["held_strict_defect_charts"]=held;
        } catch(const std::exception& error) {writer_failure("output_hygiene",error,record);}
        report["output_hygiene"]=record;
        report["retained"]=run_state_summary(kept);report["provenance"]=kept.provenance;
    }
    if(outcome.adopted) {
        state.cs=kept.surface;state.uv=kept.uv;
        for(std::size_t e=0;e<ne;++e)state.seams[e]=kept.cuts[e] && !t.edge_boundary[e];
    }
    // Effective (committed), solved (surface) and reported (seams) cuts are one set; the
    // surface always cuts boundary and nonmanifold edges, which the seam output never lists.
    bool consistent=state.cs->tables().cut==kept.cuts;
    for(std::size_t e=0;e<ne;++e) {
        const bool fact=t.edge_boundary[e] || !t.edge_manifold[e] || !t.edge_consistent[e] || t.edge_corners[e].size()!=2;
        consistent=consistent && (fact || (state.seams[e]!=0)==(kept.cuts[e]!=0));
    }
    outcome.feasible=kept.feasible;outcome.valid=kept.valid && consistent;outcome.admissible=kept.admissible() && consistent;
    report["cuts_consistent"]=consistent;report["feasible"]=kept.feasible;report["valid"]=kept.valid && consistent;
    report["admissible"]=outcome.admissible;
    // GB (Section 11.3 policy): candidate zero's strict-only defects, repaired or explicitly
    // reported unrepaired (then `valid` is false, the decision stays native_baseline's; GB B1: an
    // accepted outcome, named in the report and its messages, never an error diagnostic),
    // measured on the state the run returns: overlaps (strict_overlap) and GB N2 pins (pin_check).
    {
        Json stage=nullptr;
        for(const auto& round:report.value("rounds",Json::array()))if(round.value("generator","")=="strict_overlap_repair")stage=round;
        report["strict_overlap"]=strict_overlap_report(ctx,origin,kept,stage);
        report["pin_check"]=pin_check_report(ctx,origin,kept);
    }
    // T33: the measured necessity of every dead end of the reported state (trial removals
    // through the transaction, never committed); the verdict reads the classification.
    {
        const auto necessity=measure_seam_necessity(kept,ctx,*ctx.necessity);
        // GB N5: ends the deadline left unmeasured make the pass a time-budget outcome.
        if(std::any_of(necessity.arcs.begin(),necessity.arcs.end(),[](const SeamEndEvidence& a) {
               return a.classification==SeamEndClass::unmeasured && a.reason=="deadline";}))deadline_cut();
        report["seam_necessity"]=necessity.to_json();
    }
    // GB (m4): every isolated failure of the pass, its writers and the necessity trials.
    report["failures"]=ctx.diagnostics.value("failures",Json::array());
    if(report.contains("diagnostics") && report["diagnostics"].is_object() && ctx.diagnostics.contains("failures"))
        report["diagnostics"]["failures"]=ctx.diagnostics["failures"];
    // T38: the explicit symmetry policy of the pass (absent in the asymmetric mode).
    if(auto symmetry=symmetry_summary(ctx,origin,kept);!symmetry.is_null())report["symmetry"]=std::move(symmetry);
    report["counters"]=ctx.counters.to_json();
    // GB N3: the result's chart_over_budget rows of an adopted state by native_baseline's own rule
    // (read by execute() when the state is adopted); the enhanced budget stays in `retained`.
    if(outcome.adopted) {
        std::string source="not_judged";
        if(rep_.contains("repair_absorption") && rep_["repair_absorption"].is_object() &&
           rep_["repair_absorption"].value("absorbed",std::uint64_t{0})>0)source="repair_absorption";
        else if(rep_.contains("chart_gate"))source="chart_gate";
        Json rule;retained_over_budget_=native_rule_over_budget_rows(kept,origin,rep_["chart_over_budget"],source,o_,&rule);
        report["chart_over_budget_rule"]=std::move(rule);
    }
    // GB N5: every cut the deadline made (a round, a generator, the contour cleanup, the
    // classification) is the pass's `stopped`.
    const auto stopped=ctx.search->stopped;
    report["stopped"]=stopped;
    if(report.contains("search") && report["search"].is_object())report["search"]["stopped"]=stopped;
    rep_["enhanced"]=std::move(report);
    Json record={{"stage","enhanced"},{"status",outcome.adopted?"committed":"baseline_retained"},{"proposals",ctx.counters.proposals},
                 {"commits",ctx.counters.commits},{"feasible",kept.feasible},{"valid",outcome.valid},
                 {"strict_overlap",rep_["enhanced"]["strict_overlap"]["status"]}};
    if(!stopped.empty())record["stopped"]=stopped;
    stage(std::move(record));
    rep_["timings"]["enhanced"]=seconds_since(t0);rep_["cpu_timings"]["enhanced"]=cpu_seconds()-cpu0;
    return outcome;
}
PipelineResult run_pipeline_with_services(const PipelineInput& input,UnwrapBackend* backend,const ProgressSink& progress,const Control& control,const LayoutPostprocessServices& services) {
    if(!input.mesh)throw ContractError(ErrorCode::invalid_input,"mesh","missing mesh snapshot");
    // Gate A F3: the run's allowance reaches every table, solve and exposure field it
    // builds (worker tasks inherit it), and the simultaneous tables are checked first.
    MemoryBudgetScope budget(input.memory_budget);
    check_run_memory(table_counts(*input.mesh),input.memory_budget);
    const auto cpu=cpu_seconds();const auto wall=Clock::now();
    const auto scheduling=resolve_scheduling(input.scheduling);
    // A run owns a bounded cache unless its caller (a session) supplies one.
    std::optional<TopologyCache> cache;if(!active_cache() && input.cache_bytes)cache.emplace(input.cache_bytes);
    CacheScope scope(cache?&*cache:active_cache());
    // T24: the run owns its native workers (joined when the run returns or throws) unless
    // the caller has a pool active. One thread is the reference serial schedule.
    std::optional<WorkerPool> pool;if(!active_workers() && scheduling.threads>1)pool.emplace(scheduling.threads);
    WorkerScope workers(pool?&*pool:active_workers());
    check_layout_operation(input);
    auto result=input.operation=="postprocess_layout"?postprocess_layout(input,progress,control,services):Run(input,backend,progress,control).execute();
    // enhanced_v1: a run that stopped before a complete map never entered the pass.
    if(input.options.policy=="enhanced_v1" && !result.report.contains("enhanced"))
        result.report["enhanced"]={{"policy",enhanced_policy_version},{"status","not_run"},{"reason","baseline_incomplete"},{"baseline_reason",result.reason}};
    if(const auto* used=active_workers()) {
        const auto st=used->stats();
        result.scheduling["workers"]={{"threads",used->threads()},{"parallel_calls",st.parallel_calls},{"inline_calls",st.inline_calls},
                                      {"tasks",st.tasks},{"skipped",st.skipped},{"max_width",st.max_width},{"nested_solver_threads",st.nested_solver_threads}};
    } else result.scheduling["workers"]={{"threads",1},{"parallel_calls",0},{"nested_solver_threads",1}};
    result.scheduling["threads"]=scheduling.threads;result.scheduling["batch_max_charts"]=scheduling.batch_max_charts;
    result.scheduling["batch_max_triangles"]=scheduling.batch_max_triangles;
    auto& rep=result.report;
    rep["cpu_s"]=cpu_seconds()-cpu;
    double spent=0;for(const auto& [k,v]:rep["cpu_timings"].items())spent+=v.get<double>();
    rep["cpu_timings"]["other"]=std::max(0.0,rep["cpu_s"].get<double>()-spent);
    rep["wall_s"]=seconds_since(wall);
    // GB m15: a deadline that stopped the enhanced pass (its stage record's `stopped`) is a
    // time-budget outcome as well: the kept layout is the best validated one the pass reached.
    bool timed_out=false;
    for(const auto& s:rep["stages"])timed_out=timed_out || s.value("status","")=="timed_out" || s.value("stopped","")=="deadline";
    rep["timed_out"]=timed_out;
    return result;
}
PipelineResult run_pipeline(const PipelineInput& input,UnwrapBackend* backend,const ProgressSink& progress,const Control& control) {
    return run_pipeline_with_services(input,backend,progress,control,LayoutPostprocessServices{});
}

}
