// SPDX-License-Identifier: GPL-3.0-or-later
// Phase B thin-shell proposals (T37, Section 18): area-weighted PCA nominates thin
// components; opposite-side correspondence (inward rays on the geometry-revisioned BVH),
// front/back labels from the thin-axis normal alignment and the position relative to the
// locally estimated middle surface, smoothing under hard barriers and geometric speckle
// repair confirm them; the front/back border (real edges) and a visibility-biased nearby
// border are proposed as compound cuts through the ordinary transaction. There is no early
// acceptance: every piece is solved and validated like any other move. enhanced_v1 only.
#include "cleanseams/proposals.hpp"
#include "cleanseams/physical_scale.hpp"
#include "cleanseams/geometry_bvh.hpp"
#include "cleanseams/cache.hpp"
#include "cleanseams/errors.hpp"
#include <Eigen/Dense>
#include <algorithm>
#include <cmath>
#include <map>
#include <numeric>
#include <set>
namespace cleanseams {
namespace {
using Json=nlohmann::json;
constexpr const char* name="shells";
constexpr const char* version="1";
bool always_cut(const MeshTables& t,std::size_t e) {
    return t.edge_boundary[e] || !t.edge_manifold[e] || !t.edge_consistent[e] || t.edge_corners[e].size()!=2;
}
Vec3 add(const Vec3& a,const Vec3& b) {return {a[0]+b[0],a[1]+b[1],a[2]+b[2]};}
Vec3 sub(const Vec3& a,const Vec3& b) {return {a[0]-b[0],a[1]-b[1],a[2]-b[2]};}
Vec3 mul(const Vec3& a,double s) {return {a[0]*s,a[1]*s,a[2]*s};}
double dot(const Vec3& a,const Vec3& b) {return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];}
double median_weighted(std::vector<std::pair<double,double>> v) {
    if(v.empty())return 0.0;
    std::sort(v.begin(),v.end());double total=0;for(const auto& x:v)total+=x.second;
    double run=0;for(const auto& x:v) {run+=x.second;if(run>=0.5*total)return x.first;}
    return v.back().first;
}
// A component's sorted faces and the slot of a face in them (-1 outside): GB m7, every array of
// the evidence is local to the component.
struct Local {
    std::span<const MeshIndex> faces;
    std::ptrdiff_t operator()(MeshIndex f) const {const auto it=std::lower_bound(faces.begin(),faces.end(),f);return it!=faces.end() && *it==f?it-faces.begin():-1;}
};
// Front/back border of `label` (per slot) over the component: edges whose two faces are in it
// and labelled differently.
std::vector<MeshIndex> border(const MeshSnapshot& m,const Local& local,const std::vector<std::int8_t>& label) {
    const auto& t=m.tables();std::set<MeshIndex> out;
    auto at=[&](MeshIndex f) -> std::int8_t {const auto i=local(f);return i<0?std::int8_t{-1}:label[static_cast<std::size_t>(i)];};
    for(auto f:local.faces)for(auto c:m.face_corners(f)) {
        const auto e=t.corner_edge[c];const auto pair=t.edge_faces[e];
        if(pair.size()!=2 || at(pair[0])<0 || at(pair[1])<0)continue;
        if(at(pair[0])!=at(pair[1]))out.insert(e);
    }
    return {out.begin(),out.end()};
}
// Every component of `edges` is a simple closed loop; the loop count.
std::pair<bool,std::uint64_t> simple_loops(const MeshSnapshot& m,std::span<const MeshIndex> edges) {
    const auto& t=m.tables();std::map<MeshIndex,std::vector<MeshIndex>> at;
    for(auto e:edges)for(auto v:t.edge_vertices[e])at[v].push_back(e);
    std::set<MeshIndex> left(edges.begin(),edges.end());std::uint64_t loops=0;bool simple=!edges.empty();
    while(!left.empty()) {
        std::vector<MeshIndex> comp,stack{*left.begin()};left.erase(left.begin());
        while(!stack.empty()) {const auto e=stack.back();stack.pop_back();comp.push_back(e);for(auto v:t.edge_vertices[e])for(auto f:at[v])if(left.erase(f))stack.push_back(f);}
        std::sort(comp.begin(),comp.end());++loops;
        simple=simple && simple_loop_order(m,comp).size()==comp.size();
    }
    return {simple,loops};
}
}
nlohmann::json ShellEvidence::to_json() const {
    return {{"component",component},{"faces",faces},{"area",area},{"extents",extents},{"thin_axis",thin_axis},{"thin_ratio",thin_ratio},{"pca_ratio",pca_ratio},
            {"nominated",nominated},{"front_fraction",front_fraction},{"back_fraction",back_fraction},{"edge_on_fraction",edge_on_fraction},
            {"correspondence",correspondence},{"correspondence_samples",correspondence_samples},{"thickness",thickness},
            {"speckles_repaired",speckles_repaired},{"rim_edges",rim.size()},{"rim_loops",rim_loops},{"rim_simple",rim_simple},
            {"biased_rim_edges",biased_rim.size()},{"confidence",confidence},{"status",status},{"reasons",reasons},{"scale",scale}};
}
std::vector<ShellEvidence> shell_evidence(const MeshSnapshot& m,const TriangleBvh& bvh,const HardConstraints& hc,std::span<const double> routing,
                                          const PipelineOptions& options,const PhysicalScale* scale,DetectorPass* pass) {
    const auto& t=m.tables();const auto& g=m.geometry();const auto& o=options.enhanced.shells;
    if(!bvh.serves(m))throw ContractError(ErrorCode::stale_state,"shells.bvh","the BVH belongs to another geometry revision");
    const auto nf=g.face_area.size();
    DetectorPass local_pass;DetectorPass& work=pass?*pass:local_pass;
    std::vector<std::vector<MeshIndex>> components(g.n_components);
    for(std::size_t f=0;f<nf;++f)if(g.face_area[f]>0 && g.face_degenerate_kind[f]==0)components[g.face_component[f]].push_back(static_cast<MeshIndex>(f));
    auto hard=[&](std::size_t e) {return (e<hc.locked.size() && hc.locked[e]) || (e<hc.required.size() && hc.required[e]);};
    std::vector<ShellEvidence> out;
    for(MeshIndex k=0;k<g.n_components;++k) {
        if(work.stop())break;
        const auto& faces=components[k];
        const Local local{faces};const auto nc=faces.size();
        auto slot=[&](MeshIndex f) {return static_cast<std::size_t>(local(f));};
        ShellEvidence ev;ev.component=k;ev.faces=faces.size();
        // T39: the face floor is an area on the component's reference chart under the physical
        // scale; label smoothing and the nearby rim's offset get physical radii.
        std::int64_t smoothing_passes=o.smoothing_passes,rim_offset=1;
        if(scale && !faces.empty()) {
            double area=0;for(auto f:faces)area+=g.face_area[f];
            const auto floor=scale->floor(o.min_faces,faces);const auto smooth=scale->diffusion(o.smoothing_passes,faces);
            const auto offset=scale->rings(1,faces);
            smoothing_passes=smooth.steps;rim_offset=std::max<std::int64_t>(1,offset.steps);
            ev.scale={{"face_floor",floor.to_json()},{"smoothing",smooth.to_json()},{"rim_offset",offset.to_json()}};
            if(floor.below(area,faces.size())) {ev.status="rejected";ev.reasons={"too_few_faces"};out.push_back(std::move(ev));continue;}
        }
        else if(static_cast<std::int64_t>(faces.size())<o.min_faces) {ev.status="rejected";ev.reasons={"too_few_faces"};out.push_back(std::move(ev));continue;}
        // Area-weighted PCA of the face centroids.
        Vec3 mean{0,0,0};for(auto f:faces) {ev.area+=g.face_area[f];mean=add(mean,mul(g.face_center[f],g.face_area[f]));}
        mean=mul(mean,1/ev.area);
        Eigen::Matrix3d cov=Eigen::Matrix3d::Zero();
        for(auto f:faces) {const auto d=sub(g.face_center[f],mean);const Eigen::Vector3d v(d[0],d[1],d[2]);cov+=g.face_area[f]*v*v.transpose();}
        cov/=ev.area;
        const Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eig(cov);
        // Principal full extents: the range of the component's vertices along each axis.
        std::set<MeshIndex> vertices;for(auto f:faces)for(auto c:m.face_corners(f))vertices.insert(static_cast<MeshIndex>(m.input().corner_vertices[c]));
        for(int i=0;i<3;++i) {
            const Eigen::Vector3d a=eig.eigenvectors().col(i);double lo=INFINITY,hi=-INFINITY;
            for(auto v:vertices) {const auto& p=m.input().positions[v];const double x=a[0]*p[0]+a[1]*p[1]+a[2]*p[2];lo=std::min(lo,x);hi=std::max(hi,x);}
            ev.extents[static_cast<std::size_t>(i)]=hi-lo;
        }
        const double spread0=std::sqrt(std::max(0.0,eig.eigenvalues()[0])),spread1=std::sqrt(std::max(0.0,eig.eigenvalues()[1]));
        ev.pca_ratio=spread1>0?spread0/spread1:1.0;
        Eigen::Vector3d axis=eig.eigenvectors().col(0);
        int big=0;for(int i=1;i<3;++i)if(std::abs(axis[i])>std::abs(axis[big]))big=i;
        if(axis[big]<0)axis=-axis; // deterministic sign
        const Vec3 axis_v{axis[0],axis[1],axis[2]};ev.thin_axis=axis_v;
        // Thin-axis normal alignment.
        std::vector<double> s(nc,0.0);
        for(std::size_t i=0;i<nc;++i) {
            const auto f=faces[i];s[i]=dot(g.face_normal[f],axis_v);
            (s[i]>0.5?ev.front_fraction:s[i]<-0.5?ev.back_fraction:ev.edge_on_fraction)+=g.face_area[f]/ev.area;
        }
        // Opposite-side correspondence: a ray from a front or back face into the solid (-n)
        // must meet an opposite-facing face of the same component within a few thicknesses.
        const double reach=std::max(4.0*spread0,1e-9*std::max(ev.extents[2],1.0));
        std::vector<MeshIndex> sampled;for(std::size_t i=0;i<nc;++i)if(std::abs(s[i])>0.5)sampled.push_back(faces[i]);
        if(static_cast<std::int64_t>(sampled.size())>o.max_samples) {
            std::vector<MeshIndex> kept;const auto n=sampled.size(),cap=static_cast<std::size_t>(o.max_samples);
            for(std::size_t i=0;i<cap;++i)kept.push_back(sampled[i*n/cap]);
            sampled=std::move(kept);
        }
        std::map<MeshIndex,double> mid;std::vector<std::pair<double,double>> widths;double consistent=0,total=0;
        for(auto f:sampled) {
            const auto n=g.face_normal[f];const auto origin=sub(g.face_center[f],mul(n,1e-9*reach));
            const auto hit=bvh.ray(origin,mul(n,-1.0),reach,[&](MeshIndex h) {return h!=f && g.face_component[h]==k;});
            total+=g.face_area[f];++ev.correspondence_samples;
            if(hit.hit && dot(g.face_normal[hit.face],n)< -0.5) {
                consistent+=g.face_area[f];widths.push_back({hit.t,g.face_area[f]});
                const auto middle=sub(g.face_center[f],mul(n,0.5*hit.t));mid.emplace(f,dot(middle,axis_v));
            }
        }
        ev.correspondence=total>0?consistent/total:0.0;ev.thickness=median_weighted(widths);
        // The nominating thickness ratio: the local opposite-side thickness over the middle
        // principal extent (a bent paddle stays thin; global PCA alone would count its bend),
        // or the smallest extent when no opposite side answers.
        const double thick=ev.thickness>0?ev.thickness:ev.extents[0];
        ev.thin_ratio=ev.extents[1]>0?thick/ev.extents[1]:1.0;
        ev.nominated=ev.thin_ratio<o.thin_ratio;
        // Middle-surface height propagated from the corresponding faces over the face graph.
        // (GB m7: component-local slots; the adjacency across non-hard edges in corner order.)
        std::vector<std::vector<std::pair<std::size_t,MeshIndex>>> adjacent(nc);
        for(std::size_t i=0;i<nc;++i)for(auto c:m.face_corners(faces[i])) {
            const auto e=t.corner_edge[c];if(hard(e))continue;
            for(auto h:t.edge_faces[e])if(h!=faces[i]) {const auto j=local(h);if(j>=0)adjacent[i].emplace_back(static_cast<std::size_t>(j),e);}
        }
        std::vector<double> middle(nc,dot(mean,axis_v));std::vector<std::uint8_t> known(nc,0);
        std::vector<std::size_t> frontier;for(const auto& [f,x]:mid) {const auto i=slot(f);middle[i]=x;known[i]=1;frontier.push_back(i);}
        while(!frontier.empty()) {
            std::vector<std::size_t> next;
            for(auto i:frontier)for(const auto& [j,e]:adjacent[i])if(!known[j]) {known[j]=1;middle[j]=middle[i];next.push_back(j);}
            std::sort(next.begin(),next.end());frontier=std::move(next);
        }
        const double half=std::max(0.5*thick,1e-9*std::max(ev.extents[2],1e-300));
        std::vector<double> value(nc,0.0);
        for(std::size_t i=0;i<nc;++i)value[i]=0.5*s[i]+std::clamp((dot(g.face_center[faces[i]],axis_v)-middle[i])/half,-1.0,1.0);
        // Smoothing that never crosses a hard barrier.
        bool stopped=false;
        for(std::int64_t smooth=0;smooth<smoothing_passes;++smooth) {
            if(work.stop()) {stopped=true;break;}
            auto next=value;
            for(std::size_t i=0;i<nc;++i) {
                double sum=0;std::uint64_t count=0;for(const auto& [j,e]:adjacent[i]) {sum+=value[j];++count;}
                if(count)next[i]=0.5*value[i]+0.5*sum/static_cast<double>(count);
            }
            value=std::move(next);
        }
        if(stopped)break;
        std::vector<std::int8_t> label(nc,-1);for(std::size_t i=0;i<nc;++i)label[i]=value[i]>=0?1:0;
        // Speckle repair: a label region below speckle_fraction of the area joins the label it
        // shares the most border length with (geometry, not an exemption).
        for(int repair=0;repair<4;++repair) {
            bool changed=false;std::vector<std::uint8_t> seen(nc,0);
            for(std::size_t i0=0;i0<nc;++i0) {
                if(seen[i0])continue;
                std::vector<std::size_t> comp,stack{i0};seen[i0]=1;double area=0;
                while(!stack.empty()) {const auto i=stack.back();stack.pop_back();comp.push_back(i);area+=g.face_area[faces[i]];
                    for(const auto& [j,e]:adjacent[i])if(!seen[j] && label[j]==label[i0]) {seen[j]=1;stack.push_back(j);}}
                if(area>=o.speckle_fraction*ev.area)continue;
                double other=0;for(auto i:comp)for(const auto& [j,e]:adjacent[i])if(label[j]!=label[i0])other+=g.edge_length[e];
                if(other<=0)continue;
                for(auto i:comp)label[i]=static_cast<std::int8_t>(1-label[i0]);
                ++ev.speckles_repaired;changed=true;
            }
            if(!changed)break;
        }
        ev.labels.assign(label.begin(),label.end());
        ev.rim=border(m,local,label);
        const auto [simple,loops]=simple_loops(m,ev.rim);ev.rim_simple=simple;ev.rim_loops=loops;
        bool forbidden=false;for(auto e:ev.rim)forbidden=forbidden || (e<hc.forbidden.size() && hc.forbidden[e] && !always_cut(t,e));
        // A visibility-biased nearby rim: the front grown or shrunk by one face ring, when its
        // mean routing multiplier (raw exposure, explicit view, paint) is lower than the raw rim's.
        auto mean_cost=[&](std::span<const MeshIndex> edges) {
            double c=0,l=0;for(auto e:edges) {c+=routing.empty()?g.edge_length[e]:routing[e];l+=g.edge_length[e];}return l>0?c/l:0.0;
        };
        if(!routing.empty() && simple) {
            double best=mean_cost(ev.rim);
            for(int grow=0;grow<2;++grow) {
                auto shifted=label;const std::int8_t from=static_cast<std::int8_t>(grow?0:1),to=static_cast<std::int8_t>(grow?1:0);
                // One face ring per pass across the current border (rim_offset passes).
                auto front=ev.rim;
                for(std::int64_t ring=0;ring<rim_offset && !front.empty();++ring) {
                    auto next=shifted;
                    for(auto e:front) {const auto pair=t.edge_faces[e];for(auto f:pair) {const auto i=local(f);if(i>=0 && shifted[static_cast<std::size_t>(i)]==from)next[static_cast<std::size_t>(i)]=to;}}
                    shifted=std::move(next);front=border(m,local,shifted);
                }
                const auto alt=border(m,local,shifted);
                if(alt.empty() || alt==ev.rim || !simple_loops(m,alt).first)continue;
                bool alt_forbidden=false;for(auto e:alt)alt_forbidden=alt_forbidden || (e<hc.forbidden.size() && hc.forbidden[e] && !always_cut(t,e));
                if(alt_forbidden)continue;
                const double c=mean_cost(alt);
                if(c<best*(1-1e-9)) {best=c;ev.biased_rim=alt;}
            }
        }
        const double thinness=ev.nominated?std::clamp(1.0-ev.thin_ratio/o.thin_ratio,0.0,1.0):0.0;
        const double sides=std::min(1.0,std::min(ev.front_fraction,ev.back_fraction)/0.3);
        ev.confidence=std::sqrt(thinness)*ev.correspondence*sides*(simple && loops>=1?1.0:0.0);
        if(!ev.nominated)ev.reasons.push_back("not_thin");
        if(std::min(ev.front_fraction,ev.back_fraction)<0.2)ev.reasons.push_back("no_opposite_sides");
        if(ev.correspondence<o.min_correspondence)ev.reasons.push_back("low_correspondence");
        if(!simple || loops==0)ev.reasons.push_back("incoherent_rim");
        if(forbidden)ev.reasons.push_back("forbidden_rim");
        if(ev.confidence<o.min_confidence)ev.reasons.push_back("low_confidence");
        ev.status=ev.reasons.empty()?"candidate":"rejected";
        out.push_back(std::move(ev));
    }
    return out;
}
ProposalGenerator shell_generator() {
    return {name,version,[](const RunState& state,ProposalContext& ctx) {
        std::vector<Proposal> out;
        const auto& o=ctx.options.enhanced.shells;
        if(!o.enabled || !state.surface)return out;
        const auto& mesh=ctx.mesh;const auto& m=*mesh;const auto& t=m.tables();const auto& g=m.geometry();const auto nf=g.face_area.size();
        const auto& hc=ctx.constraints;
        auto hard=[&](std::size_t e) {return (e<hc.locked.size() && hc.locked[e]) || (e<hc.required.size() && hc.required[e]);};
        std::shared_ptr<const TriangleBvh> bvh;
        try {bvh=context_bvh(ctx);}
        catch(const ContractError& error) {
            if(error.code!=ErrorCode::memory_limit)throw;
            ctx.diagnostics[name]={{"version",version},{"skipped","memory_limit"},{"reason",error.what()},{"proposals",0}};
            return out;
        }
        const auto routing=context_routing_cost(ctx);
        // GB m7: the pass's control inside the evidence (per component and smoothing pass).
        DetectorPass work;work.control.cancel=ctx.cancel;work.control.deadline=ctx.deadline;
        if(work.stop()) {
            ctx.diagnostics[name]={{"version",version},{"stopped_by_deadline",true},{"components",Json::array()},{"proposals",0}};
            return out;
        }
        const auto evidence=shell_evidence(m,*bvh,hc,routing,ctx.options,&context_scale(ctx,state),&work);
        Json records=Json::array();std::uint64_t id=0;std::set<EdgeMask> emitted;
        std::map<std::string,std::uint64_t> statuses;bool stopped=work.stopped;
        // Every face of each component (degenerate ones too), built once when a candidate needs it.
        std::vector<std::vector<MeshIndex>> component_faces;
        // Reports stay bounded: the first 32 components in full, every status counted.
        auto keep=[&](Json record) {if(records.size()<32)records.push_back(std::move(record));};
        for(const auto& ev:evidence) {
            Json record={{"evidence",ev.to_json()},{"proposals",Json::array()}};
            ++statuses[ev.status];
            if(ev.status!="candidate") {keep(std::move(record));continue;}
            if(work.stop()) {record["skipped"]="deadline";keep(std::move(record));stopped=true;break;}
            // The component as one patch: its non-hard cuts are dropped (and free to reuse for
            // any companion opening); the rim is added; every side is reopened if needed.
            if(component_faces.empty()) {
                component_faces.resize(g.n_components);
                for(std::size_t f=0;f<nf;++f)component_faces[g.face_component[f]].push_back(static_cast<MeshIndex>(f));
            }
            const auto& faces=component_faces[ev.component];
            std::vector<std::uint8_t> scope(nf,0);for(auto f:faces)scope[f]=1;
            auto patch=state.cuts;auto cost=routing;std::vector<MeshIndex> dropped;
            std::vector<MeshIndex> edges;for(auto f:faces)for(auto c:m.face_corners(f))edges.push_back(t.corner_edge[c]);
            std::sort(edges.begin(),edges.end());edges.erase(std::unique(edges.begin(),edges.end()),edges.end());
            for(auto e:edges) {
                const auto pair=t.edge_faces[e];
                if(!state.cuts[e] || always_cut(t,e) || hard(e) || pair.size()!=2 || !scope[pair[0]] || !scope[pair[1]])continue;
                patch[e]=0;cost[e]=0.0;dropped.push_back(e);
            }
            std::vector<std::pair<std::string,std::vector<MeshIndex>>> rims{{"geometric_rim",ev.rim}};
            if(!ev.biased_rim.empty())rims.push_back({"visibility_biased_rim",ev.biased_rim});
            // At most max_candidates alternatives, candidate zero included.
            const auto limit=static_cast<std::size_t>(std::max<std::int64_t>(0,o.max_candidates-1));
            for(std::size_t r=0;r<rims.size();++r) {
                Json p={{"variant",rims[r].first},{"rim_edges",rims[r].second.size()}};
                if(r>=limit) {p["status"]="beyond_candidate_budget";record["proposals"].push_back(p);continue;}
                auto cuts=patch;for(auto e:rims[r].second)cuts[e]=1;
                const auto opening=open_non_disk_charts(mesh,cuts,cost,hc,ctx.exempt,scope,work.control);
                if(opening.stopped) {p["status"]="stopped";record["proposals"].push_back(p);stopped=true;break;}
                p["openings"]=opening.added;
                if(!opening.ok) {p["status"]="not_openable";p["failures"]=opening.failures;record["proposals"].push_back(p);continue;}
                const auto violations=constraint_violations(m,hc,cuts);
                if(!violations.empty()) {p["status"]="constraint";p["violations"]=violations;record["proposals"].push_back(p);continue;}
                if(cuts==state.cuts || !emitted.insert(cuts).second) {p["status"]="duplicate";record["proposals"].push_back(p);continue;}
                ++id;
                Json provenance={{"component",ev.component},{"variant",rims[r].first},{"rim",rims[r].second},{"openings",opening.added},
                                 {"dropped_interior_cuts",dropped},{"thin_ratio",ev.thin_ratio},{"correspondence",ev.correspondence},
                                 {"confidence",ev.confidence},{"speckles_repaired",ev.speckles_repaired}};
                out.push_back(cut_proposal(id,state,cuts,name,version,std::move(provenance)));
                p["status"]="proposed";p["id"]=id;record["proposals"].push_back(p);
            }
            keep(std::move(record));
            if(stopped)break;
        }
        ctx.diagnostics[name]={{"version",version},
                               {"options",{{"thin_ratio",o.thin_ratio},{"max_candidates",o.max_candidates},{"min_correspondence",o.min_correspondence},
                                           {"min_confidence",o.min_confidence},{"speckle_fraction",o.speckle_fraction},{"smoothing_passes",o.smoothing_passes}}},
                               {"component_status",statuses},{"stopped_by_deadline",stopped},{"components",records},{"proposals",out.size()}};
        return out;
    }};
}
}
