// SPDX-License-Identifier: GPL-3.0-or-later
// Phase B preference-aware region proposals (T34, Section 15). Region growth uses seam
// desirability with the correct sign: a preferred seam location is expensive for a growing
// region to cross, so neighbouring regions meet there. The primal routing cost (cheap to
// follow as a seam) and the dual crossing weight (expensive to cross) stay distinct arrays.
// Forbidden seams contract their faces into one dual node before seeding and growth;
// required and locked cuts are barriers. Used by policy enhanced_v1 only.
#include "cleanseams/proposals.hpp"
#include "cleanseams/physical_scale.hpp"
#include "cleanseams/cache.hpp"
#include "cleanseams/errors.hpp"
#include "cleanseams/detail/mesh_math.hpp"
#include <algorithm>
#include <cmath>
#include <map>
#include <numbers>
#include <numeric>
#include <set>
namespace cleanseams {
namespace {
using Json=nlohmann::json;
constexpr const char* name="region_preferences";
constexpr const char* version="1";
// GB (visual review V1): ids of the unrefined fallbacks (above every refined proposal's id).
constexpr std::uint64_t fallback_id_offset=1'000'000;
bool always_cut(const MeshTables& t,std::size_t e) {
    return t.edge_boundary[e] || !t.edge_manifold[e] || !t.edge_consistent[e] || t.edge_corners[e].size()!=2;
}
double median(std::vector<double> v) {
    if(v.empty())return 0.0;
    std::sort(v.begin(),v.end());const auto n=v.size();
    return n%2?v[n/2]:0.5*(v[n/2-1]+v[n/2]);
}
struct Target {
    std::string kind;                // "split" (over-budget chart) or "pair" (adjacent charts)
    std::vector<MeshIndex> charts;   // incumbent charts
    std::vector<MeshIndex> faces;    // ascending
    double priority=0;
};
// The node of `dual` closest to the area-weighted centre and normal of `faces` (lloyd()'s
// seed score: squared distance over the largest, plus a quarter of the squared normal gap).
MeshIndex proxy_node(const RegionDual& dual,const MeshSnapshot& m,std::span<const MeshIndex> faces,double quantum) {
    const auto& g=m.geometry();Vec3 centre{0,0,0},normal{0,0,0};double area=0;
    for(auto f:faces) {const double a=g.face_area[f];area+=a;for(std::size_t d=0;d<3;++d) {centre[d]+=g.face_center[f][d]*a;normal[d]+=g.face_normal[f][d]*a;}}
    for(auto& x:centre)x/=std::max(area,1e-300);
    const double length=std::max(detail::norm(normal),1e-300);for(auto& x:normal)x/=length;
    std::vector<MeshIndex> nodes;for(auto f:faces)nodes.push_back(dual.node(f));
    std::sort(nodes.begin(),nodes.end());nodes.erase(std::unique(nodes.begin(),nodes.end()),nodes.end());
    std::vector<double> d2(nodes.size());double biggest=-INFINITY;
    for(std::size_t i=0;i<nodes.size();++i) {const auto dx=detail::sub(dual.node_center[nodes[i]],centre);d2[i]=(dx[0]*dx[0]+dx[1]*dx[1])+dx[2]*dx[2];biggest=std::max(biggest,d2[i]);}
    std::size_t pick=0;double best=INFINITY;
    for(std::size_t i=0;i<nodes.size();++i) {
        const double angular=1-std::clamp(detail::dot(dual.node_normal[nodes[i]],normal),-1.0,1.0);
        // T39: on the score grid under the physical scale (exact symmetric ties stay ties).
        const double score=snap_to_quantum(d2[i]/std::max(biggest,1e-300)+.25*(angular*angular),quantum>0?score_quantum:0.0);
        if(score<best) {best=score;pick=i;}
    }
    return nodes[pick];
}
MeshIndex farthest(const std::vector<double>& distance) {
    MeshIndex best=0;double d=-1;
    for(std::size_t i=0;i<distance.size();++i)if(std::isfinite(distance[i]) && distance[i]>d) {d=distance[i];best=static_cast<MeshIndex>(i);}
    return best;
}
// Every connected component of the dual keeps a seed (component coverage, as the legacy
// farthest-point seeding guarantees); returns the seeds added.
std::uint64_t cover_components(const RegionDual& dual,std::vector<MeshIndex>& seeds) {
    const auto n=dual.graph.size();std::vector<std::int64_t> component(n,-1);std::int64_t count=0;
    for(std::size_t s=0;s<n;++s) {
        if(component[s]>=0)continue;
        std::vector<std::size_t> stack{s};component[s]=count;
        while(!stack.empty()) {
            const auto v=stack.back();stack.pop_back();
            for(auto j=dual.graph.indptr[v];j<dual.graph.indptr[v+1];++j) {const auto w=dual.graph.indices[j];if(component[w]<0) {component[w]=count;stack.push_back(w);}}
        }
        ++count;
    }
    std::vector<std::uint8_t> seeded(static_cast<std::size_t>(count),0);for(auto s:seeds)seeded[static_cast<std::size_t>(component[s])]=1;
    std::uint64_t added=0;
    for(std::size_t v=0;v<n;++v)if(!seeded[static_cast<std::size_t>(component[v])]) {seeded[static_cast<std::size_t>(component[v])]=1;seeds.push_back(static_cast<MeshIndex>(v));++added;}
    return added;
}
}
std::vector<double> preferred_seam_strength(const MeshSnapshot& mesh,const ResolvedPreferences& prefs,const PreferenceOptions& options,
                                            std::span<const MeshIndex> reference) {
    const auto& g=mesh.geometry();const auto ne=g.edge_length.size();
    const auto cost=enhanced_edge_cost(mesh,prefs,options);
    // The routing multiplier per edge (a degenerate edge carries no preference).
    std::vector<double> m(ne,1.0);
    for(std::size_t e=0;e<ne;++e)if(g.edge_length[e]>0)m[e]=cost[e]/g.edge_length[e];
    std::vector<double> sample;
    if(reference.empty())sample=m;
    else for(auto e:reference) {
        if(e>=ne)throw ContractError(ErrorCode::invalid_input,"regions.reference","reference edge out of range");
        sample.push_back(m[e]);
    }
    const double ref=median(std::move(sample));std::vector<double> strength(ne,0.0);
    if(ref>0 && std::isfinite(ref))for(std::size_t e=0;e<ne;++e)strength[e]=std::clamp(1.0-m[e]/ref,0.0,1.0);
    return strength;
}
ProposalGenerator region_preference_generator() {
    return {name,version,[](const RunState& state,ProposalContext& ctx) {
        std::vector<Proposal> out;
        const auto& o=ctx.options.enhanced.regions;
        if(!o.enabled || !state.surface)return out;
        const auto& mesh=ctx.mesh;const auto& m=*mesh;const auto& t=m.tables();const auto& g=m.geometry();
        const auto ne=t.edge_vertices.size(),nf=g.face_area.size();const auto& ct=state.surface->tables();
        const auto& hc=ctx.constraints;
        // GB (visual review V3): an incumbent seam on a hard feature edge (the baseline seam cost's
        // hard edges: |dihedral| >= its feature angle) bounds the growth like a required or locked
        // cut, so the regrown borders never move a seam off a sharp edge onto a smooth face.
        const double feature_angle=(ctx.options.layout.prefer_features?60.0:ctx.options.layout.feature_angle)*std::numbers::pi/180.0;
        // A painted forbidden edge is never kept: the hard constraint overrides the preference.
        auto feature_seam=[&](std::size_t e) {
            return o.keep_feature_seams && feature_angle>0 && state.cuts[e] && std::abs(g.edge_dihedral[e])>=feature_angle &&
                   !(e<hc.forbidden.size() && hc.forbidden[e]);
        };
        std::uint64_t feature_seams=0;
        auto constrained=[&](std::size_t e) {return (e<hc.locked.size() && hc.locked[e]) || (e<hc.required.size() && hc.required[e]);};
        // Barriers of the growth and of the pair targets: the hard cuts and the feature seams. A feature
        // seam stays cut where the regrown regions meet across it; one that ends inside a region
        // (a slit of it) is not forced to stay, so no dead-end cut is left in a regrown chart.
        auto hard=[&](std::size_t e) {return constrained(e) || feature_seam(e);};
        const auto prefs=context_preferences(ctx);const auto routing=context_routing_cost(ctx);
        // Strength relative to the median routing multiplier of the mesh's interior edges.
        std::vector<MeshIndex> interior;for(std::size_t e=0;e<ne;++e)if(!always_cut(t,e))interior.push_back(static_cast<MeshIndex>(e));
        const auto strength=preferred_seam_strength(m,*prefs,ctx.options.enhanced.preferences,interior);
        std::uint64_t preferred=0;for(auto s:strength)preferred+=s>0;
        // Targets: over-budget charts first (a split alternative), then adjacent chart pairs
        // whose shared seam is not entirely hard, the longest unwelcome shared seam first.
        std::vector<Target> splits,pairs;
        const auto& scale=context_scale(ctx,state);
        for(MeshIndex ch=0;ch<ct.n_charts;++ch) {
            const auto& c=state.charts[ch];
            if(c.exempt_only || !c.measured || c.within_budget)continue;
            const auto faces=ct.chart_faces[ch];
            // T39: the split floor is an area under the physical scale (the face count in hops mode).
            double area=0;for(auto f:faces)area+=g.face_area[f];
            if(scale.floor(std::max<std::int64_t>(2,ctx.options.layout.split_min_faces),faces).below(area,faces.size()))continue;
            splits.push_back({"split",{ch},std::vector<MeshIndex>(faces.begin(),faces.end()),c.budget_ratio});
        }
        for(const auto& adj:ct.adjacency) {
            const auto a=adj.charts[0],b=adj.charts[1];
            if(a==b || state.charts[a].exempt_only || state.charts[b].exempt_only)continue;
            double unwelcome=0;bool soft=false;
            for(auto e:adj.edges)if(!hard(e) && !always_cut(t,e)) {soft=true;unwelcome+=g.edge_length[e]*(1.0-strength[e]);}
            if(!soft)continue;
            std::vector<MeshIndex> faces;
            for(auto ch:{a,b})for(auto f:ct.chart_faces[ch])faces.push_back(f);
            std::sort(faces.begin(),faces.end());
            pairs.push_back({"pair",{std::min(a,b),std::max(a,b)},std::move(faces),unwelcome});
        }
        auto by_priority=[](const Target& x,const Target& y) {return x.priority!=y.priority?x.priority>y.priority:x.charts<y.charts;};
        std::stable_sort(splits.begin(),splits.end(),by_priority);std::stable_sort(pairs.begin(),pairs.end(),by_priority);
        std::vector<Target> targets=splits;for(auto& p:pairs)targets.push_back(std::move(p));
        const auto available=targets.size();
        if(static_cast<std::int64_t>(targets.size())>o.max_targets)targets.resize(static_cast<std::size_t>(o.max_targets));
        EdgeMask barriers(ne,0),weld(ne,0);
        for(std::size_t e=0;e<ne;++e) {
            barriers[e]=hard(e) || always_cut(t,e);weld[e]=e<hc.forbidden.size() && hc.forbidden[e];
            feature_seams+=!always_cut(t,e) && feature_seam(e);
        }
        Json records=Json::array();std::map<EdgeMask,std::size_t> emitted; // cut set -> index in `out`
        for(std::size_t ti=0;ti<targets.size();++ti) {
            if(ctx.deadline.expired()) {records.push_back({{"target",ti},{"skipped","deadline"}});break;}
            throw_if_cancelled(ctx.cancel);
            const auto& target=targets[ti];
            Json record={{"kind",target.kind},{"charts",target.charts},{"faces",target.faces.size()},{"priority",target.priority},{"variants",Json::array()}};
            std::vector<std::uint8_t> scope(nf,0);for(auto f:target.faces)scope[f]=1;
            for(int variant=0;variant<2;++variant) {
                const std::string variant_name=variant==0?"preference":"baseline";
                const double w_pref=variant==0?o.w_preference:0.0;
                Json v={{"variant",variant_name},{"w_preference",w_pref}};
                auto dual=region_dual(m,target.faces,barriers,weld,strength,o.w_normal,w_pref);
                // T39: crossing weights on the target's length grid (physical scale; exact in
                // hops mode), so seeds and growth never follow coordinate rounding noise.
                const double quantum=scale.length_quantum(target.faces);
                for(auto& w:dual.graph.weights)w=snap_to_quantum(w,quantum);
                std::vector<MeshIndex> seeds;
                if(target.kind=="pair") {
                    for(auto ch:target.charts) {const auto faces=ct.chart_faces[ch];seeds.push_back(proxy_node(dual,m,faces,quantum));}
                } else {
                    // Legacy split_regions seeding on this variant's weights: the farthest node
                    // from the chart's first face, then the farthest from that one.
                    std::vector<double> distance(dual.graph.size(),INFINITY);
                    update_distance(dual.graph,dual.node(target.faces.front()),distance);
                    const auto first=farthest(distance);
                    std::fill(distance.begin(),distance.end(),INFINITY);update_distance(dual.graph,first,distance);
                    seeds={first,farthest(distance)};
                }
                if(seeds.size()<2 || seeds[0]==seeds[1]) {v["status"]="refused";v["reason"]="degenerate_seeds";record["variants"].push_back(v);continue;}
                const auto covered=cover_components(dual,seeds);
                const auto grown=region_lloyd(dual,seeds,o.lloyd_steps,ctx.deadline,quantum);
                // Labels back to the original faces; the target's non-hard cuts are replaced by
                // the region borders (required/locked cuts stay: they were barriers).
                std::vector<std::int64_t> label(nf,-1);
                for(std::size_t i=0;i<dual.faces.size();++i)label[dual.faces[i]]=grown.labels[dual.node_of_face[i]];
                auto cuts=state.cuts;
                for(std::size_t e=0;e<ne;++e) {
                    const auto pair=t.edge_faces[e];
                    if(always_cut(t,e) || pair.size()!=2 || label[pair[0]]<0 || label[pair[1]]<0)continue;
                    cuts[e]=constrained(e) || label[pair[0]]!=label[pair[1]];
                }
                std::vector<MeshIndex> seed_faces;for(auto s:grown.seeds)seed_faces.push_back(dual.node_faces[s][0]);
                v["seed_faces"]=seed_faces;v["regions"]=grown.seeds.size();v["coverage_seeds"]=covered;v["welded_edges"]=dual.welded_edges;
                v["nodes"]=dual.graph.size();
                // Exact disk opening of every chart the new partition touches (companion cuts).
                const auto opening=open_non_disk_charts(mesh,cuts,routing,hc,ctx.exempt,scope,Control{ctx.cancel,ctx.deadline}); // GB m7
                v["openings"]=opening.added;
                if(!opening.ok) {v["status"]="refused";v["reason"]="not_openable";v["failures"]=opening.failures;record["variants"].push_back(v);continue;}
                const auto violations=constraint_violations(m,hc,cuts);
                if(!violations.empty()) {v["status"]="refused";v["reason"]="constraint";v["violations"]=violations;record["variants"].push_back(v);continue;}
                // GB (visual review V1): the region borders and openings follow the triangulation's
                // staircase; they go through the baseline's chain refinement inside the target (its
                // hard constraints and topology kept) before the transaction evaluates them.
                const auto raw=cuts;
                const auto refinement=refine_compound_boundaries(ctx,state,cuts,scope);
                v["boundary_refinement"]=refinement.to_json();
                if(cuts==state.cuts) {v["status"]="unchanged";record["variants"].push_back(v);continue;}
                if(const auto found=emitted.find(cuts);found!=emitted.end()) {
                    auto& same=out[found->second];same.provenance["variants"].push_back(variant_name);
                    v["status"]="duplicate";v["same_as"]=same.identity.id;record["variants"].push_back(v);continue;
                }
                const auto id=static_cast<std::uint64_t>(1+2*ti+static_cast<std::size_t>(variant));
                Json provenance={{"target",target.kind},{"charts",target.charts},{"variants",Json::array({variant_name})},
                                 {"coefficients",{{"w_normal",o.w_normal},{"w_preference",w_pref},{"lloyd_steps",o.lloyd_steps}}},
                                 {"strength_version",region_strength_version},{"seed_faces",seed_faces},{"openings",opening.added},
                                 {"boundary_refinement",refinement.to_json()}};
                emitted.emplace(cuts,out.size());
                out.push_back(cut_proposal(id,state,cuts,name,version,provenance));
                v["status"]="proposed";v["id"]=id;
                // The refinement has no distortion guard of its own (the baseline's per-route guard
                // solves every flanking chart of every tried route): the unrefined partition is
                // proposed too, after it (a higher id, so ties keep the refined one). The common
                // ordering then keeps the refined cut unless it fails a check or a budget the raw
                // one meets; a raw fallback with a higher J than a compliant refined leader is
                // refused as dominated before its solve.
                if(cuts!=raw && raw!=state.cuts && !emitted.contains(raw)) {
                    const auto fallback=id+fallback_id_offset;
                    provenance["boundary_refinement"]={{"unrefined_fallback_of",id}};
                    emitted.emplace(raw,out.size());
                    out.push_back(cut_proposal(fallback,state,raw,name,version,std::move(provenance)));
                    v["unrefined_fallback_id"]=fallback;
                }
                record["variants"].push_back(v);
            }
            records.push_back(std::move(record));
        }
        ctx.diagnostics[name]={{"version",version},{"strength",{{"version",region_strength_version},{"preferred_edges",preferred}}},
                               {"coefficients",{{"w_normal",o.w_normal},{"w_preference",o.w_preference},{"lloyd_steps",o.lloyd_steps}}},
                               {"targets_available",available},{"max_targets",o.max_targets},{"targets",records},{"proposals",out.size()},
                               {"feature_seams",{{"kept",o.keep_feature_seams},{"feature_angle_deg",feature_angle*180.0/std::numbers::pi},{"edges",feature_seams}}}};
        return out;
    }};
}
}
