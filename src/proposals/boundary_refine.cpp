// SPDX-License-Identifier: GPL-3.0-or-later
// GB (visual review V1): the baseline's seam-chain refinement applied to the new boundaries of
// compound proposals. Region growth (T34) and the level-set or label borders of rings, cavities
// and shells (T35-T37) are face-region boundaries: on a triangle mesh they follow the
// triangulation staircase, and unlike the baseline's own chains (pipeline refine step:
// refine_chain / refine_cycle) they were never refined, so the transaction judged zig-zags the
// objective's weak turn term could not tell from a straight seam. Here every new or changed chain
// of the proposal's cut set goes through that same refinement, inside the proposal's faces and
// under its hard constraints, before the proposal is emitted: the T28 transaction then evaluates
// the refined cut like any other proposal. Used by policy enhanced_v1 only.
#include "cleanseams/proposals.hpp"
#include "cleanseams/routing.hpp"
#include "cleanseams/physical_scale.hpp"
#include "cleanseams/fields.hpp"
#include "cleanseams/quality.hpp"
#include "cleanseams/cache.hpp"
#include "cleanseams/errors.hpp"
#include <algorithm>
#include <map>
#include <set>
namespace cleanseams {
namespace {
using Json=nlohmann::json;
// The transaction's exact pin rule: no sector of the surface holds two different pins.
bool pins_compatible(const CutSurface& cs,const PinMap& pins) {
    std::map<MeshIndex,Vec2> first;
    for(const auto& [c,p]:pins) {
        const auto s=cs.tables().sector_of_corner[c];const auto [at,inserted]=first.emplace(s,p);
        if(!inserted && at->second!=p)return false;
    }
    return true;
}
bool flag(const EdgeMask& mask,std::size_t e) {return e<mask.size() && mask[e];}
}
nlohmann::json BoundaryRefinement::to_json() const {
    return {{"enabled",enabled},{"chains",chains},{"refined",refined},{"deferred",deferred},{"pin_refusals",pin_refusals},{"passes",passes},
            {"stopped",stopped},{"reverted",reverted},{"corridor_rings",corridor_rings},{"corridor_max_rings",corridor_max_rings},
            {"refined_chains",{{"edges",{edges_before,edges_after}},{"length",{length_before,length_after}},{"turns",{turns_before,turns_after}}}}};
}
BoundaryRefinement refine_compound_boundaries(ProposalContext& ctx,const RunState& parent,EdgeMask& cuts,std::span<const std::uint8_t> scope) {
    BoundaryRefinement out;
    const auto& o=ctx.options;
    out.enabled=o.enhanced.boundary_refinement.enabled && o.refine_paths; // the baseline refines only with refine_paths
    if(!out.enabled || !ctx.mesh)return out;
    const auto& mesh=ctx.mesh;const auto& m=*mesh;const auto& t=m.tables();const auto& g=m.geometry();
    const auto ne=t.edge_vertices.size(),nf=g.face_area.size();
    if(cuts.size()!=ne || parent.cuts.size()!=ne)
        throw ContractError(ErrorCode::invalid_input,"enhanced.boundary_refinement","cut mask does not match the mesh edges");
    if(!scope.empty() && scope.size()!=nf)throw ContractError(ErrorCode::invalid_input,"enhanced.boundary_refinement","scope needs one flag per face");
    auto in_scope=[&](std::size_t e) {
        if(scope.empty())return true;
        for(auto f:t.edge_faces[e])if(!scope[f])return false;
        return true;
    };
    const auto& hc=ctx.constraints;
    // Internal seams (the baseline's `seams`), the parent's, the edges that never move and the edges
    // no corridor may use. Never moved (the baseline's locks: they end or lock chains): locked and
    // required cuts, and every seam the proposal keeps from its parent, so only the proposal's own
    // new boundaries are refined (the parent's seams were refined by the baseline, and a kept seam
    // on a sharp edge never slides off it). No corridor uses those, forbidden edges or any edge with
    // a face outside the proposal's faces.
    EdgeMask seams(ne,0),before(ne,0),locked(ne,0),protect(ne,0);
    for(std::size_t e=0;e<ne;++e) {
        seams[e]=cuts[e] && !t.edge_boundary[e];before[e]=parent.cuts[e] && !t.edge_boundary[e];
        locked[e]=flag(hc.locked,e) || flag(hc.required,e) || (seams[e] && before[e]);
        protect[e]=locked[e] || flag(hc.forbidden,e) || !in_scope(e);
    }
    const std::vector<TransferIndex> anchors(hc.anchors.begin(),hc.anchors.end());
    // The baseline's routing options (pipeline refine step): a triangle-dominant mesh drops the
    // flow term and doubles the turn weight; the native_baseline fixes FIX-01..03 are on.
    auto paths=o.paths;
    std::size_t quads=0;for(std::size_t f=0;f<nf;++f)quads+=m.input().face_offsets[f+1]-m.input().face_offsets[f]==4;
    if(nf>0 && static_cast<double>(quads)/static_cast<double>(nf)<.6) {paths.w_flow=0.;paths.w_turn=2*paths.w_turn;}
    // T39: the corridor radii are ring operators on the proposal faces' frozen reference chart
    // (unchanged in hops mode and with the automatic references, where unit = local length).
    std::vector<MeshIndex> faces;for(std::size_t f=0;f<nf;++f)if(scope.empty() || scope[f])faces.push_back(static_cast<MeshIndex>(f));
    if(faces.empty())return out;
    const auto& scale=context_scale(ctx,parent);
    paths.corridor_rings=std::max<std::int64_t>(1,scale.rings(std::max<std::int64_t>(1,paths.corridor_rings),faces).steps);
    paths.corridor_max_rings=std::max<std::int64_t>(1,scale.rings(std::max<std::int64_t>(1,paths.corridor_max_rings),faces).steps);
    out.corridor_rings=paths.corridor_rings;out.corridor_max_rings=paths.corridor_max_rings;
    // The baseline's seam cost field (seam_cost_field): the raw all-views exposure with its valley
    // and hard-edge (feature) discounts and the explicit viewer, clipped to [0, 1]. The raw field
    // is the pass's measured one (the same directions as the baseline's), never modified.
    std::vector<double> cost;
    if(o.layout.use_exposure && ctx.raw_exposure.size()==ne) {
        ExposureField field;field.edge=ctx.raw_exposure;
        SeamCostOptions options;options.feature_angle=o.layout.prefer_features?60.0:o.layout.feature_angle;options.hide_dir=o.layout.hide_dir;
        cost=seam_cost(m,field,options).cost;
    }
    auto changed=[&](const Chain& chain) {
        bool added=false;
        for(auto e:chain.edges) {if(!in_scope(e))return false;added=added || !before[e];}
        return added;
    };
    auto cs=cut_surface(mesh,seams);
    const bool pins_ok=pins_compatible(*cs,hc.pins);
    RoutingWorkspace workspace;FlankingCache previous;
    const auto sweeps=std::max<std::int64_t>(1,std::min<std::int64_t>(3,o.refine_passes));
    for(std::int64_t sweep=0;sweep<sweeps;++sweep) {
        throw_if_cancelled(ctx.cancel);
        if(ctx.deadline.expired()) {out.stopped=true;break;}
        ++out.passes;
        const SeamGraph graph(mesh,seams,locked,anchors,cs);
        std::vector<const Chain*> chains;for(const auto& chain:graph.chains())chains.push_back(&chain);
        std::stable_sort(chains.begin(),chains.end(),[&](const Chain* a,const Chain* b) {
            return std::pair(-a->edge_length(m),a->id)<std::pair(-b->edge_length(m),b->id);
        });
        std::set<MeshIndex> dirty;std::uint64_t optimized=0;
        for(const auto* chain:chains) {
            throw_if_cancelled(ctx.cancel);
            if(ctx.deadline.expired()) {out.stopped=true;break;}
            if(!changed(*chain))continue;
            if(sweep==0)++out.chains;
            if(std::any_of(chain->verts.begin(),chain->verts.end(),[&](auto v) {return dirty.contains(v);})) {++out.deferred;continue;}
            if(chain->edges.size()<3)continue; // the baseline keeps them (too_short)
            RefineRequest request;request.cut=cs;request.seams=seams;request.options=paths;request.protected_edges=std::span<const std::uint8_t>(protect);
            request.exposure=cost;request.evaluate=nullptr;request.previous=&previous;request.workspace=&workspace;
            request.fixes.weight_variants=true;request.fixes.corridor_bounds=true;request.fixes.full_schedule=true; // FIX-01..03
            auto outcome=chain->cycle && o.refine_cycles?refine_cycle(request,*chain):refine_chain(request,*chain);
            if(outcome.result.status!="optimized" || !outcome.result.cut)continue;
            // Never glue differently pinned sectors (the transaction would refuse the whole proposal).
            if(pins_ok && !pins_compatible(*outcome.result.cut,hc.pins)) {++out.pin_refusals;continue;}
            for(std::size_t e=0;e<ne;++e) {
                if(seams[e] && !outcome.seams[e]) {++out.edges_before;out.length_before+=g.edge_length[e];}
                if(!seams[e] && outcome.seams[e]) {++out.edges_after;out.length_after+=g.edge_length[e];}
            }
            out.turns_before+=outcome.result.turns_old.turns;
            if(outcome.result.turns_new)out.turns_after+=outcome.result.turns_new->turns;
            seams=std::move(outcome.seams);cs=outcome.result.cut;dirty.insert(outcome.result.dirty.begin(),outcome.result.dirty.end());
            ++optimized;++out.refined;
        }
        if(out.stopped || !optimized)break;
    }
    if(!out.refined)return out;
    EdgeMask refined=cuts;
    for(std::size_t e=0;e<ne;++e)if(!t.edge_boundary[e])refined[e]=seams[e];
    // Checked, never assumed: the refined cut set keeps every hard constraint.
    if(!constraint_violations(m,hc,refined).empty()) {out.reverted=true;return out;}
    cuts=std::move(refined);
    return out;
}
}
