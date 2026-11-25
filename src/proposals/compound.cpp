// SPDX-License-Identifier: GPL-3.0-or-later
// Phase B compound cuts shared by the region (T34), ring (T35), cavity (T36) and shell (T37)
// generators: admissible cheapest paths, exact disk opening of every chart a compound
// proposal touches (its companion cuts), and the single delta that carries all of them to
// the common transaction. Used by policy enhanced_v1 only.
#include "cleanseams/proposals.hpp"
#include "cleanseams/geometry_bvh.hpp"
#include "cleanseams/cache.hpp"
#include "cleanseams/errors.hpp"
#include <algorithm>
#include <map>
#include <numeric>
#include <queue>
#include <set>
namespace cleanseams {
namespace {
using Json=nlohmann::json;
bool always_cut(const MeshTables& t,std::size_t e) {
    return t.edge_boundary[e] || !t.edge_manifold[e] || !t.edge_consistent[e] || t.edge_corners[e].size()!=2;
}
bool measured_chart(const CutTables& ct,MeshIndex ch,const MeshGeometry& g,std::span<const std::uint8_t> exempt) {
    for(auto f:ct.chart_faces[ch])if(!(f<exempt.size() && exempt[f]) && g.face_degenerate_kind[f]==0)return true;
    return false;
}
struct UnionFind {
    std::vector<std::size_t> parent;
    explicit UnionFind(std::size_t n):parent(n) {std::iota(parent.begin(),parent.end(),std::size_t{0});}
    std::size_t find(std::size_t x) {while(parent[x]!=x) {parent[x]=parent[parent[x]];x=parent[x];}return x;}
    bool join(std::size_t a,std::size_t b) {a=find(a);b=find(b);if(a==b)return false;parent[std::max(a,b)]=std::min(a,b);return true;}
};
}
nlohmann::json CompanionOpening::to_json() const {return {{"ok",ok},{"stopped",stopped},{"added",added},{"failures",failures}};}
double DetectorPass::area(const MeshSnapshot& m) {
    if(!mesh_area)mesh_area=std::accumulate(m.geometry().face_area.begin(),m.geometry().face_area.end(),0.0);
    return *mesh_area;
}
std::span<const double> DetectorPass::diagonals(const MeshSnapshot& m) {
    if(component_diagonals.size()!=m.geometry().n_components)component_diagonals=cleanseams::component_diagonals(m);
    return component_diagonals;
}
std::vector<MeshIndex> cheapest_edge_path(const MeshSnapshot& mesh,std::span<const MeshIndex> sources,std::span<const std::uint8_t> targets,
                                          const std::function<bool(MeshIndex)>& usable,std::span<const double> cost) {
    const auto& t=mesh.tables();const auto nv=mesh.input().positions.size();
    if(targets.size()!=nv)throw ContractError(ErrorCode::invalid_input,"compound.path","target mask needs one flag per vertex");
    if(cost.size()!=t.edge_vertices.size())throw ContractError(ErrorCode::invalid_input,"compound.path","cost needs one value per edge");
    std::vector<double> dist(nv,std::numeric_limits<double>::infinity());std::vector<std::int64_t> via(nv,-1);
    // T39: a settled vertex keeps its predecessor. Before, an equal-distance tie across a
    // zero-cost edge (dropped cuts are priced at zero) could re-point a settled vertex at a
    // later one, so two vertices pointed at each other and the path reconstruction never ended
    // (found by the renumbering property test on F25).
    std::vector<std::uint8_t> settled(nv,0);
    using Item=std::pair<double,MeshIndex>;std::priority_queue<Item,std::vector<Item>,std::greater<>> heap;
    for(auto s:sources)if(s<nv && dist[s]>0) {dist[s]=0;heap.push({0.0,s});}
    while(!heap.empty()) {
        const auto [d,v]=heap.top();heap.pop();
        if(d>dist[v] || settled[v])continue;
        settled[v]=1;
        if(targets[v]) {
            std::vector<MeshIndex> edges;auto at=v;
            while(via[at]>=0) {
                const auto e=static_cast<MeshIndex>(via[at]);edges.push_back(e);at=mesh.edge_other(e,at);
                if(edges.size()>nv)throw ContractError(ErrorCode::numerical_failure,"compound.path","predecessor cycle");
            }
            std::reverse(edges.begin(),edges.end());
            return edges;
        }
        for(auto e:t.vert_edges[v]) {
            if(!usable(e))continue;
            const double w=cost[e];
            if(!std::isfinite(w) || w<0)throw ContractError(ErrorCode::invalid_input,"compound.path","path costs must be finite and nonnegative");
            const auto u=mesh.edge_other(e,v);if(settled[u])continue;
            const double nd=d+w;
            if(nd<dist[u] || (nd==dist[u] && via[u]>=0 && e<static_cast<MeshIndex>(via[u]))) {dist[u]=nd;via[u]=e;heap.push({nd,u});}
        }
    }
    return {};
}
CompanionOpening open_non_disk_charts(const std::shared_ptr<const MeshSnapshot>& mesh,EdgeMask& cuts,std::span<const double> cost,
                                      const HardConstraints& hc,std::span<const std::uint8_t> exempt,std::span<const std::uint8_t> scope,
                                      const Control& control) {
    const auto& t=mesh->tables();const auto& g=mesh->geometry();const auto nv=mesh->input().positions.size();
    if(cuts.size()!=t.edge_vertices.size())throw ContractError(ErrorCode::invalid_input,"compound.cuts","cut mask does not match the mesh edges");
    CompanionOpening out;std::set<MeshIndex> added,failed; // failed: charts by smallest face
    auto admissible=[&](MeshIndex e) {return !always_cut(t,e) && !(e<hc.forbidden.size() && hc.forbidden[e]);};
    for(int round=0;round<4096;++round) {
        // GB m7: every round rebuilds the whole surface; the pass's control is checked first.
        if(control.check()) {
            out.ok=false;out.stopped=true;out.failures.push_back({{"reason","deadline"},{"round",round}});
            out.added.assign(added.begin(),added.end());
            return out;
        }
        const auto surface=cut_surface(mesh,cuts);const auto& ct=surface->tables();
        std::optional<MeshIndex> target;
        for(MeshIndex ch=0;ch<ct.n_charts && !target;++ch) {
            if(ct.chart_is_disk[ch] || failed.contains(ct.chart_faces[ch][0]))continue;
            if(!measured_chart(ct,ch,g,exempt))continue;
            if(!scope.empty()) {
                bool inside=false;for(auto f:ct.chart_faces[ch])inside=inside || scope[f];
                if(!inside)continue;
            }
            if(ct.chart_pinched[ch] || ct.chart_unsupported[ch]) {
                failed.insert(ct.chart_faces[ch][0]);
                out.failures.push_back({{"chart_face",ct.chart_faces[ch][0]},{"faces",ct.chart_faces[ch].size()},
                                        {"reason",ct.chart_pinched[ch]?"pinched":"unsupported"}});
                continue;
            }
            target=ch;
        }
        if(!target)break;
        const auto ch=*target;bool changed=false;
        auto interior=[&](MeshIndex e) {
            if(cuts[e] || !admissible(e))return false;
            for(auto f:t.edge_faces[e])if(ct.chart_of_face[f]!=ch)return false;
            return true;
        };
        const auto loops=ct.chart_loops[ch];
        if(loops.size()>=2) {
            // Two boundary loops: the cheapest admissible interior path joins them.
            const auto start=surface->loop_verts(loops[0]);std::vector<std::uint8_t> goal(nv,0);
            for(std::size_t k=1;k<loops.size();++k)for(auto v:surface->loop_verts(loops[k]))goal[v]=1;
            for(auto v:start)goal[v]=0;
            for(auto e:cheapest_edge_path(*mesh,start,goal,interior,cost)) {cuts[e]=1;added.insert(e);changed=true;}
        } else if(loops.empty()) {
            // A closed chart: its cheapest admissible interior edge opens it (a sphere becomes a
            // disk; a handle is then opened by the cut graph in the next round).
            std::optional<MeshIndex> best;
            for(auto f:ct.chart_faces[ch])for(auto c:mesh->face_corners(f)) {
                const auto e=t.corner_edge[c];
                if(interior(e) && (!best || cost[e]<cost[*best] || (cost[e]==cost[*best] && e<*best)))best=e;
            }
            if(best) {cuts[*best]=1;added.insert(*best);changed=true;}
        }
        if(!changed) {
            // Cut graph: a dual spanning tree keeps forbidden edges first and expensive edges
            // whole; every other interior edge is cut, then dangling cuts are glued back.
            const auto faces=ct.chart_faces[ch];std::map<MeshIndex,std::size_t> local;
            for(std::size_t i=0;i<faces.size();++i)local.emplace(faces[i],i);
            std::vector<MeshIndex> dual;
            for(auto f:faces)for(auto c:mesh->face_corners(f)) {const auto e=t.corner_edge[c];if(!cuts[e] && !always_cut(t,e))dual.push_back(e);}
            std::sort(dual.begin(),dual.end());dual.erase(std::unique(dual.begin(),dual.end()),dual.end());
            std::vector<MeshIndex> order=dual;
            auto forbidden=[&](MeshIndex e) {return e<hc.forbidden.size() && hc.forbidden[e];};
            std::stable_sort(order.begin(),order.end(),[&](MeshIndex a,MeshIndex b) {
                if(forbidden(a)!=forbidden(b))return forbidden(a);
                return cost[a]>cost[b];
            });
            UnionFind uf(faces.size());std::vector<MeshIndex> cut_graph;
            for(auto e:order) {
                const auto pair=t.edge_faces[e];
                if(!uf.join(local.at(pair[0]),local.at(pair[1])) && !forbidden(e))cut_graph.push_back(e);
            }
            // Cut degree of every vertex: the chart's existing cut edges plus the cut graph.
            std::map<MeshIndex,std::int64_t> degree;std::set<MeshIndex> chart_cut;
            for(auto f:faces)for(auto c:mesh->face_corners(f)) {const auto e=t.corner_edge[c];if(cuts[e])chart_cut.insert(e);}
            for(auto e:chart_cut)for(auto v:t.edge_vertices[e])++degree[v];
            std::set<MeshIndex> kept(cut_graph.begin(),cut_graph.end());
            for(auto e:kept)for(auto v:t.edge_vertices[e])++degree[v];
            for(bool pruned=true;pruned;) {
                pruned=false;
                for(auto it=kept.begin();it!=kept.end();) {
                    const auto [a,b]=t.edge_vertices[*it];
                    if(degree[a]<=1 || degree[b]<=1) {--degree[a];--degree[b];it=kept.erase(it);pruned=true;}
                    else ++it;
                }
            }
            for(auto e:kept) {if(!cuts[e]) {cuts[e]=1;added.insert(e);changed=true;}}
        }
        if(!changed) {
            failed.insert(ct.chart_faces[ch][0]);
            out.failures.push_back({{"chart_face",ct.chart_faces[ch][0]},{"faces",ct.chart_faces[ch].size()},{"reason","no_admissible_opening"},
                                    {"loops",loops.size()},{"euler",ct.chart_euler[ch]}});
        }
    }
    // Final exact check: every measured chart in scope is a disk.
    const auto surface=cut_surface(mesh,cuts);const auto& ct=surface->tables();
    for(MeshIndex ch=0;ch<ct.n_charts;++ch) {
        if(ct.chart_is_disk[ch] || !measured_chart(ct,ch,g,exempt))continue;
        bool inside=scope.empty();for(auto f:ct.chart_faces[ch])inside=inside || scope[f];
        if(!inside)continue;
        out.ok=false;
        if(!failed.contains(ct.chart_faces[ch][0]))
            out.failures.push_back({{"chart_face",ct.chart_faces[ch][0]},{"faces",ct.chart_faces[ch].size()},{"reason","not_disk"}});
    }
    out.added.assign(added.begin(),added.end());
    return out;
}
Proposal cut_proposal(std::uint64_t id,const RunState& parent,std::span<const std::uint8_t> cuts,std::string generator,std::string version,
                      nlohmann::json provenance) {
    const auto& t=parent.surface->mesh().tables();
    if(cuts.size()!=parent.cuts.size())throw ContractError(ErrorCode::invalid_input,"compound.cuts","cut mask does not match the mesh edges");
    Proposal p;p.identity.id=id;p.identity.parent_revision=parent.revision;p.identity.generator=std::move(generator);
    p.identity.generator_version=std::move(version);
    for(std::size_t e=0;e<cuts.size();++e) {
        if(always_cut(t,e))continue;
        if(cuts[e] && !parent.cuts[e])p.identity.added_edges.push_back(static_cast<std::int64_t>(e));
        if(!cuts[e] && parent.cuts[e])p.identity.removed_edges.push_back(static_cast<std::int64_t>(e));
    }
    p.provenance=std::move(provenance);
    return p;
}
std::shared_ptr<const ResolvedPreferences> context_preferences(const ProposalContext& ctx) {
    if(ctx.preferences)return ctx.preferences;
    auto neutral=resolve_preferences(*ctx.mesh,{},ctx.options.enhanced.preferences,std::nullopt,{},ctx.options.layout.feature_angle,ctx.options.layout.hide_dir);
    neutral.raw_exposure=ctx.raw_exposure;
    return std::make_shared<const ResolvedPreferences>(std::move(neutral));
}
std::vector<double> context_routing_cost(const ProposalContext& ctx) {
    return enhanced_edge_cost(*ctx.mesh,*context_preferences(ctx),ctx.options.enhanced.preferences);
}
std::shared_ptr<const TriangleBvh> context_bvh(ProposalContext& ctx) {
    if(!ctx.mesh)throw ContractError(ErrorCode::invalid_input,"bvh","missing mesh");
    if(!ctx.bvh || !ctx.bvh->serves(*ctx.mesh))ctx.bvh=std::make_shared<const TriangleBvh>(ctx.mesh);
    return ctx.bvh;
}
}
