// SPDX-License-Identifier: GPL-3.0-or-later
// Phase B geodesic corridors (Section 20, T39): surface distances along a patch's own edges in
// object-space units (or edge counts for the documented hop mode), the band of a distance field
// and the deterministic widening schedule of a search. A distance never crosses empty space:
// another component, a nearby sheet or a folded layer is reached along the surface or not at
// all. The legacy hop corridor of native_baseline stays build_corridor (routing/search.cpp).
#include "cleanseams/physical_scale.hpp"
#include "cleanseams/errors.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <queue>
namespace cleanseams {
namespace {
void require(bool yes,const char* field,const char* message) {if(!yes)throw ContractError(ErrorCode::invalid_input,field,message);}
constexpr double infinity=std::numeric_limits<double>::infinity();
// Patch edges: every face of the edge is in the patch.
std::vector<std::uint8_t> patch_edges(const MeshSnapshot& m,std::span<const std::uint8_t> patch_faces) {
    const auto& t=m.tables();const auto ne=t.edge_vertices.size();
    std::vector<std::uint8_t> out(ne,0);
    for(std::size_t e=0;e<ne;++e) {
        bool inside=!t.edge_faces[e].empty();for(auto f:t.edge_faces[e])inside=inside && patch_faces[f];
        out[e]=inside;
    }
    return out;
}
std::vector<std::uint8_t> patch_vertex_mask(const MeshSnapshot& m,std::span<const std::uint8_t> patch_faces) {
    const auto nf=m.input().face_offsets.size()-1;std::vector<std::uint8_t> out(m.input().positions.size(),0);
    for(std::size_t f=0;f<nf;++f)if(patch_faces[f])for(auto c:m.face_corners(static_cast<MeshIndex>(f)))out[static_cast<std::size_t>(m.input().corner_vertices[c])]=1;
    return out;
}
}
std::string_view distance_metric_name(DistanceMetric metric) {return metric==DistanceMetric::hops?"hops":"physical";}
std::vector<double> patch_distances(const MeshSnapshot& m,std::span<const std::uint8_t> patch_faces,std::span<const MeshIndex> seeds,DistanceMetric metric) {
    const auto& t=m.tables();const auto& g=m.geometry();
    const auto nv=m.input().positions.size(),nf=m.input().face_offsets.size()-1;
    require(patch_faces.size()==nf,"routing.corridor.faces","one flag per face");
    const auto edges=patch_edges(m,patch_faces);const auto vertices=patch_vertex_mask(m,patch_faces);
    std::vector<double> distance(nv,infinity);
    // Dijkstra in (distance, vertex) order: deterministic, and exact for nonnegative lengths.
    using Item=std::pair<double,MeshIndex>;std::priority_queue<Item,std::vector<Item>,std::greater<>> heap;
    for(auto s:seeds) {
        require(s<nv && vertices[s],"routing.corridor.seeds","seeds must be patch vertices");
        if(distance[s]>0) {distance[s]=0;heap.push({0.0,s});}
    }
    while(!heap.empty()) {
        const auto [d,v]=heap.top();heap.pop();
        if(d>distance[v])continue;
        for(auto e:t.vert_edges[v]) {
            if(!edges[e])continue;
            const auto w=m.edge_other(e,v);
            const double step=metric==DistanceMetric::hops?1.0:g.edge_length[e];
            require(std::isfinite(step) && step>=0,"routing.corridor.length","edge lengths must be finite and nonnegative");
            const double nd=d+step;
            if(nd<distance[w]) {distance[w]=nd;heap.push({nd,w});}
        }
    }
    return distance;
}
PatchCorridor corridor_band(const MeshSnapshot& m,std::span<const std::uint8_t> patch_faces,std::span<const std::uint8_t> edge_usable,
                            std::span<const std::uint8_t> vertex_usable,std::span<const double> distance,std::optional<double> radius) {
    const auto& t=m.tables();
    const auto ne=t.edge_vertices.size(),nv=m.input().positions.size(),nf=m.input().face_offsets.size()-1;
    require(patch_faces.size()==nf,"routing.corridor.faces","one flag per face");
    require(edge_usable.size()==ne,"routing.corridor.edges","one flag per mesh edge");
    require(vertex_usable.size()==nv,"routing.corridor.vertices","one flag per mesh vertex");
    require(distance.size()==nv,"routing.corridor.distance","one distance per mesh vertex");
    if(radius)require(std::isfinite(*radius) && *radius>=0,"routing.corridor.radius","radius must be finite and nonnegative");
    const auto edges=patch_edges(m,patch_faces);const auto vertices=patch_vertex_mask(m,patch_faces);
    PatchCorridor out;out.distance.assign(distance.begin(),distance.end());
    if(radius)out.radius=*radius;
    std::vector<std::uint8_t> band(nv,0);
    for(std::size_t v=0;v<nv;++v) {
        out.patch_vertices+=vertices[v];
        band[v]=vertices[v] && std::isfinite(distance[v]) && (!radius || distance[v]<=*radius);
        out.band_vertices+=band[v];
    }
    out.whole_patch=out.band_vertices==out.patch_vertices;
    out.edge_allowed.assign(ne,0);out.vertex_allowed.assign(nv,0);
    for(std::size_t e=0;e<ne;++e) {
        const auto& ends=t.edge_vertices[e];
        out.edge_allowed[e]=edge_usable[e] && edges[e] && band[ends[0]] && band[ends[1]];out.band_edges+=out.edge_allowed[e];
    }
    for(std::size_t v=0;v<nv;++v)out.vertex_allowed[v]=vertex_usable[v] && band[v];
    return out;
}
nlohmann::json CorridorSchedule::to_json() const {
    nlohmann::json r=nlohmann::json::array();for(const auto& x:radii)r.push_back(x?nlohmann::json(*x):nlohmann::json(nullptr));
    return {{"metric",metric},{"radii",r},{"band_vertices",band_vertices},{"patch_vertices",patch_vertices}};
}
CorridorSchedule corridor_schedule(std::span<const double> distance,std::uint64_t patch_vertices,double first_radius,std::uint64_t vertex_budget,
                                   DistanceMetric metric) {
    require(vertex_budget>=1,"routing.corridor.vertex_budget","must be at least 1");
    CorridorSchedule out;out.metric=std::string(distance_metric_name(metric));out.patch_vertices=patch_vertices;
    if(patch_vertices<=vertex_budget) {out.radii={std::nullopt};out.band_vertices={patch_vertices};return out;}
    std::vector<double> reached;for(auto d:distance)if(std::isfinite(d))reached.push_back(d);
    std::sort(reached.begin(),reached.end());
    auto within=[&](double r) {return static_cast<std::uint64_t>(std::upper_bound(reached.begin(),reached.end(),r)-reached.begin());};
    double radius=first_radius;
    if(!(radius>0) || !std::isfinite(radius))radius=reached.empty()?1.0:std::max(reached.back(),1e-300);
    // The first band is always kept; each doubling is kept while it stays within the budget.
    // The radius only grows: a trial cut never shrinks it (Section 20).
    out.radii.push_back(radius);out.band_vertices.push_back(within(radius));
    while(out.radii.size()<64 && within(2.0*radius)<=vertex_budget) {radius*=2.0;out.radii.push_back(radius);out.band_vertices.push_back(within(radius));}
    return out;
}
}
