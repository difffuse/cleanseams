// SPDX-License-Identifier: GPL-3.0-or-later
#include "cleanseams/panels.hpp"
#include "cleanseams/cache.hpp"
#include "cleanseams/fields.hpp"
#include "cleanseams/errors.hpp"
#include "cleanseams/detail/mesh_math.hpp"
#include "cleanseams/detail/numpy_compat.hpp"
#include <algorithm>
#include <cmath>
#include <map>
#include <numbers>
#include <queue>
#include <set>
namespace cleanseams {
namespace {
void require(bool yes,const char* field,const char* message) {if(!yes)throw ContractError(ErrorCode::invalid_input,field,message);}
double loops_total(const CutTables& ct) {double n=0;for(auto k:ct.chart_n_loops)n+=static_cast<double>(k);return n;}
}
std::optional<std::vector<MeshIndex>> join_boundary_loops(const CutSurface& cs,MeshIndex chart,std::optional<std::span<const double>> cost) {
    const auto& m=cs.mesh();const auto& ct=cs.tables();const auto& t=m.tables();const auto& g=m.geometry();
    require(chart<ct.n_charts,"panels.chart","chart index out of range");
    if(cost)require(cost->size()==t.edge_vertices.size(),"panels.cost","cost does not match the mesh edges");
    const auto loops=ct.chart_loops[chart];
    if(loops.size()<2)return std::nullopt;
    MeshIndex first=loops[0];for(auto li:loops)if(ct.loops[li].size()<ct.loops[first].size())first=li;
    std::set<MeshIndex> seeds,targets;
    for(auto s:cs.loop_sectors(first))seeds.insert(s);
    for(auto li:loops)if(li!=first)for(auto s:cs.loop_sectors(li))targets.insert(s);
    struct Arc {MeshIndex other,edge;double weight;};
    std::map<MeshIndex,std::vector<Arc>> graph;
    for(std::size_t edge=0;edge<ct.cut.size();++edge) {
        if(ct.cut[edge])continue;
        const auto corner=t.edge_corners[edge][0];
        if(ct.chart_of_corner[corner]!=chart)continue;
        const auto a=ct.sector_of_corner[corner],b=ct.sector_of_corner[m.corner_next(corner)];
        double weight=g.edge_length[edge];
        if(cost)weight*=1+.5*(*cost)[edge];
        graph[a].push_back({b,static_cast<MeshIndex>(edge),weight});graph[b].push_back({a,static_cast<MeshIndex>(edge),weight});
    }
    std::map<MeshIndex,double> distance;std::map<MeshIndex,std::pair<MeshIndex,MeshIndex>> previous;
    using Item=std::pair<double,MeshIndex>;std::priority_queue<Item,std::vector<Item>,std::greater<Item>> queue;
    for(auto s:seeds) {distance[s]=0.0;queue.push({0.0,s});}
    while(!queue.empty()) {
        const auto [value,sector]=queue.top();queue.pop();
        if(value>distance.at(sector))continue;
        if(targets.contains(sector)) {
            std::vector<MeshIndex> path;auto at=sector;
            while(!seeds.contains(at)) {const auto [from,edge]=previous.at(at);path.push_back(edge);at=from;}
            if(path.empty())return std::nullopt;
            return path;
        }
        const auto found=graph.find(sector);if(found==graph.end())continue;
        for(const auto& arc:found->second) {
            const double candidate=value+arc.weight;const auto known=distance.find(arc.other);
            if(candidate<(known==distance.end()?INFINITY:known->second)) {
                distance[arc.other]=candidate;previous[arc.other]={sector,arc.edge};queue.push({candidate,arc.other});
            }
        }
    }
    return std::nullopt;
}
std::optional<CadPanelCuts> cad_panel_cuts(const std::shared_ptr<const MeshSnapshot>& mesh,std::span<const std::uint8_t> initial_cut,
                                           std::optional<std::span<const double>> cost,const Deadline& deadline,bool allow_triangles) {
    require(bool(mesh),"panels.mesh","missing mesh snapshot");const auto& m=*mesh;const auto& t=m.tables();const auto& g=m.geometry();
    const auto ne=t.edge_vertices.size(),nf=m.input().face_offsets.size()-1;
    require(initial_cut.size()==ne,"panels.cut","edge mask does not match the mesh edges");
    if(cost)require(cost->size()==ne,"panels.cost","cost does not match the mesh edges");
    std::vector<std::uint8_t> original(ne);bool user_cut=false;
    for(std::size_t e=0;e<ne;++e) {original[e]=initial_cut[e] || t.edge_boundary[e];user_cut=user_cut || (original[e] && !t.edge_boundary[e]);}
    std::size_t ngons=0,triangles=0;
    for(std::size_t f=0;f<nf;++f) {const auto n=m.input().face_offsets[f+1]-m.input().face_offsets[f];ngons+=n>4;triangles+=n==3;}
    const bool eligible=static_cast<double>(ngons)/static_cast<double>(nf)>.02 ||
                        (allow_triangles && static_cast<double>(triangles)/static_cast<double>(nf)>.9);
    if(!eligible || user_cut || deadline.expired())return std::nullopt;
    std::vector<std::uint8_t> angular(ne,1);const double coplanar=std::cos(3*(std::numbers::pi/180.0));
    for(std::size_t e=0;e<ne;++e) {
        const auto pair=t.edge_faces[e];
        if(pair.size()==2 && t.edge_consistent[e] && detail::dot(g.face_normal[pair[0]],g.face_normal[pair[1]])>coplanar)angular[e]=0;
    }
    const auto groups_owner=cut_surface(mesh,angular);const CutSurface& groups=*groups_owner;const auto& gt=groups.tables();const auto nc=gt.n_charts;
    std::vector<double> areas(nc,0.0);for(std::size_t f=0;f<nf;++f)areas[gt.chart_of_face[f]]+=g.face_area[f];
    const double total=detail::numpy_sum(areas);
    std::vector<std::uint8_t> flat(nc,0);std::vector<MeshIndex> unfenced;
    for(MeshIndex chart=0;chart<nc;++chart) {
        if(!(areas[chart]>=.01*total))continue;
        const auto faces=gt.chart_faces[chart];Vec3 normal{0,0,0};std::vector<double> weights;
        for(auto f:faces) {weights.push_back(g.face_area[f]);for(std::size_t d=0;d<3;++d)normal[d]+=g.face_normal[f][d]*g.face_area[f];}
        normal=detail::divide(normal,detail::numpy_sum(weights));
        flat[chart]=detail::norm(normal)>=.995 && !gt.chart_unsupported[chart];
        if(flat[chart] && faces.size()==1) {
            // A coarse barrel's facets are planar but only a few degrees apart: a
            // one-polygon panel needs a substantial crease fence.
            std::vector<double> lengths,sharp;
            for(auto c:m.face_corners(faces[0])) {
                const auto e=t.corner_edge[c];lengths.push_back(g.edge_length[e]);
                if(std::abs(g.edge_dihedral[e])>=45*(std::numbers::pi/180.0))sharp.push_back(g.edge_length[e]);
            }
            if(detail::numpy_sum(sharp)<.5*detail::numpy_sum(lengths))unfenced.push_back(chart);
        }
    }
    auto area_sum=[&](auto predicate) {std::vector<double> v;for(MeshIndex c=0;c<nc;++c)if(predicate(c))v.push_back(areas[c]);return detail::numpy_sum(v);};
    if(unfenced.size()>=8) {
        std::vector<double> v;for(auto c:unfenced)v.push_back(areas[c]);
        if(detail::numpy_sum(v)>.25*total)return std::nullopt;
    }
    const auto flat_count=std::count(flat.begin(),flat.end(),1);const double flat_area=area_sum([&](MeshIndex c) {return flat[c]!=0;});
    if(flat_count<4 || flat_area<.6*total)return std::nullopt;
    std::vector<std::int64_t> labels(nf);for(std::size_t f=0;f<nf;++f)labels[f]=flat[gt.chart_of_face[f]]?static_cast<std::int64_t>(gt.chart_of_face[f]):-1;
    auto boundaries=original;
    for(std::size_t e=0;e<ne;++e) {const auto pair=t.edge_faces[e];if(pair.size()==2 && labels[pair[0]]!=labels[pair[1]])boundaries[e]=1;}
    auto reassigned=reassign_bevel_bands(mesh,boundaries,std::span<const std::uint8_t>(original),true);
    if(deadline.expired())return std::nullopt;
    CadPanelCuts out;out.bevel_bands=reassigned.changed_charts.size();
    if(reassigned.changed_charts.empty()) {
        if(flat_area<.98*total)return std::nullopt;
        out.rims=boundaries;
    } else out.rims=std::move(reassigned.cut);
    out.cut.resize(ne);for(std::size_t e=0;e<ne;++e)out.cut[e]=out.rims[e] || original[e];
    auto cs=std::make_unique<CutSurface>(mesh,out.cut);std::uint64_t added=0;
    const auto rounds=static_cast<std::uint64_t>(loops_total(cs->tables()));
    for(std::uint64_t round=0;round<rounds;++round) {
        const auto& ct=cs->tables();
        if(std::all_of(ct.chart_is_disk.begin(),ct.chart_is_disk.end(),[](auto d) {return d!=0;})) {
            out.panels=ct.n_charts;out.hole_cut_edges=added;return out;
        }
        if(deadline.expired())return std::nullopt;
        MeshIndex chart=0;while(ct.chart_is_disk[chart])++chart;
        const auto path=join_boundary_loops(*cs,chart,cost);
        if(!path)return std::nullopt;
        auto candidate=out.cut;for(auto e:*path)candidate[e]=1;
        auto trial=std::make_unique<CutSurface>(mesh,candidate);
        if(trial->tables().n_charts!=ct.n_charts || loops_total(trial->tables())>=loops_total(ct))return std::nullopt;
        out.cut=std::move(candidate);cs=std::move(trial);added+=path->size();
    }
    return std::nullopt;
}
}
