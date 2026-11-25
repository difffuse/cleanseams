// SPDX-License-Identifier: GPL-3.0-or-later
#include "detail.hpp"
#include <map>
namespace cleanseams {
std::int64_t CutSurface::euler_after_gluing(std::span<const MeshIndex> edges) const {
    using namespace cut_detail;const auto& m=mesh_->tables();const auto& input=mesh_->input();
    require(edges.size()<=m.edge_vertices.size(),"glue.edges","too many edges");
    if(edges.size()>(budget_-storage_allowance_)/256)throw ContractError(ErrorCode::memory_limit,"glue.memory_budget","gluing workspace exceeds memory allowance");
    std::set<MeshIndex> unique,charts;std::map<MeshIndex,MeshIndex> parent;
    // Only touched sectors participate; no full-mesh parent array is allocated.
    auto find=[&](MeshIndex value) {
        parent.try_emplace(value,value);auto root=value;
        while(parent.at(root)!=root) root=parent.at(root);
        while(parent.at(value)!=root) {const auto next=parent.at(value);parent.at(value)=root;value=next;}
        return root;
    };
    std::int64_t unions=0;
    for(auto edge:edges) {
        require(edge<m.edge_vertices.size() && unique.insert(edge).second,"glue.edges","invalid or repeated edge");
        // Two-sided effective cuts, as the source (cut_topology.euler_after_gluing). An
        // inconsistently wound edge stays forced into the cut on a rebuild, but the source
        // still predicts its gluing (the pairing below handles either orientation) and the
        // layout callers rely on that; refusing it failed whole runs (T25).
        require(data_.cut[edge] && !m.edge_boundary[edge] && m.edge_manifold[edge] && m.edge_corners[edge].size()==2,
            "glue.edges","edge is not a two-sided effective cut");
        const auto corners=m.edge_corners[edge];const auto a=corners[0],b=corners[1];
        const auto an=mesh_->corner_next(a),bn=mesh_->corner_next(b);
        const auto sa=data_.sector_of_corner[a],sb=data_.sector_of_corner[b],san=data_.sector_of_corner[an],sbn=data_.sector_of_corner[bn];
        const auto pairs=input.corner_vertices[a]==input.corner_vertices[bn]?
            std::array<std::array<MeshIndex,2>,2>{{{sa,sbn},{san,sb}}}:
            std::array<std::array<MeshIndex,2>,2>{{{sa,sb},{san,sbn}}};
        for(const auto& pair:pairs) {
            const auto x=find(pair[0]),y=find(pair[1]);if(x!=y) {parent.at(std::max(x,y))=std::min(x,y);++unions;}
        }
        charts.insert(data_.chart_of_corner[a]);charts.insert(data_.chart_of_corner[b]);
    }
    std::int64_t result=static_cast<std::int64_t>(edges.size())-unions;
    for(auto chart:charts) result+=data_.chart_euler[chart];return result;
}
}
