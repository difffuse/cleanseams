// SPDX-License-Identifier: GPL-3.0-or-later
#include "cleanseams/seam_graph.hpp"
#include "cleanseams/cache.hpp"
#include "cleanseams/detail/numpy_compat.hpp"
#include "detail.hpp"
namespace cleanseams {
using namespace cut_detail;
double Chain::edge_length(const MeshSnapshot& mesh) const {
    // Source: float(mesh.edge_length[edges].sum()), a NumPy pairwise sum.
    std::vector<double> lengths;lengths.reserve(edges.size());
    for(auto edge:edges) lengths.push_back(mesh.geometry().edge_length.at(edge));
    return detail::numpy_sum(lengths);
}
SeamGraph::SeamGraph(std::shared_ptr<const MeshSnapshot> mesh,std::span<const std::uint8_t> seams,
    std::span<const std::uint8_t> locked,std::span<const TransferIndex> anchors,
    std::shared_ptr<const CutSurface> cut,std::span<const TransferIndex> pins,std::uint64_t budget) {
    require(bool(mesh),"seams.mesh","mesh owner is required");
    const auto& m=mesh->tables();const auto& input=mesh->input();
    const auto nv=input.positions.size(),ne=m.edge_vertices.size(),nc=input.corner_vertices.size(),nf=input.face_offsets.size()-1;
    require(seams.size()==ne,"seams.mask","edge count mismatch");
    require(locked.empty() || locked.size()==ne,"seams.locked","edge count mismatch");
    for(auto bit:seams) require(bit<=1,"seams.mask","mask must contain zero or one");
    for(auto bit:locked) require(bit<=1,"seams.locked","mask must contain zero or one");
    require(anchors.size()<=budget/sizeof(TransferIndex) && pins.size()<=budget/sizeof(TransferIndex),
        "seams.constraints","constraint input exceeds budget");
    for(auto v:anchors) require(v>=0 && static_cast<std::uint64_t>(v)<nv,"seams.anchors","vertex out of range");
    for(auto c:pins) require(c>=0 && static_cast<std::uint64_t>(c)<nc,"seams.pins","corner out of range");
    // The mesh tables plus the measured charge of seam tables, per-vertex side provenance,
    // chains (every edge a seam at worst) and traversal scratch (Gate A F3); a cut surface
    // built here checks its own tables.
    const auto allowance=validate_mesh_counts(nv,nf,nc,input.edge_vertices.size(),budget);
    const auto own=seam_table_bytes({nv,nf,nc,input.edge_vertices.size()});
    if(own>budget-allowance)throw ContractError(ErrorCode::memory_limit,"seams.memory_budget","seam network exceeds the memory budget");
    if(cut) {
        require(mesh==cut->mesh_owner() || mesh->same_geometry(cut->mesh()),"seams.cut","cut geometry mismatch");
        for(std::size_t e=0;e<ne;++e) {
            const bool effective=seams[e] || m.edge_boundary[e] || !m.edge_manifold[e] || !m.edge_consistent[e];
            require(bool(cut->tables().cut[e])==effective,"seams.cut","cut mask mismatch");
        }
        cut_=std::move(cut);
    } else cut_=cut_surface(mesh,seams,budget);
    auto& d=data_;d.seams.assign(seams.begin(),seams.end());d.locked.assign(ne,0);
    if(!locked.empty()) std::copy(locked.begin(),locked.end(),d.locked.begin());
    d.provenance.resize(nv);
    for(auto v:anchors) d.provenance[static_cast<std::size_t>(v)].artist_fixed=true;
    for(MeshIndex v=0;v<nv;++v) {
        auto& p=d.provenance[v];p.sectors=cut_->sectors_at_vertex(v);
        p.physical_boundary=bool(m.vert_boundary[v]);p.unsupported_fan=!m.fans[v];
        if(p.artist_fixed) d.manual_anchors.push_back(v);
    }
    // Duplicate input corners must not grow retained storage with input length.
    if(!pins.empty()) {
        std::vector<std::uint8_t> pinned(cut_->tables().n_sectors,0);
        for(auto c:pins) pinned[cut_->tables().sector_of_corner[static_cast<std::size_t>(c)]]=1;
        for(auto& p:d.provenance) {
            for(auto s:p.sectors) if(pinned[s]) p.pinned_sectors.push_back(s);
            p.pinned_sector=!p.pinned_sectors.empty();
        }
    }
    build();
}
void SeamGraph::build() {
    auto& d=data_;const auto& mesh=cut_->mesh();const auto& m=mesh.tables();
    const auto nv=mesh.input().positions.size(),ne=d.seams.size();
    d.degree.assign(nv,0);d.edge_charts.resize(ne);d.anchor.assign(nv,0);
    d.vert_seam_edges=rows(nv,[&](auto add){for(std::size_t e=0;e<ne;++e) if(d.seams[e])
        for(auto v:m.edge_vertices[e]) add(v,static_cast<MeshIndex>(e));});
    for(std::size_t e=0;e<ne;++e) if(d.seams[e]) {
        for(auto v:m.edge_vertices[e]) {++d.degree[v];if(d.locked[e]) d.provenance[v].locked_incident=true;}
        const auto faces=m.edge_faces[e];
        if(faces.size()==2) {
            auto a=cut_->tables().chart_of_face[faces[0]],b=cut_->tables().chart_of_face[faces[1]];
            if(a>b) std::swap(a,b);d.edge_charts[e]=std::array<MeshIndex,2>{a,b};
        }
    }
    for(std::size_t v=0;v<nv;++v) {
        auto& p=d.provenance[v];const auto degree=d.degree[v];
        p.generated_junction=degree>2;p.seam_end=degree==1;
        if(degree==2) {const auto edges=d.vert_seam_edges[v];p.flanking_change=d.edge_charts[edges[0]]!=d.edge_charts[edges[1]];}
        d.anchor[v]=p.artist_fixed || p.locked_incident || p.flanking_change ||
            (degree>0 && (degree!=2 || p.physical_boundary || p.unsupported_fan));
    }
    std::vector<std::uint8_t> used(ne,0);
    auto walk=[&](MeshIndex vertex,MeshIndex edge,bool remaining_cycle) {
        Chain chain;chain.id=static_cast<MeshIndex>(chains_.size());chain.verts.push_back(vertex);
        while(true) {
            used[edge]=1;chain.edges.push_back(edge);vertex=mesh.edge_other(edge,vertex);chain.verts.push_back(vertex);
            if(d.anchor[vertex]) break;
            MeshIndex next=0;std::size_t count=0;
            for(auto e:d.vert_seam_edges[vertex]) if(e!=edge) {next=e;++count;}
            if(count!=1 || used[next]) {chain.cycle=chain.verts.front()==vertex && !d.anchor[vertex];break;}
            edge=next;
        }
        if(remaining_cycle) chain.cycle=true;
        const auto pair=d.edge_charts[chain.edges[0]];bool same_pair=true;
        for(auto e:chain.edges) {chain.locked=chain.locked || d.locked[e];same_pair=same_pair && d.edge_charts[e]==pair;}
        if(same_pair && pair) chain.charts=std::vector<MeshIndex>{(*pair)[0],(*pair)[1]};
        if(chain.cycle) chain.reason="closed_cycle";
        else if(chain.locked) chain.reason="locked";
        else if(!chain.charts) chain.reason="mixed_flanking_charts";
        else if(chain.edges.size()<2) chain.reason="too_short";
        if(!chain.cycle) for(auto v:{chain.verts.front(),chain.verts.back()})
            chain.free_end=chain.free_end || (d.degree[v]==1 && !m.vert_boundary[v] && bool(m.fans[v]));
        chains_.push_back(std::move(chain));
    };
    for(MeshIndex v=0;v<nv;++v) if(d.anchor[v]) for(auto edge:d.vert_seam_edges[v]) if(!used[edge]) walk(v,edge,false);
    for(MeshIndex e=0;e<ne;++e) if(d.seams[e] && !used[e]) walk(m.edge_vertices[e][0],e,true);
}
std::vector<MeshIndex> SeamGraph::editable_chains() const {
    std::vector<MeshIndex> result;for(const auto& chain:chains_) if(chain.reason.empty()) result.push_back(chain.id);return result;
}
nlohmann::json SeamGraph::summary() const {
    std::size_t cycles=0,locked=0;for(const auto& chain:chains_) {cycles+=chain.cycle;locked+=chain.locked;}
    return {{"seam_edges",std::count(data_.seams.begin(),data_.seams.end(),1)},
        {"anchors",std::count(data_.anchor.begin(),data_.anchor.end(),1)},
        {"chains",chains_.size()},{"editable",editable_chains().size()},
        {"cycles",cycles},{"locked",locked},{"charts",cut_->tables().n_charts}};
}
}
