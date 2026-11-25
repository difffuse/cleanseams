// SPDX-License-Identifier: GPL-3.0-or-later
#include "detail.hpp"
#include "cleanseams/cache.hpp"
#include <map>
#include <unordered_map>
namespace cleanseams {
using namespace cut_detail;
std::uint64_t CutWorkspace::capacity_bytes() const {
    return (corners_.capacity()+faces_.capacity()+labels_.capacity())*sizeof(MeshIndex)+visited_.capacity();
}
void CutWorkspace::prepare(std::size_t nc,std::size_t nf) {
    // Release buffers when switching mesh sizes; same-mesh searches retain
    // capacity. This prevents a prior large mesh retaining an unbudgeted cache.
    if(corners_.size()!=nc || faces_.size()!=nf) {
        std::vector<MeshIndex>().swap(corners_);std::vector<MeshIndex>().swap(faces_);
        std::vector<MeshIndex>().swap(labels_);std::vector<std::uint8_t>().swap(visited_);
    }
    corners_.resize(nc);faces_.resize(nf);labels_.reserve(std::max(nc,nf));visited_.assign(nc,0);
    std::iota(corners_.begin(),corners_.end(),MeshIndex{0});std::iota(faces_.begin(),faces_.end(),MeshIndex{0});
}
CutSurface::CutSurface(std::shared_ptr<const MeshSnapshot> mesh,std::span<const std::uint8_t> cut,std::uint64_t budget):mesh_(std::move(mesh)) {
    CutWorkspace workspace;build(cut,workspace,budget);detail::note_cut_build(*mesh_,data_.cut);
}
CutSurface::CutSurface(std::shared_ptr<const MeshSnapshot> mesh,std::span<const std::uint8_t> cut,CutWorkspace& workspace,std::uint64_t budget):mesh_(std::move(mesh)) {
    build(cut,workspace,budget);detail::note_cut_build(*mesh_,data_.cut);
}
void CutSurface::build(std::span<const std::uint8_t> cut,CutWorkspace& workspace,std::uint64_t budget) {
    require(bool(mesh_),"cut.mesh","mesh owner is required");
    const auto& in=mesh_->input();const auto& m=mesh_->tables();auto& d=data_;
    const auto nc=in.corner_vertices.size(),nf=in.face_offsets.size()-1,nv=in.positions.size(),ne=m.edge_vertices.size();
    require(cut.size()==ne,"cut.mask","cut count must match canonical mesh edges");
    for(auto bit:cut) require(bit<=1,"cut.mask","cut mask must contain zero or one");
    // The mesh tables plus the measured charge of cut tables, keys, ordered adjacency and
    // arcs, validation maps and constructor scratch for any cut (Gate A F3).
    std::uint64_t allowance=validate_mesh_counts(nv,nf,nc,in.edge_vertices.size(),budget);
    const auto own=cut_table_bytes({nv,nf,nc,in.edge_vertices.size()});
    if(own>budget-allowance)throw ContractError(ErrorCode::memory_limit,"cut.memory_budget","cut tables exceed the memory budget");
    allowance+=own;
    budget_=budget;storage_allowance_=allowance;
    workspace.prepare(nc,nf);d.cut.assign(cut.begin(),cut.end());
    for(std::size_t e=0;e<ne;++e) {
        d.cut[e]=d.cut[e] || m.edge_boundary[e] || !m.edge_manifold[e] || !m.edge_consistent[e];
        if(d.cut[e]) continue;
        const auto corners=m.edge_corners[e];if(corners.size()!=2) continue;
        const auto a=corners[0],b=corners[1],an=mesh_->corner_next(a),bn=mesh_->corner_next(b);
        if(in.corner_vertices[a]==in.corner_vertices[bn]) {unite(workspace.corners_,a,bn);unite(workspace.corners_,an,b);}
        else {unite(workspace.corners_,a,b);unite(workspace.corners_,an,bn);}
        unite(workspace.faces_,m.corner_face[a],m.corner_face[b]);
    }
    d.n_sectors=labels(workspace.corners_,d.sector_of_corner,workspace.labels_);
    d.n_charts=labels(workspace.faces_,d.chart_of_face,workspace.labels_);
    d.chart_of_corner.resize(nc);d.sector_vert.resize(d.n_sectors);d.sector_chart.resize(d.n_sectors);
    for(std::size_t c=0;c<nc;++c) {
        d.chart_of_corner[c]=d.chart_of_face[m.corner_face[c]];
        d.sector_vert[d.sector_of_corner[c]]=static_cast<MeshIndex>(in.corner_vertices[c]);
        d.sector_chart[d.sector_of_corner[c]]=d.chart_of_corner[c];
        if(d.cut[m.corner_edge[c]]) d.boundary_corners.push_back(static_cast<MeshIndex>(c));
    }
    build_boundaries(workspace);
    d.chart_faces=rows(d.n_charts,[&](auto add){for(std::size_t f=0;f<nf;++f) add(d.chart_of_face[f],static_cast<MeshIndex>(f));});
    d.chart_euler.assign(d.n_charts,0);d.chart_boundary_sides.assign(d.n_charts,0);
    d.chart_n_loops.resize(d.n_charts);d.chart_pinched.assign(d.n_charts,0);d.chart_unsupported.assign(d.n_charts,0);d.chart_is_disk.resize(d.n_charts);
    for(auto ch:d.sector_chart) ++d.chart_euler[ch];
    for(auto ch:d.chart_of_face) ++d.chart_euler[ch];
    for(auto c:d.boundary_corners) {const auto ch=d.chart_of_corner[c];--d.chart_euler[ch];++d.chart_boundary_sides[ch];}
    for(std::size_t e=0;e<ne;++e) if(!d.cut[e] && m.edge_faces[e].size()==2) --d.chart_euler[d.chart_of_face[m.edge_faces[e][0]]];
    for(auto sector:d.pinched_sectors) d.chart_pinched[d.sector_chart[sector]]=1;
    for(std::size_t v=0;v<nv;++v) if(!m.fans[v]) for(auto f:m.vert_faces[v]) d.chart_unsupported[d.chart_of_face[f]]=1;
    for(MeshIndex ch=0;ch<d.n_charts;++ch) {
        d.chart_n_loops[ch]=d.chart_loops[ch].size();
        d.chart_is_disk[ch]=d.chart_n_loops[ch]==1 && d.chart_euler[ch]==1 && !d.chart_pinched[ch] && !d.chart_unsupported[ch];
    }
    // Adjacencies in first-encountered edge order; the hash map only finds a pair's entry
    // (GB performance: no ordered-map nodes per build).
    d.cut_degree.assign(nv,0);std::unordered_map<std::uint64_t,std::size_t> lookup;
    for(std::size_t e=0;e<ne;++e) if(d.cut[e]) {
        for(auto v:m.edge_vertices[e]) ++d.cut_degree[v];
        const auto faces=m.edge_faces[e];if(m.edge_boundary[e] || faces.size()!=2) continue;
        auto a=d.chart_of_face[faces[0]],b=d.chart_of_face[faces[1]];if(a>b) std::swap(a,b);
        const std::array<MeshIndex,2> key{a,b};
        const auto [found,inserted]=lookup.try_emplace((static_cast<std::uint64_t>(a)<<32)|b,d.adjacency.size());
        if(inserted) d.adjacency.push_back({key,{},{}});
        d.adjacency[found->second].edges.push_back(static_cast<MeshIndex>(e));
    }
    for(auto& adjacent:d.adjacency) adjacent.arcs=split_arcs(*mesh_,adjacent.edges,d.cut_degree);
    build_chart_keys();
}
std::vector<MeshIndex> CutSurface::loop_sectors(MeshIndex loop) const {
    std::vector<MeshIndex> result;for(auto c:data_.loops[loop]) result.push_back(data_.sector_of_corner[c]);return result;
}
std::vector<MeshIndex> CutSurface::loop_verts(MeshIndex loop) const {
    std::vector<MeshIndex> result;for(auto c:data_.loops[loop]) result.push_back(static_cast<MeshIndex>(mesh_->input().corner_vertices[c]));return result;
}
std::vector<MeshIndex> CutSurface::sectors_at_vertex(MeshIndex vertex) const {
    std::vector<MeshIndex> result;for(auto c:mesh_->tables().vert_corners[vertex]) result.push_back(data_.sector_of_corner[c]);
    std::sort(result.begin(),result.end());result.erase(std::unique(result.begin(),result.end()),result.end());return result;
}
std::set<std::array<MeshIndex,2>> CutSurface::adjacency_keys() const {
    std::set<std::array<MeshIndex,2>> result;for(const auto& item:data_.adjacency) result.insert(item.charts);return result;
}
nlohmann::json CutSurface::chart_summary(MeshIndex chart) const {
    return {{"faces",data_.chart_faces[chart].size()},{"loops",data_.chart_n_loops.at(chart)},
        {"euler",data_.chart_euler.at(chart)},{"pinched",bool(data_.chart_pinched.at(chart))},
        {"unsupported",bool(data_.chart_unsupported.at(chart))},{"disk",bool(data_.chart_is_disk.at(chart))}};
}
}
