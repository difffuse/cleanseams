// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
#include "cleanseams/mesh.hpp"
#include <memory>
#include <mutex>
#include <set>
namespace cleanseams {
// Caller-owned scratch; use one workspace per concurrent builder. Snapshots
// never retain references into it. Repeated builds on the same mesh reuse it.
class CLEANSEAMS_CORE_EXPORT CutWorkspace {
public:
    std::uint64_t capacity_bytes() const;
private:
    friend class CutSurface;
    std::vector<MeshIndex> corners_,faces_,labels_;
    std::vector<std::uint8_t> visited_;
    void prepare(std::size_t corners,std::size_t faces);
};
class CLEANSEAMS_CORE_EXPORT ChartKey {
public:
    const std::string& digest() const;
    const std::string& legacy_digest() const;
    const std::vector<MeshIndex>& faces() const {return faces_;}
    const std::vector<MeshIndex>& boundary_corners() const {return boundary_;}
    bool operator==(const ChartKey&) const;
private:
    friend class CutSurface;
    ChartKey(std::shared_ptr<const MeshSnapshot>,std::span<const MeshIndex> faces,std::vector<MeshIndex> boundary);
    // GB performance: the digests are computed on first use (most charts of a surface are
    // never keyed), once, thread-safely, and shared by copies of the key. The values are
    // exactly those the constructor computed before.
    struct Digests {std::once_flag once,legacy_once;std::string digest,legacy;};
    std::shared_ptr<const MeshSnapshot> mesh_;
    std::vector<MeshIndex> faces_,boundary_;
    std::shared_ptr<Digests> digests_;
};
struct ChartAdjacency {
    std::array<MeshIndex,2> charts{};
    std::vector<MeshIndex> edges;
    IndexRows arcs;
};
struct CutTables {
    std::vector<std::uint8_t> cut;
    std::vector<MeshIndex> sector_of_corner,chart_of_face,chart_of_corner,sector_vert,sector_chart;
    MeshIndex n_sectors=0,n_charts=0;
    std::vector<MeshIndex> boundary_corners,pinched_sectors,loop_chart;
    IndexRows loops,chart_loops,chart_faces;
    std::vector<std::int64_t> chart_euler;
    std::vector<Offset> chart_boundary_sides,chart_n_loops,cut_degree;
    std::vector<std::uint8_t> chart_pinched,chart_unsupported,chart_is_disk;
    // First encountered edge order, as in the source's insertion-ordered dict.
    std::vector<ChartAdjacency> adjacency;
};
class CLEANSEAMS_CORE_EXPORT CutSurface {
public:
    CutSurface(std::shared_ptr<const MeshSnapshot>,std::span<const std::uint8_t> cut,
        std::uint64_t memory_budget=active_memory_budget());
    CutSurface(std::shared_ptr<const MeshSnapshot>,std::span<const std::uint8_t> cut,CutWorkspace&,
        std::uint64_t memory_budget=active_memory_budget());
    const MeshSnapshot& mesh() const {return *mesh_;}
    const std::shared_ptr<const MeshSnapshot>& mesh_owner() const {return mesh_;}
    const CutTables& tables() const {return data_;}
    const ChartKey& chart_key(MeshIndex chart) const {return keys_.at(chart);}
    std::vector<MeshIndex> loop_sectors(MeshIndex loop) const;
    std::vector<MeshIndex> loop_verts(MeshIndex loop) const;
    std::vector<MeshIndex> sectors_at_vertex(MeshIndex vertex) const;
    std::set<std::array<MeshIndex,2>> adjacency_keys() const;
    nlohmann::json chart_summary(MeshIndex chart) const;
    // Sum of resulting Euler characteristics over charts touched by these
    // unique, removable effective cuts; empty input returns zero as in source.
    std::int64_t euler_after_gluing(std::span<const MeshIndex> edges) const;
private:
    std::shared_ptr<const MeshSnapshot> mesh_;
    CutTables data_;
    std::vector<ChartKey> keys_;
    std::uint64_t budget_=0,storage_allowance_=0;
    void build(std::span<const std::uint8_t>,CutWorkspace&,std::uint64_t);
    void build_boundaries(CutWorkspace&);
    void build_chart_keys();
};
// Source arc_key: the node vertices (cut degree other than two) the edges touch and
// CPython's hash(frozenset(edges)). Ordered like the source's (ends, hash) tuple; the
// sorted edge set breaks a hash tie so that distinct edge sets never compare equal.
struct ArcKey {
    std::vector<MeshIndex> ends,edges;
    std::int64_t hash=0;
    bool operator==(const ArcKey&) const=default;
    std::strong_ordering operator<=>(const ArcKey& other) const {
        if(auto c=ends<=>other.ends;c!=0)return c;
        if(auto c=hash<=>other.hash;c!=0)return c;
        return edges<=>other.edges;
    }
    nlohmann::json to_json() const {return nlohmann::json::array({ends,hash});}
};
CLEANSEAMS_CORE_EXPORT IndexRows split_arcs(const MeshSnapshot&,std::span<const MeshIndex> edges,std::span<const Offset> cut_degree);
CLEANSEAMS_CORE_EXPORT ArcKey arc_key(const CutSurface&,std::span<const MeshIndex> edges);
CLEANSEAMS_CORE_EXPORT Offset free_ends(const CutSurface&,std::span<const MeshIndex> edges);
}
