// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
#include "cleanseams/solver.hpp"
#include <span>
#include <ranges>
namespace cleanseams {
using MeshIndex=std::uint32_t;
struct IndexRows {
    std::vector<Offset> offsets{0};
    std::vector<MeshIndex> values;
    std::size_t size() const {return offsets.size()-1;}
    std::span<const MeshIndex> operator[](std::size_t row) const {
        const auto start=offsets.at(row),end=offsets.at(row+1);
        return std::span<const MeshIndex>(values).subspan(static_cast<std::size_t>(start),static_cast<std::size_t>(end-start));
    }
};
struct CLEANSEAMS_CORE_EXPORT Fan {
    std::vector<MeshIndex> edges;
    std::vector<MeshIndex> faces;
    bool closed=false;
    TransferIndex straight(MeshIndex edge) const;
};
struct MeshTables {
    std::vector<MeshIndex> corner_face,corner_pos,corner_next,corner_edge;
    std::vector<std::array<MeshIndex,2>> edge_vertices;
    IndexRows edge_corners,edge_faces,vert_corners,vert_edges,vert_faces;
    std::vector<std::uint8_t> edge_boundary,edge_manifold,edge_consistent;
    std::vector<std::uint8_t> vert_boundary,vert_regular;
    std::vector<std::optional<Fan>> fans;
    // The forward mapping follows supplied host edge rows. The reverse mapping
    // contains explicit host identities (or row indices if no identities supplied).
    // Wire edges have no surface edge and retain -1.
    std::vector<TransferIndex> host_edge_to_core,core_edge_to_host,corner_to_host_loop;
};
struct MeshGeometry {
    std::vector<double> edge_length,edge_dihedral,face_area,face_aspect,face_flatness;
    std::vector<Vec3> face_normal,face_center,vert_normal;
    std::vector<std::uint8_t> vert_normal_ok,face_degenerate_kind,face_nonplanar;
    std::vector<MeshIndex> face_component;
    MeshIndex n_components=0;
};
struct Triangulation {
    std::vector<std::array<MeshIndex,3>> corners;
    std::vector<MeshIndex> faces;
    // -1 denotes an artificial diagonal in this polygon, even if another
    // polygon elsewhere happens to have an edge with the same endpoints.
    std::vector<std::array<TransferIndex,3>> original_edges;
    std::uint64_t fallbacks=0;
};
CLEANSEAMS_CORE_EXPORT bool polygon_self_intersects(std::span<const Vec2>,double eps=1e-12);
// Gate A F3: byte charges of the native tables from measured allocations. Each is an
// upper bound (per-element strides 1.5x the measured envelope, every measured build covered
// at least 1.2x) of the peak live bytes of one build measured under a
// counting allocator on quad, triangle, n-gon, triangle-soup and host-row meshes, with
// sparse and all-edges cuts.
// Edges and triangles are bounded by corners, so the charges use these four counts.
struct TableCounts {std::uint64_t vertices=0,faces=0,corners=0,host_edges=0;};
CLEANSEAMS_CORE_EXPORT std::uint64_t mesh_table_bytes(const TableCounts&);    // MeshSnapshot (with its input copy)
CLEANSEAMS_CORE_EXPORT std::uint64_t cut_table_bytes(const TableCounts&);     // one CutSurface, any cut
CLEANSEAMS_CORE_EXPORT std::uint64_t quality_table_bytes(const TableCounts&); // QualityGeometry plus one full report
CLEANSEAMS_CORE_EXPORT std::uint64_t seam_table_bytes(const TableCounts&);    // one SeamGraph, any seams
// The largest set a pipeline run holds at once: mesh, cut surface, quality geometry and
// seam network. UINT64_MAX when the sum overflows.
CLEANSEAMS_CORE_EXPORT std::uint64_t run_table_bytes(const TableCounts&);
// A run's up-front refusal: memory_limit (field "memory_budget", with the bytes needed)
// when run_table_bytes exceeds `memory_budget`, before any work starts.
CLEANSEAMS_CORE_EXPORT void check_run_memory(const TableCounts&,std::uint64_t memory_budget);
// The budget the Blender host and the command line give a run: half the physical memory,
// and at least default_memory_budget. The frozen add-on had no limit of its own.
CLEANSEAMS_CORE_EXPORT std::uint64_t recommended_memory_budget();
// Checks internal representability, then that the mesh tables fit `memory_budget`
// (memory_limit otherwise); returns mesh_table_bytes.
CLEANSEAMS_CORE_EXPORT std::uint64_t validate_mesh_counts(std::uint64_t vertices,
    std::uint64_t faces,std::uint64_t corners,std::uint64_t host_edges,std::uint64_t memory_budget);
class CLEANSEAMS_CORE_EXPORT MeshSnapshot {
public:
    explicit MeshSnapshot(const MeshInput&,std::uint64_t memory_budget=active_memory_budget());
    const MeshInput& input() const {return input_;}
    const MeshTables& tables() const {return tables_;}
    const MeshGeometry& geometry() const {return geometry_;}
    const Triangulation& triangulation() const {return triangulation_;}
    // Content identity of the owned float64 geometry, distinct from a host
    // restore fingerprint and from the Session's monotonic revision ordinal.
    const std::string& geometry_key() const {return geometry_key_;}
    bool same_geometry(const MeshSnapshot&) const;
    std::ranges::iota_view<MeshIndex,MeshIndex> face_corners(MeshIndex face) const;
    MeshIndex corner_next(std::size_t corner) const;
    MeshIndex corner_prev(std::size_t corner) const;
    TransferIndex edge_index(TransferIndex a,TransferIndex b) const;
    TransferIndex face_across(MeshIndex face,MeshIndex edge) const;
    // Source Mesh.edge_other: the endpoint that is not `vertex`; for a vertex not on the
    // edge the first endpoint, as the source returns (the edge id is range-checked).
    MeshIndex edge_other(MeshIndex edge,MeshIndex vertex) const;
    Vec3 edge_dir(MeshIndex vertex,MeshIndex edge) const;
    Vec3 tangent_dir(MeshIndex vertex,MeshIndex edge) const;
    double mean_edge_length(std::optional<std::span<const MeshIndex>> edges=std::nullopt) const;
    IndexRows face_rings(TransferIndex rings) const;
    nlohmann::json diagnose() const;
    nlohmann::json summary() const;
private:
    MeshInput input_;
    MeshTables tables_;
    MeshGeometry geometry_;
    Triangulation triangulation_;
    std::string geometry_key_;
    std::uint64_t memory_budget_=0;
    std::uint64_t storage_allowance_=0;
    void build_tables();
    void build_identity();
    void build_geometry();
    void build_diagnostics();
    void build_triangulation();
    std::vector<Vec3> face_points(MeshIndex face) const;
    std::optional<Fan> build_fan(MeshIndex vertex) const;
};
// The counts the table charges are computed from.
CLEANSEAMS_CORE_EXPORT TableCounts table_counts(const MeshSnapshot&);
}
