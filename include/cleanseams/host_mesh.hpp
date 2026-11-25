// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
#include "cleanseams/solver.hpp"
#include <span>
namespace cleanseams {
// A synchronous view of Blender foreach_get buffers, never stored by the core.
// Fingerprinting uses the float32/int32 raw host order before any transform.
struct HostMeshView {
    std::span<const float> positions; // flat xyz triples
    std::span<const std::int32_t> edge_vertices; // flat endpoint pairs, wires included
    std::span<const std::int32_t> loop_vertices,loop_edges,face_loop_start,face_loop_total;
};
struct ObjectTransform {
    std::array<double,16> object_to_world{1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1};
    Vec3 object_scale{1,1,1}; // obj.scale, not a decomposition of matrix_world
    bool world_space=false;
};
struct HostMeshInput {
    MeshInput mesh;
    bool nonuniform_scale=false;
};
class MeshSnapshot;
// A painted attribute as read from the host (T29): values in raw host order of its domain.
struct HostPaintedAttribute {
    std::string name,domain,data_type;
    std::span<const double> values;
};
struct HostAttributes {
    // Optional masks in raw host edge/vertex/loop order; empty means absent.
    std::span<const std::uint8_t> seams,locked,anchors,pins;
    std::optional<std::span<const Vec2>> uv;
};
struct MappedAttributes {
    std::vector<std::uint8_t> seams,locked,pins;
    std::vector<TransferIndex> anchors;
    std::optional<std::vector<Vec2>> uv;
    std::vector<CornerPin> pin_constraints;
    // Painted attributes in core order of their domain (POINT, EDGE, FACE, CORNER); other
    // domains keep their name, domain and type with no values (reported by the resolver).
    std::vector<PaintedAttribute> painted;
};
CLEANSEAMS_CORE_EXPORT MappedAttributes read_host_attributes(const MeshSnapshot&,
    const HostAttributes&,std::uint64_t memory_budget=active_memory_budget(),std::span<const HostPaintedAttribute> painted={});
// Compatible with the source adapter's little-endian Blender fingerprint.
// This restore token is not the float64 native geometry content key.
CLEANSEAMS_CORE_EXPORT std::string host_geometry_fingerprint(const HostMeshView&,
    std::uint64_t memory_budget=active_memory_budget());
CLEANSEAMS_CORE_EXPORT HostMeshInput read_host_mesh(const HostMeshView&,
    const ObjectTransform& transform={},std::uint64_t memory_budget=active_memory_budget());
}
