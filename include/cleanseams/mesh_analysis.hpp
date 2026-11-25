// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
#include "cleanseams/mesh.hpp"
namespace cleanseams {
// Owned original-input IDs. Artificial triangle edges and absent host edges use -1.
struct MeshAnalysis {
    std::vector<std::array<TransferIndex,3>> triangle_corners,triangle_edges;
    std::vector<TransferIndex> triangle_faces,corner_faces,corner_edges;
    std::vector<EdgeEndpoints> edge_vertices;
    std::vector<TransferIndex> host_edge_to_core,core_edge_to_host,corner_to_host_loop;
    std::vector<double> face_area_3d;
    std::vector<std::uint8_t> face_degenerate_kind;
    std::uint64_t triangulation_fallbacks=0;
    std::string geometry_key;
};
// Conservative geometry allocation allowance: snapshot, ingestion copy and two
// owned analysis copies (native result plus transfer arrays). Saturates on overflow.
CLEANSEAMS_CORE_EXPORT std::uint64_t mesh_analysis_bytes(const TableCounts&);
CLEANSEAMS_CORE_EXPORT MeshAnalysis analyze_mesh(const MeshInput&,
    std::uint64_t memory_budget=active_memory_budget());
}
