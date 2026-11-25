// SPDX-License-Identifier: GPL-3.0-or-later
#include "cleanseams/mesh_analysis.hpp"
#include "cleanseams/detail/mesh_math.hpp"
#include "cleanseams/errors.hpp"
#include <limits>
namespace cleanseams {
std::uint64_t mesh_analysis_bytes(const TableCounts& counts) {
    constexpr auto top=std::numeric_limits<std::uint64_t>::max();
    auto bytes=mesh_table_bytes(counts);
    for(const auto& [count,stride]:{std::pair<std::uint64_t,std::uint64_t>{1,64u<<10},
        {counts.vertices,32},{counts.faces,32},{counts.corners,240},{counts.host_edges,48}}) {
        if(count>(top-bytes)/stride)return top;
        bytes+=count*stride;
    }
    return bytes;
}
MeshAnalysis analyze_mesh(const MeshInput& in,std::uint64_t budget) {
    if(in.face_offsets.empty())throw ContractError(ErrorCode::invalid_input,"mesh.face_offsets","offsets must include initial zero");
    const TableCounts counts{in.positions.size(),in.face_offsets.size()-1,in.corner_vertices.size(),in.edge_vertices.size()};
    validate_mesh_counts(counts.vertices,counts.faces,counts.corners,counts.host_edges,budget);
    const auto needed=mesh_analysis_bytes(counts);
    if(needed>budget || in.geometry_fingerprint.size()>budget-needed)
        throw ContractError(ErrorCode::memory_limit,"mesh_analysis.memory_budget","mesh analysis and owned transfer arrays exceed the memory budget");
    MemoryBudgetScope scope(budget);const MeshSnapshot mesh(in,budget);
    const auto& triangles=mesh.triangulation();const auto& tables=mesh.tables();
    MeshAnalysis out;out.geometry_key=mesh.geometry_key();out.triangulation_fallbacks=triangles.fallbacks;
    out.face_degenerate_kind=mesh.geometry().face_degenerate_kind;out.face_area_3d.assign(static_cast<std::size_t>(counts.faces),0);
    out.triangle_corners.reserve(triangles.corners.size());out.triangle_faces.assign(triangles.faces.begin(),triangles.faces.end());
    out.triangle_edges=triangles.original_edges;
    for(std::size_t t=0;t<triangles.corners.size();++t) {
        const auto corners=triangles.corners[t];out.triangle_corners.push_back({corners[0],corners[1],corners[2]});
        const auto a=in.positions[static_cast<std::size_t>(in.corner_vertices[corners[0]])];
        const auto b=in.positions[static_cast<std::size_t>(in.corner_vertices[corners[1]])];
        const auto c=in.positions[static_cast<std::size_t>(in.corner_vertices[corners[2]])];
        const auto cross=detail::cross(detail::sub(b,a),detail::sub(c,a));
        auto& area=out.face_area_3d[triangles.faces[t]];area+=.5*std::hypot(cross[0],cross[1],cross[2]);
        if(!std::isfinite(area))throw ContractError(ErrorCode::invalid_input,"mesh_analysis.face_area_3d","nonfinite source triangle area");
    }
    out.edge_vertices.reserve(tables.edge_vertices.size());
    for(const auto& edge:tables.edge_vertices)out.edge_vertices.push_back({edge[0],edge[1]});
    out.corner_faces.assign(tables.corner_face.begin(),tables.corner_face.end());out.corner_edges.assign(tables.corner_edge.begin(),tables.corner_edge.end());
    out.host_edge_to_core=tables.host_edge_to_core;out.core_edge_to_host=tables.core_edge_to_host;out.corner_to_host_loop=tables.corner_to_host_loop;
    return out;
}
}
