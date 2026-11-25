// SPDX-License-Identifier: GPL-3.0-or-later
#include "mesh_analysis.hpp"
#include "arrays.hpp"
#include "cleanseams/mesh_analysis.hpp"
#include "cleanseams/errors.hpp"
#include <set>
namespace {
using namespace cs_bindings;
py::dict analyze(const py::dict& payload,const py::object& allowance) {
    static const std::set<std::string> known{"positions","face_offsets","corner_vertices","edge_vertices","host_edge_indices","host_loop_indices"};
    for(const auto& item:payload)if(!py::isinstance<py::str>(item.first) || !known.contains(py::cast<std::string>(item.first)))
        throw py::value_error("unknown mesh analysis input field");
    py::dict limits;limits["memory_budget"]=allowance;const auto budget=dict_uint(limits,"memory_budget","mesh analysis");
    const auto positions=dict_array(payload,"positions"),offsets=dict_array(payload,"face_offsets"),corners=dict_array(payload,"corner_vertices");
    require_array<double>(positions,2,"positions");require_width(positions,3,"positions");
    require_array<std::uint64_t>(offsets,1,"face_offsets");require_array<std::int64_t>(corners,1,"corner_vertices");
    if(offsets.size()==0)throw py::value_error("face_offsets must include initial zero");
    std::optional<py::array> edges,edge_ids,loops;
    if(payload.contains("edge_vertices")){edges=dict_array(payload,"edge_vertices");require_array<std::int64_t>(*edges,2,"edge_vertices");require_width(*edges,2,"edge_vertices");}
    if(payload.contains("host_edge_indices")){edge_ids=dict_array(payload,"host_edge_indices");require_array<std::int64_t>(*edge_ids,1,"host_edge_indices");}
    if(payload.contains("host_loop_indices")){loops=dict_array(payload,"host_loop_indices");require_array<std::int64_t>(*loops,1,"host_loop_indices");}
    const cleanseams::TableCounts counts{static_cast<std::uint64_t>(positions.shape(0)),static_cast<std::uint64_t>(offsets.size()-1),
        static_cast<std::uint64_t>(corners.size()),edges?static_cast<std::uint64_t>(edges->shape(0)):0};
    if(edge_ids && edge_ids->size()!=0 && static_cast<std::uint64_t>(edge_ids->size())!=counts.host_edges)throw py::value_error("host_edge_indices count must match edge_vertices");
    if(loops && loops->size()!=0 && static_cast<std::uint64_t>(loops->size())!=counts.corners)throw py::value_error("host_loop_indices count must match corner_vertices");
    cleanseams::validate_mesh_counts(counts.vertices,counts.faces,counts.corners,counts.host_edges,budget);
    if(cleanseams::mesh_analysis_bytes(counts)>budget)
        throw cleanseams::ContractError(cleanseams::ErrorCode::memory_limit,"mesh_analysis.memory_budget","mesh analysis and owned transfer arrays exceed the memory budget");
    cleanseams::MeshInput in;const auto p=elements<double>(positions);in.positions.reserve(static_cast<std::size_t>(counts.vertices));
    for(std::size_t i=0;i<p.size();i+=3)in.positions.push_back({p[i],p[i+1],p[i+2]});
    in.face_offsets=copy_elements<cleanseams::Offset>(offsets);in.corner_vertices=copy_elements<cleanseams::TransferIndex>(corners);
    if(edges){const auto e=elements<std::int64_t>(*edges);in.edge_vertices.reserve(e.size()/2);for(std::size_t i=0;i<e.size();i+=2)in.edge_vertices.push_back({e[i],e[i+1]});}
    if(edge_ids)in.host_edge_indices=copy_elements<cleanseams::TransferIndex>(*edge_ids);
    if(loops)in.host_loop_indices=copy_elements<cleanseams::TransferIndex>(*loops);
    cleanseams::MeshAnalysis result;
    {py::gil_scoped_release release;result=cleanseams::analyze_mesh(in,budget);}
    py::dict out;out["triangle_corners"]=owned_array(result.triangle_corners);out["triangle_faces"]=owned_array(result.triangle_faces);
    out["triangle_edges"]=owned_array(result.triangle_edges);out["face_area_3d"]=owned_array(result.face_area_3d);
    out["edge_vertices"]=owned_array(result.edge_vertices);out["corner_faces"]=owned_array(result.corner_faces);out["corner_edges"]=owned_array(result.corner_edges);
    out["host_edge_to_core"]=owned_array(result.host_edge_to_core);out["core_edge_to_host"]=owned_array(result.core_edge_to_host);
    out["corner_to_host_loop"]=owned_array(result.corner_to_host_loop);out["face_degenerate_kind"]=owned_array(result.face_degenerate_kind);
    out["triangulation_fallbacks"]=result.triangulation_fallbacks;out["geometry_key"]=result.geometry_key;return out;
}
}
void register_mesh_analysis(pybind11::module_& module) {
    module.def("analyze_mesh",&analyze,pybind11::arg("mesh_payload"),pybind11::arg("memory_budget")=pybind11::int_(cleanseams::default_memory_budget));
}
