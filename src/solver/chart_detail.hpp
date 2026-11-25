// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
#include "cleanseams/chart_solver.hpp"
#include "cleanseams/errors.hpp"
namespace cleanseams::chart_detail {
inline void require(bool yes,const char* field,const char* message) {
    if(!yes)throw ContractError(ErrorCode::invalid_input,field,message);
}
inline void storage(bool yes,const char* field) {
    if(!yes)throw ContractError(ErrorCode::memory_limit,field,"chart solve exceeds memory budget");
}
// Checked chart index, policy and geometry of a request.
MeshIndex checked_chart(const ChartSolveRequest&);
std::size_t local(std::span<const MeshIndex> sectors,MeshIndex global);
// Local sector constraints of this chart's corner and sector pins.
struct Constraints {
    std::vector<std::uint8_t> fixed;
    std::vector<Vec2> values;
    std::size_t count=0;
};
Constraints resolve_pins(const ChartSolveRequest&,const ChartTriangulation&);
// Conservative preflight of count-derived chart storage; returns bytes reserved.
std::uint64_t chart_storage(const ChartSolveRequest&,std::uint64_t triangles,std::uint64_t per_triangle);
// Corner output of sector UVs and, when requested, the full chart validator.
void finish(ChartSolveResult&,const ChartSolveRequest&,const ChartTriangulation&);
void validate(ChartSolveResult&,const ChartSolveRequest&);
// Source flatten_chart acceptance of an ARAP refinement: no flipped, folded,
// nonfinite or collapsed triangles and no boundary crossing in this chart.
bool refinement_acceptable(const std::shared_ptr<const QualityGeometry>&,MeshIndex chart,std::span<const MeshIndex> corners,std::span<const Vec2> uv);
// Row-major (free, 2) solves with a reusable factor, recording diagnostics/trace.
std::vector<double> factored_solve(SpdFactor&,std::vector<double> rhs,ChartSolveResult&,LaplacianTrace*);
template<class F> auto translate_allocation(const char* field,F&& body) {
    try {return body();}
    catch(const std::bad_alloc&) {throw ContractError(ErrorCode::memory_limit,field,"native chart solver allocation failed");}
}
}
