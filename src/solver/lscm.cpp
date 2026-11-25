// SPDX-License-Identifier: GPL-3.0-or-later
#include "chart_detail.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <new>
namespace cleanseams {
namespace {
using chart_detail::require;
using chart_detail::local;
std::vector<MeshIndex> pin_candidates(const QualityGeometry& geometry,MeshIndex chart,const ChartTriangulation& data) {
    const auto& ct=geometry.cut().tables();std::vector<MeshIndex> candidates;
    for(auto loop:ct.chart_loops[chart])for(auto corner:ct.loops[loop])candidates.push_back(ct.sector_of_corner[corner]);
    if(candidates.empty())candidates=data.sectors;
    std::sort(candidates.begin(),candidates.end());candidates.erase(std::unique(candidates.begin(),candidates.end()),candidates.end());
    require(!candidates.empty(),"solver.chart","chart has no triangulated sectors");
    return candidates;
}
// First argmax of the 3D distance from `from`'s vertex, as np.argmax(norm(P-p)).
std::pair<std::size_t,double> farthest(const QualityGeometry& geometry,std::span<const MeshIndex> candidates,MeshIndex from) {
    const auto& ct=geometry.cut().tables();const auto& positions=geometry.cut().mesh().input().positions;
    const auto& p=positions[ct.sector_vert[from]];double best=-1;std::size_t which=0;
    for(std::size_t i=0;i<candidates.size();++i) {
        const auto& q=positions[ct.sector_vert[candidates[i]]];double d=0;
        for(std::size_t axis=0;axis<3;++axis)d+=(q[axis]-p[axis])*(q[axis]-p[axis]);
        d=std::sqrt(d);require(std::isfinite(d),"solver.pins","nonfinite pin-selection distance");
        if(d>best){best=d;which=i;}
    }
    return {which,best};
}
void fix(const ChartTriangulation& data,MeshIndex sector,Vec2 uv,std::vector<std::uint8_t>& fixed,std::vector<double>& values,
         std::vector<MeshIndex>& automatic) {
    const auto ns=data.sectors.size(),i=local(data.sectors,sector);
    fixed[i]=fixed[ns+i]=1;values[i]=uv[0];values[ns+i]=uv[1];automatic.push_back(sector);
}
// Source _default_pins: two farthest-point sweeps from the first candidate.
void default_pins(const QualityGeometry& geometry,MeshIndex chart,const ChartTriangulation& data,
                  std::vector<std::uint8_t>& fixed,std::vector<double>& values,std::vector<MeshIndex>& automatic) {
    const auto candidates=pin_candidates(geometry,chart,data);
    const auto first=farthest(geometry,candidates,candidates[0]).first,second=farthest(geometry,candidates,candidates[first]).first;
    fix(data,candidates[first],{0,0},fixed,values,automatic);fix(data,candidates[second],{1,0},fixed,values,automatic);
}
// PIN-01 correction: keep the user's pin and add the farthest distinct anchor at
// the same V, offset by their 3D distance (1 if coincident). Boundary candidates
// are preferred; other chart sectors are used only if no distinct one remains.
void complete_single_pin(const QualityGeometry& geometry,MeshIndex chart,const ChartTriangulation& data,MeshIndex pinned,
                         std::vector<std::uint8_t>& fixed,std::vector<double>& values,std::vector<MeshIndex>& automatic) {
    auto candidates=pin_candidates(geometry,chart,data);
    const auto without=[&](std::vector<MeshIndex> ids) {ids.erase(std::remove(ids.begin(),ids.end(),pinned),ids.end());return ids;};
    candidates=without(std::move(candidates));
    if(candidates.empty())candidates=without(data.sectors);
    if(candidates.empty())return;
    const auto [which,distance]=farthest(geometry,candidates,pinned);const auto i=local(data.sectors,pinned),ns=data.sectors.size();
    fix(data,candidates[which],{values[i]+(distance>0?distance:1.0),values[ns+i]},fixed,values,automatic);
}
LscmAssembly assemble_impl(const ChartSolveRequest& request) {
    const auto chart=chart_detail::checked_chart(request);const auto& geometry=*request.geometry;
    const auto triangle_ids=geometry.chart_triangles()[chart];const auto nt=triangle_ids.size();
    // Includes COO construction/copies, reduced arrays, chart indexing and pins.
    const auto owned_bound=chart_detail::chart_storage(request,nt,2048);
    auto data=chart_triangulation(geometry,chart);
    const auto ns=data.sectors.size();const auto pins=chart_detail::resolve_pins(request,data);
    std::vector<double> fixed_values(2*ns,0);std::vector<std::uint8_t> fixed(2*ns,0);
    for(std::size_t i=0;i<ns;++i)if(pins.fixed[i]) {fixed[i]=fixed[ns+i]=1;fixed_values[i]=pins.values[i][0];fixed_values[ns+i]=pins.values[i][1];}
    auto user_pins=pins.count;std::vector<MeshIndex> automatic;
    if(user_pins==1 && request.pin_completion==PinCompletion::legacy_010) {
        // Source flatten.py:183 discards a lone pin (PIN-01); reproduced only on request.
        std::fill(fixed.begin(),fixed.end(),0);std::fill(fixed_values.begin(),fixed_values.end(),0);user_pins=0;
    }
    if(user_pins==0)default_pins(geometry,chart,data,fixed,fixed_values,automatic);
    else if(user_pins==1) {
        std::size_t i=0;while(!fixed[i])++i;
        complete_single_pin(geometry,chart,data,data.sectors[i],fixed,fixed_values,automatic);
    }
    std::vector<TransferIndex> rows,columns;std::vector<double> values;
    rows.reserve(nt*12);columns.reserve(nt*12);values.reserve(nt*12);
    auto add=[&](std::size_t row,std::size_t column,double value) {
        rows.push_back(static_cast<TransferIndex>(row));columns.push_back(static_cast<TransferIndex>(column));values.push_back(value);
    };
    for(std::size_t t=0;t<nt;++t) {
        const auto& frame=geometry.frames()[triangle_ids[t]];if(frame.double_area<=1e-14)continue;
        const double scale=1/std::sqrt(frame.double_area);const auto& q=frame.coordinates;
        for(std::size_t k=0;k<3;++k) {
            const auto j=(k+1)%3,l=(k+2)%3;const double wr=(q[l][0]-q[j][0])*scale,wi=(q[l][1]-q[j][1])*scale;
            const auto u=static_cast<std::size_t>(data.triangles[t][k]),v=ns+u;
            add(2*t,u,wr);add(2*t,v,-wi);add(2*t+1,u,wi);add(2*t+1,v,wr);
        }
    }
    std::vector<TransferIndex> free_index(2*ns,-1);TransferIndex nfree=0;
    for(std::size_t i=0;i<2*ns;++i)if(!fixed[i])free_index[i]=nfree++;
    std::vector<double> rhs(2*nt,0),reduced_values;std::vector<TransferIndex> reduced_rows,reduced_columns;
    reduced_values.reserve(values.size());reduced_rows.reserve(rows.size());reduced_columns.reserve(columns.size());
    for(std::size_t k=0;k<values.size();++k) {
        const auto c=static_cast<std::size_t>(columns[k]),r=static_cast<std::size_t>(rows[k]);
        if(fixed[c])rhs[r]-=values[k]*fixed_values[c];
        else {reduced_rows.push_back(rows[k]);reduced_columns.push_back(free_index[c]);reduced_values.push_back(values[k]);}
    }
    for(double b:rhs)require(std::isfinite(b),"solver.lscm.rhs","nonfinite reduced RHS");
    return {std::move(data),SparseMatrix(rows,columns,values,{2*nt,2*ns},request.memory_budget),
        SparseMatrix(reduced_rows,reduced_columns,reduced_values,{2*nt,static_cast<Offset>(nfree)},request.memory_budget),
        std::move(fixed_values),std::move(rhs),std::move(free_index),std::move(automatic),owned_bound};
}
ChartSolveResult solve_impl(const ChartSolveRequest& request) {
    auto assembly=std::make_unique<LscmAssembly>(assemble_impl(request));
    auto linear=solve_least_squares(assembly->reduced,assembly->rhs,request.policy,
        request.memory_budget-assembly->storage_bound,request.retain_assembly);
    auto& x=assembly->fixed_values;
    // Preserve the assembly's fixed vector for trace comparison.
    auto solution=x;for(std::size_t i=0;i<solution.size();++i)if(assembly->free_index[i]>=0)solution[i]=linear.solution[static_cast<std::size_t>(assembly->free_index[i])];
    ChartSolveResult result;result.method="lscm";result.automatic_pins=assembly->automatic_pins;
    const auto ns=assembly->chart.sectors.size();result.sector_uv.resize(ns);
    for(std::size_t i=0;i<ns;++i)result.sector_uv[i]={solution[i],solution[ns+i]};
    result.numerics=std::move(linear.numerics);result.diagnostics=linear.diagnostics;result.solve_diagnostics={linear.diagnostics};
    result.equation_relative_residual=linear.equation_relative_residual;
    result.normal=std::move(linear.normal);result.normal_rhs=std::move(linear.normal_rhs);
    chart_detail::finish(result,request,assembly->chart);
    if(request.retain_assembly)result.assembly=std::move(assembly);
    return result;
}
}
// Public entry points translate allocation failure into the structured limit.
LscmAssembly assemble_lscm(const ChartSolveRequest& request) {
    return chart_detail::translate_allocation("solver.lscm",[&] {return assemble_impl(request);});
}
ChartSolveResult solve_lscm(const ChartSolveRequest& request) {
    return chart_detail::translate_allocation("solver.lscm",[&] {return solve_impl(request);});
}
}
