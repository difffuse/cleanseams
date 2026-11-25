// SPDX-License-Identifier: GPL-3.0-or-later
#include "chart_detail.hpp"
#include <algorithm>
#include <cmath>
namespace cleanseams {
namespace chart_detail {
MeshIndex checked_chart(const ChartSolveRequest& request) {
    require(bool(request.geometry),"solver.geometry","missing geometry snapshot");
    const auto& ct=request.geometry->cut().tables();
    require(request.chart>=0 && static_cast<std::uint64_t>(request.chart)<ct.n_charts,"solver.chart","chart index out of range");
    require(request.policy==NumericalPolicy::scipy_compatible || request.policy==NumericalPolicy::numpy_compatible,"solver.policy","unknown numerical policy");
    require(request.pin_completion==PinCompletion::legacy_010 || request.pin_completion==PinCompletion::preserve,"solver.pin_completion","unknown pin completion policy");
    require(request.weights==HarmonicWeights::mean_value || request.weights==HarmonicWeights::uniform,"solver.weights","harmonic weights must be meanvalue or uniform");
    return static_cast<MeshIndex>(request.chart);
}
std::size_t local(std::span<const MeshIndex> sectors,MeshIndex global) {
    const auto found=std::lower_bound(sectors.begin(),sectors.end(),global);
    require(found!=sectors.end() && *found==global,"solver.chart","sector missing from chart triangulation");
    return static_cast<std::size_t>(found-sectors.begin());
}
Constraints resolve_pins(const ChartSolveRequest& request,const ChartTriangulation& data) {
    const auto& ct=request.geometry->cut().tables();const auto chart=static_cast<MeshIndex>(request.chart);
    Constraints out;out.fixed.assign(data.sectors.size(),0);out.values.assign(data.sectors.size(),Vec2{0,0});
    auto add=[&](MeshIndex sector,Vec2 uv) {
        require(std::isfinite(uv[0]) && std::isfinite(uv[1]),"solver.pins","nonfinite pin coordinate");
        const auto found=std::lower_bound(data.sectors.begin(),data.sectors.end(),sector);
        if(ct.sector_chart[sector]!=chart || found==data.sectors.end() || *found!=sector)return; // source: `if s in local`
        const auto i=static_cast<std::size_t>(found-data.sectors.begin());
        if(out.fixed[i]) {
            // legacy_010 reproduces the source's per-sector dict: the last pin written wins
            // (corner pins in corner order, then sector pins), and the pipeline's final pin
            // check refuses the run with pins_changed (Gate A F7). preserve never overwrites
            // a pin: two coordinates on one sector are a constraint_conflict.
            if(request.pin_completion==PinCompletion::legacy_010) {out.values[i]=uv;return;}
            if(out.values[i]!=uv)throw ContractError(ErrorCode::constraint_conflict,"solver.pins","different pins identify one cut sector");
            return;
        }
        out.fixed[i]=1;out.values[i]=uv;++out.count;
    };
    for(const auto& pin:request.pins) {
        require(pin.corner>=0 && static_cast<std::uint64_t>(pin.corner)<ct.sector_of_corner.size(),"solver.pins","corner index out of range");
        add(ct.sector_of_corner[static_cast<std::size_t>(pin.corner)],pin.uv);
    }
    for(const auto& pin:request.sector_pins) {
        require(pin.sector<ct.n_sectors,"solver.sector_pins","sector index out of range");
        add(pin.sector,pin.uv);
    }
    return out;
}
std::uint64_t chart_storage(const ChartSolveRequest& request,std::uint64_t triangles,std::uint64_t per_triangle) {
    const auto pins=request.pins.size()+request.sector_pins.size()+request.initial_sector_uv.size();
    storage(request.memory_budget>=1024 && triangles<=(request.memory_budget-1024)/per_triangle,"solver.chart");
    const auto bound=1024+triangles*per_triangle;
    storage(pins<=(request.memory_budget-bound)/64,"solver.chart");
    return bound+pins*64;
}
void validate(ChartSolveResult& result,const ChartSolveRequest& request) {
    const std::array<MeshIndex,1> selected{static_cast<MeshIndex>(request.chart)};
    const QualityEvaluation evaluation(request.geometry,{result.corner_uv,result.corner_ids},{},std::span<const MeshIndex>(selected));
    result.valid=evaluation.chart_valid(selected[0]);result.quality=evaluation.snapshot().report(selected[0]);
}
void finish(ChartSolveResult& result,const ChartSolveRequest& request,const ChartTriangulation& data) {
    const auto& ct=request.geometry->cut().tables();
    result.sectors=data.sectors;result.corner_ids=data.corner_ids;result.corner_uv.resize(result.corner_ids.size());
    for(std::size_t i=0;i<result.corner_ids.size();++i)
        result.corner_uv[i]=result.sector_uv[local(result.sectors,ct.sector_of_corner[result.corner_ids[i]])];
    if(request.validate)validate(result,request);
}
bool refinement_acceptable(const std::shared_ptr<const QualityGeometry>& geometry,MeshIndex chart,
                           std::span<const MeshIndex> corners,std::span<const Vec2> uv) {
    const std::array<MeshIndex,1> selected{chart};const CornerUvView view{uv,corners};
    const auto report=QualitySnapshot(geometry,view,{},std::span<const MeshIndex>(selected)).report();
    if(report.flipped || report.folded || report.nonfinite || report.collapsed_uv)return false;
    return boundary_crossings(geometry->cut(),view,{},std::span<const MeshIndex>(selected)).counts[chart]==0;
}
std::vector<double> factored_solve(SpdFactor& factor,std::vector<double> rhs,ChartSolveResult& result,LaplacianTrace* trace) {
    auto solved=factor.solve(rhs,2);
    for(const auto& column:solved.columns)result.solve_diagnostics.push_back(column);
    if(!solved.columns.empty())result.diagnostics=solved.columns.front();
    if(trace) {trace->rhs.push_back(std::move(rhs));trace->solutions.push_back(solved.solution);}
    return std::move(solved.solution);
}
}
ChartTriangulation chart_triangulation(const QualityGeometry& geometry,MeshIndex chart) {
    const auto& cut=geometry.cut();const auto& ct=cut.tables();
    chart_detail::require(chart<ct.n_charts,"solver.chart","chart index out of range");
    const auto triangle_ids=geometry.chart_triangles()[chart];const auto nt=triangle_ids.size();
    ChartTriangulation data;data.triangle_ids.assign(triangle_ids.begin(),triangle_ids.end());
    data.triangle_corners.reserve(nt);data.triangles.reserve(nt);data.sectors.reserve(3*nt);data.corner_ids.reserve(3*nt);
    for(auto id:triangle_ids) {
        const auto corners=cut.mesh().triangulation().corners[id];data.triangle_corners.push_back(corners);
        for(auto c:corners) {data.sectors.push_back(ct.sector_of_corner[c]);data.corner_ids.push_back(c);}
    }
    for(auto* ids:{&data.sectors,&data.corner_ids}) {
        std::sort(ids->begin(),ids->end());ids->erase(std::unique(ids->begin(),ids->end()),ids->end());
    }
    for(const auto& corners:data.triangle_corners) {
        std::array<MeshIndex,3> tri{};
        for(std::size_t k=0;k<3;++k)tri[k]=static_cast<MeshIndex>(chart_detail::local(data.sectors,ct.sector_of_corner[corners[k]]));
        data.triangles.push_back(tri);
    }
    return data;
}
}
