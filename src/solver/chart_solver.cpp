// SPDX-License-Identifier: GPL-3.0-or-later
#include "chart_detail.hpp"
#include "cleanseams/cache.hpp"
#include "cleanseams/control.hpp"
#include "cleanseams/detail/mesh_math.hpp"
#include "cleanseams/detail/numpy_compat.hpp"
#include <algorithm>
#include <cmath>
#include <map>
namespace cleanseams {
namespace {
using namespace chart_detail;
void refuse_pins(const ChartSolveRequest& request,MeshIndex chart,const char* field,const char* message) {
    if(resolve_pins(request,chart_triangulation(*request.geometry,chart)).count)
        throw ContractError(ErrorCode::constraint_conflict,field,message);
}
std::optional<ChartSolveResult> planar_impl(const ChartSolveRequest& request) {
    const auto chart=checked_chart(request);const auto& geometry=*request.geometry;
    const auto& cut=geometry.cut();const auto& ct=cut.tables();const auto& mesh=cut.mesh();const auto& g=mesh.geometry();
    refuse_pins(request,chart,"solver.planar","the planar n-gon map cannot honor in-chart pins");
    const auto faces=ct.chart_faces[chart];
    const auto& offsets=mesh.input().face_offsets;
    const bool ngon=std::any_of(faces.begin(),faces.end(),[&](auto f) {return offsets[f+1]-offsets[f]>4;});
    if(faces.empty() || !ngon || !ct.chart_is_disk[chart])return std::nullopt;
    Vec3 normal{0,0,0};std::vector<double> areas;areas.reserve(faces.size());
    for(auto f:faces) {
        for(std::size_t d=0;d<3;++d)normal[d]+=g.face_normal[f][d]*g.face_area[f];
        areas.push_back(g.face_area[f]);
    }
    const double norm=detail::norm(normal);
    if(norm<.95*detail::numpy_sum(areas) || norm<=0)return std::nullopt;
    normal=detail::divide(normal,norm);
    std::size_t axis=0;for(std::size_t d=1;d<3;++d)if(std::abs(normal[d])<std::abs(normal[axis]))axis=d;
    Vec3 unit{0,0,0};unit[axis]=1;
    auto u=detail::cross(normal,unit);u=detail::divide(u,detail::norm(u));const auto v=detail::cross(normal,u);
    std::vector<MeshIndex> ids;
    for(auto f:faces)for(auto c=offsets[f];c<offsets[f+1];++c)ids.push_back(static_cast<MeshIndex>(c));
    std::vector<Vec3> points;points.reserve(ids.size());Vec3 mean{0,0,0};
    for(auto c:ids) {
        points.push_back(mesh.input().positions[static_cast<std::size_t>(mesh.input().corner_vertices[c])]);
        for(std::size_t d=0;d<3;++d)mean[d]+=points.back()[d];
    }
    for(auto& m:mean)m/=double(ids.size());
    std::vector<Vec2> values;values.reserve(ids.size());
    for(auto& p:points) {p=detail::sub(p,mean);values.push_back({detail::dot(p,u),detail::dot(p,v)});}
    const std::array<MeshIndex,1> selected{chart};const CornerUvView view{values,ids};
    const auto measured=QualitySnapshot(request.geometry,view,{},std::span<const MeshIndex>(selected)).report();
    if(measured.flipped || measured.folded || measured.nonfinite || measured.collapsed_uv ||
       measured.anisotropy_p95>1.08 || measured.anisotropy_max>4 || measured.area_p95>1.1 ||
       boundary_crossings(cut,view,{},std::span<const MeshIndex>(selected)).counts[chart])return std::nullopt;
    ChartSolveResult result;result.method="planar";result.corner_ids=ids;result.corner_uv=values;
    result.numerics.policy=request.policy;result.numerics.version="orthographic-projection-v1";
    std::vector<std::pair<MeshIndex,std::size_t>> first;
    for(std::size_t i=0;i<ids.size();++i)first.push_back({ct.sector_of_corner[ids[i]],i});
    std::stable_sort(first.begin(),first.end(),[](const auto& a,const auto& b) {return a.first<b.first;});
    for(std::size_t i=0;i<first.size();++i)if(!i || first[i].first!=first[i-1].first) {
        result.sectors.push_back(first[i].first);result.sector_uv.push_back(values[first[i].second]);
    }
    if(request.validate)validate(result,request);
    return result;
}
ChartSolveResult flatten_uncached(const ChartSolveRequest& request,ChartMethod method);
// T22: a chart map is reused only under its complete numerical contract.
std::string solve_contract(const ChartSolveRequest& r,ChartMethod method,MeshIndex chart) {
    const auto& cs=r.geometry->cut();const auto& ct=cs.tables();
    nlohmann::json key{{"geometry",cs.mesh().geometry_key()},{"chart",cs.chart_key(chart).digest()},{"method",static_cast<int>(method)},
        {"iterations",r.iterations},{"weights",static_cast<int>(r.weights)},{"policy",static_cast<int>(r.policy)},
        {"completion",static_cast<int>(r.pin_completion)},{"validate",r.validate},{"version","cleanseams.solver.v1"}};
    nlohmann::json pins=nlohmann::json::array();
    for(const auto& p:r.pins)if(p.corner>=0 && static_cast<std::size_t>(p.corner)<ct.chart_of_corner.size() && ct.chart_of_corner[static_cast<std::size_t>(p.corner)]==chart)
        pins.push_back({p.corner,p.uv[0],p.uv[1]});
    nlohmann::json sectors=nlohmann::json::array();
    for(const auto& p:r.sector_pins)if(p.sector<ct.sector_chart.size() && ct.sector_chart[p.sector]==chart)sectors.push_back({p.sector,p.uv[0],p.uv[1]});
    key["pins"]=pins;key["sector_pins"]=sectors;key["start"]=r.initial_sector_uv; // exact coordinates, no rounding
    return key.dump();
}
// A representative corner (the first in corner order) of each global sector in `sectors`.
// Sector ids are global to a cut surface, while a solve contract keys the chart only, so a
// cache hit relabels them through these corners for the surface that asks (Gate A F31).
std::vector<MeshIndex> sector_corners(const CutSurface& cs,std::span<const MeshIndex> corner_ids,std::span<const MeshIndex> sectors) {
    std::map<MeshIndex,MeshIndex> first;for(auto c:corner_ids)first.emplace(cs.tables().sector_of_corner[c],c);
    std::vector<MeshIndex> out;out.reserve(sectors.size());
    for(auto s:sectors) {const auto found=first.find(s);require(found!=first.end(),"solver.cache","sector without a chart corner");out.push_back(found->second);}
    return out;
}
std::vector<MeshIndex> relabel(const CutSurface& cs,std::span<const MeshIndex> corners) {
    std::vector<MeshIndex> out;out.reserve(corners.size());for(auto c:corners)out.push_back(cs.tables().sector_of_corner[c]);return out;
}
ChartSolveResult flatten_impl(const ChartSolveRequest& request,ChartMethod method) {
    auto* cache=active_cache();
    if(!cache || request.retain_assembly)return flatten_uncached(request,method);
    const auto chart=checked_chart(request);const auto contract=solve_contract(request,method,chart);
    const auto& cs=request.geometry->cut();
    if(auto hit=cache->solved(contract)) {
        ChartSolveResult r;r.method=std::move(hit->method);r.corner_ids=std::move(hit->corner_ids);
        r.sectors=relabel(cs,hit->sector_corners);r.automatic_pins=relabel(cs,hit->pin_corners);
        r.corner_uv=std::move(hit->corner_uv);r.sector_uv=std::move(hit->sector_uv);
        r.numerics=std::move(hit->numerics);r.diagnostics=hit->diagnostics;r.solve_diagnostics=std::move(hit->solve_diagnostics);
        r.iterations=hit->iterations;r.equation_relative_residual=hit->equation_relative_residual;r.valid=hit->valid;r.quality=std::move(hit->quality);
        return r;
    }
    auto result=flatten_uncached(request,method);
    cache->store(contract,{result.method,result.corner_ids,result.corner_uv,result.sector_uv,result.numerics,
                           result.diagnostics,result.solve_diagnostics,result.iterations,result.equation_relative_residual,result.valid,result.quality,
                           sector_corners(cs,result.corner_ids,result.sectors),sector_corners(cs,result.corner_ids,result.automatic_pins)});
    return result;
}
ChartSolveResult flatten_uncached(const ChartSolveRequest& request,ChartMethod method) {
    detail::note_chart_solve();
    const auto chart=checked_chart(request);
    require(method==ChartMethod::lscm || method==ChartMethod::harmonic || method==ChartMethod::arap,"solver.method","unknown chart method");
    require(request.initial_sector_uv.empty(),"solver.initial_sector_uv","flatten_chart chooses its own start map");
    if(method==ChartMethod::lscm)return solve_lscm(request);
    if(method==ChartMethod::harmonic)return solve_harmonic(request);
    if(!resolve_pins(request,chart_triangulation(*request.geometry,chart)).count) {
        if(auto planar=planar_impl(request)) {
            // Refine a valid geometric map; keep it if the solve folds or crosses.
            ChartSolveRequest refine=request;refine.initial_sector_uv=planar->sector_uv;
            require(planar->sectors==chart_triangulation(*request.geometry,chart).sectors,"solver.planar","planar sectors differ from the chart triangulation");
            auto candidate=solve_arap(refine);
            if(refinement_acceptable(request.geometry,chart,candidate.corner_ids,candidate.corner_uv))return candidate;
            return std::move(*planar);
        }
    }
    return solve_arap(request);
}
FlattenAllResult flatten_all_impl(const FlattenAllRequest& request) {
    require(bool(request.geometry),"solver.geometry","missing geometry snapshot");
    const auto& ct=request.geometry->cut().tables();const auto nc=ct.sector_of_corner.size();
    require(request.uv.empty() || request.uv.size()==nc,"solver.uv","UV array does not match the mesh corners");
    FlattenAllResult out;out.uv=request.uv.empty()?std::vector<Vec2>(nc,Vec2{0,0}):request.uv;
    std::vector<MeshIndex> charts;
    if(request.charts)charts=*request.charts;
    else for(MeshIndex ch=0;ch<ct.n_charts;++ch)charts.push_back(ch);
    std::vector<MeshIndex> solve;
    for(auto ch:charts) {
        require(ch<ct.n_charts,"solver.charts","chart index out of range");
        if(!ct.chart_faces[ch].empty())solve.push_back(ch);
    }
    // Independent charts of one immutable geometry: each solve owns its factor and
    // workspace; results are committed in the requested chart order (T24).
    std::vector<std::optional<ChartSolveResult>> slots(solve.size());
    parallel_for(solve.size(),[&](std::size_t k) {
        ChartSolveRequest chart;chart.geometry=request.geometry;chart.chart=solve[k];chart.pins=request.pins;
        chart.policy=request.policy;chart.pin_completion=request.pin_completion;chart.iterations=request.iterations;
        chart.validate=request.validate;chart.memory_budget=request.memory_budget;
        slots[k]=flatten_impl(chart,request.method);
    },&request.cancel);
    throw_if_cancelled(request.cancel);
    for(std::size_t k=0;k<solve.size();++k) {
        auto& result=*slots[k];
        for(std::size_t i=0;i<result.corner_ids.size();++i)out.uv[result.corner_ids[i]]=result.corner_uv[i];
        out.solved_charts.push_back(solve[k]);out.charts.push_back(std::move(result));
    }
    return out;
}
}
std::optional<ChartSolveResult> planar_ngon_map(const ChartSolveRequest& request) {
    return translate_allocation("solver.planar",[&] {return planar_impl(request);});
}
ChartSolveResult flatten_chart(const ChartSolveRequest& request,ChartMethod method) {
    return translate_allocation("solver.flatten_chart",[&] {return flatten_impl(request,method);});
}
FlattenAllResult flatten_all(const FlattenAllRequest& request) {
    return translate_allocation("solver.flatten_all",[&] {return flatten_all_impl(request);});
}
NativeArapBackend::NativeArapBackend(std::shared_ptr<const MeshSnapshot> mesh,Settings settings)
    :mesh_(std::move(mesh)),settings_(settings) {
    require(bool(mesh_),"backend.mesh","missing mesh snapshot");
}
NativeArapBackend::Unwrap NativeArapBackend::unwrap(std::span<const std::uint8_t> cut,std::optional<std::span<const MeshIndex>> charts,
    std::span<const CornerPin> pins,std::optional<std::uint64_t> iterations,std::optional<std::span<const Vec2>> init) const {
    return translate_allocation("backend.unwrap",[&] {
        Unwrap out;out.cut=cut_surface(mesh_,cut,settings_.memory_budget);
        out.geometry=quality_geometry(out.cut,settings_.memory_budget);
        const auto& ct=out.cut->tables();const auto nc=ct.sector_of_corner.size();
        const auto steps=iterations.value_or(settings_.iterations);
        if(init)require(init->size()==nc,"backend.init","start map does not match the mesh corners");
        std::vector<MeshIndex> selected;
        if(charts)selected.assign(charts->begin(),charts->end());
        else for(MeshIndex ch=0;ch<ct.n_charts;++ch)selected.push_back(ch);
        if(!init && steps==10) { // source: flatten_all(method="arap") at its default ten iterations
            FlattenAllRequest request;request.geometry=out.geometry;request.charts=selected;
            request.pins.assign(pins.begin(),pins.end());request.policy=settings_.policy;
            request.pin_completion=settings_.pin_completion;request.validate=settings_.validate;
            request.memory_budget=settings_.memory_budget;request.cancel=settings_.cancel;
            auto all=flatten_all(request);
            out.uv=std::move(all.uv);out.solved_charts=std::move(all.solved_charts);out.charts=std::move(all.charts);
            return out;
        }
        out.uv.assign(nc,Vec2{0,0});
        std::vector<MeshIndex> solve;
        for(auto ch:selected) {
            require(ch<ct.n_charts,"backend.charts","chart index out of range");
            if(!ct.chart_faces[ch].empty())solve.push_back(ch);
        }
        std::vector<std::optional<ChartSolveResult>> slots(solve.size());
        parallel_for(solve.size(),[&](std::size_t k) {
            const auto ch=solve[k];
            ChartSolveRequest request;request.geometry=out.geometry;request.chart=ch;request.pins.assign(pins.begin(),pins.end());
            request.policy=settings_.policy;request.pin_completion=settings_.pin_completion;request.iterations=steps;
            request.validate=settings_.validate;request.memory_budget=settings_.memory_budget;
            if(init) {
                // Source: last corner of each sector in face order supplies its start;
                // any nonfinite start falls back to the LSCM initialization.
                const auto sectors=chart_triangulation(*out.geometry,ch).sectors;std::vector<Vec2> start(sectors.size(),Vec2{0,0});
                for(auto f:ct.chart_faces[ch])for(auto c:mesh_->face_corners(f))start[local(sectors,ct.sector_of_corner[c])]=(*init)[c];
                if(std::all_of(start.begin(),start.end(),[](const Vec2& p) {return std::isfinite(p[0]) && std::isfinite(p[1]);}))
                    request.initial_sector_uv=std::move(start);
            }
            slots[k]=solve_arap(request);
        },&settings_.cancel);
        throw_if_cancelled(settings_.cancel);
        for(std::size_t k=0;k<solve.size();++k) { // committed in chart order (T24)
            auto& result=*slots[k];
            for(std::size_t i=0;i<result.corner_ids.size();++i)out.uv[result.corner_ids[i]]=result.corner_uv[i];
            out.solved_charts.push_back(solve[k]);out.charts.push_back(std::move(result));
        }
        return out;
    });
}
}
