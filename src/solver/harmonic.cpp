// SPDX-License-Identifier: GPL-3.0-or-later
#include "chart_detail.hpp"
#include "cleanseams/detail/mesh_math.hpp"
#include "cleanseams/detail/numpy_compat.hpp"
#include <algorithm>
#include <cmath>
#include <numbers>
#include <unordered_map>
namespace cleanseams {
namespace {
using namespace chart_detail;
using Pairs=std::vector<std::array<MeshIndex,2>>;
// Source mean_value_weights: at each triangle corner, tan(angle/2)/|edge| is
// added to both incident edges; keys keep dict first-insertion order.
void mean_value_weights(const QualityGeometry& geometry,const ChartTriangulation& data,Pairs& pairs,std::vector<double>& weights) {
    const auto& ct=geometry.cut().tables();const auto& positions=geometry.cut().mesh().input().positions;
    std::unordered_map<std::uint64_t,std::size_t> index;
    auto add=[&](MeshIndex a,MeshIndex b,double w) {
        const std::array<MeshIndex,2> key=a<b?std::array<MeshIndex,2>{a,b}:std::array<MeshIndex,2>{b,a};
        const auto [at,inserted]=index.try_emplace((std::uint64_t(key[0])<<32)|key[1],pairs.size());
        if(inserted) {pairs.push_back(key);weights.push_back(0.0);}
        weights[at->second]+=w;
    };
    for(const auto& s:data.triangles) {
        std::array<Vec3,3> p;for(std::size_t k=0;k<3;++k)p[k]=positions[ct.sector_vert[data.sectors[s[k]]]];
        for(std::size_t k=0;k<3;++k) {
            const auto a=detail::sub(p[(k+1)%3],p[k]),b=detail::sub(p[(k+2)%3],p[k]);
            const double la=detail::norm(a),lb=detail::norm(b);
            if(la<1e-14 || lb<1e-14)continue;
            const double c=std::max(-1.0,std::min(1.0,detail::dot(a,b)/(la*lb)));
            const double t2=std::tan(std::min(std::acos(c),std::numbers::pi-1e-6)/2.0);
            add(s[k],s[(k+1)%3],t2/la);add(s[k],s[(k+2)%3],t2/lb);
        }
    }
}
// Source uniform weights: unique sorted triangle edges, unit weight.
void uniform_weights(const ChartTriangulation& data,Pairs& pairs,std::vector<double>& weights) {
    for(const auto [a,b]:{std::pair<std::size_t,std::size_t>{0,1},{1,2},{2,0}})
        for(const auto& t:data.triangles)pairs.push_back({std::min(t[a],t[b]),std::max(t[a],t[b])});
    std::sort(pairs.begin(),pairs.end());pairs.erase(std::unique(pairs.begin(),pairs.end()),pairs.end());
    pairs.erase(std::remove_if(pairs.begin(),pairs.end(),[](const auto& p) {return p[0]==p[1];}),pairs.end());
    weights.assign(pairs.size(),1.0);
}
ChartSolveResult harmonic_impl(const ChartSolveRequest& request,const char* method) {
    const auto chart=checked_chart(request);const auto& geometry=*request.geometry;
    const auto nt=geometry.chart_triangles()[chart].size();const auto reserved=chart_storage(request,nt,1024);
    const auto data=chart_triangulation(geometry,chart);const auto ns=data.sectors.size();
    Pairs pairs;std::vector<double> weights;
    if(request.weights==HarmonicWeights::uniform)uniform_weights(data,pairs,weights);
    else mean_value_weights(geometry,data,pairs,weights);
    std::vector<TransferIndex> rows,cols;std::vector<double> values,diagonal(ns,0.0);
    rows.reserve(2*pairs.size()+ns);cols.reserve(2*pairs.size()+ns);values.reserve(2*pairs.size()+ns);
    for(std::size_t e=0;e<pairs.size();++e) {
        const auto [i,j]=pairs[e];const double w=weights[e];
        rows.insert(rows.end(),{i,j});cols.insert(cols.end(),{j,i});values.insert(values.end(),{-w,-w});
        diagonal[i]+=w;diagonal[j]+=w;
    }
    for(std::size_t i=0;i<ns;++i) {rows.push_back(TransferIndex(i));cols.push_back(TransferIndex(i));values.push_back(diagonal[i]);}
    const auto pins=resolve_pins(request,data);
    std::vector<TransferIndex> free_index(ns,-1);TransferIndex nfree=0;
    for(std::size_t i=0;i<ns;++i)if(!pins.fixed[i])free_index[i]=nfree++;
    ChartSolveResult result;result.method=method;result.sector_uv=pins.values;
    std::unique_ptr<LaplacianTrace> trace;
    if(request.retain_assembly) {
        trace=std::make_unique<LaplacianTrace>();trace->weight_pairs=pairs;trace->weights=weights;
        trace->matrix=std::make_unique<SparseMatrix>(rows,cols,values,std::array<Offset,2>{ns,ns},request.memory_budget);
        trace->fixed=pins.fixed;trace->free_index=free_index;
    }
    if(nfree) {
        // b = -sum over free rows of L[r, fixed c] * xfix[c], accumulated in COO order.
        std::vector<double> rhs(2*static_cast<std::size_t>(nfree),0.0);std::vector<TransferIndex> rr,rc;std::vector<double> rv;
        for(std::size_t k=0;k<values.size();++k) {
            const auto r=static_cast<std::size_t>(rows[k]),c=static_cast<std::size_t>(cols[k]);
            if(free_index[r]<0)continue;
            const auto row=static_cast<std::size_t>(free_index[r]);
            if(pins.fixed[c]) {rhs[2*row]+=values[k]*pins.values[c][0];rhs[2*row+1]+=values[k]*pins.values[c][1];}
            else {rr.push_back(free_index[r]);rc.push_back(free_index[c]);rv.push_back(values[k]);}
        }
        for(auto& b:rhs)b=-b;
        const SparseMatrix reduced(rr,rc,rv,{Offset(nfree),Offset(nfree)},request.memory_budget);
        if(trace)trace->reduced=std::make_unique<SparseMatrix>(reduced);
        FactorOptions options;options.policy=request.policy;options.memory_budget=request.memory_budget-reserved;
        SpdFactor factor(reduced,options);
        const auto x=factored_solve(factor,std::move(rhs),result,trace.get());
        for(std::size_t i=0;i<ns;++i)if(free_index[i]>=0) {
            const auto f=static_cast<std::size_t>(free_index[i]);result.sector_uv[i]={x[2*f],x[2*f+1]};
        }
        result.numerics=factor.metadata();
    } else {result.numerics.policy=request.policy;result.numerics.version="empty-v1";}
    for(const auto& uv:result.sector_uv)if(!std::isfinite(uv[0]) || !std::isfinite(uv[1]))
        throw ContractError(ErrorCode::numerical_failure,"solver.harmonic","nonfinite harmonic solution");
    finish(result,request,data);result.laplacian=std::move(trace);
    return result;
}
std::optional<ChartSolveResult> convex_impl(const ChartSolveRequest& request) {
    const auto chart=checked_chart(request);const auto& geometry=*request.geometry;
    const auto& cut=geometry.cut();const auto& ct=cut.tables();const auto& mesh=cut.mesh();
    if(resolve_pins(request,chart_triangulation(geometry,chart)).count)
        throw ContractError(ErrorCode::constraint_conflict,"solver.convex","the convex boundary map cannot honor in-chart pins");
    if(ct.chart_n_loops[chart]!=1 || ct.chart_euler[chart]!=1 || ct.chart_pinched[chart])return std::nullopt;
    const auto loop=ct.loops[ct.chart_loops[chart][0]];const auto n=loop.size();
    std::vector<MeshIndex> sectors;for(auto c:loop)sectors.push_back(ct.sector_of_corner[c]);
    auto unique=sectors;std::sort(unique.begin(),unique.end());
    if(n<3 || std::unique(unique.begin(),unique.end())!=unique.end())return std::nullopt;
    const auto& positions=mesh.input().positions;std::vector<double> lengths(n);
    for(std::size_t i=0;i<n;++i) {
        const auto& p=positions[static_cast<std::size_t>(mesh.input().corner_vertices[loop[i]])];
        const auto& q=positions[static_cast<std::size_t>(mesh.input().corner_vertices[loop[(i+1)%n]])];
        lengths[i]=detail::norm(detail::sub(q,p));
    }
    const double floor=std::max(detail::numpy_sum(lengths)/double(n)*.01,1e-15);
    for(auto& l:lengths)l=std::max(l,floor);
    const double total=detail::numpy_sum(lengths);
    std::vector<Vec2> boundary(n);double arc=0;
    for(std::size_t i=0;i<n;++i) {
        if(i)arc+=lengths[i-1];
        const double angle=2*std::numbers::pi*arc/total;boundary[i]={std::cos(angle),std::sin(angle)};
    }
    std::vector<double> terms(n);
    for(std::size_t i=0;i<n;++i) {const auto& b=boundary[i];const auto& c=boundary[(i+1)%n];terms[i]=b[0]*c[1]-b[1]*c[0];}
    const double polygon_area=.5*std::abs(detail::numpy_sum(terms));
    std::vector<double> areas;for(auto f:ct.chart_faces[chart])areas.push_back(mesh.geometry().face_area[f]);
    const double scale=std::sqrt(detail::numpy_sum(areas)/std::max(polygon_area,1e-300));
    ChartSolveRequest harmonic=request;harmonic.pins.clear();harmonic.sector_pins.clear();harmonic.weights=HarmonicWeights::uniform;
    for(std::size_t i=0;i<n;++i)harmonic.sector_pins.push_back({sectors[i],{boundary[i][0]*scale,boundary[i][1]*scale}});
    return harmonic_impl(harmonic,"convex");
}
}
ChartSolveResult solve_harmonic(const ChartSolveRequest& request) {
    return translate_allocation("solver.harmonic",[&] {return harmonic_impl(request,"harmonic");});
}
std::optional<ChartSolveResult> convex_boundary_map(const ChartSolveRequest& request) {
    return translate_allocation("solver.convex",[&] {return convex_impl(request);});
}
}
