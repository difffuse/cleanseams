// SPDX-License-Identifier: GPL-3.0-or-later
#include "chart_detail.hpp"
#include "cleanseams/control.hpp"
#include <algorithm>
#include <array>
#include <exception>
#include <cmath>
namespace cleanseams {
namespace {
using namespace chart_detail;
// Source arap(): cotangent Laplacian clamped to [-5, 50], LSCM (or supplied)
// start, user pins reapplied, otherwise local sector 0 fixed; one factorization
// of the free system and a two-column solve per local/global iteration.
ChartSolveResult arap_impl(const ChartSolveRequest& request) {
    const auto chart=checked_chart(request);const auto& geometry=*request.geometry;
    const auto nt=geometry.chart_triangles()[chart].size();const auto reserved=chart_storage(request,nt,2048);
    const auto data=chart_triangulation(geometry,chart);const auto ns=data.sectors.size();
    ChartSolveResult result;result.method="arap";
    if(nt==0) {result.numerics.policy=request.policy;result.numerics.version="empty-v1";return result;}
    std::vector<std::array<Vec2,3>> q(nt);std::vector<std::array<double,3>> cot(nt);
    for(std::size_t t=0;t<nt;++t) {
        const auto& frame=geometry.frames()[data.triangle_ids[t]];q[t]=frame.coordinates;
        const bool good=frame.double_area>1e-14;
        for(std::size_t k=0;k<3;++k) {
            const auto& o=q[t][k];const auto& p=q[t][(k+1)%3];const auto& r=q[t][(k+2)%3];
            const Vec2 a{p[0]-o[0],p[1]-o[1]},b{r[0]-o[0],r[1]-o[1]};
            const double cross=a[0]*b[1]-a[1]*b[0],dot=a[0]*b[0]+a[1]*b[1];
            double value=std::abs(cross)>1e-14?dot/cross:0.0;if(!good)value=0.0;
            cot[t][k]=std::min(std::max(value,-5.0),50.0);
        }
    }
    std::vector<TransferIndex> rows,cols;std::vector<double> values;const auto entries=12*nt;
    rows.reserve(entries);cols.reserve(entries);values.reserve(entries);
    for(std::size_t k=0;k<3;++k) {
        const auto ik=(k+1)%3,jk=(k+2)%3;
        for(std::size_t part=0;part<4;++part)for(std::size_t t=0;t<nt;++t) {
            const TransferIndex i=data.triangles[t][ik],j=data.triangles[t][jk];const double w=cot[t][k];
            rows.push_back(part==0 || part==2?i:j);cols.push_back(part==0 || part==3?j:part==1?i:part==2?i:j);
            values.push_back(part<2?-w:w);
        }
    }
    std::vector<Vec2> suv;
    const bool lscm_start=request.initial_sector_uv.empty();
    ChartSolveRequest start;
    if(lscm_start) {start=request;start.validate=false;start.retain_assembly=false;start.initial_sector_uv.clear();}
    else {
        require(request.initial_sector_uv.size()==ns,"solver.initial_sector_uv","initial map does not match the chart's sectors");
        suv=request.initial_sector_uv;
    }
    std::unique_ptr<LaplacianTrace> trace;
    // The free system: user pins, otherwise local sector 0 fixed (removing the translation
    // freedom), the reduced Laplacian and its factor. It does not depend on the start map.
    Constraints pins;std::vector<TransferIndex> free_index;TransferIndex nfree=0;std::vector<std::size_t> fixed_entries;
    std::unique_ptr<SpdFactor> factor;
    auto free_system=[&] {
        pins=resolve_pins(request,data);
        auto fixed=pins.fixed;if(!pins.count)fixed[0]=1;
        free_index.assign(ns,-1);
        for(std::size_t i=0;i<ns;++i)if(!fixed[i])free_index[i]=nfree++;
        std::vector<TransferIndex> rr,rc;std::vector<double> rv;
        for(std::size_t e=0;e<values.size();++e) {
            const auto r=static_cast<std::size_t>(rows[e]),c=static_cast<std::size_t>(cols[e]);
            if(free_index[r]<0)continue;
            if(free_index[c]>=0) {rr.push_back(free_index[r]);rc.push_back(free_index[c]);rv.push_back(values[e]);}
            else fixed_entries.push_back(e);
        }
        const SparseMatrix reduced(rr,rc,rv,{Offset(nfree),Offset(nfree)},request.memory_budget);
        if(trace) {trace->fixed=fixed;trace->free_index=free_index;trace->reduced=std::make_unique<SparseMatrix>(reduced);}
        FactorOptions options;options.policy=request.policy;options.memory_budget=request.memory_budget-reserved;
        factor=std::make_unique<SpdFactor>(reduced,options); // once, reused by every iteration
    };
    const auto* pool=active_workers();
    if(lscm_start && !request.retain_assembly && pool && pool->threads()>1) {
        // GB performance: the LSCM start and the free system's factorization are independent;
        // idle workers take one while this thread computes the other (both are exactly the
        // serial computations). The LSCM failure, which the serial order reaches first, wins.
        std::vector<Vec2> lscm;std::array<std::exception_ptr,2> failure;
        parallel_for_nested(2,[&](std::size_t i) {
            try {if(i==0)lscm=solve_lscm(start).sector_uv;else free_system();}
            catch(...) {failure[i]=std::current_exception();}
        });
        for(const auto& error:failure)if(error)std::rethrow_exception(error);
        suv=std::move(lscm);
    } else {
        if(lscm_start)suv=solve_lscm(start).sector_uv;
        if(request.retain_assembly) {
            trace=std::make_unique<LaplacianTrace>();trace->cotangents=cot;trace->initial_sector_uv=suv;
            trace->matrix=std::make_unique<SparseMatrix>(rows,cols,values,std::array<Offset,2>{ns,ns},request.memory_budget);
        }
        free_system();
    }
    for(std::size_t i=0;i<ns;++i)if(pins.fixed[i])suv[i]=pins.values[i];
    if(!pins.count)result.automatic_pins={data.sectors[0]};
    auto& factor_ref=*factor;
    std::vector<Vec2> rhs(ns);std::vector<std::array<Vec2,3>> rotated(nt); // every entry rewritten per iteration
    for(std::uint64_t it=0;it<request.iterations;++it) {
        std::fill(rhs.begin(),rhs.end(),Vec2{0,0});
        for(std::size_t t=0;t<nt;++t) {
            const auto& tri=data.triangles[t];std::array<Vec2,3> e,eu;
            for(std::size_t k=0;k<3;++k) {
                const auto a=(k+1)%3,b=(k+2)%3;
                e[k]={q[t][b][0]-q[t][a][0],q[t][b][1]-q[t][a][1]};
                eu[k]={suv[tri[b]][0]-suv[tri[a]][0],suv[tri[b]][1]-suv[tri[a]][1]};
            }
            // S = sum_k cot_k Eu_k E_k^T, einsum order ((cot*Eu)*E) summed over k.
            double s[2][2]={{0,0},{0,0}};
            for(std::size_t i=0;i<2;++i)for(std::size_t j=0;j<2;++j)
                for(std::size_t k=0;k<3;++k)s[i][j]+=cot[t][k]*eu[k][i]*e[k][j];
            const double a=s[0][0]+s[1][1],b=s[1][0]-s[0][1],n=std::sqrt(a*a+b*b);
            const double c=n>1e-14?a/std::max(n,1e-300):1.0,sn=n>1e-14?b/std::max(n,1e-300):0.0;
            for(std::size_t k=0;k<3;++k)rotated[t][k]={c*e[k][0]+(-sn)*e[k][1],sn*e[k][0]+c*e[k][1]};
        }
        for(std::size_t k=0;k<3;++k) {
            const auto ik=(k+1)%3,jk=(k+2)%3;
            for(std::size_t t=0;t<nt;++t)for(std::size_t d=0;d<2;++d)rhs[data.triangles[t][jk]][d]+=cot[t][k]*rotated[t][k][d];
            for(std::size_t t=0;t<nt;++t)for(std::size_t d=0;d<2;++d)rhs[data.triangles[t][ik]][d]+=-cot[t][k]*rotated[t][k][d];
        }
        if(nfree==0)break;
        std::vector<double> bfree(2*static_cast<std::size_t>(nfree)),fixed_sum(2*static_cast<std::size_t>(nfree),0.0);
        for(std::size_t i=0;i<ns;++i)if(free_index[i]>=0) {
            const auto f=static_cast<std::size_t>(free_index[i]);bfree[2*f]=rhs[i][0];bfree[2*f+1]=rhs[i][1];
        }
        for(auto e:fixed_entries) {
            const auto f=static_cast<std::size_t>(free_index[static_cast<std::size_t>(rows[e])]);const auto& x=suv[static_cast<std::size_t>(cols[e])];
            fixed_sum[2*f]+=values[e]*x[0];fixed_sum[2*f+1]+=values[e]*x[1];
        }
        if(!fixed_entries.empty())for(std::size_t i=0;i<bfree.size();++i)bfree[i]-=fixed_sum[i];
        const auto x=factored_solve(factor_ref,std::move(bfree),result,trace.get());
        for(std::size_t i=0;i<ns;++i)if(free_index[i]>=0) {
            const auto f=static_cast<std::size_t>(free_index[i]);suv[i]={x[2*f],x[2*f+1]};
        }
        ++result.iterations;
    }
    for(const auto& uv:suv)if(!std::isfinite(uv[0]) || !std::isfinite(uv[1]))
        throw ContractError(ErrorCode::numerical_failure,"solver.arap","nonfinite ARAP solution");
    result.sector_uv=std::move(suv);result.numerics=factor_ref.metadata();
    finish(result,request,data);result.laplacian=std::move(trace);
    return result;
}
}
ChartSolveResult solve_arap(const ChartSolveRequest& request) {
    return translate_allocation("solver.arap",[&] {return arap_impl(request);});
}
}
