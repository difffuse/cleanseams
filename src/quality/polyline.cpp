// SPDX-License-Identifier: GPL-3.0-or-later
#include "detail.hpp"
#include "cleanseams/detail/numpy_compat.hpp"
namespace cleanseams {
using namespace quality_detail;
namespace {
void storage(bool fits) {if(!fits) throw ContractError(ErrorCode::memory_limit,"quality.polyline","polyline exceeds memory allowance");}
void validate(std::span<const Vec2> points,std::uint64_t budget) {
    storage(points.size()<=budget/64 && points.size()<static_cast<std::size_t>(std::numeric_limits<std::ptrdiff_t>::max()/4));
    for(auto p:points) require(finite(p),"quality.polyline","points must be finite");
}
double length(Vec2 a,Vec2 b) {
    const double x=b[0]-a[0],y=b[1]-a[1],n=std::sqrt(x*x+y*y);
    require(std::isfinite(n),"quality.polyline","segment length overflow");return n;
}
std::size_t count(double n,std::uint64_t stride,std::uint64_t budget) {
    storage(std::isfinite(n) && n>=0 && n<static_cast<double>(std::numeric_limits<std::ptrdiff_t>::max()/4));
    const auto k=static_cast<std::size_t>(n);storage(k<=budget/stride);return k;
}
}
ResampledPolyline resample_polyline(std::span<const Vec2> points,double step,std::uint64_t budget) {
    validate(points,budget);require(std::isfinite(step),"quality.polyline.step","step must be finite");
    storage(std::max<std::size_t>(points.size(),1)<=budget/sizeof(double));
    std::vector<double> arc(std::max<std::size_t>(points.size(),1),0);
    for(std::size_t i=1;i<points.size();++i) {
        arc[i]=arc[i-1]+length(points[i-1],points[i]);require(std::isfinite(arc[i]),"quality.polyline","arc length overflow");
    }
    const auto total=arc.back();
    if(total<=0 || step<=0) return {{points.begin(),points.end()},std::move(arc)};
    const auto n=count(std::max(2.,std::ceil(total/step)+1),64,budget-points.size()*16);
    ResampledPolyline result;result.points.reserve(n);result.arc.reserve(n);
    const double spacing=total/static_cast<double>(n-1);
    for(std::size_t i=0;i<n;++i) {
        const double t=i==n-1?total:static_cast<double>(i)*spacing;result.arc.push_back(t);
        if(t>=total) {result.points.push_back(points.back());continue;}
        const auto upper=std::upper_bound(arc.begin(),arc.end(),t);const auto hi=static_cast<std::size_t>(upper-arc.begin()),lo=hi-1;
        if(arc[lo]==t) {result.points.push_back(points[lo]);continue;}
        Vec2 p{};for(std::size_t d=0;d<2;++d) p[d]=(points[hi][d]-points[lo][d])/(arc[hi]-arc[lo])*(t-arc[lo])+points[lo][d];
        result.points.push_back(p);
    }
    return result;
}
std::vector<Vec2> gaussian_smooth_open(std::span<const Vec2> points,double sigma,std::uint64_t budget) {
    validate(points,budget);require(std::isfinite(sigma),"quality.polyline.sigma","sigma must be finite");
    const auto n=points.size();if(n<3 || sigma<=0) return {points.begin(),points.end()};
    const auto k=count(std::ceil(3*sigma),80,budget-n*64),kk=std::min(k,n-1);
    std::vector<double> weights(2*k+1);
    for(std::size_t j=0;j<weights.size();++j) {const double x=(static_cast<double>(j)-static_cast<double>(k))/sigma;weights[j]=std::exp(-.5*x*x);}
    const double sum=detail::numpy_sum(weights); // source g /= g.sum(): NumPy pairwise sum
    for(auto& w:weights) w/=sum;
    std::vector<Vec2> padded(n+2*k);
    for(std::size_t j=0;j<padded.size();++j) {
        auto index=static_cast<std::ptrdiff_t>(j)-static_cast<std::ptrdiff_t>(k);
        index=std::clamp(index,-static_cast<std::ptrdiff_t>(kk),static_cast<std::ptrdiff_t>(n+kk-1));
        if(index<0) {
            for(std::size_t d=0;d<2;++d) padded[j][d]=2*points.front()[d]-points[static_cast<std::size_t>(-index)][d];
        } else if(static_cast<std::size_t>(index)>=n) {
            for(std::size_t d=0;d<2;++d) padded[j][d]=2*points.back()[d]-points[2*n-2-static_cast<std::size_t>(index)][d];
        } else padded[j]=points[static_cast<std::size_t>(index)];
        require(finite(padded[j]),"quality.polyline","reflection overflow");
    }
    // Source np.convolve: every output is one Accelerate cblas_ddot over unit strides, whose
    // accumulation order depends on the operands' pointer alignment (measured on the oracle
    // interpreter: a view 8 bytes into an array and its aligned copy give different sums in
    // 628 of 1120 random dots), so it is not a function of the values. Summed in kernel
    // order here; the results agree to the last ulp (Gate A F22, proposed NUM entry).
    std::vector<Vec2> result(n,Vec2{0,0});
    for(std::size_t i=0;i<n;++i) for(std::size_t j=0;j<weights.size();++j)
        for(std::size_t d=0;d<2;++d) result[i][d]+=padded[i+j][d]*weights[weights.size()-1-j];
    result.front()=points.front();result.back()=points.back();return result;
}
ArcJaggedness arc_jaggedness(std::span<const Vec2> points,std::optional<double> scale,bool straight,std::uint64_t budget) {
    validate(points,budget);require(!scale || std::isfinite(*scale),"quality.polyline.scale","scale must be finite");
    ArcJaggedness report;report.samples=points.size();if(points.size()<3) return report;
    std::vector<double> lengths;lengths.reserve(points.size()-1);
    for(std::size_t i=1;i<points.size();++i) {const auto n=length(points[i-1],points[i]);if(n>0) lengths.push_back(n);}
    if(lengths.empty()) return report;
    std::sort(lengths.begin(),lengths.end());const auto n=lengths.size();
    report.h=n%2?lengths[n/2]:(lengths[n/2-1]+lengths[n/2])*.5;
    double sum=0;std::size_t samples=0;
    auto accumulate=[&](double deviation){const double d=deviation/report.h;sum+=d*d;report.maximum=std::max(report.maximum,d);++samples;};
    if(straight) {
        const auto a=points.front(),b=points.back();const double chord=length(a,b);
        if(chord<=0) for(auto p:points) accumulate(length(a,p));
        else {
            const Vec2 normal{-(b[1]-a[1])/chord,(b[0]-a[0])/chord};
            for(auto p:points) accumulate(std::abs((p[0]-a[0])*normal[0]+(p[1]-a[1])*normal[1]));
        }
    } else {
        const auto res=resample_polyline(points,report.h,budget/3);
        const double sigma=(scale && *scale!=0?*scale:4*report.h)/report.h;
        const auto smooth=gaussian_smooth_open(res.points,sigma,budget/3);
        for(std::size_t i=0;i<res.points.size();++i) accumulate(length(res.points[i],smooth[i]));
    }
    report.rms=std::sqrt(sum/static_cast<double>(samples));return report;
}
ChainTurns chain_turns(const MeshSnapshot& mesh,std::span<const TransferIndex> vertices) {
    require(vertices.size()<static_cast<std::size_t>(std::numeric_limits<std::int64_t>::max()),"quality.chain","chain length overflow");
    for(auto v:vertices) require(v>=0 && static_cast<std::uint64_t>(v)<mesh.input().positions.size(),"quality.chain","vertex out of range");
    for(std::size_t i=1;i<vertices.size();++i) require(mesh.edge_index(vertices[i-1],vertices[i])>=0,"quality.chain","consecutive vertices must share a mesh edge");
    ChainTurns result;result.edges=static_cast<std::int64_t>(vertices.size())-1;
    const double threshold=std::cos(30*(std::acos(-1.)/180));
    auto turn=[&](TransferIndex a,TransferIndex v,TransferIndex b) {
        const auto vertex=static_cast<MeshIndex>(v),incoming=static_cast<MeshIndex>(mesh.edge_index(a,v)),outgoing=static_cast<MeshIndex>(mesh.edge_index(v,b));
        const auto di=mesh.tangent_dir(vertex,incoming),dout=mesh.tangent_dir(vertex,outgoing);
        const Vec3 negative{-di[0],-di[1],-di[2]};if(dot(negative,dout)<threshold) ++result.turns;
        if(mesh.tables().vert_regular[vertex] && mesh.tables().fans[vertex]->straight(incoming)!=outgoing) ++result.lane_changes;
    };
    if(vertices.size()>2 && vertices.front()==vertices.back()) turn(vertices[vertices.size()-2],vertices[0],vertices[1]);
    for(std::size_t i=1;i+1<vertices.size();++i) turn(vertices[i-1],vertices[i],vertices[i+1]);
    return result;
}
}
