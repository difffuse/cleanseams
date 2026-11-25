// SPDX-License-Identifier: GPL-3.0-or-later
// Source core/output.py: _hull and orient_islands.
#include "cleanseams/output.hpp"
#include "cleanseams/errors.hpp"
#include "cleanseams/detail/numpy_blas.hpp"
#include <algorithm>
#include <cmath>
#include <numbers>
#include <set>
namespace cleanseams {
namespace {
double cross(const Vec2& o,const Vec2& a,const Vec2& b) {return (a[0]-o[0])*(b[1]-o[1])-(a[1]-o[1])*(b[0]-o[0]);}
// Python/NumPy float remainder (result takes the divisor's sign).
double python_mod(double a,double b) {
    double mod=std::fmod(a,b);
    if(mod!=0) {if((b<0)!=(mod<0))mod+=b;}
    else mod=std::copysign(0.0,b);
    return mod;
}
std::array<Vec2,2> rotation(double angle) {
    const double c=std::cos(angle),s=std::sin(angle);return {Vec2{c,-s},Vec2{s,c}};
}
}
nlohmann::json IslandOrientation::info() const {
    nlohmann::json angles=nlohmann::json::array();for(const auto& [ch,a]:angles_radians)angles.push_back(nlohmann::json::array({ch,a}));
    return {{"angles_radians",angles},{"protected_charts",protected_charts}};
}
std::vector<Vec2> output_hull(std::span<const Vec2> points) {
    // sorted(set(map(tuple, points))): the first of equal points survives, then a
    // lexicographic sort (distinct points never tie).
    std::vector<Vec2> pts;
    for(const auto& p:points)if(std::none_of(pts.begin(),pts.end(),[&](const Vec2& q) {return q[0]==p[0] && q[1]==p[1];}))pts.push_back(p);
    std::sort(pts.begin(),pts.end(),[](const Vec2& a,const Vec2& b) {return a[0]!=b[0]?a[0]<b[0]:a[1]<b[1];});
    if(pts.size()<=2)return pts;
    std::vector<Vec2> lower,upper;
    for(const auto& p:pts) {
        while(lower.size()>1 && cross(lower[lower.size()-2],lower.back(),p)<=0)lower.pop_back();
        lower.push_back(p);
    }
    for(auto it=pts.rbegin();it!=pts.rend();++it) {
        while(upper.size()>1 && cross(upper[upper.size()-2],upper.back(),*it)<=0)upper.pop_back();
        upper.push_back(*it);
    }
    lower.pop_back();upper.pop_back();lower.insert(lower.end(),upper.begin(),upper.end());
    return lower;
}
IslandOrientation orient_islands(const CutSurface& cs,std::span<const Vec2> uv,std::span<const MeshIndex> pinned_corners) {
    const auto& mesh=cs.mesh();const auto& t=cs.tables();const auto nc=mesh.input().corner_vertices.size();
    if(uv.size()!=nc)throw ContractError(ErrorCode::invalid_input,"orientation.uv","one UV per original corner required");
    IslandOrientation out;out.uv.assign(uv.begin(),uv.end());
    std::set<MeshIndex> fixed;
    for(auto c:pinned_corners) {
        if(c>=nc)throw ContractError(ErrorCode::invalid_input,"orientation.pins","pinned corner out of range");
        fixed.insert(t.chart_of_corner[c]);
    }
    for(MeshIndex ch=0;ch<t.n_charts;++ch) {
        if(fixed.contains(ch))continue;
        std::vector<Vec2> boundary;bool finite=true;
        for(auto loop:t.chart_loops[ch])for(auto c:t.loops[loop]) {
            boundary.push_back(out.uv[c]);finite=finite && std::isfinite(out.uv[c][0]) && std::isfinite(out.uv[c][1]);
        }
        if(boundary.size()<3 || !finite)continue;
        const auto hull=output_hull(boundary);
        // Candidate angles: every hull edge direction modulo pi/2, sorted and distinct.
        std::vector<double> candidates;
        for(std::size_t i=0;i<hull.size();++i) {
            const auto& a=hull[i];const auto& b=hull[(i+1)%hull.size()];
            candidates.push_back(python_mod(std::atan2(b[1]-a[1],b[0]-a[0]),std::numbers::pi/2));
        }
        std::sort(candidates.begin(),candidates.end());candidates.erase(std::unique(candidates.begin(),candidates.end()),candidates.end());
        bool have=false;double best_area=0,best_angle=0;Vec2 best_extent{};
        for(double angle:candidates) {
            const auto rotated=detail::accelerate_rows_times_2x2(hull,rotation(angle));
            Vec2 lo=rotated[0],hi=rotated[0];
            for(const auto& p:rotated)for(std::size_t k=0;k<2;++k) {lo[k]=std::min(lo[k],p[k]);hi[k]=std::max(hi[k],p[k]);}
            const Vec2 extent{hi[0]-lo[0],hi[1]-lo[1]};const double area=extent[0]*extent[1];
            if(!have || area<best_area || (area==best_area && angle<best_angle)) {have=true;best_area=area;best_angle=angle;best_extent=extent;}
        }
        double angle=best_angle;
        if(best_extent[1]>best_extent[0])angle+=std::numbers::pi/2;
        std::vector<MeshIndex> ids;for(auto f:t.chart_faces[ch])for(auto c:mesh.face_corners(f))ids.push_back(c);
        Vec2 centre{0,0};for(auto c:ids)for(std::size_t k=0;k<2;++k)centre[k]+=out.uv[c][k];
        for(std::size_t k=0;k<2;++k)centre[k]/=static_cast<double>(ids.size());
        std::vector<Vec2> shifted;shifted.reserve(ids.size());
        for(auto c:ids)shifted.push_back({out.uv[c][0]-centre[0],out.uv[c][1]-centre[1]});
        const auto turned=detail::accelerate_rows_times_2x2(shifted,rotation(angle));
        for(std::size_t i=0;i<ids.size();++i)out.uv[ids[i]]={turned[i][0]+centre[0],turned[i][1]+centre[1]};
        out.angles_radians[ch]=angle;
    }
    out.protected_charts.assign(fixed.begin(),fixed.end());
    return out;
}
}
