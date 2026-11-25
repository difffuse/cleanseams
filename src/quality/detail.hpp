// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
#include "cleanseams/quality.hpp"
#include "../topology/detail.hpp"
#include <cmath>
namespace cleanseams::quality_detail {
using cut_detail::require;
inline bool finite(Vec2 p) {return std::isfinite(p[0]) && std::isfinite(p[1]);}
inline double orient(Vec2 a,Vec2 b,Vec2 c) {return (b[0]-a[0])*(c[1]-a[1])-(c[0]-a[0])*(b[1]-a[1]);}
inline double dot(Vec3 a,Vec3 b) {return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];}
inline Vec3 sub(Vec3 a,Vec3 b) {return {a[0]-b[0],a[1]-b[1],a[2]-b[2]};}
inline Vec3 div(Vec3 a,double v) {return {a[0]/v,a[1]/v,a[2]/v};}
inline Vec3 cross(Vec3 a,Vec3 b) {return {a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]};}
inline void validate_view(const MeshSnapshot& mesh,CornerUvView uv) {
    const auto nc=mesh.input().corner_vertices.size();
    if(uv.corners.empty()) require(uv.values.size()==nc,"quality.uv","dense buffer corner count mismatch");
    else {
        require(uv.values.size()==uv.corners.size(),"quality.uv","compact buffer length mismatch");
        for(std::size_t i=0;i<uv.corners.size();++i) {
            require(uv.corners[i]<nc,"quality.corners","corner out of range");
            require(i==0 || uv.corners[i-1]<uv.corners[i],"quality.corners","IDs must be strictly increasing");
        }
    }
}
inline Vec2 point(CornerUvView uv,MeshIndex c) {
    if(uv.corners.empty()) return uv.values[c];
    const auto it=std::lower_bound(uv.corners.begin(),uv.corners.end(),c);
    require(it!=uv.corners.end() && *it==c,"quality.corners","missing required chart corner");
    return uv.values[static_cast<std::size_t>(it-uv.corners.begin())];
}
inline double signed_area(const MeshSnapshot& mesh,CornerUvView uv,MeshIndex f) {
    double area=0;bool good=true;
    for(auto c:mesh.face_corners(f)) {
        const auto a=point(uv,c),b=point(uv,mesh.corner_next(c));
        area+=a[0]*b[1]-b[0]*a[1];good=good && finite(a);
    }
    return good?.5*area:quality_nan;
}
}
