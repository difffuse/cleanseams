// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
#include "cleanseams/mesh.hpp"
#include <cmath>
namespace cleanseams::detail {
inline Vec3 add(const Vec3& a,const Vec3& b) {return {a[0]+b[0],a[1]+b[1],a[2]+b[2]};}
inline Vec3 sub(const Vec3& a,const Vec3& b) {return {a[0]-b[0],a[1]-b[1],a[2]-b[2]};}
inline Vec3 mul(const Vec3& a,double s) {return {a[0]*s,a[1]*s,a[2]*s};}
inline Vec3 divide(const Vec3& a,double s) {return {a[0]/s,a[1]/s,a[2]/s};}
inline double dot(const Vec3& a,const Vec3& b) {return (a[0]*b[0]+a[1]*b[1])+a[2]*b[2];}
inline Vec3 cross(const Vec3& a,const Vec3& b) {return {a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]};}
inline double norm(const Vec3& a) {return std::sqrt(dot(a,a));}
inline double orient(const Vec2& a,const Vec2& b,const Vec2& c) {return (b[0]-a[0])*(c[1]-a[1])-(c[0]-a[0])*(b[1]-a[1]);}
std::vector<Vec2> project_polygon(std::span<const Vec3>,const Vec3&,bool quad_projection);
}
