// SPDX-License-Identifier: GPL-3.0-or-later
#include "cleanseams/detail/mesh_math.hpp"
#include "cleanseams/detail/numpy_compat.hpp"
#include <algorithm>
#include <numeric>
namespace cleanseams {
using namespace detail;
std::vector<Vec2> detail::project_polygon(std::span<const Vec3> p,const Vec3& n,bool quad) {
    // The source's vectorized quad projection (_newell_2d) uses einsum dot products;
    // its general polygon projection uses plain dots.
    auto project=[quad](const Vec3& a,const Vec3& b) {return quad?numpy_einsum_dot3(a,b):dot(a,b);};
    Vec3 center{};for(const auto& point:p) center=add(center,point);center=divide(center,static_cast<double>(p.size()));
    auto e1=sub(p[1],p[0]);e1=sub(e1,mul(n,project(e1,n)));
    if(norm(e1)<1e-300) {
        e1={1,0,0};e1=sub(e1,mul(n,project(e1,n)));
        // The source's vectorized quad projection has only the X fallback.
        // Its general polygon projection also tries Y. Keep that distinction.
        if(!quad && norm(e1)<1e-300) {e1={0,1,0};e1=sub(e1,mul(n,dot(e1,n)));}
    }
    e1=divide(e1,std::max(norm(e1),1e-300));const auto e2=cross(n,e1);
    std::vector<Vec2> result;result.reserve(p.size());
    for(const auto& point:p) {const auto r=sub(point,center);result.push_back({project(r,e1),project(r,e2)});}
    return result;
}
bool polygon_self_intersects(std::span<const Vec2> p,double eps) {
    if(p.size()<4) return false;
    for(std::size_t i=0;i<p.size();++i) for(std::size_t j=i+2;j<p.size();++j) {
        if(i==0 && j+1==p.size()) continue;
        const auto& a=p[i];const auto& b=p[(i+1)%p.size()];const auto& c=p[j];const auto& d=p[(j+1)%p.size()];
        if(orient(c,d,a)*orient(c,d,b)<-eps && orient(a,b,c)*orient(a,b,d)<-eps) return true;
    }
    return false;
}
namespace {
using Triangles=std::vector<std::array<MeshIndex,3>>;
std::optional<Triangles> ear_clip(std::span<const Vec2> p) {
    const auto n=p.size();if(n<3) return std::nullopt;if(n==3) return Triangles{{0,1,2}};
    if(polygon_self_intersects(p)) return std::nullopt;
    auto lo=p[0],hi=p[0];double forward=0.,backward=0.;
    for(std::size_t i=0;i<n;++i) {
        for(std::size_t k=0;k<2;++k) {lo[k]=std::min(lo[k],p[i][k]);hi[k]=std::max(hi[k],p[i][k]);}
        forward+=p[i][0]*p[(i+1)%n][1];backward+=p[(i+1)%n][0]*p[i][1];
    }
    const auto extent=std::max({hi[0]-lo[0],hi[1]-lo[1],1e-300});
    const auto scale=extent*extent,eps=1e-12*scale,area2=forward-backward;
    std::vector<MeshIndex> indices(n);std::iota(indices.begin(),indices.end(),MeshIndex{0});
    if(area2<0) std::reverse(indices.begin(),indices.end());Triangles triangles;triangles.reserve(n-2);bool strict=true;
    while(indices.size()>3) {
        const auto m=indices.size();bool found=false;
        for(std::size_t t=0;t<m;++t) {
            const auto i=indices[(t+m-1)%m],j=indices[t],k=indices[(t+1)%m];
            const auto& a=p[i];const auto& b=p[j];const auto& c=p[k];const auto o=orient(a,b,c);
            if(strict?o<=eps:o<-eps) continue;
            bool blocked=false;
            for(auto v:indices) if(v!=i && v!=j && v!=k) {
                const auto& point=p[v];
                if(o>eps) {
                    const bool inside=orient(a,b,point)>-eps && orient(b,c,point)>-eps && orient(c,a,point)>-eps;
                    const auto close=1e-9*std::sqrt(scale);
                    const bool coincident=std::max(std::abs(point[0]-a[0]),std::abs(point[1]-a[1]))<=close ||
                        std::max(std::abs(point[0]-c[0]),std::abs(point[1]-c[1]))<=close;
                    if(inside && !coincident) {blocked=true;break;}
                } else {
                    bool on=std::abs(orient(a,c,point))<=eps;
                    for(std::size_t axis=0;axis<2;++axis) on=on && std::min(a[axis],c[axis])-1e-9<=point[axis] && point[axis]<=std::max(a[axis],c[axis])+1e-9;
                    if(on) {blocked=true;break;}
                }
            }
            if(blocked) continue;
            triangles.push_back({i,j,k});indices.erase(indices.begin()+static_cast<std::ptrdiff_t>(t));found=true;break;
        }
        if(!found) {if(strict) {strict=false;continue;}return std::nullopt;}
        strict=true;
    }
    triangles.push_back({indices[0],indices[1],indices[2]});
    if(area2<0) for(auto& tri:triangles) std::swap(tri[0],tri[2]);
    return triangles;
}
}
void MeshSnapshot::build_triangulation() {
    auto& result=triangulation_;const auto nf=geometry_.face_area.size();
    const auto count=input_.corner_vertices.size()-2*nf;
    result.corners.reserve(count);result.faces.reserve(count);result.original_edges.reserve(count);
    for(std::size_t f=0;f<nf;++f) {
        const auto start=input_.face_offsets[f],size=input_.face_offsets[f+1]-start;
        Triangles triangles;
        if(size==3) triangles={{0,1,2}};
        else if(size==4) {
            const auto p=face_points(static_cast<MeshIndex>(f));
            const auto a0=norm(cross(sub(p[1],p[0]),sub(p[2],p[0]))),a1=norm(cross(sub(p[2],p[0]),sub(p[3],p[0])));
            const auto b0=norm(cross(sub(p[2],p[1]),sub(p[3],p[1]))),b1=norm(cross(sub(p[3],p[1]),sub(p[0],p[1])));
            bool use_b=std::min(b0,b1)>std::min(a0,a1);
            const auto q=project_polygon(p,geometry_.face_normal[f],true);
            const bool ok_a=orient(q[0],q[1],q[2])>0 && orient(q[0],q[2],q[3])>0;
            const bool ok_b=orient(q[1],q[2],q[3])>0 && orient(q[1],q[3],q[0])>0;
            if(ok_a!=ok_b) use_b=ok_b;
            triangles=use_b?Triangles{{1,2,3},{1,3,0}}:Triangles{{0,1,2},{0,2,3}};
        } else {
            std::optional<Triangles> ears;
            if(geometry_.face_area[f]>0 && geometry_.face_degenerate_kind[f]!=2)
                ears=ear_clip(project_polygon(face_points(static_cast<MeshIndex>(f)),geometry_.face_normal[f],false));
            if(ears) triangles=std::move(*ears);
            else {++result.fallbacks;for(MeshIndex j=1;static_cast<std::uint64_t>(j)+1<size;++j) triangles.push_back({0,j,j+1});}
        }
        for(auto tri:triangles) {
            for(auto& c:tri) c=static_cast<MeshIndex>(start+c);
            result.corners.push_back(tri);result.faces.push_back(static_cast<MeshIndex>(f));
            std::array<TransferIndex,3> edges{-1,-1,-1};
            for(std::size_t side=0;side<3;++side) {
                const auto a=tri[side],b=tri[(side+1)%3];
                if(corner_next(a)==b) edges[side]=tables_.corner_edge[a];
                else if(corner_next(b)==a) edges[side]=tables_.corner_edge[b];
            }
            result.original_edges.push_back(edges);
        }
    }
}
}
