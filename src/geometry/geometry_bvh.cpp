// SPDX-License-Identifier: GPL-3.0-or-later
// Bounded immutable triangle BVH over the fixed triangulation of one geometry revision
// (T36, Section 17). Built once per geometry revision and enhancement pass and shared
// read-only (ProposalContext); queries never touch host objects, so workers may call them.
// Used by the enhanced_v1 detectors only.
#include "cleanseams/geometry_bvh.hpp"
#include "cleanseams/errors.hpp"
#include <algorithm>
#include <atomic>
#include <cmath>
#include <limits>
#include <string>
#include <numeric>
namespace cleanseams {
namespace {
std::atomic<std::uint64_t> builds{0};
Vec3 sub(const Vec3& a,const Vec3& b) {return {a[0]-b[0],a[1]-b[1],a[2]-b[2]};}
Vec3 cross(const Vec3& a,const Vec3& b) {return {a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]};}
double dot(const Vec3& a,const Vec3& b) {return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];}
// Slab test of a ray against an AABB: the entry parameter in [0, t_max], or +inf when the
// ray misses the box within (0, t_max].
double enter(const std::array<Vec3,2>& box,const Vec3& o,const Vec3& inv,double t_max) {
    double t0=0,t1=t_max;
    for(std::size_t k=0;k<3;++k) {
        double a=(box[0][k]-o[k])*inv[k],b=(box[1][k]-o[k])*inv[k];
        // A zero direction component with the origin on a slab plane gives 0 * inf = NaN:
        // the slab then does not restrict the ray.
        if(std::isnan(a))a=-INFINITY;
        if(std::isnan(b))b=INFINITY;
        if(a>b)std::swap(a,b);
        t0=std::max(t0,a);t1=std::min(t1,b);
        if(t0>t1)return INFINITY;
    }
    return t0;
}
std::size_t min_leaf(std::size_t leaf_size) {return std::max<std::size_t>(1,(leaf_size+1)/2);}
// A node splits only above leaf_size triangles into halves, so every leaf holds at least
// (leaf_size + 1) / 2 of them: at most 2 * ceil(n / that) - 1 nodes.
std::uint64_t max_nodes(std::uint64_t triangles,std::size_t leaf_size) {
    const auto m=min_leaf(std::max<std::size_t>(1,leaf_size));
    const auto leaves=std::max<std::uint64_t>(1,(triangles+m-1)/m);
    return 2*leaves-1;
}
}
std::uint64_t triangle_bvh_bytes(std::uint64_t triangles,std::size_t leaf_size) {
    return max_nodes(triangles,leaf_size)*static_cast<std::uint64_t>(sizeof(Vec3)*2+4*sizeof(MeshIndex))+triangles*sizeof(MeshIndex);
}
std::uint64_t triangle_bvh_builds() {return builds.load();}
TriangleBvh::TriangleBvh(std::shared_ptr<const MeshSnapshot> mesh,std::size_t leaf_size,std::uint64_t memory_budget)
    :mesh_(std::move(mesh)),leaf_size_(std::max<std::size_t>(1,leaf_size)) {
    if(!mesh_)throw ContractError(ErrorCode::invalid_input,"bvh","missing mesh");
    const auto n=mesh_->triangulation().corners.size();
    if(n>std::numeric_limits<MeshIndex>::max()/2)throw ContractError(ErrorCode::memory_limit,"bvh","too many triangles for a bounded BVH");
    const auto bytes=triangle_bvh_bytes(n,leaf_size_);
    if(bytes>memory_budget)throw ContractError(ErrorCode::memory_limit,"bvh","the BVH needs "+std::to_string(bytes)+" bytes");
    order_.resize(n);std::iota(order_.begin(),order_.end(),MeshIndex{0});
    std::vector<Vec3> centres(n);
    for(std::size_t t=0;t<n;++t) {
        const auto c=corners(static_cast<MeshIndex>(t));
        for(std::size_t k=0;k<3;++k)centres[t][k]=(c[0][k]+c[1][k]+c[2][k])/3;
    }
    static_assert(sizeof(Node)==sizeof(Vec3)*2+4*sizeof(MeshIndex));
    nodes_.reserve(static_cast<std::size_t>(max_nodes(n,leaf_size_)));
    if(n)build(0,n,1,centres);
    ++builds;
}
std::array<Vec3,3> TriangleBvh::corners(MeshIndex triangle) const {
    const auto& in=mesh_->input();const auto& tri=mesh_->triangulation().corners[triangle];
    return {in.positions[static_cast<std::size_t>(in.corner_vertices[tri[0]])],in.positions[static_cast<std::size_t>(in.corner_vertices[tri[1]])],
            in.positions[static_cast<std::size_t>(in.corner_vertices[tri[2]])]};
}
bool TriangleBvh::serves(const MeshSnapshot& m) const {
    return &m==mesh_.get() || (m.geometry_key()==mesh_->geometry_key() && m.same_geometry(*mesh_));
}
std::size_t TriangleBvh::build(std::size_t begin,std::size_t end,std::size_t level,const std::vector<Vec3>& centres) {
    depth_=std::max(depth_,level);
    const auto index=nodes_.size();nodes_.push_back({});
    Node node;node.box={Vec3{INFINITY,INFINITY,INFINITY},Vec3{-INFINITY,-INFINITY,-INFINITY}};
    std::array<Vec3,2> spread={Vec3{INFINITY,INFINITY,INFINITY},Vec3{-INFINITY,-INFINITY,-INFINITY}};
    for(std::size_t i=begin;i<end;++i) {
        const auto t=order_[i];
        for(const auto& p:corners(t))for(std::size_t k=0;k<3;++k) {node.box[0][k]=std::min(node.box[0][k],p[k]);node.box[1][k]=std::max(node.box[1][k],p[k]);}
        for(std::size_t k=0;k<3;++k) {spread[0][k]=std::min(spread[0][k],centres[t][k]);spread[1][k]=std::max(spread[1][k],centres[t][k]);}
    }
    if(end-begin<=leaf_size_) {
        node.first=static_cast<MeshIndex>(begin);node.count=static_cast<MeshIndex>(end-begin);nodes_[index]=node;return index;
    }
    // Median split on the longest centroid axis; ties by triangle id (deterministic).
    std::size_t axis=0;for(std::size_t k=1;k<3;++k)if(spread[1][k]-spread[0][k]>spread[1][axis]-spread[0][axis])axis=k;
    const auto mid=begin+(end-begin)/2;
    std::nth_element(order_.begin()+static_cast<std::ptrdiff_t>(begin),order_.begin()+static_cast<std::ptrdiff_t>(mid),order_.begin()+static_cast<std::ptrdiff_t>(end),
                     [&](MeshIndex a,MeshIndex b) {const double x=centres[a][axis],y=centres[b][axis];return x!=y?x<y:a<b;});
    const auto left=build(begin,mid,level+1,centres);const auto right=build(mid,end,level+1,centres);
    node.left=static_cast<MeshIndex>(left);node.right=static_cast<MeshIndex>(right);node.count=0;
    nodes_[index]=node;
    return index;
}
RayHit TriangleBvh::ray(const Vec3& origin,const Vec3& direction,double t_max,const std::function<bool(MeshIndex)>& accept) const {
    RayHit best;best.t=t_max;
    if(nodes_.empty() || !(t_max>0))return best;
    const Vec3 inv{1/direction[0],1/direction[1],1/direction[2]};
    const auto& faces=mesh_->triangulation().faces;
    // The depth is at most 2 + log2(2^32): a fixed stack suffices (each pop pushes two).
    std::array<MeshIndex,128> stack{};std::size_t top=0;stack[top++]=0;
    while(top) {
        const auto& node=nodes_[stack[--top]];
        if(enter(node.box,origin,inv,best.t)==INFINITY)continue;
        if(node.count) {
            for(MeshIndex i=node.first;i<node.first+node.count;++i) {
                const auto t=order_[i];const auto c=corners(t);
                // Moller-Trumbore, double precision; hits at t <= 0 are behind the origin.
                const auto e1=sub(c[1],c[0]),e2=sub(c[2],c[0]);const auto p=cross(direction,e2);const double det=dot(e1,p);
                if(det==0)continue;
                const double inv_det=1/det;const auto s=sub(origin,c[0]);const double u=dot(s,p)*inv_det;
                if(u<0 || u>1)continue;
                const auto q=cross(s,e1);const double v=dot(direction,q)*inv_det;
                if(v<0 || u+v>1)continue;
                const double hit=dot(e2,q)*inv_det;
                // Closest first, then the lower triangle id; t_max itself is excluded.
                const bool closer=hit<best.t || (hit==best.t && best.hit && t<best.triangle);
                if(!(hit>0) || !closer)continue;
                if(accept && !accept(faces[t]))continue;
                best.hit=true;best.t=hit;best.triangle=t;best.face=faces[t];
            }
        } else {stack[top++]=node.right;stack[top++]=node.left;}
    }
    return best;
}
std::uint64_t TriangleBvh::memory_bytes() const {
    return nodes_.capacity()*sizeof(Node)+order_.capacity()*sizeof(MeshIndex);
}
}
