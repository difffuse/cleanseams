// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
#include "cleanseams/mesh.hpp"
#include <functional>
#include <memory>
namespace cleanseams {
// Bounded immutable bounding-volume hierarchy over the fixed triangulation of one geometry
// revision (T36/T37, Section 17). It keeps its mesh snapshot (immutable) alive and reads the
// triangles from it, so it never reads host objects and any worker may query it. Median
// splits on the longest centroid axis (ties by triangle id), leaves of at most `leaf_size`
// triangles: fewer than 2 * triangles nodes and a depth of at most 2 + log2(triangles).
// Enhanced detectors only; the migrated exposure field keeps its own z-buffer.
struct RayHit {bool hit=false;double t=0;MeshIndex triangle=0,face=0;};
class CLEANSEAMS_CORE_EXPORT TriangleBvh {
public:
    // Refuses (memory_limit) when triangle_bvh_bytes exceeds `memory_budget`.
    explicit TriangleBvh(std::shared_ptr<const MeshSnapshot>,std::size_t leaf_size=4,std::uint64_t memory_budget=active_memory_budget());
    const MeshSnapshot& mesh() const {return *mesh_;}
    const std::string& geometry_key() const {return mesh_->geometry_key();}
    // The geometry revision this BVH answers for: the same snapshot, or an exactly equal
    // geometry (the content key, confirmed by full equality; a hash never authorizes reuse).
    bool serves(const MeshSnapshot&) const;
    // Closest hit with 0 < t < t_max of the ray origin + t * direction (direction need not be
    // unit; t is in its units) on a triangle whose face `accept` admits (all when empty).
    // Double-precision Moller-Trumbore; equal distances resolve to the lower triangle id,
    // so the answer does not depend on the traversal order.
    RayHit ray(const Vec3& origin,const Vec3& direction,double t_max,const std::function<bool(MeshIndex face)>& accept={}) const;
    std::size_t nodes() const {return nodes_.size();}
    std::size_t triangles() const {return order_.size();}
    std::size_t depth() const {return depth_;}
    std::uint64_t memory_bytes() const;
private:
    struct Node {std::array<Vec3,2> box{};MeshIndex left=0,right=0,first=0,count=0;};
    std::shared_ptr<const MeshSnapshot> mesh_;
    std::size_t leaf_size_=4,depth_=0;
    std::vector<MeshIndex> order_;   // triangle ids, leaves own contiguous ranges
    std::vector<Node> nodes_;
    std::array<Vec3,3> corners(MeshIndex triangle) const;
    std::size_t build(std::size_t begin,std::size_t end,std::size_t level,const std::vector<Vec3>& centres);
};
// Upper bound of the bytes of a BVH over `triangles` triangles (nodes and order).
CLEANSEAMS_CORE_EXPORT std::uint64_t triangle_bvh_bytes(std::uint64_t triangles,std::size_t leaf_size=4);
// BVH builds since the process started (reuse evidence: one per geometry revision and pass).
CLEANSEAMS_CORE_EXPORT std::uint64_t triangle_bvh_builds();
}
