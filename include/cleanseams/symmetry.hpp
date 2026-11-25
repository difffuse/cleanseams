// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
#include "cleanseams/mesh.hpp"
#include <optional>
namespace cleanseams {
// core/symmetry.py: mirror planes of a mesh, vertex/edge/face correspondences and
// the mirror share of a seam set. Diagnostics only: generation does not force
// symmetric seams (atomic mirrored operations are Phase B, T38).
struct NearestPoints {
    std::vector<TransferIndex> index; // nearest P within the radius, or -1
    std::vector<double> distance;     // that distance, or +inf
};
// Source nearest_points: uniform grid with cell = radius, 27-cell search. Among equal
// distances the source keeps the last candidate in its enumeration order (neighbour
// cell offsets dx, dy, dz in -1..1, then ascending point index); so does this.
CLEANSEAMS_CORE_EXPORT NearestPoints nearest_points(std::span<const Vec3> points,std::span<const Vec3> queries,double radius);
// Source closest_point_on_triangles (Ericson 5.1.5) with NumPy's einsum dot products.
CLEANSEAMS_CORE_EXPORT std::vector<Vec3> closest_point_on_triangles(std::span<const Vec3> points,std::span<const Vec3> a,
    std::span<const Vec3> b,std::span<const Vec3> c);
struct SurfaceDistance {
    std::vector<double> distance;        // +inf when no vertex is within reach
    std::vector<TransferIndex> vertex;   // nearest used vertex, or -1
    std::vector<Vec3> foot;              // closest surface point, or the query itself
};
// Source surface_distance: nearest vertex within `radius` (two median edge lengths by
// default), then the exact distance to the fixed triangles around it.
CLEANSEAMS_CORE_EXPORT SurfaceDistance surface_distance(const MeshSnapshot&,std::span<const Vec3> queries,std::optional<double> radius=std::nullopt);
// Mirror plane {x : normal . x = offset}; score/far as in the source Plane.
struct CLEANSEAMS_CORE_EXPORT MirrorPlane {
    Vec3 normal{0,0,1};
    double offset=0,score=0;
    std::string source;
    double far=0;
    double signed_distance(const Vec3&) const;
    Vec3 reflect(const Vec3&) const;
    nlohmann::json to_json() const; // source Plane.as_dict
};
CLEANSEAMS_CORE_EXPORT double bbox_diagonal(const MeshSnapshot&);
// Principal axes of the area-weighted face-centre covariance through its centroid
// (pca0..2, ascending eigenvalues) and the world axes through the bounding-box centre.
CLEANSEAMS_CORE_EXPORT std::vector<MirrorPlane> candidate_planes(const MeshSnapshot&);
struct PlaneScore {
    double score=0,far=0;             // matched and far vertex shares
    std::vector<std::uint8_t> matched;
    std::vector<Vec3> foot;
};
CLEANSEAMS_CORE_EXPORT PlaneScore score_plane(const MeshSnapshot&,const MirrorPlane&,double tol=0.003,std::optional<double> far_tol=std::nullopt);
CLEANSEAMS_CORE_EXPORT MirrorPlane refine_plane(const MeshSnapshot&,const MirrorPlane&,std::span<const std::uint8_t> matched,std::span<const Vec3> foot);
CLEANSEAMS_CORE_EXPORT std::vector<MirrorPlane> detect_planes(const MeshSnapshot&,double tol=0.003,double accept=0.85,double max_far=0.01,bool refine=true);
CLEANSEAMS_CORE_EXPORT std::optional<MirrorPlane> detect_plane(const MeshSnapshot&,double tol=0.003,double accept=0.85,double max_far=0.01,bool refine=true);
struct MirrorMaps {std::vector<TransferIndex> vertex,edge,face;}; // -1 where no counterpart
CLEANSEAMS_CORE_EXPORT MirrorMaps mirror_maps(const MeshSnapshot&,const MirrorPlane&,double tol=0.003);
// Every element has a counterpart and the vertex map is an involution.
CLEANSEAMS_CORE_EXPORT bool is_combinatorial(const MirrorMaps&);
struct SeamSymmetry {
    double score=std::numeric_limits<double>::quiet_NaN(); // NaN without seams or plane
    std::vector<std::uint8_t> mask;                         // seam edges with a mirrored counterpart
};
CLEANSEAMS_CORE_EXPORT SeamSymmetry seam_symmetry(const MeshSnapshot&,std::span<const std::uint8_t> seams,const std::optional<MirrorPlane>& plane,
    double tol_edges=0.75,std::int64_t samples=5);
namespace detail {
// numpy.linalg.eigh (LAPACK dsyevd, lower triangle) for a 3x3 symmetric matrix,
// following reference LAPACK's n=3 path: dsytd2 Householder reduction, dsteqr
// implicit QL/QR with its splitting tests and selection sort, dorm2r back
// transformation. Eigenvalues ascend; column i of `vectors` is `vectors[.][i]`.
// The sign and the basis of a degenerate eigenspace follow that path; the BLAS
// kernels of the reference build are not bit-reproduced (last-ulp differences).
struct SymmetricEigen3 {std::array<double,3> values{};std::array<std::array<double,3>,3> vectors{};};
CLEANSEAMS_CORE_EXPORT SymmetricEigen3 lapack_syevd3(const std::array<std::array<double,3>,3>& lower);
// The area-weighted face-centre covariance and centroid that candidate_planes
// decomposes (NUM-03 judges eigenvalue degeneracy on exactly this matrix).
struct FaceCovariance {std::array<std::array<double,3>,3> matrix{};Vec3 centroid{0,0,0};};
CLEANSEAMS_CORE_EXPORT FaceCovariance face_centre_covariance(const MeshSnapshot&);
}
}
