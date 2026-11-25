// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
#include "cleanseams/fields.hpp"
#include "cleanseams/quality.hpp"
#include <any>
#include <functional>
namespace cleanseams {
// Source core/seam_metrics.py: the columns a comparison between seam methods needs.
// Seam masks have one flag per mesh edge; maps are dense per original corner (NaN
// where unavailable). NaN results stay NaN (JSON conversion is the caller's choice).
CLEANSEAMS_CORE_EXPORT double seam_length(const MeshSnapshot&,std::span<const std::uint8_t> seams);
CLEANSEAMS_CORE_EXPORT double seam_length_norm(const MeshSnapshot&,std::span<const std::uint8_t> seams);
CLEANSEAMS_CORE_EXPORT double seam_length_ratio(const MeshSnapshot&,std::span<const std::uint8_t> seams,std::span<const std::uint8_t> reference);
// Interior two-corner edges whose faces disagree (by more than eps) on a shared vertex's UV.
CLEANSEAMS_CORE_EXPORT std::vector<std::uint8_t> seams_from_uv(const MeshSnapshot&,std::span<const Vec2> uv,double eps=1e-6);
struct IslandLabels {std::vector<MeshIndex> labels;MeshIndex count=0;};
// Face labels across unseamed two-face edges, numbered by lowest face.
CLEANSEAMS_CORE_EXPORT IslandLabels islands_from_seams(const MeshSnapshot&,std::span<const std::uint8_t> seams);
// Vertex hops to the nearest seam vertex (breadth first over mesh edges), -1 if unreachable.
CLEANSEAMS_CORE_EXPORT std::vector<std::int64_t> seam_hop_distance(const MeshSnapshot&,std::span<const std::uint8_t> seams);
// Source _edge_hops(mesh, seams, target): per edge, the larger end distance to `target`
// (0 on target edges, the vertex count where unreachable). `seams` is unused, as in source.
CLEANSEAMS_CORE_EXPORT std::vector<double> edge_hops(const MeshSnapshot&,std::span<const std::uint8_t> seams,std::span<const std::uint8_t> target);
struct CLEANSEAMS_CORE_EXPORT SeamAgreement {
    double exact=0,one_ring=0,mean_hops=0,exact_of_ours=0,one_ring_of_ours=0,mean_hops_of_ours=0,island_count_ratio=0;
    std::uint64_t islands=0,reference_islands=0;
    nlohmann::json to_json() const;
};
CLEANSEAMS_CORE_EXPORT SeamAgreement seam_agreement(const MeshSnapshot&,std::span<const std::uint8_t> seams,std::span<const std::uint8_t> reference);
// Seam vertices of degree one that are not on the mesh boundary.
CLEANSEAMS_CORE_EXPORT std::uint64_t dangling_seam_ends(const MeshSnapshot&,std::span<const std::uint8_t> seams);
// Andrew's monotone chain, counter-clockwise, over the finite points (fewer than three
// finite points, or fewer than three distinct ones, are returned without a hull).
CLEANSEAMS_CORE_EXPORT std::vector<Vec2> convex_hull(std::span<const Vec2> points);
struct MinAreaRect {double width=0,height=0,angle=0;};
CLEANSEAMS_CORE_EXPORT MinAreaRect min_area_rect(std::span<const Vec2> points);
struct SkylinePlacement {double x=0,y=0;bool rotated=false;};
struct SkylinePack {std::vector<SkylinePlacement> placements;double height=0;};
// Deterministic bottom-left skyline packing, tallest first (a packing-efficiency proxy,
// not an atlas packer).
CLEANSEAMS_CORE_EXPORT SkylinePack skyline_pack(std::span<const std::array<double,2>> sizes,double width,bool allow_rotate=true);
struct IslandRect {double width=0,height=0,area=0;};
CLEANSEAMS_CORE_EXPORT std::vector<IslandRect> island_rects(const CutSurface&,std::span<const Vec2> uv);
// Share of a square atlas covered after the skyline pack of minimum-area rectangles with
// `padding`; an estimate that compares seam sets, not an actual packed atlas.
CLEANSEAMS_CORE_EXPORT double packing_efficiency(const CutSurface&,std::span<const Vec2> uv,double padding=4.0/1024.0,std::int64_t iterations=20);
// Source _triangle_distortion: per triangle of the fixed triangulation.
struct TriangleDistortion {std::vector<MeshIndex> chart;std::vector<double> area3,anisotropy,area;};
CLEANSEAMS_CORE_EXPORT TriangleDistortion triangle_distortion(const std::shared_ptr<const CutSurface>&,std::span<const Vec2> uv);
struct CLEANSEAMS_CORE_EXPORT ChartWorst {
    std::int64_t chart_worst_id=-1;
    double chart_worst_aniso_p95=quality_nan,chart_worst_area_p95=quality_nan;
    std::vector<double> per_chart_aniso_p95,per_chart_area_p95;
    nlohmann::json to_json() const;
};
CLEANSEAMS_CORE_EXPORT ChartWorst chart_worst(const std::shared_ptr<const CutSurface>&,std::span<const Vec2> uv);
// Weighted mean of a per-face value over `rings` sector-spread rings (never across a cut).
CLEANSEAMS_CORE_EXPORT std::vector<double> windowed_face_values(const CutSurface&,std::span<const double> values,
    std::span<const double> weights,std::int64_t rings=2);
struct CLEANSEAMS_CORE_EXPORT Hotspot {
    double hotspot_max=quality_nan,hotspot_area_max=quality_nan;
    std::int64_t hotspot_face=-1,hotspot_chart=-1;
    double hotspot_area_pct=quality_nan,hotspot_chart_area_pct=quality_nan;
    std::vector<double> windowed_anisotropy; // empty when the mesh has no triangle
    nlohmann::json to_json() const; // without windowed_anisotropy
};
CLEANSEAMS_CORE_EXPORT Hotspot hotspot(const std::shared_ptr<const CutSurface>&,std::span<const Vec2> uv,std::int64_t rings=2,double hot=2.0);
// The audit tool's chain and island summary (source chain_statistics).
CLEANSEAMS_CORE_EXPORT nlohmann::json chain_statistics(std::span<const std::uint8_t> seams,const std::shared_ptr<const CutSurface>&);
// ------------------------------------------------------------------ seam_report
// Mirror-plane part of seam_report. The source always calls symmetry.detect_plane(mesh,
// tol) and seam_symmetry(mesh, seams, plane). `description` is Plane.as_dict(); `native`
// carries the T14 plane object when the native service detected it.
struct DetectedPlane {double score=0;nlohmann::json description;std::any native;};
struct SymmetryService {
    std::function<std::optional<DetectedPlane>(const MeshSnapshot&,double tol)> detect_plane;
    std::function<double(const MeshSnapshot&,std::span<const std::uint8_t> seams,const DetectedPlane&)> seam_symmetry;
    bool available() const {return detect_plane && seam_symmetry;}
};
// The native T14 module (detect_plane, seam_symmetry): seam_report's default service.
CLEANSEAMS_CORE_EXPORT SymmetryService native_symmetry_service();
// The caller-kept `meta` dict: exposure field and mirror plane computed once per mesh.
// `plane` is empty until detected; an engaged empty optional means "no plane".
struct SeamReportContext {
    std::shared_ptr<const ExposureField> exposure;
    std::optional<std::optional<DetectedPlane>> plane;
};
struct SeamReportOptions {std::uint64_t K=64;double hidden_below=0.3,symmetry_tol=0.003,hot=2.0;};
// All seam metrics of one seam set (source seam_report, same keys). The mirror plane is
// detected with `symmetry` when a service is injected (tests supply recorded results),
// otherwise with native_symmetry_service(), as the source always detects it.
CLEANSEAMS_CORE_EXPORT nlohmann::json seam_report(const std::shared_ptr<const MeshSnapshot>&,std::span<const std::uint8_t> seams,
    std::shared_ptr<const CutSurface> cs=nullptr,SeamReportContext* context=nullptr,
    std::optional<std::span<const std::uint8_t>> reference=std::nullopt,std::optional<std::span<const Vec2>> uv=std::nullopt,
    const SeamReportOptions& options={},const SymmetryService& symmetry={});
// ------------------------------------------------------------------ Phase B (T32)
// Sector-aware distortion field of a measured map (Section 14.2), used by enhanced_v1's
// relief proposals. Values live on cut-surface sectors, never on original vertices: the
// two sides of a slit through one original vertex are two sectors with their own mass and
// value, and smoothing follows the actual cut adjacency (a triangle joins only the sectors
// of its own corners, so it never crosses a cut or reaches another chart).
//
//   log_area       = ln(s0*s1 / chart density)   signed: < 0 compressed, > 0 expanded.
//                    AREA units: a uniform linear checker scale k gives 2 ln k, so the
//                    threshold ln(1.25) is a 25% area change (an 11.8% linear change).
//   log_anisotropy = ln(s0 / s1) >= 0            shape units (ratio of singular values).
//
// Measured triangles are the distortion report's (finite, non-sliver, non-exempt). Each
// gives a third of its 3D area to each corner's sector (mass); a sector's value is the
// mass-weighted mean of its triangles. One smoothing step averages each triangle's three
// sectors and gives the triangle averages back with the same weights, which conserves
// sum(mass * value) exactly (up to rounding) for every chart.
inline constexpr const char* sector_distortion_version="sector_distortion.v1";
struct CLEANSEAMS_CORE_EXPORT SectorDistortionField {
    std::string version=sector_distortion_version;
    std::int64_t smoothing_steps=0;
    std::vector<MeshIndex> charts;                           // measured charts of the snapshot (sorted)
    std::vector<double> mass;                                // per sector
    std::vector<double> log_area,log_anisotropy;             // per sector, smoothed; NaN without mass
    std::vector<double> raw_log_area,raw_log_anisotropy;     // per sector, before smoothing
    std::vector<double> chart_area,chart_mass;               // per chart of the surface: measured 3D area, total sector mass
    std::uint64_t triangles=0;                               // measured triangles
    double mass_error=0;                                     // max over charts of |mass - area| / area
    double moment_error=0;                                   // max over charts and fields of the relative change of sum(mass*value)
    static nlohmann::json units();
    nlohmann::json summary() const;                          // without the per-sector arrays
};
CLEANSEAMS_CORE_EXPORT SectorDistortionField sector_distortion_field(const QualitySnapshot&,std::int64_t smoothing_steps=2);
// Phase B (T39, Section 20): normalized seam geometry for comparisons across rigid motions,
// uniform scales and tessellations. Internal seams (effective cuts minus the physical boundary)
// are sampled uniformly at most `spacing` apart and expressed as (p - centroid) / sqrt(area)
// of the measured (non-exempt, positive-area) faces, next to the chart topology. Different
// meshes are compared by these point sets, never by edge IDs.
inline constexpr const char* seam_geometry_version="seam_geometry.v1";
struct CLEANSEAMS_CORE_EXPORT SeamGeometry {
    std::string version=seam_geometry_version;
    double area=0,seam_length=0,normalized_length=0,spacing=0;
    Vec3 centroid{0,0,0};
    std::uint64_t seam_edges=0,charts=0,disks=0,boundary_loops=0;
    std::vector<Vec3> samples;                               // normalized units
    nlohmann::json summary() const;                          // without the samples
};
CLEANSEAMS_CORE_EXPORT SeamGeometry seam_geometry(const CutSurface&,std::span<const std::uint8_t> exempt={},double spacing=1.0/64.0);
// Symmetric Hausdorff distance of two sample sets in normalized units (0 when both are empty,
// +inf when only one is).
CLEANSEAMS_CORE_EXPORT double seam_hausdorff(const SeamGeometry&,const SeamGeometry&);
}
