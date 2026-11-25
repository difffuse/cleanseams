// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
#include "cleanseams/cut_surface.hpp"
#include <limits>
#include <functional>
namespace cleanseams {
inline constexpr double quality_nan=std::numeric_limits<double>::quiet_NaN();
class CLEANSEAMS_CORE_EXPORT WeightedDistribution {
public:
    WeightedDistribution(std::span<const double> values,std::span<const double> weights,
                         std::uint64_t memory_budget=active_memory_budget());
    double percentile(double q) const;
    std::vector<double> percentiles(std::span<const double> queries) const;
private:
    std::vector<double> values_,centers_;
    double total_=0;
    std::uint64_t memory_budget_=0;
};
struct TriangleFrame {std::array<Vec2,3> coordinates{};double double_area=0;};
CLEANSEAMS_CORE_EXPORT TriangleFrame triangle_frame(const std::array<Vec3,3>& positions);
// A dense original-corner buffer, or a compact buffer indexed by strictly
// increasing original corner IDs. Views are consumed synchronously, not retained.
struct CornerUvView {std::span<const Vec2> values;std::span<const MeshIndex> corners{};};
struct QualityPolicy {
    bool exempt_input_degenerate=true;
    std::vector<TransferIndex> exempt_faces;
    std::vector<std::int8_t> reference_orientation;
    Revision uv_revision=0;
};
class CLEANSEAMS_CORE_EXPORT QualityGeometry {
public:
    explicit QualityGeometry(std::shared_ptr<const CutSurface>,std::uint64_t memory_budget=active_memory_budget());
    const CutSurface& cut() const {return *cut_;}
    const std::vector<TriangleFrame>& frames() const {return frames_;}
    const IndexRows& chart_triangles() const {return chart_triangles_;}
    double sliver_threshold() const {return sliver_threshold_;}
    std::uint64_t memory_budget() const {return memory_budget_;}
private:
    std::shared_ptr<const CutSurface> cut_;
    std::vector<TriangleFrame> frames_;
    IndexRows chart_triangles_;
    double sliver_threshold_=quality_nan;
    std::uint64_t memory_budget_=0;
};
struct TriangleQuality {
    std::array<MeshIndex,3> corners{};
    MeshIndex face=0,chart=0;
    Vec2 singular_values{quality_nan,quality_nan};
    double area3=0,area2=0,raw_area2=0;
    bool finite=false,sliver=false,exempt=false;
};
struct PolygonQuality {
    MeshIndex face=0,chart=0;
    double signed_area=0,raw_signed_area=0;
    bool finite=false,flipped=false,raw_flipped=false,exempt=false,degenerate_input=false;
};
struct CLEANSEAMS_CORE_EXPORT QualityReport {
    std::uint64_t triangles=0,flipped=0,folded=0,degenerate=0,nonfinite=0,slivers=0,
                  collapsed_uv=0,exempt=0,degenerate_input=0,flipped_exempt=0;
    double anisotropy_p50=quality_nan,anisotropy_p95=quality_nan,anisotropy_p99=quality_nan,
           anisotropy_max=quality_nan,area_p50=quality_nan,area_p95=quality_nan,
           area_p99=quality_nan,area_max=quality_nan,density=0,density_spread=1;
    // Finite JSON, with nulls and explicit path/reason records for unavailable fields.
    nlohmann::json to_json() const;
};
CLEANSEAMS_CORE_EXPORT std::vector<std::int8_t> chart_reference_orientation(const QualityGeometry&,CornerUvView);
class CLEANSEAMS_CORE_EXPORT QualitySnapshot {
public:
    QualitySnapshot(std::shared_ptr<const QualityGeometry>,CornerUvView,const QualityPolicy& policy={},
                    std::optional<std::span<const MeshIndex>> charts=std::nullopt);
    const std::vector<TriangleQuality>& triangles() const {return triangles_;}
    const std::vector<PolygonQuality>& polygons() const {return polygons_;}
    const QualityPolicy& policy() const {return policy_;}
    const QualityGeometry& geometry() const {return *geometry_;}
    QualityReport report(std::optional<MeshIndex> chart=std::nullopt,
                         std::optional<double> density=std::nullopt,bool per_chart=true) const;
    // Source ordering: UV/3D density, summed double 3D area, summed absolute double UV area.
    std::array<std::vector<double>,3> chart_density() const;
    std::vector<double> chart_anisotropy_p95() const;
private:
    std::shared_ptr<const QualityGeometry> geometry_;
    QualityPolicy policy_;
    std::vector<MeshIndex> selected_charts_;
    std::vector<TriangleQuality> triangles_;
    std::vector<PolygonQuality> polygons_;
    IndexRows chart_triangles_,chart_polygons_;
};
struct BoundaryCrossingOptions {
    std::optional<double> cell;
    bool single_face_slit=true,retain_pairs=false;
    // Exact crossings and endpoints in another segment's interior. Matching
    // endpoints/whole edges remain legal: a shared affine transform preserves
    // their contact after storage rounding. Uses orientation_sign's contract.
    bool nonconforming_contacts=false;
    // True stops enumeration; an incomplete report cannot certify a boundary.
    std::function<bool()> stopped;
    std::vector<TransferIndex> ignored_faces;
    std::uint64_t memory_budget=active_memory_budget();
    // Exact streaming enumeration is used if this acceleration allowance is
    // exhausted. Zero selects that bounded-memory path directly.
    std::uint64_t workspace_bytes=16777216ULL;
};
struct BoundaryCrossingPair {
    MeshIndex chart=0,first_corner=0,second_corner=0;
    auto operator<=>(const BoundaryCrossingPair&) const = default;
};
struct BoundaryCrossingReport {
    std::vector<std::uint64_t> counts;
    std::vector<std::uint8_t> testable;
    std::vector<MeshIndex> selected_charts,crossing_faces;
    std::vector<BoundaryCrossingPair> pairs;
    std::uint64_t bucket_entries=0,candidate_pairs=0;
    bool used_exhaustive=false,complete=true;
};
// Partial collinear contacts, T-junctions and proper crossings; matching
// endpoints/whole edges are false. Same numeric domain as orientation_sign.
CLEANSEAMS_CORE_EXPORT bool segments_have_nonconforming_contact(Vec2,Vec2,Vec2,Vec2);
CLEANSEAMS_CORE_EXPORT bool segments_cross(Vec2 p,Vec2 q,Vec2 r,Vec2 s,double eps=1e-12);
CLEANSEAMS_CORE_EXPORT BoundaryCrossingReport boundary_crossings(const CutSurface&,CornerUvView,
    const BoundaryCrossingOptions& options={},std::optional<std::span<const MeshIndex>> charts=std::nullopt);
struct CLEANSEAMS_CORE_EXPORT ChartValidity {
    std::vector<std::uint64_t> flipped,nonfinite,collapsed,crossings;
    std::vector<std::uint8_t> valid,measured;
    nlohmann::json to_json() const;
};
// Builds related measurements from one UV input. Raw chart validity and
// orientation-adjusted distortion remain distinct, matching the source helpers.
class CLEANSEAMS_CORE_EXPORT QualityEvaluation {
public:
    QualityEvaluation(std::shared_ptr<const QualityGeometry>,CornerUvView,const QualityPolicy& policy={},
                      std::optional<std::span<const MeshIndex>> charts=std::nullopt);
    const QualitySnapshot& snapshot() const {return snapshot_;}
    const ChartValidity& validity() const {return validity_;}
    const BoundaryCrossingReport& validation_crossings() const {return *validation_crossings_;}
    const BoundaryCrossingReport& metric_crossings() const {return *metric_crossings_;}
    bool chart_valid(MeshIndex chart) const;
    std::vector<MeshIndex> offending_faces(MeshIndex chart) const;
    nlohmann::json metrics() const;
private:
    QualitySnapshot snapshot_;
    ChartValidity validity_;
    IndexRows offending_faces_;
    std::shared_ptr<const BoundaryCrossingReport> metric_crossings_,validation_crossings_;
    bool whole_snapshot_=false;
};
// ------------------------------------------------ strict overlap (Phase B, T28)
// Positive-area interior overlap of UV triangles, the enhanced_v1 validity check that the
// compatibility crossing predicate above does not replace (Section 11.3). It is exact:
// orientation signs come from a filtered predicate with an exact expansion fallback, and two
// triangles overlap when no edge line of either separates them (closed half-planes), so
// triangles that only share an edge, a vertex or a collinear stretch of boundary (legal slit
// contacts) never overlap. No pair is skipped because its triangles share a sector or an
// original vertex. The broad phase is a per-chart uniform grid holding every cell an
// AABB touches; each AABB-overlapping pair is tested once, so the grid never hides a pair.
// Label of this policy in reports; the compatibility crossing metrics keep their own.
inline constexpr const char* strict_overlap_policy="strict_overlap.v1";
// Exact sign of orient(a, b, c) (+1 counter-clockwise, -1 clockwise, 0 collinear) for finite
// coordinates whose products neither overflow nor underflow; `exact` reports whether the
// double filter was inconclusive and the expansion fallback decided.
CLEANSEAMS_CORE_EXPORT int orientation_sign(Vec2 a,Vec2 b,Vec2 c,bool* exact=nullptr);
// True when the open interiors of the two triangles intersect (a degenerate triangle has an
// empty interior). Exact under the orientation_sign contract.
CLEANSEAMS_CORE_EXPORT bool triangle_interiors_intersect(const std::array<Vec2,3>& first,const std::array<Vec2,3>& second,
    std::uint64_t* exact_fallbacks=nullptr);
// Area of the intersection of two triangles (either orientation), in double precision.
CLEANSEAMS_CORE_EXPORT double triangle_intersection_area(const std::array<Vec2,3>& first,const std::array<Vec2,3>& second);
struct OverlapOptions {
    // An interior intersection whose area is at most this fraction of the smaller triangle's
    // UV area is a contact (counted in `contacts`), not an overlap. Explicit and relative:
    // exact contacts have zero area and never depend on it.
    double area_tolerance=1e-10;
    // Faces left out, as the validity policy's exemptions (input-degenerate faces are
    // added when `ignore_input_degenerate`).
    std::vector<TransferIndex> ignored_faces;
    bool ignore_input_degenerate=true;
    bool retain_pairs=true;
    std::size_t max_pairs=4096;
    std::uint64_t memory_budget=active_memory_budget();
};
struct OverlapPair {
    MeshIndex first_triangle=0,second_triangle=0,first_face=0,second_face=0,first_chart=0,second_chart=0;
    double area=0;
};
struct CLEANSEAMS_CORE_EXPORT OverlapReport {
    std::string policy=strict_overlap_policy;
    bool cross_chart=false;                 // true: the packing question (pairs of different charts)
    std::vector<std::uint64_t> counts;      // overlapping pairs per chart (a cross pair counts for both)
    std::vector<std::uint8_t> tested;       // per chart: every eligible triangle had finite corners
    std::vector<MeshIndex> selected_charts,overlapping_faces;
    std::vector<OverlapPair> pairs;         // the max_pairs overlapping pairs with the smallest triangle IDs
    // Broad phase (m5: a bounding volume hierarchy of the UV boxes per tested group): `cells`
    // counts its leaves and `cell_entries` the boxes in them; `candidate_pairs` are the pairs
    // whose box interiors intersect (every pair that can overlap), each tested exactly once.
    std::uint64_t triangles=0,degenerate=0,untestable=0,cells=0,cell_entries=0,candidate_pairs=0,exact_tests=0,
                  exact_fallbacks=0,contacts=0,overlaps=0;
    bool pairs_truncated=false;
    nlohmann::json to_json() const;
};
// The strict check's broad phase (m5), exposed for tests: every pair (i < j) of the boxes whose
// open interiors intersect (low < other high on both axes), each once, sorted. A bounding
// volume hierarchy over the boxes: O(n log n + pairs) for charts of any shape (long strips
// included); never a heuristic that could leave out a pair that can overlap (Section 11.3).
struct UvBox {Vec2 low{},high{};};
CLEANSEAMS_CORE_EXPORT std::vector<std::array<std::size_t,2>> uv_box_pairs(std::span<const UvBox> boxes,
    std::uint64_t memory_budget=active_memory_budget());
// Chart-local check: pairs of triangles of the same chart (intrachart validity).
CLEANSEAMS_CORE_EXPORT OverlapReport chart_overlaps(const QualityGeometry&,CornerUvView,const OverlapOptions& options={},
    std::optional<std::span<const MeshIndex>> charts=std::nullopt);
// Pairs of triangles of different selected charts. Separately solved, unpacked charts share
// UV space, so this answers the packing question only; it is never an intrachart invalidity.
CLEANSEAMS_CORE_EXPORT OverlapReport cross_chart_overlaps(const QualityGeometry&,CornerUvView,const OverlapOptions& options={},
    std::optional<std::span<const MeshIndex>> charts=std::nullopt);
// Polyline measures (source quality.resample_polyline, gaussian_smooth_open,
// arc_jaggedness). Input contract: every point finite. The source would propagate NaN
// through NumPy; native refuses with invalid_input instead. No production caller passes
// a non-finite point (contour arcs come from validated maps; Gate A F32).
struct ResampledPolyline {std::vector<Vec2> points;std::vector<double> arc;};
struct ArcJaggedness {double rms=0,maximum=0,h=0;std::size_t samples=0;};
struct ChainTurns {std::uint64_t turns=0,lane_changes=0;std::int64_t edges=-1;};
CLEANSEAMS_CORE_EXPORT ResampledPolyline resample_polyline(std::span<const Vec2>,double step,
    std::uint64_t memory_budget=active_memory_budget());
CLEANSEAMS_CORE_EXPORT std::vector<Vec2> gaussian_smooth_open(std::span<const Vec2>,double sigma_samples,
    std::uint64_t memory_budget=active_memory_budget());
CLEANSEAMS_CORE_EXPORT ArcJaggedness arc_jaggedness(std::span<const Vec2>,std::optional<double> scale=std::nullopt,
    bool straight=false,std::uint64_t memory_budget=active_memory_budget());
CLEANSEAMS_CORE_EXPORT ChainTurns chain_turns(const MeshSnapshot&,std::span<const TransferIndex> vertices);
}
