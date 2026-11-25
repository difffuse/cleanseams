// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
// Phase B physical-distance policy and tessellation robustness (Section 20, T39). Used by
// policy enhanced_v1 only; legacy_010 and native_baseline keep their hop and face-count
// behavior and never reach any of it.
//
//   operation start (candidate zero)  -> ScaleReference: areas, median edge lengths and median
//                                        face areas per reference chart, the legacy corridor
//                                        rings, frozen for the whole pass (never recomputed as
//                                        trial cuts change the charts)
//   PhysicalScale (enhanced.scale)    -> corridor radii, smoothing steps and region floors in
//                                        object-space units (mode physical), or the legacy
//                                        hop/face-count semantics (mode hops, documented)
//   geodesic corridors (corridor.cpp) -> surface distances along the patch's own edges: never
//                                        Euclidean, so a nearby sheet, a folded layer or another
//                                        component is never a shortcut
//
// Exact topology never depends on this policy: a hole is opened whatever its edge count or
// size (open_non_disk_charts and the transaction's disk check have no size threshold).
#include "cleanseams/endpoints.hpp"
#include <cmath>
namespace cleanseams {
inline constexpr const char* physical_scale_version="physical_scale.v1";
// The scale of the run origin, frozen once per enhancement pass (Section 20: "during a
// candidate search, do not recompute that reference length as cuts change").
struct CLEANSEAMS_CORE_EXPORT ScaleReference {
    std::string geometry_key;                  // mesh geometry the reference belongs to
    Revision revision=0;                       // the state it was frozen on (candidate zero)
    std::int64_t corridor_rings=0;             // the legacy UI radius captured at the start
    double area=0;                             // non-exempt 3D area of the mesh
    double median_edge_length=0;               // over every surface edge of the mesh
    double median_face_area=0;                 // over non-exempt faces with positive area
    std::vector<MeshIndex> chart_of_face;      // reference chart of each face (candidate zero's charts)
    std::vector<double> face_area;             // per face: 3D area (area-weighted reference charts)
    std::vector<double> chart_area;            // per reference chart: non-exempt 3D area
    std::vector<double> chart_edge_length;     // per reference chart: median edge length of its faces' edges
    std::vector<double> chart_face_area;       // per reference chart: median non-exempt positive face area
    nlohmann::json to_json() const;            // without the per-face array
};
// Freezes the reference of `origin` (candidate zero) under the run's exemption policy.
CLEANSEAMS_CORE_EXPORT ScaleReference freeze_scale_reference(const CutSurface& origin,std::span<const std::uint8_t> exempt,
    const PipelineOptions&,Revision revision=0);
// A resolved length of the policy and where it came from.
struct CLEANSEAMS_CORE_EXPORT ScaleLength {
    double value=0;
    std::string source;                        // explicit, reference, hops
    MeshIndex reference_chart=0;
    nlohmann::json to_json() const;
};
// A smoothing operator's step count under the policy. Diffusion smoothers (averaging with the
// neighbours) spread like sqrt(steps) edge lengths, ring operators (growth by one face ring)
// like steps edge lengths; their default counts are defined at one step per reference length.
struct CLEANSEAMS_CORE_EXPORT SmoothingSteps {
    std::int64_t default_steps=0,steps=0;
    double radius=0;                           // physical extent: sqrt(steps) or steps local lengths
    double local_length=0,unit=0;              // the chart's frozen median edge length and the length unit
    bool bounded=false;                        // max_steps capped the count
    nlohmann::json to_json() const;
};
// A legacy face-count floor as a region-preservation limit: an area (physical) or the count (hops).
struct CLEANSEAMS_CORE_EXPORT RegionFloor {
    std::int64_t faces=0;
    double area=0;                             // physical: faces * reference_area, or * median face area * (unit / local length)^2
    bool physical=true;
    std::string source="reference";            // explicit (reference_area), reference_length or reference
    // Relative tolerance: a region of exactly `faces` equal faces is never below its floor.
    bool below(double region_area,std::uint64_t region_faces) const;
    nlohmann::json to_json() const;
};
class CLEANSEAMS_CORE_EXPORT PhysicalScale {
public:
    PhysicalScale(ScaleReference reference,ScaleOptions options);
    // Freezes the reference of `origin` with the run's options (enhanced.scale).
    static PhysicalScale freeze(const CutSurface& origin,std::span<const std::uint8_t> exempt,const PipelineOptions&,Revision revision=0);
    bool physical() const {return options_.mode!="hops";}
    const ScaleReference& reference() const {return reference_;}
    const ScaleOptions& options() const {return options_;}
    // The reference chart of a face set: the reference chart holding the largest share of its
    // 3D area (area-weighted, ties: the lowest chart); faces outside the mesh are refused.
    MeshIndex reference_chart(std::span<const MeshIndex> faces) const;
    // The frozen median edge length of the face set's reference chart (its own resolution).
    ScaleLength local_length(std::span<const MeshIndex> faces) const;
    // The length unit: enhanced.scale.reference_length, or the local length.
    ScaleLength unit(std::span<const MeshIndex> faces) const;
    // Corridor radius around a face set (Section 20): corridor_radius, or corridor_scale (0:
    // the captured paths.corridor_rings) times the unit. In hops mode: max(1, rings) hop rings.
    ScaleLength corridor_radius(std::span<const MeshIndex> faces) const;
    // Step counts of a diffusion smoother or a ring operator on the face set.
    SmoothingSteps diffusion(std::int64_t default_steps,std::span<const MeshIndex> faces) const;
    SmoothingSteps rings(std::int64_t default_passes,std::span<const MeshIndex> faces) const;
    // A legacy face-count floor as an area on the face set's reference chart.
    RegionFloor floor(std::int64_t faces,std::span<const MeshIndex> region) const;
    // The length quantum of decisions on the face set: 2^(floor(log2 unit) - 32), about 2e-10
    // of the unit (0 in hops mode: exact legacy comparisons). Lengths compared on this grid
    // (and sums of them, exact on a power-of-two grid) keep exact geometric ties exact, so the
    // rounding noise of a rigid motion or a uniform scale never decides a symmetric tie.
    double length_quantum(std::span<const MeshIndex> faces) const;
    nlohmann::json to_json() const;
private:
    ScaleReference reference_;
    ScaleOptions options_;
};
// The pass's frozen scale: ProposalContext::scale, frozen by begin_enhancement at the
// operation start, or here from `state` when a caller built the context without it.
CLEANSEAMS_CORE_EXPORT const PhysicalScale& context_scale(ProposalContext&,const RunState& state);
// report.enhanced.scale: the frozen policy and the normalized seam geometry (Section 20 units:
// lengths over sqrt(area)) of candidate zero and of the retained state.
CLEANSEAMS_CORE_EXPORT nlohmann::json physical_scale_report(const ProposalContext&,const RunState& origin,const RunState& retained);
// ------------------------------------------------------------ geodesic corridors (routing/corridor.cpp)
// Dimensionless decision scores (ratios, normalized distances) are compared on this grid.
inline constexpr double score_quantum=0x1p-32;
// `x` on the grid of `quantum` (unchanged when quantum is 0).
inline double snap_to_quantum(double x,double quantum) {return quantum>0 && std::isfinite(x)?std::nearbyint(x/quantum)*quantum:x;}
enum class DistanceMetric {physical,hops};
CLEANSEAMS_CORE_EXPORT std::string_view distance_metric_name(DistanceMetric);
// Distances from `seeds` along the edges of a face patch (every face of the edge is in the
// patch): 3D edge lengths (physical) or edge counts (hops). Surface distance only: vertices of
// another component, of a nearby or folded sheet not joined by patch edges are unreachable
// (+inf) or reached along the surface, never through space. Seeds must be patch vertices.
CLEANSEAMS_CORE_EXPORT std::vector<double> patch_distances(const MeshSnapshot&,std::span<const std::uint8_t> patch_faces,
    std::span<const MeshIndex> seeds,DistanceMetric);
// The band of a patch distance field within `radius` (none: every reached patch vertex),
// restricted to usable edges and vertices (the PatchCorridor of T30 over any metric).
CLEANSEAMS_CORE_EXPORT PatchCorridor corridor_band(const MeshSnapshot&,std::span<const std::uint8_t> patch_faces,
    std::span<const std::uint8_t> edge_usable,std::span<const std::uint8_t> vertex_usable,std::span<const double> distance,
    std::optional<double> radius);
// The deterministic widening schedule of a search (Sections 13.2 and 20): the whole patch when
// it has at most `vertex_budget` vertices; otherwise bands of radius first * 2^k (never below
// the first radius, never shrinking as cuts change), doubled while the next band still holds at
// most `vertex_budget` patch vertices (the first band is always kept).
struct CLEANSEAMS_CORE_EXPORT CorridorSchedule {
    std::string metric;                        // physical or hops
    std::vector<std::optional<double>> radii;  // none: the whole patch
    std::vector<std::uint64_t> band_vertices;
    std::uint64_t patch_vertices=0;
    std::optional<double> search_radius() const {return radii.empty()?std::nullopt:radii.back();}
    nlohmann::json to_json() const;
};
CLEANSEAMS_CORE_EXPORT CorridorSchedule corridor_schedule(std::span<const double> distance,std::uint64_t patch_vertices,
    double first_radius,std::uint64_t vertex_budget,DistanceMetric);
}
