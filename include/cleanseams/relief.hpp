// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
// Phase B (T32, Section 14): sector-aware hotspot relief slits for enhanced_v1. A relief
// slit is a dead-end cut from a chart's boundary (or a slit side) into a distorted
// region: it adds freedom without adding an island. Candidates are ordinary proposals of
// the T28 transaction; nothing here accepts, measures or commits a state on its own.
#include "cleanseams/metrics.hpp"
#include "cleanseams/proposals.hpp"
namespace cleanseams {
inline constexpr const char* hotspot_relief_version="hotspot_relief.v1";
// Explicit feature-search budgets and thresholds of the relief stage (Section 14.2 initial
// enhanced defaults). They bound this feature's own search; no legacy search is reduced.
struct CLEANSEAMS_CORE_EXPORT ReliefSettings {
    std::int64_t max_hot_regions=4;             // separated hot regions per chart
    std::int64_t max_candidates_per_pass=12;    // distinct slit paths per affected chart and pass
    std::int64_t passes=2;
    std::vector<double> reaches={0.25,0.5,0.75,1.0}; // fractions of the available physical reach
    std::int64_t targets_per_region=3;          // the worst point and candidates near it
    // Excess thresholds of the sector field. ln(1.25) on log_area is a 25% AREA change
    // (an 11.8% linear checker change); on log_anisotropy a singular value ratio of 1.25.
    double log_area_threshold=0.22314355131420976;
    double log_anisotropy_threshold=0.22314355131420976;
    std::int64_t smoothing_steps=2;             // conservative smoothing steps on the cut adjacency
    double separation=0.25;                     // hot-region separation: this fraction of sqrt(chart 3D area), by surface distance
    nlohmann::json to_json() const;
};
struct HotRegion {
    MeshIndex sector=0,vertex=0;                // centre: the sector of largest excess not yet covered
    bool boundary=false;                        // the centre lies on the chart boundary or a slit side
    double excess=0,weighted_excess=0,log_area=0,log_anisotropy=0,mass=0;
};
struct ReliefTarget {
    MeshIndex region=0,rank=0;                  // region index, rank within the region (0: the worst point)
    MeshIndex tip_sector=0,tip_vertex=0;        // interior tip of the slit
    std::optional<MeshIndex> start_sector;      // fixed boundary start (boundary regions); none: cheapest boundary sector
    double score=0,excess=0,suitability=0,seam_distance=0,depth=0;
};
struct ReliefPath {
    MeshIndex target=0;                         // index into ReliefPlan::targets
    std::vector<MeshIndex> sectors,vertices,edges; // boundary start first, tip last; original mesh edges
    double length=0,available=0;                // physical length of this slit and of the full path to the tip
    std::vector<double> reaches;                // requested reach fractions that snapped to this slit
};
struct CLEANSEAMS_CORE_EXPORT ReliefPlan {
    MeshIndex chart=0;
    nlohmann::json field=nlohmann::json::object(); // SectorDistortionField summary (units included)
    double chart_area=0,separation=0;
    std::vector<HotRegion> regions;
    std::vector<ReliefTarget> targets;
    std::vector<ReliefPath> paths;              // selected distinct slits: target by target in rank order, shortest first
    std::uint64_t generated=0,duplicates=0,shared_tips=0,truncated=0,unreachable=0;
    std::string skipped;                        // why no slit was produced ("" otherwise)
    nlohmann::json to_json() const;
};
struct ReliefInput {
    std::shared_ptr<const CutSurface> surface;
    std::span<const Vec2> uv;                   // map of `chart` on `surface` (per original corner)
    MeshIndex chart=0;
    std::span<const std::uint8_t> exempt;       // the run's exemption policy (faces)
    const HardConstraints* constraints=nullptr; // forbidden edges, locked seam tips and pins are respected
    std::span<const double> edge_cost;          // nonnegative routing heuristic per edge; empty: edge length
};
// The relief candidates of one disk chart (Section 14.2): the sector field, up to
// max_hot_regions hot regions separated by surface distance, tips (interior regions: the
// worst interior sectors near the centre, by excess x soft suitability x distance from
// existing slits; boundary regions: the deepest interior sector from the hot boundary
// sector), cheapest real-edge paths over the chart's sectors from an eligible boundary or
// slit-side sector, and the requested reaches of each path, deduplicated.
CLEANSEAMS_CORE_EXPORT ReliefPlan plan_relief(const ReliefInput&,const ReliefSettings& settings={});
// The relief stage of the Section 21 schedule: `passes` generators sharing one state. Each
// proposes relief slits for the incumbent's over-budget disk charts and, for every broad
// over-budget split of the baseline recorded in ProposalContext::baseline_splits, the
// split's welding plus relief slits on the reunited chart (relief before the split). Charts
// within budget get nothing. A target tried on an unchanged incumbent is not proposed again.
CLEANSEAMS_CORE_EXPORT std::vector<ProposalGenerator> hotspot_relief_generators(const PipelineOptions&,ReliefSettings settings={});
}
