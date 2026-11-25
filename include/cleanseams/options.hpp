// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
#include "cleanseams/types.hpp"
#include <array>
#include <cmath>
#include <optional>
#include <string>
#include <nlohmann/json.hpp>
namespace cleanseams {
using Vec2=std::array<double,2>;
using Vec3=std::array<double,3>;
using IslandBand=std::array<std::int64_t,2>;
struct LayoutOptions {
    std::string seeding = "auto";
    std::int64_t seed_density = 400;
    std::int64_t max_patches = 4000;
    std::int64_t segment_iterations = 3;
    bool relaxed_motorcycles = false;
    double feature_angle = 60.0;
    bool prefer_features = false;
    bool lock_features = false;
    double feature_lock_fraction = 0.5;
    std::int64_t min_feature_chain = 6;
    std::int64_t feature_scale = 2;
    double feature_nms_ratio = 0.6;
    std::int64_t max_extend = 8;
    std::optional<double> max_defect = std::nullopt;
    double min_island_fraction = 0.002;
    std::int64_t min_island_faces = 8;
    double strip_isoperimetric = 12.0;
    double absorb_relax = 2.0;
    bool strict_budget = false;
    double max_anisotropy_p95 = 1.6;
    double max_area_p95 = 1.8;
    double max_anisotropy_p99 = 4.0;
    double max_anisotropy_max = 25.0;
    double max_faces_fraction = 1.0;
    std::int64_t max_rounds = 60;
    std::int64_t max_merge_edges = 0;
    bool evaluate_true_distortion = true;
    std::string solver = "arap";
    std::string priority = "compact";
    double min_fusion = 0.15;
    std::int64_t max_split_passes = 6;
    std::int64_t split_min_faces = 40;
    std::int64_t arap_iterations = 6;
    bool use_exposure = true;
    std::int64_t exposure_K = 64;
    double exposure_weight = 1.0;
    std::optional<Vec3> hide_dir = std::nullopt;
    double revolution_residual = 0.04;
    double max_isoperimetric_growth = 1.1;
    std::int64_t isoperimetric_min_faces = 40;
    bool arc_moves = true;
    bool spur_removal = true;
    std::int64_t strip_max_width = 3;
};
struct PathOptions {
    std::int64_t corridor_rings = 2;
    std::int64_t corridor_max_rings = 4;
    double w_turn = 4.0;
    double w_flow = 2.0;
    double w_drift = 3.0;
    double w_convex = 0.5;
    double w_exposure = 0.5;
    double anchor_tangent_weight = 0.5;
    std::vector<Vec3> variants = {{4.,2.,3.},{2.,4.,3.},{6.,2.,1.5},{4.,2.,0.5}};
    bool require_improvement = true;
    bool require_budget = false;
};
struct BoundaryOptions {
    std::string mode = "auto";
    std::vector<double> intensities = {1.,.5,.25,.125};
    double corner_angle = 60.0;
    double straight_angle = 150.0;
    double smooth_radius = 3.0;
    double max_displacement = 1.5;
    double straight_tolerance = 1.5;
    double p95_growth = 1.05;
    double max_growth = 1.1;
    double abs_anisotropy_p95 = 2.0;
    double abs_area_p95 = 2.5;
    std::int64_t min_closed_anchors = 4;
    std::int64_t arap_iterations = 8;
};
struct PackOptions {
    std::int64_t texture_width = 2048;
    std::int64_t texture_height = 2048;
    double padding_px = 4.0;
    double border_px = 2.0;
    std::vector<std::int64_t> tiles = {1001};
    bool keep_orientation = true;
    bool allow_mirroring = false;
    bool equalize_density = true;
    double time_budget = 0.0;
    std::int64_t trials = 8;
    std::int64_t threads = 1;
};
struct SolverOptions {
    std::string method = "arap";
    std::int64_t iterations = 10;
    std::string backend = "superlu";
};
// Phase B (Section 11.2): the weights of the versioned seam objective of the common final
// ordering, J = w_length*L + w_exposure*V + w_turn*T + w_fragment*F. The defaults are the
// initial engineering weights frozen with seam_objective.v1; not learned parameters.
// Reports name the version and whether these are exactly its frozen weights.
// GB (M3): two versions of the ordering over the same J. seam_objective.v1 (T28-T40) compares
// J within the relative tie band only. seam_objective.v2 (the enhanced_v1 default) also treats
// as a tie any J difference within the sampling uncertainty of the raw exposure term
// (exposure_uncertainty.v1, include/cleanseams/proposals.hpp): the sampled all-views exposure
// no longer decides between cut states it cannot tell apart, and a tie keeps the incumbent.
inline constexpr const char* seam_objective_v1 = "seam_objective.v1";
inline constexpr const char* seam_objective_v2 = "seam_objective.v2";
inline constexpr const char* seam_objective_version = seam_objective_v2; // the default ordering
// GB (visual review V2): the practical-effect floor of seam_objective.v2. Wherever the v2
// significance margin applies, a comparison decided by J (both states admissible, the same
// strict-defect, compliance and over-budget status) adopts the challenger only when it lowers J
// by more than max(U, practical_floor * |J_incumbent|): the exposure-sampling margin U alone let
// ~0.1 % changes through (seam churn without a visible benefit). A change that only removes seam
// edges (a subset of the incumbent's cuts: spur removal) has no floor. 0 disables the floor.
// Budget, compliance and strict-defect decisions never see it.
struct SignificanceOptions {
    double practical_floor = 0.004;
};
struct ObjectiveOptions {
    std::string version = seam_objective_v2; // seam_objective.v1 reproduces T40's decisions
    double w_length = 1.0;          // L: internal seam length / sqrt(relevant 3D area)
    double w_exposure = 2.0;        // V: sum(edge length * raw exposure) / sqrt(area)
    double w_turn = 0.05;           // T: sum(local mean edge length * turn penalty) / sqrt(area)
    double w_fragment = 0.10;       // F: small-fragment area fraction + chart_count_weight * additional charts
    double chart_count_weight = 0.01;
    double tie_band = 1e-9;         // relative; scores within it keep the incumbent
    SignificanceOptions significance = {}; // GB (V2): seam_objective.v2's practical-effect floor
};
// Phase B painted preferences and the enhanced routing heuristic (Section 12, T29):
//   edge_cost = length * (1 + w_exposure*exposure^2 + w_view*visible + w_avoid*avoid)
//                      / (1 + w_allow*allow + w_crease*crease)
// with the raw measured exposure, visibility from an explicit viewer only, painted weights
// in [0, 1] and geometric crease desirability. Nonnegative weights keep the denominator at
// least 1. Initial engineering values, frozen with the objective before held-out evaluation.
struct PreferenceOptions {
    double w_exposure = 2.0;
    double w_view = 2.0;
    double w_avoid = 4.0;
    double w_allow = 4.0;
    double w_crease = 1.0;
    // Host attribute names read for the painted constraints and preferences: the documented
    // names or explicit user aliases, never inferred from other names. Empty disables one.
    std::string required_attribute = "cleanseams_required";   // EDGE BOOLEAN: must be a seam
    std::string forbidden_attribute = "cleanseams_forbidden"; // EDGE BOOLEAN: never an internal seam
    std::string allow_attribute = "cleanseams_allow";         // POINT FLOAT [0, 1]: soft preference
    std::string avoid_attribute = "cleanseams_avoid";         // POINT FLOAT [0, 1]: soft preference
};
// Sliding seam endpoints on tube-like patches (Section 13, T30/T31). A slide route costs the
// enhanced routing heuristic per edge plus w_turn * local mean edge length * turn penalty at
// each interior vertex, w_tangent * local rim edge length * the terminal turn at each end
// and w_slide * the distance an endpoint moves along its rim. Initial engineering values.
struct EndpointSlideOptions {
    bool enabled = true;                   // T40: the Section 21 sliding stage (ablation switch)
    double max_slide_distance = 0.0;       // object-space distance along a rim an endpoint may move; 0 disables the limit
    std::int64_t vertex_budget = 20000;    // whole-patch search up to this many patch vertices; beyond it a widening surface-distance band
    double w_turn = 0.5;
    double w_tangent = 0.5;
    double w_slide = 0.25;
};
// Phase B preference-aware region proposals (Section 15, T34). Dual crossing weight of the
// face pair across edge e (weighted region growth, not a min-cut):
//   W = centroid_distance * (1 + w_normal*normal_change + w_preference*strength(e)^2)
// strength(e) in [0, 1] grows where a seam is welcome (cheap to cut = expensive to cross).
// w_preference = 0 is the retained baseline growth (the legacy dual weights).
struct RegionOptions {
    // Region proposals can change the placement of visible feature seams. Opt in per run.
    bool enabled = false;
    double w_normal = 4.0;
    double w_preference = 4.0;
    std::int64_t lloyd_steps = 3;        // Lloyd updates of the grown regions
    std::int64_t max_targets = 4;        // chart pairs and over-budget charts re-partitioned per pass
    // GB (visual review V3): the incumbent's seams on hard feature edges (|dihedral| at least the
    // baseline seam cost's feature angle, the edges its routing discounts as hard; never a painted
    // forbidden edge) are barriers of the region dual, like required and locked cuts, and stay cut
    // wherever the regrown regions meet across them: a region proposal never moves a seam off a
    // sharp edge onto a smooth face. A dead end of one inside a regrown region is not forced to stay.
    bool keep_feature_seams = true;
};
// Phase B generalized ring proposals (Section 16, T35): physical distance from a patch's
// opening sampled at section_samples quantiles in [low_quantile, high_quantile] of its
// maximum; simple closed level-set loops ranked by constriction, crease, exposure, balance
// and compactness. Loops below min_confidence (max of constriction and crease evidence) are
// reported, not proposed; at most max_candidates loops per chart on max_charts charts.
struct RingOptions {
    bool enabled = true;
    std::int64_t section_samples = 12;
    double low_quantile = 0.15;
    double high_quantile = 0.90;
    std::int64_t max_candidates = 4;
    std::int64_t max_charts = 4;
    double min_confidence = 0.15;
};
// Phase B confidence-gated cavity proposals (Section 17, T36): up to max_samples
// area-distributed faces of a chart cast rays_per_sample deterministic cone rays
// (cone_half_angle_deg about the face normal) reaching reach_fraction of their connected
// component's bounding-box diagonal; faces blocked at least occlusion_threshold of the time
// seed recess cores, refined on the original mesh at level_samples height levels. The area
// fractions (of the chart), depth over opening diameter and confidence are initial detector
// filters: they only generate proposals and are never exemptions from validation. Deep
// cavities (slit_depth_ratio) need an interior slit; at most max_cavities cavities per pass
// on charts of at least min_chart_faces faces.
struct CavityOptions {
    bool enabled = true;
    std::int64_t max_samples = 20000;
    std::int64_t rays_per_sample = 8;
    double cone_half_angle_deg = 75.0;
    double reach_fraction = 0.08;
    double occlusion_threshold = 0.25;
    std::int64_t level_samples = 16;
    double min_area_fraction = 0.01;
    double max_area_fraction = 0.30;
    double min_depth_ratio = 0.2;
    double slit_depth_ratio = 0.35;
    double min_confidence = 0.25;
    std::int64_t max_cavities = 4;
    std::int64_t min_chart_faces = 16;
};
// Phase B thin-shell proposals (Section 18, T37): a connected component of at least
// min_faces faces is nominated when its local opposite-side thickness is below thin_ratio of
// its middle principal extent (area-weighted PCA); up to max_samples front/back faces test
// the opposite side, which must answer at least min_correspondence of their area. Label
// smoothing passes, speckle repair below speckle_fraction of the area, the confidence gate
// and at most max_candidates alternatives per component, candidate zero included.
struct ShellOptions {
    bool enabled = true;
    double thin_ratio = 0.15;
    std::int64_t max_candidates = 3;
    std::int64_t max_samples = 20000;
    double min_correspondence = 0.7;
    double min_confidence = 0.5;
    double speckle_fraction = 0.03;
    std::int64_t smoothing_passes = 3;
    std::int64_t min_faces = 16;
};
// Phase B operational symmetry (Section 19, T38): an explicit policy, never an environment
// variable. mode: asymmetric (view-optimal: moves are not paired), strict (every move is paired
// with its mirror image as one transaction on an exact combinatorial correspondence), or
// approximate (metric-only diagnostics; no paired moves). on_conflict: reject_pair keeps the
// incumbent, disable_for_proposal evaluates the unpaired proposal with an explicit label. The
// plane is explicit (plane_normal, plane_offset; object space) or detected (T14). Tolerances are
// fractions of the bounding-box diagonal: map_tolerance finds counterparts (T14 mirror maps),
// geometry_tolerance validates them for strict pairing (above map_tolerance it has no further
// effect: counterparts are only searched within map_tolerance). max_candidates bounds the
// symmetrization alternatives of the Section 21 symmetry stage.
struct SymmetryOptions {
    std::string mode = "asymmetric";
    std::string on_conflict = "reject_pair";
    std::optional<Vec3> plane_normal = std::nullopt;
    double plane_offset = 0.0;
    double map_tolerance = 0.003;
    double geometry_tolerance = 1e-6;
    std::int64_t max_candidates = 3;
};
// Phase B physical-distance policy (Section 20, T39). `physical` gives every enhanced
// operation explicit object-space units, frozen at the operation start from candidate zero:
// geodesic corridor bands of radius corridor_radius (or corridor_scale, 0: the legacy
// paths.corridor_rings captured at the start, times the reference length), smoothing whose
// step counts are defined at one step per reference length, and region floors measured as
// area (a legacy face-count floor times reference_area, or 0: the reference chart's median
// face area at the reference length's resolution). The reference length is reference_length,
// or 0: the frozen median edge length of the reference chart (never recomputed as trial cuts
// change; supply the same explicit values to compare different tessellations). `hops` keeps the legacy hop-ring corridors, face-count
// floors and fixed step counts, for comparisons with native_baseline. Detector distances stay
// normalized by component or chart size. Holes are always opened exactly, whatever their size.
struct ScaleOptions {
    std::string mode = "physical";         // physical | hops
    double reference_length = 0.0;         // object-space length unit; 0: frozen median edge length of the reference chart
    double reference_area = 0.0;           // object-space area of one face of a face-count floor; 0: median face area at the unit
    double corridor_radius = 0.0;          // object-space corridor radius; 0: corridor_scale * reference length
    double corridor_scale = 0.0;           // corridor radius in reference lengths; 0: paths.corridor_rings captured at the start
    std::int64_t max_steps = 64;           // bound on smoothing steps or passes derived from physical lengths
};
// T40: switches of the Section 21 stages that have no settings record of their own (their
// budgets and thresholds are versioned native records: hotspot_relief.v1, seam_necessity.v1).
struct ReliefOptions {
    bool enabled = true;                   // T32 sector-aware hotspot relief before the broad split
};
struct SpurRemovalOptions {
    bool enabled = true;                   // T33 validated unnecessary-spur removal (the legacy spur cleanup is unchanged)
};
// GB (Section 11.3 policy): the validated strict-overlap repair of candidate zero's strict-only
// defects (native record strict_overlap_repair.v1). Off, the defects are reported unrepaired.
struct StrictRepairOptions {
    bool enabled = true;
};
// GB (M5): the Section 21 step 1 constraint repair (T29) switch, for the painted ablations.
// Off, a candidate zero that violates painted hard constraints gets no repair proposal (the
// constraints are never relaxed: a run whose kept state violates them is refused).
struct ConstraintRepairOptions {
    bool enabled = true;
};
// T40 bounded alternative search (Section 21). Next to the per-chart bounds of the pass
// (max_candidates_per_chart, max_solves_per_chart) every feature (generator) has its own
// share per origin chart, so one feature cannot use up a chart's budget before the later
// stages of the schedule see it. prune_dominated refuses, before its solve, a candidate
// whose exact objective (J depends on the cuts only) proves it cannot beat a compliant,
// admissible incumbent or round leader; estimates only order the queue.
struct SearchOptions {
    std::int64_t feature_candidates_per_chart = 4;
    std::int64_t feature_solves_per_chart = 12;
    bool prune_dominated = true;
};
// GB (visual review V1): every new or changed seam chain of a compound proposal (region, ring,
// cavity and shell alternatives) goes through the baseline's chain refinement (refine_chain /
// refine_cycle: the same corridors, turn, flow, drift, convexity and crease-discounted exposure
// terms) inside the proposal's target faces, before the transaction evaluates the proposal.
// Off, the raw region borders (triangulation staircases) are proposed as before.
struct BoundaryRefinementOptions {
    bool enabled = true;
};
// Phase B enhanced policy (enhanced_v1 only; every other policy ignores this record).
struct EnhancedOptions {
    ObjectiveOptions objective = {};
    PreferenceOptions preferences = {};
    EndpointSlideOptions endpoint_slide = {};
    RegionOptions regions = {};
    RingOptions rings = {};
    CavityOptions cavities = {};
    ShellOptions shells = {};
    SymmetryOptions symmetry = {};
    ScaleOptions scale = {};
    ReliefOptions relief = {};
    SpurRemovalOptions spur_removal = {};
    StrictRepairOptions strict_repair = {};
    ConstraintRepairOptions constraint_repair = {};
    SearchOptions search = {};
    BoundaryRefinementOptions boundary_refinement = {}; // GB (V1)
    // Section 21 bounds of the enhancement pass, per affected chart; existing legacy work
    // is not counted.
    std::int64_t max_candidates_per_chart = 8;
    std::int64_t max_solves_per_chart = 24;
    // Explicit pins must hold within this absolute + relative tolerance after a solve.
    double pin_tolerance = 1e-9;
    // Strict overlap: intersections at most this fraction of the smaller triangle's UV area
    // are contacts, not overlaps (exact contacts have zero area).
    double overlap_area_tolerance = 1e-10;
};
struct UiSettings {
    std::string profile = "hard_surface";
    bool generate = true;
    bool protect_existing = true;
    double feature_angle = 60.0;
    bool refine_paths = true;
    std::string contour_mode = "auto";
    std::string budget = "normal";
    double min_island = 0.2;
    std::int64_t solver_iterations = 10;
    std::int64_t min_island_faces = 8;
    bool target_enabled = false;
    std::int64_t target_min = 3;
    std::int64_t target_max = 8;
    double deadline_s = 0.0;
    std::string hide_from = "none";
    Vec3 hide_direction = {0, -1, 0};
    bool show_advanced = false;
    double w_turn = 4.0;
    double w_flow = 2.0;
    double w_drift = 3.0;
    std::int64_t corridor_rings = 2;
    double smooth_radius = 3.0;
    double max_displacement = 1.5;
    double min_fusion = 0.15;
    std::string pack_after = "none";
    bool keep_orientation = true;
    std::int64_t pack_width = 2048;
    std::int64_t pack_height = 2048;
    double pack_padding = 4;
    double pack_border = 2;
    std::string pack_tiles = "1001";
    std::string last_status = "";
    std::string last_report = "";
    // T41 enhanced controls (Blender panel and `--ui-options`): resolve_ui_options maps them
    // into PipelineOptions::policy and ::enhanced (inert under native_baseline). The defaults
    // are the native product defaults, so a project saved before these fields existed reads
    // policy native_baseline and the enhanced record's own defaults. Metadata, sections and
    // enum items: src/engine/ui_metadata.inc (Blender stores enums by item position: append
    // new items, never reorder them).
    std::string policy = "native_baseline";
    bool show_enhanced = false;
    bool enhanced_sliding = EndpointSlideOptions{}.enabled;
    double max_slide_distance = EndpointSlideOptions{}.max_slide_distance;
    bool enhanced_relief = ReliefOptions{}.enabled;
    bool enhanced_regions = RegionOptions{}.enabled;
    bool enhanced_rings = RingOptions{}.enabled;
    std::int64_t ring_confidence = std::lround(RingOptions{}.min_confidence*100);        // percent: Blender stores floats in float32, integers exactly
    bool enhanced_cavities = CavityOptions{}.enabled;
    std::int64_t cavity_confidence = std::lround(CavityOptions{}.min_confidence*100);    // percent
    bool enhanced_shells = ShellOptions{}.enabled;
    std::int64_t shell_confidence = std::lround(ShellOptions{}.min_confidence*100);      // percent
    bool enhanced_spur_removal = SpurRemovalOptions{}.enabled;
    std::string symmetry_mode = SymmetryOptions{}.mode;
    std::string symmetry_on_conflict = SymmetryOptions{}.on_conflict;
    std::string symmetry_plane = "detect";            // detect | X | Y | Z (object axes; the plane axis = offset)
    double symmetry_offset = SymmetryOptions{}.plane_offset;
    std::string scale_mode = ScaleOptions{}.mode;
    double reference_length = ScaleOptions{}.reference_length;
    std::string required_attribute = PreferenceOptions{}.required_attribute;
    std::string forbidden_attribute = PreferenceOptions{}.forbidden_attribute;
    std::string allow_attribute = PreferenceOptions{}.allow_attribute;
    std::string avoid_attribute = PreferenceOptions{}.avoid_attribute;
    std::int64_t max_candidates_per_chart = EnhancedOptions{}.max_candidates_per_chart;
    std::int64_t max_solves_per_chart = EnhancedOptions{}.max_solves_per_chart;
    bool prune_dominated = SearchOptions{}.prune_dominated;
    std::string last_messages = "";                   // host only: the native run_diagnostics lines of the last run (JSON)
};
struct HostOptions {
    bool protect_existing = true;
    std::string pack_after = "none";
    bool show_advanced = false;
    std::string last_status = "";
    std::string last_report = "";
};
// Optional recovery allowances are total phase seconds, never per-chart renewals.
// Zero disables the corresponding operation. Nonzero allowances become supported as
// the native recovery implementations are qualified; unsupported requests fail explicitly.
struct RepairOptions {
    bool require_open_protected_cuts=false; // locked seam endpoints must differ after storage
    double relax_faces_seconds=0;
    double region_recovery_seconds=0;
    double minimum_relative_area=0; // per isolated face, fraction of original absolute UV area
    double singleton_seconds=0;
    double join_seconds=0;
    std::string join_mode="small"; // small: <=4 faces/1 s; fewer: unrestricted faces/2 s
    double host_proposal_seconds=0;
    std::string storage_precision="float64";
};
struct PipelineOptions {
    std::string profile = "custom";
    std::optional<IslandBand> target_islands = std::nullopt;
    std::optional<double> deadline_s = std::nullopt;
    bool refine_guard = true;
    bool generate = true;
    bool refine_paths = true;
    std::int64_t refine_passes = 3;
    bool refine_cycles = true;
    bool clean_contours = true;
    bool uniform_density = true;
    bool orient_islands = true;
    LayoutOptions layout = {};
    PathOptions paths = {};
    BoundaryOptions boundary = {};
    double p95_growth = 1.05;
    double max_growth = 1.1;
    double abs_anisotropy_p95 = 2.0;
    double abs_area_p95 = 2.5;
    bool accept_partial = true;
    bool gate = true;
    std::int64_t gate_max_faces = 1500000;
    std::int64_t gate_max_components = 200;
    double gate_min_main_fraction = 0.3;
    double gate_max_degenerate = 0.2;
    double passthrough_fraction = 0.001;
    std::int64_t repair_exempt_faces = 3;
    double repair_exempt_fraction = 0.01;
    double fallback_density = 0.01;
    bool final_gate = true;
    double final_gate_tolerance = 1.05;
    std::int64_t version = 1;
    // The product default preserves user pins and applies the native corrections;
    // legacy_010 reproduces the frozen 0.1.0 source and is selected explicitly by differential tests.
    std::string policy = "native_baseline";
    // The enhanced_v1 proposal policy (T28+); inert under legacy_010 and native_baseline.
    EnhancedOptions enhanced = {};
    SolverOptions solver = {};
    PackOptions packing = {};
    HostOptions host = {};
    RepairOptions repair = {};
};
struct ViewContext {
    Vec3 direction_world{0,0,1};
    std::array<double,9> object_to_world{1,0,0,0,1,0,0,0,1}; // row-major linear part
};
CLEANSEAMS_CORE_EXPORT PipelineOptions pipeline_options_from_json(const nlohmann::json&);
CLEANSEAMS_CORE_EXPORT nlohmann::json pipeline_options_json(const PipelineOptions&);
CLEANSEAMS_CORE_EXPORT UiSettings ui_settings_from_json(const nlohmann::json&);
CLEANSEAMS_CORE_EXPORT nlohmann::json ui_settings_json(const UiSettings&);
CLEANSEAMS_CORE_EXPORT nlohmann::json option_registry();
// The run options of the UI settings (source ops.build_options); packing settings keep
// their defaults and are never validated for a run.
CLEANSEAMS_CORE_EXPORT PipelineOptions resolve_ui_options(const UiSettings&, const std::optional<ViewContext>& = std::nullopt);
// The packing options of the UI settings: packing and the
// host's packer choice only, tiles parsed as Python int(); run settings keep their defaults.
CLEANSEAMS_CORE_EXPORT PipelineOptions resolve_ui_pack_options(const UiSettings&);
CLEANSEAMS_CORE_EXPORT UiSettings apply_ui_profile(UiSettings, const std::string&);
CLEANSEAMS_CORE_EXPORT PipelineOptions apply_classified_profile(PipelineOptions, const std::string&);
}
