// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
#include "cleanseams/options.hpp"
#include "cleanseams/layout.hpp"
#include <atomic>
#include <functional>
#include <set>
namespace cleanseams {
struct ExposureField; // cleanseams/fields.hpp (the measured all-views exposure)
class TriangleBvh; // cleanseams/geometry_bvh.hpp (T36/T37 detectors)
struct StrictMirrorMap; // cleanseams/symmetric_moves.hpp (T38 strict symmetry)
class PhysicalScale; // cleanseams/physical_scale.hpp (T39 physical-distance policy)
// Raw host attributes. Interpretation/clamping and proposal execution belong to Phase B.
struct ConstraintInput {
    std::vector<std::uint8_t> required_edges;
    std::vector<std::uint8_t> forbidden_edges;
    std::vector<double> allow_points;
    std::vector<double> avoid_points;
    // Painted host attributes read by name (T29); interpreted by resolve_preferences.
    std::vector<PaintedAttribute> attributes;
};
struct ProposalIdentity {
    std::uint64_t id=0;
    std::uint64_t parent_revision=0;
    std::string generator;
    std::string generator_version;
    std::vector<std::string> affected_chart_keys;
    std::vector<std::int64_t> added_edges;
    std::vector<std::int64_t> removed_edges;
};
// ============================================================================
// Phase B common proposal and acceptance model (Section 11, T28). Everything below is
// used by policy enhanced_v1 only; legacy_010 and native_baseline never reach it.
//
//   propose -> cheap constraint checks -> complete trial state -> exact topology
//   -> solve changed charts -> measure (all validity checks + strict overlap, pins,
//   budgets, run-origin envelope) -> common ordering -> atomic commit, or discard.
//
// A proposal is an edge delta against one committed RunState revision. There is no
// per-operation validation switch: every generator's proposal goes through the same
// evaluate_proposal, and commit_proposal refuses anything that evaluation did not admit
// at the trial's current revision. The unchanged starting state is always candidate zero.
// ============================================================================
constexpr const char* enhanced_policy_version="enhanced_v1.policy.1";
constexpr const char* quality_envelope_version="run_origin_envelope.v1";
// GB (Section 11.3 policy): the validated repair of candidate zero's strict-only overlaps and the
// explicit report of the defects it cannot repair (src/proposals/strict_repair.cpp).
constexpr const char* strict_repair_version="strict_overlap_repair.v1";
// Resolved hard constraints of a run (Section 11.1). Edge masks have one flag per mesh edge.
struct HardConstraints {
    EdgeMask locked;                   // must stay cut: legacy locks, protect_existing, generated protection
    EdgeMask required;                 // must be cut (T29 painted cleanseams_required)
    EdgeMask forbidden;                // must not be an internal seam (T29 cleanseams_forbidden)
    std::vector<MeshIndex> anchors;    // manual anchors (sorted): never moved by an endpoint change
    std::vector<MeshIndex> seam_anchors; // anchors on the run origin's internal seams: never abandoned
    PinMap pins;                       // explicit user pins {original corner: uv}
    Revision revision=0;               // changes whenever painted constraints change (T29)
};
struct EndpointChange {MeshIndex from=0,to=0;};
struct Proposal {
    ProposalIdentity identity;         // id, parent revision, generator/version, affected chart keys, edge delta
    // Declared seam endpoint moves (sliding, T30/T31) and the endpoint set they must stay in.
    std::vector<EndpointChange> endpoints;
    std::optional<std::vector<MeshIndex>> allowed_endpoints;
    // Solver initialization identity: "" is the backend's default solve. A nonempty
    // identity also re-solves the listed charts of the parent that the delta leaves unchanged.
    std::string initialization;
    std::vector<std::string> resolve_chart_keys;
    nlohmann::json predicted=nlohmann::json::object();  // inexpensive predictions: ordering only, never a rejection
    nlohmann::json provenance=nlohmann::json::object();
};
// Measurements of one chart of a candidate state.
struct CLEANSEAMS_CORE_EXPORT ChartOutcome {
    std::string key;                   // ChartKey digest
    std::uint64_t faces=0;
    double area=0;                     // 3D area of the chart's faces
    bool exempt_only=false;            // every face exempt: not measured, never a budget or validity failure
    bool disk=false,measured=false,valid=false,overlap_free=false,pins_kept=true,within_budget=false,tested=true;
    // GB N2: every pin of the chart passes native_baseline's pin rule (numpy.allclose, rtol 1e-5,
    // atol 1e-8: the pipeline's pins_kept); `pins_kept` is the strict rule (pin_tolerance, or the
    // exact float32 rounding of the pin on a float32-storing host).
    bool pins_native_kept=true;
    // GB: a strict-only defect of candidate zero that this chart inherits unchanged: every
    // compatibility check and the overlap test pass, and only a strict check fails: positive-area
    // intrachart overlaps (strict overlap check) and/or (GB N2) a pinned UV that native_baseline's
    // pin rule keeps and the strict pin rule does not (a host solve that moved the pin by a few
    // float32 steps). It is reported and never silently compliant, but it does not make the state
    // inadmissible: the flag is set only on candidate zero's own charts and survives only while a
    // trial leaves the chart exactly as it was; a chart a trial solves or moves is measured afresh
    // with the strict rules.
    bool strict_defect=false;
    std::uint64_t flipped=0,nonfinite=0,collapsed=0,crossings=0,overlaps=0;
    double anisotropy_p95=quality_nan,anisotropy_p99=quality_nan,anisotropy_max=quality_nan,area_p95=quality_nan;
    double budget_ratio=0;             // worst of the budgeted quantities over its limit
    std::vector<MeshIndex> ancestors;  // charts of the run origin holding this chart's faces (face ancestry)
    nlohmann::json to_json() const;
};
// Section 11.2 terms of the seam objective of one state.
struct CLEANSEAMS_CORE_EXPORT SeamObjective {
    std::string version;
    double area=0,seam_length=0,exposed_length=0,turn_length=0,fragment_area=0;
    std::uint64_t seam_edges=0,charts=0,components=0,additional_charts=0,small_fragments=0;
    bool exposure_measured=false;
    double L=0,V=0,T=0,F=0,J=0;
    nlohmann::json to_json() const;
};
enum class BudgetStatus {compliant,improved_over_budget,over_budget};
CLEANSEAMS_CORE_EXPORT std::string_view budget_status_name(BudgetStatus);
struct BudgetSummary {
    bool measured=false;               // GB N7: set by the measurement (a trial refused before it has none)
    bool compliant=true;
    double worst_ratio=0;              // max over measured charts of ChartOutcome::budget_ratio
    double total_excess=0;             // area-weighted sum of max(0, ratio - 1)
    std::uint64_t over_budget_charts=0;
};
// Face-domain comparison against the run origin (Section 11.1): the faces of every chart
// whose key the origin does not have, measured on the candidate and on the origin.
struct CLEANSEAMS_CORE_EXPORT EnvelopeCheck {
    // GB N7: `inside` is the admission decision; `checked` whether the origin and the trial were
    // actually compared (not for an empty domain, an unmeasured origin, an infeasible parent or a
    // trial refused before the envelope). The report prints `inside` null when nothing was compared.
    bool inside=true;
    bool checked=false;
    std::uint64_t faces=0;
    double area=0;
    double origin_anisotropy_p95=quality_nan,trial_anisotropy_p95=quality_nan,origin_area_p95=quality_nan,trial_area_p95=quality_nan;
    double origin_anisotropy_max=quality_nan,trial_anisotropy_max=quality_nan;
    std::string reason;
    std::vector<std::size_t> relief_rebases;   // GB m10: the staged relief rebases the comparison applied
    nlohmann::json to_json() const;
};
// Per-triangle distortion samples of a measured map, exactly as the distortion report weighs
// them (finite, non-sliver, non-exempt triangles; area relative to the chart's density).
struct DistortionSample {MeshIndex face=0;double weight=0,anisotropy=0,area=0;};
CLEANSEAMS_CORE_EXPORT std::vector<DistortionSample> distortion_samples(const QualitySnapshot&);
// The run-origin quality envelope: cumulative changes of an enhancement pass are compared
// with the state the pass started from, never with the previous step, so repeated small
// degradations cannot ratchet. Limits: P95 anisotropy and area of the affected face domain
// may grow by p95_growth, its maximum anisotropy by max_growth (PipelineOptions).
class CLEANSEAMS_CORE_EXPORT BaselineEnvelope {
public:
    BaselineEnvelope()=default;
    // `snapshot`: the origin measured on every chart with the run's exemptions; `uv` its map.
    BaselineEnvelope(Revision origin,const QualitySnapshot& snapshot,std::span<const Vec2> uv,double p95_growth,double max_growth);
    Revision origin_revision() const {return origin_;}
    const std::vector<std::string>& origin_chart_keys() const {return keys_;}
    const std::vector<MeshIndex>& origin_chart_of_face() const {return chart_of_face_;}
    bool has_key(const std::string& key) const;
    // Charts of a candidate that differ from the origin: a key the origin lacks, or any
    // corner whose UV the origin map does not have. Their faces form the compared domain.
    std::vector<MeshIndex> domain_charts(const CutSurface&,std::span<const Vec2> uv) const;
    // `faces`: the domain (sorted, unique original faces); `trial`: samples of the candidate.
    // `active`: the relief rebases (stage_rebase) that apply to this trial (GB m10).
    EnvelopeCheck check(std::span<const MeshIndex> faces,std::span<const DistortionSample> trial,std::span<const std::size_t> active={}) const;
    // T32: faces whose run origin is an earlier baseline state, the chart before the
    // ordinary gate's over-budget split (Section 21 schedules the enhanced operations before
    // that gate). Their origin samples become `samples` (of that state); chart keys and the
    // origin map that decide the compared domain stay candidate zero's. Unconditional: every
    // later trial is measured against it.
    void rebase_faces(std::span<const MeshIndex> faces,std::span<const DistortionSample> samples);
    // GB m10: the same rebase for relief only. Staged per recorded split (the returned index
    // is the split record's position), it applies to a trial only when that trial is the
    // relief proposal welding the split (hotspot_relief, kind relief_before_split, that
    // split_index) or its parent state adopted such a proposal; every other trial is measured
    // against candidate zero's split pieces.
    std::size_t stage_rebase(std::span<const MeshIndex> faces,std::span<const DistortionSample> samples);
    std::size_t staged_rebases() const {return staged_.size();}
    nlohmann::json to_json() const;
private:
    Revision origin_=0;
    std::uint64_t rebased_faces_=0;
    double p95_growth_=1.05,max_growth_=1.1;
    std::vector<std::string> keys_;          // sorted
    std::vector<MeshIndex> chart_of_face_;
    std::vector<DistortionSample> samples_;  // sorted by face
    std::vector<Vec2> uv_;                   // the origin map
    struct Staged {std::vector<MeshIndex> faces;std::vector<DistortionSample> samples;}; // faces sorted; samples of those faces, by face
    std::vector<Staged> staged_;
};
// Per-sector area mass: each triangle gives a third of its 3D area to each corner's sector.
// Split sectors keep separate masses and the masses always total the triangulated area.
CLEANSEAMS_CORE_EXPORT std::vector<double> sector_area_mass(const CutSurface&);
// Charts of `origin` holding faces of each chart of `surface` (sorted), from face identity.
CLEANSEAMS_CORE_EXPORT std::vector<std::vector<MeshIndex>> face_ancestry(const CutSurface& surface,std::span<const MeshIndex> origin_chart_of_face);
// One committed, authoritative state: effective cuts, their surface, the map solved on
// that surface, measurements and provenance (Section 6.5).
struct RunState {
    Revision revision=0;
    EdgeMask cuts;                              // effective cuts: always surface->tables().cut
    std::shared_ptr<const CutSurface> surface;
    std::vector<Vec2> uv;                       // dense, per original corner
    std::vector<ChartOutcome> charts;           // per chart of `surface`
    SeamObjective objective;
    BudgetSummary budget;
    bool feasible=false,valid=false;            // hard constraints; every validity check of every chart
    // GB (Section 11.3 policy): charts carrying a strict-only overlap of candidate zero that no
    // validated trial replaced yet (ChartOutcome::strict_defect). `valid` stays false while any
    // remain; `strict_only` marks a feasible state whose only failures are these charts: it
    // ranks as admissible, after every state with fewer of them (validity first).
    std::vector<MeshIndex> strict_defects;
    bool strict_only=false;
    bool admissible() const {return feasible && (valid || strict_only);}
    std::string invalid_reason;
    nlohmann::json violations=nlohmann::json::object();  // exact edge/anchor/chart IDs of a failure
    nlohmann::json provenance=nlohmann::json::array();   // commit records, origin first
};
// GB m10: the staged relief rebases (BaselineEnvelope::stage_rebase) that apply to a trial of
// `proposal` on `parent`: the split welded by that relief proposal and every split whose relief
// proposal the parent's provenance committed (sorted, unique).
CLEANSEAMS_CORE_EXPORT std::vector<std::size_t> relief_rebases(const RunState& parent,const Proposal& proposal);
// The chart solver of the transaction: a (nc) map holding every requested chart of the
// surface, honoring the pins. `initialization` is the proposal's identity.
using ProposalSolver=std::function<std::vector<Vec2>(const std::shared_ptr<const CutSurface>&,std::span<const MeshIndex> charts,
                                                     const PinMap& pins,const std::string& initialization)>;
// Solves and candidate states retained per affected origin chart (Section 21 bounds).
struct CLEANSEAMS_CORE_EXPORT SearchCounters {
    std::map<std::string,std::uint64_t> solves_by_chart,candidates_by_chart;
    std::map<std::string,std::map<std::string,std::uint64_t>> by_generator; // generator -> {tried, admissible, committed, <rejection>}
    std::uint64_t proposals=0,solves=0,charts_solved=0,duplicates=0,commits=0;
    nlohmann::json to_json() const;
};
// ------------------------------------------------------------ objective significance (GB, M3)
// seam_objective.v2: a J difference is a tie when it is within the sampling uncertainty of the
// raw exposure term. L, T and F are exact functions of the cuts; only V is sampled: the raw
// exposure of a vertex is the cosine-weighted visible fraction over K all-views directions,
//   x_v = sum_d w_vd b_vd / W_v,   W_v = sum_d w_vd,   b_vd in {0, 1},  w_vd = cos+ (float32),
// an edge's exposure the mean of its two vertices', and V = sum(length * exposure) / sqrt(A).
// Rotating the object moves the fixed direction set relative to it, so the directions are the
// sampling units. Between two cut states only the edges they do not share contribute:
//   dV = sum_v a_v x_v,   a_v = 1/2 sum_{e in A xor B, e ~ v} (+-length_e / sqrt(A))
// (+ for an edge of A only, - for one of B only; physical boundary edges are in both). The
// linearized (delta-method) variance of this sum of ratio estimators over shared direction
// samples, treated as K independent draws, is
//   Var(dV) ~= K/(K-1) sum_d G_d^2,   G_d = sum_v a_v w_vd (b_vd - x_v) / W_v,
// which keeps the correlation between vertices seen by the same direction. U = z w_exposure
// sigma(dV) with z = 2 (about 95 %, one-sided about 97.7 %) is the objective's uncertainty:
// |J_a - J_b| <= max(tie_band * max|J|, U) is a tie, and a tie keeps the incumbent. Budget and
// validity ordering are unchanged; identical cut sets have U = 0 (map writers compare as in v1).
// Measured calibration and limits: docs/ENHANCED_SEAMS.md, "GB — objective v2".
constexpr const char* exposure_uncertainty_version="exposure_uncertainty.v1";
constexpr double seam_objective_v2_significance=2.0;   // z of seam_objective.v2
class CLEANSEAMS_CORE_EXPORT ExposureUncertainty {
public:
    ExposureUncertainty()=default;
    // From the measured field (its per-vertex exposure and visibility) of `mesh`.
    ExposureUncertainty(std::shared_ptr<const MeshSnapshot> mesh,const ExposureField& field);
    std::uint64_t directions() const {return static_cast<std::uint64_t>(directions_.size());}
    bool usable() const {return mesh_ && directions_.size()>=2;}
    // sigma(V(a) - V(b)) for two effective cut masks; `area` is the objective's relevant area.
    double delta_v_sigma(std::span<const std::uint8_t> a,std::span<const std::uint8_t> b,double area) const;
    nlohmann::json to_json() const;
private:
    std::shared_ptr<const MeshSnapshot> mesh_;
    std::vector<Vec3> directions_;
    std::size_t words_=0;
    std::vector<std::uint64_t> visible_;   // per vertex, `words_` words of visibility bits
    std::vector<double> exposure_,weight_; // per vertex: x_v and W_v
};
struct ProposalContext;
// U of seam_objective.v2 between two cut masks: z * w_exposure * sigma(dV); 0 under
// seam_objective.v1, without a measured field or an uncertainty model.
CLEANSEAMS_CORE_EXPORT double objective_margin(const ProposalContext&,std::span<const std::uint8_t> a,std::span<const std::uint8_t> b,double area);
// Whether the context's ordering applies the v2 margin (seam_objective.v2 with a usable model).
CLEANSEAMS_CORE_EXPORT bool significance_active(const ProposalContext&);
// The record of one v2 comparison: delta J, sigma(dV), U and whether it is significant.
CLEANSEAMS_CORE_EXPORT nlohmann::json significance_record(const ProposalContext&,double challenger_J,double incumbent_J,double margin,bool ranks_first);
// GB (visual review V2): the practical-effect floor of seam_objective.v2 between a challenger's
// and an incumbent's cut masks: enhanced.objective.significance.practical_floor * |incumbent_J|
// wherever the v2 margin applies (significance_active) and the challenger adds at least one cut
// edge; 0 for a change that only removes seam edges (identical cuts included), under v1 and
// without a sampling model. The margin of an adoption decided by J is max(U, floor).
CLEANSEAMS_CORE_EXPORT double practical_floor(const ProposalContext&,std::span<const std::uint8_t> challenger,
    std::span<const std::uint8_t> incumbent,double incumbent_J);
struct CLEANSEAMS_CORE_EXPORT ObjectiveSignificance {
    double uncertainty=0;   // U: z * w_exposure * sigma(dV) (objective_margin)
    double floor=0;         // practical_floor
    double margin() const {return uncertainty>floor?uncertainty:floor;}
};
CLEANSEAMS_CORE_EXPORT ObjectiveSignificance objective_significance(const ProposalContext&,std::span<const std::uint8_t> challenger,
    std::span<const std::uint8_t> incumbent,double area,double incumbent_J);
// max(U, floor) of the pair (the margin compare_candidates receives for an adoption).
CLEANSEAMS_CORE_EXPORT double objective_margin(const ProposalContext&,std::span<const std::uint8_t> challenger,
    std::span<const std::uint8_t> incumbent,double area,double incumbent_J);
// The record with the floor: also `practical_floor`, `margin`, `within_margin` and
// `below_practical_floor` (beyond U, within the floor: the floor alone kept the incumbent).
CLEANSEAMS_CORE_EXPORT nlohmann::json significance_record(const ProposalContext&,double challenger_J,double incumbent_J,
    const ObjectiveSignificance&,bool ranks_first);
struct ResolvedPreferences;
class SeamNecessityCache;
// ------------------------------------------------------------ bounded alternative search (T40)
// Section 21: the enhancement pass is a bounded alternative search, not a claim of global
// optimality. Within a pass (run_enhancement installs a SearchLedger in the context) every
// round ranks its unevaluated proposals by a cheap exact estimate (the length and raw exposure
// part of J that the edge delta alone determines), deduplicates by effective cut state plus
// the numerical initialization contract across the whole pass, applies the per-chart bounds
// and each feature's share of them, and refuses a candidate before its solve only on an
// actual constraint/topology failure or a proven bound (DominanceBound). The incumbent is the
// best validated state; a trial never replaces it unless it ranks strictly before it, and a
// deadline or cancellation stops new work without erasing it. select_proposal without a
// ledger is the exhaustive T28 round (every representable proposal in id order).
constexpr const char* bounded_search_version="bounded_search.v1";
// A proven lower bound. J depends on the effective cuts only (never on the solved map), so
// the exact J of a trial is known once its topology is built. A trial cannot rank before a
// compliant, admissible incumbent unless it is compliant with a J strictly lower outside the
// tie band; it cannot become the round's leader (a compliant, admissible evaluated candidate)
// unless its J is strictly lower, or tied with a lower candidate id.
struct CLEANSEAMS_CORE_EXPORT DominanceBound {
    std::optional<double> incumbent_J;
    std::optional<double> leader_J;
    std::uint64_t leader_id=0;
    double tie_band=1e-9;
    // GB (seam_objective.v2): the incumbent test also needs the trial's J to be lower by more than
    // U (objective_margin of the trial's and the incumbent's cuts, computed by the transaction).
    bool significance=false;
    bool dominates(double J,std::uint64_t id,double incumbent_margin=0.0) const;
    nlohmann::json to_json() const;
};
// GB (m6): the pass-wide deduplication keys. A state is keyed by the BLAKE2b-256 digest of its
// effective cuts plus its numerical initialization contract (initialization identity, charts
// re-solved on request, parent revision), never by the cut bytes themselves; each entry keeps
// its cut delta against the ledger's reference cuts (the first incumbent), and two states are
// the same only when the digest key and that delta are equal (full comparison on a collision).
class CLEANSEAMS_CORE_EXPORT StateKeys {
public:
    struct Found {bool duplicate=false;const nlohmann::json* first=nullptr;};
    // The record of an equal state, if any (does not insert).
    Found find(const std::string& key,std::span<const std::uint8_t> cuts) const;
    // Inserts the state unless an equal one is present (then keeps the first record).
    bool insert(const std::string& key,std::span<const std::uint8_t> cuts,nlohmann::json record);
    std::size_t size() const {return size_;}
    std::uint64_t collisions() const {return collisions_;}
private:
    struct Entry {std::vector<MeshIndex> delta;nlohmann::json record;};
    std::vector<MeshIndex> delta_of(std::span<const std::uint8_t> cuts) const;
    std::vector<std::uint8_t> reference_;
    bool has_reference_=false;
    std::map<std::string,std::vector<Entry>> entries_;
    std::size_t size_=0;
    std::uint64_t collisions_=0;
};
// The digest key of a trial state (m6): cut_digest(cuts), the initialization identity, the sorted
// charts re-solved on request and the parent revision.
CLEANSEAMS_CORE_EXPORT std::string state_digest_key(std::span<const std::uint8_t> cuts,const std::string& initialization,
    std::vector<std::string> resolve_chart_keys,std::optional<Revision> parent);
// The pass-wide record of the bounded search: deduplication keys with the first candidate
// that carried them, each feature's use of every origin chart's bounds, where the pass
// stopped, and every incumbent the pass kept (candidate zero first).
struct CLEANSEAMS_CORE_EXPORT SearchLedger {
    StateKeys states;                                             // digest key -> {id, generator, round, rejection, J}
    std::map<std::string,std::map<std::string,std::uint64_t>> feature_solves,feature_candidates; // feature -> origin chart -> count
    std::map<std::string,std::map<std::string,double>> feature_seconds; // feature -> {generate, evaluate}
    std::uint64_t rounds=0,ranked=0,reordered=0,dominated=0,duplicates=0,bound_refusals=0,stopped_proposals=0;
    std::string stopped;                                          // "", "deadline" or "cancelled"
    nlohmann::json incumbents=nlohmann::json::array();            // {round, generator, revision, J, compliant, worst_ratio}
    nlohmann::json to_json() const;
};
struct ProposalContext {
    std::shared_ptr<const MeshSnapshot> mesh;
    PipelineOptions options;                    // layout budgets, growth limits and the enhanced record
    HardConstraints constraints;
    std::vector<double> raw_exposure;           // measured per edge; empty: exposure unmeasured
    // The run's resolved paint and explicit viewer (T29), read by generators' routing heuristics
    // and by those that price seams (T31, T32, T34+); null: none painted (context_preferences
    // resolves a neutral record from the raw exposure).
    std::shared_ptr<const ResolvedPreferences> preferences;
    // Generator diagnostics (detector evidence, refused partitions and why), by generator
    // name; reported as report.enhanced.diagnostics.
    nlohmann::json diagnostics=nlohmann::json::object();
    std::vector<std::uint8_t> exempt;           // the run's exemption policy (faces); proposals never change it
    ProposalSolver solve;
    BaselineEnvelope envelope;                  // run origin (set by begin_enhancement)
    SearchCounters counters;
    std::shared_ptr<std::atomic<Revision>> revisions=std::make_shared<std::atomic<Revision>>(0);
    CancellationToken cancel;
    Deadline deadline;                          // stops new work; the incumbent is never erased
    // T32: the ordinary gate's over-budget splits of candidate zero (relief before the split);
    // relief generators route with `preferences` above.
    std::vector<SplitRecord> baseline_splits;
    // T33: seam-end necessity evidence of the pass, by evidence key (null: generators keep their own).
    std::shared_ptr<SeamNecessityCache> necessity;
    // T36/T37: the triangle BVH of `mesh`'s geometry revision, built on first use by
    // context_bvh and shared by the pass's detectors (geometry never changes in a pass).
    std::shared_ptr<const TriangleBvh> bvh;
    // T38: the pass's mirror map (context_symmetry builds it on first use unless the symmetry
    // mode is asymmetric); in strict mode the transaction pairs every proposal with it.
    std::shared_ptr<const StrictMirrorMap> symmetry;
    // T39: the pass's physical scale, frozen by begin_enhancement on candidate zero (the
    // operation start) and never recomputed as cuts change (Section 20).
    std::shared_ptr<const PhysicalScale> scale;
    // T40: the bounded search of the pass (run_enhancement installs one; null: exhaustive rounds).
    std::shared_ptr<SearchLedger> search;
    // GB: the run's backend stores maps in float32 (UnwrapBackend::float32_storage, the Blender
    // host): the strict pin check also keeps a pinned coordinate that is exactly the float32
    // rounding of the pin (check_stored_uv's rule), nothing looser. Core maps keep pin_tolerance.
    // GB N2: this covers storage only. A host solve that moves a pin by a few float32 steps (Blender's
    // CONFORMAL unwrap does) fails this rule, as it fails the host's exact apply check
    // (Session::result_for_apply); on candidate zero that is a strict-only pin defect (reported,
    // repaired when a validated re-map can, never a refusal native_baseline would not make).
    bool float32_storage=false;
    // GB (M3): the sampling uncertainty of `raw_exposure` (the pipeline builds it from the same
    // measured field); null: none, and seam_objective.v2 then compares like v1.
    std::shared_ptr<const ExposureUncertainty> exposure_uncertainty;
    // GB (m3): the solver initialization identities `solve` honours. Null: every identity (a
    // caller's own solver). The pipeline's backend solver honours only "" (its default solve), so
    // a proposal naming another identity is refused before any work (unsupported_initialization)
    // instead of being solved as the default under another numerical contract.
    std::function<bool(const std::string&)> supports_initialization;
    Revision next_revision() {return ++*revisions;}
};
// GB (m3): the origin charts a proposal affects, for the per-chart bounds: the charts of the
// faces its edge delta touches and the charts it asks to re-solve (resolve_chart_keys, keys of
// `parent`'s charts). Out-of-range edges and unknown keys are skipped (the transaction names them).
CLEANSEAMS_CORE_EXPORT std::set<std::string> proposal_origin_charts(const Proposal&,const RunState& parent,const ProposalContext&);
// GB (m4): whether a failure of one generator or trial is isolated (recorded, the incumbent
// kept, the pass goes on). Cancellation stops the pass, and contract errors of the caller's
// inputs (options, versions) propagate; every other error of that stage or trial is isolated.
CLEANSEAMS_CORE_EXPORT bool isolated_failure(const std::exception&);
// The report record of an isolated failure (code, field, message).
CLEANSEAMS_CORE_EXPORT nlohmann::json failure_record(const std::exception&);
// The trial of one proposal. Every writer (the proposal itself first, then any forced
// post-generation cleanup) gives the trial a new revision; the evaluation belongs to one
// revision, and a trial is committable only when both match.
struct TrialState {
    Revision revision=0,evaluated_revision=0;
    EdgeMask cuts;
    std::shared_ptr<const CutSurface> surface;
    std::vector<Vec2> uv;
    bool uv_written=false;                      // a writer replaced the map: measured as written, not re-solved
    std::vector<std::string> writers;
};
struct CLEANSEAMS_CORE_EXPORT EvaluatedProposal {
    Proposal proposal;
    Revision parent_revision=0;
    TrialState trial;
    bool feasible=false,topology_ok=false,valid=false,inside_envelope=false;
    // GB: inherited, unchanged strict-only defect charts of candidate zero (the trial is
    // admissible with them; `valid` is false while any remain).
    std::vector<MeshIndex> strict_defects;
    std::string rejection;                      // "" when admissible
    nlohmann::json detail=nlohmann::json::object();
    std::vector<MeshIndex> changed_charts;      // trial charts solved and measured
    std::vector<ChartOutcome> charts;           // per trial chart
    SeamObjective objective;
    BudgetSummary budget;
    EnvelopeCheck envelope;
    std::uint64_t solves=0;
    double seconds=0;
    Revision constraints_revision=0;            // HardConstraints::revision the evaluation used (T29)
    bool admissible() const {return rejection.empty() && trial.evaluated_revision==trial.revision && trial.revision!=0;}
    nlohmann::json summary() const;
};
// Resolves the run's hard constraints: locks, anchors and pins as given; required and
// forbidden masks (T29 painted attributes) default to empty.
CLEANSEAMS_CORE_EXPORT HardConstraints hard_constraints(const MeshSnapshot&,const std::optional<EdgeMask>& locked,std::span<const TransferIndex> anchors,
    const PinMap& pins,const std::optional<EdgeMask>& required=std::nullopt,const std::optional<EdgeMask>& forbidden=std::nullopt);
// Objective terms of a cut state (raw exposure, never a routing cost).
CLEANSEAMS_CORE_EXPORT SeamObjective seam_objective(const CutSurface&,std::span<const double> raw_exposure,std::span<const std::uint8_t> exempt,
    const LayoutOptions& layout,const ObjectiveOptions& objective);
// Measures a complete state: feasibility, disks, validity + strict overlap, pins, budgets,
// objective. `ancestry` maps faces to the origin's charts (empty: none).
CLEANSEAMS_CORE_EXPORT RunState measure_state(ProposalContext&,std::shared_ptr<const CutSurface>,std::vector<Vec2> uv,
    Revision revision,nlohmann::json provenance);
// BLAKE2b-256 hex digest of an effective cut mask (report identity of a cut state).
CLEANSEAMS_CORE_EXPORT std::string cut_digest(std::span<const std::uint8_t> cuts);
// Report record of a state: revision, feasibility/validity, charts, cut digest, objective, budget.
CLEANSEAMS_CORE_EXPORT nlohmann::json run_state_summary(const RunState&);
// Starts an enhancement pass from a complete candidate zero: measures it and fixes the
// run-origin envelope on it.
CLEANSEAMS_CORE_EXPORT RunState begin_enhancement(ProposalContext&,std::shared_ptr<const CutSurface>,std::vector<Vec2> uv,
    nlohmann::json provenance=nlohmann::json::object());
// Fixes the run-origin envelope on `state` (candidate zero, or the first feasible state when
// candidate zero violates hard constraints).
CLEANSEAMS_CORE_EXPORT void anchor_envelope(ProposalContext&,const RunState&);
// Every hard-constraint violation of an effective cut set (exact IDs; empty when feasible).
CLEANSEAMS_CORE_EXPORT nlohmann::json constraint_violations(const MeshSnapshot&,const HardConstraints&,std::span<const std::uint8_t> cuts);
CLEANSEAMS_CORE_EXPORT EvaluatedProposal evaluate_proposal(const RunState& parent,Proposal,ProposalContext&);
// T40: the same evaluation with a proven bound checked after the exact topology and before
// the solve (`dominated`: no solve, the trial's exact objective and the bound in the detail).
CLEANSEAMS_CORE_EXPORT EvaluatedProposal evaluate_proposal(const RunState& parent,Proposal,ProposalContext&,const DominanceBound* bound);
// Forced post-generation writers: replace the trial's cuts (its changed charts are solved
// again) or map (measured as written). Either gives the trial a new revision and drops its
// evaluation; reevaluate_proposal must run before it can be committed.
CLEANSEAMS_CORE_EXPORT void write_trial_cuts(EvaluatedProposal&,EdgeMask cuts,const std::string& writer,ProposalContext&);
CLEANSEAMS_CORE_EXPORT void write_trial_uv(EvaluatedProposal&,std::vector<Vec2> uv,const std::string& writer,ProposalContext&);
CLEANSEAMS_CORE_EXPORT void reevaluate_proposal(const RunState& parent,EvaluatedProposal&,ProposalContext&);
// Atomic commit: the new state carries exactly the trial's evaluated cuts, surface and map.
// Throws ContractError(stale_state) for another parent revision or a stale evaluation, and
// ContractError(constraint_conflict) for a proposal the evaluation did not admit.
CLEANSEAMS_CORE_EXPORT RunState commit_proposal(const RunState& parent,const EvaluatedProposal&,ProposalContext&);
// A post-generation map writer on a committed state (output hygiene, T28): the written map
// is measured on the same cuts through the transaction and committed only if admissible;
// otherwise the state is returned unchanged. `record` receives the trial summary. T40:
// `no_worse` also refuses an admissible map that the common ordering ranks after `state`
// (`ranks_after_state`), e.g. a contour cleanup that loses budget compliance.
CLEANSEAMS_CORE_EXPORT RunState rewrite_state_uv(const RunState& state,std::vector<Vec2> uv,const std::string& writer,ProposalContext&,
    nlohmann::json* record=nullptr,bool no_worse=false);
// ------------------------------------------------------------ common final ordering (search.cpp)
// Lexicographic: feasible and valid first, then budget compliance; among compliant states
// the objective J, among over-budget states the worst budget ratio, the total excess, then
// J. Scores within the relative tie band are equal. Negative: `a` ranks before `b`.
// GB (strict-only defects): right after admissibility, the state with fewer unrepaired strict
// defects of candidate zero ranks first, but (GB N4, Section 11.1) only when it also keeps the
// other state's budget (keeps_budget); otherwise the defect count gives no rank and the budget
// ordering decides.
struct CandidateView {
    std::uint64_t id=0;                // candidate zero is 0; proposals use their stable IDs
    bool admissible=false;
    const SeamObjective* objective=nullptr;
    const BudgetSummary* budget=nullptr;
    // GB: unrepaired strict-only defect charts, compared right after admissibility (validity
    // first): fewer rank first. Zero for every state without such a defect.
    std::uint64_t strict_defects=0;
};
// GB (seam_objective.v2): `j_margin` is U of the pair (objective_margin of their cuts): the J
// term is a tie within max(tie band, U). 0 (the default) is the v1 comparison.
CLEANSEAMS_CORE_EXPORT int compare_candidates(const CandidateView& a,const CandidateView& b,const ObjectiveOptions&,double j_margin=0.0);
// GB N4 (Section 11.1): `trial` keeps `incumbent`'s budget: a compliant incumbent stays compliant;
// an over-budget one gets no higher worst ratio and no higher total excess (the tie band).
CLEANSEAMS_CORE_EXPORT bool keeps_budget(const BudgetSummary& trial,const BudgetSummary& incumbent,double tie_band);
// Budget label of a candidate against the incumbent it would replace.
CLEANSEAMS_CORE_EXPORT BudgetStatus budget_status(const BudgetSummary& candidate,const BudgetSummary& incumbent,const ObjectiveOptions&);
struct ProposalGenerator {
    std::string name,version;
    std::function<std::vector<Proposal>(const RunState&,ProposalContext&)> generate;
};
// T40: the cheap exact estimate that orders a round's queue: the change of w_length * L +
// w_exposure * V that the proposal's edge delta alone determines (internal edges only;
// `area` is the relevant 3D area of the pass). It orders; it never refuses.
CLEANSEAMS_CORE_EXPORT double proposal_estimate(const MeshSnapshot&,const Proposal&,std::span<const double> raw_exposure,double area,
    const ObjectiveOptions&);
struct SearchRound {
    RunState state;                            // the incumbent after the round (unchanged when nothing wins)
    bool committed=false;
    std::string stopped;                       // T40: "deadline" or "cancelled" when the round stopped early
    nlohmann::json candidates=nlohmann::json::array(); // candidate zero first, then every proposal
    nlohmann::json winner;                     // id and budget status of the retained state
};
// One round: every proposal is evaluated against `incumbent` (deduplicated by effective cuts
// plus initialization, bounded per chart), and the best admissible candidate replaces the
// incumbent only when it ranks strictly before it. The incumbent is always a candidate.
CLEANSEAMS_CORE_EXPORT SearchRound select_proposal(const RunState& incumbent,std::vector<Proposal> proposals,ProposalContext&);
// The enhancement pass over generators in schedule order (Section 21); each generator sees
// the current incumbent. Returns the retained state and the pass report.
struct EnhancementResult {RunState state;RunState origin;nlohmann::json report;};
CLEANSEAMS_CORE_EXPORT EnhancementResult run_enhancement(const RunState& origin,std::span<const ProposalGenerator> generators,ProposalContext&);
// The built-in generators of enhanced_v1 in schedule order (none before T30).
CLEANSEAMS_CORE_EXPORT std::vector<ProposalGenerator> enhanced_generators(const PipelineOptions&);
// T40: the complete Section 21 schedule of a run: every step in order with its feature name,
// whether it runs (and why not), the generators it contributes and its explicit bounds. The
// baseline steps (generation, refinement, contour cleanup, repair and the gate) run unchanged
// in candidate zero; `enhanced_generators` returns the enabled proposal steps in this order.
CLEANSEAMS_CORE_EXPORT nlohmann::json enhanced_schedule(const PipelineOptions&,bool painted_hard_constraints);
// ------------------------------------------------------------ strict-overlap repair (GB, strict_repair.cpp)
// Section 11.3 policy. The strict overlap check stays a validity check of every new candidate:
// a trial's solved or moved charts must pass it, so an invalid map never wins. Candidate zero's
// own strict-only defects (charts whose compatibility checks, pins and test pass but which have
// positive-area intrachart overlaps) do not refuse the run that native_baseline accepts: the
// pass first tries the validated repair below on exactly those charts; a defect it cannot repair
// stays in the state (ChartOutcome::strict_defect), is reported (strict_overlap_report) and only
// a validated trial may ever replace that chart.
//
// The repair ladder per chart (strict_overlap_repair.v1), each attempt one trial through the
// transaction (every validity check, strict overlap, pins, the run-origin envelope, measured
// budgets), in order of increasing change, stopping at the first attempt that is admissible,
// keeps the budget (a compliant incumbent stays compliant; an over-budget one gets no higher worst
// ratio or total excess) and ranks strictly before the incumbent:
//   untangle:1..3   the repair ladder's local untangle (mean-value harmonic re-solve of the rings
//                   around the overlapping faces, the rest of the chart fixed)
//   arap:20, arap:50  native ARAP re-solve of the chart from its LSCM map with more iterations
//   lscm            the chart's LSCM map
//   harmonic_arap:20  ARAP from the convex-boundary harmonic map (unpinned single-loop disks)
//   harmonic        the convex-boundary harmonic map itself (injective, usually over the envelope)
// then, as the last resort, the minimal splits (cut proposals solved by the run's solver; the best
// admissible one under the common ordering). The two sides of the overlapping face pairs are
// separated over the chart's dual (its own cuts are barriers, forbidden edges never cut):
// `min_cut` by the cheapest separating edge set (maximum flow; an edge costs what it adds to J's
// length and exposure terms), `separate` by the border of two regions grown from the sides;
// `bisect` grows two regions from the chart's two farthest faces (the ordinary split seeding).
// Existing cuts are kept and every non-disk piece is opened along admissible cheapest paths.
// Bounds: every attempt is one solve on the chart's origin chart, within max_solves_per_chart and
// the feature share (search.feature_solves_per_chart); a deadline or cancellation stops the stage.
// Runs as soon as the incumbent is feasible and fails only the strict check, before any other
// proposal stage (after the constraint repair when painted constraints need it).
// GB N2: a chart whose only strict defect is a moved pin (no overlapping pair) gets the re-maps
// only (every native re-map holds the pins exactly); the untangle rungs and the splits, which
// address overlapping faces, are not applicable to it.
CLEANSEAMS_CORE_EXPORT SearchRound repair_strict_overlaps(const RunState& incumbent,ProposalContext&);
// report.enhanced.strict_overlap: candidate zero's strict-only defects (charts, overlapping pairs,
// overlap area), each chart's status (repaired or unrepaired) and attempts, the overall status
// (none, repaired, unrepaired) and the defects the retained state still carries. `stage` is the
// repair stage's round record (null when it did not run).
CLEANSEAMS_CORE_EXPORT nlohmann::json strict_overlap_report(const ProposalContext&,const RunState& origin,const RunState& retained,
    const nlohmann::json& stage);
// The strict pin rule of one pinned coordinate: kept within tolerance * (1 + |pin|) or, on a
// float32-storing host, when it is exactly the float32 rounding of the pin (storage only,
// nothing looser; GB N2: a pin a host solve moved by float32 steps is not kept).
CLEANSEAMS_CORE_EXPORT bool strict_pin_kept(double got,double pin,double tolerance,bool float32_storage);
// GB N2: report.enhanced.pin_check: the strict pin rule and native_baseline's, every pin of
// candidate zero that the strict rule does not keep (corner, chart, pin, UV, deviation and its
// float32 steps), the charts carrying such a strict-only pin defect with their status (repaired
// or unrepaired) and repair attempts, and the pins the retained state still moves.
CLEANSEAMS_CORE_EXPORT nlohmann::json pin_check_report(const ProposalContext&,const RunState& origin,const RunState& retained);
// GB N3 (src/engine/pipeline.cpp): the result's chart_over_budget rows of an adopted state `kept` by
// native_baseline's own listing rule. `native_rows` are native_baseline's rows of candidate zero
// (`origin`) and `source` the rule that produced them: "chart_gate" (P95 anisotropy over
// max_anisotropy_p95 * final_gate_tolerance on charts of at least split_min_faces faces),
// "repair_absorption" (the same P95 rule on every chart) or "not_judged" (no rows); both rules add
// FIX-04's area P95 rows (same face floor) after no anisotropy listing. A kept chart with candidate
// zero's key and map keeps native_baseline's row (or absence); every other chart is judged by the
// rule. Rows [chart, P95 to 3 places] in chart order; `record` receives the rule and the charts
// inherited and measured. The enhanced budget (P99/maximum caps) is not part of this rule.
CLEANSEAMS_CORE_EXPORT nlohmann::json native_rule_over_budget_rows(const RunState& kept,const RunState& origin,const nlohmann::json& native_rows,
    const std::string& source,const PipelineOptions&,nlohmann::json* record=nullptr);
// A written map of `parent`'s cuts evaluated through the transaction (a writer trial: the charts it
// moves are measured afresh), not committed. `identity` names the generator and initialization.
CLEANSEAMS_CORE_EXPORT EvaluatedProposal evaluate_written_map(const RunState& parent,Proposal identity,std::vector<Vec2> uv,const std::string& writer,
    ProposalContext&);
// ------------------------------------------------------------ painted preferences (T29)
// Section 12: four separate fields. Nothing here overwrites the measured exposure; soft
// preferences only change the enhanced routing heuristic; hard admissibility is enforced by
// the transaction (constraint_violations) for every candidate, candidate zero included.
struct CLEANSEAMS_CORE_EXPORT ResolvedPreferences {
    std::vector<double> raw_exposure;          // RawExposure per edge, measured (empty: not measured)
    std::vector<double> feature_desirability;  // FeatureDesirability per edge in [0, 1]
    std::vector<double> allow,avoid;           // SoftPreference per vertex in [0, 1] (finite weights clamped)
    std::vector<double> allow_edge,avoid_edge; // SoftPreference per edge: mean of its two vertices
    EdgeMask required,forbidden,locked;        // HardAdmissibility (forbidden flags on physical boundaries removed)
    std::vector<MeshIndex> anchors;
    std::optional<Vec3> view_direction;        // explicit viewer only; there is no implicit front
    std::vector<double> facing;                // per edge, only with an explicit viewer
    nlohmann::json sources=nlohmann::json::array();       // every input: domain/type status, counts, clamping
    nlohmann::json inapplicable=nlohmann::json::object(); // flags on physical boundaries (topology facts)
    nlohmann::json conflicts=nlohmann::json::object();    // hard conflicts with exact edge IDs: the run is refused
    std::string digest;                        // painted inputs, aliases' results, weights, view: proposal invalidation
    bool hard_constraints() const;             // any required or forbidden edge
    nlohmann::json summary() const;
    // Why an edge is attractive, avoided or forbidden: every term of the routing cost.
    nlohmann::json explain(const MeshSnapshot&,MeshIndex edge,const PreferenceOptions&) const;
};
// Resolves typed constraint arrays and painted host attributes (matched only by the names in
// `options`) natively: wrong domains/types are reported and ignored, NaN or infinite weights
// are invalid input, finite out-of-range weights are clamped and counted, a missing
// attribute is zero preference. `view` is the explicit viewer (layout.hide_dir) or none.
CLEANSEAMS_CORE_EXPORT ResolvedPreferences resolve_preferences(const MeshSnapshot&,const ConstraintInput&,const PreferenceOptions& options,
    const std::optional<EdgeMask>& locked,std::span<const TransferIndex> anchors,double feature_angle_deg,std::optional<Vec3> view);
// The enhanced routing heuristic per edge (Section 12), nonnegative with a denominator of at
// least 1. Admissibility is separate: forbidden edges are excluded by the searches, never
// priced out with a large finite cost.
CLEANSEAMS_CORE_EXPORT std::vector<double> enhanced_edge_cost(const MeshSnapshot&,const ResolvedPreferences&,const PreferenceOptions&);
// The host attribute names a run transfers (documented names or configured aliases); none
// unless the options' policy is enhanced_v1, the only policy that reads them (GB M2).
CLEANSEAMS_CORE_EXPORT std::vector<std::string> painted_attribute_names(const PipelineOptions&);
// BLAKE2b-256 hex digest of a float64 array (report identity of the raw exposure used).
CLEANSEAMS_CORE_EXPORT std::string values_digest(std::span<const double> values);
// The first generator of the Section 21 schedule: when candidate zero violates the painted
// hard constraints (or leaves a chart open after them), it proposes a reroute of the
// forbidden cuts and a weld-and-reopen alternative, along admissible cheapest paths.
CLEANSEAMS_CORE_EXPORT ProposalGenerator constraint_repair_generator(std::shared_ptr<const ResolvedPreferences>,PreferenceOptions);
// ------------------------------------------------------------ compound cuts (T34-T37)
// Shared by the region, ring, cavity and shell generators (src/proposals/compound.cpp).
// A compound proposal carries every companion cut it needs: the generator builds the whole
// trial cut set, opens every non-disk chart it touches along admissible paths and emits one
// delta; the transaction then validates all affected charts together.
struct CLEANSEAMS_CORE_EXPORT CompanionOpening {
    bool ok=true;                               // every chart in scope is a measured-disk or exempt-only
    bool stopped=false;                         // GB m7: the control's deadline stopped the opening (not ok)
    std::vector<MeshIndex> added;               // opening edges added to the cut set (ascending)
    nlohmann::json failures=nlohmann::json::array(); // charts left open: {chart, faces, reason}
    nlohmann::json to_json() const;
};
// Cheapest path over mesh edges from any source vertex to any target vertex (Dijkstra in
// (distance, vertex) order, deterministic). `usable` filters edges; `cost` prices them.
CLEANSEAMS_CORE_EXPORT std::vector<MeshIndex> cheapest_edge_path(const MeshSnapshot&,std::span<const MeshIndex> sources,
    std::span<const std::uint8_t> targets,const std::function<bool(MeshIndex)>& usable,std::span<const double> cost);
// Exact disk opening of every measured chart of `cuts` holding a face of `scope` (every
// chart when `scope` is empty): two or more boundary loops are joined along the cheapest
// admissible interior path; closed or higher-genus charts get a cut graph from a dual
// spanning tree that keeps forbidden edges uncut and expensive edges whole, pruned to its
// cycles. Admissible: not forbidden, not a physical boundary. Charts that stay non-disk
// (pinched, unsupported, or no admissible opening) are reported, never forced. GB m7: each
// round rebuilds the surface, so `control` is checked before every round: a cancellation
// throws, an expired deadline stops with `ok` false, `stopped` and a "deadline" failure.
CLEANSEAMS_CORE_EXPORT CompanionOpening open_non_disk_charts(const std::shared_ptr<const MeshSnapshot>&,EdgeMask& cuts,std::span<const double> cost,
    const HardConstraints&,std::span<const std::uint8_t> exempt,std::span<const std::uint8_t> scope={},const Control& control={});
// GB (visual review V1): the baseline's seam-chain refinement on a compound proposal. Region
// growth and the level-set or label borders of rings, cavities and shells are face-region
// boundaries that follow the triangulation's staircase; the baseline refines its own chains
// (pipeline refine step) and these were never refined. Every new chain of `cuts` (the seams the
// parent keeps are locked like the baseline's locks, so chains split where the new boundary meets
// them) with every edge inside `scope` (all faces when empty) goes through refine_chain /
// refine_cycle exactly as the baseline's refine step runs them: the same corridors (hop rings; T39:
// ring operators under the frozen scale), weight variants, turn/flow/drift/convexity terms, the
// crease-discounted seam cost of the raw exposure (seam_cost, feature and valley discounts), the
// native_baseline routing fixes, a required improvement, up to min(3, refine_passes) sweeps with
// the dirty-vertex deferral. Hard constraints and the proposal's topology are kept: locked and
// required chains are never moved, forbidden and out-of-scope edges are outside every corridor,
// manual anchors end chains, each accepted route keeps the chart count, adjacency, disks and
// pinches (topology_preserved), a route that would glue differently pinned sectors is not taken,
// and constraint_violations is checked on the result (a violation keeps `cuts` unchanged). No
// solve and no flanking guard: the transaction evaluates the refined proposal like any other.
// Off (enhanced.boundary_refinement.enabled, or refine_paths false) `cuts` is left unchanged.
struct CLEANSEAMS_CORE_EXPORT BoundaryRefinement {
    bool enabled=false,stopped=false,reverted=false;
    std::uint64_t chains=0,refined=0,deferred=0,pin_refusals=0,passes=0;
    std::int64_t corridor_rings=0,corridor_max_rings=0;
    std::uint64_t edges_before=0,edges_after=0,turns_before=0,turns_after=0; // of the refined chains
    double length_before=0,length_after=0;
    nlohmann::json to_json() const;
};
CLEANSEAMS_CORE_EXPORT BoundaryRefinement refine_compound_boundaries(ProposalContext&,const RunState& parent,EdgeMask& cuts,
    std::span<const std::uint8_t> scope={});
// The delta from `parent.cuts` to `cuts` as one proposal (physical boundary edges never listed).
CLEANSEAMS_CORE_EXPORT Proposal cut_proposal(std::uint64_t id,const RunState& parent,std::span<const std::uint8_t> cuts,std::string generator,
    std::string version,nlohmann::json provenance);
// The context's resolved paint, or a neutral record (no paint) over its raw exposure.
CLEANSEAMS_CORE_EXPORT std::shared_ptr<const ResolvedPreferences> context_preferences(const ProposalContext&);
// The enhanced routing heuristic of the context (enhanced_edge_cost on context_preferences).
CLEANSEAMS_CORE_EXPORT std::vector<double> context_routing_cost(const ProposalContext&);
// The context's BVH (T36/T37): built once per geometry revision and pass, on first use;
// rebuilt only when `mesh` is another geometry. Throws memory_limit when it does not fit
// the run's memory budget (the detectors then skip with a diagnostic).
CLEANSEAMS_CORE_EXPORT std::shared_ptr<const TriangleBvh> context_bvh(ProposalContext&);
// ------------------------------------------------------------ preference-aware regions (T34)
// Section 15. Preferred seam strength per edge, in [0, 1], grows where a seam is welcome:
//   m(e) = enhanced_edge_cost(e) / length(e)        (the Section 12 routing multiplier)
//   strength(e) = 1 - min(1, m(e) / m_ref)          m_ref: median of m over `reference`
// A uniform field has strength 0 everywhere, so preference-aware growth then reproduces the
// baseline growth exactly. It is a dual crossing weight (expensive to cross = welcome to
// cut), never the primal routing cost.
constexpr const char* region_strength_version="region_strength.v1";
CLEANSEAMS_CORE_EXPORT std::vector<double> preferred_seam_strength(const MeshSnapshot&,const ResolvedPreferences&,const PreferenceOptions&,
    std::span<const MeshIndex> reference);
// Section 21 "preference-aware split alternatives": for up to enhanced.regions.max_targets
// targets (over-budget charts, then adjacent chart pairs whose shared seam is not hard), the
// target's faces are re-partitioned by weighted region growth over the must-weld contracted
// dual (region_dual; required/locked cuts are barriers) in two variants: `preference`
// (w_preference as configured) and `baseline` (w_preference = 0: the original farthest-seed /
// Lloyd routine). Non-hard cuts inside the target are replaced by the region borders plus exact
// disk openings; partitions that cannot be opened admissibly are refused with diagnostics
// (ctx.diagnostics["region_preferences"]). Identical variants are proposed once.
CLEANSEAMS_CORE_EXPORT ProposalGenerator region_preference_generator();
// ------------------------------------------------------------ generalized rings (T35)
// Section 16: constriction and crease loops from geometry alone (no object front, world-up or
// anatomy label). On a patch (a chart's faces, its interior cuts ignored: physical surface
// distance) the physical distance from its opening (the longest boundary component) is
// sampled at section_samples quantiles in [low_quantile, high_quantile] of the maximum; the
// borders between sublevel and superlevel face sets are real mesh edges, and only simple
// closed loops are kept. Every loop is ranked (ring_rank.v1):
//   constriction = 1 - length / min(widest section before, widest section after)   (a waist)
//   crease       = length-weighted feature desirability along the loop
//   exposure     = length-weighted raw exposure along the loop
//   balance      = smaller side area / patch area;  compactness = min(1, sqrt(smaller side area) / length)
//   confidence   = max(constriction, crease)
//   score        = constriction + 0.5 crease + 0.25 balance + 0.25 compactness - 0.5 exposure
// Ranking only orders and bounds the candidates; the common transaction judges them.
constexpr const char* ring_rank_version="ring_rank.v1";
struct CLEANSEAMS_CORE_EXPORT RingCandidate {
    std::vector<MeshIndex> edges;        // ascending
    std::vector<MeshIndex> vertices;     // in loop order
    double level=0,quantile=0,length=0;
    std::array<double,2> side_area{};
    std::array<std::uint64_t,2> side_faces{};
    double constriction=0,crease=0,exposure=0,balance=0,compactness=0,confidence=0,score=0;
    // candidate, or: non_separating, tiny_side (below the minimum-island policy), forbidden_edge,
    // low_confidence; the generator adds proposed, not_openable, constraint, beyond_candidate_budget.
    std::string status;
    nlohmann::json to_json() const;
};
struct CLEANSEAMS_CORE_EXPORT RingSearch {
    std::vector<MeshIndex> opening;      // edges of the opening
    double opening_length=0,far_length=0,max_distance=0;
    std::vector<double> levels,section_length; // per sampled level: distance and level-set length
    std::uint64_t non_simple_components=0;
    std::vector<RingCandidate> loops;    // every distinct simple loop, candidates first, best first
    nlohmann::json island_floor=nlohmann::json::object(); // T39: the minimum-island floor applied to the sides
    bool stopped=false;                  // GB m7: the pass's deadline stopped the search (partial)
    nlohmann::json to_json() const;
};
// The vertices of a simple closed loop in order, or empty when `edges` is not one.
CLEANSEAMS_CORE_EXPORT std::vector<MeshIndex> simple_loop_order(const MeshSnapshot&,std::span<const MeshIndex> edges);
// GB m7: the per-pass work of the per-chart detectors (rings, cavities, shells): mesh-level
// quantities computed once per pass rather than once per chart (every other array a detector
// uses is local to its chart or component), and the pass's control, checked inside the
// detectors' loops: a cancellation throws `cancelled`, an expired deadline stops the detector,
// which returns what it found so far with `stopped` set (the generator then stops).
struct CLEANSEAMS_CORE_EXPORT DetectorPass {
    Control control;
    std::optional<double> mesh_area;                  // total face area; computed on first use
    std::vector<double> component_diagonals;          // per connected component; computed on first use
    bool stopped=false;                               // a detector stopped at the deadline
    double area(const MeshSnapshot&);
    std::span<const double> diagonals(const MeshSnapshot&);
    // Throws cancelled when requested; true (and `stopped`) when the deadline has passed.
    bool stop() {stopped=stopped || control.check();return stopped;}
};
// `scale` (T39): the pass's physical scale; the minimum-island face floor of a side is then an
// area floor (area-weighted region preservation). Null: the legacy face count. A side is tiny
// under the complete minimum-island policy (GB m11: layout's below_island_policy: face count
// with the large n-gon exception, mesh area fraction and thin strips; face counts are area
// floors under the physical scale). `pass` (GB m7): the pass's shared work and control (null:
// computed here, no control); the search is local to `faces`.
CLEANSEAMS_CORE_EXPORT RingSearch ring_search(const MeshSnapshot&,std::span<const MeshIndex> faces,const HardConstraints&,
    std::span<const double> feature_desirability,std::span<const double> raw_exposure,const PipelineOptions&,const PhysicalScale* scale=nullptr,
    DetectorPass* pass=nullptr);
// Section 21 "confidence-gated ring alternatives": for up to enhanced.rings.max_charts charts
// (best loop first), up to max_candidates loops each become one compound proposal: the loop
// plus the chart's companion openings. The chart's non-hard interior cuts are dropped and
// every side is reopened along admissible cheapest paths on which the old openings are free,
// so an existing slit is reused where it still fits; a ring that leaves a closed or annular
// side without an admissible opening is refused (never a partial candidate).
CLEANSEAMS_CORE_EXPORT ProposalGenerator ring_generator();
// ------------------------------------------------------------ cavities (T36)
// Section 17: recessed connected regions, confidence-gated. Nomination: up to max_samples
// area-distributed faces of a chart each cast rays_per_sample deterministic cone rays
// (cone_half_angle_deg about the face normal, azimuths evenly spaced in a frame fixed by the
// face's first edge, so the rays move with the geometry) on the geometry-revisioned BVH, with
// a reach of reach_fraction of the sample's connected-component bounding-box diagonal. Only
// faces of the same connected component count as hits, so a nearby disconnected sheet never
// closes a pocket; escaped rays lower the occlusion and never prove a cavity.
// Refinement on the original mesh: faces blocked at least occlusion_threshold of the time
// (holes filled, seeds never removed) form cores; the core's area-weighted normal is the
// opening axis; sublevel sets of the height along it grow a region from the deepest core
// face, and the largest of level_samples levels whose region stays inside the chart with one
// simple border loop and at most max_area_fraction of the chart area is the rim. The seam is
// the nearest such loop inside the rim (the rim itself when there is none); the rim stays an
// alternative loop. Evidence:
//   depth_ratio = (mean rim height - lowest region height) / (rim length / pi)
//   concavity   = area share of the region facing the opening (towards the rim centroid):
//                 the normal-consistency term (inverted or inconsistent normals lower it)
//   confidence  = min(1, occlusion / 0.5) * min(1, depth_ratio / (2 min_depth_ratio)) * concavity
// A loop's interior slit runs inside the loop's region from its deepest interior vertex to the
// loop; it is needed when that region is itself deep: its depth over the loop's diameter at
// least slit_depth_ratio (and it has an interior vertex to start from).
// Detector filters only generate proposals; they are never exemptions from validation.
constexpr const char* cavity_detector_version="cavity_detector.v1";
struct CLEANSEAMS_CORE_EXPORT CavityOcclusion {
    std::vector<MeshIndex> samples;            // sampled faces, ascending
    std::map<MeshIndex,double> occlusion;      // per sampled face: fraction of blocked rays
    std::map<MeshIndex,double> reach;          // per sampled component: the ray reach used
    std::uint64_t rays=0,hits=0;
    double seconds=0;
    nlohmann::json to_json() const;
};
// Throws stale_state when the BVH does not serve `mesh`'s geometry revision. `diagonals`
// (GB m7): component_diagonals of `mesh`, computed once per pass; empty: computed here.
CLEANSEAMS_CORE_EXPORT CavityOcclusion cavity_occlusion(const MeshSnapshot&,const TriangleBvh&,std::span<const MeshIndex> faces,
    const CavityOptions&,bool same_component=true,std::span<const double> diagonals={});
struct CLEANSEAMS_CORE_EXPORT CavityRecord {
    std::vector<MeshIndex> core,region;        // faces: the occluded core and the region inside the rim
    std::vector<MeshIndex> rim,seam,slit;      // edges: rim loop, seam loop (inside the rim when possible), interior slit
    Vec3 axis{0,0,0};
    double occlusion=0,concavity=0,depth=0,opening_diameter=0,depth_ratio=0,area_fraction=0,confidence=0;
    double seam_depth_ratio=0;                 // the seam region's depth over its loop's diameter
    bool needs_slit=false,seam_inside_rim=false; // needs_slit: seam_depth_ratio >= slit_depth_ratio with an interior vertex
    // The rim itself as an alternative loop when the seam is inside it (and not forbidden).
    bool rim_alternative=false,rim_needs_slit=false;
    std::vector<MeshIndex> rim_slit;
    // candidate, or: no_opening_axis, incoherent_opening (no level with one simple border loop
    // inside the chart: open sheets, through-holes, recesses reaching a seam), too_small,
    // shallow, forbidden_edge, low_confidence
    std::string status;
    nlohmann::json to_json() const;
};
struct CLEANSEAMS_CORE_EXPORT CavityDetection {
    CavityOcclusion occlusion;
    std::uint64_t seeds=0,cores=0;
    std::vector<CavityRecord> cavities;         // at most max_cavities, strongest core first
    nlohmann::json scale=nlohmann::json::object(); // T39: the seed-filling passes under the physical scale
    bool stopped=false;                         // GB m7: the pass's deadline stopped the detection (partial)
    nlohmann::json to_json() const;
};
// Bounding-box diagonal of every connected component (the detector's scale).
CLEANSEAMS_CORE_EXPORT std::vector<double> component_diagonals(const MeshSnapshot&);
// Detection on one chart (`chart_faces`) of the cut state `cuts`; `routing` prices the slit.
// `scale` (T39): the pass's physical scale; the seed filling then has a physical radius (null:
// three face-ring passes). `pass` (GB m7): the pass's shared work (component diagonals) and
// control, checked per core (null: computed here, no control); every array is local to the chart.
CLEANSEAMS_CORE_EXPORT CavityDetection detect_cavities(const MeshSnapshot&,const TriangleBvh&,std::span<const MeshIndex> chart_faces,
    std::span<const std::uint8_t> cuts,const HardConstraints&,std::span<const double> routing,const PipelineOptions&,bool same_component=true,
    const PhysicalScale* scale=nullptr,DetectorPass* pass=nullptr);
// Section 21 "confidence-gated cavity alternatives": for up to max_cavities candidate
// cavities per pass (charts of at least min_chart_faces faces), compound proposals `seam` (the
// seam loop) and, when the seam lies inside the rim, `rim` (the rim loop itself), each with a
// `_slit` variant (the loop plus its interior slit) when the loop's region needs one. Every
// one carries the companion openings every changed chart needs (the annular outer side is
// joined to its boundary). Existing seams are never removed and no exemption is added: the transaction's
// strict overlap check refuses any overlap instead of trading it for a hidden seam.
CLEANSEAMS_CORE_EXPORT ProposalGenerator cavity_generator();
// ------------------------------------------------------------ thin shells (T37)
// Section 18. Per connected component: an area-weighted PCA of the face centroids gives the
// thin axis and the principal full extents (vertex ranges along the axes). The nominating
// thickness ratio is the local thickness (area-weighted median opposite-side distance) over
// the middle extent (the smallest extent when no opposite side answers); below thin_ratio it
// only nominates. The global PCA spread ratio is reported (it counts a bend as thickness).
// Confirmation: front and back normal alignment with the thin axis (|n.t| > 0.5 on both sides),
// opposite-side correspondence (a ray from a front or back face into the solid meets an
// opposite-facing face of the same component within four standard deviations along the thin
// axis), front/back labels
// from 0.5 n.t plus the position relative to the locally estimated middle surface (propagated
// from the corresponding faces), smoothed without crossing required or locked edges, label
// speckles below speckle_fraction of the area repaired geometrically, and a rim (the label
// border, real edges) made of simple closed loops.
//   confidence = sqrt(max(0, 1 - ratio/thin_ratio)) * correspondence * min(1, min(front, back)/0.3) * [rim coherent]
// Reasons of a rejection: not_thin, no_opposite_sides, low_correspondence, incoherent_rim,
// forbidden_rim, low_confidence, too_few_faces.
struct CLEANSEAMS_CORE_EXPORT ShellEvidence {
    MeshIndex component=0;
    std::uint64_t faces=0,correspondence_samples=0,speckles_repaired=0,rim_loops=0;
    double area=0;
    std::array<double,3> extents{};    // full vertex ranges along the thin, middle and largest principal axes
    Vec3 thin_axis{0,0,0};
    double thin_ratio=1,pca_ratio=1,front_fraction=0,back_fraction=0,edge_on_fraction=0,correspondence=0,thickness=0,confidence=0;
    bool nominated=false,rim_simple=false;
    std::vector<std::int8_t> labels;   // per face of the component (ascending face ids): 1 front, 0 back
    std::vector<MeshIndex> rim,biased_rim;
    std::string status;                // candidate or rejected
    std::vector<std::string> reasons;
    nlohmann::json scale=nlohmann::json::object(); // T39: face floor, label smoothing and rim offset under the physical scale
    nlohmann::json to_json() const;
};
// `routing` (per edge, optional) prices the visibility-biased nearby rim: the front grown or
// shrunk by one face ring, kept when its mean routing multiplier is lower than the raw rim's.
// Throws stale_state when the BVH does not serve `mesh`'s geometry revision.
// `scale` (T39): the pass's physical scale: the component face floor becomes an area floor, the
// label smoothing and the nearby rim's offset get physical radii. Null: the legacy counts.
// `pass` (GB m7): the pass's control, checked before every component and smoothing pass (a
// deadline returns the components examined so far with `pass->stopped` set); every array is
// local to the component.
CLEANSEAMS_CORE_EXPORT std::vector<ShellEvidence> shell_evidence(const MeshSnapshot&,const TriangleBvh&,const HardConstraints&,
    std::span<const double> routing,const PipelineOptions&,const PhysicalScale* scale=nullptr,DetectorPass* pass=nullptr);
// Section 21 "confidence-gated shell alternatives": at most max_candidates alternatives per
// component, candidate zero included (the geometric rim and the visibility-biased rim). The
// component's non-hard cuts are dropped (free for any companion opening), the rim added, every
// side reopened if needed; the ordinary transaction solves and validates both sides, all
// remaining pieces and every pin. There is no early acceptance.
CLEANSEAMS_CORE_EXPORT ProposalGenerator shell_generator();
}
