// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
#include "cleanseams/evaluator.hpp"
#include "cleanseams/fields.hpp"
#include <set>
namespace cleanseams {
// One validated chart map of the uv cache, with the exact chart identity it was
// measured for (lookups compare the structure, not only the digest).
struct CachedChartMap {
    ChartKey key;
    std::vector<MeshIndex> corner_ids; // strictly increasing
    std::vector<Vec2> corner_uv;
};
// Source quad_layout.LayoutReport. Seeding (T12) and merging (T13) fill it.
struct CLEANSEAMS_CORE_EXPORT LayoutReport {
    bool timed_out=false;
    std::string tier="quad";
    std::uint64_t seeds=0;
    double dijkstra_s=0;
    std::uint64_t paired_sources=0,sources=0,seeds_added=0,trace_passes=0,patches=0,charts=0,rounds=0;
    std::uint64_t merges_accepted=0,merges_rejected_proxy=0,merges_rejected_topology=0,merges_rejected_distortion=0;
    std::uint64_t merges_rejected_compactness=0,merges_rejected_feature=0,merges_rejected_strip=0,merges_rejected_isoperimetric=0;
    double seam_exposure=std::numeric_limits<double>::quiet_NaN();
    std::uint64_t seeds_by_exposure=0;
    nlohmann::json revolution_cuts=nlohmann::json::array();
    std::uint64_t shoulder_edges=0,single_arc_merges=0,self_cuts_removed=0,spurs_removed=0,spurs_kept=0,feature_edges=0;
    std::vector<std::vector<MeshIndex>> feature_chains;
    std::vector<MeshIndex> feature_panel_faces,open_features;
    std::optional<std::vector<std::uint8_t>> protected_seams;
    std::uint64_t seam_edges=0;
    double seam_length=0;
    std::uint64_t non_disk_charts=0;
    nlohmann::json notes=nlohmann::json::array(); // strings, or dicts such as the CAD panel summary
    std::vector<std::pair<MeshIndex,ChartMeasurement>> chart_metrics;
    nlohmann::json rejections=nlohmann::json::array();
    std::map<std::string,std::uint64_t> stop_reasons;
    nlohmann::json stop_pairs=nlohmann::json::array();
    std::uint64_t last_round_rejections=0,max_rejections=2000,splits=0,absorbed=0,absorbed_relaxed=0,absorbed_exempt=0;
    std::uint64_t unsupported_absorbed=0,premerged_nonquad=0;
    nlohmann::json tiny_remaining=nlohmann::json::array(); // {"chart","faces","reason"} per tiny chart left
    std::vector<MeshIndex> absorbed_faces,exempt_faces;
    std::map<std::string,CachedChartMap> uv_cache; // validated chart maps by ChartKey::digest()
    std::optional<std::vector<Vec2>> uv; // (nc) map assembled from the validated charts, NaN where missing
    std::set<ArcKey> frozen_keys; // arc keys of reverted moves
    std::map<ArcKey,std::string> frozen_reasons;
    nlohmann::json as_dict() const;
    // A validated map for this exact chart (key digest and structure), if cached.
    const CachedChartMap* cached(const ChartKey&) const;
    void cache(const ChartKey&,std::vector<MeshIndex> corner_ids,std::vector<Vec2> corner_uv);
};
// Ordered stage events in the source's call-completion order, for parity tests
// and diagnostics. Recording is optional and never changes a decision.
struct LayoutTrace {nlohmann::json events=nlohmann::json::array();};
using EdgeMask=std::vector<std::uint8_t>;
// User pins {corner: (u, v)}.
using PinMap=std::map<MeshIndex,Vec2>;
// Source progress callback: stage events such as {"stage": "generate", "status":
// "running", ...}; the pipeline's sink also checks cancellation.
using LayoutProgress=std::function<void(const nlohmann::json&)>;
// PIN-03: gluing `edges` into `cs` would join sectors whose corners are pinned to
// different coordinates. Sectors join as a rebuilt surface joins corners across the
// uncut manifold edge; edges a surface always cuts are skipped, and a conflict
// already inside one sector is not attributed to the move.
CLEANSEAMS_CORE_EXPORT bool pins_conflict_after_gluing(const CutSurface&,const PinMap& pins,std::span<const MeshIndex> edges);
using CostField=std::optional<std::span<const double>>;
// ------------------------------------------------------------ motorcycle layout
CLEANSEAMS_CORE_EXPORT EdgeMask irregular_sources(const MeshSnapshot&,std::span<const std::uint8_t> cut,bool relaxed=false,LayoutTrace* trace=nullptr);
struct PairedSources {std::vector<std::uint8_t> sources;EdgeMask cut;std::uint64_t count=0;};
CLEANSEAMS_CORE_EXPORT PairedSources pair_sources(const MeshSnapshot&,std::span<const std::uint8_t> sources,std::span<const std::uint8_t> cut,LayoutTrace* trace=nullptr);
CLEANSEAMS_CORE_EXPORT TransferIndex motorcycle_continuation(const MeshSnapshot&,MeshIndex vertex,MeshIndex incoming,bool relaxed=false,double max_angle=35.0);
struct MotorcycleTrace {EdgeMask cut;std::uint64_t motorcycles=0;};
CLEANSEAMS_CORE_EXPORT MotorcycleTrace trace_motorcycles(const MeshSnapshot&,std::span<const std::uint8_t> sources,std::span<const std::uint8_t> cut,
    bool relaxed=false,LayoutTrace* trace=nullptr);
CLEANSEAMS_CORE_EXPORT std::vector<double> angle_defect(const MeshSnapshot&);
// Source feature_cuts: detect_features with the layout's scale/ratio, stored in the report.
CLEANSEAMS_CORE_EXPORT EdgeMask feature_cuts(const std::shared_ptr<const MeshSnapshot>&,double angle_deg,std::int64_t min_chain,std::int64_t scale,
    double nms_ratio,std::int64_t max_extend,std::optional<std::span<const std::uint8_t>> barriers,const Deadline&,LayoutReport*,LayoutTrace* trace=nullptr);
struct DiskLayout {EdgeMask cut;std::shared_ptr<const CutSurface> cs;};
CLEANSEAMS_CORE_EXPORT DiskLayout make_all_disks(const std::shared_ptr<const MeshSnapshot>&,std::span<const std::uint8_t> cut,
    std::span<const std::uint8_t> sources,LayoutReport&,std::int64_t max_passes=12,CostField cost=std::nullopt,const Deadline& deadline={},
    bool relaxed=false,LayoutTrace* trace=nullptr);
CLEANSEAMS_CORE_EXPORT std::vector<MeshIndex> trace_from(const MeshSnapshot&,std::span<const std::uint8_t> cut,MeshIndex vertex,std::uint64_t max_steps=100000);
CLEANSEAMS_CORE_EXPORT std::optional<MeshIndex> cheapest_seed(const MeshSnapshot&,std::span<const std::uint8_t> cut,std::vector<MeshIndex> candidates,
    CostField cost,std::size_t max_candidates=64,LayoutTrace* trace=nullptr);
CLEANSEAMS_CORE_EXPORT std::optional<MeshIndex> seed_vertex(const CutSurface&,MeshIndex chart,std::span<const std::uint8_t> sources,CostField cost,LayoutTrace* trace=nullptr);
// {loop: {sector: uncut-edge distance to the chart's other loops}}.
CLEANSEAMS_CORE_EXPORT std::map<MeshIndex,std::map<MeshIndex,std::int64_t>> distance_from_other_loops(const CutSurface&,MeshIndex chart);
CLEANSEAMS_CORE_EXPORT std::optional<MeshIndex> interior_seed(const CutSurface&,MeshIndex chart,std::span<const std::uint8_t> sources,CostField cost,LayoutTrace* trace=nullptr);
struct SplitLayout {EdgeMask cut;std::shared_ptr<const CutSurface> cs;std::vector<std::uint8_t> sources;};
CLEANSEAMS_CORE_EXPORT SplitLayout split_over_budget(const std::shared_ptr<const MeshSnapshot>&,std::span<const std::uint8_t> cut,
    const std::shared_ptr<const CutSurface>& cs,std::span<const std::uint8_t> sources,const LayoutOptions&,LayoutReport&,const ChartEvaluator*,
    CostField cost=std::nullopt,const Deadline& deadline={},bool segmented=false,LayoutTrace* trace=nullptr);
// ------------------------------------------------------------ split records (T32, split.cpp)
// A broad over-budget split of the baseline, recorded under enhanced_v1 so that relief
// proposals can be evaluated before it (Sections 14.3 and 21). Recording is read-only: the
// split itself, and every legacy_010/native_baseline decision, is unchanged.
struct CLEANSEAMS_CORE_EXPORT SplitRecord {
    std::string stage;                          // "final_chart_gate"
    std::shared_ptr<const CutSurface> surface;  // the surface before the split
    std::vector<Vec2> uv;                       // its complete map (per original corner)
    MeshIndex chart=0;                          // the over-budget chart of `surface`
    std::vector<MeshIndex> faces;               // its faces (sorted)
    std::vector<MeshIndex> added_edges;         // seam edges the split added inside it (sorted)
    double anisotropy_p95=std::numeric_limits<double>::quiet_NaN(); // the gate's measurement before the split
    nlohmann::json summary() const;
};
// The records of one split: for every over-budget chart of `before` (the surface the split
// started from, with its map), the edges `after` cuts inside that chart and `before` did
// not. Charts the split did not change are left out.
CLEANSEAMS_CORE_EXPORT std::vector<SplitRecord> record_split(const std::string& stage,const std::shared_ptr<const CutSurface>& before,
    std::span<const Vec2> uv,const CutSurface& after,std::span<const MeshIndex> over_charts,std::span<const double> anisotropy_p95={});
// Split edges of `record` that `cuts` still cuts (the edges relief before the split welds).
CLEANSEAMS_CORE_EXPORT std::vector<MeshIndex> split_edges_present(const SplitRecord&,std::span<const std::uint8_t> cuts);
// ------------------------------------------------------------ dead-end arcs (T33, merge.cpp)
// A dead-end seam arc, by the legacy dangling-end definition (dangling_seam_ends): from a
// seam vertex of degree one that is not on the mesh boundary, along seam vertices of degree
// two, to the first junction or boundary vertex, or to a second dangling end (a floating
// slit). The ends of all arcs are exactly the dangling ends. cleanup_spurs is unchanged.
struct DeadEndArc {
    std::vector<MeshIndex> tips;                // one dangling end; two for a floating slit
    std::vector<MeshIndex> edges;               // from the first tip to the junction
    std::vector<MeshIndex> vertices;            // first tip first
    MeshIndex junction=0;                       // the last vertex (the second tip when floating)
    bool floating=false;
};
CLEANSEAMS_CORE_EXPORT std::vector<DeadEndArc> dead_end_arcs(const MeshSnapshot&,std::span<const std::uint8_t> seams);
// ------------------------------------------------------------ region growing (segment.py)
struct DualGraph {
    std::vector<Offset> indptr{0};
    std::vector<MeshIndex> indices,edges;
    std::vector<double> weights;
    std::size_t size() const {return indptr.size()-1;}
};
CLEANSEAMS_CORE_EXPORT DualGraph dual_graph(const MeshSnapshot&,std::optional<std::span<const std::uint8_t>> barriers=std::nullopt,double w_crease=4.0);
CLEANSEAMS_CORE_EXPORT void update_distance(const DualGraph&,MeshIndex seed,std::vector<double>& distance);
CLEANSEAMS_CORE_EXPORT std::vector<MeshIndex> farthest_point_seeds(const DualGraph&,std::uint64_t count,
    std::optional<std::vector<MeshIndex>> initial=std::nullopt,const Deadline& deadline={},LayoutTrace* trace=nullptr);
CLEANSEAMS_CORE_EXPORT std::vector<std::int64_t> grow_labels(const DualGraph&,std::span<const MeshIndex> seeds,
    std::optional<std::span<const Vec3>> normals=std::nullopt,std::optional<std::span<const Vec3>> proxies=std::nullopt,double w_normal=2.0,LayoutTrace* trace=nullptr);
struct LloydResult {std::vector<std::int64_t> labels;std::vector<MeshIndex> seeds;};
CLEANSEAMS_CORE_EXPORT LloydResult lloyd(const MeshSnapshot&,const DualGraph&,std::vector<MeshIndex> seeds,std::int64_t steps=3,const Deadline& deadline={},LayoutTrace* trace=nullptr);
CLEANSEAMS_CORE_EXPORT EdgeMask regions_to_cut(const MeshSnapshot&,std::span<const std::int64_t> labels,std::optional<std::span<const std::uint8_t>> initial_cut=std::nullopt);
CLEANSEAMS_CORE_EXPORT DiskLayout make_region_disks(const std::shared_ptr<const MeshSnapshot>&,std::span<const std::uint8_t> cut,CostField cost=std::nullopt,
    const Deadline& deadline={},LayoutTrace* trace=nullptr);
struct Segmentation {EdgeMask cut;std::shared_ptr<const CutSurface> cs;std::uint64_t seeds=0,regions=0;double dijkstra_s=0;};
CLEANSEAMS_CORE_EXPORT Segmentation segment_mesh(const std::shared_ptr<const MeshSnapshot>&,std::optional<std::span<const std::uint8_t>> initial_cut=std::nullopt,
    std::int64_t seed_density=400,std::int64_t max_patches=4000,std::int64_t steps=3,CostField cost=std::nullopt,const Deadline& deadline={},
    std::optional<std::vector<MeshIndex>> seeds=std::nullopt,LayoutTrace* trace=nullptr);
CLEANSEAMS_CORE_EXPORT DiskLayout split_regions(const std::shared_ptr<const MeshSnapshot>&,const CutSurface& surface,std::span<const MeshIndex> charts,
    CostField cost=std::nullopt,const Deadline& deadline={},LayoutTrace* trace=nullptr);
// ------------------------------------------------------------ preference-aware regions (T34)
// enhanced_v1 only; the legacy functions above are unchanged. The dual graph of a face set
// with must-weld contraction: faces joined across `weld` edges (painted forbidden seams) are
// one node, so no partition of the nodes can cut a welded edge; `barriers` (required and
// locked cuts, physical boundaries) are never crossed. The crossing weight of a face pair is
//   centroid_distance * (1 + w_normal*(1 - n_a.n_b) + w_preference*strength(e)^2)
// (Section 15): a preferred seam is expensive to cross, so grown regions meet there. Between
// two nodes the cheapest crossing is kept (ties: lower edge id). With w_preference = 0 and
// no welds the node weights are exactly those of dual_graph(mesh, barriers, w_normal).
struct CLEANSEAMS_CORE_EXPORT RegionDual {
    DualGraph graph;                        // nodes: must-weld groups, numbered by smallest face
    std::vector<MeshIndex> faces;           // the region's faces, ascending
    std::vector<MeshIndex> node_of_face;    // per entry of `faces`
    IndexRows node_faces;                   // faces of every node, ascending
    std::vector<double> node_area;
    std::vector<Vec3> node_center,node_normal; // area-weighted centre and unit normal
    std::uint64_t welded_edges=0;           // weld edges inside the region
    MeshIndex node(MeshIndex face) const;   // node of a region face (throws for another face)
};
CLEANSEAMS_CORE_EXPORT RegionDual region_dual(const MeshSnapshot&,std::span<const MeshIndex> faces,std::span<const std::uint8_t> barriers,
    std::span<const std::uint8_t> weld,std::span<const double> strength,double w_normal,double w_preference);
// Lloyd iterations of lloyd() on the nodes of a region dual: competing growth from `seeds`,
// then seeds move to the node closest to each region's area-weighted centre and normal
// (the same score), growth with the normal-proxy penalty, until stable or `steps`.
// `quantum` (T39, physical mode): growth distances are kept on this length grid and the seed
// scores on a 2^-32 grid, so exact geometric ties stay exact under rigid motions and uniform
// scales and are broken by (label, node); 0 is the legacy growth (grow_labels) exactly.
CLEANSEAMS_CORE_EXPORT LloydResult region_lloyd(const RegionDual&,std::vector<MeshIndex> seeds,std::int64_t steps=3,const Deadline& deadline={},
    double quantum=0.0);
// ------------------------------------------------------------ absorb helpers used before splitting
CLEANSEAMS_CORE_EXPORT bool topological_disk(const CutSurface&,MeshIndex chart);
CLEANSEAMS_CORE_EXPORT bool large_planar_ngons(const MeshSnapshot&,std::span<const MeshIndex> faces,std::int64_t minimum=8);
CLEANSEAMS_CORE_EXPORT bool small_face_count(const MeshSnapshot&,std::span<const MeshIndex> faces,std::int64_t minimum);
CLEANSEAMS_CORE_EXPORT DiskLayout premerge_single_nonquad(const std::shared_ptr<const MeshSnapshot>&,std::shared_ptr<const CutSurface> cs,LayoutReport&,
    std::optional<std::span<const std::uint8_t>> locked=std::nullopt,const Deadline& deadline={},std::int64_t minimum_faces=8,double minimum_fraction=.002,
    LayoutTrace* trace=nullptr,const PinMap& pins={});
// ------------------------------------------------------------ generation up to merging
// The first half of source generate_seams: feature/CAD/lathe/revolution cuts,
// irregular sources, motorcycle or region seeding, single n-gon premerge and
// over-budget splitting. T13's merge stage continues from this state.
struct SeededLayout {
    EdgeMask cut; // after splitting, with locked edges added (source merge input)
    std::shared_ptr<const CutSurface> cs;
    std::vector<std::uint8_t> sources;
    EdgeMask features,refine_locked,user_cut;
    std::optional<EdgeMask> locked;
    std::vector<double> cost; // empty when exposure is disabled
    bool segmented=false;
    LayoutReport report;
    ChartEvaluator evaluate; // effective evaluator (empty when none)
    std::string next_stage; // absorb_tiny (unsupported charts) or merge_patches
    // Set by the deadline-bounded batched evaluator when it stops early; merged into
    // the report (source: the wrapper sets report.timed_out).
    std::shared_ptr<bool> evaluation_timed_out=std::make_shared<bool>(false);
};
struct SeedingInput {
    std::optional<EdgeMask> initial_cut,locked;
    std::optional<std::vector<double>> cost; // computed from the options when empty
    const ChartEvaluator* evaluate=nullptr; // default_evaluator when null and enabled
    Deadline deadline;
    NumericalPolicy policy=NumericalPolicy::scipy_compatible;
    // User pins. PinCompletion::preserve (PIN-03) refuses every generation move that
    // would weld corners pinned to different coordinates; legacy_010 ignores pins
    // during generation as the source does.
    PinMap pins;
    PinCompletion pin_policy=PinCompletion::preserve;
    LayoutProgress progress; // region growth, merge rounds and absorptions
    // T24: under a deadline the evaluator runs in bounded batches (never one call per
    // chart), checking the deadline before and cancellation before and after each batch.
    // Without a deadline each evaluator call is one call, as in the source.
    SchedulingOptions scheduling;
    CancellationToken cancel;
    std::shared_ptr<BatchLog> batches; // records the evaluation batches (optional)
    PinMap effective_pins() const {return pin_policy==PinCompletion::preserve?pins:PinMap{};}
};
CLEANSEAMS_CORE_EXPORT SeededLayout seed_layout(const std::shared_ptr<const MeshSnapshot>&,const LayoutOptions&,const SeedingInput& input={},
    LayoutTrace* trace=nullptr);
// ------------------------------------------------------------ merging and absorption (T13)
// Source absorb.tiny_charts: charts below the face/area floor or thin strips,
// except a chart that is a whole connected component.
struct TinyCharts {std::vector<std::uint8_t> tiny,whole;std::vector<double> areas;std::vector<std::uint64_t> counts;};
CLEANSEAMS_CORE_EXPORT TinyCharts tiny_charts(const MeshSnapshot&,const CutSurface&,const LayoutOptions&);
// The minimum-island predicate of tiny_charts for one face set (a chart, or a region a
// proposal would make one; GB m11): fewer faces than min_island_faces unless they are large
// planar n-gons (small_face_count), an area below `area_floor` (min_island_fraction of the
// mesh area), or a thin strip (fewer than 50 faces with perimeter^2 / (4 pi area) above
// strip_isoperimetric). The whole-component exception is the caller's.
CLEANSEAMS_CORE_EXPORT bool below_island_policy(const MeshSnapshot&,std::span<const MeshIndex> faces,double area,double perimeter,double area_floor,
                                                const LayoutOptions&);
enum class MoveKind : std::uint8_t {merge=0,single_arc=1,self=2,spur=3}; // source _MOVE_RANK order
CLEANSEAMS_CORE_EXPORT std::string_view move_kind_name(MoveKind);
struct MergeMove {MoveKind kind=MoveKind::merge;std::vector<MeshIndex> edges;};
// Source _is_strip and _moves_for_pair (candidate moves before any gate; counts
// strip refusals in the report).
CLEANSEAMS_CORE_EXPORT bool is_strip(const CutSurface&,MeshIndex chart,const IndexRows& arcs,std::int64_t max_width);
CLEANSEAMS_CORE_EXPORT std::vector<MergeMove> moves_for_pair(const CutSurface&,MeshIndex a,MeshIndex b,const IndexRows& arcs,const LayoutOptions&,LayoutReport&);
struct MergeConstraints {
    std::optional<EdgeMask> locked,features; // locked edges are never removed; mostly-feature moves are refused
    std::set<std::array<MeshIndex,2>> frozen_pairs; // (min face of a, min face of b), never merged unless tiny
    CostField cost;
    Deadline deadline;
    PinMap pins; // PIN-03 guard (empty: none)
    LayoutProgress progress; // one event per round
};
// Source merge_patches: greedy rounds over arc moves with topology, proxy,
// compactness, isoperimetric and measured distortion gates. A null evaluator uses
// default_evaluator when LayoutOptions::evaluate_true_distortion is set.
CLEANSEAMS_CORE_EXPORT DiskLayout merge_patches(const std::shared_ptr<const MeshSnapshot>&,std::span<const std::uint8_t> cut,const LayoutOptions&,
    LayoutReport&,const ChartEvaluator* evaluate=nullptr,const MergeConstraints& constraints={},LayoutTrace* trace=nullptr);
struct AbsorbRequest {
    std::optional<std::vector<std::uint8_t>> exempt; // faces already exempt from measurement
    std::optional<EdgeMask> locked;
    Deadline deadline;
    bool unsupported_only=false;
    PinMap pins; // PIN-03 guard (empty: none)
    LayoutProgress progress; // one event per absorption
};
// Source absorb_tiny: glue tiny (or unsupported) charts along their shared arcs,
// judged by the strict/2x/3x/validity-only/exempt ladder.
CLEANSEAMS_CORE_EXPORT DiskLayout absorb_tiny(const std::shared_ptr<const MeshSnapshot>&,std::shared_ptr<const CutSurface> cs,const LayoutOptions&,
    LayoutReport&,const ChartEvaluator* evaluate,const AbsorbRequest& request={},LayoutTrace* trace=nullptr);
struct SpurCleanup {EdgeMask seams;std::vector<MeshIndex> removed;};
// Source cleanup_spurs: dead-end chains removed edge by edge, purely topological.
CLEANSEAMS_CORE_EXPORT SpurCleanup cleanup_spurs(const MeshSnapshot&,std::span<const std::uint8_t> seams,std::optional<std::span<const std::uint8_t>> locked=std::nullopt);
// Source _final_spur_sweep: remove dead-end chains, validated per chart by the evaluator.
CLEANSEAMS_CORE_EXPORT DiskLayout final_spur_sweep(const std::shared_ptr<const MeshSnapshot>&,std::span<const std::uint8_t> cut,
    const std::shared_ptr<const CutSurface>& cs,const LayoutOptions&,LayoutReport&,const ChartEvaluator* evaluate,
    std::optional<std::span<const std::uint8_t>> locked=std::nullopt,LayoutTrace* trace=nullptr,const PinMap& pins={});
// Source assemble_uv: cached validated maps, the evaluator once more for the rest.
CLEANSEAMS_CORE_EXPORT std::vector<Vec2> assemble_uv(const MeshSnapshot&,const std::shared_ptr<const CutSurface>&,LayoutReport&,const ChartEvaluator& evaluate,
    LayoutTrace* trace=nullptr);
// ------------------------------------------------------------ complete generation
struct GeneratedLayout {
    EdgeMask seams; // cut edges without mesh boundary edges
    std::shared_ptr<const CutSurface> cs;
    LayoutReport report;
};
// Source generate_seams: seed_layout, unsupported-chart absorption, merging, tiny
// absorption, the final spur sweep and the assembled map.
CLEANSEAMS_CORE_EXPORT GeneratedLayout generate_seams(const std::shared_ptr<const MeshSnapshot>&,const LayoutOptions&,const SeedingInput& input={},
    LayoutTrace* trace=nullptr);
// Source seam_cost_field: the legacy per-edge cost, or empty when exposure is off.
CLEANSEAMS_CORE_EXPORT std::vector<double> seam_cost_field(const MeshSnapshot&,const LayoutOptions&);
}
