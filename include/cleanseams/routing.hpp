// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
#include "cleanseams/evaluator.hpp"
#include "cleanseams/seam_graph.hpp"
namespace cleanseams {
// core/seam_paths.py. Search states are directed incident edges (prev -> cur),
// the priority tuple is (cost, prev, cur) exactly as the source heap.
struct Corridor {
    std::vector<MeshIndex> verts; // CPython set iteration order of the source
    std::vector<std::uint8_t> vertex_mask,edge_mask;
    std::vector<MeshIndex> edges; // CPython set iteration order (fixes lref's summation)
    std::vector<std::int64_t> dist; // -1 outside the corridor
    double lref=0;
    std::int64_t rings=0;
};
CLEANSEAMS_CORE_EXPORT Corridor build_corridor(const CutSurface&,const Chain&,std::int64_t rings,
    std::optional<std::span<const std::uint8_t>> protected_edges=std::nullopt);
struct AnchorTangent {std::optional<Vec3> direction;double weight=0;};
CLEANSEAMS_CORE_EXPORT AnchorTangent anchor_tangent(const CutSurface&,const Chain&,int end);
CLEANSEAMS_CORE_EXPORT double turn_cost(const Vec3& in,const Vec3& out);
using RouteWeights=Vec3; // (w_turn, w_flow, w_drift)
struct RoutingStats {
    std::uint64_t searches=0,pops=0,pushes=0,early_stops=0,tangent_lookups=0,tangent_computed=0,stamp_resets=0;
    bool nonnegative_costs=true; // asserted before the safe early stop is used
};
// Reusable per-thread buffers: best cost/predecessor per directed edge with
// generation stamps, and directed tangent vectors cached per (vertex, edge).
class CLEANSEAMS_CORE_EXPORT RoutingWorkspace {
public:
    const RoutingStats& stats() const {return stats_;}
    // Test hook: start the stamp just below overflow to exercise the reset.
    void set_generation_for_test(std::uint32_t generation) {generation_=generation;}
private:
    friend struct RoutingAccess;
    const MeshSnapshot* mesh_=nullptr;
    std::string mesh_key_;
    std::vector<std::uint32_t> stamp_,tangent_stamp_;
    std::vector<double> best_;
    std::vector<std::int64_t> pred_;
    std::vector<Vec3> tangent_;
    std::uint32_t generation_=0;
    RoutingStats stats_;
};
struct RouteSearch {
    std::optional<std::vector<MeshIndex>> verts;
    double cost=0;
    std::string reason; // no_route or route_not_simple
};
CLEANSEAMS_CORE_EXPORT RouteSearch find_path(const CutSurface&,const Chain&,const Corridor&,const RouteWeights&,const PathOptions&,
    std::span<const double> exposure,RoutingWorkspace&);
CLEANSEAMS_CORE_EXPORT double route_cost(const CutSurface&,const Chain&,std::span<const MeshIndex> verts,const Corridor&,
    const RouteWeights&,const PathOptions&,std::span<const double> exposure,RoutingWorkspace* workspace=nullptr);
CLEANSEAMS_CORE_EXPORT bool route_covers_chain(const MeshSnapshot&,const Corridor&,std::span<const MeshIndex> route,const Chain&,std::int64_t radius);
CLEANSEAMS_CORE_EXPORT std::optional<std::vector<std::uint8_t>> apply_route(const MeshSnapshot&,std::span<const std::uint8_t> seams,
    const Chain&,std::span<const MeshIndex> verts);
CLEANSEAMS_CORE_EXPORT std::vector<std::uint8_t> corridor_band_faces(const MeshSnapshot&,const Corridor&);
struct TopologyCheck {bool ok=false;std::string reason;std::shared_ptr<const CutSurface> after;};
CLEANSEAMS_CORE_EXPORT TopologyCheck topology_preserved(const std::shared_ptr<const CutSurface>& before,std::span<const std::uint8_t> seams_after,
    const Chain&,std::optional<std::span<const std::uint8_t>> band=std::nullopt);
// Before-reports of flanking charts keyed by chart key digest (source `previous`).
using FlankingCache=std::map<std::string,ChartMeasurement>;
struct FlankingCheck {bool ok=false;std::string reason;ChartMeasurements after;};
CLEANSEAMS_CORE_EXPORT FlankingCheck flanking_ok(const std::shared_ptr<const CutSurface>& before,const std::shared_ptr<const CutSurface>& after,
    const Chain&,const ChartEvaluator*,FlankingCache& previous,bool require_budget=false);
struct RouteAttempt {
    std::int64_t rings=0;
    RouteWeights weights{};
    std::string outcome; // failure reason or skip label
};
struct CLEANSEAMS_CORE_EXPORT RefineResult {
    MeshIndex chain=0;
    std::string status="kept",reason;
    std::vector<MeshIndex> old_route;
    std::optional<std::vector<MeshIndex>> new_route;
    std::optional<RouteWeights> variant;
    std::optional<double> cost_old,cost_new;
    ChainTurns turns_old;
    std::optional<ChainTurns> turns_new;
    std::optional<std::int64_t> rings;
    std::optional<std::int64_t> max_rings; // FIX-02: the resolved maximum radius
    std::vector<RouteAttempt> tried;
    std::vector<RefineResult> parts; // cycle arcs, without their attempts
    bool cycle=false; // refine_cycle result: the source's cycle dict shape
    // Accepted state: new cut surface and corridor vertices (source `state`).
    std::shared_ptr<const CutSurface> cut;
    std::vector<MeshIndex> dirty;
    nlohmann::json to_json() const;
};
struct RefineOutcome {RefineResult result;std::vector<std::uint8_t> seams;};
// T23 corrections of the source's routing, enabled by the native_baseline policy;
// all false reproduces the source (legacy_010).
struct RoutingFixes {
    bool weight_variants=false; // FIX-01: variants follow the exposed w_turn/w_flow/w_drift
    bool corridor_bounds=false; // FIX-02: search up to max(corridor_rings, corridor_max_rings)
    bool full_schedule=false;   // FIX-03: a same route at one radius does not stop wider radii
};
// FIX-01: the exposed weights times (1,1,1), (0.5,2,1), (1.5,1,0.5), (1,1,1/6); the
// defaults (4,2,3) give exactly the source's fixed variants.
CLEANSEAMS_CORE_EXPORT std::vector<RouteWeights> weight_variants(const PathOptions&);
struct RefineRequest {
    std::shared_ptr<const CutSurface> cut;
    std::span<const std::uint8_t> seams;
    PathOptions options;
    std::optional<std::span<const std::uint8_t>> protected_edges;
    std::span<const double> exposure; // seam cost field, empty for none
    const ChartEvaluator* evaluate=nullptr; // guard evaluator, null for none
    FlankingCache* previous=nullptr;
    RoutingWorkspace* workspace=nullptr;
    RoutingFixes fixes;
    // GB performance: with an active worker pool, measure the flanking charts of several
    // candidate routes (and the before-charts) in one parallel call, then judge them in
    // schedule order exactly as the serial loop does. Only for a guard evaluator that
    // measures charts independently and returns every requested chart or throws (the
    // native evaluator); needs `previous`. Results equal the serial schedule's.
    bool speculative=false;
};
CLEANSEAMS_CORE_EXPORT RefineOutcome refine_chain(const RefineRequest&,const Chain&);
CLEANSEAMS_CORE_EXPORT RefineOutcome refine_cycle(const RefineRequest&,const Chain&);
}
