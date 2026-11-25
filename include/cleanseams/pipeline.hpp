// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
#include "cleanseams/layout.hpp"
#include "cleanseams/contours.hpp"
#include "cleanseams/options.hpp"
#include "cleanseams/proposals.hpp"
#include "cleanseams/routing.hpp"
#include "cleanseams/seam_graph.hpp"
#include <functional>
namespace cleanseams {
// core/pipeline.py: gate -> baseline -> [generate] -> [refine] -> unwrap -> repair per
// chart -> per-chart budget gate -> [contours] -> validate -> output hygiene. Every
// stage works on copies; the input mesh, seams, map and pins are never modified.
using ProgressSink=std::function<void(const nlohmann::json&)>;
// Source policy label -> pin handling: legacy_010 reproduces the source (PIN-01..03
// defects); native_baseline and enhanced_v1 keep every explicit pin.
CLEANSEAMS_CORE_EXPORT PinCompletion pin_policy(const PipelineOptions&);
// The solver service (source backend protocol): unwrap(cut, charts, pins) -> (nc) map.
// `iterations` and `init` are honored by the core; a host backend may ignore `init`.
class CLEANSEAMS_CORE_EXPORT UnwrapBackend {
public:
    virtual ~UnwrapBackend();
    virtual std::vector<Vec2> unwrap(std::span<const std::uint8_t> cut,std::optional<std::span<const MeshIndex>> charts,const PinMap& pins,
        std::optional<std::uint64_t> iterations=std::nullopt,std::optional<std::span<const Vec2>> init=std::nullopt)=0;
    virtual std::uint64_t iterations() const=0;
    // The source's CoreBackend (the only backend with a guard evaluator and
    // iteration/start-map support in _unwrap_charts).
    virtual bool core() const {return false;}
    // GB: the maps it returns are stored in float32 (the Blender host's UV layer), so a pinned
    // coordinate comes back as exactly its float32 rounding (T25's check_stored_uv rule).
    virtual bool float32_storage() const {return false;}
};
// Source CoreBackend: NativeArapBackend at ten iterations with the layout's evaluators.
class CLEANSEAMS_CORE_EXPORT CoreUnwrapBackend final:public UnwrapBackend {
public:
    // `cancel` is checked after every chart solve of its unwraps and evaluators (T24).
    CoreUnwrapBackend(std::shared_ptr<const MeshSnapshot>,PipelineOptions,CancellationToken cancel={});
    std::vector<Vec2> unwrap(std::span<const std::uint8_t> cut,std::optional<std::span<const MeshIndex>> charts,const PinMap& pins,
        std::optional<std::uint64_t> iterations=std::nullopt,std::optional<std::span<const Vec2>> init=std::nullopt) override;
    std::uint64_t iterations() const override {return backend_.settings().iterations;}
    bool core() const override {return true;}
    void set_options(PipelineOptions options) {options_=std::move(options);}
    ChartEvaluator evaluator() const;       // default_evaluator(mesh, layout)
    ChartEvaluator guard_evaluator() const; // the final unwrap's iteration count
    InteriorSolver interior_solver() const; // ARAP with boundary.arap_iterations
private:
    std::shared_ptr<const MeshSnapshot> mesh_;
    PipelineOptions options_;
    CancellationToken cancel_;
    NativeArapBackend backend_;
};
// Source backend_evaluator / backend_interior_solver for a host backend.
CLEANSEAMS_CORE_EXPORT ChartEvaluator backend_evaluator(std::shared_ptr<const MeshSnapshot>,UnwrapBackend&,const PipelineOptions&);
CLEANSEAMS_CORE_EXPORT InteriorSolver backend_interior_solver(std::shared_ptr<const MeshSnapshot>,UnwrapBackend&);
struct GateResult {bool ok=true;std::string reason;nlohmann::json info;};
CLEANSEAMS_CORE_EXPORT GateResult gate(const MeshSnapshot&,const PipelineOptions&);
CLEANSEAMS_CORE_EXPORT std::vector<std::uint8_t> passthrough_faces(const MeshSnapshot&,const PipelineOptions&);
// Source _metrics: distortion report with exemptions, summed crossings, validity,
// chart and disk counts.
CLEANSEAMS_CORE_EXPORT nlohmann::json pipeline_metrics(const std::shared_ptr<const CutSurface>&,std::span<const Vec2> uv,
    std::span<const std::uint8_t> exempt);
// ------------------------------------------------------------------ repair ladder
struct RepairContext {
    std::shared_ptr<const MeshSnapshot> mesh;
    UnwrapBackend* backend=nullptr;
    const PipelineOptions* options=nullptr;
    const std::map<std::string,CachedChartMap>* uv_cache=nullptr; // validated generation maps
    const PinMap* pins=nullptr;
    std::uint64_t iterations=10;
    PinCompletion pin_completion=PinCompletion::preserve; // the pipeline sets it from PipelineOptions::policy
    // T24: the deadline-budgeted unwrap solves bounded batches of charts, checking the
    // deadline before and cancellation before and after each batch.
    SchedulingOptions scheduling;
    CancellationToken cancel;
    BatchLog* batches=nullptr; // records the unwrap batches (optional)
    // Native recovery services. Empty uses the monotonic clock and built-in LSCM.
    // Overrides still pass exact corner identity, finite/coverage/precision validation.
    const EdgeMask* recovery_locked=nullptr; // preserve these edges at their current cut state
    ClockSource recovery_clock;
    std::function<ChartSolveResult(const ChartSolveRequest&)> recovery_solver;
};
struct RepairState {
    std::vector<Vec2> uv;
    std::vector<std::uint8_t> exempt;
    EdgeMask seams;
    std::shared_ptr<const CutSurface> cs;
    // Quality geometry of `cs`, rebuilt when the surface changes.
    std::shared_ptr<const QualityGeometry> geometry() const;
private:
    mutable std::shared_ptr<const QualityGeometry> geometry_;
};
// Source repair_chart: L1 re-solve, L2 exempt tiny interior offenders, L2b local
// untangle (2- and 3-ring), L5 cached map, fallback split. Returns the record.
CLEANSEAMS_CORE_EXPORT nlohmann::json repair_chart(const RepairContext&,RepairState&,MeshIndex chart,double density=1.0);
struct RepairAllResult {nlohmann::json records=nlohmann::json::array();std::vector<MeshIndex> invalid;};
CLEANSEAMS_CORE_EXPORT RepairAllResult repair_all(const RepairContext&,RepairState&,const std::function<void(const std::string&)>& say={},
    std::size_t max_repairs=64);
// Source final_chart_gate: one split pass for charts over the P95 budget.
CLEANSEAMS_CORE_EXPORT nlohmann::json final_chart_gate(const RepairContext&,RepairState&,const std::function<void(const std::string&)>& say={});
// FIX-04 (policies other than legacy_010): [chart, area P95 rounded to 3 digits] of the
// charts (every chart by default) with at least layout.split_min_faces faces whose area
// P95, exempt faces left out, exceeds layout.max_area_p95 x final_gate_tolerance.
CLEANSEAMS_CORE_EXPORT nlohmann::json area_over_budget(const RepairState&,const PipelineOptions&,
    std::optional<std::span<const MeshIndex>> charts=std::nullopt);
// FIX-04: the chart_over_budget rows of area-only charts, [chart, anisotropy P95 rounded
// to 3 digits or null], in the given order.
CLEANSEAMS_CORE_EXPORT nlohmann::json area_only_rows(std::span<const MeshIndex> charts,std::span<const double> anisotropy_p95);
// Source _budget_unwrap: cached charts, solves, conformal/harmonic fallbacks and
// explicitly reported emergency face projections. Returns projected faces. T24: the
// charts to solve go in bounded batches (the source solved one chart per call); a batch
// that fails is retried chart by chart, so a failure stays local to its chart.
CLEANSEAMS_CORE_EXPORT std::uint64_t budget_unwrap(const RepairContext&,RepairState&,const std::optional<std::vector<Vec2>>& cached,
    const Deadline&,nlohmann::json& fallback_info);
// Source absorb_repair_fragments: repaired specks pass generation's absorption rules.
CLEANSEAMS_CORE_EXPORT nlohmann::json absorb_repair_fragments(const RepairContext&,RepairState&,const ChartEvaluator* evaluate,
    std::span<const std::uint8_t> locked,const Deadline&);
// Source local_untangle (core path): pin the chart outside the rings of `faces` and
// re-solve the ring with mean value weights. Empty when no sector is free.
CLEANSEAMS_CORE_EXPORT std::optional<std::vector<Vec2>> local_untangle(const RepairContext&,const std::shared_ptr<const QualityGeometry>&,
    MeshIndex chart,std::span<const Vec2> uv,std::span<const MeshIndex> faces,std::int64_t rings=2);
CLEANSEAMS_CORE_EXPORT std::vector<Vec2> planar_projection(const CutSurface&,std::span<const MeshIndex> charts,std::span<const Vec2> uv,double density);
// ------------------------------------------------------------------ orchestration
struct PipelineInput {
    std::shared_ptr<const MeshSnapshot> mesh;
    EdgeMask seams;                         // bool over edges
    std::optional<std::vector<Vec2>> uv;    // incoming map (baseline), per corner
    PinMap pins;                            // user pins {corner: (u, v)}
    std::optional<EdgeMask> locked;         // edges that must not move
    std::vector<TransferIndex> anchors;     // manual anchor vertices
    PipelineOptions options;
    // Topology/geometry cache for the run when the caller has none active; 0 disables.
    // Default: the lesser of 512 MiB and a quarter of the 2 GiB session budget.
    std::uint64_t cache_bytes=536870912ULL;
    // The run's memory allowance (Gate A F3; a session passes SessionConfig's). Checked
    // up front against the tables the run builds (memory_limit), then the allowance of
    // every table, solve and exposure field of the run (active_memory_budget()).
    std::uint64_t memory_budget=default_memory_budget;
    // T24: native evaluation threads (the run owns its WorkerPool unless the caller has
    // one active) and the batch bounds of deadline work. Results never depend on them.
    SchedulingOptions scheduling;
    // enhanced_v1 only (T28): proposal generators evaluated after the built-in ones
    // (enhanced_generators), through the same transaction and common ordering. Every
    // other policy ignores them.
    std::vector<ProposalGenerator> proposal_generators;
    // enhanced_v1 only (T29): painted hard constraints and soft preferences (typed arrays in
    // core order and raw host attributes), resolved natively. Every other policy ignores them.
    ConstraintInput constraints;
    std::string operation="pipeline";
    std::vector<std::uint8_t> selected_faces;
};
struct PipelineResult {
    bool accepted=false;
    std::string reason;
    EdgeMask seams;                         // without mesh boundary edges
    std::vector<Vec2> uv;
    std::shared_ptr<const CutSurface> cut;
    nlohmann::json report=nlohmann::json::object();
    std::vector<MeshIndex> invalid_charts;
    std::vector<std::uint8_t> exempt;       // faces excluded from validation
    // T24 scheduling record, kept out of the legacy report: threads and batch bounds,
    // every batch actually scheduled (stage, charts, triangles, evaluated or why not)
    // and the worker pool's counters.
    nlohmann::json scheduling=nlohmann::json::object();
    std::vector<MeshIndex> changed_faces;  // original face IDs, supplied-layout operations
};
// Source pipeline.run: one synchronous run with wall/CPU timings. A null backend uses
// CoreUnwrapBackend. Cancellation throws ContractError(cancelled) from a progress
// checkpoint and leaves every input untouched (no partial result is applied).
CLEANSEAMS_CORE_EXPORT PipelineResult run_pipeline(const PipelineInput&,UnwrapBackend* backend=nullptr,const ProgressSink& progress={},
    const Control& control={});
}
