// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
#include "cleanseams/chart_solver.hpp"
#include "cleanseams/control.hpp"
#include <functional>
#include <map>
#include <mutex>
namespace cleanseams {
// One chart flattened and measured by the source's default_evaluator: the
// distortion report of that chart alone, its boundary crossings and the merge
// gate verdict. The map it was measured on is retained.
struct ChartMeasurement {
    bool within_budget=false;
    QualityReport report;
    std::uint64_t crossings=0;
    std::string method;
    std::vector<MeshIndex> corner_ids;
    std::vector<Vec2> corner_uv;
    FactorMetadata numerics;
};
using ChartMeasurements=std::map<MeshIndex,ChartMeasurement>;
// Evaluate charts of a cut surface; `relax` loosens the budget per chart (source relax dict).
using ChartEvaluator=std::function<ChartMeasurements(const std::shared_ptr<const CutSurface>&,std::span<const MeshIndex> charts,
                                                     const std::map<MeshIndex,double>& relax)>;
// Source quad_layout.within_budget: validity counts are zero and P95/P99/max
// anisotropy and P95 area stay within the (optionally relaxed) layout budget.
CLEANSEAMS_CORE_EXPORT bool within_budget(const QualityReport&,std::uint64_t crossings,const LayoutOptions&,double relax=1.0);
struct EvaluatorSettings {
    std::optional<std::uint64_t> iterations; // default: LayoutOptions::arap_iterations
    NumericalPolicy policy=NumericalPolicy::scipy_compatible;
    PinCompletion pin_completion=PinCompletion::preserve;
    std::uint64_t memory_budget=active_memory_budget();
    // Checked after every chart solve: charts not yet started are skipped and the call
    // throws ContractError(cancelled); nothing of a cancelled call is committed.
    CancellationToken cancel{};
};
// Source quad_layout.default_evaluator: flatten_chart with the layout's solver and
// iteration count, no pins, then the chart's report, crossings and budget. The charts
// of one call are independent: with an active WorkerPool they are measured in parallel
// against the same immutable surface and committed in the requested order.
CLEANSEAMS_CORE_EXPORT ChartEvaluator default_evaluator(std::shared_ptr<const MeshSnapshot>,const LayoutOptions&,const EvaluatorSettings& settings={});
// Commit-time revalidation (T24): a measurement may be committed for `chart` only when
// every corner it maps belongs to that chart of `surface`.
CLEANSEAMS_CORE_EXPORT bool measurement_matches(const CutSurface& surface,MeshIndex chart,const ChartMeasurement&);
// ------------------------------------------------------------ bounded batches (T24)
// One evaluation or host-solve batch as it was actually scheduled.
struct BatchRecord {
    std::string stage;            // "generate" (layout evaluator) or "unwrap" (deadline-budgeted unwrap)
    std::uint64_t call=0;         // the evaluator call (or unwrap pass) the batch belongs to
    std::vector<MeshIndex> charts;
    std::uint64_t triangles=0;
    bool evaluated=false;         // false: skipped, `stop` says why
    std::string stop;             // "", "deadline" or "cancelled"
};
// Thread-safe record of a run's batches. Every batch is counted; the first
// `max_records` are kept in full with their chart lists.
class CLEANSEAMS_CORE_EXPORT BatchLog {
public:
    explicit BatchLog(std::size_t max_records=1024):max_records_(max_records) {}
    std::uint64_t begin_call();
    void add(BatchRecord);
    std::vector<BatchRecord> records() const;
    struct Totals {
        std::uint64_t calls=0,batches=0,evaluated_batches=0,skipped_batches=0,charts_evaluated=0,charts_not_evaluated=0;
        std::uint64_t singleton_batches=0,peak_batch_charts=0,peak_batch_triangles=0,records_dropped=0;
    };
    Totals totals() const;
    nlohmann::json json() const;
private:
    mutable std::mutex mutex_;
    std::size_t max_records_;
    std::vector<BatchRecord> records_;
    Totals totals_;
};
// Triangles of each chart of `surface` (the fixed triangulation the solvers use).
CLEANSEAMS_CORE_EXPORT std::vector<std::uint64_t> chart_triangle_counts(const CutSurface&);
// Consecutive batches of `charts` in the given order, each holding at most
// batch_max_charts charts and batch_max_triangles triangles; a chart above the triangle
// limit is a batch alone. The plan depends only on the charts and the limits, never on
// the thread count, so serial and parallel runs evaluate the same batches.
CLEANSEAMS_CORE_EXPORT std::vector<std::vector<MeshIndex>> plan_batches(const CutSurface&,std::span<const MeshIndex> charts,const SchedulingOptions&);
// The same plan from precomputed per-chart triangle counts (indexed by chart).
CLEANSEAMS_CORE_EXPORT std::vector<std::vector<MeshIndex>> plan_batches(std::span<const MeshIndex> charts,std::span<const std::uint64_t> triangles,
    const SchedulingOptions&);
struct BatchedEvaluation {
    SchedulingOptions scheduling;
    Deadline deadline;
    CancellationToken cancel;
    std::shared_ptr<BatchLog> log;       // optional
    std::shared_ptr<bool> timed_out;     // set when the deadline stops a call early (optional)
    std::string stage="generate";
};
// Replaces the source's deadline wrapper, which made one evaluator call per chart. The
// charts of a call are evaluated in bounded batches, in order. Cancellation is checked
// before every batch and after every active solve (ContractError(cancelled)); the deadline
// is checked before every batch. Charts of batches the deadline stops are left out of the
// result, so callers keep their validated incumbent for them and the log says so.
CLEANSEAMS_CORE_EXPORT ChartEvaluator batched_evaluator(ChartEvaluator raw,BatchedEvaluation settings);
}
