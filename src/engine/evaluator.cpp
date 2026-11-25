// SPDX-License-Identifier: GPL-3.0-or-later
#include "cleanseams/evaluator.hpp"
#include "cleanseams/cache.hpp"
#include "cleanseams/errors.hpp"
#include <cmath>
#include <string>
namespace cleanseams {
bool within_budget(const QualityReport& r,std::uint64_t crossings,const LayoutOptions& o,double relax) {
    if(r.flipped || r.nonfinite || crossings || r.collapsed_uv)return false;
    for(const auto [value,bound]:{std::pair{r.anisotropy_p95,o.max_anisotropy_p95},std::pair{r.area_p95,o.max_area_p95},
                                 std::pair{r.anisotropy_p99,o.max_anisotropy_p99}})
        if(!std::isfinite(value) || !(value<=1+(bound-1)*relax))return false;
    return std::isfinite(r.anisotropy_max) && r.anisotropy_max<=o.max_anisotropy_max;
}
bool measurement_matches(const CutSurface& surface,MeshIndex chart,const ChartMeasurement& m) {
    const auto& ct=surface.tables();
    if(chart>=ct.n_charts || m.corner_ids.size()!=m.corner_uv.size())return false;
    for(auto c:m.corner_ids)if(c>=ct.chart_of_corner.size() || ct.chart_of_corner[c]!=chart)return false;
    return true;
}
ChartEvaluator default_evaluator(std::shared_ptr<const MeshSnapshot> mesh,const LayoutOptions& options,const EvaluatorSettings& settings) {
    if(!mesh)throw ContractError(ErrorCode::invalid_input,"evaluator.mesh","missing mesh snapshot");
    // Source flatten_chart: lscm, harmonic, and ARAP for any other name; a negative ARAP
    // iteration count runs no iteration (range() is empty).
    ChartMethod method=ChartMethod::arap;
    if(options.solver=="lscm")method=ChartMethod::lscm;
    else if(options.solver=="harmonic")method=ChartMethod::harmonic;
    const auto iterations=settings.iterations.value_or(static_cast<std::uint64_t>(std::max<std::int64_t>(0,options.arap_iterations)));
    return [mesh,options,settings,method,iterations](const std::shared_ptr<const CutSurface>& cut,std::span<const MeshIndex> charts,
                                                     const std::map<MeshIndex,double>& relax) {
        if(!cut || !cut->mesh().same_geometry(*mesh))throw ContractError(ErrorCode::invalid_input,"evaluator.cut","cut surface belongs to another mesh");
        const auto geometry=quality_geometry(cut,settings.memory_budget);
        // Each chart is independent: its own factor and workspace, the shared immutable
        // surface and geometry. Results land in their own slots.
        std::vector<std::optional<ChartMeasurement>> slots(charts.size());
        parallel_for(charts.size(),[&](std::size_t i) {
            const auto ch=charts[i];
            ChartSolveRequest request;request.geometry=geometry;request.chart=ch;request.policy=settings.policy;
            request.pin_completion=settings.pin_completion;request.iterations=iterations;request.validate=false;
            request.memory_budget=settings.memory_budget;
            auto solved=flatten_chart(request,method);
            const std::array<MeshIndex,1> selected{ch};const CornerUvView view{solved.corner_uv,solved.corner_ids};
            ChartMeasurement m;m.report=QualitySnapshot(geometry,view,{},std::span<const MeshIndex>(selected)).report();
            m.crossings=boundary_crossings(*cut,view,{},std::span<const MeshIndex>(selected)).counts[ch];
            const auto found=relax.find(ch);
            m.within_budget=within_budget(m.report,m.crossings,options,found==relax.end()?1.0:found->second);
            m.method=std::move(solved.method);m.corner_ids=std::move(solved.corner_ids);m.corner_uv=std::move(solved.corner_uv);
            m.numerics=std::move(solved.numerics);slots[i]=std::move(m);
        },&settings.cancel);
        throw_if_cancelled(settings.cancel);
        // Commit in the requested order after revalidating that every result belongs to
        // its chart of this snapshot.
        ChartMeasurements out;
        for(std::size_t i=0;i<charts.size();++i) {
            if(!slots[i] || !measurement_matches(*cut,charts[i],*slots[i]))
                throw ContractError(ErrorCode::stale_state,"evaluator.commit","measurement does not belong to chart "+std::to_string(charts[i]));
            out[charts[i]]=std::move(*slots[i]);
        }
        return out;
    };
}
// ------------------------------------------------------------ bounded batches
std::uint64_t BatchLog::begin_call() {std::lock_guard lock(mutex_);return ++totals_.calls;}
void BatchLog::add(BatchRecord r) {
    std::lock_guard lock(mutex_);auto& t=totals_;
    const auto n=static_cast<std::uint64_t>(r.charts.size());
    ++t.batches;
    if(r.evaluated) {++t.evaluated_batches;t.charts_evaluated+=n;} else {++t.skipped_batches;t.charts_not_evaluated+=n;}
    t.singleton_batches+=n==1;t.peak_batch_charts=std::max(t.peak_batch_charts,n);t.peak_batch_triangles=std::max(t.peak_batch_triangles,r.triangles);
    if(records_.size()<max_records_)records_.push_back(std::move(r));else ++t.records_dropped;
}
std::vector<BatchRecord> BatchLog::records() const {std::lock_guard lock(mutex_);return records_;}
BatchLog::Totals BatchLog::totals() const {std::lock_guard lock(mutex_);return totals_;}
nlohmann::json BatchLog::json() const {
    std::lock_guard lock(mutex_);const auto& t=totals_;
    nlohmann::json rows=nlohmann::json::array();
    for(const auto& r:records_)rows.push_back({{"stage",r.stage},{"call",r.call},{"charts",r.charts},{"triangles",r.triangles},
                                              {"evaluated",r.evaluated},{"stop",r.stop.empty()?nlohmann::json(nullptr):nlohmann::json(r.stop)}});
    return {{"calls",t.calls},{"batches",t.batches},{"evaluated_batches",t.evaluated_batches},{"skipped_batches",t.skipped_batches},
            {"charts_evaluated",t.charts_evaluated},{"charts_not_evaluated",t.charts_not_evaluated},{"singleton_batches",t.singleton_batches},
            {"peak_batch_charts",t.peak_batch_charts},{"peak_batch_triangles",t.peak_batch_triangles},{"records_dropped",t.records_dropped},
            {"records",std::move(rows)}};
}
std::vector<std::uint64_t> chart_triangle_counts(const CutSurface& cs) {
    const auto& ct=cs.tables();std::vector<std::uint64_t> out(ct.n_charts,0);
    for(auto f:cs.mesh().triangulation().faces)++out[ct.chart_of_face[f]];
    return out;
}
std::vector<std::vector<MeshIndex>> plan_batches(std::span<const MeshIndex> charts,std::span<const std::uint64_t> triangles,
                                                 const SchedulingOptions& options) {
    const auto o=resolve_scheduling(options);
    std::vector<std::vector<MeshIndex>> out;
    std::uint64_t load=0;
    for(auto ch:charts) {
        if(ch>=triangles.size())throw ContractError(ErrorCode::invalid_input,"batches.charts","chart index out of range");
        const auto size=triangles[ch];
        if(out.empty() || out.back().size()>=o.batch_max_charts || load+size>o.batch_max_triangles) {out.emplace_back();load=0;}
        out.back().push_back(ch);load+=size;
    }
    return out;
}
std::vector<std::vector<MeshIndex>> plan_batches(const CutSurface& cs,std::span<const MeshIndex> charts,const SchedulingOptions& options) {
    if(charts.empty()) {resolve_scheduling(options);return {};}
    return plan_batches(charts,chart_triangle_counts(cs),options);
}
ChartEvaluator batched_evaluator(ChartEvaluator raw,BatchedEvaluation settings) {
    if(!raw)throw ContractError(ErrorCode::invalid_input,"evaluator","missing evaluator");
    settings.scheduling=resolve_scheduling(settings.scheduling);
    return [raw=std::move(raw),settings](const std::shared_ptr<const CutSurface>& surface,std::span<const MeshIndex> charts,
                                         const std::map<MeshIndex,double>& relax) {
        ChartMeasurements out;
        if(charts.empty())return out;
        if(!surface)throw ContractError(ErrorCode::invalid_input,"evaluator.cut","missing cut surface");
        const auto triangles=chart_triangle_counts(*surface);
        const auto batches=plan_batches(charts,triangles,settings.scheduling);
        const auto call=settings.log?settings.log->begin_call():0;
        auto record=[&](const std::vector<MeshIndex>& batch,bool evaluated,const char* stop) {
            if(!settings.log)return;
            BatchRecord r;r.stage=settings.stage;r.call=call;r.charts=batch;r.evaluated=evaluated;r.stop=stop;
            for(auto ch:batch)r.triangles+=triangles[ch];
            settings.log->add(std::move(r));
        };
        for(std::size_t b=0;b<batches.size();++b) {
            const auto& batch=batches[b];
            if(settings.cancel.requested()) {
                for(auto k=b;k<batches.size();++k)record(batches[k],false,"cancelled");
                throw_if_cancelled(settings.cancel);
            }
            if(settings.deadline.expired()) {
                // Partial completion: the remaining charts are not evaluated (callers keep
                // their incumbent) and the call is reported as timed out.
                if(settings.timed_out)*settings.timed_out=true;
                for(auto k=b;k<batches.size();++k)record(batches[k],false,"deadline");
                break;
            }
            auto results=raw(surface,batch,relax);
            record(batch,true,"");
            if(settings.cancel.requested()) { // after the active solve: discard, never commit
                for(auto k=b+1;k<batches.size();++k)record(batches[k],false,"cancelled");
                throw_if_cancelled(settings.cancel);
            }
            for(auto& [ch,measured]:results)out[ch]=std::move(measured);
        }
        return out;
    };
}
}
