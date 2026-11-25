// SPDX-License-Identifier: GPL-3.0-or-later
// Phase B proposal transaction (T28, Section 11): feasibility before scoring, exact trial
// topology, changed-chart solves, strict validity, pins, budgets, the run-origin quality
// envelope and atomic commits. Used by policy enhanced_v1 only.
#include "cleanseams/proposals.hpp"
#include "cleanseams/symmetric_moves.hpp"
#include "cleanseams/physical_scale.hpp"
#include "cleanseams/cache.hpp"
#include "cleanseams/errors.hpp"
#include "cleanseams/detail/numpy_compat.hpp"
#include "cleanseams/detail/blake2b.hpp"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <numeric>
#include <set>
namespace cleanseams {
namespace {
using Json=nlohmann::json;
bool always_cut(const MeshTables& t,std::size_t e) {
    return t.edge_boundary[e] || !t.edge_manifold[e] || !t.edge_consistent[e] || t.edge_corners[e].size()!=2;
}
std::vector<TransferIndex> face_ids(std::span<const std::uint8_t> mask) {
    std::vector<TransferIndex> out;for(std::size_t f=0;f<mask.size();++f)if(mask[f])out.push_back(static_cast<TransferIndex>(f));return out;
}
Json finite_number(double x) {return std::isfinite(x)?Json(x):Json(nullptr);}
std::string generator_name(const Proposal& p) {return p.identity.generator.empty()?"anonymous":p.identity.generator;}
// The chart budget of the enhanced policy: the final chart gate's P95 rule (layout P95
// budgets times final_gate_tolerance on charts of at least split_min_faces non-exempt faces,
// FIX-04 included) plus the layout's P99/maximum sanity caps on every chart.
double budget_ratio(const QualityReport& r,std::uint64_t measured_faces,const PipelineOptions& o) {
    double ratio=0;
    auto take=[&](double value,double limit) {if(std::isfinite(value) && limit>0)ratio=std::max(ratio,value/limit);};
    if(static_cast<std::int64_t>(measured_faces)>=o.layout.split_min_faces) {
        take(r.anisotropy_p95,o.layout.max_anisotropy_p95*o.final_gate_tolerance);
        take(r.area_p95,o.layout.max_area_p95*o.final_gate_tolerance);
    }
    take(r.anisotropy_p99,o.layout.max_anisotropy_p99);
    take(r.anisotropy_max,o.layout.max_anisotropy_max);
    return ratio;
}
struct Measured {
    std::vector<ChartOutcome> outcomes;          // indexed by chart; filled for `selected`
    std::vector<MeshIndex> selected;
    std::shared_ptr<const QualityGeometry> geometry;
    std::optional<QualitySnapshot> snapshot;
};
// Every validity check of the selected charts: the compatibility checks (flips, nonfinite,
// collapse, boundary crossings) under the run's exemptions, the strict overlap check, pins,
// the chart budget and face ancestry.
Measured measure_charts(const ProposalContext& ctx,const std::shared_ptr<const CutSurface>& surface,std::span<const Vec2> uv,
                        std::optional<std::span<const MeshIndex>> charts) {
    const auto& cs=*surface;const auto& ct=cs.tables();const auto& mesh=cs.mesh();const auto& g=mesh.geometry();
    Measured out;out.outcomes.resize(ct.n_charts);out.geometry=quality_geometry(surface);
    if(charts)out.selected.assign(charts->begin(),charts->end());
    else {out.selected.resize(ct.n_charts);std::iota(out.selected.begin(),out.selected.end(),MeshIndex{0});}
    std::sort(out.selected.begin(),out.selected.end());out.selected.erase(std::unique(out.selected.begin(),out.selected.end()),out.selected.end());
    if(out.selected.empty())return out;
    QualityPolicy policy;policy.exempt_faces=face_ids(ctx.exempt);
    const std::span<const MeshIndex> selection(out.selected);
    const QualityEvaluation evaluation(out.geometry,CornerUvView{uv},policy,selection);
    OverlapOptions overlap_options;overlap_options.area_tolerance=ctx.options.enhanced.overlap_area_tolerance;
    overlap_options.ignored_faces=policy.exempt_faces;overlap_options.retain_pairs=false;
    const auto overlap=chart_overlaps(*out.geometry,CornerUvView{uv},overlap_options,selection);
    const auto& validity=evaluation.validity();
    const auto ancestry=ctx.envelope.origin_chart_of_face().empty()?std::vector<std::vector<MeshIndex>>{}:
        face_ancestry(cs,ctx.envelope.origin_chart_of_face());
    const double tolerance=ctx.options.enhanced.pin_tolerance;
    // A pinned coordinate is kept within pin_tolerance or, for a float32-stored host map, when it
    // is exactly the float32 rounding of the pin (a serialization artefact, not a move).
    auto kept=[&](double got,double pin) {return strict_pin_kept(got,pin,tolerance,ctx.float32_storage);};
    for(auto ch:out.selected) {
        auto& o=out.outcomes[ch];o.key=cs.chart_key(ch).digest();const auto faces=ct.chart_faces[ch];o.faces=faces.size();
        std::uint64_t measured_faces=0;
        for(auto f:faces) {o.area+=g.face_area[f];measured_faces+=!(f<ctx.exempt.size() && ctx.exempt[f]) && g.face_degenerate_kind[f]==0;}
        o.exempt_only=measured_faces==0;o.disk=ct.chart_is_disk[ch]!=0;o.measured=true;
        if(!ancestry.empty())o.ancestors=ancestry[ch];
        if(o.exempt_only) {o.valid=o.overlap_free=o.within_budget=true;}
        else {
            o.flipped=validity.flipped[ch];o.nonfinite=validity.nonfinite[ch];o.collapsed=validity.collapsed[ch];o.crossings=validity.crossings[ch];
            o.valid=validity.valid[ch]!=0;o.overlaps=overlap.counts[ch];o.tested=overlap.tested[ch]!=0;o.overlap_free=o.overlaps==0 && o.tested;
            const auto r=evaluation.snapshot().report(ch);
            o.anisotropy_p95=r.anisotropy_p95;o.anisotropy_p99=r.anisotropy_p99;o.anisotropy_max=r.anisotropy_max;o.area_p95=r.area_p95;
            o.budget_ratio=budget_ratio(r,measured_faces,ctx.options);o.within_budget=o.budget_ratio<=1.0;
        }
        for(const auto& [corner,pin]:ctx.constraints.pins) {
            if(corner>=ct.chart_of_corner.size() || ct.chart_of_corner[corner]!=ch)continue;
            const auto& got=uv[corner];
            for(std::size_t k=0;k<2;++k)
                if(!kept(got[k],pin[k]))o.pins_kept=false;
            // GB N2: native_baseline's rule (the pipeline's pins_kept), recorded to recognise a
            // strict-only pin defect of candidate zero; it never admits a new map.
            const std::array<double,2> a{got[0],got[1]},b{pin[0],pin[1]};
            if(!detail::numpy_allclose(a,b))o.pins_native_kept=false;
        }
    }
    out.snapshot.emplace(evaluation.snapshot());
    return out;
}
BudgetSummary summarize_budget(std::span<const ChartOutcome> charts) {
    BudgetSummary b;b.measured=true;
    for(const auto& o:charts) {
        if(o.exempt_only || !o.measured)continue;
        b.worst_ratio=std::max(b.worst_ratio,o.budget_ratio);
        if(!o.within_budget) {b.compliant=false;++b.over_budget_charts;b.total_excess+=o.area*std::max(0.0,o.budget_ratio-1.0);}
    }
    return b;
}
// First failing class of a set of chart outcomes ("" when every chart is valid). GB: an
// inherited strict-only defect of candidate zero (ChartOutcome::strict_defect: overlaps, and GB N2
// a pin only the strict rule refuses) is not a failure here; it is counted by
// strict_defect_charts and keeps the state's `valid` false.
std::pair<std::string,Json> chart_failures(std::span<const ChartOutcome> charts,std::span<const MeshIndex> ids) {
    std::vector<MeshIndex> invalid,overlapping,moved,untestable;
    for(auto ch:ids) {
        const auto& o=charts[ch];if(o.exempt_only)continue;
        if(!o.valid)invalid.push_back(ch);
        if(!o.tested)untestable.push_back(ch);else if(o.overlaps && !o.strict_defect)overlapping.push_back(ch);
        if(!o.pins_kept && !o.strict_defect)moved.push_back(ch);
    }
    if(!moved.empty())return {"pin_moved",{{"charts",moved}}};
    if(!invalid.empty())return {"invalid_uv",{{"charts",invalid}}};
    if(!untestable.empty())return {"invalid_uv",{{"charts",untestable},{"untestable",true}}};
    if(!overlapping.empty())return {"overlap",{{"charts",overlapping}}};
    return {"",Json::object()};
}
// GB: the charts carrying an inherited strict-only defect of candidate zero.
std::vector<MeshIndex> strict_defect_charts(std::span<const ChartOutcome> charts) {
    std::vector<MeshIndex> out;
    for(MeshIndex ch=0;ch<charts.size();++ch)
        if(charts[ch].strict_defect && (charts[ch].overlaps || !charts[ch].pins_kept) && !charts[ch].exempt_only)out.push_back(ch);
    return out;
}
// A strict-only defect: every compatibility check and the overlap test pass, and only a strict
// check fails (Section 11.3: the source's epsilon behaviour stays in the compatibility metrics;
// the robust policy is labelled separately): the strict overlap check finds positive-area
// intrachart overlaps, and/or (GB N2) a pinned UV fails the strict pin rule while
// native_baseline's pin rule keeps it (a host solve moved it by a few float32 steps).
bool strict_only_failure(const ChartOutcome& o) {
    return !o.exempt_only && o.measured && o.valid && o.tested && o.pins_native_kept && (o.overlaps>0 || !o.pins_kept);
}
std::string violation_reason(const Json& v) {
    if(v.contains("locked_edges"))return "locked_edge";
    if(v.contains("required_edges"))return "required_edge";
    if(v.contains("forbidden_edges"))return "forbidden_edge";
    if(v.contains("anchors"))return "fixed_anchor";
    return "constraint";
}
std::vector<MeshIndex> all_charts(const CutSurface& cs) {std::vector<MeshIndex> out(cs.tables().n_charts);std::iota(out.begin(),out.end(),MeshIndex{0});return out;}
}
bool strict_pin_kept(double got,double pin,double tolerance,bool float32_storage) {
    return std::abs(got-pin)<=tolerance*(1.0+std::abs(pin)) || (float32_storage && got==static_cast<double>(static_cast<float>(pin)));
}
std::string cut_digest(std::span<const std::uint8_t> cuts) {
    // BLAKE2b-256 of the effective cut mask bytes.
    detail::Blake2b hash(32);hash.update(cuts);return hash.hex_digest();
}
nlohmann::json run_state_summary(const RunState& s) {
    std::uint64_t disks=0;std::vector<MeshIndex> non_disk;
    for(MeshIndex ch=0;ch<s.charts.size();++ch) {disks+=s.charts[ch].disk;if(!s.charts[ch].disk && !s.charts[ch].exempt_only)non_disk.push_back(ch);}
    // GB N3: `budget` is the enhanced chart budget (enhanced_chart_budget: the P95 anisotropy and
    // area rules on charts of at least split_min_faces measured faces plus the P99 and maximum
    // anisotropy caps on every chart), not native_baseline's chart_over_budget listing rule.
    return {{"revision",s.revision},{"feasible",s.feasible},{"valid",s.valid},{"strict_only",s.strict_only},{"admissible",s.admissible()},
            {"strict_defects",s.strict_defects},{"invalid_reason",s.invalid_reason},{"violations",s.violations},
            {"charts",s.charts.size()},{"disks",disks},{"non_disk_charts",non_disk},{"cut_edges",std::count(s.cuts.begin(),s.cuts.end(),std::uint8_t{1})},
            {"cut_digest",cut_digest(s.cuts)},{"objective",s.objective.to_json()},
            {"budget",{{"rule","enhanced_chart_budget"},{"compliant",s.budget.compliant},{"worst_ratio",s.budget.worst_ratio},{"total_excess",s.budget.total_excess},
                       {"over_budget_charts",s.budget.over_budget_charts}}}};
}
// ------------------------------------------------------------ records
nlohmann::json ChartOutcome::to_json() const {
    return {{"key",key},{"faces",faces},{"area",area},{"exempt_only",exempt_only},{"disk",disk},{"valid",valid},{"overlap_free",overlap_free},
            {"tested",tested},{"pins_kept",pins_kept},{"within_budget",within_budget},{"flipped",flipped},{"nonfinite",nonfinite},
            {"collapsed",collapsed},{"crossings",crossings},{"overlaps",overlaps},{"anisotropy_p95",finite_number(anisotropy_p95)},
            {"anisotropy_p99",finite_number(anisotropy_p99)},{"anisotropy_max",finite_number(anisotropy_max)},{"area_p95",finite_number(area_p95)},
            {"budget_ratio",finite_number(budget_ratio)},{"ancestors",ancestors},{"strict_defect",strict_defect},{"pins_native_kept",pins_native_kept}};
}
nlohmann::json SeamObjective::to_json() const {
    return {{"version",version},{"area",area},{"seam_length",seam_length},{"exposed_length",exposed_length},{"turn_length",turn_length},
            {"fragment_area",fragment_area},{"seam_edges",seam_edges},{"charts",charts},{"components",components},
            {"additional_charts",additional_charts},{"small_fragments",small_fragments},{"exposure_measured",exposure_measured},
            {"L",L},{"V",V},{"T",T},{"F",F},{"J",J}};
}
nlohmann::json EnvelopeCheck::to_json() const {
    // GB N7: an admission without a comparison (empty domain, unmeasured origin, infeasible parent,
    // never reached) is reported as not measured, never as `inside: true`.
    return {{"inside",checked || !inside?Json(inside):Json(nullptr)},{"checked",checked},{"faces",faces},{"area",area},
            {"origin_anisotropy_p95",finite_number(origin_anisotropy_p95)},
            {"trial_anisotropy_p95",finite_number(trial_anisotropy_p95)},{"origin_area_p95",finite_number(origin_area_p95)},
            {"trial_area_p95",finite_number(trial_area_p95)},{"origin_anisotropy_max",finite_number(origin_anisotropy_max)},
            {"trial_anisotropy_max",finite_number(trial_anisotropy_max)},{"reason",reason},{"relief_rebases",relief_rebases}};
}
nlohmann::json SearchCounters::to_json() const {
    return {{"proposals",proposals},{"solves",solves},{"charts_solved",charts_solved},{"duplicates",duplicates},{"commits",commits},
            {"solves_by_chart",solves_by_chart},{"candidates_by_chart",candidates_by_chart},{"by_generator",by_generator}};
}
nlohmann::json EvaluatedProposal::summary() const {
    // GB N7: a check that never ran is null, never its default (a trial refused before its
    // measurement has no validity, budget, envelope or J; one refused before its trial state was
    // built has no feasibility or topology either).
    const bool reached=trial.revision!=0,measured=trial.evaluated_revision!=0;
    const Json none(nullptr);
    return {{"id",proposal.identity.id},{"generator",proposal.identity.generator},{"generator_version",proposal.identity.generator_version},
            {"parent_revision",parent_revision},{"trial_revision",trial.revision},{"evaluated_revision",trial.evaluated_revision},
            {"added_edges",proposal.identity.added_edges},{"removed_edges",proposal.identity.removed_edges},
            {"initialization",proposal.initialization},{"writers",trial.writers},{"admissible",admissible()},{"rejection",rejection},
            {"detail",detail},{"feasible",reached?Json(feasible):none},{"topology_ok",reached && feasible?Json(topology_ok):none},
            {"valid",measured?Json(valid):none},{"strict_defects",measured?Json(strict_defects):none},
            {"inside_envelope",envelope.checked || !envelope.inside?Json(inside_envelope):none},
            {"changed_charts",changed_charts},{"solves",solves},{"objective",objective.version.empty()?none:objective.to_json()},
            {"budget",budget.measured?Json{{"compliant",budget.compliant},{"worst_ratio",budget.worst_ratio},{"total_excess",budget.total_excess}}:none},
            {"envelope",envelope.checked || !envelope.inside || !envelope.reason.empty()?envelope.to_json():none},{"seconds",seconds},
            {"provenance",proposal.provenance}};
}
std::string_view budget_status_name(BudgetStatus s) {
    switch(s) {case BudgetStatus::compliant:return "compliant";case BudgetStatus::improved_over_budget:return "improved_over_budget";default:return "over_budget";}
}
// ------------------------------------------------------------ measurement helpers
std::vector<DistortionSample> distortion_samples(const QualitySnapshot& snapshot) {
    // The distortion report's accumulation: per-chart density from the usable triangles in
    // triangle order, anisotropy s0/s1 and area max(r, 1/r) with the report's reciprocal
    // handling and its 1e9 replacement of nonfinite values.
    const auto& tris=snapshot.triangles();
    auto usable=[](const TriangleQuality& t) {return t.finite && !t.sliver && !t.exempt && std::isfinite(t.singular_values[0]) && std::isfinite(t.singular_values[1]);};
    std::map<MeshIndex,std::pair<double,double>> sums;
    for(const auto& t:tris)if(usable(t)) {auto& s=sums[t.chart];s.first+=t.area3;s.second+=std::abs(t.area2);}
    std::vector<DistortionSample> out;out.reserve(tris.size());
    for(const auto& t:tris) {
        if(!usable(t))continue;
        const auto& s=sums.at(t.chart);const double density=s.first>0?s.second/std::max(s.first,1e-300):1.0;
        double a=t.singular_values[0]/std::max(t.singular_values[1],1e-300);
        const double ratio=t.singular_values[0]*t.singular_values[1]/std::max(density,1e-300);
        double b=std::max(ratio,1/std::max(ratio,1e-300));
        if(!std::isfinite(a))a=1e9;if(!std::isfinite(b))b=1e9;
        out.push_back({t.face,t.area3,a,b});
    }
    return out;
}
BaselineEnvelope::BaselineEnvelope(Revision origin,const QualitySnapshot& snapshot,std::span<const Vec2> uv,double p95_growth,double max_growth):
    origin_(origin),p95_growth_(p95_growth),max_growth_(max_growth),chart_of_face_(snapshot.geometry().cut().tables().chart_of_face),uv_(uv.begin(),uv.end()) {
    const auto& cs=snapshot.geometry().cut();
    for(MeshIndex ch=0;ch<cs.tables().n_charts;++ch)keys_.push_back(cs.chart_key(ch).digest());
    std::sort(keys_.begin(),keys_.end());samples_=distortion_samples(snapshot);
    std::stable_sort(samples_.begin(),samples_.end(),[](const DistortionSample& a,const DistortionSample& b) {return a.face<b.face;});
}
bool BaselineEnvelope::has_key(const std::string& key) const {return std::binary_search(keys_.begin(),keys_.end(),key);}
std::vector<MeshIndex> BaselineEnvelope::domain_charts(const CutSurface& surface,std::span<const Vec2> uv) const {
    // A chart is unchanged since the origin when the origin has its key and it keeps the
    // origin's map on every corner; every other chart is in the domain.
    const auto& ct=surface.tables();std::vector<MeshIndex> out;
    std::vector<std::uint8_t> moved(ct.n_charts,0);
    if(uv.size()==uv_.size())for(std::size_t c=0;c<uv.size();++c)if(uv[c]!=uv_[c])moved[ct.chart_of_corner[c]]=1;
    for(MeshIndex ch=0;ch<ct.n_charts;++ch)
        if(uv.size()!=uv_.size() || moved[ch] || !has_key(surface.chart_key(ch).digest()))out.push_back(ch);
    return out;
}
EnvelopeCheck BaselineEnvelope::check(std::span<const MeshIndex> faces,std::span<const DistortionSample> trial,std::span<const std::size_t> active) const {
    EnvelopeCheck out;out.faces=faces.size();
    if(faces.empty()) {out.reason="empty_domain";return out;}
    auto in_domain=[&](MeshIndex f) {return std::binary_search(faces.begin(),faces.end(),f);};
    std::vector<double> ov,oa,ow,tv,ta,tw;
    // GB m10: the origin samples of faces under an active relief rebase are the chart's before
    // the split; merged in face order like rebase_faces, so the accumulation order is the same.
    std::vector<DistortionSample> merged;const std::vector<DistortionSample>* origin=&samples_;
    if(!active.empty()) {
        std::vector<const Staged*> rebases;
        for(auto k:active) {
            if(k>=staged_.size())throw ContractError(ErrorCode::invalid_input,"enhanced.envelope","unknown staged rebase");
            rebases.push_back(&staged_[k]);
        }
        auto rebased=[&](MeshIndex f) {
            for(const auto* r:rebases)if(std::binary_search(r->faces.begin(),r->faces.end(),f))return true;
            return false;
        };
        for(const auto& s:samples_)if(!rebased(s.face))merged.push_back(s);
        for(const auto* r:rebases)for(const auto& s:r->samples)merged.push_back(s);
        std::stable_sort(merged.begin(),merged.end(),[](const DistortionSample& a,const DistortionSample& b) {return a.face<b.face;});
        origin=&merged;
    }
    for(const auto& s:*origin)if(in_domain(s.face)) {ov.push_back(s.anisotropy);oa.push_back(s.area);ow.push_back(s.weight);out.area+=s.weight;}
    // The trial's samples keep their triangle order (the report's accumulation order).
    for(const auto& s:trial)if(in_domain(s.face)) {tv.push_back(s.anisotropy);ta.push_back(s.area);tw.push_back(s.weight);}
    if(ov.empty()) {out.reason="origin_unmeasured";return out;}
    const WeightedDistribution oan(ov,ow),oar(oa,ow);
    out.origin_anisotropy_p95=oan.percentile(95);out.origin_area_p95=oar.percentile(95);out.origin_anisotropy_max=*std::max_element(ov.begin(),ov.end());
    if(tv.empty()) {out.inside=false;out.reason="trial_unmeasured";return out;}
    const WeightedDistribution tan(tv,tw),tar(ta,tw);
    out.trial_anisotropy_p95=tan.percentile(95);out.trial_area_p95=tar.percentile(95);out.trial_anisotropy_max=*std::max_element(tv.begin(),tv.end());
    // Relative to the run origin only; the absolute budgets are separate (chart budgets).
    out.checked=true;   // GB N7: origin and trial were compared
    const double slack=1.0+1e-12;
    if(!(out.trial_anisotropy_p95<=out.origin_anisotropy_p95*p95_growth_*slack))out.reason="anisotropy_p95";
    else if(!(out.trial_area_p95<=out.origin_area_p95*p95_growth_*slack))out.reason="area_p95";
    else if(!(out.trial_anisotropy_max<=out.origin_anisotropy_max*max_growth_*slack))out.reason="anisotropy_max";
    out.inside=out.reason.empty();
    return out;
}
void BaselineEnvelope::rebase_faces(std::span<const MeshIndex> faces,std::span<const DistortionSample> samples) {
    // T32: the faces' origin is an earlier baseline state; everything else is unchanged.
    const std::set<MeshIndex> domain(faces.begin(),faces.end());
    std::vector<DistortionSample> kept;kept.reserve(samples_.size()+samples.size());
    for(const auto& s:samples_)if(!domain.contains(s.face))kept.push_back(s);
    for(const auto& s:samples)if(domain.contains(s.face))kept.push_back(s);
    std::stable_sort(kept.begin(),kept.end(),[](const DistortionSample& a,const DistortionSample& b) {return a.face<b.face;});
    samples_=std::move(kept);rebased_faces_+=domain.size();
}
std::size_t BaselineEnvelope::stage_rebase(std::span<const MeshIndex> faces,std::span<const DistortionSample> samples) {
    Staged s;s.faces.assign(faces.begin(),faces.end());std::sort(s.faces.begin(),s.faces.end());s.faces.erase(std::unique(s.faces.begin(),s.faces.end()),s.faces.end());
    for(const auto& x:samples)if(std::binary_search(s.faces.begin(),s.faces.end(),x.face))s.samples.push_back(x);
    std::stable_sort(s.samples.begin(),s.samples.end(),[](const DistortionSample& a,const DistortionSample& b) {return a.face<b.face;});
    staged_.push_back(std::move(s));
    return staged_.size()-1;
}
nlohmann::json BaselineEnvelope::to_json() const {
    std::uint64_t staged_faces=0;for(const auto& s:staged_)staged_faces+=s.faces.size();
    return {{"version",quality_envelope_version},{"origin_revision",origin_},{"p95_growth",p95_growth_},{"max_growth",max_growth_},
            {"origin_charts",keys_.size()},{"samples",samples_.size()},{"rebased_faces",rebased_faces_},
            {"relief_rebases",{{"splits",staged_.size()},{"faces",staged_faces},
                               {"applies_to","hotspot_relief relief_before_split trials of that split and states that adopted one (GB m10)"}}}};
}
std::vector<std::size_t> relief_rebases(const RunState& parent,const Proposal& proposal) {
    // The split a relief proposal welds (its provenance, from the relief generator only).
    auto split_of=[](const std::string& generator,const nlohmann::json& provenance) -> std::optional<std::size_t> {
        if(generator!="hotspot_relief" || !provenance.is_object() || provenance.value("kind","")!="relief_before_split")return std::nullopt;
        const auto found=provenance.find("split_index");
        if(found==provenance.end() || !found->is_number_integer() || found->get<std::int64_t>()<0)return std::nullopt;
        return static_cast<std::size_t>(found->get<std::int64_t>());
    };
    std::set<std::size_t> out;
    if(const auto k=split_of(proposal.identity.generator,proposal.provenance))out.insert(*k);
    for(const auto& record:parent.provenance)
        if(record.is_object() && record.contains("provenance"))
            if(const auto k=split_of(record.value("generator",""),record.at("provenance")))out.insert(*k);
    return {out.begin(),out.end()};
}
std::vector<double> sector_area_mass(const CutSurface& cs) {
    const auto& mesh=cs.mesh();const auto& tri=mesh.triangulation();const auto& ct=cs.tables();const auto& pos=mesh.input().positions;
    std::vector<double> mass(ct.n_sectors,0.0);
    for(std::size_t t=0;t<tri.corners.size();++t) {
        std::array<Vec3,3> p{};for(std::size_t k=0;k<3;++k)p[k]=pos[static_cast<std::size_t>(mesh.input().corner_vertices[tri.corners[t][k]])];
        const Vec3 u{p[1][0]-p[0][0],p[1][1]-p[0][1],p[1][2]-p[0][2]},v{p[2][0]-p[0][0],p[2][1]-p[0][1],p[2][2]-p[0][2]};
        const Vec3 n{u[1]*v[2]-u[2]*v[1],u[2]*v[0]-u[0]*v[2],u[0]*v[1]-u[1]*v[0]};
        const double third=0.5*std::sqrt(n[0]*n[0]+n[1]*n[1]+n[2]*n[2])/3.0;
        for(std::size_t k=0;k<3;++k)mass[ct.sector_of_corner[tri.corners[t][k]]]+=third;
    }
    return mass;
}
std::vector<std::vector<MeshIndex>> face_ancestry(const CutSurface& surface,std::span<const MeshIndex> origin_chart_of_face) {
    const auto& ct=surface.tables();
    if(origin_chart_of_face.size()!=ct.chart_of_face.size())throw ContractError(ErrorCode::invalid_input,"ancestry","origin face map does not match the mesh faces");
    std::vector<std::vector<MeshIndex>> out(ct.n_charts);
    for(MeshIndex ch=0;ch<ct.n_charts;++ch) {
        for(auto f:ct.chart_faces[ch])out[ch].push_back(origin_chart_of_face[f]);
        std::sort(out[ch].begin(),out[ch].end());out[ch].erase(std::unique(out[ch].begin(),out[ch].end()),out[ch].end());
    }
    return out;
}
// ------------------------------------------------------------ constraints
HardConstraints hard_constraints(const MeshSnapshot& mesh,const std::optional<EdgeMask>& locked,std::span<const TransferIndex> anchors,const PinMap& pins,
                                 const std::optional<EdgeMask>& required,const std::optional<EdgeMask>& forbidden) {
    const auto ne=mesh.tables().edge_vertices.size(),nv=mesh.input().positions.size(),nc=mesh.input().corner_vertices.size();
    auto mask=[&](const std::optional<EdgeMask>& given,const char* field) {
        if(!given)return EdgeMask(ne,0);
        if(given->size()!=ne)throw ContractError(ErrorCode::invalid_input,field,"edge mask does not match the mesh edges");
        for(auto x:*given)if(x>1)throw ContractError(ErrorCode::invalid_input,field,"edge mask values must be zero or one");
        return *given;
    };
    HardConstraints out;out.locked=mask(locked,"constraints.locked");out.required=mask(required,"constraints.required_edges");
    out.forbidden=mask(forbidden,"constraints.forbidden_edges");
    for(auto a:anchors) {
        if(a<0 || static_cast<std::uint64_t>(a)>=nv)throw ContractError(ErrorCode::invalid_input,"constraints.anchors","anchor vertex out of range");
        out.anchors.push_back(static_cast<MeshIndex>(a));
    }
    std::sort(out.anchors.begin(),out.anchors.end());out.anchors.erase(std::unique(out.anchors.begin(),out.anchors.end()),out.anchors.end());
    for(const auto& [c,p]:pins) {
        if(c>=nc || !std::isfinite(p[0]) || !std::isfinite(p[1]))throw ContractError(ErrorCode::invalid_input,"constraints.pins","pins need a valid corner and finite coordinates");
        out.pins.emplace(c,p);
    }
    return out;
}
nlohmann::json constraint_violations(const MeshSnapshot& mesh,const HardConstraints& hc,std::span<const std::uint8_t> cuts) {
    const auto& t=mesh.tables();const auto ne=t.edge_vertices.size();
    if(cuts.size()!=ne)throw ContractError(ErrorCode::invalid_input,"constraints.cuts","cut mask does not match the mesh edges");
    std::vector<MeshIndex> locked,required,forbidden,anchors;
    for(std::size_t e=0;e<ne;++e) {
        const auto id=static_cast<MeshIndex>(e);
        if(e<hc.locked.size() && hc.locked[e] && !cuts[e])locked.push_back(id);
        if(e<hc.required.size() && hc.required[e] && !cuts[e])required.push_back(id);
        // A physical boundary (or nonmanifold) cut is a topology fact, not an internal seam.
        if(e<hc.forbidden.size() && hc.forbidden[e] && cuts[e] && !always_cut(t,e))forbidden.push_back(id);
    }
    // A manual anchor on the internal seam network of the run origin stays on it: a proposal
    // may reroute the seams through it but never abandon it.
    for(auto v:hc.seam_anchors) {
        bool on_seam=false;
        for(auto e:t.vert_edges[v])on_seam=on_seam || (cuts[e] && !t.edge_boundary[e]);
        if(!on_seam)anchors.push_back(v);
    }
    Json out=Json::object();
    if(!locked.empty())out["locked_edges"]=locked;
    if(!required.empty())out["required_edges"]=required;
    if(!forbidden.empty())out["forbidden_edges"]=forbidden;
    if(!anchors.empty())out["anchors"]=anchors;
    return out;
}
// ------------------------------------------------------------ states
namespace {
// Validity of a state from its chart outcomes and its feasibility (violations of the hard
// constraints already in `s.violations`). GB: inherited strict-only defects of candidate zero
// keep `valid` false and are listed in `strict_defects`; a feasible state whose only failures
// they are is `strict_only` (admissible, ranked after states with fewer of them).
void settle_validity(RunState& s) {
    s.violations.erase("charts");s.violations.erase("strict_defects");
    const auto ids=all_charts(*s.surface);const auto [failure,detail]=chart_failures(s.charts,ids);
    s.strict_defects=strict_defect_charts(s.charts);
    s.valid=s.feasible && failure.empty() && s.strict_defects.empty();
    s.strict_only=s.feasible && failure.empty() && !s.strict_defects.empty();
    // GB N2: a strict-only pin defect is named like the check that sees it (pins first, as chart_failures).
    const bool pin_defect=std::any_of(s.strict_defects.begin(),s.strict_defects.end(),[&](MeshIndex ch) {return !s.charts[ch].pins_kept;});
    s.invalid_reason=!s.feasible?violation_reason(s.violations):!failure.empty()?failure:s.strict_defects.empty()?"":pin_defect?"pin_moved":"overlap";
    if(!failure.empty())s.violations["charts"]=detail;
    if(!s.strict_defects.empty())s.violations["strict_defects"]=s.strict_defects;
}
}
RunState measure_state(ProposalContext& ctx,std::shared_ptr<const CutSurface> surface,std::vector<Vec2> uv,Revision revision,nlohmann::json provenance) {
    if(!surface || surface->mesh_owner()!=ctx.mesh)throw ContractError(ErrorCode::invalid_input,"enhanced.surface","state surface belongs to another mesh");
    if(uv.size()!=surface->tables().sector_of_corner.size())throw ContractError(ErrorCode::invalid_input,"enhanced.uv","one UV per corner required");
    RunState s;s.revision=revision;s.cuts=surface->tables().cut;s.surface=std::move(surface);s.uv=std::move(uv);
    auto measured=measure_charts(ctx,s.surface,s.uv,std::nullopt);s.charts=std::move(measured.outcomes);
    s.objective=seam_objective(*s.surface,ctx.raw_exposure,ctx.exempt,ctx.options.layout,ctx.options.enhanced.objective);
    s.budget=summarize_budget(s.charts);
    s.violations=constraint_violations(*ctx.mesh,ctx.constraints,s.cuts);s.feasible=s.violations.empty();
    settle_validity(s);
    if(!provenance.is_null())s.provenance.push_back(std::move(provenance));
    return s;
}
RunState begin_enhancement(ProposalContext& ctx,std::shared_ptr<const CutSurface> surface,std::vector<Vec2> uv,nlohmann::json provenance) {
    if(!surface)throw ContractError(ErrorCode::invalid_input,"enhanced.surface","missing candidate zero");
    // Anchors on candidate zero's internal seams are held there for the whole pass.
    const auto& t=ctx.mesh->tables();const auto& cut=surface->tables().cut;ctx.constraints.seam_anchors.clear();
    for(auto v:ctx.constraints.anchors) {
        bool on_seam=false;for(auto e:t.vert_edges[v])on_seam=on_seam || (cut[e] && !t.edge_boundary[e]);
        if(on_seam)ctx.constraints.seam_anchors.push_back(v);
    }
    ctx.envelope=BaselineEnvelope();
    const auto revision=ctx.next_revision();
    // T39: the operation's physical scale is frozen here, on candidate zero, for the whole pass
    // (a later re-anchored envelope or committed cut never changes it; Section 20).
    ctx.scale=std::make_shared<const PhysicalScale>(PhysicalScale::freeze(*surface,ctx.exempt,ctx.options,revision));
    if(provenance.is_object() && !provenance.contains("candidate"))provenance["candidate"]=0;
    auto state=measure_state(ctx,surface,std::move(uv),revision,std::move(provenance));
    // GB (Section 11.3 policy): candidate zero's strict-only defects are its own, reported
    // charts, not a refusal native_baseline would not make; every later trial inherits the flag
    // only while it leaves the chart exactly as it is (evaluate_trial).
    bool marked=false;
    for(auto& o:state.charts)if(strict_only_failure(o)) {o.strict_defect=true;marked=true;}
    if(marked)settle_validity(state);
    anchor_envelope(ctx,state);
    for(MeshIndex ch=0;ch<state.charts.size();++ch)state.charts[ch].ancestors={ch};
    return state;
}
void anchor_envelope(ProposalContext& ctx,const RunState& state) {
    QualityPolicy policy;policy.exempt_faces=face_ids(ctx.exempt);
    const QualitySnapshot snapshot(quality_geometry(state.surface),CornerUvView{state.uv},policy);
    ctx.envelope=BaselineEnvelope(state.revision,snapshot,state.uv,ctx.options.p95_growth,ctx.options.max_growth);
}
// ------------------------------------------------------------ proposals
namespace {
// Trial topology, changed-chart solve and measurement at the trial's current revision.
// T40: `bound` (optional) refuses, after the exact topology and before the solve, a trial
// whose exact objective proves it cannot win (J depends on the cuts only).
void evaluate_trial(const RunState& parent,EvaluatedProposal& e,ProposalContext& ctx,const DominanceBound* bound=nullptr) {
    const auto started=std::chrono::steady_clock::now();
    auto& trial=e.trial;const auto& mesh=*ctx.mesh;const auto gen=generator_name(e.proposal);
    auto& tally=ctx.counters.by_generator[gen];
    e.feasible=e.topology_ok=e.valid=e.inside_envelope=false;e.rejection.clear();e.detail=Json::object();e.changed_charts.clear();e.charts.clear();
    e.strict_defects.clear();
    trial.evaluated_revision=0;e.constraints_revision=ctx.constraints.revision;
    auto reject=[&](std::string reason,Json detail) {
        e.rejection=std::move(reason);e.detail=std::move(detail);++tally[e.rejection];
        e.seconds+=std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count();
    };
    throw_if_cancelled(ctx.cancel);
    // Hard constraints hold on the complete trial cut set, before anything is solved.
    const auto violations=constraint_violations(mesh,ctx.constraints,trial.cuts);
    if(!violations.empty())return reject(violation_reason(violations),violations);
    e.feasible=true;
    // Exact topology of the rebuilt trial surface.
    trial.surface=cut_surface(ctx.mesh,trial.cuts);
    const auto& ct=trial.surface->tables();const auto& pt=parent.surface->tables();
    std::map<std::string,MeshIndex> parent_keys;for(MeshIndex ch=0;ch<pt.n_charts;++ch)parent_keys.emplace(parent.surface->chart_key(ch).digest(),ch);
    std::set<std::string> resolve(e.proposal.resolve_chart_keys.begin(),e.proposal.resolve_chart_keys.end());
    std::vector<MeshIndex> changed;std::vector<std::optional<MeshIndex>> same(ct.n_charts);
    // A writer's map is measured as written: every chart whose corners it moved is changed.
    std::vector<std::uint8_t> moved(ct.n_charts,0),kept_topology(ct.n_charts,0);
    if(trial.uv_written)for(std::size_t c=0;c<trial.uv.size();++c)if(trial.uv[c]!=parent.uv[c])moved[ct.chart_of_corner[c]]=1;
    for(MeshIndex ch=0;ch<ct.n_charts;++ch) {
        const auto& key=trial.surface->chart_key(ch);const auto found=parent_keys.find(key.digest());
        kept_topology[ch]=found!=parent_keys.end() && parent.surface->chart_key(found->second)==key && !resolve.contains(key.digest());
        if(kept_topology[ch] && !moved[ch])same[ch]=found->second;
        else changed.push_back(ch);
    }
    e.changed_charts=changed;
    std::vector<MeshIndex> not_disk;
    const auto& g=mesh.geometry();
    for(auto ch:changed) {
        // GB: a writer that only moved a chart the parent already has (same key, not re-solved)
        // keeps the parent's topology: its map is judged by the measurements below, and a
        // non-disk chart native_baseline accepted is not refused for being moved (output hygiene,
        // strict-overlap repair). Charts whose topology the trial creates or re-solves must be disks.
        if(kept_topology[ch])continue;
        bool measured=false;for(auto f:ct.chart_faces[ch])measured=measured || (!(f<ctx.exempt.size() && ctx.exempt[f]) && g.face_degenerate_kind[f]==0);
        if(measured && !ct.chart_is_disk[ch])not_disk.push_back(ch);
    }
    if(!not_disk.empty())return reject("not_disk",{{"charts",not_disk}});
    // Exact pin compatibility on the trial's sectors: never averaged, never overwritten.
    std::map<MeshIndex,Vec2> first;std::vector<MeshIndex> conflicts;
    for(const auto& [c,p]:ctx.constraints.pins) {
        const auto s=ct.sector_of_corner[c];const auto [at,inserted]=first.emplace(s,p);
        if(!inserted && at->second!=p)conflicts.push_back(s);
    }
    if(!conflicts.empty())return reject("pin_conflict",{{"sectors",conflicts}});
    e.topology_ok=true;
    if(bound) {
        // The objective is exact here: no solve can change it. The map and budget stay
        // unmeasured; the trial is not committable (no evaluated revision).
        e.objective=seam_objective(*trial.surface,ctx.raw_exposure,ctx.exempt,ctx.options.layout,ctx.options.enhanced.objective);
        // GB (seam_objective.v2): U against the incumbent (the parent) is exact here too: it
        // depends on the two cut sets only.
        // GB (visual review V2): and the practical-effect floor against the incumbent (exact too).
        const auto sig=bound->significance && bound->incumbent_J?objective_significance(ctx,trial.cuts,parent.cuts,parent.objective.area,parent.objective.J)
                                                              :ObjectiveSignificance{};
        const double margin=sig.margin();
        if(bound->dominates(e.objective.J,e.proposal.identity.id,margin)) {
            Json detail={{"J",e.objective.J},{"bound",bound->to_json()},{"version",bounded_search_version}};
            if(bound->significance) {detail["uncertainty"]=sig.uncertainty;detail["practical_floor"]=sig.floor;detail["margin"]=margin;}
            return reject("dominated",std::move(detail));
        }
    }
    // Solve the changed charts; every other chart keeps the parent's map.
    if(!trial.uv_written) {
        trial.uv=parent.uv;
        if(!changed.empty()) {
            std::vector<Vec2> solved;
            try {solved=ctx.solve(trial.surface,changed,ctx.constraints.pins,e.proposal.initialization);}
            catch(const ContractError& error) {
                if(error.code==ErrorCode::cancelled)throw;
                return reject(error.code==ErrorCode::constraint_conflict?"pin_conflict":"solver_failure",{{"message",error.what()}});
            } catch(const std::exception& error) {return reject("solver_failure",{{"message",error.what()}});}
            if(solved.size()!=trial.uv.size())return reject("solver_failure",{{"message","solver returned a map of the wrong size"}});
            const std::set<MeshIndex> changed_set(changed.begin(),changed.end());
            for(std::size_t c=0;c<trial.uv.size();++c)if(changed_set.contains(ct.chart_of_corner[c]))trial.uv[c]=solved[c];
            e.solves=1;++ctx.counters.solves;ctx.counters.charts_solved+=changed.size();
        }
    }
    // Measure: changed charts afresh, unchanged charts keep the parent's outcome.
    auto measured=measure_charts(ctx,trial.surface,trial.uv,std::span<const MeshIndex>(changed));
    e.charts.resize(ct.n_charts);
    for(MeshIndex ch=0;ch<ct.n_charts;++ch)e.charts[ch]=same[ch]?parent.charts[*same[ch]]:measured.outcomes[ch];
    std::vector<MeshIndex> ids(ct.n_charts);std::iota(ids.begin(),ids.end(),MeshIndex{0});
    const auto [failure,detail]=chart_failures(e.charts,ids);
    e.objective=seam_objective(*trial.surface,ctx.raw_exposure,ctx.exempt,ctx.options.layout,ctx.options.enhanced.objective);
    e.budget=summarize_budget(e.charts);
    trial.evaluated_revision=trial.revision;
    if(!failure.empty())return reject(failure,detail);
    // GB: charts this trial inherits unchanged with candidate zero's strict-only defect (every
    // chart it solved or moved was measured afresh above and had to pass the strict check).
    e.strict_defects=strict_defect_charts(e.charts);
    e.valid=e.strict_defects.empty();
    // Run-origin envelope on the faces of every chart the origin does not have.
    if(!parent.feasible) {
        // The parent violates hard constraints: the envelope anchors on the first feasible state.
        e.envelope=EnvelopeCheck{};e.envelope.reason="parent_infeasible";e.inside_envelope=true;++tally["admissible"];
        e.seconds+=std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count();
        return;
    }
    const auto domain_charts=ctx.envelope.domain_charts(*trial.surface,trial.uv);std::vector<MeshIndex> domain_faces;
    for(auto ch:domain_charts)for(auto f:ct.chart_faces[ch])domain_faces.push_back(f);
    std::sort(domain_faces.begin(),domain_faces.end());
    std::vector<DistortionSample> samples;
    if(!domain_charts.empty()) {
        QualityPolicy policy;policy.exempt_faces=face_ids(ctx.exempt);
        samples=distortion_samples(QualitySnapshot(quality_geometry(trial.surface),CornerUvView{trial.uv},policy,std::span<const MeshIndex>(domain_charts)));
    }
    // GB m10: the relief-only rebases of the split faces apply to the relief reunion trial and to
    // states that adopted one; every other trial is measured against candidate zero.
    std::vector<std::size_t> rebases;
    if(ctx.envelope.staged_rebases())
        for(auto k:relief_rebases(parent,e.proposal))if(k<ctx.envelope.staged_rebases())rebases.push_back(k);
    e.envelope=ctx.envelope.check(domain_faces,samples,rebases);e.envelope.relief_rebases=rebases;
    if(!e.envelope.inside)return reject("quality_envelope",e.envelope.to_json());
    e.inside_envelope=true;++tally["admissible"];
    e.seconds+=std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count();
}
}
EvaluatedProposal evaluate_proposal(const RunState& parent,Proposal proposal,ProposalContext& ctx) {
    return evaluate_proposal(parent,std::move(proposal),ctx,nullptr);
}
EvaluatedProposal evaluate_proposal(const RunState& parent,Proposal proposal,ProposalContext& ctx,const DominanceBound* bound) {
    EvaluatedProposal e;e.parent_revision=parent.revision;
    const auto gen=generator_name(proposal);auto& tally=ctx.counters.by_generator[gen];++tally["tried"];++ctx.counters.proposals;
    auto reject=[&](std::string reason,Json detail) {e.rejection=std::move(reason);e.detail=std::move(detail);++tally[e.rejection];e.proposal=std::move(proposal);return e;};
    if(!parent.surface || parent.revision==0)throw ContractError(ErrorCode::invalid_input,"enhanced.parent","proposals need a measured parent state");
    if(proposal.identity.parent_revision!=parent.revision)
        return reject("stale_parent",{{"parent_revision",proposal.identity.parent_revision},{"current_revision",parent.revision}});
    throw_if_cancelled(ctx.cancel);
    const auto& t=ctx.mesh->tables();const auto ne=t.edge_vertices.size();
    // Representable delta: distinct original surface edges that change the effective cuts.
    std::vector<MeshIndex> bad,boundary,unchanged;std::set<std::int64_t> seen;
    auto check=[&](const std::vector<std::int64_t>& edges,bool add) {
        for(auto e:edges) {
            if(e<0 || static_cast<std::uint64_t>(e)>=ne || !seen.insert(e).second) {bad.push_back(e<0?0:static_cast<MeshIndex>(std::min<std::int64_t>(e,UINT32_MAX)));continue;}
            const auto id=static_cast<std::size_t>(e);
            if(always_cut(t,id))boundary.push_back(static_cast<MeshIndex>(id));
            else if((parent.cuts[id]!=0)==add)unchanged.push_back(static_cast<MeshIndex>(id));
        }
    };
    check(proposal.identity.added_edges,true);check(proposal.identity.removed_edges,false);
    if(!bad.empty())return reject("unrepresentable_edge",{{"edges",bad},{"why","out of range or repeated"}});
    if(!boundary.empty())return reject("physical_boundary",{{"edges",boundary}});
    if(!unchanged.empty())return reject("unrepresentable_edge",{{"edges",unchanged},{"why","not a change of the effective cuts"}});
    if(proposal.identity.added_edges.empty() && proposal.identity.removed_edges.empty() && proposal.initialization.empty())
        return reject("duplicate",{{"why","no delta and no initialization: candidate zero"}});
    // GB (m3): the initialization identity is a numerical contract. A solver that cannot honour
    // it refuses the proposal before any work, never solving it as its default under another name.
    if(!proposal.initialization.empty() && ctx.supports_initialization && !ctx.supports_initialization(proposal.initialization))
        return reject("unsupported_initialization",{{"initialization",proposal.initialization},
                                                     {"why","the run's solver honours only its default initialization"}});
    // T38: under strict symmetry the proposal and its mirror image are one delta, paired before
    // anything is checked or solved. A conflict refuses the pair with its reason, or (explicit
    // user policy disable_for_proposal) the unpaired proposal is evaluated with that label.
    if(strict_symmetry(ctx)) {
        const auto map=context_symmetry(ctx);
        auto paired=pair_proposal(proposal,parent,ctx.constraints,*map);
        if(paired.ok()) {++tally["symmetry_"+paired.status];proposal=std::move(paired.proposal);}
        else {
            ++tally["symmetry_conflict:"+paired.reason];
            if(ctx.options.enhanced.symmetry.on_conflict!="disable_for_proposal")return reject("symmetry_conflict",paired.detail);
            ++tally["symmetry_disabled"];
            auto record=paired.detail;record["version"]=symmetric_moves_version;record["status"]="disabled_for_proposal";
            if(!proposal.provenance.is_object())proposal.provenance=Json{{"generator_provenance",proposal.provenance}};
            proposal.provenance["symmetry"]=std::move(record);
        }
    }
    // Declared endpoint moves: anchors never slide; targets stay in the allowed set.
    for(const auto& m:proposal.endpoints) {
        if(std::binary_search(ctx.constraints.anchors.begin(),ctx.constraints.anchors.end(),m.from))return reject("fixed_anchor",{{"anchors",{m.from}}});
        if(proposal.allowed_endpoints && std::find(proposal.allowed_endpoints->begin(),proposal.allowed_endpoints->end(),m.to)==proposal.allowed_endpoints->end())
            return reject("endpoint_not_allowed",{{"endpoints",{m.to}}});
    }
    // Cheap pin screen before any rebuild: gluing that joins differently pinned sectors.
    std::vector<MeshIndex> removed;for(auto x:proposal.identity.removed_edges)removed.push_back(static_cast<MeshIndex>(x));
    if(pins_conflict_after_gluing(*parent.surface,ctx.constraints.pins,removed))return reject("pin_conflict",{{"edges",removed}});
    // Section 21 bound: expensive solves per affected origin chart. GB (m3): the charts a
    // proposal asks to re-solve are affected too (an initialization-only proposal counts).
    const auto affected=proposal_origin_charts(proposal,parent,ctx);
    if(proposal.identity.affected_chart_keys.empty()) {
        std::set<std::string> keys(proposal.resolve_chart_keys.begin(),proposal.resolve_chart_keys.end());
        for(const auto& list:{proposal.identity.added_edges,proposal.identity.removed_edges})
            for(auto x:list)for(auto f:t.edge_faces[static_cast<std::size_t>(x)])keys.insert(parent.charts[parent.surface->tables().chart_of_face[f]].key);
        proposal.identity.affected_chart_keys.assign(keys.begin(),keys.end());
    }
    for(const auto& key:affected) {
        const auto found=ctx.counters.solves_by_chart.find(key);
        if(found!=ctx.counters.solves_by_chart.end() && found->second>=static_cast<std::uint64_t>(ctx.options.enhanced.max_solves_per_chart))
            return reject("solve_budget_exhausted",{{"chart",key},{"limit",ctx.options.enhanced.max_solves_per_chart}});
    }
    e.trial.cuts=parent.cuts;
    for(auto x:proposal.identity.added_edges)e.trial.cuts[static_cast<std::size_t>(x)]=1;
    for(auto x:proposal.identity.removed_edges)e.trial.cuts[static_cast<std::size_t>(x)]=0;
    e.trial.revision=ctx.next_revision();e.trial.writers.push_back("proposal:"+gen);
    e.proposal=std::move(proposal);
    // A trial that throws after its solve (m4: the search isolates the failure) still used it.
    try {evaluate_trial(parent,e,ctx,bound);}
    catch(...) {if(e.solves)for(const auto& key:affected)++ctx.counters.solves_by_chart[key];throw;}
    if(e.solves)for(const auto& key:affected)++ctx.counters.solves_by_chart[key];
    return e;
}
std::set<std::string> proposal_origin_charts(const Proposal& p,const RunState& parent,const ProposalContext& ctx) {
    std::set<std::string> out;
    if(!parent.surface || !ctx.mesh)return out;
    const auto& t=ctx.mesh->tables();const auto ne=t.edge_vertices.size();
    const auto& origin=ctx.envelope.origin_chart_of_face();const auto& pt=parent.surface->tables();
    auto add_face=[&](MeshIndex f) {out.insert(origin.empty()?parent.charts[pt.chart_of_face[f]].key:"origin:"+std::to_string(origin[f]));};
    for(const auto& list:{p.identity.added_edges,p.identity.removed_edges})
        for(auto x:list) {
            if(x<0 || static_cast<std::uint64_t>(x)>=ne)continue;
            for(auto f:t.edge_faces[static_cast<std::size_t>(x)])add_face(f);
        }
    // GB (m3): charts re-solved on request (keys of the parent's charts), every face of each.
    if(!p.resolve_chart_keys.empty()) {
        const std::set<std::string> wanted(p.resolve_chart_keys.begin(),p.resolve_chart_keys.end());
        for(MeshIndex ch=0;ch<parent.charts.size() && ch<pt.n_charts;++ch)
            if(wanted.contains(parent.charts[ch].key))for(auto f:pt.chart_faces[ch])add_face(f);
    }
    return out;
}
bool isolated_failure(const std::exception& error) {
    const auto* contract=dynamic_cast<const ContractError*>(&error);
    if(!contract)return true;
    return contract->code!=ErrorCode::cancelled && contract->code!=ErrorCode::invalid_option && contract->code!=ErrorCode::unsupported_version;
}
nlohmann::json failure_record(const std::exception& error) {
    const auto* contract=dynamic_cast<const ContractError*>(&error);
    return {{"code",contract?std::string(error_code_name(contract->code)):std::string("exception")},{"field",contract?contract->field:std::string()},
            {"message",error.what()}};
}
void write_trial_cuts(EvaluatedProposal& e,EdgeMask cuts,const std::string& writer,ProposalContext& ctx) {
    if(cuts.size()!=ctx.mesh->tables().edge_vertices.size())throw ContractError(ErrorCode::invalid_input,"enhanced.writer","cut mask does not match the mesh edges");
    // The surface re-adds the physical boundary; the trial records the effective cuts.
    e.trial.surface=cut_surface(ctx.mesh,cuts);e.trial.cuts=e.trial.surface->tables().cut;
    e.trial.uv_written=false;e.trial.revision=ctx.next_revision();e.trial.writers.push_back(writer);
}
void write_trial_uv(EvaluatedProposal& e,std::vector<Vec2> uv,const std::string& writer,ProposalContext& ctx) {
    if(uv.size()!=ctx.mesh->input().corner_vertices.size())throw ContractError(ErrorCode::invalid_input,"enhanced.writer","one UV per corner required");
    e.trial.uv=std::move(uv);e.trial.uv_written=true;e.trial.revision=ctx.next_revision();e.trial.writers.push_back(writer);
}
void reevaluate_proposal(const RunState& parent,EvaluatedProposal& e,ProposalContext& ctx) {
    if(e.parent_revision!=parent.revision)throw ContractError(ErrorCode::stale_state,"enhanced.parent","the trial was built on another state revision");
    ++ctx.counters.by_generator[generator_name(e.proposal)]["reevaluated"];
    const auto solves=e.solves;evaluate_trial(parent,e,ctx);e.solves+=solves;
}
RunState commit_proposal(const RunState& parent,const EvaluatedProposal& e,ProposalContext& ctx) {
    if(e.parent_revision!=parent.revision)throw ContractError(ErrorCode::stale_state,"enhanced.commit","the trial was built on another state revision");
    if(e.trial.revision==0 || e.trial.evaluated_revision!=e.trial.revision)
        throw ContractError(ErrorCode::stale_state,"enhanced.commit","a writer changed the trial after it was evaluated; evaluate it again");
    // Painted constraints changed after the evaluation (T29): the proposal is stale.
    if(e.constraints_revision!=ctx.constraints.revision)
        throw ContractError(ErrorCode::stale_state,"enhanced.commit","the painted constraints changed after the trial was evaluated; evaluate it again");
    if(!e.rejection.empty())throw ContractError(ErrorCode::constraint_conflict,"enhanced.commit","the proposal was not admitted: "+e.rejection);
    // T38: under strict symmetry no half pair is ever committed: every declared change needs its
    // mirror half in the committed cuts, unless the user's policy explicitly disabled symmetry
    // for this proposal (labelled in its provenance).
    if(strict_symmetry(ctx)) {
        const auto map=context_symmetry(ctx);
        const auto record=e.proposal.provenance.is_object()?e.proposal.provenance.value("symmetry",Json::object()):Json::object();
        const bool disabled=record.value("status","")=="disabled_for_proposal" && ctx.options.enhanced.symmetry.on_conflict=="disable_for_proposal";
        if(!disabled)
            if(const auto missing=missing_mirror_half(e.proposal,e.trial.cuts,*map);!missing.is_null())
                throw ContractError(ErrorCode::constraint_conflict,"enhanced.commit","half pair: the mirror half of the proposal is missing "+missing.dump());
    }
    // What is committed is exactly what was solved and measured.
    if(!e.trial.surface || e.trial.surface->tables().cut!=e.trial.cuts || e.charts.size()!=e.trial.surface->tables().n_charts)
        throw ContractError(ErrorCode::stale_state,"enhanced.commit","trial cuts, solved surface and measurements disagree");
    RunState s;s.revision=e.trial.revision;s.cuts=e.trial.cuts;s.surface=e.trial.surface;s.uv=e.trial.uv;s.charts=e.charts;
    s.objective=e.objective;s.budget=e.budget;s.feasible=true;s.provenance=parent.provenance;
    settle_validity(s); // valid, or strict_only with the inherited strict defects the trial left unchanged
    s.provenance.push_back({{"revision",s.revision},{"parent_revision",parent.revision},{"proposal",e.proposal.identity.id},
                            {"generator",e.proposal.identity.generator},{"generator_version",e.proposal.identity.generator_version},
                            {"added_edges",e.proposal.identity.added_edges},{"removed_edges",e.proposal.identity.removed_edges},
                            {"initialization",e.proposal.initialization},{"writers",e.trial.writers},{"changed_charts",e.changed_charts},
                            {"provenance",e.proposal.provenance}});
    ++ctx.counters.commits;++ctx.counters.by_generator[generator_name(e.proposal)]["committed"];
    return s;
}
RunState rewrite_state_uv(const RunState& state,std::vector<Vec2> uv,const std::string& writer,ProposalContext& ctx,nlohmann::json* record,bool no_worse) {
    EvaluatedProposal e;e.parent_revision=state.revision;
    e.proposal.identity.parent_revision=state.revision;e.proposal.identity.generator=writer;e.proposal.initialization="writer:"+writer;
    ++ctx.counters.by_generator[writer]["tried"];++ctx.counters.proposals;
    e.trial.cuts=state.cuts;e.trial.surface=state.surface;e.trial.writers.push_back("state:"+std::to_string(state.revision));
    write_trial_uv(e,std::move(uv),writer,ctx);
    evaluate_trial(state,e,ctx);
    bool keep=e.admissible();
    if(keep && no_worse) {
        const CandidateView trial{1,true,&e.objective,&e.budget,e.strict_defects.size()},
                            current{0,state.admissible(),&state.objective,&state.budget,state.strict_defects.size()};
        if(compare_candidates(trial,current,ctx.options.enhanced.objective)>0) {keep=false;++ctx.counters.by_generator[writer]["ranks_after_state"];}
    }
    if(record) {*record=e.summary();if(!keep && e.admissible())(*record)["rejection"]="ranks_after_state";}
    return keep?commit_proposal(state,e,ctx):state;
}
EvaluatedProposal evaluate_written_map(const RunState& parent,Proposal identity,std::vector<Vec2> uv,const std::string& writer,ProposalContext& ctx) {
    if(!parent.surface || parent.revision==0)throw ContractError(ErrorCode::invalid_input,"enhanced.parent","written maps need a measured parent state");
    EvaluatedProposal e;e.parent_revision=parent.revision;
    identity.identity.parent_revision=parent.revision;
    const auto gen=generator_name(identity);++ctx.counters.by_generator[gen]["tried"];++ctx.counters.proposals;
    e.proposal=std::move(identity);
    e.trial.cuts=parent.cuts;e.trial.surface=parent.surface;e.trial.writers.push_back("state:"+std::to_string(parent.revision));
    write_trial_uv(e,std::move(uv),writer,ctx);
    evaluate_trial(parent,e,ctx);
    return e;
}
}
