// SPDX-License-Identifier: GPL-3.0-or-later
// Source core/target.py.
#include "cleanseams/target.hpp"
#include "cleanseams/cache.hpp"
#include "cleanseams/errors.hpp"
#include <algorithm>
#include <cmath>
#include <tuple>
namespace cleanseams {
namespace {
[[noreturn]] void bad_band() {
    throw ContractError(ErrorCode::invalid_option,"target_islands","target_islands must be two positive integers in ascending order");
}
// One element of the source's np.asarray(value) as a float64 comparison value.
double band_value(const nlohmann::json& x) {
    if(x.is_boolean())return x.get<bool>()?1.0:0.0;
    if(x.is_number_integer() || x.is_number_unsigned()) {
        if(x.is_number_unsigned() && x.get<std::uint64_t>()>static_cast<std::uint64_t>(std::numeric_limits<std::int64_t>::max()))bad_band();
        return static_cast<double>(x.get<std::int64_t>());
    }
    if(x.is_number_float()) {
        const double value=x.get<double>();
        // astype(int) must give back the same value (finite, integral, in int64 range).
        if(!std::isfinite(value) || value!=std::trunc(value) || std::abs(value)>=0x1p63)bad_band();
        return value;
    }
    bad_band();
}
struct CachedMeasurement {ChartKey key;ChartMeasurement measurement;};
struct TrialCache {
    std::map<std::string,CachedMeasurement> maps;
    std::uint64_t hits=0;
    const ChartMeasurement* find(const ChartKey& key) const {
        const auto found=maps.find(key.digest());
        return found!=maps.end() && found->second.key==key?&found->second.measurement:nullptr;
    }
};
struct Candidate {
    std::tuple<bool,bool,std::size_t,std::int64_t,double,double> score;
    std::size_t index=0;
    GeneratedLayout layout;
    LayoutOptions options;
};
}
std::optional<IslandBand> island_band(const nlohmann::json& value) {
    if(value.is_null())return std::nullopt;
    if(!value.is_array() || value.size()!=2)bad_band();
    const double low=band_value(value[0]),high=band_value(value[1]);
    if(!(1<=low && low<=high))bad_band();
    return IslandBand{static_cast<std::int64_t>(low),static_cast<std::int64_t>(high)};
}
TargetSearch generate_for_target(const std::shared_ptr<const MeshSnapshot>& mesh,const LayoutOptions& opts,IslandBand target,
                                 const EvaluatorFactory& factory,const SeedingInput& input,LayoutTrace* trace) {
    if(!mesh)throw ContractError(ErrorCode::invalid_input,"target.mesh","missing mesh snapshot");
    if(!factory)throw ContractError(ErrorCode::invalid_input,"target.evaluator","missing evaluator factory");
    const auto [lo,hi]=*island_band(nlohmann::json::array({target[0],target[1]}));
    const auto& m=*mesh;
    auto cache=std::make_shared<TrialCache>();
    std::vector<Candidate> candidates;nlohmann::json attempts=nlohmann::json::array();
    double left=.6,right=1.6,scale=1.;
    for(std::size_t index=0;index<5;++index) {
        if(index && input.deadline.expired())break;
        LayoutOptions layout=opts;
        layout.max_anisotropy_p95=1+(opts.max_anisotropy_p95-1)*scale;layout.max_area_p95=1+(opts.max_area_p95-1)*scale;
        layout.strict_budget=true;layout.absorb_relax=1.;layout.min_island_fraction=opts.min_island_fraction*scale;
        const auto raw=factory(layout);
        if(!raw)throw ContractError(ErrorCode::invalid_input,"target.evaluator","the evaluator factory returned no evaluator");
        // Measured maps are shared across trials; the verdict is re-judged at this
        // trial's budget and the caller's relax multipliers.
        ChartEvaluator evaluate=[cache,raw,layout](const std::shared_ptr<const CutSurface>& cs,std::span<const MeshIndex> charts,
                                                   const std::map<MeshIndex,double>& relax) {
            std::vector<MeshIndex> missing;for(auto ch:charts)if(!cache->find(cs->chart_key(ch)))missing.push_back(ch);
            cache->hits+=charts.size()-missing.size();
            if(!missing.empty())
                for(auto& [ch,measured]:raw(cs,missing,relax)) {
                    const auto& key=cs->chart_key(ch);
                    cache->maps.insert_or_assign(key.digest(),CachedMeasurement{key,std::move(measured)});
                }
            ChartMeasurements result;
            for(auto ch:charts) {
                const auto* hit=cache->find(cs->chart_key(ch));if(!hit)continue;
                auto measured=*hit;const auto found=relax.find(ch);
                measured.within_budget=within_budget(measured.report,measured.crossings,layout,found==relax.end()?1.:found->second);
                result[ch]=std::move(measured);
            }
            return result;
        };
        SeedingInput trial=input;trial.evaluate=&evaluate;
        auto generated=generate_seams(mesh,layout,trial,trace);
        const auto& ct=generated.cs->tables();const auto& report=generated.report;
        bool uv_finite=report.uv.has_value();
        if(uv_finite)for(const auto& p:*report.uv)uv_finite=uv_finite && std::isfinite(p[0]) && std::isfinite(p[1]);
        std::shared_ptr<const QualityGeometry> geometry;
        std::vector<MeshIndex> bad;double worst=1.;
        for(MeshIndex ch=0;ch<ct.n_charts;++ch) {
            const auto faces=ct.chart_faces[ch];
            if(std::all_of(faces.begin(),faces.end(),[&](auto f) {return m.geometry().face_degenerate_kind[f]!=0;}))continue;
            std::optional<std::pair<QualityReport,std::uint64_t>> metric;
            if(const auto* hit=cache->find(generated.cs->chart_key(ch)))metric.emplace(hit->report,hit->crossings);
            else if(uv_finite) {
                if(!geometry)geometry=quality_geometry(generated.cs);
                const std::array<MeshIndex,1> one{ch};const CornerUvView view{*report.uv};
                metric.emplace(QualitySnapshot(geometry,view,{},std::span<const MeshIndex>(one)).report(),
                               boundary_crossings(*generated.cs,view,{},std::span<const MeshIndex>(one)).counts[ch]);
            }
            if(!metric || !within_budget(metric->first,metric->second,layout))bad.push_back(ch);
            if(metric && metric->first.anisotropy_p95>worst)worst=metric->first.anisotropy_p95;
        }
        const auto charts=static_cast<std::int64_t>(ct.n_charts);
        const auto distance=std::max<std::int64_t>({lo-charts,charts-hi,0});
        attempts.push_back({{"scale",scale},{"charts",charts},{"over_budget_charts",bad},{"timed_out",report.timed_out},
                            {"anisotropy_limit",layout.max_anisotropy_p95},{"area_limit",layout.max_area_p95}});
        if(trace) {
            std::vector<MeshIndex> seams;for(std::size_t e=0;e<generated.seams.size();++e)if(generated.seams[e])seams.push_back(static_cast<MeshIndex>(e));
            trace->events.push_back({{"stage","target_trial"},{"index",index},{"charts",charts},{"seams",seams},{"timed_out",report.timed_out},
                                     {"max_anisotropy_p95",layout.max_anisotropy_p95},{"max_area_p95",layout.max_area_p95},
                                     {"min_island_fraction",layout.min_island_fraction}});
        }
        const bool timed_out=report.timed_out;
        candidates.push_back({{timed_out,!bad.empty(),bad.size(),distance,std::abs(scale-1),worst},index,std::move(generated),layout});
        if(distance==0 && bad.empty() && !timed_out)break;
        // Check the reachable end of the range before bisecting. A repeated scale
        // tells us the count is outside the permitted budget range.
        double next;
        if(charts>hi) {left=scale>left?scale:left;next=index==0?right:.5*(left+right);}
        else {right=scale<right?scale:right;next=index==0?left:.5*(left+right);}
        if(std::any_of(attempts.begin(),attempts.end(),[&](const nlohmann::json& row) {return std::abs(next-row.at("scale").get<double>())<1e-6;}))break;
        scale=next;
    }
    auto best=std::min_element(candidates.begin(),candidates.end(),[](const Candidate& a,const Candidate& b) {
        return std::tie(a.score,a.index)<std::tie(b.score,b.index);
    });
    const auto& chosen=attempts.at(best->index);
    const auto charts=static_cast<std::int64_t>(best->layout.cs->tables().n_charts);
    const bool within=lo<=charts && charts<=hi,over=!chosen.at("over_budget_charts").empty(),timed_out=chosen.at("timed_out").get<bool>();
    TargetSearch out;
    out.cache_hits=cache->hits;
    out.info={{"requested",{lo,hi}},{"attempts",attempts},{"selected",best->index},{"scale",chosen.at("scale")},{"cache_hits",cache->hits},
              {"status",within && !over && !timed_out?"reached":"unreached"},
              {"reason",timed_out?"time budget":over?"distortion budget":!within?"merge constraints":""}};
    out.selected_layout=best->options;out.layout=std::move(best->layout);
    return out;
}
}
