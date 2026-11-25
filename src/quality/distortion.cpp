// SPDX-License-Identifier: GPL-3.0-or-later
#include "detail.hpp"
#include "cleanseams/detail/numpy_compat.hpp"
namespace cleanseams {
using namespace quality_detail;
nlohmann::json QualityReport::to_json() const {
    nlohmann::json result={{"triangles",triangles},{"flipped",flipped},{"folded",folded},{"degenerate",degenerate},
        {"nonfinite",nonfinite},{"slivers",slivers},{"collapsed_uv",collapsed_uv},{"exempt",exempt},
        {"degenerate_input",degenerate_input},{"flipped_exempt",flipped_exempt},{"density",density},{"density_spread",density_spread},
        {"anisotropy_p50",anisotropy_p50},{"anisotropy_p95",anisotropy_p95},{"anisotropy_p99",anisotropy_p99},{"anisotropy_max",anisotropy_max},
        {"area_p50",area_p50},{"area_p95",area_p95},{"area_p99",area_p99},{"area_max",area_max}};
    auto unavailable=nlohmann::json::array();
    for(auto it=result.begin();it!=result.end();++it) if(it->is_number_float()) {
        const double v=it->get<double>();if(std::isfinite(v)) continue;
        unavailable.push_back({{"path","/"+it.key()},{"reason",std::isnan(v)?"nan":v>0?"positive_infinity":"negative_infinity"}});it.value()=nullptr;
    }
    result["unavailable"]=std::move(unavailable);return result;
}
QualityReport QualitySnapshot::report(std::optional<MeshIndex> chart,std::optional<double> density,bool per_chart) const {
    require(!chart || std::binary_search(selected_charts_.begin(),selected_charts_.end(),*chart),"quality.chart","chart was not measured");
    require(!density || (std::isfinite(*density) && *density>=0),"quality.density","density must be finite and nonnegative");
    QualityReport r;
    std::vector<MeshIndex> all_t,all_p;
    if(!chart) {all_t.resize(triangles_.size());std::iota(all_t.begin(),all_t.end(),0);all_p.resize(polygons_.size());std::iota(all_p.begin(),all_p.end(),0);}
    const std::span<const MeshIndex> ts=chart?chart_triangles_[*chart]:std::span<const MeshIndex>(all_t);
    const std::span<const MeshIndex> ps=chart?chart_polygons_[*chart]:std::span<const MeshIndex>(all_p);
    r.triangles=ts.size();if(ts.empty()) return r;
    for(auto id:ps) {
        const auto& p=polygons_[id];r.exempt+=p.exempt;r.degenerate_input+=p.degenerate_input && p.exempt;
        r.flipped+=p.flipped && !p.exempt;r.flipped_exempt+=p.flipped && p.exempt;
    }
    // One chart report needs one accumulator, not an array over every chart.
    const auto n=chart?1:geometry_->cut().tables().n_charts;
    std::vector<double> d3(n,0),duv(n,0),ok_a3,ok_auv;ok_a3.reserve(ts.size());ok_auv.reserve(ts.size());
    std::vector<MeshIndex> ok;ok.reserve(ts.size());
    for(auto id:ts) {
        const auto& t=triangles_[id];r.nonfinite+=!t.finite && !t.exempt;r.slivers+=t.sliver;
        r.degenerate+=t.area3<=1e-14 && !t.exempt;
        if(t.finite && !t.sliver && !t.exempt) {
            r.folded+=t.area2 < -1e-12*t.area3;r.collapsed_uv+=std::abs(t.area2)<=1e-12*t.area3;
            if(finite(t.singular_values)) {
                ok.push_back(id);ok_a3.push_back(t.area3);ok_auv.push_back(std::abs(t.area2));
                const auto ch=chart?0:t.chart;d3[ch]+=t.area3;duv[ch]+=std::abs(t.area2);
            }
        }
    }
    // Source a2[ok].sum() and abs(uv_a2[ok]).sum(): NumPy pairwise sums in triangle order.
    const double a3=detail::numpy_sum(ok_a3),auv=detail::numpy_sum(ok_auv);
    const double global_density=a3>0?auv/a3:1;r.density=density.value_or(global_density);
    std::vector<double> dens(n,global_density);double low=std::numeric_limits<double>::infinity(),high=-std::numeric_limits<double>::infinity();bool present=false;
    if(!density && per_chart) {
        for(std::size_t ch=0;ch<n;++ch) if(d3[ch]>0) {
            dens[ch]=duv[ch]/std::max(d3[ch],1e-300);low=std::min(low,dens[ch]);high=std::max(high,dens[ch]);present=true;
        }
        if(present && low>0) r.density_spread=high/low;
    }
    std::vector<double> weights,aniso,area;weights.reserve(ok.size());aniso.reserve(ok.size());area.reserve(ok.size());
    for(auto id:ok) {
        const auto& t=triangles_[id];weights.push_back(t.area3);
        double a=t.singular_values[0]/std::max(t.singular_values[1],1e-300);
        const double den=density?*density:per_chart?dens[chart?0:t.chart]:global_density;
        const double ratio=t.singular_values[0]*t.singular_values[1]/std::max(den,1e-300);
        double b=std::max(ratio,1/std::max(ratio,1e-300));
        if(!std::isfinite(a)) a=1e9;if(!std::isfinite(b)) b=1e9;aniso.push_back(a);area.push_back(b);
    }
    const WeightedDistribution anisotropy(aniso,weights),area_distribution(area,weights);
    const std::array<double,3> queries{50,95,99};const auto aq=anisotropy.percentiles(queries),bq=area_distribution.percentiles(queries);
    r.anisotropy_p50=aq[0];r.anisotropy_p95=aq[1];r.anisotropy_p99=aq[2];r.area_p50=bq[0];r.area_p95=bq[1];r.area_p99=bq[2];
    if(!ok.empty()) {r.anisotropy_max=*std::max_element(aniso.begin(),aniso.end());r.area_max=*std::max_element(area.begin(),area.end());}
    return r;
}
std::array<std::vector<double>,3> QualitySnapshot::chart_density() const {
    const auto n=geometry_->cut().tables().n_charts;
    std::array<std::vector<double>,3> result;result[0].assign(n,quality_nan);result[1].assign(n,0);result[2].assign(n,0);
    for(const auto& t:triangles_) if(std::isfinite(t.area2) && finite(t.singular_values) && !t.sliver && !t.exempt) {
        result[1][t.chart]+=t.area3;result[2][t.chart]+=std::abs(t.area2);
    }
    for(std::size_t ch=0;ch<n;++ch) if(result[1][ch]>0) result[0][ch]=result[2][ch]/result[1][ch];
    return result;
}
std::vector<double> QualitySnapshot::chart_anisotropy_p95() const {
    const auto n=geometry_->cut().tables().n_charts;std::vector<double> result(n,quality_nan);
    for(auto ch:selected_charts_) {
        std::vector<double> values,weights;
        for(auto id:chart_triangles_[ch]) {
            const auto& t=triangles_[id];if(!t.finite || !finite(t.singular_values) || t.area3<=1e-14 || t.exempt) continue;
            auto a=t.singular_values[0]/std::max(t.singular_values[1],1e-300);if(!std::isfinite(a)) a=1e9;
            values.push_back(a);weights.push_back(t.area3);
        }
        result[ch]=WeightedDistribution(values,weights).percentile(95);
    }
    return result;
}
}
