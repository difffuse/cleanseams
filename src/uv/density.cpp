// SPDX-License-Identifier: GPL-3.0-or-later
// Source core/output.py: normalise_density.
#include "cleanseams/output.hpp"
#include "cleanseams/cache.hpp"
#include "cleanseams/errors.hpp"
#include "cleanseams/detail/numpy_compat.hpp"
#include <algorithm>
#include <cmath>
#include <set>
namespace cleanseams {
namespace {
nlohmann::json number(double x) {return nlohmann::json(x);}
}
nlohmann::json DensityNormalization::info() const {
    return {{"density_spread_before",number(density_spread_before)},{"target_density",number(target_density)},
            {"protected_charts",protected_charts},{"density_outliers",density_outliers}};
}
DensityNormalization normalise_density(const std::shared_ptr<const CutSurface>& cs,std::span<const Vec2> uv,
                                       std::span<const MeshIndex> pinned_corners,std::span<const std::uint8_t> exempt) {
    if(!cs)throw ContractError(ErrorCode::invalid_input,"density.cut","cut surface owner required");
    const auto& mesh=cs->mesh();const auto& t=cs->tables();
    const auto nc=mesh.input().corner_vertices.size(),nf=mesh.input().face_offsets.size()-1;
    if(uv.size()!=nc)throw ContractError(ErrorCode::invalid_input,"density.uv","one UV per original corner required");
    if(!exempt.empty() && exempt.size()!=nf)throw ContractError(ErrorCode::invalid_input,"density.exempt","one flag per face required");
    std::set<MeshIndex> fixed;
    for(auto c:pinned_corners) {
        if(c>=nc)throw ContractError(ErrorCode::invalid_input,"density.pins","pinned corner out of range");
        fixed.insert(t.chart_of_corner[c]);
    }
    // quality.chart_density: the same meaningful triangles as distortion_report.
    QualityPolicy policy;for(std::size_t f=0;f<exempt.size();++f)if(exempt[f])policy.exempt_faces.push_back(static_cast<TransferIndex>(f));
    const QualitySnapshot snapshot(quality_geometry(cs),CornerUvView{uv},policy);
    const auto [density,area3,area2]=snapshot.chart_density();
    std::vector<MeshIndex> eligible;std::vector<double> a2,a3,d;
    for(MeshIndex ch=0;ch<t.n_charts;++ch)
        if(std::isfinite(density[ch]) && density[ch]>0) {eligible.push_back(ch);a2.push_back(area2[ch]);a3.push_back(area3[ch]);d.push_back(density[ch]);}
    DensityNormalization out;out.uv.assign(uv.begin(),uv.end());
    if(!eligible.empty()) {
        out.target_density=detail::numpy_sum(a2)/detail::numpy_sum(a3);
        out.density_spread_before=*std::max_element(d.begin(),d.end())/ *std::min_element(d.begin(),d.end());
    }
    const double target=out.target_density;
    for(auto ch:eligible) {
        if(fixed.contains(ch))continue;
        std::vector<MeshIndex> ids;for(auto f:t.chart_faces[ch])for(auto c:mesh.face_corners(f))ids.push_back(c);
        Vec2 centre{0,0};for(auto c:ids)for(std::size_t k=0;k<2;++k)centre[k]+=out.uv[c][k];
        for(std::size_t k=0;k<2;++k)centre[k]/=static_cast<double>(ids.size());
        const double scale=std::sqrt(target/density[ch]);
        for(auto c:ids)for(std::size_t k=0;k<2;++k)out.uv[c][k]=(out.uv[c][k]-centre[k])*scale+centre[k];
    }
    out.protected_charts.assign(fixed.begin(),fixed.end());
    for(auto ch:eligible)if(density[ch]>target*50 || density[ch]<target/50)out.density_outliers.push_back(ch);
    return out;
}
}
