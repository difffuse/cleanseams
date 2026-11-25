// SPDX-License-Identifier: GPL-3.0-or-later
// Source core/classify.py and the auto-profile step of core/pipeline.py.
#include "cleanseams/classify.hpp"
#include "cleanseams/axis.hpp"
#include "cleanseams/errors.hpp"
#include "cleanseams/detail/mesh_math.hpp"
#include "cleanseams/detail/numpy_compat.hpp"
#include <algorithm>
#include <cmath>
#include <numbers>
namespace cleanseams {
namespace {
using detail::sub;using detail::dot;using detail::norm;
std::size_t face_count(const MeshSnapshot& m) {return m.input().face_offsets.size()-1;}
}
PlanarCoverage planar_coverage(const MeshSnapshot& m) {
    const auto& g=m.geometry();const auto nf=face_count(m);
    const auto& weights=g.face_area;const double total=detail::numpy_sum(weights);
    if(total<=0)return {};
    // Area samples: searchsorted(cumsum(area), (k + .5) * total / 192), unique.
    std::vector<double> cumulative(nf);double running=0;
    for(std::size_t f=0;f<nf;++f) {running+=weights[f];cumulative[f]=running;}
    std::vector<std::size_t> samples;
    for(int k=0;k<192;++k) {
        const double value=((static_cast<double>(k)+.5)*total)/192;
        const auto index=static_cast<std::size_t>(std::lower_bound(cumulative.begin(),cumulative.end(),value)-cumulative.begin());
        // The source indexes the face tables with it; a sample past the last face is
        // an IndexError there (unreachable for finite areas).
        if(index>=nf)throw ContractError(ErrorCode::numerical_failure,"classify.samples","area sample beyond the last face");
        samples.push_back(index);
    }
    std::sort(samples.begin(),samples.end());samples.erase(std::unique(samples.begin(),samples.end()),samples.end());
    const auto& v=m.input().positions;
    Vec3 lo=v.empty()?Vec3{0,0,0}:v[0],hi=lo;
    for(const auto& p:v)for(std::size_t d=0;d<3;++d) {lo[d]=std::min(lo[d],p[d]);hi[d]=std::max(hi[d],p[d]);}
    const double extent=norm(sub(hi,lo)),diagonal=1e-15>extent?1e-15:extent;
    const double cone=std::cos(10*(std::numbers::pi/180.0)),band=.008*diagonal;
    std::vector<std::uint8_t> remaining(nf,1);PlanarCoverage out;double covered=0;
    std::vector<std::uint8_t> near(nf),best_mask;
    for(int round=0;round<8;++round) {
        std::optional<double> best;
        for(auto sample:samples) {
            const auto& normal=g.face_normal[sample];std::vector<double> area;
            for(std::size_t f=0;f<nf;++f) {
                near[f]=remaining[f] && dot(g.face_normal[f],normal)>cone &&
                        std::abs(dot(sub(g.face_center[f],g.face_center[sample]),normal))<band;
                if(near[f])area.push_back(weights[f]);
            }
            const double sum=detail::numpy_sum(area);
            if(!best || sum>*best) {best=sum;best_mask=near;}
        }
        if(!best || *best<.035*total)break;
        covered+=*best;++out.clusters;
        for(std::size_t f=0;f<nf;++f)if(best_mask[f])remaining[f]=0;
    }
    out.coverage=covered/total;
    return out;
}
AssetClass asset_class(const MeshSnapshot& m,const AssetClassInput& input) {
    const auto& t=m.tables();const auto& g=m.geometry();const auto nf=face_count(m);const auto& v=m.input().positions;
    std::size_t quads=0;for(std::size_t f=0;f<nf;++f)quads+=m.input().face_offsets[f+1]-m.input().face_offsets[f]==4;
    const auto regular_count=static_cast<std::size_t>(std::count_if(t.vert_regular.begin(),t.vert_regular.end(),[](auto x) {return x!=0;}));
    const double quad_fraction=nf?static_cast<double>(quads)/static_cast<double>(nf):0.;
    const double regular_fraction=v.empty()?0.:static_cast<double>(regular_count)/static_cast<double>(v.size());
    const auto coverage=planar_coverage(m);
    nlohmann::json evidence{{"quad_fraction",quad_fraction},{"regular_fraction",regular_fraction},{"planar_coverage",coverage.coverage},
                            {"planar_clusters",coverage.clusters},{"symmetry_supplied",input.plane_supplied}};
    if(input.features) {
        if(input.features->size()!=g.edge_length.size())throw ContractError(ErrorCode::invalid_input,"classify.features","edge mask does not match the mesh edges");
        std::vector<double> hard;for(std::size_t e=0;e<g.edge_length.size();++e)if((*input.features)[e])hard.push_back(g.edge_length[e]);
        const double all=detail::numpy_sum(g.edge_length);
        evidence["feature_length_fraction"]=detail::numpy_sum(hard)/(1e-15>all?1e-15:all);
    }
    auto result=[&](const char* profile,const char* reason) {auto r=evidence;r["reason"]=reason;return AssetClass{profile,r};};
    if(coverage.coverage>=.65 && (coverage.clusters>=4 || coverage.coverage>=.98))return result("hard","broad planar panels");
    // A sphere has no unique revolution axis; constant-radius vertex samples are
    // profile evidence, not a sphere diagnosis.
    if(!v.empty()) {
        Vec3 lo=v[0],hi=v[0];for(const auto& p:v)for(std::size_t d=0;d<3;++d) {lo[d]=std::min(lo[d],p[d]);hi[d]=std::max(hi[d],p[d]);}
        const Vec3 centre{.5*(lo[0]+hi[0]),.5*(lo[1]+hi[1]),.5*(lo[2]+hi[2])};
        std::vector<double> radius(v.size());for(std::size_t i=0;i<v.size();++i)radius[i]=norm(sub(v[i],centre));
        const double n=static_cast<double>(radius.size()),mean=detail::numpy_sum(radius)/n;
        std::vector<double> squared(radius.size());for(std::size_t i=0;i<radius.size();++i) {const double d=radius[i]-mean;squared[i]=d*d;}
        const double spread=std::sqrt(detail::numpy_sum(squared)/n);
        if(mean>0 && spread/mean<.025) {auto r=result("hard","constant-radius vertex samples");r.rationale["radius_residual"]=spread/mean;return r;}
    }
    std::vector<double> residuals;
    if(input.axis_fit_residuals)residuals=*input.axis_fit_residuals;
    else if(auto fit=revolution_frame(m,std::nullopt,.02))residuals.push_back(fit->residual);
    if(!residuals.empty()) {auto r=result("hard","coherent surface of revolution");r.rationale["revolution_residual"]=residuals.front();return r;}
    if(quad_fraction<.6 && regular_fraction<.3)return result("tri","triangle-heavy surface without a dominant primitive");
    return result("organic","curved quad surface without dominant panels or an axis");
}
ResolvedProfile resolve_auto_profile(const MeshSnapshot& m,PipelineOptions options,const Deadline& deadline) {
    if(options.profile!="auto" && options.profile!="custom")
        throw ContractError(ErrorCode::invalid_option,"profile","profile must be custom or auto; explicit settings belong in layout and boundary");
    if(options.profile!="auto" || !options.generate || deadline.expired())return {std::move(options),nullptr};
    const auto classified=asset_class(m);
    options=apply_classified_profile(std::move(options),classified.profile);
    return {std::move(options),{{"requested","auto"},{"selected",classified.profile},{"rationale",classified.rationale}}};
}
}
