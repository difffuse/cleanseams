// SPDX-License-Identifier: GPL-3.0-or-later
// GB (M3): seam_objective.v2, the significance margin of the common ordering. The only sampled
// term of J is the raw all-views exposure V; its sampling uncertainty between two cut states is
// the linearized variance of the vertex ratio estimators over the shared direction samples
// (exposure_uncertainty.v1; derivation in include/cleanseams/proposals.hpp and
// docs/ENHANCED_SEAMS.md). Used by policy enhanced_v1 only.
#include "cleanseams/proposals.hpp"
#include "cleanseams/fields.hpp"
#include "cleanseams/errors.hpp"
#include "cleanseams/detail/mesh_math.hpp"
#include <algorithm>
#include <cmath>
#include <map>
namespace cleanseams {
namespace {
using Json=nlohmann::json;
}
ExposureUncertainty::ExposureUncertainty(std::shared_ptr<const MeshSnapshot> mesh,const ExposureField& field):mesh_(std::move(mesh)) {
    if(!mesh_)throw ContractError(ErrorCode::invalid_input,"enhanced.exposure_uncertainty","missing mesh");
    const auto nv=mesh_->input().positions.size();const auto K=field.directions.size();
    if(!field.visibility || field.vertex.size()!=nv || field.visibility->vertex.size()!=nv*K || field.visibility->vertex_weight.size()!=nv*K)
        throw ContractError(ErrorCode::invalid_input,"enhanced.exposure_uncertainty","the exposure field does not belong to this mesh");
    directions_=field.directions;words_=(K+63)/64;
    visible_.assign(nv*words_,0);exposure_=field.vertex;weight_.assign(nv,0.0);
    const auto& vis=*field.visibility;
    for(std::size_t v=0;v<nv;++v) {
        double total=0;
        for(std::size_t d=0;d<K;++d) {
            total+=static_cast<double>(vis.vertex_weight[v*K+d]);
            if(vis.vertex[v*K+d])visible_[v*words_+d/64]|=std::uint64_t{1}<<(d%64);
        }
        weight_[v]=total;
    }
}
double ExposureUncertainty::delta_v_sigma(std::span<const std::uint8_t> a,std::span<const std::uint8_t> b,double area) const {
    if(!usable())return 0.0;
    const auto& t=mesh_->tables();const auto& g=mesh_->geometry();const auto ne=t.edge_vertices.size();
    if(a.size()!=ne || b.size()!=ne)throw ContractError(ErrorCode::invalid_input,"enhanced.exposure_uncertainty","cut masks do not match the mesh edges");
    const double root=area>0?std::sqrt(area):1.0;
    // a_v: half the signed, normalized length of every internal edge one state has and the other
    // lacks (an edge's exposure is the mean of its vertices').
    std::map<MeshIndex,double> coefficient;
    for(std::size_t e=0;e<ne;++e) {
        if((a[e]!=0)==(b[e]!=0) || t.edge_boundary[e])continue;
        const double c=(a[e]?0.5:-0.5)*g.edge_length[e]/root;
        for(auto v:t.edge_vertices[e])coefficient[v]+=c;
    }
    const auto K=directions_.size();
    std::vector<double> G(K,0.0),w(K,0.0);
    const auto& corner_face=t.corner_face;
    for(const auto& [v,c]:coefficient) {
        if(c==0 || !(weight_[v]>0))continue;
        // The visibility's weights: float32 of the largest front-facing cosine of the vertex's
        // faces (zbuffer_visibility), recomputed instead of stored.
        for(std::size_t d=0;d<K;++d) {
            double cosv=0;
            for(auto corner:t.vert_corners[v])cosv=std::max(cosv,std::max(detail::dot(g.face_normal[corner_face[corner]],directions_[d]),0.0));
            w[d]=static_cast<double>(static_cast<float>(cosv));
        }
        const double x=exposure_[v],W=weight_[v];
        for(std::size_t d=0;d<K;++d) {
            const double seen=(visible_[v*words_+d/64]>>(d%64))&1U?1.0:0.0;
            G[d]+=c*w[d]*(seen-x)/W;
        }
    }
    double sum=0;for(double x:G)sum+=x*x;
    const double k=static_cast<double>(K);
    return std::sqrt(sum*k/(k-1.0));
}
nlohmann::json ExposureUncertainty::to_json() const {
    std::uint64_t mixed=0;
    if(mesh_)for(std::size_t v=0;v<exposure_.size();++v)mixed+=weight_[v]>0 && exposure_[v]>0 && exposure_[v]<1;
    return {{"version",exposure_uncertainty_version},{"directions",directions()},{"usable",usable()},{"vertices_with_mixed_visibility",mixed},
            {"model","linearized ratio-estimator variance over the shared direction samples"},{"z",seam_objective_v2_significance}};
}
bool significance_active(const ProposalContext& ctx) {
    const auto& o=ctx.options.enhanced.objective;
    return o.version==seam_objective_v2 && ctx.exposure_uncertainty && ctx.exposure_uncertainty->usable() && !ctx.raw_exposure.empty() && o.w_exposure>0;
}
double objective_margin(const ProposalContext& ctx,std::span<const std::uint8_t> a,std::span<const std::uint8_t> b,double area) {
    if(!significance_active(ctx))return 0.0;
    const auto& o=ctx.options.enhanced.objective;
    return seam_objective_v2_significance*o.w_exposure*ctx.exposure_uncertainty->delta_v_sigma(a,b,area);
}
nlohmann::json significance_record(const ProposalContext& ctx,double challenger_J,double incumbent_J,double margin,bool ranks_first) {
    const auto& o=ctx.options.enhanced.objective;
    const double sigma=o.w_exposure>0?margin/(seam_objective_v2_significance*o.w_exposure):0.0;
    const double delta=challenger_J-incumbent_J;
    return {{"version",o.version},{"model",exposure_uncertainty_version},{"z",seam_objective_v2_significance},{"delta_J",delta},
            {"sigma_V",sigma},{"uncertainty",margin},{"within_uncertainty",std::abs(delta)<=margin},{"ranks_before_incumbent",ranks_first}};
}
// GB (visual review V2): the practical-effect floor. U models only the exposure sampling, so where
// the sampling noise is tiny (sigma_V ~ 1e-4 on a cylinder) a 0.1 % change of J was significant
// and adopted although no one can see it (a slit sliding two columns). The floor asks a change
// that adds seam to improve J by a practical fraction of it; removing seam (spur removal, a
// subset of the incumbent's cuts) is a simplification and keeps U alone.
double practical_floor(const ProposalContext& ctx,std::span<const std::uint8_t> challenger,std::span<const std::uint8_t> incumbent,double incumbent_J) {
    if(!significance_active(ctx))return 0.0;
    const double fraction=ctx.options.enhanced.objective.significance.practical_floor;
    if(!(fraction>0) || !std::isfinite(incumbent_J))return 0.0;
    if(challenger.size()!=incumbent.size())throw ContractError(ErrorCode::invalid_input,"enhanced.significance","cut masks do not match");
    for(std::size_t e=0;e<challenger.size();++e)if(challenger[e] && !incumbent[e])return fraction*std::abs(incumbent_J);
    return 0.0;
}
ObjectiveSignificance objective_significance(const ProposalContext& ctx,std::span<const std::uint8_t> challenger,std::span<const std::uint8_t> incumbent,
                                             double area,double incumbent_J) {
    return {objective_margin(ctx,challenger,incumbent,area),practical_floor(ctx,challenger,incumbent,incumbent_J)};
}
double objective_margin(const ProposalContext& ctx,std::span<const std::uint8_t> challenger,std::span<const std::uint8_t> incumbent,double area,double incumbent_J) {
    return objective_significance(ctx,challenger,incumbent,area,incumbent_J).margin();
}
nlohmann::json significance_record(const ProposalContext& ctx,double challenger_J,double incumbent_J,const ObjectiveSignificance& s,bool ranks_first) {
    auto out=significance_record(ctx,challenger_J,incumbent_J,s.uncertainty,ranks_first);
    const double delta=std::abs(challenger_J-incumbent_J);
    out["practical_floor"]=s.floor;out["practical_floor_fraction"]=ctx.options.enhanced.objective.significance.practical_floor;
    out["margin"]=s.margin();out["within_margin"]=delta<=s.margin();
    out["below_practical_floor"]=challenger_J<incumbent_J && delta>s.uncertainty && delta<=s.floor; // an improvement only the floor refused
    return out;
}
}
