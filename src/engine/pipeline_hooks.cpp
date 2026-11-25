// SPDX-License-Identifier: GPL-3.0-or-later
#include "pipeline_hooks.hpp"
#include "cleanseams/errors.hpp"
#include "cleanseams/classify.hpp"
#include "cleanseams/output.hpp"
#include "cleanseams/target.hpp"
namespace cleanseams::pipeline_hooks {
std::optional<IslandBand> island_band(const std::optional<IslandBand>& band) {
    return band?cleanseams::island_band(nlohmann::json(*band)):std::nullopt;
}
TargetGeneration generate_for_target(const std::shared_ptr<const MeshSnapshot>& mesh,const PipelineOptions& options,const IslandBand& band,bool core,
                                     UnwrapBackend* backend,const SeedingInput& request) {
    const EvaluatorFactory factory=[&](const LayoutOptions& layout) {
        if(core) {EvaluatorSettings settings;settings.pin_completion=pin_policy(options);settings.cancel=request.cancel;return default_evaluator(mesh,layout,settings);}
        auto trial=options;trial.layout=layout;return backend_evaluator(mesh,*backend,trial);
    };
    auto search=cleanseams::generate_for_target(mesh,options.layout,band,factory,request);
    return {std::move(search.layout),std::move(search.info),std::move(search.selected_layout)};
}
void finish_target(nlohmann::json& target,const IslandBand& band,const RepairState& state,const PipelineOptions& o,const nlohmann::json& final_metrics) {
    target["generation_status"]=target.at("status");
    const auto& ct=state.cs->tables();target["final_charts"]=ct.n_charts;
    std::vector<MeshIndex> bad;const auto geometry=state.geometry();
    QualityPolicy policy;for(std::size_t f=0;f<state.exempt.size();++f)if(state.exempt[f])policy.exempt_faces.push_back(static_cast<TransferIndex>(f));
    BoundaryCrossingOptions ignore;ignore.ignored_faces=policy.exempt_faces;
    for(MeshIndex ch=0;ch<ct.n_charts;++ch) {
        const auto faces=ct.chart_faces[ch];
        if(std::all_of(faces.begin(),faces.end(),[&](auto f) {return state.exempt[f]!=0;}))continue;
        const std::array<MeshIndex,1> one{ch};
        const auto report=QualitySnapshot(geometry,CornerUvView{state.uv},policy,std::span<const MeshIndex>(one)).report();
        const auto crossings=boundary_crossings(geometry->cut(),CornerUvView{state.uv},ignore,std::span<const MeshIndex>(one)).counts[ch];
        if(!within_budget(report,crossings,o.layout))bad.push_back(ch);
    }
    target["final_over_budget_charts"]=bad;
    const auto n=static_cast<std::int64_t>(ct.n_charts);
    if(!bad.empty() || n<band[0] || n>band[1] || !final_metrics.at("valid").get<bool>()) {
        target["status"]="unreached";target["reason"]=!bad.empty()?"final distortion budget":"final island count";
    }
}
std::pair<std::vector<Vec2>,nlohmann::json> normalise_density(const std::shared_ptr<const QualityGeometry>& geometry,std::span<const Vec2> uv,
                                                              const PinMap& pins,std::span<const std::uint8_t> exempt) {
    std::vector<MeshIndex> corners;for(const auto& [c,p]:pins)corners.push_back(c);
    const auto cut=std::shared_ptr<const CutSurface>(geometry,&geometry->cut());
    auto result=cleanseams::normalise_density(cut,uv,corners,exempt);return {std::move(result.uv),result.info()};
}
std::pair<std::vector<Vec2>,nlohmann::json> orient_islands(const std::shared_ptr<const QualityGeometry>& geometry,std::span<const Vec2> uv,const PinMap& pins) {
    std::vector<MeshIndex> corners;for(const auto& [c,p]:pins)corners.push_back(c);
    auto result=cleanseams::orient_islands(geometry->cut(),uv,corners);return {std::move(result.uv),result.info()};
}
}
