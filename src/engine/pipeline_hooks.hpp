// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
// Pipeline entry points owned by T14 (classification, target search) and T17
// (output density and orientation). Each forwards to the native module.
#include "cleanseams/pipeline.hpp"
namespace cleanseams::pipeline_hooks {
struct TargetGeneration {GeneratedLayout generated;nlohmann::json info;LayoutOptions layout;};
// Source target.island_band: validated inclusive (minimum, maximum), or none.
std::optional<IslandBand> island_band(const std::optional<IslandBand>&);
// Source target.generate_for_target with the pipeline's evaluator factory.
TargetGeneration generate_for_target(const std::shared_ptr<const MeshSnapshot>&,const PipelineOptions&,const IslandBand&,bool core,
    UnwrapBackend* backend,const SeedingInput& request);
// Source pipeline post-processing of rep["island_target"] after validation.
void finish_target(nlohmann::json& target,const IslandBand&,const RepairState&,const PipelineOptions&,const nlohmann::json& final_metrics);
// Source output.normalise_density / orient_islands: (map, info).
std::pair<std::vector<Vec2>,nlohmann::json> normalise_density(const std::shared_ptr<const QualityGeometry>&,std::span<const Vec2> uv,
    const PinMap& pins,std::span<const std::uint8_t> exempt);
std::pair<std::vector<Vec2>,nlohmann::json> orient_islands(const std::shared_ptr<const QualityGeometry>&,std::span<const Vec2> uv,const PinMap& pins);
}
