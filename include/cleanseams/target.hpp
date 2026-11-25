// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
#include "cleanseams/layout.hpp"
#include <functional>
namespace cleanseams {
// core/target.py. Source island_band: null means no target; otherwise two
// integer-valued finite numbers with 1 <= low <= high (ContractError otherwise).
CLEANSEAMS_CORE_EXPORT std::optional<IslandBand> island_band(const nlohmann::json& value);
// Builds the evaluator for one trial's layout (the pipeline's default_evaluator or backend).
using EvaluatorFactory=std::function<ChartEvaluator(const LayoutOptions&)>;
struct TargetSearch {
    GeneratedLayout layout;      // the selected trial's seams, surface and report
    LayoutOptions selected_layout;
    nlohmann::json info;         // requested, attempts, selected, scale, cache_hits, status, reason
    std::uint64_t cache_hits=0;
};
// Source generate_for_target: at most five generate_seams trials with scaled
// distortion budgets, strict absorption and a scaled island floor, sharing one
// cache of measured chart maps (verdicts re-judged at each trial's budget). The
// best validated candidate wins; an unreachable band is reported, never faked.
// `input` carries the caller's cuts, locks, cost, deadline and pins for every trial.
CLEANSEAMS_CORE_EXPORT TargetSearch generate_for_target(const std::shared_ptr<const MeshSnapshot>&,const LayoutOptions&,IslandBand target,
    const EvaluatorFactory& factory,const SeedingInput& input={},LayoutTrace* trace=nullptr);
}
