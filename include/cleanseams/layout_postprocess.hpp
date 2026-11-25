// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
#include "cleanseams/pipeline.hpp"
namespace cleanseams {
// Optional native numerical service for synchronous C++ callers and deterministic
// solver-failure tests. Session/wire/binding calls use the native LSCM default
// and opt into fragment host services through SessionConfig.host_capabilities.
// Supplied proposals pass the same mapping, quality and storage validation.
struct LayoutPostprocessServices {
    std::function<ChartSolveResult(const ChartSolveRequest&)> region_solver;
    // Optional synchronous host operations. Projection may propose new cuts; unwrap
    // receives a scratch mesh split at those native-derived cuts. Neither callback
    // decides whether a proposal replaces the layout. nullopt means host failure.
    // The original-corner map and cumulative stop predicate are borrowed only for
    // the duration of the call; queue waits must honor the predicate.
    std::function<std::optional<std::vector<Vec2>>(const MeshInput&,std::span<const TransferIndex>,const std::function<bool()>&)> region_proposal;
    std::function<std::optional<std::vector<Vec2>>(const HostSolveRequest&,const std::function<bool()>&)> fragment_unwrap;
    std::string fragment_unwrap_method="ANGLE_BASED";
};
// Shared synchronous contract check for direct native and Session callers.
CLEANSEAMS_CORE_EXPORT void check_layout_operation(const PipelineInput&);
// Processes an existing layout; never generates seams or invokes a global unwrap.
// Optional repair stages are added here under repair_v1, using a single phase budget.
CLEANSEAMS_CORE_EXPORT PipelineResult postprocess_layout(const PipelineInput&,const ProgressSink& progress={},const Control& control={});
CLEANSEAMS_CORE_EXPORT PipelineResult postprocess_layout(const PipelineInput&,const ProgressSink&,const Control&,const LayoutPostprocessServices&);
}
