// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
#include "cleanseams/layout_postprocess.hpp"
namespace cleanseams::repair_detail {
// Twelve synchronous weighted rounds on actual cut sectors. The proposal is never
// adopted here: the caller validates it after requested storage conversion.
std::optional<std::vector<Vec2>> relax_selected_faces(const PipelineInput&,const CutSurface&,
    const std::function<bool()>& stopped,nlohmann::json& diagnostic);
struct LayoutRecovery {std::vector<Vec2> uv;EdgeMask seams;};
std::optional<LayoutRecovery> recover_selected_regions(const PipelineInput&,const std::shared_ptr<const CutSurface>&,
    std::span<const std::int8_t> face_orientation,const std::function<bool()>& stopped,nlohmann::json& diagnostic,
    const std::function<ChartSolveResult(const ChartSolveRequest&)>& solver={},bool isolate_only=false);
double face_weighted_p95(const QualitySnapshot&);
std::optional<LayoutRecovery> join_small_charts(const PipelineInput&,const std::shared_ptr<const CutSurface>&,
    const std::function<bool()>& stopped,const std::function<bool(const LayoutRecovery&,const std::shared_ptr<const CutSurface>&,
        const LayoutRecovery&,const std::shared_ptr<const CutSurface>&,std::span<const MeshIndex>)>& valid,nlohmann::json& diagnostic);
std::optional<LayoutRecovery> reconstruct_singletons(const PipelineInput&,const std::shared_ptr<const CutSurface>&,
    std::span<const std::int8_t> face_orientation,const std::function<bool()>& stopped,
    const std::function<bool(const LayoutRecovery&,const std::shared_ptr<const CutSurface>&,
        const LayoutRecovery&,const std::shared_ptr<const CutSurface>&,std::span<const MeshIndex>,MeshIndex)>& valid,nlohmann::json& diagnostic);
std::optional<LayoutRecovery> recover_fragments(const PipelineInput&,const std::shared_ptr<const CutSurface>&,
    std::span<const std::int8_t> face_orientation,const std::function<bool()>& stopped,
    const LayoutPostprocessServices&,
    const std::function<bool(const LayoutRecovery&,const std::shared_ptr<const CutSurface>&,std::span<const MeshIndex>)>& valid,
    nlohmann::json& diagnostic);
struct LocalRecovery {
    std::vector<MeshIndex> corners,isolated_faces;
    std::vector<Vec2> values;
    double anisotropy_p95=0; // ratio, stretch + 1; original face-area weights
};
std::optional<LocalRecovery> recover_mixed_chart(const RepairContext&,const RepairState&,MeshIndex,
    const Deadline& operation_deadline,const Deadline& local_deadline,nlohmann::json& diagnostic);
}
