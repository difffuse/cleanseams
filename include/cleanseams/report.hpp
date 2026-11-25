// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
#include "cleanseams/pipeline.hpp"
#include "cleanseams/metrics.hpp"
#include "cleanseams/result.hpp"
namespace cleanseams {
// Human summaries of a pipeline result (moved from the add-on's ops.py so no selection
// or formatting policy remains in UI code) and the finite JSON form of reports.
// Source ops.repairs_message: one clause on what the repair ladder did.
CLEANSEAMS_CORE_EXPORT std::string repairs_message(const nlohmann::json& report);
// Source ops.refusal_message: "reason: explanation (islands ...)".
CLEANSEAMS_CORE_EXPORT std::string refusal_message(const std::string& reason,std::span<const MeshIndex> invalid_charts);
// Source ops.summarize: the one-line status of a run.
CLEANSEAMS_CORE_EXPORT std::string summarize(const PipelineResult&);
// T41: the actionable lines a host panel or a command-line user reads for a run, derived from
// its finite report alone (policy, the enhancement outcome and why proposals were refused,
// budgets, required relief slits and painted constraints, the solver backend, the symmetry
// and distance policies, stops and time budgets, strict-validity defects, before/after seam
// measurements). Each line is {"level": info|warning|error, "topic": ..., "text": ...}. A
// native_baseline or legacy_010 report gets its policy line and the generic run lines only.
CLEANSEAMS_CORE_EXPORT nlohmann::json run_diagnostics(const nlohmann::json& report);
// T41: the lines of a packing step, never a generic success: `engine` is the host's packer
// choice ("blender"). The native engine does not verify host packing.
CLEANSEAMS_CORE_EXPORT nlohmann::json packing_diagnostics(const nlohmann::json& packing_report,const std::string& engine);
// Source ops._json_safe: every non-finite number becomes null (no NaN/Infinity tokens).
CLEANSEAMS_CORE_EXPORT nlohmann::json finite_json(const nlohmann::json&);
// The one conversion of a pipeline result into the public RunResult, shared by
// sessions, the Python bindings and the CLI: independent validity flags (a missing
// final gate never counts as within budget), exempt faces with their kind, the finite
// report with its summary line (and refusal text when refused).
CLEANSEAMS_CORE_EXPORT RunResult run_result(const PipelineResult&,const PipelineOptions&,std::span<const CornerPin> pins,const std::string& backend);
// native_symmetry_service() (seam_report's default symmetry part) is declared in metrics.hpp.
}
