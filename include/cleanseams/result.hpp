// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
#include "cleanseams/solver.hpp"
namespace cleanseams {
struct ExemptFace {
    TransferIndex face=0;
    std::string reason;
};
struct UnavailableValue {
    std::string path;
    std::string reason;
};
struct RunResult {
    std::uint64_t version=1;
    bool accepted=false; // compatibility decision; never substitutes for independent flags
    std::string reason;
    bool topology_valid=false;
    bool uv_valid=false;
    bool budget_satisfied=false;
    bool partial=false;
    bool cancelled=false;
    bool timed_out=false;
    std::string geometry_fingerprint;
    Revision state_revision=0;
    std::vector<std::uint8_t> seams;
    std::optional<std::vector<Vec2>> uv;
    std::vector<CornerPin> pins;
    std::vector<ExemptFace> exempt_faces;
    std::vector<TransferIndex> invalid_charts;
    std::vector<TransferIndex> changed_charts;
    PipelineOptions resolved_options;
    std::string backend="superlu";
    std::vector<ErrorInfo> diagnostics;
    nlohmann::json report=nlohmann::json::object();
    nlohmann::json decision_trace=nlohmann::json::array();
    std::vector<UnavailableValue> unavailable;
    // Owned original-input mappings. Chart numbers are local to this result.
    std::vector<TransferIndex> chart_of_face,chart_of_corner,changed_faces;
    std::vector<EdgeEndpoints> edge_vertices;
};
struct ProgressEvent {
    std::string stage;
    std::uint64_t done=0;
    std::uint64_t total=0;
    double elapsed_s=0.;
    std::optional<double> eta_s;
    nlohmann::json counters=nlohmann::json::object();
};
CLEANSEAMS_CORE_EXPORT RunResult run_result_from_json(const nlohmann::json&);
CLEANSEAMS_CORE_EXPORT nlohmann::json run_result_json(const RunResult&);
}
