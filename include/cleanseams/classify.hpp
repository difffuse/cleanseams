// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
#include "cleanseams/control.hpp"
#include "cleanseams/mesh.hpp"
#include "cleanseams/options.hpp"
#include <optional>
namespace cleanseams {
// core/classify.py: conservative geometry-based profile selection with evidence.
// Source _planar_coverage: surface share in up to eight nearly planar normal clusters.
struct PlanarCoverage {double coverage=0;std::uint64_t clusters=0;};
CLEANSEAMS_CORE_EXPORT PlanarCoverage planar_coverage(const MeshSnapshot&);
struct AssetClassInput {
    std::optional<std::vector<std::uint8_t>> features;       // evidence only: feature_length_fraction
    std::optional<std::vector<double>> axis_fit_residuals;  // precomputed revolution fits (none: fit here)
    bool plane_supplied=false;                              // evidence only: symmetry_supplied
};
// Source asset_class: "hard", "organic" or "tri" and the rationale dictionary.
struct AssetClass {std::string profile;nlohmann::json rationale;};
CLEANSEAMS_CORE_EXPORT AssetClass asset_class(const MeshSnapshot&,const AssetClassInput& input={});
// Source pipeline auto profile: with profile "auto", generation on and time left,
// classify and apply the layout/boundary settings; the report entry is
// {"requested": "auto", "selected", "rationale"} (null when not applied).
struct ResolvedProfile {PipelineOptions options;nlohmann::json report;};
CLEANSEAMS_CORE_EXPORT ResolvedProfile resolve_auto_profile(const MeshSnapshot&,PipelineOptions options,const Deadline& deadline={});
}
