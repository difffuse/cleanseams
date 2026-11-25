// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
#include "cleanseams/quality.hpp"
#include <map>
namespace cleanseams {
// Source core/output.py: per-island similarity transforms applied after unwrapping.
// Maps are dense per original corner (nc rows, NaN where unavailable). Charts holding
// a pinned corner are never moved.
struct CLEANSEAMS_CORE_EXPORT DensityNormalization {
    std::vector<Vec2> uv;
    double density_spread_before=1.0,target_density=1.0;
    std::vector<MeshIndex> protected_charts,density_outliers;
    nlohmann::json info() const; // the source's info dict (NaN stays a JSON float)
};
// normalise_density: scale every chart with a finite positive density (quality
// chart_density, degenerate and `exempt` faces excluded) about its corner mean to the
// mesh's area-weighted target density.
CLEANSEAMS_CORE_EXPORT DensityNormalization normalise_density(const std::shared_ptr<const CutSurface>&,std::span<const Vec2> uv,
    std::span<const MeshIndex> pinned_corners={},std::span<const std::uint8_t> exempt={});
struct CLEANSEAMS_CORE_EXPORT IslandOrientation {
    std::vector<Vec2> uv;
    std::map<MeshIndex,double> angles_radians;
    std::vector<MeshIndex> protected_charts;
    nlohmann::json info() const;
};
// Source output._hull: monotone chain over the distinct points (sorted), collinear
// points dropped; two or fewer distinct points are returned as they are.
CLEANSEAMS_CORE_EXPORT std::vector<Vec2> output_hull(std::span<const Vec2> points);
// orient_islands: rotate each unpinned chart with a finite boundary about its corner
// mean so that its minimum-area hull rectangle is axis aligned and wider than tall.
CLEANSEAMS_CORE_EXPORT IslandOrientation orient_islands(const CutSurface&,std::span<const Vec2> uv,std::span<const MeshIndex> pinned_corners={});
}
