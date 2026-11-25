// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
#include "cleanseams/chart_solver.hpp"
#include "cleanseams/control.hpp"
#include "cleanseams/grid.hpp"
#include <functional>
#include <set>
namespace cleanseams {
// core/uv_boundary.py: straighten or smooth a chart's UV contour arc by arc,
// pin the candidate contour, re-solve the interior and keep the first intensity
// that validates; otherwise return the incoming map unchanged.
struct ContourMeasurement {QualityReport report;std::uint64_t crossings=0;};
struct ContourArc {
    std::int64_t index=0;
    std::string mode; // straight, smooth, degenerate or pinned
    std::size_t n=0;
    double jag_before=0,jag_after=quality_nan;
    std::vector<std::size_t> positions; // loop positions; empty for grid sides
    std::vector<Vec2> target; // diagnostic: the arc's full-intensity target
};
struct ContourResult {
    MeshIndex chart=0;
    std::string status="kept",reason; // optimized | kept | protected | skipped
    double intensity=0;
    std::vector<ContourArc> arcs;
    std::optional<ContourMeasurement> before,after;
    std::optional<GridMap> grid; // dimensions of an accepted exact rectangle
    std::vector<std::string> flags; // loop anchor flags after arc anchors were added
    std::uint64_t grid_recognitions=0,interior_solves=0;
};
// Interior re-solve against pinned loop sectors: (chart corner ids, their UVs).
using InteriorSolver=std::function<std::pair<std::vector<MeshIndex>,std::vector<Vec2>>(const std::shared_ptr<const QualityGeometry>&,
    MeshIndex chart,std::span<const SectorPin> pins,std::span<const Vec2> uv)>;
// Candidate validation against the accepted state of the chart.
using ContourValidator=std::function<std::pair<bool,ContourMeasurement>(const std::shared_ptr<const QualityGeometry>&,
    MeshIndex chart,std::span<const Vec2> uv,const ContourMeasurement& baseline)>;
// Chart report and boundary crossings of a full per-corner UV array.
CLEANSEAMS_CORE_EXPORT ContourMeasurement measure_chart(const std::shared_ptr<const QualityGeometry>&,MeshIndex chart,std::span<const Vec2> uv);
CLEANSEAMS_CORE_EXPORT std::vector<std::int64_t> cut_degree(const CutSurface&);
CLEANSEAMS_CORE_EXPORT std::vector<std::string> loop_anchor_flags(const CutSurface&,std::span<const MeshIndex> loop,
    std::span<const std::int64_t> degree,const std::set<MeshIndex>& pinned_corners,const std::set<MeshIndex>& manual_anchors,
    const BoundaryOptions&);
CLEANSEAMS_CORE_EXPORT bool arc_is_loop_straight(const MeshSnapshot&,std::span<const MeshIndex> verts,const BoundaryOptions&);
CLEANSEAMS_CORE_EXPORT bool arc_near_chord(std::span<const Vec2>,double tol_edges);
CLEANSEAMS_CORE_EXPORT bool chord_is_clear(std::span<const Vec2> loop_uv,const Vec2& a,const Vec2& b,std::span<const std::size_t> arc_positions);
CLEANSEAMS_CORE_EXPORT std::optional<std::vector<Vec2>> straight_target(std::span<const Vec3> points,const Vec2& a,const Vec2& b);
CLEANSEAMS_CORE_EXPORT std::optional<std::vector<Vec2>> smooth_target(std::span<const Vec2>,double radius,double band);
// A pipeline run passes its own settings (Gate A F30): its memory budget, cancellation
// token and pin completion policy (PIN-01 applies to a re-solve with one pinned sector).
struct InteriorSolverSettings {
    NumericalPolicy policy=NumericalPolicy::scipy_compatible;
    std::uint64_t memory_budget=active_memory_budget();
    PinCompletion pin_completion=PinCompletion::preserve;
    CancellationToken cancel{}; // checked after every re-solve
};
// Source default_interior_solver: ARAP with `arap_iterations`, the sector pins, LSCM start.
CLEANSEAMS_CORE_EXPORT InteriorSolver default_interior_solver(const BoundaryOptions&,const InteriorSolverSettings& settings={});
// Source default_validator: valid map, P95 growth and absolute caps, max growth.
CLEANSEAMS_CORE_EXPORT ContourValidator default_validator(const BoundaryOptions&);
struct ContourRequest {
    BoundaryOptions options;
    std::set<MeshIndex> pinned_corners; // original corner IDs
    std::set<MeshIndex> manual_anchors; // original vertex IDs
    const InteriorSolver* solve=nullptr; // default_interior_solver when null
    const ContourValidator* validate=nullptr; // default_validator when null
    bool protect=false;
    // PIN-02: preserve keeps user pins on
    // sectors off the contour at their incoming UV through the re-solve;
    // legacy_010 reproduces the source, which frees them.
    PinCompletion pin_policy=PinCompletion::preserve;
    // GB performance: clean_all may clean the charts in parallel on the active worker pool
    // (without a deadline): every chart is cleaned from the incoming map and its corners are
    // committed in chart order. Only for an interior solver and validator that read and
    // write the cleaned chart alone and may run on worker threads (the native ones).
    bool parallel=false;
};
struct ContourOutcome {std::vector<Vec2> uv;ContourResult result;};
CLEANSEAMS_CORE_EXPORT ContourOutcome clean_chart(const std::shared_ptr<const QualityGeometry>&,std::span<const Vec2> uv,MeshIndex chart,
    const ContourRequest&);
struct CleanAllOutcome {std::vector<Vec2> uv;std::vector<ContourResult> charts;};
CLEANSEAMS_CORE_EXPORT CleanAllOutcome clean_all(const std::shared_ptr<const QualityGeometry>&,std::span<const Vec2> uv,const ContourRequest&,
    std::optional<std::span<const MeshIndex>> charts=std::nullopt,const Deadline& deadline={},
    const std::function<void(const nlohmann::json&)>& progress={});
}
