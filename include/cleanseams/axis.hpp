// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
#include "cleanseams/control.hpp"
#include "cleanseams/cut_surface.hpp"
namespace cleanseams {
// core/axis.py. A rotational tangent field omega x x + v fitted to (position,
// normal) samples; translation is eliminated before solving for rotation.
struct LinearComplexFit {
    Vec3 point{},direction{}; // direction: unit, largest-magnitude component positive
    double pitch=0,residual=0,scale=0;
    std::array<double,3> eigenvalues{}; // whitened reduced system, ascending
};
CLEANSEAMS_CORE_EXPORT std::optional<LinearComplexFit> fit_linear_complex(std::span<const Vec3> positions,
    std::span<const Vec3> normals,std::optional<std::span<const double>> weights=std::nullopt,
    double residual_limit=1e-3,double pitch_limit=0.02);
// A fit on two rounds of area-weighted neighbour smoothing, with per-vertex
// cylindrical coordinates about the fitted axis.
struct RevolutionFrame : LinearComplexFit {
    std::vector<double> theta,height,radius; // per mesh vertex
};
CLEANSEAMS_CORE_EXPORT std::optional<RevolutionFrame> revolution_frame(const MeshSnapshot&,
    std::optional<std::span<const MeshIndex>> faces=std::nullopt,double residual_limit=1e-3);
// Straight trace along regular fans until a cut is reached (source _trace).
struct TracePath {std::vector<MeshIndex> edges;bool complete=false;};
CLEANSEAMS_CORE_EXPORT TracePath trace_path(const MeshSnapshot&,std::span<const std::uint8_t> cut,MeshIndex vertex,MeshIndex edge);
struct DirectedSeed {
    TransferIndex vertex=0,edge=0;
    auto operator<=>(const DirectedSeed&) const=default;
};
// Traces only the requested (vertex, edge) directions; incomplete paths are discarded.
CLEANSEAMS_CORE_EXPORT std::vector<std::uint8_t> trace_directed(const MeshSnapshot&,std::span<const std::uint8_t> cut,
    std::span<const DirectedSeed> seeds);
CLEANSEAMS_CORE_EXPORT std::vector<std::uint8_t> shoulder_seed(const CutSurface&,MeshIndex chart,const RevolutionFrame&);
CLEANSEAMS_CORE_EXPORT std::optional<DirectedSeed> generatrix_seed(const CutSurface&,MeshIndex chart,const RevolutionFrame&,
    std::optional<std::span<const double>> cost=std::nullopt,std::optional<Vec3> hide_dir=std::nullopt);
// Geometry-only fits are reused within one call when the face set and limit match.
struct FitStats {std::uint64_t fits=0,reused=0;};
struct RevolutionRecord {
    DirectedSeed seed;
    std::uint64_t edges=0;
    double residual=0,pitch=0;
    Vec3 direction{},point{};
};
struct RevolutionCuts {
    std::vector<std::uint8_t> cut,shoulders;
    std::vector<RevolutionRecord> records;
    FitStats stats;
};
CLEANSEAMS_CORE_EXPORT RevolutionCuts revolution_cuts(const std::shared_ptr<const MeshSnapshot>&,std::span<const std::uint8_t> cut,
    double residual_limit=0.04,std::optional<std::span<const double>> cost=std::nullopt,std::optional<Vec3> hide_dir=std::nullopt,
    const Deadline& deadline={});
struct LatheRecord {
    std::uint64_t charts=0,rows=0;
    std::vector<std::int64_t> rim_rows;
    double residual=0;
    Vec3 direction{};
};
struct LatheProfileCuts {
    std::vector<std::uint8_t> cut,protected_rims;
    std::vector<LatheRecord> records;
    std::vector<MeshIndex> handled_faces;
    FitStats stats;
};
CLEANSEAMS_CORE_EXPORT LatheProfileCuts lathe_profile_cuts(const std::shared_ptr<const MeshSnapshot>&,std::span<const std::uint8_t> initial_cut,
    std::optional<std::span<const double>> cost=std::nullopt,std::optional<Vec3> hide_dir=std::nullopt,const Deadline& deadline={});
}
