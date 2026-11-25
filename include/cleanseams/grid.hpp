// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
#include "cleanseams/cut_surface.hpp"
namespace cleanseams {
using GridCoordinate=std::array<std::int64_t,2>;
// A recognized rectangular quad lattice on cut sectors (source grid_topology).
// It is tied to the chart key it was recognized on and reused by mapping and
// contour cleanup instead of recognizing the chart again.
struct GridDescriptor {
    MeshIndex chart=0;
    std::string chart_key; // ChartKey::digest() of the recognized chart
    std::vector<MeshIndex> faces,corners; // chart faces; four corners per face in face order
    std::vector<std::array<MeshIndex,4>> face_sectors;
    std::vector<MeshIndex> sectors; // sorted lattice sectors
    std::vector<GridCoordinate> ij; // lattice (i, j) of each entry of `sectors`
    std::uint64_t nu=0,nv=0;
    std::array<std::vector<MeshIndex>,4> arcs; // boundary corners anchor to anchor, closing anchor included
    GridCoordinate coordinate(MeshIndex sector) const;
};
struct GridStats {std::uint64_t recognitions=0;};
CLEANSEAMS_CORE_EXPORT std::optional<GridDescriptor> grid_topology(const CutSurface&,MeshIndex chart,GridStats* stats=nullptr);
CLEANSEAMS_CORE_EXPORT bool is_grid_chart(const CutSurface&,MeshIndex chart);
struct GridMap {
    std::vector<MeshIndex> corner_ids; // the descriptor's corners
    std::vector<Vec2> uv;
    double lu=0,lv=0;
    std::array<std::uint64_t,2> grid{};
};
// Source grid_uv: arc-length Coons coordinates of the lattice. Empty when a side
// has zero length. In-chart pins must equal their mapped coordinates exactly;
// otherwise the request is a constraint_conflict and the caller picks a validated
// fallback. The descriptor must belong to this cut surface's chart.
CLEANSEAMS_CORE_EXPORT std::optional<GridMap> grid_uv(const GridDescriptor&,const CutSurface&,std::span<const CornerPin> pins={});
}
