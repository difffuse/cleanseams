// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
#include "cleanseams/control.hpp"
#include "cleanseams/cut_surface.hpp"
namespace cleanseams {
// core/panel_layout.py. Shortest existing-edge path joining the chart's shortest
// boundary loop to another loop through distinct UV sectors; empty when none.
CLEANSEAMS_CORE_EXPORT std::optional<std::vector<MeshIndex>> join_boundary_loops(const CutSurface&,MeshIndex chart,
    std::optional<std::span<const double>> cost=std::nullopt);
struct CadPanelCuts {
    std::vector<std::uint8_t> cut,rims;
    std::uint64_t panels=0,hole_cut_edges=0,bevel_bands=0;
};
// Planar CAD panels seeded from coherent 3-degree groups, bevels reassigned to
// neighbouring panels, holes opened by boundary joins that keep the chart count.
// Empty whenever the recogniser declines the mesh.
CLEANSEAMS_CORE_EXPORT std::optional<CadPanelCuts> cad_panel_cuts(const std::shared_ptr<const MeshSnapshot>&,
    std::span<const std::uint8_t> initial_cut,std::optional<std::span<const double>> cost=std::nullopt,
    const Deadline& deadline={},bool allow_triangles=false);
}
