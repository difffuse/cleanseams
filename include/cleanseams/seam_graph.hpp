// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
#include "cleanseams/cut_surface.hpp"
namespace cleanseams {
using ChartPair=std::optional<std::array<MeshIndex,2>>;
struct CLEANSEAMS_CORE_EXPORT Chain {
    MeshIndex id=0;
    // The frozen source repeats the first vertex at the end of cycles too.
    std::vector<MeshIndex> verts,edges;
    bool cycle=false,locked=false;
    // SeamGraph emits pairs (including self-pairs); routing can construct a
    // cycle part with a deduplicated single-chart tuple, as the source does.
    std::optional<std::vector<MeshIndex>> charts;
    std::string reason;
    nlohmann::json tags=nlohmann::json::object();
    bool free_end=false;
    std::size_t length() const {return edges.size();}
    double edge_length(const MeshSnapshot&) const;
};
struct AnchorProvenance {
    bool artist_fixed=false,locked_incident=false,pinned_sector=false;
    bool physical_boundary=false,generated_junction=false,temporary_cycle_anchor=false;
    bool seam_end=false,unsupported_fan=false,flanking_change=false;
    // Preserve distinct cut-side identities even at one degree-two vertex.
    std::vector<MeshIndex> sectors,pinned_sectors;
};
struct SeamTables {
    std::vector<std::uint8_t> seams,locked,anchor;
    std::vector<MeshIndex> manual_anchors;
    std::vector<Offset> degree;
    // Indexed by canonical edge, meaningful for seam edges only.
    std::vector<ChartPair> edge_charts;
    IndexRows vert_seam_edges;
    std::vector<AnchorProvenance> provenance;
};
class CLEANSEAMS_CORE_EXPORT SeamGraph {
public:
    // Empty locked span means no locks. Pins identify original corners; their
    // coordinates/conflicts are solver concerns. Provenance never moves anchors.
    SeamGraph(std::shared_ptr<const MeshSnapshot>,std::span<const std::uint8_t> seams,
        std::span<const std::uint8_t> locked={},std::span<const TransferIndex> anchors={},
        std::shared_ptr<const CutSurface> cut={},std::span<const TransferIndex> pinned_corners={},
        std::uint64_t memory_budget=active_memory_budget());
    const CutSurface& cut() const {return *cut_;}
    const SeamTables& tables() const {return data_;}
    const std::vector<Chain>& chains() const {return chains_;}
    // Stable IDs into chains(), retaining source order without copying chains.
    std::vector<MeshIndex> editable_chains() const;
    nlohmann::json summary() const;
private:
    std::shared_ptr<const CutSurface> cut_;
    SeamTables data_;
    std::vector<Chain> chains_;
    void build();
};
}
