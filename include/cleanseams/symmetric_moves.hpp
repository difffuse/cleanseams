// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
#include "cleanseams/proposals.hpp"
#include "cleanseams/symmetry.hpp"
namespace cleanseams {
// ============================================================================
// Phase B operational symmetry (Section 19, T38), policy enhanced_v1 only.
//
// The user selects the policy explicitly (PipelineOptions::enhanced.symmetry.mode; there is no
// environment variable): `asymmetric` (view-optimal, moves are never paired), `strict` or
// `approximate` (metric-only diagnostics). In strict mode the transaction pairs every proposal
// with its mirror image before anything is solved:
//
//   exact correspondence on the affected patch -> mirror the added/removed edges, endpoint
//   moves and allowed endpoints -> merge the halves (self-mirrored and overlapping edges once)
//   -> conflicts of the mirror half (locks, required, forbidden, anchors, pins) -> one delta
//   evaluated on every affected chart -> commit only when every change has its mirror half.
//
// A conflict refuses the pair with its reason (reject_pair) or, when the user selected
// disable_for_proposal, evaluates the unpaired proposal with an explicit label. Painted
// constraints, anchors and pin coordinates are never mirrored or overwritten.
// ============================================================================
constexpr const char* symmetric_moves_version="symmetric_moves.v1";
// The mirror correspondence of a mesh under a plane (T14 mirror maps at the counterpart
// tolerance), validated element by element for strict pairing: a vertex is exact when its
// counterpart is an involution partner whose position is the reflection within the geometry
// tolerance; an edge when both vertices are exact, its image is an edge mapped back to it and
// both are topology facts or neither; a face when every vertex is exact and its image is a
// face mapped back to it.
struct CLEANSEAMS_CORE_EXPORT StrictMirrorMap {
    // exact (every element), partial (some exact faces), approximate (counterparts only within
    // the counterpart tolerance), noncombinatorial (no complete correspondence), no_plane.
    std::string status;
    std::optional<MirrorPlane> plane;
    std::string plane_source;                  // explicit, detected or none
    std::string geometry_key;                  // the mesh geometry the map belongs to
    SymmetryOptions options;                   // the options the map was built with
    double diagonal=0,map_radius=0,tolerance=0; // bounding-box diagonal; absolute radii
    double max_deviation=0;                    // over matched vertices: |reflect(p) - p(counterpart)|
    MirrorMaps maps;                           // -1 where there is no counterpart
    std::vector<std::uint8_t> vertex_exact,edge_exact,face_exact;
    std::vector<std::int8_t> vertex_side;      // -1, 0 (on the plane within the tolerance), +1
    std::uint64_t matched_vertices=0,exact_vertices=0,exact_edges=0,exact_faces=0,self_mirrored_edges=0;
    bool usable() const {return exact_faces>0;}
    // Image of an exact edge (nullopt otherwise).
    std::optional<MeshIndex> mirror_edge(MeshIndex e) const;
    std::optional<MeshIndex> mirror_vertex(MeshIndex v) const;
    nlohmann::json to_json() const;
};
CLEANSEAMS_CORE_EXPORT StrictMirrorMap strict_mirror_map(const MeshSnapshot&,const SymmetryOptions&);
// The result of pairing one proposal. `status`: self_mirrored (the delta is its own image),
// mirror_complete (the mirror half changes nothing the delta does not already change),
// paired (the mirror half added changes) or conflict (`reason`, exact IDs in `detail`):
// no_strict_map, approximate_correspondence, noncombinatorial, overlapping_halves,
// asymmetric_lock, asymmetric_required, asymmetric_forbidden (the paired cut set violates a
// constraint the original half alone does not), asymmetric_anchor (the mirror half moves an
// anchored endpoint), asymmetric_pin (the pair changes the cuts at a vertex pinned on one side
// only, or the mirror half glues differently pinned sectors), unrepresentable_edge (left to the
// transaction, which names it).
struct CLEANSEAMS_CORE_EXPORT SymmetricPairing {
    std::string status,reason;
    Proposal proposal;                          // the paired proposal (a conflict: the original)
    nlohmann::json detail=nlohmann::json::object();
    bool ok() const {return status!="conflict";}
    nlohmann::json to_json() const;
};
// Pairs `proposal` (a delta against `parent`) under `map`; the constraints are only read.
// Pairing a paired proposal returns it unchanged (idempotent).
CLEANSEAMS_CORE_EXPORT SymmetricPairing pair_proposal(const Proposal&,const RunState& parent,const HardConstraints&,const StrictMirrorMap&);
// Null when every declared change of `proposal` has its mirror half realized in `cuts` and every
// endpoint move its mirror move; otherwise the missing halves (edge and vertex IDs).
CLEANSEAMS_CORE_EXPORT nlohmann::json missing_mirror_half(const Proposal&,std::span<const std::uint8_t> cuts,const StrictMirrorMap&);
// The pass's mirror map (built once per geometry and symmetry options, on first use): null in
// the asymmetric mode. Throws invalid_option for an unknown mode or conflict policy.
CLEANSEAMS_CORE_EXPORT std::shared_ptr<const StrictMirrorMap> context_symmetry(ProposalContext&);
// True when the pass pairs proposals (mode strict).
CLEANSEAMS_CORE_EXPORT bool strict_symmetry(const ProposalContext&);
// Metric-only symmetry measurements of a cut state: exact seam pairing (length share of the
// internal seams on the exact correspondence whose image is cut), the T14 seam-symmetry
// score under the map's plane, unpaired edges, and the nearest-surface deviation of the
// reflected vertices (approximate diagnostics; no hard edge equality is promised).
CLEANSEAMS_CORE_EXPORT nlohmann::json symmetry_metrics(const MeshSnapshot&,const StrictMirrorMap&,std::span<const std::uint8_t> cuts);
// Report record of a pass: mode, conflict policy, map, per-generator pairing counts, and the
// metrics of candidate zero and of the retained state. Null in the asymmetric mode.
CLEANSEAMS_CORE_EXPORT nlohmann::json symmetry_summary(ProposalContext&,const RunState& origin,const RunState& retained);
// Section 21 "strict symmetry pairing where requested": in strict mode, up to
// enhanced.symmetry.max_candidates symmetrization alternatives of the incumbent's unpaired
// seams on the exact correspondence: `closure` (the images of the unpaired seams added),
// `mirror_negative` and `mirror_positive` (one side's seams copied onto the other). Edges that
// are hard-constrained, or whose image is, keep their state (asymmetric constraints are never
// mirrored); non-disk charts are reopened along admissible paths whose images are added too.
// In the approximate mode it reports metric-only diagnostics and proposes nothing.
CLEANSEAMS_CORE_EXPORT ProposalGenerator symmetry_generator();
}
