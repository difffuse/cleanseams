// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
// Phase B sliding seam endpoints on tube-like patches (Section 13; T30, T31). Used by
// policy enhanced_v1 only; legacy_010 and native_baseline never reach any of it.
//
//   boundary-side provenance (CutSurface)      -> which chart boundary sides are physical
//                                                 rims, forced cuts or seam sides
//   underlying patch (slit glued, rims)        -> PhysicalRimId of each rim of the annulus
//   endpoint eligibility (constraints)         -> explicit admissible endpoint vertex/sector
//                                                 sets on two distinct rims, or a reason
//   multi-endpoint routing over a patch corridor -> direction-aware multi-source/multi-goal
//                                                 search with displacement/terminal costs
//   whole-slit relocation (T31)                -> up to four deduplicated routes per eligible
//                                                 slit as one-delta proposals of the common
//                                                 T28 transaction
#include "cleanseams/proposals.hpp"
#include "cleanseams/routing.hpp"
#include <limits>
namespace cleanseams {
// ============================================================ boundary sides (topology/boundaries.cpp)
// Every boundary side of a cut chart is one original corner whose outgoing edge is cut. Its
// kind is the reason the edge is cut: a physical (mesh) boundary, a forced cut (nonmanifold
// or inconsistently wound edge), or an internal seam. A seam side also names its twin, the
// corner of the same edge on the other side, and the twin's chart: a slit has both sides in
// one chart, a seam between two charts does not.
enum class BoundarySideKind : std::uint8_t {physical,forced,seam};
CLEANSEAMS_CORE_EXPORT std::string_view boundary_side_kind_name(BoundarySideKind);
struct CLEANSEAMS_CORE_EXPORT BoundarySide {
    MeshIndex corner=0,edge=0,loop=0;
    BoundarySideKind kind=BoundarySideKind::physical;
    std::optional<MeshIndex> twin,twin_chart;
    nlohmann::json to_json() const;
};
// Boundary sides of `chart` in the surface's loop order.
CLEANSEAMS_CORE_EXPORT std::vector<BoundarySide> boundary_sides(const CutSurface&,MeshIndex chart);
// Identity of a rim of an underlying patch: its smallest original edge and its edge count.
// Rims of one patch never share edges, so the smallest edge identifies a rim independently of
// where a slit meets it (the ID does not change when the slit moves).
struct PhysicalRimId {
    MeshIndex min_edge=0,edge_count=0;
    auto operator<=>(const PhysicalRimId&) const=default;
    nlohmann::json to_json() const {return {{"min_edge",min_edge},{"edges",edge_count}};}
};
struct CLEANSEAMS_CORE_EXPORT PhysicalRim {
    PhysicalRimId id;
    // physical: every side a physical boundary edge; seam: every side an internal cut towards
    // another chart (a patch rim); mixed: both.
    std::string kind;
    std::vector<MeshIndex> vertices,edges,corners; // loop order; side k runs vertices[k] -> vertices[k+1]
    std::vector<double> arc;                       // cumulative length at vertices[k]
    double length=0;
    std::optional<std::size_t> position(MeshIndex vertex) const;
    // Shortest distance along the rim between two of its vertices (either direction).
    double distance(MeshIndex a,MeshIndex b) const;
    nlohmann::json to_json() const;
};
// The underlying patch of a chart before its slit: the surface with `slit_edges` glued. It is
// the intended tube when its chart is an annulus (Euler characteristic 0, two loops, no pinch,
// every fan supported) whose two loops are simple, disjoint rims without forced or slit sides.
struct CLEANSEAMS_CORE_EXPORT PatchRims {
    bool annulus=false;
    // "" for the intended annulus; otherwise closed_patch (no rim), not_tube (one loop: the
    // slit is a spur of a disk), handle (two loops, Euler characteristic below 0),
    // extra_holes (more than two loops), unsupported_patch, pinched_patch, forced_rim,
    // slit_side_rim, pinched_rim or rims_intersect.
    std::string reason;
    std::shared_ptr<const CutSurface> patch;
    MeshIndex patch_chart=0;
    std::int64_t euler=0;
    Offset loops=0;
    std::vector<PhysicalRim> rims;   // two for an annulus, ordered by PhysicalRimId
    nlohmann::json to_json() const;
};
// `slit_edges` must be internal cut edges of `chart` (both faces in it).
CLEANSEAMS_CORE_EXPORT PatchRims patch_rims(const std::shared_ptr<const CutSurface>&,MeshIndex chart,std::span<const MeshIndex> slit_edges);
// ============================================================ multi-endpoint routing (routing/search.cpp)
// A route end: a vertex, a nonnegative displacement cost of using it (for example the slide
// distance from the current endpoint) and an optional required terminal direction: a
// source's direction of departure, a goal's direction of arrival.
struct RouteEndpoint {
    MeshIndex vertex=0;
    double displacement=0;
    std::optional<Vec3> tangent;
    double tangent_weight=0;
};
// Direction-aware search over directed incident-edge states (previous -> current), like
// find_path, from every source to any goal:
//   start   displacement(s) + tangent_weight(s) * turn(tangent(s), departure) + edge_cost(e1)
//   step    edge_cost(e_out) + w_turn * (len(e_in) + len(e_out)) / 2 * turn(arrival, departure)
//   arrive  displacement(g) + tangent_weight(g) * turn(arrival, tangent(g))
// with turn = (1 - cos)/2. Every term must be finite and nonnegative (checked: invalid_input),
// which makes the lower-bound stop exact: once the cheapest open state costs at least the best
// completed route, no extension or arrival can beat it. Interior vertices are
// `vertex_allowed` and neither a source nor a goal; edges are `edge_allowed`. Ties keep the
// first route found in (cost, previous, current) order; sources and goals are taken in
// ascending vertex order, so their listing order never changes the result.
struct MultiEndpointQuery {
    std::vector<RouteEndpoint> sources,goals;
    std::span<const double> edge_cost;
    std::span<const std::uint8_t> edge_allowed,vertex_allowed;
    double w_turn=0;
};
struct CLEANSEAMS_CORE_EXPORT MultiEndpointRoute {
    std::optional<std::vector<MeshIndex>> verts;   // source ... goal
    std::vector<MeshIndex> edges;
    double cost=0;
    std::optional<MeshIndex> source,goal;
    std::string reason;                            // no_route, route_not_simple, no_endpoints
    std::uint64_t pops=0,pushes=0;
    bool early_stop=false;
    nlohmann::json to_json() const;
};
CLEANSEAMS_CORE_EXPORT MultiEndpointRoute find_multi_endpoint_path(const MeshSnapshot&,const MultiEndpointQuery&,RoutingWorkspace&);
// The cost of a given route under the query's terms (the route must start at a source and end
// at a goal); used to check reported costs.
CLEANSEAMS_CORE_EXPORT double multi_endpoint_route_cost(const MeshSnapshot&,const MultiEndpointQuery&,std::span<const MeshIndex> verts);
// A corridor over one patch. Distances follow the patch's own edges (every face of the edge
// in the patch) from the seeds, never Euclidean proximity, so a nearby limb or sheet is never
// reached. With a radius, the band is every vertex within that surface distance; without
// one, the whole patch.
struct CLEANSEAMS_CORE_EXPORT PatchCorridor {
    std::vector<std::uint8_t> edge_allowed,vertex_allowed;
    std::vector<double> distance;                  // +inf outside the patch or unreached
    double radius=std::numeric_limits<double>::infinity();
    bool whole_patch=false;
    std::uint64_t patch_vertices=0,band_vertices=0,band_edges=0;
    nlohmann::json to_json() const;
};
// `edge_usable`/`vertex_usable` are the admissibility masks (per mesh edge/vertex); the
// corridor keeps the usable ones inside the band.
CLEANSEAMS_CORE_EXPORT PatchCorridor patch_corridor(const MeshSnapshot&,std::span<const std::uint8_t> patch_faces,
    std::span<const std::uint8_t> edge_usable,std::span<const std::uint8_t> vertex_usable,std::span<const MeshIndex> seeds,
    std::optional<double> radius);
// ============================================================ endpoint eligibility (proposals/endpoint_slide.cpp)
// Why a rim vertex may (not) carry the slit's endpoint.
struct CLEANSEAMS_CORE_EXPORT EndpointClass {
    MeshIndex vertex=0;
    std::optional<MeshIndex> sector;   // the patch sector at the vertex (one per rim vertex)
    double slide_distance=0;           // along the rim from the current endpoint
    bool current=false;
    bool anchor=false;                 // manual anchor: never slides, never receives the endpoint
    bool locked_incident=false;        // touches a locked or required edge other than the rim
    bool pinned=false;                 // a pinned corner of the patch at the vertex
    bool neighbor_seam=false;          // another seam (not the rim, not the slit) meets the vertex
    bool junction=false;               // not exactly two rim cuts at the vertex besides the slit
    bool unsupported=false;            // nonmanifold/unsupported vertex fan
    bool pinched=false;                // the vertex has several patch sectors or repeats on the rim
    bool beyond_distance=false;        // farther than the configured maximum slide distance
    bool protected_vertex() const {return anchor || locked_incident || pinned || neighbor_seam || junction || unsupported || pinched;}
    bool admissible() const {return !protected_vertex() && !beyond_distance;}
    std::vector<std::string> reasons() const;
    nlohmann::json to_json() const;
};
struct CLEANSEAMS_CORE_EXPORT RimEndpoints {
    PhysicalRim rim;
    MeshIndex current=0;               // where the slit meets this rim now
    bool fixed=false;                  // the current endpoint is protected: it keeps its place
    std::string fixed_reason;
    std::vector<EndpointClass> vertices;   // every rim vertex, rim order
    std::vector<MeshIndex> admissible;     // sorted: vertices a route may end at (the current endpoint included)
    std::vector<MeshIndex> movable;        // admissible minus the current endpoint
    std::vector<MeshIndex> movable_sectors;// patch sectors of `movable`, same order
    nlohmann::json to_json() const;
};
struct CLEANSEAMS_CORE_EXPORT EndpointEligibility {
    bool eligible=false;
    // "" when eligible; otherwise the first failing rule: multiple_slits, branching_slit,
    // closed_slit, a PatchRims reason, side_provenance, slit_touches_rim, endpoint_off_rim,
    // same_rim, locked_slit, required_slit, anchor_on_slit, pinned_slit_sector,
    // both_endpoints_fixed or no_movable_endpoint.
    std::string reason;
    nlohmann::json detail=nlohmann::json::object();
    MeshIndex chart=0;
    std::string chart_key;
    std::vector<MeshIndex> slit_edges;     // sorted
    std::vector<MeshIndex> slit_path;      // vertices, from rims[0].current to rims[1].current
    Offset slit_chart_loops=0;             // boundary loops of the slit chart now (1 for a slit disk)
    std::vector<BoundarySide> sides;       // boundary-side provenance of the slit chart
    std::uint64_t rim_sides=0,slit_sides=0;// sides on the rims and on the two slit sides
    PatchRims patch;
    std::array<RimEndpoints,2> rims;
    nlohmann::json to_json() const;
};
// One record per chart holding internal slit edges (cuts with both faces in the chart), in
// chart order. Only a single unlocked, nonbranching slit joining two distinct rims of an
// annular patch is eligible; everything else gets a reason and is never forced into the tube case.
// GB m7: the patch analysis rebuilds a surface per through-slit chart, so `control` is checked
// before each: a cancellation throws, and after the deadline the remaining charts are recorded
// with reason "deadline" (never analysed, never eligible).
CLEANSEAMS_CORE_EXPORT std::vector<EndpointEligibility> endpoint_eligibility(const std::shared_ptr<const CutSurface>&,const HardConstraints&,
    const EndpointSlideOptions&,const Control& control={});
// The admissible source (rims[0]) and goal (rims[1]) sets of an eligible slit as route ends:
// displacement w_slide * slide distance, terminal direction along the patch's interior edges
// at the vertex (departure into the patch, arrival towards the rim) weighted by w_tangent
// times the mean length of the rim sides at the vertex.
struct SlideEndpoints {std::vector<RouteEndpoint> sources,goals;};
CLEANSEAMS_CORE_EXPORT SlideEndpoints slide_endpoints(const MeshSnapshot&,const EndpointEligibility&,double w_slide,double w_tangent);
// The corridor of an eligible slit: interior edges of its patch that are uncut (or the slit's
// own) and not forbidden; interior vertices that are supported, unpinned and not anchors.
// `radius` none is the whole patch.
CLEANSEAMS_CORE_EXPORT PatchCorridor slide_corridor(const CutSurface&,const EndpointEligibility&,const HardConstraints&,std::optional<double> radius);
enum class DistanceMetric;
// T39: the same corridor over a distance metric (physical surface distance or hop count).
CLEANSEAMS_CORE_EXPORT PatchCorridor slide_corridor(const CutSurface&,const EndpointEligibility&,const HardConstraints&,std::optional<double> radius,
    DistanceMetric);
// The corridor schedule of a slide search (Sections 13.2 and 20): the whole patch when it has at
// most `vertex_budget` vertices; otherwise surface-distance bands around the slit of radius
// r0 * 2^k, with r0 = max(1, corridor_rings) * the median edge length of the slit chart (captured
// once for the operation, never recomputed as cuts change), widened deterministically while the
// next band still holds at most `vertex_budget` patch vertices (the first band is always kept).
// The search uses the last, widest entry: every narrower band is a subset of it, so its optimum
// is never worse, and the budget bounds the work.
struct CLEANSEAMS_CORE_EXPORT SlideCorridorSchedule {
    std::vector<std::optional<double>> radii;
    std::vector<std::uint64_t> band_vertices;  // patch vertices within each radius
    double reference_length=0;         // median edge length of the slit chart (T39: the frozen length unit)
    std::uint64_t patch_vertices=0;
    // T39 (PhysicalScale overload): metric (physical or hops), the resolved first radius with its
    // source (explicit, reference_length, reference or hops) and the chart's current median
    // edge length, which the schedule never follows.
    std::string metric;
    nlohmann::json radius=nlohmann::json::object();
    double local_length=0;
    std::optional<double> search_radius() const {return radii.empty()?std::nullopt:radii.back();}
    nlohmann::json to_json() const;
};
CLEANSEAMS_CORE_EXPORT SlideCorridorSchedule slide_corridor_schedule(const CutSurface&,const EndpointEligibility&,const EndpointSlideOptions&,
    std::int64_t corridor_rings);
// T39: the schedule of an enhancement pass under its frozen physical scale (Section 20): the
// first radius is PhysicalScale::corridor_radius of the slit chart's faces (explicit, or the
// captured corridor rings times the frozen reference length; hop rings in hops mode).
class PhysicalScale;
CLEANSEAMS_CORE_EXPORT SlideCorridorSchedule slide_corridor_schedule(const CutSurface&,const EndpointEligibility&,const EndpointSlideOptions&,
    const PhysicalScale&);
// ============================================================ whole-slit relocation (T31)
// Section 13.2: for every eligible slit of `state`, up to four distinct routes between its two
// admissible endpoint sets, one per cost variant, searched over the corridor schedule:
//   preference  the enhanced routing heuristic (paint, explicit viewer, raw exposure; T29)
//   free_slide  the same without the slide displacement cost
//   straight    the same with four times the turn weight
//   shortest    physical edge length only
// Routes are deduplicated by their complete slit edge set (the current slit included, so the
// old path is never proposed again). Each distinct route becomes one proposal: the old slit's
// edges removed and the new path's added in one delta, the endpoint moves declared against the
// admissible sets, predicted route cost for ordering only, and provenance recording moved
// endpoints, slide distances, raw exposed length and seam length before and after. A cheap
// topology screen drops routes that would not leave one disk in place of the slit chart with
// every other chart unchanged; the transaction checks everything again. `report` receives one
// record per slit (eligibility, schedule, every variant's route and fate).
CLEANSEAMS_CORE_EXPORT std::vector<Proposal> endpoint_slide_proposals(const RunState&,ProposalContext&,nlohmann::json* report=nullptr);
// The Section 21 "sliding proposals on eligible existing tube slits" generator (name
// endpoint_slide, version 1); options from ProposalContext::options.enhanced.endpoint_slide.
CLEANSEAMS_CORE_EXPORT ProposalGenerator endpoint_slide_generator();
}
