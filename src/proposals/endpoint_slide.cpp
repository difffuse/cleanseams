// SPDX-License-Identifier: GPL-3.0-or-later
// Phase B sliding seam endpoints on tube-like patches (Section 13). T30: endpoint
// eligibility on the underlying patch, explicit admissible endpoint sets, their route ends and
// the patch corridor. T31: whole-slit relocation proposals through the common transaction.
// Used by policy enhanced_v1 only.
#include "cleanseams/endpoints.hpp"
#include "cleanseams/physical_scale.hpp"
#include "cleanseams/cache.hpp"
#include "cleanseams/errors.hpp"
#include "cleanseams/control.hpp"
#include "cleanseams/detail/mesh_math.hpp"
#include "cleanseams/detail/numpy_compat.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <map>
#include <set>
namespace cleanseams {
namespace {
using Json=nlohmann::json;
bool always_cut(const MeshTables& t,std::size_t e) {
    return t.edge_boundary[e] || !t.edge_manifold[e] || !t.edge_consistent[e] || t.edge_corners[e].size()!=2;
}
void require(bool yes,const char* field,const char* message) {if(!yes)throw ContractError(ErrorCode::invalid_input,field,message);}
bool contains(const std::vector<MeshIndex>& sorted,MeshIndex x) {return std::binary_search(sorted.begin(),sorted.end(),x);}
// Pinned corners of one chart, by vertex.
std::set<MeshIndex> pinned_vertices(const CutSurface& cs,const HardConstraints& hc,MeshIndex chart) {
    std::set<MeshIndex> out;const auto& ct=cs.tables();
    for(const auto& [c,uv]:hc.pins)if(c<ct.chart_of_corner.size() && ct.chart_of_corner[c]==chart)
        out.insert(static_cast<MeshIndex>(cs.mesh().input().corner_vertices[c]));
    return out;
}
// The slit's vertex path from `start` along the slit edges (a simple path by construction).
std::vector<MeshIndex> walk_slit(const MeshSnapshot& m,const std::map<MeshIndex,std::vector<MeshIndex>>& incident,MeshIndex start) {
    std::vector<MeshIndex> path{start};std::set<MeshIndex> used;auto at=start;
    while(true) {
        std::optional<MeshIndex> next;
        for(auto e:incident.at(at))if(!used.contains(e)) {next=e;break;}
        if(!next)break;
        used.insert(*next);at=m.edge_other(*next,at);path.push_back(at);
    }
    return path;
}
RimEndpoints classify_rim(const CutSurface& cs,const EndpointEligibility& rec,const PhysicalRim& rim,MeshIndex current,const HardConstraints& hc,
                          const EndpointSlideOptions& o,const std::set<MeshIndex>& pinned) {
    const auto& m=cs.mesh();const auto& t=m.tables();const auto& ct=cs.tables();
    const auto& pt=rec.patch.patch->tables();
    std::vector<MeshIndex> rim_edges=rim.edges;std::sort(rim_edges.begin(),rim_edges.end());
    RimEndpoints out;out.rim=rim;out.current=current;
    for(auto v:rim.vertices) {
        EndpointClass k;k.vertex=v;k.current=v==current;k.slide_distance=rim.distance(current,v);
        k.anchor=std::binary_search(hc.anchors.begin(),hc.anchors.end(),v);
        std::uint64_t cuts=0;
        for(auto e:t.vert_edges[v]) {
            if(contains(rec.slit_edges,e))continue;
            if(ct.cut[e])++cuts;
            if(contains(rim_edges,e))continue;
            if((e<hc.locked.size() && hc.locked[e]) || (e<hc.required.size() && hc.required[e]))k.locked_incident=true;
            if(ct.cut[e])k.neighbor_seam=true;
        }
        k.junction=cuts!=2;
        k.unsupported=!t.fans[v];
        k.pinned=pinned.contains(v);
        std::set<MeshIndex> sectors;
        for(auto c:t.vert_corners[v])if(pt.chart_of_corner[c]==rec.patch.patch_chart)sectors.insert(pt.sector_of_corner[c]);
        if(sectors.size()==1)k.sector=*sectors.begin();else k.pinched=true;
        k.beyond_distance=o.max_slide_distance>0 && k.slide_distance>o.max_slide_distance*(1.0+1e-12);
        out.vertices.push_back(k);
    }
    const auto& here=*std::find_if(out.vertices.begin(),out.vertices.end(),[&](const EndpointClass& k) {return k.current;});
    out.fixed=here.protected_vertex();
    if(out.fixed) {
        std::string why;for(const auto& r:here.reasons()) {if(!why.empty())why+=",";why+=r;}
        out.fixed_reason=why;out.admissible={current};
    } else for(const auto& k:out.vertices)if(k.admissible())out.admissible.push_back(k.vertex);
    std::sort(out.admissible.begin(),out.admissible.end());
    for(auto v:out.admissible) {
        if(v==current)continue;
        out.movable.push_back(v);
        out.movable_sectors.push_back(*std::find_if(out.vertices.begin(),out.vertices.end(),[&](const EndpointClass& k) {return k.vertex==v;})->sector);
    }
    return out;
}
// Direction from a rim vertex into its patch: the mean of the patch's interior edge
// directions at the vertex (tangent-plane projected).
std::optional<Vec3> inward(const MeshSnapshot& m,const EndpointEligibility& rec,const RimEndpoints& r,MeshIndex v) {
    const auto& t=m.tables();const auto& pt=rec.patch.patch->tables();
    Vec3 sum{0,0,0};std::size_t n=0;
    for(auto e:t.vert_edges[v]) {
        if(always_cut(t,e) || std::find(r.rim.edges.begin(),r.rim.edges.end(),e)!=r.rim.edges.end())continue;
        bool inside=true;for(auto f:t.edge_faces[e])inside=inside && pt.chart_of_face[f]==rec.patch.patch_chart;
        if(!inside)continue;
        sum=detail::add(sum,m.tangent_dir(v,e));++n;
    }
    const double length=detail::norm(sum);
    if(!n || !(length>1e-12))return std::nullopt;
    return detail::divide(sum,length);
}
double rim_side_length(const MeshSnapshot& m,const PhysicalRim& rim,MeshIndex v) {
    const auto& t=m.tables();double sum=0;std::size_t n=0;
    for(auto e:rim.edges)if(t.edge_vertices[e][0]==v || t.edge_vertices[e][1]==v) {sum+=m.geometry().edge_length[e];++n;}
    return n?sum/static_cast<double>(n):0.0;
}
}
std::vector<std::string> EndpointClass::reasons() const {
    std::vector<std::string> out;
    if(anchor)out.push_back("anchor");
    if(locked_incident)out.push_back("locked_incident");
    if(pinned)out.push_back("pinned");
    if(neighbor_seam)out.push_back("neighbor_seam");
    if(junction)out.push_back("junction");
    if(unsupported)out.push_back("unsupported");
    if(pinched)out.push_back("pinched");
    if(beyond_distance)out.push_back("beyond_distance");
    return out;
}
nlohmann::json EndpointClass::to_json() const {
    return {{"vertex",vertex},{"sector",sector?Json(*sector):Json(nullptr)},{"slide_distance",slide_distance},{"current",current},
            {"admissible",admissible()},{"reasons",reasons()}};
}
nlohmann::json RimEndpoints::to_json() const {
    Json excluded=Json::object();
    for(const auto& k:vertices)if(!k.admissible())excluded[std::to_string(k.vertex)]=k.reasons();
    return {{"rim",rim.to_json()},{"current",current},{"fixed",fixed},{"fixed_reason",fixed_reason},{"admissible",admissible},
            {"movable",movable},{"movable_sectors",movable_sectors},{"excluded",excluded}};
}
nlohmann::json EndpointEligibility::to_json() const {
    Json sides_by_kind=Json::object();
    for(const auto& s:sides)sides_by_kind[std::string(boundary_side_kind_name(s.kind))]=sides_by_kind.value(std::string(boundary_side_kind_name(s.kind)),0)+1;
    Json out={{"eligible",eligible},{"reason",reason},{"detail",detail},{"chart",chart},{"chart_key",chart_key},{"slit_edges",slit_edges},
              {"slit_path",slit_path},{"slit_chart_loops",slit_chart_loops},{"sides",sides_by_kind},{"rim_sides",rim_sides},{"slit_sides",slit_sides},
              {"patch",patch.to_json()}};
    if(eligible)out["rims"]={rims[0].to_json(),rims[1].to_json()};
    return out;
}
std::vector<EndpointEligibility> endpoint_eligibility(const std::shared_ptr<const CutSurface>& surface,const HardConstraints& hc,const EndpointSlideOptions& o,
                                                     const Control& control) {
    require(bool(surface),"endpoints.surface","surface is required");
    require(std::isfinite(o.max_slide_distance) && o.max_slide_distance>=0,"endpoints.max_slide_distance","must be finite and nonnegative");
    const auto& cs=*surface;const auto& ct=cs.tables();const auto& m=cs.mesh();const auto& t=m.tables();const auto ne=t.edge_vertices.size();
    // Internal slit edges per chart: effective cuts that are not topology facts, with both faces in one chart.
    std::map<MeshIndex,std::vector<MeshIndex>> slits;
    for(std::size_t e=0;e<ne;++e) {
        if(!ct.cut[e] || always_cut(t,e))continue;
        const auto a=ct.chart_of_face[t.edge_faces[e][0]],b=ct.chart_of_face[t.edge_faces[e][1]];
        if(a==b)slits[a].push_back(static_cast<MeshIndex>(e));
    }
    std::vector<EndpointEligibility> out;
    for(auto& [chart,edges]:slits) {
        EndpointEligibility rec;rec.chart=chart;rec.chart_key=cs.chart_key(chart).digest();rec.slit_edges=edges;
        rec.slit_chart_loops=ct.chart_n_loops[chart];rec.sides=boundary_sides(cs,chart);
        auto refuse=[&](const char* reason,Json detail=Json::object()) {rec.reason=reason;rec.detail=std::move(detail);out.push_back(std::move(rec));};
        // Structure: one nonbranching open slit.
        std::map<MeshIndex,std::vector<MeshIndex>> incident;
        for(auto e:edges)for(auto v:t.edge_vertices[e])incident[v].push_back(e);
        std::vector<MeshIndex> branching,ends;
        for(const auto& [v,list]:incident) {if(list.size()>2)branching.push_back(v);if(list.size()==1)ends.push_back(v);}
        std::size_t components=0;
        {
            std::set<MeshIndex> seen;
            for(const auto& [v,list]:incident) {
                if(seen.contains(v))continue;++components;std::vector<MeshIndex> stack{v};seen.insert(v);
                while(!stack.empty()) {
                    const auto x=stack.back();stack.pop_back();
                    for(auto e:incident.at(x)) {const auto y=m.edge_other(e,x);if(seen.insert(y).second)stack.push_back(y);}
                }
            }
        }
        if(components>1) {refuse("multiple_slits",{{"components",components}});continue;}
        if(!branching.empty()) {refuse("branching_slit",{{"vertices",branching}});continue;}
        if(ends.size()!=2) {refuse("closed_slit");continue;}
        // A spur of a disk chart (one free end touching no other cut): gluing it leaves the same
        // disk, so the underlying-patch analysis would answer not_tube. Answered here without
        // building the patch surface (dead-end relief slits are common); every other shape takes
        // the full analysis below.
        std::vector<MeshIndex> interior_ends;
        for(auto v:ends) {
            bool other=false;for(auto e:t.vert_edges[v])other=other || (ct.cut[e] && !contains(edges,e));
            if(!other)interior_ends.push_back(v);
        }
        if(ct.chart_is_disk[chart] && interior_ends.size()==1) {refuse("not_tube",{{"spur_end",interior_ends[0]},{"loops",1}});continue;}
        rec.slit_path=walk_slit(m,incident,ends[0]);
        // GB m7: the patch analysis below rebuilds a whole surface; the pass's control first.
        if(control.check()) {refuse("deadline");continue;}
        // The underlying patch before the slit: the intended annulus with two distinct rims.
        rec.patch=patch_rims(surface,chart,edges);
        if(!rec.patch.annulus) {refuse(rec.patch.reason.c_str(),{{"euler",rec.patch.euler},{"loops",rec.patch.loops}});continue;}
        // Boundary-side provenance: the slit disk's boundary is exactly both rims plus both sides
        // of every slit edge.
        std::uint64_t unaccounted=0;
        for(const auto& s:rec.sides) {
            if(s.kind==BoundarySideKind::seam && s.twin_chart && *s.twin_chart==chart && contains(edges,s.edge))++rec.slit_sides;
            else if(s.kind==BoundarySideKind::physical || (s.kind==BoundarySideKind::seam && s.twin_chart && *s.twin_chart!=chart))++rec.rim_sides;
            else ++unaccounted;
        }
        const auto& r0=rec.patch.rims[0];const auto& r1=rec.patch.rims[1];
        if(unaccounted || rec.slit_sides!=2*edges.size() || rec.rim_sides!=r0.edges.size()+r1.edges.size()) {
            refuse("side_provenance",{{"rim_sides",rec.rim_sides},{"slit_sides",rec.slit_sides},{"unaccounted",unaccounted}});continue;
        }
        auto rim_of=[&](MeshIndex v) -> int {return r0.position(v)?0:(r1.position(v)?1:-1);};
        std::vector<MeshIndex> touching;
        for(std::size_t k=1;k+1<rec.slit_path.size();++k)if(rim_of(rec.slit_path[k])>=0)touching.push_back(rec.slit_path[k]);
        if(!touching.empty()) {refuse("slit_touches_rim",{{"vertices",touching}});continue;}
        const int ra=rim_of(rec.slit_path.front()),rb=rim_of(rec.slit_path.back());
        if(ra<0 || rb<0) {refuse("endpoint_off_rim",{{"ends",{rec.slit_path.front(),rec.slit_path.back()}}});continue;}
        if(ra==rb) {refuse("same_rim",{{"rim",rec.patch.rims[static_cast<std::size_t>(ra)].id.to_json()}});continue;}
        if(ra==1)std::reverse(rec.slit_path.begin(),rec.slit_path.end());
        // Hard constraints on the slit itself.
        std::vector<MeshIndex> locked,required,anchored;
        for(auto e:edges) {
            if(e<hc.locked.size() && hc.locked[e])locked.push_back(e);
            if(e<hc.required.size() && hc.required[e])required.push_back(e);
        }
        for(std::size_t k=1;k+1<rec.slit_path.size();++k)
            if(std::binary_search(hc.anchors.begin(),hc.anchors.end(),rec.slit_path[k]))anchored.push_back(rec.slit_path[k]);
        if(!locked.empty()) {refuse("locked_slit",{{"edges",locked}});continue;}
        if(!required.empty()) {refuse("required_slit",{{"edges",required}});continue;}
        if(!anchored.empty()) {refuse("anchor_on_slit",{{"anchors",anchored}});continue;}
        // Moving the slit merges the sectors along it: a pin there would change its sector constraint.
        const auto pinned=pinned_vertices(cs,hc,chart);
        std::vector<MeshIndex> pinned_on_slit;for(auto v:rec.slit_path)if(pinned.contains(v))pinned_on_slit.push_back(v);
        if(!pinned_on_slit.empty()) {refuse("pinned_slit_sector",{{"vertices",pinned_on_slit}});continue;}
        rec.rims[0]=classify_rim(cs,rec,r0,rec.slit_path.front(),hc,o,pinned);
        rec.rims[1]=classify_rim(cs,rec,r1,rec.slit_path.back(),hc,o,pinned);
        if(rec.rims[0].fixed && rec.rims[1].fixed) {
            refuse("both_endpoints_fixed",{{"reasons",{rec.rims[0].fixed_reason,rec.rims[1].fixed_reason}}});continue;
        }
        if(rec.rims[0].movable.empty() && rec.rims[1].movable.empty()) {refuse("no_movable_endpoint");continue;}
        rec.eligible=true;out.push_back(std::move(rec));
    }
    return out;
}
SlideEndpoints slide_endpoints(const MeshSnapshot& m,const EndpointEligibility& rec,double w_slide,double w_tangent) {
    require(rec.eligible,"endpoints.slit","slide endpoints need an eligible slit");
    require(std::isfinite(w_slide) && w_slide>=0 && std::isfinite(w_tangent) && w_tangent>=0,"endpoints.weights","weights must be finite and nonnegative");
    SlideEndpoints out;
    for(int side=0;side<2;++side) {
        const auto& r=rec.rims[static_cast<std::size_t>(side)];
        for(auto v:r.admissible) {
            RouteEndpoint end;end.vertex=v;
            end.displacement=w_slide*r.rim.distance(r.current,v);
            if(auto d=inward(m,rec,r,v)) {
                // Sources leave into the patch; goals arrive heading out to their rim.
                end.tangent=side==0?*d:detail::mul(*d,-1.0);end.tangent_weight=w_tangent*rim_side_length(m,r.rim,v);
            }
            (side==0?out.sources:out.goals).push_back(end);
        }
    }
    return out;
}
namespace {
// The slit chart's faces and the admissibility masks of its slide corridor.
struct SlideMasks {std::vector<std::uint8_t> faces,edge_usable,vertex_usable;};
SlideMasks slide_masks(const CutSurface& cs,const EndpointEligibility& rec,const HardConstraints& hc) {
    require(rec.eligible,"endpoints.slit","a slide corridor needs an eligible slit");
    const auto& m=cs.mesh();const auto& t=m.tables();const auto& ct=cs.tables();
    const auto ne=t.edge_vertices.size(),nv=m.input().positions.size(),nf=ct.chart_of_face.size();
    SlideMasks out;auto& faces=out.faces;auto& edge_usable=out.edge_usable;auto& vertex_usable=out.vertex_usable;
    faces.assign(nf,0);for(auto f:ct.chart_faces[rec.chart])faces[f]=1;
    edge_usable.assign(ne,0);vertex_usable.assign(nv,0);
    for(std::size_t e=0;e<ne;++e) {
        if(always_cut(t,e) || (e<hc.forbidden.size() && hc.forbidden[e]))continue;
        edge_usable[e]=!ct.cut[e] || contains(rec.slit_edges,static_cast<MeshIndex>(e));
    }
    const auto pinned=pinned_vertices(cs,hc,rec.chart);
    std::vector<std::uint8_t> on_rim(nv,0);for(const auto& r:rec.rims)for(auto v:r.rim.vertices)on_rim[v]=1;
    for(auto f:ct.chart_faces[rec.chart])for(auto c:m.face_corners(f)) {
        const auto v=static_cast<MeshIndex>(m.input().corner_vertices[c]);
        vertex_usable[v]=!on_rim[v] && t.fans[v] && !pinned.contains(v) && !std::binary_search(hc.anchors.begin(),hc.anchors.end(),v);
    }
    return out;
}
}
PatchCorridor slide_corridor(const CutSurface& cs,const EndpointEligibility& rec,const HardConstraints& hc,std::optional<double> radius) {
    const auto masks=slide_masks(cs,rec,hc);
    return patch_corridor(cs.mesh(),masks.faces,masks.edge_usable,masks.vertex_usable,rec.slit_path,radius);
}
PatchCorridor slide_corridor(const CutSurface& cs,const EndpointEligibility& rec,const HardConstraints& hc,std::optional<double> radius,DistanceMetric metric) {
    const auto masks=slide_masks(cs,rec,hc);
    const auto distance=patch_distances(cs.mesh(),masks.faces,rec.slit_path,metric);
    return corridor_band(cs.mesh(),masks.faces,masks.edge_usable,masks.vertex_usable,distance,radius);
}
nlohmann::json SlideCorridorSchedule::to_json() const {
    Json r=Json::array();for(const auto& x:radii)r.push_back(x?Json(*x):Json(nullptr));
    Json out={{"radii",r},{"band_vertices",band_vertices},{"reference_length",reference_length},{"patch_vertices",patch_vertices}};
    if(!metric.empty())out.update({{"metric",metric},{"radius",radius},{"local_length",local_length}});
    return out;
}
SlideCorridorSchedule slide_corridor_schedule(const CutSurface& cs,const EndpointEligibility& rec,const EndpointSlideOptions& o,std::int64_t corridor_rings) {
    require(rec.eligible,"endpoints.slit","a slide corridor needs an eligible slit");
    require(o.vertex_budget>=1,"endpoints.vertex_budget","must be at least 1");
    const auto& m=cs.mesh();const auto& t=m.tables();const auto& ct=cs.tables();
    // The reference length of the operation: the median edge length of the slit chart.
    std::set<MeshIndex> chart_edges;
    for(auto f:ct.chart_faces[rec.chart])for(auto c:m.face_corners(f))chart_edges.insert(t.corner_edge[c]);
    std::vector<double> lengths;for(auto e:chart_edges)lengths.push_back(m.geometry().edge_length[e]);
    SlideCorridorSchedule out;out.reference_length=detail::numpy_median(lengths);
    const auto whole=patch_corridor(m,[&] {std::vector<std::uint8_t> f(ct.chart_of_face.size(),0);for(auto x:ct.chart_faces[rec.chart])f[x]=1;return f;}(),
                                    std::vector<std::uint8_t>(t.edge_vertices.size(),0),std::vector<std::uint8_t>(m.input().positions.size(),0),rec.slit_path,std::nullopt);
    out.patch_vertices=whole.patch_vertices;
    const auto budget=static_cast<std::uint64_t>(o.vertex_budget);
    if(out.patch_vertices<=budget) {out.radii={std::nullopt};out.band_vertices={out.patch_vertices};return out;}
    // Patch vertices within a radius (surface distances along the patch's own edges).
    std::vector<double> distances;for(auto d:whole.distance)if(std::isfinite(d))distances.push_back(d);
    std::sort(distances.begin(),distances.end());
    auto within=[&](double r) {return static_cast<std::uint64_t>(std::upper_bound(distances.begin(),distances.end(),r)-distances.begin());};
    double radius=static_cast<double>(std::max<std::int64_t>(1,corridor_rings))*out.reference_length;
    if(!(radius>0))radius=distances.empty()?1.0:std::max(distances.back(),1e-300);
    // The first band is always kept; each doubling is kept while it stays within the budget. The
    // whole patch exceeds the budget here, so the widening stops before it.
    out.radii.push_back(radius);out.band_vertices.push_back(within(radius));
    while(out.radii.size()<64 && within(2.0*radius)<=budget) {radius*=2.0;out.radii.push_back(radius);out.band_vertices.push_back(within(radius));}
    return out;
}
SlideCorridorSchedule slide_corridor_schedule(const CutSurface& cs,const EndpointEligibility& rec,const EndpointSlideOptions& o,const PhysicalScale& scale) {
    require(rec.eligible,"endpoints.slit","a slide corridor needs an eligible slit");
    require(o.vertex_budget>=1,"endpoints.vertex_budget","must be at least 1");
    const auto& m=cs.mesh();const auto& t=m.tables();const auto& ct=cs.tables();
    const auto faces=ct.chart_faces[rec.chart];
    // The radius of the operation (Section 20): frozen at its start on the slit's reference
    // chart (the run origin's chart holding most of the slit chart's area), or the explicit
    // object-space radius; never the median edge length of the chart as trial cuts left it.
    const auto radius=scale.corridor_radius(faces);
    const auto metric=scale.physical()?DistanceMetric::physical:DistanceMetric::hops;
    std::vector<std::uint8_t> mask(ct.chart_of_face.size(),0);for(auto f:faces)mask[f]=1;
    const auto distance=patch_distances(m,mask,rec.slit_path,metric);
    std::uint64_t patch_vertices=0;{
        std::vector<std::uint8_t> seen(m.input().positions.size(),0);
        for(auto f:faces)for(auto c:m.face_corners(f)) {const auto v=static_cast<std::size_t>(m.input().corner_vertices[c]);patch_vertices+=!seen[v];seen[v]=1;}
    }
    const auto schedule=corridor_schedule(distance,patch_vertices,radius.value,static_cast<std::uint64_t>(o.vertex_budget),metric);
    SlideCorridorSchedule out;out.radii=schedule.radii;out.band_vertices=schedule.band_vertices;out.patch_vertices=schedule.patch_vertices;
    out.metric=schedule.metric;out.radius=radius.to_json();out.reference_length=scale.physical()?scale.unit(faces).value:0.0;
    // The chart's median edge length now, reported next to the frozen reference (they differ
    // once cuts changed the chart; the schedule never follows it).
    std::set<MeshIndex> chart_edges;for(auto f:faces)for(auto c:m.face_corners(f))chart_edges.insert(t.corner_edge[c]);
    std::vector<double> lengths;for(auto e:chart_edges)lengths.push_back(m.geometry().edge_length[e]);
    out.local_length=lengths.empty()?0.0:detail::numpy_median(lengths);
    return out;
}
// ============================================================ T31: whole-slit relocation
namespace {
struct Variant {const char* name;bool length_only;double turn_factor;bool slide_cost;};
constexpr std::array<Variant,4> slide_variants{{
    {"preference",false,1.0,true},
    {"free_slide",false,1.0,false},
    {"straight",false,4.0,true},
    {"shortest",true,1.0,true},
}};
constexpr std::size_t max_slide_routes=4;
Json path_record(const MeshSnapshot& m,std::span<const double> exposure,std::span<const MeshIndex> edges) {
    double length=0,exposed=0;
    for(auto e:edges) {length+=m.geometry().edge_length[e];if(!exposure.empty())exposed+=m.geometry().edge_length[e]*exposure[e];}
    return {{"edges",std::vector<MeshIndex>(edges.begin(),edges.end())},{"seam_length",length},{"exposed_length",exposure.empty()?Json(nullptr):Json(exposed)}};
}
// The cheap topology screen of a trial: one disk where the slit chart was, every other chart
// unchanged (same chart keys), the same chart count.
std::string topology_screen(const CutSurface& parent,const CutSurface& trial,const EndpointEligibility& rec) {
    const auto& pt=parent.tables();const auto& tt=trial.tables();
    if(tt.n_charts!=pt.n_charts)return "chart_count";
    const auto ch=tt.chart_of_face[pt.chart_faces[rec.chart][0]];
    if(!tt.chart_is_disk[ch] || tt.chart_faces[ch].size()!=pt.chart_faces[rec.chart].size())return "not_disk";
    std::set<std::string> keys;for(MeshIndex k=0;k<tt.n_charts;++k)keys.insert(trial.chart_key(k).digest());
    for(MeshIndex k=0;k<pt.n_charts;++k)if(k!=rec.chart && !keys.contains(parent.chart_key(k).digest()))return "outside_changed";
    return "";
}
}
std::vector<Proposal> endpoint_slide_proposals(const RunState& state,ProposalContext& ctx,nlohmann::json* report) {
    require(bool(state.surface) && state.surface->mesh_owner()==ctx.mesh,"endpoints.state","the state belongs to another mesh");
    const auto& m=*ctx.mesh;const auto& o=ctx.options.enhanced.endpoint_slide;
    auto& tally=ctx.counters.by_generator["endpoint_slide"];
    // Routing cost: the enhanced heuristic under the pass's paint, viewer and raw exposure.
    std::shared_ptr<const ResolvedPreferences> prefs=ctx.preferences;
    if(!prefs) {
        auto p=resolve_preferences(m,{},ctx.options.enhanced.preferences,std::nullopt,{},ctx.options.layout.feature_angle,ctx.options.layout.hide_dir);
        p.raw_exposure=ctx.raw_exposure;prefs=std::make_shared<const ResolvedPreferences>(std::move(p));
    }
    const auto heuristic=enhanced_edge_cost(m,*prefs,ctx.options.enhanced.preferences);
    const auto& lengths=m.geometry().edge_length;
    std::vector<Proposal> out;Json records=Json::array();std::uint64_t next_id=1;
    RoutingWorkspace ws;
    const Control control{ctx.cancel,ctx.deadline};
    for(const auto& rec:endpoint_eligibility(state.surface,ctx.constraints,o,control)) {
        throw_if_cancelled(ctx.cancel);
        ++tally["slits_examined"];
        Json item={{"chart",rec.chart},{"chart_key",rec.chart_key},{"eligible",rec.eligible},{"reason",rec.reason},{"detail",rec.detail}};
        if(!rec.eligible) {++tally["ineligible:"+rec.reason];records.push_back(std::move(item));continue;}
        if(ctx.deadline.expired()) {++tally["deadline"];item["stopped"]="deadline";records.push_back(std::move(item));continue;}
        ++tally["eligible"];
        // T39: the corridor radius is the pass's frozen physical scale (or the hop rings in
        // the documented hops mode), never re-derived from the chart as cuts changed it.
        const auto& scale=context_scale(ctx,state);
        const auto schedule=slide_corridor_schedule(*state.surface,rec,o,scale);
        const auto band=slide_corridor(*state.surface,rec,ctx.constraints,schedule.search_radius(),
                                       scale.physical()?DistanceMetric::physical:DistanceMetric::hops);
        std::vector<MeshIndex> allowed;for(const auto& r:rec.rims)allowed.insert(allowed.end(),r.admissible.begin(),r.admissible.end());
        std::sort(allowed.begin(),allowed.end());
        // The old path is a known cut set: a route equal to it is never proposed.
        std::set<std::vector<MeshIndex>> seen{rec.slit_edges};
        Json routes=Json::array();std::size_t proposed=0;
        for(const auto& v:slide_variants) {
            if(proposed>=max_slide_routes)break;
            const auto ends=slide_endpoints(m,rec,v.slide_cost?o.w_slide:0.0,o.w_tangent);
            // The widest band of the schedule (the whole patch under the vertex budget): every
            // narrower band is a subset of it, so its optimum is never worse.
            MultiEndpointQuery q;q.sources=ends.sources;q.goals=ends.goals;q.edge_cost=v.length_only?std::span<const double>(lengths):std::span<const double>(heuristic);
            q.edge_allowed=band.edge_allowed;q.vertex_allowed=band.vertex_allowed;q.w_turn=o.w_turn*v.turn_factor;
            const auto route=find_multi_endpoint_path(m,q,ws);const auto corridor=band.to_json();
            Json r={{"variant",v.name},{"route",route.to_json()},{"corridor",corridor}};
            if(!route.verts) {++tally["no_route"];routes.push_back(std::move(r));continue;}
            auto edges=route.edges;std::sort(edges.begin(),edges.end());
            if(!seen.insert(edges).second) {
                const char* why=edges==rec.slit_edges?"unchanged_route":"duplicate_route";++tally[why];r["fate"]=why;routes.push_back(std::move(r));continue;
            }
            // The whole slit is replaced in one delta: old edges removed, new path added.
            EdgeMask cuts=state.cuts;for(auto e:rec.slit_edges)cuts[e]=0;for(auto e:edges)cuts[e]=1;
            const auto trial=cut_surface(ctx.mesh,cuts);
            if(const auto screen=topology_screen(*state.surface,*trial,rec);!screen.empty()) {
                ++tally["screen:"+screen];r["fate"]="screen:"+screen;routes.push_back(std::move(r));continue;
            }
            Proposal p;p.identity.id=next_id++;p.identity.parent_revision=state.revision;p.identity.generator="endpoint_slide";
            p.identity.generator_version="1";p.identity.affected_chart_keys={rec.chart_key};
            for(auto e:edges)if(!state.cuts[e])p.identity.added_edges.push_back(e);
            for(auto e:rec.slit_edges)if(!std::binary_search(edges.begin(),edges.end(),e))p.identity.removed_edges.push_back(e);
            const std::array<MeshIndex,2> to{route.verts->front(),route.verts->back()};
            Json moved=Json::array();
            for(std::size_t k=0;k<2;++k) {
                const auto& r2=rec.rims[k];
                if(to[k]!=r2.current)p.endpoints.push_back({r2.current,to[k]});
                moved.push_back({{"rim",r2.rim.id.to_json()},{"rim_kind",r2.rim.kind},{"from",r2.current},{"to",to[k]},
                                 {"slide_distance",r2.rim.distance(r2.current,to[k])}});
            }
            p.allowed_endpoints=allowed;
            p.predicted={{"route_cost",route.cost},{"variant",v.name}};
            p.provenance={{"operation","endpoint_slide"},{"variant",v.name},{"chart_key",rec.chart_key},{"endpoints",moved},
                          {"before",path_record(m,ctx.raw_exposure,rec.slit_edges)},{"after",path_record(m,ctx.raw_exposure,edges)},
                          {"route_cost",route.cost},{"corridor",corridor}};
            r["fate"]="proposed";r["id"]=p.identity.id;routes.push_back(std::move(r));
            ++tally["routes"];++proposed;out.push_back(std::move(p));
        }
        item["schedule"]=schedule.to_json();item["routes"]=std::move(routes);
        item["rims"]={rec.rims[0].to_json(),rec.rims[1].to_json()};
        records.push_back(std::move(item));
    }
    if(report)*report=std::move(records);
    return out;
}
ProposalGenerator endpoint_slide_generator() {
    return {"endpoint_slide","1",[](const RunState& state,ProposalContext& ctx) {return endpoint_slide_proposals(state,ctx,nullptr);}};
}
}
