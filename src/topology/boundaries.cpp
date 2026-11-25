// SPDX-License-Identifier: GPL-3.0-or-later
#include "detail.hpp"
#include "cleanseams/detail/numpy_compat.hpp"
#include "cleanseams/endpoints.hpp"
#include <map>
namespace cleanseams {
using namespace cut_detail;
void CutSurface::build_boundaries(CutWorkspace& workspace) {
    auto& d=data_;const auto& m=*mesh_;
    const auto outgoing=rows(d.n_sectors,[&](auto add){for(auto c:d.boundary_corners) add(d.sector_of_corner[c],c);});
    std::vector<Offset> incoming(d.n_sectors,0);
    for(auto c:d.boundary_corners) ++incoming[d.sector_of_corner[m.corner_next(c)]];
    for(MeshIndex s=0;s<d.n_sectors;++s) if(outgoing[s].size()>1 || incoming[s]>1) d.pinched_sectors.push_back(s);
    auto& visited=workspace.visited_;
    for(auto first:d.boundary_corners) {
        if(visited[first]) continue;
        auto current=first;d.loop_chart.push_back(d.chart_of_corner[first]);
        while(!visited[current]) {
            visited[current]=1;d.loops.values.push_back(current);
            const auto candidates=outgoing[d.sector_of_corner[m.corner_next(current)]];
            if(candidates.size()==1) current=candidates[0];
            else {
                const auto next=std::find_if(candidates.begin(),candidates.end(),[&](auto c){return !visited[c];});
                if(next==candidates.end()) break;current=*next;
            }
        }
        d.loops.offsets.push_back(d.loops.values.size());
    }
    d.chart_loops=rows(d.n_charts,[&](auto add){for(std::size_t i=0;i<d.loop_chart.size();++i) add(d.loop_chart[i],static_cast<MeshIndex>(i));});
}
IndexRows split_arcs(const MeshSnapshot& mesh,std::span<const MeshIndex> edges,std::span<const Offset> degree) {
    require(degree.size()==mesh.input().positions.size(),"cut_degree","vertex count mismatch");
    require(edges.size()<=mesh.tables().edge_vertices.size(),"arc.edges","too many edges");
    std::map<MeshIndex,std::vector<MeshIndex>> incident;std::set<MeshIndex> unique,used;IndexRows arcs;
    for(auto e:edges) {
        require(e<mesh.tables().edge_vertices.size() && unique.insert(e).second,"arc.edges","invalid or repeated edge");
        for(auto v:mesh.tables().edge_vertices[e]) incident[v].push_back(e);
    }
    auto walk=[&](MeshIndex vertex,MeshIndex edge) {
        while(true) {
            used.insert(edge);arcs.values.push_back(edge);vertex=mesh.edge_other(edge,vertex);
            if(degree[vertex]!=2) break;
            const auto found=incident.find(vertex);MeshIndex next=0;std::size_t count=0;
            if(found!=incident.end()) for(auto candidate:found->second) if(candidate!=edge && !used.contains(candidate)) {next=candidate;++count;}
            if(count!=1) break;edge=next;
        }
        arcs.offsets.push_back(arcs.values.size());
    };
    for(auto& [vertex,local]:incident) {
        if(degree[vertex]==2) continue;std::sort(local.begin(),local.end());
        for(auto e:local) if(!used.contains(e)) walk(vertex,e);
    }
    for(auto e:unique) if(!used.contains(e)) walk(mesh.tables().edge_vertices[e][0],e);
    return arcs;
}
namespace {
std::set<MeshIndex> touched_vertices(const CutSurface& cs,std::span<const MeshIndex> edges) {
    require(edges.size()<=cs.mesh().tables().edge_vertices.size(),"arc.edges","too many edges");
    std::set<MeshIndex> vertices;
    for(auto edge:edges) {
        require(edge<cs.mesh().tables().edge_vertices.size(),"arc.edges","edge out of range");
        for(auto v:cs.mesh().tables().edge_vertices[edge]) vertices.insert(v);
    }
    return vertices;
}
}
ArcKey arc_key(const CutSurface& cs,std::span<const MeshIndex> edges) {
    ArcKey result;for(auto v:touched_vertices(cs,edges)) if(cs.tables().cut_degree[v]!=2) result.ends.push_back(v);
    result.edges.assign(edges.begin(),edges.end());std::sort(result.edges.begin(),result.edges.end());
    result.edges.erase(std::unique(result.edges.begin(),result.edges.end()),result.edges.end());
    result.hash=detail::python_frozenset_hash(result.edges);return result;
}
Offset free_ends(const CutSurface& cs,std::span<const MeshIndex> edges) {
    Offset count=0;for(auto v:touched_vertices(cs,edges)) if(cs.tables().cut_degree[v]==1) ++count;return count;
}
// ============================================================ Phase B: boundary-side provenance and rims (T30)
// Section 13.1: the rims of a slit tube are recovered on the underlying patch (the slit glued),
// never from the cut chart's boundary count: a slit cylinder has one disk boundary holding both
// physical rims and the two slit sides. Every boundary side keeps why it is cut.
namespace {
bool forced_cut(const MeshTables& t,MeshIndex e) {
    return !t.edge_boundary[e] && (!t.edge_manifold[e] || !t.edge_consistent[e] || t.edge_corners[e].size()!=2);
}
}
std::string_view boundary_side_kind_name(BoundarySideKind k) {
    switch(k) {case BoundarySideKind::physical:return "physical";case BoundarySideKind::forced:return "forced";default:return "seam";}
}
nlohmann::json BoundarySide::to_json() const {
    return {{"corner",corner},{"edge",edge},{"loop",loop},{"kind",boundary_side_kind_name(kind)},
            {"twin",twin?nlohmann::json(*twin):nlohmann::json(nullptr)},{"twin_chart",twin_chart?nlohmann::json(*twin_chart):nlohmann::json(nullptr)}};
}
std::vector<BoundarySide> boundary_sides(const CutSurface& cs,MeshIndex chart) {
    const auto& ct=cs.tables();const auto& t=cs.mesh().tables();
    require(chart<ct.n_charts,"boundary.chart","chart out of range");
    std::vector<BoundarySide> out;
    for(auto loop:ct.chart_loops[chart])for(auto c:ct.loops[loop]) {
        BoundarySide s;s.corner=c;s.edge=t.corner_edge[c];s.loop=loop;
        if(t.edge_boundary[s.edge])s.kind=BoundarySideKind::physical;
        else if(forced_cut(t,s.edge))s.kind=BoundarySideKind::forced;
        else {
            s.kind=BoundarySideKind::seam;
            for(auto other:t.edge_corners[s.edge])if(other!=c) {s.twin=other;s.twin_chart=ct.chart_of_corner[other];}
        }
        out.push_back(s);
    }
    return out;
}
std::optional<std::size_t> PhysicalRim::position(MeshIndex vertex) const {
    const auto found=std::find(vertices.begin(),vertices.end(),vertex);
    if(found==vertices.end())return std::nullopt;
    return static_cast<std::size_t>(found-vertices.begin());
}
double PhysicalRim::distance(MeshIndex a,MeshIndex b) const {
    const auto i=position(a),j=position(b);
    require(i && j,"rim.vertex","vertex is not on the rim");
    const double d=std::abs(arc[*i]-arc[*j]);
    return std::min(d,length-d);
}
nlohmann::json PhysicalRim::to_json() const {
    return {{"id",id.to_json()},{"kind",kind},{"vertices",vertices},{"edges",edges},{"length",length}};
}
nlohmann::json PatchRims::to_json() const {
    nlohmann::json r=nlohmann::json::array();for(const auto& rim:rims)r.push_back(rim.to_json());
    return {{"annulus",annulus},{"reason",reason},{"patch_chart",patch_chart},{"euler",euler},{"loops",loops},{"rims",r}};
}
PatchRims patch_rims(const std::shared_ptr<const CutSurface>& surface,MeshIndex chart,std::span<const MeshIndex> slit_edges) {
    require(bool(surface),"rims.surface","surface is required");
    const auto& cs=*surface;const auto& ct=cs.tables();const auto& mesh=cs.mesh();const auto& t=mesh.tables();const auto& g=mesh.geometry();
    require(chart<ct.n_charts,"rims.chart","chart out of range");
    // The underlying patch: the same cuts with the slit glued (a temporary analysis surface;
    // the state's cuts are never modified).
    auto cut=ct.cut;
    for(auto e:slit_edges) {
        require(e<cut.size() && cut[e] && !t.edge_boundary[e] && !forced_cut(t,e),"rims.slit","slit edges must be internal cuts");
        for(auto f:t.edge_faces[e])require(ct.chart_of_face[f]==chart,"rims.slit","slit edges must lie inside the chart");
        cut[e]=0;
    }
    PatchRims out;out.patch=std::make_shared<const CutSurface>(surface->mesh_owner(),cut);
    const auto& pt=out.patch->tables();
    const auto faces=ct.chart_faces[chart];
    out.patch_chart=pt.chart_of_face[faces[0]];
    // Gluing edges inside one chart keeps its face set.
    require(pt.chart_faces[out.patch_chart].size()==faces.size(),"rims.patch","gluing the slit changed the chart's faces");
    out.euler=pt.chart_euler[out.patch_chart];out.loops=pt.chart_n_loops[out.patch_chart];
    auto fail=[&](const char* reason) {out.reason=reason;return out;};
    if(pt.chart_unsupported[out.patch_chart])return fail("unsupported_patch");
    if(pt.chart_pinched[out.patch_chart])return fail("pinched_patch");
    if(out.loops==0)return fail("closed_patch");
    if(out.loops==1)return fail("not_tube");
    if(out.loops>2)return fail("extra_holes");
    if(out.euler!=0)return fail(out.euler<0?"handle":"not_tube");
    for(auto loop:pt.chart_loops[out.patch_chart]) {
        PhysicalRim rim;bool physical=false,seam=false;
        for(auto c:pt.loops[loop]) {
            const auto e=t.corner_edge[c];
            if(forced_cut(t,e))return fail("forced_rim");
            if(t.edge_boundary[e])physical=true;
            else {
                // A seam side of a patch rim faces another chart: a side whose twin lies in the
                // patch itself is a remaining slit, not a rim.
                for(auto other:t.edge_corners[e])if(other!=c && pt.chart_of_corner[other]==out.patch_chart)return fail("slit_side_rim");
                seam=true;
            }
            rim.vertices.push_back(static_cast<MeshIndex>(mesh.input().corner_vertices[c]));rim.edges.push_back(e);rim.corners.push_back(c);
            rim.arc.push_back(rim.length);rim.length+=g.edge_length[e];
        }
        auto sorted=rim.vertices;std::sort(sorted.begin(),sorted.end());
        if(std::adjacent_find(sorted.begin(),sorted.end())!=sorted.end())return fail("pinched_rim");
        rim.kind=physical && seam?"mixed":physical?"physical":"seam";
        rim.id={*std::min_element(rim.edges.begin(),rim.edges.end()),static_cast<MeshIndex>(rim.edges.size())};
        out.rims.push_back(std::move(rim));
    }
    std::vector<MeshIndex> a=out.rims[0].vertices,b=out.rims[1].vertices;std::sort(a.begin(),a.end());std::sort(b.begin(),b.end());
    std::vector<MeshIndex> shared;std::set_intersection(a.begin(),a.end(),b.begin(),b.end(),std::back_inserter(shared));
    if(!shared.empty())return fail("rims_intersect");
    std::sort(out.rims.begin(),out.rims.end(),[](const PhysicalRim& x,const PhysicalRim& y) {return x.id<y.id;});
    out.annulus=true;
    return out;
}
}
