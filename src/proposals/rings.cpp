// SPDX-License-Identifier: GPL-3.0-or-later
// Phase B generalized ring proposals (T35, Section 16): constriction and crease loops found
// from geometry alone (no object front, world-up or anatomy label), proposed as compound
// cuts that carry every companion opening both sides need, and judged by the common
// transaction. Used by policy enhanced_v1 only.
#include "cleanseams/proposals.hpp"
#include "cleanseams/physical_scale.hpp"
#include "cleanseams/cache.hpp"
#include "cleanseams/layout.hpp"
#include "cleanseams/errors.hpp"
#include <algorithm>
#include <cmath>
#include <map>
#include <numbers>
#include <numeric>
#include <queue>
#include <set>
namespace cleanseams {
namespace {
using Json=nlohmann::json;
constexpr const char* name="rings";
constexpr const char* version="1";
bool always_cut(const MeshTables& t,std::size_t e) {
    return t.edge_boundary[e] || !t.edge_manifold[e] || !t.edge_consistent[e] || t.edge_corners[e].size()!=2;
}
// Vertex components of an edge set (ascending edges per component, components by smallest edge).
std::vector<std::vector<MeshIndex>> edge_components(const MeshSnapshot& m,std::span<const MeshIndex> edges) {
    const auto& t=m.tables();std::map<MeshIndex,std::vector<MeshIndex>> at;
    for(auto e:edges)for(auto v:t.edge_vertices[e])at[v].push_back(e);
    std::set<MeshIndex> left(edges.begin(),edges.end());std::vector<std::vector<MeshIndex>> out;
    while(!left.empty()) {
        std::vector<MeshIndex> comp,stack{*left.begin()};left.erase(left.begin());
        while(!stack.empty()) {
            const auto e=stack.back();stack.pop_back();comp.push_back(e);
            for(auto v:t.edge_vertices[e])for(auto f:at[v])if(left.erase(f))stack.push_back(f);
        }
        std::sort(comp.begin(),comp.end());out.push_back(std::move(comp));
    }
    return out;
}
double edges_length(const MeshSnapshot& m,std::span<const MeshIndex> edges) {double s=0;for(auto e:edges)s+=m.geometry().edge_length[e];return s;}
}
std::vector<MeshIndex> simple_loop_order(const MeshSnapshot& m,std::span<const MeshIndex> edges) {
    const auto& t=m.tables();if(edges.size()<3)return {};
    std::map<MeshIndex,std::vector<MeshIndex>> at;
    for(auto e:edges)for(auto v:t.edge_vertices[e])at[v].push_back(e);
    if(at.size()!=edges.size())return {};
    for(const auto& [v,list]:at)if(list.size()!=2)return {};
    std::vector<MeshIndex> order;const auto start=t.edge_vertices[edges.front()][0];auto v=start;auto e=edges.front();
    std::set<MeshIndex> used;
    while(true) {
        order.push_back(v);used.insert(e);v=m.edge_other(e,v);
        if(v==start)break;
        const auto& list=at.at(v);e=list[0]==e?list[1]:list[0];
        if(used.contains(e))return {};
    }
    return used.size()==edges.size()?order:std::vector<MeshIndex>{};
}
nlohmann::json RingCandidate::to_json() const {
    return {{"edges",edges},{"vertices",vertices.size()},{"level",level},{"quantile",quantile},{"length",length},
            {"side_area",side_area},{"side_faces",side_faces},{"constriction",constriction},{"crease",crease},{"exposure",exposure},
            {"balance",balance},{"compactness",compactness},{"confidence",confidence},{"score",score},{"status",status}};
}
nlohmann::json RingSearch::to_json() const {
    // Reports stay bounded: the 16 best loops in full, the rest counted.
    Json loops_json=Json::array();for(std::size_t k=0;k<loops.size() && k<16;++k)loops_json.push_back(loops[k].to_json());
    return {{"version",ring_rank_version},{"opening_edges",opening.size()},{"opening_length",opening_length},{"far_length",far_length},
            {"max_distance",max_distance},{"levels",levels},{"section_length",section_length},{"non_simple_components",non_simple_components},
            {"loops_found",loops.size()},{"loops",loops_json},{"island_floor",island_floor},{"stopped",stopped}};
}
RingSearch ring_search(const MeshSnapshot& m,std::span<const MeshIndex> faces,const HardConstraints& hc,std::span<const double> feature,
                       std::span<const double> raw_exposure,const PipelineOptions& options,const PhysicalScale* scale,DetectorPass* pass) {
    const auto& t=m.tables();const auto& g=m.geometry();const auto& o=options.enhanced.rings;
    const auto nf=g.face_area.size(),ne=t.edge_vertices.size();
    if(!feature.empty() && feature.size()!=ne)throw ContractError(ErrorCode::invalid_input,"rings.feature","feature desirability needs one value per edge");
    if(!raw_exposure.empty() && raw_exposure.size()!=ne)throw ContractError(ErrorCode::invalid_input,"rings.exposure","raw exposure needs one value per edge");
    DetectorPass local_pass;DetectorPass& work=pass?*pass:local_pass;
    RingSearch out;
    // GB m7: everything below is local to the patch (sorted face, vertex and edge lists with
    // order-preserving local indices), never mesh-sized.
    std::vector<MeshIndex> patch(faces.begin(),faces.end());std::sort(patch.begin(),patch.end());patch.erase(std::unique(patch.begin(),patch.end()),patch.end());
    for(auto f:patch)if(f>=nf)throw ContractError(ErrorCode::invalid_input,"rings.faces","face outside the mesh");
    auto face_slot=[&](MeshIndex f) -> std::ptrdiff_t {
        const auto it=std::lower_bound(patch.begin(),patch.end(),f);return it!=patch.end() && *it==f?it-patch.begin():-1;
    };
    // Patch edges: interior (two patch faces; existing cuts ignored: physical distance) and
    // boundary (one patch face, or a physical/nonmanifold edge).
    std::vector<MeshIndex> patch_edges;for(auto f:patch)for(auto c:m.face_corners(f))patch_edges.push_back(t.corner_edge[c]);
    std::sort(patch_edges.begin(),patch_edges.end());patch_edges.erase(std::unique(patch_edges.begin(),patch_edges.end()),patch_edges.end());
    std::vector<MeshIndex> interior,boundary;std::vector<std::array<std::size_t,2>> interior_faces;
    for(auto e:patch_edges) {
        const auto pair=t.edge_faces[e];bool both=pair.size()==2;
        std::array<std::ptrdiff_t,2> slots{-1,-1};
        for(std::size_t k=0;k<pair.size() && k<2;++k)slots[k]=face_slot(pair[k]);
        both=both && slots[0]>=0 && slots[1]>=0;
        if(both) {interior.push_back(e);interior_faces.push_back({static_cast<std::size_t>(slots[0]),static_cast<std::size_t>(slots[1])});}
        else boundary.push_back(e);
    }
    if(boundary.empty() || interior.empty())return out;
    // The opening: the longest boundary component (ties: smallest edge).
    auto comps=edge_components(m,boundary);
    std::stable_sort(comps.begin(),comps.end(),[&](const auto& a,const auto& b) {return edges_length(m,a)>edges_length(m,b);});
    out.opening=comps.front();out.opening_length=edges_length(m,out.opening);
    for(std::size_t k=1;k<comps.size();++k)out.far_length+=edges_length(m,comps[k]);
    // Physical surface distance from the opening over the patch's edges (local vertices: the
    // sorted patch vertices, so local order is global order and the search is the same).
    std::vector<MeshIndex> vertices;for(auto e:patch_edges)for(auto v:t.edge_vertices[e])vertices.push_back(v);
    std::sort(vertices.begin(),vertices.end());vertices.erase(std::unique(vertices.begin(),vertices.end()),vertices.end());
    auto vertex_slot=[&](MeshIndex v) {return static_cast<std::size_t>(std::lower_bound(vertices.begin(),vertices.end(),v)-vertices.begin());};
    std::vector<std::vector<std::pair<std::size_t,MeshIndex>>> adjacent(vertices.size()); // in vert_edges order
    for(std::size_t i=0;i<vertices.size();++i)
        for(auto e:t.vert_edges[vertices[i]])
            if(std::binary_search(patch_edges.begin(),patch_edges.end(),e))adjacent[i].emplace_back(vertex_slot(m.edge_other(e,vertices[i])),e);
    if(work.stop()) {out.stopped=true;return out;}
    std::vector<double> dist(vertices.size(),INFINITY);using Item=std::pair<double,std::size_t>;std::priority_queue<Item,std::vector<Item>,std::greater<>> heap;
    for(auto e:out.opening)for(auto v:t.edge_vertices[e]) {const auto i=vertex_slot(v);if(dist[i]>0) {dist[i]=0;heap.push({0.0,i});}}
    while(!heap.empty()) {
        const auto [d,v]=heap.top();heap.pop();if(d>dist[v])continue;
        for(const auto& [w,e]:adjacent[v]) {
            const double nd=d+g.edge_length[e];
            if(nd<dist[w]) {dist[w]=nd;heap.push({nd,w});}
        }
    }
    std::vector<double> value(patch.size(),INFINITY);double top=0;
    for(std::size_t i=0;i<patch.size();++i) {
        double s=0;std::size_t n=0;for(auto c:m.face_corners(patch[i])) {s+=dist[vertex_slot(static_cast<MeshIndex>(m.input().corner_vertices[c]))];++n;}
        value[i]=n?s/static_cast<double>(n):INFINITY;if(std::isfinite(value[i]))top=std::max(top,value[i]);
    }
    out.max_distance=top;if(!(top>0))return out;
    const auto samples=static_cast<std::size_t>(std::max<std::int64_t>(1,o.section_samples));
    std::map<std::vector<MeshIndex>,std::pair<double,double>> found; // loop -> (level, quantile)
    for(std::size_t k=0;k<samples;++k) {
        if(work.stop()) {out.stopped=true;break;}
        const double q=samples==1?o.low_quantile:o.low_quantile+(o.high_quantile-o.low_quantile)*static_cast<double>(k)/static_cast<double>(samples-1);
        const double level=q*top;out.levels.push_back(level);
        std::vector<MeshIndex> cut;
        for(std::size_t i=0;i<interior.size();++i) {const auto [a,b]=interior_faces[i];if((value[a]<level)!=(value[b]<level))cut.push_back(interior[i]);}
        out.section_length.push_back(edges_length(m,cut));
        for(auto& comp:edge_components(m,cut)) {
            if(simple_loop_order(m,comp).empty()) {++out.non_simple_components;continue;}
            found.emplace(comp,std::pair{level,q});
        }
    }
    // Metrics of every distinct simple loop.
    const auto total_area=work.area(m);
    // T39: area-weighted region preservation: the minimum-island face count is an area floor on
    // the patch's reference chart under the physical scale (the face count without one).
    RegionFloor island;island.faces=options.layout.min_island_faces;island.physical=false;
    if(scale)island=scale->floor(options.layout.min_island_faces,faces);
    out.island_floor=island.to_json();
    // GB m11: the complete minimum-island policy (absorb.tiny_charts: face count with the large
    // n-gon exception, area fraction, thin strips). Under the physical scale its face counts are
    // area floors on the patch's reference chart (the strip rule's 50 faces too).
    const double area_floor=options.layout.min_island_fraction*total_area;
    std::optional<RegionFloor> strip_floor;if(island.physical)strip_floor=scale->floor(50,faces);
    auto tiny_side=[&](std::span<const MeshIndex> side_faces,double area,double perimeter) {
        if(!island.physical)return below_island_policy(m,side_faces,area,perimeter,area_floor,options.layout);
        const double iq=perimeter*perimeter/std::max(4*std::numbers::pi*area,1e-300);
        return island.below(area,side_faces.size()) || area<area_floor || (strip_floor->below(area,side_faces.size()) && iq>options.layout.strip_isoperimetric);
    };
    out.island_floor["strip_isoperimetric"]=options.layout.strip_isoperimetric;out.island_floor["area_fraction"]=options.layout.min_island_fraction;
    for(const auto& [edges,where]:found) {
        if(work.stop()) {out.stopped=true;break;}
        RingCandidate r;r.edges=edges;r.vertices=simple_loop_order(m,edges);r.level=where.first;r.quantile=where.second;r.length=edges_length(m,edges);
        // Sides: patch faces flooded across interior edges other than the loop.
        const std::set<MeshIndex> loop(edges.begin(),edges.end());std::vector<int> side(patch.size(),-1);int sides=0;
        for(auto f0:faces) {
            const auto s0=static_cast<std::size_t>(face_slot(f0));
            if(side[s0]>=0)continue;
            std::vector<MeshIndex> stack{f0};side[s0]=sides;
            while(!stack.empty()) {
                const auto f=stack.back();stack.pop_back();
                for(auto c:m.face_corners(f)) {
                    const auto e=t.corner_edge[c];if(loop.contains(e))continue;
                    const auto pair=t.edge_faces[e];if(pair.size()!=2)continue;
                    for(auto h:pair) {
                        if(h==f)continue;
                        const auto sh=face_slot(h);
                        if(sh>=0 && side[static_cast<std::size_t>(sh)]<0) {side[static_cast<std::size_t>(sh)]=sides;stack.push_back(h);}
                    }
                }
            }
            ++sides;
        }
        if(sides!=2) {r.status="non_separating";out.loops.push_back(std::move(r));continue;}
        std::array<std::vector<MeshIndex>,2> side_list;std::array<double,2> perimeter{};
        for(std::size_t i=0;i<patch.size();++i) {
            const auto s=static_cast<std::size_t>(side[i]);const auto f=patch[i];
            r.side_area[s]+=g.face_area[f];++r.side_faces[s];side_list[s].push_back(f);
            // The side's border as a region: the loop, the patch boundary and nonmanifold edges.
            for(auto c:m.face_corners(f)) {
                const auto e=t.corner_edge[c];const auto pair=t.edge_faces[e];bool inside=pair.size()==2 && !loop.contains(e);
                if(inside)for(auto h:pair)if(h!=f) {const auto sh=face_slot(h);inside=sh>=0 && side[static_cast<std::size_t>(sh)]==side[i];}
                if(!inside)perimeter[s]+=g.edge_length[e];
            }
        }
        // Constriction: shorter than the widest section on both of its sides (a waist).
        const auto k=static_cast<std::size_t>(std::find(out.levels.begin(),out.levels.end(),r.level)-out.levels.begin());
        double before=out.opening_length,after=out.far_length;
        for(std::size_t j=0;j<out.levels.size();++j) {if(j<k)before=std::max(before,out.section_length[j]);if(j>k)after=std::max(after,out.section_length[j]);}
        const double widest=std::min(before,after);
        r.constriction=widest>0?std::clamp(1.0-r.length/widest,0.0,1.0):0.0;
        double crease=0,exposure=0;
        for(auto e:edges) {
            const double l=g.edge_length[e];
            if(!feature.empty())crease+=l*feature[e];
            if(!raw_exposure.empty())exposure+=l*raw_exposure[e];
        }
        r.crease=r.length>0?crease/r.length:0.0;r.exposure=r.length>0?exposure/r.length:0.0;
        const double small=std::min(r.side_area[0],r.side_area[1]);
        r.balance=small/std::max(r.side_area[0]+r.side_area[1],1e-300);
        r.compactness=r.length>0?std::min(1.0,std::sqrt(small)/r.length):0.0;
        r.confidence=std::max(r.constriction,r.crease);
        r.score=r.constriction+0.5*r.crease+0.25*r.balance+0.25*r.compactness-0.5*r.exposure;
        // The minimum-island policy is never relaxed: tiny sides are refused.
        bool tiny=false;
        for(std::size_t s=0;s<2;++s)tiny=tiny || tiny_side(side_list[s],r.side_area[s],perimeter[s]);
        bool forbidden=false;for(auto e:edges)forbidden=forbidden || (e<hc.forbidden.size() && hc.forbidden[e] && !always_cut(t,e));
        r.status=forbidden?"forbidden_edge":tiny?"tiny_side":r.confidence<o.min_confidence?"low_confidence":"candidate";
        out.loops.push_back(std::move(r));
    }
    // Ranked: candidates first, best score first (ties: smallest first edge).
    std::stable_sort(out.loops.begin(),out.loops.end(),[](const RingCandidate& a,const RingCandidate& b) {
        const bool ca=a.status=="candidate",cb=b.status=="candidate";
        if(ca!=cb)return ca;
        if(a.score!=b.score)return a.score>b.score;
        return a.edges<b.edges;
    });
    return out;
}
ProposalGenerator ring_generator() {
    return {name,version,[](const RunState& state,ProposalContext& ctx) {
        std::vector<Proposal> out;
        const auto& o=ctx.options.enhanced.rings;
        if(!o.enabled || !state.surface)return out;
        const auto& mesh=ctx.mesh;const auto& m=*mesh;const auto& t=m.tables();const auto nf=m.geometry().face_area.size();
        const auto& ct=state.surface->tables();const auto& hc=ctx.constraints;
        auto hard=[&](std::size_t e) {return (e<hc.locked.size() && hc.locked[e]) || (e<hc.required.size() && hc.required[e]);};
        const auto prefs=context_preferences(ctx);const auto routing=context_routing_cost(ctx);
        const auto& scale=context_scale(ctx,state);
        // GB m7: mesh-level quantities once per pass; the pass's control inside the searches.
        DetectorPass work;work.control.cancel=ctx.cancel;work.control.deadline=ctx.deadline;bool stopped=false;
        struct Chart {MeshIndex chart;RingSearch search;double best=-INFINITY;};
        std::vector<Chart> charts;Json records=Json::array();
        for(MeshIndex ch=0;ch<ct.n_charts;++ch) {
            const auto& c=state.charts[ch];
            if(c.exempt_only || !c.measured)continue;
            const auto faces=ct.chart_faces[ch];
            // Two sides of at least the minimum island each (an area floor under the physical scale).
            double area=0;for(auto f:faces)area+=m.geometry().face_area[f];
            if(scale.floor(2*std::max<std::int64_t>(1,ctx.options.layout.min_island_faces),faces).below(area,faces.size()))continue;
            if(work.stop()) {stopped=true;break;}
            Chart entry{ch,ring_search(m,faces,hc,prefs->feature_desirability,ctx.raw_exposure,ctx.options,&scale,&work)};
            for(const auto& l:entry.search.loops)if(l.status=="candidate")entry.best=std::max(entry.best,l.score);
            const bool partial=entry.search.stopped;
            charts.push_back(std::move(entry));
            if(partial) {stopped=true;break;}
        }
        std::stable_sort(charts.begin(),charts.end(),[](const Chart& a,const Chart& b) {return a.best>b.best;});
        std::uint64_t id=0,charts_used=0;std::set<EdgeMask> emitted;
        for(auto& entry:charts) {
            Json record={{"chart",entry.chart},{"faces",ct.chart_faces[entry.chart].size()},{"loops_found",entry.search.loops.size()},
                         {"best_score",std::isfinite(entry.best)?Json(entry.best):Json(nullptr)},{"proposals",Json::array()}};
            if(!std::isfinite(entry.best) || static_cast<std::int64_t>(charts_used)>=o.max_charts) {
                // Charts without a candidate keep a short record (their loops' reasons only
                // for the charts that are searched in full below, to bound the report).
                record["status"]=std::isfinite(entry.best)?"beyond_chart_budget":"no_candidate";
                if(records.size()<16) {
                    Json reasons=Json::object();for(const auto& l:entry.search.loops)reasons[l.status]=reasons.value(l.status,0)+1;
                    record["loop_status"]=reasons;
                }
                records.push_back(std::move(record));continue;
            }
            ++charts_used;
            const auto faces=ct.chart_faces[entry.chart];std::vector<std::uint8_t> scope(nf,0);for(auto f:faces)scope[f]=1;
            // The underlying patch: the chart without its non-hard interior cuts; the old
            // openings are free for the companion paths, so they are reused when they fit.
            // (GB m7: the chart's own edges only, ascending.)
            auto patch=state.cuts;auto cost=routing;std::vector<MeshIndex> dropped;
            std::vector<MeshIndex> chart_edges;for(auto f:faces)for(auto c:m.face_corners(f))chart_edges.push_back(t.corner_edge[c]);
            std::sort(chart_edges.begin(),chart_edges.end());chart_edges.erase(std::unique(chart_edges.begin(),chart_edges.end()),chart_edges.end());
            for(auto e:chart_edges) {
                const auto pair=t.edge_faces[e];
                if(!state.cuts[e] || always_cut(t,e) || hard(e) || pair.size()!=2 || !scope[pair[0]] || !scope[pair[1]])continue;
                patch[e]=0;cost[e]=0.0;dropped.push_back(e);
            }
            std::int64_t used=0;
            for(auto& loop:entry.search.loops) {
                if(loop.status!="candidate")continue;
                if(used>=o.max_candidates) {loop.status="beyond_candidate_budget";continue;}
                ++used;
                auto cuts=patch;for(auto e:loop.edges)cuts[e]=1;
                const auto opening=open_non_disk_charts(mesh,cuts,cost,hc,ctx.exempt,scope,work.control);
                if(opening.stopped) {stopped=true;loop.status="stopped";break;}
                Json p={{"loop_edges",loop.edges.size()},{"quantile",loop.quantile},{"score",loop.score},{"confidence",loop.confidence},{"openings",opening.added}};
                if(!opening.ok) {loop.status="not_openable";p["status"]="refused";p["failures"]=opening.failures;record["proposals"].push_back(p);continue;}
                const auto violations=constraint_violations(m,hc,cuts);
                if(!violations.empty()) {loop.status="constraint";p["status"]="refused";p["violations"]=violations;record["proposals"].push_back(p);continue;}
                if(cuts==state.cuts || !emitted.insert(cuts).second) {p["status"]="duplicate";record["proposals"].push_back(p);continue;}
                ++id;
                Json provenance={{"chart",entry.chart},{"loop",loop.edges},{"quantile",loop.quantile},{"level",loop.level},{"length",loop.length},
                                 {"constriction",loop.constriction},{"crease",loop.crease},{"exposure",loop.exposure},{"balance",loop.balance},
                                 {"confidence",loop.confidence},{"score",loop.score},{"rank_version",ring_rank_version},{"openings",opening.added},
                                 {"dropped_interior_cuts",dropped}};
                out.push_back(cut_proposal(id,state,cuts,name,version,std::move(provenance)));
                loop.status="proposed";p["status"]="proposed";p["id"]=id;record["proposals"].push_back(p);
            }
            record["search"]=entry.search.to_json();record["status"]="searched";
            records.push_back(std::move(record));
            if(stopped)break;
        }
        ctx.diagnostics[name]={{"version",version},{"rank_version",ring_rank_version},
                               {"options",{{"section_samples",o.section_samples},{"low_quantile",o.low_quantile},{"high_quantile",o.high_quantile},
                                           {"max_candidates",o.max_candidates},{"max_charts",o.max_charts},{"min_confidence",o.min_confidence}}},
                               {"charts",records},{"proposals",out.size()},{"stopped_by_deadline",stopped}};
        return out;
    }};
}
}
