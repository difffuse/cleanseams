// SPDX-License-Identifier: GPL-3.0-or-later
// Phase B (T32, Section 14): sector-aware hotspot relief proposals, evaluated before the
// baseline's broad over-budget split. Only candidate construction lives here; every
// candidate is an ordinary T28 proposal (topology, solve, validity, strict overlap, pins,
// budgets, run-origin envelope, common ordering, atomic commit). Used by enhanced_v1 only.
#include "cleanseams/relief.hpp"
#include "cleanseams/physical_scale.hpp"
#include "cleanseams/cache.hpp"
#include "cleanseams/errors.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <numeric>
#include <queue>
#include <set>
namespace cleanseams {
namespace {
using Json=nlohmann::json;
constexpr double infinity=std::numeric_limits<double>::infinity();
constexpr std::size_t none=std::numeric_limits<std::size_t>::max();
Json finite_or_null(double x) {return std::isfinite(x)?Json(x):Json(nullptr);}
bool always_cut(const MeshTables& t,std::size_t e) {
    return t.edge_boundary[e] || !t.edge_manifold[e] || !t.edge_consistent[e] || t.edge_corners[e].size()!=2;
}
// The sectors of one chart and its surface graph: every polygon edge of the chart's faces
// joins the sectors of its two corners (a cut edge joins the sectors of one side only), so
// walks never cross a cut to the other side or reach another chart through a shared
// original vertex.
struct ChartGraph {
    std::vector<MeshIndex> sectors;                 // local index -> sector (ascending)
    std::vector<std::vector<std::pair<std::size_t,MeshIndex>>> adjacent; // local -> (local, edge)
    std::vector<std::uint8_t> interior;             // vertex with no incident effective cut and one sector
    std::vector<std::uint8_t> seam_side;            // on a side of an internal seam (not a physical boundary)
    // Sector -> local index (none outside the chart).
    std::size_t local(MeshIndex s) const {
        const auto it=std::lower_bound(sectors.begin(),sectors.end(),s);
        return it!=sectors.end() && *it==s?static_cast<std::size_t>(it-sectors.begin()):none;
    }
};
ChartGraph chart_graph(const CutSurface& cs,MeshIndex chart) {
    const auto& m=cs.mesh();const auto& t=m.tables();const auto& ct=cs.tables();
    // GB m7: chart-local: the chart's sectors are those of its faces' corners (ascending, as a
    // scan of every sector gave them); the cut degree and the sector count of a vertex are
    // read per vertex, never over the whole mesh.
    ChartGraph g;
    for(auto f:ct.chart_faces[chart])for(auto c:m.face_corners(f))g.sectors.push_back(ct.sector_of_corner[c]);
    std::sort(g.sectors.begin(),g.sectors.end());g.sectors.erase(std::unique(g.sectors.begin(),g.sectors.end()),g.sectors.end());
    g.adjacent.resize(g.sectors.size());g.interior.assign(g.sectors.size(),0);g.seam_side.assign(g.sectors.size(),0);
    for(std::size_t i=0;i<g.sectors.size();++i) {
        const auto v=ct.sector_vert[g.sectors[i]];
        g.interior[i]=ct.cut_degree[v]==0 && cs.sectors_at_vertex(v).size()==1;
    }
    for(auto f:ct.chart_faces[chart])
        for(auto c:m.face_corners(f)) {
            const auto e=t.corner_edge[c];const auto a=g.local(ct.sector_of_corner[c]),b=g.local(ct.sector_of_corner[m.corner_next(c)]);
            if(a==none || b==none || a==b)continue;
            g.adjacent[a].emplace_back(b,e);g.adjacent[b].emplace_back(a,e);
            if(ct.cut[e] && !always_cut(t,e)) {g.seam_side[a]=1;g.seam_side[b]=1;}
        }
    for(auto& list:g.adjacent) {std::sort(list.begin(),list.end());list.erase(std::unique(list.begin(),list.end()),list.end());}
    return g;
}
// Surface distances (3D edge lengths along the chart's graph) from `sources`, up to `limit`.
std::vector<double> surface_distance(const ChartGraph& g,const MeshGeometry& geo,std::span<const std::size_t> sources,double limit=infinity) {
    std::vector<double> d(g.sectors.size(),infinity);
    using Item=std::pair<double,std::size_t>;std::priority_queue<Item,std::vector<Item>,std::greater<>> queue;
    for(auto s:sources) {d[s]=0;queue.emplace(0.0,s);}
    while(!queue.empty()) {
        const auto [du,u]=queue.top();queue.pop();
        if(du>d[u])continue;
        for(const auto& [w,e]:g.adjacent[u]) {
            const double dw=du+geo.edge_length[e];
            if(dw<d[w] && dw<=limit) {d[w]=dw;queue.emplace(dw,w);}
        }
    }
    return d;
}
struct Admissibility {
    std::vector<std::uint8_t> traversable;          // local: interior sectors a slit may pass through or end at
    std::vector<std::uint8_t> start;                // local: boundary sectors a slit may start from
    const CutTables* cut=nullptr;const MeshTables* mesh=nullptr;const HardConstraints* hc=nullptr;
    // An uncut, non-forbidden edge (GB m7: read per edge, no mesh-sized table per chart).
    bool edge_ok(MeshIndex e) const {return !cut->cut[e] && !always_cut(*mesh,e) && !(hc && e<hc->forbidden.size() && hc->forbidden[e]);}
};
// A slit never runs along a cut or a forbidden edge, never passes through or ends at a
// vertex with a pinned corner (it would detach the pinned corner's neighbours on one side),
// and never starts at the free tip of a locked seam (a relief slit does not lengthen a
// locked seam).
Admissibility admissibility(const CutSurface& cs,const ChartGraph& g,const HardConstraints* hc) {
    const auto& m=cs.mesh();const auto& t=m.tables();const auto& ct=cs.tables();
    Admissibility a;a.traversable.assign(g.sectors.size(),0);a.start.assign(g.sectors.size(),0);a.cut=&ct;a.mesh=&t;a.hc=hc;
    std::set<MeshIndex> pinned;
    if(hc)for(const auto& [corner,uv]:hc->pins)if(corner<t.corner_face.size())pinned.insert(static_cast<MeshIndex>(m.input().corner_vertices[corner]));
    for(std::size_t i=0;i<g.sectors.size();++i) {
        const auto v=ct.sector_vert[g.sectors[i]];
        if(g.interior[i]) {a.traversable[i]=!pinned.contains(v);continue;}
        bool locked_tip=false;
        if(hc && !t.vert_boundary[v]) {
            std::uint32_t cuts=0;bool locked=false;
            for(auto e:t.vert_edges[v])if(ct.cut[e]) {++cuts;locked=locked || (e<hc->locked.size() && hc->locked[e]);}
            locked_tip=cuts==1 && locked;
        }
        a.start[i]=!locked_tip;
    }
    return a;
}
struct Route {std::vector<std::size_t> nodes;std::vector<MeshIndex> edges;double cost=0;};
// Cheapest slit over admissible edges and interior sectors, from a fixed boundary start
// (`start` set) or from any eligible boundary start, ending at `tip`. Nodes boundary first.
std::optional<Route> slit_route(const ChartGraph& g,const Admissibility& a,std::span<const double> cost,std::size_t tip,std::optional<std::size_t> start) {
    const auto n=g.sectors.size();std::vector<double> d(n,infinity);std::vector<std::size_t> previous(n,none);std::vector<MeshIndex> via(n,0);
    std::vector<std::uint8_t> done(n,0);
    using Item=std::pair<double,std::size_t>;std::priority_queue<Item,std::vector<Item>,std::greater<>> queue;
    // Search from the tip outward (symmetric costs); a boundary sector is a terminal.
    d[tip]=0;queue.emplace(0.0,tip);
    std::size_t reached=none;
    while(!queue.empty()) {
        const auto [du,u]=queue.top();queue.pop();
        if(done[u] || du>d[u])continue;
        done[u]=1;
        if(!g.interior[u]) {reached=u;break;} // the first terminal popped is the cheapest
        for(const auto& [w,e]:g.adjacent[u]) {
            if(!a.edge_ok(e) || done[w])continue;
            if(g.interior[w]) {if(!a.traversable[w])continue;}
            else if(start ? w!=*start : !a.start[w])continue;
            const double dw=du+cost[e];
            if(dw<d[w] || (dw==d[w] && u<previous[w])) {d[w]=dw;previous[w]=u;via[w]=e;queue.emplace(dw,w);}
        }
    }
    if(reached==none)return std::nullopt;
    Route r;r.cost=d[reached];
    for(auto x=reached;x!=none;x=previous[x]) {r.nodes.push_back(x);if(previous[x]!=none)r.edges.push_back(via[x]);}
    return r; // nodes: boundary start ... tip; edges in the same order
}
double suitability(const MeshSnapshot& m,MeshIndex v,std::span<const double> cost) {
    const auto& g=m.geometry();double sum=0;std::size_t n=0;
    for(auto e:m.tables().vert_edges[v]) {
        const double c=cost.empty()?g.edge_length[e]:cost[e];
        if(g.edge_length[e]>0 && c>0 && std::isfinite(c)) {sum+=g.edge_length[e]/c;++n;}
    }
    return n?sum/static_cast<double>(n):1.0;
}
void check_settings(const ReliefSettings& s) {
    auto require=[](bool ok,const char* field,const char* what) {if(!ok)throw ContractError(ErrorCode::invalid_input,field,what);};
    require(s.max_hot_regions>=0 && s.max_candidates_per_pass>=0 && s.passes>=0 && s.targets_per_region>=0,"relief.settings","budgets must be nonnegative");
    for(double r:s.reaches)require(std::isfinite(r) && r>0 && r<=1,"relief.reaches","reaches are fractions in (0, 1]");
    require(std::isfinite(s.log_area_threshold) && s.log_area_threshold>=0 && std::isfinite(s.log_anisotropy_threshold) && s.log_anisotropy_threshold>=0,
            "relief.thresholds","thresholds must be finite and nonnegative");
    require(std::isfinite(s.separation) && s.separation>0,"relief.separation","separation must be positive");
    require(s.smoothing_steps>=0,"relief.smoothing_steps","smoothing steps must be nonnegative");
}
}
nlohmann::json ReliefSettings::to_json() const {
    return {{"version",hotspot_relief_version},{"max_hot_regions",max_hot_regions},{"max_candidates_per_pass",max_candidates_per_pass},{"passes",passes},
            {"reaches",reaches},{"targets_per_region",targets_per_region},{"log_area_threshold",log_area_threshold},
            {"log_anisotropy_threshold",log_anisotropy_threshold},{"smoothing_steps",smoothing_steps},{"separation",separation},
            {"units",{{"log_area_threshold","area (natural log of the area ratio)"},{"log_anisotropy_threshold","shape (natural log of s0/s1)"},
                      {"separation","fraction of sqrt(chart 3D area), measured as surface distance"},{"reaches","fraction of the physical path length"}}}};
}
nlohmann::json ReliefPlan::to_json() const {
    Json regions_json=Json::array(),targets_json=Json::array(),paths_json=Json::array();
    for(const auto& r:regions)regions_json.push_back({{"sector",r.sector},{"vertex",r.vertex},{"boundary",r.boundary},{"excess",r.excess},
        {"weighted_excess",r.weighted_excess},{"log_area",finite_or_null(r.log_area)},{"log_anisotropy",finite_or_null(r.log_anisotropy)},{"mass",r.mass}});
    for(const auto& t:targets)targets_json.push_back({{"region",t.region},{"rank",t.rank},{"tip_sector",t.tip_sector},{"tip_vertex",t.tip_vertex},
        {"start_sector",t.start_sector?Json(*t.start_sector):Json(nullptr)},{"score",finite_or_null(t.score)},{"excess",finite_or_null(t.excess)},
        {"suitability",finite_or_null(t.suitability)},{"seam_distance",finite_or_null(t.seam_distance)},{"depth",finite_or_null(t.depth)}});
    for(const auto& p:paths)paths_json.push_back({{"target",p.target},{"edges",p.edges},{"vertices",p.vertices},{"length",p.length},
        {"available",p.available},{"reaches",p.reaches}});
    return {{"chart",chart},{"field",field},{"chart_area",chart_area},{"separation",separation},{"regions",regions_json},{"targets",targets_json},
            {"paths",paths_json},{"generated",generated},{"duplicates",duplicates},{"shared_tips",shared_tips},{"truncated",truncated},
            {"unreachable",unreachable},{"skipped",skipped}};
}
ReliefPlan plan_relief(const ReliefInput& in,const ReliefSettings& settings) {
    check_settings(settings);
    if(!in.surface)throw ContractError(ErrorCode::invalid_input,"relief.surface","missing cut surface");
    const auto& cs=*in.surface;const auto& m=cs.mesh();const auto& t=m.tables();const auto& geo=m.geometry();const auto& ct=cs.tables();
    const auto ne=t.edge_vertices.size();
    if(in.chart>=ct.n_charts)throw ContractError(ErrorCode::invalid_input,"relief.chart","chart out of range");
    if(in.uv.size()!=t.corner_face.size())throw ContractError(ErrorCode::invalid_input,"relief.uv","one UV per corner required");
    if(!in.edge_cost.empty() && in.edge_cost.size()!=ne)throw ContractError(ErrorCode::invalid_input,"relief.edge_cost","one cost per edge required");
    // GB m7: only the chart's edges are routed, so only their costs are checked.
    if(!in.edge_cost.empty())
        for(auto f:ct.chart_faces[in.chart])for(auto c:m.face_corners(f))
            if(!(in.edge_cost[t.corner_edge[c]]>=0))throw ContractError(ErrorCode::invalid_input,"relief.edge_cost","costs must be nonnegative");
    ReliefPlan plan;plan.chart=in.chart;
    for(auto f:ct.chart_faces[in.chart])plan.chart_area+=geo.face_area[f];
    plan.separation=settings.separation*std::sqrt(plan.chart_area);
    if(!ct.chart_is_disk[in.chart]) {plan.skipped="not_disk";return plan;}
    // Sector field of this chart, with the run's exemptions.
    // (GB m7: the chart's own exempt faces; the snapshot measures this chart only.)
    QualityPolicy policy;for(auto f:ct.chart_faces[in.chart])if(f<in.exempt.size() && in.exempt[f])policy.exempt_faces.push_back(static_cast<TransferIndex>(f));
    const std::array<MeshIndex,1> one{in.chart};
    const QualitySnapshot snapshot(quality_geometry(in.surface),CornerUvView{in.uv},policy,std::span<const MeshIndex>(one));
    const auto field=sector_distortion_field(snapshot,settings.smoothing_steps);
    plan.field=field.summary();plan.field["chart_mass"]=field.chart_mass[in.chart];plan.field["chart_measured_area"]=field.chart_area[in.chart];
    if(!(field.chart_area[in.chart]>0)) {plan.skipped="unmeasured";return plan;}
    const std::span<const double> cost=in.edge_cost.empty()?std::span<const double>(geo.edge_length):in.edge_cost;
    const auto g=chart_graph(cs,in.chart);const auto a=admissibility(cs,g,in.constraints);const auto n=g.sectors.size();
    // Distortion excess per sector (area and shape parts, each over its threshold).
    std::vector<double> excess(n,0.0);
    for(std::size_t i=0;i<n;++i) {
        const auto s=g.sectors[i];if(!(field.mass[s]>0))continue;
        excess[i]=std::max(0.0,std::abs(field.log_area[s])-settings.log_area_threshold)+std::max(0.0,field.log_anisotropy[s]-settings.log_anisotropy_threshold);
    }
    std::vector<std::size_t> boundary_nodes,seam_nodes;
    for(std::size_t i=0;i<n;++i) {if(!g.interior[i])boundary_nodes.push_back(i);if(g.seam_side[i])seam_nodes.push_back(i);}
    const auto depth=surface_distance(g,geo,boundary_nodes);
    const auto seam_distance=surface_distance(g,geo,seam_nodes);
    // Hot regions: greedy by excess, separated by surface distance.
    std::vector<std::size_t> order;for(std::size_t i=0;i<n;++i)if(excess[i]>0)order.push_back(i);
    std::stable_sort(order.begin(),order.end(),[&](std::size_t x,std::size_t y) {return excess[x]!=excess[y]?excess[x]>excess[y]:g.sectors[x]<g.sectors[y];});
    std::vector<std::uint8_t> covered(n,0);std::vector<std::vector<double>> region_distance;
    for(auto c:order) {
        if(static_cast<std::int64_t>(plan.regions.size())>=settings.max_hot_regions)break;
        if(covered[c])continue;
        const std::array<std::size_t,1> source{c};
        auto within=surface_distance(g,geo,source,plan.separation);
        HotRegion r;r.sector=g.sectors[c];r.vertex=ct.sector_vert[r.sector];r.boundary=!g.interior[c];r.excess=excess[c];
        r.log_area=field.log_area[r.sector];r.log_anisotropy=field.log_anisotropy[r.sector];
        for(std::size_t i=0;i<n;++i)if(std::isfinite(within[i])) {covered[i]=1;r.weighted_excess+=field.mass[g.sectors[i]]*excess[i];r.mass+=field.mass[g.sectors[i]];}
        plan.regions.push_back(r);region_distance.push_back(std::move(within));
    }
    if(plan.regions.empty()) {plan.skipped="no_hotspot";return plan;}
    // Tips per region.
    const auto k=static_cast<std::size_t>(settings.targets_per_region);
    for(std::size_t ri=0;ri<plan.regions.size();++ri) {
        const auto& within=region_distance[ri];const auto centre=g.local(plan.regions[ri].sector);
        if(!plan.regions[ri].boundary) {
            // Interior region: the worst interior sectors near the centre.
            std::vector<std::pair<double,std::size_t>> scored;
            for(std::size_t i=0;i<n;++i) {
                if(!std::isfinite(within[i]) || !g.interior[i] || !a.traversable[i] || !(excess[i]>0))continue;
                const auto v=ct.sector_vert[g.sectors[i]];
                const double factor=std::isfinite(seam_distance[i])?std::min(1.0,seam_distance[i]/plan.separation):1.0;
                scored.emplace_back(excess[i]*suitability(m,v,cost)*factor,i);
            }
            std::stable_sort(scored.begin(),scored.end(),[&](const auto& x,const auto& y) {return x.first!=y.first?x.first>y.first:g.sectors[x.second]<g.sectors[y.second];});
            for(std::size_t j=0;j<std::min(k,scored.size());++j) {
                const auto i=scored[j].second;ReliefTarget target;target.region=static_cast<MeshIndex>(ri);target.rank=static_cast<MeshIndex>(j);
                target.tip_sector=g.sectors[i];target.tip_vertex=ct.sector_vert[target.tip_sector];target.score=scored[j].first;target.excess=excess[i];
                target.suitability=suitability(m,target.tip_vertex,cost);target.seam_distance=seam_distance[i];target.depth=depth[i];
                plan.targets.push_back(target);
            }
        } else {
            // Boundary region: slits start at the hot boundary sectors (the centre first) and
            // reach toward the deepest interior sector, measured from that start.
            std::vector<std::size_t> starts;
            if(a.start[centre])starts.push_back(centre);
            std::vector<std::size_t> others;
            for(std::size_t i=0;i<n;++i)if(i!=centre && std::isfinite(within[i]) && !g.interior[i] && a.start[i] && excess[i]>0)others.push_back(i);
            std::stable_sort(others.begin(),others.end(),[&](std::size_t x,std::size_t y) {
                return excess[x]!=excess[y]?excess[x]>excess[y]:(within[x]!=within[y]?within[x]<within[y]:g.sectors[x]<g.sectors[y]);});
            for(auto i:others)if(starts.size()<k)starts.push_back(i);
            if(starts.size()>k)starts.resize(k);
            for(std::size_t j=0;j<starts.size();++j) {
                const std::array<std::size_t,1> source{starts[j]};const auto from=surface_distance(g,geo,source);
                std::size_t tip=none;double best=-infinity;
                for(std::size_t i=0;i<n;++i) {
                    if(!g.interior[i] || !a.traversable[i] || !std::isfinite(from[i]) || !std::isfinite(depth[i]))continue;
                    const double value=depth[i]-0.5*from[i];
                    if(value>best || (tip!=none && value==best && g.sectors[i]<g.sectors[tip])) {best=value;tip=i;}
                }
                if(tip==none)continue;
                ReliefTarget target;target.region=static_cast<MeshIndex>(ri);target.rank=static_cast<MeshIndex>(j);target.tip_sector=g.sectors[tip];
                target.tip_vertex=ct.sector_vert[target.tip_sector];target.start_sector=g.sectors[starts[j]];target.score=excess[starts[j]];
                target.excess=excess[starts[j]];target.suitability=suitability(m,target.tip_vertex,cost);target.seam_distance=seam_distance[tip];target.depth=depth[tip];
                plan.targets.push_back(target);
            }
        }
    }
    // Routes and reaches.
    struct Entry {std::size_t target;std::size_t reach;std::vector<MeshIndex> sectors,vertices,edges;double length=0,available=0;};
    std::vector<Entry> entries;
    for(std::size_t ti=0;ti<plan.targets.size();++ti) {
        const auto& target=plan.targets[ti];
        const auto start=target.start_sector?std::optional<std::size_t>(g.local(*target.start_sector)):std::nullopt;
        const auto route=slit_route(g,a,cost,g.local(target.tip_sector),start);
        if(!route || route->edges.empty()) {++plan.unreachable;continue;}
        std::vector<double> cumulative{0.0};for(auto e:route->edges)cumulative.push_back(cumulative.back()+geo.edge_length[e]);
        const double total=cumulative.back();
        for(std::size_t ri=0;ri<settings.reaches.size();++ri) {
            // Snap the requested physical reach to the nearest path vertex (at least one edge).
            const double want=settings.reaches[ri]*total;std::size_t best=1;
            for(std::size_t i=1;i<cumulative.size();++i)if(std::abs(cumulative[i]-want)<std::abs(cumulative[best]-want))best=i;
            Entry entry;entry.target=ti;entry.reach=ri;entry.length=cumulative[best];entry.available=total;
            for(std::size_t i=0;i<=best;++i) {entry.sectors.push_back(g.sectors[route->nodes[i]]);entry.vertices.push_back(ct.sector_vert[g.sectors[route->nodes[i]]]);}
            entry.edges.assign(route->edges.begin(),route->edges.begin()+static_cast<std::ptrdiff_t>(best));
            entries.push_back(std::move(entry));
        }
    }
    plan.generated=entries.size();
    // Selection by rank (the worst point of every region first, then the candidates near it;
    // the requested reaches in order), deduplicated by edge set, bounded per pass.
    std::stable_sort(entries.begin(),entries.end(),[&](const Entry& x,const Entry& y) {
        const auto& tx=plan.targets[x.target];const auto& ty=plan.targets[y.target];
        return std::tuple(tx.rank,tx.region,x.reach,x.target)<std::tuple(ty.rank,ty.region,y.reach,y.target);
    });
    std::map<std::vector<MeshIndex>,std::size_t> seen;std::vector<ReliefPath> distinct;
    // A tip belongs to the best-ranked target that reaches it: other starts toward the same
    // tip (the same deep point seen from further along a hot boundary) repeat that slit.
    std::map<MeshIndex,std::size_t> tip_owner;
    for(const auto& entry:entries) {
        const auto owner=tip_owner.emplace(plan.targets[entry.target].tip_sector,entry.target).first->second;
        if(owner!=entry.target) {++plan.shared_tips;continue;}
        auto key=entry.edges;std::sort(key.begin(),key.end());
        if(const auto found=seen.find(key);found!=seen.end()) {++plan.duplicates;distinct[found->second].reaches.push_back(settings.reaches[entry.reach]);continue;}
        seen.emplace(std::move(key),distinct.size());
        ReliefPath p;p.target=static_cast<MeshIndex>(entry.target);p.sectors=entry.sectors;p.vertices=entry.vertices;p.edges=entry.edges;
        p.length=entry.length;p.available=entry.available;p.reaches={settings.reaches[entry.reach]};
        distinct.push_back(std::move(p));
    }
    const auto cap=static_cast<std::size_t>(settings.max_candidates_per_pass);
    if(distinct.size()>cap) {plan.truncated=distinct.size()-cap;distinct.resize(cap);}
    // Evaluation order: target by target in rank order, the shortest slit of each first.
    std::stable_sort(distinct.begin(),distinct.end(),[&](const ReliefPath& x,const ReliefPath& y) {
        const auto& tx=plan.targets[x.target];const auto& ty=plan.targets[y.target];
        if(std::tuple(tx.rank,tx.region)!=std::tuple(ty.rank,ty.region))return std::tuple(tx.rank,tx.region)<std::tuple(ty.rank,ty.region);
        if(x.length!=y.length)return x.length<y.length;
        auto kx=x.edges,ky=y.edges;std::sort(kx.begin(),kx.end());std::sort(ky.begin(),ky.end());return kx<ky;
    });
    plan.paths=std::move(distinct);
    if(plan.paths.empty())plan.skipped=plan.targets.empty()?"no_target":"unreachable";
    return plan;
}
// ------------------------------------------------------------ generators
namespace {
struct ReliefState {
    Revision last_revision=0;
    std::set<std::string> tried;                   // targets tried on an unchanged chart or split region
    std::optional<std::vector<double>> cost;
};
std::vector<double> routing_cost(const ProposalContext& ctx) {
    // The enhanced routing heuristic (Section 12) with the run's paint, or without paint.
    const ResolvedPreferences* prefs=ctx.preferences.get();ResolvedPreferences local;
    if(!prefs || (prefs->raw_exposure.empty() && !ctx.raw_exposure.empty())) {
        local=prefs?*prefs:resolve_preferences(*ctx.mesh,{},ctx.options.enhanced.preferences,std::nullopt,{},ctx.options.layout.feature_angle,
                                               ctx.options.layout.hide_dir);
        local.raw_exposure=ctx.raw_exposure;prefs=&local;
    }
    return enhanced_edge_cost(*ctx.mesh,*prefs,ctx.options.enhanced.preferences);
}
std::vector<Proposal> relief_round(const RunState& incumbent,ProposalContext& ctx,const ReliefSettings& settings,ReliefState& state,std::int64_t pass) {
    auto& tally=ctx.counters.by_generator["hotspot_relief"];
    // A later pass on the same incumbent would propose exactly the same slits again.
    if(pass>1 && state.last_revision==incumbent.revision) {++tally["unchanged_incumbent"];return {};}
    state.last_revision=incumbent.revision;
    const auto& mesh=*ctx.mesh;const auto& t=mesh.tables();const auto& ct=incumbent.surface->tables();
    struct Target {std::string kind,key;std::shared_ptr<const CutSurface> surface;std::span<const Vec2> uv;MeshIndex chart=0;
                   std::vector<MeshIndex> removed;const SplitRecord* record=nullptr;std::size_t split=0;};
    std::vector<Target> targets;
    // Over-budget disk charts of the incumbent: relief slits only.
    for(MeshIndex ch=0;ch<ct.n_charts;++ch) {
        const auto& o=incumbent.charts[ch];
        if(o.exempt_only || !o.measured || o.within_budget)continue;
        if(!o.disk) {++tally["not_disk"];continue;}
        Target target;target.kind="relief";target.key="chart:"+o.key;target.surface=incumbent.surface;target.uv=incumbent.uv;target.chart=ch;
        targets.push_back(std::move(target));
    }
    // Broad over-budget splits of the baseline: weld the split, then relief on the reunited chart.
    for(std::size_t k=0;k<ctx.baseline_splits.size();++k) {
        const auto& record=ctx.baseline_splits[k];
        const auto present=split_edges_present(record,incumbent.cuts);
        if(present.empty()) {++tally["split_resolved"];continue;}
        bool protected_edge=false;
        for(auto e:present)protected_edge=protected_edge || ctx.constraints.locked[e] || ctx.constraints.required[e];
        if(protected_edge) {++tally["split_edges_protected"];continue;}
        auto base=incumbent.cuts;for(auto e:present)base[e]=0;
        auto surface=cut_surface(ctx.mesh,base);
        const auto reunited=surface->tables().chart_of_face[record.faces.front()];
        if(!(surface->chart_key(reunited)==record.surface->chart_key(record.chart))) {++tally["split_region_changed"];continue;}
        Target target;target.kind="relief_before_split";target.surface=std::move(surface);target.uv=record.uv;target.chart=reunited;
        target.removed=present;target.record=&record;target.split=k;
        std::string key="split:"+std::to_string(k)+":";for(auto e:present)key+=std::to_string(e)+",";
        target.key=std::move(key);
        targets.push_back(std::move(target));
    }
    if(targets.empty())return {};
    if(!state.cost)state.cost=routing_cost(ctx);
    struct Candidate {Proposal proposal;std::size_t target=0,order=0;};
    std::vector<Candidate> candidates;
    for(std::size_t ti=0;ti<targets.size();++ti) {
        const auto& target=targets[ti];
        throw_if_cancelled(ctx.cancel);
        // A deadline stops new candidate work; the incumbent stays (Section 21).
        if(ctx.deadline.expired()) {++tally["deadline"];break;}
        if(state.tried.contains(target.key)) {++tally["already_tried"];continue;}
        state.tried.insert(target.key);
        ReliefInput input;input.surface=target.surface;input.uv=target.uv;input.chart=target.chart;input.exempt=ctx.exempt;
        input.constraints=&ctx.constraints;input.edge_cost=*state.cost;
        // T39: the field's smoothing has a physical radius (its step count is defined at one
        // step per reference length of the pass's frozen scale; hops mode keeps the count).
        auto local=settings;
        const auto smoothing=context_scale(ctx,incumbent).diffusion(settings.smoothing_steps,target.surface->tables().chart_faces[target.chart]);
        local.smoothing_steps=smoothing.steps;
        const auto plan=plan_relief(input,local);
        tally["targets"]+=1;tally["hot_regions"]+=plan.regions.size();tally["slit_candidates"]+=plan.paths.size();
        tally["candidates_truncated"]+=plan.truncated;tally["candidate_duplicates"]+=plan.duplicates;tally["unreachable_tips"]+=plan.unreachable;
        tally["shared_tips"]+=plan.shared_tips;
        if(!plan.skipped.empty())++tally["skipped_"+plan.skipped];
        const auto chart_key=target.surface->chart_key(target.chart).digest();
        for(std::size_t pi=0;pi<plan.paths.size();++pi) {
            const auto& path=plan.paths[pi];
            Candidate c;c.target=ti;c.order=pi;
            auto& p=c.proposal;p.identity.parent_revision=incumbent.revision;p.identity.generator="hotspot_relief";
            p.identity.generator_version=hotspot_relief_version;
            // The net delta against the incumbent: a split edge the slit reuses stays cut.
            const std::set<MeshIndex> slit(path.edges.begin(),path.edges.end()),welded(target.removed.begin(),target.removed.end());
            for(auto e:path.edges)if(!welded.contains(e))p.identity.added_edges.push_back(e);
            for(auto e:target.removed)if(!slit.contains(e))p.identity.removed_edges.push_back(e);
            const auto& tip=plan.targets[path.target];const auto& region=plan.regions[tip.region];
            p.predicted={{"slit_length",path.length},{"region_excess",region.excess}};
            p.provenance={{"kind",target.kind},{"pass",pass},{"chart_key",chart_key},{"region",tip.region},{"region_boundary",region.boundary},
                          {"region_excess",region.excess},{"target_rank",tip.rank},{"tip_vertex",tip.tip_vertex},{"start_vertex",path.vertices.front()},
                          {"slit_edges",path.edges},{"slit_length",path.length},{"available_reach",path.available},{"reaches",path.reaches},
                          {"field",plan.field.value("version","")},{"settings",local.to_json()},{"smoothing",smoothing.to_json()}};
            // GB m10: split_index names the staged run-origin rebase this reunion is measured against.
            if(target.record) {p.provenance["split"]=target.record->summary();p.provenance["split_index"]=target.split;}
            candidates.push_back(std::move(c));
        }
    }
    // Stable IDs in the plans' evaluation order (select_proposal orders by ID): target chart
    // by target chart, each plan's best-ranked tip first and its shortest slit first.
    std::stable_sort(candidates.begin(),candidates.end(),[](const Candidate& x,const Candidate& y) {
        return std::tie(x.target,x.order)<std::tie(y.target,y.order);
    });
    std::vector<Proposal> out;
    for(std::size_t i=0;i<candidates.size();++i) {
        auto& p=candidates[i].proposal;p.identity.id=32'000'000+static_cast<std::uint64_t>(pass)*100'000+i+1;
        out.push_back(std::move(p));
    }
    (void)t;
    return out;
}
}
std::vector<ProposalGenerator> hotspot_relief_generators(const PipelineOptions&,ReliefSettings settings) {
    check_settings(settings);
    auto state=std::make_shared<ReliefState>();
    std::vector<ProposalGenerator> out;
    for(std::int64_t pass=1;pass<=settings.passes;++pass)
        out.push_back({"hotspot_relief",hotspot_relief_version,[state,settings,pass](const RunState& incumbent,ProposalContext& ctx) {
            return relief_round(incumbent,ctx,settings,*state,pass);
        }});
    return out;
}
}
