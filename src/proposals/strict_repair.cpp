// SPDX-License-Identifier: GPL-3.0-or-later
// Gate B strict-overlap policy (Section 11.3): the validated repair of candidate zero's
// strict-only defects and the explicit report of the defects it cannot repair. The strict
// overlap check stays a validity check of every new candidate (a trial's solved or moved charts
// must pass it); candidate zero's own charts that fail only that check never turn into a refusal
// native_baseline would not make. Used by policy enhanced_v1 only.
#include "cleanseams/proposals.hpp"
#include "cleanseams/pipeline.hpp"
#include "cleanseams/chart_solver.hpp"
#include "cleanseams/cache.hpp"
#include "cleanseams/errors.hpp"
#include "cleanseams/detail/numpy_compat.hpp"
#include <algorithm>
#include <cmath>
#include <cstring>
#include <functional>
#include <limits>
#include <map>
#include <queue>
#include <set>
#include <tuple>
namespace cleanseams {
namespace {
using Json=nlohmann::json;
constexpr const char* stage_name="strict_overlap_repair";
constexpr std::size_t listed_faces=64;       // overlapping faces listed per chart in reports
bool always_cut(const MeshTables& t,std::size_t e) {
    return t.edge_boundary[e] || !t.edge_manifold[e] || !t.edge_consistent[e] || t.edge_corners[e].size()!=2;
}
// The repair must pass the budgets too (keeps_budget, search.cpp): a compliant incumbent stays
// compliant, an over-budget one gets no higher worst ratio and no higher total excess.
std::vector<TransferIndex> exempt_faces(std::span<const std::uint8_t> mask) {
    std::vector<TransferIndex> out;for(std::size_t f=0;f<mask.size();++f)if(mask[f])out.push_back(static_cast<TransferIndex>(f));return out;
}
// The strict check of one chart with its overlapping pairs retained.
OverlapReport chart_pairs(const ProposalContext& ctx,const QualityGeometry& geometry,std::span<const Vec2> uv,MeshIndex chart) {
    OverlapOptions o;o.area_tolerance=ctx.options.enhanced.overlap_area_tolerance;o.ignored_faces=exempt_faces(ctx.exempt);o.retain_pairs=true;
    const std::array<MeshIndex,1> one{chart};
    return chart_overlaps(geometry,CornerUvView{uv},o,std::span<const MeshIndex>(one));
}
struct DefectFacts {std::uint64_t pairs=0;double area=0;bool truncated=false;std::vector<MeshIndex> faces;std::vector<OverlapPair> list;};
DefectFacts defect_facts(const ProposalContext& ctx,const RunState& s,MeshIndex chart) {
    DefectFacts out;const auto report=chart_pairs(ctx,*quality_geometry(s.surface),s.uv,chart);
    out.pairs=report.counts.size()>chart?report.counts[chart]:report.overlaps;out.truncated=report.pairs_truncated;out.list=report.pairs;
    std::set<MeshIndex> faces;
    for(const auto& p:report.pairs) {out.area+=p.area;faces.insert(p.first_face);faces.insert(p.second_face);}
    out.faces.assign(faces.begin(),faces.end());
    return out;
}
Json facts_json(const DefectFacts& d) {
    std::vector<MeshIndex> listed(d.faces.begin(),d.faces.begin()+static_cast<std::ptrdiff_t>(std::min(d.faces.size(),listed_faces)));
    return {{"pairs",d.pairs},{"overlap_area",d.area},{"pairs_truncated",d.truncated},{"overlapping_faces",d.faces.size()},{"faces",listed}};
}
std::string origin_label(const ProposalContext& ctx,const RunState& s,MeshIndex chart) {
    const auto& origin=ctx.envelope.origin_chart_of_face();const auto faces=s.surface->tables().chart_faces[chart];
    if(origin.empty() || faces.empty())return s.charts[chart].key;
    return "origin:"+std::to_string(origin[faces[0]]);
}
std::optional<MeshIndex> chart_with_key(const RunState& s,const std::string& key) {
    for(MeshIndex ch=0;ch<s.charts.size();++ch)if(s.charts[ch].key==key)return ch;
    return std::nullopt;
}
// ------------------------------------------------------------ native re-maps of one chart
std::vector<CornerPin> corner_pins(const ProposalContext& ctx) {
    std::vector<CornerPin> out;for(const auto& [c,p]:ctx.constraints.pins)out.push_back({static_cast<TransferIndex>(c),p});return out;
}
ChartSolveRequest chart_request(const ProposalContext& ctx,const std::shared_ptr<const QualityGeometry>& geometry,MeshIndex chart,std::uint64_t iterations) {
    ChartSolveRequest r;r.geometry=geometry;r.chart=static_cast<TransferIndex>(chart);r.pins=corner_pins(ctx);r.pin_completion=pin_policy(ctx.options);
    r.iterations=iterations;r.validate=false;
    return r;
}
std::vector<Vec2> with_chart(std::span<const Vec2> uv,const ChartSolveResult& solved) {
    std::vector<Vec2> out(uv.begin(),uv.end());
    for(std::size_t i=0;i<solved.corner_ids.size();++i)out[solved.corner_ids[i]]=solved.corner_uv[i];
    return out;
}
// One rung of the ladder: the method's map of `chart` in a copy of the state's map, or the reason
// it does not apply (no map, no solve).
struct Remap {std::optional<std::vector<Vec2>> uv;std::string not_applicable;};
Remap remap(const ProposalContext& ctx,const RunState& s,MeshIndex chart,const std::string& method,std::span<const MeshIndex> faces) {
    const auto geometry=quality_geometry(s.surface);
    if(method.rfind("untangle:",0)==0) {
        // GB N2: a pin-only defect has no overlapping face to untangle around.
        if(faces.empty())return {std::nullopt,"no_overlapping_faces"};
        // The repair ladder's local untangle (L2b) on the rings around the overlapping faces.
        RepairContext rc;rc.mesh=ctx.mesh;rc.pins=&ctx.constraints.pins;rc.pin_completion=pin_policy(ctx.options);rc.options=&ctx.options;rc.cancel=ctx.cancel;
        auto uv=local_untangle(rc,geometry,chart,s.uv,faces,std::stoll(method.substr(9)));
        if(!uv)return {std::nullopt,"no_free_sectors"};
        return {std::move(*uv),""};
    }
    if(method.rfind("arap:",0)==0)return {with_chart(s.uv,solve_arap(chart_request(ctx,geometry,chart,std::stoull(method.substr(5))))),""};
    if(method=="lscm")return {with_chart(s.uv,solve_lscm(chart_request(ctx,geometry,chart,10))),""};
    if(method=="harmonic" || method.rfind("harmonic_arap:",0)==0) {
        const auto iterations=method=="harmonic"?10:std::stoull(method.substr(14));
        std::optional<ChartSolveResult> convex;
        try {convex=convex_boundary_map(chart_request(ctx,geometry,chart,iterations));}
        catch(const ContractError& error) {if(error.code==ErrorCode::constraint_conflict)return {std::nullopt,"pinned_chart"};throw;}
        if(!convex)return {std::nullopt,"not_a_single_loop_disk"};
        if(method=="harmonic")return {with_chart(s.uv,*convex),""};
        auto request=chart_request(ctx,geometry,chart,iterations);request.initial_sector_uv=convex->sector_uv;
        return {with_chart(s.uv,solve_arap(request)),""};
    }
    throw ContractError(ErrorCode::invalid_input,"strict_repair.method","unknown repair method "+method);
}
// ------------------------------------------------------------ minimal splits of one chart
// Competing multi-source growth over the chart's dual (the nearest seed set labels a node; ties
// go to the lower label, then the lower node): the partition of `separate` and `bisect`.
std::vector<std::int64_t> competing_labels(const DualGraph& g,const std::array<std::vector<MeshIndex>,2>& seeds) {
    const auto n=g.size();std::vector<double> dist(n,std::numeric_limits<double>::infinity());std::vector<std::int64_t> label(n,-1);
    using Item=std::tuple<double,std::int64_t,MeshIndex>;std::priority_queue<Item,std::vector<Item>,std::greater<>> queue;
    for(std::int64_t l=0;l<2;++l)for(auto s:seeds[static_cast<std::size_t>(l)])if(label[s]<0) {dist[s]=0;label[s]=l;queue.push({0.0,l,s});}
    while(!queue.empty()) {
        const auto [d,l,v]=queue.top();queue.pop();
        if(d>dist[v] || l!=label[v])continue;
        for(auto j=g.indptr[v];j<g.indptr[v+1];++j) {
            const auto w=g.indices[j];const double nd=d+g.weights[j];
            if(nd<dist[w] || (nd==dist[w] && l<label[w])) {dist[w]=nd;label[w]=l;queue.push({nd,l,w});}
        }
    }
    return label;
}
std::vector<double> distances(const DualGraph& g,std::span<const MeshIndex> sources) {
    std::vector<double> dist(g.size(),std::numeric_limits<double>::infinity());
    using Item=std::pair<double,MeshIndex>;std::priority_queue<Item,std::vector<Item>,std::greater<>> queue;
    for(auto s:sources) {dist[s]=0;queue.push({0.0,s});}
    while(!queue.empty()) {
        const auto [d,v]=queue.top();queue.pop();if(d>dist[v])continue;
        for(auto j=g.indptr[v];j<g.indptr[v+1];++j) {const auto w=g.indices[j];if(d+g.weights[j]<dist[w]) {dist[w]=d+g.weights[j];queue.push({dist[w],w});}}
    }
    return dist;
}
MeshIndex farthest(const std::vector<double>& d) {
    MeshIndex best=0;double far=-1;for(std::size_t i=0;i<d.size();++i)if(std::isfinite(d[i]) && d[i]>far) {far=d[i];best=static_cast<MeshIndex>(i);}
    return best;
}
// The two sides of the overlapping pairs: the pair with the largest overlap first; a pair with
// one labelled face gets the other face on the opposite side; a new pair is oriented by its
// distances to the sides found so far (both faces of one node, or both already labelled: skipped).
std::array<std::vector<MeshIndex>,2> overlap_sides(const RegionDual& dual,std::vector<OverlapPair> pairs) {
    std::stable_sort(pairs.begin(),pairs.end(),[](const OverlapPair& a,const OverlapPair& b) {return a.area>b.area;});
    std::vector<std::int64_t> side(dual.graph.size(),-1);std::array<std::vector<MeshIndex>,2> out;
    for(const auto& p:pairs) {
        auto a=dual.node(p.first_face),b=dual.node(p.second_face);
        if(a==b || (side[a]>=0 && side[b]>=0))continue;
        std::int64_t sa=0;
        if(side[a]>=0)sa=side[a];
        else if(side[b]>=0)sa=1-side[b];
        else if(!out[0].empty() && !out[1].empty()) {
            const auto d0=distances(dual.graph,out[0]),d1=distances(dual.graph,out[1]);
            sa=d0[a]+d1[b]<=d1[a]+d0[b]?0:1;
        }
        if(side[a]<0) {side[a]=sa;out[static_cast<std::size_t>(sa)].push_back(a);}
        if(side[b]<0) {side[b]=1-sa;out[static_cast<std::size_t>(1-sa)].push_back(b);}
    }
    return out;
}
// The cheapest separation of the two sides (a minimum cut, Dinic's max flow): every interior,
// uncut edge between two dual nodes is priced by what cutting it adds to J's length and exposure
// terms, length * (w_length + w_exposure * raw exposure) (forbidden edges are inside the nodes:
// never cut). Labels: 0 for the nodes the source side still reaches, 1 for the others.
std::vector<std::int64_t> min_cut_labels(const ProposalContext& ctx,const RegionDual& dual,std::span<const std::uint8_t> barriers,
                                         const std::array<std::vector<MeshIndex>,2>& sides) {
    const auto& m=*ctx.mesh;const auto& t=m.tables();const auto& g=m.geometry();const auto ne=t.edge_vertices.size();
    const auto& w=ctx.options.enhanced.objective;const bool exposure=ctx.raw_exposure.size()==ne;
    const auto n=dual.graph.size();std::vector<std::int64_t> local(m.geometry().face_area.size(),-1);
    for(std::size_t i=0;i<dual.faces.size();++i)local[dual.faces[i]]=static_cast<std::int64_t>(dual.node_of_face[i]);
    std::map<std::pair<MeshIndex,MeshIndex>,double> capacity;double total=0;
    for(std::size_t e=0;e<ne;++e) {
        const auto pair=t.edge_faces[e];
        if(pair.size()!=2 || barriers[e] || local[pair[0]]<0 || local[pair[1]]<0)continue;
        const auto a=static_cast<MeshIndex>(local[pair[0]]),b=static_cast<MeshIndex>(local[pair[1]]);
        if(a==b)continue;
        const double c=g.edge_length[e]*(w.w_length+(exposure?w.w_exposure*ctx.raw_exposure[e]:0.0));
        capacity[{std::min(a,b),std::max(a,b)}]+=c;total+=c;
    }
    const double floor=std::max(total,1e-300)*1e-12,infinite=2*total+1;
    struct Arc {std::size_t to;double cap;std::size_t rev;};
    const std::size_t source=n,sink=n+1;std::vector<std::vector<Arc>> graph(n+2);
    auto add=[&](std::size_t u,std::size_t v,double cu,double cv) {
        graph[u].push_back({v,cu,graph[v].size()});graph[v].push_back({u,cv,graph[u].size()-1});
    };
    for(const auto& [key,c]:capacity)add(key.first,key.second,std::max(c,floor),std::max(c,floor));
    for(auto v:sides[0])add(source,v,infinite,0);
    for(auto v:sides[1])add(v,sink,infinite,0);
    std::vector<std::int64_t> level(n+2);std::vector<std::size_t> next(n+2);
    auto bfs=[&] {
        std::fill(level.begin(),level.end(),-1);level[source]=0;std::vector<std::size_t> queue{source};
        for(std::size_t i=0;i<queue.size();++i)for(const auto& a:graph[queue[i]])
            if(a.cap>floor && level[a.to]<0) {level[a.to]=level[queue[i]]+1;queue.push_back(a.to);}
        return level[sink]>=0;
    };
    // Blocking flows on the level graph, with an explicit path stack (no recursion: charts of
    // tens of thousands of faces run on worker threads with small stacks).
    while(bfs()) {
        std::fill(next.begin(),next.end(),0);
        for(bool more=true;more;) {
            std::vector<std::pair<std::size_t,std::size_t>> path;std::size_t u=source;
            while(u!=sink) {
                bool advanced=false;
                for(;next[u]<graph[u].size();++next[u]) {
                    const auto& a=graph[u][next[u]];
                    if(a.cap>floor && level[a.to]==level[u]+1) {path.push_back({u,next[u]});u=a.to;advanced=true;break;}
                }
                if(advanced)continue;
                if(path.empty()) {more=false;break;}      // the source is exhausted: the phase ends
                level[u]=-1;u=path.back().first;path.pop_back();++next[u];
            }
            if(!more)break;
            double f=infinite;for(const auto& [v,i]:path)f=std::min(f,graph[v][i].cap);
            for(const auto& [v,i]:path) {auto& a=graph[v][i];a.cap-=f;graph[a.to][a.rev].cap+=f;}
        }
    }
    // The source side of the minimum cut: nodes still reachable in the residual graph.
    std::vector<std::int64_t> label(n,1);std::vector<std::size_t> queue{source};std::vector<std::uint8_t> seen(n+2,0);seen[source]=1;
    for(std::size_t i=0;i<queue.size();++i)for(const auto& a:graph[queue[i]])
        if(a.cap>floor && !seen[a.to]) {seen[a.to]=1;queue.push_back(a.to);if(a.to<n)label[a.to]=0;}
    return label;
}
struct SplitCandidate {std::string method;std::optional<Proposal> proposal;Json record;};
SplitCandidate split_candidate(const ProposalContext& ctx,const RunState& s,MeshIndex chart,const std::string& method,const DefectFacts& facts,
                               std::uint64_t id,const std::vector<double>& routing) {
    SplitCandidate out;out.method=method;out.record={{"method",method},{"id",id}};
    const auto& m=*ctx.mesh;const auto& t=m.tables();const auto ne=t.edge_vertices.size(),nf=m.geometry().face_area.size();
    const auto& hc=ctx.constraints;const auto faces_span=s.surface->tables().chart_faces[chart];
    const std::vector<MeshIndex> faces(faces_span.begin(),faces_span.end());
    // The chart's own cuts (its slits) are barriers: the sides of a slit are far apart on the chart.
    EdgeMask barriers(ne,0),weld(ne,0);
    for(std::size_t e=0;e<ne;++e) {
        barriers[e]=s.cuts[e] || always_cut(t,e) || (e<hc.locked.size() && hc.locked[e]) || (e<hc.required.size() && hc.required[e]);
        weld[e]=e<hc.forbidden.size() && hc.forbidden[e];
    }
    const auto dual=region_dual(m,faces,barriers,weld,{},ctx.options.enhanced.regions.w_normal,0.0);
    std::array<std::vector<MeshIndex>,2> seeds;
    if(method=="min_cut" || method=="separate")seeds=overlap_sides(dual,facts.list);
    else {
        // The ordinary split seeding: the node farthest from the chart's first face, then the
        // node farthest from that one.
        const std::array<MeshIndex,1> start{dual.node(faces.front())};
        const auto first=farthest(distances(dual.graph,start));const std::array<MeshIndex,1> from{first};
        seeds={std::vector<MeshIndex>{first},std::vector<MeshIndex>{farthest(distances(dual.graph,from))}};
    }
    if(seeds[0].empty() || seeds[1].empty() || seeds[0]==seeds[1]) {out.record["status"]="refused";out.record["reason"]="degenerate_seeds";return out;}
    const auto labels=method=="min_cut"?min_cut_labels(ctx,dual,barriers,seeds):competing_labels(dual.graph,seeds);
    std::vector<std::int64_t> label(nf,-1);
    for(std::size_t i=0;i<dual.faces.size();++i)label[dual.faces[i]]=labels[dual.node_of_face[i]];
    // Existing cuts stay; the border of the two regions is added.
    auto cuts=s.cuts;std::uint64_t border=0;
    for(std::size_t e=0;e<ne;++e) {
        const auto pair=t.edge_faces[e];
        if(always_cut(t,e) || pair.size()!=2 || cuts[e] || label[pair[0]]<0 || label[pair[1]]<0 || label[pair[0]]==label[pair[1]])continue;
        cuts[e]=1;++border;
    }
    std::vector<MeshIndex> seed_faces[2];
    for(std::size_t k=0;k<2;++k)for(auto node:seeds[k])seed_faces[k].push_back(dual.node_faces[node][0]);
    out.record["seed_faces"]={seed_faces[0],seed_faces[1]};out.record["border_edges"]=border;
    if(!border) {out.record["status"]="refused";out.record["reason"]="no_border";return out;}
    std::vector<std::uint8_t> scope(nf,0);for(auto f:faces)scope[f]=1;
    const auto opening=open_non_disk_charts(ctx.mesh,cuts,routing,hc,ctx.exempt,scope,Control{ctx.cancel,ctx.deadline}); // GB m7
    out.record["openings"]=opening.added;
    if(!opening.ok) {out.record["status"]="refused";out.record["reason"]="not_openable";out.record["failures"]=opening.failures;return out;}
    if(const auto v=constraint_violations(m,hc,cuts);!v.empty()) {out.record["status"]="refused";out.record["reason"]="constraint";out.record["violations"]=v;return out;}
    out.record["status"]="proposed";
    Json provenance={{"stage",stage_name},{"version",strict_repair_version},{"method",method},{"chart_key",s.charts[chart].key},
                     {"seed_faces",out.record["seed_faces"]},{"border_edges",border},{"openings",opening.added}};
    out.proposal=cut_proposal(id,s,cuts,stage_name,strict_repair_version,std::move(provenance));
    return out;
}
Json attempt_json(const EvaluatedProposal& e,const std::string& method) {
    // GB N7: what the evaluation never measured is null (a trial refused before its measurement).
    const bool measured=e.trial.evaluated_revision!=0;const Json none(nullptr);
    Json out={{"method",method},{"id",e.proposal.identity.id},{"admissible",e.admissible()},{"rejection",e.rejection},{"solves",e.solves},
              {"changed_charts",e.changed_charts},{"strict_defects",measured?Json(e.strict_defects.size()):none},
              {"J",e.objective.version.empty()?none:Json(e.objective.J)},
              {"budget",e.budget.measured?Json{{"compliant",e.budget.compliant},{"worst_ratio",e.budget.worst_ratio},{"total_excess",e.budget.total_excess}}:none},
              {"envelope",e.envelope.checked || !e.envelope.inside || !e.envelope.reason.empty()?e.envelope.to_json():none}};
    if(!e.rejection.empty())out["detail"]=e.detail;
    return out;
}
}
SearchRound repair_strict_overlaps(const RunState& incumbent,ProposalContext& ctx) {
    SearchRound round;round.state=incumbent;
    auto zero=run_state_summary(incumbent);zero["id"]=0;zero["candidate"]="incumbent";zero["admissible"]=incumbent.admissible();
    round.candidates.push_back(std::move(zero));
    round.winner={{"id",0},{"revision",incumbent.revision}};
    Json charts=Json::array();
    const auto& o=ctx.options.enhanced;const auto band=o.objective.tie_band;
    // The defects are identified by chart key: a commit on one chart keeps every other chart's key.
    std::vector<std::string> keys;for(auto ch:incumbent.strict_defects)keys.push_back(incumbent.charts[ch].key);
    static const std::vector<std::string> ladder{"untangle:1","untangle:2","untangle:3","arap:20","arap:50","lscm","harmonic_arap:20","harmonic"};
    const auto routing=context_routing_cost(ctx);
    std::uint64_t next_id=1;
    for(std::size_t k=0;k<keys.size();++k) {
        auto& s=round.state;
        const auto found=chart_with_key(s,keys[k]);
        if(!found || !s.charts[*found].strict_defect)continue;
        const auto chart=*found;const auto origin=origin_label(ctx,s,chart);
        const auto facts=defect_facts(ctx,s,chart);
        Json record={{"key",keys[k]},{"chart",chart},{"faces",s.surface->tables().chart_faces[chart].size()},{"origin_chart",origin},{"defect",facts_json(facts)},
                     {"attempts",Json::array()},{"status","unrepaired"}};
        bool repaired=false;
        // Bounds before any work: the pass's solve and candidate bounds on the origin chart and
        // this stage's share of them.
        auto bound_refusal=[&]() -> std::string {
            if(ctx.counters.solves_by_chart[origin]>=static_cast<std::uint64_t>(o.max_solves_per_chart))return "solve_budget_exhausted";
            if(ctx.counters.candidates_by_chart[origin]>=static_cast<std::uint64_t>(o.max_candidates_per_chart))return "candidate_budget_exhausted";
            if(ctx.search) {
                if(ctx.search->feature_solves[stage_name][origin]>=static_cast<std::uint64_t>(o.search.feature_solves_per_chart))return "feature_solve_budget_exhausted";
                if(ctx.search->feature_candidates[stage_name][origin]>=static_cast<std::uint64_t>(o.search.feature_candidates_per_chart))return "feature_candidate_budget_exhausted";
            }
            return "";
        };
        auto stop=[&]() -> bool {
            if(round.stopped.empty() && ctx.cancel.requested())round.stopped="cancelled";
            if(round.stopped.empty() && ctx.deadline.expired())round.stopped="deadline";
            return !round.stopped.empty();
        };
        // Accounting of one evaluated attempt; true when it repairs the chart (committed).
        auto settle=[&](EvaluatedProposal& e,const std::string& method,bool count_solve) {
            // A native re-map is one expensive solve of the chart, whatever its evaluation says.
            if(count_solve) {e.solves=1;++ctx.counters.solves;++ctx.counters.charts_solved;++ctx.counters.solves_by_chart[origin];}
            if(ctx.search && e.solves)++ctx.search->feature_solves[stage_name][origin];
            auto summary=attempt_json(e,method);
            bool commit=false;
            if(e.admissible()) {
                ++ctx.counters.candidates_by_chart[origin];if(ctx.search)++ctx.search->feature_candidates[stage_name][origin];
                const CandidateView trial{e.proposal.identity.id,true,&e.objective,&e.budget,e.strict_defects.size()};
                const CandidateView current{0,s.admissible(),&s.objective,&s.budget,s.strict_defects.size()};
                const bool budget_ok=keeps_budget(e.budget,s.budget,band);
                const bool ranks=compare_candidates(trial,current,o.objective)<0;
                // The chart itself must be gone from the defects (a trial that left it unchanged
                // ranks equal and is not a repair).
                const bool gone=std::none_of(e.strict_defects.begin(),e.strict_defects.end(),[&](MeshIndex c) {return e.charts[c].key==keys[k];});
                summary["keeps_budget"]=budget_ok;summary["ranks_before_incumbent"]=ranks;
                commit=budget_ok && ranks && gone;
                if(!budget_ok)summary["rejection"]="repair_over_budget";
                else if(!gone || !ranks)summary["rejection"]="not_a_repair";
            }
            summary["committed"]=commit;
            record["attempts"].push_back(summary);round.candidates.push_back(summary);
            if(!commit)return false;
            const auto before=s.revision;
            s=commit_proposal(s,e,ctx);
            round.committed=true;round.winner={{"id",e.proposal.identity.id},{"revision",s.revision},{"chart_key",keys[k]},{"method",method},{"parent_revision",before}};
            record["status"]="repaired";record["repaired_by"]={{"stage",stage_name},{"method",method},{"id",e.proposal.identity.id},{"revision",s.revision}};
            return true;
        };
        // 1. Re-maps of the chart's own cuts, least change first; the first validated repair wins.
        for(const auto& method:ladder) {
            if(repaired || stop())break;
            if(const auto refusal=bound_refusal();!refusal.empty()) {record["attempts"].push_back({{"method",method},{"rejection",refusal}});break;}
            const auto id=next_id++;
            Remap map;
            try {map=remap(ctx,s,chart,method,facts.faces);}
            catch(const ContractError& error) {
                if(error.code==ErrorCode::cancelled)throw;
                record["attempts"].push_back({{"method",method},{"id",id},{"rejection","solver_failure"},{"detail",{{"message",error.what()}}}});continue;
            } catch(const std::exception& error) {
                record["attempts"].push_back({{"method",method},{"id",id},{"rejection","solver_failure"},{"detail",{{"message",error.what()}}}});continue;
            }
            if(!map.uv) {record["attempts"].push_back({{"method",method},{"id",id},{"rejection","not_applicable"},{"reason",map.not_applicable}});continue;}
            Proposal p;p.identity.id=id;p.identity.generator=stage_name;p.identity.generator_version=strict_repair_version;
            // A writer trial: the charts whose corners the map moves are the changed charts (a map
            // that leaves the chart as it was is not a repair).
            p.identity.affected_chart_keys={keys[k]};p.initialization=std::string(strict_repair_version)+":"+method;
            p.provenance={{"stage",stage_name},{"version",strict_repair_version},{"method",method},{"chart_key",keys[k]}};
            auto e=evaluate_written_map(s,std::move(p),std::move(*map.uv),std::string(stage_name)+":"+method,ctx);
            repaired=settle(e,method,true);
        }
        // 2. The last resort: a minimal split of the chart (the best admissible candidate). GB N2: a
        // split separates overlapping faces; a pin-only defect (no overlapping pair) is never split.
        if(!repaired && facts.pairs==0 && round.stopped.empty())
            for(const auto* method:{"min_cut","separate","bisect"})
                record["attempts"].push_back({{"method",method},{"rejection","not_applicable"},{"reason","pin_deviation_only"}});
        if(!repaired && facts.pairs>0 && !stop()) {
            std::optional<EvaluatedProposal> best;std::string best_method;
            for(const auto* method:{"min_cut","separate","bisect"}) {
                if(stop())break;
                if(const auto refusal=bound_refusal();!refusal.empty()) {record["attempts"].push_back({{"method",method},{"rejection",refusal}});break;}
                const auto id=next_id++;
                auto candidate=split_candidate(ctx,s,chart,method,facts,id,routing);
                if(!candidate.proposal) {candidate.record["rejection"]=candidate.record.value("reason","refused");record["attempts"].push_back(candidate.record);continue;}
                EvaluatedProposal e;
                try {e=evaluate_proposal(s,std::move(*candidate.proposal),ctx);}
                catch(const ContractError& error) {
                    if(error.code==ErrorCode::cancelled) {round.stopped="cancelled";break;}
                    // GB (m4): a failing split trial is one refused attempt; the ladder goes on.
                    if(!isolated_failure(error))throw;
                    record["attempts"].push_back({{"method",method},{"id",id},{"rejection","trial_failure"},{"detail",failure_record(error)}});continue;
                } catch(const std::exception& error) {
                    record["attempts"].push_back({{"method",method},{"id",id},{"rejection","trial_failure"},{"detail",failure_record(error)}});continue;
                }
                if(ctx.search && e.solves)++ctx.search->feature_solves[stage_name][origin];
                auto summary=attempt_json(e,method);summary["split"]=candidate.record;
                if(e.admissible()) {
                    ++ctx.counters.candidates_by_chart[origin];if(ctx.search)++ctx.search->feature_candidates[stage_name][origin];
                    const CandidateView trial{e.proposal.identity.id,true,&e.objective,&e.budget,e.strict_defects.size()};
                    const CandidateView current{0,s.admissible(),&s.objective,&s.budget,s.strict_defects.size()};
                    // GB (seam_objective.v2): a split changes the cuts, so its J is compared with the margin.
                    const bool budget_ok=keeps_budget(e.budget,s.budget,band),
                               ranks=compare_candidates(trial,current,o.objective,objective_margin(ctx,e.trial.cuts,s.cuts,s.objective.area,s.objective.J))<0;
                    summary["keeps_budget"]=budget_ok;summary["ranks_before_incumbent"]=ranks;
                    if(!budget_ok)summary["rejection"]="repair_over_budget";
                    else if(!ranks)summary["rejection"]="not_a_repair";
                    else {
                        const CandidateView leader=best?CandidateView{best->proposal.identity.id,true,&best->objective,&best->budget,best->strict_defects.size()}:CandidateView{};
                        if(!best || compare_candidates(trial,leader,o.objective)<0) {best=e;best_method=method;}
                    }
                }
                summary["committed"]=false;
                record["attempts"].push_back(summary);round.candidates.push_back(summary);
            }
            if(best && round.stopped!="cancelled") {
                const auto before=s.revision;
                s=commit_proposal(s,*best,ctx);repaired=true;
                for(auto& a:record["attempts"])if(a.value("id",std::uint64_t{0})==best->proposal.identity.id)a["committed"]=true;
                for(auto& a:round.candidates)if(a.value("id",std::uint64_t{0})==best->proposal.identity.id && a.contains("method"))a["committed"]=true;
                round.committed=true;
                round.winner={{"id",best->proposal.identity.id},{"revision",s.revision},{"chart_key",keys[k]},{"method",best_method},{"parent_revision",before}};
                record["status"]="repaired";record["repaired_by"]={{"stage",stage_name},{"method",best_method},{"id",best->proposal.identity.id},{"revision",s.revision}};
            }
        }
        if(!round.stopped.empty() && !repaired)record["stopped"]=round.stopped;
        charts.push_back(std::move(record));
        if(!round.stopped.empty())break;
    }
    ctx.diagnostics[stage_name]={{"version",strict_repair_version},{"ladder",ladder},{"splits",{"min_cut","separate","bisect"}},{"charts",charts},
                                 {"stopped",round.stopped}};
    return round;
}
nlohmann::json strict_overlap_report(const ProposalContext& ctx,const RunState& origin,const RunState& retained,const nlohmann::json& stage) {
    Json out={{"policy",strict_overlap_policy},{"version",strict_repair_version},{"repair_enabled",ctx.options.enhanced.strict_repair.enabled}};
    const auto diagnostics=ctx.diagnostics.contains(stage_name)?ctx.diagnostics.at(stage_name):Json::object();
    std::map<std::string,Json> attempts;
    for(const auto& c:diagnostics.value("charts",Json::array()))attempts[c.at("key").get<std::string>()]=c;
    std::set<std::string> still;for(auto ch:retained.strict_defects)if(retained.charts[ch].overlaps)still.insert(retained.charts[ch].key);
    Json charts=Json::array();std::uint64_t repaired=0,unrepaired=0,pairs=0,defect_charts=0;double area=0;
    for(auto ch:origin.strict_defects) {
        if(!origin.charts[ch].overlaps)continue;   // GB N2: a pin-only defect is report.enhanced.pin_check's
        ++defect_charts;
        const auto& key=origin.charts[ch].key;const auto facts=defect_facts(ctx,origin,ch);
        pairs+=facts.pairs;area+=facts.area;
        const bool open=still.contains(key);(open?unrepaired:repaired)+=1;
        Json record={{"chart",ch},{"key",key},{"faces",origin.charts[ch].faces},{"defect",facts_json(facts)},{"status",open?"unrepaired":"repaired"}};
        if(const auto found=attempts.find(key);found!=attempts.end()) {
            record["attempts"]=found->second.value("attempts",Json::array());
            if(found->second.contains("repaired_by"))record["repaired_by"]=found->second.at("repaired_by");
            if(found->second.contains("stopped"))record["stopped"]=found->second.at("stopped");
        } else record["attempts"]=Json::array();
        // Replaced after the repair stage by another stage's validated trial.
        if(!open && !record.contains("repaired_by"))record["repaired_by"]={{"stage","validated_trial"}};
        charts.push_back(std::move(record));
    }
    Json left=Json::array();std::uint64_t left_pairs=0;double left_area=0;
    for(auto ch:retained.strict_defects) {
        if(!retained.charts[ch].overlaps)continue;
        const auto facts=defect_facts(ctx,retained,ch);left_pairs+=facts.pairs;left_area+=facts.area;
        left.push_back({{"chart",ch},{"key",retained.charts[ch].key},{"pairs",facts.pairs},{"overlap_area",facts.area}});
    }
    out["status"]=defect_charts==0?"none":unrepaired?"unrepaired":"repaired";
    out["candidate_zero"]={{"charts",defect_charts},{"pairs",pairs},{"overlap_area",area},{"strict_only",origin.strict_only}};
    out["charts"]=std::move(charts);out["repaired"]=repaired;out["unrepaired"]=unrepaired;
    out["retained"]={{"charts",left},{"pairs",left_pairs},{"overlap_area",left_area}};
    out["stage"]=stage.is_null()?Json(nullptr):Json{{"ran",!stage.contains("skipped")},{"committed",stage.value("committed",false)},
                                                    {"stopped",stage.value("stopped",stage.value("skipped",std::string()))}};
    return out;
}
namespace {
// The order of a float32 value on the number line (adjacent floats differ by one).
std::int64_t float32_order(double x) {
    const float f=static_cast<float>(x);std::uint32_t u=0;std::memcpy(&u,&f,sizeof u);
    return (u & 0x80000000u)?-static_cast<std::int64_t>(u & 0x7fffffffu):static_cast<std::int64_t>(u);
}
}
nlohmann::json pin_check_report(const ProposalContext& ctx,const RunState& origin,const RunState& retained) {
    const double tolerance=ctx.options.enhanced.pin_tolerance;
    Json out={{"version","pin_check.v1"},{"pins",ctx.constraints.pins.size()},
              {"rule",{{"strict",ctx.float32_storage?"|uv - pin| <= pin_tolerance * (1 + |pin|), or uv is exactly float32(pin) (float32-storing host)":
                                                     "|uv - pin| <= pin_tolerance * (1 + |pin|)"},
                       {"pin_tolerance",tolerance},{"float32_storage",ctx.float32_storage},
                       {"native","numpy.allclose(uv, pin): |uv - pin| <= 1e-8 + 1e-5 * |pin| per coordinate (native_baseline's pins_kept)"},
                       {"host_apply","a host writes every pinned UV exactly on its pin: within the native tolerance it is put back on the pin, beyond it nothing is applied (Session::result_for_apply, HOST-01)"}}}};
    // Every pin the strict rule does not keep (bounded list), with its deviation.
    auto deviations=[&](const RunState& s) {
        Json list=Json::array();std::uint64_t moved=0;const auto& ct=s.surface->tables();
        for(const auto& [corner,pin]:ctx.constraints.pins) {
            if(corner>=s.uv.size() || corner>=ct.chart_of_corner.size())continue;
            const auto& got=s.uv[corner];
            const bool strict=strict_pin_kept(got[0],pin[0],tolerance,ctx.float32_storage) && strict_pin_kept(got[1],pin[1],tolerance,ctx.float32_storage);
            if(strict)continue;
            ++moved;
            if(list.size()>=64)continue;
            const std::array<double,2> a{got[0],got[1]},b{pin[0],pin[1]};
            list.push_back({{"corner",corner},{"chart",ct.chart_of_corner[corner]},{"pin",{pin[0],pin[1]}},{"uv",{got[0],got[1]}},
                            {"delta",{got[0]-pin[0],got[1]-pin[1]}},
                            {"float32_steps",{float32_order(got[0])-float32_order(pin[0]),float32_order(got[1])-float32_order(pin[1])}},
                            {"native_kept",detail::numpy_allclose(a,b)}});
        }
        return Json{{"moved",moved},{"deviations",std::move(list)}};
    };
    const auto diagnostics=ctx.diagnostics.contains(stage_name)?ctx.diagnostics.at(stage_name):Json::object();
    std::map<std::string,Json> attempts;
    for(const auto& c:diagnostics.value("charts",Json::array()))attempts[c.at("key").get<std::string>()]=c;
    std::set<std::string> still;for(auto ch:retained.strict_defects)if(!retained.charts[ch].pins_kept)still.insert(retained.charts[ch].key);
    Json charts=Json::array();std::uint64_t repaired=0,unrepaired=0;
    for(auto ch:origin.strict_defects) {
        if(origin.charts[ch].pins_kept)continue;   // an overlap-only defect is report.enhanced.strict_overlap's
        const auto& key=origin.charts[ch].key;const bool open=still.contains(key);(open?unrepaired:repaired)+=1;
        Json record={{"chart",ch},{"key",key},{"faces",origin.charts[ch].faces},{"overlaps",origin.charts[ch].overlaps},{"status",open?"unrepaired":"repaired"}};
        if(const auto found=attempts.find(key);found!=attempts.end()) {
            record["attempts"]=found->second.value("attempts",Json::array());
            if(found->second.contains("repaired_by"))record["repaired_by"]=found->second.at("repaired_by");
            if(found->second.contains("stopped"))record["stopped"]=found->second.at("stopped");
        } else record["attempts"]=Json::array();
        if(!open && !record.contains("repaired_by"))record["repaired_by"]={{"stage","validated_trial"}};
        charts.push_back(std::move(record));
    }
    out["candidate_zero"]=deviations(origin);out["retained"]=deviations(retained);
    out["status"]=charts.empty()?"none":unrepaired?"unrepaired":"repaired";
    out["charts"]=std::move(charts);out["repaired"]=repaired;out["unrepaired"]=unrepaired;
    return out;
}
}
