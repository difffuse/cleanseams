// SPDX-License-Identifier: GPL-3.0-or-later
#include "cleanseams/routing.hpp"
#include "cleanseams/cache.hpp"
#include "cleanseams/errors.hpp"
#include "cleanseams/detail/pyset.hpp"
#include <algorithm>
#include <cmath>
#include <exception>
#include <set>
namespace cleanseams {
namespace {
void require(bool yes,const char* field,const char* message) {if(!yes)throw ContractError(ErrorCode::invalid_input,field,message);}
std::vector<TransferIndex> transfer(std::span<const MeshIndex> verts) {return {verts.begin(),verts.end()};}
nlohmann::json turns_json(const ChainTurns& t) {return {{"turns",t.turns},{"lane_changes",t.lane_changes},{"edges",t.edges}};}
}
nlohmann::json RefineResult::to_json() const {
    nlohmann::json j{{"chain",chain},{"status",status},{"reason",reason},{"old",old_route},
        {"new",new_route?nlohmann::json(*new_route):nlohmann::json(nullptr)}};
    if(cycle) {
        // Source refine_cycle: chain/status/reason/old/new/parts, plus costs and turns
        // once optimized; no weight variant, rings or attempts.
        auto p=nlohmann::json::array();for(const auto& part:parts) {auto x=part.to_json();x.erase("tried");p.push_back(x);}j["parts"]=p;
        if(status=="optimized") {
            j["cost_old"]=cost_old?nlohmann::json(*cost_old):nlohmann::json(nullptr);j["cost_new"]=cost_new?nlohmann::json(*cost_new):nlohmann::json(nullptr);
            j["turns_old"]=turns_json(turns_old);j["turns_new"]=turns_new?turns_json(*turns_new):nlohmann::json(nullptr);
        }
        return j;
    }
    j["variant"]=variant?nlohmann::json(*variant):nlohmann::json(nullptr);
    j["cost_old"]=cost_old?nlohmann::json(*cost_old):nlohmann::json(nullptr);
    j["cost_new"]=cost_new?nlohmann::json(*cost_new):nlohmann::json(nullptr);
    j["turns_old"]=turns_json(turns_old);j["turns_new"]=turns_new?turns_json(*turns_new):nlohmann::json(nullptr);
    if(rings)j["rings"]=*rings;
    if(max_rings)j["max_rings"]=*max_rings;
    auto tried_json=nlohmann::json::array();
    for(const auto& a:tried)tried_json.push_back({a.rings,a.weights,a.outcome});
    j["tried"]=tried_json;
    if(!parts.empty()) {auto p=nlohmann::json::array();for(const auto& part:parts) {auto x=part.to_json();x.erase("tried");p.push_back(x);}j["parts"]=p;}
    return j;
}
TopologyCheck topology_preserved(const std::shared_ptr<const CutSurface>& before,std::span<const std::uint8_t> seams_after,const Chain& chain,
                                 std::optional<std::span<const std::uint8_t>> band) {
    require(bool(before),"routing.cut","missing cut surface");
    const auto& b=before->tables();const auto& m=before->mesh();
    TopologyCheck out;out.after=cut_surface(before->mesh_owner(),seams_after);const auto& a=out.after->tables();
    auto fail=[&](const char* reason) {out.ok=false;out.reason=reason;return out;};
    if(a.n_charts!=b.n_charts)return fail("chart_count_changed");
    if(band) {
        require(band->size()==m.input().face_offsets.size()-1,"routing.band","band mask does not match the mesh faces");
        std::set<MeshIndex> owner;
        for(MeshIndex ch=0;ch<b.n_charts;++ch) {
            const auto faces=b.chart_faces[ch];if(faces.empty())continue;
            std::vector<std::size_t> counts(a.n_charts,0);for(auto f:faces)++counts[a.chart_of_face[f]];
            const auto nid=static_cast<MeshIndex>(std::max_element(counts.begin(),counts.end())-counts.begin());
            if(static_cast<double>(counts[nid])<0.6*static_cast<double>(faces.size()))return fail("chart_lost_majority");
            for(auto f:faces)if(!(*band)[f] && a.chart_of_face[f]!=nid)return fail("chart_split_outside_corridor");
            if(owner.contains(nid))return fail("charts_merged_outside_corridor");
            owner.insert(nid);
        }
    }
    if(out.after->adjacency_keys()!=before->adjacency_keys()) {
        // Charts are relabelled by first face; compare as pairs of lowest faces.
        auto pairs=[](const CutTables& t) {
            std::set<std::array<MeshIndex,2>> s;
            for(const auto& adj:t.adjacency) {
                const auto x=t.chart_faces[adj.charts[0]][0],y=t.chart_faces[adj.charts[1]][0];s.insert({std::min(x,y),std::max(x,y)});
            }
            return s;
        };
        if(pairs(b)!=pairs(a))return fail("adjacency_changed");
    }
    if(chain.charts)for(auto ch:*chain.charts) {
        const auto f=b.chart_faces[ch][0];const auto ch2=a.chart_of_face[f];
        if(b.chart_is_disk[ch] && !a.chart_is_disk[ch2])return fail("chart_not_disk");
    }
    if(!a.pinched_sectors.empty() && b.pinched_sectors.empty())return fail("pinch_created");
    out.ok=true;return out;
}
namespace {
// The after-surface charts the chain's charts become (new chart -> old chart; later set
// entries win, as the source's dict) and their ids in ascending order.
struct FlankingPlan {std::map<MeshIndex,MeshIndex> after_ids;std::vector<MeshIndex> ids;};
FlankingPlan flanking_plan(const CutSurface& before,const CutSurface& after,const Chain& chain) {
    const auto& b=before.tables();const auto& a=after.tables();
    detail::PyIntSet charts;for(auto ch:*chain.charts)charts.add(ch);
    FlankingPlan plan;
    for(auto ch:charts.order())plan.after_ids[a.chart_of_face[b.chart_faces[static_cast<std::size_t>(ch)][0]]]=static_cast<MeshIndex>(ch);
    for(const auto& [nid,ch]:plan.after_ids)plan.ids.push_back(nid);
    return plan;
}
// The chain's before-charts whose report `previous` lacks, sorted.
std::vector<MeshIndex> missing_before(const CutSurface& before,const Chain& chain,const FlankingCache& previous) {
    detail::PyIntSet charts;for(auto ch:*chain.charts)charts.add(ch);
    std::vector<MeshIndex> need;
    for(auto ch:charts.order())if(!previous.contains(before.chart_key(static_cast<MeshIndex>(ch)).digest()))need.push_back(static_cast<MeshIndex>(ch));
    std::sort(need.begin(),need.end());
    return need;
}
bool all_measured(const FlankingPlan& plan,const ChartMeasurements& after) {
    return std::all_of(plan.ids.begin(),plan.ids.end(),[&](MeshIndex nid) {return after.contains(nid);});
}
// The source's judgement of the measured after-charts against the before-reports.
FlankingCheck flanking_judge(const CutSurface& before,const FlankingPlan& plan,ChartMeasurements after,const FlankingCache& previous,bool require_budget) {
    FlankingCheck out;out.after=std::move(after);
    auto fail=[&](std::string reason) {out.ok=false;out.reason=std::move(reason);return out;};
    for(const auto& [nid,measured]:out.after) {
        const auto found=plan.after_ids.find(nid);if(found==plan.after_ids.end())continue;
        if(require_budget && !measured.within_budget)return fail("route_exceeds_budget");
        const auto key=before.chart_key(found->second).digest();const auto prev=previous.find(key);
        const auto& r=measured.report;
        if(r.flipped || measured.crossings || r.nonfinite || r.collapsed_uv)return fail("route_flips_or_crossings");
        for(const auto* metric:{"anisotropy_p95","area_p95"}) {
            const bool aniso=std::string(metric)=="anisotropy_p95";
            const double x=aniso?r.anisotropy_p95:r.area_p95;
            if(prev==previous.end())return fail("route_missing_distortion");
            const double y=aniso?prev->second.report.anisotropy_p95:prev->second.report.area_p95;
            if(!std::isfinite(x) || !std::isfinite(y))return fail("route_missing_distortion");
            if(x>y*1.05+1e-9)return fail(std::string("route_raises_")+metric);
        }
    }
    out.ok=true;return out;
}
// GB performance: the flanking measurements of a window of candidate routes of one chain,
// and of the chain's before-charts missing from `previous`, taken as one parallel call of
// single-chart evaluations (the native evaluator measures every chart independently, so a
// chart measured alone is the chart measured with others). judge(k), called for k = 0, 1,
// ... in schedule order, returns exactly what flanking_ok returns for candidate k at that
// point of the serial loop: it rethrows the failure flanking_ok would raise first (the
// candidate's lowest failing after-chart, then the before-charts') and commits the
// before-reports to `previous` where flanking_ok writes them. What the serial loop never
// reaches (later candidates, their errors) is discarded.
class FlankingWindow {
public:
    FlankingWindow(const std::shared_ptr<const CutSurface>& before,const std::vector<std::shared_ptr<const CutSurface>>& afters,const Chain& chain,
                   const ChartEvaluator& evaluate,const FlankingCache& previous):before_(before) {
        const auto n=afters.size();plans_.reserve(n);after_.resize(n);after_error_.resize(n);
        struct Task {const std::shared_ptr<const CutSurface>* surface;MeshIndex chart;std::size_t owner;}; // owner n: before-charts
        std::vector<Task> tasks;
        for(std::size_t k=0;k<n;++k) {
            plans_.push_back(flanking_plan(*before,*afters[k],chain));
            for(auto nid:plans_[k].ids)tasks.push_back({&afters[k],nid,k});
            if(k==0)for(auto ch:missing_before(*before,chain,previous))tasks.push_back({&before_,ch,n});
        }
        std::vector<std::optional<ChartMeasurements>> slots(tasks.size());std::vector<std::exception_ptr> errors(tasks.size());
        parallel_for(tasks.size(),[&](std::size_t i) {
            const auto& t=tasks[i];
            try {slots[i]=evaluate(*t.surface,std::span<const MeshIndex>(&t.chart,1),{});}
            catch(...) {errors[i]=std::current_exception();}
        });
        // Per owner, the failure of its lowest chart (one evaluator call rethrows its lowest
        // failing index); the single-chart results merge into one map per call.
        for(std::size_t i=0;i<tasks.size();++i) {
            const auto owner=tasks[i].owner;auto& error=owner==n?need_error_:after_error_[owner];
            if(errors[i]) {if(!error)error=errors[i];continue;}
            if(error)continue;
            auto& target=owner==n?need_measured_:after_[owner];
            for(auto& [ch,m]:*slots[i])target[ch]=std::move(m);
        }
    }
    FlankingCheck judge(std::size_t k,FlankingCache& previous,bool require_budget) {
        if(after_error_[k])std::rethrow_exception(after_error_[k]);
        if(!all_measured(plans_[k],after_[k])) {FlankingCheck out;out.after=std::move(after_[k]);out.reason="route_not_evaluated";return out;}
        if(!need_committed_) {
            // flanking_ok's before-charts: missing from `previous` when the window opened,
            // and nothing else writes it before the first judgement that reaches here.
            if(need_error_)std::rethrow_exception(need_error_);
            for(auto& [ch,m]:need_measured_)previous[before_->chart_key(ch).digest()]=std::move(m);
            need_committed_=true;
        }
        return flanking_judge(*before_,plans_[k],std::move(after_[k]),previous,require_budget);
    }
private:
    std::shared_ptr<const CutSurface> before_;
    std::vector<FlankingPlan> plans_;
    std::vector<ChartMeasurements> after_;
    std::vector<std::exception_ptr> after_error_;
    ChartMeasurements need_measured_;
    std::exception_ptr need_error_;
    bool need_committed_=false;
};
}
FlankingCheck flanking_ok(const std::shared_ptr<const CutSurface>& before,const std::shared_ptr<const CutSurface>& after,const Chain& chain,
                          const ChartEvaluator* evaluate,FlankingCache& previous,bool require_budget) {
    FlankingCheck out;
    if(!evaluate || !chain.charts || chain.charts->empty()) {out.ok=true;return out;}
    const auto plan=flanking_plan(*before,*after,chain);
    auto measured=(*evaluate)(after,plan.ids,{});
    if(!all_measured(plan,measured)) {out.after=std::move(measured);out.ok=false;out.reason="route_not_evaluated";return out;}
    const auto need=missing_before(*before,chain,previous);
    if(!need.empty())for(auto& [ch,m]:(*evaluate)(before,need,{}))previous[before->chart_key(ch).digest()]=m;
    return flanking_judge(*before,plan,std::move(measured),previous,require_budget);
}
std::vector<RouteWeights> weight_variants(const PathOptions& o) {
    std::vector<RouteWeights> out;
    for(const auto& m:{RouteWeights{1,1,1},RouteWeights{0.5,2,1},RouteWeights{1.5,1,0.5},RouteWeights{1,1,1.0/6.0}})
        out.push_back({o.w_turn*m[0],o.w_flow*m[1],o.w_drift*m[2]});
    return out;
}
namespace {
// A route of the schedule that passed every check before the flanking judgement.
struct Candidate {
    std::int64_t rings=0;RouteWeights weights{};std::size_t slot=0; // slot: its entry in `tried`
    std::vector<MeshIndex> verts;double cost=0,old_cost=0;
    std::vector<std::uint8_t> seams;std::shared_ptr<const CutSurface> after;std::vector<MeshIndex> dirty;
};
// Parallel width of the speculative flanking judgement: 0 keeps the serial schedule (no
// pool, one thread, inside a worker task, or a request that does not allow it).
std::size_t speculation_width(const RefineRequest& q) {
    if(!q.speculative || !q.evaluate || !*q.evaluate || !q.previous || in_worker_task())return 0;
    const auto* pool=active_workers();
    return pool && pool->threads()>1?static_cast<std::size_t>(pool->threads()):0;
}
RefineOutcome refine_chain_impl(const RefineRequest& q,const Chain& chain,RoutingWorkspace& ws) {
    const auto& m=q.cut->mesh();const auto& o=q.options;
    RefineOutcome out;auto& r=out.result;r.chain=chain.id;r.old_route=chain.verts;r.turns_old=chain_turns(m,transfer(chain.verts));
    out.seams.assign(q.seams.begin(),q.seams.end());
    if(!chain.reason.empty()) {r.reason=chain.reason;return out;}
    // Source: flanking_ok(previous=None) starts from an empty dict on every call, so
    // without a caller cache no before-report is shared between candidates (Gate A F32).
    auto previous=[&q](FlankingCache& fresh)->FlankingCache& {return q.previous?*q.previous:fresh;};
    auto attempted=[&](std::int64_t rings,const RouteWeights& w,std::string outcome) {r.tried.push_back({rings,w,std::move(outcome)});};
    auto optimized=[&](Candidate& c) {
        r.status="optimized";r.reason.clear();r.new_route=c.verts;r.variant=c.weights;r.cost_old=c.old_cost;r.cost_new=c.cost;
        r.turns_new=chain_turns(m,transfer(c.verts));r.rings=c.rings;r.cut=c.after;r.dirty=std::move(c.dirty);
        out.seams=std::move(c.seams);
    };
    // GB performance: with a worker pool, routes that reach the flanking judgement wait in a
    // window (the first window holds one route, later ones the pool width) whose flanking
    // charts are measured in one parallel call; the window is then judged in schedule order.
    // A pending route keeps its place in `tried`; the first route that passes ends the call
    // exactly as in the serial loop, and every attempt recorded after it is dropped. Checks
    // before the judgement run in schedule order on this thread; a failure there waits
    // until the pending routes before it are judged, as the serial loop would reach it.
    const auto width=speculation_width(q);
    std::vector<Candidate> pending;std::size_t window=1;
    // Judges the pending routes; true when one passed (the outcome is then final).
    auto resolve=[&]() {
        if(pending.empty())return false;
        auto batch=std::move(pending);pending.clear();window=width;
        std::vector<std::shared_ptr<const CutSurface>> afters;for(const auto& c:batch)afters.push_back(c.after);
        FlankingWindow measured(q.cut,afters,chain,*q.evaluate,*q.previous);
        for(std::size_t k=0;k<batch.size();++k) {
            auto& c=batch[k];
            const auto flanking=measured.judge(k,*q.previous,o.require_budget);
            if(!flanking.ok) {r.tried[c.slot].outcome=flanking.reason;continue;}
            r.tried.resize(c.slot);optimized(c);return true;
        }
        return false;
    };
    // FIX-02: a starting radius above the maximum still searches (source: skipped).
    const auto max_rings=q.fixes.corridor_bounds?std::max(o.corridor_rings,o.corridor_max_rings):o.corridor_max_rings;
    if(q.fixes.corridor_bounds)r.max_rings=max_rings;
    try {
        for(auto rings=o.corridor_rings;rings<=max_rings;++rings) {
            const auto corridor=build_corridor(*q.cut,chain,rings,q.protected_edges);
            for(const auto& weights:q.fixes.weight_variants?weight_variants(o):o.variants) {
                const auto found=find_path(*q.cut,chain,corridor,weights,o,q.exposure,ws);
                if(!found.verts) {attempted(rings,weights,found.reason);continue;}
                const auto& verts=*found.verts;
                const double old_cost=route_cost(*q.cut,chain,chain.verts,corridor,weights,o,q.exposure,&ws);
                if(verts==chain.verts) {attempted(rings,weights,"same_route");continue;}
                if(o.require_improvement && found.cost>=old_cost-1e-9) {attempted(rings,weights,"no_improvement");continue;}
                if(!route_covers_chain(m,corridor,verts,chain,rings)) {attempted(rings,weights,"old_route_not_covered");continue;}
                auto seams=apply_route(m,q.seams,chain,verts);
                if(!seams) {attempted(rings,weights,"route_edges_missing");continue;}
                const auto band=corridor_band_faces(m,corridor);
                auto topology=topology_preserved(q.cut,*seams,chain,std::span<const std::uint8_t>(band));
                if(!topology.ok) {attempted(rings,weights,topology.reason);continue;}
                Candidate c{rings,weights,r.tried.size(),verts,found.cost,old_cost,std::move(*seams),topology.after,corridor.verts};
                if(width) {
                    attempted(rings,weights,""); // placeholder until judged
                    pending.push_back(std::move(c));
                    if(pending.size()>=window && resolve())return out;
                    continue;
                }
                FlankingCache fresh;const auto flanking=flanking_ok(q.cut,topology.after,chain,q.evaluate,previous(fresh),o.require_budget);
                if(!flanking.ok) {attempted(rings,weights,flanking.reason);continue;}
                optimized(c);return out;
            }
            // Legacy: a same-route result stops the widening. FIX-03 finishes the configured
            // schedule: the same path at a narrow radius is no proof a wider one cannot improve.
            // (Only checks before the judgement give same_route, so pending routes never matter.)
            if(!q.fixes.full_schedule && std::any_of(r.tried.begin(),r.tried.end(),[](const auto& t) {return t.outcome=="same_route";}))break;
        }
    } catch(...) {
        // The serial loop reaches this failure only after judging the routes before it.
        if(resolve())return out;
        throw;
    }
    if(resolve())return out;
    const bool clean=std::any_of(r.tried.begin(),r.tried.end(),[](const auto& t) {return t.outcome=="same_route" || t.outcome=="no_improvement";});
    r.reason=clean?"already_clean":"no_valid_route";
    return out;
}
void validate(const RefineRequest& q) {
    require(bool(q.cut),"routing.cut","missing cut surface");
    require(q.seams.size()==q.cut->mesh().tables().edge_vertices.size(),"routing.seams","seam mask does not match the mesh edges");
    // Legacy searches any radius as the source does (a negative radius is a chain-only
    // corridor, like 0); FIX-02 requires positive radii.
    if(q.fixes.corridor_bounds)
        require(q.options.corridor_rings>0 && q.options.corridor_max_rings>0,"routing.options","corridor rings must be positive integers");
}
}
RefineOutcome refine_chain(const RefineRequest& q,const Chain& chain) {
    validate(q);RoutingWorkspace local;return refine_chain_impl(q,chain,q.workspace?*q.workspace:local);
}
RefineOutcome refine_cycle(const RefineRequest& q,const Chain& chain) {
    validate(q);RoutingWorkspace local_ws;auto& ws=q.workspace?*q.workspace:local_ws;
    const auto& m=q.cut->mesh();const auto& o=q.options;
    RefineOutcome out;auto& r=out.result;r.chain=chain.id;r.reason="closed_cycle";r.old_route=chain.verts;r.cycle=true;
    r.turns_old=chain_turns(m,transfer(chain.verts));out.seams.assign(q.seams.begin(),q.seams.end());
    if(!chain.cycle || chain.locked || !chain.charts || chain.edges.size()<6) {r.reason=chain.locked?"locked":"cycle_not_editable";return out;}
    if(q.protected_edges)for(auto e:chain.edges)if((*q.protected_edges)[e]) {r.reason="locked";return out;}
    auto vertices=chain.verts;if(vertices.back()==vertices.front())vertices.pop_back();
    const auto n=vertices.size();
    std::vector<double> arc{0.0};for(auto e:chain.edges)arc.push_back(arc.back()+m.geometry().edge_length[e]);
    std::vector<std::size_t> anchors{0};
    for(const double fraction:{1.0/3.0,2.0/3.0}) {
        const double target=fraction*arc.back();
        const auto at=static_cast<std::size_t>(std::lower_bound(arc.begin(),arc.end(),target)-arc.begin());
        anchors.push_back(std::clamp<std::size_t>(at,1,n-1));
    }
    if(std::set<std::size_t>(anchors.begin(),anchors.end()).size()!=3) {r.reason="cycle_too_short";return out;}
    auto current=out.seams;auto surface=q.cut;std::set<MeshIndex> dirty;
    // Source: the arcs and the cycle judgement pass `previous` through (None stays None).
    for(std::size_t k=0;k<3;++k) {
        const auto start=anchors[k],end=k<2?anchors[k+1]:n;
        Chain part;part.id=chain.id;
        part.verts.assign(vertices.begin()+static_cast<std::ptrdiff_t>(start),vertices.begin()+static_cast<std::ptrdiff_t>(end));part.verts.push_back(vertices[end%n]);
        part.edges.assign(chain.edges.begin()+static_cast<std::ptrdiff_t>(start),chain.edges.begin()+static_cast<std::ptrdiff_t>(end));
        std::set<MeshIndex> sides;for(auto f:m.tables().edge_faces[part.edges.front()])sides.insert(surface->tables().chart_of_face[f]);
        part.charts=std::vector<MeshIndex>(sides.begin(),sides.end());
        RefineRequest sub=q;sub.cut=surface;sub.seams=current;
        auto record=refine_chain_impl(sub,part,ws);
        if(record.result.status=="optimized") {
            current=std::move(record.seams);surface=record.result.cut;dirty.insert(record.result.dirty.begin(),record.result.dirty.end());
        }
        record.result.tried.clear();r.parts.push_back(std::move(record.result));
    }
    if(std::equal(current.begin(),current.end(),q.seams.begin(),q.seams.end())) {r.reason="already_clean";return out;}
    const SeamGraph graph(q.cut->mesh_owner(),current,q.protected_edges?*q.protected_edges:std::span<const std::uint8_t>{},{},surface);
    std::set<MeshIndex> anchor_set;for(auto i:anchors)anchor_set.insert(vertices[i]);
    std::vector<const Chain*> loops;
    for(const auto& c:graph.chains()) {
        if(!c.cycle)continue;
        const std::set<MeshIndex> vs(c.verts.begin(),c.verts.end());
        if(std::includes(vs.begin(),vs.end(),anchor_set.begin(),anchor_set.end()))loops.push_back(&c);
    }
    if(loops.size()!=1) {r.reason="cycle_changed_topology";return out;}
    const auto corridor=build_corridor(*q.cut,chain,q.fixes.corridor_bounds?std::max(o.corridor_rings,o.corridor_max_rings):o.corridor_max_rings,
                                       q.protected_edges);
    const auto band=corridor_band_faces(m,corridor);
    auto topology=topology_preserved(q.cut,current,chain,std::span<const std::uint8_t>(band));
    if(!topology.ok) {r.reason=topology.reason;return out;}
    FlankingCache fresh;
    // GB performance: with a pool, the after- and before-charts are measured together.
    const auto flanking=speculation_width(q)?FlankingWindow(q.cut,{topology.after},chain,*q.evaluate,*q.previous).judge(0,*q.previous,o.require_budget):
                                             flanking_ok(q.cut,topology.after,chain,q.evaluate,q.previous?*q.previous:fresh,o.require_budget);
    if(!flanking.ok) {r.reason=flanking.reason;return out;}
    // Compare complete cycles with the same start and direction.
    std::vector<MeshIndex> route(loops[0]->verts);if(route.front()==route.back())route.pop_back();
    const auto index=static_cast<std::size_t>(std::find(route.begin(),route.end(),vertices[0])-route.begin());
    std::rotate(route.begin(),route.begin()+static_cast<std::ptrdiff_t>(index),route.end());
    auto position=[&](MeshIndex v) {return std::find(route.begin(),route.end(),v)-route.begin();};
    if(position(vertices[anchors[1]])>position(vertices[anchors[2]]))std::reverse(route.begin()+1,route.end());
    route.push_back(route.front());
    auto old_route=vertices;old_route.push_back(vertices.front());
    std::vector<std::pair<double,double>> costs,improved;
    for(const auto& w:q.fixes.weight_variants?weight_variants(o):o.variants)costs.push_back({route_cost(*q.cut,chain,old_route,corridor,w,o,q.exposure,&ws),route_cost(*q.cut,chain,route,corridor,w,o,q.exposure,&ws)});
    for(const auto& c:costs)if(c.second<c.first-1e-9)improved.push_back(c);
    if(o.require_improvement && improved.empty()) {r.reason="no_cycle_improvement";return out;}
    const auto& pool=improved.empty()?costs:improved;std::size_t best=0;
    for(std::size_t i=1;i<pool.size();++i)if(pool[i].first-pool[i].second>pool[best].first-pool[best].second)best=i;
    r.status="optimized";r.reason.clear();r.new_route=route;r.cost_old=pool[best].first;r.cost_new=pool[best].second;
    r.turns_old=chain_turns(m,transfer(old_route));r.turns_new=chain_turns(m,transfer(route));
    r.cut=topology.after;r.dirty.assign(dirty.begin(),dirty.end());out.seams=std::move(current);
    return out;
}
}
