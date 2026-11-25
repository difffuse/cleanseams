// SPDX-License-Identifier: GPL-3.0-or-later
// Phase B common scoring and candidate selection (T28, Sections 11.2 and 21): the
// versioned seam objective J, the lexicographic feasibility/budget ordering with its tie
// band, bounded and deduplicated proposal rounds, and the enhancement pass that always
// keeps candidate zero. T40: the pass is a bounded alternative search (estimate-ordered
// queues, pass-wide deduplication, per-feature shares of the per-chart bounds, proven
// bounds before solves, deadlines and cancellation that never erase the incumbent) over the
// complete Section 21 schedule. Used by policy enhanced_v1 only.
#include "cleanseams/proposals.hpp"
#include "cleanseams/relief.hpp"
#include "cleanseams/necessity.hpp"
#include "cleanseams/routing.hpp"
#include "cleanseams/endpoints.hpp"
#include "cleanseams/symmetric_moves.hpp"
#include "cleanseams/physical_scale.hpp"
#include "cleanseams/errors.hpp"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <numeric>
#include <set>
namespace cleanseams {
namespace {
using Json=nlohmann::json;
Vec3 negate(Vec3 v) {return {-v[0],-v[1],-v[2]};}
// |x - y| within the relative tie band of the larger magnitude.
bool tied(double x,double y,double band) {
    if(x==y)return true;
    if(!std::isfinite(x) || !std::isfinite(y))return false;
    return std::abs(x-y)<=band*std::max(std::abs(x),std::abs(y));
}
int order(double x,double y,double band) {return tied(x,y,band)?0:(x<y?-1:1);}
// GB (seam_objective.v2): also a tie within the absolute margin U of the pair.
int order(double x,double y,double band,double margin) {
    if(margin>0 && std::isfinite(x) && std::isfinite(y) && std::abs(x-y)<=margin)return 0;
    return order(x,y,band);
}
// GB: a feasible state whose only failures are inherited strict-only defects of candidate zero
// is admissible; its defect count ranks it after states with fewer (compare_candidates).
bool admissible_state(const RunState& s) {return s.admissible();}
CandidateView view_of(const RunState& s,std::uint64_t id) {return {id,admissible_state(s),&s.objective,&s.budget,s.strict_defects.size()};}
CandidateView view_of(const EvaluatedProposal& e) {return {e.proposal.identity.id,true,&e.objective,&e.budget,e.strict_defects.size()};}
// T40: the pass-wide key of a trial: its effective cuts plus the numerical initialization
// contract (initialization identity, charts re-solved on request, and the parent revision
// whose map every unchanged chart inherits). Equal keys give the same trial state. GB (m6): the
// cuts enter as their digest; StateKeys compares the cut delta on a digest collision.
std::string pass_key(std::span<const std::uint8_t> cuts,const Proposal& p,Revision parent) {
    return state_digest_key(cuts,p.initialization,p.resolve_chart_keys,parent);
}
std::string feature_of(const Proposal& p) {return p.identity.generator.empty()?"anonymous":p.identity.generator;}
Json incumbent_record(std::uint64_t round,const std::string& generator,const RunState& s) {
    return {{"round",round},{"generator",generator},{"revision",s.revision},{"feasible",s.feasible},{"valid",s.valid},
            {"strict_defects",s.strict_defects.size()},{"J",s.objective.J},
            {"compliant",s.budget.compliant},{"worst_ratio",s.budget.worst_ratio},{"total_excess",s.budget.total_excess},{"cut_digest",cut_digest(s.cuts)}};
}
// The weights with the objective version they are evaluated under, and whether they are
// exactly that version's frozen weights. GB (M3): the version is the ordering in use (v1, or
// v2 with its exposure-sampling significance model).
Json objective_options_json(const ObjectiveOptions& o,const ProposalContext& ctx) {
    const ObjectiveOptions frozen;
    const bool is_frozen=o.w_length==frozen.w_length && o.w_exposure==frozen.w_exposure && o.w_turn==frozen.w_turn &&
                         o.w_fragment==frozen.w_fragment && o.chart_count_weight==frozen.chart_count_weight;
    Json out={{"version",o.version},{"frozen_weights",is_frozen},{"w_length",o.w_length},{"w_exposure",o.w_exposure},{"w_turn",o.w_turn},{"w_fragment",o.w_fragment},
              {"chart_count_weight",o.chart_count_weight},{"tie_band",o.tie_band}};
    if(o.version==seam_objective_v2) {
        Json s={{"model",exposure_uncertainty_version},{"z",seam_objective_v2_significance},{"active",significance_active(ctx)},
                {"practical_floor",o.significance.practical_floor}}; // GB (V2)
        if(ctx.exposure_uncertainty)s["exposure"]=ctx.exposure_uncertainty->to_json();
        else s["reason"]=ctx.raw_exposure.empty()?"exposure not measured: V is 0 and no margin applies":"no sampling model: the tie band only";
        out["significance"]=std::move(s);
    }
    return out;
}
}
SeamObjective seam_objective(const CutSurface& cs,std::span<const double> raw_exposure,std::span<const std::uint8_t> exempt,
                             const LayoutOptions& layout,const ObjectiveOptions& objective) {
    const auto& m=cs.mesh();const auto& t=m.tables();const auto& g=m.geometry();const auto& ct=cs.tables();
    const auto ne=t.edge_vertices.size(),nf=g.face_area.size(),nv=m.input().positions.size();
    if(!raw_exposure.empty() && raw_exposure.size()!=ne)throw ContractError(ErrorCode::invalid_input,"enhanced.exposure","raw exposure needs one value per edge");
    if(!exempt.empty() && exempt.size()!=nf)throw ContractError(ErrorCode::invalid_input,"enhanced.exempt","exempt mask needs one flag per face");
    SeamObjective o;o.version=objective.version;o.exposure_measured=!raw_exposure.empty();
    auto relevant=[&](std::size_t f) {return exempt.empty() || !exempt[f];};
    for(std::size_t f=0;f<nf;++f)if(relevant(f))o.area+=g.face_area[f];
    // Internal seams: effective cuts minus the physical boundary (the reported seam set).
    std::vector<std::uint8_t> seam(ne,0);
    for(std::size_t e=0;e<ne;++e)if(ct.cut[e] && !t.edge_boundary[e]) {
        seam[e]=1;++o.seam_edges;o.seam_length+=g.edge_length[e];
        // Raw exposure: hidden length adds length, never dilutes the exposed length.
        if(o.exposure_measured)o.exposed_length+=g.edge_length[e]*raw_exposure[e];
    }
    // Unwanted turns where a seam passes through a vertex: the routing turn penalty
    // (0 straight, 1 reversal) weighted by the local mean edge length.
    for(std::size_t v=0;v<nv;++v) {
        std::array<MeshIndex,2> pair{};std::size_t count=0;
        for(auto e:t.vert_edges[v])if(seam[e]) {if(count<2)pair[count]=e;++count;}
        if(count!=2)continue;
        const auto vid=static_cast<MeshIndex>(v);
        const double turn=turn_cost(negate(m.tangent_dir(vid,pair[0])),m.tangent_dir(vid,pair[1]));
        o.turn_length+=0.5*(g.edge_length[pair[0]]+g.edge_length[pair[1]])*turn;
    }
    // Small fragments: charts below the minimum-island policy (absorb.tiny_charts).
    const auto tiny=tiny_charts(m,cs,layout);
    for(MeshIndex ch=0;ch<ct.n_charts;++ch)if(tiny.tiny[ch]) {
        ++o.small_fragments;for(auto f:ct.chart_faces[ch])if(relevant(f))o.fragment_area+=g.face_area[f];
    }
    o.charts=ct.n_charts;o.components=g.n_components;o.additional_charts=o.charts>o.components?o.charts-o.components:0;
    const double root=o.area>0?std::sqrt(o.area):1.0;
    o.L=o.seam_length/root;o.V=o.exposed_length/root;o.T=o.turn_length/root;
    o.F=(o.area>0?o.fragment_area/o.area:0.0)+objective.chart_count_weight*static_cast<double>(o.additional_charts);
    o.J=objective.w_length*o.L+objective.w_exposure*o.V+objective.w_turn*o.T+objective.w_fragment*o.F;
    return o;
}
bool keeps_budget(const BudgetSummary& trial,const BudgetSummary& incumbent,double band) {
    if(incumbent.compliant)return trial.compliant;
    if(trial.compliant)return true;
    return order(trial.worst_ratio,incumbent.worst_ratio,band)<=0 && order(trial.total_excess,incumbent.total_excess,band)<=0;
}
int compare_candidates(const CandidateView& a,const CandidateView& b,const ObjectiveOptions& o,double j_margin) {
    if(a.admissible!=b.admissible)return a.admissible?-1:1;
    if(!a.admissible)return 0;
    if(!a.objective || !a.budget || !b.objective || !b.budget)throw ContractError(ErrorCode::invalid_input,"enhanced.compare","admissible candidates need measurements");
    // GB (Section 11.3 policy): validity first. An admissible state that still carries candidate
    // zero's unrepaired strict-only defects ranks after a state with fewer of them, provided (GB N4,
    // Section 11.1: budgets hold unless the starting layout is already over them) that state keeps
    // the other's budget, as the repair stage requires (keeps_budget). Otherwise removing a defect
    // buys no rank and the budget ordering below decides (it cannot tie then).
    if(a.strict_defects!=b.strict_defects) {
        const bool a_fewer=a.strict_defects<b.strict_defects;
        const auto& fewer=a_fewer?*a.budget:*b.budget;const auto& more=a_fewer?*b.budget:*a.budget;
        if(keeps_budget(fewer,more,o.tie_band))return a_fewer?-1:1;
    }
    if(a.budget->compliant!=b.budget->compliant)return a.budget->compliant?-1:1;
    if(!a.budget->compliant) {
        if(const int c=order(a.budget->worst_ratio,b.budget->worst_ratio,o.tie_band))return c;
        if(const int c=order(a.budget->total_excess,b.budget->total_excess,o.tie_band))return c;
    }
    // GB (seam_objective.v2): J differences within the pair's exposure-sampling uncertainty are
    // ties (j_margin is 0 under v1); budgets and validity above are unchanged.
    return order(a.objective->J,b.objective->J,o.tie_band,j_margin);
}
BudgetStatus budget_status(const BudgetSummary& candidate,const BudgetSummary& incumbent,const ObjectiveOptions& o) {
    if(candidate.compliant)return BudgetStatus::compliant;
    if(incumbent.compliant)return BudgetStatus::over_budget;
    const int worst=order(candidate.worst_ratio,incumbent.worst_ratio,o.tie_band);
    if(worst<0 || (worst==0 && order(candidate.total_excess,incumbent.total_excess,o.tie_band)<0))return BudgetStatus::improved_over_budget;
    return BudgetStatus::over_budget;
}
// ------------------------------------------------------------ bounded alternative search (T40)
bool DominanceBound::dominates(double J,std::uint64_t id,double incumbent_margin) const {
    // Only the cuts decide J, so these comparisons are exact for the trial; the budget and
    // validity a solve would measure can only rank the trial lower, never higher. GB (v2): the
    // incumbent's margin U is exact too (it depends on the two cut sets only); the round leader
    // is already significantly better than the incumbent, and leaders compare by J alone.
    if(!std::isfinite(J))return false;
    if(incumbent_J && order(J,*incumbent_J,tie_band,significance?incumbent_margin:0.0)>=0)return true; // a tie keeps the incumbent
    if(leader_J) {const int c=order(J,*leader_J,tie_band);if(c>0 || (c==0 && id>leader_id))return true;}
    return false;
}
nlohmann::json DominanceBound::to_json() const {
    Json out={{"tie_band",tie_band}};
    if(significance)out["significance"]=true;
    out["incumbent_J"]=incumbent_J?Json(*incumbent_J):Json(nullptr);
    out["leader_J"]=leader_J?Json(*leader_J):Json(nullptr);
    if(leader_J)out["leader_id"]=leader_id;
    return out;
}
nlohmann::json SearchLedger::to_json() const {
    // Per feature: the most any origin chart used of each share (the bounds are per chart).
    Json features=Json::object();
    auto peak=[](const std::map<std::string,std::uint64_t>& by_chart) {std::uint64_t m=0;for(const auto& [k,v]:by_chart)m=std::max(m,v);return m;};
    for(const auto& [name,by_chart]:feature_solves)features[name]["solves_by_chart"]=by_chart,features[name]["max_solves_on_a_chart"]=peak(by_chart);
    for(const auto& [name,by_chart]:feature_candidates)features[name]["candidates_by_chart"]=by_chart,features[name]["max_candidates_on_a_chart"]=peak(by_chart);
    for(const auto& [name,seconds]:feature_seconds)features[name]["seconds"]=seconds;
    return {{"version",bounded_search_version},{"rounds",rounds},{"ranked",ranked},{"reordered",reordered},{"dominated",dominated},
            {"duplicates",duplicates},{"bound_refusals",bound_refusals},{"stopped",stopped},{"stopped_proposals",stopped_proposals},
            {"distinct_states",states.size()},{"state_keys",{{"kind","cut_digest"},{"collisions",states.collisions()}}},{"features",features},{"incumbents",incumbents}};
}
// ------------------------------------------------------------ digest-keyed states (GB, m6)
std::string state_digest_key(std::span<const std::uint8_t> cuts,const std::string& initialization,std::vector<std::string> resolve_chart_keys,
                             std::optional<Revision> parent) {
    auto key=cut_digest(cuts);key.push_back('\0');key+=initialization;
    std::sort(resolve_chart_keys.begin(),resolve_chart_keys.end());
    for(const auto& k:resolve_chart_keys) {key.push_back('\0');key+=k;}
    if(parent) {key.push_back('\1');key+=std::to_string(*parent);}
    return key;
}
std::vector<MeshIndex> StateKeys::delta_of(std::span<const std::uint8_t> cuts) const {
    // The edges where `cuts` differs from the reference: small for the states of one pass, and
    // equal deltas mean equal cut masks.
    std::vector<MeshIndex> out;
    const auto n=std::max(cuts.size(),reference_.size());
    for(std::size_t e=0;e<n;++e) {
        const bool a=e<cuts.size() && cuts[e],b=e<reference_.size() && reference_[e];
        if(a!=b)out.push_back(static_cast<MeshIndex>(e));
    }
    return out;
}
StateKeys::Found StateKeys::find(const std::string& key,std::span<const std::uint8_t> cuts) const {
    const auto found=entries_.find(key);
    if(found==entries_.end())return {};
    const auto delta=delta_of(cuts);
    for(const auto& entry:found->second)if(entry.delta==delta)return {true,&entry.record};
    return {};
}
bool StateKeys::insert(const std::string& key,std::span<const std::uint8_t> cuts,nlohmann::json record) {
    if(!has_reference_) {reference_.assign(cuts.begin(),cuts.end());has_reference_=true;}
    auto delta=delta_of(cuts);
    auto& list=entries_[key];
    for(const auto& entry:list)if(entry.delta==delta)return false;
    if(!list.empty())++collisions_; // the same digest key for another cut mask
    list.push_back({std::move(delta),std::move(record)});++size_;
    return true;
}
double proposal_estimate(const MeshSnapshot& m,const Proposal& p,std::span<const double> raw_exposure,double area,const ObjectiveOptions& o) {
    const auto& t=m.tables();const auto& g=m.geometry();const auto ne=t.edge_vertices.size();
    const bool exposure=raw_exposure.size()==ne;
    double length=0,exposed=0;
    auto add=[&](std::int64_t x,double sign) {
        if(x<0 || static_cast<std::uint64_t>(x)>=ne)return;
        const auto e=static_cast<std::size_t>(x);
        if(t.edge_boundary[e])return; // the physical boundary is never an internal seam
        length+=sign*g.edge_length[e];
        if(exposure)exposed+=sign*g.edge_length[e]*raw_exposure[e];
    };
    for(auto x:p.identity.added_edges)add(x,1.0);
    for(auto x:p.identity.removed_edges)add(x,-1.0);
    const double root=area>0?std::sqrt(area):1.0;
    return (o.w_length*length+o.w_exposure*exposed)/root;
}
SearchRound select_proposal(const RunState& incumbent,std::vector<Proposal> proposals,ProposalContext& ctx) {
    SearchRound round;round.state=incumbent;
    const auto& objective=ctx.options.enhanced.objective;
    auto zero=run_state_summary(incumbent);zero["id"]=0;zero["candidate"]="incumbent";zero["admissible"]=admissible_state(incumbent);
    round.candidates.push_back(std::move(zero));
    // Stable order: proposal IDs, then generation order.
    std::stable_sort(proposals.begin(),proposals.end(),[](const Proposal& a,const Proposal& b) {return a.identity.id<b.identity.id;});
    // T38: under strict symmetry each proposal is paired with its mirror image first, so a move
    // and its image are one candidate state below (evaluate_proposal pairs again, idempotently,
    // and refuses the conflicts left unpaired here).
    if(strict_symmetry(ctx)) {
        const auto map=context_symmetry(ctx);
        for(auto& p:proposals)if(auto paired=pair_proposal(p,incumbent,ctx.constraints,*map);paired.ok())p=std::move(paired.proposal);
    }
    // T40: within a pass the queue is ordered by the cheap exact estimate (ties: the stable
    // ids above); without a pass ledger every proposal is evaluated in id order (T28).
    SearchLedger* ledger=ctx.search.get();
    std::vector<std::size_t> queue(proposals.size());std::iota(queue.begin(),queue.end(),std::size_t{0});
    std::vector<double> estimate(proposals.size(),0.0);
    if(ledger) {
        ++ledger->rounds;
        for(std::size_t i=0;i<proposals.size();++i)estimate[i]=proposal_estimate(*ctx.mesh,proposals[i],ctx.raw_exposure,incumbent.objective.area,objective);
        std::stable_sort(queue.begin(),queue.end(),[&](std::size_t a,std::size_t b) {return estimate[a]<estimate[b];});
        ledger->ranked+=queue.size();
        for(std::size_t i=0;i<queue.size();++i)ledger->reordered+=queue[i]!=i;
        ledger->states.insert(pass_key(incumbent.cuts,Proposal{},incumbent.revision),incumbent.cuts,Json{{"id",0},{"generator","incumbent"},{"round",ledger->rounds}});
    }
    // The incumbent is the caller's input: a state that was never measured is a contract error
    // (it propagates, m4), not a failure of each trial.
    if(!proposals.empty() && (!incumbent.surface || incumbent.revision==0))
        throw ContractError(ErrorCode::invalid_input,"enhanced.parent","proposals need a measured parent state");
    const auto& search=ctx.options.enhanced.search;
    StateKeys seen;seen.insert(state_digest_key(incumbent.cuts,"",{},std::nullopt),incumbent.cuts,Json{{"id",0}});
    // GB (M3, seam_objective.v2): only candidates that rank strictly before the incumbent with
    // the exposure-sampling margin compete for the round; among them the best J wins (ties: the
    // lowest id). v1 keeps T40's round exactly (margin 0).
    const bool significance=significance_active(ctx);
    std::optional<EvaluatedProposal> best;
    double best_margin=0.0;
    std::uint64_t rank=0;
    for(const auto index:queue) {
        auto& p=proposals[index];const auto gen=feature_of(p);
        Json record={{"id",p.identity.id},{"generator",p.identity.generator}};
        if(ledger) {record["estimate"]=estimate[index];record["rank"]=rank++;}
        // A deadline or cancellation stops new work; what was evaluated stays, the incumbent stays.
        if(ledger) {
            if(round.stopped.empty() && ctx.cancel.requested())round.stopped="cancelled";
            if(round.stopped.empty() && ctx.deadline.expired())round.stopped="deadline";
            if(!round.stopped.empty()) {++ledger->stopped_proposals;record["rejection"]=round.stopped;round.candidates.push_back(std::move(record));continue;}
        } else {
            throw_if_cancelled(ctx.cancel);
            if(ctx.deadline.expired()) {record["rejection"]="deadline";round.candidates.push_back(std::move(record));continue;}
        }
        // Deduplicate by effective cut state plus solver initialization (T40: across the
        // pass, with the parent map the trial inherits). A delta that does not change the
        // effective cuts is left to evaluate_proposal, which names the unrepresentable edges.
        auto cuts=incumbent.cuts;bool representable=true;
        for(const auto& list:{std::pair(&p.identity.added_edges,std::uint8_t{1}),std::pair(&p.identity.removed_edges,std::uint8_t{0})})
            for(auto x:*list.first) {
                if(x<0 || static_cast<std::uint64_t>(x)>=cuts.size() || cuts[static_cast<std::size_t>(x)]==list.second) {representable=false;break;}
                cuts[static_cast<std::size_t>(x)]=list.second;
            }
        std::string key;
        if(representable) {
            bool duplicate=false;
            if(ledger) {
                key=pass_key(cuts,p,incumbent.revision);
                if(const auto found=ledger->states.find(key,cuts);found.duplicate) {duplicate=true;record["duplicate_of"]=*found.first;++ledger->duplicates;}
            } else duplicate=!seen.insert(state_digest_key(cuts,p.initialization,{},std::nullopt),cuts,Json::object());
            if(duplicate) {
                ++ctx.counters.duplicates;++ctx.counters.by_generator[gen]["duplicate"];
                record["rejection"]="duplicate";round.candidates.push_back(std::move(record));continue;
            }
        }
        // Section 21: at most max_candidates_per_chart complete candidate states per chart;
        // T40: and each feature's share of that bound and of the solve bound.
        const auto affected=proposal_origin_charts(p,incumbent,ctx);
        std::string refusal;Json refusal_detail;
        for(const auto& chart:affected) {
            if(ctx.counters.candidates_by_chart[chart]>=static_cast<std::uint64_t>(ctx.options.enhanced.max_candidates_per_chart)) {
                refusal="candidate_budget_exhausted";refusal_detail={{"chart",chart},{"limit",ctx.options.enhanced.max_candidates_per_chart}};break;
            }
            if(!ledger)continue;
            if(ledger->feature_candidates[gen][chart]>=static_cast<std::uint64_t>(search.feature_candidates_per_chart)) {
                refusal="feature_candidate_budget_exhausted";refusal_detail={{"chart",chart},{"feature",gen},{"limit",search.feature_candidates_per_chart}};break;
            }
            if(ledger->feature_solves[gen][chart]>=static_cast<std::uint64_t>(search.feature_solves_per_chart)) {
                refusal="feature_solve_budget_exhausted";refusal_detail={{"chart",chart},{"feature",gen},{"limit",search.feature_solves_per_chart}};break;
            }
        }
        if(!refusal.empty()) {
            ++ctx.counters.by_generator[gen][refusal];if(ledger)++ledger->bound_refusals;
            record["rejection"]=refusal;record["detail"]=refusal_detail;
            if(ledger && !key.empty())ledger->states.insert(key,cuts,Json{{"id",p.identity.id},{"generator",gen},{"round",ledger->rounds},{"rejection",refusal}});
            round.candidates.push_back(std::move(record));continue;
        }
        // T40: the proven bound of this candidate (compliant, admissible incumbent and leader).
        std::optional<DominanceBound> bound;
        if(ledger && search.prune_dominated) {
            DominanceBound b;b.tie_band=objective.tie_band;b.significance=significance;
            // GB: J proves nothing against a state with strict defects (a trial with fewer ranks
            // first whatever its J), so such an incumbent or leader gives no bound.
            if(admissible_state(incumbent) && incumbent.budget.compliant && incumbent.strict_defects.empty())b.incumbent_J=incumbent.objective.J;
            if(best && best->budget.compliant && best->strict_defects.empty()) {b.leader_J=best->objective.J;b.leader_id=best->proposal.identity.id;}
            if(b.incumbent_J || b.leader_J)bound=b;
        }
        std::optional<EvaluatedProposal> evaluated;
        // GB (m4): a trial that fails (its topology, solve or measurement throws) is refused with
        // the error recorded; the incumbent and the round go on. Cancellation and the caller's
        // contract errors propagate as before.
        auto trial_failure=[&](const std::exception& error) {
            ++ctx.counters.by_generator[gen]["trial_failure"];
            record["rejection"]="trial_failure";record["detail"]=failure_record(error);
            auto entry=record["detail"];entry["stage"]="trial";entry["generator"]=gen;entry["id"]=record.at("id");
            if(!ctx.diagnostics.contains("failures"))ctx.diagnostics["failures"]=Json::array();
            ctx.diagnostics["failures"].push_back(std::move(entry));
            if(ledger && !key.empty())ledger->states.insert(key,cuts,Json{{"id",record.at("id")},{"generator",gen},{"round",ledger->rounds},{"rejection","trial_failure"}});
            round.candidates.push_back(std::move(record));
        };
        try {evaluated=evaluate_proposal(incumbent,std::move(p),ctx,bound?&*bound:nullptr);}
        catch(const ContractError& error) {
            // A cancellation inside a trial (its solve) discards the trial, never the incumbent.
            if(error.code==ErrorCode::cancelled) {
                if(!ledger)throw;
                round.stopped="cancelled";++ledger->stopped_proposals;record["rejection"]="cancelled";round.candidates.push_back(std::move(record));continue;
            }
            if(!isolated_failure(error))throw;
            trial_failure(error);continue;
        } catch(const std::exception& error) {trial_failure(error);continue;}
        auto summary=evaluated->summary();
        if(ledger) {
            summary["estimate"]=record.at("estimate");summary["rank"]=record.at("rank");
            if(evaluated->rejection=="dominated")++ledger->dominated;
            if(evaluated->solves)for(const auto& chart:affected)ledger->feature_solves[gen][chart]+=1;
            if(!key.empty())
                ledger->states.insert(key,cuts,Json{{"id",evaluated->proposal.identity.id},{"generator",gen},{"round",ledger->rounds},
                                                    {"rejection",evaluated->rejection},{"J",evaluated->objective.J}});
        }
        if(evaluated->admissible()) {
            for(const auto& chart:affected) {++ctx.counters.candidates_by_chart[chart];if(ledger)++ledger->feature_candidates[gen][chart];}
            summary["budget_status"]=budget_status_name(budget_status(evaluated->budget,incumbent.budget,objective));
            const CandidateView challenger=view_of(*evaluated);
            bool eligible=true;double margin=0.0;
            if(significance) {
                // GB (visual review V2): max(U, the practical-effect floor) for a J-decided adoption.
                const auto sig=objective_significance(ctx,evaluated->trial.cuts,incumbent.cuts,incumbent.objective.area,incumbent.objective.J);
                margin=sig.margin();
                eligible=compare_candidates(challenger,view_of(incumbent,0),objective,margin)<0;
                summary["significance"]=significance_record(ctx,evaluated->objective.J,incumbent.objective.J,sig,eligible);
                if(!eligible)++ctx.counters.by_generator[gen][summary["significance"].value("below_practical_floor",false)?"below_practical_floor":"within_sampling_uncertainty"];
            }
            if(eligible) {
                if(!best) {best=std::move(*evaluated);best_margin=margin;}
                else {
                    const CandidateView current=view_of(*best);
                    const int c=compare_candidates(challenger,current,objective);
                    if(c<0 || (c==0 && challenger.id<current.id)) {best=std::move(*evaluated);best_margin=margin;}
                }
            }
        }
        round.candidates.push_back(std::move(summary));
    }
    // The incumbent is replaced only by a candidate that ranks strictly before it; a tie
    // within the band keeps it (Section 11.2); under seam_objective.v2 a J difference within
    // the exposure-sampling uncertainty U is a tie too. A cancelled round commits nothing.
    if(best && round.stopped!="cancelled") {
        const CandidateView challenger=view_of(*best);
        if(compare_candidates(challenger,view_of(incumbent,0),objective,best_margin)<0) {
            const auto status=budget_status(best->budget,incumbent.budget,objective);
            round.state=commit_proposal(incumbent,*best,ctx);round.committed=true;
            // A repair out of an infeasible incumbent re-anchors the run-origin envelope.
            if(!incumbent.feasible && round.state.feasible)anchor_envelope(ctx,round.state);
            round.winner={{"id",best->proposal.identity.id},{"revision",round.state.revision},{"budget_status",budget_status_name(status)}};
            return round;
        }
    }
    round.winner={{"id",0},{"revision",incumbent.revision},{"budget_status",budget_status_name(budget_status(incumbent.budget,incumbent.budget,objective))}};
    return round;
}
EnhancementResult run_enhancement(const RunState& origin,std::span<const ProposalGenerator> generators,ProposalContext& ctx) {
    EnhancementResult result;result.origin=origin;result.state=origin;
    const auto& objective=ctx.options.enhanced.objective;
    // T40: the pass is a bounded alternative search with one ledger.
    if(!ctx.search)ctx.search=std::make_shared<SearchLedger>();
    auto& ledger=*ctx.search;
    ledger.incumbents.push_back(incumbent_record(0,"candidate_zero",origin));
    Json rounds=Json::array(),names=Json::array(),strict_stage=nullptr;
    // GB (m4): an isolated failure of a stage, recorded with the trials' failures.
    auto record_failure=[&](const std::string& stage,const std::exception& error) {
        auto entry=failure_record(error);entry["stage"]="generate";entry["generator"]=stage;
        if(!ctx.diagnostics.contains("failures"))ctx.diagnostics["failures"]=Json::array();
        ctx.diagnostics["failures"].push_back(std::move(entry));
        ++ctx.counters.by_generator[stage]["generator_failure"];
    };
    // GB (Section 11.3 policy): the validated strict-overlap repair runs once, as soon as the
    // incumbent is feasible and fails only the strict check (candidate zero's strict-only
    // defects), before any other proposal stage (after the constraint repair when needed).
    auto repair_strict=[&] {
        if(!strict_stage.is_null() || !ctx.options.enhanced.strict_repair.enabled || !result.state.strict_only)return;
        if(ledger.stopped.empty() && ctx.cancel.requested())ledger.stopped="cancelled";
        if(ledger.stopped.empty() && ctx.deadline.expired())ledger.stopped="deadline";
        if(!ledger.stopped.empty()) {strict_stage={{"generator","strict_overlap_repair"},{"skipped",ledger.stopped}};rounds.push_back(strict_stage);return;}
        const auto t0=std::chrono::steady_clock::now();
        SearchRound round;
        // GB (m4): a failure of the stage keeps the incumbent (its defects stay reported).
        try {round=repair_strict_overlaps(result.state,ctx);}
        catch(const std::exception& error) {
            if(!isolated_failure(error))throw;
            strict_stage={{"generator","strict_overlap_repair"},{"failed",failure_record(error)}};rounds.push_back(strict_stage);
            record_failure("strict_overlap_repair",error);return;
        }
        const double seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-t0).count();
        ledger.feature_seconds["strict_overlap_repair"]["evaluate"]+=seconds;
        Json record={{"generator","strict_overlap_repair"},{"version",strict_repair_version},{"proposals",round.candidates.size()>0?round.candidates.size()-1:0},
                     {"committed",round.committed},{"winner",round.winner},{"candidates",round.candidates},{"generate_s",0.0},{"evaluate_s",seconds}};
        if(!round.stopped.empty()) {record["stopped"]=round.stopped;ledger.stopped=round.stopped;}
        strict_stage=record;rounds.push_back(std::move(record));
        if(round.committed)ledger.incumbents.push_back(incumbent_record(ledger.rounds,"strict_overlap_repair",round.state));
        result.state=std::move(round.state);
    };
    for(const auto& generator:generators) {
        // The constraint repair (painted hard constraints) comes first; everything else after
        // the strict-overlap repair.
        if(generator.name!="constraint_repair")repair_strict();
        names.push_back({{"name",generator.name},{"version",generator.version}});
        // New work stops at a deadline or a cancellation; the incumbent is kept.
        if(ledger.stopped.empty() && ctx.cancel.requested())ledger.stopped="cancelled";
        if(ledger.stopped.empty() && ctx.deadline.expired())ledger.stopped="deadline";
        if(!ledger.stopped.empty()) {rounds.push_back({{"generator",generator.name},{"skipped",ledger.stopped}});continue;}
        const auto t0=std::chrono::steady_clock::now();
        std::vector<Proposal> proposals;
        try {proposals=generator.generate?generator.generate(result.state,ctx):std::vector<Proposal>{};}
        catch(const ContractError& error) {
            if(error.code==ErrorCode::cancelled) {
                ledger.stopped="cancelled";rounds.push_back({{"generator",generator.name},{"skipped","cancelled"},{"during","generation"}});continue;
            }
            // GB (m4): a generator that fails contributes no proposal; the pass goes on.
            if(!isolated_failure(error))throw;
            rounds.push_back({{"generator",generator.name},{"version",generator.version},{"failed",failure_record(error)},{"proposals",0},{"committed",false}});
            record_failure(generator.name,error);continue;
        } catch(const std::exception& error) {
            rounds.push_back({{"generator",generator.name},{"version",generator.version},{"failed",failure_record(error)},{"proposals",0},{"committed",false}});
            record_failure(generator.name,error);continue;
        }
        const auto t1=std::chrono::steady_clock::now();
        // GB N5: a deadline that expired while the generator ran may have cut its detection short
        // (detectors stop with their partial result): the round is a time-budget outcome, never
        // "nothing eligible", even when it produced no proposal.
        const bool generation_cut=ctx.deadline.expired();
        for(auto& p:proposals) {
            if(p.identity.generator.empty())p.identity.generator=generator.name;
            if(p.identity.generator_version.empty())p.identity.generator_version=generator.version;
        }
        const auto count=proposals.size();
        auto round=select_proposal(result.state,std::move(proposals),ctx);
        const double generate_s=std::chrono::duration<double>(t1-t0).count(),evaluate_s=std::chrono::duration<double>(std::chrono::steady_clock::now()-t1).count();
        ledger.feature_seconds[generator.name]["generate"]+=generate_s;ledger.feature_seconds[generator.name]["evaluate"]+=evaluate_s;
        Json record={{"generator",generator.name},{"version",generator.version},{"proposals",count},{"committed",round.committed},
                     {"winner",round.winner},{"candidates",round.candidates},{"generate_s",generate_s},{"evaluate_s",evaluate_s}};
        if(round.stopped.empty() && generation_cut) {round.stopped="deadline";record["deadline_during_generation"]=true;}
        if(!round.stopped.empty()) {record["stopped"]=round.stopped;ledger.stopped=round.stopped;}
        rounds.push_back(std::move(record));
        if(round.committed)ledger.incumbents.push_back(incumbent_record(ledger.rounds,generator.name,round.state));
        result.state=std::move(round.state);
    }
    repair_strict();
    // Common final comparison with candidate zero: the retained state must rank strictly
    // before the unchanged baseline, otherwise the baseline is kept. GB (seam_objective.v2): with
    // the exposure-sampling margin of the two cut sets (0 under v1).
    Json final_significance=nullptr;
    if(result.state.revision!=origin.revision) {
        const auto sig=objective_significance(ctx,result.state.cuts,origin.cuts,origin.objective.area,origin.objective.J); // GB (V2): with the floor
        const bool first=compare_candidates(view_of(result.state,1),view_of(origin,0),objective,sig.margin())<0;
        if(significance_active(ctx))final_significance=significance_record(ctx,result.state.objective.J,origin.objective.J,sig,first);
        if(!first)result.state=origin;
    }
    const bool retained=result.state.revision==origin.revision;
    const auto& o=ctx.options.enhanced;
    result.report={{"policy",enhanced_policy_version},{"objective_options",objective_options_json(objective,ctx)},
                   {"bounds",{{"max_candidates_per_chart",o.max_candidates_per_chart},{"max_solves_per_chart",o.max_solves_per_chart},
                              {"feature_candidates_per_chart",o.search.feature_candidates_per_chart},{"feature_solves_per_chart",o.search.feature_solves_per_chart},
                              {"prune_dominated",o.search.prune_dominated}}},
                   {"envelope",ctx.envelope.to_json()},{"generators",names},{"candidate_zero",run_state_summary(origin)},
                   {"retained",run_state_summary(result.state)},{"baseline_retained",retained},
                   {"budget_status",budget_status_name(budget_status(result.state.budget,origin.budget,objective))},
                   {"rounds",rounds},{"counters",ctx.counters.to_json()},{"provenance",result.state.provenance},{"diagnostics",ctx.diagnostics},
                   {"scale",physical_scale_report(ctx,origin,result.state)},{"search",ledger.to_json()},{"stopped",ledger.stopped},
                   {"strict_overlap",strict_overlap_report(ctx,origin,result.state,strict_stage)},
                   {"pin_check",pin_check_report(ctx,origin,result.state)}};
    if(!final_significance.is_null())result.report["final_significance"]=std::move(final_significance);
    // GB (m4): stages and trials whose errors were isolated (the incumbent was kept).
    result.report["failures"]=ctx.diagnostics.value("failures",Json::array());
    return result;
}
namespace {
Json step(int index,const char* section,const char* feature,bool runs,std::string reason,Json generators,Json bounds=Json::object()) {
    return {{"step",index},{"section21",section},{"feature",feature},{"runs",runs},{"reason",std::move(reason)},{"generators",std::move(generators)},
            {"bounds",std::move(bounds)}};
}
}
nlohmann::json enhanced_schedule(const PipelineOptions& options,bool painted) {
    const auto& o=options.enhanced;const ReliefSettings relief;
    auto on=[](bool enabled,const char* option) {return enabled?std::string("enabled"):std::string("disabled by ")+option;};
    Json relief_passes=Json::array();for(std::int64_t k=0;k<relief.passes;++k)relief_passes.push_back("hotspot_relief");
    const Json per_chart={{"feature_candidates_per_chart",o.search.feature_candidates_per_chart},{"feature_solves_per_chart",o.search.feature_solves_per_chart}};
    auto with=[&](Json own) {own["pass_share"]=per_chart;return own;};
    Json out=Json::array();
    // GB M5: the repair has its own switch (enhanced.constraint_repair.enabled, an ablation switch).
    const bool repair=painted && o.constraint_repair.enabled;
    out.push_back(step(1,"resolve constraints/preferences","constraint_repair",repair,
                       !painted?"no painted hard constraint":repair?"painted hard constraints: repair proposals when candidate zero violates them":
                                "disabled by enhanced.constraint_repair.enabled: a kept state violating painted hard constraints is refused",
                       repair?Json::array({"constraint_repair"}):Json::array(),with({})));
    out.push_back(step(2,"baseline generation/refinement","candidate_zero",true,
                       "the complete native_baseline pipeline (generation, refinement, contour cleanup, repair ladder, chart gate) runs unchanged; its "
                       "accepted state is candidate zero",Json::array()));
    // GB (Section 11.3 policy): candidate zero's strict-only defects get the validated repair
    // stage first (after the constraint repair when painted constraints need it); a defect it
    // cannot repair is reported and keeps native_baseline's decision.
    out.back()["strict_overlap_repair"]={{"version",strict_repair_version},{"runs",o.strict_repair.enabled},
        {"reason",o.strict_repair.enabled?"when candidate zero fails only a strict check (overlap; or a pinned UV a host solve moved, GB N2): re-maps, then "
                                          "for overlaps a minimal split, of exactly the failing charts":
                                          "disabled by enhanced.strict_repair.enabled: strict-only defects are reported unrepaired"},
        {"bounds",{{"per_chart","one solve per attempt within max_solves_per_chart and search.feature_solves_per_chart"}}}};
    out.push_back(step(3,"sliding proposals on eligible existing tube slits","endpoint_slide",o.endpoint_slide.enabled,
                       on(o.endpoint_slide.enabled,"enhanced.endpoint_slide.enabled"),o.endpoint_slide.enabled?Json::array({"endpoint_slide"}):Json::array(),
                       with({{"routes_per_slit",4},{"vertex_budget",o.endpoint_slide.vertex_budget},{"max_slide_distance",o.endpoint_slide.max_slide_distance}})));
    out.push_back(step(4,"relief candidates for over-budget charts","hotspot_relief",o.relief.enabled,on(o.relief.enabled,"enhanced.relief.enabled"),
                       o.relief.enabled?relief_passes:Json::array(),with(relief.to_json())));
    out.push_back(step(5,"preference-aware split alternatives","region_preferences",o.regions.enabled,on(o.regions.enabled,"enhanced.regions.enabled"),
                       o.regions.enabled?Json::array({"region_preferences"}):Json::array(),with({{"max_targets",o.regions.max_targets},{"variants_per_target",2}})));
    out.push_back(step(6,"confidence-gated ring/cavity/shell alternatives","rings",o.rings.enabled,on(o.rings.enabled,"enhanced.rings.enabled"),
                       o.rings.enabled?Json::array({"rings"}):Json::array(),with({{"max_charts",o.rings.max_charts},{"max_candidates",o.rings.max_candidates},
                       {"min_confidence",o.rings.min_confidence}})));
    out.push_back(step(6,"confidence-gated ring/cavity/shell alternatives","cavities",o.cavities.enabled,on(o.cavities.enabled,"enhanced.cavities.enabled"),
                       o.cavities.enabled?Json::array({"cavities"}):Json::array(),with({{"max_cavities",o.cavities.max_cavities},{"max_samples",o.cavities.max_samples},
                       {"min_confidence",o.cavities.min_confidence}})));
    out.push_back(step(6,"confidence-gated ring/cavity/shell alternatives","shells",o.shells.enabled,on(o.shells.enabled,"enhanced.shells.enabled"),
                       o.shells.enabled?Json::array({"shells"}):Json::array(),with({{"max_candidates",o.shells.max_candidates},{"max_samples",o.shells.max_samples},
                       {"min_confidence",o.shells.min_confidence}})));
    const bool symmetry=o.symmetry.mode!="asymmetric";
    out.push_back(step(7,"strict symmetry pairing where requested","symmetry",symmetry,
                       symmetry?"enhanced.symmetry.mode "+o.symmetry.mode:std::string("enhanced.symmetry.mode asymmetric (not requested)"),
                       symmetry?Json::array({"symmetry"}):Json::array(),with({{"max_candidates",o.symmetry.max_candidates},{"mode",o.symmetry.mode}})));
    Json spur_rounds=Json::array();if(o.spur_removal.enabled)for(int k=0;k<3;++k)spur_rounds.push_back("spur_removal");
    out.push_back(step(8,"validated unnecessary-spur removal","spur_removal",o.spur_removal.enabled,
                       o.spur_removal.enabled?(options.layout.spur_removal?"enabled":"measurement only (layout.spur_removal false keeps every end)"):
                                             "disabled by enhanced.spur_removal.enabled",spur_rounds,with({{"rounds",3}})));
    out.push_back(step(9,"contour cleanup and ordinary repair/gate","transaction_gate",true,
                       "candidate zero carries the baseline's contour cleanup, repair and gate; every trial passes the same validity checks, strict "
                       "overlap, pins, P95 budgets with the P99/maximum caps and the run-origin envelope; the retained state gets the output-hygiene writer",
                       Json::array()));
    out.push_back(step(10,"common final candidate comparison","final_comparison",true,
                       "the retained state replaces candidate zero only when it ranks strictly before it (lexicographic feasibility/budget, then J)",
                       Json::array(),{{"max_candidates_per_chart",o.max_candidates_per_chart},{"max_solves_per_chart",o.max_solves_per_chart},
                                      {"prune_dominated",o.search.prune_dominated},{"search",bounded_search_version}}));
    return out;
}
std::vector<ProposalGenerator> enhanced_generators(const PipelineOptions& options) {
    // Section 21 schedule (enhanced_schedule): sliding, relief, region, ring/cavity/shell,
    // symmetry, spur removal. Candidate zero is always evaluated; the constraint repair
    // (only with painted hard constraints) runs before these.
    const auto& o=options.enhanced;
    std::vector<ProposalGenerator> out;
    if(o.endpoint_slide.enabled)out.push_back(endpoint_slide_generator());   // T31: sliding proposals on eligible existing tube slits
    // T32: relief candidates for over-budget charts, before the ordinary over-budget split.
    if(o.relief.enabled)for(auto& g:hotspot_relief_generators(options))out.push_back(std::move(g));
    // T34-T37 (chain C): preference-aware split alternatives, then the confidence-gated
    // ring, cavity and shell alternatives.
    if(o.regions.enabled)out.push_back(region_preference_generator());
    // T35: constriction/crease rings as compound cuts (confidence-gated), before spur removal.
    if(o.rings.enabled)out.push_back(ring_generator());
    // T36: confidence-gated cavity rings and internal slits.
    if(o.cavities.enabled)out.push_back(cavity_generator());
    // T37: thin-shell front/back rims.
    if(o.shells.enabled)out.push_back(shell_generator());
    // T38: strict symmetry pairing where requested (the explicit enhanced.symmetry.mode;
    // approximate: metric-only diagnostics), after the shells and before spur removal.
    if(o.symmetry.mode!="asymmetric")out.push_back(symmetry_generator());
    // T33: validated unnecessary-spur removal (measured necessity), after the other operations.
    if(o.spur_removal.enabled)for(auto& g:spur_removal_generators(options))out.push_back(std::move(g));
    return out;
}
}
