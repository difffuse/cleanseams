// SPDX-License-Identifier: GPL-3.0-or-later
// Phase B (T33, Section 14.3): measured necessity of dead-end seams and the validated
// unnecessary-spur removal of the Section 21 schedule. Every removal is a trial state of the
// T28 transaction; nothing here measures quality or commits on its own. enhanced_v1 only.
#include "cleanseams/necessity.hpp"
#include "cleanseams/symmetric_moves.hpp"
#include "cleanseams/errors.hpp"
#include "cleanseams/metrics.hpp"
#include "cleanseams/detail/blake2b.hpp"
#include <algorithm>
#include <cmath>
#include <set>
namespace cleanseams {
namespace {
using Json=nlohmann::json;
Json finite_or_null(double x) {return std::isfinite(x)?Json(x):Json(nullptr);}
std::string digest_of(const Json& j) {
    const auto bytes=j.dump();detail::Blake2b hash(32);
    hash.update({reinterpret_cast<const std::uint8_t*>(bytes.data()),bytes.size()});return hash.hex_digest();
}
// Seam mask of a state as the pipeline reports it (effective cuts minus the physical boundary).
EdgeMask reported_seams(const MeshSnapshot& m,std::span<const std::uint8_t> cuts) {
    const auto& t=m.tables();EdgeMask seams(cuts.size(),0);
    for(std::size_t e=0;e<cuts.size();++e)seams[e]=cuts[e] && !t.edge_boundary[e];
    return seams;
}
MeshIndex arc_chart(const RunState& s,const DeadEndArc& arc) {
    const auto& t=s.surface->mesh().tables();return s.surface->tables().chart_of_face[t.edge_faces[arc.edges.front()].front()];
}
// Edges a committed hotspot_relief proposal added, from the state's provenance.
std::set<MeshIndex> relief_edges(const RunState& s) {
    std::set<MeshIndex> out;
    for(const auto& record:s.provenance)
        if(record.is_object() && record.value("generator","")=="hotspot_relief" && record.contains("added_edges"))
            for(const auto& e:record.at("added_edges"))out.insert(e.get<MeshIndex>());
    return out;
}
Proposal removal(const RunState& s,const DeadEndArc& arc,const std::string& initialization,std::uint64_t id) {
    Proposal p;p.identity.id=id;p.identity.parent_revision=s.revision;p.identity.generator="seam_necessity";p.identity.generator_version=seam_necessity_version;
    for(auto e:arc.edges)p.identity.removed_edges.push_back(e);
    p.initialization=initialization;p.provenance={{"kind","dead_end_removal"},{"tips",arc.tips},{"edges",arc.edges},{"junction",arc.junction}};
    return p;
}
SeamEndEvidence classify(const RunState& s,const DeadEndArc& arc,ProposalContext& ctx,const std::string& key,const std::string& initialization,std::uint64_t id) {
    SeamEndEvidence ev;ev.arc=arc;ev.key=key;
    const auto charts=s.surface->tables().n_charts;(void)charts;
    if(ctx.deadline.expired()) {ev.classification=SeamEndClass::unmeasured;ev.reason="deadline";return ev;}
    EvaluatedProposal e;
    // GB (m4): a removal trial that fails leaves the end unmeasured, with the error recorded.
    try {e=evaluate_proposal(s,removal(s,arc,initialization,id),ctx);}
    catch(const std::exception& error) {
        if(!isolated_failure(error))throw;
        auto entry=failure_record(error);ev.measured={{"rejection","trial_failure"},{"detail",entry}};
        entry["stage"]="trial";entry["generator"]="seam_necessity";entry["id"]=id;
        if(!ctx.diagnostics.contains("failures"))ctx.diagnostics["failures"]=Json::array();
        ctx.diagnostics["failures"].push_back(std::move(entry));
        ev.classification=SeamEndClass::unmeasured;ev.reason="trial_failure";return ev;
    }
    const auto chart=arc_chart(s,arc);const auto& before=s.charts[chart];
    ev.measured={{"rejection",e.rejection},{"detail",e.detail},{"solves",e.solves},{"state_revision",s.revision},{"trial_revision",e.trial.revision},
                 {"state",{{"J",s.objective.J},{"compliant",s.budget.compliant},{"worst_ratio",s.budget.worst_ratio},
                           {"chart_budget_ratio",before.budget_ratio},{"chart_anisotropy_p95",finite_or_null(before.anisotropy_p95)}}}};
    // T38: under strict symmetry the removal is paired with its mirror image; a pair the user's
    // symmetry policy refuses (symmetry_conflict) keeps the end like a hard constraint, no solve.
    static const std::set<std::string> constraint{"locked_edge","required_edge","fixed_anchor","pin_conflict","physical_boundary","forbidden_edge","endpoint_not_allowed",
                                                  "symmetry_conflict"};
    // GB N6: a removal whose solve failed (or that named an initialization the solver does not
    // honour) measured nothing about the end: it is unmeasured, never a necessary relief end.
    static const std::set<std::string> unmeasured{"solve_budget_exhausted","candidate_budget_exhausted","deadline","stale_parent","duplicate","unrepresentable_edge",
                                                  "solver_failure","unsupported_initialization"};
    if(!e.rejection.empty()) {
        ev.reason=e.rejection;
        ev.classification=constraint.contains(e.rejection)?SeamEndClass::protected_end:unmeasured.contains(e.rejection)?SeamEndClass::unmeasured:SeamEndClass::necessary_relief_end;
        return ev;
    }
    // Admissible removal: the common ordering against the same state decides.
    const CandidateView trial{id,true,&e.objective,&e.budget,e.strict_defects.size()},state{0,s.admissible(),&s.objective,&s.budget,s.strict_defects.size()};
    // GB (seam_objective.v2): with the exposure-sampling margin, as the spur-removal commit uses.
    const int order=compare_candidates(trial,state,ctx.options.enhanced.objective,objective_margin(ctx,e.trial.cuts,s.cuts,s.objective.area,s.objective.J));
    Json changed=Json::array();
    for(auto ch:e.changed_charts)changed.push_back(e.charts[ch].to_json());
    ev.measured["trial"]={{"J",e.objective.J},{"compliant",e.budget.compliant},{"worst_ratio",e.budget.worst_ratio},{"changed_charts",changed},
                          {"envelope",e.envelope.to_json()},{"ranks_before_state",order<0}};
    if(order<0) {ev.classification=SeamEndClass::removable_spur;ev.reason="removal_ranks_first";return ev;}
    ev.classification=SeamEndClass::necessary_relief_end;
    if(!e.budget.compliant && s.budget.compliant)ev.reason="over_budget";
    else if(!e.budget.compliant)ev.reason="budget_worse";
    else ev.reason="objective";
    return ev;
}
std::uint64_t arc_id(const DeadEndArc& arc,std::uint64_t base) {return base+arc.edges.front();}
}
std::string_view seam_end_class_name(SeamEndClass c) {
    switch(c) {
        case SeamEndClass::necessary_relief_end: return "necessary_relief_end";
        case SeamEndClass::removable_spur: return "removable_spur";
        case SeamEndClass::protected_end: return "protected_end";
        default: return "unmeasured";
    }
}
nlohmann::json SeamEndEvidence::to_json() const {
    return {{"tips",arc.tips},{"edges",arc.edges},{"junction",arc.junction},{"floating",arc.floating},{"key",key},
            {"classification",seam_end_class_name(classification)},{"reason",reason},{"relief",relief},{"cached",cached},{"measured",measured}};
}
nlohmann::json SeamNecessity::to_json() const {
    Json list=Json::array();for(const auto& a:arcs)list.push_back(a.to_json());
    return {{"version",version},{"revision",revision},{"cut_digest",cut_digest},{"raw_dangling_ends",raw_dangling_ends},
            {"necessary_relief_ends",necessary_relief_ends},{"removable_spur_ends",removable_spur_ends},{"protected_ends",protected_ends},
            {"unmeasured_ends",unmeasured_ends},{"trials",trials},{"cache_hits",cache_hits},{"arcs",list}};
}
namespace {
// GB m9: the run-origin envelope domain a removal trial is compared on besides the arc's own
// chart: every chart of the state that differs from the origin (its key and its map), the
// origin's revision and the relief rebases the state adopted (m10). A removal's
// quality_envelope outcome depends on all of them, so the evidence is keyed by them. Computed
// once per classified state; the arc's chart is left out per arc (its topology and map are
// keyed on their own).
struct EnvelopeDomain {
    Revision origin=0;
    std::vector<std::size_t> rebases;
    std::vector<std::pair<MeshIndex,std::string>> charts;   // state chart, digest of (key, map)
};
EnvelopeDomain envelope_domain(const ProposalContext& ctx,const RunState& s) {
    EnvelopeDomain d;d.origin=ctx.envelope.origin_revision();
    for(auto k:relief_rebases(s,Proposal{}))if(k<ctx.envelope.staged_rebases())d.rebases.push_back(k);
    const auto& m=s.surface->mesh();const auto& ct=s.surface->tables();
    for(auto ch:ctx.envelope.domain_charts(*s.surface,s.uv)) {
        Json uv=Json::array();
        for(auto f:ct.chart_faces[ch])for(auto c:m.face_corners(f))uv.push_back({s.uv[c][0],s.uv[c][1]});
        d.charts.emplace_back(ch,digest_of(Json::array({s.surface->chart_key(ch).digest(),digest_of(uv)})));
    }
    return d;
}
Json envelope_part(const EnvelopeDomain& d,const RunState& s,MeshIndex arc_chart) {
    std::vector<std::string> charts;for(const auto& [ch,digest]:d.charts)if(ch!=arc_chart)charts.push_back(digest);
    std::sort(charts.begin(),charts.end());
    Json part={{"origin_revision",d.origin},{"relief_rebases",d.rebases},{"domain",digest_of(charts)}};
    // The common ordering compares an over-budget state's worst ratio and total excess before J.
    if(!s.budget.compliant)part["state_budget"]={s.budget.worst_ratio,s.budget.total_excess};
    return part;
}
nlohmann::json key_parts(const ProposalContext& ctx,const RunState& s,const DeadEndArc& arc,const std::string& initialization,const EnvelopeDomain& domain);
}
nlohmann::json necessity_key_parts(const ProposalContext& ctx,const RunState& s,const DeadEndArc& arc,const std::string& initialization) {
    if(!s.surface || arc.edges.empty())throw ContractError(ErrorCode::invalid_input,"necessity.arc","a dead-end arc of a measured state is required");
    return key_parts(ctx,s,arc,initialization,envelope_domain(ctx,s));
}
namespace {
nlohmann::json key_parts(const ProposalContext& ctx,const RunState& s,const DeadEndArc& arc,const std::string& initialization,const EnvelopeDomain& domain) {
    if(!s.surface || arc.edges.empty())throw ContractError(ErrorCode::invalid_input,"necessity.arc","a dead-end arc of a measured state is required");
    const auto& m=s.surface->mesh();const auto& ct=s.surface->tables();const auto chart=arc_chart(s,arc);
    const auto& o=ctx.options;
    Json pins=Json::array(),uv=Json::array(),exempt=Json::array();
    for(auto f:ct.chart_faces[chart]) {
        if(f<ctx.exempt.size() && ctx.exempt[f])exempt.push_back(f);
        for(auto c:m.face_corners(f)) {
            uv.push_back({s.uv[c][0],s.uv[c][1]});
            if(const auto found=ctx.constraints.pins.find(c);found!=ctx.constraints.pins.end())pins.push_back({c,found->second[0],found->second[1]});
        }
    }
    const auto& ob=o.enhanced.objective;
    // Hard constraints a removal of this arc meets: its edges' locked/required/forbidden
    // flags and the fixed anchors on its vertices (the only ones its removal can violate).
    Json hard=Json::array();
    for(auto e:arc.edges)hard.push_back({e,e<ctx.constraints.locked.size() && ctx.constraints.locked[e],
                                         e<ctx.constraints.required.size() && ctx.constraints.required[e],
                                         e<ctx.constraints.forbidden.size() && ctx.constraints.forbidden[e]});
    Json anchors=Json::array();
    for(auto v:arc.vertices) {
        if(std::binary_search(ctx.constraints.anchors.begin(),ctx.constraints.anchors.end(),v))anchors.push_back({v,"anchor"});
        if(std::find(ctx.constraints.seam_anchors.begin(),ctx.constraints.seam_anchors.end(),v)!=ctx.constraints.seam_anchors.end())anchors.push_back({v,"seam_anchor"});
    }
    Json key={{"version",seam_necessity_version},{"geometry",m.geometry_key()},{"constraints",{{"edges",hard},{"anchors",anchors}}},
            {"topology",{{"chart_key",s.surface->chart_key(chart).digest()},{"arc",arc.edges}}},
            {"solver",{{"policy",o.policy},{"method",o.solver.method},{"iterations",o.solver.iterations},{"backend",o.solver.backend},
                       {"arap_iterations",o.layout.arap_iterations},{"initialization",initialization}}},
            {"pins",digest_of(pins)},
            {"policy",{{"enhanced",enhanced_policy_version},{"max_anisotropy_p95",o.layout.max_anisotropy_p95},{"max_area_p95",o.layout.max_area_p95},
                       {"max_anisotropy_p99",o.layout.max_anisotropy_p99},{"max_anisotropy_max",o.layout.max_anisotropy_max},
                       {"split_min_faces",o.layout.split_min_faces},{"final_gate_tolerance",o.final_gate_tolerance},{"p95_growth",o.p95_growth},
                       {"max_growth",o.max_growth},{"objective",{ob.w_length,ob.w_exposure,ob.w_turn,ob.w_fragment,ob.chart_count_weight,ob.tie_band}},{"objective_version",ob.version},
                       {"min_island_faces",o.layout.min_island_faces},{"min_island_fraction",o.layout.min_island_fraction},
                       {"overlap_area_tolerance",o.enhanced.overlap_area_tolerance},{"pin_tolerance",o.enhanced.pin_tolerance},
                       {"constraints_revision",ctx.constraints.revision},{"exempt",digest_of(exempt)}}},
            {"map",digest_of(uv)}};
    // T38: a strict-symmetry trial also removes the mirror arc (another chart), so its evidence
    // belongs to the whole cut state and the symmetry options; other modes keep the key above.
    if(strict_symmetry(ctx)) {
        const auto& sym=o.enhanced.symmetry;
        key["symmetry"]={{"mode",sym.mode},{"state",cut_digest(s.cuts)},{"plane_normal",sym.plane_normal?Json(*sym.plane_normal):Json(nullptr)},
                         {"plane_offset",sym.plane_offset},{"map_tolerance",sym.map_tolerance},{"geometry_tolerance",sym.geometry_tolerance}};
    }
    key["envelope"]=envelope_part(domain,s,chart);
    return key;
}
}
std::string necessity_key(const ProposalContext& ctx,const RunState& s,const DeadEndArc& arc,const std::string& initialization) {
    return digest_of(necessity_key_parts(ctx,s,arc,initialization));
}
const SeamEndEvidence* SeamNecessityCache::find(const std::string& key) const {
    const auto found=entries_.find(key);
    if(found==entries_.end()) {++misses_;return nullptr;}
    ++hits_;return &found->second;
}
void SeamNecessityCache::store(const SeamEndEvidence& evidence) {
    // Unmeasured ends are not evidence: they are measured again when budgets allow.
    if(evidence.classification!=SeamEndClass::unmeasured)entries_.insert_or_assign(evidence.key,evidence);
}
SeamNecessity measure_seam_necessity(const RunState& s,ProposalContext& ctx,SeamNecessityCache& cache,const std::string& initialization) {
    if(!s.surface)throw ContractError(ErrorCode::invalid_input,"necessity.state","a measured state is required");
    const auto& m=*ctx.mesh;SeamNecessity out;out.revision=s.revision;out.cut_digest=cut_digest(s.cuts);
    const auto seams=reported_seams(m,s.cuts);out.raw_dangling_ends=dangling_seam_ends(m,seams);
    const auto relief=relief_edges(s);
    const auto domain=envelope_domain(ctx,s); // GB m9: once per classified state
    for(const auto& arc:dead_end_arcs(m,seams)) {
        throw_if_cancelled(ctx.cancel);
        const auto key=digest_of(key_parts(ctx,s,arc,initialization,domain));
        SeamEndEvidence ev;
        if(const auto* hit=cache.find(key)) {ev=*hit;ev.cached=true;++out.cache_hits;}
        else {ev=classify(s,arc,ctx,key,initialization,arc_id(arc,33'000'000));++out.trials;cache.store(ev);}
        ev.arc=arc;ev.relief=std::any_of(arc.edges.begin(),arc.edges.end(),[&](MeshIndex e) {return relief.contains(e);});
        const auto ends=static_cast<std::uint64_t>(arc.tips.size());
        switch(ev.classification) {
            case SeamEndClass::necessary_relief_end: out.necessary_relief_ends+=ends;break;
            case SeamEndClass::removable_spur: out.removable_spur_ends+=ends;break;
            case SeamEndClass::protected_end: out.protected_ends+=ends;break;
            default: out.unmeasured_ends+=ends;break;
        }
        out.arcs.push_back(std::move(ev));
    }
    return out;
}
namespace {
struct SpurState {std::shared_ptr<SeamNecessityCache> cache=std::make_shared<SeamNecessityCache>();Revision last_revision=0;};
std::vector<Proposal> spur_round(const RunState& incumbent,ProposalContext& ctx,SpurState& state,std::int64_t round) {
    auto& tally=ctx.counters.by_generator["spur_removal"];
    if(round>1 && state.last_revision==incumbent.revision) {++tally["unchanged_incumbent"];return {};}
    state.last_revision=incumbent.revision;
    auto& cache=ctx.necessity?*ctx.necessity:*state.cache;
    const auto necessity=measure_seam_necessity(incumbent,ctx,cache);
    tally["dead_end_arcs"]+=necessity.arcs.size();tally["necessity_trials"]+=necessity.trials;tally["necessity_cache_hits"]+=necessity.cache_hits;
    if(!ctx.options.layout.spur_removal) {++tally["spur_removal_disabled"];return {};}
    std::vector<Proposal> out;std::vector<MeshIndex> all;
    for(const auto& ev:necessity.arcs) {
        if(ev.classification!=SeamEndClass::removable_spur)continue;
        auto p=removal(incumbent,ev.arc,"",arc_id(ev.arc,34'000'000+static_cast<std::uint64_t>(round)*1'000'000));
        p.identity.generator="spur_removal";p.provenance["round"]=round;p.provenance["evidence"]=ev.key;
        out.push_back(std::move(p));all.insert(all.end(),ev.arc.edges.begin(),ev.arc.edges.end());
    }
    if(out.size()>=2) {
        Proposal p;p.identity.id=34'000'000+static_cast<std::uint64_t>(round)*1'000'000-1;p.identity.parent_revision=incumbent.revision;
        p.identity.generator="spur_removal";p.identity.generator_version=seam_necessity_version;
        std::sort(all.begin(),all.end());for(auto e:all)p.identity.removed_edges.push_back(e);
        p.provenance={{"kind","dead_end_removal_all"},{"round",round},{"arcs",out.size()}};
        out.push_back(std::move(p));
    }
    tally["removals_proposed"]+=out.size();
    return out;
}
}
std::vector<ProposalGenerator> spur_removal_generators(const PipelineOptions&,std::int64_t rounds) {
    if(rounds<0)throw ContractError(ErrorCode::invalid_input,"necessity.rounds","rounds must be nonnegative");
    auto state=std::make_shared<SpurState>();std::vector<ProposalGenerator> out;
    for(std::int64_t r=1;r<=rounds;++r)
        out.push_back({"spur_removal",seam_necessity_version,[state,r](const RunState& incumbent,ProposalContext& ctx) {return spur_round(incumbent,ctx,*state,r);}});
    return out;
}
}
