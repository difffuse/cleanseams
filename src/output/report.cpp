// SPDX-License-Identifier: GPL-3.0-or-later
// Source ops.py pure helpers: summarize, refusal_message, repairs_message, _json_safe.
#include "cleanseams/report.hpp"
#include "cleanseams/symmetry.hpp"
#include "cleanseams/errors.hpp"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <map>
#include <set>
namespace cleanseams {
namespace {
using Json=nlohmann::json;
// Python f"{x:.Nf}" (correctly rounded; nan/inf spelled as Python does).
std::string fixed(double x,int digits) {
    if(std::isnan(x))return "nan";
    if(std::isinf(x))return x>0?"inf":"-inf";
    char buffer[512];std::snprintf(buffer,sizeof buffer,"%.*f",digits,x);return buffer;
}
double number(const Json& j) {return j.is_null()?std::numeric_limits<double>::quiet_NaN():j.get<double>();}
std::string plural(const char* word,std::int64_t n) {return std::string(word)+(n!=1?"s":"");}
std::string integer_text(const Json& j) {return j.is_number_integer() || j.is_number_unsigned()?std::to_string(j.get<std::int64_t>()):j.dump();}
// ------------------------------------------------------------ T41 diagnostics helpers
const Json& member(const Json& j,const char* key) {
    static const Json none;
    return j.is_object() && j.contains(key)?j.at(key):none;
}
std::string text_of(const Json& j,const char* key,const std::string& fallback="") {
    const auto& v=member(j,key);return v.is_string()?v.get<std::string>():fallback;
}
std::int64_t count_of(const Json& j,const char* key) {
    const auto& v=member(j,key);
    if(v.is_number_integer() || v.is_number_unsigned())return v.get<std::int64_t>();
    return v.is_array() || v.is_object()?static_cast<std::int64_t>(v.size()):0;
}
// Three significant digits (C "%.3g"); null (a non-finite measurement) prints "n/a".
std::string short_number(const Json& j) {
    if(!j.is_number())return "n/a";
    char buffer[64];std::snprintf(buffer,sizeof buffer,"%.3g",j.get<double>());return buffer;
}
std::string id_list(const Json& ids,std::size_t limit=6) {
    std::string out;std::size_t n=0;
    if(!ids.is_array())return out;
    for(const auto& id:ids) {if(n==limit) {out+=", ...";break;}out+=(n++?", ":"")+(id.is_number()?integer_text(id):id.dump());}
    return out;
}
std::string joined(const std::vector<std::string>& parts,const char* separator=", ") {
    std::string out;for(std::size_t i=0;i<parts.size();++i)out+=(i?separator:"")+parts[i];return out;
}
// Why the transaction or the search refused a proposal, in the words a user acts on.
std::string rejection_text(const std::string& code) {
    static const std::map<std::string,std::string> text{
        {"dominated","provably worse than the kept layout"},{"duplicate","already evaluated"},
        {"candidate_budget_exhausted","island search bound reached"},{"solve_budget_exhausted","island solve bound reached"},
        {"feature_candidate_budget_exhausted","feature share of the search bound reached"},
        {"feature_solve_budget_exhausted","feature share of the solve bound reached"},
        {"deadline","time budget reached"},{"cancelled","cancelled"},
        {"quality_envelope","more distortion than the baseline allows"},{"overlap","strict UV overlap"},
        {"invalid_uv","invalid UVs (flips, collapse or crossings)"},{"not_disk","an island would not be a disk"},
        {"pin_conflict","would glue differently pinned UVs"},{"pin_moved","would move a pinned UV"},
        {"locked_edge","would remove a locked seam"},{"required_edge","would drop a required seam"},
        {"forbidden_edge","would cut a forbidden edge"},{"fixed_anchor","would move an anchor"},
        {"endpoint_not_allowed","slit end outside its allowed rim"},{"symmetry_conflict","the mirror half breaks a constraint"},
        {"solver_failure","the unwrap failed"},{"physical_boundary","touches a physical boundary"},
        {"unrepresentable_edge","not a change of the surface's seams"},{"stale_parent","stale"},{"ranked_worse","ranked worse"},
        // GB review fixes (M3, m3, m4)
        {"within_sampling_uncertainty","no better beyond the exposure sampling uncertainty (a tie keeps the layout)"},
        {"below_practical_floor","an improvement too small to be worth a change (practical-effect floor; the layout is kept)"}, // GB (V2)
        {"unsupported_initialization","an initialization the solver does not honour"},{"trial_failure","the trial failed (error recorded)"}};
    const auto found=text.find(code);return found==text.end()?code:found->second;
}
Json diagnostic(const char* level,const char* topic,std::string text) {return {{"level",level},{"topic",topic},{"text",std::move(text)}};}
// GB m15: the run met its time budget only in the enhanced pass (it stopped there with the best
// validated layout; no stage before it timed out).
bool enhanced_deadline_only(const Json& rep) {
    for(const auto& st:member(rep,"stages"))if(member(st,"status")=="timed_out")return false;
    return text_of(member(rep,"enhanced"),"stopped")=="deadline";
}
// GB M4: an enhanced state replaced candidate zero, so the report's generation, refinement,
// contour and chart gate records describe the baseline, not the result.
bool enhanced_adopted(const Json& rep) {
    const auto& e=member(rep,"enhanced");return text_of(e,"status")=="done" && member(e,"baseline_retained")==false;
}
}
std::string repairs_message(const Json& rep) {
    const auto rs=rep.contains("repairs_summary") && rep.at("repairs_summary").is_object()?rep.at("repairs_summary"):Json::object();
    const auto repaired=rs.value("repaired",std::int64_t{0})+rs.value("exempted",std::int64_t{0});
    const auto unresolved=rs.value("unresolved",std::int64_t{0}),exempt=rs.value("exempt_faces",std::int64_t{0});
    const auto over=rep.contains("chart_over_budget") && rep.at("chart_over_budget").is_array()?static_cast<std::int64_t>(rep.at("chart_over_budget").size()):0;
    if(!(repaired || unresolved || exempt || over))return "applied, every island valid";
    std::string out="applied, "+std::to_string(repaired)+" "+plural("island",repaired)+" repaired";
    if(unresolved)out+=", "+std::to_string(unresolved)+" unresolved";
    if(exempt)out+=", "+std::to_string(exempt)+" "+plural("face",exempt)+" exempt";
    if(over)out+=", "+std::to_string(over)+" "+plural("island",over)+" over budget";
    return out;
}
std::string refusal_message(const std::string& reason,std::span<const MeshIndex> invalid) {
    static const std::map<std::string,std::string> refusals{
        {"fragment_cloud","too many disconnected parts or too little surface in the main part; see the input counts in the report"},
        {"scan_debris","more than a fifth of the faces are degenerate (scan debris)"},
        {"too_many_faces","more faces than the operator handles in one run"},
        {"chart_unrepairable","an island could not be repaired or cut apart"},
        {"unwrap_invalid","the unwrap is invalid and partial results are off"},
        {"final_invalid","the final map is invalid"},
        {"degraded_beyond_budget","the result would be more distorted than the existing UVs"},
        {"pins_changed","the candidate could not preserve every pinned UV"},
        {"enhanced_invalid","the retained map fails the enhanced validity checks (strict overlap, pins or chart validity)"},
        {"constraint_conflict","a hard seam constraint cannot be satisfied; see report.enhanced for the edge and anchor IDs"}};
    const auto found=refusals.find(reason);std::string text=found!=refusals.end()?found->second:reason;
    if(!invalid.empty()) {
        text+=" ("+plural("island",static_cast<std::int64_t>(invalid.size()))+" ";
        for(std::size_t i=0;i<std::min<std::size_t>(6,invalid.size());++i)text+=(i?", ":"")+std::to_string(invalid[i]);
        text+=")";
    }
    return reason+": "+text;
}
std::string summarize(const PipelineResult& res) {
    const auto& rep=res.report;std::vector<std::string> parts;
    // GB M4: after an adopted enhanced state the layout, chain and contour clauses are the baseline's.
    const bool adopted=enhanced_adopted(rep);const std::string baseline=adopted?"baseline ":"";
    if(rep.value("timed_out",false))
        parts.push_back(enhanced_deadline_only(rep)?"time budget reached in the enhanced pass: kept the best validated layout": // GB m15
                        "time budget exhausted: validated partial layout");
    if(rep.contains("layout") && rep.at("layout").is_object() && !rep.at("layout").empty()) {
        const auto& lay=rep.at("layout");
        parts.push_back((adopted?"baseline layout: ":"")+integer_text(lay.at("charts"))+" islands from "+integer_text(lay.at("patches"))+" patches");
        std::vector<std::pair<std::string,std::int64_t>> reasons;
        if(lay.contains("stop_reasons"))for(auto& [kind,count]:lay.at("stop_reasons").items())reasons.emplace_back(kind,count.get<std::int64_t>());
        std::sort(reasons.begin(),reasons.end(),[](const auto& a,const auto& b) {return std::pair(-a.second,a.first)<std::pair(-b.second,b.first);});
        if(!reasons.empty()) {
            std::string line="merging stopped: ";
            for(std::size_t i=0;i<std::min<std::size_t>(3,reasons.size());++i) {
                auto kind=reasons[i].first;std::replace(kind.begin(),kind.end(),'_',' ');
                line+=(i?", ":"")+kind+" ("+std::to_string(reasons[i].second)+")";
            }
            parts.push_back(line);
        }
    }
    if(rep.contains("profile") && rep.at("profile").is_object() && !rep.at("profile").empty())
        parts.push_back("auto selected "+rep.at("profile").at("selected").get<std::string>());
    if(rep.contains("island_target") && rep.at("island_target").is_object() && !rep.at("island_target").empty()) {
        const auto& target=rep.at("island_target");
        parts.push_back("island target "+target.at("status").get<std::string>()+"; budget scale "+fixed(number(target.at("scale")),2));
    }
    for(const auto& st:rep.value("stages",Json::array())) {
        if(st.at("stage")=="refine")parts.push_back(baseline+"chains: "+integer_text(st.at("optimized"))+" optimized, "+integer_text(st.at("kept"))+" kept");
        if(st.at("stage")=="contours") {
            // The source's {status: count} dict is in first-seen order of the chart results.
            std::vector<std::string> order;
            for(const auto& r:rep.value("contours",Json::array())) {
                const auto status=r.at("status").get<std::string>();if(std::find(order.begin(),order.end(),status)==order.end())order.push_back(status);
            }
            std::string line=baseline+"contours "+st.at("status").get<std::string>()+": ";
            for(std::size_t i=0;i<order.size();++i)line+=(i?", ":"")+integer_text(st.at("charts").at(order[i]))+" "+order[i];
            parts.push_back(line);
        }
    }
    if(rep.contains("final") && rep.at("final").is_object() && !rep.at("final").empty()) {
        const auto& f=rep.at("final");
        parts.push_back("P95 anisotropy "+fixed(number(f.at("anisotropy_p95")),2)+", area "+fixed(number(f.at("area_p95")),2)+", flips "+
                        integer_text(f.at("flipped"))+", crossings "+integer_text(f.at("crossings")));
    }
    parts.push_back(repairs_message(rep));
    // T41 (enhanced_v1): the outcome of the proposal pass; the other policies have no pass.
    if(rep.contains("enhanced") && rep.at("enhanced").is_object()) {
        const auto& e=rep.at("enhanced");const auto status=text_of(e,"status");
        if(status=="done" && member(e,"baseline_retained")==true)parts.push_back("enhanced: native baseline kept");
        else if(status=="done") {
            std::int64_t commits=0;for(const auto& round:member(e,"rounds"))commits+=member(round,"committed")==true;
            const auto islands=count_of(member(e,"retained"),"charts"); // GB M4: the result's islands
            parts.push_back("enhanced: "+std::to_string(commits)+" validated "+plural("proposal",commits)+" applied, "+
                            std::to_string(islands)+" "+plural("island",islands));
        } else if(status=="not_run")parts.push_back("enhanced: not run");
        else if(status=="failed")parts.push_back("enhanced: failed, native baseline kept"); // GB m4
        // GB B1/N2: strict-only defects of the native baseline map that the result keeps (reported).
        if(status=="done") {
            const auto overlap_charts=count_of(member(member(e,"strict_overlap"),"retained"),"charts");
            if(text_of(member(e,"strict_overlap"),"status")=="unrepaired" && overlap_charts>0)
                parts.push_back("strict-only UV overlaps kept in "+std::to_string(overlap_charts)+" "+plural("island",overlap_charts)+" (reported)");
            const auto moved=count_of(member(member(e,"pin_check"),"retained"),"moved");
            if(text_of(member(e,"pin_check"),"status")=="unrepaired" && moved>0)
                parts.push_back(std::to_string(moved)+" pinned "+plural("UV",moved)+" moved by float32 steps kept (reported)");
        }
    }
    // T33 (enhanced_v1): measured seam ends; the other policies have no classification.
    if(rep.contains("enhanced") && rep.at("enhanced").is_object() && rep.at("enhanced").contains("seam_necessity") &&
       rep.at("enhanced").at("seam_necessity").is_object()) {
        const auto& n=rep.at("enhanced").at("seam_necessity");
        const auto necessary=n.value("necessary_relief_ends",std::int64_t{0}),removable=n.value("removable_spur_ends",std::int64_t{0});
        const auto kept=n.value("protected_ends",std::int64_t{0}),unmeasured=n.value("unmeasured_ends",std::int64_t{0});
        if(necessary || removable || kept || unmeasured) {
            std::string line="seam ends: "+std::to_string(necessary)+" necessary relief "+plural("end",necessary);
            if(kept)line+=", "+std::to_string(kept)+" held by constraints";
            if(removable)line+=", "+std::to_string(removable)+" removable";
            if(unmeasured)line+=", "+std::to_string(unmeasured)+" unmeasured";
            parts.push_back(line);
        }
    }
    double total=0;if(rep.contains("timings") && rep.at("timings").contains("total"))total=rep.at("timings").at("total").get<double>();
    parts.push_back(fixed(total,1)+"s");
    if(rep.value("nonuniform_scale",false))
        parts.push_back("object has non-uniform scale: measured in object space, apply the scale for exact texel density");
    std::string out;for(std::size_t i=0;i<parts.size();++i)out+=(i?"; ":"")+parts[i];
    return out;
}
nlohmann::json run_diagnostics(const nlohmann::json& rep) {
    Json out=Json::array();
    auto add=[&](const char* level,const char* topic,std::string text) {out.push_back(diagnostic(level,topic,std::move(text)));};
    const auto policy=text_of(rep,"policy");
    const auto& enh=member(rep,"enhanced");
    const auto status=text_of(enh,"status");
    const auto& retained=member(enh,"retained");
    const bool done=policy=="enhanced_v1" && status=="done";
    // Policy and outcome.
    if(policy=="enhanced_v1") {
        if(status=="done") {
            std::vector<std::string> committed;std::int64_t commits=0;
            for(const auto& round:member(enh,"rounds"))
                if(member(round,"committed")==true) {
                    ++commits;const auto name=text_of(round,"generator","?");
                    if(std::find(committed.begin(),committed.end(),name)==committed.end())committed.push_back(name);
                }
            const auto& counters=member(enh,"counters");
            const auto proposals=count_of(counters,"proposals"),solves=count_of(counters,"solves");
            const auto work=std::to_string(proposals)+" "+plural("proposal",proposals)+", "+std::to_string(solves)+" "+plural("solve",solves);
            if(member(enh,"baseline_retained")==true)
                add("info","outcome","Enhanced: kept the native baseline; no proposal ranked better ("+work+")");
            else
                add("info","outcome","Enhanced: applied "+std::to_string(commits)+" validated "+plural("proposal",commits)+" ("+joined(committed)+"; "+work+")");
        } else if(status=="not_run") {
            const auto base=text_of(enh,"baseline_reason");
            add("warning","outcome",std::string("Enhanced pass did not run: the native baseline ")+
                (text_of(enh,"reason")=="baseline_refused"?"was refused":"did not finish")+(base.empty()?"":" ("+base+")"));
        } else if(status=="refused") {
            add("error","outcome","Enhanced pass refused before any stage: the painted constraints conflict");
        } else if(status=="failed") {
            // GB m4: an error of the pass itself keeps native_baseline's decision and result.
            const auto& failure=member(enh,"failure");
            add("warning","outcome","Enhanced pass failed ("+text_of(failure,"code","error")+": "+text_of(failure,"message")+
                "): kept the native baseline's result (report.enhanced.failure)");
        } else add("info","policy","Policy enhanced_v1");
    } else if(policy=="legacy_010") {
        add("info","policy","Policy legacy_010: the frozen 0.1.0 behavior, for comparisons; no proposal pass");
    } else if(!policy.empty()) {
        add("info","policy","Policy "+policy+": the validated baseline layout; the enhanced proposal pass is off");
    }
    // The refusal and what to do about it.
    if(member(rep,"refusal").is_string())add("error","refusal","Kept the previous state: "+text_of(rep,"refusal"));
    // GB B1: a kept layout that is admissible with candidate zero's strict-only defects is an
    // accepted outcome, named by the strict_overlap / pins lines below, not an error line.
    if(done && member(enh,"valid")==false && member(enh,"admissible")!=true) {
        const auto reason=text_of(retained,"invalid_reason","invalid_uv");
        const auto& charts=member(member(member(retained,"violations"),"charts"),"charts");
        std::string text="Strict validity failed on the kept layout: "+rejection_text(reason);
        if(charts.is_array() && !charts.empty())text+=" in "+plural("island",static_cast<std::int64_t>(charts.size()))+" "+id_list(charts);
        if(member(enh,"cuts_consistent")==false)text+="; the reported seams differ from the solved cuts";
        else if(member(enh,"baseline_retained")==true && member(retained,"feasible")==true)
            text+="; this is the native baseline map, which the baseline checks accept: choose Native baseline to apply it";
        add("error","validity",text);
    }
    for(const std::string kind:{"required_and_forbidden","locked_and_forbidden"}) {
        const auto& ids=member(member(member(enh,"preferences"),"conflicts"),kind.c_str());
        if(ids.is_array() && !ids.empty())
            add("error","constraints",kind=="required_and_forbidden"?"Edges both required and forbidden: "+id_list(ids,12)+"; remove one of the two marks":
                "Locked seams painted forbidden (protected existing seams are locked too): "+id_list(ids,12)+"; remove the lock or the forbidden mark");
    }
    for(const char* kind:{"locked_edges","required_edges","forbidden_edges","anchors"}) {
        const auto& ids=member(member(retained,"violations"),kind);
        if(ids.is_array() && !ids.empty())add("error","constraints",std::string("Unmet hard constraint (")+kind+"): "+id_list(ids,12));
    }
    if(done) {
        // Before/after measurements: the seam_objective.v1 terms of candidate zero and of the kept state.
        const auto& zero=member(member(enh,"candidate_zero"),"objective");const auto& kept=member(retained,"objective");
        if(zero.is_object() && kept.is_object())
            add("info","measurements","Seams "+integer_text(member(zero,"seam_edges"))+" -> "+integer_text(member(kept,"seam_edges"))+" edges, length "+
                short_number(member(zero,"seam_length"))+" -> "+short_number(member(kept,"seam_length"))+", exposed "+
                short_number(member(zero,"exposed_length"))+" -> "+short_number(member(kept,"exposed_length"))+", islands "+
                integer_text(member(zero,"charts"))+" -> "+integer_text(member(kept,"charts"))+" (objective "+short_number(member(zero,"J"))+" -> "+
                short_number(member(kept,"J"))+")");
        // A baseline map that fails the strict checks (or the painted constraints) and the
        // proposals that repaired it.
        const auto& zero_state=member(enh,"candidate_zero");
        if(member(enh,"baseline_retained")==false && member(enh,"valid")==true && member(zero_state,"valid")==false) {
            const auto& listed=member(member(member(zero_state,"violations"),"charts"),"charts");
            const auto& charts=listed.is_array()?listed:member(member(zero_state,"violations"),"strict_defects"); // GB N2: strict-only defects
            std::string text="The native baseline map failed the strict checks ("+rejection_text(text_of(zero_state,"invalid_reason","invalid_uv"));
            if(charts.is_array() && !charts.empty())text+=" in "+plural("island",static_cast<std::int64_t>(charts.size()))+" "+id_list(charts);
            add("info","validity",text+"); the applied proposals pass them");
        }
        // GB B1/N2: strict-only defects of the native baseline map the kept layout still carries: an
        // accepted outcome with native_baseline's decision, reported (never an error or a refusal).
        if(const auto& so=member(enh,"strict_overlap");text_of(so,"status")=="unrepaired") {
            const auto& left=member(so,"retained");const auto charts=count_of(left,"charts"),pairs=count_of(left,"pairs");
            add("warning","strict_overlap",std::to_string(charts)+" "+plural("island",charts)+(charts==1?" keeps ":" keep ")+std::to_string(pairs)+" positive-area UV overlap "+
                plural("pair",pairs)+" (area "+short_number(member(left,"overlap_area"))+") that only the strict overlap check ("+text_of(so,"policy","strict")+
                ") finds; the validated repair did not fix them, so the layout keeps the native baseline's decision (report.enhanced.strict_overlap)");
        }
        if(const auto& pc=member(enh,"pin_check");text_of(pc,"status")!="" && text_of(pc,"status")!="none") {
            std::int64_t steps=0;
            for(const auto& d:member(member(pc,"candidate_zero"),"deviations"))
                for(const auto& s:member(d,"float32_steps"))if(s.is_number_integer())steps=std::max<std::int64_t>(steps,std::abs(s.get<std::int64_t>()));
            const auto moved=count_of(member(pc,"candidate_zero"),"moved");
            const std::string what=std::to_string(moved)+" pinned "+plural("UV",moved)+" by up to "+std::to_string(steps)+" float32 "+plural("step",steps)+
                                   " (within the native pin check, not the strict one; the host applies only exact pins)";
            if(text_of(pc,"status")=="unrepaired")
                add("warning","pins","The unwrap moved "+what+"; no validated re-map fixed it, so the layout keeps the native baseline's decision "
                    "(report.enhanced.pin_check)");
            else add("info","pins","The native baseline map moved "+what+"; re-mapped with exact pins (report.enhanced.pin_check)");
        }
        // Budgets: the worst island against its limit, before and after.
        const auto& b0=member(zero_state,"budget");const auto& b1=member(retained,"budget");
        const auto budget=text_of(enh,"budget_status");const auto over=count_of(b1,"over_budget_charts");
        if(budget=="compliant")
            add("info","budget","Every island within its distortion budget (the worst at "+short_number(member(b1,"worst_ratio"))+" of its limit"+
                (member(b0,"compliant")==false?"; the baseline's worst was "+short_number(member(b0,"worst_ratio"))+"x":std::string())+")");
        else if(budget=="improved_over_budget")
            add("warning","budget","Still over budget but improved: the worst island at "+short_number(member(b1,"worst_ratio"))+"x its limit (was "+
                short_number(member(b0,"worst_ratio"))+"x), "+std::to_string(over)+" "+plural("island",over)+" over");
        else if(!budget.empty())
            add("warning","budget","Over budget: the worst island at "+short_number(member(b1,"worst_ratio"))+"x its limit, "+std::to_string(over)+" "+
                plural("island",over)+" over");
        // Why the other proposals were not taken, per reason; stages without an eligible part.
        std::map<std::string,std::int64_t> refused;
        std::vector<std::string> idle,skipped,proposing,failed,cut;
        for(const auto& round:member(enh,"rounds")) {
            const auto name=text_of(round,"generator","?");
            if(member(round,"skipped").is_string()) {skipped.push_back(name);continue;}
            if(member(round,"failed").is_object()) {failed.push_back(name+" ("+text_of(member(round,"failed"),"code","error")+")");continue;} // GB (m4)
            // GB N5: a round the deadline cut short is a time-budget outcome, not "nothing eligible".
            if(!text_of(round,"stopped").empty()) {if(std::find(cut.begin(),cut.end(),name)==cut.end())cut.push_back(name);}
            else if(std::find(idle.begin(),idle.end(),name)==idle.end())idle.push_back(name);
            if(count_of(round,"proposals")>0)proposing.push_back(name);
            const auto& winner=member(member(round,"winner"),"id");
            for(const auto& c:member(round,"candidates")) {
                if(member(c,"candidate")=="incumbent")continue;
                const auto code=text_of(c,"rejection");
                if(!code.empty())++refused[code];
                else if(member(member(c,"significance"),"below_practical_floor")==true)++refused["below_practical_floor"]; // GB (V2)
                else if(member(member(c,"significance"),"ranks_before_incumbent")==false)++refused["within_sampling_uncertainty"]; // GB (M3)
                else if(member(c,"id")!=winner)++refused["ranked_worse"];
            }
        }
        if(!refused.empty()) {
            std::vector<std::pair<std::string,std::int64_t>> order(refused.begin(),refused.end());
            std::stable_sort(order.begin(),order.end(),[](const auto& a,const auto& b) {return a.second>b.second;});
            std::int64_t total=0;std::vector<std::string> parts;
            for(const auto& [code,n]:order) {total+=n;parts.push_back(std::to_string(n)+" "+rejection_text(code)+" ("+code+")");}
            add("info","refusals","Not taken: "+std::to_string(total)+" "+plural("proposal",total)+": "+joined(parts,"; "));
        }
        // Stages that ran and proposed nothing (no eligible part of this mesh). GB N5: approximate
        // symmetry never proposes (report only; its own line says so), and a stage a deadline cut
        // short is named by the stopped line.
        const bool approximate=text_of(member(enh,"symmetry"),"mode")=="approximate";
        std::erase_if(idle,[&](const std::string& name) {
            return std::find(proposing.begin(),proposing.end(),name)!=proposing.end() || std::find(cut.begin(),cut.end(),name)!=cut.end() ||
                   (approximate && name=="symmetry");
        });
        if(!idle.empty())add("info","eligibility","No candidate from: "+joined(idle)+" (nothing eligible on this mesh)");
        std::vector<std::string> disabled;
        for(const auto& step:member(enh,"schedule"))
            if(member(step,"runs")==false && text_of(step,"reason").rfind("disabled",0)==0)disabled.push_back(text_of(step,"feature","?"));
        if(!disabled.empty())add("info","eligibility","Turned off: "+joined(disabled));
        // GB (m4): stages whose errors were isolated; the best validated layout was kept.
        if(!failed.empty())add("warning","failures","Failed and skipped: "+joined(failed)+"; kept the best validated layout (report.enhanced.failures)");
        const auto stopped=text_of(enh,"stopped");
        if(!stopped.empty()) {
            // GB N5: every stage, writer and classification the deadline cut.
            std::vector<std::string> after;
            if(text_of(member(enh,"contour_cleanup"),"skipped")=="deadline")after.push_back("contour cleanup skipped");
            std::int64_t unmeasured=0;
            for(const auto& arc:member(member(enh,"seam_necessity"),"arcs"))
                if(member(arc,"classification")=="unmeasured" && member(arc,"reason")=="deadline")unmeasured+=count_of(arc,"tips");
            if(unmeasured)after.push_back(std::to_string(unmeasured)+" seam "+plural("end",unmeasured)+" unmeasured");
            add("warning","stopped",std::string(stopped=="deadline"?"Time budget reached during the enhanced pass":"Enhanced pass cancelled")+
                ": kept the best validated layout"+(cut.empty()?std::string():"; cut short: "+joined(cut))+
                (skipped.empty()?std::string():"; skipped "+joined(skipped))+(after.empty()?std::string():"; "+joined(after)));
        }
    }
    // Required relief slits and painted constraints.
    if(const auto& n=member(enh,"seam_necessity");n.is_object()) {
        std::vector<std::string> tips;
        for(const auto& arc:member(n,"arcs"))
            if(member(arc,"classification")=="necessary_relief_end" && tips.size()<6)tips.push_back(id_list(member(arc,"tips"),2));
        const auto necessary=count_of(n,"necessary_relief_ends"),kept=count_of(n,"protected_ends");
        const auto removable=count_of(n,"removable_spur_ends"),unmeasured=count_of(n,"unmeasured_ends");
        if(necessary || kept || removable || unmeasured) {
            std::string text="Seam ends: "+std::to_string(necessary)+" required relief "+plural("slit",necessary);
            if(!tips.empty())text+=" (ending at vertices "+joined(tips,"; ")+")";
            if(kept)text+=", "+std::to_string(kept)+" held by locks, pins or painted constraints";
            if(removable)text+=", "+std::to_string(removable)+" removable";
            if(unmeasured)text+=", "+std::to_string(unmeasured)+" unmeasured";
            add("info","seam_ends",text);
        }
    }
    if(const auto& p=member(enh,"preferences");p.is_object()) {
        const auto required=count_of(p,"required_edges"),forbidden=count_of(p,"forbidden_edges");
        const auto allow=count_of(p,"allow_points"),avoid=count_of(p,"avoid_points");
        if(required || forbidden || allow || avoid)
            add("info","painted","Painted: "+std::to_string(required)+" required and "+std::to_string(forbidden)+" forbidden "+plural("edge",required+forbidden)+
                "; seams welcome at "+std::to_string(allow)+" and avoided at "+std::to_string(avoid)+" "+plural("point",allow+avoid));
        for(const auto& s:member(p,"sources")) {
            const auto st=text_of(s,"status");const auto name=text_of(s,"name",text_of(s,"source","?"));
            if(st.rfind("ignored_",0)==0) {
                std::string text="Attribute "+name+" ignored ("+st.substr(8);
                if(member(s,"domain").is_string())text+=": "+text_of(s,"domain")+" "+text_of(s,"data_type");
                if(member(s,"expected").is_object())text+=", expected "+text_of(member(s,"expected"),"domain")+" "+text_of(member(s,"expected"),"data_type");
                add("warning","painted",text+")");
            } else if(count_of(s,"clamped")>0)
                add("info","painted","Attribute "+name+": "+std::to_string(count_of(s,"clamped"))+" "+plural("weight",count_of(s,"clamped"))+" clamped to [0, 1]");
        }
        const auto& inapplicable=member(p,"inapplicable");
        if(inapplicable.is_object())
            for(auto it=inapplicable.begin();it!=inapplicable.end();++it)
                if(it.value().is_array() && !it.value().empty())
                    add("info","painted","Not applicable ("+it.key()+"): "+std::to_string(it.value().size())+" "+plural("edge",static_cast<std::int64_t>(it.value().size())));
    }
    // Solver backend, symmetry and distance policies.
    if(const auto& paths=member(rep,"solver_paths");paths.is_object()) {
        std::string text="Solver: the host's unwrap ("+text_of(paths,"primary","?")+")";
        if(member(rep,"unwrap_calls").is_number())text+=", "+std::to_string(count_of(rep,"unwrap_calls"))+" "+plural("call",count_of(rep,"unwrap_calls"));
        if(count_of(paths,"conformal_fallbacks"))
            text+=", "+std::to_string(count_of(paths,"conformal_fallbacks"))+" conformal "+(count_of(paths,"conformal_fallbacks")==1?"retry":"retries");
        add("info","solver",text);
    } else if(!policy.empty()) add("info","solver","Solver: the native ARAP backend");
    if(done) {
        if(const auto& s=member(enh,"symmetry");s.is_object()) {
            const auto& map=member(s,"map");const auto map_status=text_of(map,"status","?");
            if(text_of(s,"mode")=="strict") {
                std::int64_t paired=0,conflicts=0;
                const auto& pairing=member(s,"pairing");
                if(pairing.is_object())
                    for(auto g=pairing.begin();g!=pairing.end();++g)
                        if(g.value().is_object())
                            for(auto k=g.value().begin();k!=g.value().end();++k) {
                                if(!k.value().is_number())continue;
                                const auto n=k.value().get<std::int64_t>();
                                if(k.key()=="symmetry_paired" || k.key()=="symmetry_mirror_complete" || k.key()=="symmetry_self_mirrored")paired+=n;
                                if(k.key().rfind("symmetry_conflict",0)==0)conflicts+=n;
                            }
                const bool moves=member(s,"moves_enabled")==true;
                add(moves?"info":"warning","symmetry","Symmetry strict: "+map_status+" mirror correspondence ("+text_of(map,"plane_source","?")+" plane); exact pairing "+
                    short_number(member(member(s,"candidate_zero"),"exact_pairing"))+" -> "+short_number(member(member(s,"retained"),"exact_pairing"))+"; "+
                    std::to_string(paired)+" paired, "+std::to_string(conflicts)+" "+plural("conflict",conflicts)+
                    (moves?std::string():"; no paired move is possible without an exact correspondence"));
                // GB m13: moves of the kept layout committed without their mirror half (on_conflict
                // disable_for_proposal evaluated them alone, labelled in their provenance).
                std::int64_t unpaired=0;
                for(const auto& commit:member(enh,"provenance"))
                    unpaired+=text_of(member(member(commit,"provenance"),"symmetry"),"status")=="disabled_for_proposal";
                if(unpaired)
                    add("warning","symmetry",std::to_string(unpaired)+" applied "+plural("proposal",unpaired)+" unpaired (your conflict policy "+
                        text_of(s,"on_conflict","disable_for_proposal")+"): the mirror half conflicted, so the layout is not mirror-symmetric there");
            } else
                add("info","symmetry","Symmetry approximate (report only; no paired moves, no hard edge equality): "+map_status+" correspondence, seam symmetry "+
                    short_number(member(member(s,"retained"),"seam_symmetry"))+", largest mirror deviation "+
                    short_number(member(member(member(s,"retained"),"nearest_surface"),"max_relative"))+" of the diagonal");
        } else add("info","symmetry","Symmetry asymmetric: view-optimal, moves are not paired");
        if(const auto& sc=member(enh,"scale");sc.is_object()) {
            if(text_of(sc,"mode")=="hops")add("info","scale","Distances: edge hops and face counts (the legacy measure)");
            else {
                const auto& given=member(member(sc,"options"),"reference_length");
                const bool explicit_unit=given.is_number() && given.get<double>()>0;
                add("info","scale","Distances: physical, unit "+short_number(explicit_unit?given:member(member(sc,"reference"),"median_edge_length"))+
                    (explicit_unit?" (the given reference length)":" (the median edge length, frozen at the start)"));
            }
        }
    }
    // Time budget and host facts.
    // GB m15: a deadline met only in the enhanced pass is named by its own "stopped" line above.
    if(member(rep,"timed_out")==true && !enhanced_deadline_only(rep))add("warning","time_budget","Time budget exhausted: a validated partial layout");
    if(member(rep,"nonuniform_scale")==true)
        add("warning","object_scale","Object has non-uniform scale: measured in object space, apply the scale for exact texel density");
    return out;
}
nlohmann::json packing_diagnostics(const nlohmann::json&,const std::string& engine) {
    if(engine!="blender")
        throw ContractError(ErrorCode::invalid_option,"packing.engine","only Blender host packing is supported");
    return Json::array({diagnostic("warning","packing",
        "Packed with Blender's packer: not verified by the engine (overlaps, texel gaps and pins are not checked)")});
}
Json finite_json(const Json& j) {
    if(j.is_object()) {Json out=Json::object();for(auto& [k,v]:j.items())out[k]=finite_json(v);return out;}
    if(j.is_array()) {Json out=Json::array();for(const auto& v:j)out.push_back(finite_json(v));return out;}
    if(j.is_number_float() && !std::isfinite(j.get<double>()))return nullptr;
    return j;
}
SymmetryService native_symmetry_service() {
    SymmetryService service;
    service.detect_plane=[](const MeshSnapshot& mesh,double tol) -> std::optional<DetectedPlane> {
        const auto plane=detect_plane(mesh,tol);if(!plane)return std::nullopt;
        return DetectedPlane{plane->score,plane->to_json(),*plane};
    };
    service.seam_symmetry=[](const MeshSnapshot& mesh,std::span<const std::uint8_t> seams,const DetectedPlane& plane) {
        // A plane detected by another service carries no native payload (Gate A F32).
        const auto* native=std::any_cast<MirrorPlane>(&plane.native);
        if(!native)throw ContractError(ErrorCode::invalid_input,"report.symmetry_plane","the plane was not detected by the native symmetry service");
        return seam_symmetry(mesh,seams,*native).score;
    };
    return service;
}
RunResult run_result(const PipelineResult& res,const PipelineOptions& options,std::span<const CornerPin> pins,const std::string& backend) {
    RunResult out;out.accepted=res.accepted;out.reason=res.reason;
    const auto& rep=res.report;
    if(res.cut) {const auto& disks=res.cut->tables().chart_is_disk;out.topology_valid=std::all_of(disks.begin(),disks.end(),[](auto d) {return d!=0;});}
    if(rep.contains("stored_validation"))out.topology_valid=rep.at("stored_validation").value("topology_valid",false);
    out.uv_valid=rep.contains("final") && rep.at("final").is_object() && rep.at("final").value("valid",false);
    // enhanced_v1 (T28): the strict policy also judges the retained map (overlap, pins).
    if(rep.contains("enhanced") && rep.at("enhanced").is_object() && rep.at("enhanced").contains("valid"))
        out.uv_valid=out.uv_valid && rep.at("enhanced").at("valid").get<bool>();
    // GB (Section 11.3 policy) and GB B1: a strict-only defect of candidate zero that the validated
    // repair could not fix keeps native_baseline's decision. It is an outcome of an accepted run,
    // not an error: it is reported in report.enhanced.strict_overlap / pin_check, `uv_valid` is
    // false, and the report's summary and `messages` name it. RunResult::diagnostics stays the
    // error channel (ErrorCode values; hosts treat an entry there as a failed run).
    out.budget_satisfied=rep.contains("chart_over_budget") && rep.at("chart_over_budget").empty();
    out.timed_out=rep.value("timed_out",false);
    out.partial=out.timed_out || (rep.contains("repairs_summary") && rep.at("repairs_summary").value("unresolved",0)>0);
    if(res.cut) {
        const auto& given=res.cut->mesh().input().geometry_fingerprint;
        out.geometry_fingerprint=given.empty()?res.cut->mesh().geometry_key():given;
        const auto& kind=res.cut->mesh().geometry().face_degenerate_kind;
        for(std::size_t f=0;f<res.exempt.size();++f)if(res.exempt[f])out.exempt_faces.push_back({static_cast<TransferIndex>(f),kind[f]?"degenerate_input":"exempt"});
    }
    if(res.cut) {
        const auto& t=res.cut->tables();
        out.chart_of_face.assign(t.chart_of_face.begin(),t.chart_of_face.end());
        out.chart_of_corner.assign(t.chart_of_corner.begin(),t.chart_of_corner.end());
        for(const auto& edge:res.cut->mesh().tables().edge_vertices)out.edge_vertices.push_back({edge[0],edge[1]});
    }
    out.seams=res.seams;
    out.changed_faces.assign(res.changed_faces.begin(),res.changed_faces.end());
    if(!res.uv.empty())out.uv=res.uv;
    out.pins.assign(pins.begin(),pins.end());
    for(auto ch:res.invalid_charts)out.invalid_charts.push_back(ch);
    out.resolved_options=options;out.backend=backend;
    out.report=finite_json(rep);out.report["summary"]=summarize(res);
    if(!res.scheduling.empty())out.report["scheduling"]=res.scheduling; // T24: batches and workers
    out.report["policy"]=options.policy; // which behavior produced this result (legacy_010 or the corrected product policy)
    if(!res.accepted)out.report["refusal"]=refusal_message(res.reason,res.invalid_charts);
    // T41: enhanced_v1 results carry their actionable lines (CLI and library callers); the
    // reports of the other policies are unchanged.
    if(out.report.contains("enhanced"))out.report["messages"]=run_diagnostics(out.report);
    return out;
}
}
