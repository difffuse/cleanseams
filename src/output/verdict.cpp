// SPDX-License-Identifier: GPL-3.0-or-later
// Source core/verdict.py (finite, measurements, rejected, verdict, compare, aggregate).
#include "cleanseams/verdict.hpp"
#include "cleanseams/errors.hpp"
#include "cleanseams/detail/numpy_compat.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <map>
namespace cleanseams::verdicts {
namespace {
using Json=nlohmann::json;
[[noreturn]] void python_error(const char* what) {throw ContractError(ErrorCode::invalid_input,"verdict",what);}
// Python truthiness of a JSON value.
bool truthy(const Json& v) {
    if(v.is_null())return false;
    if(v.is_boolean())return v.get<bool>();
    if(v.is_number_float())return v.get<double>()!=0.0 || std::isnan(v.get<double>());
    if(v.is_number())return v.get<std::int64_t>()!=0;
    if(v.is_string())return !v.get_ref<const std::string&>().empty();
    return !v.empty();
}
// dict.get(key[, default]) on a value that must be a dict.
Json get(const Json& d,const std::string& key,const Json& fallback=nullptr) {
    if(!d.is_object())python_error("'get' on a value that is not a dict");
    const auto it=d.find(key);return it==d.end()?fallback:*it;
}
Json dict_or_empty(const Json& v) {
    if(!truthy(v))return Json::object();
    if(!v.is_object())python_error("dict() of a value that is not a mapping");
    return v;
}
// Python numbers (bool is an int); arithmetic keeps int-ness.
struct Number {bool is_int=true;std::int64_t i=0;double d=0;double value() const {return is_int?static_cast<double>(i):d;}};
Number number(const Json& v) {
    if(v.is_boolean())return {true,v.get<bool>()?1:0,0};
    if(v.is_number_float())return {false,0,v.get<double>()};
    if(v.is_number())return {true,v.get<std::int64_t>(),0};
    python_error("unsupported operand type for a number");
}
Json to_json(const Number& n) {return n.is_int?Json(n.i):Json(n.d);}
Number sub(const Number& a,const Number& b) {return a.is_int && b.is_int?Number{true,a.i-b.i,0}:Number{false,0,a.value()-b.value()};}
Number mul(const Number& a,const Number& b) {return a.is_int && b.is_int?Number{true,a.i*b.i,0}:Number{false,0,a.value()*b.value()};}
Number absolute(const Number& a) {return a.is_int?Number{true,a.i<0?-a.i:a.i,0}:Number{false,0,std::abs(a.d)};}
bool less(const Number& a,const Number& b) {
    if(a.is_int && b.is_int)return a.i<b.i;
    return a.value()<b.value();
}
bool greater(const Json& a,const Json& b) {return less(number(b),number(a));}
bool less_equal(const Json& a,const Json& b) {const auto x=number(a),y=number(b);return !(less(y,x)) && !(std::isnan(x.value()) || std::isnan(y.value()));}
std::string python_string_repr(const std::string& s) {
    const bool single=s.find('\'')!=std::string::npos && s.find('"')==std::string::npos;
    const char quote=single?'"':'\'';std::string out(1,quote);
    for(const char ch:s) {const auto c=static_cast<unsigned char>(ch);
        if(c=='\\')out+="\\\\";
        else if(ch==quote) {out+='\\';out+=ch;}
        else if(c=='\n')out+="\\n";
        else if(c=='\r')out+="\\r";
        else if(c=='\t')out+="\\t";
        else if(c<0x20 || c==0x7f) {char buffer[8];std::snprintf(buffer,sizeof buffer,"\\x%02x",c);out+=buffer;}
        else out+=static_cast<char>(c);
    }
    return out+quote;
}
const std::map<std::string,double> cpu_budgets{{"T1",10.0},{"T2",60.0},{"T3",300.0},{"T4",900.0}};
const std::array<const char*,6> quality_metrics{"islands","chart_worst_aniso_p95","seam_length_norm","exposure_ratio",
                                                "lane_changes_per_100_edges","packing_efficiency"};
struct Comparison {const char* key;int direction;double tolerance;bool relative;};
const std::array<Comparison,9> comparisons{{{"islands",1,0.0,false},{"chart_worst_aniso_p95",1,0.02,false},{"chart_worst_area_p95",1,0.02,false},
    {"seam_length_norm",1,0.05,true},{"exposure_ratio",1,0.05,false},{"symmetry",-1,0.05,false},{"lane_changes_per_100_edges",1,0.5,false},
    {"packing_efficiency",-1,0.03,false},{"cpu_s",1,0.20,true}}};
}
std::string python_repr(const Json& v) {
    if(v.is_null())return "None";
    if(v.is_boolean())return v.get<bool>()?"True":"False";
    if(v.is_number_float())return detail::python_float_repr(v.get<double>());
    if(v.is_number_unsigned())return std::to_string(v.get<std::uint64_t>());
    if(v.is_number())return std::to_string(v.get<std::int64_t>());
    if(v.is_string())return python_string_repr(v.get_ref<const std::string&>());
    std::string out;
    if(v.is_array()) {
        out="[";bool first=true;for(const auto& x:v) {if(!first)out+=", ";out+=python_repr(x);first=false;}return out+"]";
    }
    out="{";bool first=true;
    for(const auto& [k,x]:v.items()) {if(!first)out+=", ";out+=python_string_repr(k)+": "+python_repr(x);first=false;}
    return out+"}";
}
std::string python_str(const Json& v) {return v.is_string()?v.get<std::string>():python_repr(v);}
bool finite(const Json& v) {
    if(v.is_number_float())return std::isfinite(v.get<double>());
    return v.is_number();
}
Json measurements(const Json& run) {
    const auto report=dict_or_empty(get(run,"report"));
    auto m=dict_or_empty(get(report,"final"));
    const auto metrics=dict_or_empty(get(run,"metrics"));
    for(const auto& [k,v]:metrics.items())m[k]=v;
    const auto s=dict_or_empty(get(run,"seams"));
    auto out=m;for(const auto& [k,v]:s.items())out[k]=v;
    out["islands"]=m.contains("charts")?m.at("charts"):get(s,"islands",get(run,"islands"));
    if(!out.contains("non_disk_charts") && finite(get(m,"charts")) && finite(get(m,"disks")))
        out["non_disk_charts"]=to_json(sub(number(m.at("charts")),number(m.at("disks"))));
    for(const char* key:{"cpu_s","quads_pct"})out[key]=get(run,key);
    // T33 (enhanced_v1): classified dead ends from the run's measured necessity evidence. The
    // legacy raw count stays in dangling_seam_ends; a run without evidence is unchanged.
    const auto necessity=dict_or_empty(get(dict_or_empty(get(report,"enhanced")),"seam_necessity"));
    if(!necessity.empty()) {
        for(const auto& [from,to]:std::initializer_list<std::pair<const char*,const char*>>{{"necessary_relief_ends","necessary_relief_ends"},
                {"removable_spur_ends","removable_spur_ends"},{"protected_ends","protected_seam_ends"},{"unmeasured_ends","unmeasured_seam_ends"}})
            out[to]=get(necessity,from);
        // The evidence describes these seams only when its raw count is the seam report's.
        const auto raw=get(out,"dangling_seam_ends"),evidence=get(necessity,"raw_dangling_ends");
        const bool consistent=finite(raw) && finite(evidence) && !less(number(raw),number(evidence)) && !less(number(evidence),number(raw)) &&
                              finite(out["necessary_relief_ends"]) && finite(out["protected_seam_ends"]);
        out["seam_necessity_consistent"]=consistent;
        if(consistent)out["dangling_seam_ends_unexplained"]=to_json(sub(sub(number(raw),number(out["necessary_relief_ends"])),number(out["protected_seam_ends"])));
    }
    return out;
}
bool rejected(const Json& run) {
    // The stage list is built first, as in the source.
    Json stages=Json::array();
    for(const auto& list:{get(run,"stages"),get(dict_or_empty(get(run,"report")),"stages")}) {
        if(!truthy(list))continue;
        if(!list.is_array())python_error("stages must be a list");
        for(const auto& s:list)stages.push_back(s);
    }
    if(!truthy(get(run,"accepted",false)))return true;
    const auto status=get(run,"status");
    if(status.is_string()) {
        const auto& s=status.get_ref<const std::string&>();
        if(s=="timeout" || s=="error" || s=="unavailable" || s=="rejected")return true;
    }
    if(truthy(get(run,"timed_out",false)))return true;
    for(const auto& s:stages)if(get(s,"status")=="timed_out")return true;
    return false;
}
Json verdict(const Json& run,const Json& meta_in) {
    auto meta=run;if(!meta.is_object())python_error("run must be a dict");
    const auto extra=dict_or_empty(meta_in);
    for(const auto& [k,v]:extra.items())meta[k]=v;
    const auto values=measurements(run);
    Json gates=Json::object(),reasons=Json::array();
    const std::map<std::string,int> severity{{"ship",0},{"fix",1},{"reject",2}};
    std::string result="ship";
    auto gate=[&](const std::string& name,const Json& value,const Json& limit,bool failed,std::string level="fix",bool required=false) {
        const bool missing=value.is_null() || (value.is_number_float() && !std::isfinite(value.get<double>()));
        const char* status=missing?"missing":failed?"fail":"pass";
        gates[name]={{"value",missing?Json(nullptr):value},{"threshold",limit},{"status",status}};
        if((missing && required) || (!missing && failed)) {
            if(missing)level="fix";
            const auto it=severity.find(level);if(it==severity.end())python_error("unknown verdict level");
            if(it->second>severity.at(result))result=level;
            reasons.push_back(name+"="+python_repr(gates[name]["value"])+"; required "+python_str(limit));
        }
    };
    const bool accepted=truthy(get(run,"accepted"));
    gate("accepted",accepted,true,!accepted,"reject");
    gate("completed",!rejected(run),true,rejected(run),"reject");
    for(const char* key:{"flipped","crossings","nonfinite","collapsed_uv","non_disk_charts"}) {
        const auto value=get(values,key);gate(key,value,"== 0",finite(value) && greater(value,0),"reject",true);
    }
    const auto n=get(values,"islands");const auto band=get(meta,"expected_islands");
    auto band_size=[&]() {
        if(band.is_array() || band.is_object() || band.is_string())return band.size();
        python_error("object has no len()");
    };
    auto band_at=[&](std::size_t i) {
        if(!band.is_array())python_error("expected_islands must be a list");
        if(i>=band.size())python_error("list index out of range");
        return band[i];
    };
    if(truthy(band) && band_size()==2) {
        const auto lo=band_at(0),hi=band_at(1);
        const bool outside=finite(n) && !(less_equal(lo,n) && less_equal(n,hi));
        gate("island_band",n,python_str(lo)+" <= islands <= "+python_str(hi),outside,"fix",true);
        const auto triple=to_json(mul(Number{true,3,0},number(hi)));
        gate("island_explosion",n,"<= "+python_str(triple),finite(n) && greater(n,triple),"reject");
    } else gate("expected_islands",nullptr,"a measured asset's expected island band",false,"fix",true);
    const auto tiny=get(values,"tiny_islands_under_0.2pct");const bool tiny_bad=finite(tiny) && greater(tiny,0);
    const bool explosive=truthy(band) && finite(n) && greater(n,band_at(1));
    gate("tiny_islands_under_0.2pct",tiny,"== 0",tiny_bad,explosive?"reject":"fix",true);
    const std::array<std::pair<const char*,Json>,5> caps{{{"chart_worst_aniso_p95",1.6},{"chart_worst_area_p95",1.8},{"exposure_ratio",1.0},
        {"dangling_seam_ends",0},{"hotspot_max",get(meta,"hotspot_gate",4.0)}}};
    // T33: with measured necessity evidence of these seams, only removable or unmeasured dead
    // ends are defects; the raw count stays a measurement. Without it the legacy gate applies.
    const bool classified=get(values,"seam_necessity_consistent",false)==Json(true);
    for(const auto& [key,cap]:caps) {
        auto v=get(values,key);std::string limit="<= "+python_str(cap);
        if(classified && std::string(key)=="dangling_seam_ends") {
            v=get(values,"dangling_seam_ends_unexplained");
            limit+=" removable or unmeasured dead ends (measured necessary relief and protected ends are not defects)";
        }
        gate(key,v,limit,finite(v) && greater(v,cap),"fix",true);
    }
    if(classified)gate("seam_necessity",get(values,"necessary_relief_ends"),"measured necessary relief ends (informational)",false,"fix",false);
    const auto density=get(values,"density_spread");
    gate("density_spread",density,"<= 10",finite(density) && greater(density,10),"reject",true);
    const auto quads=get(values,"quads_pct");
    if(finite(quads) && !less(number(quads),Number{true,90,0})) {
        const auto lane=get(values,"lane_changes_per_100_edges");
        gate("lane_changes_per_100_edges",lane,"<= 1.0",finite(lane) && greater(lane,1),"fix",true);
    }
    if(truthy(get(meta,"symmetry")) || truthy(get(meta,"symmetry_plane")) || truthy(get(meta,"symmetry_axes"))) {
        const auto sym=get(values,"symmetry");
        gate("symmetry",sym,">= 0.8",finite(sym) && less(number(sym),Number{false,0,0.8}),"fix",true);
    }
    const auto tier=get(meta,"tier");
    if(tier.is_array() || tier.is_object())python_error("unhashable type");
    if(tier.is_string()) {
        const auto it=cpu_budgets.find(tier.get<std::string>());
        if(it!=cpu_budgets.end()) {
            const auto cpu=truthy(get(run,"official"))?get(values,"cpu_s"):Json(nullptr);
            gate("cpu_s",cpu,"<= "+detail::python_float_repr(it->second)+" (official timing)",finite(cpu) && greater(cpu,it->second));
        }
    }
    const auto visual=dict_or_empty(get(run,"visual_review"));
    const auto verdict_value=get(visual,"verdict");
    if(verdict_value=="fix" || verdict_value=="reject") {
        gate("visual_review",verdict_value,"no unresolved visible defect",true,verdict_value.get<std::string>());
        reasons.back()=reasons.back().get<std::string>()+": "+python_str(get(visual,"reason","see visual review"));
    }
    return {{"verdict",result},{"reasons",reasons},{"gates",gates}};
}
Json compare(const Json& a,const Json& b,const Json& meta) {
    const auto va=measurements(a),vb=measurements(b);
    const bool bad_a=verdict(a,meta)["verdict"]=="reject",bad_b=verdict(b,meta)["verdict"]=="reject";
    Json band=get(dict_or_empty(meta),"expected_islands");
    if(!truthy(band))band=get(a,"expected_islands");
    if(!truthy(band))band=get(b,"expected_islands");
    Json out=Json::object();std::int64_t wins=0,ties=0,losses=0;
    for(const auto& c:comparisons) {
        const auto x=get(va,c.key),y=get(vb,c.key);
        const Number limit=c.relative && finite(y)?mul(Number{false,0,c.tolerance},absolute(number(y))):Number{false,0,c.tolerance};
        std::string status;Json delta=nullptr;bool severe=false;
        if(bad_a) {status=bad_b?"tie":"lose";severe=status=="lose";}
        else if(!finite(x) || !finite(y) || (std::string(c.key)=="cpu_s" && !(truthy(get(a,"official")) && truthy(get(b,"official")))))
            status="unmeasured";
        else if(bad_b)status="win";
        else {
            Number d=mul(Number{true,c.direction,0},sub(number(x),number(y)));
            if(std::string(c.key)=="islands" && truthy(band)) {
                auto at=[&](std::size_t i) {
                    if(!band.is_array())python_error("expected_islands must be a list");
                    if(i>=band.size())python_error("list index out of range");
                    return number(band[i]);
                };
                // Python max() keeps the first of equal maxima (and its type).
                auto pymax=[](std::initializer_list<Number> values) {
                    Number best=*values.begin();for(const auto& v:values)if(less(best,v))best=v;return best;
                };
                const auto nx=number(x),ny=number(y);
                const auto dx=pymax({sub(at(0),nx),Number{true,0,0},sub(nx,at(1))});
                const auto dy=pymax({sub(at(0),ny),Number{true,0,0},sub(ny,at(1))});
                if(dx.value()!=0 || dy.value()!=0)d=sub(dx,dy);
            }
            const double dv=d.value(),lv=limit.value();
            status=dv>lv+1e-12?"lose":dv<-lv-1e-12?"win":"tie";
            severe=dv>2*lv+1e-12;delta=to_json(d);
        }
        wins+=status=="win";ties+=status=="tie";losses+=status=="lose";
        out[c.key]={{"a",x},{"b",y},{"delta",delta},{"tolerance",to_json(limit)},{"result",status},{"severe_loss",severe}};
    }
    return {{"metrics",out},{"rejected",bad_a},{"reference_rejected",bad_b},{"wins",wins},{"ties",ties},{"losses",losses}};
}
Json aggregate(const Json& class_runs) {
    if(!class_runs.is_array())python_error("class_runs must be a list");
    std::int64_t winners=0,rejections=0,severe=0;bool complete=true;
    for(const auto& run:class_runs) {
        const auto& metrics=run.at("metrics");
        // The source's all(), `measured and sum(...)` and any() short-circuit: nothing after
        // the first unmeasured metric (or the first severe loss) is read (Gate A F32).
        bool measured=true;for(const char* k:quality_metrics)if(metrics.at(k).at("result")=="unmeasured") {measured=false;break;}
        complete=complete && measured;
        if(measured) {std::int64_t won=0;for(const char* k:quality_metrics)won+=metrics.at(k).at("result")=="win";winners+=won>=4;}
        rejections+=truthy(run.at("rejected"));
        bool any=false;for(const auto& [k,v]:metrics.items())if(truthy(v.at("severe_loss"))) {any=true;break;}
        severe+=any;
    }
    const auto count=static_cast<std::int64_t>(class_runs.size());
    const double fraction=count?static_cast<double>(winners)/static_cast<double>(count):0.0;
    return {{"assets",count},{"winning_assets",winners},{"winning_fraction",fraction},{"rejections",rejections},{"severe_losses",severe},
            {"complete",complete},{"claim_ok",count && complete && fraction>=0.8 && !rejections && !severe}};
}
}
