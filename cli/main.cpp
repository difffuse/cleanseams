// SPDX-License-Identifier: GPL-3.0-or-later
// cleanseams_cli (plan section 25.1). Machine output goes to --output files (or stdout
// for --version/--schemas); diagnostics go to stderr as one JSON error record.
// Exit codes: 0 done, 2 invalid input/options, 3 missing/incompatible backend,
// 4 unsuccessful validation or no accepted result, 5 external process failure,
// 6 cancellation, 7 timeout without an accepted partial result, 8 internal failure.
#include "cleanseams/types.hpp"
#include "cleanseams/build_info.hpp"
#include "cleanseams/session.hpp"
#include "cleanseams/errors.hpp"
#include "cleanseams/cache.hpp"
#include "cleanseams/io.hpp"
#include "cleanseams/pipeline.hpp"
#include "cleanseams/report.hpp"
#include "cleanseams/verdict.hpp"
#include "cleanseams/detail/numpy_compat.hpp"
#include "cs_build_info.hpp"
#include <nlohmann/json.hpp>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <string_view>
namespace {
using Json=nlohmann::json;
namespace fs=std::filesystem;
using namespace cleanseams;
constexpr std::uint64_t json_limit=268435456ULL; // 256 MiB of JSON text per input file
struct Exit {int code;std::string field,message;};
[[noreturn]] void quit(int code,const std::string& field,const std::string& message) {throw Exit{code,field,message};}
// Bounded nesting (parse_untrusted_json): results, requests and responses are untrusted.
Json read_json(const fs::path& path,const char* field) {return parse_untrusted_json(read_text_file(path,json_limit,field),field);}
// A missing or oversized configuration is invalid input (exit 2), not an internal error.
Json read_options(const char* path) {return parse_untrusted_json(read_text_file(path,1048576,"options"),"options");}
void write_json(const fs::path& path,const Json& value) {
    std::ofstream out(path,std::ios::binary|std::ios::trunc);
    if(!out)throw ContractError(ErrorCode::io_error,"output","cannot open "+path.string());
    out<<value.dump(2)<<'\n';
    if(!out)throw ContractError(ErrorCode::io_error,"output","cannot write "+path.string());
}
// "--flag value" pairs and bare switches after the command word.
struct Args {
    std::map<std::string,std::string> values;std::set<std::string> flags;
    const std::string& need(const std::string& key) const {
        const auto found=values.find(key);if(found==values.end())quit(2,key,"missing "+key);return found->second;
    }
    std::optional<std::string> get(const std::string& key) const {const auto f=values.find(key);return f==values.end()?std::nullopt:std::optional(f->second);}
};
Args parse(int argc,char** argv,int first,const std::set<std::string>& switches) {
    Args a;
    for(int i=first;i<argc;++i) {
        const std::string key=argv[i];
        if(key.rfind("--",0)!=0)quit(2,key,"unexpected argument "+key);
        if(switches.contains(key)) {a.flags.insert(key);continue;}
        if(i+1>=argc)quit(2,key,"missing value for "+key);
        if(!a.values.emplace(key,argv[++i]).second)quit(2,key,"repeated "+key);
    }
    return a;
}
// Outputs must not overwrite inputs or each other.
void distinct_outputs(const Args& a,const std::vector<std::string>& outputs,const std::vector<std::string>& inputs) {
    std::set<fs::path> seen;
    for(const auto& key:inputs)if(auto v=a.get(key))seen.insert(fs::weakly_canonical(*v));
    for(const auto& key:outputs)if(auto v=a.get(key)) {
        if(!seen.insert(fs::weakly_canonical(*v)).second)quit(2,key,key+" overlaps another input or output");
    }
}
std::shared_ptr<const MeshSnapshot> snapshot(const CaseData& c) {return std::make_shared<const MeshSnapshot>(c.mesh);}
EdgeMask mask_of(const MeshSnapshot& m,std::span<const EdgeEndpoints> pairs,const char* field) {
    EdgeMask out(m.tables().edge_vertices.size(),0);
    for(const auto& e:pairs) {
        const auto id=m.edge_index(e[0],e[1]);
        if(id<0)quit(2,field,"edge endpoint pair is not a surface edge");
        out[static_cast<std::size_t>(id)]=1;
    }
    return out;
}
PinMap pin_map(std::span<const CornerPin> pins) {PinMap out;for(const auto& p:pins)out[static_cast<MeshIndex>(p.corner)]=p.uv;return out;}
// A positive decimal count for a scheduling flag (bounds are checked by resolve_scheduling).
std::uint64_t count_flag(const Args& a,const std::string& key,std::uint64_t fallback) {
    const auto text=a.get(key);if(!text)return fallback;
    if(text->empty() || text->size()>18 || text->find_first_not_of("0123456789")!=std::string::npos)quit(2,key,key+" must be a positive integer");
    return std::stoull(*text);
}
// --threads (native evaluation threads, 1..64; default: the product default, default_native_threads()),
// --batch-charts and --batch-triangles (T24). Every thread count gives the same results.
SchedulingOptions scheduling_of(const Args& a) {
    SchedulingOptions s;s.threads=count_flag(a,"--threads",default_native_threads());s.batch_max_charts=count_flag(a,"--batch-charts",s.batch_max_charts);
    s.batch_max_triangles=count_flag(a,"--batch-triangles",s.batch_max_triangles);
    try {return resolve_scheduling(s);} catch(const ContractError& error) {quit(2,"--"+error.field,error.what());}
}
Json version_json() {return build_info();}

int command_run(const Args& a) {
    distinct_outputs(a,{"--output","--trace","--counters"},{"--input","--options"});
    // The behavior is --behavior when given, else the options' policy, whose default is the
    // product policy native_baseline; legacy_010 reproduces the frozen source for differential
    // runs and enhanced_v1 runs the Phase B proposal pass after it (T41), with the `enhanced`
    // record of the options (case file or --options) as its settings.
    auto check_behavior=[](const std::string& value,const std::string& field) {
        if(value!="legacy_010" && value!="native_baseline" && value!="enhanced_v1" && value!="repair_v1")
            quit(2,field,"behavior must be legacy_010, native_baseline, enhanced_v1 or repair_v1");
    };
    const auto behavior=a.get("--behavior");
    if(behavior)check_behavior(*behavior,"--behavior");
    if(a.get("--backend").value_or("native_arap")!="native_arap")quit(3,"--backend","only the native_arap backend is available from the CLI");
    const auto scheduling=scheduling_of(a);
    const auto c=read_case(a.need("--input"));
    // --options replaces the case file's own options; either passes whole through the
    // options contract (every field applies, unknown fields are refused).
    Json option_json=Json::object();
    if(const auto path=a.get("--options"))option_json=read_options(path->c_str());
    else if(c.extra.contains("options") && !c.extra.at("options").is_null()) {
        if(!c.extra.at("options").is_object())quit(2,"options","the case file's options must be an object");
        option_json=c.extra.at("options");
    }
    if(!option_json.is_object())quit(2,"--options","options must be a JSON object");
    if(behavior)option_json["policy"]=*behavior; // validated with the policy it will run under
    auto options=pipeline_options_from_json(option_json);
    check_behavior(options.policy,"policy");
    if(a.flags.contains("--no-deadline"))options.deadline_s.reset();
    const auto mesh=snapshot(c);
    PipelineInput in;in.operation=c.operation;in.selected_faces=c.selected_faces;in.mesh=mesh;in.seams=mask_of(*mesh,c.seam_edges,"seam_edges");
    const auto locked=mask_of(*mesh,c.locked_edges,"locked_edges");
    if(std::any_of(locked.begin(),locked.end(),[](auto x) {return x!=0;}))in.locked=locked;
    in.pins=pin_map(c.pins);in.anchors=c.anchors;in.uv=c.uv;in.options=options;in.scheduling=scheduling;
    // T41: painted constraints of the case (the run input's `constraints` record: typed masks
    // and weights in core order, raw attributes), checked as a session checks them. GB M2: only
    // enhanced_v1 reads them; under the other policies the record is neither parsed nor checked
    // (it stays an ignored extra key of the case, as at Gate A).
    if(options.policy=="enhanced_v1" && c.extra.contains("constraints") && !c.extra.at("constraints").is_null())
        in.constraints=constraint_input_from_json(c.extra.at("constraints"),mesh->tables().edge_vertices.size(),c.mesh.positions.size());
    in.memory_budget=active_memory_budget(); // recommended_memory_budget(), set in main (Gate A F3)
    Json trace=Json::array();
    const ProgressSink sink=[&](const Json& event) {trace.push_back(finite_json(event));};
    reset_work_counters();const auto started=std::chrono::steady_clock::now();
    const auto res=run_pipeline(in,nullptr,sink);
    const auto pipeline_s=std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count();
    const auto result=run_result(res,options,c.pins,"superlu");
    write_json(a.need("--output"),run_result_json(result));
    if(const auto path=a.get("--trace"))write_json(*path,trace);
    if(const auto path=a.get("--counters")) {
        // Measurement record (T27 matrix): repeated-work counters and the pipeline's own
        // wall time next to the report's stage timings; never part of the result.
        const auto w=work_counters();
        write_json(*path,{{"pipeline_wall_s",pipeline_s},{"timings",finite_json(res.report.value("timings",Json::object()))},
                          {"cpu_timings",finite_json(res.report.value("cpu_timings",Json::object()))},
                          {"work",{{"cut_builds",w.cut_builds},{"unchanged_topology_rebuilds",w.unchanged_topology_rebuilds},{"geometry_builds",w.geometry_builds},
                                   {"duplicate_geometry_builds",w.duplicate_geometry_builds},{"cache_hits",w.cache_hits},{"cache_misses",w.cache_misses},
                                   {"cache_evictions",w.cache_evictions},{"chart_solves",w.chart_solves},{"solve_cache_hits",w.solve_cache_hits}}}});
    }
    if(res.accepted)return 0;
    return result.timed_out?7:4;
}
// Independent measurement of a saved result: topology, validity and pins are recomputed.
int command_validate(const Args& a) {
    distinct_outputs(a,{"--output"},{"--input","--result"});
    const auto c=read_case(a.need("--input"));const auto result=run_result_from_json(read_json(a.need("--result"),"result"));
    const auto mesh=snapshot(c);const auto ne=mesh->tables().edge_vertices.size(),nf=c.mesh.face_offsets.size()-1;
    if(result.seams.size()!=ne)quit(2,"result.seams","seam mask does not match the mesh edges");
    if(!result.uv || result.uv->size()!=c.mesh.corner_vertices.size())quit(4,"result.uv","the result has no complete UV map");
    const auto cs=std::make_shared<const CutSurface>(mesh,result.seams);
    std::vector<std::uint8_t> exempt(nf,0);for(const auto& e:result.exempt_faces)if(e.face>=0 && static_cast<std::size_t>(e.face)<nf)exempt[static_cast<std::size_t>(e.face)]=1;
    const auto metrics=pipeline_metrics(cs,*result.uv,exempt);
    const auto& disks=cs->tables().chart_is_disk;const bool topology=std::all_of(disks.begin(),disks.end(),[](auto d) {return d!=0;});
    bool pins_kept=true;
    for(const auto& p:c.pins) {
        const auto& got=(*result.uv)[static_cast<std::size_t>(p.corner)];
        pins_kept=pins_kept && std::abs(got[0]-p.uv[0])<=1e-8+1e-5*std::abs(p.uv[0]) && std::abs(got[1]-p.uv[1])<=1e-8+1e-5*std::abs(p.uv[1]);
    }
    const bool valid=metrics.at("valid").get<bool>() && topology && pins_kept;
    write_json(a.need("--output"),{{"valid",valid},{"uv_valid",metrics.at("valid")},{"topology_valid",topology},{"pins_kept",pins_kept},
                                   {"charts",cs->tables().n_charts},{"metrics",finite_json(metrics)},{"result_accepted",result.accepted},
                                   {"agrees_with_result",result.accepted==valid}});
    return valid?0:4;
}
int command_convert(const Args& a) {
    distinct_outputs(a,{"--output"},{"--obj"});
    const auto obj=read_obj(a.need("--obj"));
    CaseData c;c.mesh.positions=obj.positions;c.mesh.face_offsets={0};
    for(const auto& f:obj.faces) {for(auto v:f)c.mesh.corner_vertices.push_back(v);c.mesh.face_offsets.push_back(c.mesh.corner_vertices.size());}
    c.uv=obj.uv;
    MeshSnapshot validated(c.mesh); // rejects meshes the engine cannot represent
    write_case(a.need("--output"),c);
    return 0;
}
int command_export(const Args& a) {
    distinct_outputs(a,{"--obj"},{"--input","--result"});
    const auto c=read_case(a.need("--input"));const auto result=run_result_from_json(read_json(a.need("--result"),"result"));
    const auto mesh=snapshot(c);const auto& t=mesh->tables();
    if(result.seams.size()!=t.edge_vertices.size())quit(2,"result.seams","seam mask does not match the mesh edges");
    std::vector<EdgeEndpoints> lines;
    for(std::size_t e=0;e<result.seams.size();++e)if(result.seams[e])lines.push_back({static_cast<TransferIndex>(t.edge_vertices[e][0]),static_cast<TransferIndex>(t.edge_vertices[e][1])});
    const auto faces=case_faces(c.mesh);
    write_obj(a.need("--obj"),c.mesh.positions,faces,result.uv?std::optional<std::span<const Vec2>>(*result.uv):std::nullopt,lines);
    return 0;
}
int command_inspect(const Args& a) {
    distinct_outputs(a,{"--output"},{"--input"});
    const auto c=read_case(a.need("--input"));const auto mesh=snapshot(c);
    const auto cs=std::make_shared<const CutSurface>(mesh,mask_of(*mesh,c.seam_edges,"seam_edges"));
    const auto& disks=cs->tables().chart_is_disk;
    write_json(a.need("--output"),{{"mesh",mesh->summary()},{"diagnose",finite_json(mesh->diagnose())},{"geometry_key",mesh->geometry_key()},
        {"seam_edges",c.seam_edges.size()},{"locked_edges",c.locked_edges.size()},{"pins",c.pins.size()},{"anchors",c.anchors.size()},
        {"uv",c.uv.has_value()},{"charts",cs->tables().n_charts},{"disks",std::count(disks.begin(),disks.end(),std::uint8_t{1})}});
    return 0;
}
// Seam measurements and the evidence verdict of a saved result: the source's
// seam_metrics.seam_report (native mirror-plane detection; exposure over --directions
// view directions, 64 as in the source's benchmark) and verdict.verdict/compare, over the
// run record the source's benchmark built (acceptance, report, seam report, quad share).
// --reference adds the agreement block and a comparison with another result of the same
// case; --meta supplies the verdict's asset facts (expected_islands, symmetry, tier).
int command_report(const Args& a) {
    distinct_outputs(a,{"--output"},{"--input","--result","--reference","--meta"});
    const auto c=read_case(a.need("--input"));const auto mesh=snapshot(c);
    const auto ne=mesh->tables().edge_vertices.size(),nc=c.mesh.corner_vertices.size(),nf=c.mesh.face_offsets.size()-1;
    auto load=[&](const std::string& key) {
        auto r=run_result_from_json(read_json(a.need(key),"result"));
        if(r.seams.size()!=ne)quit(2,key,"seam mask does not match the mesh edges");
        if(r.uv && r.uv->size()!=nc)quit(2,key,"UV map does not match the mesh corners");
        return r;
    };
    const auto result=load("--result");
    std::optional<RunResult> reference;if(a.get("--reference"))reference=load("--reference");
    const Json meta=a.get("--meta")?read_options(a.get("--meta")->c_str()):Json(nullptr);
    if(!meta.is_null() && !meta.is_object())quit(2,"--meta","verdict metadata must be a JSON object");
    SeamReportOptions options;options.K=count_flag(a,"--directions",64);
    if(options.K<1)quit(2,"--directions","--directions must be a positive integer");
    std::size_t quads=0;for(std::size_t f=0;f<nf;++f)quads+=c.mesh.face_offsets[f+1]-c.mesh.face_offsets[f]==4;
    const double quads_pct=detail::python_round(100.0*static_cast<double>(quads)/static_cast<double>(std::max<std::size_t>(nf,1)),1);
    SeamReportContext context; // exposure field and mirror plane: computed once for both results
    auto measure=[&](const RunResult& r,const std::optional<RunResult>& other) {
        const auto cs=std::make_shared<const CutSurface>(mesh,r.seams);
        const auto ref=other?std::optional<std::span<const std::uint8_t>>(other->seams):std::nullopt;
        const auto uv=r.uv?std::optional<std::span<const Vec2>>(*r.uv):std::nullopt;
        return seam_report(mesh,r.seams,cs,&context,ref,uv,options);
    };
    auto run_of=[&](const RunResult& r,const Json& seams) {
        return Json{{"accepted",r.accepted},{"timed_out",r.timed_out},{"report",r.report},{"seams",seams},{"quads_pct",quads_pct}};
    };
    const auto seams=measure(result,reference);const auto run=run_of(result,seams);
    Json out{{"seam_report",finite_json(seams)},{"verdict",finite_json(verdicts::verdict(run,meta))},{"quads_pct",quads_pct},
             {"exposure_directions",options.K}};
    if(reference) {
        const auto reference_seams=measure(*reference,std::nullopt);
        out["reference_seam_report"]=finite_json(reference_seams);
        out["reference_verdict"]=finite_json(verdicts::verdict(run_of(*reference,reference_seams),meta));
        out["compare"]=finite_json(verdicts::compare(run,run_of(*reference,reference_seams),meta));
    }
    write_json(a.need("--output"),out);
    return 0;
}
// The source's claim rule (verdict.aggregate) over a JSON array of `report --reference`
// comparisons of one asset class.
int command_aggregate(const Args& a) {
    distinct_outputs(a,{"--output"},{"--comparisons"});
    const auto comparisons=read_json(a.need("--comparisons"),"comparisons");
    if(!comparisons.is_array())quit(2,"--comparisons","comparisons must be a JSON array of compare records");
    write_json(a.need("--output"),finite_json(verdicts::aggregate(comparisons)));
    return 0;
}
int usage() {
    std::cerr<<"usage: cleanseams_cli --version --json | inspect | run | validate | report | aggregate | convert | export\n"
               "       cleanseams_cli --roundtrip-empty | --schemas | --options FILE | --ui-options FILE | --ui-pack-options FILE\n";
    return 2;
}
}
int main(int argc,char** argv) {
    // Gate A F3: every command builds its tables within half the physical memory (at least
    // the 2 GiB library default); the frozen add-on had no limit of its own.
    const MemoryBudgetScope budget(recommended_memory_budget());
    try {
        if(argc<2)return usage();
        const std::string_view command=argv[1];
        if(argc==2 && command=="--roundtrip-empty") {
            const std::vector<std::uint64_t> offsets{0};const auto mesh=copy_mesh_buffers({},offsets,{});
            std::cout<<Json{{"positions",mesh.positions},{"face_offsets",mesh.face_offsets},{"corner_vertices",mesh.corner_vertices}}.dump()<<'\n';return 0;
        }
        if(argc==2 && command=="--schemas") {std::cout<<contract_schemas().dump(2)<<'\n';return 0;}
        if(argc==3 && command=="--options") {std::cout<<pipeline_options_json(pipeline_options_from_json(read_options(argv[2]))).dump(2)<<'\n';return 0;}
        if(argc==3 && command=="--ui-options") {std::cout<<pipeline_options_json(resolve_ui_options(ui_settings_from_json(read_options(argv[2])))).dump(2)<<'\n';return 0;}
        if(argc==3 && command=="--ui-pack-options") {std::cout<<pipeline_options_json(resolve_ui_pack_options(ui_settings_from_json(read_options(argv[2])))).dump(2)<<'\n';return 0;}
        if(command=="--version") {std::cout<<version_json().dump(2)<<'\n';return 0;}
        if(command=="run")return command_run(parse(argc,argv,2,{"--no-deadline"}));
        if(command=="validate")return command_validate(parse(argc,argv,2,{}));
        if(command=="report")return command_report(parse(argc,argv,2,{}));
        if(command=="aggregate")return command_aggregate(parse(argc,argv,2,{}));
        if(command=="inspect")return command_inspect(parse(argc,argv,2,{}));
        if(command=="convert")return command_convert(parse(argc,argv,2,{}));
        if(command=="export")return command_export(parse(argc,argv,2,{}));
        return usage();
    } catch(const Exit& exit) {
        std::cerr<<Json{{"error",{{"code",exit.code},{"field",exit.field},{"message",exit.message}}}}.dump()<<'\n';return exit.code;
    } catch(const ContractError& error) {
        const int code=error.code==ErrorCode::cancelled?6:error.code==ErrorCode::host_failure?5:2;
        std::cerr<<Json{{"error",{{"code",error_code_name(error.code)},{"field",error.field},{"message",error.what()}}}}.dump()<<'\n';return code;
    } catch(const nlohmann::json::exception& error) {
        std::cerr<<Json{{"error",{{"code","invalid_input"},{"field","json"},{"message",error.what()}}}}.dump()<<'\n';return 2;
    } catch(const std::filesystem::filesystem_error& error) {
        // Unresolvable user paths (missing directories, permissions) are input errors.
        std::cerr<<Json{{"error",{{"code","io_error"},{"field","path"},{"message",error.what()}}}}.dump()<<'\n';return 2;
    } catch(const std::exception& error) {
        std::cerr<<Json{{"error",{{"code","internal"},{"field",""},{"message",error.what()}}}}.dump()<<'\n';return 8;
    }
}
