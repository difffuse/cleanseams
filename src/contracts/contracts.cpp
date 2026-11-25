// SPDX-License-Identifier: GPL-3.0-or-later
#include "cleanseams/session.hpp"
#include "cleanseams/detail/option_fields.hpp"
#include "cleanseams/detail/contract_fields.hpp"
#include "cleanseams/detail/codec.hpp"
#include "cleanseams/control.hpp"
#include <cmath>
#include <limits>
#include <map>
#include <set>
namespace cleanseams {
namespace {
using detail::Json;
void require(bool condition,const std::string& field,const char* message) {
    if(!condition) throw ContractError(ErrorCode::invalid_input,field,message);
}
void index(TransferIndex value,std::size_t size,const std::string& path) {
    require(value>=0 && static_cast<std::uint64_t>(value)<size,path,"index out of range");
}
void check_mesh(const MeshInput& mesh,const std::string& path) {
    const auto& offsets=mesh.face_offsets;
    require(!offsets.empty() && offsets.front()==0 && offsets.back()==mesh.corner_vertices.size() && std::is_sorted(offsets.begin(),offsets.end()),path+".face_offsets","offsets must span corners monotonically from zero");
    require(mesh.positions.size()<std::numeric_limits<std::uint32_t>::max() && mesh.corner_vertices.size()<std::numeric_limits<std::uint32_t>::max() && offsets.size()<=std::numeric_limits<std::uint32_t>::max(),path,"mesh exceeds internal index capacity");
    for(auto i:mesh.corner_vertices) index(i,mesh.positions.size(),path+".corner_vertices");
    for(const auto& edge:mesh.edge_vertices) for(auto i:edge) index(i,mesh.positions.size(),path+".edge_vertices");
    auto mapping=[&](const auto& values,std::size_t count,const char* field) {
        require(values.empty() || values.size()==count,path+field,"mapping length mismatch");
        std::set<TransferIndex> seen;
        for(auto v:values) {require(v>=0,path+field,"negative host identity");require(seen.insert(v).second,path+field,"duplicate host identity");}
    };
    mapping(mesh.host_loop_indices,mesh.corner_vertices.size(),".host_loop_indices");
    mapping(mesh.host_edge_indices,mesh.edge_vertices.size(),".host_edge_indices");
}
void check_uv(const std::optional<std::vector<Vec2>>& uv,std::size_t count,const std::string& path) {
    require(!uv || uv->size()==count,path,"UV count must equal corner count");
}
void check_pins(const std::vector<CornerPin>& pins,std::size_t count,const std::string& path) {
    std::set<TransferIndex> corners;
    for(const auto& pin:pins) {index(pin.corner,count,path+".corner");require(corners.insert(pin.corner).second,path,"duplicate corner pin");}
}
void check_mask(const std::vector<std::uint8_t>& mask,std::size_t count,const std::string& path) {
    require(mask.empty() || mask.size()==count,path,"mask length mismatch");
    for(auto value:mask) require(value<=1,path,"mask values must be zero or one");
}
void check_input(const RunInput& in) {
    check_mesh(in.mesh,"mesh");const auto nc=in.mesh.corner_vertices.size(),nv=in.mesh.positions.size(),ne=in.mesh.edge_vertices.size();
    check_mask(in.selected_faces,in.mesh.face_offsets.size()-1,"selected_faces");
    check_uv(in.uv,nc,"uv");
    require(in.operation=="pipeline" || in.operation=="postprocess_layout","operation","unknown operation");
    if(in.operation=="postprocess_layout") {
        require(in.uv.has_value(),"uv","postprocess_layout requires supplied UVs");
        require(in.options.policy=="repair_v1","options.policy","postprocess_layout requires repair_v1");
        for(std::size_t f=0;f+1<in.mesh.face_offsets.size();++f)for(auto c=in.mesh.face_offsets[f];c<in.mesh.face_offsets[f+1];++c) {
            const auto& p=(*in.uv)[c];
            const bool recover=in.options.repair.region_recovery_seconds>0 && !in.selected_faces.empty() && in.selected_faces[f];
            require(recover || (std::isfinite(p[0]) && std::isfinite(p[1])),"uv","nonfinite supplied UVs require selected region recovery");
        }
    }
    check_pins(in.pins,nc,"pins");
    check_mask(in.seam_mask,ne,"seam_mask");check_mask(in.lock_mask,ne,"lock_mask");
    require(in.seam_edges.empty() || in.seam_mask.empty(),"seam_edges","cannot combine endpoint pairs and a mask");
    require(in.locked_edges.empty() || in.lock_mask.empty(),"locked_edges","cannot combine endpoint pairs and a mask");
    for(const auto& edges:{in.seam_edges,in.locked_edges}) for(const auto& edge:edges) for(auto i:edge) index(i,nv,"edge_endpoints");
    for(auto i:in.anchors) index(i,nv,"anchors");
    check_mask(in.constraints.required_edges,ne,"constraints.required_edges");
    check_mask(in.constraints.forbidden_edges,ne,"constraints.forbidden_edges");
    require(in.constraints.allow_points.empty() || in.constraints.allow_points.size()==nv,"constraints.allow_points","preference length mismatch");
    require(in.constraints.avoid_points.empty() || in.constraints.avoid_points.size()==nv,"constraints.avoid_points","preference length mismatch");
    pipeline_options_json(in.options);
}
Json identity_schema() {
    auto s=detail::schema_for(HostIdentity{});s["properties"]["version"]["enum"]={1};return s;
}
Json input_schema() {
    auto s=detail::schema_for(RunInput{});
    s["properties"]["version"]["enum"]={1};s["properties"]["coordinate_space"]["enum"]={"object"};
    s["properties"]["operation"]["enum"]={"pipeline","postprocess_layout"};
    auto& required=s["required"];required.erase(std::remove(required.begin(),required.end(),"operation"),required.end());
    required.erase(std::remove(required.begin(),required.end(),"selected_faces"),required.end());
    s["properties"]["options"]=option_registry()["pipeline"];return s;
}
Json result_schema() {
    auto s=detail::schema_for(RunResult{});s["properties"]["version"]["enum"]={1};
    s["properties"]["resolved_options"]=option_registry()["pipeline"];
    s["properties"]["uv"]["items"]["items"]["type"]={"number","null"};
    s["properties"]["report"]["type"]="object";s["properties"]["decision_trace"]["type"]="array";
    s["properties"]["diagnostics"]["items"]["properties"]["code"]["enum"]=error_codes;
    for(const auto* name:{"chart_of_face","chart_of_corner","changed_faces","edge_vertices"}) {
        auto& required=s["required"];required.erase(std::remove(required.begin(),required.end(),name),required.end());
    }
    return s;
}
Json request_schema() {
    auto s=detail::schema_for(HostRequest{});s["properties"]["identity"]=identity_schema();
    s["properties"]["kind"]["enum"]={"solve","pack","region_proposal"};
    s["properties"]["solve"]["properties"]["method"]["enum"]={"MINIMUM_STRETCH","ANGLE_BASED","CONFORMAL"};
    s["properties"]["pack"]["properties"]["options"]=option_registry()["pipeline"]["properties"]["packing"];
    s["properties"]["pack"]["properties"]["backend"]["enum"]={"blender"};
    auto& required=s["required"];required.erase(std::remove(required.begin(),required.end(),"region_proposal"),required.end());
    return s;
}
Json response_schema() {
    auto s=detail::schema_for(HostResponse{});s["properties"]["identity"]=identity_schema();
    s["properties"]["kind"]["enum"]={"solve","pack","region_proposal"};
    s["properties"]["method"]["enum"]={"MINIMUM_STRETCH","ANGLE_BASED","CONFORMAL","SMART_PROJECT","blender"};
    s["properties"]["error"]["properties"]["code"]["enum"]=result_schema()["properties"]["diagnostics"]["items"]["properties"]["code"]["enum"];
    s["properties"]["elapsed_s"]["minimum"]=0.;return s;
}
Json descriptor_schema() {
    const Json unsigned_count={{"type","integer"},{"minimum",0},{"maximum",std::numeric_limits<std::uint64_t>::max()}};
    return {{"type","object"},{"additionalProperties",false},
        {"required",{"file","dtype","shape","byte_offset","byte_length","sha256"}},
        {"properties",{{"file",{{"type","string"},{"minLength",1},{"description","Relative path; T18 additionally rejects absolute paths, traversal and symlink escape."}}},
            {"dtype",{{"type","string"},{"enum",{"float64-le","int64-le","uint64-le","uint8"}}}},
            {"shape",{{"type","array"},{"minItems",1},{"items",unsigned_count}}},
            {"byte_offset",unsigned_count},{"byte_length",unsigned_count},
            {"sha256",{{"type","string"},{"pattern","^[0-9a-f]{64}$"}}}}}};
}
void option_patch_schema(Json& value) {
    if(value.is_object()) {value.erase("required");for(auto it=value.begin();it!=value.end();++it) option_patch_schema(it.value());}
    else if(value.is_array()) for(auto& v:value) option_patch_schema(v);
}
Json case_schema() {
    auto schema=input_schema();schema.erase("default");auto& props=schema["properties"];
    const auto mesh=props["mesh"]["properties"];props.erase("mesh");
    for(auto it=mesh.begin();it!=mesh.end();++it) props[it.key()=="edge_vertices"?"host_edge_vertices":it.key()]=it.value();
    props["schema"]={{"type","string"},{"enum",{"cleanseams.case"}}};
    schema["required"]={"schema","version","coordinate_space","positions","face_offsets","corner_vertices"};
    for(const auto* name:{"positions","face_offsets","corner_vertices","host_edge_vertices","host_edge_indices","host_loop_indices","seam_mask","lock_mask"}) {
        auto descriptor=descriptor_schema();
        descriptor["properties"]["dtype"]["enum"]={std::string(name)=="positions"?"float64-le":std::string(name)=="face_offsets"?"uint64-le":std::string(name)=="seam_mask" || std::string(name)=="lock_mask"?"uint8":"int64-le"};
        props[name]={{"oneOf",Json::array({props[name],descriptor})}};
    }
    auto uv_descriptor=descriptor_schema();uv_descriptor["properties"]["dtype"]["enum"]={"float64-le"};
    props["uv"]={{"oneOf",Json::array({props["uv"],uv_descriptor})}};
    option_patch_schema(props["options"]);
    // Fixture generation metadata is inert and does not enable runtime behavior.
    for(const auto* name:{"id","family"}) props[name]={{"type","string"}};
    props["pipeline_capture"]={{"type","boolean"}};
    for(const auto* name:{"parameters","fields","host_setup","packing","symmetry"}) props[name]=Json::object();
    for(const auto* name:{"host_actions","cache_trials","random_cut_states","reference_direction","scenarios"}) props[name]=Json::object();
    return schema;
}
void check_request(const HostRequest& request) {
    if(request.kind=="solve") {
        require(request.solve.has_value() && !request.pack && !request.region_proposal,"solve","solve request must contain only solve payload");
        const auto& s=*request.solve;check_mesh(s.mesh,"solve.mesh");const auto count=s.mesh.corner_vertices.size();
        require(s.scratch_to_original_corner.size()==count,"solve.scratch_to_original_corner","mapping length mismatch");
        for(auto c:s.scratch_to_original_corner) require(c>=0,"solve.scratch_to_original_corner","negative original corner");
        check_uv(s.initial_uv,count,"solve.initial_uv");check_pins(s.pins,count,"solve.pins");
    } else if(request.kind=="region_proposal") {
        require(request.region_proposal.has_value() && !request.solve && !request.pack,"region_proposal","projection request must contain only region proposal payload");
        const auto& p=*request.region_proposal;check_mesh(p.mesh,"region_proposal.mesh");
        require(p.scratch_to_original_corner.size()==p.mesh.corner_vertices.size(),"region_proposal.scratch_to_original_corner","mapping length mismatch");
        std::set<TransferIndex> seen;
        for(auto c:p.scratch_to_original_corner) {
            require(c>=0,"region_proposal.scratch_to_original_corner","negative original corner");
            require(seen.insert(c).second,"region_proposal.scratch_to_original_corner","duplicate original corner");
        }
    } else {
        require(request.pack.has_value() && !request.solve && !request.region_proposal,"pack","pack request must contain only pack payload");
        const auto& p=*request.pack;check_mesh(p.mesh,"pack.mesh");
        check_uv(p.uv,p.mesh.corner_vertices.size(),"pack.uv");check_pins(p.pins,p.mesh.corner_vertices.size(),"pack.pins");
        PipelineOptions options;options.packing=p.options;pipeline_options_json(options);
    }
}
void check_response(const HostResponse& response) {
    require(response.success ? !response.error : response.error.has_value(),"error","failure requires a structured error; success forbids it");
    if(response.success && (response.kind=="solve" || response.kind=="region_proposal"))
        require(response.uv.has_value(),"uv","successful solve or proposal requires UV output");
    if(response.success && response.kind=="region_proposal")require(response.method=="SMART_PROJECT","method","region proposal must report SMART_PROJECT");
    if(response.kind=="solve")require(response.method!="SMART_PROJECT","method","projection is not a fixed-cut solve");
}
void sanitize(Json& value,const std::string& path,std::vector<UnavailableValue>& unavailable) {
    if(value.is_number_float() && !std::isfinite(value.get<double>())) {
        const auto x=value.get<double>();unavailable.push_back({path,std::isnan(x)?"nan":x>0?"positive_infinity":"negative_infinity"});value=nullptr;
    } else if(value.is_object()) for(auto it=value.begin();it!=value.end();++it) sanitize(it.value(),(Json::json_pointer(path)/it.key()).to_string(),unavailable);
    else if(value.is_array()) for(std::size_t i=0;i<value.size();++i) sanitize(value[i],path+"/"+std::to_string(i),unavailable);
}
void check_result(const RunResult& result) {
    for(auto v:result.seams) require(v<=1,"seams","mask values must be zero or one");
    if(result.uv) check_pins(result.pins,result.uv->size(),"pins");
    if(result.uv && (result.accepted || result.uv_valid))for(const auto& p:*result.uv)
        require(std::isfinite(p[0]) && std::isfinite(p[1]),"uv","an accepted or UV-valid result must have finite coordinates");
    for(const auto& item:result.exempt_faces) require(item.face>=0 && !item.reason.empty(),"exempt_faces","face identity and reason required");
    for(const auto& item:result.unavailable) require(!item.path.empty() && !item.reason.empty(),"unavailable","path and reason required");
    pipeline_options_json(result.resolved_options);
}
template<class T> T parse(const Json& input,const Json& schema) {
    detail::validate_schema(input,schema,"",ErrorCode::invalid_input);
    return detail::decode<T>(input,"");
}
// Phase B records are optional on the wire (T28): an options object without `enhanced`, or
// with a partial one, is read with the record's defaults, so every Gate A document (fuzz
// corpus, stored results) stays valid. Every other field stays required.
Json with_enhanced_defaults(Json document,const char* key) {
    if(document.is_object() && document.contains(key) && document.at(key).is_object()) {
        auto& options=document[key];const auto defaults=detail::encode(EnhancedOptions{});
        const auto repair_defaults=detail::encode(RepairOptions{});
        options["repair"]=options.contains("repair")?detail::merge_defaults(repair_defaults,options.at("repair")):repair_defaults;
        options["enhanced"]=options.contains("enhanced")?detail::merge_defaults(defaults,options.at("enhanced")):defaults;
    }
    return document;
}
// The painted attributes of `constraints` (T29) are optional on the wire as well.
Json with_attribute_defaults(Json document) {
    if(document.is_object()) {
        if(!document.contains("operation"))document["operation"]="pipeline";
        if(!document.contains("selected_faces"))document["selected_faces"]=Json::array();
    }
    if(document.is_object() && document.contains("constraints") && document.at("constraints").is_object() &&
       !document.at("constraints").contains("attributes"))
        document["constraints"]["attributes"]=Json::array();
    return document;
}
}
RunInput run_input_from_json(const Json& input) {
    auto in=parse<RunInput>(with_attribute_defaults(with_enhanced_defaults(input,"options")),input_schema());check_input(in);return in;
}
ConstraintInput constraint_input_from_json(const Json& input,std::uint64_t edges,std::uint64_t vertices) {
    if(!input.is_object())throw ContractError(ErrorCode::invalid_input,"constraints","constraints must be an object");
    const auto complete=detail::merge_defaults(detail::encode(ConstraintInput{}),input);
    detail::validate_schema(complete,input_schema()["properties"]["constraints"],"constraints",ErrorCode::invalid_input);
    auto out=detail::decode<ConstraintInput>(complete,"constraints");
    check_constraint_input(out,edges,vertices);return out;
}
Json run_input_json(const RunInput& in) {auto j=detail::encode(in);detail::validate_schema(j,input_schema(),"",ErrorCode::invalid_input);check_input(in);return j;}
HostRequest host_request_from_json(const Json& input) {
    auto document=input;if(document.is_object() && !document.contains("region_proposal"))document["region_proposal"]=nullptr;
    auto r=parse<HostRequest>(document,request_schema());check_request(r);return r;
}
Json host_request_json(const HostRequest& r) {auto j=detail::encode(r);detail::validate_schema(j,request_schema(),"",ErrorCode::invalid_input);check_request(r);return j;}
HostResponse host_response_from_json(const Json& input) {auto r=parse<HostResponse>(input,response_schema());check_response(r);return r;}
Json host_response_json(const HostResponse& r) {auto j=detail::encode(r);detail::validate_schema(j,response_schema(),"",ErrorCode::invalid_input);check_response(r);return j;}
RunResult run_result_from_json(const Json& input) {auto document=with_enhanced_defaults(input,"resolved_options");
    if(document.is_object())for(const auto* name:{"chart_of_face","chart_of_corner","changed_faces","edge_vertices"})
        if(!document.contains(name))document[name]=Json::array();
    detail::validate_schema(document,result_schema(),"",ErrorCode::invalid_input);
    std::map<std::string,std::string> missing;
    for(const auto& record:document.at("unavailable")) {
        const auto path=record.at("path").get<std::string>();
        if(path=="/uv" || path.starts_with("/uv/"))
            require(missing.emplace(path,record.at("reason").get<std::string>()).second,"unavailable","duplicate UV unavailable path");
    }
    auto& uv=document["uv"];
    if(uv.is_array())for(std::size_t c=0;c<uv.size();++c)for(std::size_t d=0;d<2;++d)if(uv[c][d].is_null()) {
        const std::string path="/uv/"+std::to_string(c)+"/"+std::to_string(d);const auto found=missing.find(path);
        require(found!=missing.end(),"unavailable","null UV coordinate requires an unavailable reason");
        const auto& reason=found->second;
        require(reason=="nan" || reason=="positive_infinity" || reason=="negative_infinity","unavailable","unknown UV unavailable reason");
        uv[c][d]=reason=="nan"?std::numeric_limits<double>::quiet_NaN():reason=="positive_infinity"?INFINITY:-INFINITY;
        missing.erase(found);
    }
    require(missing.empty(),"unavailable","UV unavailable path must identify a null coordinate");
    // Structural/finite JSON validation precedes restoring typed unavailable values.
    auto r=detail::decode<RunResult>(document,"");check_result(r);return r;}
Json run_result_json(const RunResult& result) {
    check_result(result);auto normalized=result;
    // Regenerate UV unavailable records from the owned array, making repeated
    // encode/decode cycles idempotent without trusting stale caller metadata.
    std::erase_if(normalized.unavailable,[](const auto& value){return value.path=="/uv" || value.path.starts_with("/uv/");});
    sanitize(normalized.report,"/report",normalized.unavailable);sanitize(normalized.decision_trace,"/decision_trace",normalized.unavailable);
    auto j=detail::encode(normalized);sanitize(j["uv"],"/uv",normalized.unavailable);
    j["unavailable"]=detail::encode(normalized.unavailable);
    detail::validate_schema(j,result_schema(),"",ErrorCode::invalid_input);return j;
}
SessionConfig resolve_session_config(SessionConfig c) {
    if(c.worker_count==0)c.worker_count=default_native_threads(); // automatic (GB performance)
    require(c.worker_count>0 && c.worker_count<=max_native_threads,"worker_count","must be between 1 and 64, or 0 for automatic");
    require(c.memory_budget_bytes>0,"memory_budget_bytes","must be positive");
    require(c.batch_max_charts>0,"batch_max_charts","must be positive");require(c.batch_max_triangles>0,"batch_max_triangles","must be positive");
    // "host" solves through typed host requests (Blender); the others are native.
    require(c.backend=="superlu" || c.backend=="numpy_compat" || c.backend=="host","backend","unknown backend");
    require(c.trace_level=="off" || c.trace_level=="stages" || c.trace_level=="full","trace_level","unknown trace level");
    std::set<std::string> capabilities;
    for(const auto& capability:c.host_capabilities) {
        require(capability=="region_proposal","host_capabilities","unknown optional host capability");
        require(capabilities.insert(capability).second,"host_capabilities","duplicate optional host capability");
    }
    require(c.fragment_unwrap_method=="ANGLE_BASED" || c.fragment_unwrap_method=="CONFORMAL" || c.fragment_unwrap_method=="MINIMUM_STRETCH",
            "fragment_unwrap_method","expected a fixed-cut host unwrap method");
    require(c.host_unwrap_method=="ANGLE_BASED" || c.host_unwrap_method=="CONFORMAL" || c.host_unwrap_method=="MINIMUM_STRETCH",
            "host_unwrap_method","expected a fixed-cut host unwrap method");
    if(!c.cache_budget_bytes) c.cache_budget_bytes=std::min<std::uint64_t>(536870912ULL,c.memory_budget_bytes/4);
    require(*c.cache_budget_bytes<=c.memory_budget_bytes,"cache_budget_bytes","exceeds session budget");return c;
}
Json contract_schemas() {
    auto session_schema=detail::schema_for(SessionConfig{});
    for(const auto* field:{"host_capabilities","fragment_unwrap_method","host_unwrap_method"}) {
        auto& required=session_schema["required"];required.erase(std::remove(required.begin(),required.end(),field),required.end());
    }
    session_schema["properties"]["host_capabilities"]["items"]["enum"]={"region_proposal"};
    session_schema["properties"]["host_capabilities"]["uniqueItems"]=true;
    session_schema["properties"]["fragment_unwrap_method"]["enum"]={"ANGLE_BASED","CONFORMAL","MINIMUM_STRETCH"};
    session_schema["properties"]["host_unwrap_method"]["enum"]={"ANGLE_BASED","CONFORMAL","MINIMUM_STRETCH"};
    auto schemas=Json{{"options",option_registry()["pipeline"]},{"ui",option_registry()["ui"]},
        {"run_input",input_schema()},{"result",result_schema()},{"host_request",request_schema()},
        {"host_response",response_schema()},{"session",session_schema},{"progress",detail::schema_for(ProgressEvent{})},
        {"array_descriptor",descriptor_schema()},{"case",case_schema()}};
    for(auto it=schemas.begin();it!=schemas.end();++it) {
        it.value()["$schema"]="https://json-schema.org/draft/2020-12/schema";
        it.value()["$id"]="urn:cleanseams:"+it.key()+":1";
        it.value()["description"]="Structural wire contract; native validation also checks indices, lengths and cross-field invariants.";
    }
    return schemas;
}
}
