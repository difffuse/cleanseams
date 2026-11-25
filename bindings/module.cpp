// SPDX-License-Identifier: GPL-3.0-or-later
#include "arrays.hpp"
#include <pybind11/stl.h>
#include "mesh_analysis.hpp"
#include "cleanseams/types.hpp"
#include "cleanseams/build_info.hpp"
#include "cleanseams/options.hpp"
#include "cleanseams/session.hpp"
#include "cleanseams/host_mesh.hpp"
#include "cleanseams/mesh.hpp"
#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstring>
#include <filesystem>
#include <limits>
#include <mutex>
#include <set>
#include "cleanseams/errors.hpp"
#include "cleanseams/report.hpp"
#include "cleanseams/io.hpp"
#include "cleanseams/cut_surface.hpp"
#include "cleanseams/metrics.hpp"
namespace py = pybind11;
using namespace cs_bindings;
namespace {
cleanseams::HostMeshView host_view(const py::array& positions,const py::array& edges,
    const py::array& vertices,const py::array& loop_edges,const py::array& starts,const py::array& totals) {
    require_array<float>(positions,2,"positions");require_array<std::int32_t>(edges,2,"edge_vertices");
    require_array<std::int32_t>(vertices,1,"loop_vertices");require_array<std::int32_t>(loop_edges,1,"loop_edges");
    require_array<std::int32_t>(starts,1,"face_loop_start");require_array<std::int32_t>(totals,1,"face_loop_total");
    if(positions.shape(1)!=3 || edges.shape(1)!=2) throw py::value_error("positions/edge_vertices require shapes (n,3)/(n,2)");
    return {elements<float>(positions),elements<std::int32_t>(edges),elements<std::int32_t>(vertices),
        elements<std::int32_t>(loop_edges),elements<std::int32_t>(starts),elements<std::int32_t>(totals)};
}
py::dict read_host(const py::array& positions,const py::array& edges,const py::array& vertices,
    const py::array& loop_edges,const py::array& starts,const py::array& totals,
    const py::array& matrix,const py::array& scale,bool world_space,std::uint64_t budget,const py::dict& attributes) {
    const auto view=host_view(positions,edges,vertices,loop_edges,starts,totals);
    require_array<double>(matrix,2,"object_to_world");require_array<double>(scale,1,"object_scale");
    if(matrix.shape(0)!=4 || matrix.shape(1)!=4 || scale.size()!=3) throw py::value_error("expected a (4,4) matrix and (3,) scale");
    cleanseams::ObjectTransform transform;transform.world_space=world_space;
    std::copy_n(static_cast<const double*>(matrix.data()),16,transform.object_to_world.begin());
    std::copy_n(static_cast<const double*>(scale.data()),3,transform.object_scale.begin());
    for(const auto& item:attributes) {
        if(!py::isinstance<py::str>(item.first)) throw py::value_error("host attribute names must be strings");
        const auto key=py::cast<std::string>(item.first);
        if(key!="seams" && key!="locked" && key!="anchors" && key!="pins" && key!="uv" && key!="painted") throw py::value_error("unknown host attribute: "+key);
    }
    // Synchronous reads hold the GIL. No Python view or object enters a snapshot.
    const auto result=cleanseams::read_host_mesh(view,transform,budget);
    const cleanseams::MeshSnapshot snapshot(result.mesh,budget);
    auto mask=[&](const char* key)->std::span<const std::uint8_t> {
        if(!attributes.contains(key)) return {};
        const auto value=dict_array(attributes,key);require_array<bool>(value,1,key);
        return elements<std::uint8_t>(value);
    };
    cleanseams::HostAttributes host_attributes{mask("seams"),mask("locked"),mask("anchors"),mask("pins"),std::nullopt};
    std::vector<cleanseams::HostPaintedAttribute> host_painted;
    // Painted attributes (T29): {name, domain, data_type, values (float64, host order)}.
    std::vector<py::array> painted_values;
    if(attributes.contains("painted")) {
        for(const auto& item:py::cast<py::list>(attributes["painted"])) {
            const auto entry=py::cast<py::dict>(item);
            cleanseams::HostPaintedAttribute a;a.name=dict_str(entry,"name","painted");a.domain=dict_str(entry,"domain","painted");
            a.data_type=dict_str(entry,"data_type","painted");
            const auto values=dict_array(entry,"values");require_array<double>(values,1,"painted.values");painted_values.push_back(values);
            a.values=elements<double>(values);host_painted.push_back(std::move(a));
        }
    }
    std::vector<cleanseams::Vec2> host_uv;
    if(attributes.contains("uv")) {
        const auto value=dict_array(attributes,"uv");require_array<double>(value,2,"uv");
        if(value.shape(1)!=2 || static_cast<std::size_t>(value.shape(0))!=result.mesh.corner_vertices.size())
            throw py::value_error("uv must have shape (host_loop_count,2)");
        const auto coordinates=elements<double>(value);host_uv.reserve(coordinates.size()/2);
        for(std::size_t i=0;i<coordinates.size();i+=2) host_uv.push_back({coordinates[i],coordinates[i+1]});
        host_attributes.uv=std::span<const cleanseams::Vec2>(host_uv);
    }
    const auto mapped=cleanseams::read_host_attributes(snapshot,host_attributes,budget,host_painted);
    py::dict out;const auto& mesh=result.mesh;
    out["positions"]=owned_array(mesh.positions);out["face_offsets"]=owned_array(mesh.face_offsets);
    out["corner_vertices"]=owned_array(mesh.corner_vertices);out["edge_vertices"]=owned_array(mesh.edge_vertices);
    out["host_loop_indices"]=owned_array(mesh.host_loop_indices);
    out["host_edge_to_core"]=owned_array(snapshot.tables().host_edge_to_core);
    out["core_edge_to_host"]=owned_array(snapshot.tables().core_edge_to_host);
    out["geometry_fingerprint"]=mesh.geometry_fingerprint;out["geometry_key"]=snapshot.geometry_key();
    out["nonuniform_scale"]=result.nonuniform_scale;
    out["seams"]=owned_array(mapped.seams);out["locked"]=owned_array(mapped.locked);
    out["anchors"]=owned_array(mapped.anchors);out["pins"]=owned_array(mapped.pins);
    if(mapped.uv) out["uv"]=owned_array(*mapped.uv);else out["uv"]=py::none();
    std::vector<cleanseams::TransferIndex> pin_corners;std::vector<cleanseams::Vec2> pin_uv;
    pin_corners.reserve(mapped.pin_constraints.size());pin_uv.reserve(mapped.pin_constraints.size());
    for(const auto& pin:mapped.pin_constraints) {pin_corners.push_back(pin.corner);pin_uv.push_back(pin.uv);}
    out["pin_corners"]=owned_array(pin_corners);out["pin_uv"]=owned_array(pin_uv);
    py::list painted;
    for(const auto& a:mapped.painted) {
        py::dict entry;entry["name"]=a.name;entry["domain"]=a.domain;entry["data_type"]=a.data_type;entry["values"]=owned_array(a.values);
        painted.append(entry);
    }
    out["painted"]=painted;return out;
}
py::tuple copy_buffers(const py::array& points, const py::array& offsets, const py::array& corners) {
    require_array<double>(points,2,"positions");
    require_array<std::uint64_t>(offsets,1,"face_offsets");
    require_array<std::int64_t>(corners,1,"corner_vertices");
    if (points.shape(1) != 3) throw py::value_error("positions must have shape (n,3)");
    const auto mesh=cleanseams::copy_mesh_buffers(
        {static_cast<const double*>(points.data()), static_cast<std::size_t>(points.size())},
        {static_cast<const std::uint64_t*>(offsets.data()), static_cast<std::size_t>(offsets.size())},
        {static_cast<const std::int64_t*>(corners.data()), static_cast<std::size_t>(corners.size())});
    py::array_t<double> p(std::vector<py::ssize_t>{static_cast<py::ssize_t>(mesh.positions.size()/3),3});
    py::array_t<std::uint64_t> o(static_cast<py::ssize_t>(mesh.face_offsets.size()));
    py::array_t<std::int64_t> c(static_cast<py::ssize_t>(mesh.corner_vertices.size()));
    if(!mesh.positions.empty()) std::memcpy(p.mutable_data(),mesh.positions.data(),mesh.positions.size()*sizeof(double));
    if(!mesh.face_offsets.empty()) std::memcpy(o.mutable_data(),mesh.face_offsets.data(),mesh.face_offsets.size()*sizeof(std::uint64_t));
    if(!mesh.corner_vertices.empty()) std::memcpy(c.mutable_data(),mesh.corner_vertices.data(),mesh.corner_vertices.size()*sizeof(std::int64_t));
    return py::make_tuple(std::move(p),std::move(o),std::move(c));
}
// ------------------------------------------------------------------ sessions (T19, T25)
py::array_t<double> uv_array(const std::vector<cleanseams::Vec2>& values) {return owned_array(values);}
// Reads every buffer synchronously (GIL held) into owned native storage: no Python
// object or view reaches the worker, so caller buffers may change after start().
// Session::start then validates every index, length and value against the snapshot.
cleanseams::RunInput run_input(const py::dict& d) {
    static const std::set<std::string> known{"positions","face_offsets","corner_vertices","seam_mask","lock_mask","uv","pin_corners","pin_uv",
        "anchors","options","geometry_revision","geometry_fingerprint","painted_attributes","operation","selected_faces"};
    for(const auto& item:d) {
        if(!py::isinstance<py::str>(item.first)) throw py::value_error("run input keys must be strings");
        if(!known.contains(py::cast<std::string>(item.first))) throw py::value_error("unknown run input: "+py::cast<std::string>(item.first));
    }
    cleanseams::RunInput in;
    if(d.contains("operation"))in.operation=dict_str(d,"operation","run input");
    if(d.contains("selected_faces")) {
        const auto selected=dict_array(d,"selected_faces");require_array<bool>(selected,1,"selected_faces");
        in.selected_faces=copy_elements<std::uint8_t>(selected);
    }
    const auto positions=dict_array(d,"positions");require_array<double>(positions,2,"positions");require_width(positions,3,"positions");
    const auto p=elements<double>(positions);in.mesh.positions.reserve(p.size()/3);
    for(std::size_t i=0;i<p.size();i+=3) {
        if(!std::isfinite(p[i]) || !std::isfinite(p[i+1]) || !std::isfinite(p[i+2])) throw py::value_error("positions must be finite");
        in.mesh.positions.push_back({p[i],p[i+1],p[i+2]});
    }
    const auto offsets=dict_array(d,"face_offsets");require_array<std::uint64_t>(offsets,1,"face_offsets");in.mesh.face_offsets=copy_elements<std::uint64_t>(offsets);
    const auto corners=dict_array(d,"corner_vertices");require_array<std::int64_t>(corners,1,"corner_vertices");in.mesh.corner_vertices=copy_elements<std::int64_t>(corners);
    for(const auto* key:{"seam_mask","lock_mask"}) if(d.contains(key)) {
        const auto m=dict_array(d,key);require_array<bool>(m,1,key);
        (std::string(key)=="seam_mask"?in.seam_mask:in.lock_mask)=copy_elements<std::uint8_t>(m);
    }
    if(d.contains("uv") && !d["uv"].is_none()) in.uv=copy_uv(dict_array(d,"uv"),"uv",false); // non-finite: no baseline (F6)
    if(d.contains("pin_corners") || d.contains("pin_uv")) {
        const auto c=dict_array(d,"pin_corners");require_array<std::int64_t>(c,1,"pin_corners");const auto uv=copy_uv(dict_array(d,"pin_uv"),"pin_uv",true);
        const auto ids=elements<std::int64_t>(c);if(ids.size()!=uv.size()) throw py::value_error("pin_corners and pin_uv lengths differ");
        for(std::size_t i=0;i<ids.size();++i) in.pins.push_back({ids[i],uv[i]});
    }
    if(d.contains("anchors")) {const auto a=dict_array(d,"anchors");require_array<std::int64_t>(a,1,"anchors");in.anchors=copy_elements<std::int64_t>(a);}
    if(d.contains("options")) in.options=cleanseams::pipeline_options_from_json(parse_json(dict_str(d,"options","run input"),"options"));
    if(d.contains("geometry_revision")) in.geometry_revision=dict_uint(d,"geometry_revision","run input");
    if(d.contains("geometry_fingerprint")) in.mesh.geometry_fingerprint=dict_str(d,"geometry_fingerprint","run input");
    // Painted attributes in core order (read_host_mesh's "painted"): copied, never interpreted.
    if(d.contains("painted_attributes")) {
        for(const auto& item:py::cast<py::list>(d["painted_attributes"])) {
            const auto entry=py::cast<py::dict>(item);cleanseams::PaintedAttribute a;
            a.name=dict_str(entry,"name","painted_attributes");a.domain=dict_str(entry,"domain","painted_attributes");
            a.data_type=dict_str(entry,"data_type","painted_attributes");
            const auto values=dict_array(entry,"values");require_array<double>(values,1,"painted_attributes.values");
            a.values=copy_elements<double>(values);in.constraints.attributes.push_back(std::move(a));
        }
    }
    return in;
}
py::dict identity_dict(const cleanseams::HostIdentity& id) {
    py::dict out;out["version"]=id.version;out["session_id"]=id.session_id;out["request_id"]=id.request_id;
    out["geometry_revision"]=id.geometry_revision;out["state_revision"]=id.state_revision;out["chart_keys"]=id.chart_keys;return out;
}
py::object request_dict(const std::optional<cleanseams::HostRequest>& request) {
    if(!request) return py::none();
    py::dict out;out["identity"]=identity_dict(request->identity);out["kind"]=request->kind;
    if(request->region_proposal) {
        const auto& p=*request->region_proposal;
        out["positions"]=owned_array(p.mesh.positions);out["face_offsets"]=owned_array(p.mesh.face_offsets);
        out["corner_vertices"]=owned_array(p.mesh.corner_vertices);out["scratch_to_original_corner"]=owned_array(p.scratch_to_original_corner);
        out["method"]="SMART_PROJECT";
    }
    if(request->solve) {
        const auto& s=*request->solve;
        out["positions"]=owned_array(s.mesh.positions);out["face_offsets"]=owned_array(s.mesh.face_offsets);
        out["corner_vertices"]=owned_array(s.mesh.corner_vertices);out["scratch_to_original_corner"]=owned_array(s.scratch_to_original_corner);
        out["initial_uv"]=s.initial_uv?py::object(uv_array(*s.initial_uv)):py::none();
        std::vector<cleanseams::TransferIndex> pin_corners;std::vector<cleanseams::Vec2> pin_uv;
        for(const auto& pin:s.pins) {pin_corners.push_back(pin.corner);pin_uv.push_back(pin.uv);}
        out["pin_corners"]=owned_array(pin_corners);out["pin_uv"]=uv_array(pin_uv);out["method"]=s.method;out["iterations"]=s.iterations;
    }
    return out;
}
// A host answer: every field typed and named in the error (the session then checks the
// identity, revision, kind, method and map against the open request).
cleanseams::HostResponse response_of(const py::dict& d) {
    static const std::set<std::string> known{"identity","kind","success","method","uv","error"};
    for(const auto& item:d) {
        if(!py::isinstance<py::str>(item.first) || !known.contains(py::cast<std::string>(item.first)))
            throw py::value_error("unknown host response field: "+py::cast<std::string>(py::str(item.first)));
    }
    if(!d.contains("identity") || !py::isinstance<py::dict>(d["identity"])) throw py::value_error("host response identity must be a dict");
    cleanseams::HostResponse r;const auto id=py::cast<py::dict>(d["identity"]);
    static const std::set<std::string> identity_fields{"version","session_id","request_id","geometry_revision","state_revision","chart_keys"};
    for(const auto& item:id)if(!py::isinstance<py::str>(item.first) || !identity_fields.contains(py::cast<std::string>(item.first)))
        throw py::value_error("unknown host response identity field");
    if(id.contains("state_revision"))r.identity.state_revision=dict_uint(id,"state_revision","identity");
    if(id.contains("chart_keys"))r.identity.chart_keys=py::cast<std::vector<std::string>>(id["chart_keys"]);
    if(id.contains("version")) r.identity.version=dict_uint(id,"version","identity");
    r.identity.session_id=dict_uint(id,"session_id","identity");r.identity.request_id=dict_uint(id,"request_id","identity");
    r.identity.geometry_revision=dict_uint(id,"geometry_revision","identity");
    r.kind=d.contains("kind")?dict_str(d,"kind","host response"):"solve";r.success=dict_bool(d,"success","host response");
    if(r.kind=="region_proposal")r.method="SMART_PROJECT";
    if(d.contains("method") && !d["method"].is_none()) r.method=dict_str(d,"method","host response"); // the method the host actually used
    if(d.contains("uv") && !d["uv"].is_none()) r.uv=copy_uv(dict_array(d,"uv"),"uv",false); // ranked invalid natively (F9)
    if(d.contains("error") && !d["error"].is_none()) r.error=cleanseams::ErrorInfo{"host_failure","",dict_str(d,"error","host response")};
    return r;
}
py::dict result_dict(const cleanseams::RunResult& r) {
    py::dict out;out["accepted"]=r.accepted;out["reason"]=r.reason;out["cancelled"]=r.cancelled;out["timed_out"]=r.timed_out;
    out["partial"]=r.partial;out["uv_valid"]=r.uv_valid;out["topology_valid"]=r.topology_valid;out["budget_satisfied"]=r.budget_satisfied;
    py::array_t<bool> seams(static_cast<py::ssize_t>(r.seams.size()));
    for(std::size_t i=0;i<r.seams.size();++i) seams.mutable_data()[i]=r.seams[i]!=0;
    out["seams"]=seams;out["uv"]=r.uv?py::object(uv_array(*r.uv)):py::none();
    std::vector<cleanseams::TransferIndex> exempt;for(const auto& e:r.exempt_faces) exempt.push_back(e.face);
    out["exempt_faces"]=owned_array(exempt);out["invalid_charts"]=owned_array(r.invalid_charts);out["backend"]=r.backend;
    out["geometry_fingerprint"]=r.geometry_fingerprint;
    out["chart_of_face"]=owned_array(r.chart_of_face);out["chart_of_corner"]=owned_array(r.chart_of_corner);
    out["changed_faces"]=owned_array(r.changed_faces);
    out["edge_vertices"]=owned_array(r.edge_vertices);
    out["report_json"]=r.report.dump();
    py::list diagnostics;for(const auto& d:r.diagnostics) diagnostics.append(py::make_tuple(d.code,d.field,d.message));
    out["diagnostics"]=diagnostics;return out;
}
// Every Session the module hands out is registered, so extension shutdown and interpreter
// exit can close them (cancel, wake, join) before the module goes away: no native worker
// outlives its owner (T25). The registry is deliberately leaked so that it is usable
// during interpreter finalization, after static destructors may have run.
struct BoundSession;
struct Registry {std::mutex mutex;std::set<BoundSession*> sessions;};
Registry& registry() {static auto* instance=new Registry;return *instance;}
struct BoundSession {
    cleanseams::Session session;
    std::atomic<bool> open{true};
    explicit BoundSession(cleanseams::SessionConfig config):session(std::move(config)) {
        std::lock_guard lock(registry().mutex);registry().sessions.insert(this);
    }
    ~BoundSession() {
        {std::lock_guard lock(registry().mutex);registry().sessions.erase(this);}
        session.close();
    }
    BoundSession(const BoundSession&)=delete;
    BoundSession& operator=(const BoundSession&)=delete;
    void close() {session.close();open=false;}
};
// Called with the GIL held: no Python thread can drop a session while it is being
// closed, and workers never take the GIL, so joining them here cannot deadlock.
std::size_t close_all_sessions() {
    std::vector<BoundSession*> live;
    {
        std::lock_guard lock(registry().mutex);
        for(auto* s:registry().sessions)if(s->open)live.push_back(s);
    }
    for(auto* s:live)s->close();
    return live.size();
}
std::size_t live_sessions() {
    std::lock_guard lock(registry().mutex);
    return static_cast<std::size_t>(std::count_if(registry().sessions.begin(),registry().sessions.end(),[](auto* s) {return s->open.load();}));
}
}
PYBIND11_MODULE(cleanseams_native, module) {
    register_mesh_analysis(module);
    module.def("build_info",[] {return py::module_::import("json").attr("loads")(cleanseams::build_info().dump());});
    module.def("host_geometry_fingerprint",[](const py::array& p,const py::array& e,const py::array& v,
        const py::array& le,const py::array& s,const py::array& t,std::uint64_t budget) {
        return cleanseams::host_geometry_fingerprint(host_view(p,e,v,le,s,t),budget);
    },py::arg("positions").noconvert(),py::arg("edge_vertices").noconvert(),py::arg("loop_vertices").noconvert(),
      py::arg("loop_edges").noconvert(),py::arg("face_loop_start").noconvert(),py::arg("face_loop_total").noconvert(),
      py::arg("memory_budget")=2147483648ULL);
    module.def("read_host_mesh",&read_host,py::arg("positions").noconvert(),py::arg("edge_vertices").noconvert(),
        py::arg("loop_vertices").noconvert(),py::arg("loop_edges").noconvert(),py::arg("face_loop_start").noconvert(),
        py::arg("face_loop_total").noconvert(),py::arg("object_to_world").noconvert(),py::arg("object_scale").noconvert(),
        py::arg("world_space")=false,py::arg("memory_budget")=2147483648ULL,py::arg("attributes")=py::dict());
    // T29: the host attribute names a run transfers (documented names or configured aliases).
    module.def("painted_attribute_names",[](const std::string& value) {
        py::list names;
        for(const auto& name:cleanseams::painted_attribute_names(cleanseams::pipeline_options_from_json(parse_json(value,"options")))) names.append(name);
        return names;
    },py::arg("options_json"));
    module.def("pipeline_options_json",[](const std::string& value){return cleanseams::pipeline_options_json(cleanseams::pipeline_options_from_json(parse_json(value,"options"))).dump();});
    // `view_direction` (world) and `object_to_world` (row-major 3x3) are the host's view
    // data for hide_from="view"; the native side resolves the object-space direction.
    module.def("resolve_ui_options_json",[](const std::string& value,const py::object& view_direction,const py::object& object_to_world){
        std::optional<cleanseams::ViewContext> view;
        if(!view_direction.is_none() || !object_to_world.is_none()) {
            if(view_direction.is_none() || object_to_world.is_none()) throw py::value_error("view_direction and object_to_world go together");
            const auto direction=numbers(view_direction,"view_direction");const auto matrix=numbers(object_to_world,"object_to_world");
            if(direction.size()!=3 || matrix.size()!=9) throw py::value_error("view_direction needs 3 values and object_to_world 9");
            cleanseams::ViewContext context;std::copy(direction.begin(),direction.end(),context.direction_world.begin());
            std::copy(matrix.begin(),matrix.end(),context.object_to_world.begin());view=context;
        }
        return cleanseams::pipeline_options_json(cleanseams::resolve_ui_options(cleanseams::ui_settings_from_json(parse_json(value,"ui")),view)).dump();
    },py::arg("value"),py::arg("view_direction")=py::none(),py::arg("object_to_world")=py::none());
    // Packing settings only: no view, hide direction or island target.
    module.def("resolve_ui_pack_options_json",[](const std::string& value){
        return cleanseams::pipeline_options_json(cleanseams::resolve_ui_pack_options(cleanseams::ui_settings_from_json(parse_json(value,"ui")))).dump();
    },py::arg("value"));
    module.def("apply_ui_profile_json",[](const std::string& value,const std::string& profile){return cleanseams::ui_settings_json(cleanseams::apply_ui_profile(cleanseams::ui_settings_from_json(parse_json(value,"ui")),profile)).dump();});
    module.def("option_registry_json",[]{return cleanseams::option_registry().dump();});
    module.def("contract_schemas_json",[]{return cleanseams::contract_schemas().dump();});
    module.doc()="CleanSeams native engine: owned buffers, options and persistent sessions";
    // Every call copies at the boundary; the worker never sees a Python object and the
    // GIL is released only around pure native waits.
    py::class_<BoundSession>(module,"Session",py::module_local())
        .def(py::init([](const std::string& backend,std::uint64_t memory_budget,std::uint64_t worker_count,std::uint64_t batch_max_charts,
                         std::uint64_t batch_max_triangles,const std::vector<std::string>& host_capabilities,const std::string& fragment_unwrap_method,
                         const std::string& host_unwrap_method) {
            cleanseams::SessionConfig config;config.backend=backend;config.memory_budget_bytes=memory_budget;
            // T24: native evaluation threads (1 = reference; 0, the default, = automatic,
            // default_native_threads(): GB performance) and deadline batch bounds.
            config.worker_count=worker_count;config.batch_max_charts=batch_max_charts;config.batch_max_triangles=batch_max_triangles;
            config.host_capabilities=host_capabilities;config.fragment_unwrap_method=fragment_unwrap_method;
            config.host_unwrap_method=host_unwrap_method;
            return std::make_unique<BoundSession>(config);
        }),py::arg("backend")="superlu",py::arg("memory_budget")=2147483648ULL,py::arg("worker_count")=0,
            py::arg("batch_max_charts")=16,py::arg("batch_max_triangles")=50000,
            py::arg("host_capabilities")=std::vector<std::string>{},py::arg("fragment_unwrap_method")="ANGLE_BASED",
            py::arg("host_unwrap_method")="MINIMUM_STRETCH")
        .def("start",[](BoundSession& s,const py::dict& input) {
            if(!s.open)throw cleanseams::ContractError(cleanseams::ErrorCode::stale_state,"session","the session is closed");
            auto in=run_input(input);py::gil_scoped_release release;s.session.start(std::move(in));
        },py::arg("input"))
        .def("poll_progress",[](BoundSession& s) {
            py::list out;for(const auto& e:s.session.poll_progress()) {
                py::dict d;d["stage"]=e.stage;d["done"]=e.done;d["total"]=e.total;d["elapsed_s"]=e.elapsed_s;
                d["eta_s"]=e.eta_s?py::object(py::float_(*e.eta_s)):py::none();d["counters_json"]=e.counters.dump();out.append(d);
            }
            return out;
        })
        .def("next_host_request",[](BoundSession& s) {return request_dict(s.session.next_host_request());})
        .def("submit_host_response",[](BoundSession& s,const py::dict& response) {auto r=response_of(response);s.session.submit_host_response(std::move(r));})
        .def("ready",[](const BoundSession& s) {return s.session.ready();})
        .def("result",[](const BoundSession& s) {return result_dict(s.session.result());})
        // The result only for the geometry the run started on: the host passes the native
        // fingerprint of its live data, taken immediately before it writes anything.
        .def("result_for_apply",[](const BoundSession& s,const std::string& fingerprint) {return result_dict(s.session.result_for_apply(fingerprint));},
             py::arg("geometry_fingerprint"))
        .def("request_cancel",[](BoundSession& s) {s.session.request_cancel();})
        .def("close",[](BoundSession& s) {py::gil_scoped_release release;s.close();})
        .def_property_readonly("open",[](const BoundSession& s) {return s.open.load();});
    module.def("close_all_sessions",&close_all_sessions,"Cancel, wake and join every open session; returns how many were open.");
    module.def("live_sessions",&live_sessions,"Sessions of this module that are not closed yet.");
    module.def("recommended_memory_budget",&cleanseams::recommended_memory_budget,
        "Half the physical memory, at least 2 GiB: the budget the Blender host gives a run (Gate A F3).");
    // Interpreter exit (Blender quitting, a script ending) closes whatever is still open.
    py::module_::import("atexit").attr("register")(py::cpp_function([] {close_all_sessions();}));
    py::register_local_exception<cleanseams::ContractError>(module,"ContractError",PyExc_ValueError);
    // The status line and summaries the UI shows (source ops.summarize/refusal_message/
    // repairs_message); the host only adds the object's non-uniform scale flag.
    module.def("run_messages",[](const std::string& report_json,bool accepted,const std::string& reason,const py::array& invalid_charts,
                                 bool nonuniform_scale) {
        require_array<std::int64_t>(invalid_charts,1,"invalid_charts");
        cleanseams::PipelineResult result;result.accepted=accepted;result.reason=reason;
        result.report=parse_json(report_json,"report_json");
        if(!result.report.is_object()) throw py::value_error("report_json must be a JSON object");
        // The source adapter always records the flag (true or false) in the report.
        result.report["nonuniform_scale"]=nonuniform_scale;
        for(auto c:elements<std::int64_t>(invalid_charts)) {
            if(c<0 || c>=std::numeric_limits<cleanseams::MeshIndex>::max()) throw py::value_error("invalid_charts must be chart identities");
            result.invalid_charts.push_back(static_cast<cleanseams::MeshIndex>(c));
        }
        // A run result's own summary line was formatted from the unrounded report (the
        // report JSON writes non-finite measurements as null, which would print "nan" for
        // an infinite one); the host's non-uniform scale clause is appended as the source's
        // summarize() appends it last.
        const auto& stored=result.report.contains("summary") && result.report.at("summary").is_string()?result.report.at("summary"):nlohmann::json();
        const auto summary=stored.is_string()?stored.get<std::string>()+(nonuniform_scale?
            "; object has non-uniform scale: measured in object space, apply the scale for exact texel density":""):cleanseams::summarize(result);
        py::dict out;
        const auto refusal=cleanseams::refusal_message(reason,result.invalid_charts);
        out["summary"]=summary;out["refusal"]=accepted?py::object(py::none()):py::object(py::str(refusal));
        out["repairs"]=cleanseams::repairs_message(result.report);
        out["status"]=accepted?summary:"kept previous state: "+refusal;
        // T41: the panel's lines, from this report alone (the host draws them as given); an
        // enhanced report keeps the same lines in its own `messages`.
        const auto lines=cleanseams::run_diagnostics(cleanseams::finite_json(result.report));
        if(result.report.contains("messages")) result.report["messages"]=lines;
        py::list messages;
        for(const auto& line:lines) {
            py::dict d;d["level"]=line.at("level").get<std::string>();d["topic"]=line.at("topic").get<std::string>();d["text"]=line.at("text").get<std::string>();
            messages.append(d);
        }
        out["messages"]=messages;out["messages_json"]=lines.dump();
        out["report_json"]=cleanseams::finite_json(result.report).dump();
        return out;
    },py::arg("report_json"),py::arg("accepted"),py::arg("reason"),py::arg("invalid_charts").noconvert(),py::arg("nonuniform_scale")=false);
    // Blender packing diagnostics appended to the retained run messages.
    module.def("packing_messages",[](const std::string& report_json,const std::string& engine,const std::string& messages_json) {
        const auto report=report_json.empty()?nlohmann::json::object():parse_json(report_json,"report_json");
        auto lines=messages_json.empty()?nlohmann::json::array():parse_json(messages_json,"messages_json");
        if(!lines.is_array()) throw py::value_error("messages_json must be a JSON array");
        const auto packing=cleanseams::packing_diagnostics(report,engine);
        for(const auto& line:packing) lines.push_back(line);
        py::dict out;out["status"]=packing.at(0).at("text").get<std::string>();out["level"]=packing.at(0).at("level").get<std::string>();
        out["messages_json"]=lines.dump();
        return out;
    },py::arg("report_json"),py::arg("engine"),py::arg("messages_json")="");
    // Digest of host state a run reads and would replace (seams, locks, anchors, the active
    // UV layer and its pins); the host only gathers the arrays. Each part is (name, array
    // or None), hashed with its length so different splits never collide.
    module.def("state_digest",[](const py::list& parts) {
        std::vector<std::uint8_t> bytes;
        auto put=[&](std::uint64_t n) {for(int i=0;i<8;++i) bytes.push_back(static_cast<std::uint8_t>(n>>(8*i)));};
        for(const auto& item:parts) {
            const auto part=py::cast<py::tuple>(item);
            if(part.size()!=2) throw py::value_error("state_digest parts are (name, array or None)");
            const auto name=py::cast<std::string>(part[0]);put(name.size());bytes.insert(bytes.end(),name.begin(),name.end());
            if(part[1].is_none()) {bytes.push_back(0);continue;}
            if(!py::isinstance<py::array>(part[1])) throw py::value_error("state_digest: NumPy array or None required");
            const auto array=py::cast<py::array>(part[1]);
            if(!(array.flags() & py::array::c_style)) throw py::value_error("state_digest arrays must be C-contiguous");
            bytes.push_back(1);put(static_cast<std::uint64_t>(array.nbytes()));
            const auto* data=static_cast<const std::uint8_t*>(array.data());bytes.insert(bytes.end(),data,data+array.nbytes());
        }
        return cleanseams::sha256_hex(bytes);
    },py::arg("parts"));
    // The host's read-back of the UV layer it just wrote (float32 per loop) against the
    // validated map: the decision is native, the host only gathers the arrays (Gate A F28).
    module.def("check_stored_uv",[](const py::array& uv,const py::array& host_loop_indices,const py::array& stored) {
        require_array<double>(uv,2,"uv");require_width(uv,2,"uv");
        require_array<std::int64_t>(host_loop_indices,1,"host_loop_indices");
        require_array<float>(stored,2,"stored");require_width(stored,2,"stored");
        const auto values=elements<double>(uv);std::vector<cleanseams::Vec2> rows(values.size()/2);
        for(std::size_t i=0;i<rows.size();++i) rows[i]={values[2*i],values[2*i+1]};
        cleanseams::check_stored_uv(rows,elements<std::int64_t>(host_loop_indices),elements<float>(stored));
    },py::arg("uv").noconvert(),py::arg("host_loop_indices").noconvert(),py::arg("stored").noconvert());
    module.def("copy_mesh_buffers",&copy_buffers,py::arg("positions").noconvert(),
        py::arg("face_offsets").noconvert(),py::arg("corner_vertices").noconvert());
}
