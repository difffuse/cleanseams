// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
namespace cleanseams {
template<class F> void visit_fields(PaintedAttribute& v,F&& f) {
    f("name",v.name);
    f("domain",v.domain);
    f("data_type",v.data_type);
    f("values",v.values);
}
template<class F> void visit_fields(const PaintedAttribute& v,F&& f) {
    f("name",v.name);
    f("domain",v.domain);
    f("data_type",v.data_type);
    f("values",v.values);
}
template<class F> void visit_fields(ConstraintInput& v,F&& f) {
    f("required_edges",v.required_edges);
    f("forbidden_edges",v.forbidden_edges);
    f("allow_points",v.allow_points);
    f("avoid_points",v.avoid_points);
    f("attributes",v.attributes);
}
template<class F> void visit_fields(const ConstraintInput& v,F&& f) {
    f("required_edges",v.required_edges);
    f("forbidden_edges",v.forbidden_edges);
    f("allow_points",v.allow_points);
    f("avoid_points",v.avoid_points);
    f("attributes",v.attributes);
}
template<class F> void visit_fields(ProposalIdentity& v,F&& f) {
    f("id",v.id);
    f("parent_revision",v.parent_revision);
    f("generator",v.generator);
    f("generator_version",v.generator_version);
    f("affected_chart_keys",v.affected_chart_keys);
    f("added_edges",v.added_edges);
    f("removed_edges",v.removed_edges);
}
template<class F> void visit_fields(const ProposalIdentity& v,F&& f) {
    f("id",v.id);
    f("parent_revision",v.parent_revision);
    f("generator",v.generator);
    f("generator_version",v.generator_version);
    f("affected_chart_keys",v.affected_chart_keys);
    f("added_edges",v.added_edges);
    f("removed_edges",v.removed_edges);
}
template<class F> void visit_fields(MeshInput& v,F&& f) {
    f("positions",v.positions);
    f("face_offsets",v.face_offsets);
    f("corner_vertices",v.corner_vertices);
    f("edge_vertices",v.edge_vertices);
    f("host_edge_indices",v.host_edge_indices);
    f("host_loop_indices",v.host_loop_indices);
    f("geometry_fingerprint",v.geometry_fingerprint);
}
template<class F> void visit_fields(const MeshInput& v,F&& f) {
    f("positions",v.positions);
    f("face_offsets",v.face_offsets);
    f("corner_vertices",v.corner_vertices);
    f("edge_vertices",v.edge_vertices);
    f("host_edge_indices",v.host_edge_indices);
    f("host_loop_indices",v.host_loop_indices);
    f("geometry_fingerprint",v.geometry_fingerprint);
}
template<class F> void visit_fields(CornerPin& v,F&& f) {
    f("corner",v.corner);
    f("uv",v.uv);
}
template<class F> void visit_fields(const CornerPin& v,F&& f) {
    f("corner",v.corner);
    f("uv",v.uv);
}
template<class F> void visit_fields(ErrorInfo& v,F&& f) {
    f("code",v.code);
    f("field",v.field);
    f("message",v.message);
}
template<class F> void visit_fields(const ErrorInfo& v,F&& f) {
    f("code",v.code);
    f("field",v.field);
    f("message",v.message);
}
template<class F> void visit_fields(HostIdentity& v,F&& f) {
    f("version",v.version);
    f("session_id",v.session_id);
    f("request_id",v.request_id);
    f("geometry_revision",v.geometry_revision);
    f("state_revision",v.state_revision);
    f("chart_keys",v.chart_keys);
}
template<class F> void visit_fields(const HostIdentity& v,F&& f) {
    f("version",v.version);
    f("session_id",v.session_id);
    f("request_id",v.request_id);
    f("geometry_revision",v.geometry_revision);
    f("state_revision",v.state_revision);
    f("chart_keys",v.chart_keys);
}
template<class F> void visit_fields(HostSolveRequest& v,F&& f) {
    f("mesh",v.mesh);
    f("scratch_to_original_corner",v.scratch_to_original_corner);
    f("initial_uv",v.initial_uv);
    f("pins",v.pins);
    f("method",v.method);
    f("iterations",v.iterations);
    f("required_capabilities",v.required_capabilities);
}
template<class F> void visit_fields(const HostSolveRequest& v,F&& f) {
    f("mesh",v.mesh);
    f("scratch_to_original_corner",v.scratch_to_original_corner);
    f("initial_uv",v.initial_uv);
    f("pins",v.pins);
    f("method",v.method);
    f("iterations",v.iterations);
    f("required_capabilities",v.required_capabilities);
}
template<class F> void visit_fields(HostRegionProposalRequest& v,F&& f) {
    f("mesh",v.mesh);
    f("scratch_to_original_corner",v.scratch_to_original_corner);
}
template<class F> void visit_fields(const HostRegionProposalRequest& v,F&& f) {
    f("mesh",v.mesh);
    f("scratch_to_original_corner",v.scratch_to_original_corner);
}
template<class F> void visit_fields(HostPackRequest& v,F&& f) {
    f("mesh",v.mesh);
    f("uv",v.uv);
    f("pins",v.pins);
    f("options",v.options);
    f("backend",v.backend);
    f("payload",v.payload);
}
template<class F> void visit_fields(const HostPackRequest& v,F&& f) {
    f("mesh",v.mesh);
    f("uv",v.uv);
    f("pins",v.pins);
    f("options",v.options);
    f("backend",v.backend);
    f("payload",v.payload);
}
template<class F> void visit_fields(HostRequest& v,F&& f) {
    f("identity",v.identity);
    f("kind",v.kind);
    f("solve",v.solve);
    f("pack",v.pack);
    f("region_proposal",v.region_proposal);
}
template<class F> void visit_fields(const HostRequest& v,F&& f) {
    f("identity",v.identity);
    f("kind",v.kind);
    f("solve",v.solve);
    f("pack",v.pack);
    f("region_proposal",v.region_proposal);
}
template<class F> void visit_fields(HostResponse& v,F&& f) {
    f("identity",v.identity);
    f("kind",v.kind);
    f("success",v.success);
    f("uv",v.uv);
    f("method",v.method);
    f("iterations",v.iterations);
    f("capabilities",v.capabilities);
    f("warnings",v.warnings);
    f("elapsed_s",v.elapsed_s);
    f("error",v.error);
    f("packing_result",v.packing_result);
}
template<class F> void visit_fields(const HostResponse& v,F&& f) {
    f("identity",v.identity);
    f("kind",v.kind);
    f("success",v.success);
    f("uv",v.uv);
    f("method",v.method);
    f("iterations",v.iterations);
    f("capabilities",v.capabilities);
    f("warnings",v.warnings);
    f("elapsed_s",v.elapsed_s);
    f("error",v.error);
    f("packing_result",v.packing_result);
}
template<class F> void visit_fields(ExemptFace& v,F&& f) {
    f("face",v.face);
    f("reason",v.reason);
}
template<class F> void visit_fields(const ExemptFace& v,F&& f) {
    f("face",v.face);
    f("reason",v.reason);
}
template<class F> void visit_fields(UnavailableValue& v,F&& f) {
    f("path",v.path);
    f("reason",v.reason);
}
template<class F> void visit_fields(const UnavailableValue& v,F&& f) {
    f("path",v.path);
    f("reason",v.reason);
}
template<class F> void visit_fields(RunResult& v,F&& f) {
    f("version",v.version);
    f("accepted",v.accepted);
    f("reason",v.reason);
    f("topology_valid",v.topology_valid);
    f("uv_valid",v.uv_valid);
    f("budget_satisfied",v.budget_satisfied);
    f("partial",v.partial);
    f("cancelled",v.cancelled);
    f("timed_out",v.timed_out);
    f("geometry_fingerprint",v.geometry_fingerprint);
    f("state_revision",v.state_revision);
    f("seams",v.seams);
    f("uv",v.uv);
    f("pins",v.pins);
    f("exempt_faces",v.exempt_faces);
    f("invalid_charts",v.invalid_charts);
    f("changed_charts",v.changed_charts);
    f("resolved_options",v.resolved_options);
    f("backend",v.backend);
    f("diagnostics",v.diagnostics);
    f("report",v.report);
    f("decision_trace",v.decision_trace);
    f("unavailable",v.unavailable);
    f("chart_of_face",v.chart_of_face);
    f("chart_of_corner",v.chart_of_corner);
    f("changed_faces",v.changed_faces);
    f("edge_vertices",v.edge_vertices);
}
template<class F> void visit_fields(const RunResult& v,F&& f) {
    f("version",v.version);
    f("accepted",v.accepted);
    f("reason",v.reason);
    f("topology_valid",v.topology_valid);
    f("uv_valid",v.uv_valid);
    f("budget_satisfied",v.budget_satisfied);
    f("partial",v.partial);
    f("cancelled",v.cancelled);
    f("timed_out",v.timed_out);
    f("geometry_fingerprint",v.geometry_fingerprint);
    f("state_revision",v.state_revision);
    f("seams",v.seams);
    f("uv",v.uv);
    f("pins",v.pins);
    f("exempt_faces",v.exempt_faces);
    f("invalid_charts",v.invalid_charts);
    f("changed_charts",v.changed_charts);
    f("resolved_options",v.resolved_options);
    f("backend",v.backend);
    f("diagnostics",v.diagnostics);
    f("report",v.report);
    f("decision_trace",v.decision_trace);
    f("unavailable",v.unavailable);
    f("chart_of_face",v.chart_of_face);
    f("chart_of_corner",v.chart_of_corner);
    f("changed_faces",v.changed_faces);
    f("edge_vertices",v.edge_vertices);
}
template<class F> void visit_fields(ProgressEvent& v,F&& f) {
    f("stage",v.stage);
    f("done",v.done);
    f("total",v.total);
    f("elapsed_s",v.elapsed_s);
    f("eta_s",v.eta_s);
    f("counters",v.counters);
}
template<class F> void visit_fields(const ProgressEvent& v,F&& f) {
    f("stage",v.stage);
    f("done",v.done);
    f("total",v.total);
    f("elapsed_s",v.elapsed_s);
    f("eta_s",v.eta_s);
    f("counters",v.counters);
}
template<class F> void visit_fields(SessionConfig& v,F&& f) {
    f("worker_count",v.worker_count);
    f("memory_budget_bytes",v.memory_budget_bytes);
    f("cache_budget_bytes",v.cache_budget_bytes);
    f("trace_level",v.trace_level);
    f("backend",v.backend);
    f("host_capabilities",v.host_capabilities);
    f("fragment_unwrap_method",v.fragment_unwrap_method);
    f("host_unwrap_method",v.host_unwrap_method);
    f("batch_max_charts",v.batch_max_charts);
    f("batch_max_triangles",v.batch_max_triangles);
}
template<class F> void visit_fields(const SessionConfig& v,F&& f) {
    f("worker_count",v.worker_count);
    f("memory_budget_bytes",v.memory_budget_bytes);
    f("cache_budget_bytes",v.cache_budget_bytes);
    f("trace_level",v.trace_level);
    f("backend",v.backend);
    f("host_capabilities",v.host_capabilities);
    f("fragment_unwrap_method",v.fragment_unwrap_method);
    f("host_unwrap_method",v.host_unwrap_method);
    f("batch_max_charts",v.batch_max_charts);
    f("batch_max_triangles",v.batch_max_triangles);
}
template<class F> void visit_fields(RunInput& v,F&& f) {
    f("version",v.version);
    f("coordinate_space",v.coordinate_space);
    f("mesh",v.mesh);
    f("geometry_revision",v.geometry_revision);
    f("seam_edges",v.seam_edges);
    f("locked_edges",v.locked_edges);
    f("seam_mask",v.seam_mask);
    f("lock_mask",v.lock_mask);
    f("uv",v.uv);
    f("pins",v.pins);
    f("anchors",v.anchors);
    f("constraints",v.constraints);
    f("object_to_world",v.object_to_world);
    f("view_direction_world",v.view_direction_world);
    f("options",v.options);
    f("operation",v.operation);
    f("selected_faces",v.selected_faces);
}
template<class F> void visit_fields(const RunInput& v,F&& f) {
    f("version",v.version);
    f("coordinate_space",v.coordinate_space);
    f("mesh",v.mesh);
    f("geometry_revision",v.geometry_revision);
    f("seam_edges",v.seam_edges);
    f("locked_edges",v.locked_edges);
    f("seam_mask",v.seam_mask);
    f("lock_mask",v.lock_mask);
    f("uv",v.uv);
    f("pins",v.pins);
    f("anchors",v.anchors);
    f("constraints",v.constraints);
    f("object_to_world",v.object_to_world);
    f("view_direction_world",v.view_direction_world);
    f("options",v.options);
    f("operation",v.operation);
    f("selected_faces",v.selected_faces);
}
}
