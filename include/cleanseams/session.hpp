// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
#include "cleanseams/result.hpp"
#include "cleanseams/proposals.hpp"
#include <memory>
#include <string_view>
namespace cleanseams {
struct SessionConfig {
    // Native evaluation threads for the run (1 is the reference serial schedule, at most
    // 64; 0, the default, is automatic: default_native_threads(), resolved by
    // resolve_session_config). Every count gives the same results. Host requests are
    // always serviced by the host on its own main thread.
    std::uint64_t worker_count=0;
    std::uint64_t memory_budget_bytes=default_memory_budget;
    std::optional<std::uint64_t> cache_budget_bytes;
    std::string trace_level="stages";
    std::string backend="superlu";
    // Additive opt-in; existing host solve clients do not implicitly support projection.
    std::vector<std::string> host_capabilities;
    std::string fragment_unwrap_method="ANGLE_BASED";
    // Bounds of deadline work batches (evaluations and host solves), T24.
    std::uint64_t batch_max_charts=16;
    std::uint64_t batch_max_triangles=50000;
    // Primary generation solve; native retries still request their own method.
    std::string host_unwrap_method="MINIMUM_STRETCH";
};
struct RunInput {
    std::uint64_t version=1;
    std::string coordinate_space="object";
    MeshInput mesh;
    Revision geometry_revision=0;
    std::vector<EdgeEndpoints> seam_edges;
    std::vector<EdgeEndpoints> locked_edges;
    std::vector<std::uint8_t> seam_mask;
    std::vector<std::uint8_t> lock_mask;
    std::optional<std::vector<Vec2>> uv;
    std::vector<CornerPin> pins;
    std::vector<TransferIndex> anchors;
    ConstraintInput constraints;
    std::array<double,16> object_to_world{1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1};
    std::optional<Vec3> view_direction_world;
    PipelineOptions options;
    std::string operation="pipeline";
    std::vector<std::uint8_t> selected_faces;
};
// The host stores UVs in float32. After writing the core map `uv` (one row per core corner)
// to host loops `host_loop_indices[c]`, it reads the loops back as `stored` (two float32 per
// host loop); every mapped loop must be exactly the float32 rounding of the map, else a
// host_failure ContractError (field apply.uv_storage) and the host rolls the write back.
CLEANSEAMS_CORE_EXPORT void check_stored_uv(std::span<const Vec2> uv,std::span<const TransferIndex> host_loop_indices,std::span<const float> stored);
CLEANSEAMS_CORE_EXPORT RunInput run_input_from_json(const nlohmann::json&);
// Painted constraints (T29) against a mesh of `edges` edges and `vertices` vertices: typed
// masks of zeros and ones in core edge order, finite weights in core vertex order, finite
// raw attribute values. Sessions check them in start() and the CLI checks a case file's, for
// enhanced_v1 runs only: the other policies neither read nor check them (GB M2).
// Domains, types, clamping and conflicts are the native resolver's (enhanced_v1 only).
CLEANSEAMS_CORE_EXPORT void check_constraint_input(const ConstraintInput&,std::uint64_t edges,std::uint64_t vertices);
// T41: the `constraints` record of a case file (the run input's wire form; absent fields are
// empty), schema-checked and then checked against the mesh as above.
CLEANSEAMS_CORE_EXPORT ConstraintInput constraint_input_from_json(const nlohmann::json&,std::uint64_t edges,std::uint64_t vertices);
CLEANSEAMS_CORE_EXPORT nlohmann::json run_input_json(const RunInput&);
CLEANSEAMS_CORE_EXPORT SessionConfig resolve_session_config(SessionConfig);
CLEANSEAMS_CORE_EXPORT nlohmann::json contract_schemas();
// Interface only in T02. Stateful execution and definitions are owned by T19.
class CLEANSEAMS_CORE_EXPORT Session {
public:
    explicit Session(SessionConfig config);
    ~Session();
    Session(const Session&)=delete;
    Session& operator=(const Session&)=delete;
    void start(RunInput input);
    std::vector<ProgressEvent> poll_progress();
    std::optional<HostRequest> next_host_request();
    void submit_host_response(HostResponse response);
    bool ready() const;
    RunResult result() const;
    // The result, only when `current_geometry_fingerprint` (the host's native fingerprint
    // of its live data, taken right before writing) equals the one the run started with;
    // otherwise a stale_state ContractError and nothing may be applied (T25). An accepted
    // map must hold every input pin exactly as given (compared as stored, never after
    // aligning maps); otherwise a constraint_conflict ContractError (field apply.pins).
    RunResult result_for_apply(std::string_view current_geometry_fingerprint) const;
    void request_cancel();
    void close();
private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};
}
