// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
#include "cleanseams/options.hpp"
namespace cleanseams {
using Revision=std::uint64_t;
using RequestId=std::uint64_t;
using EdgeEndpoints=std::array<TransferIndex,2>;
struct MeshInput {
    std::vector<Vec3> positions;
    std::vector<Offset> face_offsets{0};
    std::vector<TransferIndex> corner_vertices;
    std::vector<EdgeEndpoints> edge_vertices;
    std::vector<TransferIndex> host_edge_indices;
    std::vector<TransferIndex> host_loop_indices;
    std::string geometry_fingerprint;
};
struct CornerPin {
    TransferIndex corner=0;
    Vec2 uv{};
};
// A painted host attribute as the host stores it (T29): name, Blender domain and data type,
// and for scalar types one value per element of the domain in core order (empty otherwise).
// The host only transfers it; the native resolver interprets, validates and reports it.
struct PaintedAttribute {
    std::string name;
    std::string domain;
    std::string data_type;
    std::vector<double> values;
};
struct ErrorInfo {
    std::string code="host_failure";
    std::string field;
    std::string message;
};
struct HostIdentity {
    std::uint64_t version=1;
    std::uint64_t session_id=0;
    RequestId request_id=0;
    Revision geometry_revision=0;
    Revision state_revision=0;
    std::vector<std::string> chart_keys;
};
struct HostSolveRequest {
    MeshInput mesh;
    std::vector<TransferIndex> scratch_to_original_corner;
    std::optional<std::vector<Vec2>> initial_uv;
    std::vector<CornerPin> pins;
    std::string method="MINIMUM_STRETCH";
    std::uint64_t iterations=10;
    std::vector<std::string> required_capabilities;
};
// A projection proposes cuts on original region vertices. It is not a fixed-cut solve.
struct HostRegionProposalRequest {
    MeshInput mesh;
    std::vector<TransferIndex> scratch_to_original_corner;
};
struct HostPackRequest {
    MeshInput mesh;
    std::vector<Vec2> uv;
    std::vector<CornerPin> pins;
    PackOptions options;
    std::string backend="blender";
    nlohmann::json payload=nlohmann::json::object();
};
struct HostRequest {
    HostIdentity identity;
    std::string kind="solve";
    std::optional<HostSolveRequest> solve;
    std::optional<HostPackRequest> pack;
    std::optional<HostRegionProposalRequest> region_proposal;
};
struct HostResponse {
    HostIdentity identity;
    std::string kind="solve";
    bool success=false;
    std::optional<std::vector<Vec2>> uv;
    std::string method="MINIMUM_STRETCH";
    std::uint64_t iterations=0;
    std::vector<std::string> capabilities;
    std::vector<std::string> warnings;
    double elapsed_s=0.;
    std::optional<ErrorInfo> error;
    nlohmann::json packing_result=nlohmann::json::object();
};
CLEANSEAMS_CORE_EXPORT HostRequest host_request_from_json(const nlohmann::json&);
CLEANSEAMS_CORE_EXPORT nlohmann::json host_request_json(const HostRequest&);
CLEANSEAMS_CORE_EXPORT HostResponse host_response_from_json(const nlohmann::json&);
CLEANSEAMS_CORE_EXPORT nlohmann::json host_response_json(const HostResponse&);
}
