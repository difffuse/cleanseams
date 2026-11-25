// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
#include "cleanseams/mesh.hpp"
#include "cleanseams/solver.hpp"
#include <nlohmann/json.hpp>
#include <filesystem>
#include <optional>
#include <string>
#include <string_view>
namespace cleanseams {
// ------------------------------------------------------------------ OBJ (core/io.py)
// Positions and polygons only; per-corner UVs when every face carries vt indices.
struct ObjMesh {
    std::vector<Vec3> positions;
    std::vector<std::vector<MeshIndex>> faces;
    std::optional<std::vector<Vec2>> uv; // per corner in face order
};
// Source read_obj: one-based indices, negative indices relative to the elements defined
// so far, '#' comments, "v" x y z, "vt" u v, "f" v[/vt[/vn]]...; other records ignored.
// Text mode: "\n", "\r\n" and a lone "\r" end a line; str.split() whitespace; Python
// int()/float() number syntax, including PEP 515 underscores.
CLEANSEAMS_CORE_EXPORT ObjMesh parse_obj(std::string_view text,std::uint64_t memory_budget=default_memory_budget);
CLEANSEAMS_CORE_EXPORT ObjMesh read_obj(const std::filesystem::path&,std::uint64_t memory_budget=default_memory_budget);
// Source write_obj: "%.17g" coordinates, "f v/t" with running UV indices, and one "l a b"
// line per seam edge (original vertex endpoints) when given.
CLEANSEAMS_CORE_EXPORT std::string format_obj(std::span<const Vec3> positions,const std::vector<std::vector<MeshIndex>>& faces,
    std::optional<std::span<const Vec2>> uv=std::nullopt,std::span<const EdgeEndpoints> seam_lines={});
CLEANSEAMS_CORE_EXPORT void write_obj(const std::filesystem::path&,std::span<const Vec3> positions,const std::vector<std::vector<MeshIndex>>& faces,
    std::optional<std::span<const Vec2>> uv=std::nullopt,std::span<const EdgeEndpoints> seam_lines={});
// ------------------------------------------------------------------ case files
// SHA-256 (FIPS 180-4) of a byte buffer, lowercase hex; binary-array descriptors use it.
CLEANSEAMS_CORE_EXPORT std::string sha256_hex(std::span<const std::uint8_t> bytes);
// A "cleanseams.case" file (plan section 25.2): inline arrays or checked binary-array
// descriptors. Unknown top-level fields (id, family, parameters, options, packing, ...)
// are kept in `extra` untouched.
struct CaseData {
    std::uint64_t version=1;
    std::string coordinate_space="object";
    MeshInput mesh;
    std::vector<EdgeEndpoints> seam_edges,locked_edges; // original vertex endpoint pairs
    std::vector<TransferIndex> anchors;
    std::vector<CornerPin> pins;
    std::optional<std::vector<Vec2>> uv;
    nlohmann::json extra=nlohmann::json::object();
    std::string operation="pipeline";
    std::vector<std::uint8_t> selected_faces;
};
struct CaseLimits {
    std::uint64_t memory_budget=default_memory_budget; // total bytes of JSON text plus binary arrays
};
// Reads a case; descriptors resolve against `base` (the case's directory). Rejects
// absolute paths, "..", symlink escapes, unsupported dtypes/versions, size or count
// mismatches, truncated files, hash mismatches, nonfinite values and budget overruns.
CLEANSEAMS_CORE_EXPORT CaseData case_from_json(const nlohmann::json&,const std::filesystem::path& base,const CaseLimits& limits={});
CLEANSEAMS_CORE_EXPORT CaseData read_case(const std::filesystem::path&,const CaseLimits& limits={});
// Arrays with at least `inline_limit` scalar entries are written as descriptors under
// <directory of path>/arrays/<stem>-<field>.<ext>; smaller ones inline. Existing files are
// only replaced when their content differs from what is written.
CLEANSEAMS_CORE_EXPORT nlohmann::json case_json(const CaseData&,const std::filesystem::path& path,std::size_t inline_limit=1u<<20);
CLEANSEAMS_CORE_EXPORT void write_case(const std::filesystem::path&,const CaseData&,std::size_t inline_limit=1u<<20);
// Faces of a case mesh as index lists (OBJ writing, packing requests).
CLEANSEAMS_CORE_EXPORT std::vector<std::vector<MeshIndex>> case_faces(const MeshInput&);
// Reads a whole regular file under a byte limit (JSON inputs and responses).
CLEANSEAMS_CORE_EXPORT std::string read_text_file(const std::filesystem::path&,std::uint64_t limit,const char* field);
// Untrusted JSON text (case files, CLI inputs, external packer output, host strings).
// nlohmann::json copies, dumps and compares values recursively, so the nesting depth is
// bounded before a value exists (real documents nest at most 8 levels); syntax errors,
// invalid UTF-8 and number overflow become ContractErrors on `field`.
inline constexpr std::size_t max_json_depth=128;
CLEANSEAMS_CORE_EXPORT nlohmann::json parse_untrusted_json(std::string_view text,const std::string& field,
    std::size_t max_depth=max_json_depth);
}
