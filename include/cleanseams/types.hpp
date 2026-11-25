// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
#include <cstdint>
#include <span>
#include <vector>
#include "cleanseams/export.hpp"
namespace cleanseams {
using Scalar = double;
using TransferIndex = std::int64_t;
using Offset = std::uint64_t;
// Owned transfer buffers. Mesh topology and full request DTOs follow in T02/T03.
struct MeshBuffers {
    std::vector<Scalar> positions; // xyz triples
    std::vector<Offset> face_offsets;
    std::vector<TransferIndex> corner_vertices;
};
CLEANSEAMS_CORE_EXPORT MeshBuffers copy_mesh_buffers(std::span<const Scalar> positions,
    std::span<const Offset> face_offsets, std::span<const TransferIndex> corner_vertices);
// Default native memory allowance of a session or run, in bytes.
inline constexpr std::uint64_t default_memory_budget = 2147483648ULL;
// The memory allowance of the run on the calling thread (Gate A F3). A pipeline run sets
// it from PipelineInput::memory_budget (a session from SessionConfig) for its duration,
// and its worker tasks inherit it; outside a run it is default_memory_budget. Every
// table, solver and exposure allowance not given explicitly defaults to it.
CLEANSEAMS_CORE_EXPORT std::uint64_t active_memory_budget();
class CLEANSEAMS_CORE_EXPORT MemoryBudgetScope {
public:
    explicit MemoryBudgetScope(std::uint64_t budget);
    ~MemoryBudgetScope();
    MemoryBudgetScope(const MemoryBudgetScope&) = delete;
    MemoryBudgetScope& operator=(const MemoryBudgetScope&) = delete;
private:
    std::uint64_t previous_;
};
}
