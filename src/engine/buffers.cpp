// SPDX-License-Identifier: GPL-3.0-or-later
#include "cleanseams/types.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
namespace cleanseams {
MeshBuffers copy_mesh_buffers(std::span<const Scalar> positions,
        std::span<const Offset> face_offsets, std::span<const TransferIndex> corner_vertices) {
    static_assert(sizeof(Scalar) == 8 && std::numeric_limits<Scalar>::is_iec559);
    if (positions.size() % 3 != 0)
        throw std::invalid_argument("positions must contain xyz triples");
    if (face_offsets.empty() || face_offsets.front() != 0 ||
            face_offsets.back() != corner_vertices.size() ||
            !std::is_sorted(face_offsets.begin(), face_offsets.end()))
        throw std::invalid_argument("face offsets must span corners monotonically from zero");
    if (!std::all_of(positions.begin(), positions.end(), [](Scalar v){return std::isfinite(v);}))
        throw std::invalid_argument("positions must be finite");
    const auto vertices = positions.size() / 3;
    for (const auto index : corner_vertices)
        if (index < 0 || static_cast<std::uint64_t>(index) >= vertices)
            throw std::invalid_argument("corner vertex index out of range");
    // Copy after validation; caller memory cannot become session-owned by accident.
    return {{positions.begin(), positions.end()}, {face_offsets.begin(), face_offsets.end()},
            {corner_vertices.begin(), corner_vertices.end()}};
}
}
