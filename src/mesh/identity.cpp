// SPDX-License-Identifier: GPL-3.0-or-later
#include "cleanseams/mesh.hpp"
#include "cleanseams/detail/blake2b.hpp"
#include <bit>
namespace cleanseams {
namespace {
void length(detail::Blake2b& hash,std::size_t size) {
    const auto n=static_cast<std::uint64_t>(size);detail::hash_little_endian(hash,std::span(&n,1));
}
template<class T> void array(detail::Blake2b& hash,const std::vector<T>& values) {
    length(hash,values.size());detail::hash_little_endian(hash,std::span<const T>(values));
}
}
void MeshSnapshot::build_identity() {
    detail::Blake2b hash(32);
    constexpr char domain[]="cleanseams.mesh.geometry.v1";
    hash.update({reinterpret_cast<const std::uint8_t*>(domain),sizeof(domain)});
    length(hash,input_.positions.size());
    for(const auto& point:input_.positions) detail::hash_little_endian(hash,std::span<const double>(point));
    array(hash,input_.face_offsets);array(hash,input_.corner_vertices);
    length(hash,input_.edge_vertices.size());
    for(const auto& edge:input_.edge_vertices) detail::hash_little_endian(hash,std::span<const TransferIndex>(edge));
    array(hash,input_.host_edge_indices);array(hash,input_.host_loop_indices);
    geometry_key_=hash.hex_digest();
}
bool MeshSnapshot::same_geometry(const MeshSnapshot& other) const {
    // Exact equality is mandatory after a cache hash hit, including signed zero.
    // External restore labels are deliberately not trusted as content identity.
    const auto& b=other.input_;
    if(input_.positions.size()!=b.positions.size() || input_.face_offsets!=b.face_offsets ||
        input_.corner_vertices!=b.corner_vertices || input_.edge_vertices!=b.edge_vertices ||
        input_.host_edge_indices!=b.host_edge_indices || input_.host_loop_indices!=b.host_loop_indices) return false;
    for(std::size_t v=0;v<input_.positions.size();++v) for(std::size_t d=0;d<3;++d)
        if(std::bit_cast<std::uint64_t>(input_.positions[v][d])!=std::bit_cast<std::uint64_t>(b.positions[v][d])) return false;
    return true;
}
}
