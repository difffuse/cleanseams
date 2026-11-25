// SPDX-License-Identifier: GPL-3.0-or-later
#include "detail.hpp"
#include "cleanseams/detail/blake2b.hpp"
namespace cleanseams {
namespace {
void count(detail::Blake2b& hash,std::size_t size) {
    const auto n=static_cast<std::uint64_t>(size);detail::hash_little_endian(hash,std::span(&n,1));
}
void ids(detail::Blake2b& hash,std::span<const MeshIndex> values) {
    std::array<std::int64_t,128> buffer{};
    while(!values.empty()) {
        const auto n=std::min(values.size(),buffer.size());
        std::copy_n(values.begin(),n,buffer.begin());detail::hash_little_endian(hash,std::span<const std::int64_t>(buffer).first(n));
        values=values.subspan(n);
    }
}
}
ChartKey::ChartKey(std::shared_ptr<const MeshSnapshot> mesh,std::span<const MeshIndex> faces,std::vector<MeshIndex> boundary):
    mesh_(std::move(mesh)),faces_(faces.begin(),faces.end()),boundary_(std::move(boundary)),digests_(std::make_shared<Digests>()) {}
const std::string& ChartKey::digest() const {
    std::call_once(digests_->once,[this] {
        detail::Blake2b hash(32);constexpr char domain[]="cleanseams.chart.v1";
        hash.update({reinterpret_cast<const std::uint8_t*>(domain),sizeof(domain)});
        const auto& geometry=mesh_->geometry_key();hash.update({reinterpret_cast<const std::uint8_t*>(geometry.data()),geometry.size()});
        count(hash,faces_.size());ids(hash,faces_);count(hash,boundary_.size());ids(hash,boundary_);digests_->digest=hash.hex_digest();
    });
    return digests_->digest;
}
const std::string& ChartKey::legacy_digest() const {
    std::call_once(digests_->legacy_once,[this] {
        detail::Blake2b legacy(20);count(legacy,faces_.size());ids(legacy,faces_);ids(legacy,boundary_);digests_->legacy=legacy.hex_digest();
    });
    return digests_->legacy;
}
bool ChartKey::operator==(const ChartKey& other) const {
    return digest()==other.digest() && faces_==other.faces_ && boundary_==other.boundary_ &&
        (mesh_==other.mesh_ || mesh_->same_geometry(*other.mesh_));
}
void CutSurface::build_chart_keys() {
    keys_.reserve(data_.n_charts);
    for(MeshIndex ch=0;ch<data_.n_charts;++ch) {
        std::vector<MeshIndex> boundary;boundary.reserve(static_cast<std::size_t>(data_.chart_boundary_sides[ch]));
        for(auto loop:data_.chart_loops[ch]) for(auto corner:data_.loops[loop]) boundary.push_back(corner);
        std::sort(boundary.begin(),boundary.end());
        keys_.push_back(ChartKey(mesh_,data_.chart_faces[ch],std::move(boundary)));
    }
}
}
