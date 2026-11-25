// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
#include "cleanseams/export.hpp"
#include <array>
#include <bit>
#include <cstdint>
#include <span>
#include <string>
#include <type_traits>
namespace cleanseams::detail {
// Unkeyed sequential BLAKE2b for geometry identity. Algorithm: RFC 7693,
// sections 2-3, https://www.rfc-editor.org/rfc/rfc7693.html . No platform crypto
// dependency and no heap buffer proportional to the message size.
class CLEANSEAMS_CORE_EXPORT Blake2b {
public:
    explicit Blake2b(std::size_t digest_bytes);
    void update(std::span<const std::uint8_t> bytes);
    std::string hex_digest() const;
private:
    std::array<std::uint64_t,8> state_{};
    std::array<std::uint8_t,128> block_{};
    std::uint64_t low_=0,high_=0;
    std::size_t used_=0,digest_bytes_=0;
    void compress(bool final);
};
template<class T> void hash_little_endian(Blake2b& hash,std::span<const T> values) {
    static_assert(std::is_arithmetic_v<T> && (sizeof(T)==4 || sizeof(T)==8));
    if constexpr(std::endian::native==std::endian::little) {
        hash.update({reinterpret_cast<const std::uint8_t*>(values.data()),values.size_bytes()});
    } else {
        for(const auto& value:values) {
            auto bytes=std::bit_cast<std::array<std::uint8_t,sizeof(T)>>(value);
            for(std::size_t i=0;i<bytes.size()/2;++i) std::swap(bytes[i],bytes[bytes.size()-i-1]);
            hash.update(bytes);
        }
    }
}
}
