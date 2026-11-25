// SPDX-License-Identifier: GPL-3.0-or-later
#include "cleanseams/detail/blake2b.hpp"
#include <algorithm>
#include <stdexcept>
namespace cleanseams::detail {
namespace {
constexpr std::array<std::uint64_t,8> iv{
    0x6a09e667f3bcc908ULL,0xbb67ae8584caa73bULL,0x3c6ef372fe94f82bULL,0xa54ff53a5f1d36f1ULL,
    0x510e527fade682d1ULL,0x9b05688c2b3e6c1fULL,0x1f83d9abfb41bd6bULL,0x5be0cd19137e2179ULL};
constexpr std::array<std::array<std::size_t,16>,10> schedule{{
    {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15},
    {14,10,4,8,9,15,13,6,1,12,0,2,11,7,5,3},
    {11,8,12,0,5,2,15,13,10,14,3,6,7,1,9,4},
    {7,9,3,1,13,12,11,14,2,6,5,10,4,0,15,8},
    {9,0,5,7,2,4,10,15,14,1,11,12,6,8,3,13},
    {2,12,6,10,0,11,8,3,4,13,7,5,15,14,1,9},
    {12,5,1,15,14,13,4,10,0,7,6,3,9,2,8,11},
    {13,11,7,14,12,1,3,9,5,0,15,4,8,6,2,10},
    {6,15,14,9,11,3,0,8,12,2,13,7,1,4,10,5},
    {10,2,8,4,7,6,1,5,15,11,9,14,3,12,13,0}}};
constexpr std::array<std::array<std::size_t,4>,8> lanes{{
    {0,4,8,12},{1,5,9,13},{2,6,10,14},{3,7,11,15},
    {0,5,10,15},{1,6,11,12},{2,7,8,13},{3,4,9,14}}};
}
Blake2b::Blake2b(std::size_t bytes):state_(iv),digest_bytes_(bytes) {
    if(bytes==0 || bytes>64) throw std::invalid_argument("BLAKE2b digest size must be 1..64 bytes");
    state_[0]^=0x01010000ULL|static_cast<std::uint64_t>(bytes);
}
void Blake2b::compress(bool final) {
    const auto old=low_;low_+=used_;if(low_<old) ++high_;
    std::array<std::uint64_t,16> words{},work{};
    for(std::size_t i=0;i<8;++i) {work[i]=state_[i];work[i+8]=iv[i];}
    work[12]^=low_;work[13]^=high_;if(final) work[14]=~work[14];
    for(std::size_t i=0;i<block_.size();++i) words[i/8]|=std::uint64_t(block_[i])<<(8*(i%8));
    for(std::size_t round=0;round<12;++round) for(std::size_t lane=0;lane<lanes.size();++lane) {
        const auto indices=lanes[lane];
        auto& a=work[indices[0]];auto& b=work[indices[1]];
        auto& c=work[indices[2]];auto& d=work[indices[3]];
        const auto& order=schedule[round%10];
        for(std::size_t half=0;half<2;++half) {
            a+=b+words[order[2*lane+half]];
            d=std::rotr(d^a,half==0?32:16);
            c+=d;b=std::rotr(b^c,half==0?24:63);
        }
    }
    for(std::size_t i=0;i<state_.size();++i) state_[i]^=work[i]^work[i+8];
}
void Blake2b::update(std::span<const std::uint8_t> bytes) {
    while(!bytes.empty()) {
        // Keep even an exactly full last block pending until we know whether
        // another byte follows. Empty input still receives a final block.
        if(used_==block_.size()) {compress(false);used_=0;}
        const auto count=std::min(block_.size()-used_,bytes.size());
        std::copy_n(bytes.begin(),count,block_.begin()+static_cast<std::ptrdiff_t>(used_));
        used_+=count;bytes=bytes.subspan(count);
    }
}
std::string Blake2b::hex_digest() const {
    auto copy=*this;
    std::fill(copy.block_.begin()+static_cast<std::ptrdiff_t>(copy.used_),copy.block_.end(),0);
    copy.compress(true);
    constexpr char digits[]="0123456789abcdef";
    std::string result(digest_bytes_*2,'0');
    for(std::size_t i=0;i<digest_bytes_;++i) {
        const auto byte=static_cast<std::uint8_t>(copy.state_[i/8]>>(8*(i%8)));
        result[2*i]=digits[byte>>4];result[2*i+1]=digits[byte&15];
    }
    return result;
}
}
