// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
// GB performance: the permutation of COO entries sorted by (major, minor, entry index).
//
// The sparse assembly sorted entry indices with a comparison sort whose explicit index
// tie-break makes the order total, so the permutation is unique. Two stable counting
// passes, least significant key first, produce exactly that permutation in linear time:
// the first pass orders by minor key keeping index order, the second by major key keeping
// the first pass's order. Duplicate summation therefore still runs in original COO order
// and every sum is bit-for-bit the one the comparison sort gave.
#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <numeric>
#include <span>
#include <tuple>
#include <vector>
namespace cleanseams::solver_detail {
// Keys must lie in [0, extent). An extent far larger than the entry count would make the
// counting arrays dominate; that case keeps the comparison sort (the same permutation).
template<class Key>
std::vector<std::size_t> lexicographic_order(std::span<const Key> major,std::uint64_t major_extent,std::span<const Key> minor,std::uint64_t minor_extent) {
    const auto n=major.size();
    std::vector<std::size_t> order(n);
    const std::uint64_t dense_limit=4*static_cast<std::uint64_t>(n)+1024;
    if(major_extent>dense_limit || minor_extent>dense_limit) {
        std::iota(order.begin(),order.end(),std::size_t{0});
        std::sort(order.begin(),order.end(),[&](std::size_t a,std::size_t b) {
            return std::tie(major[a],minor[a],a)<std::tie(major[b],minor[b],b);
        });
        return order;
    }
    std::vector<std::size_t> start(static_cast<std::size_t>(std::max(major_extent,minor_extent))+1),by_minor(n);
    auto counting_pass=[&](std::span<const Key> key,std::uint64_t extent,auto&& source,std::vector<std::size_t>& out) {
        std::fill(start.begin(),start.begin()+static_cast<std::ptrdiff_t>(extent)+1,std::size_t{0});
        for(std::size_t i=0;i<n;++i)++start[static_cast<std::size_t>(key[i])+1];
        std::partial_sum(start.begin(),start.begin()+static_cast<std::ptrdiff_t>(extent)+1,start.begin());
        for(std::size_t i=0;i<n;++i) {const auto k=source(i);out[start[static_cast<std::size_t>(key[k])]++]=k;}
    };
    counting_pass(minor,minor_extent,[](std::size_t i) {return i;},by_minor);
    counting_pass(major,major_extent,[&](std::size_t i) {return by_minor[i];},order);
    return order;
}
}
