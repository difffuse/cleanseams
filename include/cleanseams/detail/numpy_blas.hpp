// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
#include "cleanseams/types.hpp"
#include <cmath>
#include <cstddef>
#include <span>
#include <vector>
// Float64 products that the reference NumPy 2.4.4 routes to Apple Accelerate on arm64
// (np.show_config(): blas/lapack "accelerate"). Each model was measured bit for bit on
// the oracle interpreter; other BLAS builds
// order these operations differently, so parity claims are per measured platform.
// Every function compiles under -ffp-contract=off: fused steps are explicit std::fma.
namespace cleanseams::detail {
// np.dot(x, y) of two 1-D float64 vectors when either operand is strided (a column of
// an (n, 2) array, as in the shoelace sums of seam_metrics.island_rects): four
// accumulators by index modulo 4, later elements fused, reduced as (a0+a2)+(a1+a3).
inline double accelerate_strided_dot(std::span<const double> x,std::span<const double> y) {
    double acc[4]={0,0,0,0};bool used[4]={false,false,false,false};
    for(std::size_t i=0;i<x.size();++i) {
        const auto j=i%4;
        acc[j]=used[j]?std::fma(x[i],y[i],acc[j]):x[i]*y[i];used[j]=true;
    }
    return (acc[0]+acc[2])+(acc[1]+acc[3]);
}
// A @ B for a C-contiguous (n, 2) A and a C-contiguous (2, 2) B (orient_islands'
// rotations): rows are fused, out = fma(a1, b1j, a0*b0j), except the last row of an
// odd row count below 1600, which is p0*b0j + p1*b1j.
inline std::vector<Vec2> accelerate_rows_times_2x2(std::span<const Vec2> a,const std::array<Vec2,2>& b) {
    std::vector<Vec2> out(a.size());const auto n=a.size();
    for(std::size_t i=0;i<n;++i)
        for(std::size_t j=0;j<2;++j)
            out[i][j]=(n%2 && i==n-1 && n<1600)?a[i][0]*b[0][j]+a[i][1]*b[1][j]:std::fma(a[i][1],b[1][j],a[i][0]*b[0][j]);
    return out;
}
// A @ C.T for C-contiguous (n, 2) A and (n, 2) C (min_area_rect's projections onto
// edge directions and normals, a square product): for n < 16 columns in complete blocks
// of four are plain, the remaining columns fused except the last row of an odd row count;
// from n = 16 every entry is fused. Returns row-major (n, n).
inline std::vector<double> accelerate_square_projection(std::span<const Vec2> a,std::span<const Vec2> c) {
    const auto n=a.size(),full=n-n%4;std::vector<double> out(n*c.size());
    for(std::size_t i=0;i<n;++i)
        for(std::size_t j=0;j<c.size();++j) {
            const bool plain=n<16 && (j<full || (n%2 && i==n-1));
            out[i*c.size()+j]=plain?a[i][0]*c[j][0]+a[i][1]*c[j][1]:std::fma(a[i][1],c[j][1],a[i][0]*c[j][0]);
        }
    return out;
}
}
