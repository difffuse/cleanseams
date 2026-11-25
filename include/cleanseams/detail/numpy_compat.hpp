// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
#include <array>
#include <algorithm>
#include <charconv>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <span>
#include <string>
#include <vector>
namespace cleanseams::detail {
// numpy's contiguous float64 1-D add.reduce (ndarray.sum): pairwise summation with
// eight accumulators over blocks of at most 128 values. Axis-0 reductions of
// row-major matrices are sequential instead and do not use this helper.
inline double numpy_sum(std::span<const double> a) {
    const auto n=a.size();
    if(n<8) {double sum=0;for(double x:a)sum+=x;return sum;}
    if(n<=128) {
        double r[8];for(std::size_t j=0;j<8;++j)r[j]=a[j];
        std::size_t i=8;
        for(;i<n-(n%8);i+=8)for(std::size_t j=0;j<8;++j)r[j]+=a[i+j];
        double sum=((r[0]+r[1])+(r[2]+r[3]))+((r[4]+r[5])+(r[6]+r[7]));
        for(;i<n;++i)sum+=a[i];
        return sum;
    }
    auto half=n/2;half-=half%8;
    return numpy_sum(a.first(half))+numpy_sum(a.subspan(half));
}
// numpy.einsum("ij,ij->i") over contiguous 3-vectors on the reference build
// (arm64 NEON, two float64 lanes, fused multiply-add): lanes (a0*b0, a1*b1),
// the tail fused into lane 0, then the lanes summed. Verified bit-for-bit against
// the oracle interpreter. Other SIMD widths
// order this differently; parity claims are per measured platform.
inline double numpy_einsum_dot3(const std::array<double,3>& a,const std::array<double,3>& b) {
    return std::fma(a[2],b[2],a[0]*b[0])+a[1]*b[1];
}
// np.interp for nondecreasing xp: j = (elements <= x) - 1, end values outside,
// the node value on an exact hit, NumPy's slope formula and its NaN fallback.
inline double numpy_interp(double x,std::span<const double> xp,std::span<const double> fp) {
    const auto n=xp.size();
    if(std::isnan(x))return x;
    if(x>xp[n-1])return fp[n-1];
    if(x<xp[0])return fp[0];
    const auto j=static_cast<std::size_t>(std::upper_bound(xp.begin(),xp.end(),x)-xp.begin())-1;
    if(j==n-1 || xp[j]==x)return fp[j];
    const double slope=(fp[j+1]-fp[j])/(xp[j+1]-xp[j]);
    double y=slope*(x-xp[j])+fp[j];
    if(std::isnan(y)) {
        y=slope*(x-xp[j+1])+fp[j+1];
        if(std::isnan(y) && fp[j]==fp[j+1])y=fp[j];
    }
    return y;
}
// np.allclose(a, b) with default rtol 1e-5, atol 1e-8 and equal_nan=False.
inline bool numpy_allclose(std::span<const double> a,std::span<const double> b,double rtol=1e-5,double atol=1e-8) {
    if(a.size()!=b.size())return false;
    for(std::size_t i=0;i<a.size();++i) {
        if(std::isinf(a[i]) || std::isinf(b[i])) {if(a[i]!=b[i])return false;continue;}
        if(!(std::abs(a[i]-b[i])<=atol+rtol*std::abs(b[i])))return false;
    }
    return true;
}
// np.median of a nonempty set of values (mean of the two middle values when even).
inline double numpy_median(std::vector<double> values) {
    std::sort(values.begin(),values.end());const auto n=values.size();
    return n%2?values[n/2]:(values[n/2-1]+values[n/2])/2.0;
}
// Python repr/str of a float (shortest round-trip digits; fixed for decimal
// exponents -4 <= e < 16, otherwise d.ddde[+-]XX), as in source reason strings.
// Shortest decimal digits of a finite nonzero |x| that read back to exactly x (the digits of
// Python repr and of std::to_chars' shortest form), without floating std::to_chars, which
// needs macOS 13.3. For each precision the correctly rounded value is tried first, then its
// two decimal neighbours: when x's rounding interval is asymmetric (a power of two) the
// shortest representation can be a neighbour, and at most one of them reads back.
inline void shortest_decimal(double x,std::string& digits,int& exponent) {
    for(int precision=1;precision<=17;++precision) {
        char buffer[64];std::snprintf(buffer,sizeof buffer,"%.*e",precision-1,x);
        std::string mantissa;const char* e=buffer;
        for(;*e!='e';++e)if(*e>='0' && *e<='9')mantissa+=*e;
        const int lead=std::atoi(e+1);const auto scaled=std::stoull(mantissa);
        const int unit=lead-precision+1; // value = integer * 10^unit
        for(const auto candidate:{scaled,scaled-1,scaled+1}) {
            if(candidate==0)continue;
            const auto text=std::to_string(candidate)+"e"+std::to_string(unit);
            if(std::strtod(text.c_str(),nullptr)!=x)continue;
            digits=std::to_string(candidate);
            exponent=unit+static_cast<int>(digits.size())-1;
            while(digits.size()>1 && digits.back()=='0')digits.pop_back();
            return;
        }
    }
    char buffer[64];std::snprintf(buffer,sizeof buffer,"%.16e",x); // 17 digits always read back
    digits.clear();const char* e=buffer;for(;*e!='e';++e)if(*e>='0' && *e<='9')digits+=*e;
    exponent=std::atoi(e+1);while(digits.size()>1 && digits.back()=='0')digits.pop_back();
}
inline std::string python_float_repr(double x) {
    if(std::isnan(x))return "nan";
    if(std::isinf(x))return x>0?"inf":"-inf";
    std::string sign=std::signbit(x)?"-":"",digits;int exponent=0;
    if(x==0) {digits="0";exponent=0;}
    else shortest_decimal(std::fabs(x),digits,exponent);
    if(exponent>=-4 && exponent<16) {
        const int point=exponent+1,count=static_cast<int>(digits.size());
        if(point<=0)return sign+"0."+std::string(static_cast<std::size_t>(-point),'0')+digits;
        if(point>=count)return sign+digits+std::string(static_cast<std::size_t>(point-count),'0')+".0";
        return sign+digits.substr(0,static_cast<std::size_t>(point))+"."+digits.substr(static_cast<std::size_t>(point));
    }
    std::string out=sign+digits.substr(0,1);
    if(digits.size()>1)out+="."+digits.substr(1);
    const auto magnitude=std::abs(exponent);
    return out+(exponent<0?"e-":"e+")+(magnitude<10?"0":"")+std::to_string(magnitude);
}
// CPython 3.11 hash(frozenset(values)) for distinct nonnegative integers below
// 2**61 - 1, whose own hashes are the values (frozenset_hash with _shuffle_bits;
// empty and dummy table entries cancel). Verified against the oracle interpreter.
template<class T> std::int64_t python_frozenset_hash(std::span<const T> distinct) {
    auto shuffle=[](std::uint64_t h) {return ((h^89869747ULL)^(h<<16))*3644798167ULL;};
    std::uint64_t hash=0;
    for(auto v:distinct)hash^=shuffle(static_cast<std::uint64_t>(v));
    hash^=(static_cast<std::uint64_t>(distinct.size())+1)*1927868237ULL;
    hash^=(hash>>11)^(hash>>25);
    hash=hash*69069ULL+907133923ULL;
    if(hash==~std::uint64_t{0})hash=590923713ULL;
    return static_cast<std::int64_t>(hash);
}
template<class T> std::int64_t python_frozenset_hash(const std::vector<T>& distinct) {return python_frozenset_hash(std::span<const T>(distinct));}
// Python round(x, digits) for a float: the correctly rounded decimal (the C
// library's exact %.Nf conversion on the reference platforms), read back.
inline double python_round(double x,int digits) {
    if(!std::isfinite(x))return x;
    char buffer[512];std::snprintf(buffer,sizeof buffer,"%.*f",digits,x);return std::strtod(buffer,nullptr);
}
// round(numpy.float64, decimals) for decimals 0..8: rint(x * 10**d) / 10**d.
inline double numpy_round(double x,int decimals) {
    static constexpr double p10[]={1e0,1e1,1e2,1e3,1e4,1e5,1e6,1e7,1e8};
    const double f=p10[decimals];return std::nearbyint(x*f)/f;
}
}
