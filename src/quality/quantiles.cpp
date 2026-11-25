// SPDX-License-Identifier: GPL-3.0-or-later AND BSD-3-Clause
// Index sorting adapted from NumPy 2.4.4 aquicksort_/aheapsort_.
// Copyright (c) 2005-2025, NumPy Developers. All rights reserved.
// See LICENSES/numpy-2.4.4.txt and docs/QUALITY.md for provenance.
#include "detail.hpp"
#include <bit>
namespace cleanseams {
namespace {
// Preserve the scalar oracle's unstable equal-value order: midpoint weights
// make replacing this with stable_sort (or another quicksort) observable.
void heap_indices(std::span<std::size_t> indices,std::span<const double> values) {
    auto a=[&](std::size_t i)->std::size_t& {return indices[i-1];};
    auto n=indices.size();
    auto sift=[&](std::size_t start,std::size_t value) {
        auto i=start,j=start*2;
        while(j<=n) {
            if(j<n && values[a(j)]<values[a(j+1)]) ++j;
            if(!(values[value]<values[a(j)])) break;
            a(i)=a(j);i=j;j*=2;
        }
        a(i)=value;
    };
    for(auto l=n/2;l>0;--l) sift(l,a(l));
    while(n>1) {const auto last=a(n);a(n)=a(1);--n;sift(1,last);}
}
void sort_indices(std::vector<std::size_t>& order,std::span<const double> values) {
    if(order.size()<2) return;
    using Index=std::ptrdiff_t;
    struct Partition {Index lo,hi;int depth;};
    std::vector<Partition> stack;
    Index lo=0,hi=static_cast<Index>(order.size())-1;
    int depth=static_cast<int>(std::bit_width(order.size())-1)*2;
    auto at=[&](Index i)->std::size_t& {return order[static_cast<std::size_t>(i)];};
    auto value=[&](Index i){return values[at(i)];};
    for(;;) {
        if(depth<0) heap_indices(std::span(order).subspan(static_cast<std::size_t>(lo),static_cast<std::size_t>(hi-lo+1)),values);
        else {
            while(hi-lo>15) {
                const Index mid=lo+(hi-lo)/2;
                if(value(mid)<value(lo)) std::swap(at(mid),at(lo));
                if(value(hi)<value(mid)) std::swap(at(hi),at(mid));
                if(value(mid)<value(lo)) std::swap(at(mid),at(lo));
                const double pivot=value(mid);Index i=lo,j=hi-1;
                std::swap(at(mid),at(j));
                for(;;) {
                    do {++i;} while(value(i)<pivot);
                    do {--j;} while(pivot<value(j));
                    if(i>=j) break;std::swap(at(i),at(j));
                }
                std::swap(at(i),at(hi-1));--depth;
                if(i-lo<hi-i) {stack.push_back({i+1,hi,depth});hi=i-1;}
                else {stack.push_back({lo,i-1,depth});lo=i+1;}
            }
            for(auto i=lo+1;i<=hi;++i) {
                const auto item=at(i);auto j=i;
                while(j>lo && values[item]<value(j-1)) {at(j)=at(j-1);--j;}
                at(j)=item;
            }
        }
        if(stack.empty()) break;
        const auto next=stack.back();stack.pop_back();lo=next.lo;hi=next.hi;depth=next.depth;
    }
}
}
WeightedDistribution::WeightedDistribution(std::span<const double> values,std::span<const double> weights,std::uint64_t budget):memory_budget_(budget) {
    using quality_detail::require;
    require(values.size()==weights.size(),"quality.weights","value/weight count mismatch");
    require(values.size()<=budget/64 && values.size()<=static_cast<std::size_t>(std::numeric_limits<std::ptrdiff_t>::max()),"quality.quantiles","distribution exceeds memory allowance");
    for(auto v:values) require(std::isfinite(v),"quality.values","distribution values must be finite");
    for(auto w:weights) require(std::isfinite(w),"quality.weights","distribution weights must be finite");
    std::vector<std::size_t> order(values.size());std::iota(order.begin(),order.end(),0);sort_indices(order,values);
    values_.reserve(order.size());centers_.reserve(order.size());
    for(auto i:order) {values_.push_back(values[i]);total_+=weights[i];centers_.push_back(total_-weights[i]/2);}
    require(std::isfinite(total_),"quality.weights","total weight overflow");
    if(total_>0) for(auto w:weights) require(w>=0,"quality.weights","positive-total distribution requires nonnegative weights");
}
double WeightedDistribution::percentile(double q) const {
    quality_detail::require(std::isfinite(q) && q>=0 && q<=100,"quality.quantile","percentile must be in [0,100]");
    if(values_.empty()) return quality_nan;
    if(total_<=0) {
        const double index=(q/100)*static_cast<double>(values_.size()-1);
        const auto lo=static_cast<std::size_t>(index),hi=std::min(lo+1,values_.size()-1);
        const auto f=index-static_cast<double>(lo),diff=values_[hi]-values_[lo];
        return f>=.5?values_[hi]-diff*(1-f):values_[lo]+diff*f;
    }
    const double x=q/100*total_;
    if(x<centers_.front()) return values_.front();
    if(x>=centers_.back()) return values_.back();
    const auto it=std::upper_bound(centers_.begin(),centers_.end(),x);
    const auto hi=static_cast<std::size_t>(it-centers_.begin()),lo=hi-1;
    if(centers_[lo]==x) return values_[lo];
    const double slope=(values_[hi]-values_[lo])/(centers_[hi]-centers_[lo]);
    return slope*(x-centers_[lo])+values_[lo];
}
std::vector<double> WeightedDistribution::percentiles(std::span<const double> queries) const {
    quality_detail::require(queries.size()<=(memory_budget_-values_.size()*2*sizeof(double))/sizeof(double),"quality.quantiles","query output exceeds memory allowance");
    std::vector<double> result;result.reserve(queries.size());for(auto q:queries) result.push_back(percentile(q));return result;
}
}
