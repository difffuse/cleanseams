// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
#include <cstdint>
#include <optional>
#include <vector>
namespace cleanseams::detail {
// Iteration order of a CPython 3.11 `set` of non-negative ints (hash(i) == i),
// reproduced where the source sums or iterates such a set and the order changes
// float rounding or decisions. Supports add() and set(other_set) construction;
// validated against the oracle interpreter.
class PyIntSet {
public:
    PyIntSet():table_(8) {}
    void add(std::uint64_t key) {
        std::size_t mask=table_.size()-1,i=static_cast<std::size_t>(key)&mask;std::uint64_t perturb=key;
        while(true) {
            std::size_t probes=i+linear_probes<=mask?linear_probes:0,j=i;
            while(true) {
                auto& entry=table_[j];
                if(!entry) {
                    entry=key;++fill_;++used_;
                    if(fill_*5<mask*3)return;
                    resize(used_>50000?used_*2:used_*4);return;
                }
                if(*entry==key)return;
                ++j;
                if(probes==0)break;
                --probes;
            }
            perturb>>=perturb_shift;i=(i*5+1+static_cast<std::size_t>(perturb))&mask;
        }
    }
    // set(other): CPython's set_merge into an empty set.
    static PyIntSet copy_of(const PyIntSet& other) {
        PyIntSet out;
        if(!other.used_)return out;
        if((out.fill_+other.used_)*5>=(out.table_.size()-1)*3)out.resize((out.used_+other.used_)*2);
        if(out.table_.size()==other.table_.size()) {out.table_=other.table_;out.fill_=other.fill_;out.used_=other.used_;return out;}
        out.fill_=out.used_=other.used_;
        for(const auto& k:other.table_)if(k)insert_clean(out.table_,*k);
        return out;
    }
    std::vector<std::uint64_t> order() const {
        std::vector<std::uint64_t> out;out.reserve(used_);for(const auto& k:table_)if(k)out.push_back(*k);return out;
    }
    std::size_t size() const {return used_;}
private:
    static constexpr std::size_t linear_probes=9,perturb_shift=5;
    std::vector<std::optional<std::uint64_t>> table_;
    std::size_t fill_=0,used_=0;
    static void insert_clean(std::vector<std::optional<std::uint64_t>>& table,std::uint64_t key) {
        const std::size_t mask=table.size()-1;std::size_t i=static_cast<std::size_t>(key)&mask;std::uint64_t perturb=key;
        while(true) {
            if(!table[i]) {table[i]=key;return;}
            if(i+linear_probes<=mask)for(std::size_t j=1;j<=linear_probes;++j)if(!table[i+j]) {table[i+j]=key;return;}
            perturb>>=perturb_shift;i=(i*5+1+static_cast<std::size_t>(perturb))&mask;
        }
    }
    void resize(std::size_t minused) {
        std::size_t size=8;while(size<=minused)size<<=1;
        std::vector<std::optional<std::uint64_t>> old(size);old.swap(table_);
        for(const auto& k:old)if(k)insert_clean(table_,*k);
    }
};
}
