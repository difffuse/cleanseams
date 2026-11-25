// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
#include "cleanseams/cut_surface.hpp"
#include "cleanseams/errors.hpp"
#include <algorithm>
#include <limits>
#include <numeric>
namespace cleanseams::cut_detail {
inline void require(bool condition,const char* field,const char* message) {
    if(!condition) throw ContractError(ErrorCode::invalid_input,field,message);
}
inline MeshIndex root(std::vector<MeshIndex>& parent,MeshIndex value) {
    auto r=value;while(parent[r]!=r) r=parent[r];
    while(parent[value]!=r) {const auto next=parent[value];parent[value]=r;value=next;}
    return r;
}
inline void unite(std::vector<MeshIndex>& parent,MeshIndex a,MeshIndex b) {
    a=root(parent,a);b=root(parent,b);if(a!=b) parent[std::max(a,b)]=std::min(a,b);
}
inline MeshIndex labels(std::vector<MeshIndex>& parent,std::vector<MeshIndex>& output,std::vector<MeshIndex>& scratch) {
    scratch.assign(parent.size(),std::numeric_limits<MeshIndex>::max());output.resize(parent.size());MeshIndex count=0;
    for(std::size_t i=0;i<parent.size();++i) {
        const auto r=root(parent,static_cast<MeshIndex>(i));
        if(scratch[r]==std::numeric_limits<MeshIndex>::max()) scratch[r]=count++;
        output[i]=scratch[r];
    }
    return count;
}
template<class Enumerate> IndexRows rows(std::size_t count,Enumerate enumerate) {
    IndexRows result;result.offsets.assign(count+1,0);
    enumerate([&](MeshIndex row,MeshIndex){++result.offsets[static_cast<std::size_t>(row)+1];});
    std::partial_sum(result.offsets.begin(),result.offsets.end(),result.offsets.begin());
    result.values.resize(static_cast<std::size_t>(result.offsets.back()));auto cursor=result.offsets;
    enumerate([&](MeshIndex row,MeshIndex value){result.values[static_cast<std::size_t>(cursor[row]++)]=value;});
    return result;
}
}
