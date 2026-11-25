// SPDX-License-Identifier: GPL-3.0-or-later
#include "cleanseams/mesh.hpp"
#include "cleanseams/errors.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <initializer_list>
#include <numeric>
#include <set>
#include <string>
#if defined(_WIN32)
#include <windows.h>
#else
#include <unistd.h>
#endif
#if defined(__APPLE__)
#include <sys/sysctl.h>
#endif
namespace cleanseams {
namespace {
void require(bool condition,const std::string& field,const char* message) {
    if(!condition) throw ContractError(ErrorCode::invalid_input,field,message);
}
void validate(const MeshInput& in,std::uint64_t budget) {
    require(!in.face_offsets.empty(),"mesh.face_offsets","offsets must include initial zero");
    validate_mesh_counts(in.positions.size(),in.face_offsets.size()-1,in.corner_vertices.size(),in.edge_vertices.size(),budget);
    require(in.face_offsets.front()==0 && in.face_offsets.back()==in.corner_vertices.size() &&
        std::is_sorted(in.face_offsets.begin(),in.face_offsets.end()),"mesh.face_offsets","offsets must span corners monotonically from zero");
    for(std::size_t v=0;v<in.positions.size();++v) for(auto x:in.positions[v])
        require(std::isfinite(x),"mesh.positions["+std::to_string(v)+"]","coordinate must be finite");
    for(std::size_t c=0;c<in.corner_vertices.size();++c) {
        const auto v=in.corner_vertices[c];
        require(v>=0 && static_cast<std::uint64_t>(v)<in.positions.size(),"mesh.corner_vertices["+std::to_string(c)+"]","vertex index out of range");
    }
    for(std::size_t f=0;f+1<in.face_offsets.size();++f) {
        const auto begin=in.face_offsets[f],end=in.face_offsets[f+1];
        const auto path="mesh.faces["+std::to_string(f)+"]";
        require(end-begin>=3,path,"face must have at least three distinct vertices");
        std::set<TransferIndex> seen;
        for(auto c=begin;c<end;++c) require(seen.insert(in.corner_vertices[static_cast<std::size_t>(c)]).second,path,"repeated face vertex");
    }
    for(std::size_t e=0;e<in.edge_vertices.size();++e) for(auto v:in.edge_vertices[e])
        require(v>=0 && static_cast<std::uint64_t>(v)<in.positions.size(),"mesh.edge_vertices["+std::to_string(e)+"]","host vertex index out of range");
    auto identities=[](const auto& ids,std::size_t count,const char* field) {
        require(ids.empty() || ids.size()==count,field,"host identity count mismatch");
        std::set<TransferIndex> seen;
        for(auto id:ids) require(id>=0 && seen.insert(id).second,field,"host identities must be unique and nonnegative");
    };
    identities(in.host_edge_indices,in.edge_vertices.size(),"mesh.host_edge_indices");
    identities(in.host_loop_indices,in.corner_vertices.size(),"mesh.host_loop_indices");
}
template<class Enumerate> IndexRows make_rows(std::size_t count,Enumerate enumerate) {
    IndexRows rows;rows.offsets.assign(count+1,0);
    enumerate([&](MeshIndex row,MeshIndex){++rows.offsets[static_cast<std::size_t>(row)+1];});
    std::partial_sum(rows.offsets.begin(),rows.offsets.end(),rows.offsets.begin());
    rows.values.resize(static_cast<std::size_t>(rows.offsets.back()));
    auto cursor=rows.offsets;
    enumerate([&](MeshIndex row,MeshIndex value){rows.values[static_cast<std::size_t>(cursor[row]++)]=value;});
    return rows;
}
}
namespace {
// base + sum(count * stride), saturating at UINT64_MAX.
std::uint64_t charge(std::uint64_t base,std::initializer_list<std::pair<std::uint64_t,std::uint64_t>> terms) {
    constexpr auto top=std::numeric_limits<std::uint64_t>::max();
    for(const auto& [count,stride]:terms) {
        if(stride && count>(top-base)/stride)return top;
        base+=count*stride;
    }
    return base;
}
std::uint64_t triangles(const TableCounts& c) {return c.corners>2*c.faces?c.corners-2*c.faces:0;}
std::string mib(std::uint64_t bytes) {return std::to_string((bytes+(1u<<20)-1)>>20)+" MiB";}
}
// Per-element strides: 1.5x the envelope of peaks measured on meshes of 2 to 450 cells a
// side (macOS arm64 malloc sizes); with the 64 KiB constants every measured build, down
// to a few faces, is covered at least 1.2x. The old strides (mesh
// 256/128/512/128, cut +512 per corner, quality +1024 per triangle, seams +1024 per
// corner) were 3x to 36x the measured peaks and refused 360k-quad meshes at 2 GiB.
std::uint64_t mesh_table_bytes(const TableCounts& c) {
    return charge(64u<<10,{{c.vertices,280},{c.faces,64},{c.corners,176},{c.host_edges,72}});
}
std::uint64_t cut_table_bytes(const TableCounts& c) {return charge(64u<<10,{{c.faces,272},{c.corners,304}});}
std::uint64_t quality_table_bytes(const TableCounts& c) {return charge(64u<<10,{{c.faces,256},{triangles(c),336}});}
std::uint64_t seam_table_bytes(const TableCounts& c) {return charge(64u<<10,{{c.vertices,328},{c.faces,1120}});}
std::uint64_t run_table_bytes(const TableCounts& c) {
    return charge(0,{{mesh_table_bytes(c),1},{cut_table_bytes(c),1},{quality_table_bytes(c),1},{seam_table_bytes(c),1}});
}
void check_run_memory(const TableCounts& c,std::uint64_t budget) {
    const auto needed=run_table_bytes(c);
    if(needed>budget)
        throw ContractError(ErrorCode::memory_limit,"memory_budget","the run's mesh, cut, quality and seam tables need "+mib(needed)+
            "; the memory budget is "+mib(budget));
}
std::uint64_t recommended_memory_budget() {
    std::uint64_t physical=0;
#if defined(_WIN32)
    MEMORYSTATUSEX info{};info.dwLength=sizeof(info);
    if(GlobalMemoryStatusEx(&info))physical=info.ullTotalPhys;
#elif defined(__APPLE__)
    std::size_t size=sizeof physical;
    if(::sysctlbyname("hw.memsize",&physical,&size,nullptr,0)!=0)physical=0;
#elif defined(_SC_PHYS_PAGES) && defined(_SC_PAGE_SIZE)
    const long pages=::sysconf(_SC_PHYS_PAGES),page=::sysconf(_SC_PAGE_SIZE);
    if(pages>0 && page>0)physical=static_cast<std::uint64_t>(pages)*static_cast<std::uint64_t>(page);
#endif
    return std::max(default_memory_budget,physical/2);
}
TableCounts table_counts(const MeshSnapshot& m) {
    return {m.input().positions.size(),m.input().face_offsets.size()-1,m.input().corner_vertices.size(),m.input().edge_vertices.size()};
}
std::uint64_t validate_mesh_counts(std::uint64_t nv,std::uint64_t nf,std::uint64_t nc,std::uint64_t nh,std::uint64_t budget) {
    const auto maximum=std::numeric_limits<MeshIndex>::max();
    require(nv<maximum && nf<maximum && nc<maximum && nh<maximum,"mesh.counts","mesh exceeds internal index capacity");
    const auto bytes=mesh_table_bytes({nv,nf,nc,nh});
    if(bytes>budget || bytes>std::numeric_limits<std::size_t>::max())
        throw ContractError(ErrorCode::memory_limit,"mesh.memory_budget","mesh tables need "+mib(bytes)+"; the memory budget is "+mib(budget));
    return bytes;
}
MeshSnapshot::MeshSnapshot(const MeshInput& input,std::uint64_t budget) {
    validate(input,budget);memory_budget_=budget;
    storage_allowance_=validate_mesh_counts(input.positions.size(),input.face_offsets.size()-1,input.corner_vertices.size(),input.edge_vertices.size(),budget);
    input_=input;build_identity();build_tables();build_geometry();build_diagnostics();build_triangulation();
}
std::ranges::iota_view<MeshIndex,MeshIndex> MeshSnapshot::face_corners(MeshIndex face) const {
    const auto end=input_.face_offsets.at(static_cast<std::size_t>(face)+1);
    return std::views::iota(static_cast<MeshIndex>(input_.face_offsets[face]),static_cast<MeshIndex>(end));
}
MeshIndex MeshSnapshot::corner_next(std::size_t c) const {return tables_.corner_next.at(c);}
MeshIndex MeshSnapshot::corner_prev(std::size_t c) const {
    const auto f=tables_.corner_face.at(c);
    const auto start=input_.face_offsets[f],end=input_.face_offsets[static_cast<std::size_t>(f)+1];
    return static_cast<MeshIndex>(c==start?end-1:c-1);
}
TransferIndex MeshSnapshot::edge_index(TransferIndex a,TransferIndex b) const {
    if(a<0 || b<0 || static_cast<std::uint64_t>(a)>=input_.positions.size() || static_cast<std::uint64_t>(b)>=input_.positions.size()) return -1;
    const std::array<MeshIndex,2> key{static_cast<MeshIndex>(std::min(a,b)),static_cast<MeshIndex>(std::max(a,b))};
    const auto& edges=tables_.edge_vertices;const auto found=std::lower_bound(edges.begin(),edges.end(),key);
    return found!=edges.end() && *found==key ? static_cast<TransferIndex>(found-edges.begin()) : -1;
}
TransferIndex MeshSnapshot::face_across(MeshIndex face,MeshIndex edge) const {
    const auto faces=tables_.edge_faces[edge];
    if(faces.size()!=2) return -1;
    if(faces[0]==face) return faces[1];if(faces[1]==face) return faces[0];return -1;
}
TransferIndex Fan::straight(MeshIndex edge) const {
    if(!closed || edges.size()%2!=0 || edges.empty()) return -1;
    const auto found=std::find(edges.begin(),edges.end(),edge);if(found==edges.end()) return -1;
    return edges[(static_cast<std::size_t>(found-edges.begin())+edges.size()/2)%edges.size()];
}
void MeshSnapshot::build_tables() {
    auto& t=tables_;const auto nv=input_.positions.size(),nc=input_.corner_vertices.size(),nf=input_.face_offsets.size()-1;
    t.corner_face.resize(nc);t.corner_pos.resize(nc);t.corner_next.resize(nc);t.corner_edge.resize(nc);
    std::vector<std::array<MeshIndex,2>> corner_pairs(nc);
    for(std::size_t f=0;f<nf;++f) for(auto c=input_.face_offsets[f];c<input_.face_offsets[f+1];++c) {
        const auto i=static_cast<std::size_t>(c);
        const auto next=c+1==input_.face_offsets[f+1]?input_.face_offsets[f]:c+1;
        t.corner_face[i]=static_cast<MeshIndex>(f);t.corner_pos[i]=static_cast<MeshIndex>(c-input_.face_offsets[f]);t.corner_next[i]=static_cast<MeshIndex>(next);
        const auto a=input_.corner_vertices[i],b=input_.corner_vertices[static_cast<std::size_t>(next)];
        corner_pairs[i]={static_cast<MeshIndex>(std::min(a,b)),static_cast<MeshIndex>(std::max(a,b))};
    }
    t.edge_vertices=corner_pairs;std::sort(t.edge_vertices.begin(),t.edge_vertices.end());
    t.edge_vertices.erase(std::unique(t.edge_vertices.begin(),t.edge_vertices.end()),t.edge_vertices.end());
    const auto ne=t.edge_vertices.size();
    for(std::size_t c=0;c<nc;++c) t.corner_edge[c]=static_cast<MeshIndex>(std::lower_bound(t.edge_vertices.begin(),t.edge_vertices.end(),corner_pairs[c])-t.edge_vertices.begin());
    t.edge_corners=make_rows(ne,[&](auto add){for(std::size_t c=0;c<nc;++c) add(t.corner_edge[c],static_cast<MeshIndex>(c));});
    t.edge_faces=t.edge_corners;for(auto& c:t.edge_faces.values) c=t.corner_face[c];
    t.edge_boundary.resize(ne);t.edge_manifold.resize(ne);t.edge_consistent.assign(ne,1);t.vert_boundary.assign(nv,0);
    for(std::size_t e=0;e<ne;++e) {
        const auto cs=t.edge_corners[e];t.edge_boundary[e]=cs.size()==1;t.edge_manifold[e]=cs.size()<=2;
        if(cs.size()==2 && input_.corner_vertices[cs[0]]==input_.corner_vertices[cs[1]]) t.edge_consistent[e]=0;
        if(t.edge_boundary[e]) for(auto v:t.edge_vertices[e]) t.vert_boundary[v]=1;
    }
    t.vert_corners=make_rows(nv,[&](auto add){for(std::size_t c=0;c<nc;++c) add(static_cast<MeshIndex>(input_.corner_vertices[c]),static_cast<MeshIndex>(c));});
    t.vert_faces=t.vert_corners;for(auto& c:t.vert_faces.values) c=t.corner_face[c];
    t.vert_edges=make_rows(nv,[&](auto add){for(std::size_t e=0;e<ne;++e) for(auto v:t.edge_vertices[e]) add(v,static_cast<MeshIndex>(e));});
    t.fans.reserve(nv);t.vert_regular.assign(nv,0);
    for(std::size_t v=0;v<nv;++v) {
        auto fan=build_fan(static_cast<MeshIndex>(v));
        if(fan && fan->closed && fan->edges.size()==4)
            t.vert_regular[v]=std::all_of(fan->faces.begin(),fan->faces.end(),[&](auto f){return input_.face_offsets[static_cast<std::size_t>(f)+1]-input_.face_offsets[f]==4;});
        t.fans.push_back(std::move(fan));
    }
    t.host_edge_to_core.reserve(input_.edge_vertices.size());t.core_edge_to_host.assign(ne,-1);
    for(std::size_t h=0;h<input_.edge_vertices.size();++h) {
        const auto pair=input_.edge_vertices[h];const auto e=edge_index(pair[0],pair[1]);t.host_edge_to_core.push_back(e);
        if(e>=0) t.core_edge_to_host[static_cast<std::size_t>(e)]=input_.host_edge_indices.empty()?static_cast<TransferIndex>(h):input_.host_edge_indices[h];
    }
    t.corner_to_host_loop=input_.host_loop_indices;
    if(t.corner_to_host_loop.empty()) {t.corner_to_host_loop.resize(nc);std::iota(t.corner_to_host_loop.begin(),t.corner_to_host_loop.end(),TransferIndex{0});}
}
std::optional<Fan> MeshSnapshot::build_fan(MeshIndex v) const {
    const auto& t=tables_;const auto corners=t.vert_corners[v];if(corners.empty()) return std::nullopt;
    std::map<MeshIndex,MeshIndex> next,prev,by_next;
    for(auto c:corners) {
        const auto f=t.corner_face[c];if(!next.emplace(f,t.corner_edge[c]).second) return std::nullopt;
        prev.emplace(f,t.corner_edge[corner_prev(c)]);
    }
    for(auto e:t.vert_edges[v]) if(!t.edge_manifold[e] || !t.edge_consistent[e]) return std::nullopt;
    std::optional<MeshIndex> boundary_start;
    for(const auto& [f,e]:next) {
        if(!by_next.emplace(e,f).second) return std::nullopt;
        if(t.edge_boundary[e]) {if(boundary_start) return std::nullopt;boundary_start=f;}
    }
    auto f=boundary_start.value_or(next.begin()->first);Fan fan;fan.closed=!boundary_start;
    fan.edges.push_back(next.at(f));std::set<MeshIndex> seen;
    while(true) {
        fan.faces.push_back(f);seen.insert(f);const auto e=prev.at(f);
        if(fan.closed && e==fan.edges.front()) break;
        fan.edges.push_back(e);
        if(t.edge_boundary[e]) {if(fan.closed) return std::nullopt;break;}
        const auto found=by_next.find(e);if(found==by_next.end() || seen.contains(found->second)) return std::nullopt;f=found->second;
    }
    if(fan.faces.size()!=corners.size()) return std::nullopt;
    return fan;
}
}
