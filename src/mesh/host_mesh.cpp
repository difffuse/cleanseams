// SPDX-License-Identifier: GPL-3.0-or-later
#include "cleanseams/host_mesh.hpp"
#include "cleanseams/mesh.hpp"
#include "cleanseams/errors.hpp"
#include "cleanseams/detail/blake2b.hpp"
#include <algorithm>
#include <cmath>
namespace cleanseams {
namespace {
void require(bool condition,const std::string& field,const char* message) {
    if(!condition) throw ContractError(ErrorCode::invalid_input,"host_mesh."+field,message);
}
void validate(const HostMeshView& raw,std::uint64_t budget) {
    require(raw.positions.size()%3==0,"positions","expected xyz triples");
    require(raw.edge_vertices.size()%2==0,"edge_vertices","expected endpoint pairs");
    const auto nv=raw.positions.size()/3,ne=raw.edge_vertices.size()/2;
    const auto nf=raw.face_loop_start.size(),nl=raw.loop_vertices.size();
    validate_mesh_counts(nv,nf,nl,ne,budget);
    require(raw.loop_edges.size()==nl,"loop_edges","loop counts differ");
    require(raw.face_loop_total.size()==nf,"face_loop_total","polygon counts differ");
    for(std::size_t i=0;i<raw.positions.size();++i)
        require(std::isfinite(raw.positions[i]),"positions["+std::to_string(i)+"]","coordinate must be finite");
    auto indices=[&](auto values,std::size_t bound,const char* name) {
        for(std::size_t i=0;i<values.size();++i)
            require(values[i]>=0 && static_cast<std::uint64_t>(values[i])<bound,
                std::string(name)+"["+std::to_string(i)+"]","index out of range");
    };
    indices(raw.edge_vertices,nv,"edge_vertices");indices(raw.loop_vertices,nv,"loop_vertices");
    indices(raw.loop_edges,ne,"loop_edges");
    // Check every interval before allocating the overlap/coverage workspace.
    for(std::size_t f=0;f<nf;++f) {
        const auto start=raw.face_loop_start[f],count=raw.face_loop_total[f];
        require(start>=0 && count>=3 && static_cast<std::uint64_t>(start)<=nl &&
            static_cast<std::uint64_t>(count)<=nl-static_cast<std::uint64_t>(start),
            "polygons["+std::to_string(f)+"]","invalid polygon loop range");
    }
    std::vector<std::uint8_t> seen(nl,0);
    for(std::size_t f=0;f<nf;++f) {
        const auto start=static_cast<std::size_t>(raw.face_loop_start[f]);
        const auto end=start+static_cast<std::size_t>(raw.face_loop_total[f]);
        for(auto loop=start;loop<end;++loop) {
            require(!seen[loop],"polygons["+std::to_string(f)+"]","polygon loop ranges overlap");seen[loop]=1;
            const auto next=loop+1==end?start:loop+1;
            const auto a=raw.loop_vertices[loop],b=raw.loop_vertices[next];
            const auto edge=static_cast<std::size_t>(raw.loop_edges[loop]);
            const auto x=raw.edge_vertices[2*edge],y=raw.edge_vertices[2*edge+1];
            require((a==x && b==y)||(a==y && b==x),"loop_edges["+std::to_string(loop)+"]","edge does not match polygon side");
        }
    }
    require(std::all_of(seen.begin(),seen.end(),[](auto used){return used!=0;}),"polygons","polygon ranges must cover every loop");
}
std::string fingerprint(const HostMeshView& raw) {
    detail::Blake2b hash(12);
    const std::array<std::int64_t,4> counts{static_cast<std::int64_t>(raw.positions.size()/3),
        static_cast<std::int64_t>(raw.edge_vertices.size()/2),static_cast<std::int64_t>(raw.face_loop_start.size()),
        static_cast<std::int64_t>(raw.loop_vertices.size())};
    detail::hash_little_endian(hash,std::span<const std::int64_t>(counts));
    detail::hash_little_endian(hash,raw.positions);
    for(auto values:{raw.edge_vertices,raw.loop_vertices,raw.loop_edges,raw.face_loop_start,raw.face_loop_total})
        detail::hash_little_endian(hash,values);
    return hash.hex_digest();
}
}
std::string host_geometry_fingerprint(const HostMeshView& raw,std::uint64_t budget) {
    validate(raw,budget);return fingerprint(raw);
}
HostMeshInput read_host_mesh(const HostMeshView& raw,const ObjectTransform& transform,std::uint64_t budget) {
    validate(raw,budget);
    for(auto x:transform.object_to_world) require(std::isfinite(x),"object_to_world","matrix must be finite");
    Vec3 scales{};
    for(std::size_t i=0;i<3;++i) {
        require(std::isfinite(transform.object_scale[i]),"object_scale","scale must be finite");
        scales[i]=std::abs(transform.object_scale[i]);
    }
    HostMeshInput result;auto& mesh=result.mesh;
    const auto [minimum,maximum]=std::minmax_element(scales.begin(),scales.end());
    result.nonuniform_scale=*maximum>1.001*std::max(*minimum,1e-12);
    mesh.geometry_fingerprint=fingerprint(raw);
    mesh.positions.reserve(raw.positions.size()/3);
    for(std::size_t i=0;i<raw.positions.size();i+=3) {
        Vec3 point{raw.positions[i],raw.positions[i+1],raw.positions[i+2]};
        if(transform.world_space) {
            const auto original=point;const auto& m=transform.object_to_world;
            for(std::size_t d=0;d<3;++d) {
                point[d]=original[0]*m[4*d]+original[1]*m[4*d+1]+original[2]*m[4*d+2]+m[4*d+3];
                require(std::isfinite(point[d]),"positions["+std::to_string(i/3)+"]","world-space coordinate overflow");
            }
        }
        mesh.positions.push_back(point);
    }
    mesh.corner_vertices.reserve(raw.loop_vertices.size());mesh.host_loop_indices.reserve(raw.loop_vertices.size());
    mesh.face_offsets.reserve(raw.face_loop_start.size()+1);
    for(std::size_t f=0;f<raw.face_loop_start.size();++f) {
        const auto start=static_cast<std::size_t>(raw.face_loop_start[f]);
        const auto end=start+static_cast<std::size_t>(raw.face_loop_total[f]);
        for(auto c=start;c<end;++c) {mesh.corner_vertices.push_back(raw.loop_vertices[c]);mesh.host_loop_indices.push_back(static_cast<TransferIndex>(c));}
        mesh.face_offsets.push_back(mesh.corner_vertices.size());
    }
    mesh.edge_vertices.reserve(raw.edge_vertices.size()/2);
    for(std::size_t e=0;e<raw.edge_vertices.size();e+=2) mesh.edge_vertices.push_back({raw.edge_vertices[e],raw.edge_vertices[e+1]});
    return result;
}
MappedAttributes read_host_attributes(const MeshSnapshot& mesh,const HostAttributes& input,std::uint64_t budget,
                                      std::span<const HostPaintedAttribute> painted) {
    const auto& geometry=mesh.input();const auto& maps=mesh.tables();
    const auto nv=geometry.positions.size(),nc=geometry.corner_vertices.size(),ne=geometry.edge_vertices.size();
    validate_mesh_counts(nv,geometry.face_offsets.size()-1,nc,ne,budget);
    auto mask=[&](auto values,std::size_t count,const char* field) {
        require(values.empty() || values.size()==count,field,"mask length does not match host geometry");
        for(auto value:values) require(value<=1,field,"boolean masks must contain only zero or one");
    };
    mask(input.seams,ne,"seams");mask(input.locked,ne,"locked");
    mask(input.anchors,nv,"anchors");mask(input.pins,nc,"pins");
    // Unpinned UVs may be non-finite: the source then runs without a baseline (Gate A F6).
    // A pinned loop's UV is a constraint and must be finite.
    if(input.uv) {
        require(input.uv->size()==nc,"uv","UV count does not match host loops");
        for(std::size_t i=0;i<nc;++i) if(!input.pins.empty() && input.pins[i]) for(auto x:(*input.uv)[i])
            require(std::isfinite(x),"uv["+std::to_string(i)+"]","a pinned UV coordinate must be finite");
    } else require(std::none_of(input.pins.begin(),input.pins.end(),[](auto x){return x!=0;}),"pins","pinned loops require UV coordinates");
    for(auto loop:maps.corner_to_host_loop)
        require(loop>=0 && static_cast<std::uint64_t>(loop)<nc,"host_loop_indices","loop identity outside supplied host arrays");
    MappedAttributes result;
    result.seams.assign(maps.edge_vertices.size(),0);result.locked.assign(maps.edge_vertices.size(),0);
    result.pins.assign(nc,0);if(input.uv) result.uv=std::vector<Vec2>(nc);
    for(std::size_t edge=0;edge<ne;++edge) {
        const auto core=maps.host_edge_to_core[edge];if(core<0) continue;
        const auto index=static_cast<std::size_t>(core);
        if(!input.seams.empty()) result.seams[index]=input.seams[edge];
        if(!input.locked.empty()) result.locked[index]=input.locked[edge];
    }
    for(std::size_t v=0;v<input.anchors.size();++v) if(input.anchors[v]) result.anchors.push_back(static_cast<TransferIndex>(v));
    for(std::size_t corner=0;corner<nc;++corner) {
        const auto loop=static_cast<std::size_t>(maps.corner_to_host_loop[corner]);
        if(input.uv) (*result.uv)[corner]=(*input.uv)[loop];
        if(!input.pins.empty()) result.pins[corner]=input.pins[loop];
        // Preserve each corner separately, even when several share a geometric
        // vertex. Sector compatibility belongs to the candidate cut topology.
        if(result.pins[corner]) result.pin_constraints.push_back({static_cast<TransferIndex>(corner),(*result.uv)[corner]});
    }
    // Painted attributes (T29): mapped to core identities only; interpretation is native.
    const auto nf=geometry.face_offsets.size()-1;
    for(const auto& a:painted) {
        PaintedAttribute out;out.name=a.name;out.domain=a.domain;out.data_type=a.data_type;
        const std::size_t count=a.domain=="POINT"?nv:a.domain=="EDGE"?ne:a.domain=="FACE"?nf:a.domain=="CORNER"?nc:0;
        const auto field="painted."+a.name;
        if(count && !a.values.empty()) {
            require(a.values.size()==count,field.c_str(),"attribute length does not match its host domain");
            if(a.domain=="EDGE") {
                out.values.assign(maps.edge_vertices.size(),0.0);
                for(std::size_t edge=0;edge<ne;++edge) {const auto core=maps.host_edge_to_core[edge];if(core>=0) out.values[static_cast<std::size_t>(core)]=a.values[edge];}
            } else if(a.domain=="CORNER") {
                out.values.resize(nc);
                for(std::size_t corner=0;corner<nc;++corner) out.values[corner]=a.values[static_cast<std::size_t>(maps.corner_to_host_loop[corner])];
            } else out.values.assign(a.values.begin(),a.values.end());
        }
        result.painted.push_back(std::move(out));
    }
    return result;
}
}
