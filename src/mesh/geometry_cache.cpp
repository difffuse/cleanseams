// SPDX-License-Identifier: GPL-3.0-or-later
#include "cleanseams/detail/mesh_math.hpp"
#include "cleanseams/detail/numpy_compat.hpp"
#include "cleanseams/errors.hpp"
#include <algorithm>
#include <limits>
#include <numeric>
namespace cleanseams {
using namespace detail;
namespace {
void finite(double x,const char* field) {
    if(!std::isfinite(x)) throw ContractError(ErrorCode::invalid_input,field,"geometry arithmetic exceeds finite float64 range");
}
void finite(const Vec3& v,const char* field) {for(auto x:v) finite(x,field);}
}
std::vector<Vec3> MeshSnapshot::face_points(MeshIndex f) const {
    std::vector<Vec3> points;
    for(auto c=input_.face_offsets.at(f);c<input_.face_offsets.at(static_cast<std::size_t>(f)+1);++c)
        points.push_back(input_.positions[static_cast<std::size_t>(input_.corner_vertices[static_cast<std::size_t>(c)])]);
    return points;
}
void MeshSnapshot::build_geometry() {
    auto& g=geometry_;const auto nf=input_.face_offsets.size()-1,nv=input_.positions.size(),ne=tables_.edge_vertices.size();
    g.edge_length.reserve(ne);
    for(const auto& e:tables_.edge_vertices) {const auto length=norm(sub(input_.positions[e[1]],input_.positions[e[0]]));finite(length,"mesh.edge_length");g.edge_length.push_back(length);}
    g.face_area.resize(nf);g.face_normal.resize(nf);g.face_center.resize(nf);
    std::vector<Vec3> weighted(nf);g.vert_normal.resize(nv);g.vert_normal_ok.resize(nv);
    for(std::size_t f=0;f<nf;++f) {
        const auto points=face_points(static_cast<MeshIndex>(f));Vec3 n{},center{};
        for(std::size_t j=0;j<points.size();++j) {n=add(n,cross(points[j],points[(j+1)%points.size()]));center=add(center,points[j]);}
        weighted[f]=mul(n,.5);g.face_area[f]=norm(weighted[f]);finite(g.face_area[f],"mesh.face_area");
        g.face_normal[f]=g.face_area[f]>0?divide(weighted[f],std::max(g.face_area[f],1e-300)):Vec3{};
        g.face_center[f]=divide(center,static_cast<double>(points.size()));finite(g.face_center[f],"mesh.face_center");
    }
    for(std::size_t c=0;c<input_.corner_vertices.size();++c) {
        const auto v=static_cast<std::size_t>(input_.corner_vertices[c]);g.vert_normal[v]=add(g.vert_normal[v],weighted[tables_.corner_face[c]]);
    }
    for(std::size_t v=0;v<nv;++v) {
        const auto length=norm(g.vert_normal[v]);finite(length,"mesh.vert_normal");g.vert_normal_ok[v]=length>1e-12;
        g.vert_normal[v]=length>0?divide(g.vert_normal[v],std::max(length,1e-300)):Vec3{};
    }
    g.edge_dihedral.assign(ne,0.);
    for(std::size_t e=0;e<ne;++e) {
        const auto faces=tables_.edge_faces[e];if(faces.size()!=2 || !tables_.edge_consistent[e]) continue;
        const auto a=faces[0],b=faces[1];if(g.face_area[a]<=0 || g.face_area[b]<=0) continue;
        const auto cosine=std::clamp(dot(g.face_normal[a],g.face_normal[b]),-1.,1.);
        const auto angle=std::acos(cosine),side=dot(g.face_normal[a],sub(g.face_center[b],g.face_center[a]));
        g.edge_dihedral[e]=side<0?angle:-angle;
    }
}
MeshIndex MeshSnapshot::edge_other(MeshIndex e,MeshIndex v) const {
    const auto endpoints=tables_.edge_vertices.at(e);
    return endpoints[0]==v?endpoints[1]:endpoints[0];
}
Vec3 MeshSnapshot::edge_dir(MeshIndex v,MeshIndex e) const {
    const auto w=edge_other(e,v);const auto d=sub(input_.positions[w],input_.positions[v]);const auto length=norm(d);
    return length>0?divide(d,length):d;
}
Vec3 MeshSnapshot::tangent_dir(MeshIndex v,MeshIndex e) const {
    const auto d=edge_dir(v,e);if(!geometry_.vert_normal_ok[v]) return d;
    const auto n=geometry_.vert_normal[v],t=sub(d,mul(n,dot(d,n)));const auto length=norm(t);
    return length>1e-9?divide(t,length):d;
}
double MeshSnapshot::mean_edge_length(std::optional<std::span<const MeshIndex>> edges) const {
    // NumPy mean: pairwise sum in the given edge order, then division.
    if(!edges) {const auto n=geometry_.edge_length.size();return n?detail::numpy_sum(geometry_.edge_length)/static_cast<double>(n):0.;}
    std::vector<double> lengths;lengths.reserve(edges->size());for(auto e:*edges) lengths.push_back(geometry_.edge_length.at(e));
    return lengths.empty()?0.:detail::numpy_sum(lengths)/static_cast<double>(lengths.size());
}
void MeshSnapshot::build_diagnostics() {
    auto& g=geometry_;const auto nf=g.face_area.size();
    const auto total=detail::numpy_sum(g.face_area);finite(total,"mesh.total_area"); // source A.mean(): pairwise sum
    const auto mean=nf?total/static_cast<double>(nf):0.;auto sorted=g.face_area;std::sort(sorted.begin(),sorted.end());
    const auto median=nf?(nf%2?sorted[nf/2]:(sorted[nf/2-1]+sorted[nf/2])*.5):0.;
    g.face_aspect.resize(nf);g.face_flatness.resize(nf);g.face_degenerate_kind.resize(nf);g.face_nonplanar.assign(nf,0);
    for(std::size_t f=0;f<nf;++f) {
        const auto area=g.face_area[f];double hi=0.,lo=std::numeric_limits<double>::infinity();
        const auto start=input_.face_offsets[f],end=input_.face_offsets[f+1];
        for(auto c=start;c<end;++c) {const auto length=g.edge_length[tables_.corner_edge[static_cast<std::size_t>(c)]];hi=std::max(hi,length);lo=std::min(lo,length);}
        g.face_aspect[f]=area>0?(hi*hi)/std::max(area,1e-300):std::numeric_limits<double>::infinity();
        const auto flat=lo>0?(area/std::max(hi,1e-300))/std::max(lo,1e-300):std::numeric_limits<double>::infinity();g.face_flatness[f]=flat;
        const bool zero=area<=1e-9*std::max(mean,1e-300),collapsed=!zero && flat<.05,sliver=!zero && !collapsed && area<.02*median && flat<.25;
        bool twisted=false;
        if(!zero && end-start>=4) {
            const auto points=face_points(static_cast<MeshIndex>(f));const bool quad=points.size()==4;
            const auto q=project_polygon(points,g.face_normal[f],quad);
            if(quad) twisted=((orient(q[0],q[1],q[2])<0)!=(orient(q[0],q[2],q[3])<0)) &&
                ((orient(q[1],q[2],q[3])<0)!=(orient(q[1],q[3],q[0])<0));
            else {
                double dev=0.;for(const auto& p:points) dev=std::max(dev,std::abs(dot(sub(p,g.face_center[f]),g.face_normal[f])));
                g.face_nonplanar[f]=dev>.25*std::sqrt(area);twisted=polygon_self_intersects(q);
            }
        }
        g.face_degenerate_kind[f]=zero?1:twisted?2:collapsed?3:sliver?4:0;
    }
    std::vector<MeshIndex> parent(nf);std::iota(parent.begin(),parent.end(),MeshIndex{0});
    auto root=[&](MeshIndex v){while(parent[v]!=v) {parent[v]=parent[parent[v]];v=parent[v];}return v;};
    for(std::size_t e=0;e<tables_.edge_faces.size();++e) {
        const auto faces=tables_.edge_faces[e];for(std::size_t i=1;i<faces.size();++i) {const auto a=root(faces[i-1]),b=root(faces[i]);parent[std::max(a,b)]=std::min(a,b);}
    }
    g.face_component.resize(nf);std::vector<MeshIndex> labels(nf,std::numeric_limits<MeshIndex>::max());
    for(std::size_t f=0;f<nf;++f) {const auto r=root(static_cast<MeshIndex>(f));if(labels[r]==std::numeric_limits<MeshIndex>::max()) labels[r]=g.n_components++;g.face_component[f]=labels[r];}
}
IndexRows MeshSnapshot::face_rings(TransferIndex rings) const {
    if(rings<0) throw ContractError(ErrorCode::invalid_input,"rings","ring count must be nonnegative");
    const auto nf=geometry_.face_area.size();IndexRows result;
    // Per-call workspace and output stay inside the remaining allowance.
    const auto available=memory_budget_-storage_allowance_;
    const auto fixed_bytes=static_cast<std::uint64_t>(nf)*64+4096;
    if(fixed_bytes>available) throw ContractError(ErrorCode::memory_limit,"rings","ring workspace exceeds memory allowance");
    // Leave space for old and new allocations during vector growth.
    const auto maximum_values=(available-fixed_bytes)/(2*sizeof(MeshIndex));
    result.offsets.reserve(nf+1);
    std::vector<MeshIndex> seen(nf,std::numeric_limits<MeshIndex>::max());
    std::vector<MeshIndex> frontier,next,visited;
    for(std::size_t f=0;f<nf;++f) {
        const auto seed=static_cast<MeshIndex>(f);seen[f]=seed;
        frontier.clear();frontier.push_back(seed);visited.clear();visited.push_back(seed);
        for(TransferIndex step=0;step<rings && !frontier.empty();++step) {
            next.clear();
            for(auto face:frontier) for(auto c=input_.face_offsets[face];c<input_.face_offsets[static_cast<std::size_t>(face)+1];++c) {
                const auto across=face_across(face,tables_.corner_edge[static_cast<std::size_t>(c)]);
                if(across>=0 && seen[static_cast<std::size_t>(across)]!=seed) {
                    seen[static_cast<std::size_t>(across)]=seed;
                    next.push_back(static_cast<MeshIndex>(across));visited.push_back(static_cast<MeshIndex>(across));
                }
            }
            frontier.swap(next);
        }
        std::sort(visited.begin(),visited.end());
        for(auto face:visited) {
            if(result.values.size()>=maximum_values) throw ContractError(ErrorCode::memory_limit,"rings","ring output exceeds memory allowance");
            if(result.values.size()==result.values.capacity()) result.values.reserve(static_cast<std::size_t>(std::min(maximum_values,std::max<std::uint64_t>(1,result.values.capacity()*2ULL))));
            result.values.push_back(face);
        }
        result.offsets.push_back(result.values.size());
    }
    return result;
}
nlohmann::json MeshSnapshot::diagnose() const {
    const auto& g=geometry_;const auto nf=g.face_area.size();std::array<std::uint64_t,5> kinds{};
    std::vector<MeshIndex> first;std::vector<double> areas(g.n_components,0.);std::uint64_t nonplanar=0;
    for(std::size_t f=0;f<nf;++f) {++kinds[g.face_degenerate_kind[f]];areas[g.face_component[f]]+=g.face_area[f];nonplanar+=g.face_nonplanar[f];if(g.face_degenerate_kind[f] && first.size()<20) first.push_back(static_cast<MeshIndex>(f));}
    const auto degenerate=nf-kinds[0];const auto total=detail::numpy_sum(g.face_area); // source area.sum(): pairwise
    return {{"faces",nf},{"degenerate",degenerate},{"degenerate_fraction",nf?static_cast<double>(degenerate)/static_cast<double>(nf):0.},
        {"by_kind",{{"zero_area",kinds[1]},{"twisted",kinds[2]},{"collapsed",kinds[3]},{"sliver",kinds[4]}}},
        {"first_degenerate",first},{"components",g.n_components},{"largest_component_fraction",total>0?*std::max_element(areas.begin(),areas.end())/total:0.},
        {"component_area",g.n_components<=64?nlohmann::json(areas):nlohmann::json(nullptr)},{"nonplanar_ngons",nonplanar}};
}
nlohmann::json MeshSnapshot::summary() const {
    std::uint64_t tris=0,quads=0,ngons=0;const auto nf=geometry_.face_area.size();
    for(std::size_t f=0;f<nf;++f) {const auto n=input_.face_offsets[f+1]-input_.face_offsets[f];tris+=n==3;quads+=n==4;ngons+=n>4;}
    auto count=[](const auto& a){return std::count_if(a.begin(),a.end(),[](const auto& x){return bool(x);});};
    return {{"verts",input_.positions.size()},{"faces",nf},{"edges",tables_.edge_vertices.size()},{"corners",input_.corner_vertices.size()},
        {"quads",quads},{"tris",tris},{"ngons",ngons},{"boundary_edges",count(tables_.edge_boundary)},
        {"nonmanifold_edges",tables_.edge_vertices.size()-static_cast<std::size_t>(count(tables_.edge_manifold))},
        {"nonmanifold_verts",input_.positions.size()-static_cast<std::size_t>(count(tables_.fans))},
        {"regular_verts",count(tables_.vert_regular)},{"degenerate_faces",count(geometry_.face_degenerate_kind)},{"components",geometry_.n_components}};
}
}
