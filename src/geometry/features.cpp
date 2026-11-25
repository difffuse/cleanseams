// SPDX-License-Identifier: GPL-3.0-or-later
#include "cleanseams/fields.hpp"
#include "cleanseams/cache.hpp"
#include "cleanseams/errors.hpp"
#include "cleanseams/detail/mesh_math.hpp"
#include "cleanseams/detail/numpy_compat.hpp"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <numbers>
#include <numeric>
#include <queue>
#include <set>
#include <tuple>
namespace cleanseams {
namespace {
using detail::sub;using detail::dot;using detail::norm;
void require(bool yes,const char* field,const char* message) {if(!yes)throw ContractError(ErrorCode::invalid_input,field,message);}
std::size_t edge_count(const MeshSnapshot& m) {return m.tables().edge_vertices.size();}
void check_mask(const MeshSnapshot& m,std::span<const std::uint8_t> mask,const char* field) {
    require(mask.size()==edge_count(m),field,"edge mask does not match the mesh edges");
}
std::vector<std::uint8_t> barrier_mask(const MeshSnapshot& m,std::optional<std::span<const std::uint8_t>> barriers) {
    std::vector<std::uint8_t> out(m.tables().edge_boundary.begin(),m.tables().edge_boundary.end());
    if(barriers) {check_mask(m,*barriers,"features.barriers");for(std::size_t e=0;e<out.size();++e)out[e]=out[e] || (*barriers)[e];}
    return out;
}
double radians(double degrees) {return degrees*(std::numbers::pi/180.0);}
// Python round(x, 8): correctly rounded decimal, as CPython's dtoa mode 3.
double round8(double x) {
    if(!std::isfinite(x))return x;
    char buffer[512];std::snprintf(buffer,sizeof buffer,"%.8f",x);return std::strtod(buffer,nullptr);
}
double median(std::vector<double> values) {
    std::sort(values.begin(),values.end());const auto n=values.size();
    if(n%2)return values[n/2];
    return (values[n/2-1]+values[n/2])/2.0;
}
struct UnionFind {
    std::vector<MeshIndex> parent;
    explicit UnionFind(std::size_t n):parent(n) {std::iota(parent.begin(),parent.end(),MeshIndex(0));}
    MeshIndex find(MeshIndex x) {while(parent[x]!=x) {parent[x]=parent[parent[x]];x=parent[x];}return x;}
    void unite(MeshIndex a,MeshIndex b) {a=find(a);b=find(b);if(a!=b)parent[b]=a;}
};
std::vector<std::uint64_t> vertex_degree(const MeshSnapshot& m,std::span<const MeshIndex> edges) {
    std::vector<std::uint64_t> degree(m.input().positions.size(),0);
    for(auto e:edges)for(auto v:m.tables().edge_vertices[e])++degree[v];
    return degree;
}
std::vector<MeshIndex> selected(std::span<const std::uint8_t> mask) {
    std::vector<MeshIndex> out;for(std::size_t e=0;e<mask.size();++e)if(mask[e])out.push_back(static_cast<MeshIndex>(e));return out;
}
std::size_t face_size(const MeshSnapshot& m,MeshIndex f) {return static_cast<std::size_t>(m.input().face_offsets[f+1]-m.input().face_offsets[f]);}
}
std::vector<double> wide_dihedral(const MeshSnapshot& m,std::int64_t rings,std::optional<double> radius,
                                  std::optional<std::span<const std::uint8_t>> candidates,const Deadline& deadline) {
    const auto& t=m.tables();const auto& g=m.geometry();const auto ne=edge_count(m);
    const auto neighbourhood=m.face_rings(rings);
    std::vector<double> result(ne);for(std::size_t e=0;e<ne;++e)result[e]=std::abs(g.edge_dihedral[e]);
    std::vector<MeshIndex> edges;
    if(candidates) {check_mask(m,*candidates,"features.candidates");edges=selected(*candidates);}
    else for(std::size_t e=0;e<ne;++e)if(!t.edge_boundary[e])edges.push_back(static_cast<MeshIndex>(e));
    for(std::size_t index=0;index<edges.size();++index) {
        if(index%128==0 && deadline.expired())break;
        const auto edge=edges[index];const auto pair=t.edge_faces[edge];
        if(pair.size()!=2 || !t.edge_consistent[edge])continue;
        const auto& p0=m.input().positions[t.edge_vertices[edge][0]];const auto& p1=m.input().positions[t.edge_vertices[edge][1]];
        const Vec3 midpoint{(p0[0]+p1[0])/2.0,(p0[1]+p1[1])/2.0,(p0[2]+p1[2])/2.0};
        auto tangent=sub(p1,p0);tangent=detail::divide(tangent,std::max(norm(tangent),1e-300));
        auto across=sub(g.face_center[pair[1]],g.face_center[pair[0]]);
        const double along=dot(across,tangent);across=sub(across,detail::mul(tangent,along));
        std::array<Vec3,2> normals{};
        for(std::size_t side=0;side<2;++side) {
            const auto face=pair[side];const double sign=side?1.0:-1.0;
            const double slack=-1e-12*std::max(norm(across),1e-12);
            std::vector<MeshIndex> nearby;
            for(auto f:neighbourhood[face]) {
                bool keep=dot(sub(g.face_center[f],midpoint),across)*sign>=slack;
                if(radius)keep=keep && norm(sub(g.face_center[f],g.face_center[face]))<=*radius;
                if(keep)nearby.push_back(f);
            }
            if(nearby.empty())nearby.push_back(face);
            Vec3 normal{0,0,0};
            for(auto f:nearby)for(std::size_t d=0;d<3;++d)normal[d]+=g.face_normal[f][d]*g.face_area[f];
            normals[side]=detail::divide(normal,std::max(norm(normal),1e-300));
        }
        result[edge]=std::max(result[edge],std::acos(std::clamp(dot(normals[0],normals[1]),-1.0,1.0)));
    }
    return result;
}
std::vector<MeshIndex> parallel_edges(const MeshSnapshot& m,MeshIndex edge,std::int64_t steps) {
    const auto& t=m.tables();require(edge<edge_count(m),"features.edge","edge index out of range");
    std::set<MeshIndex> seen{edge};
    for(auto start:t.edge_faces[edge]) {
        TransferIndex face=start;MeshIndex current=edge;
        for(std::int64_t s=0;s<steps;++s) {
            if(face_size(m,static_cast<MeshIndex>(face))!=4)break;
            std::array<MeshIndex,4> edges{};std::size_t k=0;
            for(auto c:m.face_corners(static_cast<MeshIndex>(face)))edges[k++]=t.corner_edge[c];
            const auto at=static_cast<std::size_t>(std::find(edges.begin(),edges.end(),current)-edges.begin());
            const auto opposite=edges[(at+2)%4];
            if(seen.contains(opposite))break;
            seen.insert(opposite);
            face=m.face_across(static_cast<MeshIndex>(face),opposite);current=opposite;
            if(face<0)break;
        }
    }
    return {seen.begin(),seen.end()};
}
std::vector<std::uint8_t> suppress_parallel(const MeshSnapshot& m,std::span<const std::uint8_t> candidates,std::span<const double> strength,
                                            double ratio,std::int64_t rings) {
    check_mask(m,candidates,"features.candidates");require(strength.size()==candidates.size(),"features.strength","strength does not match the mesh edges");
    const auto& g=m.geometry();std::vector<std::uint8_t> keep(candidates.begin(),candidates.end());
    auto local=[&](MeshIndex e) {return std::abs(g.edge_dihedral[e]);};
    for(auto edge:selected(candidates)) {
        const auto parallel=parallel_edges(m,edge,rings);
        double maximum=-INFINITY;std::vector<double> values;
        for(auto p:parallel) {maximum=std::max(maximum,local(p));values.push_back(local(p));}
        if(parallel.size()>=3 && local(edge)<1.1*median(values)) {keep[edge]=0;continue;}
        if(local(edge)+1e-12<ratio*maximum) {keep[edge]=0;continue;}
        // Prefer the centre of the broad crease; an exact plateau has one stable winner.
        std::optional<MeshIndex> best;std::tuple<double,double,std::int64_t> best_key{};
        for(auto p:parallel) {
            if(!candidates[p] || !(local(p)>=maximum*ratio))continue;
            const std::tuple<double,double,std::int64_t> key{round8(strength[p]),round8(local(p)),-static_cast<std::int64_t>(p)};
            if(!best || key>best_key) {best=p;best_key=key;}
        }
        require(best.has_value(),"features.suppress","no parallel peer");
        const auto b=*best;
        if(strength[edge]+1e-8<strength[b] || (std::abs(strength[edge]-strength[b])<=1e-8 && std::abs(local(edge)-local(b))<=1e-8 && edge!=b))
            keep[edge]=0;
    }
    return keep;
}
std::vector<std::vector<MeshIndex>> feature_chains(const MeshSnapshot& m,std::span<const std::uint8_t> mask) {
    check_mask(m,mask,"features.mask");const auto& t=m.tables();
    UnionFind uf(m.input().positions.size());const auto edges=selected(mask);
    for(auto e:edges)uf.unite(t.edge_vertices[e][0],t.edge_vertices[e][1]);
    std::vector<std::vector<MeshIndex>> out;std::vector<std::int64_t> slot(m.input().positions.size(),-1);
    for(auto e:edges) {
        const auto root=uf.find(t.edge_vertices[e][0]);
        if(slot[root]<0) {slot[root]=static_cast<std::int64_t>(out.size());out.emplace_back();}
        out[static_cast<std::size_t>(slot[root])].push_back(e);
    }
    return out;
}
TerminatedChains terminate_chains(const MeshSnapshot& m,std::span<const std::uint8_t> mask,std::int64_t min_chain,std::int64_t max_extend,
                                  std::optional<std::span<const std::uint8_t>> barriers_in) {
    check_mask(m,mask,"features.mask");const auto& t=m.tables();
    std::vector<std::uint8_t> result(mask.begin(),mask.end());const auto barriers=barrier_mask(m,barriers_in);
    for(const auto& edges:feature_chains(m,result))if(static_cast<std::int64_t>(edges.size())<min_chain)for(auto e:edges)result[e]=0;
    const auto degree=vertex_degree(m,selected(result));
    auto touches=[&](MeshIndex v,const auto& predicate) {for(auto e:t.vert_edges[v])if(predicate(e))return true;return false;};
    for(MeshIndex vertex=0;vertex<degree.size();++vertex) {
        if(degree[vertex]!=1)continue;
        if(touches(vertex,[&](MeshIndex e) {return barriers[e]!=0;}))continue;
        MeshIndex incoming=0;for(auto e:t.vert_edges[vertex])if(result[e]) {incoming=e;break;}
        std::vector<MeshIndex> path;std::set<MeshIndex> seen{vertex};MeshIndex current=vertex;
        for(std::int64_t s=0;s<max_extend;++s) {
            const auto& fan=t.fans[current];const TransferIndex edge=fan?fan->straight(incoming):-1;
            if(edge<0 || result[static_cast<std::size_t>(edge)] || barriers[static_cast<std::size_t>(edge)])break;
            path.push_back(static_cast<MeshIndex>(edge));current=m.edge_other(static_cast<MeshIndex>(edge),current);
            if(seen.contains(current))break;
            seen.insert(current);
            if(touches(current,[&](MeshIndex e) {return result[e] || barriers[e];})) {for(auto p:path)result[p]=1;break;}
            incoming=static_cast<MeshIndex>(edge);
        }
    }
    TerminatedChains out;
    for(const auto& edges:feature_chains(m,result)) {
        const auto d=vertex_degree(m,edges);bool free=false;
        for(MeshIndex v=0;v<d.size() && !free;++v)if(d[v]==1 && !touches(v,[&](MeshIndex e) {return barriers[e]!=0;}))free=true;
        if(free && static_cast<std::int64_t>(edges.size())<2*min_chain)for(auto e:edges)result[e]=0;
        else if(free)out.open_edges.insert(out.open_edges.end(),edges.begin(),edges.end());
    }
    std::sort(out.open_edges.begin(),out.open_edges.end());out.mask=std::move(result);
    return out;
}
BevelReassignment reassign_bevel_bands(const std::shared_ptr<const MeshSnapshot>& mesh,std::span<const std::uint8_t> cut,
                                       std::optional<std::span<const std::uint8_t>> barriers,bool allow_holes) {
    require(bool(mesh),"features.mesh","missing mesh snapshot");const auto& m=*mesh;const auto& t=m.tables();const auto& g=m.geometry();
    check_mask(m,cut,"features.cut");const auto ne=edge_count(m);const auto nf=m.input().face_offsets.size()-1;
    const auto protect=barrier_mask(m,barriers);
    std::vector<std::uint8_t> initial(ne);for(std::size_t e=0;e<ne;++e)initial[e]=cut[e] || protect[e];
    const auto cs_owner=cut_surface(mesh,initial);const CutSurface& cs=*cs_owner;const auto& ct=cs.tables();const auto nc=ct.n_charts;
    std::vector<double> areas(nc,0.0);std::vector<Vec3> normals(nc,Vec3{0,0,0});
    for(std::size_t f=0;f<nf;++f) {
        const auto ch=ct.chart_of_face[f];areas[ch]+=g.face_area[f];
        for(std::size_t d=0;d<3;++d)normals[ch][d]+=g.face_normal[f][d]*g.face_area[f];
    }
    std::vector<double> coherence(nc);
    for(std::size_t ch=0;ch<nc;++ch) {
        for(auto& x:normals[ch])x/=std::max(areas[ch],1e-300);
        coherence[ch]=norm(normals[ch]);
    }
    const double total=detail::numpy_sum(areas);
    std::vector<std::uint8_t> flat(nc);
    for(std::size_t ch=0;ch<nc;++ch) {
        const bool topology=allow_holes?(!ct.chart_unsupported[ch] && ct.chart_euler[ch]==2-static_cast<std::int64_t>(ct.chart_n_loops[ch]) && ct.chart_n_loops[ch]>0)
                                       :ct.chart_is_disk[ch]!=0;
        flat[ch]=coherence[ch]>=.98 && areas[ch]>=.01*total && topology;
        for(auto& x:normals[ch])x/=std::max(coherence[ch],1e-300);
    }
    std::vector<std::int64_t> labels(ct.chart_of_face.begin(),ct.chart_of_face.end());
    std::vector<MeshIndex> changed;std::set<std::int64_t> panels;
    auto cut_of=[&](const std::vector<std::int64_t>& proposal) {
        auto out=protect;
        for(std::size_t e=0;e<ne;++e) {const auto faces=t.edge_faces[e];if(faces.size()==2 && proposal[faces[0]]!=proposal[faces[1]])out[e]=1;}
        return out;
    };
    for(MeshIndex chart=0;chart<nc;++chart) {
        if(flat[chart] || areas[chart]>(allow_holes?.25:.12)*total)continue;
        std::set<std::int64_t> neighbours;
        for(const auto& a:ct.adjacency) {
            if(a.charts[0]==a.charts[1] || (a.charts[0]!=chart && a.charts[1]!=chart))continue;
            neighbours.insert(a.charts[0]==chart?a.charts[1]:a.charts[0]);
        }
        if(neighbours.size()<3 || !std::all_of(neighbours.begin(),neighbours.end(),[&](auto n) {return flat[static_cast<std::size_t>(n)]!=0;}))continue;
        auto topology_ok=[&](const CutSurface& surface,const std::vector<std::int64_t>& proposal) {
            const auto& st=surface.tables();std::set<MeshIndex> affected;
            for(std::size_t f=0;f<nf;++f)if(neighbours.contains(proposal[f]))affected.insert(st.chart_of_face[f]);
            if(affected.size()!=neighbours.size())return false;
            if(!allow_holes)return std::all_of(affected.begin(),affected.end(),[&](auto c) {return st.chart_is_disk[c]!=0;});
            for(auto owner:neighbours) {
                std::set<MeshIndex> target;for(std::size_t f=0;f<nf;++f)if(proposal[f]==owner)target.insert(st.chart_of_face[f]);
                if(target.size()!=1)return false;
                const auto part=*target.begin();
                if(st.chart_unsupported[part] || st.chart_n_loops[part]>ct.chart_n_loops[static_cast<std::size_t>(owner)] ||
                   st.chart_euler[part]!=2-static_cast<std::int64_t>(st.chart_n_loops[part]))return false;
            }
            return true;
        };
        const auto faces=ct.chart_faces[chart];
        std::vector<std::uint8_t> chosen(nf,0);for(auto f:faces)chosen[f]=1;
        auto proposed=labels;for(auto f:faces)proposed[f]=-1;
        std::vector<double> distance(nf,INFINITY);
        using Item=std::tuple<double,std::int64_t,std::int64_t>;
        std::priority_queue<Item,std::vector<Item>,std::greater<Item>> queue;
        auto push=[&](MeshIndex face,std::int64_t owner,MeshIndex source,double cost) {
            const double normal_error=std::max(0.0,1-dot(g.face_normal[face],normals[static_cast<std::size_t>(owner)]));
            const double travel=norm(sub(g.face_center[face],g.face_center[source]));
            const double value=cost+travel*(1+8*normal_error);
            if(value<distance[face]) {distance[face]=value;queue.push({value,owner,face});}
        };
        for(auto face:faces)for(auto corner:m.face_corners(face)) {
            const auto edge=t.corner_edge[corner];if(protect[edge])continue;
            const auto other=m.face_across(face,edge);
            if(other>=0 && neighbours.contains(labels[static_cast<std::size_t>(other)]))push(face,labels[static_cast<std::size_t>(other)],static_cast<MeshIndex>(other),0.0);
        }
        while(!queue.empty()) {
            const auto [value,owner,f]=queue.top();queue.pop();const auto face=static_cast<MeshIndex>(f);
            if(proposed[face]>=0 || value>distance[face])continue;
            proposed[face]=owner;
            for(auto corner:m.face_corners(face)) {
                const auto edge=t.corner_edge[corner];if(protect[edge])continue;
                const auto other=m.face_across(face,edge);
                if(other>=0 && chosen[static_cast<std::size_t>(other)] && proposed[static_cast<std::size_t>(other)]<0)push(static_cast<MeshIndex>(other),owner,face,value);
            }
        }
        if(std::any_of(faces.begin(),faces.end(),[&](auto f) {return proposed[f]<0;}))continue;
        // A constant-curvature strip is equally close to its two panels: resolve the
        // whole strip by panel normal (rounded to four decimals, first maximum).
        const std::vector<std::int64_t> owners(neighbours.begin(),neighbours.end());
        auto normal_proposal=labels;
        for(auto face:faces) {
            std::size_t best=0;double best_score=-INFINITY;
            for(std::size_t k=0;k<owners.size();++k) {
                const double score=std::nearbyint(dot(g.face_normal[face],normals[static_cast<std::size_t>(owners[k])])*1e4)/1e4;
                if(score>best_score) {best_score=score;best=k;}
            }
            normal_proposal[face]=owners[best];
        }
        if(topology_ok(*cut_surface(mesh,cut_of(normal_proposal)),normal_proposal))proposed=normal_proposal;
        if(!topology_ok(*cut_surface(mesh,cut_of(proposed)),proposed))continue;
        labels=std::move(proposed);changed.push_back(chart);panels.insert(neighbours.begin(),neighbours.end());
    }
    BevelReassignment out;out.changed_charts=changed;
    if(changed.empty()) {out.cut.assign(cut.begin(),cut.end());return out;}
    out.cut=cut_of(labels);
    for(std::size_t e=0;e<ne;++e)if(t.edge_boundary[e])out.cut[e]=0;
    for(std::size_t f=0;f<nf;++f)if(panels.contains(labels[f]))out.panel_faces.push_back(static_cast<MeshIndex>(f));
    return out;
}
FeatureField detect_features(const std::shared_ptr<const MeshSnapshot>& mesh,const FeatureOptions& o,
                             std::optional<std::span<const std::uint8_t>> barriers,const Deadline& deadline) {
    require(bool(mesh),"features.mesh","missing mesh snapshot");const auto& m=*mesh;const auto& t=m.tables();const auto& g=m.geometry();
    // A negative min_chain or max_extend acts as 0, as in the source (plain comparisons and
    // loops that run no step); a negative ring count raises in the source's face_rings.
    require(std::isfinite(o.angle_deg) && std::isfinite(o.ratio) && o.rings>=0,"features.options","invalid feature options");
    const auto ne=edge_count(m);const auto nf=m.input().face_offsets.size()-1;FeatureField out;
    if(o.angle_deg<=0 || deadline.expired()) {out.cut.assign(ne,0);return out;}
    const double threshold=radians(o.angle_deg);
    std::vector<double> local(ne);for(std::size_t e=0;e<ne;++e)local[e]=std::abs(g.edge_dihedral[e]);
    std::vector<std::uint8_t> candidate(ne);
    for(std::size_t e=0;e<ne;++e)candidate[e]=local[e]>=threshold/static_cast<double>(std::max<std::int64_t>(2*o.rings+1,1))*.4 && !t.edge_boundary[e];
    std::size_t quads=0,ngons=0;for(MeshIndex f=0;f<nf;++f) {quads+=face_size(m,f)==4;ngons+=face_size(m,f)>4;}
    std::optional<double> radius;
    if(!(static_cast<double>(quads)/static_cast<double>(nf)>=.9)) {
        Vec3 low{INFINITY,INFINITY,INFINITY},high{-INFINITY,-INFINITY,-INFINITY};
        for(const auto& p:m.input().positions)for(std::size_t d=0;d<3;++d) {low[d]=std::min(low[d],p[d]);high[d]=std::max(high[d],p[d]);}
        radius=.015*norm(sub(high,low));
    }
    // Large n-gons span several sides of a CAD part: no broad fitting there.
    if(static_cast<double>(ngons)/static_cast<double>(nf)>.02)out.wide_dihedral=local;
    else {out.wide_dihedral=wide_dihedral(m,o.rings,radius,std::span<const std::uint8_t>(candidate),deadline);out.wide_measurement=true;}
    const auto& strength=out.wide_dihedral;
    std::vector<std::uint8_t> strong(ne);
    for(std::size_t e=0;e<ne;++e) {strong[e]=candidate[e] && strength[e]>=threshold;candidate[e]=candidate[e] && strength[e]>=.9*threshold;}
    const auto components=feature_chains(m,candidate);
    auto any_strong=[&](const std::vector<MeshIndex>& c) {return std::any_of(c.begin(),c.end(),[&](auto e) {return strong[e]!=0;});};
    for(const auto& c:components)if(!any_strong(c))for(auto e:c)candidate[e]=0;
    auto thinned=suppress_parallel(m,candidate,strength,o.ratio,o.rings);out.thinned=thinned;
    // A broad crease may end in a plateau around a rounded corner: keep complete rims.
    for(const auto& c:components) {
        if(!any_strong(c))continue;
        const auto degree=vertex_degree(m,c);
        if(std::all_of(degree.begin(),degree.end(),[](auto d) {return d==0 || d==2;}))for(auto e:c)thinned[e]=1;
    }
    candidate=std::move(thinned);
    for(std::size_t e=0;e<ne;++e)if(local[e]>=threshold && !t.edge_boundary[e])candidate[e]=1;
    for(const auto& c:feature_chains(m,candidate)) {
        std::vector<double> positive_lengths,lengths;
        for(auto e:c) {lengths.push_back(g.edge_length[e]);if(g.edge_dihedral[e]>=0)positive_lengths.push_back(g.edge_length[e]);}
        const double positive=detail::numpy_sum(positive_lengths),total=detail::numpy_sum(lengths);
        if(std::max(positive,total-positive)<.6*total)for(auto e:c)candidate[e]=0;
    }
    auto terminated=terminate_chains(m,candidate,o.min_chain,o.max_extend,barriers);
    out.terminated=terminated.mask;out.terminated_open=terminated.open_edges;
    out.cut=std::move(terminated.mask);out.open_features=std::move(terminated.open_edges);
    if(!deadline.expired()) {
        auto reassigned=reassign_bevel_bands(mesh,out.cut,barriers);
        out.reassigned_charts=reassigned.changed_charts;
        out.cut=std::move(reassigned.cut);out.panel_faces=std::move(reassigned.panel_faces);
        if(!out.reassigned_charts.empty()) {
            const auto degree=vertex_degree(m,selected(out.cut));out.open_features.clear();
            for(auto e:selected(out.cut))if(degree[t.edge_vertices[e][0]]==1 || degree[t.edge_vertices[e][1]]==1)out.open_features.push_back(e);
        }
    }
    out.chains=feature_chains(m,out.cut);out.measured=true;
    return out;
}
}
