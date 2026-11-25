// SPDX-License-Identifier: GPL-3.0-or-later
#include "cleanseams/layout.hpp"
#include "cleanseams/cache.hpp"
#include "cleanseams/panels.hpp"
#include "cleanseams/errors.hpp"
#include "cleanseams/detail/mesh_math.hpp"
#include "cleanseams/detail/numpy_compat.hpp"
#include "cleanseams/detail/pyset.hpp"
#include <algorithm>
#include <chrono>
#include <map>
#include <cmath>
#include <numeric>
#include <queue>
#include <set>
#include <tuple>
namespace cleanseams {
namespace {
using detail::sub;using detail::dot;using detail::norm;
void require(bool yes,const char* field,const char* message) {if(!yes)throw ContractError(ErrorCode::invalid_input,field,message);}
std::vector<MeshIndex> ones(std::span<const std::uint8_t> mask) {
    std::vector<MeshIndex> out;for(std::size_t i=0;i<mask.size();++i)if(mask[i])out.push_back(static_cast<MeshIndex>(i));return out;
}
void emit(LayoutTrace* trace,nlohmann::json event) {if(trace)trace->events.push_back(std::move(event));}
bool disk(const CutTables& t,MeshIndex ch) {return t.chart_n_loops[ch]==1 && t.chart_euler[ch]==1 && !t.chart_pinched[ch];}
}
DualGraph dual_graph(const MeshSnapshot& m,std::optional<std::span<const std::uint8_t>> barriers,double w_crease) {
    const auto& t=m.tables();const auto& g=m.geometry();const auto ne=t.edge_vertices.size(),nf=m.input().face_offsets.size()-1;
    if(barriers)require(barriers->size()==ne,"segment.barriers","barrier mask does not match the mesh edges");
    std::vector<MeshIndex> edge_ids,a,b;
    for(MeshIndex e=0;e<ne;++e) {
        const auto pair=t.edge_faces[e];
        if(pair.size()==2 && t.edge_consistent[e] && !(barriers && (*barriers)[e])) {edge_ids.push_back(e);a.push_back(pair[0]);b.push_back(pair[1]);}
    }
    const double floor=std::max(m.mean_edge_length()*1e-6,1e-300);
    std::vector<double> weight(edge_ids.size());
    for(std::size_t i=0;i<edge_ids.size();++i) {
        const double length=std::max(norm(sub(g.face_center[a[i]],g.face_center[b[i]])),floor);
        const double bend=1-std::clamp(detail::numpy_einsum_dot3(g.face_normal[a[i]],g.face_normal[b[i]]),-1.0,1.0);
        weight[i]=length*(1+w_crease*bend);
    }
    // Rows then columns, stable: np.lexsort((cols, rows)) over r_[a, b] / r_[b, a].
    const auto k=edge_ids.size();std::vector<std::size_t> order(2*k);std::iota(order.begin(),order.end(),0);
    auto row=[&](std::size_t i) {return i<k?a[i]:b[i-k];};auto col=[&](std::size_t i) {return i<k?b[i]:a[i-k];};
    std::stable_sort(order.begin(),order.end(),[&](auto x,auto y) {return std::pair(row(x),col(x))<std::pair(row(y),col(y));});
    DualGraph out;out.indptr.assign(nf+1,0);
    for(std::size_t i=0;i<2*k;++i)++out.indptr[row(i)+1];
    for(std::size_t f=0;f<nf;++f)out.indptr[f+1]+=out.indptr[f];
    for(auto i:order) {out.indices.push_back(col(i));out.weights.push_back(weight[i%k]);out.edges.push_back(edge_ids[i%k]);}
    return out;
}
void update_distance(const DualGraph& graph,MeshIndex seed,std::vector<double>& distance) {
    distance[seed]=0.0;
    using Item=std::pair<double,MeshIndex>;std::priority_queue<Item,std::vector<Item>,std::greater<Item>> queue;queue.push({0.0,seed});
    while(!queue.empty()) {
        const auto [value,face]=queue.top();queue.pop();
        if(value>distance[face])continue;
        for(auto j=graph.indptr[face];j<graph.indptr[face+1];++j) {
            const auto other=graph.indices[j];const double candidate=value+graph.weights[j];
            if(candidate<distance[other]) {distance[other]=candidate;queue.push({candidate,other});}
        }
    }
}
std::vector<MeshIndex> farthest_point_seeds(const DualGraph& graph,std::uint64_t count,std::optional<std::vector<MeshIndex>> initial,const Deadline& deadline,LayoutTrace* trace) {
    std::vector<MeshIndex> seeds;
    if(graph.size()) {
        std::vector<double> distance(graph.size(),INFINITY);
        std::vector<MeshIndex> first=initial?*initial:std::vector<MeshIndex>{0};std::set<MeshIndex> seen;
        for(auto face:first) {
            if(!seen.insert(face).second)continue;
            require(face<graph.size(),"segment.seeds","seed outside the face array");
            seeds.push_back(face);update_distance(graph,face,distance);
        }
        count=std::min<std::uint64_t>(graph.size(),std::max<std::uint64_t>(1,count));const auto minimum=std::min<std::uint64_t>(16,count);
        while(true) {
            const auto face=static_cast<MeshIndex>(std::max_element(distance.begin(),distance.end())-distance.begin());
            const bool uncovered=!std::isfinite(distance[face]);
            if(!uncovered && (seeds.size()>=count || (seeds.size()>=minimum && deadline.expired()) || distance[face]<=0))break;
            seeds.push_back(face);update_distance(graph,face,distance);
        }
    }
    emit(trace,{{"stage","farthest_point_seeds"},{"result",seeds}});
    return seeds;
}
std::vector<std::int64_t> grow_labels(const DualGraph& graph,std::span<const MeshIndex> seeds,std::optional<std::span<const Vec3>> normals,
                                      std::optional<std::span<const Vec3>> proxies,double w_normal,LayoutTrace* trace) {
    std::set<MeshIndex> distinct(seeds.begin(),seeds.end());
    require(distinct.size()==seeds.size() && std::all_of(seeds.begin(),seeds.end(),[&](auto s) {return s<graph.size();}),"segment.seeds","seeds must be distinct valid face ids");
    std::vector<std::int64_t> labels(graph.size(),-1);std::vector<double> distance(graph.size(),INFINITY);
    for(std::size_t i=0;i<seeds.size();++i) {labels[seeds[i]]=static_cast<std::int64_t>(i);distance[seeds[i]]=0.0;}
    using Item=std::tuple<double,std::int64_t,MeshIndex>;std::priority_queue<Item,std::vector<Item>,std::greater<Item>> queue;
    auto expand=[&](MeshIndex face,std::int64_t label,double value) {
        for(auto j=graph.indptr[face];j<graph.indptr[face+1];++j) {
            const auto other=graph.indices[j];if(labels[other]>=0)continue;
            double penalty=1.0;
            if(proxies)penalty+=w_normal*std::pow(1-std::clamp(dot((*normals)[other],(*proxies)[static_cast<std::size_t>(label)]),-1.0,1.0),2.0);
            const double candidate=value+graph.weights[j]*penalty;
            if(candidate<distance[other]) {distance[other]=candidate;queue.push({candidate,label,other});}
        }
    };
    for(std::size_t i=0;i<seeds.size();++i)expand(seeds[i],static_cast<std::int64_t>(i),0.0);
    while(!queue.empty()) {
        const auto [value,label,face]=queue.top();queue.pop();
        if(labels[face]>=0 || value>distance[face])continue;
        labels[face]=label;expand(face,label,value);
    }
    require(std::all_of(labels.begin(),labels.end(),[](auto l) {return l>=0;}),"segment.seeds","every dual-graph component needs a seed");
    emit(trace,{{"stage","grow_labels"},{"labels",labels}});
    return labels;
}
LloydResult lloyd(const MeshSnapshot& m,const DualGraph& graph,std::vector<MeshIndex> seeds,std::int64_t steps,const Deadline& deadline,LayoutTrace* trace) {
    const auto& g=m.geometry();LloydResult out;out.labels=grow_labels(graph,seeds,std::nullopt,std::nullopt,2.0,trace);
    for(std::int64_t step=0;step<std::max<std::int64_t>(0,steps);++step) {
        if(deadline.expired())break;
        const auto k=seeds.size();std::vector<double> area(k,0.0);std::vector<Vec3> normals(k,Vec3{0,0,0}),centres(k,Vec3{0,0,0});
        for(std::size_t f=0;f<out.labels.size();++f) {
            const auto l=static_cast<std::size_t>(out.labels[f]);area[l]+=g.face_area[f];
            for(std::size_t d=0;d<3;++d) {normals[l][d]+=g.face_normal[f][d]*g.face_area[f];centres[l][d]+=g.face_center[f][d]*g.face_area[f];}
        }
        for(std::size_t l=0;l<k;++l) {
            const double n=std::max(norm(normals[l]),1e-300);for(auto& x:normals[l])x/=n;
            for(auto& x:centres[l])x/=std::max(area[l],1e-300);
        }
        std::vector<MeshIndex> next;
        for(std::size_t l=0;l<k;++l) {
            std::vector<MeshIndex> ids;for(std::size_t f=0;f<out.labels.size();++f)if(out.labels[f]==static_cast<std::int64_t>(l))ids.push_back(static_cast<MeshIndex>(f));
            std::vector<double> d2(ids.size());double biggest=-INFINITY;
            for(std::size_t i=0;i<ids.size();++i) {
                const auto dx=sub(g.face_center[ids[i]],centres[l]);d2[i]=(dx[0]*dx[0]+dx[1]*dx[1])+dx[2]*dx[2];biggest=std::max(biggest,d2[i]);
            }
            std::size_t best=0;double best_score=INFINITY;
            for(std::size_t i=0;i<ids.size();++i) {
                const double angular=1-std::clamp(dot(g.face_normal[ids[i]],normals[l]),-1.0,1.0);
                const double score=d2[i]/std::max(biggest,1e-300)+.25*(angular*angular);
                if(score<best_score) {best_score=score;best=i;}
            }
            next.push_back(ids[best]);
        }
        seeds=std::move(next);
        auto proposed=grow_labels(graph,seeds,std::span<const Vec3>(g.face_normal),std::span<const Vec3>(normals),2.0,trace);
        if(proposed==out.labels)break;
        out.labels=std::move(proposed);
    }
    out.seeds=seeds;
    emit(trace,{{"stage","lloyd"},{"labels",out.labels},{"seeds",out.seeds}});
    return out;
}
EdgeMask regions_to_cut(const MeshSnapshot& m,std::span<const std::int64_t> labels,std::optional<std::span<const std::uint8_t>> initial) {
    const auto& t=m.tables();const auto nf=m.input().face_offsets.size()-1;
    require(labels.size()==nf && std::all_of(labels.begin(),labels.end(),[](auto l) {return l>=0;}),"segment.labels","regions must cover every face");
    EdgeMask cut(t.edge_boundary.begin(),t.edge_boundary.end());
    if(initial)for(std::size_t e=0;e<cut.size();++e)cut[e]=cut[e] || (*initial)[e];
    for(MeshIndex e=0;e<cut.size();++e) {
        const auto pair=t.edge_faces[e];
        if(pair.size()!=2 || !t.edge_consistent[e] || labels[pair[0]]!=labels[pair[1]])cut[e]=1;
    }
    return cut;
}
DiskLayout make_region_disks(const std::shared_ptr<const MeshSnapshot>& mesh,std::span<const std::uint8_t> cut,CostField cost,const Deadline& deadline,LayoutTrace* trace) {
    auto cs=cut_surface(mesh,cut);
    auto all_disks=[](const CutSurface& s) {for(MeshIndex c=0;c<s.tables().n_charts;++c)if(!disk(s.tables(),c))return false;return true;};
    auto loops_total=[](const CutSurface& s) {std::uint64_t n=0;for(auto x:s.tables().chart_n_loops)n+=x;return n;};
    while(!all_disks(*cs) && !deadline.expired()) {
        bool changed=false;
        for(MeshIndex chart=0;chart<cs->tables().n_charts && !changed;++chart) {
            if(disk(cs->tables(),chart) || cs->tables().chart_loops[chart].size()<=1)continue;
            const auto path=join_boundary_loops(*cs,chart,cost);
            if(!path || path->empty())continue;
            auto trial=cs->tables().cut;for(auto e:*path)trial[e]=1;
            auto candidate=cut_surface(mesh,trial);
            if(candidate->tables().n_charts==cs->tables().n_charts && loops_total(*candidate)<loops_total(*cs)) {cs=candidate;changed=true;}
        }
        if(!changed)break;
    }
    DiskLayout out;
    if(all_disks(*cs)) {out.cut=cs->tables().cut;out.cs=cs;}
    else {
        // Glue polygonal disks only along a dual spanning tree: a disk, and handles open.
        const auto graph=dual_graph(*mesh,std::span<const std::uint8_t>(cs->tables().cut));auto result=cs->tables().cut;
        for(MeshIndex chart=0;chart<cs->tables().n_charts;++chart) {
            if(disk(cs->tables(),chart))continue;
            const auto faces=cs->tables().chart_faces[chart];
            std::set<MeshIndex> seen{faces[0]};std::vector<MeshIndex> queue{faces[0]};std::set<MeshIndex> tree;
            for(std::size_t q=0;q<queue.size();++q) {
                const auto face=queue[q];
                for(auto j=graph.indptr[face];j<graph.indptr[face+1];++j) {
                    const auto other=graph.indices[j];const auto edge=graph.edges[j];result[edge]=1;
                    if(seen.insert(other).second) {queue.push_back(other);tree.insert(edge);}
                }
            }
            for(auto e:tree)result[e]=0;
        }
        out.cs=cut_surface(mesh,result);out.cut=std::move(result);
    }
    emit(trace,{{"stage","make_region_disks"},{"cut",ones(out.cut)},{"charts",out.cs->tables().n_charts}});
    return out;
}
Segmentation segment_mesh(const std::shared_ptr<const MeshSnapshot>& mesh,std::optional<std::span<const std::uint8_t>> initial_cut,std::int64_t seed_density,
                          std::int64_t max_patches,std::int64_t steps,CostField cost,const Deadline& deadline,std::optional<std::vector<MeshIndex>> seeds,LayoutTrace* trace) {
    if(seed_density<=0 || max_patches<1)throw ContractError(ErrorCode::invalid_option,"layout.seed_density","seed density and patch cap must be positive");
    const auto started=std::chrono::steady_clock::now();const auto nf=mesh->input().face_offsets.size()-1;
    const auto graph=dual_graph(*mesh,initial_cut);
    const auto wanted=static_cast<std::int64_t>(std::ceil(static_cast<double>(nf)/static_cast<double>(seed_density)));
    const auto count=std::min<std::int64_t>(max_patches,std::max<std::int64_t>(16,std::min<std::int64_t>(256,wanted)));
    auto selected=farthest_point_seeds(graph,static_cast<std::uint64_t>(count),seeds,deadline,trace);
    auto grown=lloyd(*mesh,graph,selected,steps,deadline,trace);
    const auto cut=regions_to_cut(*mesh,grown.labels,initial_cut);
    auto disks=make_region_disks(mesh,cut,cost,deadline,trace);
    Segmentation out;out.cut=std::move(disks.cut);out.cs=std::move(disks.cs);out.seeds=grown.seeds.size();out.regions=out.cs->tables().n_charts;
    out.dijkstra_s=std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count();
    emit(trace,{{"stage","segment_mesh"},{"cut",ones(out.cut)},{"seeds",out.seeds},{"regions",out.regions}});
    return out;
}
DiskLayout split_regions(const std::shared_ptr<const MeshSnapshot>& mesh,const CutSurface& surface,std::span<const MeshIndex> charts,CostField cost,const Deadline& deadline,LayoutTrace* trace) {
    const auto& st=surface.tables();const auto graph=dual_graph(*mesh,std::span<const std::uint8_t>(st.cut));
    detail::PyIntSet wanted;for(auto c:charts)wanted.add(c);
    std::vector<MeshIndex> seeds;for(MeshIndex c=0;c<st.n_charts;++c)seeds.push_back(st.chart_faces[c][0]);
    std::vector<double> distance(mesh->input().face_offsets.size()-1,INFINITY);
    for(auto key:wanted.order()) {
        const auto chart=static_cast<MeshIndex>(key);const auto faces=st.chart_faces[chart];
        if(faces.size()<2 || deadline.expired())continue;
        update_distance(graph,seeds[chart],distance);
        auto farthest=[&] {MeshIndex best=faces[0];double d=-INFINITY;for(auto f:faces)if(distance[f]>d) {d=distance[f];best=f;}return best;};
        const auto first=farthest();for(auto f:faces)distance[f]=INFINITY;
        update_distance(graph,first,distance);const auto second=farthest();
        seeds[chart]=first;if(first!=second)seeds.push_back(second);
    }
    const auto labels=grow_labels(graph,seeds,std::nullopt,std::nullopt,2.0,trace);
    const auto cut=regions_to_cut(*mesh,labels,std::span<const std::uint8_t>(st.cut));
    auto out=make_region_disks(mesh,cut,cost,deadline,trace);
    emit(trace,{{"stage","split_regions"},{"charts",std::vector<MeshIndex>(charts.begin(),charts.end())},{"cut",ones(out.cut)},{"charts_after",out.cs->tables().n_charts}});
    return out;
}
// ------------------------------------------------------------ preference-aware regions (T34)
// enhanced_v1 only: nothing above calls these, so legacy_010 and native_baseline growth is
// unchanged.
MeshIndex RegionDual::node(MeshIndex face) const {
    const auto it=std::lower_bound(faces.begin(),faces.end(),face);
    if(it==faces.end() || *it!=face)throw ContractError(ErrorCode::invalid_input,"regions.face","face outside the region");
    return node_of_face[static_cast<std::size_t>(it-faces.begin())];
}
RegionDual region_dual(const MeshSnapshot& m,std::span<const MeshIndex> faces,std::span<const std::uint8_t> barriers,std::span<const std::uint8_t> weld,
                       std::span<const double> strength,double w_normal,double w_preference) {
    const auto& t=m.tables();const auto& g=m.geometry();const auto ne=t.edge_vertices.size(),nf=m.input().face_offsets.size()-1;
    require(barriers.empty() || barriers.size()==ne,"regions.barriers","barrier mask does not match the mesh edges");
    require(weld.empty() || weld.size()==ne,"regions.weld","weld mask does not match the mesh edges");
    require(strength.empty() || strength.size()==ne,"regions.strength","strength needs one value per edge");
    require(std::isfinite(w_normal) && w_normal>=0 && std::isfinite(w_preference) && w_preference>=0,"regions.weights","dual weights must be finite and nonnegative");
    for(auto s:strength)require(std::isfinite(s) && s>=0 && s<=1,"regions.strength","strength must lie in [0, 1]");
    RegionDual out;out.faces.assign(faces.begin(),faces.end());
    std::sort(out.faces.begin(),out.faces.end());out.faces.erase(std::unique(out.faces.begin(),out.faces.end()),out.faces.end());
    require(out.faces.empty() || out.faces.back()<nf,"regions.faces","face outside the mesh");
    const auto n=out.faces.size();std::vector<std::int64_t> local(nf,-1);
    for(std::size_t i=0;i<n;++i)local[out.faces[i]]=static_cast<std::int64_t>(i);
    // Must-weld contraction (union-find over the region's faces).
    std::vector<std::size_t> parent(n);std::iota(parent.begin(),parent.end(),std::size_t{0});
    auto find=[&](std::size_t x) {while(parent[x]!=x) {parent[x]=parent[parent[x]];x=parent[x];}return x;};
    struct Crossing {MeshIndex edge;std::size_t a,b;};std::vector<Crossing> crossings;
    for(MeshIndex e=0;e<ne;++e) {
        const auto pair=t.edge_faces[e];
        if(pair.size()!=2 || !t.edge_consistent[e])continue;
        const auto la=local[pair[0]],lb=local[pair[1]];
        if(la<0 || lb<0 || (!barriers.empty() && barriers[e]))continue;
        if(!weld.empty() && weld[e]) {
            ++out.welded_edges;const auto ra=find(static_cast<std::size_t>(la)),rb=find(static_cast<std::size_t>(lb));
            if(ra!=rb)parent[std::max(ra,rb)]=std::min(ra,rb);
            continue;
        }
        crossings.push_back({e,static_cast<std::size_t>(la),static_cast<std::size_t>(lb)});
    }
    // Nodes numbered by their smallest face (faces are ascending).
    std::vector<std::int64_t> node_of_root(n,-1);MeshIndex nodes=0;out.node_of_face.resize(n);
    for(std::size_t i=0;i<n;++i) {
        const auto r=find(i);if(node_of_root[r]<0)node_of_root[r]=nodes++;
        out.node_of_face[i]=static_cast<MeshIndex>(node_of_root[r]);
    }
    std::vector<std::vector<MeshIndex>> members(nodes);
    for(std::size_t i=0;i<n;++i)members[out.node_of_face[i]].push_back(out.faces[i]);
    out.node_area.assign(nodes,0.0);out.node_center.assign(nodes,Vec3{0,0,0});out.node_normal.assign(nodes,Vec3{0,0,0});
    for(MeshIndex k=0;k<nodes;++k) {
        for(auto f:members[k]) {
            out.node_faces.values.push_back(f);const double a=g.face_area[f];out.node_area[k]+=a;
            for(std::size_t d=0;d<3;++d) {out.node_center[k][d]+=g.face_center[f][d]*a;out.node_normal[k][d]+=g.face_normal[f][d]*a;}
        }
        out.node_faces.offsets.push_back(static_cast<Offset>(out.node_faces.values.size()));
        if(out.node_area[k]>0)for(auto& x:out.node_center[k])x/=out.node_area[k];
        else {out.node_center[k]=g.face_center[members[k][0]];}
        // A welded group whose normals cancel (a whole ring of a tube) has no normal: the
        // growth's normal-proxy term then treats it uniformly instead of amplifying noise.
        const double length=norm(out.node_normal[k]);
        if(length>1e-9*out.node_area[k] && length>0)for(auto& x:out.node_normal[k])x/=length;
        else out.node_normal[k]=Vec3{0,0,0};
    }
    // Crossing weights between nodes: the cheapest face pair (ties: lower edge id).
    const double floor=std::max(m.mean_edge_length()*1e-6,1e-300);
    std::map<std::pair<MeshIndex,MeshIndex>,std::pair<double,MeshIndex>> best;
    for(const auto& c:crossings) {
        const auto u=out.node_of_face[c.a],v=out.node_of_face[c.b];if(u==v)continue;
        const auto fa=out.faces[c.a],fb=out.faces[c.b];
        const double length=std::max(norm(sub(g.face_center[fa],g.face_center[fb])),floor);
        const double bend=1-std::clamp(detail::numpy_einsum_dot3(g.face_normal[fa],g.face_normal[fb]),-1.0,1.0);
        double inner=1+w_normal*bend;
        const double s=strength.empty()?0.0:strength[c.edge];
        if(w_preference>0 && s>0)inner+=w_preference*s*s;
        const double w=length*inner;const auto key=std::minmax(u,v);
        const auto found=best.find(key);
        if(found==best.end() || w<found->second.first || (w==found->second.first && c.edge<found->second.second))best[key]={w,c.edge};
    }
    struct Entry {MeshIndex row,col;double w;MeshIndex edge;};std::vector<Entry> entries;
    for(const auto& [key,value]:best) {entries.push_back({key.first,key.second,value.first,value.second});entries.push_back({key.second,key.first,value.first,value.second});}
    std::sort(entries.begin(),entries.end(),[](const Entry& a,const Entry& b) {return std::pair(a.row,a.col)<std::pair(b.row,b.col);});
    out.graph.indptr.assign(nodes+1,0);
    for(const auto& x:entries) {++out.graph.indptr[x.row+1];out.graph.indices.push_back(x.col);out.graph.weights.push_back(x.w);out.graph.edges.push_back(x.edge);}
    for(std::size_t k=0;k<nodes;++k)out.graph.indptr[k+1]+=out.graph.indptr[k];
    return out;
}
namespace {
double on_grid(double x,double q) {return q>0 && std::isfinite(x)?std::nearbyint(x/q)*q:x;}
// grow_labels with every tentative distance on the length grid `q` (T39): the same
// (distance, label, node) order, with exact ties where the geometry is tied.
std::vector<std::int64_t> grow_on_grid(const DualGraph& graph,std::span<const MeshIndex> seeds,std::optional<std::span<const Vec3>> normals,
                                       std::optional<std::span<const Vec3>> proxies,double w_normal,double q) {
    std::set<MeshIndex> distinct(seeds.begin(),seeds.end());
    require(distinct.size()==seeds.size() && std::all_of(seeds.begin(),seeds.end(),[&](auto s) {return s<graph.size();}),"segment.seeds","seeds must be distinct valid face ids");
    std::vector<std::int64_t> labels(graph.size(),-1);std::vector<double> distance(graph.size(),INFINITY);
    for(std::size_t i=0;i<seeds.size();++i) {labels[seeds[i]]=static_cast<std::int64_t>(i);distance[seeds[i]]=0.0;}
    using Item=std::tuple<double,std::int64_t,MeshIndex>;std::priority_queue<Item,std::vector<Item>,std::greater<Item>> queue;
    auto expand=[&](MeshIndex face,std::int64_t label,double value) {
        for(auto j=graph.indptr[face];j<graph.indptr[face+1];++j) {
            const auto other=graph.indices[j];if(labels[other]>=0)continue;
            double penalty=1.0;
            if(proxies)penalty+=w_normal*std::pow(1-std::clamp(dot((*normals)[other],(*proxies)[static_cast<std::size_t>(label)]),-1.0,1.0),2.0);
            const double candidate=on_grid(value+graph.weights[j]*on_grid(penalty,0x1p-32),q);
            if(candidate<distance[other]) {distance[other]=candidate;queue.push({candidate,label,other});}
        }
    };
    for(std::size_t i=0;i<seeds.size();++i)expand(seeds[i],static_cast<std::int64_t>(i),0.0);
    while(!queue.empty()) {
        const auto [value,label,face]=queue.top();queue.pop();
        if(labels[face]>=0 || value>distance[face])continue;
        labels[face]=label;expand(face,label,value);
    }
    require(std::all_of(labels.begin(),labels.end(),[](auto l) {return l>=0;}),"segment.seeds","every dual-graph component needs a seed");
    return labels;
}
}
LloydResult region_lloyd(const RegionDual& d,std::vector<MeshIndex> seeds,std::int64_t steps,const Deadline& deadline,double quantum) {
    auto grow=[&](std::span<const MeshIndex> s,std::optional<std::span<const Vec3>> normals,std::optional<std::span<const Vec3>> proxies) {
        return quantum>0?grow_on_grid(d.graph,s,normals,proxies,2.0,quantum):grow_labels(d.graph,s,normals,proxies,2.0,nullptr);
    };
    LloydResult out;out.labels=grow(seeds,std::nullopt,std::nullopt);
    for(std::int64_t step=0;step<std::max<std::int64_t>(0,steps);++step) {
        if(deadline.expired())break;
        const auto k=seeds.size();std::vector<double> area(k,0.0);std::vector<Vec3> normals(k,Vec3{0,0,0}),centres(k,Vec3{0,0,0});
        for(std::size_t v=0;v<out.labels.size();++v) {
            const auto l=static_cast<std::size_t>(out.labels[v]);const double a=d.node_area[v];area[l]+=a;
            for(std::size_t x=0;x<3;++x) {normals[l][x]+=d.node_normal[v][x]*a;centres[l][x]+=d.node_center[v][x]*a;}
        }
        for(std::size_t l=0;l<k;++l) {
            const double length=std::max(norm(normals[l]),1e-300);for(auto& x:normals[l])x/=length;
            for(auto& x:centres[l])x/=std::max(area[l],1e-300);
        }
        std::vector<MeshIndex> next;
        for(std::size_t l=0;l<k;++l) {
            std::vector<MeshIndex> ids;for(std::size_t v=0;v<out.labels.size();++v)if(out.labels[v]==static_cast<std::int64_t>(l))ids.push_back(static_cast<MeshIndex>(v));
            std::vector<double> d2(ids.size());double biggest=-INFINITY;
            for(std::size_t i=0;i<ids.size();++i) {
                const auto dx=sub(d.node_center[ids[i]],centres[l]);d2[i]=(dx[0]*dx[0]+dx[1]*dx[1])+dx[2]*dx[2];biggest=std::max(biggest,d2[i]);
            }
            std::size_t pick=0;double best=INFINITY;
            for(std::size_t i=0;i<ids.size();++i) {
                const double angular=1-std::clamp(dot(d.node_normal[ids[i]],normals[l]),-1.0,1.0);
                const double raw=d2[i]/std::max(biggest,1e-300)+.25*(angular*angular);
                const double score=quantum>0?on_grid(raw,0x1p-32):raw;
                if(score<best) {best=score;pick=i;}
            }
            next.push_back(ids[pick]);
        }
        seeds=std::move(next);
        auto proposed=grow(seeds,std::span<const Vec3>(d.node_normal),std::span<const Vec3>(normals));
        if(proposed==out.labels)break;
        out.labels=std::move(proposed);
    }
    out.seeds=std::move(seeds);
    return out;
}
}
