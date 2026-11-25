// SPDX-License-Identifier: GPL-3.0-or-later
#include "cleanseams/layout.hpp"
#include "cleanseams/cache.hpp"
#include "cleanseams/errors.hpp"
#include "cleanseams/detail/mesh_math.hpp"
#include "cleanseams/detail/numpy_compat.hpp"
#include <algorithm>
#include <cmath>
#include <deque>
#include <numbers>
#include <set>
#include <tuple>
namespace cleanseams {
namespace {
using detail::sub;using detail::dot;using detail::norm;
void require(bool yes,const char* field,const char* message) {if(!yes)throw ContractError(ErrorCode::invalid_input,field,message);}
std::size_t nv(const MeshSnapshot& m) {return m.input().positions.size();}
std::size_t ne(const MeshSnapshot& m) {return m.tables().edge_vertices.size();}
std::size_t face_size(const MeshSnapshot& m,MeshIndex f) {return static_cast<std::size_t>(m.input().face_offsets[f+1]-m.input().face_offsets[f]);}
std::vector<MeshIndex> ones(std::span<const std::uint8_t> mask) {
    std::vector<MeshIndex> out;for(std::size_t i=0;i<mask.size();++i)if(mask[i])out.push_back(static_cast<MeshIndex>(i));return out;
}
nlohmann::json optional_json(std::optional<MeshIndex> v) {return v?nlohmann::json(*v):nlohmann::json(nullptr);}
void emit(LayoutTrace* trace,nlohmann::json event) {if(trace)trace->events.push_back(std::move(event));}
std::size_t valence(const MeshSnapshot& m,MeshIndex v) {return m.tables().vert_edges[v].size();}
std::vector<MeshIndex> face_vertices(const MeshSnapshot& m,MeshIndex f) {
    std::vector<MeshIndex> out;for(auto c:m.face_corners(f))out.push_back(static_cast<MeshIndex>(m.input().corner_vertices[c]));return out;
}
}
EdgeMask irregular_sources(const MeshSnapshot& m,std::span<const std::uint8_t> cut,bool relaxed,LayoutTrace* trace) {
    const auto& t=m.tables();require(cut.size()==ne(m),"layout.cut","cut mask does not match the mesh edges");
    std::vector<std::uint8_t> src(nv(m),0);std::vector<std::int64_t> degree(nv(m),0);
    for(std::size_t e=0;e<cut.size();++e)if(cut[e] && !t.edge_boundary[e])for(auto v:t.edge_vertices[e])++degree[v];
    for(MeshIndex v=0;v<nv(m);++v) {
        const auto& fan=t.fans[v];if(!fan)continue;
        const bool nonquad=!relaxed && std::any_of(fan->faces.begin(),fan->faces.end(),[&](auto f) {return face_size(m,f)!=4;});
        if(t.vert_boundary[v]) {if(fan->edges.size()!=3 || nonquad)src[v]=1;}
        else if(fan->edges.size()!=4 || nonquad)src[v]=1;
        if(degree[v]!=0 && degree[v]!=2)src[v]=1;
    }
    emit(trace,{{"stage","irregular_sources"},{"result",ones(src)}});
    return src;
}
PairedSources pair_sources(const MeshSnapshot& m,std::span<const std::uint8_t> sources,std::span<const std::uint8_t> cut,LayoutTrace* trace) {
    const auto& t=m.tables();const auto& g=m.geometry();
    PairedSources out;out.sources.assign(sources.begin(),sources.end());out.cut.assign(cut.begin(),cut.end());
    using Candidate=std::tuple<double,MeshIndex,MeshIndex,std::vector<MeshIndex>>;std::vector<Candidate> candidates;
    for(MeshIndex vertex=0;vertex<nv(m);++vertex) {
        const auto val=valence(m,vertex);
        if(!out.sources[vertex] || t.vert_boundary[vertex] || (val!=3 && val!=5))continue;
        for(auto edge:t.vert_edges[vertex]) {
            if(out.cut[edge] || !t.edge_manifold[edge])continue;
            const auto middle=m.edge_other(edge,vertex);
            std::vector<std::pair<MeshIndex,std::vector<MeshIndex>>> paths{{middle,{edge}}};
            if(!out.sources[middle] && t.fans[middle] && !t.vert_boundary[middle])
                for(auto e:t.vert_edges[middle])if(e!=edge && !out.cut[e] && t.edge_manifold[e])paths.push_back({m.edge_other(e,middle),{edge,e}});
            for(const auto& [end,path]:paths) {
                if(end>vertex && out.sources[end] && !t.vert_boundary[end] && valence(m,end)+val==8) {
                    std::vector<double> lengths;for(auto e:path)lengths.push_back(g.edge_length[e]);
                    candidates.push_back({detail::numpy_sum(lengths),vertex,end,path});
                }
            }
        }
    }
    std::sort(candidates.begin(),candidates.end());
    for(const auto& [length,a,b,path]:candidates) {
        if(out.sources[a] && out.sources[b]) {for(auto e:path)out.cut[e]=1;out.sources[a]=out.sources[b]=0;++out.count;}
    }
    emit(trace,{{"stage","pair_sources"},{"sources",ones(out.sources)},{"cut",ones(out.cut)},{"count",out.count}});
    return out;
}
TransferIndex motorcycle_continuation(const MeshSnapshot& m,MeshIndex vertex,MeshIndex incoming,bool relaxed,double max_angle) {
    const auto& t=m.tables();const auto& g=m.geometry();
    if(t.vert_regular[vertex])return t.fans[vertex]->straight(incoming);
    if(!relaxed || !t.fans[vertex] || t.vert_boundary[vertex] || valence(m,vertex)>=6)return -1;
    const auto& normal=g.vert_normal[vertex];const auto& p=m.input().positions[vertex];
    auto travel=sub(p,m.input().positions[m.edge_other(incoming,vertex)]);travel=sub(travel,detail::mul(normal,dot(travel,normal)));
    const double length=norm(travel);if(length<1e-15)return -1;
    travel=detail::divide(travel,length);
    std::optional<std::pair<double,std::int64_t>> best;
    for(auto edge:t.vert_edges[vertex]) {
        if(edge==incoming || !t.edge_manifold[edge] || !t.edge_consistent[edge])continue;
        auto direction=sub(m.input().positions[m.edge_other(edge,vertex)],p);direction=sub(direction,detail::mul(normal,dot(direction,normal)));
        const double n=norm(direction);
        if(n>1e-15) {const std::pair<double,std::int64_t> choice{dot(direction,travel)/n,-static_cast<std::int64_t>(edge)};if(!best || choice>*best)best=choice;}
    }
    if(!best)return -1;
    return best->first>std::cos(max_angle*(std::numbers::pi/180.0))?-best->second:-1;
}
MotorcycleTrace trace_motorcycles(const MeshSnapshot& m,std::span<const std::uint8_t> sources,std::span<const std::uint8_t> cut_in,bool relaxed,LayoutTrace* trace) {
    const auto& t=m.tables();require(sources.size()==nv(m) && cut_in.size()==ne(m),"layout.trace","sources or cut do not match the mesh");
    MotorcycleTrace out;out.cut.assign(cut_in.begin(),cut_in.end());auto& cut=out.cut;
    std::vector<std::uint8_t> visited(nv(m),0);
    for(MeshIndex v=0;v<nv(m);++v)visited[v]=sources[v] || t.vert_boundary[v] || !t.fans[v];
    for(std::size_t e=0;e<cut.size();++e)if(cut[e])for(auto v:t.edge_vertices[e])visited[v]=1;
    std::vector<std::pair<MeshIndex,MeshIndex>> active;
    for(MeshIndex v=0;v<nv(m);++v) {
        if(!sources[v] || !t.fans[v])continue;
        std::vector<MeshIndex> edges(t.vert_edges[v].begin(),t.vert_edges[v].end());std::sort(edges.begin(),edges.end());
        for(auto e:edges)if(!t.edge_boundary[e] && !cut[e])active.push_back({v,e});
    }
    out.motorcycles=active.size();
    std::vector<std::uint32_t> on_edge(ne(m),0),arrived(nv(m),0);
    while(!active.empty()) {
        for(const auto& [v,e]:active)++on_edge[e];
        struct Arrival {MeshIndex w,e;bool head_on;};std::vector<Arrival> arrivals;
        for(const auto& [v,e]:active) {cut[e]=1;arrivals.push_back({m.edge_other(e,v),e,on_edge[e]>1});}
        for(const auto& a:arrivals)++arrived[a.w];
        std::vector<std::pair<MeshIndex,MeshIndex>> next;
        for(const auto& a:arrivals) {
            if(a.head_on || visited[a.w] || arrived[a.w]>1)continue;
            const auto e2=motorcycle_continuation(m,a.w,a.e,relaxed);
            if(e2<0 || cut[static_cast<std::size_t>(e2)])continue;
            next.push_back({a.w,static_cast<MeshIndex>(e2)});
        }
        for(const auto& a:arrivals) {visited[a.w]=1;arrived[a.w]=0;}
        for(const auto& [v,e]:active)on_edge[e]=0;
        active=std::move(next);
    }
    emit(trace,{{"stage","trace_motorcycles"},{"cut",ones(cut)},{"count",out.motorcycles}});
    return out;
}
std::vector<double> angle_defect(const MeshSnapshot& m) {
    const auto& t=m.tables();const auto& V=m.input().positions;std::vector<double> sums(nv(m),0.0);
    for(std::size_t c=0;c<t.corner_face.size();++c) {
        const auto v=static_cast<std::size_t>(m.input().corner_vertices[c]);
        const auto a=sub(V[static_cast<std::size_t>(m.input().corner_vertices[m.corner_next(c)])],V[v]);
        const auto b=sub(V[static_cast<std::size_t>(m.input().corner_vertices[m.corner_prev(c)])],V[v]);
        const double la=norm(a),lb=norm(b);if(la<1e-14 || lb<1e-14)continue;
        sums[v]+=std::acos(std::max(-1.0,std::min(1.0,dot(a,b)/(la*lb))));
    }
    std::vector<double> k(nv(m));
    for(MeshIndex v=0;v<nv(m);++v)k[v]=(t.vert_boundary[v] || !t.fans[v])?0.0:2*std::numbers::pi-sums[v];
    return k;
}
EdgeMask feature_cuts(const std::shared_ptr<const MeshSnapshot>& mesh,double angle_deg,std::int64_t min_chain,std::int64_t scale,double nms_ratio,
                      std::int64_t max_extend,std::optional<std::span<const std::uint8_t>> barriers,const Deadline& deadline,LayoutReport* report,LayoutTrace* trace) {
    FeatureOptions o;o.angle_deg=angle_deg;o.min_chain=min_chain;o.rings=scale;o.ratio=nms_ratio;o.max_extend=max_extend;
    auto field=detect_features(mesh,o,barriers,deadline);
    if(report) {report->feature_chains=field.chains;report->open_features=field.open_features;report->feature_panel_faces=field.panel_faces;}
    emit(trace,{{"stage","feature_cuts"},{"cut",ones(field.cut)}});
    return std::move(field.cut);
}
std::vector<MeshIndex> trace_from(const MeshSnapshot& m,std::span<const std::uint8_t> cut,MeshIndex v,std::uint64_t max_steps) {
    const auto& t=m.tables();std::vector<MeshIndex> out;std::set<MeshIndex> seen{v};
    for(auto e0:t.vert_edges[v]) {
        if(cut[e0] || t.edge_boundary[e0])continue;
        MeshIndex cur=v,e=e0;std::uint64_t steps=0;
        while(steps<max_steps) {
            out.push_back(e);const auto w=m.edge_other(e,cur);++steps;
            if(seen.contains(w) || !t.vert_regular[w])break;
            seen.insert(w);
            const auto e2=t.fans[w]->straight(e);
            if(e2<0 || cut[static_cast<std::size_t>(e2)])break;
            cur=w;e=static_cast<MeshIndex>(e2);
        }
    }
    return out;
}
std::optional<MeshIndex> cheapest_seed(const MeshSnapshot& m,std::span<const std::uint8_t> cut,std::vector<MeshIndex> cand,CostField cost,std::size_t max_cand,LayoutTrace* trace) {
    const auto input=cand;std::optional<MeshIndex> best;
    if(!cost || cand.empty())best=cand.empty()?std::nullopt:std::optional<MeshIndex>(cand.front());
    else {
        if(cand.size()>max_cand) {
            const double step=static_cast<double>(cand.size())/static_cast<double>(max_cand);std::vector<MeshIndex> spread;
            for(std::size_t i=0;i<max_cand;++i)spread.push_back(cand[static_cast<std::size_t>(static_cast<double>(i)*step)]);
            cand=std::move(spread);
        }
        best=cand.front();double best_c=INFINITY;const auto& L=m.geometry().edge_length;
        for(auto v:cand) {
            const auto edges=trace_from(m,cut,v);if(edges.empty())continue;
            std::vector<double> weighted,lengths;
            for(auto e:edges) {weighted.push_back((*cost)[e]*L[e]);lengths.push_back(L[e]);}
            const double c=detail::numpy_sum(weighted)/std::max(detail::numpy_sum(lengths),1e-12);
            if(c<best_c-1e-12) {best=v;best_c=c;}
        }
    }
    emit(trace,{{"stage","_cheapest"},{"cand",input},{"result",optional_json(best)}});
    return best;
}
std::map<MeshIndex,std::map<MeshIndex,std::int64_t>> distance_from_other_loops(const CutSurface& cs,MeshIndex chart) {
    const auto& m=cs.mesh();const auto& t=m.tables();const auto& ct=cs.tables();
    std::map<MeshIndex,std::vector<MeshIndex>> adj;
    for(MeshIndex e=0;e<ne(m);++e) {
        if(ct.cut[e])continue;
        const auto c0=t.edge_corners[e][0];if(ct.chart_of_corner[c0]!=chart)continue;
        const auto s0=ct.sector_of_corner[c0],s1=ct.sector_of_corner[m.corner_next(c0)];adj[s0].push_back(s1);adj[s1].push_back(s0);
    }
    std::map<MeshIndex,std::set<MeshIndex>> loop_sectors;
    for(auto li:ct.chart_loops[chart]) {const auto s=cs.loop_sectors(li);loop_sectors[li]=std::set<MeshIndex>(s.begin(),s.end());}
    std::map<MeshIndex,std::map<MeshIndex,std::int64_t>> out;
    for(const auto& [li,own]:loop_sectors) {
        std::set<MeshIndex> seeds;for(const auto& [lj,secs]:loop_sectors)if(lj!=li)seeds.insert(secs.begin(),secs.end());
        std::map<MeshIndex,std::int64_t> dist;for(auto s:seeds)dist[s]=0;
        std::vector<MeshIndex> frontier(seeds.begin(),seeds.end());
        while(!frontier.empty()) {
            std::vector<MeshIndex> next;
            for(auto s:frontier) {const auto found=adj.find(s);if(found==adj.end())continue;
                for(auto x:found->second)if(!dist.contains(x)) {dist[x]=dist[s]+1;next.push_back(x);}}
            frontier=std::move(next);
        }
        out[li]=std::move(dist);
    }
    return out;
}
std::optional<MeshIndex> seed_vertex(const CutSurface& cs,MeshIndex chart,std::span<const std::uint8_t> sources,CostField cost,LayoutTrace* trace) {
    const auto& m=cs.mesh();const auto& t=m.tables();const auto& ct=cs.tables();
    const auto loops=ct.chart_loops[chart];std::optional<MeshIndex> result;
    std::map<MeshIndex,std::pair<MeshIndex,MeshIndex>> cand; // vertex -> (loop, sector) of first sighting
    for(auto li:loops)for(auto c:ct.loops[li]) {
        const auto v=static_cast<MeshIndex>(m.input().corner_vertices[c]);
        if(!sources[v] && !cand.contains(v))cand[v]={li,ct.sector_of_corner[c]};
    }
    if(!cand.empty()) {
        std::vector<MeshIndex> good;
        for(const auto& [v,info]:cand)
            if(t.fans[v] && std::any_of(t.vert_edges[v].begin(),t.vert_edges[v].end(),[&](auto e) {return !ct.cut[e];}))good.push_back(v);
        if(!good.empty()) {
            if(loops.size()<2)result=cheapest_seed(m,ct.cut,good,cost,64,trace);
            else {
                const auto dist=distance_from_other_loops(cs,chart);std::vector<std::pair<std::int64_t,MeshIndex>> scored;
                for(auto v:good) {
                    const auto [li,sec]=cand.at(v);const auto& d=dist.at(li);const auto found=d.find(sec);
                    scored.push_back({found==d.end()?0:found->second,v});
                }
                std::int64_t best=scored.front().first;for(const auto& [d,v]:scored)best=std::max(best,d);
                std::vector<MeshIndex> wide;for(const auto& [d,v]:scored)if(static_cast<double>(d)>=0.7*static_cast<double>(best))wide.push_back(v);
                result=cheapest_seed(m,ct.cut,wide,cost,64,trace);
            }
        } else result=cand.begin()->first;
    } else {
        std::set<MeshIndex> verts;for(auto f:ct.chart_faces[chart])for(auto v:face_vertices(m,f))verts.insert(v);
        for(auto v:verts)if(t.vert_regular[v] && !sources[v]) {result=v;break;}
        if(!result)for(auto v:verts)if(!sources[v]) {result=v;break;}
    }
    emit(trace,{{"stage","_seed_vertex"},{"chart",chart},{"result",optional_json(result)}});
    return result;
}
std::optional<MeshIndex> interior_seed(const CutSurface& cs,MeshIndex chart,std::span<const std::uint8_t> sources,CostField cost,LayoutTrace* trace) {
    const auto& m=cs.mesh();const auto& t=m.tables();const auto& ct=cs.tables();
    std::set<MeshIndex> verts;for(auto f:ct.chart_faces[chart])for(auto v:face_vertices(m,f))verts.insert(v);
    std::map<MeshIndex,std::int64_t> dist;std::vector<MeshIndex> frontier;
    for(auto v:verts)if(std::any_of(t.vert_edges[v].begin(),t.vert_edges[v].end(),[&](auto e) {return ct.cut[e]!=0;})) {dist[v]=0;frontier.push_back(v);}
    while(!frontier.empty()) {
        std::vector<MeshIndex> next;
        for(auto v:frontier)for(auto e:t.vert_edges[v]) {
            if(ct.cut[e])continue;
            const auto w=m.edge_other(e,v);if(verts.contains(w) && !dist.contains(w)) {dist[w]=dist[v]+1;next.push_back(w);}
        }
        frontier=std::move(next);
    }
    std::vector<std::pair<std::int64_t,MeshIndex>> scored;
    for(auto v:verts) {if(sources[v] || !t.vert_regular[v])continue;const auto f=dist.find(v);scored.push_back({f==dist.end()?0:f->second,v});}
    std::optional<MeshIndex> result;
    if(!scored.empty()) {
        std::int64_t best=scored.front().first;for(const auto& [d,v]:scored)best=std::max(best,d);
        if(!cost) {for(const auto& [d,v]:scored)if(d==best) {result=v;break;}}
        else {
            std::vector<MeshIndex> far;for(const auto& [d,v]:scored)if(static_cast<double>(d)>=0.7*static_cast<double>(best))far.push_back(v);
            result=cheapest_seed(m,ct.cut,far,cost,64,trace);
        }
    }
    emit(trace,{{"stage","_interior_seed"},{"chart",chart},{"result",optional_json(result)}});
    return result;
}
DiskLayout make_all_disks(const std::shared_ptr<const MeshSnapshot>& mesh,std::span<const std::uint8_t> cut,std::span<const std::uint8_t> sources_in,
                          LayoutReport& report,std::int64_t max_passes,CostField cost,const Deadline& deadline,bool relaxed,LayoutTrace* trace) {
    const auto& m=*mesh;std::vector<std::uint8_t> sources(sources_in.begin(),sources_in.end());DiskLayout out;
    auto finish=[&](EdgeMask c,std::shared_ptr<const CutSurface> cs) {
        out.cut=std::move(c);out.cs=std::move(cs);
        emit(trace,{{"stage","make_all_disks"},{"cut",ones(out.cut)},{"charts",out.cs->tables().n_charts}});
        return out;
    };
    for(std::int64_t pass=0;pass<max_passes;++pass) {
        if(deadline.expired()) {report.timed_out=true;return finish(EdgeMask(cut.begin(),cut.end()),cut_surface(mesh,cut));}
        auto traced=trace_motorcycles(m,sources,cut,relaxed,trace);++report.trace_passes;
        auto cs=cut_surface(mesh,traced.cut);const auto& ct=cs->tables();
        std::vector<MeshIndex> bad;for(MeshIndex c=0;c<ct.n_charts;++c)if(!ct.chart_is_disk[c] && !ct.chart_unsupported[c])bad.push_back(c);
        if(bad.empty())return finish(std::move(traced.cut),cs);
        for(auto c:bad) {
            const auto seed=seed_vertex(*cs,c,sources,cost,trace);
            if(cost)++report.seeds_by_exposure;
            if(!seed) {report.notes.push_back("chart "+std::to_string(c)+": no seed vertex, left as is");continue;}
            sources[*seed]=1;++report.seeds_added;
        }
    }
    auto traced=trace_motorcycles(m,sources,cut,relaxed,trace);
    auto cs=cut_surface(mesh,traced.cut);
    return finish(std::move(traced.cut),cs);
}
SplitLayout split_over_budget(const std::shared_ptr<const MeshSnapshot>& mesh,std::span<const std::uint8_t> cut_in,const std::shared_ptr<const CutSurface>& cs_in,
                              std::span<const std::uint8_t> sources_in,const LayoutOptions& o,LayoutReport& report,const ChartEvaluator* evaluate,
                              CostField cost,const Deadline& deadline,bool segmented,LayoutTrace* trace) {
    SplitLayout out;out.cut.assign(cut_in.begin(),cut_in.end());out.cs=cs_in;out.sources.assign(sources_in.begin(),sources_in.end());
    auto finish=[&] {
        emit(trace,{{"stage","split_over_budget"},{"cut",ones(out.cut)},{"sources",ones(out.sources)},{"charts",out.cs->tables().n_charts}});
        return out;
    };
    if(!evaluate || !*evaluate)return finish();
    for(std::int64_t pass=0;pass<o.max_split_passes;++pass) {
        if(deadline.expired()) {report.timed_out=true;break;}
        const auto& ct=out.cs->tables();std::vector<MeshIndex> charts;
        for(MeshIndex c=0;c<ct.n_charts;++c) {
            const bool disk=segmented?topological_disk(*out.cs,c):ct.chart_is_disk[c]!=0;
            if(disk && static_cast<std::int64_t>(ct.chart_faces[c].size())>=o.split_min_faces &&
               (!segmented || !report.cached(out.cs->chart_key(c))))charts.push_back(c);
        }
        const auto results=(*evaluate)(out.cs,charts,{});
        std::vector<MeshIndex> bad;report.chart_metrics.clear();
        for(const auto& [c,measured]:results) {
            if(!measured.within_budget)bad.push_back(c);
            report.chart_metrics.push_back({c,measured});
            if(measured.within_budget)report.cache(out.cs->chart_key(c),measured.corner_ids,measured.corner_uv);
        }
        if(bad.empty())return finish();
        if(segmented) {
            const auto old_count=ct.n_charts;auto split=split_regions(mesh,*out.cs,bad,cost,deadline,trace);
            out.cut=std::move(split.cut);out.cs=std::move(split.cs);report.splits+=out.cs->tables().n_charts-old_count;
            if(out.cs->tables().n_charts==old_count)return finish();
            continue;
        }
        std::uint64_t added=0;
        for(auto c:bad) {const auto v=interior_seed(*out.cs,c,out.sources,cost,trace);if(v) {out.sources[*v]=1;++added;}}
        if(!added) {report.notes.push_back("split: no seed available for over-budget charts");return finish();}
        report.seeds_added+=added;report.splits+=added;
        auto disks=make_all_disks(mesh,out.cut,out.sources,report,12,cost,deadline,o.relaxed_motorcycles,trace);
        out.cut=std::move(disks.cut);out.cs=std::move(disks.cs);
    }
    return finish();
}
}
