// SPDX-License-Identifier: GPL-3.0-or-later
// Source core/absorb.py: tiny_charts, _valid, _within and absorb_tiny.
#include "cleanseams/layout.hpp"
#include "cleanseams/cache.hpp"
#include "cleanseams/errors.hpp"
#include "cleanseams/detail/numpy_compat.hpp"
#include <algorithm>
#include <cmath>
#include <numbers>
#include <set>
#include <tuple>
namespace cleanseams {
namespace {
void emit(LayoutTrace* trace,nlohmann::json event) {if(trace)trace->events.push_back(std::move(event));}
std::vector<MeshIndex> ones(std::span<const std::uint8_t> mask) {
    std::vector<MeshIndex> out;for(std::size_t i=0;i<mask.size();++i)if(mask[i])out.push_back(static_cast<MeshIndex>(i));return out;
}
// Source _valid: every validity count present and zero.
bool valid(const QualityReport& r,std::uint64_t crossings) {return !r.flipped && !r.nonfinite && !crossings && !r.collapsed_uv;}
// Source _within: valid and inside the budget scaled by `relax` (NaN never passes).
bool within(const QualityReport& r,std::uint64_t crossings,const LayoutOptions& o,double relax) {
    if(!valid(r,crossings))return false;
    for(const auto [value,limit]:{std::pair{r.anisotropy_p95,o.max_anisotropy_p95},std::pair{r.area_p95,o.max_area_p95},
                                 std::pair{r.anisotropy_p99,o.max_anisotropy_p99}})
        if(!(value<=1+(limit-1)*relax))return false;
    return r.anisotropy_max<=o.max_anisotropy_max;
}
std::vector<TransferIndex> exempt_list(const std::vector<std::uint8_t>& exempt) {
    std::vector<TransferIndex> out;for(std::size_t f=0;f<exempt.size();++f)if(exempt[f])out.push_back(static_cast<TransferIndex>(f));return out;
}
}
TinyCharts tiny_charts(const MeshSnapshot& m,const CutSurface& cs,const LayoutOptions& o) {
    const auto& t=m.tables();const auto& g=m.geometry();const auto& ct=cs.tables();const auto n=ct.n_charts;
    TinyCharts out;out.tiny.assign(n,0);out.whole.assign(n,0);out.areas.assign(n,0.0);out.counts.assign(n,0);
    std::vector<double> perimeters(n,0.0);
    for(std::size_t f=0;f<ct.chart_of_face.size();++f)out.areas[ct.chart_of_face[f]]+=g.face_area[f];
    for(auto c:ct.boundary_corners)perimeters[ct.chart_of_corner[c]]+=g.edge_length[t.corner_edge[c]];
    std::vector<std::uint64_t> component_size(g.n_components,0);for(auto c:g.face_component)++component_size[c];
    const double area_floor=o.min_island_fraction*detail::numpy_sum(g.face_area);
    for(MeshIndex c=0;c<n;++c) {
        const auto fs=ct.chart_faces[c];out.counts[c]=fs.size();
        out.whole[c]=!fs.empty() && fs.size()==component_size[g.face_component[fs[0]]];
        out.tiny[c]=below_island_policy(m,fs,out.areas[c],perimeters[c],area_floor,o) && !out.whole[c];
    }
    return out;
}
bool below_island_policy(const MeshSnapshot& m,std::span<const MeshIndex> faces,double area,double perimeter,double area_floor,const LayoutOptions& o) {
    const double iq=perimeter*perimeter/std::max(4*std::numbers::pi*area,1e-300);
    return small_face_count(m,faces,o.min_island_faces) || area<area_floor || (faces.size()<50 && iq>o.strip_isoperimetric);
}
DiskLayout absorb_tiny(const std::shared_ptr<const MeshSnapshot>& mesh,std::shared_ptr<const CutSurface> cs,const LayoutOptions& o,LayoutReport& report,
                       const ChartEvaluator* evaluate,const AbsorbRequest& request,LayoutTrace* trace) {
    if(!mesh || !cs)throw ContractError(ErrorCode::invalid_input,"absorb.mesh","missing mesh snapshot or cut surface");
    const auto& m=*mesh;const auto& t=m.tables();const auto& g=m.geometry();
    const auto ne=t.edge_vertices.size(),nf=m.input().face_offsets.size()-1;
    std::vector<std::uint8_t> locked(ne,0),exempt(nf,0);
    if(request.locked) {
        if(request.locked->size()!=ne)throw ContractError(ErrorCode::invalid_input,"absorb.locked","edge mask does not match the mesh edges");
        locked=*request.locked;
    }
    if(request.exempt) {
        if(request.exempt->size()!=nf)throw ContractError(ErrorCode::invalid_input,"absorb.exempt","face mask does not match the mesh faces");
        exempt=*request.exempt;
    }
    for(auto f:report.exempt_faces) {
        if(f>=nf)throw ContractError(ErrorCode::invalid_input,"absorb.exempt_faces","exempt face out of range");
        exempt[f]=1;
    }
    emit(trace,{{"stage","absorb_tiny"},{"unsupported_only",request.unsupported_only}});
    const auto& deadline=request.deadline;
    std::set<std::tuple<std::string,std::string,std::vector<MeshIndex>>> attempted;
    std::set<MeshIndex> unsupported_done;
    while(cs->tables().n_charts>1) {
        if(deadline.expired()) {report.timed_out=true;break;}
        const auto& ct=cs->tables();const auto info=tiny_charts(m,*cs,o);
        std::vector<MeshIndex> order;
        for(MeshIndex c=0;c<ct.n_charts;++c)
            if(request.unsupported_only?(ct.chart_unsupported[c] && !info.whole[c]):info.tiny[c]!=0)order.push_back(c);
        std::stable_sort(order.begin(),order.end(),[&](auto x,auto y) {
            return std::tuple(info.counts[x],info.areas[x],x)<std::tuple(info.counts[y],info.areas[y],y);
        });
        bool changed=false;
        for(auto a:order) {
            const auto a_faces=ct.chart_faces[a];
            if(request.unsupported_only && std::any_of(a_faces.begin(),a_faces.end(),[&](auto f) {return unsupported_done.contains(f);}))continue;
            struct Option {double first,second;MeshIndex b;std::vector<MeshIndex> edges;};
            std::vector<Option> options;
            for(const auto& adjacency:ct.adjacency) {
                const auto x=adjacency.charts[0],y=adjacency.charts[1];
                if(x==y || (a!=x && a!=y))continue;
                const auto b=a==x?y:x;
                for(std::size_t i=0;i<adjacency.arcs.size();++i) {
                    const auto arc=adjacency.arcs[i];std::vector<double> lengths;
                    for(auto e:arc)lengths.push_back(g.edge_length[e]);
                    const double length=detail::numpy_sum(lengths);
                    options.push_back(request.unsupported_only?Option{-static_cast<double>(info.counts[b]),-length,b,{arc.begin(),arc.end()}}:
                                                                Option{-length,-info.areas[b],b,{arc.begin(),arc.end()}});
                }
            }
            std::stable_sort(options.begin(),options.end(),[](const Option& x,const Option& y) {
                return std::tie(x.first,x.second,x.b,x.edges)<std::tie(y.first,y.second,y.b,y.edges);
            });
            for(const auto& option:options) {
                const auto b=option.b;const auto& edges=option.edges;
                auto sorted_edges=edges;std::sort(sorted_edges.begin(),sorted_edges.end());
                auto key=std::tuple(cs->chart_key(a).digest(),cs->chart_key(b).digest(),sorted_edges);
                if(attempted.contains(key) || std::any_of(edges.begin(),edges.end(),[&](auto e) {return locked[e]!=0;}))continue;
                attempted.insert(std::move(key));
                if(cs->euler_after_gluing(edges)!=1)continue;
                if(pins_conflict_after_gluing(*cs,request.pins,edges))continue;
                auto cut=ct.cut;for(auto e:edges)cut[e]=0;
                auto candidate=cut_surface(mesh,cut);
                const auto union_chart=candidate->tables().chart_of_face[a_faces[0]];
                if(!topological_disk(*candidate,union_chart))continue;
                if(deadline.expired()) {report.timed_out=true;break;}
                std::string rung="topology";auto new_exempt=exempt;
                if(evaluate) {
                    const std::array<MeshIndex,1> selected{union_chart};
                    auto results=(*evaluate)(candidate,selected,{});
                    const auto found=results.find(union_chart);
                    if(found==results.end())continue;
                    const auto& measured=found->second;
                    const CornerUvView view{measured.corner_uv,measured.corner_ids};
                    // Exempting a face may relax its deformation, but must never erase a
                    // real crossing from the island's visible outline.
                    if(boundary_crossings(*candidate,view,{},std::span<const MeshIndex>(selected)).counts[union_chart])continue;
                    QualityReport rep=measured.report;std::uint64_t crossings=measured.crossings;
                    const bool any_exempt=std::any_of(exempt.begin(),exempt.end(),[](auto x) {return x!=0;});
                    std::shared_ptr<const QualityGeometry> geometry;
                    auto measure=[&](const std::vector<std::uint8_t>& faces) {
                        if(!geometry)geometry=quality_geometry(candidate);
                        QualityPolicy policy;policy.exempt_faces=exempt_list(faces);
                        BoundaryCrossingOptions ignore;ignore.ignored_faces=policy.exempt_faces;
                        return std::pair(QualitySnapshot(geometry,view,policy,std::span<const MeshIndex>(selected)).report(),
                                         boundary_crossings(*candidate,view,ignore,std::span<const MeshIndex>(selected)).counts[union_chart]);
                    };
                    if(any_exempt)std::tie(rep,crossings)=measure(exempt);
                    if(o.strict_budget) {
                        if(!within(rep,crossings,o,1))continue;
                        rung="strict_budget";
                    } else if(within(rep,crossings,o,2))rung="relax_2";
                    else if(within(rep,crossings,o,3))rung="relax_3";
                    else if(valid(rep,crossings))rung="validity_only";
                    else {
                        for(auto f:a_faces)new_exempt[f]=1;
                        const auto [check,check_crossings]=measure(new_exempt);
                        const bool finite=std::all_of(measured.corner_uv.begin(),measured.corner_uv.end(),
                                                      [](const Vec2& p) {return std::isfinite(p[0]) && std::isfinite(p[1]);});
                        if(!finite || !valid(check,check_crossings))continue;
                        rung="exempt";
                    }
                    report.cache(candidate->chart_key(union_chart),measured.corner_ids,measured.corner_uv);
                }
                if(candidate->tables().chart_unsupported[union_chart])++report.unsupported_absorbed;
                ++report.absorbed;
                report.absorbed_relaxed+=rung=="relax_3" || rung=="validity_only";
                report.absorbed_exempt+=rung=="exempt";
                report.absorbed_faces.insert(report.absorbed_faces.end(),a_faces.begin(),a_faces.end());
                for(auto f:candidate->tables().chart_faces[union_chart])unsupported_done.insert(f);
                exempt=std::move(new_exempt);
                report.exempt_faces.clear();for(MeshIndex f=0;f<nf;++f)if(exempt[f])report.exempt_faces.push_back(f);
                report.notes.push_back({{"absorption",rung},{"faces",info.counts[a]},{"into_faces",info.counts[b]},{"shared_edges",edges.size()}});
                cs=std::move(candidate);changed=true;
                if(request.progress)
                    request.progress({{"stage","generate"},{"status","running"},{"msg","absorbing small islands"},
                                      {"metrics",{{"charts",cs->tables().n_charts},{"absorbed",report.absorbed}}}});
                break;
            }
            if(changed || deadline.expired())break;
        }
        if(!changed)break;
    }
    const auto info=tiny_charts(m,*cs,o);const auto& ct=cs->tables();
    report.tiny_remaining=nlohmann::json::array();
    for(MeshIndex ch=0;ch<ct.n_charts;++ch) {
        if(!info.tiny[ch])continue;
        bool removable=false;
        for(const auto& adjacency:ct.adjacency) {
            if(adjacency.charts[0]==adjacency.charts[1] || (adjacency.charts[0]!=ch && adjacency.charts[1]!=ch))continue;
            for(std::size_t i=0;i<adjacency.arcs.size() && !removable;++i) {
                const auto arc=adjacency.arcs[i];
                removable=std::none_of(arc.begin(),arc.end(),[&](auto e) {return locked[e]!=0;});
            }
        }
        const char* reason=deadline.expired()?"time_budget":removable?"no_valid_union":"protected_or_boundary";
        report.tiny_remaining.push_back({{"chart",ch},{"faces",info.counts[ch]},{"reason",reason}});
    }
    std::sort(report.absorbed_faces.begin(),report.absorbed_faces.end());
    report.absorbed_faces.erase(std::unique(report.absorbed_faces.begin(),report.absorbed_faces.end()),report.absorbed_faces.end());
    DiskLayout out;out.cut=ct.cut;out.cs=cs;
    emit(trace,{{"stage","absorb_tiny_done"},{"cut",ones(out.cut)},{"charts",ct.n_charts}});
    return out;
}
}
