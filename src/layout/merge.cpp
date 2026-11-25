// SPDX-License-Identifier: GPL-3.0-or-later
// Source core/quad_layout.py: merge_patches, _moves_for_pair, _is_strip,
// _note_rejection, cleanup_spurs, _final_spur_sweep and assemble_uv.
#include "cleanseams/layout.hpp"
#include "cleanseams/cache.hpp"
#include "cleanseams/errors.hpp"
#include "cleanseams/detail/blake2b.hpp"
#include "cleanseams/detail/numpy_compat.hpp"
#include <algorithm>
#include <cmath>
#include <cstring>
#include <numbers>
#include <set>
#include <tuple>
namespace cleanseams {
namespace {
constexpr double four_pi=4*std::numbers::pi;
void emit(LayoutTrace* trace,nlohmann::json event) {if(trace)trace->events.push_back(std::move(event));}
// T24 commit revalidation. The charts of one evaluator call are independent moves of one
// immutable surface (at most one move per chart), so they may be measured in parallel;
// a result is committed only for the chart of this surface it was computed for, and a
// result for anything else is refused, never applied. Missing: not evaluated.
const ChartMeasurement* committed(const ChartMeasurements& results,const CutSurface& cs,MeshIndex ch,const char* field) {
    const auto found=results.find(ch);
    if(found==results.end())return nullptr;
    if(!measurement_matches(cs,ch,found->second))
        throw ContractError(ErrorCode::stale_state,field,"evaluation result does not belong to chart "+std::to_string(ch)+" of this surface");
    return &found->second;
}
std::vector<MeshIndex> ones(std::span<const std::uint8_t> mask) {
    std::vector<MeshIndex> out;for(std::size_t i=0;i<mask.size();++i)if(mask[i])out.push_back(static_cast<MeshIndex>(i));return out;
}
double edge_sum(const MeshSnapshot& m,std::span<const MeshIndex> edges) {
    std::vector<double> lengths;lengths.reserve(edges.size());for(auto e:edges)lengths.push_back(m.geometry().edge_length[e]);
    return detail::numpy_sum(lengths);
}
bool any_of(std::span<const std::uint8_t> mask,std::span<const MeshIndex> edges) {
    return std::any_of(edges.begin(),edges.end(),[&](auto e) {return mask[e]!=0;});
}
std::vector<std::uint8_t> mask_or_zero(const std::optional<EdgeMask>& mask,std::size_t size,const char* field) {
    if(!mask)return std::vector<std::uint8_t>(size,0);
    if(mask->size()!=size)throw ContractError(ErrorCode::invalid_input,field,"edge mask does not match the mesh edges");
    return *mask;
}
// Source _note_rejection.
void note_rejection(LayoutReport& report,std::int64_t round,MeshIndex a,MeshIndex b,const std::string& kind,const ArcKey* key,
                    nlohmann::json info=nlohmann::json::object()) {
    if(key && (kind=="distortion" || kind=="topology"))report.frozen_reasons[*key]=kind;
    ++report.stop_reasons[kind];++report.last_round_rejections;
    nlohmann::json item{{"round",round},{"a",a},{"b",b},{"kind",kind}};
    for(auto& [name,value]:info.items())item[name]=value;
    report.stop_pairs.push_back(item);
    auto rows=report.stop_pairs.get<std::vector<nlohmann::json>>();
    auto edges=[](const nlohmann::json& row) {return row.contains("edges")?row.at("edges").get<std::int64_t>():std::int64_t{0};};
    std::stable_sort(rows.begin(),rows.end(),[&](const nlohmann::json& x,const nlohmann::json& y) {
        return std::make_tuple(-edges(x),x.at("a").get<MeshIndex>(),x.at("b").get<MeshIndex>(),x.at("kind").get<std::string>())<
               std::make_tuple(-edges(y),y.at("a").get<MeshIndex>(),y.at("b").get<MeshIndex>(),y.at("kind").get<std::string>());
    });
    if(rows.size()>10)rows.resize(10);
    report.stop_pairs=rows;
    if(report.rejections.size()<report.max_rejections) {
        nlohmann::json row{{"round",round},{"a",a},{"b",b},{"key",key?key->to_json():nlohmann::json(nullptr)},{"kind",kind}};
        for(auto& [name,value]:info.items())row[name]=value;
        report.rejections.push_back(std::move(row));
    }
}
struct Candidate {
    std::int64_t feature_first=0,rank=0;
    std::vector<double> prio;
    ArcKey key;
    MeshIndex a=0,b=0;
    std::vector<MeshIndex> edges;
    MoveKind kind=MoveKind::merge;
    bool tiny=false;
};
bool candidate_less(const Candidate& x,const Candidate& y) {
    if(x.feature_first!=y.feature_first)return x.feature_first<y.feature_first;
    if(x.rank!=y.rank)return x.rank<y.rank;
    for(std::size_t i=0;i<std::min(x.prio.size(),y.prio.size());++i)if(x.prio[i]!=y.prio[i])return x.prio[i]<y.prio[i];
    if(x.prio.size()!=y.prio.size())return x.prio.size()<y.prio.size();
    if(auto c=x.key<=>y.key;c!=0)return c<0;
    return std::pair(x.a,x.b)<std::pair(y.a,y.b);
}
// The capture tool's candidate_bytes: little-endian int64/float64 words.
std::string candidates_digest(const std::vector<Candidate>& candidates) {
    detail::Blake2b hash(20);std::vector<std::uint8_t> bytes;
    auto word=[&](auto value) {std::uint8_t raw[8];std::memcpy(raw,&value,8);bytes.insert(bytes.end(),raw,raw+8);};
    for(const auto& c:candidates) {
        bytes.clear();word(c.feature_first);word(c.rank);for(double x:c.prio)word(x);
        word(static_cast<std::int64_t>(c.key.ends.size()));for(auto v:c.key.ends)word(static_cast<std::int64_t>(v));
        word(c.key.hash);word(static_cast<std::int64_t>(c.a));word(static_cast<std::int64_t>(c.b));word(static_cast<std::int64_t>(c.kind));
        word(static_cast<std::int64_t>(c.tiny));word(static_cast<std::int64_t>(c.edges.size()));for(auto e:c.edges)word(static_cast<std::int64_t>(e));
        hash.update(bytes);
    }
    return hash.hex_digest();
}
// Source default_evaluator verdict and rounded fields of a distortion rejection.
nlohmann::json distortion_info(const ChartMeasurement& m) {
    const auto& r=m.report;auto value=[](double x) {return std::isfinite(x)?nlohmann::json(detail::python_round(x,3)):nlohmann::json(nullptr);};
    return {{"anisotropy_p95",value(r.anisotropy_p95)},{"area_p95",value(r.area_p95)},{"anisotropy_max",value(r.anisotropy_max)},
            {"flipped",r.flipped},{"crossings",m.crossings}};
}
}
bool pins_conflict_after_gluing(const CutSurface& cs,const PinMap& pins,std::span<const MeshIndex> edges) {
    if(pins.empty())return false;
    const auto& m=cs.mesh();const auto& t=m.tables();const auto& sector=cs.tables().sector_of_corner;
    const auto& corner_vertex=m.input().corner_vertices;
    std::map<MeshIndex,MeshIndex> parent;
    auto find=[&](MeshIndex s) {for(auto it=parent.find(s);it!=parent.end() && it->second!=s;it=parent.find(s))s=it->second;return s;};
    for(auto e:edges) {
        if(e>=t.edge_vertices.size())throw ContractError(ErrorCode::invalid_input,"pins.edges","edge out of range");
        if(t.edge_boundary[e] || !t.edge_manifold[e] || !t.edge_consistent[e])continue;
        const auto corners=t.edge_corners[e];if(corners.size()!=2)continue;
        const auto c0=corners[0],c1=corners[1],n0=t.corner_next[c0],n1=t.corner_next[c1];
        const bool opposite=corner_vertex[c0]==corner_vertex[n1];
        for(const auto [x,y]:{opposite?std::pair(c0,n1):std::pair(c0,c1),opposite?std::pair(n0,c1):std::pair(n0,n1)}) {
            const auto rx=find(sector[x]),ry=find(sector[y]);
            if(rx!=ry)parent[std::max(rx,ry)]=std::min(rx,ry);
        }
    }
    std::map<MeshIndex,std::vector<std::pair<MeshIndex,Vec2>>> groups;
    for(const auto& [c,uv]:pins) {
        if(c>=sector.size())throw ContractError(ErrorCode::invalid_input,"pins.corner","pinned corner out of range");
        groups[find(sector[c])].emplace_back(sector[c],uv);
    }
    for(const auto& [root,rows]:groups)
        for(const auto& [s,p]:rows)for(const auto& [u,q]:rows)if(s!=u && p!=q)return true;
    return false;
}
std::string_view move_kind_name(MoveKind kind) {
    switch(kind) {
        case MoveKind::merge: return "merge";
        case MoveKind::single_arc: return "single_arc";
        case MoveKind::self: return "self";
        case MoveKind::spur: return "spur";
    }
    return "merge";
}
const CachedChartMap* LayoutReport::cached(const ChartKey& key) const {
    const auto found=uv_cache.find(key.digest());
    return found!=uv_cache.end() && found->second.key==key?&found->second:nullptr;
}
void LayoutReport::cache(const ChartKey& key,std::vector<MeshIndex> corner_ids,std::vector<Vec2> corner_uv) {
    uv_cache.insert_or_assign(key.digest(),CachedChartMap{key,std::move(corner_ids),std::move(corner_uv)});
}
bool is_strip(const CutSurface& cs,MeshIndex chart,const IndexRows& arcs,std::int64_t max_width) {
    if(arcs.size()<2)return false;
    std::vector<std::int64_t> lengths;for(std::size_t i=0;i<arcs.size();++i)lengths.push_back(static_cast<std::int64_t>(arcs[i].size()));
    std::sort(lengths.begin(),lengths.end(),std::greater<>());
    const auto rest=static_cast<std::int64_t>(cs.tables().chart_boundary_sides[chart])-lengths[0]-lengths[1];
    return rest<=2*max_width && lengths[1]>max_width;
}
std::vector<MergeMove> moves_for_pair(const CutSurface& cs,MeshIndex a,MeshIndex b,const IndexRows& arcs,const LayoutOptions& o,LayoutReport& report) {
    std::vector<MergeMove> out;
    auto arc=[&](std::size_t i) {auto r=arcs[i];return std::vector<MeshIndex>(r.begin(),r.end());};
    if(a==b) {
        for(std::size_t i=0;i<arcs.size();++i) {
            const auto n=free_ends(cs,arcs[i]);
            if(n==1 && o.spur_removal)out.push_back({MoveKind::spur,arc(i)});
            else if(n==2 && o.arc_moves)out.push_back({MoveKind::self,arc(i)});
        }
        return out;
    }
    if(arcs.size()==1)return {{MoveKind::merge,arc(0)}};
    const auto& t=cs.tables();const auto small=t.chart_faces[a].size()<=t.chart_faces[b].size()?a:b;
    if(o.arc_moves && is_strip(cs,small,arcs,o.strip_max_width))++report.merges_rejected_strip;
    else out.push_back({MoveKind::merge,arcs.values});
    if(o.arc_moves)for(std::size_t i=0;i<arcs.size();++i)out.push_back({MoveKind::single_arc,arc(i)});
    return out;
}
DiskLayout merge_patches(const std::shared_ptr<const MeshSnapshot>& mesh,std::span<const std::uint8_t> cut,const LayoutOptions& o,LayoutReport& report,
                         const ChartEvaluator* evaluate,const MergeConstraints& constraints,LayoutTrace* trace) {
    if(!mesh)throw ContractError(ErrorCode::invalid_input,"merge.mesh","missing mesh snapshot");
    const auto& m=*mesh;const auto& t=m.tables();const auto& g=m.geometry();
    const auto ne=t.edge_vertices.size(),nv=m.input().positions.size(),nf=m.input().face_offsets.size()-1;
    if(cut.size()!=ne)throw ContractError(ErrorCode::invalid_input,"merge.cut","edge mask does not match the mesh edges");
    if(constraints.cost && constraints.cost->size()!=ne)throw ContractError(ErrorCode::invalid_input,"merge.cost","cost does not match the mesh edges");
    emit(trace,{{"stage","merge_patches"}});
    const auto locked=mask_or_zero(constraints.locked,ne,"merge.locked");const auto features=mask_or_zero(constraints.features,ne,"merge.features");
    const bool features_any=std::any_of(features.begin(),features.end(),[](auto x) {return x!=0;});
    auto K=angle_defect(m);for(auto& k:K)k=std::abs(k);
    ChartEvaluator fallback;
    if(!evaluate && o.evaluate_true_distortion) {fallback=default_evaluator(mesh,o);evaluate=&fallback;}
    // Source: max(1, int(fraction * nf)). Two charts never hold more than nf faces, so any cap
    // at or above nf behaves as nf; clamping keeps the float-to-integer conversion defined.
    const double scaled=o.max_faces_fraction*static_cast<double>(nf);
    const auto max_faces=scaled>=static_cast<double>(nf)?std::max<std::int64_t>(1,static_cast<std::int64_t>(nf)):
                         scaled<1.?std::int64_t{1}:static_cast<std::int64_t>(scaled);
    auto cs=cut_surface(mesh,cut);
    const double total_area=detail::numpy_sum(g.face_area);
    std::set<std::tuple<std::string,std::string,ArcKey>> attempted_tiny;
    const double max_defect=o.max_defect.value_or(std::numeric_limits<double>::infinity());
    for(std::int64_t rnd=0;rnd<o.max_rounds;++rnd) {
        if(constraints.deadline.expired()) {report.timed_out=true;break;}
        report.rounds=static_cast<std::uint64_t>(rnd+1);report.stop_reasons.clear();report.stop_pairs=nlohmann::json::array();report.last_round_rejections=0;
        const auto& ct=cs->tables();const auto n_charts=ct.n_charts;
        std::vector<std::uint8_t> interior(nv,1);
        for(std::size_t e=0;e<ne;++e)if(ct.cut[e])for(auto v:t.edge_vertices[e])interior[v]=0;
        // flatness proxy per chart: sum of |K| over its interior vertices
        std::vector<double> proxy(n_charts,0.0);
        for(MeshIndex v=0;v<nv;++v)if(interior[v] && !t.vert_faces[v].empty())proxy[ct.chart_of_face[t.vert_faces[v][0]]]+=K[v];
        std::vector<MeshIndex> chart_min_face(n_charts);std::vector<double> chart_area(n_charts,0.0),chart_perim(n_charts,0.0);
        std::vector<std::uint8_t> face_tiny(n_charts);
        for(MeshIndex c=0;c<n_charts;++c) {
            const auto fs=ct.chart_faces[c];chart_min_face[c]=*std::min_element(fs.begin(),fs.end());
            face_tiny[c]=small_face_count(m,fs,o.min_island_faces);
        }
        for(MeshIndex f=0;f<nf;++f)chart_area[ct.chart_of_face[f]]+=g.face_area[f];
        // removable perimeter per chart: boundary sides a move could glue
        for(auto c:ct.boundary_corners) {
            const auto e=t.corner_edge[c];
            if(t.edge_boundary[e] || locked[e] || (features_any && features[e]))continue;
            chart_perim[ct.chart_of_corner[c]]+=g.edge_length[e];
        }
        std::vector<Candidate> candidates;
        for(const auto& adjacency:ct.adjacency) {
            const auto a=adjacency.charts[0],b=adjacency.charts[1];const auto& arcs=adjacency.arcs;
            const bool pair=a!=b;
            const double small=std::min(chart_area[a],chart_area[b])/std::max(total_area,1e-300);
            const bool tiny=pair && (small<o.min_island_fraction || face_tiny[a] || face_tiny[b]);
            if(pair) {
                auto pkey=std::array<MeshIndex,2>{chart_min_face[a],chart_min_face[b]};if(pkey[0]>pkey[1])std::swap(pkey[0],pkey[1]);
                if(constraints.frozen_pairs.contains(pkey) && !tiny) {note_rejection(report,rnd,a,b,"protected_pair",nullptr);continue;}
                if(!tiny && static_cast<std::int64_t>(ct.chart_faces[a].size()+ct.chart_faces[b].size())>max_faces) {
                    note_rejection(report,rnd,a,b,"size",nullptr);continue;
                }
            }
            if((ct.chart_unsupported[a] || ct.chart_unsupported[b]) && report.tier!="triangle_regions") {
                note_rejection(report,rnd,a,b,"unsupported",nullptr);continue;
            }
            double shared=0.0;
            if(pair)for(std::size_t i=0;i<arcs.size();++i)shared+=edge_sum(m,arcs[i]);
            const bool fenced=tiny && chart_perim[chart_area[a]<=chart_area[b]?a:b]<=1e-12;
            for(auto& move:moves_for_pair(*cs,a,b,arcs,o,report)) {
                const auto& edges=move.edges;const auto kind=std::string(move_kind_name(move.kind));
                const auto n_edges=static_cast<std::int64_t>(edges.size());
                auto key=arc_key(*cs,edges);
                auto refused=[&](const std::string& why,nlohmann::json info) {note_rejection(report,rnd,a,b,why,&key,std::move(info));};
                if(report.frozen_keys.contains(key) && !tiny) {
                    const auto found=report.frozen_reasons.find(key);
                    refused(found!=report.frozen_reasons.end()?found->second:"previously_rejected",
                            {{"move",kind},{"edges",n_edges},{"cached_refusal",true}});
                    continue;
                }
                if(tiny && attempted_tiny.contains({cs->chart_key(a).digest(),cs->chart_key(b).digest(),key})) {
                    const auto found=report.frozen_reasons.find(key);
                    refused(found!=report.frozen_reasons.end()?found->second:"previously_rejected",
                            {{"move",kind},{"edges",n_edges},{"cached_refusal",true}});
                    continue;
                }
                if(o.max_merge_edges && n_edges>o.max_merge_edges) {refused("size",{{"move",kind},{"edges",n_edges}});continue;}
                if(any_of(locked,edges)) {refused("protected",{{"move",kind},{"edges",n_edges}});continue;}
                if(features_any && (!tiny || !fenced)) {
                    std::size_t hard=0;for(auto e:edges)hard+=features[e]!=0;
                    if(static_cast<double>(hard)/static_cast<double>(edges.size())>=o.feature_lock_fraction) {
                        ++report.merges_rejected_feature;refused("feature",{{"move",kind},{"edges",n_edges}});continue;
                    }
                }
                const auto euler=cs->euler_after_gluing(edges);
                if(euler!=1) {++report.merges_rejected_topology;refused("topology",{{"move",kind},{"euler",euler},{"edges",n_edges}});continue;}
                // PIN-03: a move never welds corners pinned to different coordinates.
                if(pins_conflict_after_gluing(*cs,constraints.pins,edges)) {refused("pin_conflict",{{"move",kind},{"edges",n_edges}});continue;}
                const std::set<MeshIndex> edge_set(edges.begin(),edges.end());
                std::set<MeshIndex> verts;for(auto e:edges)for(auto v:t.edge_vertices[e])verts.insert(v);
                double newly_interior=0.0;
                for(auto v:verts) {
                    if(t.vert_boundary[v])continue;const auto around=t.vert_edges[v];
                    if(std::all_of(around.begin(),around.end(),[&](auto e) {return !ct.cut[e] || edge_set.contains(e);}))newly_interior+=K[v];
                }
                const double merged_proxy=proxy[a]+(pair?proxy[b]:0.0)+newly_interior;
                if(merged_proxy>max_defect) {
                    ++report.merges_rejected_proxy;refused("proxy",{{"move",kind},{"proxy",detail::python_round(merged_proxy,3)}});continue;
                }
                const double length=edge_sum(m,edges);
                // compactness: the share of the smaller chart's perimeter the move removes
                const double least=std::min(chart_perim[a],chart_perim[b]);const bool numpy_scalar=!(1e-12>least);
                const double denominator=std::max(least,1e-12);
                const double fusion=length/denominator,outline=std::max(length,shared)/denominator;
                if(pair && outline<o.min_fusion && !tiny) {
                    ++report.merges_rejected_compactness;
                    refused("compactness",{{"move",kind},{"edges",n_edges},{"fusion",numpy_scalar?detail::numpy_round(outline,4):detail::python_round(outline,4)}});
                    continue;
                }
                if(pair && !tiny && o.max_isoperimetric_growth>0 &&
                   static_cast<std::int64_t>(std::min(ct.chart_faces[a].size(),ct.chart_faces[b].size()))>=o.isoperimetric_min_faces) {
                    const double p_union=chart_perim[a]+chart_perim[b]-2.0*std::max(length,shared),a_union=chart_area[a]+chart_area[b];
                    const double iq_union=p_union*p_union/(four_pi*std::max(a_union,1e-300));
                    const double iq_a=std::pow(chart_perim[a],2.0)/(four_pi*std::max(chart_area[a],1e-300));
                    const double iq_b=std::pow(chart_perim[b],2.0)/(four_pi*std::max(chart_area[b],1e-300));
                    if(iq_union>o.max_isoperimetric_growth*std::max({iq_a,iq_b,2.0})) {
                        ++report.merges_rejected_isoperimetric;refused("isoperimetric",{{"move",kind},{"iq",detail::python_round(iq_union,2)}});continue;
                    }
                }
                // mean seam cost of what the move removes: expensive (exposed) seams go first
                double expo=0.0;
                if(constraints.cost) {
                    std::vector<double> weighted;weighted.reserve(edges.size());
                    for(auto e:edges)weighted.push_back((*constraints.cost)[e]*g.edge_length[e]);
                    expo=detail::numpy_sum(weighted)/std::max(length,1e-12);
                }
                const double score=fusion*(1.0+o.exposure_weight*expo);
                Candidate c;
                if(o.priority=="small_first")c.prio={detail::python_round(std::log10(std::max(small,1e-12)),1),-score,-length};
                else if(o.priority=="compact")c.prio={-score,-length};
                else c.prio={-length/(1.0+newly_interior)};
                c.rank=static_cast<std::int64_t>(move.kind);
                c.feature_first=tiny && any_of(features,edges);
                c.key=std::move(key);c.a=a;c.b=b;c.edges=edges;c.kind=move.kind;c.tiny=tiny;
                candidates.push_back(std::move(c));
            }
        }
        if(candidates.empty())break;
        std::stable_sort(candidates.begin(),candidates.end(),candidate_less);
        std::set<MeshIndex> used;std::vector<const Candidate*> chosen;
        for(const auto& c:candidates) {
            if(used.contains(c.a) || used.contains(c.b))continue;
            used.insert(c.a);used.insert(c.b);chosen.push_back(&c);
        }
        if(trace) {
            nlohmann::json rows=nlohmann::json::array();
            for(const auto* c:chosen)rows.push_back({c->a,c->b,static_cast<int>(c->kind),static_cast<int>(c->tiny),c->key.hash});
            emit(trace,{{"stage","merge_round"},{"round",rnd},{"candidates",candidates.size()},{"digest",candidates_digest(candidates)},{"chosen",rows}});
        }
        auto new_cut=ct.cut;
        for(const auto* c:chosen) {
            if(c->tiny)attempted_tiny.insert({cs->chart_key(c->a).digest(),cs->chart_key(c->b).digest(),c->key});
            for(auto e:c->edges)new_cut[e]=0;
        }
        auto new_cs=cut_surface(mesh,new_cut);
        // charts of the new surface that came from a move: identified by a face of `a`
        std::vector<std::pair<MeshIndex,const Candidate*>> moved;std::map<MeshIndex,double> relax;
        for(const auto* c:chosen) {
            const auto ch=new_cs->tables().chart_of_face[ct.chart_faces[c->a][0]];
            moved.emplace_back(ch,c);if(c->tiny)relax[ch]=o.absorb_relax;
        }
        std::set<MeshIndex> reverted;
        for(const auto& [ch,c]:moved) {
            const bool disk=report.tier=="triangle_regions"?topological_disk(*new_cs,ch):new_cs->tables().chart_is_disk[ch]!=0;
            if(!disk) {
                // the Euler count said disk; a pinch or an unsupported vertex says otherwise
                for(auto e:c->edges)new_cut[e]=1;
                report.frozen_keys.insert(c->key);reverted.insert(ch);++report.merges_rejected_topology;
                auto info=new_cs->chart_summary(ch);info["move"]=std::string(move_kind_name(c->kind));
                note_rejection(report,rnd,c->a,c->b,"topology",&c->key,std::move(info));
            }
        }
        if(evaluate) {
            std::vector<MeshIndex> to_check;for(const auto& [ch,c]:moved)if(!reverted.contains(ch))to_check.push_back(ch);
            std::sort(to_check.begin(),to_check.end());
            auto results=(*evaluate)(new_cs,to_check,relax);
            // A budgeted evaluator can leave candidates unevaluated. They must revert.
            for(auto ch:to_check) {
                const auto* c=std::find_if(moved.begin(),moved.end(),[&](const auto& x) {return x.first==ch;})->second;
                const auto* measured=committed(results,*new_cs,ch,"merge.evaluate");const bool ok=measured && measured->within_budget;
                if(ok)report.cache(new_cs->chart_key(ch),measured->corner_ids,measured->corner_uv);
                else {
                    for(auto e:c->edges)new_cut[e]=1;
                    report.frozen_keys.insert(c->key);reverted.insert(ch);++report.merges_rejected_distortion;
                    if(c->kind==MoveKind::spur)++report.spurs_kept;
                    auto info=measured?distortion_info(*measured):nlohmann::json::object();
                    info["move"]=std::string(move_kind_name(c->kind));
                    note_rejection(report,rnd,c->a,c->b,"distortion",&c->key,std::move(info));
                }
            }
        }
        for(const auto& [ch,c]:moved) {
            if(reverted.contains(ch))continue;
            if(c->kind==MoveKind::merge || c->kind==MoveKind::single_arc) {
                ++report.merges_accepted;
                if(relax.contains(ch)) {
                    ++report.absorbed;++report.absorbed_relaxed;
                    const auto small_chart=chart_area[c->b]<chart_area[c->a]?c->b:c->a;
                    const auto fs=ct.chart_faces[small_chart];report.absorbed_faces.insert(report.absorbed_faces.end(),fs.begin(),fs.end());
                }
            }
            if(c->kind==MoveKind::single_arc)++report.single_arc_merges;
            else if(c->kind==MoveKind::self)++report.self_cuts_removed;
            else if(c->kind==MoveKind::spur)++report.spurs_removed;
        }
        if(!reverted.empty())new_cs=cut_surface(mesh,new_cut);
        cs=std::move(new_cs);
        if(constraints.progress)
            constraints.progress({{"stage","generate"},{"status","running"},{"done",rnd+1},{"total",o.max_rounds},{"metrics",{{"charts",cs->tables().n_charts}}}});
    }
    DiskLayout out;out.cut=cs->tables().cut;out.cs=cs;
    emit(trace,{{"stage","merge_patches_done"},{"cut",ones(out.cut)},{"charts",cs->tables().n_charts}});
    return out;
}
SpurCleanup cleanup_spurs(const MeshSnapshot& m,std::span<const std::uint8_t> seams_in,std::optional<std::span<const std::uint8_t>> locked) {
    const auto& t=m.tables();const auto ne=t.edge_vertices.size(),nv=m.input().positions.size();
    if(seams_in.size()!=ne || (locked && locked->size()!=ne))throw ContractError(ErrorCode::invalid_input,"spurs.seams","edge mask does not match the mesh edges");
    SpurCleanup out;out.seams.assign(seams_in.begin(),seams_in.end());
    std::vector<std::uint8_t> cut(ne),removable(ne);std::vector<std::int64_t> degree(nv,0);
    for(std::size_t e=0;e<ne;++e) {
        cut[e]=out.seams[e] || t.edge_boundary[e] || !t.edge_manifold[e] || !t.edge_consistent[e];
        removable[e]=out.seams[e] && !(locked && (*locked)[e]) && !t.edge_boundary[e] && t.edge_manifold[e] && t.edge_consistent[e];
        if(cut[e])for(auto v:t.edge_vertices[e])++degree[v];
    }
    std::vector<MeshIndex> stack;for(MeshIndex v=0;v<nv;++v)if(degree[v]==1)stack.push_back(v);
    while(!stack.empty()) {
        const auto v=stack.back();stack.pop_back();
        if(degree[v]!=1)continue;
        MeshIndex e=0;for(auto x:t.vert_edges[v])if(cut[x]) {e=x;break;}
        if(!removable[e])continue;
        cut[e]=0;out.seams[e]=0;out.removed.push_back(e);
        const auto w=m.edge_other(e,v);--degree[v];--degree[w];
        if(degree[w]==1)stack.push_back(w);
    }
    std::sort(out.removed.begin(),out.removed.end());
    return out;
}
std::vector<DeadEndArc> dead_end_arcs(const MeshSnapshot& m,std::span<const std::uint8_t> seams) {
    // T33: the dead ends dangling_seam_ends counts, grouped into arcs (read-only).
    const auto& t=m.tables();const auto ne=t.edge_vertices.size(),nv=m.input().positions.size();
    if(seams.size()!=ne)throw ContractError(ErrorCode::invalid_input,"spurs.seams","edge mask does not match the mesh edges");
    std::vector<std::uint32_t> degree(nv,0);
    for(std::size_t e=0;e<ne;++e)if(seams[e])for(auto v:t.edge_vertices[e])++degree[v];
    auto dangling=[&](MeshIndex v) {return degree[v]==1 && !t.vert_boundary[v];};
    std::vector<std::uint8_t> used(ne,0);std::vector<DeadEndArc> out;
    for(MeshIndex v=0;v<nv;++v) {
        if(!dangling(v))continue;
        DeadEndArc arc;arc.tips.push_back(v);arc.vertices.push_back(v);
        MeshIndex at=v;
        while(true) {
            std::optional<MeshIndex> next;
            for(auto e:t.vert_edges[at])if(seams[e] && !used[e]) {next=e;break;}
            if(!next)break;
            used[*next]=1;arc.edges.push_back(*next);at=m.edge_other(*next,at);arc.vertices.push_back(at);
            if(dangling(at)) {arc.floating=true;arc.tips.push_back(at);break;}
            if(degree[at]!=2 || t.vert_boundary[at])break;
        }
        // The other end of a floating slit already reached this tip.
        if(arc.edges.empty())continue;
        arc.junction=at;out.push_back(std::move(arc));
    }
    return out;
}
DiskLayout final_spur_sweep(const std::shared_ptr<const MeshSnapshot>& mesh,std::span<const std::uint8_t> cut,const std::shared_ptr<const CutSurface>& cs,
                            const LayoutOptions&,LayoutReport& report,const ChartEvaluator* evaluate,std::optional<std::span<const std::uint8_t>> locked,
                            LayoutTrace* trace,const PinMap& pins) {
    if(!mesh || !cs)throw ContractError(ErrorCode::invalid_input,"spurs.mesh","missing mesh snapshot or cut surface");
    const auto& m=*mesh;const auto& t=m.tables();const auto ne=t.edge_vertices.size();
    if(cut.size()!=ne)throw ContractError(ErrorCode::invalid_input,"spurs.cut","edge mask does not match the mesh edges");
    emit(trace,{{"stage","final_spur_sweep"}});
    auto finish=[&](DiskLayout out) {emit(trace,{{"stage","final_spur_sweep_done"},{"cut",ones(out.cut)},{"charts",out.cs->tables().n_charts}});return out;};
    EdgeMask seams(ne);for(std::size_t e=0;e<ne;++e)seams[e]=cut[e] && !t.edge_boundary[e];
    auto removed=cleanup_spurs(m,seams,locked).removed;
    if(removed.empty())return finish({EdgeMask(cut.begin(),cut.end()),cs});
    std::set<MeshIndex> removed_set(removed.begin(),removed.end());
    const auto& ct=cs->tables();
    if(!evaluate) {
        if(pins_conflict_after_gluing(*cs,pins,removed)) {
            // PIN-03: keep whole chains whose removal would weld conflicting pins.
            std::vector<MeshIndex> keep;
            for(const auto& adjacency:ct.adjacency) {
                if(adjacency.charts[0]!=adjacency.charts[1])continue;
                for(std::size_t i=0;i<adjacency.arcs.size();++i) {
                    const auto arc=adjacency.arcs[i];
                    if(!std::all_of(arc.begin(),arc.end(),[&](auto e) {return removed_set.contains(e);}))continue;
                    auto trial=keep;trial.insert(trial.end(),arc.begin(),arc.end());
                    if(!pins_conflict_after_gluing(*cs,pins,trial))keep=std::move(trial);
                }
            }
            removed=std::move(keep);removed_set=std::set<MeshIndex>(removed.begin(),removed.end());
        }
        EdgeMask new_cut(cut.begin(),cut.end());for(auto e:removed)new_cut[e]=0;
        for(const auto& adjacency:ct.adjacency) {
            if(adjacency.charts[0]!=adjacency.charts[1])continue;
            for(std::size_t i=0;i<adjacency.arcs.size();++i)if(removed_set.contains(adjacency.arcs[i][0]))++report.spurs_removed;
        }
        auto surface=cut_surface(mesh,new_cut);
        return finish({std::move(new_cut),std::move(surface)});
    }
    std::vector<std::pair<MeshIndex,std::vector<std::vector<MeshIndex>>>> by_chart;
    for(const auto& adjacency:ct.adjacency) {
        const auto a=adjacency.charts[0];if(a!=adjacency.charts[1])continue;
        for(std::size_t i=0;i<adjacency.arcs.size();++i) {
            const auto arc=adjacency.arcs[i];
            if(!std::all_of(arc.begin(),arc.end(),[&](auto e) {return removed_set.contains(e);}))continue;
            if(report.frozen_keys.contains(arc_key(*cs,arc))) {++report.spurs_kept;continue;}
            auto slot=std::find_if(by_chart.begin(),by_chart.end(),[&](const auto& x) {return x.first==a;});
            if(!pins.empty()) {
                std::vector<MeshIndex> trial;
                if(slot!=by_chart.end())for(const auto& kept:slot->second)trial.insert(trial.end(),kept.begin(),kept.end());
                trial.insert(trial.end(),arc.begin(),arc.end());
                if(pins_conflict_after_gluing(*cs,pins,trial)) {++report.spurs_kept;continue;}
            }
            if(slot==by_chart.end()) {by_chart.push_back({a,{}});slot=by_chart.end()-1;}
            slot->second.emplace_back(arc.begin(),arc.end());
        }
    }
    if(by_chart.empty())return finish({EdgeMask(cut.begin(),cut.end()),cs});
    EdgeMask new_cut(cut.begin(),cut.end());
    for(const auto& [a,arcs]:by_chart)for(const auto& arc:arcs)for(auto e:arc)new_cut[e]=0;
    auto new_cs=cut_surface(mesh,new_cut);
    std::vector<std::pair<MeshIndex,std::size_t>> moved;std::vector<MeshIndex> charts;
    for(std::size_t i=0;i<by_chart.size();++i) {
        const auto ch=new_cs->tables().chart_of_face[ct.chart_faces[by_chart[i].first][0]];
        auto slot=std::find_if(moved.begin(),moved.end(),[&](const auto& x) {return x.first==ch;});
        if(slot==moved.end()) {moved.emplace_back(ch,i);charts.push_back(ch);} else slot->second=i;
    }
    std::sort(charts.begin(),charts.end());
    auto results=(*evaluate)(new_cs,charts,{});
    bool reverted=false;
    for(const auto& [ch,index]:moved) {
        const auto& arcs=by_chart[index].second;const auto* measured=committed(results,*new_cs,ch,"spurs.evaluate");
        if(measured && measured->within_budget) {
            report.cache(new_cs->chart_key(ch),measured->corner_ids,measured->corner_uv);report.spurs_removed+=arcs.size();
        } else {
            for(const auto& arc:arcs) {for(auto e:arc)new_cut[e]=1;report.frozen_keys.insert(arc_key(*cs,arc));}
            report.spurs_kept+=arcs.size();reverted=true;
        }
    }
    if(reverted)new_cs=cut_surface(mesh,new_cut);
    return finish({std::move(new_cut),std::move(new_cs)});
}
std::vector<Vec2> assemble_uv(const MeshSnapshot& m,const std::shared_ptr<const CutSurface>& cs,LayoutReport& report,const ChartEvaluator& evaluate,
                              LayoutTrace* trace) {
    if(!cs || !evaluate)throw ContractError(ErrorCode::invalid_input,"assemble.evaluate","missing cut surface or evaluator");
    emit(trace,{{"stage","assemble_uv"}});
    const auto nan=std::numeric_limits<double>::quiet_NaN();
    std::vector<Vec2> uv(m.tables().corner_face.size(),Vec2{nan,nan});std::vector<MeshIndex> missing;
    const auto& ct=cs->tables();
    for(MeshIndex ch=0;ch<ct.n_charts;++ch) {
        if(ct.chart_faces[ch].empty())continue;
        if(const auto* hit=report.cached(cs->chart_key(ch))) {for(std::size_t i=0;i<hit->corner_ids.size();++i)uv[hit->corner_ids[i]]=hit->corner_uv[i];}
        else missing.push_back(ch);
    }
    if(!missing.empty()) {
        auto results=evaluate(cs,missing,{});
        for(const auto& [ch,measured]:results)committed(results,*cs,ch,"assemble.evaluate");
        for(auto& [ch,measured]:results) {
            for(std::size_t i=0;i<measured.corner_ids.size();++i)uv[measured.corner_ids[i]]=measured.corner_uv[i];
            report.chart_metrics.emplace_back(ch,std::move(measured));
        }
    }
    return uv;
}
}
