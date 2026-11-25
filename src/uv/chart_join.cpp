// SPDX-License-Identifier: GPL-3.0-or-later
#include "layout_postprocess_internal.hpp"
#include "cleanseams/cache.hpp"
#include "cleanseams/metrics.hpp"
#include <algorithm>
#include <cmath>
#include <numeric>
namespace cleanseams::repair_detail {
namespace {
std::array<MeshIndex,2> endpoints(const MeshSnapshot& mesh,MeshIndex corner) {
    std::array<MeshIndex,2> pair{corner,mesh.corner_next(corner)};
    if(mesh.input().corner_vertices[pair[0]]>mesh.input().corner_vertices[pair[1]])std::swap(pair[0],pair[1]);
    return pair;
}
}
std::optional<LayoutRecovery> join_small_charts(const PipelineInput& in,const std::shared_ptr<const CutSurface>& original,
    const std::function<bool()>& stopped,const std::function<bool(const LayoutRecovery&,const std::shared_ptr<const CutSurface>&,
        const LayoutRecovery&,const std::shared_ptr<const CutSurface>&,std::span<const MeshIndex>)>& valid,nlohmann::json& diagnostic) {
    const auto& mesh=*in.mesh;const auto& mt=mesh.tables();const auto& ct=original->tables();
    diagnostic["before"]=ct.n_charts;diagnostic["after"]=ct.n_charts;diagnostic["joined"]=0;
    diagnostic["attempts"]=0;diagnostic["reason"]="no_safe_join";
    diagnostic["rejections"]=nlohmann::json::object();
    diagnostic["proposal_records"]=nlohmann::json::array();
    diagnostic["work_seconds"]=nlohmann::json::object();
    nlohmann::json* record=nullptr;
    auto reject=[&](const char* reason){auto& counts=diagnostic["rejections"];counts[reason]=counts.value(reason,0)+1;if(record)(*record)["reason"]=reason;};
    auto measured=[&](const char* stage,auto&& work) {
        const auto started=MonotonicClock::now();auto value=work();
        auto& times=diagnostic["work_seconds"];
        times[stage]=times.value(stage,0.0)+std::chrono::duration<double>(MonotonicClock::now()-started).count();
        return value;
    };
    if(ct.n_charts<2 || stopped())return std::nullopt;
    if(!std::all_of(in.uv->begin(),in.uv->end(),[](auto p){return std::isfinite(p[0]) && std::isfinite(p[1]);}))return std::nullopt;
    LayoutRecovery result{*in.uv,in.seams};auto labels=ct.chart_of_face;
    const auto metric_geometry=quality_geometry(original);
    const double old_ratio=measured("quality",[&]{return face_weighted_p95(QualitySnapshot(metric_geometry,CornerUvView{result.uv}));});
    diagnostic["anisotropy_p95_before"]=old_ratio;diagnostic["anisotropy_p95_after"]=old_ratio;
    std::vector<std::vector<MeshIndex>> members(ct.n_charts);
    for(MeshIndex f=0;f<labels.size();++f)members[labels[f]].push_back(f);
    std::vector<MeshIndex> ordered(ct.n_charts);std::iota(ordered.begin(),ordered.end(),0);
    std::stable_sort(ordered.begin(),ordered.end(),[&](auto a,auto b){return members[a].size()<members[b].size();});
    struct Neighbour {MeshIndex other,edge;std::array<MeshIndex,2> source,target;};
    std::vector<std::vector<Neighbour>> neighbours(labels.size());
    for(MeshIndex e=0;e<mt.edge_corners.size();++e) {
        const auto edge=mt.edge_corners[e];if(edge.size()!=2)continue;
        if(in.locked && (*in.locked)[e] && in.seams[e])continue;
        const auto a=edge[0],b=edge[1];const auto fa=mt.corner_face[a],fb=mt.corner_face[b];
        neighbours[fa].push_back({fb,e,endpoints(mesh,a),endpoints(mesh,b)});
        neighbours[fb].push_back({fa,e,endpoints(mesh,b),endpoints(mesh,a)});
    }
    auto current=original;
    for(const auto label:ordered) {
        if(stopped())break;
        const auto group=members[label];if(group.empty() || (in.options.repair.join_mode=="small" && group.size()>4))continue;
        std::vector<MeshIndex> corners;for(auto f:group)for(auto c:mesh.face_corners(f))corners.push_back(c);
        struct Choice {double length;MeshIndex target;Neighbour edge;};std::vector<Choice> choices;
        for(auto f:group)for(const auto& neighbour:neighbours[f]) {
            const auto target=labels[neighbour.other];if(target==label)continue;
            if(members[target].size()<group.size() || (members[target].size()==group.size() && target>label))continue;
            const auto a=result.uv[neighbour.target[0]],b=result.uv[neighbour.target[1]];
            choices.push_back({std::hypot(b[0]-a[0],b[1]-a[1]),target,neighbour});
        }
        std::stable_sort(choices.begin(),choices.end(),[](const auto& a,const auto& b){return a.length>b.length;});
        if(choices.size()>8)choices.resize(8);
        for(const auto& choice:choices) {
            if(stopped())break;
            diagnostic["attempts"]=diagnostic["attempts"].get<int>()+1;
            const auto source=choice.edge.source,target=choice.edge.target;
            record=nullptr;
            if(diagnostic["proposal_records"].size()<256) {
                diagnostic["proposal_records"].push_back({{"moving_first_face",group.front()},{"moving_face_count",group.size()},
                    {"source_edge",source},{"target_edge",target},{"accepted",false},{"reason","pending"}});
                record=&diagnostic["proposal_records"].back();
            }
            const auto a=result.uv[source[0]],b=result.uv[source[1]],c=result.uv[target[0]],d=result.uv[target[1]];
            const Vec2 v{b[0]-a[0],b[1]-a[1]},w{d[0]-c[0],d[1]-c[1]};
            const double denominator=v[0]*v[0]+v[1]*v[1];
            if(!(denominator>1e-20) || !(w[0]*w[0]+w[1]*w[1]>1e-20)){reject("short_edge");continue;}
            const double cosine=(v[0]*w[0]+v[1]*w[1])/denominator,sine=(v[0]*w[1]-v[1]*w[0])/denominator;
            auto trial=result;
            for(auto id:corners) {
                const double x=result.uv[id][0]-a[0],y=result.uv[id][1]-a[1];
                trial.uv[id]={x*cosine-y*sine+c[0],x*sine+y*cosine+c[1]};
            }
            // Matching source sectors share the exact target endpoint values; this
            // prevents arithmetic noise from manufacturing a new UV discontinuity.
            for(auto id:corners)for(std::size_t k=0;k<2;++k)
                if(current->tables().sector_of_corner[id]==current->tables().sector_of_corner[source[k]])trial.uv[id]=result.uv[target[k]];
            if(in.options.repair.storage_precision=="float32")for(auto id:corners)for(auto& x:trial.uv[id])
                x=std::abs(x)>std::numeric_limits<float>::max()?NAN:static_cast<double>(static_cast<float>(x));
            bool pins=true;for(const auto& [id,value]:in.pins)if(trial.uv[id]!=value)pins=false;
            if(!pins){reject("pin");continue;}
            // Assess adjacency at requested storage precision, including unchanged
            // target corners, while returning those healthy corner values exactly.
            auto stored=trial.uv;
            if(in.options.repair.storage_precision=="float32")for(auto& p:stored)for(auto& x:p)
                x=std::abs(x)>std::numeric_limits<float>::max()?NAN:static_cast<double>(static_cast<float>(x));
            if(!std::all_of(stored.begin(),stored.end(),[](auto p){return std::isfinite(p[0]) && std::isfinite(p[1]);})){reject("nonfinite_storage");continue;}
            const auto discontinuities=seams_from_uv(mesh,stored,0);
            bool cuts=true;
            for(std::size_t e=0;e<trial.seams.size();++e) {
                if(!result.seams[e] && discontinuities[e])cuts=false;
                if(in.locked && (*in.locked)[e]) {
                    // An explicitly open protected cut cannot be closed indirectly
                    // by matching a different edge between the same two charts.
                    if(in.seams[e] && mt.edge_corners[e].size()==2 && !discontinuities[e]) {
                        const auto pair=mt.edge_corners[e];const auto left=endpoints(mesh,pair[0]),right=endpoints(mesh,pair[1]);
                        if(result.uv[left[0]]!=result.uv[right[0]] || result.uv[left[1]]!=result.uv[right[1]])cuts=false;
                    }
                } else if(!discontinuities[e] && mt.edge_faces[e].size()==2) {
                    const auto faces=mt.edge_faces[e];
                    if(labels[faces[0]]!=labels[faces[1]] && (labels[faces[0]]==label || labels[faces[1]]==label))trial.seams[e]=0;
                }
            }
            if(!cuts){reject("cut");continue;}
            auto joined=measured("cut_surface",[&]{return cut_surface(in.mesh,trial.seams);});
            if(joined->tables().n_charts>=current->tables().n_charts){reject("no_reduction");continue;}
            const double new_ratio=measured("quality",[&]{return face_weighted_p95(QualitySnapshot(metric_geometry,CornerUvView{stored}));});
            if(!std::isfinite(new_ratio) || new_ratio>old_ratio+1e-5){reject("stretch");continue;}
            if(!measured("validation",[&]{return valid(result,current,trial,joined,group);})){reject("validation");continue;}
            // A similarity attachment must not introduce a partial boundary
            // contact between its old charts. Different endpoints on the same
            // line can round to crossing segments during a later rigid pack.
            // Existing intra-chart slits remain the caller's unchanged input.
            const auto joined_chart=joined->tables().chart_of_face[group.front()];
            const std::array<MeshIndex,1> selected{joined_chart};
            BoundaryCrossingOptions contact_options;contact_options.nonconforming_contacts=true;contact_options.retain_pairs=true;contact_options.stopped=stopped;
            const auto contacts=measured("boundary_contacts",[&]{return boundary_crossings(*joined,CornerUvView{stored},contact_options,std::span<const MeshIndex>(selected));});
            if(!contacts.complete || stopped()){reject("deadline");break;}
            auto incumbent_point=[&](MeshIndex corner) {
                auto p=result.uv[corner];
                if(in.options.repair.storage_precision=="float32")for(auto& x:p)x=static_cast<double>(static_cast<float>(x));
                return p;
            };
            bool new_contact=false;
            for(const auto& pair:contacts.pairs) {
                const auto first_corner=pair.first_corner,second_corner=pair.second_corner;
                // Old charts may coincide before packing. Only a contact already
                // inside one old chart is an inherited slit, not a new contact.
                if(current->tables().chart_of_corner[first_corner]!=current->tables().chart_of_corner[second_corner] ||
                   !segments_have_nonconforming_contact(incumbent_point(first_corner),incumbent_point(mesh.corner_next(first_corner)),
                       incumbent_point(second_corner),incumbent_point(mesh.corner_next(second_corner)))) {new_contact=true;break;}
                if(stopped())break;
            }
            if(new_contact){reject("boundary_contact");continue;}
            if(stopped()){reject("deadline");continue;}
            if(record){(*record)["accepted"]=true;(*record)["reason"]="reduced";}
            result=std::move(trial);current=std::move(joined);
            // Refresh membership from the actual cut graph: exact endpoint matches
            // can join more than the edge selected for the proposal.
            std::vector<MeshIndex> representative(current->tables().n_charts,std::numeric_limits<MeshIndex>::max());
            for(MeshIndex f=0;f<labels.size();++f){auto& rep=representative[current->tables().chart_of_face[f]];rep=std::min(rep,labels[f]==label?choice.target:labels[f]);}
            for(auto& member:members)member.clear();
            for(MeshIndex f=0;f<labels.size();++f){labels[f]=representative[current->tables().chart_of_face[f]];members[labels[f]].push_back(f);}
            diagnostic["joined"]=ct.n_charts-current->tables().n_charts;diagnostic["after"]=current->tables().n_charts;
            diagnostic["reason"]="reduced";diagnostic["anisotropy_p95_after"]=new_ratio;break;
        }
    }
    if(diagnostic["joined"]==0)return std::nullopt;
    return result;
}
}
