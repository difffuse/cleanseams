// SPDX-License-Identifier: GPL-3.0-or-later
#include "layout_postprocess_internal.hpp"
#include "cleanseams/cache.hpp"
#include "cleanseams/metrics.hpp"
#include "cleanseams/errors.hpp"
#include <algorithm>
#include <cmath>
#include <map>
#include <set>
namespace cleanseams::repair_detail {
namespace {
using Faces=std::set<MeshIndex>;
std::vector<MeshIndex> singles(const CutSurface& surface) {
    const auto& ct=surface.tables();std::vector<MeshIndex> result;
    for(MeshIndex f=0;f<ct.chart_of_face.size();++f)if(ct.chart_faces[ct.chart_of_face[f]].size()==1)result.push_back(f);
    return result;
}
std::array<MeshIndex,2> endpoints(const MeshSnapshot& mesh,MeshIndex c) {
    std::array<MeshIndex,2> result{c,mesh.corner_next(c)};
    if(mesh.input().corner_vertices[result[0]]>mesh.input().corner_vertices[result[1]])std::swap(result[0],result[1]);
    return result;
}
bool finite(Vec2 p){return std::isfinite(p[0]) && std::isfinite(p[1]);}
double uv_area(const QualitySnapshot& q) {
    double result=0;for(const auto& t:q.triangles())result+=.5*std::abs(t.raw_area2);return result;
}
}
std::optional<LayoutRecovery> reconstruct_singletons(const PipelineInput& in,const std::shared_ptr<const CutSurface>& original,
    std::span<const std::int8_t> orientation,const std::function<bool()>& stopped,
    const std::function<bool(const LayoutRecovery&,const std::shared_ptr<const CutSurface>&,
        const LayoutRecovery&,const std::shared_ptr<const CutSurface>&,std::span<const MeshIndex>,MeshIndex)>& valid,
    nlohmann::json& diagnostic) {
    const auto& mesh=*in.mesh;const auto& mt=mesh.tables();auto current=original;
    LayoutRecovery result{*in.uv,in.seams};auto remaining=singles(*current);
    diagnostic["before"]=remaining.size();diagnostic["after"]=remaining.size();diagnostic["attached"]=0;
    diagnostic["rebuilt_patches"]=0;diagnostic["attempts"]=0;diagnostic["reason"]="no_safe_reconstruction";
    diagnostic["remaining_faces"]=remaining;diagnostic["unavoidable_faces"]=nlohmann::json::array();diagnostic["rejections"]=nlohmann::json::object();
    auto reject=[&](const char* reason){auto& counts=diagnostic["rejections"];counts[reason]=counts.value(reason,0)+1;};
    if(remaining.empty() || stopped())return std::nullopt;
    if(!std::all_of(result.uv.begin(),result.uv.end(),finite)){diagnostic["reason"]="nonfinite_input";return std::nullopt;}
    auto rounded=[&](Vec2 p) {
        if(in.options.repair.storage_precision=="float32")for(auto& x:p)
            x=std::abs(x)>std::numeric_limits<float>::max()?NAN:static_cast<double>(static_cast<float>(x));
        return p;
    };
    auto round=[&](std::vector<Vec2>& uv){for(auto& p:uv)p=rounded(p);};
    // Numerical constraints all describe the same stored layout. Keep the original
    // values separately so unchanged corners and user pins return byte-for-byte exact.
    auto working=result.uv;round(working);
    if(!std::all_of(working.begin(),working.end(),finite)){diagnostic["reason"]="nonfinite_storage";return std::nullopt;}
    // Geometry-exempt faces have no reliable 3D distortion metric. If their
    // existing stored UVs have positive simple coverage, preserve those sectors
    // as solve constraints instead of letting a neighbouring attachment flatten
    // them invisibly to the percentile metric.
    auto covered_exempt=[&]() {
        auto fixed=mesh.geometry().face_degenerate_kind;
        const auto& triangles=mesh.triangulation();
        for(std::size_t t=0;t<triangles.corners.size();++t) {
            const auto f=triangles.faces[t];if(!fixed[f])continue;
            const auto ids=triangles.corners[t];const auto a=working[ids[0]],b=working[ids[1]],c=working[ids[2]];
            const double area=(b[0]-a[0])*(c[1]-a[1])-(b[1]-a[1])*(c[0]-a[0]);
            if(!std::isfinite(area) || !(area*orientation[f]>0))fixed[f]=0;
        }
        for(MeshIndex f=0;f<fixed.size();++f)if(fixed[f]) {
            std::vector<Vec2> polygon;for(auto c:mesh.face_corners(f))polygon.push_back(working[c]);
            if(polygon_self_intersects(polygon,0))fixed[f]=0;
        }
        return fixed;
    };
    auto preserved=covered_exempt();
    const auto metric=quality_geometry(original);
    // Ratio units: source max(old stretch, .25)+1e-5 becomes max(ratio,1.25)+1e-5.
    // The singleton source compares each proposal with the last accepted incumbent.
    // This differs from the fixed phase-original similarity-joining envelope.
    double old_ratio=face_weighted_p95(QualitySnapshot(metric,CornerUvView{working}));
    diagnostic["anisotropy_p95_before"]=old_ratio;diagnostic["anisotropy_p95_after"]=old_ratio;
    struct Neighbour {MeshIndex face,edge;std::array<MeshIndex,2> source,target;};
    std::vector<std::vector<Neighbour>> neighbours(mt.corner_face.empty()?0:mesh.input().face_offsets.size()-1);
    for(MeshIndex e=0;e<mt.edge_corners.size();++e) {
        const auto corners=mt.edge_corners[e];if(corners.size()!=2 || (in.locked && (*in.locked)[e] && in.seams[e]))continue;
        const auto a=corners[0],b=corners[1];const auto fa=mt.corner_face[a],fb=mt.corner_face[b];
        neighbours[fa].push_back({fb,e,endpoints(mesh,a),endpoints(mesh,b)});
        neighbours[fb].push_back({fa,e,endpoints(mesh,b),endpoints(mesh,a)});
    }
    for(auto f:remaining)if(neighbours[f].empty())diagnostic["unavoidable_faces"].push_back(f);
    auto accept=[&](LayoutRecovery trial,MeshIndex face,const char* method) {
        if(stopped())return false;
        auto stored=trial.uv;round(stored);
        for(std::size_t c=0;c<trial.uv.size();++c)trial.uv[c]=stored[c]==working[c]?result.uv[c]:stored[c];
        for(const auto& [c,p]:in.pins) {
            if(stored[c]!=rounded(p)){reject("pin");return false;}
            trial.uv[c]=p;
        }
        if(!std::all_of(stored.begin(),stored.end(),finite)){reject("nonfinite_storage");return false;}
        const auto discontinuities=seams_from_uv(mesh,stored,0);
        const auto before_cuts=seams_from_uv(mesh,working,0);
        for(std::size_t e=0;e<trial.seams.size();++e) {
            trial.seams[e]=trial.seams[e] || discontinuities[e];
            if(in.locked && (*in.locked)[e] &&
               (trial.seams[e]!=in.seams[e] || (in.seams[e] && before_cuts[e] && !discontinuities[e]))) {
                reject("protected_cut");return false;
            }
        }
        auto joined=cut_surface(in.mesh,trial.seams);const auto next=singles(*joined);
        if(next.size()>=remaining.size() || std::binary_search(next.begin(),next.end(),face)){reject("no_reduction");return false;}
        const double ratio=face_weighted_p95(QualitySnapshot(metric,CornerUvView{stored}));
        if(!std::isfinite(ratio) || ratio>std::max(old_ratio,1.25)+1e-5){reject("stretch");return false;}
        std::vector<MeshIndex> affected;
        for(MeshIndex f=0;f<neighbours.size();++f)for(auto c:mesh.face_corners(f))if(trial.uv[c]!=result.uv[c]){affected.push_back(f);break;}
        // Seam-only changes also need coverage and topology checks on both sides.
        for(MeshIndex e=0;e<trial.seams.size();++e)if(trial.seams[e]!=result.seams[e])for(auto f:mt.edge_faces[e])affected.push_back(f);
        std::sort(affected.begin(),affected.end());affected.erase(std::unique(affected.begin(),affected.end()),affected.end());
        if(!valid(result,current,trial,joined,affected,face)){reject("validation");return false;}
        if(stopped())return false;
        result=std::move(trial);current=std::move(joined);remaining=next;old_ratio=ratio;working=std::move(stored);preserved=covered_exempt();
        diagnostic[method]=diagnostic[method].get<int>()+1;diagnostic["after"]=remaining.size();
        diagnostic["anisotropy_p95_after"]=ratio;diagnostic["reason"]="reduced";
        return true;
    };
    auto solve=[&](const ChartSolveRequest& request,const char* method)->std::optional<ChartSolveResult> {
        if(stopped())return std::nullopt;
        diagnostic["attempts"]=diagnostic["attempts"].get<int>()+1;
        try {
            std::optional<ChartSolveResult> mapped;
            if(std::string_view(method)=="harmonic")mapped=solve_harmonic(request);
            else if(std::string_view(method)=="lscm")mapped=solve_lscm(request);
            else mapped=convex_boundary_map(request);
            if(stopped())return std::nullopt;
            return mapped;
        } catch(const ContractError& error) {
            if(error.code==ErrorCode::cancelled || error.code==ErrorCode::memory_limit)throw;
            reject(error.code==ErrorCode::constraint_conflict?"constraint_conflict":"solve_failed");
        } catch(const std::runtime_error&){reject("solve_failed");}
        return std::nullopt;
    };
    auto corners_of=[&](const Faces& faces) {
        std::vector<MeshIndex> ids;for(auto f:faces)for(auto c:mesh.face_corners(f))ids.push_back(c);return ids;
    };
    auto mapped_trial=[&](const ChartSolveResult& mapped,const std::vector<MeshIndex>& ids,const EdgeMask& seams)->std::optional<LayoutRecovery> {
        if(mapped.corner_ids!=ids || mapped.corner_uv.size()!=ids.size() || !std::all_of(mapped.corner_uv.begin(),mapped.corner_uv.end(),finite)) {
            reject("invalid_mapping");return std::nullopt;
        }
        auto trial=result;trial.seams=seams;for(std::size_t i=0;i<ids.size();++i)trial.uv[ids[i]]=mapped.corner_uv[i];return trial;
    };
    auto expand=[&](Faces& region,const CutSurface* joined,MeshIndex chart) {
        auto next=region;for(auto f:region)for(const auto& n:neighbours[f])if(!joined || joined->tables().chart_of_face[n.face]==chart)next.insert(n.face);
        region=std::move(next);
    };
    while(!remaining.empty() && !stopped()) {
        bool improved=false;const auto candidates=remaining;
        for(auto face:candidates) {
            if(stopped())break;
            const auto& ct=current->tables();if(ct.chart_faces[ct.chart_of_face[face]].size()!=1)continue;
            std::vector<MeshIndex> targets;for(const auto& n:neighbours[face])if(ct.chart_of_face[n.face]!=ct.chart_of_face[face])targets.push_back(ct.chart_of_face[n.face]);
            std::sort(targets.begin(),targets.end());targets.erase(std::unique(targets.begin(),targets.end()),targets.end());
            std::stable_sort(targets.begin(),targets.end(),[&](auto a,auto b){return ct.chart_faces[a].size()>ct.chart_faces[b].size();});
            bool accepted=false;
            // Uniform harmonic attachment. Four expanding rings, at most 128 faces;
            // sectors touching the rest of the chart stay fixed at the incumbent UVs.
            auto attach=[&](EdgeMask seams,const std::map<MeshIndex,Vec2>& leaf_pins,bool leaf) {
                auto joined=cut_surface(in.mesh,seams);const auto& jt=joined->tables();const auto chart=jt.chart_of_face[face];
                Faces faces(jt.chart_faces[chart].begin(),jt.chart_faces[chart].end());if(faces.size()<2)return false;
                const auto ids=corners_of(faces);Faces region{face};
                for(int ring=0;ring<4;++ring) {
                    if(stopped())break;expand(region,joined.get(),chart);if(region.size()>128)break;
                    std::vector<MeshIndex> outside;for(auto c:ids)if(!region.contains(mt.corner_face[c]))outside.push_back(c);
                    if(outside.empty() && !leaf)break;
                    if(outside.empty())for(auto c:ids)if(mt.corner_face[c]!=face)outside.push_back(c);
                    ChartSolveRequest request;request.geometry=quality_geometry(joined);request.chart=chart;request.validate=false;request.weights=HarmonicWeights::uniform;
                    std::map<MeshIndex,Vec2> pins;bool conflict=false;
                    auto pin=[&](MeshIndex sector,Vec2 value){const auto [it,inserted]=pins.emplace(sector,value);if(!inserted && it->second!=value)conflict=true;};
                    for(auto c:outside)pin(jt.sector_of_corner[c],working[c]);
                    // Harmonic coordinates solve the interior of a fixed boundary.
                    // Letting established exterior sectors float pulls them toward
                    // their neighbours and can flatten healthy boundary triangles.
                    for(auto c:jt.boundary_corners)if(jt.chart_of_corner[c]==chart && mt.corner_face[c]!=face) {
                        pin(jt.sector_of_corner[c],working[c]);
                        const auto next=mesh.corner_next(c);pin(jt.sector_of_corner[next],working[next]);
                    }
                    for(auto c:ids)if(preserved[mt.corner_face[c]])pin(jt.sector_of_corner[c],working[c]);
                    for(const auto& [sector,value]:leaf_pins)pin(sector,value);
                    for(const auto& [c,value]:in.pins)if(jt.chart_of_corner[c]==chart)pin(jt.sector_of_corner[c],rounded(value));
                    if(conflict){reject("constraint_conflict");continue;}
                    bool boundary_fixed=true;
                    for(auto c:jt.boundary_corners)if(jt.chart_of_corner[c]==chart)
                        boundary_fixed=boundary_fixed && pins.contains(jt.sector_of_corner[c]) &&
                            pins.contains(jt.sector_of_corner[mesh.corner_next(c)]);
                    if(!boundary_fixed){reject("unconstrained_boundary");continue;}
                    for(const auto& [sector,value]:pins)request.sector_pins.push_back({sector,value});
                    auto mapped=solve(request,"harmonic");if(!mapped)continue;
                    auto trial=mapped_trial(*mapped,ids,seams);if(trial && accept(std::move(*trial),face,"attached"))return true;
                }
                return false;
            };
            for(auto target:targets) {
                if(stopped())break;auto seams=result.seams;
                for(const auto& n:neighbours[face])if(ct.chart_of_face[n.face]==target)seams[n.edge]=0;
                if(attach(std::move(seams),{},false)){accepted=true;break;}
            }
            if(accepted){improved=true;continue;}
            // Boundary leaves keep their free-corner shape. Shared sectors remain
            // governed by the healthy chart; only the leaf-only sectors are pinned.
            for(const auto& n:neighbours[face]) {
                if(stopped())break;auto seams=result.seams;seams[n.edge]=0;
                auto joined=cut_surface(in.mesh,seams);const auto& jt=joined->tables();const auto chart=jt.chart_of_face[face];
                std::set<MeshIndex> healthy;for(auto f:jt.chart_faces[chart])if(f!=face)for(auto c:mesh.face_corners(f))healthy.insert(jt.sector_of_corner[c]);
                const auto a=working[n.source[0]],b=working[n.source[1]],c=working[n.target[0]],d=working[n.target[1]];
                const Vec2 v{b[0]-a[0],b[1]-a[1]},w{d[0]-c[0],d[1]-c[1]};const double length=v[0]*v[0]+v[1]*v[1];
                if(length<1e-20)continue;
                const double co=(v[0]*w[0]+v[1]*w[1])/length,si=(v[0]*w[1]-v[1]*w[0])/length;
                std::map<MeshIndex,Vec2> pins;
                for(auto id:mesh.face_corners(face))if(!healthy.contains(jt.sector_of_corner[id])) {
                    const double x=working[id][0]-a[0],y=working[id][1]-a[1];pins[jt.sector_of_corner[id]]={x*co-y*si+c[0],x*si+y*co+c[1]};
                }
                if(attach(std::move(seams),pins,true)){accepted=true;break;}
            }
            if(accepted){improved=true;continue;}
            // Reconstruct connected groups: complete small target charts, neighbour
            // pairs, then three rings. Never strand a partner by hiding a new singleton.
            std::vector<Faces> groups;
            for(auto target:targets)if(ct.chart_faces[target].size()<64){Faces group{face};group.insert(ct.chart_faces[target].begin(),ct.chart_faces[target].end());groups.push_back(std::move(group));}
            for(const auto& n:neighbours[face])groups.push_back({face,n.face});
            Faces region{face};for(int ring=0;ring<3;++ring){expand(region,nullptr,0);if(region.size()>=2 && region.size()<=64)groups.push_back(region);}
            for(const auto& group:groups) {
                if(stopped())break;if(group.size()<2 || group.size()>64)continue;
                bool pinned=false;for(const auto& [c,p]:in.pins)if(group.contains(mt.corner_face[c]))pinned=true;
                if(pinned){reject("pinned_patch");continue;}
                if(std::any_of(group.begin(),group.end(),[&](auto f){return preserved[f]!=0;})) {
                    reject("covered_exempt_patch");continue;
                }
                auto seams=result.seams;bool conflict=false;
                for(MeshIndex e=0;e<seams.size();++e) {
                    const auto faces=mt.edge_faces[e];if(faces.size()!=2)continue;
                    const bool a=group.contains(faces[0]),b=group.contains(faces[1]);if(!a && !b)continue;
                    const auto required=static_cast<std::uint8_t>(a!=b || (in.locked && (*in.locked)[e] && in.seams[e]));
                    if(in.locked && (*in.locked)[e] && required!=in.seams[e])conflict=true;
                    seams[e]=required;
                }
                if(conflict){reject("protected_cut");continue;}
                auto joined=cut_surface(in.mesh,seams);const auto& jt=joined->tables();const auto chart=jt.chart_of_face[face];
                if(jt.chart_faces[chart].size()!=group.size() || jt.chart_n_loops[chart]!=1 || jt.chart_euler[chart]!=1 || jt.chart_pinched[chart])continue;
                const auto ids=corners_of(group);const auto geometry=quality_geometry(joined);const std::vector<MeshIndex> one{chart};
                const QualitySnapshot before(geometry,CornerUvView{working},{},std::span<const MeshIndex>(one));
                const double local_limit=std::max(face_weighted_p95(before),2.0)+1e-5;
                for(const auto* method:{"lscm","convex"}) {
                    ChartSolveRequest request;request.geometry=geometry;request.chart=chart;request.validate=false;
                    auto mapped=solve(request,method);if(!mapped)continue;
                    auto trial=mapped_trial(*mapped,ids,seams);if(!trial)continue;
                    const QualitySnapshot raw(geometry,CornerUvView{trial->uv},{},std::span<const MeshIndex>(one));
                    const double ratio=face_weighted_p95(raw);if(!std::isfinite(ratio) || ratio>local_limit){reject("patch_stretch");continue;}
                    const double scale=std::sqrt(std::max(uv_area(before),1e-12)/std::max(uv_area(raw),1e-20));
                    Vec2 low{INFINITY,INFINITY};for(auto id:ids){trial->uv[id][0]*=scale;trial->uv[id][1]*=scale*orientation[face];for(std::size_t d=0;d<2;++d)low[d]=std::min(low[d],trial->uv[id][d]);}
                    double cursor=-INFINITY;for(const auto& p:working)cursor=std::max(cursor,p[0]);cursor+=1;
                    for(auto id:ids){trial->uv[id][0]+=cursor-low[0];trial->uv[id][1]-=low[1];}
                    if(accept(std::move(*trial),face,"rebuilt_patches")){accepted=true;break;}
                }
                if(accepted)break;
            }
            improved=improved || accepted;
        }
        if(!improved)break;
    }
    diagnostic["remaining_faces"]=remaining;
    if(remaining.size()==diagnostic["before"].get<std::size_t>())return std::nullopt;
    return result;
}
}
