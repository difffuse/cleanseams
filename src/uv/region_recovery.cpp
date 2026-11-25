// SPDX-License-Identifier: GPL-3.0-or-later
#include "layout_postprocess_internal.hpp"
#include "cleanseams/cache.hpp"
#include "cleanseams/errors.hpp"
#include "cleanseams/metrics.hpp"
#include <Eigen/SVD>
#include <algorithm>
#include <cmath>
#include <numbers>
#include <numeric>
#include <set>
namespace cleanseams::repair_detail {
// The recovery stretch percentile is a step CDF with ORIGINAL FACE area repeated for
// each of its fixed triangles. The engine's interpolated triangle-area percentile
// is a different statistic and cannot replace this migration acceptance envelope.
double face_weighted_p95(const QualitySnapshot& snapshot) {
    std::map<MeshIndex,double> areas;
    for(const auto& t:snapshot.triangles())areas[t.face]+=.5*t.area3;
    std::vector<std::pair<double,double>> rows;double total=0;
    for(const auto& t:snapshot.triangles()) {
        if(!(t.singular_values[1]>1e-12))continue;
        const auto ratio=t.singular_values[0]/t.singular_values[1];
        if(!std::isfinite(ratio))return std::numeric_limits<double>::infinity();
        rows.emplace_back(ratio,areas[t.face]);total+=areas[t.face];
    }
    if(rows.empty() || !(total>0))return std::numeric_limits<double>::infinity();
    std::sort(rows.begin(),rows.end());double accumulated=0;
    for(const auto& [ratio,area]:rows){accumulated+=area;if(accumulated>=.95*total)return ratio;}
    return rows.back().first;
}
namespace {
bool finite(Vec2 p) {return std::isfinite(p[0]) && std::isfinite(p[1]);}
double cross(Vec2 a,Vec2 b) {return a[0]*b[1]-a[1]*b[0];}
Vec2 minus(Vec2 a,Vec2 b) {return {a[0]-b[0],a[1]-b[1]};}
std::set<MeshIndex> defects(const std::shared_ptr<const QualityGeometry>& geometry,const std::vector<Vec2>& uv,
                          const std::vector<MeshIndex>& charts) {
    const auto selected=std::span<const MeshIndex>(charts);
    const QualityEvaluation q(geometry,CornerUvView{uv},{},selected);
    std::set<MeshIndex> bad;
    for(auto ch:charts)for(auto f:q.offending_faces(ch))bad.insert(f);
    // Fixed source triangles expose folds even when the polygon's signed area is positive.
    for(const auto& t:q.snapshot().triangles())if(!t.finite || t.raw_area2<=2e-14)bad.insert(t.face);
    OverlapOptions strict;strict.area_tolerance=0;strict.ignore_input_degenerate=false;
    for(const auto& r:{chart_overlaps(*geometry,CornerUvView{uv},strict,selected),
                       cross_chart_overlaps(*geometry,CornerUvView{uv},strict,selected)})
        bad.insert(r.overlapping_faces.begin(),r.overlapping_faces.end());
    return bad;
}
// The fragment source chooses a greedy face cover of overlap pairs. Isolating
// every participant needlessly fragments both sides of an otherwise repairable
// contact. Keep strict native overlap predicates, but move only the chosen cover.
std::optional<std::set<MeshIndex>> fragment_damage(const std::shared_ptr<const QualityGeometry>& geometry,
    const std::vector<Vec2>& uv,const std::vector<MeshIndex>& charts) {
    const auto selected=std::span<const MeshIndex>(charts);
    const QualityEvaluation q(geometry,CornerUvView{uv},{},selected);std::set<MeshIndex> bad;
    for(auto ch:charts)for(auto f:q.offending_faces(ch))bad.insert(f);
    for(const auto& t:q.snapshot().triangles())if(!t.finite || t.raw_area2<=2e-14)bad.insert(t.face);
    OverlapOptions strict;strict.area_tolerance=0;strict.ignore_input_degenerate=false;
    std::set<std::pair<MeshIndex,MeshIndex>> edges;
    for(const auto& report:{chart_overlaps(*geometry,CornerUvView{uv},strict,selected),
                            cross_chart_overlaps(*geometry,CornerUvView{uv},strict,selected)}) {
        // Never treat a truncated witness list as the complete overlap graph.
        // The bounded optional proposal is rejected; the incumbent stays intact.
        if(report.pairs_truncated)return std::nullopt;
        for(const auto& pair:report.pairs) {
            if(pair.first_face==pair.second_face)bad.insert(pair.first_face);
            else edges.emplace(std::min(pair.first_face,pair.second_face),std::max(pair.first_face,pair.second_face));
        }
    }
    while(!edges.empty()) {
        std::map<MeshIndex,std::size_t> degree;
        for(const auto& [a,b]:edges){++degree[a];++degree[b];}
        auto chosen=degree.begin();for(auto it=degree.begin();it!=degree.end();++it)if(it->second>chosen->second)chosen=it;
        const auto face=chosen->first;bad.insert(face); // ties: lowest original region-face index
        std::erase_if(edges,[&](const auto& edge){return edge.first==face || edge.second==face;});
    }
    return bad;
}
bool isolate(const MeshSnapshot& mesh,std::vector<Vec2>& uv,const std::set<MeshIndex>& faces,
             const QualitySnapshot& before,const std::vector<MeshIndex>& corners,
             std::span<const std::int8_t> orientation={},double minimum_relative_area=0,const std::function<bool()>& stopped={}) {
    std::map<MeshIndex,double> areas;double area3=0,area2=0;
    for(const auto& t:before.triangles()) {areas[t.face]+=.5*t.area3;area3+=.5*t.area3;area2+=.5*std::abs(t.raw_area2);}
    const double density2=std::max(area2,1e-8)/std::max(area3,1e-12);
    double cursor=1;for(auto c:corners)for(auto x:uv[c])if(std::isfinite(x))cursor=std::max(cursor,x);cursor+=1;
    for(auto f:faces) {
        if(stopped && stopped())return false;
        const auto range=mesh.face_corners(f);const std::vector<MeshIndex> ids(range.begin(),range.end());const auto n=ids.size();
        Eigen::MatrixXd points(static_cast<Eigen::Index>(n),3);
        for(std::size_t i=0;i<n;++i)for(std::size_t d=0;d<3;++d)
            points(static_cast<Eigen::Index>(i),static_cast<Eigen::Index>(d))=mesh.input().positions[static_cast<std::size_t>(mesh.input().corner_vertices[ids[i]])][d];
        points.rowwise()-=points.colwise().mean().eval();
        const Eigen::JacobiSVD<Eigen::MatrixXd> svd(points,Eigen::ComputeThinV);
        const Eigen::MatrixXd xy=points*svd.matrixV().leftCols(2);
        std::vector<Vec2> chart(n);for(std::size_t i=0;i<n;++i)chart[i]={xy(static_cast<Eigen::Index>(i),0),xy(static_cast<Eigen::Index>(i),1)};
        bool positive=true,negative=true;
        for(std::size_t i=0;i<n;++i) {
            const auto turn=cross(minus(chart[(i+1)%n],chart[i]),minus(chart[(i+2)%n],chart[(i+1)%n]));
            positive=positive && turn>1e-14;negative=negative && turn< -1e-14;
        }
        for(std::size_t i=1;i+1<n;++i) {
            const auto a=cross(minus(chart[i],chart[0]),minus(chart[i+1],chart[0]));
            positive=positive && a>1e-14;negative=negative && a< -1e-14;
        }
        if(!positive && !negative)for(std::size_t i=0;i<n;++i) {
            const double angle=2*std::numbers::pi*double(i)/double(n);chart[i]={std::cos(angle),std::sin(angle)};
        }
        double area=0;for(std::size_t i=0;i<n;++i)area+=.5*cross(chart[i],chart[(i+1)%n]);
        const int sign=orientation.empty()?1:orientation[f];
        if(area*sign<0)for(auto& p:chart)p[1]*=-1;
        const double target=std::max(density2*std::max(areas[f],1e-12),minimum_relative_area*std::max(area2,1e-8));
        const double scale=std::sqrt(target/std::max(std::abs(area),1e-12));
        Vec2 low{INFINITY,INFINITY};for(auto& p:chart)for(std::size_t d=0;d<2;++d){p[d]*=scale;low[d]=std::min(low[d],p[d]);}
        double high=cursor;
        for(std::size_t i=0;i<n;++i){uv[ids[i]]={chart[i][0]-low[0]+cursor,chart[i][1]-low[1]};high=std::max(high,uv[ids[i]][0]);}
        cursor=high+.01;
    }
    return true;
}
}
std::optional<LayoutRecovery> recover_selected_regions(const PipelineInput& in,const std::shared_ptr<const CutSurface>& original,
    std::span<const std::int8_t> orientation,const std::function<bool()>& stopped,nlohmann::json& diagnostic,
    const std::function<ChartSolveResult(const ChartSolveRequest&)>& solver,bool isolate_only) {
    diagnostic["reason"]="no_selected_faces";diagnostic["connected_groups"]=0;
    diagnostic["attempts"]=0;diagnostic["groups"]=nlohmann::json::array();
    const auto& mesh=*in.mesh;const auto& mt=mesh.tables();std::set<MeshIndex> selected;
    for(MeshIndex f=0;f<in.selected_faces.size();++f)if(in.selected_faces[f])selected.insert(f);
    if(selected.empty() || stopped())return std::nullopt;
    // A pinned selected region cannot be parked or rescaled. Keep the incumbent,
    // rather than treating user pins as disposable automatic LSCM anchors.
    for(const auto& [c,p]:in.pins)if(selected.contains(mt.corner_face[c])) {
        diagnostic["reason"]="selected_pin";return std::nullopt;
    }
    auto safe=*in.uv;for(auto& p:safe)for(auto& x:p)if(!std::isfinite(x))x=0;
    const auto geometry=quality_geometry(original);
    const QualitySnapshot before(geometry,CornerUvView{safe});
    std::vector<MeshIndex> all_corners(safe.size());std::iota(all_corners.begin(),all_corners.end(),0);
    LayoutRecovery result{*in.uv,in.seams};
    if(!isolate(mesh,result.uv,selected,before,all_corners,orientation,in.options.repair.minimum_relative_area,stopped))return std::nullopt;
    diagnostic["isolated_faces"]=selected.size();
    EdgeMask group_seams=in.seams;
    for(std::size_t e=0;e<group_seams.size();++e) {
        const auto faces=mt.edge_faces[e];if(faces.size()!=2)continue;
        const bool a=selected.contains(faces[0]),b=selected.contains(faces[1]);
        if(a || b)result.seams[e]=1;
        if(a!=b)group_seams[e]=1;
    }
    if(stopped())return std::nullopt;
    const auto groups=cut_surface(in.mesh,group_seams);const auto& ct=groups->tables();
    std::vector<MeshIndex> ordered;
    for(MeshIndex ch=0;ch<ct.n_charts;++ch) {
        const auto faces=ct.chart_faces[ch];
        if(!isolate_only && faces.size()>=2 && faces.size()<=256 && selected.contains(faces[0]))ordered.push_back(ch);
    }
    std::stable_sort(ordered.begin(),ordered.end(),[&](auto a,auto b){return ct.chart_faces[a].size()>ct.chart_faces[b].size();});
    if(ordered.size()>8)ordered.resize(8);
    double cursor=1;for(const auto& p:result.uv)for(auto x:p)if(std::isfinite(x))cursor=std::max(cursor,x);cursor+=1;
    const auto group_geometry=quality_geometry(groups);
    for(auto ch:ordered) {
        if(stopped())return std::nullopt;
        const auto faces=ct.chart_faces[ch];std::vector<MeshIndex> corners;
        for(auto f:faces)for(auto c:mesh.face_corners(f))corners.push_back(c);
        std::sort(corners.begin(),corners.end());
        nlohmann::json record={{"faces",std::vector<MeshIndex>(faces.begin(),faces.end())},{"accepted",false},{"reason","ineligible_chart"}};
        auto reject=[&](const char* reason){record["reason"]=reason;diagnostic["groups"].push_back(record);};
        if(corners.size()>1024 || ct.chart_n_loops[ch]!=1 || ct.chart_euler[ch]!=1 || ct.chart_pinched[ch]){reject("ineligible_chart");continue;}
        diagnostic["attempts"]=diagnostic["attempts"].get<int>()+1;
        ChartSolveRequest request;request.geometry=group_geometry;request.chart=ch;request.validate=false;
        ChartSolveResult mapped;
        try {mapped=solver?solver(request):solve_lscm(request);}
        catch(const ContractError& error) {
            if(error.code==ErrorCode::cancelled || error.code==ErrorCode::constraint_conflict || error.code==ErrorCode::memory_limit)throw;
            reject("solve_failed");continue;
        } catch(const std::runtime_error&) {reject("solve_failed");continue;}
        if(stopped())return std::nullopt;
        if(mapped.corner_ids!=corners || mapped.corner_uv.size()!=corners.size() ||
            !std::all_of(mapped.corner_uv.begin(),mapped.corner_uv.end(),finite)){reject("invalid_mapping");continue;}
        auto trial=safe;for(std::size_t i=0;i<corners.size();++i)trial[corners[i]]=mapped.corner_uv[i];
        const std::vector<MeshIndex> one{ch};
        if(!defects(group_geometry,trial,one).empty()){reject("invalid_mapping");continue;}
        const QualitySnapshot measured(group_geometry,CornerUvView{trial},{},std::span<const MeshIndex>(one));
        const double stretch=face_weighted_p95(measured);record["anisotropy_p95"]=stretch;
        if(!std::isfinite(stretch) || stretch>2){reject("stretch_budget");continue;}
        double mapped_area=0,target_area=0;
        for(const auto& t:measured.triangles())mapped_area+=.5*std::abs(t.raw_area2);
        // Match the sum of the selected faces' isolated coverage, using the same
        // immutable source triangles throughout. No retriangulation from UVs.
        for(auto t:group_geometry->chart_triangles()[ch]) {
            const auto ids=mesh.triangulation().corners[t];
            target_area+=.5*std::abs(cross(minus(result.uv[ids[1]],result.uv[ids[0]]),minus(result.uv[ids[2]],result.uv[ids[0]])));
        }
        const double scale=std::sqrt(target_area/std::max(mapped_area,1e-30));
        Vec2 low{INFINITY,INFINITY};const int sign=orientation[faces[0]];
        for(auto& p:mapped.corner_uv){p[0]*=scale;p[1]*=scale*sign;for(std::size_t d=0;d<2;++d)low[d]=std::min(low[d],p[d]);}
        double high=cursor;
        for(std::size_t i=0;i<corners.size();++i) {
            trial[corners[i]]={mapped.corner_uv[i][0]-low[0]+cursor,mapped.corner_uv[i][1]-low[1]};
            if(in.options.repair.storage_precision=="float32")for(auto& x:trial[corners[i]])
                x=std::abs(x)>std::numeric_limits<float>::max()?NAN:static_cast<double>(static_cast<float>(x));
            high=std::max(high,trial[corners[i]][0]);
        }
        auto normalized=trial;for(auto c:corners)normalized[c][1]*=sign;
        if(!defects(group_geometry,normalized,one).empty()){reject("invalid_after_storage");continue;}
        const double stored_stretch=face_weighted_p95(QualitySnapshot(group_geometry,CornerUvView{normalized},{},std::span<const MeshIndex>(one)));
        if(!std::isfinite(stored_stretch) || stored_stretch>2){reject("stored_stretch_budget");continue;}
        if(stopped())return std::nullopt;
        for(auto c:corners)result.uv[c]=trial[c];
        for(auto f:faces)for(auto c:mesh.face_corners(f)) {
            const auto e=mt.corner_edge[c];if(!group_seams[e] && mt.edge_faces[e].size()==2)result.seams[e]=in.seams[e];
        }
        cursor=high+1;record["accepted"]=true;record["reason"]="connected";
        diagnostic["groups"].push_back(record);diagnostic["connected_groups"]=diagnostic["connected_groups"].get<int>()+1;
        diagnostic["isolated_faces"]=diagnostic["isolated_faces"].get<std::size_t>()-faces.size();
    }
    for(auto f:selected)for(auto c:mesh.face_corners(f))if(in.options.repair.storage_precision=="float32")for(auto& x:result.uv[c])
        x=std::abs(x)>std::numeric_limits<float>::max()?NAN:static_cast<double>(static_cast<float>(x));
    for(std::size_t e=0;e<result.seams.size();++e)if(in.locked && (*in.locked)[e] && result.seams[e]!=in.seams[e]) {
        diagnostic["reason"]="protected_cut";return std::nullopt;
    }
    if(stopped())return std::nullopt;
    return result;
}
std::optional<LayoutRecovery> recover_fragments(const PipelineInput& in,const std::shared_ptr<const CutSurface>& original,
    std::span<const std::int8_t> orientation,const std::function<bool()>& stopped,const LayoutPostprocessServices& services,
    const std::function<bool(const LayoutRecovery&,const std::shared_ptr<const CutSurface>&,std::span<const MeshIndex>)>& valid,
    nlohmann::json& diagnostic) {
    diagnostic["before"]=original->tables().n_charts;diagnostic["after"]=original->tables().n_charts;
    diagnostic["recovered_regions"]=0;diagnostic["recovered_faces"]=0;diagnostic["rejected_regions"]=0;
    diagnostic["proposal_calls"]=0;diagnostic["unwrap_calls"]=0;diagnostic["groups"]=nlohmann::json::array();
    if(!services.region_proposal || !services.fragment_unwrap) {
        diagnostic["status"]="skipped";diagnostic["reason"]="host_capability_unavailable";return std::nullopt;
    }
    const auto& method=services.fragment_unwrap_method;
    if(method!="ANGLE_BASED" && method!="CONFORMAL" && method!="MINIMUM_STRETCH")
        throw ContractError(ErrorCode::invalid_option,"fragment_unwrap_method","expected a fixed-cut host unwrap method");
    diagnostic["reason"]="no_eligible_regions";
    if(stopped())return std::nullopt;
    const auto& mesh=*in.mesh;const auto& mt=mesh.tables();const auto nf=mesh.input().face_offsets.size()-1;
    // Eligibility follows actual UV islands. Explicit unprotected seam flags do not
    // turn a continuous supplied layout into accidental singleton fragmentation.
    const auto islands=cut_surface(in.mesh,seams_from_uv(mesh,*in.uv,0.0));
    std::vector<std::uint8_t> eligible(nf,0);std::vector<MeshIndex> parent(nf);
    std::iota(parent.begin(),parent.end(),0);
    for(MeshIndex f=0;f<nf;++f)eligible[f]=islands->tables().chart_faces[islands->tables().chart_of_face[f]].size()==1;
    for(MeshIndex e=0;e<in.seams.size();++e)if(in.locked && (*in.locked)[e] && in.seams[e])
        for(auto f:mt.edge_faces[e])eligible[f]=0;
    auto root=[&](MeshIndex f) {while(parent[f]!=f){parent[f]=parent[parent[f]];f=parent[f];}return f;};
    for(MeshIndex e=0;e<in.seams.size();++e) {
        if((e&1023)==0 && stopped())return std::nullopt;
        const auto faces=mt.edge_faces[e];if(faces.size()!=2 || !eligible[faces[0]] || !eligible[faces[1]])continue;
        const auto a=root(faces[0]),b=root(faces[1]);parent[std::max(a,b)]=std::min(a,b);
    }
    std::map<MeshIndex,std::vector<MeshIndex>> components;
    for(MeshIndex f=0;f<nf;++f)if(eligible[f])components[root(f)].push_back(f);
    std::vector<std::vector<MeshIndex>> groups;
    for(auto& [id,faces]:components)if(faces.size()>=4 && faces.size()<=512)groups.push_back(std::move(faces));
    std::stable_sort(groups.begin(),groups.end(),[](const auto& a,const auto& b){return a.size()>b.size();});
    if(groups.size()>8)groups.resize(8);
    LayoutRecovery result{*in.uv,in.seams};bool improved=false;
    for(const auto& faces:groups) {
        if(stopped())break;
        nlohmann::json record={{"faces",faces},{"accepted",false}};
        auto reject=[&](const char* reason) {record["reason"]=reason;diagnostic["groups"].push_back(record);
            diagnostic["rejected_regions"]=diagnostic["rejected_regions"].get<int>()+1;diagnostic["reason"]=reason;};
        bool pinned=false;for(const auto& [c,p]:in.pins)if(std::binary_search(faces.begin(),faces.end(),mt.corner_face[c]))pinned=true;
        if(pinned){reject("selected_pin");continue;}
        MeshInput sub;std::map<TransferIndex,TransferIndex> vertices;std::vector<TransferIndex> corners;
        for(auto f:faces)for(auto c:mesh.face_corners(f)){vertices.emplace(mesh.input().corner_vertices[c],0);corners.push_back(c);}
        for(auto& [v,index]:vertices){index=static_cast<TransferIndex>(sub.positions.size());sub.positions.push_back(mesh.input().positions[static_cast<std::size_t>(v)]);}
        for(auto f:faces){for(auto c:mesh.face_corners(f))sub.corner_vertices.push_back(vertices.at(mesh.input().corner_vertices[c]));sub.face_offsets.push_back(sub.corner_vertices.size());}
        std::optional<std::vector<Vec2>> proposed;
        auto call=[&](const auto& operation) {
            try {proposed=operation();}
            catch(const ContractError& error) {if(error.code!=ErrorCode::host_failure)throw;proposed.reset();}
            catch(const std::runtime_error&) {proposed.reset();}
            return proposed && proposed->size()==corners.size() && std::all_of(proposed->begin(),proposed->end(),finite);
        };
        diagnostic["proposal_calls"]=diagnostic["proposal_calls"].get<int>()+1;
        const bool projection_ok=call([&]{return services.region_proposal(sub,corners,stopped);});
        if(stopped())break;
        if(!projection_ok){reject("invalid_projection");continue;}
        const auto local_mesh=std::make_shared<const MeshSnapshot>(sub,in.memory_budget);
        auto local_cut=cut_surface(local_mesh,seams_from_uv(*local_mesh,*proposed,0.0));
        HostSolveRequest request;request.method=method;request.initial_uv=*proposed;request.scratch_to_original_corner=corners;
        request.mesh.face_offsets=sub.face_offsets;
        for(auto v:local_cut->tables().sector_vert)request.mesh.positions.push_back(sub.positions[v]);
        for(auto sector:local_cut->tables().sector_of_corner)request.mesh.corner_vertices.push_back(sector);
        if(stopped())break;
        diagnostic["unwrap_calls"]=diagnostic["unwrap_calls"].get<int>()+1;
        const bool unwrap_ok=call([&]{return services.fragment_unwrap(request,stopped);});
        if(stopped())break;
        if(!unwrap_ok){reject("invalid_unwrap");continue;}
        auto values=std::move(*proposed);
        // Normalize the host map's global orientation for defect/isolation checks;
        // restore the incumbent's reference sign only when parking the proposal.
        std::int64_t vote=0;
        for(const auto& ids:local_mesh->triangulation().corners) {
            const auto area=cross(minus(values[ids[1]],values[ids[0]]),minus(values[ids[2]],values[ids[0]]));vote+=(area>0)-(area<0);
        }
        if(vote<0)for(auto& p:values)p[1]*=-1;
        local_cut=cut_surface(local_mesh,seams_from_uv(*local_mesh,values,0.0));
        auto geometry=quality_geometry(local_cut);
        std::vector<MeshIndex> charts(local_cut->tables().n_charts);std::iota(charts.begin(),charts.end(),0);
        const auto damage=fragment_damage(geometry,values,charts);
        if(!damage){reject("overlap_pair_budget");continue;}
        const auto& damaged=*damage;
        std::vector<MeshIndex> local_corners(corners.size());std::iota(local_corners.begin(),local_corners.end(),0);
        if(!damaged.empty() && !isolate(*local_mesh,values,damaged,QualitySnapshot(geometry,CornerUvView{values}),local_corners,{},0,stopped))break;
        record["isolated_faces"]=damaged;
        local_cut=cut_surface(local_mesh,seams_from_uv(*local_mesh,values,0.0));geometry=quality_geometry(local_cut);
        charts.resize(local_cut->tables().n_charts);std::iota(charts.begin(),charts.end(),0);
        if(local_cut->tables().n_charts>=faces.size() || !defects(geometry,values,charts).empty()){reject("invalid_or_unimproved");continue;}
        const QualitySnapshot measured(geometry,CornerUvView{values});const double stretch=face_weighted_p95(measured);
        record["anisotropy_p95"]=std::isfinite(stretch)?nlohmann::json(stretch):nlohmann::json(nullptr);
        if(!std::isfinite(stretch) || stretch>1.25){reject("stretch_budget");continue;}
        double old_area=0,new_area=0;
        for(const auto& ids:local_mesh->triangulation().corners) {
            old_area+=.5*std::abs(cross(minus(result.uv[static_cast<std::size_t>(corners[ids[1]])],result.uv[static_cast<std::size_t>(corners[ids[0]])]),
                                        minus(result.uv[static_cast<std::size_t>(corners[ids[2]])],result.uv[static_cast<std::size_t>(corners[ids[0]])])));
            new_area+=.5*std::abs(cross(minus(values[ids[1]],values[ids[0]]),minus(values[ids[2]],values[ids[0]])));
        }
        const double scale=std::sqrt(std::max(old_area,1e-12)/std::max(new_area,1e-12));
        Vec2 low{INFINITY,INFINITY};double cursor=-INFINITY;for(auto p:result.uv)cursor=std::max(cursor,p[0]);cursor+=1;
        for(auto& p:values){p[0]*=scale;p[1]*=scale*orientation[faces[0]];for(std::size_t d=0;d<2;++d)low[d]=std::min(low[d],p[d]);}
        for(auto& p:values) {
            p={p[0]-low[0]+cursor,p[1]-low[1]};
            if(in.options.repair.storage_precision=="float32")for(auto& x:p)
                x=std::abs(x)>std::numeric_limits<float>::max()?NAN:static_cast<double>(static_cast<float>(x));
        }
        auto normalized=values;for(auto& p:normalized)p[1]*=orientation[faces[0]];
        local_cut=cut_surface(local_mesh,seams_from_uv(*local_mesh,values,0.0));geometry=quality_geometry(local_cut);
        charts.resize(local_cut->tables().n_charts);std::iota(charts.begin(),charts.end(),0);
        if(local_cut->tables().n_charts>=faces.size() || !defects(geometry,normalized,charts).empty()){reject("invalid_after_storage");continue;}
        const double stored_stretch=face_weighted_p95(QualitySnapshot(geometry,CornerUvView{normalized}));
        record["stored_anisotropy_p95"]=std::isfinite(stored_stretch)?nlohmann::json(stored_stretch):nlohmann::json(nullptr);
        if(!std::isfinite(stored_stretch) || stored_stretch>1.25){reject("stored_stretch_budget");continue;}
        auto trial=result;for(std::size_t c=0;c<corners.size();++c)trial.uv[static_cast<std::size_t>(corners[c])]=values[c];
        const auto discontinuities=seams_from_uv(mesh,trial.uv,0.0);
        for(auto f:faces)for(auto c:mesh.face_corners(f))trial.seams[mt.corner_edge[c]]=discontinuities[mt.corner_edge[c]];
        const auto joined=cut_surface(in.mesh,trial.seams);
        if(stopped())break;
        if(!valid(trial,joined,faces)){reject("final_validation");continue;}
        if(stopped())break;
        result=std::move(trial);improved=true;record["accepted"]=true;record["reason"]="recovered";
        diagnostic["groups"].push_back(record);diagnostic["reason"]="recovered";
        diagnostic["recovered_regions"]=diagnostic["recovered_regions"].get<int>()+1;
        diagnostic["recovered_faces"]=diagnostic["recovered_faces"].get<std::size_t>()+faces.size();
        diagnostic["after"]=joined->tables().n_charts;
    }
    return improved?std::optional(std::move(result)):std::nullopt;
}
std::optional<LocalRecovery> recover_mixed_chart(const RepairContext& ctx,const RepairState& state,MeshIndex chart,
                                                const Deadline& operation,const Deadline& local,nlohmann::json& diagnostic) {
    diagnostic={{"chart",chart},{"accepted",false},{"reason","ineligible_chart"}};
    auto stopped=[&] {throw_if_cancelled(ctx.cancel);const bool expired=operation.expired() || local.expired();if(expired)diagnostic["reason"]="deadline";return expired;};
    const auto& ct=state.cs->tables();const auto faces=ct.chart_faces[chart];
    diagnostic["topology"]=state.cs->chart_summary(chart);
    // Match the source recovery eligibility: a disk of CUT SECTORS (one loop,
    // Euler one, no pinched sector). The mesh-level unsupported-fan flag is still
    // reported, but does not prove this cut-sector parameterization is invalid.
    if(faces.size()<2 || faces.size()>256 || ct.chart_n_loops[chart]!=1 || ct.chart_euler[chart]!=1 || ct.chart_pinched[chart] || stopped())return std::nullopt;
    LocalRecovery result;
    for(auto f:faces)for(auto c:ctx.mesh->face_corners(f))result.corners.push_back(c);
    std::sort(result.corners.begin(),result.corners.end());if(result.corners.size()>1024)return std::nullopt;
    if(ctx.pins)for(const auto& [c,p]:*ctx.pins)if(ct.chart_of_corner[c]==chart)return std::nullopt;
    const auto geometry=quality_geometry(state.cs);
    ChartSolveRequest request;request.geometry=geometry;request.chart=chart;request.validate=false;request.pin_completion=ctx.pin_completion;
    diagnostic["reason"]="solve_failed";
    ChartSolveResult mapped;
    try {mapped=ctx.recovery_solver?ctx.recovery_solver(request):solve_lscm(request);}
    catch(const ContractError& error) {
        if(error.code==ErrorCode::cancelled || error.code==ErrorCode::constraint_conflict || error.code==ErrorCode::memory_limit)throw;
        return std::nullopt;
    } catch(const std::runtime_error&) {return std::nullopt;}
    diagnostic["reason"]="invalid_mapping";
    if(stopped() || mapped.corner_ids!=result.corners || mapped.corner_uv.size()!=result.corners.size() ||
        !std::all_of(mapped.corner_uv.begin(),mapped.corner_uv.end(),finite))return std::nullopt;
    std::vector<Vec2> trial(ctx.mesh->input().corner_vertices.size(),Vec2{NAN,NAN});
    for(std::size_t i=0;i<result.corners.size();++i)trial[result.corners[i]]=mapped.corner_uv[i];
    const std::vector<MeshIndex> one{chart};
    const auto damaged=defects(geometry,trial,one);
    diagnostic["damaged_faces"]=damaged;diagnostic["reason"]="too_many_invalid_faces";
    if(damaged.size()>std::max<std::size_t>(3,faces.size()/10))return std::nullopt;
    const QualitySnapshot before(geometry,CornerUvView{trial},{},std::span<const MeshIndex>(one));
    if(!damaged.empty())isolate(*ctx.mesh,trial,damaged,before,result.corners);
    if(ctx.options->repair.storage_precision=="float32")for(auto c:result.corners)for(auto& x:trial[c]) {
        if(std::abs(x)>std::numeric_limits<float>::max()){diagnostic["reason"]="storage_overflow";return std::nullopt;}
        x=static_cast<double>(static_cast<float>(x));
    }
    auto seams=state.seams;
    for(auto f:damaged)for(auto c:ctx.mesh->face_corners(f))seams[ctx.mesh->tables().corner_edge[c]]=1;
    for(std::size_t e=0;e<seams.size();++e) {
        if(ctx.mesh->tables().edge_boundary[e])seams[e]=0;
        if(ctx.recovery_locked && (*ctx.recovery_locked)[e] && seams[e]!=state.seams[e]){diagnostic["reason"]="protected_cut";return std::nullopt;}
    }
    const auto cut=cut_surface(ctx.mesh,seams);const auto after_geometry=quality_geometry(cut);
    std::set<MeshIndex> chosen;for(auto f:faces)chosen.insert(cut->tables().chart_of_face[f]);
    const std::vector<MeshIndex> charts(chosen.begin(),chosen.end());
    const auto remaining=defects(after_geometry,trial,charts);
    diagnostic["reason"]="invalid_after_isolation";diagnostic["remaining_invalid_faces"]=remaining;
    if(stopped() || !remaining.empty())return std::nullopt;
    const QualitySnapshot after(after_geometry,CornerUvView{trial},{},std::span<const MeshIndex>(charts));
    result.anisotropy_p95=face_weighted_p95(after);
    diagnostic["reason"]="stretch_budget";diagnostic["anisotropy_p95"]=result.anisotropy_p95;
    if(!std::isfinite(result.anisotropy_p95) || result.anisotropy_p95>2.0 || stopped())return std::nullopt;
    for(auto c:result.corners)result.values.push_back(trial[c]);
    result.isolated_faces.assign(damaged.begin(),damaged.end());diagnostic["accepted"]=true;diagnostic["reason"]="recovered";return result;
}
}
