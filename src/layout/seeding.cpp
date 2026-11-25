// SPDX-License-Identifier: GPL-3.0-or-later
#include "cleanseams/layout.hpp"
#include "cleanseams/cache.hpp"
#include "cleanseams/axis.hpp"
#include "cleanseams/panels.hpp"
#include "cleanseams/errors.hpp"
#include "cleanseams/detail/numpy_compat.hpp"
#include <algorithm>
#include <cmath>
#include <set>
namespace cleanseams {
namespace {
std::vector<MeshIndex> ones(std::span<const std::uint8_t> mask) {
    std::vector<MeshIndex> out;for(std::size_t i=0;i<mask.size();++i)if(mask[i])out.push_back(static_cast<MeshIndex>(i));return out;
}
void emit(LayoutTrace* trace,nlohmann::json event) {if(trace)trace->events.push_back(std::move(event));}
std::size_t face_size(const MeshSnapshot& m,MeshIndex f) {return static_cast<std::size_t>(m.input().face_offsets[f+1]-m.input().face_offsets[f]);}
nlohmann::json vec(const Vec3& v) {return nlohmann::json::array({v[0],v[1],v[2]});}
}
nlohmann::json LayoutReport::as_dict() const {
    nlohmann::json j{{"timed_out",timed_out},{"tier",tier},{"seeds",seeds},{"dijkstra_s",dijkstra_s},{"paired_sources",paired_sources},
        {"sources",sources},{"seeds_added",seeds_added},{"trace_passes",trace_passes},{"patches",patches},{"charts",charts},{"rounds",rounds},
        {"merges_accepted",merges_accepted},{"merges_rejected_proxy",merges_rejected_proxy},{"merges_rejected_topology",merges_rejected_topology},
        {"merges_rejected_distortion",merges_rejected_distortion},{"merges_rejected_compactness",merges_rejected_compactness},
        {"merges_rejected_feature",merges_rejected_feature},{"merges_rejected_strip",merges_rejected_strip},
        {"merges_rejected_isoperimetric",merges_rejected_isoperimetric},
        {"seam_exposure",std::isfinite(seam_exposure)?nlohmann::json(seam_exposure):nlohmann::json(nullptr)},
        {"seeds_by_exposure",seeds_by_exposure},{"revolution_cuts",revolution_cuts},{"shoulder_edges",shoulder_edges},
        {"single_arc_merges",single_arc_merges},{"self_cuts_removed",self_cuts_removed},{"spurs_removed",spurs_removed},{"spurs_kept",spurs_kept},
        {"feature_edges",feature_edges},{"feature_chains",feature_chains},{"open_features",open_features},{"seam_edges",seam_edges},
        {"seam_length",seam_length},{"non_disk_charts",non_disk_charts},{"notes",notes},{"stop_reasons",stop_reasons},{"stop_pairs",stop_pairs},
        {"last_round_rejections",last_round_rejections},{"splits",splits},{"absorbed",absorbed},{"absorbed_relaxed",absorbed_relaxed},
        {"absorbed_exempt",absorbed_exempt},{"unsupported_absorbed",unsupported_absorbed},{"premerged_nonquad",premerged_nonquad},
        {"tiny_remaining",tiny_remaining},{"absorbed_faces",absorbed_faces},{"exempt_faces",exempt_faces}};
    return j;
}
bool topological_disk(const CutSurface& cs,MeshIndex chart) {
    const auto& t=cs.tables();return t.chart_euler[chart]==1 && t.chart_n_loops[chart]==1 && !t.chart_pinched[chart];
}
bool large_planar_ngons(const MeshSnapshot& m,std::span<const MeshIndex> faces,std::int64_t minimum) {
    if(faces.empty())return false;
    std::int64_t complexity=0;
    for(auto f:faces) {
        if(face_size(m,f)<=4 || m.geometry().face_nonplanar[f])return false;
        complexity+=static_cast<std::int64_t>(face_size(m,f))-2;
    }
    return complexity>=std::max<std::int64_t>(8,minimum);
}
bool small_face_count(const MeshSnapshot& m,std::span<const MeshIndex> faces,std::int64_t minimum) {
    return static_cast<std::int64_t>(faces.size())<minimum && !large_planar_ngons(m,faces,minimum);
}
DiskLayout premerge_single_nonquad(const std::shared_ptr<const MeshSnapshot>& mesh,std::shared_ptr<const CutSurface> cs,LayoutReport& report,
                                   std::optional<std::span<const std::uint8_t>> locked_in,const Deadline& deadline,std::int64_t minimum_faces,
                                   double minimum_fraction,LayoutTrace* trace,const PinMap& pins) {
    const auto& m=*mesh;const auto& t=m.tables();const auto& g=m.geometry();
    const auto ne=t.edge_vertices.size();std::vector<std::uint8_t> locked(ne,0);
    if(locked_in)locked.assign(locked_in->begin(),locked_in->end());
    std::vector<MeshIndex> singletons;
    for(MeshIndex c=0;c<cs->tables().n_charts;++c) {const auto fs=cs->tables().chart_faces[c];if(fs.size()==1 && face_size(m,fs[0])!=4)singletons.push_back(fs[0]);}
    const double total_area=detail::numpy_sum(g.face_area);
    for(auto face:singletons) {
        if(deadline.expired()) {report.timed_out=true;break;}
        const std::array<MeshIndex,1> one{face};
        if(large_planar_ngons(m,one,minimum_faces) && g.face_area[face]>=minimum_fraction*total_area)continue;
        const auto ch=cs->tables().chart_of_face[face];if(cs->tables().chart_faces[ch].size()!=1)continue;
        std::vector<MeshIndex> edges;for(auto c:m.face_corners(face))edges.push_back(t.corner_edge[c]);
        std::stable_sort(edges.begin(),edges.end(),[&](auto a,auto b) {return std::pair(-g.edge_length[a],a)<std::pair(-g.edge_length[b],b);});
        for(auto e:edges) {
            if(locked[e] || !t.edge_manifold[e] || t.edge_faces[e].size()!=2)continue;
            const std::array<MeshIndex,1> glue{e};
            if(!cs->tables().cut[e] || cs->euler_after_gluing(glue)!=1)continue;
            if(pins_conflict_after_gluing(*cs,pins,glue))continue;
            auto cut=cs->tables().cut;cut[e]=0;
            auto candidate=cut_surface(mesh,cut);
            if(candidate->tables().chart_is_disk[candidate->tables().chart_of_face[face]]) {cs=candidate;++report.premerged_nonquad;break;}
        }
    }
    DiskLayout out;out.cut=cs->tables().cut;out.cs=cs;
    emit(trace,{{"stage","premerge_single_nonquad"},{"cut",ones(out.cut)},{"charts",out.cs->tables().n_charts}});
    return out;
}
std::vector<double> seam_cost_field(const MeshSnapshot& m,const LayoutOptions& o) {
    if(!o.use_exposure)return {};
    SeamCostOptions options;options.feature_angle=o.prefer_features?60.0:o.feature_angle;options.hide_dir=o.hide_dir;
    const auto field=exposure_field(m,static_cast<std::uint64_t>(o.exposure_K));
    return seam_cost(m,field,options).cost;
}
SeededLayout seed_layout(const std::shared_ptr<const MeshSnapshot>& mesh,const LayoutOptions& o,const SeedingInput& input,LayoutTrace* trace) {
    if(!mesh)throw ContractError(ErrorCode::invalid_input,"layout.mesh","missing mesh snapshot");
    if(o.seeding!="auto" && o.seeding!="motorcycle" && o.seeding!="voronoi")
        throw ContractError(ErrorCode::invalid_option,"layout.seeding","seeding must be auto, motorcycle or voronoi");
    const auto& m=*mesh;const auto& t=m.tables();const auto ne=t.edge_vertices.size(),nf=m.input().face_offsets.size()-1,nv=m.input().positions.size();
    const auto& deadline=input.deadline;
    SeededLayout s;auto& report=s.report;
    for(const auto* mask:{input.initial_cut?&*input.initial_cut:nullptr,input.locked?&*input.locked:nullptr})
        if(mask && mask->size()!=ne)throw ContractError(ErrorCode::invalid_input,"layout.cut","edge mask does not match the mesh edges");
    s.refine_locked=input.locked?*input.locked:EdgeMask(ne,0);s.locked=input.locked;
    if(input.cost)s.cost=*input.cost;
    else if(!deadline.expired())s.cost=seam_cost_field(m,o);
    const CostField cost=s.cost.empty()?CostField():CostField(std::span<const double>(s.cost));
    EdgeMask cut=input.initial_cut?*input.initial_cut:EdgeMask(ne,0);
    s.user_cut=cut;for(std::size_t e=0;e<ne;++e)s.user_cut[e]=s.user_cut[e] || t.edge_boundary[e] || (input.locked && (*input.locked)[e]);
    auto feats=feature_cuts(mesh,o.feature_angle,o.min_feature_chain,o.feature_scale,o.feature_nms_ratio,o.max_extend,
                            std::span<const std::uint8_t>(cut),deadline,&report,trace);
    report.feature_edges=static_cast<std::uint64_t>(std::count(feats.begin(),feats.end(),1));
    bool cad=false;
    if(o.feature_angle>0 && !deadline.expired()) {
        if(auto proposal=cad_panel_cuts(mesh,cut,cost,deadline,o.seeding!="motorcycle")) {
            cut=std::move(proposal->cut);feats=std::move(proposal->rims);cad=true;
            if(!s.locked)s.locked=feats;else for(std::size_t e=0;e<ne;++e)(*s.locked)[e]=(*s.locked)[e] || feats[e];
            report.notes.push_back({{"kind","cad_panels"},{"panels",proposal->panels},{"hole_cut_edges",proposal->hole_cut_edges},{"bevel_bands",proposal->bevel_bands}});
            report.notes.push_back("coarse CAD panels keep their boundaries despite the minimum face count");
        }
    }
    std::vector<MeshIndex> revolution_faces;
    if(o.revolution_residual>0 && !deadline.expired()) {
        auto lathe=lathe_profile_cuts(mesh,cut,cost,o.hide_dir,deadline);
        if(!lathe.handled_faces.empty()) {
            revolution_faces=lathe.handled_faces;
            for(auto f:revolution_faces)for(auto c:m.face_corners(f))feats[t.corner_edge[c]]=0;
            for(std::size_t e=0;e<ne;++e) {feats[e]=feats[e] || lathe.protected_rims[e];cut[e]=cut[e] || lathe.cut[e];}
            for(const auto& r:lathe.records)
                report.revolution_cuts.push_back({{"kind","lathe_profile"},{"charts",r.charts},{"rows",r.rows},{"rim_rows",r.rim_rows},
                    {"residual",r.residual},{"direction",vec(r.direction)}});
        }
    }
    for(std::size_t e=0;e<ne;++e)cut[e]=cut[e] || feats[e] || t.edge_boundary[e];
    if(o.revolution_residual>0 && !deadline.expired()) {
        auto revolution=revolution_cuts(mesh,cut,o.revolution_residual,cost,o.hide_dir,deadline);
        cut=std::move(revolution.cut);
        for(const auto& r:revolution.records)
            report.revolution_cuts.push_back({{"seed",{r.seed.vertex,r.seed.edge}},{"edges",r.edges},{"residual",r.residual},{"pitch",r.pitch},
                {"direction",vec(r.direction)},{"point",vec(r.point)}});
        for(std::size_t e=0;e<ne;++e)feats[e]=feats[e] || revolution.shoulders[e];
        report.shoulder_edges=static_cast<std::uint64_t>(std::count(revolution.shoulders.begin(),revolution.shoulders.end(),1));
    }
    if(o.lock_features && std::any_of(feats.begin(),feats.end(),[](auto x) {return x!=0;})) {
        if(!s.locked)s.locked=feats;else for(std::size_t e=0;e<ne;++e)(*s.locked)[e]=(*s.locked)[e] || feats[e];
        for(std::size_t e=0;e<ne;++e)s.refine_locked[e]=s.refine_locked[e] || feats[e];
    }
    auto sources=irregular_sources(m,cut,o.relaxed_motorcycles,trace);
    if(!report.feature_panel_faces.empty() || !revolution_faces.empty() || cad) {
        std::vector<std::uint8_t> panel(nf,0);
        for(auto f:report.feature_panel_faces)panel[f]=1;for(auto f:revolution_faces)panel[f]=1;
        if(cad)std::fill(panel.begin(),panel.end(),1);
        for(MeshIndex v=0;v<nv;++v) {
            if(!sources[v])continue;const auto faces=t.vert_faces[v];
            if(std::all_of(faces.begin(),faces.end(),[&](auto f) {return panel[f]!=0;}))sources[v]=0;
        }
        report.notes.push_back("recognised panels seed the layout without extra motorcycle cuts");
    }
    bool all_quads=true;std::size_t quads=0;for(MeshIndex f=0;f<nf;++f) {const bool q=face_size(m,f)==4;all_quads=all_quads && q;quads+=q;}
    if(!all_quads) {
        report.tier="quad_dominant";
        if(o.relaxed_motorcycles) {auto paired=pair_sources(m,sources,cut,trace);sources=std::move(paired.sources);cut=std::move(paired.cut);report.paired_sources=paired.count;}
    }
    report.sources=static_cast<std::uint64_t>(std::count(sources.begin(),sources.end(),1));
    const auto regular=static_cast<double>(std::count(t.vert_regular.begin(),t.vert_regular.end(),1));
    bool segmented=!cad && (o.seeding=="voronoi" || (o.seeding=="auto" && static_cast<double>(quads)/static_cast<double>(nf)<.6 &&
                                                     regular/static_cast<double>(nv)<.3));
    const auto before_tracing=cut;std::shared_ptr<const CutSurface> cs;
    if(!segmented) {
        auto disks=make_all_disks(mesh,cut,sources,report,12,cost,deadline,o.relaxed_motorcycles,trace);
        cut=std::move(disks.cut);cs=std::move(disks.cs);
        segmented=o.seeding=="auto" && static_cast<std::int64_t>(cs->tables().n_charts)>o.max_patches && !all_quads;
    }
    if(segmented) {
        const auto& barriers=o.lock_features?before_tracing:s.user_cut;
        if(input.progress)input.progress({{"stage","generate"},{"status","running"},{"msg","growing connected face regions"}});
        auto regions=segment_mesh(mesh,std::span<const std::uint8_t>(barriers),o.seed_density,o.max_patches,o.segment_iterations,cost,deadline,std::nullopt,trace);
        cut=std::move(regions.cut);cs=std::move(regions.cs);report.tier="triangle_regions";report.seeds=regions.seeds;report.dijkstra_s=regions.dijkstra_s;
        report.notes.push_back("connected dual-graph regions replace motorcycle tracing");
    }
    const auto locked_span=s.locked?std::optional<std::span<const std::uint8_t>>(std::span<const std::uint8_t>(*s.locked)):std::nullopt;
    const auto pins=input.effective_pins();
    for(const auto& [c,uv]:pins)
        if(c>=t.corner_face.size())throw ContractError(ErrorCode::invalid_input,"layout.pins","pinned corner out of range");
    auto premerged=premerge_single_nonquad(mesh,cs,report,locked_span,deadline,o.min_island_faces,o.min_island_fraction,trace,pins);
    cut=std::move(premerged.cut);cs=std::move(premerged.cs);
    if(input.evaluate && *input.evaluate)s.evaluate=*input.evaluate;
    else if(o.evaluate_true_distortion) {EvaluatorSettings settings;settings.policy=input.policy;settings.cancel=input.cancel;s.evaluate=default_evaluator(mesh,o,settings);}
    if(s.evaluate && deadline.at()) {
        // T24: the source wrapper evaluated one chart per call until the deadline passed.
        // Bounded batches keep the deadline and cancellation checks between batches
        // without forcing singletons; charts of stopped batches stay unevaluated.
        BatchedEvaluation batched;batched.scheduling=input.scheduling;batched.deadline=deadline;batched.cancel=input.cancel;
        batched.log=input.batches;batched.timed_out=s.evaluation_timed_out;batched.stage="generate";
        s.evaluate=batched_evaluator(std::move(s.evaluate),std::move(batched));
    }
    auto split=split_over_budget(mesh,cut,cs,sources,o,report,s.evaluate?&s.evaluate:nullptr,cost,deadline,segmented,trace);
    s.cut=std::move(split.cut);s.cs=std::move(split.cs);s.sources=std::move(split.sources);
    report.patches=s.cs->tables().n_charts;
    if(s.locked)for(std::size_t e=0;e<ne;++e)s.cut[e]=s.cut[e] || (*s.locked)[e];
    const auto& ct=s.cs->tables();
    s.next_stage=std::any_of(ct.chart_unsupported.begin(),ct.chart_unsupported.end(),[](auto x) {return x!=0;}) && !segmented?"absorb_tiny":"merge_patches";
    s.features=std::move(feats);s.segmented=segmented;
    report.timed_out=report.timed_out || *s.evaluation_timed_out;
    return s;
}
GeneratedLayout generate_seams(const std::shared_ptr<const MeshSnapshot>& mesh,const LayoutOptions& o,const SeedingInput& input,LayoutTrace* trace) {
    auto s=seed_layout(mesh,o,input,trace);
    const auto& m=*mesh;const auto& t=m.tables();const auto& g=m.geometry();const auto ne=t.edge_vertices.size();
    GeneratedLayout out;auto& report=out.report;report=std::move(s.report);
    const auto& deadline=input.deadline;const ChartEvaluator* evaluate=s.evaluate?&s.evaluate:nullptr;
    const CostField cost=s.cost.empty()?CostField():CostField(std::span<const double>(s.cost));
    auto cut=std::move(s.cut);auto cs=std::move(s.cs);
    if(s.next_stage=="absorb_tiny") {
        AbsorbRequest request;request.locked=s.locked;request.deadline=deadline;request.unsupported_only=true;request.pins=input.effective_pins();
        request.progress=input.progress;
        auto absorbed=absorb_tiny(mesh,cs,o,report,evaluate,request,trace);cut=std::move(absorbed.cut);cs=std::move(absorbed.cs);
    }
    MergeConstraints constraints;constraints.locked=s.locked;constraints.cost=cost;constraints.deadline=deadline;constraints.pins=input.effective_pins();
    constraints.progress=input.progress;
    if(!s.segmented || o.lock_features)constraints.features=s.features;
    auto merged=merge_patches(mesh,cut,o,report,evaluate,constraints,trace);
    AbsorbRequest request;request.locked=s.locked;request.deadline=deadline;request.pins=input.effective_pins();request.progress=input.progress;
    auto absorbed=absorb_tiny(mesh,merged.cs,o,report,evaluate,request,trace);cut=std::move(absorbed.cut);cs=std::move(absorbed.cs);
    if(o.spur_removal && !deadline.expired()) {
        const auto locked=s.locked?std::optional<std::span<const std::uint8_t>>(std::span<const std::uint8_t>(*s.locked)):std::nullopt;
        auto swept=final_spur_sweep(mesh,cut,cs,o,report,evaluate,locked,trace,input.effective_pins());cut=std::move(swept.cut);cs=std::move(swept.cs);
    }
    out.seams.assign(ne,0);for(std::size_t e=0;e<ne;++e)out.seams[e]=cut[e] && !t.edge_boundary[e];
    // CAD panel borders block merging, but their interior route may still be improved.
    // Only artist locks and an explicit lock-features request freeze refinement edges.
    report.protected_seams=s.refine_locked;
    report.charts=cs->tables().n_charts;
    if(evaluate)report.uv=assemble_uv(m,cs,report,*evaluate,trace);
    report.timed_out=report.timed_out || *s.evaluation_timed_out || deadline.expired();
    std::vector<double> lengths,weighted;
    for(std::size_t e=0;e<ne;++e)if(out.seams[e]) {lengths.push_back(g.edge_length[e]);if(cost)weighted.push_back((*cost)[e]*g.edge_length[e]);}
    report.seam_edges=lengths.size();report.seam_length=detail::numpy_sum(lengths);
    const auto& disks=cs->tables().chart_is_disk;report.non_disk_charts=static_cast<std::uint64_t>(std::count(disks.begin(),disks.end(),0));
    if(cost && !lengths.empty())report.seam_exposure=detail::numpy_sum(weighted)/std::max(detail::numpy_sum(lengths),1e-12);
    out.cs=std::move(cs);
    return out;
}
}
