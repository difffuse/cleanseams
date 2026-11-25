// SPDX-License-Identifier: GPL-3.0-or-later
#include "cleanseams/layout_postprocess.hpp"
#include "cleanseams/cache.hpp"
#include "cleanseams/metrics.hpp"
#include "cleanseams/errors.hpp"
#include "layout_postprocess_internal.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
namespace cleanseams {
void check_layout_operation(const PipelineInput& in) {
    auto require=[](bool ok,const char* field,const char* why) {
        if(!ok)throw ContractError(ErrorCode::invalid_input,field,why);
    };
    require(in.operation=="pipeline" || in.operation=="postprocess_layout","operation","unknown operation");
    require(bool(in.mesh),"mesh","missing mesh snapshot");
    const auto nf=in.mesh->input().face_offsets.size()-1,nc=in.mesh->input().corner_vertices.size();
    require(in.selected_faces.empty() || in.selected_faces.size()==nf,"selected_faces","mask must match original faces");
    for(auto v:in.selected_faces)require(v<=1,"selected_faces","mask must contain only zero or one");
    const auto& r=in.options.repair;
    if((r.relax_faces_seconds!=0 && (in.operation!="postprocess_layout" || in.options.policy!="repair_v1")) ||
       (r.region_recovery_seconds!=0 && in.options.policy!="repair_v1") ||
       (r.singleton_seconds!=0 && (in.operation!="postprocess_layout" || in.options.policy!="repair_v1")) ||
       (r.join_seconds!=0 && (in.operation!="postprocess_layout" || in.options.policy!="repair_v1")) ||
       (r.host_proposal_seconds!=0 && (in.operation!="postprocess_layout" || in.options.policy!="repair_v1")))
        throw ContractError(ErrorCode::invalid_option,"repair","nonzero recovery allowances are not supported by this candidate");
    if(r.require_open_protected_cuts && (in.operation!="postprocess_layout" || in.options.policy!="repair_v1"))
        throw ContractError(ErrorCode::invalid_option,"repair.require_open_protected_cuts","requires postprocess_layout and repair_v1");
    if(in.operation!="postprocess_layout")return;
    pipeline_options_json(in.options);
    require(in.options.policy=="repair_v1","options.policy","postprocess_layout requires repair_v1");
    require(in.uv.has_value(),"uv","postprocess_layout requires supplied UVs");
    require(in.uv->size()==nc,"uv","UV count must equal corner count");
    for(std::size_t f=0;f<nf;++f)for(auto c:in.mesh->face_corners(static_cast<MeshIndex>(f))) {
        const auto& p=(*in.uv)[c];
        const bool recover=r.region_recovery_seconds>0 && !in.selected_faces.empty() && in.selected_faces[f];
        require(recover || (std::isfinite(p[0]) && std::isfinite(p[1])),"uv","nonfinite supplied UVs require selected region recovery");
    }
    const auto ne=in.mesh->tables().edge_vertices.size();
    require(in.seams.size()==ne,"seams","mask must match mesh edges");
    for(auto v:in.seams)require(v<=1,"seams","mask must contain only zero or one");
    if(in.locked) {
        require(in.locked->size()==ne,"locked","mask must match mesh edges");
        for(auto v:*in.locked)require(v<=1,"locked","mask must contain only zero or one");
    }
    for(const auto& [c,p]:in.pins) {
        require(c<nc,"pins","corner out of range");
        require(std::isfinite(p[0]) && std::isfinite(p[1]),"pins","pin coordinates must be finite");
    }
}
PipelineResult postprocess_layout(const PipelineInput& in,const ProgressSink& progress,const Control& control,const LayoutPostprocessServices& services) {
    MemoryBudgetScope budget(in.memory_budget);
    check_layout_operation(in);
    check_run_memory(table_counts(*in.mesh),in.memory_budget);
    const auto deadline=Deadline::after(in.options.deadline_s);
    bool timed_out=control.check() || deadline.expired();
    if(progress)progress({{"stage","postprocess_layout"},{"status","started"}});
    PipelineResult out;out.uv=*in.uv;out.seams=in.seams;
    // A supplied cut remains a cut even when endpoint coordinates coincide. Additional
    // UV discontinuities use exact source-edge correspondence, never coordinate welding.
    auto discontinuities=seams_from_uv(*in.mesh,out.uv,0.0);
    // Unlike the legacy metric's NaN comparison, unknown endpoints cannot prove
    // continuity for local recovery. Preserve them as separate sectors.
    for(std::size_t e=0;e<discontinuities.size();++e)if(in.mesh->tables().edge_corners[e].size()==2)
        for(auto c:in.mesh->tables().edge_corners[e])for(auto v:{c,in.mesh->corner_next(c)})
            if(!std::isfinite(out.uv[v][0]) || !std::isfinite(out.uv[v][1]))discontinuities[e]=1;
    bool protected_valid=true;
    for(std::size_t e=0;e<out.seams.size();++e) {
        out.seams[e]=out.seams[e] || discontinuities[e];
        if(in.locked && (*in.locked)[e] && out.seams[e]!=in.seams[e])protected_valid=false;
    }
    out.cut=cut_surface(in.mesh,out.seams);
    const bool single=in.options.repair.storage_precision=="float32";
    auto rounded=[&](double x) {
        if(!single)return x;
        if(std::abs(x)>std::numeric_limits<float>::max())return std::numeric_limits<double>::quiet_NaN();
        return static_cast<double>(static_cast<float>(x));
    };
    const auto& ct=out.cut->tables();
    std::int64_t global_vote=0;
    auto orientation_uv=out.uv;for(auto& p:orientation_uv)for(auto& x:p)if(!std::isfinite(x))x=0;
    for(std::size_t t=0;t<in.mesh->triangulation().corners.size();++t) {
        const auto ids=in.mesh->triangulation().corners[t];const auto a=orientation_uv[ids[0]],b=orientation_uv[ids[1]],c=orientation_uv[ids[2]];
        const double area=(b[0]-a[0])*(c[1]-a[1])-(b[1]-a[1])*(c[0]-a[0]);
        if(std::isfinite(area))global_vote+=(area>0)-(area<0);
    }
    // Invalid-face classification uses one global signed-triangle vote,
    // with positive tie breaking. Keep it immutable across all partition changes.
    const std::vector<std::int8_t> orientation(ct.chart_of_face.size(),global_vote<0?-1:1);
    struct Validation {
        bool valid=false;
        std::vector<MeshIndex> invalid,defect_faces;
        std::vector<QualityReport> quality;
        nlohmann::json report;
    };
    auto validate=[&](const std::vector<Vec2>& uv,const std::shared_ptr<const CutSurface>& surface,bool summarize=true,
                      std::optional<std::span<const MeshIndex>> selected=std::nullopt,bool include_exempt_overlap=false) {
        control.check();auto stored=uv;
        std::vector<std::uint8_t> changed(orientation.size(),0);std::vector<MeshIndex> coverage_faces;
        for(MeshIndex f=0;f<changed.size();++f)for(auto c:in.mesh->face_corners(f)) {
            for(std::size_t d=0;d<2;++d)if(uv[c][d]!=(*in.uv)[c][d] && !(std::isnan(uv[c][d]) && std::isnan((*in.uv)[c][d])))changed[f]=1;
            if(changed[f])break;
        }
        for(MeshIndex f=0;f<changed.size();++f)if(changed[f])coverage_faces.push_back(f);
        for(auto& p:stored)for(auto& x:p)x=rounded(x);
        bool pins_valid=true;
        for(const auto& [c,p]:in.pins)for(std::size_t d=0;d<2;++d)
            pins_valid=pins_valid && uv[c][d]==p[d] && stored[c][d]==rounded(p[d]);
        auto normalized=stored;
        for(std::size_t c=0;c<normalized.size();++c)normalized[c][1]*=orientation[in.mesh->tables().corner_face[c]];
        const auto geometry=quality_geometry(surface);const auto& charts=surface->tables();
        const bool input_topology_supported=std::none_of(charts.chart_unsupported.begin(),charts.chart_unsupported.end(),[](auto x){return x!=0;});
        const bool cut_sectors_valid=std::none_of(charts.chart_pinched.begin(),charts.chart_pinched.end(),[](auto x){return x!=0;});
        // A split original vertex fan is still unsupported input topology, but
        // separate unpinched sectors can carry a strictly valid supplied UV map.
        const bool topology_valid=input_topology_supported && cut_sectors_valid;
        bool cuts_valid=protected_valid;
        for(std::size_t e=0;e<in.seams.size();++e)if(in.locked && (*in.locked)[e] && charts.cut[e]!=bool(in.seams[e] || in.mesh->tables().edge_boundary[e]))cuts_valid=false;
        const bool logical_cuts_valid=cuts_valid;
        std::vector<MeshIndex> closed_edges,recovery_faces;
        if(in.options.repair.require_open_protected_cuts && in.locked) {
            const auto& mt=in.mesh->tables();std::uint64_t pairs=0;
            std::vector<std::uint8_t> recover(orientation.size(),0);
            for(MeshIndex e=0;e<in.seams.size();++e)if(in.seams[e] && (*in.locked)[e]) {
                // Stable original-corner ordering matches the source's second-side
                // isolation rule. Do not materialize the quadratic all-pairs array.
                const auto corners=mt.edge_corners[e];bool closed_edge=false;
                for(std::size_t i=0;i<corners.size();++i)for(std::size_t j=i+1;j<corners.size();++j) {
                    if((pairs++ & 1023)==0)control.check();
                    auto endpoints=[&](MeshIndex c) {
                        std::array<MeshIndex,2> ids{c,in.mesh->corner_next(c)};
                        if(in.mesh->input().corner_vertices[ids[0]]>in.mesh->input().corner_vertices[ids[1]])std::swap(ids[0],ids[1]);
                        return ids;
                    };
                    const auto a=endpoints(corners[i]),b=endpoints(corners[j]);bool closed=true;
                    for(std::size_t k=0;k<2;++k)for(std::size_t d=0;d<2;++d)
                        closed=closed && std::abs(stored[a[k]][d]-stored[b[k]][d])<=1e-7;
                    if(closed) {closed_edge=true;recover[mt.corner_face[corners[j]]]=1;}
                }
                if(closed_edge)closed_edges.push_back(e);
            }
            for(MeshIndex f=0;f<recover.size();++f)if(recover[f])recovery_faces.push_back(f);
        }
        cuts_valid=cuts_valid && closed_edges.empty();
        const QualityEvaluation quality(geometry,CornerUvView{normalized},{},selected);
        OverlapOptions overlap_options;overlap_options.area_tolerance=0;overlap_options.ignore_input_degenerate=false;
        for(MeshIndex f=0;f<changed.size();++f)if(!include_exempt_overlap && !changed[f] && in.mesh->geometry().face_degenerate_kind[f])overlap_options.ignored_faces.push_back(f);
        const auto overlap=chart_overlaps(*geometry,CornerUvView{stored},overlap_options,selected);
        Validation result;std::vector<std::uint8_t> defects(changed.size(),0);
        std::vector<MeshIndex> checked_charts;
        if(selected)checked_charts.assign(selected->begin(),selected->end());
        else for(MeshIndex ch=0;ch<charts.n_charts;++ch)checked_charts.push_back(ch);
        // Selected snapshots retain original IDs but have compact measurement arrays.
        // Never index those arrays with an original face or triangle index.
        std::vector<double> polygon_area(changed.size(),NAN);
        for(const auto& polygon:quality.snapshot().polygons())polygon_area[polygon.face]=polygon.raw_signed_area;
        std::vector<std::uint8_t> triangles_valid(charts.n_charts,1);
        for(const auto& triangle:quality.snapshot().triangles()) {
            if((!triangle.exempt || changed[triangle.face]) && (!triangle.finite || triangle.raw_area2<=0))triangles_valid[triangle.chart]=0;
            if((include_exempt_overlap || !triangle.exempt || changed[triangle.face]) && (!triangle.finite || triangle.raw_area2<=0))defects[triangle.face]=1;
        }
        for(auto f:overlap.overlapping_faces)defects[f]=1;
        for(auto ch:checked_charts) {
            bool valid=quality.chart_valid(ch) && !overlap.counts[ch] && overlap.tested[ch];
            for(auto f:quality.offending_faces(ch))defects[f]=1;
            for(auto f:charts.chart_faces[ch])for(auto c:in.mesh->face_corners(f))
                if(!std::isfinite(stored[c][0]) || !std::isfinite(stored[c][1])){valid=false;defects[f]=1;}
            for(auto f:charts.chart_faces[ch])if(in.mesh->geometry().face_degenerate_kind[f]) {
                std::vector<Vec2> polygon;for(auto c:in.mesh->face_corners(f))polygon.push_back(normalized[c]);
                // Unreliable 3D distortion does not exempt a crossed UV boundary.
                // Report inherited crossings too, so local recovery can select
                // them before this layout is submitted for packing.
                if(polygon_self_intersects(polygon,0) || (changed[f] && !(polygon_area[f]>0))){valid=false;defects[f]=1;}
            }
            // Polygon signed area alone cannot certify the fixed interior triangles.
            valid=valid && triangles_valid[ch];
            if(!valid)result.invalid.push_back(ch);
            if(summarize)result.quality.push_back(quality.snapshot().report(ch));
        }
        for(MeshIndex f=0;f<defects.size();++f)if(defects[f])result.defect_faces.push_back(f);
        result.valid=result.invalid.empty() && pins_valid && cuts_valid && cut_sectors_valid;
        result.report=summarize?quality.metrics():nlohmann::json::object();result.report["valid"]=result.valid;
        result.report["storage_precision"]=in.options.repair.storage_precision;
        result.report["validation_scope"]=selected?"selected_charts":"whole_layout";
        result.report["validated_face_count"]=quality.snapshot().polygons().size();
        result.report["topology_valid"]=topology_valid;result.report["pins_valid"]=pins_valid;
        result.report["input_topology_supported"]=input_topology_supported;
        result.report["cut_sectors_valid"]=cut_sectors_valid;
        result.report["uv_coverage_valid"]=result.invalid.empty();
        std::vector<MeshIndex> unsupported_vertices;
        const auto& topology=in.mesh->tables();
        for(MeshIndex v=0;v<topology.fans.size();++v)if(!topology.fans[v] && !topology.vert_faces[v].empty())unsupported_vertices.push_back(v);
        result.report["unsupported_input_vertices"]=unsupported_vertices;
        result.report["pinched_sectors"]=charts.pinched_sectors;
        result.report["logical_protected_cuts_valid"]=logical_cuts_valid;
        result.report["protected_uv_cuts_required"]=in.options.repair.require_open_protected_cuts;
        result.report["protected_uv_cuts_valid"]=closed_edges.empty();
        result.report["closed_protected_uv_edges"]=closed_edges;
        result.report["protected_cut_recovery_faces"]=recovery_faces;
        result.report["protected_cuts_valid"]=cuts_valid;result.report["strict_overlap"]=overlap.to_json();
        result.report["reference_orientation_of_face"]=orientation;
        result.report["coverage_checked_faces"]=coverage_faces;
        result.report["invalid_charts"]=result.invalid;
        result.report["coverage_defect_faces"]=result.defect_faces;
        return result;
    };
    auto validation=validate(out.uv,out.cut);
    nlohmann::json relax={{"accepted",false},{"status","skipped"},{"reason","zero_optional_budget"}};
    const bool requested=in.options.repair.relax_faces_seconds>0;
    if(requested) {
        const auto local=Deadline::after(in.options.repair.relax_faces_seconds);
        bool phase_expired=false;
        auto stopped=[&] {
            const bool expired=control.check() || deadline.expired() || local.expired() || phase_expired;
            if(expired){phase_expired=true;timed_out=true;relax["status"]="timed_out";relax["reason"]="deadline";}
            return expired;
        };
        relax["status"]="unchanged";
        if(progress)progress({{"stage","relax_faces"},{"status","started"}});
        const bool finite_input=std::all_of(out.uv.begin(),out.uv.end(),[](auto p){return std::isfinite(p[0]) && std::isfinite(p[1]);});
        if(!finite_input)relax["reason"]="nonfinite_input";
        if(!stopped() && finite_input) {
            auto trial=repair_detail::relax_selected_faces(in,*out.cut,stopped,relax);
            if(trial && !stopped()) {
                std::vector<std::uint8_t> affected(ct.n_charts,0);
                for(std::size_t c=0;c<trial->size();++c) {
                    if((*trial)[c]!=out.uv[c])for(auto& x:(*trial)[c])x=rounded(x);
                    if((*trial)[c]!=out.uv[c])affected[ct.chart_of_corner[c]]=1;
                }
                auto candidate=validate(*trial,out.cut);
                bool quality_ok=candidate.valid;
                // A valid chart cannot buy improvement elsewhere by degrading. For an
                // invalid chart, bound the repaired P95 ratio instead of comparing NaNs.
                for(MeshIndex ch=0;ch<ct.n_charts && quality_ok;++ch) {
                    if(!affected[ch])continue;
                    const auto& before=validation.quality[ch];const auto& after=candidate.quality[ch];
                    const bool was_invalid=std::binary_search(validation.invalid.begin(),validation.invalid.end(),ch);
                    const double p95_limit=was_invalid?2.0:before.anisotropy_p95*(1+1e-6);
                    quality_ok=std::isfinite(after.anisotropy_p95) && after.anisotropy_p95<=p95_limit;
                    if(!was_invalid)quality_ok=quality_ok && after.anisotropy_max<=before.anisotropy_max*(1+1e-6);
                }
                relax["status"]="rejected";relax["reason"]=candidate.valid?"distortion_regression":"invalid_proposal";
                if(quality_ok && !stopped()) {
                    std::vector<MeshIndex> changed;
                    for(MeshIndex f=0;f+1<in.mesh->input().face_offsets.size();++f) {
                        const auto corners=in.mesh->face_corners(f);
                        if(std::any_of(corners.begin(),corners.end(),[&](auto c){return (*trial)[c]!=out.uv[c];}))changed.push_back(f);
                    }
                    if(stopped()) { /* retain the original layout and its empty change set */ }
                    else if(!changed.empty()) {
                        out.changed_faces=std::move(changed);
                        out.uv=std::move(*trial);validation=std::move(candidate);
                        relax["accepted"]=true;relax["status"]="accepted";relax["reason"]="recovered";
                    } else {relax["status"]="unchanged";relax["reason"]="no_coordinate_change";}
                }
                // Expiry after a rejected proposal must also remain visible.
                stopped();
            }
        }
        if(progress)progress({{"stage","relax_faces"},{"status",relax["status"]}});
    }
    nlohmann::json protected_recovery={{"accepted",false},{"status","skipped"},{"reason","zero_optional_budget"}};
    nlohmann::json regions={{"accepted",false},{"status","skipped"},{"reason","zero_optional_budget"}};
    if(in.options.repair.region_recovery_seconds>0) {
        const auto local=Deadline::after(in.options.repair.region_recovery_seconds);
        bool phase_expired=false;
        auto stopped=[&] {
            const bool expired=control.check() || deadline.expired() || local.expired() || phase_expired;
            if(expired){phase_expired=true;timed_out=true;regions["status"]="timed_out";regions["reason"]="deadline";}
            return expired;
        };
        regions["status"]="unchanged";
        if(progress)progress({{"stage","repair_regions"},{"status","started"}});
        protected_recovery["reason"]="no_closed_cuts";
        const auto closed_faces=validation.report.at("protected_cut_recovery_faces").get<std::vector<MeshIndex>>();
        if(!closed_faces.empty() && !stopped()) {
            auto phase=in;phase.uv=out.uv;phase.seams=out.seams;phase.selected_faces.assign(orientation.size(),0);
            for(auto f:closed_faces)phase.selected_faces[f]=1;
            const auto previous=validate(out.uv,out.cut,false,std::nullopt,true);
            protected_recovery["status"]="unchanged";
            auto trial=repair_detail::recover_selected_regions(phase,out.cut,orientation,stopped,protected_recovery,{},true);
            if(trial && !stopped()) {
                auto surface=cut_surface(in.mesh,trial->seams);
                auto candidate=validate(trial->uv,surface);const auto coverage=validate(trial->uv,surface,false,std::nullopt,true);
                bool safe=candidate.report.at("pins_valid").get<bool>() && candidate.report.at("protected_cuts_valid").get<bool>() &&
                    candidate.report.at("cut_sectors_valid").get<bool>() &&
                    std::includes(previous.defect_faces.begin(),previous.defect_faces.end(),coverage.defect_faces.begin(),coverage.defect_faces.end());
                for(auto f:closed_faces)if(std::binary_search(coverage.defect_faces.begin(),coverage.defect_faces.end(),f))safe=false;
                protected_recovery["status"]="rejected";protected_recovery["reason"]="invalid_proposal";
                protected_recovery["coverage_defects_before"]=previous.defect_faces;protected_recovery["coverage_defects_after"]=coverage.defect_faces;
                if(safe && !stopped()) {
                    out.uv=std::move(trial->uv);out.seams=std::move(trial->seams);out.cut=std::move(surface);validation=std::move(candidate);
                    out.changed_faces.clear();
                    for(MeshIndex f=0;f<orientation.size();++f)for(auto c:in.mesh->face_corners(f))if(out.uv[c]!=(*in.uv)[c]){out.changed_faces.push_back(f);break;}
                    protected_recovery["accepted"]=true;protected_recovery["status"]="accepted";protected_recovery["reason"]="recovered";
                }
            }
            stopped();
        }
        if(!closed_faces.empty() && phase_expired) {
            protected_recovery["status"]="timed_out";protected_recovery["reason"]="deadline";
        }
        if(!stopped()) {
            auto phase=in;phase.uv=out.uv;phase.seams=out.seams;
            auto trial=repair_detail::recover_selected_regions(phase,out.cut,orientation,stopped,regions,services.region_solver);
            if(trial && !stopped()) {
                const auto surface=cut_surface(in.mesh,trial->seams);auto candidate=validate(trial->uv,surface);
                regions["status"]="rejected";regions["reason"]="invalid_proposal";
                if(candidate.valid && !stopped()) {
                    out.uv=std::move(trial->uv);out.seams=std::move(trial->seams);out.cut=surface;validation=std::move(candidate);
                    out.changed_faces.clear();
                    for(MeshIndex f=0;f<orientation.size();++f)for(auto c:in.mesh->face_corners(f))if(out.uv[c]!=(*in.uv)[c]){out.changed_faces.push_back(f);break;}
                    regions["accepted"]=true;regions["status"]="accepted";regions["reason"]="recovered";
                }
                stopped();
            }
        }
        if(progress)progress({{"stage","repair_regions"},{"status",regions["status"]}});
    }
    nlohmann::json fragments={{"accepted",false},{"status","skipped"},{"reason","zero_optional_budget"}};
    if(in.options.repair.host_proposal_seconds>0) {
        const auto local=Deadline::after(in.options.repair.host_proposal_seconds);bool phase_expired=false;
        fragments["status"]="unchanged";fragments["phase_budget_seconds"]=in.options.repair.host_proposal_seconds;
        auto stopped=[&] {
            const bool expired=control.check() || deadline.expired() || local.expired() || phase_expired;
            if(expired){phase_expired=true;timed_out=true;fragments["status"]="timed_out";fragments["reason"]="deadline";}return expired;
        };
        if(progress)progress({{"stage","recover_fragments"},{"status","started"}});
        if(!stopped()) {
            auto phase=in;phase.uv=out.uv;phase.seams=out.seams;
            const auto previous=validate(out.uv,out.cut,false,std::nullopt,true);
            auto valid=[&](const repair_detail::LayoutRecovery& trial,const std::shared_ptr<const CutSurface>& joined,std::span<const MeshIndex> moved) {
                const auto checked=validate(trial.uv,joined,false,std::nullopt,true);
                if(!checked.report.at("pins_valid").get<bool>() || !checked.report.at("protected_cuts_valid").get<bool>() ||
                   !checked.report.at("cut_sectors_valid").get<bool>())return false;
                for(auto f:checked.defect_faces)if(std::binary_search(moved.begin(),moved.end(),f))return false;
                return std::includes(previous.defect_faces.begin(),previous.defect_faces.end(),checked.defect_faces.begin(),checked.defect_faces.end());
            };
            auto trial=repair_detail::recover_fragments(phase,out.cut,orientation,stopped,services,valid,fragments);
            if(trial) {
                out.uv=std::move(trial->uv);out.seams=std::move(trial->seams);out.cut=cut_surface(in.mesh,out.seams);
                validation=validate(out.uv,out.cut);out.changed_faces.clear();
                for(MeshIndex f=0;f<orientation.size();++f)for(auto c:in.mesh->face_corners(f))if(out.uv[c]!=(*in.uv)[c]){out.changed_faces.push_back(f);break;}
                fragments["accepted"]=true;if(!phase_expired)fragments["status"]="accepted";
            }
        }
        if(progress)progress({{"stage","recover_fragments"},{"status",fragments["status"]}});
    }
    nlohmann::json singleton={{"accepted",false},{"status","skipped"},{"reason","zero_optional_budget"}};
    if(in.options.repair.singleton_seconds>0) {
        const double allowance=std::min(45.0,in.options.repair.singleton_seconds);
        const auto local=Deadline::after(allowance);bool phase_expired=false;
        singleton["phase_budget_seconds"]=allowance;singleton["status"]="unchanged";
        auto stopped=[&] {
            const bool expired=control.check() || deadline.expired() || local.expired() || phase_expired;
            if(expired){phase_expired=true;timed_out=true;singleton["status"]="timed_out";}
            return expired;
        };
        if(progress)progress({{"stage","repair_singletons"},{"status","started"}});
        if(!stopped()) {
            auto phase=in;phase.uv=out.uv;phase.seams=out.seams;
            std::shared_ptr<const CutSurface> measured_incumbent;
            std::vector<std::uint8_t> old_bad(orientation.size(),0);
            auto valid=[&](const repair_detail::LayoutRecovery& incumbent,const std::shared_ptr<const CutSurface>& incumbent_cut,
                           const repair_detail::LayoutRecovery& trial,const std::shared_ptr<const CutSurface>& joined,
                           std::span<const MeshIndex> faces,MeshIndex required_face) {
                // Compare with the last accepted layout, never a rejected proposal or
                // the phase-original defects after a later proposal repaired them.
                if(measured_incumbent!=incumbent_cut) {
                    const auto previous=validate(incumbent.uv,incumbent_cut,false,std::nullopt,true);
                    std::fill(old_bad.begin(),old_bad.end(),0);for(auto f:previous.defect_faces)old_bad[f]=1;
                    measured_incumbent=incumbent_cut;
                    if(!singleton.contains("coverage_defects_before"))singleton["coverage_defects_before"]=previous.defect_faces;
                }
                std::vector<MeshIndex> charts;for(auto f:faces)charts.push_back(joined->tables().chart_of_face[f]);
                std::sort(charts.begin(),charts.end());charts.erase(std::unique(charts.begin(),charts.end()),charts.end());
                const auto checked=validate(trial.uv,joined,false,std::span<const MeshIndex>(charts),true);
                if(!checked.report.at("pins_valid").get<bool>() || !checked.report.at("protected_cuts_valid").get<bool>() ||
                   !checked.report.at("cut_sectors_valid").get<bool>())return false;
                for(auto f:checked.defect_faces)if(f==required_face || !old_bad[f])return false;
                return true;
            };
            auto trial=repair_detail::reconstruct_singletons(phase,out.cut,orientation,stopped,valid,singleton);
            if(trial) {
                out.uv=std::move(trial->uv);out.seams=std::move(trial->seams);out.cut=cut_surface(in.mesh,out.seams);
                validation=validate(out.uv,out.cut);out.changed_faces.clear();
                for(MeshIndex f=0;f<orientation.size();++f)for(auto c:in.mesh->face_corners(f))if(out.uv[c]!=(*in.uv)[c]){out.changed_faces.push_back(f);break;}
                singleton["accepted"]=true;singleton["remaining_invalid_charts"]=validation.invalid;
                singleton["coverage_defects_after"]=validate(out.uv,out.cut,false,std::nullopt,true).defect_faces;
                if(!phase_expired)singleton["status"]="accepted";
            }
            stopped();
        }
        if(progress)progress({{"stage","repair_singletons"},{"status",singleton["status"]}});
    }
    nlohmann::json joining={{"accepted",false},{"status","skipped"},{"reason","zero_optional_budget"}};
    if(in.options.repair.join_seconds>0) {
        const double allowance=std::min(in.options.repair.join_mode=="fewer"?2.0:1.0,in.options.repair.join_seconds);
        const auto local=Deadline::after(allowance);
        joining["phase_budget_seconds"]=allowance;joining["mode"]=in.options.repair.join_mode;
        bool phase_expired=false;
        auto stopped=[&] {
            const bool expired=control.check() || deadline.expired() || local.expired() || phase_expired;
            if(expired){phase_expired=true;timed_out=true;joining["status"]="timed_out";}
            return expired;
        };
        joining["status"]="unchanged";
        if(progress)progress({{"stage","join_small_islands"},{"status","started"}});
        if(!stopped()) {
            auto phase=in;phase.uv=out.uv;phase.seams=out.seams;
            joining["validation_rejections"]=nlohmann::json::array();
            const auto phase_coverage=validate(out.uv,out.cut,false,std::nullopt,true);
            joining["coverage_defects_before"]=phase_coverage.defect_faces;
            auto no_new_defects=[](const Validation& before,const Validation& after) {
                return after.report.at("pins_valid").get<bool>() && after.report.at("protected_cuts_valid").get<bool>() &&
                    after.report.at("cut_sectors_valid").get<bool>() &&
                    std::includes(before.defect_faces.begin(),before.defect_faces.end(),after.defect_faces.begin(),after.defect_faces.end());
            };
            auto valid=[&](const repair_detail::LayoutRecovery& incumbent,const std::shared_ptr<const CutSurface>& incumbent_cut,
                           const repair_detail::LayoutRecovery& trial,const std::shared_ptr<const CutSurface>& joined,std::span<const MeshIndex> moved){
                std::vector<MeshIndex> affected;
                for(auto f:moved)affected.push_back(joined->tables().chart_of_face[f]);
                std::sort(affected.begin(),affected.end());affected.erase(std::unique(affected.begin(),affected.end()),affected.end());
                // Measure every old chart participating in this merge, including
                // distant/exempt coverage, against the actual last accepted layout.
                std::vector<MeshIndex> previous_charts;
                for(auto ch:affected)for(auto f:joined->tables().chart_faces[ch])previous_charts.push_back(incumbent_cut->tables().chart_of_face[f]);
                std::sort(previous_charts.begin(),previous_charts.end());previous_charts.erase(std::unique(previous_charts.begin(),previous_charts.end()),previous_charts.end());
                const auto previous=validate(incumbent.uv,incumbent_cut,false,std::span<const MeshIndex>(previous_charts),true);
                const auto checked=validate(trial.uv,joined,false,std::span<const MeshIndex>(affected),true);
                bool accepted=no_new_defects(previous,checked);
                // Joining must not merely relocate the moving chart's existing defect.
                for(auto f:moved)if(std::binary_search(checked.defect_faces.begin(),checked.defect_faces.end(),f))accepted=false;
                joining["baseline_face_visits"]=joining.value("baseline_face_visits",std::uint64_t{0})+previous.report.at("validated_face_count").get<std::uint64_t>();
                joining["validated_face_visits"]=joining.value("validated_face_visits",std::uint64_t{0})+checked.report.at("validated_face_count").get<std::uint64_t>();
                if(!accepted && joining["validation_rejections"].size()<8) {
                    std::vector<MeshIndex> changed_exempt;
                    for(auto f:checked.report.at("coverage_checked_faces").get<std::vector<MeshIndex>>())
                        if(in.mesh->geometry().face_degenerate_kind[f])changed_exempt.push_back(f);
                    joining["validation_rejections"].push_back({{"invalid_charts",checked.invalid},
                        {"coverage_defects_before",previous.defect_faces},{"coverage_defects_after",checked.defect_faces},
                        {"changed_exempt_faces",changed_exempt},{"overlapping_faces",checked.report.at("strict_overlap").at("overlapping_faces")},
                        {"pins_valid",checked.report.at("pins_valid")},{"protected_cuts_valid",checked.report.at("protected_cuts_valid")},
                        {"cut_sectors_valid",checked.report.at("cut_sectors_valid")}});
                }
                return accepted;
            };
            auto trial=repair_detail::join_small_charts(phase,out.cut,stopped,valid,joining);
            if(trial) {
                // Accepted improvements survive optional timeout. Recheck the full
                // no-new-defect contract before adoption; final validity stays separate.
                auto candidate_cut=cut_surface(in.mesh,trial->seams);
                auto candidate_validation=validate(trial->uv,candidate_cut);
                const auto candidate_coverage=validate(trial->uv,candidate_cut,false,std::nullopt,true);
                joining["coverage_defects_after"]=candidate_coverage.defect_faces;
                if(no_new_defects(phase_coverage,candidate_coverage)) {
                    out.uv=std::move(trial->uv);out.seams=std::move(trial->seams);out.cut=std::move(candidate_cut);
                    validation=std::move(candidate_validation);out.changed_faces.clear();
                    for(MeshIndex f=0;f<orientation.size();++f)for(auto c:in.mesh->face_corners(f))if(out.uv[c]!=(*in.uv)[c]){out.changed_faces.push_back(f);break;}
                    joining["accepted"]=true;joining["remaining_invalid_charts"]=validation.invalid;
                    if(!phase_expired)joining["status"]="accepted";
                } else {
                    joining["status"]="rejected";joining["reason"]="final_validation";
                    joining["proposed_coverage_defects"]=candidate_coverage.defect_faces;joining["coverage_defects_after"]=phase_coverage.defect_faces;
                    joining["proposed_after"]=joining["after"];joining["after"]=joining["before"];joining["joined"]=0;
                    joining["proposed_anisotropy_p95"]=joining["anisotropy_p95_after"];joining["anisotropy_p95_after"]=joining["anisotropy_p95_before"];
                }
            }
            stopped();
        }
        if(progress)progress({{"stage","join_small_islands"},{"status",joining["status"]}});
    }
    out.invalid_charts=validation.invalid;
    out.exempt.resize(in.mesh->geometry().face_degenerate_kind.size());
    for(std::size_t f=0;f<out.exempt.size();++f)out.exempt[f]=in.mesh->geometry().face_degenerate_kind[f]!=0;
    out.accepted=validation.valid;
    const bool repaired=!out.changed_faces.empty() || joining["accepted"].get<bool>() || singleton["accepted"].get<bool>();
    out.reason=out.accepted?(repaired?"repaired":"unchanged"):"invalid_layout";
    timed_out=control.check() || deadline.expired() || timed_out;
    out.report={{"operation","postprocess_layout"},{"raw_generation",{{"status","skipped"}}},
        {"postprocessing",{{"status",timed_out?"timed_out":repaired?"repaired":"unchanged"},
            {"reason",timed_out?nlohmann::json("deadline"):in.options.repair.join_seconds>0?joining["reason"]:in.options.repair.singleton_seconds>0?singleton["reason"]:protected_recovery["accepted"].get<bool>() && !regions["accepted"].get<bool>()?protected_recovery["reason"]:in.options.repair.host_proposal_seconds>0?fragments["reason"]:in.options.repair.region_recovery_seconds>0?regions["reason"]:relax["reason"]}}},
        {"relax_faces",relax},{"recover_fragments",fragments},{"repair_regions",regions},{"recover_protected_cuts",protected_recovery},{"join_small_islands",joining},{"repair_singletons",singleton},{"stored_validation",validation.report},{"final",validation.report},{"chart_over_budget",nlohmann::json::array()},
        {"cpu_timings",nlohmann::json::object()},{"timed_out",timed_out},
        {"stages",nlohmann::json::array({{{"stage","postprocess_layout"},{"status",timed_out?"timed_out":"complete"}}})}};
    control.check();
    if(progress)progress({{"stage","postprocess_layout"},{"status","complete"}});
    return out;
}
PipelineResult postprocess_layout(const PipelineInput& in,const ProgressSink& progress,const Control& control) {
    return postprocess_layout(in,progress,control,LayoutPostprocessServices{});
}
}
