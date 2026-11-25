// SPDX-License-Identifier: GPL-3.0-or-later
// Source core/pipeline.py: _chart_corners, _unwrap_charts, _face_rings, local_untangle,
// planar_projection, _chart_valid_now, repair_chart, repair_all, final_chart_gate,
// _budget_unwrap and absorb_repair_fragments.
#include "cleanseams/pipeline.hpp"
#include "layout_postprocess_internal.hpp"
#include "cleanseams/cache.hpp"
#include "cleanseams/errors.hpp"
#include "cleanseams/detail/numpy_compat.hpp"
#include <algorithm>
#include <cmath>
#include <set>
namespace cleanseams {
namespace {
using Json=nlohmann::json;
constexpr double nan_value=std::numeric_limits<double>::quiet_NaN();
bool finite(const Vec2& p) {return std::isfinite(p[0]) && std::isfinite(p[1]);}
std::vector<TransferIndex> face_list(std::span<const std::uint8_t> mask) {
    std::vector<TransferIndex> out;for(std::size_t f=0;f<mask.size();++f)if(mask[f])out.push_back(static_cast<TransferIndex>(f));return out;
}
QualityPolicy exempt_policy(std::span<const std::uint8_t> exempt) {QualityPolicy p;p.exempt_faces=face_list(exempt);return p;}
// Source chart_validity over all charts (per-chart counts do not depend on others).
ChartValidity validity(const std::shared_ptr<const QualityGeometry>& geometry,std::span<const Vec2> uv,std::span<const std::uint8_t> exempt,
                       std::optional<std::span<const MeshIndex>> charts=std::nullopt) {
    return QualityEvaluation(geometry,CornerUvView{uv},exempt_policy(exempt),charts).validity();
}
// Source _chart_valid_now: (valid, {flipped, nonfinite, collapsed, crossings}).
std::pair<bool,Json> valid_now(const std::shared_ptr<const QualityGeometry>& geometry,MeshIndex chart,std::span<const Vec2> uv,
                               std::span<const std::uint8_t> exempt) {
    const std::array<MeshIndex,1> one{chart};const auto v=validity(geometry,uv,exempt,std::span<const MeshIndex>(one));
    return {v.valid[chart]!=0,{{"flipped",v.flipped[chart]},{"nonfinite",v.nonfinite[chart]},{"collapsed",v.collapsed[chart]},{"crossings",v.crossings[chart]}}};
}
std::string error_text(const std::exception& error) {
    if(const auto* contract=dynamic_cast<const ContractError*>(&error))return "ContractError: "+std::string(contract->what());
    return std::string("RuntimeError: ")+error.what();
}
// A rung or the gate's split records a solver failure and moves on (source `except
// Exception`), but a cancellation and a pin conflict end the run (Gate A F29). The source
// raises neither inside these calls: its Cancelled comes only from progress events, and
// its per-sector pin dict keeps the last pin; natively legacy_010 keeps the last pin too
// and preserve refuses conflicting sectors before the unwrap stage (F7).
void rethrow_control(const std::exception& error) {
    if(const auto* contract=dynamic_cast<const ContractError*>(&error))
        if(contract->code==ErrorCode::cancelled || contract->code==ErrorCode::constraint_conflict)throw *contract;
}
const CachedChartMap* cached_map(const RepairContext& ctx,const ChartKey& key) {
    if(!ctx.uv_cache)return nullptr;
    const auto found=ctx.uv_cache->find(key.digest());
    return found!=ctx.uv_cache->end() && found->second.key==key?&found->second:nullptr;
}
std::shared_ptr<const QualityGeometry> geometry_of(const std::shared_ptr<const CutSurface>& cs) {return quality_geometry(cs);}
}
std::shared_ptr<const QualityGeometry> RepairState::geometry() const {
    if(!geometry_ || &geometry_->cut()!=cs.get())geometry_=geometry_of(cs);
    return geometry_;
}
namespace repair_detail {
// Source _chart_corners: corners of the given charts in corner order.
std::vector<MeshIndex> chart_corners(const CutSurface& cs,std::span<const MeshIndex> charts) {
    const auto& ct=cs.tables();std::vector<std::uint8_t> wanted(ct.n_charts,0);for(auto c:charts)wanted.at(c)=1;
    std::vector<MeshIndex> out;const auto& corner_face=cs.mesh().tables().corner_face;
    for(std::size_t c=0;c<corner_face.size();++c)if(wanted[ct.chart_of_face[corner_face[c]]])out.push_back(static_cast<MeshIndex>(c));
    return out;
}
// Source _unwrap_charts: solve only `charts` and copy their corners into a copy of uv.
std::vector<Vec2> unwrap_charts(const RepairContext& ctx,const CutSurface& cs,std::span<const MeshIndex> charts,std::span<const Vec2> uv,
                                std::optional<std::uint64_t> iterations=std::nullopt,std::optional<std::span<const Vec2>> init=std::nullopt) {
    PinMap sub;const std::set<MeshIndex> wanted(charts.begin(),charts.end());
    if(ctx.pins)for(const auto& [c,p]:*ctx.pins)if(wanted.contains(cs.tables().chart_of_corner[c]))sub[c]=p;
    const auto out=ctx.backend->core()?ctx.backend->unwrap(cs.tables().cut,charts,sub,iterations,init)
                                      :ctx.backend->unwrap(cs.tables().cut,charts,sub,iterations);
    std::vector<Vec2> result(uv.begin(),uv.end());
    for(auto c:chart_corners(cs,charts))result[c]=out.at(c);
    return result;
}
// Source _face_rings: faces of `chart` within `rings` steps across uncut edges.
std::set<MeshIndex> face_rings(const CutSurface& cs,MeshIndex chart,std::span<const MeshIndex> faces,std::int64_t rings) {
    const auto& m=cs.mesh();const auto& t=m.tables();const auto& ct=cs.tables();
    std::set<MeshIndex> region(faces.begin(),faces.end());std::vector<MeshIndex> frontier(region.begin(),region.end());
    for(std::int64_t r=0;r<rings;++r) {
        std::vector<MeshIndex> next;
        for(auto f:frontier)for(auto c:m.face_corners(f)) {
            const auto e=t.corner_edge[c];if(ct.cut[e])continue;
            for(auto g:t.edge_faces[e])if(!region.contains(g) && ct.chart_of_face[g]==chart) {region.insert(g);next.push_back(g);}
        }
        frontier=std::move(next);
    }
    return region;
}
}
using namespace repair_detail;
std::vector<Vec2> planar_projection(const CutSurface& cs,std::span<const MeshIndex> charts,std::span<const Vec2> uv,double density) {
    const auto& m=cs.mesh();const auto& g=m.geometry();const auto& in=m.input();
    std::vector<Vec2> out(uv.begin(),uv.end());const double scale=std::sqrt(std::max(density,1e-12));
    auto norm=[](const Vec3& v) {return std::sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);};
    for(auto ch:charts) {
        const auto faces=cs.tables().chart_faces[ch];if(faces.empty())continue;
        Vec3 n{0,0,0};for(std::size_t d=0;d<3;++d)for(auto f:faces)n[d]+=g.face_normal[f][d]*g.face_area[f];
        const double ln=norm(n);
        if(ln>1e-300)for(auto& x:n)x/=ln;else n={0,0,1};
        Vec3 e1=std::abs(n[0])<0.9?Vec3{1,0,0}:Vec3{0,1,0};
        const double along=e1[0]*n[0]+e1[1]*n[1]+e1[2]*n[2];for(std::size_t d=0;d<3;++d)e1[d]-=along*n[d];
        const double l1=std::max(norm(e1),1e-300);for(auto& x:e1)x/=l1;
        const Vec3 e2{n[1]*e1[2]-n[2]*e1[1],n[2]*e1[0]-n[0]*e1[2],n[0]*e1[1]-n[1]*e1[0]};
        std::vector<MeshIndex> corners;for(auto f:faces)for(auto c:m.face_corners(f))corners.push_back(c);
        std::vector<Vec3> P;for(auto c:corners)P.push_back(in.positions[static_cast<std::size_t>(in.corner_vertices[c])]);
        Vec3 mean{0,0,0};for(const auto& p:P)for(std::size_t d=0;d<3;++d)mean[d]+=p[d];
        for(auto& x:mean)x/=static_cast<double>(P.size());
        Vec2 centre{0,0};bool all_finite=true;for(auto c:corners)all_finite=all_finite && finite(out[c]);
        if(all_finite) {for(auto c:corners)for(std::size_t d=0;d<2;++d)centre[d]+=out[c][d];for(auto& x:centre)x/=static_cast<double>(corners.size());}
        for(std::size_t i=0;i<corners.size();++i) {
            Vec3 p;for(std::size_t d=0;d<3;++d)p[d]=P[i][d]-mean[d];
            out[corners[i]]={(p[0]*e1[0]+p[1]*e1[1]+p[2]*e1[2])*scale+centre[0],(p[0]*e2[0]+p[1]*e2[1]+p[2]*e2[2])*scale+centre[1]};
        }
    }
    return out;
}
std::optional<std::vector<Vec2>> local_untangle(const RepairContext& ctx,const std::shared_ptr<const QualityGeometry>& geometry,MeshIndex chart,
                                                std::span<const Vec2> uv,std::span<const MeshIndex> faces,std::int64_t rings) {
    const auto& cs=geometry->cut();const auto& m=cs.mesh();const auto& ct=cs.tables();
    const auto region=face_rings(cs,chart,faces,rings);
    const auto chart_faces=ct.chart_faces[chart];
    std::set<MeshIndex> region_sectors,rim;
    for(auto f:region)for(auto c:m.face_corners(f))region_sectors.insert(ct.sector_of_corner[c]);
    for(auto f:chart_faces)if(!region.contains(f))for(auto c:m.face_corners(f))if(region_sectors.contains(ct.sector_of_corner[c]))rim.insert(ct.sector_of_corner[c]);
    std::set<MeshIndex> free;for(auto s:region_sectors)if(!rim.contains(s))free.insert(s);
    if(ctx.pins)for(const auto& [c,p]:*ctx.pins)if(ct.chart_of_corner[c]==chart)free.erase(ct.sector_of_corner[c]);
    if(free.empty())return std::nullopt;
    // Sector pins in first-encounter order, user pins overriding their sectors.
    std::vector<std::pair<MeshIndex,Vec2>> spins;std::map<MeshIndex,std::size_t> slot;
    for(auto f:chart_faces)for(auto c:m.face_corners(f)) {
        const auto s=ct.sector_of_corner[c];
        if(!free.contains(s) && !slot.contains(s)) {slot[s]=spins.size();spins.emplace_back(s,uv[c]);}
    }
    if(ctx.pins)for(const auto& [c,p]:*ctx.pins) {
        if(ct.chart_of_corner[c]!=chart)continue;
        const auto s=ct.sector_of_corner[c];
        if(slot.contains(s))spins[slot[s]].second=p;else {slot[s]=spins.size();spins.emplace_back(s,p);}
    }
    if(spins.empty())return std::nullopt;
    // Source: the core path solves mean-value harmonic weights (no backend given).
    ChartSolveRequest request;request.geometry=geometry;
    request.chart=chart;for(const auto& [s,p]:spins)request.sector_pins.push_back({s,p});
    request.pin_completion=ctx.pin_completion;request.validate=false;
    auto solved=solve_harmonic(request);
    std::vector<Vec2> out(uv.begin(),uv.end());
    for(std::size_t i=0;i<solved.corner_ids.size();++i)out[solved.corner_ids[i]]=solved.corner_uv[i];
    return out;
}
Json repair_chart(const RepairContext& ctx,RepairState& state,MeshIndex chart,double density) {
    const auto& m=*ctx.mesh;const auto& o=*ctx.options;const auto& g=m.geometry();
    auto cs=state.cs;const auto& ct=cs->tables();auto geometry=state.geometry();
    auto exempt=state.exempt;auto uv=state.uv;
    const auto key=cs->chart_key(chart).legacy_digest();
    Json rec{{"chart",chart},{"faces",ct.chart_faces[chart].size()},{"key",key},{"levels",Json::array()},{"level",""},{"status","unrepairable"},
             {"exempted",Json::array()}};
    auto [ok_before,before]=valid_now(geometry,chart,uv,exempt);(void)ok_before;
    rec["before"]=before;
    const auto* cached=cached_map(ctx,cs->chart_key(chart));
    auto finish=[&](std::vector<Vec2> new_uv,std::vector<std::uint8_t> new_exempt) {state.uv=std::move(new_uv);state.exempt=std::move(new_exempt);return rec;};
    // L1: re-solve alone with twice the iterations from the last validated map of this face set.
    std::optional<std::vector<Vec2>> init;
    if(cached) {init=std::vector<Vec2>(uv.size(),Vec2{nan_value,nan_value});for(std::size_t i=0;i<cached->corner_ids.size();++i)(*init)[cached->corner_ids[i]]=cached->corner_uv[i];}
    std::vector<Vec2> uv1;bool ok=false;Json after;
    try {
        const std::array<MeshIndex,1> one{chart};
        uv1=unwrap_charts(ctx,*cs,one,uv,2*ctx.iterations,init?std::optional<std::span<const Vec2>>(*init):std::nullopt);
        std::tie(ok,after)=valid_now(geometry,chart,uv1,exempt);
    } catch(const std::exception& error) {rethrow_control(error);ok=false;after={{"error",error_text(error)}};uv1=uv;}
    rec["levels"].push_back({{"level","L1"},{"ok",ok},{"after",after}});
    if(ok) {rec["level"]="L1";rec["status"]="repaired";return finish(std::move(uv1),std::move(exempt));}
    auto count=[](const Json& a,const char* k) {return a.contains(k)?a.at(k).get<double>():1e9;};
    if(count(after,"flipped")+count(after,"crossings")<before.at("flipped").get<double>()+before.at("crossings").get<double>())uv=uv1;
    // offenders
    const std::array<MeshIndex,1> one{chart};
    const auto off_ids=QualityEvaluation(geometry,CornerUvView{uv},exempt_policy(exempt),std::span<const MeshIndex>(one)).offending_faces(chart);
    const std::set<MeshIndex> off(off_ids.begin(),off_ids.end());
    const auto chart_faces=ct.chart_faces[chart];
    std::vector<double> areas;for(auto f:chart_faces)areas.push_back(g.face_area[f]);
    const double chart_area=detail::numpy_sum(areas);
    bool any_boundary=false,all_small=true;
    for(auto f:off_ids) {
        for(auto c:m.face_corners(f))any_boundary=any_boundary || ct.cut[m.tables().corner_edge[c]];
        all_small=all_small && g.face_area[f]<o.repair_exempt_fraction*chart_area;
    }
    // L2: exempt a few tiny interior offenders
    if(!off_ids.empty() && static_cast<std::int64_t>(off_ids.size())<=o.repair_exempt_faces && all_small && !any_boundary) {
        auto ex2=exempt;for(auto f:off_ids)ex2[f]=1;
        auto [ok2,after2]=valid_now(geometry,chart,uv,ex2);
        rec["levels"].push_back({{"level","L2"},{"ok",ok2},{"faces",off_ids},{"after",after2}});
        if(ok2) {rec["level"]="L2";rec["status"]="exempted";rec["exempted"]=off_ids;return finish(std::move(uv),std::move(ex2));}
    }
    // L2b: mean value weights on the 2-ring, then on the 3-ring
    if(!off_ids.empty())for(std::int64_t rings:{2,3}) {
        const auto label="L2b-"+std::to_string(rings);std::optional<std::vector<Vec2>> uv2;
        try {uv2=local_untangle(ctx,geometry,chart,uv,off_ids,rings);}
        catch(const std::exception& error) {rethrow_control(error);rec["levels"].push_back({{"level",label},{"ok",false},{"error",error_text(error)}});continue;}
        if(!uv2) {rec["levels"].push_back({{"level",label},{"ok",false},{"after","no_free_sectors"}});continue;}
        auto [ok2,after2]=valid_now(geometry,chart,*uv2,exempt);
        rec["levels"].push_back({{"level",label},{"ok",ok2},{"after",after2}});
        if(ok2) {rec["level"]="L2b";rec["status"]="repaired";return finish(std::move(*uv2),std::move(exempt));}
    }
    // L5: the last validated map of this face set
    if(cached) {
        auto uv5=uv;for(std::size_t i=0;i<cached->corner_ids.size();++i)uv5[cached->corner_ids[i]]=cached->corner_uv[i];
        if(ctx.pins)for(const auto& [c,p]:*ctx.pins)if(ct.chart_of_corner[c]==chart)uv5[c]=p;
        auto [ok5,after5]=valid_now(geometry,chart,uv5,exempt);
        rec["levels"].push_back({{"level","L5"},{"ok",ok5},{"after",after5}});
        if(ok5) {rec["level"]="L5";rec["status"]="repaired";return finish(std::move(uv5),std::move(exempt));}
    }
    // fallback: cut the offenders off
    if(off_ids.empty() || off_ids.size()==chart_faces.size()) {
        if(off_ids.size()==chart_faces.size() && chart_faces.size()==1) {
            for(auto f:off_ids)exempt[f]=1;
            auto uvp=planar_projection(*cs,one,uv,density*o.fallback_density);
            rec["levels"].push_back({{"level","fallback"},{"ok",true},{"faces",off_ids}});
            rec["level"]="fallback";rec["status"]="unresolved";rec["exempted"]=off_ids;
            return finish(std::move(uvp),std::move(exempt));
        }
        rec["levels"].push_back({{"level","fallback"},{"ok",false},{"after","no_offender"}});
        return finish(std::move(uv),std::move(exempt));
    }
    auto seams2=state.seams;
    for(auto f:off_ids)for(auto c:m.face_corners(f)) {const auto e=m.tables().corner_edge[c];if(!m.tables().edge_boundary[e])seams2[e]=1;}
    auto cs2=cut_surface(ctx.mesh,seams2);
    for(auto f:off_ids)exempt[f]=1;
    std::set<MeshIndex> cut_off,remainder;
    for(auto f:off_ids)cut_off.insert(cs2->tables().chart_of_face[f]);
    for(auto f:chart_faces)if(!off.contains(f))remainder.insert(cs2->tables().chart_of_face[f]);
    const std::vector<MeshIndex> cut_off_list(cut_off.begin(),cut_off.end()),remainder_list(remainder.begin(),remainder.end());
    auto uv3=planar_projection(*cs2,cut_off_list,uv,density*o.fallback_density);
    bool ok3=false;Json after3;
    try {
        uv3=unwrap_charts(ctx,*cs2,remainder_list,uv3,ctx.iterations);
        const auto v=validity(geometry_of(cs2),uv3,exempt);ok3=true;
        std::uint64_t f=0,n=0,c=0,x=0;
        for(auto ch:remainder_list) {ok3=ok3 && v.valid[ch];f+=v.flipped[ch];n+=v.nonfinite[ch];c+=v.collapsed[ch];x+=v.crossings[ch];}
        after3={{"flipped",f},{"nonfinite",n},{"collapsed",c},{"crossings",x}};
    } catch(const std::exception& error) {rethrow_control(error);ok3=false;after3={{"error",error_text(error)}};}
    rec["levels"].push_back({{"level","fallback"},{"ok",ok3},{"faces",off_ids},{"remainder_charts",remainder_list},{"after",after3}});
    rec["level"]="fallback";rec["status"]=ok3?"unresolved":"split";rec["exempted"]=off_ids;
    state.seams=std::move(seams2);state.cs=std::move(cs2);
    return finish(std::move(uv3),std::move(exempt));
}
RepairAllResult repair_all(const RepairContext& ctx,RepairState& state,const std::function<void(const std::string&)>& say,std::size_t max_repairs) {
    RepairAllResult out;std::set<std::string> tried;
    double density=QualitySnapshot(state.geometry(),CornerUvView{state.uv},exempt_policy(state.exempt)).report().density;
    if(density==0)density=1.0; // source: report.get("density", 1.0) or 1.0
    auto invalid_charts=[&]() {
        const auto v=validity(state.geometry(),state.uv,state.exempt);std::vector<MeshIndex> bad;
        for(MeshIndex ch=0;ch<state.cs->tables().n_charts;++ch)if(!v.valid[ch] && !state.cs->tables().chart_faces[ch].empty())bad.push_back(ch);
        return bad;
    };
    for(std::size_t i=0;i<max_repairs;++i) {
        std::vector<MeshIndex> bad;
        for(auto ch:invalid_charts())if(!tried.contains(state.cs->chart_key(ch).digest()))bad.push_back(ch);
        if(bad.empty())break;
        const auto ch=bad.front();tried.insert(state.cs->chart_key(ch).digest());
        if(say)say("repairing island "+std::to_string(ch)+" ("+std::to_string(state.cs->tables().chart_faces[ch].size())+" faces)");
        auto rec=repair_chart(ctx,state,ch,density);
        const bool stop=rec.at("status")=="unrepairable";out.records.push_back(std::move(rec));
        if(stop)break;
    }
    out.invalid=invalid_charts();
    return out;
}
Json area_over_budget(const RepairState& state,const PipelineOptions& o,std::optional<std::span<const MeshIndex>> charts) {
    const auto& ct=state.cs->tables();const double budget=o.layout.max_area_p95*o.final_gate_tolerance;
    std::vector<MeshIndex> all;
    if(!charts)for(MeshIndex ch=0;ch<ct.n_charts;++ch)all.push_back(ch);
    const QualitySnapshot snapshot(state.geometry(),CornerUvView{state.uv},exempt_policy(state.exempt));
    Json out=Json::array();
    for(auto ch:charts?*charts:std::span<const MeshIndex>(all)) {
        if(static_cast<std::int64_t>(ct.chart_faces[ch].size())<o.layout.split_min_faces)continue;
        const double a=snapshot.report(ch).area_p95;
        if(std::isfinite(a) && a>budget)out.push_back(Json::array({ch,detail::python_round(a,3)}));
    }
    return out;
}
Json area_only_rows(std::span<const MeshIndex> charts,std::span<const double> anisotropy_p95) {
    Json out=Json::array();
    for(auto ch:charts)out.push_back(Json::array({ch,std::isfinite(anisotropy_p95[ch])?Json(detail::python_round(anisotropy_p95[ch],3)):Json(nullptr)}));
    return out;
}
Json final_chart_gate(const RepairContext& ctx,RepairState& state,const std::function<void(const std::string&)>& say) {
    const auto& m=*ctx.mesh;const auto& o=*ctx.options;const auto& t=m.tables();
    const double budget=o.layout.max_anisotropy_p95*o.final_gate_tolerance;
    const auto p95=QualitySnapshot(state.geometry(),CornerUvView{state.uv},exempt_policy(state.exempt)).chart_anisotropy_p95();
    const auto before=state.cs; // keeps the gated surface alive after a split replaces it
    const auto& ct=before->tables();std::vector<MeshIndex> over;Json over_json=Json::array();
    for(MeshIndex ch=0;ch<ct.n_charts;++ch)
        if(std::isfinite(p95[ch]) && p95[ch]>budget && static_cast<std::int64_t>(ct.chart_faces[ch].size())>=o.layout.split_min_faces) {
            over.push_back(ch);over_json.push_back(Json::array({ch,detail::python_round(p95[ch],3)}));
        }
    Json info{{"budget",budget},{"over",over_json},{"split",Json::array()},{"chart_over_budget",Json::array()}};
    // FIX-04 (native_baseline): the area P95 constraint is judged too, with the same
    // tolerance, face floor and exemptions. Area-only violators are reported over budget;
    // the split pass remains the anisotropy one.
    std::vector<MeshIndex> area_only;
    const bool area_gate=o.policy!="legacy_010";
    if(area_gate) {
        info["area_budget"]=o.layout.max_area_p95*o.final_gate_tolerance;info["area_over"]=area_over_budget(state,o);
        for(const auto& row:info["area_over"])
            if(const auto ch=row[0].get<MeshIndex>();std::find(over.begin(),over.end(),ch)==over.end())area_only.push_back(ch);
    }
    // Every exit that keeps the gated surface lists the area-only charts after `over`.
    auto add_area_only=[&](Json& listed) {for(auto& row:area_only_rows(area_only,p95))listed.push_back(std::move(row));};
    if(over.empty()) {add_area_only(info["chart_over_budget"]);return info;}
    if(say)say(std::to_string(over.size())+" island(s) over budget: one split pass");
    const auto nf=m.input().face_offsets.size()-1,nv=m.input().positions.size();
    std::size_t quads=0;for(MeshIndex f=0;f<nf;++f)quads+=m.input().face_offsets[f+1]-m.input().face_offsets[f]==4;
    const auto regular=static_cast<double>(std::count(t.vert_regular.begin(),t.vert_regular.end(),1));
    const bool segmented=o.layout.seeding=="voronoi" || (o.layout.seeding=="auto" && static_cast<double>(quads)/static_cast<double>(nf)<.6 &&
                                                         regular/static_cast<double>(nv)<.3);
    DiskLayout split;
    if(segmented)split=split_regions(ctx.mesh,*state.cs,over);
    else {
        std::vector<std::uint8_t> sources(nv,0);
        for(auto ch:over)if(auto v=interior_seed(*state.cs,ch,sources,std::nullopt))sources[*v]=1;
        if(std::none_of(sources.begin(),sources.end(),[](auto x) {return x!=0;})) {
            info["chart_over_budget"]=info["over"];add_area_only(info["chart_over_budget"]);return info;
        }
        LayoutReport scratch;split=make_all_disks(ctx.mesh,state.cs->tables().cut,sources,scratch);
    }
    EdgeMask seams2(split.cut.size());for(std::size_t e=0;e<seams2.size();++e)seams2[e]=split.cut[e] && !t.edge_boundary[e];
    const std::set<MeshIndex> over_set(over.begin(),over.end());std::set<MeshIndex> fresh;
    for(MeshIndex f=0;f<nf;++f)if(over_set.contains(ct.chart_of_face[f]))fresh.insert(split.cs->tables().chart_of_face[f]);
    const std::vector<MeshIndex> new_charts(fresh.begin(),fresh.end());
    bool ok=false,improved=false;double worst_new=nan_value,worst_old=nan_value;std::vector<Vec2> uv2;std::vector<double> p2;
    try {
        uv2=unwrap_charts(ctx,*split.cs,new_charts,state.uv,ctx.iterations);
        const auto geometry=geometry_of(split.cs);
        const auto v2=validity(geometry,uv2,state.exempt);
        p2=QualitySnapshot(geometry,CornerUvView{uv2},exempt_policy(state.exempt)).chart_anisotropy_p95();
        ok=std::all_of(new_charts.begin(),new_charts.end(),[&](auto ch) {return v2.valid[ch]!=0;});
        for(auto ch:new_charts)if(!std::isnan(p2[ch]))worst_new=std::isnan(worst_new)?p2[ch]:std::max(worst_new,p2[ch]);
        worst_old=p95[over.front()];for(auto ch:over)worst_old=std::max(worst_old,p95[ch]);
        improved=std::isfinite(worst_new) && worst_new<worst_old-1e-6;
    } catch(const std::exception& error) {rethrow_control(error);ok=false;improved=false;worst_new=nan_value;worst_old=nan_value;info["error"]=error_text(error);}
    auto rounded=[](double x) {return std::isfinite(x)?Json(detail::python_round(x,3)):Json(nullptr);};
    info["split"]=Json::array({{{"charts_before",over},{"charts_after",new_charts.size()},{"valid",ok},{"worst_before",rounded(worst_old)},
                                {"worst_after",rounded(worst_new)}}});
    if(ok && improved) {
        Json still=Json::array();
        for(auto ch:new_charts)if(std::isfinite(p2[ch]) && p2[ch]>budget)still.push_back(Json::array({ch,detail::python_round(p2[ch],3)}));
        info["chart_over_budget"]=still;
        state.uv=std::move(uv2);state.seams=std::move(seams2);state.cs=std::move(split.cs);
        if(area_gate) {
            // Area-only charts were not split (their faces keep one chart of the new surface);
            // the pieces of the split are judged on area too.
            std::set<MeshIndex> extra;for(auto ch:area_only)extra.insert(state.cs->tables().chart_of_face[ct.chart_faces[ch][0]]);
            std::set<MeshIndex> listed;for(const auto& row:still)listed.insert(row[0].get<MeshIndex>());
            for(const auto& row:area_over_budget(state,o,std::span<const MeshIndex>(new_charts)))
                if(const auto ch=row[0].get<MeshIndex>();!listed.contains(ch))extra.insert(ch);
            const std::vector<MeshIndex> rows(extra.begin(),extra.end());
            for(auto& row:area_only_rows(rows,p2))info["chart_over_budget"].push_back(std::move(row));
        }
        return info;
    }
    info["chart_over_budget"]=info["over"];add_area_only(info["chart_over_budget"]);
    return info;
}
std::uint64_t budget_unwrap(const RepairContext& ctx,RepairState& state,const std::optional<std::vector<Vec2>>& cached,const Deadline& deadline,
                            Json& fallback_info) {
    const auto& m=*ctx.mesh;const auto& t=m.tables();const auto& g=m.geometry();const auto& in=m.input();
    const auto nc=t.corner_face.size(),nf=in.face_offsets.size()-1;
    auto uv=cached?*cached:std::vector<Vec2>(nc,Vec2{nan_value,nan_value});
    auto cs=state.cs;const auto& ct=cs->tables();auto geometry=state.geometry();
    std::set<MeshIndex> pinned_charts;if(ctx.pins)for(const auto& [c,p]:*ctx.pins)pinned_charts.insert(ct.chart_of_corner[c]);
    fallback_info={{"conformal_charts",0},{"conformal_faces",0},{"harmonic_charts",0},{"harmonic_faces",0}};
    const bool recovery=ctx.options && ctx.options->policy=="repair_v1";
    if(recovery) {
        pipeline_options_json(*ctx.options);
        if(ctx.recovery_locked && (ctx.recovery_locked->size()!=t.edge_vertices.size() ||
           std::any_of(ctx.recovery_locked->begin(),ctx.recovery_locked->end(),[](auto v){return v>1;})))
            throw ContractError(ErrorCode::invalid_input,"repair.locked","expected a boolean mask in original edge order");
    }
    double local_remaining=recovery?ctx.options->repair.region_recovery_seconds:0;
    if(!std::isfinite(local_remaining) || local_remaining<0)throw ContractError(ErrorCode::invalid_option,"repair.region_recovery_seconds","must be finite and nonnegative");
    local_remaining=std::min(1.0,local_remaining);std::size_t local_attempts=0;
    EdgeMask local_cuts(state.seams.size(),0);
    const auto clock=std::make_shared<const ClockSource>(ctx.recovery_clock?ctx.recovery_clock:ClockSource{[]{return MonotonicClock::now();}});
    if(recovery) {
        for(const auto* key:{"local_charts","local_faces","local_isolated_faces","local_attempts","deadline_expired_charts","rejected_charts"})fallback_info[key]=0;
        fallback_info["local_records"]=Json::array();
    }
    const std::vector<std::uint8_t> no_exempt(nf,0);
    auto valid=[&](MeshIndex ch,const std::vector<Vec2>& map) {return valid_now(geometry,ch,map,no_exempt).first;};
    // Charts whose cached map is finite, keeps its pins and validates are kept. A chart's
    // check reads only its own corners, which no other chart's solve or fallback writes.
    std::vector<MeshIndex> needs;
    for(MeshIndex ch=0;ch<ct.n_charts;++ch) {
        const std::array<MeshIndex,1> one{ch};
        bool kept=true;for(auto c:chart_corners(*cs,one))kept=kept && finite(uv[c]);
        if(kept && ctx.pins)for(const auto& [c,p]:*ctx.pins)if(ct.chart_of_corner[c]==ch && uv[c]!=p)kept=false;
        if(!(kept && valid(ch,uv)))needs.push_back(ch);
    }
    // Copy a solved chart's corners into the map when that chart validates.
    auto commit=[&](MeshIndex ch,const std::vector<Vec2>& trial) {
        if(!valid(ch,trial))return false;
        const std::array<MeshIndex,1> one{ch};for(auto c:chart_corners(*cs,one))uv[c]=trial[c];
        return true;
    };
    auto solve_one=[&](MeshIndex ch) {
        if(deadline.expired())return false;
        try {
            const std::array<MeshIndex,1> one{ch};
            return commit(ch,unwrap_charts(ctx,*cs,one,uv));
        } catch(const ContractError& error) {
            // As in the rungs (rethrow_control, Gate A F29): the source's `except Exception:
            // pass` here never saw a Cancelled (raised only from progress events), and a
            // preserved pin is never overwritten.
            if(error.code==ErrorCode::constraint_conflict || error.code==ErrorCode::cancelled)throw;
        } catch(const std::exception&) {}
        return false;
    };
    std::vector<MeshIndex> missing;
    auto fallback=[&](MeshIndex ch) {
        if(recovery) {
            auto& count=fallback_info[deadline.expired()?"deadline_expired_charts":"rejected_charts"];count=count.get<std::int64_t>()+1;
            throw_if_cancelled(ctx.cancel);
        }
        const auto faces=ct.chart_faces[ch];
        const bool triangles=std::all_of(faces.begin(),faces.end(),[&](auto f) {return in.face_offsets[f+1]-in.face_offsets[f]==3;});
        if(ct.n_charts>1 && !pinned_charts.contains(ch) && triangles) {
            // A single conformal solve keeps small regions' shape; very large unsolved
            // regions use the positive-weight map.
            std::vector<std::string> methods;if(faces.size()<=50000)methods.push_back("conformal");methods.push_back("harmonic");
            for(const auto& kind:methods) {
                std::optional<ChartSolveResult> mapped;
                ChartSolveRequest request;request.geometry=geometry;request.chart=ch;request.validate=false;request.pin_completion=ctx.pin_completion;
                try {mapped=kind=="conformal"?std::optional(solve_lscm(request)):convex_boundary_map(request);}
                catch(const ContractError& error) {
                    if(recovery && (error.code==ErrorCode::cancelled || error.code==ErrorCode::memory_limit || error.code==ErrorCode::constraint_conflict))throw;
                    mapped.reset();
                }
                if(!mapped)continue;
                auto trial=uv;for(std::size_t i=0;i<mapped->corner_ids.size();++i)trial[mapped->corner_ids[i]]=mapped->corner_uv[i];
                if(valid(ch,trial)) {
                    uv=std::move(trial);
                    fallback_info[kind+"_charts"]=fallback_info[kind+"_charts"].get<std::int64_t>()+1;
                    fallback_info[kind+"_faces"]=fallback_info[kind+"_faces"].get<std::int64_t>()+static_cast<std::int64_t>(faces.size());
                    return;
                }
            }
        }
        if(recovery && !triangles && !pinned_charts.contains(ch) && local_attempts<8 && local_remaining>0 && !deadline.expired()) {
            ++local_attempts;fallback_info["local_attempts"]=local_attempts;
            const auto started=(*clock)();
            auto until=started+std::chrono::duration_cast<MonotonicClock::duration>(std::chrono::duration<double>(local_remaining));
            // Each deadline keeps its own clock; the operation is checked separately by
            // recovery as well, including injected check-count deadlines.
            const Deadline local(until,clock);
            std::optional<LocalRecovery> mapped;Json diagnostic;
            try {mapped=recover_mixed_chart(ctx,state,ch,deadline,local,diagnostic);}
            catch(...) {local_remaining-=std::max(0.0,std::chrono::duration<double>((*clock)()-started).count());throw;}
            local_remaining-=std::max(0.0,std::chrono::duration<double>((*clock)()-started).count());
            fallback_info["local_records"].push_back(diagnostic);
            if(mapped) {
                for(std::size_t i=0;i<mapped->corners.size();++i)uv[mapped->corners[i]]=mapped->values[i];
                for(auto f:mapped->isolated_faces)for(auto c:m.face_corners(f))local_cuts[t.corner_edge[c]]=1;
                fallback_info["local_charts"]=fallback_info["local_charts"].get<std::size_t>()+1;
                fallback_info["local_faces"]=fallback_info["local_faces"].get<std::size_t>()+faces.size();
                fallback_info["local_isolated_faces"]=fallback_info["local_isolated_faces"].get<std::size_t>()+mapped->isolated_faces.size();
                fallback_info["local_records"].back().update({{"chart",ch},{"isolated_faces",mapped->isolated_faces},
                    {"anisotropy_p95",mapped->anisotropy_p95},{"stretch_p95",mapped->anisotropy_p95-1.0},
                    {"weighting","original_face_area_step_cdf"},{"triangulation","mesh_snapshot_fixed"},{"storage_precision",ctx.options->repair.storage_precision}});
                return;
            }
        }
        missing.push_back(ch);
    };
    // T24: the source solved one chart per call until the deadline passed. Bounded batches
    // of compatible charts (one surface, solver, iteration count and start policy) are
    // solved together; the deadline is checked before each batch and cancellation before
    // and after it. Charts of a batch the deadline stops, and charts whose map does not
    // validate, take the fallbacks right away, as the source did chart by chart.
    const auto tri=chart_triangle_counts(*cs);
    const auto batches=plan_batches(needs,tri,ctx.scheduling);
    const auto call=ctx.batches && !batches.empty()?ctx.batches->begin_call():0;
    auto record=[&](const std::vector<MeshIndex>& batch,bool evaluated,const char* stop) {
        if(!ctx.batches)return;
        BatchRecord r;r.stage="unwrap";r.call=call;r.charts=batch;r.evaluated=evaluated;r.stop=stop;
        for(auto ch:batch)r.triangles+=tri[ch];
        ctx.batches->add(std::move(r));
    };
    for(std::size_t b=0;b<batches.size();++b) {
        const auto& batch=batches[b];
        if(ctx.cancel.requested()) {for(auto k=b;k<batches.size();++k)record(batches[k],false,"cancelled");throw_if_cancelled(ctx.cancel);}
        if(deadline.expired()) {
            record(batch,false,"deadline");
            for(auto ch:batch)fallback(ch);
            continue;
        }
        std::optional<std::vector<Vec2>> trial;
        try {trial=unwrap_charts(ctx,*cs,batch,uv);}
        catch(const ContractError& error) {
            if(error.code==ErrorCode::constraint_conflict || error.code==ErrorCode::cancelled)throw;
        } catch(const std::exception&) {}
        record(batch,true,"");
        throw_if_cancelled(ctx.cancel); // after the active solve: nothing is committed
        for(auto ch:batch) {
            // A failed batch is retried chart by chart, so a failure stays with its chart.
            const bool solved=trial?commit(ch,*trial):solve_one(ch);
            if(!solved)fallback(ch);
        }
    }
    if(missing.empty()) {
        state.uv=std::move(uv);
        if(std::any_of(local_cuts.begin(),local_cuts.end(),[](auto x){return x!=0;})) {
            for(std::size_t e=0;e<state.seams.size();++e)state.seams[e]=(state.seams[e] || local_cuts[e]) && !t.edge_boundary[e];
            state.cs=cut_surface(ctx.mesh,state.seams);
        }
        return 0;
    }
    // Emergency face cuts only add seams and keep user pins; they are reported.
    const std::set<MeshIndex> missing_set(missing.begin(),missing.end());
    std::vector<std::uint8_t> face_missing(nf,0);std::uint64_t projected=0;
    for(MeshIndex f=0;f<nf;++f)if(missing_set.contains(ct.chart_of_face[f])) {face_missing[f]=1;++projected;}
    auto seams=state.seams;
    for(std::size_t e=0;e<seams.size();++e)seams[e]=seams[e] || local_cuts[e];
    std::vector<MeshIndex> corners;for(std::size_t c=0;c<nc;++c)if(face_missing[t.corner_face[c]]) {
        const auto e=t.corner_edge[c];
        if(recovery && ctx.recovery_locked && (*ctx.recovery_locked)[e] && !state.seams[e] && !t.edge_boundary[e])
            throw ContractError(ErrorCode::constraint_conflict,"repair.locked","emergency isolation would cut a locked continuous edge");
        corners.push_back(static_cast<MeshIndex>(c));seams[e]=1;
    }
    for(std::size_t e=0;e<seams.size();++e)if(t.edge_boundary[e])seams[e]=0;
    auto new_cs=cut_surface(ctx.mesh,seams);
    for(auto c:corners) {
        const auto f=t.corner_face[c];auto n=g.face_normal[f];
        if(std::sqrt(n[0]*n[0]+n[1]*n[1]+n[2]*n[2])<1e-12)n={0,0,1};
        std::size_t axis=0;for(std::size_t d=1;d<3;++d)if(std::abs(n[d])<std::abs(n[axis]))axis=d;
        Vec3 a{0,0,0};a[axis]=1;
        Vec3 u{a[1]*n[2]-a[2]*n[1],a[2]*n[0]-a[0]*n[2],a[0]*n[1]-a[1]*n[0]};
        const double lu=std::max(std::sqrt(u[0]*u[0]+u[1]*u[1]+u[2]*u[2]),1e-300);for(auto& x:u)x/=lu;
        const Vec3 v{n[1]*u[2]-n[2]*u[1],n[2]*u[0]-n[0]*u[2],n[0]*u[1]-n[1]*u[0]};
        const auto& P=in.positions[static_cast<std::size_t>(in.corner_vertices[c])];const auto& C=g.face_center[f];
        const std::array<double,3> p{P[0]-C[0],P[1]-C[1],P[2]-C[2]};
        uv[c]={detail::numpy_einsum_dot3(p,u),detail::numpy_einsum_dot3(p,v)};
    }
    if(ctx.pins)for(const auto& [c,p]:*ctx.pins)uv[c]=p;
    state.uv=std::move(uv);state.seams=std::move(seams);state.cs=std::move(new_cs);
    return projected;
}
Json absorb_repair_fragments(const RepairContext& ctx,RepairState& state,const ChartEvaluator* evaluate,std::span<const std::uint8_t> locked,
                             const Deadline& deadline) {
    const auto& m=*ctx.mesh;const auto& o=*ctx.options;
    const auto info_tiny=tiny_charts(m,*state.cs,o.layout);
    if(std::none_of(info_tiny.tiny.begin(),info_tiny.tiny.end(),[](auto x) {return x!=0;}) || deadline.expired())return {{"absorbed",0}};
    const auto& ct=state.cs->tables();
    EdgeMask protected_edges(locked.begin(),locked.end());
    std::set<MeshIndex> pinned_charts;if(ctx.pins)for(const auto& [c,p]:*ctx.pins)pinned_charts.insert(ct.chart_of_corner[c]);
    if(!pinned_charts.empty())
        for(auto c:ct.boundary_corners)if(pinned_charts.contains(ct.chart_of_corner[c]))protected_edges[m.tables().corner_edge[c]]=1;
    LayoutReport report;for(std::size_t f=0;f<state.exempt.size();++f)if(state.exempt[f])report.exempt_faces.push_back(static_cast<MeshIndex>(f));
    AbsorbRequest request;request.exempt=state.exempt;request.locked=protected_edges;request.deadline=deadline;
    auto trial=absorb_tiny(ctx.mesh,state.cs,o.layout,report,evaluate,request).cs;
    Json info{{"absorbed",report.absorbed},{"absorbed_relaxed",report.absorbed_relaxed},{"absorbed_exempt",report.absorbed_exempt},
              {"tiny_remaining",report.tiny_remaining},{"absorbed_faces",report.absorbed_faces}};
    if(!report.absorbed)return info;
    auto candidate=state.uv;
    for(MeshIndex ch=0;ch<trial->tables().n_charts;++ch)
        if(const auto* hit=report.cached(trial->chart_key(ch)))for(std::size_t i=0;i<hit->corner_ids.size();++i)candidate[hit->corner_ids[i]]=hit->corner_uv[i];
    auto new_exempt=state.exempt;for(auto f:report.exempt_faces)new_exempt[f]=1;
    const auto metrics=pipeline_metrics(trial,candidate,new_exempt);
    // A masked face cannot conceal a crossing of the actual visible outline.
    const auto crossings=boundary_crossings(*trial,CornerUvView{candidate}).counts;
    bool pins_moved=false;if(ctx.pins)for(const auto& [c,p]:*ctx.pins)pins_moved=pins_moved || candidate[c]!=p;
    if(!metrics.at("valid").get<bool>() || std::any_of(crossings.begin(),crossings.end(),[](auto x) {return x!=0;}) || pins_moved) {
        info["absorbed"]=0;info["reverted"]="final_validation";return info;
    }
    state.cs=std::move(trial);state.uv=std::move(candidate);state.exempt=std::move(new_exempt);
    state.seams.assign(state.cs->tables().cut.begin(),state.cs->tables().cut.end());
    for(std::size_t e=0;e<state.seams.size();++e)if(m.tables().edge_boundary[e])state.seams[e]=0;
    return info;
}
}
