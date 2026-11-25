// SPDX-License-Identifier: GPL-3.0-or-later
// Phase B physical-distance policy (Section 20, T39): the scale reference frozen at the
// operation start and the policy that turns the enhanced operations' corridors, smoothing and
// region floors into object-space units. Used by policy enhanced_v1 only.
#include "cleanseams/physical_scale.hpp"
#include "cleanseams/metrics.hpp"
#include "cleanseams/errors.hpp"
#include "cleanseams/detail/numpy_compat.hpp"
#include <algorithm>
#include <cmath>
#include <map>
#include <set>
namespace cleanseams {
namespace {
using Json=nlohmann::json;
void require(bool yes,const char* field,const char* message) {if(!yes)throw ContractError(ErrorCode::invalid_input,field,message);}
double median_or(std::vector<double> values,double fallback) {
    if(values.empty())return fallback;
    const double m=detail::numpy_median(std::move(values));
    return std::isfinite(m) && m>0?m:fallback;
}
Json finite_or_null(double x) {return std::isfinite(x)?Json(x):Json(nullptr);}
}
nlohmann::json ScaleReference::to_json() const {
    // Bounded: the first 64 reference charts in full, all of them counted.
    Json charts=Json::array();
    for(std::size_t k=0;k<chart_area.size() && k<64;++k)
        charts.push_back({{"chart",k},{"area",chart_area[k]},{"median_edge_length",chart_edge_length[k]},{"median_face_area",chart_face_area[k]}});
    return {{"geometry_key",geometry_key},{"revision",revision},{"corridor_rings",corridor_rings},{"area",area},
            {"median_edge_length",median_edge_length},{"median_face_area",median_face_area},{"reference_charts",chart_area.size()},{"charts",charts}};
}
ScaleReference freeze_scale_reference(const CutSurface& origin,std::span<const std::uint8_t> exempt,const PipelineOptions& options,Revision revision) {
    const auto& m=origin.mesh();const auto& t=m.tables();const auto& g=m.geometry();const auto& ct=origin.tables();
    const auto nf=g.face_area.size(),ne=t.edge_vertices.size();
    require(exempt.empty() || exempt.size()==nf,"enhanced.scale.exempt","exempt mask needs one flag per face");
    auto relevant=[&](std::size_t f) {return exempt.empty() || !exempt[f];};
    ScaleReference r;r.geometry_key=m.geometry_key();r.revision=revision;r.corridor_rings=options.paths.corridor_rings;
    std::vector<double> edges,faces;
    for(std::size_t e=0;e<ne;++e)if(!t.edge_faces[e].empty() && std::isfinite(g.edge_length[e]) && g.edge_length[e]>0)edges.push_back(g.edge_length[e]);
    for(std::size_t f=0;f<nf;++f)if(relevant(f)) {r.area+=g.face_area[f];if(g.face_area[f]>0)faces.push_back(g.face_area[f]);}
    r.median_edge_length=median_or(std::move(edges),0.0);r.median_face_area=median_or(std::move(faces),0.0);
    r.chart_of_face=ct.chart_of_face;r.face_area=g.face_area;
    const auto nc=static_cast<std::size_t>(ct.n_charts);
    r.chart_area.assign(nc,0.0);r.chart_edge_length.assign(nc,r.median_edge_length);r.chart_face_area.assign(nc,r.median_face_area);
    for(std::size_t ch=0;ch<nc;++ch) {
        // The chart's own resolution: the median length of the edges of its faces (the
        // reference T31 captured for its corridor) and the median area of its faces.
        std::set<MeshIndex> chart_edges;std::vector<double> areas;
        for(auto f:ct.chart_faces[ch]) {
            for(auto c:m.face_corners(f))chart_edges.insert(t.corner_edge[c]);
            if(relevant(f)) {r.chart_area[ch]+=g.face_area[f];if(g.face_area[f]>0)areas.push_back(g.face_area[f]);}
        }
        std::vector<double> lengths;for(auto e:chart_edges)lengths.push_back(g.edge_length[e]);
        r.chart_edge_length[ch]=median_or(std::move(lengths),r.median_edge_length);
        r.chart_face_area[ch]=median_or(std::move(areas),r.median_face_area);
    }
    return r;
}
nlohmann::json ScaleLength::to_json() const {return {{"value",value},{"source",source},{"reference_chart",reference_chart}};}
nlohmann::json SmoothingSteps::to_json() const {
    return {{"default_steps",default_steps},{"steps",steps},{"radius",finite_or_null(radius)},{"local_length",local_length},{"unit",unit},{"bounded",bounded}};
}
bool RegionFloor::below(double region_area,std::uint64_t region_faces) const {
    if(!physical)return static_cast<std::int64_t>(region_faces)<faces;
    return region_area<area*(1.0-1e-9);
}
nlohmann::json RegionFloor::to_json() const {
    return {{"faces",faces},{"area",physical?Json(area):Json(nullptr)},{"measure",physical?"area":"faces"},{"source",physical?source:"hops"}};
}
PhysicalScale::PhysicalScale(ScaleReference reference,ScaleOptions options):reference_(std::move(reference)),options_(std::move(options)) {
    require(options_.mode=="physical" || options_.mode=="hops","enhanced.scale.mode","mode must be physical or hops");
    for(const auto x:{options_.reference_length,options_.reference_area,options_.corridor_radius,options_.corridor_scale})
        require(std::isfinite(x) && x>=0,"enhanced.scale","explicit lengths must be finite and nonnegative");
    require(options_.max_steps>=1,"enhanced.scale.max_steps","must be at least 1");
}
PhysicalScale PhysicalScale::freeze(const CutSurface& origin,std::span<const std::uint8_t> exempt,const PipelineOptions& options,Revision revision) {
    return PhysicalScale(freeze_scale_reference(origin,exempt,options,revision),options.enhanced.scale);
}
MeshIndex PhysicalScale::reference_chart(std::span<const MeshIndex> faces) const {
    require(!faces.empty(),"enhanced.scale.faces","a reference chart needs faces");
    require(!reference_.chart_area.empty(),"enhanced.scale.reference","the scale reference has no charts");
    // Area-weighted majority: the reference chart holding the largest share of the faces' 3D area.
    std::map<MeshIndex,std::pair<double,std::uint64_t>> share;
    for(auto f:faces) {
        require(f<reference_.chart_of_face.size(),"enhanced.scale.faces","face outside the reference mesh");
        auto& s=share[reference_.chart_of_face[f]];s.first+=reference_.face_area[f];++s.second;
    }
    return std::max_element(share.begin(),share.end(),[](const auto& a,const auto& b) {
        if(a.second.first!=b.second.first)return a.second.first<b.second.first;
        if(a.second.second!=b.second.second)return a.second.second<b.second.second;
        return a.first>b.first;
    })->first;
}
ScaleLength PhysicalScale::local_length(std::span<const MeshIndex> faces) const {
    const auto ch=reference_chart(faces);
    return {reference_.chart_edge_length[ch],"reference",ch};
}
ScaleLength PhysicalScale::unit(std::span<const MeshIndex> faces) const {
    auto local=local_length(faces);
    if(physical() && options_.reference_length>0) {local.value=options_.reference_length;local.source="explicit";}
    return local;
}
ScaleLength PhysicalScale::corridor_radius(std::span<const MeshIndex> faces) const {
    const double rings=static_cast<double>(std::max<std::int64_t>(1,reference_.corridor_rings));
    const auto ch=reference_chart(faces);
    if(!physical())return {rings,"hops",ch};
    if(options_.corridor_radius>0)return {options_.corridor_radius,"explicit",ch};
    const auto u=unit(faces);
    return {(options_.corridor_scale>0?options_.corridor_scale:rings)*u.value,u.source=="explicit"?"reference_length":"reference",ch};
}
namespace {
std::int64_t scaled_steps(double steps,std::int64_t bound,bool& bounded) {
    if(!std::isfinite(steps) || steps>static_cast<double>(bound)) {bounded=true;return bound;}
    return static_cast<std::int64_t>(std::llround(std::max(0.0,steps)));
}
}
SmoothingSteps PhysicalScale::diffusion(std::int64_t default_steps,std::span<const MeshIndex> faces) const {
    require(default_steps>=0,"enhanced.scale.steps","steps must be nonnegative");
    SmoothingSteps s;s.default_steps=default_steps;s.steps=default_steps;
    s.local_length=local_length(faces).value;s.unit=unit(faces).value;
    // A diffusion smoother spreads like sqrt(steps) edge lengths: its default count, defined at
    // one step per unit, takes (unit / local)^2 times as many steps at the local resolution.
    if(physical() && s.unit!=s.local_length && s.local_length>0) {
        const double ratio=s.unit/s.local_length;
        s.steps=scaled_steps(static_cast<double>(default_steps)*ratio*ratio,options_.max_steps,s.bounded);
    }
    s.radius=std::sqrt(static_cast<double>(s.steps))*s.local_length;
    return s;
}
SmoothingSteps PhysicalScale::rings(std::int64_t default_passes,std::span<const MeshIndex> faces) const {
    require(default_passes>=0,"enhanced.scale.steps","passes must be nonnegative");
    SmoothingSteps s;s.default_steps=default_passes;s.steps=default_passes;
    s.local_length=local_length(faces).value;s.unit=unit(faces).value;
    // A ring operator grows by one face ring (about one edge length) per pass.
    if(physical() && s.unit!=s.local_length && s.local_length>0)
        s.steps=scaled_steps(static_cast<double>(default_passes)*s.unit/s.local_length,options_.max_steps,s.bounded);
    s.radius=static_cast<double>(s.steps)*s.local_length;
    return s;
}
RegionFloor PhysicalScale::floor(std::int64_t faces,std::span<const MeshIndex> region) const {
    RegionFloor out;out.faces=faces;out.physical=physical();
    if(!out.physical)return out;
    const auto count=static_cast<double>(std::max<std::int64_t>(0,faces));
    if(options_.reference_area>0) {out.area=count*options_.reference_area;out.source="explicit";return out;}
    const auto ch=reference_chart(region);
    const double local=reference_.chart_edge_length[ch],u=unit(region).value;
    const double ratio=local>0?u/local:1.0;
    // Area-weighted region preservation: the legacy count times the reference chart's median
    // face area, at the unit's resolution (a linear subdivision quarters both, so the floor is
    // the same object-space area whenever the same explicit reference length is supplied;
    // reference_area fixes it for tessellations of other shapes).
    out.area=count*reference_.chart_face_area[ch]*ratio*ratio;out.source=ratio!=1.0?"reference_length":"reference";
    return out;
}
double PhysicalScale::length_quantum(std::span<const MeshIndex> faces) const {
    if(!physical())return 0.0;
    const double u=unit(faces).value;
    if(!(u>0) || !std::isfinite(u))return 0.0;
    return std::ldexp(1.0,std::ilogb(u)-32);
}
nlohmann::json PhysicalScale::to_json() const {
    return {{"version",physical_scale_version},{"mode",options_.mode},
            {"options",{{"reference_length",options_.reference_length},{"reference_area",options_.reference_area},{"corridor_radius",options_.corridor_radius},
                        {"corridor_scale",options_.corridor_scale},{"max_steps",options_.max_steps}}},
            {"distance","surface (geodesic graph) distance along patch edges; never Euclidean"},
            {"decision_quantum","lengths on a 2^(floor(log2 unit) - 32) grid, dimensionless scores on a 2^-32 grid (physical mode)"},
            {"detectors","component- or chart-size-normalized (relief separation, cavity reach, shell extents)"},
            {"reference",reference_.to_json()}};
}
const PhysicalScale& context_scale(ProposalContext& ctx,const RunState& state) {
    require(bool(state.surface),"enhanced.scale.state","the state has no surface");
    if(!ctx.scale)ctx.scale=std::make_shared<const PhysicalScale>(PhysicalScale::freeze(*state.surface,ctx.exempt,ctx.options,state.revision));
    if(ctx.scale->reference().geometry_key!=state.surface->mesh().geometry_key() ||
       ctx.scale->reference().chart_of_face.size()!=state.surface->tables().chart_of_face.size())
        throw ContractError(ErrorCode::stale_state,"enhanced.scale","the frozen scale belongs to another geometry");
    return *ctx.scale;
}
nlohmann::json physical_scale_report(const ProposalContext& ctx,const RunState& origin,const RunState& retained) {
    Json out=ctx.scale?ctx.scale->to_json():Json{{"version",physical_scale_version},{"mode",ctx.options.enhanced.scale.mode},{"frozen",false}};
    if(origin.surface)out["candidate_zero"]=seam_geometry(*origin.surface,ctx.exempt).summary();
    if(retained.surface)out["retained"]=seam_geometry(*retained.surface,ctx.exempt).summary();
    return out;
}
}
