// SPDX-License-Identifier: GPL-3.0-or-later
#include "cleanseams/options.hpp"
#include "cleanseams/detail/option_fields.hpp"
#include "cleanseams/detail/codec.hpp"
#include <limits>
#include <map>
#include <set>
#include <string_view>
#include <utility>
#include <cmath>
namespace cleanseams {
namespace {
using detail::Json;
// GB M1: the largest strict-validity tolerance (enhanced.pin_tolerance, enhanced.overlap_area_tolerance).
constexpr double max_strict_tolerance=1e-6;
// GB m8: the largest per-chart candidate and solve bounds (the pass's and each feature's share;
// the UI metadata of the two panel controls carries the same maxima).
constexpr std::int64_t max_candidates_per_chart=1024;
constexpr std::int64_t max_solves_per_chart=4096;
Json& property(Json& schema,const std::string& path) {
    auto* result=&schema;std::size_t start=0;
    while(true) {const auto end=path.find('.',start);result=&(*result)["properties"].at(path.substr(start,end-start));if(end==std::string::npos) return *result;start=end+1;}
}
// Gate A F11: the structural checks (types, finite numbers, int64 range) hold for every
// field. Value ranges are kept only where the source itself refuses the value on every
// path it is used, or where the source's result is ill-defined; everywhere else any finite
// value runs exactly as in the source (threshold comparisons, loops that run zero times,
// fallbacks for unknown names). The kept ranges are the INPUT-02 refusals (docs/CONTRACTS.md).
Json pipeline_schema() {
    auto schema=detail::schema_for(PipelineOptions{});
    auto& required=schema["required"];required.erase(std::remove(required.begin(),required.end(),"repair"),required.end());
    auto choice=[&](const char* path,Json values){property(schema,path)["enum"]=std::move(values);};
    // The source raises for any other profile at the start of every run (pipeline.py:801);
    // policy, version, solver.* and host.* are native fields.
    choice("profile",{"auto","custom"});choice("policy",{"legacy_010","native_baseline","enhanced_v1","repair_v1"});
    choice("repair.storage_precision",{"float32","float64"});
    choice("repair.join_mode",{"small","fewer"});
    property(schema,"repair.minimum_relative_area")["minimum"]=0;
    for(const auto* name:{"repair.relax_faces_seconds","repair.region_recovery_seconds","repair.singleton_seconds","repair.join_seconds","repair.host_proposal_seconds"})property(schema,name)["minimum"]=0;
    choice("solver.method",{"arap","lscm","harmonic","blender"});choice("solver.backend",{"superlu","numpy_compat"});
    choice("host.pack_after",{"none","blender"});choice("version",{1});
    // Nonnegative where a negative value is refused or ill-defined in the source: a negative
    // deadline raises (control.py:32); negative routing weights break the Dijkstra search
    // (a negative cycle never ends); a negative feature scale raises in face_rings whenever
    // features run; a negative displacement limit divides by zero. Host iterations and the
    // packing sizes are native/packer fields with the packer's own limits (packing.py:42-44).
    // A negative exposure_K or layout.arap_iterations acts as 0 in the source; they are
    // refused rather than given a second spelling of 0 (0 itself is accepted).
    for(const auto* path:{"deadline_s","layout.feature_scale","layout.exposure_K","layout.arap_iterations","boundary.max_displacement","paths.w_turn","paths.w_flow","paths.w_drift",
        "paths.w_convex","paths.w_exposure","solver.iterations","packing.padding_px","packing.border_px","packing.time_budget"})
        property(schema,path)["minimum"]=0;
    property(schema,"paths.variants")["items"]["items"]["minimum"]=0;
    // Empty variants make the source's refine_cycle index an empty list.
    property(schema,"paths.variants")["minItems"]=1;
    // island_band refuses these at the start of a run (target.py:17-19).
    property(schema,"target_islands")["items"]["minimum"]=1;
    for(const auto* path:{"packing.texture_width","packing.texture_height","packing.trials","packing.threads"}) property(schema,path)["minimum"]=1;
    property(schema,"packing.tiles")["minItems"]=1;property(schema,"packing.tiles")["items"]["minimum"]=1001;
    property(schema,"packing.tiles")["uniqueItems"]=true;
    property(schema,"layout.max_defect")["description"]="null means unbounded legacy angle-defect proxy; no nonfinite number is accepted";
    property(schema,"paths.corridor_rings")["description"]="any integer under legacy_010 (the source searches from it); at least 1 under other policies (FIX-02)";
    property(schema,"paths.corridor_max_rings")["description"]="any integer under legacy_010; at least 1 under other policies (FIX-02)";
    property(schema,"paths.variants")["description"]="legacy_010 only: the fixed weight variants the source searches; other policies derive them from w_turn, w_flow and w_drift (FIX-01)";
    property(schema,"layout.solver")["description"]="lscm, harmonic, or anything else for ARAP (source flatten_chart)";
    property(schema,"layout.priority")["description"]="compact, small_first, or anything else for longest-arc order (source merge_patches)";
    property(schema,"boundary.mode")["description"]="auto, panels, or anything else to smooth every arc (source clean_chart)";
    property(schema,"layout.seeding")["description"]="auto, motorcycle or voronoi; any other value is refused when seams are generated (source generate_seams)";
    property(schema,"solver.method")["description"]="recorded, not a selector: the core backend always runs ARAP (source CoreBackend) and a host session the host's unwrap";
    property(schema,"solver.backend")["description"]="recorded, not a selector: the pipeline's sparse solves use the SciPy-compatible numerical policy";
    property(schema,"solver.iterations")["description"]="iterations a host session asks the host's unwrap for; the core backend uses 10 (source CoreBackend)";
    property(schema,"host")["description"]="Blender settings carried with the options for the host; the engine and the command line do not read them";
    // Phase B enhanced policy record (T28): read by enhanced_v1 only. Weights and tolerances
    // are nonnegative, bounds positive; the objective's version is a native constant.
    for(const auto* path:{"enhanced.objective.w_length","enhanced.objective.w_exposure","enhanced.objective.w_turn","enhanced.objective.w_fragment",
        "enhanced.objective.chart_count_weight","enhanced.objective.tie_band","enhanced.pin_tolerance","enhanced.overlap_area_tolerance"})
        property(schema,path)["minimum"]=0;
    property(schema,"enhanced.objective.tie_band")["maximum"]=0.5;
    // GB M1: no value of the two strict-validity tolerances may switch a strict check off. An
    // intersection never exceeds the smaller triangle, so an overlap tolerance of 1 or more
    // accepts every overlap; a pin tolerance of order 1 accepts a moved pin. Both are capped at
    // max_strict_tolerance (1e-6): overlaps above a millionth of the smaller triangle's UV area
    // stay overlaps, and a pin may move at most 1e-6 * (1 + |pin|) per coordinate, which still
    // admits float32 storage of a host map (half an ulp is at most 2^-24 * |pin|).
    for(const auto* path:{"enhanced.pin_tolerance","enhanced.overlap_area_tolerance"}) property(schema,path)["maximum"]=max_strict_tolerance;
    for(const auto* path:{"enhanced.preferences.w_exposure","enhanced.preferences.w_view","enhanced.preferences.w_avoid",
        "enhanced.preferences.w_allow","enhanced.preferences.w_crease"})
        property(schema,path)["minimum"]=0;
    property(schema,"enhanced.preferences")["description"]="enhanced_v1 painted preferences (T29): routing heuristic weights and the host "
        "attribute names read for required/forbidden edges and allow/avoid points (documented names or explicit aliases; empty disables one)";
    for(const auto* path:{"enhanced.max_candidates_per_chart","enhanced.max_solves_per_chart"}) property(schema,path)["minimum"]=1;
    for(const auto* path:{"enhanced.endpoint_slide.max_slide_distance","enhanced.endpoint_slide.w_turn","enhanced.endpoint_slide.w_tangent",
        "enhanced.endpoint_slide.w_slide"})
        property(schema,path)["minimum"]=0;
    property(schema,"enhanced.endpoint_slide.vertex_budget")["minimum"]=1;
    property(schema,"enhanced.endpoint_slide")["description"]="enhanced_v1 sliding seam endpoints on tube-like patches (T30/T31): optional "
        "maximum slide distance along a rim (0 disables it), whole-patch search vertex budget and route cost weights";
    // T34 region proposals: nonnegative dual weights, nonnegative Lloyd steps, positive target count.
    for(const auto* path:{"enhanced.regions.w_normal","enhanced.regions.w_preference","enhanced.regions.lloyd_steps"}) property(schema,path)["minimum"]=0;
    property(schema,"enhanced.regions.max_targets")["minimum"]=1;
    // T35 rings: at least one section sample, quantiles in [0, 1], positive bounds.
    for(const auto* path:{"enhanced.rings.section_samples","enhanced.rings.max_candidates","enhanced.rings.max_charts"}) property(schema,path)["minimum"]=1;
    for(const auto* path:{"enhanced.rings.low_quantile","enhanced.rings.high_quantile","enhanced.rings.min_confidence"}) {
        property(schema,path)["minimum"]=0;property(schema,path)["maximum"]=1;
    }
    // T36 cavities: positive sample, ray, level, cavity and chart bounds; fractions and
    // thresholds in [0, 1]; a cone strictly inside the tangent plane.
    for(const auto* path:{"enhanced.cavities.max_samples","enhanced.cavities.rays_per_sample","enhanced.cavities.max_cavities",
        "enhanced.cavities.min_chart_faces"}) property(schema,path)["minimum"]=1;
    property(schema,"enhanced.cavities.level_samples")["minimum"]=2;
    for(const auto* path:{"enhanced.cavities.reach_fraction","enhanced.cavities.occlusion_threshold","enhanced.cavities.min_area_fraction",
        "enhanced.cavities.max_area_fraction","enhanced.cavities.min_confidence"}) {
        property(schema,path)["minimum"]=0;property(schema,path)["maximum"]=1;
    }
    for(const auto* path:{"enhanced.cavities.min_depth_ratio","enhanced.cavities.slit_depth_ratio"}) property(schema,path)["minimum"]=0;
    property(schema,"enhanced.cavities.cone_half_angle_deg")["minimum"]=0;property(schema,"enhanced.cavities.cone_half_angle_deg")["maximum"]=89;
    property(schema,"enhanced.cavities")["description"]="enhanced_v1 confidence-gated cavity proposals (T36): same-component cone rays on the "
        "geometry BVH (samples, rays, cone, reach as a fraction of the component diagonal), refinement levels, detector filters and bounds";
    // T37 shells: ratios, fractions and gates in [0, 1]; positive bounds; nonnegative passes.
    for(const auto* path:{"enhanced.shells.thin_ratio","enhanced.shells.min_correspondence","enhanced.shells.min_confidence",
        "enhanced.shells.speckle_fraction"}) {
        property(schema,path)["minimum"]=0;property(schema,path)["maximum"]=1;
    }
    for(const auto* path:{"enhanced.shells.max_candidates","enhanced.shells.max_samples","enhanced.shells.min_faces"}) property(schema,path)["minimum"]=1;
    property(schema,"enhanced.shells.smoothing_passes")["minimum"]=0;
    property(schema,"enhanced.shells")["description"]="enhanced_v1 thin-shell front/back rim proposals (T37): nominating thickness ratio, "
        "opposite-side correspondence and confidence gates, label smoothing and speckle repair, alternatives per component (baseline included)";
    // T38 symmetry: explicit policy values, nonnegative tolerances and candidate bound.
    choice("enhanced.symmetry.mode",{"asymmetric","strict","approximate"});
    choice("enhanced.symmetry.on_conflict",{"reject_pair","disable_for_proposal"});
    for(const auto* path:{"enhanced.symmetry.map_tolerance","enhanced.symmetry.geometry_tolerance","enhanced.symmetry.max_candidates"})
        property(schema,path)["minimum"]=0;
    property(schema,"enhanced.symmetry")["description"]="enhanced_v1 operational symmetry (T38): asymmetric (view-optimal, unpaired moves), strict "
        "(every move paired with its mirror image as one transaction on an exact combinatorial correspondence) or approximate (metric-only "
        "diagnostics); the conflict policy, an explicit or detected plane, the counterpart and strict geometry tolerances (fractions of the "
        "bounding-box diagonal) and the symmetrization alternatives per pass";
    // T39 physical-distance policy: explicit object-space lengths are nonnegative (0 selects
    // the frozen reference), the mode is physical or hops, and the step bound is positive.
    choice("enhanced.scale.mode",{"physical","hops"});
    for(const auto* path:{"enhanced.scale.reference_length","enhanced.scale.reference_area","enhanced.scale.corridor_radius","enhanced.scale.corridor_scale"})
        property(schema,path)["minimum"]=0;
    property(schema,"enhanced.scale.max_steps")["minimum"]=1;
    property(schema,"enhanced.scale")["description"]="enhanced_v1 physical-distance policy (T39): physical (geodesic corridors, area floors and "
        "smoothing in object-space units frozen at the operation start; explicit reference_length, reference_area and corridor_radius for cross-tessellation "
        "work) or hops (the legacy hop-ring corridors, face-count floors and step counts, for comparisons with native_baseline)";
    // T40 bounded search: per-feature shares of the per-chart bounds are positive.
    for(const auto* path:{"enhanced.search.feature_candidates_per_chart","enhanced.search.feature_solves_per_chart"}) property(schema,path)["minimum"]=1;
    property(schema,"enhanced.search")["description"]="enhanced_v1 bounded alternative search (T40, Section 21): each feature's share of the "
        "per-chart candidate and solve bounds, and whether a candidate whose exact objective proves it cannot win is refused before its solve";
    // GB (M3): the versioned ordering over the objective; v1 reproduces T40, v2 (the default) adds
    // the exposure-sampling significance margin (exposure_uncertainty.v1).
    choice("enhanced.objective.version",{"seam_objective.v1","seam_objective.v2"});
    property(schema,"enhanced.objective.version")["description"]="the ordering over J: seam_objective.v1 compares J within the tie band only "
        "(T40); seam_objective.v2 also treats a J difference within the sampling uncertainty of the raw exposure term as a tie that keeps "
        "the incumbent";
    // GB (visual review V2): the practical-effect floor of seam_objective.v2, a fraction of the
    // incumbent's J in [0, 0.1]; 0 disables it.
    property(schema,"enhanced.objective.significance.practical_floor")["minimum"]=0;
    property(schema,"enhanced.objective.significance.practical_floor")["maximum"]=0.1;
    property(schema,"enhanced.objective.significance")["description"]="seam_objective.v2 significance rule (enhanced_v1): where the "
        "exposure-sampling margin U applies, a change decided by J is adopted only when it lowers J by more than "
        "max(U, practical_floor * |J of the incumbent|); a change that only removes seam edges has no floor, and budget, compliance and "
        "strict-defect decisions never see it";
    // GB (visual review V1): the baseline chain refinement of compound proposals' new boundaries.
    property(schema,"enhanced.boundary_refinement")["description"]="enhanced_v1: every new or changed seam chain of a region, ring, "
        "cavity or shell proposal goes through the baseline's chain refinement (corridors, turn, flow, drift, convexity and the "
        "crease-discounted exposure) inside the proposal's faces before the transaction evaluates it (requires refine_paths)";
    property(schema,"enhanced.relief")["description"]="enhanced_v1 sector-aware hotspot relief (T32; settings: native record hotspot_relief.v1)";
    property(schema,"enhanced.spur_removal")["description"]="enhanced_v1 validated unnecessary-spur removal (T33; the legacy spur cleanup is unchanged)";
    property(schema,"enhanced.strict_repair")["description"]="enhanced_v1 validated strict-overlap repair (Section 11.3 policy; native record "
        "strict_overlap_repair.v1): candidate zero's charts that fail only the strict overlap check are re-mapped or minimally split through the "
        "transaction; a defect left unrepaired is reported and keeps native_baseline's decision";
    property(schema,"enhanced.constraint_repair")["description"]="enhanced_v1 constraint repair of painted hard constraints (T29, Section 21 "
        "step 1; an ablation switch): off, a candidate zero violating required or forbidden edges gets no repair proposal and such a run is refused "
        "(constraint_conflict); the constraints themselves are never relaxed";
    property(schema,"enhanced.rings")["description"]="enhanced_v1 constriction/crease ring proposals (T35): distance-level sections from a "
        "patch opening, the confidence gate and the candidate bounds per chart and pass";
    property(schema,"enhanced.regions")["description"]="enhanced_v1 preference-aware region proposals (T34): dual crossing weight "
        "centroid_distance * (1 + w_normal*normal_change + w_preference*strength^2), Lloyd steps and targets per pass";
    property(schema,"enhanced.regions.keep_feature_seams")["description"]="GB (visual review V3): the incumbent's seams on hard feature edges "
        "(|dihedral| at least the baseline seam cost's feature angle; never a painted forbidden edge) bound the region growth like required "
        "and locked cuts and stay cut wherever the regrown regions meet across them";
    property(schema,"enhanced")["description"]="enhanced_v1 proposal policy: objective weights (Section 11.2), per-chart search bounds (Section 21) and "
        "strict-validity tolerances; legacy_010 and native_baseline ignore it";
    // GB m8 (G08): every count and step option of the enhanced record has an upper bound, so no
    // option value can make a detector, a smoothing loop or the search allocate or loop without
    // limit. The maxima are generous (at least 16x the default and above every tested value).
    for(const auto& [path,maximum]:std::initializer_list<std::pair<const char*,std::int64_t>>{
            {"enhanced.max_candidates_per_chart",max_candidates_per_chart},{"enhanced.max_solves_per_chart",max_solves_per_chart},
            {"enhanced.search.feature_candidates_per_chart",max_candidates_per_chart},{"enhanced.search.feature_solves_per_chart",max_solves_per_chart},
            {"enhanced.endpoint_slide.vertex_budget",10'000'000},{"enhanced.regions.lloyd_steps",256},{"enhanced.regions.max_targets",1024},
            {"enhanced.rings.section_samples",256},{"enhanced.rings.max_candidates",256},{"enhanced.rings.max_charts",1024},
            {"enhanced.cavities.max_samples",1'000'000},{"enhanced.cavities.rays_per_sample",256},{"enhanced.cavities.level_samples",256},
            {"enhanced.cavities.max_cavities",256},{"enhanced.shells.max_candidates",256},{"enhanced.shells.max_samples",1'000'000},
            {"enhanced.shells.smoothing_passes",256},{"enhanced.symmetry.max_candidates",256},{"enhanced.scale.max_steps",4096}})
        property(schema,path)["maximum"]=maximum;
    return schema;
}
Json ui_metadata() {return Json::parse(
#include "ui_metadata.inc"
);}
// T41: the UiSettings fields added after the frozen source's (metadata group "enhanced"), in
// declaration order: the host declares their properties from this metadata and the defaults
// of the UI schema, and transfers their values; nothing about them is spelled in the host.
Json ui_enhanced_fields() {
    const auto metadata=ui_metadata();Json out=Json::array();
    visit_fields(UiSettings{},[&](const char* name,const auto&) {
        if(metadata.contains(name) && metadata.at(name).value("group","")=="enhanced") out.push_back(name);
    });
    return out;
}
// T41: the panel sections of the enhanced controls, in drawing order. `requires` names the
// policy under which the section applies (the host greys it out otherwise).
Json ui_sections() {
    return Json::array({
        {{"id","mode"},{"label","Seam policy"},{"fields",{"policy"}}},
        {{"id","features"},{"label","Enhanced features"},{"requires","enhanced_v1"},
         {"fields",{"enhanced_sliding","max_slide_distance","enhanced_relief","enhanced_regions","enhanced_rings","ring_confidence",
                    "enhanced_cavities","cavity_confidence","enhanced_shells","shell_confidence","enhanced_spur_removal"}}},
        {{"id","symmetry"},{"label","Symmetry"},{"requires","enhanced_v1"},
         {"fields",{"symmetry_mode","symmetry_on_conflict","symmetry_plane","symmetry_offset"}}},
        {{"id","scale"},{"label","Distances"},{"requires","enhanced_v1"},{"fields",{"scale_mode","reference_length"}}},
        {{"id","painted"},{"label","Painted attributes"},{"requires","enhanced_v1"},
         {"fields",{"required_attribute","forbidden_attribute","allow_attribute","avoid_attribute"}}},
        {{"id","search"},{"label","Search bounds"},{"requires","enhanced_v1"},
         {"fields",{"max_candidates_per_chart","max_solves_per_chart","prune_dominated"}}}});
}
Json ui_schema() {
    auto schema=detail::schema_for(UiSettings{});
    const auto metadata=ui_metadata();
    for(auto it=metadata.begin();it!=metadata.end();++it) {
        auto& field=schema["properties"][it.key()];const auto& info=it.value();
        if(info.contains("min")) field["minimum"]=info["min"];
        if(info.contains("max")) field["maximum"]=info["max"];
        if(info.contains("items")) {field["enum"]=Json::array();for(const auto& item:info["items"]) field["enum"].push_back(item[0]);}
        if(info.contains("description")) field["description"]=info["description"];
        if(info.contains("name")) field["title"]=info["name"];
    }
    return schema;
}
void check_options(const PipelineOptions& o) {
    if(o.enhanced.rings.low_quantile>o.enhanced.rings.high_quantile)
        throw ContractError(ErrorCode::invalid_option,"enhanced.rings.low_quantile","the low quantile exceeds the high quantile");
    if(o.enhanced.cavities.min_area_fraction>o.enhanced.cavities.max_area_fraction)
        throw ContractError(ErrorCode::invalid_option,"enhanced.cavities.min_area_fraction","the minimum area fraction exceeds the maximum");
    if(o.enhanced.symmetry.plane_normal) {
        const auto& n=*o.enhanced.symmetry.plane_normal;
        if(!(n[0]*n[0]+n[1]*n[1]+n[2]*n[2]>0))
            throw ContractError(ErrorCode::invalid_option,"enhanced.symmetry.plane_normal","the plane normal must be a nonzero vector");
    }
    // Review m2: one host attribute cannot carry two painted roles (it would be resolved by
    // priority). Only enhanced_v1 reads the names; the other policies ignore them (M2).
    if(o.policy=="enhanced_v1") {
        const auto& p=o.enhanced.preferences;std::map<std::string,std::string> roles;
        for(const auto& [field,name]:{std::pair{"enhanced.preferences.required_attribute",&p.required_attribute},
                                      std::pair{"enhanced.preferences.forbidden_attribute",&p.forbidden_attribute},
                                      std::pair{"enhanced.preferences.allow_attribute",&p.allow_attribute},
                                      std::pair{"enhanced.preferences.avoid_attribute",&p.avoid_attribute}})
            if(!name->empty() && !roles.emplace(*name,field).second)
                throw ContractError(ErrorCode::invalid_option,field,"attribute \""+*name+"\" is already configured as "+roles.at(*name));
    }
    if(o.target_islands && o.target_islands->at(0)>o.target_islands->at(1))
        throw ContractError(ErrorCode::invalid_option,"target_islands","minimum exceeds maximum");
    // FIX-02: the corrected routing searches radii up to max(corridor_rings, corridor_max_rings)
    // and requires both positive; legacy_010 keeps the source's search of any radius.
    if(o.policy!="legacy_010")
        for(const auto& [field,value]:{std::pair{"paths.corridor_rings",o.paths.corridor_rings},std::pair{"paths.corridor_max_rings",o.paths.corridor_max_rings}})
            if(value<1) throw ContractError(ErrorCode::invalid_option,field,"corridor radii must be positive integers under "+o.policy+" (FIX-02)");
    std::set<std::int64_t> tiles(o.packing.tiles.begin(),o.packing.tiles.end());
    if(tiles.size()!=o.packing.tiles.size()) throw ContractError(ErrorCode::invalid_option,"packing.tiles","duplicate tile");
    if(o.layout.hide_dir) {
        const auto& d=*o.layout.hide_dir;
        if(std::hypot(d[0],d[1],d[2])==0) throw ContractError(ErrorCode::invalid_option,"layout.hide_dir","direction must be nonzero");
    }
}
// Source exposure.view_direction: vector / np.linalg.norm(vector), the norm being
// sqrt(x.dot(x)) summed in order (the squares of Blender's float32 values are exact).
Vec3 normalized(Vec3 d) {
    const double squares=(d[0]*d[0]+d[1]*d[1])+d[2]*d[2];
    const auto norm=std::sqrt(squares);
    if(!std::isfinite(norm) || norm<=0) throw ContractError(ErrorCode::invalid_option,"hide_direction","finite nonzero direction required");
    for(auto& v:d) v/=norm;return d;
}
// Source ops.build_options for the current view, in mathutils float32 (Gate A F21):
// obj.matrix_world.inverted().to_3x3() @ view. mathutils inverts a 4x4 matrix as
// Blender's adjoint_m4_m4 / determinant_m4 (classic cofactor formulas, one float
// rounding per operation) and multiplies Matrix @ Vector with float products summed in a
// double. The host passes the world matrix's linear part and the view vector as Blender's
// float32 values; the translation of an affine matrix meets only exact zeros in the
// linear block, so it is taken as zero. Bit-exact with Blender 5.2.1 on random
// transforms.
Vec3 mathutils_object_direction(const ViewContext& view) {
    float m[4][4]={}; // BLI layout m[column][row]
    for(int r=0;r<3;++r)for(int c=0;c<3;++c)m[c][r]=static_cast<float>(view.object_to_world[static_cast<std::size_t>(r*3+c)]);
    m[3][3]=1.f;
    float v[3];for(std::size_t i=0;i<3;++i)v[i]=static_cast<float>(view.direction_world[i]);
    for(const auto& column:m)for(const float x:column)
        if(!std::isfinite(x)) throw ContractError(ErrorCode::invalid_option,"view","finite invertible object transform required");
    for(const float x:v) if(!std::isfinite(x)) throw ContractError(ErrorCode::invalid_option,"view","finite view direction required");
    auto det2=[](float a,float b,float c,float d) {const float ad=a*d;const float bc=b*c;return ad-bc;};
    auto det3=[&](float a1,float a2,float a3,float b1,float b2,float b3,float c1,float c2,float c3) {
        const float t1=a1*det2(b2,b3,c2,c3);const float t2=b1*det2(a2,a3,c2,c3);const float t3=c1*det2(a2,a3,b2,b3);
        const float t12=t1-t2;return t12+t3;
    };
    const float a1=m[0][0],b1=m[0][1],c1=m[0][2],d1=m[0][3],a2=m[1][0],b2=m[1][1],c2=m[1][2],d2=m[1][3];
    const float a3=m[2][0],b3=m[2][1],c3=m[2][2],d3=m[2][3],a4=m[3][0],b4=m[3][1],c4=m[3][2],d4=m[3][3];
    const float x1=a1*det3(b2,b3,b4,c2,c3,c4,d2,d3,d4);const float x2=b1*det3(a2,a3,a4,c2,c3,c4,d2,d3,d4);
    const float x3=c1*det3(a2,a3,a4,b2,b3,b4,d2,d3,d4);const float x4=d1*det3(a2,a3,a4,b2,b3,b4,c2,c3,c4);
    const float x12=x1-x2;const float x123=x12+x3;const float det=x123-x4;
    if(det==0.f || !std::isfinite(det)) throw ContractError(ErrorCode::invalid_option,"view","finite invertible object transform required");
    float adj[3][3]; // adjoint_m4_m4's linear block, BLI layout adj[column][row]
    adj[0][0]=det3(b2,b3,b4,c2,c3,c4,d2,d3,d4);adj[1][0]=-det3(a2,a3,a4,c2,c3,c4,d2,d3,d4);adj[2][0]=det3(a2,a3,a4,b2,b3,b4,d2,d3,d4);
    adj[0][1]=-det3(b1,b3,b4,c1,c3,c4,d1,d3,d4);adj[1][1]=det3(a1,a3,a4,c1,c3,c4,d1,d3,d4);adj[2][1]=-det3(a1,a3,a4,b1,b3,b4,d1,d3,d4);
    adj[0][2]=det3(b1,b2,b4,c1,c2,c4,d1,d2,d4);adj[1][2]=-det3(a1,a2,a4,c1,c2,c4,d1,d2,d4);adj[2][2]=det3(a1,a2,a4,b1,b2,b4,d1,d2,d4);
    Vec3 out{};
    for(std::size_t row=0;row<3;++row) {
        double dot=0.0;
        for(std::size_t col=0;col<3;++col) {const float inverse=adj[col][row]/det;const float product=inverse*v[col];dot+=static_cast<double>(product);}
        out[row]=static_cast<double>(static_cast<float>(dot));
    }
    return out;
}
// Python int() of one tile text (int(tile.strip())): surrounding
// whitespace (str.isspace in ASCII), an optional sign, decimal digits with single underscores
// between them (PEP 515) and leading zeros. Non-ASCII digits and whitespace are refused.
std::int64_t python_int(std::string_view text) {
    auto bad=[](const char* message)->std::int64_t {throw ContractError(ErrorCode::invalid_option,"pack_tiles",message);};
    auto space=[](char c) {return c==' ' || (c>='\t' && c<='\r') || (c>='\x1c' && c<='\x1f');};
    while(!text.empty() && space(text.front())) text.remove_prefix(1);
    while(!text.empty() && space(text.back())) text.remove_suffix(1);
    if(text.empty()) return bad("empty tile");
    const bool negative=text.front()=='-';
    if(text.front()=='+' || text.front()=='-') text.remove_prefix(1);
    std::int64_t value=0;bool digits=false,underscore=false;
    for(const char c:text) {
        if(c=='_') {if(!digits || underscore) return bad("invalid tile integer");underscore=true;continue;}
        if(c<'0' || c>'9') return bad("invalid tile integer");
        const auto d=c-'0';
        if(value>(std::numeric_limits<std::int64_t>::max()-d)/10) return bad("tile integer out of range");
        value=value*10+d;digits=true;underscore=false;
    }
    if(!digits || underscore) return bad("invalid tile integer");
    return negative?-value:value;
}
std::vector<std::int64_t> python_tiles(const std::string& text) {
    std::vector<std::int64_t> tiles;std::size_t start=0;
    while(true) {
        const auto end=text.find(',',start);
        tiles.push_back(python_int(std::string_view(text).substr(start,end==std::string::npos?std::string::npos:end-start)));
        if(end==std::string::npos) return tiles;
        start=end+1;
    }
}
}
PipelineOptions pipeline_options_from_json(const Json& input) {
    const auto complete=detail::merge_defaults(detail::encode(PipelineOptions{}),input);
    detail::validate_schema(complete,pipeline_schema(),"");
    auto options=detail::decode<PipelineOptions>(complete,"");check_options(options);return options;
}
Json pipeline_options_json(const PipelineOptions& options) {
    auto value=detail::encode(options);detail::validate_schema(value,pipeline_schema(),"");check_options(options);return value;
}
UiSettings ui_settings_from_json(const Json& input) {
    const auto complete=detail::merge_defaults(detail::encode(UiSettings{}),input);
    detail::validate_schema(complete,ui_schema(),"ui");return detail::decode<UiSettings>(complete,"ui");
}
Json ui_settings_json(const UiSettings& settings) {
    auto value=detail::encode(settings);detail::validate_schema(value,ui_schema(),"ui");return value;
}
Json option_registry() {
    return {{"version",1},{"pipeline",pipeline_schema()},{"ui",ui_schema()},{"ui_metadata",ui_metadata()},
        {"ui_enhanced_fields",ui_enhanced_fields()},{"ui_sections",ui_sections()},
        {"policy_status",{{"repair_v1","development: supplied-layout validation; optional recovery algorithms pending"},{"legacy_010","reference compatibility"},{"native_baseline","corrections implemented by T23; the product default"},
                          {"enhanced_v1","Gate B: the native_baseline result is candidate zero of the bounded proposal pass (T28-T40); selectable in "
                           "the CLI (--behavior, options) and the Blender panel (T41); the default stays native_baseline until Gate B (T42)"}}}};
}
UiSettings apply_ui_profile(UiSettings settings,const std::string& profile) {
    settings.profile=profile;
    if(profile=="hard_surface") {settings.feature_angle=60.;settings.contour_mode="auto";}
    else if(profile=="organic") {settings.feature_angle=0.;settings.contour_mode="organic";}
    ui_settings_json(settings);return settings;
}
PipelineOptions apply_classified_profile(PipelineOptions options,const std::string& selected) {
    // Source asset_class returns hard, organic or tri; tri keeps the hard layout
    // settings with organic contours (pipeline.py auto profile).
    if(selected!="hard" && selected!="organic" && selected!="tri") throw ContractError(ErrorCode::invalid_option,"selected_profile","unknown classification");
    if(options.profile=="auto" && options.generate) {
        options.layout.feature_angle=selected=="organic"?0.:60.;options.layout.prefer_features=selected=="organic";
        options.boundary.mode=selected=="hard"?"auto":"organic";
    }
    pipeline_options_json(options);return options;
}
PipelineOptions resolve_ui_options(const UiSettings& st,const std::optional<ViewContext>& view) {
    ui_settings_json(st);
    PipelineOptions o;
    const auto aniso=st.budget=="tight"?1.3:st.budget=="normal"?1.6:2.;
    const auto area=st.budget=="tight"?1.4:st.budget=="normal"?1.8:2.4;
    o.profile=st.profile=="auto"?"auto":"custom";o.generate=st.generate;o.refine_paths=st.refine_paths;
    if(st.target_enabled) o.target_islands=IslandBand{st.target_min,st.target_max};
    if(st.deadline_s!=0) o.deadline_s=st.deadline_s;
    o.clean_contours=st.contour_mode!="off";
    o.layout.feature_angle=st.feature_angle;o.layout.max_anisotropy_p95=aniso;o.layout.max_area_p95=area;
    o.layout.max_anisotropy_p99=std::max(4.,aniso*2.5);o.layout.min_island_fraction=st.min_island/100.;
    o.layout.min_island_faces=st.min_island_faces;o.layout.min_fusion=st.min_fusion;
    o.paths.corridor_rings=st.corridor_rings;o.paths.w_turn=st.w_turn;o.paths.w_flow=st.w_flow;o.paths.w_drift=st.w_drift;
    o.boundary.mode=st.contour_mode=="off"?"auto":st.contour_mode;o.boundary.smooth_radius=st.smooth_radius;
    o.boundary.max_displacement=st.max_displacement;
    o.abs_anisotropy_p95=o.boundary.abs_anisotropy_p95=aniso*1.25;o.abs_area_p95=o.boundary.abs_area_p95=area*1.25;
    if(st.hide_from=="custom") o.layout.hide_dir=normalized(st.hide_direction);
    else if(st.hide_from=="view") {
        if(!view) throw ContractError(ErrorCode::invalid_option,"hide_from","current view requires host view data");
        o.layout.hide_dir=normalized(mathutils_object_direction(*view));
    } else if(st.hide_from!="none") {
        Vec3 d{};const auto axis=std::string("XYZ").find(st.hide_from[1]);d[axis]=st.hide_from[0]=='+'?1.:-1.;o.layout.hide_dir=d;
    }
    o.solver.method="blender";o.solver.iterations=st.solver_iterations;
    o.host.protect_existing=st.protect_existing;o.host.pack_after=st.pack_after;
    o.host.show_advanced=st.show_advanced;o.host.last_status=st.last_status;o.host.last_report=st.last_report;
    // T41: the policy and the enhanced controls. The record is copied under every policy and
    // read only by enhanced_v1; the UI defaults are the record's own, so settings saved before
    // these fields existed resolve exactly as before (native_baseline, default record).
    o.policy=st.policy;
    auto& e=o.enhanced;
    e.endpoint_slide.enabled=st.enhanced_sliding;e.endpoint_slide.max_slide_distance=st.max_slide_distance;
    e.relief.enabled=st.enhanced_relief;e.regions.enabled=st.enhanced_regions;
    e.rings.enabled=st.enhanced_rings;e.rings.min_confidence=static_cast<double>(st.ring_confidence)/100.;
    e.cavities.enabled=st.enhanced_cavities;e.cavities.min_confidence=static_cast<double>(st.cavity_confidence)/100.;
    e.shells.enabled=st.enhanced_shells;e.shells.min_confidence=static_cast<double>(st.shell_confidence)/100.;
    e.spur_removal.enabled=st.enhanced_spur_removal;
    e.symmetry.mode=st.symmetry_mode;e.symmetry.on_conflict=st.symmetry_on_conflict;
    if(st.symmetry_plane!="detect") {
        // An explicit object-axis plane {x : x[axis] = offset}; "detect" leaves the plane to T14.
        Vec3 normal{};normal[std::string("XYZ").find(st.symmetry_plane)]=1.;
        e.symmetry.plane_normal=normal;e.symmetry.plane_offset=st.symmetry_offset;
    }
    e.scale.mode=st.scale_mode;e.scale.reference_length=st.reference_length;
    e.preferences.required_attribute=st.required_attribute;e.preferences.forbidden_attribute=st.forbidden_attribute;
    e.preferences.allow_attribute=st.allow_attribute;e.preferences.avoid_attribute=st.avoid_attribute;
    e.max_candidates_per_chart=st.max_candidates_per_chart;e.max_solves_per_chart=st.max_solves_per_chart;
    e.search.prune_dominated=st.prune_dominated;
    // Packing settings are resolved only for packing (resolve_ui_pack_options): the source's
    // build_options never reads them, so a malformed tile list cannot block a run (Gate A F8).
    pipeline_options_json(o);return o;
}
PipelineOptions resolve_ui_pack_options(const UiSettings& st) {
    // Resolve only packing settings; the hide direction, island
    // target and budgets are run settings and are neither resolved nor validated here.
    ui_settings_json(st);
    PipelineOptions o;
    o.host.pack_after=st.pack_after;
    o.packing.texture_width=st.pack_width;o.packing.texture_height=st.pack_height;
    o.packing.padding_px=st.pack_padding;o.packing.border_px=st.pack_border;o.packing.keep_orientation=st.keep_orientation;
    o.packing.tiles=python_tiles(st.pack_tiles);
    pipeline_options_json(o);return o; // PackOptions checks tile ranges.
}
}
