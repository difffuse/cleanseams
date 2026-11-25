// SPDX-License-Identifier: GPL-3.0-or-later
// Phase B confidence-gated cavity proposals (T36, Section 17): recessed regions nominated by
// short same-component cone rays on the geometry-revisioned BVH, refined on the original
// mesh into one coherent opening with measured depth, and proposed as compound cuts (a seam
// loop just inside the rim, optionally an interior slit from the deepest point, plus every
// companion opening) through the common transaction. Low-confidence, open, thin or
// through-hole cases keep the baseline and report why. Used by policy enhanced_v1 only.
#include "cleanseams/proposals.hpp"
#include "cleanseams/physical_scale.hpp"
#include "cleanseams/geometry_bvh.hpp"
#include "cleanseams/cache.hpp"
#include "cleanseams/errors.hpp"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <numbers>
#include <numeric>
#include <queue>
#include <set>
namespace cleanseams {
namespace {
using Json=nlohmann::json;
constexpr const char* name="cavities";
constexpr const char* version="1";
bool always_cut(const MeshTables& t,std::size_t e) {
    return t.edge_boundary[e] || !t.edge_manifold[e] || !t.edge_consistent[e] || t.edge_corners[e].size()!=2;
}
Vec3 add(const Vec3& a,const Vec3& b) {return {a[0]+b[0],a[1]+b[1],a[2]+b[2]};}
Vec3 sub(const Vec3& a,const Vec3& b) {return {a[0]-b[0],a[1]-b[1],a[2]-b[2]};}
Vec3 mul(const Vec3& a,double s) {return {a[0]*s,a[1]*s,a[2]*s};}
double dot(const Vec3& a,const Vec3& b) {return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];}
Vec3 cross(const Vec3& a,const Vec3& b) {return {a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]};}
double norm(const Vec3& a) {return std::sqrt(dot(a,a));}
Vec3 unit(const Vec3& a) {const double n=norm(a);return n>0?mul(a,1/n):Vec3{0,0,0};}
double length_of(const MeshSnapshot& m,std::span<const MeshIndex> edges) {double s=0;for(auto e:edges)s+=m.geometry().edge_length[e];return s;}
Json edge_list(std::span<const MeshIndex> edges) {return Json(std::vector<MeshIndex>(edges.begin(),edges.end()));}
// GB m7: cheapest_edge_path (compound.cpp) restricted to the sorted vertex set `vertices`
// that every usable edge stays in: the same (distance, vertex) order, relaxation order,
// settled rule and tie-break on equal distances, with region-sized arrays instead of mesh-sized.
std::vector<MeshIndex> local_cheapest_path(const MeshSnapshot& mesh,std::span<const MeshIndex> vertices,MeshIndex source,std::span<const std::uint8_t> targets,
                                           const std::function<bool(MeshIndex)>& usable,std::span<const double> cost) {
    const auto& t=mesh.tables();
    auto slot=[&](MeshIndex v) {
        const auto it=std::lower_bound(vertices.begin(),vertices.end(),v);
        if(it==vertices.end() || *it!=v)throw ContractError(ErrorCode::invalid_input,"cavities.path","a usable edge leaves the region");
        return static_cast<std::size_t>(it-vertices.begin());
    };
    std::vector<double> dist(vertices.size(),std::numeric_limits<double>::infinity());std::vector<std::int64_t> via(vertices.size(),-1);
    std::vector<std::uint8_t> settled(vertices.size(),0);
    using Item=std::pair<double,MeshIndex>;std::priority_queue<Item,std::vector<Item>,std::greater<>> heap;
    dist[slot(source)]=0;heap.push({0.0,source});
    while(!heap.empty()) {
        const auto [d,v]=heap.top();heap.pop();const auto iv=slot(v);
        if(d>dist[iv] || settled[iv])continue;
        settled[iv]=1;
        if(targets[iv]) {
            std::vector<MeshIndex> edges;auto at=v;
            while(via[slot(at)]>=0) {
                const auto e=static_cast<MeshIndex>(via[slot(at)]);edges.push_back(e);at=mesh.edge_other(e,at);
                if(edges.size()>vertices.size())throw ContractError(ErrorCode::numerical_failure,"cavities.path","predecessor cycle");
            }
            std::reverse(edges.begin(),edges.end());
            return edges;
        }
        for(auto e:t.vert_edges[v]) {
            if(!usable(e))continue;
            const double w=cost[e];
            if(!std::isfinite(w) || w<0)throw ContractError(ErrorCode::invalid_input,"compound.path","path costs must be finite and nonnegative");
            const auto u=mesh.edge_other(e,v);const auto iu=slot(u);if(settled[iu])continue;
            const double nd=d+w;
            if(nd<dist[iu] || (nd==dist[iu] && via[iu]>=0 && e<static_cast<MeshIndex>(via[iu]))) {dist[iu]=nd;via[iu]=e;heap.push({nd,u});}
        }
    }
    return {};
}
}
std::vector<double> component_diagonals(const MeshSnapshot& m) {
    const auto& g=m.geometry();const auto& in=m.input();
    std::vector<std::array<Vec3,2>> box(g.n_components,{Vec3{INFINITY,INFINITY,INFINITY},Vec3{-INFINITY,-INFINITY,-INFINITY}});
    for(std::size_t f=0;f<g.face_area.size();++f)for(auto c:m.face_corners(static_cast<MeshIndex>(f))) {
        const auto& p=in.positions[static_cast<std::size_t>(in.corner_vertices[c])];auto& b=box[g.face_component[f]];
        for(std::size_t k=0;k<3;++k) {b[0][k]=std::min(b[0][k],p[k]);b[1][k]=std::max(b[1][k],p[k]);}
    }
    std::vector<double> out;for(const auto& b:box)out.push_back(b[0][0]<=b[1][0]?norm(sub(b[1],b[0])):0.0);
    return out;
}
nlohmann::json CavityOcclusion::to_json() const {
    Json reaches=Json::object();for(const auto& [k,r]:reach)reaches[std::to_string(k)]=r;
    return {{"samples",samples.size()},{"rays",rays},{"hits",hits},{"escaped",rays-hits},{"reach",reaches},{"seconds",seconds}};
}
CavityOcclusion cavity_occlusion(const MeshSnapshot& m,const TriangleBvh& bvh,std::span<const MeshIndex> faces,const CavityOptions& o,bool same_component,
                                 std::span<const double> given_diagonals) {
    const auto& g=m.geometry();const auto& t=m.tables();const auto& in=m.input();
    if(!bvh.serves(m))throw ContractError(ErrorCode::stale_state,"cavities.bvh","the BVH belongs to another geometry revision");
    CavityOcclusion out;const auto started=std::chrono::steady_clock::now();
    // GB m7: the component diagonals are a pass-wide quantity (computed here only without one).
    if(!given_diagonals.empty() && given_diagonals.size()!=g.n_components)
        throw ContractError(ErrorCode::invalid_input,"cavities.diagonals","one diagonal per connected component required");
    std::vector<double> computed;if(given_diagonals.empty())computed=component_diagonals(m);
    const std::span<const double> diagonals=given_diagonals.empty()?std::span<const double>(computed):given_diagonals;
    // Area-distributed samples: every usable face when there are at most max_samples,
    // otherwise systematic sampling of the cumulative area (deterministic).
    std::vector<MeshIndex> usable;
    for(auto f:faces) {
        if(f>=g.face_area.size())throw ContractError(ErrorCode::invalid_input,"cavities.faces","face outside the mesh");
        if(g.face_area[f]>0 && g.face_degenerate_kind[f]==0)usable.push_back(f);
    }
    std::sort(usable.begin(),usable.end());usable.erase(std::unique(usable.begin(),usable.end()),usable.end());
    const auto cap=static_cast<std::size_t>(std::max<std::int64_t>(1,o.max_samples));
    if(usable.size()<=cap)out.samples=usable;
    else {
        std::vector<double> cumulative;double total=0;for(auto f:usable) {total+=g.face_area[f];cumulative.push_back(total);}
        for(std::size_t k=0;k<cap;++k) {
            const double at=(static_cast<double>(k)+0.5)*total/static_cast<double>(cap);
            const auto i=static_cast<std::size_t>(std::lower_bound(cumulative.begin(),cumulative.end(),at)-cumulative.begin());
            const auto f=usable[std::min(i,usable.size()-1)];if(out.samples.empty() || out.samples.back()!=f)out.samples.push_back(f);
        }
    }
    const auto rays=static_cast<std::size_t>(std::max<std::int64_t>(1,o.rays_per_sample));
    const double angle=o.cone_half_angle_deg*std::numbers::pi/180.0;
    for(auto f:out.samples) {
        const auto n=g.face_normal[f];const auto component=g.face_component[f];
        const double reach=o.reach_fraction*diagonals[component];out.reach.emplace(component,reach);
        // Tangent frame from the face's first edge (deterministic, moves with the geometry).
        const auto c0=m.face_corners(f).front();
        const auto& a=in.positions[static_cast<std::size_t>(in.corner_vertices[c0])];
        const auto& b=in.positions[static_cast<std::size_t>(in.corner_vertices[t.corner_next[c0]])];
        auto u=sub(b,a);u=unit(sub(u,mul(n,dot(u,n))));
        if(norm(u)==0)u=unit(cross(n,std::abs(n[0])<0.9?Vec3{1,0,0}:Vec3{0,1,0}));
        const auto w=cross(n,u);
        const auto origin=add(g.face_center[f],mul(n,1e-9*reach));
        std::uint64_t hits=0;
        for(std::size_t k=0;k<rays;++k) {
            const double phi=2*std::numbers::pi*static_cast<double>(k)/static_cast<double>(rays);
            const auto d=unit(add(mul(n,std::cos(angle)),mul(add(mul(u,std::cos(phi)),mul(w,std::sin(phi))),std::sin(angle))));
            // Only the sample's own connected component can close a pocket around it.
            const auto hit=bvh.ray(origin,d,reach,[&](MeshIndex h) {return h!=f && (!same_component || g.face_component[h]==component);});
            ++out.rays;if(hit.hit) {++hits;++out.hits;}
        }
        out.occlusion.emplace(f,static_cast<double>(hits)/static_cast<double>(rays));
    }
    out.seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count();
    return out;
}
nlohmann::json CavityRecord::to_json() const {
    return {{"status",status},{"core_faces",core.size()},{"region_faces",region.size()},{"rim_edges",rim.size()},{"seam",seam},{"slit",slit},
            {"axis",axis},{"occlusion",occlusion},{"concavity",concavity},{"depth",depth},
            {"opening_diameter",opening_diameter},{"depth_ratio",depth_ratio},{"area_fraction",area_fraction},{"confidence",confidence},
            {"seam_depth_ratio",seam_depth_ratio},{"needs_slit",needs_slit},{"seam_inside_rim",seam_inside_rim},
            {"rim_alternative",rim_alternative},{"rim_slit",rim_slit},{"rim_needs_slit",rim_needs_slit}};
}
nlohmann::json CavityDetection::to_json() const {
    Json list=Json::array();for(const auto& c:cavities)list.push_back(c.to_json());
    return {{"version",cavity_detector_version},{"occlusion",occlusion.to_json()},{"seeds",seeds},{"cores",cores},{"cavities",list},{"scale",scale},{"stopped",stopped}};
}
CavityDetection detect_cavities(const MeshSnapshot& m,const TriangleBvh& bvh,std::span<const MeshIndex> chart_faces,std::span<const std::uint8_t> cuts,
                                const HardConstraints& hc,std::span<const double> routing,const PipelineOptions& options,bool same_component,
                                const PhysicalScale* scale,DetectorPass* pass) {
    const auto& t=m.tables();const auto& g=m.geometry();const auto& in=m.input();const auto& o=options.enhanced.cavities;
    const auto nf=g.face_area.size(),ne=t.edge_vertices.size();
    if(cuts.size()!=ne)throw ContractError(ErrorCode::invalid_input,"cavities.cuts","cut mask does not match the mesh edges");
    if(!routing.empty() && routing.size()!=ne)throw ContractError(ErrorCode::invalid_input,"cavities.routing","routing cost needs one value per edge");
    DetectorPass local_pass;DetectorPass& work=pass?*pass:local_pass;
    CavityDetection out;
    // GB m7: every array below is local to the chart (its sorted faces; `local` maps a face to
    // its slot, -1 outside), never mesh-sized: a chart's work is proportional to the chart.
    std::vector<MeshIndex> ordered(chart_faces.begin(),chart_faces.end());std::sort(ordered.begin(),ordered.end());
    ordered.erase(std::unique(ordered.begin(),ordered.end()),ordered.end());
    for(auto f:ordered)if(f>=nf)throw ContractError(ErrorCode::invalid_input,"cavities.faces","face outside the mesh");
    const auto n=ordered.size();
    auto local=[&](MeshIndex f) -> std::ptrdiff_t {const auto it=std::lower_bound(ordered.begin(),ordered.end(),f);return it!=ordered.end() && *it==f?it-ordered.begin():-1;};
    auto slot=[&](MeshIndex f) {return static_cast<std::size_t>(local(f));};
    double chart_area=0;for(auto f:chart_faces)chart_area+=g.face_area[f];
    // Chart adjacency across uncut interior edges (in corner order); faces touching the chart's boundary.
    std::vector<std::vector<std::pair<std::size_t,MeshIndex>>> adjacent(n);std::vector<std::uint8_t> on_boundary(n,0);
    for(std::size_t i=0;i<n;++i) {
        const auto f=ordered[i];
        for(auto c:m.face_corners(f)) {
            const auto e=t.corner_edge[c];const bool open=!cuts[e] && !always_cut(t,e);bool inner=open;
            for(auto h:t.edge_faces[e]) {
                const auto j=local(h);inner=inner && j>=0;
                if(open && h!=f && j>=0)adjacent[i].emplace_back(static_cast<std::size_t>(j),e);
            }
            if(!inner)on_boundary[i]=1;
        }
    }
    if(work.stop()) {out.stopped=true;return out;}
    out.occlusion=cavity_occlusion(m,bvh,chart_faces,o,same_component,work.diagonals(m));
    // Seeds: sampled faces whose cone rays are blocked at least occlusion_threshold of the
    // time; a face surrounded mostly by seeds joins them (holes filled, seeds never removed).
    std::vector<std::uint8_t> seed(n,0);
    for(const auto& [f,x]:out.occlusion.occlusion)if(x>0 && x>=o.occlusion_threshold)seed[slot(f)]=1;
    // T39: three face-ring passes at the reference length (a physical filling radius).
    std::int64_t fill_passes=3;
    if(scale && !chart_faces.empty()) {const auto fill=scale->rings(3,chart_faces);fill_passes=fill.steps;out.scale={{"fill",fill.to_json()}};}
    for(std::int64_t fill=0;fill<fill_passes;++fill) {
        std::vector<std::size_t> grow;
        for(std::size_t i=0;i<n;++i) {
            if(seed[i])continue;std::uint64_t yes=0,all=0;
            for(const auto& [j,e]:adjacent[i]) {++all;yes+=seed[j];}
            if(all && 2*yes>all)grow.push_back(i);
        }
        if(grow.empty())break;for(auto i:grow)seed[i]=1;
    }
    for(auto f:chart_faces)out.seeds+=seed[slot(f)];
    // Cores: connected seed components, strongest occlusion mass first (ties: smallest face).
    std::vector<std::vector<MeshIndex>> cores;std::vector<std::uint8_t> seen(n,0);
    for(std::size_t i0=0;i0<n;++i0) {
        if(!seed[i0] || seen[i0])continue;
        std::vector<MeshIndex> comp;std::vector<std::size_t> stack{i0};seen[i0]=1;
        while(!stack.empty()) {
            const auto i=stack.back();stack.pop_back();comp.push_back(ordered[i]);
            for(const auto& [j,e]:adjacent[i])if(seed[j] && !seen[j]) {seen[j]=1;stack.push_back(j);}
        }
        std::sort(comp.begin(),comp.end());cores.push_back(std::move(comp));
    }
    auto mass=[&](const std::vector<MeshIndex>& core) {
        double s=0;for(auto f:core) {const auto it=out.occlusion.occlusion.find(f);s+=g.face_area[f]*(it==out.occlusion.occlusion.end()?0.0:it->second);}
        return s;
    };
    std::stable_sort(cores.begin(),cores.end(),[&](const auto& a,const auto& b) {return mass(a)>mass(b);});
    out.cores=cores.size();
    auto forbidden=[&](std::span<const MeshIndex> edges) {
        for(auto e:edges)if(e<hc.forbidden.size() && hc.forbidden[e])return true;
        return false;
    };
    for(const auto& core:cores) {
        if(static_cast<std::int64_t>(out.cavities.size())>=o.max_cavities)break;
        if(work.stop()) {out.stopped=true;break;}
        CavityRecord r;r.core=core;
        Vec3 normal{0,0,0},centre{0,0,0};double area=0,occ=0,occ_area=0;
        for(auto f:core) {
            const double a=g.face_area[f];area+=a;normal=add(normal,mul(g.face_normal[f],a));centre=add(centre,mul(g.face_center[f],a));
            if(const auto it=out.occlusion.occlusion.find(f);it!=out.occlusion.occlusion.end()) {occ+=a*it->second;occ_area+=a;}
        }
        r.occlusion=occ_area>0?occ/occ_area:0.0;
        if(!(area>0) || norm(normal)<=1e-9*area) {r.status="no_opening_axis";out.cavities.push_back(std::move(r));continue;}
        r.axis=unit(normal);centre=mul(centre,1/area);
        // Height along the opening axis (normals point out of the recess, towards its opening).
        auto height=[&](const Vec3& p) {return dot(sub(p,centre),r.axis);};
        std::vector<double> h(n,INFINITY);double h_top=-INFINITY;
        for(auto f:chart_faces) {h[slot(f)]=height(g.face_center[f]);h_top=std::max(h_top,h[slot(f)]);}
        // Regions grow from the deepest core face, so loops inside the rim exist too.
        MeshIndex deepest_face=core.front();for(auto f:core)if(h[slot(f)]<h[slot(deepest_face)])deepest_face=f;
        const double h_low=h[slot(deepest_face)];
        struct Level {double t=0;std::vector<MeshIndex> region,border;bool valid=false;double area=0;};
        std::vector<Level> levels;
        const auto samples=static_cast<std::size_t>(std::max<std::int64_t>(2,o.level_samples));
        for(std::size_t k=0;k<samples;++k) {
            Level L;L.t=h_low+(h_top-h_low)*static_cast<double>(k)/static_cast<double>(samples-1);
            std::vector<std::uint8_t> mark(n,0);std::vector<MeshIndex> stack{deepest_face};mark[slot(deepest_face)]=1;
            while(!stack.empty()) {
                const auto f=stack.back();stack.pop_back();L.region.push_back(f);
                for(const auto& [j,e]:adjacent[slot(f)])if(!mark[j] && h[j]<=L.t) {mark[j]=1;stack.push_back(ordered[j]);}
            }
            std::sort(L.region.begin(),L.region.end());bool touches=false;
            for(auto f:L.region) {
                L.area+=g.face_area[f];touches=touches || on_boundary[slot(f)];
                for(const auto& [j,e]:adjacent[slot(f)])if(!mark[j])L.border.push_back(e);
            }
            std::sort(L.border.begin(),L.border.end());L.border.erase(std::unique(L.border.begin(),L.border.end()),L.border.end());
            L.valid=!touches && !L.border.empty() && simple_loop_order(m,L.border).size()==L.border.size() && L.area<=o.max_area_fraction*chart_area;
            levels.push_back(std::move(L));
        }
        std::optional<std::size_t> best;for(std::size_t k=0;k<levels.size();++k)if(levels[k].valid)best=k;
        if(!best) {r.status="incoherent_opening";out.cavities.push_back(std::move(r));continue;}
        const auto& rim=levels[*best];r.rim=rim.border;r.region=rim.region;
        r.area_fraction=chart_area>0?rim.area/chart_area:0.0;
        // Depth: the rim's mean height over the lowest region vertex, against the opening
        // diameter (rim length / pi).
        std::set<MeshIndex> rim_vertices;for(auto e:rim.border)for(auto v:t.edge_vertices[e])rim_vertices.insert(v);
        double rim_height=0;for(auto v:rim_vertices)rim_height+=height(in.positions[v]);rim_height/=static_cast<double>(rim_vertices.size());
        double lowest=INFINITY;
        for(auto f:rim.region)for(auto c:m.face_corners(f))lowest=std::min(lowest,height(in.positions[static_cast<std::size_t>(in.corner_vertices[c])]));
        r.depth=std::max(0.0,rim_height-lowest);r.opening_diameter=length_of(m,rim.border)/std::numbers::pi;
        r.depth_ratio=r.opening_diameter>0?r.depth/r.opening_diameter:0.0;
        // Concavity (normal consistency): the area share of the region whose faces look towards
        // the opening (into the recess's free space). Inverted or inconsistent normals lower it;
        // the cancelling walls of a deep pocket do not. (A chart never spans an inconsistently
        // oriented edge: those are always cut.)
        {
            Vec3 opening{0,0,0};double weight=0;
            for(auto e:rim.border) {const auto [a,b]=t.edge_vertices[e];const double l=g.edge_length[e];
                opening=add(opening,mul(add(in.positions[a],in.positions[b]),0.5*l));weight+=l;}
            opening=mul(opening,1/std::max(weight,1e-300));
            double good=0,all=0;
            for(auto f:rim.region) {const double a=g.face_area[f];all+=a;if(dot(g.face_normal[f],sub(opening,g.face_center[f]))>0)good+=a;}
            r.concavity=all>0?good/all:0.0;
        }
        r.confidence=std::min(1.0,r.occlusion/0.5)*std::min(1.0,r.depth_ratio/std::max(2*o.min_depth_ratio,1e-300))*r.concavity;
        // The seam slightly inside the rim: the nearest valid level below it with another
        // loop and no forbidden edge (the rim itself otherwise). The rim stays an alternative.
        const Level* seam=nullptr;
        for(std::size_t k=*best;k-->0;) {
            const auto& L=levels[k];
            if(L.valid && L.border!=rim.border && L.area>=o.min_area_fraction*chart_area && !forbidden(L.border)) {seam=&L;break;}
        }
        const bool rim_allowed=!forbidden(rim.border);
        if(!seam && rim_allowed)seam=&rim;
        if(seam) {r.seam=seam->border;r.seam_inside_rim=seam!=&rim;}
        // Interior slit of a loop: the cheapest admissible path inside the loop's region from
        // its deepest interior vertex (the farthest from the loop among equally deep ones) to
        // the loop. It is needed when that region is itself deep: its depth over the loop's
        // diameter at least slit_depth_ratio (and it has an interior vertex to start from).
        struct Slit {std::vector<MeshIndex> edges;double depth_ratio=0;bool needed=false;};
        auto slit_of=[&](const Level& L) {
            Slit out_slit;
            // GB m7: region-local vertices (sorted: local order is global order) and flags.
            std::vector<std::uint8_t> region(n,0);for(auto f:L.region)region[slot(f)]=1;
            std::vector<MeshIndex> rv;for(auto f:L.region)for(auto c:m.face_corners(f))rv.push_back(static_cast<MeshIndex>(in.corner_vertices[c]));
            std::sort(rv.begin(),rv.end());rv.erase(std::unique(rv.begin(),rv.end()),rv.end());
            auto vslot=[&](MeshIndex v) {return static_cast<std::size_t>(std::lower_bound(rv.begin(),rv.end(),v)-rv.begin());};
            std::vector<std::uint8_t> on_loop(rv.size(),0);for(auto e:L.border)for(auto v:t.edge_vertices[e])on_loop[vslot(v)]=1;
            auto usable=[&](MeshIndex e) {
                if(cuts[e] || always_cut(t,e) || (e<hc.forbidden.size() && hc.forbidden[e]))return false;
                for(auto f:t.edge_faces[e]) {const auto j=local(f);if(j<0 || !region[static_cast<std::size_t>(j)])return false;}
                return true;
            };
            // Surface distance from the loop inside the region (the tie-break of equal depths).
            std::vector<double> from_loop(rv.size(),INFINITY);using Item=std::pair<double,MeshIndex>;std::priority_queue<Item,std::vector<Item>,std::greater<>> heap;
            for(std::size_t i=0;i<rv.size();++i)if(on_loop[i]) {from_loop[i]=0;heap.push({0.0,rv[i]});}
            while(!heap.empty()) {
                const auto [d,v]=heap.top();heap.pop();if(d>from_loop[vslot(v)])continue;
                for(auto e:t.vert_edges[v]) {
                    if(!usable(e))continue;
                    const auto w=m.edge_other(e,v);const double nd=d+g.edge_length[e];const auto iw=vslot(w);
                    if(nd<from_loop[iw]) {from_loop[iw]=nd;heap.push({nd,w});}
                }
            }
            std::set<MeshIndex> interior;
            for(auto f:L.region)for(auto c:m.face_corners(f)) {
                const auto v=static_cast<MeshIndex>(in.corner_vertices[c]);const auto i=vslot(v);
                if(!on_loop[i] && std::isfinite(from_loop[i]))interior.insert(v);
            }
            std::set<MeshIndex> loop_vertices;for(auto e:L.border)for(auto v:t.edge_vertices[e])loop_vertices.insert(v);
            double loop_height=0;for(auto v:loop_vertices)loop_height+=height(in.positions[v]);loop_height/=static_cast<double>(loop_vertices.size());
            double lowest_here=INFINITY;
            for(auto f:L.region)for(auto c:m.face_corners(f))lowest_here=std::min(lowest_here,height(in.positions[static_cast<std::size_t>(in.corner_vertices[c])]));
            const double diameter=length_of(m,L.border)/std::numbers::pi;
            out_slit.depth_ratio=diameter>0?std::max(0.0,loop_height-lowest_here)/diameter:0.0;
            out_slit.needed=!interior.empty() && out_slit.depth_ratio>=o.slit_depth_ratio;
            double deepest_h=INFINITY;for(auto v:interior)deepest_h=std::min(deepest_h,height(in.positions[v]));
            const double tolerance=1e-9*std::max(r.opening_diameter,1e-300);
            std::optional<MeshIndex> deepest;
            for(auto v:interior) {
                if(height(in.positions[v])>deepest_h+tolerance)continue;
                if(!deepest || from_loop[vslot(v)]>from_loop[vslot(*deepest)])deepest=v;
            }
            if(deepest && !routing.empty())out_slit.edges=local_cheapest_path(m,rv,*deepest,on_loop,usable,routing);
            return out_slit;
        };
        if(seam) {auto sl=slit_of(*seam);r.slit=std::move(sl.edges);r.seam_depth_ratio=sl.depth_ratio;r.needs_slit=sl.needed;}
        if(r.seam_inside_rim && rim_allowed) {
            auto sl=slit_of(rim);r.rim_alternative=true;r.rim_slit=std::move(sl.edges);r.rim_needs_slit=sl.needed;
        }
        if(r.area_fraction<o.min_area_fraction)r.status="too_small";
        else if(r.depth_ratio<o.min_depth_ratio)r.status="shallow";
        else if(!seam)r.status="forbidden_edge";
        else if(r.confidence<o.min_confidence)r.status="low_confidence";
        else r.status="candidate";
        out.cavities.push_back(std::move(r));
    }
    return out;
}
ProposalGenerator cavity_generator() {
    return {name,version,[](const RunState& state,ProposalContext& ctx) {
        std::vector<Proposal> out;
        const auto& o=ctx.options.enhanced.cavities;
        if(!o.enabled || !state.surface)return out;
        const auto started=std::chrono::steady_clock::now();
        const auto& mesh=ctx.mesh;const auto& m=*mesh;const auto nf=m.geometry().face_area.size();
        const auto& ct=state.surface->tables();const auto& hc=ctx.constraints;
        const auto builds_before=triangle_bvh_builds();
        std::shared_ptr<const TriangleBvh> bvh;
        try {bvh=context_bvh(ctx);}
        catch(const ContractError& error) {
            if(error.code!=ErrorCode::memory_limit)throw;
            ctx.diagnostics[name]={{"version",version},{"skipped","memory_limit"},{"reason",error.what()},{"proposals",0}};
            return out;
        }
        const double bvh_seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count();
        const auto routing=context_routing_cost(ctx);
        const auto& scale=context_scale(ctx,state);
        // GB m7: the component diagonals once per pass; the pass's control inside the detection.
        DetectorPass work;work.control.cancel=ctx.cancel;work.control.deadline=ctx.deadline;
        Json records=Json::array();std::uint64_t id=0,proposed_cavities=0,rays=0,hits=0,charts_searched=0,charts_with_cores=0;
        std::map<std::string,std::uint64_t> statuses;std::set<EdgeMask> emitted;bool stopped=false;
        for(MeshIndex ch=0;ch<ct.n_charts && !stopped;++ch) {
            const auto& c=state.charts[ch];
            if(c.exempt_only || !c.measured)continue;
            const auto faces=ct.chart_faces[ch];
            // T39: the chart floor is an area under the physical scale (the face count in hops mode).
            double area=0;for(auto f:faces)area+=m.geometry().face_area[f];
            if(scale.floor(o.min_chart_faces,faces).below(area,faces.size()))continue;
            if(work.stop()) {stopped=true;break;}
            const auto detection=detect_cavities(m,*bvh,faces,state.cuts,hc,routing,ctx.options,true,&scale,&work);
            stopped=detection.stopped;
            ++charts_searched;rays+=detection.occlusion.rays;hits+=detection.occlusion.hits;
            for(const auto& cav:detection.cavities)++statuses[cav.status];
            if(detection.cores==0)continue;       // no occluded samples: counted, not listed
            ++charts_with_cores;
            Json record={{"chart",ch},{"faces",faces.size()},{"detection",detection.to_json()},{"proposals",Json::array()}};
            std::vector<std::uint8_t> scope;
            for(std::size_t k=0;k<detection.cavities.size() && !stopped;++k) {
                const auto& cav=detection.cavities[k];
                if(cav.status!="candidate")continue;
                if(scope.empty()) {scope.assign(nf,0);for(auto f:faces)scope[f]=1;}
                if(static_cast<std::int64_t>(proposed_cavities)>=o.max_cavities) {record["proposals"].push_back({{"cavity",k},{"status","beyond_cavity_budget"}});continue;}
                ++proposed_cavities;
                // The seam (inside the rim when a nearby loop exists), then the rim itself; each
                // loop also with its interior slit when its region needs one (deep against its
                // opening; a slit in a shallow region is only a longer seam).
                struct Variant {std::string label;const std::vector<MeshIndex>* loop;const std::vector<MeshIndex>* slit;bool needs_slit;};
                std::vector<Variant> variants{{"seam",&cav.seam,nullptr,cav.needs_slit}};
                if(cav.needs_slit)variants.push_back({"seam_slit",&cav.seam,&cav.slit,true});
                if(cav.rim_alternative) {
                    variants.push_back({"rim",&cav.rim,nullptr,cav.rim_needs_slit});
                    if(cav.rim_needs_slit)variants.push_back({"rim_slit",&cav.rim,&cav.rim_slit,true});
                }
                for(const auto& variant:variants) {
                    Json p={{"cavity",k},{"variant",variant.label}};
                    if(variant.slit && variant.slit->empty()) {p["status"]="no_slit_path";record["proposals"].push_back(p);continue;}
                    auto cuts=state.cuts;for(auto e:*variant.loop)cuts[e]=1;
                    if(variant.slit)for(auto e:*variant.slit)cuts[e]=1;
                    const auto opening=open_non_disk_charts(mesh,cuts,routing,hc,ctx.exempt,scope,work.control);
                    if(opening.stopped) {stopped=true;p["status"]="stopped";record["proposals"].push_back(p);break;}
                    p["openings"]=opening.added;
                    if(!opening.ok) {p["status"]="not_openable";p["failures"]=opening.failures;record["proposals"].push_back(p);continue;}
                    const auto violations=constraint_violations(m,hc,cuts);
                    if(!violations.empty()) {p["status"]="constraint";p["violations"]=violations;record["proposals"].push_back(p);continue;}
                    if(cuts==state.cuts || !emitted.insert(cuts).second) {p["status"]="duplicate";record["proposals"].push_back(p);continue;}
                    ++id;
                    Json provenance={{"chart",ch},{"cavity",k},{"variant",variant.label},{"loop",edge_list(*variant.loop)},
                                     {"slit",variant.slit?edge_list(*variant.slit):Json::array()},{"rim",edge_list(cav.rim)},{"openings",opening.added},
                                     {"confidence",cav.confidence},{"occlusion",cav.occlusion},{"depth_ratio",cav.depth_ratio},
                                     {"area_fraction",cav.area_fraction},{"needs_slit",variant.needs_slit},{"loop_inside_rim",variant.loop!=&cav.rim && cav.seam_inside_rim},
                                     {"detector_version",cavity_detector_version}};
                    out.push_back(cut_proposal(id,state,cuts,name,version,std::move(provenance)));
                    p["status"]="proposed";p["id"]=id;record["proposals"].push_back(p);
                }
            }
            // Reports stay bounded: the first 32 charts with occluded cores in full.
            if(records.size()<32)records.push_back(std::move(record));
        }
        const auto bvh_json=Json{{"geometry_key",bvh->geometry_key()},{"nodes",bvh->nodes()},{"triangles",bvh->triangles()},{"depth",bvh->depth()},
                                 {"memory_bytes",bvh->memory_bytes()},{"built_now",triangle_bvh_builds()>builds_before},{"seconds",bvh_seconds}};
        ctx.diagnostics[name]={{"version",version},{"detector_version",cavity_detector_version},
                               {"options",{{"max_samples",o.max_samples},{"rays_per_sample",o.rays_per_sample},{"cone_half_angle_deg",o.cone_half_angle_deg},
                                           {"reach_fraction",o.reach_fraction},{"occlusion_threshold",o.occlusion_threshold},{"level_samples",o.level_samples},
                                           {"min_area_fraction",o.min_area_fraction},{"max_area_fraction",o.max_area_fraction},
                                           {"min_depth_ratio",o.min_depth_ratio},{"slit_depth_ratio",o.slit_depth_ratio},
                                           {"min_confidence",o.min_confidence},{"max_cavities",o.max_cavities},{"min_chart_faces",o.min_chart_faces}}},
                               {"bvh",bvh_json},{"rays",rays},{"hits",hits},{"charts_searched",charts_searched},{"charts_with_cores",charts_with_cores},
                               {"cavity_status",statuses},{"stopped_by_deadline",stopped},
                               {"detector_seconds",std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count()},
                               {"charts",records},{"proposals",out.size()}};
        return out;
    }};
}
}
