// SPDX-License-Identifier: GPL-3.0-or-later
// Phase B explicit painted preferences and hard constraint composition (T29, Section 12).
// The four fields stay separate: raw exposure (measured, never modified), geometric
// feature desirability, soft artist preference (allow/avoid) and hard admissibility
// (required/forbidden/locked). Painted attributes are interpreted here, natively; the host
// only transfers them. Used by policy enhanced_v1 only.
#include "cleanseams/proposals.hpp"
#include "cleanseams/cache.hpp"
#include "cleanseams/errors.hpp"
#include "cleanseams/detail/blake2b.hpp"
#include <algorithm>
#include <cmath>
#include <cstring>
#include <numbers>
#include <queue>
#include <set>
namespace cleanseams {
namespace {
using Json=nlohmann::json;
bool always_cut(const MeshTables& t,std::size_t e) {
    return t.edge_boundary[e] || !t.edge_manifold[e] || !t.edge_consistent[e] || t.edge_corners[e].size()!=2;
}
std::vector<MeshIndex> ids(std::span<const std::uint8_t> mask) {
    std::vector<MeshIndex> out;for(std::size_t i=0;i<mask.size();++i)if(mask[i])out.push_back(static_cast<MeshIndex>(i));return out;
}
enum class Role {required,forbidden,allow,avoid,none};
struct Contract {const char* domain;const char* data_type;};
Contract contract_of(Role r) {
    switch(r) {
    case Role::required: case Role::forbidden: return {"EDGE","BOOLEAN"};
    default: return {"POINT","FLOAT"};
    }
}
const char* role_name(Role r) {
    switch(r) {case Role::required:return "required";case Role::forbidden:return "forbidden";case Role::allow:return "allow";case Role::avoid:return "avoid";default:return "none";}
}
// Clamps finite weights into [0, 1]; a NaN or infinite weight is invalid input.
std::uint64_t clamp_weights(std::vector<double>& w,const std::string& field) {
    std::uint64_t clamped=0;
    for(auto& x:w) {
        if(!std::isfinite(x))throw ContractError(ErrorCode::invalid_input,field,"preference weights must be finite (NaN or infinity refused)");
        if(x<0 || x>1) {x=std::clamp(x,0.0,1.0);++clamped;}
    }
    return clamped;
}
void check_mask(std::span<const double> values,const std::string& field) {
    for(auto x:values)if(x!=0 && x!=1)throw ContractError(ErrorCode::invalid_input,field,"boolean attribute values must be 0 or 1");
}
void add_bytes(detail::Blake2b& h,const void* data,std::size_t size) {h.update({static_cast<const std::uint8_t*>(data),size});}
void add_doubles(detail::Blake2b& h,std::span<const double> v) {std::uint64_t n=v.size();add_bytes(h,&n,sizeof n);if(n)add_bytes(h,v.data(),v.size_bytes());}
void add_mask(detail::Blake2b& h,std::span<const std::uint8_t> v) {std::uint64_t n=v.size();add_bytes(h,&n,sizeof n);if(n)add_bytes(h,v.data(),v.size());}
void add_string(detail::Blake2b& h,const std::string& s) {std::uint64_t n=s.size();add_bytes(h,&n,sizeof n);add_bytes(h,s.data(),s.size());}
}
std::string values_digest(std::span<const double> values) {
    detail::Blake2b h(32);add_doubles(h,values);return h.hex_digest();
}
std::vector<std::string> painted_attribute_names(const PipelineOptions& options) {
    const auto& p=options.enhanced.preferences;
    std::vector<std::string> out;
    // GB M2: only enhanced_v1 reads painted attributes. Under legacy_010 and native_baseline no
    // name is returned, so a host transfers nothing and its stale-state digest covers exactly
    // what it covered at Gate A.
    if(options.policy!="enhanced_v1")return out;
    for(const auto* name:{&p.required_attribute,&p.forbidden_attribute,&p.allow_attribute,&p.avoid_attribute})
        if(!name->empty() && std::find(out.begin(),out.end(),*name)==out.end())out.push_back(*name);
    return out;
}
bool ResolvedPreferences::hard_constraints() const {
    return std::any_of(required.begin(),required.end(),[](auto x) {return x!=0;}) || std::any_of(forbidden.begin(),forbidden.end(),[](auto x) {return x!=0;});
}
ResolvedPreferences resolve_preferences(const MeshSnapshot& mesh,const ConstraintInput& input,const PreferenceOptions& o,
                                        const std::optional<EdgeMask>& locked,std::span<const TransferIndex> anchors,double feature_angle_deg,
                                        std::optional<Vec3> view) {
    const auto& t=mesh.tables();const auto& g=mesh.geometry();
    const auto ne=t.edge_vertices.size(),nv=mesh.input().positions.size();
    for(const double w:{o.w_exposure,o.w_view,o.w_avoid,o.w_allow,o.w_crease})
        if(!std::isfinite(w) || w<0)throw ContractError(ErrorCode::invalid_option,"enhanced.preferences","routing weights must be finite and nonnegative");
    ResolvedPreferences r;r.required.assign(ne,0);r.forbidden.assign(ne,0);r.allow.assign(nv,0.0);r.avoid.assign(nv,0.0);
    if(locked) {
        if(locked->size()!=ne)throw ContractError(ErrorCode::invalid_input,"locked","edge mask does not match the mesh edges");
        r.locked=*locked;
    } else r.locked.assign(ne,0);
    for(auto a:anchors) {
        if(a<0 || static_cast<std::uint64_t>(a)>=nv)throw ContractError(ErrorCode::invalid_input,"anchors","anchor vertex out of range");
        r.anchors.push_back(static_cast<MeshIndex>(a));
    }
    std::sort(r.anchors.begin(),r.anchors.end());r.anchors.erase(std::unique(r.anchors.begin(),r.anchors.end()),r.anchors.end());
    Json sources=Json::array();
    // Typed inputs (already in their domains: case files, library callers).
    auto typed_mask=[&](const std::vector<std::uint8_t>& mask,EdgeMask& into,const char* field) {
        if(mask.empty())return;
        if(mask.size()!=ne)throw ContractError(ErrorCode::invalid_input,field,"edge mask does not match the mesh edges");
        std::uint64_t n=0;for(std::size_t e=0;e<ne;++e) {
            if(mask[e]>1)throw ContractError(ErrorCode::invalid_input,field,"edge mask values must be zero or one");
            if(mask[e]) {into[e]=1;++n;}
        }
        sources.push_back({{"source",field},{"status","applied"},{"count",n}});
    };
    auto typed_weights=[&](const std::vector<double>& w,std::vector<double>& into,const char* field) {
        if(w.empty())return;
        if(w.size()!=nv)throw ContractError(ErrorCode::invalid_input,field,"preference length mismatch");
        auto values=w;const auto clamped=clamp_weights(values,field);
        std::uint64_t n=0;for(std::size_t v=0;v<nv;++v) {into[v]=std::max(into[v],values[v]);n+=values[v]>0;}
        sources.push_back({{"source",field},{"status","applied"},{"count",n},{"clamped",clamped}});
    };
    typed_mask(input.required_edges,r.required,"constraints.required_edges");
    typed_mask(input.forbidden_edges,r.forbidden,"constraints.forbidden_edges");
    typed_weights(input.allow_points,r.allow,"constraints.allow_points");
    typed_weights(input.avoid_points,r.avoid,"constraints.avoid_points");
    // Painted host attributes, matched only by the configured names (documented names or
    // explicit aliases); never inferred from other names.
    std::set<std::string> seen;
    for(const auto& a:input.attributes) {
        Role role=Role::none;
        if(!o.required_attribute.empty() && a.name==o.required_attribute)role=Role::required;
        else if(!o.forbidden_attribute.empty() && a.name==o.forbidden_attribute)role=Role::forbidden;
        else if(!o.allow_attribute.empty() && a.name==o.allow_attribute)role=Role::allow;
        else if(!o.avoid_attribute.empty() && a.name==o.avoid_attribute)role=Role::avoid;
        Json record={{"source","attribute"},{"name",a.name},{"domain",a.domain},{"data_type",a.data_type},{"role",role_name(role)}};
        if(role==Role::none) {record["status"]="ignored_unconfigured_name";sources.push_back(record);continue;}
        if(!seen.insert(a.name).second) {record["status"]="ignored_duplicate";sources.push_back(record);continue;}
        const auto c=contract_of(role);
        record["expected"]={{"domain",c.domain},{"data_type",c.data_type}};
        // Wrong domains/types are reported and ignored, never guessed into another domain.
        if(a.domain!=c.domain) {record["status"]="ignored_wrong_domain";sources.push_back(record);continue;}
        if(a.data_type!=c.data_type) {record["status"]="ignored_wrong_type";sources.push_back(record);continue;}
        const std::string field="constraints.attributes."+a.name;
        const auto count=std::string(c.domain)=="EDGE"?ne:nv;
        if(a.values.size()!=count)throw ContractError(ErrorCode::invalid_input,field,"attribute length does not match its domain");
        std::uint64_t n=0,clamped=0;
        if(role==Role::required || role==Role::forbidden) {
            check_mask(a.values,field);auto& into=role==Role::required?r.required:r.forbidden;
            for(std::size_t e=0;e<ne;++e)if(a.values[e]!=0) {into[e]=1;++n;}
        } else {
            auto values=a.values;clamped=clamp_weights(values,field);auto& into=role==Role::allow?r.allow:r.avoid;
            for(std::size_t v=0;v<nv;++v) {into[v]=std::max(into[v],values[v]);n+=values[v]>0;}
        }
        record["status"]="applied";record["count"]=n;record["clamped"]=clamped;sources.push_back(record);
    }
    // Hard rules: physical boundaries are topology facts (a forbidden flag there is
    // inapplicable, a required flag already met); required/locked against forbidden is a
    // conflict refused with the exact IDs, never settled by priority.
    std::vector<MeshIndex> inapplicable,already_cut,required_and_forbidden,locked_and_forbidden;
    for(std::size_t e=0;e<ne;++e) {
        const auto id=static_cast<MeshIndex>(e);
        if(always_cut(t,e)) {
            if(r.forbidden[e]) {inapplicable.push_back(id);r.forbidden[e]=0;}
            if(r.required[e])already_cut.push_back(id);
            continue;
        }
        if(r.required[e] && r.forbidden[e])required_and_forbidden.push_back(id);
        if(r.locked[e] && r.forbidden[e])locked_and_forbidden.push_back(id);
    }
    r.inapplicable={{"forbidden_on_physical_boundary",inapplicable},{"required_on_physical_boundary",already_cut}};
    r.conflicts=Json::object();
    if(!required_and_forbidden.empty())r.conflicts["required_and_forbidden"]=required_and_forbidden;
    if(!locked_and_forbidden.empty())r.conflicts["locked_and_forbidden"]=locked_and_forbidden;
    r.sources=std::move(sources);
    // Geometric feature desirability: sharpness relative to the feature angle, in [0, 1].
    r.feature_desirability.assign(ne,0.0);
    if(feature_angle_deg>0)
        for(std::size_t e=0;e<ne;++e)r.feature_desirability[e]=std::clamp(std::abs(g.edge_dihedral[e])/(feature_angle_deg*std::numbers::pi/180.0),0.0,1.0);
    // Soft preference per edge: the mean of its end points.
    r.allow_edge.resize(ne);r.avoid_edge.resize(ne);
    for(std::size_t e=0;e<ne;++e) {
        const auto [a,b]=t.edge_vertices[e];
        r.allow_edge[e]=0.5*(r.allow[a]+r.allow[b]);r.avoid_edge[e]=0.5*(r.avoid[a]+r.avoid[b]);
    }
    // Explicit viewer only: no implicit front. Facing as the legacy seam cost measures it.
    if(view) {
        const auto d=view_direction(*view);r.view_direction=d;r.facing.resize(ne);
        for(std::size_t e=0;e<ne;++e) {
            const auto& na=g.vert_normal[t.edge_vertices[e][0]];const auto& nb=g.vert_normal[t.edge_vertices[e][1]];
            const Vec3 n{0.5*(na[0]+nb[0]),0.5*(na[1]+nb[1]),0.5*(na[2]+nb[2])};
            r.facing[e]=std::clamp(n[0]*d[0]+n[1]*d[1]+n[2]*d[2],-1.0,1.0);
        }
    }
    // Content identity of everything painted and configured: proposals evaluated under
    // another digest are stale; geometry and topology caches never depend on it.
    detail::Blake2b h(32);add_string(h,"preferences.v1");add_mask(h,r.required);add_mask(h,r.forbidden);add_mask(h,r.locked);
    add_doubles(h,r.allow);add_doubles(h,r.avoid);
    for(const double w:{o.w_exposure,o.w_view,o.w_avoid,o.w_allow,o.w_crease,feature_angle_deg})add_bytes(h,&w,sizeof w);
    if(view)add_bytes(h,r.view_direction->data(),sizeof(Vec3));
    std::vector<std::uint64_t> anchor_ids(r.anchors.begin(),r.anchors.end());std::uint64_t na=anchor_ids.size();add_bytes(h,&na,sizeof na);
    if(na)add_bytes(h,anchor_ids.data(),na*sizeof(std::uint64_t));
    r.digest=h.hex_digest();
    return r;
}
std::vector<double> enhanced_edge_cost(const MeshSnapshot& mesh,const ResolvedPreferences& p,const PreferenceOptions& o) {
    const auto& g=mesh.geometry();const auto ne=g.edge_length.size();
    if(p.allow_edge.size()!=ne || p.feature_desirability.size()!=ne)throw ContractError(ErrorCode::invalid_input,"preferences","preferences belong to another mesh");
    if(!p.raw_exposure.empty() && p.raw_exposure.size()!=ne)throw ContractError(ErrorCode::invalid_input,"preferences.raw_exposure","exposure needs one value per edge");
    std::vector<double> cost(ne);
    for(std::size_t e=0;e<ne;++e) {
        const double x=p.raw_exposure.empty()?0.0:p.raw_exposure[e];
        const double seen=p.facing.empty()?0.0:std::max(0.0,p.facing[e]);
        const double numerator=1.0+o.w_exposure*x*x+o.w_view*seen+o.w_avoid*p.avoid_edge[e];
        const double denominator=1.0+o.w_allow*p.allow_edge[e]+o.w_crease*p.feature_desirability[e]; // >= 1
        cost[e]=g.edge_length[e]*numerator/denominator;
    }
    return cost;
}
nlohmann::json ResolvedPreferences::explain(const MeshSnapshot& mesh,MeshIndex edge,const PreferenceOptions& o) const {
    const auto& t=mesh.tables();const auto& g=mesh.geometry();
    if(edge>=t.edge_vertices.size())throw ContractError(ErrorCode::invalid_input,"preferences.edge","edge out of range");
    const double x=raw_exposure.empty()?0.0:raw_exposure[edge];
    const double seen=facing.empty()?0.0:std::max(0.0,facing[edge]);
    const double numerator=1.0+o.w_exposure*x*x+o.w_view*seen+o.w_avoid*avoid_edge[edge];
    const double denominator=1.0+o.w_allow*allow_edge[edge]+o.w_crease*feature_desirability[edge];
    std::string hard="free";
    if(always_cut(t,edge))hard="physical_boundary";
    else if(forbidden[edge])hard="forbidden";
    else if(required[edge] || locked[edge])hard=required[edge]?"required":"locked";
    return {{"edge",edge},{"vertices",{t.edge_vertices[edge][0],t.edge_vertices[edge][1]}},{"length",g.edge_length[edge]},
            {"raw_exposure",raw_exposure.empty()?Json(nullptr):Json(x)},{"facing",facing.empty()?Json(nullptr):Json(facing[edge])},
            {"allow",allow_edge[edge]},{"avoid",avoid_edge[edge]},{"feature_desirability",feature_desirability[edge]},
            {"terms",{{"exposure",o.w_exposure*x*x},{"view",o.w_view*seen},{"avoid",o.w_avoid*avoid_edge[edge]},
                      {"allow",o.w_allow*allow_edge[edge]},{"crease",o.w_crease*feature_desirability[edge]}}},
            {"cost",g.edge_length[edge]*numerator/denominator},{"hard",hard}};
}
nlohmann::json ResolvedPreferences::summary() const {
    std::uint64_t allowed=0,avoided=0;for(auto x:allow)allowed+=x>0;for(auto x:avoid)avoided+=x>0;
    return {{"version","preferences.v1"},{"digest",digest},{"sources",sources},{"required_edges",ids(required)},{"forbidden_edges",ids(forbidden)},
            {"allow_points",allowed},{"avoid_points",avoided},{"inapplicable",inapplicable},{"conflicts",conflicts},
            {"view_direction",view_direction?Json(*view_direction):Json(nullptr)},{"raw_exposure_measured",!raw_exposure.empty()}};
}
// ------------------------------------------------------------ constraint repair
namespace {
struct Search {
    const MeshSnapshot& mesh;
    std::span<const double> cost;
    // Cheapest path over mesh edges from any source to any target: `usable(e)` filters
    // edges, `price(e)` is the added cost. Deterministic: (distance, vertex) order.
    template<class Usable,class Price>
    std::vector<MeshIndex> path(std::span<const MeshIndex> sources,const std::vector<std::uint8_t>& target,Usable usable,Price price) const {
        const auto& t=mesh.tables();const auto nv=mesh.input().positions.size();
        std::vector<double> dist(nv,std::numeric_limits<double>::infinity());std::vector<std::int64_t> via(nv,-1);
        using Item=std::pair<double,MeshIndex>;std::priority_queue<Item,std::vector<Item>,std::greater<>> heap;
        for(auto s:sources) {if(dist[s]>0) {dist[s]=0;heap.push({0.0,s});}}
        while(!heap.empty()) {
            const auto [d,v]=heap.top();heap.pop();
            if(d>dist[v])continue;
            if(target[v]) {
                std::vector<MeshIndex> edges;auto at=v;
                while(via[at]>=0) {const auto e=static_cast<MeshIndex>(via[at]);edges.push_back(e);at=mesh.edge_other(e,at);}
                return edges;
            }
            for(auto e:t.vert_edges[v]) {
                if(!usable(e))continue;
                const auto w=mesh.edge_other(e,v);const double nd=d+price(e);
                if(nd<dist[w]) {dist[w]=nd;via[w]=e;heap.push({nd,w});}
            }
        }
        return {};
    }
};
// Opens every measured non-disk chart it can: two or more boundary loops are joined by the
// cheapest admissible interior path; a closed chart gets its cheapest admissible edge.
// GB m7: every round rebuilds the whole surface; `control` is checked first (cancellation
// throws; false when the deadline stopped it, the cut set then left incomplete).
bool open_charts(const std::shared_ptr<const MeshSnapshot>& mesh,EdgeMask& cuts,const ResolvedPreferences& p,std::span<const double> cost,
                 std::span<const std::uint8_t> exempt,const Control& control) {
    const auto& t=mesh->tables();const auto& g=mesh->geometry();const Search search{*mesh,cost};
    for(int round=0;round<256;++round) {
        if(control.check())return false;
        const auto surface=cut_surface(mesh,cuts);const auto& ct=surface->tables();bool changed=false;
        for(MeshIndex ch=0;ch<ct.n_charts && !changed;++ch) {
            if(ct.chart_is_disk[ch] || ct.chart_pinched[ch] || ct.chart_unsupported[ch])continue;
            bool measured=false;for(auto f:ct.chart_faces[ch])measured=measured || (!(f<exempt.size() && exempt[f]) && g.face_degenerate_kind[f]==0);
            if(!measured)continue;
            auto interior=[&](MeshIndex e) {
                if(cuts[e] || always_cut(t,e) || p.forbidden[e])return false;
                for(auto f:t.edge_faces[e])if(ct.chart_of_face[f]!=ch)return false;
                return true;
            };
            const auto loops=ct.chart_loops[ch];
            if(loops.size()>=2) {
                const auto start=surface->loop_verts(loops[0]);std::vector<std::uint8_t> target(mesh->input().positions.size(),0);
                for(std::size_t k=1;k<loops.size();++k)for(auto v:surface->loop_verts(loops[k]))target[v]=1;
                for(auto v:start)target[v]=0;
                const auto edges=search.path(start,target,interior,[&](MeshIndex e) {return cost[e];});
                for(auto e:edges) {cuts[e]=1;changed=true;}
            } else if(loops.empty()) {
                std::optional<MeshIndex> best;
                for(auto f:ct.chart_faces[ch])for(auto c:mesh->face_corners(f)) {
                    const auto e=t.corner_edge[c];if(interior(e) && (!best || cost[e]<cost[*best] || (cost[e]==cost[*best] && e<*best)))best=e;
                }
                if(best) {cuts[*best]=1;changed=true;}
            }
        }
        if(!changed)return true;
    }
    return true;
}
Proposal delta(std::uint64_t id,const RunState& s,const EdgeMask& cuts,const char* strategy) {
    Proposal p;p.identity.id=id;p.identity.parent_revision=s.revision;p.identity.generator="constraint_repair";p.identity.generator_version="1";
    for(std::size_t e=0;e<cuts.size();++e) {
        if(cuts[e] && !s.cuts[e])p.identity.added_edges.push_back(static_cast<std::int64_t>(e));
        if(!cuts[e] && s.cuts[e])p.identity.removed_edges.push_back(static_cast<std::int64_t>(e));
    }
    p.provenance={{"strategy",strategy}};
    return p;
}
}
ProposalGenerator constraint_repair_generator(std::shared_ptr<const ResolvedPreferences> preferences,PreferenceOptions options) {
    return {"constraint_repair","1",[preferences,options](const RunState& state,ProposalContext& ctx) {
        std::vector<Proposal> out;
        const auto& mesh=ctx.mesh;const auto& t=mesh->tables();const auto ne=t.edge_vertices.size();const auto& p=*preferences;
        const auto violations=constraint_violations(*mesh,ctx.constraints,state.cuts);
        bool open_needed=false;for(const auto& o:state.charts)open_needed=open_needed || (!o.exempt_only && !o.disk);
        if(violations.empty() && !open_needed)return out;
        const auto cost=enhanced_edge_cost(*mesh,p,options);const Search search{*mesh,cost};
        EdgeMask base=state.cuts;for(std::size_t e=0;e<ne;++e)if(p.required[e] || ctx.constraints.locked[e])base[e]=1;
        std::vector<MeshIndex> forbidden_cut;for(std::size_t e=0;e<ne;++e)if(base[e] && p.forbidden[e] && !always_cut(t,e))forbidden_cut.push_back(static_cast<MeshIndex>(e));
        // Reroute: each group of forbidden cut edges is replaced by the cheapest admissible
        // connection of the vertices where it met the rest of the seam network (existing
        // cuts and physical boundaries cost nothing to follow).
        EdgeMask reroute=base;for(auto e:forbidden_cut)reroute[e]=0;
        if(!forbidden_cut.empty()) {
            std::vector<std::uint8_t> in_group(mesh->input().positions.size(),0),terminal(mesh->input().positions.size(),0);
            std::vector<MeshIndex> terminals;
            for(auto e:forbidden_cut)for(auto v:t.edge_vertices[e]) {
                in_group[v]=1;
                bool touches=false;for(auto f:t.vert_edges[v])touches=touches || reroute[f];
                if(touches && !terminal[v]) {terminal[v]=1;terminals.push_back(v);}
            }
            std::sort(terminals.begin(),terminals.end());
            if(terminals.size()>=2) {
                auto usable=[&](MeshIndex e) {return !p.forbidden[e] || always_cut(t,e);};
                auto price=[&](MeshIndex e) {return reroute[e]?0.0:cost[e];};
                std::vector<MeshIndex> tree{terminals[0]};std::vector<std::uint8_t> remaining(terminal);remaining[terminals[0]]=0;
                for(std::size_t k=1;k<terminals.size();++k) {
                    const auto edges=search.path(tree,remaining,usable,price);
                    if(edges.empty())break;
                    for(auto e:edges) {reroute[e]=1;for(auto v:t.edge_vertices[e]) {if(remaining[v])remaining[v]=0;tree.push_back(v);}}
                }
            }
        }
        // GB m7: a deadline during an opening leaves that repair incomplete: it is not proposed.
        const Control control{ctx.cancel,ctx.deadline};auto& tally=ctx.counters.by_generator["constraint_repair"];
        if(!open_charts(mesh,reroute,p,cost,ctx.exempt,control)) {++tally["deadline"];return out;}
        out.push_back(delta(1,state,reroute,"reroute"));
        // Weld: the forbidden cut edges are closed and the charts that stop being disks are
        // reopened along admissible paths.
        EdgeMask weld=base;for(auto e:forbidden_cut)weld[e]=0;
        if(!open_charts(mesh,weld,p,cost,ctx.exempt,control)) {++tally["deadline"];return out;}
        if(weld!=reroute)out.push_back(delta(2,state,weld,"weld_and_open"));
        return out;
    }};
}
}
