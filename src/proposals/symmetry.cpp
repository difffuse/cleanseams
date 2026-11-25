// SPDX-License-Identifier: GPL-3.0-or-later
// Phase B operational symmetry (Section 19, T38): the strict mirror correspondence of a mesh,
// atomic pairing of proposals with their mirror images (the transaction calls pair_proposal
// before anything is solved and refuses to commit a half pair), the Section 21 symmetrization
// stage and metric-only diagnostics. Used by policy enhanced_v1 only; the policy is the explicit
// option enhanced.symmetry, never an environment variable.
#include "cleanseams/symmetric_moves.hpp"
#include "cleanseams/cache.hpp"
#include "cleanseams/errors.hpp"
#include "cleanseams/layout.hpp"
#include <algorithm>
#include <cmath>
#include <map>
#include <set>
namespace cleanseams {
namespace {
using Json=nlohmann::json;
constexpr std::uint64_t symmetry_id_base=38'000'000;
constexpr const char* generator_name="symmetry";
constexpr const char* generator_version="1";
bool topology_fact(const MeshTables& t,std::size_t e) {
    return t.edge_boundary[e] || !t.edge_manifold[e] || !t.edge_consistent[e] || t.edge_corners[e].size()!=2;
}
std::size_t face_total(const MeshSnapshot& m) {return m.input().face_offsets.empty()?0:m.input().face_offsets.size()-1;}
double distance(const Vec3& a,const Vec3& b) {
    const double x=a[0]-b[0],y=a[1]-b[1],z=a[2]-b[2];return std::sqrt(x*x+y*y+z*z);
}
Json finite_or_null(double x) {return std::isfinite(x)?Json(x):Json(nullptr);}
Json plane_json(const std::optional<MirrorPlane>& p) {
    if(!p)return nullptr;
    return {{"normal",p->normal},{"offset",p->offset},{"source",p->source},{"score",p->score},{"far",p->far}};
}
template<class T> std::vector<T> sorted(std::set<T> s) {return {s.begin(),s.end()};}
std::vector<std::int64_t> as_ids(const std::vector<MeshIndex>& v) {return {v.begin(),v.end()};}
bool same_map_options(const SymmetryOptions& a,const SymmetryOptions& b) {
    return a.plane_normal==b.plane_normal && a.plane_offset==b.plane_offset && a.map_tolerance==b.map_tolerance &&
           a.geometry_tolerance==b.geometry_tolerance;
}
// Side of an edge relative to the plane: -1 or +1 when both ends are on that side (or one is on
// the plane), 0 when both lie on the plane or the edge crosses it.
int edge_side(const StrictMirrorMap& map,const MeshTables& t,MeshIndex e) {
    const int s=map.vertex_side[t.edge_vertices[e][0]]+map.vertex_side[t.edge_vertices[e][1]];
    return s<0?-1:(s>0?1:0);
}
// A chart is measured when it holds a non-exempt, nondegenerate face (the transaction's rule).
bool measured_chart(const MeshSnapshot& m,const CutSurface& cs,MeshIndex ch,std::span<const std::uint8_t> exempt) {
    for(auto f:cs.tables().chart_faces[ch])if(!(f<exempt.size() && exempt[f]) && m.geometry().face_degenerate_kind[f]==0)return true;
    return false;
}
const Json* symmetry_record(const Proposal& p) {
    if(!p.provenance.is_object())return nullptr;
    const auto found=p.provenance.find("symmetry");
    return found!=p.provenance.end() && found->is_object()?&*found:nullptr;
}
void set_symmetry_record(Proposal& p,Json record) {
    if(!p.provenance.is_object())p.provenance=Json{{"generator_provenance",p.provenance}};
    p.provenance["symmetry"]=std::move(record);
}
}
// ------------------------------------------------------------ the strict correspondence
std::optional<MeshIndex> StrictMirrorMap::mirror_edge(MeshIndex e) const {
    if(e>=edge_exact.size() || !edge_exact[e])return std::nullopt;
    return static_cast<MeshIndex>(maps.edge[e]);
}
std::optional<MeshIndex> StrictMirrorMap::mirror_vertex(MeshIndex v) const {
    if(v>=vertex_exact.size() || !vertex_exact[v])return std::nullopt;
    return static_cast<MeshIndex>(maps.vertex[v]);
}
nlohmann::json StrictMirrorMap::to_json() const {
    return {{"version",symmetric_moves_version},{"status",status},{"plane",plane_json(plane)},{"plane_source",plane_source},
            {"diagonal",diagonal},{"map_tolerance",options.map_tolerance},{"geometry_tolerance",options.geometry_tolerance},
            {"map_radius",map_radius},{"tolerance",tolerance},{"max_deviation",max_deviation},
            {"vertices",vertex_exact.size()},{"edges",edge_exact.size()},{"faces",face_exact.size()},
            {"matched_vertices",matched_vertices},{"exact_vertices",exact_vertices},{"exact_edges",exact_edges},{"exact_faces",exact_faces},
            {"self_mirrored_edges",self_mirrored_edges}};
}
StrictMirrorMap strict_mirror_map(const MeshSnapshot& m,const SymmetryOptions& o) {
    StrictMirrorMap out;out.options=o;out.geometry_key=m.geometry_key();
    const auto& t=m.tables();const auto& pos=m.input().positions;
    const auto nv=pos.size(),ne=t.edge_vertices.size(),nf=face_total(m);
    out.vertex_exact.assign(nv,0);out.edge_exact.assign(ne,0);out.face_exact.assign(nf,0);out.vertex_side.assign(nv,0);
    out.diagonal=bbox_diagonal(m);out.map_radius=o.map_tolerance*out.diagonal;out.tolerance=o.geometry_tolerance*out.diagonal;
    if(o.plane_normal) {
        const auto& n=*o.plane_normal;const double length=std::sqrt(n[0]*n[0]+n[1]*n[1]+n[2]*n[2]);
        if(!(length>0) || !std::isfinite(length) || !std::isfinite(o.plane_offset))
            throw ContractError(ErrorCode::invalid_option,"enhanced.symmetry.plane_normal","the plane needs a finite nonzero normal and a finite offset");
        MirrorPlane p;p.normal={n[0]/length,n[1]/length,n[2]/length};p.offset=o.plane_offset/length;p.source="explicit";
        out.plane=p;out.plane_source="explicit";
    } else if(auto detected=detect_plane(m)) {out.plane=*detected;out.plane_source="detected";}
    else out.plane_source="none";
    if(!out.plane) {
        out.status="no_plane";out.maps.vertex.assign(nv,-1);out.maps.edge.assign(ne,-1);out.maps.face.assign(nf,-1);
        return out;
    }
    const auto& plane=*out.plane;
    out.maps=mirror_maps(m,plane,o.map_tolerance);
    for(std::size_t v=0;v<nv;++v) {
        const double d=plane.signed_distance(pos[v]);
        out.vertex_side[v]=static_cast<std::int8_t>(std::abs(d)<=out.tolerance?0:(d<0?-1:1));
        const auto w=out.maps.vertex[v];
        if(w<0)continue;
        ++out.matched_vertices;
        const double deviation=distance(plane.reflect(pos[v]),pos[static_cast<std::size_t>(w)]);
        out.max_deviation=std::max(out.max_deviation,deviation);
        if(out.maps.vertex[static_cast<std::size_t>(w)]==static_cast<TransferIndex>(v) && deviation<=out.tolerance) {out.vertex_exact[v]=1;++out.exact_vertices;}
    }
    for(std::size_t e=0;e<ne;++e) {
        const auto a=t.edge_vertices[e][0],b=t.edge_vertices[e][1];const auto image=out.maps.edge[e];
        if(!out.vertex_exact[a] || !out.vertex_exact[b] || image<0)continue;
        const auto img=static_cast<std::size_t>(image);
        if(out.maps.edge[img]!=static_cast<TransferIndex>(e) || topology_fact(t,e)!=topology_fact(t,img))continue;
        out.edge_exact[e]=1;++out.exact_edges;out.self_mirrored_edges+=img==e;
    }
    for(std::size_t f=0;f<nf;++f) {
        const auto image=out.maps.face[f];if(image<0 || out.maps.face[static_cast<std::size_t>(image)]!=static_cast<TransferIndex>(f))continue;
        bool exact=true;
        for(auto c:m.face_corners(static_cast<MeshIndex>(f)))
            exact=exact && out.vertex_exact[static_cast<std::size_t>(m.input().corner_vertices[c])] && out.edge_exact[t.corner_edge[c]];
        if(exact) {out.face_exact[f]=1;++out.exact_faces;}
    }
    if(out.exact_faces==nf && out.exact_edges==ne && out.exact_vertices==nv)out.status="exact";
    else if(out.exact_faces>0)out.status="partial";
    else if(is_combinatorial(out.maps))out.status="approximate";
    else out.status="noncombinatorial";
    return out;
}
// ------------------------------------------------------------ pairing
nlohmann::json SymmetricPairing::to_json() const {
    Json out={{"status",status},{"detail",detail}};
    if(!reason.empty())out["reason"]=reason;
    return out;
}
SymmetricPairing pair_proposal(const Proposal& p,const RunState& parent,const HardConstraints& hc,const StrictMirrorMap& map) {
    SymmetricPairing out;out.proposal=p;
    auto conflict=[&](std::string reason,Json detail) {
        out.status="conflict";out.reason=std::move(reason);out.detail=std::move(detail);out.detail["reason"]=out.reason;out.detail["half"]="mirror";
        return out;
    };
    if(!parent.surface)throw ContractError(ErrorCode::invalid_input,"symmetry.parent","pairing needs a measured parent state");
    const auto& m=parent.surface->mesh();const auto& t=m.tables();const auto& ct=parent.surface->tables();const auto ne=t.edge_vertices.size();
    if(map.geometry_key!=m.geometry_key() || map.edge_exact.size()!=ne)
        throw ContractError(ErrorCode::stale_state,"symmetry.map","the mirror map belongs to another geometry");
    if(!map.usable())
        return conflict(map.status=="approximate"?"approximate_correspondence":map.status=="noncombinatorial"?"noncombinatorial":"no_strict_map",
                        {{"map_status",map.status}});
    // The original half must be a representable delta; otherwise the transaction names it.
    std::vector<MeshIndex> A,R;std::set<std::int64_t> seen;
    for(const auto& [list,add]:{std::pair(&p.identity.added_edges,true),std::pair(&p.identity.removed_edges,false)})
        for(auto x:*list) {
            if(x<0 || static_cast<std::uint64_t>(x)>=ne || !seen.insert(x).second || topology_fact(t,static_cast<std::size_t>(x)) ||
               (parent.cuts[static_cast<std::size_t>(x)]!=0)==add)
                return conflict("unrepresentable_edge",{{"edges",Json::array({x})}});
            (add?A:R).push_back(static_cast<MeshIndex>(x));
        }
    const auto nv=m.input().positions.size();
    for(const auto& move:p.endpoints)
        if(move.from>=nv || move.to>=nv)return conflict("unrepresentable_edge",{{"endpoints",Json::array({move.from,move.to})}});
    // Exact correspondence on the affected patch: the edges and endpoints of the original half,
    // then every face of the parent charts the delta or its image touches. A failing element is
    // approximate when it has an involutive counterpart within the counterpart tolerance only,
    // noncombinatorial otherwise (no counterpart, or a triangulation without a mirror image).
    auto vertex_failure=[&](MeshIndex v) {
        if(map.vertex_exact[v])return 0;
        const auto w=map.maps.vertex[v];
        return w>=0 && map.maps.vertex[static_cast<std::size_t>(w)]==static_cast<TransferIndex>(v)?1:2;
    };
    auto edge_failure=[&](MeshIndex e) {
        if(map.edge_exact[e])return 0;
        const auto img=map.maps.edge[e];
        int worst=img>=0 && map.maps.edge[static_cast<std::size_t>(img)]==static_cast<TransferIndex>(e) &&
                  topology_fact(t,e)==topology_fact(t,static_cast<std::size_t>(img))?1:2;
        return std::max({worst,vertex_failure(t.edge_vertices[e][0]),vertex_failure(t.edge_vertices[e][1])});
    };
    auto face_failure=[&](MeshIndex f) {
        if(map.face_exact[f])return 0;
        const auto img=map.maps.face[f];
        int worst=img>=0 && map.maps.face[static_cast<std::size_t>(img)]==static_cast<TransferIndex>(f)?1:2;
        for(auto c:m.face_corners(f))
            worst=std::max({worst,vertex_failure(static_cast<MeshIndex>(m.input().corner_vertices[c])),edge_failure(t.corner_edge[c])});
        return worst;
    };
    int failure=0;std::vector<MeshIndex> failing_edges,failing_vertices,failing_faces;std::uint64_t non_exact_faces=0;
    for(const auto* list:{&A,&R})for(auto e:*list)if(const int f=edge_failure(e)) {failure=std::max(failure,f);failing_edges.push_back(e);}
    for(const auto& move:p.endpoints)for(auto v:{move.from,move.to})if(const int f=vertex_failure(v)) {failure=std::max(failure,f);failing_vertices.push_back(v);}
    auto correspondence_conflict=[&]() {
        return conflict(failure==1?"approximate_correspondence":"noncombinatorial",
                        {{"edges",failing_edges},{"vertices",failing_vertices},{"faces",failing_faces},{"non_exact_faces",non_exact_faces},
                         {"max_deviation",map.max_deviation},{"tolerance",map.tolerance},{"map_status",map.status}});
    };
    if(failure)return correspondence_conflict();
    auto image=[&](MeshIndex e) {return *map.mirror_edge(e);};
    auto vimage=[&](MeshIndex v) {return *map.mirror_vertex(v);};
    std::set<MeshIndex> charts;
    auto touch_edge=[&](MeshIndex e) {for(auto f:t.edge_faces[e])charts.insert(ct.chart_of_face[f]);};
    auto touch_vertex=[&](MeshIndex v) {for(auto f:t.vert_faces[v])charts.insert(ct.chart_of_face[f]);};
    for(auto e:A) {touch_edge(e);touch_edge(image(e));}
    for(auto e:R) {touch_edge(e);touch_edge(image(e));}
    for(const auto& move:p.endpoints) {touch_vertex(move.from);touch_vertex(move.to);touch_vertex(vimage(move.from));touch_vertex(vimage(move.to));}
    std::uint64_t patch_faces=0;
    for(auto ch:charts)for(auto f:ct.chart_faces[ch]) {
        ++patch_faces;
        if(const int fail=face_failure(f)) {failure=std::max(failure,fail);++non_exact_faces;if(failing_faces.size()<32)failing_faces.push_back(f);}
    }
    if(failure)return correspondence_conflict();
    // Merge the halves: every edge once (self-mirrored and overlapping edges included).
    std::set<MeshIndex> add(A.begin(),A.end()),rem(R.begin(),R.end());
    std::uint64_t deduplicated=0,self_mirrored=0;
    for(auto e:A) {const auto i=image(e);deduplicated+=add.contains(i);self_mirrored+=i==e;}
    for(auto e:R) {const auto i=image(e);deduplicated+=rem.contains(i);self_mirrored+=i==e;}
    for(auto e:A)add.insert(image(e));
    for(auto e:R)rem.insert(image(e));
    std::vector<MeshIndex> overlap;for(auto e:add)if(rem.contains(e))overlap.push_back(e);
    if(!overlap.empty())return conflict("overlapping_halves",{{"edges",overlap}});
    std::vector<MeshIndex> merged_add,merged_rem,mirror_added,mirror_removed;std::uint64_t unchanged=0;
    const std::set<MeshIndex> a0(A.begin(),A.end()),r0(R.begin(),R.end());
    for(auto e:add) {if(parent.cuts[e]) {++unchanged;continue;}merged_add.push_back(e);if(!a0.contains(e))mirror_added.push_back(e);}
    for(auto e:rem) {if(!parent.cuts[e]) {++unchanged;continue;}merged_rem.push_back(e);if(!r0.contains(e))mirror_removed.push_back(e);}
    // Endpoint moves and the allowed endpoint set with their images.
    std::vector<EndpointChange> moves=p.endpoints,mirror_moves;
    auto has_move=[&](const std::vector<EndpointChange>& list,MeshIndex from,MeshIndex to) {
        return std::any_of(list.begin(),list.end(),[&](const EndpointChange& c) {return c.from==from && c.to==to;});
    };
    for(const auto& move:p.endpoints) {
        const auto f=vimage(move.from),to=vimage(move.to);
        if(!has_move(moves,f,to)) {moves.push_back({f,to});mirror_moves.push_back({f,to});}
    }
    std::optional<std::vector<MeshIndex>> allowed=p.allowed_endpoints;
    if(allowed) {
        std::set<MeshIndex> all(allowed->begin(),allowed->end());
        for(auto v:*allowed)if(v<nv && map.vertex_exact[v])all.insert(vimage(v));
        allowed=sorted(all);
    }
    // Conflicts of the mirror half, before anything is solved: hard constraints the paired cut
    // set violates while the original half alone does not, anchors the mirror half moves, and
    // pins the mirror half's gluing joins with different coordinates. Nothing is mirrored.
    EdgeMask original=parent.cuts,paired=parent.cuts;
    for(auto e:A)original[e]=1;
    for(auto e:R)original[e]=0;
    for(auto e:merged_add)paired[e]=1;
    for(auto e:merged_rem)paired[e]=0;
    const auto before=constraint_violations(m,hc,original),after=constraint_violations(m,hc,paired);
    for(const auto& [key,reason]:{std::pair("locked_edges","asymmetric_lock"),std::pair("required_edges","asymmetric_required"),
                                  std::pair("forbidden_edges","asymmetric_forbidden"),std::pair("anchors","asymmetric_anchor")}) {
        if(!after.contains(key))continue;
        std::set<MeshIndex> known;if(before.contains(key))for(const auto& x:before.at(key))known.insert(x.get<MeshIndex>());
        std::vector<MeshIndex> added;for(const auto& x:after.at(key))if(!known.contains(x.get<MeshIndex>()))added.push_back(x.get<MeshIndex>());
        if(!added.empty())return conflict(reason,{{std::string(key)=="anchors"?"anchors":"edges",added},{"constraint",key}});
    }
    for(const auto& move:mirror_moves)
        if(std::binary_search(hc.anchors.begin(),hc.anchors.end(),move.from))return conflict("asymmetric_anchor",{{"anchors",Json::array({move.from})},{"constraint","anchors"}});
    // Pins are never mirrored: a pair that changes the cuts at a vertex pinned on one side only
    // would change that side's pinned sectors without a counterpart (whichever half is the
    // original), and a mirror half that glues differently pinned sectors cannot hold them.
    std::set<MeshIndex> pinned;
    for(const auto& [corner,uv]:hc.pins)if(corner<m.input().corner_vertices.size())pinned.insert(static_cast<MeshIndex>(m.input().corner_vertices[corner]));
    std::set<MeshIndex> one_sided;
    for(const auto* list:{&merged_add,&merged_rem})for(auto e:*list)for(auto v:t.edge_vertices[e])if(pinned.contains(v) && !pinned.contains(vimage(v)))one_sided.insert(v);
    if(!one_sided.empty())return conflict("asymmetric_pin",{{"vertices",sorted(one_sided)},{"constraint","pins"}});
    if(!mirror_removed.empty() && pins_conflict_after_gluing(*parent.surface,hc.pins,merged_rem) && !pins_conflict_after_gluing(*parent.surface,hc.pins,R))
        return conflict("asymmetric_pin",{{"edges",mirror_removed},{"constraint","pins"}});
    // The paired proposal: one delta over every affected chart.
    auto& q=out.proposal;
    const bool closed_delta=mirror_added.empty() && mirror_removed.empty() && mirror_moves.empty();
    bool self=closed_delta;
    for(auto e:A)self=self && a0.contains(image(e));
    for(auto e:R)self=self && r0.contains(image(e));
    out.status=self?"self_mirrored":closed_delta?"mirror_complete":"paired";
    // Idempotent: a proposal paired before keeps its record.
    if(const auto* record=symmetry_record(p);record && closed_delta && record->value("version","")==symmetric_moves_version &&
                                               record->value("status","")!="disabled_for_proposal") {
        out.status=record->value("status",out.status);return out;
    }
    q.identity.added_edges=as_ids(merged_add);q.identity.removed_edges=as_ids(merged_rem);
    q.endpoints=moves;q.allowed_endpoints=allowed;
    std::set<std::string> keys;
    for(const auto& list:{merged_add,merged_rem})for(auto e:list)for(auto f:t.edge_faces[e])keys.insert(parent.surface->chart_key(ct.chart_of_face[f]).digest());
    q.identity.affected_chart_keys=sorted(keys);
    if(!p.resolve_chart_keys.empty()) {
        std::set<std::string> resolve(p.resolve_chart_keys.begin(),p.resolve_chart_keys.end());
        for(MeshIndex ch=0;ch<ct.n_charts;++ch) {
            if(!resolve.contains(parent.surface->chart_key(ch).digest()))continue;
            for(auto f:ct.chart_faces[ch])
                if(map.face_exact[f])resolve.insert(parent.surface->chart_key(ct.chart_of_face[static_cast<std::size_t>(map.maps.face[f])]).digest());
        }
        q.resolve_chart_keys=sorted(resolve);
    }
    Json moves_json=Json::array();for(const auto& mv:mirror_moves)moves_json.push_back({{"from",mv.from},{"to",mv.to}});
    Json original_moves=Json::array();for(const auto& mv:p.endpoints)original_moves.push_back({{"from",mv.from},{"to",mv.to}});
    set_symmetry_record(q,{{"version",symmetric_moves_version},{"status",out.status},{"plane",plane_json(map.plane)},
                           {"original",{{"added",A},{"removed",R},{"endpoints",original_moves}}},
                           {"mirror_added",mirror_added},{"mirror_removed",mirror_removed},{"mirror_endpoints",moves_json},
                           {"deduplicated_edges",deduplicated},{"self_mirrored_edges",self_mirrored},{"unchanged_mirror_edges",unchanged},
                           {"patch_charts",sorted(charts)},{"patch_faces",patch_faces}});
    out.detail={{"mirror_added",mirror_added.size()},{"mirror_removed",mirror_removed.size()},{"deduplicated_edges",deduplicated}};
    return out;
}
nlohmann::json missing_mirror_half(const Proposal& p,std::span<const std::uint8_t> cuts,const StrictMirrorMap& map) {
    if(cuts.size()!=map.edge_exact.size())throw ContractError(ErrorCode::invalid_input,"symmetry.cuts","cut mask does not match the mirror map");
    std::vector<std::int64_t> added,removed,unmapped;std::vector<MeshIndex> vertices;
    for(auto x:p.identity.added_edges) {
        const auto img=x>=0?map.mirror_edge(static_cast<MeshIndex>(x)):std::nullopt;
        if(!img)unmapped.push_back(x);else if(!cuts[*img])added.push_back(x);
    }
    for(auto x:p.identity.removed_edges) {
        const auto img=x>=0?map.mirror_edge(static_cast<MeshIndex>(x)):std::nullopt;
        if(!img)unmapped.push_back(x);else if(cuts[*img])removed.push_back(x);
    }
    for(const auto& move:p.endpoints) {
        const auto f=map.mirror_vertex(move.from),to=map.mirror_vertex(move.to);
        const bool found=f && to && std::any_of(p.endpoints.begin(),p.endpoints.end(),[&](const EndpointChange& c) {return c.from==*f && c.to==*to;});
        if(!found)vertices.push_back(move.from);
    }
    if(added.empty() && removed.empty() && unmapped.empty() && vertices.empty())return nullptr;
    return {{"added_without_image",added},{"removed_without_image",removed},{"unmapped_edges",unmapped},{"endpoint_moves_without_image",vertices}};
}
// ------------------------------------------------------------ context
bool strict_symmetry(const ProposalContext& ctx) {return ctx.options.enhanced.symmetry.mode=="strict";}
std::shared_ptr<const StrictMirrorMap> context_symmetry(ProposalContext& ctx) {
    const auto& o=ctx.options.enhanced.symmetry;
    if(o.mode=="asymmetric")return nullptr;
    if(o.mode!="strict" && o.mode!="approximate")
        throw ContractError(ErrorCode::invalid_option,"enhanced.symmetry.mode","the symmetry mode must be asymmetric, strict or approximate");
    if(o.on_conflict!="reject_pair" && o.on_conflict!="disable_for_proposal")
        throw ContractError(ErrorCode::invalid_option,"enhanced.symmetry.on_conflict","the conflict policy must be reject_pair or disable_for_proposal");
    if(!ctx.mesh)throw ContractError(ErrorCode::invalid_input,"symmetry","missing mesh");
    if(!ctx.symmetry || ctx.symmetry->geometry_key!=ctx.mesh->geometry_key() || !same_map_options(ctx.symmetry->options,o))
        ctx.symmetry=std::make_shared<const StrictMirrorMap>(strict_mirror_map(*ctx.mesh,o));
    return ctx.symmetry;
}
// ------------------------------------------------------------ diagnostics
nlohmann::json symmetry_metrics(const MeshSnapshot& m,const StrictMirrorMap& map,std::span<const std::uint8_t> cuts) {
    const auto& t=m.tables();const auto& g=m.geometry();const auto ne=t.edge_vertices.size();
    if(cuts.size()!=ne)throw ContractError(ErrorCode::invalid_input,"symmetry.cuts","cut mask does not match the mesh edges");
    double exact_length=0,paired_length=0;std::uint64_t seams=0,outside=0;std::vector<MeshIndex> unpaired;
    EdgeMask internal(ne,0);
    for(std::size_t e=0;e<ne;++e) {
        if(!cuts[e] || topology_fact(t,e))continue;
        internal[e]=1;++seams;
        const auto img=map.mirror_edge(static_cast<MeshIndex>(e));
        if(!img) {++outside;continue;}
        exact_length+=g.edge_length[e];
        if(cuts[*img])paired_length+=g.edge_length[e];else unpaired.push_back(static_cast<MeshIndex>(e));
    }
    Json out={{"seam_edges",seams},{"outside_correspondence",outside},{"unpaired_edges",unpaired.size()},
              {"exact_pairing",exact_length>0?Json(paired_length/exact_length):Json(nullptr)}};
    std::vector<MeshIndex> listed(unpaired.begin(),unpaired.begin()+static_cast<std::ptrdiff_t>(std::min<std::size_t>(unpaired.size(),64)));
    out["unpaired"]=listed;
    // T14 seam symmetry (sampled, tolerant): the approximate, metric-only view.
    out["seam_symmetry"]=map.plane && seams?finite_or_null(seam_symmetry(m,internal,map.plane).score):Json(nullptr);
    // Nearest-surface deviation of the reflected vertices (no hard correspondence promised).
    if(map.plane && !m.input().positions.empty() && map.diagonal>0) {
        std::vector<Vec3> reflected;reflected.reserve(m.input().positions.size());
        for(const auto& p:m.input().positions)reflected.push_back(map.plane->reflect(p));
        const auto near=surface_distance(m,reflected);
        double worst=0,sum=0;std::uint64_t reached=0,within=0;
        for(double d:near.distance) {if(!std::isfinite(d))continue;++reached;sum+=d;worst=std::max(worst,d);within+=d<=map.map_radius;}
        out["nearest_surface"]={{"vertices",near.distance.size()},{"reached",reached},{"within_map_radius",within},
                                {"max_relative",reached?Json(worst/map.diagonal):Json(nullptr)},{"mean_relative",reached?Json(sum/static_cast<double>(reached)/map.diagonal):Json(nullptr)}};
    } else out["nearest_surface"]=nullptr;
    return out;
}
nlohmann::json symmetry_summary(ProposalContext& ctx,const RunState& origin,const RunState& retained) {
    const auto map=context_symmetry(ctx);
    if(!map)return nullptr;
    const auto& o=ctx.options.enhanced.symmetry;
    Json pairing=Json::object();
    // GB (m12): a conflict is tallied once per reason ("symmetry_conflict:<reason>", under both
    // conflict policies); the transaction's rejection count ("symmetry_conflict", reject_pair
    // only) counts the same proposals again and stays with the rejections, not in the pairing.
    for(const auto& [generator,tally]:ctx.counters.by_generator)
        for(const auto& [key,count]:tally)if(key.rfind("symmetry",0)==0 && key!="symmetry_conflict")pairing[generator][key]=count;
    const bool strict=o.mode=="strict";
    return {{"version",symmetric_moves_version},{"mode",o.mode},{"on_conflict",o.on_conflict},{"moves_enabled",strict && map->usable()},
            {"label",strict?"strict: every committed move has its mirror half on an exact correspondence":
                            "approximate: metric-only diagnostics; no paired moves and no hard edge equality"},
            {"map",map->to_json()},{"pairing",pairing},{"candidate_zero",symmetry_metrics(*ctx.mesh,*map,origin.cuts)},
            {"retained",symmetry_metrics(*ctx.mesh,*map,retained.cuts)}};
}
// ------------------------------------------------------------ the Section 21 symmetry stage
ProposalGenerator symmetry_generator() {
    return {generator_name,generator_version,[](const RunState& s,ProposalContext& ctx) {
        std::vector<Proposal> out;
        const auto& o=ctx.options.enhanced.symmetry;
        Json diag={{"version",symmetric_moves_version},{"mode",o.mode},{"on_conflict",o.on_conflict}};
        const auto map=context_symmetry(ctx);
        if(!map) {diag["status"]="off";ctx.diagnostics[generator_name]=diag;return out;}
        diag["map"]=map->to_json();diag["incumbent"]=symmetry_metrics(*ctx.mesh,*map,s.cuts);
        if(o.mode!="strict") {
            diag["status"]="diagnostics_only";diag["moves_enabled"]=false;
            diag["label"]="approximate: metric-only diagnostics; no paired moves and no hard edge equality";
            ctx.diagnostics[generator_name]=diag;return out;
        }
        diag["moves_enabled"]=map->usable();
        if(!map->usable()) {diag["status"]="no_strict_map";ctx.diagnostics[generator_name]=diag;return out;}
        const auto& m=*ctx.mesh;const auto& t=m.tables();const auto ne=t.edge_vertices.size();const auto& hc=ctx.constraints;
        std::vector<MeshIndex> unpaired;
        for(std::size_t e=0;e<ne;++e) {
            if(!s.cuts[e] || topology_fact(t,e))continue;
            const auto img=map->mirror_edge(static_cast<MeshIndex>(e));
            if(img && !s.cuts[*img])unpaired.push_back(static_cast<MeshIndex>(e));
        }
        if(unpaired.empty()) {diag["status"]="already_symmetric";ctx.diagnostics[generator_name]=diag;return out;}
        // An edge keeps its state when it is hard-constrained, or when changing it would make the
        // mirror half violate its image's constraint (adding the image of a forbidden edge,
        // removing the image of a locked or required one). Constraints are never mirrored.
        auto kept_cut=[&](MeshIndex e) {return (e<hc.locked.size() && hc.locked[e]) || (e<hc.required.size() && hc.required[e]);};
        auto forbidden=[&](MeshIndex e) {return e<hc.forbidden.size() && hc.forbidden[e];};
        std::set<MeshIndex> pinned;
        for(const auto& [corner,uv]:hc.pins)if(corner<m.input().corner_vertices.size())pinned.insert(static_cast<MeshIndex>(m.input().corner_vertices[corner]));
        auto one_sided_pin=[&](MeshIndex e) {
            for(auto v:t.edge_vertices[e]) {const auto w=map->mirror_vertex(v);if(pinned.contains(v) && (!w || !pinned.contains(*w)))return true;}
            return false;
        };
        auto frozen=[&](MeshIndex e,bool cut) {
            if(kept_cut(e) || forbidden(e))return true;
            const auto img=map->mirror_edge(e);
            return !img || (cut?forbidden(*img):kept_cut(*img)) || one_sided_pin(e) || one_sided_pin(*img);
        };
        const auto routing=context_routing_cost(ctx);
        std::set<EdgeMask> seen{s.cuts};Json variants=Json::array();std::int64_t proposed=0;
        const std::array<std::pair<const char*,int>,3> names{{{"closure",0},{"mirror_negative",-1},{"mirror_positive",1}}};
        for(std::size_t k=0;k<names.size();++k) {
            const auto [name,side]=names[k];
            Json record={{"variant",name}};
            if(proposed>=o.max_candidates) {record["status"]="beyond_candidate_budget";variants.push_back(record);continue;}
            EdgeMask target=s.cuts;std::vector<MeshIndex> residual;
            for(std::size_t e=0;e<ne;++e) {
                if(topology_fact(t,e))continue;
                const auto id=static_cast<MeshIndex>(e);const auto img=map->mirror_edge(id);
                if(!img)continue;
                bool want=s.cuts[e]!=0;
                if(side==0) want=want || s.cuts[*img];
                else {
                    const int es=edge_side(*map,t,id);
                    if(es==-side)want=s.cuts[*img]!=0;          // the overwritten side copies its image
                    else if(es==0 && *img!=id)want=want || s.cuts[*img];
                }
                if(want==(s.cuts[e]!=0))continue;
                if(frozen(id,want)) {residual.push_back(id);continue;}  // asymmetric constraints are never mirrored
                target[e]=want;
            }
            // Companion openings for charts the variant leaves open, with their images.
            auto opening=open_non_disk_charts(ctx.mesh,target,routing,hc,ctx.exempt);
            std::vector<MeshIndex> mirrored;bool refused=false;
            for(auto e:opening.added) {
                const auto img=map->mirror_edge(e);
                if(!img || target[*img])continue;
                if(frozen(*img,true)) {refused=true;residual.push_back(*img);continue;}
                target[*img]=1;mirrored.push_back(*img);
            }
            record["openings"]=opening.added;record["mirrored_openings"]=mirrored;record["residual_edges"]=residual;
            std::vector<MeshIndex> open;
            if(!opening.ok || refused)record["status"]="not_openable";
            else {
                const auto surface=cut_surface(ctx.mesh,target);
                for(MeshIndex ch=0;ch<surface->tables().n_charts;++ch)
                    if(!surface->tables().chart_is_disk[ch] && measured_chart(m,*surface,ch,ctx.exempt))open.push_back(ch);
                if(!open.empty()) {record["status"]="not_openable_symmetric";record["open_charts"]=open;}
                else if(target==s.cuts)record["status"]="unchanged";
                else if(!seen.insert(target).second)record["status"]="duplicate";
                else {
                    std::uint64_t added=0,removed=0;
                    for(std::size_t e=0;e<ne;++e) {added+=target[e] && !s.cuts[e];removed+=!target[e] && s.cuts[e];}
                    auto p=cut_proposal(symmetry_id_base+k,s,target,generator_name,generator_version,
                                        {{"operation","symmetrize"},{"variant",name},{"unpaired_edges",unpaired},{"residual_edges",residual},
                                         {"openings",opening.added},{"mirrored_openings",mirrored},{"added",added},{"removed",removed}});
                    record["status"]="proposed";record["id"]=p.identity.id;record["added"]=added;record["removed"]=removed;
                    out.push_back(std::move(p));++proposed;
                }
            }
            variants.push_back(record);
        }
        diag["status"]=out.empty()?"no_candidate":"proposed";diag["unpaired_edges"]=unpaired.size();diag["variants"]=variants;
        ctx.diagnostics[generator_name]=diag;
        return out;
    }};
}
}
