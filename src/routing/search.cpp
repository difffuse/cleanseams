// SPDX-License-Identifier: GPL-3.0-or-later
#include "cleanseams/routing.hpp"
#include "cleanseams/endpoints.hpp"
#include "cleanseams/errors.hpp"
#include "cleanseams/detail/mesh_math.hpp"
#include "cleanseams/detail/pyset.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <numbers>
#include <queue>
#include <set>
#include <tuple>
namespace cleanseams {
namespace {
void require(bool yes,const char* field,const char* message) {if(!yes)throw ContractError(ErrorCode::invalid_input,field,message);}
}
struct RoutingAccess {
    static void prepare(RoutingWorkspace& w,const MeshSnapshot& m) {
        const auto directed=2*m.tables().edge_vertices.size();
        if(w.mesh_!=&m || w.mesh_key_!=m.geometry_key() || w.stamp_.size()!=directed) {
            w.mesh_=&m;w.mesh_key_=m.geometry_key();w.stamp_.assign(directed,0);w.best_.assign(directed,0);w.pred_.assign(directed,-1);
            w.tangent_stamp_.assign(directed,0);w.tangent_.assign(directed,Vec3{0,0,0});
        }
    }
    static std::uint32_t next_generation(RoutingWorkspace& w) {
        if(w.generation_==std::numeric_limits<std::uint32_t>::max()) {
            std::fill(w.stamp_.begin(),w.stamp_.end(),0);w.generation_=0;++w.stats_.stamp_resets;
        }
        return ++w.generation_;
    }
    static std::size_t directed(const MeshSnapshot& m,MeshIndex from,MeshIndex edge) {
        return 2*static_cast<std::size_t>(edge)+(m.tables().edge_vertices[edge][0]==from?0:1);
    }
    // Directed tangent of `edge` leaving `vertex`, computed once per mesh.
    static const Vec3& tangent(RoutingWorkspace& w,const MeshSnapshot& m,MeshIndex vertex,MeshIndex edge) {
        const auto d=directed(m,vertex,edge);++w.stats_.tangent_lookups;
        if(!w.tangent_stamp_[d]) {w.tangent_[d]=m.tangent_dir(vertex,edge);w.tangent_stamp_[d]=1;++w.stats_.tangent_computed;}
        return w.tangent_[d];
    }
    static RoutingStats& stats(RoutingWorkspace& w) {return w.stats_;}
    static std::vector<std::uint32_t>& stamps(RoutingWorkspace& w) {return w.stamp_;}
    static std::vector<double>& best(RoutingWorkspace& w) {return w.best_;}
    static std::vector<std::int64_t>& pred(RoutingWorkspace& w) {return w.pred_;}
};
Corridor build_corridor(const CutSurface& cs,const Chain& chain,std::int64_t rings,std::optional<std::span<const std::uint8_t>> protect) {
    const auto& m=cs.mesh();const auto& t=m.tables();const auto& ct=cs.tables();
    const auto nv=m.input().positions.size(),ne=t.edge_vertices.size();
    require(!chain.verts.empty(),"routing.chain","empty chain");
    if(protect)require(protect->size()==ne,"routing.protected","protected mask does not match the mesh edges");
    std::set<MeshIndex> allowed;if(chain.charts)allowed.insert(chain.charts->begin(),chain.charts->end());
    std::vector<std::uint8_t> chain_edge(ne,0);for(auto e:chain.edges)chain_edge.at(e)=1;
    const auto a=chain.verts.front(),b=chain.verts.back();
    auto edge_ok=[&](MeshIndex e) {
        if(t.edge_boundary[e] || !t.edge_manifold[e])return false;
        if(protect && (*protect)[e])return false;
        if(ct.cut[e] && !chain_edge[e])return false;
        for(auto f:t.edge_faces[e])if(!allowed.contains(ct.chart_of_face[f]))return false;
        return true;
    };
    auto vert_ok=[&](MeshIndex v) {
        if(v==a || v==b)return true;
        if(!t.fans[v] || t.vert_boundary[v])return false;
        for(auto e:t.vert_edges[v])if(ct.cut[e] && !chain_edge[e])return false;
        return true;
    };
    detail::PyIntSet chain_verts;for(auto v:chain.verts)chain_verts.add(v);
    Corridor out;out.rings=rings;out.dist.assign(nv,-1);
    std::vector<MeshIndex> frontier;for(auto v:chain_verts.order()) {frontier.push_back(static_cast<MeshIndex>(v));out.dist[v]=0;}
    auto verts=detail::PyIntSet::copy_of(chain_verts);
    for(std::int64_t r=1;r<=rings;++r) {
        std::vector<MeshIndex> next;
        for(auto v:frontier)for(auto e:t.vert_edges[v]) {
            if(!edge_ok(e))continue;
            const auto w=m.edge_other(e,v);
            if(out.dist[w]>=0 || !vert_ok(w))continue;
            out.dist[w]=r;verts.add(w);next.push_back(w);
        }
        frontier=std::move(next);
    }
    out.vertex_mask.assign(nv,0);for(auto v:verts.order()) {out.verts.push_back(static_cast<MeshIndex>(v));out.vertex_mask[v]=1;}
    detail::PyIntSet edges;
    for(auto v:out.verts)for(auto e:t.vert_edges[v])if(edge_ok(e) && out.vertex_mask[m.edge_other(e,v)])edges.add(e);
    out.edge_mask.assign(ne,0);for(auto e:edges.order()) {out.edges.push_back(static_cast<MeshIndex>(e));out.edge_mask[e]=1;}
    out.lref=out.edges.empty()?m.mean_edge_length():m.mean_edge_length(std::span<const MeshIndex>(out.edges));
    return out;
}
AnchorTangent anchor_tangent(const CutSurface& cs,const Chain& chain,int end) {
    const auto& m=cs.mesh();const auto& t=m.tables();const auto& ct=cs.tables();
    require(!chain.verts.empty() && !chain.edges.empty(),"routing.chain","chain without edges");
    const auto v=end==0?chain.verts.front():chain.verts.back();
    const std::set<MeshIndex> chain_edges(chain.edges.begin(),chain.edges.end());
    std::vector<MeshIndex> others;
    for(auto e:t.vert_edges[v])if(ct.cut[e] && !chain_edges.contains(e) && !t.edge_boundary[e])others.push_back(e);
    if(others.size()==1) {
        const auto d=m.tangent_dir(v,others[0]);
        return {end==0?detail::mul(d,-1.0):d,1.0};
    }
    const auto d=m.tangent_dir(v,end==0?chain.edges.front():chain.edges.back());
    return {end==0?d:detail::mul(d,-1.0),0.5};
}
double turn_cost(const Vec3& in,const Vec3& out) {
    const double c=detail::dot(in,out);
    return (1.0-std::max(-1.0,std::min(1.0,c)))/2.0;
}
namespace {
struct Pricing {
    const CutSurface& cs;const Corridor& corridor;const RouteWeights& weights;const PathOptions& opts;
    std::span<const double> exposure;RoutingWorkspace* ws;double lref;std::int64_t radius;
    Vec3 tangent(MeshIndex v,MeshIndex e) const {return ws?RoutingAccess::tangent(*ws,cs.mesh(),v,e):cs.mesh().tangent_dir(v,e);}
    // Source step_cost; `e_in` is ignored without a previous vertex.
    double step(bool has_prev,MeshIndex cur,MeshIndex next,MeshIndex e_in,MeshIndex e_out) const {
        const auto& m=cs.mesh();const auto& t=m.tables();const auto& g=m.geometry();
        double turn=0.0,flow=0.0;
        if(has_prev) {
            turn=turn_cost(detail::mul(tangent(cur,e_in),-1.0),tangent(cur,e_out));
            if(t.vert_regular[cur] && t.fans[cur]->straight(e_in)!=static_cast<TransferIndex>(e_out))flow=1.0;
        }
        const auto d=corridor.dist[next]>=0?corridor.dist[next]:radius;
        const double drift=std::min(1.0,std::pow(static_cast<double>(d)/static_cast<double>(radius),2.0));
        const double convex=std::max(0.0,std::min(1.0,g.edge_dihedral[e_out]/std::numbers::pi));
        const double expo=exposure.empty()?0.0:exposure[e_out];
        const double length=g.edge_length[e_out]/lref;
        return length*(1.0+weights[2]*drift+opts.w_convex*convex+opts.w_exposure*expo)+weights[0]*turn+weights[1]*flow;
    }
    double terminal(double weight,const Vec3& a,const Vec3& b) const {
        return weights[0]*weight*opts.anchor_tangent_weight*2.0*turn_cost(a,b);
    }
};
}
RouteSearch find_path(const CutSurface& cs,const Chain& chain,const Corridor& corridor,const RouteWeights& weights,const PathOptions& opts,
                      std::span<const double> exposure,RoutingWorkspace& ws) {
    const auto& m=cs.mesh();const auto& t=m.tables();const auto ne=t.edge_vertices.size();
    require(chain.verts.size()>=2,"routing.chain","chain needs two anchors");
    require(corridor.edge_mask.size()==ne && corridor.dist.size()==m.input().positions.size(),"routing.corridor","corridor belongs to another mesh");
    require(exposure.empty() || exposure.size()==ne,"routing.exposure","cost field does not match the mesh edges");
    RoutingAccess::prepare(ws,m);auto& stats=RoutingAccess::stats(ws);++stats.searches;
    const auto generation=RoutingAccess::next_generation(ws);
    auto& stamp=RoutingAccess::stamps(ws);auto& best=RoutingAccess::best(ws);auto& pred=RoutingAccess::pred(ws);
    const auto a=chain.verts.front(),b=chain.verts.back();
    const auto start=anchor_tangent(cs,chain,0),finish=anchor_tangent(cs,chain,1);
    const Pricing price{cs,corridor,weights,opts,exposure,&ws,corridor.lref!=0?corridor.lref:1.0,std::max<std::int64_t>(1,corridor.rings)};
    // Safe early stop needs every increment and terminal cost to be nonnegative.
    bool nonnegative=weights[0]>=0 && weights[1]>=0 && weights[2]>=0 && opts.w_convex>=0 && opts.w_exposure>=0 && opts.anchor_tangent_weight>=0;
    for(double x:exposure)nonnegative=nonnegative && x>=0;
    stats.nonnegative_costs=stats.nonnegative_costs && nonnegative;
    constexpr double inf=std::numeric_limits<double>::infinity();
    auto cost_of=[&](std::size_t d) {return stamp[d]==generation?best[d]:inf;};
    using Item=std::tuple<double,std::int64_t,std::int64_t>;
    std::priority_queue<Item,std::vector<Item>,std::greater<Item>> heap;
    for(auto e:t.vert_edges[a]) {
        if(!corridor.edge_mask[e])continue;
        const auto w=m.edge_other(e,a);double c=price.step(false,a,w,0,e);
        if(start.direction)c+=price.terminal(start.weight,*start.direction,price.tangent(a,e));
        const auto state=RoutingAccess::directed(m,a,e);
        if(c<cost_of(state)) {stamp[state]=generation;best[state]=c;pred[state]=-1;heap.push({c,a,w});++stats.pushes;}
    }
    double goal=inf;std::int64_t goal_state=-1;
    while(!heap.empty()) {
        const auto [c,prev_i,cur_i]=heap.top();heap.pop();++stats.pops;
        if(nonnegative && c>=goal) {++stats.early_stops;break;} // nothing left can improve the goal
        const auto prev=static_cast<MeshIndex>(prev_i),cur=static_cast<MeshIndex>(cur_i);
        const auto e_in=static_cast<MeshIndex>(m.edge_index(prev,cur));const auto state=RoutingAccess::directed(m,prev,e_in);
        if(c>cost_of(state))continue;
        if(cur==b) {
            double total=c;
            if(finish.direction)total+=price.terminal(finish.weight,detail::mul(price.tangent(cur,e_in),-1.0),*finish.direction);
            if(total<goal) {goal=total;goal_state=static_cast<std::int64_t>(state);}
            continue;
        }
        if(cur==a)continue;
        for(auto e_out:t.vert_edges[cur]) {
            if(!corridor.edge_mask[e_out] || e_out==e_in)continue;
            const auto next=m.edge_other(e_out,cur);
            if(next==prev || next==a)continue;
            const double nc=c+price.step(true,cur,next,e_in,e_out);
            if(nc>=goal)continue;
            const auto target=RoutingAccess::directed(m,cur,e_out);
            if(nc<cost_of(target)) {stamp[target]=generation;best[target]=nc;pred[target]=static_cast<std::int64_t>(state);heap.push({nc,cur,next});++stats.pushes;}
        }
    }
    RouteSearch out;
    if(goal_state<0) {out.reason="no_route";return out;}
    std::vector<MeshIndex> verts;
    for(auto s=goal_state;s>=0;s=pred[static_cast<std::size_t>(s)]) {
        const auto e=static_cast<std::size_t>(s)/2;const auto& ends=t.edge_vertices[e];
        verts.push_back(s%2?ends[0]:ends[1]); // the state's current vertex
    }
    verts.push_back(a);std::reverse(verts.begin(),verts.end());
    auto unique=verts;std::sort(unique.begin(),unique.end());
    if(std::adjacent_find(unique.begin(),unique.end())!=unique.end()) {out.reason="route_not_simple";return out;}
    out.verts=std::move(verts);out.cost=goal;
    return out;
}
double route_cost(const CutSurface& cs,const Chain& chain,std::span<const MeshIndex> verts,const Corridor& corridor,
                  const RouteWeights& weights,const PathOptions& opts,std::span<const double> exposure,RoutingWorkspace* ws) {
    const auto& m=cs.mesh();
    if(ws)RoutingAccess::prepare(*ws,m);
    const Pricing price{cs,corridor,weights,opts,exposure,ws,corridor.lref!=0?corridor.lref:1.0,std::max<std::int64_t>(1,corridor.rings)};
    const auto start=anchor_tangent(cs,chain,0),finish=anchor_tangent(cs,chain,1);
    auto edge=[&](MeshIndex x,MeshIndex y) {const auto e=m.edge_index(x,y);require(e>=0,"routing.route","route step is not a mesh edge");return static_cast<MeshIndex>(e);};
    double total=0.0;
    for(std::size_t k=0;k+1<verts.size();++k) {
        const auto cur=verts[k],next=verts[k+1],e_out=edge(cur,next);
        if(k==0 && start.direction)total+=price.terminal(start.weight,*start.direction,price.tangent(cur,e_out));
        total+=price.step(k>0,cur,next,k>0?edge(verts[k-1],cur):0,e_out);
    }
    if(finish.direction && verts.size()>1) {
        const auto e_in=edge(verts[verts.size()-2],verts.back());
        total+=price.terminal(finish.weight,detail::mul(price.tangent(verts.back(),e_in),-1.0),*finish.direction);
    }
    return total;
}
bool route_covers_chain(const MeshSnapshot& m,const Corridor& corridor,std::span<const MeshIndex> route,const Chain& chain,std::int64_t radius) {
    std::vector<std::int64_t> dist(m.input().positions.size(),-1);std::vector<MeshIndex> frontier;
    for(auto v:route)if(dist[v]<0) {dist[v]=0;frontier.push_back(v);}
    for(std::int64_t r=1;r<=radius;++r) {
        std::vector<MeshIndex> next;
        for(auto v:frontier)for(auto e:m.tables().vert_edges[v]) {
            if(!corridor.edge_mask[e])continue;
            const auto w=m.edge_other(e,v);if(dist[w]<0) {dist[w]=r;next.push_back(w);}
        }
        frontier=std::move(next);
    }
    return std::all_of(chain.verts.begin(),chain.verts.end(),[&](auto v) {return dist[v]>=0;});
}
std::optional<std::vector<std::uint8_t>> apply_route(const MeshSnapshot& m,std::span<const std::uint8_t> seams,const Chain& chain,std::span<const MeshIndex> verts) {
    require(seams.size()==m.tables().edge_vertices.size(),"routing.seams","seam mask does not match the mesh edges");
    std::vector<std::uint8_t> out(seams.begin(),seams.end());
    for(auto e:chain.edges)out.at(e)=0;
    for(std::size_t k=0;k+1<verts.size();++k) {
        const auto e=m.edge_index(verts[k],verts[k+1]);if(e<0)return std::nullopt;out[static_cast<std::size_t>(e)]=1;
    }
    return out;
}
std::vector<std::uint8_t> corridor_band_faces(const MeshSnapshot& m,const Corridor& corridor) {
    std::vector<std::uint8_t> band(m.input().face_offsets.size()-1,0);
    for(auto v:corridor.verts)for(auto f:m.tables().vert_faces[v])band[f]=1;
    return band;
}
// ============================================================ Phase B multi-endpoint routing (T30)
// Section 13.2: direction-aware routing from several sources to several goals, for sliding
// seam endpoints. Legacy routing above is untouched; only enhanced_v1 generators call this.
namespace {
void require_cost(double x,const char* field) {
    if(!std::isfinite(x) || x<0)throw ContractError(ErrorCode::invalid_input,field,"route costs must be finite and nonnegative");
}
struct Ends {
    std::vector<RouteEndpoint> list;          // ascending vertex
    std::vector<std::int64_t> index;          // per vertex: position in list or -1
};
Ends prepare_ends(const std::vector<RouteEndpoint>& given,std::size_t nv,const char* field) {
    Ends out;out.list=given;out.index.assign(nv,-1);
    std::stable_sort(out.list.begin(),out.list.end(),[](const RouteEndpoint& a,const RouteEndpoint& b) {return a.vertex<b.vertex;});
    for(std::size_t k=0;k<out.list.size();++k) {
        const auto& e=out.list[k];
        require(e.vertex<nv,field,"route endpoint out of range");
        require(out.index[e.vertex]<0,field,"repeated route endpoint");
        require_cost(e.displacement,field);require_cost(e.tangent_weight,field);
        if(e.tangent)require(std::isfinite((*e.tangent)[0]) && std::isfinite((*e.tangent)[1]) && std::isfinite((*e.tangent)[2]),field,"terminal direction must be finite");
        out.index[e.vertex]=static_cast<std::int64_t>(k);
    }
    return out;
}
double terminal_cost(const RouteEndpoint& end,const Vec3& a,const Vec3& b) {
    return end.tangent && end.tangent_weight>0?end.tangent_weight*turn_cost(a,b):0.0;
}
void check_query(const MeshSnapshot& m,const MultiEndpointQuery& q) {
    const auto ne=m.tables().edge_vertices.size(),nv=m.input().positions.size();
    require(q.edge_cost.size()==ne,"routing.multi.edge_cost","one cost per mesh edge");
    require(q.edge_allowed.size()==ne,"routing.multi.edge_allowed","one flag per mesh edge");
    require(q.vertex_allowed.size()==nv,"routing.multi.vertex_allowed","one flag per mesh vertex");
    require_cost(q.w_turn,"routing.multi.w_turn");
    for(std::size_t e=0;e<ne;++e)if(q.edge_allowed[e])require_cost(q.edge_cost[e],"routing.multi.edge_cost");
}
}
nlohmann::json MultiEndpointRoute::to_json() const {
    return {{"verts",verts?nlohmann::json(*verts):nlohmann::json(nullptr)},{"edges",edges},{"cost",cost},
            {"source",source?nlohmann::json(*source):nlohmann::json(nullptr)},{"goal",goal?nlohmann::json(*goal):nlohmann::json(nullptr)},
            {"reason",reason},{"pops",pops},{"pushes",pushes},{"early_stop",early_stop}};
}
MultiEndpointRoute find_multi_endpoint_path(const MeshSnapshot& m,const MultiEndpointQuery& q,RoutingWorkspace& ws) {
    const auto& t=m.tables();const auto& g=m.geometry();const auto nv=m.input().positions.size();
    check_query(m,q);
    const auto sources=prepare_ends(q.sources,nv,"routing.multi.sources"),goals=prepare_ends(q.goals,nv,"routing.multi.goals");
    MultiEndpointRoute out;
    if(sources.list.empty() || goals.list.empty()) {out.reason="no_endpoints";return out;}
    for(const auto& s:sources.list)require(goals.index[s.vertex]<0,"routing.multi.goals","a vertex cannot be both a source and a goal");
    RoutingAccess::prepare(ws,m);auto& stats=RoutingAccess::stats(ws);++stats.searches;
    const auto generation=RoutingAccess::next_generation(ws);
    auto& stamp=RoutingAccess::stamps(ws);auto& best=RoutingAccess::best(ws);auto& pred=RoutingAccess::pred(ws);
    constexpr double inf=std::numeric_limits<double>::infinity();
    auto cost_of=[&](std::size_t d) {return stamp[d]==generation?best[d]:inf;};
    auto tangent=[&](MeshIndex v,MeshIndex e) -> const Vec3& {return RoutingAccess::tangent(ws,m,v,e);};
    // A vertex may be entered when it is a goal (terminal) or an allowed interior vertex that is
    // not a source: sources only start routes, goals only end them.
    auto enterable=[&](MeshIndex v) {return sources.index[v]<0 && (goals.index[v]>=0 || q.vertex_allowed[v]);};
    using Item=std::tuple<double,std::int64_t,std::int64_t>;
    std::priority_queue<Item,std::vector<Item>,std::greater<Item>> heap;
    for(const auto& s:sources.list)for(auto e:t.vert_edges[s.vertex]) {
        if(!q.edge_allowed[e])continue;
        const auto w=m.edge_other(e,s.vertex);if(!enterable(w))continue;
        const double c=s.displacement+(s.tangent?terminal_cost(s,*s.tangent,tangent(s.vertex,e)):0.0)+q.edge_cost[e];
        const auto state=RoutingAccess::directed(m,s.vertex,e);
        if(c<cost_of(state)) {stamp[state]=generation;best[state]=c;pred[state]=-1;heap.push({c,s.vertex,w});++out.pushes;++stats.pushes;}
    }
    double goal=inf;std::int64_t goal_state=-1;
    while(!heap.empty()) {
        const auto [c,prev_i,cur_i]=heap.top();heap.pop();++out.pops;++stats.pops;
        // Lower-bound stop: every remaining state costs at least c, and extensions and arrivals
        // only add nonnegative terms (checked above), so nothing can beat the goal.
        if(c>=goal) {out.early_stop=true;++stats.early_stops;break;}
        const auto prev=static_cast<MeshIndex>(prev_i),cur=static_cast<MeshIndex>(cur_i);
        const auto e_in=static_cast<MeshIndex>(m.edge_index(prev,cur));const auto state=RoutingAccess::directed(m,prev,e_in);
        if(c>cost_of(state))continue;
        const auto arrival=detail::mul(tangent(cur,e_in),-1.0);
        if(goals.index[cur]>=0) {
            const auto& end=goals.list[static_cast<std::size_t>(goals.index[cur])];
            const double total=c+end.displacement+(end.tangent?terminal_cost(end,arrival,*end.tangent):0.0);
            if(total<goal) {goal=total;goal_state=static_cast<std::int64_t>(state);}
            continue;
        }
        for(auto e_out:t.vert_edges[cur]) {
            if(!q.edge_allowed[e_out] || e_out==e_in)continue;
            const auto next=m.edge_other(e_out,cur);
            if(next==prev || !enterable(next))continue;
            const double turn=q.w_turn>0?q.w_turn*0.5*(g.edge_length[e_in]+g.edge_length[e_out])*turn_cost(arrival,tangent(cur,e_out)):0.0;
            const double nc=c+q.edge_cost[e_out]+turn;
            if(nc>=goal)continue;
            const auto target=RoutingAccess::directed(m,cur,e_out);
            if(nc<cost_of(target)) {stamp[target]=generation;best[target]=nc;pred[target]=static_cast<std::int64_t>(state);heap.push({nc,cur,next});++out.pushes;++stats.pushes;}
        }
    }
    if(goal_state<0) {out.reason="no_route";return out;}
    std::vector<MeshIndex> verts;std::int64_t first=goal_state;
    for(auto s=goal_state;s>=0;s=pred[static_cast<std::size_t>(s)]) {
        const auto e=static_cast<std::size_t>(s)/2;const auto& ends=t.edge_vertices[e];
        verts.push_back(s%2?ends[0]:ends[1]);first=s;
    }
    {const auto e=static_cast<std::size_t>(first)/2;const auto& ends=t.edge_vertices[e];verts.push_back(first%2?ends[1]:ends[0]);}
    std::reverse(verts.begin(),verts.end());
    auto unique=verts;std::sort(unique.begin(),unique.end());
    if(std::adjacent_find(unique.begin(),unique.end())!=unique.end()) {out.reason="route_not_simple";return out;}
    for(std::size_t k=0;k+1<verts.size();++k)out.edges.push_back(static_cast<MeshIndex>(m.edge_index(verts[k],verts[k+1])));
    out.source=verts.front();out.goal=verts.back();out.cost=goal;out.verts=std::move(verts);
    return out;
}
double multi_endpoint_route_cost(const MeshSnapshot& m,const MultiEndpointQuery& q,std::span<const MeshIndex> verts) {
    const auto& g=m.geometry();const auto nv=m.input().positions.size();
    check_query(m,q);
    require(verts.size()>=2,"routing.multi.route","a route needs two vertices");
    const auto sources=prepare_ends(q.sources,nv,"routing.multi.sources"),goals=prepare_ends(q.goals,nv,"routing.multi.goals");
    require(sources.index[verts.front()]>=0 && goals.index[verts.back()]>=0,"routing.multi.route","a route runs from a source to a goal");
    auto edge=[&](MeshIndex x,MeshIndex y) {const auto e=m.edge_index(x,y);require(e>=0,"routing.multi.route","route step is not a mesh edge");return static_cast<MeshIndex>(e);};
    const auto& s=sources.list[static_cast<std::size_t>(sources.index[verts.front()])];
    const auto& end=goals.list[static_cast<std::size_t>(goals.index[verts.back()])];
    double total=s.displacement;
    for(std::size_t k=0;k+1<verts.size();++k) {
        const auto e_out=edge(verts[k],verts[k+1]);
        if(k==0) {if(s.tangent)total+=terminal_cost(s,*s.tangent,m.tangent_dir(verts[0],e_out));}
        else {
            const auto e_in=edge(verts[k-1],verts[k]);
            if(q.w_turn>0)total+=q.w_turn*0.5*(g.edge_length[e_in]+g.edge_length[e_out])*turn_cost(detail::mul(m.tangent_dir(verts[k],e_in),-1.0),m.tangent_dir(verts[k],e_out));
        }
        total+=q.edge_cost[e_out];
    }
    const auto e_last=edge(verts[verts.size()-2],verts.back());
    total+=end.displacement;
    if(end.tangent)total+=terminal_cost(end,detail::mul(m.tangent_dir(verts.back(),e_last),-1.0),*end.tangent);
    return total;
}
nlohmann::json PatchCorridor::to_json() const {
    return {{"radius",std::isfinite(radius)?nlohmann::json(radius):nlohmann::json(nullptr)},{"whole_patch",whole_patch},
            {"patch_vertices",patch_vertices},{"band_vertices",band_vertices},{"band_edges",band_edges}};
}
PatchCorridor patch_corridor(const MeshSnapshot& m,std::span<const std::uint8_t> patch_faces,std::span<const std::uint8_t> edge_usable,
                             std::span<const std::uint8_t> vertex_usable,std::span<const MeshIndex> seeds,std::optional<double> radius) {
    const auto& t=m.tables();const auto& g=m.geometry();
    const auto ne=t.edge_vertices.size(),nv=m.input().positions.size(),nf=m.input().face_offsets.size()-1;
    require(patch_faces.size()==nf,"routing.patch.faces","one flag per face");
    require(edge_usable.size()==ne,"routing.patch.edges","one flag per mesh edge");
    require(vertex_usable.size()==nv,"routing.patch.vertices","one flag per mesh vertex");
    if(radius)require(std::isfinite(*radius) && *radius>=0,"routing.patch.radius","radius must be finite and nonnegative");
    // Patch edges: every face they touch is in the patch. Distances follow them only.
    std::vector<std::uint8_t> patch_edge(ne,0),patch_vertex(nv,0);
    for(std::size_t e=0;e<ne;++e) {
        bool inside=!t.edge_faces[e].empty();for(auto f:t.edge_faces[e])inside=inside && patch_faces[f];
        patch_edge[e]=inside;
    }
    for(std::size_t f=0;f<nf;++f)if(patch_faces[f])for(auto c:m.face_corners(static_cast<MeshIndex>(f)))patch_vertex[static_cast<std::size_t>(m.input().corner_vertices[c])]=1;
    PatchCorridor out;out.distance.assign(nv,std::numeric_limits<double>::infinity());
    using Item=std::pair<double,MeshIndex>;std::priority_queue<Item,std::vector<Item>,std::greater<>> heap;
    for(auto s:seeds) {
        require(s<nv && patch_vertex[s],"routing.patch.seeds","seeds must be patch vertices");
        if(out.distance[s]>0) {out.distance[s]=0;heap.push({0.0,s});}
    }
    while(!heap.empty()) {
        const auto [d,v]=heap.top();heap.pop();
        if(d>out.distance[v])continue;
        for(auto e:t.vert_edges[v]) {
            if(!patch_edge[e])continue;
            const auto w=m.edge_other(e,v);const double nd=d+g.edge_length[e];
            if(nd<out.distance[w]) {out.distance[w]=nd;heap.push({nd,w});}
        }
    }
    if(radius)out.radius=*radius;
    std::vector<std::uint8_t> band(nv,0);
    for(std::size_t v=0;v<nv;++v) {
        out.patch_vertices+=patch_vertex[v];
        band[v]=patch_vertex[v] && std::isfinite(out.distance[v]) && (!radius || out.distance[v]<=*radius);
        out.band_vertices+=band[v];
    }
    out.whole_patch=out.band_vertices==out.patch_vertices;
    out.edge_allowed.assign(ne,0);out.vertex_allowed.assign(nv,0);
    for(std::size_t e=0;e<ne;++e) {
        const auto& ends=t.edge_vertices[e];
        out.edge_allowed[e]=edge_usable[e] && patch_edge[e] && band[ends[0]] && band[ends[1]];out.band_edges+=out.edge_allowed[e];
    }
    for(std::size_t v=0;v<nv;++v)out.vertex_allowed[v]=vertex_usable[v] && band[v];
    return out;
}
}
