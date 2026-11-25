// SPDX-License-Identifier: GPL-3.0-or-later
#include "cleanseams/cache.hpp"
#include "cleanseams/errors.hpp"
#include <atomic>
#include <set>
namespace cleanseams {
namespace {
struct Counters {
    std::atomic<std::uint64_t> cut_builds{0},unchanged{0},geometry_builds{0},duplicate_geometry{0},hits{0},misses{0},evictions{0},solves{0},solve_hits{0};
    std::mutex mutex;
    // Identities built in the counting window: (mesh geometry key, effective cut hash).
    std::set<std::pair<std::string,std::uint64_t>> cuts,geometries;
};
Counters& counters() {static Counters c;return c;}
// The identity sets are bounded (Gate A F31): a long session that never resets the
// counters starts a new identity window instead of growing without limit.
constexpr std::size_t identity_window=std::size_t{1}<<16;
bool first_in_window(std::set<std::pair<std::string,std::uint64_t>>& seen,std::pair<std::string,std::uint64_t> identity) {
    if(seen.size()>=identity_window && !seen.contains(identity))seen.clear();
    return seen.insert(std::move(identity)).second;
}
thread_local TopologyCache* active=nullptr;
std::uint64_t surface_bytes(const CutSurface& cs) {
    const auto& t=cs.tables();std::uint64_t n=0;
    auto add=[&](const auto& v) {n+=static_cast<std::uint64_t>(v.size()*sizeof(typename std::decay_t<decltype(v)>::value_type));};
    add(t.cut);add(t.sector_of_corner);add(t.chart_of_face);add(t.chart_of_corner);add(t.sector_vert);add(t.sector_chart);add(t.boundary_corners);
    add(t.loops.values);add(t.chart_loops.values);add(t.chart_faces.values);add(t.cut_degree);
    for(const auto& a:t.adjacency) {add(a.edges);add(a.arcs.values);}
    return n;
}
std::uint64_t geometry_bytes(const QualityGeometry& g) {
    return static_cast<std::uint64_t>(g.frames().size()*sizeof(TriangleFrame)+g.chart_triangles().values.size()*sizeof(MeshIndex));
}
}
namespace detail {
void note_chart_solve() {++counters().solves;}
std::uint64_t cut_hash(std::span<const std::uint8_t> cut) {
    std::uint64_t h=1469598103934665603ULL;for(std::size_t i=0;i<cut.size();++i)if(cut[i]) {h^=i+1;h*=1099511628211ULL;}
    return h^cut.size();
}
void note_cut_build(const MeshSnapshot& mesh,std::span<const std::uint8_t> cut) {
    auto& c=counters();++c.cut_builds;std::lock_guard lock(c.mutex);
    if(!first_in_window(c.cuts,{mesh.geometry_key(),cut_hash(cut)}))++c.unchanged;
}
void note_geometry_build(const CutSurface& cs) {
    auto& c=counters();++c.geometry_builds;std::lock_guard lock(c.mutex);
    if(!first_in_window(c.geometries,{cs.mesh().geometry_key(),cut_hash(cs.tables().cut)}))++c.duplicate_geometry;
}
}
WorkCounters work_counters() {
    auto& c=counters();
    return {c.cut_builds.load(),c.unchanged.load(),c.geometry_builds.load(),c.duplicate_geometry.load(),c.hits.load(),c.misses.load(),c.evictions.load(),
            c.solves.load(),c.solve_hits.load()};
}
void reset_work_counters() {
    auto& c=counters();std::lock_guard lock(c.mutex);
    c.cut_builds=0;c.unchanged=0;c.geometry_builds=0;c.duplicate_geometry=0;c.hits=0;c.misses=0;c.evictions=0;c.solves=0;c.solve_hits=0;
    c.cuts.clear();c.geometries.clear();
}
TopologyCache::TopologyCache(std::uint64_t budget):budget_(budget) {}
std::uint64_t TopologyCache::bytes() const {std::lock_guard lock(mutex_);return bytes_;}
std::uint64_t TopologyCache::peak_bytes() const {std::lock_guard lock(mutex_);return peak_;}
std::uint64_t TopologyCache::SolvedMap::bytes(const std::string& contract) const {
    return static_cast<std::uint64_t>(corner_uv.size()*sizeof(Vec2)+sector_uv.size()*sizeof(Vec2)+
                                      (corner_ids.size()+sector_corners.size()+pin_corners.size())*sizeof(MeshIndex)+contract.size());
}
void TopologyCache::evict_locked() {
    // GB performance: least recently used surfaces go first, then solved maps. A surface and
    // its geometry are large and rebuilt in milliseconds, while a solved map is small and
    // costs a chart solve (evicting maps first re-solved hundreds of charts on 12k-face
    // meshes). Cached and recomputed values are identical, so only the work changes. A map
    // is released with exactly the bytes store() charged for it.
    while(bytes_>budget_ && !entries_.empty()) {bytes_-=entries_.back().bytes;entries_.pop_back();++counters().evictions;}
    while(bytes_>budget_ && !maps_.empty()) {
        const auto size=maps_.back().second.bytes(maps_.back().first);
        bytes_-=size;map_bytes_-=size;map_index_.erase(maps_.back().first);maps_.pop_back();++counters().evictions;
    }
}
std::optional<TopologyCache::SolvedMap> TopologyCache::solved(const std::string& contract) {
    std::lock_guard lock(mutex_);const auto found=map_index_.find(contract);
    if(found==map_index_.end())return std::nullopt;
    maps_.splice(maps_.begin(),maps_,found->second);++counters().solve_hits;return found->second->second;
}
void TopologyCache::store(const std::string& contract,SolvedMap map) {
    const auto size=map.bytes(contract);
    std::lock_guard lock(mutex_);
    if(size>budget_ || map_index_.contains(contract))return;
    maps_.emplace_front(contract,std::move(map));map_index_[contract]=maps_.begin();bytes_+=size;map_bytes_+=size;
    evict_locked();peak_=std::max(peak_,bytes_);
}
std::shared_ptr<const CutSurface> TopologyCache::surface(const std::shared_ptr<const MeshSnapshot>& mesh,std::span<const std::uint8_t> cut,
                                                        std::uint64_t memory_budget) {
    const auto& t=mesh->tables();
    if(cut.size()!=t.edge_vertices.size())throw ContractError(ErrorCode::invalid_input,"cache.cut","edge mask does not match the mesh edges");
    // The effective cut, as CutSurface builds it, identifies the topology exactly.
    std::vector<std::uint8_t> effective(cut.size());
    for(std::size_t e=0;e<cut.size();++e)effective[e]=cut[e] || t.edge_boundary[e] || !t.edge_manifold[e] || !t.edge_consistent[e];
    const auto hash=detail::cut_hash(effective);
    {
        std::lock_guard lock(mutex_);
        for(auto it=entries_.begin();it!=entries_.end();++it)
            if(it->mesh==mesh.get() && it->hash==hash && it->cut==effective) {
                entries_.splice(entries_.begin(),entries_,it);++counters().hits;return it->surface;
            }
    }
    ++counters().misses;
    auto surface=std::make_shared<const CutSurface>(mesh,cut,memory_budget);
    const auto size=surface_bytes(*surface)+effective.size();
    std::lock_guard lock(mutex_);
    if(size<=budget_) {
        entries_.push_front({mesh.get(),hash,std::move(effective),surface,nullptr,size});bytes_+=size;evict_locked();peak_=std::max(peak_,bytes_);
    }
    return surface;
}
std::shared_ptr<const QualityGeometry> TopologyCache::geometry(const std::shared_ptr<const CutSurface>& cs,std::uint64_t memory_budget) {
    {
        std::lock_guard lock(mutex_);
        for(auto& entry:entries_)if(entry.surface==cs && entry.geometry) {++counters().hits;return entry.geometry;}
    }
    ++counters().misses;
    auto geometry=std::make_shared<const QualityGeometry>(cs,memory_budget);
    std::lock_guard lock(mutex_);
    for(auto it=entries_.begin();it!=entries_.end();++it)
        if(it->surface==cs) {
            if(!it->geometry) {const auto extra=geometry_bytes(*geometry);it->geometry=geometry;it->bytes+=extra;bytes_+=extra;}
            entries_.splice(entries_.begin(),entries_,it);evict_locked();peak_=std::max(peak_,bytes_);break;
        }
    return geometry;
}
CacheScope::CacheScope(TopologyCache* cache):previous_(active) {active=cache;}
CacheScope::~CacheScope() {active=previous_;}
namespace {thread_local std::uint64_t run_budget=default_memory_budget;}
std::uint64_t active_memory_budget() {return run_budget;}
MemoryBudgetScope::MemoryBudgetScope(std::uint64_t budget):previous_(run_budget) {run_budget=budget;}
MemoryBudgetScope::~MemoryBudgetScope() {run_budget=previous_;}
std::shared_ptr<const CutSurface> cut_surface(const std::shared_ptr<const MeshSnapshot>& mesh,std::span<const std::uint8_t> cut,std::uint64_t budget) {
    return active?active->surface(mesh,cut,budget):std::make_shared<const CutSurface>(mesh,cut,budget);
}
std::shared_ptr<const QualityGeometry> quality_geometry(const std::shared_ptr<const CutSurface>& cs,std::uint64_t budget) {
    return active?active->geometry(cs,budget):std::make_shared<const QualityGeometry>(cs,budget);
}
TopologyCache* active_cache() {return active;}
}
