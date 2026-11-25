// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
#include "cleanseams/quality.hpp"
#include "cleanseams/factorization.hpp"
#include <list>
#include <map>
#include <mutex>
namespace cleanseams {
// Repeated-work counters (T21). Every CutSurface/QualityGeometry construction is
// recorded with its exact identity (mesh + effective cut); building one whose identity
// was already built during the counting window is an unchanged rebuild. Counting is
// per process and cheap; reset it around the work being measured. A window also ends
// after 65536 distinct identities, so memory stays bounded when nothing resets it.
struct WorkCounters {
    std::uint64_t cut_builds=0,unchanged_topology_rebuilds=0,geometry_builds=0,duplicate_geometry_builds=0,
                  cache_hits=0,cache_misses=0,cache_evictions=0,chart_solves=0,solve_cache_hits=0;
};
CLEANSEAMS_CORE_EXPORT WorkCounters work_counters();
CLEANSEAMS_CORE_EXPORT void reset_work_counters();
// Bounded cache of cut surfaces and their quality geometry, keyed by the exact
// effective cut of one mesh snapshot (hash lookup, then full equality). Owned by a
// pipeline run or session, never global. Eviction is least-recently-used by bytes.
class CLEANSEAMS_CORE_EXPORT TopologyCache {
public:
    explicit TopologyCache(std::uint64_t budget_bytes);
    // A miss builds with `memory_budget`; a hit returns the retained object.
    std::shared_ptr<const CutSurface> surface(const std::shared_ptr<const MeshSnapshot>&,std::span<const std::uint8_t> cut,
        std::uint64_t memory_budget=active_memory_budget());
    std::shared_ptr<const QualityGeometry> geometry(const std::shared_ptr<const CutSurface>&,std::uint64_t memory_budget=active_memory_budget());
    // Solved chart maps (T22) under their complete numerical contract: geometry, chart
    // key, method, iterations, weights, policy, pin completion, in-chart pin
    // coordinates, start map and validation. Anything else is a miss.
    // Sector ids are global to one cut while the contract keys the chart, so the sectors
    // (in sector_uv order) and the automatic pins are kept as one corner of each and are
    // relabelled for the surface that asks.
    struct SolvedMap {
        std::string method;std::vector<MeshIndex> corner_ids;std::vector<Vec2> corner_uv,sector_uv;
        FactorMetadata numerics;ColumnSolveDiagnostics diagnostics;std::vector<ColumnSolveDiagnostics> solve_diagnostics;
        std::uint64_t iterations=0;double equation_relative_residual=0;bool valid=false;std::optional<QualityReport> quality;
        std::vector<MeshIndex> sector_corners,pin_corners;
        // Bytes charged to the cache for this map under `contract`.
        CLEANSEAMS_CORE_EXPORT std::uint64_t bytes(const std::string& contract) const;
    };
    std::optional<SolvedMap> solved(const std::string& contract);
    void store(const std::string& contract,SolvedMap map);
    std::uint64_t bytes() const;
    std::uint64_t peak_bytes() const;
    std::uint64_t budget() const {return budget_;}
private:
    struct Entry {
        const MeshSnapshot* mesh=nullptr;std::uint64_t hash=0;std::vector<std::uint8_t> cut;
        std::shared_ptr<const CutSurface> surface;std::shared_ptr<const QualityGeometry> geometry;std::uint64_t bytes=0;
    };
    mutable std::mutex mutex_;
    std::list<Entry> entries_; // most recent first
    std::list<std::pair<std::string,SolvedMap>> maps_; // most recent first
    std::map<std::string,std::list<std::pair<std::string,SolvedMap>>::iterator> map_index_;
    std::uint64_t map_bytes_=0;
    std::uint64_t budget_=0,bytes_=0,peak_=0;
    void evict_locked();
};
// Makes `cache` the calling thread's active cache for the scope (nesting restores).
class CLEANSEAMS_CORE_EXPORT CacheScope {
public:
    explicit CacheScope(TopologyCache* cache);
    ~CacheScope();
    CacheScope(const CacheScope&)=delete;
    CacheScope& operator=(const CacheScope&)=delete;
private:
    TopologyCache* previous_=nullptr;
};
// Construction through the active cache when there is one, otherwise a plain build.
CLEANSEAMS_CORE_EXPORT std::shared_ptr<const CutSurface> cut_surface(const std::shared_ptr<const MeshSnapshot>&,std::span<const std::uint8_t> cut,
    std::uint64_t memory_budget=active_memory_budget());
CLEANSEAMS_CORE_EXPORT std::shared_ptr<const QualityGeometry> quality_geometry(const std::shared_ptr<const CutSurface>&,
    std::uint64_t memory_budget=active_memory_budget());
CLEANSEAMS_CORE_EXPORT TopologyCache* active_cache();
namespace detail {
void note_cut_build(const MeshSnapshot&,std::span<const std::uint8_t> effective_cut);
void note_geometry_build(const CutSurface&);
void note_chart_solve();
std::uint64_t cut_hash(std::span<const std::uint8_t> cut);
}
}
