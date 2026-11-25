// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
// Phase B (T33, Section 14.3): measured seam-end necessity for enhanced_v1. A dead-end seam
// (a relief slit or any other spur) is removed only in a trial state that goes through the
// T28 transaction (topology, solve, validity, strict overlap, pins, budgets, envelope, common
// ordering); if the removal breaks the budget or validity the end is kept and the measured
// reason recorded. A relief tag grants no immunity: the evidence belongs to one state and is
// invalidated when topology, solver, pins, quality policy, geometry or the chart's map change.
#include "cleanseams/proposals.hpp"
namespace cleanseams {
inline constexpr const char* seam_necessity_version="seam_necessity.v1";
enum class SeamEndClass : std::uint8_t {
    necessary_relief_end=0, // removal trial refused or ranked after the state for a measured quality/validity reason
    removable_spur=1,       // removal trial admissible and ranking strictly before the state
    protected_end=2,        // a hard constraint keeps it (locked/required edge, fixed anchor, pins, topology fact)
    unmeasured=3            // no trial: budget exhausted, deadline or cancellation
};
CLEANSEAMS_CORE_EXPORT std::string_view seam_end_class_name(SeamEndClass);
struct CLEANSEAMS_CORE_EXPORT SeamEndEvidence {
    DeadEndArc arc;
    std::string key;                          // necessity_key: the evidence's identity
    SeamEndClass classification=SeamEndClass::unmeasured;
    std::string reason;                       // measured reason (a rejection, "over_budget", "objective", "removal_ranks_first", ...)
    bool relief=false;                        // the arc was added by a committed hotspot_relief proposal (provenance only)
    bool cached=false;                        // served from evidence with the same key
    nlohmann::json measured=nlohmann::json::object(); // the removal trial and the comparison with the state
    nlohmann::json to_json() const;
};
struct CLEANSEAMS_CORE_EXPORT SeamNecessity {
    std::string version=seam_necessity_version;
    Revision revision=0;                      // the classified state
    std::string cut_digest;
    std::uint64_t raw_dangling_ends=0;        // legacy dangling_seam_ends of the state's seams (diagnostic)
    std::uint64_t necessary_relief_ends=0,removable_spur_ends=0,protected_ends=0,unmeasured_ends=0; // per dangling end
    std::uint64_t trials=0,cache_hits=0;
    std::vector<SeamEndEvidence> arcs;
    nlohmann::json to_json() const;
};
// Evidence identity of one dead-end arc of `state`: the mesh geometry, the arc and its chart's
// topology (chart key), the solver identity (policy, solver record, trial initialization), the
// pins of the chart, the quality policy (budgets, caps, gate tolerance, envelope growth,
// objective, enhanced policy version, exemptions of the chart), the chart's map and (GB m9)
// the run-origin envelope domain of the removal trial: the origin revision, the relief
// rebases the state adopted and the key and map of every other chart that differs from the
// origin; for an over-budget state also its worst ratio and total excess (the ordering).
CLEANSEAMS_CORE_EXPORT nlohmann::json necessity_key_parts(const ProposalContext&,const RunState&,const DeadEndArc&,const std::string& initialization="");
CLEANSEAMS_CORE_EXPORT std::string necessity_key(const ProposalContext&,const RunState&,const DeadEndArc&,const std::string& initialization="");
// Evidence by key. A lookup under any other key misses, so stale evidence is never used.
class CLEANSEAMS_CORE_EXPORT SeamNecessityCache {
public:
    const SeamEndEvidence* find(const std::string& key) const;
    void store(const SeamEndEvidence&);
    std::uint64_t hits() const {return hits_;}
    std::uint64_t misses() const {return misses_;}
    std::size_t size() const {return entries_.size();}
private:
    std::map<std::string,SeamEndEvidence> entries_;
    mutable std::uint64_t hits_=0,misses_=0;
};
// Classifies every dead-end arc of `state` by a trial removal through the transaction (or
// the cache's evidence under the same key). Nothing is committed.
CLEANSEAMS_CORE_EXPORT SeamNecessity measure_seam_necessity(const RunState& state,ProposalContext&,SeamNecessityCache&,const std::string& initialization="");
// Section 21 "validated unnecessary-spur removal": `rounds` generators sharing the cache
// (ProposalContext::necessity when set). Each measures the incumbent's dead ends and proposes
// the removals that rank strictly before it (and all of them together); the round commits the
// best. Removals whose trial was refused on the same state are never proposed again.
// layout.spur_removal false keeps every end (measurement only).
CLEANSEAMS_CORE_EXPORT std::vector<ProposalGenerator> spur_removal_generators(const PipelineOptions&,std::int64_t rounds=3);
}
