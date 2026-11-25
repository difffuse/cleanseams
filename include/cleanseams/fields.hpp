// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
#include "cleanseams/control.hpp"
#include "cleanseams/cut_surface.hpp"
#include <map>
#include <mutex>
namespace cleanseams {
// ------------------------------------------------------------------ features
// core/features.py. Masks are per edge (0/1). Chains are feature networks in the
// source's order: first appearance of a component, edges ascending.
struct FeatureOptions {
    double angle_deg=60.0;
    std::int64_t min_chain=6,rings=2;
    double ratio=0.6;
    std::int64_t max_extend=8;
};
struct TerminatedChains {std::vector<std::uint8_t> mask;std::vector<MeshIndex> open_edges;};
struct BevelReassignment {std::vector<std::uint8_t> cut;std::vector<MeshIndex> changed_charts,panel_faces;};
struct FeatureField {
    std::vector<std::uint8_t> cut;
    std::vector<std::vector<MeshIndex>> chains;
    std::vector<MeshIndex> open_features,panel_faces;
    std::vector<double> wide_dihedral; // strength; empty when detection did not run
    bool measured=false;
    // Stage records retained for diagnostics and parity tests: parallel suppression
    // output (before complete rims are restored) and terminated chains.
    std::vector<std::uint8_t> thinned,terminated;
    std::vector<MeshIndex> terminated_open,reassigned_charts;
    bool wide_measurement=false; // true when broad fitting ran (not n-gon-dominated meshes)
};
CLEANSEAMS_CORE_EXPORT std::vector<double> wide_dihedral(const MeshSnapshot&,std::int64_t rings=2,std::optional<double> radius=std::nullopt,
    std::optional<std::span<const std::uint8_t>> candidates=std::nullopt,const Deadline& deadline={});
CLEANSEAMS_CORE_EXPORT std::vector<MeshIndex> parallel_edges(const MeshSnapshot&,MeshIndex edge,std::int64_t steps);
CLEANSEAMS_CORE_EXPORT std::vector<std::uint8_t> suppress_parallel(const MeshSnapshot&,std::span<const std::uint8_t> candidates,
    std::span<const double> strength,double ratio=0.6,std::int64_t rings=2);
CLEANSEAMS_CORE_EXPORT std::vector<std::vector<MeshIndex>> feature_chains(const MeshSnapshot&,std::span<const std::uint8_t> mask);
CLEANSEAMS_CORE_EXPORT TerminatedChains terminate_chains(const MeshSnapshot&,std::span<const std::uint8_t> mask,std::int64_t min_chain=6,
    std::int64_t max_extend=8,std::optional<std::span<const std::uint8_t>> barriers=std::nullopt);
CLEANSEAMS_CORE_EXPORT BevelReassignment reassign_bevel_bands(const std::shared_ptr<const MeshSnapshot>&,std::span<const std::uint8_t> cut,
    std::optional<std::span<const std::uint8_t>> barriers=std::nullopt,bool allow_holes=false);
CLEANSEAMS_CORE_EXPORT FeatureField detect_features(const std::shared_ptr<const MeshSnapshot>&,const FeatureOptions& options={},
    std::optional<std::span<const std::uint8_t>> barriers=std::nullopt,const Deadline& deadline={});
// ------------------------------------------------------------------ exposure
// core/exposure.py: orthographic z-buffer visibility over Fibonacci directions.
CLEANSEAMS_CORE_EXPORT std::vector<Vec3> fibonacci_directions(std::uint64_t count);
struct ScreenFrames {std::vector<Vec3> u,w;};
CLEANSEAMS_CORE_EXPORT ScreenFrames screen_frames(std::span<const Vec3> directions);
struct SurfaceSamples {
    std::vector<Vec3> points;
    std::vector<TransferIndex> face_a,face_b,vertex; // vertex is -1 away from corners
};
CLEANSEAMS_CORE_EXPORT SurfaceSamples surface_samples(const MeshSnapshot&,double spacing,std::uint64_t memory_budget=active_memory_budget());
struct VisibilityOptions {
    double cell_factor=1.0,base=0.25,tan_max=5.7,spacing=0.5;
    std::uint64_t memory_budget=active_memory_budget();
    auto operator<=>(const VisibilityOptions&) const=default;
};
// Row-major (element, direction) visibility with float32 cosine weights, as the source.
struct Visibility {
    std::vector<Vec3> directions;
    std::vector<std::uint8_t> vertex,face;
    std::vector<float> vertex_weight,face_weight;
};
CLEANSEAMS_CORE_EXPORT Visibility zbuffer_visibility(const MeshSnapshot&,std::span<const Vec3> directions,
    const VisibilityOptions& options={},const SurfaceSamples* samples=nullptr);
// Raw exposure per vertex/face/edge; independent of any seam desirability discount.
struct ExposureField {
    std::vector<Vec3> directions;
    std::vector<double> vertex,face,edge;
    std::shared_ptr<const Visibility> visibility;
    std::string geometry_key;
};
CLEANSEAMS_CORE_EXPORT ExposureField exposure_field(const MeshSnapshot&,std::uint64_t count=64,
    std::optional<std::span<const Vec3>> directions=std::nullopt,const VisibilityOptions& options={});
// Bounded, session-owned reuse of geometry-only exposure. The key is the exact
// geometry content key plus directions and sampling options; an edited mesh
// never matches an old entry.
class CLEANSEAMS_CORE_EXPORT ExposureCache {
public:
    explicit ExposureCache(std::size_t max_entries=8):max_entries_(max_entries) {}
    std::shared_ptr<const ExposureField> get(const MeshSnapshot&,std::uint64_t count=64,const VisibilityOptions& options={});
    std::uint64_t hits() const {return hits_;}
    std::uint64_t misses() const {return misses_;}
    std::size_t size() const;
private:
    struct Key {std::string geometry;std::uint64_t count;VisibilityOptions options;auto operator<=>(const Key&) const=default;};
    mutable std::mutex mutex_;
    std::map<Key,std::shared_ptr<const ExposureField>> entries_;
    std::vector<Key> order_;
    std::size_t max_entries_;
    std::uint64_t hits_=0,misses_=0;
};
CLEANSEAMS_CORE_EXPORT double mesh_exposure(const MeshSnapshot&,const ExposureField&);
struct SeamExposure {double seam_exposure=0,mesh_exposure=0,exposure_ratio=0,hidden_seam_pct=0;};
CLEANSEAMS_CORE_EXPORT SeamExposure seam_exposure(const MeshSnapshot&,std::span<const std::uint8_t> seams,
    const ExposureField&,double hidden_below=0.3);
// Signed axis ("+X".."-Z") or finite nonzero vector, normalized. There is no default axis.
CLEANSEAMS_CORE_EXPORT Vec3 view_direction(std::string_view axis);
CLEANSEAMS_CORE_EXPORT Vec3 view_direction(const Vec3& vector);
CLEANSEAMS_CORE_EXPORT double seam_on_front_pct(const MeshSnapshot&,std::span<const std::uint8_t> seams,std::optional<Vec3> front);
struct SeamCostOptions {
    double feature_angle=60.0,valley_angle=20.0,valley_factor=0.5,hard_factor=0.1;
    std::optional<Vec3> hide_dir; // explicit only; never an implicit front axis
    double hide_strength=1.0;
};
// Legacy composed seam cost with its separately retained raw measurements.
struct SeamCostField {
    std::vector<double> exposure;          // raw edge exposure
    std::vector<std::uint8_t> valley,hard; // discount masks
    std::vector<double> facing;            // empty without hide_dir
    std::vector<double> cost;              // legacy composed, clipped to [0, 1]
};
CLEANSEAMS_CORE_EXPORT SeamCostField seam_cost(const MeshSnapshot&,const ExposureField&,const SeamCostOptions& options={});
}
