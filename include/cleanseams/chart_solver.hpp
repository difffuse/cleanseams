// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
#include "cleanseams/quality.hpp"
#include "cleanseams/factorization.hpp"
#include "cleanseams/control.hpp"
namespace cleanseams {
// PIN-01: legacy_010 reproduces the source,
// which replaces a single in-chart user pin by the automatic pair; preserve keeps
// that pin and adds the farthest distinct anchor at its 3D distance along +U.
enum class PinCompletion { legacy_010, preserve };
enum class HarmonicWeights { mean_value, uniform };
enum class ChartMethod { lscm, harmonic, arap };
// A constraint on a global cut sector, as the source's {sector: (u, v)} pins.
struct SectorPin {MeshIndex sector=0;Vec2 uv{};};
struct ChartSolveRequest {
    std::shared_ptr<const QualityGeometry> geometry;
    TransferIndex chart=0;
    // Corner and sector pins resolve to this chart's cut sectors; pins elsewhere
    // are ignored. Identical repeats are one constraint; different coordinates on
    // one sector are a constraint_conflict under preserve, never averaged or
    // overwritten; legacy_010 keeps the last one written, as the source (Gate A F7).
    std::vector<CornerPin> pins;
    std::vector<SectorPin> sector_pins;
    NumericalPolicy policy=NumericalPolicy::scipy_compatible;
    PinCompletion pin_completion=PinCompletion::preserve;
    HarmonicWeights weights=HarmonicWeights::mean_value;
    std::uint64_t iterations=10; // ARAP local/global steps: source 6 generation, 10 final
    std::vector<Vec2> initial_sector_uv; // ARAP start in local sector order; empty uses LSCM
    bool validate=true; // full chart validity and distortion evaluation of the map
    std::uint64_t memory_budget=active_memory_budget();
    bool retain_assembly=false;
};
struct ChartTriangulation {
    std::vector<std::array<MeshIndex,3>> triangles,triangle_corners;
    std::vector<MeshIndex> triangle_ids,sectors,corner_ids;
};
struct LscmAssembly {
    ChartTriangulation chart;
    SparseMatrix matrix,reduced;
    std::vector<double> fixed_values,rhs;
    std::vector<TransferIndex> free_index;
    std::vector<MeshIndex> automatic_pins; // global sectors fixed by completion, not the user
    std::uint64_t storage_bound=0;
};
// Harmonic/ARAP system in source order, retained only on request.
struct LaplacianTrace {
    std::vector<std::array<MeshIndex,2>> weight_pairs; // harmonic weights, dict insertion order
    std::vector<double> weights;
    std::vector<std::array<double,3>> cotangents; // ARAP, per triangle corner after clamping
    std::unique_ptr<SparseMatrix> matrix,reduced;
    std::vector<std::uint8_t> fixed;
    std::vector<TransferIndex> free_index;
    std::vector<Vec2> initial_sector_uv; // ARAP start before pins are reapplied
    std::vector<std::vector<double>> rhs,solutions; // each factored solve, row-major (free, 2)
};
struct LinearLeastSquaresResult {
    std::vector<double> solution;
    FactorMetadata numerics;
    ColumnSolveDiagnostics diagnostics;
    double equation_relative_residual=0;
    std::unique_ptr<SparseMatrix> normal; // regularized, retained only on request
    std::vector<double> normal_rhs;
};
struct ChartSolveResult {
    std::string method; // lscm, harmonic, arap, convex or planar
    std::vector<MeshIndex> corner_ids,sectors,automatic_pins;
    std::vector<Vec2> corner_uv,sector_uv;
    FactorMetadata numerics;
    ColumnSolveDiagnostics diagnostics; // last factored solve, first column
    std::vector<ColumnSolveDiagnostics> solve_diagnostics; // every factored solve and column
    std::uint64_t iterations=0; // ARAP steps completed
    double equation_relative_residual=0;
    bool valid=false; // measured by the chart validator; never inferred from finiteness
    std::optional<QualityReport> quality;
    std::unique_ptr<LscmAssembly> assembly;
    std::unique_ptr<LaplacianTrace> laplacian;
    std::unique_ptr<SparseMatrix> normal;
    std::vector<double> normal_rhs;
};
// Source chart_triangles: the chart's triangles of the shared fixed triangulation,
// its sectors sorted by global ID and the triangulated corners.
CLEANSEAMS_CORE_EXPORT ChartTriangulation chart_triangulation(const QualityGeometry&,MeshIndex chart);
CLEANSEAMS_CORE_EXPORT LscmAssembly assemble_lscm(const ChartSolveRequest&);
CLEANSEAMS_CORE_EXPORT LinearLeastSquaresResult solve_least_squares(const SparseMatrix&,std::span<const double> rhs,
    NumericalPolicy policy=NumericalPolicy::scipy_compatible,std::uint64_t memory_budget=active_memory_budget(),bool retain_system=false);
CLEANSEAMS_CORE_EXPORT ChartSolveResult solve_lscm(const ChartSolveRequest&);
CLEANSEAMS_CORE_EXPORT ChartSolveResult solve_harmonic(const ChartSolveRequest&);
// One factorization of the free Laplacian serves every iteration's two-column solve.
CLEANSEAMS_CORE_EXPORT ChartSolveResult solve_arap(const ChartSolveRequest&);
// Coarse positive-weight fallback for single-loop Euler-one unpinched charts; empty
// when the source declines the chart. It cannot honor pins: in-chart pins are a
// constraint_conflict (the source's callers never pass pinned charts).
CLEANSEAMS_CORE_EXPORT std::optional<ChartSolveResult> convex_boundary_map(const ChartSolveRequest&);
// Orthographic map of a near-planar disk containing an n-gon, if it passes the
// source's geometry and quality gates. In-chart pins are a constraint_conflict.
CLEANSEAMS_CORE_EXPORT std::optional<ChartSolveResult> planar_ngon_map(const ChartSolveRequest&);
// Source flatten_chart. For ARAP without in-chart pins a valid planar map is
// refined by ARAP and kept when the refinement flips, folds or crosses.
CLEANSEAMS_CORE_EXPORT ChartSolveResult flatten_chart(const ChartSolveRequest&,ChartMethod method=ChartMethod::arap);
struct FlattenAllRequest {
    std::shared_ptr<const QualityGeometry> geometry;
    std::optional<std::vector<MeshIndex>> charts; // all charts when empty
    std::vector<Vec2> uv; // entries for other charts are kept; zeros when empty
    std::vector<CornerPin> pins;
    ChartMethod method=ChartMethod::arap;
    NumericalPolicy policy=NumericalPolicy::scipy_compatible;
    PinCompletion pin_completion=PinCompletion::preserve;
    std::uint64_t iterations=10;
    bool validate=true;
    std::uint64_t memory_budget=active_memory_budget();
    CancellationToken cancel{}; // charts not yet started are skipped; the call then throws cancelled
};
struct FlattenAllResult {
    std::vector<Vec2> uv; // original corners
    std::vector<MeshIndex> solved_charts;
    std::vector<ChartSolveResult> charts; // parallel to solved_charts
};
CLEANSEAMS_CORE_EXPORT FlattenAllResult flatten_all(const FlattenAllRequest&);
// Native replacement of the source's CoreBackend.unwrap (the solver used when no
// Blender is present): cut topology, then flatten_all at the default ten
// iterations without a start map, or per-chart ARAP with the requested
// iteration count and optional start map. Blender remains a separate host backend.
class CLEANSEAMS_CORE_EXPORT NativeArapBackend {
public:
    struct Settings {
        std::uint64_t iterations=10;
        NumericalPolicy policy=NumericalPolicy::scipy_compatible;
        PinCompletion pin_completion=PinCompletion::preserve;
        bool validate=false; // callers measure the assembled map with their own gates
        std::uint64_t memory_budget=active_memory_budget();
        CancellationToken cancel{}; // checked after every chart solve (T24)
    };
    struct Unwrap {
        std::shared_ptr<const CutSurface> cut;
        std::shared_ptr<const QualityGeometry> geometry;
        std::vector<Vec2> uv;
        std::vector<MeshIndex> solved_charts;
        std::vector<ChartSolveResult> charts;
    };
    explicit NativeArapBackend(std::shared_ptr<const MeshSnapshot>,Settings settings);
    explicit NativeArapBackend(std::shared_ptr<const MeshSnapshot> mesh):NativeArapBackend(std::move(mesh),Settings{}) {}
    const Settings& settings() const {return settings_;}
    Unwrap unwrap(std::span<const std::uint8_t> cut,std::optional<std::span<const MeshIndex>> charts=std::nullopt,
        std::span<const CornerPin> pins={},std::optional<std::uint64_t> iterations=std::nullopt,
        std::optional<std::span<const Vec2>> init=std::nullopt) const;
private:
    std::shared_ptr<const MeshSnapshot> mesh_;
    Settings settings_;
};
}
