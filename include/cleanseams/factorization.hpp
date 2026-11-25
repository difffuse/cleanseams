// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
#include "cleanseams/sparse_system.hpp"
#include <memory>
#include <optional>
#include <string>
namespace cleanseams {
enum class NumericalPolicy { scipy_compatible, numpy_compatible };
enum class FactorBackend { empty, superlu, dense_inverse, jacobi_cg, dense_lu, dense_least_squares };
enum class FactorSystem { laplacian, least_squares_normal };
struct FactorOptions {
    FactorSystem system=FactorSystem::laplacian;
    NumericalPolicy policy=NumericalPolicy::scipy_compatible;
    std::uint64_t memory_budget=active_memory_budget();
};
struct FactorMetadata {
    FactorSystem system=FactorSystem::laplacian;
    std::string fallback_reason;
    NumericalPolicy policy{};
    FactorBackend backend=FactorBackend::empty;
    std::string version;
    double regularization=0;
    std::optional<double> reciprocal_condition;
    std::uint64_t factorizations=0,solve_calls=0,peak_bytes=0;
};
struct ColumnSolveDiagnostics {
    std::uint64_t iterations=0;
    double relative_residual=0,original_relative_residual=0;
    CgTermination termination=CgTermination::converged;
};
struct FactorSolveResult {
    std::vector<double> solution; // row-major (unknowns, rhs_columns)
    std::vector<ColumnSolveDiagnostics> columns;
};
// Owns matrix and reusable factors. Independent instances can run concurrently;
// calls on a single instance are serialized by its chart/session owner.
class CLEANSEAMS_CORE_EXPORT SpdFactor {
public:
    explicit SpdFactor(const SparseMatrix&,const FactorOptions& options={});
    ~SpdFactor();
    SpdFactor(const SpdFactor&)=delete;
    SpdFactor& operator=(const SpdFactor&)=delete;
    FactorSolveResult solve(std::span<const double> rhs,Offset rhs_columns);
    const FactorMetadata& metadata() const;
private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};
}
