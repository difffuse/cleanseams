// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
#include "cleanseams/types.hpp"
#include <array>
namespace cleanseams {
struct CompressedRows {
    std::array<Offset,2> shape{};
    std::vector<Offset> row_offsets;
    std::vector<TransferIndex> columns;
    std::vector<double> values;
};
// Owned COO entries retain assembly order and duplicates. Compression merges
// duplicates, while matvec and the legacy normal preconditioner use COO order.
class CLEANSEAMS_CORE_EXPORT SparseMatrix {
public:
    SparseMatrix(std::span<const TransferIndex> rows,std::span<const TransferIndex> columns,
        std::span<const double> values,std::array<Offset,2> shape,
        std::uint64_t memory_budget=active_memory_budget());
    const std::array<Offset,2>& shape() const {return shape_;}
    std::span<const TransferIndex> rows() const {return rows_;}
    std::span<const TransferIndex> columns() const {return columns_;}
    std::span<const double> values() const {return values_;}
    std::vector<double> matvec(std::span<const double>) const;
    std::vector<double> rmatvec(std::span<const double>) const;
    std::vector<double> dense() const; // row-major
    std::vector<double> diag_of_normal() const;
    CompressedRows csr() const;
private:
    void output_allowance(std::uint64_t count,std::uint64_t stride) const;
    std::array<Offset,2> shape_;
    std::vector<TransferIndex> rows_,columns_;
    std::vector<double> values_;
    std::uint64_t budget_,owned_bytes_=0;
};
enum class CgTermination { converged, iteration_limit, non_positive_curvature };
struct CgOptions {
    double tolerance=1e-10;
    std::uint64_t max_iterations=5000;
    bool normal_equations=false; // apply A^T A without assembling it
    std::uint64_t memory_budget=active_memory_budget();
};
struct CgResult {
    std::vector<double> solution;
    std::uint64_t iterations=0;
    double relative_residual=0;
    CgTermination termination=CgTermination::converged;
};
CLEANSEAMS_CORE_EXPORT CgResult conjugate_gradient(const SparseMatrix&,std::span<const double> rhs,
    const CgOptions& options={},std::span<const double> initial={},std::span<const double> diagonal={});
}
