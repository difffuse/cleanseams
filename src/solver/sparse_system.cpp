// SPDX-License-Identifier: GPL-3.0-or-later
#include "cleanseams/sparse_system.hpp"
#include "cleanseams/errors.hpp"
#include "sparse_order.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <tuple>
namespace cleanseams {
namespace {
// GB performance: the check inlines into the per-entry validation loops; only the failure
// path (an out-of-line, cold call) builds the error.
[[noreturn,gnu::noinline,gnu::cold]] void invalid(const char* field,const char* message) {
    throw ContractError(ErrorCode::invalid_input,field,message);
}
inline void require(bool condition,const char* field,const char* message) {
    if(!condition) [[unlikely]] invalid(field,message);
}
void storage(bool condition) {
    if(!condition) throw ContractError(ErrorCode::memory_limit,"solver.sparse","sparse storage exceeds memory allowance");
}
void finite_result(std::span<const double> values) {
    for(auto x:values) if(!std::isfinite(x))
        throw ContractError(ErrorCode::numerical_failure,"solver.sparse","nonfinite arithmetic result");
}
void vector_input(std::span<const double> values,Offset expected) {
    require(values.size()==expected,"solver.sparse.vector","vector length mismatch");
    for(auto x:values) require(std::isfinite(x),"solver.sparse.vector","vector must be finite");
}
}
SparseMatrix::SparseMatrix(std::span<const TransferIndex> rows,std::span<const TransferIndex> columns,
    std::span<const double> values,std::array<Offset,2> shape,std::uint64_t budget)
    :shape_(shape),budget_(std::min<std::uint64_t>(budget,std::numeric_limits<std::size_t>::max())) {
    require(rows.size()==columns.size() && rows.size()==values.size(),"solver.sparse.entries","COO lengths must match");
    // Even a sparse extent must be representable by this operation's vector
    // outputs. Check signed indexing and byte arithmetic before any allocation.
    for(auto n:shape) {
        storage(n<=budget_/sizeof(double));
        storage(n<static_cast<std::uint64_t>(std::numeric_limits<std::ptrdiff_t>::max())/sizeof(double));
    }
    constexpr auto stride=2*sizeof(TransferIndex)+sizeof(double);
    storage(values.size()<=budget_/stride);owned_bytes_=values.size()*stride;
    for(std::size_t i=0;i<values.size();++i) {
        require(rows[i]>=0 && static_cast<Offset>(rows[i])<shape[0],"solver.sparse.rows","row out of range");
        require(columns[i]>=0 && static_cast<Offset>(columns[i])<shape[1],"solver.sparse.columns","column out of range");
        require(std::isfinite(values[i]),"solver.sparse.values","coefficient must be finite");
    }
    rows_.assign(rows.begin(),rows.end());columns_.assign(columns.begin(),columns.end());values_.assign(values.begin(),values.end());
}
void SparseMatrix::output_allowance(std::uint64_t count,std::uint64_t stride) const {
    storage(count<=(budget_-owned_bytes_)/stride);
}
std::vector<double> SparseMatrix::matvec(std::span<const double> x) const {
    vector_input(x,shape_[1]);output_allowance(shape_[0],sizeof(double));
    std::vector<double> out(static_cast<std::size_t>(shape_[0]),0);
    for(std::size_t i=0;i<values_.size();++i) out[static_cast<std::size_t>(rows_[i])]+=values_[i]*x[static_cast<std::size_t>(columns_[i])];
    finite_result(out);return out;
}
std::vector<double> SparseMatrix::rmatvec(std::span<const double> x) const {
    vector_input(x,shape_[0]);output_allowance(shape_[1],sizeof(double));
    std::vector<double> out(static_cast<std::size_t>(shape_[1]),0);
    for(std::size_t i=0;i<values_.size();++i) out[static_cast<std::size_t>(columns_[i])]+=values_[i]*x[static_cast<std::size_t>(rows_[i])];
    finite_result(out);return out;
}
std::vector<double> SparseMatrix::dense() const {
    const auto available=(budget_-owned_bytes_)/sizeof(double);
    storage(shape_[0]==0 || shape_[1]<=available/shape_[0]);
    std::vector<double> out(static_cast<std::size_t>(shape_[0]*shape_[1]),0);
    for(std::size_t i=0;i<values_.size();++i)
        out[static_cast<std::size_t>(static_cast<Offset>(rows_[i])*shape_[1]+static_cast<Offset>(columns_[i]))]+=values_[i];
    finite_result(out);return out;
}
std::vector<double> SparseMatrix::diag_of_normal() const {
    output_allowance(shape_[1],sizeof(double));std::vector<double> out(static_cast<std::size_t>(shape_[1]),0);
    // Match source Sparse.diag_of_normal, including duplicate COO semantics.
    for(std::size_t i=0;i<values_.size();++i) out[static_cast<std::size_t>(columns_[i])]+=values_[i]*values_[i];
    finite_result(out);return out;
}
CompressedRows SparseMatrix::csr() const {
    output_allowance(shape_[0]+1,sizeof(Offset));
    const auto offsets_bytes=(shape_[0]+1)*sizeof(Offset);
    // Permutation scratch plus worst-case output entries, with no growth copies.
    storage(values_.size()<=(budget_-owned_bytes_-offsets_bytes)/(sizeof(std::size_t)+sizeof(TransferIndex)+sizeof(double)));
    // Entries ordered by (row, column, original index): duplicate summation retains the
    // original COO order. Linear-time counting passes give the same permutation as the
    // comparison sort with the index tie-break (GB performance, sparse_order.hpp).
    const auto order=solver_detail::lexicographic_order(std::span<const TransferIndex>(rows_),shape_[0],std::span<const TransferIndex>(columns_),shape_[1]);
    CompressedRows out;out.shape=shape_;out.row_offsets.assign(static_cast<std::size_t>(shape_[0]+1),0);
    out.columns.reserve(values_.size());out.values.reserve(values_.size());
    for(std::size_t start=0;start<order.size();) {
        const auto first=order[start];auto end=start+1;double sum=values_[first];
        while(end<order.size() && rows_[order[end]]==rows_[first] && columns_[order[end]]==columns_[first]) sum+=values_[order[end++]];
        out.columns.push_back(columns_[first]);out.values.push_back(sum);
        ++out.row_offsets[static_cast<std::size_t>(rows_[first])+1];start=end;
    }
    std::partial_sum(out.row_offsets.begin(),out.row_offsets.end(),out.row_offsets.begin());
    finite_result(out.values);return out;
}
CgResult conjugate_gradient(const SparseMatrix& matrix,std::span<const double> rhs,
    const CgOptions& options,std::span<const double> initial,std::span<const double> diagonal) {
    require(std::isfinite(options.tolerance) && options.tolerance>=0,"solver.cg.tolerance","tolerance must be finite and nonnegative");
    require(options.normal_equations || matrix.shape()[0]==matrix.shape()[1],"solver.cg.matrix","direct CG requires a square matrix");
    vector_input(rhs,matrix.shape()[1]);
    if(!initial.empty()) vector_input(initial,rhs.size());
    if(!diagonal.empty()) vector_input(diagonal,rhs.size());
    const auto budget=std::min<std::uint64_t>(options.memory_budget,std::numeric_limits<std::size_t>::max());
    const auto n=rhs.size();storage(n<=budget/(8*sizeof(double)));
    const auto workspace_left=budget-n*8*sizeof(double);
    storage(!options.normal_equations || matrix.shape()[0]<=workspace_left/(2*sizeof(double)));
    auto product=[&](std::span<const double> x) {
        auto ax=matrix.matvec(x);return options.normal_equations?matrix.rmatvec(ax):std::move(ax);
    };
    auto dot=[](std::span<const double> a,std::span<const double> b) {
        double sum=0;for(std::size_t i=0;i<a.size();++i) sum+=a[i]*b[i];
        finite_result(std::array<double,1>{sum});return sum;
    };
    CgResult out;out.solution.assign(n,0);
    if(!initial.empty()) std::copy(initial.begin(),initial.end(),out.solution.begin());
    auto ap=product(out.solution);
    std::vector<double> r(n),inverse(n,1),z(n),p(n);
    for(std::size_t i=0;i<n;++i) {
        r[i]=rhs[i]-ap[i];
        if(!diagonal.empty()) inverse[i]=1/std::max(diagonal[i],1e-300);
        z[i]=inverse[i]*r[i];
    }
    finite_result(z);p=z;double rz=dot(r,z);
    double bnorm=std::sqrt(dot(rhs,rhs));if(bnorm==0) bnorm=1;
    out.relative_residual=std::sqrt(dot(r,r))/bnorm;
    while(out.relative_residual>options.tolerance && out.iterations<options.max_iterations) {
        ap=product(p);const double curvature=dot(p,ap);
        if(curvature<=0) {out.termination=CgTermination::non_positive_curvature;return out;}
        const double alpha=rz/curvature;
        for(std::size_t i=0;i<n;++i) {out.solution[i]+=alpha*p[i];r[i]-=alpha*ap[i];}
        finite_result(out.solution);out.relative_residual=std::sqrt(dot(r,r))/bnorm;
        for(std::size_t i=0;i<n;++i) z[i]=inverse[i]*r[i];
        const double next_rz=dot(r,z),beta=next_rz/rz;
        for(std::size_t i=0;i<n;++i) p[i]=z[i]+beta*p[i];
        finite_result(p);rz=next_rz;++out.iterations;
    }
    out.termination=out.relative_residual<=options.tolerance?CgTermination::converged:CgTermination::iteration_limit;
    return out;
}
}
