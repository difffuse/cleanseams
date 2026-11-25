// SPDX-License-Identifier: GPL-3.0-or-later
#include "cleanseams/factorization.hpp"
#include "cleanseams/errors.hpp"
#include "superlu_bridge.h"
#include "sparse_order.hpp"
#include <Eigen/LU>
#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <tuple>
namespace cleanseams {
namespace {
void require(bool yes,const char* message) {
    if(!yes)throw ContractError(ErrorCode::invalid_input,"solver.factor",message);
}
void memory(bool yes) {
    if(!yes)throw ContractError(ErrorCode::memory_limit,"solver.factor","factor or solve exceeds memory budget");
}
void finite(double value) {
    if(!std::isfinite(value))throw ContractError(ErrorCode::numerical_failure,"solver.factor","nonfinite factorization or solve arithmetic");
}
void lu_check(int result,const cs_lu* factor) {
    if(result!=CS_LU_OK)throw ContractError(result==CS_LU_MEMORY?ErrorCode::memory_limit:
        (result==CS_LU_ARGUMENT?ErrorCode::invalid_input:ErrorCode::numerical_failure),"solver.superlu",cs_lu_error(factor));
}
std::uint64_t bytes(std::uint64_t count,std::uint64_t stride,std::uint64_t budget) {
    memory(count<=budget/stride);return count*stride;
}
struct Entry {TransferIndex column,row;std::size_t order;double value;};
using LuOwner=std::unique_ptr<cs_lu,decltype(&cs_lu_destroy)>;
}
struct SpdFactor::Impl {
    FactorOptions options;
    FactorMetadata meta;
    std::unique_ptr<SparseMatrix> matrix;
    Eigen::MatrixXd inverse;
    std::unique_ptr<Eigen::PartialPivLU<Eigen::MatrixXd>> dense_lu;
    std::vector<double> diagonal;
    LuOwner lu{nullptr,cs_lu_destroy};
    std::uint64_t base_bytes=sizeof(Impl)+sizeof(SparseMatrix)+256,retained_bytes=0;
    std::size_t n=0;
    explicit Impl(const SparseMatrix& input,FactorOptions config):options(config) {
        require(input.shape()[0]==input.shape()[1],"matrix must be square");
        require(options.policy==NumericalPolicy::scipy_compatible || options.policy==NumericalPolicy::numpy_compatible,"unknown numerical policy");
        options.memory_budget=std::min<std::uint64_t>(options.memory_budget,std::numeric_limits<std::size_t>::max());
        memory(base_bytes<=options.memory_budget);
        base_bytes+=bytes(input.values().size(),24,options.memory_budget-base_bytes);
        n=static_cast<std::size_t>(input.shape()[0]);meta.policy=options.policy;meta.system=options.system;
        require(options.system==FactorSystem::laplacian || options.system==FactorSystem::least_squares_normal,"unknown factor system");
        const bool normal=options.system==FactorSystem::least_squares_normal;
        require(!normal || options.policy!=NumericalPolicy::numpy_compatible || n<=1800,"large NumPy least squares requires the matrix-free normal operator");
        matrix=std::make_unique<SparseMatrix>(input.rows(),input.columns(),input.values(),input.shape(),options.memory_budget);
        retained_bytes=base_bytes;meta.peak_bytes=base_bytes;
        if(!n) {meta.version="empty-v1";return;}
        if(options.policy==NumericalPolicy::scipy_compatible) {
            meta.backend=FactorBackend::superlu;meta.version="SuperLU-7.0.0/COLAMD/pivot=1/serial";meta.regularization=normal?1e-12:1e-10;
            sparse_factor();
        } else if(n<=(normal?1800U:900U)) {
            meta.backend=normal?FactorBackend::dense_lu:FactorBackend::dense_inverse;
            meta.version=normal?"Eigen-3.4.0/PartialPivLU/solve":"Eigen-3.4.0/PartialPivLU/inverse";meta.regularization=normal?1e-12:1e-9;
            dense_factor();
        } else {
            meta.backend=FactorBackend::jacobi_cg;meta.version="CleanSeams-CG-v1/tol=1e-8/maxiter=6000";
            retained_bytes+=bytes(n,sizeof(double),options.memory_budget-retained_bytes);
            diagonal.assign(n,0);
            for(std::size_t k=0;k<input.values().size();++k)if(input.rows()[k]==input.columns()[k])
                diagonal[static_cast<std::size_t>(input.rows()[k])]+=input.values()[k];
            for(double d:diagonal)finite(d);
            meta.peak_bytes=retained_bytes;
        }
    }
    void dense_factor() {
        // Simultaneous matrix, LU, inverse and conservative Eigen solve scratch.
        const auto matrix_bytes=bytes(n,n*sizeof(double),options.memory_budget-base_bytes);
        memory(matrix_bytes<=(options.memory_budget-base_bytes)/5);
        auto peak=base_bytes+5*matrix_bytes;
        peak+=bytes(n,64,options.memory_budget-peak);
        meta.peak_bytes=peak;
        Eigen::MatrixXd dense=Eigen::MatrixXd::Zero(static_cast<Eigen::Index>(n),static_cast<Eigen::Index>(n));
        for(std::size_t k=0;k<matrix->values().size();++k)dense(matrix->rows()[k],matrix->columns()[k])+=matrix->values()[k];
        dense.diagonal().array()+=meta.regularization;
        require(dense.rows()>0,"empty dense factor");
        if(!dense.allFinite())finite(std::numeric_limits<double>::infinity());
        Eigen::PartialPivLU<Eigen::MatrixXd> decomposition(dense);
        if((decomposition.matrixLU().diagonal().array()==0).any())
            throw ContractError(ErrorCode::numerical_failure,meta.backend==FactorBackend::dense_lu?"solver.dense_solve.singular":"solver.dense_inverse","singular pivot after regularization");
        meta.reciprocal_condition=decomposition.rcond();finite(*meta.reciprocal_condition);
        if(meta.backend==FactorBackend::dense_lu) dense_lu=std::make_unique<Eigen::PartialPivLU<Eigen::MatrixXd>>(std::move(decomposition));
        else {inverse=decomposition.inverse();if(!inverse.allFinite())finite(std::numeric_limits<double>::infinity());}
        retained_bytes=base_bytes+matrix_bytes+n*64;meta.factorizations=1;
    }
    void sparse_factor() {
        require(n<=static_cast<std::size_t>(std::numeric_limits<int>::max()),"SuperLU dimension exceeds index representation");
        const auto count=matrix->values().size()+n;
        require(count<=static_cast<std::size_t>(std::numeric_limits<int>::max()),"SuperLU nonzero count exceeds index representation");
        auto temporary=bytes(count,sizeof(Entry)+sizeof(int)+sizeof(double),options.memory_budget-base_bytes);
        temporary+=bytes(n+1,sizeof(int),options.memory_budget-base_bytes-temporary);
        // The matrix entries, then the regularization diagonal, ordered by (column, row,
        // entry index): duplicates merge in their original order. The linear-time counting
        // order is the same permutation as the comparison sort with that tie-break (GB
        // performance, sparse_order.hpp).
        const auto nnz=matrix->values().size();
        std::vector<TransferIndex> entry_rows(count),entry_columns(count);std::vector<double> entry_values(count);
        std::copy(matrix->rows().begin(),matrix->rows().end(),entry_rows.begin());
        std::copy(matrix->columns().begin(),matrix->columns().end(),entry_columns.begin());
        std::copy(matrix->values().begin(),matrix->values().end(),entry_values.begin());
        for(std::size_t i=0;i<n;++i) {entry_rows[nnz+i]=entry_columns[nnz+i]=static_cast<TransferIndex>(i);entry_values[nnz+i]=meta.regularization;}
        const auto order=solver_detail::lexicographic_order(std::span<const TransferIndex>(entry_columns),n,std::span<const TransferIndex>(entry_rows),n);
        std::vector<int> offsets(n+1,0),rows;std::vector<double> values;rows.reserve(count);values.reserve(count);
        double norm=0,column_sum=0;TransferIndex column=-1;
        for(std::size_t k=0;k<order.size();) {
            const auto head=order[k];const Entry first{entry_columns[head],entry_rows[head],head,entry_values[head]};double sum=first.value;++k;
            while(k<order.size() && entry_columns[order[k]]==first.column && entry_rows[order[k]]==first.row)sum+=entry_values[order[k++]];
            finite(sum);
            if(sum==0)continue; // scipy sparse addition removes exact zeros
            if(first.column!=column) {norm=std::max(norm,column_sum);column_sum=0;column=first.column;}
            column_sum+=std::abs(sum);finite(column_sum);
            rows.push_back(static_cast<int>(first.row));values.push_back(sum);++offsets[static_cast<std::size_t>(first.column)+1];
        }
        norm=std::max(norm,column_sum);
        std::partial_sum(offsets.begin(),offsets.end(),offsets.begin());
        const auto overhead=base_bytes+temporary;
        lu.reset(cs_lu_create(static_cast<std::size_t>(options.memory_budget-overhead)));memory(bool(lu));
        lu_check(cs_lu_factor(lu.get(),static_cast<int>(n),static_cast<int>(values.size()),offsets.data(),rows.data(),values.data(),norm),lu.get());
        meta.peak_bytes=overhead+cs_lu_peak(lu.get());
        meta.reciprocal_condition=cs_lu_rcond(lu.get());finite(*meta.reciprocal_condition);
        retained_bytes=base_bytes+cs_lu_bytes(lu.get());meta.factorizations=1;
    }
    FactorSolveResult solve(std::span<const double> rhs,Offset columns) {
        require(columns>0 && columns<=static_cast<Offset>(std::numeric_limits<int>::max()),"invalid number of RHS columns");
        require(!n || columns<=std::numeric_limits<std::size_t>::max()/n,"RHS extent overflow");
        require(rhs.size()==n*columns,"RHS must be row-major (unknowns, columns)");
        for(double b:rhs)require(std::isfinite(b),"nonfinite RHS");
        // Solution, column-major staging, per-column residual/CG work, metadata.
        auto temporary=bytes(rhs.size(),2*sizeof(double),options.memory_budget-retained_bytes);
        temporary+=bytes(n,12*sizeof(double),options.memory_budget-retained_bytes-temporary);
        temporary+=bytes(columns,sizeof(ColumnSolveDiagnostics),options.memory_budget-retained_bytes-temporary);
        const auto peak=retained_bytes+temporary;
        FactorSolveResult out;out.solution.resize(rhs.size());out.columns.resize(static_cast<std::size_t>(columns));
        if(meta.backend==FactorBackend::superlu) {
            std::vector<double> staging(rhs.size());
            for(std::size_t i=0;i<n;++i)for(std::size_t c=0;c<columns;++c)staging[c*n+i]=rhs[i*columns+c];
            lu_check(cs_lu_solve(lu.get(),static_cast<int>(columns),staging.data(),static_cast<std::size_t>(options.memory_budget-base_bytes-temporary)),lu.get());
            for(std::size_t i=0;i<n;++i)for(std::size_t c=0;c<columns;++c)out.solution[i*columns+c]=staging[c*n+i];
            meta.peak_bytes=std::max(meta.peak_bytes,base_bytes+temporary+cs_lu_peak(lu.get()));
        } else if(meta.backend==FactorBackend::dense_inverse) {
            for(std::size_t i=0;i<n;++i)for(std::size_t c=0;c<columns;++c) {
                double sum=0;for(std::size_t j=0;j<n;++j)sum+=inverse(static_cast<Eigen::Index>(i),static_cast<Eigen::Index>(j))*rhs[j*columns+c];
                out.solution[i*columns+c]=sum;
            }
        } else if(meta.backend==FactorBackend::dense_lu) {
            Eigen::VectorXd b(static_cast<Eigen::Index>(n));
            for(std::size_t c=0;c<columns;++c) {
                for(std::size_t i=0;i<n;++i)b[static_cast<Eigen::Index>(i)]=rhs[i*columns+c];
                const Eigen::VectorXd x=dense_lu->solve(b);
                for(std::size_t i=0;i<n;++i)out.solution[i*columns+c]=x[static_cast<Eigen::Index>(i)];
            }
        } else if(meta.backend==FactorBackend::jacobi_cg) {
            std::vector<double> b(n);CgOptions cg;cg.tolerance=1e-8;cg.max_iterations=6000;cg.memory_budget=n*8*sizeof(double);
            for(std::size_t c=0;c<columns;++c) {
                for(std::size_t i=0;i<n;++i)b[i]=rhs[i*columns+c];
                const auto answer=conjugate_gradient(*matrix,b,cg,{},diagonal);
                for(std::size_t i=0;i<n;++i)out.solution[i*columns+c]=answer.solution[i];
                out.columns[c].iterations=answer.iterations;out.columns[c].termination=answer.termination;
            }
        }
        for(double x:out.solution)finite(x);
        // Recompute both true residuals; CG recurrence residual is not a
        // substitute for checking the original and regularized equations.
        std::vector<double> x(n);
        for(std::size_t c=0;c<columns;++c) {
            for(std::size_t i=0;i<n;++i)x[i]=out.solution[i*columns+c];
            const auto ax=matrix->matvec(x);double bnorm=0,original=0,regularized=0;
            for(std::size_t i=0;i<n;++i) {
                const double b=rhs[i*columns+c],r=ax[i]-b;
                bnorm=std::hypot(bnorm,b);original=std::hypot(original,r);
                regularized=std::hypot(regularized,r+meta.regularization*x[i]);
            }
            if(bnorm==0)bnorm=1;
            out.columns[c].relative_residual=regularized/bnorm;out.columns[c].original_relative_residual=original/bnorm;
            finite(out.columns[c].relative_residual);finite(out.columns[c].original_relative_residual);
        }
        meta.peak_bytes=std::max(meta.peak_bytes,peak);++meta.solve_calls;return out;
    }
};
SpdFactor::SpdFactor(const SparseMatrix& matrix,const FactorOptions& options) try :impl_(std::make_unique<Impl>(matrix,options)) {}
catch(const std::bad_alloc&) {throw ContractError(ErrorCode::memory_limit,"solver.factor","native factor allocation failed");}
SpdFactor::~SpdFactor()=default;
const FactorMetadata& SpdFactor::metadata() const {return impl_->meta;}
FactorSolveResult SpdFactor::solve(std::span<const double> rhs,Offset columns) {
    try {return impl_->solve(rhs,columns);}
    catch(const std::bad_alloc&) {throw ContractError(ErrorCode::memory_limit,"solver.factor","native solve allocation failed");}
}
}
