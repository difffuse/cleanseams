// SPDX-License-Identifier: GPL-3.0-or-later
#include "cleanseams/chart_solver.hpp"
#include "cleanseams/errors.hpp"
#include <Eigen/SVD>
#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
namespace cleanseams {
namespace {
void storage(bool yes) {if(!yes)throw ContractError(ErrorCode::memory_limit,"solver.least_squares","least-squares workspace exceeds memory budget");}
void finite(double value) {if(!std::isfinite(value))throw ContractError(ErrorCode::numerical_failure,"solver.least_squares","nonfinite least-squares arithmetic");}
void reserve_bound(std::uint64_t& total,std::uint64_t count,std::uint64_t stride,std::uint64_t budget) {
    storage(total<=budget && count<=(budget-total)/stride);total+=count*stride;
}
SparseMatrix normal_matrix(const SparseMatrix& matrix,std::uint64_t budget,std::uint64_t base,std::uint64_t& peak) {
    const auto csr=matrix.csr();const auto n=matrix.shape()[1];std::uint64_t products=0;
    for(std::size_t row=0;row<matrix.shape()[0];++row) {
        const auto count=csr.row_offsets[row+1]-csr.row_offsets[row];
        storage(!count || count<=((budget-base)/128-products)/count);products+=count*count;
    }
    peak=base+products*128;
    // The source forms every in-row product csr[a]*csr[b] as a COO entry (column of a,
    // column of b) in (row, a, b) order and compresses that matrix, summing the duplicates
    // of one entry in COO order, i.e. over the rows of M in ascending order, starting from
    // the first product. GB performance: the same sums are formed directly, row by row of
    // the normal matrix (Gustavson): the rows of M holding column i are visited in
    // ascending order, and each product is the same operand pair. The products and the
    // sums are therefore bit-for-bit those of the COO compression, without the product
    // matrix and its sort. Error order is kept: a nonfinite product is reported before a
    // nonfinite sum. (The COO matrix's own budget checks are implied by the products check
    // above: 128 bytes per product plus 64 per column are within the allowance.)
    const auto m=static_cast<std::size_t>(matrix.shape()[0]);
    std::vector<std::size_t> column_start(static_cast<std::size_t>(n)+1,0);
    for(auto c:csr.columns)++column_start[static_cast<std::size_t>(c)+1];
    std::partial_sum(column_start.begin(),column_start.end(),column_start.begin());
    std::vector<std::size_t> fill(column_start.begin(),column_start.end()-1),row_of(csr.columns.size()),entry_of(csr.columns.size());
    for(std::size_t row=0;row<m;++row)for(auto a=csr.row_offsets[row];a<csr.row_offsets[row+1];++a) {
        const auto at=fill[static_cast<std::size_t>(csr.columns[a])]++;row_of[at]=row;entry_of[at]=static_cast<std::size_t>(a);
    }
    std::vector<TransferIndex> rows,cols;std::vector<double> values;std::vector<double> sum(static_cast<std::size_t>(n),0.0);
    std::vector<std::size_t> seen(static_cast<std::size_t>(n),static_cast<std::size_t>(-1));std::vector<TransferIndex> touched;
    bool nonfinite_product=false;
    for(std::size_t i=0;i<n;++i) {
        touched.clear();
        for(auto at=column_start[i];at<column_start[i+1];++at) {
            const auto row=row_of[at];const double left=csr.values[entry_of[at]];
            for(auto b=csr.row_offsets[row];b<csr.row_offsets[row+1];++b) {
                const auto j=static_cast<std::size_t>(csr.columns[b]);const double product=left*csr.values[b];
                nonfinite_product=nonfinite_product || !std::isfinite(product);
                if(seen[j]!=i) {seen[j]=i;sum[j]=product;touched.push_back(static_cast<TransferIndex>(j));}
                else sum[j]+=product;
            }
        }
        std::sort(touched.begin(),touched.end());
        for(auto j:touched) {rows.push_back(static_cast<TransferIndex>(i));cols.push_back(j);values.push_back(sum[static_cast<std::size_t>(j)]);}
    }
    if(nonfinite_product)finite(std::numeric_limits<double>::infinity());
    for(double v:values)if(!std::isfinite(v))throw ContractError(ErrorCode::numerical_failure,"solver.sparse","nonfinite arithmetic result");
    return SparseMatrix(rows,cols,values,{n,n},budget);
}
std::unique_ptr<SparseMatrix> regularized_trace(const SparseMatrix& normal,std::uint64_t budget) {
    std::vector<TransferIndex> rows(normal.rows().begin(),normal.rows().end()),cols(normal.columns().begin(),normal.columns().end());
    std::vector<double> values(normal.values().begin(),normal.values().end());
    const auto count=normal.values().size()+normal.shape()[0];rows.reserve(count);cols.reserve(count);values.reserve(count);
    for(std::size_t i=0;i<normal.shape()[0];++i) {rows.push_back(static_cast<TransferIndex>(i));cols.push_back(static_cast<TransferIndex>(i));values.push_back(1e-12);}
    return std::make_unique<SparseMatrix>(rows,cols,values,normal.shape(),budget);
}
void svd_fallback(LinearLeastSquaresResult& result,const SparseMatrix& matrix,std::span<const double> rhs,
                  std::uint64_t budget,std::uint64_t retained) {
    const auto m=matrix.shape()[0],n=matrix.shape()[1];
    storage(!m || n<=(budget-retained)/(12*sizeof(double))/m);
    result.numerics.peak_bytes=std::max(result.numerics.peak_bytes,retained+m*n*12*sizeof(double));
    Eigen::MatrixXd dense=Eigen::MatrixXd::Zero(static_cast<Eigen::Index>(m),static_cast<Eigen::Index>(n));
    for(std::size_t k=0;k<matrix.values().size();++k)dense(matrix.rows()[k],matrix.columns()[k])+=matrix.values()[k];
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(dense,Eigen::ComputeThinU|Eigen::ComputeThinV);
    svd.setThreshold(std::numeric_limits<double>::epsilon()*static_cast<double>(std::max(m,n)));
    const Eigen::Map<const Eigen::VectorXd> b(rhs.data(),static_cast<Eigen::Index>(rhs.size()));
    const Eigen::VectorXd x=svd.solve(b);result.solution.assign(x.data(),x.data()+x.size());
    result.numerics.backend=FactorBackend::dense_least_squares;result.numerics.version="Eigen-3.4.0/JacobiSVD/numpy-rcond";
    result.numerics.regularization=0;result.numerics.factorizations=1;result.numerics.solve_calls=1;
    result.numerics.fallback_reason="singular dense normal system after 1e-12 regularization";
    if(svd.singularValues().size()) {
        const auto& s=svd.singularValues();result.numerics.reciprocal_condition=s[0]>0?s[s.size()-1]/s[0]:0;
    }
}
LinearLeastSquaresResult solve_impl(const SparseMatrix& matrix,std::span<const double> rhs,NumericalPolicy policy,std::uint64_t budget,bool retain) {
    if(rhs.size()!=matrix.shape()[0])throw ContractError(ErrorCode::invalid_input,"solver.least_squares.rhs","RHS row count differs from matrix");
    for(double value:rhs)if(!std::isfinite(value))throw ContractError(ErrorCode::invalid_input,"solver.least_squares.rhs","nonfinite RHS");
    if(policy!=NumericalPolicy::scipy_compatible && policy!=NumericalPolicy::numpy_compatible)
        throw ContractError(ErrorCode::invalid_input,"solver.least_squares.policy","unknown numerical policy");
    const auto m=matrix.shape()[0],n=matrix.shape()[1];std::uint64_t base=1024;
    reserve_bound(base,m+n,64,budget);reserve_bound(base,matrix.values().size(),64,budget);
    LinearLeastSquaresResult result;result.numerics.policy=policy;result.numerics.system=FactorSystem::least_squares_normal;
    result.numerics.peak_bytes=base;
    if(n==0) {result.numerics.version="empty-v1";result.numerics.solve_calls=1;}
    else if(policy==NumericalPolicy::numpy_compatible && n>1800) {
        reserve_bound(base,n,12*sizeof(double),budget);result.numerics.peak_bytes=base;
        const auto normal_rhs=matrix.rmatvec(rhs),diagonal=matrix.diag_of_normal();
        CgOptions options;options.normal_equations=true;options.tolerance=1e-9;options.max_iterations=20000;
        options.memory_budget=budget-(base-n*12*sizeof(double));
        const auto cg=conjugate_gradient(matrix,normal_rhs,options,{},diagonal);result.solution=cg.solution;
        result.diagnostics.iterations=cg.iterations;result.diagnostics.relative_residual=cg.relative_residual;
        result.diagnostics.original_relative_residual=cg.relative_residual;result.diagnostics.termination=cg.termination;
        result.numerics.backend=FactorBackend::jacobi_cg;result.numerics.version="CleanSeams-CG-v1/normal/tol=1e-9/maxiter=20000";
        result.numerics.solve_calls=1;
        if(retain)result.normal_rhs=normal_rhs;
    } else {
        std::uint64_t product_peak=0;const auto normal=normal_matrix(matrix,budget,base,product_peak);
        const auto normal_rhs=matrix.rmatvec(rhs);auto retained=base;
        reserve_bound(retained,normal.values().size(),24,budget);
        if(retain)reserve_bound(retained,normal.values().size()+n,72,budget);
        FactorOptions options;options.policy=policy;options.system=FactorSystem::least_squares_normal;options.memory_budget=budget-retained;
        try {
            SpdFactor factor(normal,options);const auto solved=factor.solve(normal_rhs,1);
            result.solution=solved.solution;result.diagnostics=solved.columns[0];result.numerics=factor.metadata();
            result.numerics.peak_bytes+=retained;
        } catch(const ContractError& error) {
            if(policy!=NumericalPolicy::numpy_compatible || error.code!=ErrorCode::numerical_failure || error.field!="solver.dense_solve.singular")throw;
            svd_fallback(result,matrix,rhs,budget,retained);
            const auto ax=normal.matvec(result.solution);double norm=0,residual=0;
            for(std::size_t i=0;i<n;++i) {norm=std::hypot(norm,normal_rhs[i]);residual=std::hypot(residual,ax[i]-normal_rhs[i]);}
            result.diagnostics.relative_residual=result.diagnostics.original_relative_residual=residual/(norm!=0?norm:1);
        }
        result.numerics.peak_bytes=std::max(result.numerics.peak_bytes,product_peak);
        if(retain) {result.normal=regularized_trace(normal,budget);result.normal_rhs=normal_rhs;}
    }
    for(double x:result.solution)finite(x);
    const auto ax=matrix.matvec(result.solution);double norm=0,residual=0;
    for(std::size_t i=0;i<m;++i) {norm=std::hypot(norm,rhs[i]);residual=std::hypot(residual,ax[i]-rhs[i]);}
    result.equation_relative_residual=residual/(norm!=0?norm:1);finite(result.equation_relative_residual);
    finite(result.diagnostics.relative_residual);return result;
}
}
LinearLeastSquaresResult solve_least_squares(const SparseMatrix& matrix,std::span<const double> rhs,NumericalPolicy policy,std::uint64_t budget,bool retain) {
    try {return solve_impl(matrix,rhs,policy,budget,retain);}
    catch(const std::bad_alloc&) {throw ContractError(ErrorCode::memory_limit,"solver.least_squares","native least-squares allocation failed");}
}
}
