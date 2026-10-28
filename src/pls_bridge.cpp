// Bridges between R and the PLS C++ implementations for BPCells matrices.
// kernelpls_cpp was extracted from BPCells r/src/matrix_utils.cpp; predict_pls_cpp
// computes Y_hat = X_new %*% B + B0 without materialising X_new into memory.

#include <Rcpp.h>
#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]

#include "R_interrupts.h"
#include "R_xptr_wrapper.h"
#include "pls.h"

using namespace Rcpp;
using namespace BPCells;

// [[Rcpp::export]]
List cppls_cpp(SEXP matrix_X, Eigen::Map<Eigen::MatrixXd> Y,
               Eigen::Map<Eigen::MatrixXd> Y_add,
               int ncomp, bool center, bool stripped, bool transpose,
               double w_tol, double X_tol) {
    auto X = take_unique_xptr<MatrixLoader<double>>(matrix_X);

    // Copy Y and Y_add so the background thread owns the data
    Eigen::MatrixXd Y_copy = Y;
    Eigen::MatrixXd Y_add_copy = Y_add;

    CPPLSResult res = run_with_R_interrupt_check(
        cppls, X.get(), Y_copy, Y_add_copy, ncomp, center, stripped, transpose,
        w_tol, X_tol
    );

    List result = List::create(
        Named("coefficients") = res.coefficients,
        Named("Xmeans") = res.Xmeans,
        Named("Ymeans") = res.Ymeans,
        Named("gammas") = res.gammas,
        Named("canonical.correlations") = res.canonical_correlations,
        Named("A") = res.A_weights
    );

    // Convert smallNorm from 0-indexed to 1-indexed for R
    if (!res.smallNorm.empty()) {
        Rcpp::IntegerVector sn(res.smallNorm.begin(), res.smallNorm.end());
        for (int i = 0; i < sn.size(); i++) sn[i] += 1;
        result["smallNorm"] = sn;
    } else {
        result["smallNorm"] = Rcpp::IntegerVector(0);
    }

    if (!stripped) {
        result["scores"] = res.scores;
        result["loadings"] = res.loadings;
        result["loading.weights"] = res.loading_weights;
        result["Yscores"] = res.Yscores;
        result["Yloadings"] = res.Yloadings;
        result["projection"] = res.projection;
        result["fitted.values"] = res.fitted_values;
        result["residuals"] = res.residuals;
        result["Xvar"] = res.Xvar;
        result["Xtotvar"] = res.Xtotvar;
    }

    return result;
}

// [[Rcpp::export]]
List kernelpls_cpp(SEXP matrix_X, Eigen::Map<Eigen::MatrixXd> Y, int ncomp, bool center, bool stripped, bool transpose) {
    auto X = take_unique_xptr<MatrixLoader<double>>(matrix_X);

    // Copy Y into an Eigen::MatrixXd so the background thread owns the data
    // (the Map points to R-managed memory which must not be accessed from a
    // background thread without a copy).
    Eigen::MatrixXd Y_copy = Y;

    KernelPLSResult res = run_with_R_interrupt_check(
        kernelpls, X.get(), Y_copy, ncomp, center, stripped, transpose
    );

    List result = List::create(
        Named("coefficients") = res.coefficients,
        Named("Xmeans") = res.Xmeans,
        Named("Ymeans") = res.Ymeans
    );

    if (!stripped) {
        result["scores"] = res.scores;
        result["loadings"] = res.loadings;
        result["loading.weights"] = res.loading_weights;
        result["Yscores"] = res.Yscores;
        result["Yloadings"] = res.Yloadings;
        result["projection"] = res.projection;
        result["fitted.values"] = res.fitted_values;
        result["residuals"] = res.residuals;
        result["Xvar"] = res.Xvar;
        result["Xtotvar"] = res.Xtotvar;
    }

    return result;
}

// Core implementation run on background thread.
// X is the on-disk matrix (p x n when transpose=true, n x p when transpose=false).
// B is p x q coefficient matrix, B0 is q-length intercept.
// Returns n x q prediction matrix.
Eigen::MatrixXd predict_pls_impl(
    MatrixLoader<double> *X,
    const Eigen::MatrixXd &B,
    const Eigen::VectorXd &B0,
    bool transpose,
    std::atomic<bool> *user_interrupt
) {
    Eigen::MatrixXd XB;
    if (transpose) {
        // Stored p x n: X_stored^T * B^T gives (n x p) * (p x q) ... but we use
        // denseMultiplyLeft(Bt) which computes Bt * X_stored = (q x p)(p x n) = q x n
        // then transpose to get n x q.
        Eigen::MatrixXd Bt = B.transpose();  // q x p
        Eigen::Map<Eigen::MatrixXd> Bt_map(Bt.data(), Bt.rows(), Bt.cols());
        XB = X->denseMultiplyLeft(Bt_map, user_interrupt).transpose();  // n x q
    } else {
        // Stored n x p: X_stored * B gives n x q directly.
        Eigen::Map<Eigen::MatrixXd> B_map(
            const_cast<double*>(B.data()), B.rows(), B.cols());
        XB = X->denseMultiplyRight(B_map, user_interrupt);  // n x q
    }

    // Add intercept: each row gets + B0^T
    XB.rowwise() += B0.transpose();
    return XB;
}

// [[Rcpp::export]]
Eigen::MatrixXd predict_pls_cpp(SEXP matrix_X,
                                 Eigen::Map<Eigen::MatrixXd> B,
                                 Eigen::Map<Eigen::VectorXd> B0,
                                 bool transpose) {
    auto X = take_unique_xptr<MatrixLoader<double>>(matrix_X);

    // Copy R-managed Eigen data so the background thread owns the memory.
    Eigen::MatrixXd B_copy = B;
    Eigen::VectorXd B0_copy = B0;

    return run_with_R_interrupt_check(
        predict_pls_impl, X.get(), B_copy, B0_copy, transpose
    );
}
