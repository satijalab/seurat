// Bridge between R and BPCells CPPLS C++ implementation.

#include <Rcpp.h>
#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]

#include "R_interrupts.h"
#include "R_xptr_wrapper.h"
#include "bpcells/matrixIterators/CPPLS.h"

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
