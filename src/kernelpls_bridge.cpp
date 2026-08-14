// Bridge between R and BPCells kernelpls C++ implementation.
// Extracted from BPCells r/src/matrix_utils.cpp (kernelpls_cpp function).

#include <Rcpp.h>
#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]

#include "R_interrupts.h"
#include "R_xptr_wrapper.h"
#include "bpcells/matrixIterators/KernelPLS.h"

using namespace Rcpp;
using namespace BPCells;

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
