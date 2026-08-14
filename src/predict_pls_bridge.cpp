// Bridge between R and BPCells for PLS prediction on IterableMatrix objects.
// Computes Y_hat = X_new %*% B + B0 without materialising X_new into memory.

#include <Rcpp.h>
#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]

#include "R_interrupts.h"
#include "R_xptr_wrapper.h"
#include "bpcells/matrixIterators/MatrixIterator.h"

using namespace Rcpp;
using namespace BPCells;

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
