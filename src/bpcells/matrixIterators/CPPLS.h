// Copyright 2024 BPCells contributors
//
// Licensed under the Apache License, Version 2.0
// https://www.apache.org/licenses/LICENSE-2.0 or the MIT license
// https://opensource.org/licenses/MIT, at your option. This file
// may not be copied, modified, or distributed except according to
// those terms.

#pragma once

#include <atomic>
#include <vector>

#include <Eigen/Core>

#include "MatrixIterator.h"

namespace BPCells {

struct CCAResult {
    Eigen::MatrixXd A;   // weight vectors [dz × dxy]
    double r;            // first canonical correlation
};

struct CPPLSResult {
    Eigen::MatrixXd coefficients;    // p × (q_prim*ncomp), cumulative coefficients
    Eigen::MatrixXd scores;          // n × ncomp (T)
    Eigen::MatrixXd loadings;        // p × ncomp (P)
    Eigen::MatrixXd loading_weights; // p × ncomp (W)
    Eigen::MatrixXd Yscores;         // n × ncomp (U)
    Eigen::MatrixXd Yloadings;       // q_prim × ncomp (Q, primary only)
    Eigen::MatrixXd projection;      // p × ncomp (R)
    Eigen::VectorXd Xmeans;          // p × 1
    Eigen::VectorXd Ymeans;          // q_prim × 1
    Eigen::MatrixXd fitted_values;   // n × (q_prim*ncomp), cumulative fitted values
    Eigen::MatrixXd residuals;       // n × (q_prim*ncomp), cumulative residuals
    Eigen::VectorXd Xvar;            // ncomp × 1
    double Xtotvar;                  // scalar
    // CPLS-specific fields:
    Eigen::VectorXd gammas;                  // power values (0.5 for all components)
    Eigen::VectorXd canonical_correlations;  // squared CC per component
    Eigen::MatrixXd A_weights;               // CCA weight vectors [q_full × ncomp]
    std::vector<int> smallNorm;              // indices of near-zero-norm variables
};

// Canonical correlation analysis helper (in-memory).
// Z: [n × dz], Yprim: [n × dy]
// Returns: CCA weight matrix A [dz × dxy] and first canonical correlation r.
// Z is taken by value so the caller can std::move() it in, avoiding a copy.
CCAResult cancorr(Eigen::MatrixXd Z, const Eigen::MatrixXd &Yprim);

// CPLS (Canonical PLS) algorithm for on-disk BPCells matrices.
//
// Parameters:
//   X            - predictor matrix accessed via MatrixLoader (stays on disk).
//                  When transpose=false, X is n x p (samples x features).
//                  When transpose=true,  X is p x n (features x samples, the
//                  BPCells convention) and operations are swapped internally.
//   Y            - primary response matrix (n x q_prim), must fit in memory
//   Y_add        - additional response matrix (n x q_add), can be 0-column
//   ncomp        - number of PLS components to compute
//   center       - whether to center X and Y (subtract column means)
//   stripped     - if true, only compute coefficients/means (skip scores, fitted, etc.)
//   transpose    - if true, treat X as transposed (features x samples layout)
//   w_tol        - threshold for zeroing small loading weights
//   X_tol        - threshold for detecting small-norm variables
//   user_interrupt - pointer to atomic bool for cooperative cancellation
CPPLSResult cppls(
    MatrixLoader<double> *X,
    const Eigen::MatrixXd &Y,
    const Eigen::MatrixXd &Y_add,
    int ncomp,
    bool center,
    bool stripped,
    bool transpose,
    double w_tol,
    double X_tol,
    std::atomic<bool> *user_interrupt);

} // namespace BPCells
