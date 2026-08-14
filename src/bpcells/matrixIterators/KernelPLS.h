// Copyright 2024 BPCells contributors
//
// Licensed under the Apache License, Version 2.0
// https://www.apache.org/licenses/LICENSE-2.0 or the MIT license
// https://opensource.org/licenses/MIT, at your option. This file
// may not be copied, modified, or distributed except according to
// those terms.

#pragma once

#include <atomic>

#include <Eigen/Core>

#include "MatrixIterator.h"

namespace BPCells {

struct KernelPLSResult {
    Eigen::MatrixXd coefficients;    // p × (q*ncomp), cumulative coefficients
    Eigen::MatrixXd scores;          // n × ncomp (T)
    Eigen::MatrixXd loadings;        // p × ncomp (P)
    Eigen::MatrixXd loading_weights; // p × ncomp (W)
    Eigen::MatrixXd Yscores;         // n × ncomp (U)
    Eigen::MatrixXd Yloadings;       // q × ncomp
    Eigen::MatrixXd projection;      // p × ncomp (R)
    Eigen::VectorXd Xmeans;          // p × 1
    Eigen::VectorXd Ymeans;          // q × 1
    Eigen::MatrixXd fitted_values;   // n × (q*ncomp), cumulative fitted values
    Eigen::MatrixXd residuals;       // n × (q*ncomp), cumulative residuals
    Eigen::VectorXd Xvar;            // ncomp × 1
    double Xtotvar;                  // scalar
};

// Kernel PLS algorithm (Dayal & MacGregor, 1997)
//
// Parameters:
//   X            - predictor matrix accessed via MatrixLoader (stays on disk).
//                  When transpose=false, X is n x p (samples x features).
//                  When transpose=true,  X is p x n (features x samples, the
//                  BPCells convention) and operations are swapped internally.
//   Y            - response matrix (n x q), must fit in memory
//   ncomp        - number of PLS components to compute
//   center       - whether to center X and Y (subtract column means)
//   stripped     - if true, only compute coefficients/means (skip scores, fitted, etc.)
//   transpose    - if true, treat X as transposed (features x samples layout)
//   user_interrupt - pointer to atomic bool for cooperative cancellation
KernelPLSResult kernelpls(
    MatrixLoader<double> *X,
    const Eigen::MatrixXd &Y,
    int ncomp,
    bool center,
    bool stripped,
    bool transpose,
    std::atomic<bool> *user_interrupt);

} // namespace BPCells
