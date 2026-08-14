// Copyright 2024 BPCells contributors
//
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

#include "KernelPLS.h"

#include <Eigen/Eigenvalues>

namespace BPCells {

KernelPLSResult kernelpls(
    MatrixLoader<double> *X,
    const Eigen::MatrixXd &Y,
    int ncomp,
    bool center,
    bool stripped,
    bool transpose,
    std::atomic<bool> *user_interrupt
) {
    KernelPLSResult result;

    // When transpose=false, X is n x p (samples x features) — standard PLS layout.
    // When transpose=true,  X is p x n (features x samples) — BPCells layout.
    // We define nobj (n) and npred (p) according to the conceptual n x p orientation.
    const int nobj  = transpose ? X->cols() : X->rows();
    const int npred = transpose ? X->rows() : X->cols();
    const int nresp = Y.cols();

    // Helper lambdas to abstract away the transpose logic.
    // "conceptual X" is always n x p regardless of stored layout.
    //
    // conceptualXtY: compute X^T * Y  (p x q)
    //   non-transposed: denseMultiplyLeft(Y^T)^T  = (Y^T * X)^T = X^T * Y
    //   transposed:     denseMultiplyRight(Y)      = X_stored * Y  (p x n * n x q = p x q)
    auto conceptualXtY = [&](const Eigen::MatrixXd &M) -> Eigen::MatrixXd {
        if (transpose) {
            Eigen::Map<Eigen::MatrixXd> M_map(
                const_cast<double*>(M.data()), M.rows(), M.cols());
            return X->denseMultiplyRight(M_map, user_interrupt);
        } else {
            Eigen::MatrixXd Mt = M.transpose();
            Eigen::Map<Eigen::MatrixXd> Mt_map(Mt.data(), Mt.rows(), Mt.cols());
            return X->denseMultiplyLeft(Mt_map, user_interrupt).transpose();
        }
    };

    // conceptualXr: compute X * r  (n x 1)  where r is p x 1
    //   non-transposed: vecMultiplyRight(r)   = X * r
    //   transposed:     vecMultiplyLeft(r)    = X_stored^T * r  (n x p * p x 1 = n x 1)
    auto conceptualXr = [&](Eigen::VectorXd &r) -> Eigen::VectorXd {
        Eigen::Map<Eigen::VectorXd> r_map(r.data(), r.size());
        if (transpose) {
            return X->vecMultiplyLeft(r_map, user_interrupt);
        } else {
            return X->vecMultiplyRight(r_map, user_interrupt);
        }
    };

    // conceptualXtv: compute X^T * v  (p x 1)  where v is n x 1
    //   non-transposed: vecMultiplyLeft(v)    = X^T * v
    //   transposed:     vecMultiplyRight(v)   = X_stored * v  (p x n * n x 1 = p x 1)
    auto conceptualXtv = [&](Eigen::VectorXd &v) -> Eigen::VectorXd {
        Eigen::Map<Eigen::VectorXd> v_map(v.data(), v.size());
        if (transpose) {
            return X->vecMultiplyRight(v_map, user_interrupt);
        } else {
            return X->vecMultiplyLeft(v_map, user_interrupt);
        }
    };

    // 1. Center X and Y
    Eigen::MatrixXd Y_centered = Y;
    if (center) {
        // Compute X column means (conceptual columns = features, p means).
        // When transpose=false: features are stored columns → col stats
        // When transpose=true:  features are stored rows   → row stats
        StatsResult X_stats = transpose
            ? X->computeMatrixStats(Stats::Mean, Stats::None, user_interrupt)
            : X->computeMatrixStats(Stats::None, Stats::Mean, user_interrupt);
        // Extract the mean vector (p × 1)
        if (transpose) {
            result.Xmeans = X_stats.rowMean().row(0).matrix().transpose();
        } else {
            result.Xmeans = X_stats.colMean().row(0).matrix().transpose();
        }

        // Compute and subtract Y column means
        result.Ymeans = Y.colwise().mean();
        Y_centered.rowwise() -= result.Ymeans.transpose();
    } else {
        result.Xmeans = Eigen::VectorXd::Zero(npred);
        result.Ymeans = Eigen::VectorXd::Zero(nresp);
    }

    // 2. Compute X^T * Y_centered (p x q).
    // When center=true, we want (X - Xmeans)^T * Y_centered.  Because Y is
    // already centered (column sums are zero), the centering correction vanishes:
    //   (X - mu)^T * Y_c = X^T * Y_c - mu * 1^T * Y_c = X^T * Y_c
    Eigen::MatrixXd XtY = conceptualXtY(Y_centered);

    // 3. Initialize component matrices
    Eigen::MatrixXd R_mat(npred, ncomp);  // Projection matrix
    Eigen::MatrixXd P(npred, ncomp);      // X-loadings
    Eigen::MatrixXd tQ(ncomp, nresp);     // Y-loadings (transposed)

    Eigen::VectorXd tsqs(ncomp); // store t'*t for each component

    if (!stripped) {
        result.scores = Eigen::MatrixXd(nobj, ncomp);
        result.Yscores = Eigen::MatrixXd(nobj, ncomp);
        result.loading_weights = Eigen::MatrixXd(npred, ncomp);
        result.Xvar = Eigen::VectorXd(ncomp);
    }

    // 4. Iterate over components
    for (int a = 0; a < ncomp; a++) {
        if (user_interrupt != nullptr && *user_interrupt) break;

        // 4.1 Compute loading weight w_a
        Eigen::VectorXd w_a;
        if (nresp == 1) {
            w_a = XtY.col(0);
            double norm = w_a.norm();
            if (norm > 0) w_a /= norm;
        } else if (nresp < npred) {
            Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eig(XtY.transpose() * XtY);
            Eigen::VectorXd q_vec = eig.eigenvectors().col(nresp - 1);
            w_a = XtY * q_vec;
            double norm = w_a.norm();
            if (norm > 0) w_a /= norm;
        } else {
            Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eig(XtY * XtY.transpose());
            w_a = eig.eigenvectors().col(npred - 1);
        }

        // 4.2 Compute projection r_a (orthogonalize w_a against previous P loadings)
        // r_a = w_a - sum_j (P_j' * w_a) * R_j  =  w_a - R * (P' * w_a)
        Eigen::VectorXd r_a = w_a;
        if (a > 5) {
            r_a -= R_mat.leftCols(a) * (P.leftCols(a).transpose() * w_a);
        } else if (a > 0) {
            for (int j = 0; j < a; j++) {
                r_a -= (P.col(j).dot(w_a)) * R_mat.col(j);
            }
        }

        // 4.3 Compute scores: t_a = (X - Xmeans) * r_a  (n × 1)
        Eigen::VectorXd t_a = conceptualXr(r_a);
        if (center) {
            t_a.array() -= result.Xmeans.dot(r_a);
        }

        double tsq = t_a.squaredNorm();

        // 4.4 Compute X-loadings: p_a = (X - Xmeans)^T * t_a / tsq  (p × 1)
        Eigen::VectorXd p_a = conceptualXtv(t_a);
        p_a /= tsq;
        if (center) {
            p_a -= result.Xmeans * (t_a.sum() / tsq);
        }

        // 4.5 Compute Y-loading: q_a = XtY^T * r_a / tsq
        Eigen::VectorXd q_a = XtY.transpose() * r_a / tsq;

        // 4.6 Deflate XtY
        XtY -= (tsq * p_a) * q_a.transpose();

        // 4.7 Store components
        R_mat.col(a) = r_a;
        P.col(a) = p_a;
        tQ.row(a) = q_a.transpose();
        tsqs(a) = tsq;

        if (!stripped) {
            result.scores.col(a) = t_a;
            result.loading_weights.col(a) = w_a;

            Eigen::VectorXd u_a = Y_centered * q_a / q_a.squaredNorm();
            if (a > 0) {
                Eigen::VectorXd proj = result.scores.leftCols(a).transpose() * u_a;
                for (int j = 0; j < a; j++) {
                    proj(j) /= tsqs(j);
                }
                u_a -= result.scores.leftCols(a) * proj;
            }
            result.Yscores.col(a) = u_a;

            result.Xvar(a) = p_a.squaredNorm() * tsq;
        }
    }

    // 5. Compute cumulative coefficients B[,,a] = R[,1:a] * tQ[1:a,]
    result.coefficients.resize(npred, nresp * ncomp);
    for (int a = 0; a < ncomp; a++) {
        result.coefficients.block(0, a * nresp, npred, nresp) =
            R_mat.leftCols(a + 1) * tQ.topRows(a + 1);
    }

    // 6. Compute fitted values, residuals, and remaining outputs
    if (!stripped) {
        result.loadings = P;
        result.Yloadings = tQ.transpose();
        result.projection = R_mat;

        result.fitted_values.resize(nobj, nresp * ncomp);
        result.residuals.resize(nobj, nresp * ncomp);
        for (int a = 0; a < ncomp; a++) {
            Eigen::MatrixXd fitted_a =
                result.scores.leftCols(a + 1) * tQ.topRows(a + 1);
            result.residuals.block(0, a * nresp, nobj, nresp) =
                Y_centered - fitted_a;
            if (center) {
                fitted_a.rowwise() += result.Ymeans.transpose();
            }
            result.fitted_values.block(0, a * nresp, nobj, nresp) = fitted_a;
        }

        // Total X variance: use conceptual column (feature) variances
        // When transpose=false: feature variances are column variances
        // When transpose=true:  feature variances are row variances
        StatsResult var_stats = transpose
            ? X->computeMatrixStats(Stats::Variance, Stats::None, user_interrupt)
            : X->computeMatrixStats(Stats::None, Stats::Variance, user_interrupt);
        if (transpose) {
            result.Xtotvar = var_stats.rowVariance().sum() * (nobj - 1);
        } else {
            result.Xtotvar = var_stats.colVariance().sum() * (nobj - 1);
        }
    }

    return result;
}

} // namespace BPCells
