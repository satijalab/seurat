#include <RcppEigen.h>
#include <algorithm>
#include <cmath>

// [[Rcpp::depends(RcppEigen)]]

using namespace Rcpp;

// Exact PCA via the Gram matrix and a symmetric eigendecomposition, using
// Eigen's own (BLAS-independent) kernels. `object` is the feature-by-cell
// scaled matrix X. Since X = U D V', the eigenvectors of X X' are the feature
// loadings U and the eigenvalues are D^2, so the top npcs components are
// obtained without an iterative SVD or a dense transpose. Cell embeddings are
// X' U = V D (weighted by variance) or V (unweighted). This matches the output
// of RunPCA's approx path up to sign, and is exact for the leading components
// (it loses precision only for singular values near the numerical noise floor,
// which the top npcs of scaled data are not).
//
// [[Rcpp::export(rng = false)]]
List EigenGramPCA(const Eigen::Map<Eigen::MatrixXd> object,
                  int npcs,
                  bool weight_by_var) {
  const int nfeatures = object.rows();
  const int ncells = object.cols();
  npcs = std::min(npcs, nfeatures);

  // Gram matrix X X' (features x features); only the lower triangle is formed.
  Eigen::MatrixXd XtX = Eigen::MatrixXd::Zero(nfeatures, nfeatures);
  XtX.selfadjointView<Eigen::Lower>().rankUpdate(object);

  // Symmetric eigensolver returns eigenvalues in ascending order.
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(XtX);
  const Eigen::VectorXd& evalues = es.eigenvalues();
  const Eigen::MatrixXd& evectors = es.eigenvectors();

  // Take the top npcs (largest eigenvalues) in descending order.
  Eigen::MatrixXd loadings(nfeatures, npcs);
  Eigen::VectorXd d(npcs);
  for (int j = 0; j < npcs; ++j) {
    const int src = nfeatures - 1 - j;
    loadings.col(j) = evectors.col(src);
    d(j) = std::sqrt(std::max(0.0, evalues(src)));
  }

  Eigen::MatrixXd embeddings = object.transpose() * loadings;  // X' U = V D
  if (!weight_by_var) {
    for (int j = 0; j < npcs; ++j) {
      if (d(j) > 0) {
        embeddings.col(j) /= d(j);
      }
    }
  }

  Eigen::VectorXd sdev = d / std::sqrt(static_cast<double>(std::max(1, ncells - 1)));

  return List::create(
    _["loadings"] = loadings,
    _["embeddings"] = embeddings,
    _["sdev"] = sdev
  );
}
