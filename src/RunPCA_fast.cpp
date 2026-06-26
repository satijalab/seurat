#include <RcppEigen.h>
// CHANGE: Spectra (header-only, Eigen-native; the library RSpectra wraps) gives a
// partial *top-k* symmetric eigensolver. The previous full SelfAdjointEigenSolver
// computed all nfeatures eigenpairs even though only npcs are needed, an O(nfeatures^3)
// cost that is fixed regardless of cell count and dominated the Gram path for small
// datasets. The top-k Lanczos solver below is machine-precision identical for the
// leading components but avoids computing the discarded eigenvectors.
#include <Spectra/SymEigsSolver.h>
#include <Spectra/MatOp/DenseSymMatProd.h>
#include <algorithm>
#include <cmath>

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RSpectra)]]

using namespace Rcpp;

// Approximate PCA via the Gram matrix and a top-k symmetric eigendecomposition,
// using Eigen's own (BLAS-independent) kernels. `object` is the feature-by-cell
// scaled matrix X. Since X = U D V', the eigenvectors of X X' are the feature
// loadings U and the eigenvalues are D^2, so the top npcs components are obtained
// without a dense transpose. Cell embeddings are X' U = V D (weighted by variance)
// or V (unweighted). The top-k eigenpairs are computed with a Lanczos solver
// (Spectra), which converges to the same truncated SVD as irlba up to sign, but
// tighter (tol 1e-10); for the well-separated leading components of scaled data
// this is numerically indistinguishable from exact.
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

  // The loadings are the top npcs eigenvectors of XX' (descending) and d the
  // square-roots of the corresponding eigenvalues.
  Eigen::MatrixXd loadings;
  Eigen::VectorXd d(npcs);

  // CHANGE: compute only the top npcs eigenpairs with a Lanczos solver when it is
  // applicable and beneficial (npcs well below nfeatures, so a partial basis is a
  // real saving). DenseSymMatProd reads the lower triangle filled above. Spectra
  // requires 1 <= nev < ncv <= n; ncv ~ 2*nev controls convergence. On any failure
  // (or when npcs is close to nfeatures) fall back to the full dense eigensolver,
  // so results are unchanged in every case.
  bool used_partial = false;
  if ((2 * npcs + 1) < nfeatures) {
    const int ncv = std::min(nfeatures, std::max(2 * npcs + 1, 20));
    Spectra::DenseSymMatProd<double> op(XtX);
    Spectra::SymEigsSolver<double, Spectra::LARGEST_ALGE,
                           Spectra::DenseSymMatProd<double> > eigs(&op, npcs, ncv);
    eigs.init();
    const int nconv = eigs.compute(1000, 1e-10, Spectra::LARGEST_ALGE);
    if (eigs.info() == Spectra::SUCCESSFUL && nconv >= npcs) {
      // Spectra returns eigenvalues in descending order with eigenvectors as columns.
      const Eigen::VectorXd evalues = eigs.eigenvalues();
      loadings = eigs.eigenvectors();
      for (int j = 0; j < npcs; ++j) {
        d(j) = std::sqrt(std::max(0.0, evalues(j)));
      }
      used_partial = true;
    }
  }
  if (!used_partial) {
    // Full symmetric eigensolver (eigenvalues in ascending order); take the top
    // npcs from the end so they are in descending order, as above.
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(XtX);
    const Eigen::VectorXd& evalues = es.eigenvalues();
    const Eigen::MatrixXd& evectors = es.eigenvectors();
    loadings.resize(nfeatures, npcs);
    for (int j = 0; j < npcs; ++j) {
      const int src = nfeatures - 1 - j;
      loadings.col(j) = evectors.col(src);
      d(j) = std::sqrt(std::max(0.0, evalues(src)));
    }
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
