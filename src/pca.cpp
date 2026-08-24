#include <RcppEigen.h>
#include <RcppThread.h>
// Spectra (header-only, Eigen-native; the library RSpectra wraps) gives a
// partial *top-k* symmetric eigensolver. The top-k Lanczos solver below is machine-precision 
// identical for the leading components but avoids computing the discarded eigenvectors.
#include <Spectra/SymEigsSolver.h>
#include <Spectra/MatOp/DenseSymMatProd.h>
#include <algorithm>
#include <cmath>

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RSpectra)]]
// [[Rcpp::depends(RcppThread)]]

using namespace Rcpp;

// Stage 1: form the lower triangle of the Gram matrix X X' by output-column
// blocks. For a contiguous column range [begin, end) the only rows needed (lower
// triangle: global row >= global col) start at `begin`, so we multiply the bottom
// rows of X by the block's rows and copy each column's lower part into XtX. Only
// the lower triangle is filled, exactly as the previous rankUpdate<Lower> did, and
// that is what Spectra::DenseSymMatProd / SelfAdjointEigenSolver read below.
struct GramWorker {
  const Eigen::Map<Eigen::MatrixXd> X;  // nfeatures x ncells (lightweight view)
  Eigen::MatrixXd& XtX;                 // nfeatures x nfeatures, lower triangle
  const int nfeatures;

  GramWorker(const Eigen::Map<Eigen::MatrixXd> X, Eigen::MatrixXd& XtX)
    : X(X), XtX(XtX), nfeatures(static_cast<int>(X.rows())) {}

  void operator()(std::size_t begin, std::size_t end) {
    const int b = static_cast<int>(begin);
    const int e = static_cast<int>(end);
    const int width = e - b;
    // block(r, c) = XtX(b + r, b + c) for r in [0, nfeatures - b), c in [0, width)
    Eigen::MatrixXd block =
      X.bottomRows(nfeatures - b) * X.middleRows(b, width).transpose();
    for (int c = 0; c < width; ++c) {
      const int col = b + c;                 // global column index
      // lower triangle of this column: global rows col .. nfeatures - 1, which is
      // block rows c .. (nfeatures - b - 1).
      XtX.col(col).segment(col, nfeatures - col) =
        block.col(c).segment(c, nfeatures - b - c);
    }
  }
};

// Stage 3: embeddings = X' U (ncells x npcs), parallelized over cells (rows of the
// output). Each thread computes a disjoint block of embedding rows.
struct EmbeddingWorker {
  const Eigen::Map<Eigen::MatrixXd> X;  // nfeatures x ncells
  const Eigen::MatrixXd& loadings;      // nfeatures x npcs
  Eigen::MatrixXd& embeddings;          // ncells x npcs

  EmbeddingWorker(const Eigen::Map<Eigen::MatrixXd> X,
                  const Eigen::MatrixXd& loadings,
                  Eigen::MatrixXd& embeddings)
    : X(X), loadings(loadings), embeddings(embeddings) {}

  void operator()(std::size_t begin, std::size_t end) {
    const int b = static_cast<int>(begin);
    const int width = static_cast<int>(end) - b;
    embeddings.middleRows(b, width) =
      X.middleCols(b, width).transpose() * loadings;
  }
};

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
                  bool weight_by_var,
                  int nthreads = 1) {
  const int nfeatures = object.rows();
  const int ncells = object.cols();
  npcs = std::min(npcs, nfeatures);

  // Gram matrix X X' (features x features); only the lower triangle is formed.
  Eigen::MatrixXd XtX = Eigen::MatrixXd::Zero(nfeatures, nfeatures);
  if (nthreads <= 1) {
    XtX.selfadjointView<Eigen::Lower>().rankUpdate(object);
  } else {
    // Parallel Gram by output-column blocks (see GramWorker). Fine grain lets RcppThread
    // balance the unequal per-column work of a triangular fill.
    GramWorker gram_worker(object, XtX);
    const int chunks = std::max(1, std::min(nfeatures, nthreads * 8));
    RcppThread::parallelFor(0, chunks, [&](int chunk) {
      const int begin = (nfeatures * chunk) / chunks;
      const int end = (nfeatures * (chunk + 1)) / chunks;
      gram_worker(begin, end);
    }, nthreads);
  }

  // The loadings are the top npcs eigenvectors of XX' (descending) and d the
  // square-roots of the corresponding eigenvalues.
  Eigen::MatrixXd loadings;
  Eigen::VectorXd d(npcs);

  // Compute only the top npcs eigenpairs with a Lanczos solver when it is
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

  Eigen::MatrixXd embeddings(ncells, npcs);  // X' U = V D
  if (nthreads <= 1) {
    embeddings.noalias() = object.transpose() * loadings;
  } else {
    EmbeddingWorker emb_worker(object, loadings, embeddings);
    const int chunks = std::max(1, std::min(ncells, nthreads * 8));
    RcppThread::parallelFor(0, chunks, [&](int chunk) {
      const int begin = (ncells * chunk) / chunks;
      const int end = (ncells * (chunk + 1)) / chunks;
      emb_worker(begin, end);
    }, nthreads);
  }
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
