#include <RcppEigen.h>
#include <RcppThread.h>
// Spectra (header-only) provides the top-k symmetric eigensolver used below
#include <Spectra/SymEigsSolver.h>
#include <Spectra/MatOp/DenseSymMatProd.h>
#include <algorithm>
#include <cmath>
#include <limits>

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RSpectra)]]
// [[Rcpp::depends(RcppThread)]]

using namespace Rcpp;

// Fill the lower triangle of the Gram matrix X X' by output-column blocks.
// Spectra::DenseSymMatProd and SelfAdjointEigenSolver read this triangle below.
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

// Approximate PCA via the Gram matrix and a top-k symmetric eigendecomposition.
// `object` is the feature-by-cell scaled matrix X. The eigenvectors of X X' are
// feature loadings. Cell embeddings are X' U = V D, or V when weight_by_var is false.
//
// [[Rcpp::export(rng = false)]]
List EigenGramPCA(const Eigen::Map<Eigen::MatrixXd> object,
                  int npcs,
                  bool weight_by_var,
                  int nthreads = 1) {
  const int nfeatures = object.rows();
  const int ncells = object.cols();
  if (npcs < 1 || npcs >= std::min(nfeatures, ncells)) {
    stop("npcs must be positive and strictly less than both matrix dimensions.");
  }

  // Gram matrix X X' (features x features); only the lower triangle is formed.
  Eigen::MatrixXd XtX = Eigen::MatrixXd::Zero(nfeatures, nfeatures);
  if (nthreads <= 1) {
    XtX.selfadjointView<Eigen::Lower>().rankUpdate(object);
  } else {
    // Parallelize the triangular Gram fill by output-column blocks.
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

  // Use the partial Lanczos solver when it can compute fewer eigenpairs than
  // the full dense eigensolver; otherwise use the full solver.
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
    // Use SVD to recover orthonormal scores for near-zero singular values (R will catch the error and fall back to SVD)
    const double min_singular = d.maxCoeff() * std::sqrt(std::numeric_limits<double>::epsilon());
    if (d.minCoeff() <= min_singular) {
      stop("Unweighted PCA requires SVD for numerically zero singular values.");
    }
    for (int j = 0; j < npcs; ++j) {
      embeddings.col(j) /= d(j);
    }
  }

  Eigen::VectorXd sdev = d / std::sqrt(static_cast<double>(std::max(1, ncells - 1)));

  return List::create(
    _["loadings"] = loadings,
    _["embeddings"] = embeddings,
    _["sdev"] = sdev
  );
}
