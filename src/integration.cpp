#include <RcppEigen.h>
#include <RcppThread.h>
#include <R_ext/BLAS.h>
#include <progress.hpp>
#include <algorithm>
#include <cmath>
#include <numeric>
#include <unordered_map>
#include <vector>
#include "data_manipulation.h"

using namespace Rcpp;
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppProgress)]]
// [[Rcpp::depends(RcppThread)]]

// Matrix-vector product used by the implicit CCA SVD path. Given two
// feature-by-cell matrices A and B and a vector x over B's cells, compute
// A' * (B * x) without returning the intermediate feature vector to R.
//
// Reference for calling BLAS from C via Fortran: https://cran.r-project.org/doc/manuals/r-release/R-exts.html#Calling-C-from-Fortran-and-vice-versa-1
//
// [[Rcpp::export(rng = false)]]
NumericVector CcaCrossprodMultiply(const Eigen::Map<Eigen::MatrixXd> left,
                                   const Eigen::Map<Eigen::MatrixXd> right,
                                   const Eigen::Map<Eigen::VectorXd> x) {
  const int nfeatures = left.rows();
  const int ncells = left.cols();
  if (right.rows() != nfeatures) {
    stop("CCA matrices must have the same number of rows.");
  }
  if (right.cols() != x.size()) {
    stop("CCA multiplier length must match the number of columns in the right matrix.");
  }

  const double one = 1.0;
  const double zero = 0.0;
  const int inc = 1;
  const int right_cells = right.cols();

  NumericVector feature_scores(no_init(nfeatures));
  NumericVector out(no_init(ncells));

  // BLAS is called through Fortran, so pass scalar args as pointers
  // feature_scores = right %*% x
  F77_CALL(dgemv)(
    "N",
    &nfeatures,
    &right_cells,
    &one,
    right.data(),
    &nfeatures,
    x.data(),
    &inc,
    &zero,
    REAL(feature_scores),
    &inc FCONE // required when passing char* to Fortran, see R-exts manual
  );

  // out = t(left) %*% feature_scores
  F77_CALL(dgemv)(
    "T",
    &nfeatures,
    &ncells,
    &one,
    left.data(),
    &nfeatures,
    REAL(feature_scores),
    &inc,
    &zero,
    REAL(out),
    &inc FCONE
  );

  return out;
}

typedef Eigen::Triplet<double> T;

// Use at most one work chunk per thread
static inline int integration_thread_chunk_count(const int nitems, const int nthreads) {
  return std::max(1, std::min(nitems, nthreads));
}

// Count how many unique nearest neighbors are shared by each anchor pair.
//
// Neighbor tables and anchor ids come from R as 1-based indices. B-neighbor ids
// are shifted by offset into a combined A+B id space
// [[Rcpp::export(rng = false)]]
IntegerVector CountAnchorSharedNeighbors(
    IntegerMatrix indices_aa,
    IntegerMatrix indices_ab,
    IntegerMatrix indices_ba,
    IntegerMatrix indices_bb,
    IntegerVector anchor_cell1,
    IntegerVector anchor_cell2,
    int offset,
    int k_score,
    int nthreads = 1
) {
  const int nanchors = anchor_cell1.size();
  std::vector<int> scores(nanchors);
  const int max_neighbor = offset + indices_bb.nrow();
  const int aa_rows = indices_aa.nrow();
  const int ab_rows = indices_ab.nrow();
  const int ba_rows = indices_ba.nrow();
  const int bb_rows = indices_bb.nrow();

  // Move R inputs into C++ storage before threaded work
  std::vector<int> aa(indices_aa.begin(), indices_aa.end());
  std::vector<int> ab(indices_ab.begin(), indices_ab.end());
  std::vector<int> ba(indices_ba.begin(), indices_ba.end());
  std::vector<int> bb(indices_bb.begin(), indices_bb.end());
  std::vector<int> cell1(anchor_cell1.begin(), anchor_cell1.end());
  std::vector<int> cell2(anchor_cell2.begin(), anchor_cell2.end());

  auto score_range = [&](std::size_t begin, std::size_t end) {
    // Reuse marker vectors within each worker
    std::vector<int> neighbors_a(max_neighbor + 1, 0);
    std::vector<int> seen_b(max_neighbor + 1, 0);

    for (std::size_t anchor = begin; anchor < end; ++anchor) {
      // The stamp marks entries that belong to the current anchor
      const int stamp = static_cast<int>(anchor) + 1;

      const int anchor_cell_a = cell1[anchor] - 1;
      const int anchor_cell_b = cell2[anchor] - 1;
      // Mark neighbors for the first anchor cell across both datasets
      for (int k = 0; k < k_score; ++k) {
        neighbors_a[aa[anchor_cell_a + static_cast<std::size_t>(k) * aa_rows]] = stamp;
        neighbors_a[ab[anchor_cell_a + static_cast<std::size_t>(k) * ab_rows] + offset] = stamp;
      }

      int shared = 0;
      // Count unique neighbors shared with the second anchor cell
      // seen_b ensures duplicate neighbor ids count only once per anchor
      for (int k = 0; k < k_score; ++k) {
        const int neighbor_ba = ba[anchor_cell_b + static_cast<std::size_t>(k) * ba_rows];
        if (seen_b[neighbor_ba] != stamp) {
          seen_b[neighbor_ba] = stamp;
          if (neighbors_a[neighbor_ba] == stamp) {
            ++shared;
          }
        }
        const int neighbor_bb = bb[anchor_cell_b + static_cast<std::size_t>(k) * bb_rows] + offset;
        if (seen_b[neighbor_bb] != stamp) {
          seen_b[neighbor_bb] = stamp;
          if (neighbors_a[neighbor_bb] == stamp) {
            ++shared;
          }
        }
      }
      scores[anchor] = shared;
    }
  };

  if (nthreads <= 1 || nanchors < 2) {
    score_range(0, nanchors);
  } else {
    const int chunks = integration_thread_chunk_count(nanchors, nthreads);
    RcppThread::parallelFor(0, chunks, [&](int chunk) {
      const int begin = (nanchors * chunk) / chunks;
      const int end = (nanchors * (chunk + 1)) / chunks;
      score_range(begin, end);
    }, nthreads);
  }
  return wrap(scores);
}

// Construct the anchor-by-cell weight matrix used during integration.
// [[Rcpp::export(rng = false)]]
Eigen::SparseMatrix<double> FindWeightsC(
    NumericVector cells2,
    Eigen::MatrixXd distances,
    std::vector<std::string> anchor_cells2,
    std::vector<std::string> integration_matrix_rownames,
    Eigen::MatrixXd cell_index,
    Eigen::VectorXd anchor_score,
    double min_dist,
    double sd,
    bool display_progress,
    int nthreads = 1
) {
  std::vector<T> tripletList;
  tripletList.reserve(anchor_cells2.size() * 10);
  std::vector<std::vector<int>> cell_map(anchor_cells2.size());

  const bool run_parallel_cells = nthreads > 1 && cells2.size() >= 2;
  Progress p(anchor_cells2.size() + (!run_parallel_cells ? cells2.size() : 0), display_progress);

  // Duplicate row names are expected when multiple anchors map to the same cell
  // Index the integration rows by anchor cell id for fast lookup
  std::unordered_map<std::string, std::vector<int>> integration_rows;
  integration_rows.reserve(integration_matrix_rownames.size());
  for (int i = 0; i < integration_matrix_rownames.size(); ++i) {
    integration_rows[integration_matrix_rownames[i]].push_back(i);
  }
  // build map from anchor_cells2 to integration_matrix rows
  for (int i = 0; i < anchor_cells2.size(); ++i) {
    const auto iter = integration_rows.find(anchor_cells2[i]);
    if (iter != integration_rows.end()) {
      cell_map[i] = iter->second;
    }
    p.increment();
  }

  std::vector<int> cells(cells2.size());
  for (int i = 0; i < cells2.size(); ++i) {
    cells[i] = static_cast<int>(cells2[i]);
  }
  const double scale = std::pow(2 / sd, 2);

  // Construct sparse distance-weight entries per cell
  auto weights_range = [&](std::size_t begin, std::size_t end, std::vector<T>& triplets, Progress* progress = NULL, RcppThread::ProgressBar* bar = NULL) {
    for (std::size_t cell_pos = begin; cell_pos < end; ++cell_pos) {
      const int cell = cells[cell_pos];
      int k = 0;
      for (int i = 0; i < cell_index.cols() && k < cell_index.cols(); ++i) {
        const int anchor_idx = static_cast<int>(cell_index(cell, i)) - 1;
        const std::vector<int>& mnn_idx = cell_map[anchor_idx];
        for (int j = 0; j < mnn_idx.size() && k < cell_index.cols(); ++j) {
          double to_add = 1 - std::exp(-1 * distances(cell, i) * anchor_score[mnn_idx[j]] / scale);
          triplets.push_back(T(mnn_idx[j], cell, to_add));
          k++;
        }
      }
      if (bar != NULL) {
        (*bar)++;
      }
      if (progress != NULL) {
        progress->increment();
      }
    }
  };

  if (!run_parallel_cells) {
    weights_range(0, cells.size(), tripletList, &p);
  } else { // Each parallel worker writes triplets into its own vector
    const int chunks = integration_thread_chunk_count(static_cast<int>(cells.size()), nthreads);
    std::vector<std::vector<T>> chunk_triplets(chunks);
    if (display_progress) {
      RcppThread::ProgressBar bar(cells.size(), 1);
      RcppThread::parallelFor(0, chunks, [&](int chunk) {
        const int begin = (static_cast<int>(cells.size()) * chunk) / chunks;
        const int end = (static_cast<int>(cells.size()) * (chunk + 1)) / chunks;
        chunk_triplets[chunk].reserve(static_cast<size_t>(end - begin) * cell_index.cols());
        weights_range(begin, end, chunk_triplets[chunk], NULL, &bar);
      }, nthreads);
    } else {
      RcppThread::parallelFor(0, chunks, [&](int chunk) {
        const int begin = (static_cast<int>(cells.size()) * chunk) / chunks;
        const int end = (static_cast<int>(cells.size()) * (chunk + 1)) / chunks;
        chunk_triplets[chunk].reserve(static_cast<size_t>(end - begin) * cell_index.cols());
        weights_range(begin, end, chunk_triplets[chunk]);
      }, nthreads);
    }

    size_t triplet_count = 0;
    for (const auto& chunk : chunk_triplets) {
      triplet_count += chunk.size();
    }
    tripletList.reserve(triplet_count);
    // Merge per-chunk triplets after workers finish
    for (auto& chunk : chunk_triplets) {
      tripletList.insert(tripletList.end(), chunk.begin(), chunk.end());
    }
  }
  Eigen::SparseMatrix<double> return_mat;
  if(min_dist == 0){
    Eigen::SparseMatrix<double> dist_weights(integration_matrix_rownames.size(), cells2.size());
    dist_weights.setFromTriplets(tripletList.begin(), tripletList.end(), [] (const double&, const double &b) { return b; });
    Eigen::VectorXd colSums = dist_weights.transpose() * Eigen::VectorXd::Ones(dist_weights.rows());
    for (int k=0; k < dist_weights.outerSize(); ++k){
      for (Eigen::SparseMatrix<double>::InnerIterator it(dist_weights, k); it; ++it){
        it.valueRef() = it.value()/colSums[k];
      }
    }
    return_mat = dist_weights;
  } else {
    Eigen::MatrixXd dist_weights = Eigen::MatrixXd::Constant(integration_matrix_rownames.size(), cells2.size(), min_dist);
    // Start with dense baseline weights from the minimum distance
    auto dense_range = [&](std::size_t begin, std::size_t end) {
      for (std::size_t i = begin; i < end; ++i) {
        for (int j = 0; j < dist_weights.rows(); ++j) {
          dist_weights(j, i) = 1 - std::exp(-1 * dist_weights(j, i) * anchor_score[j] / scale);
        }
      }
    };
    if (nthreads <= 1 || dist_weights.cols() < 2) {
      dense_range(0, dist_weights.cols());
    } else {
      const int chunks = integration_thread_chunk_count(static_cast<int>(dist_weights.cols()), nthreads);
      RcppThread::parallelFor(0, chunks, [&](int chunk) {
        const int begin = (dist_weights.cols() * chunk) / chunks;
        const int end = (dist_weights.cols() * (chunk + 1)) / chunks;
        dense_range(begin, end);
      }, nthreads);
    }
    for(auto const &weight : tripletList){
      dist_weights(weight.row(), weight.col()) = weight.value();
    }
    Eigen::VectorXd colSums = dist_weights.colwise().sum();
    auto normalize_range = [&](std::size_t begin, std::size_t end) {
      for (std::size_t i = begin; i < end; ++i) {
        for (int j = 0; j < dist_weights.rows(); ++j) {
          dist_weights(j, i) = dist_weights(j, i) / colSums[i];
        }
      }
    };
    if (nthreads <= 1 || dist_weights.cols() < 2) {
      normalize_range(0, dist_weights.cols());
    } else {
      const int chunks = integration_thread_chunk_count(static_cast<int>(dist_weights.cols()), nthreads);
      RcppThread::parallelFor(0, chunks, [&](int chunk) {
        const int begin = (dist_weights.cols() * chunk) / chunks;
        const int end = (dist_weights.cols() * (chunk + 1)) / chunks;
        normalize_range(begin, end);
      }, nthreads);
    }
    return_mat = dist_weights.sparseView();
  }
  return(return_mat);
}

// Apply the integration correction to the query expression matrix
// [[Rcpp::export(rng = false)]]
Eigen::SparseMatrix<double> IntegrateDataC(
    Eigen::SparseMatrix<double> integration_matrix,
    Eigen::SparseMatrix<double> weights,
    Eigen::SparseMatrix<double> expression_cells2
) {
  Eigen::SparseMatrix<double> corrected = expression_cells2 - weights.transpose() * integration_matrix;
  return(corrected);
}

// Return stable ascending indices for a vector without modifying the vector
template <typename S>
std::vector<size_t> sort_indexes(const std::vector<S> &v) {
  std::vector<size_t> idx(v.size());
  std::iota(idx.begin(), idx.end(), 0);
  std::stable_sort(idx.begin(), idx.end(),
                   [&v](size_t i1, size_t i2) {return v[i1] < v[i2];});
  return idx;
}

// Compute mapping scores for each query cell after integration
//
// query_pca is dimensions x cells
// query_dists and corrected_nns are cells x neighbor rank
// corrected_nns entries are 1-based ids (from R)
// [[Rcpp::export]]
std::vector<double> ScoreHelper(
    Eigen::SparseMatrix<double> snn,
    Eigen::MatrixXd query_pca,
    Eigen::MatrixXd query_dists,
    Eigen::MatrixXd corrected_nns,
    int k_snn,
    bool subtract_first_nn,
    bool display_progress,
    int nthreads = 1
) {
  std::vector<double> scores(snn.outerSize());

  const bool run_parallel = nthreads > 1 && snn.outerSize() >= 2;
  Progress p(snn.outerSize(), display_progress && !run_parallel);
  auto score_range = [&](std::size_t begin, std::size_t end, Progress* progress = NULL, RcppThread::ProgressBar* bar = NULL) {
    for (std::size_t i = begin; i < end; ++i) {
      // create vectors to store the nonzero snn elements and their indices
      std::vector<double> nonzero;
      std::vector<size_t> nonzero_idx;
      for (Eigen::SparseMatrix<double>::InnerIterator it(snn, i); it; ++it) {
        nonzero.push_back(it.value());
        nonzero_idx.push_back(it.row());
      }
      // find the k_snn cells with the smallest non-zero edge weights to use in
      // computing the transition probability bandwidth
      std::vector<size_t> nonzero_order = sort_indexes(nonzero);
      std::vector<double> bw_dists;
      int k_snn_i = k_snn;
      if (k_snn_i > nonzero_order.size()) k_snn_i = nonzero_order.size();
      for (int j = 0; j < nonzero_order.size(); ++j) {
        // compute euclidean distances to cells with small edge weights
        size_t cell = nonzero_idx[nonzero_order[j]];
        if(bw_dists.size() < k_snn_i || nonzero[nonzero_order[j]] == nonzero[nonzero_order[k_snn_i-1]]) {
          double res = (query_pca.col(cell) - query_pca.col(i)).norm();
          bw_dists.push_back(res);
        } else {
          break;
        }
      }
      // compute bandwidth as the mean distance of the farthest k_snn cells
      double bw;
      if (bw_dists.size() > k_snn_i) {
        std::sort(bw_dists.rbegin(), bw_dists.rend());
        bw = std::accumulate(bw_dists.begin(), bw_dists.begin() + k_snn_i, 0.0) / k_snn_i;
      } else {
        bw = std::accumulate(bw_dists.begin(), bw_dists.end(), 0.0) / bw_dists.size();
      }
      // compute transition probabilites
      double first_neighbor_dist;
      // subtract off distance to first neighbor?
      if (subtract_first_nn) {
        first_neighbor_dist = query_dists(i, 1);
      } else {
        first_neighbor_dist = 0;
      }
      bw = bw - first_neighbor_dist;
      double q_tps = 0;
      for(int j = 0; j < query_dists.cols(); ++j) {
        q_tps += std::exp(-1 * (query_dists(i, j) - first_neighbor_dist) / bw);
      }
      q_tps = q_tps/(query_dists.cols());
      double c_tps = 0;
      for(int j = 0; j < corrected_nns.cols(); ++j) {
        c_tps += exp(-1 * ((query_pca.col(i) - query_pca.col(corrected_nns(i, j)-1)).norm() - first_neighbor_dist) / bw);
      }
      c_tps = c_tps/(corrected_nns.cols());
      scores[i] = c_tps / q_tps;
      if (bar != NULL) {
        (*bar)++;
      }
      if (progress != NULL) {
        progress->increment();
      }
    }
  };

  if (!run_parallel) {
    score_range(0, snn.outerSize(), &p);
  } else {
    const int chunks = integration_thread_chunk_count(static_cast<int>(snn.outerSize()), nthreads);
    if (display_progress) {
      RcppThread::ProgressBar bar(snn.outerSize(), 1);
      RcppThread::parallelFor(0, chunks, [&](int chunk) {
        const int begin = (snn.outerSize() * chunk) / chunks;
        const int end = (snn.outerSize() * (chunk + 1)) / chunks;
        score_range(begin, end, NULL, &bar);
      }, nthreads);
    } else {
      RcppThread::parallelFor(0, chunks, [&](int chunk) {
        const int begin = (snn.outerSize() * chunk) / chunks;
        const int end = (snn.outerSize() * (chunk + 1)) / chunks;
        score_range(begin, end);
      }, nthreads);
    }
  }
  return(scores);
}
