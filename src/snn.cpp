#include <RcppEigen.h>
#include <RcppThread.h>
#include <algorithm>
#include <limits>
#include <numeric>
#include <thread>
#include <utility>

using namespace Rcpp;

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppThread)]]

typedef Eigen::Triplet<double> T;

struct SNNOverlapWorker {
  const std::vector<int>& nn_buf;
  const std::vector<int>& offsets;
  const std::vector<int>& postings;
  const int n;
  const int k;
  const double k_d;
  const double prune;
  std::vector<int> overlap;
  std::vector<int> touched;
  std::vector<T> triplets;

  SNNOverlapWorker(
    const std::vector<int>& nn_buf,
    const std::vector<int>& offsets,
    const std::vector<int>& postings,
    const int n,
    const int k,
    const double prune
  ) : nn_buf(nn_buf), offsets(offsets), postings(postings), n(n), k(k), k_d(static_cast<double>(k)), prune(prune), overlap(n, 0) {
    touched.reserve(static_cast<size_t>(k) * k);
    triplets.reserve(static_cast<size_t>(k) * 64);
  }

  void operator()(std::size_t begin, std::size_t end) {
    for (std::size_t row = begin; row < end; ++row) {
      for (int col = 0; col < k; ++col) {
        const int neighbor = nn_buf[row + static_cast<std::size_t>(col) * n];

        for (int idx = offsets[neighbor]; idx < offsets[neighbor + 1]; ++idx) {
          const int other = postings[idx];

          // Remember the candidate the first time it is touched so it can be reset cheaply.
          if (overlap[other] == 0) {
            touched.push_back(other);
          }

          ++overlap[other];
        }
      }

      for (const int other : touched) {
        const double shared = static_cast<double>(overlap[other]);
        const double value = shared / (k_d + (k_d - shared));

        if (value >= prune) {
          triplets.emplace_back(static_cast<int>(row), other, value);
        }

        // Reset only the candidates touched for this row instead of clearing overlap[0:n].
        overlap[other] = 0;
      }

      touched.clear();
    }
  }

  void join(const SNNOverlapWorker& other) {
    triplets.insert(triplets.end(), other.triplets.begin(), other.triplets.end());
  }
};

// Compute a Jaccard-style SNN graph from a 1-based nearest-neighbor rank matrix.
// [[Rcpp::export(rng = false)]]
S4 ComputeSNN(IntegerMatrix nn_ranked, double prune, int nthreads) {
  const int n = nn_ranked.nrow();
  const int k = nn_ranked.ncol();

  // Allocate a plain C++ buffer so parallel workers never touch R's SEXP memory.
  std::vector<int> nn_buf(static_cast<size_t>(n) * k);

  // Copy the matrix to native memory and convert R's 1-based indices to 0-based.
  const int* src = INTEGER(nn_ranked);
  for (size_t i = 0; i < nn_buf.size(); ++i) {
    nn_buf[i] = src[i] - 1;
  }

  // Validate neighbor indices before entering parallel code so Rcpp::stop() is safe.
  for (size_t i = 0; i < nn_buf.size(); ++i) {
    if (nn_buf[i] < 0 || nn_buf[i] >= n) {
      stop("nn_ranked contains an index outside the valid cell range");
    }
  }

  // Build an inverted index from neighbor id to rows containing that neighbor.
  std::vector<int> offsets(n + 1, 0);
  for (int col = 0; col < k; ++col) {
    for (int row = 0; row < n; ++row) {
      ++offsets[nn_buf[row + static_cast<size_t>(col) * n] + 1];
    }
  }
  for (int cell = 0; cell < n; ++cell) {
    offsets[cell + 1] += offsets[cell];
  }
  std::vector<int> postings(static_cast<size_t>(n) * k);
  {
    std::vector<int> cursor(offsets.begin(), offsets.end());
    for (int col = 0; col < k; ++col) {
      for (int row = 0; row < n; ++row) {
        const int neighbor = nn_buf[row + static_cast<size_t>(col) * n];
        postings[cursor[neighbor]++] = row;
      }
    }
  }

  std::vector<T> triplets;

  // Generate overlap triplets, using chunked workers only for multi-threading.
  if (nthreads == 1) {
    SNNOverlapWorker worker(nn_buf, offsets, postings, n, k, prune);
    worker(0, n);
    if (worker.triplets.size() > static_cast<size_t>(std::numeric_limits<int>::max())) {
      stop("SNN graph has too many non-zero entries for a dgCMatrix");
    }
    triplets = std::move(worker.triplets);
  } else {
    const int chunks = std::max(1, std::min(n, nthreads * 8));
    std::vector<SNNOverlapWorker> workers;
    workers.reserve(chunks);
    for (int chunk = 0; chunk < chunks; ++chunk) {
      workers.emplace_back(nn_buf, offsets, postings, n, k, prune);
    }

    RcppThread::parallelFor(0, chunks, [&](int chunk) {
      const int begin = (n * chunk) / chunks;
      const int end = (n * (chunk + 1)) / chunks;
      workers[chunk](begin, end);
    }, nthreads);

    size_t triplet_count = 0;
    for (int chunk = 0; chunk < chunks; ++chunk) {
      triplet_count += workers[chunk].triplets.size();
    }
    if (triplet_count > static_cast<size_t>(std::numeric_limits<int>::max())) {
      stop("SNN graph has too many non-zero entries for a dgCMatrix");
    }
    triplets.reserve(triplet_count);
    for (int chunk = 0; chunk < chunks; ++chunk) {
      triplets.insert(triplets.end(), workers[chunk].triplets.begin(), workers[chunk].triplets.end());
    }
  }

  // Sort triplets before building the sparse matrix slots
  std::sort(
    triplets.begin(),
    triplets.end(),
    [](const T& lhs, const T& rhs) {
      return lhs.col() == rhs.col() ? lhs.row() < rhs.row() : lhs.col() < rhs.col();
    }
  );

  const int nnz = static_cast<int>(triplets.size());
  IntegerVector i(nnz);
  IntegerVector p(n + 1);
  NumericVector x(nnz);

  for (const T& entry : triplets) {
    ++p[entry.col() + 1];
  }
  for (int col = 0; col < n; ++col) {
    p[col + 1] += p[col];
  }

  std::vector<int> cursor(p.begin(), p.end());
  for (const T& entry : triplets) {
    const int offset = cursor[entry.col()]++;
    i[offset] = entry.row();
    x[offset] = entry.value();
  }

  S4 SNN("dgCMatrix");
  SNN.slot("i") = i;
  SNN.slot("p") = p;
  SNN.slot("Dim") = IntegerVector::create(n, n);
  SNN.slot("Dimnames") = List::create(R_NilValue, R_NilValue);
  SNN.slot("x") = x;
  SNN.slot("factors") = List::create();

  return SNN;
}

// [[Rcpp::export]]
std::vector<double> SNN_SmallestNonzero_Dist(
    Eigen::SparseMatrix<double> snn,
    Eigen::MatrixXd mat,
    int n,
    std::vector<double> nearest_dist
) {
  std::vector<double> results;
  results.reserve(snn.outerSize());
  std::vector<std::pair<double, size_t>> nonzero;
  std::vector<double> dists;
  for (int i=0; i < snn.outerSize(); ++i){
    // Store nonzero SNN edge weights with their row indices.
    nonzero.clear();
    for (Eigen::SparseMatrix<double>::InnerIterator it(snn, i); it; ++it) {
      nonzero.emplace_back(it.value(), it.row());
    }
    int n_i = n;
    if (n_i > nonzero.size()) n_i = nonzero.size();
    std::stable_sort(
      nonzero.begin(),
      nonzero.end(),
      [](const std::pair<double, size_t>& lhs, const std::pair<double, size_t>& rhs) {
        return lhs.first < rhs.first;
      }
    );
    dists.clear();
    for (size_t j = 0; j < nonzero.size(); ++j) {
      // compute euclidean distances to cells with small edge weights
      // if multiple entries have same value as nth element, calc dist to all
      size_t cell = nonzero[j].second;
      if(dists.size() < n_i  || nonzero[j].first == nonzero[n_i-1].first) {
        double res = (mat.row(cell) - mat.row(i)).norm();
        if (nearest_dist[i] > 0) {
          res = res - nearest_dist[i];
          if (res < 0) res = 0;
        }
        dists.push_back(res);
      } else {
        break;
      }
    }
    double avg_dist;
    if (dists.size() > n_i) {
      std::sort(dists.rbegin(), dists.rend());
      avg_dist = std::accumulate(dists.begin(), dists.begin() + n_i, 0.0) / n_i;
    } else {
      avg_dist = std::accumulate(dists.begin(), dists.end(), 0.0) / dists.size();
    }
    results.push_back(avg_dist);
  }
  return results;
}
