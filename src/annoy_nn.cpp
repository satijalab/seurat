// Threaded search of an Annoy index built by RcppAnnoy.
#include <RcppAnnoy.h>
#include <RcppThread.h>

#include <algorithm>
#include <atomic>
#include <cstdlib>
#include <cstdint>
#include <string>
#include <type_traits>
#include <vector>

typedef std::remove_const<decltype(Kiss64Random::default_seed)>::type SeedType;

// The saved index is loaded through Annoy's public interface
// Hamming indexes use uint64_t coordinates; the other supported metrics use float
template <typename Coord>
using AnnoyInterface = Annoy::AnnoyIndexInterface<int32_t, Coord, SeedType>;

template <typename Coord>
static Rcpp::List SearchIndex(AnnoyInterface<Coord>* index, const Rcpp::NumericMatrix& query, int k, int search_k, bool include_distance, bool angular, int nthreads) {
  const int n = query.nrow(), f = query.ncol();

  // Pack the column-major R matrix into row-major query vectors for Annoy.
  std::vector<Coord> points(static_cast<std::size_t>(n) * f);
  for (int j = 0; j < f; ++j) {
    for (int i = 0; i < n; ++i) {
      points[static_cast<std::size_t>(i) * f + j] =
        static_cast<Coord>(query[static_cast<std::size_t>(j) * n + i]);
    }
  }

  // Allocate R objects before starting worker threads. The worker loop writes
  // into raw matrix storage and does not call back into the R API.
  Rcpp::NumericMatrix nn_idx(n, k);
  Rcpp::NumericMatrix nn_dist(n, k);
  double* idx_out = nn_idx.begin();
  double* dist_out = nn_dist.begin();
  const std::size_t rows = static_cast<std::size_t>(n);
  std::atomic<int> answered(0);

  RcppThread::parallelFor(0, n, [&](int i) {
    std::vector<int32_t> items;
    std::vector<Coord> distances;
    index->get_nns_by_vector(
      &points[static_cast<std::size_t>(i) * f], k, search_k, &items,
      include_distance ? &distances : NULL
    );
    // Fill any missing neighbor slots with NA.
    const int found = std::min(static_cast<int>(items.size()), k);
    for (int j = 0; j < found; ++j) {
      idx_out[static_cast<std::size_t>(j) * rows + i] = items[j] + 1;
      if (include_distance) {
        // Convert Annoy angular distance to cosine distance.
        const double d = static_cast<double>(distances[j]);
        dist_out[static_cast<std::size_t>(j) * rows + i] = angular ? 0.5 * (d * d) : d;
      } else {
        dist_out[static_cast<std::size_t>(j) * rows + i] = NA_REAL;
      }
    }
    for (int j = found; j < k; ++j) {
      idx_out[static_cast<std::size_t>(j) * rows + i] = NA_REAL;
      dist_out[static_cast<std::size_t>(j) * rows + i] = NA_REAL;
    }
    answered.fetch_add(1, std::memory_order_relaxed);
  }, nthreads);

  // If an interrupt stops the parallel loop early, don't return partial results
  if (answered.load(std::memory_order_relaxed) != n) {
    Rcpp::stop("Annoy neighbor search was interrupted");
  }

  return Rcpp::List::create(
    Rcpp::_["nn.idx"] = nn_idx,
    Rcpp::_["nn.dists"] = nn_dist
  );
}

// Load an Annoy index file and search it in parallel
/*
* Credit for file-backed parallel search idea: @jlmelville / uwot package
*/
template <typename Distance, typename Coord>
static SEXP LoadAndSearchIndex(const std::string& index_path, const Rcpp::NumericMatrix& query, int k, int search_k, bool include_distance, bool angular, int nthreads) {
  // Construct the index with the query dimensionality
  Annoy::AnnoyIndex<int32_t, Coord, Distance, Kiss64Random, RcppAnnoyIndexThreadPolicy> index(query.ncol());
  char* error = NULL;
  const bool loaded = index.load(index_path.c_str(), false, &error);
  if (!loaded) {
    std::string message = "Unable to load Annoy index";
    if (error != NULL) {
      message += ": ";
      message += error;
      std::free(error);
    }
    Rcpp::stop(message);
  }
  if (index.get_f() != query.ncol()) {
    Rcpp::stop("Loaded Annoy index dimensionality does not match query dimensionality");
  }
  if (index.get_n_items() <= 0) {
    Rcpp::stop("Loaded Annoy index contains no items");
  }
  return SearchIndex<Coord>(&index, query, k, search_k, include_distance, angular, nthreads);
}

// [[Rcpp::export(rng = false)]]
SEXP AnnoySearchCpp(std::string index_path, Rcpp::NumericMatrix query, int k, int search_k, bool include_distance, std::string metric, int nthreads) {
  // Dispatch to the Annoy distance type that matches the index file written by RcppAnnoy
  if (metric == "euclidean") {
    return LoadAndSearchIndex<Annoy::Euclidean, float>(index_path, query, k, search_k, include_distance, false, nthreads);
  } else if (metric == "cosine") {  // SearchIndex converts Annoy angular distance to cosine distance
    return LoadAndSearchIndex<Annoy::Angular, float>(index_path, query, k, search_k, include_distance, true, nthreads);
  } else if (metric == "manhattan") {
    return LoadAndSearchIndex<Annoy::Manhattan, float>(index_path, query, k, search_k, include_distance, false, nthreads);
  } else if (metric == "hamming") {
    return LoadAndSearchIndex<Annoy::Hamming, uint64_t>(index_path, query, k, search_k, include_distance, false, nthreads);
  }
  Rcpp::stop("Unsupported Annoy metric");
}
