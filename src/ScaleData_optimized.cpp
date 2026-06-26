#include <Rcpp.h>
#include <RcppParallel.h>
#include <progress.hpp>
#include <cmath>
#include <cstring>
#include <vector>

// [[Rcpp::depends(RcppProgress)]]
// [[Rcpp::depends(RcppParallel)]]

using namespace Rcpp;

inline double clip_scale_value(double value, double scale_max) {
  if (value > scale_max) {
    return scale_max;
  }
  return value;
}

// ---------------------------------------------------------------------------
// Sparse row scaling
// ---------------------------------------------------------------------------

// Scale a single column of the sparse output in place. Shared by the serial and
// parallel paths so the arithmetic is guaranteed to be identical regardless of
// threading.
static inline void scale_sparse_column(
  const int col,
  double* out,
  const int* i, const int* p, const double* x,
  const int* gmap,
  const double* mu, const double* inv_sigma, const double* zero_value,
  const char* valid,
  const int n_sel, const bool clip, const double scale_max
) {
  double* col_ptr = out + static_cast<R_xlen_t>(col) * n_sel;
  // Start every cell at the scaled value of a structural zero, then overwrite
  // the stored (non-zero) entries below.
  std::memcpy(col_ptr, zero_value, static_cast<size_t>(n_sel) * sizeof(double));
  for (int ptr = p[col]; ptr < p[col + 1]; ++ptr) {
    const int row = gmap[i[ptr]];
    // row < 0 marks a feature that was not requested; valid == 0 marks a row
    // with zero/NA sigma that must be left at zero (matches prior behaviour).
    if (row >= 0 && valid[row]) {
      double value = (x[ptr] - mu[row]) * inv_sigma[row];
      if (clip) {
        value = clip_scale_value(value, scale_max);
      }
      col_ptr[row] = value;
    }
  }
}

// RcppParallel worker: scales a contiguous range of columns. Only reads raw
// C pointers / std::vector data captured below, so it never touches the R API
// from a worker thread.
struct SparseScaleWorker : public RcppParallel::Worker {
  double* out;
  const int* i; const int* p; const double* x;
  const int* gmap;
  const double* mu; const double* inv_sigma; const double* zero_value;
  const char* valid;
  const int n_sel; const bool clip; const double scale_max;

  SparseScaleWorker(
    double* out, const int* i, const int* p, const double* x, const int* gmap,
    const double* mu, const double* inv_sigma, const double* zero_value,
    const char* valid, int n_sel, bool clip, double scale_max
  ) : out(out), i(i), p(p), x(x), gmap(gmap), mu(mu), inv_sigma(inv_sigma),
      zero_value(zero_value), valid(valid), n_sel(n_sel), clip(clip),
      scale_max(scale_max) {}

  void operator()(std::size_t begin, std::size_t end) {
    for (std::size_t col = begin; col < end; ++col) {
      scale_sparse_column(static_cast<int>(col), out, i, p, x, gmap, mu,
                          inv_sigma, zero_value, valid, n_sel, clip, scale_max);
    }
  }
};

// [[Rcpp::export(rng = false)]]
NumericMatrix FastSparseRowScale_optimized(
  NumericVector x,
  IntegerVector i,
  IntegerVector p,
  int rows,
  int cols,
  IntegerVector features = IntegerVector::create(),
  bool scale = true,
  bool center = true,
  double scale_max = 10,
  int nthreads = 1,
  bool display_progress = false
) {
  // `features` holds the 0-based row indices (of the full matrix) to scale, in
  // output order. When empty, every row is scaled in its original order; this
  // keeps callers that pass an already-subset matrix working unchanged.
  const bool all_rows = (features.size() == 0);
  const int n_sel = all_rows ? rows : static_cast<int>(features.size());

  // Map each full-matrix row to its output row (-1 = not requested). Selecting
  // features here lets the caller pass the full matrix and avoid an R-level
  // sparse-matrix subset copy.
  std::vector<int> gmap(rows, -1);
  if (all_rows) {
    for (int r = 0; r < rows; ++r) {
      gmap[r] = r;
    }
  } else {
    for (int g = 0; g < n_sel; ++g) {
      gmap[features[g]] = g;
    }
  }

  // First pass: per-feature sum and sum-of-squares over all stored values.
  std::vector<double> row_sum(n_sel, 0.0);
  std::vector<double> row_sq_sum(n_sel, 0.0);
  for (R_xlen_t idx = 0; idx < x.size(); ++idx) {
    const int row = gmap[i[idx]];
    if (row >= 0) {
      const double value = x[idx];
      row_sum[row] += value;
      row_sq_sum[row] += value * value;
    }
  }

  // Per-feature mean, inverse standard deviation and the scaled value of a
  // structural zero. `valid` records the rows we actually scale; rows with a
  // zero or NA sigma stay at zero, exactly as the original implementation did.
  std::vector<double> mu(n_sel);
  std::vector<double> inv_sigma(n_sel, 0.0);
  std::vector<double> zero_value(n_sel, 0.0);
  std::vector<char> valid(n_sel, 0);
  const bool clip = (scale_max != R_PosInf);
  for (int row = 0; row < n_sel; ++row) {
    const double mean = row_sum[row] / static_cast<double>(cols);
    mu[row] = center ? mean : 0.0;
    double sigma;
    if (scale) {
      const double variance_numerator = center
        ? row_sq_sum[row] - (row_sum[row] * row_sum[row] / static_cast<double>(cols))
        : row_sq_sum[row];
      sigma = std::sqrt(variance_numerator / static_cast<double>(cols - 1));
    } else {
      sigma = 1.0;
    }
    if (sigma > 0.0 && !R_IsNA(sigma)) {
      valid[row] = 1;
      inv_sigma[row] = 1.0 / sigma;  // multiply by reciprocal instead of divide
      zero_value[row] = (0.0 - mu[row]) * inv_sigma[row];
      if (clip) {
        zero_value[row] = clip_scale_value(zero_value[row], scale_max);
      }
    }
  }

  NumericMatrix out = no_init_matrix(n_sel, cols);
  double* out_ptr = REAL(out);
  const int* ip = INTEGER(i);
  const int* pp = INTEGER(p);
  const double* xp = REAL(x);

  if (nthreads <= 1) {
    // Serial path: show an RcppProgress bar (one tick per cell), matching the
    // progress bar ScaleData has historically printed.
    Progress progress(cols, display_progress);
    for (int col = 0; col < cols; ++col) {
      if (Progress::check_abort()) {
        return out;
      }
      scale_sparse_column(col, out_ptr, ip, pp, xp, gmap.data(), mu.data(),
                          inv_sigma.data(), zero_value.data(), valid.data(),
                          n_sel, clip, scale_max);
      progress.increment();
    }
  } else {
    // Parallel path (RcppParallel): no progress bar, because RcppProgress cannot
    // safely tick from worker threads. Columns are independent, so we simply
    // split the cell range across the requested number of threads.
    SparseScaleWorker worker(out_ptr, ip, pp, xp, gmap.data(), mu.data(),
                             inv_sigma.data(), zero_value.data(), valid.data(),
                             n_sel, clip, scale_max);
    const std::size_t grain = std::max<std::size_t>(
      1, static_cast<std::size_t>(cols) / (static_cast<std::size_t>(nthreads) * 8));
    RcppParallel::parallelFor(0, cols, worker, grain, nthreads);
  }
  return out;
}

// ---------------------------------------------------------------------------
// Dense row scaling
// ---------------------------------------------------------------------------

// Scale a single column of the dense output in place. `sel[row]` gives the
// source row in `mat` for output row `row`. Shared by serial and parallel paths.
static inline void scale_dense_column(
  const int col,
  double* out, const double* mat,
  const int* sel,
  const double* mu, const double* inv_sigma, const char* valid,
  const int n_sel, const int full_rows, const bool clip, const double scale_max
) {
  double* out_col = out + static_cast<R_xlen_t>(col) * n_sel;
  const double* mat_col = mat + static_cast<R_xlen_t>(col) * full_rows;
  for (int row = 0; row < n_sel; ++row) {
    double value = 0.0;
    if (valid[row]) {
      value = (mat_col[sel[row]] - mu[row]) * inv_sigma[row];
      if (clip) {
        value = clip_scale_value(value, scale_max);
      }
    }
    out_col[row] = value;
  }
}

struct DenseScaleWorker : public RcppParallel::Worker {
  double* out; const double* mat; const int* sel;
  const double* mu; const double* inv_sigma; const char* valid;
  const int n_sel; const int full_rows; const bool clip; const double scale_max;

  DenseScaleWorker(
    double* out, const double* mat, const int* sel, const double* mu,
    const double* inv_sigma, const char* valid, int n_sel, int full_rows,
    bool clip, double scale_max
  ) : out(out), mat(mat), sel(sel), mu(mu), inv_sigma(inv_sigma), valid(valid),
      n_sel(n_sel), full_rows(full_rows), clip(clip), scale_max(scale_max) {}

  void operator()(std::size_t begin, std::size_t end) {
    for (std::size_t col = begin; col < end; ++col) {
      scale_dense_column(static_cast<int>(col), out, mat, sel, mu, inv_sigma,
                         valid, n_sel, full_rows, clip, scale_max);
    }
  }
};

// [[Rcpp::export(rng = false)]]
NumericMatrix FastDenseRowScale_optimized(
  NumericMatrix mat,
  IntegerVector features = IntegerVector::create(),
  bool scale = true,
  bool center = true,
  double scale_max = 10,
  int nthreads = 1,
  bool display_progress = false
) {
  const int full_rows = mat.nrow();
  const int cols = mat.ncol();
  // `features` holds the 0-based rows of `mat` to scale, in output order. When
  // empty, every row is scaled in its original order.
  const bool all_rows = (features.size() == 0);
  const int n_sel = all_rows ? full_rows : static_cast<int>(features.size());
  const double* mat_ptr = REAL(mat);

  // Forward map output row -> source row in `mat`.
  std::vector<int> sel(n_sel);
  if (all_rows) {
    for (int r = 0; r < n_sel; ++r) {
      sel[r] = r;
    }
  } else {
    for (int g = 0; g < n_sel; ++g) {
      sel[g] = features[g];
    }
  }

  std::vector<double> row_sum(n_sel, 0.0);
  std::vector<double> row_sq_sum(n_sel, 0.0);
  for (int col = 0; col < cols; ++col) {
    const double* mat_col = mat_ptr + static_cast<R_xlen_t>(col) * full_rows;
    for (int row = 0; row < n_sel; ++row) {
      const double value = mat_col[sel[row]];
      row_sum[row] += value;
      row_sq_sum[row] += value * value;
    }
  }

  std::vector<double> mu(n_sel);
  std::vector<double> inv_sigma(n_sel, 0.0);
  std::vector<char> valid(n_sel, 0);
  const bool clip = (scale_max != R_PosInf);
  for (int row = 0; row < n_sel; ++row) {
    const double sum = row_sum[row];
    const double mean = sum / static_cast<double>(cols);
    mu[row] = center ? mean : 0.0;
    double sigma;
    if (scale) {
      const double numerator = center
        ? row_sq_sum[row] - (sum * sum / static_cast<double>(cols))
        : row_sq_sum[row];
      sigma = std::sqrt(numerator / static_cast<double>(cols - 1));
    } else {
      sigma = 1.0;
    }
    if (sigma > 0.0 && !R_IsNA(sigma)) {
      valid[row] = 1;
      inv_sigma[row] = 1.0 / sigma;  // multiply by reciprocal instead of divide
    }
  }

  NumericMatrix out = no_init_matrix(n_sel, cols);
  double* out_ptr = REAL(out);

  if (nthreads <= 1) {
    // Serial path: show an RcppProgress bar (one tick per cell).
    Progress progress(cols, display_progress);
    for (int col = 0; col < cols; ++col) {
      if (Progress::check_abort()) {
        return out;
      }
      scale_dense_column(col, out_ptr, mat_ptr, sel.data(), mu.data(),
                         inv_sigma.data(), valid.data(), n_sel, full_rows,
                         clip, scale_max);
      progress.increment();
    }
  } else {
    // Parallel path (RcppParallel): no progress bar; split columns across threads.
    DenseScaleWorker worker(out_ptr, mat_ptr, sel.data(), mu.data(),
                            inv_sigma.data(), valid.data(), n_sel, full_rows,
                            clip, scale_max);
    const std::size_t grain = std::max<std::size_t>(
      1, static_cast<std::size_t>(cols) / (static_cast<std::size_t>(nthreads) * 8));
    RcppParallel::parallelFor(0, cols, worker, grain, nthreads);
  }
  return out;
}
