#include <Rcpp.h>
#include <RcppThread.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

using namespace Rcpp;

// [[Rcpp::depends(RcppThread)]]

inline double sct_model_var(double mu, double theta) {
  if (R_finite(theta)) {
    return mu + (mu * mu / theta);
  }
  return mu;
}

inline double sct_clip(double value, double lower, double upper) {
  if (value < lower) {
    return lower;
  }
  if (value > upper) {
    return upper;
  }
  return value;
}

inline double sct_round0(double value) {
  return std::nearbyint(value);
}

inline S4 sct_corrected_matrix(
  const std::vector<int>& corrected_i,
  const std::vector<int>& corrected_p,
  const std::vector<double>& corrected_x,
  int rows,
  int cols
) {
  S4 corrected("dgCMatrix");
  corrected.slot("i") = wrap(corrected_i);
  corrected.slot("p") = wrap(corrected_p);
  corrected.slot("x") = wrap(corrected_x);
  corrected.slot("Dim") = IntegerVector::create(rows, cols);
  corrected.slot("Dimnames") = List::create(R_NilValue, R_NilValue);
  return corrected;
}

// Reducer for per-gene residual mean/variance across a range of cells (columns).
// Writes no shared output; accumulates per-row sums that are merged in join().
struct SCTStatsReducer {
  const int* p;
  const int* i;
  const double* x;
  const int rows;
  const double* log_umi_ptr;
  const double* theta_ptr;
  const double* slope_ptr;
  const std::vector<double>& exp_intercept;
  const bool common_slope;
  const double first_slope;
  const double min_var;
  const double residual_clip_min;
  const double residual_clip_max;
  std::vector<double> residual_sum;
  std::vector<double> residual_sq_sum;

  SCTStatsReducer(
    const int* p,
    const int* i,
    const double* x,
    int rows,
    const double* log_umi_ptr,
    const double* theta_ptr,
    const double* slope_ptr,
    const std::vector<double>& exp_intercept,
    bool common_slope,
    double first_slope,
    double min_var,
    double residual_clip_min,
    double residual_clip_max
  ) : p(p), i(i), x(x), rows(rows), log_umi_ptr(log_umi_ptr),
      theta_ptr(theta_ptr), slope_ptr(slope_ptr), exp_intercept(exp_intercept),
      common_slope(common_slope), first_slope(first_slope), min_var(min_var),
      residual_clip_min(residual_clip_min), residual_clip_max(residual_clip_max),
      residual_sum(rows, 0.0), residual_sq_sum(rows, 0.0) {}

  void operator()(std::size_t begin, std::size_t end) {
    for (std::size_t col = begin; col < end; ++col) {
      int ptr = p[col];
      const int ptr_end = p[col + 1];
      const double log_umi_col = log_umi_ptr[col];
      const double common_factor = common_slope ? std::exp(first_slope * log_umi_col) : 0.0;
      for (int row = 0; row < rows; ++row) {
        double y_value = 0.0;
        if (ptr < ptr_end && i[ptr] == row) {
          y_value = x[ptr];
          ++ptr;
        }
        const double mu_original = common_slope ?
          exp_intercept[row] * common_factor :
          exp_intercept[row] * std::exp(slope_ptr[row] * log_umi_col);
        double variance_original = sct_model_var(mu_original, theta_ptr[row]);
        if (variance_original < min_var) {
          variance_original = min_var;
        }
        const double residual = (y_value - mu_original) / std::sqrt(variance_original);
        const double clipped_residual = sct_clip(residual, residual_clip_min, residual_clip_max);
        residual_sum[row] += clipped_residual;
        residual_sq_sum[row] += clipped_residual * clipped_residual;
      }
    }
  }

  void join(const SCTStatsReducer& rhs) {
    for (int row = 0; row < rows; ++row) {
      residual_sum[row] += rhs.residual_sum[row];
      residual_sq_sum[row] += rhs.residual_sq_sum[row];
    }
  }
};

static inline int sct_thread_chunk_count(const int nitems, const int nthreads) {
  return std::max(1, std::min(nitems, nthreads));
}

// [[Rcpp::export(rng = false)]]
List SCTResidualStatsAndCorrected(
  NumericVector x,
  IntegerVector i,
  IntegerVector p,
  int rows,
  int cols,
  NumericVector theta,
  NumericVector intercept,
  NumericVector slope,
  NumericVector log_umi,
  double target_log_umi,
  double min_var,
  double residual_clip_min,
  double residual_clip_max,
  int n_threads = 1,
  bool compute_corrected = true
) {
  NumericVector residual_mean(rows);
  NumericVector residual_variance(rows);

  const double* theta_ptr = REAL(theta);
  const double* intercept_ptr = REAL(intercept);
  const double* slope_ptr = REAL(slope);
  const double* log_umi_ptr = REAL(log_umi);
  const int* i_ptr = INTEGER(i);
  const int* p_ptr = INTEGER(p);
  const double* x_ptr = REAL(x);
  std::vector<double> exp_intercept(rows);
  std::vector<double> target_mu(rows);
  std::vector<double> target_sqrt_var(rows);
  bool common_slope = true;
  const double first_slope = slope_ptr[0];
  for (int row = 0; row < rows; ++row) {
    exp_intercept[row] = std::exp(intercept_ptr[row]);
    if (std::abs(slope_ptr[row] - first_slope) > 1e-12) {
      common_slope = false;
    }
  }
  if (common_slope) {
    const double target_factor = std::exp(first_slope * target_log_umi);
    for (int row = 0; row < rows; ++row) {
      target_mu[row] = exp_intercept[row] * target_factor;
      target_sqrt_var[row] = std::sqrt(sct_model_var(target_mu[row], theta_ptr[row]));
    }
  } else {
    for (int row = 0; row < rows; ++row) {
      target_mu[row] = exp_intercept[row] * std::exp(slope_ptr[row] * target_log_umi);
      target_sqrt_var[row] = std::sqrt(sct_model_var(target_mu[row], theta_ptr[row]));
    }
  }

  // Residual mean/variance: parallel reduction over cells when requested.
  SCTStatsReducer reducer(
    p_ptr, i_ptr, x_ptr, rows, log_umi_ptr, theta_ptr, slope_ptr,
    exp_intercept, common_slope, first_slope,
    min_var, residual_clip_min, residual_clip_max
  );
  if (n_threads > 1 && cols > 1) {
    const int chunks = sct_thread_chunk_count(cols, n_threads);
    std::vector<SCTStatsReducer> reducers;
    reducers.reserve(chunks);
    for (int chunk = 0; chunk < chunks; ++chunk) {
      reducers.emplace_back(
        p_ptr, i_ptr, x_ptr, rows, log_umi_ptr, theta_ptr, slope_ptr,
        exp_intercept, common_slope, first_slope,
        min_var, residual_clip_min, residual_clip_max
      );
    }
    RcppThread::parallelFor(0, chunks, [&](int chunk) {
      const int begin = (cols * chunk) / chunks;
      const int end = (cols * (chunk + 1)) / chunks;
      reducers[chunk](begin, end);
    }, chunks);
    for (int chunk = 0; chunk < chunks; ++chunk) {
      reducer.join(reducers[chunk]);
    }
  } else {
    reducer(0, cols);
  }

  const double cols_d = static_cast<double>(cols);
  const double denom = static_cast<double>(cols - 1);
  for (int row = 0; row < rows; ++row) {
    residual_mean[row] = reducer.residual_sum[row] / cols_d;
    residual_variance[row] = (
      reducer.residual_sq_sum[row] - reducer.residual_sum[row] * reducer.residual_sum[row] / cols_d
    ) / denom;
  }

  List out = List::create(
    _["residual_mean"] = residual_mean,
    _["residual_variance"] = residual_variance
  );

  // Corrected counts are assembled serially to preserve CSC column ordering.
  if (compute_corrected) {
    std::vector<int> corrected_i;
    std::vector<double> corrected_x;
    std::vector<int> corrected_p(cols + 1, 0);
    corrected_i.reserve(x.size());
    corrected_x.reserve(x.size());

    for (int col = 0; col < cols; ++col) {
      corrected_p[col] = static_cast<int>(corrected_i.size());
      int ptr = p_ptr[col];
      const int ptr_end = p_ptr[col + 1];
      const double log_umi_col = log_umi_ptr[col];
      const double common_factor = common_slope ? std::exp(first_slope * log_umi_col) : 0.0;
      for (int row = 0; row < rows; ++row) {
        double y_value = 0.0;
        if (ptr < ptr_end && i_ptr[ptr] == row) {
          y_value = x_ptr[ptr];
          ++ptr;
        }
        const double mu_original = common_slope ?
          exp_intercept[row] * common_factor :
          exp_intercept[row] * std::exp(slope_ptr[row] * log_umi_col);
        double variance_original = sct_model_var(mu_original, theta_ptr[row]);
        if (variance_original < min_var) {
          variance_original = min_var;
        }
        const double residual = (y_value - mu_original) / std::sqrt(variance_original);
        double corrected = target_mu[row] + residual * target_sqrt_var[row];
        corrected = sct_round0(corrected);
        if (corrected < 0.0) {
          corrected = 0.0;
        }
        if (corrected != 0.0) {
          corrected_i.push_back(row);
          corrected_x.push_back(corrected);
        }
      }
    }
    corrected_p[cols] = static_cast<int>(corrected_i.size());
    out["corrected"] = sct_corrected_matrix(corrected_i, corrected_p, corrected_x, rows, cols);
  }

  return out;
}

// Reducer that fills the dense residual matrix (disjoint columns) and, as a
// side reduction, accumulates the per-row sums used for optional centering.
struct SCTResidualMatrixReducer {
  double* out_ptr;
  const int selected;
  const int* p;
  const int* i;
  const double* x;
  const double* log_umi_ptr;
  const std::vector<double>& selected_exp_intercept;
  const std::vector<double>& selected_slope;
  const std::vector<double>& selected_theta;
  const std::vector<double>& selected_min_var;
  const std::vector<int>& row_to_selected;
  const bool common_slope;
  const double first_slope;
  const double clip_min;
  const double clip_max;
  std::vector<double> local_y;
  std::vector<int> local_touched;
  std::vector<double> row_sum;

  SCTResidualMatrixReducer(
    double* out_ptr,
    int selected,
    const int* p,
    const int* i,
    const double* x,
    const double* log_umi_ptr,
    const std::vector<double>& selected_exp_intercept,
    const std::vector<double>& selected_slope,
    const std::vector<double>& selected_theta,
    const std::vector<double>& selected_min_var,
    const std::vector<int>& row_to_selected,
    bool common_slope,
    double first_slope,
    double clip_min,
    double clip_max
  ) : out_ptr(out_ptr), selected(selected), p(p), i(i), x(x),
      log_umi_ptr(log_umi_ptr), selected_exp_intercept(selected_exp_intercept),
      selected_slope(selected_slope), selected_theta(selected_theta),
      selected_min_var(selected_min_var), row_to_selected(row_to_selected),
      common_slope(common_slope), first_slope(first_slope),
      clip_min(clip_min), clip_max(clip_max),
      local_y(selected, 0.0), row_sum(selected, 0.0) {
    local_touched.reserve(selected);
  }

  void operator()(std::size_t begin, std::size_t end) {
    for (std::size_t col = begin; col < end; ++col) {
      local_touched.clear();
      for (int ptr = p[col]; ptr < p[col + 1]; ++ptr) {
        const int selected_row = row_to_selected[i[ptr]];
        if (selected_row >= 0) {
          local_y[selected_row] = x[ptr];
          local_touched.push_back(selected_row);
        }
      }

      double* out_col = out_ptr + static_cast<R_xlen_t>(col) * selected;
      const double log_umi_col = log_umi_ptr[col];
      const double common_factor = common_slope ? std::exp(first_slope * log_umi_col) : 0.0;
      for (int idx = 0; idx < selected; ++idx) {
        const double mu = common_slope ?
          selected_exp_intercept[idx] * common_factor :
          selected_exp_intercept[idx] * std::exp(selected_slope[idx] * log_umi_col);
        double variance = sct_model_var(mu, selected_theta[idx]);
        if (variance < selected_min_var[idx]) {
          variance = selected_min_var[idx];
        }
        const double residual = (local_y[idx] - mu) / std::sqrt(variance);
        const double clipped = sct_clip(residual, clip_min, clip_max);
        out_col[idx] = clipped;
        row_sum[idx] += clipped;
      }

      for (std::vector<int>::const_iterator it = local_touched.begin(); it != local_touched.end(); ++it) {
        local_y[*it] = 0.0;
      }
    }
  }

  void join(const SCTResidualMatrixReducer& rhs) {
    for (int idx = 0; idx < selected; ++idx) {
      row_sum[idx] += rhs.row_sum[idx];
    }
  }
};

// Worker that subtracts each selected gene's residual mean (centering).
struct SCTCenterWorker {
  double* out_ptr;
  const int selected;
  const std::vector<double>& row_sum;
  const double inv_cols;

  SCTCenterWorker(double* out_ptr, int selected, const std::vector<double>& row_sum, double inv_cols)
    : out_ptr(out_ptr), selected(selected), row_sum(row_sum), inv_cols(inv_cols) {}

  void operator()(std::size_t begin, std::size_t end) {
    for (std::size_t col = begin; col < end; ++col) {
      double* out_col = out_ptr + static_cast<R_xlen_t>(col) * selected;
      for (int idx = 0; idx < selected; ++idx) {
        out_col[idx] -= row_sum[idx] * inv_cols;
      }
    }
  }
};

// [[Rcpp::export(rng = false)]]
NumericMatrix SCTPearsonResidualMatrix(
  NumericVector x,
  IntegerVector i,
  IntegerVector p,
  int rows,
  int cols,
  NumericVector theta,
  NumericVector intercept,
  NumericVector slope,
  NumericVector log_umi,
  IntegerVector feature_index,
  NumericVector min_var,
  double clip_min,
  double clip_max,
  bool do_center = true,
  int n_threads = 1
) {
  const int selected = feature_index.size();
  NumericMatrix out = no_init_matrix(selected, cols);
  double* out_ptr = REAL(out);

  std::vector<int> row_to_selected(rows, -1);
  for (int idx = 0; idx < selected; ++idx) {
    row_to_selected[feature_index[idx]] = idx;
  }

  const double* theta_ptr = REAL(theta);
  const double* intercept_ptr = REAL(intercept);
  const double* slope_ptr = REAL(slope);
  const double* log_umi_ptr = REAL(log_umi);
  const int* i_ptr = INTEGER(i);
  const int* p_ptr = INTEGER(p);
  const double* x_ptr = REAL(x);
  std::vector<double> selected_theta(selected);
  std::vector<double> selected_exp_intercept(selected);
  std::vector<double> selected_slope(selected);
  std::vector<double> selected_min_var(selected);
  const bool scalar_min_var = min_var.size() == 1;
  if (!scalar_min_var && min_var.size() != selected) {
    stop("min_var must have length 1 or match feature_index");
  }
  bool common_slope = true;
  const double first_slope = selected > 0 ? slope_ptr[feature_index[0]] : 0.0;
  for (int idx = 0; idx < selected; ++idx) {
    const int row = feature_index[idx];
    selected_theta[idx] = theta_ptr[row];
    selected_exp_intercept[idx] = std::exp(intercept_ptr[row]);
    selected_slope[idx] = slope_ptr[row];
    selected_min_var[idx] = scalar_min_var ? min_var[0] : min_var[idx];
    if (std::abs(selected_slope[idx] - first_slope) > 1e-12) {
      common_slope = false;
    }
  }

  const bool run_parallel = (n_threads > 1 && cols > 1);

  // Fill residuals (disjoint columns) and accumulate per-gene sums for centering.
  SCTResidualMatrixReducer reducer(
    out_ptr, selected, p_ptr, i_ptr, x_ptr, log_umi_ptr,
    selected_exp_intercept, selected_slope, selected_theta, selected_min_var,
    row_to_selected, common_slope, first_slope, clip_min, clip_max
  );
  if (run_parallel) {
    const int chunks = sct_thread_chunk_count(cols, n_threads);
    std::vector<SCTResidualMatrixReducer> reducers;
    reducers.reserve(chunks);
    for (int chunk = 0; chunk < chunks; ++chunk) {
      reducers.emplace_back(
        out_ptr, selected, p_ptr, i_ptr, x_ptr, log_umi_ptr,
        selected_exp_intercept, selected_slope, selected_theta, selected_min_var,
        row_to_selected, common_slope, first_slope, clip_min, clip_max
      );
    }
    RcppThread::parallelFor(0, chunks, [&](int chunk) {
      const int begin = (cols * chunk) / chunks;
      const int end = (cols * (chunk + 1)) / chunks;
      reducers[chunk](begin, end);
    }, chunks);
    for (int chunk = 0; chunk < chunks; ++chunk) {
      reducer.join(reducers[chunk]);
    }
  } else {
    reducer(0, cols);
  }

  if (do_center) {
    SCTCenterWorker center(out_ptr, selected, reducer.row_sum, 1.0 / static_cast<double>(cols));
    if (run_parallel) {
      RcppThread::parallelFor(0, cols, [&](int col) {
        center(col, col + 1);
      }, n_threads);
    } else {
      center(0, cols);
    }
  }

  return out;
}