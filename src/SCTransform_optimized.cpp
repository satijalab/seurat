#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <thread>
#include <vector>

using namespace Rcpp;

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

inline int sct_thread_count(int requested, int cols) {
  if (cols <= 1) {
    return 1;
  }
  int threads = requested;
  if (threads < 1) {
    threads = static_cast<int>(std::thread::hardware_concurrency());
  }
  if (threads < 1) {
    threads = 1;
  }
  return std::min(threads, cols);
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

struct SCTStatsWorkerResult {
  int start = 0;
  int end = 0;
  std::vector<double> residual_sum;
  std::vector<double> residual_sq_sum;
  std::vector<int> corrected_i;
  std::vector<double> corrected_x;
  std::vector<int> corrected_p;
};

// [[Rcpp::export(rng = false)]]
List SCTResidualStatsAndCorrected_optimized(
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

  const int threads = sct_thread_count(n_threads, cols);
  if (threads == 1) {
    std::vector<double> residual_sum(rows, 0.0);
    std::vector<double> residual_sq_sum(rows, 0.0);
    std::vector<int> corrected_i;
    std::vector<double> corrected_x;
    std::vector<int> corrected_p(cols + 1, 0);
    if (compute_corrected) {
      corrected_i.reserve(x.size());
      corrected_x.reserve(x.size());
    }

    for (int col = 0; col < cols; ++col) {
      if (compute_corrected) {
        corrected_p[col] = static_cast<int>(corrected_i.size());
      }
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

        if (compute_corrected) {
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
    }
    if (compute_corrected) {
      corrected_p[cols] = static_cast<int>(corrected_i.size());
    }

    const double cols_d = static_cast<double>(cols);
    const double denom = static_cast<double>(cols - 1);
    for (int row = 0; row < rows; ++row) {
      residual_mean[row] = residual_sum[row] / cols_d;
      residual_variance[row] = (
        residual_sq_sum[row] - residual_sum[row] * residual_sum[row] / cols_d
      ) / denom;
    }

    List out = List::create(
      _["residual_mean"] = residual_mean,
      _["residual_variance"] = residual_variance
    );
    if (compute_corrected) {
      S4 corrected = sct_corrected_matrix(corrected_i, corrected_p, corrected_x, rows, cols);
      out["corrected"] = corrected;
    }
    return out;
  }

  std::vector<SCTStatsWorkerResult> results(threads);
  std::vector<std::thread> workers;
  workers.reserve(threads);
  const int chunk = (cols + threads - 1) / threads;

  for (int thread = 0; thread < threads; ++thread) {
    const int start = thread * chunk;
    const int end = std::min(cols, start + chunk);
    results[thread].start = start;
    results[thread].end = end;
    workers.emplace_back([&, thread, start, end]() {
      SCTStatsWorkerResult& result = results[thread];
      result.residual_sum.assign(rows, 0.0);
      result.residual_sq_sum.assign(rows, 0.0);
      if (compute_corrected) {
        result.corrected_p.assign(end - start + 1, 0);
        const R_xlen_t start_nnz = p[start];
        const R_xlen_t end_nnz = p[end];
        const R_xlen_t reserve_nnz = std::max<R_xlen_t>(end_nnz - start_nnz, 1);
        result.corrected_i.reserve(static_cast<size_t>(reserve_nnz));
        result.corrected_x.reserve(static_cast<size_t>(reserve_nnz));
      }

      for (int col = start; col < end; ++col) {
        if (compute_corrected) {
          result.corrected_p[col - start] = result.corrected_i.size();
        }
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
          result.residual_sum[row] += clipped_residual;
          result.residual_sq_sum[row] += clipped_residual * clipped_residual;

          if (compute_corrected) {
            double corrected = target_mu[row] + residual * target_sqrt_var[row];
            corrected = sct_round0(corrected);
            if (corrected < 0.0) {
              corrected = 0.0;
            }
            if (corrected != 0.0) {
              result.corrected_i.push_back(row);
              result.corrected_x.push_back(corrected);
            }
          }
        }
      }
      if (compute_corrected) {
        result.corrected_p[end - start] = result.corrected_i.size();
      }
    });
  }
  for (std::thread& worker : workers) {
    worker.join();
  }

  std::vector<double> residual_sum(rows, 0.0);
  std::vector<double> residual_sq_sum(rows, 0.0);
  size_t corrected_nnz = 0;
  for (int thread = 0; thread < threads; ++thread) {
    if (compute_corrected) {
      corrected_nnz += results[thread].corrected_i.size();
    }
    for (int row = 0; row < rows; ++row) {
      residual_sum[row] += results[thread].residual_sum[row];
      residual_sq_sum[row] += results[thread].residual_sq_sum[row];
    }
  }

  std::vector<int> corrected_i;
  std::vector<double> corrected_x;
  std::vector<int> corrected_p(cols + 1, 0);
  size_t offset = 0;
  if (compute_corrected) {
    corrected_i.reserve(corrected_nnz);
    corrected_x.reserve(corrected_nnz);
    for (int thread = 0; thread < threads; ++thread) {
      const SCTStatsWorkerResult& result = results[thread];
      const int local_cols = result.end - result.start;
      for (int local_col = 0; local_col < local_cols; ++local_col) {
        corrected_p[result.start + local_col] = static_cast<int>(offset + result.corrected_p[local_col]);
      }
      corrected_i.insert(corrected_i.end(), result.corrected_i.begin(), result.corrected_i.end());
      corrected_x.insert(corrected_x.end(), result.corrected_x.begin(), result.corrected_x.end());
      offset += result.corrected_i.size();
    }
    corrected_p[cols] = static_cast<int>(offset);
  }

  const double cols_d = static_cast<double>(cols);
  const double denom = static_cast<double>(cols - 1);
  for (int row = 0; row < rows; ++row) {
    residual_mean[row] = residual_sum[row] / cols_d;
    residual_variance[row] = (
      residual_sq_sum[row] - residual_sum[row] * residual_sum[row] / cols_d
    ) / denom;
  }

  List out = List::create(
    _["residual_mean"] = residual_mean,
    _["residual_variance"] = residual_variance
  );
  if (compute_corrected) {
    S4 corrected = sct_corrected_matrix(corrected_i, corrected_p, corrected_x, rows, cols);
    out["corrected"] = corrected;
  }
  return out;
}

// [[Rcpp::export(rng = false)]]
NumericMatrix SCTPearsonResidualMatrix_optimized(
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
  const int threads = sct_thread_count(n_threads, cols);

  std::vector<int> row_to_selected(rows, -1);
  for (int idx = 0; idx < selected; ++idx) {
    row_to_selected[feature_index[idx]] = idx;
  }

  const double* theta_ptr = REAL(theta);
  const double* intercept_ptr = REAL(intercept);
  const double* slope_ptr = REAL(slope);
  const double* log_umi_ptr = REAL(log_umi);
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
  std::vector<double> row_sum(selected, 0.0);
  std::vector<std::vector<double> > thread_row_sums(threads, std::vector<double>(selected, 0.0));
  std::vector<std::thread> workers;
  workers.reserve(threads);
  const int chunk = (cols + threads - 1) / threads;

  for (int thread = 0; thread < threads; ++thread) {
    const int start = thread * chunk;
    const int end = std::min(cols, start + chunk);
    workers.emplace_back([&, thread, start, end]() {
      std::vector<double> local_y(selected, 0.0);
      std::vector<int> local_touched;
      local_touched.reserve(selected);
      std::vector<double>& local_row_sum = thread_row_sums[thread];

      for (int col = start; col < end; ++col) {
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
          local_row_sum[idx] += clipped;
        }

        for (std::vector<int>::const_iterator it = local_touched.begin(); it != local_touched.end(); ++it) {
          local_y[*it] = 0.0;
        }
      }
    });
  }
  for (std::thread& worker : workers) {
    worker.join();
  }

  for (int thread = 0; thread < threads; ++thread) {
    for (int idx = 0; idx < selected; ++idx) {
      row_sum[idx] += thread_row_sums[thread][idx];
    }
  }

  if (do_center) {
    workers.clear();
    for (int thread = 0; thread < threads; ++thread) {
      const int start = thread * chunk;
      const int end = std::min(cols, start + chunk);
      workers.emplace_back([&, start, end]() {
        for (int col = start; col < end; ++col) {
          double* out_col = out_ptr + static_cast<R_xlen_t>(col) * selected;
          for (int idx = 0; idx < selected; ++idx) {
            out_col[idx] -= row_sum[idx] / static_cast<double>(cols);
          }
        }
      });
    }
    for (std::thread& worker : workers) {
      worker.join();
    }
  }

  return out;
}
