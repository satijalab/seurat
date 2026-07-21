#include <RcppEigen.h>
#include <RcppParallel.h>
#include <progress.hpp>
#include <algorithm>
#include <cmath>
#include <unordered_map>
#include <fstream>
#include <string>
#include <Rinternals.h>

using namespace Rcpp;
using namespace RcppParallel;

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppProgress)]]
// [[Rcpp::depends(RcppParallel)]]

// [[Rcpp::export]]
Eigen::SparseMatrix<double> RunUMISampling(Eigen::SparseMatrix<double> data, int sample_val, bool upsample = false, bool display_progress=true){
    Progress p(data.outerSize(), display_progress);
    Eigen::VectorXd colSums = data.transpose() * Eigen::VectorXd::Ones(data.rows());
    for (int k=0; k < data.outerSize(); ++k){
      p.increment();
      for (Eigen::SparseMatrix<double>::InnerIterator it(data, k); it; ++it){
        double entry = it.value();
        if( (upsample) || (colSums[k] > sample_val)){
          entry = entry * double(sample_val) / colSums[k];
          if (fmod(entry, 1) != 0){
            double rn = R::runif(0,1);
            if(fmod(entry, 1) <= rn){
              it.valueRef() = floor(entry);
            }
            else{
              it.valueRef() = ceil(entry);
            }
          }
          else{
            it.valueRef() = entry;
          }
        }
      }
    }
  return(data);
}

// [[Rcpp::export]]
Eigen::SparseMatrix<double> RunUMISamplingPerCell(Eigen::SparseMatrix<double> data, NumericVector sample_val, bool upsample = false, bool display_progress=true){
  Progress p(data.outerSize(), display_progress);
  Eigen::VectorXd colSums = data.transpose() * Eigen::VectorXd::Ones(data.rows());
  for (int k=0; k < data.outerSize(); ++k){
    p.increment();
    for (Eigen::SparseMatrix<double>::InnerIterator it(data, k); it; ++it){
      double entry = it.value();
      if( (upsample) || (colSums[k] > sample_val[k])){
        entry = entry * double(sample_val[k]) / colSums[k];
        if (fmod(entry, 1) != 0){
          double rn = R::runif(0,1);
          if(fmod(entry, 1) <= rn){
            it.valueRef() = floor(entry);
          }
          else{
            it.valueRef() = ceil(entry);
          }
        }
        else{
          it.valueRef() = entry;
        }
      }
    }
  }
  return(data);
}


typedef Eigen::Triplet<double> T;
// [[Rcpp::export(rng = false)]]
Eigen::SparseMatrix<double> RowMergeMatrices(Eigen::SparseMatrix<double, Eigen::RowMajor> mat1, Eigen::SparseMatrix<double, Eigen::RowMajor> mat2, std::vector< std::string > mat1_rownames,
                                             std::vector< std::string > mat2_rownames, std::vector< std::string > all_rownames){


  // Set up hash maps for rowname based lookup
  std::unordered_map<std::string, int> mat1_map;
  for(unsigned int i = 0; i < mat1_rownames.size(); i++){
    mat1_map[mat1_rownames[i]] = i;
  }
  std::unordered_map<std::string, int> mat2_map;
  for(unsigned int i = 0; i < mat2_rownames.size(); i++){
    mat2_map[mat2_rownames[i]] = i;
  }

  // set up tripletList for new matrix creation
  std::vector<T> tripletList;
  int num_rows = all_rownames.size();
  int num_col1 = mat1.cols();
  int num_col2 = mat2.cols();


  tripletList.reserve(mat1.nonZeros() + mat2.nonZeros());
  for(int i = 0; i < num_rows; i++){
    std::string key = all_rownames[i];
    if (mat1_map.count(key)){
      for(Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it1(mat1, mat1_map[key]); it1; ++it1){
        tripletList.emplace_back(i, it1.col(), it1.value());
      }
    }
    if (mat2_map.count(key)){
      for(Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it2(mat2, mat2_map[key]); it2; ++it2){
        tripletList.emplace_back(i, num_col1 + it2.col(), it2.value());
      }
    }
  }
  Eigen::SparseMatrix<double> combined_mat(num_rows, num_col1 + num_col2);
  combined_mat.setFromTriplets(tripletList.begin(), tripletList.end());
  return combined_mat;
}

// log-normalize a given column of a sparse matrix
inline void LogNormColumn(const int col, const int *ip, const double *rx, double *ro, const int scale_factor) {
  double col_sum = 0;
  const int col_start = ip[col];
  const int col_end = ip[col + 1];
  for (int j = col_start; j < col_end; j++) {
    col_sum += rx[j];
  }
  // scale factor and column sum are loop-invariant here - compute once for reuse
  const double mult = scale_factor / col_sum;
  for (int j = col_start; j < col_end; j++) {
    ro[j] = log1p(rx[j] * mult);
  }
}

// worker struct for parallel processing
struct LogNormWorker : public RcppParallel::Worker {
  const int *ip;
  const double *rx;
  double *ro;
  const int scale_factor;

  LogNormWorker(const int *ip, const double *rx, double *ro, const int scale_factor)
    : ip(ip), rx(rx), ro(ro), scale_factor(scale_factor) {}

  // each worker is given a range of columns to process
  void operator()(std::size_t begin, std::size_t end) {
    for (std::size_t i = begin; i < end; i++) {
      LogNormColumn(i, ip, rx, ro, scale_factor);
    }
  }
};

inline void LogNormSerial(const int num_cols, const int *ip, const double *rx, double *ro, const int scale_factor, const bool display_progress) {
  Progress prog(num_cols, display_progress);
  // compute col sums and do normalization in one pass
  for(int i = 0; i < num_cols; i++){
    LogNormColumn(i, ip, rx, ro, scale_factor);
    prog.increment(1);
  }
}

inline void LogNormParallel(const int num_cols, const int *ip, const double *rx, double *ro, const int scale_factor, const int nthreads, const bool display_progress) {
  if (display_progress) {
    Progress prog(num_cols, true);
    LogNormWorker worker(ip, rx, ro, scale_factor);
    // show ~100 progress increments
    // use at least 512 columns per block to avoid too many small calls
    const int block_size = std::max(512, num_cols / 100);
    for (int i = 0; i < num_cols; i += block_size) {
      const int end = std::min(i + block_size, num_cols);
      RcppParallel::parallelFor(i, end, worker, 1, nthreads);
      prog.increment(end - i);
    }
  } else {
    LogNormWorker worker(ip, rx, ro, scale_factor);
    RcppParallel::parallelFor(0, num_cols, worker, 1, nthreads);
  }
}

// log normalize that uses the x and p slots from the sparse matrix class
// x consists of the actual non-zero values and p consists of offsets for the start of each matrix column
// [[Rcpp::export(rng = false)]]
NumericVector LogNorm(NumericVector x, IntegerVector p, int scale_factor, int nthreads, bool display_progress) {
  NumericVector out(no_init(x.size()));

  // we use vector accessor functions to get pointers to the underlying data of the Rcpp vectors
  // R0 indicates an accessor that returns pointers as const for read-only access
  const int *ip = INTEGER_RO(p);
  const double *rx = REAL_RO(x);
  double *ro = REAL(out);

  const int num_cols = p.size() - 1;
  if (nthreads > 1) {
    LogNormParallel(num_cols, ip, rx, ro, scale_factor, nthreads, display_progress);
  } else {
    LogNormSerial(num_cols, ip, rx, ro, scale_factor, display_progress);
  }
  return(out);
}

/* Performs column scaling and/or centering. Equivalent to using scale(mat, TRUE, apply(x,2,sd)) in R.
 Note: Doesn't handle NA/NaNs in the same way the R implementation does, */

// [[Rcpp::export(rng = false)]]
NumericMatrix Standardize(Eigen::Map<Eigen::MatrixXd> mat, bool display_progress = true){
  Progress p(mat.cols(), display_progress);
  NumericMatrix std_mat(mat.rows(), mat.cols());
  for(int i=0; i < mat.cols(); ++i){
    p.increment();
    Eigen::ArrayXd r = mat.col(i).array();
    double colMean = r.mean();
    double colSdev = sqrt((r - colMean).square().sum() / (mat.rows() - 1));
    NumericMatrix::Column new_col = std_mat(_, i);
    for(int j=0; j < new_col.size(); j++) {
      new_col[j] = (r[j] - colMean) / colSdev;
    }
  }
  return std_mat;
}

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
NumericMatrix FastSparseRowScale(NumericVector x,
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
NumericMatrix FastDenseRowScale(
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

// [[Rcpp::export(rng = false)]]
Eigen::MatrixXd FastSparseRowScaleWithKnownStats(Eigen::SparseMatrix<double> mat, NumericVector mu, NumericVector sigma, bool scale = true, bool center = true,
                                   double scale_max = 10, bool display_progress = true){
    mat = mat.transpose();
    Progress p(mat.outerSize(), display_progress);
    Eigen::MatrixXd scaled_mat(mat.rows(), mat.cols());
    for (int k=0; k<mat.outerSize(); ++k){
        p.increment();
        double colMean = 0;
        double colSdev = 1;
        if (scale == true){
            colSdev = sigma[k];
        }
        if(center == true){
            colMean = mu[k];
        }
        Eigen::VectorXd col = Eigen::VectorXd(mat.col(k));
        scaled_mat.col(k) = (col.array() - colMean) / colSdev;
        for(int s=0; s<scaled_mat.col(k).size(); ++s){
            if(scaled_mat(s,k) > scale_max){
                scaled_mat(s,k) = scale_max;
            }
        }
    }
    return scaled_mat.transpose();
}

/* Note: May not handle NA/NaNs in the same way the R implementation does, */

// [[Rcpp::export(rng = false)]]
Eigen::MatrixXd FastCov(Eigen::MatrixXd mat, bool center = true){
  if (center) {
    mat = mat.rowwise() - mat.colwise().mean();
  }
  Eigen::MatrixXd cov = (mat.adjoint() * mat) / double(mat.rows() - 1);
  return(cov);
}

// [[Rcpp::export(rng = false)]]
Eigen::MatrixXd FastCovMats(Eigen::MatrixXd mat1, Eigen::MatrixXd mat2, bool center = true){
  if(center){
    mat1 = mat1.rowwise() - mat1.colwise().mean();
    mat2 = mat2.rowwise() - mat2.colwise().mean();
  }
  Eigen::MatrixXd cov = (mat1.adjoint() * mat2) / double(mat1.rows() - 1);
  return(cov);
}

/* Note: Faster than the R implementation but is not in-place */
// [[Rcpp::export(rng = false)]]
Eigen::MatrixXd FastRBind(Eigen::MatrixXd mat1, Eigen::MatrixXd mat2){
  Eigen::MatrixXd mat3(mat1.rows() + mat2.rows(), mat1.cols());
  mat3 << mat1, mat2;
  return(mat3);
}

/* Calculates the row means of the logged values in non-log space */
// [[Rcpp::export(rng = false)]]
Eigen::VectorXd FastExpMean(Eigen::SparseMatrix<double> mat, bool display_progress){
  int ncols = mat.cols();
  Eigen::VectorXd rowmeans(mat.rows());
  mat = mat.transpose();
  if(display_progress == true){
    Rcpp::Rcerr << "Calculating gene means" << std::endl;
  }
  Progress p(mat.outerSize(), display_progress);
  for (int k=0; k<mat.outerSize(); ++k){
    p.increment();
    double rm = 0;
    for (Eigen::SparseMatrix<double>::InnerIterator it(mat,k); it; ++it){
      rm += expm1(it.value());
    }
    rm = rm / ncols;
    rowmeans[k] = log1p(rm);
  }
  return(rowmeans);
}


struct SparseRowMeanVarWorker : public Worker {
  const int* p;
  const int* i;
  const double* x;
  const int rows;
  std::vector<double> sum;
  std::vector<double> sumsq;

  SparseRowMeanVarWorker(
    const int* p,
    const int* i,
    const double* x,
    const int rows
  ) : p(p), i(i), x(x), rows(rows), sum(rows, 0.0), sumsq(rows, 0.0) {}

  SparseRowMeanVarWorker(const SparseRowMeanVarWorker& other, Split)
    : p(other.p), i(other.i), x(other.x), rows(other.rows),
      sum(rows, 0.0), sumsq(rows, 0.0) {}

  void operator()(std::size_t begin, std::size_t end) {
    for (std::size_t col = begin; col < end; ++col) {
      for (int idx = p[col]; idx < p[col + 1]; ++idx) {
        const int row = i[idx];
        const double value = x[idx];
        sum[row] += value;
        sumsq[row] += value * value;
      }
    }
  }

  void join(const SparseRowMeanVarWorker& rhs) {
    for (int row = 0; row < rows; ++row) {
      sum[row] += rhs.sum[row];
      sumsq[row] += rhs.sumsq[row];
    }
  }
};

struct SparseRowVarStdWorker : public Worker {
  const int* p;
  const int* i;
  const double* x;
  const double* mu;
  const double* inv_sd;
  const double* clip_threshold;
  const double vmax;
  const int rows;
  std::vector<double> sumsq;
  std::vector<int> nnz;

  SparseRowVarStdWorker(
    const int* p,
    const int* i,
    const double* x,
    const double* mu,
    const double* inv_sd,
    const double* clip_threshold,
    const double vmax,
    const int rows
  ) : p(p), i(i), x(x), mu(mu), inv_sd(inv_sd), clip_threshold(clip_threshold),
      vmax(vmax), rows(rows), sumsq(rows, 0.0), nnz(rows, 0) {}

  SparseRowVarStdWorker(const SparseRowVarStdWorker& other, Split)
    : p(other.p), i(other.i), x(other.x), mu(other.mu), inv_sd(other.inv_sd),
      clip_threshold(other.clip_threshold), vmax(other.vmax), rows(other.rows),
      sumsq(rows, 0.0), nnz(rows, 0) {}

  void operator()(std::size_t begin, std::size_t end) {
    for (std::size_t col = begin; col < end; ++col) {
      for (int idx = p[col]; idx < p[col + 1]; ++idx) {
        const int row = i[idx];
        if (inv_sd[row] == 0.0) {
          continue;
        }
        double value_sq;
        if (x[idx] > clip_threshold[row]) {
          value_sq = vmax * vmax;
        } else {
          const double value = (x[idx] - mu[row]) * inv_sd[row];
          value_sq = value * value;
        }
        sumsq[row] += value_sq;
        nnz[row] += 1;
      }
    }
  }

  void join(const SparseRowVarStdWorker& rhs) {
    for (int row = 0; row < rows; ++row) {
      sumsq[row] += rhs.sumsq[row];
      nnz[row] += rhs.nnz[row];
    }
  }
};

inline void SparseRowMeanVarParallel(
  NumericVector means,
  NumericVector vars,
  const int* p_ptr,
  const int* i_ptr,
  const double* x_ptr,
  const int rows,
  const int cols,
  const int nthreads,
  const bool display_progress
) {
  const double cols_d = static_cast<double>(cols);
  const double denom = static_cast<double>(cols - 1);
  const int block_size = std::max(512, cols / 100);

  if (display_progress) {
    Progress prog(cols + rows, true);
    for (int begin = 0; begin < cols; begin += block_size) {
      const int end = std::min(begin + block_size, cols);
      SparseRowMeanVarWorker worker(p_ptr, i_ptr, x_ptr, rows);
      parallelReduce(begin, end, worker, 1, nthreads);
      for (int row = 0; row < rows; ++row) {
        REAL(means)[row] += worker.sum[row];
        REAL(vars)[row] += worker.sumsq[row];
      }
      prog.increment(end - begin);
    }
    for (int row = 0; row < rows; ++row) {
      prog.increment();
      const double sum = REAL(means)[row];
      REAL(means)[row] = sum / cols_d;
      REAL(vars)[row] = (REAL(vars)[row] - (sum * sum / cols_d)) / denom;
    }
  } else {
    SparseRowMeanVarWorker worker(p_ptr, i_ptr, x_ptr, rows);
    parallelReduce(0, cols, worker, 1, nthreads);
    for (int row = 0; row < rows; ++row) {
      const double sum = worker.sum[row];
      REAL(means)[row] = sum / cols_d;
      REAL(vars)[row] = (worker.sumsq[row] - (sum * sum / cols_d)) / denom;
    }
  }
}

inline void SparseRowVarStdParallel(
  NumericVector vars,
  IntegerVector nnz,
  const int* p_ptr,
  const int* i_ptr,
  const double* x_ptr,
  const double* mu_ptr,
  const double* inv_sd_ptr,
  const double* zero_value_sq_ptr,
  const double* clip_threshold_ptr,
  const int rows,
  const int cols,
  const double vmax,
  const int nthreads,
  const bool display_progress
) {
  const double denom = static_cast<double>(cols - 1);
  const int block_size = std::max(512, cols / 100);

  if (display_progress) {
    Progress prog(rows + cols + rows, true);
    for (int row = 0; row < rows; ++row) {
      prog.increment();
    }
    for (int begin = 0; begin < cols; begin += block_size) {
      const int end = std::min(begin + block_size, cols);
      SparseRowVarStdWorker worker(
        p_ptr, i_ptr, x_ptr, mu_ptr, inv_sd_ptr, clip_threshold_ptr, vmax, rows
      );
      parallelReduce(begin, end, worker, 1, nthreads);
      for (int row = 0; row < rows; ++row) {
        REAL(vars)[row] += worker.sumsq[row];
        INTEGER(nnz)[row] += worker.nnz[row];
      }
      prog.increment(end - begin);
    }
    for (int row = 0; row < rows; ++row) {
      prog.increment();
      if (inv_sd_ptr[row] == 0.0) {
        continue;
      }
      const int nzero = cols - INTEGER(nnz)[row];
      REAL(vars)[row] = (REAL(vars)[row] + (zero_value_sq_ptr[row] * nzero)) / denom;
    }
  } else {
    SparseRowVarStdWorker worker(
      p_ptr, i_ptr, x_ptr, mu_ptr, inv_sd_ptr, clip_threshold_ptr, vmax, rows
    );
    parallelReduce(0, cols, worker, 1, nthreads);
    for (int row = 0; row < rows; ++row) {
      if (inv_sd_ptr[row] == 0.0) {
        continue;
      }
      const int nzero = cols - worker.nnz[row];
      REAL(vars)[row] = (worker.sumsq[row] + (zero_value_sq_ptr[row] * nzero)) / denom;
    }
  }
}

// [[Rcpp::export(rng = false)]]
List SparseRowMeanVar(
  NumericVector x,
  IntegerVector i,
  IntegerVector p,
  int rows,
  int cols,
  int nthreads,
  bool display_progress
) {
  NumericVector means(rows);
  NumericVector vars(rows);
  const double* x_ptr = REAL(x);
  const int* i_ptr = INTEGER(i);
  const int* p_ptr = INTEGER(p);

  if (nthreads > 1) {
    SparseRowMeanVarParallel(
      means, vars, p_ptr, i_ptr, x_ptr, rows, cols, nthreads, display_progress
    );
  } else {
    double* means_ptr = REAL(means);
    double* vars_ptr = REAL(vars);
    const R_xlen_t x_size = x.size();
    const double cols_d = static_cast<double>(cols);
    const double denom = static_cast<double>(cols - 1);
    Progress prog(x_size + rows, display_progress);

    for (R_xlen_t idx = 0; idx < x_size; ++idx) {
      prog.increment();
      const int row = i_ptr[idx];
      const double value = x_ptr[idx];
      means_ptr[row] += value;
      vars_ptr[row] += value * value;
    }
    for (int row = 0; row < rows; ++row) {
      prog.increment();
      const double sum = means_ptr[row];
      means_ptr[row] = sum / cols_d;
      vars_ptr[row] = (vars_ptr[row] - (sum * sum / cols_d)) / denom;
    }
  }
  return List::create(
    _["mean"] = means,
    _["variance"] = vars
  );
}

/* standardize matrix rows using given mean and standard deviation,
   clip values larger than vmax to vmax,
   then return variance for each row */
// [[Rcpp::export(rng = false)]]
NumericVector SparseRowVarStd(
  NumericVector x,
  IntegerVector i,
  IntegerVector p,
  NumericVector mu,
  NumericVector sd,
  int rows,
  int cols,
  double vmax,
  int nthreads,
  bool display_progress
) {
  NumericVector vars(rows);
  std::vector<double> inv_sd(rows, 0.0);
  std::vector<double> zero_value_sq(rows, 0.0);
  std::vector<double> clip_threshold(rows, 0.0);
  const double* x_ptr = REAL(x);
  const int* i_ptr = INTEGER(i);
  const int* p_ptr = INTEGER(p);
  const double* mu_ptr = REAL(mu);
  const double* sd_ptr = REAL(sd);
  double* vars_ptr = REAL(vars);

  for (int row = 0; row < rows; ++row) {
    if (sd_ptr[row] == 0.0 || R_IsNA(sd_ptr[row])) {
      continue;
    }
    inv_sd[row] = 1.0 / sd_ptr[row];
    const double zero_value = -mu_ptr[row] * inv_sd[row];
    zero_value_sq[row] = zero_value * zero_value;
    clip_threshold[row] = mu_ptr[row] + (vmax * sd_ptr[row]);
  }

  if (nthreads > 1) {
    IntegerVector nnz(rows);
    SparseRowVarStdParallel(
      vars,
      nnz,
      p_ptr,
      i_ptr,
      x_ptr,
      mu_ptr,
      inv_sd.data(),
      zero_value_sq.data(),
      clip_threshold.data(),
      rows,
      cols,
      vmax,
      nthreads,
      display_progress
    );
  } else {
    const R_xlen_t x_size = x.size();
    const double denom = static_cast<double>(cols - 1);
    Progress prog(rows + x_size + rows, display_progress);

    for (int row = 0; row < rows; ++row) {
      prog.increment();
      vars_ptr[row] = zero_value_sq[row] * static_cast<double>(cols);
    }

    for (R_xlen_t idx = 0; idx < x_size; ++idx) {
      prog.increment();
      const int row = i_ptr[idx];
      if (inv_sd[row] == 0.0) {
        continue;
      }
      double value_sq;
      if (x_ptr[idx] > clip_threshold[row]) {
        value_sq = vmax * vmax;
      } else {
        const double value = (x_ptr[idx] - mu_ptr[row]) * inv_sd[row];
        value_sq = value * value;
      }
      vars_ptr[row] += value_sq - zero_value_sq[row];
    }
    for (int row = 0; row < rows; ++row) {
      prog.increment();
      if (inv_sd[row] == 0.0) {
        continue;
      }
      vars_ptr[row] /= denom;
    }
  }
  return vars;
}

/* Calculate the variance to mean ratio (VMR) in non-logspace (return answer in
log-space) */
// [[Rcpp::export(rng = false)]]
Eigen::VectorXd FastLogVMR(Eigen::SparseMatrix<double> mat,  bool display_progress){
  int ncols = mat.cols();
  Eigen::VectorXd rowdisp(mat.rows());
  mat = mat.transpose();
  if(display_progress == true){
    Rcpp::Rcerr << "Calculating gene variance to mean ratios" << std::endl;
  }
  Progress p(mat.outerSize(), display_progress);
  for (int k=0; k<mat.outerSize(); ++k){
    p.increment();
    double rm = 0;
    double v = 0;
    int nnZero = 0;
    for (Eigen::SparseMatrix<double>::InnerIterator it(mat,k); it; ++it){
      rm += expm1(it.value());
    }
    rm = rm / ncols;
    for (Eigen::SparseMatrix<double>::InnerIterator it(mat,k); it; ++it){
      v += pow(expm1(it.value()) - rm, 2);
      nnZero += 1;
    }
    v = (v + (ncols - nnZero) * pow(rm, 2)) / (ncols - 1);
    rowdisp[k] = log(v/rm);

  }
  return(rowdisp);
}

/* Calculates the variance of rows of a matrix */
// [[Rcpp::export(rng = false)]]
NumericVector RowVar(Eigen::Map<Eigen::MatrixXd> x){
  NumericVector out(x.rows());
  for(int i=0; i < x.rows(); ++i){
    Eigen::ArrayXd r = x.row(i).array();
    double rowMean = r.mean();
    out[i] = (r - rowMean).square().sum() / (x.cols() - 1);
  }
  return out;
}

/* Calculate the variance in non-logspace (return answer in non-logspace) */
// [[Rcpp::export(rng = false)]]
Eigen::VectorXd SparseRowVar(Eigen::SparseMatrix<double> mat, bool display_progress){
  int ncols = mat.cols();
  Eigen::VectorXd rowdisp(mat.rows());
  mat = mat.transpose();
  if(display_progress == true){
    Rcpp::Rcerr << "Calculating gene variances" << std::endl;
  }
  Progress p(mat.outerSize(), display_progress);
  for (int k=0; k<mat.outerSize(); ++k){
    p.increment();
    double rm = 0;
    double v = 0;
    int nnZero = 0;
    for (Eigen::SparseMatrix<double>::InnerIterator it(mat,k); it; ++it){
      rm += (it.value());
    }
    rm = rm / ncols;
    for (Eigen::SparseMatrix<double>::InnerIterator it(mat,k); it; ++it){
      v += pow((it.value()) - rm, 2);
      nnZero += 1;
    }
    v = (v + (ncols - nnZero) * pow(rm, 2)) / (ncols - 1);
    rowdisp[k] = v;
  }
  return(rowdisp);
}

//cols_idx should be 0-indexed
// [[Rcpp::export(rng = false)]]
Eigen::SparseMatrix<double> ReplaceColsC(Eigen::SparseMatrix<double> mat, NumericVector col_idx, Eigen::SparseMatrix<double> replacement){
  int rep_idx = 0;
  for(auto const &ci : col_idx){
    mat.col(ci) = replacement.col(rep_idx);
    rep_idx += 1;
  }
  return(mat);
}

template <typename S>
std::vector<size_t> sort_indexes(const std::vector<S> &v) {
  // initialize original index locations
  std::vector<size_t> idx(v.size());
  std::iota(idx.begin(), idx.end(), 0);
  std::stable_sort(idx.begin(), idx.end(),
                   [&v](size_t i1, size_t i2) {return v[i1] < v[i2];});
  return idx;
}

// [[Rcpp::export(rng = false)]]
List GraphToNeighborHelper(Eigen::SparseMatrix<double> mat) {
  mat = mat.transpose();
  //determine the number of neighbors
  int n = 0;
  for(Eigen::SparseMatrix<double>::InnerIterator it(mat, 0); it; ++it) {
    n += 1;
  }
  Eigen::MatrixXd nn_idx(mat.rows(), n);
  Eigen::MatrixXd nn_dist(mat.rows(), n);

  for (int k=0; k<mat.outerSize(); ++k){
    int n_k = 0;
    std::vector<double> row_idx;
    std::vector<double> row_dist;
    row_idx.reserve(n);
    row_dist.reserve(n);
    for (Eigen::SparseMatrix<double>::InnerIterator it(mat,k); it; ++it) {
      if (n_k > (n-1)) {
        Rcpp::stop("Not all cells have an equal number of neighbors.");
      }
      row_idx.push_back(it.row() + 1);
      row_dist.push_back(it.value());
      n_k += 1;
    }
    if (n_k != n) {
      Rcpp::Rcout << n << ":::" << n_k << std::endl;
      Rcpp::stop("Not all cells have an equal number of neighbors.");
    }
    //order the idx based on dist
    std::vector<size_t> idx_order = sort_indexes(row_dist);
    for(int i = 0; i < n; ++i) {
      nn_idx(k, i) = row_idx[idx_order[i]];
      nn_dist(k, i) = row_dist[idx_order[i]];
    }
  }
  List neighbors = List::create(nn_idx, nn_dist);
  return(neighbors);
}
