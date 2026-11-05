// Vendored BPCells sources used by the scSLIDE on-disk PLS code (pls.h,
// pls.cpp, pls_bridge.cpp).
//
// The PLS bridges receive a MatrixLoader built by the installed BPCells R
// package (via BPCells' iterate_matrix() XPtr) and make virtual calls on it
// through the class declarations below. Those declarations must therefore
// match the installed BPCells' MatrixLoader layout exactly: if a BPCells
// release adds, removes or reorders a virtual method, the calls here would
// silently dispatch to the wrong method instead of failing at build time.
// Re-check this file against BPCells whenever BPCells is upgraded, and do not
// trim unused members from the classes. Built and tested against BPCells 0.3.1
// (commit c293a0a3).
//
// Each section below is copied from the BPCells file named in its header,
// with only the #pragma once / local #include lines removed and the
// StatsResult getters (MatrixStats.cpp) made inline so this stays
// header-only.
//
// Copyright 2021-2023 BPCells contributors
//
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

#pragma once

#include <atomic>
#include <chrono>
#include <cstdint>
#include <future>
#include <memory>
#include <stdexcept>
#include <vector>

#define RCPP_NO_RTTI
#define RCPP_NO_SUGAR
#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]

// ---- from BPCells MatrixStats.h ----

namespace BPCells {

// Enum to request stats to calculate in MatrixLoader::computeMatrixStats
enum class Stats { None = 0, NonZeroCount = 1, Mean = 2, Variance = 3 };

// Result class for row + column stats.
// Each stat is a row in the matrix
class StatsResult {
  public:
    Eigen::ArrayXXd row_stats;
    Eigen::ArrayXXd col_stats;

    Eigen::ArrayXXd rowNonzeros();
    Eigen::ArrayXXd rowMean();
    Eigen::ArrayXXd rowVariance();

    Eigen::ArrayXXd colNonzeros();
    Eigen::ArrayXXd colMean();
    Eigen::ArrayXXd colVariance();

    StatsResult transpose();
};

} // end namespace BPCells

// ---- from BPCells MatrixStats.cpp ----

namespace BPCells {

inline Eigen::ArrayXXd StatsResult::rowNonzeros() {
    if (row_stats.rows() < 1)
        throw std::runtime_error("Nonzero not calculated in this StatsResult");
    return row_stats.row(0);
}
inline Eigen::ArrayXXd StatsResult::rowMean() {
    if (row_stats.rows() < 2) throw std::runtime_error("Mean not calculated in this StatsResult");
    return row_stats.row(1);
}
inline Eigen::ArrayXXd StatsResult::rowVariance() {
    if (row_stats.rows() < 3)
        throw std::runtime_error("Variance not calculated in this StatsResult");
    return row_stats.row(2);
}

inline Eigen::ArrayXXd StatsResult::colNonzeros() {
    if (col_stats.rows() < 1)
        throw std::runtime_error("Nonzero not calculated in this StatsResult");
    return col_stats.row(0);
}
inline Eigen::ArrayXXd StatsResult::colMean() {
    if (col_stats.rows() < 2) throw std::runtime_error("Mean not calculated in this StatsResult");
    return col_stats.row(1);
}
inline Eigen::ArrayXXd StatsResult::colVariance() {
    if (col_stats.rows() < 3)
        throw std::runtime_error("Variance not calculated in this StatsResult");
    return col_stats.row(2);
}

inline StatsResult StatsResult::transpose() { return StatsResult{col_stats, row_stats}; }

} // end namespace BPCells

// ---- from BPCells MatrixIterator.h ----

namespace BPCells {

// Base class to load entries from (sparse) matrices
// Subclasses implement reading from various in-memory and disk formats,
// calculation of matrices from fragments, as well as transformations like
// normalization, and re-indexing.
//
// For simplicity, matrix loaders are purely in column-major order, similar to
// Eigen and R's default compressed sparse column layouts.
// Transposition can be performed only with a MatrixWriter,
// since the re-ordering requires storing all entries in an intermediate matrix.
//
// For flexibility with indexing, matrix iterators need only be grouped by column,
// and neither column or row ordering matters, so long as all entries for a single column
// are consecutive
// To implement a new MatrixIterator:
// 1. Implement rows() and cols() methods to return dimension
// 2. Implement load() method to load the next chunk from the current column.
//    this should return 0 repeatedly at the end of a column until nextCol is called.
// 3. Implement nextCol() method to advance to the next available column.
// 4. Implement restart() method to restart the iterator from the beginning
template <typename T> class MatrixLoader {
  public:
    virtual ~MatrixLoader(){};

    // Return the count of rows and columns
    virtual uint32_t rows() const = 0;
    virtual uint32_t cols() const = 0;

    // Return name for a given row or column.
    // If a matrix doesn't have assigned names this will return NULL
    virtual const char *rowNames(uint32_t row) = 0;
    virtual const char *colNames(uint32_t col) = 0;

    // Reset the iterator to start from the beginning
    virtual void restart() = 0;

    // Seek to a specific column without reading data
    // Next call should be to load(). col must be < cols()
    virtual void seekCol(uint32_t col) = 0;

    // Advance to the next column, or return false if there
    // are no more columns
    virtual bool nextCol() = 0;

    // Return the index of the current column
    virtual uint32_t currentCol() const = 0;

    // Return false if there are no more entries to load
    virtual bool load() = 0;

    // Number of loaded entries available
    virtual uint32_t capacity() const = 0;

    // Pointers to the loaded entries
    virtual uint32_t *rowData() = 0;
    virtual T *valData() = 0;

    // Matrix math operations (implemented in MatrixOps.cpp and MatrixStats.cpp)
    // These operations can be overloaded by matrix transform operations

    // Calculate matrix-matrix product A*B where A (this) is sparse and B is a dense matrix.
    virtual Eigen::MatrixXd denseMultiplyRight(
        const Eigen::Map<Eigen::MatrixXd> B, std::atomic<bool> *user_interrupt = NULL
    );
    virtual Eigen::MatrixXd denseMultiplyLeft(
        const Eigen::Map<Eigen::MatrixXd> B, std::atomic<bool> *user_interrupt = NULL
    );
    // Calculate matrix-vector product A*v where A (this) is sparse and B is a dense matrix.
    virtual Eigen::VectorXd
    vecMultiplyRight(const Eigen::Map<Eigen::VectorXd> v, std::atomic<bool> *user_interrupt = NULL);
    virtual Eigen::VectorXd
    vecMultiplyLeft(const Eigen::Map<Eigen::VectorXd> v, std::atomic<bool> *user_interrupt = NULL);

    // Calculate row/column sums of the matrix
    virtual std::vector<T> colSums(std::atomic<bool> *user_interrupt = NULL);
    virtual std::vector<T> rowSums(std::atomic<bool> *user_interrupt = NULL);

    // Calculate stats on the rows or columns of a matrix in a single pass.
    // For each of rows and columns, the user can choose to from the following
    // Available statistics:
    // 1. None
    // 2. NonZeroCount
    // 3. Mean
    // 4. Variance
    // Each of the later statistics is always calculated simultaneously to the earlier statistics
    // Outputs results to matrices row_output and col_output, which are col-major matrices,
    // with one column per # rows or # columns as appropriate, and one row per output statistic
    virtual StatsResult
    computeMatrixStats(Stats row_stats, Stats col_stats, std::atomic<bool> *user_interrupt = NULL);
};

// Usually downstream users will want MatrixLoaderWrapper.
// This base type is provided in rare cases where we want a wrapper that might
// not match the output type of its input type (e.g. colwise rank)
template <typename Tin, typename Tout = Tin> class MatrixConverterLoaderWrapper : public MatrixLoader<Tout> {
  protected:
    std::unique_ptr<MatrixLoader<Tin>> loader;
    uint32_t current_col = UINT32_MAX - 1;
    bool take_ownership = true;

  public:
    MatrixConverterLoaderWrapper(std::unique_ptr<MatrixLoader<Tin>> &&loader)
        : loader(std::move(loader)){};

    ~MatrixConverterLoaderWrapper() {
        if (!take_ownership) loader.release();
    }
    MatrixConverterLoaderWrapper() = default;
    MatrixConverterLoaderWrapper(MatrixConverterLoaderWrapper &&) = default;
    MatrixConverterLoaderWrapper &operator=(MatrixConverterLoaderWrapper &&) = default;

    // Set the object so that the inner loader will be preserved
    // rather than calling the destructor when this loader is destructed
    void preserve_input_loader() { take_ownership = false; }

    uint32_t rows() const override { return loader->rows(); }
    uint32_t cols() const override { return loader->cols(); }

    const char *rowNames(uint32_t row) override { return loader->rowNames(row); }
    const char *colNames(uint32_t col) override { return loader->colNames(col); }

    void restart() override { loader->restart(); }
    void seekCol(uint32_t col) override { loader->seekCol(col); }

    bool nextCol() override {
        if (!loader->nextCol()) return false;
        current_col = loader->currentCol();
        return true;
    }

    uint32_t currentCol() const override { return loader->currentCol(); }

    bool load() override { return loader->load(); }

    uint32_t capacity() const override { return loader->capacity(); }

    uint32_t *rowData() override { return loader->rowData(); }
};

template <typename T> class MatrixLoaderWrapper : public MatrixConverterLoaderWrapper<T, T> {

  public:
    MatrixLoaderWrapper(std::unique_ptr<MatrixLoader<T>> &&loader)
        : MatrixConverterLoaderWrapper<T,T>(std::move(loader)){};

    MatrixLoaderWrapper() = default;
    MatrixLoaderWrapper(MatrixLoaderWrapper &&) = default;
    MatrixLoaderWrapper &operator=(MatrixLoaderWrapper &&) = default;
    T *valData() override { return this->loader->valData(); }
};

template <typename T> class MatrixIterator : public MatrixLoaderWrapper<T> {
  private:
    uint32_t idx = UINT32_MAX;
    uint32_t current_col;
    uint32_t current_capacity = 0;
    uint32_t *current_row;
    T *current_val;
    bool take_ownership;

  public:
    MatrixIterator(std::unique_ptr<MatrixLoader<T>> &&loader)
        : MatrixLoaderWrapper<T>(std::move(loader)) {}

    // Reset the iterator to start from the beginning
    void restart() override {
        this->loader->restart();
        idx = UINT32_MAX;
        current_capacity = 0;
    }

    void seekCol(uint32_t col) override {
        this->loader->seekCol(col);
        idx = UINT32_MAX;
        current_capacity = 0;
        current_col = this->loader->currentCol();
    }

    // Return false if there isn't another column to access
    bool nextCol() override {
        bool res = this->loader->nextCol();
        if (res) current_col = this->loader->currentCol();
        idx = UINT32_MAX;
        current_capacity = 0;
        return res;
    }

    // Return false if there isn't another entry in the current column
    inline bool nextValue() {
        idx += 1;
        if (idx >= current_capacity) {
            return load();
        }
        return true;
    }
    // Access current row, column, and value
    inline uint32_t row() const { return current_row[idx]; };
    inline uint32_t col() const { return current_col; };
    inline T val() const { return current_val[idx]; };

    bool load() override {
        if (!this->loader->load()) {
            current_capacity = 0;
            return false;
        }
        idx = 0;
        current_capacity = this->loader->capacity();
        current_row = this->loader->rowData();
        current_val = this->loader->valData();
        return true;
    };
    uint32_t capacity() const override { return current_capacity; }
};

template <typename T> class MatrixWriter {
  public:
    virtual ~MatrixWriter(){};
    virtual void write(MatrixLoader<T> &mat, std::atomic<bool> *user_interrupt = NULL) = 0;
};

template <typename Tin, typename Tout> class MatrixConverterLoader : public MatrixLoader<Tout> {
  private:
    std::unique_ptr<MatrixLoader<Tin>> loader;
    std::vector<Tout> vals;

  public:
    MatrixConverterLoader(std::unique_ptr<MatrixLoader<Tin>> &&loader)
        : loader(std::move(loader)) {}

    uint32_t rows() const override { return loader->rows(); }
    uint32_t cols() const override { return loader->cols(); }

    const char *rowNames(uint32_t row) override { return loader->rowNames(row); }
    const char *colNames(uint32_t col) override { return loader->colNames(col); }

    void restart() override { loader->restart(); }
    void seekCol(uint32_t col) override { loader->seekCol(col); }

    bool nextCol() override { return loader->nextCol(); }
    uint32_t currentCol() const override { return loader->currentCol(); }

    bool load() override {
        if (!loader->load()) return false;

        uint32_t capacity = loader->capacity();
        vals.resize(std::max(vals.size(), (size_t)capacity));

        Tin *val_in = loader->valData();

        for (uint32_t i = 0; i < capacity; i++) {
            vals[i] = val_in[i];
        }
        return true;
    };

    uint32_t capacity() const override { return loader->capacity(); }
    uint32_t *rowData() override { return loader->rowData(); }
    Tout *valData() override { return vals.data(); }
};

} // end namespace BPCells

// ---- from BPCells MatrixOps.h ----

namespace BPCells {
// Calculate matrix-matrix product A*B where A (this) is sparse and B is a dense matrix.
template <typename T>
Eigen::MatrixXd MatrixLoader<T>::denseMultiplyRight(
    const Eigen::Map<Eigen::MatrixXd> B, std::atomic<bool> *user_interrupt
) {
    // Use transposed output so that results write contiguously in memory.
    // Note that left multiply will have better performance properties due
    // to better cache locality on the input matrix.
    // So if the performance is key and the size of the input is small,
    // it may be worth paying the price of additional transpose operations
    if (cols() != B.rows()) throw std::runtime_error("Incompatible dimensions for matrix multiply");
    Eigen::MatrixXd res(B.cols(), rows());
    Eigen::VectorXd B_row;
    res.setZero();
    restart();
    while (nextCol()) {
        const uint32_t col = currentCol();
        B_row = B.row(col);
        if (user_interrupt != NULL && *user_interrupt) return res;
        while (load()) {
            const T *val_data = valData();
            const uint32_t *row_data = rowData();
            const uint32_t count = capacity();
            for (uint32_t i = 0; i < count; i++) {
                res.col(row_data[i]) += ((double)val_data[i]) * B_row;
            }
        }
    }
    return res.transpose();
}
template <typename T>
Eigen::MatrixXd MatrixLoader<T>::denseMultiplyLeft(
    const Eigen::Map<Eigen::MatrixXd> B, std::atomic<bool> *user_interrupt
) {
    if (rows() != B.cols()) throw std::runtime_error("Incompatible dimensions for matrix multiply");
    Eigen::MatrixXd res(B.rows(), cols());
    res.setZero();
    restart();
    while (nextCol()) {
        const uint32_t col = currentCol();
        if (user_interrupt != NULL && *user_interrupt) return res;
        while (load()) {
            const T *val_data = valData();
            const uint32_t *row_data = rowData();
            const uint32_t count = capacity();
            for (uint32_t i = 0; i < count; i++) {
                res.col(col) += ((double)val_data[i]) * B.col(row_data[i]);
            }
        }
    }
    return res;
}

// Calculate matrix-vector product A*v where A (this) is sparse and B is a dense matrix.
template <typename T>
Eigen::VectorXd MatrixLoader<T>::vecMultiplyRight(
    const Eigen::Map<Eigen::VectorXd> v, std::atomic<bool> *user_interrupt
) {
    if (cols() != v.rows()) throw std::runtime_error("Incompatible dimensions for vector multiply");
    Eigen::VectorXd res(rows());
    res.setZero();
    restart();
    while (nextCol()) {
        const uint32_t col = currentCol();
        double v_col = v(col);
        if (user_interrupt != NULL && *user_interrupt) return res;
        while (load()) {
            const T *val_data = valData();
            const uint32_t *row_data = rowData();
            const uint32_t count = capacity();
            // I considered manual loop unrolling here too, but it was <10% performance gain for the
            // function so I discarded that idea
            for (uint32_t i = 0; i < count; i++) {
                res(row_data[i]) += ((double)val_data[i]) * v_col;
            }
        }
    }
    return res;
}
template <typename T>
Eigen::VectorXd MatrixLoader<T>::vecMultiplyLeft(
    const Eigen::Map<Eigen::VectorXd> v, std::atomic<bool> *user_interrupt
) {
    if (rows() != v.rows()) throw std::runtime_error("Incompatible dimensions for vector multiply");
    Eigen::VectorXd res(cols());
    res.setZero();
    restart();
    while (nextCol()) {
        const uint32_t col = currentCol();
        if (user_interrupt != NULL && *user_interrupt) return res;
        while (load()) {
            const T *val_data = valData();
            const uint32_t *row_data = rowData();
            const uint32_t count = capacity();
            // Use unrolling trick from Eigen:
            // https://gitlab.com/libeigen/eigen/-/blob/master/Eigen/src/SparseCore/SparseDenseProduct.h#L66-80
            // I tried 4x unrolling and the benefit was marginal so here we are
            uint32_t i;
            double tmp1 = 0;
            double tmp2 = 0;
            for (i = 0; i + 2 <= count; i += 2) {
                tmp1 += ((double)val_data[i]) * v(row_data[i]);
                tmp2 += ((double)val_data[i + 1]) * v(row_data[i + 1]);
            }
            for (; i < count; i++) {
                tmp1 += ((double)val_data[i]) * v(row_data[i]);
            }
            res(col) += tmp1 + tmp2;
        }
    }
    return res;
}

// Calculate row/column sums of the matrix
template <typename T> std::vector<T> MatrixLoader<T>::colSums(std::atomic<bool> *user_interrupt) {
    std::vector<T> sums(cols());
    restart();
    while (nextCol()) {
        const uint32_t col = currentCol();
        if (user_interrupt != NULL && *user_interrupt) return sums;
        while (load()) {
            const T *val_data = valData();
            const uint32_t count = capacity();
            for (uint32_t i = 0; i < count; i++) {
                sums[col] += val_data[i];
            }
        }
    }
    return sums;
}
template <typename T> std::vector<T> MatrixLoader<T>::rowSums(std::atomic<bool> *user_interrupt) {
    std::vector<T> sums(rows());
    restart();
    while (nextCol()) {
        if (user_interrupt != NULL && *user_interrupt) return sums;
        while (load()) {
            const uint32_t *row_data = rowData();
            const T *val_data = valData();
            const uint32_t count = capacity();
            for (uint32_t i = 0; i < count; i++) {
                sums[row_data[i]] += val_data[i];
            }
        }
    }
    return sums;
}

// Helpers for calculating running mean + running variance
inline void update_mean(double val, double &count, double &running_mean) {
    count += 1;
    running_mean += (val - running_mean) / count;
}

inline void update_mean_and_variance(
    double val, double &count, double &running_mean, double &running_variance
) {
    count += 1;
    double delta = val - running_mean;
    running_mean += delta / count;
    double delta2 = val - running_mean;
    running_variance += delta * delta2;
}

// Calculate stats on the rows or columns of a matrix in a single pass.
// For each of rows and columns, the user can choose to from the following
// Available statistics:
// 1. None
// 2. NonZeroCount
// 3. Mean
// 4. Variance
// Each of the later statistics is always calculated simultaneously to the earlier statistics
// Outputs results to matrices row_output and col_output, which are col-major matrices,
// with one column per # rows or # columns as appropriate, and one row per output statistic
template <typename T>
StatsResult MatrixLoader<T>::computeMatrixStats(
    Stats row_stats, Stats col_stats, std::atomic<bool> *user_interrupt
) {
    restart();
    StatsResult res{
        Eigen::ArrayXXd((int)row_stats, rows()), Eigen::ArrayXXd((int)col_stats, cols())};
    res.row_stats.setZero();
    res.col_stats.setZero();

    while (nextCol()) {
        const uint32_t col = currentCol();
        if (user_interrupt != NULL && *user_interrupt) return res;
        while (load()) {
            const uint32_t *row_data = rowData();
            const T *val_data = valData();
            const uint32_t count = capacity();

            // Update row stats if needed
            if (row_stats == Stats::NonZeroCount) {
                for (size_t i = 0; i < count; i++) {
                    res.row_stats(0, row_data[i]) += 1;
                }
            } else if (row_stats == Stats::Mean) {
                for (size_t i = 0; i < count; i++) {
                    update_mean(
                        val_data[i], res.row_stats(0, row_data[i]), res.row_stats(1, row_data[i])
                    );
                }
            } else if (row_stats == Stats::Variance) {
                for (size_t i = 0; i < count; i++) {
                    update_mean_and_variance(
                        val_data[i],
                        res.row_stats(0, row_data[i]),
                        res.row_stats(1, row_data[i]),
                        res.row_stats(2, row_data[i])
                    );
                }
            }
            // Update col stats if needed
            if (col_stats == Stats::NonZeroCount) {
                for (size_t i = 0; i < count; i++) {
                    res.col_stats(0, col) += 1;
                }
            } else if (col_stats == Stats::Mean) {
                for (size_t i = 0; i < count; i++) {
                    update_mean(val_data[i], res.col_stats(0, col), res.col_stats(1, col));
                }
            } else if (col_stats == Stats::Variance) {
                for (size_t i = 0; i < count; i++) {
                    update_mean_and_variance(
                        val_data[i],
                        res.col_stats(0, col),
                        res.col_stats(1, col),
                        res.col_stats(2, col)
                    );
                }
            }
        }
    }
    uint32_t n_rows = rows();
    uint32_t n_cols = cols();
    // fix up mean + variance accounting for the zero entries
    if (col_stats == Stats::Mean) {
        auto nnz = res.col_stats.row(0);
        auto mean = res.col_stats.row(1);
        res.col_stats.row(1) = mean * nnz / n_rows;
    } else if (col_stats == Stats::Variance) {
        // m2 += mean^2 * nnz*(nrow - nnz)/nrow;
        auto nnz = res.col_stats.row(0);
        auto mean = res.col_stats.row(1);
        auto m2 = res.col_stats.row(2);
        res.col_stats.row(2) =
            (m2 + mean.square() * nnz * ((double)n_rows - nnz) / n_rows) / (n_rows - 1);
        res.col_stats.row(1) = mean * nnz / n_rows;
    }

    if (row_stats == Stats::Mean) {
        auto nnz = res.row_stats.row(0);
        auto mean = res.row_stats.row(1);
        res.row_stats.row(1) = mean * nnz / n_cols;
    } else if (row_stats == Stats::Variance) {
        // m2 += mean^2 * nnz*(ncol - nnz)/ncol;
        auto nnz = res.row_stats.row(0);
        auto mean = res.row_stats.row(1);
        auto m2 = res.row_stats.row(2);
        res.row_stats.row(2) =
            (m2 + mean.square() * nnz * ((double)n_cols - nnz) / n_cols) / (n_cols - 1);
        res.row_stats.row(1) = mean * nnz / n_cols;
    }

    return res;
}

} // end namespace BPCells

// ---- from BPCells R_interrupts.h ----

namespace {
inline void myCheckInterruptFn(void * /*dummy*/) { R_CheckUserInterrupt(); }
} // anonymous namespace

// Check for interrupts and throw the sentinel exception if one is pending
inline bool hasUserInterrupt() { return (R_ToplevelExec(myCheckInterruptFn, NULL) == FALSE); }

// Wrap a function call such that we will check for R user interrupts
// The function itself is run in a background thread, and the main thread checks for an
// R user interrupt every 100ms.
// The function must take a pointer to as std::atomic<bool> as the last argument, and exit
// early if the bool gets set to true
// NOTE: It is EXTREMELY IMPORTANT that no R objects are created/destroyed inside the spawned
// thread, which includes destructors that mess with R's GC protection.
template <class F, class... Args>
std::invoke_result_t<F, Args..., std::atomic<bool> *>
run_with_R_interrupt_check(F &&f, Args &&...args) {
    std::atomic<bool> interrupt(false);
    auto job =
        std::async(std::launch::async, std::forward<F>(f), std::forward<Args>(args)..., &interrupt);
    while (job.wait_for(std::chrono::milliseconds(100)) == std::future_status::timeout) {
        if (hasUserInterrupt()) {
            interrupt = true;
        }
    }
    if (interrupt) {
        throw Rcpp::internal::InterruptedException();
    }
    return job.get();
}

// Wrap and parallelize a function call while checking for R interrupts
// Arguments:
//  Function f:
//     parameters are: (int) task_start, (int) task_end, and pointer to std::atomic<bool>
//     Will run tasks with IDs [task_start, task_end).
//     Return value if any is ignored.
//     Should terminate early if the std::atomic<bool> becomes true
//  int n_tasks:
//     Number of tasks to run
//  int n_threads:
//     Number of threads to run with
// Parallelization assigns a contigous range of tasks to each thread
template <class F>
void run_parallel_with_R_interrupt_check(F &&f, size_t n_tasks, size_t n_threads) {
    std::vector<std::future<void>> task_vec;
    std::atomic<bool> interrupt = false;
    size_t idx = 0;
    n_threads = std::max<size_t>(1, n_threads);

    for (size_t i = 0; i < n_threads; i++) {
        size_t items = (n_tasks - idx) / (n_threads - i);
        task_vec.push_back(std::async(std::launch::async, std::forward<F>(f), idx, idx+items, &interrupt));
        idx += items;
    }
    // Wait for work to finish while checking for interrupts
    for (size_t i = 0; i < n_threads; i++) {
        while(task_vec[i].wait_for(std::chrono::milliseconds(100)) == std::future_status::timeout) {
            if (hasUserInterrupt()) {
                interrupt = true;
            }
        }
    }
    if (interrupt) throw Rcpp::internal::InterruptedException();
}

// ---- from BPCells R_xptr_wrapper.h ----

// Design problem:
// - Once Rcpp creates an object in an XPtr, it permanently takes
//   ownership and *will* eventually call the destructor
// - Since MatrixLoaders and FragmentLoaders now take ownership of their inputs via unique_ptr,
//   they will also eventually call the destructor for their inputs
// - This could result in double-free, where both R and C++ try to free an input loader
// Solution:
// - Rather than using XPtr<T>, use XPtr<unique_ptr<T>>. This works because unique_ptr is capable of
//   giving up ownership of an object. So then when R calls the destructor on the unique_ptr, the unique
//   ptr will only call the destructor on contents that it owns.
// - This pair of functions allows creating an object and giving ownership to R, then taking back ownership
//   later when passed the object from R.

// Analog of make_unique, but wrapped in an XPtr
template<class T, class... Args>
SEXP make_unique_xptr(Args&&... args) {
  return Rcpp::wrap(Rcpp::XPtr<std::unique_ptr<T>>(new std::unique_ptr<T>(new T(std::forward<Args>(args)...))));
}

// Construct the XPtr wrapper if we already have a unique_ptr object
template<class T>
SEXP unique_xptr(std::unique_ptr<T> &&ptr) {
  return Rcpp::wrap(Rcpp::XPtr<std::unique_ptr<T>>(new std::unique_ptr<T>(std::move(ptr))));
}

// Take ownership of the unique_ptr from the XPtr
template<class T>
std::unique_ptr<T> take_unique_xptr(SEXP &sexp) {
  Rcpp::XPtr<std::unique_ptr<T>> ptr(sexp);
  return std::unique_ptr<T>(ptr->release());
}

// Get a pointer to the object within the XPtr without assuming ownership
template<class T>
T* peek_unique_xptr(SEXP &sexp) {
  Rcpp::XPtr<std::unique_ptr<T>> ptr(sexp);
  return ptr->get();
}
