// Copyright 2021 BPCells contributors
//
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

#pragma once

#include <atomic>
#include <cstdint>
#include <memory>
#include <stdexcept>
#include <vector>

#define RCPP_NO_RTTI
#define RCPP_NO_SUGAR
#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]

#include "MatrixStats.h"

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

#include "MatrixOps.h"
