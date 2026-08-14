// Copyright 2022 BPCells contributors
//
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

#pragma once

#define RCPP_NO_RTTI
#define RCPP_NO_SUGAR
#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]

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
