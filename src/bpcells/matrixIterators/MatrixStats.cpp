// Copyright 2022 BPCells contributors
// 
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

#include "MatrixStats.h"

namespace BPCells {

Eigen::ArrayXXd StatsResult::rowNonzeros() {
    if (row_stats.rows() < 1)
        throw std::runtime_error("Nonzero not calculated in this StatsResult");
    return row_stats.row(0);
}
Eigen::ArrayXXd StatsResult::rowMean() {
    if (row_stats.rows() < 2) throw std::runtime_error("Mean not calculated in this StatsResult");
    return row_stats.row(1);
}
Eigen::ArrayXXd StatsResult::rowVariance() {
    if (row_stats.rows() < 3)
        throw std::runtime_error("Variance not calculated in this StatsResult");
    return row_stats.row(2);
}

Eigen::ArrayXXd StatsResult::colNonzeros() {
    if (col_stats.rows() < 1)
        throw std::runtime_error("Nonzero not calculated in this StatsResult");
    return col_stats.row(0);
}
Eigen::ArrayXXd StatsResult::colMean() {
    if (col_stats.rows() < 2) throw std::runtime_error("Mean not calculated in this StatsResult");
    return col_stats.row(1);
}
Eigen::ArrayXXd StatsResult::colVariance() {
    if (col_stats.rows() < 3)
        throw std::runtime_error("Variance not calculated in this StatsResult");
    return col_stats.row(2);
}

StatsResult StatsResult::transpose() { return StatsResult{col_stats, row_stats}; }

} // end namespace BPCells
