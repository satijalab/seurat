#ifndef CORRECT_EXPRESSION
#define CORRECT_EXPRESSION

#include <RcppEigen.h>
#include <unordered_map>

using namespace Rcpp;

//----------------------------------------------------
Eigen::SparseMatrix<double> FindWeightsC(
  NumericVector cells2,
  Eigen::MatrixXd distances,
  std::vector<std::string> anchor_cells2,
  std::vector<std::string> integration_matrix_rownames,
  Eigen::MatrixXd cell_index,
  Eigen::VectorXd anchor_score,
  double min_dist,
  double sd,
  bool display_progress
);
Eigen::SparseMatrix<double> IntegrateDataC(
  Eigen::SparseMatrix<double> integration_matrix,
  Eigen::SparseMatrix<double> weights,
  Eigen::SparseMatrix<double> expression_cells2
);
std::vector<double> ScoreHelper(
    Eigen::SparseMatrix<double> snn,
    Eigen::MatrixXd query_pca,
    Eigen::MatrixXd query_dists,
    Eigen::MatrixXd corrected_nns,
    int k_snn,
    bool subtract_first_nn,
    bool display_progress
);
//----------------------------------------------------

#endif//CORRECT_EXPRESSION
