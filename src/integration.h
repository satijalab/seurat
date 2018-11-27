#ifndef CORRECT_EXPRESSION
#define CORRECT_EXPRESSION

#include <RcppEigen.h>
#include <unordered_map>

using namespace Rcpp;

//----------------------------------------------------
Eigen::SparseMatrix<double> FindWeightsC(
  Eigen::SparseMatrix<double> integration_matrix, 
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
Eigen::SparseMatrix<double> SNNAnchor(
  Eigen::SparseMatrix<double> k_matrix, 
  Eigen::SparseMatrix<double> anchor_only
);
//----------------------------------------------------

#endif//CORRECT_EXPRESSION
