#include <RcppEigen.h>
#include <progress.hpp>
#include <unordered_map>


using namespace Rcpp;
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppProgress)]]

typedef Eigen::Triplet<double> T;

// [[Rcpp::export]]
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
) {
  std::vector<T> tripletList;
  tripletList.reserve(anchor_cells2.size() * 10);
  std::unordered_map<int, std::vector<int>> cell_map;
  Progress p(anchor_cells2.size() + cells2.size() , display_progress);
  // build map from anchor_cells2 to integration_matrix rows
  for(int i=0; i<anchor_cells2.size(); ++i){
    std::vector<int> matches;
    std::vector<std::string>::iterator iter = integration_matrix_rownames.begin();
    while ((iter = std::find(iter, integration_matrix_rownames.end(), anchor_cells2[i])) != integration_matrix_rownames.end()) {
      int idx = std::distance(integration_matrix_rownames.begin(), iter);
      matches.push_back(idx);
      iter++;
    }
    cell_map[i] = matches;
    p.increment();
  }

  // Construct dist_weights matrix
  for(auto const &cell : cells2){
    Eigen::VectorXd dist = distances.row(cell);
    Eigen::VectorXd indices = cell_index.row(cell);
    for(int i=0; i<indices.size(); ++i){
      std::vector<int> mnn_idx = cell_map[indices[i]-1];
      for(int j=0; j<mnn_idx.size(); ++j){
        double to_add = 1 - exp(-1 * dist[i] * anchor_score[mnn_idx[j]]/2 * pow(1/sd, 2));
        tripletList.push_back(T(mnn_idx[j], cell, to_add));
      }
    }
    p.increment();
  }
  Eigen::SparseMatrix<double> return_mat;
  if(min_dist == 0){
    Eigen::SparseMatrix<double> dist_weights(integration_matrix.rows(), cells2.size());
    dist_weights.setFromTriplets(tripletList.begin(), tripletList.end(), [] (const double&, const double &b) { return b; });
    Eigen::VectorXd colSums = dist_weights.transpose() * Eigen::VectorXd::Ones(dist_weights.rows());
    for (int k=0; k < dist_weights.outerSize(); ++k){
      for (Eigen::SparseMatrix<double>::InnerIterator it(dist_weights, k); it; ++it){
        it.valueRef() = it.value()/colSums[k];
      }
    }
    return_mat = dist_weights;
  } else {
    Eigen::MatrixXd dist_weights = Eigen::MatrixXd::Constant(integration_matrix.rows(), cells2.size(), min_dist);
    for(int i = 0; i < dist_weights.cols(); ++i){
      for(int j = 0; j < dist_weights.rows(); ++j){
        dist_weights(j, i) = 1 - exp(-1 * dist_weights(j, i) * anchor_score[j]/2 * pow(1/sd, 2));
      }
    }
    for(auto const &weight : tripletList){
      dist_weights(weight.row(), weight.col()) = weight.value();
    }
    Eigen::VectorXd colSums = dist_weights.colwise().sum();
    for(int i = 0; i < dist_weights.cols(); ++i){
      for(int j = 0; j < dist_weights.rows(); ++j){
        dist_weights(j, i) = dist_weights(j, i) / colSums[i];
      }
    }
    return_mat = dist_weights.sparseView();
  }
  return(return_mat);
}


// [[Rcpp::export]]
Eigen::SparseMatrix<double> IntegrateDataC(
  Eigen::SparseMatrix<double> integration_matrix,
  Eigen::SparseMatrix<double> weights,
  Eigen::SparseMatrix<double> expression_cells2
) {
  Eigen::SparseMatrix<double> corrected = expression_cells2 - weights.transpose() * integration_matrix;
  return(corrected);
}


//[[Rcpp::export]]
Eigen::SparseMatrix<double> SNNAnchor(
  Eigen::SparseMatrix<double> k_matrix, 
  Eigen::SparseMatrix<double> anchor_only
) {
  typedef Eigen::SparseMatrix<double, 0, std::ptrdiff_t> SpMat;
  SpMat mat2 = k_matrix;
  SpMat mat3 = mat2 * mat2.transpose();
  for (int k=0; k<anchor_only.outerSize(); ++k){
    for (Eigen::SparseMatrix<double>::InnerIterator it(anchor_only,k); it; ++it){
      it.valueRef() = mat3.coeff(it.row(), it.col());
    }
  }
  return(anchor_only);
}
