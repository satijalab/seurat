#include <RcppEigen.h>
#include <progress.hpp>
#include <unordered_map>
#include "data_manipulation.h"

using namespace Rcpp;
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppProgress)]]

typedef Eigen::Triplet<double> T;

// [[Rcpp::export(rng = false)]]
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
    int k=0; //number of anchors used so far; a cell in the neighbor list may contribute to multiple anchors
    for(int i=0; i<indices.size() && k<indices.size(); ++i){ //index in neighbor list
      std::vector<int> mnn_idx = cell_map[indices[i]-1];
      for(int j=0; j<mnn_idx.size() && k<indices.size(); ++j){
        double to_add = 1 - exp(-1 * dist[i] * anchor_score[mnn_idx[j]]/ pow(2/sd, 2));
        tripletList.push_back(T(mnn_idx[j], cell, to_add));
        k++;
      }
    }
    p.increment();
  }
  Eigen::SparseMatrix<double> return_mat;
  if(min_dist == 0){
    Eigen::SparseMatrix<double> dist_weights(integration_matrix_rownames.size(), cells2.size());
    dist_weights.setFromTriplets(tripletList.begin(), tripletList.end(), [] (const double&, const double &b) { return b; });
    Eigen::VectorXd colSums = dist_weights.transpose() * Eigen::VectorXd::Ones(dist_weights.rows());
    for (int k=0; k < dist_weights.outerSize(); ++k){
      for (Eigen::SparseMatrix<double>::InnerIterator it(dist_weights, k); it; ++it){
        it.valueRef() = it.value()/colSums[k];
      }
    }
    return_mat = dist_weights;
  } else {
    Eigen::MatrixXd dist_weights = Eigen::MatrixXd::Constant(integration_matrix_rownames.size(), cells2.size(), min_dist);
    for(int i = 0; i < dist_weights.cols(); ++i){
      for(int j = 0; j < dist_weights.rows(); ++j){
        dist_weights(j, i) = 1 - exp(-1 * dist_weights(j, i) * anchor_score[j]/ pow(2/sd, 2) );
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

// [[Rcpp::export(rng = false)]]
Eigen::SparseMatrix<double> IntegrateDataC(
    Eigen::SparseMatrix<double> integration_matrix,
    Eigen::SparseMatrix<double> weights,
    Eigen::SparseMatrix<double> expression_cells2
) {
  Eigen::SparseMatrix<double> corrected = expression_cells2 - weights.transpose() * integration_matrix;
  return(corrected);
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

// [[Rcpp::export]]
std::vector<double> ScoreHelper(
    Eigen::SparseMatrix<double> snn,
    Eigen::MatrixXd query_pca,
    Eigen::MatrixXd query_dists,
    Eigen::MatrixXd corrected_nns,
    int k_snn,
    bool subtract_first_nn,
    bool display_progress
) {
  std::vector<double> scores;
  // Loop over all query cells
  Progress p(snn.outerSize(), display_progress);
  for (int i=0; i < snn.outerSize(); ++i){
    p.increment();
    // create vectors to store the nonzero snn elements and their indices
    std::vector<double> nonzero;
    std::vector<size_t> nonzero_idx;
    for (Eigen::SparseMatrix<double>::InnerIterator it(snn, i); it; ++it) {
      nonzero.push_back(it.value());
      nonzero_idx.push_back(it.row());
    }
    // find the k_snn cells with the smallest non-zero edge weights to use in
    // computing the transition probability bandwidth
    std::vector<size_t> nonzero_order = sort_indexes(nonzero);
    std::vector<double> bw_dists;
    int k_snn_i = k_snn;
    if (k_snn_i > nonzero_order.size()) k_snn_i = nonzero_order.size();
    for (int j = 0; j < nonzero_order.size(); ++j) {
      // compute euclidean distances to cells with small edge weights
      size_t cell = nonzero_idx[nonzero_order[j]];
      if(bw_dists.size() < k_snn_i || nonzero[nonzero_order[j]] == nonzero[nonzero_order[k_snn_i-1]]) {
        double res = (query_pca.col(cell) - query_pca.col(i)).norm();
        bw_dists.push_back(res);
      } else {
        break;
      }
    }
    // compute bandwidth as the mean distance of the farthest k_snn cells
    double bw;
    if (bw_dists.size() > k_snn_i) {
      std::sort(bw_dists.rbegin(), bw_dists.rend());
      bw = std::accumulate(bw_dists.begin(), bw_dists.begin() + k_snn_i, 0.0) / k_snn_i;
    } else {
      bw = std::accumulate(bw_dists.begin(), bw_dists.end(), 0.0) / bw_dists.size();
    }
    // compute transition probabilites
    double first_neighbor_dist;
    // subtract off distance to first neighbor?
    if (subtract_first_nn) {
      first_neighbor_dist = query_dists(i, 1);
    } else {
      first_neighbor_dist = 0;
    }
    bw = bw - first_neighbor_dist;
    double q_tps = 0;
    for(int j = 0; j < query_dists.cols(); ++j) {
      q_tps += std::exp(-1 * (query_dists(i, j) - first_neighbor_dist) / bw);
    }
    q_tps = q_tps/(query_dists.cols());
    double c_tps = 0;
    for(int j = 0; j < corrected_nns.cols(); ++j) {
      c_tps += exp(-1 * ((query_pca.col(i) - query_pca.col(corrected_nns(i, j)-1)).norm() - first_neighbor_dist) / bw);
    }
    c_tps = c_tps/(corrected_nns.cols());
    scores.push_back(c_tps/q_tps);
  }
  return(scores);
}

