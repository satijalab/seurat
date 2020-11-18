#include <RcppEigen.h>
#include "data_manipulation.h"
#include <progress.hpp>
#include <cmath>
#include <unordered_map>
#include <fstream>
#include <string>
#include <iomanip>

using namespace Rcpp;
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppProgress)]]

typedef Eigen::Triplet<double> T;
// [[Rcpp::export(rng = false)]]
Eigen::SparseMatrix<double> ComputeSNN(Eigen::MatrixXd nn_ranked, double prune) {
  std::vector<T> tripletList;
  int k = nn_ranked.cols();
  tripletList.reserve(nn_ranked.rows() * nn_ranked.cols());
  for(int j=0; j<nn_ranked.cols(); ++j){
    for(int i=0; i<nn_ranked.rows(); ++i) {
      tripletList.push_back(T(i, nn_ranked(i, j) - 1, 1));
    }
  }
  Eigen::SparseMatrix<double> SNN(nn_ranked.rows(), nn_ranked.rows());
  SNN.setFromTriplets(tripletList.begin(), tripletList.end());
  SNN = SNN * (SNN.transpose());
  for (int i=0; i < SNN.outerSize(); ++i){
    for (Eigen::SparseMatrix<double>::InnerIterator it(SNN, i); it; ++it){
      it.valueRef() = it.value()/(k + (k - it.value()));
      if(it.value() < prune){
        it.valueRef() = 0;
      }
    }
  }
  SNN.prune(0.0); // actually remove pruned values
  return SNN;
}

// [[Rcpp::export(rng = false)]]
void WriteEdgeFile(Eigen::SparseMatrix<double> snn, String filename, bool display_progress){
  if (display_progress == true) {
    Rcpp::Rcerr << "Writing SNN as edge file" << std::endl;
  }
  // Write out lower triangle
  std::ofstream output;
  output.open(filename);
  Progress p(snn.outerSize(), display_progress);
  for (int k=0; k < snn.outerSize(); ++k){
    p.increment();
    for (Eigen::SparseMatrix<double>::InnerIterator it(snn, k); it; ++it){
      if(it.col() >= it.row()){
        continue;
      }
      output << std::setprecision(15) << it.col() << "\t" << it.row() << "\t" << it.value() << "\n";
    }
  }
  output.close();
}

// Wrapper function so that we don't have to go back into R before writing to file
// [[Rcpp::export(rng = false)]]
Eigen::SparseMatrix<double> DirectSNNToFile(Eigen::MatrixXd nn_ranked,
                                            double prune, bool display_progress,
                                            String filename) {
  Eigen::SparseMatrix<double> SNN = ComputeSNN(nn_ranked, prune);
  WriteEdgeFile(SNN, filename, display_progress);
  return SNN;
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
std::vector<double> SNN_SmallestNonzero_Dist(
    Eigen::SparseMatrix<double> snn,
    Eigen::MatrixXd mat,
    int n,
    std::vector<double> nearest_dist
) {
  std::vector<double> results;
  for (int i=0; i < snn.outerSize(); ++i){
    // create vectors to store the nonzero snn elements and their indices
    std::vector<double> nonzero;
    std::vector<size_t> nonzero_idx;
    for (Eigen::SparseMatrix<double>::InnerIterator it(snn, i); it; ++it) {
      nonzero.push_back(it.value());
      nonzero_idx.push_back(it.row());
    }
    std::vector<size_t> nonzero_order = sort_indexes(nonzero);
    int n_i = n;
    if (n_i > nonzero_order.size()) n_i = nonzero_order.size();
    std::vector<double> dists;
    for (int j = 0; j < nonzero_order.size(); ++j) {
      // compute euclidean distances to cells with small edge weights
      // if multiple entries have same value as nth element, calc dist to all
      size_t cell = nonzero_idx[nonzero_order[j]];
      if(dists.size() < n_i  || nonzero[nonzero_order[j]] == nonzero[nonzero_order[n_i-1]]) {
        double res = (mat.row(cell) - mat.row(i)).norm();
        if (nearest_dist[i] > 0) {
          res = res - nearest_dist[i];
          if (res < 0) res = 0;
        }
        dists.push_back(res);
      } else {
        break;
      }
    }
    double avg_dist;
    if (dists.size() > n_i) {
      std::sort(dists.rbegin(), dists.rend());
      avg_dist = std::accumulate(dists.begin(), dists.begin() + n_i, 0.0) / n_i;
    } else {
      avg_dist = std::accumulate(dists.begin(), dists.end(), 0.0) / dists.size();
    }
    results.push_back(avg_dist);
  }
  return results;
}
