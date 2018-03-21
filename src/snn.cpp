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
//[[Rcpp::export]]
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

//[[Rcpp::export]]
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
//[[Rcpp::export]]
Eigen::SparseMatrix<double> DirectSNNToFile(Eigen::MatrixXd nn_ranked,
                                            double prune, bool display_progress,
                                            String filename) {
  Eigen::SparseMatrix<double> SNN = ComputeSNN(nn_ranked, prune);
  WriteEdgeFile(SNN, filename, display_progress);
  return SNN;
}
