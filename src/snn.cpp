#include <RcppEigen.h>
#include "data_manipulation.h"
#include <progress.hpp>
#include <cmath>
#include <unordered_map>
#include <fstream>
#include <string>

using namespace Rcpp;
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppProgress)]]

typedef Eigen::Triplet<double> T;
//[[Rcpp::export]]
Eigen::SparseMatrix<double> ComputeSNN(Eigen::MatrixXd nn_large, Eigen::MatrixXd nn_ranked,
                                       double prune, bool display_progress) {
  Progress p(nn_large.rows(), display_progress);
  std::vector<T> tripletList;
  tripletList.reserve(nn_large.rows() * nn_large.cols()/2);
  for(int i=0; i<nn_large.rows(); ++i){
    p.increment();
    for(int j=0; j<nn_large.cols(); ++j){
      std::vector<int> cell1 = ToVector(nn_ranked.row(i));
      std::vector<int> cell2 = ToVector(nn_ranked.row(nn_large(i, j) - 1));
      int s = IntersectLength(cell1, cell2);
      int u = UnionLength(cell1, cell2, s);
      double e = double(s) / u;
      if(e > prune){
        tripletList.push_back(T(i, nn_large(i, j) - 1, e));
      }
    }
    tripletList.push_back(T(i, i, 1));
  }
  Eigen::SparseMatrix<double> SNN(nn_large.rows(), nn_large.rows());
  SNN.setFromTriplets(tripletList.begin(), tripletList.end());
  return SNN;
}

//[[Rcpp::export]]
void WriteEdgeFile(Eigen::SparseMatrix<double> snn, String filename, bool display_progress){
  if (display_progress == true) {
    Rcpp::Rcerr << "Writing SNN as edge file" << std::endl;
  }
  Progress p(snn.outerSize(), display_progress);
  std::ofstream output;
  output.open(filename);
  for (int k=0; k < snn.outerSize(); ++k){
    p.increment();
    for (Eigen::SparseMatrix<double>::InnerIterator it(snn, k); it; ++it){
      if(it.col() >= it.row()){
        continue;
      }
      output << it.col() << "\t" << it.row() << "\t" << it.value() << "\n";
    }
  }
  output.close();
}

//[[Rcpp::export]]
Eigen::SparseMatrix<double> DirectSNNToFile(Eigen::MatrixXd nn_large, Eigen::MatrixXd nn_ranked,
                                         double prune, bool display_progress, String filename) {
  Eigen::SparseMatrix<double> SNN = ComputeSNN(nn_large, nn_ranked, prune, display_progress);
  WriteEdgeFile(SNN, filename, display_progress);
  return SNN;
}
