#include <RcppEigen.h>
#include <progress.hpp>
#include <cmath>
using namespace Rcpp;


// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppProgress)]]



// [[Rcpp::export]]
Eigen::SparseMatrix<double> RunUMISampling(Eigen::SparseMatrix<double> data, int sample_val, bool upsample = false, bool display_progress=true){
    Progress p(data.outerSize(), display_progress);
    Eigen::VectorXd colSums = data.transpose() * Eigen::VectorXd::Ones(data.cols());
    for (int k=0; k < data.outerSize(); ++k){
      p.increment();
      for (Eigen::SparseMatrix<double>::InnerIterator it(data, k); it; ++it){
        double entry = it.value();
        if( (upsample) || (colSums[k] > sample_val)){
          entry = entry * double(sample_val) / colSums[k];
          if (fmod(entry, 1) != 0){
            double rn = runif(1)[0];
            if(fmod(entry, 1) <= rn){
              data.coeffRef(it.row(), it.col()) = floor(entry);
            }
            else{
              data.coeffRef(it.row(), it.col()) = ceil(entry);
            }
          }
          else{
            data.coeffRef(it.row(), it.col()) = entry;
          }
        }
      }
    } 
  return(data);
}

// [[Rcpp::export]]
Eigen::SparseMatrix<double> RunUMISamplingPerCell(Eigen::SparseMatrix<double> data, NumericVector sample_val, bool upsample = false, bool display_progress=true){
  Progress p(data.outerSize(), display_progress);
  Eigen::VectorXd colSums = data.transpose() * Eigen::VectorXd::Ones(data.cols());
  for (int k=0; k < data.outerSize(); ++k){
    p.increment();
    for (Eigen::SparseMatrix<double>::InnerIterator it(data, k); it; ++it){
      double entry = it.value();
      if( (upsample) || (colSums[k] > sample_val[k])){
        entry = entry * double(sample_val[k]) / colSums[k];
        if (fmod(entry, 1) != 0){
          double rn = runif(1)[0];
          if(fmod(entry, 1) <= rn){
            data.coeffRef(it.row(), it.col()) = floor(entry);
          }
          else{
            data.coeffRef(it.row(), it.col()) = ceil(entry);
          }
        }
        else{
          data.coeffRef(it.row(), it.col()) = entry;
        }
      }
    }
  } 
  return(data);
}