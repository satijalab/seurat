#include <RcppEigen.h>
#include "data_manipulation.h"

using namespace Rcpp;


// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppProgress)]]


// [[Rcpp::export]]
Eigen::MatrixXd CanonCor(Eigen::MatrixXd mat1, Eigen::MatrixXd mat2, bool center = true){
  Eigen::MatrixXd cmat1 = FastCov(mat1);
  Eigen::MatrixXd cmat2 = FastCov(mat2);
  Eigen::LDLT<Eigen::MatrixXd> cr1;
  Eigen::LDLT<Eigen::MatrixXd> cr2;
  
  cr1.compute(mat1);
  cr2.compute(mat2);
  
  return(mat1);
}