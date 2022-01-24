#include <Rcpp.h>

using namespace Rcpp;

// the following code in-parts taken from sparseMatrixStats (http://www.bioconductor.org/packages/release/bioc/html/sparseMatrixStats.html).
// [[Rcpp::export]]
NumericVector row_sum_dgcmatrix(NumericVector &x, IntegerVector &i, int rows, int cols) {
  NumericVector rowsum(rows, 0.0);
  int x_length = x.length();
  for (int k=0; k<x_length; ++k) {
    rowsum[i[k]] += x[k];
  }
  
  return rowsum;
}

// [[Rcpp::export]]
NumericVector row_mean_dgcmatrix(NumericVector &x, IntegerVector &i, int rows, int cols) {
  NumericVector rowmean = row_sum_dgcmatrix(x, i, rows, cols);
  for (int k=0; k<rows; ++k) {
    rowmean[k] /= cols;
  }

  return rowmean;
}

// [[Rcpp::export]]
NumericVector row_var_dgcmatrix(NumericVector &x, IntegerVector &i, int rows, int cols) {
  NumericVector rowmean = row_mean_dgcmatrix(x, i, rows, cols);

  int x_length = x.length();
  NumericVector rowvar(rows, 0.0);
  IntegerVector nzero(rows, cols);
  for (int k=0; k<x_length; ++k) {
    rowvar[i[k]] += pow(x[k] - rowmean[i[k]], 2);
    nzero[i[k]] -= 1;
  }
  for (int k=0; k<rows; ++k) {
    rowvar[k] = (rowvar[k] + (pow(rowmean[k], 2) * nzero[k])) / (cols - 1);
  }
  return rowvar;
}
