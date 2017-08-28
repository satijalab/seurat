#ifndef DATA_MANIPULATION
#define DATA_MANIPULATION

#include <RcppEigen.h>
#include <progress.hpp>
#include <cmath>
#include <unordered_map>

using namespace Rcpp;

//----------------------------------------------------
Eigen::SparseMatrix<double> RunUMISampling(Eigen::SparseMatrix<double> data, int sample_val,
                                           bool upsample, bool display_progress);
Eigen::SparseMatrix<double> RunUMISamplingPerCell(Eigen::SparseMatrix<double> data,
                                                  NumericVector sample_val, bool upsample,
                                                  bool display_progress);
Eigen::SparseMatrix<double> RowMergeMatrices(Eigen::SparseMatrix<double, Eigen::RowMajor> mat1,
                                             Eigen::SparseMatrix<double, Eigen::RowMajor> mat2,
                                             std::vector< std::string > mat1_rownames,
                                             std::vector< std::string > mat2_rownames,
                                             std::vector< std::string > all_rownames);
Eigen::SparseMatrix<double> LogNorm(Eigen::SparseMatrix<double> data, int scale_factor,
                                    bool display_progress );
Eigen::MatrixXd FastMatMult(Eigen::MatrixXd m1, Eigen::MatrixXd m2);
Eigen::MatrixXd FastRowScale(Eigen::MatrixXd mat, bool scale, bool center, double scale_max,
                             bool display_progress);
Eigen::MatrixXd FastSparseRowScale(Eigen::SparseMatrix<double> mat, bool scale, bool center,
                                   double scale_max, bool display_progress);
Eigen::MatrixXd FastCov(Eigen::MatrixXd mat, bool center);
Eigen::MatrixXd FastCovMats(Eigen::MatrixXd mat1, Eigen::MatrixXd mat2, bool center);
Eigen::MatrixXd Standardize(Eigen::MatrixXd mat, bool display_progress);
Eigen::MatrixXd FastRBind(Eigen::MatrixXd mat1, Eigen::MatrixXd mat2);
Eigen::VectorXd ExpMeanFast(Eigen::MatrixXd mat, bool display_progress));
Eigen::VectorXd FastLogVMR(Eigen::SparseMatrix<double> mat, bool display_progress);
//----------------------------------------------------

#endif//DATA_MANIPULATION
