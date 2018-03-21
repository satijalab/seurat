#ifndef SNN
#define SNN

#include <RcppEigen.h>
#include "data_manipulation.h"
#include <progress.hpp>
#include <cmath>
#include <unordered_map>
#include <fstream>
#include <string>
#include <iomanip>

using namespace Rcpp;

//----------------------------------------------------
Eigen::SparseMatrix<double> ComputeSNN(Eigen::MatrixXd nn_ranked);
void WriteEdgeFile(Eigen::SparseMatrix<double> snn, String filename, bool display_progress);
Eigen::SparseMatrix<double> DirectSNNToFile(Eigen::MatrixXd nn_ranked, double prune, bool display_progress, String filename);
//----------------------------------------------------

#endif//SNN
