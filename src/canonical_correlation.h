#ifndef CANONICAL_CORRELATION
#define CANONICAL_CORRELATION

#include <RcppEigen.h>
#include <progress.hpp>
#include <cmath>
#include <unordered_map>
#include <data_manipulation.h>
//----------------------------------------------------
List CanonCor(Eigen::MatrixXd mat1, Eigen::MatrixXd mat2, bool center);
List SparseCanonCor(Eigen::MatrixXd mat1, Eigen::MatrixXd mat2, bool standardize, int k){
//----------------------------------------------------

#endif//CANONICAL_CORRELATION