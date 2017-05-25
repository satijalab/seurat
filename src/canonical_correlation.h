#ifndef CANONICAL_CORRELATION
#define CANONICAL_CORRELATION

#include <RcppEigen.h>
#include <progress.hpp>
#include <cmath>
#include <unordered_map>
#include <data_manipulation.h>
//----------------------------------------------------
List CanonCor(Eigen::MatrixXd mat1, Eigen::MatrixXd mat2, bool center);
List CalcPartialCCA(Eigen::MatrixXd mat1, Eigen::MatrixXd mat2, int k, Eigen::MatrixXd v_init);
//----------------------------------------------------

#endif//CANONICAL_CORRELATION