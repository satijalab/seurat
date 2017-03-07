#ifndef CANONICAL_CORRELATION
#define CANONICAL_CORRELATION

#include <RcppEigen.h>
#include <progress.hpp>
#include <cmath>
#include <unordered_map>
#include <data_manipulation.h>

//----------------------------------------------------
Eigen::MatrixXd CanonCor(Eigen::MatrixXd mat1, Eigen::MatrixXd mat2, bool center);
//----------------------------------------------------

#endif//CANONICAL_CORRELATION