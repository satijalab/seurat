#include <RcppEigen.h>
#include <progress.hpp>
#include <cmath>
#include <unordered_map>
#include <fstream>
#include <string>

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


typedef Eigen::Triplet<double> T;
// [[Rcpp::export]]
Eigen::SparseMatrix<double> RowMergeMatrices(Eigen::SparseMatrix<double, Eigen::RowMajor> mat1, Eigen::SparseMatrix<double, Eigen::RowMajor> mat2, std::vector< std::string > mat1_rownames,
                                             std::vector< std::string > mat2_rownames, std::vector< std::string > all_rownames){


  // Set up hash maps for rowname based lookup
  std::unordered_map<std::string, int> mat1_map;
  for(int i = 0; i < mat1_rownames.size(); i++){
    mat1_map[mat1_rownames[i]] = i;
  }
  std::unordered_map<std::string, int> mat2_map;
  for(int i = 0; i < mat2_rownames.size(); i++){
    mat2_map[mat2_rownames[i]] = i;
  }

  // set up tripletList for new matrix creation
  std::vector<T> tripletList;
  int num_rows = all_rownames.size();
  int num_col1 = mat1.cols();
  int num_col2 = mat2.cols();


  tripletList.reserve(mat1.nonZeros() + mat2.nonZeros());
  for(int i = 0; i < num_rows; i++){
    std::string key = all_rownames[i];
    if (mat1_map.count(key)){
      for(Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it1(mat1, mat1_map[key]); it1; ++it1){
        tripletList.push_back(T(i, it1.col(), it1.value()));
      }
    }
    if (mat2_map.count(key)){
      for(Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it2(mat2, mat2_map[key]); it2; ++it2){
        tripletList.push_back(T(i, num_col1 + it2.col(), it2.value()));
      }
    }
  }
  Eigen::SparseMatrix<double> combined_mat(num_rows, num_col1 + num_col2);
  combined_mat.setFromTriplets(tripletList.begin(), tripletList.end());
  return combined_mat;
}

// [[Rcpp::export]]
Eigen::SparseMatrix<double> LogNorm(Eigen::SparseMatrix<double> data, int scale_factor, bool display_progress = true){
  Progress p(data.outerSize(), display_progress);
  Eigen::VectorXd colSums = data.transpose() * Eigen::VectorXd::Ones(data.cols());
  for (int k=0; k < data.outerSize(); ++k){
    p.increment();
    for (Eigen::SparseMatrix<double>::InnerIterator it(data, k); it; ++it){
      data.coeffRef(it.row(), it.col()) = log1p(double(it.value()) / colSums[k] * scale_factor);
    }
  }
  return data;
}

// [[Rcpp::export]]
Eigen::MatrixXd FastMatMult(Eigen::MatrixXd m1, Eigen::MatrixXd m2){
  Eigen::MatrixXd m3 = m1 * m2;
  return(m3);
}


/* Performs row scaling and/or centering. Equivalent to using t(scale(t(mat))) in R.
   Note: Doesn't handle NA/NaNs in the same way the R implementation does, */

// [[Rcpp::export]]
Eigen::MatrixXd FastRowScale(Eigen::MatrixXd mat, bool scale = true, bool center = true,
                             double scale_max = 10, bool display_progress = true){
  Progress p(mat.rows(), display_progress);
  Eigen::MatrixXd scaled_mat(mat.rows(), mat.cols());
  for(int i=0; i < mat.rows(); ++i){
    p.increment();
    Eigen::ArrayXd r = mat.row(i).array();
    double rowMean = r.mean();
    double rowSdev = 1;
    if(scale == true){
      if(center == true){
        rowSdev = sqrt((r - rowMean).square().sum() / (mat.cols() - 1));
      }
      else{
        rowSdev = sqrt(r.square().sum() / (mat.cols() - 1));
      }
    }
    if(center == false){
      rowMean = 0;
    }
    scaled_mat.row(i) = (r - rowMean) / rowSdev;
    for(int s=0; s<scaled_mat.row(i).size(); ++s){
      if(scaled_mat(i, s) > scale_max){
        scaled_mat(i, s) = scale_max;
      }
    }
  }
  return scaled_mat;
}

/* Performs column scaling and/or centering. Equivalent to using scale(mat, TRUE, apply(x,2,sd)) in R.
 Note: Doesn't handle NA/NaNs in the same way the R implementation does, */

// [[Rcpp::export]]
Eigen::MatrixXd Standardize(Eigen::MatrixXd mat, bool display_progress = true){
  Progress p(mat.cols(), display_progress);
  Eigen::MatrixXd std_mat(mat.rows(), mat.cols());
  for(int i=0; i < mat.cols(); ++i){
    p.increment();
    Eigen::ArrayXd r = mat.col(i).array();
    double colMean = r.mean();
    double colSdev = sqrt((r - colMean).square().sum() / (mat.rows() - 1));
    std_mat.col(i) = (r - colMean) / colSdev;
  }
  return std_mat;
}

// [[Rcpp::export]]
Eigen::MatrixXd FastSparseRowScale(Eigen::SparseMatrix<double> mat, bool scale = true, bool center = true,
                                   double scale_max = 10, bool display_progress = true){
  mat = mat.transpose();
  Progress p(mat.outerSize(), display_progress);
  Eigen::MatrixXd scaled_mat(mat.rows(), mat.cols());
  for (int k=0; k<mat.outerSize(); ++k){
    p.increment();
    double colMean = 0;
    double colSdev = 0;
    for (Eigen::SparseMatrix<double>::InnerIterator it(mat,k); it; ++it)
    {
      colMean += it.value();
    }
    colMean = colMean / mat.rows();
    if (scale == true){
      int nnZero = 0;
      if(center == true){
        for (Eigen::SparseMatrix<double>::InnerIterator it(mat,k); it; ++it)
        {
          nnZero += 1;
          colSdev += pow((it.value() - colMean), 2);
        }
        colSdev += pow(colMean, 2) * (mat.rows() - nnZero);
      }
      else{
        for (Eigen::SparseMatrix<double>::InnerIterator it(mat,k); it; ++it)
        {
          colSdev += pow(it.value(), 2);
        }
      }
      colSdev = sqrt(colSdev / (mat.rows() - 1));
    }
    else{
      colSdev = 1;
    }
    if(center == false){
      colMean = 0;
    }
    Eigen::VectorXd col = Eigen::VectorXd(mat.col(k));
    scaled_mat.col(k) = (col.array() - colMean) / colSdev;
    for(int s=0; s<scaled_mat.col(k).size(); ++s){
      if(scaled_mat(s,k) > scale_max){
        scaled_mat(s,k) = scale_max;
      }
    }
  }
  return scaled_mat.transpose();
}

/* Note: May not handle NA/NaNs in the same way the R implementation does, */

// [[Rcpp::export]]
Eigen::MatrixXd FastCov(Eigen::MatrixXd mat, bool center = true){
  if (center) {
    mat = mat.rowwise() - mat.colwise().mean();
  }
  Eigen::MatrixXd cov = (mat.adjoint() * mat) / double(mat.rows() - 1);
  return(cov);
}

// [[Rcpp::export]]
Eigen::MatrixXd FastCovMats(Eigen::MatrixXd mat1, Eigen::MatrixXd mat2, bool center = true){
  if(center){
    mat1 = mat1.rowwise() - mat1.colwise().mean();
    mat2 = mat2.rowwise() - mat2.colwise().mean();
  }
  Eigen::MatrixXd cov = (mat1.adjoint() * mat2) / double(mat1.rows() - 1);
  return(cov);
}

/* Note: Faster than the R implementation but is not in-place */
//[[Rcpp::export]]
Eigen::MatrixXd FastRBind(Eigen::MatrixXd mat1, Eigen::MatrixXd mat2){
  Eigen::MatrixXd mat3(mat1.rows() + mat2.rows(), mat1.cols());
  mat3 << mat1, mat2;
  return(mat3);
}

/* Calculates the row means of the logged values in non-log space */
//[[Rcpp::export]]
Eigen::VectorXd FastExpMean(Eigen::SparseMatrix<double> mat, bool display_progress){
  int ncols = mat.cols();
  Eigen::VectorXd rowmeans(mat.rows());
  mat = mat.transpose();
  if(display_progress == true){
    Rcpp::Rcerr << "Calculating gene means" << std::endl;
  }
  Progress p(mat.outerSize(), display_progress);
  for (int k=0; k<mat.outerSize(); ++k){
    p.increment();
    double rm = 0;
    for (Eigen::SparseMatrix<double>::InnerIterator it(mat,k); it; ++it){
      rm += expm1(it.value());
    }
    rm = rm / ncols;
    rowmeans[k] = log1p(rm);
  }
  return(rowmeans);
}

/* Calculate the variance to mean ratio (VMR) in non-logspace (return answer in
log-space) */
//[[Rcpp::export]]
Eigen::VectorXd FastLogVMR(Eigen::SparseMatrix<double> mat,  bool display_progress){
  int ncols = mat.cols();
  Eigen::VectorXd rowdisp(mat.rows());
  mat = mat.transpose();
  if(display_progress == true){
    Rcpp::Rcerr << "Calculating gene variance to mean ratios" << std::endl;
  }
  Progress p(mat.outerSize(), display_progress);
  for (int k=0; k<mat.outerSize(); ++k){
    p.increment();
    double rm = 0;
    double v = 0;
    int nnZero = 0;
    for (Eigen::SparseMatrix<double>::InnerIterator it(mat,k); it; ++it){
      rm += expm1(it.value());
    }
    rm = rm / ncols;
    for (Eigen::SparseMatrix<double>::InnerIterator it(mat,k); it; ++it){
      v += pow(expm1(it.value()) - rm, 2);
      nnZero += 1;
    }
    v = (v + (ncols - nnZero) * pow(rm, 2)) / (ncols - 1);
    rowdisp[k] = log(v/rm);

  }
  return(rowdisp);
}

int IntersectLength(std::vector<int> a, std::vector<int> b){
  std::unordered_set<int> s(a.begin(), a.end());
  int intersect = count_if(b.begin(), b.end(), [&](int k) {return s.find(k) != s.end();});
  return(intersect);
}

// IMPORTANT: assumes that a and b are vectors with non-duplicated elements
int UnionLength(std::vector<int> a, std::vector<int> b, int intersect_length) {
  return(a.size() + b.size() - intersect_length);
}

std::vector<int> ToVector(Eigen::VectorXd v1){
  std::vector<int> v2(v1.data(), v1.data() + v1.size());
  return(v2);
}
