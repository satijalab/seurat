#include <RcppEigen.h>
#include <progress.hpp>
#include <cmath>
#include <unordered_map>
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

