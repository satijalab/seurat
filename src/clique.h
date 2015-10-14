#ifndef CLIQUE_H
#define CLIQUE_H

#include <string>
#include <vector>
#include <ctime>
#include <queue>
#include <iostream>
#include <algorithm>
#include <RcppArmadillo.h>

//----------------------------------------------------
class Clique{
public:
  Clique(Rcpp::IntegerVector m, int i);
  Clique(Rcpp::IntegerVector m);
  Clique();
  Rcpp::IntegerVector getMembers() const;
  int getSize() const;
  void setIdx(int i);
  int getIdx();
  void emptyMemberVector();
  void setMembers(Rcpp::IntegerVector m);
  void removeMember(int m);
  bool cmpI (Clique* c1, Clique* c2);
  bool cmpD (Clique* c1,  Clique* c2);
  bool operator ()(const Clique* c1, const Clique* c2) const;
private:
  Rcpp::IntegerVector members;
  int idx;
};

Rcpp::List r_wrapper(Rcpp::NumericMatrix adj_mat, arma::sp_mat adj_mat_sp, double r_param, double m_param, double q, double qup, double update, int min_cluster_size, bool do_sparse);
std::vector<Clique*> findQuasiCliques(Rcpp::NumericMatrix adj_mat, arma::sp_mat adj_mat_sp, double r_param,double update, bool do_sparse);
std::vector<Clique*> mergeCliques(Rcpp::NumericMatrix adj_mat, arma::sp_mat adj_mat_sp, std::vector<Clique*> cliqueList, double m_param, double q, double qup, double update, int min_cluster_size, bool do_sparse);

Rcpp::IntegerVector whichNotZero(Rcpp::NumericVector x);
Rcpp::NumericMatrix subsetMatrix(Rcpp::NumericMatrix m, Rcpp::IntegerVector rows, Rcpp::IntegerVector cols);
Rcpp::IntegerVector removeNode(Rcpp::IntegerVector x, int y);
bool removeRedundantClique(Rcpp::IntegerVector x, Rcpp::IntegerVector y);
void constructIndex(std::vector<Clique*> cliqueList);
std::vector<Clique*> findLargestCliques(std::vector<Clique*> cliqueList, int n );
int sizeCliqueIntersection(Rcpp::IntegerVector x, Rcpp::IntegerVector y);
Rcpp::IntegerVector cliqueUnion(Rcpp::IntegerVector x, Rcpp::IntegerVector y);
double scoreCluster(Rcpp::IntegerVector, int, Rcpp::NumericMatrix );
Rcpp::NumericMatrix setRow(Rcpp::NumericMatrix m, int r, int n);
Rcpp::NumericMatrix setCol(Rcpp::NumericMatrix m, int c, int n);

//----------------------------------------------------

#endif//CLIQUE_H