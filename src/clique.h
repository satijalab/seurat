#ifndef CLIQUE_H
#define CLIQUE_H

#include <string>
#include <vector>
#include <ctime>
#include <queue>
#include <iostream>
#include <algorithm>
#include <Rcpp.h>

//----------------------------------------------------
class Clique{
public:
  Clique(Rcpp::IntegerVector members);
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



std::vector<Clique*> findQuasiCliques(Rcpp::NumericMatrix adj_mat, double r_param,double update);
Rcpp::NumericMatrix subsetMatrix(Rcpp::NumericMatrix m, Rcpp::IntegerVector rows, Rcpp::IntegerVector cols);
std::string labelClique(int cell);
Rcpp::IntegerVector whichNotZero(Rcpp::NumericVector x);
Rcpp::IntegerVector removeNode(Rcpp::IntegerVector x, int y);
bool removeRedundantClique(Rcpp::IntegerVector x, Rcpp::IntegerVector y);
std::vector<Clique*> mergeCliques(Rcpp::NumericMatrix adj_mat, std::vector<Clique*> cliqueList, double m_param=0.5, double q=0.01, double qup=0.1, double update=0.1);
std::string intToString (int Number );


Rcpp::List r_wrapper(Rcpp::NumericMatrix adj_mat, double r_param, double m_param=0.5, double q=0.1, double qup=0.1, double update=0.1 );
int sizeCliqueIntersection(Rcpp::IntegerVector x, Rcpp::IntegerVector y);
Rcpp::IntegerVector cliqueUnion(Rcpp::IntegerVector x, Rcpp::IntegerVector y);
double scoreCluster(Rcpp::IntegerVector, int, Rcpp::NumericMatrix );

Rcpp::NumericMatrix setRow(Rcpp::NumericMatrix m, int r, int n);
Rcpp::NumericMatrix setCol(Rcpp::NumericMatrix m, int c, int n);
Rcpp::NumericMatrix delRowCol(Rcpp::NumericMatrix m, int i);
std::vector<Clique*> findLargestCliques(std::vector<Clique*> cliqueList, int n );
//----------------------------------------------------

#endif//CLIQUE_H