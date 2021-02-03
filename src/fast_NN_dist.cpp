#include <Rcpp.h>
using namespace Rcpp;

// code in-parts taken from https://gallery.rcpp.org/articles/parallel-distance-matrix/
// Assumption: the end position of vector2 is implied by the end position of vector1
// generic function for euclidean distance
template <typename InputIterator1, typename InputIterator2>
inline double euclidean_distance(InputIterator1 begin1, 
                                 InputIterator1 end1, 
                                 InputIterator2 begin2) {
  // value to return
  double rval = 0;
  
  // set iterators to beginning of ranges
  InputIterator1 it1 = begin1;
  InputIterator2 it2 = begin2;
  
  // for each input item
  while (it1 != end1) {
    
    // take the value and increment the iterator
    double d1 = *it1++;
    double d2 = *it2++;
    
    // update the distance
    rval += pow(d1-d2, 2);
  }
  
  return sqrt(rval);
}


// [[Rcpp::export]]
List fast_dist(NumericMatrix x, NumericMatrix y, List n) {
  // extracting the number of element in the knn graph
  size_t ngraph_size = n.size();
  if (x.nrow() != ngraph_size) { return List(); }
  List distances_list = clone(n);
  
  // looping over the neigbors
  for (size_t i=0; i<ngraph_size; ++i) {
    // extracting the neighbors vector
    SEXP ll = n[i];
    NumericVector neighbors(ll);
    
    // prepopulating the distance vector
    NumericVector distances(neighbors.size());
    
    // extracting the row iterator of x
    NumericMatrix::Row row1 = x.row(i);
    for (size_t j=0; j<neighbors.size(); ++j) {
      size_t n_idx = neighbors[j] - 1;
      
      // extracting the row iterator of y
      NumericMatrix::Row row2 = y.row(n_idx);
      // extracting the distance
      double distance = euclidean_distance(row1.begin(),
                                           row1.end(),
                                           row2.begin());
      
      if (distance == -1) { return List(); }
      distances[j] = distance;
    }
    
    // updating the distance
    distances_list[i] = distances;
  }
  
  return distances_list;
}