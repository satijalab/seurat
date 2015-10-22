#include "clique.h"
using namespace std;

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

Clique::Clique(IntegerVector m, int i){
  members = m;
  idx = i;
}


Clique::Clique(IntegerVector m){
  members = m;
  idx = 0;
}

Clique::Clique(){
  IntegerVector m;
  members = m;
  idx = 0;
}

void Clique::emptyMemberVector(){
  IntegerVector m;
  members = m;
}

IntegerVector Clique::getMembers() const{
  return members;
}

int Clique::getSize() const {
  return members.size();
}


void Clique::setMembers(IntegerVector m){
  members = m;
}

void Clique::setIdx(int i){
  idx = i;
}

int Clique::getIdx(){
  return idx;
}

void Clique::removeMember(int m){
  for(int i = 0; i<members.size();++i){
    if(m == members[i]) members.erase(members.begin()+i);
  }
}

bool cmpI (Clique* c1,  Clique* c2) {
  return c1->getMembers().size() < c2->getMembers().size();
}

bool cmpD (Clique* c1,  Clique* c2) {
  return c1->getMembers().size() > c2->getMembers().size();
}

bool Clique::operator ()(const Clique* c1, const Clique* c2) const
{
  return c1->getMembers().size() > c2->getMembers().size();
}

struct cp{
  Clique* c1;
  Clique* c2;
  int combined_size;
};

bool cmpCP(cp pair1, cp pair2)
{
  return pair1.combined_size < pair2.combined_size;
}

// [[Rcpp::export]]

IntegerVector whichNotZero(NumericVector x) {
  IntegerVector v = Rcpp::seq(0, x.size()-1);
  return v[x!=0];
}

// [[Rcpp::export]]

NumericMatrix subsetMatrix(NumericMatrix m, IntegerVector rows, IntegerVector cols){
  NumericMatrix subMat(rows.size(),cols.size());
  for(int r=0; r<rows.size(); ++r){
    for(int c=0; c<cols.size(); ++c){
      subMat(r,c)= m(rows(r),cols(c));
    }
  }
  return subMat;
}

// [[Rcpp::export]]

IntegerVector removeNode(IntegerVector x, int y){
  as<std::vector<int> >(x);
  x.erase(y);
  return wrap(x);
}

// [[Rcpp::export]]
bool removeRedundantClique(IntegerVector x, IntegerVector y){
  as<std::vector<int> >(x);
  as<std::vector<int> >(y);
  std::sort(x.begin(), x.end());
  std::sort(y.begin(), y.end());
  vector<int> z;
  set_intersection(x.begin(),x.end(),y.begin(),y.end(),back_inserter(z));
  return z.size() == x.size();
}

void constructIndex(vector<Clique*> cliqueList){
  for(int i=0;i<cliqueList.size();++i){
    cliqueList[i]->setIdx(i);
  }
}

std::vector<Clique*> findLargestCliques(vector<Clique*> cliqueList, int n ){
  vector<Clique*> sortedCliques = cliqueList;
  vector<Clique*> largestCliques(n);
  /*only need to sort the largest n cliques*/
  nth_element(sortedCliques.begin(), sortedCliques.begin()+n, sortedCliques.end(), cmpD );
  partial_sort_copy(sortedCliques.begin(), sortedCliques.begin()+n, largestCliques.begin(), largestCliques.end(), cmpD);
  return largestCliques;
}

// [[Rcpp::export]]
int sizeCliqueIntersection(IntegerVector x, IntegerVector y){
  as<std::vector<int> >(x);
  as<std::vector<int> >(y);
  std::sort(x.begin(), x.end());
  std::sort(y.begin(), y.end());
  vector<int> z;
  set_intersection(x.begin(),x.end(),y.begin(),y.end(),back_inserter(z));
  return z.size();
}

IntegerVector cliqueUnion(IntegerVector x, IntegerVector y){
  as<std::vector<int> >(x);
  as<std::vector<int> >(y);
  std::sort(x.begin(), x.end());
  std::sort(y.begin(), y.end());
  vector<int> z;
  set_union(x.begin(),x.end(),y.begin(),y.end(),back_inserter(z));
  return wrap(z);
}


// [[Rcpp::export]]

NumericMatrix setRow(NumericMatrix m, int r, int n){
  for(int i=0; i<m.ncol();++i){
    m(r,i) = n;
  }
  return m;
}

// [[Rcpp::export]]
NumericMatrix setCol(NumericMatrix m, int c, int n){
  for(int i=0; i<m.ncol();++i){
    m(i,c) = n;
  }
  return m;
}

double scoreCluster(IntegerVector cluster, int cell, NumericMatrix adj_mat){
  double score;
  for(int i=0;i<cluster.size(); ++i){
    score = score + adj_mat(cell, cluster[i]);
  }
  score = (score - 1)/cluster.size(); 
  return score;
}

/*Functions for dealing with sparse matrix*/
NumericVector subviewToVec(arma::SpSubview<double> sbv, int n){
  NumericVector v(n);
  for(int i=0; i<n; ++i){
    v[i] = sbv[i];
  }
  return v;
}

NumericMatrix subsetMatrix(sp_mat m, IntegerVector rows, IntegerVector cols){
  NumericMatrix subMat(rows.size(),cols.size());
  for(int r=0; r<rows.size(); ++r){
    for(int c=0; c<cols.size(); ++c){
      subMat(r,c)= m(rows(r),cols(c));
    }
  }
  return subMat;
}

double scoreCluster(IntegerVector cluster, int cell, sp_mat adj_mat){
  double score;
  for(int i=0;i<cluster.size(); ++i){
    score = score + adj_mat(cell, cluster[i]);
  }
  score = (score - 1)/cluster.size(); 
  return score;
}

//' @export 
// [[Rcpp::export]]
List r_wrapper(NumericMatrix adj_mat, arma::sp_mat adj_mat_sp, double r_param, double m_param, double q, double qup, double update, int min_cluster_size, bool do_sparse){
  vector<Clique*> cliqueList = findQuasiCliques(adj_mat, adj_mat_sp, r_param, update, do_sparse);
  cliqueList = mergeCliques(adj_mat, adj_mat_sp, cliqueList, m_param, q, qup, update, min_cluster_size, do_sparse);
  list<IntegerVector > cliqueListFinal;
  list<int> cliqueSizesFinal;
  list<bool> unassigned;
  for(int i =0; i<cliqueList.size(); ++i){
    cliqueListFinal.push_back(cliqueList[i]->getMembers()+1);
    cliqueSizesFinal.push_back(cliqueList[i]->getSize());
  }
  if(cliqueList[cliqueList.size()]->getIdx()==-1){
    unassigned.push_back(true);
  }
  return List::create(cliqueListFinal, cliqueSizesFinal, unassigned);
}


vector<Clique*> findQuasiCliques(NumericMatrix adj_mat, arma::sp_mat adj_mat_sp, double r_param, double update, bool do_sparse){
  int num_cells;
  if(do_sparse)  num_cells = adj_mat_sp.n_rows;
  else num_cells = adj_mat.nrow();
  if(num_cells==0) Rcpp::stop("error: adjacency matrix required");
  /*Creates a vector to hold all the Clique objects,  
  storing as a vector of pointers prevents object slicing*/
  vector<Clique*> cliqueList;
  int counter=1;
  for (int i = 0; i<num_cells; ++i){
    if ((update>0) &&(i==int(num_cells*update*counter))) {
      Rcout << "Find Quasi-Cliques : Processed " << i << " cells" <<endl;
      counter++;
    }
    bool in_clique = false;
    
    NumericVector r;
    if(do_sparse){
      SpSubview<double> spr = adj_mat_sp.row(i);
       r = subviewToVec(spr, num_cells);
    }
    else  r = adj_mat.row(i);
    /* possible improvement here is to limit the number of nodes to compare to only the nearest X */
    IntegerVector nodes = whichNotZero(r);
    NumericMatrix sub_adj_mat(nodes.size(), nodes.size()); 
    if(nodes.size()>1){   
      /* Subsetting matrices isn't well supported. */
      if(do_sparse) sub_adj_mat = subsetMatrix(adj_mat_sp, nodes, nodes);
      else sub_adj_mat = subsetMatrix(adj_mat, nodes, nodes);
    }
    else{
      in_clique = true;
    }
    while(!in_clique){
      if(sub_adj_mat.nrow() >= 3){
        NumericVector rsums(sub_adj_mat.nrow());
        for(int r = 0; r<sub_adj_mat.nrow();r++){
          NumericVector rowR = sub_adj_mat.row(r);
          rowR = rowR[rowR!=0];
          rsums(r) = rowR.size()-1;
        }
        int least_connected = which_min(rsums);
        if(rsums(least_connected)/sub_adj_mat.nrow() < r_param){
          nodes = removeNode(nodes, least_connected);
          if(do_sparse) sub_adj_mat = subsetMatrix(adj_mat_sp, nodes, nodes);
          else sub_adj_mat = subsetMatrix(adj_mat, nodes, nodes);
        }
        else {
          in_clique = true;
          cliqueList.push_back(new Clique(nodes));
        }
      }
      else{
        in_clique = true;
      }
    }
  }
  /* cmpI is a custom comparator for Clique objects - sorts in increasing order */
  sort(cliqueList.begin(), cliqueList.end(), cmpI); 
  for(int i = 0; i<cliqueList.size(); ++i){ 
    for(int j = i+1; j<cliqueList.size(); ++j){
      if (removeRedundantClique(cliqueList[i]->getMembers(), cliqueList[j]->getMembers())){ 
        cliqueList.erase(cliqueList.begin()+j);
        --j;
      }
    } 
  } 
  return cliqueList;
}


/*Function to merge quasi-cliques identified in findQuasiClique() 
 *based on clique overlap
 */

vector<Clique*> mergeCliques(NumericMatrix adj_mat, arma::sp_mat adj_mat_sp, vector<Clique*> cliqueList, double m_param, double q, double qup, double update, int min_cluster_size, bool do_sparse){
  int num_cells;
  if(do_sparse)  num_cells = adj_mat_sp.n_rows;
  else num_cells = adj_mat.nrow();
  /* count singletons */
  int num_singletons = 0;
  int counter=1;
  for(int i =0; i<cliqueList.size(); ++i){
    if (cliqueList[i]->getSize() <= 1) num_singletons++;
  }
  constructIndex(cliqueList);
  NumericMatrix allComparisons(cliqueList.size(), cliqueList.size());

  while(q<=1){
    bool allChecked = false;
    while(!allChecked){
      /* determine cliques to compare */
      vector<Clique*> largestCliques;
      largestCliques = findLargestCliques(cliqueList, (cliqueList.size()-num_singletons)*q);
      
      /*create vector of pairs of cliques that could possibly merge */
      vector<cp> cliquePairs;
      for(int i =0; i<largestCliques.size()-1; ++i){
        int idx1 = largestCliques[i]->getIdx();
        for(int j =i+1; j<largestCliques.size(); ++j){
          int idx2 = largestCliques[j]->getIdx();
          if(allComparisons(idx1, idx2)==0){
            double overlap = sizeCliqueIntersection(largestCliques[i]->getMembers(), largestCliques[j]->getMembers())/double(min(largestCliques[i]->getSize(), largestCliques[j]->getSize()));
            if(overlap > m_param){
              cliquePairs.push_back(cp());
              cliquePairs[cliquePairs.size()-1].c1 = largestCliques[i]; cliquePairs[cliquePairs.size()-1].c2 = largestCliques[j];
              cliquePairs[cliquePairs.size()-1].combined_size = largestCliques[i]->getSize() + largestCliques[j]->getSize();
            }
            else{
              allComparisons(idx1, idx2) = 1;
              allComparisons(idx2, idx1) = 1;  
            }
          }
        }
      }
      /* just need to find max if made overlap comp already */
      if (cliquePairs.size() > 0){
        cp cliquesToMerge = *max_element(cliquePairs.begin(), cliquePairs.end(), cmpCP);
        int idx1 = cliquesToMerge.c1->getIdx();
        allComparisons=setRow(allComparisons,idx1, 0); allComparisons=setCol(allComparisons, idx1, 0); 
        cliquesToMerge.c1->setMembers(cliqueUnion(cliquesToMerge.c1->getMembers(),cliquesToMerge.c2->getMembers()));
        for(int i = 0; i<cliqueList.size(); ++i){
          if(cliqueList[i] == cliquesToMerge.c2){
            cliqueList.erase(cliqueList.begin()+i);
            break;
          }
        }
      }
      else allChecked = true;
      int processedCliques=num_cells-cliqueList.size();
      if ((update>0) &&(processedCliques>int(num_cells*update*counter))) {
        Rcout << "Merge Quasi-Cliques :  " << cliqueList.size() << " cliques remain" <<endl;
        counter++;
      }
    }
    q=q+qup;
  }

  /* remove clusters that are smaller than minimum requirement*/
  for(int i=0; i<cliqueList.size(); ++i){
    if(cliqueList[i]->getSize() < min_cluster_size){
      cliqueList.erase(cliqueList.begin()+i);
      --i;
    }
  }
  
  /* ensure cells are only assigned to one cluster 
   * It is possible cells to belong to no cluster, 
   * assign all of those to a separate 'ignore' cluster
   */
  
  IntegerVector ignore;
  for(int i=0;i<num_cells;++i){
    typedef std::multimap<double, int, greater<double> > MapType;
    MapType clusterRepeats;
    for(int j=0;j<cliqueList.size();++j){
      if(sizeCliqueIntersection(IntegerVector::create(i),cliqueList[j]->getMembers()) == 1){
        if(do_sparse) clusterRepeats.insert(MapType::value_type(scoreCluster(cliqueList[j]->getMembers(),i,adj_mat_sp),j));
        else clusterRepeats.insert(MapType::value_type(scoreCluster(cliqueList[j]->getMembers(),i,adj_mat),j));
      }
    }
    if(clusterRepeats.size()>1){
      clusterRepeats.erase(clusterRepeats.begin());
      for (multimap<double, int>::iterator it= clusterRepeats.begin(); it != clusterRepeats.end(); ++it){
        cliqueList[it->second]->removeMember(i);
      }
    }
    if(clusterRepeats.size()<1){
      ignore.push_back(i);
    }
  }
  if(ignore.size()>0) cliqueList.push_back(new Clique(ignore,-1));
  if(update>0) Rcout << "Merge Quasi-Cliques :  " << cliqueList.size() << " cliques remain" <<endl;
  return cliqueList;
}