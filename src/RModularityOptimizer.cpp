#include <chrono>
#include <exception>
#include <fstream>
#include <iostream>
#include <limits>
#include <memory>
#include <sstream>

#include <RcppEigen.h>
#include <Rcpp.h>
#include <progress.hpp>

#include "ModularityOptimizer.h"

using namespace ModularityOptimizer;
using namespace std::chrono;
using namespace Rcpp;


// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppProgress)]]
// [[Rcpp::export]]
IntegerVector RunModularityClusteringCpp(Eigen::SparseMatrix<double> SNN,
    int modularityFunction,
    double resolution,
    int algorithm,
    int nRandomStarts,
    int nIterations,
    int randomSeed,
    bool printOutput,
    std::string edgefilename) {

  // validate arguments
  if(modularityFunction != 1 && modularityFunction != 2)
    stop("Modularity parameter must be equal to 1 or 2.");
  if(algorithm != 1 && algorithm !=2 && algorithm !=3 && algorithm !=4)
    stop("Algorithm for modularity optimization must be 1, 2, 3, or 4");
  if(nRandomStarts < 1)
    stop("Have to have at least one start");
  if(nIterations < 1)
    stop("Need at least one interation");
  if (modularityFunction == 2 && resolution > 1.0)
    stop("error: resolution<1 for alternative modularity");
  try {
  bool update;
  double modularity, maxModularity, resolution2;
  int i, j;

  std::string msg = "Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck";
  if (printOutput)
    Rcout << msg << std::endl << std::endl;

  // Load netwrok
  std::shared_ptr<Network> network;
  if(edgefilename != "") {
    if (printOutput)
      Rcout << "Reading input file..." << std::endl << std::endl;
    try{
      network = readInputFile(edgefilename, modularityFunction);
    } catch(...) {
      stop("Could not parse edge file.");
    }
  } else {
    // Load lower triangle
    int network_size = (SNN.nonZeros() / 2) + 3;
    IVector node1;
    IVector node2;
    DVector edgeweights;
    node1.reserve(network_size);
    node2.reserve(network_size);
    edgeweights.reserve(network_size);
    for (int k=0; k < SNN.outerSize(); ++k){
      for (Eigen::SparseMatrix<double>::InnerIterator it(SNN, k); it; ++it){
        if(it.col() >= it.row()){
          continue;
        }
        node1.emplace_back(it.col());
        node2.emplace_back(it.row());
        edgeweights.emplace_back(it.value());
      }
    }
    if(node1.size() == 0) {
      stop("Matrix contained no network data.  Check format.");
    }
    int nNodes = std::max(SNN.cols(), SNN.rows());
    network = matrixToNetwork(node1, node2, edgeweights, modularityFunction, nNodes);
    Rcpp::checkUserInterrupt();
  }

  if (printOutput)
  {
    Rprintf("Number of nodes: %d\n", network->getNNodes());
    Rprintf("Number of edges: %d\n", network->getNEdges());
    Rcout << std::endl;
    Rcout << "Running " <<  ((algorithm == 1) ? "Louvain algorithm" : ((algorithm == 2) ? "Louvain algorithm with multilevel refinement" : "smart local moving algorithm")) << "...";
    Rcout << std::endl;
  }

  resolution2 = ((modularityFunction == 1) ? (resolution / (2 * network->getTotalEdgeWeight() + network->getTotalEdgeWeightSelfLinks())) : resolution);

  auto beginTime = duration_cast<milliseconds>(system_clock::now().time_since_epoch());
  std::shared_ptr<Clustering> clustering;
  maxModularity = -std::numeric_limits<double>::infinity();
  JavaRandom random(randomSeed);

  Progress p(nRandomStarts, printOutput);
  for (i = 0; i < nRandomStarts; i++)
  {
    //if (printOutput && (nRandomStarts > 1))
    //Rprintf("Random start: %d\n", i + 1);

    VOSClusteringTechnique vosClusteringTechnique(network, resolution2);

    j = 0;
    update = true;
    do
    {
      /*if (printOutput && (nIterations > 1))
        Rprintf("Iteration: %d\n", j + 1);
      */
      if (algorithm == 1)
        update = vosClusteringTechnique.runLouvainAlgorithm(random);
      else if (algorithm == 2)
        update = vosClusteringTechnique.runLouvainAlgorithmWithMultilevelRefinement(random);
      else if (algorithm == 3)
        vosClusteringTechnique.runSmartLocalMovingAlgorithm(random);
      j++;

      modularity = vosClusteringTechnique.calcQualityFunction();

      //if (printOutput && (nIterations > 1))
      //  Rprintf("Modularity: %.4f\n", modularity);
      Rcpp::checkUserInterrupt();
    }
    while ((j < nIterations) && update);

    if (modularity > maxModularity)
    {
      clustering = vosClusteringTechnique.getClustering();
      maxModularity = modularity;
    }

    /*if (printOutput && (nRandomStarts > 1))
    {
      if (nIterations == 1)
        Rprintf("Modularity: %.4f\n", modularity);
      Rcout << std::endl;
    }*/
    p.increment();
  }
  auto endTime = duration_cast<milliseconds>(system_clock::now().time_since_epoch());
  if(clustering == nullptr) {
    stop("Clustering step failed.");
  }
  if (printOutput)
  {
    if (nRandomStarts == 1)
    {
      if (nIterations > 1)
        Rcout << std::endl;
      Rprintf("Modularity: %.4f\n", maxModularity);
    }
    else
      Rprintf("Maximum modularity in %d random starts: %.4f\n", nRandomStarts, maxModularity);
      Rprintf("Number of communities: %d\n", clustering->getNClusters());
      Rprintf("Elapsed time: %d seconds\n", static_cast<int>((endTime - beginTime).count() / 1000.0));
  }

  // Return results
  clustering->orderClustersByNNodes();
  IntegerVector iv(clustering->cluster.cbegin(), clustering->cluster.cend());
  return iv;
  } catch(std::exception &ex) {
    forward_exception_to_r(ex);
  } catch(...) {
    ::Rf_error("c++ exception (unknown reason)");
  }
  return IntegerVector(1);
}
