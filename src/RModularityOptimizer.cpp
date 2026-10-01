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
#include <RcppThread.h>

#include "ModularityOptimizer.h"

using namespace ModularityOptimizer;
using namespace std::chrono;
using namespace Rcpp;

// Convert an SNN adjacency matrix into the edge-list representation consumed by
// ModularityOptimizer. Sparse Graph/dgCMatrix inputs are traversed column-wise
// through their compressed slots; dense numeric matrices are scanned over the
// lower triangle. Only lower-triangle nonzero edges are emitted because the SNN
// graph is undirected and stored symmetrically.
std::shared_ptr<Network> matrixToNetworkRcpp(SEXP SNN, int modularityFunction) {
  IVector node1;
  IVector node2;
  DVector edgeweights;
  int nNodes = 0;

  if (Rf_isS4(SNN)) {
    S4 mat(SNN);
    IntegerVector dims = mat.slot("Dim");
    if (dims.size() < 2) {
      stop("SNN matrix must have two dimensions.");
    }
    nNodes = std::max(dims[0], dims[1]);
    if (!Rf_inherits(SNN, "dgCMatrix") && !Rf_inherits(SNN, "Graph")) {
      stop("RunModularityClusteringCpp supports dgCMatrix/Graph sparse matrices or base numeric matrices.");
    }
    IntegerVector i = mat.slot("i");
    IntegerVector p = mat.slot("p");
    NumericVector x = mat.slot("x");
    node1.reserve((x.size() / 2) + 3);
    node2.reserve((x.size() / 2) + 3);
    edgeweights.reserve((x.size() / 2) + 3);
    for (int col = 0; col < dims[1]; ++col) {
      for (int idx = p[col]; idx < p[col + 1]; ++idx) {
        int row = i[idx];
        if (col >= row) {
          continue;
        }
        node1.emplace_back(col);
        node2.emplace_back(row);
        edgeweights.emplace_back(x[idx]);
      }
    }
  } else {
    NumericMatrix mat(SNN);
    nNodes = std::max(mat.nrow(), mat.ncol());
    int reserveSize = ((mat.nrow() * mat.ncol()) / 2) + 3;
    node1.reserve(reserveSize);
    node2.reserve(reserveSize);
    edgeweights.reserve(reserveSize);
    for (int col = 0; col < mat.ncol(); ++col) {
      for (int row = col + 1; row < mat.nrow(); ++row) {
        double value = mat(row, col);
        if (value == 0) {
          continue;
        }
        node1.emplace_back(col);
        node2.emplace_back(row);
        edgeweights.emplace_back(value);
      }
    }
  }

  if (node1.empty()) {
    stop("Matrix contained no network data.  Check format.");
  }
  return matrixToNetwork(node1, node2, edgeweights, modularityFunction, nNodes);
}

struct ModularityRestartWorker {
  std::shared_ptr<Network> network;
  double resolution;
  int algorithm;
  int nIterations;
  int randomSeed;
  std::vector<double>& startModularity;
  std::vector<std::shared_ptr<Clustering>>& startClustering;

  ModularityRestartWorker(std::shared_ptr<Network> network, double resolution,
                              int algorithm, int nIterations, int randomSeed,
                              std::vector<double>& startModularity,
                              std::vector<std::shared_ptr<Clustering>>& startClustering)
    : network(network), resolution(resolution), algorithm(algorithm),
      nIterations(nIterations), randomSeed(randomSeed),
      startModularity(startModularity), startClustering(startClustering) {}

  void operator()(std::size_t begin, std::size_t end, RcppThread::ProgressBar* bar = NULL) {
    // Each restart owns its VOSClusteringTechnique and JavaRandom stream, and
    // writes to restart-indexed output slots. The shared Network is read-only.
    for (std::size_t s = begin; s < end; ++s) {
      JavaRandom random(static_cast<uint64_t>(randomSeed) + static_cast<uint64_t>(s));
      VOSClusteringTechnique vos(network, resolution);
      int it = 0;
      bool upd = true;
      double mod = -std::numeric_limits<double>::infinity();
      do {
        if (algorithm == 1) upd = vos.runLouvainAlgorithm(random);
        else if (algorithm == 2) upd = vos.runLouvainAlgorithmWithMultilevelRefinement(random);
        else if (algorithm == 3) vos.runSmartLocalMovingAlgorithm(random);
        it++;
        mod = vos.calcQualityFunction();
      } while ((it < nIterations) && upd);
      startModularity[s] = mod;
      startClustering[s] = vos.getClustering();
      if (bar != NULL) {
        (*bar)++;
      }
    }
  }
};

IntegerVector RunModularityClusteringCppOnNetwork(std::shared_ptr<Network> network,
    double resolution,
    int modularityFunction,
    int algorithm,
    int nRandomStarts,
    int nIterations,
    int randomSeed,
    bool printOutput,
    int nThreads) {

  if (printOutput) {
    Rprintf("Number of nodes: %d\n", network->getNNodes());
    Rprintf("Number of edges: %d\n", network->getNEdges());
    Rcout << std::endl;
    Rcout << "Running " <<  ((algorithm == 1) ? "Louvain algorithm" : ((algorithm == 2) ? "Louvain algorithm with multilevel refinement" : "smart local moving algorithm")) << "...";
    Rcout << std::endl;
  }

  double resolution2 = ((modularityFunction == 1) ? (resolution / (2 * network->getTotalEdgeWeight() + network->getTotalEdgeWeightSelfLinks())) : resolution);
  // Avoid paying extra overhead if running single-threaded or with a single random start
  const bool serialRestarts = nThreads <= 1 || nRandomStarts <= 1;

  auto beginTime = duration_cast<milliseconds>(system_clock::now().time_since_epoch());
  std::shared_ptr<Clustering> clustering;
  double maxModularity = -std::numeric_limits<double>::infinity();

  if (serialRestarts) {
    JavaRandom random(randomSeed);
    Progress p(nRandomStarts, printOutput);
    for (int i = 0; i < nRandomStarts; i++) {
      VOSClusteringTechnique vosClusteringTechnique(network, resolution2);
      int j = 0;
      bool update = true;
      double modularity = -std::numeric_limits<double>::infinity();
      do {
        if (algorithm == 1)
          update = vosClusteringTechnique.runLouvainAlgorithm(random);
        else if (algorithm == 2)
          update = vosClusteringTechnique.runLouvainAlgorithmWithMultilevelRefinement(random);
        else if (algorithm == 3)
          vosClusteringTechnique.runSmartLocalMovingAlgorithm(random);
        j++;
        modularity = vosClusteringTechnique.calcQualityFunction();
        Rcpp::checkUserInterrupt();
      } while ((j < nIterations) && update);
      if (modularity > maxModularity) {
        clustering = vosClusteringTechnique.getClustering();
        maxModularity = modularity;
      }
      p.increment();
    }
  } else {
    int nWorkers = std::min(nThreads, nRandomStarts);
    std::vector<double> startModularity(nRandomStarts, -std::numeric_limits<double>::infinity());
    std::vector<std::shared_ptr<Clustering>> startClustering(nRandomStarts);
    ModularityRestartWorker worker(network, resolution2, algorithm, nIterations, randomSeed,
                                       startModularity, startClustering);
    if (printOutput) {
      RcppThread::ProgressBar bar(nRandomStarts, 1);
      RcppThread::parallelFor(0, nRandomStarts, [&](int restart) {
        worker(restart, restart + 1, &bar);
      }, nWorkers);
    } else {
      RcppThread::parallelFor(0, nRandomStarts, [&](int restart) {
        worker(restart, restart + 1);
      }, nWorkers);
    }
    for (int i = 0; i < nRandomStarts; i++) {
      if (startModularity[i] > maxModularity) {
        maxModularity = startModularity[i];
        clustering = startClustering[i];
      }
    }
  }

  auto endTime = duration_cast<milliseconds>(system_clock::now().time_since_epoch());
  if(clustering == nullptr) {
    stop("Clustering step failed.");
  }
  if (printOutput) {
    if (nRandomStarts == 1) {
      if (nIterations > 1)
        Rcout << std::endl;
      Rprintf("Modularity: %.4f\n", maxModularity);
    } else {
      Rprintf("Maximum modularity in %d random starts: %.4f\n", nRandomStarts, maxModularity);
    }
    Rprintf("Number of communities: %d\n", clustering->getNClusters());
    Rprintf("Elapsed time: %d seconds\n", static_cast<int>((endTime - beginTime).count() / 1000.0));
  }

  clustering->orderClustersByNNodes();
  IntegerVector iv(clustering->cluster.cbegin(), clustering->cluster.cend());
  return iv;
}

// [[Rcpp::depends(RcppThread)]]
// [[Rcpp::export]]
IntegerVector RunModularityClusteringCpp(SEXP SNN,
    int modularityFunction,
    double resolution,
    int algorithm,
    int nRandomStarts,
    int nIterations,
    int randomSeed,
    bool printOutput,
    std::string edgefilename = "",
    int nThreads = 1) {

  // Validate optimizer arguments, build a Network from the SNN graph, run one or
  // more random starts, and return cluster IDs ordered by community size. The
  // serial branch uses one JavaRandom stream across starts; the parallel branch
  // assigns a deterministic stream to each restart and reduces the best
  // modularity in restart order.
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

  if (printOutput) {
    Rcout << "Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck"
          << std::endl << std::endl;
  }

  std::shared_ptr<Network> network;
  if (edgefilename != "") {
    if (printOutput) {
      Rcout << "Reading input file..." << std::endl << std::endl;
    }
    try {
      network = readInputFile(edgefilename, modularityFunction);
    } catch (...) {
      stop("Could not parse edge file.");
    }
  } else {
    network = matrixToNetworkRcpp(SNN, modularityFunction);
  }
  Rcpp::checkUserInterrupt();

  return RunModularityClusteringCppOnNetwork(
    network,
    resolution,
    modularityFunction,
    algorithm,
    nRandomStarts,
    nIterations,
    randomSeed,
    printOutput,
    nThreads
  );
}

// [[Rcpp::depends(RcppThread)]]
// [[Rcpp::export]]
List RunModularityClusteringCpp_multi(SEXP SNN,
    int modularityFunction,
    NumericVector resolutions,
    int algorithm,
    int nRandomStarts,
    int nIterations,
    int randomSeed,
    bool printOutput,
    std::string edgefilename = "",
    int nThreads = 1) {

  if(modularityFunction != 1 && modularityFunction != 2)
    stop("Modularity parameter must be equal to 1 or 2.");
  if(algorithm != 1 && algorithm !=2 && algorithm !=3 && algorithm !=4)
    stop("Algorithm for modularity optimization must be 1, 2, 3, or 4");
  if(nRandomStarts < 1)
    stop("Have to have at least one start");
  if(nIterations < 1)
    stop("Need at least one interation");
  if (modularityFunction == 2) {
    for (int i = 0; i < resolutions.size(); i++) {
      if (resolutions[i] > 1.0)
        stop("error: resolution<1 for alternative modularity");
    }
  }

  if (printOutput) {
    Rcout << "Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck"
          << std::endl << std::endl;
  }

  std::shared_ptr<Network> network;
  if (edgefilename != "") {
    if (printOutput) {
      Rcout << "Reading input file..." << std::endl << std::endl;
    }
    try {
      network = readInputFile(edgefilename, modularityFunction);
    } catch (...) {
      stop("Could not parse edge file.");
    }
  } else {
    network = matrixToNetworkRcpp(SNN, modularityFunction);
  }
  Rcpp::checkUserInterrupt();

  List results(resolutions.size());
  for (int i = 0; i < resolutions.size(); i++) {
    results[i] = RunModularityClusteringCppOnNetwork(
      network,
      resolutions[i],
      modularityFunction,
      algorithm,
      nRandomStarts,
      nIterations,
      randomSeed,
      printOutput,
      nThreads
    );
  }
  return results;
}
