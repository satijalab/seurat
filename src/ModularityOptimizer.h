#pragma once

#include <chrono>
#include <exception>
#include <fstream>
#include <limits>
#include <memory>
#include <sstream>
#include <iostream>
#include <vector>

typedef std::vector<int> IVector;
typedef std::vector<double> DVector;

namespace ModularityOptimizer {

class JavaRandom {
private:
  uint64_t seed;
  int next(int bits);
public:
  JavaRandom(uint64_t seed);
  int nextInt(int n);
  void setSeed(uint64_t seed);
};

namespace Arrays2 {
  IVector generateRandomPermutation(int nElements);
  IVector generateRandomPermutation(int nElements, JavaRandom& random);
}


class Clustering {
private:
  int nNodes;
public:
  // Note: These two variables were "protected" in java, which means it is accessible to the whole package/public.
  // Although we could have used friend classes, this allows for better mirroring of the original code.
  int nClusters;
  IVector cluster;
  
  Clustering(int nNodes);
  Clustering(IVector cluster);
  int getNNodes() const {return nNodes;};
  int getNClusters() const {return nClusters;};
  IVector getClusters() const {return cluster;};
  int getCluster(int node) const {return cluster[node];};
  IVector getNNodesPerCluster() const;
  std::vector<IVector> getNodesPerCluster() const;
  void setCluster(int node, int cluster);
  void initSingletonClusters();
  void orderClustersByNNodes();
  void mergeClusters(const Clustering& clustering);

};

class Network {
  friend class VOSClusteringTechnique;
protected:
  int nNodes;
  int nEdges;
  DVector nodeWeight;
  IVector firstNeighborIndex;
  IVector neighbor;
  DVector edgeWeight;
  double totalEdgeWeightSelfLinks;
public:
  Network();
  Network(int nNodes, DVector* nodeWeight, std::vector<IVector>& edge, DVector* edgeWeight);
  Network(int nNodes, std::vector<IVector>& edge) : 
    Network(nNodes, nullptr, edge, nullptr) { };
  Network(int nNodes, DVector* nodeWeight, std::vector<IVector> edge) : 
    Network(nNodes, nodeWeight, edge, nullptr) {};
  Network(int nNodes, std::vector<IVector>& edge, DVector* edgeWeight) :
    Network(nNodes, nullptr, edge, edgeWeight) {};
  
  Network(int nNodes, DVector* nodeWeight, IVector& firstNeighborIndex, IVector& neighbor, DVector* edgeWeight);
  Network(int nNodes, IVector& firstNeighborIndex, IVector& neighbor) : 
    Network(nNodes, nullptr, firstNeighborIndex, neighbor, nullptr) {};
  
  Network(int nNodes, DVector* nodeWeight, IVector& firstNeighborIndex, IVector& neighbor) :
    Network(nNodes, nodeWeight, firstNeighborIndex, neighbor, nullptr){};
  
  Network(int nNodes, IVector& firstNeighborIndex, IVector& neighbor, DVector* edgeWeight) :
    Network(nNodes, nullptr, firstNeighborIndex, neighbor, edgeWeight) {};
  

 int getNNodes() {return nNodes;};
 double getTotalNodeWeight();
 DVector getNodeWeights();
 double getNodeWeight(int node) { return nodeWeight.at(node);};
 int getNEdges() {return nEdges / 2;};
 int getNEdges(int node) {return firstNeighborIndex.at(node + 1) - firstNeighborIndex.at(node);};
 IVector getNEdgesPerNode();
 std::vector<IVector> getEdges();
 IVector getEdges(int node);
 std::vector<IVector> getEdgesPerNode();
 double getTotalEdgeWeight();
 double getTotalEdgeWeight(int node);
 DVector getTotalEdgeWeightPerNode();
 DVector getEdgeWeights() {return edgeWeight;};
 DVector getEdgeWeights(int node);
 std::vector<DVector> getEdgeWeightsPerNode();
 double getTotalEdgeWeightSelfLinks()
 {
   return totalEdgeWeightSelfLinks;
 };
 // Added these to avoid making these values public
 int getFirstNeighborIndexValue(int i) const {
   return firstNeighborIndex.at(i);
 };
 int getNeighborValue(int index) const {
   return neighbor.at(index);
 }

 std::vector<Network> createSubnetworks(Clustering clustering) const;
 Network createReducedNetwork(const Clustering& clustering) const;
 Clustering identifyComponents();
private:
  double generateRandomNumber(int node1, int node2, const IVector& nodePermutation);
  Network createSubnetwork(const Clustering& clustering, int cluster, IVector& node, 
                           IVector& subnetworkNode, IVector& subnetworkNeighbor, DVector& subnetworkEdgeWeight) const;
};


class VOSClusteringTechnique {
private:
  std::shared_ptr<Network> network;
  std::shared_ptr<Clustering> clustering;
  double resolution;

public:
  VOSClusteringTechnique(std::shared_ptr<Network> network, double resolution);
  VOSClusteringTechnique(std::shared_ptr<Network> network, std::shared_ptr<Clustering> clustering, double resolution);
  std::shared_ptr<Network> getNetwork() { return network;}
  std::shared_ptr<Clustering> getClustering()  { return clustering; }
  double getResolution() {return resolution; }
  void setNetwork(std::shared_ptr<Network> network) {this->network = network;}
  void setClustering(std::shared_ptr<Clustering> clustering) {this->clustering = clustering;}
  void setResolution(double resolution) {this->resolution = resolution;}
  double calcQualityFunction();

  bool runLocalMovingAlgorithm(JavaRandom& random);
  bool runLouvainAlgorithm(JavaRandom& random);
  bool runIteratedLouvainAlgorithm(int maxNIterations, JavaRandom& random);
  bool runLouvainAlgorithmWithMultilevelRefinement(JavaRandom& random);
  bool runIteratedLouvainAlgorithmWithMultilevelRefinement(int maxNIterations, JavaRandom& random);
  bool runSmartLocalMovingAlgorithm(JavaRandom& random);
  bool runIteratedSmartLocalMovingAlgorithm(int nIterations, JavaRandom& random);

  int removeCluster(int cluster);
  void removeSmallClusters(int minNNodesPerCluster);
};


std::shared_ptr<Network> matrixToNetwork(IVector& node1, IVector& node2, DVector& edgeWeight1, int modularityFunction, int nNodes);
std::shared_ptr<Network> readInputFile(std::string fname, int modularityFunction);
std::vector<std::string> split(const std::string& s, char delimiter);
};
