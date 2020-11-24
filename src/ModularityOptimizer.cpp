// This code is a translation of the Java package, ModularityOptimizer.jar,
// available from http://www.ludowaltman.nl/slm/ that performs clustering and was
// used by an earlier version of Seurat.
//
// In translating the code, the interface was maintained, such that the although
// the programming languages changed, the results are identical. As a
// consequence, rather than rewriting the program to be more idiomatic Rcpp/C++,
// both the output and the code was kept as close to the original Java as
// possible. For example, in order to maintain compatibility, the Java random
// number generator was implemented in C++, and std::stable_sort was used in
// place of std::sort to match the output results of the Java code exactly.
//
// In order to test and verify that the C++ and Java code return the same
// results, the C++ code added also includes a mirror of the original command
// line interface. This version of the program can be compiled with a command
// such as:
//
//    clang++ -O3 -std=c++11 -DSTANDALONE -Wall -g ModularityOptimizer.cpp
//
// And can then be called in an identical fashion to the original Java version to
// verify compatibility or identify any issues in the input/output.

#include "ModularityOptimizer.h"

#include <algorithm>
#include <exception>
#include <functional>
#include <numeric>
#include <stdexcept>

using namespace ModularityOptimizer;
using namespace std::chrono;

JavaRandom::JavaRandom(uint64_t seed) {
  setSeed(seed);
}

void JavaRandom::setSeed(uint64_t seed) {

  this->seed = (seed ^ uint64_t(0x5DEECE66D)) & ((uint64_t(1) << 48) - 1);
}

int JavaRandom::next(int bits) {
  // Only 31 bits ever used.
  seed = (seed * uint64_t(0x5DEECE66D) + uint64_t(0xB)) & ((uint64_t(1) << 48) - 1);
  return static_cast<int>(seed >> (48 - bits));
}


int JavaRandom::nextInt(int n) {
  if (n <= 0)
    throw std::out_of_range("n must be positive");
  if ((n & -n) == n) // i.e., n is a power of 2
    return static_cast<int>((static_cast<uint64_t>(n) * static_cast<uint64_t>(next(31))) >> 31);
  int bits, val;
  do
  {
    bits = next(31);
    val = bits % n;
  }
  while (bits - val + (n - 1) < 0);
  return val;
}


IVector Arrays2::generateRandomPermutation(int nElements, JavaRandom& random)
{
    IVector permutation(nElements, 0);
    for (int i = 0; i < nElements; i++)
      permutation[i] = i;
    for (int i = 0; i < nElements; i++)
    {
      int j = random.nextInt(nElements);
      int k = permutation[i];
      permutation[i] = permutation[j];
      permutation[j] = k;
    }
    return permutation;
}


Clustering::Clustering(int nNodes): 
  nNodes(nNodes),
  nClusters(1),
  cluster(nNodes)
{};

Clustering::Clustering(IVector cluster) :
  nNodes(cluster.size()),
  cluster(cluster.cbegin(), cluster.cend()) 
  {
    nClusters = *std::max_element(cluster.cbegin(), cluster.cend()) + 1;  
  }


IVector Clustering::getNNodesPerCluster() const {
  IVector nNodesPerCluster(nClusters, 0);
  for(const int& clust: cluster) {
    nNodesPerCluster.at(clust)++;
  }
  return nNodesPerCluster;
}

std::vector<IVector> Clustering::getNodesPerCluster() const {
  std::vector<IVector> nodePerCluster(nClusters);
  IVector nNodesPerCluster = getNNodesPerCluster();
  for(int i =0; i < nClusters; i++)
  {
    const int cnt = nNodesPerCluster.at(i);
    nodePerCluster.at(i).reserve(cnt);
  }
  for(int i=0; i< nNodes; i++) {
    nodePerCluster.at(cluster.at(i)).push_back(i);
  }
  return nodePerCluster;
}

void Clustering::setCluster(int node, int cluster) {
  this->cluster.at(node) = cluster;
  nClusters = std::max(nClusters, cluster+1);
}

void Clustering::initSingletonClusters() {
  for(int i=0; i < nNodes; i++) {
    cluster.at(i) = i;
  }
  nClusters = nNodes;
}

void Clustering::orderClustersByNNodes() {
  typedef std::pair<int, int> ipair; // holds numNodes, cluster
  std::vector<ipair> clusterNNodes;
  clusterNNodes.reserve(nClusters);
  IVector nNodesPerCluster = getNNodesPerCluster();
  for(int i=0; i<nClusters; i++) {
    clusterNNodes.push_back(std::make_pair(nNodesPerCluster.at(i), i));
  }
  // Note order is descending
  stable_sort(clusterNNodes.begin(), clusterNNodes.end(), 
       [](const std::pair<int, int>&a, const std::pair<int, int>& b) {
         return b.first < a.first;
       });
       //std::greater<ipair>());
  
  // now make a map from old to new names
  IVector newCluster(nClusters, 0);
  int i=0;
  do {
    newCluster[clusterNNodes[i].second] = i;
    i++;
  } while (i < nClusters && clusterNNodes[i].first > 0);
  nClusters = i;
  for(int i=0; i<nNodes; i++) {
    cluster[i] = newCluster[cluster[i]];
  }
}

void Clustering::mergeClusters(const Clustering& clustering) {
  for (int i = 0; i < nNodes; i++)
    cluster.at(i) = clustering.cluster.at(cluster.at(i));
  nClusters = clustering.nClusters;
}


Network::Network() {};

Network::Network(int nNodes, DVector* nodeWeight, IVector& firstNeighborIndex, IVector& neighbor, DVector* edgeWeight) :
  nNodes(nNodes),
  nEdges(neighbor.size()),
  nodeWeight(nNodes),
  firstNeighborIndex(firstNeighborIndex.cbegin(), firstNeighborIndex.cend()),
  neighbor(neighbor.cbegin(), neighbor.cend()),
  edgeWeight(nEdges, 1.0),
  totalEdgeWeightSelfLinks(0)
  {
  
  if (edgeWeight != nullptr)
    std::copy(edgeWeight->cbegin(), edgeWeight->cend(), this->edgeWeight.begin());

  if (nodeWeight != nullptr) {
    std::copy(nodeWeight->cbegin(), nodeWeight->cend(), this->nodeWeight.begin());
  } else {
    this->nodeWeight = getTotalEdgeWeightPerNode();
  }
}


Network::Network(int nNodes, DVector* nodeWeight, std::vector<IVector>& edge, DVector* edgeWeight) :
  nNodes(nNodes),
  nEdges(0),
  nodeWeight(),
  firstNeighborIndex(nNodes + 1, 0),
  neighbor(),
  edgeWeight(),
  totalEdgeWeightSelfLinks(0)
{
  if(edge.size() != 2 || edge[0].size() != edge[1].size()) {
    throw std::length_error("Edge was supposed to be an array with 2 columns of equal size.");
  }
  IVector neighbor(edge.at(0).size(), 0);
  DVector edgeWeight2(edge.at(0).size(), 0.0);
  
  int i = 1;
  for (size_t j = 0; j < edge[0].size(); j++)
    if (edge[0][j] != edge[1][j])
    {
      if (edge[0][j] >= i)
        for (; i <= edge[0][j]; i++)
          firstNeighborIndex.at(i) = nEdges;
      neighbor[nEdges] = edge[1][j];
      edgeWeight2[nEdges] = (edgeWeight != nullptr) ? (*edgeWeight)[j] : 1.0;
      nEdges++;
    }
    else
      totalEdgeWeightSelfLinks += (edgeWeight != nullptr) ? (*edgeWeight)[j] : 1;
    for (; i <= nNodes; i++)
      firstNeighborIndex.at(i) = nEdges;
    
    this->neighbor.resize(nEdges);
    std::copy(neighbor.begin(), neighbor.begin() + nEdges, this->neighbor.begin());
    this->edgeWeight.resize(nEdges);
    std::copy(edgeWeight2.begin(), edgeWeight2.begin() + nEdges, this->edgeWeight.begin());
    if(nodeWeight == nullptr) {
      this->nodeWeight = getTotalEdgeWeightPerNode();
    } else {
      this->nodeWeight = *nodeWeight;
    }
}

double Network::getTotalNodeWeight() {
  return std::accumulate(nodeWeight.cbegin(), nodeWeight.cend(), 0.0);
}
DVector Network::getNodeWeights() {
  return nodeWeight;
}
IVector Network::getNEdgesPerNode() {
  IVector nEdgesPerNode(nNodes, 0);
  for(int i=0; i< nNodes; i++) {
    nEdgesPerNode.at(i) = firstNeighborIndex.at(i + 1) - firstNeighborIndex.at(i);
  }
  return nEdgesPerNode;
}
std::vector<IVector> Network::getEdges() {
  std::vector<IVector> edge(2);
  edge[0].resize(nEdges);
  for(int i=0; i < nNodes; i++) {
    std::fill(edge[0].begin() + firstNeighborIndex.at(i), edge[0].begin() + firstNeighborIndex.at(i + 1), i);
  }
  edge.at(1) = neighbor;
  return edge;
}
IVector Network::getEdges(int node) {
  return IVector(neighbor.begin() + firstNeighborIndex.at(node),
                 neighbor.begin() + firstNeighborIndex.at(node + 1));
}

std::vector<IVector> Network::getEdgesPerNode() {
  std::vector<IVector> edgePerNode(nNodes);
  for (int i = 0; i < nNodes; i++) {
    edgePerNode[i] = IVector(neighbor.begin() + firstNeighborIndex.at(i),
                             neighbor.begin() + firstNeighborIndex.at(i + 1));
  }
  return edgePerNode;
}

double Network::getTotalEdgeWeight() {
  return std::accumulate(edgeWeight.cbegin(), edgeWeight.cend(), 0.0) / 2.0;
}

double Network::getTotalEdgeWeight(int node) {
  return std::accumulate(edgeWeight.cbegin() + firstNeighborIndex.at(node),
                         edgeWeight.cbegin() + firstNeighborIndex.at(node + 1),
                         0.0);
}

DVector Network::getTotalEdgeWeightPerNode() {
  DVector totalEdgeWeightPerNode(nNodes, 0.0);
  for (int i = 0; i < nNodes; i++) {
    totalEdgeWeightPerNode[i] = getTotalEdgeWeight(i);
  }
  return totalEdgeWeightPerNode;
}

DVector Network::getEdgeWeights(int node) {
  return DVector(edgeWeight.cbegin() + firstNeighborIndex.at(node),
                 edgeWeight.cbegin() + firstNeighborIndex.at(node+1));
}

std::vector<DVector> Network::getEdgeWeightsPerNode() {
  std::vector<DVector> edgeWeightPerNode(nNodes);
  for (int i = 0; i < nNodes; i++)
    edgeWeightPerNode[i] = getEdgeWeights(i);
  return edgeWeightPerNode;
}


// Skipping unused Network creators
// Network createNetworkWithoutNodeWeights()
// Network createNetworkWithoutEdgeWeights()
// Network createNetworkWithoutNodeAndEdgeWeights()
// Network createNormalizedNetwork1()
// Network createNormalizedNetwork2()
// Network createPrunedNetwork(int nEdges)
// Network createPrunedNetwork(int nEdges, Random random)
// Network createSubnetwork(int[] node)
// Network createSubnetwork(boolean[] nodeInSubnetwork)
// Network createSubnetwork(Clustering clustering, int cluster)
std::vector<Network> Network::createSubnetworks(Clustering clustering) const {
  std::vector<Network> subnetwork(clustering.nClusters);
  auto nodePerCluster = clustering.getNodesPerCluster();
  IVector subnetworkNode(nNodes);
  IVector subnetworkNeighbor(nEdges);
  DVector subnetworkEdgeWeight(nEdges);
  for (int i = 0; i < clustering.nClusters; i++)
    subnetwork[i] = createSubnetwork(clustering, i, nodePerCluster[i], subnetworkNode, subnetworkNeighbor, subnetworkEdgeWeight);
  return subnetwork;
}
// Network createSubnetworkLargestComponent()
// Network createReducedNetwork(Clustering clustering)
Network Network::createReducedNetwork(const Clustering& clustering) const {
  Network reducedNetwork;
  reducedNetwork.nNodes = clustering.nClusters;
  
  reducedNetwork.nEdges = 0;
  reducedNetwork.nodeWeight = DVector(clustering.nClusters);
  reducedNetwork.firstNeighborIndex = IVector(clustering.nClusters + 1);
  reducedNetwork.totalEdgeWeightSelfLinks = totalEdgeWeightSelfLinks;
  IVector reducedNetworkNeighbor1(nEdges);
  DVector reducedNetworkEdgeWeight1(nEdges);
  IVector reducedNetworkNeighbor2(clustering.nClusters - 1);
  DVector reducedNetworkEdgeWeight2(clustering.nClusters);
  std::vector<IVector> nodePerCluster = clustering.getNodesPerCluster();
  for (int i = 0; i < clustering.nClusters; i++)
  {
    int j = 0;
    for (size_t k = 0; k < nodePerCluster[i].size(); k++)
    {
      int l = nodePerCluster[i][k];
      
      reducedNetwork.nodeWeight[i] += nodeWeight[l];
      
      for (int m = firstNeighborIndex[l]; m < firstNeighborIndex[l + 1]; m++)
      {
        int n = clustering.cluster[neighbor[m]];
        if (n != i)
        {
          if (reducedNetworkEdgeWeight2[n] == 0)
          {
            reducedNetworkNeighbor2[j] = n;
            j++;
          }
          reducedNetworkEdgeWeight2[n] += edgeWeight[m];
        }
        else
          reducedNetwork.totalEdgeWeightSelfLinks += edgeWeight[m];
      }
    }
    
    for (int k = 0; k < j; k++)
    {
      reducedNetworkNeighbor1[reducedNetwork.nEdges + k] = reducedNetworkNeighbor2[k];
      reducedNetworkEdgeWeight1[reducedNetwork.nEdges + k] = reducedNetworkEdgeWeight2[reducedNetworkNeighbor2[k]];
      reducedNetworkEdgeWeight2[reducedNetworkNeighbor2[k]] = 0;
    }
    reducedNetwork.nEdges += j;
    reducedNetwork.firstNeighborIndex[i + 1] = reducedNetwork.nEdges;
  }
  reducedNetwork.neighbor = IVector(reducedNetworkNeighbor1.cbegin(), reducedNetworkNeighbor1.cbegin() + reducedNetwork.nEdges);
  reducedNetwork.edgeWeight = DVector(reducedNetworkEdgeWeight1.cbegin(), reducedNetworkEdgeWeight1.cbegin() + reducedNetwork.nEdges);
  return reducedNetwork;
}
Clustering Network::identifyComponents() {
  std::vector<bool> nodeVisited(nNodes, false);
  IVector node(nNodes);
  Clustering clustering(nNodes);
  clustering.nClusters = 0;
  for (int i = 0; i < nNodes; i++)
    if (!nodeVisited[i])
    {
      clustering.cluster[i] = clustering.nClusters;
      nodeVisited[i] = true;
      node[0] = i;
      int j = 1;
      int k = 0;
      do
      {
        for (int l = firstNeighborIndex[node[k]]; l < firstNeighborIndex[node[k] + 1]; l++)
          if (!nodeVisited[neighbor[l]])
          {
            clustering.cluster[neighbor[l]] = clustering.nClusters;
            nodeVisited[neighbor[l]] = true;
            node[j] = neighbor[l];
            j++;
          }
          k++;
      } while (k < j);
      clustering.nClusters++;
    }
    clustering.orderClustersByNNodes();
    return clustering;
}
// private:
  // double generateRandomNumber(int node1, int node2, const IVector& nodePermutation);
Network Network::createSubnetwork(const Clustering& clustering, int cluster, IVector& node, IVector& subnetworkNode, 
                                  IVector& subnetworkNeighbor, DVector& subnetworkEdgeWeight) const {
  Network subnetwork;
  subnetwork.nNodes = node.size();
  
  if (subnetwork.nNodes == 1)
  {
    subnetwork.nEdges = 0;
    subnetwork.nodeWeight = DVector(1, nodeWeight[node[0]]);
    subnetwork.firstNeighborIndex = IVector(2);
    subnetwork.neighbor = IVector(0);
    subnetwork.edgeWeight = DVector(0);
  }
  else
  {
    for (size_t i = 0; i < node.size(); i++)
      subnetworkNode[node[i]] = i;
    
    subnetwork.nEdges = 0;
    subnetwork.nodeWeight = DVector(subnetwork.nNodes, 0);
    subnetwork.firstNeighborIndex = IVector(subnetwork.nNodes + 1);
    for (int i = 0; i < subnetwork.nNodes; i++)
    {
      int j = node[i];
      subnetwork.nodeWeight[i] = nodeWeight[j];
      for (int k = firstNeighborIndex[j]; k < firstNeighborIndex[j + 1]; k++)
        if (clustering.cluster[neighbor[k]] == cluster)
        {
          subnetworkNeighbor[subnetwork.nEdges] = subnetworkNode[neighbor[k]];
          subnetworkEdgeWeight[subnetwork.nEdges] = edgeWeight[k];
          subnetwork.nEdges++;
        }
        subnetwork.firstNeighborIndex[i + 1] = subnetwork.nEdges;
    }
    subnetwork.neighbor = IVector(subnetworkNeighbor.cbegin(), subnetworkNeighbor.cbegin() + subnetwork.nEdges);
    subnetwork.edgeWeight = DVector(subnetworkEdgeWeight.cbegin(), subnetworkEdgeWeight.cbegin() + subnetwork.nEdges);
  }
  
  subnetwork.totalEdgeWeightSelfLinks = 0;
  
  return subnetwork;
}

VOSClusteringTechnique::VOSClusteringTechnique(std::shared_ptr<Network> network, double resolution) :
  network(network),
  resolution(resolution)
  {
    clustering = std::make_shared<Clustering>(network->getNNodes());
    clustering->initSingletonClusters();
  };

VOSClusteringTechnique::VOSClusteringTechnique(std::shared_ptr<Network> network, std::shared_ptr<Clustering> clustering, double resolution) :
  network(network),
  clustering(clustering),
  resolution(resolution){};

double VOSClusteringTechnique::calcQualityFunction() {
  double qualityFunction = 0.0;
  for (int i = 0; i < network->getNNodes(); i++)
  {
    int j = clustering->cluster[i];
    for (int k = network->getFirstNeighborIndexValue(i); k < network->getFirstNeighborIndexValue(i + 1); k++)
      if (clustering->cluster[network->getNeighborValue(k)] == j)
        qualityFunction += network->edgeWeight[k];
  }
  qualityFunction += network->totalEdgeWeightSelfLinks;
  
  DVector clusterWeight(clustering->nClusters);
  for (int i = 0; i < network->nNodes; i++)
    clusterWeight[clustering->cluster[i]] += network->nodeWeight[i];
  for (int i = 0; i < clustering->nClusters; i++)
    qualityFunction -= clusterWeight[i] * clusterWeight[i] * resolution;
  
  qualityFunction /= 2 * network->getTotalEdgeWeight() + network->totalEdgeWeightSelfLinks;
  
  return qualityFunction;
}

bool VOSClusteringTechnique::runLocalMovingAlgorithm(JavaRandom& random){
  bool update = false;
  double maxQualityFunction, qualityFunction;
  DVector clusterWeight(network->getNNodes(), 0); 
  IVector nNodesPerCluster(network->getNNodes(), 0);
  
  int bestCluster, j, k, l, nNeighboringClusters, nStableNodes;
  if (network->getNNodes() == 1)
    return false;
  
  for (int i = 0; i < network->getNNodes(); i++)
  {
    clusterWeight[clustering->cluster[i]] += network->nodeWeight[i];
    nNodesPerCluster[clustering->cluster[i]]++;
  }
  
  int nUnusedClusters = 0;
  IVector unusedCluster(network->getNNodes(), 0);
  for (int i = 0; i < network->getNNodes(); i++) {
    if (nNodesPerCluster[i] == 0)
    {
      unusedCluster[nUnusedClusters] = i;
      nUnusedClusters++;
    }
  }
    
  IVector nodePermutation = Arrays2::generateRandomPermutation(network->nNodes, random);
  DVector edgeWeightPerCluster(network->getNNodes(), 0.0);
  IVector neighboringCluster(network->getNNodes() - 1, 0);
  nStableNodes = 0;
  int i = 0;
  do {
    j = nodePermutation[i];
    nNeighboringClusters = 0;
    for (k = network->firstNeighborIndex.at(j); k < network->firstNeighborIndex.at(j + 1); k++)
      {
        l = clustering->cluster[network->neighbor[k]];
        if (edgeWeightPerCluster[l] == 0)
        {
          neighboringCluster[nNeighboringClusters] = l;
          nNeighboringClusters++;
        }
        edgeWeightPerCluster[l] += network->edgeWeight[k];
      }
      
      clusterWeight[clustering->cluster[j]] -= network->nodeWeight[j];
      nNodesPerCluster[clustering->cluster[j]]--;
      if (nNodesPerCluster[clustering->cluster[j]] == 0)
      {
        unusedCluster[nUnusedClusters] = clustering->cluster[j];
        nUnusedClusters++;
      }
      
      bestCluster = -1;
      maxQualityFunction = 0;
      for (k = 0; k < nNeighboringClusters; k++)
      {
        l = neighboringCluster[k];
        qualityFunction = edgeWeightPerCluster[l] - network->nodeWeight[j] * clusterWeight[l] * resolution;
        if ((qualityFunction > maxQualityFunction) || ((qualityFunction == maxQualityFunction) && (l < bestCluster)))
        {
          bestCluster = l;
          maxQualityFunction = qualityFunction;
        }
        edgeWeightPerCluster[l] = 0;
      }
      if (maxQualityFunction == 0)
      {
        bestCluster = unusedCluster[nUnusedClusters - 1];
        nUnusedClusters--;
      }
      
      clusterWeight[bestCluster] += network->nodeWeight[j];
      nNodesPerCluster[bestCluster]++;
      if (bestCluster == clustering->cluster[j])
        nStableNodes++;
      else
      {
        clustering->cluster[j] = bestCluster;
        nStableNodes = 1;
        update = true;
      }
      
      i = (i < network->nNodes - 1) ? (i + 1) : 0;
    }
    while (nStableNodes < network->nNodes);
    
    IVector newCluster(network->getNNodes());
    clustering->nClusters = 0;
    for (i = 0; i < network->nNodes; i++)
      if (nNodesPerCluster[i] > 0)
      {
        newCluster[i] = clustering->nClusters;
        clustering->nClusters++;
      }
    for (i = 0; i < network->nNodes; i++)
        clustering->cluster[i] = newCluster[clustering->cluster[i]];
      
      return update;
}

bool VOSClusteringTechnique::runLouvainAlgorithm(JavaRandom& random) {

  if (network->nNodes == 1)
    return false;
  bool update = runLocalMovingAlgorithm(random);
  if (clustering->nClusters < network->nNodes)
  {
    VOSClusteringTechnique vosClusteringTechnique(std::make_shared<Network>(network->createReducedNetwork(*clustering)), resolution);
    
    bool update2 = vosClusteringTechnique.runLouvainAlgorithm(random);
    
    if (update2)
    {
      update = true;
      
      clustering->mergeClusters(*vosClusteringTechnique.clustering);
    }
  }
  return update;
}

bool VOSClusteringTechnique::runIteratedLouvainAlgorithm(int maxNIterations, JavaRandom& random) {
  bool update;
  int i = 0;
  do
  {
    update = runLouvainAlgorithm(random);
    i++;
  }
  while ((i < maxNIterations) && update);
  return ((i > 1) || update);
}

bool VOSClusteringTechnique::runLouvainAlgorithmWithMultilevelRefinement(JavaRandom& random) {
  if (network->nNodes == 1)
    return false;
  
  bool update = runLocalMovingAlgorithm(random);
  
  if (clustering->nClusters < network->nNodes)
  {
    VOSClusteringTechnique vosClusteringTechnique(std::make_shared<Network>(network->createReducedNetwork(*clustering)), resolution);
    
    bool update2 = vosClusteringTechnique.runLouvainAlgorithmWithMultilevelRefinement(random);
    
    if (update2)
    {
      update = true;
      clustering->mergeClusters(*vosClusteringTechnique.clustering);
      runLocalMovingAlgorithm(random);
    }
  }
  return update;}

bool VOSClusteringTechnique::runIteratedLouvainAlgorithmWithMultilevelRefinement(int maxNIterations, JavaRandom& random) {
  bool update;
  int i = 0;
  do
  {
    update = runLouvainAlgorithmWithMultilevelRefinement(random);
    i++;
  }
  while ((i < maxNIterations) && update);
  return ((i > 1) || update);
}

bool VOSClusteringTechnique::runSmartLocalMovingAlgorithm(JavaRandom& random) {
  if (network->nNodes == 1)
    return false;
  
  bool update = runLocalMovingAlgorithm(random);
  
  if (clustering->nClusters < network->nNodes)
  {
    std::vector<Network> subnetwork = network->createSubnetworks(*clustering);
    auto nodePerCluster = clustering->getNodesPerCluster();
    clustering->nClusters = 0;
    IVector nNodesPerClusterReducedNetwork(subnetwork.size());
    for (size_t i = 0; i < subnetwork.size(); i++)
    {
      VOSClusteringTechnique vosClusteringTechnique(std::make_shared<Network>(subnetwork[i]), resolution);
      vosClusteringTechnique.runLocalMovingAlgorithm(random);
      for (int j = 0; j < subnetwork[i].nNodes; j++)
        clustering->cluster[nodePerCluster[i][j]] = clustering->nClusters + vosClusteringTechnique.clustering->cluster[j];
      clustering->nClusters += vosClusteringTechnique.clustering->nClusters;
      nNodesPerClusterReducedNetwork[i] = vosClusteringTechnique.clustering->nClusters;
    }
    
    VOSClusteringTechnique vosClusteringTechnique2(std::make_shared<Network>(network->createReducedNetwork(*clustering)), resolution);
    
    int i = 0;
    for (size_t j = 0; j < nNodesPerClusterReducedNetwork.size(); j++) {
      for (int k = 0; k < nNodesPerClusterReducedNetwork[j]; k++)
      {
        vosClusteringTechnique2.clustering->cluster[i] = static_cast<int>(j);
        i++;
      }
    }
    vosClusteringTechnique2.clustering->nClusters = nNodesPerClusterReducedNetwork.size();
    
    update |= vosClusteringTechnique2.runSmartLocalMovingAlgorithm(random);
    
    clustering->mergeClusters(*vosClusteringTechnique2.clustering);
  }
  return update;
}

bool VOSClusteringTechnique::runIteratedSmartLocalMovingAlgorithm(int nIterations, JavaRandom& random) {
  bool update = false;
  for (int i = 0; i < nIterations; i++)
    update |= runSmartLocalMovingAlgorithm(random);
  return update;
}

int VOSClusteringTechnique::removeCluster(int cluster) {
  DVector clusterWeight(clustering->nClusters);
  DVector totalEdgeWeightPerCluster(clustering->nClusters);
  for (int i = 0; i < network->nNodes; i++)
  {
    clusterWeight[clustering->cluster[i]] += network->nodeWeight[i];
    if (clustering->cluster[i] == cluster)
      for (int j = network->firstNeighborIndex[i]; j < network->firstNeighborIndex[i + 1]; j++)
        totalEdgeWeightPerCluster[clustering->cluster[network->neighbor[j]]] += network->edgeWeight[j];
  }
  
  int i = -1;
  double maxQualityFunction = 0;
  for (int j = 0; j < clustering->nClusters; j++)
    if ((j != cluster) && (clusterWeight[j] > 0))
    {
      double qualityFunction = totalEdgeWeightPerCluster[j] / clusterWeight[j];
      if (qualityFunction > maxQualityFunction)
      {
        i = j;
        maxQualityFunction = qualityFunction;
      }
    }
    
    if (i >= 0)
    {
      for (int j = 0; j < network->nNodes; j++)
        if (clustering->cluster[j] == cluster)
          clustering->cluster[j] = i;
        if (cluster == clustering->nClusters - 1)
          clustering->nClusters = *std::max_element(clustering->cluster.cbegin(), clustering->cluster.cend()) + 1;
    }
    return i;
}

void VOSClusteringTechnique::removeSmallClusters(int minNNodesPerCluster) {
  VOSClusteringTechnique vosClusteringTechnique(std::make_shared<Network>(network->createReducedNetwork(*clustering)), resolution);
  auto nNodesPerCluster = clustering->getNNodesPerCluster();
  int i;
  do
  {
    i = -1;
    int j = minNNodesPerCluster;
    for (int k = 0; k < vosClusteringTechnique.clustering->nClusters; k++)
      if ((nNodesPerCluster[k] > 0) && (nNodesPerCluster[k] < j))
      {
        i = k;
        j = nNodesPerCluster[k];
      }
      
      if (i >= 0)
      {
        j = vosClusteringTechnique.removeCluster(i);
        if (j >= 0)
          nNodesPerCluster[j] += nNodesPerCluster[i];
        nNodesPerCluster[i] = 0;
      }
  }
  while (i >= 0);
  clustering->mergeClusters(*vosClusteringTechnique.clustering);
}

std::shared_ptr<Network> ModularityOptimizer::matrixToNetwork(IVector& node1, IVector& node2, DVector& edgeWeight1, int modularityFunction, int nNodes) {
  
  int n1_max = *std::max_element(node1.cbegin(), node1.cend());
  int n2_max = *std::max_element(node2.cbegin(), node2.cend());
  IVector nNeighbors(nNodes);
  for (size_t i = 0; i < node1.size(); i++)
    if (node1[i] < node2[i])
    {
      nNeighbors[node1[i]]++;
      nNeighbors[node2[i]]++;
    }
    
    IVector firstNeighborIndex(nNodes + 1);
    int nEdges = 0;
    for (int i = 0; i < nNodes; i++)
    {
      firstNeighborIndex[i] = nEdges;
      nEdges += nNeighbors[i];
    }
    firstNeighborIndex[nNodes] = nEdges;
    
    IVector neighbor(nEdges);
    DVector edgeWeight2(nEdges);
    std::fill(nNeighbors.begin(), nNeighbors.end(), 0);
    for (size_t i = 0; i < node1.size(); i++)
      if (node1[i] < node2[i])
      {
        int j = firstNeighborIndex[node1[i]] + nNeighbors[node1[i]];
        neighbor[j] = node2[i];
        edgeWeight2[j] = edgeWeight1[i];
        nNeighbors[node1[i]]++;
        j = firstNeighborIndex[node2[i]] + nNeighbors[node2[i]];
        neighbor[j] = node1[i];
        edgeWeight2[j] = edgeWeight1[i];
        nNeighbors[node2[i]]++;
      }
      
      if (modularityFunction == 1)
        return std::make_shared<Network>(nNodes, firstNeighborIndex, neighbor, &edgeWeight2);
      else
      {
        DVector nodeWeight(nNodes, 1.0);
        return std::make_shared<Network>(nNodes, &nodeWeight, firstNeighborIndex, neighbor, &edgeWeight2);
      }
}


std::shared_ptr<Network> ModularityOptimizer::readInputFile(std::string fname, int modularityFunction) {
  std::ifstream f;
  f.open(fname, std::ios::in);
  if(!f) {
    throw std::runtime_error("File could not be opened.");
  }
  std::string line;
  int nLines = 0;
  while(std::getline(f, line)) {
    nLines++;
  }
  f.clear();
  f.seekg(0, std::ios::beg);
  
  IVector node1(nLines);
  IVector node2(nLines);
  DVector edgeWeight1(nLines, 1.0);
  
  for (int j = 0; j < nLines; j++)
  {
    std::getline(f, line);
    auto splittedLine = split(line, '\t');
    node1[j] =std::stoi(splittedLine[0]);
    node2[j] = std::stoi(splittedLine[1]);
    if(splittedLine.size() > 2) {
      edgeWeight1[j] = std::stod(splittedLine[2]);
    }
  }
  int n1_max = *std::max_element(node1.cbegin(), node1.cend());
  int n2_max = *std::max_element(node2.cbegin(), node2.cend());
  int nNodes = std::max(n1_max, n2_max) + 1;
  return matrixToNetwork(node1, node2, edgeWeight1, modularityFunction, nNodes);
}

std::vector<std::string> ModularityOptimizer::split(const std::string& s, char delimiter)
{
  std::vector<std::string> tokens;
  std::string token;
  std::istringstream tokenStream(s);
  while (std::getline(tokenStream, token, delimiter))
  {
    tokens.push_back(token);
  }
  return tokens;
}
#ifdef STANDALONE

void writeOutputFile(std::string fname, Clustering& clustering) {
  int nNodes = clustering.getNNodes();
  clustering.orderClustersByNNodes();
  std::ofstream f(fname, std::ios::out);
  for(int i=0; i < nNodes; i++)
    f << clustering.getCluster(i) << std::endl;
  f.close();
}




template<typename T>
void input(std::string msg, T& value) {
  std::cout << msg << std::endl << std::endl;
  std::cin >> value;
}


int main(int argc, char* argv[]) {
  
  std::string msg = "Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck";
  std::vector<std::string> args;
  std::string inputFileName, outputFileName;
  
  bool printOutput, update;
  double modularity, maxModularity, resolution, resolution2;
  int algorithm, i, j, modularityFunction, nIterations, nRandomStarts;

  unsigned long long int randomSeed;
  for(int i=0; i<argc; i++) {
    args.emplace_back(std::string(argv[i]));
  }

  try {
    if (args.size() == 10)
    {
      inputFileName = args[1];
      outputFileName = args[2];
      modularityFunction = stoi(args[3]);
      resolution = stod(args[4]);
      algorithm = stoi(args[5]);
      nRandomStarts = stoi(args[6]);
      nIterations = stoi(args[7]);
      randomSeed = stoull(args[8]);
      printOutput = (stoi(args[9]) > 0);
      
      if (printOutput)
      {
        std::cout << msg << std::endl << std::endl;
        
      }
    }
    else
    {
      std::cout << msg << std::endl << std::endl;
      input<std::string>("Input file name: ", inputFileName);
      input<std::string>("Output file name: ", outputFileName);
      input<int>("Modularity function (1 = standard; 2 = alternative): ", modularityFunction);
      input<double>("Resolution parameter (e.g., 1.0): ", resolution);
      input<int>("Algorithm (1 = Louvain; 2 = Louvain with multilevel refinement; 3 = smart local moving): ", algorithm);
      input<int>("Number of random starts (e.g., 10): ", nRandomStarts);
      input<int>("Number of iterations (e.g., 10): ",nIterations);
      input<unsigned long long int>("Random seed (e.g., 0): ", randomSeed);
      int tmp;
      input<int>("Print output (0 = no; 1 = yes): ",tmp);
      printOutput = tmp > 0;
      std::cout << std::endl;
    }
    
    
    if (printOutput)
    {
      std::cout << "Reading input file..." << std::endl << std::endl;
    }
    
    std::shared_ptr<Network> network = readInputFile(inputFileName, modularityFunction);
    
    if (printOutput)
    {
      std::printf("Number of nodes: %d\n", network->getNNodes());
      std::printf("Number of edges: %d\n", network->getNEdges());
      std::cout << std::endl;
      std::cout << "Running " <<  ((algorithm == 1) ? "Louvain algorithm" : ((algorithm == 2) ? "Louvain algorithm with multilevel refinement" : "smart local moving algorithm")) << "...";
      std::cout << std::endl;
    }
    
    resolution2 = ((modularityFunction == 1) ? (resolution / (2 * network->getTotalEdgeWeight() + network->getTotalEdgeWeightSelfLinks())) : resolution);
    
    auto beginTime = duration_cast<milliseconds>(system_clock::now().time_since_epoch());
    std::shared_ptr<Clustering> clustering;
    maxModularity = -std::numeric_limits<double>::infinity();
    JavaRandom random(randomSeed);
    for (i = 0; i < nRandomStarts; i++)
    {
      if (printOutput && (nRandomStarts > 1))
        std::printf("Random start: %d\n", i + 1);
      
      VOSClusteringTechnique vosClusteringTechnique(network, resolution2);
      
      j = 0;
      update = true;
      do
      {
        if (printOutput && (nIterations > 1))
          std::printf("Iteration: %d\n", j + 1);
        
        if (algorithm == 1)
          update = vosClusteringTechnique.runLouvainAlgorithm(random);
        else if (algorithm == 2)
          update = vosClusteringTechnique.runLouvainAlgorithmWithMultilevelRefinement(random);
        else if (algorithm == 3)
          vosClusteringTechnique.runSmartLocalMovingAlgorithm(random);
        j++;
        
        modularity = vosClusteringTechnique.calcQualityFunction();
        
        if (printOutput && (nIterations > 1))
          std::printf("Modularity: %.4f\n", modularity);
      }
      while ((j < nIterations) && update);
      
      if (modularity > maxModularity)
      {
        clustering = vosClusteringTechnique.getClustering();
        maxModularity = modularity;
      }
      
      if (printOutput && (nRandomStarts > 1))
      {
        if (nIterations == 1)
          std::printf("Modularity: %.4f\n", modularity);
        std::cout << std::endl;
      }
    }
    auto endTime = duration_cast<milliseconds>(system_clock::now().time_since_epoch());
    
    if (printOutput)
    {
      if (nRandomStarts == 1)
      {
        if (nIterations > 1)
          std::cout << std::endl;
        std::printf("Modularity: %.4f\n", maxModularity);
      }
      else
      std::printf("Maximum modularity in %d random starts: %.4f\n", nRandomStarts, maxModularity);
      std::printf("Number of communities: %d\n", clustering->getNClusters());
      std::printf("Elapsed time: %d seconds\n", static_cast<int>((endTime - beginTime).count() / 1000.0));
      std::cout << std::endl << "Writing output file..." << std::endl;
    }
    
    writeOutputFile(outputFileName, *clustering);
  } catch (std::exception a) {
    std::cout << a.what() << std::endl;
  }
  return 0;
};

#endif
