/**
 * Network
 *
 * @author Ludo Waltman
 * @author Nees Jan van Eck
 * @version 1.3.1, 08/30/15
 */

import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.Serializable;
import java.util.Arrays;
import java.util.Random;

public class Network implements Serializable
{
    private static final long serialVersionUID = 1;

    protected int nNodes;
    protected int nEdges;
    protected double[] nodeWeight;
    protected int[] firstNeighborIndex;
    protected int[] neighbor;
    protected double[] edgeWeight;
    protected double totalEdgeWeightSelfLinks;

    public static Network load(String fileName) throws ClassNotFoundException, IOException
    {
        Network network;
        ObjectInputStream objectInputStream;

        objectInputStream = new ObjectInputStream(new FileInputStream(fileName));

        network = (Network)objectInputStream.readObject();

        objectInputStream.close();

        return network;
    }

    public Network(int nNodes, int[][] edge)
    {
        this(nNodes, null, edge, null);
    }

    public Network(int nNodes, double[] nodeWeight, int[][] edge)
    {
        this(nNodes, nodeWeight, edge, null);
    }

    public Network(int nNodes, int[][] edge, double[] edgeWeight)
    {
        this(nNodes, null, edge, edgeWeight);
    }

    public Network(int nNodes, double[] nodeWeight, int[][] edge, double[] edgeWeight)
    {
        double[] edgeWeight2;
        int i, j;
        int[] neighbor;

        this.nNodes = nNodes;

        nEdges = 0;
        firstNeighborIndex = new int[nNodes + 1];
        neighbor = new int[edge[0].length];
        edgeWeight2 = new double[edge[0].length];
        totalEdgeWeightSelfLinks = 0;
        i = 1;
        for (j = 0; j < edge[0].length; j++)
            if (edge[0][j] != edge[1][j])
            {
                if (edge[0][j] >= i)
                    for (; i <= edge[0][j]; i++)
                        firstNeighborIndex[i] = nEdges;
                neighbor[nEdges] = edge[1][j];
                edgeWeight2[nEdges] = (edgeWeight != null) ? edgeWeight[j] : 1;
                nEdges++;
            }
            else
                totalEdgeWeightSelfLinks += (edgeWeight != null) ? edgeWeight[j] : 1;
        for (; i <= nNodes; i++)
            firstNeighborIndex[i] = nEdges;
        this.neighbor = Arrays.copyOfRange(neighbor, 0, nEdges);
        this.edgeWeight = Arrays.copyOfRange(edgeWeight2, 0, nEdges);

        this.nodeWeight = (nodeWeight != null) ? (double[])nodeWeight.clone() : getTotalEdgeWeightPerNode();
    }

    public Network(int nNodes, int[] firstNeighborIndex, int[] neighbor)
    {
        this(nNodes, null, firstNeighborIndex, neighbor, null);
    }

    public Network(int nNodes, double[] nodeWeight, int[] firstNeighborIndex, int[] neighbor)
    {
        this(nNodes, nodeWeight, firstNeighborIndex, neighbor, null);
    }

    public Network(int nNodes, int[] firstNeighborIndex, int[] neighbor, double[] edgeWeight)
    {
        this(nNodes, null, firstNeighborIndex, neighbor, edgeWeight);
    }

    public Network(int nNodes, double[] nodeWeight, int[] firstNeighborIndex, int[] neighbor, double[] edgeWeight)
    {
        this.nNodes = nNodes;

        nEdges = neighbor.length;
        this.firstNeighborIndex = (int[])firstNeighborIndex.clone();
        this.neighbor = (int[])neighbor.clone();
        if (edgeWeight != null)
            this.edgeWeight = (double[])edgeWeight.clone();
        else
        {
            this.edgeWeight = new double[nEdges];
            Arrays.fill(this.edgeWeight, 1);
        }
        totalEdgeWeightSelfLinks = 0;

        this.nodeWeight = (nodeWeight != null) ? (double[])nodeWeight.clone() : getTotalEdgeWeightPerNode();
    }

    public void save(String fileName) throws IOException
    {
        ObjectOutputStream objectOutputStream;

        objectOutputStream = new ObjectOutputStream(new FileOutputStream(fileName));

        objectOutputStream.writeObject(this);

        objectOutputStream.close();
    }

    public int getNNodes()
    {
        return nNodes;
    }

    public double getTotalNodeWeight()
    {
        return Arrays2.calcSum(nodeWeight);
    }

    public double[] getNodeWeights()
    {
        return (double[])nodeWeight.clone();
    }

    public double getNodeWeight(int node)
    {
        return nodeWeight[node];
    }

    public int getNEdges()
    {
        return nEdges / 2;
    }

    public int getNEdges(int node)
    {
        return firstNeighborIndex[node + 1] - firstNeighborIndex[node];
    }

    public int[] getNEdgesPerNode()
    {
        int i;
        int[] nEdgesPerNode;

        nEdgesPerNode = new int[nNodes];
        for (i = 0; i < nNodes; i++)
            nEdgesPerNode[i] = firstNeighborIndex[i + 1] - firstNeighborIndex[i];
        return nEdgesPerNode;
    }

    public int[][] getEdges()
    {
        int i;
        int[][] edge;

        edge = new int[2][];
        edge[0] = new int[nEdges];
        for (i = 0; i < nNodes; i++)
            Arrays.fill(edge[0], firstNeighborIndex[i], firstNeighborIndex[i + 1], i);
        edge[1] = (int[])neighbor.clone();
        return edge;
    }

    public int[] getEdges(int node)
    {
        return Arrays.copyOfRange(neighbor, firstNeighborIndex[node], firstNeighborIndex[node + 1]);
    }

    public int[][] getEdgesPerNode()
    {
        int i;
        int[][] edgePerNode;

        edgePerNode = new int[nNodes][];
        for (i = 0; i < nNodes; i++)
            edgePerNode[i] = Arrays.copyOfRange(neighbor, firstNeighborIndex[i], firstNeighborIndex[i + 1]);
        return edgePerNode;
    }

    public double getTotalEdgeWeight()
    {
        return Arrays2.calcSum(edgeWeight) / 2;
    }

    public double getTotalEdgeWeight(int node)
    {
        return Arrays2.calcSum(edgeWeight, firstNeighborIndex[node], firstNeighborIndex[node + 1]);
    }

    public double[] getTotalEdgeWeightPerNode()
    {
        double[] totalEdgeWeightPerNode;
        int i;

        totalEdgeWeightPerNode = new double[nNodes];
        for (i = 0; i < nNodes; i++)
            totalEdgeWeightPerNode[i] = Arrays2.calcSum(edgeWeight, firstNeighborIndex[i], firstNeighborIndex[i + 1]);
        return totalEdgeWeightPerNode;
    }

    public double[] getEdgeWeights()
    {
        return (double[])edgeWeight.clone();
    }

    public double[] getEdgeWeights(int node)
    {
        return Arrays.copyOfRange(edgeWeight, firstNeighborIndex[node], firstNeighborIndex[node + 1]);
    }

    public double[][] getEdgeWeightsPerNode()
    {
        double[][] edgeWeightPerNode;
        int i;

        edgeWeightPerNode = new double[nNodes][];
        for (i = 0; i < nNodes; i++)
            edgeWeightPerNode[i] = Arrays.copyOfRange(edgeWeight, firstNeighborIndex[i], firstNeighborIndex[i + 1]);
        return edgeWeightPerNode;
    }

    public double getTotalEdgeWeightSelfLinks()
    {
        return totalEdgeWeightSelfLinks;
    }

    public Network createNetworkWithoutNodeWeights()
    {
        Network networkWithoutNodeWeights;

        networkWithoutNodeWeights = new Network();
        networkWithoutNodeWeights.nNodes = nNodes;
        networkWithoutNodeWeights.nEdges = nEdges;
        networkWithoutNodeWeights.nodeWeight = new double[nNodes];
        Arrays.fill(networkWithoutNodeWeights.nodeWeight, 1);
        networkWithoutNodeWeights.firstNeighborIndex = firstNeighborIndex;
        networkWithoutNodeWeights.neighbor = neighbor;
        networkWithoutNodeWeights.edgeWeight = edgeWeight;
        networkWithoutNodeWeights.totalEdgeWeightSelfLinks = totalEdgeWeightSelfLinks;
        return networkWithoutNodeWeights;
    }

    public Network createNetworkWithoutEdgeWeights()
    {
        Network networkWithoutEdgeWeights;

        networkWithoutEdgeWeights = new Network();
        networkWithoutEdgeWeights.nNodes = nNodes;
        networkWithoutEdgeWeights.nEdges = nEdges;
        networkWithoutEdgeWeights.nodeWeight = nodeWeight;
        networkWithoutEdgeWeights.firstNeighborIndex = firstNeighborIndex;
        networkWithoutEdgeWeights.neighbor = neighbor;
        networkWithoutEdgeWeights.edgeWeight = new double[nEdges];
        Arrays.fill(networkWithoutEdgeWeights.edgeWeight, 1);
        networkWithoutEdgeWeights.totalEdgeWeightSelfLinks = 0;
        return networkWithoutEdgeWeights;
    }

    public Network createNetworkWithoutNodeAndEdgeWeights()
    {
        Network networkWithoutNodeAndEdgeWeights;

        networkWithoutNodeAndEdgeWeights = new Network();
        networkWithoutNodeAndEdgeWeights.nNodes = nNodes;
        networkWithoutNodeAndEdgeWeights.nEdges = nEdges;
        networkWithoutNodeAndEdgeWeights.nodeWeight = new double[nNodes];
        Arrays.fill(networkWithoutNodeAndEdgeWeights.nodeWeight, 1);
        networkWithoutNodeAndEdgeWeights.firstNeighborIndex = firstNeighborIndex;
        networkWithoutNodeAndEdgeWeights.neighbor = neighbor;
        networkWithoutNodeAndEdgeWeights.edgeWeight = new double[nEdges];
        Arrays.fill(networkWithoutNodeAndEdgeWeights.edgeWeight, 1);
        networkWithoutNodeAndEdgeWeights.totalEdgeWeightSelfLinks = 0;
        return networkWithoutNodeAndEdgeWeights;
    }

    public Network createNormalizedNetwork1()
    {
        double totalNodeWeight;
        int i, j;
        Network normalizedNetwork;

        normalizedNetwork = new Network();

        normalizedNetwork.nNodes = nNodes;
        normalizedNetwork.nEdges = nEdges;
        normalizedNetwork.nodeWeight = new double[nNodes];
        Arrays.fill(normalizedNetwork.nodeWeight, 1);
        normalizedNetwork.firstNeighborIndex = firstNeighborIndex;
        normalizedNetwork.neighbor = neighbor;

        normalizedNetwork.edgeWeight = new double[nEdges];
        totalNodeWeight = getTotalNodeWeight();
        for (i = 0; i < nNodes; i++)
            for (j = firstNeighborIndex[i]; j < firstNeighborIndex[i + 1]; j++)
                normalizedNetwork.edgeWeight[j] = edgeWeight[j] / ((nodeWeight[i] * nodeWeight[neighbor[j]]) / totalNodeWeight);

        normalizedNetwork.totalEdgeWeightSelfLinks = 0;

        return normalizedNetwork;
    }

    public Network createNormalizedNetwork2()
    {
        int i, j;
        Network normalizedNetwork;

        normalizedNetwork = new Network();

        normalizedNetwork.nNodes = nNodes;
        normalizedNetwork.nEdges = nEdges;
        normalizedNetwork.nodeWeight = new double[nNodes];
        Arrays.fill(normalizedNetwork.nodeWeight, 1);
        normalizedNetwork.firstNeighborIndex = firstNeighborIndex;
        normalizedNetwork.neighbor = neighbor;

        normalizedNetwork.edgeWeight = new double[nEdges];
        for (i = 0; i < nNodes; i++)
            for (j = firstNeighborIndex[i]; j < firstNeighborIndex[i + 1]; j++)
                normalizedNetwork.edgeWeight[j] = edgeWeight[j] / (2 / (nNodes / nodeWeight[i] + nNodes / nodeWeight[neighbor[j]]));

        normalizedNetwork.totalEdgeWeightSelfLinks = 0;

        return normalizedNetwork;
    }

    public Network createPrunedNetwork(int nEdges)
    {
        return createPrunedNetwork(nEdges, new Random());
    }

    public Network createPrunedNetwork(int nEdges, Random random)
    {
        double edgeWeightThreshold, randomNumberThreshold;
        double[] edgeWeight, randomNumber;
        int i, j, k, nEdgesAboveThreshold, nEdgesAtThreshold;
        int[] nodePermutation;
        Network prunedNetwork;

        nEdges *= 2;

        if (nEdges >= this.nEdges)
            return this;

        edgeWeight = new double[this.nEdges / 2];
        i = 0;
        for (j = 0; j < nNodes; j++)
            for (k = firstNeighborIndex[j]; k < firstNeighborIndex[j + 1]; k++)
                if (neighbor[k] < j)
                {
                    edgeWeight[i] = this.edgeWeight[k];
                    i++;
                }
        Arrays.sort(edgeWeight);
        edgeWeightThreshold = edgeWeight[(this.nEdges - nEdges) / 2];

        nEdgesAboveThreshold = 0;
        while (edgeWeight[this.nEdges / 2 - nEdgesAboveThreshold - 1] > edgeWeightThreshold)
            nEdgesAboveThreshold++;
        nEdgesAtThreshold = 0;
        while ((nEdgesAboveThreshold + nEdgesAtThreshold < this.nEdges / 2) && (edgeWeight[this.nEdges / 2 - nEdgesAboveThreshold - nEdgesAtThreshold - 1] == edgeWeightThreshold))
            nEdgesAtThreshold++;

        nodePermutation = Arrays2.generateRandomPermutation(nNodes, random);

        randomNumber = new double[nEdgesAtThreshold];
        i = 0;
        for (j = 0; j < nNodes; j++)
            for (k = firstNeighborIndex[j]; k < firstNeighborIndex[j + 1]; k++)
                if ((neighbor[k] < j) && (this.edgeWeight[k] == edgeWeightThreshold))
                {
                    randomNumber[i] = generateRandomNumber(j, neighbor[k], nodePermutation);
                    i++;
                }
        Arrays.sort(randomNumber);
        randomNumberThreshold = randomNumber[nEdgesAboveThreshold + nEdgesAtThreshold - nEdges / 2];

        prunedNetwork = new Network();

        prunedNetwork.nNodes = nNodes;
        prunedNetwork.nEdges = nEdges;
        prunedNetwork.nodeWeight = nodeWeight;

        prunedNetwork.firstNeighborIndex = new int[nNodes + 1];
        prunedNetwork.neighbor = new int[nEdges];
        prunedNetwork.edgeWeight = new double[nEdges];
        i = 0;
        for (j = 0; j < nNodes; j++)
        {
            for (k = firstNeighborIndex[j]; k < firstNeighborIndex[j + 1]; k++)
                if ((this.edgeWeight[k] > edgeWeightThreshold) || ((this.edgeWeight[k] == edgeWeightThreshold) && (generateRandomNumber(j, neighbor[k], nodePermutation) >= randomNumberThreshold)))
                {
                    prunedNetwork.neighbor[i] = neighbor[k];
                    prunedNetwork.edgeWeight[i] = this.edgeWeight[k];
                    i++;
                }
            prunedNetwork.firstNeighborIndex[j + 1] = i;
        }

        prunedNetwork.totalEdgeWeightSelfLinks = totalEdgeWeightSelfLinks;

        return prunedNetwork;
    }

    public Network createSubnetwork(int[] node)
    {
        double[] subnetworkEdgeWeight;
        int i, j, k;
        int[] subnetworkNode, subnetworkNeighbor;
        Network subnetwork;

        subnetwork = new Network();

        subnetwork.nNodes = node.length;

        if (subnetwork.nNodes == 1)
        {
            subnetwork.nEdges = 0;
            subnetwork.nodeWeight = new double[] {nodeWeight[node[0]]};
            subnetwork.firstNeighborIndex = new int[2];
            subnetwork.neighbor = new int[0];
            subnetwork.edgeWeight = new double[0];
        }
        else
        {
            subnetworkNode = new int[nNodes];
            Arrays.fill(subnetworkNode, -1);
            for (i = 0; i < node.length; i++)
                subnetworkNode[node[i]] = i;

            subnetwork.nEdges = 0;
            subnetwork.nodeWeight = new double[subnetwork.nNodes];
            subnetwork.firstNeighborIndex = new int[subnetwork.nNodes + 1];
            subnetworkNeighbor = new int[nEdges];
            subnetworkEdgeWeight = new double[nEdges];
            for (i = 0; i < subnetwork.nNodes; i++)
            {
                j = node[i];
                subnetwork.nodeWeight[i] = nodeWeight[j];
                for (k = firstNeighborIndex[j]; k < firstNeighborIndex[j + 1]; k++)
                    if (subnetworkNode[neighbor[k]] >= 0)
                    {
                        subnetworkNeighbor[subnetwork.nEdges] = subnetworkNode[neighbor[k]];
                        subnetworkEdgeWeight[subnetwork.nEdges] = edgeWeight[k];
                        subnetwork.nEdges++;
                    }
                subnetwork.firstNeighborIndex[i + 1] = subnetwork.nEdges;
            }
            subnetwork.neighbor = Arrays.copyOfRange(subnetworkNeighbor, 0, subnetwork.nEdges);
            subnetwork.edgeWeight = Arrays.copyOfRange(subnetworkEdgeWeight, 0, subnetwork.nEdges);
        }

        subnetwork.totalEdgeWeightSelfLinks = 0;

        return subnetwork;
    }

    public Network createSubnetwork(boolean[] nodeInSubnetwork)
    {
        int i, j;
        int[] node;

        i = 0;
        for (j = 0; j < nNodes; j++)
            if (nodeInSubnetwork[j])
                i++;
        node = new int[i];
        i = 0;
        for (j = 0; j < nNodes; j++)
            if (nodeInSubnetwork[j])
            {
                node[i] = j;
                i++;
            }
        return createSubnetwork(node);
    }

    public Network createSubnetwork(Clustering clustering, int cluster)
    {
        double[] subnetworkEdgeWeight;
        int[] subnetworkNeighbor, subnetworkNode;
        int[][] nodePerCluster;
        Network subnetwork;

        nodePerCluster = clustering.getNodesPerCluster();
        subnetworkNode = new int[nNodes];
        subnetworkNeighbor = new int[nEdges];
        subnetworkEdgeWeight = new double[nEdges];
        subnetwork = createSubnetwork(clustering, cluster, nodePerCluster[cluster], subnetworkNode, subnetworkNeighbor, subnetworkEdgeWeight);
        return subnetwork;
    }

    public Network[] createSubnetworks(Clustering clustering)
    {
        double[] subnetworkEdgeWeight;
        int i;
        int[] subnetworkNeighbor, subnetworkNode;
        int[][] nodePerCluster;
        Network[] subnetwork;

        subnetwork = new Network[clustering.nClusters];
        nodePerCluster = clustering.getNodesPerCluster();
        subnetworkNode = new int[nNodes];
        subnetworkNeighbor = new int[nEdges];
        subnetworkEdgeWeight = new double[nEdges];
        for (i = 0; i < clustering.nClusters; i++)
            subnetwork[i] = createSubnetwork(clustering, i, nodePerCluster[i], subnetworkNode, subnetworkNeighbor, subnetworkEdgeWeight);
        return subnetwork;
    }

    public Network createSubnetworkLargestComponent()
    {
        return createSubnetwork(identifyComponents(), 0);
    }

    public Network createReducedNetwork(Clustering clustering)
    {
        double[] reducedNetworkEdgeWeight1, reducedNetworkEdgeWeight2;
        int i, j, k, l, m, n;
        int[] reducedNetworkNeighbor1, reducedNetworkNeighbor2;
        int[][] nodePerCluster;
        Network reducedNetwork;

        reducedNetwork = new Network();

        reducedNetwork.nNodes = clustering.nClusters;

        reducedNetwork.nEdges = 0;
        reducedNetwork.nodeWeight = new double[clustering.nClusters];
        reducedNetwork.firstNeighborIndex = new int[clustering.nClusters + 1];
        reducedNetwork.totalEdgeWeightSelfLinks = totalEdgeWeightSelfLinks;
        reducedNetworkNeighbor1 = new int[nEdges];
        reducedNetworkEdgeWeight1 = new double[nEdges];
        reducedNetworkNeighbor2 = new int[clustering.nClusters - 1];
        reducedNetworkEdgeWeight2 = new double[clustering.nClusters];
        nodePerCluster = clustering.getNodesPerCluster();
        for (i = 0; i < clustering.nClusters; i++)
        {
            j = 0;
            for (k = 0; k < nodePerCluster[i].length; k++)
            {
                l = nodePerCluster[i][k];

                reducedNetwork.nodeWeight[i] += nodeWeight[l];

                for (m = firstNeighborIndex[l]; m < firstNeighborIndex[l + 1]; m++)
                {
                    n = clustering.cluster[neighbor[m]];
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

            for (k = 0; k < j; k++)
            {
                reducedNetworkNeighbor1[reducedNetwork.nEdges + k] = reducedNetworkNeighbor2[k];
                reducedNetworkEdgeWeight1[reducedNetwork.nEdges + k] = reducedNetworkEdgeWeight2[reducedNetworkNeighbor2[k]];
                reducedNetworkEdgeWeight2[reducedNetworkNeighbor2[k]] = 0;
            }
            reducedNetwork.nEdges += j;
            reducedNetwork.firstNeighborIndex[i + 1] = reducedNetwork.nEdges;
        }
        reducedNetwork.neighbor = Arrays.copyOfRange(reducedNetworkNeighbor1, 0, reducedNetwork.nEdges);
        reducedNetwork.edgeWeight = Arrays.copyOfRange(reducedNetworkEdgeWeight1, 0, reducedNetwork.nEdges);

        return reducedNetwork;
    }

    public Clustering identifyComponents()
    {
        boolean[] nodeVisited;
        Clustering clustering;
        int i, j, k, l;
        int[] node;

        clustering = new Clustering(nNodes);

        clustering.nClusters = 0;
        nodeVisited = new boolean[nNodes];
        node = new int[nNodes];
        for (i = 0; i < nNodes; i++)
            if (!nodeVisited[i])
            {
                clustering.cluster[i] = clustering.nClusters;
                nodeVisited[i] = true;
                node[0] = i;
                j = 1;
                k = 0;
                do
                {
                    for (l = firstNeighborIndex[node[k]]; l < firstNeighborIndex[node[k] + 1]; l++)
                        if (!nodeVisited[neighbor[l]])
                        {
                            clustering.cluster[neighbor[l]] = clustering.nClusters;
                            nodeVisited[neighbor[l]] = true;
                            node[j] = neighbor[l];
                            j++;
                        }
                    k++;
                }
                while (k < j);

                clustering.nClusters++;
            }

        clustering.orderClustersByNNodes();

        return clustering;
    }

    private Network()
    {
    }

    private double generateRandomNumber(int node1, int node2, int[] nodePermutation)
    {
        int i, j;
        Random random;

        if (node1 < node2)
        {
            i = node1;
            j = node2;
        }
        else
        {
            i = node2;
            j = node1;
        }
        random = new Random(nodePermutation[i] * nNodes + nodePermutation[j]);
        return random.nextDouble();
    }

    private Network createSubnetwork(Clustering clustering, int cluster, int[] node, int[] subnetworkNode, int[] subnetworkNeighbor, double[] subnetworkEdgeWeight)
    {
        int i, j, k;
        Network subnetwork;

        subnetwork = new Network();

        subnetwork.nNodes = node.length;

        if (subnetwork.nNodes == 1)
        {
            subnetwork.nEdges = 0;
            subnetwork.nodeWeight = new double[] {nodeWeight[node[0]]};
            subnetwork.firstNeighborIndex = new int[2];
            subnetwork.neighbor = new int[0];
            subnetwork.edgeWeight = new double[0];
        }
        else
        {
            for (i = 0; i < node.length; i++)
                subnetworkNode[node[i]] = i;

            subnetwork.nEdges = 0;
            subnetwork.nodeWeight = new double[subnetwork.nNodes];
            subnetwork.firstNeighborIndex = new int[subnetwork.nNodes + 1];
            for (i = 0; i < subnetwork.nNodes; i++)
            {
                j = node[i];
                subnetwork.nodeWeight[i] = nodeWeight[j];
                for (k = firstNeighborIndex[j]; k < firstNeighborIndex[j + 1]; k++)
                    if (clustering.cluster[neighbor[k]] == cluster)
                    {
                        subnetworkNeighbor[subnetwork.nEdges] = subnetworkNode[neighbor[k]];
                        subnetworkEdgeWeight[subnetwork.nEdges] = edgeWeight[k];
                        subnetwork.nEdges++;
                    }
                subnetwork.firstNeighborIndex[i + 1] = subnetwork.nEdges;
            }
            subnetwork.neighbor = Arrays.copyOfRange(subnetworkNeighbor, 0, subnetwork.nEdges);
            subnetwork.edgeWeight = Arrays.copyOfRange(subnetworkEdgeWeight, 0, subnetwork.nEdges);
        }

        subnetwork.totalEdgeWeightSelfLinks = 0;

        return subnetwork;
    }
}
