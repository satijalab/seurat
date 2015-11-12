/**
 * Clustering
 *
 * @author Ludo Waltman
 * @author Nees Jan van Eck
 * @version 1.3.1 11/17/14
 */

import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.Serializable;
import java.util.Arrays;

public class Clustering implements Cloneable, Serializable
{
    private static final long serialVersionUID = 1;

    protected int nNodes;
    protected int nClusters;
    protected int[] cluster;

    public static Clustering load(String fileName) throws ClassNotFoundException, IOException
    {
        Clustering clustering;
        ObjectInputStream objectInputStream;

        objectInputStream = new ObjectInputStream(new FileInputStream(fileName));

        clustering = (Clustering)objectInputStream.readObject();

        objectInputStream.close();

        return clustering;
    }

    public Clustering(int nNodes)
    {
        this.nNodes = nNodes;
        cluster = new int[nNodes];
        nClusters = 1;
    }

    public Clustering(int[] cluster)
    {
        nNodes = cluster.length;
        this.cluster = (int[])cluster.clone();
        nClusters = Arrays2.calcMaximum(cluster) + 1;
    }

    public Object clone()
    {
        Clustering clonedClustering;

        try
        {
            clonedClustering = (Clustering)super.clone();
            clonedClustering.cluster = getClusters();
            return clonedClustering;
        }
        catch (CloneNotSupportedException e)
        {
            return null;
        }
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

    public int getNClusters()
    {
        return nClusters;
    }

    public int[] getClusters()
    {
        return (int[])cluster.clone();
    }

    public int getCluster(int node)
    {
        return cluster[node];
    }

    public int[] getNNodesPerCluster()
    {
        int i;
        int[] nNodesPerCluster;

        nNodesPerCluster = new int[nClusters];
        for (i = 0; i < nNodes; i++)
            nNodesPerCluster[cluster[i]]++;
        return nNodesPerCluster;
    }

    public int[][] getNodesPerCluster()
    {
        int i;
        int[] nNodesPerCluster;
        int[][] nodePerCluster;

        nodePerCluster = new int[nClusters][];
        nNodesPerCluster = getNNodesPerCluster();
        for (i = 0; i < nClusters; i++)
        {
            nodePerCluster[i] = new int[nNodesPerCluster[i]];
            nNodesPerCluster[i] = 0;
        }
        for (i = 0; i < nNodes; i++)
        {
            nodePerCluster[cluster[i]][nNodesPerCluster[cluster[i]]] = i;
            nNodesPerCluster[cluster[i]]++;
        }
        return nodePerCluster;
    }

    public void setCluster(int node, int cluster)
    {
        this.cluster[node] = cluster;
        nClusters = Math.max(nClusters, cluster + 1);
    }

    public void initSingletonClusters()
    {
        int i;

        for (i = 0; i < nNodes; i++)
            cluster[i] = i;
        nClusters = nNodes;
    }

    public void orderClustersByNNodes()
    {
        class ClusterNNodes implements Comparable<ClusterNNodes>
        {
            public int cluster;
            public int nNodes;

            public ClusterNNodes(int cluster, int nNodes)
            {
                this.cluster = cluster;
                this.nNodes = nNodes;
            }

            public int compareTo(ClusterNNodes clusterNNodes)
            {
                return (clusterNNodes.nNodes > nNodes) ? 1 : ((clusterNNodes.nNodes < nNodes) ? -1 : 0);
            }
        }

        ClusterNNodes[] clusterNNodes;
        int i;
        int[] newCluster, nNodesPerCluster;

        nNodesPerCluster = getNNodesPerCluster();
        clusterNNodes = new ClusterNNodes[nClusters];
        for (i = 0; i < nClusters; i++)
            clusterNNodes[i] = new ClusterNNodes(i, nNodesPerCluster[i]);

        Arrays.sort(clusterNNodes);

        newCluster = new int[nClusters];
        i = 0;
        do
        {
            newCluster[clusterNNodes[i].cluster] = i;
            i++;
        }
        while ((i < nClusters) && (clusterNNodes[i].nNodes > 0));
        nClusters = i;
        for (i = 0; i < nNodes; i++)
            cluster[i] = newCluster[cluster[i]];
    }

    public void mergeClusters(Clustering clustering)
    {
        int i;

        for (i = 0; i < nNodes; i++)
            cluster[i] = clustering.cluster[cluster[i]];
        nClusters = clustering.nClusters;
    }
}
