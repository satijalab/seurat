#' @include seurat.R
NULL

#' SNN Graph Construction
#'
#' Constructs a Shared Nearest Neighbor (SNN) Graph for a given dataset. We
#' first determine the k-nearest neighbors of each cell (defined by k.param *
#' k.scale). We use this knn graph to construct the SNN graph by calculating the
#' neighborhood overlap (Jaccard distance) between every cell and its k.param *
#' k.scale nearest neighbors (defining the neighborhood for each cell as the
#' k.param nearest neighbors).
#'
#' @param object Seurat object
#' @param genes.use A vector of gene names to use in construction of SNN graph
#' if building directly based on expression data rather than a dimensionally
#' reduced representation (i.e. PCs).
#' @param reduction.type Name of dimensional reduction technique to use in
#' construction of SNN graph. (e.g. "pca", "ica")
#' @param dims.use A vector of the dimensions to use in construction of the SNN
#' graph (e.g. To use the first 10 PCs, pass 1:10)
#' @param k.param Defines k for the k-nearest neighbor algorithm
#' @param k.scale Granularity option for k.param
#' @param plot.SNN Plot the SNN graph
#' @param prune.SNN Sets the cutoff for acceptable Jaccard distances when
#' computing the neighborhood overlap for the SNN construction. Any edges with
#' values less than or equal to this will be set to 0 and removed from the SNN
#' graph. Essentially sets the strigency of pruning (0 --- no pruning, 1 ---
#' prune everything).
#' @param print.output Whether or not to print output to the console
#' @param distance.matrix Build SNN from distance matrix (experimental)
#' @param force.recalc Force recalculation of SNN.
#' @param filename Write SNN directly to file named here as an edge list compatible with FindClusters
#' @param save.SNN Default behavior is to store the SNN in object@@snn. Setting to FALSE can be used
#' together with a provided filename to only write the SNN out as an edge file to disk.
#' @importFrom RANN nn2
#' @importFrom igraph plot.igraph graph.adjlist graph.adjacency E
#' @importFrom Matrix sparseMatrix
#'
#' @return Returns the object with object@@snn filled
#'
#' @rdname BuildSNN
#' @export BuildSNN
#'
#' @examples
#'
#' pbmc_small
#' # Compute an SNN on the gene expression level
#' pbmc_small <- BuildSNN(pbmc_small, genes.use = pbmc_small@var.genes)
#'
#' # More commonly, we build the SNN on a dimensionally reduced form of the data
#' # such as the first 10 principle components.
#'
#' pbmc_small <- BuildSNN(pbmc_small, reduction.type = "pca", dims.use = 1:10)
#'
BuildSNN <- function(object, ...) {
  UseMethod(generic = 'BuildSNN', object = object)
}
