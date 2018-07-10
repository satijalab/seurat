#' SNN Graph Construction
#'
#' Constructs a Shared Nearest Neighbor (SNN) Graph for a given dataset. We
#' first determine the k-nearest neighbors of each cell. We use this knn graph
#' to construct the SNN graph by calculating the neighborhood overlap
#' (Jaccard index) between every cell and its k.param nearest neighbors.
#'
#' @param object Seurat object
#' @param k.param Defines k for the k-nearest neighbor algorithm
#' @param prune.SNN Sets the cutoff for acceptable Jaccard index when
#' computing the neighborhood overlap for the SNN construction. Any edges with
#' values less than or equal to this will be set to 0 and removed from the SNN
#' graph. Essentially sets the strigency of pruning (0 --- no pruning, 1 ---
#' prune everything).
#' @param nn.eps Error bound when performing nearest neighbor seach using RANN;
#' default of 0.0 implies exact nearest neighbor search
#' @param verbose Whether or not to print output to the console
#' @param force.recalc Force recalculation of SNN.
#'
#' @importFrom RANN nn2
#' @importFrom igraph plot.igraph graph.adjlist graph.adjacency E
#' @importFrom Matrix sparseMatrix
#' @return Returns the object with object@@snn filled
#' @export
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
#' @rdname BuildSNN
#' @export BuildSNN
#'
BuildSNN <- function(
  object,
  k.param,
  prune.SNN,
  nn.eps,
  verbose,
  force.recalc,
  ...
) {
  UseMethod(generic = 'BuildSNN', object = object)
}
