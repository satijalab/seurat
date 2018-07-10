#' @include seurat.R
NULL

#' Cluster Determination
#'
#' Identify clusters of cells by a shared nearest neighbor (SNN) modularity
#' optimization based clustering algorithm. First calculate k-nearest neighbors
#' and construct the SNN graph. Then optimize the modularity function to
#' determine clusters. For a full description of the algorithms, see Waltman and
#' van Eck (2013) \emph{The European Physical Journal B}.
#'
#' @param object Seurat object
#' @param modularity.fxn Modularity function (1 = standard; 2 = alternative).
#' @param resolution Value of the resolution parameter, use a value above
#' (below) 1.0 if you want to obtain a larger (smaller) number of communities.
#' @param algorithm Algorithm for modularity optimization (1 = original Louvain
#' algorithm; 2 = Louvain algorithm with multilevel refinement; 3 = SLM
#' algorithm).
#' @param n.start Number of random starts.
#' @param n.iter Maximal number of iterations per random start.
#' @param random.seed Seed of the random number generator.
#' @param temp.file.location Directory where intermediate files will be written.
#' Specify the ABSOLUTE path.
#' @param edge.file.name Edge file to use as input for modularity optimizer jar.
#' @param verbose Print output
#'
#' @importFrom Matrix sparseMatrix
#' @importFrom methods .hasSlot
#' @importFrom igraph plot.igraph graph.adjlist
#'
#' @return Returns a Seurat object and optionally the SNN matrix,
#'         object idents have been updated with new cluster info
#'
#' @export
#'
#' @rdname FindClusters
#' @export FindClusters
#'
FindClusters <- function(
  object,
  modularity.fxn,
  resolution,
  algorithm,
  n.start,
  n.iter,
  random.seed,
  temp.file.location,
  edge.file.name,
  verbose,
  ...
) {
  UseMethod(generic = 'FindClusters', object = object)
}
