# Runs the modularity optimizer (C++ port of java program ModularityOptimizer.jar)
#
#
# @param SNN SNN              matrix to use as input for the clustering
#                             algorithms
# @param modularity           Modularity function to use in clustering (1 =
#                             standard; 2 = alternative).
# @param resolution           Value of the resolution parameter, use a value
#                             above (below) 1.0 if you want to obtain a larger
#                             (smaller) number of communities.
# @param algorithm            Algorithm for modularity optimization (1 =
#                             original Louvain algorithm; 2 = Louvain algorithm
#                             with multilevel refinement; 3 = SLM algorithm)
# @param n.start              Number of random starts.
# @param n.iter               Maximal number of iterations per random start.
# @param random.seed          Seed of the random number generator
# @param print.output         Whether or not to print output to the console
# @param temp.file.location   Deprecated and no longer used.
# @param edge.file.name       Path to edge file to use
# @return                     Seurat object with identities set to the results
#                             of the clustering procedure.
#
#' @importFrom utils read.table write.table
#
RunModularityClustering <- function(
  SNN = matrix(),
  modularity = 1,
  resolution = 0.8,
  algorithm = 1,
  n.start = 100,
  n.iter = 10,
  random.seed = 0,
  print.output = TRUE,
  temp.file.location = NULL,
  edge.file.name = NULL
) {

  edge_file = SetIfNull(x = edge.file.name, default="")
  clusters = RunModularityClusteringCpp(SNN, modularity, resolution, algorithm,
                             n.start, n.iter, random.seed, print.output, edge_file)
  return(clusters)
}

# Group single cells that make up their own cluster in with the cluster they are
# most connected to.
#
# @param ids     Named vector of cluster ids
# @param SNN     SNN graph used in clustering
# @return        Returns Seurat object with all singletons merged with most
#                connected cluster

GroupSingletons <- function(ids, SNN, verbose) {
  # identify singletons
  singletons <- c()
  singletons <- names(which(table(ids) == 1))
  singletons <- intersect(unique(ids), singletons)
  # calculate connectivity of singletons to other clusters, add singleton
  # to cluster it is most connected to
  cluster_names <- as.character(unique(x = ids))
  cluster_names <- setdiff(x = cluster_names, y = singletons)
  connectivity <- vector(mode = "numeric", length = length(x = cluster_names))
  names(x = connectivity) <- cluster_names
  new.ids <- ids
  for (i in singletons) {
    i.cells <- names(which(ids == i))
    for (j in cluster_names) {
      j.cells <- names(which(ids == j))
      subSNN <- SNN[i.cells, j.cells]
      set.seed(1) # to match previous behavior, random seed being set in WhichCells
      if (is.object(x = subSNN)) {
        connectivity[j] <- sum(subSNN) / (nrow(x = subSNN) * ncol(x = subSNN))
      } else {
        connectivity[j] <- mean(x = subSNN)
      }
    }
    m <- max(connectivity, na.rm = T)
    mi <- which(x = connectivity == m, arr.ind = TRUE)
    closest_cluster <- sample(x = names(x = connectivity[mi]), 1)
    ids[i.cells] <- closest_cluster
  }
  if (length(x = singletons) > 0 && verbose) {
    message(paste(
      length(x = singletons),
      "singletons identified.",
      length(x = unique(x = ids)),
      "final clusters."
    ))
  }
  return(ids)
}


# Set up kmeans class
# This is an infrequently used slot, but some people still find it very useful to do kmeans clustering
# and in particular, to do so at the gene level
# potential to be updated in the future

kmeans.info <- setClass(
  Class = "kmeans.info",
  slots = list(
    gene.kmeans.obj = "ANY",
    cell.kmeans.obj = "ANY"
  )
)
