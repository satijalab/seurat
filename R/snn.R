#' @include seurat.R
NULL
#' SNN Graph Construction
#'
#' Constructs a Shared Nearest Neighbor (SNN) Graph for a given dataset. We
#' first determine the k-nearest neighbors of each cell. We use this knn graph
#' to construct the SNN graph by calculating the neighborhood overlap
#' (Jaccard index) between every cell and its k.param nearest neighbors.
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
#' @param plot.SNN Plot the SNN graph
#' @param prune.SNN Sets the cutoff for acceptable Jaccard index when
#' computing the neighborhood overlap for the SNN construction. Any edges with
#' values less than or equal to this will be set to 0 and removed from the SNN
#' graph. Essentially sets the strigency of pruning (0 --- no pruning, 1 ---
#' prune everything).
#' @param print.output Whether or not to print output to the console
#' @param distance.matrix Build SNN from distance matrix (experimental)
#' @param force.recalc Force recalculation of SNN.
#' @param filename Write SNN directly to file named here as an edge list
#' compatible with FindClusters
#' @param save.SNN Default behavior is to store the SNN in object@@snn. Setting
#' to FALSE can be used together with a provided filename to only write the SNN
#' out as an edge file to disk.
#' @param nn.eps Error bound when performing nearest neighbor seach using RANN;
#' default of 0.0 implies exact nearest neighbor search
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
BuildSNN <- function(
  object,
  genes.use = NULL,
  reduction.type = "pca",
  dims.use = NULL,
  k.param = 10,
  plot.SNN = FALSE,
  prune.SNN = 1/15,
  print.output = TRUE,
  distance.matrix = NULL,
  force.recalc = FALSE,
  filename = NULL,
  save.SNN = TRUE,
  nn.eps = 0
) {
  if (! is.null(x = distance.matrix)) {
    data.use <- distance.matrix
    force.recalc <- TRUE
  } else if (is.null(x = dims.use)) {
    genes.use <- SetIfNull(x = genes.use, default = object@var.genes)
    data.use <- t(x = as.matrix(x = object@data[genes.use, ]))
  } else {
    data.use <- GetCellEmbeddings(object = object,
                                  reduction.type = reduction.type,
                                  dims.use = dims.use)
  }
  parameters.to.store <- as.list(environment(), all = TRUE)[names(formals("BuildSNN"))]
  parameters.to.store$object <- NULL
  parameters.to.store$distance.matrix <- NULL
  parameters.to.store$print.output <- NULL
  if (CalcInfoExists(object, "BuildSNN") && ! force.recalc){
    old.parameters <- GetAllCalcParam(object, "BuildSNN")
    old.parameters$time <- NULL
    old.parameters$print.output <- NULL
    if(all(all.equal(old.parameters, parameters.to.store) == TRUE)){
      warning("Build parameters exactly match those of already computed and stored SNN. To force recalculation, set force.recalc to TRUE.")
      return(object)
    }
  }
  object <- SetCalcParams(object = object,
                          calculation = "BuildSNN",
                          ... = parameters.to.store)
  n.cells <- nrow(x = data.use)
  if (n.cells < k.param) {
    warning("k.param set larger than number of cells. Setting k.param to number of cells - 1.")
    k.param <- n.cells - 1
  }
  # find the k-nearest neighbors for each single cell
  if (is.null(x = distance.matrix)) {
    if (print.output) {
      cat("Computing nearest neighbor graph\n", file = stderr())
    }
    my.knn <- nn2(
        data = data.use,
        k = k.param,
        searchtype = 'standard',
        eps = nn.eps)
    nn.ranked <- my.knn$nn.idx
  } else {
    if (print.output) {
      cat("Building SNN based on a provided distance matrix\n", file = stderr())
    }
    n <- nrow(x = distance.matrix)
    k.for.nn <- k.param
    knn.mat <- matrix(data = 0, ncol = k.for.nn, nrow = n)
    knd.mat <- knn.mat
    for (i in 1:n){
      knn.mat[i, ] <- order(data.use[i, ])[1:k.for.nn]
      knd.mat[i, ] <- data.use[i, knn.mat[i, ]]
    }
    nn.ranked <- knn.mat[, 1:k.param]
  }
  if (print.output) {
    cat("Computing SNN\n", file = stderr())
  }
  if (save.SNN | is.null(filename)) {
    object@snn <- ComputeSNN(nn_ranked = nn.ranked,
                             prune = prune.SNN)
    rownames(object@snn) <- object@cell.names
    colnames(object@snn) <- object@cell.names
    if (!is.null(filename)) {
      WriteEdgeFile(snn = object@snn, filename = filename, display_progress = print.output)
    }
  } else {
    DirectSNNToFile(nn_ranked = nn.ranked, prune = prune.SNN,
                    display_progress = print.output, filename = filename)
  }
  if (plot.SNN & save.SNN) {
    if(!"tsne" %in% names(object@dr)) {
      warning("Please compute a tSNE for SNN visualization. See RunTSNE().")
    } else {
      if (nrow(object@dr$tsne@cell.embeddings) != length(x = object@cell.names)) {
        warning("Please compute a tSNE for SNN visualization. See RunTSNE().")
      } else {
        net <- graph.adjacency(
          adjmatrix = as.matrix(object@snn),
          mode = "undirected",
          weighted = TRUE,
          diag = FALSE
        )
        plot.igraph(
          x = net,
          layout = as.matrix(x = object@dr$tsne@cell.embeddings),
          edge.width = E(graph = net)$weight,
          vertex.label = NA,
          vertex.size = 0
        )
      }
    }
  }
  return(object)
}

# This function calculates the pairwise connectivity of clusters.

# @param object  Seurat object containing the snn graph and cluster assignments
# @return        matrix with all pairwise connectivities calculated

CalcConnectivity <- function(object) {
  SNN <- object@snn
  cluster.names <- unique(x = object@ident)
  num.clusters <- length(x = cluster.names)
  connectivity <- matrix(data = 0, nrow = num.clusters, ncol = num.clusters)
  rownames(x = connectivity) <- cluster.names
  colnames(x = connectivity) <- cluster.names
  n <- 1
  for (i in cluster.names) {
    for (j in cluster.names[-(1:n)]) {
      subSNN <- SNN[
        match(x = WhichCells(object = object, ident = i), colnames(x = SNN)),
        match(x = WhichCells(object = object, ident = j), rownames(x = SNN))
      ]
      if (is.object(x = subSNN)) {
        connectivity[i, j] <- sum(subSNN) / (nrow(x = subSNN) * ncol(x = subSNN))
      } else {
        connectivity[i, j] <- mean(x = subSNN)
      }
    }
    n <- n + 1
  }
  return(connectivity)
}
