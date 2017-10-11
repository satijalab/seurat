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
#' @importFrom FNN get.knn
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
  k.scale = 10,
  plot.SNN = FALSE,
  prune.SNN = 1/15,
  print.output = TRUE,
  distance.matrix = NULL,
  force.recalc = FALSE
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
    my.knn <- get.knn(
      data <- as.matrix(x = data.use),
      k = min(k.scale * k.param, n.cells - 1)
    )
    nn.ranked <- cbind(1:n.cells, my.knn$nn.index[, 1:(k.param-1)])
    nn.large <- my.knn$nn.index
  } else {
    cat("Building SNN based on a provided distance matrix\n", file = stderr())
    n <- nrow(x = distance.matrix)
    k.for.nn <- k.param * k.scale
    knn.mat <- matrix(data = 0, ncol = k.for.nn, nrow = n)
    knd.mat <- knn.mat
    for (i in 1:n){
      knn.mat[i, ] <- order(data.use[i, ])[1:k.for.nn]
      knd.mat[i, ] <- data.use[i, knn.mat[i, ]]
    }
    nn.large <- knn.mat[, 2:(min(n, k.for.nn))]
    nn.ranked <- knn.mat[, 1:k.param]
  }
  w <- CalcSNNSparse(
    cell.names = object@cell.names,
    k.param = k.param,
    nn.large = nn.large,
    nn.ranked = nn.ranked,
    prune.SNN = prune.SNN,
    print.output = print.output
  )
  object@snn <- w
  if (plot.SNN) {
    if(!"tsne" %in% names(object@dr)) {
      warning("Please compute a tSNE for SNN visualization. See RunTSNE().")
    } else {
      if (nrow(object@dr$tsne@cell.embeddings) != length(x = object@cell.names)) {
        warning("Please compute a tSNE for SNN visualization. See RunTSNE().")
      } else {
        net <- graph.adjacency(
          adjmatrix = as.matrix(w),
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

# Function to convert the knn graph into the snn graph. Stored in a sparse
# representation.

# @param cell.names    A vector of cell names which will correspond to the row/
#                      column names of the SNN
# @param k.param       Defines nearest neighborhood when computing NN graph
# @param nn.large      Full KNN graph (computed with get.knn with k set to
#                      k.param * k.scale)
# @param nn.ranked     Subset of full KNN graph that only contains the first
#                      k.param nearest neighbors. Used to define Jaccard
#                      distances between any two cells
# @param prune.snn     Sets the cutoff for acceptable Jaccard distances when
#                      computing the neighborhood overlap for the SNN
#                      construction. Any edges with values less than or equal to
#                      this will be set to 0 and removed from the SNN graph.
#                      Essentially sets the strigency of pruning (0 --- no
#                      pruning, 1 --- prune everything).
# @param print.output  Whether or not to print output to the console
# @return              Returns an adjacency matrix representation of the SNN
#                      graph
#
#' @importFrom utils txtProgressBar setTxtProgressBar
#
CalcSNNSparse <- function(
  cell.names,
  k.param,
  nn.large,
  nn.ranked,
  prune.SNN,
  print.output
) {
  n.cells <- length(cell.names)
  counter <- 1
  idx1 <- vector(mode = "integer", length = n.cells ^ 2 / k.param)
  idx2 <- vector(mode = "integer", length = n.cells ^ 2 / k.param)
  edge.weight <- vector(mode = "double", length = n.cells ^ 2 / k.param)
  id <- 1
  # fill out the adjacency matrix w with edge weights only between your target
  # cell and its k.scale*k.param-nearest neighbors
  # speed things up (don't have to calculate all pairwise distances)
  # define the edge weights with Jaccard distance
  if (print.output) {
    print("Constructing SNN")
    pb <- txtProgressBar(min = 0, max = n.cells, style = 3)
  }
  for (i in 1:n.cells) {
    for (j in 1:ncol(x = nn.large)) {
      s <- intersect(x = nn.ranked[i, ], y = nn.ranked[nn.large[i, j], ])
      u <- union(nn.ranked[i, ], nn.ranked[nn.large[i, j], ])
      e <- length(x = s) / length(x = u)
      if (e > prune.SNN) {
        idx1[id] <- i
        idx2[id] <- nn.large[i, j]
        edge.weight[id] <- e
        id <- id + 1
      }
    }
    if (print.output) {
      setTxtProgressBar(pb = pb, value = i)
    }
  }
  if (print.output) {
    close(con = pb)
  }
  idx1 <- idx1[! is.na(x = idx1) & idx1 != 0]
  idx2 <- idx2[! is.na(x = idx2) & idx2 != 0]
  edge.weight <- edge.weight[! is.na(x = edge.weight) & edge.weight != 0]
  w <- sparseMatrix(
    i = idx1,
    j = idx2,
    x = edge.weight,
    dims = c(n.cells, n.cells)
  )
  diag(x = w) <- 1
  rownames(x = w) <- cell.names
  colnames(x = w) <- cell.names
  return(w)
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
