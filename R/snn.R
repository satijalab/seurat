
#' @param distance.matrix Boolean value of whether the provided matrix is a
#' distance matrix
#' @describeIn BuildSNN Build an SNN on a given matrix
#' @export
#' @method BuildSNN default
#'
BuildSNN.default <- function(
  object,
  distance.matrix = FALSE,
  k.param = 10,
  prune.SNN = 1/15,
  nn.eps = 0,
  verbose = TRUE,
  force.recalc = FALSE
) {
  n.cells <- nrow(x = object)
  if (n.cells < k.param) {
    warning("k.param set larger than number of cells. Setting k.param to number of cells - 1.")
    k.param <- n.cells - 1
  }

  # find the k-nearest neighbors for each single cell
  if (! distance.matrix) {
    if (verbose) {
      message("Computing nearest neighbor graph")
    }
    my.knn <- nn2(
      data = object,
      k = k.param,
      searchtype = 'standard',
      eps = nn.eps)
    nn.ranked <- my.knn$nn.idx
  } else {
    if (verbose) {
      message("Building SNN based on a provided distance matrix")
    }
    knn.mat <- matrix(data = 0, ncol = k.param, nrow = n.cells)
    knd.mat <- knn.mat
    for (i in 1:n.cells){
      knn.mat[i, ] <- order(data.use[i, ])[1:k.param]
      knd.mat[i, ] <- data.use[i, knn.mat[i, ]]
    }
    nn.ranked <- knn.mat[, 1:k.param]
  }
  if (verbose) {
    message("Computing SNN")
  }
  snn.matrix <- ComputeSNN(nn_ranked = nn.ranked,
                           prune = prune.SNN)
  rownames(snn.matrix) <- rownames(object)
  colnames(snn.matrix) <- rownames(object)
  snn.matrix <- as(object = snn.matrix, Class = "Graph")
  return(snn.matrix)
}

#' @param features.use A vector of feature names to use in construction of SNN
#' graph if building directly based on data rather than a dimensionally reduced
#' representation (i.e. PCs).
#' @param reduction.use Name of dimensional reduction technique to use in
#' construction of SNN graph. (e.g. "pca", "ica")
#' @param dims.use A vector of the dimensions to use in construction of the SNN
#' graph (e.g. To use the first 10 PCs, pass 1:10)
#'
#'
#' @describeIn BuildSNN Build an SNN on an Assay object
#' @export
#' @method BuildSNN Assay
#'
BuildSNN.Assay <- function(
  object,
  features.use = NULL,
  reduction.use = "pca",
  dims.use = NULL,
  k.param = 10,
  prune.SNN = 1/15,
  nn.eps = 0,
  verbose = TRUE,
  force.recalc = FALSE
) {
  if (is.null(dims.use)) {
    features.use <- features.use %||% VariableFeatures(object = object)
    data.use <- t(GetAssayData(object = object, slot = "data")[features.use, ])
  } else {
    data.use <- Embeddings(object = object)[, dims.use]
  }
  snn.matrix <- BuildSNN(object = data.use,
                         k.param = k.param,
                         prune.SNN = prune.SNN,
                         nn.eps = nn.eps,
                         verbose = verbose,
                         force.recalc = force.recalc)
  return(snn.matrix)
}

#' @param assay.use Assay to use in construction of SNN
#' @param features.use Features to use as input for building the SNN
#' @param reduction.use Reduction to use as input for building the SNN
#' @param dims.use Dimensions of reduction to use as input
#' @param do.plot Plot SNN graph on tSNE coordinates
#' @param graph.name Optional naming parameter for stored SNN graph. Default is
#' assay.name_snn.
#'
#' @describeIn BuildSNN Build an SNN on a Seurat object
#' @export
#' @method BuildSNN Seurat
#'
BuildSNN.Seurat <- function(
  object,
  assay.use = NULL,
  features.use = NULL,
  reduction.use = "pca",
  dims.use = NULL,
  k.param = 10,
  prune.SNN = 1/15,
  nn.eps = 0,
  verbose = TRUE,
  force.recalc = FALSE,
  do.plot = FALSE,
  graph.name = NULL
) {
  assay.use <- assay.use %||% DefaultAssay(object = object)
  assay.data <- GetAssay(object = object, assay.use = assay.use)
  snn.matrix <- BuildSNN(object = assay.data,
                         features.use = features.use,
                         reduction.use = reduction.use,
                         dims.use = dims.use,
                         k.param = k.param,
                         prune.SNN = prune.SNN,
                         nn.eps = nn.eps,
                         verbose = verbose,
                         force.recalc = force.recalc)

  graph.name <- graph.name %||% paste0(assay.use, "_snn")
  object[[graph.name]] <- snn.matrix

  if (do.plot) {
    if(!"tsne" %in% names(object@reductions)) {
      warning("Please compute a tSNE for SNN visualization. See RunTSNE().")
    } else {
      if (nrow(Embeddings(object = object[["tsne"]])) != ncol(object)) {
        warning("Please compute a tSNE for SNN visualization. See RunTSNE().")
      } else {
        net <- graph.adjacency(
          adjmatrix = as.matrix(snn.matrix),
          mode = "undirected",
          weighted = TRUE,
          diag = FALSE
        )
        plot.igraph(
          x = net,
          layout = as.matrix(x = Embeddings(object = object[["tsne"]])),
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
