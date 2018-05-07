#' @include seurat.R
#' @include snn_generics.R
NULL

#' @describeIn BuildSNN ...
#' @export BuildSNN.seurat
#' @method BuildSNN seurat
#'
BuildSNN.seurat <- function(
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
  if (!is.null(x = distance.matrix)) {
    data.use <- distance.matrix
    force.recalc <- TRUE
  } else if (is.null(x = dims.use)) {
    genes.use <- SetIfNull(x = genes.use, default = object@var.genes)
    data.use <- t(x = as.matrix(x = object@data[genes.use, ]))
  } else {
    data.use <- GetCellEmbeddings(
      object = object,
      reduction.type = reduction.type,
      dims.use = dims.use
    )
  }
  parameters.to.store <- as.list(environment(), all = TRUE)[names(formals())]
  parameters.to.store$object <- NULL
  parameters.to.store$distance.matrix <- NULL
  parameters.to.store$print.output <- NULL
  if (CalcInfoExists(object, "BuildSNN") && !force.recalc) {
    old.parameters <- GetAllCalcParam(object, "BuildSNN")
    old.parameters$time <- NULL
    old.parameters$print.output <- NULL
    if (all(all.equal(old.parameters, parameters.to.store) == TRUE)) {
      warning("Build parameters exactly match those of already computed and stored SNN. To force recalculation, set force.recalc to TRUE.")
      return(object)
    }
  }
  object <- SetCalcParams(
    object = object,
    calculation = "BuildSNN",
    ... = parameters.to.store
  )
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
    for (i in 1:n) {
      knn.mat[i, ] <- order(data.use[i, ])[1:k.for.nn]
      knd.mat[i, ] <- data.use[i, knn.mat[i, ]]
    }
    nn.ranked <- knn.mat[, 1:k.param]
  }
  if (print.output) {
    cat("Computing SNN\n", file = stderr())
  }
  if (save.SNN | is.null(filename)) {
    object@snn <- ComputeSNN(
      nn_ranked = nn.ranked,
      prune = prune.SNN
    )
    rownames(object@snn) <- object@cell.names
    colnames(object@snn) <- object@cell.names
    if (!is.null(filename)) {
      WriteEdgeFile(
        snn = object@snn,
        filename = filename,
        display_progress = print.output
      )
    }
  } else {
    DirectSNNToFile(
      nn_ranked = nn.ranked,
      prune = prune.SNN,
      display_progress = print.output,
      filename = filename
    )
  }
  if (plot.SNN & save.SNN) {
    if (!"tsne" %in% names(object@dr)) {
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

#' @param overwrite Overwrite existing SNN graph?
#'
#' @describeIn BuildSNN ...
#' @export
#' @method BuildSNN loom
#'
BuildSNN.loom <- function(
  object,
  genes.use = NULL,
  reduction.type  = 'pca',
  dims.use = NULL,
  k.param = 10,
  k.scale = 10,
  prune.SNN = 1/15,
  print.output = TRUE,
  distance.matrix = NULL,
  filename = NULL,
  save.SNN = TRUE,
  nn.eps = 0,
  overwrite = FALSE
) {
  if (!is.null(x = distance.matrix)) {
    data.use <- distance.matrix
  } else if (is.null(x = dims.use)) {
    genes.use <- SetIfNull(x = genes.use, default = which(object[['row_attrs/var_genes']][]))
    data.use <- t(GetAssayData.loom(object = object, genes.use = genes.use))
  } else {
    # data.use <- as.matrix(object[[reduction.data]])
    data.use <- GetDimReduction(
      object = object,
      reduction.type = reduction.type,
      slot = 'cell_embeddings'
    )
    data.use <- data.use[, dims.use]
  }
  # needs code for storing calculation parameters
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
      eps = nn.eps
    )
    nn.ranked <- my.knn$nn.idx
  } else {
    if (print.output) {
      cat("Building SNN based on a provided distance matrix\n", file = stderr())
    }
    n <- nrow(x = distance.matrix)
    k.for.nn <- k.param
    knn.mat <- matrix(data = 0, ncol = k.for.nn, nrow = n)
    knd.mat <- knn.mat
    for (i in 1:n) {
      knn.mat[i, ] <- order(data.use[i, ])[1:k.for.nn]
      knd.mat[i, ] <- data.use[i, knn.mat[i, ]]
    }
    nn.ranked <- knn.mat[, 1:k.param]
  }
  if (print.output) {
    cat("Computing SNN\n", file = stderr())
  }
  # needs option to store SNN matrix in loom object - should be supported in
  # loom2, for now just write out edge file
  if (save.SNN | is.null(filename)) {
    if ('SNN' %in% names(x = object[['col_graphs']])) {
      if (overwrite) {
        object[['col_graphs']]$link_delete(name = 'SNN')
      } else {
        stop("SNN already exists!")
      }
    }
    snn <- ComputeSNN(
      nn_ranked = nn.ranked,
      prune = prune.SNN
    )
    object$add.graph.matrix(mat = snn, name = 'SNN')
    if (!is.null(filename)) {
      WriteEdgeFile(
        snn = snn,
        filename = filename,
        display_progress = print.output
      )
    }
  } else {
    DirectSNNToFile(
      nn_ranked = nn.ranked,
      prune = prune.SNN,
      display_progress = print.output,
      filename = filename
    )
  }
  object$flush()
  gc(verbose = FALSE)
  invisible(x = object)
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
