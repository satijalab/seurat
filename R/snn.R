#' @include seurat.R
NULL
#' SNN Graph Construction
#'
#' Construct a Shared Nearest Neighbor (SNN) Graph for a given
#' dataset.
#'
#'
#' @param object Seurat object
#' @param genes.use Gene expression data
#' @param pc.use Which PCs to use for construction of the SNN graph
#' @param k.param Defines k for the k-nearest neighbor algorithm
#' @param k.scale granularity option for k.param
#' @param plot.SNN Plot the SNN graph
#' @param prune.SNN Stringency of pruning for the SNN graph (0 - no pruning,
#'        1 - prune everything)
#' @param do.sparse Whether to compute and return the SNN graph as a sparse
#' matrix or not
#' @param print.output Whether or not to print output to the console
#' @param distance.matrix Build SNN from distance matrix (experimental)
#' @importFrom FNN get.knn
#' @importFrom igraph plot.igraph graph.adjlist graph.adjacency E
#' @importFrom Matrix sparseMatrix
#' @return Returns the object with object@@snn.k and either
#' object@@snn.dense or object@@snn.sparse filled depending on the option
#' set
#' @export
BuildSNN <- function(
  object,
  genes.use = NULL,
  reduction.type="pca",
  pc.use = NULL,
  k.param = 10,
  k.scale = 10,
  plot.SNN = FALSE,
  prune.SNN = 1/15,
  do.sparse = FALSE,
  print.output = TRUE,
  distance.matrix = NULL
) {
  if (! is.null(x = distance.matrix)) {
    data.use <- distance.matrix
  } else if (is.null(x = genes.use) && is.null(x = pc.use)) {
    genes.use <- object@var.genes
    data.use <- t(x = as.matrix(x = object@data[genes.use, ]))
  } else if (! is.null(x = pc.use)) {
    #data.use <- as.matrix(object@pca.rot[, pc.use])
    dim_scores <- eval(expr = parse(
      text = paste0("object@dr$", reduction.type, "@rotation")
    ))
    data.use <- dim_scores[, pc.use]
  } else if (!is.null(genes.use) && is.null(pc.use)) {
    data.use <- t(x = as.matrix(x = object@data[genes.use, ]))
  } else {
      stop("Data error!")
  }
  n.cells <- nrow(x = data.use)
  if (n.cells < k.param) {
    warning("k.param set larger than number of cells. Setting k.param to number of cells - 1.")
    k.param <- n.cells - 1
  }
  #find the k-nearest neighbors for each single cell
  if (is.null(x = distance.matrix)) {
    my.knn <- get.knn(
      data <- as.matrix(x = data.use),
      k = min(k.scale * k.param, n.cells - 1)
    )
    nn.ranked <- cbind(1:n.cells, my.knn$nn.index[, 1:(k.param-1)])
    nn.large <- my.knn$nn.index
  } else {
    warning("Building SNN based on a provided distance matrix")
    n <- nrow(x = distance.matrix)
    k.for.nn <- k.param * k.scale
    knn.mat <- matrix(data = 0, ncol = k.for.nn, nrow = n)
    knd.mat <- knn.mat
    for (i in 1:n){
      knn.mat[i, ] <- order(data.use[i, ])[1:k.for.nn]
      knd.mat[i, ] <- data.use[i, knn.mat[i, ]]
    }
    nn.large <- knn.mat
    nn.ranked <- knn.mat[, 2:k.param]
  }
  if (do.sparse) {
    w <- CalcSNNSparse(
      object = object,
      n.cells = n.cells,
      k.param = k.param,
      nn.large = nn.large,
      nn.ranked = nn.ranked,
      prune.SNN = prune.SNN,
      print.output = print.output
    )
  } else {
      w <- CalcSNNDense(
        object = object,
        n.cells = n.cells,
        nn.large = nn.large,
        nn.ranked = nn.ranked,
        prune.SNN = prune.SNN,
        print.output = print.output
      )
  }
  if (plot.SNN) {
    if (length(x = object@dr$tsne@rotation) < 1) {
      warning("Please compute a tSNE for SNN visualization. See RunTSNE().")
    } else {
      net <- graph.adjacency(
        adjmatrix = w,
        mode = "undirected",
        weighted = TRUE,
        diag = FALSE
      )
      plot.igraph(
        x = net,
        layout = as.matrix(x = object@dr$tsne@rotation),
        edge.width = E(graph = net)$weight,
        vertex.label = NA,
        vertex.size = 0
      )
    }
  }
  #only allow one of the snn matrix slots to be filled
  object@snn.k <- k.param
  if (do.sparse) {
    object@snn.sparse <- w
    object@snn.dense <- matrix()
  } else {
    object@snn.dense <- w
    object@snn.sparse <- sparseMatrix(i = 1, j = 1, x = 1)
  }
  return(object)
}

CalcSNNSparse <- function(
  object,
  n.cells,
  k.param,
  nn.large,
  nn.ranked,
  prune.SNN,
  print.output
) {
  counter <- 1
  idx1 <- vector(mode = "integer", length = n.cells ^ 2 / k.param)
  idx2 <- vector(mode = "integer", length = n.cells ^ 2 / k.param)
  edge.weight <- vector(mode = "double", length = n.cells ^ 2 / k.param)
  id <- 1
  #fill out the adjacency matrix w with edge weights only between your target
  #cell and its 10*k.param-nearest neighbors
  #speed things up (don't have to calculate all pairwise distances)
  #define the edge weights with Jaccard distance
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
  rownames(x = w) <- object@cell.names
  colnames(x = w) <- object@cell.names
  return(w)
}

CalcSNNDense <- function(
  object,
  n.cells,
  nn.large,
  nn.ranked,
  prune.SNN,
  print.output = TRUE
) {
  counter <- 1
  w <- matrix(data = 0, nrow = n.cells, ncol = n.cells)
  rownames(x = w) <- object@cell.names
  colnames(x = w) <- object@cell.names
  diag(x = w) <- 1
  #fill out the adjacency matrix w with edge weights only between your target
  #cell and its 10*k.param-nearest neighbors
  #speed things up (don't have to calculate all pairwise distances)
  if (print.output){
    print("Constructing SNN")
    pb <- txtProgressBar(min = 0, max = n.cells, style = 3)
  }
  for (i in 1:n.cells) {
    for (j in 1:ncol(x = nn.large)) {
      s <- intersect(x = nn.ranked[i, ], y = nn.ranked[nn.large[i, j], ])
      u <- union(nn.ranked[i, ], nn.ranked[nn.large[i, j], ])
      e <- length(x = s) / length(x = u)
      if (e > prune.SNN) {
        w[i, nn.large[i, j]] <- e
      } else {
        w[i,nn.large[i, j]] <- 0
      }
    }
    if (print.output) {
      setTxtProgressBar(pb = pb, value = i)
    }
  }
  if (print.output) {
    close(con = pb)
  }
  return(w)
}

CalcConnectivity <- function(object, SNN) {
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
