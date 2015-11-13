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
#' @param update Adjust how verbose the output is
#' @importFrom FNN get.knn
#' @importFrom igraph plot.igraph graph.adjlist
#' @importFrom Matrix sparseMatrix
#' @return Returns the object with object@@snn.k and either 
#' object@@snn.dense or object@@snn.sparse filled depending on the option
#' set
#' @export
setGeneric("BuildSNN", function(object, genes.use = NULL, pc.use = NULL, 
                                k.param = 10, k.scale = 10, plot.SNN = FALSE, 
                                prune.SNN = 0.1, do.sparse = FALSE, 
                                update = 0.25)  
standardGeneric("BuildSNN"))
#' @export
setMethod("BuildSNN", signature = "seurat",
          function(object, genes.use = NULL, pc.use = NULL, k.param = 10, 
                   k.scale = 10, plot.SNN = FALSE, prune.SNN = 0.1,
                   do.sparse = FALSE, update = 0.25) {

  if (is.null(genes.use) && is.null(pc.use)) {
    genes.use <- object@var.genes
    data.use <- t(as.matrix(object@data[genes.use, ]))
  } else if (!is.null(pc.use)) {
      data.use <- as.matrix(object@pca.rot[, pc.use])
  } else if (!is.null(genes.use) && is.null(pc.use)) {
      data.use <- t(as.matrix(object@data[genes.use, ]))
  } else {
      stop("Data error!")
  }
                
  n.cells <- nrow(data.use)
  if (k.param * k.scale > n.cells){
    stop("k.scale x k.param can't be larger than the number of cells")
  }

  #find the k-nearest neighbors for each single cell
  my.knn <- get.knn(as.matrix(data.use), k = min(k.scale * k.param, n.cells))
  nn.ranked <- cbind(1:n.cells, my.knn$nn.index[, 1:(k.param-1)])
  nn.large <- my.knn$nn.index
  if (do.sparse){
    w <- CalcSNNSparse(object, n.cells, k.param, nn.large, nn.ranked, prune.SNN, 
                       update)
  } else {
      w <- CalcSNNDense(object, n.cells, nn.large, nn.ranked, prune.SNN, update)
      if (plot.SNN==TRUE) {
        net <- graph.adjacency(w, mode = "undirected", weighted = TRUE, 
                               diag = FALSE)
        plot.igraph(net, layout = as.matrix(object@tsne.rot), 
                    edge.width = E(net)$weight, vertex.label = NA, 
                    vertex.size = 0)
      }
  }
  
  #only allow one of the snn matrix slots to be filled
  object@snn.k <- k.param
  if (do.sparse == TRUE) {
    object@snn.sparse <- w
    object@snn.dense <- matrix()
  } else {
    object@snn.dense <- w
    object@snn.sparse <- sparseMatrix(1, 1, x = 1)
  }
  return(object)
})

CalcSNNSparse <- function(object, n.cells, k.param, nn.large, nn.ranked, 
                         prune.SNN, update) {
  counter <- 1
  idx1 <- vector(mode = "integer", length = n.cells ^ 2 / k.param)
  idx2 <- vector(mode = "integer", length = n.cells ^ 2 / k.param)
  edge.weight <- vector(mode = "double", length = n.cells ^ 2 / k.param)
  id <- 1
  
  #fill out the adjacency matrix w with edge weights only between your target 
  #cell and its 10*k.param-nearest neighbors
  #speed things up (don't have to calculate all pairwise distances)
  #define the edge weights with Jaccard distance
  for (i in 1:n.cells) {
    for (j in 1:ncol(nn.large)) {
      s <- intersect(nn.ranked[i, ], nn.ranked[nn.large[i, j], ])
      u <- union(nn.ranked[i, ], nn.ranked[nn.large[i, j], ])
      e <- length(s) / length(u)
      if (e > prune.SNN) {
        idx1[id] <- i
        idx2[id] <- nn.large[i, j]
        edge.weight[id] <- e
        id <- id + 1
      }
    }
    if (i == round(counter * n.cells * update)) {
      print(paste("SNN : processed ", i, " cells", sep = ""))
      counter <- counter + 1;
    }
  }
  idx1 <- idx1[!is.na(idx1) & idx1 != 0]
  idx2 <- idx2[!is.na(idx2) & idx2 != 0]
  edge.weight <- edge.weight[!is.na(edge.weight) & edge.weight != 0]
  w <- sparseMatrix(i = idx1, j = idx2, x = edge.weight, 
                    dims = c(n.cells, n.cells))
  diag(w) <- 1
  rownames(w) <- object@cell.names
  colnames(w) <- object@cell.names
  return(w)
}

CalcSNNDense <- function(object, n.cells, nn.large, nn.ranked, prune.SNN, 
                         update) {
  counter <- 1
  w <- matrix(0, n.cells, n.cells)
  rownames(w) <- object@cell.names
  colnames(w) <- object@cell.names
  diag(w) <- 1
  #fill out the adjacency matrix w with edge weights only between your target 
  #cell and its 10*k.param-nearest neighbors
  #speed things up (don't have to calculate all pairwise distances)
  for (i in 1:n.cells) {
    for (j in 1:ncol(nn.large)) {
      s <- intersect(nn.ranked[i, ], nn.ranked[nn.large[i, j], ])
      u <- union(nn.ranked[i, ], nn.ranked[nn.large[i, j], ])
      e <- length(s) / length(u)
      if (e > prune.SNN) {
        w[i, nn.large[i, j]] <- e
      } else {
        w[i,nn.large[i, j]] <- 0 
      }
    }
    if (i == round(counter * n.cells * update)) {
      print(paste("SNN : processed ", i, " cells", sep = ""))
      counter <- counter + 1;
    }
  }
  return(w)
}

CalcConnectivity <- function(object, SNN) {
  cluster.names <- unique(object@ident)
  num.clusters <- length(cluster.names)
  connectivity <- matrix(0, nrow = num.clusters, ncol = num.clusters)
  rownames(connectivity) = cluster.names
  colnames(connectivity) = cluster.names
  n <- 1
  for (i in cluster.names) {
    for (j in cluster.names[-(1:n)]) {
      subSNN <- SNN[match(which.cells(object, i), colnames(SNN)), 
                    match(which.cells(object, j), rownames(SNN))]
      if (is.object(subSNN)) {
        connectivity[i, j] <- sum(subSNN) / (nrow(subSNN) * ncol(subSNN))
      } else {
        connectivity[i, j] <- mean(subSNN)
      }
    }
    n <- n + 1
  }
  return(connectivity)
}