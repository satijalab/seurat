#' @include seurat.R
NULL
#' Cluster Determination
#'
#' Identify clusters of cells by a shared nearest neighbor (SNN) quasi-clique
#' based clustering algorithm. First calculate k-nearest neighbors and construct 
#' the SNN graph. Then determine the quasi-cliques associated with each cell. 
#' Finally, merge the quasi-cliques into clusters. For a full description of the
#' algorithm, see Xu and Su (2015) \emph{Bioinformatics}. 
#' 
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
#' @param save.SNN Whether to return the SNN matrix or not. If true, returns a 
#'        list with the object as the first item
#'         and the SNN matrix as the second item.
#' @param r.param r defines the connectivity for the quasi-cliques. 
#'        Higher r gives a more compact subgraph
#' @param m.param m is the threshold for merging two quasi-cliques. 
#'        Higher m results in less merging
#' @param q Defines the percentage of quasi-cliques to examine for merging each 
#'        iteration
#' @param qup Determines how to change q once all possible merges have been made
#' @param update Adjust how verbose the output is
#' @param min.cluster.size Smallest allowed size for a cluster
#' @param do.sparse Option to store and use SNN matrix as a sparse matrix. 
#'        May be necessary datasets containing a large number of cells.
#' @param do.modularity Option to use modularity optimization for single cell 
#'        clustering.
#' @param modularity Modularity function (1 = standard; 2 = alternative).
#' @param resolution Value of the resolution parameter, use a value above 
#'        (below) 1.0 if you want to obtain a larger (smaller) number of 
#'        communities.
#' @param algorithm Algorithm for modularity optimization (1 = original Louvain 
#'        algorithm; 2 = Louvain algorithm with multilevel refinement; 
#'        3 = SLM algorithm).
#' @param n.start Number of random starts.
#' @param n.iter Maximal number of iterations per random start.
#' @param random.seed Seed of the random number generator.
#' @param print.output Whether or not to print output to the console (0 = no; 
#'        1 = yes).
#' @importFrom FNN get.knn
#' @importFrom igraph plot.igraph graph.adjlist
#' @importFrom Matrix sparseMatrix
#' @return Returns a Seurat object and optionally the SNN matrix, 
#'         object@@ident has been updated with new cluster info
#' @export
setGeneric("FindClusters", function(object, genes.use = NULL, pc.use = NULL, 
                                     k.param = 10, k.scale = 10,
                                     plot.SNN = FALSE, prune.SNN = 0.1, 
                                     save.SNN = FALSE, r.param = 0.7, 
                                     m.param = NULL, q = 0.1, qup = 0.1, 
                                     update = 0.25, min.cluster.size = 1, 
                                     do.sparse = FALSE, do.modularity = TRUE, 
                                     modularity = 1, resolution = 0.8, 
                                     algorithm = 1, n.start = 100, 
                                     n.iter = 10, random.seed = 0, 
                                     print.output = 1)  
standardGeneric("FindClusters"))
#' @export
setMethod("FindClusters", signature = "seurat",
          function(object, genes.use = NULL, pc.use = NULL, k.param = 10, 
                   k.scale = 10, plot.SNN = FALSE, prune.SNN = 0.1, 
                   save.SNN = FALSE, r.param = 0.7, m.param = NULL, 
                   q = 0.1, qup = 0.1, update = 0.25, min.cluster.size = 1, 
                   do.sparse = FALSE, do.modularity = TRUE, modularity = 1, 
                   resolution = 0.8, algorithm = 1, n.start = 100, n.iter = 10, 
                   random.seed = 0, print.output = 1){

  # if any SNN building parameters are provided, build a new SNN
  if (k.param != 10 | k.scale != 10) {
    object <- BuildSNN(object, genes.use, pc.use, k.param, k.scale, 
                        plot.SNN, prune.SNN, do.sparse, update)
  }
  
  # if the SNN hasn't been built yet, build it
  snn.built = FALSE
  if (.hasSlot(object, "snn.dense")) {
    if (length(object@snn.dense) > 1) {
      snn.built = TRUE
    }
  }
  if (.hasSlot(object, "snn.sparse")) {
    if (length(object@snn.sparse) > 1) {
      snn.built = TRUE
    }
  }
  if (!snn.built) {
    object <- BuildSNN(object, genes.use, pc.use, k.param, k.scale, 
                       plot.SNN, prune.SNN, do.sparse, update)
  }
  # deal with sparse SNNs
  # this part should be refactored given new slots to make it cleaner
  # (will require changing called functions as well)
  if (length(object@snn.sparse) > 1) {
    SNN.sp <- object@snn.sparse
    SNN.use <- matrix()
    do.sparse <- TRUE
  } else {
    SNN.use <- object@snn.dense
    SNN.sp <- sparseMatrix(1, 1, x = 1)
    do.sparse <- FALSE
  }

  if (do.modularity) {
    if (do.sparse) {
      object <- RunModularityClustering(object, SNN.sp, modularity, resolution, 
                                   algorithm, n.start, n.iter, random.seed, 
                                   print.output)
      object <- GroupSingletons(object, SNN.sp)
    } else {
      object <- RunModularityClustering(object, SNN.use, modularity, resolution, 
                                   algorithm, n.start, n.iter, random.seed, 
                                   print.output)
      object <- GroupSingletons(object, SNN.use)
    }
  } else {
    if (is.null(m.param)) {
      clusters <- r_wrapper(SNN.use, SNN.sp, r.param, m.param <- r.param, q, 
                            qup, update, min.cluster.size, do.sparse)
    } else {
      clusters <- r_wrapper(SNN.use, SNN.sp, r.param, m.param, q, qup, update, 
                            min.cluster.size, do.sparse)
    }
    clusters.list <- rep(1:length(clusters[[2]]), clusters[[2]])
    if (!is.null(clusters[[3]])) {
      clusters.list <- replace(clusters.list, 
                               seq(length(clusters.list)-tail(clusters[[2]],1),
                               length(clusters.list)), 0)
    }
    cells.use <- object@cell.names[unlist(clusters[[1]])]
    ident.use <- clusters.list
    object <- set.ident(object, cells.use, ident.use)
  }
            
  if (!save.SNN) {
    object@snn.sparse <- sparseMatrix(1, 1, x = 1)
    object@snn.dense <- matrix()
    object@snn.k <- integer()
  }
  return(object)
})

#' @export 
setGeneric("GetClusters", function(object) standardGeneric("GetClusters"))
setMethod("GetClusters", signature="seurat",
          function(object){
            return(data.frame(object@ident))
          }
)

#' @export 
setGeneric("SetClusters", function(object, clusters=NULL) 
  standardGeneric("SetClusters"))
setMethod("SetClusters", signature="seurat",
          function(object, clusters = NULL){
            cells.use <- rownames(clusters)
            ident.use <- as.numeric(clusters[, 1])
            object <- set.ident(object, cells.use, ident.use)
            return(object)
          }
)

#' @export 
setGeneric("SaveClusters", function(object, file) 
  standardGeneric("SaveClusters"))
setMethod("SaveClusters", signature="seurat",
          function(object, file){
            my.clusters <- GetClusters(object)
            write.table(my.clusters, file = file, sep="\t", quote = FALSE)
          }
)

#' @export 
setGeneric("NumberClusters", function(object) 
  standardGeneric("NumberClusters"))
setMethod("NumberClusters", signature="seurat",
          function(object){
            clusters <- unique(object@ident)
            if (typeof(clusters) == "integer") {
              n <- as.numeric(max(clusters)) + 1
              for (i in clusters) {
                object <- set.ident(object, cells.use = which.cells(object, i), 
                                    ident.use = n)
                n <- n + 1
              }
              clusters <- unique(object@ident)
            }
            n <- 1
            for (i in clusters) {
              object <- set.ident(object, cells.use = which.cells(object, i), 
                                  ident.use = n)
              n <- n + 1
            }
            return(object)
          }
)


RunModularityClustering <- function(object, SNN = matrix(), modularity = 1, 
                                    resolution = 0.8, algorithm = 1, 
                                    n.start = 100, n.iter = 10, random.seed = 0,
                                    print.output = 1){

  ModularityJarFile <- paste(system.file(package="Seurat"), 
                             "/java/ModularityOptimizer.jar", sep = "")
  diag(SNN) <- 0
  if (is.object(SNN)) {
    SNN <- as(SNN, "dgTMatrix")
    edge <- cbind(i = SNN@i, j = SNN@j, x = SNN@x)
  } else {
    edge <- cbind((which(SNN != 0, arr.ind = TRUE) - 1), 
                  SNN[which(SNN != 0, arr.ind = TRUE)])
  }
  rownames(edge) <- NULL
  colnames(edge) <- NULL
  
  unique_ID <- sample(10000 : 99999, 1)
  edge_file <- paste("edge_", unique_ID, ".txt", sep = "")
  output_file <- paste("output_", unique_ID, ".txt", sep = "")
  while (file.exists(edge_file)) {
    unique_ID <- sample(10000 : 99999, 1)
    edge_file <- paste("edge_", unique_ID, ".txt", sep = "")
    output_file <- paste("output", unique_ID, ".txt", sep = "")
  }
  
  
  write.table(x = edge, file = edge_file, sep = "\t", row.names = FALSE, 
              col.names = FALSE)
  if (modularity == 2 && resolution > 1){
    stop("error: resolution<1 for alternative modularity")
  }
  command <- paste("java -jar", ModularityJarFile, edge_file, output_file, 
                   modularity, resolution, algorithm, n.start, n.iter, 
                   random.seed, print.output, sep = " ")
  system(command, wait = TRUE)
  ident.use <- read.table(file = output_file, header = FALSE, sep = "\t")[, 1]
  object <- set.ident(object, object@cell.names, ident.use)
  file.remove(edge_file)
  file.remove(output_file)
  return (object)
}


GroupSingletons <- function(object, SNN){
  # identify singletons
  singletons <- c()
  for (cluster in unique(object@ident)) {
    if (length(which.cells(object, cluster)) == 1) {
      singletons <- append(singletons, cluster)
    }
  }
  # calculate connectivity of singletons to other clusters, add singleton
  # to cluster it is most connected to
  cluster.names <- unique(object@ident)
  cluster.names <- setdiff(cluster.names, singletons)
  while (length(singletons) > 0) {
    connectivity <- matrix(0, ncol = length(cluster.names), 
                           nrow = length(singletons))
    rownames(connectivity) <- singletons
    colnames(connectivity) <- cluster.names
    for (i in singletons) {
      for (j in cluster.names) {
        subSNN = SNN[which.cells(object, i), match(which.cells(object, j), colnames(SNN))]
        if (is.object(subSNN)) {
          connectivity[i, j] <- sum(subSNN) / (nrow(subSNN) * ncol(subSNN))
        } else {
          connectivity[i, j] <- mean(subSNN)
        }
      }
    }
    m <- max(connectivity, na.rm = T)
    mi <- which(connectivity == m, arr.ind = TRUE)
    c1 <- rownames(connectivity)[mi[1, 1]]
    c2 <- colnames(connectivity)[mi[1, 2]]
    object <- set.ident(object, cells.use = which.cells(object, c1), 
                        ident.use = c2)
    singletons <- singletons[singletons != c1]
  }
  return(object)
}