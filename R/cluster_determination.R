#' @include seurat.R
NULL
#' Cluster Determination
#'
#' Identify clusters of cells by a shared nearest neighbor (SNN) modularity optimization
#' based clustering algorithm. First calculate k-nearest neighbors and construct
#' the SNN graph. Then optimize the modularity function to determine clusters.
#' For a full description of the algorithms, see Waltman and van Eck (2013) 
#' \emph{The European Physical Journal B}.
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
#' @param save.SNN Whether to save the SNN in an object slot
#' @param reuse.SNN Force utilization of stored SNN. If none store, this will
#' throw an error.
#' @param do.sparse Option to store and use SNN matrix as a sparse matrix.
#'        May be necessary datasets containing a large number of cells.
#' @param modularity.fxn Modularity function (1 = standard; 2 = alternative).
#' @param resolution Value of the resolution parameter, use a value above
#'        (below) 1.0 if you want to obtain a larger (smaller) number of
#'        communities.
#' @param algorithm Algorithm for modularity optimization (1 = original Louvain
#'        algorithm; 2 = Louvain algorithm with multilevel refinement;
#'        3 = SLM algorithm).
#' @param n.start Number of random starts.
#' @param n.iter Maximal number of iterations per random start.
#' @param random.seed Seed of the random number generator.
#' @param print.output Whether or not to print output to the console
#' @importFrom FNN get.knn
#' @importFrom igraph plot.igraph graph.adjlist
#' @importFrom Matrix sparseMatrix
#' @return Returns a Seurat object and optionally the SNN matrix,
#'         object@@ident has been updated with new cluster info
#' @export
setGeneric("FindClusters", function(object, genes.use = NULL, pc.use = NULL,
                                     k.param = 30, k.scale = 25,
                                     plot.SNN = FALSE, prune.SNN = 1/15,
                                     save.SNN = FALSE, reuse.SNN = FALSE,
                                     do.sparse = FALSE, modularity.fxn = 1, 
                                     resolution = 0.8, algorithm = 1, 
                                     n.start = 100, n.iter = 10, random.seed = 0,
                                     print.output = TRUE)
standardGeneric("FindClusters"))
#' @export
setMethod("FindClusters", signature = "seurat",
          function(object, genes.use = NULL, pc.use = NULL, k.param = 30,
                   k.scale = 25, plot.SNN = FALSE, prune.SNN = 1/15,
                   save.SNN = FALSE, reuse.SNN = FALSE, do.sparse = FALSE, 
                   modularity.fxn = 1, resolution = 0.8, algorithm = 1,
                   n.start = 100, n.iter = 10, random.seed = 0, print.output = TRUE){

  # for older objects without the snn.k slot
  if(typeof(validObject(object, test = T)) == "character"){
    object@snn.k <- numeric()
  }
  
  snn.built <- FALSE
  if (.hasSlot(object, "snn.dense")) {
    if (length(object@snn.dense) > 1) {
      snn.built <- TRUE
    }
  }
  if (.hasSlot(object, "snn.sparse")) {
    if (length(object@snn.sparse) > 1) {
      snn.built <- TRUE
    }
  }
  
  if((missing(genes.use) && missing(pc.use) && missing(k.param) && missing(k.scale) && 
     missing(prune.SNN) && snn.built) || reuse.SNN){
    save.SNN <- TRUE
    if (reuse.SNN && !snn.built){
      stop("No SNN stored to reuse.")
    }
    if (reuse.SNN && (!missing(genes.use) || !missing(pc.use) || !missing(k.param) || 
                      !missing(k.scale) || !missing(prune.SNN))){
      warning("SNN was not be rebuilt with new parameters. Continued with stored SNN. To suppress this
              warning, remove all SNN building parameters.")
    }
  }
  # if any SNN building parameters are provided or it hasn't been built, build a new SNN
  else{
    object <- BuildSNN(object, genes.use, pc.use, k.param, k.scale,
                       plot.SNN, prune.SNN, do.sparse, print.output)
  }
  
  # deal with sparse SNNs
  if (length(object@snn.sparse) > 1) {
    SNN.use <- object@snn.sparse
  } else {
    SNN.use <- object@snn.dense
  }
  for (r in resolution) {
    object <- RunModularityClustering(object, SNN.use, modularity.fxn, r,
                                      algorithm, n.start, n.iter, random.seed,
                                      print.output)
    object <- GroupSingletons(object, SNN.use)
    name <- paste("res.", r, sep = "")
    object <- StashIdent(object, name)
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
            object <- SetIdent(object, cells.use, ident.use)
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
                object <- SetIdent(object, cells.use = WhichCells(object, i),
                                    ident.use = n)
                n <- n + 1
              }
              clusters <- unique(object@ident)
            }
            n <- 1
            for (i in clusters) {
              object <- SetIdent(object, cells.use = WhichCells(object, i),
                                  ident.use = n)
              n <- n + 1
            }
            return(object)
          }
)


RunModularityClustering <- function(object, SNN = matrix(), modularity = 1,
                                    resolution = 0.8, algorithm = 1,
                                    n.start = 100, n.iter = 10, random.seed = 0,
                                    print.output = TRUE){

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
  if (print.output) {
    print.output <- 1
  }
  else {
    print.output <- 0
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

  object <- SetIdent(object, object@cell.names, ident.use)
  file.remove(edge_file)
  file.remove(output_file)
  return (object)
}

GroupSingletons <- function(object, SNN){
  # identify singletons
  singletons <- c()
  for (cluster in unique(object@ident)) {
    if (length(WhichCells(object, cluster)) == 1) {
      singletons <- append(singletons, cluster)
    }
  }
  # calculate connectivity of singletons to other clusters, add singleton
  # to cluster it is most connected to
  cluster_names <- unique(object@ident)
  cluster_names <- setdiff(cluster_names, singletons)
  connectivity <- vector(mode="numeric", length = length(cluster_names))
  names(connectivity) <- cluster_names
  for (i in singletons) {
    for (j in cluster_names) {
      subSNN = SNN[WhichCells(object, i), match(WhichCells(object, j), colnames(SNN))]
      if (is.object(subSNN)) {
        connectivity[j] <- sum(subSNN) / (nrow(subSNN) * ncol(subSNN))
      } else {
        connectivity[j] <- mean(subSNN)
      }
    }
    m <- max(connectivity, na.rm = T)
    mi <- which(connectivity == m, arr.ind = TRUE)
    closest_cluster <- sample(names(connectivity[mi]), 1)
    object <- SetIdent(object, cells.use = WhichCells(object,i), 
                        ident.use = closest_cluster)
    
  }
  if (length(singletons) > 0){
    print(paste(length(singletons), "singletons identified.", length(unique(object@ident)), "final clusters."))
  }
  return(object)
}