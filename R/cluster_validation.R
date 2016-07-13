#' @include seurat.R
NULL
#' Cluster Validation
#'
#' Methods for validating the legitimacy of clusters using
#' classification. SVMs are used as the basis for the classification.
#' Merging is done based on the connectivity from an SNN graph.
#'
#'
#' @param object Seurat object
#' @param pc.use Which PCs to use for model construction
#' @param top.genes Use the top X genes for model construction
#' @param min.connectivity Threshold of connectedness for comparison
#' of two clusters
#' @param acc.cutoff Accuracy cutoff for classifier
#' @param verbose Controls whether to display progress and merge results
#' @importFrom caret trainControl train
#' @return Returns a Seurat object, object@@ident has been updated with
#' new cluster info
#' @export
setGeneric("ValidateClusters", function(object, pc.use = NULL, top.genes = 30,
                                        min.connectivity = 0.01,
                                        acc.cutoff = 0.9, verbose = TRUE)
standardGeneric("ValidateClusters"))
#' @export
setMethod("ValidateClusters", signature = "seurat",
  function(object, pc.use = NULL, top.genes = 30, min.connectivity = 0.01,
           acc.cutoff = 0.9, verbose = TRUE){

  #probably should refactor to make cleaner
  if (length(object@snn.dense) > 1) {
    SNN.use <- object@snn.dense
  } else if (length(object@snn.sparse) > 1){
    SNN.use <- object@snn.sparse
  } else {
    stop("SNN matrix required. Please run BuildSNN() to save the SNN matrix in the object slot")
  }

  num.clusters.orig <- length(unique(object@ident))
  still_merging <- TRUE
  if (verbose) {
    connectivity <- CalcConnectivity(object, SNN.use)
    end <- length(connectivity[connectivity > min.connectivity])
    progress <- end
    status <- 0
  }
  # find connectedness of every two clusters
  while (still_merging) {
    connectivity <- CalcConnectivity(object, SNN.use)
    merge.done <- FALSE
    while (!merge.done) {
      m <- max(connectivity, na.rm = TRUE)
      mi <- which(connectivity == m, arr.ind = TRUE)
      c1 <- rownames(connectivity)[mi[, 1]]
      c2 <- rownames(connectivity)[mi[, 2]]
      if (m > min.connectivity) {
        acc <- RunClassifier(object, c1, c2, pc.use, top.genes)
        # if classifier can't classify them well enough, merge clusters
        if (acc < acc.cutoff) {
          object <- SetIdent(object, cells.use = which.cells(object, c1),
                              ident.use = c2)
          if (verbose) {
            progress <- length(connectivity[connectivity > min.connectivity])
            print(paste(sprintf("%3.0f", (1 - progress / end) * 100),
                  "% complete --- merge clusters ",c1, " and ", c2,
                  ", classification accuracy of ",
                  sprintf("%1.4f", acc), sep = ""))
          }
          merge.done <- TRUE
        } else {
          if (verbose & status == 5) {
            print(paste(sprintf("%3.0f", (1 - progress / end) * 100),
                  "% complete --- Last 5 cluster comparisons failed to merge",
                  " ,still checking possible merges ...", sep = ""))
            status <- 0
          }
          status <- status + 1
          connectivity[c1, c2] <- 0
          connectivity[c2, c1] <- 0
        }
      } else {
        still_merging <- FALSE
        break
      }
    }
  }
  if (verbose) {
    print(paste("100% complete --- started with ", num.clusters.orig,
          " clusters, ", length(unique(object@ident)), " clusters remaining",
          sep = "" ))
  }
  return(object)
})



#' Specific Cluster Validation
#'
#' Methods for validating the legitimacy of two specific clusters using
#' classification. SVMs are used as the basis for the classification.
#' Merging is done based on the connectivity from an SNN graph.
#'
#'
#' @param object Seurat object
#' @param cluster1 First cluster to check classification
#' @param cluster2 Second cluster to check with classification
#' @param pc.use Which PCs to use for model construction
#' @param top.genes Use the top X genes for model construction
#' @param acc.cutoff Accuracy cutoff for classifier
#' @importFrom caret trainControl train
#' @return Returns a Seurat object, object@@ident has been updated with
#' new cluster info
#' @export
setGeneric("ValidateSpecificClusters", function(object, cluster1 = NULL,
                                                  cluster2 = 1, pc.use=2,
                                                  top.genes = 30,
                                                  acc.cutoff = 0.9)
standardGeneric("ValidateSpecificClusters"))
#' @export
setMethod("ValidateSpecificClusters", signature = "seurat",
  function(object, cluster1 = NULL, cluster2 = 1, pc.use = 2, top.genes = 30,
           acc.cutoff = 0.9){
  acc <- RunClassifier(object, cluster1, cluster2, pc.use, top.genes)
  print(paste("Comparing cluster ", cluster1, " and ", cluster2, ": Acc = ",
        acc, sep = ""))
  if (acc < acc.cutoff) {
    object <- SetIdent(object, cells.use = which.cells(object, cluster1),
                        ident.use = cluster2)
    print(paste("merge cluster ", cluster1, " and ", cluster2))
    merge.done <- TRUE
  }
  return(object)
})


RunClassifier <- function(object, group1, group2, pcs, num.genes) {
  d1 <- which.cells(object, group1)
  d2 <- which.cells(object, group2)
  y  <- as.numeric(object@ident[c(d1, d2)]) - 1
  x  <- data.frame(t(object@data[pcTopGenes(object, pcs, num.genes),
                                            c(d1, d2)]));
  xv <- apply(x, 2, var)
  x  <- x[, names(xv > 0)]
  # run k-fold cross validation
  ctrl <- trainControl(method = "repeatedcv", repeats = 5)
  set.seed(1500)
  model <- train(as.factor(y)~., data = x, method = "svmLinear",
                 trControl = ctrl)
  acc <- model$results[, 2]
  return(acc)
}

