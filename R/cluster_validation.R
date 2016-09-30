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
          object <- SetIdent(object, cells.use = WhichCells(object, c1),
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
    object <- SetIdent(object, cells.use = WhichCells(object, cluster1),
                        ident.use = cluster2)
    print(paste("merge cluster ", cluster1, " and ", cluster2))
    merge.done <- TRUE
  }
  return(object)
})


RunClassifier <- function(object, group1, group2, pcs, num.genes) {
  d1 <- WhichCells(object, group1)
  d2 <- WhichCells(object, group2)
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


#' Assess Internal Nodes
#' 
#' Method for automating assessment of tree splits over all internal nodes,
#' or a provided list of internal nodes. Uses AssessSplit() for calculation
#' of Out of Bag error (proxy for confidence in split).
#' 
#' @param object Seurat object
#' @param node.list List of internal nodes to assess and return 
#' @param all.below If single node provided in node.list, assess all splits below (and including)
#' provided node.
#' @return Returns the Out of Bag error for a random forest classifiers trained on 
#' each internal node split or each split provided in the node list.
#' @export 
setGeneric("AssessNodes", function(object, node.list, all.below = FALSE) standardGeneric("AssessNodes"))
#' @export
setMethod("AssessNodes", signature = "seurat", 
          function(object, node.list, all.below = FALSE){
            tree <- object@cluster.tree[[1]]
            if(missing(node.list)){
              node.list <- GetAllInternalNodes(tree)
            }
            else{
              possible.nodes <- GetAllInternalNodes(tree)
              if(any(!node.list %in% possible.nodes)){
                stop(paste0(node.list[!(node.list %in% possible.nodes)], "not valid internal nodes", sep = ""))
              }
              if(length(node.list == 1) && all.below){
                node.list <- c(node.list, DFT(tree, node.list))
              }
            }
          oobe <- pbsapply(node.list, function(x) AssessSplit(object, node = x, print.output = F, verbose = F))
          return(data.frame(node = node.list, oobe))
          })


#' Assess Cluster Split
#'
#' Method for determining confidence in specific bifurcations in
#' the cluster tree. Use the Out of Bag (OOB) error of a random
#' forest classifier to judge confidence.
#' 
#'
#' @param object Seurat object
#' @param node Node in the cluster tree in question
#' @param cluster1 First cluster to compare
#' @param cluster2 Second cluster to compare 
#' @param print.output Print the OOB error for the classifier
#' @param ... additional parameters to pass to BuildRFClassifier
#' @return Returns the Out of Bag error for a random forest classifier 
#' trained on the split from the given node
#' @export
setGeneric("AssessSplit", function(object, node, cluster1, cluster2, print.output = T, ...)
  standardGeneric("AssessSplit"))
#' @export
setMethod("AssessSplit", signature = "seurat",
          function(object, node, cluster1, cluster2, print.output = T, ...){
            tree <- object@cluster.tree[[1]]
            if(!missing(node)){
              if(!missing(cluster1) || !missing(cluster2)){
                warning("Both node and cluster IDs provided. Defaulting to using node ID")
              }
              possible.nodes <- c(DFT(tree, tree$edge[,1][1]), tree$edge[,1][1])
              if(! node %in% possible.nodes){
                stop("Not a valid node")
              }
              split <- tree$edge[which(tree$edge[,1] == node), ][,2]
              group1 <- DFT(tree, node = split[1], only.children = T)
              group2 <- DFT(tree, node = split[2], only.children = T)
              if(any(is.na(group1))) group1 <- split[1]
              if(any(is.na(group2))) group2 <- split[2]
            }
            else{
              group1 <- cluster1
              group2 <- cluster2
            }
            group1.cells <- WhichCells(object, ident = group1)
            group2.cells <- WhichCells(object, ident = group2)
            assess.data <- SubsetData(object, cells.use = c(group1.cells, group2.cells))
            assess.data <- SetIdent(assess.data, cells.use = group1.cells, ident.use = "g1")
            assess.data <- SetIdent(assess.data, cells.use = group2.cells, ident.use = "g2")
            rfc <- BuildRFClassifier(object = assess.data, training.genes = assess.data@var.genes, training.classes = assess.data@ident, ...)
            oobe <- rfc$prediction.error
            if(print.output){
              print(paste("Out of Bag Error: ", round(oobe, 4) * 100, "%", sep= ""))
            }
            return(oobe)
          }
)

#' Color tSNE Plot Based on Split
#' 
#' Returns a tSNE plot colored based on whether the cells fall in clusters
#' to the left or to the right of a node split in the cluster tree.
#' 
#' @param object Seurat object
#' @param node Node in cluster tree on which to base the split
#' @param color1 Color for the left side of the split
#' @param color2 Color for the right side of the split
#' @param color3 Color for all other cells
#' @param ... Additional parameters for TSNEPlot()
#' @return Returns a tSNE plot
#' @export
setGeneric("ColorTSNESplit", function(object, node, color1 = "red", color2 = "blue", color3 = "gray")
  standardGeneric("ColorTSNESplit"))
#' @export
setMethod("ColorTSNESplit", signature = "seurat",
          function(object, node, color1 = "red", color2 = "blue", color3 = "gray"){
            tree <- object@cluster.tree[[1]]
            split <- tree$edge[which(tree$edge[,1] == node), ][,2]
            all.children <- DFT(tree, node = tree$edge[,1][1], only.children = T)
            left.group <- DFT(tree, node = split[1], only.children = T)
            right.group <- DFT(tree, node = split[2], only.children = T)
            if(any(is.na(left.group))) left.group <- split[1]
            if(any(is.na(right.group))) right.group <- split[2]

            remaining.group <- setdiff(all.children, c(left.group, right.group))
            left.cells <- WhichCells(object, left.group)
            right.cells <- WhichCells(object, right.group)
            remaining.cells <- WhichCells(object, remaining.group)
            object <- SetIdent(object, cells.use = left.cells, ident.use = "Left Split")
            object <- SetIdent(object, cells.use = right.cells, ident.use = "Right Split")
            object <- SetIdent(object, cells.use = remaining.cells, ident.use = "Not in Split")
            colors.use = c(color1, color3, color2)
          
            return(TSNEPlot(object, colors.use = colors.use))
          }
)