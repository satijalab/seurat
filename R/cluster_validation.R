#' @include seurat.R
NULL
#' Cluster Validation
#'
#' Methods for validating the legitimacy of clusters using classification. SVMs
#' are used as the basis for the classification. Merging is done based on the
#' connectivity from an SNN graph.
#'
#' @param object Seurat object
#' @param pc.use Which PCs to use to define genes in model construction
#' @param top.genes Use the top X genes for each PC in model construction
#' @param min.connectivity Threshold of connectedness for comparison of two
#' clusters
#' @param acc.cutoff Accuracy cutoff for classifier
#' @param verbose Controls whether to display progress and merging results
#' @importFrom caret trainControl train
#' @return Returns a Seurat object, object@@ident has been updated with new
#' cluster info
#' @export
#'
#' @examples
#' pbmc_small
#' # May throw warnings when cluster sizes are particularly small
#' \dontrun{
#' pbmc_small <- FindClusters(object = pbmc_small, reduction.type = "pca",
#'                            dims.use = 1:10, resolution = 1.1, save.SNN = TRUE)
#' pbmc_small <- ValidateClusters(pbmc_small, pc.use = 1:10)
#'}
#'
ValidateClusters <- function(
  object,
  pc.use = NULL,
  top.genes = 30,
  min.connectivity = 0.01,
  acc.cutoff = 0.9,
  verbose = TRUE
) {
  # probably should refactor to make cleaner
  if (length(x = object@snn) > 1) {
    SNN.use <- object@snn
  } else {
    stop("SNN matrix required. Please run BuildSNN() to save the SNN matrix in
         the object slot")
  }
  if (is.null(pc.use)){
    stop("pc.use not set. Please choose PCs.")
  }
  num.clusters.orig <- length(x = unique(x = object@ident))
  still_merging <- TRUE
  if (verbose) {
    connectivity <- CalcConnectivity(object = object)
    end <- length(x = connectivity[connectivity > min.connectivity])
    progress <- end
    status <- 0
  }
  # find connectedness of every two clusters
  while (still_merging) {
    connectivity <- CalcConnectivity(object = object)
    merge.done <- FALSE
    while (! merge.done) {
      m <- max(connectivity, na.rm = TRUE)
      mi <- which(x = connectivity == m, arr.ind = TRUE)
      c1 <- rownames(x = connectivity)[mi[, 1]]
      c2 <- rownames(x = connectivity)[mi[, 2]]
      if (m > min.connectivity) {
        acc <- RunClassifier(
          object = object,
          group1 = c1,
          group2 = c2,
          pcs = pc.use,
          num.genes = top.genes
        )
        # if classifier can't classify them well enough, merge clusters
        if (acc < acc.cutoff) {
          object <- SetIdent(
            object = object,
            cells.use = WhichCells(object = object, ident = c1),
            ident.use = c2
          )
          if (verbose) {
            progress <- length(x = connectivity[connectivity > min.connectivity])
            print(paste0(
              sprintf("%3.0f", (1 - progress / end) * 100),
              "% complete --- merge clusters ",
              c1,
              " and ",
              c2,
              ", classification accuracy of ",
              sprintf("%1.4f", acc)
            ))
          }
          merge.done <- TRUE
        } else {
          if (verbose & status == 5) {
            print(paste0(
              sprintf("%3.0f", (1 - progress / end) * 100),
              "% complete --- Last 5 cluster comparisons failed to merge, ",
              "still checking possible merges ..."
            ))
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
    print(paste0(
      "100% complete --- started with ",
      num.clusters.orig,
      " clusters, ",
      length(x = unique(x = object@ident)),
      " clusters remaining"
    ))
  }
  return(object)
}

#' Specific Cluster Validation
#'
#' Methods for validating the legitimacy of two specific clusters using
#' classification. SVMs are used as the basis for the classification.
#' Merging is done based on the connectivity from an SNN graph.
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
#'
#' @examples
#' \dontrun{
#' pbmc_small
#' pbmc_small <- FindClusters(object = pbmc_small, reduction.type = "pca",
#'                            dims.use = 1:10, resolution = 1.1, save.SNN = TRUE)
#' pbmc_small <- ValidateSpecificClusters(pbmc_small, cluster1 = 1,
#'                                        cluster2 = 2,  pc.use = 1:10)
#'}
#'
ValidateSpecificClusters <- function(
  object,
  cluster1 = NULL,
  cluster2 = 1,
  pc.use = 2,
  top.genes = 30,
  acc.cutoff = 0.9
) {
  acc <- RunClassifier(
    object = object,
    group1 = cluster1,
    group2 = cluster2,
    pcs = pc.use,
    num.genes = top.genes
  )
  print(paste0(
    "Comparing cluster ",
    cluster1,
    " and ",
    cluster2,
    ": Acc = ",
    acc
  ))
  if (acc < acc.cutoff) {
    object <- SetIdent(
      object = object,
      cells.use = WhichCells(object = object, ident = cluster1),
      ident.use = cluster2
    )
    print(paste("merge cluster", cluster1, "and", cluster2))
    merge.done <- TRUE
  }
  return(object)
}

# Train an SVM classifier and return the accuracy after 5 fold CV
#
# @param object     Seurat object
# @param group1     One identity to train classifier on
# @param group2     Second identity to train classifier on
# @param pcs        Vector of PCs on which to base genes to train classifier on.
#                   Pulls top num.genes genes associated with these PCs
# @param num.genes  Number of genes to pull for each PC
# @return           Returns the accuracy of the classifier after CV

RunClassifier <- function(object, group1, group2, pcs, num.genes) {
  d1 <- WhichCells(object = object, ident = group1)
  d2 <- WhichCells(object = object, ident = group2)
  y  <- as.numeric(x = object@ident[c(d1, d2)]) - 1
  x  <- data.frame(as.matrix(t(
    x = object@data[PCTopGenes(object = object, pc.use = pcs, num.genes =
                                 num.genes), c(d1, d2)]
    )))
  xv <- apply(X = x, MARGIN = 2, FUN = var)
  x  <- x[, names(x = which(xv > 0))]
  # run k-fold cross validation
  ctrl <- trainControl(method = "repeatedcv", repeats = 5)
  set.seed(seed = 1500)
  model <- train(
    x = x,
    y = as.factor(x = y),
    formula = as.factor(x = y) ~ .,
    method = "svmLinear",
    trControl = ctrl
  )
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
#' provided node
#' .
#' @return Returns the Out of Bag error for a random forest classifiers trained on
#' each internal node split or each split provided in the node list.
#'
#' @export
#'
#' @examples
#' pbmc_small
#' pbmc_small <- FindClusters(object = pbmc_small, reduction.type = "pca",
#'                            dims.use = 1:10, resolution = 1.1, save.SNN = TRUE)
#' pbmc_small <- BuildClusterTree(pbmc_small, reorder.numeric = TRUE, do.reorder = TRUE)
#' AssessNodes(pbmc_small)
#'
AssessNodes <- function(object, node.list, all.below = FALSE) {
  tree <- object@cluster.tree[[1]]
  if (missing(x = node.list)) {
    node.list <- GetAllInternalNodes(tree = tree)
  } else {
    possible.nodes <- GetAllInternalNodes(tree = tree)
    if (any(! node.list %in% possible.nodes)){
      stop(paste(
        node.list[!(node.list %in% possible.nodes)],
        "not valid internal nodes"
      ))
    }
    if (length(x = node.list == 1) && all.below) {
      node.list <- c(node.list, DFT(tree = tree, node = node.list))
    }
  }
  oobe <- pbsapply(
    X = node.list,
    FUN = function(x) {
      return(AssessSplit(
        object = object,
        node = x,
        print.output = FALSE,
        verbose = FALSE
      ))
    }
  )
  return(data.frame(node = node.list, oobe))
}

#' Assess Cluster Split
#'
#' Method for determining confidence in specific bifurcations in
#' the cluster tree. Use the Out of Bag (OOB) error of a random
#' forest classifier to judge confidence.
#'
#' @param object Seurat object
#' @param node Node in the cluster tree in question
#' @param cluster1 First cluster to compare
#' @param cluster2 Second cluster to compare
#' @param print.output Print the OOB error for the classifier
#' @inheritDotParams BuildRFClassifier -object
#' @return Returns the Out of Bag error for a random forest classifier
#' trained on the split from the given node
#' @export
#'
#' @examples
#' pbmc_small
#' pbmc_small <- FindClusters(object = pbmc_small, reduction.type = "pca",
#'                            dims.use = 1:10, resolution = 1.1, save.SNN = TRUE)
#' pbmc_small <- BuildClusterTree(pbmc_small, reorder.numeric = TRUE, do.reorder = TRUE)
#' # Assess based on a given node
#' AssessSplit(pbmc_small, node = 11)
#' # Asses based on two given clusters (or vectors of clusters)
#' AssessSplit(pbmc_small, cluster1 = 5, cluster2 = 6)
#' AssessSplit(pbmc_small, cluster1 = 4, cluster2 = c(5, 6))
#'
AssessSplit <- function(
  object,
  node,
  cluster1,
  cluster2,
  print.output = TRUE,
  ...
) {
  tree <- object@cluster.tree[[1]]
  if (! missing(x = node)){
    if (! missing(x = cluster1) || ! missing(x = cluster2)) {
      warning("Both node and cluster IDs provided. Defaulting to using node ID")
    }
    possible.nodes <- c(
      DFT(tree = tree, node = tree$edge[,1][1]),
      tree$edge[,1][1]
    )
    if (! node %in% possible.nodes) {
      stop("Not a valid node")
    }
    split <- tree$edge[which(x = tree$edge[,1] == node), ][,2]
    group1 <- DFT(tree = tree, node = split[1], only.children = TRUE)
    group2 <- DFT(tree = tree, node = split[2], only.children = TRUE)
    if (any(is.na(x = group1))) {
      group1 <- split[1]
    }
    if (any(is.na(x = group2))) {
      group2 <- split[2]
    }
  } else {
    group1 <- cluster1
    group2 <- cluster2
  }
  group1.cells <- WhichCells(object = object, ident = group1)
  group2.cells <- WhichCells(object = object, ident = group2)
  assess.data <- SubsetData(
    object = object,
    cells.use = c(group1.cells, group2.cells)
  )
  assess.data <- SetIdent(
    object = assess.data,
    cells.use = group1.cells,
    ident.use = "g1"
  )
  assess.data <- SetIdent(
    object = assess.data,
    cells.use = group2.cells,
    ident.use = "g2"
  )
  rfc <- BuildRFClassifier(
    object = assess.data,
    training.genes = assess.data@var.genes,
    training.classes = assess.data@ident,
    ...
  )
  oobe <- rfc$prediction.error
  if (print.output) {
    print(paste0("Out of Bag Error: ", round(x = oobe, digits = 4) * 100, "%"))
  }
  return(oobe)
}
