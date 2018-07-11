FindClusters.default <- function(
  object,
  modularity.fxn = 1,
  resolution = 0.8,
  algorithm = 1,
  n.start = 100,
  n.iter = 10,
  random.seed = 0,
  temp.file.location = NULL,
  edge.file.name = NULL,
  verbose = TRUE
){
  if(is.null(object)) {
    stop("Please provide an SNN graph")
  }
  clustering.results <- data.frame(row.names = colnames(object))
  for(r in resolution){
    ids <- RunModularityClustering(
      SNN = object,
      modularity = modularity.fxn,
      resolution = r,
      algorithm = algorithm,
      n.start = n.start,
      n.iter = n.iter,
      random.seed = random.seed,
      print.output = verbose,
      temp.file.location = temp.file.location,
      edge.file.name = edge.file.name)
    names(ids) <- colnames(object)
    ids <- GroupSingletons(ids = ids, SNN = object, verbose = verbose)
    clustering.results[, paste0("res.", r)] <- ids
  }
  return(clustering.results)
}

#' @param graph.name Name of graph to use for the clustering algorithm
#'
#' @describeIn FindClusters FindClusters on a Seurat object
#' @export
#' @method FindClusters Seurat
#'
FindClusters.Seurat <- function(
  object,
  graph.name = NULL,
  modularity.fxn = 1,
  resolution = 0.8,
  algorithm = 1,
  n.start = 100,
  n.iter = 10,
  random.seed = 0,
  temp.file.location = NULL,
  edge.file.name = NULL,
  verbose = TRUE
){
  graph.name <- graph.name %||% paste0(DefaultAssay(object), "_snn")
  if(! graph.name %in% names(object)) {
    stop("Provided graph.name not present in Seurat object")
  }
  if(! is(object = object[[graph.name]], class2 = "Graph")){
    stop("Provided graph.name does not correspond to a graph object.")
  }
  clustering.results <- FindClusters(
    object = object[[graph.name]],
    modularity.fxn = modularity.fxn,
    resolution = resolution,
    algorithm = algorithm,
    n.start = n.start,
    n.iter = n.iter,
    random.seed = random.seed,
    temp.file.location = temp.file.location,
    edge.file.name = edge.file.name,
    verbose = verbose
  )
  colnames(clustering.results) <- paste0(graph.name, "_", colnames(clustering.results))
  object <- AddMetaData(object = object, metadata = clustering.results)
  Idents(object = object) <- colnames(clustering.results)[ncol(clustering.results)]
  return(object)
}

#' Get Cluster Assignments
#'
#' Retrieve cluster IDs as a dataframe. First column will be the cell name,
#' second column will be the current cluster identity (pulled from object@ident).

#' @param object Seurat object with cluster assignments
#' @return Returns a dataframe with cell names and cluster assignments
#' @export
#'
#'@examples
#' pbmc_small
#' clusters <- GetClusters(object = pbmc_small)
#' head(clusters)
#'
GetClusters <- function(object) {
  clusters <- data.frame(cell.name = names(object@ident), cluster = object@ident)
  rownames(clusters) <- NULL
  clusters$cell.name <- as.character(clusters$cell.name)
  return(clusters)
}

#' Set Cluster Assignments
#'
#' Easily set the cluster assignments using the output of GetClusters() ---
#' a dataframe with cell names as the first column and cluster assignments as
#' the second.
#'
#' @param object Seurat object
#' @param clusters A dataframe containing the cell names and cluster assignments
#' to set for the object.
#' @return Returns a Seurat object with the identities set to the cluster
#' assignments that were passed.
#' @export
#'
#'@examples
#' pbmc_small
#' # Get clusters as a dataframe with GetClusters.
#' clusters <- GetClusters(object = pbmc_small)
#' # Use SetClusters to set cluster IDs
#' pbmc_small <- SetClusters(object = pbmc_small, clusters = clusters)
#'
SetClusters <- function(object, clusters = NULL) {
  if(!(all(c("cell.name", "cluster") %in% colnames(clusters)))){
    stop("The clusters parameter must be the output from GetClusters (i.e.
         Columns must be cell.name and cluster)")
  }
  cells.use <- clusters$cell.name
  ident.use <- clusters$cluster
  object <- SetIdent(
    object = object,
    cells.use = cells.use,
    ident.use = ident.use
  )
  return(object)
  }

#' Save cluster assignments to a TSV file
#'
#' @param object Seurat object with cluster assignments
#' @param file Path to file to write cluster assignments to
#'
#' @return No return value. Writes clusters assignments to specified file.
#'
#' @importFrom utils write.table
#'
#' @export
#'
#' @examples
#' \dontrun{
#' pbmc_small
#' file.loc <- "~/Desktop/cluster_assignments.tsv"
#' SaveClusters(object = pbmc_small, file = file.loc)
#' }
#'
SaveClusters <- function(object, file) {
  my.clusters <- GetClusters(object = object)
  write.table(my.clusters, file = file, sep="\t", quote = FALSE, row.names = F)
}

#' Convert the cluster labels to a numeric representation
#'
#' @param object Seurat object
#' @return Returns a Seurat object with the identities relabeled numerically
#' starting from 1.
#'
#' @export
#'
#' @examples
#' # Append "Cluster_" to cluster IDs to demonstrate numerical conversion
#' new.cluster.labels <- paste0("Cluster_", pbmc_small@ident)
#' pbmc_small <- SetIdent(
#'   object = pbmc_small,
#'   cells.use = pbmc_small@cell.names,
#'   ident.use = new.cluster.labels
#' )
#' unique(pbmc_small@ident)
#' # Now relabel the IDs numerically starting from 1
#' pbmc_small <- NumberClusters(pbmc_small)
#' unique(pbmc_small@ident)
#'
NumberClusters <- function(object) {
  clusters <- unique(x = object@ident)
  if(any(sapply(X = clusters,
                FUN = function(x) { !grepl("\\D", x) }))
     ) {
    n <- as.numeric(x = max(clusters)) + 1
    for (i in clusters) {
      object <- SetIdent(
        object = object,
        cells.use = WhichCells(object = object, ident = i),
        ident.use = n
      )
      n <- n + 1
    }
    clusters <- unique(x = object@ident)
  }
  n <- 1
  for (i in clusters) {
    object <- SetIdent(
      object,
      cells.use = WhichCells(object = object, ident = i),
      ident.use = n
    )
    n <- n + 1
  }
  return(object)
}

#' Classify New Data
#'
#' Classify new data based on the cluster information of the provided object.
#' Random Forests are used as the basis of the classification.
#'
#' @param object Seurat object on which to train the classifier
#' @param classifier Random Forest classifier from BuildRFClassifier. If not provided,
#' it will be built from the training data provided.
#' @param training.genes Vector of genes to build the classifier on
#' @param training.classes Vector of classes to build the classifier on
#' @param new.data New data to classify
#' @param ... additional parameters passed to ranger
#'
#' @return Vector of cluster ids
#'
#' @import Matrix
#' @importFrom stats predict
#'
#' @export
#'
#' @examples
#' pbmc_small
#' # take the first 10 cells as test data and train on the remaining 70 cells
#' test.pbmc <- SubsetData(object = pbmc_small, cells.use = pbmc_small@cell.names[1:10])
#' train.pbmc <- SubsetData(object = pbmc_small, cells.use = pbmc_small@cell.names[11:80])
#' predicted.classes <- ClassifyCells(
#'   object = train.pbmc,
#'   training.classes = train.pbmc@ident,
#'   new.data = test.pbmc@data
#' )
#'
ClassifyCells <- function(
  object,
  classifier,
  training.genes = NULL,
  training.classes = NULL,
  new.data = NULL,
  ...
) {
  PackageCheck('ranger')
  # build the classifier
  if (missing(classifier)){
    classifier <- BuildRFClassifier(
      object = object,
      training.genes = training.genes,
      training.classes = training.classes,
      ...
    )
  }
  # run the classifier on the new data
  features <- classifier$forest$independent.variable.names
  genes.to.add <- setdiff(x = features, y = rownames(x = new.data))
  data.to.add <- matrix(
    data = 0,
    nrow = length(x = genes.to.add),
    ncol = ncol(x = new.data)
  )
  rownames(x = data.to.add) <- genes.to.add
  new.data <- rbind(new.data, data.to.add)
  new.data <- new.data[features, ]
  new.data <- as.matrix(x = t(x = new.data))
  message("Running Classifier ...")
  prediction <- predict(classifier, new.data)
  new.classes <- prediction$predictions
  return(new.classes)
}

#' Build Random Forest Classifier
#'
#' Train the random forest classifier
#'
#'
#' @param object Seurat object on which to train the classifier
#' @param training.genes Vector of genes to build the classifier on
#' @param training.classes Vector of classes to build the classifier on
#' @param verbose Additional progress print statements
#' @param ... additional parameters passed to ranger
#'
#' @return Returns the random forest classifier
#'
#' @import Matrix
#'
#' @export
#'
#' @examples
#' pbmc_small
#' # Builds the random forest classifier to be used with ClassifyCells
#' # Useful if you want to use the same classifier with several sets of new data
#' classifier <- BuildRFClassifier(pbmc_small, training.classes = pbmc_small@ident)
#'
BuildRFClassifier <- function(
  object,
  training.genes = NULL,
  training.classes = NULL,
  verbose = TRUE,
  ...
) {
  PackageCheck('ranger')
  training.classes <- as.vector(x = training.classes)
  training.genes <- SetIfNull(
    x = training.genes,
    default = rownames(x = object@data)
  )
  training.data <- as.data.frame(
    x = as.matrix(
      x = t(
        x = object@data[training.genes, ]
      )
    )
  )
  training.data$class <- factor(x = training.classes)
  if (verbose) {
    message("Training Classifier ...")
  }
  classifier <- ranger::ranger(
    data = training.data,
    dependent.variable.name = "class",
    classification = TRUE,
    write.forest = TRUE,
    ...
  )
  return(classifier)
}

#' K-Means Clustering
#'
#' Perform k-means clustering on both genes and single cells
#'
#' K-means and heatmap are calculated on object@@scale.data
#'
#' @param object Seurat object
#' @param genes.use Genes to use for clustering
#' @param k.genes K value to use for clustering genes
#' @param k.cells K value to use for clustering cells (default is NULL, cells
#' are not clustered)
#' @param k.seed Random seed
#' @param do.plot Draw heatmap of clustered genes/cells (default is FALSE).
#' @param data.cut Clip all z-scores to have an absolute value below this.
#' Reduces the effect of huge outliers in the data.
#' @param k.cols Color palette for heatmap
#' @param set.ident If clustering cells (so k.cells>0), set the cell identity
#' class to its K-means cluster (default is TRUE)
#' @param do.constrained FALSE by default. If TRUE, use the constrained K-means function implemented in the tclust package.
#' @param assay.type Type of data to normalize for (default is RNA), but can be changed for multimodal analyses.
#' @param \dots Additional parameters passed to kmeans (or tkmeans)
#'
#' @importFrom methods new
#' @importFrom stats kmeans
#'
#' @return Seurat object where the k-means results for genes is stored in
#' object@@kmeans.obj[[1]], and the k-means results for cells is stored in
#' object@@kmeans.col[[1]]. The cluster for each cell is stored in object@@meta.data[,"kmeans.ident"]
#' and also object@@ident (if set.ident=TRUE)
#'
#' @export
#'
#' @examples
#' pbmc_small
#' # Cluster on genes only
#' pbmc_small <- DoKMeans(pbmc_small, k.genes = 3)
#' # Cluster on genes and cell
#' pbmc_small <- DoKMeans(pbmc_small, k.genes = 3, k.cells = 3)
#'
DoKMeans <- function(
  object,
  genes.use = NULL,
  k.genes = NULL,
  k.cells = 0,
  k.seed = 1,
  do.plot = FALSE,
  data.cut = 2.5,
  k.cols = PurpleAndYellow(),
  set.ident = TRUE,
  do.constrained = FALSE,
  assay.type="RNA",
  ...
) {
  data.use.orig <- GetAssayData(
    object = object,
    assay.type = assay.type,
    slot = "scale.data"
  )
  data.use <- MinMax(data = data.use.orig, min = data.cut * (-1), max = data.cut)
  genes.use <- SetIfNull(x = genes.use, default = object@var.genes)
  genes.use <- genes.use[genes.use %in% rownames(x = data.use)]
  cells.use <- object@cell.names
  kmeans.data <- data.use[genes.use, cells.use]
  if (do.constrained) {
    set.seed(seed = k.seed)
    PackageCheck('tclust')
    kmeans.obj <- tclust::tkmeans(x = kmeans.data, k = k.genes, ...)
  } else {
    set.seed(seed = k.seed)
    kmeans.obj <- kmeans(x = kmeans.data, centers = k.genes, ...)
  }
  names(x = kmeans.obj$cluster) <- genes.use
  #if we are going to k-means cluster cells in addition to genes
  kmeans.col <- c()
  if (k.cells > 0) {
    kmeans.col <- kmeans(x = t(x = kmeans.data), centers = k.cells)
    names(x = kmeans.col$cluster) <- cells.use
  }
  object.kmeans <- new(
    Class = "kmeans.info",
    gene.kmeans.obj = kmeans.obj,
    cell.kmeans.obj = kmeans.col
  )
  object@kmeans <- object.kmeans
  if (k.cells > 0) {
    kmeans.code=paste("kmeans",k.cells,"ident",sep=".")
    object@meta.data[names(x = kmeans.col$cluster), kmeans.code] <- kmeans.col$cluster
  }
  if (set.ident && (k.cells > 0)) {
    object <- SetIdent(
      object = object,
      cells.use = names(x = kmeans.col$cluster),
      ident.use = kmeans.col$cluster
    )
  }
  if (do.plot) {
    KMeansHeatmap(object = object)
  }
  return(object)
}

globalVariables(
  names = 'WeightedEuclideanDist',
  package = 'Seurat',
  add = TRUE
)
#' Phylogenetic Analysis of Identity Classes
#'
#' Constructs a phylogenetic tree relating the 'average' cell from each
#' identity class. Tree is estimated based on a distance matrix constructed in
#' either gene expression space or PCA space.
#'
#' Note that the tree is calculated for an 'average' cell, so gene expression
#' or PC scores are averaged across all cells in an identity class before the
#' tree is constructed.
#'
#' @param object Seurat object
#' @param genes.use Genes to use for the analysis. Default is the set of
#' variable genes (object@@var.genes). Assumes pcs.use=NULL (tree calculated in
#' gene expression space)
#' @param pcs.use If set, tree is calculated in PCA space.
#' @param SNN.use If SNN is passed, build tree based on SNN graph connectivity between clusters
#' @param do.plot Plot the resulting phylogenetic tree
#' @param do.reorder Re-order identity classes (factor ordering), according to
#' position on the tree. This groups similar classes together which can be
#' helpful, for example, when drawing violin plots.
#' @param reorder.numeric Re-order identity classes according to position on
#' the tree, assigning a numeric value ('1' is the leftmost node)
#' @param show.progress Show progress updates
#'
#' @return A Seurat object where the cluster tree is stored in
#' object@@cluster.tree[[1]]
#'
#' @importFrom ape as.phylo
#' @importFrom stats dist hclust
#'
#' @export
#'
#' @examples
#' pbmc_small
#' pbmc_small <- BuildClusterTree(pbmc_small, do.plot = FALSE)
#'
BuildClusterTree <- function(
  object,
  genes.use = NULL,
  pcs.use = NULL,
  SNN.use = NULL,
  do.plot = TRUE,
  do.reorder = FALSE,
  reorder.numeric = FALSE,
  show.progress = TRUE
) {
  genes.use <- SetIfNull(x = genes.use, default = object@var.genes)
  ident.names <- as.character(x = unique(x = object@ident))
  if (! is.null(x = genes.use)) {
    genes.use <- intersect(x = genes.use, y = rownames(x = object@data))
    data.avg <- AverageExpression(
      object = object,
      genes.use = genes.use,
      show.progress = show.progress
    )
    data.dist <- dist(t(x = data.avg[genes.use, ]))
  }
  if (! is.null(x = pcs.use)) {
    data.pca <- AveragePCA(object = object)
    data.dist <- dist(t(x = data.pca[pcs.use,]))
  }
  if (! is.null(x = SNN.use)) {
    num.clusters <- length(x = ident.names)
    data.dist <- matrix(data = 0, nrow = num.clusters, ncol = num.clusters)
    rownames(data.dist) <- ident.names
    colnames(data.dist) <- ident.names
    for (i in 1:(num.clusters - 1)) {
      for (j in (i + 1):num.clusters) {
        subSNN <- SNN.use[
          match(
            x = WhichCells(object = object, ident = ident.names[i]),
            table = colnames(x = SNN.use)
          ), # Row
          match(
            x = WhichCells(object = object, ident = ident.names[j]),
            table = rownames(x = SNN.use)
          ) # Column
          ]
        d <- mean(subSNN)
        if (is.na(x = d)) {
          data.dist[i, j] <- 0
        } else {
          data.dist[i, j] <- d
        }
      }
    }
    diag(x = data.dist) <- 1
    data.dist <- dist(data.dist)
  }
  data.tree <- as.phylo(x = hclust(d = data.dist))
  object@cluster.tree[[1]] <- data.tree
  if (do.reorder) {
    old.ident.order <- sort(x = unique(x = object@ident))
    data.tree <- object@cluster.tree[[1]]
    all.desc <- GetDescendants(tree = data.tree, node = (data.tree$Nnode + 2))
    all.desc <- old.ident.order[all.desc[all.desc <= (data.tree$Nnode + 1)]]
    object@ident <- factor(x = object@ident, levels = all.desc, ordered = TRUE)
    if (reorder.numeric) {
      object <- SetIdent(
        object = object,
        cells.use = object@cell.names,
        ident.use = as.integer(x = object@ident)
      )
      object@meta.data[object@cell.names, "tree.ident"] <- as.integer(x = object@ident)
    }
    object <- BuildClusterTree(
      object = object,
      genes.use = genes.use,
      pcs.use = pcs.use,
      do.plot = FALSE,
      do.reorder = FALSE,
      show.progress = show.progress
    )
  }
  if (do.plot) {
    PlotClusterTree(object)
  }
  return(object)
}


#' Perform spectral density clustering on single cells
#'
#' Find point clounds single cells in a two-dimensional space using density clustering (DBSCAN).
#'
#' @param object Seurat object
#' @param dim.1 First dimension to use
#' @param dim.2 second dimension to use
#' @param reduction.use Which dimensional reduction to use (either 'pca' or 'ica')
#' @param G.use Parameter for the density clustering. Lower value to get more fine-scale clustering
#' @param set.ident TRUE by default. Set identity class to the results of the density clustering.
#' Unassigned cells (cells that cannot be assigned a cluster) are placed in cluster 1, if there are any.
#' @param seed.use Random seed for the dbscan function
#' @param ... Additional arguments to be passed to the dbscan function
#'
#' @export
#'
#' @examples
#' pbmc_small
#' # Density based clustering on the first two tSNE dimensions
#' pbmc_small <- DBClustDimension(pbmc_small)
#'
DBClustDimension <- function(
  object,
  dim.1 = 1,
  dim.2 = 2,
  reduction.use = "tsne",
  G.use = NULL,
  set.ident = TRUE,
  seed.use = 1,
  ...
) {
  dim.code <- GetDimReduction(
    object = object,
    reduction.type = reduction.use,
    slot = 'key'
  )
  dim.codes <- paste0(dim.code, c(dim.1, dim.2))
  data.plot <- FetchData(object = object, vars.all = dim.codes)
  x1 <- paste0(dim.code, dim.1)
  x2 <- paste0(dim.code, dim.2)
  data.plot$x <- data.plot[, x1]
  data.plot$y <- data.plot[, x2]
  set.seed(seed = seed.use)
  data.mclust <- ds <- dbscan(data = data.plot[, c("x", "y")], eps = G.use, ...)
  to.set <- as.numeric(x = data.mclust$cluster + 1)
  data.names <- names(x = object@ident)
  object@meta.data[data.names, "DBclust.ident"] <- to.set
  if (set.ident) {
    object@ident <- factor(x = to.set)
    names(x = object@ident) <- data.names
  }
  return(object)
}

#' Perform spectral k-means clustering on single cells
#'
#' Find point clounds single cells in a low-dimensional space using k-means clustering.
#' Can be useful for smaller datasets, where graph-based clustering can perform poorly
#'
#' @param object A Seurat object
#' @param dims.use Dimensions to use for clustering
#' @param reduction.use Dimmensional Reduction to use for k-means clustering
#' @param k.use Number of clusters
#' @param set.ident Set identity of Seurat object
#' @param seed.use Random seed to use
#'
#' @return Object with clustering information
#'
#' @importFrom stats kmeans
#'
#' @export
#'
#' @examples
#' pbmc_small
#' # K-means clustering on the first two tSNE dimensions
#' pbmc_small <- KClustDimension(pbmc_small)
#'
KClustDimension <- function(
  object,
  dims.use = c(1,2),
  reduction.use = "tsne",
  k.use = 5,
  set.ident = TRUE,
  seed.use = 1
) {
  dim.code <- GetDimReduction(
    object = object,
    reduction.type = reduction.use,
    slot = 'key'
  )
  dim.codes <- paste0(dim.code, dims.use)
  data.plot <- FetchData(object = object, vars.all = dim.codes)
  set.seed(seed = seed.use)
  data.mclust <- ds <- kmeans(x = data.plot, centers = k.use)
  to.set <- as.numeric(x = data.mclust$cluster)
  data.names <- names(x = object@ident)
  object@meta.data[data.names, "kdimension.ident"] <- to.set
  if (set.ident) {
    object@ident <- factor(x = to.set)
    names(x = object@ident) <- data.names
  }
  return(object)
}
