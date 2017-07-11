#' @include seurat.R
NULL
#' Cluster Determination
#'
#' Identify clusters of cells by a shared nearest neighbor (SNN) modularity
#' optimization based clustering algorithm. First calculate k-nearest neighbors
#' and construct the SNN graph. Then optimize the modularity function to
#' determine clusters. For a full description of the algorithms, see Waltman and
#' van Eck (2013) \emph{The European Physical Journal B}.
#'
#' @param object Seurat object
#' @param genes.use A vector of gene names to use in construction of SNN graph
#' if building directly based on expression data rather than a dimensionally
#' reduced representation (i.e. PCs).
#' @param reduction.type Name of dimensional reduction technique to use in
#' construction of SNN graph. (e.g. "pca", "ica")
#' @param dims.use A vector of the dimensions to use in construction of the SNN
#' graph (e.g. To use the first 10 PCs, pass 1:10)
#' @param k.param Defines k for the k-nearest neighbor algorithm
#' @param k.scale Granularity option for k.param
#' @param plot.SNN Plot the SNN graph
#' @param prune.SNN Sets the cutoff for acceptable Jaccard distances when
#' computing the neighborhood overlap for the SNN construction. Any edges with
#' values less than or equal to this will be set to 0 and removed from the SNN
#' graph. Essentially sets the strigency of pruning (0 --- no pruning, 1 ---
#' prune everything).
#' @param print.output Whether or not to print output to the console
#' @param distance.matrix Build SNN from distance matrix (experimental)
#' @param save.SNN Saves the SNN matrix associated with the calculation in
#' object@@snn
#' @param reuse.SNN Force utilization of stored SNN. If none store, this will
#' throw an error.
#' @param modularity.fxn Modularity function (1 = standard; 2 = alternative).
#' @param resolution Value of the resolution parameter, use a value above
#' (below) 1.0 if you want to obtain a larger (smaller) number of communities.
#' @param algorithm Algorithm for modularity optimization (1 = original Louvain
#' algorithm; 2 = Louvain algorithm with multilevel refinement; 3 = SLM
#' algorithm).
#' @param n.start Number of random starts.
#' @param n.iter Maximal number of iterations per random start.
#' @param random.seed Seed of the random number generator.
#' @param temp.file.location Directory where intermediate files will be written.
#' Specify the ABSOLUTE path.
#' @importFrom FNN get.knn
#' @importFrom igraph plot.igraph graph.adjlist
#' @importFrom Matrix sparseMatrix
#' @return Returns a Seurat object and optionally the SNN matrix,
#'         object@@ident has been updated with new cluster info
#'
#' @export
#'
FindClusters <- function(
  object,
  genes.use = NULL,
  reduction.type = "pca",
  dims.use = NULL,
  k.param = 30,
  k.scale = 25,
  plot.SNN = FALSE,
  prune.SNN = 1/15,
  print.output = TRUE,
  distance.matrix = NULL,
  save.SNN = FALSE,
  reuse.SNN = FALSE,
  modularity.fxn = 1,
  resolution = 0.8,
  algorithm = 1,
  n.start = 100,
  n.iter = 10,
  random.seed = 0,
  temp.file.location = NULL
) {
  # for older objects without the snn.k slot
  if(typeof(x = validObject(object = object, test = TRUE)) == "character") {
    object@snn.k <- numeric()
  }
  snn.built <- FALSE
  if (.hasSlot(object = object, name = "snn")) {
    if (length(x = object@snn) > 1) {
      snn.built <- TRUE
    }
  }
  if ((
    missing(x = genes.use) && missing(x = dims.use) && missing(x = k.param) &&
    missing(x = k.scale) && missing(x = prune.SNN) && snn.built
  ) || reuse.SNN) {
    save.SNN <- TRUE
    if (reuse.SNN && !snn.built) {
      stop("No SNN stored to reuse.")
    }
    if (reuse.SNN && (
      ! missing(x = genes.use) || ! missing(x = dims.use) || ! missing(x = k.param)
      || ! missing(x = k.scale) || ! missing(x = prune.SNN)
    )) {
      warning("SNN was not be rebuilt with new parameters. Continued with stored
               SNN. To suppress this warning, remove all SNN building parameters.")
    }
  } else {
    # if any SNN building parameters are provided or it hasn't been built, build
    # a new SNN
    object <- BuildSNN(
      object = object,
      genes.use = genes.use,
      reduction.type = reduction.type,
      dims.use = dims.use,
      k.param = k.param,
      k.scale = k.scale,
      plot.SNN = plot.SNN,
      prune.SNN = prune.SNN,
      print.output = print.output,
      distance.matrix = distance.matrix
    )
  }
  for (r in resolution) {
    object <- RunModularityClustering(
      object = object,
      SNN = object@snn,
      modularity = modularity.fxn,
      resolution = r,
      algorithm = algorithm,
      n.start = n.start,
      n.iter = n.iter,
      random.seed = random.seed,
      print.output = print.output,
      temp.file.location = temp.file.location
    )
    object <- GroupSingletons(object = object, SNN = object@snn)
    name <- paste0("res.", r)
    object <- StashIdent(object = object, save.name = name)
  }
  if (!save.SNN) {
    object@snn <- sparseMatrix(1, 1, x = 1)
    object@snn.k <- integer()
  }
  return(object)
}

# Runs the modularity optimizer java program (ModularityOptimizer.jar)
#
#
# @param object               Seurat object
# @param SNN SNN              matrix to use as input for the clustering
#                             algorithms
# @param modularity           Modularity function to use in clustering (1 =
#                             standard; 2 = alternative).
# @param resolution           Value of the resolution parameter, use a value
#                             above (below) 1.0 if you want to obtain a larger
#                             (smaller) number of communities.
# @param algorithm            Algorithm for modularity optimization (1 =
#                             original Louvain algorithm; 2 = Louvain algorithm
#                             with multilevel refinement; 3 = SLM algorithm)
# @param n.start              Number of random starts.
# @param n.iter               Maximal number of iterations per random start.
# @param random.seed          Seed of the random number generator
# @param print.output         Whether or not to print output to the console
# @param temp.file.location   Directory where intermediate files will be written.
# @return                     Seurat object with identities set to the results
#                             of the clustering procedure.

RunModularityClustering <- function(
  object,
  SNN = matrix(),
  modularity = 1,
  resolution = 0.8,
  algorithm = 1,
  n.start = 100,
  n.iter = 10,
  random.seed = 0,
  print.output = TRUE,
  temp.file.location = NULL
) {
  seurat.dir <- system.file(package = "Seurat")
  ModularityJarFile <- paste0(seurat.dir, "/java/ModularityOptimizer.jar")
  seurat.dir.base <- strsplit(x = seurat.dir, split = "/")[[1]]
  seurat.dir <- paste0(
    seurat.dir.base[0:(length(x = seurat.dir.base) - 1)],
    collapse = "/"
  )
  seurat.dir <- paste0(seurat.dir, "/")
  diag(x = SNN) <- 0
  if (is.object(x = SNN)) {
    SNN <- as(object = SNN, Class = "dgTMatrix")
    edge <- cbind(i = SNN@j, j = SNN@i, x = SNN@x)
  } else {
    swap <- which(x = SNN != 0, arr.ind = TRUE) - 1
    temp <- swap[, 1]
    swap[, 1] <- swap[, 2]
    swap[, 2] <- temp
    edge <- cbind(swap, SNN[which(x = SNN != 0, arr.ind = TRUE)])
  }
  rownames(x = edge) <- NULL
  colnames(x = edge) <- NULL
  edge <- edge[! duplicated(x = edge[, 1:2]), ]
  temp.file.location <- SetIfNull(x = temp.file.location, default = seurat.dir)
  unique_ID <- sample(x = 10000:99999, size = 1)
  edge_file <- paste0(temp.file.location, "edge_", unique_ID, ".txt")
  output_file <- paste0(temp.file.location, "output_", unique_ID, ".txt")
  while (file.exists(edge_file)) {
    unique_ID <- sample(x = 10000:99999, size = 1)
    edge_file <- paste0(temp.file.location, "edge_", unique_ID, ".txt")
    output_file <- paste0(temp.file.location, "output", unique_ID, ".txt")
  }
  if (print.output) {
    print.output <- 1
  } else {
    print.output <- 0
  }
  write.table(
    x = edge,
    file = edge_file,
    sep = "\t",
    row.names = FALSE,
    col.names = FALSE
  )
  if (modularity == 2 && resolution > 1) {
    stop("error: resolution<1 for alternative modularity")
  }
  command <- paste(
    "java -jar",
    shQuote(string = ModularityJarFile),
    shQuote(string = edge_file),
    shQuote(string = output_file),
    modularity,
    resolution,
    algorithm,
    n.start,
    n.iter,
    random.seed,
    print.output
  )
  system(command, wait = TRUE)
  ident.use <- read.table(file = output_file, header = FALSE, sep = "\t")[, 1]

  object <- SetIdent(
    object = object,
    cells.use = object@cell.names,
    ident.use = ident.use
  )
  file.remove(edge_file)
  file.remove(output_file)
  return (object)
}

# Group single cells that make up their own cluster in with the cluster they are
# most connected to.
#
# @param object  Seurat object
# @param SNN     SNN graph used in clustering
# @return        Returns Seurat object with all singletons merged with most
#                connected cluster

GroupSingletons <- function(object, SNN) {
  # identify singletons
  singletons <- c()
  for (cluster in unique(x = object@ident)) {
    if (length(x = WhichCells(object = object, ident = cluster)) == 1) {
      singletons <- append(x = singletons, values = cluster)
    }
  }
  # calculate connectivity of singletons to other clusters, add singleton
  # to cluster it is most connected to
  cluster_names <- unique(x = object@ident)
  cluster_names <- setdiff(x = cluster_names, y = singletons)
  connectivity <- vector(mode="numeric", length = length(x = cluster_names))
  names(x = connectivity) <- cluster_names
  for (i in singletons) {
    for (j in cluster_names) {
      subSNN = SNN[
        WhichCells(object = object, ident = i), # Row
        match(
          x = WhichCells(object = object, ident = j),
          table = colnames(x = SNN)
        )
      ]
      if (is.object(x = subSNN)) {
        connectivity[j] <- sum(subSNN) / (nrow(x = subSNN) * ncol(x = subSNN))
      } else {
        connectivity[j] <- mean(x = subSNN)
      }
    }
    m <- max(connectivity, na.rm = T)
    mi <- which(x = connectivity == m, arr.ind = TRUE)
    closest_cluster <- sample(x = names(x = connectivity[mi]), 1)
    object <- SetIdent(
      object = object,
      cells.use = WhichCells(object = object, ident = i),
      ident.use = closest_cluster
    )
  }
  if (length(x = singletons) > 0) {
    print(paste(
      length(x = singletons),
      "singletons identified.",
      length(x = unique(object@ident)),
      "final clusters."
    ))
  }
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
#' @return No return value. Writes clusters assignments to specified file.
#' @export
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
NumberClusters <- function(object) {
  clusters <- unique(x = object@ident)
  if (typeof(x = clusters) == "integer") {
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
#' @importFrom ranger ranger
#'
#' @export
#'
ClassifyCells <- function(
  object,
  classifier,
  training.genes = NULL,
  training.classes = NULL,
  new.data = NULL,
  ...
) {
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
  print("Running Classifier ...")
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
#' @importFrom ranger ranger
#'
#' @export
#'
BuildRFClassifier <- function(
  object,
  training.genes = NULL,
  training.classes = NULL,
  verbose = TRUE,
  ...
) {
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
    print("Training Classifier ...")
  }
  classifier <- ranger(
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
#' Perform k=means clustering on both genes and single cells
#'
#' K-means and heatmap are calculated on object@@scale.data
#'
#' @param object Seurat object
#' @param genes.use Genes to use for clustering
#' @param k.genes K value to use for clustering genes
#' @param k.cells K value to use for clustering cells (default is NULL, cells
#' are not clustered)
#' @param k.seed Random seed
#' @param do.plot Draw heatmap of clustered genes/cells (default is TRUE)
#' @param data.cut Clip all z-scores to have an absolute value below this.
#' Reduces the effect of huge outliers in the data.
#' @param k.cols Color palette for heatmap
#' @param pc.row.order Order gene clusters based on the average PC score within
#' a cluster. Can be useful if you want to visualize clusters, for example,
#' based on their average score for PC1.
#' @param pc.col.order Order cell clusters based on the average PC score within
#' a cluster
#' @param rev.pc.order Use the reverse PC ordering for gene and cell clusters
#' (since the sign of a PC is arbitrary)
#' @param use.imputed Cluster imputed values (default is FALSE)
#' @param set.ident If clustering cells (so k.cells>0), set the cell identity
#' class to its K-means cluster (default is TRUE)
#' @param do.constrained FALSE by default. If TRUE, use the constrained K-means function implemented in the tclust package.
#' @param \dots Additional parameters passed to kmeans (or tkmeans)
#'
#' @importFrom tclust tkmeans
#'
#' @return Seurat object where the k-means results for genes is stored in
#'
#' object@@kmeans.obj[[1]], and the k-means results for cells is stored in
#' object@@kmeans.col[[1]]. The cluster for each cell is stored in object@@data.info[,"kmeans.ident"]
#' and also object@@ident (if set.ident=TRUE)
#'
#' @export
DoKMeans <- function(
  object,
  genes.use = NULL,
  k.genes = NULL,
  k.cells = 0,
  k.seed = 1,
  do.plot = TRUE,
  data.cut = 2.5,
  k.cols = pyCols,
  pc.row.order = NULL,
  pc.col.order = NULL,
  rev.pc.order = FALSE,
  use.imputed = FALSE,
  set.ident = TRUE,
  do.constrained = FALSE,
  ...
) {
  if (use.imputed) {
    data.use.orig <- data.frame(t(x = scale(x = t(x = object@imputed))))
  } else {
    data.use.orig <- object@scale.data
  }
  data.use <- minmax(data = data.use.orig, min = data.cut * (-1), max = data.cut)
  if (rev.pc.order) {
    revFxn <- function(x) {
      max(x) + 1 - x
    }
  } else {
    revFxn <- same
  }
  kmeans.col <- NULL
  genes.use <- SetIfNull(x = genes.use, default = object@var.genes)
  genes.use <- genes.use[genes.use %in% rownames(x = data.use)]
  cells.use <- object@cell.names
  kmeans.data <- data.use[genes.use, cells.use]
  if (do.constrained) {
    set.seed(seed = k.seed)
    kmeans.obj <- tkmeans(objectkmeans.data, k = k.genes, ...)
  } else {
    set.seed(seed = k.seed)
    kmeans.obj <- kmeans(x = kmeans.data, centers = k.genes, ...)
  }
  if (! is.null(x = pc.row.order)) {
    if (nrow(x = object@pca.x.full > 0)) {
      pcx.use <- object@pca.x.full
      pc.genes <- ainb(a = genes.use, b = rownames(x = pcx.use))
    } else {
      pcx.use <- object@pca.x
      pc.genes <- ainb(a = genes.use, b = rownames(x = pcx.use))
    }
    kmeans.obj$cluster <- as.numeric(
      x = revFxn(
        rank(
          x = tapply(
            X = pcx.use[genes.use, pc.row.order],
            INDEX = as.numeric(x = kmeans.obj$cluster),
            FUN = mean
          )
        )
      )[as.numeric(x = kmeans.obj$cluster)]
    )
  }
  names(x = kmeans.obj$cluster) <- genes.use
  if (k.cells > 0) {
    kmeans.col <- kmeans(x = t(x = kmeans.data), centers = k.cells)
    if (! is.null(x = pc.col.order)) {
      kmeans.col$cluster <- as.numeric(
        x = revFxn(
          rank(
            x = tapply(
              X = object@pca.rot[cells.use, pc.col.order],
              INDEX = kmeans.col$cluster,
              FUN = mean
            )
          )
        )[as.numeric(x = kmeans.col$cluster)]
      )
    }
    names(x = kmeans.col$cluster) <- cells.use
  }
  object@kmeans.obj <- list(kmeans.obj)
  object@kmeans.col <- list(kmeans.col)
  kmeans.obj <- object@kmeans.obj[[1]]
  kmeans.col <- object@kmeans.col[[1]]
  if (k.cells > 0) {
    object@data.info[names(x = kmeans.col$cluster), "kmeans.ident"] <- kmeans.col$cluster
  }
  if (set.ident && (k.cells > 0)) {
    object <- SetIdent(
      object = object,
      cells.use = names(x = kmeans.col$cluster),
      ident.use = kmeans.col$cluster
    )
  }
  if (do.plot) {
    disp.data <- minmax(
      data = kmeans.data[order(kmeans.obj$cluster[genes.use]), ],
      min = data.cut * (-1),
      max = data.cut
    )
    #DoHeatmap(object,object@cell.names,names(sort(kmeans.obj$cluster)),data.cut*(-1),data.cut,col.use = k.cols,...)
    KMeansHeatmap(object = object)
  }
  return(object)
}
