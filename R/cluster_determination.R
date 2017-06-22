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
#' @param temp.file.location Directory where intermediate files will be written. Specify the
#'        ABSOLUTE path.
#' @param distance.matrix Build SNN from distance matrix (experimental)
#'
#' @importFrom FNN get.knn
#' @importFrom igraph plot.igraph graph.adjlist
#' @importFrom Matrix sparseMatrix
#'
#' @return Returns a Seurat object and optionally the SNN matrix,
#'         object@@ident has been updated with new cluster info
#'
#' @export
#'
FindClusters <- function(
  object,
  genes.use = NULL,
  reduction.type = "pca",
  pc.use = NULL,
  k.param = 30,
  k.scale = 25,
  plot.SNN = FALSE,
  prune.SNN = 1/15,
  save.SNN = FALSE,
  reuse.SNN = FALSE,
  do.sparse = FALSE,
  modularity.fxn = 1,
  resolution = 0.8,
  algorithm = 1,
  n.start = 100,
  n.iter = 10,
  random.seed = 0,
  print.output = TRUE,
  temp.file.location = NULL,
  distance.matrix = NULL
) {
  # for older objects without the snn.k slot
  if(typeof(x = validObject(object = object, test = TRUE)) == "character") {
    object@snn.k <- numeric()
  }
  snn.built <- FALSE
  if (.hasSlot(object = object, name = "snn.dense")) {
    if (length(x = object@snn.dense) > 1) {
      snn.built <- TRUE
    }
  }
  if (.hasSlot(object = object, name = "snn.sparse")) {
    if (length(x = object@snn.sparse) > 1) {
      snn.built <- TRUE
    }
  }
  if ((
    missing(x = genes.use) && missing(x = pc.use) && missing(x = k.param) &&
    missing(x = k.scale) && missing(x = prune.SNN) && snn.built
  ) || reuse.SNN) {
    save.SNN <- TRUE
    if (reuse.SNN && !snn.built) {
      stop("No SNN stored to reuse.")
    }
    if (reuse.SNN && (
      ! missing(x = genes.use) || ! missing(x = pc.use) || ! missing(x = k.param)
      || ! missing(x = k.scale) || ! missing(x = prune.SNN)
    )) {
      warning("SNN was not be rebuilt with new parameters. Continued with stored SNN. To suppress this
                      warning, remove all SNN building parameters.")
    }
  } else { # if any SNN building parameters are provided or it hasn't been built, build a new SNN
    object <- BuildSNN(
      object = object,
      genes.use = genes.use,
      reduction.type = reduction.type,
      pc.use = pc.use,
      k.param = k.param,
      k.scale = k.scale,
      plot.SNN = plot.SNN,
      prune.SNN = prune.SNN,
      do.sparse = do.sparse,
      print.output = print.output,
      distance.matrix = distance.matrix
    )
  }
  # deal with sparse SNNs
  if (length(x = object@snn.sparse) > 1) {
    SNN.use <- object@snn.sparse
  } else {
    SNN.use <- object@snn.dense
  }
  for (r in resolution) {
    object <- RunModularityClustering(
      object = object,
      SNN = SNN.use,
      modularity = modularity.fxn,
      resolution = r,
      algorithm = algorithm,
      n.start = n.start,
      n.iter = n.iter,
      random.seed = random.seed,
      print.output = print.output,
      temp.file.location = temp.file.location
    )
    object <- GroupSingletons(object = object, SNN = SNN.use)
    name <- paste0("res.", r)
    object <- StashIdent(object = object, save.name = name)
  }
  if (!save.SNN) {
    object@snn.sparse <- sparseMatrix(1, 1, x = 1)
    object@snn.dense <- matrix()
    object@snn.k <- integer()
  }
  return(object)
}

# Documentation
###############
#' @export
#'
GetClusters <- function(object) {
  return(data.frame(object@ident))
}

# Documentation
###############
#' @export
#'
SetClusters <- function(object, clusters = NULL) {
  cells.use <- rownames(x = clusters)
  ident.use <- as.numeric(x = clusters[, 1])
  object <- SetIdent(
    object = object,
    cells.use = cells.use,
    ident.use = ident.use
  )
  return(object)
}

# Documentation
###############
#' @export
#'
SaveClusters <- function(object, file) {
  my.clusters <- GetClusters(object = object)
  write.table(my.clusters, file = file, sep="\t", quote = FALSE)
}


# Documentation
###############
#' @export
#'
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

# Documentation
###############
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
  temp.file.location <- set.ifnull(x = temp.file.location, y = seurat.dir)
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

# Documentation
###############
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
        ) # Column
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
