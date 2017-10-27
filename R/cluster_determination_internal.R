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
#
#' @importFrom utils read.table write.table
#
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
  edge <- edge[!duplicated(x = edge[, 1:2]), ]
  temp.file.location <- SetIfNull(x = temp.file.location, default = tempfile())
  unique_ID <- sample(x = 10000:99999, size = 1)
  edge_file <- paste0(temp.file.location, "_edge_", unique_ID, ".txt")
  output_file <- paste0(temp.file.location, "_output_", unique_ID, ".txt")
  while (file.exists(edge_file)) {
    unique_ID <- sample(x = 10000:99999, size = 1)
    edge_file <- paste0(temp.file.location, "_edge_", unique_ID, ".txt")
    output_file <- paste0(temp.file.location, "_output", unique_ID, ".txt")
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
  return(object)
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


# Set up kmeans class
# This is an infrequently used slot, but some people still find it very useful to do kmeans clustering
# and in particular, to do so at the gene level
# potential to be updated in the future

kmeans.info <- setClass(
  Class = "kmeans.info",
  slots = list(
    gene.kmeans.obj = "ANY",
    cell.kmeans.obj = "ANY"
  )
)

