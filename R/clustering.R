#' @include generics.R
#'
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Methods for Seurat-defined generics
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @importFrom pbapply pblapply
#' @importFrom future.apply future_lapply
#'
#' @param modularity.fxn Modularity function (1 = standard; 2 = alternative).
#' @param initial.membership,weights,node.sizes Parameters to pass to the Python leidenalg function.
#' @param resolution Value of the resolution parameter, use a value above
#' (below) 1.0 if you want to obtain a larger (smaller) number of communities.
#' @param algorithm Algorithm for modularity optimization (1 = original Louvain
#' algorithm; 2 = Louvain algorithm with multilevel refinement; 3 = SLM
#' algorithm; 4 = Leiden algorithm). Leiden requires the leidenalg python.
#' @param n.start Number of random starts.
#' @param n.iter Maximal number of iterations per random start.
#' @param random.seed Seed of the random number generator.
#' @param temp.file.location Directory where intermediate files will be written.
#' Specify the ABSOLUTE path.
#' @param edge.file.name Edge file to use as input for modularity optimizer jar.
#' @param verbose Print output
#'
#' @rdname FindClusters
#' @export
#'
FindClusters.default <- function(
  object,
  modularity.fxn = 1,
  initial.membership = NULL,
  weights = NULL,
  node.sizes = NULL,
  resolution = 0.8,
  algorithm = 1,
  n.start = 10,
  n.iter = 10,
  random.seed = 0,
  temp.file.location = NULL,
  edge.file.name = NULL,
  verbose = TRUE,
  ...
) {
  if (is.null(x = object)) {
    stop("Please provide an SNN graph")
  }
  if (tolower(x = algorithm) == "louvain") {
    algorithm <- 1
  }
  if (tolower(x = algorithm) == "leiden") {
    algorithm <- 4
  }
  if (PlanThreads() > 1) {
    clustering.results <- future_lapply(
      X = resolution,
      FUN = function(r) {
        if (algorithm %in% c(1:3)) {
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
            edge.file.name = edge.file.name
          )
        } else if (algorithm == 4) {
          ids <- RunLeiden(
            adj_mat = object,
            partition.type = "RBConfigurationVertexPartition",
            initial.membership = initial.membership,
            weights = weights,
            node.sizes = node.sizes,
            resolution.parameter = resolution
          )
        } else {
          stop("algorithm not recognised, please specify as an integer or string")
        }
        names(x = ids) <- colnames(x = object)
        ids <- GroupSingletons(ids = ids, SNN = object, verbose = verbose)
        results <- list(factor(x = ids))
        names(x = results) <- paste0('res.', r)
        return(results)
      }
    )
    clustering.results <- as.data.frame(x = clustering.results)
  } else {
    clustering.results <- data.frame(row.names = colnames(x = object))
    for (r in resolution) {
      if (algorithm %in% c(1:3)) {
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
      } else if (algorithm == 4) {
        ids <- RunLeiden(
          adj_mat = object,
          partition.type = "RBConfigurationVertexPartition",
          initial.membership = initial.membership,
          weights = weights,
          node.sizes = node.sizes,
          resolution.parameter = r
        )
      } else {
        stop("algorithm not recognised, please specify as an integer or string")
      }
      names(x = ids) <- colnames(x = object)
      ids <- GroupSingletons(ids = ids, SNN = object, verbose = verbose)
      clustering.results[, paste0("res.", r)] <- factor(x = ids)
    }
  }
  return(clustering.results)
}

#' @importFrom methods is
#'
#' @param graph.name Name of graph to use for the clustering algorithm
#'
#' @rdname FindClusters
#' @export
#' @method FindClusters Seurat
#'
FindClusters.Seurat <- function(
  object,
  graph.name = NULL,
  modularity.fxn = 1,
  initial.membership = NULL,
  weights = NULL,
  node.sizes = NULL,
  resolution = 0.8,
  algorithm = 1,
  n.start = 10,
  n.iter = 10,
  random.seed = 0,
  temp.file.location = NULL,
  edge.file.name = NULL,
  verbose = TRUE,
  ...
) {
  graph.name <- graph.name %||% paste0(DefaultAssay(object = object), "_snn")
  if (!graph.name %in% names(x = object)) {
    stop("Provided graph.name not present in Seurat object")
  }
  if (!is(object = object[[graph.name]], class2 = "Graph")) {
    stop("Provided graph.name does not correspond to a graph object.")
  }
  clustering.results <- FindClusters(
    object = object[[graph.name]],
    modularity.fxn = modularity.fxn,
    initial.membership = initial.membership,
    weights = weights,
    node.sizes = node.sizes,
    resolution = resolution,
    algorithm = algorithm,
    n.start = n.start,
    n.iter = n.iter,
    random.seed = random.seed,
    temp.file.location = temp.file.location,
    edge.file.name = edge.file.name,
    verbose = verbose
  )
  colnames(x = clustering.results) <- paste0(graph.name, "_", colnames(x = clustering.results))
  object <- AddMetaData(object = object, metadata = clustering.results)
  Idents(object = object) <- colnames(x = clustering.results)[ncol(x = clustering.results)]
  object[['seurat_clusters']] <- Idents(object = object)
  object <- LogSeuratCommand(object)
  return(object)
}

#' @param distance.matrix Boolean value of whether the provided matrix is a
#' distance matrix; note, for objects of class \code{dist}, this parameter will
#' be set automatically
#' @param k.param Defines k for the k-nearest neighbor algorithm
#' @param compute.SNN also compute the shared nearest neighbor graph
#' @param prune.SNN Sets the cutoff for acceptable Jaccard index when
#' computing the neighborhood overlap for the SNN construction. Any edges with
#' values less than or equal to this will be set to 0 and removed from the SNN
#' graph. Essentially sets the strigency of pruning (0 --- no pruning, 1 ---
#' prune everything).
#' @param nn.eps Error bound when performing nearest neighbor seach using RANN;
#' default of 0.0 implies exact nearest neighbor search
#' @param verbose Whether or not to print output to the console
#' @param force.recalc Force recalculation of SNN.
#'
#' @importFrom RANN nn2
#' @importFrom methods as
#'
#' @rdname FindNeighbors
#' @export
#' @method FindNeighbors default
#'
FindNeighbors.default <- function(
  object,
  distance.matrix = FALSE,
  k.param = 20,
  compute.SNN = TRUE,
  prune.SNN = 1/15,
  nn.eps = 0,
  verbose = TRUE,
  force.recalc = FALSE,
  ...
) {
  if (is.null(x = dim(x = object))) {
    warning(
      "Object should have two dimensions, attempting to coerce to matrix",
      call. = FALSE
    )
    object <- as.matrix(x = object)
  }
  if (is.null(rownames(x = object))) {
    stop("Please provide rownames (cell names) with the input object")
  }
  n.cells <- nrow(x = object)
  if (n.cells < k.param) {
    warning(
      "k.param set larger than number of cells. Setting k.param to number of cells - 1.",
      call. = FALSE
    )
    k.param <- n.cells - 1
  }
  # find the k-nearest neighbors for each single cell
  if (!distance.matrix) {
    if (verbose) {
      message("Computing nearest neighbor graph")
    }
    my.knn <- nn2(
      data = object,
      k = k.param,
      searchtype = 'standard',
      eps = nn.eps)
    nn.ranked <- my.knn$nn.idx
  } else {
    if (verbose) {
      message("Building SNN based on a provided distance matrix")
    }
    knn.mat <- matrix(data = 0, ncol = k.param, nrow = n.cells)
    knd.mat <- knn.mat
    for (i in 1:n.cells) {
      knn.mat[i, ] <- order(object[i, ])[1:k.param]
      knd.mat[i, ] <- object[i, knn.mat[i, ]]
    }
    nn.ranked <- knn.mat[, 1:k.param]
  }
  # convert nn.ranked into a Graph
  j <- as.numeric(x = t(x = nn.ranked))
  i <- ((1:length(x = j)) - 1) %/% k.param + 1
  nn.matrix <- as(object = sparseMatrix(i = i, j = j, x = 1), Class = "Graph")
  rownames(x = nn.matrix) <- rownames(x = object)
  colnames(x = nn.matrix) <- rownames(x = object)
  neighbor.graphs <- list(nn = nn.matrix)
  if (compute.SNN) {
    if (verbose) {
      message("Computing SNN")
    }
    snn.matrix <- ComputeSNN(
      nn_ranked = nn.ranked,
      prune = prune.SNN
    )
    rownames(x = snn.matrix) <- rownames(x = object)
    colnames(x = snn.matrix) <- rownames(x = object)
    snn.matrix <- as.Graph(x = snn.matrix)
    neighbor.graphs[["snn"]] <- snn.matrix
  }
  return(neighbor.graphs)
}

#' @rdname FindNeighbors
#' @export
#' @method FindNeighbors Assay
#'
FindNeighbors.Assay <- function(
  object,
  features = NULL,
  k.param = 20,
  compute.SNN = TRUE,
  prune.SNN = 1/15,
  nn.eps = 0,
  verbose = TRUE,
  force.recalc = FALSE,
  ...
) {
  features <- features %||% VariableFeatures(object = object)
  data.use <- t(x = GetAssayData(object = object, slot = "data")[features, ])
  neighbor.graphs <- FindNeighbors(
    object = data.use,
    k.param = k.param,
    compute.SNN = compute.SNN,
    prune.SNN = prune.SNN,
    nn.eps = nn.eps,
    verbose = verbose,
    force.recalc = force.recalc
  )
  return(neighbor.graphs)
}

#' @rdname FindNeighbors
#' @export
#' @method FindNeighbors dist
#'
FindNeighbors.dist <- function(
  object,
  k.param = 20,
  compute.SNN = TRUE,
  prune.SNN = 1/15,
  nn.eps = 0,
  verbose = TRUE,
  force.recalc = FALSE,
  ...
) {
  return(FindNeighbors(
    object = as.matrix(x = object),
    distance.matrix = TRUE,
    k.param = k.param,
    compute.SNN = compute.SNN,
    prune.SNN = prune.SNN,
    nn.eps = nn.eps,
    verbose = verbose,
    force.recalc = force.recalc,
    ...
  ))
}

#' @param assay Assay to use in construction of SNN
#' @param features Features to use as input for building the SNN
#' @param reduction Reduction to use as input for building the SNN
#' @param dims Dimensions of reduction to use as input
#' @param do.plot Plot SNN graph on tSNE coordinates
#' @param graph.name Optional naming parameter for stored SNN graph. Default is
#' assay.name_snn.
#'
#' @importFrom igraph graph.adjacency plot.igraph E
#'
#' @rdname FindNeighbors
#' @export
#' @method FindNeighbors Seurat
#'
FindNeighbors.Seurat <- function(
  object,
  reduction = "pca",
  dims = 1:10,
  assay = NULL,
  features = NULL,
  k.param = 20,
  compute.SNN = TRUE,
  prune.SNN = 1/15,
  nn.eps = 0,
  verbose = TRUE,
  force.recalc = FALSE,
  do.plot = FALSE,
  graph.name = NULL,
  ...
) {
  if (!is.null(x = dims)) {
    assay <- assay %||% DefaultAssay(object = object)
    data.use <- Embeddings(object = object[[reduction]])
    if (max(dims) > ncol(x = data.use)) {
      stop("More dimensions specified in dims than have been computed")
    }
    data.use <- data.use[, dims]
    neighbor.graphs <- FindNeighbors(
      object = data.use,
      k.param = k.param,
      compute.SNN = compute.SNN,
      prune.SNN = prune.SNN,
      nn.eps = nn.eps,
      verbose = verbose,
      force.recalc = force.recalc
    )
  } else {
    assay <- assay %||% DefaultAssay(object = object)
    data.use <- GetAssay(object = object, assay = assay)
    neighbor.graphs <- FindNeighbors(
      object = data.use,
      features = features,
      k.param = k.param,
      compute.SNN = compute.SNN,
      prune.SNN = prune.SNN,
      nn.eps = nn.eps,
      verbose = verbose,
      force.recalc = force.recalc
    )
  }
  graph.name <- graph.name %||% paste0(assay, "_", names(x = neighbor.graphs))
  for (ii in 1:length(x = graph.name)) {
    object[[graph.name[[ii]]]] <- neighbor.graphs[[ii]]
  }
  if (do.plot) {
    if (!"tsne" %in% names(x = object@reductions)) {
      warning("Please compute a tSNE for SNN visualization. See RunTSNE().")
    } else {
      if (nrow(x = Embeddings(object = object[["tsne"]])) != ncol(x = object)) {
        warning("Please compute a tSNE for SNN visualization. See RunTSNE().")
      } else {
        net <- graph.adjacency(
          adjmatrix = as.matrix(x = neighbor.graphs[[2]]),
          mode = "undirected",
          weighted = TRUE,
          diag = FALSE
        )
        plot.igraph(
          x = net,
          layout = as.matrix(x = Embeddings(object = object[["tsne"]])),
          edge.width = E(graph = net)$weight,
          vertex.label = NA,
          vertex.size = 0
        )
      }
    }
  }
  object <- LogSeuratCommand(object = object)
  return(object)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Internal
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Group single cells that make up their own cluster in with the cluster they are
# most connected to.
#
# @param ids Named vector of cluster ids
# @param SNN SNN graph used in clustering
#
# @return Returns Seurat object with all singletons merged with most connected cluster
#
GroupSingletons <- function(ids, SNN, verbose) {
  # identify singletons
  singletons <- c()
  singletons <- names(which(table(ids) == 1))
  singletons <- intersect(unique(ids), singletons)
  # calculate connectivity of singletons to other clusters, add singleton
  # to cluster it is most connected to
  cluster_names <- as.character(unique(x = ids))
  cluster_names <- setdiff(x = cluster_names, y = singletons)
  connectivity <- vector(mode = "numeric", length = length(x = cluster_names))
  names(x = connectivity) <- cluster_names
  new.ids <- ids
  for (i in singletons) {
    i.cells <- names(which(ids == i))
    for (j in cluster_names) {
      j.cells <- names(which(ids == j))
      subSNN <- SNN[i.cells, j.cells]
      set.seed(1) # to match previous behavior, random seed being set in WhichCells
      if (is.object(x = subSNN)) {
        connectivity[j] <- sum(subSNN) / (nrow(x = subSNN) * ncol(x = subSNN))
      } else {
        connectivity[j] <- mean(x = subSNN)
      }
    }
    m <- max(connectivity, na.rm = T)
    mi <- which(x = connectivity == m, arr.ind = TRUE)
    closest_cluster <- sample(x = names(x = connectivity[mi]), 1)
    ids[i.cells] <- closest_cluster
  }
  if (length(x = singletons) > 0 && verbose) {
    message(paste(
      length(x = singletons),
      "singletons identified.",
      length(x = unique(x = ids)),
      "final clusters."
    ))
  }
  return(ids)
}

# Run Leiden clustering algorithm
#
# Implements the Leiden clustering algorithm in R using reticulate
# to run the Python version. Requires the python "leidenalg" and "igraph" modules
# to be installed. Returns a vector of partition indices.
#
# @param adj_mat An adjacency matrix or SNN matrix
# @param partition.type Type of partition to use for Leiden algorithm.
# Defaults to RBConfigurationVertexPartition. Options include: ModularityVertexPartition,
# RBERVertexPartition, CPMVertexPartition, MutableVertexPartition,
# SignificanceVertexPartition, SurpriseVertexPartition (see the Leiden python
# module documentation for more details)
# @param initial.membership,weights,node.sizes Parameters to pass to the Python leidenalg function.
# @param resolution.parameter A parameter controlling the coarseness of the clusters
# for Leiden algorithm. Higher values lead to more clusters. (defaults to 1.0 for
# partition types that accept a resolution parameter)
#
# @keywords graph network igraph mvtnorm simulation
#
#' @importFrom reticulate py_module_available import r_to_py
#
# @author Tom Kelly
#
# @export
#
RunLeiden <- function(
  adj_mat,
  partition.type = c(
    'RBConfigurationVertexPartition',
    'ModularityVertexPartition',
    'RBERVertexPartition',
    'CPMVertexPartition',
    'MutableVertexPartition',
    'SignificanceVertexPartition',
    'SurpriseVertexPartition'
  ),
  initial.membership = NULL,
  weights = NULL,
  node.sizes = NULL,
  resolution.parameter = 1
) {
  if (!py_module_available(module = 'leidenalg')) {
    stop("Cannot find Leiden algorithm, please install through pip (e.g. pip install leidenalg).")
  }
  # import python modules with reticulate
  leidenalg <- import(module = "leidenalg")
  ig <- import(module = "igraph")
  # convert matrix input
  adj_mat <- as.matrix(x = ceiling(x = adj_mat))
  # convert to python numpy.ndarray, then a list
  adj_mat_py <- r_to_py(x = adj_mat)
  adj_mat_py <- adj_mat_py$tolist()
  # convert graph structure to a Python compatible object
  snn_graph <- ig$Graph$Adjacency(adj_mat_py)
  # compute partitions
  part <- switch(
    EXPR = partition.type,
    'RBConfigurationVertexPartition' = leidenalg$find_partition(
      snn_graph,
      leidenalg$RBConfigurationVertexPartition,
      initial_membership = initial.membership,
      weights = weights,
      resolution_parameter = resolution.parameter
    ),
    'ModularityVertexPartition' = leidenalg$find_partition(
      snn_graph,
      leidenalg$ModularityVertexPartition,
      initial_membership = initial.membership,
      weights = weights
    ),
    'RBERVertexPartition' = leidenalg$find_partition(
      snn_graph,
      leidenalg$RBERVertexPartition,
      initial_membership = initial.membership,
      weights = weights,
      node_sizes = node.sizes,
      resolution_parameter = resolution.parameter
    ),
    'CPMVertexPartition' = leidenalg$find_partition(
      snn_graph,
      leidenalg$CPMVertexPartition,
      initial_membership = initial.membership,
      weights = weights,
      node_sizes = node.sizes,
      resolution_parameter = resolution.parameter
    ),
    'MutableVertexPartition' = leidenalg$find_partition(
      snn_graph,
      leidenalg$MutableVertexPartition,
      initial_membership = initial.membership
    ),
    'SignificanceVertexPartition' = leidenalg$find_partition(
      snn_graph,
      leidenalg$SignificanceVertexPartition,
      initial_membership = initial.membership,
      node_sizes = node.sizes,
      resolution_parameter = resolution.parameter
    ),
    'SurpriseVertexPartition' = leidenalg$find_partition(
      snn_graph,
      leidenalg$SurpriseVertexPartition,
      initial_membership = initial.membership,
      weights = weights,
      node_sizes = node.sizes
    ),
    stop("please specify a partition type as a string out of those documented")
  )
  return(part$membership + 1)
}

# Runs the modularity optimizer (C++ port of java program ModularityOptimizer.jar)
#
# @param SNN SNN matrix to use as input for the clustering algorithms
# @param modularity Modularity function to use in clustering (1 = standard; 2 = alternative)
# @param resolution Value of the resolution parameter, use a value above (below) 1.0 if you want to obtain a larger (smaller) number of communities
# @param algorithm Algorithm for modularity optimization (1 = original Louvain algorithm; 2 = Louvain algorithm with multilevel refinement; 3 = SLM algorithm; 4 = Leiden algorithm). Leiden requires the leidenalg python module.
# @param n.start Number of random starts
# @param n.iter Maximal number of iterations per random start
# @param random.seed Seed of the random number generator
# @param print.output Whether or not to print output to the console
# @param temp.file.location Deprecated and no longer used
# @param edge.file.name Path to edge file to use
#
# @return Seurat object with identities set to the results of the clustering procedure
#
#' @importFrom utils read.table write.table
#
RunModularityClustering <- function(
  SNN = matrix(),
  modularity = 1,
  resolution = 0.8,
  algorithm = 1,
  n.start = 10,
  n.iter = 10,
  random.seed = 0,
  print.output = TRUE,
  temp.file.location = NULL,
  edge.file.name = NULL
) {
  edge_file <- edge.file.name %||% ''
  clusters <- RunModularityClusteringCpp(
    SNN,
    modularity,
    resolution,
    algorithm,
    n.start,
    n.iter,
    random.seed,
    print.output,
    edge_file
  )
  return(clusters)
}
