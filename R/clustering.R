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
#' @importFrom future nbrOfWorkers
#'
#' @param modularity.fxn Modularity function (1 = standard; 2 = alternative).
#' @param initial.membership,node.sizes Parameters to pass to the Python leidenalg function.
#' @param resolution Value of the resolution parameter, use a value above
#' (below) 1.0 if you want to obtain a larger (smaller) number of communities.
#' @param algorithm Algorithm for modularity optimization (1 = original Louvain
#' algorithm; 2 = Louvain algorithm with multilevel refinement; 3 = SLM
#' algorithm; 4 = Leiden algorithm). Leiden requires the leidenalg python.
#' @param method Method for running leiden (defaults to matrix which is fast for small datasets).
#' Enable method = "igraph" to avoid casting large data to a dense matrix.
#' @param n.start Number of random starts.
#' @param n.iter Maximal number of iterations per random start.
#' @param random.seed Seed of the random number generator.
#' @param group.singletons Group singletons into nearest cluster. If FALSE, assign all singletons to
#' a "singleton" group
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
  node.sizes = NULL,
  resolution = 0.8,
  method = "matrix",
  algorithm = 1,
  n.start = 10,
  n.iter = 10,
  random.seed = 0,
  group.singletons = TRUE,
  temp.file.location = NULL,
  edge.file.name = NULL,
  verbose = TRUE,
  ...
) {
  CheckDots(...)
  if (is.null(x = object)) {
    stop("Please provide an SNN graph")
  }
  if (tolower(x = algorithm) == "louvain") {
    algorithm <- 1
  }
  if (tolower(x = algorithm) == "leiden") {
    algorithm <- 4
  }
  if (nbrOfWorkers() > 1) {
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
            object = object,
            method = method,
            partition.type = "RBConfigurationVertexPartition",
            initial.membership = initial.membership,
            node.sizes = node.sizes,
            resolution.parameter = r,
            random.seed = random.seed,
            n.iter = n.iter
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
          object = object,
          method = method,
          partition.type = "RBConfigurationVertexPartition",
          initial.membership = initial.membership,
          node.sizes = node.sizes,
          resolution.parameter = r,
          random.seed = random.seed,
          n.iter = n.iter
        )
      } else {
        stop("algorithm not recognised, please specify as an integer or string")
      }
      names(x = ids) <- colnames(x = object)
      ids <- GroupSingletons(ids = ids, SNN = object, group.singletons = group.singletons, verbose = verbose)
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
  node.sizes = NULL,
  resolution = 0.8,
  method = "matrix",
  algorithm = 1,
  n.start = 10,
  n.iter = 10,
  random.seed = 0,
  group.singletons = TRUE,
  temp.file.location = NULL,
  edge.file.name = NULL,
  verbose = TRUE,
  ...
) {
  CheckDots(...)
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
    node.sizes = node.sizes,
    resolution = resolution,
    method = method,
    algorithm = algorithm,
    n.start = n.start,
    n.iter = n.iter,
    random.seed = random.seed,
    group.singletons = group.singletons,
    temp.file.location = temp.file.location,
    edge.file.name = edge.file.name,
    verbose = verbose,
    ...
  )
  colnames(x = clustering.results) <- paste0(graph.name, "_", colnames(x = clustering.results))
  object <- AddMetaData(object = object, metadata = clustering.results)
  Idents(object = object) <- colnames(x = clustering.results)[ncol(x = clustering.results)]
  levels <- levels(x = object)
  levels <- tryCatch(
    expr = as.numeric(x = levels),
    warning = function(...) {
      return(levels)
    },
    error = function(...) {
      return(levels)
    }
  )
  Idents(object = object) <- factor(x = Idents(object = object), levels = sort(x = levels))
  object[['seurat_clusters']] <- Idents(object = object)
  cmd <- LogSeuratCommand(object = object, return.command = TRUE)
  slot(object = cmd, name = 'assay.used') <- DefaultAssay(object = object[[graph.name]])
  object[[slot(object = cmd, name = 'name')]] <- cmd
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
#' @param nn.method Method for nearest neighbor finding. Options include: rann,
#' annoy
#' @param annoy.metric Distance metric for annoy. Options include: euclidean,
#' cosine, manhattan, and hamming
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
  nn.method = 'rann',
  annoy.metric = "euclidean",
  nn.eps = 0,
  verbose = TRUE,
  force.recalc = FALSE,
  ...
) {
  CheckDots(...)
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
    nn.ranked <- NNHelper(
      data = object,
      k = k.param,
      method = nn.method,
      searchtype = "standard",
      eps = nn.eps,
      metric = annoy.metric)
    nn.ranked <- nn.ranked$nn.idx
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
  nn.matrix <- as(object = sparseMatrix(i = i, j = j, x = 1, dims = c(nrow(x = object), nrow(x = object))), Class = "Graph")
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
  nn.method = 'rann',
  annoy.metric = "euclidean",
  nn.eps = 0,
  verbose = TRUE,
  force.recalc = FALSE,
  ...
) {
  CheckDots(...)
  features <- features %||% VariableFeatures(object = object)
  data.use <- t(x = GetAssayData(object = object, slot = "data")[features, ])
  neighbor.graphs <- FindNeighbors(
    object = data.use,
    k.param = k.param,
    compute.SNN = compute.SNN,
    prune.SNN = prune.SNN,
    nn.method = nn.method,
    annoy.metric = annoy.metric,
    nn.eps = nn.eps,
    verbose = verbose,
    force.recalc = force.recalc,
    ...
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
  nn.method = "rann",
  annoy.metric = "euclidean",
  nn.eps = 0,
  verbose = TRUE,
  force.recalc = FALSE,
  ...
) {
  CheckDots(...)
  return(FindNeighbors(
    object = as.matrix(x = object),
    distance.matrix = TRUE,
    k.param = k.param,
    compute.SNN = compute.SNN,
    prune.SNN = prune.SNN,
    nn.eps = nn.eps,
    nn.method = nn.method,
    annoy.metric = annoy.metric,
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
  nn.method = "rann",
  annoy.metric = "euclidean",
  nn.eps = 0,
  verbose = TRUE,
  force.recalc = FALSE,
  do.plot = FALSE,
  graph.name = NULL,
  ...
) {
  CheckDots(...)
  if (!is.null(x = dims)) {
    # assay <- assay %||% DefaultAssay(object = object)
    assay <- DefaultAssay(object = object[[reduction]])
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
      nn.method = nn.method,
      annoy.metric = annoy.metric,
      nn.eps = nn.eps,
      verbose = verbose,
      force.recalc = force.recalc,
      ...
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
      nn.method = nn.method,
      annoy.metric = annoy.metric,
      nn.eps = nn.eps,
      verbose = verbose,
      force.recalc = force.recalc,
      ...
    )
  }
  graph.name <- graph.name %||% paste0(assay, "_", names(x = neighbor.graphs))
  for (ii in 1:length(x = graph.name)) {
    DefaultAssay(object = neighbor.graphs[[ii]]) <- assay
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

# Run annoy
#
# @param data Data to build the index with
# @param query A set of data to be queried against data
# @param metric Distance metric; can be one of "euclidean", "cosine", "manhattan",
# "hamming"
# @param n.trees More trees gives higher precision when querying
# @param k Number of neighbors
# @param search.k During the query it will inspect up to search_k nodes which
# gives you a run-time tradeoff between better accuracy and speed.
# @ param include.distance Include the corresponding distances
#
AnnoyNN <- function(data, query = data, metric = "euclidean", n.trees = 50, k,
                    search.k = -1, include.distance = TRUE) {
  idx <- AnnoyBuildIndex(
    data = data,
    metric = metric,
    n.trees = n.trees)
  nn <- AnnoySearch(
    index = idx,
    query = query,
    k = k,
    search.k = search.k,
    include.distance = include.distance)
  return(nn)
}

# Build the annoy index
#
# @param data Data to build the index with
# @param metric Distance metric; can be one of "euclidean", "cosine", "manhattan",
# "hamming"
# @param n.trees More trees gives higher precision when querying
#' @importFrom RcppAnnoy AnnoyEuclidean AnnoyAngular AnnoyManhattan AnnoyHamming
#
AnnoyBuildIndex <- function(data, metric = "euclidean", n.trees = 50) {
  f <- ncol(x = data)
  a <- switch(
    EXPR = metric,
    "euclidean" =  new(Class = RcppAnnoy::AnnoyEuclidean, f),
    "cosine" = new(Class = RcppAnnoy::AnnoyAngular, f),
    "manhattan" = new(Class = RcppAnnoy::AnnoyManhattan, f),
    "hamming" = new(Class = RcppAnnoy::AnnoyHamming, f),
    stop ("Invalid metric")
  )
  for (ii in seq(nrow(x = data))) {
    a$addItem(ii - 1, data[ii, ])
  }
  a$build(n.trees)
  return(a)
}

# Search the annoy index
#
# @param Annoy index, build with AnnoyBuildIndex
# @param query A set of data to be queried against the index
# @param k Number of neighbors
# @param search.k During the query it will inspect up to search_k nodes which
# gives you a run-time tradeoff between better accuracy and speed.
# @ param include.distance Include the corresponding distances
#
AnnoySearch <- function(index, query, k, search.k = -1, include.distance = TRUE) {
  n <- nrow(x = query)
  idx <- matrix(nrow = n,  ncol = k)
  dist <- matrix(nrow = n, ncol = k)
  convert <- methods::is(index, "Rcpp_AnnoyAngular")
  res <- future_lapply(X = 1:n, FUN = function(x) {
    res <- index$getNNsByVectorList(query[x, ], k, search.k, include.distance)
    # Convert from Angular to Cosine distance
    if (convert) {
      res$dist <- 0.5 * (res$dist * res$dist)
    }
    list(res$item + 1, res$distance)
  })
  for (i in 1:n) {
    idx[i, ] <- res[[i]][[1]]
    if (include.distance) {
      dist[i, ] <- res[[i]][[2]]
    }
  }
  return(list(nn.idx = idx, nn.dists = dist))
}

# Group single cells that make up their own cluster in with the cluster they are
# most connected to.
#
# @param ids Named vector of cluster ids
# @param SNN SNN graph used in clustering
# @param group.singletons Group singletons into nearest cluster. If FALSE, assign all singletons to
# a "singleton" group
#
# @return Returns Seurat object with all singletons merged with most connected cluster
#
GroupSingletons <- function(ids, SNN, group.singletons = TRUE, verbose = TRUE) {
  # identify singletons
  singletons <- c()
  singletons <- names(x = which(x = table(ids) == 1))
  singletons <- intersect(x = unique(x = ids), singletons)
  if (!group.singletons) {
    ids[which(ids %in% singletons)] <- "singleton"
    return(ids)
  }
  # calculate connectivity of singletons to other clusters, add singleton
  # to cluster it is most connected to
  cluster_names <- as.character(x = unique(x = ids))
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

# Internal helper function to dispatch to various neighbor finding methods
#
# @param data Input data
# @param query Data to query against data
# @param k Number of nearest neighbors to compute
# @param method Nearest neighbor method to use: "rann", "annoy"
# @param ... additional parameters to specific neighbor finding method
#
NNHelper <- function(data, query = data, k, method, ...) {
  args <- as.list(x = sys.frame(which = sys.nframe()))
  args <- c(args, list(...))
  return(
    switch(
      EXPR = method,
      "rann" = {
        args <- args[intersect(x = names(x = args), y = names(x = formals(fun = nn2)))]
        do.call(what = 'nn2', args = args)
      },
      "annoy" = {
        args <- args[intersect(x = names(x = args), y = names(x = formals(fun = AnnoyNN)))]
        do.call(what = 'AnnoyNN', args = args)
      },
      stop("Invalid method. Please choose one of 'rann', 'annoy'")
    )
  )
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
# @param initial.membership,node.sizes Parameters to pass to the Python leidenalg function.
# @param resolution.parameter A parameter controlling the coarseness of the clusters
# for Leiden algorithm. Higher values lead to more clusters. (defaults to 1.0 for
# partition types that accept a resolution parameter)
# @param random.seed Seed of the random number generator
# @param n.iter Maximal number of iterations per random start
#
# @keywords graph network igraph mvtnorm simulation
#
#' @importFrom leiden leiden
#' @importFrom reticulate py_module_available
#' @importFrom igraph graph_from_adjacency_matrix graph_from_adj_list
#
# @author Tom Kelly
#
# @export
#
RunLeiden <- function(
  object,
  method = c("matrix", "igraph"),
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
  node.sizes = NULL,
  resolution.parameter = 1,
  random.seed = 0,
  n.iter = 10
) {
  if (!py_module_available(module = 'leidenalg')) {
    stop(
      "Cannot find Leiden algorithm, please install through pip (e.g. pip install leidenalg).",
      call. = FALSE
    )
  }
  switch(
    EXPR = method,
    "matrix" = {
      input <- as(object = object, Class = "matrix")
      },
    "igraph" = {
      input <- if (inherits(x = object, what = 'list')) {
        graph_from_adj_list(adjlist = object)
      } else if (inherits(x = object, what = c('dgCMatrix', 'matrix', 'Matrix'))) {
        if (inherits(x = object, what = 'Graph')) {
          object <- as(object = object, Class = "dgCMatrix")
        }
          graph_from_adjacency_matrix(adjmatrix = object, weighted = TRUE)
      } else if (inherits(x = object, what = 'igraph')) {
        object
      } else {
        stop(
          "Method for Leiden not found for class", class(x = object), 
           call. = FALSE
        )
      }
    },
    stop("Method for Leiden must be either 'matrix' or igraph'")
  )
  #run leiden from CRAN package (calls python with reticulate)
  partition <- leiden(
    object = input,
    partition_type = partition.type,
    initial_membership = initial.membership,
    weights = NULL,
    node_sizes = node.sizes,
    resolution_parameter = resolution.parameter,
    seed = random.seed,
    n_iterations = n.iter
  )
  return(partition)
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

#' Predict expression value from knn
#'
#' @param object The object used to calculate knn
#' @param nn.idx k near neighbour indices. A cells x k matrix.
#' @param assay Assay used for prediction
#' @param reduction Cell embedding of the reduction used for prediction
#' @param dims Number of dimensions of cell embedding
#' @param slot slot used for prediction
#' @param features features used for prediction
#' @param mean.type the type of mean,	arithmetic mean (amean) or ExpMean
#'
#' @importFrom pbapply pbapply
#' @importFrom future.apply future_apply
#' @importFrom future nbrOfWorkers
#'
#' @return return an assay containing predicted expression value in the data slot
#' @export
#' 
PredictAssay <- function(
  object,
  nn.idx,
  assay,
  reduction = NULL,
  dims = NULL, 
  return.assay = TRUE,
  cells = NULL,
  slot = "scale.data",
  features = NULL,
  mean.type = NULL,
  verbose = TRUE
){
  if (is.null(x = mean.type)) {
    mean.type <- "amean"
  }
  if (is.null(x = reduction)) {
    reference.data <- GetAssayData(object = object,
                                   assay = assay,
                                   slot = slot)
    if (is.null(x = features)) {
      features <- VariableFeatures(object[[assay]])
      if (length(x = features) == 0) {
        message("VariableFeatures are empty in the ",assay," assay, features in the ", slot, " slot will be used" )
        features <- rownames(reference.data)
        if (length(x = features) == 0) {
          stop("No features in the ",slot, " slot of the assay ",assay )
        }
      }
    }
    reference.data <- reference.data[features, ,drop = FALSE]
  } else {
    if (is.null(x = dims)) {
      stop("dims is empty")
    }
    reference.data <- t(Embeddings(object = object, reduction = reduction)[, dims])
  }
  if (nrow(x = nn.idx) > 1) {
    if (all(nn.idx[1:10,1] == 1:10)) {
      if(verbose){
        message("The nearest neighbor is the query cell itself, and it will not be used for prediction")
      }
      nn.idx <- nn.idx[,-1]
    }
    if (mean.type == "ExpMean") {
      predicted <- apply(X = nn.idx,
                         MARGIN = 1,
                         FUN = function(x) FastExpMean(mat = reference.data[,x], 
                                                       display_progress = FALSE))
    } else {
      predicted <- apply(X = nn.idx,
                         MARGIN = 1,
                         FUN = function(x) rowMeans(x = reference.data[,x], 
                                              na.rm = TRUE))
    } 
  } else if (nrow(x = nn.idx) ==1) {
    if (verbose) {
      message("Only one cell's neighbors are given, and the first one will not be used for prediction")
    }
    nn.idx <- nn.idx[,-1]
    predicted <-  rowMeans(x = reference.data[,nn.idx],  na.rm = T)
    predicted <- as.matrix(x = predicted)
  }
  if (is.null(x = cells)) {
    cells <-  Cells(object)
  }
  colnames(x = predicted) <- cells
  if (return.assay) {
    predicted.assay <- CreateAssayObject(data = predicted)
    return (predicted.assay)
  } else {
    return (predicted)
  }
}

# Calculate NN distance for the given nn.idx
#' @param nn.idx
#' @param redunction.embedding
#' @param metric
#' @param query.reduction.embedding
#' @param nearest.dist
#' 
NNdist <- function( nn.idx, 
                    redunction.embedding, 
                    metric = "euclidean",
                    query.reduction.embedding = NULL, 
                    nearest.dist = NULL){
  if (!is.list(x = nn.idx)) {
    nn.idx <- lapply(X = 1:nrow(x = nn.idx), FUN = function(x) nn.idx[x,])
  }
  if (is.null(x = query.reduction.embedding)) {
    query.reduction.embedding <- redunction.embedding
  }
  nn.dist <- fast_dist(x = query.reduction.embedding,
                       y = redunction.embedding,
                       n = nn.idx)
  if (!is.null(x = nearest.dist)) {
    nn.dist <- lapply(X = 1:nrow(x = query.reduction.embedding), 
                      FUN = function(x) {
                        r_dist = nn.dist[[x]] - nearest.dist[x]  
                        r_dist[r_dist < 0] <- 0
                        return (r_dist)
                      })
  }
  return (nn.dist)
}


# Find multi-model neighbors 
#
#' @param object The object used to calculate knn
#' @param k.nn .Number of nearest multi-model neighbors to compute
#' @param reduction.list A list of reduction name 
#' @param dims.list A list of dimentions used for the reduction
#' @param knn.range The number of approximate neighbors to compute
#' @param modality.weight the cell-specific modeality weights
#' @param verbose Whether or not to print output to the console
#'
#' @return return a list containing nn index and nn multi-model distance
#' @export
#' 
MultiModalNN <- function(object, 
                         query = NULL,
                         modality.weight = NULL,
                         k.nn =  modality.weight$params$k.nn, 
                         reduction.list =  modality.weight$params$reduction.list,
                         dims.list = modality.weight$params$dims.list,
                         knn.range = 200,
                         kernel.power = 1, 
                         nearest.dist = modality.weight$params$nearest.dist,
                         sigma.list = modality.weight$params$sigma.list,
                         l2.norm =  modality.weight$params$l2.norm, 
                         verbose = TRUE
){
  modality.weight.value <- list(modality.weight$first.modality.weight, 1 - modality.weight$first.modality.weight)
  names(x = modality.weight.value) <- unlist(x = reduction.list)
  if (class(x = object)[1] == "Seurat") {
    redunction_embedding <- lapply( X = 1:length(x = reduction.list), 
                                    FUN = function(x) {
                                      Embeddings(object = object, 
                                                 reduction = reduction.list[[x]] )[ ,dims.list[[x]] ]
                                    })
  } else {
    redunction_embedding <- object
  }
  if (is.null(x = query)) {
    query.redunction_embedding <- redunction_embedding
    query <- object
  } else {
    if (class(x = query)[1] == "Seurat") {
      query.redunction_embedding <- lapply(X = 1:length(x = reduction.list), 
                                            FUN = function(x) {
                                              Embeddings(object = query,
                                                         reduction = reduction.list[[x]] )[ ,dims.list[[x]] ]
                                            })
    } else {
      query.redunction_embedding <- query
    }
  }
  if (l2.norm) {
    query.redunction_embedding <- lapply( X = query.redunction_embedding,
                                          FUN =  function(x)  L2Norm(mat = x))
    redunction_embedding <- lapply(X = redunction_embedding, FUN = function(x) L2Norm(mat = x))
  }
  query.cell.num <- nrow(x = query.redunction_embedding[[1]])
  reduction.num <- length(x = query.redunction_embedding)
  if (verbose) {
    message("Finding multi-modal nearest neighbors")
    pb <- txtProgressBar(min = 0, max = reduction.num, style = 3)
  }
  redunction_nn <- lapply(X = 1:reduction.num, 
                           FUN = function(x) {
                             nn_x <- NNHelper(data = redunction_embedding[[x]], 
                                              query = query.redunction_embedding[[x]],
                                              k = knn.range,
                                              method = 'annoy',
                                              metric = "euclidean") 
                             if (verbose) {
                               setTxtProgressBar(pb = pb, value = x)
                             }
                             return (nn_x)
                           })
  if (verbose) {
    close(pb)
  }
  # union of rna and adt nn, remove itself from neighobors
  redunction_nn <- lapply(X = redunction_nn , 
                          FUN = function(x)  x$nn.idx[, -1]  )
  nn_idx <- lapply( X = 1:query.cell.num , 
                    FUN = function(x)  Reduce(f = union, 
                                              x = lapply(X = redunction_nn, 
                                                         FUN = function(y) y[x,] )))
  # calculate euclidean distance of all neighbors
  nn_dist <- lapply(X = 1:reduction.num,  
                    FUN = function(r) {
                      NNdist(nn.idx = nn_idx,
                             redunction.embedding = redunction_embedding[[r]], 
                             query.reduction.embedding = query.redunction_embedding[[r]], 
                             nearest.dist = nearest.dist[[r]])
                    })
  # modality weighted distance
  if (length(x = sigma.list[[1]]) == 1) {
    sigma.list <- lapply(X = sigma.list, FUN = function(x) rep(x = x, ncol(x = object) ))
  }
  nn_weighted_dist <- lapply(X = 1:reduction.num,  
                             FUN = function(r){
                               lapply(X = 1:query.cell.num,
                                        FUN = function(x) { 
                                          exp(-1*(nn_dist[[r]][[x]] / sigma.list[[r]][x] )**
                                                kernel.power) * 
                                            modality.weight.value[[r]][x] })
                             })
  nn_weighted_dist <- sapply(X = 1:query.cell.num, 
                             FUN =  function(x) { 
                               Reduce(f = "+", 
                                      x = lapply( X = 1:reduction.num, 
                                              FUN = function(r) nn_weighted_dist[[r]][[x]] )) 
                             })
  # select k nearest joint neighbors
  select_order <- lapply( X = nn_weighted_dist,
                         FUN = function(dist) {
                           order(dist, decreasing = TRUE)
                         })
  select_nn <- t(sapply( X = 1:query.cell.num,
                         FUN = function(x) nn_idx[[x]][select_order[[x]]][ 1:k.nn ])
  )
  select_dist <- t(sapply( X = 1:query.cell.num, 
                           FUN = function(x) nn_weighted_dist[[x]][select_order[[x]]][ 1:k.nn ])
  )
  select_dist <- 1 - select_dist
  rownames(x = select_nn) <- rownames(x = select_dist) <- Cells(query)
  joint.nn <- list(select_nn, select_dist)
  names(x = joint.nn) <- c("nn.idx", "nn.dists")
  return (joint.nn)
}


#' Multimodel KNN and SNN Graph Construction
#' 
#' @param object The object used to calculate knn
#' @param modality.weight the cell-specific modeality weights 
#' @param prune.SNN .Number of nearest multi-model neighbors to compute
#' @param knn.graph.name The name of multi-model knn graph
#' @param snn.graph.name The name of multi-model snn graph
#' @param joint.nn.name The name of multi-model neighbors
#' @param modality.weight.name The name of first modality weights
#' @param knn.range The number of approximate neighbors to compute
#' @param verbose Whether or not to print output to the console
#' 
#' @return return an object containing multi-model KNN, SNN and neighbors
#' @export

FindMultiModelNeighbors  <- function(object, 
                                     modality.weight = NULL,
                                     k.nn = modality.weight$params$k.nn,
                                     prune.SNN = 1/15, 
                                     knn.graph.name = "jknn",
                                     snn.graph.name = "jsnn",
                                     joint.nn.name = "joint.nn",
                                     modality.weight.name = "first.modality.weight",
                                     knn.range = 200,
                                     weighted.graph = FALSE,
                                     verbose = TRUE
){
  joint.nn <- MultiModalNN(object = object, 
                           modality.weight = modality.weight,
                           knn.range = knn.range, 
                           verbose = verbose )
  select_nn <- joint.nn$nn.idx
  select_nn_dist <- joint.nn$nn.dists 
  if (weighted.graph) {
    if(verbose){
      message("Constructing joint weighted knn graph")
    }
    joint.nn$nn.dists <- t(apply(X = joint.nn$nn.dists,
                                 MARGIN = 1, 
                                 FUN = function(x) log2(k.nn)*x/sum(x) ))
    nn.matrix <- sparseMatrix(i = 1:ncol(x = object),
                              j = 1:ncol(x = object),
                              x = 1)
    for (i in 1:ncol(x = object)) {
      nn.matrix[i, select_nn[i,]] <- joint.nn$nn.dists[i, ]
    }
  } else {
    if (verbose) {
      message("Constructing joint KNN graph")
    }
    j <- as.numeric(x = t(x = select_nn ))
    i <- ((1:length(x = j)) - 1) %/% k.nn + 1
    nn.matrix <- sparseMatrix(i = i,
                              j = j,
                              x = 1, 
                              dims = c(ncol(x = object), ncol(x = object)))
    diag(x = nn.matrix) <- 1
  }
  rownames(x = nn.matrix) <-  colnames(x = nn.matrix) <- colnames(x = object)
  nn.matrix <- nn.matrix + t(nn.matrix) - t(nn.matrix)*nn.matrix
  nn.matrix <- as.Graph(x = nn.matrix)
  suppressWarnings(object[[knn.graph.name]] <- nn.matrix)
  if (verbose) {
    message("Constructing joint SNN graph")
  }
  snn.matrix <- ComputeSNN(nn_ranked = select_nn, prune = prune.SNN )
  rownames(x = snn.matrix) <- colnames(x = snn.matrix) <- Cells(object)
  suppressWarnings(object[[snn.graph.name]] <- as(object = snn.matrix, Class = "Graph"))
  object@neighbors[[joint.nn.name]] <- joint.nn
  object@meta.data[, modality.weight.name] <- modality.weight$first.modality.weight
  return (object)
}


#' Calculate modality weights
#'
#' @param object A Seurat object
#' @param reduction.list A list of name of dimension reduction 
#' @param dims.list A list of number of dimensions to use
#' @param k.nn A list of number of nearest neighbors to use
#' @param snn.far.nn
#' @param s.nn 
#' @param prune.snn
#' @param l2.norm A list of features for different modality
#' @param sd.scale  A list of slots 
#' @param query distance metric for finding neighbors
#' @param cross.contant.list the minimal cross modality prediction similarity
#' @param sigma.idx the maximal modality score
#' @param smooth seed for bootstrapping
#' @param verbose Display messages
#'
#' @return Returns a list of similarities
#' @export
#' 
FindModalityWeights  <- function(object, 
                                    reduction.list, 
                                    dims.list, 
                                    k.nn = 20, 
                                    snn.far.nn = TRUE, 
                                    s.nn = NULL, 
                                    prune.snn = 0, 
                                    l2.norm = TRUE, 
                                    sd.scale = 1, 
                                    query = NULL, 
                                    cross.contant.list = NULL, 
                                    sigma.idx = NULL,
                                    smooth = FALSE, 
                                    verbose = TRUE
){
  if (is.null(x = s.nn)) {
    s.nn <- k.nn
  }
  if (is.null(x = sigma.idx)) {
    sigma.idx <- k.nn
  }
  if (is.null(x = cross.contant.list)) {
    cross.contant.list <- list(1e-4, 1e-4)
  }
  reduction.set <- unlist(x = reduction.list)
  names(x = reduction.list) <- names(x = dims.list) <- 
    names(x = cross.contant.list) <- reduction.set
  embeddings.list <- lapply(X = reduction.list, 
                             FUN = function(r) Embeddings(object = object, 
                                                          reduction = r)[, dims.list[[r]]])
  if (l2.norm) {
    embeddings.list.norm <- lapply(X = embeddings.list,
                                    FUN = function(embeddings) L2Norm(mat = embeddings)) 
  } else {
    embeddings.list.norm <- embeddings.list
  }
  if (is.null(x = query)) {
    query.embeddings.list.norm <- embeddings.list.norm
    query <- object
  } else { 
    if (snn.far.nn) {
      stop("query does not support to use snn to find distant neighbors")
    }
    query.embeddings.list <- lapply( X = reduction.list, 
                                     FUN = function(r) {
                                       Embeddings(object = query, reduction = r)[, dims.list[[r]]]
                                     })
    if (l2.norm) {
      query.embeddings.list <- lapply( X = query.embeddings.list,
                                            FUN = function(embeddings) L2Norm(mat = embeddings)) 
    }  
      query.embeddings.list.norm <- query.embeddings.list
  }
  if (verbose) {
    message("Finding ",k.nn ," nearest neighrbos for each modal") 
    pb <- txtProgressBar(min = 0, max = length(x = reduction.list) , style = 3)
  }
  nn.list <- lapply(X = reduction.list, 
                    FUN = function(r){
                      nn.r <- NNHelper(data = embeddings.list.norm[[r]],
                                        query = query.embeddings.list.norm[[r]],
                                        k = max(k.nn, sigma.idx, s.nn), 
                                        method = "annoy", 
                                        metric = "euclidean")
                      rownames(x = nn.r$nn.idx) <- Cells(query)
                      if (verbose) {
                        setTxtProgressBar(pb = pb, value = which(reduction.list == r))
                      }
                      return (nn.r)
                    }
  )
  sigma.nn.list <- nn.list
  if (verbose) {
    close(con = pb)
  }
  if (sigma.idx > k.nn || s.nn > k.nn) {
    nn.list <- lapply(X = nn.list, 
                      FUN = function(nn){
                        nn$nn.idx <- nn$nn.idx[, 1:k.nn]
                        nn$nn.dists <- nn$nn.dists[, 1:k.nn]
                        return (nn)
                      })
  }
  nearest_dist <-  lapply(X = reduction.list, FUN = function(r) nn.list[[r]]$nn.dists[,2])
  within_impute <- list()
  cross_impute <- list()
  # Calculating within and cross modality distance
  for (r in reduction.set) {
    reduction.norm <- paste0(r, ".norm")
    object[[ reduction.norm ]] <- CreateDimReducObject(embeddings = embeddings.list.norm[[r]],
                                                       key = paste0("norm", object[[r]]@key), 
                                                       assay = object[[r]]@assay.used )
    within_impute[[r]] <- PredictAssay(object = object, 
                                       nn.idx =  nn.list[[r]]$nn.idx,
                                       reduction = reduction.norm,
                                       dims = 1:ncol(x = embeddings.list.norm[[r]]), 
                                       verbose = FALSE,
                                       cells = Cells(query),
                                       return.assay = FALSE )
    cross_impute[[r]] <- PredictAssay(object = object,
                                     nn.idx = nn.list[[setdiff(x = reduction.set, y = r )]]$nn.idx,
                                     reduction = reduction.norm, 
                                     dims = 1:ncol(x = embeddings.list.norm[[r]]), 
                                     verbose = FALSE,
                                     cells = Cells(query),
                                     return.assay = FALSE )
  }
  
  within_impute_dist <- lapply( X = reduction.list, 
                                FUN = function(r) {
                                 r_dist <- sqrt(rowSums((query.embeddings.list.norm[[r]] -
                                                           t(within_impute[[r]]))**2))
                                 r_dist <- r_dist -  nearest_dist[[r]]
                                 r_dist[r_dist < 0] <- 0
                                 return (r_dist)
                                })
  cross_impute_dist <- lapply( X = reduction.list, 
                               FUN = function(r) {
                                 r_dist <-  sqrt(rowSums((query.embeddings.list.norm[[r]] - 
                                                            t(cross_impute[[r]]))**2))
                                 r_dist <- r_dist - nearest_dist[[r]]
                                 r_dist[r_dist < 0] <-0
                                 return(r_dist)
                               })
  # calculate kernel width
  if (snn.far.nn) {
    if (verbose) {
      message("Constructing SNN graphs for each modality by ", s.nn, " nearest neighbors") 
    }
    snn.graph.list <- lapply(X = sigma.nn.list,
                             FUN = function(nn){
                            snn.matrix <- ComputeSNN(
                            nn_ranked =  nn$nn.idx[, 1:s.nn],
                            prune = prune.snn
                               )
                           colnames(x = snn.matrix) <- rownames(x = snn.matrix) <- Cells(object)
                             return (snn.matrix)
                             })
    if (verbose) {
      message("Finding ", k.nn ," distant neighbors from snn graph") 
      pb <- txtProgressBar(min = 0, max = length(reduction.list) , style = 3)
    }
    farthest_nn_dist <- lapply(X = 1:length(x = snn.graph.list),
                       FUN = function(s) {
                         distant_nn <- ComputeSNNwidth(snn.graph = snn.graph.list[[s]],
                                                       k.nn = k.nn, 
                                                       l2.norm = FALSE,
                                                       embeddings =  embeddings.list.norm[[s]],
                                                       nearest.dist = nearest_dist[[s]] )
                         if (verbose) {
                           setTxtProgressBar(pb = pb, value = s)
                         }
                         return (distant_nn)
                       })
    names(x = farthest_nn_dist) <- unlist(x = reduction.list)
    if (verbose) {
      close(con = pb)
    }
    modality_sd.list <- lapply( X = farthest_nn_dist, 
                                FUN =  function(sd)  sd*sd.scale)
  } else {
    if (verbose) {
      message("Calculating sigma by ", sigma.idx, "th neighbor") 
    }
    modality_sd.list <- lapply(X = reduction.list , 
                                FUN =  function(r) {
      rdist <- sigma.nn.list[[r]]$nn.dists[, sigma.idx] - nearest_dist[[r]]
      rdist <- rdist * sd.scale
      return (rdist)
    })  
  }
  # Calculating within and cross modality kernel, and modalit weights
  within_impute_kernel <- lapply(X = reduction.list,
                                  FUN = function(r) {
                                    exp(-1*( within_impute_dist[[r]]/modality_sd.list[[r]] )**1) 
                                  })
  cross_impute_kernel <- lapply(X = reduction.list,
                                FUN = function(r) {
                                  exp(-1*( cross_impute_dist[[r]]/modality_sd.list[[r]] )**1) 
                                })
  params <- list( reduction.list,
                  dims.list,
                  l2.norm,
                  k.nn, 
                  sigma.idx,
                  snn.far.nn ,
                  modality_sd.list, 
                  nearest_dist)
  names(x = params) <- c("reduction.list", "dims.list", "l2.norm", "k.nn" , 
                     "sigma.idx", "snn.far.nn", "sigma.list", "nearest.dist")
  modality_score <-  lapply( X = reduction.list,
                             FUN = function(r) {
                               score = within_impute_kernel[[r]] / 
                                 ( cross_impute_kernel[[r]] + cross.contant.list[[r]] )
                               score = MinMax(data = score, min = 0, max = 200)
                             })
  if (smooth) {
    modality_score <- lapply( X = reduction.list, 
                              FUN = function(r) {
                                apply( X = nn.list[[r]]$nn.idx,
                                       MARGIN = 1, 
                                       FUN = function(nn)  mean(x = modality_score[[r]][ nn[-1]])
                                       ) 
                              })
  }
  modality1.weight <- exp(modality_score[[1]])/(exp(modality_score[[1]]) + exp(modality_score[[2]]))
  score.mat<- cbind(Reduce(f = cbind, x = within_impute_dist), 
                    Reduce(f = cbind, x = cross_impute_dist), 
                    Reduce(f = cbind, x = within_impute_kernel), 
                    Reduce(f = cbind, x = cross_impute_kernel), 
                    Reduce(f = cbind, x = modality_score))
  colnames(x = score.mat) <- c( "modality1_nn1", "modality2_nn2", 
                            "modality1_nn2",  "modality2_nn1", 
                            "modality1_nn1_kernel",  "modality2_nn2_kernel",
                            "modality1_nn2_kernel",  "modality2_nn1_kernel",
                            "modality1_score", "modality2_score")
  score.mat <- as.data.frame(x = score.mat)
  weight.list <- list(modality1.weight, params, score.mat)
  names(x = weight.list) <- c("first.modality.weight", "params", "score.matrix")
  return (weight.list)
}


#' Find farthest jarccard neighbors from snn graph
#' @param snn.graph a SNN graph
#' @param k.nn the number of neighbors is calculated from the SNN graph
#' @param far.nn determine if it find farthest neighbors in SNN graph
#' 
#' @importFrom Matrix summary
#' @importFrom dplyr group_by arrange slice select group_split
#' 
#' 
snn_nn <- function(snn.graph,
                   k.nn, 
                   far.nn = TRUE
                   ){
  if (far.nn) {
    direction <- 1
  } else {
    direction <- (-1)
  }
  edge <- summary(snn.graph)
  edge$x <- edge$x * direction
  nn.idx.snn <- edge %>% 
    group_by(j) %>% 
    arrange(x,  .by_group = TRUE)%>%
    slice(1:k.nn) %>% 
    select( -x) %>%
    group_split( .keep = FALSE) %>%
    Reduce(f = cbind, .) %>%
    t()
  rownames(x = nn.idx.snn) <- colnames(x = snn.graph)
  return (nn.idx.snn)
}

#' Find calculate mean of the farthest neighbors width from SNN graph
#' @param snn.graph
#' @param embeddings
#' @param k.nn
#' @param l2.norm
#' @param nearest.dist
#' 
#' @importFrom rdist cdist
ComputeSNNwidth <- function(snn.graph,
                            embeddings, 
                            k.nn, 
                            l2.norm = TRUE, 
                            nearest.dist = NULL) {
  if (l2.norm) {
    embeddings <- L2Norm(mat = embeddings)
  }
    # find farthest jarccard neighbors
    snn.idx <- snn_nn(snn.graph = snn.graph, k.nn = k.nn)
 if (any(is.na(x = snn.idx))) {
   warning("NA is detected, and SNN graph may be pruned")
 }
  snn.width <- t(sapply(X = 1:nrow(x = embeddings),
                                 FUN = function(x) {
                                   dist <- cdist(X = embeddings[x, , drop = FALSE], 
                                                Y = embeddings[snn.idx[x,], ],
                                                metric = "euclidean")
                                   return (dist)
                                 }))
  if (!is.null(x = nearest.dist)) {
    snn.width <- apply(X = snn.width,
           MARGIN = 2,
           FUN =  function(x) {
           dist <- x - nearest.dist
           dist[ dist < 0 ] <- 0
           return (dist)
    })
  }
  snn.width <- rowMeans(x = snn.width, na.rm = TRUE) 
  return (snn.width)
}