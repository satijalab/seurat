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

#' @param query Matrix of data to query against object. If missing, defaults to
#' object.
#' @param distance.matrix Boolean value of whether the provided matrix is a
#' distance matrix; note, for objects of class \code{dist}, this parameter will
#' be set automatically
#' @param k.param Defines k for the k-nearest neighbor algorithm
#' @param return.neighbor Return result as \code{\link{Neighbor}} object. Not
#' used with distance matrix input.
#' @param compute.SNN also compute the shared nearest neighbor graph
#' @param prune.SNN Sets the cutoff for acceptable Jaccard index when
#' computing the neighborhood overlap for the SNN construction. Any edges with
#' values less than or equal to this will be set to 0 and removed from the SNN
#' graph. Essentially sets the stringency of pruning (0 --- no pruning, 1 ---
#' prune everything).
#' @param nn.method Method for nearest neighbor finding. Options include: rann,
#' annoy
#' @param annoy.metric Distance metric for annoy. Options include: euclidean,
#' cosine, manhattan, and hamming
#' @param nn.eps Error bound when performing nearest neighbor seach using RANN;
#' default of 0.0 implies exact nearest neighbor search
#' @param verbose Whether or not to print output to the console
#' @param force.recalc Force recalculation of (S)NN.
#' @param l2.norm Take L2Norm of the data
#' @param cache.index Include cached index in returned Neighbor object
#' (only relevant if return.neighbor = TRUE)
#' @param index Precomputed index. Useful if querying new data against existing
#' index to avoid recomputing.
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
  query = NULL,
  distance.matrix = FALSE,
  k.param = 20,
  return.neighbor = FALSE,
  compute.SNN = !return.neighbor,
  prune.SNN = 1/15,
  nn.method = "annoy",
  annoy.metric = "euclidean",
  nn.eps = 0,
  verbose = TRUE,
  force.recalc = FALSE,
  l2.norm = FALSE,
  cache.index = FALSE,
  index = NULL,
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
  if (l2.norm) {
    object <- L2Norm(mat = object)
    query <- query %iff% L2Norm(mat = query)
  }
  query <- query %||% object
  # find the k-nearest neighbors for each single cell
  if (!distance.matrix) {
    if (verbose) {
      message("Computing nearest neighbor graph")
    }
    nn.ranked <- NNHelper(
      data = object,
      query = query,
      k = k.param,
      method = nn.method,
      searchtype = "standard",
      eps = nn.eps,
      metric = annoy.metric,
      cache.index = cache.index,
      index = index
    )
    if (return.neighbor) {
      if (compute.SNN) {
        warning("The SNN graph is not computed if return.neighbor is TRUE.", call. = FALSE)
      }
      return(nn.ranked)
    }
    nn.ranked <- Indices(object = nn.ranked)
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
  return.neighbor = FALSE,
  compute.SNN = !return.neighbor,
  prune.SNN = 1/15,
  nn.method = "annoy",
  annoy.metric = "euclidean",
  nn.eps = 0,
  verbose = TRUE,
  force.recalc = FALSE,
  l2.norm = FALSE,
  cache.index = FALSE,
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
    l2.norm = l2.norm,
    return.neighbor = return.neighbor,
    cache.index = cache.index,
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
  return.neighbor = FALSE,
  compute.SNN = !return.neighbor,
  prune.SNN = 1/15,
  nn.method = "annoy",
  annoy.metric = "euclidean",
  nn.eps = 0,
  verbose = TRUE,
  force.recalc = FALSE,
  l2.norm = FALSE,
  cache.index = FALSE,
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
    l2.norm = l2.norm,
    return.neighbor = return.neighbor,
    cache.index = cache.index,
    ...
  ))
}

#' @param assay Assay to use in construction of (S)NN; used only when \code{dims}
#' is \code{NULL}
#' @param features Features to use as input for building the (S)NN; used only when
#' \code{dims} is \code{NULL}
#' @param reduction Reduction to use as input for building the (S)NN
#' @param dims Dimensions of reduction to use as input
#' @param do.plot Plot SNN graph on tSNE coordinates
#' @param graph.name Optional naming parameter for stored (S)NN graph
#' (or Neighbor object, if return.neighbor = TRUE). Default is assay.name_(s)nn.
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
  return.neighbor = FALSE,
  compute.SNN = !return.neighbor,
  prune.SNN = 1/15,
  nn.method = "annoy",
  annoy.metric = "euclidean",
  nn.eps = 0,
  verbose = TRUE,
  force.recalc = FALSE,
  do.plot = FALSE,
  graph.name = NULL,
  l2.norm = FALSE,
  cache.index = FALSE,
  ...
) {
  CheckDots(...)
  if (!is.null(x = dims)) {
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
      l2.norm = l2.norm,
      return.neighbor = return.neighbor,
      cache.index = cache.index,
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
      l2.norm = l2.norm,
      return.neighbor = return.neighbor,
      cache.index = cache.index,
      ...
    )
  }
  if (length(x = neighbor.graphs) == 1) {
    neighbor.graphs <- list(nn = neighbor.graphs)
  }
  graph.name <- graph.name %||% paste0(assay, "_", names(x = neighbor.graphs))
  for (ii in 1:length(x = graph.name)) {
    if (inherits(x = neighbor.graphs[[ii]], what = "Graph")) {
      DefaultAssay(object = neighbor.graphs[[ii]]) <- assay
    }
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
# @param include.distance Include the corresponding distances
# @param index optional index object, will be recomputed if not provided
#
AnnoyNN <- function(data,
                    query = data,
                    metric = "euclidean",
                    n.trees = 50,
                    k,
                    search.k = -1,
                    include.distance = TRUE,
                    index = NULL
                    ) {
  idx <- index %||% AnnoyBuildIndex(
    data = data,
    metric = metric,
    n.trees = n.trees)
  nn <- AnnoySearch(
    index = idx,
    query = query,
    k = k,
    search.k = search.k,
    include.distance = include.distance)
  nn$idx <- idx
  nn$alg.info <- list(metric = metric, ndim = ncol(x = data))
  return(nn)
}

# Build the annoy index
#
# @param data Data to build the index with
# @param metric Distance metric; can be one of "euclidean", "cosine", "manhattan",
# "hamming"
# @param n.trees More trees gives higher precision when querying
#
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

# Search an Annoy approximate nearest neighbor index
#
# @param Annoy index, built with AnnoyBuildIndex
# @param query A set of data to be queried against the index
# @param k Number of neighbors
# @param search.k During the query it will inspect up to search_k nodes which
# gives you a run-time tradeoff between better accuracy and speed.
# @param include.distance Include the corresponding distances in the result
#
# @return A list with 'nn.idx' (for each element in 'query', the index of the
# nearest k elements in the index) and 'nn.dists' (the distances of the nearest
# k elements)
#
#' @importFrom future plan
#' @importFrom future.apply future_lapply
#
AnnoySearch <- function(index, query, k, search.k = -1, include.distance = TRUE) {
  n <- nrow(x = query)
  idx <- matrix(nrow = n,  ncol = k)
  dist <- matrix(nrow = n, ncol = k)
  convert <- methods::is(index, "Rcpp_AnnoyAngular")
  if (!inherits(x = plan(), what = "multicore")) {
    oplan <- plan(strategy = "sequential")
    on.exit(plan(oplan), add = TRUE)
  }
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

# Create an Annoy index
#
# @note Function exists because it's not exported from \pkg{uwot}
#
# @param name Distance metric name
# @param ndim Number of dimensions
#
# @return An nn index object
#
#' @importFrom methods new
#' @importFrom RcppAnnoy AnnoyAngular AnnoyManhattan AnnoyEuclidean AnnoyHamming
#
CreateAnn <- function(name, ndim) {
  return(switch(
    EXPR = name,
    cosine = new(Class = AnnoyAngular, ndim),
    manhattan = new(Class = AnnoyManhattan, ndim),
    euclidean = new(Class = AnnoyEuclidean, ndim),
    hamming = new(Class = AnnoyHamming, ndim),
    stop("BUG: unknown Annoy metric '", name, "'")
  ))
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
# @param cache.index Store algorithm index with results for reuse
# @param ... additional parameters to specific neighbor finding method
#
NNHelper <- function(data, query = data, k, method, cache.index = FALSE, ...) {
  args <- as.list(x = sys.frame(which = sys.nframe()))
  args <- c(args, list(...))
  results <- (
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
  n.ob <- Neighbor(
    nn.idx = results$nn.idx,
    nn.dist = results$nn.dists,
    alg.info = results$alg.info %||% list(),
    cell.names = rownames(x = query)
  )
  if (isTRUE(x = cache.index) && !is.null(x = results$idx)) {
    slot(object = n.ob, name = "alg.idx") <- results$idx
  }
  return(n.ob)
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

#' Find subclusters under one cluster
#'
#' @inheritParams FindClusters
#' @param cluster the cluster to be sub-clustered
#' @param subcluster.name the name of sub cluster added in the meta.data
#'
#' @return return a object with sub cluster labels in the sub-cluster.name variable
#' @export
#'
FindSubCluster <- function(
  object,
  cluster,
  graph.name,
  subcluster.name = "sub.cluster",
  resolution = 0.5,
  algorithm = 1
) {
  sub.cell <- WhichCells(object = object, idents = cluster)
  sub.graph <- as.Graph(x = object[[graph.name]][sub.cell, sub.cell])
  sub.clusters <- FindClusters(
    object = sub.graph,
    resolution = resolution,
    algorithm = algorithm
  )
  sub.clusters[, 1] <- paste(cluster,  sub.clusters[, 1], sep = "_")
  object[[subcluster.name]] <- as.character(x = Idents(object = object))
  object[[subcluster.name]][sub.cell, ] <- sub.clusters[, 1]
  return(object)
}


#' Predict value from nearest neighbors
#'
#' This function will predict expression or cell embeddings from its k nearest
#' neighbors index. For each cell, it will average its k neighbors value to get
#' its new imputed value. It can average expression value in assays and cell
#' embeddings from dimensional reductions.
#'
#' @param object The object used to calculate knn
#' @param nn.idx k near neighbour indices. A cells x k matrix.
#' @param assay Assay used for prediction
#' @param reduction Cell embedding of the reduction used for prediction
#' @param dims Number of dimensions of cell embedding
#' @param return.assay Return an assay or a predicted matrix
#' @param slot slot used for prediction
#' @param features features used for prediction
#' @param mean.function the function used to calculate row mean
#' @param seed Sets the random seed to check if the nearest neighbor is query
#' cell
#' @param verbose Print progress
#'
#' @return return an assay containing predicted expression value in the data
#' slot
#' @export
#'
PredictAssay <- function(
  object,
  nn.idx,
  assay,
  reduction = NULL,
  dims = NULL,
  return.assay = TRUE,
  slot = "scale.data",
  features = NULL,
  mean.function = rowMeans,
  seed = 4273,
  verbose = TRUE
){
  if (!inherits(x = mean.function, what = 'function')) {
    stop("'mean.function' must be a function")
  }
  if (is.null(x = reduction)) {
    reference.data <- GetAssayData(
      object = object,
      assay = assay,
      slot = slot
    )
    features <- features %||% VariableFeatures(object = object[[assay]])
    if (length(x = features) == 0) {
      features <- rownames(x = reference.data)
      if (verbose) {
        message("VariableFeatures are empty in the ", assay,
                " assay, features in the ", slot, " slot will be used" )
      }
    }
    reference.data <- reference.data[features, , drop = FALSE]
  } else {
    if (is.null(x = dims)) {
      stop("dims is empty")
    }
    reference.data <- t(x = Embeddings(object = object, reduction = reduction)[, dims])
  }
  set.seed(seed = seed)
  nn.check <- sample(x = 1:nrow(x = nn.idx), size = min(50, nrow(x = nn.idx)))
  if (all(nn.idx[nn.check, 1] == nn.check)) {
    if(verbose){
      message("The nearest neighbor is the query cell itself, and it will not be used for prediction")
    }
    nn.idx <- nn.idx[,-1]
  }
  predicted <- apply(
    X = nn.idx,
    MARGIN = 1,
    FUN = function(x) mean.function(reference.data[, x] )
  )
  colnames(x = predicted) <- Cells(x = object)
  if (return.assay) {
    predicted.assay <- CreateAssayObject(data = predicted)
    return (predicted.assay)
  } else {
    return (predicted)
  }
}

# Calculate NN distance for the given nn.idx
# @param nn.idx The nearest neighbors position index
# @param embeddings cell embeddings
# @param metric distance metric
# @param query.embeddings query cell embeddings
# @param nearest.dist The list of distance to the nearest neighbors
#
NNdist <- function(
  nn.idx,
  embeddings,
  metric = "euclidean",
  query.embeddings = NULL,
  nearest.dist = NULL
) {
  if (!is.list(x = nn.idx)) {
    nn.idx <- lapply(X = 1:nrow(x = nn.idx), FUN = function(x) nn.idx[x, ])
  }
  query.embeddings <- query.embeddings %||% embeddings
  nn.dist <- fast_dist(
    x = query.embeddings,
    y = embeddings,
    n = nn.idx
  )
  if (!is.null(x = nearest.dist)) {
    nn.dist <- lapply(
      X = 1:nrow(x = query.embeddings),
      FUN = function(x) {
        r_dist = nn.dist[[x]] - nearest.dist[x]
        r_dist[r_dist < 0] <- 0
        return(r_dist)
      }
    )
  }
  return(nn.dist)
}


# Find multimodal neighbors
#
# @param object The object used to calculate knn
# @param query The query object when query and reference are different
# @param modality.weight A \code{\link{ModalityWeights}} object generated by
# \code{\link{FindModalityWeights}}
# @param k.nn .Number of nearest multi-model neighbors to compute
# @param reduction.list A list of reduction name
# @param dims.list A list of dimentions used for the reduction
# @param knn.range The number of approximate neighbors to compute
# @param kernel.power The power for the exponential kernel
# @param nearest.dist The list of distance to the nearest neighbors
# @param sigma.list The list of kernel width
# @param l2.norm Perform L2 normalization on the cell embeddings after
# dimensional reduction
# @param verbose Print output to the console
# @importFrom pbapply pblapply
# @return return a list containing nn index and nn multimodal distance
#
MultiModalNN <- function(
  object,
  query = NULL,
  modality.weight = NULL,
  k.nn = NULL,
  reduction.list = NULL,
  dims.list = NULL,
  knn.range = 200,
  kernel.power = 1,
  nearest.dist = NULL,
  sigma.list = NULL,
  l2.norm =  NULL,
  verbose = TRUE
){
  my.lapply <- ifelse(
    test = verbose,
    yes = pblapply,
    no = lapply
  )
  k.nn <-  k.nn %||% slot(object = modality.weight, name = "params")$k.nn
  reduction.list <- reduction.list %||%
    slot(object = modality.weight, name = "params")$reduction.list
  dims.list = dims.list %||%
    slot(object = modality.weight, name = "params")$dims.list
  nearest.dist = nearest.dist %||%
    slot(object = modality.weight, name = "params")$nearest.dist
  sigma.list =sigma.list %||%
    slot(object = modality.weight, name = "params")$sigma.list
  l2.norm = l2.norm %||%
    slot(object = modality.weight, name = "params")$l2.norm
  fmw <- slot(object = modality.weight, name = "first.modality.weight")
  modality.weight.value <- list(fmw, 1 - fmw)
  names(x = modality.weight.value) <- unlist(x = reduction.list)
  if (inherits(x = object, what = "Seurat")) {
    reduction_embedding <- lapply(
      X = 1:length(x = reduction.list),
      FUN = function(x) {
        Embeddings(object = object, reduction = reduction.list[[x]])[, dims.list[[x]]]
      }
    )
  } else {
    reduction_embedding <- object
  }
  if (is.null(x = query)) {
    query.reduction_embedding <- reduction_embedding
    query <- object
  } else {
    if (inherits(x = object, what = "Seurat")) {
      query.reduction_embedding <- lapply(
        X = 1:length(x = reduction.list),
        FUN = function(x) {
          Embeddings(object = query, reduction = reduction.list[[x]] )[, dims.list[[x]]]
        }
      )
    } else {
      query.reduction_embedding <- query
    }
  }
  if (l2.norm) {
    query.reduction_embedding <- lapply(
      X = query.reduction_embedding,
      FUN = function(x)  L2Norm(mat = x)
    )
    reduction_embedding <- lapply(
      X = reduction_embedding,
      FUN = function(x) L2Norm(mat = x)
    )
  }
  query.cell.num <- nrow(x = query.reduction_embedding[[1]])
  reduction.num <- length(x = query.reduction_embedding)
  if (verbose) {
    message("Finding multi-modal neighbors")
  }
  redunction_nn <- my.lapply(
    X = 1:reduction.num,
    FUN = function(x) {
      nn_x <- NNHelper(
        data = reduction_embedding[[x]],
        query = query.reduction_embedding[[x]],
        k = knn.range,
        method = 'annoy',
        metric = "euclidean"
      )
      return (nn_x)
    }
  )
  # union of rna and adt nn, remove itself from neighobors
  redunction_nn <- lapply(
    X = redunction_nn,
    FUN = function(x)  Indices(object = x)[, -1]
  )
  nn_idx <- lapply(
    X = 1:query.cell.num ,
    FUN = function(x) {
      Reduce(
        f = union,
        x = lapply(
          X = redunction_nn,
          FUN = function(y) y[x, ]
        )
      )
    }
  )
  if (verbose) {
    message("Calculating distance of multi-modal neighbors")
  }
  # calculate euclidean distance of all neighbors
  nn_dist <- my.lapply(
    X = 1:reduction.num,
    FUN = function(r) {
      nndist <- NNdist(
        nn.idx = nn_idx,
        embeddings = reduction_embedding[[r]],
        query.embeddings = query.reduction_embedding[[r]],
        nearest.dist = nearest.dist[[r]]
      )
      return(nndist)
   }
  )
  # modality weighted distance
  if (length(x = sigma.list[[1]]) == 1) {
    sigma.list <- lapply(X = sigma.list, FUN = function(x) rep(x = x, ncol(x = object)))
  }
  nn_weighted_dist <- lapply(
    X = 1:reduction.num,
    FUN = function(r) {
      lapply(
        X = 1:query.cell.num,
        FUN = function(x) {
          exp(-1*(nn_dist[[r]][[x]] / sigma.list[[r]][x] ) ** kernel.power) * modality.weight.value[[r]][x]
        }
      )
    }
  )
  nn_weighted_dist <- sapply(
    X = 1:query.cell.num,
    FUN =  function(x) {
      Reduce(
        f = "+",
        x = lapply(
          X = 1:reduction.num,
          FUN = function(r) nn_weighted_dist[[r]][[x]]
        )
      )
    }
  )
  # select k nearest joint neighbors
  select_order <- lapply(
    X = nn_weighted_dist,
    FUN = function(dist) {
      order(dist, decreasing = TRUE)
  })
  select_nn <- t(x = sapply(
    X = 1:query.cell.num,
    FUN = function(x) nn_idx[[x]][select_order[[x]]][1:k.nn]
    )
  )
  select_dist <- t(x = sapply(
    X = 1:query.cell.num,
    FUN = function(x) nn_weighted_dist[[x]][select_order[[x]]][1:k.nn])
  )
  select_dist <- sqrt(x = (1 - select_dist) / 2)
  weighted.nn <- Neighbor(
    nn.idx = select_nn,
    nn.dist = select_dist,
    alg.info = list(),
    cell.names = Cells(x = query)
  )
  return(weighted.nn)
}


#' Construct multimodal neighbors, KNN and SNN Graph
#'
#' This function will construct multimodal neighbors, Kth Nearest Neighbors
#' (KNN) and Shared Nearest Neighbor (SNN) Graphs. According to the input
#' \code{modality.weight}, it constructs a cell-specific weighted joint kernel.
#' Then, for each cell, it will find \code{knn.range}s individual modal
#' neighbors, and get the union of those neighbors. Next, it will find
#' \code{k.nn} multimodal neighbors by the weighted joint kernel. Given the set
#' of multimodal neighbors, we construct its KNN and SNN Graph.
#'
#' @param object A Seurat object
#' @param reduction.list A list of name of dimension reduction
#' @param dims.list A list of number of dimensions to use
#' @param modality.weight.name The variable name of first modality weights
#' stored in the meta.data.
#' @param k.nn the number of multi-modal neighbors computed
#' @param l2.norm Perform L2 normalization on the cell embeddings after
#' dimensional reduction
#' @param sd.scale  The scaling factor for kernel width, and the default is 1.
#' @param cross.contant.list the minimal cross-modality prediction similarity
#' used in the modality score calculation.
#' @param smooth Smoothing modality score across each individual modality
#' neighbors.
#' @param prune.SNN Sets the cutoff for acceptable Jaccard index when computing
#' the neighborhood overlap for the SNN construction
#' @param knn.graph.name The name of multimodal knn graph
#' @param snn.graph.name The name of multimodal snn graph
#' @param weighted.nn.name The name of multimodal neighbors

#' @param knn.range The number of approximate neighbors to compute
#' @param modality.weight A \code{\link{ModalityWeights}} object generated by
#' \code{FindModalityWeights}
#' @param weighted.graph Add consider neighbor distance as the edges to
#' construct KNN graph
#' @param return.intermediate Return intermediate results in the FindModalityWeights.
#' It will be stored in the misc as a ModalityClass.
#' @param verbose Print progress bars and output
#'
#' @return return an object containing multimodal KNN, SNN and neighbors
#' @export

FindMultiModalNeighbors  <- function(
  object,
  reduction.list,
  dims.list,
  modality.weight.name = "first.modality.weight",
  k.nn = 20,
  l2.norm = TRUE,
  sd.scale = 1,
  cross.contant.list = list(1e-4, 1e-4),
  smooth = FALSE,
  knn.range = 200,
  knn.graph.name = "wknn",
  snn.graph.name = "wsnn",
  weighted.nn.name = "weighted.nn",
  modality.weight = NULL,
  prune.SNN = 1/15,
  weighted.graph = FALSE,
  return.intermediate = FALSE,
  verbose = TRUE
) {
  if (is.null(x = modality.weight)) {
    modality.weight <- FindModalityWeights(
      object = object,
      reduction.list = reduction.list,
      dims.list = dims.list,
      k.nn = k.nn,
      sd.scale = sd.scale,
      l2.norm = l2.norm,
      cross.contant.list = cross.contant.list,
      smooth = smooth,
      verbose = verbose
   )
  }
  k.nn <- k.nn %||% slot(object = modality.weight, name = "params")$k.nn
  first.assay <- slot(object = modality.weight, name = "modality.assay")[1]
  weighted.nn <- MultiModalNN(
    object = object,
    k.nn = k.nn,
    modality.weight = modality.weight,
    knn.range = knn.range,
    verbose = verbose
  )
  select_nn <- Indices(object = weighted.nn)
  select_nn_dist <- Distances(object = weighted.nn)
  # compute KNN graph
  if (weighted.graph) {
    if (verbose) {
      message("Constructing joint weighted knn graph")
    }
    select_nn_dist <- t(x = apply(
      X = select_nn_dist,
      MARGIN = 1,
      FUN = function(x) log2(k.nn) * x / sum(x))
    )
    nn.matrix <- sparseMatrix(
      i = 1:ncol(x = object),
      j = 1:ncol(x = object),
      x = 1
    )
    for (i in 1:ncol(x = object)) {
      nn.matrix[i, select_nn[i, ]] <- select_nn_dist[i, ]
    slot(object = weighted.nn, name = "nn.dist") <- select_nn_dist
    }
  } else {
    if (verbose) {
      message("Constructing multi-modal KNN graph")
    }
    j <- as.numeric(x = t(x = select_nn ))
    i <- ((1:length(x = j)) - 1) %/% k.nn + 1
    nn.matrix <- sparseMatrix(
      i = i,
      j = j,
      x = 1,
      dims = c(ncol(x = object), ncol(x = object))
    )
    diag(x = nn.matrix) <- 1
  }
  rownames(x = nn.matrix) <-  colnames(x = nn.matrix) <- colnames(x = object)
  nn.matrix <- nn.matrix + t(x = nn.matrix) - t(x = nn.matrix) * nn.matrix
  nn.matrix <- as.Graph(x = nn.matrix)
  slot(object = nn.matrix, name = "assay.used") <- first.assay
  object[[knn.graph.name]] <- nn.matrix

  # compute SNN graph
  if (verbose) {
    message("Constructing multi-modal SNN graph")
  }
  snn.matrix <- ComputeSNN(nn_ranked = select_nn, prune = prune.SNN)
  rownames(x = snn.matrix) <- colnames(x = snn.matrix) <- Cells(x = object)
  snn.matrix <- as.Graph(x = snn.matrix )
  slot(object = snn.matrix, name = "assay.used") <- first.assay
  object[[snn.graph.name]] <- snn.matrix

  # add neighbors and modality weights
  object[[weighted.nn.name]] <- weighted.nn
  object[[modality.weight.name]] <- slot(object = modality.weight, name = "first.modality.weight")

  # add command log
   modality.weight.command <- slot(object = modality.weight, name = "command")
   slot(object = modality.weight.command, name = "assay.used") <- first.assay
   modality.weight.command.name <- slot(object = modality.weight.command, name = "name")
   object[[modality.weight.command.name]] <- modality.weight.command
   command <- LogSeuratCommand(object = object, return.command = TRUE)
   slot(object = command, name = "params")$modality.weight <- NULL
   slot(object = command, name = "assay.used") <- first.assay
   command.name <- slot(object = command, name = "name")
   object[[command.name]] <- command

   if (return.intermediate) {
     Misc(object = object, slot = "modality.weight") <- modality.weight
   }
   return (object)
}


# Calculate modality weights
#
# This function calculates cell-specific modality weights which are used to
# construct the multimodal kernel to find multimodal neighbors. It finds
# neighbors from each modality and performs within- and cross- modality
# prediction to calculate modality weights.
#
# @param object A Seurat object
# @param reduction.list A list of name of dimension reduction
# @param dims.list A list of number of dimensions to use
# @param k.nn How many neighbors (k) to use
# @param snn.far.nn Use SNN to find farthest neighbors to calculate
# the kernel width
# @param s.nn How many neighbors (k) to use from the SNN graph
# @param prune.SNN Sets the cutoff for acceptable Jaccard index when computing
# the neighborhood overlap for the SNN construction.
# @param l2.norm Perform L2 normalization on the cell embeddings after
# dimensional reduction
# @param sd.scale  The scaling factor for kernel width, and the default is 1.
# @param query A Seurat object used as the query when query and reference
# objects are different. snn.far.nn does not support for query object.
# @param cross.contant.list the minimal cross-modality prediction similarity
# used in the modality score calculation.
# @param sigma.idx Use sigma.idx-th neighbor's distance as the kernel width.
# When snn.far.nn is TRUE, this parameter is not used
# @param smooth Smoothing modality score across each individual modality
# neighbors.
# @param verbose Display messages
# @importFrom pbapply pblapply
# @return Returns a \code{ModalityWeights} object that can be used as input to
# \code{\link{FindMultiModalNeighbors}}
#
#' @importFrom pbapply pblapply
#
FindModalityWeights  <- function(
  object,
  reduction.list,
  dims.list,
  k.nn = 20,
  snn.far.nn = TRUE,
  s.nn = k.nn,
  prune.SNN = 0,
  l2.norm = TRUE,
  sd.scale = 1,
  query = NULL,
  cross.contant.list = list(1e-4, 1e-4),
  sigma.idx = k.nn,
  smooth = FALSE,
  verbose = TRUE
) {
  my.lapply <- ifelse(
    test = verbose,
    yes = pblapply,
    no = lapply
  )
  reduction.set <- unlist(x = reduction.list)
  names(x = reduction.list) <- names(x = dims.list) <-
    names(x = cross.contant.list) <- reduction.set
  embeddings.list <- lapply(
    X = reduction.list,
    FUN = function(r) Embeddings(object = object, reduction = r)[, dims.list[[r]]]
  )
  if (l2.norm) {
    embeddings.list.norm <- lapply(
      X = embeddings.list,
      FUN = function(embeddings) L2Norm(mat = embeddings)
    )
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
    query.embeddings.list <- lapply(
      X = reduction.list,
      FUN = function(r) {
        Embeddings(object = query, reduction = r)[, dims.list[[r]]]
      }
    )
    if (l2.norm) {
      query.embeddings.list <- lapply(
        X = query.embeddings.list,
        FUN = function(embeddings) L2Norm(mat = embeddings)
      )
    }
    query.embeddings.list.norm <- query.embeddings.list
  }
  if (verbose) {
    message("Finding ", k.nn, " nearest neighbors for each modality.")
  }
  nn.list <- my.lapply(
    X = reduction.list,
    FUN = function(r) {
      nn.r <- NNHelper(
        data = embeddings.list.norm[[r]],
        query = query.embeddings.list.norm[[r]],
        k = max(k.nn, sigma.idx, s.nn),
        method = "annoy",
        metric = "euclidean"
      )
      return(nn.r)
    }
  )
  sigma.nn.list <- nn.list

  if (sigma.idx > k.nn || s.nn > k.nn) {
    nn.list <- lapply(
      X = nn.list,
      FUN = function(nn){
        slot(object = nn, name = "nn.idx") <- Indices(object = nn)[, 1:k.nn]
        slot(object = nn, name = "nn.dists") <- Distances(object = nn)[, 1:k.nn]
        return(nn)
      }
    )
  }
  nearest_dist <- lapply(X = reduction.list, FUN = function(r) Distances(object = nn.list[[r]])[, 2])
  within_impute <- list()
  cross_impute <- list()
  # Calculating within and cross modality distance
  for (r in reduction.set) {
    reduction.norm <- paste0(r, ".norm")
    object[[ reduction.norm ]] <- CreateDimReducObject(
      embeddings = embeddings.list.norm[[r]],
      key = paste0("norm", Key(object = object[[r]])),
      assay = DefaultAssay(object = object[[r]])
    )
    within_impute[[r]] <- PredictAssay(
      object = object,
      nn.idx =  Indices(object = nn.list[[r]]),
      reduction = reduction.norm,
      dims = 1:ncol(x = embeddings.list.norm[[r]]),
      verbose = FALSE,
      return.assay = FALSE
    )
    cross_impute[[r]] <- PredictAssay(
      object = object,
      nn.idx = Indices(object = nn.list[[setdiff(x = reduction.set, y = r )]]),
      reduction = reduction.norm,
      dims = 1:ncol(x = embeddings.list.norm[[r]]),
      verbose = FALSE,
      return.assay = FALSE
    )
  }
  within_impute_dist <- lapply(
    X = reduction.list,
    FUN = function(r) {
     r_dist <- sqrt(x = rowSums(x = (query.embeddings.list.norm[[r]] - t(x = within_impute[[r]])) ** 2))
     r_dist <- r_dist -  nearest_dist[[r]]
     r_dist[r_dist < 0] <- 0
     return(r_dist)
    }
  )
  cross_impute_dist <- lapply(
    X = reduction.list,
    FUN = function(r) {
      r_dist <-  sqrt(x = rowSums(x = (query.embeddings.list.norm[[r]] - t(x = cross_impute[[r]])) ** 2))
      r_dist <- r_dist - nearest_dist[[r]]
      r_dist[r_dist < 0] <-0
      return(r_dist)
    }
  )
  # calculate kernel width
  if (snn.far.nn) {
    if (verbose) {
      message("Constructing SNN graphs for each modality by ", s.nn, " nearest neighbors")
    }
    snn.graph.list <- lapply(
      X = sigma.nn.list,
      FUN = function(nn) {
        snn.matrix <- ComputeSNN(
          nn_ranked =  Indices(object = nn)[, 1:s.nn],
          prune = prune.SNN
        )
        colnames(x = snn.matrix) <- rownames(x = snn.matrix) <- Cells(x = object)
        return (snn.matrix)
      }
    )
    if (verbose) {
      message("Finding ", k.nn, " distant neighbors from snn graph")
    }
    farthest_nn_dist <- my.lapply(
      X = 1:length(x = snn.graph.list),
      FUN = function(s) {
        distant_nn <- ComputeSNNwidth(
          snn.graph = snn.graph.list[[s]],
          k.nn = k.nn,
          l2.norm = FALSE,
          embeddings =  embeddings.list.norm[[s]],
          nearest.dist = nearest_dist[[s]]
        )
        return (distant_nn)
      }
    )
    names(x = farthest_nn_dist) <- unlist(x = reduction.list)
    modality_sd.list <- lapply(
      X = farthest_nn_dist,
      FUN =  function(sd)  sd * sd.scale
    )
  } else {
    if (verbose) {
      message("Calculating sigma by ", sigma.idx, "th neighbor")
    }
    modality_sd.list <- lapply(
      X = reduction.list ,
      FUN =  function(r) {
        rdist <- Distances(object = sigma.nn.list[[r]])[, sigma.idx] - nearest_dist[[r]]
        rdist <- rdist * sd.scale
        return (rdist)
      }
    )
  }
  # Calculating within and cross modality kernel, and modalit weights
  within_impute_kernel <- lapply(
    X = reduction.list,
    FUN = function(r) {
      exp(-1 * (within_impute_dist[[r]] / modality_sd.list[[r]]) )
    }
  )
  cross_impute_kernel <- lapply(
    X = reduction.list,
    FUN = function(r) {
      exp(-1 * (cross_impute_dist[[r]] / modality_sd.list[[r]]) )
    }
  )
  params <- list(
    "reduction.list" = reduction.list,
    "dims.list" = dims.list,
    "l2.norm" = l2.norm,
    "k.nn" = k.nn,
    "sigma.idx" = sigma.idx,
    "snn.far.nn" = snn.far.nn ,
    "sigma.list" = modality_sd.list,
    "nearest.dist" = nearest_dist
  )
  modality_score <-  lapply(
    X = reduction.list,
    FUN = function(r) {
      score <- within_impute_kernel[[r]] / (cross_impute_kernel[[r]] + cross.contant.list[[r]])
      score <- MinMax(data = score, min = 0, max = 200)
    }
  )
  if (smooth) {
    modality_score <- lapply(
      X = reduction.list,
      FUN = function(r) {
        apply(
          X = Indices(object = nn.list[[r]]),
          MARGIN = 1,
          FUN = function(nn)  mean(x = modality_score[[r]][nn[-1]])
        )
      }
    )
  }
  modality1.weight <- exp(x = modality_score[[1]]) / (exp(x = modality_score[[1]]) + exp(x = modality_score[[2]]))
  score.mat<- cbind(
    Reduce(f = cbind, x = within_impute_dist),
    Reduce(f = cbind, x = cross_impute_dist),
    Reduce(f = cbind, x = within_impute_kernel),
    Reduce(f = cbind, x = cross_impute_kernel),
    Reduce(f = cbind, x = modality_score)
  )
  colnames(x = score.mat) <- c(
    "modality1_nn1", "modality2_nn2", "modality1_nn2",  "modality2_nn1",
    "modality1_nn1_kernel", "modality2_nn2_kernel", "modality1_nn2_kernel",
    "modality2_nn1_kernel", "modality1_score", "modality2_score"
  )
  score.mat <- as.data.frame(x = score.mat)
  # unlist the input parameters
  command <- LogSeuratCommand(object = object, return.command = TRUE)
  command@params <- lapply(X =  command@params , FUN = function (l) unlist(x = l))
  modality.assay <- sapply(
    X = reduction.list ,
    FUN = function (r) slot(object[[r]], name = "assay.used")
  )
  modality.weights <- new(
    Class = "ModalityWeights",
    first.modality.weight = modality1.weight,
    modality.assay = modality.assay,
    params = params,
    score.matrix = score.mat,
    command = command
  )
  return (modality.weights)
}

# Calculate mean distance of the farthest neighbors from SNN graph
#
# This function will compute the average distance of the farthest k.nn
# neighbors with the lowest nonzero SNN edge weight. First, for each cell it
# finds the k.nn neighbors with the smallest edge weight. If there are multiple
# cells with the same edge weight at the k.nn-th index, consider all of those
# cells in the next step. Next, it computes the euclidean distance to all k.nn
# cells in the space defined by the embeddings matrix and returns the average
# distance to the farthest k.nn cells.
#
# @param snn.graph An SNN graph
# @param embeddings The cell embeddings used to calculate neighbor distances
# @param k.nn The number of neighbors to calculate
# @param l2.norm Perform L2 normalization on the cell embeddings
# @param nearest.dist The vector of distance to the nearest neighbors to
# subtract off from distance calculations
#
#
ComputeSNNwidth <- function(
  snn.graph,
  embeddings,
  k.nn,
  l2.norm = TRUE,
  nearest.dist = NULL
) {
  if (l2.norm) {
    embeddings <- L2Norm(mat = embeddings)
  }
  nearest.dist <- nearest.dist %||% rep(x = 0, times = ncol(x = snn.graph))
  if (length(x = nearest.dist) != ncol(x = snn.graph)) {
    stop("Please provide a vector for nearest.dist that has as many elements as",
         " there are columns in the snn.graph (", ncol(x = snn.graph), ").")
  }
  snn.width <- SNN_SmallestNonzero_Dist(
    snn = snn.graph,
    mat = embeddings,
    n = k.nn,
    nearest_dist = nearest.dist
  )
  return (snn.width)
}
