#' @include zzz.R
#' @include generics.R
#' @importFrom rlang enquo is_quosure quo_get_env quo_get_expr
#'
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Generics
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' Sketch Data
#'
#' This function uses sketching methods to downsample high-dimensional single-cell RNA expression data,
#' which can help with scalability for large datasets.
#'
#' @param object A Seurat object.
#' @param assay Assay name. Default is NULL, in which case the default assay of the object is used.
#' @param ncells A positive integer indicating the number of cells to sample for the sketching. Default is 5000.
#' @param sketched.assay Sketched assay name. A  sketch assay is created or overwrite with the sketch data. Default is 'sketch'.
#' @param method  Sketching method to use. Can be 'LeverageScore' or 'Uniform'.
#'               Default is 'LeverageScore'.
#' @param var.name A metadata column name to store the leverage scores. Default is 'leverage.score'.
#' @param over.write whether to overwrite existing column in the metadata. Default is FALSE.
#' @param seed A positive integer for the seed of the random number generator. Default is 123.
#' @param cast The type to cast the resulting assay to. Default is 'dgCMatrix'.
#' @param verbose Print progress and diagnostic messages
#' @param ... Arguments passed to other methods
#'
#' @return A Seurat object with the sketched data added as a new assay.
#'
#' @importFrom SeuratObject CastAssay Key Key<- Layers
#'
#' @export
#' @concept sketching
#'
SketchData <- function(
  object,
  assay = NULL,
  ncells = 5000L,
  sketched.assay = 'sketch',
  method = c('LeverageScore', 'Uniform'),
  var.name = "leverage.score",
  over.write = FALSE,
  seed = 123L,
  cast = 'dgCMatrix',
  verbose = TRUE,
  features = NULL,
  ...
) {
  assay <- assay[1L] %||% DefaultAssay(object = object)
  assay <- match.arg(arg = assay, choices = Assays(object = object))
  method <- match.arg(arg = method)
  if (sketched.assay == assay) {
    abort(message = "Cannot overwrite existing assays")
  }
  if (sketched.assay %in% Assays(object = object)) {
    if (sketched.assay == DefaultAssay(object = object)) {
      DefaultAssay(object = object) <- assay
    }
    object[[sketched.assay]] <- NULL
  }
  if (!over.write) {
    var.name <- CheckMetaVarName(object = object, var.name = var.name)
  }

  if (method == 'LeverageScore') {
    if (verbose) {
      message("Calcuating Leverage Score")
    }
    object <- LeverageScore(
      object = object,
      assay = assay,
      var.name = var.name,
      over.write = over.write,
      seed = seed,
      verbose = FALSE,
      features = features,
      ...
    )
  } else if (method == 'Uniform') {
    if (verbose) {
      message("Uniformly sampling")
    }
    object[[var.name]] <- 1
  }
  leverage.score <- object[[var.name]]
  layers.data <- Layers(object = object[[assay]], search = 'data')
  cells <- list()
  for (i in seq_along(layers.data)){
    set.seed(seed = seed) # does this need to be set in the forloop? is it getting updated somehow 
    lyr <- layers.data[i]
    if (length(ncells) == 1) { # use the same number of cells per layer 
      ncells.lyr <- ncells
    } else {
      ncells.lyr <- ncells[i]
    }
    lcells <- Cells(x = object[[assay]], layer = lyr)
    if (length(x = lcells) < ncells.lyr) {
      cells[[i]] <- lcells
    } else {
      cells[[i]] <- sample(
        x = lcells,
        size = ncells.lyr,
        prob = leverage.score[lcells,]
      )
    }
    seed = seed
  }
  # cells <- lapply(
  #   X = seq_along(along.with = layers.data),
  #   FUN = function(i, seed) {
  #     set.seed(seed = seed)
  #     lcells <- Cells(x = object[[assay]], layer = layers.data[i])
  #     if (length(x = lcells) < ncells) {
  #       return(lcells)
  #     }
  #     return(sample(
  #       x = lcells,
  #       size = ncells,
  #       prob = leverage.score[lcells,]
  #     ))
  #   },
  #   seed = seed
  # )
  sketched <- suppressWarnings(expr = subset(
    x = object[[assay]],
    cells = unlist(cells),
    layers = Layers(object = object[[assay]], search = c('counts', 'data'))
  ))
  for (lyr in layers.data) {
    try(
      expr = VariableFeatures(object = sketched, method = "sketch", layer = lyr) <-
        VariableFeatures(object = object[[assay]], layer = lyr),
      silent = TRUE
    )
  }
  if (!is.null(x = cast) && inherits(x = sketched, what = 'Assay5')) {
    sketched <- CastAssay(object = sketched, to = cast, ...)
  }
  Key(object = sketched) <- Key(object = sketched.assay, quiet = TRUE)
  object[[sketched.assay]] <- sketched
  DefaultAssay(object = object) <- sketched.assay
  return(object)
}


#' Project full data to the sketch assay
#'
#'
#' This function allows projection of high-dimensional single-cell RNA expression data from a full dataset
#' onto the lower-dimensional embedding of the sketch of the dataset.
#'
#' @param object A Seurat object.
#' @param assay Assay name for the full data. Default is 'RNA'.
#' @param sketched.assay Sketched assay name to project onto. Default is 'sketch'.
#' @param sketched.reduction Dimensional reduction results of the sketched assay to project onto.
#' @param full.reduction Dimensional reduction name for the projected full dataset.
#' @param dims Dimensions to include in the projection.
#' @param normalization.method Normalization method to use. Can be 'LogNormalize' or 'SCT'.
#'               Default is 'LogNormalize'.
#' @param refdata An optional list for label transfer from sketch to full data. Default is NULL.
#'        Similar to refdata in `MapQuery`
#' @param k.weight Number of neighbors to consider when weighting labels for transfer. Default is 50.
#' @param umap.model An optional pre-computed UMAP model. Default is NULL.
#' @param recompute.neighbors Whether to recompute the neighbors for label transfer. Default is FALSE.
#' @param recompute.weights Whether to recompute the weights for label transfer. Default is FALSE.
#' @param verbose Print progress and diagnostic messages.
#'
#' @return A Seurat object with the full data projected onto the sketched dimensional reduction results.
#' The projected data are stored in the specified full reduction.
#'
#' @export
#' @concept sketching
#'
ProjectData <- function(
  object,
  assay = 'RNA',
  sketched.assay = 'sketch',
  sketched.reduction,
  full.reduction,
  dims,
  normalization.method = c("LogNormalize", "SCT"),
  refdata = NULL,
  k.weight = 50,
  umap.model = NULL,
  recompute.neighbors = FALSE,
  recompute.weights = FALSE,
  verbose = TRUE
) {
  if (!full.reduction %in% Reductions(object)) {
    if (verbose) {
      message(full.reduction, ' is not in the object.'
              ,' Data from all cells will be projected to ', sketched.reduction)
    }
    proj.emb <- ProjectCellEmbeddings(
      query = object,
      reference = object,
      query.assay = assay,
      dims = dims,
      normalization.method = normalization.method,
      reference.assay = sketched.assay,
      reduction = sketched.reduction,
      verbose = verbose)
    object[[full.reduction]] <- CreateDimReducObject(
      embeddings = proj.emb,
      assay = assay,
      key = Key(object = full.reduction, quiet = TRUE)
    )
  }
  object <- TransferSketchLabels(
    object = object,
    sketched.assay = sketched.assay,
    reduction = full.reduction,
    dims = dims,
    k = k.weight,
    refdata = refdata,
    reduction.model = umap.model,
    recompute.neighbors = recompute.neighbors,
    recompute.weights = recompute.weights,
    verbose = verbose)
  return(object)
}


#' Transfer data from sketch data to full data
#'
#' This function transfers cell type labels from a sketched dataset to a full dataset
#' based on the similarities in the lower dimensional space.
#'
#' @param object A Seurat object.
#' @param sketched.assay Sketched assay name. Default is 'sketch'.
#' @param reduction Dimensional reduction name to use for label transfer.
#' @param dims An integer vector indicating which dimensions to use for label transfer.
#' @param refdata A list of character strings indicating the metadata columns containing labels to transfer. Default is NULL.
#'                Similar to refdata in `MapQuery`
#' @param k Number of neighbors to use for label transfer. Default is 50.
#' @param reduction.model Dimensional reduction model to use for label transfer. Default is NULL.
#' @param neighbors An object storing the neighbors found during the sketching process. Default is NULL.
#' @param recompute.neighbors Whether to recompute the neighbors for label transfer. Default is FALSE.
#' @param recompute.weights Whether to recompute the weights for label transfer. Default is FALSE.
#' @param verbose Print progress and diagnostic messages
#'
#' @return A Seurat object with transferred labels stored in the metadata. If a UMAP model is provided,
#' the full data are also projected onto the UMAP space, with the results stored in a new reduction, full.`reduction.model`
#'
#' @export
#' @concept sketching
#'
TransferSketchLabels <- function(
  object,
  sketched.assay = 'sketch',
  reduction,
  dims,
  refdata = NULL,
  k = 50,
  reduction.model = NULL,
  neighbors = NULL,
  recompute.neighbors = FALSE,
  recompute.weights = FALSE,
  verbose = TRUE
){
  full_sketch.nn <- neighbors %||% Tool(
    object = object,
    slot = 'TransferSketchLabels'
  )$full_sketch.nn
  full_sketch.weight <- Tool(
    object = object,
    slot = 'TransferSketchLabels'
  )$full_sketch.weight

  compute.neighbors <- is.null(x = full_sketch.nn) ||
    !all(Cells(full_sketch.nn) == Cells(object[[reduction]])) ||
    max(Indices(full_sketch.nn)) >  ncol(object[[sketched.assay]]) ||
    !identical(x = full_sketch.nn@alg.info$dims, y =  dims) ||
    !identical(x = full_sketch.nn@alg.info$reduction, y =  reduction) ||
    recompute.neighbors

  compute.weights <- is.null(x = full_sketch.weight) ||
    !all(colnames(full_sketch.weight) == Cells(object[[reduction]])) ||
    !all(rownames(full_sketch.weight) == colnames(object[[sketched.assay]]))  ||
    recompute.weights ||
    recompute.neighbors

  if (compute.neighbors) {
    if (verbose) {
      message("Finding sketch neighbors")
    }
    full_sketch.nn <- NNHelper(
      query = Embeddings(object[[reduction]])[, dims],
      data = Embeddings(object[[reduction]])[colnames(object[[sketched.assay]]), dims],
      k = k,
      method = "annoy"
    )
    slot(object = full_sketch.nn, name = 'alg.info')$dims <- dims
    slot(object = full_sketch.nn, name = 'alg.info')$reduction <- reduction
  }
  if (compute.weights) {
    if (verbose) {
      message("Finding sketch weight matrix")
    }
    full_sketch.weight <- FindWeightsNN(
      nn.obj = full_sketch.nn,
      query.cells = Cells(object[[reduction]]),
      reference.cells = colnames(object[[sketched.assay]]),
      verbose = verbose)
    rownames(full_sketch.weight) <- colnames(object[[sketched.assay]])
    colnames(full_sketch.weight) <- Cells(object[[reduction]])
  }
  slot(
    object = object, name = 'tools'
    )$TransferSketchLabels$full_sketch.nn <- full_sketch.nn
  slot(
    object = object, name = 'tools'
    )$TransferSketchLabels$full_sketch.weight <- full_sketch.weight

  if (!is.null(refdata)) {
    if (length(refdata) == 1  & is.character(refdata)) {
      refdata <- list(refdata)
      names(refdata) <- unlist(refdata)
    }
    if (verbose) {
      message("Transfering refdata from sketch")
    }
    for (rd in 1:length(x = refdata)) {
      if (isFALSE(x = refdata[[rd]])) {
        transfer.results[[rd]] <- NULL
        next
      }
      rd.name <- names(x = refdata)[rd]
      label.rd <- refdata[[rd]]
      ## FetchData not work
      if (!label.rd %in% colnames( object[[]])) {
        stop(label.rd, ' is not in the meta.data')
      }
      reference.labels <- object[[]][colnames(object[[sketched.assay]]), label.rd]
      predicted.labels.list <- TransferLablesNN(
        reference.labels = reference.labels,
        weight.matrix = full_sketch.weight)
      object[[paste0(rd.name)]] <- predicted.labels.list$labels
      object[[paste0(rd.name, '.score')]] <- predicted.labels.list$scores
    }
  }
  if (!is.null(reduction.model)) {
    umap.model <- Misc(object = object[[reduction.model]], slot = 'model')
    if (is.null(umap.model)) {
      warning(reduction.model, ' does not have a stored umap model')
      return(object)
    }
    if (verbose) {
      message("Projection to sketch umap")
    }
    if (ncol(full_sketch.nn) > umap.model$n_neighbors) {
      full_sketch.nn@nn.idx <- full_sketch.nn@nn.idx[, 1:umap.model$n_neighbors]
      full_sketch.nn@nn.dist <- full_sketch.nn@nn.dist[, 1:umap.model$n_neighbors]
    }
    proj.umap <- RunUMAP(
      object = full_sketch.nn,
      reduction.model = object[[reduction.model]],
      verbose = verbose,
      assay =  slot(object = object[[reduction]], name = 'assay.used')
    )
    full.umap.reduction <- rev(
      x = make.unique(
        names = c(
          Reductions(object = object),
          paste0('full.',reduction.model)
          )
        )
      )[1]
    Key(object = proj.umap) <- Key(object = full.umap.reduction)
    object[[full.umap.reduction ]] <- proj.umap
  }
  return(object)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Methods for Seurat-defined generics
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @param nsketch A positive integer. The number of sketches to be used in the approximation.
#'                Default is 5000.
#' @param ndims A positive integer or NULL. The number of dimensions to use. If NULL, the number
#'              of dimensions will default to the number of columns in the object.
#' @param method The sketching method to use, defaults to CountSketch.
#' @param eps A numeric. The error tolerance for the approximation in Johnson–Lindenstrauss embeddings,
#'            defaults to 0.5.
#' @param seed A positive integer. The seed for the random number generator, defaults to 123.
#' @param verbose Print progress and diagnostic messages
#' @importFrom Matrix qrR t
#' @importFrom irlba irlba
#'
#' @rdname LeverageScore
#' @method LeverageScore default
#'
#' @export
#' @concept sketching
#'
LeverageScore.default <- function(
  object,
  nsketch = 5000L,
  ndims = NULL,
  method = CountSketch,
  eps = 0.5,
  seed = 123L,
  verbose = TRUE,
  ...
) {
  # Check the dimensions of the object, nsketch, and ndims
  ncells <- ncol(x = object)
  if (ncells < nsketch * 1.5) {
    nv <- ifelse(nrow(x = object) < 50, nrow(x = object) - 1, 50)
    Z <- irlba(A = object, nv = 50, nu = 0, verbose = FALSE)$v
    return(rowSums(x = Z ^ 2))
  }
  if (nrow(x = object) > 5000L) {
    abort(message = "too slow")
  } else if (nrow(x = object) > (ncells / 1.1)) {
    abort(message = "too square")
  }
  ndims <- ndims %||% ncells
  if (nsketch < (1.1 * nrow(x = object))) {
    nsketch <- 1.1 * nrow(x = object)
    warning(
      "'nsketch' is too close to the number of features, setting to ",
      round(x = nsketch, digits = 2L),
      call. = FALSE,
      immediate. = TRUE
    )
  }
  nsketch <- min(nsketch, ndims)
  # Check the method
  if (is_quosure(x = method)) {
    method <- eval(
      expr = quo_get_expr(quo = method),
      envir = quo_get_env(quo = method)
    )
  }
  if (is.character(x = method)) {
    method <- match.fun(FUN = method)
  }
  stopifnot(is.function(x = method))
  # Run the sketching
  if (isTRUE(x = verbose)) {
    message("sampling ", nsketch, " cells for random sketching")
  }
  S <- method(nsketch = nsketch, ncells = ncells, seed = seed, ...)
  object <- t(x = object)
  if (isTRUE(x = verbose)) {
    message("Performing QR decomposition")
  }
  if (inherits(x = object, what = 'IterableMatrix')) {
    temp <- tempdir()
    object.gene_index <- BPCells::transpose_storage_order(matrix = object, tmpdir = temp)
    sa <- as(object = S %*% object.gene_index, Class = 'dgCMatrix')
    rm(object.gene_index)
    unlink(list.files(path = temp, full.names = TRUE))
  } else {
    sa <- S %*% object
  }
  if (!inherits(x = sa, what = 'dgCMatrix')) {
    sa <- as(object = sa, Class = 'dgCMatrix')
  }
  qr.sa <- base::qr(x = sa)
  R <- if (inherits(x = qr.sa, what = 'sparseQR')) {
    qrR(qr = qr.sa)
  } else {
    base::qr.R(qr = qr.sa)
  }
  R.inv <- as.sparse(x = backsolve(r = R, x = diag(x = ncol(x = R))))
  if (isTRUE(x = verbose)) {
    message("Performing random projection")
  }
  JL <- as.sparse(x = JLEmbed(
    nrow = ncol(x = R.inv),
    ncol = ndims,
    eps = eps,
    seed = seed
  ))
  Z <- object %*% (R.inv %*% JL)
  if (inherits(x = Z, what = 'IterableMatrix')) {
    Z.score <- BPCells::matrix_stats(matrix = Z ^ 2, row_stats = 'mean'
                            )$row_stats['mean',]*ncol(x = Z)
    } else {
    Z.score <- rowSums(x = Z ^ 2)
  }
  return(Z.score)
}

#' @param nsketch A positive integer. The number of sketches to be used in the approximation.
#'                Default is 5000.
#' @param ndims A positive integer or NULL. The number of dimensions to use. If NULL, the number
#'              of dimensions will default to the number of columns in the object.
#' @param method The sketching method to use, defaults to CountSketch.
#' @param vf.method VariableFeatures method
#' @param layer layer to use
#' @param eps A numeric. The error tolerance for the approximation in Johnson–Lindenstrauss embeddings,
#'            defaults to 0.5.
#' @param seed A positive integer. The seed for the random number generator, defaults to 123.
#' @param verbose Print progress and diagnostic messages
#'
#' @importFrom SeuratObject EmptyDF
#'
#' @rdname LeverageScore
#' @method LeverageScore StdAssay
#'
#' @export
#' @concept sketching
#'
LeverageScore.StdAssay <- function(
  object,
  nsketch = 5000L,
  ndims = NULL,
  method = CountSketch,
  vf.method = NULL,
  layer = 'data',
  eps = 0.5,
  seed = 123L,
  verbose = TRUE,
  features = NULL,
  ...
) {
  layer <- unique(x = layer) %||% DefaultLayer(object = object)
  layer <- Layers(object = object, search = layer)
  if (!is_quosure(x = method)) {
    method <- enquo(arg = method)
  }
  scores <- EmptyDF(n = ncol(x = object))
  row.names(x = scores) <- colnames(x = object)
  scores[, 1] <- NA_real_
  for (i in seq_along(along.with = layer)) {
    l <- layer[i]
    if (isTRUE(x = verbose)) {
      message("Running LeverageScore for layer ", l)
    }
    features <- features %||% tryCatch({
      VariableFeatures(
        object = object,
        method = vf.method,
        layer = l
      )
    }, error = function(e) {
      stop("Unable to get Variable Features from layer ", l, ". Try providing `features` argument instead.")
    })
    
    scores[Cells(x = object, layer = l), 1] <- LeverageScore(
      object = LayerData(
        object = object,
        layer = l,
        features = features,
        fast = TRUE
      ),
      nsketch = nsketch,
      ndims = ndims %||% ncol(x = object),
      method = method,
      eps = eps,
      seed = seed,
      verbose = verbose,
      ...
    )
  }
  return(scores)
}

#' @rdname LeverageScore
#' @method LeverageScore Assay
#' @export
#'
LeverageScore.Assay <- LeverageScore.StdAssay


#' @param assay assay to use
#' @param nsketch A positive integer. The number of sketches to be used in the approximation.
#'                Default is 5000.
#' @param ndims A positive integer or NULL. The number of dimensions to use. If NULL, the number
#'              of dimensions will default to the number of columns in the object.
#' @param method The sketching method to use, defaults to CountSketch.
#' @param var.name name of slot to store leverage scores
#' @param over.write whether to overwrite slot that currently stores leverage scores. Defaults
#' to FALSE, in which case the 'var.name' is modified if it already exists in the object
#'
#' @rdname LeverageScore
#' @method LeverageScore Seurat
#'
#' @export
#' @concept sketching
#'
LeverageScore.Seurat <- function(
  object,
  assay = NULL,
  nsketch = 5000L,
  ndims = NULL,
  var.name = 'leverage.score',
  over.write = FALSE,
  method = CountSketch,
  vf.method = NULL,
  layer = 'data',
  eps = 0.5,
  seed = 123L,
  verbose = TRUE,
  features = NULL,
  ...
) {
  if (!over.write) {
    var.name <- CheckMetaVarName(object = object, var.name = var.name)
  }
  assay <- assay[1L] %||% DefaultAssay(object = object)
  assay <- match.arg(arg = assay, choices = Assays(object = object))
  method <- enquo(arg = method)
  scores <- LeverageScore(
    object = object[[assay]],
    nsketch = nsketch,
    ndims = ndims,
    method = method,
    vf.method = vf.method,
    layer = layer,
    eps = eps,
    seed = seed,
    verbose = verbose,
    features = features,
    ...
  )
  names(x = scores) <- var.name
  object[[]] <- scores
  return(object)
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Methods for R-defined generics
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Internal
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Generate CountSketch random matrix
#'
#' @inheritParams base::set.seed
#' @param nsketch Number of sketching random cells
#' @param ncells Number of cells in the original data
#' @param ... Ignored
#'
#' @return ...
#'
#' @importFrom Matrix sparseMatrix
#'
#' @export
#' @concept sketching
#'
#' @keywords internal
#'
#' @references Clarkson, KL. & Woodruff, DP.
#' Low-rank approximation and regression in input sparsity time.
#' Journal of the ACM (JACM). 2017 Jan 30;63(6):1-45.
#' \url{https://dl.acm.org/doi/abs/10.1145/3019134};

CountSketch <- function(nsketch, ncells, seed = NA_integer_, ...) {
  if (!is.na(x = seed)) {
    set.seed(seed = seed)
  }
  iv <- xv <- vector(mode = "numeric", length = ncells)
  jv <- seq_len(length.out = ncells)
  for (i in jv) {
    iv[i] <- sample(x = seq_len(length.out = nsketch), size = 1L)
    xv[i] <- sample(x = c(-1L, 1L), size = 1L)
  }
  return(sparseMatrix(
    i = iv,
    j = jv,
    x = xv,
    dims = c(nsketch, ncells)
  ))
}

#' Gaussian sketching
#'
#' @inheritParams CountSketch
#'
#' @return ...
#'
#' @export
#' @concept sketching
#'
#' @keywords internal
#'
GaussianSketch <- function(nsketch, ncells, seed = NA_integer_, ...) {
  if (!is.na(x = seed)) {
    set.seed(seed = seed)
  }
  return(matrix(
    data = rnorm(n = nsketch * ncells, mean = 0L, sd = 1 / (ncells ^ 2)),
    nrow = nsketch,
    ncol = ncells
  ))
}

#' Generate JL random projection embeddings
#'
#' @keywords internal
#'
#' @references Aghila G and Siddharth R (2020).
#' RandPro: Random Projection with Classification. R package version 0.2.2.
#' \url{https://CRAN.R-project.org/package=RandPro}
#'
#' @noRd
#
JLEmbed <- function(nrow, ncol, eps = 0.1, seed = NA_integer_, method = "li") {
  if (!is.na(x = seed)) {
    set.seed(seed = seed)
  }
  method <- method[1L]
  method <- match.arg(arg = method)
  if (!is.null(x = eps)) {
    if (eps > 1 || eps <= 0) {
      stop("'eps' must be 0 < eps <= 1")
    }
    ncol <- floor(x = 4 * log(x = ncol) / ((eps ^ 2) / 2 - (eps ^ 3 / 3)))
  }
  m <- switch(
    EXPR = method,
    "li" = {
      s <- ceiling(x = sqrt(x = ncol))
      prob <- c(
        1 / (2 * s),
        1 - (1 / s),
        1 / (2 * s)
      )
      matrix(
        data = sample(
          x = seq.int(from = -1L, to = 1L),
          size = nrow * ncol,
          replace = TRUE,
          prob = prob
        ),
        nrow = nrow
      )
    }
  )
  return(m)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# S4 Methods
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
