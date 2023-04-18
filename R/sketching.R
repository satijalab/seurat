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


#' @importFrom SeuratObject CastAssay Key Key<- Layers
#'
#' @export
#'
#'
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
      verbose = verbose,
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
  cells <- lapply(
    X = seq_along(along.with = layers.data),
    FUN = function(i, seed) {
      set.seed(seed = seed)
      lcells <- Cells(x = object[[assay]], layer = layers.data[i])
      if (length(x = lcells) < ncells) {
        return(lcells)
      }
      return(sample(
        x = lcells,
        size = ncells,
        prob = leverage.score[lcells,]
      ))
    },
    seed = seed
  )
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
  if (!is.null(x = cast)) {
    sketched <- CastAssay(object = sketched, to = cast, ...)
  }
  Key(object = sketched) <- Key(object = sketched.assay, quiet = TRUE)
  object[[sketched.assay]] <- sketched
  DefaultAssay(object = object) <- sketched.assay
  return(object)
}


#' Project full data to the sketch assay
#'
#' @export
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
    proj.emb <- ProjectCellEmbeddings(query = object,
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
  
  object <- TransferSketchLabels(object = object,
                                 atoms = sketched.assay,
                                 reduction = full.reduction,
                                 dims = dims,
                                 k = k.weight,
                                 refdata = refdata,
                                 reduction.model = umap.model,
                                 recompute.neighbors = recompute.neighbors,
                                 recompute.weights = recompute.weights,
                                 verbose = verbose
  )
  return(object)
}


#' Transfer data from sketch data to full data
#' @export
#'
TransferSketchLabels <- function(
  object,
  atoms = 'sketch',
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
    max(Indices(full_sketch.nn)) >  ncol(object[[atoms]]) ||
    !identical(x = full_sketch.nn@alg.info$dims, y =  dims) ||
    !identical(x = full_sketch.nn@alg.info$reduction, y =  reduction) ||
    recompute.neighbors
  
  compute.weights <- is.null(x = full_sketch.weight) ||
    !all(colnames(full_sketch.weight) == Cells(object[[reduction]])) ||
    !all(rownames(full_sketch.weight) == colnames(object[[atoms]]))  ||
    recompute.weights || 
    recompute.neighbors
  
  if (compute.neighbors) {
    if (verbose) {
      message("Finding sketch neighbors")
    }
    full_sketch.nn <- NNHelper(
      query = Embeddings(object[[reduction]])[, dims],
      data = Embeddings(object[[reduction]])[colnames(object[[atoms]]), dims],
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
    full_sketch.weight <- FindWeightsNN(nn.obj = full_sketch.nn,
                                        query.cells = Cells(object[[reduction]]),
                                        reference = colnames(object[[atoms]]),
                                        verbose = verbose)
    rownames(full_sketch.weight) <- colnames(object[[atoms]])
    colnames(full_sketch.weight) <- Cells(object[[reduction]])
  }
  slot(object = object, name = 'tools')$TransferSketchLabels$full_sketch.nn <- full_sketch.nn
  slot(object = object, name = 'tools')$TransferSketchLabels$full_sketch.weight <- full_sketch.weight
  
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
      reference.labels <- object[[]][colnames(object[[atoms]]), label.rd]
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
    Key(proj.umap) <- paste0('ref', Key(proj.umap))
    object[[paste0('ref.',reduction.model )]] <- proj.umap
  }
  return(object)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Methods for Seurat-defined generics
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @importFrom Matrix qrR t
#' @importFrom irlba irlba
#'
#' @method LeverageScore default
#' @export
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
  if (ncells < nsketch*1.5) {
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
    message("sampling ", nsketch, " cells")
  }
  S <- method(nsketch = nsketch, ncells = ncells, seed = seed, ...)
  object <- t(x = object)
  if (isTRUE(x = verbose)) {
    message("Performing QR decomposition")
  }
  if (inherits(x = object, what = 'IterableMatrix')) {
    temp <- tempdir()
    object.gene_index <- transpose_storage_order(matrix = object, tmpdir = temp)
    sa <- as(object = S %*% object, Class = 'dgCMatrix')
    unlink(x = temp, recursive = TRUE)
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
    Z.score <- matrix_stats(matrix = Z ^ 2, row_stats = 'mean'
                            )$row_stats['mean',]*ncol(x = Z)
    } else {
    Z.score <- rowSums(x = Z ^ 2)
  }
  return(Z.score)
}

#' @importFrom Matrix qrR t
#' @method LeverageScore DelayedMatrix
#' @export
#'
LeverageScore.DelayedMatrix <- function(
  object,
  nsketch = 5000L,
  ndims = NULL,
  method = CountSketch,
  eps = 0.5,
  seed = 123L,
  block.size = 1e8,
  verbose = TRUE,
  ...
) {
  check_installed(
    pkg = 'DelayedArray',
    reason = 'for working with delayed matrices'
  )
  if (!is_quosure(x = method)) {
    method <- enquo(arg = method)
  }
  sa <- SketchMatrixProd(object = object,
                         block.size = block.size,
                         nsketch = nsketch,
                         method = method,
                         ...)
  qr.sa <- base::qr(x = sa)
  R <- if (inherits(x = qr.sa, what = 'sparseQR')) {
    qrR(qr = qr.sa)
  } else {
    base::qr.R(qr = qr.sa)
  }
  if (length(x = which(x = diag(x = R) == 0))> 0) {
    warning("not all features are variable features")
    var.index <- which(x = diag(x = R) != 0)
    R <- R[var.index, var.index]
  }
  R.inv <- as.sparse(x = backsolve(r = R, x = diag(x = ncol(x = R))))
  JL <- as.sparse(x = JLEmbed(
    nrow = ncol(x = R.inv),
    ncol = ndims,
    eps = eps,
    seed = seed
  ))
  RP.mat <- R.inv %*% JL
  sparse <- DelayedArray::is_sparse(x = object)
  suppressMessages(setAutoBlockSize(size = block.size))
  cells.grid <- DelayedArray::colAutoGrid(x = object)
  norm.list <- list()
  for (i in seq_len(length.out = length(x = cells.grid))) {
    vp <- cells.grid[[i]]
    block <- DelayedArray::read_block(x = object, viewport = vp, as.sparse = sparse)
    if (sparse) {
      block <- as(object = block, Class = 'dgCMatrix')
    } else {
      block <- as(object = block, Class = 'Matrix')
    }
    norm.list[[i]] <- colSums(x = as.matrix(t(RP.mat) %*% block[1:ncol(R),]) ^ 2)
  }
 scores <- unlist(norm.list)
  return(scores)
}


#' @method LeverageScore StdAssay
#' 
#' @export
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
  ...
) {
  layer <- unique(x = layer) %||% DefaultLayer(object = object)
  layer <- Layers(object = object, search = layer)
  if (!is_quosure(x = method)) {
    method <- enquo(arg = method)
  }
  scores <- SeuratObject:::EmptyDF(n = ncol(x = object))
  row.names(x = scores) <- colnames(x = object)
  scores[, 1] <- NA_real_
  for (i in seq_along(along.with = layer)) {
    l <- layer[i]
    if (isTRUE(x = verbose)) {
      message("Running LeverageScore for layer ", l)
    }
    scores[Cells(x = object, layer = l), 1] <- LeverageScore(
      object = LayerData(
        object = object,
        layer = l,
        features = VariableFeatures(
          object = object,
          method = vf.method,
          layer = l
        ),
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

#' @method LeverageScore Assay
#' @export
#'
LeverageScore.Assay <- LeverageScore.StdAssay

#' @method LeverageScore Seurat
#' @export
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



SketchMatrixProd <- function(
    object,
    method = CountSketch,
    block.size = 1e9,
    nsketch = 5000L,
    seed = 123L,
    ...) {

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
  sparse <- DelayedArray::is_sparse(x = object)
  suppressMessages(setAutoBlockSize(size = block.size))
  cells.grid <- DelayedArray::colAutoGrid(x = object)
  SA.mat <- matrix(data = 0, nrow = nsketch, ncol = nrow(object))
  for (i in seq_len(length.out = length(x = cells.grid))) {
    vp <- cells.grid[[i]]
    block <- DelayedArray::read_block(x = object, viewport = vp, as.sparse = sparse)

    if (sparse) {
      block <- as(object = block, Class = 'dgCMatrix')
    } else {
      block <- as(object = block, Class = 'Matrix')
    }
    ncells.block <- ncol(block)
    S.block <- method(nsketch = nsketch, ncells = ncells.block, seed = seed, ...)
    SA.mat <- SA.mat + as.matrix(S.block %*% t(block))
  }
  return(SA.mat)
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# S4 Methods
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
