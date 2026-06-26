#' Fast PCA dimensional reduction
#'
#' Faster variant of \code{\link{RunPCA}}. For the common case (a dense scaled
#' feature-by-cell matrix) the default \code{approx = TRUE} computes PCA via the
#' Gram matrix and a Lanczos top-k eigendecomposition (\code{EigenGramPCA}), using
#' Eigen's own BLAS-independent kernels. This replaces \code{irlba} on dense input:
#' it converges to the same truncated SVD, but tighter (tol 1e-10) and faster.
#'
#' Dispatch:
#' \itemize{
#'   \item \code{approx = TRUE} (default), dense input: Gram + Lanczos top-k
#'     eigendecomposition, with \code{prcomp} retained only as a safety net if the
#'     eigen path errors.
#'   \item \code{approx = FALSE}, dense input: exact PCA via \code{prcomp} (as in
#'     base \code{\link{RunPCA}}).
#'   \item Sparse (\code{dgCMatrix}) or on-disk (\code{IterableMatrix}) input:
#'     \code{approx} is forced to \code{TRUE} and the existing \code{irlba} /
#'     BPCells SVD path is used (the dense Gram path does not apply).
#'   \item \code{rev.pca = TRUE}: unchanged from \code{\link{RunPCA}}.
#' }
#' The Gram path costs \eqn{O(\mathrm{nfeatures}^2 \times \mathrm{ncells})} to form
#' the Gram matrix; the eigendecomposition computes only the top \code{npcs} pairs
#' via a Lanczos solver (\code{Spectra}).
#'
#' A second, behavior-neutral optimization avoids computing per-feature variances
#' twice: \code{PrepDR}/\code{PrepDR5} already compute them to drop zero-variance
#' features, so those values are carried through (as a \code{"feature.var"}
#' attribute) and reused for the \code{total.variance} statistic instead of a
#' second \code{RowVar} pass.
#'
#' @param object An object
#' @param ... Arguments passed to other methods
#'
#' @export
#'
RunPCA_fast <- function(object, ...) {
  UseMethod(generic = 'RunPCA_fast', object = object)
}

#' @inheritParams RunPCA
#'
#' @importFrom irlba irlba
#' @importFrom stats prcomp
#' @importFrom utils capture.output
#'
#' @rdname RunPCA_fast
#' @concept dimensional_reduction
#' @export
#'
RunPCA_fast.default <- function(
  object,
  assay = NULL,
  npcs = 50,
  rev.pca = FALSE,
  weight.by.var = TRUE,
  verbose = TRUE,
  ndims.print = 1:5,
  nfeatures.print = 30,
  reduction.key = "PC_",
  seed.use = 42,
  approx = TRUE,
  ...
) {
  if (!is.null(x = seed.use)) {
    set.seed(seed = seed.use)
  }
  # CHANGE: per-feature variances precomputed during data prep (PrepDR/PrepDR5),
  # carried on the matrix so the non-reversed total.variance can reuse them
  # rather than running a second full-matrix RowVar pass. NULL when called
  # directly on a matrix, in which case we fall back to the original behavior.
  feature.var <- attr(x = object, which = 'feature.var')
 if (inherits(x = object, what = 'matrix')) {
   RowVar.function <- RowVar
   svd.function <- irlba
 } else if (inherits(x = object, what = 'dgCMatrix')) {
   RowVar.function <- RowVarSparse
   svd.function <- irlba
 } else if (inherits(x = object, what = 'IterableMatrix')) {
   RowVar.function <- function(x) {
     return(BPCells::matrix_stats(
       matrix = x,
       row_stats = 'variance'
     )$row_stats['variance',])
    }
    svd.function <- function(A, nv, ...) BPCells::svds(A=A, k = nv)
 }
  # CHANGE: sparse (dgCMatrix) and on-disk (IterableMatrix) inputs cannot use the
  # dense Gram path nor an exact prcomp; force approx = TRUE so they take the
  # irlba / BPCells SVD path below.
  if (!inherits(x = object, what = 'matrix')) {
    approx <- TRUE
  }
  if (rev.pca) {
    npcs <- min(npcs, ncol(x = object) - 1)
    pca.results <- svd.function(A = object, nv = npcs, ...)
    total.variance <- sum(RowVar.function(x = t(x = object)))
    sdev <- pca.results$d/sqrt(max(1, nrow(x = object) - 1))
    if (weight.by.var) {
      feature.loadings <- pca.results$u %*% diag(pca.results$d)
    } else{
      feature.loadings <- pca.results$u
    }
    cell.embeddings <- pca.results$v
  }
  else {
    # CHANGE: reuse the precomputed per-feature variances (see top of function)
    # for total.variance instead of a second RowVar pass over the matrix;
    # numerically identical. Original was: sum(RowVar.function(x = object)).
    total.variance <- if (is.null(x = feature.var)) {
      sum(RowVar.function(x = object))
    } else {
      sum(feature.var)
    }
    if (approx && inherits(x = object, what = 'matrix')) {
      # CHANGE: dense approximate PCA (the default) goes through the Gram matrix +
      # Lanczos top-k eigendecomposition (EigenGramPCA), using Eigen's own
      # BLAS-independent kernels, instead of irlba. It converges to the same
      # truncated SVD but tighter (tol 1e-10) and is faster across dataset sizes
      # (irlba paid a per-cell iterative cost; the top-k solver also avoids the
      # full-spectrum eigendecomposition). prcomp is retained purely as a safety
      # net if the eigen path errors.
      npcs <- min(npcs, nrow(x = object))
      gram.results <- tryCatch(
        expr = EigenGramPCA(
          object = object,
          npcs = npcs,
          weight_by_var = weight.by.var
        ),
        error = function(e) NULL
      )
      if (!is.null(x = gram.results)) {
        feature.loadings <- gram.results$loadings
        cell.embeddings <- gram.results$embeddings
        sdev <- as.numeric(x = gram.results$sdev)
      } else {
        pca.results <- prcomp(x = t(object), rank. = npcs, ...)
        feature.loadings <- pca.results$rotation
        sdev <- pca.results$sdev
        if (weight.by.var) {
          cell.embeddings <- pca.results$x
        } else {
          cell.embeddings <- pca.results$x / (pca.results$sdev[1:npcs] * sqrt(x = ncol(x = object) - 1))
        }
      }
    } else if (approx) {
      # Sparse (dgCMatrix) / on-disk (IterableMatrix): unchanged irlba / BPCells
      # SVD path (the dense Gram path does not apply to these inputs).
      npcs <- min(npcs, nrow(x = object) - 1)
      pca.results <- svd.function(A = t(x = object), nv = npcs, ...)
      feature.loadings <- pca.results$v
      sdev <- pca.results$d/sqrt(max(1, ncol(object) - 1))
      if (weight.by.var) {
        cell.embeddings <- pca.results$u %*% diag(pca.results$d)
      } else {
        cell.embeddings <- pca.results$u
      }
    } else {
      npcs <- min(npcs, nrow(x = object))
      pca.results <- prcomp(x = t(object), rank. = npcs, ...)
      feature.loadings <- pca.results$rotation
      sdev <- pca.results$sdev
      if (weight.by.var) {
        cell.embeddings <- pca.results$x
      } else {
        cell.embeddings <- pca.results$x / (pca.results$sdev[1:npcs] * sqrt(x = ncol(x = object) - 1))
      }
    }
  }
  rownames(x = feature.loadings) <- rownames(x = object)
  colnames(x = feature.loadings) <- paste0(reduction.key, 1:npcs)
  rownames(x = cell.embeddings) <- colnames(x = object)
  colnames(x = cell.embeddings) <- colnames(x = feature.loadings)
  reduction.data <- CreateDimReducObject(
    embeddings = cell.embeddings,
    loadings = feature.loadings,
    assay = assay,
    stdev = sdev,
    key = reduction.key,
    misc = list(total.variance = total.variance)
  )
  if (verbose) {
    msg <- capture.output(print(
      x = reduction.data,
      dims = ndims.print,
      nfeatures = nfeatures.print
    ))
    message(paste(msg, collapse = '\n'))
  }
  return(reduction.data)
}

#' @inheritParams RunPCA
#'
#' @rdname RunPCA_fast
#' @concept dimensional_reduction
#' @export
#' @method RunPCA_fast Assay
#'
RunPCA_fast.Assay <- function(
  object,
  assay = NULL,
  features = NULL,
  npcs = 50,
  rev.pca = FALSE,
  weight.by.var = TRUE,
  verbose = TRUE,
  ndims.print = 1:5,
  nfeatures.print = 30,
  reduction.key = "PC_",
  seed.use = 42,
  ...
) {
  data.use <- PrepDR_fast(
    object = object,
    features = features,
    verbose = verbose
  )
  reduction.data <- RunPCA_fast(
    object = data.use,
    assay = assay,
    npcs = npcs,
    rev.pca = rev.pca,
    weight.by.var = weight.by.var,
    verbose = verbose,
    ndims.print = ndims.print,
    nfeatures.print = nfeatures.print,
    reduction.key = reduction.key,
    seed.use = seed.use,
    ...
  )
  return(reduction.data)
}

#' @rdname RunPCA_fast
#' @concept dimensional_reduction
#' @export
#' @method RunPCA_fast StdAssay
#'
RunPCA_fast.StdAssay <- function(
  object,
  assay = NULL,
  features = NULL,
  layer = 'scale.data',
  npcs = 50,
  rev.pca = FALSE,
  weight.by.var = TRUE,
  verbose = TRUE,
  ndims.print = 1:5,
  nfeatures.print = 30,
  reduction.key = "PC_",
  seed.use = 42,
  ...
) {
  data.use <- PrepDR5_fast(
    object = object,
    features = features,
    layer = layer,
    verbose = verbose
  )
  return(RunPCA_fast(
    object = data.use,
    assay = assay,
    npcs = npcs,
    rev.pca = rev.pca,
    weight.by.var = weight.by.var,
    verbose = verbose,
    ndims.print = ndims.print,
    nfeatures.print = nfeatures.print,
    reduction.key = reduction.key,
    seed.use = seed.use,
    ...
  ))
}

#' @inheritParams RunPCA
#'
#' @rdname RunPCA_fast
#' @concept dimensional_reduction
#' @export
#' @method RunPCA_fast Seurat
#'
RunPCA_fast.Seurat <- function(
  object,
  assay = NULL,
  features = NULL,
  npcs = 50,
  rev.pca = FALSE,
  weight.by.var = TRUE,
  verbose = TRUE,
  ndims.print = 1:5,
  nfeatures.print = 30,
  reduction.name = "pca",
  reduction.key = "PC_",
  seed.use = 42,
  ...
) {
  assay <- assay %||% DefaultAssay(object = object)
  reduction.data <- RunPCA_fast(
    object = object[[assay]],
    assay = assay,
    features = features,
    npcs = npcs,
    rev.pca = rev.pca,
    weight.by.var = weight.by.var,
    verbose = verbose,
    ndims.print = ndims.print,
    nfeatures.print = nfeatures.print,
    reduction.key = reduction.key,
    seed.use = seed.use,
    ...
  )
  object[[reduction.name]] <- reduction.data
  object <- LogSeuratCommand(object = object)
  return(object)
}

#' @method RunPCA_fast Seurat5
#' @export
#'
RunPCA_fast.Seurat5 <- function(
  object,
  assay = NULL,
  features = NULL,
  npcs = 50,
  rev.pca = FALSE,
  weight.by.var = TRUE,
  verbose = TRUE,
  ndims.print = 1:5,
  nfeatures.print = 30,
  reduction.name = "pca",
  reduction.key = "PC_",
  seed.use = 42,
  ...
) {
  assay <- assay %||% DefaultAssay(object = object)
  reduction.data <- RunPCA_fast(
    object = object[[assay]],
    assay = assay,
    features = features,
    npcs = npcs,
    rev.pca = rev.pca,
    weight.by.var = weight.by.var,
    verbose = verbose,
    ndims.print = ndims.print,
    nfeatures.print = nfeatures.print,
    reduction.key = reduction.key,
    seed.use = seed.use,
    ...
  )
  object[[reduction.name]] <- reduction.data
  # object <- LogSeuratCommand(object = object)
  return(object)
}

# Copy of PrepDR that additionally carries the per-feature variances it already
# computes (for the kept features, in returned-row order) on the result as a
# "feature.var" attribute, so RunPCA_fast.default can reuse them for
# total.variance instead of recomputing. Variance values and feature selection
# are otherwise identical to PrepDR.
PrepDR_fast <- function(
  object,
  features = NULL,
  slot = 'scale.data',
  verbose = TRUE
) {
  if (length(x = VariableFeatures(object = object)) == 0 && is.null(x = features)) {
    stop("Variable features haven't been set. Run FindVariableFeatures() or provide a vector of feature names.")
  }
  data.use <- GetAssayData(object = object, layer = slot)
  if (nrow(x = data.use ) == 0 && slot == "scale.data") {
    stop("Data has not been scaled. Please run ScaleData and retry")
  }
  features <- features %||% VariableFeatures(object = object)
  features.keep <- unique(x = features[features %in% rownames(x = data.use)])
  if (length(x = features.keep) < length(x = features)) {
    features.exclude <- setdiff(x = features, y = features.keep)
    if (verbose) {
      warning(paste0("The following ", length(x = features.exclude), " features requested have not been scaled (running reduction without them): ", paste0(features.exclude, collapse = ", ")))
    }
  }
  features <- features.keep

  if (inherits(x = data.use, what = 'dgCMatrix')) {
    features.var <- RowVarSparse(mat = data.use[features, ])
  }
  else {
    features.var <- RowVar(x = data.use[features, ])
  }
  names(x = features.var) <- features  # CHANGE: name variances so they can be carried through
  features.keep <- features[features.var > 0]
  if (length(x = features.keep) < length(x = features)) {
    features.exclude <- setdiff(x = features, y = features.keep)
    if (verbose) {
      warning(paste0("The following ", length(x = features.exclude), " features requested have zero variance (running reduction without them): ", paste0(features.exclude, collapse = ", ")))
    }
  }
  features <- features.keep
  features <- features[!is.na(x = features)]
  data.use <- data.use[features, ]
  # CHANGE: carry the kept-feature variances for reuse in RunPCA_fast.default.
  attr(x = data.use, which = 'feature.var') <- features.var[features]
  return(data.use)
}

# Copy of PrepDR5 that additionally carries the per-feature variances it already
# computes (for the kept features, in returned-row order) on the result as a
# "feature.var" attribute, so RunPCA_fast.default can reuse them for
# total.variance. Variance method and feature selection are identical to PrepDR5.
PrepDR5_fast <- function(object, features = NULL, layer = 'scale.data', verbose = TRUE) {
  layer <- layer[1L]
  olayer <- layer
  layer <- Layers(object = object, search = layer)
  if (is.null(layer)) {
    abort(paste0("No layer matching pattern '", olayer, "' not found. Please run ScaleData and retry"))
  }
  data.use <- LayerData(object = object, layer = layer)
  features <- features %||% VariableFeatures(object = object)
  if (!length(x = features)) {
    stop("No variable features, run FindVariableFeatures() or provide a vector of features", call. = FALSE)
  }
  if (is(data.use, "IterableMatrix")) {
    features.var <- BPCells::matrix_stats(matrix=data.use, row_stats="variance")$row_stats["variance",]
  } else {
    features.var <- apply(X = data.use, MARGIN = 1L, FUN = var)
  }
  features.keep <- features[features.var > 0]
  if (!length(x = features.keep)) {
    stop("None of the requested features have any variance", call. = FALSE)
  } else if (length(x = features.keep) < length(x = features)) {
    exclude <- setdiff(x = features, y = features.keep)
    if (isTRUE(x = verbose)) {
      warning(
        "The following ",
        length(x = exclude),
        " features requested have zero variance; running reduction without them: ",
        paste(exclude, collapse = ', '),
        call. = FALSE,
        immediate. = TRUE
      )
    }
  }
  features <- features.keep
  features <- features[!is.na(x = features)]
  features.use <- features[features %in% rownames(data.use)]
  if(!isTRUE(all.equal(features, features.use))) {
    missing_features <- setdiff(features, features.use)
    if(length(missing_features) > 0) {
    warning_message <- paste("The following features were not available: ",
                             paste(missing_features, collapse = ", "),
                             ".", sep = "")
    warning(warning_message, immediate. = TRUE)
    }
  }
  data.use <- data.use[features.use, ]
  # CHANGE: carry the kept-feature variances for reuse in RunPCA_fast.default.
  attr(x = data.use, which = 'feature.var') <- features.var[features.use]
  return(data.use)
}
