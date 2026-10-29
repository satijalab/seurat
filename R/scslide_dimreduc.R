#' Run Diffusion Map
#'
#' Run diffusion map dimensionality reduction on single-cell data.
#'
#' @param object An object to run diffusion map on
#' @param assay Name of Assay diffusion map is being run on
#' @param ncomp Number of components to compute
#' @param metric This determines the choice of metric used to measure
#' distance in the input space. Check DiffusionMap() for details. If a distance.matrix is provided,
#' this will be ignored.
#' @param distance.matrix If set, runs DiffusionMap on the given distance matrix
#' instead of data matrix.
#' @param reduction.key dimensional reduction key, specifies the string before
#' the number for the dimension names.
#' @param reduction.name the name of the DimReduc object
#' @param weight.by.var Weight the cell embeddings by the eigenvalues of the diffusion map
#' @param verbose To print diagnostic messages. Default is TRUE.
#' @param seed.use Set a random seed.  Setting NULL will not set a seed.
#' @param features Features to compute Diffusion Map
#' @param layer The layer in \code{assay} to use when running diffusion map analysis.
#' @param ... Additional arguments to be passed to DiffusionMap()
#'
#' @return Returns a DimReduc object with diffusion map results
#' @export
#' @concept dimensional_reduction
#'
#' @examples
#' \dontrun{
#' # Run diffusion map on a Seurat object
#' seurat_obj <- RunDiffusionMap(seurat_obj, ncomp = 10)
#' }
#'
#' @importFrom SeuratObject CreateDimReducObject DefaultAssay Assays
#'
RunDiffusionMap <- function(object, ...) {
  UseMethod(generic = 'RunDiffusionMap', object = object)
}

#' @rdname RunDiffusionMap
#' @method RunDiffusionMap default
#' @export
RunDiffusionMap.default <- function(
    object,
    assay = NULL,
    metric = "cosine",
    distance.matrix = NULL,
    ncomp = 20,
    reduction.name = "DiffMap",
    reduction.key = "DC_",
    weight.by.var = TRUE,
    verbose = TRUE,
    seed.use = 101,
    ...
) {
  # Check if destiny is available
  if (!requireNamespace("destiny", quietly = TRUE)) {
    stop(
      "The 'destiny' package is required for RunDiffusionMap() but is not installed.\n\n",
      "To install destiny, please run:\n",
      "  if (!require('BiocManager', quietly = TRUE))\n",
      "    install.packages('BiocManager')\n",
      "  BiocManager::install('destiny')\n\n",
      call. = FALSE
    )
  }

  # Check for knn.covertree with a warning instead of error
  if (!requireNamespace("knn.covertree", quietly = TRUE)) {
    warning(
      "Package 'knn.covertree' is not installed. ",
      "Some advanced distance metrics may not be available.\n",
      "To install: install.packages('knn.covertree')",
      call. = FALSE
    )
  }
  #
  if (!is.null(x = seed.use)) {
    set.seed(seed = seed.use)
  }
  ncomp <- min(ncomp, ncol(x = object))

  # run Diffusion Map
  if (is.null(distance.matrix) & !is.null(metric)){
    dm_results <- destiny::DiffusionMap(data = t(object), n_eigs = ncomp, distance = metric)
  } else if (!is.null(distance.matrix)){
    dm_results <- destiny::DiffusionMap(distance = distance.matrix, n_eigs = ncomp)
  } else {
    stop("You must provide either a 'metric' or a 'distance.matrix' to DiffusionMap().")
  }
  #
  if(isTRUE(x = weight.by.var)){
    stdev <- dm_results@eigenvalues
    stdev[which(stdev < 0)] <- 0
    diag_mat <- diag(sqrt(abs(stdev)))
    cell.embeddings <- dm_results@eigenvectors %*% diag_mat
  } else {
    cell.embeddings <- dm_results@eigenvectors
  }
  #
  rownames(cell.embeddings) <- colnames(object)
  reduction.data <- CreateDimReducObject(
    embeddings = cell.embeddings,
    assay = assay,
    stdev = dm_results@eigenvalues,
    key = reduction.key
  )
  #
  return(reduction.data)
}

#' @rdname RunDiffusionMap
#' @method RunDiffusionMap Assay
#' @export
RunDiffusionMap.Assay <- function(
    object,
    assay = NULL,
    features = NULL,
    ncomp = 20,
    metric = "cosine",
    distance.matrix = NULL,
    reduction.name = "DiffMap",
    reduction.key = "DC_",
    weight.by.var = TRUE,
    verbose = TRUE,
    seed.use = 101,
    ...
) {
  # Get internal function from Seurat
  
  data.use <- PrepDR(
    object = object,
    features = features,
    verbose = verbose
  )
  reduction.data <- RunDiffusionMap(
    object = data.use,
    assay = assay,
    ncomp = ncomp,
    metric = metric,
    distance.matrix = distance.matrix,
    reduction.name = reduction.name,
    reduction.key = reduction.key,
    weight.by.var = weight.by.var,
    verbose = verbose,
    seed.use = seed.use,
    ...
  )
  return(reduction.data)
}

#' @rdname RunDiffusionMap
#' @method RunDiffusionMap StdAssay
#' @export
RunDiffusionMap.StdAssay <- function(
    object,
    assay = NULL,
    features = NULL,
    layer = 'scale.data',
    ncomp = 20,
    metric = "cosine",
    distance.matrix = NULL,
    reduction.name = "DiffMap",
    reduction.key = "DC_",
    weight.by.var = TRUE,
    verbose = TRUE,
    seed.use = 101,
    ...
) {
  # Get internal function from Seurat
  
  if (is.null(distance.matrix)) {
    data.use <- PrepDR5(
      object = object,
      features = features,
      layer = layer,
      verbose = verbose
    )
  } else {
    data.use <- PrepDR5(
      object = object,
      features = rownames(object),
      layer = layer,
      verbose = verbose
    )
  }
  reduction.data <- RunDiffusionMap(
    object = data.use,
    assay = assay,
    ncomp = ncomp,
    metric = metric,
    distance.matrix = distance.matrix,
    reduction.name = reduction.name,
    reduction.key = reduction.key,
    weight.by.var = weight.by.var,
    verbose = verbose,
    seed.use = seed.use,
    ...
  )
  return(reduction.data)
}

#' @rdname RunDiffusionMap
#' @method RunDiffusionMap Seurat
#' @export
RunDiffusionMap.Seurat <- function(
    object,
    assay = NULL,
    features = NULL,
    layer = 'scale.data',
    ncomp = 20,
    metric = "cosine",
    distance.matrix = NULL,
    reduction.name = "DiffMap",
    reduction.key = "DC_",
    weight.by.var = TRUE,
    verbose = TRUE,
    seed.use = 101,
    ...
) {
  assay <- assay %||% DefaultAssay(object = object)
  assay <- match.arg(arg = assay, choices = Assays(object = object))

  reduction.data <- RunDiffusionMap(
    object = object[[assay]],
    assay = assay,
    features = features,
    layer = layer,
    ncomp = ncomp,
    metric = metric,
    distance.matrix = distance.matrix,
    reduction.name = reduction.name,
    reduction.key = reduction.key,
    weight.by.var = weight.by.var,
    verbose = verbose,
    seed.use = seed.use,
    ...
  )
  object[[reduction.name]] <- reduction.data
  object <- LogSeuratCommand(object = object)
  return(object)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Partial Least Squares (PLS) dimensional reduction
# (moved from scslide_pls.R; PLS is a dimensional reduction method)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Resolve Y (or Y.add) into a numeric matrix.
#
# Accepts a numeric matrix (returned as-is), a bare numeric vector
# (promoted to single-column matrix), or a data.frame (converted via
# model.matrix(~. + 0, data = ...)).  Character vectors are NOT
# handled here; the Seurat method resolves those before calling this helper.
#
# @param Y Input response object.
# @param expected_nrow Expected number of rows (for validation); NULL skips.
# @param label Label used in error messages (e.g. "Y", "Y.add").
# @return A numeric matrix, or NULL if Y is NULL.
# @keywords internal
.resolve_Y_matrix <- function(Y, expected_nrow = NULL, label = "Y") {
  if (is.null(Y)) return(NULL)
  if (is.matrix(Y) && is.numeric(Y)) {
    Y_mat <- Y
  } else if (is.numeric(Y) && is.null(dim(Y))) {
    # bare numeric vector -> single-column matrix
    Y_mat <- matrix(Y, ncol = 1)
  } else if (is.data.frame(Y)) {
    Y_mat <- model.matrix(~. + 0, data = Y)
  } else {
    stop(sprintf(
      "'%s' must be a numeric matrix, a data.frame, or a character vector of metadata column names. Got: %s",
      label, paste(class(Y), collapse = ", ")
    ))
  }
  if (!is.null(expected_nrow) && nrow(Y_mat) != expected_nrow) {
    stop(sprintf(
      "nrow(%s) is %d but expected %d (number of cells/samples).",
      label, nrow(Y_mat), expected_nrow
    ))
  }
  Y_mat
}

# Compute R2 and RMSEP from raw Y and fitted values (no pls model object needed).
#
# @param Y n x q response matrix.
# @param fitted_values n x q x ncomp array of cumulative fitted values.
# @return A list with \code{$R2} and \code{$RMSEP}, each an \code{mvrVal}-class
#   object matching the structure of \code{pls::R2()} / \code{pls::RMSEP()}.
# @keywords internal
.compute_r2_rmsep <- function(Y, fitted_values) {
  Y <- as.matrix(Y)
  n <- nrow(Y)
  q <- ncol(Y)
  ncomp <- dim(fitted_values)[3]

  # SST per response (intercept-only model = column means)
  Y_mean <- matrix(colMeans(Y), nrow = n, ncol = q, byrow = TRUE)
  SST <- colSums((Y - Y_mean)^2)  # length q

  # R2$val layout from pls: 1 x (ncomp+1) x q
  # [1,1,] = intercept (0 comps), [1,a+1,] = a comps
  R2_val <- array(NA_real_, dim = c(1, ncomp + 1, q))
  RMSEP_val <- array(NA_real_, dim = c(1, ncomp + 1, q))

  # Intercept-only model
  R2_val[1, 1, ] <- 0
  RMSEP_val[1, 1, ] <- sqrt(SST / n)

  for (a in seq_len(ncomp)) {
    fitted_a <- fitted_values[, , a, drop = FALSE]
    dim(fitted_a) <- c(n, q)
    SSE <- colSums((Y - fitted_a)^2)
    R2_val[1, a + 1, ] <- 1 - SSE / SST
    RMSEP_val[1, a + 1, ] <- sqrt(SSE / n)
  }

  # Set dimnames matching pls convention
  comp_names <- c("(Intercept)", paste(seq_len(ncomp), "comps"))
  resp_names <- colnames(Y)
  if (is.null(resp_names)) resp_names <- paste0("Y", seq_len(q))

  dimnames(R2_val) <- list("train", comp_names, resp_names)
  dimnames(RMSEP_val) <- list("train", comp_names, resp_names)

  R2_out <- list(val = R2_val, type = "train")
  class(R2_out) <- "mvrVal"

  RMSEP_out <- list(val = RMSEP_val, type = "train")
  class(RMSEP_out) <- "mvrVal"

  list(R2 = R2_out, RMSEP = RMSEP_out)
}

#' Run Partial Least Squares (PLS) on Seurat Objects
#'
#' Performs Partial Least Squares regression analysis on single-cell data.
#' This function provides methods for plsr, spls, and cppls from the pls and spls packages.
#'
#' @param object An object to run PLS on
#' @param assay Name of Assay PLS is being run on
#' @param ncomp Number of components to compute
#' @param Y The response variable(s) for PLS regression.
#'   Accepted formats:
#'   \itemize{
#'     \item A character vector of column names in the Seurat object metadata
#'           (Seurat method only).
#'     \item A \code{data.frame}, which is converted via
#'           \code{model.matrix(~. + 0, data = Y)}.
#'     \item A numeric matrix (used as-is).
#'     \item A numeric vector (promoted to a single-column matrix).
#'   }
#'   The number of rows must match the number of cells/samples.
#' @param Y.add Additional response variable(s) containing relevant
#'   information about the observations. Only used for cppls. Accepts the same
#'   formats as \code{Y}.
#' @param pls.function PLS function from pls package to run (options: plsr, spls, cppls)
#' @param verbose Print progress messages, diagnostics (e.g. R2), and the top
#' genes associated with high/low loadings for the components
#' @param ndims.print components to print genes for
#' @param nfeatures.print Number of genes to print for each component
#' @param reduction.name the name of the DimReduc object
#' @param reduction.key dimensional reduction key, specifies the string before
#' the number for the dimension names.
#' @param seed.use Set a random seed.  Setting NULL will not set a seed.
#' @param eta Thresholding parameter that controls the sparsity of the spls method (larger --> sparser). eta should be between 0 and 1.
#' @param features Features to compute PLS on
#' @param save.model Logical; if TRUE, save model components (coefficients, Xmeans, Ymeans, feature.mean, feature.sd) into the misc slot for later prediction. The per-gene \code{feature.mean} and \code{feature.sd} are computed from the training \code{data} layer so that \code{PredictPLS} can scale new data into the same coordinate system. Default FALSE.
#' @param threads Integer; number of threads for BPCells parallel computation
#'   (IterableMatrix path only). \code{NULL} (default) uses the package-wide
#'   thread count set by \code{\link{setThreads}} (\code{1} unless changed);
#'   \code{0} = auto-detect via \code{parallel::detectCores(logical = FALSE)};
#'   values \eqn{\ge 1} use that many threads. Ignored for in-memory data.
#' @param layer The layer in \code{assay} to use when running PLS analysis.
#' @param materialize Logical; for on-disk (BPCells) assays only. If \code{TRUE},
#'   the feature-subsetted matrix is written once to a single contiguous BPCells
#'   store before PLS, instead of re-streaming a fragmented lazy chain on each of
#'   the ~\code{2*ncomp+3} passes. This can substantially speed up PLS on a
#'   CellAnova-corrected assay whose on-disk root is a per-sample
#'   \code{ColBindMatrices}. Ignored for in-memory data. Default \code{FALSE}.
#' @param materialize.dir Directory for the materialized matrix (only used when
#'   \code{materialize = TRUE}). \code{NULL} (default) writes to a temporary
#'   directory that is removed after PLS completes; supplying a path keeps the
#'   matrix there for reuse by downstream steps.
#' @param overwrite Logical; only used when \code{materialize = TRUE} and a
#'   \code{materialize.dir} is supplied. Passed to
#'   \code{BPCells::write_matrix_dir()} to control whether an existing directory
#'   may be overwritten. Default \code{FALSE} so that a non-empty
#'   \code{materialize.dir} is never clobbered without explicit consent.
#' @param ... Additional arguments to be passed to the PLS function
#'
#' @return Returns a DimReduc object with PLS results
#' @export
#' @concept dimensional_reduction
#'
#' @examples
#' \dontrun{
#' # Run PLS on a Seurat object using metadata column names
#' seurat_obj <- RunPLS(seurat_obj, Y = "condition", ncomp = 10)
#'
#' # Run PLS with a pre-built numeric design matrix
#' design_mat <- model.matrix(~ age + sex, data = seurat_obj[[]])
#' seurat_obj <- RunPLS(seurat_obj, Y = design_mat, ncomp = 10)
#' }
#'
#' @import pls
#' @importFrom stats model.matrix sd
#' @importFrom SeuratObject CreateDimReducObject DefaultAssay Assays Cells
#' @importFrom utils getFromNamespace
#' @rdname RunPLS
#'
RunPLS <- function(object, ...) {
  UseMethod(generic = 'RunPLS', object = object)
}

#' @rdname RunPLS
#' @method RunPLS default
#' @export
RunPLS.default <- function(
    object,
    assay = NULL,
    ncomp = 20,
    Y = NULL,
    Y.add = NULL,
    pls.function = c("plsr", "spls", "cppls"),
    verbose = TRUE,
    ndims.print = 1:5,
    nfeatures.print = 30,
    reduction.name = "pls",
    reduction.key = "PLS_",
    seed.use = 42,
    eta = 0.5,
    save.model = FALSE,
    ...
) {
  pls.function <- match.arg(arg = pls.function)
  if (!is.null(x = seed.use)) {
    set.seed(seed = seed.use)
  }
  ncomp <- min(ncomp, ncol(x = object))
  pls.fxn <- if (pls.function == "spls") {
    if (!requireNamespace("spls", quietly = TRUE)) {
      stop("The 'spls' package is required for pls.function = 'spls' but is not installed. Please install it with install.packages('spls').")
    }
    spls::spls
  } else {
    eval(expr = parse(text = pls.function))
  }

  Y_mat <- .resolve_Y_matrix(Y, expected_nrow = ncol(object), label = "Y")
  # run PLS using the selected method
  if (pls.function == "cppls" & !is.null(Y.add)){
    requireNamespace("pls", quietly = FALSE)
    message("Fitting Y.add in CPPLS...")
    Y.add <- .resolve_Y_matrix(Y.add, expected_nrow = ncol(object), label = "Y.add")
    pls.results <- pls.fxn(Y_mat ~ t(object), ncomp = ncomp, validation = "none", Y.add = Y.add, ...)
  } else if (pls.function == "spls") {
    pls.results <- pls.fxn(x = t(object), y = Y_mat, K = ncomp, fit = "kernelpls", eta = eta, ...)
  } else {
    requireNamespace("pls", quietly = FALSE)
    pls.results <- pls.fxn(Y_mat ~ t(object), ncomp = ncomp, validation = "none", ...)
  }
  #
  if (pls.function == "spls"){
    feature.loadings <- pls.results$projection
    colnames(feature.loadings) <- paste0(reduction.key, 1:ncol(feature.loadings))
    rownames(feature.loadings) <- rownames(object)[pls.results$A]
    cell.embeddings <- t(as.matrix(object[pls.results$A, ])) %*% feature.loadings
    colnames(cell.embeddings) <- paste0(reduction.key, 1:ncol(cell.embeddings))
    stdev <- numeric()
    misc <- list()
    if (save.model) {
      warning("save.model is not supported for spls; model components will not be saved.")
    }
  } else {
    feature.loadings <- unclass(pls.results$projection)
    colnames(feature.loadings) <- paste0(reduction.key, 1:ncol(feature.loadings))
    rownames(feature.loadings) <- rownames(object)
    cell.embeddings <- unclass(pls.results$scores)
    colnames(cell.embeddings) <- paste0(reduction.key, 1:ncol(cell.embeddings))
    # if we are using plsr or cppls, we are able to retrieve the var_explained for each comp
    stdev <- pls.results$Xvar
    # other miscellaneous data might be useful for model selection
    misc <- list(R2 = pls::R2(pls.results),
                 RMSEP = pls::RMSEP(pls.results))
    # a diagnosis component
    R2 <- matrix(pls::R2(pls.results)$val, byrow = T, ncol = dim(pls::R2(pls.results)$val)[2])
    mean_final_R2 <- mean(R2[nrow(R2), ])
    message("The average R2 of the PLS model is ", mean_final_R2)
    if (save.model) {
      misc$model <- list(
        coefficients = pls.results$coefficients,
        Xmeans = pls.results$Xmeans,
        Ymeans = pls.results$Ymeans
      )
    }
  }
  #
  reduction.data <- CreateDimReducObject(
    embeddings = cell.embeddings,
    loadings = feature.loadings,
    stdev = stdev,
    assay = assay,
    key = reduction.key,
    misc = misc
  )
  if (verbose) {
    print(x = reduction.data, dims = ndims.print, nfeatures = nfeatures.print)
  }
  return(reduction.data)
}

#' @rdname RunPLS
#' @method RunPLS IterableMatrix
#' @export
RunPLS.IterableMatrix <- function(
    object,
    assay = NULL,
    ncomp = 20,
    Y = NULL,
    Y.add = NULL,
    pls.function = c("plsr", "spls", "cppls"),
    verbose = TRUE,
    ndims.print = 1:5,
    nfeatures.print = 30,
    reduction.name = "pls",
    reduction.key = "PLS_",
    seed.use = 42,
    eta = 0.5,
    save.model = FALSE,
    threads = NULL,
    ...
) {
  pls.function <- match.arg(arg = pls.function)

  if (!is.null(x = seed.use)) {
    set.seed(seed = seed.use)
  }

  n.samples <- ncol(object)
  n.features <- nrow(object)
  max.ncomp <- min(n.samples, n.features) - 1L
  if (max.ncomp < 1L) {
    rlang::abort("IterableMatrix must have at least 2 samples/features to run PLS.")
  }
  ncomp <- min(ncomp, max.ncomp)

  Y_mat <- .resolve_Y_matrix(Y, expected_nrow = n.samples, label = "Y")

  if (pls.function == "cppls") {
    pls.results <- cppls_ondisk(X = object, Y = Y_mat, Y.add = Y.add,
                                ncomp = ncomp, threads = threads, ...)
  } else if (pls.function == "plsr") {
    if (!is.null(Y.add)) {
      warning("Y.add is not supported for plsr/kernelpls and will be ignored")
    }
    pls.results <- kernelpls(X = object, Y = Y_mat, ncomp = ncomp,
                             threads = threads, ...)
  } else {
    rlang::abort("Only pls.function = 'plsr' or 'cppls' is supported for IterableMatrix currently")
  }

  feature.loadings <- unclass(pls.results$projection)
  colnames(feature.loadings) <- paste0(reduction.key, seq_len(ncol(feature.loadings)))
  feature.names <- rownames(object)
  if (!is.null(feature.names) && length(feature.names) == nrow(feature.loadings)) {
    rownames(feature.loadings) <- feature.names
  }

  cell.embeddings <- unclass(pls.results$scores)
  colnames(cell.embeddings) <- paste0(reduction.key, seq_len(ncol(cell.embeddings)))
  sample.names <- colnames(object)
  if (!is.null(sample.names) && length(sample.names) == nrow(cell.embeddings)) {
    rownames(cell.embeddings) <- sample.names
  }

  stdev <- pls.results$Xvar
  r2_rmsep <- .compute_r2_rmsep(Y_mat, pls.results$fitted.values)
  misc <- list(R2 = r2_rmsep$R2,
               RMSEP = r2_rmsep$RMSEP)
  # diagnostic message matching RunPLS.default()
  R2 <- matrix(r2_rmsep$R2$val, byrow = TRUE, ncol = dim(r2_rmsep$R2$val)[2])
  mean_final_R2 <- mean(R2[nrow(R2), ])
  message("The average R2 of the PLS model is ", mean_final_R2)
  if (save.model) {
    misc$model <- list(
      coefficients = pls.results$coefficients,
      Xmeans = pls.results$Xmeans,
      Ymeans = pls.results$Ymeans
    )
  }

  reduction.data <- CreateDimReducObject(
    embeddings = cell.embeddings,
    loadings = feature.loadings,
    stdev = stdev,
    assay = assay,
    key = reduction.key,
    misc = misc
  )
  if (verbose) {
    print(x = reduction.data, dims = ndims.print, nfeatures = nfeatures.print)
  }
  return(reduction.data)
}

#' @rdname RunPLS
#' @method RunPLS Assay
#' @export
RunPLS.Assay <- function(
    object,
    assay = NULL,
    features = NULL,
    ncomp = 20,
    Y = NULL,
    Y.add = NULL,
    pls.function = c("plsr", "spls", "cppls"),
    verbose = TRUE,
    ndims.print = 1:5,
    nfeatures.print = 30,
    reduction.name = "pls",
    reduction.key = "PLS_",
    seed.use = 42,
    eta = 0.5,
    save.model = FALSE,
    threads = NULL,
    ...
) {
  # Get internal function from Seurat
  
  data.use <- PrepDR(
    object = object,
    features = features,
    verbose = verbose
  )
  reduction.data <- RunPLS(
    object = data.use,
    assay = assay,
    ncomp = ncomp,
    Y = Y,
    Y.add = Y.add,
    pls.function = pls.function,
    verbose = verbose,
    ndims.print = ndims.print,
    nfeatures.print = nfeatures.print,
    reduction.key = reduction.key,
    seed.use = seed.use,
    eta = eta,
    save.model = save.model,
    threads = threads,
    ...

  )
  return(reduction.data)
}

#' @rdname RunPLS
#' @method RunPLS StdAssay
#' @export
RunPLS.StdAssay <- function(
    object,
    assay = NULL,
    features = NULL,
    layer = 'scale.data',
    ncomp = 20,
    Y = NULL,
    Y.add = NULL,
    pls.function = c("plsr", "spls", "cppls"),
    verbose = TRUE,
    ndims.print = 1:5,
    nfeatures.print = 30,
    reduction.name = "pls",
    reduction.key = "PLS_",
    seed.use = 42,
    eta = 0.5,
    save.model = FALSE,
    threads = NULL,
    materialize = FALSE,
    materialize.dir = NULL,
    overwrite = FALSE,
    ...
) {
  # Get internal function from Seurat

  data.use <- PrepDR5(
    object = object,
    features = features,
    layer = layer,
    verbose = verbose
  )

  # Optional: checkpoint the (feature-subsetted) on-disk matrix to a single
  # contiguous BPCells store before streaming. RunPLS streams the matrix
  # ~2*ncomp+3 times; when the source is a fragmented lazy chain (e.g. a
  # per-sample ColBindMatrices root, as produced by CellAnova correction),
  # materializing once is far cheaper than re-streaming the fragments each pass.
  # materialize.dir = NULL writes to a temp dir removed after RunPLS; supplying a
  # path keeps the matrix for reuse by downstream steps.
  if (isTRUE(materialize) && inherits(data.use, "IterableMatrix")) {
    if (!requireNamespace("BPCells", quietly = TRUE)) {
      rlang::abort("'materialize = TRUE' requires the BPCells package.")
    }
    dir <- materialize.dir %||% tempfile(pattern = "scSLIDE_materialize_")
    if (verbose) {
      message("Materializing on-disk layer to a contiguous store: ", dir)
    }
    data.use <- BPCells::write_matrix_dir(mat = data.use, dir = dir, overwrite = overwrite)
    if (is.null(materialize.dir)) {
      on.exit(unlink(dir, recursive = TRUE), add = TRUE)
    }
  }

  reduction.data <- RunPLS(
    object = data.use,
    assay = assay,
    ncomp = ncomp,
    Y = Y,
    Y.add = Y.add,
    pls.function = pls.function,
    verbose = verbose,
    ndims.print = ndims.print,
    nfeatures.print = nfeatures.print,
    reduction.key = reduction.key,
    seed.use = seed.use,
    eta = eta,
    save.model = save.model,
    threads = threads,
    ...

  )
  return(reduction.data)
}

#' @rdname RunPLS
#' @method RunPLS Seurat
#' @export
RunPLS.Seurat <- function(
    object,
    assay = NULL,
    features = NULL,
    ncomp = 20,
    Y = NULL,
    Y.add = NULL,
    pls.function = c("plsr", "spls", "cppls"),
    verbose = TRUE,
    ndims.print = 1:5,
    nfeatures.print = 30,
    reduction.name = "pls",
    reduction.key = "PLS_",
    seed.use = 42,
    eta = 0.5,
    save.model = FALSE,
    threads = NULL,
    materialize = FALSE,
    materialize.dir = NULL,
    overwrite = FALSE,
    ...
) {
  assay <- assay %||% DefaultAssay(object = object)
  assay <- match.arg(arg = assay, choices = Assays(object = object))
  pls.function <- match.arg(arg = pls.function)

  # Resolve Y: character -> metadata lookup; matrix/data.frame -> use directly
  cells <- Cells(object[[assay]])
  if (is.character(Y)) {
    Y <- object[[Y]][cells, , drop = FALSE]
  } else if (is.matrix(Y) || is.data.frame(Y) || (is.numeric(Y) && is.null(dim(Y)))) {
    Y <- .resolve_Y_matrix(Y, expected_nrow = length(cells), label = "Y")
    # If rownames exist, reorder to match assay cell order
    if (!is.null(rownames(Y))) {
      missing <- setdiff(cells, rownames(Y))
      if (length(missing) > 0) {
        stop(sprintf("Y is missing %d cells present in the assay.", length(missing)))
      }
      Y <- Y[cells, , drop = FALSE]
    }
  } else {
    stop("'Y' must be a character vector of metadata column names, a numeric matrix/vector, or a data.frame.")
  }

  # Resolve Y.add with the same logic
  if (!is.null(Y.add)) {
    if (is.character(Y.add)) {
      Y.add <- object[[Y.add]][cells, , drop = FALSE]
    } else if (is.matrix(Y.add) || is.data.frame(Y.add) || (is.numeric(Y.add) && is.null(dim(Y.add)))) {
      Y.add <- .resolve_Y_matrix(Y.add, expected_nrow = length(cells), label = "Y.add")
      if (!is.null(rownames(Y.add))) {
        missing <- setdiff(cells, rownames(Y.add))
        if (length(missing) > 0) {
          stop(sprintf("Y.add is missing %d cells present in the assay.", length(missing)))
        }
        Y.add <- Y.add[cells, , drop = FALSE]
      }
    } else {
      stop("'Y.add' must be a character vector of metadata column names, a numeric matrix/vector, or a data.frame.")
    }
  }

  reduction.data <- RunPLS(
    object = object[[assay]],
    assay = assay,
    features = features,
    ncomp = ncomp,
    Y = Y,
    Y.add = Y.add,
    pls.function = pls.function,
    verbose = verbose,
    ndims.print = ndims.print,
    nfeatures.print = nfeatures.print,
    reduction.key = reduction.key,
    seed.use = seed.use,
    eta = eta,
    save.model = save.model,
    threads = threads,
    materialize = materialize,
    materialize.dir = materialize.dir,
    overwrite = overwrite,
    ...
  )

  # Store per-gene scaling parameters from the training data layer so that
  # PredictPLS.Seurat can transform new (unscaled) data into training-scaled
  # space before applying PLS coefficients.
  if (save.model && !is.null(reduction.data@misc$model)) {
    pls.features <- rownames(Loadings(reduction.data))
    train.data <- LayerData(object[[assay]], layer = "data")
    train.data <- train.data[pls.features, , drop = FALSE]

    if (inherits(train.data, "IterableMatrix")) {
      stats <- BPCells::matrix_stats(train.data, row_stats = "variance")
      feat.mean <- stats$row_stats["mean", ]
      feat.sd <- sqrt(stats$row_stats["variance", ])
    } else {
      feat.mean <- rowMeans(train.data)
      feat.sd <- apply(train.data, 1, sd)
    }

    # Edge cases: NA means -> 0, zero/NA SDs -> 1
    feat.mean[is.na(feat.mean)] <- 0
    feat.sd[is.na(feat.sd) | feat.sd == 0] <- 1

    names(feat.mean) <- pls.features
    names(feat.sd) <- pls.features
    reduction.data@misc$model$feature.mean <- feat.mean
    reduction.data@misc$model$feature.sd <- feat.sd
  }

  object[[reduction.name]] <- reduction.data
  object <- LogSeuratCommand(object = object)
  return(object)
}


#' Kernel Partial Least Squares
#'
#' Generic function for kernel-PLS. The \code{IterableMatrix} method computes
#' kernel-PLS for BPCells on-disk matrices.
#'
#' @param X Predictor matrix (or an IterableMatrix).
#' @param ... Arguments passed to methods.
#' @export
kernelpls <- function(X, ...) UseMethod("kernelpls")

#' Kernel Partial Least Squares for IterableMatrix
#'
#' @description
#' Compute kernel-PLS decomposition for large on-disk BPCells matrices using
#' the algorithm of Dayal and MacGregor (1997).
#'
#' @param X An IterableMatrix (features x samples, the BPCells convention).
#'   The transpose is handled internally.
#' @param Y A matrix or vector of responses (samples x q).
#' @param ncomp Integer; number of PLS components to compute.
#' @param center Logical; whether to center X and Y (default \code{TRUE}).
#' @param stripped Logical; if \code{TRUE} return only coefficients, Xmeans, and
#'   Ymeans for speed (default \code{FALSE}).
#' @param threads Integer; number of threads for BPCells parallel computation.
#'   \code{NULL} (default) uses the package-wide thread count set by
#'   \code{\link{setThreads}} (\code{1} unless changed);
#'   \code{0} = auto-detect via \code{parallel::detectCores(logical = FALSE)};
#'   values \eqn{\ge 1} use that many threads.
#' @param ... Additional arguments (currently unused).
#'
#' @return A list with components:
#' \describe{
#'   \item{coefficients}{Regression coefficients (p x q x ncomp array).}
#'   \item{scores}{X-scores (n x ncomp matrix).}
#'   \item{loadings}{X-loadings (p x ncomp matrix).}
#'   \item{loading.weights}{Loading weights (p x ncomp matrix).}
#'   \item{Yscores}{Y-scores (n x ncomp matrix).}
#'   \item{Yloadings}{Y-loadings (q x ncomp matrix).}
#'   \item{projection}{Projection matrix (p x ncomp matrix).}
#'   \item{Xmeans}{Column means of X (length p).}
#'   \item{Ymeans}{Column means of Y (length q).}
#'   \item{fitted.values}{Fitted Y values (n x q x ncomp array).}
#'   \item{residuals}{Residuals (n x q x ncomp array).}
#'   \item{Xvar}{Variance explained in X per component (length ncomp).}
#'   \item{Xtotvar}{Total variance in X (scalar).}
#' }
#' If \code{stripped = TRUE}, only \code{coefficients}, \code{Xmeans}, and
#' \code{Ymeans} are returned.
#'
#' @details
#' Implements the kernel-PLS algorithm (Dayal & MacGregor, 1997) for BPCells
#' on-disk matrices. The large predictor matrix X is accessed from disk
#' column-by-column; the response matrix Y is kept in memory.
#'
#' @references
#' Dayal, B. S. and MacGregor, J. F. (1997) Improved PLS algorithms.
#' \emph{Journal of Chemometrics}, \bold{11}, 73--85.
#'
#' @export
#' @method kernelpls IterableMatrix
kernelpls.IterableMatrix <- function(X, Y, ncomp, center = TRUE,
                                      stripped = FALSE, threads = NULL, ...) {
  if (!requireNamespace("BPCells", quietly = TRUE)) {
    rlang::abort("Package 'BPCells' is required for kernelpls.IterableMatrix but is not installed.")
  }

  if (!is.numeric(ncomp) || length(ncomp) != 1 || ncomp != as.integer(ncomp)) {
    rlang::abort("ncomp must be a single whole number")
  }
  ncomp <- as.integer(ncomp)
  if (ncomp <= 0) {
    rlang::abort("ncomp must be a positive integer")
  }

  threads <- as.integer(threads %||% getThreads())
  if (threads < 0L) rlang::abort("threads must be a non-negative integer")

  # Convert factor/character Y to a dummy indicator matrix
  if (is.factor(Y) || is.character(Y)) {
    Y <- stats::model.matrix(~ 0 + factor(Y))
    colnames(Y) <- sub("^factor\\(Y\\)", "", colnames(Y))
  }

  Y <- as.matrix(Y)

  if (!is.numeric(Y)) {
    rlang::abort("Y must be numeric (or a factor/character vector for classification)")
  }

  # BPCells' nrow()/ncol() already account for @transpose, so they always
  # return the user-facing dimensions: nrow = features, ncol = samples.
  n <- ncol(X)
  p <- nrow(X)
  q <- ncol(Y)

  if (n != nrow(Y)) {
    rlang::abort(sprintf("Number of samples in X (%d) must match nrow(Y) (%d)", n, nrow(Y)))
  }
  if (ncomp >= min(n, p)) {
    rlang::abort(sprintf("ncomp must be less than min(n, p) = %d", min(n, p)))
  }

  # The transpose flag tells C++ whether to swap row/col semantics.
  # BPCells default (features x samples) needs transpose=TRUE so that C++
  # treats rows as features and columns as samples.
  is_transposed <- !X@transpose

  if (threads == 0L) {
    threads <- max(parallel::detectCores(logical = FALSE), 1L, na.rm = TRUE)
  }

  n_splits <- min(threads * 4L, n)

  parallel_split <- getFromNamespace("parallel_split", "BPCells")
  iterate_matrix <- getFromNamespace("iterate_matrix", "BPCells")
  it <- X %>%
    BPCells::convert_matrix_type("double") %>%
    parallel_split(threads, n_splits) %>%
    iterate_matrix()

  result <- kernelpls_cpp(it, Y, as.integer(ncomp), center, stripped, is_transposed)

  # Reshape 2D coefficient matrix to 3D array (p x q x ncomp)
  comp_names <- paste(seq_len(ncomp), "comps")
  result$coefficients <- array(result$coefficients, dim = c(p, q, ncomp),
                               dimnames = list(NULL, NULL, comp_names))

  if (!stripped) {
    # Reshape fitted.values and residuals to 3D arrays (n x q x ncomp)
    result$fitted.values <- array(result$fitted.values, dim = c(n, q, ncomp),
                                  dimnames = list(NULL, NULL, comp_names))
    result$residuals     <- array(result$residuals, dim = c(n, q, ncomp),
                                  dimnames = list(NULL, NULL, comp_names))

    # Add class attributes for compatibility with pls package conventions
    class(result$scores)  <- class(result$Yscores) <- "scores"
    class(result$loadings) <- class(result$loading.weights) <-
      class(result$Yloadings) <- "loadings"
  }

  result
}


#' CPLS (Canonical PLS) for On-Disk Matrices
#'
#' Generic function for CPLS. The \code{IterableMatrix} method computes
#' CPLS for BPCells on-disk matrices.
#'
#' @param X Predictor matrix (or an IterableMatrix).
#' @param ... Arguments passed to methods.
#' @export
cppls_ondisk <- function(X, ...) UseMethod("cppls_ondisk")

#' CPLS (Canonical PLS) for IterableMatrix
#'
#' @description
#' Compute CPLS decomposition for large on-disk BPCells matrices.
#' CPLS uses canonical correlation analysis (CCA) to determine loading weights,
#' which can improve performance with multi-response or mixed-type response data.
#'
#' @param X An IterableMatrix (features x samples, the BPCells convention).
#'   The transpose is handled internally.
#' @param Y A matrix or vector of primary responses (samples x q).
#' @param Y.add A matrix of additional responses (samples x q_add), or NULL.
#'   Y.add provides information that guides loading weights but does not
#'   enter the regression coefficients.
#' @param ncomp Integer; number of PLS components to compute.
#' @param center Logical; whether to center X and Y (default \code{TRUE}).
#' @param stripped Logical; if \code{TRUE} return only coefficients, Xmeans, and
#'   Ymeans for speed (default \code{FALSE}).
#' @param w.tol Numeric; threshold for zeroing small loading weights (default 0).
#' @param X.tol Numeric; threshold for small-norm variable detection (default 1e-12).
#' @param threads Integer; number of threads for BPCells parallel computation.
#'   \code{NULL} (default) uses the package-wide thread count set by
#'   \code{\link{setThreads}} (\code{1} unless changed);
#'   \code{0} = auto-detect via \code{parallel::detectCores(logical = FALSE)};
#'   values \eqn{\ge 1} use that many threads.
#' @param ... Additional arguments (currently unused).
#'
#' @return A list with components:
#' \describe{
#'   \item{coefficients}{Regression coefficients (p x q x ncomp array).}
#'   \item{scores}{X-scores (n x ncomp matrix).}
#'   \item{loadings}{X-loadings (p x ncomp matrix).}
#'   \item{loading.weights}{Loading weights (p x ncomp matrix).}
#'   \item{Yscores}{Y-scores (n x ncomp matrix).}
#'   \item{Yloadings}{Y-loadings (q x ncomp matrix).}
#'   \item{projection}{Projection matrix (p x ncomp matrix).}
#'   \item{Xmeans}{Column means of X (length p).}
#'   \item{Ymeans}{Column means of Y (length q).}
#'   \item{fitted.values}{Fitted Y values (n x q x ncomp array).}
#'   \item{residuals}{Residuals (n x q x ncomp array).}
#'   \item{Xvar}{Variance explained in X per component (length ncomp).}
#'   \item{Xtotvar}{Total variance in X (scalar).}
#'   \item{gammas}{Power values per component (0.5 for all in this implementation).}
#'   \item{canonical.correlations}{Squared canonical correlations per component.}
#'   \item{A}{CCA weight vectors (q_full x ncomp matrix).}
#'   \item{smallNorm}{Indices of near-zero-norm variables.}
#' }
#' If \code{stripped = TRUE}, only \code{coefficients}, \code{Xmeans},
#' \code{Ymeans}, \code{gammas}, \code{canonical.correlations}, \code{A},
#' and \code{smallNorm} are returned.
#'
#' @export
#' @method cppls_ondisk IterableMatrix
cppls_ondisk.IterableMatrix <- function(X, Y, Y.add = NULL, ncomp,
                                         center = TRUE, stripped = FALSE,
                                         w.tol = 0, X.tol = 1e-12,
                                         threads = NULL, ...) {
  if (!requireNamespace("BPCells", quietly = TRUE)) {
    rlang::abort("Package 'BPCells' is required for cppls_ondisk.IterableMatrix but is not installed.")
  }

  if (!is.numeric(ncomp) || length(ncomp) != 1 || ncomp != as.integer(ncomp)) {
    rlang::abort("ncomp must be a single whole number")
  }
  ncomp <- as.integer(ncomp)
  if (ncomp <= 0) {
    rlang::abort("ncomp must be a positive integer")
  }

  threads <- as.integer(threads %||% getThreads())
  if (threads < 0L) rlang::abort("threads must be a non-negative integer")

  # Convert factor/character Y to a dummy indicator matrix
  if (is.factor(Y) || is.character(Y)) {
    Y <- stats::model.matrix(~ 0 + factor(Y))
    colnames(Y) <- sub("^factor\\(Y\\)", "", colnames(Y))
  }

  Y <- as.matrix(Y)

  if (!is.numeric(Y)) {
    rlang::abort("Y must be numeric (or a factor/character vector for classification)")
  }

  # Handle Y.add
  if (is.null(Y.add)) {
    Y_add_mat <- matrix(0.0, nrow = nrow(Y), ncol = 0)
  } else {
    if (is.factor(Y.add) || is.character(Y.add)) {
      Y.add <- stats::model.matrix(~ 0 + factor(Y.add))
      colnames(Y.add) <- sub("^factor\\(Y\\.add\\)", "", colnames(Y.add))
    }
    Y_add_mat <- as.matrix(Y.add)
    if (!is.numeric(Y_add_mat)) {
      rlang::abort("Y.add must be numeric (or a factor/character vector)")
    }
    if (nrow(Y_add_mat) != nrow(Y)) {
      rlang::abort(sprintf("nrow(Y.add) (%d) must match nrow(Y) (%d)",
                           nrow(Y_add_mat), nrow(Y)))
    }
  }

  # BPCells' nrow()/ncol() already account for @transpose, so they always
  # return the user-facing dimensions: nrow = features, ncol = samples.
  n <- ncol(X)
  p <- nrow(X)
  q <- ncol(Y)

  if (n != nrow(Y)) {
    rlang::abort(sprintf("Number of samples in X (%d) must match nrow(Y) (%d)", n, nrow(Y)))
  }
  if (ncomp >= min(n, p)) {
    rlang::abort(sprintf("ncomp must be less than min(n, p) = %d", min(n, p)))
  }

  is_transposed <- !X@transpose

  if (threads == 0L) {
    threads <- max(parallel::detectCores(logical = FALSE), 1L, na.rm = TRUE)
  }

  n_splits <- min(threads * 4L, n)

  parallel_split <- getFromNamespace("parallel_split", "BPCells")
  iterate_matrix <- getFromNamespace("iterate_matrix", "BPCells")
  it <- X %>%
    BPCells::convert_matrix_type("double") %>%
    parallel_split(threads, n_splits) %>%
    iterate_matrix()

  result <- cppls_cpp(it, Y, Y_add_mat, as.integer(ncomp), center, stripped,
                       is_transposed, as.double(w.tol), as.double(X.tol))

  # Reshape 2D coefficient matrix to 3D array (p x q x ncomp)
  comp_names <- paste(seq_len(ncomp), "comps")
  result$coefficients <- array(result$coefficients, dim = c(p, q, ncomp),
                               dimnames = list(NULL, NULL, comp_names))

  if (!stripped) {
    # Reshape fitted.values and residuals to 3D arrays (n x q x ncomp)
    result$fitted.values <- array(result$fitted.values, dim = c(n, q, ncomp),
                                  dimnames = list(NULL, NULL, comp_names))
    result$residuals     <- array(result$residuals, dim = c(n, q, ncomp),
                                  dimnames = list(NULL, NULL, comp_names))

    # Add class attributes for compatibility with pls package conventions
    class(result$scores)  <- class(result$Yscores) <- "scores"
    class(result$loadings) <- class(result$loading.weights) <-
      class(result$Yloadings) <- "loadings"
  }

  result
}

#' Predict Responses from a PLS Model
#'
#' Performs response prediction on new data using a PLS model previously
#' fitted with \code{\link{RunPLS}(..., save.model = TRUE)}. Mirrors the
#' core behaviour of \code{predict.mvr()} from the \pkg{pls} package.
#'
#' @param object A \code{DimReduc} object or \code{Seurat} object containing
#'   a PLS model fitted with \code{save.model = TRUE}.
#' @param newdata New observations: a cells x features numeric matrix (for
#'   the DimReduc method) or a Seurat object (for the Seurat method).
#' @param reduction (Seurat method only) Name of the PLS reduction slot.
#'   Default \code{"pls"}.
#' @param ncomp Integer vector of component counts to predict with. Each
#'   value produces a separate set of predictions using the cumulative
#'   coefficients for that many components. Default: maximum available.
#' @param assay (Seurat method only) Assay to use from \code{newdata}.
#'   Default: \code{DefaultAssay(newdata)}.
#' @param layer (Seurat method only) Layer to use from the assay. Default
#'   \code{"data"}.
#' @param ... Additional arguments (currently unused).
#'
#' @return A matrix (n x q) if a single \code{ncomp} value, or a 3D array
#'   (n x q x length(ncomp)) if multiple values are given.
#'
#' @details
#' When called on a Seurat object, \code{PredictPLS} automatically scales the
#' new data using per-gene mean and standard deviation stored from the training
#' \code{data} layer (saved by \code{RunPLS(..., save.model = TRUE)}). This
#' ensures that new observations are projected into the same coordinate system
#' used to fit the PLS model. If the model was saved without scaling parameters
#' (e.g., from an older version), a warning is emitted and prediction proceeds
#' without scaling.
#'
#' @export
#' @rdname PredictPLS
#'
#' @importFrom SeuratObject Loadings DefaultAssay LayerData
#' @importFrom utils head
#'
#' @examples
#' \dontrun{
#' # Fit PLS model with save.model = TRUE
#' obj <- RunPLS(obj, Y = "response", save.model = TRUE)
#'
#' # Predict on new Seurat object
#' predictions <- PredictPLS(obj, newdata = new_obj)
#' }
PredictPLS <- function(object, newdata, ...) {
  UseMethod("PredictPLS")
}

#' @rdname PredictPLS
#' @method PredictPLS DimReduc
#' @export
PredictPLS.DimReduc <- function(
    object,
    newdata,
    ncomp = NULL,
    ...
) {
  # Validate model exists

if (is.null(object@misc$model)) {
    stop("No saved model found in the DimReduc object. ",
         "Re-run RunPLS() with save.model = TRUE.")
  }

  # Dispatch to on-disk path if newdata is a BPCells IterableMatrix
  if (inherits(newdata, "IterableMatrix")) {
    return(.predict_pls_iterable(object, newdata, ncomp))
  }

  # Extract model components
  coefficients <- object@misc$model$coefficients
  Xmeans <- object@misc$model$Xmeans
  Ymeans <- object@misc$model$Ymeans

  # Get feature names from loadings
  model.features <- rownames(Loadings(object))
  p <- dim(coefficients)[1]

  # Coerce newdata to matrix if needed
  if (!is.matrix(newdata)) {
    newdata <- as.matrix(newdata)
  }

  # Feature alignment
  if (length(model.features) > 0 && !is.null(colnames(newdata))) {
    missing.features <- setdiff(model.features, colnames(newdata))
    if (length(missing.features) > 0) {
      stop("newdata is missing ", length(missing.features),
           " features required by the model: ",
           paste(head(missing.features, 5), collapse = ", "),
           if (length(missing.features) > 5) ", ..." else "")
    }
    newdata <- newdata[, model.features, drop = FALSE]
  } else {
    if (ncol(newdata) != p) {
      stop("newdata has ", ncol(newdata),
           " columns but model expects ", p, " features")
    }
  }

  # Validate ncomp
  ncomp_max <- dim(coefficients)[3]
  if (is.null(ncomp)) {
    ncomp <- ncomp_max
  }
  if (any(ncomp < 1) || any(ncomp > ncomp_max)) {
    stop("ncomp values must be between 1 and ", ncomp_max)
  }

  nobs <- nrow(newdata)
  q <- length(Ymeans)

  if (length(ncomp) == 1) {
    B <- coefficients[, , ncomp, drop = FALSE]
    dim(B) <- dim(B)[1:2]
    B0 <- Ymeans - crossprod(Xmeans, B)
    pred <- newdata %*% B + rep(B0, each = nobs)
    # Set dimnames
    rownames(pred) <- rownames(newdata)
    colnames(pred) <- names(Ymeans)
  } else {
    pred <- array(NA_real_, dim = c(nobs, q, length(ncomp)))
    for (i in seq_along(ncomp)) {
      k <- ncomp[i]
      B <- coefficients[, , k, drop = FALSE]
      dim(B) <- dim(B)[1:2]
      B0 <- Ymeans - crossprod(Xmeans, B)
      pred[, , i] <- newdata %*% B + rep(B0, each = nobs)
    }
    # Set dimnames
    comp.names <- dimnames(coefficients)[[3]]
    if (!is.null(comp.names)) {
      comp.labels <- comp.names[ncomp]
    } else {
      comp.labels <- paste(ncomp, "comps")
    }
    dimnames(pred) <- list(
      rownames(newdata),
      names(Ymeans),
      comp.labels
    )
  }

  return(pred)
}

#' @rdname PredictPLS
#' @method PredictPLS Seurat
#' @export
PredictPLS.Seurat <- function(
    object,
    newdata,
    reduction = "pls",
    ncomp = NULL,
    assay = NULL,
    layer = "data",
    ...
) {
  # Extract the DimReduc object
  reduction.obj <- object[[reduction]]
  if (!inherits(reduction.obj, "DimReduc")) {
    stop("reduction '", reduction, "' not found or is not a DimReduc object")
  }

  # Validate model exists
  if (is.null(reduction.obj@misc$model)) {
    stop("No saved model found in reduction '", reduction, "'. ",
         "Re-run RunPLS() with save.model = TRUE.")
  }

  # Get model features
  features <- rownames(Loadings(reduction.obj))

  # Resolve assay/layer in newdata
  assay <- assay %||% DefaultAssay(newdata)

  # Extract data matrix for matching features
  layer.data <- LayerData(newdata[[assay]], layer = layer)
  available.features <- rownames(layer.data)
  missing.features <- setdiff(features, available.features)
  if (length(missing.features) > 0) {
    stop("newdata is missing ", length(missing.features),
         " features required by the model: ",
         paste(head(missing.features, 5), collapse = ", "),
         if (length(missing.features) > 5) ", ..." else "")
  }

  # Apply training scaling parameters if available
  feat.mean <- reduction.obj@misc$model$feature.mean
  feat.sd <- reduction.obj@misc$model$feature.sd
  has.scaling <- !is.null(feat.mean) && !is.null(feat.sd)
  if (!has.scaling) {
    warning("PLS model does not contain feature.mean/feature.sd scaling ",
            "parameters. Prediction will proceed without scaling. ",
            "Re-run RunPLS() with save.model = TRUE to store scaling parameters.")
  }

  if (inherits(layer.data, "IterableMatrix")) {
    # Subset features lazily (features x samples layout), pass through
    new.mat <- layer.data[features, ]
    # Apply row-wise centering and scaling via BPCells lazy arithmetic
    if (has.scaling) {
      new.mat <- (new.mat - feat.mean[features]) / feat.sd[features]
    }
    return(PredictPLS.DimReduc(
      object = reduction.obj,
      newdata = new.mat,
      ncomp = ncomp,
      ...
    ))
  }

  # Dense path: materialise and transpose to cells x features
  new.mat <- t(as.matrix(layer.data[features, ]))
  if (has.scaling) {
    new.mat <- scale(new.mat, center = feat.mean[features], scale = feat.sd[features])
    attr(new.mat, "scaled:center") <- NULL
    attr(new.mat, "scaled:scale") <- NULL
  }

  # Delegate to DimReduc method
  PredictPLS.DimReduc(
    object = reduction.obj,
    newdata = new.mat,
    ncomp = ncomp,
    ...
  )
}

# Internal helper: predict from PLS model using BPCells IterableMatrix (on-disk)
# newdata is an IterableMatrix; prediction uses C++ bridge to avoid materialisation.
.predict_pls_iterable <- function(object, newdata, ncomp) {
  if (!requireNamespace("BPCells", quietly = TRUE)) {
    stop("Package 'BPCells' is required for IterableMatrix prediction but is not installed.")
  }

  # Extract model components
  coefficients <- object@misc$model$coefficients
  Xmeans <- object@misc$model$Xmeans
  Ymeans <- object@misc$model$Ymeans
  model.features <- rownames(Loadings(object))
  p <- dim(coefficients)[1]

  # Feature alignment
  # BPCells' rownames()/colnames() already account for @transpose.
  feat.names <- rownames(newdata)

  if (length(model.features) > 0 && !is.null(feat.names)) {
    missing.features <- setdiff(model.features, feat.names)
    if (length(missing.features) > 0) {
      stop("newdata is missing ", length(missing.features),
           " features required by the model: ",
           paste(head(missing.features, 5), collapse = ", "),
           if (length(missing.features) > 5) ", ..." else "")
    }
    # Subset/reorder lazily
    newdata <- newdata[model.features, ]
  } else {
    n.feat <- nrow(newdata)
    if (n.feat != p) {
      stop("newdata has ", n.feat,
           " features but model expects ", p, " features")
    }
  }

  # Determine transpose flag (same convention as kernelpls.IterableMatrix)
  is_transposed <- !newdata@transpose

  # Get sample count and names
  n.samples <- ncol(newdata)
  sample.names <- colnames(newdata)

  # Validate ncomp
  ncomp_max <- dim(coefficients)[3]
  if (is.null(ncomp)) {
    ncomp <- ncomp_max
  }
  if (any(ncomp < 1) || any(ncomp > ncomp_max)) {
    stop("ncomp values must be between 1 and ", ncomp_max)
  }

  q <- length(Ymeans)
  iterate_matrix <- getFromNamespace("iterate_matrix", "BPCells")

  if (length(ncomp) == 1) {
    B <- coefficients[, , ncomp, drop = FALSE]
    dim(B) <- dim(B)[1:2]
    B0 <- as.numeric(Ymeans - crossprod(Xmeans, B))

    it <- newdata %>%
      BPCells::convert_matrix_type("double") %>%
      iterate_matrix()
    pred <- predict_pls_cpp(it, B, B0, is_transposed)

    rownames(pred) <- sample.names
    colnames(pred) <- names(Ymeans)
  } else {
    pred <- array(NA_real_, dim = c(n.samples, q, length(ncomp)))
    for (i in seq_along(ncomp)) {
      k <- ncomp[i]
      B <- coefficients[, , k, drop = FALSE]
      dim(B) <- dim(B)[1:2]
      B0 <- as.numeric(Ymeans - crossprod(Xmeans, B))

      # Fresh iterator each time (take_unique_xptr consumes the previous one)
      it <- newdata %>%
        BPCells::convert_matrix_type("double") %>%
        iterate_matrix()
      pred[, , i] <- predict_pls_cpp(it, B, B0, is_transposed)
    }
    comp.names <- dimnames(coefficients)[[3]]
    if (!is.null(comp.names)) {
      comp.labels <- comp.names[ncomp]
    } else {
      comp.labels <- paste(ncomp, "comps")
    }
    dimnames(pred) <- list(
      sample.names,
      names(Ymeans),
      comp.labels
    )
  }

  return(pred)
}
