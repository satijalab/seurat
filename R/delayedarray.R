#%% DelayedArray / HDF5Array backend %%#
#
# Support for DelayedMatrix-backed layers (DelayedArray, HDF5Array, TENxMatrix).
# These backends process matrices block-wise and use 64-bit indexing, so they
# bypass the 2^31-1 element limit of dgCMatrix. The methods below mirror the
# existing BPCells IterableMatrix backend so the standard preprocessing
# pipeline (NormalizeData -> FindVariableFeatures -> ScaleData) works on
# DelayedMatrix layers without ever materializing the full matrix in memory.
NULL

#' Coerce to a DelayedMatrix
#'
#' Wrap a matrix-like object in a \code{\link[DelayedArray]{DelayedMatrix}} so it
#' can be used as an out-of-memory Seurat layer. Existing \code{DelayedMatrix}
#' objects are returned unchanged; other inputs (\code{matrix},
#' \code{\link[Matrix]{dgCMatrix}}, \code{IterableMatrix}) are wrapped, preserving
#' dimnames.
#'
#' @param x A matrix-like object
#' @param ... Ignored
#'
#' @return A \code{\link[DelayedArray]{DelayedMatrix}}
#'
#' @export
#' @concept utilities
#'
#' @examples
#' \dontrun{
#' mat <- Matrix::rsparsematrix(100, 20, 0.1)
#' rownames(mat) <- paste0("g", 1:100)
#' colnames(mat) <- paste0("c", 1:20)
#' d <- as.DelayedMatrix(mat)
#' }
#'
as.DelayedMatrix <- function(x, ...) {
  if (!requireNamespace('DelayedArray', quietly = TRUE)) {
    stop("Package 'DelayedArray' must be installed to create a DelayedMatrix.",
         call. = FALSE)
  }
  if (inherits(x = x, what = 'DelayedMatrix')) {
    return(x)
  }
  if (inherits(x = x, what = 'IterableMatrix')) {
    x <- as.sparse(x = x)
  }
  return(DelayedArray::DelayedArray(seed = x))
}

#' @method as.sparse DelayedMatrix
#' @export
#'
as.sparse.DelayedMatrix <- function(x, ...) {
  return(as(object = x, Class = 'dgCMatrix'))
}

#' @importFrom SeuratObject .CalcN
#'
#' @method .CalcN DelayedMatrix
#' @export
#'
.CalcN.DelayedMatrix <- function(object, ...) {
  return(list(
    nCount = DelayedMatrixStats::colSums2(x = object),
    nFeature = nrow(x = object) -
      DelayedMatrixStats::colCounts(x = object, value = 0)
  ))
}

#' @method LogNormalize DelayedMatrix
#' @export
#'
LogNormalize.DelayedMatrix <- function(
    data,
    scale.factor = 1e4,
    margin = 2L,
    verbose = TRUE,
    ...
) {
  # Divide each cell (column) by its total counts, then log1p-transform; all
  # operations stay delayed (block-wise), so the matrix is never materialized.
  col.sums <- DelayedMatrixStats::colSums2(x = data)
  col.sums[col.sums == 0] <- 1
  data <- DelayedArray::sweep(x = data, MARGIN = 2L, STATS = col.sums, FUN = '/')
  data <- log1p(x = data * scale.factor)
  return(data)
}

#' @rdname VST
#' @method VST DelayedMatrix
#' @importFrom SeuratObject EmptyDF
#' @export
#'
VST.DelayedMatrix <- function(
    data,
    margin = 1L,
    nselect = 2000L,
    span = 0.3,
    clip = NULL,
    verbose = TRUE,
    ...
) {
  nfeatures <- nrow(x = data)
  hvf.info <- EmptyDF(n = nfeatures)
  hvf.info$mean <- DelayedMatrixStats::rowMeans2(x = data)
  hvf.info$variance <- DelayedMatrixStats::rowVars(x = data)
  hvf.info$variance.expected <- 0L
  not.const <- hvf.info$variance > 0
  fit <- loess(
    formula = log10(x = variance) ~ log10(x = mean),
    data = hvf.info[not.const, , drop = TRUE],
    span = span
  )
  hvf.info$variance.expected[not.const] <- 10 ^ fit$fitted
  feature.mean <- hvf.info$mean
  feature.sd <- sqrt(x = hvf.info$variance.expected)
  feature.sd[feature.sd == 0] <- 1
  standard.max <- clip %||% sqrt(x = ncol(x = data))
  # Standardize every entry (including the implicit zeros), clip at standard.max,
  # and take the per-feature variance of the clipped standardized values.
  standardized <- DelayedArray::sweep(
    x = data, MARGIN = 1L, STATS = feature.mean, FUN = '-'
  )
  standardized <- DelayedArray::sweep(
    x = standardized, MARGIN = 1L, STATS = feature.sd, FUN = '/'
  )
  standardized <- DelayedArray::pmin2(e1 = standardized, e2 = standard.max)
  hvf.info$variance.standardized <- DelayedMatrixStats::rowVars(x = standardized)
  hvf.info$variable <- FALSE
  hvf.info$rank <- NA
  vf <- head(
    x = order(hvf.info$variance.standardized, decreasing = TRUE),
    n = nselect
  )
  hvf.info$variable[vf] <- TRUE
  hvf.info$rank[vf] <- seq_along(along.with = vf)
  rownames(x = hvf.info) <- rownames(x = data)
  return(hvf.info)
}
