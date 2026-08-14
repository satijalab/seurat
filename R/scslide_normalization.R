#' Apply Chi-squared normalization to a Seurat object
#'
#' This function applies Chi-squared normalization to a (sample-level) Seurat object. It calculates the expected density (count of cells) for each
#' landmark-sample pair assuming random distribution, then normalizes observed counts using:
#' (observed - expected) / sqrt(expected). The resulting values show how much each landmark's density deviates
#' from the expected baseline density values. Positive values indicate higher-than-expected density,
#' negative values indicate lower-than-expected density.
#'
#' @param object a count matrix (a landmark-sample density matrix), Assay, StdAssay, or Seurat object
#' @param assay Name of assay to use (for Seurat objects)
#' @param layer Layer to use (for StdAssay objects)
#' @param save Layer to save results to (for StdAssay objects)
#' @param verbose display progress bar for normalization procedure
#' @param ... Additional arguments
#'
#' @return Normalized object of the same class as input
#' @export
#' @concept normalization
#'
#' @examples
#' \dontrun{
#' # Normalize a sample-level Seurat object
#' sample_obj <- NormalizeChiSquared(sample_obj)
#' }
#'
#' @importFrom SeuratObject LayerData LayerData<-
#'
NormalizeChiSquared <- function(object, ...) {
  UseMethod(generic = 'NormalizeChiSquared', object = object)
}

#' @rdname NormalizeChiSquared
#' @method NormalizeChiSquared default
#' @export
NormalizeChiSquared.default <- function(
    object = NULL,
    verbose = TRUE,
    ...
){
  if (verbose) {
    cat("Performing Chi-squared-normalization\n")
  }
  #
  ct_matrix <- object
  # calculate the row sums and column sums
  col_Sums_ct <- SeuratObject::colSums(ct_matrix)
  row_Sums_ct <- SeuratObject::rowSums(ct_matrix)
  total_ct <- sum(ct_matrix)
  # the expected count matrix
  expect_ct_matrix <- t(outer(col_Sums_ct, row_Sums_ct)/total_ct)
  # calculate the Pearson residual
  relative_matrix <- (ct_matrix - expect_ct_matrix) / sqrt(expect_ct_matrix)
  # there are some 0(NA) in the relative_matrix, need to convert them
  relative_matrix[is.na(relative_matrix)] <- 0

  # a standard matrix object should be fine since the matrix will always be dense
  relative_matrix <- as.matrix(relative_matrix)
  #
  return(relative_matrix)
}

#' @rdname NormalizeChiSquared
#' @method NormalizeChiSquared Assay
#' @export
NormalizeChiSquared.Assay <- function(
    object,
    verbose = TRUE,
    ...
) {
  LayerData(object, layer = "data") <- NormalizeChiSquared(
    object = LayerData(object = object, layer = "counts"),
    verbose = verbose
  )
  return(object)
}

#' @rdname NormalizeChiSquared
#' @method NormalizeChiSquared StdAssay
#' @export
NormalizeChiSquared.StdAssay <- function(
    object,
    layer = 'counts',
    save = 'data',
    verbose = TRUE,
    ...
) {
  # TODO: need to implement the real assay5-based normalization (allows chisqNorm for multiple counts layers)
  LayerData(object, layer = "data") <- NormalizeChiSquared(
    object = LayerData(object = object, layer = "counts"),
    verbose = verbose
  )
  return(object)
}

#' @rdname NormalizeChiSquared
#' @method NormalizeChiSquared Seurat
#' @export
NormalizeChiSquared.Seurat <- function(
    object,
    assay = NULL,
    verbose = TRUE,
    ...
) {
  assay <- assay %||% DefaultAssay(object = object)
  assay <- match.arg(arg = assay, choices = Assays(object = object))

  assay.data <- NormalizeChiSquared(
    object = object[[assay]],
    verbose = verbose
  )
  object[[assay]] <- assay.data
  object <- LogSeuratCommand(object = object)
  return(object)
}