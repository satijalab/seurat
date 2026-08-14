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
