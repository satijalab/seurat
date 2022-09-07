#' @include generics.R
#' @include preprocessing.R
#' @importFrom stats loess
#' @importFrom methods slot
#' @importFrom SeuratObject .MARGIN .SparseSlots
#' @importFrom utils txtProgressBar setTxtProgressBar
#'
NULL

hvf.methods <- list()

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Generics
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Methods for Seurat-defined generics
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @importFrom rlang is_quosure quo_get_env quo_get_expr
#' @method FindVariableFeatures default
#' @export
#'
FindVariableFeatures.default <- function(
  object,
  method = VST,
  nselect = 2000L,
  verbose = TRUE,
  ...
) {
  if (is_quosure(x = method)) {
    method <- eval(
      expr = quo_get_expr(quo = method),
      envir = quo_get_env(quo = method)
    )
  }
  if (is.character(x = method)) {
    method <- get(x = method)
  }
  if (!is.function(x = method)) {
    stop(
      "'method' must be a function for calculating highly variable features",
      call. = FALSE
    )
  }
  return(method(
    data = object,
    nselect = nselect,
    verbose = verbose,
    ...
  ))
}

g <- function(x, method = VST) {
  method <- enquo(arg = method)
  FindVariableFeatures(object = x, method = method, layer = 'counts')
}

#' @importFrom rlang as_name enquo is_quosure
#' @importFrom SeuratObject DefaultLayer Features Key Layers
#'
#' @method FindVariableFeatures StdAssay
#' @export
#'
FindVariableFeatures.StdAssay <- function(
  object,
  method = VST,
  nselect = 2000L,
  layer = NULL,
  span = 0.3,
  clip = NULL,
  key = NULL,
  verbose = TRUE,
  ...
) {
  layer <- unique(x = layer) %||% DefaultLayer(object = object)
  layer <- Layers(object = object, search = layer)
  if (is.null(x = key)) {
    false <- function(...) {
      return(FALSE)
    }
    key <- if (tryCatch(expr = is_quosure(x = method), error = false)) {
      method
    } else if (is.function(x = method)) {
      substitute(expr = method)
    } else if (is.call(x = enquo(arg = method))) {
      enquo(arg = method)
    } else if (is.character(x = method)) {
      method
    } else {
      parse(text = method)
    }
    key <- .Abbrv(x = as_name(x = key))
  }
  warn.var <- warn.rank <- TRUE
  for (i in seq_along(along.with = layer)) {
    if (isTRUE(x = verbose)) {
      message("Finding variable features for layer ", layer[i])
    }
    data <- LayerData(object = object, layer = layer[i], fast = TRUE)
    f <- if (inherits(x = data, what = 'V3Matrix')) {
      FindVariableFeatures.default
    } else {
      FindVariableFeatures
    }
    hvf.info <- f(
      object = data,
      method = method,
      nselect = nselect,
      span = span,
      clip = clip,
      verbose = verbose,
      ...
    )
    if (warn.var) {
      if (!'variable' %in% colnames(x = hvf.info) || !is.logical(x = hvf.info$variable)) {
        warning(
          "No variable feature indication in HVF info for method ",
          key,
          ", `VariableFeatures` will not work",
          call. = FALSE,
          immediate. = TRUE
        )
        warn.var <- FALSE
      }
    } else if (warn.rank && !'rank' %in% colnames(x = hvf.info)) {
      warning(
        "No variable feature rank in HVF info for method ",
        key,
        ", `VariableFeatures` will return variable features in assay order",
        call. = FALSE,
        immediate. = TRUE
      )
      warn.rank <- FALSE
    }
    colnames(x = hvf.info) <- paste(
      'vf',
      key,
      layer[i],
      colnames(x = hvf.info),
      sep = '_'
    )
    rownames(x = hvf.info) <- Features(x = object, layer = layer[i])
    object[[colnames(x = hvf.info)]] <- hvf.info
  }
  return(object)
}

#' @importFrom rlang enquo
#' @method FindVariableFeatures Seurat5
#' @export
#'
FindVariableFeatures.Seurat5 <- function(
  object,
  assay = NULL,
  method = VST,
  nselect = 2000L,
  layer = NULL,
  span = 0.3,
  clip = NULL,
  key = NULL,
  verbose = TRUE,
  ...
) {
  assay <- assay[1L] %||% DefaultAssay(object = object)
  assay <- match.arg(arg = assay, choices = Assays(object = object))
  method <- enquo(arg = method)
  object[[assay]] <- FindVariableFeatures(
    object = object[[assay]],
    method = method,
    nselect = nselect,
    layer = layer,
    span = span,
    clip = clip,
    key = key,
    verbose = verbose,
    ...
  )
  return(object)
}

#' @rdname LogNormalize
#' @method LogNormalize default
#' @export
#'
LogNormalize.default <- function(
  data,
  scale.factor = 1e4,
  margin = 2L,
  verbose = TRUE,
  ...
) {
  margin <- SeuratObject:::.CheckFmargin(fmargin = margin)
  ncells <- dim(x = data)[margin]
  if (isTRUE(x = verbose)) {
    pb <- txtProgressBar(file = stderr(), style = 3)
  }
  for (i in seq_len(length.out = ncells)) {
    x <- if (margin == 1L) {
      data[i, ]
    } else {
      data[, i]
    }
    xnorm <- log1p(x = x / sum(x) * scale.factor)
    if (margin == 1L) {
      data[i, ] <- xnorm
    } else {
      data[, i] <- xnorm
    }
    if (isTRUE(x = verbose)) {
      setTxtProgressBar(pb = pb, value = i / ncells)
    }
  }
  if (isTRUE(x = verbose)) {
    close(con = pb)
  }
  return(data)
}

#' @method LogNormalize DelayedMatrix
#' @export
#'
LogNormalize.DelayedMatrix <- function(
  data,
  scale.factor = 1e4,
  margin = 2L,
  verbose = TRUE,
  sink = NULL,
  ...
) {
  check_installed(
    pkg = 'DelayedArray',
    reason = 'for working with delayed matrices'
  )
  if (is.null(x = sink)) {
    sink <- DelayedArray::AutoRealizationSink(
      dim = dim(x = data),
      dimnames = dimnames(x = data),
      as.sparse = DelayedArray::is_sparse(x = data)
    )
  }
  if (!inherits(x = sink, what = 'RealizationSink')) {
    abort(message = "'sink' must be a RealizationSink")
  } else if (!all(dim(x = sink) == dim(x = data))) {
    abort(message = "'sink' must be the same size as 'data'")
  }
  if (!margin %in% c(1L, 2L)) {
    abort(message = "'margin' must be 1 or 2")
  }
  grid <- if (margin == 1L) {
    DelayedArray::rowAutoGrid(x = data)
  } else {
    DelayedArray::colAutoGrid(x = data)
  }
  sparse <- DelayedArray::is_sparse(x = data)
  if (isTRUE(x = verbose)) {
    pb <- txtProgressBar(file = stderr(), style = 3)
  }
  for (i in seq_len(length.out = length(x = grid))) {
    vp <- grid[[i]]
    x <- DelayedArray::read_block(x = data, viewport = vp, as.sparse = sparse)
    x <- LogNormalize(
      data = x,
      scale.factor = scale.factor,
      margin = margin,
      verbose = FALSE,
      ...
    )
    DelayedArray::write_block(sink = sink, viewport = vp, block = x)
    if (isTRUE(x = verbose)) {
      setTxtProgressBar(pb = pb, value = i / length(x = grid))
    }
  }
  if (isTRUE(x = verbose)) {
    close(con = pb)
  }
  DelayedArray::close(con = sink)
  return(as(object = sink, Class = "DelayedArray"))
}

#' @method LogNormalize H5ADMatrix
#' @export
#'
LogNormalize.H5ADMatrix <- function(
  data,
  scale.factor = 1e4,
  margin = 2L,
  verbose = TRUE,
  layer = 'data',
  ...
) {
  results <- LogNormalize.HDF5Matrix(
    data = data,
    scale.factor = scale.factor,
    margin = margin,
    verbose = verbose,
    layer = file.path('layers', layer, fsep = '/'),
    ...
  )
  rpath <- slot(object = slot(object = results, name = 'seed'), name = 'filepath')
  return(HDF5Array::H5ADMatrix(filepath = rpath, layer = layer))
}

#' @method LogNormalize HDF5Matrix
#' @export
#'
LogNormalize.HDF5Matrix <- function(
  data,
  scale.factor = 1e4,
  margin = 2L,
  verbose = TRUE,
  layer = 'data',
  ...
) {
  check_installed(pkg = 'HDF5Array', reason = 'for working with HDF5 matrices')
  fpath <- slot(object = slot(object = data, name = 'seed'), name = 'filepath')
  if (.DelayedH5DExists(object = data, path = layer)) {
    rhdf5::h5delete(file = fpath, name = layer)
    dpath <- file.path(
      dirname(path = layer),
      paste0('.', basename(layer), '_dimnames'),
      fsep = '/'
    )
    rhdf5::h5delete(file = fpath, name = dpath)
  }
  sink <- HDF5Array::HDF5RealizationSink(
    dim = dim(x = data),
    dimnames = dimnames(x = data),
    as.sparse = DelayedArray::is_sparse(x = data),
    filepath = fpath,
    name = layer
  )
  return(LogNormalize.DelayedMatrix(
    data = data,
    scale.factor = scale.factor,
    margin = margin,
    verbose = verbose,
    sink = sink,
    ...
  ))
}

#' @method LogNormalize SparseArraySeed
#' @export
#'
LogNormalize.SparseArraySeed <- function(
  data,
  scale.factor = 1e4,
  margin = 2L,
  return.seed = TRUE,
  verbose= TRUE,
  ...
) {
  check_installed(
    pkg = 'DelayedArray',
    reason = 'for working with SparseArraySeeds'
  )
  data <- LogNormalize(
    data = as(object = data, Class = 'CsparseMatrix'),
    scale.factor = scale.factor,
    margin = margin,
    verbose = verbose,
    ...
  )
  if (!isFALSE(x = return.seed)) {
    data <- as(object = data, Class = 'SparseArraySeed')
  }
  return(data)
}

#' @importFrom SeuratObject IsSparse
#'
#' @method NormalizeData default
#' @export
#'
NormalizeData.default <- function(
  object,
  method = c('LogNormalize'),
  scale.factor = 1e4,
  cmargin = 2L,
  margin = 1L,
  verbose = TRUE,
  ...
) {
  method <- method[1L]
  method <- match.arg(arg = method)
  # TODO: enable parallelization via future
  normalized <- switch(
    EXPR = method,
    'LogNormalize' = {
      if (IsSparse(x = object) && .MARGIN(object = object) == cmargin) {
        .SparseNormalize(
          data = object,
          scale.factor = scale.factor,
          verbose = verbose
        )
      } else {
        LogNormalize(
          data = object,
          scale.factor = scale.factor,
          margin = cmargin,
          verbose = verbose,
          ...
        )
      }
    }
  )
  return(normalized)
}

.DelayedH5DExists <- function(object, path) {
  check_installed(pkg = 'HDF5Array', reason = 'for working with HDF5 files')
  if (!inherits(x = object, what = c('HDF5Array', 'H5ADMatrix'))) {
    abort(message = "'object' must be an HDF5Array or H5ADMatrix")
  }
  on.exit(expr = rhdf5::h5closeAll(), add = TRUE)
  fpath <- slot(object = slot(object = object, name = 'seed'), name = 'filepath')
  h5loc <- rhdf5::H5Fopen(
    name = fpath,
    flags = 'H5F_ACC_RDWR',
    fapl = NULL,
    native = FALSE
  )
  return(rhdf5::H5Lexists(h5loc = h5loc, name = path))
}

# #' @method NormalizeData DelayedArray
# #' @export
# #'
# NormalizeData.DelayedArray <- function(
#   object,
#   method = c('LogNormalize'),
#   scale.factor = 1e4,
#   cmargin = 2L,
#   margin = 1L,
#   layer = 'data',
#   verbose = TRUE,
#   ...
# ) {
#   method <- arg_match(arg = method)
#   normalized <- switch(
#     EXPR = method,
#     LogNormalize = LogNormalize(
#       data = object,
#       scale.factor = scale.factor,
#       margin = 2L,
#       verbose = TRUE,
#       layer = layer,
#       ...
#     )
#   )
#   return(normalized)
# }

#' @importFrom SeuratObject Cells DefaultLayer DefaultLayer<- Features
#' LayerData LayerData<-
#'
#' @method NormalizeData StdAssay
#' @export
#'
NormalizeData.StdAssay <- function(
  object,
  method = 'LogNormalize',
  scale.factor = 1e4,
  margin = 1L,
  layer = NULL,
  save = 'data',
  default = TRUE,
  verbose = TRUE,
  ...
) {
  olayer <- layer <- unique(x = layer) %||% DefaultLayer(object = object)
  layer <- Layers(object = object, search = layer)
  if (save == DefaultLayer(object = object)) {
    default <- FALSE
  }
  if (length(x = save) != length(x = layer)) {
    save <- make.unique(names = gsub(
      pattern = olayer,
      replacement = save,
      x = layer
    ))
  }
  for (i in seq_along(along.with = layer)) {
    l <- layer[i]
    if (isTRUE(x = verbose)) {
      message("Normalizing layer: ", l)
    }
    LayerData(
      object = object,
      layer = save[i],
      features = Features(x = object, layer = l),
      cells = Cells(x = object, layer = l)
    ) <- NormalizeData(
      object = LayerData(object = object, layer = l, fast = NA),
      method = method,
      scale.factor = scale.factor,
      margin = margin,
      verbose = verbose,
      layer = save,
      ...
    )
  }
  if (isTRUE(x = default)) {
    DefaultLayer(object = object) <- save
  }
  gc(verbose = FALSE)
  return(object)
}

#' @importFrom SeuratObject DefaultAssay
#'
#' @method NormalizeData Seurat5
#' @export
#'
NormalizeData.Seurat5 <- function(
  object,
  assay = NULL,
  method = 'LogNormalize',
  scale.factor = 1e4,
  margin = 1L,
  layer = NULL,
  save = 'data',
  default = TRUE,
  verbose = TRUE,
  ...
) {
  assay <- assay[1L] %||% DefaultAssay(object = object)
  assay <- match.arg(arg = assay, choices = Assays(object = object))
  object[[assay]] <- NormalizeData(
    object = object[[assay]],
    method = method,
    scale.factor = scale.factor,
    margin = margin,
    layer = layer,
    save = save,
    default = default,
    verbose = verbose,
    ...
  )
  return(object)
}

#' @importFrom SeuratObject StitchMatrix
#'
#' @method ScaleData StdAssay
#' @export
#'
ScaleData.StdAssay <- function(
  object,
  features = NULL,
  layer = NULL,
  vars.to.regress = NULL,
  latent.data = NULL,
  split.by = NULL,
  model.use = 'linear',
  use.umi = FALSE,
  do.scale= TRUE,
  do.center = TRUE,
  scale.max = 10,
  block.size = 1000,
  min.cells.to.block = 3000,
  save = 'scale.data',
  verbose = TRUE,
  ...
) {
  use.umi <- ifelse(test = model.use != 'linear', yes = TRUE, no = use.umi)
  layer <- Layers(object = object, search = layer)
  if (isTRUE(x = use.umi)) {
    message("'use.umi' is TRUE, please make sure 'layer' specifies raw counts")
  }
  features <- features %||% Reduce(
    f = union,
    x = lapply(
      X = layer,
      FUN = function(x) {
        return(VariableFeatures(object = object, layer = layer))
      }
    )
  )
  if (!length(x = features)) {
    features <- Reduce(f = union, x = lapply(X = layer, FUN = Features, x = object))
  }
  ldata <- if (length(x = layer) > 1L) {
    StitchMatrix(
      x = LayerData(object = object, layer = layer[1L], features = features),
      y = lapply(
        X = layer[2:length(x = layer)],
        FUN = LayerData,
        object = object,
        features = features
      ),
      rowmap = slot(object = object, name = 'features')[features, layer],
      colmap = slot(object = object, name = 'cells')[, layer]
    )
  } else {
    LayerData(object = object, layer = layer, features = features)
  }
  LayerData(object = object, layer = save, features = features) <- ScaleData(
    object = ldata,
    features = features,
    vars.to.regress = vars.to.regress,
    latent.data = latent.data,
    split.by = split.by,
    model.use = model.use,
    use.umi = use.umi,
    do.scale = do.scale,
    do.center = do.center,
    scale.max = scale.max,
    block.size = block.size,
    min.cells.to.block = min.cells.to.block,
    verbose = verbose,
    ...
  )
  return(object)
}

#' @importFrom rlang is_scalar_character
#'
#' @method ScaleData Seurat5
#' @export
#'
ScaleData.Seurat5 <- function(
  object,
  features = NULL,
  assay = NULL,
  layer = NULL,
  vars.to.regress = NULL,
  split.by = NULL,
  model.use = 'linear',
  use.umi = FALSE,
  do.scale = TRUE,
  do.center = TRUE,
  scale.max = 10,
  block.size = 1000,
  min.cells.to.block = 3000,
  verbose = TRUE,
  ...
) {
  assay <- assay %||% DefaultAssay(object = object)
  if (!is.null(x = vars.to.regress)) {
    vars.to.regress <- intersect(x = vars.to.regress, y = names(x = object[[]]))
  }
  latent.data <- if (length(x = vars.to.regress)) {
    object[[vars.to.regress]]
  } else {
    NULL
  }
  if (is_scalar_character(x = split.by)) {
    split.by <- object[[split.by]]
  }
  object[[assay]] <- ScaleData(
    object = object[[assay]],
    features = features,
    layer = layer,
    vars.to.regress = vars.to.regress,
    latent.data = latent.data,
    split.by = split.by,
    model.use = model.use,
    use.umi = use.umi,
    do.scale = do.scale,
    do.center = do.center,
    scale.max = scale.max,
    min.cells.to.block = min.cells.to.block,
    verbose = verbose,
    ...
  )
  return(object)
}

#' @rdname VST
#' @method VST default
#' @export
#'
VST.default <- function(
  data,
  margin = 1L,
  nselect = 2000L,
  span = 0.3,
  clip = NULL,
  ...
) {
  .NotYetImplemented()
}

#' @method VST DelayedMatrix
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
  check_installed(
    pkg = 'DelayedArray',
    reason = 'for working with delayed matrices'
  )
  if (!margin %in% c(1L, 2L)) {
    abort(message = "'margin' must be 1 or 2")
  }
  grid <- if (margin == 1L) {
    DelayedArray::rowAutoGrid(x = data)
  } else {
    DelayedArray::colAutoGrid(x = data)
  }
  nfeatures <- dim(x = data)[margin]
  ncells <- dim(x = data)[-margin]
  # hvf.info <- SeuratObject::EmptyDF(n = nfeatures)
  hvf.info <- vector(mode = 'list', length = length(x = grid))
  sparse <- DelayedArray::is_sparse(x = data)
  # Calculate feature means
  # if (isTRUE(x = verbose)) {
  #   inform(message = "Calculating feature means")
  # }
  # hvf.info$mean <- if (margin == 1L) {
  #   DelayedArray::rowMeans(x = data)
  # } else {
  #   DelayedArray::colMeans(x = data)
  # }
  # Calculate variance
  # hvf.info$variance <- NA_real_
  if (isTRUE(x = verbose)) {
    # inform(message = "Calculating feature variances")
    inform(message = "Identifying variable features")
    pb <- txtProgressBar(style = 3L, file = stderr())
  }
  for (i in seq_len(length.out = length(x = grid))) {
    vp <- grid[[i]]
    idx <- seq.int(
      from = IRanges::start(x = slot(object = vp, name = 'ranges')[margin]),
      to = IRanges::end(x = slot(object = vp, name = 'ranges')[margin])
    )
    x <- DelayedArray::read_block(x = data, viewport = vp, as.sparse = sparse)
    if (isTRUE(x = sparse)) {
      x <- as(object = x, Class = "CsparseMatrix")
    }
    hvf.info[[i]] <- VST(
      data = x,
      margin = margin,
      nselect = floor(x = nselect / length(x = grid)),
      span = span,
      clip = clip,
      verbose = FALSE,
      ...
    )
  #   if (margin == 2L) {
  #     x <- t(x = x)
  #   }
  #   mu <- hvf.info$mean[idx]
  #   hvf.info$variance[idx] <- rowSums(x = ((x - mu) ^ 2) / (ncells - 1L))
  #   # hvf.info$variance[idx] <- vapply(
  #   #   X = seq_along(along.with = mu),
  #   #   FUN = function(j) {
  #   #     y <- if (margin == 1L) {
  #   #       x[j, ]
  #   #     } else {
  #   #       x[, j]
  #   #     }
  #   #     y <- y - mu[j]
  #   #     return(sum(y ^ 2) / (ncells - 1L))
  #   #   },
  #   #   FUN.VALUE = numeric(length = 1L)
  #   # )
  #   if (isTRUE(x = verbose)) {
  #     setTxtProgressBar(pb = pb, value = i / length(x = grid))
  #   }
  # }
  # if (isTRUE(x = verbose)) {
  #   close(con = pb)
  # }
  # hvf.info$variance.expected <- 0
  # not.const <- hvf.info$variance > 0
  # fit <- loess(
  #   formula = log10(x = variance) ~ log10(x = mean),
  #   data = hvf.info[not.const, , drop = FALSE],
  #   span = span
  # )
  # hvf.info$variance.expected[not.const] <- 10 ^ fit$fitted
  # # Calculate standardized variance
  # hvf.info$variance.standardized <- NA_real_
  # if (isTRUE(x = verbose)) {
  #   inform(
  #     message = "Calculating feature variances of standardized and clipped values"
  #   )
  #   pb <- txtProgressBar(style = 3L, file = stderr())
  # }
  # clip <- clip %||% sqrt(x = ncells)
  # for (i in seq_len(length.out = length(x = grid))) {
  #   vp <- grid[[i]]
  #   idx <- seq.int(
  #     from = IRanges::start(x = slot(object = vp, name = 'ranges')[margin]),
  #     to = IRanges::end(x = slot(object = vp, name = 'ranges')[margin])
  #   )
  #   x <- DelayedArray::read_block(x = data, viewport = vp, as.sparse = sparse)
  #   if (isTRUE(x = sparse)) {
  #     x <- as(object = x, Class = "CsparseMatrix")
  #   }
  #   if (margin == 2L) {
  #     x <- t(x = x)
  #   }
  #   mu <- hvf.info$mean[idx]
  #   sd <- sqrt(x = hvf.info$variance.expected[idx])
  #   hvf.info$variance.standardized[idx] <- 0
  #   sdn <- which(x = sd != 0)
  #   hvf.info$variance.standardized[idx[sdn]] <- rowSums(x = (((x[sdn, ] - mu[sdn]) / sd[sdn]) ^ 2) / (ncells - 1L))
  #   # hvf.info$variance.standardized[idx] <- vapply(
  #   #   X = seq_along(along.with = mu),
  #   #   FUN = function(j) {
  #   #     if (sd[j] == 0) {
  #   #       return(0)
  #   #     }
  #   #     y <- if (margin == 1L) {
  #   #       x[j, ]
  #   #     } else {
  #   #       x[, j]
  #   #     }
  #   #     y <- y - mu[j]
  #   #     y <- y / sd[j]
  #   #     y[y > clip] <- clip
  #   #     return(sum(y ^ 2) / (ncells - 1L))
  #   #   },
  #   #   FUN.VALUE = numeric(length = 1L)
  #   # )
    if (isTRUE(x = verbose)) {
      setTxtProgressBar(pb = pb, value = i / length(x = grid))
    }
  }
  if (isTRUE(x = verbose)) {
    close(con = pb)
  }
  hvf.info <- do.call(what = 'rbind', args = hvf.info)
  # Set variable status
  hvf.info$variable <- FALSE
  hvf.info$rank <- NA_integer_
  vs <- hvf.info$variance.standardized
  vs[vs == 0] <- NA
  vf <- head(
    x = order(vs, decreasing = TRUE),
    n = nselect
  )
  hvf.info$variable[vf] <- TRUE
  hvf.info$rank[vf] <- seq_along(along.with = vf)
  return(hvf.info)
}

#' @importFrom Matrix rowMeans
#'
#' @rdname VST
#' @method VST dgCMatrix
#' @export
#'
VST.dgCMatrix <- function(
  data,
  margin = 1L,
  nselect = 2000L,
  span = 0.3,
  clip = NULL,
  verbose = TRUE,
  ...
) {
  nfeatures <- nrow(x = data)
  hvf.info <- SeuratObject:::EmptyDF(n = nfeatures)
  # Calculate feature means
  hvf.info$mean <- Matrix::rowMeans(x = data)
  # Calculate feature variance
  hvf.info$variance <- SparseRowVar2(
    mat = data,
    mu = hvf.info$mean,
    display_progress = verbose
  )
  hvf.info$variance.expected <- 0L
  not.const <- hvf.info$variance > 0
  fit <- loess(
    formula = log10(x = variance) ~ log10(x = mean),
    data = hvf.info[not.const, , drop = TRUE],
    span = span
  )
  hvf.info$variance.expected[not.const] <- 10 ^ fit$fitted
  hvf.info$variance.standardized <- SparseRowVarStd(
    mat = data,
    mu = hvf.info$mean,
    sd = sqrt(x = hvf.info$variance.expected),
    vmax = clip %||% sqrt(x = ncol(x = data)),
    display_progress = verbose
  )
  # Set variable features
  hvf.info$variable <- FALSE
  hvf.info$rank <- NA
  vf <- head(
    x = order(hvf.info$variance.standardized, decreasing = TRUE),
    n = nselect
  )
  hvf.info$variable[vf] <- TRUE
  hvf.info$rank[vf] <- seq_along(along.with = vf)
  return(hvf.info)
}

#' @rdname VST
#' @method VST matrix
#' @export
#'
VST.matrix <- function(
  data,
  margin = 1L,
  nselect = 2000L,
  span = 0.3,
  clip = NULL,
  ...
) {
  return(VST(
    data = as.sparse(x = data),
    margin = margin,
    nselect = nselect,
    span = span,
    clip = clip,
    ...
  ))
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Methods for R-defined generics
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Internal
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

.FeatureVar <- function(
  data,
  mu,
  fmargin = 1L,
  standardize = FALSE,
  sd = NULL,
  clip = NULL,
  verbose = TRUE
) {
  fmargin <- SeuratObject:::.CheckFmargin(fmargin = fmargin)
  ncells <- dim(x = data)[-fmargin]
  nfeatures <- dim(x = data)[fmargin]
  fvars <- vector(mode = 'numeric', length = nfeatures)
  if (length(x = mu) != nfeatures) {
    stop("Wrong number of feature means provided")
  }
  if (isTRUE(x = standardize)) {
    clip <- clip %||% sqrt(x = ncells)
    if (length(x = sd) != nfeatures) {
      stop("Wrong number of standard deviations")
    }
  }
  if (isTRUE(x = verbose)) {
    msg <- 'Calculating feature variances'
    if (isTRUE(x = standardize)) {
      msg <- paste(msg, 'of standardized and clipped values')
    }
    message(msg)
    pb <- txtProgressBar(style = 3, file = stderr())
  }
  for (i in seq_len(length.out = nfeatures)) {
    if (isTRUE(x = standardize) && sd[i] == 0) {
      if (isTRUE(x = verbose)) {
        setTxtProgressBar(pb = pb, value = i / nfeatures)
      }
      next
    }
    x <- if (fmargin == 1L) {
      data[i, , drop = TRUE]
    } else {
      data[, i, drop = TRUE]
    }
    x <- x - mu[i]
    if (isTRUE(x = standardize)) {
      x <- x / sd[i]
      x[x > clip] <- clip
    }
    fvars[i] <- sum(x ^ 2) / (ncells - 1L)
    if (isTRUE(x = verbose)) {
      setTxtProgressBar(pb = pb, value = i / nfeatures)
    }
  }
  if (isTRUE(x = verbose)) {
    close(con = pb)
  }
  return(fvars)
}

.Mean <- function(data, margin = 1L) {
  nout <- dim(x = data)[margin]
  nobs <- dim(x = data)[-margin]
  means <- vector(mode = 'numeric', length = nout)
  for (i in seq_len(length.out = nout)) {
    x <- if (margin == 1L) {
      data[i, , drop = TRUE]
    } else {
      data[, i, drop = TRUE]
    }
    means[i] <- sum(x) / nobs
  }
  return(means)
}

.SparseNormalize <- function(data, scale.factor = 1e4, verbose = TRUE) {
  entryname <- .SparseSlots(x = data, type = 'entries')
  p <- slot(object = data, name = .SparseSlots(x = data, type = 'pointers'))
  if (p[1L] == 0) {
    p <- p + 1L
  }
  np <- length(x = p) - 1L
  if (isTRUE(x = verbose)) {
    pb <- txtProgressBar(style = 3L, file = stderr())
  }
  for (i in seq_len(length.out = np)) {
    idx <- seq.int(from = p[i], to = p[i + 1] - 1L)
    xidx <- slot(object = data, name = entryname)[idx]
    slot(object = data, name = entryname)[idx] <- log1p(
      x = xidx / sum(xidx) * scale.factor
    )
    if (isTRUE(x = verbose)) {
      setTxtProgressBar(pb = pb, value = i / np)
    }
  }
  if (isTRUE(x = verbose)) {
    close(con = pb)
  }
  return(data)
}

#' @param data A sparse matrix
#' @param mu A vector of feature means
#' @param fmargin Feature margin
#' @param standardize Standardize matrix rows prior to calculating variances
#' @param sd If standardizing, a vector of standard deviations to
#' standardize with
#' @param clip Set upper bound for standardized variances; defaults to the
#' square root of the number of cells
#' @param verbose Show progress updates
#'
#' @keywords internal
#'
#' @noRd
#'
.SparseFeatureVar <- function(
  data,
  mu,
  fmargin = 1L,
  standardize = FALSE,
  sd = NULL,
  clip = NULL,
  verbose = TRUE
) {
  fmargin <- SeuratObject:::.CheckFmargin(fmargin = fmargin)
  if (fmargin != .MARGIN(object = data)) {
    data <- t(x = data)
    fmargin <- .MARGIN(object = data)
  }
  entryname <- .SparseSlots(x = data, type = 'entries')
  p <- slot(object = data, name = .SparseSlots(x = data, type = 'pointers'))
  if (p[1L] == 0) {
    p <- p + 1L
  }
  np <- length(x = p) - 1L
  ncells <- dim(x = data)[-fmargin]
  fvars <- vector(mode = 'numeric', length = np)
  if (length(x = mu) != np) {
    stop("Wrong number of feature means provided")
  }
  if (isTRUE(x = standardize)) {
    clip <- clip %||% sqrt(x = ncells)
    if (length(x = sd) != np) {
      stop("Wrong number of standard deviations provided")
    }
  }
  if (isTRUE(x = verbose)) {
    msg <- 'Calculating feature variances'
    if (isTRUE(x = standardize)) {
      msg <- paste(msg, 'of standardized and clipped values')
    }
    message(msg)
    pb <- txtProgressBar(style = 3, file = stderr())
  }
  for (i in seq_len(length.out = np)) {
    if (isTRUE(x = standardize) && sd[i] == 0) {
      if (isTRUE(x = verbose)) {
        setTxtProgressBar(pb = pb, value = i / np)
      }
      next
    }
    idx <- seq.int(from = p[i], to = p[i + 1L] - 1L)
    xidx <- slot(object = data, name = entryname)[idx] - mu[i]
    nzero <- ncells - length(x = xidx)
    csum <- nzero * ifelse(
      test = isTRUE(x = standardize),
      yes = ((0 - mu[i]) / sd[i]) ^ 2,
      no = mu[i] ^ 2
    )
    if (isTRUE(x = standardize)) {
      xidx <- xidx / sd[i]
      xidx[xidx > clip] <- clip
    }
    fsum <- sum(xidx ^ 2) + csum
    fvars[i] <- fsum / (ncells - 1L)
    if (isTRUE(x = verbose)) {
      setTxtProgressBar(pb = pb, value = i / np)
    }
  }
  if (isTRUE(x = verbose)) {
    close(con = pb)
  }
  return(fvars)
}

.SparseMean <- function(data, margin = 1L) {
  margin <- SeuratObject:::.CheckFmargin(fmargin = margin)
  if (margin != .MARGIN(object = data)) {
    data <- t(x = data)
    margin <- .MARGIN(object = data)
  }
  entryname <- .SparseSlots(x = data, type = 'entries')
  p <- slot(object = data, name = .SparseSlots(x = data, type = 'pointers'))
  if (p[1L] == 0) {
    p <- p + 1L
  }
  np <- length(x = p) - 1L
  nobs <- dim(x = data)[-margin]
  means <- vector(mode = 'numeric', length = np)
  for (i in seq_len(length.out = np)) {
    idx <- seq.int(from = p[i], to = p[i + 1L] - 1L)
    means[i] <- sum(slot(object = data, name = entryname)[idx]) / nobs
  }
  return(means)
}

#' @inheritParams stats::loess
#' @param data A matrix
#' @param fmargin Feature margin
#' @param nselect Number of features to select
#' @param clip After standardization values larger than \code{clip} will be set
#' to \code{clip}; default is \code{NULL} which sets this value to the square
#' root of the number of cells
#'
#' @importFrom Matrix rowMeans
#'
#' @keywords internal
#'
#' @noRd
#'
.VST <- function(
  data,
  fmargin = 1L,
  nselect = 2000L,
  span = 0.3,
  clip = NULL,
  verbose = TRUE,
  ...
) {
  fmargin <- SeuratObject:::.CheckFmargin(fmargin = fmargin)
  nfeatures <- dim(x = data)[fmargin]
  # TODO: Support transposed matrices
  # nfeatures <- nrow(x = data)
  if (IsSparse(x = data)) {
    mean.func <- .SparseMean
    var.func <- .SparseFeatureVar
  } else {
    mean.func <- .Mean
    var.func <- .FeatureVar
  }
  hvf.info <- SeuratObject:::EmptyDF(n = nfeatures)
  # hvf.info$mean <- mean.func(data = data, margin = fmargin)
  hvf.info$mean <- rowMeans(x = data)
  hvf.info$variance <- var.func(
    data = data,
    mu = hvf.info$mean,
    fmargin = fmargin,
    verbose = verbose
  )
  hvf.info$variance.expected <- 0L
  not.const <- hvf.info$variance > 0
  fit <- loess(
    formula = log10(x = variance) ~ log10(x = mean),
    data = hvf.info[not.const, , drop = TRUE],
    span = span
  )
  hvf.info$variance.expected[not.const] <- 10 ^ fit$fitted
  hvf.info$variance.standardized <- var.func(
    data = data,
    mu = hvf.info$mean,
    standardize = TRUE,
    sd = sqrt(x = hvf.info$variance.expected),
    clip = clip,
    verbose = verbose
  )
  hvf.info$variable <- FALSE
  hvf.info$rank <- NA
  vs <- hvf.info$variance.standardized
  vs[vs == 0] <- NA
  vf <- head(
    x = order(hvf.info$variance.standardized, decreasing = TRUE),
    n = nselect
  )
  hvf.info$variable[vf] <- TRUE
  hvf.info$rank[vf] <- seq_along(along.with = vf)
  # colnames(x = hvf.info) <- paste0('vst.', colnames(x = hvf.info))
  return(hvf.info)
}

# hvf.methods$vst <- VST

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# S4 Methods
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




################################################################################
################################# SCTransform ##################################
################################################################################


#' @importFrom SeuratObject Cells DefaultLayer DefaultLayer<- Features
#' LayerData LayerData<-
#'
#' @method SCTransform StdAssay
#' @export
#'
SCTransform.StdAssay <- function(
    object,
    layer = 'counts',
    cell.attr = NULL,
    reference.SCT.model = NULL,
    do.correct.umi = TRUE,
    ncells = 5000,
    residual.features = NULL,
    variable.features.n = 3000,
    variable.features.rv.th = 1.3,
    vars.to.regress = NULL,
    do.scale = FALSE,
    do.center = TRUE,
    clip.range = c(-sqrt(x = ncol(x = object) / 30), sqrt(x = ncol(x = object) / 30)),
    conserve.memory = FALSE,
    return.only.var.genes = TRUE,
    seed.use = 1448145,
    verbose = TRUE,
    ...
) {
  olayer <- layer <- unique(x = layer) %||% DefaultLayer(object = object)
  layers <- Layers(object = object, search = layer)
  dataset.names <- gsub(pattern = paste0(layer, "."), replacement = "", x = layers)
  sct.assay.list <- list()
  for (i in seq_along(along.with = layers)) {
    l <- layers[i]
    if (isTRUE(x = verbose)) {
      message("Running SCTransform on layer: ", l)
    }
    counts <- LayerData(
      object = object,
      layer = l,
      features = Features(x = object, layer = l),
      cells = Cells(x = object, layer = l)
      )
    feature.grid <- DelayedArray::rowAutoGrid(x = counts)
    ##TODO: handle this later for boundary conditions
    cells.grid <- DelayedArray::colAutoGrid(x = counts, ncol = 2000)

    ##TODO: handle this later for boundary conditions
    vp <- cells.grid[[1L]]

    # Read a block from a delayed matrix
    sparse <- DelayedArray::is_sparse(x = counts) # TRUE
    block <- DelayedArray::read_block(x = counts, viewport = vp, as.sparse = sparse)

    counts <- as(object = block, Class = 'dgCMatrix')
    cell.attr.object <- cell.attr[colnames(x = counts),]

    if (!identical(rownames(cell.attr.object), colnames(counts))) {
      print(length(setdiff(rownames(cell.attr.object), colnames(counts))))
      print(length(setdiff(colnames(counts),rownames(cell.attr.object))))
      stop("cell attribute row names must match column names of count matrix")
    }
    vst.out <- SCTransform(object = counts,
                           cell.attr = cell.attr.object,
                           reference.SCT.model = reference.SCT.model,
                           do.correct.umi = do.correct.umi,
                           ncells = ncells,
                           residual.features = residual.features,
                           variable.features.n = variable.features.n,
                           variable.features.rv.th = variable.features.rv.th,
                           vars.to.regress = vars.to.regress,
                           do.scale = do.scale,
                           do.center = do.center,
                           clip.range = clip.range,
                           conserve.memory = conserve.memory,
                           return.only.var.genes = return.only.var.genes,
                           seed.use = seed.use,
                           verbose = verbose,
                           ...)

    residual.type <- vst.out[['residual_type']] %||% 'pearson'
    sct.method <- vst.out[['sct.method']]
    # create output assay and put (corrected) umi counts in count slot
    if (do.correct.umi & residual.type == 'pearson') {
      if (verbose) {
        message('Place corrected count matrix in counts slot')
      }
      assay.out <- CreateAssayObject(counts = vst.out$umi_corrected)
      vst.out$umi_corrected <- NULL
    } else {
      # TODO: restore once check.matrix is in SeuratObject
      # assay.out <- CreateAssayObject(counts = umi, check.matrix = FALSE)
      assay.out <- CreateAssayObject(counts = counts)
    }
    # set the variable genes
    VariableFeatures(object = assay.out) <- vst.out$variable_features
    # put log1p transformed counts in data
    assay.out <- SetAssayData(
      object = assay.out,
      slot = 'data',
      new.data = log1p(x = GetAssayData(object = assay.out, slot = 'counts'))
    )
    scale.data <- vst.out$y
    assay.out <- SetAssayData(
      object = assay.out,
      slot = 'scale.data',
      new.data = scale.data
    )
    vst.out$y <- NULL
    # save clip.range into vst model
    vst.out$arguments$sct.clip.range <- clip.range
    vst.out$arguments$sct.method <- sct.method
    Misc(object = assay.out, slot = 'vst.out') <- vst.out
    assay.out <- as(object = assay.out, Class = "SCTAssay")

    sct.assay.list[[dataset.names[i]]] <- assay.out
  }

  # Return array by merging everythin
  if (length(x = sct.assay.list)>1){
    merged.assay <- merge.SCTAssay(x = sct.assay.list[[1]], y = sct.assay.list[2:length(sct.assay.list)])
    # set the names of SCTmodels to be layer names
    models <- slot(object = merged.assay, name="SCTModel.list")
    names(models) <- names(x = sct.assay.list)
    slot(object = merged.assay, name="SCTModel.list") <- models
  } else {
    return (sct.assay.list[[1]])
  }
  gc(verbose = FALSE)
  return(merged.assay)
}

#' @importFrom SeuratObject DefaultAssay
#'
#' @method SCTransform Seurat5
#' @export
#'
SCTransform.Seurat5 <- function(
    object,
    assay = NULL,
    reference.SCT.model = NULL,
    do.correct.umi = TRUE,
    ncells = 5000,
    residual.features = NULL,
    variable.features.n = 3000,
    variable.features.rv.th = 1.3,
    vars.to.regress = NULL,
    do.scale = FALSE,
    do.center = TRUE,
    clip.range = c(-sqrt(x = ncol(x = object) / 30), sqrt(x = ncol(x = object) / 30)),
    conserve.memory = FALSE,
    return.only.var.genes = TRUE,
    seed.use = 1448145,
    save.data = 'data',
    save.scaledata = 'scale.data',
    verbose = TRUE,
    ...
) {
  assay <- assay[1L] %||% DefaultAssay(object = object)
  assay <- match.arg(arg = assay, choices = Assays(object = object))
  cell.attr.list <- slot(object = object, name = 'meta.data')

  object[[assay]] <- SCTransform(object = object[[assay]],
                                 cell.attr.list = cell.attr.list,
                                 reference.SCT.model = reference.SCT.model,
                                 do.correct.umi = do.correct.umi,
                                 ncells = ncells,
                                 residual.features = residual.features,
                                 variable.features.n = variable.features.n,
                                 variable.features.rv.th = variable.features.rv.th,
                                 vars.to.regress = vars.to.regress,
                                 do.scale = do.scale,
                                 do.center = do.center,
                                 clip.range = clip.range,
                                 conserve.memory = conserve.memory,
                                 return.only.var.genes = return.only.var.genes,
                                 seed.use = seed.use,
                                 verbose = verbose,
                                 ...)
  return(object)
}


#' Calculate pearson residuals of features not in the scale.data
#'
#' This function calls sctransform::get_residuals.
#'
#' @param object A seurat object
#' @param features Name of features to add into the scale.data
#' @param assay Name of the assay of the seurat object generated by SCTransform
#' @param layer Name (prefix) of the layer to pull counts from
#' @param umi.assay Name of the assay of the seurat object containing UMI matrix
#' and the default is RNA
#' @param clip.range Numeric of length two specifying the min and max values the
#' Pearson residual will be clipped to
#' @param replace.value Recalculate residuals for all features, even if they are
#' already present. Useful if you want to change the clip.range.
#' @param na.rm For features where there is no feature model stored, return NA
#' for residual value in scale.data when na.rm = FALSE. When na.rm is TRUE, only
#' return residuals for features with a model stored for all cells.
#' @param verbose Whether to print messages and progress bars
#'
#' @return Returns a Seurat object containing Pearson residuals of added
#' features in its scale.data
#'
#' @importFrom sctransform get_residuals
#' @importFrom matrixStats rowAnyNAs
#'
#' @export
#' @concept preprocessing
#'
#' @seealso \code{\link[sctransform]{get_residuals}}
#'
#' @examples
#' data("pbmc_small")
#' pbmc_small <- SCTransform(object = pbmc_small, variable.features.n = 20)
#' pbmc_small <- GetResidual(object = pbmc_small, features = c('MS4A1', 'TCL1A'))
#'
GetResidual.V5 <- function(
    object,
    features,
    assay = NULL,
    umi.assay = "RNA",
    layer = "counts",
    clip.range = NULL,
    replace.value = FALSE,
    na.rm = TRUE,
    verbose = TRUE
) {
  assay <- assay %||% DefaultAssay(object = object)
  if (IsSCT(assay = object[[assay]])) {
    object[[assay]] <- as(object[[assay]], 'SCTAssay')
  }
  if (!inherits(x = object[[assay]], what = "SCTAssay")) {
    stop(assay, " assay was not generated by SCTransform")
  }
  sct.models <- levels(x = object[[assay]])
  if (length(x = sct.models) == 0) {
    warning("SCT model not present in assay", call. = FALSE, immediate. = TRUE)
    return(object)
  }
  possible.features <- unique(x = unlist(x = lapply(X = sct.models, FUN = function(x) {
    rownames(x = SCTResults(object = object[[assay]], slot = "feature.attributes", model = x))
  }
  )))
  bad.features <- setdiff(x = features, y = possible.features)
  if (length(x = bad.features) > 0) {
    warning("The following requested features are not present in any models: ",
            paste(bad.features, collapse = ", "), call. = FALSE)
    features <- intersect(x = features, y = possible.features)
  }
  features.orig <- features
  if (na.rm) {
    # only compute residuals when feature model info is present in all
    features <- names(x = which(x = table(unlist(x = lapply(
      X = sct.models,
      FUN = function(x) {
        rownames(x = SCTResults(object = object[[assay]], slot = "feature.attributes", model = x))
      }
    ))) == length(x = sct.models)))
    if (length(x = features) == 0) {
      return(object)
    }
  }
  features <- intersect(x = features.orig, y = features)
  if (length(x = sct.models) > 1 & verbose) {
    message("This SCTAssay contains multiple SCT models. Computing residuals for cells using")
  }
  cat(sct.models)

  new.residuals <- lapply(
    X = sct.models,
    FUN = function(x) {
      GetResidualSCTModel.V5(
        object = object,
        umi.assay = umi.assay,
        assay = assay,
        layer = layer,
        SCTModel = x,
        new_features = features,
        replace.value = replace.value,
        clip.range = clip.range,
        verbose = verbose
      )
    }
  )
  existing.data <- GetAssayData(object = object, slot = 'scale.data', assay = assay)
  all.features <- union(x = rownames(x = existing.data), y = features)
  new.scale <- matrix(
    data = NA,
    nrow = length(x = all.features),
    ncol = ncol(x = object),
    dimnames = list(all.features, Cells(x = object))
  )
  if (nrow(x = existing.data) > 0){
    new.scale[1:nrow(x = existing.data), ] <- existing.data
  }
  if (length(x = new.residuals) == 1 & is.list(x = new.residuals)) {
    new.residuals <- new.residuals[[1]]
  } else {
    new.residuals <- Reduce(cbind, new.residuals)
  }
  new.scale[rownames(x = new.residuals), colnames(x = new.residuals)] <- new.residuals
  if (na.rm) {
    new.scale <- new.scale[!rowAnyNAs(x = new.scale), ]
  }
  object <- SetAssayData(
    object = object,
    assay = assay,
    slot = "scale.data",
    new.data = new.scale
  )
  if (any(!features.orig %in% rownames(x = new.scale))) {
    bad.features <- features.orig[which(!features.orig %in% rownames(x = new.scale))]
    warning("Residuals not computed for the following requested features: ",
            paste(bad.features, collapse = ", "), call. = FALSE)
  }
  return(object)
}


# Calculate pearson residuals of features not in the scale.data
# This function is the secondary function under GetResidual
#
# @param object A seurat object
# @param features Name of features to add into the scale.data
# @param assay Name of the assay of the seurat object generated by SCTransform
# @param vst_out The SCT parameter list
# @param clip.range Numeric of length two specifying the min and max values the Pearson residual
# will be clipped to
# Useful if you want to change the clip.range.
# @param verbose Whether to print messages and progress bars
#
# @return Returns a matrix containing not-centered pearson residuals of added features
#
#' @importFrom sctransform get_residuals
#
GetResidualSCTModel.V5 <- function(
    object,
    assay,
    umi.assay,
    layer,
    SCTModel,
    new_features,
    clip.range,
    replace.value,
    verbose
) {
  clip.range <- clip.range %||% SCTResults(object = object[[assay]], slot = "clips", model = SCTModel)$sct
  model.features <- rownames(x = SCTResults(object = object[[assay]], slot = "feature.attributes", model = SCTModel))
  #umi.assay <- SCTResults(object = object[[assay]], slot = "umi.assay", model = SCTModel)
  model.cells <- Cells(x = slot(object = object[[assay]], name = "SCTModel.list")[[SCTModel]])
  sct.method <-  SCTResults(object = object[[assay]], slot = "arguments", model = SCTModel)$sct.method %||% "default"
  scale.data.cells <- colnames(x = GetAssayData(object = object, assay = assay, slot = "scale.data"))
  if (length(x = setdiff(x = model.cells, y =  scale.data.cells)) == 0) {
    existing_features <- names(x = which(x = ! apply(
      X = GetAssayData(object = object, assay = assay, slot = "scale.data")[, model.cells],
      MARGIN = 1,
      FUN = anyNA)
    ))
  } else {
    existing_features <- character()
  }
  if (replace.value) {
    features_to_compute <- new_features
  } else {
    features_to_compute <- setdiff(x = new_features, y = existing_features)
  }
  if (sct.method == "reference.model") {
    if (verbose) {
      message("sct.model ", SCTModel, " is from reference, so no residuals will be recalculated")
    }
    features_to_compute <- character()
  }
  if (!umi.assay %in% Assays(object = object)) {
    warning("The umi assay (", umi.assay, ") is not present in the object. ",
            "Cannot compute additional residuals.", call. = FALSE, immediate. = TRUE)
    return(NULL)
  }
  diff_features <- setdiff(x = features_to_compute, y = model.features)
  intersect_features <- intersect(x = features_to_compute, y = model.features)
  if (length(x = diff_features) == 0) {
    #umi <- GetAssayData(object = object, assay = umi.assay, slot = "counts" )[features_to_compute, model.cells, drop = FALSE]
    print(object[[umi.assay]])
    print(layer)
    #layer <- LayerData(object = object[[umi.assay]], layer = layer)

    #olayer <- layer <- unique(x = layer) %||% DefaultLayer(object = object[[umi.assay]])
    layers <- Layers(object = object[[umi.assay]], search = layer)
    dataset.names <- gsub(pattern = paste0(layer, "."), replacement = "", x = layers)
    vst_out <- SCTModel_to_vst(SCTModel = slot(object = object[[assay]], name = "SCTModel.list")[[SCTModel]])

    for (i in seq_along(along.with = layers)) {
      l <- layers[i]
      if (isTRUE(x = verbose)) {
        message("Running SCTransform on layer: ", l)
      }
      counts <- LayerData(
        object = object,
        layer = l,
        features = Features(x = object, layer = l),
        cells = Cells(x = object, layer = l)
      )
      feature.grid <- DelayedArray::rowAutoGrid(x = counts)
      ##TODO: handle this later for boundary conditions
      cells.grid <- DelayedArray::colAutoGrid(x = counts, ncol = 2000)


      new_residuals <- list()
      for (i in seq_len(length.out = length(x = cells.grid))) {
        vp <- cells.grid[[i]]
        block <- DelayedArray::read_block(x = counts, viewport = vp, as.sparse = TRUE)
        umi <- as(object = block, Class = "dgCMatrix")
        if (verbose) {
          message("sct.model: ", SCTModel)
        }
        new_residual <- get_residuals(
          vst_out = vst_out,
          umi = umi,
          residual_type = "pearson",
          res_clip_range = c(clip.min, clip.max),
          verbosity = as.numeric(x = verbose) * 2
        )
        new_residual <- as.matrix(x = new_residual)
        # centered data
        new_residuals[[i]] <- new_residual
      }

      new_residual <- do.call(what = cbind, args = new_residuals)
      new_residual <- new_residual - rowMeans(x = new_residual)

    }

  } else {
    warning(
      "In the SCTModel ", SCTModel, ", the following ", length(x = diff_features),
      " features do not exist in the counts slot: ", paste(diff_features, collapse = ", ")
    )
    if (length(x = intersect_features) == 0) {
      return(matrix(
        data = NA,
        nrow = length(x = features_to_compute),
        ncol = length(x = model.cells),
        dimnames = list(features_to_compute, model.cells)
      ))
    }
    umi <- GetAssayData(object = object, assay = umi.assay, slot = "counts")[intersect_features, model.cells, drop = FALSE]
  }
  clip.max <- max(clip.range)
  clip.min <- min(clip.range)
  if (nrow(x = umi) > 0) {
    vst_out <- SCTModel_to_vst(SCTModel = slot(object = object[[assay]], name = "SCTModel.list")[[SCTModel]])
    if (verbose) {
      message("sct.model: ", SCTModel)
    }
    new_residual <- get_residuals(
      vst_out = vst_out,
      umi = umi,
      residual_type = "pearson",
      res_clip_range = c(clip.min, clip.max),
      verbosity = as.numeric(x = verbose) * 2
    )
    new_residual <- as.matrix(x = new_residual)
    # centered data
    new_residual <- new_residual - rowMeans(x = new_residual)
  } else {
    new_residual <- matrix(data = NA, nrow = 0, ncol = length(x = model.cells), dimnames = list(c(), model.cells))
  }
  old.features <- setdiff(x = new_features, y = features_to_compute)
  if (length(x = old.features) > 0) {
    old_residuals <- GetAssayData(object = object[[assay]], slot = "scale.data")[old.features, model.cells, drop = FALSE]
    new_residual <- rbind(new_residual, old_residuals)[new_features, ]
  }
  return(new_residual)
}

