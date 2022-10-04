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
  # fpath <- slot(object = slot(object = data, name = 'seed'), name = 'filepath')
  fpath <- DelayedArray::path(object = data)
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
  verbose = TRUE,
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

#' @method LogNormalize TileDBMatrix
#' @export
#'
LogNormalize.TileDBMatrix <- function(
  data,
  scale.factor = 1e4,
  margin = 2L,
  return.seed = TRUE,
  verbose= TRUE,
  layer = 'data',
  ...
) {
  check_installed(
    pkg = 'TileDBArray',
    reason = 'for working with TileDB matrices'
  )
  odir <- c(
    dirname(path = DelayedArray::path(object = data)),
    getwd(),
    tempdir(check = TRUE)
  )
  # file.access returns 0 (FALSE) for true and -1 (TRUE) for false
  idx <- which(x = !file.access(names = odir, mode = 2L))[1L]
  if (rlang::is_na(x = odir)) {
    abort(message = "Unable to find a directory to write normalized TileDB matrix")
  }
  out <- file.path(odir[idx], layer)
  if (!file.access(names = out, mode = 0L)) {
    warn(message = paste(sQuote(x = out), "exists, overwriting"))
    unlink(x = out, recursive = TRUE, force = TRUE)
  }
  sink <- TileDBArray::TileDBRealizationSink(
    dim = dim(x = data),
    dimnames = dimnames(x = data),
    type = BiocGenerics::type(x = data),
    path = out,
    attr = layer,
    sparse = DelayedArray::is_sparse(x = data)
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
