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
  var.gene.ouput <- method(
    data = object,
    nselect = nselect,
    verbose = verbose,
    ...
  )
  rownames(x = var.gene.ouput) <- rownames(x = object)
  return(var.gene.ouput)
}


#' @importFrom SeuratObject DefaultLayer Features Key Layers
#'
#' @method FindVariableFeatures StdAssay
#' @export
#'
FindVariableFeatures.StdAssay <- function(
  object,
  method = NULL,
  nselect = 2000L,
  layer = NULL,
  span = 0.3,
  clip = NULL,
  key = NULL,
  verbose = TRUE,
  selection.method = 'vst',
  ...
) {
  if (selection.method == 'vst') {
    layer <- layer%||%'counts'
    method <- VST
    key <- 'vst'
  } else if (selection.method %in% c('mean.var.plot', 'mvp')) {
    layer <- layer%||%'data'
    method <- MVP
    key <- 'mvp'
  } else if (selection.method %in% c('dispersion', 'disp')) {
    layer <- layer%||%'data'
    method <- DISP
    key <- 'disp'
  } else if (is.null(x = method) || is.null(x = layer)){
    stop('Custome functions and layers are both required')
  } else {
    key <- NULL
  }
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
    hvf.function <- if (inherits(x = data, what = 'V3Matrix')) {
      FindVariableFeatures.default
    } else {
      FindVariableFeatures
    }
    hvf.info <- hvf.function(
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
    object[colnames(x = hvf.info)] <- hvf.info
  }
  var.name <- paste(
    'vf',
    key,
    layer[i],
    'variable',
    sep = '_'
  )
  VariableFeatures(object = object) <- rownames(hvf.info)[hvf.info[,var.name]]
  return(object)
}

#' @param layer Layer in the Assay5 to pull data from
#' @param features If provided, only compute on given features. Otherwise,
#' compute for all features.
#' @param nfeatures Number of features to mark as the top spatially variable.
#'
#' @method FindSpatiallyVariableFeatures StdAssay
#' @rdname FindSpatiallyVariableFeatures
#' @concept preprocessing
#' @concept spatial
#' @export
#'
FindSpatiallyVariableFeatures.StdAssay <- function(
  object,
  layer = "scale.data",
  spatial.location,
  selection.method = c('markvariogram', 'moransi'),
  features = NULL,
  r.metric = 5,
  x.cuts = NULL,
  y.cuts = NULL,
  nfeatures = nfeatures,
  verbose = TRUE,
  ...
) {
  features <- features %||% rownames(x = object)
  if (selection.method == "markvariogram" && "markvariogram" %in% names(x = Misc(object = object))) {
    features.computed <- names(x = Misc(object = object, slot = "markvariogram"))
    features <- features[! features %in% features.computed]
  }
  data <- GetAssayData(object = object, layer = layer)
  data <- as.matrix(x = data[features, ])
  data <- data[RowVar(x = data) > 0, ]
  if (nrow(x = data) != 0) {
    svf.info <- FindSpatiallyVariableFeatures(
      object = data,
      spatial.location = spatial.location,
      selection.method = selection.method,
      r.metric = r.metric,
      x.cuts = x.cuts,
      y.cuts = y.cuts,
      verbose = verbose,
      ...
    )
  } else {
    svf.info <- c()
  }
  if (selection.method == "markvariogram") {
    if ("markvariogram" %in% names(x = Misc(object = object))) {
      svf.info <- c(svf.info, Misc(object = object, slot = "markvariogram"))
    }
    suppressWarnings(expr = Misc(object = object, slot = "markvariogram") <- svf.info)
    svf.info <- ComputeRMetric(mv = svf.info, r.metric)
    svf.info <- svf.info[order(svf.info[, 1]), , drop = FALSE]
  }
  if (selection.method == "moransi") {
    colnames(x = svf.info) <- paste0("MoransI_", colnames(x = svf.info))
    svf.info <- svf.info[order(svf.info[, 2], -abs(svf.info[, 1])), , drop = FALSE]
  }
  var.name <- paste0(selection.method, ".spatially.variable")
  var.name.rank <- paste0(var.name, ".rank")
  svf.info[[var.name]] <- FALSE
  svf.info[[var.name]][1:(min(nrow(x = svf.info), nfeatures))] <- TRUE
  svf.info[[var.name.rank]] <- 1:nrow(x = svf.info)
  object[names(x = svf.info)] <- svf.info
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
  } else if (inherits(x = sink, what = 'arrayRealizationSink')) {
    # arrayRealizationSinks do not support write_block with rowAutoGrid or colAutoGrid
    # Because of course they don't
    abort(message = "Array RealizationSinks are not supported due to issues with {DelayedArray}")
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
    if (isTRUE(x = sparse)) {
      x <- as(object = x, Class = 'dgCMatrix')
    }
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


#' @method LogNormalize IterableMatrix
#' @export
#'
LogNormalize.IterableMatrix <- function(
    data,
    scale.factor = 1e4,
    margin = 2L,
    verbose = TRUE,
    ...
) {
  data <- BPCells::t(BPCells::t(data) / colSums(data))
  # Log normalization
  data <- log1p(data * scale.factor)
  return(data)
}
#' @method LogNormalize TileDBMatrix
#' @export
#'
LogNormalize.TileDBMatrix <- function(
  data,
  scale.factor = 1e4,
  margin = 2L,
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
  if (rlang::is_na(x = idx)) {
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
  return(NextMethod(
    generic = 'LogNormalize',
    object = data,
    scale.factor = scale.factor,
    margin = margin,
    verbose = verbose,
    sink = sink,
    ...
  ))
}

#' @method LogNormalize TENxMatrix
#' @export
#'
LogNormalize.TENxMatrix <- function(
  data,
  scale.factor = 1e4,
  margin = 2L,
  verbose= TRUE,
  layer = 'data',
  ...
) {
  check_installed(pkg = 'HDF5Array', reason = 'for working with HDF5 matrices')
  fpath <- DelayedArray::path(object = data)
  if (.DelayedH5DExists(object = data, path = layer)) {
    rhdf5::h5delete(file = fpath, name = layer)
  }
  sink <- HDF5Array::TENxRealizationSink(
    dim = dim(x = data),
    dimnames = dimnames(x = data),
    filepath = fpath,
    group = layer
  )
  return(NextMethod(
    generic = 'LogNormalize',
    object = data,
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
  normalization.method = c('LogNormalize', 'CLR', 'RC'),
  scale.factor = 1e4,
  cmargin = 2L,
  margin = 1L,
  verbose = TRUE,
  ...
) {
  normalization.method <- normalization.method[1L]
  normalization.method <- match.arg(arg = normalization.method)
  # TODO: enable parallelization via future
  normalized <- switch(
    EXPR = normalization.method,
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
    },
    'CLR' = {
      if (!inherits(x = object, what = 'dgCMatrix') &&
          !inherits(x = object, what = 'matrix')) {
        stop('CLR normalization only supports for dense and dgCMatrix')
      }
      CustomNormalize(
        data = object,
        custom_function = function(x) {
          return(log1p(x = x/(exp(x = sum(log1p(x = x[x > 0]), na.rm = TRUE)/length(x = x)))))
        },
        margin = margin,
        verbose = verbose
      )
    },
    'RC' = {
      if (!inherits(x = object, what = 'dgCMatrix') &&
          !inherits(x = object, what = 'matrix')) {
        stop('RC normalization only supports for dense and dgCMatrix')
      }
      RelativeCounts(data = object,
                     scale.factor = scale.factor,
                     verbose = verbose)
    }
  )
  return(normalized)
}

.DelayedH5DExists <- function(object, path) {
  check_installed(pkg = 'HDF5Array', reason = 'for working with HDF5 files')
  if (!inherits(x = object, what = c('HDF5Array', 'H5ADMatrix', 'TENxMatrix'))) {
    abort(message = "'object' must be an HDF5Array or H5ADMatrix")
  }
  on.exit(expr = rhdf5::h5closeAll(), add = TRUE)
  fpath <- DelayedArray::path(object = object)
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
  normalization.method = 'LogNormalize',
  scale.factor = 1e4,
  margin = 1L,
  layer = 'counts',
  save = 'data',
  verbose = TRUE,
  ...
) {
  olayer <- layer <- unique(x = layer)
  layer <- Layers(object = object, search = layer)
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
      normalization.method = normalization.method,
      scale.factor = scale.factor,
      margin = margin,
      verbose = verbose,
      layer = save,
      ...
    )
  }
  gc(verbose = FALSE)
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
  layer = 'data',
  vars.to.regress = NULL,
  latent.data = NULL,
  by.layer = FALSE,
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
  olayer <- layer <- unique(x = layer)
  layer <- Layers(object = object, search = layer)
  if (isTRUE(x = use.umi)) {
    inform(
      message = "'use.umi' is TRUE, please make sure 'layer' specifies raw counts"
    )
  }
  features <- features %||% VariableFeatures(object = object)
  if (!length(x = features)) {
    features <- Features(x = object, layer = layer)
  }
  if (isTRUE(x = by.layer)) {
    if (length(x = save) != length(x = layer)) {
      save <- make.unique(names = gsub(
        pattern = olayer,
        replacement = save,
        x = layer
      ))
    }
    for (i in seq_along(along.with = layer)) {
      lyr <- layer[i]
      if (isTRUE(x = verbose)) {
        inform(message = paste("Scaling data for layer", sQuote(x = lyr)))
      }
      LayerData(object = object, layer = save[i], ...) <- ScaleData(
        object = LayerData(
          object = object,
          layer = lyr,
          features = features,
          fast = NA
        ),
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
    }
  } else {
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
  }
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

#' @rdname VST
#' @method VST IterableMatrix
#' @export
#'
VST.IterableMatrix <-function(
    data,
    nselect = 2000L,
    span = 0.3,
    clip = NULL,
    verbose = TRUE,
    ...
) {
  nfeatures <- nrow(x = data)
  hvf.info <- SeuratObject:::EmptyDF(n = nfeatures)
  hvf.stats <- BPCells::matrix_stats(
    matrix = data,
    row_stats = 'variance')$row_stats
  # Calculate feature means
  hvf.info$mean <- hvf.stats['mean',]
  # Calculate feature variance
  hvf.info$variance <- hvf.stats['variance',]
  hvf.info$variance.expected <- 0L
  not.const <- hvf.info$variance > 0
  fit <- loess(
    formula = log10(x = variance) ~ log10(x = mean),
    data = hvf.info[not.const, , drop = TRUE],
    span = span
  )
  hvf.info$variance.expected[not.const] <- 10 ^ fit$fitted
  feature.mean <- hvf.info$mean
  feature.sd <-  sqrt(x = hvf.info$variance.expected)
  standard.max <- clip %||% sqrt(x = ncol(x = data))
  feature.mean[feature.mean == 0] <- 0.1
  data <- BPCells::min_by_row(mat = data, vals = standard.max*feature.sd + feature.mean)
  data.standard <- (data - feature.mean)/feature.sd
  hvf.info$variance.standardized <- BPCells::matrix_stats(
    matrix = data.standard,
    row_stats = 'variance'
    )$row_stats['variance',]
  # Set variable features
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

#' @method VST DelayedMatrix
#' @export
#'
VST.DelayedMatrix <- function(
  data,
  margin = 2L,
  nselect = 2000L,
  span = 0.3,
  clip = NULL,
  verbose = TRUE,
  block.size = 1e8,
  ...
) {
  check_installed(
    pkg = 'DelayedArray',
    reason = 'for working with delayed matrices'
  )
  if (!margin %in% c(1L, 2L)) {
    abort(message = "'margin' must be 1 or 2")
  }
  nfeatures <- dim(x = data)[-margin]
  ncells <- dim(x = data)[margin]
  hvf.info <- SeuratObject::EmptyDF(n = nfeatures)
  hvf.info$mean <- RowMeanDelayedAssay(x = data, block.size = block.size)
  # Calculate feature variance
  hvf.info$variance <- RowVarDelayedAssay(x = data, block.size = block.size)
  hvf.info$variance.expected <- 0L
  not.const <- hvf.info$variance > 0
  fit <- loess(
    formula = log10(x = variance) ~ log10(x = mean),
    data = hvf.info[not.const, , drop = TRUE],
    span = span
  )
  hvf.info$variance.expected[not.const] <- 10 ^ fit$fitted

  suppressMessages(setAutoBlockSize(size = block.size))
  grid <- if (margin == 1L) {
    DelayedArray::rowAutoGrid(x = data)
  } else {
    DelayedArray::colAutoGrid(x = data)
  }
  sparse <- DelayedArray::is_sparse(x = data)
  if (sparse) {
    sweep.func <- SweepSparse
    rowsum.func <- RowSumSparse
  } else {
    sweep.func <- sweep
    rowsum.func <- rowSums2
  }
  var_stand.list <- list()
  for (i in seq_len(length.out = length(x = grid))) {
    vp <- grid[[i]]
    block <- DelayedArray::read_block(x = data, viewport = vp, as.sparse = sparse)
    block <- as(object = block, Class = 'dgCMatrix')
    block.stat <- SparseRowVarStd(mat = block,
                                  mu = hvf.info$mean,
                                  sd = sqrt(hvf.info$variance.expected),
                                  vmax =  clip %||% sqrt(x = ncol(x = data)),
                                  display_progress = FALSE)

    var_stand.list[[i]] <- block.stat * (ncol(block) - 1)
  }
  hvf.info$variance.standardized <- Reduce(f = '+', x = var_stand.list)/
    (ncol(data) - 1)
  # Set variable features
  hvf.info$variable <- FALSE
  hvf.info$rank <- NA
  vf <- head(
    x = order(hvf.info$variance.standardized, decreasing = TRUE),
    n = nselect
  )
  hvf.info$variable[vf] <- TRUE
  hvf.info$rank[vf] <- seq_along(along.with = vf)
 rownames(hvf.info) <- rownames(data)

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


#' Calculate dispersion of features
#'
#'
CalcDispersion <- function(
  object,
  mean.function = FastExpMean,
  dispersion.function = FastLogVMR,
  num.bin = 20,
  binning.method = "equal_width",
  verbose = TRUE,
  ...
) {
  if (!inherits(x = object, what = c('dgCMatrix', 'matrix'))) {
    stop('mean.var.plot and dispersion methods only support dense and sparse matrix input')
  }
  if (inherits(x = object, what =  'matrix')) {
    object <- as.sparse(x = object)
  }
  feature.mean <- mean.function(object, verbose)
  feature.dispersion <- dispersion.function(object, verbose)

  names(x = feature.mean) <- names(x = feature.dispersion) <- rownames(x = object)
  feature.dispersion[is.na(x = feature.dispersion)] <- 0
  feature.mean[is.na(x = feature.mean)] <- 0
  data.x.breaks <- switch(
    EXPR = binning.method,
    'equal_width' = num.bin,
    'equal_frequency' = c(
      quantile(
        x = feature.mean[feature.mean > 0],
        probs = seq.int(from = 0, to = 1, length.out = num.bin)
      )
    ),
    stop("Unknown binning method: ", binning.method)
  )
  data.x.bin <- cut(x = feature.mean, breaks = data.x.breaks,
                    include.lowest = TRUE)
  names(x = data.x.bin) <- names(x = feature.mean)
  mean.y <- tapply(X = feature.dispersion, INDEX = data.x.bin, FUN = mean)
  sd.y <- tapply(X = feature.dispersion, INDEX = data.x.bin, FUN = sd)
  feature.dispersion.scaled <- (feature.dispersion - mean.y[as.numeric(x = data.x.bin)]) /
    sd.y[as.numeric(x = data.x.bin)]
  names(x = feature.dispersion.scaled) <- names(x = feature.mean)
  hvf.info <- data.frame(feature.mean, feature.dispersion, feature.dispersion.scaled)
  rownames(x = hvf.info) <- rownames(x = object)
  colnames(x = hvf.info) <- paste0('mvp.', c('mean', 'dispersion', 'dispersion.scaled'))
  return(hvf.info)
}


#' @importFrom SeuratObject .CalcN
#'
CalcN <- function(object) {
  return(.CalcN(object))
}

#' @method .CalcN IterableMatrix
#' @export
#'
.CalcN.IterableMatrix <- function(object) {
  col_stat <- BPCells::matrix_stats(matrix = object, col_stats = 'mean')$col_stats
  return(list(
    nCount = round(col_stat['mean',] *nrow(object)),
    nFeature = col_stat['nonzero',]
  ))
}

#' Find variable features based on dispersion
#'
DISP <- function(
  data,
  nselect = 2000L,
  verbose = TRUE,
  ...
) {
  hvf.info <- CalcDispersion(object = data, verbose = verbose, ...)
  hvf.info$variable <- FALSE
  hvf.info$rank <- NA
  vf <- head(
    x = order(hvf.info$mvp.dispersion, decreasing = TRUE),
    n = nselect
  )
  hvf.info$variable[vf] <- TRUE
  hvf.info$rank[vf] <- seq_along(along.with = vf)
  return(hvf.info)
}

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
  if (!is.null(reference.SCT.model)){
    do.correct.umi <- FALSE
    do.center <- FALSE
  }
  olayer <- layer <- unique(x = layer)
  layers <- Layers(object = object, search = layer)
  dataset.names <- gsub(pattern = paste0(layer, "."), replacement = "", x = layers)
  sct.assay.list <- list()
  for (dataset.index in seq_along(along.with = layers)) {
    l <- layers[dataset.index]
    if (isTRUE(x = verbose)) {
      message("Running SCTransform on layer: ", l)
    }
    all_cells <-  Cells(x = object, layer = l)
    all_features <- Features(x = object, layer = l)
    counts <- LayerData(
      object = object,
      layer = l,
      features = all_features,
      cells = all_cells
    )
    sparse <- DelayedArray::is_sparse(x = counts)
    ## Sample  cells
    cells.grid <- DelayedArray::colAutoGrid(x = counts, ncol = min(ncells, ncol(counts)))
    # if there is no reference model we randomly select a subset of cells
    # TODO: randomize this set of cells
    variable.feature.list <- list()
    GetSCT.Chunked <- function(vp, reference.SCT.model = NULL, do.correct.umi = TRUE){
      # counts here is global
      block <- DelayedArray::read_block(x = counts,
                                        viewport = vp,
                                        as.sparse = sparse)
      counts.chunk <- as(object = block, Class = 'dgCMatrix')
      cell.attr.object <- cell.attr[colnames(x = counts.chunk),, drop=FALSE]

      if (!identical(rownames(cell.attr.object), colnames(counts.chunk))) {
        stop("cell attribute row names must match column names of count matrix")
      }
      vst.out <- SCTransform(object = counts.chunk,
                             cell.attr = cell.attr.object,
                             reference.SCT.model = reference.SCT.model,
                             do.correct.umi = do.correct.umi,
                             ncells = ncells,
                             residual.features = residual.features,
                             variable.features.n = variable.features.n,
                             variable.features.rv.th = variable.features.rv.th,
                             vars.to.regress = vars.to.regress,
                             do.scale = FALSE,
                             do.center = FALSE,
                             clip.range = clip.range,
                             conserve.memory = conserve.memory,
                             return.only.var.genes = return.only.var.genes,
                             seed.use = seed.use,
                             verbose = FALSE,
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
        assay.out <- CreateAssayObject(counts = counts.chunk)
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
      # does not like character(0) keys being merged
      return (assay.out)
    }
    local.reference.SCT.model <- NULL
    if (is.null(reference.SCT.model)){
      # No reference model so just select the some block of cells
      set.seed(seed = seed.use)
      selected.block <-  sample(x = seq.int(from = 1, to = length(cells.grid)), size = 1)
      if (verbose){
        message("Using block ", selected.block, " from ", dataset.names[[dataset.index]], " to learn model.")
      }
      vp <- cells.grid[[selected.block]]

      do.correct.umi.chunk <- FALSE
      # correct umi if only single chunk
      if (length(x = cells.grid) == 1) {
        do.correct.umi.chunk <- TRUE
      }
      assay.out <- GetSCT.Chunked(vp = vp, do.correct.umi = do.correct.umi.chunk)
      local.reference.SCT.model <- assay.out@SCTModel.list[[1]]
      variable.features <- VariableFeatures(assay.out)
      # once we have the model, just calculate residuals for all
      # cells
      vst_out.reference <- SCTModel_to_vst(SCTModel = local.reference.SCT.model)
      min_var <- vst_out.reference$arguments$min_variance
      if (min_var == "umi_median"){
        block <- DelayedArray::read_block(x = counts,
                                          viewport = vp,
                                          as.sparse = TRUE)

        counts.x <- as(object = block, Class = 'dgCMatrix')
        min_var <- (median(counts.x@x)/5)^2
      }
      res_clip_range <-  vst_out.reference$arguments$res_clip_range
      residuals <- list()
      corrected_counts <- list()
      cell_attrs <- list()
      if (length(x = cells.grid) == 1){
        merged.assay <- assay.out
        corrected_counts[[1]] <- GetAssayData(object = assay.out, slot="data")
        residuals[[1]] <- GetAssayData(object = assay.out, slot="scale.data")
        cell_attrs[[1]] <- vst_out.reference$cell_attr
        sct.assay.list[[dataset.names[dataset.index]]] <- assay.out
      } else {
        # iterate over chunks to get residuals
        for (i in seq_len(length.out = length(x = cells.grid))) {
          vp <- cells.grid[[i]]
          if (verbose){
            message("Getting residuals for block ", i, "(of ", length(cells.grid), ") for ", dataset.names[[dataset.index]], " dataset")
          }
          block <- DelayedArray::read_block(x = counts,
                                            viewport = vp,
                                            as.sparse = TRUE)

          counts.vp <- as(object = block, Class = 'dgCMatrix')
          cell.attr.object <- cell.attr[colnames(x = counts.vp),, drop=FALSE]
          vst_out <- vst_out.reference
          cell_attr <- data.frame(
            umi = colSums(counts.vp),
            log_umi = log10(x = colSums(counts.vp))
          )
          rownames(cell_attr) <- colnames(counts.vp)
          vst_out$cell_attr <- cell_attr
          vst_out$gene_attr <- vst_out$gene_attr[variable.features,]
          if (return.only.var.genes){
            new_residual <- get_residuals(
              vst_out = vst_out,
              umi = counts.vp[variable.features,],
              residual_type = "pearson",
              min_variance = min_var,
              res_clip_range = res_clip_range,
              verbosity =  FALSE#as.numeric(x = verbose) * 2
              )
          } else {
            new_residual <- get_residuals(
              vst_out = vst_out,
              umi = counts.vp[all.features,],
              residual_type = "pearson",
              min_variance = min_var,
              res_clip_range = res_clip_range,
              verbosity =  FALSE#as.numeric(x = verbose) * 2
            )
          }
          vst_out$y <- new_residual
          corrected_counts[[i]] <- correct_counts(
            x = vst_out,
            umi = counts.vp[all_features,],
            verbosity = FALSE# as.numeric(x = verbose) * 2
          )
          residuals[[i]] <- new_residual
          cell_attrs[[i]] <- cell_attr
        }
        new.residuals <- Reduce(cbind, residuals)

        corrected_counts <- Reduce(cbind, corrected_counts)
        cell_attrs <- Reduce(rbind, cell_attrs)

        vst_out.reference$cell_attr <- cell_attrs[colnames(new.residuals),]
        SCTModel.list <- PrepVSTResults(
          vst.res = vst_out.reference,
          cell.names = all_cells
        )
        SCTModel.list <- list(model1 = SCTModel.list)

        # scale data here as do.center and do.scale are set to FALSE inside
        new.residuals <- ScaleData(
          new.residuals,
          features = NULL,
          #vars.to.regress = vars.to.regress,
          #latent.data = cell.attr[, vars.to.regress, drop = FALSE],
          model.use = 'linear',
          use.umi = FALSE,
          do.scale = do.scale,
          do.center = do.center,
          scale.max = Inf,
          block.size = 750,
          min.cells.to.block = 3000,
          verbose = verbose
        )
        assay.out <- CreateSCTAssayObject(
          counts = corrected_counts,
          scale.data = new.residuals,
          SCTModel.list = SCTModel.list
          )
        assay.out$data <- log1p(x = corrected_counts)
        VariableFeatures(assay.out) <- variable.features
        # one assay per dataset
        if (verbose){
          message("Finished calculating residuals for ", dataset.names[dataset.index])
        }
        sct.assay.list[[dataset.names[dataset.index]]] <- assay.out
        variable.feature.list[[dataset.names[dataset.index]]] <- VariableFeatures(assay.out)
      }
    } else { ### With reference model
      sct.assay.list.temp <- list()
      for (i in seq_len(length.out = length(x = cells.grid))) {
        vp <- cells.grid[[i]]

        if (verbose){
          message("Getting residuals for block ", i, "(of ", length(cells.grid), ") for ", dataset.names[[dataset.index]], " dataset")
        }

        assay.out <- GetSCT.Chunked(vp = vp,
                                    reference.SCT.model = reference.SCT.model,
                                    do.correct.umi = do.correct.umi)
        sct.assay.list.temp[[paste0("chunk", i)]] <- assay.out
        }
      if (length(sct.assay.list.temp)>1){
        # this currently fails in merge.StdAssay step
        # assignment of an object of class “list” is not valid for
        # slot ‘key’ in an object of class “Assay”; is(value, "character") is not TRUE
        assay.out <- merge(x = sct.assay.list.temp[[1]],
                           y = sct.assay.list.temp[2:length(sct.assay.list.temp)])

        } else {
          assay.out <- sct.assay.list.temp[[1]]
        }
      ## DoScaling
      scale.data <- GetAssayData(object = assay.out, slot = "scale.data")
      # scale data here as do.center and do.scale are set to FALSE inside
      scale.data <- ScaleData(
        scale.data,
        features = NULL,
        #vars.to.regress = vars.to.regress,
        #latent.data = cell.attr[, vars.to.regress, drop = FALSE],
        model.use = 'linear',
        use.umi = FALSE,
        do.scale = do.scale,
        do.center = do.center,
        scale.max = Inf,
        block.size = 750,
        min.cells.to.block = 3000,
        verbose = verbose
      )
      assay.out <- SetAssayData(object = assay.out, slot = "scale.data", new.data = scale.data)
      if (verbose){
        message("Finished calculating residuals for ", dataset.names[dataset.index])
      }
      sct.assay.list[[dataset.names[dataset.index]]] <- assay.out
      variable.feature.list[[dataset.names[dataset.index]]] <- rownames(assay.out)
    }
  }
# Return array by merging everythin
  if (length(x = sct.assay.list) > 1){
    vf.list <- lapply(X  = sct.assay.list, FUN = function(object.i) VariableFeatures(object = object.i))
    variable.features.union <- Reduce(f = union, x = vf.list)
    var.features.sorted <- sort(
      x = table(unlist(x = vf.list, use.names = FALSE)),
      decreasing = TRUE
    )
    # idx <- which(x = var.features == length(x = sct.assay.list))
    # select top ranking features
    #var.features <- names(x = var.features.sorted[1:variable.features.n])
    # calculate residuals for union of features
    var.features <- variable.features.union
    for (layer.name in names(sct.assay.list)){
      vst_out <- SCTModel_to_vst(SCTModel = slot(object = sct.assay.list[[layer.name]], name = "SCTModel.list")[[1]])
      all_cells <-  Cells(x = object, layer = paste0(layer, ".", layer.name))
      all_features <- Features(x = object, layer = paste0(layer, ".", layer.name))
      variable.features.target <- intersect(x = rownames(x = vst_out$model_pars_fit), y = var.features)
      variable.features.target <- setdiff(x = variable.features.target, y = VariableFeatures(sct.assay.list[[layer.name]]))
      if (length(variable.features.target )<1){
        next
      }
      counts <- LayerData(
        object = object,
        layer = paste0(layer, ".", layer.name),
        cells = all_cells
      )
      cells.grid <- DelayedArray::colAutoGrid(x = counts, ncol = ncol(counts))
      vp <- cells.grid[[1L]]
      block <- DelayedArray::read_block(x = counts, viewport = vp, as.sparse = TRUE)
      counts.vp <- as(object = block, Class = 'dgCMatrix')

      if (vst_out$arguments$min_var == "umi_median"){
        nz_median <- median(counts.vp@x)
        min_var_custom <- (nz_median / 5)^2
      } else {
        min_var_custom <- vst_out$arguments$min_var
      }
      vst_out$cell_attr <- vst_out$cell_attr[, c("log_umi"), drop=FALSE]
      vst_out$model_pars_fit <- vst_out$model_pars_fit[variable.features.target,,drop=FALSE]

      new_residual <- get_residuals(
        vst_out = vst_out,
        umi = counts.vp[variable.features.target,],
        residual_type = "pearson",
        min_variance = min_var_custom,
        verbosity =  FALSE
      )
      old_residual <- GetAssayData(object = sct.assay.list[[layer.name]], slot = 'scale.data')
      merged_residual <- rbind(old_residual, new_residual)
      sct.assay.list[[layer.name]] <- SetAssayData(object = sct.assay.list[[layer.name]], slot = 'scale.data', new.data = merged_residual)
      VariableFeatures(sct.assay.list[[layer.name]]) <- rownames(x = merged_residual)
    }
    merged.assay <- merge(x = sct.assay.list[[1]], y = sct.assay.list[2:length(sct.assay.list)])
    VariableFeatures(object = merged.assay) <- VariableFeatures(
      object = merged.assay,
      use.var.features = FALSE,
      nfeatures = variable.features.n
      )
    # set the names of SCTmodels to be layer names
    models <- slot(object = merged.assay, name="SCTModel.list")
    names(models) <- names(x = sct.assay.list)
    slot(object = merged.assay, name="SCTModel.list") <- models
  } else {
    merged.assay <- sct.assay.list[[1]]
  }
  gc(verbose = FALSE)
  return(merged.assay)
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
FetchResiduals <- function(object,
                           features,
                           assay = NULL,
                           umi.assay = "RNA",
                           layer = "counts",
                           clip.range = NULL,
                           reference.SCT.model = NULL,
                           replace.value = FALSE,
                           na.rm = TRUE,
                           verbose = TRUE) {
  assay <- assay %||% DefaultAssay(object = object)
  if (IsSCT(assay = object[[assay]])) {
    object[[assay]] <- as(object[[assay]], "SCTAssay")
  }
  if (!inherits(x = object[[assay]], what = "SCTAssay")) {
    stop(assay, " assay was not generated by SCTransform")
  }
  sct.models <- levels(x = object[[assay]])
  if (length(sct.models)==1){
    sct.models <- list(sct.models)
  }
  if (length(x = sct.models) == 0) {
    warning("SCT model not present in assay", call. = FALSE, immediate. = TRUE)
    return(object)
  }
  possible.features <- Reduce(f = union, x = lapply(X = sct.models, FUN = function(x) {
    rownames(x = SCTResults(object = object[[assay]], slot = "feature.attributes", model = x))
  }))
  bad.features <- setdiff(x = features, y = possible.features)
  if (length(x = bad.features) > 0) {
    warning("The following requested features are not present in any models: ",
            paste(bad.features, collapse = ", "),
            call. = FALSE
    )
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
  if (length(features) < 1){
    warning("The following requested features are not present in all the models: ",
            paste(features.orig, collapse = ", "),
            call. = FALSE
    )
    return (NULL)
  }  #if (length(x = sct.models) > 1 & verbose) {
  #  message("This SCTAssay contains multiple SCT models. Computing residuals for cells using")
  #}

  # Get all (count) layers
  layers <- Layers(object = object[[umi.assay]], search = layer)

  # iterate over layer running sct model for each of the object names
  new.residuals <- list()
  total_cells <- 0
  all_cells <- c()
  if (!is.null(x = reference.SCT.model)) {
    if (inherits(x = reference.SCT.model, what = "SCTModel")) {
      reference.SCT.model <- SCTModel_to_vst(SCTModel = reference.SCT.model)
    }
    if (is.list(x = reference.SCT.model) & inherits(x = reference.SCT.model[[1]], what = "SCTModel")) {
      stop("reference.SCT.model must be one SCTModel rather than a list of SCTModel")
    }
    if (reference.SCT.model$model_str != "y ~ log_umi") {
      stop("reference.SCT.model must be derived using default SCT regression formula, `y ~ log_umi`")
    }
  }
  for (i in seq_along(along.with = layers)) {
    l <- layers[i]
    sct_model <- sct.models[[i]]
    # these cells belong to this layer
    layer_cells <- Cells(x = object[[umi.assay]], layer = l)
    all_cells <- c(all_cells, layer_cells)
    total_cells <- total_cells + length(layer_cells)
    # calculate residual using this model and these cells
    new.residuals[[i]] <- FetchResidualSCTModel(
      object = object,
      umi.assay = umi.assay,
      assay = assay,
      layer = l,
      layer.cells = layer_cells,
      SCTModel = sct_model,
      reference.SCT.model = reference.SCT.model,
      new_features = features,
      replace.value = replace.value,
      clip.range = clip.range,
      verbose = verbose
    )
  }

  existing.data <- GetAssayData(object = object, slot = "scale.data", assay = assay)
  all.features <- union(x = rownames(x = existing.data), y = features)
  new.scale <- matrix(
    data = NA,
    nrow = length(x = all.features),
    ncol = total_cells,
    dimnames = list(all.features, all_cells)
  )
  common_cells <- intersect(colnames(new.scale), colnames(existing.data))
  if (nrow(x = existing.data) > 0) {
    new.scale[rownames(x = existing.data), common_cells] <- existing.data[, common_cells]
  }
  if (length(x = new.residuals) == 1 & is.list(x = new.residuals)) {
    new.residuals <- new.residuals[[1]]
  } else {
    new.residuals <- Reduce(cbind, new.residuals)
    #new.residuals <- matrix(data = unlist(new.residuals), nrow = nrow(new.scale) , ncol = ncol(new.scale))
    #colnames(new.residuals) <- colnames(new.scale)
    #rownames(new.residuals) <- rownames(new.scale)
  }
  new.scale[rownames(x = new.residuals), colnames(x = new.residuals)] <- new.residuals

  if (na.rm) {
    new.scale <- new.scale[!rowAnyNAs(x = new.scale), ]
  }

  return(new.scale[features,])
}


# Calculate pearson residuals of features not in the scale.data
# This function is the secondary function under FetchResiduals
#
# @param object A seurat object
# @param assay Name of the assay of the seurat object generated by SCTransform. Default
# is "SCT"
# @param umi.assay Name of the assay of the seurat object to fetch UMIs from. Default
# is "RNA"
# @param layer Name of the layer under `umi.assay` to fetch UMIs from. Default is
# "counts"
# @param layer.cells Vector of cells to calculate the residual for. Default is NULL
# which uses all cells in the layer
# @param SCTModel Which SCTmodel to use from the object for calculating the residual.
# Will be ignored if reference.SCT.model is set
# @param reference.SCT.model If a reference SCT model should be used for calculating
# the residuals. When set to not NULL, ignores the `SCTModel` paramater.
# @param new_features A vector of features to calculate the residuals for
# @param clip.range Numeric of length two specifying the min and max values the Pearson residual will be clipped to. Useful if you want to change the clip.range.
# @param replace.value Whether to replace the value of residuals if it already exists
# @param verbose Whether to print messages and progress bars
#
# @return Returns a matrix containing centered pearson residuals of added features
#
#' @importFrom sctransform get_residuals
#' @importFrom Matrix colSums
#
FetchResidualSCTModel <- function(object,
                                  assay = "SCT",
                                  umi.assay = "RNA",
                                  layer = "counts",
                                  layer.cells = NULL,
                                  SCTModel = NULL,
                                  reference.SCT.model = NULL,
                                  new_features = NULL,
                                  clip.range = NULL,
                                  replace.value = FALSE,
                                  verbose = FALSE) {
  model.cells <- character()
  model.features <- Features(x = object, assay = assay)
  if (is.null(x = reference.SCT.model)){
    clip.range <- clip.range %||% SCTResults(object = object[[assay]], slot = "clips", model = SCTModel)$sct
    model.features <- rownames(x = SCTResults(object = object[[assay]], slot = "feature.attributes", model = SCTModel))
    model.cells <- Cells(x = slot(object = object[[assay]], name = "SCTModel.list")[[SCTModel]])
    sct.method <- SCTResults(object = object[[assay]], slot = "arguments", model = SCTModel)$sct.method %||% "default"
  }

  layer.cells <- layer.cells %||% Cells(x = object[[umi.assay]], layer = layer)
  if (!is.null(reference.SCT.model)) {
    # use reference SCT model
    sct.method <- "reference"
  }
  existing.scale.data <- NULL
  if (is.null(x=reference.SCT.model)){
    existing.scale.data <- suppressWarnings(GetAssayData(object = object, assay = assay, slot = "scale.data"))
  }
  scale.data.cells <- colnames(x = existing.scale.data)
  scale.data.cells.common <- intersect(scale.data.cells, layer.cells)
  scale.data.cells <- intersect(x = scale.data.cells, y = scale.data.cells.common)
  if (length(x = setdiff(x = layer.cells, y = scale.data.cells)) == 0) {
    # existing.scale.data <- suppressWarnings(GetAssayData(object = object, assay = assay, slot = "scale.data"))
    #full.scale.data <- matrix(data = NA, nrow = nrow(x = existing.scale.data),
    #                          ncol = length(x = layer.cells), dimnames = list(rownames(x = existing.scale.data), layer.cells))
    #full.scale.data[rownames(x = existing.scale.data), colnames(x = existing.scale.data)] <- existing.scale.data
    #existing_features <- names(x = which(x = !apply(
    #  X = full.scale.data,
    #  MARGIN = 1,
    #  FUN = anyNA
    #)))
    existing_features <- rownames(x = existing.scale.data)
  } else {
    existing_features <- character()
  }
  if (replace.value) {
    features_to_compute <- new_features
  } else {
    features_to_compute <- setdiff(x = new_features, y = existing_features)
  }
  if (length(features_to_compute)<1){
    return (existing.scale.data[intersect(x = rownames(x = scale.data.cells), y = new_features),,drop=FALSE])
  }

  if (is.null(x = reference.SCT.model) & length(x = setdiff(x = model.cells, y =  scale.data.cells)) == 0) {
    existing_features <- names(x = which(x = ! apply(
      X = GetAssayData(object = object, assay = assay, slot = "scale.data")[, model.cells],
      MARGIN = 1,
      FUN = anyNA)
    ))
  } else {
    existing_features <- character()
  }
  if (sct.method == "reference.model") {
    if (verbose) {
      message("sct.model ", SCTModel, " is from reference, so no residuals will be recalculated")
    }
    features_to_compute <- character()
  }
  if (!umi.assay %in% Assays(object = object)) {
    warning("The umi assay (", umi.assay, ") is not present in the object. ",
            "Cannot compute additional residuals.",
            call. = FALSE, immediate. = TRUE
    )
    return(NULL)
  }
  # these features do not have feature attriutes
  diff_features <- setdiff(x = features_to_compute, y = model.features)
  intersect_features <- intersect(x = features_to_compute, y = model.features)
  if (sct.method == "reference") {
    vst_out <- SCTModel_to_vst(SCTModel = reference.SCT.model)

    # override clip.range
    clip.range <- vst_out$arguments$sct.clip.range
    umi.field <- paste0("nCount_", assay)
    # get rid of the cell attributes
    vst_out$cell_attr <- NULL
    all.features <- intersect(
      x = rownames(x = vst_out$gene_attr),
      y = features_to_compute
    )
    vst_out$gene_attr <- vst_out$gene_attr[all.features, , drop = FALSE]
    vst_out$model_pars_fit <- vst_out$model_pars_fit[all.features, , drop = FALSE]
  } else {
    vst_out <- SCTModel_to_vst(SCTModel = slot(object = object[[assay]], name = "SCTModel.list")[[SCTModel]])
    clip.range <- vst_out$arguments$sct.clip.range
  }
  clip.max <- max(clip.range)
  clip.min <- min(clip.range)

  layer.cells <- layer.cells %||% Cells(x = object[[umi.assay]], layer = layer)
  if (length(x = diff_features) == 0) {
    counts <- LayerData(
      object = object[[umi.assay]],
      layer = layer,
      cells = layer.cells
    )

    # iterate over 2k cells at once
    #cells.grid <- DelayedArray::colAutoGrid(x = counts, ncol = min(2000, length(x = layer.cells)))
    cells.grid <- DelayedArray::colAutoGrid(x = counts, ncol = length(x = layer.cells))
    new_residuals <- list()

    for (i in seq_len(length.out = length(x = cells.grid))) {
      vp <- cells.grid[[i]]
      block <- DelayedArray::read_block(x = counts, viewport = vp, as.sparse = TRUE)
      ## TODO: Maybe read only interesting genes
      umi.all <- as(object = block, Class = "dgCMatrix")

      # calculate min_variance for get_residuals
      # required when vst_out$arguments$min_variance == "umi_median"
      # only calculated once
      if (i==1){
        nz_median <- median(umi.all@x)
        min_var_custom <- (nz_median / 5)^2
      }
      umi <- umi.all[features_to_compute, , drop = FALSE]

      ## Add cell_attr for missing cells
      cell_attr <- data.frame(
        umi = colSums(umi.all),
        log_umi = log10(x = colSums(umi.all))
      )
      rownames(cell_attr) <- colnames(umi.all)
      if (sct.method %in% c("reference.model", "reference")) {
        vst_out$cell_attr <- cell_attr[colnames(umi.all), ,drop=FALSE]
      } else {
        cell_attr_existing <- vst_out$cell_attr
        cells_missing <- setdiff(rownames(cell_attr), rownames(cell_attr_existing))
        if (length(cells_missing)>0){
          cell_attr_missing <- cell_attr[cells_missing, ,drop=FALSE]
          missing_cols <- setdiff(x = colnames(x = cell_attr_existing),
                                  y = colnames(x = cell_attr_missing))

          if (length(x = missing_cols) > 0) {
            cell_attr_missing[, missing_cols] <- NA
          }
          vst_out$cell_attr <- rbind(cell_attr_existing,
                                     cell_attr_missing)
          vst_out$cell_attr <- vst_out$cell_attr[colnames(umi), , drop=FALSE]
        }
      }
      if (verbose) {
        if (sct.method == "reference.model") {
          message("using reference sct model")
        } else {
          message("sct.model: ", SCTModel, " on ", ncol(x = umi), " cells: ",
                  colnames(x = umi.all)[1], " .. ", colnames(x = umi.all)[ncol(umi.all)])
        }
      }

      if (vst_out$arguments$min_variance == "umi_median"){
        min_var <- min_var_custom
      } else {
        min_var <- vst_out$arguments$min_variance
      }
      if (nrow(umi)>0){
        new_residual <- get_residuals(
          vst_out = vst_out,
          umi = umi,
          residual_type = "pearson",
          min_variance = min_var,
          res_clip_range = c(clip.min, clip.max),
          verbosity = as.numeric(x = verbose) * 2
        )
      } else {
        return(matrix(
          data = NA,
          nrow = length(x = features_to_compute),
          ncol = length(x = colnames(umi.all)),
          dimnames = list(features_to_compute, colnames(umi.all))
        ))
      }
      new_residual <- as.matrix(x = new_residual)
      new_residuals[[i]] <- new_residual
    }
    new_residual <- do.call(what = cbind, args = new_residuals)
    # centered data if no reference model is provided
    if (is.null(x = reference.SCT.model)){
      new_residual <- new_residual - rowMeans(x = new_residual)
    } else {
      # subtract residual mean from reference model
      if (verbose){
        message("Using residual mean from reference for centering")
      }
      vst_out <- SCTModel_to_vst(SCTModel = reference.SCT.model)
      ref.residuals.mean <- vst_out$gene_attr[rownames(x = new_residual),"residual_mean"]
      new_residual <- sweep(
        x = new_residual,
        MARGIN = 1,
        STATS = ref.residuals.mean,
        FUN = "-"
      )
    }
    # return (new_residuals)
  } else {
    #  Some features do not exist
    warning(
      "In the SCTModel ", SCTModel, ", the following ", length(x = diff_features),
      " features do not exist in the counts slot: ", paste(diff_features, collapse = ", ")
    )
    if (length(x = intersect_features) == 0) {
      # No features exist
      return(matrix(
        data = NA,
        nrow = length(x = features_to_compute),
        ncol = length(x = model.cells),
        dimnames = list(features_to_compute, model.cells)
      ))
    }
  }
  old.features <- setdiff(x = new_features, y = features_to_compute)
  if (length(x = old.features) > 0) {
    old_residuals <- GetAssayData(object = object[[assay]], slot = "scale.data")[old.features, model.cells, drop = FALSE]
    new_residual <- rbind(new_residual, old_residuals)[new_features, ]
  }
  return(new_residual)
}

#' temporal function to get residuals from reference
#' @importFrom sctransform get_residuals
#' @importFrom Matrix colSums
#'

FetchResiduals_reference <- function(object,
                                     reference.SCT.model = NULL,
                                     features = NULL,
                                     nCount_UMI = NULL,
                                     verbose = FALSE) {
  ## Add cell_attr for missing cells
  nCount_UMI <- nCount_UMI %||% colSums(object)
  cell_attr <- data.frame(
    umi = nCount_UMI,
    log_umi = log10(x = nCount_UMI)
  )
  features_to_compute <- features
  features_to_compute <- intersect(features_to_compute, rownames(object))
  vst_out <- SCTModel_to_vst(SCTModel = reference.SCT.model)

  # override clip.range
  clip.range <- vst_out$arguments$sct.clip.range
  # get rid of the cell attributes
  vst_out$cell_attr <- NULL
  all.features <- intersect(
    x = rownames(x = vst_out$gene_attr),
    y = features_to_compute
  )
  vst_out$gene_attr <- vst_out$gene_attr[all.features, , drop = FALSE]
  vst_out$model_pars_fit <- vst_out$model_pars_fit[all.features, , drop = FALSE]

  clip.max <- max(clip.range)
  clip.min <- min(clip.range)


  umi <- object[features_to_compute, , drop = FALSE]


  rownames(cell_attr) <- colnames(object)
  vst_out$cell_attr <- cell_attr

  if (verbose) {
    message("using reference sct model")
  }

  if (vst_out$arguments$min_variance == "umi_median"){
    nz_median <- 1
    min_var_custom <- (nz_median / 5)^2
    min_var <- min_var_custom
  } else {
    min_var <- vst_out$arguments$min_variance
  }
  new_residual <- get_residuals(
    vst_out = vst_out,
    umi = umi,
    residual_type = "pearson",
    min_variance = min_var,
    verbosity = as.numeric(x = verbose) * 2
  )

  ref.residuals.mean <- vst_out$gene_attr[rownames(x = new_residual),"residual_mean"]
  new_residual <- sweep(
    x = new_residual,
    MARGIN = 1,
    STATS = ref.residuals.mean,
    FUN = "-"
  )
  new_residual <- MinMax(data = new_residual, min = clip.min, max = clip.max)
  return(new_residual)
}

#' Find variable features based on mean.var.plot
#'
MVP <- function(
  data,
  verbose = TRUE,
  nselect = 2000L,
  mean.cutoff = c(0.1, 8),
  dispersion.cutoff = c(1, Inf),
  ...
) {
  hvf.info <- DISP(data = data, nselect = nselect, verbose = verbose)
  hvf.info$variable <- FALSE
  means.use <- (hvf.info[, 1] > mean.cutoff[1]) & (hvf.info[, 1] < mean.cutoff[2])
  dispersions.use <- (hvf.info[, 3] > dispersion.cutoff[1]) & (hvf.info[, 3] < dispersion.cutoff[2])
  hvf.info[which(x = means.use & dispersions.use), 'variable'] <- TRUE
  hvf.info[hvf.info$variable,'rank'] <- rank(x = hvf.info[hvf.info$variable,'rank'])
  hvf.info[!hvf.info$variable,'rank'] <- NA
  return(hvf.info)
}
