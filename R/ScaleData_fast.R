#' Fast ScaleData implementation
#'
#' @param object An object
#' @param ... Arguments passed to other methods
#'
#' @export
#'
ScaleData_fast <- function(object, ...) {
  UseMethod(generic = 'ScaleData_fast', object = object)
}

#' @importFrom future nbrOfWorkers
#' @importFrom future.apply future_lapply
#'
#' @inheritParams ScaleData
#' @param nthreads Number of threads to use for the single-call fast path
#' (centering/scaling without regression, splitting, or multiple workers).
#' Values greater than 1 enable RcppParallel threading; the progress bar is
#' only shown when \code{nthreads = 1}.
#'
#' @rdname ScaleData_fast
#' @concept preprocessing
#' @export
#'
ScaleData_fast.default <- function(
  object,
  features = NULL,
  vars.to.regress = NULL,
  latent.data = NULL,
  split.by = NULL,
  model.use = 'linear',
  use.umi = FALSE,
  do.scale = TRUE,
  do.center = TRUE,
  scale.max = 10,
  block.size = 1000,
  min.cells.to.block = 3000,
  nthreads = 1,
  verbose = TRUE,
  ...
) {
  CheckDots(...)
  features <- features %||% rownames(x = object)
  features <- as.vector(x = intersect(x = features, y = rownames(x = object)))
  # CHANGE: Defer the row-subset of `object` (previously performed here) so the
  # fast path below can pass the full matrix plus feature indices into C++ and
  # select the requested rows there, avoiding an extra sparse-matrix copy. The
  # general path performs the subset further down, exactly as before.
  object.names <- list(features, colnames(x = object))
  min.cells.to.block <- min(min.cells.to.block, ncol(x = object))
  suppressWarnings(expr = Parenting(
    parent.find = "ScaleData_fast.Assay",
    features = features,
    min.cells.to.block = min.cells.to.block
  ))
  split.by <- split.by %||% TRUE
  split.cells <- split(x = colnames(x = object), f = split.by)
  CheckGC()
  # CHANGE: Fast path. When no regression is requested, the data are not split,
  # and only one worker is available, compute the row statistics and materialise
  # the final dense matrix in a single C++ call. This avoids the per-block
  # subsetting, sparse transpose and dense intermediate copies that dominate on
  # large data. Feature selection happens inside C++ (via `features`), so we
  # operate on the full (un-subset) matrix here. `nthreads > 1` enables
  # RcppParallel threading; the progress bar is shown only on the serial path.
  if (
    is.null(x = vars.to.regress) &&
    is.null(x = latent.data) &&
    length(x = split.cells) == 1 &&
    nbrOfWorkers() == 1
  ) {
    if (verbose && (do.scale || do.center)) {
      msg <- paste(
        na.omit(object = c(
          ifelse(test = do.center, yes = 'centering', no = NA_character_),
          ifelse(test = do.scale, yes = 'scaling', no = NA_character_)
        )),
        collapse = ' and '
      )
      msg <- paste0(
        toupper(x = substr(x = msg, start = 1, stop = 1)),
        substr(x = msg, start = 2, stop = nchar(x = msg)),
        ' data matrix'
      )
      message(msg)
    }
    # 0-based row indices of the requested features within the full matrix.
    feature.idx <- match(x = features, table = rownames(x = object)) - 1L
    if (inherits(x = object, what = 'dgTMatrix')) {
      object <- as(object = object, Class = 'dgCMatrix')
    }
    if (is(object = object, class2 = 'dgCMatrix')) {
      scaled.data <- FastSparseRowScale_optimized(
        x = object@x,
        i = object@i,
        p = object@p,
        rows = nrow(x = object),
        cols = ncol(x = object),
        features = feature.idx,
        scale = do.scale,
        center = do.center,
        scale_max = scale.max,
        nthreads = nthreads,
        display_progress = verbose
      )
    } else {
      scaled.data <- FastDenseRowScale_optimized(
        mat = as.matrix(x = object),
        features = feature.idx,
        scale = do.scale,
        center = do.center,
        scale_max = scale.max,
        nthreads = nthreads,
        display_progress = verbose
      )
    }
    dimnames(x = scaled.data) <- object.names
    return(scaled.data)
  }
  # General path: perform the row-subset the remaining (regression / split /
  # multi-worker) code paths expect.
  object <- object[features, , drop = FALSE]
  if (!is.null(x = vars.to.regress)) {
    if (is.null(x = latent.data)) {
      latent.data <- data.frame(row.names = colnames(x = object))
    } else {
      latent.data <- latent.data[colnames(x = object), , drop = FALSE]
      rownames(x = latent.data) <- colnames(x = object)
    }
    if (any(vars.to.regress %in% rownames(x = object))) {
      latent.data <- cbind(
        latent.data,
        t(x = object[vars.to.regress[vars.to.regress %in% rownames(x = object)], , drop=FALSE])
      )
    }
    # Currently, RegressOutMatrix will do nothing if latent.data = NULL
    notfound <- setdiff(x = vars.to.regress, y = colnames(x = latent.data))
    if (length(x = notfound) == length(x = vars.to.regress)) {
      stop(
        "None of the requested variables to regress are present in the object.",
        call. = FALSE
      )
    } else if (length(x = notfound) > 0) {
      warning(
        "Requested variables to regress not in object: ",
        paste(notfound, collapse = ", "),
        call. = FALSE,
        immediate. = TRUE
      )
      vars.to.regress <- colnames(x = latent.data)
    }
    if (verbose) {
      message("Regressing out ", paste(vars.to.regress, collapse = ', '))
    }
    chunk.points <- ChunkPoints(dsize = nrow(x = object), csize = block.size)
    if (nbrOfWorkers() > 1) { # TODO: lapply
      chunks <- expand.grid(
        names(x = split.cells),
        1:ncol(x = chunk.points),
        stringsAsFactors = FALSE
      )
      object <- future_lapply(
        X = 1:nrow(x = chunks),
        FUN = function(i) {
          row <- chunks[i, ]
          group <- row[[1]]
          index <- as.numeric(x = row[[2]])
          return(RegressOutMatrix(
            data.expr = object[chunk.points[1, index]:chunk.points[2, index], split.cells[[group]], drop = FALSE],
            latent.data = latent.data[split.cells[[group]], , drop = FALSE],
            features.regress = features,
            model.use = model.use,
            use.umi = use.umi,
            verbose = FALSE
          ))
        }
      )
      if (length(x = split.cells) > 1) {
        merge.indices <- lapply(
          X = 1:length(x = split.cells),
          FUN = seq.int,
          to = length(x = object),
          by = length(x = split.cells)
        )
        object <- lapply(
          X = merge.indices,
          FUN = function(x) {
            return(do.call(what = 'rbind', args = object[x]))
          }
        )
        object <- do.call(what = 'cbind', args = object)
      } else {
        object <- do.call(what = 'rbind', args = object)
      }
    } else {
      object <- lapply(
        X = names(x = split.cells),
        FUN = function(x) {
          if (verbose && length(x = split.cells) > 1) {
            message("Regressing out variables from split ", x)
          }
          return(RegressOutMatrix(
            data.expr = object[, split.cells[[x]], drop = FALSE],
            latent.data = latent.data[split.cells[[x]], , drop = FALSE],
            features.regress = features,
            model.use = model.use,
            use.umi = use.umi,
            verbose = verbose
          ))
        }
      )
      object <- do.call(what = 'cbind', args = object)
    }
    dimnames(x = object) <- object.names
    CheckGC()
  }
  if (verbose && (do.scale || do.center)) {
    msg <- paste(
      na.omit(object = c(
        ifelse(test = do.center, yes = 'centering', no = NA_character_),
        ifelse(test = do.scale, yes = 'scaling', no = NA_character_)
      )),
      collapse = ' and '
    )
    msg <- paste0(
      toupper(x = substr(x = msg, start = 1, stop = 1)),
      substr(x = msg, start = 2, stop = nchar(x = msg)),
      ' data matrix'
    )
    message(msg)
  }
  # NOTE: The single-call fast path is handled earlier (before the row-subset),
  # so the remaining code only runs for the regression / split / multi-worker
  # cases that still require R-level feature blocking.
  if (inherits(x = object, what = c('dgCMatrix', 'dgTMatrix'))) {
    scale.function <- FastSparseRowScale
  } else {
    object <- as.matrix(x = object)
    scale.function <- FastRowScale
  }
  if (nbrOfWorkers() > 1) {
    blocks <- ChunkPoints(dsize = length(x = features), csize = block.size)
    chunks <- expand.grid(
      names(x = split.cells),
      1:ncol(x = blocks),
      stringsAsFactors = FALSE
    )
    scaled.data <- future_lapply(
      X = 1:nrow(x = chunks),
      FUN = function(index) {
        row <- chunks[index, ]
        group <- row[[1]]
        block <- as.vector(x = blocks[, as.numeric(x = row[[2]])])
        arg.list <- list(
          mat = object[features[block[1]:block[2]], split.cells[[group]], drop = FALSE],
          scale = do.scale,
          center = do.center,
          scale_max = scale.max,
          display_progress = FALSE
        )
        arg.list <- arg.list[intersect(x = names(x = arg.list), y = names(x = formals(fun = scale.function)))]
        data.scale <- do.call(what = scale.function, args = arg.list)
        dimnames(x = data.scale) <- dimnames(x = object[features[block[1]:block[2]], split.cells[[group]]])
        suppressWarnings(expr = data.scale[is.na(x = data.scale)] <- 0)
        CheckGC()
        return(data.scale)
      }
    )
    if (length(x = split.cells) > 1) {
      merge.indices <- lapply(
        X = 1:length(x = split.cells),
        FUN = seq.int,
        to = length(x = scaled.data),
        by = length(x = split.cells)
      )
      scaled.data <- lapply(
        X = merge.indices,
        FUN = function(x) {
          return(suppressWarnings(expr = do.call(what = 'rbind', args = scaled.data[x])))
        }
      )
      scaled.data <- suppressWarnings(expr = do.call(what = 'cbind', args = scaled.data))
    } else {
      suppressWarnings(expr = scaled.data <- do.call(what = 'rbind', args = scaled.data))
    }
  } else {
    scaled.data <- matrix(
      data = NA_real_,
      nrow = nrow(x = object),
      ncol = ncol(x = object),
      dimnames = object.names
    )
    max.block <- ceiling(x = length(x = features) / block.size)
    for (x in names(x = split.cells)) {
      if (verbose) {
        if (length(x = split.cells) > 1 && (do.scale || do.center)) {
          message(gsub(pattern = 'matrix', replacement = 'from split ', x = msg), x)
        }
        pb <- txtProgressBar(min = 0, max = max.block, style = 3, file = stderr())
      }
      for (i in 1:max.block) {
        my.inds <- ((block.size * (i - 1)):(block.size * i - 1)) + 1
        my.inds <- my.inds[my.inds <= length(x = features)]
        arg.list <- list(
          mat = object[features[my.inds], split.cells[[x]], drop = FALSE],
          scale = do.scale,
          center = do.center,
          scale_max = scale.max,
          display_progress = FALSE
        )
        arg.list <- arg.list[intersect(x = names(x = arg.list), y = names(x = formals(fun = scale.function)))]
        data.scale <- do.call(what = scale.function, args = arg.list)
        dimnames(x = data.scale) <- dimnames(x = object[features[my.inds], split.cells[[x]]])
        scaled.data[features[my.inds], split.cells[[x]]] <- data.scale
        rm(data.scale)
        CheckGC()
        if (verbose) {
          setTxtProgressBar(pb = pb, value = i)
        }
      }
      if (verbose) {
        close(con = pb)
      }
    }
  }
  dimnames(x = scaled.data) <- object.names
  scaled.data[is.na(x = scaled.data)] <- 0
  CheckGC()
  return(scaled.data)
}

#' @rdname ScaleData_fast
#' @concept preprocessing
#' @export
#' @method ScaleData_fast Assay
#'
ScaleData_fast.Assay <- function(
  object,
  features = NULL,
  vars.to.regress = NULL,
  latent.data = NULL,
  split.by = NULL,
  model.use = 'linear',
  use.umi = FALSE,
  do.scale = TRUE,
  do.center = TRUE,
  scale.max = 10,
  block.size = 1000,
  min.cells.to.block = 3000,
  nthreads = 1,
  verbose = TRUE,
  ...
) {
  use.umi <- ifelse(test = model.use != 'linear', yes = TRUE, no = use.umi)
  slot.use <- ifelse(test = use.umi, yes = 'counts', no = 'data')
  features <- features %||% VariableFeatures(object)
  if (length(x = features) == 0) {
    features <- rownames(x = GetAssayData(object = object, layer = slot.use))
  }
  object <- SetAssayData(
    object = object,
    layer = 'scale.data',
    new.data = ScaleData_fast(
      object = GetAssayData(object = object, layer = slot.use),
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
      nthreads = nthreads,
      verbose = verbose,
      ...
    )
  )
  suppressWarnings(expr = Parenting(
    parent.find = "ScaleData_fast.Seurat",
    features = features,
    min.cells.to.block = min.cells.to.block,
    use.umi = use.umi
  ))
  return(object)
}

#' @importFrom SeuratObject StitchMatrix
#'
#' @rdname ScaleData_fast
#' @concept preprocessing
#' @export
#' @method ScaleData_fast StdAssay
#'
ScaleData_fast.StdAssay <- function(
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
  nthreads = 1,
  save = 'scale.data',
  verbose = TRUE,
  ...
) {
  use.umi <- ifelse(test = model.use != 'linear', yes = TRUE, no = use.umi)
  olayer <- layer <- unique(x = layer)
  layer <- Layers(object = object, search = layer)
  if (is.null(layer)) {
    abort(paste0("No layer matching pattern '", olayer, "' found. Please run NormalizeData and retry"))
  }
  if (isTRUE(x = use.umi)) {
    layer <- "counts"
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
      LayerData(object = object, layer = save[i], ...) <- ScaleData_fast(
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
        nthreads = nthreads,
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
    ldata <- ScaleData_fast(
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
      nthreads = nthreads,
      verbose = verbose,
      ...
    )
    LayerData(object = object, layer = save, features = rownames(ldata)) <- ldata
  }
  return(object)
}

#' @param assay Name of Assay to scale
#'
#' @rdname ScaleData_fast
#' @concept preprocessing
#' @export
#' @method ScaleData_fast Seurat
#'
ScaleData_fast.Seurat <- function(
  object,
  features = NULL,
  assay = NULL,
  vars.to.regress = NULL,
  split.by = NULL,
  model.use = 'linear',
  use.umi = FALSE,
  do.scale = TRUE,
  do.center = TRUE,
  scale.max = 10,
  block.size = 1000,
  min.cells.to.block = 3000,
  nthreads = 1,
  verbose = TRUE,
  ...
) {
  assay <- assay[1L] %||% DefaultAssay(object = object)
  assay <- match.arg(arg = assay, choices = Assays(object = object))
  if (any(vars.to.regress %in% colnames(x = object[[]]))) {
    latent.data <- object[[vars.to.regress[vars.to.regress %in% colnames(x = object[[]])]]]
  } else {
    latent.data <- NULL
  }
  if (is.character(x = split.by) && length(x = split.by) == 1) {
    split.by <- object[[split.by]]
  }
  assay.data <- ScaleData_fast(
    object = object[[assay]],
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
    nthreads = nthreads,
    verbose = verbose,
    ...
  )
  object[[assay]] <- assay.data
  object <- LogSeuratCommand(object = object)
  return(object)
}
