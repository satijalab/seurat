#%% Single-file HDF5 Seurat container %%#
#
# SaveSeuratH5()/LoadSeuratH5() store a Seurat object in a single HDF5 file, the
# closest analogue to an AnnData .h5ad: assay layers are written as on-disk
# sparse (10x-style CSC) HDF5 datasets and read back lazily as DelayedMatrix, so
# the file is portable, cross-tool readable, and supports matrices larger than
# the 2^31 dgCMatrix limit. Cell metadata and dimensional reductions are stored
# alongside. Version 1 persists assays/layers, meta.data and reductions; graphs,
# images, commands and misc are not yet written (use SaveSeurat() to bundle a
# full object). Split layers are joined per assay on save.
NULL

# Write a character/numeric vector to an HDF5 file.
.H5Put <- function(file, name, value) {
  rhdf5::h5write(obj = value, file = file, name = name)
  invisible(x = NULL)
}

# Ensure a (possibly nested) group path exists in an HDF5 file.
.H5Group <- function(file, path) {
  parts <- strsplit(x = path, split = '/', fixed = TRUE)[[1L]]
  parts <- parts[nzchar(x = parts)]
  cur <- ''
  for (p in parts) {
    cur <- if (nzchar(x = cur)) paste(cur, p, sep = '/') else p
    suppressWarnings(expr = rhdf5::h5createGroup(file = file, group = cur))
  }
  invisible(x = NULL)
}

#' Save a Seurat object to a single HDF5 file
#'
#' Write a Seurat object to one HDF5 (\code{.h5}) container -- the closest
#' Seurat analogue to an AnnData \code{.h5ad}. Each assay layer is stored as an
#' on-disk 10x-style sparse HDF5 dataset and read back lazily as a
#' \code{\link[DelayedArray]{DelayedMatrix}} by \code{\link{LoadSeuratH5}}, so
#' the file is portable, cross-tool readable, and not bound by the
#' \code{dgCMatrix} \eqn{2^{31}} limit. Cell metadata (\code{meta.data}) and
#' dimensional reductions are stored alongside.
#'
#' Version 1 persists assays, layers, \code{meta.data} and reductions; graphs,
#' spatial images, commands and \code{misc} are not written -- use
#' \code{\link{SaveSeurat}} to bundle a complete object. Split layers are joined
#' per assay before writing.
#'
#' @param object A \code{\link[SeuratObject]{Seurat}} object
#' @param file Path to the output \code{.h5} file
#' @param verbose Print progress messages (default \code{TRUE})
#'
#' @return Invisibly, \code{file}
#'
#' @seealso \code{\link{LoadSeuratH5}}, \code{\link{SaveSeurat}}
#'
#' @export
#' @concept utilities
#'
#' @examples
#' \dontrun{
#' SaveSeuratH5(object = obj, file = "atlas.h5")
#' obj <- LoadSeuratH5("atlas.h5")
#' }
#'
SaveSeuratH5 <- function(object, file, verbose = TRUE) {
  for (pkg in c('rhdf5', 'HDF5Array', 'jsonlite')) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop("Package '", pkg, "' is required by SaveSeuratH5().", call. = FALSE)
    }
  }
  if (file.exists(file)) {
    file.remove(file)
  }
  rhdf5::h5createFile(file = file)
  on.exit(expr = rhdf5::h5closeAll(), add = TRUE)
  assays <- Assays(object = object)
  .H5Put(file = file, name = 'format_version', value = '1')
  .H5Put(file = file, name = 'default_assay', value = DefaultAssay(object = object))
  .H5Put(file = file, name = 'assay_names', value = assays)
  .H5Put(file = file, name = 'cell_names', value = colnames(x = object))
  # Cell metadata: one dataset per column plus a JSON column->class map.
  .H5Group(file = file, path = 'meta')
  meta <- object[[]]
  classes <- vapply(X = meta, FUN = function(x) class(x = x)[1L], FUN.VALUE = '')
  for (col in colnames(x = meta)) {
    value <- meta[[col]]
    if (is.factor(x = value)) {
      .H5Put(file = file, name = paste0('meta/', col), value = as.character(x = value))
    } else if (is.logical(x = value)) {
      .H5Put(file = file, name = paste0('meta/', col), value = as.integer(x = value))
    } else {
      .H5Put(file = file, name = paste0('meta/', col), value = value)
    }
  }
  .H5Put(file = file, name = 'meta_columns', value = colnames(x = meta))
  .H5Put(file = file, name = 'meta_classes', value = unname(obj = classes))
  # Assays and layers (split layers joined first).
  .H5Group(file = file, path = 'assays')
  for (assay in assays) {
    object.assay <- object[[assay]]
    layer.names <- Layers(object = object.assay)
    if (anyDuplicated(x = gsub(pattern = '\\.[^.]+$', replacement = '', x = layer.names)) > 0) {
      object.assay <- tryCatch(
        expr = JoinLayers(object = object.assay),
        error = function(e) object.assay
      )
      layer.names <- Layers(object = object.assay)
    }
    .H5Group(file = file, path = file.path('assays', assay))
    .H5Group(file = file, path = file.path('assays', assay, 'layers'))
    .H5Put(file = file, name = file.path('assays', assay, 'features'),
           value = rownames(x = object.assay))
    .H5Put(file = file, name = file.path('assays', assay, 'layer_names'),
           value = layer.names)
    for (lyr in layer.names) {
      if (isTRUE(x = verbose)) {
        message('Writing layer "', lyr, '" of assay "', assay, '"')
      }
      ldat <- as.DelayedMatrix(x = LayerData(object = object.assay, layer = lyr))
      HDF5Array::writeTENxMatrix(
        x = ldat,
        filepath = file,
        group = file.path('assays', assay, 'layers', lyr),
        verbose = FALSE
      )
    }
  }
  # Dimensional reductions (embeddings, and loadings when present).
  reductions <- Reductions(object = object)
  .H5Put(file = file, name = 'reduction_names', value = reductions %||% character(length = 0L))
  if (length(x = reductions) > 0) {
    .H5Group(file = file, path = 'reductions')
    for (rd in reductions) {
      dr <- object[[rd]]
      .H5Group(file = file, path = file.path('reductions', rd))
      emb <- Embeddings(object = dr)
      .H5Put(file = file, name = file.path('reductions', rd, 'embeddings'), value = emb)
      .H5Put(file = file, name = file.path('reductions', rd, 'embeddings_cols'),
             value = colnames(x = emb))
      load.mat <- Loadings(object = dr)
      if (length(x = load.mat) > 0 && nrow(x = load.mat) > 0) {
        .H5Put(file = file, name = file.path('reductions', rd, 'loadings'), value = load.mat)
        .H5Put(file = file, name = file.path('reductions', rd, 'loadings_rows'),
               value = rownames(x = load.mat))
      }
      .H5Put(file = file, name = file.path('reductions', rd, 'key'), value = Key(object = dr))
      .H5Put(file = file, name = file.path('reductions', rd, 'assay_used'),
             value = DefaultAssay(object = dr) %||% '')
    }
  }
  rhdf5::h5closeAll()
  if (isTRUE(x = verbose)) {
    message('Saved Seurat HDF5 container to ', file)
  }
  return(invisible(x = file))
}

#' Load a Seurat object from an HDF5 container
#'
#' Read a Seurat object written by \code{\link{SaveSeuratH5}}. Assay layers are
#' opened lazily as \code{\link[DelayedArray]{DelayedMatrix}} (on-disk, backed by
#' the HDF5 file), and \code{meta.data} and reductions are restored. The HDF5
#' file backs the returned object's layers, so keep it available for the
#' object's lifetime (or call \code{\link{AsInMemory}} to detach it).
#'
#' @param file Path to a \code{.h5} file written by \code{\link{SaveSeuratH5}}
#' @param verbose Print progress messages (default \code{TRUE})
#'
#' @return A \code{\link[SeuratObject]{Seurat}} object with on-disk
#'   \code{DelayedMatrix} layers
#'
#' @seealso \code{\link{SaveSeuratH5}}
#'
#' @export
#' @concept utilities
#'
#' @examples
#' \dontrun{
#' obj <- LoadSeuratH5("atlas.h5")
#' }
#'
LoadSeuratH5 <- function(file, verbose = TRUE) {
  for (pkg in c('rhdf5', 'HDF5Array')) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop("Package '", pkg, "' is required by LoadSeuratH5().", call. = FALSE)
    }
  }
  on.exit(expr = rhdf5::h5closeAll(), add = TRUE)
  assay.names <- as.character(x = rhdf5::h5read(file = file, name = 'assay_names'))
  default.assay <- as.character(x = rhdf5::h5read(file = file, name = 'default_assay'))[1L]
  cell.names <- as.character(x = rhdf5::h5read(file = file, name = 'cell_names'))
  open.layers <- function(assay) {
    layer.names <- as.character(x = rhdf5::h5read(
      file = file, name = file.path('assays', assay, 'layer_names')
    ))
    mats <- lapply(X = layer.names, FUN = function(lyr) {
      HDF5Array::TENxMatrix(filepath = file, group = file.path('assays', assay, 'layers', lyr))
    })
    names(x = mats) <- layer.names
    return(mats)
  }
  # Build each assay from its layers, then assemble the Seurat object.
  build.assay <- function(assay) {
    mats <- open.layers(assay = assay)
    seed <- if ('counts' %in% names(x = mats)) 'counts' else names(x = mats)[1L]
    assay.obj <- CreateAssay5Object(counts = mats[[seed]])
    for (lyr in names(x = mats)) {
      LayerData(object = assay.obj, layer = lyr) <- mats[[lyr]]
    }
    if (!('counts' %in% names(x = mats))) {
      suppressWarnings(expr = LayerData(object = assay.obj, layer = 'counts') <- NULL)
    }
    return(assay.obj)
  }
  if (isTRUE(x = verbose)) {
    message('Reading assay "', default.assay, '"')
  }
  object <- CreateSeuratObject(counts = build.assay(assay = default.assay), assay = default.assay)
  for (assay in setdiff(x = assay.names, y = default.assay)) {
    if (isTRUE(x = verbose)) {
      message('Reading assay "', assay, '"')
    }
    object[[assay]] <- build.assay(assay = assay)
  }
  DefaultAssay(object = object) <- default.assay
  # Restore cell metadata.
  meta.columns <- as.character(x = rhdf5::h5read(file = file, name = 'meta_columns'))
  meta.classes <- as.character(x = rhdf5::h5read(file = file, name = 'meta_classes'))
  for (i in seq_along(along.with = meta.columns)) {
    col <- meta.columns[i]
    if (col %in% c('orig.ident', 'nCount_RNA', 'nFeature_RNA')) {
      next
    }
    value <- as.vector(x = rhdf5::h5read(file = file, name = paste0('meta/', col)))
    value <- switch(
      EXPR = meta.classes[i],
      factor = factor(x = as.character(x = value)),
      logical = as.logical(x = value),
      value
    )
    object[[col]] <- value
  }
  # Restore reductions.
  reduction.names <- tryCatch(
    expr = as.character(x = rhdf5::h5read(file = file, name = 'reduction_names')),
    error = function(e) character(length = 0L)
  )
  for (rd in reduction.names) {
    if (!nzchar(x = rd)) {
      next
    }
    emb <- as.matrix(x = rhdf5::h5read(file = file, name = file.path('reductions', rd, 'embeddings')))
    rownames(x = emb) <- cell.names
    colnames(x = emb) <- as.character(x = rhdf5::h5read(
      file = file, name = file.path('reductions', rd, 'embeddings_cols')
    ))
    key <- as.character(x = rhdf5::h5read(file = file, name = file.path('reductions', rd, 'key')))[1L]
    assay.used <- as.character(x = rhdf5::h5read(
      file = file, name = file.path('reductions', rd, 'assay_used')
    ))[1L]
    loadings <- tryCatch(
      expr = {
        lm <- as.matrix(x = rhdf5::h5read(file = file, name = file.path('reductions', rd, 'loadings')))
        rownames(x = lm) <- as.character(x = rhdf5::h5read(
          file = file, name = file.path('reductions', rd, 'loadings_rows')
        ))
        colnames(x = lm) <- colnames(x = emb)
        lm
      },
      error = function(e) new(Class = 'matrix')
    )
    object[[rd]] <- CreateDimReducObject(
      embeddings = emb,
      loadings = loadings,
      key = key,
      assay = if (nzchar(x = assay.used)) assay.used else default.assay
    )
  }
  rhdf5::h5closeAll()
  return(object)
}
