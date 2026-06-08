#%% Portable single-file Seurat bundles %%#
#
# SaveSeurat()/LoadSeurat() write a Seurat object -- including any out-of-memory
# layers (BPCells IterableMatrix, DelayedMatrix/HDF5) -- to a single, portable
# archive file. On-disk matrix stores are copied into the archive and their
# paths are reconnected on load, so the bundle can be moved or shared like an
# AnnData .h5ad without worrying about absolute paths. Fully in-memory objects
# are simply serialized inside the archive.
NULL

# List the out-of-memory layers of a Seurat object.
#
# @param object A Seurat object
# @return A data.frame with columns assay, layer, type ("bpcells"/"delayed")
#
.OnDiskLayers <- function(object) {
  out <- data.frame(
    assay = character(),
    layer = character(),
    type = character(),
    stringsAsFactors = FALSE
  )
  for (assay in Assays(object = object)) {
    for (lyr in Layers(object = object, assay = assay)) {
      ldat <- LayerData(object = object, layer = lyr, assay = assay)
      type <- if (inherits(x = ldat, what = 'IterableMatrix')) {
        'bpcells'
      } else if (inherits(x = ldat, what = 'DelayedMatrix')) {
        'delayed'
      } else {
        next
      }
      out <- rbind(out, data.frame(
        assay = assay, layer = lyr, type = type, stringsAsFactors = FALSE
      ))
    }
  }
  return(out)
}

# Write a single on-disk layer to a canonical store directory.
#
# @param ldat An IterableMatrix or DelayedMatrix
# @param dir Destination store directory (created if needed)
# @param type "bpcells" or "delayed"
# @return Invisibly, the reopened on-disk matrix pointing at \code{dir}
#
.WriteLayerStore <- function(ldat, dir, type) {
  dir.create(path = dir, recursive = TRUE, showWarnings = FALSE)
  if (type == 'bpcells') {
    if (!requireNamespace('BPCells', quietly = TRUE)) {
      stop("Package 'BPCells' is required to bundle BPCells layers.", call. = FALSE)
    }
    BPCells::write_matrix_dir(mat = ldat, dir = dir, overwrite = TRUE)
    return(invisible(x = BPCells::open_matrix_dir(dir = dir)))
  }
  if (!requireNamespace('HDF5Array', quietly = TRUE)) {
    stop("Package 'HDF5Array' is required to bundle DelayedMatrix layers.", call. = FALSE)
  }
  fp <- file.path(dir, 'matrix.h5')
  if (file.exists(fp)) {
    file.remove(fp)
  }
  HDF5Array::writeTENxMatrix(x = ldat, filepath = fp, group = 'matrix', verbose = FALSE)
  return(invisible(x = HDF5Array::TENxMatrix(filepath = fp, group = 'matrix')))
}

# Reopen an on-disk layer store from a bundle directory.
#
# @param dir Store directory
# @param type "bpcells" or "delayed"
# @return The reopened on-disk matrix
#
.OpenLayerStore <- function(dir, type) {
  if (type == 'bpcells') {
    return(BPCells::open_matrix_dir(dir = dir))
  }
  return(HDF5Array::TENxMatrix(filepath = file.path(dir, 'matrix.h5'), group = 'matrix'))
}

#' Save a Seurat object to a single portable file
#'
#' Write a Seurat object to one self-contained archive that bundles the object
#' together with any out-of-memory layers (BPCells \code{IterableMatrix},
#' \code{DelayedMatrix}/HDF5). The on-disk matrix stores are copied into the
#' archive and reconnected by \code{\link{LoadSeurat}}, so the file can be moved,
#' copied or shared without worrying about the absolute paths of the original
#' on-disk directories -- the same portability as an AnnData \code{.h5ad}. A
#' fully in-memory object is simply serialized inside the archive.
#'
#' @param object A \code{\link[SeuratObject]{Seurat}} object
#' @param file Path to the output archive (a gzip-compressed tarball; the
#'   \code{.seurat} extension is conventional)
#' @param verbose Print progress messages (default \code{TRUE})
#'
#' @return Invisibly, \code{file}
#'
#' @seealso \code{\link{LoadSeurat}}, \code{\link{SaveSeuratH5}},
#'   \code{\link{AsInMemory}}
#'
#' @export
#' @concept utilities
#'
#' @examples
#' \dontrun{
#' SaveSeurat(object = obj, file = "atlas.seurat")
#' obj <- LoadSeurat("atlas.seurat")
#' }
#'
SaveSeurat <- function(object, file, verbose = TRUE) {
  if (!requireNamespace('jsonlite', quietly = TRUE)) {
    stop("Package 'jsonlite' is required by SaveSeurat().", call. = FALSE)
  }
  stage <- tempfile(pattern = 'seurat-bundle')
  dir.create(path = stage, recursive = TRUE, showWarnings = FALSE)
  on.exit(expr = unlink(x = stage, recursive = TRUE), add = TRUE)
  ondisk <- .OnDiskLayers(object = object)
  manifest <- list(format_version = '1', layers = list())
  if (nrow(x = ondisk) > 0) {
    for (i in seq_len(length.out = nrow(x = ondisk))) {
      assay <- ondisk$assay[i]
      lyr <- ondisk$layer[i]
      type <- ondisk$type[i]
      store.rel <- file.path('stores', paste0(assay, '__', lyr))
      if (isTRUE(x = verbose)) {
        message('Bundling on-disk layer "', lyr, '" of assay "', assay, '"')
      }
      reopened <- .WriteLayerStore(
        ldat = LayerData(object = object, layer = lyr, assay = assay),
        dir = file.path(stage, store.rel),
        type = type
      )
      # Point the in-object layer at the staged store so object.rds stays small;
      # LoadSeurat reconnects from the manifest regardless of this path.
      LayerData(object = object, layer = lyr, assay = assay) <- reopened
      manifest$layers[[length(x = manifest$layers) + 1L]] <- list(
        assay = assay, layer = lyr, type = type, store = store.rel
      )
    }
  }
  saveRDS(object = object, file = file.path(stage, 'object.rds'))
  writeLines(
    text = jsonlite::toJSON(x = manifest, auto_unbox = TRUE, pretty = TRUE),
    con = file.path(stage, 'manifest.json')
  )
  file.abs <- normalizePath(path = file, mustWork = FALSE)
  old.wd <- setwd(dir = stage)
  on.exit(expr = setwd(dir = old.wd), add = TRUE)
  utils::tar(
    tarfile = file.abs,
    files = c('object.rds', 'manifest.json', if (nrow(x = ondisk) > 0) 'stores'),
    compression = 'gzip'
  )
  if (isTRUE(x = verbose)) {
    message('Saved Seurat bundle to ', file)
  }
  return(invisible(x = file))
}

#' Load a Seurat object from a portable bundle
#'
#' Read a Seurat object written by \code{\link{SaveSeurat}}, extracting the
#' archive and reconnecting any bundled out-of-memory layers to their extracted
#' on-disk stores. The object is usable immediately, regardless of where the
#' bundle was created.
#'
#' @param file Path to a bundle written by \code{\link{SaveSeurat}}
#' @param dir Directory into which the bundle is extracted; the on-disk layer
#'   stores live here for the lifetime of the returned object (default: a fresh
#'   temporary directory)
#' @param verbose Print progress messages (default \code{TRUE})
#'
#' @return A \code{\link[SeuratObject]{Seurat}} object
#'
#' @seealso \code{\link{SaveSeurat}}
#'
#' @export
#' @concept utilities
#'
#' @examples
#' \dontrun{
#' obj <- LoadSeurat("atlas.seurat")
#' }
#'
LoadSeurat <- function(file, dir = tempfile(pattern = 'seurat-bundle'), verbose = TRUE) {
  if (!requireNamespace('jsonlite', quietly = TRUE)) {
    stop("Package 'jsonlite' is required by LoadSeurat().", call. = FALSE)
  }
  dir.create(path = dir, recursive = TRUE, showWarnings = FALSE)
  utils::untar(tarfile = normalizePath(path = file), exdir = dir)
  object <- readRDS(file = file.path(dir, 'object.rds'))
  manifest <- jsonlite::fromJSON(
    txt = file.path(dir, 'manifest.json'),
    simplifyVector = FALSE
  )
  for (entry in manifest$layers) {
    if (isTRUE(x = verbose)) {
      message('Reconnecting on-disk layer "', entry$layer, '" of assay "', entry$assay, '"')
    }
    LayerData(object = object, layer = entry$layer, assay = entry$assay) <-
      .OpenLayerStore(dir = file.path(dir, entry$store), type = entry$type)
  }
  return(object)
}
