#' @include generics.R
#' @include visualization.R
#'
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @param fov Name to store FOV as
#' @param assay Name to store expression matrix as
#' @param ... Ignored
#'
#' @return \code{LoadAkoya}: A \code{\link[SeuratObject]{Seurat}} object
#'
#' @importFrom SeuratObject Cells CreateFOV CreateSeuratObject
#'
#' @export
#'
#' @rdname ReadAkoya
#'
LoadAkoya <- function(
  filename,
  type = c('inform', 'processor', 'qupath'),
  fov,
  assay = 'Akoya',
  ...
) {
  # read in matrix and centroids
  data <- ReadAkoya(filename = filename, type = type)
  # convert centroids into coords object
  coords <- suppressWarnings(expr = CreateFOV(
    coords = data$centroids,
    type = 'centroids',
    key = 'fov',
    assay = assay
  ))
  colnames(x = data$metadata) <- suppressWarnings(
    expr = make.names(names = colnames(x = data$metadata))
  )
  # build Seurat object from matrix
  obj <- CreateSeuratObject(
    counts = data$matrix,
    assay = assay,
    meta.data = data$metadata
  )
  # make sure coords only contain cells in seurat object
  coords <- subset(x = coords, cells = Cells(x = obj))
  suppressWarnings(expr = obj[[fov]] <- coords) # add image to seurat object
  # Add additional assays
  for (i in setdiff(x = names(x = data), y = c('matrix', 'centroids', 'metadata'))) {
    suppressWarnings(expr = obj[[i]] <- CreateAssayObject(counts = data[[i]]))
  }
  return(obj)
}

#' @inheritParams ReadAkoya
#' @param data.dir Path to a directory containing Vitessce cells
#' and clusters JSONs
#'
#' @return \code{LoadHuBMAPCODEX}: A \code{\link[SeuratObject]{Seurat}} object
#'
#' @importFrom SeuratObject Cells CreateFOV CreateSeuratObject
#'
#' @export
#'
#' @rdname ReadVitessce
#'
LoadHuBMAPCODEX <- function(data.dir, fov, assay = 'CODEX') {
  data <- ReadVitessce(
    counts = file.path(data.dir, "reg1_stitched_expressions.clusters.json"),
    coords = file.path(data.dir, "reg1_stitched_expressions.cells.json"),
    type = "segmentations"
  )
  # Create spatial and Seurat objects
  coords <- CreateFOV(
    coords = data$segmentations,
    molecules = data$molecules,
    assay = assay
  )
  obj <- CreateSeuratObject(counts = data$counts, assay = assay)
  # make sure spatial coords only contain cells in seurat object
  coords <- subset(x = coords, cells = Cells(x = obj))
  obj[[fov]] <- coords
  return(obj)
}

#' @inheritParams ReadAkoya
#' @param data.dir Path to folder containing Nanostring SMI outputs
#'
#' @return \code{LoadNanostring}: A \code{\link[SeuratObject]{Seurat}} object
#'
#' @importFrom SeuratObject Cells CreateCentroids CreateFOV
#' CreateSegmentation CreateSeuratObject
#'
#' @export
#'
#' @rdname ReadNanostring
#'
LoadNanostring <- function(data.dir, fov, assay = 'Nanostring') {
  data <- ReadNanostring(
    data.dir = data.dir,
    type = c("centroids", "segmentations")
  )
  segs <- CreateSegmentation(data$segmentations)
  cents <- CreateCentroids(data$centroids)
  segmentations.data <- list(
    "centroids" = cents,
    "segmentation" = segs
  )
  coords <- CreateFOV(
    coords = segmentations.data,
    type = c("segmentation", "centroids"),
    molecules = data$pixels,
    assay = assay
  )
  obj <- CreateSeuratObject(counts = data$matrix, assay = assay)

  # subset both object and coords based on the cells shared by both
  cells <- intersect(
    Cells(x = coords, boundary = "segmentation"),
    Cells(x = coords, boundary = "centroids")
  )
  cells <- intersect(Cells(obj), cells)
  coords <- subset(x = coords, cells = cells)
  obj[[fov]] <- coords
  return(obj)
}

#' @return \code{LoadVizgen}: A \code{\link[SeuratObject]{Seurat}} object
#'
#' @importFrom SeuratObject Cells CreateCentroids CreateFOV
#' CreateSegmentation CreateSeuratObject
#'
#' @export
#'
#' @rdname ReadVizgen
#'
LoadVizgen <- function(data.dir, fov, assay = 'Vizgen', z = 3L) {
  data <- ReadVizgen(
    data.dir = data.dir,
    filter = "^Blank-",
    type = c("centroids", "segmentations"),
    z = z
  )
  segs <- CreateSegmentation(data$segmentations)
  cents <- CreateCentroids(data$centroids)
  segmentations.data <- list(
    "centroids" = cents,
    "segmentation" = segs
  )
  coords <- CreateFOV(
    coords = segmentations.data,
    type = c("segmentation", "centroids"),
    molecules = data$microns,
    assay = assay
  )
  obj <- CreateSeuratObject(counts = data$transcripts, assay = assay)
  # only consider the cells we have counts and a segmentation for
  # Cells which don't have a segmentation are probably found in other z slices.
  coords <- subset(
    x = coords,
    cells = intersect(
      x = Cells(x = coords[["segmentation"]]),
      y = Cells(x = obj)
    )
  )
  # add coords to seurat object
  obj[[fov]] <- coords
  return(obj)
}

#' @return \code{LoadXenium}: A \code{\link[SeuratObject]{Seurat}} object
#'
#' @param data.dir Path to folder containing Xenium outputs
#' @param fov FOV name
#' @param assay Assay name
#' @param mols.qv.threshold Remove transcript molecules with
#' a QV less than this threshold. QV >= 20 is the standard threshold
#' used to construct the cell x gene count matrix.
#' @param cell.centroids Whether or not to load cell centroids
#' @param molecule.coordinates Whether or not to load molecule pixel coordinates
#' @param segmentations One of "cell", "nucleus" or NULL (to load either cell
#' segmentations, nucleus segmentations or neither)
#' @param flip.xy Whether or not to flip the x/y coordinates of the Xenium outputs
#' to match what is displayed in Xenium Explorer, or fit on your screen better.
#'
#' @importFrom SeuratObject Cells CreateCentroids CreateFOV
#' CreateSegmentation CreateSeuratObject CreateMolecules
#'
#' @export
#'
#' @rdname ReadXenium
#'
LoadXenium <- function(
  data.dir,
  fov = 'fov',
  assay = 'Xenium',
  mols.qv.threshold = 20,
  cell.centroids = TRUE,
  molecule.coordinates = TRUE,
  segmentations = NULL,
  flip.xy = FALSE
) {
  if(!is.null(segmentations) && !(segmentations %in% c('nucleus', 'cell'))) {
    stop('segmentations must be NULL or one of "nucleus", "cell"')
  }

  if(!cell.centroids && is.null(segmentations)) {
    stop(
      "Must load either centroids or cell/nucleus segmentations"
    )
  }

  data <- ReadXenium(
    data.dir = data.dir,
    type = c("centroids", "segmentations", "nucleus_segmentations")[
      c(cell.centroids, isTRUE(segmentations == 'cell'), isTRUE(segmentations == 'nucleus'))
    ],
    outs = c("segmentation_method", "matrix", "microns")[
      c(cell.centroids || isTRUE(segmentations != 'nucleus'), TRUE, molecule.coordinates && (cell.centroids || !is.null(segmentations)))
    ],
    mols.qv.threshold = mols.qv.threshold,
    flip.xy = flip.xy
  )

  segmentations <- intersect(c("segmentations", "nucleus_segmentations"), names(data))

  segmentations.data <- Filter(Negate(is.null), list(
    centroids = if(is.null(data$centroids)) {
      NULL
    } else {
      CreateCentroids(data$centroids)
    },
    segmentations = if(length(segmentations) > 0) {
      CreateSegmentation(
        data[[segmentations]]
      )
    } else {
      NULL
    }
  ))

  coords <- if(length(segmentations.data) > 0) {
    CreateFOV(
      segmentations.data,
      assay = assay,
      molecules = if(is.null(data$microns)) {
        NULL
      } else {
        CreateMolecules(data$microns)
      }
    )
  } else {
    NULL
  }

  slot.map <- c(
    `Blank Codeword` = 'BlankCodeword',
    `Unassigned Codeword` = 'BlankCodeword',
    `Negative Control Codeword` = 'ControlCodeword',
    `Negative Control Probe` = 'ControlProbe',
    `Genomic Control` = 'GenomicControl',
    `Protein Expression` = 'ProteinExpression'
  )

  xenium.obj <- CreateSeuratObject(counts = data$matrix[["Gene Expression"]], assay = assay)

  if(!is.null(data$metadata)) {
    Misc(xenium.obj, 'run_metadata') <- data$metadata
  }

  if(!is.null(data$segmentation_method)) {
    xenium.obj <- AddMetaData(xenium.obj, data$segmentation_method)
  }

  for(name in intersect(names(slot.map), names(data$matrix))) {
    xenium.obj[[slot.map[name]]] <- CreateAssayObject(counts = data$matrix[[name]])
  }

  xenium.obj[[fov]] <- coords
  return(xenium.obj)
}

#' @return \code{LoadAtera}: A \code{\link[SeuratObject]{Seurat}} object
#'
#' @param data.dir Path to folder containing Atera outputs
#' @param fov FOV name
#' @param assay Assay name
#' @param mols.qv.threshold Remove transcript molecules with a calibrated
#' Q-score less than this threshold
#' @param cell.centroids Whether or not to load cell centroids
#' @param molecule.coordinates Whether or not to load transcript molecule
#' coordinates
#' @param segmentations Which segmentation boundary polygons to load, one
#' of \dQuote{cell} or \dQuote{nucleus}, or \code{NULL} to load none
#' @param genes Optional character vector of gene names to restrict the
#' initial \code{molecule.coordinates} load to. When \code{molecule.coordinates
#' = TRUE}, a lazy, on-disk handle onto the transcript table is always kept
#' on the returned object (regardless of whether \code{genes} is set), so
#' additional genes can be fetched later with \code{\link{LoadAteraMolecules}}
#' without re-reading genes already fetched. See \code{\link{ReadAtera}}
#' @param feature.types Optional character vector of feature types to
#' restrict the counts matrix to (eg \dQuote{Gene Expression}); \code{NULL}
#' (the default) loads all feature types, matching current behavior.
#' \dQuote{Gene Expression} is always included, since it backs the object's
#' primary assay. See \code{\link{ReadAtera}}
#' @param bpcells.dir Path to a directory used to cache the counts matrix on
#' disk via \code{BPCells}. \code{NULL} (the default) auto-selects a
#' session-scoped cache directory so the returned object is BPCells-backed
#' and lazy by default (even for small bundles) whenever \code{BPCells} is
#' installed; pass \code{FALSE} to always get plain in-memory matrices
#' instead. See \code{\link{ReadAtera}}
#' @param morphology.image Whether to load a single channel plane of a
#' \dQuote{morphology_2d} image (eg the DAPI stain), by default at the
#' lowest-resolution pyramid level, for use as a background/overview image
#' in plots. Requires the \code{RBioFormats} package. Defaults to
#' \code{FALSE} since it is an optional dependency and not needed for
#' downstream analysis. See \code{\link{ReadAtera}}
#' @param morphology.channel Channel to load when \code{morphology.image =
#' TRUE}: either a channel name (matched case-insensitively as a substring,
#' eg \dQuote{dapi}) or a 0-based integer channel index. See
#' \code{\link{ReadAtera}}
#' @param morphology.resolution Pyramid resolution level to load when
#' \code{morphology.image = TRUE}, where \code{1} is full resolution.
#' Defaults to \code{NULL}, which loads the lowest-resolution (smallest,
#' fastest) level available. See \code{\link{ReadAtera}}
#' @param morphology.region Optional list with \code{x}/\code{y} elements
#' (each a length-2 micron range) to crop the loaded morphology image to a
#' rectangular window instead of reading the whole plane; recommended when
#' using \code{morphology.resolution = 1} (full resolution) on real-world
#' datasets, since these whole-slide images can be tens of thousands of
#' pixels per side. Defaults to \code{NULL}, which reads the whole plane.
#' See \code{\link{ReadAtera}}
#'
#' @importFrom SeuratObject Cells CreateCentroids CreateFOV CreateSegmentation
#' CreateSeuratObject CreateAssay5Object CreateMolecules
#'
#' @export
#'
#' @rdname ReadAtera
#'
LoadAtera <- function(
  data.dir,
  fov = 'fov',
  assay = 'Atera',
  mols.qv.threshold = 20,
  cell.centroids = TRUE,
  molecule.coordinates = FALSE,
  segmentations = NULL,
  genes = NULL,
  feature.types = NULL,
  bpcells.dir = NULL,
  morphology.image = FALSE,
  morphology.channel = "dapi",
  morphology.resolution = NULL,
  morphology.region = NULL
) {
  if (!is.null(segmentations) && !(segmentations %in% c('nucleus', 'cell'))) {
    stop('segmentations must be NULL or one of "nucleus", "cell"')
  }

  if (!cell.centroids && is.null(segmentations)) {
    stop("Must load either centroids or cell/nucleus segmentations")
  }

  if (!is.null(feature.types)) {
    feature.types <- union(feature.types, "Gene Expression")
  }

  data <- ReadAtera(
    data.dir = data.dir,
    outs = c("matrix", "centroids", "segmentations", "nucleus_segmentations", "morphology")[
      c(TRUE, cell.centroids, isTRUE(segmentations == 'cell'), isTRUE(segmentations == 'nucleus'), morphology.image)
    ],
    mols.qv.threshold = mols.qv.threshold,
    feature.types = feature.types,
    bpcells.dir = bpcells.dir,
    morphology.channel = morphology.channel,
    morphology.resolution = morphology.resolution,
    morphology.region = morphology.region
  )
  mols.handle <- NULL
  if (molecule.coordinates) {
    .AteraCheckDeps()
    mols.handle <- .AteraMoleculesHandle(data.dir = data.dir, mols.qv.threshold = mols.qv.threshold)
    if (!is.null(genes)) {
      data$microns <- .AteraFetchMolecules(handle = mols.handle, genes = genes)
    }
  }

  segmentations.key <- intersect(c("segmentations", "nucleus_segmentations"), names(data))

  segmentations.data <- Filter(Negate(is.null), list(
    centroids = if (is.null(data$centroids)) {
      NULL
    } else {
      CreateCentroids(data$centroids)
    },
    segmentations = if (length(segmentations.key) > 0) {
      CreateSegmentation(data[[segmentations.key]])
    } else {
      NULL
    }
  ))

  coords <- CreateFOV(
    segmentations.data,
    assay = assay,
    molecules = if (is.null(data$microns)) {
      NULL
    } else {
      CreateMolecules(data$microns)
    }
  )

  slot.map <- c(
    `Deprecated Codeword` = 'DeprecatedCodeword',
    `Unassigned Codeword` = 'BlankCodeword',
    `Negative Control Codeword` = 'ControlCodeword',
    `Negative Control Probe` = 'ControlProbe',
    `Genomic Control` = 'GenomicControl'
  )

  atera.obj <- CreateSeuratObject(counts = data$matrix[["Gene Expression"]], assay = assay)

  if (!is.null(data$metadata)) {
    Misc(atera.obj, 'run_metadata') <- data$metadata
  }

  extra.types <- intersect(names(slot.map), names(data$matrix))
  single.feature.types <- Filter(function(name) nrow(data$matrix[[name]]) == 1L, extra.types)
  for (name in setdiff(extra.types, single.feature.types)) {
    atera.obj[[slot.map[name]]] <- CreateAssay5Object(counts = data$matrix[[name]])
  }

  atera.obj <- subset(atera.obj, cells = intersect(Cells(atera.obj), Cells(coords)))
  coords <- subset(coords, cells = intersect(Cells(atera.obj), Cells(coords)))

  atera.obj[[fov]] <- coords

  for (name in single.feature.types) {
    # CreateAssay5Object errors on single-feature layers ("Layers must be
    # two-dimensional objects"), which would otherwise force falling back to
    # the legacy v3 Assay (CreateAssayObject) here -- but that class's
    # LayerData.Assay() validates its `cells` argument via
    # rlang::arg_match(..., multiple = TRUE) every time per-assay stats are
    # (re)computed (on creation, on every later subset(), etc.), and that
    # call doesn't scale: it's effectively unusable past a few thousand
    # cells, let alone the millions typical of a whole-tissue Atera run.
    # A single-feature "assay" isn't useful for standard per-feature
    # analysis anyway, so it's stashed as a plain matrix in Misc instead,
    # subset to the cells actually retained above.
    Misc(atera.obj, paste0('atera.', slot.map[name])) <- data$matrix[[name]][, Cells(atera.obj), drop = FALSE]
  }

  if (!is.null(mols.handle)) {
    Misc(atera.obj, 'atera.molecules') <- list(handle = mols.handle, cache = data$microns, fov = fov)
  }

  if (!is.null(data$morphology)) {
    Misc(atera.obj, 'atera.morphology') <- data$morphology
  }

  return(atera.obj)
}

#' Fetch additional Atera transcript molecules for genes not yet loaded
#'
#' Fetches transcript molecule coordinates for \code{genes} and adds them to
#' the molecule layer of a Seurat object created by
#' \code{\link{LoadAtera}(..., molecule.coordinates = TRUE)}, without
#' re-reading genes that have already been fetched (by this call or by
#' \code{genes} at load time). This is meant for interactive exploration,
#' where different genes are plotted one at a time and reloading the full
#' transcript table for each one would be wasteful.
#'
#' @param object A Seurat object created by \code{LoadAtera(..., molecule.coordinates = TRUE)}
#' @param genes Character vector of gene names to fetch
#'
#' @return \code{object}, with \code{genes}' transcript coordinates added to
#' its molecule layer (in addition to any genes fetched by a previous call)
#'
#' @importFrom SeuratObject CreateMolecules
#'
#' @export
#'
LoadAteraMolecules <- function(object, genes) {
  handle.info <- Misc(object, slot = 'atera.molecules')
  if (is.null(handle.info)) {
    stop(
      "'object' was not loaded with LoadAtera(..., molecule.coordinates = TRUE)",
      call. = FALSE
    )
  }

  cached.genes <- if (is.null(handle.info$cache)) character(0) else unique(handle.info$cache$gene)
  missing.genes <- setdiff(genes, cached.genes)
  if (length(missing.genes) > 0) {
    new.df <- .AteraFetchMolecules(handle = handle.info$handle, genes = missing.genes)
    handle.info$cache <- if (is.null(handle.info$cache)) {
      new.df
    } else {
      rbind(handle.info$cache, new.df)
    }
  }

  object[[handle.info$fov]]@molecules <- list(
    molecules = CreateMolecules(handle.info$cache, key = 'mols_')
  )
  Misc(object, 'atera.molecules') <- handle.info

  return(object)
}

#' Add segmentation polygons to an Atera field of view, for only the cells
#' currently present in it
#'
#' \code{LoadAtera(..., segmentations = "cell")} builds a \code{Polygons}
#' object (via \code{sp::SpatialPolygons()}) for every cell up front, which
#' doesn't scale to bundles with many hundreds of thousands of cells -- that
#' machinery gets impractically slow, and can overflow R's protection stack,
#' well before reaching cell counts typical of a whole-tissue Atera run. This
#' function instead builds segmentation polygons only for the cells already
#' present in \code{object[[fov]]}, so a bundle too large to load
#' segmentations for up front (\code{LoadAtera(..., segmentations = NULL)})
#' can still get polygons for a \code{\link[SeuratObject]{Crop}}ped
#' region-of-interest FOV, which typically has orders of magnitude fewer
#' cells.
#'
#' @param object A Seurat object created by \code{\link{LoadAtera}}
#' @param data.dir The same \code{data.dir} passed to the original
#' \code{\link{LoadAtera}} call (the Atera output bundle directory)
#' @param fov Name of the FOV (typically already cropped via
#' \code{\link[SeuratObject]{Crop}}) to add segmentation polygons to
#' @param segmentations One of \code{"cell"} or \code{"nucleus"}
#' @param max.cells Safety limit: errors rather than attempting to build
#' polygons for more than this many cells, since \code{sp::SpatialPolygons()}
#' becomes impractically slow/crash-prone well before this many
#'
#' @return \code{object}, with a \code{"segmentations"} boundary added to
#' \code{object[[fov]]}, built only for the cells currently in that FOV
#'
#' @importFrom SeuratObject Cells CreateSegmentation
#'
#' @export
#'
LoadAteraSegmentations <- function(object, data.dir, fov, segmentations = "cell", max.cells = 50000) {
  if (!isTRUE(segmentations %in% c("cell", "nucleus"))) {
    stop("segmentations must be one of \"cell\", \"nucleus\"", call. = FALSE)
  }
  cell.ids <- Cells(object[[fov]])
  if (length(cell.ids) > max.cells) {
    stop(
      "FOV '", fov, "' has ", length(cell.ids), " cells, over `max.cells` (", max.cells, "). ",
      "Building per-cell segmentation polygons (via sp::SpatialPolygons()) doesn't scale this ",
      "far and will be extremely slow or crash. Crop to a smaller region first, or raise ",
      "`max.cells` if you're sure.",
      call. = FALSE
    )
  }

  zip.file <- file.path(data.dir, "cells.zarr.zip")
  zidx <- .AteraZipIndex(zip.file)
  con <- file(zip.file, "rb")
  on.exit(close(con), add = TRUE)

  all.cell.ids <- .AteraFormatId(.AteraDecodePackedId(.AteraReadArray(con, zidx, "cell_id")))
  set.idx <- if (segmentations == "cell") 1L else 0L
  df <- .AteraReadPolygonSet(con, zidx, set.idx, all.cell.ids)
  df <- df[df$cell %in% cell.ids, , drop = FALSE]

  object[[fov]][["segmentations"]] <- CreateSegmentation(df)
  return(object)
}

#' Load per-codeword transcript breakdowns for specific genes
#'
#' Atera's per-transcript \code{codeword_identity} (which specific barcode/
#' codeword called each detection) is collapsed away by
#' \code{\link{LoadAteraMolecules}}, which reports one combined position per
#' gene. A gene can map to more than one distinct codeword (more than one
#' probe/barcode design targeting the same gene); this instead returns one
#' \code{"<gene>_codeword_<id>"} pseudo-feature per codeword, which is useful
#' for a common QC question: whether different codewords for the same gene
#' cluster differently in space (which would suggest a probe-specific
#' artifact rather than true expression).
#'
#' @param object A Seurat object created by \code{\link{LoadAtera}}
#' @param data.dir The same \code{data.dir} passed to the original
#' \code{\link{LoadAtera}} call (the Atera output bundle directory)
#' @param genes Character vector of one or more gene names to break down by
#' codeword
#' @param fov Name of the FOV to attach the resulting molecules to; defaults
#' to the object's default FOV
#'
#' @return \code{object}, with one \code{"<gene>_codeword_<id>"} molecule
#' entry per codeword mapping to each requested gene, merged into
#' \code{object[[fov]]}'s existing molecule layer (regular genes loaded via
#' \code{\link{LoadAteraMolecules}} are preserved, not overwritten; genes
#' already broken down by a previous call are not re-fetched)
#'
#' @importFrom SeuratObject DefaultFOV Misc CreateMolecules
#'
#' @export
#'
LoadAteraCodewords <- function(object, data.dir, genes, fov = NULL) {
  fov <- fov %||% DefaultFOV(object = object)

  handle <- .AteraMoleculesHandle(data.dir = data.dir)
  unknown <- setdiff(genes, handle$gene.names)
  if (length(unknown) > 0) {
    stop("Unknown Atera gene(s): ", paste(unknown, collapse = ", "), call. = FALSE)
  }

  handle.info <- Misc(object, slot = 'atera.molecules')
  existing.cache <- if (is.null(handle.info)) NULL else handle.info$cache
  cached.genes <- if (is.null(existing.cache)) character(0) else unique(existing.cache$gene)
  already.loaded <- vapply(genes, function(g) {
    any(grepl(paste0("^", g, "_codeword_"), cached.genes))
  }, logical(1L))
  missing.genes <- genes[!already.loaded]

  if (length(missing.genes) > 0) {
    zip.file <- file.path(data.dir, "transcripts.zarr.zip")
    zidx <- .AteraZipIndex(zip.file)
    nthreads <- .AteraThreads()

    gene.dfs <- lapply(missing.genes, function(gene_name) {
      gene.idx <- match(gene_name, handle$gene.names)
      # A gene's rows are typically spread across most of the grid tiles (eg
      # ~85-88 of 102 for genes checked on a real whole-transcriptome bundle);
      # each tile's read is independent, so this is parallelized across
      # tiles rather than reading them one at a time. Unlike .AteraReadArray's
      # chunk-level parallelism (which shares one connection, reading
      # sequentially before decoding in parallel), tile reads here are full,
      # separate .AteraReadArray calls -- each forked worker opens its own
      # connection, since a connection's file position can't safely be
      # shared/seeked concurrently across forked processes.
      read.tile <- function(tile) {
        go <- tile$gene_offset
        n <- go[gene.idx, 2] - go[gene.idx, 1]
        if (n == 0) {
          return(NULL)
        }
        rows <- c(go[gene.idx, 1] + 1L, go[gene.idx, 2])
        con <- file(zip.file, "rb")
        on.exit(close(con))
        loc <- .AteraReadArray(con, zidx, paste0(tile$dir, "/location"), row.range = rows)
        cw <- .AteraReadArray(con, zidx, paste0(tile$dir, "/codeword_identity"), row.range = rows)
        data.frame(x = loc[, 1], y = loc[, 2], gene = paste0(gene_name, "_codeword_", cw))
      }
      tile.dfs <- if (nthreads > 1L) {
        parallel::mclapply(handle$tiles, read.tile, mc.cores = nthreads)
      } else {
        lapply(handle$tiles, read.tile)
      }
      do.call(rbind, tile.dfs)
    })
    new.df <- do.call(rbind, gene.dfs)
    existing.cache <- if (is.null(existing.cache)) new.df else rbind(existing.cache, new.df)
  }

  object[[fov]]@molecules <- list(molecules = CreateMolecules(existing.cache, key = 'mols_'))
  if (!is.null(handle.info)) {
    handle.info$cache <- existing.cache
    Misc(object, 'atera.molecules') <- handle.info
  }

  return(object)
}

#' @param ... Extra parameters passed to \code{DimHeatmap}
#'
#' @rdname DimHeatmap
#' @concept convenience
#' @export
#'
PCHeatmap <- function(object, ...) {
  args <- list('object' = object)
  args <- c(args, list(...))
  args$reduction <- "pca"
  return(do.call(what = 'DimHeatmap', args = args))
}

#' @param ... Extra parameters passed to \code{DimPlot}
#'
#' @rdname DimPlot
#' @concept convenience
#' @export
#'
PCAPlot <- function(object, ...) {
  return(SpecificDimPlot(object = object, ...))
}

#' @rdname SpatialPlot
#' @concept convenience
#' @concept spatial
#' @export
#'
SpatialDimPlot <- function(
  object,
  group.by = NULL,
  images = NULL,
  cols = NULL,
  crop = TRUE,
  cells.highlight = NULL,
  cols.highlight = c('#DE2D26', 'grey50'),
  facet.highlight = FALSE,
  label = FALSE,
  label.size = 7,
  label.color = 'white',
  repel = FALSE,
  ncol = NULL,
  combine = TRUE,
  pt.size.factor = 1.6,
  alpha = c(1, 1),
  image.alpha = 1,
  image.scale = "lowres",
  shape = 21,
  stroke = NA,
  stroke.alpha = NA,
  label.box = TRUE,
  interactive = FALSE,
  information = NULL,
  plot_segmentations = FALSE
) {
  return(SpatialPlot(
    object = object,
    group.by = group.by,
    images = images,
    cols = cols,
    crop = crop,
    cells.highlight = cells.highlight,
    cols.highlight = cols.highlight,
    facet.highlight = facet.highlight,
    label = label,
    label.size = label.size,
    label.color = label.color,
    repel = repel,
    ncol = ncol,
    combine = combine,
    pt.size.factor = pt.size.factor,
    alpha = alpha,
    image.alpha = image.alpha,
    image.scale = image.scale,
    shape = shape,
    stroke = stroke,
    stroke.alpha = stroke.alpha,
    label.box = label.box,
    interactive = interactive,
    information = information,
    plot_segmentations = plot_segmentations
  ))
}

#' @rdname SpatialPlot
#' @concept convenience
#' @concept spatial
#' @export
#'
SpatialFeaturePlot <- function(
  object,
  features,
  images = NULL,
  crop = TRUE,
  slot = 'data',
  keep.scale = "feature",
  min.cutoff = NA,
  max.cutoff = NA,
  ncol = NULL,
  combine = TRUE,
  pt.size.factor = 1.6,
  alpha = c(1, 1),
  image.alpha = 1,
  image.scale = "lowres",
  shape = 21,
  stroke = NA,
  stroke.alpha = NA,
  interactive = FALSE,
  information = NULL,
  plot_segmentations = FALSE
) {
  return(SpatialPlot(
    object = object,
    features = features,
    images = images,
    crop = crop,
    slot = slot,
    keep.scale = keep.scale,
    min.cutoff = min.cutoff,
    max.cutoff = max.cutoff,
    ncol = ncol,
    combine = combine,
    pt.size.factor = pt.size.factor,
    alpha = alpha,
    image.alpha = image.alpha,
    image.scale = image.scale,
    shape = shape,
    stroke = stroke,
    stroke.alpha = stroke.alpha,
    interactive = interactive,
    information = information,
    plot_segmentations = plot_segmentations
  ))
}

#' @rdname DimPlot
#' @concept convenience
#' @export
#'
TSNEPlot <- function(object, ...) {
  return(SpecificDimPlot(object = object, ...))
}

#' @rdname DimPlot
#' @concept convenience
#' @export
#'
UMAPPlot <- function(object, ...) {
  return(SpecificDimPlot(object = object, ...))
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Methods for Seurat-defined generics
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Methods for R-defined generics
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Internal
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# @rdname DimPlot
#
SpecificDimPlot <- function(object, ...) {
  funs <- sys.calls()
  name <- as.character(x = funs[[length(x = funs) - 1]])[1]
  name <- tolower(x = gsub(pattern = 'Plot', replacement = '', x = name))
  args <- list('object' = object)
  args <- c(args, list(...))
  reduc <- grep(
    pattern = name,
    x = names(x = object),
    value = TRUE,
    ignore.case = TRUE
  )
  reduc <- grep(pattern = DefaultAssay(object = object), x = reduc, value = TRUE)
  args$reduction <- ifelse(test = length(x = reduc) == 1, yes = reduc, no = name)
  tryCatch(
    expr = return(do.call(what = 'DimPlot', args = args)),
    error = function(e) {
      stop(e)
    }
  )
}

#' Read output from Parse Biosciences
#'
#' @param data.dir Directory containing the data files
#' @param ... Extra parameters passed to \code{\link{ReadMtx}}
#' @concept convenience
#' @export
#'
ReadParseBio <- function(data.dir, ...) {
  file.dir <- list.files(path = data.dir, pattern = ".mtx")
  mtx <- file.path(data.dir, file.dir)
  cells <- file.path(data.dir, "cell_metadata.csv")
  features <- file.path(data.dir, "all_genes.csv")
  return(ReadMtx(
    mtx = mtx,
    cells = cells,
    features = features,
    cell.column = 1,
    feature.column = 2,
    cell.sep = ",",
    feature.sep = ",",
    skip.cell = 1,
    skip.feature = 1,
    mtx.transpose = TRUE
  ))
}

#' Read output from STARsolo
#'
#' @param data.dir Directory containing the data files
#' @param ... Extra parameters passed to \code{\link{ReadMtx}}
#'
#' @rdname ReadSTARsolo
#' @concept convenience
#' @export
#'
ReadSTARsolo <- function(data.dir, ... ) {
  mtx <- file.path(data.dir, "matrix.mtx")
  cells <- file.path(data.dir, "barcodes.tsv")
  features <- file.path(data.dir, "features.tsv")
  return(ReadMtx(mtx = mtx, cells = cells, features = features, ...))
}
