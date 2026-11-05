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
#' @param molecule.coordinates Whether or not to load molecule coordinates
#' immediately. For Xenium, \code{FALSE} skips transcript reads.
#' @param segmentations One of "cell", "nucleus" or NULL (to load either cell
#' segmentations, nucleus segmentations or neither)
#' @param flip.xy Swap x and y in the loaded coordinates. Deprecated for
#' \code{LoadXenium}, where \code{flip.xy = TRUE} requires SeuratObject >= 5.3.0.
#' To change only the plot, use \code{flip_xy = TRUE} in
#' \code{ImageDimPlot} or \code{ImageFeaturePlot}.
#'
#' @details By default, \code{LoadXenium} keeps the original coordinates.
#' \code{ImageDimPlot} and \code{ImageFeaturePlot} display x horizontally and
#' y increasing downward, including for older objects. Use \code{flip_xy = TRUE}
#' in either plotting function to reproduce the old plot layout without
#' changing stored coordinates.
#'
#' The deprecated \code{LoadXenium(flip.xy = TRUE)} requires SeuratObject >= 5.3.0
#' to record the axis orientation. With that version or later, it swaps stored
#' coordinates and retains its previous plotting behavior.
#' To crop using stored x/y ranges, use
#' \code{\link[SeuratObject]{Crop}(..., coords = "tissue")}.
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
  flip.xy = FALSE,
  genes = NULL,
  region = NULL
) {
  if (!missing(flip.xy)) {
    deprecate_soft(
      when = "6.0.0",
      what = "LoadXenium(flip.xy = )",
      details = paste(
        "Xenium plots now use the correct orientation without needing to swap stored coordinates.",
        "For the old layout, set `flip_xy = TRUE` in `ImageDimPlot()` or `ImageFeaturePlot()`."
      )
    )
  }
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
    outs = c("segmentation_method", "matrix")[
      c(cell.centroids || isTRUE(segmentations != 'nucleus'), TRUE)
    ],
    mols.qv.threshold = mols.qv.threshold,
    flip.xy = flip.xy,
    genes = genes,
    region = region
  )

  slot.map <- c(
    `Blank Codeword` = 'BlankCodeword',
    `Unassigned Codeword` = 'BlankCodeword',
    `Negative Control Codeword` = 'ControlCodeword',
    `Negative Control Probe` = 'ControlProbe',
    `Genomic Control` = 'GenomicControl',
    `Protein Expression` = 'ProteinExpression'
  )

  xenium.obj <- CreateSeuratObject(counts = data$matrix[["Gene Expression"]], assay = assay)
  # Match transcript labels to the feature names assigned to the expression assay
  feature.map <- stats::setNames(rownames(xenium.obj[[assay]]),
    rownames(data$matrix[["Gene Expression"]]))

  if(!is.null(data$metadata)) {
    Misc(xenium.obj, 'run_metadata') <- data$metadata
  }

  if(!is.null(data$segmentation_method)) {
    xenium.obj <- AddMetaData(xenium.obj, data$segmentation_method)
  }

  for(name in intersect(names(slot.map), names(data$matrix))) {
    control.assay <- CreateAssayObject(counts = data$matrix[[name]])
    feature.map <- c(feature.map, stats::setNames(
      rownames(control.assay), rownames(data$matrix[[name]])))
    xenium.obj[[slot.map[name]]] <- control.assay
  }

  # Transcript filters still use names from the source files.
  feature.map <- feature.map[!duplicated(names(feature.map))]
  source.genes <- .SpatialSourceGenes(genes, feature.map)
  if (molecule.coordinates) {
    data$microns <- .SpatialMapMolecules(.XeniumReadMolecules(data.dir,
      genes = source.genes, mols.qv.threshold = mols.qv.threshold,
      flip.xy = flip.xy, region = region), feature.map)
  }

  segmentations.key <- intersect(c("segmentations", "nucleus_segmentations"), names(data))

  segmentations.data <- Filter(Negate(is.null), list(
    centroids = if(is.null(data$centroids)) {
      NULL
    } else {
      CreateCentroids(data$centroids)
    },
    segmentations = if(length(segmentations.key) > 0) {
      CreateSegmentation(data[[segmentations.key]])
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

  # preserve legacy plotting only for explicitly flipped Xenium loads
  # this ensures that plotting functions interpret the flipped coordinates correctly
  # should be removed once flip.xy is fully deprecated
  if (!is.null(coords) && flip.xy) {
    if (!.hasSlot(coords, "coords_x_orientation")) {
      stop("LoadXenium(flip.xy = TRUE) requires SeuratObject >= 5.3.0. Please update SeuratObject and restart R.",
        call. = FALSE)
    }
    slot(coords, "coords_x_orientation") <- "vertical"
  }

  xenium.obj[[fov]] <- coords
  return(xenium.obj)
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

# Map requested assay names back to source names for transcript reads.
# feature.map maps source names to assay names; source names are also accepted.
.SpatialSourceGenes <- function(genes, feature.map) {
  if (is.null(genes)) return(NULL)
  if (!is.character(genes) || anyNA(genes)) stop("genes must be a character vector without missing values")
  if (!length(genes)) return(character())
  source.names <- names(feature.map)
  source.genes <- lapply(genes, function(gene) {
                      if (gene %in% source.names) return(gene)
                      source <- source.names[feature.map == gene]
                      if (length(source)) source else gene
                  })
  return(unique(unlist(source.genes, use.names = FALSE)))
}

# Map molecule labels to assay feature names, leaving unmatched names unchanged.
.SpatialMapMolecules <- function(molecules, feature.map) {
  index <- match(molecules$gene, names(feature.map))
  matched <- !is.na(index)
  molecules$gene[matched] <- unname(feature.map[index[matched]])
  molecules
}
