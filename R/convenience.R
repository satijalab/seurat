#' @include generics.R
#' @include visualization.R
#'
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @param fov Name to store FOV as
#' @param assay Name to store expression matrix as
#' @inheritDotParams ReadAkoya
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
#' @param fov Name to store FOV as
#' @param assay Name to store expression matrix as
#' @param add.zIndex If to add \code{z} slice index to a cell
#' @param update.object If to update final object, default to TRUE
#' @param add.molecules If to add \code{molecules} coordinates to FOV of the object, 
#'  default to TRUE
#' @param ... Arguments passed to \code{ReadVizgen}
#'
#' @importFrom SeuratObject Cells CreateCentroids CreateFOV
#' CreateSegmentation CreateSeuratObject
#' @import dplyr
#'
#' @export
#'
#' @rdname ReadVizgen
#'
LoadVizgen <- function(
    data.dir, 
    fov = 'vz', 
    assay = 'Vizgen',
    mol.type = 'microns',
    filter = '^Blank-',
    z = 3L,
    add.zIndex = TRUE, 
    update.object = TRUE,
    add.molecules = TRUE,
    min.area = 5,
    verbose,
    ...)
{
  # reading data..
  data <- ReadVizgen(data.dir = data.dir,
                     mol.type = mol.type,
                     filter = filter,
                     z = z,
                     min.area = min.area,
                     verbose = verbose,
                     ...)
  
  if (verbose) { message("Creating Seurat object..") }  
  obj <- CreateSeuratObject(counts = data[["transcripts"]], assay = assay)
  
  # in case no segmentation is present, use boxes
  if (!"segmentations" %in% names(data)) {
    if ("boxes" %in% names(data)) {
      bound.boxes <- CreateSegmentation(data[["boxes"]])
      cents <- CreateCentroids(data[["centroids"]])
      bound.boxes.data <- list(centroids = cents, 
                               boxes = bound.boxes)
      if (verbose) { 
        message("Creating FOVs..", "\n", 
                if (!add.molecules) { ">>> `molecules` coordidates will be skipped" }, 
                "\n",
                ">>> using box coordinates instead of segmentations") 
      }
      coords <- 
        CreateFOV(coords = bound.boxes.data,
                  type = c("boxes", "centroids"),
                  molecules = 
                    if (add.molecules) {
                      data[[mol.type]] } else { NULL }, 
                  assay = assay) %>%
        subset(x = .,
               cells = intersect(x = Cells(x = .[["boxes"]]),
                                 y = Cells(x = obj)))
    } else { 
      # in case no segmentation & no boxes are present, use centroids only
      cents <- CreateCentroids(data[["centroids"]])
      if (verbose) { 
        message("Creating FOVs..", "\n", 
                if (!add.molecules) { ">>> `molecules` coordidates will be skipped" }, 
                "\n", 
                ">>> using only centroids") 
      }
      coords <-
        CreateFOV(coords = list(centroids = cents),
                  type = c("centroids"),
                  molecules = 
                    if (add.molecules) {
                      data[[mol.type]] } else { NULL }, 
                  assay = assay) %>%
        subset(x = ., 
               cells = intersect(x = Cells(x = .[["centroids"]]),
                                 y = Cells(x = obj)))
    }
  } else if ("segmentations" %in% names(data)) {
    segs <- CreateSegmentation(data[["segmentations"]])
    cents <- CreateCentroids(data[["centroids"]])
    segmentations.data <- list(centroids = cents, segmentation = segs)
    if (verbose) { 
      message("Creating FOVs..", "\n", 
              if (!add.molecules) { ">>> `molecules` coordidates will be skipped" }, 
              "\n", 
              ">>> using segmentations") 
    }
    coords <-
      CreateFOV(coords = segmentations.data, 
                type = c("segmentation", "centroids"), 
                molecules = 
                  if (add.molecules) {
                    data[[mol.type]] } else { NULL }, 
                assay = assay) %>%
      # only consider the cells we have counts and a segmentation.
      # Cells which don't have a segmentation are probably found in other z slices.
      subset(x = .,
             cells = intersect(x = Cells(x = .[["segmentation"]]),
                               y = Cells(x = obj)))
  }
  
  # add z-stack index for cells
  if (add.zIndex) { obj$z <- data$zIndex %>% pull(z) }
  
  # add metadata vars
  if (verbose) { message(">>> adding metadata infos") }
  if (c("metadata" %in% names(data))) {
    metadata <- match.arg(arg = "metadata", choices = names(data), several.ok = TRUE)
    meta.vars <- names(data[[metadata]])
    for (i in meta.vars %>% seq) {
      obj %<>% AddMetaData(metadata = data[[metadata]][[meta.vars[i]]], 
                           col.name = meta.vars[i])
    }
  }
  
  # sanity on fov name
  fov %<>% gsub("_|-", ".", .)
  
  if (verbose) { message(">>> adding FOV") }
  obj[[fov]] <- coords
  
  ## filter - keep cells with counts > 0
  # helper function to return metadata
  callmeta <- function (object = NULL) { return(object@meta.data) }
  nCount <- grep("nCount", callmeta(obj) %>% names, value = TRUE)
  if (any(obj[[nCount]] == 0)) {
    if (verbose) { message(">>> filtering object - keeping cells with counts > 0") }
    obj %<>% subset(subset = !!base::as.symbol(nCount) > 0)
  } else { if (verbose) { message(">>> all counts are > 0") } }
  
  if (update.object) { 
    if (verbose) { message("Updating object:") 
      obj %<>% UpdateSeuratObject()
    } else { 
      obj %<>% 
        UpdateSeuratObject() %>% 
        suppressMessages() } }
  
  if (verbose) { message("Object is ready!") } 
  return(obj)
}

#' @return \code{LoadXenium}: A \code{\link[SeuratObject]{Seurat}} object
#'
#' @param data.dir Path to folder containing Nanostring SMI outputs
#' @param fov FOV name
#' @param assay Assay name
#'
#' @importFrom SeuratObject Cells CreateCentroids CreateFOV
#' CreateSegmentation CreateSeuratObject
#'
#' @export
#'
#' @rdname ReadXenium
#'
LoadXenium <- function(data.dir, fov = 'fov', assay = 'Xenium') {
  data <- ReadXenium(
    data.dir = data.dir,
    type = c("centroids", "segmentations"),
  )
  
segmentations.data <- list(
    "centroids" = CreateCentroids(data$centroids),
    "segmentation" = CreateSegmentation(data$segmentations)
  )
  coords <- CreateFOV(
    coords = segmentations.data,
    type = c("segmentation", "centroids"),
    molecules = data$microns,
    assay = assay
  )
  
xenium.obj <- CreateSeuratObject(counts = data$matrix[["Gene Expression"]], assay = assay)
  if("Blank Codeword" %in% names(data$matrix))
    xenium.obj[["BlankCodeword"]] <- CreateAssayObject(counts = data$matrix[["Blank Codeword"]])
  else
    xenium.obj[["BlankCodeword"]] <- CreateAssayObject(counts = data$matrix[["Unassigned Codeword"]])
  xenium.obj[["ControlCodeword"]] <- CreateAssayObject(counts = data$matrix[["Negative Control Codeword"]])
  xenium.obj[["ControlProbe"]] <- CreateAssayObject(counts = data$matrix[["Negative Control Probe"]])
  
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
  stroke = 0.25,
  label.box = TRUE,
  interactive = FALSE,
  information = NULL
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
    stroke = stroke,
    label.box = label.box,
    interactive = interactive,
    information = information
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
  stroke = 0.25,
  interactive = FALSE,
  information = NULL
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
    stroke = stroke,
    interactive = interactive,
    information = information
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
  mtx <- file.path(data.dir, "DGE.mtx")
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
