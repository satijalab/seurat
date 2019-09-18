#' @include objects.R
#' @include generics.R
#' @importFrom methods setClass setOldClass setClassUnion slot<-
#' setMethod new
#'
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Class definitions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @importFrom tibble tibble
#'
setOldClass(Classes = 'tbl_df')

SpatialAssay <- setClass(
  Class = 'SpatialAssay',
  contains = 'Assay',
  slots = list(
    'image' = 'tbl_df',
    'scale.factors' = 'list',
    'tissue.positions' = 'data.frame'
  )
)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Read10xSpatial <- function(
  outs_path,
  filtered_feature_bc_matrix_path = outs_path,
  image_path = outs_path,
  scale_factors_path = outs_path,
  tissue_positions_path = outs_path,
  filter_matrix = TRUE,
  assay = "Spatial"
) {
  a <- if (missing(x = filtered_feature_bc_matrix_path)) {
    Read10X_h5(filename = file.path(outs_path, "filtered_feature_bc_matrix.h5"))
  } else {
    Read10X_h5(filename = filtered_feature_bc_matrix_path)
  }
  if (missing(x = image_path)) {
   image <- png::readPNG(source = file.path(outs_path, "spatial", "tissue_lowres_image.jpg"))
  } else {
   image <- png::readPNG(source = image_path)
  }
  height <- nrow(x = image)
  width <- ncol(x = image)
  image <- list(grid::rasterGrob(
    image = image,
    width = unit(x = 1, units = "npc"),
    height = unit(x = 1, units = "npc")
  ))
  image <- tibble(grob = image)
  if (missing(x = scale_factors_path)) {
    scale.factors <- jsonlite::fromJSON(txt = file.path(outs_path, "spatial", "scalefactors_json.json"))
  } else {
    scale.factors <- jsonlite::fromJSON(txt = scale_factors_path)
  }
  if (missing(x = tissue_positions_path)) {
    tissue.positions <- read.csv(
      file = file.path(outs_path, "spatial", "tissue_positions_list.txt"),
      col.names = c("barcode", "tissue", "row", "col", "imagerow", "imagecol"),
      header = FALSE
    )
  } else {
    tissue.positions <- read.csv(
      file = tissue_positions_path,
      col.names = c("barcode", "tissue", "row", "col", "imagerow", "imagecol"),
      header = FALSE
    )
  }
  tissue.positions$scaled_imagerow <- tissue.positions$imagerow * scale.factors$tissue_lowres_scalef
  tissue.positions$scaled_imagecol <- tissue.positions$imagecol * scale.factors$tissue_lowres_scalef
  image$image_height <- height
  image$image_width <- width
  # Filter for only spots under tissue if TRUE (default)
  if (filter_matrix) {
    tissue.positions <- tissue.positions %>% filter(tissue == 1)
  }
  assays <- list(CreateAssayObject(
    counts = a,
    image = image,
    scale.factors = scale.factors,
    tissue.positions = tissue.positions
  ))
  names(x = assays) <- assay
  Key(object = assays[[assay]]) <- suppressWarnings(expr = UpdateKey(key = tolower(x = assay)))
  object <- new(
    Class = "Seurat",
    assays = assays,
    meta.data = data.frame(row.names = colnames(x = assays[[assay]])),
    version = packageVersion(pkg = "Seurat"),
    project.name = "SeuratProject"
  )
 DefaultAssay(object = object) <- assay
 n.calc <- CalcN(object = assays[[1]])
 if (!is.null(x = n.calc)) {
   names(x = n.calc) <- paste(names(x = n.calc), assay, sep = '_')
   object[[names(x = n.calc)]] <- n.calc
 }
 Idents(object = object) <- "SeuratProject"
 return(object)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Methods for Seurat-defined generics
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Methods for R-defined generics
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# S4 methods
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Internal
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
