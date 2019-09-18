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

Read10xSpatial <- function(outs_path, filtered_feature_bc_matrix_path = outs_path, image_path=outs_path, scale_factors_path=outs_path, tissue_positions_path=outs_path, filter_matrix=TRUE, assay = "Spatial") {
  
  if(missing(filtered_feature_bc_matrix_path)) {
    a <- Seurat::Read10X_h5(paste(outs_path, "/filtered_feature_bc_matrix.h5", sep = ""))
  } else {
    a <- Seurat::Read10X_h5(filtered_feature_bc_matrix_path)
  }
    
  if(missing(image_path)) {
   image <- png::readPNG(paste(outs_path, "/spatial/tissue_lowres_image.jpg", sep = ""))
  } else {
   image <- png::readPNG(image_path)
    }
    height <- nrow(image)
    width <- ncol(image)
    image <- list(grid::rasterGrob(image, width=unit(1,"npc"), height=unit(1,"npc")))
    image <- tibble(grob=image)
  
  if(missing(scale_factors_path)) {
    scale.factors <- jsonlite::fromJSON(txt = paste(outs_path, "/spatial/scalefactors_json.json", sep = ""))
  } else {
    scale.factors <- jsonlite::fromJSON(txt = scale_factors_path)
  }
    
  if(missing(tissue_positions_path)) {  
    tissue.positions <- read.csv(paste(outs_path,"/spatial/tissue_positions_list.txt", sep = ""),col.names=c("barcode","tissue","row","col","imagerow","imagecol"), header = FALSE)
  } else {
    tissue.positions <- read.csv(tissue_positions_path,col.names=c("barcode","tissue","row","col","imagerow","imagecol"), header = FALSE)
  }
    
    tissue.positions$scaled_imagerow <- tissue.positions$imagerow * scale.factors$tissue_lowres_scalef
    tissue.positions$scaled_imagecol <- tissue.positions$imagecol * scale.factors$tissue_lowres_scalef
    
   image$image_height <- height
   image$image_width <- width
    
    # Filter for only spots under tissue if TRUE (default)
   if(filter_matrix) {
    tissue.positions <- tissue.positions %>%
      filter(tissue == 1)
   }
  
   
  CreateAssayObject(counts = a, image = image, scale.factors = scale.factors, tissue.positions = tissue.positions)
  
  assays <- list(CreateAssayObject(counts = a, image = image, scale.factors = scale.factors, tissue.positions = tissue.positions))
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
