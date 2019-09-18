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
    tissue.positions <- tissue.positions %>% 
		dplyr::filter(tissue == 1)
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

# for plotting the tissue image
geom_spatial <-  function(mapping = NULL,
                         data = NULL,
                         stat = "identity",
                         position = "identity",
                         na.rm = FALSE,
                         show.legend = NA,
                         inherit.aes = FALSE,
                         ...) {
  
  GeomCustom <- ggproto(
    "GeomCustom",
    Geom,
    setup_data = function(self, data, params) {
      data <- ggproto_parent(Geom, self)$setup_data(data, params)
      data
    },
    
    draw_group = function(data, panel_scales, coord) {
      vp <- grid::viewport(x=data$x, y=data$y)
      g <- grid::editGrob(data$grob[[1]], vp=vp)
      ggplot2:::ggname("geom_spatial", g)
    },
    
    required_aes = c("grob","x","y")
    
  )
  
  layer(
    geom = GeomCustom,
    mapping = mapping,
    data = data,
    stat = stat,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(na.rm = na.rm, ...)
  )
}


SpatialColors <- colorRampPalette(rev(RColorBrewer::brewer.pal(11, "Spectral")))

SingleSpatialPlot <- function(
  data,
  image.tibble,
  pt.size = NULL,
  col.by = NULL,
  na.value = 'grey50'
) {
  pt.size <- pt.size %||% AutoPointSize(data = data)
  
  if (!is.data.frame(x = data)) {
    data <- as.data.frame(x = data)
  }
  if (!is.null(x = col.by) && !col.by %in% colnames(x = data)) {
    warning("Cannot find '", col.by, "' in data, not coloring", call. = FALSE, immediate. = TRUE)
    col.by <- NULL
  }
  
  plot <- ggplot(data = data) +
    geom_spatial(data=image.tibble, aes(grob=grob), x=0.5, y=0.5)+
    geom_point(
      mapping = aes_string(
        x = colnames(x = data)[2],
        y = colnames(x = data)[1],
        color = col.by %iff% paste0("`", col.by, "`")
        ),
      size = pt.size
    ) +
    xlim(0,image.tibble$image_width)+
    ylim(image.tibble$image_height,0)+
    coord_cartesian(expand=FALSE)+
    # guides(color = guide_legend(override.aes = list(size = 3))) +
    labs(color = NULL)

  #plot <- plot + theme_cowplot()
   plot <- plot + theme_void()

  return(plot)
}

SpatialFeaturePlot <- function(
  object,
  features,
  slot = 'data',
  min.cutoff = NA,
  max.cutoff = NA,
  ncol = NULL,
  combine = TRUE
) {
 
  data <- FetchData(object = object,
                    vars = features,
                    slot = slot)
  
  features <- colnames(x = data)
  
  # # Determine cutoffs
  # min.cutoff <- mapply(
  #   FUN = function(cutoff, feature) {
  #     return(ifelse(
  #       test = is.na(x = cutoff),
  #       yes = min(data[, feature]),
  #       no = cutoff
  #     ))
  #   },
  #   cutoff = min.cutoff,
  #   feature = features
  # )
  # max.cutoff <- mapply(
  #   FUN = function(cutoff, feature) {
  #     return(ifelse(
  #       test = is.na(x = cutoff),
  #       yes = max(data[, feature]),
  #       no = cutoff
  #     ))
  #   },
  #   cutoff = max.cutoff,
  #   feature = features
  # )
  # check.lengths <- unique(x = vapply(
  #   X = list(features, min.cutoff, max.cutoff),
  #   FUN = length,
  #   FUN.VALUE = numeric(length = 1)
  # ))
  # if (length(x = check.lengths) != 1) {
  #   stop("There must be the same number of minimum and maximum cuttoffs as there are features")
  # }
  # # brewer.gran <- ifelse(
  # #   test = length(x = cols) == 1,
  # #   yes = brewer.pal.info[cols, ]$maxcolors,
  # #   no = length(x = cols)
  # # )
  # brewer.gran <- 100
  # # Apply cutoffs
  # data <- sapply(
  #   X = 1:ncol(x = data),
  #   FUN = function(index) {
  #     data.feature <- as.vector(x = data[, index])
  #     min.use <- SetQuantile(cutoff = min.cutoff[index - 3], data.feature)
  #     max.use <- SetQuantile(cutoff = max.cutoff[index - 3], data.feature)
  #     data.feature[data.feature < min.use] <- min.use
  #     data.feature[data.feature > max.use] <- max.use
  #     if (brewer.gran == 2) {
  #       return(data.feature)
  #     }
  #     data.cut <- if (all(data.feature == 0)) {
  #       0
  #     }
  #     else {
  #       as.numeric(x = as.factor(x = cut(
  #         x = as.numeric(x = data.feature),
  #         breaks = brewer.gran
  #       )))
  #     }
  #     return(data.cut)
  #   }
  # )
  colnames(x = data) <- features
  rownames(x = data) <- colnames(object)
  
  plots <- vector(
    mode = "list",
    length(features)
  )
  
  for (i in 1:length(x = features)) {
      feature <- features[i]
      plot <- SingleSpatialPlot(
        data = cbind(GetTissueCoordinates(object = object), data[, features, drop = FALSE]),
        image.tibble = GetImage(object = object),
        #pt.size = pt.size,
        col.by = feature
      )
      plot <- plot + scale_color_gradientn(colours = SpatialColors(100))
      plots[[i]] <- plot
  }
  if (combine) {
    plots <- CombinePlots(plots = plots, ncol = ncol)
  }
  return(plots)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Methods for Seurat-defined generics
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @rdname GetImage
#' @method GetImage SpatialAssay
#' @export
#'
GetImage.SpatialAssay <- function(object, ...) {
  return(slot(object = object, name = 'image'))
}

#' @rdname GetImage
#' @method GetImage Seurat
#' @export
#'
GetImage.Seurat <- function(object, assay = NULL, ...) {
  assay <- assay %||% DefaultAssay(object = object)
  if (!inherits(x = object[[assay]], what = 'SpatialAssay')) {
    stop("GetImage works only with SpatialAssay objects", call. = FALSE)
  }
  return(GetImage(object = object[[assay]]))
}

#' @rdname GetTissueCoordinates
#' @method GetTissueCoordinates SpatialAssay
#' @export
#'
GetTissueCoordinates.SpatialAssay <- function(object, ...) {
  cols.use <- c('scaled_imagerow', 'scaled_imagecol')
  return(slot(object = object, name = 'tissue.positions')[, cols.use, drop = FALSE])
}

#' @rdname GetTissueCoordinates
#' @method GetTissueCoordinates Seurat
#' @export
#'
GetTissueCoordinates.Seurat <- function(object, assay = NULL, ...) {
  assay <- assay %||% DefaultAssay(object = object)
  if (!inherits(x = object[[assay]], what = 'SpatialAssay')) {
    stop("GetTissueCoordinates works only with SpatialAssay objects", call. = FALSE)
  }
  return(GetTissueCoordinates(object = object[[assay]]))
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Methods for R-defined generics
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# S4 methods
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Internal
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
