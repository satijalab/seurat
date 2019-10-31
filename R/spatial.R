#' @include objects.R
#' @include generics.R
#' @importFrom methods setClass setOldClass slot<- setAs setMethod new
#'
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Class definitions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @note \code{scalefactors} objects can be created with \code{scalefactors()}
#'
#' @param spot Spot full resolution scale factor
#' @param fiducial Fiducial full resolution scale factor
#' @param hires High resolutoin scale factor
#' @param lowres Low resolution scale factor
#'
#' @rdname ScaleFactors
#' @export
#'
scalefactors <- function(spot, fiducial, hires, lowres) {
  object <- list(
    spot = spot,
    fiducial = fiducial,
    hires = hires,
    lowres = lowres
  )
  object <- sapply(X = object, FUN = as.numeric, simplify = FALSE, USE.NAMES = TRUE)
  return(structure(.Data = object, class = 'scalefactors'))
}

setOldClass(Classes = c('scalefactors'))

#' The SpatialImage class
#'
#' The SpatialImage class is a virtual class representing spatial information for
#' Seurat. All spatial image information must inherit from this class for use with
#' \code{Seurat} objects
#'
#' @slot assay Name of assay to associate image data with; will give this image
#' priority for visualization when the assay is set as the active/default assay
#' in a \code{Seurat} object
#'
#' @section Provided methods:
#' These methods are defined on the \code{SpatialImage} object and should not be
#' overwritten without careful thought
#' \itemize{
#'   \item \code{\link{DefaultAssay}} and \code{\link{DefaultAssay<-}}
#'   \item \code{\link{IsGlobal}}
#' }
#'
#' @section Required methods:
#' All subclasses of the \code{SpatialImage} class must define the following methods;
#' simply relying on the \code{SpatialImage} method will result in errors
#' \describe{
#'   \item{\code{\link{Cells}}}{Return the cell/spot barcodes associated with each position}
#'   \item{\code{\link{dim}}}{...}
#'   \item{\code{\link{GetImage}}}{Return image data; by default, must return a grob object}
#'   \item{\code{\link{GetTissueCoordinates}}}{Return tissue coordinates; by default,
#'   must return a two-column data.frame with x-coordinates in the first column and y-coordiantes
#'   in the second}
#'   \item{\code{\link{RenameCells}}}{Rename the cell/spot barcodes for this image}
#'   \item{\code{\link{subset}} and \code{[}}{Subset the image data by cells/spots;
#'   \code{[} should only take \code{i} for subsetting by cells/spots}
#' }
#' These methods are used throughout Seurat, so defining them and setting the proper
#' defaults will allow subclasses of \code{SpatialImage} to work seamlessly
#'
#' @name SpatialImage-class
#' @rdname SpatialImage-class
#' @exportClass SpatialImage
#'
SpatialImage <- setClass(
  Class = 'SpatialImage',
  contains = 'VIRTUAL',
  slots = list(
    'assay' = 'character'
  )
)

#' The SliceImage class
#'
#' The SliceImage class represents spatial information from the 10X Genomics Visium
#' platform
#'
#' @slot image A three-dimensional array with PNG image data, see
#' \code{\link[png]{readPNG}} for more details
#' @slot scale.factors An object of class \code{\link{scalefactors}}; see
#' \code{\link{scalefactors}} for more information
#' @slot coordinates A data frame with tissue coordinate information
#'
#' @name SliceImage-class
#' @rdname SliceImage-class
#' @exportClass SliceImage
#'
SliceImage <- setClass(
  Class = 'SliceImage',
  contains = 'SpatialImage',
  slots = list(
    'image' = 'array',
    'scale.factors' = 'scalefactors',
    'coordinates' = 'data.frame'
  )
)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Pull spatial image names
#'
#' List the names of \code{SpatialImage} objects present in a \code{Seurat} object.
#' If \code{assay} is provided, limits search to images associated with that assay
#'
#' @param object A \code{Seurat} object
#' @param assay Name of assay to limit search to
#'
#' @return A list of image names
#'
#' @export
#'
Images <- function(object, assay = NULL) {
  object <- UpdateSlots(object = object)
  images <- names(x = slot(object = object, name = 'images'))
  if (!is.null(x = assay)) {
    assays <- c(assay, DefaultAssay(object = object[[assay]]))
    images <- Filter(
      f = function(x) {
        return(DefaultAssay(object = object[[x]]) %in% assays)
      },
      x = images
    )
  }
  return(images)
}

#' Load a 10X Genomics Visium Image
#'
#' @inheritParams Read10X
#' @param ... Ignored for now
#'
#' @return A \code{\link{SliceImage}} object
#'
#' @importFrom png readPNG
#' @importFrom grid rasterGrob
#' @importFrom jsonlite fromJSON
#'
#' @seealso \code{\link{SliceImage}} \code{\link{Load10X_Spatial}}
#'
#' @export
#'
Read10X_Image <- function(data.dir, filter.matrix = TRUE, ...) {
  image <- readPNG(source = file.path(data.dir, 'spatial', 'tissue_lowres_image.png'))
  scale.factors <- fromJSON(txt = file.path(data.dir, 'spatial', 'scalefactors_json.json'))
  tissue.positions <- read.csv(
    file = file.path(data.dir, 'spatial', 'tissue_positions_list.txt'),
    col.names = c('barcodes', 'tissue', 'row', 'col', 'imagerow', 'imagecol'),
    header = FALSE,
    as.is = TRUE,
    row.names = 1
  )
  if (filter.matrix) {
    tissue.positions <- tissue.positions[which(x = tissue.positions$tissue == 1), , drop = FALSE]
  }
  return(new(
    Class = 'SliceImage',
    image = image,
    scale.factors = scalefactors(
      spot = scale.factors$tissue_hires_scalef,
      fiducial = scale.factors$fiducial_diameter_fullres,
      hires = scale.factors$tissue_hires_scalef,
      scale.factors$tissue_lowres_scalef
    ),
    coordinates = tissue.positions
  ))
}


#' @importFrom crosstalk SharedData bscols
#' @importFrom plotly layout plot_ly highlight
#'
LinkPlots <- function(
  plot1,
  plot2,
  plot1.labels = TRUE,
  plot2.labels = TRUE,
  information = NULL,
  plot1.layout = list(),
  plot2.layout = list(),
  ...
) {
  # Get plot builds for each plot
  plot1.build <- GGpointToPlotlyBuild(plot = plot1, information = information)
  plot1.labels <- if (plot1.labels) {
    GetXYAesthetics(plot = plot1)
  } else {
    list()
  }
  plot2.build <- GGpointToPlotlyBuild(plot = plot2, information = information)
  plot2.labels <- if (plot2.labels) {
    GetXYAesthetics(plot = plot2)
  } else {
    list()
  }
  # Generate the full build
  plot.build <- merge(x = plot1.build, plot2.build, by = 0)
  rownames(x = plot.build) <- plot.build$Row.names
  plot.build <- plot.build[, which(x = colnames(x = plot.build) != 'Row.names'), drop = FALSE]
  plot.build <- SharedData$new(data = plot.build)
  # Build the plotly plots
  Axis <- function(title = NULL) {
    return(list(
      title = title %||% '',
      showgrid = FALSE,
      zeroline = FALSE,
      showline = TRUE
    ))
  }
  plot1.layout$xaxis <- c(Axis(title = plot1.labels[['x']]), plot1.layout$xaxis)
  plot1.layout$yaxis <- c(Axis(title = plot1.labels[['y']]), plot1.layout$yaxis)
  plot1.layout <- c(
    list(p = plot_ly(
      data = plot.build,
      x = ~x.x,
      y = ~y.x,
      type = 'scatter',
      mode = 'markers',
      color = ~I(color.x),
      # size = ~I(cex.x),
      # symbol = ~I(pch.x),
      hoverinfo = 'text',
      text = ~feature.x
    )),
    plot1.layout
  )
  plot1.plotly <- do.call(what = 'layout', args = plot1.layout)
  plot1.plotly <- highlight(p = plot1.plotly, on = 'plotly_selected')
  plot2.layout$xaxis <- c(Axis(title = plot2.labels[['x']]), plot2.layout$xaxis)
  plot2.layout$yaxis <- c(Axis(title = plot2.labels[['y']]), plot2.layout$yaxis)
  plot2.layout <- c(
    list(p = plot_ly(
      data = plot.build,
      x = ~x.y,
      y = ~y.y,
      type = 'scatter',
      mode = 'markers',
      color = ~I(color.y),
      # size = ~I(cex.y),
      # symbol = ~I(pch.y),
      hoverinfo = 'text',
      text = ~feature.y
    )),
    plot2.layout
  )
  plot2.plotly <- do.call(what = 'layout', args = plot2.layout)
  plot2.plotly <- highlight(p = plot2.plotly, on = 'plotly_selected')
  return(bscols(plot1.plotly, plot2.plotly))
}

#' @importFrom plotly raster2uri
#'
LinkedFeaturePlot <- function(
  object,
  feature,
  dims = 1:2,
  reduction = NULL,
  image = NULL,
  slot = 'data',
  min.cutoff = NA,
  max.cutoff = NA
) {
  if (length(x = feature) > 1) {
    stop("'LinkedFeaturePlot' currently only supports one feature", call. = FALSE)
  }
  expression.data <- FetchData(object = object, vars = feature)
  image <- image %||% DefaultImage(object = object)
  image.plotly <- list(
    source = raster2uri(r = GetImage(object = object[[image]], mode = 'raster')),
    xref = 'x',
    yref = 'y',
    x = -7,
    y = -7,
    sizex = ncol(x = object[[image]]),
    sizey = nrow(x = object[[image]]),
    sizing = 'stretch',
    opacity = 1,
    layer = 'below'
  )
  spatial.plot <- SpatialFeaturePlot(
    object = object,
    features = feature,
    images = image,
    slot = 'data',
    min.cutoff = min.cutoff,
    max.cutoff = max.cutoff
  )
  feature.plot <- FeaturePlot(
    object = object,
    features = feature,
    dims = dims,
    reduction = reduction,
    min.cutoff = min.cutoff,
    max.cutoff = max.cutoff
  )
  return(LinkPlots(
    plot1 = spatial.plot,
    plot2 = feature.plot,
    plot1.labels = FALSE,
    information = expression.data,
    plot1.layout = list('images' = image.plotly)
  ))
}

#' @importFrom plotly raster2uri
#'
LinkedDimPlot <- function(
  object,
  dims = 1:2,
  reduction = NULL,
  image = NULL,
  group.by = NULL
) {
  image <- image %||% DefaultImage(object = object)
  image.plotly <- list(
    source = raster2uri(r = GetImage(object = object[[image]], mode = 'raster')),
    xref = 'x',
    yref = 'y',
    x = -7,
    y = -7,
    sizex = ncol(x = object[[image]]),
    sizey = nrow(x = object[[image]]),
    sizing = 'stretch',
    opacity = 1,
    layer = 'below'
  )
  spatial.plot <- SpatialDimPlot(
    object = object,
    group.by = group.by,
    images = image
  )
  dim.plot <- DimPlot(
    object = object,
    group.by = group.by,
    dims = dims,
    reduction = reduction
  )
  return(LinkPlots(
    plot1 = spatial.plot,
    plot2 = dim.plot,
    plot1.labels = FALSE,
    plot1.layout = list(
      'images' = image.plotly
      # 'xaxis' = list('range' = c(0, ncol(x = object[[image]]))),
      # 'yaxis' = list('range' = c(nrow(x = object[[image]]), 0))
    )
  ))
}

#' Load a 10X Genomics Visium Spatial Experiment into a \code{Seurat} object
#'
#' @inheritParams Read10X
#' @param ... Arguments passed to \code{\link{Read10X_h5}}
#'
#' @return A \code{Seurat} object
#'
#' @importFrom png readPNG
#' @importFrom grid rasterGrob
#' @importFrom jsonlite fromJSON
#'
#' @export
#'
Load10X_Spatial <- function(
  data.dir,
  assay = 'Spatial',
  slice = 'slice1',
  filter.matrix = TRUE,
  to.upper = FALSE,
  ...
) {
  filename <- file.path(data.dir, 'filtered_feature_bc_matrix.h5')
  data <- Read10X_h5(filename = filename, ...)
  if (to.upper) {
    rownames(x = data) <- toupper(x = rownames(x = data))
  }
  object <- CreateSeuratObject(counts = data, assay = assay)
  image <- Read10X_Image(data.dir = data.dir, filter.matrix = filter.matrix)
  image <- image[Cells(x = object)]
  DefaultAssay(object = image) <- assay
  object[[slice]] <- image
  return(object)
}

# for plotting the tissue image
#' @importFrom ggplot2 ggproto Geom ggproto_parent layer
#'
geom_spatial <-  function(
  mapping = NULL,
  data = NULL,
  stat = "identity",
  position = "identity",
  na.rm = FALSE,
  show.legend = NA,
  inherit.aes = FALSE,
  ...
) {
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

#' @importFrom grDevices colorRampPalette
#'
SpatialColors <- colorRampPalette(rev(RColorBrewer::brewer.pal(11, "Spectral")))

#' @importFrom ggplot2 ggplot geom_point aes_string xlim ylim
#' coord_cartesian labs theme_void
#'
SingleSpatialPlot <- function(
  data,
  image.tibble,
  pt.size = NULL,
  alpha = 1,
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
    geom_spatial(data = image.tibble, aes_string(grob = 'grob'), x = 0.5, y = 0.5) +
    geom_point(
      mapping = aes_string(
        x = colnames(x = data)[2],
        y = colnames(x = data)[1],
        color = col.by %iff% paste0("`", col.by, "`")
      ),
      size = pt.size,
      alpha = alpha
    ) +
    xlim(0, image.tibble$image_width) +
    ylim(image.tibble$image_height, 0) +
    coord_cartesian(expand = FALSE) +
    labs(color = NULL)
  plot <- plot + theme_void()
  return(plot)
}

#' @importFrom tibble tibble
#' @importFrom ggplot2 theme element_text
#'
SpatialDimPlot <- function(
  object,
  assay = NULL,
  group.by = NULL,
  images = NULL,
  pt.size = NULL,
  alpha = 1,
  combine = TRUE,
  ...
) {
  assay <- assay %||% DefaultAssay(object = object)
  DefaultAssay(object = object) <- assay
  object[['ident']] <- Idents(object = object)
  group.by <- group.by %||% 'ident'
  data <- object[[group.by]]
  for (group in group.by) {
    if (!is.factor(x = data[, group])) {
      data[, group] <- factor(x = data[, group])
    }
  }
  # TODO: Replace this once SCTransform assigns assay.orig
  # images <- images %||% Images(object = object, assay = DefaultAssay(object = object))
  images <- images %||% DefaultImage(object = object)
  plots <- vector(
    mode = "list",
    length = length(x = group.by)
  )
  for (i in 1:length(x = group.by)) {
    group <- group.by[i]
    slice.plots <- vector(mode = "list",length = length(x = images))
    for (j in 1:length(x = images)) {
      image.use <- object[[images[[j]]]]
      coordinates <- GetTissueCoordinates(object = image.use)
      plot <- SingleSpatialPlot(
        data = cbind(
          coordinates,
          data[rownames(x = coordinates), group, drop = FALSE]
        ),
        col.by = group,
        image.tibble = tibble(
          grob = list(GetImage(object = image.use)),
          image_width = ncol(x = image.use),
          image_height = nrow(x = image.use)
        ),
        pt.size = pt.size,
        alpha = alpha
      )
      if (i == 1) {
        plot <- plot + ggtitle(label = images[j])  + theme(plot.title = element_text(hjust = 0.5))
      }
      slice.plots[[j]] <- plot
    }
    plots[[i]] <- CombinePlots(plots = slice.plots, ncol = j)
  }
  if (combine) {
    plots <- CombinePlots(plots = plots, ncol = 1)
  }
  return(plots)
}

#' @importFrom tibble tibble
#' @importFrom ggplot2 scale_color_gradientn ggtitle theme element_text
#'
SpatialFeaturePlot <- function(
  object,
  features,
  images = NULL,
  slot = 'data',
  min.cutoff = NA,
  max.cutoff = NA,
  ncol = NULL,
  combine = TRUE
) {
  data <- FetchData(
    object = object,
    vars = features,
    slot = slot
  )
  # TODO: Replace this once SCTransform assigns assay.orig
  # images <- images %||% Images(object = object, assay = DefaultAssay(object = object))
  images <- images %||% DefaultImage(object = object)
  if (length(x = images) < 1) {
    stop("Could not find any spatial image information")
  }
  features <- colnames(x = data)
  colnames(x = data) <- features
  rownames(x = data) <- colnames(object)
  plots <- vector(
    mode = "list",
    length = length(x = features)
  )
  for (i in 1:length(x = features)) {
    feature <- features[i]
    slice.plots <- vector(mode = "list",length = length(x = images))
    for (j in 1:length(x = images)) {
      image.use <- object[[images[[j]]]]
      coordinates <- GetTissueCoordinates(object = image.use)
      plot <- SingleSpatialPlot(
        data = cbind(
          coordinates,
          data[rownames(x = coordinates), features, drop = FALSE]
        ),
        image.tibble = tibble(
          grob = list(GetImage(object = image.use)),
          image_width = ncol(x = image.use),
          image_height = nrow(x = image.use)
        ),
        col.by = feature
      )
      plot <- plot + scale_color_gradientn(name = feature, colours = SpatialColors(100))
      if (i == 1) {
        plot <- plot + ggtitle(label = images[[j]]) + theme(plot.title = element_text(hjust = 0.5))
      }
      slice.plots[[j]] <- plot
    }
    plots[[i]] <- CombinePlots(plots = slice.plots, ncol = j)
  }
  if (combine) {
    plots <- CombinePlots(plots = plots, ncol = 1)
  }
  return(plots)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Methods for Seurat-defined generics
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @rdname Cells
#' @method Cells SliceImage
#' @export
#'
Cells.SliceImage <- function(x) {
  return(rownames(x = GetTissueCoordinates(object = x, scale = NULL)))
}

#' @rdname DefaultAssay
#' @method DefaultAssay SpatialImage
#' @export
#'
DefaultAssay.SpatialImage <- function(object, ...) {
  CheckDots(...)
  return(slot(object = object, name = 'assay'))
}

#' @rdname DefaultAssay
#' @method DefaultAssay<- SpatialImage
#' @export
#'
"DefaultAssay<-.SpatialImage" <- function(object, ..., value) {
  CheckDots(...)
  slot(object = object, name = 'assay') <- value
  return(object)
}

#' @param image Name of \code{SpatialImage} object to pull image data for; if
#' \code{NULL}, will attempt to select an image automatically
#'
#' @rdname GetImage
#' @method GetImage Seurat
#' @export
#'
GetImage.Seurat <- function(object, mode = 'grob', image = NULL, ...) {
  image <- image %||% DefaultImage(object = object)
  if (is.null(x = image)) {
    stop("No images present in this Seurat object", call. = FALSE)
  }
  return(GetImage(object = object[[image]], mode = mode, ...))
}

#' @importFrom grDevices as.raster
#' @importFrom grid rasterGrob unit
#'
#' @rdname GetImage
#' @method GetImage SliceImage
#' @export
#'
GetImage.SliceImage <- function(object, mode = 'grob', ...) {
  image <- slot(object = object, name = 'image')
  image <- switch(
    EXPR = mode,
    'grob' = rasterGrob(
      image = image,
      width = unit(x = 1, units = 'npc'),
      height = unit(x = 1, units = 'npc')
    ),
    'raster' = as.raster(x = image),
    image
  )
  return(image)
}

#' @rdname GetImage
#' @method GetImage SpatialImage
#' @export
#'
GetImage.SpatialImage <- function(object, mode = 'grob', ...) {
  stop(
    "'GetImage' must be overridden for all sublcasses of 'SpatialImage'",
    call. = FALSE
  )
}

#' @param image Name of \code{SpatialImage} object to get coordinates for; if
#' \code{NULL}, will attempt to select an image automatically
#'
#' @rdname GetTissueCoordinates
#' @method GetTissueCoordinates Seurat
#' @export
#'
GetTissueCoordinates.Seurat <- function(object, image = NULL, ...) {
  image <- image %||% DefaultImage(object = object)
  if (is.null(x = image)) {
    stop("No images present in this Seurat object", call. = FALSE)
  }
  return(GetTissueCoordinates(object = object[[image]], ...))
}

#' @param scale A factor to scale the coordinates by; choose from: 'tissue',
#' 'fiducial', 'hires', 'lowres', or \code{NULL} for no scaling
#' @param cols Columns of tissue coordinates data.frame to pull
#'
#' @rdname GetTissueCoordinates
#' @method GetTissueCoordinates SliceImage
#' @export
#'
GetTissueCoordinates.SliceImage <- function(
  object,
  scale = 'lowres',
  cols = c('imagerow', 'imagecol'),
  ...
) {
  cols <- cols %||% colnames(x = slot(object = object, name = 'coordinates'))
  if (!is.null(x = scale)) {
    coordinates <- slot(object = object, name = 'coordinates')[, c('imagerow', 'imagecol')]
    scale <- match.arg(arg = scale, choices = c('spot', 'fiducial', 'hires', 'lowres'))
    scale.use <- ScaleFactors(object = object)[[scale]]
    coordinates <- coordinates * scale.use
  } else {
    coordinates <- slot(object = object, name = 'coordinates')[, cols]
  }
  return(coordinates)
}

#' @rdname GetTissueCoordinates
#' @method GetTissueCoordinates SpatialImage
#' @export
#'
GetTissueCoordinates.SpatialImage <- function(object, ...) {
  stop(
    "'GetTissueCoordinates' must be overridden for all sublcasses of 'SpatialImage'",
    call. = FALSE
  )
}

#' @rdname IsGlobal
#' @method IsGlobal SpatialImage
#' @export
#'
IsGlobal.SpatialImage <- function(object) {
  return(TRUE)
}

#' @rdname RenameCells
#' @method RenameCells SliceImage
#' @export
#'
RenameCells.SliceImage <- function(object, new.names = NULL, ...) {
  if (is.null(x = new.names)) {
    return(object)
  } else if (length(x = new.names) != nrow(x = object)) {
    stop("Wrong number of cell/spot names", call. = FALSE)
  }
  names(x = new.names) <- Cells(x = object)
  coordinates <- GetTissueCoordinates(object = object, scale = NULL, cols = NULL)
  rownames(x = coordinates) <- new.names[rownames(x = coordinates)]
  slot(object = object, name = 'coordinates') <- coordinates
  return(object)
}

#' @rdname RenameCells
#' @method RenameCells SpatialImage
#' @export
#'
RenameCells.SpatialImage <- function(object, ...) {
  stop(
    "'RenameCells' must be overwritten for all subclasses of 'SpatialImage'",
    call. = FALSE
  )
}

#' @rdname ScaleFactors
#' @method ScaleFactors SliceImage
#' @export
#'
ScaleFactors.SliceImage <- function(object, ...) {
  return(slot(object = object, name = 'scale.factors'))
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Methods for R-defined generics
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @method [ SliceImage
#' @export
#'
"[.SliceImage" <- function(x, i, ...) {
  return(subset(x = x, cells = i))
}

#' @method [ SpatialImage
#' @export
#'
"[.SpatialImage" <- function(x, ...) {
  stop(
    "'[' must be overwritten for all subclasses of 'SpatialImage'",
    call. = FALSE
  )
}

#' @method dim SliceImage
#' @export
#'
dim.SliceImage <- function(x) {
  return(dim(x = GetImage(object = x)$raster))
}

#' @method dim SpatialImage
#' @export
#'
dim.SpatialImage <- function(x) {
  stop(
    "'dim' must be overwritten for all subclasses of 'SpatialImage'",
    call. = FALSE
  )
}

#' @method subset SliceImage
#' @export
#'
subset.SliceImage <- function(x, cells, ...) {
  coordinates <- GetTissueCoordinates(object = x, scale = NULL, cols = NULL)
  coordinates <- coordinates[cells, ]
  slot(object = x, name = 'coordinates') <- coordinates
  return(x)
}

#' @method subset SpatialImage
#' @export
#'
subset.SpatialImage <- function(x, ...) {
  stop("'subset' must be overwritten for all subclasses of 'SpatialImage'")
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# S4 methods
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Internal
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DefaultImage <- function(object) {
  object <- UpdateSlots(object = object)
  images <- Images(object = object, assay = DefaultAssay(object = object))
  if (length(x = images) < 1) {
    images <- Images(object = object)
  }
  return(images[[1]])
}

#' @importFrom ggplot2 ggplot_build
#'
GGpointToPlotlyBuild <- function(plot, information = NULL) {
  plot.build <- GGpointToBase(plot = plot, do.plot = FALSE)
  data <- ggplot_build(plot = plot)$plot$data
  rownames(x = plot.build) <- rownames(data)
  # Reset the names to 'x' and 'y'
  names(x = plot.build) <- c(
    'x',
    'y',
    names(x = plot.build)[3:length(x = plot.build)]
  )
  # Add the hover information we're looking for
  if (is.null(x = information)) {
    plot.build$feature <- rownames(x = data)
  } else {
    info <- apply(
      X = information,
      MARGIN = 1,
      FUN = function(x, names) {
        return(paste0(names, ': ', x, collapse = '<br>'))
      },
      names = colnames(x = information)
    )
    data.info <- data.frame(
      feature = paste(rownames(x = information), info, sep = '<br>'),
      row.names = rownames(x = information)
    )
    plot.build <- merge(x = plot.build, y = data.info, by = 0)
    rownames(x = plot.build) <- plot.build$Row.names
    plot.build <- plot.build[, which(x = colnames(x = plot.build) != 'Row.names'), drop = FALSE]
  }
  return(plot.build)
}
