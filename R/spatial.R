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
#' @slot key Key for the image
#'
#' @section Provided methods:
#' These methods are defined on the \code{SpatialImage} object and should not be
#' overwritten without careful thought
#' \itemize{
#'   \item \code{\link{DefaultAssay}} and \code{\link{DefaultAssay<-}}
#'   \item \code{\link{Key}} and \code{\link{Key<-}}
#'   \item \code{\link{IsGlobal}}
#'   \item \code{\link{Radius}}; this method \emph{can} be overridden to provide
#'   a spot radius for image objects
#' }
#'
#' @section Required methods:
#' All subclasses of the \code{SpatialImage} class must define the following methods;
#' simply relying on the \code{SpatialImage} method will result in errors. For required
#' parameters and their values, see the \code{Usage} and \code{Arguments} sections
#' \describe{
#'   \item{\code{\link{Cells}}}{Return the cell/spot barcodes associated with each position}
#'   \item{\code{\link{dim}}}{Return the dimensions of the image for plotting in \code{(Y, X)} format}
#'   \item{\code{\link{GetImage}}}{Return image data; by default, must return a grob object}
#'   \item{\code{\link{GetTissueCoordinates}}}{Return tissue coordinates; by default,
#'   must return a two-column data.frame with x-coordinates in the first column and y-coordiantes
#'   in the second}
#'   \item{\code{\link{Radius}}}{Return the spot radius; returns \code{NULL} by
#'   default for use with non-spot image technologies}
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
    'assay' = 'character',
    'key' = 'character'
  )
)

#' The SlideSeq class
#'
#' The SlideSeq class represents spatial information from the Slide-seq platform
#'
#' @inheritSection SpatialImage Slots
#' @slot coordinates ...
#' @slot ...
#'
SlideSeq <- setClass(
  Class = 'SlideSeq',
  contains = 'SpatialImage',
  slots = list(
    'coordinates' = 'data.frame'
  )
)

#' The STARmap class
#'
#' The STARmap class represents spatial information from the STARmap platform
#'
#' @inheritSection SpatialImage Slots
#' @slot ...
#'
STARmap <- setClass(
  Class = 'STARmap',
  contains = 'SpatialImage',
  slots = list(
    'coordinates' = 'data.frame',
    'qhulls' = 'data.frame'
  )
)

#' The VisiumV1 class
#'
#' The VisiumV1 class represents spatial information from the 10X Genomics Visium
#' platform
#'
#' @slot image A three-dimensional array with PNG image data, see
#' \code{\link[png]{readPNG}} for more details
#' @slot scale.factors An object of class \code{\link{scalefactors}}; see
#' \code{\link{scalefactors}} for more information
#' @slot coordinates A data frame with tissue coordinate information
#' @slot spot.radius Single numeric value giving the radius of the spots
#'
#' @name VisiumV1-class
#' @rdname VisiumV1-class
#' @exportClass VisiumV1
#'
VisiumV1 <- setClass(
  Class = 'VisiumV1',
  contains = 'SpatialImage',
  slots = list(
    'image' = 'array',
    'scale.factors' = 'scalefactors',
    'coordinates' = 'data.frame',
    'spot.radius' = 'numeric'
  )
)

setClass(Class = 'SliceImage', contains = 'VisiumV1')

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Get a vector of cell names associated with an image (or set of images)
#'
#' @param object Seurat object
#' @param images Vector of image names
#' @param unlist Return as a single vector of cell names as opposed to a list,
#' named by image name.
#'
#' @return A vector of cell names
#'
#' @examples
#' \dontrun{
#' CellsByImage(object = object, images = "slice1")
#' }
#'
CellsByImage <- function(object, images = NULL, unlist = FALSE) {
  images <- images %||% Images(object = object)
  cells <- sapply(
    X = images,
    FUN = function(x) {
      Cells(x = object[[x]])
    },
    simplify = FALSE,
    USE.NAMES = TRUE
  )
  if (unlist) {
    cells <- unname(obj = unlist(x = cells))
  }
  return(cells)
}

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
#' @param image.dir Path to directory with 10X Genomics visium image data;
#' should include files \code{tissue_lowres_iamge.png},
#' \code{scalefactors_json.json} and \code{tissue_positions_list.csv}
#' @param filter.matrix Filter spot/feature matrix to only include spots that
#' have been determined to be over tissue.
#' @param ... Ignored for now
#'
#' @return A \code{\link{VisiumV1}} object
#'
#' @importFrom png readPNG
#' @importFrom jsonlite fromJSON
#'
#' @seealso \code{\link{VisiumV1}} \code{\link{Load10X_Spatial}}
#'
#' @export
#'
Read10X_Image <- function(image.dir, filter.matrix = TRUE, ...) {
  image <- readPNG(source = file.path(image.dir, 'tissue_lowres_image.png'))
  scale.factors <- fromJSON(txt = file.path(image.dir, 'scalefactors_json.json'))
  tissue.positions <- read.csv(
    file = file.path(image.dir, 'tissue_positions_list.csv'),
    col.names = c('barcodes', 'tissue', 'row', 'col', 'imagerow', 'imagecol'),
    header = FALSE,
    as.is = TRUE,
    row.names = 1
  )
  if (filter.matrix) {
    tissue.positions <- tissue.positions[which(x = tissue.positions$tissue == 1), , drop = FALSE]
  }
  unnormalized.radius <- scale.factors$fiducial_diameter_fullres * scale.factors$tissue_lowres_scalef
  spot.radius <-  unnormalized.radius / max(dim(x = image))
  return(new(
    Class = 'VisiumV1',
    image = image,
    scale.factors = scalefactors(
      spot = scale.factors$tissue_hires_scalef,
      fiducial = scale.factors$fiducial_diameter_fullres,
      hires = scale.factors$tissue_hires_scalef,
      scale.factors$tissue_lowres_scalef
    ),
    coordinates = tissue.positions,
    spot.radius = spot.radius
  ))
}

#' Load Slide-seq spatial data
#'
#' @param coord.file Path to csv file containing bead coordinate positions
#' @param assay Name of assay to associate image to
#'
#' @return A \code{\link{SlideSeq}} object
#'
#' @importFrom utils read.csv
#'
#' @seealso \code{\link{SlideSeq}}
#'
#' @export
#'
ReadSlideSeq <- function(coord.file, assay = 'Spatial') {
  if (!file.exists(paths = coord.file)) {
    stop("Cannot find coord file ", coord.file, call. = FALSE)
  }
  slide.seq <- new(
    Class = 'SlideSeq',
    assay = assay,
    coordinates = read.csv(
      file = coord.file,
      header = TRUE,
      as.is = TRUE,
      row.names = 1
    )
  )
  return(slide.seq)
}


#' Load STARmap data
#'
#' @param data.dir location of data directory that contains the counts matrix,
#' gene name, qhull, and centroid files. 
#' @param counts.file name of file containing the counts matrix (csv)
#' @param gene.file name of file containing the gene names (csv)
#' @param qhull.file name of file containing the hull coordinates (tsv)
#' @param centroidlfile name of file containing the centroid positions (tsv)
#' @param assay Name of assay to associate spatial data to
#' @param image Name of "image" object storing spatial coordinates
#'
#' @return A \code{\link{Seurat}} object
#'
#' @importFrom utils read.csv read.table
#'
#' @seealso \code{\link{STARmap}}
#'
#' @export
#'
LoadSTARmap <- function(
  data.dir, 
  counts.file = "cell_barcode_count.csv", 
  gene.file = "genes.csv",
  qhull.file = "qhulls.tsv",
  centroid.file = "centroids.tsv",
  assay = "Spatial",
  image = "image"
) {
  if (!dir.exists(paths = data.dir)) {
    stop("Cannot find directory ", data.dir, call. = FALSE)
  }
  counts <- read.csv(
    file = file.path(data.dir, counts.file),
    as.is = TRUE,
    header = FALSE
  )
  gene.names <- read.csv(
    file = file.path(data.dir, gene.file),
    as.is = TRUE,
    header = FALSE
  )
  qhulls <- read.table(
    file = file.path(data.dir, qhull.file),
    sep = '\t',
    col.names = c('cell', 'y', 'x'),
    as.is = TRUE
  )
  centroids <- read.table(
    file = file.path(data.dir, centroid.file),
    sep = '\t',
    as.is = TRUE,
    col.names = c('y', 'x')
  )
  colnames(x = counts) <- gene.names[, 1]
  rownames(x = counts) <- paste0('starmap', seq(1:nrow(x = counts)))
  counts <- as.matrix(x = counts)
  rownames(x = centroids) <- rownames(x = counts)
  qhulls$cell <- paste0('starmap', qhulls$cell)
  centroids <- as.matrix(x = centroids)
  starmap <- CreateSeuratObject(counts = t(x = counts), assay = assay)
  starmap[[image]] <- new(
    Class = 'STARmap',
    assay = assay,
    coordinates = as.data.frame(x = centroids),
    qhulls = qhulls
  )
  return(starmap)
}

#' Link two ggplot plots together
#'
#' @inheritParams HoverLocator
#' @param plot1,plot2 \code{ggplot} objects
#' @param plot1.labels,plot2.labels Logical values to inherit X/Y labels from
#' component \code{ggplot} objects
#' @param pt.size Point size in px
#' @param plot1.cols,plot2.cols A named vector of column names to pull. Vector
#' names must be 'x', 'y', 'colour', 'shape', and/or 'size'; vector values must
#' be the names of columns in plot data that correspond to these values. May
#' pass only values that differ from the default
#' (eg. \code{cols = c('size' = 'point.size.factor')})
#' @param plot1.layout,plot2.layout Extra information for \code{\link[plotly]{layout}}
#'
#' @return A \code{\link[htmltools]{browsable}} HTML element with the two plots
#' linked and side-by-side
#'
#' @importFrom crosstalk SharedData bscols
#' @importFrom plotly layout plot_ly highlight
#'
#' @export
#'
#' @seealso \code{\link[plotly]{layout}}
#'
LinkPlots <- function(
  plot1,
  plot2,
  plot1.labels = TRUE,
  plot2.labels = TRUE,
  information = NULL,
  pt.size = 6,
  plot1.cols = eval(formals(fun = GGpointToBase)$cols),
  plot2.cols = eval(formals(fun = GGpointToBase)$cols),
  plot1.layout = list(),
  plot2.layout = list()
) {
  # Get plot builds for each plot
  plot1.build <- GGpointToPlotlyBuild(
    plot = plot1,
    information = information,
    cols = plot1.cols
  )
  plot1.labels <- if (plot1.labels) {
    GetXYAesthetics(plot = plot1)
  } else {
    list()
  }
  plot2.build <- GGpointToPlotlyBuild(
    plot = plot2,
    information = information,
    cols = plot2.cols
  )
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
      size = ~I(pt.size),
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
      size = ~I(pt.size),
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

#' Visualize spatial and clustering (dimensional reduction) data in a linked,
#' interactive framework
#'
#' @inheritParams FeaturePlot
#' @inheritParams SpatialPlot
#' @param feature Feature to visualize
#' @param image Which image to use
#'
#' @return A \code{\link[htmltools]{browsable}} HTML element with the spatial plot
#' on the left and the dimensional-reduction plot on the right
#'
#' @rdname LinkedPlots
#' @name LinkedPlots
#'
#' @aliases LinkedPlot
#'
#' @export
#'
#' @examples
#' \dontrun{
#' LinkedFeaturePlot(seurat.object, feature = 'Hpca')
#' LinkedDimPlot(seurat.object)
#' }
#'
LinkedFeaturePlot <- function(
  object,
  feature,
  dims = 1:2,
  reduction = NULL,
  pt.size = 6,
  image = NULL,
  slot = 'data',
  min.cutoff = NA,
  max.cutoff = NA
) {
  if (length(x = feature) > 1) {
    stop("'LinkedFeaturePlot' currently only supports one feature", call. = FALSE)
  }
  image <- image %||% DefaultImage(object = object)
  expression.data <- FetchData(
    object = object,
    vars = feature,
    cells = Cells(x = object[[image]])
  )
  coords <- GetTissueCoordinates(object = object[[image]])
  spatial.plot <- SingleSpatialPlot(
    data = cbind(coords, expression.data),
    image = object[[image]],
    col.by = feature,
    geom = 'interactive'
  )
  spatial.plot <- spatial.plot + scale_color_gradientn(
    colours = SpatialColors(n = 100)
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
    pt.size = pt.size,
    plot1.layout = list(
      'images' = GetImage(object = object[[image]], mode = 'plotly'),
      'xaxis' = list('visible' = FALSE),
      'yaxis' = list('visible' = FALSE)
    )
  ))
}

#' @rdname LinkedPlots
#'
#' @export
#'
LinkedDimPlot <- function(
  object,
  dims = 1:2,
  reduction = NULL,
  pt.size = 6,
  image = NULL,
  group.by = NULL
) {
  image <- image %||% DefaultImage(object = object)
  group.by <- group.by %||% 'ident'
  group.data <- FetchData(
    object = object,
    vars = group.by,
    cells = Cells(x = object[[image]])
  )
  coords <- GetTissueCoordinates(object = object[[image]])
  spatial.plot <- SingleSpatialPlot(
    data = cbind(coords, group.data),
    image = object[[image]],
    col.by = group.by,
    geom = 'interactive'
  )
  dim.plot <- DimPlot(
    object = object,
    group.by = group.by,
    dims = dims,
    reduction = reduction
  )
  linked.plot <- LinkPlots(
    plot1 = spatial.plot,
    plot2 = dim.plot,
    information = FetchData(object = object, vars = group.by),
    plot1.labels = FALSE,
    pt.size = pt.size,
    plot1.layout = list(
      'images' = GetImage(object = object[[image]], mode = 'plotly'),
      'xaxis' = list('visible' = FALSE),
      'yaxis' = list('visible' = FALSE)
    )
  )
  suppressMessages(expr = suppressWarnings(expr = print(
    x = linked.plot
  )))
  return(linked.plot)
}

#' Load a 10x Genomics Visium Spatial Experiment into a \code{Seurat} object
#'
#' @inheritParams Read10X
#' @inheritParams CreateSeuratObject
#' @param slice Name for the stored image of the tissue slice
#' @param filter.matrix Only keep spots that have been determined to be over
#' tissue
#' @param to.upper Converts all feature names to upper case. Can be useful when
#' analyses require comparisons between human and mouse gene names for example.
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
#' @examples
#' \dontrun{
#' data_dir <- 'path/to/data/directory'
#' list.files(data_dir) # Should show filtered_feature_bc_matrix.h5
#' Load10X_Spatial(data.dir = data_dir)
#' }
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
  image <- Read10X_Image(
    image.dir = file.path(data.dir, 'spatial'),
    filter.matrix = filter.matrix
  )
  image <- image[Cells(x = object)]
  DefaultAssay(object = image) <- assay
  object[[slice]] <- image
  return(object)
}

#' @importFrom grid viewport editGrob grobName
#' @importFrom ggplot2 ggproto Geom ggproto_parent
#
GeomSpatialInteractive <- ggproto(
  "GeomSpatialInteractive",
  Geom,
  setup_data = function(self, data, params) {
    data <- ggproto_parent(parent = Geom, self = self)$setup_data(data, params)
    data
  },
  draw_group = function(data, panel_scales, coord) {
    vp <- viewport(x = data$x, y = data$y)
    g <- editGrob(grob = data$grob[[1]], vp = vp)
    # Replacement for ggname
    g$name <- grobName(grob = g, prefix = 'geom_spatial_interactive')
    return(g)
    # return(ggname(prefix = "geom_spatial", grob = g))
  },
  required_aes = c("grob","x","y")
)

#' @importFrom ggplot2 layer
#
geom_spatial_interactive <-  function(
  mapping = NULL,
  data = NULL,
  stat = "identity",
  position = "identity",
  na.rm = FALSE,
  show.legend = NA,
  inherit.aes = FALSE,
  ...
) {
  layer(
    geom = GeomSpatialInteractive,
    mapping = mapping,
    data = data,
    stat = stat,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(na.rm = na.rm, ...)
  )
}

# For plotting the tissue image
#' @importFrom ggplot2 ggproto Geom aes ggproto_parent alpha draw_key_point
#' @importFrom grid unit gpar editGrob pointsGrob viewport gTree addGrob grobName
#' @export
#'
GeomSpatial <- ggproto(
  "GeomSpatial",
  Geom,
  required_aes = c("x", "y"),
  extra_params = c("na.rm", "image", "crop"),
  default_aes = aes(
    shape = 21,
    colour = "black",
    point.size.factor = 1.0,
    fill = NA,
    alpha = NA,
    stroke = 0.25
  ),
  setup_data = function(self, data, params) {
    data <- ggproto_parent(Geom, self)$setup_data(data, params)
    # We need to flip the image as the Y coordinates are reversed
    data$y = max(data$y) - data$y + min(data$y)
    data
  },
  draw_key = draw_key_point,
  draw_panel = function(data, panel_scales, coord, image, crop) {
    # This should be in native units, where
    # Locations and sizes are relative to the x- and yscales for the current viewport.
    if (!crop) {
      y.transform <- c(0, nrow(x = image)) - panel_scales$y.range
      data$y <- data$y + sum(y.transform)
      panel_scales$y.range <- c(0, nrow(x = image))
      panel_scales$x.range <- c(0, ncol(x = image))
    }
    z <- coord$transform(
      data.frame(x = c(0, ncol(x = image)), y = c(0, nrow(x = image))),
      panel_scales
    )
    # Flip Y axis for image
    z$y <- -rev(z$y) + 1
    wdth <- z$x[2] - z$x[1]
    hgth <- z$y[2] - z$y[1]
    vp <- viewport(
      x = unit(x = z$x[1], units = "npc"),
      y = unit(x = z$y[1], units = "npc"),
      width = unit(x = wdth, units = "npc"),
      height = unit(x = hgth, units = "npc"),
      just = c("left", "bottom")
    )
    img.grob <- GetImage(object = image)
    img <- editGrob(grob = img.grob, vp = vp)
    # spot.size <- slot(object = image, name = "spot.radius")
    spot.size <- Radius(object = image)
    coords <- coord$transform(data, panel_scales)
    pts <- pointsGrob(
      x = coords$x,
      y = coords$y,
      pch = data$shape,
      size = unit(spot.size, "npc") * data$point.size.factor,
      gp = gpar(
        col = alpha(coords$colour, coords$alpha),
        fill = alpha(coords$fill, coords$alpha),
        lwd = coords$stroke)
    )
    vp <- viewport()
    gt <- gTree(vp = vp)
    gt <- addGrob(gTree = gt, child = img)
    gt <- addGrob(gTree = gt, child = pts)
    # Replacement for ggname
    gt$name <- grobName(grob = gt, prefix = 'geom_spatial')
    return(gt)
    # ggplot2:::ggname("geom_spatial", gt)
  }
)

# influenced by: https://stackoverflow.com/questions/49475201/adding-tables-to-ggplot2-with-facet-wrap-in-r
# https://ggplot2.tidyverse.org/articles/extending-ggplot2.html
#' @importFrom ggplot2 layer
#'
#' @export
#'
geom_spatial <-  function(
  mapping = NULL,
  data = NULL,
  image = image,
  crop = crop,
  stat = "identity",
  position = "identity",
  na.rm = FALSE,
  show.legend = NA,
  inherit.aes = TRUE,
  ...
) {
  layer(
    geom = GeomSpatial,
    mapping = mapping,
    data = data,
    stat = stat,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(na.rm = na.rm, image = image, crop = crop, ...)
  )
}

#' Run the mark variogram computation on a given Seurat object
#'
#' Wraps the functionality of markvario from the spatstat package.
#'
#' @param spatial.location A 2 column matrix giving the spatial locations of
#' each of the data points also in data
#' @param data Matrix containing the data used as "marks" (e.g. gene expression)
#' @param ... Arguments passed to markvario
#'
#' @importFrom spatstat markvario ppp
#'
#' @export
#'
RunMarkVario <- function(
  spatial.location,
  data,
  ...
) {
  pp <- ppp(
    x = spatial.location[, 1],
    y = spatial.location[, 2],
    xrange = range(spatial.location[, 1]),
    yrange = range(spatial.location[, 2])
  )
  if (nbrOfWorkers() > 1) {
    chunks <- nbrOfWorkers()
    features <- rownames(x = data)
    features <- split(
      x = features,
      f = ceiling(x = seq_along(along.with = features) / (length(x = features) / chunks))
    )
    mv <- future_lapply(X = features, FUN = function(x) {
      pp[["marks"]] <- as.data.frame(x = t(x = data[x, ]))
      markvario(X = pp, normalise = TRUE, ...)
    })
    mv <- unlist(x = mv, recursive = FALSE)
    names(x = mv) <- rownames(x = data)
  } else {
    pp[["marks"]] <- as.data.frame(x = t(x = data))
    mv <- markvario(X = pp, normalise = TRUE, ...)
  }
  return(mv)
}

#' @importFrom RColorBrewer brewer.pal
#' @importFrom grDevices colorRampPalette
#'
#'
SpatialColors <- colorRampPalette(colors = rev(x = brewer.pal(n = 11, name = "Spectral")))

# Base plotting function for all Spatial plots
#
# @param data Data.frame with info to be plotted
# @param image SpatialImage object to be plotted
# @param crop Crop the plot in to focus on points plotted. Set to FALSE to show
# entire background image.
# @param pt.size.factor Sets the size of the points relative to spot.radius
# @param stroke Control the width of the border around the spots
# @param col.by Mapping variable for the point color
# @param alpha.by Mapping variable for the point alpha value
# @param cells.highlight A list of character or numeric vectors of cells to
# highlight. If only one group of cells desired, can simply pass a vector
# instead of a list. If set, colors selected cells to the color(s) in
# cols.highlight
# @param cols.highlight A vector of colors to highlight the cells as; ordered
# the same as the groups in cells.highlight; last color corresponds to
# unselected cells.
# @param geom Switch between normal spatial geom and geom to enable hover
# functionality
# @param na.value Color for spots with NA values

#' @importFrom tibble tibble
#' @importFrom ggplot2 ggplot aes_string coord_fixed geom_point xlim ylim
#' coord_cartesian labs theme_void theme
#'
SingleSpatialPlot <- function(
  data,
  image,
  crop = TRUE,
  pt.size.factor = NULL,
  stroke = 0.25,
  col.by = NULL,
  alpha.by = NULL,
  cells.highlight = NULL,
  cols.highlight = c('#DE2D26', 'grey50'),
  geom = c('spatial', 'interactive', 'poly'),
  na.value = 'grey50'
) {
  geom <- match.arg(arg = geom)
  if (!is.null(x = col.by) && !col.by %in% colnames(x = data)) {
    warning("Cannot find '", col.by, "' in data, not coloring", call. = FALSE, immediate. = TRUE)
    col.by <- NULL
  }
  col.by <- col.by %iff% paste0("`", col.by, "`")
  alpha.by <- alpha.by %iff% paste0("`", alpha.by, "`")
  if (!is.null(x = cells.highlight)) {
    highlight.info <- SetHighlight(
      cells.highlight = cells.highlight,
      cells.all = rownames(x = data),
      sizes.highlight = pt.size.factor,
      cols.highlight = cols.highlight[1],
      col.base = cols.highlight[2]
    )
    order <- highlight.info$plot.order
    data$highlight <- highlight.info$highlight
    col.by <- 'highlight'
  }
  plot <- ggplot(data = data, aes_string(
    x = colnames(x = data)[2],
    y = colnames(x = data)[1],
    fill = col.by,
    alpha = alpha.by
  ))
  plot <- switch(
    EXPR = geom,
    'spatial' = {
      plot + geom_spatial(
        point.size.factor = pt.size.factor,
        data = data,
        image = image,
        crop = crop,
        stroke = stroke
      ) + coord_fixed()
    },
    'interactive' = {
      plot + geom_spatial_interactive(
        data = tibble(grob = list(GetImage(object = image, mode = 'grob'))),
        mapping = aes_string(grob = 'grob'),
        x = 0.5,
        y = 0.5
      ) +
        geom_point(mapping = aes_string(color = col.by)) +
        xlim(0, ncol(x = image)) +
        ylim(nrow(x = image), 0) +
        coord_cartesian(expand = FALSE)
    },
    'poly' = {
      data$cell <- rownames(x = data)
      data[, c('x', 'y')] <- NULL
      data <- merge(
        x = data, 
        y = GetTissueCoordinates(object = image, qhulls = TRUE), 
        by = "cell"
      )
      plot + geom_polygon(
        data = data,
        mapping = aes_string(fill = col.by, group = 'cell')
      ) + coord_fixed() + theme_cowplot()
        
    },
    stop("Unknown geom, choose from 'spatial' or 'interactive'", call. = FALSE)
  )
  if (!is.null(x = cells.highlight)) {
    plot <- plot + scale_fill_manual(values = cols.highlight)
  }
  plot <- plot + theme_void()
  return(plot)
}

#' Visualize spatial clustering and expression data.
#'
#' SpatialPlot plots a feature or discrete grouping (e.g. cluster assignments) as
#' spots over the image that was collected. We also provide SpatialFeaturePlot
#' and SpatialDimPlot as wrapper functions around SpatialPlot for a consistent
#' naming framework.
#'
#' @inheritParams HoverLocator
#' @param object A Seurat object
#' @param group.by Name of meta.data column to group the data by
#' @param features Name of the feature to visualize. Provide either group.by OR
#' features, not both.
#' @param images Name of the images to use in the plot(s)
#' @param crop Crop the plot in to focus on points plotted. Set to FALSE to show
#' entire background image.
#' @param slot If plotting a feature, which data slot to pull from (counts,
#' data, or scale.data)
#' @param min.cutoff,max.cutoff Vector of minimum and maximum cutoff
#' values for each feature, may specify quantile in the form of 'q##' where '##'
#' is the quantile (eg, 'q1', 'q10')
#' @param cells.highlight A list of character or numeric vectors of cells to
#' highlight. If only one group of cells desired, can simply pass a vector
#' instead of a list. If set, colors selected cells to the color(s) in
#' cols.highlight
#' @param cols.highlight A vector of colors to highlight the cells as; ordered
#' the same as the groups in cells.highlight; last color corresponds to
#' unselected cells.
#' @param facet.highlight When highlighting certain groups of cells, split each
#' group into its own plot
#' @param label Whether to label the clusters
#' @param label.size Sets the size of the labels
#' @param label.color Sets the color of the label text
#' @param label.box Whether to put a box around the label text (geom_text vs
#' geom_label)
#' @param repel Repels the labels to prevent overlap
#' @param ncol Number of columns if plotting multiple plots
#' @param combine Combine plots into a single gg object; note that if TRUE;
#' themeing will not work when plotting multiple features/groupings
#' @param pt.size.factor Scale the size of the spots.
#' @param alpha Controls opacity of spots. Provide as a vector specifying the
#' min and max
#' @param stroke Control the width of the border around the spots
#' @param do.identify Select points on the plot to highlight them
#' @param identify.ident An optional identity class to set the selected cells to;
#' will return \code{object} with the identity class updated instead of a vector
#' of selected cells
#' @param do.hover Return a plotly view of the data to get info when hovering
#' over points
#'
#' @return If \code{do.identify}, either a vector of cells selected or the object
#' with selected cells set to the value of \code{identify.ident} (if set). Else,
#' if \code{do.hover}, a plotly object with interactive graphics. Else, a ggplot
#' object
#'
#' @importFrom ggplot2 scale_fill_gradientn ggtitle theme element_text scale_alpha
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # For functionality analagous to FeaturePlot
#' SpatialPlot(seurat.object, features = "MS4A1")
#' SpatialFeaturePlot(seurat.object, features = "MS4A1")
#'
#' # For functionality analagous to DimPlot
#' SpatialPlot(seurat.object, group.by = "clusters")
#' SpatialDimPlot(seurat.object, group.by = "clusters")
#' }
#'
SpatialPlot <- function(
  object,
  group.by = NULL,
  features = NULL,
  images = NULL,
  crop = TRUE,
  slot = 'data',
  min.cutoff = NA,
  max.cutoff = NA,
  cells.highlight = NULL,
  cols.highlight = c('#DE2D26', 'grey50'),
  facet.highlight = FALSE,
  label = FALSE,
  label.size = 5,
  label.color = 'white',
  label.box = TRUE,
  repel = FALSE,
  ncol = NULL,
  combine = TRUE,
  pt.size.factor = 1.6,
  alpha = c(1, 1),
  stroke = 0.25,
  do.identify = FALSE,
  identify.ident = NULL,
  do.hover = FALSE,
  information = NULL
) {
  if (!is.null(x = group.by) & !is.null(x = features)) {
    stop("Please specific either group.by or features, not both.")
  }
  if (is.null(x = features)) {
    group.by <- group.by %||% 'ident'
    object[['ident']] <- Idents(object = object)
    data <- object[[group.by]]
    for (group in group.by) {
      if (!is.factor(x = data[, group])) {
        data[, group] <- factor(x = data[, group])
      }
    }
  } else {
    data <- FetchData(
      object = object,
      vars = features,
      slot = slot
    )
    features <- colnames(x = data)
    # Determine cutoffs
    min.cutoff <- mapply(
      FUN = function(cutoff, feature) {
        return(ifelse(
          test = is.na(x = cutoff),
          yes = min(data[, feature]),
          no = cutoff
        ))
      },
      cutoff = min.cutoff,
      feature = features
    )
    max.cutoff <- mapply(
      FUN = function(cutoff, feature) {
        return(ifelse(
          test = is.na(x = cutoff),
          yes = max(data[, feature]),
          no = cutoff
        ))
      },
      cutoff = max.cutoff,
      feature = features
    )
    check.lengths <- unique(x = vapply(
      X = list(features, min.cutoff, max.cutoff),
      FUN = length,
      FUN.VALUE = numeric(length = 1)
    ))
    if (length(x = check.lengths) != 1) {
      stop("There must be the same number of minimum and maximum cuttoffs as there are features")
    }
    # Apply cutoffs
    data <- sapply(
      X = 1:ncol(x = data),
      FUN = function(index) {
        data.feature <- as.vector(x = data[, index])
        min.use <- SetQuantile(cutoff = min.cutoff[index], data.feature)
        max.use <- SetQuantile(cutoff = max.cutoff[index], data.feature)
        data.feature[data.feature < min.use] <- min.use
        data.feature[data.feature > max.use] <- max.use
        return(data.feature)
      }
    )
    colnames(x = data) <- features
    rownames(x = data) <- Cells(x = object)
  }
  images <- images %||% Images(object = object, assay = DefaultAssay(object = object))
  if (length(x = images) == 0) {
    images <- Images(object = object)
  }
  if (length(x = images) < 1) {
    stop("Could not find any spatial image information")
  }
  features <- colnames(x = data)
  colnames(x = data) <- features
  rownames(x = data) <- colnames(x = object)
  facet.highlight <- facet.highlight && (!is.null(x = cells.highlight) && is.list(x = cells.highlight))
  if (do.hover) {
    if (length(x = images) > 1) {
      images <- images[1]
      warning(
        "'do.hover' requires only one image, using image ",
        images,
        call. = FALSE,
        immediate. = TRUE
      )
    }
    if (length(x = features) > 1) {
      features <- features[1]
      type <- ifelse(test = is.null(x = group.by), yes = 'feature', no = 'grouping')
      warning(
        "'do.hover' requires only one ",
        type,
        ", using ",
        features,
        call. = FALSE,
        immediate. = TRUE
      )
    }
    if (facet.highlight) {
      warning(
        "'do.hover' requires no faceting highlighted cells",
        call. = FALSE,
        immediate. = TRUE
      )
      facet.highlight <- FALSE
    }
  }
  if (facet.highlight) {
    if (length(x = images) > 1) {
      images <- images[1]
      warning(
        "Faceting the highlight only works with a single image, using image ",
        images,
        call. = FALSE,
        immediate. = TRUE
      )
    }
    ncols <- length(x = cells.highlight)
  } else {
    ncols <- length(x = images)
  }
  plots <- vector(
    mode = "list",
    length = length(x = features) * ncols
  )
  for (i in 1:ncols) {
    plot.idx <- i
    image.idx <- ifelse(test = facet.highlight, yes = 1, no = i)
    image.use <- object[[images[[image.idx]]]]
    coordinates <- GetTissueCoordinates(object = image.use)
    highlight.use <- if (facet.highlight) {
      cells.highlight[i]
    } else {
      cells.highlight
    }
    for (j in 1:length(x = features)) {
      plot <- SingleSpatialPlot(
        data = cbind(
          coordinates,
          data[rownames(x = coordinates), features[j], drop = FALSE]
        ),
        image = image.use,
        col.by = features[j],
        alpha.by = if (is.null(x = group.by)) {
          features[j]
        } else {
          NULL
        },
        geom = if (do.hover || do.identify) {
          'interactive'
        } else if (inherits(x = image.use, what = "STARmap")) {
          'poly'
        } else {
          'spatial'
        },
        cells.highlight = highlight.use,
        cols.highlight = cols.highlight,
        pt.size.factor = pt.size.factor,
        stroke = stroke,
        crop = crop
      )
      if (is.null(x = group.by)) {
        plot <- plot +
          scale_fill_gradientn(
            name = features[j],
            colours = SpatialColors(n = 100)
          ) +
          theme(legend.position = 'top') +
          scale_alpha(range = alpha) +
          guides(alpha = FALSE)
      } else if (label) {
        plot <- LabelClusters(
          plot = plot,
          id = ifelse(
            test = is.null(x = cells.highlight),
            yes = features[j],
            no = 'highlight'
          ),
          geom = if (inherits(x = image.use, what = "STARmap")) {
            'GeomPolygon'
          } else {
            'GeomSpatial'
          },
          repel = repel,
          size = label.size,
          color = label.color,
          box = label.box,
          position = "nearest"
        )
      }
      if (j == 1 && length(x = images) > 1 && !facet.highlight) {
        plot <- plot +
          ggtitle(label = images[[image.idx]]) +
          theme(plot.title = element_text(hjust = 0.5))
      }
      if (facet.highlight) {
        plot <- plot +
          ggtitle(label = names(x = cells.highlight)[i]) +
          theme(plot.title = element_text(hjust = 0.5)) +
          NoLegend()
      }
      plots[[plot.idx]] <- plot
      plot.idx <- plot.idx + ncols
    }
  }
  if (do.identify) {
    return(CellSelector(
      plot = plot,
      object = identify.ident %iff% object,
      ident = identify.ident
    ))
  } else if (do.hover) {
    return(HoverLocator(
      plot = plots[[1]],
      information = information %||% data[, features, drop = FALSE],
      axes = FALSE,
      # cols = c('size' = 'point.size.factor', 'colour' = 'fill'),
      images = GetImage(object = object, mode = 'plotly', image = images)
    ))
  }
  if (length(x = images) > 1 && combine) {
    plots <- CombinePlots(plots = plots, ncol = length(x = images))
  } else if (length(x = images == 1) && combine) {
    plots <- CombinePlots(plots = plots, ncol = ncol)
  }
  return(plots)
}

# TODO: move to convienence.R
#' @export
#' @rdname SpatialPlot
#'
SpatialDimPlot <- function(
  object,
  group.by = NULL,
  images = NULL,
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
  stroke = 0.25,
  label.box = TRUE,
  do.identify = FALSE,
  identify.ident = NULL,
  do.hover = FALSE,
  information = NULL
) {
  return(SpatialPlot(
    object = object,
    group.by = group.by,
    images = images,
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
    stroke = stroke,
    label.box = label.box,
    do.identify = do.identify,
    identify.ident = identify.ident,
    do.hover = do.hover,
    information = information
  ))
}

# TODO: move to convienence.R
#' @export
#' @rdname SpatialPlot
#'
SpatialFeaturePlot <- function(
  object,
  features,
  images = NULL,
  crop = TRUE,
  slot = 'data',
  min.cutoff = NA,
  max.cutoff = NA,
  ncol = NULL,
  combine = TRUE,
  pt.size.factor = 1.6,
  alpha = c(1, 1),
  stroke = 0.25,
  do.hover = FALSE,
  do.identify = FALSE,
  identify.ident = NULL,
  information = NULL
) {
  return(SpatialPlot(
    object = object,
    features = features,
    images = images,
    crop = crop,
    slot = slot,
    min.cutoff = min.cutoff,
    max.cutoff = max.cutoff,
    ncol = ncol,
    combine = combine,
    pt.size.factor = pt.size.factor,
    alpha = alpha,
    stroke = stroke,
    do.identify = do.identify,
    identify.ident = identify.ident,
    do.hover = do.hover,
    information = information
  ))
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Methods for Seurat-defined generics
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @method Cells SlideSeq
#' @export
#'
Cells.SlideSeq <- function(x) {
  return(rownames(x = GetTissueCoordinates(object = x)))
}

#' @param x,object An object inheriting from \code{SpatialImage}
#'
#' @rdname SpatialImage-class
#' @name SpatialImage-class
#'
#' @method Cells SpatialImage
#' @export
#'
Cells.SpatialImage <- function(x) {
  stop(
    "'Cells' must be overridden for all subclasses of 'SpatialImage'",
    call. = FALSE
  )
}

#' @method Cells STARmap
#' @export
#'
Cells.STARmap <- function(x) {
  return(rownames(x = GetTissueCoordinates(object = x)))
}

#' @rdname Cells
#' @method Cells VisiumV1
#' @export
#'
Cells.VisiumV1 <- function(x) {
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

#' @method FindSpatiallyVariableFeatures Assay
#' @rdname FindSpatiallyVariableFeatures
#' @export
#'
#'
FindSpatiallyVariableFeatures.default <- function(
  object,
  spatial.location,
  selection.method = 'markvariogram',
  r.metric = 5,
  ...
) {
  # error check dimensions
  if (ncol(x = object) != nrow(x = spatial.location)) {
    stop("Please provide the same number of observations as spatial locations.")
  }
  if (selection.method == "markvariogram") {
    svf.info <- RunMarkVario(
      spatial.location = spatial.location,
      data = object
    )
  }
  return(svf.info)
}

#' @param spatial.location
#'
#' @method FindSpatiallyVariableFeatures Assay
#' @rdname FindSpatiallyVariableFeatures
#' @export
#'
FindSpatiallyVariableFeatures.Assay <- function(
  object,
  slot = "scale.data",
  features = NULL,
  spatial.location,
  selection.method = 'markvariogram',
  r.metric = 5,
  nfeatures = nfeatures,
  ...
) {
  features <- features %||% rownames(x = object)
  if (selection.method == "markvariogram" && "markvariogram" %in% names(x = Misc(object = object))) {
    features.computed <- names(x = Misc(object = object, slot = "markvariogram"))
    features <- features[! features %in% features.computed]
  }
  data <- GetAssayData(object = object, slot = slot)
  data <- as.matrix(x = data[features, ])
  data <- data[RowVar(x = data) > 0, ]
  if (nrow(x = data) != 0) {
    svf.info <- FindSpatiallyVariableFeatures(
      object = data,
      spatial.location = spatial.location,
      selection.method = selection.method,
      r.metric = r.metric,
      ...
    )
  } else {
    svf.info <- c()
  }
  if (selection.method == "markvariogram") {
    if ("markvariogram" %in% names(x = Misc(object = object))) {
      svf.info <- c(svf.info, Misc(object = object, slot = "markvariogram"))
    }
    suppressWarnings(expr = Misc(object = object, slot = "markvariogram") <- svf.info)
    svf.info <- ComputeRMetric(mv = svf.info, r.metric)
    svf.info <- svf.info[order(svf.info[, 1]), , drop = FALSE]
    svf.info$markvariogram.spatially.variable <- FALSE
    svf.info$markvariogram.spatially.variable[1:(min(nrow(x = svf.info), nfeatures))] <- TRUE
    svf.info$markvariogram.spatially.variable.rank <- 1:nrow(x = svf.info)
    object[[names(x = svf.info)]] <- svf.info
  }
  return(object)
}

#' @param object A Seurat object
#' @param assay Assay to pull the features (marks) from
#' @param slot Slot in the Assay to pull data from
#' @param features If provided, only compute on given features. Otherwise,
#' compute for all features.
#' @param image Name of image to pull the coordinates from
#' @param selection.method Method for selecting spatially variable features.
#' Only 'markvariogram' method is currently implemented. See
#' \code{\link{RunMarkVario}} for method.
#' @param r.metric r value at which to report the "trans" value of the mark
#' variogram
#'
#' @method FindSpatiallyVariableFeatures Seurat
#' @rdname FindSpatiallyVariableFeatures
#' @export
#'
FindSpatiallyVariableFeatures.Seurat <- function(
  object,
  assay = NULL,
  slot = "scale.data",
  features = NULL,
  image = NULL,
  selection.method = c("markvariogram"),
  r.metric = 5,
  nfeatures = 2000,
  verbose = TRUE,
  ...
) {
  assay <- assay %||% DefaultAssay(object = object)
  features <- features %||% rownames(x = object[[assay]])
  image <- image %||% DefaultImage(object = object)
  tc <- GetTissueCoordinates(object = object[[image]])
  # check if markvariogram has been run on necessary features
  # only run for new ones
  object[[assay]] <- FindSpatiallyVariableFeatures(
    object = object[[assay]],
    slot = slot,
    features = features,
    spatial.location = tc,
    selection.method = selection.method,
    r.metric = r.metric,
    nfeatures = nfeatures,
    ...
  )
  object <- LogSeuratCommand(object = object)
}

#' @param image Name of \code{SpatialImage} object to pull image data for; if
#' \code{NULL}, will attempt to select an image automatically
#'
#' @rdname GetImage
#' @method GetImage Seurat
#' @export
#'
GetImage.Seurat <- function(
  object,
  mode = c('grob', 'raster', 'plotly', 'raw'),
  image = NULL,
  ...
) {
  mode <- match.arg(arg = mode)
  image <- image %||% DefaultImage(object = object)
  if (is.null(x = image)) {
    stop("No images present in this Seurat object", call. = FALSE)
  }
  return(GetImage(object = object[[image]], mode = mode, ...))
}

#' @method GetImage SlideSeq
#' @export
#'
GetImage.SlideSeq <- function(
  object,
  mode = c('grob', 'raster', 'plotly', 'raw'),
  ...
) {
  mode <- match.arg(arg = mode)
  return(NullImage(mode = mode))
}

#' @inheritParams GetImage
#'
#' @rdname SpatialImage-class
#' @name SpatialImage-class
#'
#' @method GetImage SpatialImage
#' @export
#'
GetImage.SpatialImage <- function(
  object,
  mode = c('grob', 'raster', 'plotly', 'raw'),
  ...
) {
  mode <- match.arg(arg = mode)
  stop(
    "'GetImage' must be overridden for all sublcasses of 'SpatialImage'",
    call. = FALSE
  )
}

#' @method GetImage STARmap
#' @export
#'
GetImage.STARmap <- function(
  object,
  mode = c('grob', 'raster', 'plotly', 'raw'),
  ...
) {
  mode <- match.arg(arg = mode)
  return(NullImage(mode = mode))
}

#' @importFrom plotly raster2uri
#' @importFrom grDevices as.raster
#' @importFrom grid rasterGrob unit
#'
#' @rdname GetImage
#' @method GetImage VisiumV1
#' @export
#'
GetImage.VisiumV1 <- function(
  object,
  mode = c('grob', 'raster', 'plotly', 'raw'),
  ...
) {
  mode <- match.arg(arg = mode)
  image <- slot(object = object, name = 'image')
  image <- switch(
    EXPR = mode,
    'grob' = rasterGrob(
      image = image,
      width = unit(x = 1, units = 'npc'),
      height = unit(x = 1, units = 'npc')
    ),
    'raster' = as.raster(x = image),
    'plotly' = list(
      source = raster2uri(r = GetImage(object = object, mode = 'raster')),
      xref = 'x',
      yref = 'y',
      # x = -7,
      # y = -7,
      sizex = ncol(x = object),
      sizey = nrow(x = object),
      sizing = 'stretch',
      opacity = 1,
      layer = 'below'
    ),
    'raw' = image,
    stop("Unknown image mode: ", mode, call. = FALSE)
  )
  return(image)
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

#' @method GetTissueCoordinates SlideSeq
#' @export
#'
GetTissueCoordinates.SlideSeq <- function(object, ...) {
  coords <- slot(object = object, name = 'coordinates')
  colnames(x = coords) <- c('x', 'y')
  # coords$y <- -rev(x = coords$y) + 1
  # coords$y <- FlipCoords(x = coords$y)
  coords$cells <- rownames(x = coords)
  return(coords)
}

#' @inheritParams GetTissueCoordinates
#'
#' @rdname SpatialImage-class
#' @name SpatialImage-class
#'
#' @method GetTissueCoordinates SpatialImage
#' @export
#'
GetTissueCoordinates.SpatialImage <- function(object, ...) {
  stop(
    "'GetTissueCoordinates' must be overridden for all sublcasses of 'SpatialImage'",
    call. = FALSE
  )
}

#' @param qhulls return qhulls instead of centroids
#' @method GetTissueCoordinates STARmap
#' @export
#'
GetTissueCoordinates.STARmap <- function(object, qhulls = FALSE, ...) {
  if (qhulls) {
    return(slot(object = object, name = 'qhulls'))
  } 
  return(slot(object = object, name = 'coordinates'))
}

#' @param scale A factor to scale the coordinates by; choose from: 'tissue',
#' 'fiducial', 'hires', 'lowres', or \code{NULL} for no scaling
#' @param cols Columns of tissue coordinates data.frame to pull
#'
#' @rdname GetTissueCoordinates
#' @method GetTissueCoordinates VisiumV1
#' @export
#'
GetTissueCoordinates.VisiumV1 <- function(
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

#' @rdname IsGlobal
#' @method IsGlobal SpatialImage
#' @export
#'
IsGlobal.SpatialImage <- function(object) {
  return(TRUE)
}

#' @rdname Key
#' @method Key SpatialImage
#' @export
#'
Key.SpatialImage <- function(object, ...) {
  CheckDots(...)
  object <- UpdateSlots(object = object)
  return(slot(object = object, name = 'key'))
}

#' @rdname Key
#' @method Key<- SpatialImage
#' @export
#'
"Key<-.SpatialImage" <- function(object, ..., value) {
  CheckDots(...)
  object <- UpdateSlots(object = object)
  value <- UpdateKey(key = value)
  slot(object = object, name = 'key') <- value
  return(object)
}

#' @rdname Radius
#' @method Radius SlideSeq
#' @export
#'
Radius.SlideSeq <- function(object) {
  return(0.005)
}

#' @rdname SpatialImage-class
#' @name SpatialImage-class
#'
#' @method Radius SpatialImage
#' @export
#'
Radius.SpatialImage <- function(object) {
  return(NULL)
}

#' @rdname Radius
#' @method Radius VisiumV1
#' @export
#'
Radius.VisiumV1 <- function(object) {
  return(slot(object = object, name = 'spot.radius'))
}

#' @method RenameCells SlideSeq
#' @export
#'
RenameCells.SlideSeq <- function(object, new.names = NULL, ...) {
  return(RenameCells.VisiumV1(object = object, new.names = new.names))
}

#' @inheritParams  RenameCells
#'
#' @rdname SpatialImage-class
#' @name SpatialImage-class
#'
#' @method RenameCells SpatialImage
#' @export
#'
RenameCells.SpatialImage <- function(object, new.names = NULL, ...) {
  stop(
    "'RenameCells' must be overwritten for all subclasses of 'SpatialImage'",
    call. = FALSE
  )
}

#' @method RenameCells STARmap
#' @export
#'
RenameCells.STARmap <- function(object, new.names = NULL, ...) {
  names(x = new.names) <- Cells(x = object)
  object <- RenameCells.VisiumV1(object = object, new.names = new.names)
  qhulls <- GetTissueCoordinates(object = object, qhull = TRUE)
  qhulls$cell <- new.names[qhulls$cell]
  slot(object = object, name = "qhulls") <- qhulls
  return(object)
}

#' @rdname RenameCells
#' @method RenameCells VisiumV1
#' @export
#'
RenameCells.VisiumV1 <- function(object, new.names = NULL, ...) {
  if (is.null(x = new.names)) {
    return(object)
  } else if (length(x = new.names) != length(x = Cells(x = object))) {
    stop("Wrong number of cell/spot names", call. = FALSE)
  }
  names(x = new.names) <- Cells(x = object)
  coordinates <- GetTissueCoordinates(object = object, scale = NULL, cols = NULL)
  rownames(x = coordinates) <- new.names[rownames(x = coordinates)]
  slot(object = object, name = 'coordinates') <- coordinates
  return(object)
}

#' @rdname ScaleFactors
#' @method ScaleFactors VisiumV1
#' @export
#'
ScaleFactors.VisiumV1 <- function(object, ...) {
  return(slot(object = object, name = 'scale.factors'))
}

#' @rdname SpatiallyVariableFeatures
#' @export
#' @method SpatiallyVariableFeatures Assay
#'
SpatiallyVariableFeatures.Assay <- function(
  object,
  selection.method = "markvariogram",
  decreasing = TRUE,
  ...
) {
  CheckDots(...)
  vf <- SVFInfo(object = object, selection.method = selection.method, status = TRUE)
  vf <- vf[rownames(x = vf)[which(x = vf[, "variable"][, 1])], ]
  if (!is.null(x = decreasing)) {
    vf <- vf[order(x = vf[, "rank"], decreasing = !decreasing), ]
  }
  return(rownames(x = vf)[which(x = vf[, "variable"][, 1])])
}

#' @param assay Name of assay to pull spatially variable features for
#'
#' @rdname SpatiallyVariableFeatures
#' @export
#' @method SpatiallyVariableFeatures Seurat
#'
SpatiallyVariableFeatures.Seurat <- function(
  object,
  assay = NULL,
  selection.method = "markvariogram",
  decreasing = TRUE,
  ...
) {
  CheckDots(...)
  assay <- assay %||% DefaultAssay(object = object)
  return(SpatiallyVariableFeatures(object = object[[assay]], selection.method = selection.method, decreasing = decreasing))
}

#' @param selection.method Which method to pull. Options: markvariogram
#' @param status Add variable status to the resulting data.frame
#'
#' @rdname SVFInfo
#' @export
#' @method SVFInfo Assay
#'
SVFInfo.Assay <- function(object, selection.method = "markvariogram", status = FALSE, ...) {
  CheckDots(...)
  vars <- switch(
    EXPR = selection.method,
    'markvariogram' = grep(pattern = "r.metric", x = colnames(x = object[[]]), value = TRUE),
    stop("Unknown method: '", selection.method, "'", call. = FALSE)
  )
  tryCatch(
    expr = svf.info <- object[[vars]],
    error = function(e) {
      stop(
        "Unable to find highly variable feature information for method '",
        selection.method,
        "'",
        call. = FALSE
      )
    }
  )
  colnames(x = svf.info) <- vars
  if (status) {
    svf.info$variable <- object[[paste0(selection.method, '.spatially.variable')]]
    svf.info$rank <- object[[paste0(selection.method, '.spatially.variable.rank')]]
  }
  return(svf.info)
}

#' @param assay Name of assay to pull highly variable feature information for
#'
#' @importFrom tools file_path_sans_ext
#'
#' @rdname HVFInfo
#' @export
#' @method HVFInfo Seurat
#'
#' @examples
#' # Get the HVF info from a specific Assay in a Seurat object
#' HVFInfo(object = pbmc_small, assay = "RNA")[1:5, ]
#'
SVFInfo.Seurat <- function(
  object,
  selection.method = "markvariogram",
  assay = NULL,
  status = FALSE,
  ...
) {
  CheckDots(...)
  assay <- assay %||% DefaultAssay(object = object)
  return(SVFInfo(
    object = GetAssay(object = object, assay = assay),
    selection.method = selection.method,
    status = status
  ))
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Methods for R-defined generics
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @method [ SlideSeq
#' @export
#'
"[.SlideSeq" <- function(x, i, ...) {
  return(subset(x = x, cells = i, ...))
}

#' @param i,cells A vector of cells to keep
#'
#' @rdname SpatialImage-class
#' @name SpatialImage-class
#'
#' @method [ SpatialImage
#' @export
#'
"[.SpatialImage" <- function(x, i, ...) {
  stop(
    "'[' must be overwritten for all subclasses of 'SpatialImage'",
    call. = FALSE
  )
}

#' @method [ VisiumV1
#' @export
#'
"[.VisiumV1" <- function(x, i, ...) {
  return(subset(x = x, cells = i))
}

#' @method dim SlideSeq
#' @export
#'
dim.SlideSeq <- function(x) {
  # return(dim(x = GetImage(object = x, mode = 'raw')))
  return(c(599, 600))
}

#' @rdname SpatialImage-class
#' @name SpatialImage-class
#'
#' @method dim SpatialImage
#' @export
#'
dim.SpatialImage <- function(x) {
  stop(
    "'dim' must be overwritten for all subclasses of 'SpatialImage'",
    call. = FALSE
  )
}

dim.STARmap <- function(x) {
  coords <- GetTissueCoordinates(object = x)
  return(c(
    max(coords[, 1]) - min(coords[, 1]), 
    max(coords[, 2]) - min(coords[, 2])
  ))
}

#' @method dim VisiumV1
#' @export
#'
dim.VisiumV1 <- function(x) {
  return(dim(x = GetImage(object = x)$raster))
}

#' @method subset SlideSeq
#' @export
#'
subset.SlideSeq <- function(x, cells, ...) {
  .NotYetImplemented()
}

#' @method subset STARmap
#' @export
#'
subset.STARmap <- function(x, cells, ...) {
  x <- subset.VisiumV1(x = x, cells = cells, ...)
  qhulls <- GetTissueCoordinates(object = x, qhulls = TRUE)
  qhulls <- qhulls[qhulls$cell %in% cells, ]
  slot(object = x, name = 'qhulls') <- qhulls
  return(x)
}

#' @rdname SpatialImage-class
#' @name SpatialImage-class
#'
#' @method subset SpatialImage
#' @export
#'
subset.SpatialImage <- function(x, cells, ...) {
  stop("'subset' must be overwritten for all subclasses of 'SpatialImage'")
}

#' @method subset VisiumV1
#' @export
#'
subset.VisiumV1 <- function(x, cells, ...) {
  coordinates <- GetTissueCoordinates(object = x, scale = NULL, cols = NULL)
  cells <- cells[cells %in% rownames(x = coordinates)]
  coordinates <- coordinates[cells, ]
  slot(object = x, name = 'coordinates') <- coordinates
  return(x)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# S4 methods
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Internal
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Computes the metric at a given r (radius) value and stores in meta.features
#
# @param mv Results of running markvario
# @param r.metric r value at which to report the "trans" value of the mark
# variogram
#
# @return Returns a data.frame with r.metric values
#
#
ComputeRMetric <- function(mv, r.metric = 5) {
  r.metric.results <- unlist(x = lapply(
    X = mv,
    FUN = function(x) {
      x$trans[which.min(x = abs(x = x$r - r.metric))]
    }
  ))
  r.metric.results <- as.data.frame(x = r.metric.results)
  colnames(r.metric.results) <- paste0("r.metric.", r.metric)
  return(r.metric.results)
}

# Splits features into groups based on log expression levels
#
# @param object Seurat object
# @param assay Assay for expression data
# @param min.cells Only compute for features in at least this many cells
# @param ngroups Number of groups to split into
#
# @return A Seurat object with the feature group stored as a factor in
# metafeatures
#
#' @importFrom Matrix rowMeans rowSums
#
GetFeatureGroups <- function(object, assay, min.cells = 5, ngroups = 6) {
  cm <- GetAssayData(object = object[[assay]], slot = "counts")
  # subset to keep only genes detected in at least min.cells cells
  cm <- cm[rowSums(cm > 0) >= min.cells, ]
  # use the geometric mean of the features to group them
  # (using the arithmetic mean would usually not change things much)
  # could use sctransform:::row_gmean here but not exported
  feature.gmean <- exp(x = rowMeans(log1p(x = cm))) - 1
  feature.grp.breaks <- seq(
    from = min(log10(x = feature.gmean)) - 10*.Machine$double.eps,
    to = max(log10(x = feature.gmean)),
    length.out = ngroups + 1
  )
  feature.grp <- cut(
    x = log10(x = feature.gmean),
    breaks = feature.grp.breaks,
    ordered_result = TRUE
  )
  feature.grp <- factor(
    x = feature.grp,
    levels = rev(x = levels(x = feature.grp)),
    ordered = TRUE
  )
  names(x = feature.grp) <- names(x = feature.gmean)
  return(feature.grp)
}

#' Compute the correlation of features broken down by groups with another
#' covariate
#'
#' @param object Seurat object
#' @param assay Assay to pull the data from
#' @param slot Slot in the assay to pull feature expression data from (counts,
#' data, or scale.data)
#' @param var Variable with which to correlate the features
#' @param group.assay Compute the gene groups based off the data in this assay.
#' @param min.cells Only compute for genes in at least this many cells
#' @param ngroups Number of groups to split into
#'
#' @return A Seurat object with the correlation stored in metafeatures
#'
#' @export
#'
GroupCorrelation <- function(
  object,
  assay = NULL,
  slot = "scale.data",
  var = NULL,
  group.assay = NULL,
  min.cells = 5,
  ngroups = 6,
  do.plot = TRUE
) {
  assay <- assay %||% DefaultAssay(object = object)
  group.assay <- group.assay %||% assay
  var <- var %||% paste0("nCount_", group.assay)
  gene.grp <- GetFeatureGroups(
    object = object,
    assay = group.assay,
    min.cells = min.cells,
    ngroups = ngroups
  )
  data <- as.matrix(x = GetAssayData(object = object[[assay]], slot = slot))
  data <- data[rowMeans(x = data) != 0, ]
  grp.cors <- apply(
    X = data,
    MARGIN = 1,
    FUN = function(x) {
      cor(x = x, y = object[[var]])
    }
  )
  grp.cors <- grp.cors[names(x = gene.grp)]
  grp.cors <- as.data.frame(x = grp.cors[which(x = !is.na(x = grp.cors))])
  grp.cors$gene_grp <- gene.grp[rownames(x = grp.cors)]
  colnames(x = grp.cors) <- c("cor", "feature_grp")
  object[[assay]][["feature.grp"]] <- grp.cors[, "feature_grp", drop = FALSE]
  object[[assay]][[paste0(var, "_cor")]] <- grp.cors[, "cor", drop = FALSE]
  if (do.plot) {
    print(GroupCorrelationPlot(
      object = object,
      assay = assay,
      feature.group = "feature.grp",
      cor = paste0(var, "_cor")
    ))
  }
  return(object)
}

#' Boxplot of correlation of a variable (e.g. number of UMIs) with expression
#' data
#'
#' @param object Seurat object
#' @param assay Assay where the feature grouping info and correlations are
#' stored
#' @param feature.group Name of the column in meta.features where the feature
#' grouping info is stored
#' @param cor Name of the column in meta.features where correlation info is
#' stored
#'
#' @return Returns a ggplot boxplot of correlations split by group
#'
#' @importFrom ggplot2 geom_boxplot scale_fill_manual geom_hline
#' @importFrom cowplot theme_cowplot
#' @importFrom scales brewer_pal
#'
#' @export
#'
GroupCorrelationPlot <- function(
  object,
  assay = NULL,
  feature.group = "feature.grp",
  cor = "nCount_RNA_cor"
) {
  assay <- assay %||% DefaultAssay(object = object)
  data <- object[[assay]][[c(feature.group, cor)]]
  data <- data[complete.cases(data), ]
  colnames(x = data) <- c('grp', 'cor')
  plot <- ggplot(data = data, aes_string(x = "grp", y = "cor", fill = "grp")) +
    geom_boxplot() +
    theme_cowplot() +
    scale_fill_manual(values = rev(x = brewer_pal(palette = 'YlOrRd')(n = 7))) +
    ylab(paste(
      "Correlation with",
      gsub(x = cor, pattern = "_cor", replacement = "")
    )) +
    geom_hline(yintercept = 0) +
    NoLegend() +
    theme(
      axis.line.x = element_blank(),
      axis.title.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.x = element_blank()
    )
  return(plot)
}

# Get the default image of an object
#
# Attempts to find all images associated with the default assay of the object.
# If none present, finds all images present in the object. Returns the name of
# the first image
#
# @param object A Seurat object
#
# @return The name of the default image
#
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
GGpointToPlotlyBuild <- function(
  plot,
  information = NULL,
  cols = eval(expr = formals(fun = GGpointToBase)$cols),
  ...
) {
  CheckDots(...)
  plot.build <- GGpointToBase(plot = plot, do.plot = FALSE, cols = cols)
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

# Return a null image
#
# @param mode Image representation to return
# see \code{\link{GetImage}} for more details
#
#' @importFrom grid nullGrob
#' @importFrom grDevices as.raster
#
NullImage <- function(mode) {
  image <- switch(
    EXPR = mode,
    'grob' = nullGrob(),
    'raster' = as.raster(x = new(Class = 'matrix')),
    'plotly' = list('visible' = FALSE),
    'raw' = NULL,
    stop("Unknown image mode: ", mode, call. = FALSE)
  )
  return(image)
}
