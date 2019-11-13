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
    'coordinates' = 'data.frame',
    'spot.radius' = 'numeric'
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
  unnormalized.radius <- scale.factors$fiducial_diameter_fullres * scale.factors$tissue_lowres_scalef
  spot.radius <-  unnormalized.radius / max(dim(x = image))
  return(new(
    Class = 'SliceImage',
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

#' Link two ggplot plots together
#'
#' @param plot1,plot2 \code{ggplot} objects
#' @param plot1.labels,plot2.labels Logical values to inherit X/Y labels from
#' component \code{ggplot} objects
#' @param information ...
#' @param pt.size Point size in px
#' @param plot1.cols,plot2.cols A named vector of column names to pull. Vector
#' names must be 'x', 'y', 'colour', 'shape', and/or 'size'; vector values must
#' be the names of columns in plot data that correspond to these values. May
#' pass only values that differ from the default
#' (eg. \code{cols = c('size' = 'point.size.factor')})
#' @param plot1.layout,plot2.layout Extra information for \code{\link[plotly]{layout}}
#' @param ... ...
#'
#' @return A \code{link[htmltools]{browsable}} HTML element with the two plots
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
  plot2.layout = list(),
  ...
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

#' @export
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
    geom = 'custom'
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
    geom = 'custom'
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
    information = FetchData(object = object, vars = group.by),
    plot1.labels = FALSE,
    pt.size = pt.size,
    plot1.layout = list(
      'images' = GetImage(object = object[[image]], mode = 'plotly'),
      'xaxis' = list('visible' = FALSE),
      'yaxis' = list('visible' = FALSE)
    )
  ))
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
  image <- Read10X_Image(data.dir = data.dir, filter.matrix = filter.matrix)
  image <- image[Cells(x = object)]
  DefaultAssay(object = image) <- assay
  object[[slice]] <- image
  return(object)
}

#' @importFrom grid viewport editGrob grobName
#' @importFrom ggplot2 ggproto Geom ggproto_parent
GeomCustom <- ggproto(
  "GeomCustom",
  Geom,
  setup_data = function(self, data, params) {
    data <- ggproto_parent(parent = Geom, self = self)$setup_data(data, params)
    data
  },
  draw_group = function(data, panel_scales, coord) {
    vp <- viewport(x = data$x, y = data$y)
    g <- editGrob(grob = data$grob[[1]], vp = vp)
    # Replacement for ggname
    g$name <- grobName(grob = g, prefix = 'geom_custom')
    return(g)
    # return(ggname(prefix = "geom_spatial", grob = g))
  },
  required_aes = c("grob","x","y")
)

#' @importFrom ggplot2 layer
#'
geom_custom <-  function(
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

# For plotting the tissue image
#' @importFrom ggplot2 ggproto Geom aes ggproto_parent alpha draw_key_point
#' @importFrom grid unit gpar editGrob pointsGrob viewport gTree addGrob grobName
#' @export
#'
GeomSpatial <- ggproto(
  "GeomSpatial",
  Geom,
  required_aes = c("x", "y"),
  extra_params = c("na.rm", "image"),
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
  draw_group = function(data, panel_scales, coord, image) {
    # This should be in native units, where
    # Locations and sizes are relative to the x- and yscales for the current viewport.
    z = coord$transform(
      data.frame(x = c(0, ncol(x = image)), y = c(0, nrow(x = image))),
      panel_scales
    )
    # Flip Y axis for image
    z$y = -rev(z$y) + 1
    wdth = z$x[2] - z$x[1]
    hgth = z$y[2] - z$y[1]
    vp <- grid::viewport(
      x = unit(x = z$x[1], units = "npc"),
      y = unit(x = z$y[1], units = "npc"),
      width = unit(x = wdth, units = "npc"),
      height = unit(x = hgth, units = "npc"),
      just = c("left", "bottom")
    )
    img.grob <- GetImage(object = image)
    img <- editGrob(grob = img.grob, vp = vp)
    spot.size <- slot(object = image, name = "spot.radius")
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
    if (abs(x = data$group[1]) == 1) {
      gt <- addGrob(gTree = gt, child = img)
    }
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
    params = list(na.rm = na.rm, image = image, ...)
  )
}

#' @importFrom RColorBrewer brewer.pal
#' @importFrom grDevices colorRampPalette
#'
SpatialColors <- colorRampPalette(colors = rev(x = brewer.pal(n = 11, name = "Spectral")))

#' @importFrom tibble tibble
#' @importFrom ggplot2 ggplot aes_string coord_fixed geom_point xlim ylim
#' coord_cartesian labs theme_void
#'
SingleSpatialPlot <- function(
  data,
  image,
  pt.size.factor = NULL,
  alpha = 1,
  stroke = 0.25,
  col.by = NULL,
  cells.highlight = NULL,
  cols.highlight = c('#DE2D26', 'grey50'),
  geom = c('spatial', 'custom'),
  na.value = 'grey50'
) {
  geom <- match.arg(arg = geom)
  if (!is.null(x = col.by) && !col.by %in% colnames(x = data)) {
    warning("Cannot find '", col.by, "' in data, not coloring", call. = FALSE, immediate. = TRUE)
    col.by <- NULL
  }
  col.by <- col.by %iff% paste0("`", col.by, "`")
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
    fill = col.by
  ))
  plot <- switch(
    EXPR = geom,
    'spatial' = {
      plot + geom_spatial(
        point.size.factor = pt.size.factor,
        alpha = alpha,
        data = data,
        image = image,
        stroke = stroke
      ) + coord_fixed()
    },
    'custom' = {
      plot + geom_custom(
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
    stop("Unknown geom, choose from 'spatial' or 'custom'", call. = FALSE)
  )
  if (!is.null(x = cells.highlight)) {
    plot <- plot + scale_fill_manual(values = cols.highlight)
  }
  plot <- plot + theme_void()
  return(plot)
}

#' @export
#' @rdname SpatialPlot
#'
SpatialDimPlot <- function(
  object,
  group.by = NULL,
  images = NULL,
  cells.highlight = NULL,
  cols.highlight = c('#DE2D26', 'grey50'),
  label = FALSE,
  label.size = 7,
  label.color = 'white',
  repel = FALSE,
  ncol = NULL,
  combine = TRUE,
  pt.size.factor = 1,
  alpha = 1,
  stroke = 0.25,
  box = TRUE,
  do.hover = FALSE
) {
  return(SpatialPlot(
    object = object,
    group.by = group.by,
    images = images,
    cells.highlight = cells.highlight,
    cols.highlight = cols.highlight,
    label = label,
    label.size = label.size,
    label.color = label.color,
    repel = repel,
    ncol = ncol,
    combine = combine,
    pt.size.factor = pt.size.factor,
    alpha = alpha,
    stroke = stroke,
    box = box,
    do.hover = do.hover
  ))
}

#' @export
#' @rdname SpatialPlot
#'
SpatialFeaturePlot <- function(
  object,
  features,
  images = NULL,
  slot = 'data',
  min.cutoff = NA,
  max.cutoff = NA,
  ncol = NULL,
  combine = TRUE,
  pt.size.factor = 1,
  alpha = 1,
  stroke = 0.25,
  do.hover = FALSE
) {
  return(SpatialPlot(
    object = object,
    features = features,
    images = images,
    slot = slot,
    min.cutoff = min.cutoff,
    max.cutoff = max.cutoff,
    ncol = ncol,
    combine = combine,
    pt.size.factor = pt.size.factor,
    alpha = alpha,
    stroke = stroke,
    do.hover = do.hover
  ))
}

#' Visualize spatial clustering and expression data
#'
#' @return A ggplot object
#'
#' @importFrom ggplot2 scale_fill_gradientn ggtitle theme element_text
#'
#' @export
#'
SpatialPlot <- function(
  object,
  group.by = NULL,
  features = NULL,
  images = NULL,
  slot = 'data',
  min.cutoff = NA,
  max.cutoff = NA,
  cells.highlight = NULL,
  cols.highlight = c('#DE2D26', 'grey50'),
  label = FALSE,
  label.size = 5,
  label.color = 'white',
  repel = FALSE,
  box = TRUE,
  ncol = NULL,
  combine = TRUE,
  pt.size.factor = 1,
  alpha = 1,
  stroke = 0.25,
  do.hover = FALSE
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
  plots <- vector(
    mode = "list",
    length = length(x = images) * length(x = features)
  )
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
  }
  for (i in 1:length(x = images)) {
    plot.idx <- i
    image.use <- object[[images[[i]]]]
    coordinates <- GetTissueCoordinates(object = image.use)
    slice.plots <- vector(mode = "list", length(x = features))
    for (j in 1:length(x = features)) {
      plot <- SingleSpatialPlot(
        data = cbind(
          coordinates,
          data[rownames(x = coordinates), features[j], drop = FALSE]
        ),
        image = image.use,
        col.by = features[j],
        geom = ifelse(test = do.hover, yes = 'custom', no = 'spatial'),
        cells.highlight = cells.highlight,
        cols.highlight = cols.highlight,
        pt.size.factor = pt.size.factor,
        alpha = alpha,
        stroke = stroke
      )
      if (is.null(x = group.by)) {
        plot <- plot + scale_fill_gradientn(name = features[j], colours = SpatialColors(100))
      } else if (label) {
        plot <- LabelClusters(
          plot = plot,
          id = features[j],
          geom = 'GeomSpatial',
          repel = repel,
          size = label.size,
          color = label.color,
          box = box
        )
      }
      if (j == 1 | length(x = images) == 1) {
        plot <- plot + ggtitle(label = images[[i]]) + theme(plot.title = element_text(hjust = 0.5))
      }
      plots[[plot.idx]] <- plot
      plot.idx <- plot.idx + length(x = images)
    }
  }
  if (do.hover) {
    return(HoverLocator(
      plot = plots[[1]],
      information = data[, features, drop = FALSE],
      axes = FALSE,
      # cols = c('size' = 'point.size.factor', 'colour' = 'fill'),
      images = GetImage(object = object, mode = 'plotly', image = images)
    ))
  }
  if (length(x = images) > 1 & combine) {
    plots <- CombinePlots(plots = plots, ncol = length(x = images))
  } else if (length(x = images == 1) & combine) {
    plots <- CombinePlots(plots = plots, ncol = ncol)
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

#' @importFrom plotly raster2uri
#' @importFrom grDevices as.raster
#' @importFrom grid rasterGrob unit
#'
#' @rdname GetImage
#' @method GetImage SliceImage
#' @export
#'
GetImage.SliceImage <- function(
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

#' @rdname GetImage
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
  } else if (length(x = new.names) != length(x = Cells(x = object))) {
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
