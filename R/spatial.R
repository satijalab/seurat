#' @include objects.R
#' @include generics.R
#' @importFrom methods setClass setOldClass slot<- setAs
#' setMethod new
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

setOldClass(Classes = 'tbl_df')
setOldClass(Classes = c('rastergrob', 'scalefactors'))

SliceImage <- setClass(
  Class = 'SliceImage',
  slots = list(
    'image' = 'rastergrob',
    'scale.factors' = 'scalefactors',
    'coordinates' = 'data.frame'
  )
)

#' @rdname Cells
#' @method Cells SliceImage
#' @export
#'
Cells.SliceImage <- function(x) {
  return(rownames(x = GetTissueCoordinates(object = x, scale = NULL)))
}

#' @method dim SliceImage
#' @export
#'
dim.SliceImage <- function(x) {
  return(dim(x = GetImage(object = x)$raster))
}

#' @method dimnames SliceImage
#' @export
#'
dimnames.SliceImage <- function(x) {
  ''
}

SpatialAssay <- setClass(
  Class = 'SpatialAssay',
  contains = 'Assay',
  slots = list(
    'images' = 'list'
  )
)

setAs(
  from = 'Assay',
  to = 'SpatialAssay',
  def = function(from) {
    object.list <- sapply(
      X = slotNames(x = from),
      FUN = slot,
      object = from,
      simplify = FALSE,
      USE.NAMES = TRUE
    )
    object.list <- c(
      list(
        'Class' = 'SpatialAssay'
        # 'image' = tibble()
      ),
      object.list
    )
    return(do.call(what = 'new', args = object.list))
  }
)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @importFrom crosstalk SharedData bscols
#' @importFrom plotly layout plot_ly highlight
#'
LinkedFeaturePlot <- function(
  object,
  feature,
  dims = 1:2,
  reduction = NULL,
  spatial = NULL,
  slice = NULL,
  slot = 'data',
  min.cutoff = NA,
  max.cutoff = NA
) {
  if (length(x = feature) > 1) {
    stop("'LinkedFeaturePlot' currently only supports one feature", call. = FALSE)
  }
  expression.data <- FetchData(object = object, vars = feature)
  # Generate the SpatialFeaturePlot and pull plotting data
  spatial.plot <- SpatialFeaturePlot(
    object = object,
    features = feature,
    spatial = spatial,
    slice = slice,
    slot = 'data',
    min.cutoff = min.cutoff,
    max.cutoff = max.cutoff
  )
  spatial.build <- GGpointToPlotlyBuild(
    plot = spatial.plot,
    information = expression.data
  )
  # Generate the normal FeaturePlot and pull plotting data
  feature.plot <- FeaturePlot(
    object = object,
    features = feature,
    dims = dims,
    reduction = reduction,
    min.cutoff = min.cutoff,
    max.cutoff = max.cutoff
  )
  feature.build <- GGpointToPlotlyBuild(
    plot = feature.plot,
    information = expression.data
  )
  # Generate full build
  plot.build <- merge(x = spatial.build, y = feature.build, by = 0)
  rownames(x = plot.build) <- plot.build$Row.names
  plot.build <- plot.build[, which(x = colnames(x = plot.build) != 'Row.names'), drop = FALSE]
  # Build the plotly plots
  plot.build <- SharedData$new(data = plot.build)
  spatial.plotly <- layout(
    p = plot_ly(
      data = plot.build,
      x = ~x.x,
      y = ~y.x,
      type = 'scatter',
      mode = 'markers',
      color = ~I(color.x),
      hoverinfo = 'text',
      text = ~feature.x
    )
  )
  spatial.plotly <- highlight(p = spatial.plotly, on = 'plotly_selected')
  feature.plotly <- layout(
    p = plot_ly(
      data = plot.build,
      x = ~y.x,
      y = ~y.y,
      type = 'scatter',
      mode = 'markers',
      color = ~I(color.y),
      hoverinfo = 'text',
      text = ~feature.y
    )
  )
  feature.plotly <- highlight(p = feature.plotly, on = 'plotly_selected')
  bscols(spatial.plotly, feature.plotly)
}


#' @importFrom crosstalk SharedData bscols
#' @importFrom plotly layout plot_ly highlight
#'
LinkedDimPlot <- function(
  object,
  dims = 1:2,
  reduction = NULL,
  spatial = NULL,
  slice = NULL,
  group.by = NULL
) {
  # Generate the SpatialFeaturePlot and pull plotting data
  spatial.plot <- SpatialDimPlot(
    object = object,
    group.by = group.by,
    spatial = spatial,
    slice = slice
  )
  spatial.build <- GGpointToPlotlyBuild(
    plot = spatial.plot
  )
  # Generate the normal FeaturePlot and pull plotting data
  feature.plot <- DimPlot(
    object = object,
    group.by = group.by,
    dims = dims,
    reduction = reduction
  )
  feature.build <- GGpointToPlotlyBuild(
    plot = feature.plot
  )
  # Generate full build
  plot.build <- merge(x = spatial.build, y = feature.build, by = 0)
  rownames(x = plot.build) <- plot.build$Row.names
  plot.build <- plot.build[, which(x = colnames(x = plot.build) != 'Row.names'), drop = FALSE]
  # Build the plotly plots
  plot.build <- SharedData$new(data = plot.build)
  spatial.plotly <- layout(
    p = plot_ly(
      data = plot.build,
      x = ~x.x,
      y = ~y.x,
      type = 'scatter',
      mode = 'markers',
      color = ~I(color.x),
      hoverinfo = 'text',
      text = ~feature.x
    )
  )
  spatial.plotly <- highlight(p = spatial.plotly, on = 'plotly_selected')
  feature.plotly <- layout(
    p = plot_ly(
      data = plot.build,
      x = ~y.x,
      y = ~y.y,
      type = 'scatter',
      mode = 'markers',
      color = ~I(color.y),
      hoverinfo = 'text',
      text = ~feature.y
    )
  )
  feature.plotly <- highlight(p = feature.plotly, on = 'plotly_selected')
  bscols(spatial.plotly, feature.plotly)
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
  slice = NULL,
  filter.matrix = TRUE,
  ...
) {
  data <- Read10X_h5(
    filename = file.path(data.dir, 'filtered_feature_bc_matrix.h5'),
    ...
  )
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
  object <- CreateSeuratObject(counts = data, assay = assay)
  object <- SetSlice(
    object = object,
    image = new(
      Class = 'SliceImage',
      image = rasterGrob(
        image = image,
        width = unit(x = 1, units = 'npc'),
        height = unit(x = 1, units = 'npc')
      ),
      scale.factors = scalefactors(
        spot = scale.factors$tissue_hires_scalef,
        fiducial = scale.factors$fiducial_diameter_fullres,
        hires = scale.factors$tissue_hires_scalef,
        scale.factors$tissue_lowres_scalef
      ),
      coordinates = tissue.positions
    ),
    slice = slice
  )
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

SpatialColors <- colorRampPalette(rev(RColorBrewer::brewer.pal(11, "Spectral")))

#' @importFrom ggplot2 ggplot aes aes_string geom_point xlim ylim
#' coord_cartesian labs theme_void
#'
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
    geom_spatial(data = image.tibble, aes(grob = grob), x = 0.5, y = 0.5) +
    geom_point(
      mapping = aes_string(
        x = colnames(x = data)[2],
        y = colnames(x = data)[1],
        color = col.by %iff% paste0("`", col.by, "`")
        ),
      size = pt.size
    ) +
    xlim(0, image.tibble$image_width) +
    ylim(image.tibble$image_height, 0) +
    coord_cartesian(expand = FALSE) +
    labs(color = NULL)
   plot <- plot + theme_void()
  return(plot)
}

SpatialDimPlot <- function(
  object,
  spatial = NULL,
  group.by = NULL,
  slice = NULL,
  combine = TRUE,
  ...
) {
  spatial = spatial %||% FilterObjects(object = object, classes.keep = 'SpatialAssay')
  if (length(x = spatial) != 1) {
    stop("Could not unabiguously find a SpatialAssay, please provide", call. = FALSE)
  }
  object[['ident']] <- Idents(object = object)
  group.by <- group.by %||% 'ident'
  data <- object[[group.by]]
  for (group in group.by) {
    if (!is.factor(x = data[, group])) {
      data[, group] <- factor(x = data[, group])
    }
  }
  slice <- slice %||% Slices(object = object[[spatial]])
  plots <- vector(
    mode = "list",
    length = length(x = group.by)
  )
  for (i in 1:length(x = group.by)) {
    group <- group.by[i]
    slice.plots <- vector(mode = "list",length = length(x = slice))
    for (s in 1:length(x = slice)) {
      coordinates <- GetTissueCoordinates(object = object, assay = spatial, slice = slice[s])
      plot.slice <- GetSlice(object = object[[spatial]], slice = slice[s])
      plot <- SingleSpatialPlot(
        data = cbind(
          coordinates,
          data[rownames(x = coordinates), group, drop = FALSE]
        ),
        col.by = group,
        image.tibble = tibble(
          grob = list(GetImage(object = plot.slice)),
          image_width = ncol(x = plot.slice),
          image_height = nrow(x = plot.slice)
        )
      )
      if (i == 1) {
        plot <- plot + ggtitle(label = slice[s]) +
          theme(plot.title = element_text(hjust = 0.5))
      }
      slice.plots[[s]] <- plot
    }
    plots[[i]] <- CombinePlots(plots = slice.plots, ncol = s)
  }
  if (combine) {
    plots <- CombinePlots(plots = plots, ncol = 1)
  }
  return(plots)
}

#' @importFrom tibble tibble
#'
SpatialFeaturePlot <- function(
  object,
  features,
  spatial = NULL,
  slice = NULL,
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
  spatial <- spatial %||% FilterObjects(object = object, classes.keep = 'SpatialAssay')
  if (length(x = spatial) != 1) {
    stop("Could not unabiguously find a SpatialAssay, please provide", call. = FALSE)
  }
  if (!inherits(x = object[[spatial]], what = 'SpatialAssay')) {
    stop("Assay ", spatial, " is not a SpatialAssay", call. = FALSE)
  }
  slice <- slice %||% Slices(object = object[[spatial]])
  features <- colnames(x = data)
  colnames(x = data) <- features
  rownames(x = data) <- colnames(object)
  plots <- vector(
    mode = "list",
    length = length(x = features)
  )
  for (i in 1:length(x = features)) {
    feature <- features[i]
    slice.plots <- vector(mode = "list",length = length(x = slice))
    for (s in 1:length(x = slice)) {
      coordinates <- GetTissueCoordinates(
        object = object,
        assay = spatial,
        slice = slice[s]
      )
      plot.slice <- GetSlice(object = object[[spatial]], slice = slice[s])
      plot <- SingleSpatialPlot(
        data = cbind(
          coordinates,
          data[rownames(x = coordinates), features, drop = FALSE]
        ),
        image.tibble = tibble(
          grob = list(GetImage(object = plot.slice)),
          image_width = ncol(x = plot.slice),
          image_height = nrow(x = plot.slice)
        ),
        col.by = feature
      )
      plot <- plot + scale_color_gradientn(name = feature, colours = SpatialColors(100))
      if (i == 1) {
        plot <- plot + ggtitle(label = slice[s]) + theme(plot.title = element_text(hjust = 0.5))
      }
      slice.plots[[s]] <- plot
    }
    plots[[i]] <- CombinePlots(plots = slice.plots, ncol = s)
  }
  if (combine) {
    plots <- CombinePlots(plots = plots, ncol = 1)
  }
  return(plots)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Methods for Seurat-defined generics
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @rdname GetImage
#' @method GetImage Seurat
#' @export
#'
GetImage.Seurat <- function(object, slice = NULL, assay = NULL, ...) {
  assay <- assay %||% DefaultAssay(object = object)
  if (!inherits(x = object[[assay]], what = 'SpatialAssay')) {
    stop("GetImage works only with SpatialAssay objects", call. = FALSE)
  }
  return(GetImage(object = object[[assay]], slice = slice))
}

#' @rdname GetImage
#' @method GetImage SliceImage
#' @export
#'
GetImage.SliceImage <- function(object, ...) {
  return(slot(object = object, name = 'image'))
}

#' @rdname GetImage
#' @method GetImage SpatialAssay
#' @export
#'
GetImage.SpatialAssay <- function(object, slice = NULL, ...) {
  slice <- slice %||% names(x = GetSlice(object = object))[1]
  return(GetImage(object = GetSlice(object = object, slice = slice)))
}

#' @rdname ImageData
#' @method GetSlice Seurat
#' @export
#'
GetSlice.Seurat <- function(object, slice = NULL, assay = NULL, ...) {
  assay <- assay %||% DefaultAssay(object = object)
  if (!inherits(x = object[[assay]], what = 'SpatialAssay')) {
    stop("'GetSlice' needs to be run on assays with spatial information", call. = FALSE)
  }
  return(GetSlice(object = object[[assay]], slice = slice))
}

#' @rdname ImageData
#' @method GetSlice SpatialAssay
#' @export
#'
GetSlice.SpatialAssay <- function(object, slice = NULL, ...) {
  images <- slot(object = object, name = 'images')
  if (is.null(x = slice)) {
    return(images)
  }
  return(images[[slice]])
}

#' @rdname GetTissueCoordinates
#' @method GetTissueCoordinates Seurat
#' @export
#'
GetTissueCoordinates.Seurat <- function(
  object,
  scale = 'lowres',
  slice = NULL,
  assay = NULL,
  ...
) {
  assay <- assay %||% DefaultAssay(object = object)
  if (!inherits(x = object[[assay]], what = 'SpatialAssay')) {
    stop("GetTissueCoordinates works only with SpatialAssay objects", call. = FALSE)
  }
  return(GetTissueCoordinates(object = object[[assay]], scale = scale, slice = slice, ...))
}

#' @rdname GetTissueCoordinates
#' @method GetTissueCoordinates SliceImage
#' @export
#'
GetTissueCoordinates.SliceImage <- function(object, scale = 'lowres', cols = c('imagerow', 'imagecol'), ...) {
  cols <- cols %||% colnames(x = slot(object = object, name = 'coordinates'))
  if (!is.null(x = scale)) {
    coordinates <- slot(object = object, name = 'coordinates')[, c('imagerow', 'imagecol')]
    scale <- match.arg(arg = scale, choices = c('spot', 'fiducial', 'hires', 'lowres'))
    scale.use <- slot(object = object, name = 'scale.factors')[[scale]]
    coordinates <- coordinates * scale.use
  } else {
    coordinates <- slot(object = object, name = 'coordinates')[, cols]
  }
  return(coordinates)
}

#' @rdname GetTissueCoordinates
#' @method GetTissueCoordinates SpatialAssay
#' @export
#'
GetTissueCoordinates.SpatialAssay <- function(
  object,
  scale = 'lowres',
  slice = NULL,
  ...
) {
  slice <- slice %||% names(x = GetSlice(object = object))[1]
  return(GetTissueCoordinates(
    object = GetSlice(object = object, slice = slice),
    scale = scale,
    ...
  ))
}

#' @rdname merge.Seurat
#' @export
#' @method merge SpatialAssay
#'
merge.SpatialAssay <- function(
  x = NULL,
  y = NULL,
  add.cell.ids = NULL,
  merge.data = TRUE,
  ...
) {
  assays <- c(x, y)
  merged.assay <- merge.Assay(
    x = x,
    y = y,
    add.cell.ids = add.cell.ids,
    merge.data = merge.data,
    ...
  )
  merged.assay <- as(object = merged.assay, Class = "SpatialAssay")
  cell.idx <- 1
  for(i in 1:length(x = assays)) {
    if (!isTRUE(x = all.equal(colnames(x = assays[[i]]), colnames(x = merged.assay)[cell.idx:(cell.idx + ncol(x = assays[[i]]) - 1)]))) {
      merged.assay <- RenameCells(object = assays[[i]], new.names = colnames(x = merged.assay)[cell.idx:(cell.idx + ncol(x = assays[[i]]) - 1)])
    }
    slices <- GetSlice(object = assays[[i]])
    for(s in 1:length(x = slices)) {
      if (names(x = slices)[s] %in% Slices(object = merged.assay)) {
        names(x = slices)[s] <- make.unique(names = c(Slices(object = merged.assay), names(x = slices)[s]))[length(x = Slices(object = merged.assay)) + 1]
      }
      merged.assay <- SetSlice(object = merged.assay, image = slices[[s]], slice = names(x = slices)[s])
    }
    cell.idx <- cell.idx + ncol(x = assays[[i]])
  }
  return(merged.assay)
}

RenameCells.SpatialAssay <- function(
  object,
  new.names = NULL,
  ...
) {
  old.names <- Cells(x = object)
  renamed.assay <- RenameCells.Assay(object = object, new.names = new.names)
  names(x = new.names) <- old.names
  for(s in Slices(object = object)) {
    coordinates <- GetTissueCoordinates(object = renamed.assay, slice = s, scale = NULL, cols = NULL)
    rownames(x = coordinates) <- new.names[rownames(x = coordinates)]
    renamed.assay <- SetTissueCoordinates(object = renamed.assay, coordinates = coordinates, slice = s)
  }
  return(renamed.assay)
}

#' @method ScaleFactors SliceImage
#' @export
#'
ScaleFactors.SliceImage <- function(object, ...) {
  return(slot(object = object, name = 'scale.factors'))
}

#' @rdname ImageData
#' @method SetSlice Assay
#' @export
#'
SetSlice.Assay <- function(object, image, slice = NULL, ...) {
  object <- as(object = object, Class = 'SpatialAssay')
  object <- SetSlice(object = object, image = image, slice = slice)
  return(object)
}

#' @rdname ImageData
#' @method SetSlice Seurat
#' @export
#'
SetSlice.Seurat <- function(object, image, slice = NULL, assay = NULL, ...) {
  spatial <- FilterObjects(object = object, classes.keep = 'SpatialAssay')
  assay <- assay %||% ifelse(
    test = length(x = spatial) > 0,
    yes = spatial[1],
    no = DefaultAssay(object = object)
  )
  object[[assay]] <- SetSlice(object = object[[assay]], image = image)
  return(object)
}

#' @rdname ImageData
#' @method SetSlice SpatialAssay
#' @export
#'
SetSlice.SpatialAssay <- function(object, image, slice = NULL, ...) {
  if (!all(Cells(x = image) %in% Cells(x = object))) {
    stop("Extra cells")
  }
  slice <- slice %||% paste0('slice', length(x = GetSlice(object = object)) + 1)
  slot(object = object, name = 'images')[[slice]] <- image
  return(object)
}

SetTissueCoordinates.Seurat <- function(object, coordinates, slice = NULL, assay = NULL, ...) {
  if (!inherits(x = object[[assay]], what = "SpatialAssay")) {
    stop("Assay must be a SpatialAssay.")
  }
  object[[assay]] <- SetTissueCoordinates(object = object[[assay]], slice = slice, coordinates = coordinates)
  return(object)
}

SetTissueCoordinates.SpatialAssay <- function(object, coordinates, slice = NULL, ...) {
  updated.slice <- GetSlice(object = object, slice = slice)
  updated.slice <- SetTissueCoordinates(object = updated.slice, coordinates = coordinates, ... )
  object <- SetSlice(object = object, image = updated.slice, slice = slice)
  return(object)
}

SetTissueCoordinates.SliceImage <- function(object, coordinates, ...) {
  slot(object = object, name = "coordinates") <- coordinates
  return(object)
}

Slices <- function(object) {
  return(names(x = slot(object = object, name = "images")))
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
