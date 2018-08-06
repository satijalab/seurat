#' Feature expression heatmap
#'
#' Draws a heatmap of single cell feature expression using superheat
#'
#' @param object Seurat object
#' @param features.use A vector of features to plot, defaults to \code{VariableFeatures(object = object)}
#' @param cells.use A vector of cells to plot
#' @param disp.min Minimum display value (all values below are clipped)
#' @param disp.max Maximum display value (all values above are clipped)
#' @param group.by Name of variable to group cells by
#' @param group.bar Add a color bar showing group status for cells
#' @param slot.use Data slot to use, choose from 'raw.data', 'data', or 'scale.data'
#' @param assay.use Assay to pull from
# @param check.plot Check that plotting will finish in a reasonable amount of time
#' @param ... Ignored for now
#'
#' @return Invisbly returns the final grob
#'
#' @importFrom scales hue_pal
#' @importFrom ggplot2 annotation_raster coord_cartesian ggplot_build
#' @export
#'
#' @examples
#' DoHeatmap(object = pbmc_small)
#'
DoHeatmap <- function(
  object,
  features.use = NULL,
  cells.use = NULL,
  disp.min = -2.5,
  disp.max = 2.5,
  group.by = "ident",
  slot.use = 'scale.data',
  group.bar = TRUE,
  # group.order = NULL,
  # draw.line = TRUE,
  assay.use = NULL,
  # check.plot = FALSE,
  ...
) {
  cells.use <- cells.use %||% colnames(x = object)
  assay.use <- assay.use %||% DefaultAssay(object = object)
  DefaultAssay(object = object) <- assay.use
  features.use <- features.use %||% VariableFeatures(object = object)
  disp.max <- ifelse(
    test = slot.use != 'scale.data',
    yes = max(disp.max, 10),
    no = disp.max
  )
  data.plot <- FetchData(
    object = object,
    vars.fetch = rev(x = features.use),
    cells.use = cells.use,
    slot = slot.use
  )
  group.by <- group.by %||% 'ident'
  group.use <- switch(
    EXPR = group.by,
    'ident' = Idents(object = object),
    object[group.by, drop = TRUE]
  )
  group.use <- factor(x = group.use[cells.use])
  p <- SingleRasterMap(
    data.plot = data.plot,
    disp.min = disp.min,
    disp.max = disp.max,
    cell.order = names(x = sort(x = group.use)),
    group.by = group.use
  )
  if (group.bar) {
    # TODO: Change group.bar to annotation.bar
    pbuild <- ggplot_build(plot = p)
    cols <- hue_pal()(length(x = levels(x = group.use)))
    names(x = cols) <- levels(x = group.use)
    y.pos <- max(pbuild$layout$panel_params[[1]]$y.range) + 0.25
    y.max <- y.pos + 0.5
    p <- p + annotation_raster(
      raster = t(x = cols[sort(x = group.use)]),
      xmin = -Inf,
      xmax = Inf,
      ymin = y.pos,
      ymax = y.max
    ) +
      coord_cartesian(ylim = c(0, y.max), clip = 'off')
  }
  p <- p + theme(line = element_blank())
  return(p)
}

#' Dimensional reduction heatmap
#'
#' Draws a heatmap focusing on a principal component. Both cells and genes are sorted by their
#' principal component scores. Allows for nice visualization of sources of heterogeneity in the dataset.
#'
#' @inheritParams DoHeatmap
#' @param assay.use A vector of assays to pull data from
#' @param reduction.use Which dimmensional reduction to use
#' @param dims.use Dimensions to plot
#' @param cells.use A list of cells to plot. If numeric, just plots the top cells.
#' @param num.features NUmber of genes to plot
#' @param do.balanced Plot an equal number of genes with both + and - scores.
#' @param num.col Number of columns to plot
#' @param combine.plots Combine plots into a single gg object; note that if TRUE; themeing will not work when plotting multiple dimensions
#' @param plot.method Either "quick" for fast exploratory visualization or "gg" for slower
#' plotting but returns ggplot objects that can be customized.
#'
#' @return Invisbly returns the final grob
#'
#' @importFrom utils menu
#' @export
#'
#' @seealso \code{\link{superheat::superheat}}
#'
#' @examples
#' DimHeatmap(object = pbmc_small)
#'
DimHeatmap <- function(
  object,
  assay.use = NULL,
  reduction.use = "pca",
  dims.use = 1,
  cells.use = NULL,
  num.features = 30,
  use.full = FALSE,
  disp.min = -2.5,
  disp.max = 2.5,
  slot.use = 'scale.data',
  do.balanced = FALSE,
  check.plot = TRUE,
  num.col = NULL,
  combine.plots = TRUE,
  plot.method = c("quick", "gg"),
  ...
) {
  plot.method <- match.arg(arg = plot.method)
  if (is.null(x = num.col)) {
    if (length(x = dims.use) > 2) {
      num.col <- 3
    } else {
      num.col <- length(x = dims.use)
    }
  }
  plot.list <- vector(mode = 'list', length = length(x = dims.use))
  assay.use <- assay.use %||% DefaultAssay(object = object)
  if (!DefaultAssay(object = object[[reduction.use]]) %in% assay.use) {
    warning("assay")
  }
  if (is.numeric(x = cells.use)) {
    cells.use <- lapply(
      X = dims.use,
      FUN = function(x){
        cells <- TopCells(
          object = object[[reduction.use]],
          dim.use = x,
          num.cells = cells.use,
          do.balanced = do.balanced
        )
        if (do.balanced) {
          cells$negative <- rev(x = cells$negative)
        }
        cells <- unlist(x = unname(obj = cells))
      }
    )
  }
  cells.use <- cells.use %||% colnames(x = object)
  if (!is.list(x = cells.use)) {
    cells.use <- lapply(X = dims.use, FUN = function(...) cells.use)
  }
  features <- lapply(
    X = dims.use,
    FUN = TopFeatures,
    object = object[[reduction.use]],
    num.features = num.features,
    do.balanced = do.balanced
  )
  features.all <- unique(x = unlist(x = features))
  if (length(x = assay.use) > 1) {
    features.keyed <- lapply(
      X = assay.use,
      FUN = function(assay.name) {
        features.use <- features.all[features.all %in% rownames(x = object[[assay.name]])]
        if (length(x = features.use) > 0) {
          return(paste0(Key(object = object[[assay.name]]), features.use))
        }
      }
    )
    features.keyed <- Filter(f = Negate(f = is.null), x = features.keyed)
    features.keyed <- unlist(x = features.keyed)
  } else {
    features.keyed <- features.all
  }
  data.all <- FetchData(
    object = object,
    vars.fetch = features.keyed,
    cells.use = unique(x = unlist(x = cells.use)),
    slot = slot.use
  )
  data.all <- MinMax(data = data.all, min = disp.min, max = disp.max)
  data.limits <- c(min(data.all), max(data.all))
  # if (check.plot && any(c(length(x = features.keyed), length(x = cells.use[[1]])) > 700)) {
  #   choice <- menu(c("Continue with plotting", "Quit"), title = "Plot(s) requested will likely take a while to plot.")
  #   if (choice != 1) {
  #     return(invisible(x = NULL))
  #   }
  # }
  if (plot.method == "quick") {
    num.row <- floor(x = length(x = dims.use) / 3.01) + 1
    orig_par <- par()$mfrow
    par(mfrow = c(num.row, min(length(x = dims.use), 3)))
  }
  for (i in 1:length(x = dims.use)) {
    dim.features <- unname(obj = unlist(x = rev(x = features[[i]])))
    dim.features <- rev(x = unlist(x = lapply(
      X = dim.features,
      FUN = function(feat) {
        return(grep(pattern = paste0(feat, '$'), x = features.keyed, value = TRUE))
      }
    )))
    data.plot <- data.all[cells.use[[i]], dim.features]
    if (plot.method == "quick") {
      SingleImageMap(data.plot = data.plot, plot.title = paste0(Key(object = object[[reduction.use]]), dims.use[i]))
    } else {
      plot.list[[i]] <- SingleRasterMap(data.plot = data.plot, limits = data.limits)
    }
  }
  if (plot.method == "quick") {
    return(invisible(x = NULL))
  }
  if (combine.plots) {
    plot.list <- CombinePlots(
      plot.list = plot.list,
      num.col = num.col,
      legend.position = 'right'
    )
  }
  return(plot.list)
}

#' Single cell ridge plot
#'
#' Draws a ridge plot of single cell data (gene expression, metrics, PC
#' scores, etc.)
#'
#' @param object Seurat object
#' @param features.plot Features to plot (gene expression, metrics, PC scores,
#' anything that can be retreived by FetchData)
#' @param ident.include Which classes to include in the plot (default is all)
#' @param num.col Number of columns if multiple plots are displayed
#' @param do.sort Sort identity classes (on the x-axis) by the average
#' expression of the attribute being potted
#' @param y.max Maximum y axis value
#' @param same.y.lims Set all the y-axis limits to the same values
#' @param cols.use Colors to use for plotting
#' @param group.by Group (color) cells in different ways (for example, orig.ident)
#' @param y.log plot Y axis on log scale
#' @param combine.plots Combine plots into a single gg object; note that if TRUE; themeing will not work when plotting multiple features
#' @param plot.counts Use non-normalized counts data for plotting
#' @param \dots additional parameters to pass to FetchData (for example, use.imputed, use.scaled, use.raw)
#'
#' @return A ggplot object
#'
#' @export
#'
#' @examples
#' RidgePlot(object = pbmc_small, features.plot = 'PC1')
#'
RidgePlot <- function(
  object,
  features.plot,
  ident.include = NULL,
  num.col = NULL,
  do.sort = FALSE,
  y.max = NULL,
  same.y.lims = FALSE,
  cols.use = NULL,
  group.by = NULL,
  y.log = FALSE,
  combine.plots = TRUE,
  plot.counts = FALSE,
  ...
) {
  return(ExIPlot(
    object = object,
    plot.type = 'ridge',
    features.plot = features.plot,
    ident.include = ident.include,
    num.col = num.col,
    do.sort = do.sort,
    y.max = y.max,
    same.y.lims = same.y.lims,
    cols.use = cols.use,
    group.by = group.by,
    y.log = y.log,
    combine.plots = combine.plots,
    plot.counts = plot.counts,
    ...
  ))
}

#' Single cell violin plot
#'
#' Draws a violin plot of single cell data (gene expression, metrics, PC
#' scores, etc.)
#'
#' @inheritParams RidgePlot
#' @param adjust.use Adjust parameter for geom_violin
#' @param point.size.use Point size for geom_violin
#'
#' @return A ggplot object
#'
#' @export
#'
#' @examples
#' VlnPlot(object = pbmc_small, features.plot = 'PC1')
#'
VlnPlot <- function(
  object,
  features.plot,
  ident.include = NULL,
  num.col = NULL,
  do.sort = FALSE,
  y.max = NULL,
  same.y.lims = FALSE,
  adjust.use = 1,
  point.size.use = 1,
  cols.use = NULL,
  group.by = NULL,
  y.log = FALSE,
  combine.plots = TRUE,
  plot.counts = FALSE,
  ...
) {
  return(ExIPlot(
    object = object,
    plot.type = 'violin',
    features.plot = features.plot,
    ident.include = ident.include,
    num.col = num.col,
    do.sort = do.sort,
    y.max = y.max,
    same.y.lims = same.y.lims,
    adjust.use = adjust.use,
    point.size.use = point.size.use,
    cols.use = cols.use,
    group.by = group.by,
    y.log = y.log,
    combine.plots = combine.plots,
    plot.counts = plot.counts,
    ...
  ))
}

#' Dot plot visualization
#'
#' Intuitive way of visualizing how feature expression changes across different
#' identity classes (clusters). The size of the dot encodes the percentage of
#' cells within a class, while the color encodes the AverageExpression level of
#' cells within a class (blue is high).
#'
#' @param object Seurat object
#' @param features.plot Input vector of features
#' @param cols.use Colors to plot, can pass a single character giving the name of
#' a palette from \code{RColorBrewer::brewer.pal.info}
#' @param col.min Minimum scaled average expression threshold (everything smaller
#'  will be set to this)
#' @param col.max Maximum scaled average expression threshold (everything larger
#' will be set to this)
#' @param dot.min The fraction of cells at which to draw the smallest dot
#' (default is 0). All cell groups with less than this expressing the given
#' gene will have no dot drawn.
#' @param dot.scale Scale the size of the points, similar to cex
#' @param scale.by Scale the size of the points by 'size' or by 'radius'
#' @param scale.min Set lower limit for scaling, use NA for default
#' @param scale.max Set upper limit for scaling, use NA for default
#' @param group.by Factor to group the cells by
#' @param split.by Factor to split the groups by (replicates the functionality of the old SplitDotPlotGG)
#'
#' @return A ggplot object
#'
#' @importFrom ggplot2 ggplot aes_string scale_size scale_radius geom_point theme element_blank
#' scale_color_identity scale_color_distiller scale_color_gradient guides guide_legend guide_colorbar
#' @export
#'
#' @aliases SplitDotPlotGG
#' @seealso \code{RColorBrewer::brewer.pal.info}
#'
#' @examples
#' cd_genes <- c("CD247", "CD3E", "CD9")
#' DotPlot(object = pbmc_small, features.plot = cd_genes)
#' pbmc_small['groups'] <- sample(x = c('g1', 'g2', size = ncol(x = pbmc_small), replace = TRUE))
#' DotPlot(object = pbmc_small, features.plot = cd_genes, split.by = 'groups')
#'
DotPlot <- function(
  object,
  features.plot,
  cols.use = c("lightgrey", "blue"),
  col.min = -2.5,
  col.max = 2.5,
  dot.min = 0,
  dot.scale = 6,
  scale.by = 'radius',
  scale.min = NA,
  scale.max = NA,
  group.by = NULL,
  split.by = NULL
) {
  scale.func <- switch(
    EXPR = scale.by,
    'size' = scale_size,
    'radius' = scale_radius,
    stop("'scale.by' must be either 'size' or 'radius'")
  )
  data.features <- FetchData(object = object, vars.fetch = features.plot)
  data.features$id <- if (is.null(x = group.by)) {
    Idents(object = object)
  } else {
    object[group.by, drop = TRUE]
  }
  data.features$id <- as.vector(x = data.features$id)
  if (!is.null(x = split.by)) {
    splits <- object[split.by, drop = TRUE]
    if (length(x = unique(x = splits)) > length(x = cols.use)) {
      stop("Not enought colors for the number of groups")
    }
    cols.use <- cols.use[1:length(x = unique(x = splits))]
    names(x = cols.use) <- unique(x = splits)
    data.features$id <- paste(data.features$id, splits, sep = '_')
  }
  data.plot <- lapply(
    X = unique(x = data.features$id),
    FUN = function(ident) {
      data.use <- data.features[data.features$id == ident, 1:(ncol(x = data.features) - 1), drop = FALSE]
      avg.exp <- apply(
        X = data.use,
        MARGIN = 2,
        FUN = function(x) {
          return(mean(x = expm1(x = x)))
        }
      )
      pct.exp <- apply(X = data.use, MARGIN = 2, FUN = PercentAbove, threshold = 0)
      return(list(avg.exp = avg.exp, pct.exp = pct.exp))
    }
  )
  names(x = data.plot) <- unique(x = data.features$id)
  data.plot <- lapply(
    X = names(x = data.plot),
    FUN = function(x) {
      data.use <- as.data.frame(x = data.plot[[x]])
      data.use$features.plot <- rownames(x = data.use)
      data.use$id <- x
      return(data.use)
    }
  )
  data.plot <- do.call(what = 'rbind', args = data.plot)
  avg.exp.scaled <- sapply(
    X = unique(x = data.plot$features.plot),
    FUN = function(x) {
      data.use <- data.plot[data.plot$features.plot == x, 'avg.exp']
      data.use <- scale(x = data.use)
      data.use <- MinMax(data = data.use, min = col.min, max = col.max)
      return(data.use)
    }
  )
  avg.exp.scaled <- as.vector(x = t(x = avg.exp.scaled))
  if (!is.null(x = split.by)) {
    avg.exp.scaled <- as.numeric(x = cut(x = avg.exp.scaled, breaks = 20))
  }
  data.plot$avg.exp.scaled <- avg.exp.scaled
  data.plot$features.plot <- factor(
    x = data.plot$features.plot,
    levels = rev(x = features.plot)
  )
  data.plot$pct.exp[data.plot$pct.exp < dot.min] <- NA
  data.plot$pct.exp <- data.plot$pct.exp * 100
  if (!is.null(x = split.by)) {
    splits.use <- vapply(
      X = strsplit(x = data.plot$id, split = '_'),
      FUN = '[[',
      FUN.VALUE = character(length = 1L),
      2
    )
    data.plot$colors <- mapply(
      FUN = function(color, value) {
        return(colorRampPalette(colors = c('grey', color))(20)[value])
      },
      color = cols.use[splits.use],
      value = avg.exp.scaled
    )
  }
  color.by <- ifelse(test = is.null(x = split.by), yes = 'avg.exp.scaled', no = 'colors')

  p <- ggplot(data = data.plot, mapping = aes_string(x = 'features.plot', y = 'id')) +
    geom_point(mapping = aes_string(size = 'pct.exp', color = color.by)) +
    scale.func(range = c(0, dot.scale), limits = c(scale.min, scale.max)) +
    theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
    guides(size = guide_legend(title = 'Percent Expressed'))
  if (!is.null(x = split.by)) {
    p <- p + scale_color_identity()
  } else if (length(x = cols.use) == 1) {
    p <- p + scale_color_distiller(palette = cols.use)
  } else {
    p <- p + scale_color_gradient(low = cols.use[1], high = cols.use[2])
  }
  if (is.null(x = split.by)) {
    p <- p + guides(color = guide_colorbar(title = 'Average Expression'))
  }
  return(p)
}

#' Dimensional reduction plot
#'
#' Graphs the output of a dimensional reduction technique (PCA by default).
#' Cells are colored by their identity class.
#'
#' @param object Seurat object
#' @param reduction.use Which dimensionality reduction to use. Default is
#' "pca", can also be "tsne", or "ica", assuming these are precomputed.
#' @param dims.use Dimensions to plot, must be a two-length numeric vector specifying x- and y-dimensions
#' @param cells.use Vector of cells to plot (default is all cells)
#' @param pt.size Adjust point size for plotting
#' @param cols.use Vector of colors, each color corresponds to an identity class. By default, ggplot2 assigns colors.
#' @param group.by Group (color) cells in different ways (for example, orig.ident)
#' @param shape.by If NULL, all points are circles (default). You can specify any
#' cell attribute (that can be pulled with FetchData) allowing for both
#' different colors and different shapes on cells.
# @param do.hover Enable hovering over points to view information
# @param data.hover Data to add to the hover, pass a character vector of
# features to add. Defaults to cell name and ident. Pass 'NULL' to clear extra
# information.
# @param do.identify Opens a locator session to identify clusters of cells.
#' @param plot.order Specify the order of plotting for the idents. This can be
#' useful for crowded plots if points of interest are being buried. Provide
#' either a full list of valid idents or a subset to be plotted last (on top).
#' @param do.label Whether to label the clusters
#' @param label.size Sets size of labels
#' @param cells.highlight A list of character or numeric vectors of cells to
#' highlight. If only one group of cells desired, can simply
#' pass a vector instead of a list. If set, colors selected cells to the color(s)
#' in \code{cols.highlight} and other cells black (white if dark.theme = TRUE);
#'  will also resize to the size(s) passed to \code{sizes.highlight}
#' @param cols.highlight A vector of colors to highlight the cells as; will
#' repeat to the length groups in cells.highlight
#' @param sizes.highlight Size of highlighted cells; will repeat to the length
#' groups in cells.highlight
# @param vector.friendly FALSE by default. If TRUE, points are flattened into
# a PNG, while axes/labels retain full vector resolution. Useful for producing
# AI-friendly plots with large numbers of cells.
# @param png.file Used only if vector.friendly is TRUE. Location for temporary
# PNG file.
# @param png.arguments Used only if vector.friendly is TRUE. Vector of three
# elements (PNG width, PNG height, PNG DPI) to be used for temporary PNG.
# Default is c(10,10,100)
#' @param na.value Color value for NA points when using custom scale.
# @param ... Extra parameters to FeatureLocator for do.identify = TRUE
#'
#' @return A ggplot object
#'
#' @export
#'
#' @seealso \code{\link{FeaturePlot}} \code{\link{FeatureMap}}
#'
#' @examples
#' DimPlot(object = pbmc_small)
#'
DimPlot <- function(
  object,
  reduction.use = 'pca',
  dims.use = c(1, 2),
  group.by = NULL,
  cols.use = NULL,
  cells.use = NULL,
  pt.size = 1,
  shape.by = NULL,
  plot.order = NULL,
  do.label = FALSE,
  label.size = 4,
  cells.highlight = NULL,
  cols.highlight = 'red',
  sizes.highlight = 1,
  na.value = 'grey50',
  # ...,
  do.hover = FALSE,
  data.hover = 'ident',
  do.identify = FALSE,
  no.legend = FALSE,
  vector.friendly = FALSE,
  png.file = NULL,
  png.arguments = c(10,10, 100),
  ...
) {
  ReadPlotParams(object)
  if (length(x = dims.use) != 2) {
    stop("'dims.use' must be a two-length vector")
  }
  cells.use <- cells.use %||% colnames(x = object)
  data.plot <- Embeddings(object = object[[reduction.use]])[cells.use, dims.use]
  data.plot <- as.data.frame(x = data.plot)
  dims.use <- paste0(Key(object = object[[reduction.use]]), dims.use)
  group.by <- group.by %||% 'ident'
  data.plot[, group.by] <- switch(
    EXPR = group.by,
    'ident' = Idents(object = object),
    object[group.by, drop = TRUE]
  )
  if (!is.null(x = shape.by)) {
    data.plot[, shape.by] <- object[shape.by, drop = TRUE]
  }
  plot <- SingleDimPlot(
    data.plot = data.plot,
    dims.use = dims.use,
    col.by = group.by,
    cols.use = cols.use,
    pt.size = pt.size,
    shape.by = shape.by,
    plot.order = plot.order,
    do.label = do.label,
    label.size = label.size,
    cells.highlight = cells.highlight,
    cols.highlight = cols.highlight,
    sizes.highlight = sizes.highlight
  )
  return(plot)
}

#' Visualize 'features' on a dimensional reduction plot
#'
#' Colors single cells on a dimensional reduction plot according to a 'feature'
#' (i.e. gene expression, PC scores, number of genes detected, etc.)
#'
#' @inheritParams DimPlot
#' @param features.plot Vector of features to plot
#' @param min.cutoff Vector of minimum cutoff values for each feature, may specify quantile in the form of 'q##' where '##' is the quantile (eg, 1, 10)
#' @param max.cutoff Vector of maximum cutoff values for each feature, may specify quantile in the form of 'q##' where '##' is the quantile (eg, 1, 10)
#' @param cols.use The two colors to form the gradient over. Provide as string vector with
#' the first color corresponding to low values, the second to high. Also accepts a Brewer
#' color scale or vector of colors. Note: this will bin the data into number of colors provided.
# @param pch.use Pch for plotting
# @param overlay Plot two features overlayed one on top of the other
#' @param reduction.use Which dimensionality reduction to use. Default is
#' "tsne", can also be "pca", or "ica", assuming these are precomputed.
#' @param combine.plots Combine plots into a single gg object; note that if TRUE; themeing will not work when plotting multiple features
#'
#' @return A ggplot object
#'
#' @importFrom RColorBrewer brewer.pal.info
#' @importFrom ggplot2 scale_color_gradientn labs
#' @export
#'
#' @seealso \code{\link{DimPlot}} \code{\link{FeatureMap}}
#'
#' @examples
#' FeaturePlot(object = pbmc_small, features.plot = 'PC1')
#'
FeaturePlot <- function(
  object,
  features.plot,
  min.cutoff = NA,
  max.cutoff = NA,
  reduction.use = 'tsne',
  dims.use = c(1, 2),
  group.by = NULL,
  cols.use = c("yellow", "red"),
  cells.use = NULL,
  pt.size = 1,
  shape.by = NULL,
  plot.order = NULL,
  do.label = FALSE,
  label.size = 4,
  na.value = 'grey50',
  combine.plots = TRUE
) {
  ReadPlotParams(object)
  if (length(x = dims.use) != 2 || !is.numeric(x = dims.use)) {
    stop("'dims.use' must be a two-length integer vector")
  }
  dims.use <- paste0(Key(object = object[[reduction.use]]), dims.use)
  num.col <- NULL
  if (is.null(x = num.col)) {
    num.col <- 2
    if (length(x = features.plot) == 1) {
      num.col <- 1
    }
    if (length(x = features.plot) > 6) {
      num.col <- 3
    }
    if (length(x = features.plot) > 9) {
      num.col <- 4
    }
  }
  cells.use <- cells.use %||% colnames(object)
  data.features <- FetchData(
    object = object,
    vars.fetch = c(dims.use, features.plot),
    cells.use = cells.use
  )
  features.plot <- colnames(x = data.features)[3:ncol(x = data.features)]
  min.cutoff <- mapply(
    FUN = function(cutoff, feature) {
      ifelse(
        test = is.na(x = cutoff),
        yes = min(data.features[, feature]),
        no = cutoff
      )
    },
    cutoff = min.cutoff,
    feature = features.plot
  )
  max.cutoff <- mapply(
    FUN = function(cutoff, feature) {
      ifelse(
        test = is.na(x = cutoff),
        yes = max(data.features[, feature]),
        no = cutoff
      )
    },
    cutoff = max.cutoff,
    feature = features.plot
  )
  check.lengths <- unique(x = vapply(
    X = list(features.plot, min.cutoff, max.cutoff),
    FUN = length,
    FUN.VALUE = numeric(length = 1)
  ))
  if (length(x = check.lengths) != 1) {
    stop('There must be the same number of minimum and maximum cuttoffs as there are features')
  }
  brewer.gran <- ifelse(
    test = length(x = cols.use) == 1,
    yes = brewer.pal.info[cols.use, ]$maxcolors,
    no = length(x = cols.use)
  )
  data.features[, 3:ncol(x = data.features)] <- sapply(
    X = 3:ncol(x = data.features),
    FUN = function(index) {
      data.feature <- as.vector(x = data.features[, index])
      min.use <- SetQuantile(cutoff = min.cutoff[index - 2], data.feature)
      max.use <- SetQuantile(cutoff = max.cutoff[index - 2], data.feature)
      data.feature[data.feature < min.use] <- min.use
      data.feature[data.feature > max.use] <- max.use
      if (brewer.gran == 2) {
        return(data.feature)
      }
      data.cut <- if (all(data.feature == 0)) {
        0
      } else {
        as.numeric(x = as.factor(x = cut(
          x = as.numeric(x = data.feature),
          breaks = brewer.gran
        )))
      }
      return(data.cut)
    }
  )
  colnames(x = data.features)[3:ncol(x = data.features)] <- features.plot
  rownames(x = data.features) <- cells.use
  plots <- vector(mode = 'list', length = length(x = features.plot))
  for (feature in features.plot) {
    p <- SingleDimPlot(
      data.plot = data.features[, c(dims.use, feature)],
      dims.use = dims.use,
      col.by = feature,
      cols.use = cols.use,
      pt.size = pt.size,
      do.label = do.label,
      label.size = label.size
    ) + labs(title = feature)
    if (brewer.gran == 2) {
      suppressMessages(expr = p <- p + scale_color_gradientn(
        colors = cols.use,
        guide = 'colorbar'
      ))
    }
    plots[[which(x = features.plot == feature)]] <- p
  }
  if (combine.plots) {
    plots <- CombinePlots(plot.list = plots, num.col = num.col, legend.position = 'none')
  }
  return(plots)
}

#' Vizualization of multiple features
#'
#' Similar to FeaturePlot, however, also splits the plot by visualizing each
#' identity class separately.
#'
#' Particularly useful for seeing if the same groups of cells co-exhibit a
#' common feature (i.e. co-express a gene), even within an identity class. Best
#' understood by example.
#'
#' @inheritParams FeaturePlot
#' @param group.display Which identity classes to display (default is all identity classes)
#' @param slot.use Dataset to use for plotting, choose from 'raw.data', 'data', or 'scale.data'
#' @param scale.group Scale each group separately. Default is FALSE.
#' @param horizontal rotate the plot such that the features are columns, groups are the rows
#'
#' @return A ggplot object
#'
#' @importFrom rlang !! sym
#' @importFrom ggplot2 facet_grid vars theme_bw scale_color_gradient guide_colorbar
#' @export
#'
#' @aliases FeatureHeatmap
#' @seealso \code{\link{FeaturePlot}} \code{\link{DimPlot}}
#'
#' @examples
#' pbmc_small
#' FeatureMap(object = pbmc_small, features.plot = "PC1")
#'
FeatureMap <- function(
  object,
  features.plot,
  min.cutoff = NA,
  max.cutoff = NA,
  reduction.use = "tsne",
  dims.use = c(1, 2),
  group.by = NULL,
  cols.use = c("grey", "red"),
  pt.size = 2,
  group.display = NULL,
  slot.use = 'data',
  scale.group = FALSE,
  horizontal = FALSE
) {
  if (length(x = dims.use) != 2 || !is.numeric(x = dims.use)) {
    stop("'dims.use' must be a two-length integer vector")
  }
  dims.use <- paste0(Key(object = object[[reduction.use]]), dims.use)
  idents.use <- if (is.null(x = group.by)) {
    Idents(object = object)
  } else {
    object[group.by, drop = TRUE]
  }
  idents.use <- as.vector(x = idents.use)
  names(x = idents.use) <- colnames(x = object)
  group.display <- group.display %||% unique(x = idents.use)
  group.display <- intersect(x = group.display, y = idents.use)
  if (length(x = group.display) == 0) {
    stop("No groups selected to display present in the object")
  }
  idents.use <- Filter(f = function(x) x %in% group.display, x = idents.use)
  data.plot <- FetchData(
    object = object,
    vars.fetch = features.plot,
    slot = slot.use,
    cells.use = names(x = idents.use)
  )
  data.plot <- data.frame(
    expression = unlist(x = data.plot, use.names = FALSE),
    feature = unlist(x = lapply(X = features.plot, FUN = rep.int, times = length(x = idents.use))),
    cell = rep.int(x = names(x = idents.use), times = length(x = features.plot)),
    stringsAsFactors = FALSE
  )
  data.plot$ident <- idents.use[data.plot$cell]
  data.plot$scaled.expression <- if (scale.group) {
    unlist(x = apply(
      X = expand.grid(unique(x = data.plot$feature), unique(x = data.plot$ident)),
      MARGIN = 1,
      FUN = function(x) {
        return(scale(x = data.plot[data.plot$feature == x[1] & data.plot$ident == x[2], 'expression']))
      }
    ))
  } else {
    unlist(x = lapply(
      X = unique(x = data.plot$feature),
      FUN = function(x) {
        return(scale(x = data.plot[data.plot$feature == x, 'expression']))
      }
    ))
  }
  data.plot <- suppressWarnings(expr = cbind(
    data.plot,
    Embeddings(object = object[[reduction.use]])[data.plot$cell, dims.use]
  ))
  min.cutoff <- ifelse(
    test = is.na(x = min.cutoff),
    yes = min(data.plot$scaled.expression),
    no = min.cutoff
  )
  max.cutoff <- ifelse(
    test = is.na(x = max.cutoff),
    yes = max(data.plot$scaled.expression),
    no = max.cutoff
  )
  min.cutoff <- SetQuantile(cutoff = min.cutoff, data = data.plot$scaled.expression)
  max.cutoff <- SetQuantile(cutoff = max.cutoff, data = data.plot$scaled.expression)
  data.plot$scaled.expression <- MinMax(
    data = data.plot$scaled.expression,
    min = min.cutoff,
    max = max.cutoff
  )
  plot <- SingleDimPlot(
    data.plot = data.plot,
    dims.use = dims.use,
    col.by = 'scaled.expression',
    pt.size = pt.size
  )
  if (horizontal) {
    row <- 'feature'
    col <- 'ident'
  } else {
    row <- 'ident'
    col <- 'feature'
  }
  plot <- plot +
    facet_grid(rows = vars(!!sym(x = row)), cols = vars(!!sym(x = col))) +
    theme_bw() +
    NoGrid() +
    scale_color_gradient(
      low = cols.use[1],
      high = cols.use[2],
      guide = guide_colorbar(title = 'Scaled Expression')
    )
  return(plot)
}

#' Scatter plot of single cell data
#'
#' Creates a scatter plot of two features (typically gene expression), across a
#' set of single cells. Cells are colored by their identity class. Pearson
#' correlation between the two features is displayed above the plot.
#'
#' @param object Seurat object
#' @inheritParams FetchData
#' @param feature1 First feature to plot. Typically gene expression but can also
#' be metrics, PC scores, etc. - anything that can be retreived with FetchData
#' @param feature2 Second feature to plot.
#' @param cells.use Cells to include on the scatter plot.
#' @param cols.use Colors to use for identity class plotting.
#' @param pt.size Size of the points on the plot
#' @param shape.by Ignored for now
#' @param span Spline span in loess function call, if \code{NULL}, no spline added
#' @param plot.counts Use non-normalized counts data for plotting
#' @param ... Ignored for now
#'
#' @return A ggplot object
#'
#' @importFrom ggplot2 geom_smooth aes_string
#' @importFrom cowplot theme_cowplot
#' @export
#'
#' @aliases GenePlot
#'
#' @examples
#' FeatureScatter(object = pbmc_small, feature1 = 'CD9', feature2 = 'CD3E')
#'
FeatureScatter <- function(
  object,
  feature1,
  feature2,
  cells.use = NULL,
  cols.use = NULL,
  pt.size = 1,
  shape.by = NULL,
  span = NULL,
  plot.counts = FALSE,
  ...
) {
  slot.use = "data"
  if (plot.counts == TRUE) {
    slot.use = "counts"
  }
  cells.use <- cells.use %||% colnames(x = object)
  p <- SingleCorPlot(
    data.plot = FetchData(
      object = object,
      vars.fetch = c(feature1, feature2),
      cells.use = cells.use,
      slot = slot.use
    ),
    col.by = Idents(object = object)[cells.use],
    cols.use = cols.use,
    pt.size = pt.size,
    legend.title = 'Identity'
  )
  if (!is.null(x = span)) {
    p <- p + geom_smooth(
      mapping = aes_string(x = feature1, y = feature2),
      method = 'loess',
      span = span
    )
  }
  p <- p + theme_cowplot()
  return(p)
}

#' Cell-cell scatter plot
#'
#' Creates a plot of scatter plot of features across two single cells. Pearson
#' correlation between the two cells is displayed above the plot.
#'
#' @inheritParams FeatureScatter
#' @param cell1 Cell 1 name
#' @param cell2 Cell 2 name
#' @param features.use Features to plot (default, all features)
#' @param \dots Extra parameters passed to kde2d
#'
#' @return A ggplot object
#'
#' @export
#' @seealso \code{\link{MASS::kde2d}}
#'
#' @aliases CellPlot
#'
#' @examples
#' CellScatter(object = pbmc_small, cell1 = 'ATAGGAGAAACAGA', cell2 = 'CATCAGGATGCACA')
#'
CellScatter <- function(
  object,
  cell1,
  cell2,
  features.use = NULL,
  cols.use = NULL,
  pt.size = 1,
  ...
) {
  features.use <- features.use %||% rownames(x = object)
  data.plot <- FetchData(
    object = object,
    vars.fetch = features.use,
    cells.use = c(cell1, cell2)
  )
  data.plot <- as.data.frame(x = t(x = data.plot))
  p <- SingleCorPlot(
    data.plot = data.plot,
    cols.use = cols.use,
    pt.size = pt.size,
    smooth = TRUE,
    ...
  )
  return(p)
}

#' View variable features
#'
#' @inheritParams FeatureScatter
#' @param cols.use Colors to specify non-variable/variable status
#' @param assay.use Assay to pull variable features from
#'
#' @importFrom ggplot2 labs scale_color_manual
#' @export
#'
#' @aliases VariableGenePlot
#'
#' @examples
#' VariableFeaturePlot(object = pbmc_small)
#'
VariableFeaturePlot <- function(
  object,
  pt.size = 1,
  cols.use = c('black', 'red'),
  assay.use = NULL
) {
  if (length(x = cols.use) != 2) {
    stop("'cols.use' must be of length 2")
  }
  hvf.info <- HVFInfo(object = object, assay.use = assay.use)[, c('mean', 'dispersion')]
  var.features <- VariableFeatures(object = object, assay.use = assay.use)
  var.status <- ifelse(
    test = rownames(x = hvf.info) %in% var.features,
    yes = 'yes',
    no = 'no'
  )
  p <- SingleCorPlot(
    data.plot = hvf.info,
    col.by = var.status,
    pt.size = pt.size
  )
  p <- p +
    labs(title = NULL, x = 'Average Expression', y = 'Dispersion') +
    scale_color_manual(
      labels = paste(c('Non-variable', 'Variable'), 'count:', table(var.status)),
      values = cols.use
    )
  return(p)
}

#' Visualize Dimensional Reduction genes
#'
#' Visualize top genes associated with reduction components
#'
#' @param object Seurat object
#' @param reduction.use Reduction technique to visualize results for
#' @param dims.use Number of dimensions to display
#' @param num.features Number of genes to display
#' @param pt.color Color of points to use
#' @param projected Use reduction values for full dataset (i.e. projected dimensional reduction values)
#' @param num.col Number of columns to display
#' @param do.balanced Return an equal number of genes with + and - scores. If FALSE (default), returns
#' the top genes ranked by the scores absolute values
#' @param combine.plots Combine plots into a single gg object; note that if TRUE; themeing will not work when plotting multiple features
#' @param ... Ignored
#'
#' @return A ggplot object
#'
#' @importFrom ggplot2 ggplot aes_string geom_point labs
#' @importFrom cowplot theme_cowplot
#' @export
#'
#' @examples
#' VizDimReduction(object = pbmc_small)
#'
VizDimReduction <- function(
  object,
  reduction.use = "pca",
  dims.use = 1:5,
  num.features = 30,
  pt.color = 'blue',
  projected = FALSE,
  num.col = NULL,
  do.balanced = FALSE,
  combine.plots = TRUE,
  ...
) {
  if (is.null(x = num.col)) {
    num.col <- 2
    if (length(x = dims.use) == 1) {
      num.col <- 1
    }
    if (length(x = dims.use) > 6) {
      num.col <- 3
    }
    if (length(x = dims.use) > 9) {
      num.col <- 4
    }
  }
  loadings <- Loadings(object = object[[reduction.use]], projected = projected)
  features.use <- lapply(
    X = dims.use,
    FUN = TopFeatures,
    object = object[[reduction.use]],
    num.features = num.features,
    projected = projected,
    do.balanced = do.balanced
  )
  features.use <- lapply(
    X = features.use,
    FUN = unlist,
    use.names = FALSE
  )
  loadings <- loadings[unlist(x = features.use), dims.use, drop = FALSE]
  names(x = features.use) <- colnames(x = loadings) <- as.character(x = dims.use)
  plots <- lapply(
    X = as.character(x = dims.use),
    FUN = function(i) {
      data.plot <- as.data.frame(x = loadings[features.use[[i]], i, drop = FALSE])
      colnames(x = data.plot) <- paste0(Key(object = object[[reduction.use]]), i)
      data.plot$feature <- factor(x = rownames(x = data.plot), levels = rownames(x = data.plot))
      plot <- ggplot(
        data = data.plot,
        mapping = aes_string(x = colnames(x = data.plot)[1], y = 'feature')
      ) +
        geom_point(col = pt.color) +
        labs(y = NULL) + theme_cowplot()
      return(plot)
    }
  )
  if (combine.plots) {
    plots <- CombinePlots(plot.list = plots, num.col = num.col, legend.position = NULL)
  }
  return(plots)
}

#' Quickly Pick Relevant Dimensions
#'
#' Plots the standard deviations (or approximate singular values if running PCAFast)
#' of the principle components for easy identification of an elbow in the graph.
#' This elbow often corresponds well with the significant dims and is much faster to run than
#' Jackstraw
#'
#' @param object Seurat object
#' @param reduction.use Reduction technique to plot standard deviation for
#' @param dims.plot Number of dimensions to plot standard deviation for
#'
#' @return A ggplot object
#'
#' @importFrom ggplot2 ggplot aes_string geom_point labs element_line
#' @importFrom cowplot theme_cowplot
#' @export
#'
#' @examples
#' ElbowPlot(object = pbmc_small)
#'
ElbowPlot <- function(
  object,
  reduction.use = "pca",
  dims.plot = 20
) {
  data.use <- Stdev(object = object, reduction.use = reduction.use)
  if (length(data.use) == 0) {
    stop(paste("No standard deviation info stored for", reduction.use))
  }
  if (dims.plot > length(x = data.use)) {
    warning("The object only has information for ", length(x = data.use), " reductions")
    dims.plot <- length(x = data.use)
  }
  stdev <- 'Standard Deviation'
  p <- ggplot(data = data.frame(dims = 1:dims.plot, stdev = data.use[1:dims.plot])) +
    geom_point(mapping = aes_string(x = 'dims', y = 'stdev')) +
    labs(x = Key(object = object[[reduction.use]]), y = stdev) +
    theme_cowplot()
  return(p)
}

#' JackStraw Plot
#'
#' Plots the results of the JackStraw analysis for PCA significance. For each
#' PC, plots a QQ-plot comparing the distribution of p-values for all genes
#' across each PC, compared with a uniform distribution. Also determines a
#' p-value for the overall significance of each PC (see Details).
#'
#' Significant PCs should show a p-value distribution (black curve) that is
#' strongly skewed to the left compared to the null distribution (dashed line)
#' The p-value for each PC is based on a proportion test comparing the number
#' of genes with a p-value below a particular threshold (score.thresh), compared with the
#' proportion of genes expected under a uniform distribution of p-values.
#'
#' @param object Seurat object
#' @param reduction.use reduction to pull jackstraw info from
#' @param dims Dims to plot
#' @param plot.x.lim X-axis maximum on each QQ plot.
#' @param plot.y.lim Y-axis maximum on each QQ plot.
#'
#' @return A ggplot object
#'
#' @author Thanks to Omri Wurtzel for integrating with ggplot
#'
#' @importFrom SeuratObject JS
#' @importFrom ggplot2 ggplot aes_string stat_qq labs xlim ylim
#' coord_flip geom_abline guides guide_legend
#' @importFrom cowplot theme_cowplot
#'
#' @export
#'
#' @examples
#' JackStrawPlot(object = pbmc_small)
#'
JackStrawPlot <- function(
  object,
  reduction.use = "pca",
  dims = 1:5,
  plot.x.lim = 0.1,
  plot.y.lim = 0.3
) {
  pAll <- JS(object = object[[reduction.use]], slot = 'empirical')
  if (max(dims) > ncol(x = pAll)) {
    stop("Max dimension is ", ncol(pAll), ".")
  }
  pAll <- pAll[, dims, drop = FALSE]
  pAll <- as.data.frame(x = pAll)
  data.plot <- Melt(x = pAll)
  colnames(x = data.plot) <- c("Contig", "PC", "Value")
  score.df <- JS(object = object[[reduction.use]], slot = 'overall')
  if (nrow(x = score.df) < max(dims)) {
    stop("Jackstraw procedure not scored for all the provided dims. Please run ScoreJackStraw.")
  }
  score.df <- score.df[dims, ]
  if (nrow(score.df) == 0) {
    stop(paste0("JackStraw hasn't been scored. Please run ScoreJackStraw before plotting."))
  }
  data.plot$PC.Score <- rep(
    x = paste0("PC ", score.df[ ,"PC"], ": ", sprintf("%1.3g", score.df[ ,"Score"])),
    each = length(x = unique(x = data.plot$Contig))
  )
  data.plot$PC.Score <- factor(
    x = data.plot$PC.Score,
    levels = paste0("PC ", score.df[, "PC"], ": ", sprintf("%1.3g", score.df[, "Score"]))
  )
  gp <- ggplot(data = data.plot, mapping = aes_string(sample = 'Value', color = 'PC.Score')) +
    stat_qq(distribution = qunif) +
    labs(x = "Theoretical [runif(1000)]", y = "Empirical") +
    xlim(0, plot.y.lim) +
    ylim(0, plot.x.lim) +
    coord_flip() +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", na.rm = TRUE) +
    guides(color = guide_legend(title = "PC: p-value")) +
    theme_cowplot()
  return(gp)
}

#' Feature Locator
#'
#' Select points on a scatterplot and get information about them
#'
#' @param plot A ggplot2 plot
#' @param data.plot The oridinal data that went into the ggplot2 plot
#' @param ... Extra parameters, such as dark.theme, recolor, or smooth for using a dark theme,
#' recoloring based on selected cells, or using a smooth scatterplot, respectively
#'
#' @return The names of the points selected
#'
#' @seealso \code{locator}
#' @seealso \code{ggplot2::ggplot_build}
#' @export
#'
#' @examples
#' \dontrun{
#' df <- data.frame(x = rnorm(n = 100, mean = 20, sd = 2), y = rbinom(n = 100, size = 100, prob = 0.2))
#' p <- ggplot(data = df, mapping = aes(x = x, y = y)) + geom_point(mapping = aes(color = 'red'))
#' FeatureLocator(plot = p, data.plot = df)
#' }
#'
FeatureLocator <- function(plot, data.plot, ...) {
  points.located <- PointLocator(plot = plot, ...)
  #   The rownames for points.located correspond to the row indecies
  #   of data.plot thanks to the way the ggplot object was made
  selected <- data.plot[as.numeric(x = rownames(x = points.located)), ]
  return(rownames(x = selected))
}

#' Hover Locator
#'
#' Get quick information from a scatterplot by hovering over points
#'
#' @param plot A ggplot2 plot
#' @param data.plot The oridinal data that went into the ggplot2 plot
#' @param features.info An optional dataframe or matrix of extra information to be displayed on hover
#' @param dark.theme Plot using a dark theme?
#' @param ... Extra parameters to be passed to plotly::layout
#'
#' @seealso \code{plotly::layout}
#' @seealso \code{ggplot2::ggplot_build}
#' @export
#'
#' @examples
#' \dontrun{
#' df <- data.frame(x = rnorm(n = 100, mean = 20, sd = 2), y = rbinom(n = 100, size = 100, prob = 0.2))
#' p <- ggplot(data = df, mapping = aes(x = x, y = y)) + geom_point(mapping = aes(color = 'red'))
#' HoverLocator(plot = p, data.plot = df)
#' }
#'
HoverLocator <- function(
  plot,
  data.plot,
  features.info = NULL,
  dark.theme = FALSE,
  ...
) {
  #   Use GGpointToBase because we already have ggplot objects
  #   with colors (which are annoying in plotly)
  plot.build <- GGpointToBase(plot = plot, do.plot = FALSE)
  rownames(x = plot.build) <- rownames(data.plot)
  #   Reset the names to 'x' and 'y'
  names(x = plot.build) <- c(
    'x',
    'y',
    names(x = plot.build)[3:length(x = plot.build)]
  )
  #   Add the names we're looking for (eg. cell name, gene name)
  if (is.null(x = features.info)) {
    plot.build$feature <- rownames(x = data.plot)
  } else {
    info <- apply(
      X = features.info,
      MARGIN = 1,
      FUN = function(x, names) {
        return(paste0(names, ': ', x, collapse = '<br>'))
      },
      names = colnames(x = features.info)
    )
    data.info <- data.frame(
      feature = paste(rownames(x = features.info), info, sep = '<br>'),
      row.names = rownames(x = features.info)
    )
    plot.build <- merge(x = plot.build, y = data.info, by = 0)
  }
  #   Set up axis labels here
  #   Also, a bunch of stuff to get axis lines done properly
  xaxis <- list(
    title = names(x = data.plot)[1],
    showgrid = FALSE,
    zeroline = FALSE,
    showline = TRUE
  )
  yaxis <- list(
    title = names(x = data.plot)[2],
    showgrid = FALSE,
    zeroline = FALSE,
    showline = TRUE
  )
  #   Check for dark theme
  if (dark.theme) {
    title <- list(color = 'white')
    xaxis <- c(xaxis, color = 'white')
    yaxis <- c(yaxis, color = 'white')
    plotbg <- 'black'
  } else {
    title = list(color = 'black')
    plotbg = 'white'
  }
  #   Start plotly and pipe it into layout for axis modifications
  #   The `~' means pull from the data passed (this is why we reset the names)
  #   Use I() to get plotly to accept the colors from the data as is
  #   Set hoverinfo to 'text' to override the default hover information
  #   rather than append to it
  plotly::plot_ly(
    data = plot.build,
    x = ~x,
    y = ~y,
    type = 'scatter',
    mode = 'markers',
    color = ~I(color),
    hoverinfo = 'text',
    text = ~feature
  ) %>% plotly::layout(
    xaxis = xaxis,
    yaxis = yaxis,
    titlefont = title,
    paper_bgcolor = plotbg,
    plot_bgcolor = plotbg,
    ...
  )
}

#' Create a custom color palette
#'
#' Creates a custom color palette based on low, middle, and high color values
#'
#' @param low low color
#' @param high high color
#' @param mid middle color. Optional.
#' @param k number of steps (colors levels) to include between low and high values
#'
#' @return A color palette for plotting
#'
#' @importFrom grDevices col2rgb rgb
#' @export
#'
#' @rdname CustomPalette
#' @examples
#' myPalette <- CustomPalette()
#' myPalette
#'
CustomPalette <- function(
  low = "white",
  high = "red",
  mid = NULL,
  k = 50
) {
  low <- col2rgb(col = low) / 255
  high <- col2rgb(col = high) / 255
  if (is.null(x = mid)) {
    r <- seq(from = low[1], to = high[1], len = k)
    g <- seq(from = low[2], to = high[2], len = k)
    b <- seq(from = low[3], to = high[3], len = k)
  } else {
    k2 <- round(x = k / 2)
    mid <- col2rgb(col = mid) / 255
    r <- c(
      seq(from = low[1], to = mid[1], len = k2),
      seq(from = mid[1], to = high[1], len = k2)
    )
    g <- c(
      seq(from = low[2], to = mid[2], len = k2),
      seq(from = mid[2], to = high[2],len = k2)
    )
    b <- c(
      seq(from = low[3], to = mid[3], len = k2),
      seq(from = mid[3], to = high[3], len = k2)
    )
  }
  return(rgb(red = r, green = g, blue = b))
}

#' @inheritParams CustomPalette
#'
#' @export
#'
#' @rdname CustomPalette
#' @aliases BlackAndWhite
#'
#' @examples
#' df <- data.frame(x = rnorm(n = 100, mean = 20, sd = 2), y = rbinom(n = 100, size = 100, prob = 0.2))
#' plot(df, col = BlackAndWhite())
#'
BlackAndWhite <- function(mid = NULL, k = 50) {
  return(CustomPalette(low = "white", high = "black", mid = mid, k = k))
}

#' @inheritParams CustomPalette
#'
#' @export
#'
#' @rdname CustomPalette
#' @aliases PurpleAndYellow
#'
#' @examples
#' df <- data.frame(x = rnorm(n = 100, mean = 20, sd = 2), y = rbinom(n = 100, size = 100, prob = 0.2))
#' plot(df, col = PurpleAndYellow())
#'
PurpleAndYellow <- function(k = 50) {
  return(CustomPalette(low = "magenta", high = "yellow", mid = "black", k = k))
}

#' @inheritParams CustomPalette
#'
#' @export
#'
#' @rdname CustomPalette
#' @aliases BlueAndRed
#'
#' @examples
#' df <- data.frame(x = rnorm(n = 100, mean = 20, sd = 2), y = rbinom(n = 100, size = 100, prob = 0.2))
#' plot(df, col = BlueAndRed())
#'
BlueAndRed <- function(k = 50) {
  return(CustomPalette(low = "#313695" , high = "#A50026", mid = "#FFFFBF", k = k))
}

# NEED TO TREAT PNG PACKAGE PROPERLY
#' Augments ggplot2 scatterplot with a PNG image.
#'
#' Used in to creating vector friendly plots. Exported as it may be useful to others more broadly
#'
#' @param plot1 ggplot2 scatterplot. Typically will have only labeled axes and no points
#' @param imgFile location of a PNG file that contains the points to overlay onto the scatterplot.
#'
#' @return ggplot2 scatterplot that includes the original axes but also the PNG file
#'
#' @importFrom png readPNG
#' @importFrom ggplot2 annotation_raster
#'
#' @export
#'
#' @examples
#' \dontrun{
#' data("pbmc_small")
#' p <- PCAPlot(pbmc_small, do.return = TRUE)
#' ggsave(filename = 'pcaplot.png', plot = p, device = png)
#' pmod <- AugmentPlot(plot1 = p, imgFile = 'pcaplot.png')
#' pmod
#' }
AugmentPlot <- function(plot1, imgFile) {
  range.values <- c(
    ggplot_build(plot = plot1)$layout$panel_params[[1]]$x.range,
    ggplot_build(plot = plot1)$layout$panel_params[[1]]$y.range
  )
  img <- readPNG(source = imgFile)
  p1mod <- plot1 + annotation_raster(
    img,
    xmin = range.values[1],
    xmax = range.values[2],
    ymin = range.values[3],
    ymax = range.values[4]
  )
  return(p1mod)
}

#' Seurat Themes
#'
#' Various themes to be applied to ggplot2-based plots
#' \describe{
#'   \item{\code{SeuratTheme}}{The curated Seurat theme, consists of ...}
#'   \item{\code{DarkTheme}}{A dark theme, axes and text turn to white, the background becomes black}
#'   \item{\code{NoAxes}}{Removes axis lines, text, and ticks}
#'   \item{\code{NoLegend}}{Removes the legend}
#'   \item{\code{NoGrid}}{Removes grid lines}
#'   \item{\code{SeuratAxes}}{Set Seurat-style axes}
#'   \item{\code{BarePlot}}{Remove all extraneous features}
#'   \item{\code{RotatedAxis}}{Rotate X axis text 45 degrees}
#'   \item{\code{BoldTitle}}{Enlarges and emphasizes the title}
#' }
#'
#' @param ... Extra parameters to be passed to \code{theme}
#'
#' @return A ggplot2 theme object
#'
#' @export
#'
#' @rdname SeuratThemes
#' @seealso \code{\link{ggplot2::theme}}
#' @aliases SeuratTheme
#'
SeuratTheme <- function() {
  return(DarkTheme() + NoLegend() + NoGrid() + SeuratAxes())
}

#' @inheritParams SeuratThemes
#'
#' @importFrom ggplot2 theme element_rect element_text element_line margin
#' @export
#'
#' @rdname SeuratThemes
#' @aliases DarkTheme
#'
#' @examples
#' # Generate a plot with a dark theme
#' df <- data.frame(x = rnorm(n = 100, mean = 20, sd = 2), y = rbinom(n = 100, size = 100, prob = 0.2))
#' p <- ggplot(data = df, mapping = aes(x = x, y = y)) + geom_point(mapping = aes(color = 'red'))
#' p + DarkTheme(legend.position = 'none')
#'
DarkTheme <- function(...) {
  #   Some constants for easier changing in the future
  black.background <- element_rect(fill = 'black')
  black.background.no.border <- element_rect(fill = 'black', size = 0)
  font.margin <- 4
  white.text <- element_text(
    colour = 'white',
    margin = margin(
      t = font.margin,
      r = font.margin,
      b = font.margin,
      l = font.margin
    )
  )
  white.line <- element_line(colour = 'white', size = 1)
  no.line <- element_line(size = 0)
  #   Create the dark theme
  dark.theme <- theme(
    #   Set background colors
    plot.background = black.background,
    panel.background = black.background,
    legend.background = black.background,
    legend.box.background = black.background.no.border,
    legend.key = black.background.no.border,
    #   Set text colors
    plot.title = white.text,
    plot.subtitle = white.text,
    axis.title = white.text,
    axis.text = white.text,
    legend.title = white.text,
    legend.text = white.text,
    #   Set line colors
    axis.line.x = white.line,
    axis.line.y = white.line,
    panel.grid = no.line,
    panel.grid.minor = no.line,
    #   Validate the theme
    validate = TRUE,
    #   Extra parameters
    ...
  )
  return(dark.theme)
}

#' @inheritParams SeuratThemes
#'
#' @importFrom ggplot2 theme element_blank
#' @export
#'
#' @rdname SeuratThemes
#' @aliases NoAxes
#'
#' @examples
#' # Generate a plot with no axes
#' df <- data.frame(x = rnorm(n = 100, mean = 20, sd = 2), y = rbinom(n = 100, size = 100, prob = 0.2))
#' p <- ggplot(data = df, mapping = aes(x = x, y = y)) + geom_point(mapping = aes(color = 'red'))
#' p + NoAxes()
#'
NoAxes <- function(...) {
  blank <- element_blank()
  no.axes.theme <- theme(
    # Remove the axis elements
    axis.line.x = blank,
    axis.line.y = blank,
    axis.text.x = blank,
    axis.text.y = blank,
    axis.ticks.x = blank,
    axis.ticks.y = blank,
    axis.title.x = blank,
    axis.title.y = blank,
    # Validate the theme
    validate = TRUE,
    ...
  )
  return(no.axes.theme)
}

#' @inheritParams SeuratThemes
#'
#' @importFrom ggplot2 theme
#' @export
#'
#' @rdname SeuratThemes
#' @aliases NoLegend
#'
#' @examples
#' # Generate a plot with no legend
#' df <- data.frame(x = rnorm(n = 100, mean = 20, sd = 2), y = rbinom(n = 100, size = 100, prob = 0.2))
#' p <- ggplot(data = df, mapping = aes(x = x, y = y)) + geom_point(mapping = aes(color = 'red'))
#' p + NoLegend()
#'
NoLegend <- function(...) {
  no.legend.theme <- theme(
    # Remove the legend
    legend.position = 'none',
    # Validate the theme
    validate = TRUE,
    ...
  )
  return(no.legend.theme)
}

#' @inheritParams SeuratThemes
#'
#' @importFrom ggplot2 theme element_blank
#' @export
#'
#' @rdname SeuratThemes
#' @aliases NoGrid
#'
#' @examples
#' # Generate a plot with no grid lines
#' df <- data.frame(x = rnorm(n = 100, mean = 20, sd = 2), y = rbinom(n = 100, size = 100, prob = 0.2))
#' p <- ggplot(data = df, mapping = aes(x = x, y = y)) + geom_point(mapping = aes(color = 'red'))
#' p + NoGrid()
#'
NoGrid <- function(...) {
  no.grid.theme <- theme(
    # Set grid lines to blank
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    # Validate the theme
    validate = TRUE,
    ...
  )
  return(no.grid.theme)
}

#' @inheritParams SeuratThemes
#'
#' @importFrom ggplot2 theme element_text
#' @export
#'
#' @rdname SeuratThemes
#' @aliases SeuratAxes
#'
SeuratAxes <- function(...) {
  axes.theme <- theme(
    # Set axis things
    axis.title = element_text(face = 'bold', color = '#990000', size = 16),
    axis.text = element_text(vjust = 0.5, size = 12),
    # Validate the theme
    validate = TRUE,
    ...
  )
  return(axes.theme)
}

#' @export
#'
#' @rdname SeuratThemes
#' @aliases BarePlot
#'
BarePlot <- function() {
  return(NoLegend() + NoAxes() + NoGrid())
}

#' @inheritParams SeuratThemes
#'
#' @importFrom ggplot2 theme element_text
#' @export
#'
#' @rdname SeuratThemes
#' @aliases RotatedAxis
#'
RotatedAxis <- function(...) {
  rotated.theme <- theme(
    # Rotate X axis text
    axis.text.x = element_text(angle = 45, hjust = 1),
    # Validate the theme
    validate = TRUE,
    ...
  )
  return(rotated.theme)
}

#' @inheritParams SeuratThemes
#'
#' @importFrom ggplot2 theme element_text
#' @export
#'
#' @rdname SeuratThemes
#' @aliases BoldTitle
#'
BoldTitle <- function(...) {
  bold.theme <- theme(
    # Make the title bold
    plot.title = element_text(size = 20, face = 'bold'),
    # Validate the theme
    validate = TRUE,
    ...
  )
  return(bold.theme)
}

#' @inheritParams SeuratThemes
#'
#' @importFrom ggplot2 theme element_rect
#' @export
#'
#' @rdname SeuratThemes
#' @aliases WhiteBackground
#'
WhiteBackground <- function(...) {
  white.rect = element_rect(fill = 'white')
  white.theme <- theme(
    # Make the plot, panel, and legend key backgrounds white
    plot.background = white.rect,
    panel.background = white.rect,
    legend.key = white.rect,
    # Validate the theme
    validate = TRUE,
    ...
  )
  return(white.theme)
}

# Internal

# Create a scatterplot with data from a ggplot2 scatterplot
#
# @param plot.data The original ggplot2 scatterplot data
# This is taken from ggplot2::ggplot_build
# @param dark.theme Plot using a dark theme
# @param smooth Use a smooth scatterplot instead of a standard scatterplot
# @param ... Extra parameters passed to graphics::plot or graphics::smoothScatter
#
#' @importFrom graphics axis plot smoothScatter
#
PlotBuild <- function(plot.data, dark.theme = FALSE, smooth = FALSE, ...) {
  #   Do we use a smooth scatterplot?
  #   Take advantage of functions as first class objects
  #   to dynamically choose normal vs smooth scatterplot
  myplot <- ifelse(test = smooth, yes = smoothScatter, no = plot)
  if (dark.theme) {
    par(bg = 'black')
    axes = FALSE
    col.lab = 'white'
  } else {
    axes = 'TRUE'
    col.lab = 'black'
  }
  myplot(
    plot.data[, c(1, 2)],
    col = plot.data$color,
    pch = plot.data$pch,
    cex = vapply(
      X = plot.data$cex,
      FUN = function(x) {
        return(max(x / 2, 0.5))
      },
      FUN.VALUE = numeric(1)
    ),
    axes = axes,
    col.lab = col.lab,
    col.main = col.lab,
    ...
  )
  if (dark.theme) {
    axis(
      side = 1,
      at = NULL,
      labels = TRUE,
      col.axis = col.lab,
      col = col.lab
    )
    axis(
      side = 2,
      at = NULL,
      labels = TRUE,
      col.axis = col.lab,
      col = col.lab
    )
  }
}

# Convert a ggplot2 scatterplot to base R graphics
#
# @param plot A ggplot2 scatterplot
# @param do.plot Create the plot with base R graphics
# @param ... Extra parameters passed to PlotBuild
#
# @return A dataframe with the data that created the ggplot2 scatterplot
#
#' @importFrom ggplot2 ggplot_build
#
GGpointToBase <- function(plot, do.plot = TRUE, ...) {
  plot.build <- ggplot_build(plot = plot)
  build.data <- plot.build$data[[1]]
  plot.data <- build.data[, c('x', 'y', 'colour', 'shape', 'size')]
  names(x = plot.data) <- c(
    plot.build$plot$labels$x,
    plot.build$plot$labels$y,
    'color',
    'pch',
    'cex'
  )
  if (do.plot) {
    PlotBuild(plot.data = plot.data, ...)
  }
  return(plot.data)
}

# Locate points on a plot and return them
#
# @param plot A ggplot2 plot
# @param recolor Do we recolor the plot to highlight selected points?
# @param dark.theme Plot using a dark theme
# @param ... Exptra parameters to PlotBuild
#
# @return A dataframe of x and y coordinates for points selected
#
#' @importFrom graphics locator
#' @importFrom SDMTools pnt.in.poly
#
PointLocator <- function(plot, recolor = TRUE, dark.theme = FALSE, ...) {
  #   Convert the ggplot object to a data.frame
  plot.data <- GGpointToBase(plot = plot, dark.theme = dark.theme, ...)
  npoints <- nrow(x = plot.data)
  cat("Click around the cluster of points you wish to select\n")
  cat("ie. select the vertecies of a shape around the cluster you\n")
  cat("are interested in. Press <Esc> when finished (right click for R-terminal users)\n\n")
  polygon <- locator(n = npoints, type = 'l')
  polygon <- data.frame(polygon)
  #   pnt.in.poly returns a data.frame of points
  points.all <- pnt.in.poly(
    pnts = plot.data[, c(1, 2)],
    poly.pnts = polygon
  )
  #   Find the located points
  points.located <- points.all[which(x = points.all$pip == 1), ]
  #   If we're recoloring, do the recolor
  if (recolor) {
    no <- ifelse(test = dark.theme, yes = 'white', no = 'black')
    points.all$color <- ifelse(test = points.all$pip == 1, yes = 'red', no = no)
    plot.data$color <- points.all$color
    PlotBuild(plot.data = plot.data, dark.theme = dark.theme, ...)
  }
  return(points.located[, c(1, 2)])
}

# TODO Integrate BlendPlot into SingleDimPlot
# Blend two feature plots together
#
# @param data.use The data regarding the feature
# @param features.plot The features to plot
# @param data.plot The data to be plotted
# @param pt.size Size of each point
# @param pch.use Shape of each point
# @param cols.use Colors to plot
# @param dim.codes Codes for the dimensions to plot in
# @param min.cutoff Minimum cutoff for data
# @param max.cutoff Maximum cutoff for data
# @param coord.fixed Use a fixed scale coordinate system (for spatial coordinates)
# @param no.axes Remove axes from plot
# @param no.legend Remove legend from plot
# @param dark.theme Plot in dark theme
#
#' @import RColorBrewer
#' @importFrom grDevices colors
#' @importFrom utils globalVariables
#
# @return A blended ggplot2 scatterplot
#
globalVariables(names = c('x', 'y'), package = 'Seurat', add = TRUE)
BlendPlot <- function(
  data.use,
  features.plot,
  data.plot,
  pt.size,
  pch.use,
  cols.use,
  dim.codes,
  min.cutoff,
  max.cutoff,
  coord.fixed,
  no.axes,
  no.legend,
  dark.theme
) {
  num.cols <- length(x = cols.use)
  #   Create a vector of colors that weren't provided
  cols.not.provided <- colors(distinct = TRUE)
  cols.not.provided <- cols.not.provided[!(grepl(
    pattern = paste(cols.use, collapse = '|'),
    x = cols.not.provided,
    ignore.case = TRUE
  ))]
  if (num.cols > 4) {
    #   If provided more than four colors, take only the first four
    cols.use <- cols.use[c(1:4)]
  } else if ((num.cols == 2) || (num.cols == 3)) {
    #   If two or three colors, use the last two as high values for blending
    #   and add to our vector of colors
    blend <- BlendColors(cols.use[c(num.cols - 1, num.cols)])
    cols.use <- c(cols.use, blend)
    if (num.cols == 2) {
      #   If two colors, provided,
      #   we still need a low color
      cols.use <- c(sample(x = cols.not.provided, size = 1), cols.use)
    }
  } else if ((num.cols == 1)) {
    #   If only one color provided
    if (cols.use %in% rownames(x = brewer.pal.info)) {
      #   Was it a palette from RColorBrewer? If so, create
      #   our colors based on the palette
      palette <- brewer.pal(n = 3, name = cols.use)
      cols.use <- c(palette, BlendColors(palette[c(2, 3)]))
    } else {
      #   If not, randomly create our colors
      cols.high <- sample(x = cols.not.provided, size = 2, replace = FALSE)
      cols.use <- c(cols.use, cols.high, BlendColors(cols.high))
    }
  } else if (num.cols <= 0) {
    cols.use <- c('yellow','red', 'blue', BlendColors('red', 'blue'))
  }
  names(x = cols.use) <- c('low', 'high1', 'high2', 'highboth')
  length.check <- vapply(
    X = list(features.plot, min.cutoff, max.cutoff),
    FUN = function(x) {
      return(length(x = x) != 2)
    },
    FUN.VALUE = logical(length = 1)
  )
  if (any(length.check)) {
    stop("An overlayed FeaturePlot only works with two features and requires two minimum and maximum cutoffs")
  }
  #   Check for quantiles
  min.cutoff <- c(
    SetQuantile(cutoff = min.cutoff[1], data = data.gene[features.plot[1], ]),
    SetQuantile(cutoff = min.cutoff[2], data = data.gene[features.plot[2], ])
  )
  max.cutoff <- c(
    SetQuantile(cutoff = max.cutoff[1], data = data.gene[features.plot[1], ]),
    SetQuantile(cutoff = max.cutoff[2], data = data.gene[features.plot[2], ])
  )
  data.gene <- stats::na.omit(object = data.frame(data.use[features.plot, ]))
  cell.names <- colnames(x = data.gene)
  #   Minimum and maximum masking
  data.gene <- matrix(
    data = vapply(
      X = data.gene,
      FUN = function(x) ifelse(test = x < min.cutoff, yes = min.cutoff, no = x),
      FUN.VALUE = c(1, 1)
    ),
    nrow = 2
  )
  data.gene <- matrix(
    data = vapply(
      X = as.data.frame(x = data.gene),
      FUN = function(x) ifelse(test = x > max.cutoff, yes = max.cutoff, no = x),
      FUN.VALUE = c(1, 1)
    ),
    nrow = 2
  )
  data.gene <- as.data.frame(x = data.gene)
  rownames(x = data.gene) <- features.plot
  colnames(x = data.gene) <- cell.names
  #   Stuff for break points
  if(all(data.gene ==0)) {
    data.cut <- 0
  } else {
    #   Cut the expression of both features
    cuts <- apply(
      X = data.gene,
      MARGIN = 1,
      FUN = cut,
      breaks = 2,
      labels = FALSE
    )
    cuts.dim <- dim(x = cuts)
    if (cuts.dim[1] > cuts.dim[2]){
      cuts <- t(x = cuts)
    }
    #   Apply colors dependent on if the cell expresses
    #   none, one, or both features
    data.cut = apply(
      X = cuts,
      MARGIN = 2,
      FUN = function(x) {
        return(if ((x[1] == 1) && (x[2] == 2)) { # Expression in 2
          'high2'
        } else if ((x[1] == 2) && (x[2] == 1)) { # Expression in 1
          'high1'
        } else if ((x[1] == 2) && (x[2] == 2)) { # Expression in both
          'highboth'
        } else { # Expression in neither
          'low'
        })
      }
    )
    data.cut <- as.factor(x = data.cut)
  }
  data.plot$colors <- data.cut
  #   Start plotting
  legend.names <- c(
    'high1' = paste('High', features.plot[1]),
    'high2' = paste('High', features.plot[2]),
    'highboth' = 'High both'
  )
  title <- paste0(features.plot, collapse = ' x ')
  p <- ggplot(data = data.plot, mapping = aes(x = x, y = y))
  p <- p + geom_point(
    mapping = aes(color = colors),
    size = pt.size,
    shape = pch.use
  )
  p <- p + scale_color_manual(
    values = cols.use,
    limits = c('high1', 'high2', 'highboth'),
    labels = legend.names,
    guide = guide_legend(title = NULL, override.aes = list(size = 2))
  )
  #   Deal with axes and legends
  if (no.axes) {
    p <- p + labs(title = title, x ="", y="") + theme(
      axis.line = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank()
    )
  } else {
    p <- p + labs(title = title, x = dim.codes[1], y = dim.codes[2])
  }
  if (no.legend){
    p <- p + theme(legend.position = 'none')
  }
  if (dark.theme) {
    p <- p + DarkTheme()
  }
  if(coord.fixed) {
    p <- p + coord_fixed()
  }
  return(p)
}

# Blend two or more colors together
#
# @param ... Two or more colors to blend together
# These can be in a vector or standalone
# @param as.rgb Return in RGB form, otherwise return in hexadecimal form
#
# @return The blended color in RGB form (1 x 3 matrix) or hexadecimal form
#
#' @importFrom grDevices rgb col2rgb
#
BlendColors <- function(..., as.rgb = FALSE) {
  #   Assemble the arguments passed into a character vector
  colors <- as.character(x = c(...))
  if (length(x = colors) < 2) {
    stop("Please provide two or more colors to blend")
  }
  #   Check for hexadecimal values for automatic alpha blending
  alpha.value <- 255
  if (sum(sapply(X = colors, FUN = grepl, pattern = '^#')) != 0) {
    hex <- colors[which(x = grepl(pattern = '^#', x = colors))]
    hex.length <- sapply(X = hex, FUN = nchar)
    #   9-character hexadecimal values specify alpha levels
    if (9 %in% hex.length) {
      hex.alpha <- hex[which(x = hex.length == 9)]
      hex.vals <- sapply(X = hex.alpha, FUN = substr, start = 8, stop = 9)
      dec.vals <- sapply(X = hex.vals, FUN = strtoi, base = 16)
      dec.vals <- dec.vals / 255 # Convert to 0:1 scale for calculations
      alpha.value <- dec.vals[1]
      #   Blend alpha levels, going top-down
      for (val in dec.vals[-1]) {
        alpha.value <- alpha.value + (val * (1 - alpha.value))
      }
      alpha.value <- alpha.value * 255 # Convert back to 0:255 scale
    }
  }
  #   Convert to a 3 by `length(colors)` matrix of RGB values
  rgb.vals <- sapply(X = colors, FUN = col2rgb)
  if (nrow(x = rgb.vals) != 3) {
    rgb.vals <- t(x = rgb.vals)
  }
  #   Blend together using the additive method
  #   Basically, resulting colors are the mean of the component colors
  blend <- as.matrix(x = rowMeans(x = rgb.vals))
  dimnames(x = blend) <- list(c('red', 'green', 'blue'), 'blend')
  #   If we're returning RGB values, convert to matrix, just like col2rgb
  #   Otherwise, return as hexadecimal; can be used directly for plotting
  if (!as.rgb) {
    blend <- rgb(t(x = blend), alpha = alpha.value, maxColorValue = 255)
  }
  return(blend)
}

# Find the quantile of a data
#
# Converts a quantile in character form to a number regarding some data
# String form for a quantile is represented as a number prefixed with 'q'
# For example, 10th quantile is 'q10' while 2nd quantile is 'q2'
#
# Will only take a quantile of non-zero data values
#
# @param cutoff The cutoff to turn into a quantile
# @param data The data to turn find the quantile of
#
# @return The numerical representation of the quantile
#
#' @importFrom stats quantile
#
SetQuantile <- function(cutoff, data) {
  if (grepl(pattern = '^q[0-9]{1,2}$', x = as.character(x = cutoff), perl = TRUE)) {
    this.quantile <- as.numeric(x = sub(
      pattern = 'q',
      replacement = '',
      x = as.character(x = cutoff)
    )) / 100
    data <- unlist(x = data)
    data <- data[data > 0]
    cutoff <- quantile(x = data, probs = this.quantile)
  }
  return(as.numeric(x = cutoff))
}

# Set highlight information
#
# @param cells.highlight Cells to highlight
# @param cells.all A character vector of all cell names
# @param sizes.highlight Sizes of cells to highlight
# @param cols.highlight Colors to highlight cells as
# @param col.base Base color to use for unselected cells
# @param pt.size Size of unselected cells
#
# @return A list will cell highlight information
# \describe{
#   \item{plot.order}{An order to plot cells in}
#   \item{highlight}{A vector giving group information for each cell}
#   \item{size}{A vector giving size information for each cell}
#   \item{color}{Colors for highlighting in the order of plot.order}
# }
#
SetHighlight <- function(
  cells.highlight,
  cells.all,
  sizes.highlight,
  cols.highlight,
  col.base = 'black',
  pt.size = 1
) {
  if (is.character(x = cells.highlight)) {
    cells.highlight <- list(cells.highlight)
  } else if (is.data.frame(x = cells.highlight) || !is.list(x = cells.highlight)) {
    cells.highlight <- as.list(x = cells.highlight)
  }
  cells.highlight <- lapply(
    X = cells.highlight,
    FUN = function(cells) {
      cells.return <- if (is.character(x = cells)) {
        cells[cells %in% cells.all]
      } else {
        cells <- as.numeric(x = cells)
        cells <- cells[cells <= length(x = cells.all)]
        cells.all[cells]
      }
      return(cells.return)
    }
  )
  cells.highlight <- Filter(f = length, x = cells.highlight)
  names.highlight <- if (is.null(x = names(x = cells.highlight))) {
    paste0('Group_', 1L:length(x = cells.highlight))
  } else {
    names(x = cells.highlight)
  }
  sizes.highlight <- rep_len(
    x = sizes.highlight,
    length.out = length(x = cells.highlight)
  )
  cols.highlight <- c(
    col.base,
    rep_len(x = cols.highlight, length.out = length(x = cells.highlight))
  )
  size <- rep_len(x = pt.size, length.out = length(x = cells.all))
  highlight <- rep_len(x = NA_character_, length.out = length(x = cells.all))
  if (length(x = cells.highlight) > 0) {
    for (i in 1:length(x = cells.highlight)) {
      cells.check <- cells.highlight[[i]]
      index.check <- match(x = cells.check, cells.all)
      highlight[index.check] <- names.highlight[i]
      size[index.check] <- sizes.highlight[i]
    }
  }
  plot.order <- sort(x = unique(x = highlight), na.last = TRUE)
  plot.order[is.na(x = plot.order)] <- 'Unselected'
  highlight[is.na(x = highlight)] <- 'Unselected'
  highlight <- as.factor(x = highlight)
  return(list(
    plot.order = plot.order,
    highlight = highlight,
    size = size,
    color = cols.highlight
  ))
}

# Make label information for ggplot2-based scatter plots
#
# @param data.plot A three-column data frame (accessed with \code{plot$data})
# The first column should be the X axis, the second the Y, and the third should be grouping information
#
# @return A dataframe with three columns: centers along the X axis, centers along the Y axis, and group information
#
MakeLabels <- function(data.plot) {
  data.labels <- lapply(
    X = unique(x = data.plot[, 3]),
    FUN = function(group) {
      data.use <- data.plot[data.plot[, 3] == group, 1:2]
      return(apply(X = data.use, MARGIN = 2, FUN = median))
    }
  )
  names(x = data.labels) <- as.character(unique(x = data.plot[, 3]))
  data.labels <- as.data.frame(x = t(x = as.data.frame(x = data.labels)))
  data.labels[, colnames(x = data.plot)[3]] <- as.character(unique(x = data.plot[, 3]))
  return(data.labels)
}

# Combine ggplot2-based plots into a single plot
#
# @param plot.list A list of gg objects
# @param num.col Number of columns
# @param legend.position Combine legends into a single legend
# choose from 'right' or 'bottom'; pass 'none' to remove legends, or \code{NULL}
# to leave legends as they are
#
# @return A combined plot
#
#' @importFrom cowplot plot_grid get_legend
#
CombinePlots <- function(plot.list, num.col, legend.position = NULL) {
  plots.combined <- if (length(x = plot.list) > 1) {
    if (!is.null(x = legend.position)) {
      if (legend.position != 'none') {
        legend <- get_legend(plot = plot.list[[1]] + theme(legend.position = legend.position))
      }
      plot.list <- lapply(
        X = plot.list,
        FUN = function(x) {
          return(x + NoLegend())
        }
      )
    }
    plots.combined <- plot_grid(plotlist = plot.list, ncol = num.col, align = 'v')
    if (!is.null(x = legend.position)) {
      plots.combined <- switch(
        EXPR = legend.position,
        'bottom' = plot_grid(
          plots.combined,
          legend,
          ncol = 1,
          rel_heights = c(1, 0.2)
        ),
        'right' = plot_grid(
          plots.combined,
          legend,
          rel_widths = c(3, 0.3)
        ),
        plots.combined
      )
    }
    plots.combined
  } else {
    plot.list[[1]]
  }
  return(plots.combined)
}

# Plot feature expression by identity
#
# Basically combines the codebase for VlnPlot and RidgePlot
#
# @param object Seurat object
# @param features.plot Features to plot (gene expression, metrics, PC scores,
# anything that can be retreived by FetchData)
# @param ident.include Which classes to include in the plot (default is all)
# @param num.col Number of columns if multiple plots are displayed
# @param do.sort Sort identity classes (on the x-axis) by the average
# expression of the attribute being potted
# @param y.max Maximum y axis value
# @param same.y.lims Set all the y-axis limits to the same values
# @param size.x.use X axis title font size
# @param size.y.use Y axis title font size
# @param adjust.use Adjust parameter for geom_violin
# @param point.size.use Point size for geom_violin
# @param cols.use Colors to use for plotting
# @param group.by Group (color) cells in different ways (for example, orig.ident)
# @param y.log plot Y axis on log scale
# @param combine.plots Combine plots using cowplot::plot_grid
# @param plot.counts Use non-normalized counts data for plotting
# @param ... Ignores
#
#
ExIPlot <- function(
  object,
  plot.type = 'violin',
  features.plot,
  ident.include = NULL,
  num.col = NULL,
  do.sort = FALSE,
  y.max = NULL,
  same.y.lims = FALSE,
  adjust.use = 1,
  cols.use = NULL,
  group.by = NULL,
  y.log = FALSE,
  combine.plots = TRUE,
  plot.counts = FALSE,
  ...
) {
  if (is.null(x = num.col)) {
    if (length(x = features.plot) > 9) {
      num.col <- 4
    } else {
      num.col <- min(length(x = features.plot), 3)
    }
  }
  slot.use = "data"
  if (plot.counts == TRUE) {
    slot.use = "counts"
  }
  data.use <- FetchData(object = object, vars.fetch = features.plot,slot = slot.use)
  if (is.null(x = ident.include)) {
    cells.use <- colnames(x = object)
  } else {
    cells.use <- Idents(object = object)[Idents(object = object) == ident.include]
  }
  data.use <- data.use[cells.use, ,drop = FALSE]
  ident.use <- if (is.null(x = group.by)) {
    Idents(object = object)[cells.use]
  } else {
    object[group.by][cells.use, , drop = TRUE]
  }
  feature.names <- colnames(x = data.use)[colnames(x = data.use) %in% rownames(x = object)]
  if (same.y.lims && is.null(x = y.max)) {
    y.max <- max(data.use)
  }
  plots <- lapply(
    X = features.plot,
    FUN = function(x) {
      return(SingleExIPlot(
        feature = x,
        plot.type = plot.type,
        data = data.use[, x, drop = FALSE],
        cell.ident = ident.use,
        do.sort = do.sort, y.max = y.max,
        size.title.use = size.title.use,
        adjust.use = adjust.use,
        cols.use = cols.use,
        feature.names = feature.names,
        y.log = y.log
      ))
    }
  )
  if (combine.plots) {
    plots <- CombinePlots(plot.list = plots, num.col = num.col, legend.position = 'none')
  }
  return(plots)
}

# Plot a single expression by identity on a plot
#
# @param feature Feature to plot
# @param plot.type Make either a 'ridge' or 'violin' plot
# @param data Data to plot
# @param cell.ident Idents to use
# @param do.sort Sort identity classes (on the x-axis) by the average
# expression of the attribute being potted
# @param y.max Maximum Y value to plot
# @param size.x.use X axis title font size
# @param size.y.use Y axis title font size
# @param adjust.use Adjust parameter for geom_violin
# @param point.size.use Point size for geom_violin
# @param cols.use Colors to use for plotting
# @param gene.names
# @param y.log plot Y axis on log scale
# @param ... Ignored
#
# @return A ggplot-based Expression-by-Identity plot
#
# @import ggplot2
#' @importFrom stats rnorm
#' @importFrom utils globalVariables
#' @importFrom ggridges geom_density_ridges theme_ridges
#' @importFrom ggplot2 ggplot aes_string guides guide_legend theme labs geom_violin ylim
#' scale_fill_manual scale_y_discrete scale_y_log10 scale_x_log10 scale_x_continuous element_line
#' @importFrom cowplot theme_cowplot
#'
SingleExIPlot <- function(
  feature,
  plot.type = 'violin',
  data,
  cell.ident,
  do.sort,
  y.max,
  size.title.use,
  adjust.use,
  point.size.use,
  cols.use,
  feature.names,
  y.log,
  ...
) {
  set.seed(seed = 42)
  feature.name <- colnames(x = data)
  feature <- colnames(x = data) <- "feature"
  data$ident <- cell.ident
  if (do.sort) {
    data$ident <- factor(
      x = data$ident,
      levels = names(x = rev(x = sort(x = tapply(
        X = data[, feature],
        INDEX = data$ident,
        FUN = mean
      ))))
    )
  }
  if (y.log) {
    noise <- rnorm(n = length(x = data[, feature])) / 200
    data[, feature] <- data[, feature] + 1
  } else {
    noise <- rnorm(n = length(x = data[, feature])) / 100000
  }
  if (all(data[, feature] == data[, feature][1])) {
    warning(paste0("All cells have the same value of ", feature, "."))
  } else{
    data[, feature] <- data[, feature] + noise
  }
  axis.label <- if (feature %in% feature.names) {
    if (y.log) {
      "Log Expression level"
    } else {
      "Expression level"
    }
  } else {
    ""
  }
  y.max <- y.max %||% max(data[, feature])
  p <- ggplot(data = data, mapping = aes_string(fill = 'ident')) +
    guides(fill = guide_legend(title = NULL)) +
    NoGrid() +
    labs(title = feature.name)
  p <- switch(
    EXPR = plot.type,
    'violin' = {
      p <- p +
        geom_violin(
          scale = 'width',
          adjust = adjust.use,
          trim = TRUE,
          mapping = aes_string(x = 'ident', y = 'feature')
        ) +
        labs(x = 'Identity', y = axis.label)
      p <- p + if (y.log) {
        scale_y_log10()
      } else {
        ylim(min(data[, feature]), y.max)
      }
      p
    },
    'ridge' = {
      p <- p +
        geom_density_ridges(
          scale = 4,
          mapping = aes_string(x = 'feature', y = 'ident')
        ) +
        theme_ridges() +
        labs(x = axis.label, y = 'Identity') + #, title = feature.name) +
        scale_y_discrete(expand = c(0.01, 0)) +   # will generally have to set the `expand` option
        scale_x_continuous(expand = c(0, 0))      # for both axes to remove unneeded padding
      if (y.log) {
        p <- p + scale_x_log10()
      }
      p
    },
    stop("Unknown plot type: ", plot.type)
  )
  if (!is.null(x = cols.use)) {
    p <- p + scale_fill_manual(values = cols.use)
  }
  p <- p + theme_cowplot()
  return(p)
}

# Plot a single dimension
#
# @param data.plot Data to plot
# @param dims.use A two-length numeric vector with dimensions to use
# @param pt.size Adjust point size for plotting
# @param col.by ...
# @param cols.use Vector of colors, each color corresponds to an identity class. By default, ggplot assigns colors.
# @param shape.by If NULL, all points are circles (default). You can specify any
# cell attribute (that can be pulled with FetchData) allowing for both
# different colors and different shapes on cells.
# @param plot.order Specify the order of plotting for the idents. This can be
# useful for crowded plots if points of interest are being buried. Provide
# either a full list of valid idents or a subset to be plotted last (on top).
# @param do.label Whether to label the clusters
# @param label.size Sets size of labels
# @param cells.highlight A list of character or numeric vectors of cells to
# highlight. If only one group of cells desired, can simply
# pass a vector instead of a list. If set, colors selected cells to the color(s)
# in \code{cols.highlight} and other cells black (white if dark.theme = TRUE);
#  will also resize to the size(s) passed to \code{sizes.highlight}
# @param cols.highlight A vector of colors to highlight the cells as; will
# repeat to the length groups in cells.highlight
# @param sizes.highlight Size of highlighted cells; will repeat to the length
# groups in cells.highlight
# @param na.value Color value for NA points when using custom scale.
# @param ... Ignored for now
#
#' @importFrom ggplot2 ggplot aes_string labs geom_text
#' scale_color_brewer scale_color_manual element_rect
#' @importFrom cowplot theme_cowplot
#'
SingleDimPlot <- function(
  data.plot,
  dims.use,
  col.by = NULL,
  cols.use = NULL,
  pt.size = 1,
  shape.by = NULL,
  plot.order = NULL,
  do.label = FALSE,
  label.size = 4,
  cells.highlight = NULL,
  cols.highlight = 'red',
  sizes.highlight = 1,
  na.value = 'grey50',
  ...
) {
  if (length(x = dims.use) != 2) {
    stop("'dims.use' must be a two-length vector")
  }
  if (!is.data.frame(x = data.plot)) {
    data.plot <- as.data.frame(x = data.plot)
  }
  if (is.character(x = dims.use) && !all(dims.use %in% colnames(x = data.plot))) {
    stop("")
  } else if (is.numeric(x = dims.use)) {
    dims.use <- colnames(x = data.plot)[dims.use]
  }
  if (!is.null(x = cells.highlight)) {
    highlight.info <- SetHighlight(
      cells.highlight = cells.highlight,
      cells.all = colnames(x = object),
      sizes.highlight = sizes.highlight,
      cols.highlight = cols.highlight,
      col.base = cols.use[1] %||% 'black',
      pt.size = pt.size
    )
    plot.order <- highlight.info$plot.order
    data.plot$highlight <- highlight.info$highlight
    col.by <- 'highlight'
    pt.size <- highlight.info$size
    cols.use <- highlight.info$color
  }
  if (!is.null(x = plot.order) && !is.null(x = col.by)) {
    plot.order <- rev(x = c(
      plot.order,
      setdiff(x = unique(x = data.plot[, col.by]), y = plot.order)
    ))
    data.plot[, col.by] <- factor(x = data.plot[, col.by], levels = plot.order)
  }
  if (!is.null(x = col.by) && !col.by %in% colnames(x = data.plot)) {
    warning("Cannot find ", col.by, " in plotting data, not coloring plot")
    col.by <- NULL
  }
  if (!is.null(x = shape.by) && !shape.by %in% colnames(x = data.plot)) {
    warning("Cannot find ", shape.by, " in plotting data, not shaping plot")
  }
  p <- ggplot(data = data.plot) +
    geom_point(
      mapping = aes_string(
        x = dims.use[1],
        y = dims.use[2],
        color = col.by,
        shape = shape.by
      ),
      size = pt.size
    ) + labs(color = NULL)
  if (do.label && !is.null(x = col.by)) {
    labels <- MakeLabels(data.plot = p$data[, c(dims.use, col.by)])
    p <- p +
      geom_point(
        data = labels,
        mapping = aes_string(
          x = colnames(x = labels)[1],
          y = colnames(x = labels)[2]
        ),
        size = 0,
        alpha = 0
      ) +
      geom_text(
        data = labels,
        mapping = aes_string(
          x = colnames(x = labels)[1],
          y = colnames(x = labels)[2],
          label = colnames(x = labels)[3]
        ),
        size = label.size
      )
  }
  if (!is.null(x = cols.use)) {
    cols.use <- if (length(x = cols.use) == 1) {
      scale_color_brewer(palette = cols.use)
    } else {
      scale_color_manual(values = cols.use, na.value = na.value)
    }
    p <- p + cols.use
  }
  p <- p + theme_cowplot()
  return(p)
}

# A single correlation plot
#
# @param data.plot A data frame with two columns to be plotted
# @param col.by A vector or factor of values to color the plot by
# @param cols.use An optional vector of colors to use
# @param pt.size Point size for the plot
# @param smooth Make a smoothed scatter plot
# @param legend.title Optional legend title
# @param ... Extra parameters to MASS::kde2d
#
#' @importFrom stats cor
#' @importFrom MASS kde2d
#' @importFrom ggplot2 ggplot geom_point aes_string labs scale_color_brewer
#' scale_color_manual geom_tile guides element_rect
#'
SingleCorPlot <- function(
  data.plot,
  col.by = NULL,
  cols.use = NULL,
  pt.size = 1,
  smooth = FALSE,
  legend.title = NULL,
  na.value = 'grey50',
  ...
) {
  names.plot <- colnames(x = data.plot)
  plot.cor <- round(x = cor(x = data.plot[, 1], y = data.plot[, 2]), digits = 2)
  if (!is.null(x = col.by)) {
    data.plot$colors <- col.by
  }
  p <- ggplot(
    data = data.plot,
    mapping = aes_string(x = names.plot[1], y = names.plot[2])
  ) +
    labs(x = names.plot[1], y = names.plot[2], title = plot.cor, color = legend.title)
  if (!is.null(x = col.by)) {
    p <- p + geom_point(mapping = aes_string(color = 'colors'), size = pt.size)
  } else {
    p <- p + geom_point(size = pt.size)
  }
  if (smooth) {
    density <- kde2d(x = data.plot[, 1], y = data.plot[, 2], ...)
    density <- data.frame(
      x = unlist(x = lapply(X = density$x, FUN = rep.int, times = length(x = density$x))),
      y = rep.int(x = density$y, times = length(x = density$y)),
      z = unlist(x = as.data.frame(x = density$z))
    )
    p <- p + geom_tile(mapping = aes_string(x = 'x', y = 'y', fill = 'z'), data = density) +
      guides(fill = FALSE)
  }
  if (!is.null(x = cols.use)) {
    cols.scale <- if (length(x = cols.use) == 1 &&
                    cols.use %in% rownames(RColorBrewer::brewer.pal.info)) {
        scale_color_brewer(palette = cols.use)
    } else {
      scale_color_manual(values = cols.use, na.value = na.value)
    }
    p <- p + cols.scale
  }
  p <- p + WhiteBackground(panel.border = element_rect(fill = NA, colour = 'black'))
  return(p)
}

# A single heatmap from ggplot2 using geom_raster
#
# @param data.plot A matrix or data frame with data to plot
# @param cell.order ...
# @param feature.order ...
# @param colors A vector of colors to use
# @param disp.min Minimum display value (all values below are clipped)
# @param disp.max Maximum display value (all values above are clipped)
# @param limits A two-length numeric vector with the limits for colors on the plot
# @param group.by A vector to group cells by, should be one grouping identity per cell
# @param ... Ignored
#
#' @importFrom ggplot2 ggplot aes_string geom_raster scale_fill_gradient
#' scale_fill_gradientn theme element_blank labs geom_point guides guide_legend
#
SingleRasterMap <- function(
  data.plot,
  cell.order = NULL,
  feature.order = NULL,
  colors = PurpleAndYellow(),
  disp.min = -2.5,
  disp.max = 2.5,
  limits = NULL,
  group.by = NULL,
  ...
) {
  data.plot <- MinMax(data = data.plot, min = disp.min, max = disp.max)
  data.plot <- Melt(x = t(x = data.plot))
  colnames(x = data.plot) <- c('Feature', 'Cell', 'Expression')
  if (!is.null(x = feature.order)) {
    data.plot$Feature <- factor(x = data.plot$Feature, levels = unique(x = feature.order))
  }
  if (!is.null(x = cell.order)) {
    data.plot$Cell <- factor(x = data.plot$Cell, levels = unique(x = cell.order))
  }
  if (!is.null(x = group.by)) {
    data.plot$Identity <- group.by[data.plot$Cell]
  }
  limits <- limits %||% c(min(data.plot$Expression), max(data.plot$Expression))
  if (length(x = limits) != 2 || !is.numeric(x = limits)) {
    stop("limits' must be a two-length numeric vector")
  }
  p <- ggplot(data = data.plot) +
    geom_raster(mapping = aes_string(x = 'Cell', y = 'Feature', fill = 'Expression')) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
    scale_fill_gradientn(limits = limits, colors = colors) +
    labs(x = NULL, y = NULL, fill = group.by %iff% 'Expression') +
    WhiteBackground()
  if (!is.null(x = group.by)) {
    p <- p + geom_point(
      mapping = aes_string(x = 'Cell', y = 'Feature', color = 'Identity'),
      alpha = 0
    ) +
      guides(color = guide_legend(override.aes = list(alpha = 1)))
  }
  return(p)
}

# A single heatmap from base R using image
#
# @param data.plot matrix of data to plot
# @param cell.order optional vector of cell names to specify order in plot
# @param plot.title Title for plot
#
SingleImageMap <- function(
  data.plot,
  cell.order = NULL,
  plot.title = NULL
) {
  if (!is.null(x = cell.order)) {
    data.plot <- data.plot[cell.order, ]
  }
  par(mar = c(1,1,3,3))
  plot.new()
  image(
    x = as.matrix(x = data.plot),
    axes = FALSE,
    add = TRUE,
    col = PurpleAndYellow()
  )
  axis(
    side = 4,
    at = seq(from = 0, to = 1, length = ncol(x = data.plot)),
    labels = colnames(x = data.plot),
    las = 1,
    tick = FALSE,
    mgp = c(0, -0.5, 0),
    cex.axis = 0.75
  )
  title(main = plot.title)
}

globalVariables(
  names = 'cc',
  package = 'Seurat',
  add = TRUE
)
# Evaluate CCs
#
# Looks at the biweight midcorrelation of the Xth gene across the specified CCs
# for each group in the grouping.var.
#
# @param object A Seurat object
# @param grouping.var Grouping variable specified in alignment procedure
# @param dims.eval dimensions to evalutate the bicor for
# @param gene.num Xth gene to look at bicor for
# @param num.possible.genes Number of possible genes to search when choosing
# genes for the metagene. Set to 2000 by default. Lowering will decrease runtime
# but may result in metagenes constructed on fewer than num.genes genes.
# @param display.progress Show progress bar
#
EvaluateCCs <- function(
  object,
  grouping.var,
  dims.eval,
  gene.num,
  num.possible.genes,
  display.progress
) {
  reduction.type <-  "cca"
  ident.orig <- object@ident
  object <- SetAllIdent(object = object, id = grouping.var)
  levels.split <- names(x = sort(x = table(object@ident), decreasing = T))
  num.groups <- length(levels.split)
  objects <- list()
  for (i in 1:num.groups) {
    objects[[i]] <- SubsetData(object = object, ident.use = levels.split[i])
  }
  object@ident <- ident.orig
  cc.loadings <- list()
  scaled.data <- list()
  cc.embeds <- list()
  for (i in 1:num.groups) {
    cat(paste0("Rescaling group ", i, "\n"), file = stderr())
    objects[[i]] <- ScaleData(
      object = objects[[i]],
      block.size = 5000,
      display.progress = display.progress
    )
    objects[[i]] <- ProjectDim(
      object = objects[[i]],
      reduction.type = reduction.type,
      do.print = FALSE
    )
    cc.loadings[[i]] <- GetGeneLoadings(
      object = objects[[i]],
      reduction.type = reduction.type,
      use.full = TRUE
    )
    cc.embeds[[i]] <- GetCellEmbeddings(
      object = objects[[i]],
      reduction.type = reduction.type
    )
    scaled.data[[i]] <- objects[[i]]@scale.data
  }
  bc.gene <- matrix(ncol = num.groups, nrow = length(dims.eval))
  if (display.progress) {
    cat(paste0("Evaluating dims: ", paste(dims.eval, collapse = " "),  "\n"), file = stderr())
    pb <- txtProgressBar(min = 0, max = length(dims.eval) * (num.groups - 1), style = 3)
    pb.idx <- 0
  }
  for (cc.use in dims.eval) {
    bc.gene.g1 <- c()
    for (g in 2:num.groups) {
      if (display.progress) {
        pb.idx <- pb.idx + 1
        setTxtProgressBar(pb, pb.idx)
      }
      genes.rank <- data.frame(
        rank(x = abs(x = cc.loadings[[1]][, cc.use])),
        rank(x = abs(x = cc.loadings[[g]][, cc.use])),
        cc.loadings[[1]][, cc.use],
        cc.loadings[[g]][, cc.use]
      )
      genes.rank$min <- apply(X = genes.rank[,1:2], MARGIN = 1, FUN = min)
      genes.rank <- genes.rank[order(genes.rank$min, decreasing = TRUE), ]
      genes.top <- rownames(x = genes.rank)[1:min(num.possible.genes, nrow(genes.rank))]
      bicors <- list()
      for (i in c(1, g)) {
        cc.vals <- cc.embeds[[i]][, cc.use]
        bicors[[i]] <- sapply(
          X = genes.top,
          FUN = function(x) {
            return(BiweightMidcor(x = cc.vals, y = scaled.data[[i]][x, ]))
          }
        )
      }
      genes.rank <- data.frame(
        rank(x = abs(x = bicors[[1]])),
        rank(x = abs(x = bicors[[g]])),
        bicors[[1]],
        bicors[[g]]
      )
      genes.rank$min <- apply(X = abs(x = genes.rank[, 1:2]), MARGIN = 1, FUN = min)
      # genes must be correlated in same direction in both datasets
      genes.rank <- genes.rank[sign(genes.rank[,3]) == sign(genes.rank[,4]), ]
      genes.rank <- genes.rank[order(genes.rank$min, decreasing = TRUE), ]
      genes.rank <- genes.rank[order(genes.rank$min, decreasing = TRUE), ]
      bc.gene[cc.use, g] <- mean(abs(genes.rank[1:gene.num, 4]))
      bc.gene.g1 <- c(bc.gene.g1, mean(abs(genes.rank[1:gene.num, 3])))
    }
    bc.gene[cc.use, 1] <- abs(mean(bc.gene.g1))
  }
  if (display.progress) {
    close(pb)
  }
  colnames(bc.gene) <- levels.split
  bc.gene <- as.data.frame(bc.gene)
  bc.gene$cc <- 1:nrow(bc.gene)
  bc.gene <- gather(bc.gene, key = "Group",  value = "bicor", -cc)
  return(bc.gene)
}

# Documentation needed
#Sets the parameter value for a plotting function as default for an object
# @param object ...
# @param function.name ...
# @param param.name ...
# @param param.value ...
#
# Quick example :   pbmc <- FixPlotParam(pbmc,"FeaturePlot","pt.size","0.1")
#
FixPlotParam <- function(object, function.name, param.name, param.value) {
  if ("PlotParams" %in% names(object@misc)) {
    object@misc[["PlotParams"]][paste(function.name, param.name, sep=":")] <- param.value
  }
  else {
    plotParamsList=list()
    plotParamsList[paste(function.name, param.name, sep=":")] <- param.value
    object@misc[["PlotParams"]] <- plotParamsList
  }
  return(object)
}

# Documentation needed
#
# @param object ...
# @param workflow.name ...
# @param depth ...
#
ReadPlotParams <- function(object, workflow.name, depth = 1) {
  if (!("PlotParams" %in% names(object@misc))) {
    return()
  }
  command.name <- as.character(deparse(sys.calls()[[sys.nframe()-depth]]))
  command.name <- gsub(pattern = ".Seurat",replacement = "",x = command.name)
  command.name <- ExtractField(string = command.name, field = 1, delim = "\\(")
  param.list <- names(formals(fun = sys.function(sys.parent(depth))))
  paramsList <- slot(object = object, name = "misc")[["PlotParams"]]
  p.env <- parent.frame(depth)
  for(i in names(paramsList)) {
    split_arr <- unlist(strsplit(x = i,split = ":"))
    function.name <- split_arr[1]
    param.name <- split_arr[2]
    param.value <- paramsList[[i]]
    if ((function.name == command.name) && (param.name %in% param.list)) {
      assign(x = param.name,
           value = param.value,
           envir = p.env)
    }
  }
  # overwrite any arguments passed in on the command line
  argnames <- sys.call(which = depth)
  argList <- as.list(argnames[-1])
  args_ignore <- c("", "object", "workflow.name")
  args_use <- setdiff(x = names(argList), y = args_ignore)
  for(i in args_use) {
    if(as.character(unlist(argList[i])[[1]]) == "F") {
      arg.val <- FALSE
    } else if(as.character(unlist(argList[i])[[1]]) == "T") {
      arg.val <- TRUE
    } else {
      arg.val <- unlist(argList[i])[[1]]
    }
    assign(x = i, value = arg.val, envir = p.env)
  }
}
