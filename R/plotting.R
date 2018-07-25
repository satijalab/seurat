#' Feature expression heatmap
#'
#' Draws a heatmap of single cell feature expression using superheat
#'
#' @param object Seurat object
#' @param cells.use A vector of cells to plot
#' @param features.use A vector of features to plot, defaults to \code{VariableFeatures(object = object)}
#' @param disp.min Minimum display value (all values below are clipped)
#' @param disp.max Maximum display value (all values above are clipped)
# @param group.by Name of variable to group cells by
#' @param slot.use Data slot to use, choose from 'raw.data', 'data', or 'scale.data'
#' @param assay.use Assay to pull from
#' @param check.plot Check that plotting will finish in a reasonable amount of time
#' @param ... Ignored for now
#'
#' @return Invisbly returns the final grob
#'
#' @export
#'
#' @seealso \code{\link{superheat::superheat}}
#'
#' @examples
#' DoHeatmap(object = pbmc_small)
#'
DoHeatmap <- function(
  object,
  cells.use = NULL,
  features.use = NULL,
  disp.min = -2.5,
  disp.max = 2.5,
  # group.by = "ident",
  slot.use = 'data',
  # group.order = NULL,
  # draw.line = TRUE,
  assay.use = NULL,
  check.plot = FALSE,
  ...
) {
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
    vars.fetch = features.use,
    cells.use = cells.use,
    slot = slot.use
  )
  return(SingleRasterMap(data.plot = data.plot))
  # plot.grob <- SingleHeatmap(
  #   object = object,
  #   cells.use = cells.use,
  #   features.use = features.use,
  #   group.by = group.by,
  #   disp.min = disp.min,
  #   disp.max = disp.max,
  #   slot.use = slot.use,
  #   check.plot = check.plot,
  #   ...
  # )
  # grid.newpage()
  # grid.draw(x = plot.grob)
  # invisible(x = plot.grob)
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
        cells$negative <- rev(x = cells$negative)
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
      SingleImageMap(data.plot = data.plot, plot.title = paste0(Key(object = object[[reduction.use]]), i))
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
      nCol = num.col,
      legend.position = 'right'
    )
  }
  return(plot.list)
  # for (i in 1:length(x = dims.use)) {
  #   dim.features <- unname(unlist(rev(features[[i]])))
  #   dim.features <- unlist(x = lapply(
  #     X = dim.features,
  #     FUN = function(feat) {
  #       return(grep(pattern = paste0(feat, '$'), x = features.keyed, value = TRUE))
  #     }
  #   ))
  #   plot.list[[i]] <- SingleHeatmap(
  #     object = object,
  #     cells.use = unname(unlist(cells.use[[i]])),
  #     features.use = dim.features,
  #     disp.min = disp.min,
  #     disp.max = disp.max,
  #     slot.use = slot.use,
  #     title = paste0(Key(object = object[[reduction.use]]), dims.use[i]),
  #     ...
  #   )
  # }
  # plot.grob <- arrangeGrob(
  #   grobs = plot.list,
  #   ncol = num.col
  # )
  # grid.newpage()
  # grid.draw(x = plot.grob)
  # invisible(x = plot.grob)
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
#' @param nCol Number of columns if multiple plots are displayed
#' @param do.sort Sort identity classes (on the x-axis) by the average
#' expression of the attribute being potted
#' @param y.max Maximum y axis value
#' @param same.y.lims Set all the y-axis limits to the same values
#' @param cols.use Colors to use for plotting
#' @param group.by Group (color) cells in different ways (for example, orig.ident)
#' @param y.log plot Y axis on log scale
#' @param combine.plots Combine plots into a single gg object; note that if TRUE; themeing will not work when plotting multiple features
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
  nCol = NULL,
  do.sort = FALSE,
  y.max = NULL,
  same.y.lims = FALSE,
  cols.use = NULL,
  group.by = NULL,
  y.log = FALSE,
  combine.plots = TRUE,
  ...
) {
  return(ExIPlot(
    object = object,
    plot.type = 'ridge',
    features.plot = features.plot,
    ident.include = ident.include,
    nCol = nCol,
    do.sort = do.sort,
    y.max = y.max,
    same.y.lims = same.y.lims,
    cols.use = cols.use,
    group.by = group.by,
    y.log = y.log,
    combine.plots = combine.plots,
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
  nCol = NULL,
  do.sort = FALSE,
  y.max = NULL,
  same.y.lims = FALSE,
  adjust.use = 1,
  point.size.use = 1,
  cols.use = NULL,
  group.by = NULL,
  y.log = FALSE,
  combine.plots = TRUE,
  ...
) {
  return(ExIPlot(
    object = object,
    plot.type = 'violin',
    features.plot = features.plot,
    ident.include = ident.include,
    nCol = nCol,
    do.sort = do.sort,
    y.max = y.max,
    same.y.lims = same.y.lims,
    adjust.use = adjust.use,
    point.size.use = point.size.use,
    cols.use = cols.use,
    group.by = group.by,
    y.log = y.log,
    combine.plots = combine.plots,
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
#' @importFrom ggplot2 scale_color_gradientn
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
  if (length(x = dims.use) != 2 || !is.numeric(x = dims.use)) {
    stop("'dims.use' must be a two-length integer vector")
  }
  dims.use <- paste0(Key(object = object[[reduction.use]]), dims.use)
  nCol <- NULL
  if (is.null(x = nCol)) {
    nCol <- 2
    if (length(x = features.plot) == 1) {
      nCol <- 1
    }
    if (length(x = features.plot) > 6) {
      nCol <- 3
    }
    if (length(x = features.plot) > 9) {
      nCol <- 4
    }
  }
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
  rownames(x = data.features) <- colnames(x = object)
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
    )
    if (brewer.gran == 2) {
      suppressMessages(expr = p <- p + scale_color_gradientn(
        colors = cols.use,
        guide = 'colorbar'
      ))
    }
    plots[[which(x = features.plot == feature)]] <- p
  }
  if (combine.plots) {
    plots <- CombinePlots(plot.list = plots, nCol = nCol, legend.position = NULL)
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
  dims.use = c(1, 2),
  pt.size = 2,
  cols.use = c("grey", "red"),
  reduction.use = "tsne",
  group.by = NULL,
  group.display = NULL,
  slot.use = 'data',
  scale.group = FALSE,
  min.cutoff = NA,
  max.cutoff = NA,
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
#' @param gene1 First feature to plot. Typically gene expression but can also
#' be metrics, PC scores, etc. - anything that can be retreived with FetchData
#' @param gene2 Second feature to plot.
#' @param cells.use Cells to include on the scatter plot.
#' @param cols.use Colors to use for identity class plotting.
#' @param pt.size Size of the points on the plot
#' @param shape.by Ignored for now
#' @param span Spline span in loess function call, if \code{NULL}, no spline added
#' @param ... Ignored for now
#'
#' @return A ggplot object
#'
#' @importFrom ggplot2 geom_smooth aes_string
#' @export
#'
#' @examples
#' GenePlot(object = pbmc_small, gene1 = 'CD9', gene2 = 'CD3E')
#'
GenePlot <- function(
  object,
  gene1,
  gene2,
  cells.use = NULL,
  cols.use = NULL,
  slot.use = 'data',
  pt.size = 1,
  shape.by = NULL,
  span = NULL,
  ...
) {
  cells.use <- cells.use %||% colnames(x = object)
  p <- SingleCorPlot(
    data.plot = FetchData(
      object = object,
      vars.fetch = c(gene1, gene2),
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
      mapping = aes_string(x = gene1, y = gene2),
      method = 'loess',
      span = span
    )
  }
  return(p)
}

#' Cell-cell scatter plot
#'
#' Creates a plot of scatter plot of genes across two single cells. Pearson
#' correlation between the two cells is displayed above the plot.
#'
#' @inheritParams GenePlot
#' @param cell1 Cell 1 name
#' @param cell2 Cell 2 name
#' @param features.use Genes to plot (default, all genes)
#' @param \dots Extra parameters passed to kde2d
#'
#' @return A ggplot object
#'
#' @export
#' @seealso \code{\link{MASS::kde2d}}
#'
#' @examples
#' CellPlot(object = pbmc_small, cell1 = 'ATAGGAGAAACAGA', cell2 = 'CATCAGGATGCACA')
#'
CellPlot <- function(
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
#' @param nCol Number of columns to display
#' @param do.balanced Return an equal number of genes with + and - scores. If FALSE (default), returns
#' the top genes ranked by the scores absolute values
#' @param combine.plots Combine plots into a single gg object; note that if TRUE; themeing will not work when plotting multiple features
#' @param ... Ignored
#'
#' @return A ggplot object
#'
#' @importFrom ggplot2 ggplot aes_string geom_point labs
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
  nCol = NULL,
  do.balanced = FALSE,
  combine.plots = TRUE,
  ...
) {
  if (is.null(x = nCol)) {
    nCol <- 2
    if (length(x = dims.use) == 1) {
      nCol <- 1
    }
    if (length(x = dims.use) > 6) {
      nCol <- 3
    }
    if (length(x = dims.use) > 9) {
      nCol <- 4
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
      data.plot$feature <- rownames(x = data.plot)
      plot <- ggplot(
        data = data.plot,
        mapping = aes_string(x = colnames(x = data.plot)[1], y = 'feature')
      ) +
        geom_point(col = pt.color) +
        labs(y = NULL)
      return(plot)
    }
  )
  if (combine.plots) {
    plots <- CombinePlots(plot.list = plots, nCol = nCol, legend.position = NULL)
  }
  return(plots)
}

globalVariables(names = 'Value', package = 'Seurat', add = TRUE)
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
#' @return Returns a Seurat object where object@@dr$pca@@jackstraw@@overall.p.values
#' represents p-values for each PC and object@@dr$pca@@misc$jackstraw.plot
#' stores the ggplot2 plot.
#'
#' @author Thanks to Omri Wurtzel for integrating with ggplot
#'
#' @importFrom reshape2 melt
#' @importFrom stats qqplot runif prop.test qunif
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
  pAll <- GetJS(object = GetDimReduc(object = object[[reduction.use]], slot = "jackstraw"), slot = "empirical.p.values")
  if (max(dims) > ncol(pAll)) {
    stop("Max dimension is ", ncol(pAll), ".")
  }
  pAll <- pAll[, dims, drop = FALSE]
  pAll <- as.data.frame(pAll)
  pAll$Contig <- rownames(x = pAll)
  pAll.l <- reshape2::melt(data = pAll, id.vars = "Contig")
  colnames(x = pAll.l) <- c("Contig", "PC", "Value")
  score.df <- GetJS(object = GetDimReduc(object = object[[reduction.use]], slot = "jackstraw"), slot = "overall.p.values")[dims, ]
  if (nrow(score.df) == 0) {
    stop(paste0("JackStraw hasn't been scored. Please run ScoreJackStraw before plotting."))
  }
  pAll.l$PC.Score <- rep(
    x = paste0("PC ", score.df[ ,"PC"], " ", sprintf("%1.3g", score.df[ ,"Score"])),
    each = length(x = unique(x = pAll.l$Contig))
  )
  pAll.l$PC.Score <- factor(
    x = pAll.l$PC.Score,
    levels = paste0("PC ", score.df[, "PC"], " ", sprintf("%1.3g", score.df[, "Score"]))
  )
  gp <- ggplot(data = pAll.l, mapping = aes(sample=Value)) +
    stat_qq(distribution = qunif) +
    facet_wrap("PC.Score") +
    labs(x = "Theoretical [runif(1000)]", y = "Empirical") +
    xlim(0, plot.y.lim) +
    ylim(0, plot.x.lim) +
    coord_flip() +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", na.rm = TRUE) +
    theme_bw()
  return(gp)
}
