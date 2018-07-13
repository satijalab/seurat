globalVariables(names = c('cell', 'gene'), package = 'Seurat', add = TRUE)
#' Gene expression heatmap
#'
#' Draws a heatmap of single cell gene expression using ggplot2.
#'
#' @param object Seurat object
#' @param data.use Option to pass in data to use in the heatmap. Default will pick from either
#' object@@data or object@@scale.data depending on use.scaled parameter. Should have cells as columns
#' and genes as rows.
#' @param use.scaled Whether to use the data or scaled data if data.use is NULL
#' @param cells.use Cells to include in the heatmap (default is all cells)
#' @param genes.use Genes to include in the heatmap (ordered)
#' @param disp.min Minimum display value (all values below are clipped)
#' @param disp.max Maximum display value (all values above are clipped)
#' @param group.by Groups cells by this variable. Default is object@@ident
#' @param group.order Order of groups from left to right in heatmap.
#' @param draw.line Draw vertical lines delineating different groups
#' @param col.low Color for lowest expression value
#' @param col.mid Color for mid expression value
#' @param col.high Color for highest expression value
#' @param slim.col.label display only the identity class name once for each group
#' @param remove.key Removes the color key from the plot.
#' @param rotate.key Rotate color scale horizantally
#' @param title Title for plot
#' @param cex.col Controls size of column labels (cells)
#' @param cex.row Controls size of row labels (genes)
#' @param group.label.loc Place group labels on bottom or top of plot.
#' @param group.label.rot Whether to rotate the group label.
#' @param group.cex Size of group label text
#' @param group.spacing Controls amount of space between columns.
#' @param assay.type to plot heatmap for (default is RNA)
#'
#' @return Returns a ggplot2 plot object
#'
#' @importFrom dplyr %>%
#' @importFrom reshape2 melt
#'
#' @export
#'
#' @examples
#' DoHeatmap(object = pbmc_small)
#'
DoHeatmap <- function(
  object,
  data.use = NULL,
  use.scaled = TRUE,
  cells.use = NULL,
  genes.use = NULL,
  disp.min = -2.5,
  disp.max = 2.5,
  group.by = "ident",
  group.order = NULL,
  draw.line = TRUE,
  col.low = "#FF00FF",
  col.mid = "#000000",
  col.high = "#FFFF00",
  slim.col.label = FALSE,
  remove.key = FALSE,
  rotate.key = FALSE,
  title = NULL,
  cex.col = 10,
  cex.row = 10,
  group.label.loc = "bottom",
  group.label.rot = FALSE,
  group.cex = 15,
  group.spacing = 0.15,
  assay.use = NULL
) {
  data.use <- data.use %||% GetAssayData(
    object = object,
    assay.use = assay.use,
    slot = ifelse(test = use.scaled, yes = 'scale.data', no = 'data'
    ))
  cells.use <- cells.use %||% colnames(x = object)
  cells.use <- intersect(x = cells.use, y = colnames(x = object))
  if (length(x = cells.use) == 0) {
    stop("No cells given to cells.use present in object")
  }
  genes.use <- genes.use %||% rownames(x = object)
  genes.use <- intersect(x = genes.use, y = rownames(x = object))
  if (length(x = genes.use) == 0) {
    stop("No genes given to genes.use present in object")
  }
  if (is.null(x = group.by) || group.by == "ident") {
    cells.ident <- Idents(object = object)[cells.use]
  } else {
    cells.ident <- factor(x = unlist(x = object[group.by], use.names = FALSE))
    names(x = cells.ident) <- cells.use
  }
  cells.ident <- factor(
    x = cells.ident,
    labels = intersect(x = levels(x = cells.ident), y = cells.ident)
  )
  data.use <- data.use[genes.use, cells.use, drop = FALSE]
  if (!use.scaled) {
    data.use <- as.matrix(x = data.use)
    disp.max <- ifelse(test = disp.max == 2.5, yes = 10, no = disp.max)
  }
  data.use <- MinMax(data = data.use, min = disp.min, max = disp.max)
  data.use <- as.data.frame(x = t(x = data.use))
  data.use$cell <- rownames(x = data.use)
  colnames(x = data.use) <- make.unique(names = colnames(x = data.use))
  data.use %>% melt(id.vars = "cell") -> data.use
  names(x = data.use)[names(x = data.use) == 'variable'] <- 'gene'
  names(x = data.use)[names(x = data.use) == 'value'] <- 'expression'
  data.use$ident <- cells.ident[data.use$cell]
  if (!is.null(x = group.order)) {
    if (length(group.order) == length(levels(data.use$ident)) && all(group.order %in% levels(data.use$ident))) {
      data.use$ident <- factor(data.use$ident, levels = group.order)
    }
    else {
      stop("Invalid group.order")
    }
  }
  data.use$gene <- with(
    data = data.use,
    expr = factor(x = gene, levels = rev(x = unique(x = data.use$gene)))
  )
  data.use$cell <- with(
    data = data.use,
    expr = factor(x = cell, levels = cells.use)
  )
  if (rotate.key) {
    key.direction <- "horizontal"
    key.title.pos <- "top"
  } else {
    key.direction <- "vertical"
    key.title.pos <- "left"
  }
  heatmap <- ggplot(
    data = data.use,
    mapping = aes(x = cell, y = gene, fill = expression)
  ) +
    geom_tile() +
    scale_fill_gradient2(
      low = col.low,
      mid = col.mid,
      high = col.high,
      name = "Expression",
      guide = guide_colorbar(
        direction = key.direction,
        title.position = key.title.pos
      )
    ) +
    scale_y_discrete(position = "right", labels = rev(genes.use)) +
    theme(
      axis.line = element_blank(),
      axis.title.y = element_blank(),
      axis.ticks.y = element_blank(),
      strip.text.x = element_text(size = group.cex),
      axis.text.y = element_text(size = cex.row),
      axis.text.x = element_text(size = cex.col),
      axis.title.x = element_blank()
    )
  if (slim.col.label) {
    heatmap <- heatmap +
      theme(
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank()
      )
  } else {
    heatmap <- heatmap + theme(axis.text.x = element_text(angle = 90))
  }
  if (!is.null(x = group.by)) {
    if (group.label.loc == "top") {
      switch <- NULL
    } else {
      switch <- 'x'
    }
    heatmap <- heatmap +
      facet_grid(
        facets = ~ident,
        drop = TRUE,
        space = "free",
        scales = "free",
        switch = switch
      ) +
      scale_x_discrete(expand = c(0, 0), drop = TRUE)
    if (draw.line) {
      panel.spacing <- unit(x = group.spacing, units = 'lines')
    } else {
      panel.spacing <- unit(x = 0, units = 'lines')
    }
    heatmap <- heatmap +
      theme(strip.background = element_blank(), panel.spacing = panel.spacing)
    if (group.label.rot) {
      heatmap <- heatmap + theme(strip.text.x = element_text(angle = 90))
    }
  }
  if (remove.key) {
    heatmap <- heatmap + theme(legend.position = "none")
  }
  if (!is.null(x = title)) {
    heatmap <- heatmap + labs(title = title)
  }
  return(heatmap)
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
#' @return A list of ggplot objects
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
#' @return A list of ggplot objects
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

globalVariables(
  names = c('cell', 'id', 'avg.exp', 'avg.exp.scale', 'pct.exp'),
  package = 'Seurat',
  add = TRUE
)
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
#' @param plot.legend plots the legends
#' @param x.lab.rot Rotate x-axis labels
#'
#' @return default, no return, only graphical output. If do.return=TRUE, returns a ggplot2 object
#'
#' @importFrom tidyr gather
#' @importFrom dplyr %>% group_by summarize mutate ungroup
#'
#' @export
#' @seealso \code{RColorBrewer::brewer.pal.info}
#'
#' @examples
#' cd_genes <- c("CD247", "CD3E", "CD9")
#' DotPlot(object = pbmc_small, genes.plot = cd_genes)
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
  plot.legend = FALSE,
  x.lab.rot = FALSE
) {
  scale.func <- switch(
    EXPR = scale.by,
    'size' = scale_size,
    'radius' = scale_radius,
    stop("'scale.by' must be either 'size' or 'radius'")
  )
  data.to.plot <- FetchData(object = object, vars.fetch = features.plot)
  colnames(x = data.to.plot) <- features.plot
  data.to.plot$cell <- rownames(x = data.to.plot)
  data.to.plot$id <- if (is.null(x = group.by)) {
    Idents(object = object)
  } else {
    object[group.by]
  }
  data.to.plot %>% gather(
    key = features.plot,
    value = expression,
    -c(cell, id)
  ) -> data.to.plot
  data.to.plot %>%
    group_by(id, features.plot) %>%
    summarize(
      avg.exp = mean(expm1(x = expression)),
      pct.exp = PercentAbove(x = expression, threshold = 0)
    ) -> data.to.plot
  data.to.plot %>%
    ungroup() %>%
    group_by(features.plot) %>%
    mutate(avg.exp.scale = scale(x = avg.exp)) %>%
    mutate(avg.exp.scale = MinMax(
      data = avg.exp.scale,
      max = col.max,
      min = col.min
    )) ->  data.to.plot
  data.to.plot$features.plot <- factor(
    x = data.to.plot$features.plot,
    levels = rev(x = features.plot)
  )
  data.to.plot$pct.exp[data.to.plot$pct.exp < dot.min] <- NA
  data.to.plot$pct.exp <- data.to.plot$pct.exp * 100
  p <- ggplot(data = data.to.plot, mapping = aes(x = features.plot, y = id)) +
    geom_point(mapping = aes(size = pct.exp, color = avg.exp.scale)) +
    scale.func(range = c(0, dot.scale), limits = c(scale.min, scale.max)) +
    theme(axis.title.x = element_blank(), axis.title.y = element_blank())
  if (length(x = cols.use) == 1) {
    p <- p + scale_color_distiller(palette = cols.use)
  } else {
    p <- p + scale_color_gradient(low = cols.use[1], high = cols.use[2])
  }
  if (!plot.legend) {
    p <- p + theme(legend.position = "none")
  }
  if (x.lab.rot) {
    p <- p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
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
# @param do.return Return a ggplot2 object (default : FALSE)
# @param do.bare Do only minimal formatting (default : FALSE)
#' @param cols.use Vector of colors, each color corresponds to an identity
#' class. By default, ggplot assigns colors.
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
#'in \code{cols.highlight} and other cells black (white if dark.theme = TRUE);
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
#' @return If do.return==TRUE, returns a ggplot2 object. Otherwise, only
#' graphical output.
#'
# @seealso \code{\link{FeatureLocator}}
#'
#' @import SDMTools
#' @importFrom stats median
#' @importFrom dplyr summarize group_by
#' @importFrom png readPNG
#'
#' @export
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
  group.by <- group.by %||% 'ident'
  col.by <- switch(
    EXPR = group.by,
    'ident' = Idents(object = object),
    object[group.by, drop = TRUE]
  )
  if (!is.null(x = shape.by)) {
    shape.by <- object[shape.by, drop = TRUE]
  }
  plot <- SingleDimPlot(
    object = object[[reduction.use]],
    dims.use = dims.use,
    col.by = col.by,
    cols.use = cols.use,
    cells.use = cells.use,
    pt.size = pt.size,
    shape.by = shape.by,
    plot.order = plot.order,
    do.label = do.label,
    label.size = label.size,
    cells.highlight = cells.highlight,
    cols.highlight = cols.highlight,
    sizes.highlight = sizes.highlight,
    legend.title = group.by
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
#'
#' @importFrom RColorBrewer brewer.pal.info
#' @importFrom ggplot2 scale_color_gradientn
#' @importFrom cowplot plot_grid
#'
#' @return No return value, only a graphical output
#'
#' @export
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
  na.value = 'grey50'
) {
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
  data.features <- FetchData(object = object, vars.fetch = features.plot)
  features.plot <- colnames(x = data.features)
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
  data.features <- sapply(
    X = 1:ncol(x = data.features),
    FUN = function(index) {
      data.feature <- as.vector(x = data.features[, index])
      min.use <- SetQuantile(cutoff = min.cutoff[index], data.feature)
      max.use <- SetQuantile(cutoff = max.cutoff[index], data.feature)
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
  colnames(x = data.features) <- features.plot
  rownames(x = data.features) <- colnames(x = object)
  plot.list <- vector(mode = 'list', length = length(x = features.plot))
  for (i in 1:length(x = features.plot)) {
    p <- SingleDimPlot(
      object = object[[reduction.use]],
      dims.use = dims.use,
      col.by = as.vector(x = data.features[, i]),
      cols.use = cols.use,
      cells.use = cells.use,
      pt.size = pt.size,
      do.label = do.label,
      label.size = label.size,
      legend.title = features.plot[i]
    )
    if (brewer.gran == 2) {
      suppressMessages(expr = p <- p + scale_color_gradientn(
        colors = cols.use,
        guide = 'colorbar'
      ))
    }
    plot.list[[i]] <- p
  }
  plots.combined <- if (length(x = plot.list) > 1) {
    plot_grid(plotlist = plot.list, ncol = nCol)
  } else {
    plot.list[[1]]
  }
  return(plots.combined)
}
