#' Gene expression heatmap
#'
#' Draws a heatmap of single cell gene expression.
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
#' @param draw.line Draw vertical lines delineating different groups
#' @param col.low Color for lowest expression value
#' @param col.mid Color for mid expression value
#' @param col.high Color for highest expression value
#' @param slim.col.label display only the identity class name once for each group
#' @param remove.key Removes the color key from the plot.
#' @param rotate.key Rotate color scale horizantally
#' @param cex.col Controls size of column labels (cells)
#' @param cex.row Controls size of row labels (genes)
#' @param group.label.loc Place group labels on bottom or top of plot.
#' @param group.label.rot Whether to rotate the group label.
#' @param group.cex Size of group label text
#' @param group.spacing Controls amount of space between columns.
#' @param do.plot Whether to display the plot.
#' @return Returns a ggplot2 plot object
#' @importFrom reshape2 melt
#' @importFrom dplyr %>%
#' @export
DoHeatmapGG <- function(
  object,
  data.use = NULL,
  use.scaled = TRUE,
  cells.use = NULL,
  genes.use = NULL,
  disp.min = -2.5,
  disp.max = 2.5,
  group.by = "ident",
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
  do.plot = TRUE,
  ...
) {
  if (is.null(x = data.use)) {
    if (use.scaled) {
      data.use <- object@scale.data
    } else {
      data.use <- object@data
    }
  }
  # note: data.use should have cells as column names, genes as row names
  cells.use <- set.ifnull(x = cells.use, y = object@cell.names)
  cells.use <- ainb(a = cells.use, b = colnames(x = data.use))
  if (length(x = cells.use) == 0) {
    stop("No cells given to cells.use present in object")
  }
  genes.use <- set.ifnull(x = genes.use, rownames(y = data.use))
  genes.use <- ainb(a = genes.use, b = rownames(x = data.use))
  if (length(x = genes.use) == 0) {
    stop("No genes given to genes.use present in object")
  }
  if (is.null(x = group.by) || group.by == "ident") {
    cells.ident <- object@ident[cells.use]
  } else {
    cells.ident <- factor(x = FetchData(
      object = object,
      cells.use = cells.use,
      vars.all = group.by
    )[, 1])
    names(x = cells.ident) <- cells.use
  }
  cells.ident <- factor(
    x = cells.ident,
    labels = ainb(a = levels(x = cells.ident), b = cells.ident)
  )
  data.use <- data.use[genes.use, cells.use]
  if (use.scaled) {
    data.use <- minmax(data = data.use, min = disp.min, max = disp.max)
  }
  data.use <- as.data.frame(x = t(x = data.use))
  data.use$cell <- rownames(x = data.use)
  colnames(x = data.use) <- make.unique(names = colnames(x = data.use))
  data.use %>% melt(id.vars = "cell") -> data.use
  names(x = data.use)[names(x = data.use) == 'variable'] <- 'gene'
  names(x = data.use)[names(x = data.use) == 'value'] <- 'expression'
  data.use$ident <- cells.ident[data.use$cell]
  breaks <- seq(
    from = min(data.use$expression),
    to = max(data.use$expression),
    length = length(x = pyCols) + 1
  )
  data.use$gene <- with(
    data = data.use,
    expr = factor(x = gene, levels = rev(x = unique(x = data.use$gene)))
  )
  data.use$cell <- with(
    data = data.use,
    expr = factor(x = cell, levels = cells.use)
  )
  # might be a solution if we want discrete interval units, makes the legend clunky though
  #data.use$expression <- cut(data.use$expression, breaks = breaks, include.lowest = T)
  #heatmap <- ggplot(data.use, aes(x = cell, y = gene, fill = expression)) + geom_tile() +
  #                  scale_fill_manual(values = pyCols, name= "Expression") +
  #                  scale_y_discrete(position = "right", labels = rev(genes.use)) +
  #                  theme(axis.line=element_blank(), axis.title.y=element_blank(),
  #                        axis.ticks.y = element_blank())
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
      name= "Expression",
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
  if (! is.null(x = group.by)) {
    if (group.label.loc == "top") {
      switch <- NULL
      # heatmap <- heatmap +
      #   facet_grid(
      #     facets = ~ident,
      #     drop = TRUE,
      #     space = "free",
      #     scales = "free"
      #   ) +
      #   scale_x_discrete(expand = c(0, 0), drop = TRUE)
    } else {
      switch <- 'x'
      # heatmap <- heatmap +
      #   facet_grid(
      #     facets = ~ident,
      #     drop = TRUE,
      #     space = "free",
      #     scales = "free",
      #     switch = "x"
      #   ) +
      #   scale_x_discrete(expand = c(0, 0), drop = TRUE)
    }
    heatmap <- heatmap +
      facet_grid(
        facets = ~ident,
        drop = TRUE,
        space = "free",
        scales = "free",
        switch = switch,
      ) +
      scale_x_discrete(expand = c(0, 0), drop = TRUE)
    if (draw.line) {
      panel.spacing <- unit(x = group.spacing, units = 'lines')
      # heatmap <- heatmap + theme(strip.background = element_blank(), panel.spacing = unit(group.spacing, "lines"))
    } else {
      panel.spacing <- unit(x = 0, units = 'lines')
      #
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
  if (! is.null(x = title)) {
    heatmap <- heatmap + labs(title = title)
  }
  if (do.plot) {
    heatmap
  }
  return(heatmap)
}

#' Single cell violin plot
#'
#' Draws a violin plot of single cell data (gene expression, metrics, PC
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
#' @param size.x.use X axis title font size
#' @param size.y.use Y axis title font size
#' @param size.title.use Main title font size
#' @param adjust.use Adjust parameter for geom_violin
#' @param point.size.use Point size for geom_violin
#' @param cols.use Colors to use for plotting
#' @param group.by Group (color) cells in different ways (for example, orig.ident)
#' @param y.log plot Y axis on log scale
#' @param x.lab.rot Rotate x-axis labels
#' @param y.lab.rot Rotate y-axis labels
#' @param legend.position Position the legend for the plot
#' @param single.legend Consolidate legend the legend for all plots
#' @param remove.legend Remove the legend from the plot
#' @param do.return Return a ggplot2 object (default : FALSE)
#' @param return.plotlist Return the list of individual plots instead of compiled plot.
#' @param \dots additional parameters to pass to FetchData (for example, use.imputed, use.scaled, use.raw)
#' @import ggplot2
#' @importFrom cowplot plot_grid
#' @return By default, no return, only graphical output. If do.return=TRUE,
#' returns a list of ggplot objects.
#' @export
VlnPlot <- function(
  object,
  features.plot,
  ident.include = NULL,
  nCol = NULL,
  do.sort = FALSE,
  y.max = NULL,
  same.y.lims = FALSE,
  size.x.use = 16,
  size.y.use = 16,
  size.title.use = 20,
  adjust.use = 1,
  point.size.use = 1,
  cols.use = NULL,
  group.by = NULL,
  y.log = FALSE,
  x.lab.rot = FALSE,
  y.lab.rot = FALSE,
  legend.position = "right",
  single.legend = TRUE,
  remove.legend = FALSE,
  do.return = FALSE,
  return.plotlist = FALSE,
  ...
) {
  if (is.null(x = nCol)) {
    if (length(x = features.plot) > 9) {
      nCol <- 4
    } else {
      nCol <- min(length(x = features.plot), 3)
    }
  }
  data.use <- data.frame(FetchData(object = object, vars.all = features.plot, ...))
  if (is.null(x = ident.include)) {
    cells.to.include <- object@cell.names
  } else {
    cells.to.include <- WhichCells(object = object, ident = ident.include)
  }
  data.use <- data.use[cells.to.include, ,drop = FALSE]
  if (!is.null(x = group.by)) {
    ident.use <- as.factor(x = FetchData(
      object = object,
      vars.all = group.by
    )[cells.to.include, 1])
  } else {
    ident.use <- object@ident[cells.to.include]
  }
  gene.names <- colnames(x = data.use)[colnames(x = data.use) %in% rownames(x = object@data)]
  if (single.legend) {
    remove.legend <- TRUE
  }
  if (same.y.lims && is.null(x = y.max)) {
    y.max <- max(data.use)
  }
  plots <- lapply(
    X = features.plot,
    FUN = function(x) {
      return(PlotVln(
        feature = x,
        data = data.use[, x, drop = FALSE],
        cell.ident = ident.use,
        do.sort = do.sort, y.max = y.max,
        size.x.use = size.x.use,
        size.y.use = size.y.use,
        size.title.use = size.title.use,
        adjust.use = adjust.use,
        point.size.use = point.size.use,
        cols.use = cols.use,
        gene.names = gene.names,
        y.log = y.log,
        x.lab.rot = x.lab.rot,
        y.lab.rot = y.lab.rot,
        legend.position = legend.position,
        remove.legend = remove.legend
      ))
    }
  )
  if (length(x = features.plot) > 1) {
    plots.combined <- plot_grid(plotlist = plots, ncol = nCol)
    if (single.legend && !remove.legend) {
      legend <- get_legend(
        plot = plots[[1]] + theme(legend.position = legend.position)
      )
      if (legend.position == "bottom") {
        plots.combined <- plot_grid(
          plots.combined,
          legend,
          ncol = 1,
          rel_heights = c(1, .2)
        )
      } else if (legend.position == "right") {
        plots.combined <- plot_grid(
          plots.combined,
          legend,
          rel_widths = c(3, .3)
        )
      } else {
        warning("Shared legends must be at the bottom or right of the plot")
      }
    }
  } else {
    plots.combined <- plots[[1]]
  }
  if (do.return) {
    if (return.plotlist) {
      return(plots)
    } else {
      return(plots.combined)
    }
  } else {
    if (length(x = plots.combined) > 1) {
      plots.combined
    }
    else {
      invisible(x = lapply(X = plots.combined, FUN = print))
    }
  }
}

#' Dot plot visualization
#'
#' Intuitive way of visualizing how gene expression changes across different identity classes (clusters).
#' The size of the dot encodes the percentage of cells within a class, while the color encodes the
#' AverageExpression level of 'expressing' cells (green is high).
#'
#' @param object Seurat object
#' @param genes.plot Input vector of genes
#' @param cex.use Scaling factor for the dots (scales all dot sizes)
#' @param cols.use colors to plot
#' @param thresh.col The raw data value which corresponds to a red dot (lowest expression)
#' @param dot.min The fraction of cells at which to draw the smallest dot (default is 0.05)
#' @param group.by Factor to group the cells by
#' @return Only graphical output
#' @export
DotPlot <- function(
  object,
  genes.plot,
  cex.use = 2,
  cols.use = NULL,
  thresh.col = 2.5,
  dot.min = 0.05,
  group.by = NULL,
  ...
) {
  if (! is.null(x = group.by)) {
    object <- SetAllIdent(object = object, id = group.by)
  }
  #object@data=object@data[genes.plot,]
  object@data <- data.frame(t(x = FetchData(object = object, vars.all = genes.plot)))
  #this line is in case there is a '-' in the cell name
  colnames(x = object@data) <- object@cell.names
  avg.exp <- AverageExpression(object = object)
  avg.alpha <- ClusterAlpha(object = object)
  cols.use <- set.ifnull(x = cols.use, y = myPalette(low = "red", high = "green"))
  exp.scale <- t(x = scale(x = t(x = avg.exp)))
  exp.scale <- minmax(data = exp.scale, max = thresh.col, min = (-1) * thresh.col)
  n.col <- length(x = cols.use)
  data.y <- rep(x = 1:ncol(x = avg.exp), nrow(x = avg.exp))
  data.x <- unlist(x = lapply(X = 1:nrow(x = avg.exp), FUN = rep, ncol(x = avg.exp)))
  data.avg <- unlist(x = lapply(
    X = 1:length(x = data.y),
    FUN = function(x) {
      return(exp.scale[data.x[x], data.y[x]])
    }
  ))
  exp.col <- cols.use[floor(
    x = n.col * (data.avg + thresh.col) / (2 * thresh.col) + .5
  )]
  data.cex <- unlist(x = lapply(
    X = 1:length(x = data.y),
    FUN = function(x) {
      return(avg.alpha[data.x[x], data.y[x]])
    }
  )) * cex.use + dot.min
  plot(
    x = data.x,
    y = data.y,
    cex = data.cex,
    pch = 16,
    col = exp.col,
    xaxt = "n",
    xlab = "",
    ylab = "",
    yaxt = "n"
  )
  axis(side = 1, at = 1:length(x = genes.plot), labels = genes.plot)
  axis(side = 2, at = 1:ncol(x = avg.alpha), colnames(x = avg.alpha), las = 1)
}

#' Dot plot visualization
#'
#' Intuitive way of visualizing how gene expression changes across different identity classes (clusters).
#' The size of the dot encodes the percentage of cells within a class, while the color encodes the
#' AverageExpression level of 'expressing' cells (green is high).
#'
#' @param object Seurat object
#' @param genes.plot Input vector of genes
#' @param cols.use colors to plot
#' @param col.min Minimum scaled average expression threshold (everything smaller will be set to this)
#' @param col.max Maximum scaled average expression threshold (everything larger will be set to this)
#' @param dot.min The fraction of cells at which to draw the smallest dot (default is 0.05).
#' @param dot.scale Scale the size of the points, similar to cex
#' @param group.by Factor to group the cells by
#' @param plot.legend plots the legends
#' @param x.lab.rot Rotate x-axis labels
#' @param do.return Return ggplot2 object
#' @return default, no return, only graphical output. If do.return=TRUE, returns a ggplot2 object
#' @importFrom dplyr %>% group_by summarize_each mutate ungroup
#' @importFrom tidyr gather
#' @export
DotPlotGG <- function(
  object,
  genes.plot,
  cols.use = c("green", "red"),
  col.min = -2.5,
  col.max = 2.5,
  dot.min = 0,
  dot.scale = 6,
  group.by,
  plot.legend = FALSE,
  do.return = FALSE,
  x.lab.rot = FALSE
) {
  if (! missing(x = group.by)) {
    object <- SetAllIdent(object = object, id = group.by)
  }
  data.to.plot <- data.frame(FetchData(object = object, vars.all = genes.plot))
  data.to.plot$cell <- rownames(x = data.to.plot)
  data.to.plot$id <- object@ident
  data.to.plot %>% gather(
    key = genes.plot,
    value = expression,
    -c(cell, id)
  ) -> data.to.plot
  data.to.plot %>%
    group_by(id, genes.plot) %>%
    summarize(
      avg.exp = expMean(x = expression),
      pct.exp = PercentAbove(x = expression, threshold = 0)
    ) -> data.to.plot
  data.to.plot %>%
    ungroup() %>%
    group_by(genes.plot) %>%
    mutate(avg.exp = as.numeric(x = scale(center = avg.exp))) %>%
    mutate(avg.exp.scale = minmax(
      data = avg.exp,
      max = col.max,
      min = col.min
    )) ->  data.to.plot
  data.to.plot$genes.plot <- factor(
    x = data.to.plot$genes.plot,
    levels = rev(x = sub(pattern = "-", replacement = ".", x = genes.plot))
  )
  data.to.plot$pct.exp[data.to.plot$pct.exp < dot.min] <- NA
  p <- ggplot(data = data.to.plot, mapping = aes(x = genes.plot, y = id)) +
    geom_point(mapping = aes(size = pct.exp, color = avg.exp.scale)) +
    scale_radius(range = c(0, dot.scale)) +
    scale_color_gradient(low = cols.use[1], high = cols.use[2]) +
    theme(axis.title.x = element_blank(), axis.title.y = element_blank())
  if (! plot.legend) {
    p <- p + theme(legend.position = "none")
  }
  if (x.lab.rot) {
    p <- p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
  }
  suppressWarnings(print(p))
  if (do.return) {
    return(p)
  }
}


#' Split Dot plot visualization
#'
#' Intuitive way of visualizing how gene expression changes across different identity classes (clusters).
#' The size of the dot encodes the percentage of cells within a class, while the color encodes the
#' AverageExpression level of 'expressing' cells (green is high). Splits the cells into two groups based on a
#' grouping variable.
#' Still in BETA
#'
#' @param object Seurat object
#' @param grouping.var Grouping variable for splitting the dataset
#' @param genes.plot Input vector of genes
#' @param cols.use colors to plot
#' @param col.min Minimum scaled average expression threshold (everything smaller will be set to this)
#' @param col.max Maximum scaled average expression threshold (everything larger will be set to this)
#' @param dot.min The fraction of cells at which to draw the smallest dot (default is 0.05).
#' @param dot.scale Scale the size of the points, similar to cex
#' @param group.by Factor to group the cells by
#' @param plot.legend plots the legends
#' @param x.lab.rot Rotate x-axis labels
#' @param do.return Return ggplot2 object
#' @param gene.groups Add labeling bars to the top of the plot
#' @return default, no return, only graphical output. If do.return=TRUE, returns a ggplot2 object
#' @importFrom dplyr %>% group_by summarize_each mutate ungroup
#' @importFrom tidyr gather
#' @export
SplitDotPlotGG <- function(
  object,
  grouping.var,
  genes.plot,
  gene.groups,
  cols.use = c("green", "red"),
  col.min = -2.5,
  col.max = 2.5,
  dot.min = 0,
  dot.scale = 6,
  group.by,
  plot.legend = FALSE,
  do.return = FALSE,
  x.lab.rot = FALSE
) {
  if (! missing(x = group.by)) {
    object <- SetAllIdent(object = object, id = group.by)
  }
  grouping.data <- FetchData(
    object = object,
    vars.all = grouping.var
  )[names(x = object@ident), 1]
  idents.old <- levels(x = object@ident)
  object@ident <- paste(object@ident, grouping.data, sep="_")
  object@ident <- factor(
    x = object@ident,
    levels = unlist(x = lapply(
      X = idents.old,
      FUN = function(x) {
        return(c(
          paste(x, unique(x = grouping.data)[1], sep="_"),
          paste(x, unique(x = grouping.data)[2], sep="_")
        ))
      }
    )),
    ordered = TRUE
  )
  data.to.plot <- data.frame(FetchData(object = object, vars.all = genes.plot))
  data.to.plot$cell <- rownames(x = data.to.plot)
  data.to.plot$id <- object@ident
  data.to.plot %>%
    gather(key = genes.plot, value = expression, -c(cell, id)) -> data.to.plot
  data.to.plot %>%
    group_by(id, genes.plot) %>%
    summarize(
      avg.exp = expMean(x = expression),
      pct.exp = PercentAbove(x = expression, threshold = 0)
    ) -> data.to.plot
  ids.2 <- paste(
    idents.old,
    as.character(x = unique(x = grouping.data)[2]),
    sep = "_"
  )
  vals.2 <- which(x = data.to.plot$id %in% ids.2)
  ids.1 <- paste(
    idents.old,
    as.character(x = unique(x = grouping.data)[1]),
    sep = "_"
  )
  vals.1 <- which(x = data.to.plot$id %in% ids.1)
  #data.to.plot[vals.2,3]=-1*data.to.plot[vals.2,3]
  data.to.plot %>%
    ungroup() %>%
    group_by(genes.plot) %>%
    mutate(avg.exp = scale(avg.exp)) %>%
    mutate(avg.exp.scale = as.numeric(x = cut(
      x = minmax(data = avg.exp, max = col.max, min = col.min),
      breaks = 20
    ))) ->  data.to.plot
  data.to.plot$genes.plot <- factor(
    x = data.to.plot$genes.plot,
    levels = rev(x = sub(pattern = "-", replacement = ".", x = genes.plot))
  )
  data.to.plot$pct.exp[data.to.plot$pct.exp < dot.min] <- NA
  palette.1 <- myPalette(low = "grey", high = "blue", k = 20)
  palette.2 <- myPalette(low = "grey", high = "red", k = 20)
  data.to.plot$ptcolor <- "grey"
  data.to.plot[vals.1, "ptcolor"] <- palette.1[as.matrix(
    x = data.to.plot[vals.1, "avg.exp.scale"]
  )[, 1]]
  data.to.plot[vals.2, "ptcolor"] <- palette.2[as.matrix(
    x = data.to.plot[vals.2, "avg.exp.scale"]
  )[, 1]]
  if (! missing(x = gene.groups)) {
    names(x = gene.groups) <- genes.plot
    data.to.plot %>%
      mutate(gene.groups = gene.groups[genes.plot]) -> data.to.plot
  }
  p <- ggplot(data = data.to.plot, mapping = aes(x = genes.plot, y = id)) +
    geom_point(mapping = aes(size = pct.exp, color = ptcolor)) +
    scale_radius(range = c(0, dot.scale)) +
    scale_color_identity() +
    theme(axis.title.x = element_blank(), axis.title.y = element_blank())
  if (! missing(x = gene.groups)) {
    p <- p +
      facet_grid(
        facets = ~gene.groups,
        scales = "free_x",
        space = "free_x",
        switch = "y"
      ) +
      theme(
        panel.spacing = unit(x = 1, units = "lines"),
        strip.background = element_blank(),
        strip.placement = "outside"
      )
  }
  if (! plot.legend) {
    p <- p + theme(legend.position = "none")
  }
  if (x.lab.rot) {
    p <- p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
  }
  suppressWarnings(print(p))
  if (do.return) {
    return(p)
  }
}

#' Dark Theme
#'
#' Add a dark theme to ggplot objects
#'
#' @param ... Extra parameters to be passed to theme()
#' @import ggplot2
#' @return A ggplot2 theme object
#' @seealso \code{\link{theme}}
#' @export
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
    #   Make this a complete theme and validate it
    complete = TRUE,
    validate = TRUE,
    #   Extra parameters
    ...
  )
  return(dark.theme)
}

#   Functions for converting ggplot2 objects
#   to standard plots for use with locator
PlotBuild <- function(plot.data, dark.theme = FALSE, smooth = FALSE, ...) {
  #   Do we use a smooth scatterplot?
  #   Take advantage of functions as first class objects
  #   to dynamically choose normal vs smooth scatterplot
  if (smooth) {
    myplot <- smoothScatter
  } else {
    myplot <- plot
  }
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

GGpointToBase <- function(plot, do.plot = TRUE, ...) {
  plot.build <- ggplot2::ggplot_build(plot = plot)
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

#   Locate points on a plot and return them
PointLocator <- function(plot, recolor=TRUE, dark.theme = FALSE, ...) {
  #   Convert the ggplot object to a data.frame
  plot.data <- GGpointToBase(plot = plot, dark.theme = dark.theme, ...)
  npoints <- nrow(x = plot.data)
  cat("Click around the cluster of points you wish to select\n")
  cat("ie. select the vertecies of a shape around the cluster you\n")
  cat("are interested in. Press <Esc> when finished (right click for R-terminal users)\n\n")
  polygon <- locator(n = npoints, type = 'l')
  polygon <- data.frame(polygon)
  #   pnt.in.poly returns a data.frame of points
  points.all <- SDMTools::pnt.in.poly(
    pnts = plot.data[, c(1, 2)],
    poly.pnts = polygon
  )
  #   Find the located points
  points.located <- points.all[which(x = points.all$pip == 1), ]
  #   If we're recoloring, do the recolor
  if(recolor) {
    if (dark.theme) {
      no = 'white'
    } else {
      no = 'black'
    }
    points.all$color <- ifelse(test = points.all$pip == 1, yes = 'red', no = no)
    plot.data$color <- points.all$color
    PlotBuild(plot.data = plot.data, dark.theme = dark.theme, ...)
  }
  return(points.located[, c(1, 2)])
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
#' @seealso \code{\link{locator}}
#' @seealso \code{\link{ggplot2::ggplot_build}}
#' @export
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
#' @seealso \code{\link{plotly::layout}}
#' @seealso \code{\link{ggplot2::ggplot_build}}
#' @export
#'
HoverLocator <- function(plot, data.plot, features.info = NULL, dark.theme = FALSE, ...) {
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
      FUN = paste,
      collapse = '</br>'
    )
    data.info <- data.frame(
      feature = paste(rownames(x = features.info), info, sep = '</br>'),
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

PlotVln <- function(
  feature,
  data,
  cell.ident,
  do.sort,
  y.max,
  size.x.use,
  size.y.use,
  size.title.use,
  adjust.use,
  point.size.use,
  cols.use,
  gene.names,
  y.log,
  x.lab.rot,
  y.lab.rot,
  legend.position,
  remove.legend
) {
  set.seed(seed = 42)
  data$ident <- cell.ident
  if(do.sort) {
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
  data[, feature] <- data[, feature] + noise
  y.max <- set.ifnull(x = y.max, y = max(data[, feature]))
  plot <- ggplot(
    data = data,
    mapping = aes(
      x = factor(x = ident),
      y = eval(expr = parse(text = feature))
    )
  ) +
    geom_violin(
      scale = "width",
      adjust = adjust.use,
      trim = TRUE,
      mapping = aes(fill = factor(x = ident))
    ) +
    theme(
      legend.position = legend.position,
      axis.title.x = element_text(
        face = "bold",
        colour = "#990000",
        size = size.x.use
      ),
      axis.title.y = element_text(
        face = "bold",
        colour = "#990000",
        size = size.y.use
      )
    ) +
    guides(fill = guide_legend(title = NULL)) +
    geom_jitter(height = 0, size = point.size.use) +
    xlab("Cell Type") +
    nogrid +
    ggtitle(feature) +
    theme(plot.title = element_text(size = size.title.use, face = "bold"))
  if (y.log) {
    plot <- plot + scale_y_log10()
  } else {
    plot <- plot + ylim(min(data[, feature]), y.max)
  }
  if (feature %in% gene.names) {
    if (y.log) {
      plot <- plot + ylab(label = "Log Expression level")
    } else {
      plot <- plot + ylab(label = "Expression level")
    }
  } else {
    plot <- plot + ylab(label = "")
  }
  if (! is.null(x = cols.use)) {
    plot <- plot + scale_fill_manual(values = cols.use)
  }
  if (x.lab.rot) {
    plot <- plot + theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
  }
  if (y.lab.rot) {
    plot <- plot + theme(axis.text.x = element_text(angle = 90))
  }
  if (remove.legend) {
    plot <- plot + theme(legend.position = "none")
  }
  return(plot)
}
