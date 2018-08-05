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
  ...
) {
  if (is.null(x = num.col)) {
    if (length(x = features.plot) > 9) {
      num.col <- 4
    } else {
      num.col <- min(length(x = features.plot), 3)
    }
  }
  data.use <- FetchData(object = object, vars.fetch = features.plot)
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