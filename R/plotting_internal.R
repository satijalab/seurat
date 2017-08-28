# Create a scatterplot with data from a ggplot2 scatterplot
#
# @param plot.data The original ggplot2 scatterplot data
# This is taken from ggplot2::ggplot_build
# @param dark.theme Plot using a dark theme
# @param smooth Use a smooth scatterplot instead of a standard scatterplot
# @param ... Extra parameters passed to graphics::plot or graphics::smoothScatter
#
#' @importFrom graphics axis
#
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

# Convert a ggplot2 scatterplot to base R graphics
#
# @param plot A ggplot2 scatterplot
# @param do.plot Create the plot with base R graphics
# @param ... Extra parameters passed to PlotBuild
#
# @return A dataframe with the data that created the ggplot2 scatterplot
#
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
#
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

# Plot a single feature
#
# @param data.use The data regarding the feature
# @param feature The feature to plot
# @param data.plot The data to be plotted
# @param pt.size Size of each point
# @param pch.use Shape of each point
# @param cols.use Colors to plot
# @param dim.codes Codes for the dimensions to plot in
# @param min.cutoff Minimum cutoff for data
# @param max.cutoff Maximum cutoff for data
# @param no.axes Remove axes from plot
# @param no.legend Remove legend from plot
# @param dark.theme Plot in dark theme
#
# @return A ggplot2 scatterplot
#
#' @importFrom stats na.omit
#' @importFrom utils globalVariables
#
globalVariables(names = c('x', 'y', 'gene'), package = 'Seurat', add = TRUE)
SingleFeaturePlot <- function(
  data.use,
  feature,
  data.plot,
  pt.size,
  pch.use,
  cols.use,
  dim.codes,
  min.cutoff,
  max.cutoff,
  no.axes,
  no.legend,
  dark.theme
) {
  data.gene <- na.omit(object = data.frame(data.use[feature, ]))
  #   Check for quantiles
  min.cutoff <- SetQuantile(cutoff = min.cutoff, data = data.gene)
  max.cutoff <- SetQuantile(cutoff = max.cutoff, data = data.gene)
  #   Mask any values below the minimum and above the maximum values
  data.gene <- sapply(
    X = data.gene,
    FUN = function(x) {
      return(ifelse(test = x < min.cutoff, yes = min.cutoff, no = x))
    }
  )
  data.gene <- sapply(
    X = data.gene,
    FUN = function(x) {
      return(ifelse(test = x > max.cutoff, yes = max.cutoff, no = x))
    }
  )
  data.plot$gene <- data.gene
  #   Stuff for break points
  if (length(x = cols.use) == 1) {
    brewer.gran <- brewer.pal.info[cols.use, ]$maxcolors
  } else {
    brewer.gran <- length(x = cols.use)
  }
  #   Cut points
  if (all(data.gene == 0)) {
    data.cut <- 0
  } else {
    data.cut <- as.numeric(x = as.factor(x = cut(
      x = as.numeric(x = data.gene),
      breaks = brewer.gran
    )))
  }
  data.plot$col <- as.factor(x = data.cut)
  #   Start plotting
  p <- ggplot(data = data.plot, mapping = aes(x = x, y = y))
  if (brewer.gran != 2) {
    if (length(x = cols.use) == 1) {
      p <- p + geom_point(
        mapping = aes(color = col),
        size = pt.size,
        shape = pch.use
      ) + scale_color_brewer(palette = cols.use)
    } else {
      p <- p + geom_point(
        mapping = aes(color = col),
        size = pt.size,
        shape = pch.use
      ) + scale_color_manual(values = cols.use)
    }
  } else {
    if (all(data.plot$gene == data.plot$gene[1])) {
      warning(paste0("All cells have the same value of ", feature, "."))
      p <- p + geom_point(color = cols.use[1], size = pt.size, shape = pch.use)
    } else {
      p <- p + geom_point(
        mapping = aes(color = gene),
        size = pt.size,
        shape = pch.use
      ) + scale_color_gradientn(
        colors = cols.use,
        guide = guide_colorbar(title = feature)
      )
    }
  }
  if (no.axes) {
    p <- p + labs(title = feature, x ="", y="") + theme(
      axis.line = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank()
    )
  } else {
    p <- p + labs(title = feature, x = dim.codes[1], y = dim.codes[2])
  }
  if (no.legend) {
    p <- p + theme(legend.position = 'none')
  }
  if (dark.theme) {
    p <- p + DarkTheme()
  }
  return(p)
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
  rgb.vals <- sapply(X = colors, FUN = grDevices::col2rgb)
  if (nrow(x = rgb.vals) != 3) {
    rgb.vals <- t(x = rgb.vals)
  }
  #   Blend together using the additive method
  #   Basically, resulting colors are the mean of the component colors
  blend <- apply(
    X = rgb.vals,
    MARGIN = 1,
    FUN = mean
  )
  #   If we're returning RGB values, convert to matrix, just like col2rgb
  #   Otherwise, return as hexadecimal; can be used directly for plotting
  if (as.rgb) {
    result <- matrix(
      data = blend,
      nrow = 3,
      dimnames = list(c('red', 'green', 'blue'), 'blend')
    )
  } else {
    result <- grDevices::rgb(
      matrix(data = blend, ncol = 3),
      alpha = alpha.value,
      maxColorValue = 255
    )
  }
  return(result)
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

# No Grid
#
# Remove the grid lines from a ggplot2 plot
#
# @param ... Extra parameters to be passed to theme()
# @import ggplot2
# @return A ggplot2 theme object
# @seealso \code{theme}
# @import ggplot2
# @export
#
NoGrid <- function(...) {
  no.grid <- theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    ...
  )
  return(no.grid)
}

# Reset Par
#
# Reset the graphing space to
# mfrow = c(1, 1)
#
# @param ... Extra parameters for par
#
ResetPar <- function(...) {
  par(mfrow = c(1, 1), ...)
}

# Plot a single feature on a violin plot
#
# @param feature Feature to plot
# @param data Data to plot
# @param cell.ident Idents to use
# @param do.sort Sort identity classes (on the x-axis) by the average
# expression of the attribute being potted
# @param y.max Maximum Y value to plot
# @param size.x.use X axis title font size
# @param size.y.use Y axis title font size
# @param size.title.use Main title font size
# @param adjust.use Adjust parameter for geom_violin
# @param point.size.use Point size for geom_violin
# @param cols.use Colors to use for plotting
# @param gene.names
# @param y.log plot Y axis on log scale
# @param x.lab.rot Rotate x-axis labels
# @param y.lab.rot Rotate y-axis labels
# @param legend.position Position the legend for the plot
# @param remove.legend Remove the legend from the plot
#
# @return A ggplot-based violin plot
#
#' @importFrom stats rnorm
#' @importFrom utils globalVariables
#
globalVariables(names = 'ident', package = 'Seurat', add = TRUE)
SingleVlnPlot <- function(
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
  feature.name <- colnames(data)
  colnames(data) <- "feature"
  feature <- "feature"
  set.seed(seed = 42)
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
  data[, feature] <- data[, feature] + noise
  y.max <- SetIfNull(x = y.max, default = max(data[, feature]))
  plot <- ggplot(
    data = data,
    mapping = aes(
      x = factor(x = ident),
      y = feature
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
    xlab("Identity") +
    NoGrid() +
    ggtitle(feature) +
    theme(plot.title = element_text(size = size.title.use, face = "bold"))
  plot <- plot + ggtitle(feature.name)
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

# Plot a single feature on a joy plot
#
# @param feature Feature to plot
# @param data Data to plot
# @param cell.ident Idents to use
# @param do.sort Sort identity classes (on the x-axis) by the average
# expression of the attribute being potted
# @param y.max Maximum Y value to plot
# @param size.x.use X axis title font size
# @param size.y.use Y axis title font size
# @param size.title.use Main title font size
# @param cols.use Colors to use for plotting
# @param gene.names
# @param y.log plot Y axis on log scale
# @param x.lab.rot Rotate x-axis labels
# @param y.lab.rot Rotate y-axis labels
# @param legend.position Position the legend for the plot
# @param remove.legend Remove the legend from the plot
#
# @return A ggplot-based violin plot
#
#' @importFrom stats rnorm
#' @importFrom utils globalVariables
#' @importFrom ggjoy geom_joy theme_joy
#
globalVariables(names = 'ident', package = 'Seurat', add = TRUE)
SingleJoyPlot <- function(
  feature,
  data,
  cell.ident,
  do.sort,
  y.max,
  size.x.use,
  size.y.use,
  size.title.use,
  cols.use,
  gene.names,
  y.log,
  x.lab.rot,
  y.lab.rot,
  legend.position,
  remove.legend
) {
  set.seed(seed = 42)
  feature.name <- colnames(data)
  colnames(data) <- "feature"
  feature <- "feature"
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
  data[, feature] <- data[, feature] + noise
  y.max <- SetIfNull(x = y.max, default = max(data[, feature]))
  plot <- ggplot(
    data = data,
    mapping = aes(
        x = feature,
        y = factor(ident)
    )
  ) +
    geom_joy(scale = 4, mapping = aes(fill = factor(x = ident))) + theme_joy() +
  scale_y_discrete(expand = c(0.01, 0)) +   # will generally have to set the `expand` option
  scale_x_continuous(expand = c(0, 0))      # for both axes to remove unneeded padding
  plot <- plot+theme(
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
    #geom_jitter(height = 0, size = point.size.use) +
    ylab("Identity") +
    NoGrid() +
    ggtitle(feature) +
    theme(plot.title = element_text(size = size.title.use, face = "bold"))
  plot <- plot + ggtitle(feature.name)
  if (y.log) {
    plot <- plot + scale_x_log10()
  } else {
    #plot <- plot + xlim(min(data[, feature]), y.max)
  }
  if (feature %in% gene.names) {
    if (y.log) {
      plot <- plot + xlab(label = "Log Expression level")
    } else {
      plot <- plot + xlab(label = "Expression level")
    }
  } else {
    plot <- plot + xlab(label = "")
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


#remove legend title
no.legend.title <- theme(legend.title = element_blank())

#set legend text
SetLegendTextGG <- function(x = 12, y = "bold") {
  return(theme(legend.text = element_text(size = x, face = y)))
}

#set legend point size
SetLegendPointsGG <- function(x = 6) {
  return(guides(colour = guide_legend(override.aes = list(size = x))))
}

#set x axis features
SetXAxisGG <- function(x = 16, y = "#990000", z = "bold", x2 = 12) {
  return(theme(
    axis.title.x = element_text(face = z, colour = y, size = x),
    axis.text.x = element_text(angle = 90, vjust = 0.5, size = x2)
  ))
}

#set y axis features
SetYAxisGG <- function(x = 16, y = "#990000", z = "bold", x2 = 12) {
  return(theme(
    axis.title.y = element_text(face = z, colour = y, size = x),
    axis.text.y = element_text(angle = 90, vjust = 0.5, size = x2)
  ))
}

#heatmap.2, but does not draw a key.
#unclear if this is necessary, but valuable to have the function coded in for modifications
#
#' @importFrom graphics axis mtext rect abline text title hist
#' @importFrom stats median order.dendrogram as.dendrogram reorder
#
heatmap2NoKey <- function (
  x,
  Rowv = TRUE,
  Colv = if (symm) "Rowv" else TRUE,
  distfun = dist,
  hclustfun = hclust,
  dendrogram = c("both", "row", "column", "none"),
  symm = FALSE,
  scale = c("none", "row", "column"),
  na.rm = TRUE,
  revC = identical(x = Colv, y = "Rowv"),
  add.expr,
  breaks,
  symbreaks = min(x < 0, na.rm = TRUE) || scale != "none",
  col = "heat.colors",
  colsep,
  rowsep,
  sepcolor = "white",
  sepwidth = c(0.05, 0.05),
  cellnote,
  notecex = 1,
  notecol = "cyan",
  na.color = par("bg"),
  trace = c("column", "row", "both", "none"),
  tracecol = "cyan",
  hline = median(breaks),
  vline = median(breaks),
  linecol = tracecol,
  margins = c(5, 5),
  ColSideColors,
  RowSideColors,
  cexRow = 0.2 + 1 / log10(x = nr),
  cexCol = 0.2 + 1 / log10(x = nc),
  labRow = NULL,
  labCol = NULL,
  key = TRUE,
  keysize = 1.5,
  density.info = c("histogram", "density", "none"),
  denscol = tracecol,
  symkey = min(x < 0, na.rm = TRUE) || symbreaks,
  densadj = 0.25,
  main = NULL,
  xlab = NULL,
  ylab = NULL,
  lmat = NULL,
  lhei = NULL,
  axRowCol="black",
  lwid = NULL,
  dimTitle = NULL,
  pc_title = NULL,
  ...
) {
  scale01 <- function(x, low = min(x), high = max(x)) {
    return((x - low) / (high - low))
  }
  retval <- list()
  scale <- if (symm && missing(x = scale)) {
    "none"
  } else {
    match.arg(arg = scale)
  }
  dendrogram <- match.arg(arg = dendrogram)
  trace <- match.arg(arg = trace)
  density.info <- match.arg(density.info)
  if (length(x = col) == 1 && is.character(x = col)) {
    col <- get(col, mode = "function")
  }
  if (! missing(x = breaks) && (scale != "none")) {
    warning(
      "Using scale=\"row\" or scale=\"column\" when breaks are",
      "specified can produce unpredictable results.",
      "Please consider using only one or the other."
    )
  }
  if (is.null(x = Rowv) || is.na(x = Rowv)) {
    Rowv <- FALSE
  }
  if (is.null(x = Colv) || is.na(x = Colv)) {
    Colv <- FALSE
  } else if (Colv == "Rowv" && !isTRUE(x = Rowv)) {
    Colv <- FALSE
  }
  if (length(x = di <- dim(x = x)) != 2 || !is.numeric(x = x)) {
    stop("`x' must be a numeric matrix")
  }
  nr <- di[1]
  nc <- di[2]
  if (nr <= 1 || nc <= 1) {
    stop("`x' must have at least 2 rows and 2 columns")
  }
  if (! is.numeric(x = margins) || length(x = margins) != 2) {
    stop("`margins' must be a numeric vector of length 2")
  }
  if (missing(x = cellnote)) {
    cellnote <- matrix(data = "", ncol = ncol(x = x), nrow = nrow(x = x))
  }
  if (! inherits(x = Rowv, what = "dendrogram")) {
    if (((! isTRUE(x = Rowv)) || (is.null(x = Rowv))) &&
        (dendrogram %in% c("both", "row"))) {
      if (is.logical(x = Colv) && (Colv)) {
        dendrogram <- "column"
      } else {
        dedrogram <- "none"
      }
    }
  }
  if (! inherits(x = Colv, what = "dendrogram")) {
    if (((!isTRUE(x = Colv)) || (is.null(x = Colv))) &&
        (dendrogram %in% c("both", "column"))) {
      if (is.logical(x = Rowv) && (Rowv)) {
        dendrogram <- "row"
      } else {
        dendrogram <- "none"
      }
    }
  }
  if (inherits(x = Rowv, what = "dendrogram")) {
    ddr <- Rowv
    rowInd <- order.dendrogram(x = ddr)
  } else if (is.integer(x = Rowv)) {
    hcr <- hclustfun(distfun(x))
    ddr <- as.dendrogram(object = hcr)
    ddr <- reorder(x = ddr, X = Rowv)
    rowInd <- order.dendrogram(x = ddr)
    if (nr != length(x = rowInd)) {
      stop("row dendrogram ordering gave index of wrong length")
    }
  } else if (isTRUE(x = Rowv)) {
    Rowv <- rowMeans(x = x, na.rm = na.rm)
    hcr <- hclustfun(distfun(x))
    ddr <- as.dendrogram(object = hcr)
    ddr <- reorder(x = ddr, X = Rowv)
    rowInd <- order.dendrogram(x = ddr)
    if (nr != length(x = rowInd)) {
      stop("row dendrogram ordering gave index of wrong length")
    }
  } else {
    rowInd <- nr:1
  }
  if (inherits(x = Colv, what = "dendrogram")) {
    ddc <- Colv
    colInd <- order.dendrogram(x = ddc)
  } else if (identical(x = Colv, y = "Rowv")) {
    if (nr != nc) {
      stop("Colv = \"Rowv\" but nrow(x) != ncol(x)")
    }
    if (exists(x = "ddr")) {
      ddc <- ddr
      colInd <- order.dendrogram(x = ddc)
    } else {
      colInd <- rowInd
    }
  } else if (is.integer(x = Colv)) {
    hcc <- hclustfun(distfun(if (symm) {
      x
    } else {
      t(x = x)
    }))
    ddc <- as.dendrogram(object = hcc)
    ddc <- reorder(x = ddc, X = Colv)
    colInd <- order.dendrogram(x = ddc)
    if (nc != length(x = colInd)) {
      stop("column dendrogram ordering gave index of wrong length")
    }
  } else if (isTRUE(x = Colv)) {
    Colv <- colMeans(x = x, na.rm = na.rm)
    hcc <- hclustfun(distfun(if (symm) {
      x
    } else {
      t(x = x)
    }))
    ddc <- as.dendrogram(object = hcc)
    ddc <- reorder(x = ddc, X = Colv)
    colInd <- order.dendrogram(x = ddc)
    if (nc != length(x = colInd)) {
      stop("column dendrogram ordering gave index of wrong length")
    }
  } else {
    colInd <- 1:nc
  }
  retval$rowInd <- rowInd
  retval$colInd <- colInd
  retval$call <- match.call()
  x <- x[rowInd, colInd]
  x.unscaled <- x
  cellnote <- cellnote[rowInd, colInd]
  if (is.null(x = labRow)) {
    labRow <- if (is.null(rownames(x = x))) {
      (1:nr)[rowInd]
    } else {
      rownames(x = x)
    }
  } else {
    labRow <- labRow[rowInd]
  }
  if (is.null(x = labCol)) {
    labCol <- if (is.null(x = colnames(x = x))) {
      (1:nc)[colInd]
    } else {
      colnames(x = x)
    }
  } else {
    labCol <- labCol[colInd]
  }
  if (scale == "row") {
    retval$rowMeans <- rm <- rowMeans(x = x, na.rm = na.rm)
    x <- sweep(x = x, MARGIN = 1, STATS = rm)
    retval$rowSDs <- sx <- apply(X = x, MARGIN = 1, FUN = sd, na.rm = na.rm)
    x <- sweep(x = x, MARGIN = 1, STATS = sx, FUN = "/")
  } else if (scale == "column") {
    retval$colMeans <- rm <- colMeans(x = x, na.rm = na.rm)
    x <- sweep(x = x, MARGIN = 2, STATS = rm)
    retval$colSDs <- sx <- apply(X = x, MARGIN = 2, FUN = sd, na.rm = na.rm)
    x <- sweep(x = x, MARGIN = 2, STATS = sx, FUN = "/")
  }
  if (missing(x = breaks) || is.null(x = breaks) || length(x = breaks) <1) {
    if (missing(x = col) || is.function(x = col)) {
      breaks <- 16
    } else {
      breaks <- length(x = col) + 1
    }
  }
  if (length(x = breaks) == 1) {
    if (! symbreaks) {
      breaks <- seq(
        from = min(x, na.rm = na.rm),
        to = max(x, na.rm = na.rm),
        length = breaks
      )
    } else {
      extreme <- max(abs(x = x), na.rm = TRUE)
      breaks <- seq(from = -extreme, to = extreme, length = breaks)
    }
  }
  nbr <- length(x = breaks)
  ncol <- length(x = breaks) - 1
  if (class(x = col) == "function") {
    col <- col(x = ncol)
  }
  min.breaks <- min(breaks)
  max.breaks <- max(breaks)
  x[x < min.breaks] <- min.breaks
  x[x > max.breaks] <- max.breaks
  if (! missing(x = RowSideColors)) {
    par(mar = c(margins[1], 0, 0, 0.5))
    image(x = rbind(1:nr), col = RowSideColors[rowInd], axes = FALSE)
  }
  if (! missing(x = ColSideColors)) {
    par(mar = c(0.5, 0, 0, margins[2]))
    image(x = cbind(1:nc), col = ColSideColors[colInd], axes = FALSE)
  }
  oldMar <- par()$mar
  if (labCol[1] == "") {
    par(mar = c(margins[1]-3, margins[2]-2, margins[1]-3, margins[2]))
  } else {
    par(mar = c(margins[1], margins[2], margins[1], margins[2]))
  }
  x <- t(x = x)
  cellnote <- t(x = cellnote)
  if (revC) {
    iy <- nr:1
    if (exists(x = "ddr")) {
      ddr <- rev(x = ddr)
    }
    x <- x[, iy]
    cellnote <- cellnote[, iy]
  } else {
    iy <- 1:nr
  }
  # add pc number as title if plotting pc heatmaps
  if(is.null(x = dimTitle)) {
    dimTitle <- ""
  }
  if (is.null(x = pc_title)) {
    pc_title <- ''
  }
  #print(dimTitle)
  image(
    x = 1:nc,
    y = 1:nr,
    z = x,
    xlim = 0.5 + c(0, nc),
    ylim = 0.5 + c(0, nr),
    axes = FALSE,
    xlab = "",
    ylab = "",
    main = dimTitle,
    col = col,
    breaks = breaks,
    ...
  )
  retval$carpet <- x
  if (exists(x = "ddr")) {
    retval$rowDendrogram <- ddr
  }
  if (exists(x = "ddc")) {
    retval$colDendrogram <- ddc
  }
  retval$breaks <- breaks
  retval$col <- col
  if (any(is.na(x = x))) {
    mmat <- ifelse(test = is.na(x = x), yes = 1, no = NA)
    image(
      x = 1:nc,
      y = 1:nr,
      z = mmat,
      axes = FALSE,
      xlab = "",
      ylab = "",
      main = pc_title,
      col = na.color,
      add = TRUE
    )
  }
  axis(
    side = 1,
    at = 1:nc,
    labels = labCol,
    las = 2,
    line = -0.5,
    tick = 0,
    cex.axis = cexCol
  )
  if (! is.null(x = xlab)) {
    mtext(text = xlab, side = 1, line = margins[1] - 1.25)
  }
  axis(
    side = 4,
    at = iy,
    labels = labRow,
    las = 2,
    line = -0.5,
    tick = 0,
    cex.axis = cexRow,
    col = axRowCol
  )
  if (! is.null(x = ylab)) {
    mtext(text = ylab, side = 4, line = margins[2] - 1.25)
  }
  if (! missing(x = add.expr)) {
    eval(expr = substitute(expr = add.expr))
  }
  if (! missing(x = colsep)) {
    for (csep in colsep) {
      rect(
        xleft = csep + 0.5,
        ybottom = rep(x = 0, length(x = csep)),
        xright = csep + 0.5 + sepwidth[1],
        ytop = rep(x = ncol(x = x) + 1, csep),
        lty = 1,
        lwd = 1,
        col = sepcolor,
        border = sepcolor
      )
    }
  }
  if (! missing(x = rowsep)) {
    for (rsep in rowsep) {
      rect(
        xleft = 0,
        ybottom = (ncol(x = x) + 1 - rsep) - 0.5,
        xright = nrow(x = x) + 1,
        ytop = (ncol(x) + 1 - rsep) - 0.5 - sepwidth[2],
        lty = 1,
        lwd = 1,
        col = sepcolor,
        border = sepcolor
      )
    }
  }
  min.scale <- min(breaks)
  max.scale <- max(breaks)
  x.scaled <- scale01(x = t(x = x), low = min.scale, high = max.scale)
  if (trace %in% c("both", "column")) {
    retval$vline <- vline
    vline.vals <- scale01(x = vline, low = min.scale, high = max.scale)
    for (i in colInd) {
      if (! is.null(x = vline)) {
        abline(
          v = i - 0.5 + vline.vals,
          col = linecol,
          lty = 2
        )
      }
      xv <- rep(x = i, nrow(x = x.scaled)) + x.scaled[, i] - 0.5
      xv <- c(xv[1], xv)
      yv <- 1:length(x = xv) - 0.5
      ##lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
    }
  }
  if (trace %in% c("both", "row")) {
    retval$hline <- hline
    hline.vals <- scale01(x = hline, low = min.scale, high = max.scale)
    for (i in rowInd) {
      if (! is.null(x = hline)) {
        abline(h = i + hline, col = linecol, lty = 2)
      }
      yv <- rep(x = i, ncol(x = x.scaled)) + x.scaled[i, ] - 0.5
      yv <- rev(x = c(yv[1], yv))
      xv <- length(x = yv):1 - 0.5
      ##lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
    }
  }
  if (! missing(x = cellnote)) {
    text(
      x = c(row(x = cellnote)),
      y = c(col(x = cellnote)),
      labels = c(cellnote),
      col = notecol,
      cex = notecex
    )
  }
  #par(mar = c(margins[1], 0, 0, 0))
  if (dendrogram %in% c("both", "row")) {
    ##plot(ddr, horiz = TRUE, axes = FALSE, yaxs = "i", leaflab = "none")
  }
  ##else plot.new()
  #par(mar = c(0, 0, if (!is.null(main)) 5 else 0, margins[2]))
  if (dendrogram %in% c("both", "column")) {
    ##plot(ddc, axes = FALSE, xaxs = "i", leaflab = "none")
  }
  ##else plot.new()
  key <- FALSE
  if (! is.null(x = main))
    title(main = main, cex.main = 1.5 * par()[["cex.main"]])
  if (key) {
    par(mar = c(5, 4, 2, 1), cex = 0.75)
    tmpbreaks <- breaks
    if (symkey) {
      max.raw <- max(abs(x = c(x, breaks)), na.rm = TRUE)
      min.raw <- -max.raw
      tmpbreaks[1] <- -max(abs(x = x), na.rm = TRUE)
      tmpbreaks[length(x = tmpbreaks)] <- max(abs(x = x), na.rm = TRUE)
    } else {
      min.raw <- min(x, na.rm = TRUE)
      max.raw <- max(x, na.rm = TRUE)
    }
    z <- seq(from = min.raw, to = max.raw, length = length(x = col))
    #image(z = matrix(z, ncol = 1), col = col, breaks = tmpbreaks,
    #      xaxt = "n", yaxt = "n")
    par(usr = c(0, 1, 0, 1))
    lv <- pretty(x = breaks)
    xv <- scale01(x = as.numeric(x = lv), low = min.raw, high = max.raw)
    axis(side = 1, at = xv, labels = lv)
    if (scale == "row") {
      mtext(side = 1, "Row Z-Score", line = 2)
    }
    else if (scale == "column") {
      mtext(side = 1, "Column Z-Score", line = 2)
    }
    else {
      mtext(side = 1, "Value", line = 2)
    }
    if (density.info == "density") {
      dens <- density(x = x, adjust = densadj, na.rm = TRUE)
      omit <- dens$x < min(breaks) | dens$x > max(breaks)
      dens$x <- dens$x[-omit]
      dens$y <- dens$y[-omit]
      dens$x <- scale01(x = dens$x, low = min.raw, high = max.raw)
      lines(
        x = dens$x,
        y = dens$y / max(dens$y) * 0.95,
        col = denscol,
        lwd = 1
      )
      axis(
        side = 2,
        at = pretty(x = dens$y) / max(dens$y) * 0.95,
        pretty(dens$y)
      )
      title(main = "Color Key\nand Density Plot")
      par(cex = 0.5)
      mtext(side = 2, "Density", line = 2)
    } else if (density.info == "histogram") {
      h <- hist(x = x, plot = FALSE, breaks = breaks)
      hx <- scale01(x = breaks, low = min.raw, high = max.raw)
      hy <- c(h$counts, h$counts[length(h$counts)])
      lines(
        x = hx,
        y = hy / max(hy) * 0.95,
        lwd = 1,
        type = "s",
        col = denscol
      )
      axis(
        side = 2,
        at = pretty(x = hy) / max(hy) * 0.95,
        pretty(x = hy)
      )
      title(main = "Color Key\nand Histogram")
      par(cex = 0.5)
      mtext(side = 2, "Count", line = 2)
    } else {
      title(main = "Color Key")
    }
  }
  ##else plot.new()
  retval$colorTable <- data.frame(
    low = retval$breaks[-length(x = retval$breaks)],
    high = retval$breaks[-1], color = retval$col
  )
  invisible(x = retval)
  par(mar = oldMar)
}

# Documentation
###############
PlotDim <- function(
  ndim,
  object,
  reduction.type,
  use.scaled,
  use.full,
  cells.use,
  num.genes,
  group.by,
  disp.min,
  disp.max,
  col.low,
  col.mid,
  col.high,
  slim.col.label,
  do.balanced,
  remove.key,
  cex.col,
  cex.row,
  group.label.loc,
  group.label.rot,
  group.cex,
  group.spacing
) {
  if (is.numeric(x = (cells.use))) {
    cells.use <- DimTopCells(
      object = object,
      dim.use = ndim,
      reduction.type = reduction.type,
      num.cells = cells.use,
      do.balanced = do.balanced
    )
  } else {
    cells.use <- SetIfNull(x = cells.use, default = object@cell.names)
  }
  genes.use <- rev(x = DimTopGenes(
    object = object,
    dim.use = ndim,
    reduction.type = reduction.type,
    num.genes = num.genes,
    use.full = use.full,
    do.balanced = do.balanced
  ))
  dim.scores <- GetDimReduction(
    object = object,
    reduction.type = reduction.type,
    slot = "cell.embeddings"
  )
  dim.key <- GetDimReduction(
    object = object,
    reduction.type = reduction.type,
    slot = "key"
  )
  cells.ordered <- cells.use[order(dim.scores[cells.use, paste0(dim.key, ndim)])]
  if (! use.scaled) {
    data.use <- as.matrix(object@data[genes.use, cells.ordered])
  } else {
    data.use <- object@scale.data[genes.use, cells.ordered]
    data.use <- MinMax(data = data.use, min = disp.min, max = disp.max)
  }
  return(DoHeatmap(
    object = object,
    data.use = data.use,
    cells.use = cells.use,
    genes.use = genes.use,
    group.by = group.by,
    disp.min = disp.min,
    disp.max = disp.max,
    col.low = col.low,
    col.mid = col.mid,
    col.high = col.high,
    slim.col.label = slim.col.label,
    remove.key = remove.key,
    cex.col = cex.col,
    cex.row = cex.row,
    group.label.loc = group.label.loc,
    group.label.rot = group.label.rot,
    group.cex = group.cex,
    group.spacing = group.spacing,
    title = paste0(dim.key, ndim),
    do.plot = FALSE
  ))
}

# Posterior Plot
#
# @param object A Seurat object
# @param name Spatial code
#
# @seealso \code{SubsetColumn}
# @seealso \code{VlnPlot}
#
PosteriorPlot <- function(object, name) {
  post.names <- colnames(x = SubsetColumn(data = object@spatial@mix.probs, code = name))
  VlnPlot(
    object = object,
    features.plot = post.names,
    inc.first = TRUE,
    inc.final = TRUE,
    by.k = TRUE
  )
}
