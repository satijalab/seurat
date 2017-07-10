# Create a scatterplot with data from a ggplot2 scatterplot
#
# @param plot.data The original ggplot2 scatterplot data
# This is taken from ggplot2::ggplot_build
# @param dark.theme Plot using a dark theme
# @param smooth Use a smooth scatterplot instead of a standard scatterplot
# @param ... Extra parameters passed to graphics::plot or graphics::smoothScatter
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
# @return A blended ggplot2 scatterplot
#
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
  data.gene <- na.omit(object = data.frame(data.use[features.plot, ]))
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
  rgb.vals <- sapply(X = colors, FUN = col2rgb)
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
    result <- rgb(
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
# @param cutoff The cutoff to turn into a quantile
# @param data The data to turn find the quantile of
#
# @return The numerical representation of the quantile
#
SetQuantile <- function(cutoff, data) {
  if (grepl(
    pattern = '^q[0-9]{1,2}$',
    x = as.character(x = cutoff),
    perl = TRUE
  )) {
    this.quantile <- as.numeric(x = sub(
      pattern = 'q',
      replacement = '',
      x = as.character(x = cutoff)
    )) / 100
    cutoff <- quantile(x = unlist(x = data), probs = this.quantile)
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
# @seealso \code{\link{theme}}
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
ResetPar <- function() {
  par(mfrow = c(1, 1))
}
