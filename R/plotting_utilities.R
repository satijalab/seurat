#' Dark Theme
#'
#' Add a dark theme to ggplot objects
#'
#' @param ... Extra parameters to be passed to theme()
#' @import ggplot2
#' @return A ggplot2 theme object
#' @seealso \code{theme}
#' @import ggplot2
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
#'
#' @export
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

#' A black and white color palette
#'
#' @param ... Extra parameters to CustomPalette
#'
#' @return A color palette
#'
#' @seealso \code{CustomPalette}
#'
#' @export
#'
BlackAndWhite <- function(...) {
  return(CustomPalette(low = "white", high="black", ...))
}

#' A purple and yellow color palette
#'
#' @param ... Extra parameters to CustomPalette
#'
#' @return A color palette
#'
#' @seealso \code{CustomPalette}
#'
#' @export
#'
PurpleAndYellow <- function(...) {
  return(CustomPalette(low = "magenta", high = "yellow", mid = "black", ...))
}
