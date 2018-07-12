#' Seurat Themes
#'
#' Various themes to be applied to ggplot2-based plots
#'
#' @param ... Extra parameters to be passed to \code{theme}
#'
#' @return A ggplot2 theme object
#'
#' @importFrom ggplot2 theme element_rect element_text element_line margin
#' @export
#'
#' @rdname SeuratThemes
#' @seealso \code{\link{ggplot2::theme}}
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
    axis.line = blank,
    axis.text = blank,
    axis.ticks = blank,
    axis.title = blank,
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
