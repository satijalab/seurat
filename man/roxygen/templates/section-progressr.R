#' @section Progress Updates with \pkg{progressr}:
#' This function uses
#' \href{https://cran.r-project.org/package=progressr}{\pkg{progressr}} to
#' render status updates and progress bars. To enable progress updates, wrap
#' the function call in \code{\link[progressr]{with_progress}} or run
#' \code{\link[progressr:handlers]{handlers(global = TRUE)}} before running
#' this function. For more details about \pkg{progressr}, please read
#' \href{https://progressr.futureverse.org/articles/progressr-01-intro.html}{\code{vignette("progressr-intro")}}
