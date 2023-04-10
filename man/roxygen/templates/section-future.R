#' @section Parallelization with \pkg{future}:
#' This function uses
#' \href{https://cran.r-project.org/package=future}{\pkg{future}} to enable
#' parallelization. Parallelization strategies can be set using
#' \code{\link[future]{plan}}. Common plans include \dQuote{\code{sequential}}
#' for non-parallelized processing or \dQuote{\code{multisession}} for parallel
#' evaluation using multiple \R sessions; for other plans, see the
#' \dQuote{Implemented evaluation strategies} section of
#' \code{\link[future:plan]{?future::plan}}. For a more thorough introduction
#' to \pkg{future}, see
#' \href{https://future.futureverse.org/articles/future-1-overview.html}{\code{vignette("future-1-overview")}}
#'
#' @concept future
