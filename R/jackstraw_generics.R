#' Compute Jackstraw scores significance.
#'
#' Significant PCs should show a p-value distribution that is
#' strongly skewed to the left compared to the null distribution.
#' The p-value for each PC is based on a proportion test comparing the number
#' of genes with a p-value below a particular threshold (score.thresh), compared with the
#' proportion of genes expected under a uniform distribution of p-values.
#'
#' @param object An object
#' @param dims Which dimensions to examine
#' @param score.thresh Threshold to use for the proportion test of PC
#' significance (see Details)
#'
#' @return Returns a Seurat object
#'
#' @author Thanks to Omri Wurtzel for integrating with ggplot
#'
#' @rdname ScoreJackStraw
#' @export ScoreJackStraw
#'
ScoreJackStraw <- function(object, dims, score.thresh, ...) {
  UseMethod(generic = 'ScoreJackStraw', object = object)
}
