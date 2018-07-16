#' Compute Jackstraw scores significance.
#'
#' Significant PCs should show a p-value distribution that is
#' strongly skewed to the left compared to the null distribution.
#' The p-value for each PC is based on a proportion test comparing the number
#' of genes with a p-value below a particular threshold (score.thresh), compared with the
#' proportion of genes expected under a uniform distribution of p-values.
#'
#' @param object Seurat plot
#' @param PCs Which PCs to examine
#' @param nCol Number of columns
#' @param score.thresh Threshold to use for the proportion test of PC
#' significance (see Details)
#' @param do.plot Show plot
#'
#' @return Returns a Seurat object
#'
#' @author Thanks to Omri Wurtzel for integrating with ggplot
#'
#' @importFrom reshape2 melt
#' @importFrom stats qqplot runif prop.test qunif
#'
#' @export
#' @rdname ScoreJackStraw
#' @export ScoreJackStraw
#'
ScoreJackStraw <- function(object, slot, ...) {
  UseMethod(generic = 'ScoreJackStraw', object = object)
}
