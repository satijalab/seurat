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

#' Get JackStrawData slot
#'
#' Accessor method for JackStrawData
#'
#' @param object JackStrawData object
#' @param slot Name of slot to get
#' @return Returns contents of specified slot

#' @rdname GetJS
#' @export GetJS
#'
GetJS <- function(object, slot, ...) {
  UseMethod(generic = 'GetJS', object = object)
}

#' Set JackStrawData slot
#'
#' Mutator method for JackStrawData
#'
#' @param object JackStrawData object
#' @param slot Name of slot to set
#' @param new.data New data to replace slot with
#' @return Returns JackStrawData object with modified slot
#'
#' @rdname SetJS
#' @export SetJS
#'
SetJS <- function(object, slot, new.data, ...) {
  UseMethod(generic = 'SetJS', object = object)
}

