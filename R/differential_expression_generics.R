#' Gene expression markers of identity classes
#'
#' Finds markers (differentially expressed genes) for identity classes
#'
#' @param object Seurat object
#' @param features.use Genes to test. Default is to use all genes
#' @param logfc.threshold Limit testing to genes which show, on average, at least
#' X-fold difference (log-scale) between the two groups of cells. Default is 0.25
#' Increasing logfc.threshold speeds up the function, but can miss weaker signals.
#' @param test.use Denotes which test to use. Available options are:
##' \itemize{
##'  \item{"wilcox"} : Identifies differentially expressed genes between two
##'  groups of cells using a Wilcoxon Rank Sum test (default)
##'  \item{"bimod"} : Likelihood-ratio test for single cell gene expression,
##'  (McDavid et al., Bioinformatics, 2013)
##'  \item{"roc"} : Identifies 'markers' of gene expression using ROC analysis.
##'  For each gene, evaluates (using AUC) a classifier built on that gene alone,
##'  to classify between two groups of cells. An AUC value of 1 means that
##'  expression values for this gene alone can perfectly classify the two
##'  groupings (i.e. Each of the cells in cells.1 exhibit a higher level than
##'  each of the cells in cells.2). An AUC value of 0 also means there is perfect
##'  classification, but in the other direction. A value of 0.5 implies that
##'  the gene has no predictive power to classify the two groups. Returns a
##'  'predictive power' (abs(AUC-0.5)) ranked matrix of putative differentially
##'  expressed genes.
##'  \item{"t"} : Identify differentially expressed genes between two groups of
##'  cells using the Student's t-test.
##'  \item{"tobit"} : Identifies differentially expressed genes between two
##'  groups of cells using Tobit models, as proposed in Trapnell et al.,
##'  Nature Biotechnology, 2014
##'  \item{"negbinom"} : Identifies differentially expressed genes between two
##'   groups of cells using a negative binomial generalized linear model.
##'   Use only for UMI-based datasets
##'  \item{"poisson"} : Identifies differentially expressed genes between two
##'   groups of cells using a poisson generalized linear model.
##'   Use only for UMI-based datasets
##'  \item{"MAST} : Identifies differentially expressed genes between two groups
##'  of cells using a hurdle model tailored to scRNA-seq data. Utilizes the MAST
##'  package to run the DE testing.
##'  \item{"DESeq2} : DE based on a model using the negative binomial
##'  distribution (Love et al, Genome Biology, 2014)
##' }
#' @param min.pct  only test genes that are detected in a minimum fraction of
#' min.pct cells in either of the two populations. Meant to speed up the function
#' by not testing genes that are very infrequently expressed. Default is 0.1
#' @param min.diff.pct  only test genes that show a minimum difference in the
#' fraction of detection between the two groups. Set to -Inf by default
#' @param only.pos Only return positive markers (FALSE by default)
#' @param print.bar Print a progress bar once expression testing begins (uses
#' pbapply to do this)
#' @param max.cells.per.ident Down sample each identity class to a max number.
#' Default is no downsampling. Not activated by default (set to Inf)
#' @param random.seed Random seed for downsampling
#' @param latent.vars Variables to test, used only when \code{test.use} is one of
#' 'negbinom', 'poisson', or 'MAST'
#' @param min.cells.features Minimum number of cells expressing the feature in at least one
#' of the two groups, currently only used for poisson and negative binomial tests
#' @param min.cells.group Minimum number of cells in one of the groups
#' @param pseudocount.use Pseudocount to add to averaged expression values when
#' calculating logFC. 1 by default.
#' @param \dots Additional parameters to pass to specific DE functions
#' @seealso \code{\link{MASTDETest}}, and \code{\link{DESeq2DETest}} for more information on these methods
#' @return Matrix containing a ranked list of putative markers, and associated
#' statistics (p-values, ROC score, etc.)
#' @details p-value adjustment is performed using bonferroni correction based on
#' the total number of genes in the dataset. Other correction methods are not
#' recommended, as Seurat pre-filters genes using the arguments above, reducing
#' the number of tests performed. Lastly, as Aaron Lun has pointed out, p-values
#' should be interpreted cautiously, as the genes used for clustering are the
#' same genes tested for differential expression.
#'
#' @references Andrew McDavid, Greg Finak and Masanao Yajima (2017). MAST: Model-based
#' Analysis of Single Cell Transcriptomics. R package version 1.2.1.
#' https://github.com/RGLab/MAST/
#' @import pbapply
#' @importFrom lmtest lrtest
#'
#' @seealso \code{\link{NegBinomDETest}}
#'
#' @export
#'
#' @examples
#' markers <- FindMarkers(object = pbmc_small, ident.1 = 3)
#' head(markers)
#'
#' @rdname FindMarkers
#' @export FindMarkers
#'
FindMarkers <- function(
  object,
  features.use,
  logfc.threshold,
  test.use,
  min.pct,
  min.diff.pct,
  verbose,
  only.pos,
  max.cells.per.ident,
  random.seed,
  latent.vars,
  min.cells.feature,
  min.cells.group,
  pseudocount.use,
  ...
) {
  UseMethod(generic = 'FindMarkers', object = object)
}
