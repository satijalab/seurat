#' Run Principal Component Analysis on gene expression using IRLBA
#'
#' Run a PCA dimensionality reduction. For details about stored PCA calculation
#' parameters, see \code{PrintPCAParams}.
#'
#' @param object Seurat object
#' @param assay.use Name of Assay PCA is being run on
#' @param features.use Features to use as input for PCA. Defaults to variable features
#' @param pcs.compute Total Number of PCs to compute and store (20 by default)
#' @param rev.pca By default computes the PCA on the cell x gene matrix. Setting
#' to true will compute it on gene x cell matrix.
#' @param weight.by.var Weight the cell embeddings by the variance of each PC
#' (weights the gene loadings if rev.pca is TRUE)
#' @param verbose Print the top genes associated with high/low loadings for
#' the PCs
#' @param pcs.print PCs to print genes for
#' @param features.print Number of genes to print for each PC
#' @param reduction.name dimensional reduction name,  pca by default
#' @param reduction.key dimensional reduction key, specifies the string before
#' the number for the dimension names. PC by default
#' @param seed.use Set a random seed. By default, sets the seed to 42. Setting
#' NULL will not set a seed.
#' @param \dots Additional arguments to be passed to IRLBA
#'
#' @importFrom irlba irlba
#' @importFrom methods new
#'
#' @return Returns Seurat object with the PCA calculation stored in the reductions slot
#'
#' @importFrom irlba irlba
#'
#' @export
#'
#' @rdname RunPCA
#' @export RunPCA
#'
RunPCA <- function(
  object,
  features.use,
  pcs.compute,
  rev.pca,
  weight.by.var,
  verbose,
  pcs.print,
  features.print,
  reduction.name,
  reduction.key,
  seed.use,
  ...
) {
  UseMethod(generic = 'RunPCA', object = object)
}
