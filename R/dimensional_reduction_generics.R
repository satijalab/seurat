#' Run Principal Component Analysis on gene expression using IRLBA
#'
#' Run a PCA dimensionality reduction. For details about stored PCA calculation
#' parameters, see \code{PrintPCAParams}.
#'
#' @param object Seurat object
#' @param assay.use Name of Assay PCA is being run on
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

#' Run t-distributed Stochastic Neighbor Embedding
#'
#' Run t-SNE dimensionality reduction on selected features. Has the option of
#' running in a reduced dimensional space (i.e. spectral tSNE, recommended),
#' or running based on a set of genes. For details about stored TSNE calculation
#' parameters, see \code{PrintTSNEParams}.
#'
#' @param object Seurat object
#' @param reduction.use Which dimensional reduction (e.g. PCA, ICA) to use for
#' the tSNE. Default is PCA
#' @param cells.use Which cells to analyze (default, all cells)
#' @param seed.use Random seed for the t-SNE
#' @param tsne.method Select the method to use to compute the tSNE. Available
#' methods are:
#' \itemize{
#' \item{Rtsne: }{Use the Rtsne package Barnes-Hut implementation of tSNE (default)}
#' \item{tsne: }{standard tsne - not recommended for large datasets}
#' \item{FIt-SNE: }{Use the FFT-accelerated Interpolation-based t-SNE. Based on
#' Kluger Lab code found here: https://github.com/KlugerLab/FIt-SNE}
#' }
#' @param add.iter If an existing tSNE has already been computed, uses the
#' current tSNE to seed the algorithm and then adds additional iterations on top
#' of this
#' @param dim.embed The dimensional space of the resulting tSNE embedding
#' (default is 2). For example, set to 3 for a 3d tSNE
#' @param \dots Additional arguments to the tSNE call. Most commonly used is
#' perplexity (expected number of neighbors default is 30)
#' @param reduction.key dimensional reduction key, specifies the string before the number for the dimension names. tSNE_ by default
#'
#' @rdname RunTSNE
#' @export RunTSNE
#'
RunTSNE <- function(
  object,
  reduction.use,
  cells.use,
  seed.use,
  tsne.method,
  add.iter,
  dim.embed,
  reduction.key,
  ...
) {
  UseMethod(generic = 'RunTSNE', object = object)
}
