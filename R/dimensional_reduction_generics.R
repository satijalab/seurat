#' Run Principal Component Analysis on gene expression using IRLBA
#'
#' Run a PCA dimensionality reduction. For details about stored PCA calculation
#' parameters, see \code{PrintPCAParams}.
#'
#' @param object Seurat object
#' @param assay.use Name of Assay PCA is being run on
#' @param compute.dims Total Number of PCs to compute and store (20 by default)
#' @param rev.pca By default computes the PCA on the cell x gene matrix. Setting
#' to true will compute it on gene x cell matrix.
#' @param weight.by.var Weight the cell embeddings by the variance of each PC
#' (weights the gene loadings if rev.pca is TRUE)
#' @param verbose Print the top genes associated with high/low loadings for
#' the PCs
#' @param print.dims PCs to print genes for
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



#' Run UMAP
#'
#' Runs the Uniform Manifold Approximation and Projection (UMAP) dimensional
#' reduction technique. To run, you must first install the umap-learn python
#' package (e.g. via pip install umap-learn). Details on this package can be
#' found here: \url{https://github.com/lmcinnes/umap}. For a more in depth
#' discussion of the mathematics underlying UMAP, see the ArXiv paper here:
#' \url{https://arxiv.org/abs/1802.03426}.
#'
#' @param object Seurat object
#' @param cells.use Which cells to analyze (default, all cells)
#' @param dims.use Which dimensions to use as input features, used only if
#' \code{genes.use} is NULL
#' @param reduction.use Which dimensional reduction (PCA or ICA) to use for the
#' UMAP input. Default is PCA
#' @param genes.use If set, run UMAP on this subset of genes (instead of running on a
#' set of reduced dimensions). Not set (NULL) by default
#' @param assay.use Assay to pull data for when using \code{genes.use}
#' @param max.dim Max dimension to keep from UMAP procedure.
#' @param reduction.name dimensional reduction name, specifies the position in
#' the object$dr list. umap by default
#' @param reduction.key dimensional reduction key, specifies the string before
#' the number for the dimension names. UMAP by default
#' @param n_neighbors This determines the number of neighboring points used in
#' local approximations of manifold structure. Larger values will result in more
#' global structure being preserved at the loss of detailed local structure. In
#' general this parameter should often be in the range 5 to 50.
#' @param min_dist min_dist: This controls how tightly the embedding is allowed
#' compress points together. Larger values ensure embedded points are more
#' evenly distributed, while smaller values allow the algorithm to optimise more
#' accurately with regard to local structure. Sensible values are in the range
#' 0.001 to 0.5.
#' @param metric metric: This determines the choice of metric used to measure
#' distance in the input space. A wide variety of metrics are already coded, and
#' a user defined function can be passed as long as it has been JITd by numba.
#' @param seed.use Set a random seed. By default, sets the seed to 42. Setting
#' NULL will not set a seed.
#' @param ... Additional arguments to the umap
#'
#' @return Returns a Seurat object containing a UMAP representation
#'
#' @references McInnes, L, Healy, J, UMAP: Uniform Manifold Approximation and Projection for Dimension Reduction, ArXiv e-prints 1802.03426, 2018
#'
#' @importFrom reticulate import py_module_available py_set_seed
#'
#' @export
#'
#' @examples
#' \dontrun{
#' pbmc_small
#' # Run UMAP map on first 5 PCs
#' pbmc_small <- RunUMAP(object = pbmc_small, dims.use = 1:5)
#' # Plot results
#' DimPlot(object = pbmc_small, reduction.use = 'umap')
#' }
#'
#' @rdname RunUMAP
#' @export RunUMAP
#'
RunUMAP <- function(
  object,
  cells.use,
  dims.use,
  reduction.use,
  genes.use,
  assay.use,
  max.dim ,
  reduction.name,
  reduction.key,
  n_neighbors ,
  min_dist,
  metric,
  seed.use,
  ...
) {
  UseMethod(generic = 'RunUMAP', object = object)
}

#' Perform Canonical Correlation Analysis
#'
#' Runs a canonical correlation analysis using a diagonal implementation of CCA.
#' For details about stored CCA calculation parameters, see
#' \code{PrintCCAParams}.
#' @param object1 First Seurat object
#' @param object2 Second Seurat object.
#' @param num.cc Number of canonical vectors to calculate
#'
#' @return Returns a combined Seurat object with the CCA results stored.
#'
#' @seealso \code{MergeSeurat}
#'
#' @export
#' @importFrom irlba irlba
#'
#' @examples
#' pbmc_small
#' # As CCA requires two datasets, we will split our test object into two just for this example
#' pbmc1 <- SubsetData(pbmc_small,cells.use = pbmc_small@cell.names[1:40])
#' pbmc2 <- SubsetData(pbmc_small,cells.use = pbmc_small@cell.names[41:80])
#' pbmc1@meta.data$group <- "group1"
#' pbmc2@meta.data$group <- "group2"
#' pbmc_cca <- RunCCA(pbmc1,pbmc2)
#' # Print results
#' PrintDim(pbmc_cca,reduction.type = 'cca')
#'
#' @rdname RunCCA
#' @export RunCCA
#'
RunCCA <- function(
  object1,
  object2,
  num.cc,
  ...
) {
  UseMethod(generic = 'RunCCA', object = object1)
}
