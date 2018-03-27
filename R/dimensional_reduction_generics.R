#' @include seurat.R
NULL

#' Run PCA
#'
#' Run a PCA dimensionality reduction. For details about stored PCA calculation
#' parameters, see \code{PrintPCAParams}.
#'
#' @param object An object to run PCA on
#'
#' @rdname RunPCA
#' @export RunPCA
#'
RunPCA <- function(object, ...) {
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
#' @param dims.use Which dimensions to use as input features
#' @param genes.use If set, run the tSNE on this subset of genes
#' (instead of running on a set of reduced dimensions). Not set (NULL) by default
#' @param seed.use Random seed for the t-SNE
#' @param tsne.method Select the method to use to compute the tSNE. Available
#' methods are:
#' \itemize{
#' \item{Rtsne: }{Use the Rtsne package Barnes-Hut implementation of tSNE (default)}
#' \item{tsne: }{standard tsne - not recommended for large datasets}
#' \item{FIt-SNE: }{Use the FFT-accelerated Interpolation-based t-SNE. Based on
#' Kluger Lab code found here: https://github.com/ChristophH/FIt-SNE}
#' }
#' @param add.iter If an existing tSNE has already been computed, uses the
#' current tSNE to seed the algorithm and then adds additional iterations on top
#' of this
#' @param dim.embed The dimensional space of the resulting tSNE embedding
#' (default is 2). For example, set to 3 for a 3d tSNE
#' @param \dots Additional arguments to the tSNE call. Most commonly used is
#' perplexity (expected number of neighbors default is 30)
#' @param distance.matrix If set, runs tSNE on the given distance matrix
#' instead of data matrix (experimental)
#' @param reduction.name dimensional reduction name, specifies the position in the object$dr list. tsne by default
#' @param reduction.key dimensional reduction key, specifies the string before the number for the dimension names. tSNE_ by default
#'
#' @return Returns a Seurat object with a tSNE embedding in
#' object@@dr$tsne@cell.embeddings
#'
#' @importFrom Rtsne Rtsne
#' @importFrom tsne tsne
#'
#' @rdname RunTSNE
#' @export RunTSNE
#'
#' @examples
#' pbmc_small
#' # Run tSNE on first five PCs, note that for test dataset (only 80 cells)
#' # we can't use default perplexity of 30
#' pbmc_small <- RunTSNE(pbmc_small, reduction.use = "pca", dims.use = 1:5, perplexity=10)
#' # Run tSNE on first five independent components from ICA
#' pbmc_small <- RunICA(pbmc_small,ics.compute=5)
#' pbmc_small <- RunTSNE(pbmc_small, reduction.use = "ica", dims.use = 1:5, perplexity=10)
#' # Plot results
#' TSNEPlot(pbmc_small)
#'
RunTSNE <- function(object, ...) {
  UseMethod(generic = 'RunTSNE', object = object)
}

################################################################################
############################## Utilities Generics ##############################
################################################################################

#' Dimensional Reduction Accessor Function
#'
#' General accessor function for dimensional reduction objects. Pulls slot
#' contents for specified stored dimensional reduction analysis.
#'
#' @param object An object
#' @param reduction.type Type of dimensional reduction to fetch (default is PCA)
#' @param slot Specific information to pull (must be one of the following:
#'  "cell.embeddings", "gene.loadings", "gene.loadings.full", "sdev", "key", "misc")
#'  (replace the '.' with a '_' for loom objects)
#'
#' @return Returns specified slot results from given reduction technique
#'
#' @rdname GetDimReduction
#' @export GetDimReduction
#'
#' @examples
#' pbmc_small
#' # Get the PCA cell embeddings and print the top left corner
#' GetDimReduction(object = pbmc_small, reduction.type = "pca",
#'                 slot = "cell.embeddings")[1:5, 1:5]
#' # Get the standard deviation of each PC
#' GetDimReduction(object = pbmc_small, reduction.type = "pca", slot = "sdev")
#'
GetDimReduction <- function(object, ...) {
  UseMethod(generic = 'GetDimReduction', object = object)
}

#' Dimensional Reduction Mutator Function
#'
#' Set information for specified stored dimensional reduction analysis
#'
#' @param object Seurat object
#' @param reduction.type Type of dimensional reduction to set
#' @param slot Specific information to set (must be one of the following:
#' "cell.embeddings", "gene.loadings", "gene.loadings.full", "sdev", "key",
#' "misc")
#' @param new.data New data to set
#' @return Seurat object with updated slot
#'
#' @rdname SetDimReduction
#' @export SetDimReduction
#'
#' @examples
#' pbmc_small
#' # Simulate adding a new dimensional reduction
#' new.cell.embeddings <- GetCellEmbeddings(object = pbmc_small, reduction.type = "pca")
#' new.gene.loadings <- GetGeneLoadings(object = pbmc_small, reduction.type = "pca")
#' SetDimReduction(
#'   object = pbmc_small,
#'   reduction.type = "new.pca",
#'   slot = "cell.embeddings",
#'   new.data = new.cell.embeddings
#' )
#' SetDimReduction(
#'   object = pbmc_small,
#'   reduction.type = "new.pca",
#'   slot = "gene.loadings",
#'   new.data = new.gene.loadings
#' )
#'
SetDimReduction <- function(object, ...) {
  UseMethod(generic = 'SetDimReduction', object = object)
}
