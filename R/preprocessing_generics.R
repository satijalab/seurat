#' @include seurat.R
#' @importFrom methods setGeneric
NULL

#' Normalize Assay Data
#'
#' Normalize data for a given assay
#'
#' @param object Seurat object
#' @param assay.type Type of assay to normalize for (default is RNA), but can be
#' changed for multimodal analyses.
#' @param normalization.method Method for normalization. Default is
#' log-normalization (LogNormalize). More methods to be added very shortly.
#' @param scale.factor Sets the scale factor for cell-level normalization
#' @param chunk.size Chunk size to iterate over
#' @param name Basename to store results matrix in 'layers'
#' @param dataset.use Dataset to normalize, defaults to 'matrix'
#' @param display.progress Display progress bar for scaling procedure
#' @param overwrite Overwrite existing dataset with name 'layers/\code{name}'
#'
#' @return Returns object after normalization. Normalized data is stored in data
#' or scale.data slot, depending on the method
#'
#' @rdname NormalizeData
#' @export NormalizeData
#'
#' @examples
#' pbmc_small
#' pmbc_small <- NormalizeData(object = pbmc_small)
#'
setGeneric(
  name = 'NormalizeData',
  def = function(
    object,
    ...
  ) {
    return(standardGeneric(f = 'NormalizeData'))
  }
)

#' Identify variable genes
#'
#' Identifies genes that are outliers on a 'mean variability plot'. First, uses
#' a function to calculate average expression (mean.function) and dispersion (dispersion.function)
#' for each gene. Next, divides genes into num.bin (deafult 20) bins based on
#' their average expression, and calculates z-scores for dispersion within each
#' bin. The purpose of this is to identify variable genes while controlling for
#' the strong relationship between variability and average expression.
#'
#' Exact parameter settings may vary empirically from dataset to dataset, and
#' based on visual inspection of the plot.
#' Setting the y.cutoff parameter to 2 identifies genes that are more than two standard
#' deviations away from the average dispersion within a bin. The default X-axis function
#' is the mean expression level, and for Y-axis it is the log(Variance/mean). All mean/variance
#' calculations are not performed in log-space, but the results are reported in log-space -
#' see relevant functions for exact details.
#'
#' @param object Seurat object
#' @param mean.function Function to compute x-axis value (average expression). Default
#' is to take the mean of the detected (i.e. non-zero) values
#' @param dispersion.function Function to compute y-axis value (dispersion). Default is to
#' take the standard deviation of all values/
#' @param do.plot Plot the average/dispersion relationship
#' @param set.var.genes Set object@@var.genes to the identified variable genes
#' (default is TRUE)
#' @param x.low.cutoff Bottom cutoff on x-axis for identifying variable genes
#' @param x.high.cutoff Top cutoff on x-axis for identifying variable genes
#' @param y.cutoff Bottom cutoff on y-axis for identifying variable genes
#' @param y.high.cutoff Top cutoff on y-axis for identifying variable genes
#' @param num.bin Total number of bins to use in the scaled analysis (default
#' is 20)
#' @param do.recalc TRUE by default. If FALSE, plots and selects variable genes without recalculating statistics for each gene.
#' @param sort.results If TRUE (by default), sort results in object@hvg.info in decreasing order of dispersion
#' @param do.cpp Run c++ version of mean.function and dispersion.function if they
#' exist.
#' @param display.progress show progress bar for calculations
#' @param ... Extra parameters to VariableGenePlot
#' @inheritParams VariableGenePlot
#'
#' @importFrom stats sd
#' @importFrom MASS kde2d
#' @importFrom utils txtProgressBar setTxtProgressBar
#'
#' @return Returns a Seurat object, placing variable genes in object@@var.genes.
#' The result of all analysis is stored in object@@hvg.info
#'
#' @seealso \code{VariableGenePlot}
#'
#' @rdname FindVariableGenes
#' @export FindVariableGenes
#'
#' @examples
#' pbmc_small <- FindVariableGenes(object = pbmc_small, do.plot = FALSE)
#' pbmc_small@var.genes
#'
setGeneric(
  name = 'FindVariableGenes',
  def = function(
    object,
    ...
  ) {
    return(standardGeneric(f = 'FindVariableGenes'))
  }
)
