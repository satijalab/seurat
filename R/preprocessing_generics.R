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
  def = function(object, ...) {
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
#' @param chunk.size Chunk size to iterate over
#' @param normalized.data Full path to normalized data in loom file
#' @param display.progress Show progress bar for calculations
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
  def = function(object, ...) {
    return(standardGeneric(f = 'FindVariableGenes'))
  }
)

#' Scale and center the data.
#'
#' Scales and centers genes in the dataset. If variables are provided in vars.to.regress,
#' they are individually regressed against each gene, and the resulting residuals are
#' then scaled and centered.
#'
#' ScaleData now incorporates the functionality of the function formerly known
#' as RegressOut (which regressed out given the effects of provided variables
#' and then scaled the residuals). To make use of the regression functionality,
#' simply pass the variables you want to remove to the vars.to.regress parameter.
#'
#' Setting center to TRUE will center the expression for each gene by subtracting
#' the average expression for that gene. Setting scale to TRUE will scale the
#' expression level for each gene by dividing the centered gene expression
#' levels by their standard deviations if center is TRUE and by their root mean
#' square otherwise.
#'
#'
#' @param object Seurat object
#' @param genes.use Vector of gene names to scale/center. Default is all genes
#' in object@@data.
#' @param data.use Can optionally pass a matrix of data to scale, default is
#' object@data[genes.use, ]
#' @param vars.to.regress Variables to regress out (previously latent.vars in
#' RegressOut). For example, nUMI, or percent.mito.
#' @param model.use Use a linear model or generalized linear model
#' (poisson, negative binomial) for the regression. Options are 'linear'
#' (default), 'poisson', and 'negbinom'
#' @param use.umi Regress on UMI count data. Default is FALSE for linear
#' modeling, but automatically set to TRUE if model.use is 'negbinom' or 'poisson'
#' @param do.scale Whether to scale the data.
#' @param do.center Whether to center the data.
#' @param scale.max Max value to return for scaled data. The default is 10.
#' Setting this can help reduce the effects of genes that are only expressed in
#' a very small number of cells. If regressing out latent variables and using a
#' non-linear model, the default is 50.
#' @param block.size Default size for number of genes to scale at in a single
#' computation. Increasing block.size may speed up calculations but at an
#' additional memory cost.
#' @param min.cells.to.block If object contains fewer than this number of cells,
#' don't block for scaling calculations.
#' @param display.progress Displays a progress bar for scaling procedure
#' @param assay.type Assay to scale data for. Default is RNA. Can be changed for
#' multimodal analyses.
#' @param do.cpp By default (TRUE), most of the heavy lifting is done in c++.
#' We've maintained support for our previous implementation in R for
#' reproducibility (set this to FALSE) as results can change slightly due to
#' differences in numerical precision which could affect downstream calculations.
#' @param check.for.norm Check to see if data has been normalized, if not,
#' output a warning (TRUE by default)
#' @param chunk.size Chunk size to iterate over
#' @param normalized.data Full path to normalized data in loom file
#' @param Overwrite Overwrite existing dataset with name 'layers/\code{name}'
#'
#' @return Returns a seurat object with object@@scale.data updated with scaled
#' and/or centered data.
#'
#' @rdname ScaleData
#' @export ScaleData
#'
#' @examples
#' pbmc_small <- ScaleData(object = pbmc_small)
#' \dontrun{
#' # To regress out certain effects
#' pbmc_small = ScaleData(object = pbmc_small, vars.to.regress = effects_list)
#' }
#'
setGeneric(
  name = 'ScaleData',
  def = function(object, ...) {
    return(standardGeneric(f = 'ScaleData'))
  }
)
