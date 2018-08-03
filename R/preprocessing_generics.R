#' Normalize Data
#'
#' Normalize the count data present in a given assay.
#'
#' @param object An object
#' @param normalization.method Method for normalization.
#'  \itemize{
#'   \item{LogNormalize: }{Feature counts for each cell are divided by the total
#'   counts for that cell and multiplied by the scale.factor. This is then
#'   natural-log transformed using log1p.}
#'   \item{CLR: }{Applies a centered log ratio transformation}
#' }
#' More methods to be added.
#' @param scale.factor Sets the scale factor for cell-level normalization
#' @param verbose display progress bar for normalization procedure.
#'
#' @return Returns object after normalization
#'
#' @rdname NormalizeData
#' @export NormalizeData
#'
NormalizeData <- function(
  object,
  normalization.method,
  scale.factor,
  verbose
) {
  UseMethod(generic = 'NormalizeData', object = object)
}

#' Find variable features
#'
#' Identifies features that are outliers on a 'mean variability plot'.
#'
#' For the mean.var.plot method:
#' Exact parameter settings may vary empirically from dataset to dataset, and
#' based on visual inspection of the plot. Setting the y.cutoff parameter to 2
#' identifies features that are more than two standard deviations away from the
#' average dispersion within a bin. The default X-axis function is the mean
#' expression level, and for Y-axis it is the log(Variance/mean). All mean/variance
#' calculations are not performed in log-space, but the results are reported in
#' log-space - see relevant functions for exact details.
#'
#' @param object An object
#' @param selection.method How to choose top variable features. Choose one of :
#' \itemize{
#'   \item{vst:}{ First, fits a line to the relationship of log(variance) and
#'   log(mean) using local polynomial regression (loess). Then standardizes the
#'   feature values using the observed mean and expected variance (given by the
#'   fitted line). Feature variance is the calculated on the standardized values
#'   after clipping to a maximum (default is 50).}
#'   \item{mean.var.plot:}{ First, uses a function to calculate average
#'   expression (mean.function) and dispersion (dispersion.function) for each
#'   feature. Next, divides features into num.bin (deafult 20) bins based on
#'   their average expression, and calculates z-scores for dispersion within
#'   each bin. The purpose of this is to identify variable features while
#'   controlling for the strong relationship between variability and average
#'   expression.}
#'   \item{dispersion:}{ selects the genes with the highest dispersion values}
#' }
#' @param loess.span (vst method) Loess span parameter used when fitting the
#' variance-mean relationship
#' @param clip.max (vst method) After standardization values larger than
#' clip.max will be set to clip.max; default is 'auto' which sets this value to
#' the square root of the number of cells
#' @param mean.function Function to compute x-axis value (average expression).
#'  Default is to take the mean of the detected (i.e. non-zero) values
#' @param dispersion.function Function to compute y-axis value (dispersion).
#' Default is to take the standard deviation of all values
#' @param num.bin Total number of bins to use in the scaled analysis (default
#' is 20)
#' @param binning.method Specifies how the bins should be computed. Available
#' methods are:
#' \itemize{
#'   \item{equal_width:}{ each bin is of equal width along the x-axis [default]}
#'   \item{equal_frequency:}{ each bin contains an equal number of features (can
#'   increase statistical power to detect overdispersed features at high
#'   expression values, at the cost of reduced resolution along the x-axis)}
#' }
#' @param verbose show progress bar for calculations
#'
#' @rdname FindVariableFeatures
#' @export FindVariableFeatures
#'
#' @aliases FindVariableGenes
#'
FindVariableFeatures <- function(
  object,
  selection.method,
  mean.function,
  dispersion.function,
  num.bin,
  binning.method,
  verbose,
  ...
) {
  UseMethod(generic = 'FindVariableFeatures', object = object)
}

#' Scale and center the data.
#'
#' Scales and centers features in the dataset. If variables are provided in vars.to.regress,
#' they are individually regressed against each feautre, and the resulting residuals are
#' then scaled and centered.
#'
#' ScaleData now incorporates the functionality of the function formerly known
#' as RegressOut (which regressed out given the effects of provided variables
#' and then scaled the residuals). To make use of the regression functionality,
#' simply pass the variables you want to remove to the vars.to.regress parameter.
#'
#' Setting center to TRUE will center the expression for each feautre by subtracting
#' the average expression for that feautre. Setting scale to TRUE will scale the
#' expression level for each feautre by dividing the centered feautre expression
#' levels by their standard deviations if center is TRUE and by their root mean
#' square otherwise.
#'
#' @param object An object
#' @param features.use Vector of features names to scale/center. Default is all features
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
#' Setting this can help reduce the effects of feautres that are only expressed in
#' a very small number of cells. If regressing out latent variables and using a
#' non-linear model, the default is 50.
#' @param block.size Default size for number of feautres to scale at in a single
#' computation. Increasing block.size may speed up calculations but at an
#' additional memory cost.
#' @param min.cells.to.block If object contains fewer than this number of cells,
#' don't block for scaling calculations.
#' @param verbose Displays a progress bar for scaling procedure
#'
#' @rdname ScaleData
#' @export ScaleData
#'
ScaleData <- function(
  object,
  features.use,
  vars.to.regress,
  model.use,
  use.umi,
  do.scale,
  do.center,
  scale.max,
  block.size,
  min.cells.to.block,
  verbose,
  ...
) {
  UseMethod(generic = 'ScaleData', object = object)
}
