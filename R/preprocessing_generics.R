#' Normalize Assay Data
#'
#' @param object An object
#' @param normalization.method Method for normalization. Default is
#' log-normalization (LogNormalize). More methods to be added very shortly.
#' @param scale.factor Sets the scale factor for cell-level normalization
#' @param display.progress display progress bar for scaling procedure.
#' @param ... Arguments passed to other methods
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
  display.progress,
  ...
) {
  UseMethod(generic = 'NormalizeData', object = object)
}

#' Find variable features
#'
#' Identifies features that are outliers on a 'mean variability plot'. First, uses
#' a function to calculate average expression (mean.function) and dispersion (dispersion.function)
#' for each feautre Next, divides features into num.bin (deafult 20) bins based on
#' their average expression, and calculates z-scores for dispersion within each
#' bin. The purpose of this is to identify variable features while controlling for
#' the strong relationship between variability and average expression.
#'
#' Exact parameter settings may vary empirically from dataset to dataset, and
#' based on visual inspection of the plot.
#' Setting the y.cutoff parameter to 2 identifies features that are more than two standard
#' deviations away from the average dispersion within a bin. The default X-axis function
#' is the mean expression level, and for Y-axis it is the log(Variance/mean). All mean/variance
#' calculations are not performed in log-space, but the results are reported in log-space -
#' see relevant functions for exact details.
#'
#' @param object An object
#' @param mean.function Function to compute x-axis value (average expression). Default
#' is to take the mean of the detected (i.e. non-zero) values
#' @param dispersion.function Function to compute y-axis value (dispersion). Default is to
#' take the standard deviation of all values
#' @param num.bin Total number of bins to use in the scaled analysis (default
#' is 20)
#' @param binning.method Specifies how the bins should be computed. Available methods are:
#' \itemize{
#'   \item{equal_width:}{ each bin is of equal width along the x-axis [default]}
#'   \item{equal_frequency:}{ each bin contains an equal number of features (can increase
#' statistical power to detect overdispersed features at high expression values, at
#' the cost of reduced resolution along the x-axis)}
#' }
#' @param display.progress show progress bar for calculations
#' @param ... Arguments passed to other methods
#'
FindVariableFeatures <- function(
  object,
  mean.function,
  dispersion.function,
  num.bin,
  binning.method,
  display.progress,
  ...
) {
  UseMethod(generic = 'FindVariableFeatures', object = object)
}

#' Get the top variable Features
#'
#' @param object An object
#' @param num.features Number of top variable features to get
#' @param mean.cutoff Cutoff for feature means, must be a two-length numeric for low and high values
#' @param dispersion.cutoff Cutoff for feature dispersion values, must be a two-length numeric for low and high values
#' @param selection.method Specifies how to select the features to store in the object
#' \itemize{
#'   \item{mean.var.plot: }{Default method, placing cutoffs on the mean variablility plot}
#'   \item{dispersion: }{Choose the top.features with the highest dispersion}
#' }
#' @param ... Arguments passed to other methods
#'
#' @rdname GetVariableFeatures
#' @export GetVariableFeatures
#'
GetVariableFeatures <- function(
  object,
  num.features,
  mean.cutoff,
  dispersion.cutoff,
  selection.method,
  ...
) {
  UseMethod(generic = 'GetVariableFeatures', object = object)
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
#' @param display.progress Displays a progress bar for scaling procedure
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
  display.progress,
  ...
) {
  UseMethod(generic = 'ScaleData', object = object)
}
