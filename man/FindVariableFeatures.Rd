% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generics.R, R/preprocessing.R
\name{FindVariableFeatures}
\alias{FindVariableFeatures}
\alias{FindVariableGenes}
\alias{FindVariableFeatures.V3Matrix}
\alias{FindVariableFeatures.Assay}
\alias{FindVariableFeatures.SCTAssay}
\alias{FindVariableFeatures.Seurat}
\title{Find variable features}
\usage{
FindVariableFeatures(object, ...)

\method{FindVariableFeatures}{V3Matrix}(
  object,
  selection.method = "vst",
  loess.span = 0.3,
  clip.max = "auto",
  mean.function = FastExpMean,
  dispersion.function = FastLogVMR,
  num.bin = 20,
  binning.method = "equal_width",
  verbose = TRUE,
  ...
)

\method{FindVariableFeatures}{Assay}(
  object,
  selection.method = "vst",
  loess.span = 0.3,
  clip.max = "auto",
  mean.function = FastExpMean,
  dispersion.function = FastLogVMR,
  num.bin = 20,
  binning.method = "equal_width",
  nfeatures = 2000,
  mean.cutoff = c(0.1, 8),
  dispersion.cutoff = c(1, Inf),
  verbose = TRUE,
  ...
)

\method{FindVariableFeatures}{SCTAssay}(object, nfeatures = 2000, ...)

\method{FindVariableFeatures}{Seurat}(
  object,
  assay = NULL,
  selection.method = "vst",
  loess.span = 0.3,
  clip.max = "auto",
  mean.function = FastExpMean,
  dispersion.function = FastLogVMR,
  num.bin = 20,
  binning.method = "equal_width",
  nfeatures = 2000,
  mean.cutoff = c(0.1, 8),
  dispersion.cutoff = c(1, Inf),
  verbose = TRUE,
  ...
)
}
\arguments{
\item{object}{An object}

\item{...}{Arguments passed to other methods}

\item{selection.method}{How to choose top variable features. Choose one of :
\itemize{
  \item \dQuote{\code{vst}}:  First, fits a line to the relationship of
    log(variance) and log(mean) using local polynomial regression (loess).
    Then standardizes the feature values using the observed mean and
    expected variance (given by the fitted line). Feature variance is then
    calculated on the standardized values
    after clipping to a maximum (see clip.max parameter).
  \item \dQuote{\code{mean.var.plot}} (mvp): First, uses a function to
    calculate average expression (mean.function) and dispersion
    (dispersion.function) for each feature. Next, divides features into
    \code{num.bin} (default 20) bins based on their average expression,
    and calculates z-scores for dispersion within each bin. The purpose of
    this is to identify variable features while controlling for the
    strong relationship between variability and average expression
  \item \dQuote{\code{dispersion}} (disp): selects the genes with the
    highest dispersion values
}}

\item{loess.span}{(vst method) Loess span parameter used when fitting the
variance-mean relationship}

\item{clip.max}{(vst method) After standardization values larger than
clip.max will be set to clip.max; default is 'auto' which sets this value to
the square root of the number of cells}

\item{mean.function}{Function to compute x-axis value (average expression).
Default is to take the mean of the detected (i.e. non-zero) values}

\item{dispersion.function}{Function to compute y-axis value (dispersion).
Default is to take the standard deviation of all values}

\item{num.bin}{Total number of bins to use in the scaled analysis (default
is 20)}

\item{binning.method}{Specifies how the bins should be computed. Available
methods are:
\itemize{
  \item \dQuote{\code{equal_width}}: each bin is of equal width along the
    x-axis (default)
  \item \dQuote{\code{equal_frequency}}: each bin contains an equal number
    of features (can increase statistical power to detect overdispersed
    features at high expression values, at the cost of reduced resolution
    along the x-axis)
}}

\item{verbose}{show progress bar for calculations}

\item{nfeatures}{Number of features to select as top variable features;
only used when \code{selection.method} is set to \code{'dispersion'} or
\code{'vst'}}

\item{mean.cutoff}{A two-length numeric vector with low- and high-cutoffs for
feature means}

\item{dispersion.cutoff}{A two-length numeric vector with low- and high-cutoffs for
feature dispersions}

\item{assay}{Assay to use}
}
\description{
Identifies features that are outliers on a 'mean variability plot'.
}
\details{
For the mean.var.plot method:
Exact parameter settings may vary empirically from dataset to dataset, and
based on visual inspection of the plot. Setting the y.cutoff parameter to 2
identifies features that are more than two standard deviations away from the
average dispersion within a bin. The default X-axis function is the mean
expression level, and for Y-axis it is the log(Variance/mean). All mean/variance
calculations are not performed in log-space, but the results are reported in
log-space - see relevant functions for exact details.
}
\concept{preprocessing}
