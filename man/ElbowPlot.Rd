% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/visualization.R
\name{ElbowPlot}
\alias{ElbowPlot}
\title{Quickly Pick Relevant Dimensions}
\usage{
ElbowPlot(object, ndims = 20, reduction = "pca")
}
\arguments{
\item{object}{Seurat object}

\item{ndims}{Number of dimensions to plot standard deviation for}

\item{reduction}{Reduction technique to plot standard deviation for}
}
\value{
A ggplot object
}
\description{
Plots the standard deviations (or approximate singular values if running PCAFast)
of the principle components for easy identification of an elbow in the graph.
This elbow often corresponds well with the significant dims and is much faster to run than
Jackstraw
}
\examples{
data("pbmc_small")
ElbowPlot(object = pbmc_small)

}
\concept{visualization}
