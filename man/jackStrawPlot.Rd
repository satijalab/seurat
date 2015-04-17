\name{jackStrawPlot}
\alias{jackStrawPlot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
JackStraw Plot
}
\description{
Plots the results of the JackStraw analysis for PCA significance. For each PC, plots a QQ-plot comparing the distribution of p-values for all genes across each PC, compared with a uniform distribution. Also determines a p-value for the overall significance of each PC.
}
\usage{
jackStrawPlot(object, plot.lim = 0.8, PCs = 1:5, nCol = 3, score.thresh = 1e-05, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
Seurat plot
}
  \item{plot.lim}{
Y-axis maximum on each QQ plot.
}
  \item{PCs}{
Which PCs to examine
}
  \item{nCol}{
Number of columns
}
  \item{score.thresh}{
Threshold to use for the proportion test of PC significance (see Details)
}
}
\details{
The p-value for each PC is based on a proportion test comparing the number of genes with a p-value below a particular threshold, compared with the proportion of genes expected under a uniform distribution of p-values.}
\value{
No value returned, just the PC plots.
}

\author{
Thanks to Omri Wurtzel for integrating with ggplot
}
