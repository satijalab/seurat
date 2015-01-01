\name{feature.plot}
\alias{feature.plot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Visualize 'features' on a dimensional reduction plot
}
\description{
Colors single cells on a dimensional reduction plot according to a 'feature' (i.e. gene expression, PC scores, number of genes detected, etc.) 
}
\usage{
feature.plot(object, features.plot, pc.1 = 1, pc.2 = 2, cells.use = NULL, pt.size = 1, cols.use = heat.colors(10), pch.use = 16, reduction.use = "tsne", nCol = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
Seurat object
}
  \item{features.plot}{
Vector of features to plot
}
  \item{pc.1}{
Dimension for x-axis (default 1)
}
  \item{pc.2}{
Dimension for y-axis (default 2)
}
  \item{cells.use}{
Vector of cells to plot (default is all cells)
}
  \item{pt.size}{
Adjust point size for plotting
}
  \item{cols.use}{
Ordered vector of colors to use for plotting. Default is heat.colors(10). 
}
  \item{pch.use}{
Pch for plotting
}
  \item{reduction.use}{
Which dimensionality reduction to use. Default is "tsne", can also be "pca", or "ica", assuming these are precomputed.
}
  \item{nCol}{
Number of columns to use when plotting multiple features.
}
}
\details{
To determine the color, the feature values across all cells are placed into discrete bins, and then assigned a color based on cols.use. The number of bins is determined by the number of colors in cols.use
}
\value{
No return value, only a graphical output
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
structure(function (object, genes.plot, pc.1 = 1, pc.2 = 2, cells.use = NULL, 
    pt.size = 1, do.return = FALSE, do.bare = FALSE, cols.use = heat.colors(10), 
    pch.use = 16, reduction.use = "tsne", nCol = NULL) 
standardGeneric("feature.plot"), generic = structure("feature.plot", package = "Seurat"), package = "Seurat", group = list(), valueClass = character(0), signature = c("object", 
"genes.plot", "pc.1", "pc.2", "cells.use", "pt.size", "do.return", 
"do.bare", "cols.use", "pch.use", "reduction.use", "nCol"), default = `\001NULL\001`, skeleton = (function (object, 
    genes.plot, pc.1 = 1, pc.2 = 2, cells.use = NULL, pt.size = 1, 
    do.return = FALSE, do.bare = FALSE, cols.use = heat.colors(10), 
    pch.use = 16, reduction.use = "tsne", nCol = NULL) 
stop("invalid call in method dispatch to 'feature.plot' (no default method)", 
    domain = NA))(object, genes.plot, pc.1, pc.2, cells.use, 
    pt.size, do.return, do.bare, cols.use, pch.use, reduction.use, 
    nCol), class = structure("standardGeneric", package = "methods"))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
