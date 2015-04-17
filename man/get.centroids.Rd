\name{get.centroids}
\alias{get.centroids}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Get cell centroids
}
\description{
Calculate the spatial mapping centroids for each cell, based on previously calculated mapping probabilities for each bin.
}
\usage{
get.centroids(object, cells.use = NULL, get.exact = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
Seurat object
}
  \item{cells.use}{
Cells to calculate centroids for (default is all cells)
}
  \item{get.exact}{
Get exact centroid (Default is TRUE). If FALSE, identify the single closest bin.
}
}
\details{
Currently, Seurat assumes that the tissue of interest has an 8x8 bin structure. This will be broadened in a future release.
}
\value{
Data frame containing the x and y coordinates for each cell centroid.
}
