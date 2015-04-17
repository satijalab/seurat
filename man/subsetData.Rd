\name{subsetData}
\alias{subsetData}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Return a subset of the Seurat object
}
\description{
Creates a Seurat object containing only a subset of the cells in the original object. Takes either a list of cells to use as a subset, or a parameter (for example, a gene), to subset on.
}
\usage{
subsetData(object, subset.name = NULL, cells.use = NULL, accept.low = -Inf, accept.high = Inf, do.center = TRUE, do.scale = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
Seurat object
}
  \item{cells.use}{
  A vector of cell names to use as a subset. If NULL (default), then this list will be computed based on the next three arguments.
}
  \item{subset.name}{
Parameter to subset on. Eg, the name of a gene, PC1, a column name in object@data.info, etc. Any argument that can be retreived using fetch.data
}
  \item{accept.low}{
Low cutoff for the parameter (default is -Inf)
}
  \item{accept.high}{
High cutoff for the parameter (default is Inf)
}
  \item{do.center}{
Recenter the new object@scale.data
}
  \item{do.scale}{
Rescale the new object@scale.data
}
  \item{\dots}{
Additional arguments to be passed to fetch.data (for example, use.imputed=TRUE)
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Returns a Seurat object containing only the relevant subset of cells
}