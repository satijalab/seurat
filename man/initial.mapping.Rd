\name{initial.mapping}
\alias{initial.mapping}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Infer spatial origins for single cells
}
\description{
Probabilistically maps single cells based on (imputed) gene expression estimates, a set of mixture models, and an in situ spatial reference map.  
}
\usage{
initial.mapping(object, cells.use = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
Seurat object
}
  \item{cells.use}{
Which cells to map
}
}

\value{
Seurat object, where mapping probabilities for each bin are stored in object@final.prob
}
