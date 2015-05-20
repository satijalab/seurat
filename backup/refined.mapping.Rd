\name{refined.mapping}
\alias{refined.mapping}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Quantitative refinement of spatial inferences
}
\description{
Refines the initial mapping with more complex models that allow gene expression to vary quantitatively across bins (instead of 'on' or 'off'), and that also considers the covariance structure between genes.
}
\usage{
refined.mapping(object, genes.use)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
Seurat object
}
  \item{genes.use}{
Genes to use to drive the refinement procedure. 
}
}
\details{
Full details given in spatial mapping manuscript.
}
\value{
Seurat object, where mapping probabilities for each bin are stored in object@final.prob
}