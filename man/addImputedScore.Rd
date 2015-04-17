\name{addImputedScore}
\alias{addImputedScore}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate imputed expression values
}
\description{
Uses L1-constrained linear models (LASSO) to impute single cell gene expression values.
}
\usage{
addImputedScore(object, genes.use = NULL, genes.fit = NULL, s.use = 20, do.print = FALSE, gram = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
Seurat object
}
  \item{genes.use}{
A vector of genes (predictors) that can be used for building the LASSO models.
}
  \item{genes.fit}{
A vector of genes to impute values for
}
  \item{s.use}{
Maximum number of steps taken by the algorithm (lower values indicate a greater degree of smoothing)
}
  \item{do.print}{
Print progress (output the name of each gene after it has been imputed).
}
  \item{gram}{
The use.gram argument passed to lars
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Returns a Seurat object where the imputed values have been added to object@imputed
}
