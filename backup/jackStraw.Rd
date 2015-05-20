\name{jackStraw}
\alias{jackStraw}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Determine statistical significance of PCA scores.
}
\description{
Randomly permutes a subset of data, and calculates projected PCA scores for these 'random' genes. Then compares the PCA scores for the 'random' genes with the observed PCA scores to determine statistical signifance. End result is a p-value for each gene's association with each principal component.
}
\usage{
jackStraw(object, num.pc = 10, num.replicate = 100, prop.freq = 0.01, do.print = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
Seurat object
}
  \item{num.pc}{
Number of PCs to compute significance for
}
  \item{num.replicate}{
Number of replicate samplings to perform
}
  \item{prop.freq}{
Proportion of the data to randomly permute for each replicate
}
  \item{do.print}{
Print the number of replicates that have been processed.
}
}
\value{
Returns a Seurat object where object@jackStraw.empP represents p-values for each gene in the PCA analysis. If project.pca is subsequently run, object@jackStraw.empP.full then represents p-values for all genes.
}
\references{
Inspired by Chung et al, Bioinformatics (2014)
}
