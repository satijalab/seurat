\name{pca.sig.genes}
\alias{pca.sig.genes}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Significant genes from a PCA
}
\description{
Returns a set of genes, based on the jackStraw analysis, that have statistically significant associations with a set of PCs.
}
\usage{
pca.sig.genes(object, pcs.use, pval.cut = 0.1, use.full = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
Seurat object
}
  \item{pcs.use}{
PCS to use.
}
  \item{pval.cut}{
P-value cutoff
}
  \item{use.full}{
Use the full list of genes (from the projected PCA). Assumes that project.pca has been run. Default is TRUE
}
}
\value{
A vector of genes whose p-values are statistically significant for at least one of the given PCs.
}
