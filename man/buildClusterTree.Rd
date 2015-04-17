\name{buildClusterTree}
\alias{buildClusterTree}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Phylogenetic Analysis of Identity Classes
}
\description{
Constructs a phylogenetic tree relating the 'average' cell from each identity class. Tree is estimated based on a distance matrix constructed in either gene expression space or PCA space.

Note that the tree is calculated for an 'average' cell, so gene expression or PC scores are averaged across all cells in an identity class before the tree is constructed.

}
\usage{
buildClusterTree(object, genes.use = NULL, pcs.use = NULL, do.plot = TRUE, do.reorder = FALSE, reorder.numeric = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
Seurat object
}
  \item{genes.use}{
Genes to use for the analysis. Default is the set of variable genes (object@var.genes). Assumes pcs.use=NULL (tree calculated in gene expression space)
}
  \item{pcs.use}{
If set, tree is calculated in PCA space, using the eigenvalue-weighted euclidean distance across these PC scores.
}
  \item{do.plot}{
Plot the resulting phylogenetic tree
}
  \item{do.reorder}{
Re-order identity classes (factor ordering), according to position on the tree. This groups similar classes together, which can be helpful for example in violin plots.
}
  \item{reorder.numeric}{
Re-order identity classes according to position on the tree, assigning a numeric value ('1' is the leftmost node)
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A Seurat object where the cluster tree is stored in object@cluster.tree[[1]]
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