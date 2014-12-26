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
structure(function (object, genes.use = NULL, pcs.use = NULL, 
    do.plot = TRUE, do.reorder = FALSE, reorder.numeric = FALSE) 
standardGeneric("buildClusterTree"), generic = structure("buildClusterTree", package = "Seurat"), package = "Seurat", group = list(), valueClass = character(0), signature = c("object", 
"genes.use", "pcs.use", "do.plot", "do.reorder", "reorder.numeric"
), default = `\001NULL\001`, skeleton = (function (object, genes.use = NULL, 
    pcs.use = NULL, do.plot = TRUE, do.reorder = FALSE, reorder.numeric = FALSE) 
stop("invalid call in method dispatch to 'buildClusterTree' (no default method)", 
    domain = NA))(object, genes.use, pcs.use, do.plot, do.reorder, 
    reorder.numeric), class = structure("standardGeneric", package = "methods"))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
