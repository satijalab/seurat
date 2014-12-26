\name{project.pca}
\alias{project.pca}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Project Principal Components Analysis onto full dataset
}
\description{
Takes a pre-computed PCA (typically calculated on a subset of genes) and projects this onto the entire dataset (all genes). Note that the cell loadings (PCA rotation matrices) remains unchanged, but now there are gene scores for all genes.
}
\usage{
project.pca(object, do.print = TRUE, pcs.print = 5, pcs.store = 20, genes.print = 30, replace.pc = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
  Seurat object
}
  \item{do.print}{
  Print top genes associated with the projected PCs
}
  \item{pcs.print}{
  Number of PCs to print genes for
}
  \item{pcs.store}{
  Number of PCs to store (default is 20)
}
  \item{genes.print}{
  Number of genes with highest/lowest loadings to print for each PC
}
  \item{replace.pc}{
  Replace the existing PCA (overwite object@pca.x), not done by default.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Returns Seurat object with the projected PCA values in object@pca.x.full
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
structure(function (object, pc.genes = NULL, do.print = TRUE, 
    pcs.print = 5, pcs.store = 20, genes.print = 30, replace.pc = FALSE) 
standardGeneric("project.pca"), generic = structure("project.pca", package = "Seurat"), package = "Seurat", group = list(), valueClass = character(0), signature = c("object", 
"pc.genes", "do.print", "pcs.print", "pcs.store", "genes.print", 
"replace.pc"), default = `\001NULL\001`, skeleton = (function (object, 
    pc.genes = NULL, do.print = TRUE, pcs.print = 5, pcs.store = 20, 
    genes.print = 30, replace.pc = FALSE) 
stop("invalid call in method dispatch to 'project.pca' (no default method)", 
    domain = NA))(object, pc.genes, do.print, pcs.print, pcs.store, 
    genes.print, replace.pc), class = structure("standardGeneric", package = "methods"))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
