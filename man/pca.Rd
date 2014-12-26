\name{pca}
\alias{pca}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Run Principal Component Analysis on gene expression
}
\description{
Run prcomp for PCA dimensionality reduction
}
\usage{
pca(object, pc.genes = NULL, do.print = TRUE, pcs.print = 5, pcs.store = 40, genes.print = 30, use.imputed = FALSE,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
Seurat object
}
  \item{pc.genes}{
Genes to use as input for PCA. Default is object@var.genes  
}
 \item{do.print}{
Print the top genes associated with high/low loadings for the PCs
}
  \item{ics.print}{
Number of PCs to print genes for
}
  \item{ics.store}{
Number of PCs to store
}
  \item{genes.print}{
Number of genes to print for each PC
}
  \item{use.imputed}{
Run PCA on imputed values (FALSE by default)
}
  \item{\dots} {
Additional arguments to be passed to prcomp
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Returns Seurat object with an PCA embedding (object@pca.rot) and gene projection matrix (object@pca.x). The PCA object itself is stored in object@pca.obj[[1]]
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
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
