\name{ica}
\alias{ica}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Run Independent Component Analysis on gene expression
}
\description{
Run fastICA algorithm for ICA dimensionality reduction
}
\usage{
ica(object, ic.genes = NULL, do.print = TRUE, ics.print = 5, ics.store = 10, genes.print = 30, use.imputed = FALSE,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
Seurat object
}
  \item{ic.genes}{
Genes to use as input for ICA. Default is object@var.genes  
}
  \item{do.print}{
Print the top genes associated with high/low loadings for the ICs
}
  \item{ics.print}{
Number of ICs to print genes for
}
  \item{ics.store}{
Number of ICs to store
}
  \item{genes.print}{
Number of genes to print for each IC
}
  \item{use.imputed}{
Run ICA on imputed values (FALSE by default)
}
  \item{\dots}{
  Additional arguments to be passed to fastICA  
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Returns Seurat object with an ICA embedding (object@ica.rot) and gene projection matrix (object@ica.x). The ICA object itself is stored in object@ica.obj[[1]]
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
