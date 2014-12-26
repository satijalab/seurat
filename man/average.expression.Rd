\name{average.expression}
\alias{average.expression}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Averaged gene expression by identity class
}
\description{
Returns gene expression for an 'average' single cell in each identity class
}
\usage{
average.expression(object, genes.use = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
  Seurat object
}
  \item{genes.use}{
Genes to analyze. Default is all genes.
}
}
\details{
Output is in log-space, but averaging is done in non-log space.
}
\value{
Returns a matrix with genes as rows, identity classes as columns.

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
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
