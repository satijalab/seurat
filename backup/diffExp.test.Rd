\name{diffExp.test}
\alias{diffExp.test}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Likelihood ratio test for zero-inflated data
}
\description{
Identifies differentially expressed genes between two groups of cells using the LRT model proposed in Mcdavid et al, Bioinformatics, 2011
}
\usage{
diffExp.test(object, cells.1, cells.2, genes.use = NULL, thresh.use = log(2))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
Seurat object
}
  \item{cells.1}{
Group 1 cells
}
  \item{cells.2}{
Group 2 cells
}
  \item{genes.use}{
Genes to test. Default is to use all genes.
}
  \item{thresh.use}{
Limit testing to genes which show, on average, at least X-fold difference (log-scale) between the two groups of cells.

Increasing thresh.use speeds up the function, but can miss weaker signals. 
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Returns a p-value ranked matrix of putative differentially expressed genes. 
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
