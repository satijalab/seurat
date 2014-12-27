\name{marker.test}
\alias{marker.test}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
ROC-based marker discovery
}
\description{
Identifies 'markers' of gene expression using ROC analysis. For each gene, evaluates (using AUC) a classifier built on that gene alone, to classify between two groups of cells. 
}
\usage{
marker.test(object, cells.1, cells.2, genes.use = NULL, thresh.use = log(2))
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
An AUC value of 1 means that expression values for this gene alone can perfectly classify the two groupings (i.e. Each of the  cells in cells.1 exhibit a higher level than each of the cells in cells.2). An AUC value of 0 also means there is perfect classification, but in the other direction. A value of 0.5 implies that the gene has no predictive power to classify the two groups.
}
\value{
Returns a 'predictive power' (abs(AUC-0.5)) ranked matrix of putative differentially expressed genes. 
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

