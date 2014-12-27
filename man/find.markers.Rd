\name{find.markers}
\alias{find.markers}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Gene expression markers of identity classes
}
\description{
Finds markers (differentially expressed genes) for identity classes
}
\usage{
find.markers(object, stat.1, stat.2 = NULL, genes.use = NULL, thresh.use = log(2), test.use = "bimod")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
Seurat object
}
  \item{stat.1}{
Identity class to define markers for
}
  \item{stat.2}{
A second identity class for comparison. If NULL (default) - use all other cells for comparison.
}
  \item{genes.use}{
Genes to test. Default is to use all genes.
}
  \item{thresh.use}{
Limit testing to genes which show, on average, at least X-fold difference (log-scale) between the two groups of cells.

Increasing thresh.use speeds up the function, but can miss weaker signals. 
}
  \item{test.use}{
Denotes which test to use. Seurat currently implements "bimod" (likelihood-ratio test for single cell gene expression, McDavid et al., Bioinformatics, 2011, default), "roc" (standard AUC classifier), "t" (Students t-test), and "tobit" (Tobit-test for differential gene expression, as in Trapnell et al., Nature Biotech, 2014)
}
}
\details{
}
\value{
Matrix containing a ranked list of putative markers, and associated statistics (p-values, ROC score, etc.)
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
