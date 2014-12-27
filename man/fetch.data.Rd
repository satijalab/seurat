\name{fetch.data}
\alias{fetch.data}
%% ~
\title{
Access cellular data
}
\description{
Retreives data (gene expression, PCA scores, etc, metrics, etc.) for a set of cells in a Seurat object
}
\usage{
fetch.data(object, vars.all = NULL, cells.use = NULL, use.imputed = FALSE, use.scaled = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
Seurat object
}
  \item{vars.all}{
List of all variables to fetch
}
  \item{cells.use}{
Cells to collect data for (default is all cells)
}
  \item{use.imputed}{
For gene expression, use imputed values
}
  \item{use.scaled}{
For gene expression, use scaled values
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
Hi
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
