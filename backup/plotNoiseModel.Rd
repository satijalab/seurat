\name{plotNoiseModel}
\alias{plotNoiseModel}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Visualize expression/dropout curve
}
\description{
Plot the probability of detection vs average expression of a gene.

Assumes that this 'noise' model has been precomputed
}
\usage{
plotNoiseModel(object, cell.ids = c(1, 2), col.use = "black", lwd.use = 2, do.new = TRUE, x.lim = 10, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
Seurat object
}
  \item{cell.ids}{
Cells to use
}
  \item{col.use}{
Color code or name
}
  \item{lwd.use}{
Line width for curve
}
  \item{do.new}{
Create a new plot (default) or add to existing
}
  \item{x.lim}{
Maximum value for X axis
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{

Returns no value, displays a plot
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
