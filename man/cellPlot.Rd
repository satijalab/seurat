\name{cellPlot}
\alias{cellPlot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Cell-cell scatter plot 
}
\description{
Creates a plot of scatter plot of genes across two single cells 
}
\usage{
cellPlot(object, cell1, cell2, gene.ids = NULL, col.use = "black", nrpoints.use = Inf, pch.use = 16, cex.use = 0.5, do.ident = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
Seurat object
}
  \item{cell1}{
Cell 1 name (can also be a number, representing the position in object@cell.names)
}
  \item{cell2}{
Cell 2 name (can also be a number, representing the position in object@cell.names)
}
  \item{gene.ids}{
Genes to plot (default, all genes)
}
  \item{col.use}{
Colors to use for the points
}
  \item{nrpoints.use}{
Parameter for smoothScatter
}
  \item{pch.use}{
Point symbol to use
}
  \item{cex.use}{
Point size
}
  \item{do.ident}{
FALSE by default. If TRUE, allows you to click on individual points to reveal gene names (hit ESC to stop)
}
  \item{\dots}{
Additional arguments to pass to smoothScatter
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
No return value (plots a scatter plot)
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