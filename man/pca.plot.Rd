\name{pca.plot}
\alias{pca.plot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Dimensional reduction plot
}
\description{
Graphs the output of a dimensional reduction technique (PCA by default). Cells are colored by their identity class.
}
\usage{
pca.plot(object, pc.1 = 1, pc.2 = 2, cells.use = NULL, pt.size = 4, do.return = FALSE, do.bare = FALSE, by.k = FALSE, cols.use = NULL, reduction.use = "pca")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
Seurat object
}
  \item{pc.1}{
Dimension for x-axis (default 1)
}
  \item{pc.2}{
Dimension for y-axis (default 2)
}
  \item{cells.use}{
Vector of cells to plot (default is all cells)
}
  \item{pt.size}{
Adjust point size for plotting
}
  \item{do.return}{
Return a ggplot2 object (default : FALSE)
}
  \item{do.bare}{
Do only minimal formatting (default : FALSE)
}
  \item{cols.use}{
Vector of colors, each color corresponds to an identity class. By default, ggplot assigns
}
  \item{reduction.use}{
Which dimensionality reduction to use. Default is "pca", can also be "tsne", or "ica", assuming these are precomputed.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
If do.return==TRUE, returns a ggplot2 object. Otherwise, only graphical output.
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
