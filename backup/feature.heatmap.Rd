\name{feature.heatmap}
\alias{feature.heatmap}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Vizualization of multiple features 
}
\description{
Similar to feature.plot, however, also splits the plot by visualizing each identity class separately.

Particularly useful for seeing if the same groups of cells co-exhibit a common feature (i.e. co-express a gene), even within an identity class. Best understood by example.
}
\usage{
feature.heatmap(object, features.plot, pc.1 = 1, pc.2 = 2, idents.use = NULL, pt.size = 2, cols.use = heat.colors(10), pch.use = 16, reduction.use = "tsne", margin.dist = 5, font.size = 0.7)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
Seurat object
}
  \item{features.plot}{
Vector of features to plot
}
  \item{pc.1}{
Dimension for x-axis (default 1)
}
  \item{pc.2}{
Dimension for y-axis (default 2)
}
  \item{idents.use}{
Which identity classes to display (default is all identity classes)
}
  \item{pt.size}{
Adjust point size for plotting
}
  \item{cols.use}{
Ordered vector of colors to use for plotting. Default is heat.colors(10). 
}
  \item{pch.use}{
Pch for plotting
}
  \item{reduction.use}{
Which dimensionality reduction to use. Default is "tsne", can also be "pca", or "ica", assuming these are precomputed.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
No return value, only a graphical output


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

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
