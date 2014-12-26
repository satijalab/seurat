\name{run_tsne}
\alias{run_tsne}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Run t-distributed Stochastic Neighbor Embedding
}
\description{
Run t-SNE dimensionality reduction on selected features (PCA space). 
}
\usage{
run_tsne(object, cells.use = NULL, pcs.use = 1:10, do.plot = TRUE, k.seed = 1, do.fast = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
  Seurat object
}
  \item{cells.use}{
  Which cells to analyze (default, all cells)
}
  \item{pcs.use}{
  Which PCs to use as features
}
  \item{do.plot}{
  Plot resulting tSNE map
}
  \item{k.seed}{
  Random seed for the t-SNE
}
  \item{do.fast}{
  If TRUE, uses the Barnes-hut implementation, which runs faster, but is less flexible
}
  \item{\dots}{
  Additional arguments to the tSNE call. Most commonly used is perplexity (expected number of neighbors)
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Returns a Seurat object with a tSNE embedding in object@tsne_rot
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
