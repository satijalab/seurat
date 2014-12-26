\name{viz.pca}
\alias{viz.pca}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Visualize PCA genes
}
\description{
Visualize top genes associated with principal components
}
\usage{
viz.pca(object, pcs.use = 1:5, num.genes = 15, use.full = FALSE, font.size = 0.5, nCol = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
Seurat object
}
  \item{pcs.use}{
Number of PCs to display
}
  \item{num.genes}{
Number of genes to display
}
  \item{use.full}{
Use full PCA (i.e. the projected PCA, by default FALSE)
}
  \item{font.size}{
Font size
}
  \item{nCol}{
Number of columns to display
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Displays the PCA gene plot, returns no value
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
viz.pca(nbt)
viz.pca(nbt,pcs.use=1:6,num.genes=10,font.size=0.65,nCol=4)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
