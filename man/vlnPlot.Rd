\name{vlnPlot}
\alias{vlnPlot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Single cell violin plot
}
\description{
Draws a violin plot of single cell data (gene expression, metrics, PC scores, etc.)
}
\usage{
vlnPlot(object, nCol = NULL, ylab.max = 12, do.ret = TRUE, do.sort = FALSE, size.x.use = 16, size.y.use = 16, size.title.use = 20, use.imputed = FALSE, adjust.use = 1, size.use = 1, cols.use = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
Seurat object
}
  \item{features.plot}{
Features to plot (gene expression, metrics, PC scores, anything that can be retreived by fetch.data)
}
  \item{nCol}{
Number of columns if multiple plots are displayed
}
  \item{ylab.max}{
Maximum ylab value
}
  \item{do.ret}{
FALSE by default. If TRUE, returns a list of ggplot objects.
}
  \item{do.sort}{
Sort identity classes (on the x-axis) by the average expression of the attribute being potted 
}
  \item{size.x.use}{
X axis font size
}
  \item{size.y.use}{
Y axis font size
}
  \item{size.title.use}{
Title font size
}
  \item{use.imputed}{
Use imputed values for gene expression (default is FALSE)
}
  \item{adjust.use}{
Adjust parameter for geom_violin
}
  \item{size.use}{
Point size for geom_violin
}
  \item{cols.use}{
Colors to use for plotting
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
By default, no return, only graphical output. If do.return==TRUE, returns a list of ggplot objects.
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
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
structure(function (object, genes.plot, nCol = NULL, ylab.max = 12, 
    do.ret = TRUE, do.sort = FALSE, size.x.use = 16, size.y.use = 16, 
    size.title.use = 20, by.k = FALSE, use.imputed = FALSE, adjust.use = 1, 
    size.use = 1, cols.use = NULL, ...) 
standardGeneric("vlnPlot"), generic = structure("vlnPlot", package = "Seurat"), package = "Seurat", group = list(), valueClass = character(0), signature = c("object", 
"genes.plot", "nCol", "ylab.max", "do.ret", "do.sort", "size.x.use", 
"size.y.use", "size.title.use", "by.k", "use.imputed", "adjust.use", 
"size.use", "cols.use"), default = `\001NULL\001`, skeleton = (function (object, 
    genes.plot, nCol = NULL, ylab.max = 12, do.ret = TRUE, do.sort = FALSE, 
    size.x.use = 16, size.y.use = 16, size.title.use = 20, by.k = FALSE, 
    use.imputed = FALSE, adjust.use = 1, size.use = 1, cols.use = NULL, 
    ...) 
stop("invalid call in method dispatch to 'vlnPlot' (no default method)", 
    domain = NA))(object, genes.plot, nCol, ylab.max, do.ret, 
    do.sort, size.x.use, size.y.use, size.title.use, by.k, use.imputed, 
    adjust.use, size.use, cols.use, ...), class = structure("standardGeneric", package = "methods"))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
