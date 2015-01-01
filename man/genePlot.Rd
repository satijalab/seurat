\name{genePlot}
\alias{genePlot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Scatter plot of single cell data
}
\description{
Creates a scatter plot of two features (typically gene expression), across a set of single cells. Cells are colored by their identity class.
}
\usage{
genePlot(object, gene1, gene2, cell.ids = NULL, col.use = NULL, pch.use = 16, cex.use = 2, no.col = FALSE, inc.first = FALSE, use.imputed = FALSE, by.k = FALSE, do.ident = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
Seurat object
}
  \item{gene1}{
First feature to plot. Typically gene expression but can also be metrics, PC scores, etc. - anything that can be retreived with fetch.data
}
  \item{gene2}{
Second feature to plot.
}
  \item{cell.ids}{
Cells to include on the scatter plot.
}
  \item{col.use}{
Colors to use for identity class plotting.
}
  \item{pch.use}{
Pch argument for plotting
}
  \item{cex.use}{
Cex argument for plotting
}
  \item{use.imputed}{
Use imputed values for gene expression (Default is FALSE)
}
  \item{do.ident}{
False by default. If TRUE, 
}
  \item{\dots}{
Additional arguments to be passed to plot.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
No return, only graphical output
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
structure(function (object, gene1, gene2, cell.ids = NULL, col.use = NULL, 
    pch.use = 16, cex.use = 2, no.col = FALSE, inc.first = FALSE, 
    use.imputed = FALSE, by.k = FALSE, do.ident = FALSE, ...) 
standardGeneric("genePlot"), generic = structure("genePlot", package = "Seurat"), package = "Seurat", group = list(), valueClass = character(0), signature = c("object", 
"gene1", "gene2", "cell.ids", "col.use", "pch.use", "cex.use", 
"no.col", "inc.first", "use.imputed", "by.k", "do.ident"), default = `\001NULL\001`, skeleton = (function (object, 
    gene1, gene2, cell.ids = NULL, col.use = NULL, pch.use = 16, 
    cex.use = 2, no.col = FALSE, inc.first = FALSE, use.imputed = FALSE, 
    by.k = FALSE, do.ident = FALSE, ...) 
stop("invalid call in method dispatch to 'genePlot' (no default method)", 
    domain = NA))(object, gene1, gene2, cell.ids, col.use, pch.use, 
    cex.use, no.col, inc.first, use.imputed, by.k, do.ident, 
    ...), class = structure("standardGeneric", package = "methods"))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
