\name{doHeatMap}
\alias{doHeatMap}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
To add
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
doHeatMap(object, cells.use = NULL, genes.use = NULL, disp.min = -2.5, disp.max = 2.5, draw.line = TRUE, do.return = FALSE, order.by.ident = TRUE, col.use = pyCols, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
%%     ~~Describe \code{object} here~~
}
  \item{cells.use}{
%%     ~~Describe \code{cells.use} here~~
}
  \item{genes.use}{
%%     ~~Describe \code{genes.use} here~~
}
  \item{disp.min}{
%%     ~~Describe \code{disp.min} here~~
}
  \item{disp.max}{
%%     ~~Describe \code{disp.max} here~~
}
  \item{draw.line}{
%%     ~~Describe \code{draw.line} here~~
}
  \item{do.return}{
%%     ~~Describe \code{do.return} here~~
}
  \item{order.by.ident}{
%%     ~~Describe \code{order.by.ident} here~~
}
  \item{col.use}{
%%     ~~Describe \code{col.use} here~~
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
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
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
structure(function (object, cells.use = NULL, genes.use = NULL, 
    disp.min = -2.5, disp.max = 2.5, draw.line = TRUE, do.return = FALSE, 
    order.by.ident = TRUE, col.use = pyCols, ...) 
standardGeneric("doHeatMap"), generic = structure("doHeatMap", package = "Seurat"), package = "Seurat", group = list(), valueClass = character(0), signature = c("object", 
"cells.use", "genes.use", "disp.min", "disp.max", "draw.line", 
"do.return", "order.by.ident", "col.use"), default = `\001NULL\001`, skeleton = (function (object, 
    cells.use = NULL, genes.use = NULL, disp.min = -2.5, disp.max = 2.5, 
    draw.line = TRUE, do.return = FALSE, order.by.ident = TRUE, 
    col.use = pyCols, ...) 
stop("invalid call in method dispatch to 'doHeatMap' (no default method)", 
    domain = NA))(object, cells.use, genes.use, disp.min, disp.max, 
    draw.line, do.return, order.by.ident, col.use, ...), class = structure("standardGeneric", package = "methods"))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
