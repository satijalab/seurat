\name{plotClusterTree}
\alias{plotClusterTree}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot phylogenetic tree
}
\description{
Plots previously computed phylogenetic tree (from buildClusterTree)
}
\usage{
plotClusterTree(object)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
Seurat object
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Plots dendogram, returns no value
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
structure(function (object) 
standardGeneric("plotClusterTree"), generic = structure("plotClusterTree", package = "Seurat"), package = "Seurat", group = list(), valueClass = character(0), signature = "object", default = `\001NULL\001`, skeleton = (function (object) 
stop("invalid call in method dispatch to 'plotClusterTree' (no default method)", 
    domain = NA))(object), class = structure("standardGeneric", package = "methods"))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
