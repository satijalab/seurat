\name{reorder.by.tree}
\alias{reorder.by.tree}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Re-order identity classes based on phylogeny
}
\description{
Re-order (factor ordering) identity classes based on their position in the phylogenetic tree. This places similar classes together, which can be helpful for visualization (violin plots etc.)
}
\usage{
reorder.by.tree(object, genes.use = NULL, do.numeric = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
%%     ~~Describe \code{object} here~~
}
  \item{genes.use}{
Not needed, should be removed
}
  \item{do.numeric}{
Also relabels the identity classes as numbers (1 as the leftmost class)
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Returns Seurat object with a modified object@cluster.tree 
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
structure(function (object, genes.use = NULL, do.numeric = FALSE, 
    ...) 
standardGeneric("reorder.by.tree"), generic = structure("reorder.by.tree", package = "Seurat"), package = "Seurat", group = list(), valueClass = character(0), signature = c("object", 
"genes.use", "do.numeric"), default = `\001NULL\001`, skeleton = (function (object, 
    genes.use = NULL, do.numeric = FALSE, ...) 
stop("invalid call in method dispatch to 'reorder.by.tree' (no default method)", 
    domain = NA))(object, genes.use, do.numeric, ...), class = structure("standardGeneric", package = "methods"))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
