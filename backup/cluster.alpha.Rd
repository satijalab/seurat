\name{cluster.alpha}
\alias{cluster.alpha}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Probability of detection by identity class
}
\description{
For each gene, calculates the probability of detection for each identity class.
}
\usage{
cluster.alpha(object, thresh.min = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
Seurat object
}
  \item{thresh.min}{
Minimum threshold to define 'detected' (log-scale)
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Returns a matrix with genes as rows, identity classes as columns.
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
