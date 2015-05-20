\name{set.ident}
\alias{set.ident}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Set identity class information
}
\description{
Sets the identity class value for a subset (or all) cells
}
\usage{
set.ident(object, cells.use = NULL, ident.use = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
Seurat object
}
  \item{cells.use}{
Vector of cells to set identity class info for (default is all cells)
}
  \item{ident.use}{
Vector of identity class values to assign (character vector) 
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A Seurat object where object@ident has been appropriately modified
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