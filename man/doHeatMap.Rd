\name{doHeatMap}
\alias{doHeatMap}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Gene expression heatmap
}
\description{
Draws a heatmap of single cell gene expression using the heatmap.2 function.
}
\usage{
doHeatMap(object, cells.use = NULL, genes.use = NULL, disp.min = -2.5, disp.max = 2.5, draw.line = TRUE, do.return = FALSE, order.by.ident = TRUE, col.use = pyCols, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
  Seurat object
}
  \item{cells.use}{
  Cells to include in the heatmap (default is all cells)
}
  \item{genes.use}{
  Genes to include in the heatmap (ordered)
}
  \item{disp.min}{
  Minimum display value (all values below are clipped)
}
  \item{disp.max}{
  Maximum display value (all values above are clipped)
}
  \item{draw.line}{
  Draw vertical lines delineating cells in different identity classes.
}
  \item{do.return}{
Default is FALSE. If TRUE, return a matrix of scaled values which would be passed to heatmap.2
}
  \item{order.by.ident}{
Order cells in the heatmap by identity class (default is TRUE). If FALSE, cells are ordered based on their order in cells.use
}
  \item{col.use}{
Color palette to use
}
  \item{\dots}{
Additional parameters to heatmap.2. Common examples are cexRow and cexCol, which set row and column text sizes}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
If do.return==TRUE, a matrix of scaled values which would be passed to heatmap.2. Otherwise, no return value, only a graphical output
}