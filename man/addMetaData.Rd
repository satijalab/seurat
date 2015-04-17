\name{addMetaData}
\alias{addMetaData}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Add Metadata
}
\description{
Adds additional data for single cells to the Seurat object. Can be any piece of information associated with a cell (examples include read depth, alignment rate, experimental batch, or subpopulation identity). The advantage of adding it to the Seurat object is so that it can be analyzed/visualized using fetch.data, vlnPlot, genePlot, subsetData, etc.
}
\usage{
addMetaData(object, metadata)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
Seurat object
}
  \item{metadata}{
Data frame where the row names are cell names (note : these must correspond exactly to the items in object@cell.names), and the columns are additional metadata items.
}
}

\value{
Seurat object where the additional metadata has been added as columns in object@data.info
}
