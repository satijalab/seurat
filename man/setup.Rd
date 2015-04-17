\name{setup}
\alias{setup}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Setup Seurat object
}
\description{
Setup and initialize basic parameters of the Seurat object
}
\usage{
setup(object, project, min.cells = 3, min.genes = 2500, is.expr = 1, do.scale = TRUE, do.center = TRUE, names.field = 1, names.delim = "_", meta.data = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
  Seurat object
}
  \item{project}{
  Project name (string)
}
  \item{min.cells}{
  Include genes with detected expression in at least this many cells
}
  \item{min.genes}{
  Include cells where at least this many genes are detected
}
  \item{is.expr}{
  Expression threshold for 'detected' gene 
}
  \item{do.scale}{
  In object@scale.data, perform row-scaling (gene-based z-score)
}
  \item{do.center}{
  In object@scale.data, perform row-centering (gene-based centering)
}

  \item{names.field}{
  For the initial identity class for each cell, choose this field from the cell's column name
}
  \item{names.delim}{
  For the initial identity class for each cell, choose this delimiter from the cell's column name
}
  \item{meta.data}{
  Additional metadata to add to the Seurat object. Should be a data frame where the rows are cell names, and the columns are additional metadata fields
}
  \item{\dots}{
  Additional arguemnts, not currently used.
}
}

\value{
Seurat object. Fields modified include object@data, object@scale.data, object@data.info, object@data.ngene, object@ident 
}
