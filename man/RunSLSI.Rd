% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generics.R, R/dimensional_reduction.R
\name{RunSLSI}
\alias{RunSLSI}
\alias{RunSLSI.default}
\alias{RunSLSI.Assay}
\alias{RunSLSI.StdAssay}
\alias{RunSLSI.Seurat}
\title{Run Supervised Latent Semantic Indexing}
\usage{
RunSLSI(object, ...)

\method{RunSLSI}{default}(
  object,
  assay = NULL,
  n = 50,
  reduction.key = "SLSI_",
  graph = NULL,
  verbose = TRUE,
  seed.use = 42,
  ...
)

\method{RunSLSI}{Assay}(
  object,
  assay = NULL,
  features = NULL,
  n = 50,
  reduction.key = "SLSI_",
  graph = NULL,
  verbose = TRUE,
  seed.use = 42,
  ...
)

\method{RunSLSI}{StdAssay}(
  object,
  assay = NULL,
  features = NULL,
  n = 50,
  reduction.key = "SLSI_",
  graph = NULL,
  layer = "data",
  verbose = TRUE,
  seed.use = 42,
  ...
)

\method{RunSLSI}{Seurat}(
  object,
  assay = NULL,
  features = NULL,
  n = 50,
  reduction.name = "slsi",
  reduction.key = "SLSI_",
  graph = NULL,
  verbose = TRUE,
  seed.use = 42,
  ...
)
}
\arguments{
\item{object}{An object}

\item{...}{Arguments passed to IRLBA irlba}

\item{assay}{Name of Assay SLSI is being run on}

\item{n}{Total Number of SLSI components to compute and store}

\item{reduction.key}{dimensional reduction key, specifies the string before
the number for the dimension names}

\item{graph}{Graph used supervised by SLSI}

\item{verbose}{Display messages}

\item{seed.use}{Set a random seed. Setting NULL will not set a seed.}

\item{features}{Features to compute SLSI on. If features=NULL, SLSI will be run
using the variable features for the Assay5.}

\item{layer}{Layer to run SLSI on}

\item{reduction.name}{dimensional reduction name}
}
\value{
Returns Seurat object with the SLSI calculation stored in the
reductions slot
}
\description{
Run a supervised LSI (SLSI) dimensionality reduction supervised by a
cell-cell kernel. SLSI is used to capture a linear transformation of peaks
that maximizes its dependency to the given cell-cell kernel.
}
\concept{dimensional_reduction}
