#' Deprecated function(s) in the Seurat package
#'
#' These functions are provided for compatibility with older version of the Seurat package.  They may eventually be completely removed.
#' @rdname Seurat-deprecated
#' @name Seurat-deprecated
#' @param ... Parameters to be passed to the modern version of the function
#' @export setup vlnPlot subsetData mean.var.plot pca project.pca print.pca viz.pca set.ident pca.plot pcHeatmap jackStraw jackStrawPlot run_tsne tsne.plot find.markers find_all_markers genePlot feature.plot tsne.plot
#' @aliases setup vlnPlot subsetData mean.var.plot pca project.pca print.pca viz.pca set.ident pca.plot pcHeatmap jackStraw jackStrawPlot run_tsne tsne.plot find.markers find_all_markers genePlot feature.plot tnse.plot
#' @section Details:
#' \tabular{rl}{
#'   \code{setup} \tab now a synonym for \code{\link{Setup}}\cr
#'   \code{vlnPlot} \tab now a synonym for \code{\link{VlnPlot}}\cr
#'   \code{subsetData} \tab now a synonym for \code{\link{SubsetData}}\cr
#'   \code{mean.var.plot} \tab now a synonym for \code{\link{MeanVarPlot}}\cr
#'   \code{pca} \tab now a synonym for \code{\link{PCA}}\cr
#'   \code{project.pca} \tab now a synonym for \code{\link{ProjectPCA}}\cr
#'   \code{print.pca} \tab now a synonym for \code{\link{PrintPCA}}\cr
#'   \code{viz.pca} \tab now a synonym for \code{\link{VizPCA}}\cr
#'   \code{set.ident} \tab now a synonym for \code{\link{SetIdent}}\cr
#'   \code{pca.plot} \tab now a synonym for \code{\link{PCAPlot}}\cr
#'   \code{pcHeatmap} \tab now a synonym for \code{\link{PCHeatmap}}\cr
#'   \code{jackStraw} \tab now a synonym for \code{\link{JackStraw}}\cr
#'   \code{jackStrawPlot} \tab now a synonym for \code{\link{JackStrawPlot}}\cr
#'   \code{run_tsne} \tab now a synonym for \code{\link{RunTSNE}}\cr
#'   \code{tsne.plot} \tab now a synonym for \code{\link{TSNEPlot}}\cr
#'   \code{find.markers} \tab now a synonym for \code{\link{FindMarkers}}\cr
#'   \code{find_all_markers} \tab now a synonym for \code{\link{FindAllMarkers}}\cr
#'   \code{genePlot} \tab now a synonym for \code{\link{GenePlot}}\cr
#'   \code{feature.plot} \tab now a synonym for \code{\link{FeaturePlot}}\cr
#' }
#'

setup <- function(...) {
    .Deprecated("Setup", package="Seurat")
    Setup(...)
}

vlnPlot <- function(...) {
    .Deprecated("VlnPlot", package="Seurat")
    VlnPlot(...)
}

subsetData <- function(...) {
    .Deprecated("SubsetData", package="Seurat")
    SubsetData(...)
}

mean.var.plot <- function(...) {
    .Deprecated("MeanVarPlot", package="Seurat")
    MeanVarPlot(...)
}

pca <- function(...) {
    .Deprecated("PCA", package="Seurat")
    PCA(...)
}

project.pca <- function(...) {
    .Deprecated("ProjectPCA", package="Seurat")
    ProjectPCA(...)
}

print.pca <- function(...) {
    .Deprecated("PrintPCA", package="Seurat")
    PrintPCA(...)
}

viz.pca <- function(...) {
    .Deprecated("VizPCA", package="Seurat")
    VizPCA(...)
}

set.ident <- function(...) {
    .Deprecated("SetIdent", package="Seurat")
    SetIdent(...)
}

pca.plot <- function(...) {
    .Deprecated("PCAPlot", package="Seurat")
    PCAPlot(...)
}

pcHeatmap <- function(...) {
    .Deprecated("PCHeatmap", package="Seurat")
    PCHeatmap(...)
}

jackStraw <- function(...) {
    .Deprecated("jackStraw", package="Seurat")
    JackStraw(...)
}

jackStrawPlot <- function(...) {
    .Deprecated("JackStrawPlot", package="Seurat")
    JackStrawPlot(...)
}

run_tsne <- function(...) {
    .Deprecated("RunTSNE", package="Seurat")
    RunTSNE(...)
}

tsne.plot <- function(...) {
    .Deprecated("TSNEPlot", package="Seurat")
    TSNEPlot(...)
}

find.markers <- function(...) {
    .Deprecated("FindMarkers", package="Seurat")
    FindMarkers(...)
}

find_all_markers <- function(...) {
    .Deprecated("FindAllMarkers", package="Seurat")
    FindAllMarkers(...)
}

genePlot <- function(...) {
    .Deprecated("GenePlot", package="Seurat")
    GenePlot(...)
}

feature.plot <- function(...) {
    .Deprecated("FeaturePlot", package="Seurat")
    FeaturePlot(...)
}
