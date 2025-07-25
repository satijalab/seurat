#' Support for Ultra-Large Matrices with spam
#'
#' @name spam-integration
#' @rdname spam-integration
#'
#' @description
#' Seurat now provides experimental support for ultra-large sparse matrices 
#' using the \code{spam} package. This enables analysis of datasets that exceed 
#' the 2^31 element limit of standard \code{dgCMatrix} objects.
#'
#' @details
#' The \code{spam} package provides sparse matrix classes with 64-bit integer 
#' indexing support, allowing matrices with more than 2.1 billion elements. 
#' This is particularly useful for:
#' \itemize{
#'   \item Large atlas-scale single-cell datasets
#'   \item Multi-sample integration with hundreds of thousands of cells
#'   \item Spatial transcriptomics datasets with millions of spatial locations
#'   \item Multi-modal datasets (e.g., CITE-seq, spatial + scRNA-seq)
#' }
#'
#' @section Usage:
#' To use spam matrices with Seurat:
#' \enumerate{
#'   \item Install the spam package: \code{install.packages("spam")}
#'   \item Create spam matrices using \code{spam::spam()} or convert existing 
#'         matrices with \code{as.spam()}
#'   \item Use spam matrices directly with \code{CreateSeuratObject()} and 
#'         other Seurat functions
#' }
#'
#' @section Supported Functions:
#' The following Seurat functions have been enhanced to support spam matrices:
#' \itemize{
#'   \item \code{CreateSeuratObject()} - Create Seurat objects with spam matrices
#'   \item \code{NormalizeData()} - CLR and RC normalization methods
#'   \item \code{FindVariableFeatures()} - Variable feature selection using VST
#'   \item \code{ScaleData()} - Data scaling (with memory considerations)
#'   \item Matrix operations: \code{rowSums}, \code{rowMeans}, variance calculations
#' }
#'
#' @section Memory Considerations:
#' While spam matrices can handle larger datasets, they still require careful 
#' memory management:
#' \itemize{
#'   \item Very large spam matrices (>2GB) may require conversion to dgCMatrix 
#'         for some operations
#'   \item Consider using subsampling or sketching for exploratory analysis
#'   \item Monitor memory usage during intensive computations
#' }
#'
#' @section Example:
#' \preformatted{
#' # Install spam package
#' if (!requireNamespace("spam", quietly = TRUE)) {
#'   install.packages("spam")
#' }
#' 
#' library(spam)
#' library(Seurat)
#' 
#' # Create a large sparse matrix using spam
#' # This example shows a matrix that would exceed dgCMatrix limits
#' large_matrix <- spam(0, nrow = 50000, ncol = 50000)  # 2.5B elements
#' 
#' # Add some non-zero values (simulating gene expression)
#' for (i in 1:1000) {
#'   large_matrix[sample(50000, 1), sample(50000, 1)] <- rpois(1, lambda = 2)
#' }
#' 
#' # Create Seurat object (this would fail with dgCMatrix)
#' seurat_obj <- CreateSeuratObject(counts = large_matrix, 
#'                                  min.cells = 1, min.features = 1)
#' 
#' # Standard Seurat workflow
#' seurat_obj <- NormalizeData(seurat_obj, normalization.method = "CLR")
#' seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst")
#' }
#'
#' @section Performance Tips:
#' \itemize{
#'   \item Use \code{CheckMatrixSize()} to verify if spam matrices are needed
#'   \item Consider BPCells for very large on-disk matrices as an alternative
#'   \item Profile memory usage with \code{object.size()} before operations
#'   \item For integration workflows, subset to variable features first
#' }
#'
#' @seealso 
#' \code{\link[spam]{spam}}, \code{\link{CreateSeuratObject}}, 
#' \code{\link{CheckMatrixSize}}
#'
#' @concept spam
#' @concept large-matrices
#' @concept ultra-large
#' @concept integration
#'
NULL
