#' DelayedArray Integration for Memory-Efficient Operations
#'
#' Support for DelayedArray and HDF5Array matrices in Seurat for memory-efficient
#' analysis of large single-cell datasets. DelayedArray allows working with
#' matrices larger than RAM by storing data on disk and loading only required
#' portions into memory.
#'
#' @section Supported Operations:
#' \itemize{
#'   \item \code{\link{VST}} - Variance stabilizing transformation with block processing
#'   \item \code{\link{FindVariableFeatures}} - Feature selection with memory efficiency
#'   \item \code{\link{NormalizeData}} - Data normalization (via conversion)
#'   \item \code{\link{ScaleData}} - Data scaling (via conversion)
#'   \item Row-wise statistics: means, sums, variances
#' }
#'
#' @section Memory Considerations:
#' DelayedArray matrices are designed for datasets that exceed available RAM.
#' Key considerations:
#' \itemize{
#'   \item Data remains on disk until explicitly accessed
#'   \item Operations are performed in blocks to minimize memory usage
#'   \item Some operations may require conversion to in-memory formats
#'   \item Performance depends on disk I/O speed and block size
#' }
#'
#' @section Creating DelayedArray Matrices:
#' \preformatted{
#' # From HDF5 file
#' library(HDF5Array)
#' hdf5_matrix <- HDF5Array("path/to/file.h5", "dataset_name")
#' 
#' # From existing matrix
#' library(DelayedArray)
#' delayed_matrix <- DelayedArray(your_matrix)
#' 
#' # Convert to DelayedMatrix for Seurat
#' seurat_ready <- as.DelayedMatrix(your_matrix)
#' }
#'
#' @section Performance Tips:
#' \itemize{
#'   \item Install DelayedMatrixStats for optimized operations
#'   \item Use appropriate block sizes for your system
#'   \item Store frequently accessed data in faster storage (SSD)
#'   \item Consider data layout (row-major vs column-major) for operations
#' }
#'
#' @section Example Usage:
#' \preformatted{
#' # Load large dataset as DelayedArray
#' library(HDF5Array)
#' library(Seurat)
#' 
#' # Assume you have an HDF5 file with expression data
#' expr_data <- HDF5Array("large_dataset.h5", "expression")
#' 
#' # Create Seurat object (will work with DelayedArray)
#' seurat_obj <- CreateSeuratObject(counts = expr_data)
#' 
#' # Find variable features (memory-efficient)
#' seurat_obj <- FindVariableFeatures(seurat_obj, method = "vst")
#' 
#' # Normalize data (will convert to sparse for efficiency)
#' seurat_obj <- NormalizeData(seurat_obj)
#' }
#'
#' @name delayed-integration
#' @rdname delayed-integration
#' @concept delayed
#' @concept memory
#' @concept large-data
#'
NULL

#' Check DelayedArray Support
#'
#' Check if DelayedArray packages are available and provide installation guidance
#'
#' @param verbose Logical; whether to print detailed information
#' @return Logical indicating if DelayedArray support is available
#' @export
#' @concept utilities
#'
CheckDelayedArraySupport <- function(verbose = TRUE) {
  has_delayed <- requireNamespace("DelayedArray", quietly = !verbose)
  has_hdf5 <- requireNamespace("HDF5Array", quietly = !verbose)
  has_stats <- requireNamespace("DelayedMatrixStats", quietly = !verbose)
  
  if (verbose) {
    cat("DelayedArray Support Status:\n")
    cat("  DelayedArray:", if (has_delayed) "✓ Available" else "✗ Not installed", "\n")
    cat("  HDF5Array:", if (has_hdf5) "✓ Available" else "✗ Not installed", "\n")
    cat("  DelayedMatrixStats:", if (has_stats) "✓ Available (recommended)" else "✗ Not installed (recommended)", "\n")
    
    if (!has_delayed || !has_hdf5) {
      cat("\nTo install missing packages:\n")
      if (!has_delayed) cat("  BiocManager::install('DelayedArray')\n")
      if (!has_hdf5) cat("  BiocManager::install('HDF5Array')\n")
      if (!has_stats) cat("  BiocManager::install('DelayedMatrixStats')\n")
    }
  }
  
  return(has_delayed && has_hdf5)
}

#' Convert Matrix to Memory-Efficient Format
#'
#' Automatically choose the most memory-efficient matrix format based on
#' data characteristics and available packages
#'
#' @param x Input matrix
#' @param prefer_format Preferred format: "auto", "sparse", "delayed", or "spam"
#' @param ... Additional arguments passed to conversion functions
#' @return Matrix in optimized format
#' @export
#' @concept utilities
#'
OptimizeMatrixFormat <- function(x, prefer_format = "auto", ...) {
  # Get matrix dimensions and sparsity
  dims <- dim(x)
  total_elements <- prod(dims)
  
  # Estimate sparsity if possible
  sparsity <- 0
  if (inherits(x, "dgCMatrix")) {
    sparsity <- 1 - (length(x@x) / total_elements)
  } else if (is.matrix(x)) {
    sparsity <- sum(x == 0) / total_elements
  }
  
  # Auto-select format if requested
  if (prefer_format == "auto") {
    if (total_elements > 2^31 - 1) {
      # Ultra-large: use spam if available
      if (requireNamespace("spam", quietly = TRUE)) {
        prefer_format <- "spam"
      } else if (requireNamespace("DelayedArray", quietly = TRUE)) {
        prefer_format <- "delayed"
      } else {
        stop("Matrix too large for standard formats. Install 'spam' or 'DelayedArray' packages.")
      }
    } else if (total_elements > 1e8) {
      # Very large: prefer DelayedArray for memory efficiency
      if (requireNamespace("DelayedArray", quietly = TRUE)) {
        prefer_format <- "delayed"
      } else {
        prefer_format <- "sparse"
      }
    } else if (sparsity > 0.9) {
      # Highly sparse: use sparse format
      prefer_format <- "sparse"
    } else {
      # Default: keep as is or convert to sparse if beneficial
      prefer_format <- "sparse"
    }
  }
  
  # Apply conversion
  switch(prefer_format,
    "sparse" = {
      if (inherits(x, "dgCMatrix")) return(x)
      return(as(x, "dgCMatrix"))
    },
    "delayed" = {
      if (!requireNamespace("DelayedArray", quietly = TRUE)) {
        stop("DelayedArray package required for delayed format")
      }
      return(as.DelayedMatrix(x, ...))
    },
    "spam" = {
      if (!requireNamespace("spam", quietly = TRUE)) {
        stop("spam package required for spam format")
      }
      return(spam::as.spam(x))
    },
    {
      # Default: return as is
      return(x)
    }
  )
}
