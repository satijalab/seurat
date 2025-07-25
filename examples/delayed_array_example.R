# Example: Using Seurat with DelayedArray/HDF5Array for memory-efficient analysis
# 
# This example demonstrates how to use DelayedArray matrices with Seurat
# for memory-efficient single-cell analysis

library(Seurat)

# Check if DelayedArray support is available
cat("Checking DelayedArray support...\n")
delayed_support <- CheckDelayedArraySupport(verbose = TRUE)

if (delayed_support && requireNamespace("DelayedArray", quietly = TRUE)) {
  cat("DelayedArray support is available!\n\n")
  
  # Create a simulated large dataset
  cat("Creating simulated dataset...\n")
  n_genes <- 2000
  n_cells <- 500
  
  # Simulate count data (sparse)  
  set.seed(42)
  counts <- matrix(rpois(n_genes * n_cells, lambda = 0.5), 
                   nrow = n_genes, ncol = n_cells)
  rownames(counts) <- paste0("Gene_", 1:n_genes)
  colnames(counts) <- paste0("Cell_", 1:n_cells)
  
  # Convert to DelayedMatrix for memory efficiency
  cat("Converting to DelayedMatrix...\n")
  delayed_counts <- as.DelayedMatrix(counts)
  
  cat("Matrix class:", class(delayed_counts), "\n")
  cat("Matrix dimensions:", dim(delayed_counts), "\n")
  cat("Memory footprint comparison:\n")
  cat("  Original matrix:", format(object.size(counts), units = "MB"), "\n")
  cat("  DelayedMatrix:", format(object.size(delayed_counts), units = "MB"), "\n\n")
  
  # Create Seurat object with DelayedMatrix
  cat("Creating Seurat object with DelayedMatrix...\n")
  seurat_obj <- CreateSeuratObject(counts = delayed_counts, 
                                   project = "DelayedExample")
  
  cat("Seurat object created successfully!\n")
  cat("Number of features:", nrow(seurat_obj), "\n")
  cat("Number of cells:", ncol(seurat_obj), "\n")
  cat("Assay data class:", class(GetAssayData(seurat_obj, slot = "counts")), "\n\n")
  
  # Test VST with DelayedMatrix
  cat("Running VST feature selection...\n")
  hvf_result <- FindVariableFeatures(seurat_obj, 
                                     selection.method = "vst",
                                     nfeatures = 100,
                                     verbose = FALSE)
  
  top_hvf <- head(VariableFeatures(hvf_result), 10)
  cat("Top 10 highly variable features:\n")
  print(top_hvf)
  cat("\n")
  
  # Test matrix optimization
  cat("Testing matrix format optimization...\n")
  
  # Optimize for sparse format
  sparse_matrix <- OptimizeMatrixFormat(counts, prefer_format = "sparse")
  cat("Sparse optimization result:", class(sparse_matrix), "\n")
  
  # Optimize for DelayedMatrix
  delayed_matrix <- OptimizeMatrixFormat(counts, prefer_format = "delayed")  
  cat("DelayedMatrix optimization result:", class(delayed_matrix), "\n")
  
  # Auto-optimization based on matrix properties
  auto_matrix <- OptimizeMatrixFormat(counts, prefer_format = "auto")
  cat("Auto optimization result:", class(auto_matrix), "\n\n")
  
  # Memory usage summary
  cat("=== Memory Usage Summary ===\n")
  cat("DelayedArray support allows handling of matrices that don't fit in RAM\n")
  cat("by processing data in chunks and using disk-backed storage.\n")
  cat("This is particularly useful for:\n")
  cat("- Very large single-cell datasets (>100k cells)\n")
  cat("- Systems with limited RAM\n")
  cat("- HDF5-based data storage\n")
  
} else {
  cat("DelayedArray support not available.\n")
  cat("To enable DelayedArray support, install required packages:\n")
  cat("install.packages(c('DelayedArray', 'DelayedMatrixStats'))\n")
}

cat("\nExample completed!\n")
