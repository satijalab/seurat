# Example: Processing Multiple 10X Samples with Seurat
# 
# This example demonstrates how to use the new ProcessMultiple10X function
# to efficiently process multiple 10X Genomics datasets in parallel

library(Seurat)

# Example 1: Auto-detect and process all samples in a directory
# This is the most common use case - CellRanger outputs in subdirectories
data_directory <- "/path/to/cellranger/outputs/"

# Auto-detect all valid 10X samples
samples <- Detect10XSamples(data_directory)
print(samples)

# Process all detected samples in parallel
seurat_list <- ProcessMultiple10X(
  data_dir = data_directory,
  sample_dirs = samples,  # or set to NULL for auto-detection
  matrix_format = "sparse",  # or "delayed" for memory efficiency
  mc.cores = 4,
  save_file = "all_samples.qs",
  verbose = TRUE
)

# Example 2: Process specific samples with custom settings
specific_samples <- c("sample1", "sample2", "sample3")

seurat_list <- ProcessMultiple10X(
  data_dir = "/path/to/data/",
  sample_dirs = specific_samples,
  matrix_subdir = "filtered_feature_bc_matrix",  # default CellRanger output
  matrix_format = "auto",  # automatically choose best format
  mc.cores = 8,
  save_file = "processed_samples.qs",
  verbose = TRUE
)

# Example 3: For very large datasets, use DelayedArray for memory efficiency
large_seurat_list <- ProcessMultiple10X(
  data_dir = "/path/to/large/datasets/",
  sample_dirs = NULL,  # auto-detect
  matrix_format = "delayed",  # memory-efficient format
  mc.cores = 16,
  verbose = TRUE
)

# Example 4: Load previously processed data
loaded_samples <- Load10XData("all_samples.qs", nthreads = 10)

# Access individual samples
sample1 <- loaded_samples[["sample1"]]
print(sample1)

# Example 5: Traditional approach (what the function replaces)
# This is what you would do manually:
# data_dir <- "/path/to/data/"
# list_data <- c("sample1", "sample2", "sample3")
# paths <- paste0(data_dir, list_data, "/filtered_feature_bc_matrix/")
# 
# NeuroProx <- mclapply(X = paths, 
#                      FUN = function(paths){
#                        CreateSeuratObject(counts = Read10X(data.dir = paths), 
#                                          project = str_split(paths, pattern = "/")[[1]][7])
#                      },
#                      mc.cores = detectCores()
#                      )
# names(NeuroProx) <- str_split(paths, pattern = "/", simplify = T)[, 7]
# qsave(x = NeuroProx, file = "20250528_NeuroProx_full.qs", nthreads = 10)

# The new approach is much simpler and more robust:
NeuroProx <- ProcessMultiple10X(
  data_dir = data_dir,
  sample_dirs = list_data,
  save_file = "20250528_NeuroProx_full.qs",
  mc.cores = parallel::detectCores(),
  verbose = TRUE
)

# Load later:
NeuroProx <- Load10XData("20250528_NeuroProx_full.qs", nthreads = 10)

cat("Processing complete! Objects ready for downstream analysis.\n")
