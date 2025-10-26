#!/usr/bin/env Rscript

# Test for equality of Seurat objects from do.call vs direct call

cat("Testing equality of SCTransform results\n")
cat("========================================\n\n")

suppressPackageStartupMessages({
  devtools::load_all(".", quiet = TRUE)
})

data("pbmc_small")

cat("Test 1: Direct SCTransform call\n")
set.seed(42)
so_direct <- SCTransform(pbmc_small, verbose = FALSE, vst.flavor = "v2", seed.use = 42)
cat("  Completed\n\n")

cat("Test 2: do.call SCTransform\n")
data("pbmc_small")  # Reload to ensure same starting state
set.seed(42)
so_docall <- do.call(SCTransform, list(object = pbmc_small, verbose = FALSE, vst.flavor = "v2", seed.use = 42))
cat("  Completed\n\n")

cat("========================================\n")
cat("Comparing Objects\n")
cat("========================================\n\n")

# 1. Check dimensions
cat("1. Dimensions:\n")
cat(sprintf("   Direct: %d x %d\n", nrow(so_direct), ncol(so_direct)))
cat(sprintf("   do.call: %d x %d\n", nrow(so_docall), ncol(so_docall)))
cat(sprintf("   Equal: %s\n\n", identical(dim(so_direct), dim(so_docall))))

# 2. Check assays
cat("2. Assays:\n")
assays_direct <- Assays(so_direct)
assays_docall <- Assays(so_docall)
cat(sprintf("   Direct: %s\n", paste(assays_direct, collapse = ", ")))
cat(sprintf("   do.call: %s\n", paste(assays_docall, collapse = ", ")))
cat(sprintf("   Equal: %s\n\n", identical(assays_direct, assays_docall)))

# 3. Check SCT assay specifically
if ("SCT" %in% assays_direct && "SCT" %in% assays_docall) {
  cat("3. SCT Assay Data:\n")
  
  # Get layers
  layers_direct <- Layers(so_direct[["SCT"]])
  layers_docall <- Layers(so_docall[["SCT"]])
  cat(sprintf("   Direct layers: %s\n", paste(layers_direct, collapse = ", ")))
  cat(sprintf("   do.call layers: %s\n", paste(layers_docall, collapse = ", ")))
  
  # Compare counts
  cat("\n   a) Counts layer:\n")
  counts_direct <- LayerData(so_direct, assay = "SCT", layer = "counts")
  counts_docall <- LayerData(so_docall, assay = "SCT", layer = "counts")
  cat(sprintf("      Dimensions: direct %s, do.call %s\n", 
              paste(dim(counts_direct), collapse = "x"),
              paste(dim(counts_docall), collapse = "x")))
  counts_equal <- all.equal(counts_direct, counts_docall)
  cat(sprintf("      all.equal: %s\n", 
              if(isTRUE(counts_equal)) "TRUE" else as.character(counts_equal)))
  cat(sprintf("      identical: %s\n", identical(counts_direct, counts_docall)))
  
  # Compare data
  cat("\n   b) Data layer:\n")
  data_direct <- LayerData(so_direct, assay = "SCT", layer = "data")
  data_docall <- LayerData(so_docall, assay = "SCT", layer = "data")
  cat(sprintf("      Dimensions: direct %s, do.call %s\n", 
              paste(dim(data_direct), collapse = "x"),
              paste(dim(data_docall), collapse = "x")))
  data_equal <- all.equal(data_direct, data_docall, tolerance = 1e-10)
  cat(sprintf("      all.equal: %s\n", 
              if(isTRUE(data_equal)) "TRUE" else as.character(data_equal)))
  cat(sprintf("      identical: %s\n", identical(data_direct, data_docall)))
  
  # Compare scale.data
  cat("\n   c) Scale.data layer:\n")
  scale_direct <- LayerData(so_direct, assay = "SCT", layer = "scale.data")
  scale_docall <- LayerData(so_docall, assay = "SCT", layer = "scale.data")
  cat(sprintf("      Dimensions: direct %s, do.call %s\n", 
              paste(dim(scale_direct), collapse = "x"),
              paste(dim(scale_docall), collapse = "x")))
  scale_equal <- all.equal(scale_direct, scale_docall, tolerance = 1e-10)
  cat(sprintf("      all.equal: %s\n", 
              if(isTRUE(scale_equal)) "TRUE" else as.character(scale_equal)))
  cat(sprintf("      identical: %s\n", identical(scale_direct, scale_docall)))
}

# 4. Check variable features
cat("\n4. Variable Features:\n")
var_direct <- VariableFeatures(so_direct)
var_docall <- VariableFeatures(so_docall)
cat(sprintf("   Direct: %d features\n", length(var_direct)))
cat(sprintf("   do.call: %d features\n", length(var_docall)))
cat(sprintf("   identical: %s\n", identical(var_direct, var_docall)))
if (!identical(var_direct, var_docall)) {
  overlap <- length(intersect(var_direct, var_docall))
  cat(sprintf("   Overlap: %d features (%.1f%%)\n", 
              overlap, 100 * overlap / length(var_direct)))
  cat("\n   First 10 from direct:\n")
  print(head(var_direct, 10))
  cat("\n   First 10 from do.call:\n")
  print(head(var_docall, 10))
}

# 5. Sample some actual values
cat("\n5. Sample Values Check:\n")
if ("SCT" %in% assays_direct) {
  # Check a specific gene in counts
  test_gene <- rownames(so_direct)[1]
  val_direct <- LayerData(so_direct, assay = "SCT", layer = "counts")[test_gene, 1:5]
  val_docall <- LayerData(so_docall, assay = "SCT", layer = "counts")[test_gene, 1:5]
  cat(sprintf("   Gene: %s (counts, first 5 cells)\n", test_gene))
  cat("   Direct: ", paste(as.numeric(val_direct), collapse = ", "), "\n")
  cat("   do.call:", paste(as.numeric(val_docall), collapse = ", "), "\n")
  cat(sprintf("   identical: %s\n", identical(as.numeric(val_direct), as.numeric(val_docall))))
  
  # Check a specific gene in scale.data
  if (test_gene %in% rownames(LayerData(so_direct, assay = "SCT", layer = "scale.data"))) {
    val_direct_scale <- LayerData(so_direct, assay = "SCT", layer = "scale.data")[test_gene, 1:5]
    val_docall_scale <- LayerData(so_docall, assay = "SCT", layer = "scale.data")[test_gene, 1:5]
    cat(sprintf("\n   Gene: %s (scale.data, first 5 cells)\n", test_gene))
    cat("   Direct: ", paste(sprintf("%.6f", val_direct_scale), collapse = ", "), "\n")
    cat("   do.call:", paste(sprintf("%.6f", val_docall_scale), collapse = ", "), "\n")
    cat(sprintf("   identical: %s\n", identical(val_direct_scale, val_docall_scale)))
    if (!identical(val_direct_scale, val_docall_scale)) {
      cat("   Max difference:", max(abs(val_direct_scale - val_docall_scale)), "\n")
    }
  }
}

cat("\n========================================\n")
cat("CONCLUSION\n")
cat("========================================\n")

all_identical <- (
  identical(dim(so_direct), dim(so_docall)) &&
  identical(assays_direct, assays_docall) &&
  identical(var_direct, var_docall) &&
  isTRUE(counts_equal) &&
  isTRUE(data_equal) &&
  isTRUE(scale_equal)
)

if (all_identical) {
  cat("✅ Objects are IDENTICAL\n")
  cat("The results from direct call and do.call are exactly the same.\n")
} else {
  cat("⚠️  Objects are NOT IDENTICAL\n")
  cat("There are differences between direct call and do.call results.\n")
  cat("This needs investigation!\n")
}
