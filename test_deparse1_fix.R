#!/usr/bin/env Rscript

# Test the corrected Parenting function with deparse1

cat("Testing corrected Parenting function with deparse1\n\n")

devtools::load_all(".", quiet = TRUE)

# Test 1: Simple function name
cat("Test 1: Simple function call\n")
test_func1 <- function() {
  child_func1()
}
child_func1 <- function() {
  calls <- vapply(sys.calls(), function(call) {
    if (is.call(call) && length(call) > 0) {
      func <- call[[1]]
      if (is.name(func)) {
        return(as.character(func))
      } else if (is.call(func)) {
        return(deparse1(func))
      }
    }
    return("")
  }, character(1))
  cat("  Calls extracted:", paste(calls, collapse=", "), "\n")
  return(calls)
}
result1 <- test_func1()
cat("\n")

# Test 2: Package-qualified function (simulated)
cat("Test 2: Simulated package-qualified call\n")
test_func2 <- function() {
  # Simulate what would happen with base::mean
  call_sim <- quote(base::mean(x))
  func <- call_sim[[1]]
  
  if (is.call(func)) {
    cat("  func is a call to:", as.character(func)[1], "\n")
    cat("  as.character(func):", paste(as.character(func), collapse=", "), "\n")
    cat("  deparse1(func):", deparse1(func), "\n")
    
    # Test the extraction
    result_deparse <- deparse1(func)
    cat("  Extracted name:", result_deparse, "\n")
    
    # Can we grep for 'mean'?
    cat("  grep('mean', result):", grepl("mean", result_deparse), "\n")
  }
}
test_func2()
cat("\n")

# Test 3: Actual SCTransform test
cat("Test 3: SCTransform with do.call\n")
data("pbmc_small")

cat("  Running direct call...\n")
start_time <- Sys.time()
pbmc_direct <- SCTransform(pbmc_small, verbose = FALSE, vst.flavor = "v2", seed.use = 42)
direct_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
cat(sprintf("  Completed in %.2f seconds\n", direct_time))

cat("  Running do.call...\n")
data("pbmc_small")
start_time <- Sys.time()
pbmc_docall <- do.call(SCTransform, list(object = pbmc_small, verbose = FALSE, vst.flavor = "v2", seed.use = 42))
docall_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
cat(sprintf("  Completed in %.2f seconds\n", docall_time))

cat(sprintf("\n  Performance ratio: %.2fx\n", docall_time / direct_time))

if (docall_time < direct_time * 2) {
  cat("\n✅ SUCCESS: do.call performance is good!\n")
} else {
  cat("\n❌ FAILURE: do.call is still slow!\n")
}

# Test 4: Verify results are consistent
cat("\nTest 4: Result consistency\n")
direct_genes <- VariableFeatures(pbmc_direct)
docall_genes <- VariableFeatures(pbmc_docall)
overlap <- length(intersect(direct_genes, docall_genes))
cat(sprintf("  Variable features overlap: %.1f%%\n", 100 * overlap / length(direct_genes)))

if (overlap / length(direct_genes) > 0.9) {
  cat("✅ Results are consistent\n")
} else {
  cat("❌ Results differ\n")
}
