#!/usr/bin/env Rscript

# Multiple timing tests to understand the performance difference

cat("Multiple timing tests for SCTransform\n")
cat("======================================\n\n")

suppressPackageStartupMessages({
  devtools::load_all(".", quiet = TRUE)
})

n_runs <- 5

direct_times <- numeric(n_runs)
docall_times <- numeric(n_runs)

for (i in 1:n_runs) {
  cat(sprintf("Run %d/%d:\n", i, n_runs))
  
  # Direct call
  data("pbmc_small")
  set.seed(42 + i)
  start_time <- Sys.time()
  so_direct <- SCTransform(pbmc_small, verbose = FALSE, vst.flavor = "v2", seed.use = 42 + i)
  direct_times[i] <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  cat(sprintf("  Direct: %.2f seconds\n", direct_times[i]))
  
  # do.call
  data("pbmc_small")
  set.seed(42 + i)
  start_time <- Sys.time()
  so_docall <- do.call(SCTransform, list(object = pbmc_small, verbose = FALSE, vst.flavor = "v2", seed.use = 42 + i))
  docall_times[i] <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  cat(sprintf("  do.call: %.2f seconds\n", docall_times[i]))
  
  cat(sprintf("  Ratio: %.2fx\n\n", docall_times[i] / direct_times[i]))
}

cat("======================================\n")
cat("Summary Statistics\n")
cat("======================================\n\n")

cat("Direct call:\n")
cat(sprintf("  Mean: %.2f seconds\n", mean(direct_times)))
cat(sprintf("  Median: %.2f seconds\n", median(direct_times)))
cat(sprintf("  Range: %.2f - %.2f seconds\n", min(direct_times), max(direct_times)))

cat("\ndo.call:\n")
cat(sprintf("  Mean: %.2f seconds\n", mean(docall_times)))
cat(sprintf("  Median: %.2f seconds\n", median(docall_times)))
cat(sprintf("  Range: %.2f - %.2f seconds\n", min(docall_times), max(docall_times)))

cat("\nRatio (do.call / direct):\n")
ratios <- docall_times / direct_times
cat(sprintf("  Mean: %.2fx\n", mean(ratios)))
cat(sprintf("  Median: %.2fx\n", median(ratios)))
cat(sprintf("  Range: %.2fx - %.2fx\n", min(ratios), max(ratios)))

cat("\n======================================\n")
cat("Analysis\n")
cat("======================================\n\n")

if (mean(ratios) < 0.8) {
  cat("⚠️  do.call is consistently FASTER than direct call\n")
  cat("This is unexpected and may indicate:\n")
  cat("  - Different promise evaluation timing\n")
  cat("  - First call loads/compiles code, second call benefits\n")
  cat("  - Need to test with alternating order\n")
} else if (mean(ratios) > 1.2) {
  cat("❌ do.call is consistently SLOWER than direct call\n")
  cat("This would indicate a problem with the fix.\n")
} else {
  cat("✅ Performance is comparable (within 20%)\n")
  cat("Both methods have similar performance.\n")
}

cat("\nLet's test with REVERSED order (do.call first):\n")
cat("----------------------------------------------\n")

# Test with reversed order
cat("\nRun with do.call FIRST:\n")
data("pbmc_small")
set.seed(99)
start_time <- Sys.time()
so_docall_first <- do.call(SCTransform, list(object = pbmc_small, verbose = FALSE, vst.flavor = "v2", seed.use = 99))
docall_first_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
cat(sprintf("  do.call (first): %.2f seconds\n", docall_first_time))

data("pbmc_small")
set.seed(99)
start_time <- Sys.time()
so_direct_second <- SCTransform(pbmc_small, verbose = FALSE, vst.flavor = "v2", seed.use = 99)
direct_second_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
cat(sprintf("  Direct (second): %.2f seconds\n", direct_second_time))
cat(sprintf("  Ratio: %.2fx\n", docall_first_time / direct_second_time))

if (docall_first_time > direct_second_time * 1.2) {
  cat("\n⚠️  When do.call goes first, it's slower!\n")
  cat("This suggests the first call pays a performance penalty.\n")
} else if (direct_second_time > docall_first_time * 1.2) {
  cat("\n⚠️  When direct goes second, it's slower!\n")
  cat("This also suggests order matters.\n")
} else {
  cat("\n✅ Order doesn't significantly affect results.\n")
}
