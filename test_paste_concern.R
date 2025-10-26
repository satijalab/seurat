#!/usr/bin/env Rscript

# Test case to verify the PR comment about paste() with collapse=""

cat("Testing the concern raised in PR comment #2463898138\n\n")

# Test what happens with different call structures
test_cases <- list(
  # Normal function
  quote(my_function()),
  # Package qualified function
  quote(pkg::myfunc()),
  # Subsetting operator
  quote(base::`[[`(x, 1)),
  # Double colon operator itself
  quote(`::`(pkg, func)),
  # Other operators
  quote(`+`(1, 2)),
  quote(`[`(x, 1))
)

cat("Testing current implementation (collapse=''):\n")
cat("================================================\n\n")

for (i in seq_along(test_cases)) {
  call <- test_cases[[i]]
  cat(sprintf("Test %d: %s\n", i, deparse(call)))
  
  if (is.call(call) && length(call) > 0) {
    func <- call[[1]]
    
    if (is.name(func)) {
      result <- as.character(func)
      cat(sprintf("  is.name: TRUE\n"))
      cat(sprintf("  Result: '%s'\n", result))
    } else if (is.call(func)) {
      cat(sprintf("  is.call: TRUE\n"))
      cat(sprintf("  as.character(func): %s\n", paste(shQuote(as.character(func)), collapse=", ")))
      
      # Current implementation
      result_current <- paste(as.character(func), collapse = "")
      cat(sprintf("  collapse='': '%s'\n", result_current))
      
      # Suggested implementation
      result_suggested <- paste(as.character(func), collapse = "::")
      cat(sprintf("  collapse='::': '%s'\n", result_suggested))
      
      # Alternative: deparse1
      result_deparse <- deparse1(func)
      cat(sprintf("  deparse1: '%s'\n", result_deparse))
    }
  }
  cat("\n")
}

cat("\n================================================\n")
cat("Testing with actual call stack scenarios:\n")
cat("================================================\n\n")

# Create a test function that uses :: operator
test_with_qualified <- function() {
  # Simulate what sys.calls() would return
  calls <- list(
    quote(base::mean(x)),
    quote(SeuratObject::GetAssayData(object)),
    quote(object[[assay]]),
    quote(`::`(pkg, func))
  )
  
  for (call in calls) {
    cat(sprintf("Call: %s\n", deparse(call)))
    if (is.call(call) && length(call) > 0) {
      func <- call[[1]]
      if (is.call(func)) {
        cat(sprintf("  Current (collapse=''): '%s'\n", 
                    paste(as.character(func), collapse = "")))
        cat(sprintf("  Suggested (collapse='::'): '%s'\n", 
                    paste(as.character(func), collapse = "::")))
        cat(sprintf("  deparse1: '%s'\n", deparse1(func)))
      }
    }
    cat("\n")
  }
}

test_with_qualified()

cat("\n================================================\n")
cat("Analysis:\n")
cat("================================================\n\n")

cat("The concern in the PR comment is VALID if:\n")
cat("- We have operators like `[[` that get represented as multi-element vectors\n")
cat("- Using collapse='' would incorrectly concatenate them\n\n")

cat("However, let's verify what actually happens with :: operator:\n\n")

func_colons <- quote(`::`(pkg, func))[[1]]
cat(sprintf("For `::` operator:\n"))
cat(sprintf("  as.character: %s\n", paste(shQuote(as.character(func_colons)), collapse=", ")))
cat(sprintf("  length: %d\n", length(as.character(func_colons))))
cat(sprintf("  Result with collapse='': '%s'\n", paste(as.character(func_colons), collapse="")))
cat(sprintf("  Result with collapse='::': '%s'\n", paste(as.character(func_colons), collapse="::")))

cat("\n")
func_bracket <- quote(`[[`(x, 1))[[1]]
cat(sprintf("For `[[` operator:\n"))
cat(sprintf("  as.character: %s\n", paste(shQuote(as.character(func_bracket)), collapse=", ")))
cat(sprintf("  length: %d\n", length(as.character(func_bracket))))
cat(sprintf("  Result: '%s'\n", paste(as.character(func_bracket), collapse="")))
