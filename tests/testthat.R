if (Sys.getenv("TRAVIS")=="true") {
  library(reticulate)
  use_condaenv("seurat_test_env")
}

library(testthat)
library(Seurat)

test_check("Seurat")
