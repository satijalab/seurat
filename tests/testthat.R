if (Sys.getenv("CI")=="true") {
  library(reticulate)
  conda_create("seurat_test_env")
  conda_install("seurat_test_env", c("python=3.6", "numpy", "seaborn", "scikit-learn", "statsmodels", "numba", "cython"), forge=F)
  if (Sys.info()[["sysname"]] == "Linux") {
    conda_install("seurat_test_env", c("scanpy", "phate"), pip=T)
  } else if (Sys.info()[["sysname"]] == "Darwin") {
    conda_install("seurat_test_env", c("scanpy", "phate"), pip=T, pip_ignore_installed=F)
  }
  use_condaenv("seurat_test_env")
}

library(testthat)
library(Seurat)

test_check("Seurat")
