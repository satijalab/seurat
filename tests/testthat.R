if (Sys.getenv("CI")=="true") {
  library(reticulate)
  conda_create("rettest")
  conda_install("rettest", c("python=3.6", "numpy", "seaborn", "scikit-learn", "statsmodels", "numba", "cython"), forge=F)
  conda_install("rettest", c("scanpy", "phate"), pip=T, pip_ignore_installed=F)
  use_condaenv("rettest")
}

library(testthat)
library(Seurat)

test_check("Seurat")
