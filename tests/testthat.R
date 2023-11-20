library(testthat)
library(Seurat)
#
# # Run tests for 'v3'
# message('Run tests for v3 assay')
# options(Seurat.object.assay.version = 'v3')
# test_check("Seurat")

# Run tests for 'v5'
message('Run tests for v5 assay')
options(Seurat.object.assay.version = 'v5')
test_check("Seurat")


