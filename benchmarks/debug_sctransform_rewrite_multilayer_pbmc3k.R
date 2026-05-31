#!/usr/bin/env Rscript

Sys.setenv(
  OMP_NUM_THREADS = 1,
  OPENBLAS_NUM_THREADS = 1,
  MKL_NUM_THREADS = 1,
  VECLIB_MAXIMUM_THREADS = 1
)
options(future.globals.maxSize = 8 * 1024^3, mc.cores = 1)
if (requireNamespace("future", quietly = TRUE)) {
  future::plan(future::sequential)
}

pkgload::load_all(
  path = "/Users/rahulsatija/optimizer/seurat-optim",
  quiet = TRUE,
  export_all = FALSE,
  helpers = FALSE
)
library(SeuratData)

obj <- SeuratData::LoadData("pbmc3k")
if (is.list(obj) && length(obj) == 1L) {
  obj <- obj[[1L]]
}
obj <- UpdateSeuratObject(obj)

split_group <- rep(c("split1", "split2"), length.out = ncol(x = obj))
names(x = split_group) <- colnames(x = obj)
obj[["RNA"]] <- split(obj[["RNA"]], f = split_group)

cat("Input count layers:\n")
print(Layers(object = obj[["RNA"]], search = "counts"))

cat("\nInput layer dimensions:\n")
print(vapply(
  X = Layers(object = obj[["RNA"]], search = "counts"),
  FUN = function(layer_name) {
    paste(dim(x = LayerData(obj[["RNA"]], layer = layer_name)), collapse = "x")
  },
  FUN.VALUE = character(length = 1L)
))

cat("\nRunning SCTransform_rewrite...\n")
timing <- system.time({
  result <- tryCatch(
    expr = SCTransform_rewrite(
      object = obj,
      assay = "RNA",
      new.assay.name = "SCT",
      verbose = TRUE
    ),
    error = function(e) {
      cat("\nERROR:\n")
      print(e)
      cat("\nTraceback:\n")
      traceback()
      return(e)
    }
  )
})

cat("\nTiming:\n")
print(timing)

if (inherits(x = result, what = "error")) {
  quit(save = "no", status = 1)
}

scale.data <- GetAssayData(object = result[["SCT"]], layer = "scale.data")

cat("\nOutput SCT summary:\n")
cat("scale.data dim:", paste(dim(x = scale.data), collapse = "x"), "\n")
cat("variable features:", length(x = VariableFeatures(object = result[["SCT"]])), "\n")
cat("any NA in scale.data:", anyNA(x = scale.data), "\n")
cat("models:", paste(levels(x = result[["SCT"]]), collapse = ", "), "\n")

cat("\nModel feature counts:\n")
print(vapply(
  X = levels(x = result[["SCT"]]),
  FUN = function(model_name) {
    nrow(x = SCTResults(object = result[["SCT"]], slot = "feature.attributes", model = model_name))
  },
  FUN.VALUE = integer(length = 1L)
))
