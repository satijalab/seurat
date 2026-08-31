set.seed(42)

# A single layer has nothing to integrate across. Subsetting a cell type that
# happens to come from one sample produces exactly that, which is why this hits
# small subsets and not large ones.
build_object <- function(nbatch, ncell = 120L, nfeat = 120L) {
  counts <- as.sparse(matrix(rpois(nfeat * ncell, lambda = 2), nrow = nfeat))
  dimnames(counts) <- list(
    paste0("gene", seq_len(nfeat)),
    paste0("c", seq_len(ncell))
  )
  object <- suppressWarnings(CreateSeuratObject(counts = counts))
  object$batch <- rep(paste0("b", seq_len(nbatch)), length.out = ncell)
  object[["RNA"]] <- split(object[["RNA"]], f = object$batch)
  object <- suppressWarnings(NormalizeData(object, verbose = FALSE))
  object <- suppressWarnings(
    FindVariableFeatures(object, nfeatures = 60, verbose = FALSE)
  )
  object <- suppressWarnings(ScaleData(object, verbose = FALSE))
  suppressWarnings(RunPCA(object, npcs = 10, verbose = FALSE))
}

test_that("integrating a single layer explains itself", {
  object <- build_object(nbatch = 1L)
  expect_length(Layers(object, search = "data"), 1L)
  # previously: no applicable method for 'Assays' applied to an object of
  # class "NULL", from deep inside the integration method
  for (method in list(CCAIntegration, RPCAIntegration, HarmonyIntegration)) {
    expect_error(
      suppressWarnings(IntegrateLayers(
        object, method = method, orig.reduction = "pca",
        new.reduction = "integrated", verbose = FALSE
      )),
      "at least two layers"
    )
  }
})

test_that("the message names the assay and how to split it", {
  object <- build_object(nbatch = 1L)
  expect_error(
    suppressWarnings(IntegrateLayers(
      object, method = CCAIntegration, orig.reduction = "pca",
      new.reduction = "integrated", verbose = FALSE
    )),
    "split\\(object"
  )
})

test_that("two layers integrate as before", {
  object <- build_object(nbatch = 2L)
  expect_length(Layers(object, search = "data"), 2L)
  integrated <- suppressWarnings(IntegrateLayers(
    object, method = CCAIntegration, orig.reduction = "pca",
    new.reduction = "integrated.cca", dims = 1:10, k.weight = 40,
    verbose = FALSE
  ))
  expect_true("integrated.cca" %in% Reductions(integrated))
  expect_equal(nrow(Embeddings(integrated, "integrated.cca")), ncol(object))
})

test_that("the multi-umi-assay message is not mangled and names the assays", {
  skip_if_not_installed("glmGamPoi")
  set.seed(1)
  counts <- as.sparse(matrix(rpois(100 * 80, lambda = 3), nrow = 100))
  dimnames(counts) <- list(paste0("gene", 1:100), paste0("c", 1:80))
  object <- suppressWarnings(CreateSeuratObject(counts = counts))
  object <- suppressWarnings(
    SCTransform(object, method = "glmGamPoi", verbose = FALSE, ncells = 80)
  )
  assay <- object[["SCT"]]
  models <- slot(assay, "SCTModel.list")
  models <- c(models, models)
  names(models) <- c("model1", "model2")
  slot(models[[2]], "umi.assay") <- "OTHER"
  slot(assay, "SCTModel.list") <- models
  object[["SCT"]] <- assay

  message <- tryCatch(
    suppressWarnings(FindTransferAnchors(
      reference = object, query = object,
      reference.assay = "SCT", query.assay = "SCT",
      normalization.method = "SCT", verbose = FALSE
    )),
    error = function(e) conditionMessage(e)
  )
  # `call = FALSE` (no dot) was pasted into the message as "assaysFALSE"
  expect_false(grepl("FALSE", message, fixed = TRUE))
  expect_match(message, "multiple different")
  expect_match(message, "OTHER")
})
