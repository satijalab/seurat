# Tests for functions in objects.R

# Tests for interacting with the meta.data slot
# ------------------------------------------------------------------------------
context("Metadata")

cluster_letters <- LETTERS[Idents(object = pbmc_small)]
names(cluster_letters) <- colnames(x = pbmc_small)
cluster_letters_shuffled <- sample(x = cluster_letters)

test_that("AddMetaData adds in cell-level vector properly ", {
  pbmc_small <- AddMetaData(object = pbmc_small, metadata = cluster_letters, col.name = 'letter.idents')
  expect_equal(pbmc_small$letter.idents, cluster_letters)
  pbmc_small <- AddMetaData(object = pbmc_small, metadata = cluster_letters_shuffled, col.name = 'letter.idents.shuffled')
  expect_equal(pbmc_small$letter.idents, pbmc_small$letter.idents.shuffled)
})

cluster_letters_df <- data.frame(A = cluster_letters, B = cluster_letters_shuffled)
test_that("AddMetaData adds in data frame properly for cell-level metadata", {
  pbmc_small <- AddMetaData(object = pbmc_small, metadata = cluster_letters_df)
  expect_equal(pbmc_small[[c("A", "B")]], cluster_letters_df)
})

feature_letters <- sample(x = LETTERS, size = nrow(x = pbmc_small[["RNA"]]), replace = TRUE)
names(feature_letters) <- rownames(x = pbmc_small[["RNA"]])
feature_letters_shuffled <- sample(x = feature_letters)

test_that("AddMetaData adds feature level metadata", {
  pbmc_small[["RNA"]] <- AddMetaData(object = pbmc_small[["RNA"]], metadata = feature_letters, col.name = 'feature_letters')
  expect_equal(pbmc_small[["RNA"]][["feature_letters", drop = TRUE]], feature_letters)
  pbmc_small[["RNA"]] <- AddMetaData(object = pbmc_small[["RNA"]], metadata = feature_letters_shuffled, col.name = 'feature_letters_shuffled')
  expect_equal(pbmc_small[["RNA"]][["feature_letters", drop = TRUE]], pbmc_small[["RNA"]][["feature_letters_shuffled", drop = TRUE]])
})

feature_letters_df <- data.frame(A = feature_letters, B = feature_letters_shuffled)
test_that("AddMetaData adds in data frame properly for Assays", {
  pbmc_small[["RNA"]] <- AddMetaData(object = pbmc_small[["RNA"]], metadata = feature_letters_df)
  expect_equal(pbmc_small[["RNA"]][[c("A", "B")]], feature_letters_df)
})

test_that("AddMetaData errors", {
  expect_error(AddMetaData(object = pbmc_small, metadata = cluster_letters, col.name = "RNA"))
  expect_error(AddMetaData(object = pbmc_small, metadata = c(unname(cluster_letters), "A"), col.name = "letter.idents"))
  expect_error(AddMetaData(object = pbmc_small, metadata = feature_letters, col.name = "letter.idents"))
  expect_error(AddMetaData(object = pbmc_small[["RNA"]], metadata = cluster_letters, col.name = "letter.idents"))
})

# Tests for creating an Assay object
# ------------------------------------------------------------------------------
context("CreateAssayObject")

pbmc.raw <- GetAssayData(object = pbmc_small[["RNA"]], slot = "counts")
rna.assay <- CreateAssayObject(counts = pbmc.raw)
rna.assay2 <- CreateAssayObject(data = pbmc.raw)

test_that("CreateAssayObject works as expected", {
  expect_equal(dim(x = rna.assay), c(230, 80))
  expect_equal(rownames(x = rna.assay), rownames(x = pbmc.raw))
  expect_equal(colnames(x = rna.assay), colnames(x = pbmc.raw))
  expect_equal(GetAssayData(object = rna.assay, slot = "counts"), pbmc.raw)
  expect_equal(GetAssayData(object = rna.assay, slot = "data"), pbmc.raw)
  expect_equal(GetAssayData(object = rna.assay, slot = "scale.data"), new(Class = "matrix"))
  expect_equal(dim(rna.assay[[]]), c(230, 0))
  expect_equal(rownames(x = rna.assay[[]]), rownames(x = rna.assay))
  expect_equal(VariableFeatures(object = rna.assay), vector())
  expect_equal(rna.assay@misc, NULL)
  expect_equal(GetAssayData(object = rna.assay2, slot = "counts"), new(Class = "matrix"))
})

rna.assay2 <- CreateAssayObject(counts = pbmc.raw, min.cells = 10, min.features = 30)
test_that("CreateAssayObject filtering works", {
  expect_equal(dim(x = rna.assay2), c(162, 75))
  expect_true(all(rowSums(GetAssayData(object = rna.assay2, slot = "counts")) >= 10))
  expect_true(all(colSums(GetAssayData(object = rna.assay2, slot = "counts")) >= 30))
})

test_that("CreateAssayObject catches improper input", {
  expect_error(CreateAssayObject())
  expect_error(CreateAssayObject(counts = pbmc.raw, data = pbmc.raw))
  pbmc.raw2 <- cbind(pbmc.raw[, 1:10], pbmc.raw[, 1:10])
  expect_error(CreateAssayObject(counts = pbmc.raw2))
  expect_error(CreateAssayObject(data = pbmc.raw2))
  pbmc.raw2 <- rbind(pbmc.raw[1:10, ], pbmc.raw[1:10, ])
  expect_error(CreateAssayObject(counts = pbmc.raw2))
  expect_error(CreateAssayObject(data = pbmc.raw2))
  pbmc.raw2 <- pbmc.raw
  colnames(x = pbmc.raw2) <- c()
  expect_error(CreateAssayObject(counts = pbmc.raw2))
  expect_error(CreateAssayObject(data = pbmc.raw2))
  pbmc.raw2 <- pbmc.raw
  rownames(x = pbmc.raw2) <- c()
  expect_error(CreateAssayObject(counts = pbmc.raw2))
  expect_error(CreateAssayObject(data = pbmc.raw2))
})
