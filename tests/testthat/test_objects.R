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

# Tests for creating an DimReduc object
# ------------------------------------------------------------------------------
context("CreateDimReducObject")

pca <- pbmc_small[["pca"]]

test_that("CreateDimReducObject works", {
  pca.dr <- CreateDimReducObject(
    embeddings = Embeddings(object = pca),
    loadings = Loadings(object = pca),
    projected = Loadings(object = pca, projected = TRUE),
    key = "PC",
    assay = "RNA"
  )
  expect_equal(Embeddings(object = pca.dr), Embeddings(object = pca))
  expect_equal(Loadings(object = pca.dr), Loadings(object = pca))
  expect_equal(Loadings(object = pca.dr, projected = TRUE), Loadings(object = pca, projected = TRUE))
  expect_equal(Key(object = pca.dr), "PC")
  expect_equal(pca.dr@assay.used, "RNA")
})

test_that("CreateDimReducObject catches improper input", {
  bad.embeddings <- Embeddings(object = pca)
  colnames(x = bad.embeddings) <- paste0("PCA", 1:ncol(x = bad.embeddings))
  expect_error(CreateDimReducObject(embeddings = bad.embeddings, key = "PC"))
  colnames(x = bad.embeddings) <- paste0("PC", 1:ncol(x = bad.embeddings), "X")
  expect_error(CreateDimReducObject(embeddings = bad.embeddings, key = "PC"))
  expect_error(CreateDimReducObject(embeddings = bad.embeddings))
})

# Tests for creating a Seurat object
# ------------------------------------------------------------------------------
context("CreateSeuratObject")

colnames(x = pbmc.raw) <- paste0(colnames(x = pbmc.raw), "-", pbmc_small$groups)
metadata.test <- pbmc_small[[]][, 5:7]
rownames(x = metadata.test) <- colnames(x = pbmc.raw)

test_that("CreateSeuratObject works", {
  seurat.object <- CreateSeuratObject(
    counts = pbmc.raw,
    project = "TESTING",
    assay = "RNA.TEST",
    names.field = 2,
    names.delim = "-",
    meta.data = metadata.test
  )
  expect_equal(seurat.object[[]][, 4:6], metadata.test)
  expect_equal(seurat.object@project.name, "TESTING")
  expect_equal(names(x = seurat.object), "RNA.TEST")
  expect_equal(as.vector(x = unname(obj = Idents(object = seurat.object))), unname(pbmc_small$groups))
})

test_that("CreateSeuratObject handles bad names.field/names.delim", {
  expect_warning(seurat.object <- CreateSeuratObject(
    counts = pbmc.raw[1:5,1:5],
    names.field = 3,
    names.delim = ":",
    meta.data = metadata.test
  ))
})
