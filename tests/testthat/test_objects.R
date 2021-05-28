# Tests for functions in objects.R

# Tests for interacting with the meta.data slot
# ------------------------------------------------------------------------------
context("Metadata")

data("pbmc_small")

pbmc_small <- suppressWarnings(suppressMessages(UpdateSeuratObject(pbmc_small)))
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
  expect_equal(rna.assay@misc, list())
  expect_equal(GetAssayData(object = rna.assay2, slot = "counts"), new(Class = "matrix"))
})

rna.assay2 <- CreateAssayObject(counts = pbmc.raw, min.cells = 10, min.features = 30)
test_that("CreateAssayObject filtering works", {
  expect_equal(dim(x = rna.assay2), c(163, 77))
  expect_true(all(rowSums(GetAssayData(object = rna.assay2, slot = "counts")) >= 10))
  expect_true(all(colSums(GetAssayData(object = rna.assay2, slot = "counts")) >= 30))
})

test_that("CreateAssayObject catches improper input", {
  expect_error(CreateAssayObject())
  expect_error(CreateAssayObject(counts = pbmc.raw, data = pbmc.raw))
  pbmc.raw2 <- cbind(pbmc.raw[, 1:10], pbmc.raw[, 1:10])
  expect_warning(CreateAssayObject(counts = pbmc.raw2))
  expect_warning(CreateAssayObject(data = pbmc.raw2))
  pbmc.raw2 <- rbind(pbmc.raw[1:10, ], pbmc.raw[1:10, ])
  expect_warning(CreateAssayObject(counts = pbmc.raw2))
  expect_warning(CreateAssayObject(data = pbmc.raw2))
  pbmc.raw2 <- pbmc.raw
  colnames(x = pbmc.raw2) <- c()
  expect_error(CreateAssayObject(counts = pbmc.raw2))
  expect_error(CreateAssayObject(data = pbmc.raw2))
  pbmc.raw2 <- pbmc.raw
  rownames(x = pbmc.raw2) <- c()
  expect_error(CreateAssayObject(counts = pbmc.raw2))
  expect_error(CreateAssayObject(data = pbmc.raw2))
  pbmc.raw.mat <- as.matrix(x = pbmc.raw)
  pbmc.raw.df <- as.data.frame(x = pbmc.raw.mat)
  rna.assay3 <- CreateAssayObject(counts = pbmc.raw.df)
  rna.assay4 <- CreateAssayObject(counts = pbmc.raw.mat)
  expect_is(object = GetAssayData(object = rna.assay3, slot = "counts"), class = "dgCMatrix")
  expect_is(object = GetAssayData(object = rna.assay4, slot = "counts"), class = "dgCMatrix")
  pbmc.raw.underscores <- pbmc.raw
  rownames(pbmc.raw.underscores) <- gsub(pattern = "-", replacement = "_", x = rownames(pbmc.raw.underscores))
  expect_warning(CreateAssayObject(counts = pbmc.raw.underscores))
})

# Tests for creating an DimReduc object
# ------------------------------------------------------------------------------
context("CreateDimReducObject")

pca <- pbmc_small[["pca"]]
Key(object = pca) <- 'PC_'

test_that("CreateDimReducObject works", {
  pca.dr <- CreateDimReducObject(
    embeddings = Embeddings(object = pca),
    loadings = Loadings(object = pca),
    projected = Loadings(object = pca, projected = TRUE),
    assay = "RNA"
  )
  expect_equal(Embeddings(object = pca.dr), Embeddings(object = pca))
  expect_equal(Loadings(object = pca.dr), Loadings(object = pca))
  expect_equal(Loadings(object = pca.dr, projected = TRUE), Loadings(object = pca, projected = TRUE))
  expect_equal(Key(object = pca.dr), "PC_")
  expect_equal(pca.dr@assay.used, "RNA")
})

test_that("CreateDimReducObject catches improper input", {
  bad.embeddings <- Embeddings(object = pca)
  colnames(x = bad.embeddings) <- paste0("PCA", 1:ncol(x = bad.embeddings))
  expect_warning(CreateDimReducObject(embeddings = bad.embeddings, key = "PC"))
  colnames(x = bad.embeddings) <- paste0("PC", 1:ncol(x = bad.embeddings), "X")
  suppressWarnings(expect_error(CreateDimReducObject(embeddings = bad.embeddings, key = "PC")))
  suppressWarnings(expect_error(CreateDimReducObject(embeddings = bad.embeddings)))
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

# Tests for creating a Seurat object
# ------------------------------------------------------------------------------
context("Merging")

pbmc.assay <- pbmc_small[["RNA"]]
x <- merge(x = pbmc.assay, y = pbmc.assay)

test_that("Merging Assays works properly", {
  expect_equal(dim(GetAssayData(object = x, slot = "counts")), c(230, 160))
  expect_equal(dim(GetAssayData(object = x, slot = "data")), c(230, 160))
  expect_equal(GetAssayData(object = x, slot = "scale.data"), new(Class = "matrix"))
  expect_equal(Key(object = x), "rna_")
  expect_equal(VariableFeatures(object = x), vector())
  expect_equal(x[[]], data.frame(row.names = rownames(x = pbmc.assay)))
})

pbmc.assay2 <- pbmc.assay
pbmc.assay2@counts <- new("dgCMatrix")
test_that("Merging Assays handles case when counts not present", {
  y <- merge(x = pbmc.assay2, y = pbmc.assay)
  expect_equal(unname(colSums(x = GetAssayData(object = y, slot = "counts"))[1:80]), rep.int(x = 0, times = 80))
  z <- merge(x = pbmc.assay2, pbmc.assay2)
  expect_equal(nnzero(x = GetAssayData(object = z, slot = "counts")), 0)
})

pbmc.assay2 <- pbmc.assay
pbmc.assay2@data <- new("dgCMatrix")
test_that("Merging Assays handles case when data not present", {
  y <- merge(x = pbmc.assay2, y = pbmc.assay, merge.data = TRUE)
  expect_equal(unname(colSums(x = GetAssayData(object = y, slot = "data"))[1:80]), rep.int(x = 0, times = 80))
  z <- merge(x = pbmc.assay2, y = pbmc.assay2, merge.data = TRUE)
  expect_equal(nnzero(x = GetAssayData(object = z, slot = "data")), 0)
})

# Tests for Neighbor object
# ------------------------------------------------------------------------------
context("Neighbor")

# converting to Graph and back

n.rann.ob <- NNHelper(
  data = Embeddings(object = pbmc_small[["pca"]]),
  query = Embeddings(object = pbmc_small[["pca"]]),
  k = 10,
  method = "rann")

test_that("Neighbor object methods work", {
  expect_equal(dim(x = Indices(object = n.rann.ob)), c(80, 10))
  expect_equal(dim(x = n.rann.ob), c(80, 10))
  expect_equal(as.numeric(Indices(object = n.rann.ob)[1, 7]), 45, )
  expect_equal(dim(x = Distances(object = n.rann.ob)), c(80, 10))
  expect_equal(as.numeric(Distances(object = n.rann.ob)[2, 2]), 2.643759, tolerance = 1e-6)
  expect_equal(length(x = Cells(x = n.rann.ob)), 80)
  expect_equal(Cells(x = n.rann.ob)[c(1, 20, 80)], c("ATGCCAGAACGACT", "TACATCACGCTAAC", "CTTGATTGATCTTC"))
  pbmc_small[["n.ob"]] <- n.rann.ob
  pbmc_small <- RenameCells(object = pbmc_small, add.cell.id = "test")
  expect_equal(Cells(x = pbmc_small[['n.ob']])[1], c("test_ATGCCAGAACGACT"))
  expect_equal(TopNeighbors(object = n.rann.ob, cell = "ATGCCAGAACGACT", n = 5)[5], "GATATAACACGCAT")
  expect_equal(length(TopNeighbors(object = n.rann.ob, cell = "ATGCCAGAACGACT", n = 7)), 7)
  nrg <- as.Graph(x = n.rann.ob)
  expect_true(inherits(x = nrg, what = "Graph"))
  expect_equal(as.numeric(Distances(object = n.rann.ob)[2, 3]), nrg[2, Indices(object = n.rann.ob)[2, 3]])
  nro2 <- as.Neighbor(x = nrg)
  expect_true(inherits(x = nro2, what = "Neighbor"))
  expect_equal(Distances(object = n.rann.ob)[2, 3], Distances(object = nro2)[2, 3])
  expect_equal(Indices(object = n.rann.ob)[1, 6], Indices(object = nro2)[1, 6])
})

n.annoy.ob <- NNHelper(
  data = Embeddings(object = pbmc_small[["pca"]]),
  query = Embeddings(object = pbmc_small[["pca"]]),
  k = 10,
  method = "annoy",
  cache.index = TRUE)
idx.file <-  tempfile()
SaveAnnoyIndex(object = n.annoy.ob, file = idx.file)
nao2 <- LoadAnnoyIndex(object = n.annoy.ob, file = idx.file)

test_that("Saving/Loading annoy index", {
  expect_error(SaveAnnoyIndex(object = n.rann.ob, file = idx.file))
  expect_equal(head(Indices(n.annoy.ob)), head(Indices(nao2)))
  expect_equal(head(Distances(n.annoy.ob)), head(Distances(nao2)))
  expect_false(is.null(x = Index(nao2)))
})

# Tests for FetchData
# ------------------------------------------------------------------------------
context("FetchData")

# Features to test:
# able to pull cell embeddings, data, metadata
# subset of cells

test_that("Fetching a subset of cells works", {
  x <- FetchData(object = pbmc_small, cells = colnames(x = pbmc_small)[1:10], vars = rownames(x = pbmc_small)[1])
  expect_equal(rownames(x = x), colnames(x = pbmc_small)[1:10])
  random.cells <- sample(x = colnames(x = pbmc_small), size = 10)
  x <- FetchData(object = pbmc_small, cells = random.cells, vars = rownames(x = pbmc_small)[1])
  expect_equal(rownames(x = x), random.cells)
  x <- FetchData(object = pbmc_small, cells = 1:10, vars = rownames(x = pbmc_small)[1])
  expect_equal(rownames(x = x), colnames(x = pbmc_small)[1:10])
})

suppressWarnings(pbmc_small[["RNA2"]] <- pbmc_small[["RNA"]])
Key(pbmc_small[["RNA2"]]) <- "rna2_"

test_that("Fetching keyed variables works", {
  x <- FetchData(object = pbmc_small, vars = c(paste0("rna_", rownames(x = pbmc_small)[1:5]), paste0("rna2_", rownames(x = pbmc_small)[1:5])))
  expect_equal(colnames(x = x), c(paste0("rna_", rownames(x = pbmc_small)[1:5]), paste0("rna2_", rownames(x = pbmc_small)[1:5])))
  x <- FetchData(object = pbmc_small, vars = c(paste0("rna_", rownames(x = pbmc_small)[1:5]), paste0("PC_", 1:5)))
  expect_equal(colnames(x = x), c(paste0("rna_", rownames(x = pbmc_small)[1:5]), paste0("PC_", 1:5)))
})

test_that("Fetching embeddings/loadings not present returns warning or errors", {
  expect_warning(FetchData(object = pbmc_small, vars = c("PC_1", "PC_100")))
  expect_error(FetchData(object = pbmc_small, vars = "PC_100"))
})

bad.gene <- GetAssayData(object = pbmc_small[["RNA"]], slot = "data")
rownames(x = bad.gene)[1] <- paste0("rna_", rownames(x = bad.gene)[1])
pbmc_small[["RNA"]]@data <- bad.gene

# Tests for WhichCells
# ------------------------------------------------------------------------------

test_that("Specifying cells works", {
  test.cells <- Cells(x = pbmc_small)[1:10]
  expect_equal(WhichCells(object = pbmc_small, cells = test.cells), test.cells)
  expect_equal(WhichCells(object = pbmc_small, cells = test.cells, invert = TRUE), setdiff(Cells(x = pbmc_small), test.cells))
})

test_that("Specifying idents works", {
  c12 <- WhichCells(object = pbmc_small, idents = c(1, 2))
  expect_equal(length(x = c12), 44)
  expect_equal(c12[44], "CTTGATTGATCTTC")
  expect_equal(c12, WhichCells(object = pbmc_small, idents = 0, invert = TRUE))
})

test_that("downsample works", {
  expect_equal(length(x = WhichCells(object = pbmc_small, downsample = 5)), 15)
  expect_equal(length(x = WhichCells(object = pbmc_small, downsample = 100)), 80)
})

test_that("passing an expression works", {
  lyz.pos <- WhichCells(object = pbmc_small, expression = LYZ > 1)
  expect_true(all(GetAssayData(object = pbmc_small, slot = "data")["LYZ", lyz.pos] > 1))
  # multiple values in expression
  lyz.pos <- WhichCells(object = pbmc_small, expression = LYZ > 1 & groups == "g1")
  expect_equal(length(x = lyz.pos), 30)
  expect_equal(lyz.pos[30], "CTTGATTGATCTTC")
})

# Tests for small other functions
# ------------------------------------------------------------------------------
test_that("Top works", {
  dat <- Embeddings(object = pbmc_small[['pca']])[, 1, drop = FALSE]
  expect_warning(Top(data = dat, num = 1000, balanced = FALSE))
  tpc1 <- Top(data = dat, num = 20, balanced = FALSE)
  expect_equal(length(x = tpc1), 20)
  expect_equal(tpc1[1], "ACGTGATGCCATGA")
  expect_equal(tpc1[20], "GTCATACTTCGCCT")
  tpc1b <- Top(data = dat, num = 20, balanced = TRUE)
  expect_equal(length(x = tpc1b), 2)
  expect_equal(names(tpc1b), c("positive", "negative"))
  expect_equal(length(tpc1b[[1]]), 10)
  expect_equal(length(tpc1b[[2]]), 10)
  expect_equal(tpc1b[[1]][1], "GTCATACTTCGCCT")
  expect_equal(tpc1b[[1]][10], "CTTGATTGATCTTC")
  expect_equal(tpc1b[[2]][1], "ACGTGATGCCATGA")
  expect_equal(tpc1b[[2]][10], "ATTGTAGATTCCCG")
  tpc1.sub <- Top(data = dat[1:79, , drop = FALSE], num = 79, balanced = TRUE)
  expect_equal(length(tpc1.sub[[1]]), 40)
  expect_equal(length(tpc1.sub[[2]]), 39)
})


# Tests for SCE conversion
# ------------------------------------------------------------------------------
test_that("as.SingleCellExperiment works", {
  skip_on_cran()
  if (requireNamespace('SingleCellExperiment', quietly = TRUE)) {
    mat <- matrix(1:100, ncol = 10)
    colnames(mat) <- LETTERS[1:10]
    rownames(mat) <- LETTERS[1:10]
    seuratObj <- Seurat::CreateSeuratObject(mat)
    sce <- as.SingleCellExperiment(seuratObj)

    expect_equal(ncol(sce), 10)
    expect_equal(nrow(sce), 10)
    # expect_equal(length(SingleCellExperiment::altExps(sce)), 0)
    # expect_equal(SingleCellExperiment::mainExpName(sce), 'RNA')

    seuratObj <- Seurat::CreateSeuratObject(mat)
    seuratObj[['ADT']] <- CreateAssayObject(mat)
    sce <- as.SingleCellExperiment(seuratObj)
    expect_equal(ncol(sce), 10)
    expect_equal(nrow(sce), 10)
    # expect_equal(names(SingleCellExperiment::altExps(sce)), 'ADT')
    # expect_equal(SingleCellExperiment::mainExpName(sce), 'RNA')
  }
})
