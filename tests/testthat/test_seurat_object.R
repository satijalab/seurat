# Tests for functions dependent on a seurat object
set.seed(42)

# load a minimal example data set (subset of nbt dataset)
load("../testdata/nbt_small.Rdata")
nbt.small <- log(nbt.small + 1)


# Test Initial Normalization
expect_equal(LogNormalize(matrix(1:16, nrow = 4))[1,1], 6.908755, tolerance = 1e-6)

# Tests for object creation (via new/Setup)
# --------------------------------------------------------------------------------
context("Object creation")

# Generate Seurat object
min.cells <- 3
project.name <- "nbt.test"
names.field <- 2
names.delim <- "_"
min.genes <- 1000
expression.thresh <- 1

nbt.test <- new("seurat", raw.data = nbt.small)

nbt.test <- CreateSeuratObject(raw.data = nbt.small,
                               project = project.name,
                               min.cells = min.cells,
                               names.field = names.field,
                               names.delim = names.delim,
                               min.genes = min.genes,
                               is.expr = expression.thresh,
                               do.scale = T,
                               do.center = T)

test_that("object initialization creates seurat object", {
  expect_is(nbt.test, "seurat")
})

test_that("entered parameters set correctly", {
  expect_match(project.name, nbt.test@project.name)
  expect_equal(expression.thresh, nbt.test@is.expr)
})
test_that("correct cells are used",{
  gene.count <- nbt.test@meta.data$nGene
  expect_equal(min(gene.count), 2405)
  expect_true(all(gene.count >= min.genes))
})

test_that("correct genes are used", {
  useable.genes <- rowSums(nbt.test@raw.data > expression.thresh)
  useable.genes <- useable.genes[useable.genes >= min.cells]
  used.genes <- rownames(nbt.test@data)

  expect_true(length(useable.genes) > 0)
  expect_equal(length(useable.genes), length(used.genes))
})

test_that("names and IDs set correctly", {
  expect_true(length(colnames(nbt.test@raw.data)) > 0)
  expect_equal(nbt.test@cell.names, colnames(nbt.test@raw.data))

  expected.cluster.ids = c("GW21.2", "GW16", "GW21")
  expect_equal(as.vector(unique(nbt.test@ident)), expected.cluster.ids)
  expect_equal(as.vector(unique(nbt.test@ident)), as.vector(unique(nbt.test@meta.data$orig.ident)))

})

test_that("scaling done correctly", {
  expect_equal(nbt.test@scale.data["AACS", "Hi_GW21.2_3"], 1.6771640694)
  expect_equal(nbt.test@scale.data["ZYX", "Hi_GW16_1"], -0.61829233)
})

# Test dimensional reduction
# --------------------------------------------------------------------------------
context("PCA dimensional reduction")

nbt.test <- FindVariableGenes(
  nbt.test,
  y.cutoff = 2,
  x.low.cutoff = 2,
  mean.function = ExpMean,
  dispersion.function = LogVMR
)

pcs.compute <- 3
nbt.test <- RunPCA(nbt.test, pcs.compute = pcs.compute, do.print = FALSE, weight.by.var = F)

test_that("PCA returns expected data when not scaling", {
  expect_equal(abs(nbt.test@dr$pca@cell.embeddings[1,1]), 0.26627994, tolerance = 1e-6)
  expect_equal(abs(nbt.test@dr$pca@gene.loadings[1,1]), 0.5261299, tolerance = 1e-6)
  expect_equal(ncol(nbt.test@dr$pca@gene.loadings), pcs.compute)
  expect_equal(ncol(nbt.test@dr$pca@cell.embeddings), pcs.compute)

})

nbt.test <- RunPCA(nbt.test, pcs.compute = pcs.compute, do.print = FALSE)
test_that("PCA returns expected data when scaling by variance explained", {
  expect_true(nrow(nbt.test@dr$pca@cell.embeddings) == ncol(nbt.test@data))
  expect_true(nrow(nbt.test@dr$pca@gene.loadings) == length(nbt.test@var.genes))
  expect_equal(abs(nbt.test@dr$pca@cell.embeddings[1,1]), 1.423131, tolerance = 1e-6)
  expect_equal(abs(nbt.test@dr$pca@gene.loadings[1,1]), 0.5261299, tolerance = 1e-6 )
})

# Tests for tSNE
# --------------------------------------------------------------------------------
context("tSNE")
nbt.test <- RunTSNE(nbt.test, dims.use = 1:2, do.fast = T, perplexity = 4)
test_that("tSNE is run correctly", {
  expect_equal(nrow(nbt.test@dr$tsne@cell.embeddings), ncol(nbt.test@data))
})

test_that("tSNE plots correctly", {
  p <- TSNEPlot(object = nbt.test)
  expect_is(object = p, class = c('list', 'ggplot'))
  num.groups <- if (inherits(x = p, what = 'list')) {
    p[[1]][[1]]$group
  } else {
    p$data$ident
  }
  num.groups <- length(x = unique(x = num.groups))
  expect_equal(object = num.groups, expected = 3)
})

# Tests for plotting functionality (via Setup)
# --------------------------------------------------------------------------------
context("Plotting/Visualization")

test_that("Violin plots (VlnPlot() ) return as expected", {
  expect_is(VlnPlot(nbt.test, "ZYX", do.ret = T)$layers[[1]]$geom, "GeomViolin" )
  expect_equal(length(VlnPlot(nbt.test, c("ZYX", "AACS"), do.return = T, return.plotlist = T)), 2)

})

test_that("CellPlots return as expected", {
  expect_equal(CellPlot(nbt.test, nbt.test@cell.names[1], nbt.test@cell.names[2]), NULL)
})

test_that("GenePlots return as expected", {
  expect_equal(GenePlot(nbt.test,"DLX1","DLX2"), NULL)
})

test_that("FeaturePlot works as expected", {
  expect_is(FeaturePlot(nbt.test, "DLX1"), "NULL")
  plot <- FeaturePlot(nbt.test, c("DLX1", "nGene"), do.return = T)
  expect_is(plot, "list")
  expect_equal(length(plot), 2)
  expect_is(FeaturePlot(nbt.test, "DLX1", cols.use = "Purples"), "NULL")
})


# Tests for clustering related functions
# --------------------------------------------------------------------------------
context("Clustering Functions")

test_that("SNN calculations are correct and handled properly", {
  expect_true(length(nbt.test@snn) == 0)

  nbt.test <- FindClusters(nbt.test, dims.use = 1:2, print.output = 0, k.param = 4, save.SNN = T)
  expect_true(length(nbt.test@snn) > 1)
  expect_equal(nbt.test@snn[2,9], 1/3)

  nbt.test <- FindClusters(nbt.test, resolution = 1, print.output = 0)

  expect_warning(FindClusters(nbt.test, k.param = 4, reuse.SNN = T, resolution = 1, n.iter = 1, n.start = 1, print.output = 0))
  nbt.test@snn <- sparseMatrix(1, 1, x = 1)
  expect_error(FindClusters(nbt.test, resolution = 1, reuse.SNN = T))

})
nbt.test <- FindClusters(nbt.test, k.param = 4, resolution = seq(1,2,0.1), print.output = 0, n.iter = 1,
                         n.start = 1)
test_that("Clustering over multiple resolution values handled correctly", {
  expect_equal(length(nbt.test@meta.data$res.1), ncol(nbt.test@data))
  expect_equal(length(nbt.test@meta.data$res.2), ncol(nbt.test@data))
  expect_equal(length(nbt.test@snn), 1)
})

# Test subsetting functionality
# --------------------------------------------------------------------------------
context("Cell Subsetting")
test_that("WhichCells subsets properly", {
  expect_equal(length(WhichCells(nbt.test, 1)), 3)
  expect_equal(length(WhichCells(nbt.test, c(1,2))), 6)
  expect_error(WhichCells(nbt.test, 10))
  expect_equal(WhichCells(nbt.test)[1], "Hi_GW21.2_3")
  expect_equal(WhichCells(nbt.test, subset.name = "nGene", accept.high = 3000, accept.low = 2500), "Hi_GW16_23")
  expect_equal(WhichCells(nbt.test, subset.name = "PC1", accept.high = -1.4, accept.low = -1.5), "Hi_GW21.2_3")

  expect_equal(length(WhichCells(nbt.test, max.cells.per.ident = 1)), length(unique(nbt.test@ident)))
  expect_equal(length(WhichCells(nbt.test, c(1,2), max.cells.per.ident = 1)), 2)
  expect_equal(length(WhichCells(nbt.test, subset.name = "nGene", max.cells.per.ident = 1)), length(unique(nbt.test@ident)))
})

test_that("SubsetData works properly", {
  nbt.test@dr <- list()
  count <- length(WhichCells(nbt.test, 1))
  nbt.test.subset <- SubsetData(nbt.test, ident.use = 1)
  expect_equal(length(nbt.test.subset@ident), count)
})


test_that("SubsetByPredicate subsets properly", {
  expect_equal(SubsetByPredicate(nbt.test, vars.use = "nGene", predicate = "2500 < nGene & nGene < 3000")@cell.names, "Hi_GW16_23")
  expect_equal(SubsetByPredicate(nbt.test, vars.use = "PC1",   predicate = "-1.5 < PC1 & PC1 < -1.4"    )@cell.names, "Hi_GW21.2_3")

  nbt.test@dr <- list()
  count <- length(WhichCells(nbt.test, 1))
  nbt.test.subset <- SubsetByPredicate(nbt.test, "ident", "ident == '1'")
  expect_equal(length(nbt.test.subset@ident), count)

})

# Test CCA procedure
# --------------------------------------------------------------------------------
context("CCA Alignment")
scrambled.cells <- sample(nbt.test@cell.names)
c1 <- SubsetData(nbt.test, cells.use = scrambled.cells[1:7])
c2 <- SubsetData(nbt.test, cells.use = scrambled.cells[8:14])
c3 <- RunCCA(c1, c2, genes.use = c1@var.genes, num.cc = 3)
nbt.test <- SetIdent(nbt.test, cells.use = c1@cell.names, ident.use = "g1")
nbt.test <- SetIdent(nbt.test, cells.use = c2@cell.names, ident.use = "g2")
c4 <- RunCCA(nbt.test, group1 = "g1", group2 = "g2", genes.use = c1@var.genes, num.cc = 3)
c4@dr$cca@cell.embeddings <- c4@dr$cca@cell.embeddings[c3@cell.names, ]

test_that("CCA returns the expected cell.embeddings matrix values when run on two objects", {
  expect_equal(nrow(c3@dr$cca@cell.embeddings), 14)
  expect_equal(ncol(c3@dr$cca@cell.embeddings), 3)
  expect_equal(abs(unname(c3@dr$cca@cell.embeddings[1,1])), 0.3108733, tolerance = 1e-6 )
  expect_equal(abs(unname(c3@dr$cca@cell.embeddings[14,3])), 0.6064297, tolerance = 1e-6 )
})
test_that("CCA returns the expected cell.embeddings matrix values when run on single object", {
  expect_equal(nrow(c4@dr$cca@cell.embeddings), 14)
  expect_equal(ncol(c4@dr$cca@cell.embeddings), 3)
  expect_equal(abs(unname(c4@dr$cca@cell.embeddings[1,1])), 0.3108733, tolerance = 1e-6 )
  expect_equal(abs(unname(c4@dr$cca@cell.embeddings[14,3])), 0.6064297, tolerance = 1e-6 )
})

