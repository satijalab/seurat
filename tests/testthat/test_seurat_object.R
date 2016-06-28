# Tests for functions dependent on a seurat object

# load a minimal example data set (subset of nbt dataset)
load("../testdata/nbt_small.Rdata")


# Tests for object creation (via setup)
# --------------------------------------------------------------------------------
context("Object creation")

# Generate Seurat object
min_cells <- 3
project_name <- "NBT_TEST"
names_field <- 2
names_delim <- "_"
min_genes <- 1000
expression_thresh <- 1

nbt_test <- new("seurat", raw.data = nbt_small)

test_that("object initialization creates seurat object", {
  expect_is(nbt_test, "seurat")
})

nbt_test <- setup(nbt_test, project = project_name, min.cells = min_cells, names.field = names_field,  
                  names.delim = names_delim, min.genes = min_genes, is.expr = expression_thresh, 
                  large.object = T ) 

test_that("entered parameters set correctly", {
  expect_match(project_name, nbt_test@project.name)
  expect_equal(expression_thresh, nbt_test@is.expr)
  
  
})

test_that("correct cells are used",{
  gene_count <- unname(findNGene(nbt_test@raw.data, nbt_test@is.expr))
  expect_equal(min(gene_count), 2814)
  expect_true(all(gene_count >= min_genes))
})

test_that("correct genes are used", {
  usuable_genes <- rowSums(nbt_test@raw.data > expression_thresh)
  usuable_genes <- usuable_genes[usuable_genes >= min_cells]
  used_genes <- rownames(nbt_test@data)
  
  expect_true(length(usuable_genes) > 0)
  expect_equal(length(usuable_genes), length(used_genes))
})

test_that("names and IDs set correctly", {
  expect_true(length(colnames(nbt_test@raw.data)) > 0)
  expect_equal(nbt_test@cell.names, colnames(nbt_test@raw.data))
  
  expected_cluster_ids = c("GW21.2", "GW16", "GW21")
  expect_equal(as.vector(unique(nbt_test@ident)), expected_cluster_ids)
  expect_equal(as.vector(unique(nbt_test@ident)), as.vector(unique(nbt_test@data.info$orig.ident)))
  
})

test_that("scaling done correctly", {
  expect_equal(nbt_test@scale.data["AAAS", "Hi_GW21.2_3"], 3.28266251317083)
  expect_equal(nbt_test@scale.data["ZYX", "Hi_GW16_1"], -0.380777117739444)
})

test_that("nGene calculations are consistent" , {
  gene_count <- unname(findNGene(nbt_test@raw.data, nbt_test@is.expr))
  expect_equal(nbt_test@mix.probs[, 1], gene_count)
  expect_equal(nbt_test@gene.scores[, 1], gene_count)
  
})
