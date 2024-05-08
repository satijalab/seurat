PATH_TO_DATA <- system.file("extdata", "pbmc_raw.txt", package = "Seurat")

TOLERANCE = 1.0e-4

setup_test_data <- function(cluster_resolution) {
    counts <- as.sparse(
        as.matrix(read.table(PATH_TO_DATA, sep = "\t", row.names = 1))
    )
    test.data <- CreateSeuratObject(counts)
 
    test.data <- NormalizeData(test.data)
    test.data <- FindVariableFeatures(test.data)
    test.data <- ScaleData(test.data)

    test.data <- RunPCA(test.data, npcs = 20)
    test.data <- FindNeighbors(test.data, dims = 1:20)
    test.data <- FindClusters(test.data, resolution = cluster_resolution)
    
    return(test.data)
}

context("BuildClusterTree")

test_that("BuildClusterTree works as expected", {
    skip_if_not_installed("ape")

    # TODO: this test doesn't really check that the phylogenic trees being
    # generated are biologically valid but that would be more useful
    # for testing purposes we'll over-cluster our dataset to introduce
    # some structure into the resulting phylogeny
    test_case <- setup_test_data(cluster_resolution = 5)

    result <- BuildClusterTree(test_case, slot = "counts")
    tree <- Tool(result, slot = "BuildClusterTree")
    # check that the tree contains the expected number of leaf nodes
    expect_equal(length(tree$tip.label), 6)
    # check that the tree contains the expected number of edges
    expect_equal(length(tree$edge), 20)
    # spot check a few edge lengths
    expect_equal(tree$edge.length[[1]], 1127.60199, tolerance = TOLERANCE)
    expect_equal(tree$edge.length[[5]], 182.42944, tolerance = TOLERANCE)
    expect_equal(tree$edge.length[[10]], 59.24525, tolerance = TOLERANCE)

    result <- BuildClusterTree(test_case, slot = "data")
    tree <- Tool(result, slot = "BuildClusterTree")
    # check that the tree contains the expected number of leaf nodes
    expect_equal(length(tree$tip.label), 6)
    # check that the tree contains the expected number of edges
    expect_equal(length(tree$edge), 20)
    # spot check a few edge lengths
    expect_equal(tree$edge.length[[1]], 1206.4312, tolerance = TOLERANCE)
    expect_equal(tree$edge.length[[5]], 657.0097, tolerance = TOLERANCE)
    expect_equal(tree$edge.length[[10]], 982.9843, tolerance = TOLERANCE)
})
