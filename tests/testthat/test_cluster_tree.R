set.seed(42)

# hclust() needs two things to join, and said so in its own terms:
# "must have n >= 2 objects to cluster", which does not mention identities
data("pbmc_small", package = "SeuratObject", envir = environment())

test_that("one identity is reported as such", {
  one <- subset(pbmc_small, idents = 0)
  expect_error(
    suppressWarnings(BuildClusterTree(one, verbose = FALSE)),
    "Cannot build a tree from 1 identity"
  )
  expect_error(
    suppressWarnings(BuildClusterTree(one, verbose = FALSE)),
    "at least two are needed"
  )
})

test_that("two identities still build a tree", {
  skip_if_not_installed("ape")
  two <- subset(pbmc_small, idents = c(0, 1))
  built <- suppressWarnings(BuildClusterTree(two, verbose = FALSE))
  expect_s3_class(Tool(built, slot = "BuildClusterTree"), "phylo")
})
