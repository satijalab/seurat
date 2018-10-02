# Tests for functions in objects.R

# Tests for interacting with the meta.data slot
# --------------------------------------------------------------------------------
context("Metadata")

test_that("AddMetaData adds in vector properly ", {
  cluster_letters <- LETTERS[Idents(object = pbmc_small)]
  names(cluster_letters) <- colnames(x = pbmc_small)
  pbmc_small <- AddMetaData(object = pbmc_small, metadata = cluster_letters, col.name = 'letter.idents')
  expect_equal(pbmc_small$letter.idents, cluster_letters)
  cluster_letters_shuffled <- sample(x = cluster_letters)
  pbmc_small <- AddMetaData(object = pbmc_small, metadata = cluster_letters_shuffled, col.name = 'letter.idents.shuffled')
  expect_equal(pbmc_small$letter.idents, pbmc_small$letter.idents.shuffled)
})

test_that("AddMetaData adds in data frame properly", {
  cluster_letters_df <- data.frame(A=cluster_letters, B=cluster_letters_shuffled)
  pbmc_small <- AddMetaData(object = pbmc_small, metadata = cluster_letters_df)
  expect_equal(pbmc_small[[c("A", "B")]], cluster_letters_df)
})

test_that("AddMetaData errors", {
  expect_error(AddMetaData(object = pbmc_small, metadata = cluster_letters, col.name = "RNA"))
  expect_error(AddMetaData(object = pbmc_small, metadata = c(unname(cluster_letters), "A"), col.name = "letter.idents"))
})