set.seed(42)

# A feature named twice, which happens when marker lists are combined, became a
# duplicated factor level; cells with no identity came back from
# CellsByIdentities as NA and reached FetchData
data("pbmc_small", package = "SeuratObject", envir = environment())
features <- rownames(pbmc_small)[1:3]

test_that("a feature named twice is plotted once", {
  plot <- DotPlot(pbmc_small, features = c(features, features[1]))
  expect_s3_class(plot, "ggplot")
  expect_setequal(levels(plot$data$features.plot), features)
})

test_that("a feature repeated across named groups is plotted", {
  plot <- DotPlot(pbmc_small, features = list(a = features, b = features[1]))
  expect_s3_class(plot, "ggplot")
})

test_that("cells with no identity are left out and reported", {
  object <- pbmc_small
  identities <- as.character(Idents(object))
  identities[1:3] <- NA
  Idents(object) <- factor(identities)

  plot <- expect_warning(
    DotPlot(object, features = features),
    "3 cells have no identity"
  )
  expect_s3_class(plot, "ggplot")
  expect_false(anyNA(as.character(plot$data$id)))
  expect_setequal(as.character(unique(plot$data$id)), setdiff(unique(identities), NA))
})

test_that("the usual case is unchanged", {
  plot <- DotPlot(pbmc_small, features = features)
  expect_setequal(levels(plot$data$features.plot), features)
  expect_setequal(as.character(unique(plot$data$id)), levels(Idents(pbmc_small)))
})
