set.seed(42)

# The function that finds anchors for a pair used to be a closure created
# inside FindIntegrationAnchors(), so `future` exported every object being
# integrated to every worker, whichever pair that worker was given
make_object <- function(tag, ncell = 120L, nfeat = 200L) {
  counts <- as.sparse(matrix(
    rpois(nfeat * ncell, lambda = 3),
    nrow = nfeat,
    dimnames = list(paste0("g", seq_len(nfeat)), paste0(tag, seq_len(ncell)))
  ))
  object <- suppressWarnings(CreateSeuratObject(counts = counts))
  object <- suppressWarnings(NormalizeData(object, verbose = FALSE))
  object <- suppressWarnings(FindVariableFeatures(object, nfeatures = 100, verbose = FALSE))
  suppressWarnings(ScaleData(object, verbose = FALSE))
}

anchors_for <- function(objects, features) {
  set.seed(3)
  slot(
    suppressWarnings(FindIntegrationAnchors(
      objects,
      anchor.features = features,
      dims = 1:10,
      k.filter = 50,
      verbose = FALSE
    )),
    name = "anchors"
  )
}

test_that("parallel anchors match sequential ones", {
  skip_on_os("windows")
  skip_if_not(future::supportsMulticore())
  objects <- list(make_object("a"), make_object("b"), make_object("c"))
  features <- suppressWarnings(SelectIntegrationFeatures(objects, nfeatures = 80, verbose = FALSE))

  future::plan("sequential")
  sequential <- anchors_for(objects, features)

  old.plan <- future::plan()
  on.exit(future::plan(old.plan), add = TRUE)
  future::plan("multicore", workers = 2)
  parallel <- anchors_for(objects, features)
  future::plan("sequential")

  expect_equal(parallel, sequential)
})

test_that("the work sent to a worker is a package function", {
  # a closure would carry the frame it was made in, and that frame holds every
  # object being integrated
  expect_identical(environment(Seurat:::.AnchorsForPair), asNamespace("Seurat"))
})
