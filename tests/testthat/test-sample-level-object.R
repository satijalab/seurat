# Ported from scSLIDE. Uses testthat edition 3 (local_mocked_bindings etc.).
testthat::local_edition(3)

# Tests for GenerateSampleObject() argument handling and matrix construction.

# Helper: build a minimal Seurat object carrying a LANDMARK assay and a Neighbor
# object, i.e. the shape PrepareSampleObject() hands to GenerateSampleObject().
# The first n.landmark cells are the landmarks; every cell gets k neighbours.
make_sample_fixture <- function(n.genes = 50, n.cells = 60, n.landmark = 8,
                                k = 10, seed = 42) {
  set.seed(seed)
  cts <- matrix(
    rpois(n = n.genes * n.cells, lambda = 5),
    nrow = n.genes, ncol = n.cells,
    dimnames = list(paste0("g", seq_len(n.genes)), paste0("cell", seq_len(n.cells)))
  )
  obj <- suppressWarnings(CreateSeuratObject(counts = cts))
  obj$donor     <- rep(paste0("D", 1:6), each = n.cells / 6)
  obj$sex       <- rep(c("M", "F"), each = n.cells / 2)
  obj$region    <- rep(c("R1", "R2"), each = n.cells / 2)
  obj$batch     <- rep(c("b1", "b2"), each = n.cells / 2)
  obj$replicate <- rep(c("r1", "r2"), times = n.cells / 2)
  obj$only_one  <- "SAME"
  obj[["LANDMARK"]] <- suppressWarnings(
    CreateAssayObject(counts = cts[, seq_len(n.landmark), drop = FALSE])
  )
  obj[["weighted.nn"]] <- methods::new(
    Class = "Neighbor",
    nn.idx = t(sapply(
      X = seq_len(n.cells),
      FUN = function(i) sample(x = seq_len(n.landmark), size = k, replace = TRUE)
    )),
    nn.dist = matrix(runif(n = n.cells * k), nrow = n.cells, ncol = k),
    cell.names = colnames(obj),
    alg.info = list()
  )
  obj
}

# Helper: call GenerateSampleObject() quietly on the fixture's NN object. The
# normalization step prints to stdout regardless of verbose, hence capture.output().
gen <- function(object, ...) {
  invisible(capture.output(
    res <- suppressWarnings(suppressMessages(
      GenerateSampleObject(object = object, nn.name = "weighted.nn", ...)
    ))
  ))
  res
}


test_that("an explicitly named but absent Neighbor object errors", {
  obj <- make_sample_fixture()

  expect_error(
    suppressMessages(GenerateSampleObject(object = obj, nn.name = "absent",
                                          group.by = "donor")),
    "correct name of the NN object"
  )
})


test_that("nn.name = NULL auto-detects the lone default Neighbor object", {
  obj <- make_sample_fixture()  # carries a single 'weighted.nn'

  # NULL must resolve to the one present default rather than erroring
  res <- invisible(capture.output(
    out <- suppressWarnings(suppressMessages(
      GenerateSampleObject(object = obj, nn.name = NULL, group.by = "donor")
    ))
  ))
  expect_s4_class(out, "Seurat")
})


test_that("nn.name = NULL errors when no default Neighbor object is present", {
  obj <- make_sample_fixture()
  obj[["weighted.nn"]] <- NULL  # remove the only default

  expect_error(
    suppressMessages(GenerateSampleObject(object = obj, nn.name = NULL,
                                          group.by = "donor")),
    "No default NN object"
  )
})


test_that("nn.name = NULL errors when both default Neighbor objects are present", {
  obj <- make_sample_fixture()
  obj[["single.nn"]] <- obj[["weighted.nn"]]  # now both defaults exist

  expect_error(
    suppressMessages(GenerateSampleObject(object = obj, nn.name = NULL,
                                          group.by = "donor")),
    "Multiple default NN objects"
  )
})


test_that("k.nn cannot exceed the neighbours stored in the NN object", {
  obj <- make_sample_fixture(k = 10)

  expect_error(
    gen(obj, group.by = "donor", k.nn = 30),
    "exceeds the number of neighbors"
  )
  expect_s4_class(gen(obj, group.by = "donor", k.nn = 5), "Seurat")
  expect_s4_class(gen(obj, group.by = "donor", k.nn = 10), "Seurat")
})


test_that("normalization.method is validated rather than silently skipped", {
  obj <- make_sample_fixture()

  expect_error(
    gen(obj, group.by = "donor", normalization.method = "chisquared"),
    "should be one of"
  )
  expect_s4_class(
    gen(obj, group.by = "donor", normalization.method = "ChiSquared"), "Seurat"
  )
  expect_s4_class(
    gen(obj, group.by = "donor", normalization.method = "LogNormalize"), "Seurat"
  )
})


test_that("landmark cells are never counted towards their own sample", {
  n.cells <- 60
  n.landmark <- 8
  k.nn <- 5
  obj <- make_sample_fixture(n.cells = n.cells, n.landmark = n.landmark)

  mat <- gen(obj, group.by = "donor", k.nn = k.nn, return.seurat = FALSE)

  expect_equal(dim(mat), c(n.landmark, 6L))
  # every non-landmark cell contributes exactly k.nn counts; landmarks contribute none
  expect_equal(sum(mat), k.nn * (n.cells - n.landmark))

  # the toggle that used to gate this behaviour is gone
  expect_false("remove.sketch.cell.from.col" %in% names(formals(GenerateSampleObject)))
})


test_that("a sample with a single non-landmark cell does not drop to a vector", {
  obj <- make_sample_fixture()
  # cell9 is the only non-landmark cell assigned to sample "SOLO"
  obj$donor2 <- obj$donor
  obj$donor2[9] <- "SOLO"

  mat <- gen(obj, group.by = "donor2", add.meta.data = FALSE, return.seurat = FALSE)

  expect_equal(nrow(mat), 8L)
  expect_true("SOLO" %in% colnames(mat))
  expect_equal(sum(mat[, "SOLO"]), 5)
})


test_that("rename.group.by accepts several columns and uses the first", {
  obj <- make_sample_fixture()

  one <- gen(obj, group.by = "donor", rename.group.by = "sex", return.seurat = FALSE)
  two <- gen(obj, group.by = "donor", rename.group.by = c("sex", "region"),
             return.seurat = FALSE)

  # the second column must be ignored, not pasted into the names
  expect_identical(rownames(one), rownames(two))
  expect_true(all(grepl(pattern = "^(M|F)_LM[0-9]+$", x = rownames(two))))
})


test_that("rename.group.by falls back to landmark IDs when no column resolves", {
  obj <- make_sample_fixture()

  # FetchData() errors when every requested var is missing, so the only way to
  # empty the intersect() is a column coming back under a different name than
  # was asked for (FetchData rewrites vars that resolve through a key).
  local_mocked_bindings(
    FetchData = function(object, vars, ...) {
      df <- SeuratObject::FetchData(object = object, vars = vars, ...)
      if (identical(vars, "sex")) colnames(df) <- "renamed_by_fetchdata"
      df
    }
  )

  expect_warning(
    mat <- GenerateSampleObject(object = obj, nn.name = "weighted.nn",
                                group.by = "donor", rename.group.by = "sex",
                                return.seurat = FALSE),
    "Cannot find the correct meta-data columns"
  )
  # names must fall back to the landmark cell IDs, not become character(0)
  expect_equal(nrow(mat), 8L)
  expect_true(all(grepl(pattern = "^cell[0-9]+_LM[0-9]+$", x = rownames(mat))))
})


test_that("group.by must be a single meta-data column", {
  obj <- make_sample_fixture()

  expect_error(
    gen(obj, group.by = c("batch", "replicate")),
    "must be a single meta-data column"
  )

  # the documented workaround: pre-merge the columns
  obj$batch_rep <- paste(obj$batch, obj$replicate, sep = "_")
  res <- gen(obj, group.by = "batch_rep")

  expect_equal(ncol(res), 4L)
  # sample-level metadata must actually attach, not come back all-NA
  expect_false(any(is.na(res$batch_rep)))
  expect_setequal(res$batch_rep, c("b1_r1", "b1_r2", "b2_r1", "b2_r2"))
})


test_that("a group.by column with a single level is rejected", {
  obj <- make_sample_fixture()

  expect_error(gen(obj, group.by = "only_one"), "has only one value")
})
