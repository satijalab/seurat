# Tests for RenameCells.SCTAssay in objects.R
suppressWarnings(RNGversion(vstr = "3.5.3"))
set.seed(seed = 42)
is_not_cran_submission <- isTRUE(as.logical(Sys.getenv("NOT_CRAN")))

context("RenameCells.SCTAssay")

# An SCT model records the cells it was fit on, and that record is meant to stay
# a subset of the assay's cells -- `subset.SCTAssay` prunes it, and
# `CreateSCTAssayObject` intersects it. No current Seurat operation is known to
# break that invariant, so these tests write the stale cells into the model slot
# directly to reproduce what objects saved by older versions of Seurat carry.
# Without the prune, the stale cells look up as `NA` and `rownames<-` rejects
# them: one stale cell reports "missing values in 'row.names' are not allowed",
# two or more report "duplicate 'row.names' are not allowed", because
# `.rowNamesDF<-` tests `anyDuplicated` before `anyNA`.

ModelAttributes <- function(object, model = 1L, assay = "SCT") {
  slot(
    object = slot(object = object[[assay]], name = "SCTModel.list")[[model]],
    name = "cell.attributes"
  )
}

`ModelAttributes<-` <- function(object, model = 1L, assay = "SCT", value) {
  slot(
    object = slot(object = object[[assay]], name = "SCTModel.list")[[model]],
    name = "cell.attributes"
  ) <- value
  return(object)
}

# Append rows named `cells` to a model's record of the cells it was fit on
AddStaleCells <- function(object, cells, model = 1L) {
  attributes <- ModelAttributes(object = object, model = model)
  stale <- attributes[seq_along(along.with = cells), , drop = FALSE]
  rownames(x = stale) <- cells
  ModelAttributes(object = object, model = model) <- rbind(attributes, stale)
  return(object)
}

if (is_not_cran_submission) {
  sct.obj <- suppressWarnings(SCTransform(
    object = pbmc_small,
    variable.features.n = 20,
    vst.flavor = "v1",
    verbose = FALSE
  ))
  new.names <- paste0(Cells(x = sct.obj), "_query")

  test_that("renaming is unchanged when the model matches the assay", {
    renamed <- expect_no_warning(
      RenameCells(object = sct.obj, new.names = new.names)
    )
    expect_identical(object = Cells(x = renamed), expected = new.names)
    expect_setequal(
      object = rownames(x = ModelAttributes(object = renamed)),
      expected = new.names
    )
  })

  test_that("a single stale model cell is dropped rather than renamed to NA", {
    object <- AddStaleCells(object = sct.obj, cells = "stale_cell")
    expect_warning(
      object = renamed <- RenameCells(object = object, new.names = new.names),
      regexp = "Dropping 1 cell"
    )
    attributes <- ModelAttributes(object = renamed)
    expect_false(object = anyNA(x = rownames(x = attributes)))
    expect_setequal(object = rownames(x = attributes), expected = new.names)
  })

  test_that("several stale model cells are dropped", {
    object <- AddStaleCells(object = sct.obj, cells = c("stale_one", "stale_two"))
    expect_warning(
      object = renamed <- RenameCells(object = object, new.names = new.names),
      regexp = "Dropping 2 cell"
    )
    expect_setequal(
      object = rownames(x = ModelAttributes(object = renamed)),
      expected = new.names
    )
  })

  test_that("a stale cell already carrying a new name does not collide", {
    object <- AddStaleCells(object = sct.obj, cells = new.names[5])
    expect_warning(
      object = renamed <- RenameCells(object = object, new.names = new.names),
      regexp = "Dropping 1 cell"
    )
    attributes <- ModelAttributes(object = renamed)
    expect_false(object = anyDuplicated(x = rownames(x = attributes)) > 0)
    expect_setequal(object = rownames(x = attributes), expected = new.names)
  })

  test_that("the renamed object is left consistent for downstream indexing", {
    object <- AddStaleCells(object = sct.obj, cells = "stale_cell")
    renamed <- suppressWarnings(
      RenameCells(object = object, new.names = new.names)
    )
    # `PrepSCTFindMarkers` indexes the counts assay by the model's cells
    model.cells <- rownames(x = ModelAttributes(object = renamed))
    counts <- GetAssayData(object = renamed, assay = "RNA", layer = "counts")
    expect_length(object = setdiff(x = model.cells, y = Cells(x = renamed)), n = 0)
    expect_no_error(object = counts[, model.cells])
  })

  test_that("models of a merged assay are renamed independently", {
    first <- sct.obj[, 1:30]
    second <- RenameCells(
      object = sct.obj[, 31:60],
      new.names = paste0(Cells(x = sct.obj)[31:60], "_second")
    )
    merged <- merge(x = first, y = second)
    merged.names <- paste0(Cells(x = merged), "_query")
    renamed <- expect_no_warning(
      RenameCells(object = merged, new.names = merged.names)
    )
    expect_identical(object = Cells(x = renamed), expected = merged.names)
    for (model in levels(x = renamed[["SCT"]])) {
      cells <- rownames(x = SCTResults(
        object = renamed[["SCT"]],
        slot = "cell.attributes",
        model = model
      ))
      expect_length(object = setdiff(x = cells, y = merged.names), n = 0)
    }
  })
}
