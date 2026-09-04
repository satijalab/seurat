set.seed(42)

# An SCT model records the cells it was fit on, and subsetting does not prune
# that record, so a model can name cells the assay no longer has. Those have no
# new name, and assigning NA left the frame with missing row names:
#   Error in `.rowNamesDF<-`(x, value = value): missing values in 'row.names'
build_sct <- function() {
  data("pbmc_small", package = "SeuratObject", envir = environment())
  suppressWarnings(SCTransform(
    pbmc_small,
    variable.features.n = 20,
    vst.flavor = "v1",
    verbose = FALSE
  ))
}

add_ghost_cell <- function(object, model = 1L) {
  attributes <- slot(slot(object[["SCT"]], name = "SCTModel.list")[[model]], name = "cell.attributes")
  ghost <- attributes[1, , drop = FALSE]
  rownames(ghost) <- "ghost_cell"
  slot(slot(object[["SCT"]], name = "SCTModel.list")[[model]], name = "cell.attributes") <-
    rbind(attributes, ghost)
  object
}

cell_attributes <- function(object, model = 1L) {
  slot(slot(object[["SCT"]], name = "SCTModel.list")[[model]], name = "cell.attributes")
}

test_that("renaming works when a model names a cell the assay lost", {
  object <- add_ghost_cell(build_sct())
  renamed <- expect_no_error(RenameCells(object, new.names = paste0(Cells(object), "_query")))
  expect_identical(colnames(renamed), paste0(colnames(object), "_query"))

  attributes <- cell_attributes(renamed)
  expect_false(anyNA(rownames(attributes)))
  # the cells the assay has are renamed, the one it does not is left alone
  expect_true(all(paste0(colnames(object), "_query")[1:5] %in% rownames(attributes)))
  expect_true("ghost_cell" %in% rownames(attributes))
})

test_that("ordinary renaming is unchanged", {
  object <- build_sct()
  renamed <- RenameCells(object, new.names = paste0(Cells(object), "_query"))
  expect_identical(colnames(renamed), paste0(colnames(object), "_query"))
  expect_setequal(rownames(cell_attributes(renamed)), colnames(renamed))
})
