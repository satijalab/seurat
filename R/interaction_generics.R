#' Merge Seurat Objects
#'
#' Merge Seurat objects
#'
#' @param objects Vector of objects to merge
#' @param project New project name (string)
#' @param min.cells Include genes with detected expression in at least this
#' many cells
#' @param min.genes Include cells where at least this many genes are detected
#' @param is.expr Expression threshold for 'detected' gene
#' @param do.normalize Normalize the data after merging. Default is TRUE.
#' If set, will perform the same normalization strategy as stored for the first object
#' @param scale.factor If normalizing on the cell level, this sets the scale factor.
#' @param do.scale In object@@scale.data, perform row-scaling (gene-based
#' z-score). FALSE by default, so run ScaleData after merging.
#' @param do.center In object@@scale.data, perform row-centering (gene-based
#' centering). FALSE by default
#' @param names.field For the initial identity class for each cell, choose this
#' field from the cell's column name
#' @param names.delim For the initial identity class for each cell, choose this
#' delimiter from the cell's column name
#' @param add.cell.id1 String passed to \code{\link{RenameCells}} for object1
#' @param add.cell.id2 String passed to \code{\link{RenameCells}} for object1
#'
#' @return Merged Seurat object
#'
#' @import Matrix
#' @importFrom dplyr full_join filter
#'
#' @export
#'
#' @examples
#' # Split pbmc_small for this example
#' pbmc1 <- SubsetData(object = pbmc_small, cells.use = pbmc_small@cell.names[1:40])
#' pbmc1
#' pbmc2 <- SubsetData(object = pbmc_small, cells.use = pbmc_small@cell.names[41:80])
#' pbmc2
#' # Merge pbmc1 and pbmc2 into one Seurat object
#' pbmc_merged <- MergeSeurat(object1 = pbmc1, object2 = pbmc2)
#' pbmc_merged
#'
