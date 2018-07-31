#' Add samples into existing Seurat object.
#'
#' @param object Seurat object
#' @param project Project name (string)
#' @param new.data Data matrix for samples to be added
#' @param min.cells Include genes with detected expression in at least this
#' many cells
#' @param min.genes Include cells where at least this many genes are detected
#' @param is.expr Expression threshold for 'detected' gene
#' @param do.normalize Normalize the data after merging. Default is TRUE.
#' If set, will perform the same normalization strategy as stored in the object
#' @param scale.factor scale factor in the log normalization
#' @param do.scale In object@@scale.data, perform row-scaling (gene-based z-score)
#' @param do.center In object@@scale.data, perform row-centering (gene-based
#' centering)
#' @param names.field For the initial identity class for each cell, choose this
#' field from the cell's column name
#' @param names.delim For the initial identity class for each cell, choose this
#' delimiter from the cell's column name
#' @param meta.data Additional metadata to add to the Seurat object. Should be
#' a data frame where the rows are cell names, and the columns are additional
#' metadata fields
#' @param add.cell.id String to be appended to the names of all cells in
#' new.data. E.g. if add.cell.id = "rep1", "cell1" becomes "cell1.rep1"
#'
#' @import Matrix
#' @importFrom dplyr full_join
#'
#' @export
#'
#' @examples
#' pbmc1 <- SubsetData(object = pbmc_small, cells.use = pbmc_small@cell.names[1:40])
#' pbmc1
#' pbmc2 <- SubsetData(object = pbmc_small, cells.use = pbmc_small@cell.names[41:80])
#' pbmc2_data <- pbmc2@data
#' dim(pbmc2_data)
#' pbmc_added <- AddSamples(object = pbmc1, new.data = pbmc2_data)
#' pbmc_added
#'
AddSamples <- function(
  object,
  new.data,
  project = NULL,
  min.cells = 0,
  min.genes = 0,
  is.expr = 0,
  do.normalize = TRUE,
  scale.factor = 1e4,
  do.scale = FALSE,
  do.center = FALSE,
  names.field = 1,
  names.delim = "_",
  meta.data = NULL,
  add.cell.id = NULL
) {
  if (length(x = object@raw.data) < 2) {
    stop("Object provided has an empty raw.data slot. Adding/Merging performed on raw count data.")
  }
  if (!missing(x = add.cell.id)) {
    colnames(x = new.data) <- paste(
      add.cell.id,
      colnames(x = new.data),
      sep = "_"
    )
  }
  if (any(colnames(x = new.data) %in% object@cell.names)) {
    stop("Duplicate cell names, please provide 'add.cell.id' for unique names")
  }
  combined.data <- RowMergeSparseMatrices(
    mat1 = object@raw.data[, object@cell.names],
    mat2 = new.data
  )
  if (is.null(x = meta.data)) {
    filler <- matrix(NA, nrow = ncol(new.data), ncol = ncol(object@meta.data))
    rownames(filler) <- colnames(new.data)
    colnames(filler) <- colnames(object@meta.data)
    filler <- as.data.frame(filler)
    combined.meta.data <- rbind(object@meta.data, filler)
  } else {
    combined.meta.data <- suppressMessages(
      suppressWarnings(
        full_join(x = object@meta.data, y = meta.data)
      )
    )
  }
  combined.meta.data$nGene <- NULL
  combined.meta.data$nUMI <- NULL
  if (!is.null(x = add.cell.id)) {
    combined.meta.data$orig.ident <- factor(
      x = combined.meta.data$orig.ident,
      levels = c(levels(x = combined.meta.data$orig.ident), add.cell.id)
    )
    combined.meta.data[colnames(new.data), ] <- add.cell.id
  }
  project <- SetIfNull(x = project, default = object@project.name)
  new.object <- CreateSeuratObject(
    raw.data = combined.data,
    project = project,
    min.cells = min.cells,
    min.genes = min.genes,
    is.expr = is.expr,
    scale.factor = scale.factor,
    do.scale = F,
    do.center = F,
    names.field = names.field,
    names.delim = names.delim
  )
  if (do.normalize) {
    normalization.method.use = GetCalcParam(
      object = object,
      calculation = "NormalizeData",
      parameter = "normalization.method"
    )
    scale.factor.use = GetCalcParam(
      object = object,
      calculation = "NormalizeData",
      parameter = "scale.factor"
    )
    if (is.null(x = normalization.method.use)) {
      normalization.method.use <- "LogNormalize"
      scale.factor.use <- 10000
    }
    new.object <- NormalizeData(
      object = new.object,
      assay.type = "RNA",
      normalization.method = normalization.method.use,
      scale.factor = scale.factor.use
    )
  }
  if (do.scale | do.center) {
    new.object <- ScaleData(
      object = new.object,
      do.scale = do.scale,
      do.center = do.center
    )
  }
  new.object@meta.data$orig.ident <- NULL
  new.object@meta.data <- cbind(new.object@meta.data, combined.meta.data)
  return(new.object)
}

# #' Reorder identity classes
# #'
# #' Re-assigns the identity classes according to the average expression of a
# #' particular feature (i.e, gene expression, or PC score)
# #' Very useful after clustering, to re-order cells, for example, based on PC
# #' scores
# #'
# #' @param object Seurat object
# #' @param feature Feature to reorder on. Default is PC1
# #' @param rev Reverse ordering (default is FALSE)
# #' @param aggregate.fxn Function to evaluate each identity class based on
# #' (default is mean)
# #' @param reorder.numeric Rename all identity classes to be increasing numbers
# #' starting from 1 (default is FALSE)
# #' @param \dots additional arguemnts (i.e. use.imputed=TRUE)
# #'
# #' @return A seurat object where the identity have been re-oredered based on the
# #' average.
# #'
# #' @export
# #'
# #' @examples
# #' head(x = pbmc_small@ident)
# #' pbmc_small <- ReorderIdent(object = pbmc_small)
# #' head(x = pbmc_small@ident)
# #'
# ReorderIdent <- function(
#   object,
#   feature = "PC1",
#   rev = FALSE,
#   aggregate.fxn = mean,
#   reorder.numeric = FALSE,
#   ...
# ) {
#   ident.use <- object@ident
#   data.use <- FetchData(object = object, vars.all = feature, ...)[, 1]
#   revFxn <- Same
#   if (rev) {
#     revFxn <- function(x) {
#       return(max(x) + 1 - x)
#     }
#   }
#   names.sort <- names(
#     x = revFxn(
#       sort(
#         x = tapply(
#           X = data.use,
#           INDEX = (ident.use),
#           FUN = aggregate.fxn
#         )
#       )
#     )
#   )
#   ident.new <- factor(x = ident.use, levels = names.sort, ordered = TRUE)
#   if (reorder.numeric) {
#     ident.new <- factor(
#       x = revFxn(
#         rank(
#           tapply(
#             X = data.use,
#             INDEX = as.numeric(x = ident.new),
#             FUN = mean
#           )
#         )
#       )[as.numeric(ident.new)],
#       levels = 1:length(x = levels(x = ident.new)),
#       ordered = TRUE
#     )
#   }
#   names(x = ident.new) <- names(x = ident.use)
#   object@ident <- ident.new
#   return(object)
# }

# #' FastWhichCells
# #' Identify cells matching certain criteria (limited to character values)
# #' @param object Seurat object
# #' @param group.by Group cells in different ways (for example, orig.ident).
# #' Should be a column name in object@meta.data
# #' @param subset.value  Return cells matching this value
# #' @param invert invert cells to return.FALSE by default
# #'
# #' @export
# #'
# #' @examples
# #' FastWhichCells(object = pbmc_small, group.by = 'res.1', subset.value = 1)
# #'
# FastWhichCells <- function(object, group.by, subset.value, invert = FALSE) {
#   object <- SetAllIdent(object = object, id = group.by)
#   cells.return <- WhichCells(object = object, ident = subset.value)
#   if (invert) {
#     cells.return <- setdiff(x = object@cell.names, y = cells.return)
#   }
#   return(cells.return)
# }

# #' Switch identity class definition to another variable
# #'
# #' @param object Seurat object
# #' @param id Variable to switch identity class to (for example, 'DBclust.ident',
# #' the output of density clustering) Default is orig.ident - the original
# #' annotation pulled from the cell name.
# #'
# #' @return A Seurat object where object@@ident has been appropriately modified
# #'
# #' @export
# #'
# #' @examples
# #' head(x = pbmc_small@ident)
# #' pbmc_small <- SetAllIdent(object = pbmc_small, id = 'orig.ident')
# #' head(x = pbmc_small@ident)
# #'
# SetAllIdent <- function(object, id = NULL) {
#   id <- SetIfNull(x = id, default = "orig.ident")
#   if (id %in% colnames(x = object@meta.data)) {
#     cells.use <- rownames(x = object@meta.data)
#     ident.use <- object@meta.data[, id]
#     object <- SetIdent(
#       object = object,
#       cells.use = cells.use,
#       ident.use = ident.use
#     )
#   }
#   return(object)
# }

# #' Rename one identity class to another
# #'
# #' Can also be used to join identity classes together (for example, to merge
# #' clusters).
# #'
# #' @param object Seurat object
# #' @param old.ident.name The old identity class (to be renamed)
# #' @param new.ident.name The new name to apply
# #'
# #' @return A Seurat object where object@@ident has been appropriately modified
# #'
# #' @export
# #'
# #' @examples
# #' head(x = pbmc_small@ident)
# #' pbmc_small <- RenameIdent(
# #'   object = pbmc_small,
# #'   old.ident.name = 0,
# #'   new.ident.name = 'cluster_0'
# #' )
# #' head(x = pbmc_small@ident)
# #'
# RenameIdent <- function(object, old.ident.name = NULL, new.ident.name = NULL) {
#   if (!old.ident.name %in% object@ident) {
#     stop(paste("Error : ", old.ident.name, " is not a current identity class"))
#   }
#   new.levels <- old.levels <- levels(x = object@ident)
#   # new.levels <- old.levels
#   if (new.ident.name %in% old.levels) {
#     new.levels <- new.levels[new.levels != old.ident.name]
#   }
#   if (!(new.ident.name %in% old.levels)) {
#     new.levels[new.levels == old.ident.name] <- new.ident.name
#   }
#   ident.vector <- as.character(x = object@ident)
#   names(x = ident.vector) <- names(object@ident)
#   ident.vector[WhichCells(object = object, ident = old.ident.name)] <- new.ident.name
#   object@ident <- factor(x = ident.vector, levels = new.levels)
#   return(object)
# }

# #' Transfer identity class information (or meta data) from one object to another
# #'
# #' Transfers identity class information (or meta data) from one object to
# #' another, assuming the same cell barcode names are in each. Can be very useful
# #' if you have multiple Seurat objects that share a subset of underlying data.
# #'
# #' @param object.from Seurat object to transfer information from
# #' @param object.to Seurat object to transfer information onto
# #' @param data.to.transfer What data should be transferred over? Default is the
# #' identity class ("ident"), but can also include any column in
# #' object.from@@meta.data
# #' @param keep.existing For cells in object.to that are not present in
# #' object.from, keep existing data? TRUE by default. If FALSE, set to NA.
# #' @param add.cell.id1 Prefix to add (followed by an underscore) to cells in
# #'  object.from. NULL by default, in which case no prefix is added.
# #'
# #' @return A Seurat object where object@@ident or object@@meta.data has been
# #' appropriately modified
# #'
# #' @export
# #'
# #' @examples
# #' # Duplicate the test object and assign random new idents to transfer
# #' pbmc_small@@ident
# #' pbmc_small2 <- SetIdent(object = pbmc_small, cells.use = pbmc_small@@cell.names,
# #'  ident.use = sample(pbmc_small@@ident))
# #' pbmc_small2@@ident
# #' pbmc_small <- TransferIdent(object.from = pbmc_small2, object.to = pbmc_small)
# #' pbmc_small@@ident
# #'
# TransferIdent <- function(object.from, object.to, data.to.transfer = "ident", keep.existing = TRUE, add.cell.id1 = NULL) {
#   old_data <- as.character(FetchData(object = object.from, vars.all = data.to.transfer)[, 1])
#   names(old_data) <- object.from@cell.names
#   if (data.to.transfer %in% c("ident", colnames(object.to@meta.data))) {
#     new_data <- FetchData(object = object.to, vars.all = data.to.transfer)
#     if (!keep.existing) {
#       new_data[, 1] <- "NA"
#     }
#     new_data <- as.character(new_data[, 1])
#   }
#   else {
#     new_data <- rep("NA", length(object.to@cell.names))
#   }
#   names(new_data) <- object.to@cell.names
#   if (!is.null(add.cell.id1)) {
#     names(old_data) <- paste(names(old_data), add.cell.id1, sep = "_")
#   }
#   new_data[names(old_data)] <- old_data
#   if (data.to.transfer == "ident") {
#     object.to <- SetIdent(object.to, cells.use = names(new_data), ident.use = new_data)
#   }
#   else {
#     object.to <- AddMetaData(object = object.to, metadata = new_data,col.name = data.to.transfer)
#   }
#   return(object.to)
# }

# #' Sets identity class information to be a combination of two object attributes
# #'
# #' Combined two attributes to define identity classes. Very useful if, for
# #' example, you have multiple cell types and multiple replicates, and you want
# #' to group cells based on combinations of both.
# #'
# #' @param object Seurat object
# #' @param attribute.1 First attribute for combination. Default is "ident"
# #' @param attribute.2 Second attribute for combination. Default is "orig.ident"
# #' @return A Seurat object where object@@ident has been appropriately modified
# #'
# #' @export
# #'
# #' @examples
# #' groups <- sample(c("group1", "group2", "group3"), size = 80, replace = TRUE)
# #' celltype <- sample(c("celltype1", "celltype2", "celltype3"), size = 80, replace = TRUE)
# #' new.metadata <- data.frame(groups = groups, celltype = celltype)
# #' rownames(new.metadata) <- pbmc_small@@cell.names
# #' pbmc_small <- AddMetaData(object = pbmc_small, metadata = new.metadata)
# #' pbmc_small <- CombineIdent(object = pbmc_small, attribute.1 = "celltype", attribute.2 = "groups")
# #' pbmc_small@@ident
# #'
# CombineIdent <- function(object, attribute.1 = "ident", attribute.2 = "orig.ident") {
#   old_data <- FetchData(object = object, vars.all = c(attribute.1, attribute.2))
#   new_ids <- sapply(X = 1:nrow(old_data), FUN = function(x){
#     paste(as.character(old_data[x, 1]), as.character(old_data[x, 2]), sep = "_")
#     })
#   object <- SetIdent(object = object,cells.use = object@cell.names, ident.use = new_ids)
# }
