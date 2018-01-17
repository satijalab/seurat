#' @include seurat.R
NULL

globalVariables(names = 'cell.name', package = 'Seurat', add = TRUE)
#' Merge Seurat Objects
#'
#' Merge two Seurat objects
#'
#' @param object1 First Seurat object to merge
#' @param object2 Second Seurat object to merge
#' @param project Project name (string)
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
#' @param add.cell.id1 String to be appended to the names of all cells in object1
#' @param add.cell.id2 String to be appended to the names of all cells in object2
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
MergeSeurat <- function(
  object1,
  object2,
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
  add.cell.id1 = NULL,
  add.cell.id2 = NULL
) {
  if (length(x = object1@raw.data) < 2) {
    stop("First object provided has an empty raw.data slot. Adding/Merging performed on raw count data.")
  }
  if (length(x = object2@raw.data) < 2) {
    stop("Second object provided has an empty raw.data slot. Adding/Merging performed on raw count data.")
  }
  if (!missing(x = add.cell.id1)) {
    object1@cell.names <- paste(add.cell.id1,object1@cell.names, sep = "_")
    colnames(x = object1@raw.data) <- paste(
      add.cell.id1,
      colnames(x = object1@raw.data),
      sep = "_"
    )
    rownames(x = object1@meta.data) <- paste(
      add.cell.id1,
      rownames(x = object1@meta.data),
      sep = "_"
    )
  }
  if (!missing(x = add.cell.id2)) {
  object2@cell.names <- paste(add.cell.id2,object2@cell.names, sep = "_")
    colnames(x = object2@raw.data) <- paste(
      add.cell.id2,
      colnames(x = object2@raw.data),
      sep = "_"
    )
    rownames(x = object2@meta.data) <- paste(
      add.cell.id2,
      rownames(x = object2@meta.data),
      sep = "_"
    )
  }
  if (any(object1@cell.names %in% object2@cell.names)) {
    stop("Duplicate cell names, please provide 'add.cell.id1' and/or 'add.cell.id2' for unique names")
  }
  merged.raw.data <- RowMergeSparseMatrices(
    mat1 = object1@raw.data[,object1@cell.names],
    mat2 = object2@raw.data[,object2@cell.names]
  )
  object1@meta.data <- object1@meta.data[object1@cell.names, ]
  object2@meta.data <- object2@meta.data[object2@cell.names, ]
  project <- SetIfNull(x = project, default = object1@project.name)
  object1@meta.data$cell.name <- rownames(x = object1@meta.data)
  object2@meta.data$cell.name <- rownames(x = object2@meta.data)
  merged.meta.data <- suppressMessages(
    suppressWarnings(
      full_join(x = object1@meta.data, y = object2@meta.data)
    )
  )
  merged.object <- CreateSeuratObject(
    raw.data = merged.raw.data,
    project = project,
    min.cells = min.cells,
    min.genes = min.genes,
    is.expr = is.expr,
    normalization.method = NULL,
    scale.factor = scale.factor,
    do.scale = FALSE,
    do.center = FALSE,
    names.field = names.field,
    names.delim = names.delim
  )

  if (do.normalize) {
    normalization.method.use = GetCalcParam(
      object = object1,
      calculation = "NormalizeData",
      parameter = "normalization.method"
    )
    scale.factor.use = GetCalcParam(
      object = object1,
      calculation = "NormalizeData",
      parameter = "scale.factor"
    )

    if (is.null(normalization.method.use)) {
      normalization.method.use="LogNormalize"
      scale.factor.use=10000
    }
    merged.object <- NormalizeData(
      object = merged.object,
      assay.type = "RNA",
      normalization.method=normalization.method.use,
      scale.factor=scale.factor.use

    )
  }

  if (do.scale | do.center) {
    merged.object <- ScaleData(
      object = merged.object,
      do.scale = do.scale,
      do.center = do.center
    )
  }

  merged.meta.data %>% filter(
    cell.name %in% merged.object@cell.names
  ) -> merged.meta.data
  rownames(x= merged.meta.data) <- merged.object@cell.names
  merged.meta.data$cell.name <- NULL
  merged.object@meta.data <- merged.meta.data
  return(merged.object)
}

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
#' @param add.cell.id String to be appended to the names of all cells in new.data. E.g. if add.cell.id = "rep1",
#' "cell1" becomes "cell1.rep1"
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

#' Return a subset of the Seurat object
#'
#' Creates a Seurat object containing only a subset of the cells in the
#' original object. Takes either a list of cells to use as a subset, or a
#' parameter (for example, a gene), to subset on.
#'
#' @param object Seurat object
#' @param cells.use A vector of cell names to use as a subset. If NULL
#' (default), then this list will be computed based on the next three
#' arguments. Otherwise, will return an object consissting only of these cells
#' @param subset.name Parameter to subset on. Eg, the name of a gene, PC1, a
#' column name in object@@meta.data, etc. Any argument that can be retreived
#' using FetchData
#' @param ident.use Create a cell subset based on the provided identity classes
#' @param ident.remove Subtract out cells from these identity classes (used for filtration)
#' @param accept.low Low cutoff for the parameter (default is -Inf)
#' @param accept.high High cutoff for the parameter (default is Inf)
#' @param do.center Recenter the new object@@scale.data
#' @param do.scale Rescale the new object@@scale.data. FALSE by default
#' @param max.cells.per.ident Can be used to downsample the data to a certain max per cell ident. Default is inf.
#' @param random.seed Random seed for downsampling
#' @param do.clean Only keep object@@raw.data and object@@data. Cleans out most
#' other slots. Can be useful if you want to start a fresh analysis on just a
#' subset of the data. Also clears out stored clustering results in
#' object@@meta.data (any columns containing "res"). Will by default subset the
#' raw.data slot.
#' @param subset.raw Also subset object@@raw.data
#' @param \dots Additional arguments to be passed to FetchData (for example,
#' use.imputed=TRUE)
#'
#' @return Returns a Seurat object containing only the relevant subset of cells
#'
#' @importFrom stats complete.cases
#'
#' @export
#'
#' @examples
#' pbmc1 <- SubsetData(object = pbmc_small, cells.use = pbmc_small@cell.names[1:40])
#' pbmc1
#'
SubsetData <- function(
  object,
  cells.use = NULL,
  subset.name = NULL,
  ident.use = NULL,
  ident.remove = NULL,
  accept.low = -Inf,
  accept.high = Inf,
  do.center = FALSE,
  do.scale = FALSE,
  max.cells.per.ident = Inf,
  random.seed = 1,
  do.clean = FALSE,
  subset.raw,
  ...
) {
  data.use <- NULL
  cells.use <- SetIfNull(x = cells.use, default = object@cell.names)
  if (!is.null(x = ident.use)) {
    ident.use <- setdiff(ident.use, ident.remove)
    cells.use <- WhichCells(object, ident.use)
  }
  if ((is.null(x = ident.use)) && ! is.null(x = ident.remove)) {
    ident.use <- setdiff(unique(object@ident), ident.remove)
    cells.use <- WhichCells(object, ident.use)
  }
  if (! is.null(x = subset.name)) {
    data.use <- FetchData(object, subset.name, ...)
    if (length(x = data.use) == 0) {
      return(object)
    }
    subset.data <- data.use[, subset.name]
    pass.inds <- which(x = (subset.data > accept.low) & (subset.data < accept.high))
    cells.use <- rownames(data.use)[pass.inds]
  }
  cells.use <- intersect(x = cells.use, y = object@cell.names)
  cells.use <-  WhichCells(
    object = object,
    cells.use = cells.use,
    max.cells.per.ident = max.cells.per.ident,
    random.seed = random.seed
  )
  object@cell.names <- cells.use
  object@data <- object@data[, cells.use]
  if(! is.null(x = object@scale.data)) {
    if (length(x = colnames(x = object@scale.data) > 0)) {
      object@scale.data[, cells.use]
      object@scale.data <- object@scale.data[
        complete.cases(object@scale.data), # Row
        cells.use # Columns
        ]
    }
  }
  if (do.scale) {
    object <- ScaleData(
      object = object,
      do.scale = do.scale,
      do.center = do.center
    )
    object@scale.data <- object@scale.data[
      complete.cases(object@scale.data), # Row
      cells.use # Column
      ]
  }
  object@ident <- drop.levels(x = object@ident[cells.use])
  if (length(x = object@dr) > 0) {
    for (i in 1:length(object@dr)) {
      if(length(object@dr[[i]]@cell.embeddings) > 0){
        object@dr[[i]]@cell.embeddings <- object@dr[[i]]@cell.embeddings[cells.use, ,drop = FALSE]
      }
    }
  }
  # handle multimodal casess
  if (! .hasSlot(object = object, name = "assay")) {
    object@assay <- list()
  }
  if (length(object@assay) > 0) {
    for(i in 1:length(object@assay)) {
      if ((! is.null(x = object@assay[[i]]@raw.data)) && (ncol(x = object@assay[[i]]@raw.data) > 1)) {
        object@assay[[i]]@raw.data <- object@assay[[i]]@raw.data[, cells.use]
      }
      if ((! is.null(x = object@assay[[i]]@data)) && (ncol(x = object@assay[[i]]@data) > 1)) {
        object@assay[[i]]@data <- object@assay[[i]]@data[, cells.use]
      }
      if ((! is.null(x = object@assay[[i]]@scale.data)) && (ncol(x = object@assay[[i]]@scale.data) > 1)) {
        object@assay[[i]]@scale.data <- object@assay[[i]]@scale.data[, cells.use]
      }
    }
  }
  #object@tsne.rot=object@tsne.rot[cells.use, ]
  # object@gene.scores <- data.frame(object@gene.scores[cells.use,])
  # colnames(x = object@gene.scores)[1] <- "nGene"
  # rownames(x = object@gene.scores) <- colnames(x = object@data)
  object@meta.data <- data.frame(object@meta.data[cells.use,])
  #object@mix.probs=data.frame(object@mix.probs[cells.use,]); colnames(object@mix.probs)[1]="nGene"; rownames(object@mix.probs)=colnames(object@data)
  if (do.clean){
    calcs.to.keep <- c("CreateSeuratObject", "NormalizeData", "ScaleData")
    object@calc.params <- object@calc.params[calcs.to.keep]
    object@var.genes <- vector()
    object@hvg.info <- data.frame()
    object@cluster.tree <- list()
    object@snn <- as(matrix(), 'dgCMatrix')
    object@scale.data <- matrix()
    object@misc <- NULL
    object@kmeans <- NULL
    object@dr <- list()
    object@meta.data
    if(missing(subset.raw)) {
      subset.raw <- TRUE
    }
    object@meta.data[, sapply(colnames(object@meta.data), function(x){grepl("res", x)})] <- NULL
  }
  if(!missing(subset.raw)){
    if(subset.raw){
      object@raw.data <- object@raw.data[, cells.use]
    }
  }
  return(object)
}

#' Reorder identity classes
#'
#' Re-assigns the identity classes according to the average expression of a particular feature (i.e, gene expression, or PC score)
#' Very useful after clustering, to re-order cells, for example, based on PC scores
#'
#' @param object Seurat object
#' @param feature Feature to reorder on. Default is PC1
#' @param rev Reverse ordering (default is FALSE)
#' @param aggregate.fxn Function to evaluate each identity class based on (default is mean)
#' @param reorder.numeric Rename all identity classes to be increasing numbers starting from 1 (default is FALSE)
#' @param \dots additional arguemnts (i.e. use.imputed=TRUE)
#'
#' @return A seurat object where the identity have been re-oredered based on the average.
#'
#' @export
#'
#' @examples
#' head(x = pbmc_small@ident)
#' pbmc_small <- ReorderIdent(object = pbmc_small)
#' head(x = pbmc_small@ident)
#'
ReorderIdent <- function(
  object,
  feature = "PC1",
  rev = FALSE,
  aggregate.fxn = mean,
  reorder.numeric = FALSE,
  ...
) {
  ident.use <- object@ident
  data.use <- FetchData(object = object, vars.all = feature, ...)[, 1]
  revFxn <- Same
  if (rev) {
    revFxn <- function(x) {
      return(max(x) + 1 - x)
    }
  }
  names.sort <- names(
    x = revFxn(
      sort(
        x = tapply(
          X = data.use,
          INDEX = (ident.use),
          FUN = aggregate.fxn
        )
      )
    )
  )
  ident.new <- factor(x = ident.use, levels = names.sort, ordered = TRUE)
  if (reorder.numeric) {
    ident.new <- factor(
      x = revFxn(
        rank(
          tapply(
            X = data.use,
            INDEX = as.numeric(x = ident.new),
            FUN = mean
          )
        )
      )[as.numeric(ident.new)],
      levels = 1:length(x = levels(x = ident.new)),
      ordered = TRUE
    )
  }
  names(x = ident.new) <- names(x = ident.use)
  object@ident <- ident.new
  return(object)
}

#' Access cellular data
#'
#' Retreives data (gene expression, PCA scores, etc, metrics, etc.) for a set
#' of cells in a Seurat object
#'
#' @param object Seurat object
#' @param vars.all List of all variables to fetch
#' @param cells.use Cells to collect data for (default is all cells)
#' @param use.imputed For gene expression, use imputed values. Default is FALSE
#' @param use.scaled For gene expression, use scaled values. Default is FALSE
#' @param use.raw For gene expression, use raw values. Default is FALSE
#'
#' @return A data frame with cells as rows and cellular data as columns
#'
#' @export
#'
#' @examples
#' pc1 <- FetchData(object = pbmc_small, vars.all = 'PC1')
#' head(x = pc1)
#'
FetchData <- function(
  object,
  vars.all = NULL,
  cells.use = NULL,
  use.imputed = FALSE,
  use.scaled = FALSE,
  use.raw = FALSE
) {
  cells.use <- SetIfNull(x = cells.use, default = object@cell.names)
  data.return <- data.frame(row.names = cells.use)
  data.expression <- as.matrix(x = data.frame(row.names = cells.use))
  if (length(which(c(use.imputed, use.scaled, use.raw))) > 1) {
    stop("Can only set one of the following to TRUE: use.imputed, use.scaled, use.raw")
  }
  slot.use <- "data"
  # if any vars passed are genes, subset expression data
  gene.check <- vars.all %in% rownames(object@data)
  if (use.scaled) {
    slot.use <- "scale.data"
    gene.check <- vars.all %in% rownames(object@scale.data)
  }
  if (use.raw) {
    slot.use <- "raw.data"
  }
  if (any(gene.check)) {
    if (use.imputed) {
      gene.check <- vars.all %in% rownames(object@imputed)
      if (length(object@imputed) == 0) {
        stop ("Imputed expression values not calculated yet.")
      }
      data.expression <- t(object@imputed[vars.all[gene.check], cells.use, drop = FALSE])
    } else {
      data.expression <- GetAssayData(object, assay.type = "RNA", slot = slot.use)
      data.expression <- t(data.expression[vars.all[gene.check], cells.use, drop = FALSE])
    }
    if (all(gene.check)) {
      return(as.matrix(x = data.expression))
    }
  }
  # now check for multimodal data
  if (length(x = object@assay) > 0) {
    data.types <- names(x = object@assay)
    for (data.type in data.types) {
      all_data <- (GetAssayData(
        object = object,
        assay.type = data.type,
        slot = slot.use
      ))
      genes.include <- intersect(x = vars.all, y = rownames(x = all_data))
      data.expression <- cbind(
        data.expression,
        t(x = all_data[genes.include, , drop = FALSE])
      )
    }
  }
  var.options <- c("meta.data", "mix.probs", "gene.scores")
  if (length(x = names(x = object@dr)) > 0) {
    dr.options <- names(x = object@dr)
    dr.names <- paste0("dr$", names(x = object@dr), "@key")
    dr.names <- sapply(
      X = dr.names,
      FUN = function(x) {
        return(eval(expr = parse(text = paste0("object@", x))))
      }
    )
    names(x = dr.names) <- dr.options
    var.options <- c(var.options, dr.names)
  }
  object@meta.data[,"ident"] <- object@ident[rownames(x = object@meta.data)]
  for (my.var in vars.all) {
    data.use=data.frame()
    if (my.var %in% colnames(data.expression)) {
      data.use <- data.expression
    } else {
      for(i in var.options) {
        if (all(unlist(x = strsplit(x = my.var, split = "[0-9]+")) == i)) {
          eval(
            expr = parse(
              text = paste0(
                "data.use <- object@dr$",
                names(x = var.options[which(i == var.options)]),
                "@cell.embeddings"
              )
            )
          )
          colnames(x = data.use) <- paste0(i, 1:ncol(x = data.use))
          break
        }
      }
    }
    if (my.var %in% colnames(object@meta.data)) {
      data.use <- object@meta.data[, my.var, drop = FALSE]
    }
    if (ncol(x = data.use) == 0) {
      stop(paste("Error:", my.var, "not found"))
    }
    cells.use <- intersect(x = cells.use, y = rownames(x = data.use))
    if (! my.var %in% colnames(x = data.use)) {
      stop(paste("Error:", my.var, "not found"))
    }
    data.add <- data.use[cells.use, my.var]
    if (is.null(x = data.add)) {
      stop(paste("Error:", my.var, "not found"))
    }
    data.return <- cbind(data.return, data.add)
  }
  colnames(x = data.return) <- vars.all
  rownames(x = data.return) <- cells.use
  return(data.return)
}

#' FastWhichCells
#' Identify cells matching certain criteria (limited to character values)
#' @param object Seurat object
#' @param group.by Group cells in different ways (for example, orig.ident). Should be a column name in object@meta.data
#' @param subset.value  Return cells matching this value
#' @param invert invert cells to return.FALSE by default
#'
#' @export
#'
#' @examples
#' FastWhichCells(object = pbmc_small, group.by = 'res.1', subset.value = 1)
#'
FastWhichCells <- function(object, group.by, subset.value, invert = FALSE) {
  object <- SetAllIdent(object = object, id = group.by)
  cells.return <- WhichCells(object = object, ident = subset.value)
  if (invert) {
    cells.return <- setdiff(x = object@cell.names, y = cells.return)
  }
  return(cells.return)
}

#' Identify cells matching certain criteria
#'
#' Returns a list of cells that match a particular set of criteria such as
#' identity class, high/low values for particular PCs, ect..
#'
#' @param object Seurat object
#' @param ident Identity classes to subset. Default is all identities.
#' @param ident.remove Indentity classes to remove. Default is NULL.
#' @param cells.use Subset of cell names
#' @param subset.name Parameter to subset on. Eg, the name of a gene, PC1, a
#' column name in object@@meta.data, etc. Any argument that can be retreived
#' using FetchData
#' @param accept.low Low cutoff for the parameter (default is -Inf)
#' @param accept.high High cutoff for the parameter (default is Inf)
#' @param accept.value Returns all cells with the subset name equal to this value
#' @param max.cells.per.ident Can be used to downsample the data to a certain max per cell ident. Default is inf.
#' @param random.seed Random seed for downsampling
#'
#' @return A vector of cell names
#'
#' @export
#'
#' @examples
#' WhichCells(object = pbmc_small, ident = 2)
#'
WhichCells <- function(
  object,
  ident = NULL,
  ident.remove = NULL,
  cells.use = NULL,
  subset.name = NULL,
  accept.low = -Inf,
  accept.high = Inf,
  accept.value = NULL,
  max.cells.per.ident = Inf,
  random.seed = 1
) {
  set.seed(seed = random.seed)
  cells.use <- SetIfNull(x = cells.use, default = object@cell.names)
  ident <- SetIfNull(x = ident, default = unique(x = object@ident))
  ident <- setdiff(x = ident, y = ident.remove)
  if (! all(ident %in% unique(x = object@ident))) {
    bad.idents <- ident[! (ident %in% unique(x = object@ident))]
    stop(paste("Identity :", bad.idents, "not found.   "))
  }
  cells.to.use <- character()
  for (id in ident) {
    cells.in.ident <- object@ident[cells.use]
    cells.in.ident <- names(x = cells.in.ident[cells.in.ident == id])
    cells.in.ident <- cells.in.ident[! is.na(x = cells.in.ident)]
    if (length(x = cells.in.ident) > max.cells.per.ident) {
      cells.in.ident <- sample(x = cells.in.ident, size = max.cells.per.ident)
    }
    cells.to.use <- c(cells.to.use, cells.in.ident)
  }
  cells.use <- cells.to.use
  if (! is.null(x = subset.name)){
    subset.name <- as.character(subset.name)
    data.use <- FetchData(
      object = object,
      vars.all = subset.name,
      cells.use = cells.use
    )
    if (length(x = data.use) == 0) {
      stop(paste("Error : ", id, " not found"))
    }
    subset.data <- data.use[, subset.name, drop = F]
    if(! is.null(x = accept.value)) {
      pass.inds <- which(x = subset.data == accept.value)
    } else {
      pass.inds <- which(x = (subset.data > accept.low) & (subset.data < accept.high))
    }
    cells.use <- rownames(x = data.use)[pass.inds]
  }
  return(cells.use)
}

#' Switch identity class definition to another variable
#'
#' @param object Seurat object
#' @param id Variable to switch identity class to (for example, 'DBclust.ident', the output
#' of density clustering) Default is orig.ident - the original annotation pulled from the cell name.
#'
#' @return A Seurat object where object@@ident has been appropriately modified
#'
#' @export
#'
#' @examples
#' head(x = pbmc_small@ident)
#' pbmc_small <- SetAllIdent(object = pbmc_small, id = 'orig.ident')
#' head(x = pbmc_small@ident)
#'
SetAllIdent <- function(object, id = NULL) {
  id <- SetIfNull(x = id, default = "orig.ident")
  if (id %in% colnames(x = object@meta.data)) {
    cells.use <- rownames(x = object@meta.data)
    ident.use <- object@meta.data[, id]
    object <- SetIdent(
      object = object,
      cells.use = cells.use,
      ident.use = ident.use
    )
  }
  return(object)
}

#' Rename one identity class to another
#'
#' Can also be used to join identity classes together (for example, to merge clusters).
#'
#' @param object Seurat object
#' @param old.ident.name The old identity class (to be renamed)
#' @param new.ident.name The new name to apply
#'
#' @return A Seurat object where object@@ident has been appropriately modified
#'
#' @export
#'
#' @examples
#' head(x = pbmc_small@ident)
#' pbmc_small <- RenameIdent(
#'   object = pbmc_small,
#'   old.ident.name = 0,
#'   new.ident.name = 'cluster_0'
#' )
#' head(x = pbmc_small@ident)
#'
RenameIdent <- function(object, old.ident.name = NULL, new.ident.name = NULL) {
  if (! old.ident.name %in% object@ident) {
    stop(paste("Error : ", old.ident.name, " is not a current identity class"))
  }
  new.levels <- old.levels <- levels(x = object@ident)
  # new.levels <- old.levels
  if (new.ident.name %in% old.levels) {
    new.levels <- new.levels[new.levels != old.ident.name]
  }
  if(! (new.ident.name %in% old.levels)) {
    new.levels[new.levels == old.ident.name] <- new.ident.name
  }
  ident.vector <- as.character(x = object@ident)
  names(x = ident.vector) <- names(object@ident)
  ident.vector[WhichCells(object = object, ident = old.ident.name)] <- new.ident.name
  object@ident <- factor(x = ident.vector, levels = new.levels)
  return(object)
}

#' Set identity class information
#'
#' Stashes the identity in data.info to be retrieved later. Useful if, for example, testing multiple clustering parameters
#'
#' @param object Seurat object
#' @param save.name Store current object@@ident under this column name in object@@meta.data. Can be easily retrived with SetAllIdent
#'
#' @return A Seurat object where object@@ident has been appropriately modified
#'
#' @export
#'
#' @examples
#' head(x = pbmc_small@meta.data)
#' pbmc_small <- StashIdent(object = pbmc_small, save.name = 'cluster.ident')
#' head(x = pbmc_small@meta.data)
#'
StashIdent <- function(object, save.name = "oldIdent") {
  object@meta.data[, save.name] <- as.character(x = object@ident)
  return(object)
}

#' Set identity class information
#'
#' Sets the identity class value for a subset (or all) cells
#'
#' @param object Seurat object
#' @param cells.use Vector of cells to set identity class info for (default is
#' all cells)
#' @param ident.use Vector of identity class values to assign (character
#' vector)
#'
#' @return A Seurat object where object@@ident has been appropriately modified
#'
#' @importFrom gdata drop.levels
#'
#' @export
#'
#' @examples
#' cluster2 <- WhichCells(object = pbmc_small, ident = 2)
#' pbmc_small@ident[cluster2]
#' pbmc_small <- SetIdent(
#'   object = pbmc_small,
#'   cells.use = cluster2,
#'   ident.use = 'cluster_2'
#' )
#' pbmc_small@ident[cluster2]
#'
SetIdent <- function(object, cells.use = NULL, ident.use = NULL) {
  cells.use <- SetIfNull(x = cells.use, default = object@cell.names)
  if (length(x = setdiff(x = cells.use, y = object@cell.names) > 0)) {
    stop(paste(
      "ERROR : Cannot find cells ",
      setdiff(x = cells.use, y = object@cell.names)
    ))
  }
  ident.new <- setdiff(x = ident.use, y = levels(x = object@ident))
  object@ident <- factor(
    x = object@ident,
    levels = unique(
      x = c(
        as.character(x = object@ident),
        as.character(x = ident.new)
      )
    )
  )
  object@ident[cells.use] <- ident.use
  object@ident <- drop.levels(x = object@ident)
  return(object)
}

#' Add Metadata
#'
#' Adds additional data for single cells to the Seurat object. Can be any piece
#' of information associated with a cell (examples include read depth,
#' alignment rate, experimental batch, or subpopulation identity). The
#' advantage of adding it to the Seurat object is so that it can be
#' analyzed/visualized using FetchData, VlnPlot, GenePlot, SubsetData, etc.
#'
#' @param object Seurat object
#' @param metadata Data frame where the row names are cell names (note : these
#' must correspond exactly to the items in object@@cell.names), and the columns
#' are additional metadata items.
#' @param col.name Name for metadata if passing in single vector of information
#'
#' @return Seurat object where the additional metadata has been added as
#' columns in object@@meta.data
#'
#' @export
#'
#' @examples
#' cluster_letters <- LETTERS[pbmc_small@ident]
#' pbmc_small <- AddMetaData(
#'   object = pbmc_small,
#'   metadata = cluster_letters,
#'   col.name = 'letter.idents'
#' )
#' head(x = pbmc_small@meta.data)
#'
AddMetaData <- function(object, metadata, col.name = NULL) {
  if (typeof(x = metadata) != "list") {
    metadata <- as.data.frame(x = metadata)
    if (is.null(x = col.name)) {
      stop("Please provide a name for provided metadata")
    }
    colnames(x = metadata) <- col.name
  }
  cols.add <- colnames(x = metadata)
  meta.add <- metadata[rownames(x = object@meta.data), cols.add]
  if (all(is.null(x = meta.add))) {
    stop("Metadata provided doesn't match the cells in this object")
  }
  object@meta.data[, cols.add] <- meta.add
  return(object)
}
