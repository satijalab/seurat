#' @include seurat.R
#' @include interaction_generics.R
#' @importFrom methods setMethod
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



#' Return a subset of the Seurat object.
#'
#' Creates a Seurat object containing only a subset of the cells in the
#' original object. Forms a dataframe by fetching the variables in \code{vars.use}, then
#' subsets it using \code{base::subset} with \code{predicate} as the filter.
#' Returns the corresponding subset of the Seurat object.
#'
#' @param object Seurat object
#' @param vars.use Variables to fetch for use in base::subset. Character vector.
#' @param predicate String to be parsed into an R expression and evaluated as an input to base::subset.
#'
#' @export
#'
#' @examples
#' pbmc1 <- SubsetByPredicate(object = pbmc_small,
#'                       vars.use = c("nUMI", "res.1"),
#'                       predicate = "nUMI < 200 & res.1=='3'")
#' pbmc1
#'
SubsetByPredicate = function(
  object,
  vars.use,
  predicate
){
  if( typeof(vars.use) != "character"){
    stop("predicate should be a character vector. It will be parsed in `subset` as an R expression.")
  }
  if( typeof(predicate) != "character"){
    stop("vars.use should be a character vector. These variables will be passed to FetchData.")
  }
  df <- FetchData(object, vars.use) %>% as.data.frame
  cu <- df %>% subset(eval(parse(text=predicate))) %>% rownames
  object <- SubsetData(object, cells.use = cu)
  return( object )
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
#' @param ident.remove Subtract out cells from these identity classes (used for
#' filtration)
#' @param accept.low Low cutoff for the parameter (default is -Inf)
#' @param accept.high High cutoff for the parameter (default is Inf)
#' @param accept.value Returns cells with the subset name equal to this value
#' @param do.center Recenter the new object@@scale.data
#' @param do.scale Rescale the new object@@scale.data. FALSE by default
#' @param max.cells.per.ident Can be used to downsample the data to a certain
#'  max per cell ident. Default is INF.
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
  accept.value = NULL,
  do.center = FALSE,
  do.scale = FALSE,
  max.cells.per.ident = Inf,
  random.seed = 1,
  do.clean = FALSE,
  subset.raw,
  ...
) {
  data.use <- NULL
  cells.use <- WhichCells(object = object,
                          ident = ident.use,
                          ident.remove = ident.remove,
                          cells.use = cells.use,
                          subset.name = subset.name,
                          accept.low = accept.low,
                          accept.high = accept.high,
                          accept.value = accept.value,
                          max.cells.per.ident = max.cells.per.ident,
                          random.seed = random.seed,
                          ... = ...)
  object@cell.names <- cells.use
  object@data <- object@data[, cells.use]
  if(! is.null(x = object@scale.data)) {
    if (length(x = colnames(x = object@scale.data) > 0)) {
      object@scale.data[, cells.use]
      object@scale.data <- object@scale.data[, cells.use]
    }
  }
  if (do.scale) {
    object <- ScaleData(
      object = object,
      do.scale = do.scale,
      do.center = do.center
    )
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
#' Re-assigns the identity classes according to the average expression of a
#' particular feature (i.e, gene expression, or PC score)
#' Very useful after clustering, to re-order cells, for example, based on PC
#' scores
#'
#' @param object Seurat object
#' @param feature Feature to reorder on. Default is PC1
#' @param rev Reverse ordering (default is FALSE)
#' @param aggregate.fxn Function to evaluate each identity class based on
#' (default is mean)
#' @param reorder.numeric Rename all identity classes to be increasing numbers
#' starting from 1 (default is FALSE)
#' @param \dots additional arguemnts (i.e. use.imputed=TRUE)
#'
#' @return A seurat object where the identity have been re-oredered based on the
#' average.
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
#' @param group.by Group cells in different ways (for example, orig.ident).
#' Should be a column name in object@meta.data
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
#' @param max.cells.per.ident Can be used to downsample the data to a certain
#' max per cell ident. Default is INF.
#' @param random.seed Random seed for downsampling
#' @param \dots Additional arguments to be passed to FetchData (for example,
#' use.imputed=TRUE)
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
  random.seed = 1,
  ...
) {
  # input checking
  if(length(subset.name) > 1) {
    stop("subset.name must be a single parameter")
  }
  if(length(accept.low) > 1 | length(accept.high) > 1) {
    stop("Multiple values passed to accept.low or accept.high")
  }
  if(accept.low >= accept.high) {
    stop("accept.low greater than or equal to accept.high")
  }
  if (!is.na(x = random.seed)) {
    set.seed(seed = random.seed)
  }
  cells.use <- SetIfNull(x = cells.use, default = object@cell.names)
  ident <- SetIfNull(x = ident, default = unique(x = object@ident))
  bad.remove.idents <- ident.remove[! (ident.remove %in% unique(x = object@ident))]
  if(length(bad.remove.idents) > 0) {
    stop(paste("Identity :", bad.remove.idents, "not found.   "))
  }
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
      cells.use = cells.use,
      ... = ...
    )
    if (length(x = data.use) == 0) {
      stop(paste("Error : ", id, " not found"))
    }
    subset.data <- data.use[, subset.name, drop = F]
    if(! is.null(x = accept.value)) {
      if (! all(accept.value %in% unique(x = subset.data[, 1]))) {
        bad.vals <- accept.value[! (accept.value %in% unique(x = subset.data[, 1]))]
        stop(paste("Identity :", bad.vals, "not found.   "))
      }
      pass.inds <- which(apply(subset.data, MARGIN = 1, function(x) x %in% accept.value))
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
#' @param id Variable to switch identity class to (for example, 'DBclust.ident',
#' the output of density clustering) Default is orig.ident - the original
#' annotation pulled from the cell name.
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
#' Can also be used to join identity classes together (for example, to merge
#' clusters).
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
#' Stashes the identity in data.info to be retrieved later. Useful if, for
#' example, testing multiple clustering parameters
#'
#' @param object Seurat object
#' @param save.name Store current object@@ident under this column name in
#' object@@meta.data. Can be easily retrived with SetAllIdent
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


#' Transfer identity class information (or meta data) from one object to another
#'
#' Transfers identity class information (or meta data) from one object to
#' another, assuming the same cell barcode names are in each. Can be very useful
#' if you have multiple Seurat objects that share a subset of underlying data.
#'
#' @param object.from Seurat object to transfer information from
#' @param object.to Seurat object to transfer information onto
#' @param data.to.transfer What data should be transferred over? Default is the
#' identity class ("ident"), but can also include any column in
#' object.from@@meta.data
#' @param keep.existing For cells in object.to that are not present in
#' object.from, keep existing data? TRUE by default. If FALSE, set to NA.
#' @param add.cell.id1 Prefix to add (followed by an underscore) to cells in
#'  object.from. NULL by default, in which case no prefix is added.
#'
#' @return A Seurat object where object@@ident or object@@meta.data has been
#' appropriately modified
#'
#' @export
#'
#' @examples
#' # Duplicate the test object and assign random new idents to transfer
#' pbmc_small@@ident
#' pbmc_small2 <- SetIdent(object = pbmc_small, cells.use = pbmc_small@@cell.names,
#'  ident.use = sample(pbmc_small@@ident))
#' pbmc_small2@@ident
#' pbmc_small <- TransferIdent(object.from = pbmc_small2, object.to = pbmc_small)
#' pbmc_small@@ident
#'
TransferIdent <- function(object.from, object.to, data.to.transfer = "ident", keep.existing = TRUE, add.cell.id1 = NULL) {
  old_data <- as.character(FetchData(object = object.from, vars.all = data.to.transfer)[, 1])
  names(old_data) <- object.from@cell.names
  if (data.to.transfer %in% c("ident", colnames(object.to@meta.data))) {
    new_data <- FetchData(object = object.to, vars.all = data.to.transfer)
    if (!keep.existing) {
      new_data[, 1] <- "NA"
    }
    new_data <- as.character(new_data[, 1])
  }
  else {
    new_data <- rep("NA", length(object.to@cell.names))
  }
  names(new_data) <- object.to@cell.names
  if (!is.null(add.cell.id1)) {
    names(old_data) <- paste(names(old_data), add.cell.id1, sep = "_")
  }
  new_data[names(old_data)] <- old_data
  if (data.to.transfer == "ident") {
    object.to <- SetIdent(object.to, cells.use = names(new_data), ident.use = new_data)
  }
  else {
    object.to <- AddMetaData(object = object.to, metadata = new_data,col.name = data.to.transfer)
  }
  return(object.to)
}

#' Splits object into a list of subsetted objects.
#'
#' Splits object based on a single attribute into a list of subsetted objects,
#' one for each level of the attribute. For example, useful for taking an object
#' that contains cells from many patients, and subdividing it into
#' patient-specific objects.
#'
#' @param object Seurat object
#' @param attribute.1 Attribute for splitting. Default is "ident". Currently
#' only supported for class-level (i.e. non-quantitative) attributes.
#' @param \dots Additional parameters to pass to SubsetData
#' @return A named list of Seurat objects, each containing a subset of cells
#' from the original object.
#'
#' @export
#'
#' @examples
#' # Assign the test object a three level attribute
#' groups <- sample(c("group1", "group2", "group3"), size = 80, replace = TRUE)
#' names(groups) <- pbmc_small@@cell.names
#' pbmc_small <- AddMetaData(object = pbmc_small, metadata = groups, col.name = "group")
#' obj.list <- SplitObject(pbmc_small, attribute.1 = "group")
#'
SplitObject <- function(object,
                        attribute.1 = "ident",
                        ...) {
  old_data <- FetchData(object = object, vars.all = attribute.1)[, 1]
  old_levels <- unique(as.character(old_data))
  to_return <- list()
  for (i in old_levels) {
    if (attribute.1=="ident") {
      to_return[[i]] <- SubsetData(object = object, ident.use = i, ...)
    }
    else {
      to_return[[i]] <- SubsetData(object = object,
                                   subset.name = attribute.1,
                                   accept.value = i, ...)
    }
  }
  return(to_return)
}

#' Sets identity class information to be a combination of two object attributes
#'
#' Combined two attributes to define identity classes. Very useful if, for
#' example, you have multiple cell types and multiple replicates, and you want
#' to group cells based on combinations of both.
#'
#' @param object Seurat object
#' @param attribute.1 First attribute for combination. Default is "ident"
#' @param attribute.2 Second attribute for combination. Default is "orig.ident"
#' @return A Seurat object where object@@ident has been appropriately modified
#'
#' @export
#'
#' @examples
#' groups <- sample(c("group1", "group2", "group3"), size = 80, replace = TRUE)
#' celltype <- sample(c("celltype1", "celltype2", "celltype3"), size = 80, replace = TRUE)
#' new.metadata <- data.frame(groups = groups, celltype = celltype)
#' rownames(new.metadata) <- pbmc_small@@cell.names
#' pbmc_small <- AddMetaData(object = pbmc_small, metadata = new.metadata)
#' pbmc_small <- CombineIdent(object = pbmc_small, attribute.1 = "celltype", attribute.2 = "groups")
#' pbmc_small@@ident
#'
CombineIdent <- function(object, attribute.1 = "ident",
                         attribute.2 = "orig.ident") {
  old_data <- FetchData(object = object,vars.all = c(attribute.1, attribute.2))
  new_ids <- sapply(X = 1:nrow(old_data), FUN = function(x){
    paste(as.character(old_data[x, 1]), as.character(old_data[x, 2]), sep = "_")
    })
  object <- SetIdent(object = object,cells.use = object@cell.names, ident.use = new_ids)
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

#' @rdname DownsampleSeurat
#' @exportMethod DownsampleSeurat
#'
setMethod(
  f = 'DownsampleSeurat',
  signature = c('object' = 'seurat'),
  definition = function(object,
                        size,
                        dims.use = 1:10,
                        return.type = "seurat",
                        filename,
                        overwrite = FALSE,
                        ...) {
    if (! "pca" %in% names(object@dr)) {
      stop("PCA not found")
    }
    input.matrix <- GetCellEmbeddings(object = object, reduction.type = "pca")[, dims.use]
    cells.to.keep <- DownsampleMatrix(mat = input.matrix,
                                      size = size,
                                      ...)
    object <- SubsetData(object = object, cells.use = cells.to.keep, subset.raw = TRUE, random.seed = NA)
    if (return.type == "seurat") {
      return(object)
    } else {
      return(Convert(from = object, to = "loom", filename = filename, ...))
    }
  }
)

#' @rdname DownsampleSeurat
#' @exportMethod DownsampleSeurat
#'
setMethod(
  f = 'DownsampleSeurat',
  signature = c('object' = 'loom'),
  definition = function(object,
                        size,
                        dims.use = 1:10,
                        cell.names = "col_attrs/cell_names",
                        keep.layers = FALSE,
                        ...) {
    input.matrix <- GetDimReduction(object = object, reduction.type = 'pca', slot = "cell_embeddings")[, dims.use]
    cells.to.keep <- DownsampleMatrix(mat = input.matrix,
                                      size = size,
                                      ...)
    new.object <- SubsetSeurat(object = object, cells = cells.to.keep, cell.names = cell.names, keep.layers = keep.layers)
    return(new.object)
  }
)

#' @param return.type Return as either a "seurat" or "loom" object
#' @param filename file name/path for new loom object if returning loom
#'
#' @importFrom FNN get.knnx
#' @importFrom Matrix sparseMatrix
#' @importFrom class knn
#'
#' @describeIn ProjectSeurat ...
#' @export ProjectSeurat.Matrix
#' @method ProjectSeurat Matrix
#'
ProjectSeurat.Matrix <- function(
  object,
  template,
  is.scaled.data = FALSE,
  cluster.label = 'res.0.8',
  tsne.proj.k = 11,
  class.k = 11,
  return.type = "seurat",
  filename,
  display.progress = TRUE,
  ...
) {
  if (any(colnames(x = object) %in% template@cell.names)) {
    stop("One or more cell names in UMI matrix match existing cells")
  }
  if (!is.scaled.data) {
    umi.mat <- object
    # Should we filter the UMI matrix here and potentially remove cells?
    # for now assume it has already been filtered and use as is
    
    # todo (maybe): add raw data of new cells to existing raw data
    
    # make a couple of assumptions:
    # LogNormalization was used with factor 1e4
    # genes were centered and scaled
    # pca was used as dimensionality reduction
    # get some information from the existing template
    template.vg <- rownames(template@scale.data)
    template.gene.mean <- apply(template@data[template.vg, ], 1, mean)
    template.gene.sd <- apply(template@data[template.vg, ], 1, sd)
    # log-norm UMI matrix
    n.gene <- apply(umi.mat > 0, 2, sum)
    n.umi <- apply(umi.mat, 2, sum)
    expr <- log1p(x = sweep(
      x = as.matrix(x = umi.mat[template.vg, ]),
      MARGIN = 2,
      STATS = n.umi,
      FUN = '/'
    ) * 10000)
    # scale
    expr <- (expr - template.gene.mean) / template.gene.sd
    # append to scale.data (maybe we don't want to update this slot and will remove this in the future)
    template@scale.data <- cbind(template@scale.data, expr)
    # update meta data
    md <- data.frame(
      matrix(
        data = NA,
        nrow = ncol(x = expr),
        ncol = ncol(template@meta.data),
        dimnames = list(colnames(x = expr), colnames(x = template@meta.data))
      )
    )
    md$nGene <- n.gene
    md$nUMI <- n.umi
    md$projected <- TRUE
  } else {
    expr <- as.matrix(object)
    # update meta data
    md <- data.frame(matrix(NA, ncol(expr), ncol(template@meta.data), dimnames=list(colnames(expr), colnames(template@meta.data))))
    md$projected <- TRUE
  }
  # append to scale.data (maybe we don't want to update this slot and will remove this in the future)
  template@scale.data <- cbind(template@scale.data, expr)
  template@meta.data$projected <- FALSE
  template@meta.data <- rbind(template@meta.data, md)
  template@cell.names <- rownames(x = template@meta.data)
  # apply PCA transformation (blindly assume PCA was run before on template)
  v <- template@dr[['pca']]@gene.loadings
  d <- template@dr[['pca']]@sdev * sqrt(x = nrow(x = template@dr[['pca']]@cell.embeddings) - 1)
  pca.genes <- rownames(x = v)
  template@dr[['pca']]@cell.embeddings <- t(x = template@scale.data[pca.genes, ]) %*% v
  if (!template@calc.params$RunPCA$weight.by.var) {
    template@dr[['pca']]@cell.embeddings <- sweep(
      x = template@dr[['pca']]@cell.embeddings,
      MARGIN = 2,
      STATS = d,
      FUN = '/'
    )
  }
  # map all new cells onto tSNE space
  dims.use <- template@calc.params$RunTSNE$dims.use
  k <- tsne.proj.k
  new.cell <- template@meta.data$projected
  ce <- template@dr[['pca']]@cell.embeddings
  knn.out <- get.knnx(
    data = ce[!new.cell, dims.use],
    query = ce[new.cell, dims.use],
    k = k
  )
  # use weighted average of k nearest neighbors for tsne mapping
  j <- as.numeric(x = t(x = knn.out$nn.index))
  i <- ((1:length(x = j)) - 1) %/% k + 1
  rbf.gamma <- 20
  w <- exp(x = -rbf.gamma * knn.out$nn.dist / knn.out$nn.dist[, k])
  w <- w / rowSums(x = w)
  nn.mat <- sparseMatrix(
    i = i,
    j = j,
    x = as.numeric(x = t(x = w)),
    dims = c(sum(new.cell), sum(!new.cell)),
    giveCsparse = TRUE,
    dimnames = list(template@cell.names[new.cell], template@cell.names[!new.cell])
  )
  tsne.proj <- nn.mat %*% template@dr[['tsne']]@cell.embeddings
  template@dr[['tsne']]@cell.embeddings <- rbind(
    template@dr[['tsne']]@cell.embeddings,
    as.matrix(x = tsne.proj)
  )
  # assign cluster IDs to new cells using knn
  dims.use <- template@calc.params[[sprintf('FindClusters.%s', cluster.label)]]$dims.use
  knn.cluster <- knn(
    train = ce[!new.cell, dims.use],
    test = ce[new.cell, dims.use],
    cl = template@meta.data[[cluster.label]][!new.cell],
    k = class.k,
    l = 0,
    prob = TRUE,
    use.all = TRUE
  )
  template@meta.data[[cluster.label]][new.cell] <- as.character(knn.cluster)
  template@meta.data[[sprintf('%s.proj.prob', cluster.label)]] <- NA
  template@meta.data[[sprintf('%s.proj.prob', cluster.label)]][new.cell] <- attr(x = knn.cluster, which = 'prob')
  template@ident <- factor(template@meta.data[[cluster.label]], ordered=TRUE,
                           levels=as.character(sort(as.numeric(unique(template@meta.data[[cluster.label]])))))
  names(template@ident) <- template@cell.names
  if (return.type == "seurat") {
    return(template)
  } else {
    return(Convert(
      from = template,
      to = "loom",
      filename = filename,
      overwrite = overwrite,
      display.progress = display.progress,
      ...
    ))
  }
}

#' @param gene.names Path to gene names dataset in loom file
#' @param cell.names Path to cell names dataset in loom file
#' @param norm.data Path to normalized dataset in loom file
#' @param scale.data Path to scaled dataset in loom file
#' @param chunk.size ...
#' @param overwrite overwrite loom file if it already exists
#' @param tsne.k ...
#' @param tsne.rbf.gamma ...
#' @param cluster.k ...
#' @param cluster.res ...
#'
#' @importFrom FNN get.knnx
#' @importFrom Matrix sparseMatrix
#' @importFrom class knn
#'
#' @describeIn ProjectSeurat Project to a loom-based dataset
#' @export ProjectSeurat.loom
#' @method ProjectSeurat loom
#'
ProjectSeurat.loom <- function(
  object,
  template,
  gene.names = "row_attrs/gene_names",
  cell.names = "col_attrs/cell_names",
  norm.data = 'layers/norm_data',
  scale.data = 'layers/scale_data',
  chunk.size = 1000,
  overwrite = FALSE,
  tsne.k = 31,
  tsne.rbf.gamma = 20,
  cluster.k = 1,
  cluster.res = 'res.0.8',
  display.progress = TRUE,
  ...
) {
  umi.mat <- object[['matrix']]
  cell.names <- object[[cell.names]][]
  gene.names <- object[[gene.names]][]
  if (display.progress) {
    message('ProjectSeurat: object containing the raw data is loom; seurat object contains instructions for projection')
    message('               raw data loom object will be edited in place')
    message("=============================================")
    message('seurat template contains          ', ncol(x = template@scale.data), ' cells')
    message('cells to be projected             ', umi.mat$dims[1], " cells")
    common.cells <- intersect(x = cell.names, colnames(x = template@scale.data))
    message('cells common to both data sets    ', length(x = common.cells), " cells")
    message("=============================================")
  }
  # log-norm input data
  if (display.progress) {
    message('Normalize data')
  }
  object <- NormalizeData(
    object,
    scale.factor = template@calc.params$NormalizeData$scale.factor,
    overwrite = overwrite
  )
  # scale
  if (display.progress) {
    message('Scale data')
  }
  # get parameters from the seurat template
  pars <- template@calc.params$ScaleData
  #genes.use <- pars$genes.use #rownames(template@scale.data)
  genes.use <- gene.names
  if (pars$do.center) {
    gene.mean <- apply(template@data[genes.use, ], 1, mean)
  } else {
    gene.mean <- rep(0, length(genes.use))
  }
  if (pars$do.scale) {
    gene.sd <- apply(template@data[genes.use, ], 1, sd)
  } else {
    gene.sd <- rep(1, length(genes.use))
  }
  object$apply(
    name = scale.data,
    FUN = function(mat) {
      chunk.scaled <- sweep(sweep(mat, 2, gene.mean, FUN = '-'), 2, gene.sd, FUN = '/')
      chunk.scaled[chunk.scaled > pars$scale.max] <- pars$scale.max
      return(chunk.scaled)
    },
    MARGIN = 2,
    chunk.size = chunk.size,
    dataset.use = norm.data,
    overwrite = overwrite,
    display.progress = display.progress
  )
  # apply PCA transformation saved in template
  if (display.progress) {
    message('Apply PCA transformation')
  }
  v <- template@dr[['pca']]@gene.loadings
  pca.genes <- gene.names %in% rownames(x = v)
  cell.embeddings <- object$map(
    FUN = function(mat) {
      return(mat[, pca.genes] %*% v)
    },
    MARGIN = 2,
    chunk.size = chunk.size,
    dataset.use = scale.data,
    display.progress = display.progress
  )
  if (!template@calc.params$RunPCA$weight.by.var) {
    d <- template@dr[['pca']]@sdev * sqrt(x = nrow(x = template@dr[['pca']]@cell.embeddings) - 1)
    cell.embeddings <- sweep(x = cell.embeddings, MARGIN = 2, STATS = d, FUN = '/')
  }
  colnames(x = cell.embeddings) <- paste0("PC", 1:ncol(x = cell.embeddings))
  cell.embeddings <- as.data.frame(x = cell.embeddings)
  object$add.col.attribute(
    attribute = list('pca_projected_embeddings' = cell.embeddings),
    overwrite = overwrite
  )
  # map all new cells onto tSNE space
  if (display.progress) {
    message('Apply nearest neighbor tSNE mapping; tsne.k is ', tsne.k, '; rsne.rbf.gamma is ', tsne.rbf.gamma)
  }
  dims.use <- template@calc.params$RunTSNE$dims.use
  knn.out <- get.knnx(
    data = template@dr[['pca']]@cell.embeddings[, dims.use],
    query = cell.embeddings[, dims.use],
    k = tsne.k
  )
  # use weighted average of k nearest neighbors for tsne mapping
  j <- as.numeric(x = t(x = knn.out$nn.index))
  i <- ((1:length(x = j)) - 1) %/% tsne.k + 1
  w <- exp(x = -tsne.rbf.gamma * knn.out$nn.dist / knn.out$nn.dist[, tsne.k])
  w <- w / rowSums(x = w)
  nn.mat <- sparseMatrix(
    i = i,
    j = j,
    x = as.numeric(x = t(x = w)),
    dims = c(nrow(x = cell.embeddings), nrow(x = template@dr[['pca']]@cell.embeddings)),
    giveCsparse = TRUE,
    dimnames = list(cell.names, template@cell.names)
  )
  tsne.proj <- as.matrix(x = nn.mat %*% template@dr[['tsne']]@cell.embeddings)
  colnames(x = tsne.proj) <- paste0("tSNE_", 1:ncol(tsne.proj))
  tsne.proj <- as.data.frame(x = tsne.proj)
  object$add.col.attribute(
    attribute = list('tnse_projection' = tsne.proj),
    overwrite = overwrite
  )
  rm(knn.out, j, i, w, nn.mat)
  gc()
  # assign cluster IDs to new cells using knn
  if (display.progress) {
    message('Apply nearest neighbor cluster classification; cluster.k is ', cluster.k, '; cluster.res is ', cluster.res)
  }
  dims.use <- template@calc.params[[paste0('FindClusters.', cluster.res)]]$dims.use
  knn.cluster <- knn(
    train = template@dr[['pca']]@cell.embeddings[, dims.use],
    test = cell.embeddings[, dims.use],
    cl = template@meta.data[[cluster.res]],
    k = cluster.k,
    l = 0,
    prob = FALSE,
    use.all = TRUE
  )
  knn.cluster <- factor(
    x = knn.cluster,
    ordered = TRUE,
    levels = as.character(x = sort(x = as.numeric(x = levels(x = knn.cluster))))
  )
  clustering <- list(ident = knn.cluster)
  clustering[[cluster.res]] <- knn.cluster
  object$add.col.attribute(attribute = clustering, overwrite = overwrite)
  invisible(x = object)
}

#' @rdname SubsetSeurat
#' @exportMethod SubsetSeurat
#'
setMethod(
  f = 'SubsetSeurat',
  signature = c('object' = 'loom'),
  definition = function(object,
                        cells,
                        return.type = "seurat",
                        cell.names = "col_attrs/cell_names",
                        gene.names = "row_attrs/gene_names",
                        filename = NULL,
                        remove.existing = FALSE,
                        display.progress = TRUE,
                        chunk.size = 1000,
                        keep.layers = FALSE,
                        norm.data = 'layers/norm_data',
                        scale.data = 'layers/scale_data',
                        ...) {
    cell.names <- object[[cell.names]][]
    gene.names <- object[[gene.names]][]
    cells.use <- which(cell.names %in% cells)
    n.cells <- length(cells.use)
    n.genes <- length(gene.names)
    cat('Subsetting loom object; will create matrix of size', n.cells, 'x', n.genes, '(cells x genes) in memory\n')
    cat('Note that original ordering of cells will be preserved regardless of order of cell names in second parameter\n')
    mat <- GetAssayData.loom(object=object, slot='raw.data', cells.use=cells.use, chunk.size=chunk.size)
    
    # get meta data
    meta.data <- object$get.attribute.df(attribute.layer = 'col')[cells.use, ]
    gene.attrs <- object$get.attribute.df(attribute.layer = 'row')
    
    if (return.type == "seurat") {
      new.object <- CreateSeuratObject(raw.data = as(mat, "dgCMatrix"))
      duplicate.cols <- colnames(meta.data) %in% colnames(new.object@meta.data)
      colnames(meta.data)[duplicate.cols] <- paste(colnames(meta.data)[duplicate.cols], 'loom', sep='.')
      new.object <- AddMetaData(object = new.object, metadata = meta.data)
      if (!keep.layers) {
        cat('Note: SubsetSeurat is converting from loom to seurat; keep.layers=FALSE, so only raw data and cell attributes are kept, no layers\n')
      } else {
        cat('Also subsetting normalized data and scaled data\n')
        new.object@data <- GetAssayData.loom(object=object, slot=norm.data, cells.use=cells.use, chunk.size=chunk.size)
        new.object@scale.data <- GetAssayData.loom(object=object, slot=scale.data, cells.use=cells.use, chunk.size=chunk.size)
      }
      return(new.object)
    } else {
      if (file.exists(filename)) {
        if (remove.existing) {
          file.remove(filename)
        } else {
          stop('file ', filename, ' exists and remove.existing is set to FALSE\n')
        }
      }
      loomfile <- create(
        filename = filename,
        data = mat,
        cell.attrs = meta.data[, setdiff(colnames(meta.data), 'cell_names')],
        gene.attrs = gene.attrs[, setdiff(colnames(gene.attrs), 'gene_names')], ...
      )
      if (keep.layers) {
        for (layer.name in names(object[['layers']])) {
          layer.path <- paste('layers', layer.name, sep='/')
          cat('layer', layer.path, '\n')
          mat <- GetAssayData.loom(object=object, slot=layer.path, cells.use=cells.use, chunk.size=chunk.size)
          lst <- list()
          lst[[layer.name]] <- mat
          loomfile$add.layer(lst)
        }
      }
      return(loomfile)
    }

  }
)
