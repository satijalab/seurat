######################## Accessor/Mutator Functions ############################

#' Dimensional Reduction Accessor Function
#'
#' General accessor function for dimensional reduction objects. Pulls slot
#' contents for specified stored dimensional reduction analysis.
#'
#' @param object Seurat object
#' @param reduction.type Type of dimensional reduction to fetch (default is PCA)
#' @param slot Specific information to pull (must be one of the following:
#'  "cell.embeddings", "gene.loadings", "gene.loadings.full", "sdev", "key", "misc")
#'
#' @return Returns specified slot results from given reduction technique
#'
#' @importFrom methods slotNames
#'
#' @export
#'
#' @examples
#' pbmc_small
#' # Get the PCA cell embeddings and print the top left corner
#' GetDimReduction(object = pbmc_small, reduction.type = "pca",
#'                 slot = "cell.embeddings")[1:5, 1:5]
#' # Get the standard deviation of each PC
#' GetDimReduction(object = pbmc_small, reduction.type = "pca", slot = "sdev")
#'
GetDimReduction <- function(
  object,
  reduction.type = "pca",
  slot = "gene.loadings"
) {
  if (! (reduction.type %in% names(object@dr))) {
    stop(paste(reduction.type, " dimensional reduction has not been computed"))
  }
  reduction <- paste0("object@dr$", reduction.type)
  reduction.slots <- slotNames(x = eval(expr = parse(text = reduction)))
  if (! (slot %in% reduction.slots)) {
    stop(paste0(slot, " slot doesn't exist"))
  }
  return(eval(expr = parse(text = paste0(reduction, "@", slot))))
}

#' Dimensional Reduction Cell Embeddings Accessor Function
#'
#' Pull cell embeddings matrix for specified stored dimensional reduction
#' analysis
#'
#' @param object Seurat object
#' @param reduction.type Type of dimensional reduction to fetch (default is PCA)
#' @param dims.use Dimensions to include (default is all stored dims)
#' @param cells.use Cells to include (default is all cells)
#'
#' @return Cell embedding matrix for given reduction, cells, and dimensions
#'
#' @export
#'
#' @examples
#' pbmc_small
#' # Examine the head of the first 5 PC cell embeddings
#' head(GetCellEmbeddings(object = pbmc_small, reduction.type = "pca", dims.use = 1:5))
#'
GetCellEmbeddings <- function(
  object,
  reduction.type = "pca",
  dims.use = NULL,
  cells.use = NULL
) {
  object.embed <- GetDimReduction(
    object = object,
    reduction.type = reduction.type,
    slot = "cell.embeddings"
  )
  if (length(x = object.embed) == 0) {
    stop(paste0("Cell embeddings slot for ", reduction.type, " is empty."))
  }
  cells.use <- SetIfNull(x = cells.use, default = rownames(x = object.embed))
  if (any(! cells.use %in% rownames(x = object.embed))) {
    missing.cells <- paste0(
      cells.use[which(x = ! cells.use %in% rownames(x = object.embed))],
      collapse = ", "
    )
    warning(paste0("Could not find the following cell names: ", missing.cells))
    cells.use <- intersect(x = cells.use, y = rownames(x = object.embed))
  }
  dims.use <- SetIfNull(x = dims.use, default = 1:ncol(x = object.embed))
  if (any(!dims.use %in% 1:ncol(x = object.embed))) {
    missing.dims <- paste0(
      dims.use[which(x = ! dims.use %in% 1:ncol(x = object.embed))],
      collapse = ", "
    )
    stop(paste0("Could not find the following dimensions: ", missing.dims))
  }
  object.embed <- object.embed[cells.use, dims.use, drop = FALSE]
  object.key <- GetDimReduction(
    object = object,
    reduction.type = reduction.type,
    slot = "key"
  )
  if (length(x = object.key) == 0) {
    colnames(x = object.embed) <- NULL
  } else {
    colnames(x = object.embed) <- paste0(object.key, dims.use)
  }
  return(object.embed)
}

#' Dimensional Reduction Gene Loadings Accessor Function
#'
#' Pull gene loadings matrix for specified stored dimensional reduction analysis.
#'
#' @param object Seurat object
#' @param reduction.type Type of dimensional reduction to fetch (default is PCA)
#' @param dims.use Dimensions to include (default is all stored dims)
#' @param genes.use Genes to include (default is all genes)
#' @param use.full Return projected gene loadings (default is FALSE)
#' @return Gene loading matrix for given reduction, cells, and genes
#' @export
#'
#' @examples
#' pbmc_small
#' # Examine the head of the first 5 PC gene loadings
#' head(GetGeneLoadings(object = pbmc_small, reduction.type = "pca", dims.use = 1:5))
#'
GetGeneLoadings <- function(
  object,
  reduction.type = "pca",
  dims.use = NULL,
  genes.use = NULL,
  use.full = FALSE
) {
  if (use.full) {
    gene.loadings <- GetDimReduction(
      object = object,
      reduction.type = reduction.type,
      slot = "gene.loadings.full"
    )
  } else {
    gene.loadings <- GetDimReduction(
      object = object,
      reduction.type = reduction.type,
      slot = "gene.loadings"
    )
  }
  if (length(x = gene.loadings) == 0) {
    stop(paste("gene loadings slot for", reduction.type, "is empty."))
  }
  genes.use <- SetIfNull(x = genes.use, default = rownames(x = gene.loadings))
  if (any(! genes.use %in% rownames(x = gene.loadings))) {
    missing.genes <- paste0(
      genes.use[which(x = ! genes.use %in% rownames(x = gene.loadings))],
      collapse = ", "
    )
    warning(paste("Could not find the following gene names:", missing.genes))
    genes.use <- intersect(x = genes.use, y = rownames(x = gene.loadings))
  }
  dims.use <- SetIfNull(x = dims.use, default = 1:ncol(x = gene.loadings))
  if (any(! dims.use %in% 1:ncol(x = gene.loadings))) {
    missing.dims <- paste0(
      dims.use[which(x = ! dims.use %in% 1:ncol(x = gene.loadings))],
      collapse = ", "
    )
    stop(paste("Could not find the following dimensions:", missing.dims))
  }
  gene.loadings <- gene.loadings[genes.use, dims.use, drop = FALSE]
  object.key <- GetDimReduction(
    object = object,
    reduction.type = reduction.type,
    slot = "key"
  )
  if (length(x = object.key) == 0) {
    colnames(x = gene.loadings) <- NULL
  } else {
    colnames(x = gene.loadings) <- paste0(object.key, dims.use)
  }
  return(gene.loadings)
}

#' Dimensional Reduction Mutator Function
#'
#' Set information for specified stored dimensional reduction analysis
#'
#' @param object Seurat object
#' @param reduction.type Type of dimensional reduction to set
#' @param slot Specific information to set (must be one of the following:
#' "cell.embeddings", "gene.loadings", "gene.loadings.full", "sdev", "key",
#' "misc")
#' @param new.data New data to set
#'
#' @return Seurat object with updated slot
#'
#' @importFrom methods new
#'
#' @export
#'
#' @examples
#' pbmc_small
#' # Simulate adding a new dimensional reduction
#' new.cell.embeddings <- GetCellEmbeddings(object = pbmc_small, reduction.type = "pca")
#' new.gene.loadings <- GetGeneLoadings(object = pbmc_small, reduction.type = "pca")
#' SetDimReduction(
#'   object = pbmc_small,
#'   reduction.type = "new.pca",
#'   slot = "cell.embeddings",
#'   new.data = new.cell.embeddings
#' )
#' SetDimReduction(
#'   object = pbmc_small,
#'   reduction.type = "new.pca",
#'   slot = "gene.loadings",
#'   new.data = new.gene.loadings
#' )
#'
SetDimReduction <- function(
  object,
  reduction.type,
  slot,
  new.data
) {
  if (reduction.type %in% names(x = object@dr)) {
    eval(expr = parse(text = paste0(
      "object@dr$",
      reduction.type,
      "@",
      slot,
      "<- new.data"
    )))
    if (slot == "key") {
      cell.embeddings=GetCellEmbeddings(object = object,reduction.type = reduction.type)
      colnames(cell.embeddings)=paste(new.data,1:ncol(cell.embeddings),sep="")
      object=SetDimReduction(object,reduction.type = reduction.type, slot = "cell.embeddings", new.data = cell.embeddings)

      gene.loadings <- GetDimReduction(
        object = object,
        reduction.type = reduction.type,
        slot = "gene.loadings"
      )
      gene.loadings.full <- GetDimReduction(
        object = object,
        reduction.type = reduction.type,
        slot = "gene.loadings.full"
      )
      if (length(gene.loadings > 0)) {
        colnames(gene.loadings)=paste(new.data,1:ncol(gene.loadings),sep="")
        object=SetDimReduction(object,reduction.type = reduction.type, slot = "gene.loadings", new.data = gene.loadings)
      }
      if (length(gene.loadings.full > 0)) {
        colnames(gene.loadings.full)=paste(new.data,1:ncol(gene.loadings.full),sep="")
        object=SetDimReduction(object,reduction.type = reduction.type, slot = "gene.loadings.full", new.data = gene.loadings.full)
      }
    }
  } else {
    new.dr <- new(Class = "dim.reduction")
    eval(expr = parse(text = paste0("new.dr@", slot, "<- new.data")))
    eval(expr = parse(text = paste0("object@dr$", reduction.type, "<- new.dr")))
  }
  return(object)
}

################### Convienence functions for easy interaction #################

#' Diffusion Maps Cell Embeddings Accessor Function
#'
#' Pull Diffusion maps cell embedding matrix
#'
#' @param object Seurat object
#' @param dims.use Dimensions to include (default is all stored dims)
#' @param cells.use Cells to include (default is all cells)
#'
#' @return Diffusion maps embedding matrix for given cells and DMs
#'
#' @export
#'
#' @examples
#' \dontrun{
#' pbmc_small
#' pbmc_small <- RunDiffusion(pbmc_small, genes.use = pbmc_small@var.genes)
#' head(DMEmbed(object = pbmc_small))
#' }
#'
DMEmbed <- function(
  object,
  dims.use = NULL,
  cells.use = NULL
) {
  return(GetCellEmbeddings(
    object = object,
    reduction.type = "dm",
    dims.use = dims.use,
    cells.use = cells.use
  ))
}

#' PCA Cell Embeddings Accessor Function
#'
#' Pull PCA cell embedding matrix
#'
#' @param object Seurat object
#' @param dims.use Dimensions to include (default is all stored dims)
#' @param cells.use Cells to include (default is all cells)
#'
#' @return PCA cell embedding matrix for given cells and PCs
#'
#' @export
#'
#' @examples
#' pbmc_small
#' head(PCAEmbed(pbmc_small))
#' # Optionally, you can specify subsets of dims or cells to use
#' PCAEmbed(pbmc_small, dims.use = 1:5, cells.use = pbmc_small@cell.names[1:5])
#'
PCAEmbed <- function(
  object,
  dims.use = NULL,
  cells.use = NULL
) {
  return(GetCellEmbeddings(
    object = object,
    reduction.type = "pca",
    dims.use = dims.use,
    cells.use = cells.use
  ))
}

#' ICA Cell Embeddings Accessor Function
#'
#' Pull ICA cell embeddings matrix
#'
#' @param object Seurat object
#' @param dims.use Dimensions to include (default is all stored dims)
#' @param cells.use Cells to include (default is all cells)
#'
#' @return ICA cell embeddings matrix for given cells and ICs
#'
#' @export
#'
#' @examples
#' pbmc_small
#' pbmc_small <- RunICA(pbmc_small, ics.compute = 10, ics.print = 0)
#' head(ICAEmbed(pbmc_small))
#' # Optionally, you can specify subsets of dims or cells to use
#' ICAEmbed(pbmc_small, dims.use = 1:5, cells.use = pbmc_small@cell.names[1:5])
#'
ICAEmbed <- function(
  object,
  dims.use = NULL,
  cells.use = NULL
) {
  return(GetCellEmbeddings(
    object = object,
    reduction.type = "ica",
    dims.use = dims.use,
    cells.use = cells.use
  ))
}

#' PCA Gene Loadings Accessor Function
#'
#' Pull the PCA gene loadings matrix
#'
#' @param object Seurat object
#' @param dims.use Dimensions to include (default is all stored dims)
#' @param genes.use Genes to include (default is all genes)
#' @param use.full Return projected gene loadings (default is FALSE)
#'
#' @return PCA gene loading matrix for given genes and PCs
#'
#' @export
#'
#' @examples
#' pbmc_small
#' head(PCALoad(pbmc_small))
#' # Optionally, you can specify subsets of dims or genes to use
#' PCALoad(pbmc_small, dims.use = 1:5, genes.use = pbmc_small@var.genes[1:5])
#'
PCALoad <- function(
  object,
  dims.use = NULL,
  genes.use = NULL,
  use.full = FALSE
) {
  return(GetGeneLoadings(
    object = object,
    reduction.type = "pca",
    dims.use = dims.use,
    genes.use = genes.use,
    use.full = use.full
  ))
}

#' ICA Gene Loadings Accessor Function
#'
#' Pull the ICA gene loadings matrix
#'
#' @param object Seurat object
#' @param dims.use Dimensions to include (default is all stored dims)
#' @param genes.use Genes to include (default is all)
#' @param use.full Return projected gene loadings (default is FALSE)
#'
#' @return ICA gene loading matrix for given genes and ICs
#'
#' @export
#'
#' @examples
#' pbmc_small
#' pbmc_small <- RunICA(pbmc_small, ics.compute = 10, ics.print = 0)
#' head(ICALoad(pbmc_small))
#' # Optionally, you can specify subsets of dims or cells to use
#' ICALoad(pbmc_small, dims.use = 1:5, genes.use = pbmc_small@var.genes[1:5])
#'
ICALoad <- function(
  object,
  dims.use = NULL,
  genes.use = NULL,
  use.full = FALSE
) {
  return(GetGeneLoadings(
    object = object,
    reduction.type = "ica",
    dims.use = dims.use,
    genes.use = genes.use,
    use.full = use.full
  ))
}
