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
#' pbmc_small
#' pbmc_small <- RunDiffusion(pbmc_small, genes.use = pbmc_small@var.genes)
#' head(DMEmbed(object = pbmc_small))
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

################### Top Genes/Cells Related Functions ##########################

#' Find genes with highest scores for a given dimensional reduction technique
#'
#' Return a list of genes with the strongest contribution to a set of components
#'
#' @param object Seurat object
#' @param dim.use Dimension to use
#' @param reduction.type Dimensional reduction to find the highest score for
#' @param num.genes Number of genes to return
#' @param use.full Use the full PCA (projected PCA). Default i s FALSE
#' @param do.balanced Return an equal number of genes with both + and - scores.
#'
#' @return Returns a vector of genes
#'
#' @export
#'
#' @examples
#' pbmc_small
#' DimTopGenes(object = pbmc_small, dim.use = 1, reduction.type = "pca")
#' # After projection:
#' DimTopGenes(object = pbmc_small, dim.use = 1, reduction.type = "pca", use.full = TRUE)
#'
DimTopGenes <- function(
  object,
  dim.use = 1,
  reduction.type = "pca",
  num.genes = 30,
  use.full = FALSE,
  do.balanced = FALSE
) {
  #note that we use GetTopGenes, but it still works
  #error checking
  if (! reduction.type %in% names(x = object@dr)) {
    stop(paste(reduction.type, "dimensional reduction has not been computed"))
  }
  dim.scores <- GetDimReduction(
    object = object,
    reduction.type = reduction.type,
    slot = "gene.loadings"
  )
  if (use.full) {
    dim.scores <- GetDimReduction(
      object = object,
      reduction.type = reduction.type,
      slot = "gene.loadings.full"
    )
  }
  if ((is.null(x = dim.scores)) || (ncol(x = dim.scores) < 2)) {
    stop(paste0(
      "Gene loadings for ",
      reduction.type,
      " with use.full=",
      use.full,
      " have not been computed"
    ))
  }
  i <- dim.use
  num.genes <- min(num.genes, length(x = rownames(x = dim.scores)))
  key <- GetDimReduction(
    object = object,
    reduction.type = reduction.type,
    slot = "key"
  )
  dim.top.genes <- unique(x = unlist(x = lapply(
    X = i,
    FUN = GetTopGenes,
    dim.scores = dim.scores,
    do.balanced = do.balanced,
    num.genes = num.genes
  )))
  return(dim.top.genes)
}

#' Find genes with highest PCA scores
#'
#' Return a list of genes with the strongest contribution to a set of principal
#' components
#'
#' @param object Seurat object
#' @param pc.use Principal components to use
#' @param num.genes Number of genes to return
#' @param use.full Use the full PCA (projected PCA). Default is FALSE
#' @param do.balanced Return an equal number of genes with both + and - PC scores.
#'
#' @return Returns a vector of genes
#'
#' @export
#'
#' @examples
#' pbmc_small
#' PCTopGenes(object = pbmc_small, pc.use = 1)
#' # After projection:
#' PCTopGenes(object = pbmc_small, pc.use = 1, use.full = TRUE)
#'
PCTopGenes <- function(
  object,
  pc.use = 1,
  num.genes = 30,
  use.full = FALSE,
  do.balanced = FALSE
) {
  return(DimTopGenes(
    object = object,
    dim.use = pc.use,
    reduction.type = "pca",
    num.genes = num.genes,
    use.full = use.full,
    do.balanced = do.balanced
  ))
}

#' Find genes with highest ICA scores
#'
#' Return a list of genes with the strongest contribution to a set of
#' indepdendent components
#'
#' @param object Seurat object
#' @param ic.use Independent components to use
#' @param num.genes Number of genes to return
#' @param use.full Use the full ICA (projected ICA), default is FALSE
#' @param do.balanced Return an equal number of genes with both + and - IC scores.
#'
#' @return Returns a vector of genes
#'
#' @export
#'
#' @examples
#' pbmc_small
#' pbmc_small <- RunICA(object = pbmc_small, ics.compute = 10, ics.print = 0)
#' pbmc_small <- ProjectDim(object = pbmc_small, reduction.type = "ica", do.print = FALSE)
#' ICTopGenes(object = pbmc_small, ic.use = 1)
#' # After projection:
#' ICTopGenes(object = pbmc_small, ic.use = 1, use.full = TRUE)
#'
ICTopGenes <- function(
  object,
  ic.use = 1,
  num.genes = 30,
  use.full = FALSE,
  do.balanced = FALSE
) {
  return(DimTopGenes(
    object = object,
    dim.use = ic.use,
    reduction.type = "ica",
    use.full = use.full,
    num.genes = num.genes,
    do.balanced = do.balanced
  ))
}

#' Find cells with highest scores for a given dimensional reduction technique
#'
#' Return a list of genes with the strongest contribution to a set of components
#'
#' @param object Seurat object
#' @param reduction.type Dimensional reduction to find the highest score for
#' @param dim.use Components to use
#' @param num.cells Number of cells to return
#' @param do.balanced Return an equal number of cells with both + and - scores.
#'
#' @return Returns a vector of cells
#'
#' @export
#'
#' @examples
#' pbmc_small
#' head(DimTopCells(object = pbmc_small, reduction.type = "pca"))
#' # Can specify which dimension and how many cells to return
#' DimTopCells(object = pbmc_small, reduction.type = "pca", dim.use = 2, num.cells = 5)
#'
DimTopCells <- function(
  object,
  dim.use = 1,
  reduction.type = "pca",
  num.cells = NULL,
  do.balanced = FALSE
) {
  #note that we use GetTopGenes, but it still works
  #error checking
  if (! reduction.type %in% names(x = object@dr)) {
    stop(paste(reduction.type, "dimensional reduction has not been computed"))
  }
  num.cells <- SetIfNull(x = num.cells, default = length(x = object@cell.names))
  dim.scores <- GetDimReduction(
    object = object,
    reduction.type = reduction.type,
    slot = "cell.embeddings"
  )
  key <- GetDimReduction(
    object = object,
    reduction.type = reduction.type,
    slot = "key"
  )
  i <- dim.use
  dim.top.cells <- unique(x = unlist(x = lapply(
    X = i,
    FUN = GetTopGenes,
    dim.scores = dim.scores,
    do.balanced = do.balanced,
    num.genes = num.cells
  )))
  return(dim.top.cells)
}

#' Find cells with highest PCA scores
#'
#' Return a list of genes with the strongest contribution to a set of principal components
#'
#' @param object Seurat object
#' @param pc.use Principal component to use
#' @param num.cells Number of cells to return
#' @param do.balanced Return an equal number of cells with both + and - PC scores.
#'
#' @return Returns a vector of cells
#'
#' @export
#'
#' @examples
#' pbmc_small
#' head(PCTopCells(object = pbmc_small))
#' # Can specify which dimension and how many cells to return
#' DimTopCells(object = pbmc_small, dim.use = 2, num.cells = 5)
#'
PCTopCells <- function(
  object,
  pc.use = 1,
  num.cells = NULL,
  do.balanced = FALSE
) {
  return(DimTopCells(
    object = object,
    dim.use = pc.use,
    reduction.type = "pca",
    num.cells = num.cells,
    do.balanced = do.balanced
  ))
}

#' Find cells with highest ICA scores
#'
#' Return a list of genes with the strongest contribution to a set of principal
#' components
#'
#' @param object Seurat object
#' @param ic.use Independent component to use
#' @param num.cells Number of cells to return
#' @param do.balanced Return an equal number of cells with both + and - PC scores.
#'
#' @return Returns a vector of cells
#'
#' @export
#'
#' @examples
#' pbmc_small
#' pbmc_small <- RunICA(object = pbmc_small, ics.compute = 10, ics.print = 0)
#' pbmc_small <- ProjectDim(object = pbmc_small, reduction.type = "ica", do.print = FALSE)
#' ICTopCells(object = pbmc_small)
#' # Can specify which dimension and how many cells to return
#' ICTopCells(object = pbmc_small, ic.use = 2, num.cells = 5)
#'
ICTopCells <- function(
  object,
  ic.use = 1,
  num.cells = NULL,
  do.balanced = FALSE
) {
  return(DimTopCells(
    object = object,
    dim.use = ic.use,
    reduction.type = "ica",
    num.cells = num.cells,
    do.balanced = do.balanced
  ))
}

##################### Printing results #########################################

#' Print the results of a dimensional reduction analysis
#'
#' Prints a set of genes that most strongly define a set of components
#'
#' @param object Seurat object
#' @param reduction.type Reduction technique to print results for
#' @param dims.print Number of dimensions to display
#' @param genes.print Number of genes to display
#' @param use.full Use full PCA (i.e. the projected PCA, by default FALSE)
#'
#' @return Set of genes defining the components
#'
#' @export
#'
#' @examples
#' pbmc_small
#' PrintDim(object = pbmc_small, reduction.type = "pca")
#' # Options for how many dimensions and how many genes to print
#' PrintDim(object = pbmc_small, reduction.type = "pca", dims.print = 1:2, genes.print = 5)
#' # Can also print for the projected PCA
#' PrintDim(object = pbmc_small, reduction.type = "pca", use.full = TRUE)
#'
PrintDim <- function(
  object,
  reduction.type = "pca",
  dims.print = 1:5,
  genes.print = 30,
  use.full = FALSE
) {
  if (use.full) {
    slot.use <- "gene.loadings.full"
  } else {
    slot.use <- "gene.loadings"
  }
  dim.scores <- GetDimReduction(
    object = object,
    reduction.type = reduction.type,
    slot = slot.use
  )
  dim.prefix <- GetDimReduction(
    object = object,
    reduction.type = reduction.type,
    slot = "key"
  )
  dim.codes.exist <- colnames(x = dim.scores)
  dim.codes.input <- paste0(dim.prefix, dims.print)
  dims.print <- dims.print[which(x = dim.codes.input %in% dim.codes.exist)]
  genes.print <- min(genes.print, nrow(x = dim.scores))
  if (length(x = dim.scores) == 0 && use.full) {
    warning("Dimensions have not been projected. Setting use.full = FALSE")
    use.full <- FALSE
  }
  for (i in dims.print) {
    code <- paste0(
      GetDimReduction(
        object = object,
        reduction.type = reduction.type,
        slot = "key"
      ),
      i
    )
    sx <- DimTopGenes(
      object = object,
      dim.use = i,
      reduction.type = reduction.type,
      num.genes = genes.print * 2,
      use.full = use.full,
      do.balanced = TRUE
    )
    print(code)
    print((sx[1:genes.print]))
    print ("")
    print(rev(x = (sx[(length(x = sx) - genes.print + 1):length(x = sx)])))
    print ("")
    print ("")
  }
}

#' Print the results of a ICA analysis
#'
#' Prints a set of genes that most strongly define a set of independent components
#'
#' @inheritParams VizPCA
#' @param ics.print Set of ICs to print genes for
#' @param genes.print Number of genes to print for each PC
#'
#' @return Only text output
#'
#' @export
#'
#' @examples
#' pbmc_small
#' pbmc_small <- RunICA(object = pbmc_small, ics.compute = 10, ics.print = 0)
#' pbmc_small <- ProjectDim(object = pbmc_small, reduction.type = "ica", do.print = FALSE)
#' PrintICA(object = pbmc_small)
#' # Options for how many dimensions and how many genes to print
#' PrintICA(object = pbmc_small, ics.print = 1:2, genes.print = 5)
#' # Can also print for the projected PCA
#' PrintICA(object = pbmc_small, use.full = TRUE)
#'
PrintICA <- function(
  object,
  ics.print = 1:5,
  genes.print = 30,
  use.full = FALSE
) {
  PrintDim(
    object = object,
    reduction.type = "ica",
    dims.print = ics.print,
    genes.print = genes.print,
    use.full = use.full
  )
}

#' Print the results of a PCA analysis
#'
#' Prints a set of genes that most strongly define a set of principal components
#'
#' @inheritParams VizPCA
#' @param pcs.print Set of PCs to print genes for
#' @param genes.print Number of genes to print for each PC
#'
#' @return Only text output
#'
#' @export
#'
#' @examples
#' pbmc_small
#' PrintPCA(object = pbmc_small)
#' # Options for how many dimensions and how many genes to print
#' PrintPCA(object = pbmc_small, pcs.print = 1:2, genes.print = 5)
#' # Can also print for the projected PCA
#' PrintPCA(object = pbmc_small, use.full = TRUE)
#'
PrintPCA <- function(
  object,
  pcs.print = 1:5,
  genes.print = 30,
  use.full = FALSE
) {
  PrintDim(
    object = object,
    reduction.type = "pca",
    dims.print = pcs.print,
    genes.print = genes.print,
    use.full = use.full
  )
}
