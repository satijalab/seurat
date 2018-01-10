#' Run Principal Component Analysis on gene expression using IRLBA
#'
#' Run a PCA dimensionality reduction. For details about stored PCA calculation
#' parameters, see \code{PrintPCAParams}.
#'
#' @param object Seurat object
#' @param pc.genes Genes to use as input for PCA. Default is object@@var.genes
#' @param pcs.compute Total Number of PCs to compute and store (20 by default)
#' @param use.imputed Run PCA on imputed values (FALSE by default)
#' @param rev.pca By default computes the PCA on the cell x gene matrix. Setting
#' to true will compute it on gene x cell matrix.
#' @param weight.by.var Weight the cell embeddings by the variance of each PC
#' (weights the gene loadings if rev.pca is TRUE)
#' @param do.print Print the top genes associated with high/low loadings for
#' the PCs
#' @param pcs.print PCs to print genes for
#' @param genes.print Number of genes to print for each PC
#' @param reduction.name dimensional reduction name, specifies the position in the object$dr list. pca by default
#' @param reduction.key dimensional reduction key, specifies the string before the number for the dimension names. PC by default
#' @param assay.type Data type, RNA by default. Can be changed for multimodal
#' @param seed.use Set a random seed. By default, sets the seed to 42. Setting
#' NULL will not set a seed.
#' @param \dots Additional arguments to be passed to IRLBA
#'
#'@importFrom irlba irlba
#'
#' @return Returns Seurat object with the PCA calculation stored in
#' object@@dr$pca.
#'
#' @importFrom irlba irlba
#'
#' @export
#'
#' @examples
#' pbmc_small
#' # Run PCA on variable genes (default)
#' pbmc_small <- RunPCA(pbmc_small)
#' # Run PCA on different gene set (in this case all genes)
#' pbmc_small=RunPCA(pbmc_small,pc.genes = rownames(pbmc_small@data))
#' # Run PCA but compute more than 20 dimensions
#' pbmc_small=RunPCA(pbmc_small,pcs.compute=30)
#' # Plot results
#' PCAPlot(pbmc_small)
#'
RunPCA <- function(
  object,
  pc.genes = NULL,
  pcs.compute = 20,
  use.imputed = FALSE,
  rev.pca = FALSE,
  weight.by.var = TRUE,
  do.print = TRUE,
  pcs.print = 1:5,
  genes.print = 30,
  reduction.name = "pca",
  reduction.key = "PC",
  assay.type="RNA",
  seed.use = 42,
  ...
) {
  if (!is.null(seed.use)) {
    set.seed(seed = seed.use)
  }
  data.use <- PrepDR(
    object = object,
    genes.use = pc.genes,
    use.imputed = use.imputed,
    assay.type = assay.type)
  if (rev.pca) {
    pcs.compute <- min(pcs.compute, ncol(x = data.use)-1)
    pca.results <- irlba(A = data.use, nv = pcs.compute, ...)
    sdev <- pca.results$d/sqrt(max(1, nrow(data.use) - 1))
    if(weight.by.var){
      gene.loadings <- pca.results$u %*% diag(pca.results$d)
    } else{
      gene.loadings <- pca.results$u
    }
    cell.embeddings <- pca.results$v
  }
  else {
    pcs.compute <- min(pcs.compute, nrow(x = data.use)-1)
    pca.results <- irlba(A = t(x = data.use), nv = pcs.compute, ...)
    gene.loadings <- pca.results$v
    sdev <- pca.results$d/sqrt(max(1, ncol(data.use) - 1))
    if(weight.by.var){
      cell.embeddings <- pca.results$u %*% diag(pca.results$d)
    } else {
      cell.embeddings <- pca.results$u
    }
  }
  rownames(x = gene.loadings) <- rownames(x = data.use)
  colnames(x = gene.loadings) <- paste0(reduction.key, 1:pcs.compute)
  rownames(x = cell.embeddings) <- colnames(x = data.use)
  colnames(x = cell.embeddings) <- colnames(x = gene.loadings)
  pca.obj <- new(
    Class = "dim.reduction",
    gene.loadings = gene.loadings,
    cell.embeddings = cell.embeddings,
    sdev = sdev,
    key = reduction.key
  )
  #object@dr[reduction.name] <- pca.obj
  eval(expr = parse(text = paste0("object@dr$", reduction.name, "<- pca.obj")))

  parameters.to.store <- as.list(environment(), all = TRUE)[names(formals("RunPCA"))]
  object <- SetCalcParams(object = object, calculation = "RunPCA", ... = parameters.to.store)
  if(is.null(object@calc.params$RunPCA$pc.genes)){
    object@calc.params$RunPCA$pc.genes <- rownames(data.use)
  }
  if(do.print){
    PrintDim(object = object, dims.print = pcs.print, genes.print = genes.print,reduction.type = reduction.name)
  }
  return(object)
}

#' Run Independent Component Analysis on gene expression
#'
#' Run fastica algorithm from the ica package for ICA dimensionality reduction.
#' For details about stored ICA calculation parameters, see
#' \code{PrintICAParams}.
#'
#' @param object Seurat object
#' @param ic.genes Genes to use as input for ICA. Default is object@@var.genes
#' @param ics.compute Number of ICs to compute
#' @param use.imputed Run ICA on imputed values (FALSE by default)
#' @param rev.ica By default, computes the dimensional reduction on the cell x
#' gene matrix. Setting to true will compute it on the transpose (gene x cell
#' matrix).
#' @param print.results Print the top genes associated with each dimension
#' @param ics.print ICs to print genes for
#' @param genes.print Number of genes to print for each IC
#' @param ica.function ICA function from ica package to run (options: icafast,
#' icaimax, icajade)
#' @param seed.use Random seed to use for fastica
#' @param reduction.name dimensional reduction name, specifies the position in the object$dr list. ica by default
#' @param reduction.key dimensional reduction key, specifies the string before the number for the dimension names. IC by default
#' @param \dots Additional arguments to be passed to fastica
#'
#' @importFrom ica icafast icaimax icajade
#'
#' @return Returns Seurat object with an ICA calculation stored in
#' object@@dr$ica
#'
#' @export
#'
#' @examples
#' pbmc_small
#' # Run ICA on variable genes (default)
#' pbmc_small <- RunICA(pbmc_small, ics.compute=5)
#' # Run ICA on different gene set (in this case all genes)
#' pbmc_small <- RunICA(pbmc_small, ic.genes = rownames(pbmc_small@data))
#' # Plot results
#' ICAPlot(pbmc_small)
#'
RunICA <- function(
  object,
  ic.genes = NULL,
  ics.compute = 50,
  use.imputed = FALSE,
  rev.ica = FALSE,
  print.results = TRUE,
  ics.print = 1:5,
  genes.print = 50,
  ica.function = "icafast",
  seed.use = 1,
  reduction.name = "ica",
  reduction.key = "IC",
  ...
) {
  data.use <- PrepDR(
    object = object,
    genes.use = ic.genes,
    use.imputed = use.imputed)
  set.seed(seed = seed.use)
  ics.compute <- min(ics.compute, ncol(x = data.use))
  ica.fxn <- eval(parse(text = ica.function))
  if (rev.ica) {
    ica.results <- ica.fxn(data.use, nc = ics.compute,...)
    cell.embeddings <- ica.results$M
  } else {
    ica.results <- ica.fxn(t(x = data.use), nc = ics.compute,...)
    cell.embeddings <- ica.results$S
  }
  gene.loadings <- (as.matrix(x = data.use ) %*% as.matrix(x = cell.embeddings))
  colnames(x = gene.loadings) <- paste0(reduction.key, 1:ncol(x = gene.loadings))
  colnames(x = cell.embeddings) <- paste0(reduction.key, 1:ncol(x = cell.embeddings))
  ica.obj <- new(
    Class = "dim.reduction",
    gene.loadings = gene.loadings,
    cell.embeddings = cell.embeddings,
    sdev = sqrt(x = ica.results$vafs),
    key = "IC"
  )

  eval(expr = parse(text = paste0("object@dr$", reduction.name, "<- ica.obj")))
  parameters.to.store <- as.list(environment(), all = TRUE)[names(formals("ICA"))]
  object <- SetCalcParams(object = object, calculation = "ICA", ... = parameters.to.store)
  if(is.null(object@calc.params$ICA$ic.genes)){
    object@calc.params$ICA$ic.genes <- rownames(data.use)
  }
  if(print.results){
    PrintDim(object = object, dims.print = ics.print, genes.print = genes.print,reduction.type = reduction.name)
  }
  return(object)
}

#' Run t-distributed Stochastic Neighbor Embedding
#'
#' Run t-SNE dimensionality reduction on selected features. Has the option of
#' running in a reduced dimensional space (i.e. spectral tSNE, recommended),
#' or running based on a set of genes. For details about stored TSNE calculation
#' parameters, see \code{PrintTSNEParams}.
#'
#' @param object Seurat object
#' @param reduction.use Which dimensional reduction (e.g. PCA, ICA) to use for
#' the tSNE. Default is PCA
#' @param cells.use Which cells to analyze (default, all cells)
#' @param dims.use Which dimensions to use as input features
#' @param genes.use If set, run the tSNE on this subset of genes
#' (instead of running on a set of reduced dimensions). Not set (NULL) by default
#' @param seed.use Random seed for the t-SNE
#' @param do.fast If TRUE, uses the Barnes-hut implementation, which runs
#' faster, but is less flexible. TRUE by default.
#' @param add.iter If an existing tSNE has already been computed, uses the
#' current tSNE to seed the algorithm and then adds additional iterations on top
#' of this
#' @param dim.embed The dimensional space of the resulting tSNE embedding
#' (default is 2). For example, set to 3 for a 3d tSNE
#' @param \dots Additional arguments to the tSNE call. Most commonly used is
#' perplexity (expected number of neighbors default is 30)
#' @param distance.matrix If set, runs tSNE on the given distance matrix
#' instead of data matrix (experimental)
#' @param reduction.name dimensional reduction name, specifies the position in the object$dr list. tsne by default
#' @param reduction.key dimensional reduction key, specifies the string before the number for the dimension names. tSNE_ by default
#'
#' @return Returns a Seurat object with a tSNE embedding in
#' object@@dr$tsne@cell.embeddings
#'
#' @importFrom Rtsne Rtsne
#' @importFrom tsne tsne
#'
#' @export
#'
#' @examples
#' pbmc_small
#' # Run tSNE on first five PCs, note that for test dataset (only 80 cells)
#' # we can't use default perplexity of 30
#' pbmc_small <- RunTSNE(pbmc_small, reduction.use = "pca", dims.use = 1:5, perplexity=10)
#' # Run tSNE on first five independent components from ICA
#' pbmc_small <- RunICA(pbmc_small,ics.compute=5)
#' pbmc_small <- RunTSNE(pbmc_small, reduction.use = "ica", dims.use = 1:5, perplexity=10)
#' # Plot results
#' TSNEPlot(pbmc_small)
#'
RunTSNE <- function(
  object,
  reduction.use = "pca",
  cells.use = NULL,
  dims.use = 1:5,
  genes.use = NULL,
  seed.use = 1,
  do.fast = TRUE,
  add.iter = 0,
  dim.embed = 2,
  distance.matrix = NULL,
  reduction.name = "tsne",
  reduction.key = "tSNE_",
  ...
) {
  if (! is.null(x = distance.matrix)) {
    genes.use <- rownames(x = object@data)
  }
  if (is.null(x = genes.use)) {
    data.use <- GetDimReduction(
      object = object,
      reduction.type = reduction.use,
      slot = "cell.embeddings"
    )[, dims.use]
  }
  if (! is.null(x = genes.use)) {
    data.use <- t(PrepDR(
      object = object,
      genes.use = genes.use))
  }
  set.seed(seed = seed.use)
  if (do.fast) {
    if (is.null(x = distance.matrix)) {
      data.tsne <- Rtsne(X = as.matrix(x = data.use), dims = dim.embed, ...)
    } else {
      data.tsne <- Rtsne(
        X = as.matrix(x = distance.matrix),
        dims = dim.embed,
        is_distance=TRUE
      )
    }
    data.tsne <- data.tsne$Y
  } else {
    data.tsne <- tsne(X = data.use, k = dim.embed, ...)
  }
  if (add.iter > 0) {
    data.tsne <- tsne(
      X = data.use,
      initial_config = as.matrix(x = data.tsne),
      max_iter = add.iter,
      ...
    )
  }
  colnames(x = data.tsne) <- paste0(reduction.key, 1:ncol(x = data.tsne))
  rownames(x = data.tsne) <- rownames(x = data.use)
  object <- SetDimReduction(
    object = object,
    reduction.type = reduction.name,
    slot = "cell.embeddings",
    new.data = data.tsne
  )
  object <- SetDimReduction(
    object = object,
    reduction.type = reduction.name,
    slot = "key",
    new.data = reduction.key
  )
  parameters.to.store <- as.list(environment(), all = TRUE)[names(formals("RunTSNE"))]
  object <- SetCalcParams(object = object, calculation = "RunTSNE", ... = parameters.to.store)
  if(!is.null(GetCalcParam(object = object, calculation = "RunTSNE", parameter = "genes.use"))){
    object@calc.params$RunTSNE$genes.use <- colnames(data.use)
    object@calc.params$RunTSNE$cells.use <- rownames(data.use)
  }
  return(object)
}

#' Project Dimensional reduction onto full dataset
#'
#' Takes a pre-computed dimensional reduction (typically calculated on a subset
#' of genes) and projects this onto the entire dataset (all genes). Note that
#' the cell loadings will remain unchanged, but now there are gene loadings for
#' all genes.
#'
#'
#' @param object Seurat object
#' @param reduction.type Reduction to use
#' @param dims.print Number of dims to print genes for
#' @param dims.store Number of dims to store (default is 30)
#' @param genes.print Number of genes with highest/lowest loadings to print for
#' each PC
#' @param replace.dim Replace the existing data (overwrite
#' object@@dr$XXX@gene.loadings), not done by default.
#' @param do.center Center the dataset prior to projection (should be set to TRUE)
#' @param do.print Print top genes associated with the projected dimensions
#' @param assay.type Data type, RNA by default. Can be changed for multimodal
#' datasets (i.e. project a PCA done on RNA, onto CITE-seq data)
#'
#' @return Returns Seurat object with the projected values in
#' object@@dr$XXX@gene.loadings.full
#'
#' @export
#'
#' @examples
#' pbmc_small
#' pbmc_small <- ProjectDim(pbmc_small, reduction.type = "pca")
#' # Vizualize top projected genes in heatmap
#' DimHeatmap(pbmc_small,pc.use = 1,use.full = TRUE,do.balanced = TRUE,reduction.type = "pca")
#'
ProjectDim <- function(
  object,
  reduction.type = "pca",
  dims.print = 1:5,
  dims.store = 30,
  genes.print = 30,
  replace.dim = FALSE,
  do.center = FALSE,
  do.print = TRUE,
  assay.type = "RNA"
) {
  if (! reduction.type %in% names(x = object@dr)) {
    stop(paste(reduction.type, "dimensional reduction has not been computed"))
  }
  data.use <- GetAssayData(
    object = object,
    assay.type = assay.type,
    slot = "scale.data"
  )
  if (do.center) {
    data.use <- scale(x = as.matrix(x = data.use), center = TRUE, scale = FALSE)
  }
  cell.embeddings <- GetDimReduction(
    object = object,
    reduction.type = reduction.type,
    slot = "cell.embeddings"
  )
  new.gene.loadings.full <- FastMatMult(m1 = data.use, m2 = cell.embeddings)
  rownames(x = new.gene.loadings.full) <- rownames(x = data.use)
  colnames(x = new.gene.loadings.full) <- colnames(x = cell.embeddings)
  object <- SetDimReduction(
    object = object,
    reduction.type = reduction.type,
    slot = "gene.loadings.full",
    new.data = new.gene.loadings.full
  )
  if (replace.dim) {
    object <- SetDimReduction(
      object = object,
      reduction.type = reduction.type,
      slot = "gene.loadings",
      new.data = new.gene.loadings.full
    )
  }
  if (do.print) {
    PrintDim(
      object = object,
      reduction.type = reduction.type,
      genes.print = genes.print,
      use.full = TRUE,
      dims.print = dims.print
    )
  }
  return(object)
}


#' Project Principal Components Analysis onto full dataset
#'
#' Takes a pre-computed PCA (typically calculated on a subset of genes) and
#' projects this onto the entire dataset (all genes). Note that the cell
#' loadings remains unchanged, but now there are gene loading scores for all
#' genes.
#'
#' @param object Seurat object
#' @param do.print Print top genes associated with the projected PCs
#' @param pcs.print Number of PCs to print genes for
#' @param pcs.store Number of PCs to store (default is 30)
#' @param genes.print Number of genes with highest/lowest loadings to print for
#' each PC
#' @param replace.pc Replace the existing PCA (overwite
#' object@@dr$pca@gene.loadings), not done by default.
#' @param do.center Center the dataset prior to projection (should be set to TRUE)
#'
#' @return Returns Seurat object with the projected PCA values in
#' object@@dr$pca@gene.loadings.full
#'
#' @export
#'
#' @examples
#' pbmc_small
#' pbmc_small <- ProjectPCA(pbmc_small)
#' # Vizualize top projected genes in heatmap
#' PCHeatmap(pbmc_small,pc.use = 1,use.full = TRUE,do.balanced = TRUE)
#'
ProjectPCA <- function(
  object,
  do.print = TRUE,
  pcs.print = 1:5,
  pcs.store = 30,
  genes.print = 30,
  replace.pc = FALSE,
  do.center = FALSE
) {
  return(ProjectDim(
    object,
    reduction.type = "pca",
    dims.print = pcs.print,
    genes.print = 30,
    replace.dim = replace.pc,
    do.center = do.center,
    do.print = do.print,
    dims.store = pcs.store
  ))
}

#' Perform Canonical Correlation Analysis
#'
#' Runs a canonical correlation analysis using a diagonal implementation of CCA.
#' For details about stored CCA calculation parameters, see
#' \code{PrintCCAParams}.
#'
#' @param object Seurat object
#' @param object2 Optional second object. If object2 is passed, object1 will be
#' considered as group1 and object2 as group2.
#' @param group1 First set of cells (or IDs) for CCA
#' @param group2 Second set of cells (or IDs) for CCA
#' @param group.by Factor to group by (column vector stored in object@@meta.data)
#' @param num.cc Number of canonical vectors to calculate
#' @param genes.use Set of genes to use in CCA. Default is object@@var.genes. If
#' two objects are given, the default is the union of both variable gene sets
#' that are also present in both objects.
#' @param scale.data Use the scaled data from the object
#' @param rescale.groups Rescale each set of cells independently
#' @param ... Extra parameters (passed onto MergeSeurat in case with two objects
#' passed, passed onto ScaleData in case with single object and rescale.groups
#' set to TRUE)
#'
#' @return Returns Seurat object with the CCA stored in the @@dr$cca slot. If
#' one object is passed, the same object is returned. If two are passed, a
#' combined object is returned.
#'
#' @seealso \code{MergeSeurat}
#'
#' @export
#'
#' @examples
#' pbmc_small
#' # As CCA requires two datasets, we will split our test object into two just for this example
#' pbmc1 <- SubsetData(pbmc_small,cells.use = pbmc_small@cell.names[1:40])
#' pbmc2 <- SubsetData(pbmc_small,cells.use = pbmc_small@cell.names[41:80])
#' pbmc1@meta.data$group <- "group1"
#' pbmc2@meta.data$group <- "group2"
#' pbmc_cca <- RunCCA(pbmc1,pbmc2)
#' # Print results
#' PrintDim(pbmc_cca,reduction.type = 'cca')
#'
RunCCA <- function(
  object,
  object2,
  group1,
  group2,
  group.by,
  num.cc = 20,
  genes.use,
  scale.data = TRUE,
  rescale.groups = FALSE,
  ...
) {
  if (! missing(x = object2) && (! missing(x = group1) || ! missing(x = group2))) {
    warning("Both object2 and group set. Continuing with objects defining the groups")
  }
  if (! missing(x = object2)) {
    if (missing(x = genes.use)) {
      genes.use <- union(x = object@var.genes, y = object2@var.genes)
      if (length(x = genes.use) == 0) {
        stop("No variable genes present. Run MeanVarPlot and retry")
      }
    }
    if (scale.data) {
      possible.genes <- intersect(
        x = rownames(x = object@scale.data),
        y = rownames(x = object2@scale.data)
      )
      genes.use <- genes.use[genes.use %in% possible.genes]
      data.use1 <- object@scale.data[genes.use, ]
      data.use2 <- object2@scale.data[genes.use, ]
    } else {
      possible.genes <- intersect(
        x = rownames(object@data),
        y = rownames(object2@data)
      )
      genes.use <- genes.use[genes.use %in% possible.genes]
      data.use1 <- object@data[genes.use, ]
      data.use2 <- object2@data[genes.use, ]
    }
    if (length(x = genes.use) == 0) {
      stop("0 valid genes in genes.use")
    }
  } else {
    if (missing(x = group1)) {
      stop("group1 not set")
    }
    if (missing(x = group2)) {
      stop("group2 not set")
    }
    if (! missing(x = group.by)) {
      if (! group.by %in% colnames(x = object@meta.data)) {
        stop("invalid group.by parameter")
      }
    }
    if (missing(x = genes.use)) {
      genes.use <- object@var.genes
      if (length(x = genes.use) == 0) {
        stop("No variable genes present. Run MeanVarPlot and retry")
      }
    }
    if (missing(x = group.by)) {
      cells.1 <- CheckGroup(object = object, group = group1, group.id = "group1")
      cells.2 <- CheckGroup(object = object, group = group2, group.id = "group2")
    } else {
      object.current.ids <- object@ident
      object <- SetAllIdent(object = object, id = group.by)
      cells.1 <- CheckGroup(object = object, group = group1, group.id = "group1")
      cells.2 <- CheckGroup(object = object, group = group2, group.id = "group2")
      object <- SetIdent(
        object = object,
        cells.use = object@cell.names,
        ident.use = object.current.ids
      )
    }
    if (scale.data) {
      if (rescale.groups) {
        data.use1 <- ScaleData(
          object = object,
          data.use = object@data[genes.use, cells.1],
          ...
        )
        data.use1 <- data.use1@scale.data
        data.use2 <- ScaleData(
          object = object,
          data.use = object@data[genes.use, cells.2],
          ...
        )
        data.use2 <- data.use2@scale.data
      } else {
        data.use1 <- object@scale.data[genes.use, cells.1]
        data.use2 <- object@scale.data[genes.use, cells.2]
      }
    } else {
      data.use1 <- object@data[genes.use, cells.1]
      data.use2 <- object@data[genes.use, cells.2]
    }
  }
  genes.use <- CheckGenes(data.use = data.use1, genes.use = genes.use)
  genes.use <- CheckGenes(data.use = data.use2, genes.use = genes.use)
  data.use1 <- data.use1[genes.use, ]
  data.use2 <- data.use2[genes.use, ]

  cat("Running CCA\n", file = stderr())

  cca.results <- CanonCor(
    mat1 = data.use1,
    mat2 = data.use2,
    standardize = TRUE,
    k = num.cc
  )
  cca.data <- rbind(cca.results$u, cca.results$v)
  colnames(x = cca.data) <- paste0("CC", 1:num.cc)
  rownames(cca.data) <- c(colnames(data.use1), colnames(data.use2))
  cca.data <- apply(cca.data, MARGIN = 2, function(x){
    if(sign(x[1]) == -1) {
      x <- x * -1
    }
    return(x)
  })
  # wipe old CCA slot
  object@dr$cca <- NULL
  if (! missing(x = object2)) {
    cat("Merging objects\n", file = stderr())
    combined.object <- MergeSeurat(
      object1 = object,
      object2 = object2,
      do.scale = FALSE,
      do.center = FALSE,
      ...
    )
    # to improve, to pull the same normalization and scale params as previously used
    combined.object <- ScaleData(object = combined.object)
    combined.object@scale.data[is.na(x = combined.object@scale.data)] <- 0
    combined.object@var.genes <- genes.use
    if("add.cell.id1" %in% names(list(...)) && "add.cell.id2" %in% names(list(...))) {
      o1.idx <- 1:length(object@cell.names)
      o2.idx <- (length(object@cell.names) + 1):(length(object@cell.names) + length(object2@cell.names))
      rownames(cca.data)[o1.idx] <-
        paste0(list(...)$add.cell.id1, "_", rownames(cca.data)[o1.idx])
      rownames(cca.data)[o2.idx] <-
        paste0(list(...)$add.cell.id2, "_", rownames(cca.data)[o2.idx])
    }
    combined.object <- SetDimReduction(
      object = combined.object,
      reduction.type = "cca",
      slot = "cell.embeddings",
      new.data = cca.data
    )
    combined.object <- SetDimReduction(
      object = combined.object,
      reduction.type = "cca",
      slot = "key",
      new.data = "CC"
    )
    combined.object <- ProjectDim(
      object = combined.object,
      reduction.type = "cca",
      do.print = FALSE
    )
    combined.object <- SetDimReduction(
      object = combined.object,
      reduction.type = "cca",
      slot = "gene.loadings",
      new.data = GetGeneLoadings(
        object = combined.object,
        reduction.type = "cca",
        use.full = TRUE,
        genes.use = genes.use
      )
    )
    parameters.to.store <- as.list(environment(), all = TRUE)[names(formals("RunCCA"))]
    combined.object <- SetCalcParams(
      object = combined.object,
      calculation = "RunCCA",
      ... = parameters.to.store
    )
    combined.object <- SetSingleCalcParam(
      object = combined.object,
      calculation = "RunCCA",
      parameter = "object.project",
      value = object@project.name
    )
    combined.object <- SetSingleCalcParam(
      object = combined.object,
      calculation = "RunCCA",
      parameter = "object2.project",
      value = object2@project.name
    )
    return(combined.object)
  } else {
    cca.data <- cca.data[object@cell.names, ]
    object <- SetDimReduction(
      object = object,
      reduction.type = "cca",
      slot = "cell.embeddings",
      new.data = cca.data
    )
    object <- SetDimReduction(
      object = object,
      reduction.type = "cca",
      slot = "key",
      new.data = "CC"
    )
    object <- ProjectDim(
      object = object,
      reduction.type = "cca",
      do.print = FALSE
    )
    object <- SetDimReduction(
      object = object,
      reduction.type = "cca",
      slot = "gene.loadings",
      new.data = GetGeneLoadings(
        object = object,
        reduction.type = "cca",
        use.full = TRUE,
        genes.use = genes.use
      )
    )
    object@scale.data[is.na(x = object@scale.data)] <- 0
    parameters.to.store <- as.list(environment(), all = TRUE)[names(formals("RunCCA"))]
    object <- SetCalcParams(
      object = object,
      calculation = "RunCCA",
      ... = parameters.to.store
    )
    return(object)
  }
}

#' Perform Canonical Correlation Analysis with more than two groups
#'
#' Runs a canonical correlation analysis
#'
#' @param object.list List of Seurat objects
#' @param genes.use Genes to use in mCCA.
#' @param niter Number of iterations to perform. Set by default to 25.
#' @param num.ccs Number of canonical vectors to calculate
#' @param standardize standardize scale.data matrices to be centered (mean zero)
#' and scaled to have a standard deviation of 1.
#' @return Returns a combined Seurat object with the CCA stored in the @@dr$cca slot.
#' @export
RunMultiCCA <- function(object.list, genes.use, niter = 25, num.ccs = 1, standardize = TRUE){
  set.seed(42)
  if(length(object.list) < 3){
    stop("Must give at least 3 objects/matrices for MultiCCA")
  }
  mat.list <- list()
  if(class(object.list[[1]]) == "seurat"){
    if (missing(x = genes.use)) {
      genes.use <- c()
      for(obj in object.list){
        genes.use <- c(genes.use, obj@var.genes)
      }
      genes.use <- unique(genes.use)
      if (length(x = genes.use) == 0) {
        stop("No variable genes present. Run MeanVarPlot and retry")
      }
    }
    for(obj in object.list) {
      genes.use <- CheckGenes(data.use = obj@scale.data, genes.use = genes.use)
    }
    for(i in 1:length(object.list)){
      mat.list[[i]] <- object.list[[i]]@scale.data[genes.use, ]
    }
  }
  else{
    stop("input data not Seurat objects")
  }
  num.sets <- length(mat.list)
  if(standardize){
    for (i in 1:num.sets){
      mat.list[[i]] <- Standardize(mat.list[[i]], display_progress = F)
    }
  }
  ws <- list()
  for (i in 1:num.sets){
    ws[[i]] <- irlba(mat.list[[i]], nv = num.ccs)$v[, 1:num.ccs, drop = F]
  }
  ws.init <- ws
  ws.final <- list()
  cors <- NULL
  for(i in 1:length(ws)){
    ws.final[[i]] <- matrix(0, nrow=ncol(mat.list[[i]]), ncol=num.ccs)
  }
  for (cc in 1:num.ccs){
    print(paste0("Computing CC ", cc))
    ws <- list()
    for (i in 1:length(ws.init)){
      ws[[i]] <- ws.init[[i]][, cc]
    }
    cur.iter <- 1
    crit.old <- -10
    crit <- -20
    storecrits <- NULL
    while(cur.iter <= niter && abs(crit.old - crit)/abs(crit.old) > 0.001 && crit.old !=0){
      crit.old <- crit
      crit <- GetCrit(mat.list, ws, num.sets)
      storecrits <- c(storecrits, crit)
      cur.iter <- cur.iter + 1
      for(i in 1:num.sets){
        ws[[i]] <- UpdateW(mat.list, i, num.sets, ws, ws.final)
      }
    }
    for(i in 1:length(ws)){
      ws.final[[i]][, cc] <- ws[[i]]
    }
    cors <- c(cors, GetCors(mat.list, ws, num.sets))
  }
  results <- list(ws=ws.final, ws.init=ws.init, num.sets = num.sets, cors=cors)
  combined.object <- object.list[[1]]
  for(i in 2:length(object.list)){
    combined.object <- MergeSeurat(object1 = combined.object, object2 = object.list[[i]], do.scale = F, do.center = F, do.normalize = F)
  }
  combined.object <- NormalizeData(combined.object)
  combined.object@meta.data$orig.ident <- sapply(combined.object@cell.names, ExtractField, 1)
  combined.object <- ScaleData(object = combined.object)
  combined.object@scale.data[is.na(x = combined.object@scale.data)] <- 0
  combined.object@var.genes <- genes.use
  cca.data <- results$ws[[1]]
  for(i in 2:length(object.list)){
    cca.data <- rbind(cca.data, results$ws[[i]])
  }
  rownames(cca.data) <- colnames(combined.object@data)
  cca.data <- apply(cca.data, MARGIN = 2, function(x){
    if(sign(x[1]) == -1) {
      x <- x * -1
    }
    return(x)
  })
  combined.object <- SetDimReduction(
    object = combined.object,
    reduction.type = "cca",
    slot = "cell.embeddings",
    new.data = cca.data
  )
  combined.object <- SetDimReduction(
    object = combined.object,
    reduction.type = "cca",
    slot = "key",
    new.data = "CC"
  )
  combined.object <- ProjectDim(
    object = combined.object,
    reduction.type = "cca",
    do.print = FALSE
  )
  combined.object <- SetDimReduction(
    object = combined.object,
    reduction.type = "cca",
    slot = "gene.loadings",
    new.data = GetGeneLoadings(
      object = combined.object,
      reduction.type = "cca",
      use.full = TRUE,
      genes.use = genes.use
    )
  )
  parameters.to.store <- as.list(environment(), all = TRUE)[names(formals("RunMultiCCA"))]
  parameters.to.store$object.list <- NULL
  combined.object <- SetCalcParams(object = combined.object,
                                   calculation = "RunMultiCCA",
                                   ... = parameters.to.store
  )
  return(combined.object)
}

#' Run diffusion map
#'
#' @param object Seurat object
#' @param cells.use Which cells to analyze (default, all cells)
#' @param dims.use Which dimensions to use as input features
#' @param genes.use If set, run the diffusion map procedure on this subset of
#' genes (instead of running on a set of reduced dimensions). Not set (NULL) by
#' default
#' @param reduction.use Which dimensional reduction (PCA or ICA) to use for the
#' diffusion map input. Default is PCA
#' @param q.use Quantile to clip diffusion map components at. This addresses an issue where 1-2 cells will have extreme values that obscure all other points. 0.01 by default
#' @param max.dim Max dimension to keep from diffusion calculation
#' @param scale.clip Max/min value for scaled data. Default is 3
#' @param reduction.name dimensional reduction name, specifies the position in the object$dr list. dm by default
#' @param reduction.key dimensional reduction key, specifies the string before the number for the dimension names. DM by default
#' @param ... Additional arguments to the diffuse call
#'
#' @return Returns a Seurat object with a diffusion map
#'
#' @import diffusionMap
#' @importFrom stats dist quantile
#'
#' @export
#'
#' @examples
#' pbmc_small
#' # Run Diffusion on variable genes
#' pbmc_small <- RunDiffusion(pbmc_small,genes.use = pbmc_small@var.genes)
#' # Run Diffusion map on first 10 PCs
#' pbmc_small <- RunDiffusion(pbmc_small,genes.use = pbmc_small@var.genes)
#' # Plot results
#' DMPlot(pbmc_small)
#'
RunDiffusion <- function(
  object,
  cells.use = NULL,
  dims.use = 1:5,
  genes.use = NULL,
  reduction.use = 'pca',
  q.use = 0.01,
  max.dim = 2,
  scale.clip = 10,
  reduction.name = "dm",
  reduction.key = "DM",
  ...
) {
  cells.use <- SetIfNull(x = cells.use, default = colnames(x = object@data))
  if (is.null(x = genes.use)) {
    dim.code <- GetDimReduction(
      object = object,
      reduction.type = reduction.use,
      slot = 'key'
    )
    dim.codes <- paste0(dim.code, dims.use)
    data.use <- FetchData(object = object, vars.all = dim.codes)
  }
  if (! is.null(x = genes.use)) {
    genes.use <- intersect(x = genes.use, y = rownames(x = object@scale.data))
    data.use <- MinMax(
      data = t(x = object@data[genes.use, cells.use]),
      min = -1 * scale.clip,
      max = scale.clip
    )
  }
  parameters.to.store <- as.list(environment(), all = TRUE)[names(formals("RunDiffusion"))]
  object <- SetCalcParams(object = object,
                          calculation = "RunDiffusion",
                          ... = parameters.to.store)
  data.dist <- dist(data.use)
  data.diffusion <- data.frame(
    diffuse(
      D = data.dist,
      neigen = max.dim,
      maxdim = max.dim,
      ...
    )$X
  )
  colnames(x = data.diffusion) <- paste0(reduction.key, 1:ncol(x = data.diffusion))
  rownames(x = data.diffusion) <- cells.use
  for (i in 1:max.dim) {
    x <- data.diffusion[,i]
    x <- MinMax(
      data = x,
      min = quantile(x = x, probs = q.use),
      quantile(x = x, probs = 1-q.use)
    )
    data.diffusion[, i] <- x
  }
  object <- SetDimReduction(
    object = object,
    reduction.type = reduction.name,
    slot = "cell.embeddings",
    new.data = as.matrix(x = data.diffusion)
  )
  object <- SetDimReduction(
    object = object,
    reduction.type = reduction.name,
    slot = "key",
    new.data = "DM"
  )
  return(object)
}
