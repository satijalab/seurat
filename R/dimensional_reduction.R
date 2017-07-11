#' Run Principal Component Analysis on gene expression using IRLBA
#'
#' Run a PCA dimensionality reduction. For details about stored PCA calculation
#' parameters, see \code{\link{PrintPCAParams}}.
#'
#' @param object Seurat object
#' @param pc.genes Genes to use as input for PCA. Default is object@@var.genes
#' @param pcs.compute Total Number of PCs to compute and store
#' @param use.imputed Run PCA on imputed values (FALSE by default)
#' @param rev.pca By default computes the PCA on the cell x gene matrix. Setting
#' to true will compute it on gene x cell matrix.
#' @param scale.by.varexp Scale the gene loadings by the variance explained by
#' each PC (scales the cell embeddings if rev.pca is TRUE)
#' @param do.print Print the top genes associated with high/low loadings for
#' the PCs
#' @param pcs.print PCs to print genes for
#' @param genes.print Number of genes to print for each PC
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
PCA <- function(
  object,
  pc.genes = NULL,
  pcs.compute = 20,
  use.imputed = FALSE,
  rev.pca = FALSE,
  scale.by.varexp = TRUE,
  do.print = TRUE,
  pcs.print = 1:5,
  genes.print = 30,
  ...
) {
  data.use <- PrepDR(
    object = object,
    genes.use = pc.genes,
    use.imputed = use.imputed)
  pcs.compute <- min(pcs.compute, ncol(x = data.use))
  if (rev.pca) {
    pca.results <- irlba(A = data.use, nv = pcs.compute, ...)
    if(scale.by.varexp){
      gene.loadings <- pca.results$u %*% diag(pca.results$d)
    } else{
      gene.loadings <- pca.results$u
    }
    cell.embeddings <- pca.results$v
  }
  else {
    pca.results <- irlba(A = t(x = data.use), nv = pcs.compute, ...)
    gene.loadings <- pca.results$v
    if(scale.by.varexp){
      cell.embeddings <- pca.results$u %*% diag(pca.results$d)
    } else {
      cell.embeddings <- pca.results$u
    }
  }
  rownames(x = gene.loadings) <- rownames(x = data.use)
  colnames(x = gene.loadings) <- paste0("PC", 1:pcs.compute)
  rownames(x = cell.embeddings) <- colnames(x = data.use)
  colnames(x = cell.embeddings) <- colnames(x = gene.loadings)
  pca.obj <- new(
    Class = "dim.reduction",
    gene.loadings = gene.loadings,
    cell.embeddings = cell.embeddings,
    sdev = pca.results$d,
    key = "PC"
  )
  object@dr$pca <- pca.obj
  parameters.to.store <- as.list(environment(), all = TRUE)[names(formals("PCA"))]
  object <- SetCalcParams(object = object, calculation = "PCA", ... = parameters.to.store)
  if(is.null(object@calc.params$PCA$pc.genes)){
    object@calc.params$PCA$pc.genes <- rownames(data.use)
  }
  return(object)
}


#' Run Independent Component Analysis on gene expression
#'
#' Run fastica algorithm from the ica package for ICA dimensionality reduction.
#' For details about stored ICA calculation parameters, see
#' \code{\link{PrintICAParams}}.
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
#' @param \dots Additional arguments to be passed to fastica
#'
#' @importFrom ica icafast icaimax icajade
#'
#' @return Returns Seurat object with an ICA calculation stored in
#' object@@dr$ica
#'
#' @export
#'
ICA <- function(
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
  colnames(x = gene.loadings) <- paste0("IC", 1:ncol(x = gene.loadings))
  colnames(x = cell.embeddings) <- paste0("IC", 1:ncol(x = cell.embeddings))
  ica.obj <- new(
    Class = "dim.reduction",
    gene.loadings = gene.loadings,
    cell.embeddings = cell.embeddings,
    sdev = sqrt(x = ica.results$vafs),
    key = "IC"
  )
  object@dr$ica <- ica.obj
  parameters.to.store <- as.list(environment(), all = TRUE)[names(formals("ICA"))]
  object <- SetCalcParams(object = object, calculation = "ICA", ... = parameters.to.store)
  if(is.null(object@calc.params$ICA$ic.genes)){
    object@calc.params$ICA$ic.genes <- rownames(data.use)
  }
  return(object)
}


#' Run t-distributed Stochastic Neighbor Embedding
#'
#' Run t-SNE dimensionality reduction on selected features. Has the option of
#' running in a reduced dimensional space (i.e. spectral tSNE, recommended),
#' or running based on a set of genes. For details about stored TSNE calculation
#' parameters, see \code{\link{PrintTSNEParams}}.
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
#' faster, but is less flexible
#' @param add.iter If an existing tSNE has already been computed, uses the
#' current tSNE to seed the algorithm and then adds additional iterations on top
#' of this
#' @param dim.embed The dimensional space of the resulting tSNE embedding
#' (default is 2). For example, set to 3 for a 3d tSNE
#' @param \dots Additional arguments to the tSNE call. Most commonly used is
#' perplexity (expected number of neighbors default is 30)
#' @param distance.matrix If set, tuns tSNE on the given distance matrix
#' instead of data matrix (experimental)
#'
#' @return Returns a Seurat object with a tSNE embedding in
#' object@@dr$tsne@cell.embeddings
#'
#' @importFrom Rtsne Rtsne
#' @importFrom tsne tsne
#'
#' @export
#'
RunTSNE <- function(
  object,
  reduction.use = "pca",
  cells.use = NULL,
  dims.use = 1:5,
  genes.use = NULL,
  seed.use = 1,
  do.fast = FALSE,
  add.iter = 0,
  dim.embed = 2,
  distance.matrix = NULL,
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
    if (length(x = object@scale.data) == 0) {
      stop("Object@scale.data has not been set. Run ScaleData() and then retry.")
    }
    cells.use <- SetIfNull(x = cells.use, default = colnames(x = object@scale.data))
    genes.use <- ainb(a = genes.use, b = rownames(x = object@scale.data))
    data.use <- t(x = object@scale.data[genes.use, cells.use])
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
      x = data.use,
      initial_config = as.matrix(x = data.tsne),
      max_iter = add.iter,
      ...
    )
  }
  colnames(x = data.tsne) <- paste0("tSNE_", 1:ncol(x = data.tsne))
  rownames(x = data.tsne) <- rownames(x = data.use)
  object <- SetDimReduction(
    object = object,
    reduction.type = "tsne",
    slot = "cell.embeddings",
    new.data = data.tsne
  )
  object <- SetDimReduction(
    object = object,
    reduction.type = "tsne",
    slot = "key",
    new.data = "tSNE_"
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
#' \code{\link{PrintCCAParams}}.
#'
#' @param object Seurat object
#' @param object2 Optional second object. If object2 is passed, object1 will be
#' considered as group1 and object2 as group2.
#' @param group1 First set of cells (or IDs) for CCA
#' @param group2 Second set of cells (or IDs) for CCA
#' @param group.by Factor to group by (column vector stored in object@@data.info)
#' @param num.cc Number of canonical vectors to calculate
#' @param genes.use Set of genes to use in CCA. Default is object@@var.genes. If
#' two objects are given, the default is the union of both variable gene sets
#' that are also present in both objects.
#' @param scale.data Use the scaled data from the object
#' @param rescale.groups Rescale each set of cells independently
#' @return Returns Seurat object with the CCA stored in the @@dr$cca slot. If
#' one object is passed, the same object is returned. If two are passed, a
#' combined object is returned.
#' @export
RunCCA <- function(
  object,
  object2,
  group1,
  group2,
  group.by,
  num.cc = 20,
  genes.use,
  scale.data = TRUE,
  rescale.groups = FALSE
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
      if (! group.by %in% colnames(x = object@data.info)) {
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
          data.use = object@data[genes.use, cells.1]
        )
        data.use1 <- data.use1@scale.data
        data.use2 <- ScaleData(
          object = object,
          data.use = object@data[genes.use, cells.2]
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
  cca.results <- CanonCor(
      mat1 = data.use1,
      mat2 = data.use2,
      standardize = TRUE,
      k = num.cc
  )
  cca.data <- rbind(cca.results$u, cca.results$v)
  rownames(x = cca.data) <- c(colnames(x = data.use1), colnames(x = data.use2))
  colnames(x = cca.data) <- paste0("CC", 1:num.cc)
  if (! missing(x = object2)) {
    cat("Merging objects\n", file = stderr())
    combined.object <- MergeSeurat(
      object1 = object,
      object2 = object2,
      do.scale = FALSE,
      do.center = FALSE
    )
    combined.object@scale.data[which(x = is.na(x = combined.object@scale.data))] <- 0
    combined.object@var.genes <- genes.use
    combined.object <- FastScaleData(object = combined.object)
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
    combined.object <- SetCalcParams(object = combined.object,
                                     calculation = "RunCCA",
                                     ... = parameters.to.store)
    combined.object <- SetSingleCalcParam(object = combined.object,
                                          calculation = "RunCCA",
                                          parameter = "object.project",
                                          value = object@project.name)
    combined.object <- SetSingleCalcParam(object = combined.object,
                                          calculation = "RunCCA",
                                          parameter = "object2.project",
                                          value = object2@project.name)
    return(combined.object)
  } else {
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

    object <- ProjectDim(object = object,
                         reduction.type = "cca",
                         do.print = FALSE)
    object@scale.data[is.na(x = object@scale.data)] <- 0
    parameters.to.store <- as.list(environment(), all = TRUE)[names(formals("RunCCA"))]
    object <- SetCalcParams(object = object,
                                     calculation = "RunCCA",
                                     ... = parameters.to.store)
    return(object)
  }
}

#' Calculate the ratio of variance explained by ICA or PCA to CCA
#'
#' @param object Seurat object
#' @param reduction.type type of dimensional reduction to compare to CCA (pca,
#' pcafast, ica)
#' @param grouping.var variable to group by
#' @param dims.use Vector of dimensions to project onto (default is the 1:number
#'  stored for cca)
#'
#' @return Returns Seurat object with ratio of variance explained stored in
#' object@@data.info$var.ratio
#' @export
#'
CalcVarExpRatio <- function(
  object,
  reduction.type = "pca",
  grouping.var,
  dims.use
) {
  if (missing(x = grouping.var)) {
    stop("Need to provide grouping variable")
  }
  if (missing(x = dims.use)) {
    dims.use <- 1:ncol(x = GetCellEmbeddings(object = object, reduction.type = "cca"))
  }
  groups <- as.vector(x = unique(x = FetchData(
    object = object,
    vars.all = grouping.var
  )[, 1]))
  genes.use <- rownames(x = GetGeneLoadings(object = object, reduction.type = "cca"))
  var.ratio <- data.frame()
  for (group in groups) {
    cat(paste("Calculating for", group, "\n"), file = stderr())
    group.cells <- WhichCells(
      object = object,
      subset.name = grouping.var,
      accept.value = group
    )
    cat(paste("\t Separating", group, "cells\n"), file = stderr())
    group.object <- SubsetData(object = object, cells.use = group.cells)
    cat("\t Running Dimensional Reduction \n", file = stderr())
    ldp.cca <- CalcLDProj(
      object = group.object,
      reduction.type = "cca",
      dims.use = dims.use,
      genes.use = genes.use
    )
    group.object <- CalcProjectedVar(
      object = group.object,
      low.dim.data = ldp.cca,
      reduction.type = "cca",
      dims.use = dims.use,
      genes.use = genes.use
    )
    if (reduction.type == "pca") {
      group.object <- PCA(
        object = group.object,
        pc.genes = genes.use,
        do.print = FALSE
      )
      ldp.pca <- CalcLDProj(
        object = group.object,
        reduction.type = "pca",
        dims.use = dims.use,
        genes.use = genes.use
      )
      group.object <- CalcProjectedVar(
        object = group.object,
        low.dim.data = ldp.pca,
        reduction.type = "pca",
        dims.use = dims.use,
        genes.use = genes.use
      )
      group.var.ratio <- group.object@data.info[, "cca.var", drop = FALSE] /
        group.object@data.info[, "pca.var", drop = FALSE]
    } else if (reduction.type == "ica") {
      group.object <- ICA(
        object = group.object,
        ic.genes = genes.use,
        print.results = FALSE
      )
      ldp.ica <- CalcLDProj(
        object = group.object,
        reduction.type = "ica",
        dims.use = dims.use,
        genes.use = genes.use
      )
      group.object <- CalcProjectedVar(
        object = group.object,
        low.dim.data = ldp.ica,
        reduction.type = "ica",
        dims.use = dims.use,
        genes.use = genes.use
      )
      group.var.ratio <- group.object@data.info[, "cca.var", drop = FALSE] /
        group.object@data.info[, "ica.var", drop = FALSE]
    } else if (reduction.type == "pcafast") {
      group.object <- PCAFast(
        object = group.object,
        ic.genes = genes.use,
        do.print = FALSE
      )
      ldp.pca <- CalcLDProj(
        object = group.object,
        reduction.type = "pca",
        dims.use = dims.use,
        genes.use = genes.use
      )
      group.object <- CalcProjectedVar(
        object = group.object,
        low.dim.data = ldp.pca,
        reduction.type = "pca",
        dims.use = dims.use,
        genes.use = genes.use
      )
      group.var.ratio <- group.object@data.info[, "cca.var", drop = FALSE] /
        group.object@data.info[, "pca.var", drop = FALSE]
    } else {
      stop(paste("reduction.type", reduction.type, "not supported"))
    }
    var.ratio <- rbind(var.ratio, group.var.ratio)
  }
  var.ratio$cell.name <- rownames(x = var.ratio)
  eval(expr = parse(text = paste0(
    "object@data.info$var.ratio.",
    reduction.type,
    "<- NULL"
  )))
  colnames(x = var.ratio) <- c(
    paste0("var.ratio.", reduction.type),
    "cell.name"
  )
  object@data.info$cell.name <- rownames(x = object@data.info)
  object@data.info <- merge(x = object@data.info, y = var.ratio, by = "cell.name")
  rownames(x = object@data.info) <- object@data.info$cell.name
  object@data.info$cell.name <- NULL
  return(object)
}

#' Align subspaces using dynamic time warping (DTW)
#'
#' Aligns subspaces so that they line up across grouping variable (only
#' implemented for case with 2 categories in grouping.var)
#'
#'
#' @param object Seurat object
#' @param reduction.type reduction to align scores for
#' @param group.var Name of the grouping variable for which to align the scores
#' @param dims.align Dims to align, default is all
#' @param num.genes Number of genes to use in construction of "metagene"
#' @param show.plots show debugging plots
#'
#' @return Returns Seurat object with the dims aligned, stored in
#'  object@@dr$reduction.type.aligned
#'
#' @importFrom dtw dtw
#' @importFrom WGCNA bicor
#' @importFrom pbapply pbapply
#'
#' @export
#'
AlignSubspace <- function(
  object,
  reduction.type,
  grouping.var,
  dims.align,
  num.genes = 30,
  show.plots = FALSE
) {
  ident.orig <- object@ident
  object <- SetAllIdent(object = object, id = grouping.var)
  levels.split <- names(x = sort(x = table(object@ident)))
  if (length(x = levels.split) != 2) {
    stop(paste0(
      "There are not two options for ",
      grouping.var,
      ". \n Current groups include: ",
      paste(levels.split, collapse = ", ")
    ))
  }
  objects <- list(
    SubsetData(object = object, ident.use = levels.split[1]),
    SubsetData(object = object, ident.use = levels.split[2])
  )
  object@ident <- ident.orig
  cc.loadings <- list()
  scaled.data <- list()
  cc.embeds <- list()
  for (i in 1:2) {
    cat(paste0("Rescaling group ", i, "\n"), file = stderr())
    objects[[i]] <- FastScaleData(object = objects[[i]])
    objects[[i]]@scale.data[is.na(x = objects[[i]]@scale.data)] <- 0
    objects[[i]] <- ProjectDim(
      object = objects[[i]],
      reduction.type = reduction.type,
      do.print = FALSE
    )
    cc.loadings[[i]] <- GetGeneLoadings(
      object = objects[[i]],
      reduction.type = reduction.type,
      use.full = TRUE
    )
    cc.embeds[[i]] <- GetCellEmbeddings(
      object = objects[[i]],
      reduction.type = reduction.type
    )
    scaled.data[[i]] <- objects[[i]]@scale.data
  }
  cc.embeds.both <- GetCellEmbeddings(object = object, reduction.type = reduction.type)
  colnames(cc.embeds.both) <- paste0("A", colnames(x = cc.embeds.both))
  cc.embeds.orig <- cc.embeds.both
  for (cc.use in dims.align) {
    cat(paste0("Aligning dimension ", cc.use, "\n"), file = stderr())
    genes.rank <- data.frame(
      rank(x = abs(x = cc.loadings[[1]][, cc.use])),
      rank(x = abs(x = cc.loadings[[2]][, cc.use])),
      cc.loadings[[1]][, cc.use],
      cc.loadings[[2]][, cc.use]
    )
    genes.rank$min <- apply(X = genes.rank[,1:2], MARGIN = 1, FUN = min)
    genes.rank <- genes.rank[order(genes.rank$min, decreasing = TRUE), ]
    genes.top <- rownames(x = genes.rank)[1:200]
    bicors <- list()
    for (i in 1:2) {
      cc.vals <- cc.embeds[[i]][, cc.use]
      bicors[[i]] <- pbsapply(
        X = genes.top,
        FUN = function(x) {
          return(suppressWarnings(expr = bicor(x = cc.vals, scaled.data[[i]][x, ])))
        }
      )
    }
    genes.rank <- data.frame(
      rank(x = abs(x = bicors[[1]])),
      rank(x = abs(x = bicors[[2]])),
      bicors[[1]],
      bicors[[2]]
    )
    genes.rank$min <- apply(X = abs(x = genes.rank[, 1:2]), MARGIN = 1, FUN = min)
    genes.rank <- genes.rank[order(genes.rank$min, decreasing = TRUE), ]
    genes.use <- rownames(x = genes.rank)[1:num.genes]
    metagenes <- list()
    multvar.data <- list()
    for (i in 1:2) {
      scaled.use <- sweep(
        x = scaled.data[[i]][genes.use, ],
        MARGIN = 1,
        STATS = sign(x = genes.rank[genes.use, i + 2]),
        FUN = "*"
      )
      scaled.use <- scaled.use[, names(x = sort(x = cc.embeds[[i]][, cc.use]))]
      metagenes[[i]] <- apply(
        X = scaled.use[genes.use, ],
        MARGIN = 2,
        FUN = mean,
        remove.na = TRUE
      )
      metagenes[[i]] <- (
        cc.loadings[[i]][genes.use, cc.use] %*% scaled.data[[i]][genes.use, ]
      )[1, colnames(x = scaled.use)]
    }

    mean.difference <- mean(x = ReferenceRange(x = metagenes[[1]])) -
      mean(x = ReferenceRange(x = metagenes[[2]]))
    metric.use <- "Euclidean"
    align.1 <- ReferenceRange(x = metagenes[[1]])
    align.2 <- ReferenceRange(x = metagenes[[2]])
    a1q <- sapply(
      X = seq(from = 0, to = 1, by = 0.001),
      FUN = function(x) {
        return(quantile(x = align.1, probs = x))
      }
    )
    a2q <- sapply(
      X = seq(from = 0, to = 1, by = 0.001),
      FUN = function(x) {
        quantile(x = align.2, probs = x)
      }
    )
    iqr <- (a1q - a2q)[100:900]
    iqr.x <- which.min(x = abs(x = iqr))
    iqrmin <- iqr[iqr.x]
    if (show.plots) {
      print(iqrmin)
    }
    align.2 <- align.2 + iqrmin
    alignment <- dtw(
      x = align.1,
      y = align.2,
      keep = TRUE,
      dist.method = metric.use
    )
    alignment.map <- data.frame(alignment$index1, alignment$index2)
    alignment.map$cc_data1 <- sort(cc.embeds[[1]][, cc.use])[alignment$index1]
    alignment.map$cc_data2 <- sort(cc.embeds[[2]][, cc.use])[alignment$index2]
    alignment.map.orig <- alignment.map
    alignment.map <- alignment.map[! duplicated(x = alignment.map$alignment.index1), ]
    cc.embeds.both[names(x = sort(x = cc.embeds[[1]][, cc.use])), cc.use] <- alignment.map$cc_data2
    if (show.plots) {
      par(mfrow = c(3, 2))
      plot(x = ReferenceRange(x = metagenes[[1]]), main = cc.use)
      plot(x = ReferenceRange(x = metagenes[[2]]))
      plot(
        x = ReferenceRange(x = metagenes[[1]])[(alignment.map.orig$alignment.index1)],
        pch = 16
      )
      points(
        x = ReferenceRange(metagenes[[2]])[(alignment.map.orig$alignment.index2)],
        col = "red",
        pch = 16,
        cex = 0.4
      )
      plot(x = density(x = alignment.map$cc_data2))
      lines(x = density(x = sort(x = cc.embeds[[2]][, cc.use])), col = "red")
      plot(x = alignment.map.orig$cc_data1)
      points(x = alignment.map.orig$cc_data2, col = "red")
    }
  }
  new.type <- paste0(reduction.type, ".aligned")
  new.key <- paste0(
    "A",
    GetDimReduction(
      object = object,
      reduction.type = reduction.type,
      slot = "key"
    )
  )
  object <- SetDimReduction(
    object = object,
    reduction.type = new.type,
    slot = "cell.embeddings",
    new.data = scale(x = cc.embeds.both)
  )
  object <- SetDimReduction(
    object = object,
    reduction.type = new.type,
    slot = "key",
    new.data = new.key
  )
  return(object)
}

#' Run diffusion map
#'
#' @param object Seurat object
#' @param cells.use Which cells to analyze (default, all cells)
#' @param dims.use Which dimensions to use as input features
#' @param genes.use If set, run the tSNE on this subset of genes
#' (instead of running on a set of reduced dimensions). Not set (NULL) by default
#' @param reduction.use Which dimensional reduction (PCA or ICA) to use for the tSNE. Default is PCA
#' @param q.use Quantile to use
#' @param max.dim Max dimension to keep from diffusion calculation
#' @param scale.clip Max/min value for scaled data. Default is 3
#' @param ... Additional arguments to the diffuse call
#'
#' @return Returns a Seurat object with a diffusion map
#'
#' @import diffusionMap
#'
#' @export
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
    genes.use <- ainb(genes.use, rownames(x = object@scale.data))
    data.use <- minmax(
      data = t(x = object@data[genes.use, cells.use]),
      min = -1 * scale.clip,
      max = scale.clip
    )
  }
  data.dist <- dist(data.use)
  data.diffusion <- data.frame(
    diffuse( # Where is diffuse?
      D = data.dist,
      neigen = max.dim,
      maxdim = max.dim,
      ...
    )$X
  )
  colnames(x = data.diffusion) <- paste0("DM", 1:ncol(x = data.diffusion))
  rownames(x = data.diffusion) <- cells.use
  for (i in 1:max.dim) {
    x <- data.diffusion[,i]
    x <- minmax(
      data = x,
      min = quantile(x = x, probs = q.use),
      quantile(x = x, probs = 1-q.use)
    )
    data.diffusion[, i] <- x
  }
  object <- SetDimReduction(
    object = object,
    reduction.type = "dm",
    slot = "cell.embeddings",
    new.data = as.matrix(x = data.diffusion)
  )
  object <- SetDimReduction(
    object = object,
    reduction.type = "dm",
    slot = "key",
    new.data = "DM"
  )
  return(object)
}
