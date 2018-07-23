#' @export
#'
RunPCA.default <- function(
  object,
  assay.use = NULL,
  features.use = NULL,
  compute.dims = 20,
  rev.pca = FALSE,
  weight.by.var = TRUE,
  verbose = TRUE,
  print.dims = 1:5,
  features.print = 30,
  reduction.name = "pca",
  reduction.key = "PC",
  seed.use = 42,
  ...
) {
  if (!is.null(x = seed.use)) {
    set.seed(seed = seed.use)
  }
  if (rev.pca) {
    compute.dims <- min(compute.dims, ncol(x = object) - 1)
    pca.results <- irlba(A = object, nv = compute.dims, ...)
    sdev <- pca.results$d/sqrt(max(1, nrow(x = object) - 1))
    if (weight.by.var) {
      feature.loadings <- pca.results$u %*% diag(pca.results$d)
    } else{
      feature.loadings <- pca.results$u
    }
    cell.embeddings <- pca.results$v
  }
  else {
    compute.dims <- min(compute.dims, nrow(x = object) - 1)
    pca.results <- irlba(A = t(x = object), nv = compute.dims, ...)
    feature.loadings <- pca.results$v
    sdev <- pca.results$d/sqrt(max(1, ncol(object) - 1))
    if (weight.by.var) {
      cell.embeddings <- pca.results$u %*% diag(pca.results$d)
    } else {
      cell.embeddings <- pca.results$u
    }
  }
  rownames(x = feature.loadings) <- rownames(x = object)
  colnames(x = feature.loadings) <- paste0(reduction.key, 1:compute.dims)
  rownames(x = cell.embeddings) <- colnames(x = object)
  colnames(x = cell.embeddings) <- colnames(x = feature.loadings)
  reduction.data <- MakeDimReducObject(
    cell.embeddings = cell.embeddings,
    feature.loadings = feature.loadings,
    assay.used = assay.use,
    stdev = sdev,
    key = reduction.key
  )
  return(reduction.data)
}

#' @describeIn RunPCA Run a PCA on an Assay object
#' @export
#' @method RunPCA Assay
#'
RunPCA.Assay <- function(
  object,
  assay.use = NULL,
  features.use = NULL,
  compute.dims = 20,
  rev.pca = FALSE,
  weight.by.var = TRUE,
  verbose = TRUE,
  print.dims = 1:5,
  features.print = 30,
  reduction.name = "pca",
  reduction.key = "PC",
  seed.use = 42,
  ...
) {
  data.use <- PrepDR(
    object = object,
    features.use = features.use
  )
  reduction.data <- RunPCA(
    object = data.use,
    assay.use = assay.use,
    pc.features = features.use,
    compute.dims = compute.dims,
    rev.pca = rev.pca,
    weight.by.var = weight.by.var,
    verbose = verbose,
    print.dims = print.dims,
    features.print = features.print,
    reduction.name = reduction.name,
    reduction.key = reduction.key,
    seed.use = seed.use,
    ...

  )
  return(reduction.data)
}

#' @param workflow.name Name of workflow
#'
#' @describeIn RunPCA Run a PCA on a Seurat object
#' @export
#' @method RunPCA Seurat
#'
RunPCA.Seurat <- function(
  object,
  assay.use = NULL,
  features.use = NULL,
  compute.dims = 20,
  rev.pca = FALSE,
  weight.by.var = TRUE,
  verbose = TRUE,
  print.dims = 1:5,
  features.print = 30,
  reduction.name = "pca",
  reduction.key = "PC",
  seed.use = 42,
  workflow.name = NULL,
  ...
) {
  if (!is.null(workflow.name)) {
    object <- PrepareWorkflow(object = object, workflow.name = workflow.name)
  }
  assay.use <- assay.use %||% DefaultAssay(object = object)
  assay.data <- GetAssay(object = object, assay.use = assay.use)
  reduction.data <- RunPCA(
    object = assay.data,
    assay.use = assay.use,
    features.use = features.use,
    compute.dims = compute.dims,
    rev.pca = rev.pca,
    weight.by.var = weight.by.var,
    verbose = verbose,
    print.dims = print.dims,
    feature.print = features.print,
    reduction.name = reduction.name,
    reduction.key = reduction.key,
    seed.use = seed.use,
    ...
  )
  object[[reduction.name]] <- reduction.data
  object <- LogSeuratCommand(object = object)
  if (!is.null(workflow.name)) {
    object <- UpdateWorkflow(object = object, workflow.name = workflow.name)
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
#' @importFrom methods new
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
#' @param tsne.method Select the method to use to compute the tSNE. Available
#' methods are:
#' \itemize{
#' \item{Rtsne: }{Use the Rtsne package Barnes-Hut implementation of tSNE (default)}
#' \item{tsne: }{standard tsne - not recommended for large datasets}
#' \item{FIt-SNE: }{Use the FFT-accelerated Interpolation-based t-SNE. Based on
#' Kluger Lab code found here: https://github.com/KlugerLab/FIt-SNE}
#' }
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
#' \dontrun{
#' pbmc_small
#' # Run tSNE on first five PCs, note that for test dataset (only 80 cells)
#' # we can't use default perplexity of 30
#' pbmc_small <- RunTSNE(pbmc_small, reduction.use = "pca", dims.use = 1:5, perplexity=10)
#' # Run tSNE on first five independent components from ICA
#' pbmc_small <- RunICA(pbmc_small,ics.compute=5)
#' pbmc_small <- RunTSNE(pbmc_small, reduction.use = "ica", dims.use = 1:5, perplexity=10)
#' # Plot results
#' TSNEPlot(pbmc_small)
#' }
#'
RunTSNEOld <- function(
  object,
  reduction.use = "pca",
  # cells.use = NULL,
  dims.use = 1:5,
  features.use = NULL,
  seed.use = 1,
  tsne.method = "Rtsne",
  add.iter = 0,
  dim.embed = 2,
  distance.matrix = NULL,
  reduction.name = "tsne",
  reduction.key = "tSNE_",
  ...
) {
  if (length(x = dims.use) < 2) {
    stop("Cannot perform tSNE on only one dimension, please provide two or more dimensions")
  }
  if (!is.null(x = distance.matrix)) {
    genes.use <- rownames(x = object@data)
  }
  if (is.null(x = genes.use)) {
    data.use <- GetDimReduction(
      object = object,
      reduction.type = reduction.use,
      slot = "cell.embeddings"
    )[, dims.use]
  }
  if (!is.null(x = genes.use)) {
    data.use <- t(PrepDR(
      object = object,
      features.use = features.use))
  }
  set.seed(seed = seed.use)
  if (tsne.method == "Rtsne") {
    if (is.null(x = distance.matrix)) {
      data.tsne <- Rtsne(
        X = as.matrix(x = data.use),
        dims = dim.embed,
        pca = FALSE, ...
      )
    } else {
      data.tsne <- Rtsne(
        X = as.matrix(x = distance.matrix),
        dims = dim.embed,
        is_distance = TRUE,
        ...
      )
    }
    data.tsne <- data.tsne$Y
  } else if (tsne.method == "FIt-SNE" & is.null(x = distance.matrix)) {
    data.tsne <- fftRtsne(X = as.matrix(x = data.use), dims = dim.embed, rand_seed = seed.use, ...)
  } else if (tsne.method == "tsne") {
    data.tsne <- tsne(X = data.use, k = dim.embed, ...)
  } else {
    stop("Invalid tsne.method: Please select from Rtsne, tsne, or FIt-SNE")
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
  if (!is.null(GetCalcParam(object = object, calculation = "RunTSNE", parameter = "genes.use"))) {
    object@calc.params$RunTSNE$genes.use <- colnames(data.use)
    object@calc.params$RunTSNE$cells.use <- rownames(data.use)
  }
  return(object)
}

#' @export
#' @method RunTSNE matrix
#'
RunTSNE.matrix <- function(
  object,
  assay.use = NULL,
  seed.use = 1,
  tsne.method = "Rtsne",
  add.iter = 0,
  dim.embed = 2,
  distance.matrix = NULL,
  reduction.name = "tsne",
  reduction.key = "tSNE_",
  ...
) {
  set.seed(seed = seed.use)
  tsne.data <- switch(
    EXPR = tsne.method,
    'Rtsne' = Rtsne(
      X = object,
      dims = dim.embed,
      ... # PCA/is_distance
    )$Y,
    'FIt-SNE' = fftRtsne(X = object, dims = dim.embed, rand_seed = seed.use, ...),
    'tsne' = tsne(X = object, k = dim.embed, ...),
    stop("Invalid tSNE method: please choose from 'Rtsne', 'FIt-SNE', or 'tsne'")
  )
  if (add.iter > 0) {
    tsne.data <- tsne(
      X = object,
      intial_config = as.matrix(x = tsne.data),
      max_iter = add.iter,
      ...
    )
  }
  colnames(x = tsne.data) <- paste0(reduction.key, 1:ncol(x = tsne.data))
  rownames(x = tsne.data) <- rownames(x = object)
  tsne.reduction <- MakeDimReducObject(
    cell.embeddings = tsne.data,
    key = reduction.key,
    assay.used = assay.use
  )
  return(tsne.reduction)
}

#' @param dims.use Which dimensions to use as input features
#'
#' @describeIn RunTSNE Run tSNE on a DimReduc object
#' @export
#' @method RunTSNE DimReduc
#'
RunTSNE.DimReduc <- function(
  object,
  cells.use = NULL,
  dims.use = 1:5,
  seed.use = 1,
  tsne.method = "Rtsne",
  add.iter = 0,
  dim.embed = 2,
  reduction.key = "tSNE_",
  ...
) {
  tsne.reduction <- RunTSNE(
    object = object[[, dims.use]],
    assay.use = DefaultAssay(object = object),
    seed.use = seed.use,
    tsne.method = tsne.method,
    add.iter = add.iter,
    dim.embed = dim.embed,
    reduction.key = reduction.key,
    ...
  )
  return(tsne.reduction)
}

#' @param distance.matrix If set, runs tSNE on the given distance matrix
#' instead of data matrix (experimental)
#' @param features.use If set, run the tSNE on this subset of features
#' (instead of running on a set of reduced dimensions). Not set (NULL) by default
#' @param reduction.name dimensional reduction name, specifies the position in the object$dr list. tsne by default
#' @param workflow.name Name of workflow
#'
#' @describeIn RunTSNE Run tSNE on a Seurat object
#' @export
#' @method RunTSNE Seurat
#'
RunTSNE.Seurat <- function(
  object,
  reduction.use = "pca",
  cells.use = NULL,
  dims.use = 1:5,
  features.use = NULL,
  seed.use = 1,
  tsne.method = "Rtsne",
  add.iter = 0,
  dim.embed = 2,
  distance.matrix = NULL,
  reduction.name = "tsne",
  reduction.key = "tSNE_",
  workflow.name = NULL,
  ...
) {
  if (!is.null(workflow.name)) {
    object <- PrepareWorkflow(object = object, workflow.name = workflow.name)
  }
  if (length(dims.use) == 1 && !is.null(workflow.name)) {
    dims.use <- 1:dims.use
  }
  tsne.reduction <- if (tsne.method == 'Rtsne') {
    if (!is.null(x = distance.matrix)) {
      RunTSNE(
        object = as.matrix(x = distance.matrix),
        assay.use = DefaultAssay(object = object),
        seed.use = seed.use,
        tsne.method = tsne.method,
        add.iter = add.iter,
        dim.embed = dim.embed,
        reduction.key = reduction.key,
        is_distance = TRUE,
        ...
      )
    } else if (!is.null(x = features.use)) {
      RunTSNE(
        object = as.matrix(x = GetAssayData(object = object)[features.use, ]),
        assay.use = DefaultAssay(object = object),
        seed.use = seed.use,
        tsne.method = tsne.method,
        add.iter = add.iter,
        dim.embed = dim.embed,
        reduction.key = reduction.key,
        pca = FALSE,
        ...
      )
    } else {
      RunTSNE(
        object = object[[reduction.use]],
        dims.use = dims.use,
        seed.use = seed.use,
        tsne.method = tsne.method,
        add.iter = add.iter,
        dim.embed = dim.embed,
        reduction.key = reduction.key,
        pca = FALSE,
        ...
      )
    }
  } else {
    RunTSNE(
      object = object[[reduction.use]],
      dims.use = dims.use,
      seed.use = seed.use,
      tsne.method = tsne.method,
      add.iter = add.iter,
      dim.embed = dim.embed,
      reduction.key = reduction.key,
      ...
    )
  }
  object[[reduction.name]] <- tsne.reduction
  object <- LogSeuratCommand(object = object)
  if (!is.null(workflow.name)) {
    command.name <- LogSeuratCommand(object = object, return.command = TRUE)
    object <- UpdateWorkflow(
      object = object,
      workflow.name = workflow.name,
      command.name = command.name)
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
#' @param reduction.use Reduction to use
#' @param assay.use Assay to use
#' @param dims.print Number of dims to print features for
#' @param features.print Number of features with highest/lowest loadings to print for
#' each dimension
#' @param overwrite Replace the existing data in feature.loadings
#' @param do.center Center the dataset prior to projection (should be set to TRUE)
#' @param verbose Print top genes associated with the projected dimensions
#'
#' @return Returns Seurat object with the projected values
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
  reduction.use = "pca",
  assay.use = NULL,
  dims.print = 1:5,
  features.print = 20,
  overwrite = FALSE,
  do.center = FALSE,
  verbose = TRUE
) {
  reduction <- object[[reduction.use]]
  assay.use <- assay.use %||% GetDimReduc(object = reduction, slot = "assay.used")
  data.use <- GetAssayData(
    object = object,
    assay.us = assay.use,
    slot = "scale.data"
  )
  if (do.center) {
    data.use <- scale(x = as.matrix(x = data.use), center = TRUE, scale = FALSE)
  }
  cell.embeddings <- GetDimReduc(
    object = reduction,
    slot = "cell.embeddings"
  )
  new.feature.loadings.full <- FastMatMult(m1 = data.use, m2 = cell.embeddings)
  rownames(x = new.feature.loadings.full) <- rownames(x = data.use)
  colnames(x = new.feature.loadings.full) <- colnames(x = cell.embeddings)
  reduction <- SetDimReduc(
    object = reduction,
    slot = "feature.loadings.projected",
    new.data = new.feature.loadings.full
  )
  if (overwrite) {
    reduction <- SetDimReduc(
      object = reduction,
      slot = "feature.loadings",
      new.data = new.feature.loadings.full
    )
  }
  object[[reduction.use]] <- reduction
  if (verbose) {
    Print(
      object = reduction,
      dims = dims.print,
      num.features = features.print,
      projected = TRUE
    )
  }
  object <- LogSeuratCommand(object = object)
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

#' @param standardize Standardize matrices - scales columns to have unit variance
#' and mean 0
#'
#' @describeIn RunCCA Run diagonal CCA on matrices
#' @export
#' @method RunCCA default
#'
RunCCA.default <- function(
  object1,
  object2,
  standardize = TRUE,
  num.cc = 20,
  verbose = FALSE,
  use.cpp=TRUE
) {
  set.seed(seed = 42)
  cells1 <- colnames(object1)
  cells2 <- colnames(object2)
  if (standardize) {
    object1 <- Standardize(mat = object1, display_progress = FALSE)
    object2 <- Standardize(mat = object2, display_progress = FALSE)
  }
  if (use.cpp == TRUE) {
    mat3 <- FastMatMult(m1 = t(x = object1), m2 = object2)
  }
  else {
    mat3 <- crossprod(x = object1, y = object2)
  }
  cca.svd <- irlba(A = mat3, nv = num.cc)
  cca.data <- rbind(cca.svd$u, cca.svd$v)
  colnames(x = cca.data) <- paste0("CC", 1:num.cc)
  rownames(cca.data) <- c(cells1, cells2)
  cca.data <- apply(cca.data, MARGIN = 2, function(x){
    if(sign(x[1]) == -1) {
      x <- x * -1
    }
    return(x)
  })

  return(list(ccv = cca.data, d = cca.svd$d))
}

#' @param assay.use1 Assay to pull from in the first object
#' @param assay.use1 Assay to pull from in the second object
#' @param features.use Set of genes to use in CCA. Default is the union of both
#' the variable features sets present in both objects.
#' @param renormlize Renormalize raw data after merging the objects. If FALSE,
#' merge the data matrices also.
#' @param rescale Rescale the datasets prior to CCA. If FALSE, uses existing data in the scale data slots.
#' @param compute.gene.loadings Also compute the gene loadings. NOTE - this will
#' scale every gene in the dataset which may impose a high memory cost.
#' @param ... Extra parameters (passed onto MergeSeurat in case with two objects
#' passed, passed onto ScaleData in case with single object and rescale.groups
#' set to TRUE)
#' @describeIn RunCCA Run CCA on a Seurat object
#' @export
#' @method RunCCA Seurat
#'
RunCCA.Seurat <- function(
  object1,
  object2,
  assay.use1 = NULL,
  assay.use2 = NULL,
  num.cc = 20,
  features.use = NULL,
  renormalize = FALSE,
  rescale = FALSE,
  compute.gene.loadings = TRUE,
  add.cell.id1 = NULL,
  add.cell.id2 = NULL,
  verbose = TRUE,
  use.cpp = TRUE,
  ...
) {
  assay.use1 <- assay.use1 %||% DefaultAssay(object = object1)
  assay.use2 <- assay.use2 %||% DefaultAssay(object = object2)

  if(assay.use1 != assay.use2) {
    stop("assay.use1 != assay.use2. Only unimodal data currently supported.")
  }

  if (is.null(x = features.use)) {
    if (length(x = VariableFeatures(object = object1, assay.use = assay.use1)) == 0) {
      stop(paste0("VariableFeatures not computed for the ", assay.use1, " assay in object1"))
    }
    if (length(x = VariableFeatures(object = object2, assay.use = assay.use2)) == 0) {
      stop(paste0("VariableFeatures not computed for the ", assay.use2, " assay in object2"))
    }
    features.use <- union(x = VariableFeatures(object = object1), y = VariableFeatures(object = object2))
    if (length(x = features.use) == 0) {
      stop("Zero features in the union of the VariableFeature sets ")
    }
  }
  data.use1 <- GetAssayData(object = object1, assay.use = assay.use1, slot = "scale.data")
  data.use2 <- GetAssayData(object = object2, assay.use = assay.use2, slot = "scale.data")
  features.use <- CheckFeatures(data.use = data.use1, features.use = features.use, object.name = "object1")
  features.use <- CheckFeatures(data.use = data.use2, features.use = features.use, object.name = "object2")


  if (length(x = features.use) < 50) {
    warning("Fewer than 50 features used as input for CCA.")
  }
  if (verbose) {
    message("Running CCA")
  }
  data1 <- data.use1[features.use,]
  data2 <- data.use2[features.use,]
  if (rescale) {
    if (verbose) message("Rescaling groups")
    data1 <- FastRowScale(data1)
    dimnames(data1) <- list(features.use, colnames(x = object1))
    data2 <- FastRowScale(data2)
    dimnames(data2) <- list(features.use, colnames(x = object2))
  }
  cca.results <- RunCCA(
    object1 = data1,
    object2 = data2,
    standardize = TRUE,
    num.cc = num.cc,
    verbose = verbose,
    use.cpp = use.cpp
  )
  if (verbose) {
   message("Merging objects")
  }
  combined.object <- merge(
    x = object1,
    y = object2,
    ...
  )

  combined.object[['cca']] <- MakeDimReducObject(
    cell.embeddings = cca.results$ccv,
    assay.used = assay.use1,
    key = "CC"
  )
  combined.scale <- cbind(data1,data2)
  combined.object <- SetAssayData(object = combined.object,new.data = combined.scale, assay.use = assay.use1,slot = "scale.data")
  if (renormalize) {
    combined.object <- NormalizeData(
      object = combined.object,
      assay.use = assay.use1,
      normalization.method = object1[[paste0("NormalizeData.", assay.use1)]]$normalization.method,
      scale.factor = object1[[paste0("NormalizeData.", assay.use1)]]$scale.factor
    )
  }
  if (compute.gene.loadings) {
    combined.object <- ProjectDim(
      object = combined.object,
      reduction.use = "cca",
      verbose = FALSE,
      overwrite = TRUE)
  }
  return(combined.object)
}

#' @describeIn RunMultiCCA Run mCCA on a list of matrices
#' @export
#' @method RunMultiCCA default
#'
RunMultiCCA.default <- function(
  object.list,
  niter = 25,
  num.ccs = 1,
  standardize = TRUE,
  verbose = TRUE
) {
  cell.names <- c()
  set.seed(seed = 42)
  for(object in object.list) {
    cell.names <- c(cell.names, colnames(x = object))
    if (!class(x = object) %in% c("matrix", "dgCMatrix")) {
      stop("Not all objects in object.list are matrices")
    }
  }
  num.sets <- length(x = object.list)
  if (standardize){
    for (i in 1:num.sets){
      object.list[[i]] <- Standardize(object.list[[i]], display_progress = FALSE)
    }
  }
  ws <- list()
  for (i in 1:num.sets){
    ws[[i]] <- irlba(object.list[[i]], nv = num.ccs)$v[, 1:num.ccs, drop = FALSE]
  }
  ws.init <- ws
  ws.final <- list()
  cors <- NULL
  for(i in 1:length(x = ws)) {
    ws.final[[i]] <- matrix(0, nrow = ncol(object.list[[i]]), ncol = num.ccs)
  }
  for (cc in 1:num.ccs){
    if (verbose) {
      message("Computing CC", cc)
    }
    ws <- list()
    for (i in 1:length(x = ws.init)){
      ws[[i]] <- ws.init[[i]][, cc]
    }
    cur.iter <- 1
    crit.old <- -10
    crit <- -20
    storecrits <- NULL
    while(cur.iter <= niter && abs(crit.old - crit)/abs(x = crit.old) > 0.001 && crit.old != 0){
      crit.old <- crit
      crit <- GetCrit(mat.list = object.list, ws = ws, num.sets = num.sets)
      storecrits <- c(storecrits, crit)
      cur.iter <- cur.iter + 1
      for(i in 1:num.sets){
        ws[[i]] <- UpdateW(mat.list = object.list, i = i, num.sets = num.sets, ws = ws, ws.final = ws.final)
      }
    }
    for(i in 1:length(x = ws)){
      ws.final[[i]][, cc] <- ws[[i]]
    }
    cors <- c(cors, GetCors(mat.list = object.list, ws = ws, num.sets = num.sets))
  }
  cca.data <- ws.final[[1]]
  for(i in 2:length(x = object.list)){
    cca.data <- rbind(cca.data, ws.final[[i]])
  }
  rownames(cca.data) <- cell.names
  cca.data <- apply(cca.data, MARGIN = 2, function(x){
    if(sign(x[1]) == -1) {
      x <- x * -1
    }
    return(x)
  })
  results <- list(
    ccv = cca.data,
    cors = cors
  )
  return(results)
}

#' @param assay.use Assay to use
#' @param features.use Set of genes to use in CCA. Default is the union of both
#' the variable features sets present in both objects.
#' @param renormlize Renormalize raw data after merging the objects. If FALSE,
#' merge the data matrices also.
#' @param compute.gene.loadings Also compute the gene loadings. NOTE - this will
#' scale every gene in the dataset which may impose a high memory cost.
#' @describeIn RunMultiCCA Run mCCA on a Seurat object
#' @export
#' @method RunMultiCCA Seurat
#'
RunMultiCCA.Seurat <- function(
  object.list,
  assay.use = NULL,
  features.use = NULL,
  add.cell.ids = NULL,
  niter = 25,
  num.ccs = 1,
  standardize = TRUE,
  renormalize = TRUE,
  compute.gene.loadings = TRUE,
  verbose = TRUE
) {
  set.seed(seed = 42)
  for(object in object.list) {
    if (class(x = object) != "Seurat") {
      stop("Not all objects in object.list are Seurat objects")
    }
  }
  assay.use <- assay.use %||% DefaultAssay(object = object.list[[1]])
  for(object in object.list) {
    assays <- FilterObjects(object = object, classes.keep = "Assay")
    if (!assay.use %in% assays) {
      stop(paste0(assay.use, " not present in all objects."))
    }
  }
  if (is.null(x = features.use)) {
    features.use <- c()
    for(object in object.list) {
      features.use <- c(features.use, VariableFeatures(object = object))
    }
    features.use <- unique(x = features.use)
  }
  for(i in 1:length(x = object.list)) {
    features.use <- CheckFeatures(
      data.use = GetAssayData(object = object.list[[i]], assay.use = assay.use, slot = "scale.data"),
      features.use = features.use,
      object.name = paste0("object", i)
    )
  }
  if (!is.null(x = add.cell.ids)) {
    if (length(x = add.cell.ids) != length(x = object.list)) {
      stop("add.cell.ids must have the same length as object.list")
    }
    object.list <- lapply(
      X = 1:length(x = object.list),
      FUN = function(x) {
        RenameCells(object = object.list[[x]], add.cell.id = add.cell.ids[x])
      }
    )
  }
  all.cell.names <- unlist(lapply(
    X = 1:length(x = object.list),
    FUN = function(x){
      colnames(x = object.list[[x]])
    })
  )
  if (anyDuplicated(x = all.cell.names)) {
    stop("Duplicate cell names detected, please set 'add.cell.ids'")
  }

  mat.list <- list()
  for(i in 1:length(x = object.list)) {
    mat.list[[i]] <- GetAssayData(
      object = object.list[[i]],
      assay.use = assay.use,
      slot = "scale.data"
    )[features.use, ]
  }
  cca.results <- RunMultiCCA(
    object.list = mat.list,
    niter = niter,
    num.ccs = num.ccs,
    standardize = standardize,
    verbose = verbose
  )
  combined.object <- merge(
    x = object.list[[1]],
    y = object.list[2:length(x = object.list)]
  )
  combined.object[['cca']] <- MakeDimReducObject(
    cell.embeddings = cca.results$ccv,
    assay.used = assay.use,
    key = "CC"
  )
  if (renormalize) {
    combined.object <- NormalizeData(
      object = combined.object,
      assay.use = assay.use,
      normalization.method = object.list[[1]][[paste0("NormalizeData.", assay.use)]]$normalization.method,
      scale.factor = object.list[[1]][[paste0("NormalizeData.", assay.use)]]$scale.factor,
      verbose = verbose
    )
    if (compute.gene.loadings) {
      combined.object <- ScaleData(object = combined.object, verbose = verbose)
      combined.object <- ProjectDim(
        object = combined.object,
        reduction.use = "cca",
        verbose = FALSE,
        overwrite = TRUE)
    }
  }
  combined.object <- LogSeuratCommand(object = combined.object)
  return(combined.object)
}

#' Run diffusion map
#'
#' NOTE: Prior to v2.3.4, this function used the R package diffusionMap to compute
#' the diffusion map components. This package was being archived and thus
#' RunDiffusion now uses the destiny package for the diffusion computations.
#' Please be aware that this will result in different default values as the two
#' underlying package implementations are different.
#'
#' @param object Seurat object
#' @param cells.use Which cells to analyze (default, all cells)
#' @param dims.use Which dimensions to use as input features
#' @param genes.use If set, run the diffusion map procedure on this subset of
#' genes (instead of running on a set of reduced dimensions). Not set (NULL) by
#' default
#' @param reduction.use Which dimensional reduction (PCA or ICA) to use for the
#' diffusion map input. Default is PCA
#' @param q.use Quantile to clip diffusion map components at. This addresses an
#' issue where 1-2 cells will have extreme values that obscure all other points.
#' 0.01 by default
#' @param max.dim Max dimension to keep from diffusion calculation
#' @param scale.clip Max/min value for scaled data. Default is 3
#' @param reduction.name dimensional reduction name, specifies the position in
#' the object$dr list. dm by default
#' @param reduction.key dimensional reduction key, specifies the string before
#' the number for the dimension names. DM by default
#' @param ... Additional arguments to the DiffusionMap call
#'
#' @return Returns a Seurat object with a diffusion map
#'
#' @importFrom utils installed.packages
#' @importFrom stats dist quantile
#'
#' @export
#'
#' @examples
#' \dontrun{
#' pbmc_small
#' # Run Diffusion on variable genes
#' pbmc_small <- RunDiffusion(pbmc_small,genes.use = pbmc_small@var.genes)
#' # Run Diffusion map on first 10 PCs
#' pbmc_small <- RunDiffusion(pbmc_small,genes.use = pbmc_small@var.genes)
#' # Plot results
#' DMPlot(pbmc_small)
#' }
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
  # Check for destiny
  if (!'destiny' %in% rownames(x = installed.packages())) {
    stop("Please install destiny - learn more at https://bioconductor.org/packages/release/bioc/html/destiny.html")
  }
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
    destiny::DiffusionMap(data = as.matrix(data.dist),
                          n_eigs = max.dim, ...)@eigenvectors
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

#' Run PHATE
#'
#' PHATE is a data reduction method specifically designed for visualizing
#' **high** dimensional data in **low** dimensional spaces.
#' To run, you must first install the `phate` python
#' package (e.g. via pip install phate). Details on this package can be
#' found here: \url{https://github.com/KrishnaswamyLab/PHATE}. For a more in depth
#' discussion of the mathematics underlying PHATE, see the bioRxiv paper here:
#' \url{https://www.biorxiv.org/content/early/2017/12/01/120378}.
#'
#' @param object Seurat object
#' @param cells.use Which cells to analyze (default, all cells)
#' @param genes.use If set, run PHATE on this subset of genes.
#' Not set (NULL) by default
#' @param assay.type Assay to pull data for (default: 'RNA')
#' @param max.dim Total number of dimensions to embed in PHATE.
#' @param k int, optional, default: 15
#' number of nearest neighbors on which to build kernel
#' @param alpha int, optional, default: 10
#' sets decay rate of kernel tails.
#' If NA, alpha decaying kernel is not used
#' @param use.alpha boolean, default: NA
#' forces the use of alpha decaying kernel
#' If NA, alpha decaying kernel is used for small inputs
#' (n_samples < n_landmark) and not used otherwise
#' @param n.landmark int, optional, default: 2000
#' number of landmarks to use in fast PHATE
#' @param potential.method string, optional, default: 'log'
#' choose from 'log' and 'sqrt'
#' which transformation of the diffusional operator is used
#' to compute the diffusion potential
#' @param t int, optional, default: 'auto'
#' power to which the diffusion operator is powered
#' sets the level of diffusion
#' @param knn.dist.method string, optional, default: 'euclidean'.
#' The desired distance function for calculating pairwise distances on the data.
#' If 'precomputed', `data` is treated as a
#' (n_samples, n_samples) distance or affinity matrix
#' @param mds.method string, optional, default: 'metric'
#' choose from 'classic', 'metric', and 'nonmetric'
#' which MDS algorithm is used for dimensionality reduction
#' @param mds.dist.method string, optional, default: 'euclidean'
#' recommended values: 'euclidean' and 'cosine'
#' @param t.max int, optional, default: 100.
#' Maximum value of t to test for automatic t selection.
#' @param npca int, optional, default: 100
#' Number of principal components to use for calculating
#' neighborhoods. For extremely large datasets, using
#' n_pca < 20 allows neighborhoods to be calculated in
#' log(n_samples) time.
#' @param plot.optimal.t boolean, optional, default: FALSE
#' If TRUE, produce a plot showing the Von Neumann Entropy
#' curve for automatic t selection.
#' @param verbose `int` or `boolean`, optional (default : 1)
#' If `TRUE` or `> 0`, print verbose updates.
#' @param n.jobs `int`, optional (default: 1)
#' The number of jobs to use for the computation.
#' If -1 all CPUs are used. If 1 is given, no parallel computing code is
#' used at all, which is useful for debugging.
#' For n_jobs below -1, (n.cpus + 1 + n.jobs) are used. Thus for
#' n_jobs = -2, all CPUs but one are used
#' @param seed.use int or `NA`, random state (default: `NA`)
#' @param reduction.name dimensional reduction name, specifies the position in
#' the object$dr list. phate by default
#' @param reduction.key dimensional reduction key, specifies the string before
#' the number for the dimension names. PHATE by default
#' @param ... Additional arguments for `phateR::phate`
#'
#' @return Returns a Seurat object containing a PHATE representation
#'
#' @importFrom utils installed.packages
#' @export
#'
#' @references Moon K, van Dijk D, Wang Z, Burkhardt D, Chen W, van den Elzen A,
#' Hirn M, Coifman R, Ivanova N, Wolf G and Krishnaswamy S (2017).
#' "Visualizing Transitions and Structure for High Dimensional Data
#' Exploration." _bioRxiv_, pp. 120378. doi: 10.1101/120378
#' (URL: http://doi.org/10.1101/120378),
#' <URL: https://www.biorxiv.org/content/early/2017/12/01/120378>.
#' @examples
#' if (reticulate::py_module_available("phate")) {
#'
#' # Load data
#' pbmc_small
#'
#' # Run PHATE with default parameters
#' pbmc_small <- RunPHATE(object = pbmc_small)
#' # Plot results
#' DimPlot(object = pbmc_small, reduction.use = 'phate')
#'
#' # Try smaller `k` for a small dataset, and larger `t` for a noisy embedding
#' pbmc_small <- RunPHATE(object = pbmc_small, k = 4, t = 12)
#' # Plot results
#' DimPlot(object = pbmc_small, reduction.use = 'phate')
#'1
#' # For increased emphasis on local structure, use sqrt potential
#' pbmc_small <- RunPHATE(object = pbmc_small, potential.method='sqrt')
#' # Plot results
#' DimPlot(object = pbmc_small, reduction.use = 'phate')
#' }
#'
RunPHATE <- function(
  object,
  cells.use = NULL,
  genes.use = NULL,
  assay.type = 'RNA',
  max.dim = 2L,
  k = 15,
  alpha = 10,
  use.alpha = NA,
  n.landmark = 2000,
  potential.method = "log",
  t = "auto",
  knn.dist.method = "euclidean",
  mds.method = "metric",
  mds.dist.method = "euclidean",
  t.max = 100,
  npca = 100,
  plot.optimal.t = FALSE,
  verbose = 1,
  n.jobs = 1,
  seed.use = NA,
  reduction.name = "phate",
  reduction.key = "PHATE",
  ...
) {
  if (!'phateR' %in% rownames(x = installed.packages())) {
    stop("Please install phateR")
  }
  data.use <- GetAssayData(object, assay.type = assay.type, slot = "scale.data")
  if (!is.null(x = cells.use)) {
    data.use <- data.use[, cells.use]
  }
  if (!is.null(x = genes.use)) {
    data.use <- data.use[genes.use, ]
  }
  data.use <- t(x = data.use)
  parameters.to.store <- as.list(x = environment(), all = TRUE)[names(x = formals(fun = "RunPHATE"))]
  object <- SetCalcParams(
    object = object,
    calculation = "RunPHATE",
    ... = parameters.to.store
  )
  phate_output <- phateR::phate(
    data.use,
    ndim = max.dim,
    k = k,
    alpha = alpha,
    use.alpha = alpha,
    n.landmark = n.landmark,
    potential.method = potential.method,
    t = t,
    knn.dist.method = knn.dist.method,
    init = NULL,
    mds.method = mds.method,
    mds.dist.method = mds.dist.method,
    t.max = t.max,
    npca = npca,
    plot.optimal.t = plot.optimal.t,
    verbose = verbose,
    n.jobs = n.jobs,
    seed = seed.use,
    ...
  )
  phate_output <- as.matrix(x = phate_output)
  colnames(x = phate_output) <- paste0(reduction.key, 1:ncol(x = phate_output))
  object <- SetDimReduction(
    object = object,
    reduction.type = reduction.name,
    slot = "cell.embeddings",
    new.data = phate_output
  )
  object <- SetDimReduction(
    object = object,
    reduction.type = reduction.name,
    slot = "key",
    new.data = reduction.key
  )
  return(object)
}

#' @describeIn RunUMAP Run a UMAP on a Seurat object
#' @export
#' @method RunUMAP Seurat
#'
RunUMAP.Seurat <- function(
  object,
  cells.use = NULL,
  dims.use = 1:5,
  reduction.use = 'pca',
  genes.use = NULL,
  assay.use = 'RNA',
  max.dim = 2L,
  reduction.name = "umap",
  reduction.key = "UMAP",
  n_neighbors = 30L,
  min_dist = 0.3,
  metric = "correlation",
  seed.use = 42,
  ...
) {
  if (!py_module_available(module = 'umap')) {
    stop("Cannot find UMAP, please install through pip (e.g. pip install umap-learn).")
  }
  if (!is.null(x = seed.use)) {
    set.seed(seed = seed.use)
    py_set_seed(seed = seed.use)
  }
  cells.use <- SetIfNull(x = cells.use, default = colnames(x = object@data))
  if (is.null(x = genes.use)) {
    dim.code <- GetDimReduction(
      object = object,
      reduction.type = reduction.use,
      slot = 'key'
    )
    dim.codes <- paste0(dim.code, dims.use)
    data.use <- GetDimReduction(
      object = object,
      reduction.type = reduction.use,
      slot = 'cell.embeddings'
    )
    data.use <- data.use[cells.use, dim.codes, drop = FALSE]
  } else {
    data.use <- GetAssayData(object = object, assay.type = assay.use, slot = 'scale.data')
    genes.use <- intersect(x = genes.use, y = rownames(x = data.use))
    if (!length(x = genes.use)) {
      stop("No genes found in the scale.data slot of assay ", assay.use)
    }
    data.use <- data.use[genes.use, cells.use, drop = FALSE]
    data.use <- t(x = data.use)
  }
  parameters.to.store <- as.list(x = environment(), all = TRUE)[names(formals("RunUMAP"))]
  object <- SetCalcParams(
    object = object,
    calculation = "RunUMAP",
    ... = parameters.to.store
  )
  umap_import <- import(module = "umap", delay_load = TRUE)
  umap <- umap_import$UMAP(
    n_neighbors = as.integer(x = n_neighbors),
    n_components = as.integer(x = max.dim),
    metric = metric,
    min_dist = min_dist
  )
  umap_output <- umap$fit_transform(as.matrix(x = data.use))
  colnames(x = umap_output) <- paste0(reduction.key, 1:ncol(x = umap_output))
  rownames(x = umap_output) <- cells.use
  object <- SetDimReduction(
    object = object,
    reduction.type = reduction.name,
    slot = "cell.embeddings",
    new.data = as.matrix(x = umap_output)
  )
  object <- SetDimReduction(
    object = object,
    reduction.type = reduction.name,
    slot = "key",
    new.data = reduction.key
  )
  return(object)
}
