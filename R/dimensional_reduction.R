#' @include generics.R
#'
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Determine statistical significance of PCA scores.
#'
#' Randomly permutes a subset of data, and calculates projected PCA scores for
#' these 'random' genes. Then compares the PCA scores for the 'random' genes
#' with the observed PCA scores to determine statistical signifance. End result
#' is a p-value for each gene's association with each principal component.
#'
#' @param object Seurat object
#' @param reduction DimReduc to use. ONLY PCA CURRENTLY SUPPORTED.
#' @param assay Assay used to calculate reduction.
#' @param dims Number of PCs to compute significance for
#' @param num.replicate Number of replicate samplings to perform
#' @param prop.freq Proportion of the data to randomly permute for each
#' replicate
#' @param verbose Print progress bar showing the number of replicates
#' that have been processed.
#' @param maxit maximum number of iterations to be performed by the irlba function of RunPCA
#'
#' @return Returns a Seurat object where JS(object = object[['pca']], slot = 'empirical')
#' represents p-values for each gene in the PCA analysis. If ProjectPCA is
#' subsequently run, JS(object = object[['pca']], slot = 'full') then
#' represents p-values for all genes.
#'
#' @importFrom methods new
#' @importFrom pbapply pblapply pbsapply
#' @importFrom future.apply future_lapply future_sapply
#' @importFrom future nbrOfWorkers
#'
#' @references Inspired by Chung et al, Bioinformatics (2014)
#' @concept dimensional_reduction
#'
#' @export
#'
#' @examples
#' \dontrun{
#' data("pbmc_small")
#' pbmc_small = suppressWarnings(JackStraw(pbmc_small))
#' head(JS(object = pbmc_small[['pca']], slot = 'empirical'))
#' }
#'
JackStraw <- function(
  object,
  reduction = "pca",
  assay = NULL,
  dims = 20,
  num.replicate = 100,
  prop.freq = 0.01,
  verbose = TRUE,
  maxit = 1000
) {
  if (reduction != "pca") {
    stop("Only pca for reduction is currently supported")
  }
  if (verbose && nbrOfWorkers() == 1) {
    my.lapply <- pblapply
    my.sapply <- pbsapply
  } else {
    my.lapply <- future_lapply
    my.sapply <- future_sapply
  }
  assay <- assay %||% DefaultAssay(object = object)
  if (IsSCT(assay = object[[assay]])) {
    stop("JackStraw cannot be run on SCTransform-normalized data.
         Please supply a non-SCT assay.")
  }
  if (dims > length(x = object[[reduction]])) {
    dims <- length(x = object[[reduction]])
    warning("Number of dimensions specified is greater than those available. Setting dims to ", dims, " and continuing", immediate. = TRUE)
  }
  if (dims > nrow(x = object)) {
    dims <- nrow(x = object)
    warning("Number of dimensions specified is greater than the number of cells. Setting dims to ", dims, " and continuing", immediate. = TRUE)
  }
  loadings <- Loadings(object = object[[reduction]], projected = FALSE)
  reduc.features <- rownames(x = loadings)
  if (length(x = reduc.features) < 3) {
    stop("Too few features")
  }
  if (length(x = reduc.features) * prop.freq < 3) {
    warning(
      "Number of variable genes given ",
      prop.freq,
      " as the prop.freq is low. Consider including more variable genes and/or increasing prop.freq. ",
      "Continuing with 3 genes in every random sampling."
    )
  }
  data.use <- GetAssayData(object = object, assay = assay, slot = "scale.data")[reduc.features, ]
  rev.pca <- object[[paste0('RunPCA.', assay)]]$rev.pca
  weight.by.var <- object[[paste0('RunPCA.', assay)]]$weight.by.var
  fake.vals.raw <- my.lapply(
    X = 1:num.replicate,
    FUN = JackRandom,
    scaled.data = data.use,
    prop.use = prop.freq,
    r1.use = 1,
    r2.use = dims,
    rev.pca = rev.pca,
    weight.by.var = weight.by.var,
    maxit = maxit
  )
  fake.vals <- sapply(
    X = 1:dims,
    FUN = function(x) {
      return(as.numeric(x = unlist(x = lapply(
        X = 1:num.replicate,
        FUN = function(y) {
          return(fake.vals.raw[[y]][, x])
        }
      ))))
    }
  )
  fake.vals <- as.matrix(x = fake.vals)
  jackStraw.empP <- as.matrix(
    my.sapply(
      X = 1:dims,
      FUN = function(x) {
        return(unlist(x = lapply(
          X = abs(loadings[, x]),
          FUN = EmpiricalP,
          nullval = abs(fake.vals[,x])
        )))
      }
    )
  )
  colnames(x = jackStraw.empP) <- paste0("PC", 1:ncol(x = jackStraw.empP))
  jackstraw.obj <- new(
    Class = "JackStrawData",
    empirical.p.values  = jackStraw.empP,
    fake.reduction.scores = fake.vals,
    empirical.p.values.full = matrix()
  )
  JS(object = object[[reduction]]) <- jackstraw.obj
  object <- LogSeuratCommand(object = object)
  return(object)
}

#' L2-normalization
#'
#' Perform l2 normalization on given dimensional reduction
#'
#' @param object Seurat object
#' @param reduction Dimensional reduction to normalize
#' @param new.dr name of new dimensional reduction to store
#' (default is olddr.l2)
#' @param new.key name of key for new dimensional reduction
#'
#' @return Returns a \code{\link{Seurat}} object
#' @concept dimensional_reduction
#'
#' @export
#'
L2Dim <- function(object, reduction, new.dr = NULL, new.key = NULL) {
  l2.norm <- L2Norm(mat = Embeddings(object[[reduction]]))
  if(is.null(new.dr)){
    new.dr <- paste0(reduction, ".l2")
  }
  if(is.null(new.key)){
    new.key <- paste0("L2", Key(object[[reduction]]))
  }
  colnames(x = l2.norm) <- paste0(new.key, 1:ncol(x = l2.norm))
  l2.dr <- CreateDimReducObject(
    embeddings = l2.norm,
    loadings = Loadings(object = object[[reduction]], projected = FALSE),
    projected = Loadings(object = object[[reduction]], projected = TRUE),
    assay = DefaultAssay(object = object),
    stdev = slot(object = object[[reduction]], name = 'stdev'),
    key = new.key,
    jackstraw = slot(object = object[[reduction]], name = 'jackstraw'),
    misc = slot(object = object[[reduction]], name = 'misc')
  )
  object[[new.dr]] <- l2.dr
  return(object)
}

#' L2-Normalize CCA
#'
#' Perform l2 normalization on CCs
#'
#' @param object Seurat object
#' @param \dots Additional parameters to L2Dim.
#' @concept dimensional_reduction
#'
#' @export
#'
L2CCA <- function(object, ...){
  CheckDots(..., fxns = 'L2Dim')
  return(L2Dim(object = object, reduction = "cca", ...))
}

#' Significant genes from a PCA
#'
#' Returns a set of genes, based on the JackStraw analysis, that have
#' statistically significant associations with a set of PCs.
#'
#' @param object Seurat object
#' @param pcs.use PCS to use.
#' @param pval.cut P-value cutoff
#' @param use.full Use the full list of genes (from the projected PCA). Assumes
#' that \code{ProjectDim} has been run. Currently, must be set to FALSE.
#' @param max.per.pc Maximum number of genes to return per PC. Used to avoid genes from one PC dominating the entire analysis.
#'
#' @return A vector of genes whose p-values are statistically significant for
#' at least one of the given PCs.
#'
#' @export
#' @concept dimensional_reduction
#'
#' @seealso \code{\link{ProjectDim}} \code{\link{JackStraw}}
#'
#' @examples
#' data("pbmc_small")
#' PCASigGenes(pbmc_small, pcs.use = 1:2)
#'
PCASigGenes <- function(
  object,
  pcs.use,
  pval.cut = 0.1,
  use.full = FALSE,
  max.per.pc = NULL
) {
  # pvals.use <- GetDimReduction(object,reduction.type = "pca",slot = "jackstraw")@empirical.p.values
  empirical.use <- ifelse(test = use.full, yes = 'full', no = 'empirical')
  pvals.use <- JS(object = object[['pca']], slot = empirical.use)
  if (length(x = pcs.use) == 1) {
    pvals.min <- pvals.use[, pcs.use]
  }
  if (length(x = pcs.use) > 1) {
    pvals.min <- apply(X = pvals.use[, pcs.use], MARGIN = 1, FUN = min)
  }
  names(x = pvals.min) <- rownames(x = pvals.use)
  features <- names(x = pvals.min)[pvals.min < pval.cut]
  if (!is.null(x = max.per.pc)) {
    top.features <- TopFeatures(
      object = object[['pca']],
      dim = pcs.use,
      nfeatures = max.per.pc,
      projected = use.full,
      balanced = FALSE
    )
    features <- intersect(x = top.features, y = features)
  }
  return(features)
}

#' Project Dimensional reduction onto full dataset
#'
#' Takes a pre-computed dimensional reduction (typically calculated on a subset
#' of genes) and projects this onto the entire dataset (all genes). Note that
#' the cell loadings will remain unchanged, but now there are gene loadings for
#' all genes.
#'
#' @param object Seurat object
#' @param reduction Reduction to use
#' @param assay Assay to use
#' @param dims.print Number of dims to print features for
#' @param nfeatures.print Number of features with highest/lowest loadings to print for
#' each dimension
#' @param overwrite Replace the existing data in feature.loadings
#' @param do.center Center the dataset prior to projection (should be set to TRUE)
#' @param verbose Print top genes associated with the projected dimensions
#'
#' @return Returns Seurat object with the projected values
#'
#' @export
#' @concept dimensional_reduction
#'
#' @examples
#' data("pbmc_small")
#' pbmc_small
#' pbmc_small <- ProjectDim(object = pbmc_small, reduction = "pca")
#' # Vizualize top projected genes in heatmap
#' DimHeatmap(object = pbmc_small, reduction = "pca", dims = 1, balanced = TRUE)
#'
ProjectDim <- function(
  object,
  reduction = "pca",
  assay = NULL,
  dims.print = 1:5,
  nfeatures.print = 20,
  overwrite = FALSE,
  do.center = FALSE,
  verbose = TRUE
) {
  redeuc <- object[[reduction]]
  assay <- assay %||% DefaultAssay(object = redeuc)
  data.use <- GetAssayData(
    object = object[[assay]],
    slot = "scale.data"
  )
  if (do.center) {
    data.use <- scale(x = as.matrix(x = data.use), center = TRUE, scale = FALSE)
  }
  cell.embeddings <- Embeddings(object = redeuc)
  new.feature.loadings.full <- data.use %*% cell.embeddings
  rownames(x = new.feature.loadings.full) <- rownames(x = data.use)
  colnames(x = new.feature.loadings.full) <- colnames(x = cell.embeddings)
  Loadings(object = redeuc, projected = TRUE) <- new.feature.loadings.full
  if (overwrite) {
    Loadings(object = redeuc, projected = FALSE) <- new.feature.loadings.full
  }
  object[[reduction]] <- redeuc
  if (verbose) {
    print(
      x = redeuc,
      dims = dims.print,
      nfeatures = nfeatures.print,
      projected = TRUE
    )
  }
  object <- LogSeuratCommand(object = object)
  return(object)
}

#' @param query.dims Dimensions (columns) to use from query
#' @param reference.dims Dimensions (columns) to use from reference
#' @param ... Additional parameters to \code{\link{RunUMAP}}
#'
#' @inheritParams FindNeighbors
#' @inheritParams RunUMAP
#'
#' @rdname ProjectUMAP
#' @concept dimensional_reduction
#' @export
#'
ProjectUMAP.default <- function(
  query,
  query.dims = NULL,
  reference,
  reference.dims = NULL,
  k.param = 30,
  nn.method = "annoy",
  n.trees = 50,
  annoy.metric = "cosine",
  l2.norm = FALSE,
  cache.index = TRUE,
  index = NULL,
  neighbor.name = "query_ref.nn",
  reduction.model,
  ...
) {
  query.dims <- query.dims %||% 1:ncol(x = query)
  reference.dims <- reference.dims %||% query.dims
 if (length(x = reference.dims) != length(x = query.dims)) {
    stop("Length of Reference and Query number of dimensions are not equal")
   }
  if (any(reference.dims > ncol(x = reference))) {
    stop("Reference dims is larger than the number of dimensions present.", call. = FALSE)
  }
  if (any(query.dims > ncol(x = query))) {
    stop("Query dims is larger than the number of dimensions present.", call. = FALSE)
  }
  if (length(x = Misc(object = reduction.model, slot = 'model')) == 0) {
    stop(
      "The provided reduction.model does not have a model stored. Please try ",
      "running umot-learn on the object first", call. = FALSE
    )
  }
  query.neighbor <- FindNeighbors(
    object = reference[, reference.dims],
    query = query[, query.dims],
    k.param = k.param,
    nn.method = nn.method,
    n.trees = n.trees,
    annoy.metric = annoy.metric,
    cache.index = cache.index,
    index = index,
    return.neighbor = TRUE,
    l2.norm = l2.norm
  )
  proj.umap <- RunUMAP(object = query.neighbor, reduction.model = reduction.model, ...)
  return(list(proj.umap = proj.umap, query.neighbor = query.neighbor))
}

#' @rdname ProjectUMAP
#' @concept dimensional_reduction
#' @export
#' @method ProjectUMAP DimReduc
#'
ProjectUMAP.DimReduc <- function(
  query,
  query.dims = NULL,
  reference,
  reference.dims = NULL,
  k.param = 30,
  nn.method = "annoy",
  n.trees = 50,
  annoy.metric = "cosine",
  l2.norm = FALSE,
  cache.index = TRUE,
  index = NULL,
  neighbor.name = "query_ref.nn",
  reduction.model,
  ...
) {
  proj.umap <- ProjectUMAP(
    query = Embeddings(object = query),
    query.dims = query.dims,
    reference = Embeddings(object = reference),
    reference.dims = reference.dims,
    k.param = k.param,
    nn.method = nn.method,
    n.trees = 50,
    annoy.metric = annoy.metric,
    l2.norm = l2.norm,
    cache.index = cache.index,
    index = index,
    neighbor.name = neighbor.name,
    reduction.model = reduction.model,
    ...
  )
  return(proj.umap)
}

#' @param reference Reference dataset
#' @param query.reduction Name of reduction to use from the query for neighbor
#' finding
#' @param reference.reduction Name of reduction to use from the reference for
#' neighbor finding
#' @param neighbor.name Name to store neighbor information in the query
#' @param reduction.name Name of projected UMAP to store in the query
#' @param reduction.key Value for the projected UMAP key
#' @rdname ProjectUMAP
#' @concept dimensional_reduction
#' @export
#' @method ProjectUMAP Seurat
#'
ProjectUMAP.Seurat <- function(
  query,
  query.reduction,
  query.dims = NULL,
  reference,
  reference.reduction,
  reference.dims = NULL,
  k.param = 30,
  nn.method = "annoy",
  n.trees = 50,
  annoy.metric = "cosine",
  l2.norm = FALSE,
  cache.index = TRUE,
  index = NULL,
  neighbor.name = "query_ref.nn",
  reduction.model,
  reduction.name = "ref.umap",
  reduction.key = "refUMAP_",
  ...
) {
  if (!query.reduction %in% Reductions(object = query)) {
    stop("The query.reduction (", query.reduction, ") is not present in the ",
         "provided query", call. = FALSE)
  }
  if (!reference.reduction %in% Reductions(object = reference)) {
    stop("The reference.reduction (", reference.reduction, ") is not present in the ",
         "provided reference.", call. = FALSE)
  }
  if (!reduction.model %in% Reductions(object = reference)) {
    stop("The reduction.model (", reduction.model, ") is not present in the ",
         "provided reference.", call. = FALSE)
  }
  proj.umap <- ProjectUMAP(
    query = query[[query.reduction]],
    query.dims = query.dims,
    reference = reference[[reference.reduction]],
    reference.dims = reference.dims,
    k.param = k.param,
    nn.method = nn.method,
    n.trees = n.trees,
    annoy.metric = annoy.metric,
    l2.norm = l2.norm,
    cache.index = cache.index,
    index = index,
    neighbor.name = neighbor.name,
    reduction.model = reference[[reduction.model]],
    reduction.key = reduction.key,
    assay = DefaultAssay(query),
    ...
  )
  query[[reduction.name]] <- proj.umap$proj.umap
  query[[neighbor.name]] <- proj.umap$query.neighbor
  return(query)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Methods for Seurat-defined generics
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @param standardize Standardize matrices - scales columns to have unit variance
#' and mean 0
#' @param num.cc Number of canonical vectors to calculate
#' @param seed.use Random seed to set. If NULL, does not set a seed
#' @param verbose Show progress messages
#'
#' @importFrom irlba irlba
#'
#' @rdname RunCCA
#' @concept dimensional_reduction
#' @export
#'
RunCCA.default <- function(
  object1,
  object2,
  standardize = TRUE,
  num.cc = 20,
  seed.use = 42,
  verbose = FALSE,
  ...
) {
  if (!is.null(x = seed.use)) {
    set.seed(seed = seed.use)
  }
  cells1 <- colnames(x = object1)
  cells2 <- colnames(x = object2)
  if (standardize) {
    object1 <- Standardize(mat = object1, display_progress = FALSE)
    object2 <- Standardize(mat = object2, display_progress = FALSE)
  }
  mat3 <- crossprod(x = object1, y = object2)
  cca.svd <- irlba(A = mat3, nv = num.cc)
  cca.data <- rbind(cca.svd$u, cca.svd$v)
  colnames(x = cca.data) <- paste0("CC", 1:num.cc)
  rownames(cca.data) <- c(cells1, cells2)
  cca.data <- apply(
    X = cca.data,
    MARGIN = 2,
    FUN = function(x) {
      if (sign(x[1]) == -1) {
        x <- x * -1
      }
      return(x)
    }
  )
  return(list(ccv = cca.data, d = cca.svd$d))
}

#' @param assay1,assay2 Assays to pull from in the first and second objects, respectively
#' @param features Set of genes to use in CCA. Default is the union of both
#' the variable features sets present in both objects.
#' @param renormalize Renormalize raw data after merging the objects. If FALSE,
#' merge the data matrices also.
#' @param rescale Rescale the datasets prior to CCA. If FALSE, uses existing data in the scale data slots.
#' @param compute.gene.loadings Also compute the gene loadings. NOTE - this will
#' scale every gene in the dataset which may impose a high memory cost.
#' @param add.cell.id1,add.cell.id2 Add ...
#' @param ... Extra parameters (passed onto MergeSeurat in case with two objects
#' passed, passed onto ScaleData in case with single object and rescale.groups
#' set to TRUE)
#'
#' @rdname RunCCA
#' @concept dimensional_reduction
#' @export
#' @method RunCCA Seurat
#'
RunCCA.Seurat <- function(
  object1,
  object2,
  assay1 = NULL,
  assay2 = NULL,
  num.cc = 20,
  features = NULL,
  renormalize = FALSE,
  rescale = FALSE,
  compute.gene.loadings = TRUE,
  add.cell.id1 = NULL,
  add.cell.id2 = NULL,
  verbose = TRUE,
  ...
) {
  assay1 <- assay1 %||% DefaultAssay(object = object1)
  assay2 <- assay2 %||% DefaultAssay(object = object2)
  if (assay1 != assay2) {
    warning("Running CCA on different assays")
  }
  if (is.null(x = features)) {
    if (length(x = VariableFeatures(object = object1, assay = assay1)) == 0) {
      stop(paste0("VariableFeatures not computed for the ", assay1, " assay in object1"))
    }
    if (length(x = VariableFeatures(object = object2, assay = assay2)) == 0) {
      stop(paste0("VariableFeatures not computed for the ", assay2, " assay in object2"))
    }
    features <- union(x = VariableFeatures(object = object1), y = VariableFeatures(object = object2))
    if (length(x = features) == 0) {
      stop("Zero features in the union of the VariableFeature sets ")
    }
  }
  nfeatures <- length(x = features)
  if (!(rescale)) {
    data.use1 <- GetAssayData(object = object1, assay = assay1, slot = "scale.data")
    data.use2 <- GetAssayData(object = object2, assay = assay2, slot = "scale.data")
    features <- CheckFeatures(data.use = data.use1, features = features, object.name = "object1", verbose = FALSE)
    features <- CheckFeatures(data.use = data.use2, features = features, object.name = "object2", verbose = FALSE)
    data1 <- data.use1[features, ]
    data2 <- data.use2[features, ]
  }
  if (rescale) {
    data.use1 <- GetAssayData(object = object1, assay = assay1, slot = "data")
    data.use2 <- GetAssayData(object = object2, assay = assay2, slot = "data")
    features <- CheckFeatures(data.use = data.use1, features = features, object.name = "object1", verbose = FALSE)
    features <- CheckFeatures(data.use = data.use2, features = features, object.name = "object2", verbose = FALSE)
    data1 <- data.use1[features,]
    data2 <- data.use2[features,]
    if (verbose) message("Rescaling groups")
    data1 <- FastRowScale(as.matrix(data1))
    dimnames(data1) <- list(features, colnames(x = object1))
    data2 <- FastRowScale(as.matrix(data2))
    dimnames(data2) <- list(features, colnames(x = object2))
  }
  if (length(x = features) / nfeatures < 0.1 & verbose) {
    warning("More than 10% of provided features filtered out. Please check that the given features are present in the scale.data slot for both the assays provided here and that they have non-zero variance.")
  }
  if (length(x = features) < 50) {
    warning("Fewer than 50 features used as input for CCA.")
  }
  if (verbose) {
    message("Running CCA")
  }
  cca.results <- RunCCA(
    object1 = data1,
    object2 = data2,
    standardize = TRUE,
    num.cc = num.cc,
    verbose = verbose,
  )
  if (verbose) {
    message("Merging objects")
  }
  combined.object <- merge(
    x = object1,
    y = object2,
    merge.data = TRUE,
    ...
  )
  rownames(x = cca.results$ccv) <- Cells(x = combined.object)
  colnames(x = data1) <- Cells(x = combined.object)[1:ncol(x = data1)]
  colnames(x = data2) <- Cells(x = combined.object)[(ncol(x = data1) + 1):length(x = Cells(x = combined.object))]
  combined.object[['cca']] <- CreateDimReducObject(
    embeddings = cca.results$ccv[colnames(combined.object), ],
    assay = assay1,
    key = "CC_"
  )
  combined.object[['cca']]@assay.used <- DefaultAssay(combined.object)
  if (ncol(combined.object) != (ncol(object1) + ncol(object2))) {
    warning("Some cells removed after object merge due to minimum feature count cutoff")
  }
  combined.scale <- cbind(data1,data2)
  combined.object <- SetAssayData(object = combined.object,new.data = combined.scale, slot = "scale.data")
  if (renormalize) {
    combined.object <- NormalizeData(
      object = combined.object,
      assay = assay1,
      normalization.method = object1[[paste0("NormalizeData.", assay1)]]$normalization.method,
      scale.factor = object1[[paste0("NormalizeData.", assay1)]]$scale.factor
    )
  }
  if (compute.gene.loadings) {
    combined.object <- ProjectDim(
      object = combined.object,
      reduction = "cca",
      verbose = FALSE,
      overwrite = TRUE)
  }
  return(combined.object)
}

#' @param assay Name of Assay ICA is being run on
#' @param nics Number of ICs to compute
#' @param rev.ica By default, computes the dimensional reduction on the cell x
#' feature matrix. Setting to true will compute it on the transpose (feature x cell
#' matrix).
#' @param ica.function ICA function from ica package to run (options: icafast,
#' icaimax, icajade)
#' @param verbose Print the top genes associated with high/low loadings for
#' the ICs
#' @param ndims.print ICs to print genes for
#' @param nfeatures.print Number of genes to print for each IC
#' @param reduction.key dimensional reduction key, specifies the string before
#' the number for the dimension names.
#' @param seed.use Set a random seed.  Setting NULL will not set a seed.
#' @param \dots Additional arguments to be passed to fastica
#'
#' @importFrom ica icafast icaimax icajade
#'
#' @rdname RunICA
#' @concept dimensional_reduction
#' @export
#' @method RunICA default
#'
RunICA.default <- function(
  object,
  assay = NULL,
  nics = 50,
  rev.ica = FALSE,
  ica.function = "icafast",
  verbose = TRUE,
  ndims.print = 1:5,
  nfeatures.print = 30,
  reduction.name = "ica",
  reduction.key = "ica_",
  seed.use = 42,
  ...
) {
  CheckDots(..., fxns = ica.function)
  if (!is.null(x = seed.use)) {
    set.seed(seed = seed.use)
  }
  nics <- min(nics, ncol(x = object))
  ica.fxn <- eval(expr = parse(text = ica.function))
  if (rev.ica) {
    ica.results <- ica.fxn(object, nc = nics,...)
    cell.embeddings <- ica.results$M
  } else {
    ica.results <- ica.fxn(t(x = object), nc = nics,...)
    cell.embeddings <- ica.results$S
  }
  feature.loadings <- (as.matrix(x = object ) %*% as.matrix(x = cell.embeddings))
  colnames(x = feature.loadings) <- paste0(reduction.key, 1:ncol(x = feature.loadings))
  colnames(x = cell.embeddings) <- paste0(reduction.key, 1:ncol(x = cell.embeddings))
  reduction.data <- CreateDimReducObject(
    embeddings = cell.embeddings,
    loadings = feature.loadings,
    assay = assay,
    key = reduction.key
  )
  if (verbose) {
    print(x = reduction.data, dims = ndims.print, nfeatures = nfeatures.print)
  }
  return(reduction.data)
}

#' @param features Features to compute ICA on
#'
#' @rdname RunICA
#' @concept dimensional_reduction
#' @export
#' @method RunICA Assay
#'
RunICA.Assay <- function(
  object,
  assay = NULL,
  features = NULL,
  nics = 50,
  rev.ica = FALSE,
  ica.function = "icafast",
  verbose = TRUE,
  ndims.print = 1:5,
  nfeatures.print = 30,
  reduction.name = "ica",
  reduction.key = "ica_",
  seed.use = 42,
  ...
) {
  data.use <- PrepDR(
    object = object,
    features = features,
    verbose = verbose
  )
  reduction.data <- RunICA(
    object = data.use,
    assay = assay,
    nics = nics,
    rev.ica = rev.ica,
    ica.function = ica.function,
    verbose = verbose,
    ndims.print = ndims.print,
    nfeatures.print = nfeatures.print,
    reduction.key = reduction.key,
    seed.use = seed.use,
    ...

  )
  return(reduction.data)
}

#' @param reduction.name dimensional reduction name
#'
#' @rdname RunICA
#' @concept dimensional_reduction
#' @method RunICA Seurat
#' @export
#'
RunICA.Seurat <- function(
  object,
  assay = NULL,
  features = NULL,
  nics = 50,
  rev.ica = FALSE,
  ica.function = "icafast",
  verbose = TRUE,
  ndims.print = 1:5,
  nfeatures.print = 30,
  reduction.name = "ica",
  reduction.key = "IC_",
  seed.use = 42,
  ...
) {
  assay <- assay %||% DefaultAssay(object = object)
  reduction.data <- RunICA(
    object = object[[assay]],
    assay = assay,
    features = features,
    nics = nics,
    rev.ica = rev.ica,
    ica.function = ica.function,
    verbose = verbose,
    ndims.print = ndims.print,
    nfeatures.print = nfeatures.print,
    reduction.key = reduction.key,
    seed.use = seed.use,
    ...
  )
  object[[reduction.name]] <- reduction.data
  object <- LogSeuratCommand(object = object)
  return(object)
}

#' @param assay Name of Assay PCA is being run on
#' @param npcs Total Number of PCs to compute and store (50 by default)
#' @param rev.pca By default computes the PCA on the cell x gene matrix. Setting
#' to true will compute it on gene x cell matrix.
#' @param weight.by.var Weight the cell embeddings by the variance of each PC
#' (weights the gene loadings if rev.pca is TRUE)
#' @param verbose Print the top genes associated with high/low loadings for
#' the PCs
#' @param ndims.print PCs to print genes for
#' @param nfeatures.print Number of genes to print for each PC
#' @param reduction.key dimensional reduction key, specifies the string before
#' the number for the dimension names. PC by default
#' @param seed.use Set a random seed. By default, sets the seed to 42. Setting
#' NULL will not set a seed.
#' @param approx Use truncated singular value decomposition to approximate PCA
#'
#' @importFrom irlba irlba
#' @importFrom stats prcomp
#' @importFrom utils capture.output
#'
#' @rdname RunPCA
#' @concept dimensional_reduction
#' @export
#'
RunPCA.default <- function(
  object,
  assay = NULL,
  npcs = 50,
  rev.pca = FALSE,
  weight.by.var = TRUE,
  verbose = TRUE,
  ndims.print = 1:5,
  nfeatures.print = 30,
  reduction.key = "PC_",
  seed.use = 42,
  approx = TRUE,
  ...
) {
  if (!is.null(x = seed.use)) {
    set.seed(seed = seed.use)
  }
 if (inherits(x = object, what = 'matrix')) {
   RowVar.function <- RowVar
 } else if (inherits(x = object, what = 'dgCMatrix')) {
   RowVar.function <- RowVarSparse
 } else if (inherits(x = object, what = 'IterableMatrix')) {
   RowVar.function <- function(x) {
     return(BPCells::matrix_stats(
       matrix = x,
       row_stats = 'variance'
     )$row_stats['variance',])
     }
 }
  if (rev.pca) {
    npcs <- min(npcs, ncol(x = object) - 1)
    pca.results <- irlba(A = object, nv = npcs, ...)
    total.variance <- sum(RowVar.function(x = t(x = object)))
    sdev <- pca.results$d/sqrt(max(1, nrow(x = object) - 1))
    if (weight.by.var) {
      feature.loadings <- pca.results$u %*% diag(pca.results$d)
    } else{
      feature.loadings <- pca.results$u
    }
    cell.embeddings <- pca.results$v
  }
  else {
    total.variance <- sum(RowVar.function(x = object))
    if (approx) {
      npcs <- min(npcs, nrow(x = object) - 1)
      pca.results <- irlba(A = t(x = object), nv = npcs, ...)
      feature.loadings <- pca.results$v
      sdev <- pca.results$d/sqrt(max(1, ncol(object) - 1))
      if (weight.by.var) {
        cell.embeddings <- pca.results$u %*% diag(pca.results$d)
      } else {
        cell.embeddings <- pca.results$u
      }
    } else {
      npcs <- min(npcs, nrow(x = object))
      pca.results <- prcomp(x = t(object), rank. = npcs, ...)
      feature.loadings <- pca.results$rotation
      sdev <- pca.results$sdev
      if (weight.by.var) {
        cell.embeddings <- pca.results$x
      } else {
        cell.embeddings <- pca.results$x / (pca.results$sdev[1:npcs] * sqrt(x = ncol(x = object) - 1))
      }
    }
  }
  rownames(x = feature.loadings) <- rownames(x = object)
  colnames(x = feature.loadings) <- paste0(reduction.key, 1:npcs)
  rownames(x = cell.embeddings) <- colnames(x = object)
  colnames(x = cell.embeddings) <- colnames(x = feature.loadings)
  reduction.data <- CreateDimReducObject(
    embeddings = cell.embeddings,
    loadings = feature.loadings,
    assay = assay,
    stdev = sdev,
    key = reduction.key,
    misc = list(total.variance = total.variance)
  )
  if (verbose) {
    msg <- capture.output(print(
      x = reduction.data,
      dims = ndims.print,
      nfeatures = nfeatures.print
    ))
    message(paste(msg, collapse = '\n'))
  }
  return(reduction.data)
}

#' @param features Features to compute PCA on. If features=NULL, PCA will be run
#' using the variable features for the Assay. Note that the features must be present
#' in the scaled data. Any requested features that are not scaled or have 0 variance
#' will be dropped, and the PCA will be run using the remaining features.
#'
#' @rdname RunPCA
#' @concept dimensional_reduction
#' @export
#' @method RunPCA Assay
#'
RunPCA.Assay <- function(
  object,
  assay = NULL,
  features = NULL,
  npcs = 50,
  rev.pca = FALSE,
  weight.by.var = TRUE,
  verbose = TRUE,
  ndims.print = 1:5,
  nfeatures.print = 30,
  reduction.key = "PC_",
  seed.use = 42,
  ...
) {
  data.use <- PrepDR(
    object = object,
    features = features,
    verbose = verbose
  )
  reduction.data <- RunPCA(
    object = data.use,
    assay = assay,
    npcs = npcs,
    rev.pca = rev.pca,
    weight.by.var = weight.by.var,
    verbose = verbose,
    ndims.print = ndims.print,
    nfeatures.print = nfeatures.print,
    reduction.key = reduction.key,
    seed.use = seed.use,
    ...

  )
  return(reduction.data)
}

#' @method RunPCA StdAssay
#' @export
#'
RunPCA.StdAssay <- function(
  object,
  assay = NULL,
  features = NULL,
  layer = 'scale.data',
  npcs = 50,
  rev.pca = FALSE,
  weight.by.var = TRUE,
  verbose = TRUE,
  ndims.print = 1:5,
  nfeatures.print = 30,
  reduction.key = "PC_",
  seed.use = 42,
  ...
) {
  data.use <- PrepDR5(
    object = object,
    features = features,
    layer = layer,
    verbose = verbose
  )
  return(RunPCA(
    object = data.use,
    assay = assay,
    npcs = npcs,
    rev.pca = rev.pca,
    weight.by.var = weight.by.var,
    verbose = verbose,
    ndims.print = ndims.print,
    nfeatures.print = nfeatures.print,
    reduction.key = reduction.key,
    seed.use = seed.use,
    ...
  ))
}

#' @param reduction.name dimensional reduction name,  pca by default
#'
#' @rdname RunPCA
#' @concept dimensional_reduction
#' @export
#' @method RunPCA Seurat
#'
RunPCA.Seurat <- function(
  object,
  assay = NULL,
  features = NULL,
  npcs = 50,
  rev.pca = FALSE,
  weight.by.var = TRUE,
  verbose = TRUE,
  ndims.print = 1:5,
  nfeatures.print = 30,
  reduction.name = "pca",
  reduction.key = "PC_",
  seed.use = 42,
  ...
) {
  assay <- assay %||% DefaultAssay(object = object)
  reduction.data <- RunPCA(
    object = object[[assay]],
    assay = assay,
    features = features,
    npcs = npcs,
    rev.pca = rev.pca,
    weight.by.var = weight.by.var,
    verbose = verbose,
    ndims.print = ndims.print,
    nfeatures.print = nfeatures.print,
    reduction.key = reduction.key,
    seed.use = seed.use,
    ...
  )
  object[[reduction.name]] <- reduction.data
  object <- LogSeuratCommand(object = object)
  return(object)
}

#' @method RunPCA Seurat5
#' @export
#'
RunPCA.Seurat5 <- function(
  object,
  assay = NULL,
  features = NULL,
  npcs = 50,
  rev.pca = FALSE,
  weight.by.var = TRUE,
  verbose = TRUE,
  ndims.print = 1:5,
  nfeatures.print = 30,
  reduction.name = "pca",
  reduction.key = "PC_",
  seed.use = 42,
  ...
) {
  assay <- assay %||% DefaultAssay(object = object)
  reduction.data <- RunPCA(
    object = object[[assay]],
    assay = assay,
    features = features,
    npcs = npcs,
    rev.pca = rev.pca,
    weight.by.var = weight.by.var,
    verbose = verbose,
    ndims.print = ndims.print,
    nfeatures.print = nfeatures.print,
    reduction.key = reduction.key,
    seed.use = seed.use,
    ...
  )
  object[[reduction.name]] <- reduction.data
  # object <- LogSeuratCommand(object = object)
  return(object)
}

#' @param assay Name of assay that that t-SNE is being run on
#' @param seed.use Random seed for the t-SNE. If NULL, does not set the seed
#' @param tsne.method Select the method to use to compute the tSNE. Available
#' methods are:
#' \itemize{
#' \item{Rtsne: }{Use the Rtsne package Barnes-Hut implementation of tSNE (default)}
# \item{tsne: }{standard tsne - not recommended for large datasets}
#' \item{FIt-SNE: }{Use the FFT-accelerated Interpolation-based t-SNE. Based on
#' Kluger Lab code found here: https://github.com/KlugerLab/FIt-SNE}
#' }
#' @param dim.embed The dimensional space of the resulting tSNE embedding
#' (default is 2). For example, set to 3 for a 3d tSNE
#' @param reduction.key dimensional reduction key, specifies the string before the number for the dimension names. tSNE_ by default
#'
#' @importFrom Rtsne Rtsne
#'
#' @rdname RunTSNE
#' @concept dimensional_reduction
#' @export
#' @method RunTSNE matrix
#'
RunTSNE.matrix <- function(
  object,
  assay = NULL,
  seed.use = 1,
  tsne.method = "Rtsne",
  dim.embed = 2,
  reduction.key = "tSNE_",
  ...
) {
  if (!is.null(x = seed.use)) {
    set.seed(seed = seed.use)
  }
  tsne.data <- switch(
    EXPR = tsne.method,
    'Rtsne' = Rtsne(
      X = object,
      dims = dim.embed,
      pca = FALSE,
      ... # PCA/is_distance
    )$Y,
    'FIt-SNE' = fftRtsne(X = object, dims = dim.embed, rand_seed = seed.use, ...),
    stop("Invalid tSNE method: please choose from 'Rtsne' or 'FIt-SNE'")
  )
  colnames(x = tsne.data) <- paste0(reduction.key, 1:ncol(x = tsne.data))
  rownames(x = tsne.data) <- rownames(x = object)
  tsne.reduction <- CreateDimReducObject(
    embeddings = tsne.data,
    key = reduction.key,
    assay = assay,
    global = TRUE
  )
  return(tsne.reduction)
}

#' @param cells Which cells to analyze (default, all cells)
#' @param dims Which dimensions to use as input features
#'
#' @rdname RunTSNE
#' @concept dimensional_reduction
#' @export
#' @method RunTSNE DimReduc
#'
RunTSNE.DimReduc <- function(
  object,
  cells = NULL,
  dims = 1:5,
  seed.use = 1,
  tsne.method = "Rtsne",
  dim.embed = 2,
  reduction.key = "tSNE_",
  ...
) {
  args <- as.list(x = sys.frame(which = sys.nframe()))
  args <- c(args, list(...))
  args$object <- args$object[[cells, args$dims]]
  args$dims <- NULL
  args$cells <- NULL
  args$assay <- DefaultAssay(object = object)
  return(do.call(what = 'RunTSNE', args = args))
}

#' @rdname RunTSNE
#' @concept dimensional_reduction
#' @export
#' @method RunTSNE dist
#'
RunTSNE.dist <- function(
  object,
  assay = NULL,
  seed.use = 1,
  tsne.method = "Rtsne",
  dim.embed = 2,
  reduction.key = "tSNE_",
  ...
) {
  args <- as.list(x = sys.frame(which = sys.nframe()))
  args <- c(args, list(...))
  args$object <- as.matrix(x = args$object)
  args$is_distance <- TRUE
  return(do.call(what = 'RunTSNE', args = args))
}

#' @param reduction Which dimensional reduction (e.g. PCA, ICA) to use for
#' the tSNE. Default is PCA
#' @param features If set, run the tSNE on this subset of features
#' (instead of running on a set of reduced dimensions). Not set (NULL) by default;
#' \code{dims} must be NULL to run on features
#' @param distance.matrix If set, runs tSNE on the given distance matrix
#' instead of data matrix (experimental)
#' @param reduction.name dimensional reduction name, specifies the position in the object$dr list. tsne by default
#'
#' @rdname RunTSNE
#' @concept dimensional_reduction
#' @export
#' @method RunTSNE Seurat
#'
RunTSNE.Seurat <- function(
  object,
  reduction = "pca",
  cells = NULL,
  dims = 1:5,
  features = NULL,
  seed.use = 1,
  tsne.method = "Rtsne",
  dim.embed = 2,
  distance.matrix = NULL,
  reduction.name = "tsne",
  reduction.key = "tSNE_",
  ...
) {
  cells <- cells %||% Cells(x = object)
  tsne.reduction <- if (!is.null(x = distance.matrix)) {
    RunTSNE(
      object = distance.matrix,
      assay = DefaultAssay(object = object),
      seed.use = seed.use,
      tsne.method = tsne.method,
      dim.embed = dim.embed,
      reduction.key = reduction.key,
      is_distance = TRUE,
      ...
    )
  } else if (!is.null(x = dims)) {
    RunTSNE(
      object = object[[reduction]],
      cells = cells,
      dims = dims,
      seed.use = seed.use,
      tsne.method = tsne.method,
      dim.embed = dim.embed,
      reduction.key = reduction.key,
      ...
    )
  } else if (!is.null(x = features)) {
    RunTSNE(
      object = t(x = as.matrix(x = GetAssayData(object = object)[features, cells])),
      assay = DefaultAssay(object = object),
      seed.use = seed.use,
      tsne.method = tsne.method,
      dim.embed = dim.embed,
      reduction.key = reduction.key,
      ...
    )
  } else {
    stop("Unknown way of running tSNE")
  }
  object[[reduction.name]] <- tsne.reduction
  object <- LogSeuratCommand(object = object)
  return(object)
}

#' @importFrom reticulate py_module_available py_set_seed import
#' @importFrom uwot umap umap_transform
#' @importFrom future nbrOfWorkers
#'
#' @rdname RunUMAP
#' @concept dimensional_reduction
#' @method RunUMAP default
#' @export
#'
RunUMAP.default <- function(
  object,
  reduction.key = 'UMAP_',
  assay = NULL,
  reduction.model = NULL,
  return.model = FALSE,
  umap.method = 'uwot',
  n.neighbors = 30L,
  n.components = 2L,
  metric = 'cosine',
  n.epochs = NULL,
  learning.rate = 1.0,
  min.dist = 0.3,
  spread = 1.0,
  set.op.mix.ratio = 1.0,
  local.connectivity = 1L,
  repulsion.strength = 1,
  negative.sample.rate = 5,
  a = NULL,
  b = NULL,
  uwot.sgd = FALSE,
  seed.use = 42,
  metric.kwds = NULL,
  angular.rp.forest = FALSE,
  densmap = FALSE,
  dens.lambda = 2,
  dens.frac = 0.3,
  dens.var.shift = 0.1,
  verbose = TRUE,
  ...
) {
  CheckDots(...)
  if (!is.null(x = seed.use)) {
    set.seed(seed = seed.use)
  }
  if (umap.method != 'umap-learn' && getOption('Seurat.warn.umap.uwot', TRUE)) {
    warning(
      "The default method for RunUMAP has changed from calling Python UMAP via reticulate to the R-native UWOT using the cosine metric",
      "\nTo use Python UMAP via reticulate, set umap.method to 'umap-learn' and metric to 'correlation'",
      "\nThis message will be shown once per session",
      call. = FALSE,
      immediate. = TRUE
    )
    options(Seurat.warn.umap.uwot = FALSE)
  }
  if (umap.method == 'uwot-learn') {
    warning("'uwot-learn' is deprecated. Set umap.method = 'uwot' and return.model = TRUE")
    umap.method <- "uwot"
    return.model <- TRUE
  }
  if (densmap && umap.method != 'umap-learn'){
    warning("densmap is only supported by umap-learn method. Method is changed to 'umap-learn'")
    umap.method <- 'umap-learn'
  }
  if (return.model) {
    if (verbose) {
      message("UMAP will return its model")
    }
    umap.method = "uwot"
  }
  if (inherits(x = object, what = "Neighbor")) {
    object <- list( idx = Indices(object),
                    dist = Distances(object) )
  }
  if (!is.null(x = reduction.model)) {
    if (verbose) {
      message("Running UMAP projection")
    }
    umap.method <- "uwot-predict"
  }
  umap.output <- switch(
    EXPR = umap.method,
    'umap-learn' = {
      if (!py_module_available(module = 'umap')) {
        stop("Cannot find UMAP, please install through pip (e.g. pip install umap-learn).")
      }
      if (!py_module_available(module = 'sklearn')) {
        stop("Cannot find sklearn, please install through pip (e.g. pip install scikit-learn).")
      }
      if (!is.null(x = seed.use)) {
        py_set_seed(seed = seed.use)
      }
      if (typeof(x = n.epochs) == "double") {
        n.epochs <- as.integer(x = n.epochs)
      }
      umap_import <- import(module = "umap", delay_load = TRUE)
      sklearn <- import("sklearn", delay_load = TRUE)
      if (densmap &&
          numeric_version(x = umap_import$pkg_resources$get_distribution("umap-learn")$version) <
          numeric_version(x = "0.5.0")) {
        stop("densmap is only supported by versions >= 0.5.0 of umap-learn. Upgrade umap-learn (e.g. pip install --upgrade umap-learn).")
      }
      random.state <- sklearn$utils$check_random_state(seed = as.integer(x = seed.use))
      umap.args <- list(
        n_neighbors = as.integer(x = n.neighbors),
        n_components = as.integer(x = n.components),
        metric = metric,
        n_epochs = n.epochs,
        learning_rate = learning.rate,
        min_dist = min.dist,
        spread = spread,
        set_op_mix_ratio = set.op.mix.ratio,
        local_connectivity = local.connectivity,
        repulsion_strength = repulsion.strength,
        negative_sample_rate = negative.sample.rate,
        random_state = random.state,
        a = a,
        b = b,
        metric_kwds = metric.kwds,
        angular_rp_forest = angular.rp.forest,
        verbose = verbose
      )
      if (numeric_version(x = umap_import$pkg_resources$get_distribution("umap-learn")$version) >=
          numeric_version(x = "0.5.0")) {
        umap.args <- c(umap.args, list(
          densmap = densmap,
          dens_lambda = dens.lambda,
          dens_frac = dens.frac,
          dens_var_shift = dens.var.shift,
          output_dens = FALSE
        ))
      }
      umap <- do.call(what = umap_import$UMAP, args = umap.args)
      umap$fit_transform(as.matrix(x = object))
    },
    'uwot' = {
      if (is.list(x = object)) {
        umap(
          X = NULL,
          nn_method = object,
          n_threads = nbrOfWorkers(),
          n_components = as.integer(x = n.components),
          metric = metric,
          n_epochs = n.epochs,
          learning_rate = learning.rate,
          min_dist = min.dist,
          spread = spread,
          set_op_mix_ratio = set.op.mix.ratio,
          local_connectivity = local.connectivity,
          repulsion_strength = repulsion.strength,
          negative_sample_rate = negative.sample.rate,
          a = a,
          b = b,
          fast_sgd = uwot.sgd,
          verbose = verbose,
          ret_model = return.model
        )
      } else {
        umap(
          X = object,
          n_threads = nbrOfWorkers(),
          n_neighbors = as.integer(x = n.neighbors),
          n_components = as.integer(x = n.components),
          metric = metric,
          n_epochs = n.epochs,
          learning_rate = learning.rate,
          min_dist = min.dist,
          spread = spread,
          set_op_mix_ratio = set.op.mix.ratio,
          local_connectivity = local.connectivity,
          repulsion_strength = repulsion.strength,
          negative_sample_rate = negative.sample.rate,
          a = a,
          b = b,
          fast_sgd = uwot.sgd,
          verbose = verbose,
          ret_model = return.model
        )
      }
    },
    'uwot-predict' = {
      if (metric == 'correlation') {
        warning(
          "UWOT does not implement the correlation metric, using cosine instead",
          call. = FALSE,
          immediate. = TRUE
        )
        metric <- 'cosine'
      }
      if (is.null(x = reduction.model) || !inherits(x = reduction.model, what = 'DimReduc')) {
        stop(
          "If running projection UMAP, please pass a DimReduc object with the model stored to reduction.model.",
          call. = FALSE
        )
      }
      model <- Misc(
        object = reduction.model,
        slot = "model"
      )
      # add num_precomputed_nns to <v0.1.13 uwot models to prevent errors with newer versions of uwot
      if (!"num_precomputed_nns" %in% names(model)) {
        model$num_precomputed_nns <- 1
      }
      if (length(x = model) == 0) {
        stop(
          "The provided reduction.model does not have a model stored. Please try running umot-learn on the object first",
          call. = FALSE
        )
      }
      if (!"num_precomputed_nns" %in% names(x = model)) {
        model$num_precomputed_nns <- 0
      }
      if (is.list(x = object)) {
        if (ncol(object$idx) != model$n_neighbors) {
          warning("Number of neighbors between query and reference ",
          "is not equal to the number of neighbors within reference")
          model$n_neighbors <- ncol(object$idx)
        }
       umap_transform(
          X = NULL,
          nn_method = object,
          model = model,
          n_threads = nbrOfWorkers(),
          n_epochs = n.epochs,
          verbose = verbose
        )
      } else {
        umap_transform(
          X = object,
          model = model,
          n_threads = nbrOfWorkers(),
          n_epochs = n.epochs,
          verbose = verbose
        )
      }
    },
    stop("Unknown umap method: ", umap.method, call. = FALSE)
  )
  if (return.model) {
    umap.output$nn_index <- NULL
    umap.model <- umap.output
    umap.output <- umap.output$embedding
  }
  colnames(x = umap.output) <- paste0(reduction.key, 1:ncol(x = umap.output))
  if (inherits(x = object, what = 'dist')) {
    rownames(x = umap.output) <- attr(x = object, "Labels")
  } else if (is.list(x = object)) {
    rownames(x = umap.output) <- rownames(x = object$idx)
  } else {
    rownames(x = umap.output) <- rownames(x = object)
  }
  umap.reduction <- CreateDimReducObject(
    embeddings = umap.output,
    key = reduction.key,
    assay = assay,
    global = TRUE
  )
  if (return.model) {
    Misc(umap.reduction, slot = "model") <- umap.model
  }
  return(umap.reduction)
}

#' @importFrom reticulate py_module_available import
#'
#' @rdname RunUMAP
#' @concept dimensional_reduction
#' @method RunUMAP Graph
#' @export
#'
RunUMAP.Graph <- function(
  object,
  assay = NULL,
  umap.method = 'umap-learn',
  n.components = 2L,
  metric = 'correlation',
  n.epochs = 0L,
  learning.rate = 1,
  min.dist = 0.3,
  spread = 1,
  repulsion.strength = 1,
  negative.sample.rate = 5L,
  a = NULL,
  b = NULL,
  uwot.sgd = FALSE,
  seed.use = 42L,
  metric.kwds = NULL,
  densmap = FALSE,
  densmap.kwds = NULL,
  verbose = TRUE,
  reduction.key = 'UMAP_',
  ...
) {
  #CheckDots(...)
  if (umap.method != 'umap-learn') {
    warning(
      "Running UMAP on Graph objects is only supported using the umap-learn method",
      call. = FALSE,
      immediate. = TRUE
    )
  }
  if (!py_module_available(module = 'umap')) {
    stop("Cannot find UMAP, please install through pip (e.g. pip install umap-learn).")
  }
  if (!py_module_available(module = 'numpy')) {
    stop("Cannot find numpy, please install through pip (e.g. pip install numpy).")
  }
  if (!py_module_available(module = 'sklearn')) {
    stop("Cannot find sklearn, please install through pip (e.g. pip install scikit-learn).")
  }
  if (!py_module_available(module = 'scipy')) {
    stop("Cannot find scipy, please install through pip (e.g. pip install scipy).")
  }
  np <- import("numpy", delay_load = TRUE)
  sp <- import("scipy", delay_load = TRUE)
  sklearn <- import("sklearn", delay_load = TRUE)
  umap <- import("umap", delay_load = TRUE)
  diag(x = object) <- 0
  data <- object
  object <- sp$sparse$coo_matrix(arg1 = object)
  ab.params <- umap$umap_$find_ab_params(spread = spread, min_dist = min.dist)
  a <- a %||% ab.params[[1]]
  b <- b %||% ab.params[[2]]
  n.epochs <- n.epochs %||% 0L
  random.state <- sklearn$utils$check_random_state(seed = as.integer(x = seed.use))
  umap.args <- list(
    data = data,
    graph = object,
    n_components = n.components,
    initial_alpha = learning.rate,
    a = a,
    b = b,
    gamma = repulsion.strength,
    negative_sample_rate = negative.sample.rate,
    n_epochs = as.integer(x = n.epochs),
    random_state = random.state,
    init = "spectral",
    metric = metric,
    metric_kwds = metric.kwds,
    verbose = verbose
  )
  if (numeric_version(x = umap$pkg_resources$get_distribution("umap-learn")$version) >=
      numeric_version(x = "0.5.0")) {
    umap.args <- c(umap.args, list(
      densmap = densmap,
      densmap_kwds = densmap.kwds,
      output_dens = FALSE
    ))
  }
  embeddings <- do.call(what = umap$umap_$simplicial_set_embedding, args = umap.args)
  if (length(x = embeddings) == 2) {
    embeddings <- embeddings[[1]]
  }
  rownames(x = embeddings) <- colnames(x = data)
  colnames(x = embeddings) <- paste0("UMAP_", 1:n.components)
  # center the embeddings on zero
  embeddings <- scale(x = embeddings, scale = FALSE)
  umap <- CreateDimReducObject(
    embeddings = embeddings,
    key = reduction.key,
    assay = assay,
    global = TRUE
  )
  return(umap)
}

#' @rdname RunUMAP
#' @concept dimensional_reduction
#' @method RunUMAP Neighbor
#' @export
#'
RunUMAP.Neighbor <- function(
  object,
  reduction.model,
  ...
) {
  neighborlist <- list("idx" = Indices(object),
                       "dist" = Distances(object))
  RunUMAP(
    object = neighborlist,
    reduction.model = reduction.model,
    ...
  )
}


#' @param reduction.model \code{DimReduc} object that contains the umap model
#' @param dims Which dimensions to use as input features, used only if
#' \code{features} is NULL
#' @param reduction Which dimensional reduction (PCA or ICA) to use for the
#' UMAP input. Default is PCA
#' @param features If set, run UMAP on this subset of features (instead of running on a
#' set of reduced dimensions). Not set (NULL) by default; \code{dims} must be NULL to run
#' on features
#' @param graph Name of graph on which to run UMAP
#' @param assay Assay to pull data for when using \code{features}, or assay used to construct Graph
#' if running UMAP on a Graph
#' @param nn.name Name of knn output on which to run UMAP
#' @param slot The slot used to pull data for when using \code{features}. data slot is by default.
#' @param umap.method UMAP implementation to run. Can be
#' \describe{
#'   \item{\code{uwot}:}{Runs umap via the uwot R package}
#'   \item{\code{uwot-learn}:}{Runs umap via the uwot R package and return the learned umap model}
#'   \item{\code{umap-learn}:}{Run the Seurat wrapper of the python umap-learn package}
#' }
#' @param n.neighbors This determines the number of neighboring points used in
#' local approximations of manifold structure. Larger values will result in more
#' global structure being preserved at the loss of detailed local structure. In
#' general this parameter should often be in the range 5 to 50.
#' @param n.components The dimension of the space to embed into.
#' @param metric metric: This determines the choice of metric used to measure
#' distance in the input space. A wide variety of metrics are already coded, and
#' a user defined function can be passed as long as it has been JITd by numba.
#' @param n.epochs he number of training epochs to be used in optimizing the low dimensional
#' embedding. Larger values result in more accurate embeddings. If NULL is specified, a value will
#' be selected based on the size of the input dataset (200 for large datasets, 500 for small).
#' @param learning.rate The initial learning rate for the embedding optimization.
#' @param min.dist This controls how tightly the embedding is allowed compress points together.
#' Larger values ensure embedded points are moreevenly distributed, while smaller values allow the
#' algorithm to optimise more accurately with regard to local structure. Sensible values are in
#' the range 0.001 to 0.5.
#' @param spread The effective scale of embedded points. In combination with min.dist this
#' determines how clustered/clumped the embedded points are.
#' @param set.op.mix.ratio Interpolate between (fuzzy) union and intersection as the set operation
#' used to combine local fuzzy simplicial sets to obtain a global fuzzy simplicial sets. Both fuzzy
#' set operations use the product t-norm. The value of this parameter should be between 0.0 and
#' 1.0; a value of 1.0 will use a pure fuzzy union, while 0.0 will use a pure fuzzy intersection.
#' @param local.connectivity The local connectivity required - i.e. the number of nearest neighbors
#' that should be assumed to be connected at a local level. The higher this value the more connected
#' the manifold becomes locally. In practice this should be not more than the local intrinsic
#' dimension of the manifold.
#' @param repulsion.strength Weighting applied to negative samples in low dimensional embedding
#' optimization. Values higher than one will result in greater weight being given to negative
#' samples.
#' @param negative.sample.rate The number of negative samples to select per positive sample in the
#' optimization process. Increasing this value will result in greater repulsive force being applied,
#' greater optimization cost, but slightly more accuracy.
#' @param a More specific parameters controlling the embedding. If NULL, these values are set
#' automatically as determined by min. dist and spread. Parameter of differentiable approximation of
#' right adjoint functor.
#' @param b More specific parameters controlling the embedding. If NULL, these values are set
#' automatically as determined by min. dist and spread. Parameter of differentiable approximation of
#' right adjoint functor.
#' @param uwot.sgd Set \code{uwot::umap(fast_sgd = TRUE)}; see \code{\link[uwot]{umap}} for more details
#' @param metric.kwds A dictionary of arguments to pass on to the metric, such as the p value for
#' Minkowski distance. If NULL then no arguments are passed on.
#' @param angular.rp.forest Whether to use an angular random projection forest to initialise the
#' approximate nearest neighbor search. This can be faster, but is mostly on useful for metric that
#' use an angular style distance such as cosine, correlation etc. In the case of those metrics
#' angular forests will be chosen automatically.
#' @param densmap Whether to use the density-augmented objective of densMAP.
#' Turning on this option generates an embedding where the local densities
#' are encouraged to be correlated with those in the original space.
#' Parameters below with the prefix ‘dens’ further control the behavior
#' of this extension. Default is FALSE. Only compatible with 'umap-learn' method
#' and version of umap-learn >= 0.5.0
#' @param densmap.kwds A dictionary of arguments to pass on to the densMAP optimization.
#' @param dens.lambda Specific parameter which controls the regularization weight
#' of the density correlation term in densMAP. Higher values prioritize density
#' preservation over the UMAP objective, and vice versa for values closer to zero.
#' Setting this parameter to zero is equivalent to running the original UMAP algorithm.
#' Default value is 2.
#' @param dens.frac Specific parameter which controls the fraction of epochs
#' (between 0 and 1) where the density-augmented objective is used in densMAP.
#' The first (1 - dens_frac) fraction of epochs optimize the original UMAP
#' objective before introducing the density correlation term. Default is 0.3.
#' @param dens.var.shift Specific parameter which specifies a small constant
#' added to the variance of local radii in the embedding when calculating
#' the density correlation objective to prevent numerical instability from
#' dividing by a small number. Default is 0.1.
#' @param reduction.name Name to store dimensional reduction under in the Seurat object
#' @param reduction.key dimensional reduction key, specifies the string before
#' the number for the dimension names. UMAP by default
#' @param return.model whether UMAP will return the uwot model
#' @param seed.use Set a random seed. By default, sets the seed to 42. Setting
#' NULL will not set a seed
#' @param verbose Controls verbosity
#'
#' @rdname RunUMAP
#' @concept dimensional_reduction
#' @export
#' @method RunUMAP Seurat
#'
RunUMAP.Seurat <- function(
  object,
  dims = NULL,
  reduction = 'pca',
  features = NULL,
  graph = NULL,
  assay = DefaultAssay(object = object),
  nn.name = NULL,
  slot = 'data',
  umap.method = 'uwot',
  reduction.model = NULL,
  return.model = FALSE,
  n.neighbors = 30L,
  n.components = 2L,
  metric = 'cosine',
  n.epochs = NULL,
  learning.rate = 1,
  min.dist = 0.3,
  spread = 1,
  set.op.mix.ratio = 1,
  local.connectivity = 1L,
  repulsion.strength = 1,
  negative.sample.rate = 5L,
  a = NULL,
  b = NULL,
  uwot.sgd = FALSE,
  seed.use = 42L,
  metric.kwds = NULL,
  angular.rp.forest = FALSE,
  densmap = FALSE,
  dens.lambda = 2,
  dens.frac = 0.3,
  dens.var.shift = 0.1,
  verbose = TRUE,
  reduction.name = 'umap',
  reduction.key = NULL,
  ...
) {
  CheckDots(...)
  if (sum(c(is.null(x = dims), is.null(x = features), is.null(x = graph))) < 2) {
      stop("Please specify only one of the following arguments: dims, features, or graph")
  }
  if (!is.null(x = features)) {
    data.use <- as.matrix(x = t(x = GetAssayData(object = object, slot = slot, assay = assay)[features, , drop = FALSE]))
    if (ncol(x = data.use) < n.components) {
      stop(
        "Please provide as many or more features than n.components: ",
        length(x = features),
        " features provided, ",
        n.components,
        " UMAP components requested",
        call. = FALSE
      )
    }
  } else if (!is.null(x = dims)) {
    data.use <- Embeddings(object[[reduction]])[, dims]
    assay <- DefaultAssay(object = object[[reduction]])
    if (length(x = dims) < n.components) {
      stop(
        "Please provide as many or more dims than n.components: ",
        length(x = dims),
        " dims provided, ",
        n.components,
        " UMAP components requested",
        call. = FALSE
      )
    }
  }  else if (!is.null(x = nn.name)) {
    if (!inherits(x = object[[nn.name]], what = "Neighbor")) {
      stop(
        "Please specify a Neighbor object name, ",
        "instead of the name of a ",
        class(object[[nn.name]]),
        " object",
        call. = FALSE
      )
    }
    data.use <- object[[nn.name]]
  } else if (!is.null(x = graph)) {
    if (!inherits(x = object[[graph]], what = "Graph")) {
      stop(
        "Please specify a Graph object name, ",
        "instead of the name of a ",
        class(object[[graph]]),
        " object",
        call. = FALSE
      )
    }
    data.use <- object[[graph]]
  } else {
    stop("Please specify one of dims, features, or graph")
  }
  object[[reduction.name]] <- RunUMAP(
    object = data.use,
    reduction.model = reduction.model,
    return.model = return.model,
    assay = assay,
    umap.method = umap.method,
    n.neighbors = n.neighbors,
    n.components = n.components,
    metric = metric,
    n.epochs = n.epochs,
    learning.rate = learning.rate,
    min.dist = min.dist,
    spread = spread,
    set.op.mix.ratio = set.op.mix.ratio,
    local.connectivity = local.connectivity,
    repulsion.strength = repulsion.strength,
    negative.sample.rate = negative.sample.rate,
    a = a,
    b = b,
    uwot.sgd = uwot.sgd,
    seed.use = seed.use,
    metric.kwds = metric.kwds,
    angular.rp.forest = angular.rp.forest,
    densmap = densmap,
    dens.lambda = dens.lambda,
    dens.frac = dens.frac,
    dens.var.shift = dens.var.shift,
    reduction.key = reduction.key %||% Key(object = reduction.name, quiet = TRUE),
    verbose = verbose
  )
  object <- LogSeuratCommand(object = object)
  return(object)
}

#' @param dims Which dimensions to examine
#' @param score.thresh Threshold to use for the proportion test of PC
#' significance (see Details)
#'
#' @importFrom stats prop.test
#'
#' @rdname ScoreJackStraw
#' @concept dimensional_reduction
#' @export
#' @method ScoreJackStraw JackStrawData
#'
ScoreJackStraw.JackStrawData <- function(
  object,
  dims = 1:5,
  score.thresh = 1e-5,
  ...
) {
  CheckDots(...)
  pAll <- JS(object = object, slot = "empirical.p.values")
  pAll <- pAll[, dims, drop = FALSE]
  pAll <- as.data.frame(pAll)
  pAll$Contig <- rownames(x = pAll)
  score.df <- NULL
  for (i in dims) {
    pc.score <- suppressWarnings(prop.test(
      x = c(
        length(x = which(x = pAll[, i] <= score.thresh)),
        floor(x = nrow(x = pAll) * score.thresh)
      ),
      n = c(nrow(pAll), nrow(pAll))
    )$p.val)
    if (length(x = which(x = pAll[, i] <= score.thresh)) == 0) {
      pc.score <- 1
    }
    if (is.null(x = score.df)) {
      score.df <- data.frame(PC = paste0("PC", i), Score = pc.score)
    } else {
      score.df <- rbind(score.df, data.frame(PC = paste0("PC", i), Score = pc.score))
    }
  }
  score.df$PC <- dims
  score.df <- as.matrix(score.df)
  JS(object = object, slot = 'overall') <- score.df
  return(object)
}

#' @rdname ScoreJackStraw
#' @concept dimensional_reduction
#' @export
#' @method ScoreJackStraw DimReduc
#'
ScoreJackStraw.DimReduc <- function(object, dims = 1:5, score.thresh = 1e-5, ...) {
  JS(object = object) <- ScoreJackStraw(
    object = JS(object = object),
    dims = dims,
    score.thresh = score.thresh,
    ...
  )
  return(object)
}

#' @param reduction Reduction associated with JackStraw to score
#' @param do.plot Show plot. To return ggplot object, use \code{JackStrawPlot} after
#' running ScoreJackStraw.
#'
#' @seealso \code{\link{JackStrawPlot}}
#'
#' @rdname ScoreJackStraw
#' @concept dimensional_reduction
#' @export
#' @method ScoreJackStraw Seurat
#'
ScoreJackStraw.Seurat <- function(
  object,
  reduction = "pca",
  dims = 1:5,
  score.thresh = 1e-5,
  do.plot = FALSE,
  ...
) {
  object[[reduction]] <- ScoreJackStraw(
    object = object[[reduction]],
    dims = dims,
    score.thresh = score.thresh,
    ...
  )
  if (do.plot) {
    CheckDots(..., fxns = 'JackStrawPlot')
    suppressWarnings(expr = print(JackStrawPlot(
      object = object,
      reduction = reduction,
      dims = dims,
      ...
    )))
  }
  object <- LogSeuratCommand(object = object)
  return(object)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Methods for R-defined generics
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Internal
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Check that features are present and have non-zero variance
#
# @param data.use Feature matrix (features are rows)
# @param features Features to check
# @param object.name Name of object for message printing
# @param verbose Print warnings
#
# @return Returns a vector of features that is the subset of features
# that have non-zero variance
#
CheckFeatures <- function(
  data.use,
  features,
  object.name,
  verbose = TRUE
) {
  if (any(!features %in% rownames(x = data.use))) {
    missing.features <- features[!features %in% rownames(x = data.use)]
    features <- setdiff(x = features, y = missing.features)
    if (verbose){
      warning(
        paste0(
          "The following ", length(x = missing.features),
           " features are not scaled in ",
           object.name,
           ": ",
           paste0(missing.features, collapse = ", ")
        ))
    }
  }
  if (inherits(x = data.use, what = 'dgCMatrix')) {
    features.var <- SparseRowVar(mat = data.use[features, ], display_progress = F)
  }
  else {
    features.var <- RowVar(x = data.use[features, ])
  }
  no.var.features <- features[features.var == 0]
  if (length(x = no.var.features) > 0 && verbose) {
    warning(
     paste0(
       "The following features have zero variance in ",
       object.name,
       ": ",
       paste0(no.var.features, collapse = ", ")
    ))
  }
  features <- setdiff(x = features, y = no.var.features)
  features <- features[!is.na(x = features)]
  return(features)
}

#internal
EmpiricalP <- function(x, nullval) {
  return(sum(nullval > x) / length(x = nullval))
}

# FIt-SNE helper function for calling fast_tsne from R
#
# Based on Kluger Lab FIt-SNE v1.2.1 code on https://github.com/KlugerLab/FIt-SNE/blob/master/fast_tsne.R
# commit 601608ed42e4be2765970910927da20f0b0bf9b9 on June 25, 2020
#
#' @importFrom utils file_test
#
fftRtsne <- function(X,
  dims = 2,
  perplexity = 30,
  theta = 0.5,
  max_iter = 750,
  fft_not_bh = TRUE,
  ann_not_vptree = TRUE,
  stop_early_exag_iter = 250,
  exaggeration_factor = 12.0,
  no_momentum_during_exag = FALSE,
  start_late_exag_iter = -1,
  late_exag_coeff = 1.0,
  mom_switch_iter = 250,
  momentum = 0.5,
  final_momentum = 0.8,
  learning_rate = 'auto',
  n_trees = 50,
  search_k = -1,
  rand_seed = -1,
  nterms = 3,
  intervals_per_integer = 1,
  min_num_intervals = 50,
  K = -1,
  sigma = -30,
  initialization = 'pca',
  max_step_norm = 5,
  data_path = NULL,
  result_path = NULL,
  load_affinities = NULL,
  fast_tsne_path = NULL,
  nthreads = getOption('mc.cores', default = 1),
  perplexity_list = NULL,
  get_costs = FALSE,
  df = 1.0,
  ...
) {
  CheckDots(...)
  if (is.null(x = data_path)) {
    data_path <- tempfile(pattern = 'fftRtsne_data_', fileext = '.dat')
  }
  if (is.null(x = result_path)) {
    result_path <- tempfile(pattern = 'fftRtsne_result_', fileext = '.dat')
  }
  if (is.null(x = fast_tsne_path)) {
    # suppressWarnings(expr = fast_tsne_path <- system2(command = 'which', args = 'fast_tsne', stdout = TRUE))
    fast_tsne_path <- SysExec(progs = ifelse(
      test = .Platform$OS.type == 'windows',
      yes = 'FItSNE.exe',
      no = 'fast_tsne'
      ))
    if (length(x = fast_tsne_path) == 0) {
      stop("no fast_tsne_path specified and fast_tsne binary is not in the search path")
    }
  }
  fast_tsne_path <- normalizePath(path = fast_tsne_path)
  if (!file_test(op = '-x', x = fast_tsne_path)) {
    stop("fast_tsne_path '", fast_tsne_path, "' does not exist or is not executable")
  }
  # check fast_tsne version
  ft.out <- suppressWarnings(expr = system2(command = fast_tsne_path, stdout = TRUE))
  version_number <- regmatches(ft.out[1], regexpr('= t-SNE v[0-9.]+', ft.out[1]))
  if (is.null(version_number)){
    message("First line of fast_tsne output is")
    message(ft.out[1])
    stop("Our FIt-SNE wrapper requires FIt-SNE v1.0+, please install the appropriate version from github.com/KlugerLab/FIt-SNE and have fast_tsne_path point to it if it's not in your path")
  } else {
    version_number <- gsub('= t-SNE v', '', version_number)
  }

  is.wholenumber <- function(x, tol = .Machine$double.eps ^ 0.5) {
    return(abs(x = x - round(x = x)) < tol)
  }
  if (version_number == '1.0.0' && df != 1.0) {
    stop("This version of FIt-SNE does not support df!=1. Please install the appropriate version from github.com/KlugerLab/FIt-SNE")
  }
  if (!is.numeric(x = theta) || (theta < 0.0) || (theta > 1.0) ) {
    stop("Incorrect theta.")
  }
  if (nrow(x = X) - 1 < 3 * perplexity) {
    stop("Perplexity is too large.")
  }
  if (!is.matrix(x = X)) {
    stop("Input X is not a matrix")
  }
  if (!(max_iter > 0)) {
    stop("Incorrect number of iterations.")
  }
  if (!is.wholenumber(x = stop_early_exag_iter) || stop_early_exag_iter < 0) {
    stop("stop_early_exag_iter should be a positive integer")
  }
  if (!is.numeric(x = exaggeration_factor)) {
    stop("exaggeration_factor should be numeric")
  }
  if (!is.numeric(df)) {
    stop("df should be numeric")
  }
  if (!is.wholenumber(x = dims) || dims <= 0) {
    stop("Incorrect dimensionality.")
  }
  if (search_k == -1) {
    if (perplexity > 0) {
      search_k <- n_trees * perplexity * 3
    } else if (perplexity == 0) {
      search_k <- n_trees * max(perplexity_list) * 3
    } else {
      search_k <- n_trees * K
    }
  }

  if (is.character(learning_rate) && learning_rate =='auto') {
    learning_rate = max(200, nrow(X)/exaggeration_factor)
  }
  if (is.character(start_late_exag_iter) && start_late_exag_iter =='auto') {
    if (late_exag_coeff > 0) {
      start_late_exag_iter = stop_early_exag_iter
    } else {
      start_late_exag_iter = -1
    }
  }

  if (is.character(initialization) && initialization == 'pca') {
    if (rand_seed != -1)  {
      set.seed(rand_seed)
    }
    if (requireNamespace("rsvd", quietly = TRUE)) {
      message('Using rsvd() to compute the top PCs for initialization.')
      X_c <- scale(x = X, center = TRUE, scale = FALSE)
      rsvd_out <- rsvd::rsvd(A = X_c, k = dims)
      X_top_pcs <- rsvd_out$u %*% diag(x = rsvd_out$d, nrow = dims)
    } else if (requireNamespace("irlba", quietly = TRUE)) {
      message('Using irlba() to compute the top PCs for initialization.')
      X_colmeans <- colMeans(x = X)
      irlba_out <- irlba::irlba(A = X, nv = dims, center = X_colmeans)
      X_top_pcs <- irlba_out$u %*% diag(x = irlba_out$d, nrow = dims)
    } else {
      stop(
        "By default, FIt-SNE initializes the embedding with the
        top PCs. We use either rsvd or irlba for fast computation.
        To use this functionality, please install the rsvd package
        with install.packages('rsvd') or the irlba package with
        install.packages('ilrba').  Otherwise, set initialization
        to NULL for random initialization, or any N by dims matrix
        for custom initialization."
      )
    }
    initialization <- 0.0001*(X_top_pcs/sd(X_top_pcs[,1]))

  } else if (is.character(x = initialization) && initialization == 'random') {
    message('Random initialization')
    initialization = NULL
  }
  nbody_algo <- ifelse(test = fft_not_bh, yes = 2, no = 1)

  if (is.null(load_affinities)) {
    load_affinities <- 0
  } else {
    if (load_affinities == 'load') {
      load_affinities <- 1
    } else if (load_affinities == 'save') {
      load_affinities <- 2
    } else {
      load_affinities <- 0
    }
  }

  knn_algo <- ifelse(test = ann_not_vptree, yes = 1, no = 2)
  tX <- as.numeric(t(X))

  f <- file(description = data_path, open = "wb")
  n = nrow(x = X)
  D = ncol(x = X)
  writeBin(object = as.integer(x = n), con = f, size = 4)
  writeBin(object = as.integer(x = D), con = f, size = 4)
  writeBin(object = as.numeric(x = theta), con = f, size = 8)
  writeBin(object = as.numeric(x = perplexity), con = f, size = 8)

  if (perplexity == 0) {
    writeBin(object = as.integer(x = length(x = perplexity_list)), con = f, size = 4)
    writeBin(object = perplexity_list, con = f)
  }

  writeBin(object = as.integer(x = dims), con = f, size = 4) #theta
  writeBin(object = as.integer(x = max_iter), con = f, size = 4)
  writeBin(object = as.integer(x = stop_early_exag_iter), con = f, size = 4)
  writeBin(object = as.integer(x = mom_switch_iter), con = f, size = 4)
  writeBin(object = as.numeric(x = momentum), con = f, size = 8)
  writeBin(object = as.numeric(x = final_momentum), con = f, size = 8)
  writeBin(object = as.numeric(x = learning_rate), con = f, size = 8)
  if (!(version_number %in% c('1.1.0', '1.0.0'))) {
    writeBin(object = as.numeric(x = max_step_norm), f, size = 8)
  }
  writeBin(object = as.integer(x = K), con = f, size = 4) #K
  writeBin(object = as.numeric(x = sigma), con = f, size = 8) #sigma
  writeBin(object = as.integer(x = nbody_algo), con = f, size = 4)  #not barnes hut
  writeBin(object = as.integer(x = knn_algo), con = f, size = 4)
  writeBin(object = as.numeric(x = exaggeration_factor), con = f, size = 8) #compexag
  writeBin(object = as.integer(x = no_momentum_during_exag), con = f, size = 4)
  writeBin(object = as.integer(x = n_trees), con = f, size = 4)
  writeBin(object = as.integer(x = search_k), con = f, size = 4)
  writeBin(object = as.integer(x = start_late_exag_iter), con = f, size = 4)
  writeBin(object = as.numeric(x = late_exag_coeff), con = f, size = 8)
  writeBin(object = as.integer(x = nterms), con = f, size = 4)
  writeBin(object = as.numeric(x = intervals_per_integer), con = f, size = 8)
  writeBin(object = as.integer(x = min_num_intervals), con = f, size = 4)
  writeBin(object = tX, con = f)
  writeBin(object = as.integer(x = rand_seed), con = f, size = 4)
  if (version_number != "1.0.0") {
    writeBin(object = as.numeric(x = df), con = f, size = 8)
  }
  writeBin(object = as.integer(x = load_affinities), con = f, size = 4)
  if (!is.null(x = initialization)) {
    writeBin(object = c(t(x = initialization)), con = f)
  }
  close(con = f)

  if (version_number == "1.0.0") {
    flag <- system2(
      command = fast_tsne_path,
      args = c(data_path, result_path, nthreads)
    )
  } else {
    flag <- system2(
      command = fast_tsne_path,
      args = c(version_number, data_path, result_path, nthreads)
    )
  }

  if (flag != 0) {
    stop('tsne call failed')
  }

  f <- file(description = result_path, open = "rb")
  n <- readBin(con = f, what = integer(), n = 1, size = 4)
  d <- readBin(con = f, what = integer(), n = 1, size = 4)
  Y <- readBin(con = f, what = numeric(), n = n * d)
  Y <- t(x = matrix(Y, nrow = d))
  if (get_costs) {
    tmp <- readBin(con = f, what = integer(), n = 1, size = 4)
    costs <- readBin(con = f, what = numeric(), n = max_iter, size = 8)
    Yout <- list(Y = Y, costs = costs)
  } else {
    Yout <- Y
  }
  close(con = f)
  file.remove(data_path)
  file.remove(result_path)
  return(Yout)
}

#internal
#
JackRandom <- function(
  scaled.data,
  prop.use = 0.01,
  r1.use = 1,
  r2.use = 5,
  seed.use = 1,
  rev.pca = FALSE,
  weight.by.var = weight.by.var,
  maxit = 1000
) {
  if (!is.null(x = seed.use)) {
    set.seed(seed = seed.use)
  }
  rand.genes <- sample(
    x = rownames(x = scaled.data),
    size = nrow(x = scaled.data) * prop.use
  )
  # make sure that rand.genes is at least 3
  if (length(x = rand.genes) < 3) {
    rand.genes <- sample(x = rownames(x = scaled.data), size = 3)
  }
  data.mod <- scaled.data
  data.mod[rand.genes, ] <- MatrixRowShuffle(x = scaled.data[rand.genes, ])
  temp.object <- RunPCA(
    object = data.mod,
    assay = "temp",
    npcs = r2.use,
    features = rownames(x = data.mod),
    rev.pca = rev.pca,
    weight.by.var = weight.by.var,
    verbose = FALSE,
    maxit = maxit
  )
  return(Loadings(temp.object)[rand.genes, r1.use:r2.use])
}

# Calculates the l2-norm of a vector
#
# Modified from PMA package
# @references Witten, Tibshirani, and Hastie, Biostatistics 2009
# @references \url{https://github.com/cran/PMA/blob/master/R/PMD.R}
#
# @param vec numeric vector
#
# @return returns the l2-norm.
#
L2Norm <- function(vec) {
  a <- sqrt(x = sum(vec ^ 2))
  if (a == 0) {
    a <- .05
  }
  return(a)
}

# Prep data for dimensional reduction
#
# Common checks and preparatory steps before running certain dimensional
# reduction techniques
#
# @param object        Assay object
# @param features  Features to use as input for the dimensional reduction technique.
#                      Default is variable features
# @ param verbose   Print messages and warnings
#
#
PrepDR <- function(
  object,
  features = NULL,
  slot = 'scale.data',
  verbose = TRUE
) {
  if (length(x = VariableFeatures(object = object)) == 0 && is.null(x = features)) {
    stop("Variable features haven't been set. Run FindVariableFeatures() or provide a vector of feature names.")
  }
  data.use <- GetAssayData(object = object, slot = slot)
  if (nrow(x = data.use ) == 0 && slot == "scale.data") {
    stop("Data has not been scaled. Please run ScaleData and retry")
  }
  features <- features %||% VariableFeatures(object = object)
  features.keep <- unique(x = features[features %in% rownames(x = data.use)])
  if (length(x = features.keep) < length(x = features)) {
    features.exclude <- setdiff(x = features, y = features.keep)
    if (verbose) {
      warning(paste0("The following ", length(x = features.exclude), " features requested have not been scaled (running reduction without them): ", paste0(features.exclude, collapse = ", ")))
    }
  }
  features <- features.keep

  if (inherits(x = data.use, what = 'dgCMatrix')) {
    features.var <- RowVarSparse(mat = data.use[features, ])
  }
  else {
    features.var <- RowVar(x = data.use[features, ])
  }
  features.keep <- features[features.var > 0]
  if (length(x = features.keep) < length(x = features)) {
    features.exclude <- setdiff(x = features, y = features.keep)
    if (verbose) {
      warning(paste0("The following ", length(x = features.exclude), " features requested have zero variance (running reduction without them): ", paste0(features.exclude, collapse = ", ")))
    }
  }
  features <- features.keep
  features <- features[!is.na(x = features)]
  data.use <- data.use[features, ]
  return(data.use)
}

PrepDR5 <- function(object, features = NULL, layer = 'scale.data', verbose = TRUE) {
  layer <- layer[1L]
  layer <- match.arg(arg = layer, choices = Layers(object = object))
  features <- features %||% VariableFeatures(object = object)
  if (!length(x = features)) {
    stop("No variable features, run FindVariableFeatures() or provide a vector of features", call. = FALSE)
  }
  data.use <- LayerData(object = object, layer = layer, features = features)
  features.var <- apply(X = data.use, MARGIN = 1L, FUN = var)
  features.keep <- features[features.var > 0]
  if (!length(x = features.keep)) {
    stop("None of the requested features have any variance", call. = FALSE)
  } else if (length(x = features.keep) < length(x = features)) {
    exclude <- setdiff(x = features, y = features.keep)
    if (isTRUE(x = verbose)) {
      warning(
        "The following ",
        length(x = exclude),
        " features requested have zero variance; running reduction without them: ",
        paste(exclude, collapse = ', '),
        call. = FALSE,
        immediate. = TRUE
      )
    }
  }
  # features <- features.keep
  # features <- features[!is.na(x = features)]
  return(LayerData(object = object, layer = layer, features = features.keep))
}

#' @param assay Name of Assay SPCA is being run on
#' @param npcs Total Number of SPCs to compute and store (50 by default)
#' @param verbose Print the top genes associated with high/low loadings for
#' the SPCs
#' @param reduction.key dimensional reduction key, specifies the string before
#' the number for the dimension names. SPC by default
#' @param graph Graph used supervised by SPCA
#' @param seed.use Set a random seed. By default, sets the seed to 42. Setting
#' NULL will not set a seed.
#'
#' @importFrom irlba irlba
#'
#' @concept dimensional_reduction
#' @rdname RunSPCA
#' @export
RunSPCA.default <- function(
  object,
  assay = NULL,
  npcs = 50,
  reduction.key = "SPC_",
  graph = NULL,
  verbose = FALSE,
  seed.use = 42,
  ...
) {
  if (!is.null(x = seed.use)) {
    set.seed(seed = seed.use)
  }
  npcs <- min(npcs, nrow(x = object) - 1)
  if (verbose) {
    message("Computing sPCA transformation")
  }
  HSIC <- object %*% graph %*% t(x = object)
  pca.results <- irlba(A = HSIC, nv = npcs)
  feature.loadings <- pca.results$u
  rownames(x = feature.loadings) <- rownames(x = object)
  cell.embeddings <- t(object) %*% feature.loadings
  colnames(x = cell.embeddings) <- colnames(x = feature.loadings) <-
    paste0(reduction.key, 1:ncol(x = cell.embeddings))
  sdev <- pca.results$d / sqrt(max(1, nrow(x = HSIC) - 1))
  reduction.data <- CreateDimReducObject(
    embeddings = cell.embeddings,
    loadings = feature.loadings,
    assay = assay,
    stdev = sdev,
    key = reduction.key
  )
  return(reduction.data)
}

#' @param features Features to compute SPCA on. If features=NULL, SPCA will be run
#' using the variable features for the Assay.
#'
#' @rdname RunSPCA
#' @concept dimensional_reduction
#' @export
#' @method RunSPCA Assay
#'
RunSPCA.Assay <- function(
  object,
  assay = NULL,
  features = NULL,
  npcs = 50,
  reduction.key = "SPC_",
  graph = NULL,
  verbose = TRUE,
  seed.use = 42,
  ...
) {
  data.use <- PrepDR(
    object = object,
    features = features,
    verbose = verbose
  )
  reduction.data <- RunSPCA(
    object = data.use,
    assay = assay,
    npcs = npcs,
    reduction.key = reduction.key,
    graph = graph,
    verbose = verbose,
    seed.use = seed.use,
    ...
  )
  return(reduction.data)
}

#' @param features Features to compute SPCA on. If features=NULL, SPCA will be run
#' using the variable features for the Assay.
#'
#' @rdname RunSPCA
#' @concept dimensional_reduction
#' @export
#' @method RunSPCA Assay5
#'
RunSPCA.Assay5 <- function(
  object,
  assay = NULL,
  features = NULL,
  npcs = 50,
  reduction.key = "SPC_",
  graph = NULL,
  verbose = TRUE,
  seed.use = 42,
  layer = 'scale.data',
  ...
) {
  data.use <- PrepDR5(
    object = object,
    features = features,
    layer = layer,
    verbose = verbose
  )
  reduction.data <- RunSPCA(
    object = data.use,
    assay = assay,
    npcs = npcs,
    reduction.key = reduction.key,
    graph = graph,
    verbose = verbose,
    seed.use = seed.use,
    ...
  )
  return(reduction.data)
}

#' @param reduction.name dimensional reduction name, spca by default
#' @rdname RunSPCA
#' @concept dimensional_reduction
#' @export
#' @method RunSPCA Seurat
#'
RunSPCA.Seurat <- function(
  object,
  assay = NULL,
  features = NULL,
  npcs = 50,
  reduction.name = "spca",
  reduction.key = "SPC_",
  graph = NULL,
  verbose = TRUE,
  seed.use = 42,
  ...
) {
  assay <- assay %||% DefaultAssay(object = object)
  if (is.null(x = graph)) {
    stop("Graph is not provided")
  } else if (is.character(x = graph)) {
    graph <- object[[graph]]
  }
  reduction.data <- RunSPCA(
    object = object[[assay]],
    assay = assay,
    features = features,
    npcs = npcs,
    reduction.name = reduction.name,
    reduction.key = reduction.key,
    graph = graph,
    verbose = verbose,
    seed.use = seed.use,
    ...
  )
  object[[reduction.name]] <- reduction.data
  object <- LogSeuratCommand(object = object)
  return(object)
}

#' @param assay Name of Assay SLSI is being run on
#' @param n Total Number of SLSI components to compute and store
#' @param verbose Display messages
#' @param reduction.key dimensional reduction key, specifies the string before
#' the number for the dimension names
#' @param graph Graph used supervised by SLSI
#' @param seed.use Set a random seed. Setting NULL will not set a seed.
#'
#' @importFrom irlba irlba
#' @importMethodsFrom Matrix t
#'
#' @concept dimensional_reduction
#' @rdname RunSLSI
#' @export
RunSLSI.default <- function(
  object,
  assay = NULL,
  n = 50,
  reduction.key = "SLSI_",
  graph = NULL,
  verbose = TRUE,
  seed.use = 42,
  ...
) {
  if (!is.null(x = seed.use)) {
    set.seed(seed = seed.use)
  }
  n <- min(n, nrow(x = object) - 1)

  if (verbose) {
    message("Smoothing peaks matrix")
  }
  object.smooth <- t(x = graph) %*% (t(x = object) %*% object) %*% graph
  if (verbose) {
    message("Performing eigendecomposition")
  }
  svd.V <- irlba(A = object.smooth, nv = n, nu = n, ...)
  sigma <- sqrt(x = svd.V$d)
  feature.loadings <- object %*% (graph %*% svd.V$u) %*% diag(x = 1/sigma)
  feature.loadings <- as.matrix(x = feature.loadings)
  cell.embeddings <- t(x = object) %*% feature.loadings %*% diag(x = 1/sigma)
  cell.embeddings <- as.matrix(x = cell.embeddings)

  # construct svd list stored in misc for LSI projection
  svd.lsi <- svd.V
  svd.lsi$d <- sigma
  svd.lsi$u <- feature.loadings
  svd.lsi$v <- cell.embeddings

  colnames(x = cell.embeddings) <- paste0(reduction.key, 1:ncol(cell.embeddings))
  reduction.data <- CreateDimReducObject(
    embeddings = cell.embeddings,
    loadings = feature.loadings,
    key = reduction.key,
    assay = assay,
    misc = svd.lsi
  )
  return(reduction.data)
}

#' @param features Features to compute SLSI on. If NULL, SLSI will be run
#' using the variable features for the Assay.
#'
#' @rdname RunSLSI
#' @concept dimensional_reduction
#' @export
#' @method RunSLSI Assay
#'
RunSLSI.Assay <- function(
  object,
  assay = NULL,
  features = NULL,
  n = 50,
  reduction.key = "SLSI_",
  graph = NULL,
  verbose = TRUE,
  seed.use = 42,
  ...
) {
  data.use <- PrepDR(
    object = object,
    features = features,
    slot = "data",
    verbose = verbose
  )
  reduction.data <- RunSLSI(
    object = data.use,
    assay = assay,
    n = n,
    reduction.key = reduction.key,
    graph = graph,
    verbose = verbose,
    seed.use = seed.use,
    ...
  )
  return(reduction.data)
}

#' @param reduction.name dimensional reduction name
#' @rdname RunSLSI
#' @concept dimensional_reduction
#' @export
#' @method RunSLSI Seurat
#'
RunSLSI.Seurat <- function(
  object,
  assay = NULL,
  features = NULL,
  n = 50,
  reduction.name = "slsi",
  reduction.key = "SLSI_",
  graph = NULL,
  verbose = TRUE,
  seed.use = 42,
  ...
) {
  assay <- assay %||% DefaultAssay(object = object)
  assay.data <- GetAssay(object = object, assay = assay)
  if (is.null(x = graph)) {
    stop("Graph is not provided")
  } else if (is.character(x = graph)) {
    graph <- object[[graph]]
  }
  reduction.data <- RunSLSI(
    object = assay.data,
    assay = assay,
    features = features,
    n = n,
    reduction.name = reduction.name,
    reduction.key = reduction.key,
    graph = graph,
    verbose = verbose,
    seed.use = seed.use,
    ...
  )
  object[[reduction.name]] <- reduction.data
  object <- LogSeuratCommand(object = object)
  return(object)
}
