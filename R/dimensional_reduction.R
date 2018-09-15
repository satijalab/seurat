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
#' @return Returns a Seurat object where object@@dr$pca@@jackstraw@@emperical.p.value
#' represents p-values for each gene in the PCA analysis. If ProjectPCA is
#' subsequently run, object@dr$pca@jackstraw@emperical.p.value.full then
#' represents p-values for all genes.
#'
#' @import doSNOW
#' @importFrom methods new
#' @importFrom utils txtProgressBar setTxtProgressBar
#'
#' @references Inspired by Chung et al, Bioinformatics (2014)
#'
#' @export
#'
#' @examples
#' \dontrun{
#' pbmc_small = suppressWarnings(JackStraw(pbmc_small))
#' head(pbmc_small@dr$pca@jackstraw@emperical.p.value)
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
  assay <- assay %||% DefaultAssay(object = object)
  # embeddings <- Embeddings(object = object[[reduction]])
  if (dims > length(x = object[[reduction]])) {
    dims <- length(x = object[[reduction]])
    warning("Number of dimensions specified is greater than those available. Setting dims to ", dims, " and continuing", immediate. = TRUE)
  }
  if (dims > ncol(x = object)) {
    dims <- ncol(x = object)
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
  ## TODO: Parallelization
  fake.vals.raw <- lapply(
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
    sapply(
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

#' Canonical Correlation Analysis between modalities
#'
#' Runs a canonical correlation analysis between two assays (for example, RNA and ADT from CITE-seq)
#'
#' @param object Seurat object. Assumes that scale.data exist for both assays.
#' @param assay.1 The first assay for CCA.
#' @param assay.2 The second assay for CCA.
#' @param features.1 Features to use for the first assay, default is all the features (use object@var.genes if this is RNA).
#' @param features.2 Features to use for the second assay, default is all the features.
#' @param num.cc Minimal number of CCs to return.
#' @param normalize.variance Whether to scale the embeddings. Default is TRUE.
#'
#' @return Returns Seurat object with two CCA results stored (for example, object@dr$RNACCA and object@dr$ADTCCA).
#'
#' @export
#'
#' @examples
#' \dontrun{
#' object <- MultiModal_CCA(object,assay.1 = "RNA",assay.2 = "ADT")
#' }
#'
MultiModal_CCA <- function(
  object,
  assay.1 = "RNA",
  assay.2 = "ADT",
  features.1 = NULL,
  features.2 = NULL,
  num.cc = 20,
  normalize.variance = TRUE
) {
  #first pull out data, define features
  data.1 <- GetAssayData(
    object = object,
    assay.type = assay.1,
    slot = "scale.data"
  )
  data.2 <- GetAssayData(
    object = object,
    assay.type = assay.2,
    slot = "scale.data"
  )
  if (is.null(x = features.1)) {
    if ((assay.1 == "RNA") && length(x = object@var.genes) > 0) {
      features.1 <- object@var.genes
    } else {
      features.1 <- rownames(x = data.1)
    }
  }
  features.2 <- SetIfNull(x = features.2, default = rownames(x = data.2))
  data.1 <- t(x = data.1[features.1, ])
  data.2 <- t(x = data.2[features.2, ])
  #data.1.var=apply(data.1,2,var)
  #data.2.var=apply(data.2,2,var)
  data.1 <- data.1[,apply(data.1,2,var)>0]
  data.2 <- data.2[,apply(data.2,2,var)>0]

  num.cc <- max(20, min(ncol(data.1), ncol(data.2)))
  cca.data <- list(data.1, data.2)
  names(x = cca.data) <- c(assay.1, assay.2)
  # now run CCA
  out <- CCA(
    x = cca.data[[1]],
    z = cca.data[[2]],
    typex = "standard",
    typez = "standard",
    K = num.cc,
    penaltyz = 1,
    penaltyx = 1
  )
  cca.output <- list(out$u, out$v)
  embeddings.cca <- list()
  for (i in 1:length(x = cca.data)) {
    assay <- names(x = cca.data)[i]
    rownames(x = cca.output[[i]]) <- colnames(x = cca.data[[i]])
    embeddings.cca[[i]] <- cca.data[[i]] %*% cca.output[[i]]
    colnames(x = embeddings.cca[[i]]) <- paste0(
      assay,
      "CC",
      1:ncol(x = embeddings.cca[[i]])
    )
    colnames(x = cca.output[[i]]) <- colnames(x = embeddings.cca[[i]])
    if (normalize.variance) {
      embeddings.cca[[i]] <- scale(x = embeddings.cca[[i]])
    }
    object <- SetDimReduction(
      object = object,
      reduction.type = paste0(assay, "CCA"),
      slot = "cell.embeddings",
      new.data = embeddings.cca[[i]]
    )
    object <- SetDimReduction(
      object = object,
      reduction.type = paste0(assay, "CCA"),
      slot = "key",
      new.data = paste0(assay, "CC")
    )
    object <- SetDimReduction(
      object = object,
      reduction.type = paste0(assay, "CCA"),
      slot = "gene.loadings",
      new.data =  cca.output[[i]]
    )
  }
  return(object)
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
#'
#' @seealso \code{\link{ProjectDim}} \code{\link{JackStraw}}
#'
#' @examples
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
#'
#' @examples
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
    object = object,
    assay.us = assay,
    slot = "scale.data"
  )
  if (do.center) {
    data.use <- scale(x = as.matrix(x = data.use), center = TRUE, scale = FALSE)
  }
  cell.embeddings <- Embeddings(object = redeuc)
  new.feature.loadings.full <- FastMatMult(m1 = data.use, m2 = cell.embeddings)
  rownames(x = new.feature.loadings.full) <- rownames(x = data.use)
  colnames(x = new.feature.loadings.full) <- colnames(x = cell.embeddings)
  Loadings(object = redeuc, projected = TRUE) <- new.feature.loadings.full
  if (overwrite) {
    Loadings(object = redeuc, projected = FALSE) <- new.feature.loadings.full
  }
  object[[reduction]] <- redeuc
  if (verbose) {
    Print(
      object = redeuc,
      dims = dims.print,
      nfeatures = nfeatures.print,
      projected = TRUE
    )
  }
  object <- LogSeuratCommand(object = object)
  return(object)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Methods for Seurat-defined generics
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @param standardize Standardize matrices - scales columns to have unit variance
#' and mean 0
#' @param num.cc Number of canonical vectors to calculate
#' @param verbose ...
#' @param use.cpp ...
#'
#' @importFrom irlba irlba
#'
#' @rdname RunCCA
#' @export
#'
RunCCA.default <- function(
  object1,
  object2,
  standardize = TRUE,
  num.cc = 20,
  verbose = FALSE,
  use.cpp = TRUE,
  ...
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

#' @param assay1,assay2 Assays to pull from in the first and second objects, respectively
#' @param features Set of genes to use in CCA. Default is the union of both
#' the variable features sets present in both objects.
#' @param renormalize Renormalize raw data after merging the objects. If FALSE,
#' merge the data matrices also.
#' @param rescale Rescale the datasets prior to CCA. If FALSE, uses existing data in the scale data slots.
#' @param compute.gene.loadings Also compute the gene loadings. NOTE - this will
#' scale every gene in the dataset which may impose a high memory cost.
#' @param add.cell.id1,add.cell.id2 Add ...
#'
#' @rdname RunCCA
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
  use.cpp = TRUE,
  ...
) {
  assay1 <- assay1 %||% DefaultAssay(object = object1)
  assay2 <- assay2 %||% DefaultAssay(object = object2)

  if(assay1 != assay2) {
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
  data.use1 <- GetAssayData(object = object1, assay = assay1, slot = "scale.data")
  data.use2 <- GetAssayData(object = object2, assay = assay2, slot = "scale.data")
  features <- CheckFeatures(data.use = data.use1, features = features, object.name = "object1")
  features <- CheckFeatures(data.use = data.use2, features = features, object.name = "object2")


  if (length(x = features) < 50) {
    warning("Fewer than 50 featuresd as input for CCA.")
  }
  if (verbose) {
    message("Running CCA")
  }
  data1 <- data.use1[features,]
  data2 <- data.use2[features,]
  if (rescale) {
    if (verbose) message("Rescaling groups")
    data1 <- FastRowScale(data1)
    dimnames(data1) <- list(features, colnames(x = object1))
    data2 <- FastRowScale(data2)
    dimnames(data2) <- list(features, colnames(x = object2))
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
    merge.data = TRUE,
    ...
  )

  combined.object[['cca']] <- CreateDimReducObject(
    embeddings = cca.results$ccv[colnames(combined.object), ],
    assay = assay1,
    key = "CC"
  )
  if (ncol(combined.object) != (ncol(object1) + ncol(object2))) {
    warning("Some cells removed after object merge due to minimum feature count cutoff")
  }
  combined.scale <- cbind(data1,data2)
  combined.object <- SetAssayData(object = combined.object,new.data = combined.scale, assay = assay1, slot = "scale.data")
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

#' @param niter Number of iterations to perform. Set by default to 25.
#' @param num.ccs Number of canonical vectors to calculate
#' @param standardize standardize scale.data matrices to be centered (mean zero)
#' and scaled to have a standard deviation of 1.
#' @param verbose ...
#'
#' @importFrom irlba irlba
#'
#' @rdname RunMultiCCA
#' @export
#'
RunMultiCCA.default <- function(
  object.list,
  niter = 25,
  num.ccs = 1,
  standardize = TRUE,
  verbose = TRUE,
  ...
) {
  cell.names <- c()
  set.seed(seed = 42)
  for (object in object.list) {
    cell.names <- c(cell.names, colnames(x = object))
    if (!class(x = object) %in% c("matrix", "dgCMatrix")) {
      stop("Not all objects in object.list are matrices")
    }
  }
  num.sets <- length(x = object.list)
  if (standardize) {
    for (i in 1:num.sets){
      object.list[[i]] <- Standardize(object.list[[i]], display_progress = FALSE)
    }
  }
  ws <- list()
  for (i in 1:num.sets) {
    ws[[i]] <- irlba(object.list[[i]], nv = num.ccs)$v[, 1:num.ccs, drop = FALSE]
  }
  ws.init <- ws
  ws.final <- list()
  cors <- NULL
  for (i in 1:length(x = ws)) {
    ws.final[[i]] <- matrix(0, nrow = ncol(object.list[[i]]), ncol = num.ccs)
  }
  for (cc in 1:num.ccs ) {
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
  for (i in 2:length(x = object.list)) {
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

#' @param assay Assay to use
#' @param features Set of genes to use in CCA. Default is the union of both
#' the variable features sets present in both objects.
#' @param add.cell.ids ...
#' @param renormalize Renormalize raw data after merging the objects. If FALSE,
#' merge the data matrices also.
#' @param compute.gene.loadings Also compute the gene loadings. NOTE - this will
#' scale every gene in the dataset which may impose a high memory cost.
#'
#' @rdname RunMultiCCA
#' @export
#' @method RunMultiCCA Seurat
#'
RunMultiCCA.Seurat <- function(
  object.list,
  assay = NULL,
  features = NULL,
  add.cell.ids = NULL,
  niter = 25,
  num.ccs = 1,
  standardize = TRUE,
  renormalize = TRUE,
  compute.gene.loadings = TRUE,
  verbose = TRUE,
  ...
) {
  set.seed(seed = 42)
  for(object in object.list) {
    if (class(x = object) != "Seurat") {
      stop("Not all objects in object.list are Seurat objects")
    }
  }
  assay <- assay %||% DefaultAssay(object = object.list[[1]])
  for(object in object.list) {
    assays <- FilterObjects(object = object, classes.keep = "Assay")
    if (!assay %in% assays) {
      stop(paste0(assay, " not present in all objects."))
    }
  }
  if (is.null(x = features)) {
    features <- c()
    for(object in object.list) {
      features <- c(features, VariableFeatures(object = object))
    }
    features <- unique(x = features)
  }
  for(i in 1:length(x = object.list)) {
    features <- CheckFeatures(
      data.use = GetAssayData(object = object.list[[i]], assay = assay, slot = "scale.data"),
      features = features,
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
      assay = assay,
      slot = "scale.data"
    )[features, ]
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
  combined.object[['cca']] <- CreateDimReducObject(
    embeddings = cca.results$ccv,
    assay = assay,
    key = "CC"
  )
  if (renormalize) {
    combined.object <- NormalizeData(
      object = combined.object,
      assay = assay,
      normalization.method = object.list[[1]][[paste0("NormalizeData.", assay)]]$normalization.method,
      scale.factor = object.list[[1]][[paste0("NormalizeData.", assay)]]$scale.factor,
      verbose = verbose
    )
    if (compute.gene.loadings) {
      combined.object <- ScaleData(object = combined.object, verbose = verbose)
      combined.object <- ProjectDim(
        object = combined.object,
        reduction = "cca",
        verbose = FALSE,
        overwrite = TRUE)
    }
  }
  combined.object <- LogSeuratCommand(object = combined.object)
  return(combined.object)
}

#' @param assay Name of Assay PCA is being run on
#' @param npcs Total Number of PCs to compute and store (20 by default)
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
#'
#' @importFrom irlba irlba
#'
#' @rdname RunPCA
#' @export
#'
RunPCA.default <- function(
  object,
  assay = NULL,
  npcs = 20,
  rev.pca = FALSE,
  weight.by.var = TRUE,
  verbose = TRUE,
  ndims.print = 1:5,
  nfeatures.print = 30,
  reduction.key = "PC",
  seed.use = 42,
  ...
) {
  if (!is.null(x = seed.use)) {
    set.seed(seed = seed.use)
  }
  if (rev.pca) {
    npcs <- min(npcs, ncol(x = object) - 1)
    pca.results <- irlba(A = object, nv = npcs, ...)
    sdev <- pca.results$d/sqrt(max(1, nrow(x = object) - 1))
    if (weight.by.var) {
      feature.loadings <- pca.results$u %*% diag(pca.results$d)
    } else{
      feature.loadings <- pca.results$u
    }
    cell.embeddings <- pca.results$v
  }
  else {
    npcs <- min(npcs, nrow(x = object) - 1)
    pca.results <- irlba(A = t(x = object), nv = npcs, ...)
    feature.loadings <- pca.results$v
    sdev <- pca.results$d/sqrt(max(1, ncol(object) - 1))
    if (weight.by.var) {
      cell.embeddings <- pca.results$u %*% diag(pca.results$d)
    } else {
      cell.embeddings <- pca.results$u
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
    key = reduction.key
  )
  if (verbose) {
    Print(object = reduction.data, dims = ndims.print, nfeatures = nfeatures.print)
  }
  return(reduction.data)
}

#' @param features Features to compute PCA on
#'
#' @rdname RunPCA
#' @export
#' @method RunPCA Assay
#'
RunPCA.Assay <- function(
  object,
  assay = NULL,
  features = NULL,
  npcs = 20,
  rev.pca = FALSE,
  weight.by.var = TRUE,
  verbose = TRUE,
  ndims.print = 1:5,
  nfeatures.print = 30,
  reduction.key = "PC",
  seed.use = 42,
  ...
) {
  data.use <- PrepDR(
    object = object,
    features = features
  )
  reduction.data <- RunPCA(
    object = data.use,
    assay = assay,
    pc.features = features,
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

#' @param reduction.name dimensional reduction name,  pca by default
#' @param workflow.name Name of workflow
#'
#' @rdname RunPCA
#' @export
#' @method RunPCA Seurat
#'
RunPCA.Seurat <- function(
  object,
  assay = NULL,
  features = NULL,
  npcs = 20,
  rev.pca = FALSE,
  weight.by.var = TRUE,
  verbose = TRUE,
  ndims.print = 1:5,
  nfeatures.print = 30,
  reduction.name = "pca",
  reduction.key = "PC",
  seed.use = 42,
  workflow.name = NULL,
  ...
) {
  if (!is.null(workflow.name)) {
    object <- PrepareWorkflow(object = object, workflow.name = workflow.name)
  }
  assay <- assay %||% DefaultAssay(object = object)
  assay.data <- GetAssay(object = object, assay = assay)
  reduction.data <- RunPCA(
    object = assay.data,
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
  if (!is.null(workflow.name)) {
    object <- UpdateWorkflow(object = object, workflow.name = workflow.name)
  }
  return(object)
}

#' @param assay Name of assay that that t-SNE is being run on
#' @param seed.use Random seed for the t-SNE
#' @param tsne.method Select the method to use to compute the tSNE. Available
#' methods are:
#' \itemize{
#' \item{Rtsne: }{Use the Rtsne package Barnes-Hut implementation of tSNE (default)}
# \item{tsne: }{standard tsne - not recommended for large datasets}
#' \item{FIt-SNE: }{Use the FFT-accelerated Interpolation-based t-SNE. Based on
#' Kluger Lab code found here: https://github.com/KlugerLab/FIt-SNE}
#' }
#' @param add.iter If an existing tSNE has already been computed, uses the
#' current tSNE to seed the algorithm and then adds additional iterations on top
#' of this
#' @param dim.embed The dimensional space of the resulting tSNE embedding
#' (default is 2). For example, set to 3 for a 3d tSNE
#' @param reduction.key dimensional reduction key, specifies the string before the number for the dimension names. tSNE_ by default
#'
#' @importFrom tsne tsne
#' @importFrom Rtsne Rtsne
#'
#' @rdname RunTSNE
#' @export
#' @method RunTSNE matrix
#'
RunTSNE.matrix <- function(
  object,
  assay = NULL,
  seed.use = 1,
  tsne.method = "Rtsne",
  add.iter = 0,
  dim.embed = 2,
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
    stop("Invalid tSNE method: please choose from 'Rtsne' or 'FIt-SNE'")
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
  tsne.reduction <- CreateDimReducObject(
    embeddings = tsne.data,
    key = reduction.key,
    assay = assay
  )
  return(tsne.reduction)
}

#' @param cells Which cells to analyze (default, all cells)
#' @param dims Which dimensions to use as input features
#'
#' @rdname RunTSNE
#' @export
#' @method RunTSNE DimReduc
#'
RunTSNE.DimReduc <- function(
  object,
  cells = NULL,
  dims = 1:5,
  seed.use = 1,
  tsne.method = "Rtsne",
  add.iter = 0,
  dim.embed = 2,
  reduction.key = "tSNE_",
  ...
) {
  tsne.reduction <- RunTSNE(
    object = object[[, dims]],
    assay = DefaultAssay(object = object),
    seed.use = seed.use,
    tsne.method = tsne.method,
    add.iter = add.iter,
    dim.embed = dim.embed,
    reduction.key = reduction.key,
    ...
  )
  return(tsne.reduction)
}

#' @param reduction Which dimensional reduction (e.g. PCA, ICA) to use for
#' the tSNE. Default is PCA
#' @param features If set, run the tSNE on this subset of features
#' (instead of running on a set of reduced dimensions). Not set (NULL) by default
#' @param distance.matrix If set, runs tSNE on the given distance matrix
#' instead of data matrix (experimental)
#' @param reduction.name dimensional reduction name, specifies the position in the object$dr list. tsne by default
#' @param workflow.name Name of workflow
#'
#' @rdname RunTSNE
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
  if (length(dims) == 1 && !is.null(workflow.name)) {
    dims <- 1:dims
  }
  tsne.reduction <- if (!is.null(x = distance.matrix)) {
    RunTSNE(
      object = as.matrix(x = distance.matrix),
      assay = DefaultAssay(object = object),
      seed.use = seed.use,
      tsne.method = tsne.method,
      add.iter = add.iter,
      dim.embed = dim.embed,
      reduction.key = reduction.key,
      is_distance = TRUE,
      ...
    )
  } else if (!is.null(x = features)) {
    RunTSNE(
      object = as.matrix(x = GetAssayData(object = object)[features, ]),
      assay = DefaultAssay(object = object),
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
      object = object[[reduction]],
      dims = dims,
      seed.use = seed.use,
      tsne.method = tsne.method,
      add.iter = add.iter,
      dim.embed = dim.embed,
      reduction.key = reduction.key,
      pca = FALSE,
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

#' @param assay Assay to pull data for when using \code{genes.use}
#' @param nneighbors This determines the number of neighboring points used in
#' local approximations of manifold structure. Larger values will result in more
#' global structure being preserved at the loss of detailed local structure. In
#' general this parameter should often be in the range 5 to 50.
#' @param max.dim Max dimension to keep from UMAP procedure.
#' @param metric metric: This determines the choice of metric used to measure
#' distance in the input space. A wide variety of metrics are already coded, and
#' a user defined function can be passed as long as it has been JITd by numba.
#' @param min.dist min_dist: This controls how tightly the embedding is allowed
#' compress points together. Larger values ensure embedded points are more
#' evenly distributed, while smaller values allow the algorithm to optimise more
#' accurately with regard to local structure. Sensible values are in the range
#' 0.001 to 0.5.
#' @param reduction.key dimensional reduction key, specifies the string before
#' the number for the dimension names. UMAP by default
#' @param seed.use Set a random seed. By default, sets the seed to 42. Setting
#' NULL will not set a seed
#'
#' @importFrom reticulate py_module_available py_set_seed import
#'
#' @rdname RunUMAP
#' @export
#'
RunUMAP.default <- function(
  object,
  assay = NULL,
  nneighbors = 30L,
  max.dim = 2L,
  metric = 'correlation',
  min.dist = 0.3,
  reduction.key = 'umap',
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
  umap_import <- import(module = "umap", delay_load = TRUE)
  umap <- umap_import$UMAP(
    n_neighbors = as.integer(x = nneighbors),
    n_components = as.integer(x = max.dim),
    metric = metric,
    min_dist = min.dist
  )
  umap_output <- umap$fit_transform(as.matrix(x = object))
  colnames(x = umap_output) <- paste0(reduction.key, 1:ncol(x = umap_output))
  rownames(x = umap_output) <- colnames(object)
  umap.reduction <- CreateDimReducObject(
    embeddings = umap_output,
    key = reduction.key,
    assay = assay
  )
  return(umap.reduction)
}

#' @param dims Which dimensions to use as input features, used only if
#' \code{genes.use} is NULL
#' @param reduction Which dimensional reduction (PCA or ICA) to use for the
#' UMAP input. Default is PCA
#' @param features If set, run UMAP on this subset of features (instead of running on a
#' set of reduced dimensions). Not set (NULL) by default
#' @param reduction.name dimensional reduction name, specifies the position in
#' the object$dr list. umap by default
#'
#' @rdname RunUMAP
#' @export
#' @method RunUMAP Seurat
#'
RunUMAP.Seurat <- function(
  object,
  dims = 1:5,
  reduction = 'pca',
  features = NULL,
  assay = 'RNA',
  nneighbors = 30L,
  max.dim = 2L,
  min.dist = 0.3,
  reduction.name = "umap",
  reduction.key = "UMAP",
  metric = "correlation",
  seed.use = 42,
  ...
) {
  data.use <- if (is.null(x = features)) {
    Embeddings(object[[reduction]])[, dims]
    assay <- DefaultAssay(object = object[['reduction']])
  } else {
    t(x = GetAssayData(object = object, slot = 'data', assay = assay)[features, ])
  }
  object[[reduction.name]] <- RunUMAP(
    object = data.use,
    assay = assay,
    nneighbors = nneighbors,
    max.dim = max.dim,
    metric = metric,
    min.dist = min.dist,
    reduction.key = reduction.key,
    seed.use = seed.use,
    ...
  )
  return(object)
}

#' @param dims Which dimensions to examine
#' @param score.thresh Threshold to use for the proportion test of PC
#' significance (see Details)
#'
#' @importFrom stats prop.test
#'
#' @rdname ScoreJackStraw
#' @export
#' @method ScoreJackStraw JackStrawData
#'
ScoreJackStraw.JackStrawData <- function(
  object,
  dims = 1:5,
  score.thresh = 1e-5,
  ...
) {
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
#' @export
#' @method ScoreJackStraw DimReduc
#'
ScoreJackStraw.DimReduc <- function(object, dims = 1:5, score.thresh = 1e-5, ...) {
  JS(object = object) <- ScoreJackStraw(
    object = JS(object = object),
    dims = dims,
    score.thresh = 1e-5,
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
    ...
  )
  if (do.plot) {
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
# Internal
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Calculate a low dimensional projection of the data. First forms an orthonormal
# basis of the gene loadings via QR decomposition, projects the data onto that
# basis, and reconstructs the data using on the dimensions specified.
#
# @param object          Seurat object
# @param reduction.type  Type of dimensional reduction to use
# @param dims        Dimensions to use in calculation
# @param genes.use       Genes to consider when calculating
#
# @return                Returns a matrix with the low dimensional reconstruction
#
CalcLDProj <- function(object, reduction.type, dims, genes.use) {
  if (missing(x = dims)) {
    dims <- 1:ncol(x = Embeddings(
      object = object,
      reduction.type = reduction.type
    ))
  }
  x.vec <- Loadings(
    object = object,
    reduction.type = reduction.type,
    dims = dims,
    genes.use = genes.use
  )
  # form orthonormal basis via QR
  x.norm <- qr.Q(qr = qr(x = x.vec))
  data.use <- object@scale.data[rownames(x.vec), ]
  # project data onto othronormal basis
  projected.data <- t(x = data.use) %*% x.norm
  # reconstruct data using only dims specified
  low.dim.data <- x.norm %*% t(x = projected.data)
  return(low.dim.data)
}

# Calculate percent variance explained
#
# Projects dataset onto the orthonormal space defined by some dimensional
# reduction technique (e.g. PCA, CCA) and calculates the percent of the
# variance in gene expression explained by each cell in that lower dimensional
# space.
#
# @param object          Seurat object
# @param reduction.type  Name of the reduction to use for the projection
# @param dims        Vector of dimensions to project onto (default is the
#                        1:number stored for given technique)
# @param genes.use       vector of genes to use in calculation
#
# @return                Returns a Seurat object wih the variance in gene
#                        expression explained by each cell in a low dimensional
#                        space stored as metadata.
#
CalcProjectedVar <- function(
  object,
  low.dim.data,
  reduction.type = "pca",
  dims,
  genes.use
) {
  if (missing(x = low.dim.data)) {
    low.dim.data <- CalcLDProj(
      object = object,
      reduction.type = reduction.type,
      dims = dims,
      genes.use = genes.use
    )
  }
  projected.var <- apply(X = low.dim.data, MARGIN = 2, FUN = var)
  calc.name <- paste0(reduction.type, ".var")
  object <- AddMetaData(
    object = object,
    metadata = projected.var,
    col.name = calc.name
  )
  return(object)
}

# Check that features are present and have non-zero variance
#
# @param data.use      Feature matrix (features are rows)
# @param features  Features to check
# @param object.name   Name of object for message printing
#
# @return           Returns a vector of features that is the subset of features
#                   that have non-zero variance
#
CheckFeatures <- function(data.use, features, object.name) {
  if (any(!features %in% rownames(data.use))) {
    missing.features <- features[!features %in% rownames(data.use)]
    stop(
      "Following features are not scaled in ",
      object.name,
      ": ",
      paste0(missing.features, collapse = ", ")
    )
  }
  features.var <- apply(X = data.use[features, ], MARGIN = 1, FUN = var)
  no.var.features <- features[features.var == 0]
  if (length(no.var.features) > 0) {
    warning(
      "The following features have zero variance in ",
      object.name,
      ": ",
      paste0(no.var.features, collapse = ", ")
    )
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
# Based on Kluger Lab code on https://github.com/ChristophH/FIt-SNE
# commit ec25f1b36598a2d21869d10a258ac366a12f0b05
#
#' @importFrom utils file_test
#
fftRtsne <- function(
  X,
  dims = 2,
  perplexity = 30,
  theta = 0.5,
  check_duplicates = TRUE,
  max_iter = 1000,
  fft_not_bh = TRUE,
  ann_not_vptree = TRUE,
  stop_lying_iter = 250,
  exaggeration_factor = 12.0,
  no_momentum_during_exag = FALSE,
  start_late_exag_iter = -1.0,
  late_exag_coeff = 1.0,
  n_trees = 50,
  search_k = -1,
  rand_seed = -1,
  nterms = 3,
  intervals_per_integer = 1,
  min_num_intervals = 50,
  data_path = NULL,
  result_path = NULL,
  fast_tsne_path = NULL,
  nthreads = getOption('mc.cores', default = 1),
  ...
) {
  if (is.null(x = data_path)) {
    data_path <- tempfile(pattern = 'fftRtsne_data_', fileext = '.dat')
  }
  if (is.null(x = result_path)) {
    result_path <- tempfile(pattern = 'fftRtsne_result_', fileext = '.dat')
  }
  if (is.null(x = fast_tsne_path)) {
    fast_tsne_path <- system2(command = 'which', args = 'fast_tsne', stdout = TRUE)
    if (length(x = fast_tsne_path) == 0) {
      stop("no fast_tsne_path specified and fast_tsne binary is not in the search path")
    }
  }
  fast_tsne_path <- normalizePath(path = fast_tsne_path)
  if (!file_test(op = '-x', x = fast_tsne_path)) {
    stop("fast_tsne_path '", fast_tsne_path, "' does not exist or is not executable")
  }
  is.wholenumber <- function(x, tol = .Machine$double.eps ^ 0.5) {
    return(abs(x = x - round(x = x)) < tol)
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
  if (!is.wholenumber(x = stop_lying_iter) || stop_lying_iter < 0) {
    stop("stop_lying_iter should be a positive integer")
  }
  if (!is.numeric(x = exaggeration_factor)) {
    stop("exaggeration_factor should be numeric")
  }
  if (!is.wholenumber(x = dims) || dims <= 0) {
    stop("Incorrect dimensionality.")
  }
  if (search_k == -1) {
    search_k = n_trees * perplexity * 3
  }
  # if (fft_not_bh) {
  #   nbody_algo <- 2
  # } else {
  #   nbody_algo <- 1
  # }
  nbody_algo <- ifelse(test = fft_not_bh, yes = 2, no = 1)
  # if (ann_not_vptree) {
  #   knn_algo <- 1
  # }else{
  #   knn_algo <- 2
  # }
  knn_algo <- ifelse(test = ann_not_vptree, yes = 1, no = 2)
  tX = c(t(x = X))
  f <- file(data_path, "wb")
  n = nrow(x = X)
  D = ncol(x = X)
  writeBin(object = as.integer(x = n), con = f, size = 4)
  writeBin(object = as.integer(x = D), con = f, size = 4)
  writeBin(object = as.numeric(x = 0.5), con = f, size = 8) #theta
  writeBin(object = as.numeric(x = perplexity), con = f, size = 8) #theta
  writeBin(object = as.integer(x = dims), con = f, size = 4) #theta
  writeBin(object = as.integer(x = max_iter), con = f, size = 4)
  writeBin(object = as.integer(x = stop_lying_iter), con = f, size = 4)
  writeBin(object = as.integer(x = -1), con = f, size = 4) #K
  writeBin(object = as.numeric(x = -30.0), con = f, size = 8) #sigma
  writeBin(object = as.integer(x = nbody_algo), con = f, size = 4)  #not barnes hut
  writeBin(object = as.integer(x = knn_algo), con = f, size = 4)
  writeBin(object = as.numeric(x = exaggeration_factor), con = f, size = 8) #compexag
  writeBin(object = as.integer(x = no_momentum_during_exag), con = f, size = 4)
  writeBin(object = as.integer(x = n_trees), con = f, size = 4)
  writeBin(object = as.integer(x = search_k), con = f, size = 4)
  writeBin(object = as.integer(x = start_late_exag_iter), con = f, size = 4)
  writeBin(object = as.numeric(x = late_exag_coeff), con = f, size = 8)
  writeBin(object = as.integer(x = nterms), con =  f, size = 4)
  writeBin(object = as.numeric(x = intervals_per_integer), con =  f, size = 8)
  writeBin(object = as.integer(x = min_num_intervals), con =  f, size = 4)
  tX = c(t(x = X))
  writeBin(object = tX, con = f)
  writeBin(object = as.integer(x = rand_seed), con = f, size = 4)
  close(f)
  flag <- system2(command = fast_tsne_path, args = c(data_path, result_path, nthreads))
  if (flag != 0) {
    stop('tsne call failed');
  }
  f <- file(description = result_path, open = "rb")
  initialError <- readBin(f, integer(), n = 1, size = 8)
  n <- readBin(con = f, what = integer(), n = 1, size = 4)
  d <- readBin(con = f, what = integer(), n = 1, size = 4)
  Y <- readBin(con = f, what = numeric(), n = n * d)
  Yout <- t(x = matrix(data = Y, nrow = d))
  close(f)
  file.remove(data_path)
  file.remove(result_path)
  return(Yout)
}

# MultiCCA helper function - calculates correlation
#
# Modified from PMA package
# @references Witten, Tibshirani, and Hastie, Biostatistics 2009
# @references \url{https://github.com/cran/PMA/blob/master/R/MultiCCA.R}
#
# @param mat.list list of matrices to calculate correlation
# @param ws vector of projection vectors
# @param num.sets number of datasets
#
# @return total correlation
#
GetCors <- function(mat.list, ws, num.sets){
  cors <- 0
  for (i in 2:num.sets) {
    for (j in 1:(i - 1)) {
      thiscor  <-  cor(x = mat.list[[i]] %*% ws[[i]], y = mat.list[[j]] %*% ws[[j]])
      if (is.na(x = thiscor)) {
        thiscor <- 0
      }
      cors <- cors + thiscor
    }
  }
  return(cors)
}

# MultiCCA helper function - calculates critical value (when to stop iterating
# in the while loop)
#
# Modified from PMA package
# @references Witten, Tibshirani, and Hastie, Biostatistics 2009
# @references \url{https://github.com/cran/PMA/blob/master/R/MultiCCA.R}
#
# @param mat.list list of matrices
# @param ws vector of projection vectors
# @param num.sets number of datasets
#
# @return returns updated critical value
#
GetCrit <- function(mat.list, ws, num.sets){
  crit <- 0
  for (i in 2:num.sets) {
    for (j in 1:(i - 1)) {
      crit <- crit + t(ws[[i]]) %*% t(x = mat.list[[i]]) %*% mat.list[[j]] %*% ws[[j]]
    }
  }
  return(crit)
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
  set.seed(seed = seed.use)
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
    pcs.compute = r2.use,
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
L2Norm <- function(vec){
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
#
#
PrepDR <- function(
  object,
  features = NULL
) {
  if (length(x = VariableFeatures(object = object)) == 0 && is.null(x = features)) {
    stop("Variable features haven't been set. Run FindVariableFeatures() or provide a vector of genes names in genes.use and retry.")
  }
  data.use <- GetAssayData(object = object, slot = "scale.data")
  if (nrow(x = data.use ) == 0) {
    stop("Data has not been scaled. Please run ScaleData and retry")
  }
  features <- features %||% VariableFeatures(object = object)
  features.keep <- unique(x = features[features %in% rownames(x = data.use)])
  if (length(x = features.keep) < length(x = features)) {
    features.exclude <- setdiff(x = features, y = features.keep)
    warning(paste0("The following ", length(x = features.exclude), " features requested have not been scaled (running reduction without them): ", paste0(features.exclude, collapse = ", ")))
  }
  features <- features.keep
  features.var <- apply(X = data.use[features, ], MARGIN = 1, FUN = var)
  features.keep <- features[features.var > 0]
  if (length(x = features.keep) < length(x = features)) {
    features.exclude <- setdiff(x = features, y = features.keep)
    warning(paste0("The following ", length(x = features.exclude), " features requested have zero variance (running reduction without them): ", paste0(features.exclude, collapse = ", ")))
  }
  features <- features.keep
  features <- features[!is.na(x = features)]
  data.use <- data.use[features, ]
  return(data.use)
}

# MultiCCA helper function - updates W
#
# Modified from PMA package
# @references Witten, Tibshirani, and Hastie, Biostatistics 2009
# @references \url{https://github.com/cran/PMA/blob/master/R/MultiCCA.R}
#
# @param mat.list list of matrices
# @param i index of current matrix
# @param num.sets number of datasets
# @param ws initial vector of projection vectors
# @param ws.final final vector of projection vectors
#
# @return returns updated w value
#
UpdateW <- function(mat.list, i, num.sets, ws, ws.final){
  tots <- 0
  for (j in (1:num.sets)[-i]) {
    diagmat <- (t(x = ws.final[[i]]) %*% t(x = mat.list[[i]])) %*% (mat.list[[j]] %*% ws.final[[j]])
    diagmat[row(x = diagmat) != col(x = diagmat)] <- 0
    tots <- tots + t(x = mat.list[[i]]) %*% (mat.list[[j]] %*% ws[[j]]) - ws.final[[i]] %*% (diagmat %*% (t(x = ws.final[[j]]) %*% ws[[j]]))
  }
  w <- tots / L2Norm(vec = tots)
  return(w)
}
