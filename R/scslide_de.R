#' Trajectory Differential Expression Test
#'
#' Perform trajectory-based differential expression analysis using negative binomial regression.
#' This function tests for genes that change expression along a continuous trajectory variable.
#'
#' @param object Expression data matrix, Assay, StdAssay, or Seurat object
#' @param traj.var For the Seurat method: a single character string naming a column in the object
#'   metadata. For the Assay/default methods: a data frame (one column) containing the trajectory
#'   variable of the samples to be tested against.
#' @param latent.vars a data frame containing the latent variables (e.g., covariates that might affect the gene expression) to include in regression.
#' @param features Genes to test. Default is to use all genes (after QC)
#' @param fc.results fc.results calculated by FoldChange(); if not null QC will be performed based on this
#' @param verbose Print a progress bar once expression testing begins
#' @param assay the assay to be used (for Seurat objects)
#' @param layer the data layer to be used for trajectory DE test. Currently only 'counts' is supported.
#' @param samples the cells/samples to be included in the DE test
#' @param logfc.threshold Limit testing to genes which show, on average, at least X-fold difference (log-scale)
#' between the top and bottom groups of samples (see prob.break.point for details). Default is 0 (i.e., no filtering).
#' @param prob.break.point a numeric vector of probability break points with values in (0, 1). It will be used to calculate 2 quantiles
#' along the trajectory and to calculate the logfc.
#' @param min.pct only test genes that are detected in a minimum proportion of cells (passed as
#'   \code{min.prop} to the internal \code{filterByExpr()} filter). Default is 0.1.
#' @param min.count Minimum count threshold for gene filtering
#' @param pseudocount.use Pseudocount to add to averaged expression values when calculating logFC. 1 by default.
#' @param complete.results Whether to return complete results or summary
#' @param ... additional parameters to be passed to glmGamPoi::glm_gp or filterByExpr
#'
#' @return A data frame with differential expression results
#' @export
#' @concept scSLIDE_DE
#'
#' @examples
#' \dontrun{
#' # Run trajectory DE test on a Seurat object
#' de_results <- TrajDETest(seurat_obj, traj.var = "pseudotime")
#' }
#'
#' @importFrom stats model.matrix quantile hat median
#' @importFrom SeuratObject FetchData DefaultAssay Assays LayerData Layers
#'

TrajDETest <- function(object, ...) {
  UseMethod(generic = 'TrajDETest', object = object)
}

#' @rdname TrajDETest
#' @method TrajDETest default
#' @export
TrajDETest.default <- function(
    object,
    traj.var = NULL,
    latent.vars = NULL,
    features = NULL,
    fc.results = NULL,
    verbose = TRUE,
    ...
){
  if (!requireNamespace("glmGamPoi", quietly = TRUE)) {
    stop(
      "The 'glmGamPoi' package is required for TrajDETest() but is not installed.\n\n",
      "To install glmGamPoi, please run:\n",
      "  if (!require('BiocManager', quietly = TRUE))\n",
      "    install.packages('BiocManager')\n",
      "  BiocManager::install('glmGamPoi')\n"
    )
  }
  # first get the covariates model.matrix
  if(!is.null(x = latent.vars)){
    mat_all <- cbind(traj.var, latent.vars)
    colnames(mat_all)[1] <- "Traj"
    full_model_mat <- model.matrix(~., data = mat_all)
    # contrast_model_mat <- full_model_mat[,-2, drop = F]
  } else {
    mat_all <- traj.var
    colnames(mat_all)[1] <- "Traj"
    full_model_mat <- model.matrix(~., data = mat_all)
    # contrast_model_mat <- full_model_mat[,-2, drop = F]
  }

  # perform negative binomial regression
  if(is.null(x = features)) {
    features <- rownames(object)
  }
  fit <- glmGamPoi::glm_gp(data = object[features, ],
                design = full_model_mat,
                size_factors = F,
                on_disk = FALSE,
                verbose = verbose)
  # after the fitting, here we will calculate the DE p-val for the traj.var
  beta = fit$Beta
  p = glmGamPoi::test_de(fit = fit,
              contrast = "Traj",
              verbose = verbose)
  names(p)[2] = "p_Traj"
  beta = as.data.frame(beta)
  names(beta) = paste0("beta_", names(beta))
  p$gene_ID = p$name
  res = cbind(beta, p)

  # append the fc.results
  if(!is.null(x = fc.results)){
    fc.results <- fc.results[features, ]
    res <- cbind(res, fc.results)
  }
  #
  return(res)
}

#' @rdname TrajDETest
#' @method TrajDETest Assay
#' @export
TrajDETest.Assay <- function(
    object,
    layer = "counts",
    traj.var = NULL,
    latent.vars = NULL,
    samples = NULL,
    features = NULL,
    logfc.threshold = 0,
    prob.break.point = c(1/3, 2/3),
    min.pct = 0.1,
    min.count = 10,
    pseudocount.use = 1,
    verbose = TRUE,
    ...
){
  layer <- "counts"
  if (length(x = Layers(object = object, search = layer)) > 1) {
    stop(layer, " layers are not joined. Please run JoinLayers")
  }

  # get the expression matrix
  data.use <- LayerData(object = object, layer = layer, cells = samples, features = features)

  # use filterByExpr to remove very lowly-expressed genes
  idx_for_DE <- filterByExpr(data.use, min.count = min.count, min.prop = min.pct)
  idx_for_DE <- names(idx_for_DE)[which(idx_for_DE)]

  # calculate the log Fold change between the top and bottom groups along the trajectory
  if(length(prob.break.point) != 2){
    stop("prob.break.point must be a numeric vector of length 2.")
  }
  # get the top and bottom groups
  quantiles <- quantile(traj.var[, 1], probs = prob.break.point)
  idx_group_bottom <- rownames(traj.var)[traj.var[, 1] <= quantiles[1]]
  idx_group_top <- rownames(traj.var)[traj.var[, 1] > quantiles[2]]
  # calculate the logFC
  fc.results <- FoldChange(object = LayerData(object = object, layer = 'data'),
                           cells.1 = idx_group_top,
                           cells.2 = idx_group_bottom,
                           mean.fxn = function(x) log(x = (SeuratObject::rowSums(x = expm1(x = x)) + pseudocount.use)/NCOL(x), base = 2),
                           fc.name = "avg_log2FC",
                           features = idx_for_DE)
  # update idx_for_DE
  if(logfc.threshold != 0){
    message("Further removing genes based on logfc.threshold =", logfc.threshold)
    idx_for_DE <- rownames(fc.results)[abs(fc.results$avg_log2FC) >= logfc.threshold]
  }

  # pass it to the next
  de.results <- TrajDETest(
    object = data.use,
    traj.var = traj.var,
    latent.vars = latent.vars,
    features = idx_for_DE,
    fc.results = fc.results,
    verbose = verbose,
    ...
  )
  #
  return(de.results)
}

#' @rdname TrajDETest
#' @method TrajDETest StdAssay
#' @export
TrajDETest.StdAssay <- TrajDETest.Assay

#' @rdname TrajDETest
#' @method TrajDETest Seurat
#' @export
TrajDETest.Seurat <- function(
    object,
    assay = NULL,
    layer = "counts",
    traj.var = NULL,
    latent.vars = NULL,
    samples = NULL,
    features = NULL,
    logfc.threshold = 0,
    prob.break.point = c(1/3, 2/3),
    min.pct = 0.1,
    min.count = 10,
    pseudocount.use = 1,
    complete.results = FALSE,
    verbose = TRUE,
    ...
) {
  #
  assay <- assay %||% DefaultAssay(object = object)
  assay <- match.arg(arg = assay, choices = Assays(object = object))
  message("Running DE using the ", assay, " assay...")

  # get cells/samples and features
  features <- features %||% rownames(x = object)
  samples <- samples %||% colnames(x = object)

  # fetch the target variable
  if (is.null(x = traj.var) | length(x = traj.var) != 1) {
    stop("When setting 'traj.var', one and only one variable is allowed.")
  } else {
    traj.var <- FetchData(
      object = object,
      vars = traj.var,
      cells = samples
    )
  }

  # fetch latent.vars
  if (!is.null(x = latent.vars)) {
    latent.vars <- FetchData(
      object = object,
      vars = latent.vars,
      cells = samples
    )
  }

  # fetch the expression assay
  data.use <- object[[assay]]

  # pass it to the next
  de.results <- TrajDETest(
    object = data.use,
    layer = layer,
    samples = samples,
    traj.var = traj.var,
    latent.vars = latent.vars,
    features = features,
    logfc.threshold = logfc.threshold,
    prob.break.point = prob.break.point,
    min.pct = min.pct,
    min.count = min.count,
    pseudocount.use = pseudocount.use,
    verbose = verbose,
    ...
  )
  #
  de.results <- de.results[order(de.results$p_Traj), ]
  if(isTRUE(x = complete.results)){
    return(de.results)
  } else {
    return(de.results[, c("gene_ID", "beta_Traj", "p_Traj", "adj_pval", "avg_log2FC", "pct.1", "pct.2")])
  }
}

# Helper function: simplified cpm calculation
cpm <- function(y, lib.size=NULL, log=FALSE, prior.count=2) {
    if(is.null(lib.size)) lib.size <- colSums(y)
    y <- as.matrix(y)
    
    if(log) {
        t(log2(t(y + prior.count)/lib.size * 1e6))
    } else {
        t(t(y)/lib.size * 1e6)
    }
}

# Simplified filterByExpr function (internal)
filterByExpr <- function(y, design=NULL, group=NULL, lib.size=NULL, min.count=10, min.total.count=15, large.n=10, min.prop=0.7) {
    y <- as.matrix(y)
    if(mode(y) != "numeric") stop("y is not a numeric matrix")
    if(is.null(lib.size)) lib.size <- colSums(y)
    
    # Minimum effective sample size for any of the coefficients
    if(is.null(group)) {
        if(is.null(design)) {
            message("No group or design set. Assuming all samples belong to one group.")
            MinSampleSize <- ncol(y)
        } else {
            h <- hat(design)
            MinSampleSize <- 1/max(h)
        }
    } else {
        group <- as.factor(group)
        n <- tabulate(group)
        MinSampleSize <- min(n[n > 0L])
    }
    if(MinSampleSize > large.n) MinSampleSize <- large.n + (MinSampleSize-large.n)*min.prop
    
    # CPM cutoff
    MedianLibSize <- median(lib.size)
    CPM.Cutoff <- min.count/MedianLibSize*1e6
    CPM <- cpm(y, lib.size=lib.size)
    tol <- 1e-14
    keep.CPM <- rowSums(CPM >= CPM.Cutoff) >= (MinSampleSize - tol)
    
    # Total count cutoff
    keep.TotalCount <- (rowSums(y) >= min.total.count - tol)
    
    keep.CPM & keep.TotalCount
}

