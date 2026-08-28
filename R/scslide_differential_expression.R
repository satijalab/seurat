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


#' Perform cor.test between each gene and a response variable given a Seurat object
#'
#' This function take a Seurat object's gene expression matrix and a meta-data column (a response variable), and
#' run cor.test between each gene and the response variable and keep track of the top-correlated genes and return
#' that gene list to user. It is intended to give users a more response-informative gene list for downstream supervised
#' dimension reduction analyses like partial least squares.
#'
#' @param object A Seurat object
#' @param assay The name of the assay to retrieve the expression matrix
#' @param layer The name of the layer to use for cor.test
#' @param group.by A metadata column name to group cells by before cor.test. If NULL,
#'   falls back to cor.test without grouping. Default is NULL.
#' @param Y A metadata column name of responses.
#' @param cor.cutoff The correlation coefficient cutoff to select the top genes. Default is NULL.
#' @param top.prop The proportion of genes to be selected as the top genes, if cor.cutoff is not provided.
#'   Can be used together with cor.cutoff, and they cannot be NULL at the same time. Default is 0.01, i.e. the most correlated 1\% of genes will be selected.
#' @param min.cells the minimum number of cells per group to perform the cor.test. Default is 100. groups with < 100 will be skipped.
#' @param min.pct only test genes that are detected in a minimum fraction of min.pct cells in each group. Default is 0.1
#' @param verbose Print progress bars and output
#' @param unlist whether to return a unlist object (vector) or a named list (each element is a vector of DEGs for a category in group.by)
#' @param ... Additional arguments
#'
#' @return return a list of gene names for the genes that highly correlated with the responses
#' @export
#' @concept scSLIDE_DE
#'
#' @examples
#' \dontrun{
#' # Find genes correlated with pseudotime
#' cor_genes <- QuickCorTest(seurat_obj, Y = "pseudotime")
#' }
#'
#' @importFrom stats model.matrix sd
#' @importFrom SeuratObject FetchData DefaultAssay Assays LayerData Layers CellsByIdentities Idents<- Cells
#'
QuickCorTest <- function(
    object = NULL,
    assay = NULL,
    layer = "data",
    group.by = NULL,
    Y = NULL,
    cor.cutoff = NULL,
    top.prop = 0.01,
    min.pct = 0.1,
    min.cells = 100,
    verbose = TRUE,
    unlist = TRUE,
    ...
){
  assay <- assay %||% DefaultAssay(object = object)
  assay <- match.arg(arg = assay, choices = Assays(object = object))
  # fetch the target variable
  if (is.null(x = Y) | length(x = Y) != 1) {
    stop("When setting 'Y', one and only one variable is allowed.")
  } else {
    Y <- FetchData(
      object = object,
      vars = Y
    )
  }
  for(i in 1:ncol(Y)){
    if(is.character(Y[,i])) Y[,i] <- as.factor(Y[,i])
  }
  # fetch the expression assay
  data.use <- object[[assay]]
  if (length(x = Layers(object = data.use, search = layer)) > 1) {
    stop(layer, " layers are not joined. Please run JoinLayers")
  }
  # get cell IDs within each category of the "group.by"
  if (!is.null(x = group.by) && group.by %in% colnames(object@meta.data)) {
    Idents(object) <- group.by
    cell_celltype_list <- CellsByIdentities(object, cells = Cells(data.use))
  } else {
    cell_celltype_list <- list(all = Cells(data.use))
  }
  # check
  if(is.null(x = cor.cutoff) & is.null(x = top.prop)) stop("Must provide at least one of 'top.prop' and 'cor.cutoff'.")
  #### build each group's response matrix first, so that groups which are too small,
  #### or whose responses turn out to be constant, never reach the expression data
  responses <- lapply(names(cell_celltype_list), function(CELLTYPE) {
    idx <- cell_celltype_list[[CELLTYPE]]
    # Skip if the number of cells is less than min.cells
    if(length(idx) < min.cells){
      if(verbose) message(paste0("Skipping ", CELLTYPE, " since the number of cells is less than ", min.cells))
      return(NULL)  # Return NULL for skipped cell types
    } else {
      if(verbose) message(paste0("Processing ", CELLTYPE, "..."))
    }
    # get Y
    Y2 <- Y[idx, , drop = F]
    Y2 <- model.matrix(~. + 0, data = Y2)
    rm_idx = vector()
    for(i in 1:ncol(Y2)){
      if(length(unique(Y2[, i])) <= 1) rm_idx <- c(rm_idx, i)
    }
    if(length(rm_idx) > 0) Y2 <- Y2[, -rm_idx, drop = F]
    if(ncol(Y2) == 0) return(NULL)
    return(Y2)
  })
  names(responses) <- names(cell_celltype_list)
  tested <- names(responses)[!sapply(responses, is.null)]
  slct_genes_list <- vector(mode = "list", length = length(x = cell_celltype_list))
  names(slct_genes_list) <- names(cell_celltype_list)

  #### read the expression matrix once and derive every group's statistics and
  #### cross-products from that, instead of subsetting cells group by group
  if(length(tested) > 0){
    responses <- responses[tested]
    df_data <- LayerData(object = data.use, assay = assay, layer = layer)
    index <- lapply(cell_celltype_list[tested], match, table = colnames(df_data))
    groups <- factor(rep_len(NA_character_, ncol(df_data)), levels = tested)
    for(CELLTYPE in tested) groups[index[[CELLTYPE]]] <- CELLTYPE
    # cells outside every tested group still need a level to be counted under
    if(anyNA(groups)){
      levels(groups) <- c(tested, ".unused")
      groups[is.na(groups)] <- ".unused"
    }
    # Detection rate, mean and SD of every gene in every group, in a single pass
    group.stats <- GroupStats(x = df_data, groups = groups)
    # Cross-products of every gene against every group's responses, in a second pass
    xy <- GroupCrossprod(x = df_data, responses = responses, index = index)

    ####  loop through the cell types to get the top DEGs (separately for each column in Y)
    slct_genes_list[tested] <- lapply(tested, function(CELLTYPE) {
      Y2 <- responses[[CELLTYPE]]
      n <- nrow(x = Y2)
      idx_rowCount <- group.stats$nonzero[, CELLTYPE] >= min.pct*n
      features <- rownames(df_data)[idx_rowCount]
      # Calculate correlation more efficiently
      cor_res <- CorFromCrossprod(
        xy = xy[[CELLTYPE]][idx_rowCount, , drop = FALSE],
        n = n,
        x.mean = group.stats$mean[idx_rowCount, CELLTYPE],
        x.sd = group.stats$sd[idx_rowCount, CELLTYPE],
        y = Y2
      )
      # a single response collapses to the named vector that apply() used to return
      if (nrow(x = cor_res) == 1) {
        cor_res <- cor_res[1, ]
      }

      # Handle both vector (continuous Y) and matrix (categorical Y) cases
      if (is.vector(cor_res)) {
        # For continuous variables, cor_res is a vector
        n_top <- ceiling(top.prop * length(cor_res))
        abs_cor <- abs(cor_res)
        idx_slct <- order(abs_cor, decreasing = TRUE)[seq_len(length.out = n_top)]
        if(!is.null(x = cor.cutoff)){
          cutoff_idx <- which(abs_cor >= cor.cutoff)
          idx_slct <- intersect(idx_slct, cutoff_idx)
        }
      } else {
        # For categorical variables, cor_res is a matrix
        n_top <- ceiling(top.prop * ncol(cor_res))
        idx_slct <- apply(abs(cor_res), 1, function(x) {
          a <- order(x, decreasing = TRUE)[seq_len(length.out = n_top)]
          if(!is.null(x = cor.cutoff)){
            b <- which(x >= cor.cutoff)
            return(intersect(a, b))
          } else {
            return(a)
          }
        })
        idx_slct <- sort(unique(as.vector(unlist(idx_slct))))
      }

      # Return the selected genes
      features[idx_slct]
    })
  }
  # Remove NULL entries (skipped cell types)
  if(isTRUE(x = unlist)){
    slct_genes_list <- slct_genes_list[!sapply(slct_genes_list, is.null)]
    all_genes <- Reduce(union, slct_genes_list)
    return(all_genes)
  } else {
    return(slct_genes_list)
  }
}

# Helper function: nonzero count, mean and standard deviation of every gene within
# every group, as features x groups matrices. BPCells::pseudobulk_matrix() returns
# all three for all groups from a single streamed pass, so an on-disk matrix is read
# once no matter how many groups there are -- subsetting its cells group by group
# would instead cost one random-access pass each.
GroupStats <- function(x, groups) {
  if (inherits(x = x, what = 'IterableMatrix')) {
    pb <- BPCells::pseudobulk_matrix(
      mat = x,
      cell_groups = groups,
      method = c('nonzeros', 'mean', 'variance')
    )
    return(list(nonzero = pb$nonzeros, mean = pb$mean, sd = sqrt(x = pb$variance)))
  }
  n <- tabulate(bin = groups, nbins = nlevels(x = groups))
  indicator <- Matrix::sparseMatrix(
    i = seq_along(along.with = groups),
    j = as.integer(x = groups),
    x = 1,
    dims = c(length(x = groups), nlevels(x = groups)),
    dimnames = list(colnames(x = x), levels(x = groups))
  )
  # `x > 0` builds a logical copy of the whole matrix that the product below then
  # coerces straight back to double. Reusing the sparsity pattern instead skips both
  # steps; it counts stored entries rather than positive ones, which is the same
  # thing for the counts and data layers this runs on
  if (inherits(x = x, what = 'dgCMatrix')) {
    pattern <- x
    pattern@x <- rep.int(x = 1, times = length(x = x@x))
  } else {
    pattern <- x > 0
  }
  group.mean <- t(x = t(x = as.matrix(x = x %*% indicator)) / n)
  group.sumsq <- as.matrix(x = (x ^ 2) %*% indicator)
  group.var <- t(x = (t(x = group.sumsq) - n * t(x = group.mean) ^ 2) / (n - 1))
  return(list(
    nonzero = as.matrix(x = pattern %*% indicator),
    mean = group.mean,
    sd = sqrt(x = pmax(group.var, 0))
  ))
}

# Helper function: sum_i x_gi * y_ik over each group's cells, for every gene g,
# response k and group. Every group's responses are packed side by side into one
# right-hand side, zero outside that group's cells, so a single x %*% y covers all
# of them in one pass. Returns one features x responses matrix per group.
#
# Since each cell carries responses for its own group only, that packed right-hand
# side is mostly structural zeros once there is more than one group. Against a
# sparse x a dense right-hand side therefore costs nnz * total responses, growing
# linearly with the number of groups, while a sparse one stays flat (11.6x quicker
# at 40 groups). BPCells is the other way round -- it answers sparse products with
# a lazy MatrixMultiply but has an eager dense kernel -- so on-disk matrices, and
# dense in-memory ones, keep the dense right-hand side.
GroupCrossprod <- function(x, responses, index, max.bytes = 2.5e8) {
  n <- ncol(x = x)
  widths <- vapply(X = responses, FUN = ncol, FUN.VALUE = integer(length = 1L))
  sparse <- length(x = responses) > 1L && inherits(x = x, what = 'sparseMatrix')
  # a sparse right-hand side holds only the responses themselves, so every group
  # fits in one batch; a dense one is batched to keep its size bounded
  batch.size <- if (sparse) sum(widths) else max(1, floor(x = max.bytes / (8 * n)))
  res <- vector(mode = 'list', length = length(x = responses))
  names(x = res) <- names(x = responses)
  for (batch in split(x = seq_along(along.with = responses), f = ceiling(x = cumsum(x = widths) / batch.size))) {
    offsets <- cumsum(x = c(0, widths[batch]))
    if (isTRUE(x = sparse)) {
      rhs <- Matrix::sparseMatrix(
        i = unlist(x = lapply(
          X = seq_along(along.with = batch),
          FUN = function(i) rep(x = index[[batch[i]]], times = widths[batch[i]])
        )),
        j = unlist(x = lapply(
          X = seq_along(along.with = batch),
          FUN = function(i) rep(
            x = offsets[i] + seq_len(length.out = widths[batch[i]]),
            each = length(x = index[[batch[i]]])
          )
        )),
        x = unlist(x = lapply(X = responses[batch], FUN = as.vector)),
        dims = c(n, sum(widths[batch]))
      )
    } else {
      rhs <- matrix(data = 0, nrow = n, ncol = sum(widths[batch]))
      for (i in seq_along(along.with = batch)) {
        rhs[index[[batch[i]]], offsets[i] + seq_len(length.out = widths[batch[i]])] <- responses[[batch[i]]]
      }
    }
    product <- as.matrix(x = x %*% rhs)
    for (i in seq_along(along.with = batch)) {
      res[[batch[i]]] <- product[, offsets[i] + seq_len(length.out = widths[batch[i]]), drop = FALSE]
    }
  }
  return(res)
}

# Helper function: Pearson correlations from precomputed cross-products. Given
# xy[g, k] = sum_i x_gi * y_ik over the group's n cells, the correlation follows in
# closed form from the gene and response means/SDs, with no further access to the
# expression data. Returns a responses x features matrix, matching the orientation
# apply(x, 1, cor, y) produced.
CorFromCrossprod <- function(xy, n, x.mean, x.sd, y) {
  y.mean <- colMeans(x = y)
  y.sd <- apply(X = y, MARGIN = 2, FUN = sd)
  # cov(x_i, y_j) = (sum_k x_ik * y_kj - n * mean(x_i) * mean(y_j)) / (n - 1)
  covariance <- (xy - n * outer(X = x.mean, Y = y.mean)) / (n - 1)
  res <- t(x = covariance / outer(X = x.sd, Y = y.sd))
  # constant genes have no correlation, as cor() would also report
  res[!is.finite(x = res)] <- NA
  return(res)
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

