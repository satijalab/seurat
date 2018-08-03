#' @include package.R
#'
NULL

#' Load in data from 10X
#'
#' Enables easy loading of sparse data matrices provided by 10X genomics.
#'
#' @param data.dir Directory containing the matrix.mtx, genes.tsv, and barcodes.tsv
#' files provided by 10X. A vector or named vector can be given in order to load
#' several data directories. If a named vector is given, the cell barcode names
#' will be prefixed with the name.
#'
#' @return Returns a sparse matrix with rows and columns labeled
#'
#' @importFrom Matrix readMM
#'
#' @export
#'
#' @examples
#' \dontrun{
#' data_dir <- 'path/to/data/directory'
#' list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
#' expression_matrix <- Read10X(data.dir = data_dir)
#' seurat_object = CreateSeuratObject(counts = expression_matrix)
#' }
#'
Read10X <- function(data.dir = NULL){
  full.data <- list()
  for (i in seq_along(data.dir)) {
    run <- data.dir[i]
    if (! dir.exists(run)){
      stop("Directory provided does not exist")
    }
    if(!grepl("\\/$", run)){
      run <- paste(run, "/", sep = "")
    }
    barcode.loc <- paste0(run, "barcodes.tsv")
    gene.loc <- paste0(run, "genes.tsv")
    matrix.loc <- paste0(run, "matrix.mtx")
    if (!file.exists(barcode.loc)){
      stop("Barcode file missing")
    }
    if (! file.exists(gene.loc)){
      stop("Gene name file missing")
    }
    if (! file.exists(matrix.loc)){
      stop("Expression matrix file missing")
    }
    data <- readMM(file = matrix.loc)
    cell.names <- readLines(barcode.loc)
    gene.names <- readLines(gene.loc)
    if (all(grepl(pattern = "\\-1$", x = cell.names))) {
      cell.names <- as.vector(
        x = as.character(
          x = sapply(
            X = cell.names,
            FUN = ExtractField, field = 1,
            delim = "-"
          )
        )
      )
    }
    rownames(x = data) <- make.unique(
      names = as.character(
        x = sapply(
          X = gene.names,
          FUN = ExtractField,
          field = 2,
          delim = "\\t"
        )
      )
    )
    if (is.null(x = names(x = data.dir))) {
      if(i < 2){
        colnames(x = data) <- cell.names
      }
      else {
        colnames(x = data) <- paste0(i, "_", cell.names)
      }
    } else {
      colnames(x = data) <- paste0(names(x = data.dir)[i], "_", cell.names)
    }
    full.data <- append(x = full.data, values = data)
  }
  full.data <- do.call(cbind, full.data)
  return(full.data)
}

#' Read 10X hdf5 file
#'
#' Read gene expression matrix from 10X CellRanger hdf5 file
#'
#' @param filename Path to h5 file
#' @param ensg.names Label row names with ENSG names rather than unique gene names
#'
#' @return Returns a sparse matrix with rows and columns labeled. If multiple genomes are present,
#' returns a list of sparse matrices (one per genome).
#'
#' @importFrom hdf5r H5File
#'
#' @export
#'
Read10X_h5 <- function(filename, ensg.names = FALSE){
  if(!file.exists(filename)){
    stop("File not found")
  }
  infile <- H5File$new(filename)
  genomes <- names(infile)
  output <- list()
  for(genome in genomes){
    counts <- infile[[paste0(genome, '/data')]]
    indices <- infile[[paste0(genome, '/indices')]]
    indptr <- infile[[paste0(genome, '/indptr')]]
    shp <- infile[[paste0(genome, '/shape')]]
    if(ensg.names){
      gene_names <- infile[[paste0(genome, '/genes')]][]
    } else {
      gene_names <- make.unique(infile[[paste0(genome, '/gene_names')]][])
    }
    barcodes <- infile[[paste0(genome, '/barcodes')]]
    sparse.mat <- sparseMatrix(i = indices[] + 1, p = indptr[],
                               x = as.numeric(counts[]),
                               dims = shp[], giveCsparse = FALSE)
    rownames(sparse.mat) <- gene_names
    colnames(sparse.mat) <- barcodes[]
    output[[genome]] <- sparse.mat
  }
  infile$close_all()
  if(length(output) == 1) {
    return(output[[genome]])
  } else{
    return(output)
  }
}

#' @export
#'
NormalizeData.default <- function(
  object,
  normalization.method = "LogNormalize",
  scale.factor = 1e4,
  verbose = TRUE
) {
  if (is.null(x = normalization.method)) {
    return(object)
  }
  normalized.data <- switch(
    EXPR = normalization.method,
    'LogNormalize' = LogNormalize(
      data = object,
      scale.factor = scale.factor,
      verbose = verbose
    ),
    'CLR' = CustomNormalize(
      data = object,
      custom_function = function(x) {
        return(log1p(x = x / (exp(x = sum(log1p(x = x[x > 0]), na.rm = TRUE) / length(x = x + 1)))))
      },
      across = 'features'
    ),
    stop("Unkown normalization method: ", normalization.method)
  )
  return(normalized.data)
}

#' @describeIn NormalizeData Normalize data in an Assay object
#' @export
#' @method NormalizeData Assay
#'
NormalizeData.Assay <- function(
  object,
  normalization.method = "LogNormalize",
  scale.factor = 1e4,
  verbose = TRUE
) {
  object <- SetAssayData(
    object = object,
    slot = 'data',
    new.data = NormalizeData(
      object = GetAssayData(object = object, slot = 'counts'),
      normalization.method = normalization.method,
      scale.factor = scale.factor,
      verbose = verbose
    )
  )
  return(object)
}

#' @param assay.use Name of assay to use
#'
#' @describeIn NormalizeData Normalize data in a Seurat object
#' @export
#' @method NormalizeData Seurat
#'
#' @examples
#' \dontrun{
#' pbmc_small
#' pmbc_small <- NormalizeData(object = pbmc_small)
#' }
#'
NormalizeData.Seurat <- function(
  object,
  assay.use = NULL,
  normalization.method = "LogNormalize",
  scale.factor = 1e4,
  verbose = TRUE,
  workflow.name = NULL
) {
  assay.use <- assay.use %||% DefaultAssay(object = object)
  if (!is.null(workflow.name)) {
    object <- PrepareWorkflow(object = object, workflow.name = workflow.name)
  }
  assay.data <- GetAssay(object = object, assay.use = assay.use)
  assay.data <- NormalizeData(
    object = assay.data,
    normalization.method = normalization.method,
    scale.factor = scale.factor,
    verbose = verbose
  )
  object[[assay.use]] <- assay.data
  object <- LogSeuratCommand(object = object)
  if (!is.null(workflow.name)) {
    object <- UpdateWorkflow(object = object, workflow.name = workflow.name)
  }
  return(object)
}

#' Use regularized negative binomial regression to normalize UMI count data
#'
#' This function calls sctransform::vst. The sctransform package is available at
#' https://github.com/ChristophH/sctransform.
#' Use this function as an alternative to the NormalizeData, FindVariableFeatures, ScaleData workflow.
#' Results are saved in the assay's data and scale.data slot, and sctransform::vst ntermediate
#' results are saved in misc slot of seurat object.
#'
#' @param object A seurat object
#' @param assay.use Name of assay to use
#' @param do.correct.umi Place corrected UMI matrix in assay data slot
#' @param variable.features.zscore Z-score threshold for calling features highly variable;
#' z-scores are based on variances of regression model pearson residuals of all features
#' @param variable.features.n Use this many features as variable features after ranking by variance
#' @param do.scale Whether to scale residuals to have unit variance
#' @param do.center Whether to center residuals to have mean zero
#' @param scale.max Max value after scaling and/or centering
#' @param verbose Whether to print messages and progress bars
#' @param ... Additional parameters passed to \code{sctransform::vst}
#'
#' @export
RegressRegNB <- function(
  object,
  assay.use = NULL,
  do.correct.umi = FALSE,
  variable.features.zscore = 1,
  variable.features.n = NULL,
  do.scale = FALSE,
  do.center = FALSE,
  scale.max = .Machine$double.xmax,
  verbose = TRUE,
  ...
) {
  if (!requireNamespace('sctransform')) {
    stop('Install sctransform package from https://github.com/ChristophH/sctransform to use regularized negative binomial regression models.')
  }
  assay.use <- assay.use %||% DefaultAssay(object = object)
  assay.obj <- GetAssay(object = object, assay.use = assay.use)
  umi <- GetAssayData(object = assay.obj, slot = 'counts')

  vst.out <- sctransform::vst(umi, show_progress = verbose, return_cell_attr = TRUE, ...)
  # cell_attr = NULL,
  # latent_var = c('log_umi_per_gene'),
  # batch_var = NULL,
  # n_genes = 2000,
  # n_cells = NULL,
  # method = 'poisson',
  # res_clip_range = c(-50, 50),
  # bin_size = 256,
  # min_cells = 5,
  # return_cell_attr = FALSE,
  # return_gene_attr = FALSE)

  # put corrected umi counts in data slot
  if (do.correct.umi) {
    if (verbose) {
      message('Calculate corrected UMI matrix and place in data slot')
    }
    umi.corrected <- sctransform::denoise(x = vst.out)
    umi.corrected <- as(umi.corrected, 'dgCMatrix')
    # skip SetAssayData.Assay because of restrictive dimension checks there
    slot(object = assay.obj, name = 'data') <- umi.corrected
    # assay.obj <- SetAssayData(
    #   object = assay.obj,
    #   slot = 'data',
    #   new.data = umi.corrected
    # )
  }

  # set variable features
  if (verbose) {
    message('Determine variable features')
  }
  if ('residual_variance' %in% names(vst.out$gene_attr)) {
    feature.variance <- setNames(vst.out$gene_attr$residual_variance, rownames(vst.out$gene_attr))
  } else {
    feature.variance <- RowVar(vst.out$y)
    names(feature.variance) <- rownames(vst.out$y)
  }
  feature.variance <- sort(x = feature.variance, decreasing = TRUE)
  if (!is.null(variable.features.n)) {
    top.features <- names(feature.variance)[1:variable.features.n]
  } else {
    feature.variance <- scale(feature.variance)[, 1]
    top.features <- names(feature.variance)[feature.variance > variable.features.zscore]
  }
  VariableFeatures(object = assay.obj) <- top.features
  if (verbose) {
    message('Set ', length(top.features), ' variable features')
  }

  scale.data <- vst.out$y
  # re-scale the residuals
  if (do.scale || do.center) {
    if (verbose) {
      message('Re-scale residuals')
    }
    scale.data <- FastRowScale(
      mat = vst.out$y,
      scale = do.scale,
      center = do.center,
      scale_max = scale.max,
      display_progress = FALSE
    )
    dimnames(scale.data) <- dimnames(vst.out$y)
  }

  assay.obj <- SetAssayData(
    object = assay.obj,
    slot = 'scale.data',
    new.data = scale.data
  )
  object[[assay.use]] <- assay.obj

  # save vst output (except y) in @misc slot
  vst.out$y <- NULL
  object@misc[['vst.out']] <- vst.out
  return(object)
}

#' @export
#'
ScaleData.default <- function(
  object,
  features.use = NULL,
  vars.to.regress = NULL,
  latent.data = NULL,
  model.use = 'linear',
  use.umi = FALSE,
  do.scale = TRUE,
  do.center = TRUE,
  scale.max = 10,
  block.size = 1000,
  min.cells.to.block = 3000,
  verbose = TRUE,
  ...
) {
  features.use <- features.use %||% rownames(x = object)
  features.use <- as.vector(x = intersect(x = features.use, y = rownames(x = object)))
  object <- object[features.use, , drop = FALSE]
  scaled.data <- matrix(data = NA, nrow = nrow(x = object), ncol = ncol(x = object))
  dimnames(x = scaled.data) <- dimnames(x = object)
  min.cells.to.block <- min(min.cells.to.block, ncol(x = object))
  gc(verbose = FALSE)
  if (!is.null(x = vars.to.regress)) {
    if (is.null(x = latent.data)) {
      latent.data <- data.frame(row.names = colnames(x = object))
    } else {
      latent.data <- latent.data[colnames(x = object), , drop = FALSE]
      rownames(x = latent.data) <- colnames(x = object)
    }
    if (any(vars.to.regress %in% rownames(x = object))) {
      latent.data <- cbind(
        latent.data,
        t(x = object[vars.to.regress[vars.to.regress %in% rownames(x = object)], ])
      )
    }
    # Currently, RegressOutMatrix will do nothing if latent.data = NULL
    # TODO: implement parallelization for RegressOutMatrix
    if (verbose) {
      message("Regressing out ", paste(vars.to.regress, collapse = ', '))
    }
    object <- RegressOutMatrix(
      data.expr = object,
      latent.data = latent.data,
      features.regress = features.use,
      model.use = model.use,
      use.umi = use.umi,
      verbose = verbose
    )
    gc(verbose = FALSE)
  }
  max.block <- ceiling(x = length(x = features.use) / block.size)
  if (verbose) {
    message("Scaling data matrix")
    pb <- txtProgressBar(min = 0, max = max.block, style = 3, file = stderr())
  }
  for (i in 1:max.block) {
    my.inds <- ((block.size * (i - 1)):(block.size * i - 1)) + 1
    my.inds <- my.inds[my.inds <= length(x = features.use)]
    if (inherits(x = object, what = c('dgCMatrix', 'dgTMatrix'))) {
      scale.function <- FastSparseRowScale
    } else {
      object <- as.matrix(x = object)
      scale.function <- FastRowScale
    }
    data.scale <- scale.function(
      mat = object[features.use[my.inds], , drop = FALSE],
      scale = do.scale,
      center = do.center,
      scale_max = scale.max,
      display_progress = FALSE
    )
    dimnames(x = data.scale) <- dimnames(x = object[features.use[my.inds], ])
    scaled.data[features.use[my.inds], ] <- data.scale
    rm(data.scale)
    gc(verbose = FALSE)
    if (verbose) {
      setTxtProgressBar(pb = pb, value = i)
    }
  }
  if (verbose) {
    close(con = pb)
  }
  scaled.data[is.na(x = scaled.data)] <- 0
  gc(verbose = FALSE)
  return(scaled.data)
}

#' @param latent.data Extra data to regress out, should be cells x latent data
#'
#' @describeIn ScaleData Scale an Assay object
#' @export
#' @method ScaleData Assay
#'
ScaleData.Assay <- function(
  object,
  features.use = NULL,
  vars.to.regress = NULL,
  latent.data = NULL,
  model.use = 'linear',
  use.umi = FALSE,
  do.scale = TRUE,
  do.center = TRUE,
  scale.max = 10,
  block.size = 1000,
  min.cells.to.block = 3000,
  verbose = TRUE,
  ...
) {
  use.umi <- ifelse(test = model.use != 'linear', yes = TRUE, no = use.umi)
  slot.use <- ifelse(test = use.umi, yes = 'counts', no = 'data')
  features.use <- features.use %||% VariableFeatures(object)
  if (length(features.use) == 0) features.use <- rownames(GetAssayData(object = object, slot = slot.use))
  object <- SetAssayData(
    object = object,
    slot = 'scale.data',
    new.data = ScaleData(
      object = GetAssayData(object = object, slot = slot.use),
      features.use = features.use,
      vars.to.regress = vars.to.regress,
      latent.data = latent.data,
      model.use = model.use,
      use.umi = use.umi,
      do.scale = do.scale,
      do.center = do.center,
      scale.max = scale.max,
      block.size = block.size,
      min.cells.to.block = min.cells.to.block,
      verbose = verbose,
      ...
    )
  )
  return(object)
}

#' @param assay.use Name of Assay to scale
#' @param workflow.name Name of workflow
#'
#' @describeIn ScaleData Scale a Seurat object
#' @export
#' @method ScaleData Seurat
#'
ScaleData.Seurat <- function(
  object,
  features.use = NULL,
  assay.use = NULL,
  vars.to.regress = NULL,
  model.use = 'linear',
  use.umi = FALSE,
  do.scale = TRUE,
  do.center = TRUE,
  scale.max = 10,
  block.size = 1000,
  min.cells.to.block = 3000,
  verbose = TRUE,
  workflow.name = NULL,
  ...
) {
  if (!is.null(workflow.name)) {
    object <- PrepareWorkflow(object = object, workflow.name = workflow.name)
  }
  assay.use <- assay.use %||% DefaultAssay(object = object)
  assay.data <- GetAssay(object = object, assay.use = assay.use)
  if (any(vars.to.regress %in% colnames(x = object[]))) {
    latent.data <- object[vars.to.regress[vars.to.regress %in% colnames(x = object[])]]
  } else {
    latent.data <- NULL
  }
  assay.data <- ScaleData(
    object = assay.data,
    features.use = features.use,
    vars.to.regress = vars.to.regress,
    latent.data = latent.data,
    model.use = model.use,
    use.umi = use.umi,
    do.scale = do.scale,
    do.center = do.center,
    scale.max = scale.max,
    block.size = block.size,
    min.cells.to.block = min.cells.to.block,
    verbose = verbose,
    ...
  )
  object[[assay.use]] <- assay.data
  object <- LogSeuratCommand(object = object)
  if (!is.null(workflow.name)) {
    object <- UpdateWorkflow(object = object, workflow.name = workflow.name)
  }
  return(object)
}

#' Normalize raw data
#'
#' Normalize count data per cell and transform to log scale
#'
#' @param data Matrix with the raw count data
#' @param scale.factor Scale the data. Default is 1e4
#' @param verbose Print progress
#'
#' @return Returns a matrix with the normalize and log transformed data
#'
#' @import Matrix
#' @importFrom methods as
#'
#' @export
#'
#' @examples
#' mat <- matrix(data = rbinom(n = 25, size = 5, prob = 0.2), nrow = 5)
#' mat
#' mat_norm <- LogNormalize(data = mat)
#' mat_norm
#'
LogNormalize <- function(data, scale.factor = 1e4, verbose = TRUE) {
  if (class(x = data) == "data.frame") {
    data <- as.matrix(x = data)
  }
  if (class(x = data) != "dgCMatrix") {
    data <- as(object = data, Class = "dgCMatrix")
  }
  # call Rcpp function to normalize
  if (verbose) {
    cat("Performing log-normalization\n", file = stderr())
  }
  norm.data <- LogNorm(data, scale_factor = scale.factor, display_progress = verbose)
  colnames(x = norm.data) <- colnames(x = data)
  rownames(x = norm.data) <- rownames(x = data)
  return(norm.data)
}

#' Sample UMI
#'
#' Downsample each cell to a specified number of UMIs. Includes
#' an option to upsample cells below specified UMI as well.
#'
#' @param data Matrix with the raw count data
#' @param max.umi Number of UMIs to sample to
#' @param upsample Upsamples all cells with fewer than max.umi
#' @param progress.bar Display the progress bar
#'
#' @import Matrix
#' @importFrom methods as
#'
#' @return Matrix with downsampled data
#'
#' @export
#'
#' @examples
#' raw_data = as.matrix(x = pbmc_small@raw.data)
#' downsampled = SampleUMI(data = raw_data)
#' head(x = downsampled)
#'
SampleUMI <- function(
  data,
  max.umi = 1000,
  upsample = FALSE,
  progress.bar = FALSE
) {
  data <- as(data, "dgCMatrix")
  if (length(x = max.umi) == 1) {
    return(
      RunUMISampling(
        data = data,
        sample_val = max.umi,
        upsample = upsample,
        display_progress = progress.bar
      )
    )
  } else if (length(x = max.umi) != ncol(x = data)) {
    stop("max.umi vector not equal to number of cells")
  }
  new_data = RunUMISamplingPerCell(
    data = data,
    sample_val = max.umi,
    upsample = upsample,
    display_progress = progress.bar
  )
  dimnames(new_data) <- dimnames(data)
  return(new_data)
}

#' @export
#'
FindVariableFeatures.default <- function(
  object,
  selection.method = "vst",
  loess.span = 0.3,
  clip.max = 'auto',
  mean.function = FastExpMean,
  dispersion.function = FastLogVMR,
  num.bin = 20,
  binning.method = "equal_width",
  verbose = TRUE
) {
  if (!inherits(x = object, 'Matrix')) {
    object <- as(object = as.matrix(x = object), Class = 'Matrix')
  }
  if (!inherits(x = object, what = 'dgCMatrix')) {
    object <- as(object = object, Class = 'dgCMatrix')
  }
  if (selection.method == "vst") {
    if (clip.max == 'auto') {
      clip.max <- sqrt(ncol(object))
    }
    hvf.info <- data.frame(mean = rowMeans(x = object))
    hvf.info$variance <- SparseRowVar2(mat = object, mu = hvf.info$mean, display_progress = verbose)
    hvf.info$variance.expected <- 0
    hvf.info$variance.standardized <- 0

    not.const <- hvf.info$variance > 0
    fit <- loess(log10(variance) ~ log10(mean), data = hvf.info[not.const, ], span = loess.span)
    hvf.info$variance.expected[not.const] <- 10 ^ fit$fitted
    # use c function to get variance after feature standardization
    hvf.info$variance.standardized <- SparseRowVarStd(mat = object,
                                                      mu = hvf.info$mean,
                                                      sd = sqrt(hvf.info$variance.expected),
                                                      vmax = clip.max,
                                                      display_progress = verbose)
  } else {
    if (!inherits(x = mean.function, what = 'function')) {
      stop("'mean.function' must be a function")
    }
    if (!inherits(x = dispersion.function, what = 'function')) {
      stop("'dispersion.function' must be a function")
    }
    feature.mean <- mean.function(object, verbose)
    feature.dispersion <- dispersion.function(object, verbose)
    names(x = feature.mean) <- names(x = feature.dispersion) <- rownames(x = object)
    feature.dispersion[is.na(x = feature.dispersion)] <- 0
    feature.mean[is.na(x = feature.mean)] <- 0
    data.x.breaks <- switch(
      EXPR = binning.method,
      'equal_width' = num.bin,
      'equal_frequency' = c(
        -1,
        quantile(
          x = feature.mean[feature.mean > 0],
          probs = seq.int(from = 0, to = 1, length.out = num.bin)
        )
      ),
      stop("Unknown binning method: ", binning.method)
    )
    data.x.bin <- cut(x = feature.mean, breaks = data.x.breaks)
    names(x = data.x.bin) <- names(x = feature.mean)
    mean.y <- tapply(X = feature.dispersion, INDEX = data.x.bin, FUN = mean)
    sd.y <- tapply(X = feature.dispersion, INDEX = data.x.bin, FUN = sd)
    feature.dispersion.scaled <- (feature.dispersion - mean.y[as.numeric(x = data.x.bin)]) /
      sd.y[as.numeric(x = data.x.bin)]
    names(x = feature.dispersion.scaled) <- names(x = feature.mean)
    hvf.info <- data.frame(feature.mean, feature.dispersion, feature.dispersion.scaled)
    rownames(x = hvf.info) <- rownames(x = object)
    colnames(x = hvf.info) <- c('mean', 'dispersion', 'dispersion.scaled')
  }
  return(hvf.info)
}

#' @param num.features Number of features to select as top variable features;
#' only used when \code{selection.method = 'dispersion'}
#' @param mean.cutoff A two-length numeric vector with low- and high-cutoffs for
#' feature means
#' @param dispersion.cutoff A two-length numeric vector with low- and high-cutoffs for
#' feature dispersions
#'
#' @describeIn FindVariableFeatures Find variable features in an Assay object
#' @export
#' @method FindVariableFeatures Assay
#'
FindVariableFeatures.Assay <- function(
  object,
  selection.method = "vst",
  loess.span = 0.3,
  clip.max = 'auto',
  mean.function = FastExpMean,
  dispersion.function = FastLogVMR,
  num.bin = 20,
  binning.method = "equal_width",
  num.features = 1000,
  mean.cutoff = c(0.1, 8),
  dispersion.cutoff = c(1, Inf),
  verbose = TRUE
) {
  if (length(x = mean.cutoff) != 2 || length(x = dispersion.cutoff) != 2) {
    stop("Both 'mean.cutoff' and 'dispersion.cutoff' must be two numbers")
  }
  slot <- "data"
  if (selection.method == "vst") {
    slot <- "counts"
  }
  hvf.info <- FindVariableFeatures(
    object = GetAssayData(object = object, slot = slot),
    loess.span = loess.span,
    clip.max = clip.max,
    mean.function = mean.function,
    dispersion.function = dispersion.function,
    num.bin = num.bin,
    binning.method = binning.method,
    verbose = verbose
  )
  object[[names(x = hvf.info)]] <- hvf.info
  if (selection.method == "vst"){
    hvf.info <- hvf.info[order(hvf.info$variance.standardized, decreasing = TRUE), , drop = FALSE]
  } else {
    hvf.info <- hvf.info[order(hvf.info$dispersion, decreasing = TRUE), , drop = FALSE]
  }
  top.features <- switch(
    EXPR = selection.method,
    'mean.var.plot' = {
      means.use <- (hvf.info[, 1] > mean.cutoff[1]) & (hvf.info[, 1] < mean.cutoff[2])
      dispersions.use <- (hvf.info[, 3] > dispersion.cutoff[1]) & (hvf.info[, 3] < dispersion.cutoff[2])
      rownames(x = hvf.info)[which(x = means.use & dispersions.use)]
    },
    'dispersion' = head(x = rownames(x = hvf.info), n = num.features),
    'vst' = head(x = rownames(x = hvf.info), n = num.features),
    stop("Unkown selection method: ", selection.method)
  )
  VariableFeatures(object = object) <- top.features
  return(object)
}

#' @inheritParams FindVariableFeatures.Assay
#' @param assay.use Assay to use
#' @param workflow.name Name of workflow
#'
#' @describeIn FindVariableFeatures Find variable features in a Seurat object
#' @export
#' @method FindVariableFeatures Seurat
#'
FindVariableFeatures.Seurat <- function(
  object,
  assay.use = NULL,
  selection.method = "vst",
  loess.span = 0.3,
  clip.max = 'auto',
  mean.function = FastExpMean,
  dispersion.function = FastLogVMR,
  num.bin = 20,
  binning.method = "equal_width",
  num.features = 1000,
  mean.cutoff = c(0.1, 8),
  dispersion.cutoff = c(1, Inf),
  verbose = TRUE,
  workflow.name = NULL
) {
  if (!is.null(workflow.name)) {
    object <- PrepareWorkflow(object = object, workflow.name = workflow.name)
  }
  assay.use <- assay.use %||% DefaultAssay(object = object)
  assay.data <- GetAssay(object = object, assay.use = assay.use)
  assay.data <- FindVariableFeatures(
    object = assay.data,
    selection.method = selection.method,
    loess.span = loess.span,
    clip.max = clip.max,
    mean.function = mean.function,
    dispersion.function = dispersion.function,
    num.bin = num.bin,
    binning.method = binning.method,
    num.features = num.features,
    mean.cutoff = mean.cutoff,
    dispersion.cutoff = dispersion.cutoff,
    verbose = verbose
  )
  object[[assay.use]] <- assay.data
  object <- LogSeuratCommand(object = object)
  if (!is.null(workflow.name)) {
    object <- UpdateWorkflow(object = object, workflow.name = workflow.name)
  }
  return(object)
}

#' Return a subset of the Seurat object
#'
#' Creates a Seurat object containing only a subset of the cells in the
#' original object. Takes either a list of cells to use as a subset, or a
#' parameter (for example, a gene), to subset on.
#'
#' @param object Seurat object
#' @param subset.names Parameters to subset on. Eg, the name of a gene, PC1, a
#' column name in object@@meta.data, etc. Any argument that can be retreived
#' using FetchData
#' @param low.thresholds Low cutoffs for the parameters (default is -Inf)
#' @param high.thresholds High cutoffs for the parameters (default is Inf)
#' @param cells.use A vector of cell names to use as a subset
#'
#' @return Returns a Seurat object containing only the relevant subset of cells
#'
#' @export
#'
#' @examples
#' head(x = FetchData(object = pbmc_small, vars.all = 'LTB'))
#' pbmc_filtered <- FilterCells(
#'   object = pbmc_small,
#'   subset.names = 'LTB',
#'   high.thresholds = 6
#' )
#' head(x = FetchData(object = pbmc_filtered, vars.all = 'LTB'))
#'
FilterCells <- function(
  object,
  subset.names,
  low.thresholds,
  high.thresholds,
  cells.use = NULL
) {
  if (missing(x = low.thresholds)) {
    low.thresholds <- replicate(n = length(x = subset.names), expr = -Inf)
  }
  if (missing(x = high.thresholds)) {
    high.thresholds <- replicate(n = length(x = subset.names), expr = Inf)
  }
  length.check <- sapply(
    X = list(subset.names, low.thresholds, high.thresholds),
    FUN = length
  )
  if (length(x = unique(x = length.check)) != 1) {
    stop("'subset.names', 'low.thresholds', and 'high.thresholds' must all have the same length")
  }
  data.subsets <- data.frame(subset.names, low.thresholds, high.thresholds)
  cells.use <- cells.use %||% colnames(x = object)
  for (i in seq(nrow(data.subsets))) {
    cells.use <- tryCatch(
      expr = WhichCells(
        object = object,
        cells.use = cells.use,
        subset.name = data.subsets[i, 1],
        low.threshold = data.subsets[i, 2],
        high.threshold = data.subsets[i, 3]
      ),
      error = function(e) {
        warning(e)
        cells.use
      }
    )
  }
  object <- SubsetData(object, cells.use = cells.use)
  object <- LogSeuratCommand(object)
  return(object)
}
