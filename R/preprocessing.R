#' @include generics.R
#'
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Filter out cells from a Seurat object
#'
#' Creates a Seurat object containing only a subset of the cells in the
#' original object. Takes a vector of variables on which to subset along with
#' corresponding low and high threshold values for each variable, returning an
#' object contain only cells that pass the specified filters. You can also
#' restrict the filtering to only a subset of cells using cells (i.e. only
#' keep cells that BOTH pass the filters and are in the cells list).
#'
#' @param object Seurat object
#' @param subset.names Parameters to subset on. Eg, the name of a gene, PC1, a
#' column name in object@@meta.data, etc. Any argument that can be retreived
#' using FetchData
#' @param low.thresholds Low cutoffs for the parameters (default is -Inf)
#' @param high.thresholds High cutoffs for the parameters (default is Inf)
#' @param cells A vector of cell names to use as a subset
#'
#' @return Returns a Seurat object containing only the relevant subset of cells
#'
#' @export
#'
#' @examples
#' head(x = FetchData(object = pbmc_small, vars = 'LTB'))
#' pbmc_filtered <- FilterCells(
#'   object = pbmc_small,
#'   subset.names = 'LTB',
#'   high.thresholds = 6
#' )
#' head(x = FetchData(object = pbmc_filtered, vars = 'LTB'))
#'
FilterCells <- function(
  object,
  subset.names,
  low.thresholds,
  high.thresholds,
  cells = NULL
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
  cells <- cells %||% colnames(x = object)
  for (i in seq(nrow(data.subsets))) {
    cells <- tryCatch(
      expr = WhichCells(
        object = object,
        cells = cells,
        subset.name = data.subsets[i, 1],
        low.threshold = data.subsets[i, 2],
        high.threshold = data.subsets[i, 3]
      ),
      error = function(e) {
        warning(e)
        cells
      }
    )
  }
  object <- SubsetData(object = object, cells = cells)
  object <- LogSeuratCommand(object = object)
  return(object)
}

#' Demultiplex samples based on data from cell 'hashing'
#'
#' Assign sample-of-origin for each cell, annotate doublets.
#'
#' @param object Seurat object. Assumes that the hash tag oligo (HTO) data has been added and normalized.
#' @param assay Name of the Hashtag assay (HTO by default)
#' @param positive_quantile The quantile of inferred 'negative' distribution for each hashtag - over which the cell is considered 'positive'. Default is 0.99
#' @param init_centers Initial number of clusters for hashtags. Default is the # of hashtag oligo names + 1 (to account for negatives)
#' @param k_function Clustering function for initial hashtag grouping. Default is "clara" for fast k-medoids clustering on large applications, also support "kmeans" for kmeans clustering
#' @param nsamples Number of samples to be drawn from the dataset used for clustering, for k_function = "clara"
#' @param cluster_nstarts nstarts value for k-means clustering (for k_function = "kmeans"). 100 by default
#' @param verbose Prints the output
#'
#' @return The Seurat object with the following demultiplexed information stored in the meta data:
#' \item{hash_maxID}{Name of hashtag with the highest signal}
#' \item{hash_secondID}{Name of hashtag with the second highest signal}
#' \item{hash_margin}{The difference between signals for hash_maxID and hash_secondID}
#' \item{hto_classification}{Classification result, with doublets/multiplets named by the top two highest hashtags}
#' \item{hto_classification_global}{Global classification result (singlet, doublet or negative)}
#' \item{hash_ID}{Classification result where doublet IDs are collapsed}
#'
#' @importFrom cluster clara
#' @importFrom Matrix colSums
#' @importFrom fitdistrplus fitdist
#' @importFrom stats pnbinom kmeans
#' @export
#'
#' @examples
#' \dontrun{
#' object <- HTODemux(object)
#' }
#'
HTODemux <- function(
  object,
  assay = "HTO",
  positive_quantile = 0.99,
  init_centers = NULL,
  cluster_nstarts = 100,
  k_function = "clara",
  nsamples = 100,
  verbose = TRUE
) {
  #initial clustering
  assay <- assay %||% DefaultAssay(object = object)
  hash_data <- GetAssayData(object = object, assay = assay)
  hash_raw_data <- GetAssayData(
    object = object,
    assay = assay,
    slot = 'counts'
  )[, colnames(x = object)]
  hash_raw_data <- as.matrix(x = hash_raw_data)
  ncenters <- init_centers %||% nrow(x = hash_data) + 1
  if (k_function == "kmeans") {
    hto_init_clusters <- kmeans(
      x = t(x = GetAssayData(object = object, assay = assay)),
      centers = ncenters,
      nstart = cluster_nstarts
    )
    #identify positive and negative signals for all HTO
    Idents(object = object, cells = names(x = hto_init_clusters$cluster)) <- hto_init_clusters$cluster
  } else {
    #use fast k-medoid clustering
    hto_init_clusters <- clara(
      x = t(x = GetAssayData(object = object, assay = assay)),
      k = ncenters,
      samples = nsamples
    )
    #identify positive and negative signals for all HTO
    Idents(object = object, cells = names(x = hto_init_clusters$clustering)) <- hto_init_clusters$clustering
  }
  #average hto signals per cluster
  #work around so we don't average all the RNA levels which takes time
  average_hto <- AverageExpression(
    object = object,
    assay = assay,
    verbose = FALSE
  )[[assay]]
  #create a matrix to store classification result
  hto_discrete <- GetAssayData(object = object, assay = assay)
  hto_discrete[hto_discrete > 0] <- 0
  # for each HTO, we will use the minimum cluster for fitting
  for (hto_iter in rownames(x = hash_data)) {
    hto_values <- hash_raw_data[hto_iter, colnames(object)]
    #commented out if we take all but the top cluster as background
    #hto_values_negative=hto_values[setdiff(object@cell.names,WhichCells(object,which.max(average_hto[hto_iter,])))]
    hto_values_use <- hto_values[WhichCells(object = object, ident.keep = levels(Idents(object))[[which.min(x = average_hto[hto_iter, ])]])]
    hto_fit <- suppressWarnings(fitdist(hto_values_use, "nbinom"))
    hto_cutoff <- as.numeric(x = quantile(x = hto_fit, probs = positive_quantile)$quantiles[1])
    hto_discrete[hto_iter, names(x = which(x = hto_values > hto_cutoff))] <- 1
    if (verbose) {
      message(paste0("Cutoff for ", hto_iter, " : ", hto_cutoff, " reads"))
    }
  }
  # now assign cells to HTO based on discretized values
  num_hto_positive <- colSums(x = hto_discrete)
  hto_classification_global <- num_hto_positive
  hto_classification_global[num_hto_positive == 0] <- "Negative"
  hto_classification_global[num_hto_positive == 1] <- "Singlet"
  hto_classification_global[num_hto_positive > 1] <- "Doublet"
  donor_id = rownames(x = hash_data)
  hash_max <- apply(X = hash_data, MARGIN = 2, FUN = max)
  hash_maxID <- apply(X = hash_data, MARGIN = 2, FUN = which.max)
  hash_second <- apply(X = hash_data, MARGIN = 2, FUN = MaxN, N = 2)
  hash_maxID <- as.character(x = donor_id[sapply(
    X = 1:ncol(x = hash_data),
    FUN = function(x) {
      return(which(x = hash_data[, x] == hash_max[x])[1])
    }
  )])
  hash_secondID <- as.character(x = donor_id[sapply(
    X = 1:ncol(x = hash_data),
    FUN = function(x) {
      return(which(x = hash_data[, x] == hash_second[x])[1])
    }
  )])
  hash_margin <- hash_max - hash_second
  doublet_id <- sapply(
    X = 1:length(x = hash_maxID),
    FUN = function(x) {
      return(paste(sort(c(hash_maxID[x], hash_secondID[x])), collapse = "_"))
    }
  )
  # doublet_names <- names(x = table(doublet_id))[-1] # Not used
  hto_classification <- hto_classification_global
  hto_classification[hto_classification_global == "Negative"] <- "Negative"
  hto_classification[hto_classification_global == "Singlet"] <- hash_maxID[which(x = hto_classification_global == "Singlet")]
  hto_classification[hto_classification_global == "Doublet"] <- doublet_id[which(x = hto_classification_global == "Doublet")]
  classification_metadata <- data.frame(
    hash_maxID,
    hash_secondID, hash_margin,
    hto_classification,
    hto_classification_global
  )
  object <- AddMetaData(object = object, metadata = classification_metadata)
  Idents(object) <- "hto_classification"
  Idents(object, cells = rownames(object@meta.data[object@meta.data$hto_classification_global == "Doublet", ])) <- "Doublet"
  object@meta.data$hash_ID <- Idents(object)
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
# @importFrom hdf5r H5File
#'
#' @export
#'
Read10X_h5 <- function(filename, ensg.names = FALSE) {
  if (!require(hdf5r)) {
    stop("Please install hdf5r to read HDF5 files")
  }
  if (!file.exists(filename)) {
    stop("File not found")
  }
  infile <- H5File$new(filename)
  genomes <- names(infile)
  output <- list()
  for (genome in genomes) {
    counts <- infile[[paste0(genome, '/data')]]
    indices <- infile[[paste0(genome, '/indices')]]
    indptr <- infile[[paste0(genome, '/indptr')]]
    shp <- infile[[paste0(genome, '/shape')]]
    if (ensg.names) {
      gene_names <- infile[[paste0(genome, '/genes')]][]
    } else {
      gene_names <- make.unique(infile[[paste0(genome, '/gene_names')]][])
    }
    barcodes <- infile[[paste0(genome, '/barcodes')]]
    sparse.mat <- sparseMatrix(
      i = indices[] + 1, p = indptr[],
      x = as.numeric(counts[]),
      dims = shp[], giveCsparse = FALSE
    )
    rownames(sparse.mat) <- gene_names
    colnames(sparse.mat) <- barcodes[]
    output[[genome]] <- sparse.mat
  }
  infile$close_all()
  if (length(output) == 1) {
    return(output[[genome]])
  } else{
    return(output)
  }
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
#' @param assay Name of assay to use
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
  assay = NULL,
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
  assay <- assay %||% DefaultAssay(object = object)
  assay.obj <- GetAssay(object = object, assay = assay)
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
  object[[assay]] <- assay.obj

  # save vst output (except y) in @misc slot
  vst.out$y <- NULL
  object@misc[['vst.out']] <- vst.out
  return(object)
}

#' Sample UMI
#'
#' Downsample each cell to a specified number of UMIs. Includes
#' an option to upsample cells below specified UMI as well.
#'
#' @param data Matrix with the raw count data
#' @param max.umi Number of UMIs to sample to
#' @param upsample Upsamples all cells with fewer than max.umi
#' @param verbose Display the progress bar
#'
#' @import Matrix
#' @importFrom methods as
#'
#' @return Matrix with downsampled data
#'
#' @export
#'
#' @examples
#' counts = as.matrix(x = GetAssayData(object = pbmc_small, assay = "RNA", slot = "counts"))
#' downsampled = SampleUMI(data = counts)
#' head(x = downsampled)
#'
SampleUMI <- function(
  data,
  max.umi = 1000,
  upsample = FALSE,
  verbose = FALSE
) {
  data <- as(object = data, Class = "dgCMatrix")
  if (length(x = max.umi) == 1) {
    return(
      RunUMISampling(
        data = data,
        sample_val = max.umi,
        upsample = upsample,
        display_progress = verbose
      )
    )
  } else if (length(x = max.umi) != ncol(x = data)) {
    stop("max.umi vector not equal to number of cells")
  }
  new_data = RunUMISamplingPerCell(
    data = data,
    sample_val = max.umi,
    upsample = upsample,
    display_progress = verbose
  )
  dimnames(new_data) <- dimnames(data)
  return(new_data)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Methods for Seurat-defined generics
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
      clip.max <- sqrt(x = ncol(x = object))
    }
    hvf.info <- data.frame(mean = rowMeans(x = object))
    hvf.info$variance <- SparseRowVar2(
      mat = object,
      mu = hvf.info$mean,
      display_progress = verbose
    )
    hvf.info$variance.expected <- 0
    hvf.info$variance.standardized <- 0
    not.const <- hvf.info$variance > 0
    fit <- loess(
      formula = log10(x = variance) ~ log10(x = mean),
      data = hvf.info[not.const, ],
      span = loess.span
    )
    hvf.info$variance.expected[not.const] <- 10 ^ fit$fitted
    # use c function to get variance after feature standardization
    hvf.info$variance.standardized <- SparseRowVarStd(
      mat = object,
      mu = hvf.info$mean,
      sd = sqrt(hvf.info$variance.expected),
      vmax = clip.max,
      display_progress = verbose
    )
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

#' @param nfeatures Number of features to select as top variable features;
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
  nfeatures = 1000,
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
    selection.method = selection.method,
    loess.span = loess.span,
    clip.max = clip.max,
    mean.function = mean.function,
    dispersion.function = dispersion.function,
    num.bin = num.bin,
    binning.method = binning.method,
    verbose = verbose
  )
  object[[names(x = hvf.info)]] <- hvf.info
  if (selection.method == "vst") {
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
    'dispersion' = head(x = rownames(x = hvf.info), n = nfeatures),
    'vst' = head(x = rownames(x = hvf.info), n = nfeatures),
    stop("Unkown selection method: ", selection.method)
  )
  VariableFeatures(object = object) <- top.features
  return(object)
}

#' @inheritParams FindVariableFeatures.Assay
#' @param assay Assay to use
#' @param workflow.name Name of workflow
#'
#' @describeIn FindVariableFeatures Find variable features in a Seurat object
#' @export
#' @method FindVariableFeatures Seurat
#'
FindVariableFeatures.Seurat <- function(
  object,
  assay = NULL,
  selection.method = "vst",
  loess.span = 0.3,
  clip.max = 'auto',
  mean.function = FastExpMean,
  dispersion.function = FastLogVMR,
  num.bin = 20,
  binning.method = "equal_width",
  nfeatures = 1000,
  mean.cutoff = c(0.1, 8),
  dispersion.cutoff = c(1, Inf),
  verbose = TRUE,
  workflow.name = NULL
) {
  if (!is.null(workflow.name)) {
    object <- PrepareWorkflow(object = object, workflow.name = workflow.name)
  }
  assay <- assay %||% DefaultAssay(object = object)
  assay.data <- GetAssay(object = object, assay = assay)
  assay.data <- FindVariableFeatures(
    object = assay.data,
    selection.method = selection.method,
    loess.span = loess.span,
    clip.max = clip.max,
    mean.function = mean.function,
    dispersion.function = dispersion.function,
    num.bin = num.bin,
    binning.method = binning.method,
    nfeatures = nfeatures,
    mean.cutoff = mean.cutoff,
    dispersion.cutoff = dispersion.cutoff,
    verbose = verbose
  )
  object[[assay]] <- assay.data
  object <- LogSeuratCommand(object = object)
  if (!is.null(workflow.name)) {
    object <- UpdateWorkflow(object = object, workflow.name = workflow.name)
  }
  return(object)
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

#' @param assay Name of assay to use
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
  assay = NULL,
  normalization.method = "LogNormalize",
  scale.factor = 1e4,
  verbose = TRUE,
  workflow.name = NULL
) {
  assay <- assay %||% DefaultAssay(object = object)
  if (!is.null(workflow.name)) {
    object <- PrepareWorkflow(object = object, workflow.name = workflow.name)
  }
  assay.data <- GetAssay(object = object, assay = assay)
  assay.data <- NormalizeData(
    object = assay.data,
    normalization.method = normalization.method,
    scale.factor = scale.factor,
    verbose = verbose
  )
  object[[assay]] <- assay.data
  object <- LogSeuratCommand(object = object)
  if (!is.null(workflow.name)) {
    object <- UpdateWorkflow(object = object, workflow.name = workflow.name)
  }
  return(object)
}


#' @export
#'
RunALRA.default <- function(object, k = NULL, q = 10) {
  A.norm <- t(x = as.matrix(x = object))
  message("Identifying non-zero values")
  originally.nonzero <- A.norm > 0
  message("Computing Randomized SVD")
  fastDecomp.noc <- rsvd(A = A.norm, k = k, q = q)
  A.norm.rank.k <- fastDecomp.noc$u[, 1:k] %*%
    diag(x = fastDecomp.noc$d[1:k]) %*%
    t(x = fastDecomp.noc$v[,1:k])
  message("Find most negative values of each gene")
  A.norm.rank.k.mins <- abs(x = apply(X = A.norm.rank.k, MARGIN = 2, FUN = min))
  message("Thresholding by the most negative value of each gene")
  A.norm.rank.k.cor <- replace(
    x = A.norm.rank.k,
    list = A.norm.rank.k <= A.norm.rank.k.mins[col(A.norm.rank.k)],
    values = 0
  )
  sd.nonzero <- function(x) {
    return(sd(x[!x == 0]))
  }
  sigma.1 <- apply(X = A.norm.rank.k.cor, MARGIN = 2, FUN = sd.nonzero)
  sigma.2 <- apply(X = A.norm, MARGIN = 2, FUN = sd.nonzero)
  mu.1 <- colSums(x = A.norm.rank.k.cor) / colSums(x = !!A.norm.rank.k.cor)
  mu.2 <- colSums(x = A.norm) / colSums(x = !!A.norm)
  toscale <- !is.na(x = sigma.1) & !is.na(x = sigma.2)
  message(sprintf(fmt = "Scaling all except for %d columns", sum(!toscale)))
  sigma.1.2 <- sigma.2 / sigma.1
  toadd <- -1 * mu.1 * sigma.2 / sigma.1 + mu.2
  A.norm.rank.k.temp <- A.norm.rank.k.cor[, toscale]
  A.norm.rank.k.temp <- sweep(
    x = A.norm.rank.k.temp,
    MARGIN = 2,
    STATS = sigma.1.2[toscale],
    FUN = "*"
  )
  A.norm.rank.k.temp <- sweep(
    x = A.norm.rank.k.temp,
    MARGIN = 2,
    STATS = toadd[toscale],
    FUN = "+"
  )
  A.norm.rank.k.cor.sc <- A.norm.rank.k.cor
  A.norm.rank.k.cor.sc[, toscale] <- A.norm.rank.k.temp
  A.norm.rank.k.cor.sc[A.norm.rank.k.cor == 0] <- 0
  lt0 <- A.norm.rank.k.cor.sc < 0
  A.norm.rank.k.cor.sc[lt0] <- 0
  message(sprintf(
    fmt = "%.2f%% of the values became negative in the scaling process and were set to zero",
    100 * sum(lt0) / prod(dim(x = A.norm))
  ))
  A.norm.rank.k.cor.sc[originally.nonzero & A.norm.rank.k.cor.sc == 0] <-
    A.norm[originally.nonzero & A.norm.rank.k.cor.sc == 0]
  colnames(x = A.norm.rank.k) <- colnames(x = A.norm.rank.k.cor.sc) <-
    colnames(x = A.norm.rank.k.cor) <- colnames(x = A.norm)
  original.nz <- sum(A.norm > 0) / prod(dim(x = A.norm))
  completed.nz <- sum(A.norm.rank.k.cor.sc > 0) / prod(dim(x = A.norm))
  message(sprintf(
    fmt = "The matrix went from %.2f%% nonzero to %.2f%% nonzero",
    100 * original.nz,
    100 * completed.nz
  ))
  return(A.norm.rank.k.cor.sc)
}

#' @param assay Assay to use
#' @param slot slot to use
#' @param setDefaultAssay If TRUE, will set imputed results as default Assay
#' @param genes.use genes to impute
#' @param K Number of singular values to compute when choosing k. Must be less
#' than the smallest dimension of the matrix.
#' @param p.val.th  The threshold for ''significance'' when choosing k
#' @param noise.start  Index for which all smaller singular values are considered noise
#' @param q.k  Number of additional power iterations when choosing k
#' @param k.only If TRUE, only computes optimal k WITHOUT performing ALRA
#'
#' @importFrom rsvd rsvd
#' @importFrom Matrix Matrix
#' @importFrom stats pnorm sd
#'
#' @describeIn RunALRA Run ALRA on a Seurat object
#' @export
#' @method RunALRA Seurat
#'
RunALRA.Seurat <- function(
  object,
  k = NULL,
  q = 10,
  assay = NULL,
  slot = "data",
  setDefaultAssay = TRUE,
  genes.use = NULL,
  K = NULL,
  p.val.th = 1e-10,
  noise.start = NULL,
  q.k = 2,
  k.only = FALSE
) {
  if (!is.null(x = k) && k.only) {
    warning("Stop: k is already given, set k.only = FALSE or k = NULL")
  }
  genes.use <- genes.use %||% rownames(x = object)
  assay <- assay %||% DefaultAssay(object = object)
  # Add slot alra in object@tools
  if (is.null(x = object@tools[["alra"]])) {
    object@tools[["alra"]] <- list()
  }
  # Check if k is already stored
  if (is.null(x = k) & !is.null(x = object@tools[["alra"]][["k"]])) {
    k <- object@tools[["alra"]][["k"]]
    message("Using previously computed value of k")
  }
  data.used <- GetAssayData(object = object, assay = assay, slot = slot)[genes.use,]
  # Choose k with heuristics if k is not given
  if (is.null(x = k)) {
    # set K based on data dimension
    if (is.null(x = K)) {
      K <- 100
      if (K > min(dim(x = data.used))) {
        K <- min(dim(x = data.used))
        warning("For best performance, we recommend using ALRA on expression matrices larger than 100 by 100")
      }
    }
    if (K > min(dim(x = data.used))) {
      stop("For an m by n data, K must be smaller than the min(m,n)")
    }
    # set noise.start based on K
    if (is.null(x = noise.start)) {
      noise.start <- K - 20
      if (noise.start <= 0) {
        noise.start <- max(K - 5, 1)
      }
    }
    if (noise.start > K - 5) {
      stop("There need to be at least 5 singular values considered noise")
    }
    noise.svals <- noise.start:K
    rsvd.out <- rsvd(A = t(x = as.matrix(x = data.used)), k = K, q = q.k)
    diffs <- diff(x = rsvd.out$d)
    pvals <- pnorm(
      q = diffs,
      mean = mean(x = diffs[noise.svals - 1]),
      sd = sd(x = diffs[noise.svals - 1])
    )
    k <- max(which(x = pvals < p.val.th))
    object@tools[["alra"]][["d"]] <- rsvd.out$d
    object@tools[["alra"]][["k"]] <- k
    object@tools[["alra"]][["diffs"]] <- diffs
    object@tools[["alra"]][["pvals"]] <- pvals
  }
  if (k.only) {
    message("Chose rank k = ", k, ", WITHOUT performing ALRA")
    return(object)
  }
  message("Rank k = ", k)
  # Perform ALRA on data.used
  output.alra <- RunALRA(object = data.used, k = k, q = q)
  # Save ALRA data in object@assay
  data.alra <- Matrix(data = t(x = output.alra), sparse = TRUE)
  rownames(x = data.alra) <- genes.use
  colnames(x = data.alra) <- colnames(x = object)
  assay.alra <- CreateAssayObject(data = data.alra)
  object[["alra"]] <- assay.alra
  if (setDefaultAssay) {
    message("Setting default assay as alra")
    DefaultAssay(object = object) <- "alra"
  }
  return(object)
}

#' @export
#'
ScaleData.default <- function(
  object,
  features = NULL,
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
  features <- features %||% rownames(x = object)
  features <- as.vector(x = intersect(x = features, y = rownames(x = object)))
  object <- object[features, , drop = FALSE]
  scaled.data <- matrix(data = NA, nrow = nrow(x = object), ncol = ncol(x = object))
  dimnames(x = scaled.data) <- dimnames(x = object)
  min.cells.to.block <- min(min.cells.to.block, ncol(x = object))
  Parenting(
    parent.find = "ScaleData.Assay",
    features = features,
    min.cells.to.block = min.cells.to.block
  )
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
      features.regress = features,
      model.use = model.use,
      use.umi = use.umi,
      verbose = verbose
    )
    gc(verbose = FALSE)
  }
  max.block <- ceiling(x = length(x = features) / block.size)
  if (verbose) {
    message("Scaling data matrix")
    pb <- txtProgressBar(min = 0, max = max.block, style = 3, file = stderr())
  }
  for (i in 1:max.block) {
    my.inds <- ((block.size * (i - 1)):(block.size * i - 1)) + 1
    my.inds <- my.inds[my.inds <= length(x = features)]
    if (inherits(x = object, what = c('dgCMatrix', 'dgTMatrix'))) {
      scale.function <- FastSparseRowScale
    } else {
      object <- as.matrix(x = object)
      scale.function <- FastRowScale
    }
    data.scale <- scale.function(
      mat = object[features[my.inds], , drop = FALSE],
      scale = do.scale,
      center = do.center,
      scale_max = scale.max,
      display_progress = FALSE
    )
    dimnames(x = data.scale) <- dimnames(x = object[features[my.inds], ])
    scaled.data[features[my.inds], ] <- data.scale
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
  features = NULL,
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
  features <- features %||% VariableFeatures(object)
  if (length(features) == 0) features <- rownames(GetAssayData(object = object, slot = slot.use))
  object <- SetAssayData(
    object = object,
    slot = 'scale.data',
    new.data = ScaleData(
      object = GetAssayData(object = object, slot = slot.use),
      features = features,
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
  Parenting(
    parent.find = "ScaleData.Seurat",
    features = features,
    min.cells.to.block = min.cells.to.block,
    use.umi = use.umi
  )
  return(object)
}

#' @param assay Name of Assay to scale
#' @param workflow.name Name of workflow
#'
#' @describeIn ScaleData Scale a Seurat object
#' @export
#' @method ScaleData Seurat
#'
ScaleData.Seurat <- function(
  object,
  features = NULL,
  assay = NULL,
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
  assay <- assay %||% DefaultAssay(object = object)
  assay.data <- GetAssay(object = object, assay = assay)
  if (any(vars.to.regress %in% colnames(x = object[]))) {
    latent.data <- object[vars.to.regress[vars.to.regress %in% colnames(x = object[])]]
  } else {
    latent.data <- NULL
  }
  assay.data <- ScaleData(
    object = assay.data,
    features = features,
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
  object[[assay]] <- assay.data
  object <- LogSeuratCommand(object = object)
  if (!is.null(workflow.name)) {
    object <- UpdateWorkflow(object = object, workflow.name = workflow.name)
  }
  return(object)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Internal
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Normalize raw data
#
# Normalize count data per cell and transform to centered log ratio
#
# @param data Matrix with the raw count data
# @param custom_function A custom normalization function
# @parm across Which way to we normalize? Choose form 'cells' or 'features'
#
# @return Returns a matrix with the custom normalization
#
#' @importFrom methods as
# @import Matrix
#
CustomNormalize <- function(data, custom_function, across) {
  if (class(x = data) == "data.frame") {
    data <- as.matrix(x = data)
  }
  if (class(x = data) != "dgCMatrix") {
    data <- as(object = data, Class = "dgCMatrix")
  }
  margin <- switch(
    EXPR = across,
    'cells' = 2,
    'features' = 1,
    stop("'across' must be either 'cells' or 'features'")
  )
  norm.data <- apply(
    X = data,
    MARGIN = margin,
    FUN = custom_function)
  if (margin == 1) {
    norm.data = t(x = norm.data)
  }
  colnames(x = norm.data) <- colnames(x = data)
  rownames(x = norm.data) <- rownames(x = data)
  return(norm.data)
}

#
#' @importFrom stats residuals
#
NBResiduals <- function(fmla, regression.mat, gene, return.mode = FALSE) {
  fit <- 0
  try(
    fit <- glm.nb(
      formula = fmla,
      data = regression.mat
    ),
    silent = TRUE)
  if (class(fit)[1] == 'numeric') {
    message(sprintf('glm.nb failed for gene %s; falling back to scale(log(y+1))', gene))
    resid <- scale(x = log(x = regression.mat[, 'GENE'] + 1))[, 1]
    mode <- 'scale'
  } else {
    resid <- residuals(fit, type = 'pearson')
    mode = 'nbreg'
  }
  do.return <- list(resid = resid, mode = mode)
  if (return.mode) {
    return(do.return)
  } else {
    return(do.return$resid)
  }
}

# Regress out techincal effects and cell cycle from a matrix
#
# Remove unwanted effects from scale.data
#
# @parm data.expr An expression matrix to regress the effects of latent.data out of
# should be the complete expression matrix in genes x cells
# @param latent.data A matrix or data.frame of latent variables, should be cells x latent variables,
# the colnames should be the variables to regress
# @param features.regress An integer vector representing the indices of the genes to run regression on
# @param model.use Model to use, one of 'linear', 'poisson', or 'negbinom'; pass NULL to simply return data.expr
# @param use.umi Regress on UMI count data
# @param verbose Display a progress bar
#
#' @importFrom stats as.formula
#' @importFrom utils txtProgressBar setTxtProgressBar
#
RegressOutMatrix <- function(
  data.expr,
  latent.data = NULL,
  features.regress = NULL,
  model.use = NULL,
  use.umi = FALSE,
  verbose = TRUE,
  ...
) {
  # Do we bypass regression and simply return data.expr?
  bypass <- vapply(
    X = list(latent.data, model.use),
    FUN = is.null,
    FUN.VALUE = logical(length = 1L)
  )
  if (any(bypass)) {
    return(data.expr)
  }
  # Check model.use
  possible.models <- c("linear", "poisson", "negbinom")
  if (!model.use %in% possible.models) {
    stop(paste(
      model.use,
      "is not a valid model. Please use one the following:",
      paste0(possible.models, collapse = ", ")
    ))
  }
  # Check features.regress
  if (is.null(x = features.regress)) {
    features.regress <- 1:nrow(x = data.expr)
  }
  if (is.character(x = features.regress)) {
    features.regress <- intersect(x = features.regress, y = rownames(x = data.expr))
    if (length(x = features.regress) == 0) {
      stop("Cannot use features that are beyond the scope of data.expr")
    }
  } else if (max(features.regress) > nrow(x = data.expr)) {
    stop("Cannot use features that are beyond the scope of data.expr")
  }
  # Check dataset dimensions
  if (nrow(x = latent.data) != ncol(x = data.expr)) {
    stop("Uneven number of cells between latent data and expression data")
  }
  use.umi <- ifelse(test = model.use != 'linear', yes = TRUE, no = use.umi)
  # Create formula for regression
  vars.to.regress <- colnames(x = latent.data)
  fmla <- paste('GENE ~', paste(vars.to.regress, collapse = '+'))
  fmla <- as.formula(object = fmla)
  if (model.use == "linear") {
    # In this code, we'll repeatedly regress different Y against the same X
    # (latent.data) in order to calculate residuals.  Rather that repeatedly
    # call lm to do this, we'll avoid recalculating the QR decomposition for the
    # latent.data matrix each time by reusing it after calculating it once
    regression.mat <- cbind(latent.data, data.expr[1,])
    colnames(regression.mat) <- c(colnames(x = latent.data), "GENE")
    qr <- lm(fmla, data = regression.mat, qr = TRUE)$qr
    rm(regression.mat)
  }
  # Make results matrix
  data.resid <- matrix(
    nrow = nrow(x = data.expr),
    ncol = ncol(x = data.expr)
  )
  if (verbose) {
    pb <- txtProgressBar(char = '=', style = 3)
  }
  for (i in 1:length(x = features.regress)) {
    x <- features.regress[i]
    regression.mat <- cbind(latent.data, data.expr[x, ])
    colnames(x = regression.mat) <- c(vars.to.regress, 'GENE')
    regression.mat <- switch(
      EXPR = model.use,
      'linear' = qr.resid(qr = qr, y = data.expr[x,]),
      'poisson' = residuals(object = glm(
        formula = fmla,
        family = 'poisson',
        data = regression.mat),
        type = 'pearson'
      ),
      'negbinom' = NBResiduals(
        fmla = fmla,
        regression.mat = regression.mat,
        gene = x
      )
    )
    data.resid[i, ] <- regression.mat
    if (verbose) {
      setTxtProgressBar(pb = pb, value = i / length(x = features.regress))
    }
  }
  if (verbose) {
    close(con = pb)
  }
  if (use.umi) {
    data.resid <- log1p(x = sweep(
      x = data.resid,
      MARGIN = 1,
      STATS = apply(X = data.resid, MARGIN = 1, FUN = min),
      FUN = '-'
    ))
  }
  dimnames(x = data.resid) <- dimnames(x = data.expr)
  return(data.resid)
}

# Regress out technical effects and cell cycle using regularized Negative Binomial regression
#
# Remove unwanted effects from umi data and set scale.data to Pearson residuals
# Uses mclapply; you can set the number of cores it will use to n with command options(mc.cores = n)
#
# @param object Seurat object
# @param latent.vars effects to regress out
# @param genes.regress gene to run regression for (default is all genes)
# @param pr.clip.range numeric of length two specifying the min and max values the results will be clipped to
#
# @return Returns Seurat object with the scale.data (object@scale.data) genes returning the residuals fromthe regression model
#
#' @import Matrix
#' @import parallel
#' @importFrom stats glm residuals
#' @importFrom MASS theta.ml negative.binomial
#' @importFrom utils txtProgressBar setTxtProgressBar
#
RegressOutNB <- function(
  object,
  latent.vars,
  genes.regress = NULL,
  pr.clip.range = c(-30, 30),
  min.theta = 0.01
) {
  genes.regress <- SetIfNull(x = genes.regress, default = rownames(x = object@data))
  genes.regress <- intersect(x = genes.regress, y = rownames(x = object@data))
  cm <- object@raw.data[genes.regress, colnames(x = object@data), drop = FALSE]
  latent.data <- FetchData(object = object, vars = latent.vars)
  message(sprintf('Regressing out %s for %d genes\n', paste(latent.vars), length(x = genes.regress)))
  theta.fit <- RegularizedTheta(cm = cm, latent.data = latent.data, min.theta = 0.01, bin.size = 128)
  message('Second run NB regression with fixed theta')
  bin.size <- 128
  bin.ind <- ceiling(1:length(genes.regress)/bin.size)
  max.bin <- max(bin.ind)
  pb <- txtProgressBar(min = 0, max = max.bin, style = 3, file = stderr())
  pr <- c()
  for (i in 1:max.bin) {
    genes.bin.regress <- genes.regress[bin.ind == i]
    bin.pr.lst <- parallel::mclapply(
      X = genes.bin.regress,
      FUN = function(j) {
        fit <- 0
        try(
          expr = fit <- glm(
            cm[j, ] ~ .,
            data = latent.data,
            family = MASS::negative.binomial(theta = theta.fit[j])
          ),
          silent=TRUE
        )
        if (class(fit)[1] == 'numeric') {
          message(
            sprintf(
              'glm and family=negative.binomial(theta=%f) failed for gene %s; falling back to scale(log10(y+1))',
              theta.fit[j],
              j
            )
          )
          res <- scale(log10(cm[j, ] + 1))[, 1]
        } else {
          res <- residuals(fit, type = 'pearson')
        }
        return(res)
      }
    )
    pr <- rbind(pr, do.call(rbind, bin.pr.lst))
    setTxtProgressBar(pb, i)
  }
  close(pb)
  dimnames(x = pr) <- dimnames(x = cm)
  pr[pr < pr.clip.range[1]] <- pr.clip.range[1]
  pr[pr > pr.clip.range[2]] <- pr.clip.range[2]
  object@scale.data <- pr
  return(object)
}
