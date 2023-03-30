#' @include generics.R
#' @importFrom progressr progressor
#'
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

globalVariables(
  names = c('fov', 'cell_ID', 'qv'),
  package = 'Seurat',
  add = TRUE
)
#' Calculate the Barcode Distribution Inflection
#'
#' This function calculates an adaptive inflection point ("knee") of the barcode distribution
#' for each sample group. This is useful for determining a threshold for removing
#' low-quality samples.
#'
#' The function operates by calculating the slope of the barcode number vs. rank
#' distribution, and then finding the point at which the distribution changes most
#' steeply (the "knee"). Of note, this calculation often must be restricted as to the
#' range at which it performs, so `threshold` parameters are provided to restrict the
#' range of the calculation based on the rank of the barcodes. [BarcodeInflectionsPlot()]
#' is provided as a convenience function to visualize and test different thresholds and
#' thus provide more sensical end results.
#'
#' See [BarcodeInflectionsPlot()] to visualize the calculated inflection points and
#' [SubsetByBarcodeInflections()] to subsequently subset the Seurat object.
#'
#' @param object Seurat object
#' @param barcode.column Column to use as proxy for barcodes ("nCount_RNA" by default)
#' @param group.column Column to group by ("orig.ident" by default)
#' @param threshold.high Ignore barcodes of rank above thisf threshold in inflection calculation
#' @param threshold.low Ignore barcodes of rank below this threshold in inflection calculation
#'
#' @return Returns Seurat object with a new list in the `tools` slot, `CalculateBarcodeInflections` with values:
#'
#' * `barcode_distribution` - contains the full barcode distribution across the entire dataset
#' * `inflection_points` - the calculated inflection points within the thresholds
#' * `threshold_values` - the provided (or default) threshold values to search within for inflections
#' * `cells_pass` - the cells that pass the inflection point calculation
#'
#' @importFrom methods slot
#' @importFrom stats ave aggregate
#'
#' @export
#' @concept preprocessing
#'
#' @author Robert A. Amezquita, \email{robert.amezquita@fredhutch.org}
#' @seealso \code{\link{BarcodeInflectionsPlot}} \code{\link{SubsetByBarcodeInflections}}
#'
#' @examples
#' data("pbmc_small")
#' CalculateBarcodeInflections(pbmc_small, group.column = 'groups')
#'
CalculateBarcodeInflections <- function(
  object,
  barcode.column = "nCount_RNA",
  group.column = "orig.ident",
  threshold.low = NULL,
  threshold.high = NULL
) {
  ## Check that barcode.column exists in meta.data
  if (!(barcode.column %in% colnames(x = object[[]]))) {
    stop("`barcode.column` specified not present in Seurat object provided")
  }
  # Calculation of barcode distribution
  ## Append rank by grouping x umi column
  # barcode_dist <- as.data.frame(object@meta.data)[, c(group.column, barcode.column)]
  barcode_dist <- object[[c(group.column, barcode.column)]]
  barcode_dist <- barcode_dist[do.call(what = order, args = barcode_dist), ] # order by columns left to right
  barcode_dist$rank <- ave(
    x = barcode_dist[, barcode.column], barcode_dist[, group.column],
    FUN = function(x) {
      return(rev(x = order(x)))
    }
  )
  barcode_dist <- barcode_dist[order(barcode_dist[, group.column], barcode_dist[, 'rank']), ]
  ## calculate rawdiff and append per group
  top <- aggregate(
    x = barcode_dist[, barcode.column],
    by = list(barcode_dist[, group.column]),
    FUN = function(x) {
      return(c(0, diff(x = log10(x = x + 1))))
    })$x
  bot <- aggregate(
    x = barcode_dist[, 'rank'],
    by = list(barcode_dist[, group.column]),
    FUN = function(x) {
      return(c(0, diff(x = x)))
    }
  )$x
  barcode_dist$rawdiff <- unlist(x = mapply(
    FUN = function(x, y) {
      return(ifelse(test = is.na(x = x / y), yes = 0, no = x / y))
    },
    x = top,
    y = bot
  ))
  # Calculation of inflection points
  ## Set thresholds for rank of barcodes to ignore
  threshold.low <- threshold.low %||% 1
  threshold.high <- threshold.high %||% max(barcode_dist$rank)
  ## Subset the barcode distribution by thresholds
  barcode_dist_sub <- barcode_dist[barcode_dist$rank > threshold.low & barcode_dist$rank < threshold.high, ]
  ## Calculate inflection points
  ## note: if thresholds are s.t. it produces the same length across both groups,
  ## aggregate will create a data.frame with x.* columns, where * is the length
  ## using the same combine approach will yield non-symmetrical results!
  whichmin_list <- aggregate(
    x = barcode_dist_sub[, 'rawdiff'],
    by = list(barcode_dist_sub[, group.column]),
    FUN = function(x) {
      return(x == min(x))
    }
  )$x
  ## workaround for aggregate behavior noted above
  if (is.list(x = whichmin_list)) { # uneven lengths
    is_inflection <- unlist(x = whichmin_list)
  } else if (is.matrix(x = whichmin_list)) { # even lengths
    is_inflection <- as.vector(x = t(x = whichmin_list))
  }
  tmp <- cbind(barcode_dist_sub, is_inflection)
  # inflections <- tmp[tmp$is_inflection == TRUE, c(group.column, barcode.column, "rank")]
  inflections <- tmp[which(x = tmp$is_inflection), c(group.column, barcode.column, 'rank')]
  # Use inflection point for what cells to keep
  ## use the inflection points to cut the subsetted dist to what to keep
  ## keep only the barcodes above the inflection points
  keep <- unlist(x = lapply(
    X = whichmin_list,
    FUN = function(x) {
      keep <- !x
      if (sum(keep) == length(x = keep)) {
        return(keep) # prevents bug in case of keeping all cells
      }
      # toss <- which(keep == FALSE):length(x = keep) # the end cells below knee
      toss <- which(x = !keep):length(x = keep)
      keep[toss] <- FALSE
      return(keep)
    }
  ))
  barcode_dist_sub_keep <- barcode_dist_sub[keep, ]
  cells_keep <- rownames(x = barcode_dist_sub_keep)
  # Bind thresholds to keep track of where they are placed
  thresholds <- data.frame(
    threshold = c('threshold.low', 'threshold.high'),
    rank = c(threshold.low, threshold.high)
  )
  # Combine relevant info together
  ## Combine Barcode dist, inflection point, and cells to keep into list
  info <- list(
    barcode_distribution = barcode_dist,
    inflection_points = inflections,
    threshold_values = thresholds,
    cells_pass = cells_keep
  )
  # save results into object
  Tool(object = object) <- info
  return(object)
}

#' Demultiplex samples based on data from cell 'hashing'
#'
#' Assign sample-of-origin for each cell, annotate doublets.
#'
#' @param object Seurat object. Assumes that the hash tag oligo (HTO) data has been added and normalized.
#' @param assay Name of the Hashtag assay (HTO by default)
#' @param positive.quantile The quantile of inferred 'negative' distribution for each hashtag - over which the cell is considered 'positive'. Default is 0.99
#' @param init Initial number of clusters for hashtags. Default is the # of hashtag oligo names + 1 (to account for negatives)
#' @param kfunc Clustering function for initial hashtag grouping. Default is "clara" for fast k-medoids clustering on large applications, also support "kmeans" for kmeans clustering
#' @param nsamples Number of samples to be drawn from the dataset used for clustering, for kfunc = "clara"
#' @param nstarts nstarts value for k-means clustering (for kfunc = "kmeans"). 100 by default
#' @param seed Sets the random seed. If NULL, seed is not set
#' @param verbose Prints the output
#'
#' @return The Seurat object with the following demultiplexed information stored in the meta data:
#' \describe{
#'   \item{hash.maxID}{Name of hashtag with the highest signal}
#'   \item{hash.secondID}{Name of hashtag with the second highest signal}
#'   \item{hash.margin}{The difference between signals for hash.maxID and hash.secondID}
#'   \item{classification}{Classification result, with doublets/multiplets named by the top two highest hashtags}
#'   \item{classification.global}{Global classification result (singlet, doublet or negative)}
#'   \item{hash.ID}{Classification result where doublet IDs are collapsed}
#' }
#'
#' @importFrom cluster clara
#' @importFrom Matrix colSums
#' @importFrom fitdistrplus fitdist
#' @importFrom stats pnbinom kmeans
#'
#' @export
#' @concept preprocessing
#'
#' @seealso \code{\link{HTOHeatmap}}
#'
#' @examples
#' \dontrun{
#' object <- HTODemux(object)
#' }
#'
HTODemux <- function(
  object,
  assay = "HTO",
  positive.quantile = 0.99,
  init = NULL,
  nstarts = 100,
  kfunc = "clara",
  nsamples = 100,
  seed = 42,
  verbose = TRUE
) {
  if (!is.null(x = seed)) {
    set.seed(seed = seed)
  }
  #initial clustering
  assay <- assay %||% DefaultAssay(object = object)
  data <- GetAssayData(object = object, assay = assay)
  counts <- GetAssayData(
    object = object,
    assay = assay,
    slot = 'counts'
  )[, colnames(x = object)]
  counts <- as.matrix(x = counts)
  ncenters <- init %||% (nrow(x = data) + 1)
  switch(
    EXPR = kfunc,
    'kmeans' = {
      init.clusters <- kmeans(
        x = t(x = GetAssayData(object = object, assay = assay)),
        centers = ncenters,
        nstart = nstarts
      )
      #identify positive and negative signals for all HTO
      Idents(object = object, cells = names(x = init.clusters$cluster)) <- init.clusters$cluster
    },
    'clara' = {
      #use fast k-medoid clustering
      init.clusters <- clara(
        x = t(x = GetAssayData(object = object, assay = assay)),
        k = ncenters,
        samples = nsamples
      )
      #identify positive and negative signals for all HTO
      Idents(object = object, cells = names(x = init.clusters$clustering), drop = TRUE) <- init.clusters$clustering
    },
    stop("Unknown k-means function ", kfunc, ", please choose from 'kmeans' or 'clara'")
  )
  #average hto signals per cluster
  #work around so we don't average all the RNA levels which takes time
  average.expression <- AverageExpression(
    object = object,
    assays = assay,
    verbose = FALSE
  )[[assay]]
  #checking for any cluster with all zero counts for any barcode
  if (sum(average.expression == 0) > 0) {
    stop("Cells with zero counts exist as a cluster.")
  }
  #create a matrix to store classification result
  discrete <- GetAssayData(object = object, assay = assay)
  discrete[discrete > 0] <- 0
  # for each HTO, we will use the minimum cluster for fitting
  for (iter in rownames(x = data)) {
    values <- counts[iter, colnames(object)]
    #commented out if we take all but the top cluster as background
    #values_negative=values[setdiff(object@cell.names,WhichCells(object,which.max(average.expression[iter,])))]
    values.use <- values[WhichCells(
      object = object,
      idents = levels(x = Idents(object = object))[[which.min(x = average.expression[iter, ])]]
    )]
    fit <- suppressWarnings(expr = fitdist(data = values.use, distr = "nbinom"))
    cutoff <- as.numeric(x = quantile(x = fit, probs = positive.quantile)$quantiles[1])
    discrete[iter, names(x = which(x = values > cutoff))] <- 1
    if (verbose) {
      message(paste0("Cutoff for ", iter, " : ", cutoff, " reads"))
    }
  }
  # now assign cells to HTO based on discretized values
  npositive <- colSums(x = discrete)
  classification.global <- npositive
  classification.global[npositive == 0] <- "Negative"
  classification.global[npositive == 1] <- "Singlet"
  classification.global[npositive > 1] <- "Doublet"
  donor.id = rownames(x = data)
  hash.max <- apply(X = data, MARGIN = 2, FUN = max)
  hash.maxID <- apply(X = data, MARGIN = 2, FUN = which.max)
  hash.second <- apply(X = data, MARGIN = 2, FUN = MaxN, N = 2)
  hash.maxID <- as.character(x = donor.id[sapply(
    X = 1:ncol(x = data),
    FUN = function(x) {
      return(which(x = data[, x] == hash.max[x])[1])
    }
  )])
  hash.secondID <- as.character(x = donor.id[sapply(
    X = 1:ncol(x = data),
    FUN = function(x) {
      return(which(x = data[, x] == hash.second[x])[1])
    }
  )])
  hash.margin <- hash.max - hash.second
  doublet_id <- sapply(
    X = 1:length(x = hash.maxID),
    FUN = function(x) {
      return(paste(sort(x = c(hash.maxID[x], hash.secondID[x])), collapse = "_"))
    }
  )
  # doublet_names <- names(x = table(doublet_id))[-1] # Not used
  classification <- classification.global
  classification[classification.global == "Negative"] <- "Negative"
  classification[classification.global == "Singlet"] <- hash.maxID[which(x = classification.global == "Singlet")]
  classification[classification.global == "Doublet"] <- doublet_id[which(x = classification.global == "Doublet")]
  classification.metadata <- data.frame(
    hash.maxID,
    hash.secondID,
    hash.margin,
    classification,
    classification.global
  )
  colnames(x = classification.metadata) <- paste(
    assay,
    c('maxID', 'secondID', 'margin', 'classification', 'classification.global'),
    sep = '_'
  )
  object <- AddMetaData(object = object, metadata = classification.metadata)
  Idents(object) <- paste0(assay, '_classification')
  # Idents(object, cells = rownames(object@meta.data[object@meta.data$classification.global == "Doublet", ])) <- "Doublet"
  doublets <- rownames(x = object[[]])[which(object[[paste0(assay, "_classification.global")]] == "Doublet")]
  Idents(object = object, cells = doublets) <- 'Doublet'
  # object@meta.data$hash.ID <- Idents(object)
  object$hash.ID <- Idents(object = object)
  return(object)
}

#' Calculate pearson residuals of features not in the scale.data
#'
#' This function calls sctransform::get_residuals.
#'
#' @param object A seurat object
#' @param features Name of features to add into the scale.data
#' @param assay Name of the assay of the seurat object generated by SCTransform
#' @param umi.assay Name of the assay of the seurat object containing UMI matrix
#' and the default is RNA
#' @param clip.range Numeric of length two specifying the min and max values the
#' Pearson residual will be clipped to
#' @param replace.value Recalculate residuals for all features, even if they are
#' already present. Useful if you want to change the clip.range.
#' @param na.rm For features where there is no feature model stored, return NA
#' for residual value in scale.data when na.rm = FALSE. When na.rm is TRUE, only
#' return residuals for features with a model stored for all cells.
#' @param verbose Whether to print messages and progress bars
#'
#' @return Returns a Seurat object containing Pearson residuals of added
#' features in its scale.data
#'
#' @importFrom sctransform get_residuals
#' @importFrom matrixStats rowAnyNAs
#'
#' @export
#' @concept preprocessing
#'
#' @seealso \code{\link[sctransform]{get_residuals}}
#'
#' @examples
#' data("pbmc_small")
#' pbmc_small <- SCTransform(object = pbmc_small, variable.features.n = 20)
#' pbmc_small <- GetResidual(object = pbmc_small, features = c('MS4A1', 'TCL1A'))
#'
GetResidual <- function(
  object,
  features,
  assay = NULL,
  umi.assay = "RNA",
  clip.range = NULL,
  replace.value = FALSE,
  na.rm = TRUE,
  verbose = TRUE
) {
  assay <- assay %||% DefaultAssay(object = object)
  if (IsSCT(assay = object[[assay]])) {
    object[[assay]] <- as(object[[assay]], 'SCTAssay')
  }
  if (!inherits(x = object[[assay]], what = "SCTAssay")) {
    stop(assay, " assay was not generated by SCTransform")
  }
  sct.models <- levels(x = object[[assay]])
  if (length(x = sct.models) == 0) {
    warning("SCT model not present in assay", call. = FALSE, immediate. = TRUE)
    return(object)
  }
  possible.features <- unique(x = unlist(x = lapply(X = sct.models, FUN = function(x) {
    rownames(x = SCTResults(object = object[[assay]], slot = "feature.attributes", model = x))
  }
  )))
  bad.features <- setdiff(x = features, y = possible.features)
  if (length(x = bad.features) > 0) {
    warning("The following requested features are not present in any models: ",
            paste(bad.features, collapse = ", "), call. = FALSE)
    features <- intersect(x = features, y = possible.features)
  }
  features.orig <- features
  if (na.rm) {
    # only compute residuals when feature model info is present in all
    features <- names(x = which(x = table(unlist(x = lapply(
      X = sct.models,
      FUN = function(x) {
        rownames(x = SCTResults(object = object[[assay]], slot = "feature.attributes", model = x))
      }
    ))) == length(x = sct.models)))
    if (length(x = features) == 0) {
      return(object)
    }
  }
  features <- intersect(x = features.orig, y = features)
  if (length(x = sct.models) > 1 && verbose) {
    message(
      "This SCTAssay contains multiple SCT models. Computing residuals for cells using different models"
    )
  }
  if (!umi.assay %in% Assays(object = object) || 
      length(x = Layers(object = object[[umi.assay]], search = 'counts')) == 0) {
    return(object)
  }
  if (inherits(x = object[[umi.assay]], what = 'Assay')) {
    new.residuals <- lapply(
      X = sct.models,
      FUN = function(x) {
        GetResidualSCTModel(
          object = object,
          assay = assay,
          SCTModel = x,
          new_features = features,
          replace.value = replace.value,
          clip.range = clip.range,
          verbose = verbose
        )
      }
    )
  } else if (inherits(x = object[[umi.assay]], what = 'Assay5')) {
    new.residuals <- lapply(
      X = sct.models,
      FUN = function(x) {
        FetchResidualSCTModel(object = object,
                              assay = assay,
                              umi.assay = umi.assay,
                              SCTModel = x,
                              new_features = features,
                              replace.value = replace.value,
                              clip.range = clip.range,
                              verbose = verbose)
      }
    )
  }
  existing.data <- GetAssayData(object = object, slot = 'scale.data', assay = assay)
  all.features <- union(x = rownames(x = existing.data), y = features)
   new.scale <- matrix(
    data = NA,
    nrow = length(x = all.features),
    ncol = ncol(x = object),
    dimnames = list(all.features, Cells(x = object))
  )
  if (nrow(x = existing.data) > 0){
    new.scale[1:nrow(x = existing.data), ] <- existing.data
  }
  if (length(x = new.residuals) == 1 & is.list(x = new.residuals)) {
    new.residuals <- new.residuals[[1]]
  } else {
    new.residuals <- Reduce(cbind, new.residuals)
  }
  new.scale[rownames(x = new.residuals), colnames(x = new.residuals)] <- new.residuals
  if (na.rm) {
    new.scale <- new.scale[!rowAnyNAs(x = new.scale), ]
  }
  object <- SetAssayData(
    object = object,
    assay = assay,
    slot = "scale.data",
    new.data = new.scale
  )
  if (any(!features.orig %in% rownames(x = new.scale))) {
    bad.features <- features.orig[which(!features.orig %in% rownames(x = new.scale))]
    warning("Residuals not computed for the following requested features: ",
            paste(bad.features, collapse = ", "), call. = FALSE)
  }
  return(object)
}

#' Load a 10x Genomics Visium Spatial Experiment into a \code{Seurat} object
#'
#' @inheritParams Read10X
#' @inheritParams SeuratObject::CreateSeuratObject
#' @param data.dir Directory containing the H5 file specified by \code{filename}
#' and the image data in a subdirectory called \code{spatial}
#' @param filename Name of H5 file containing the feature barcode matrix
#' @param slice Name for the stored image of the tissue slice
#' @param filter.matrix Only keep spots that have been determined to be over
#' tissue
#' @param to.upper Converts all feature names to upper case. Can be useful when
#' analyses require comparisons between human and mouse gene names for example.
#' @param ... Arguments passed to \code{\link{Read10X_h5}}
#'
#' @return A \code{Seurat} object
#'
#' @importFrom png readPNG
#' @importFrom grid rasterGrob
#' @importFrom jsonlite fromJSON
#'
#' @export
#' @concept preprocessing
#'
#' @examples
#' \dontrun{
#' data_dir <- 'path/to/data/directory'
#' list.files(data_dir) # Should show filtered_feature_bc_matrix.h5
#' Load10X_Spatial(data.dir = data_dir)
#' }
#'
Load10X_Spatial <- function(
  data.dir,
  filename = 'filtered_feature_bc_matrix.h5',
  assay = 'Spatial',
  slice = 'slice1',
  filter.matrix = TRUE,
  to.upper = FALSE,
  ...
) {
  if (length(x = data.dir) > 1) {
    warning("'Load10X_Spatial' accepts only one 'data.dir'", immediate. = TRUE)
    data.dir <- data.dir[1]
  }
  data <- Read10X_h5(filename = file.path(data.dir, filename), ...)
  if (to.upper) {
    rownames(x = data) <- toupper(x = rownames(x = data))
  }
  object <- CreateSeuratObject(counts = data, assay = assay)
  image <- Read10X_Image(
    image.dir = file.path(data.dir, 'spatial'),
    filter.matrix = filter.matrix
  )
  image <- image[Cells(x = object)]
  DefaultAssay(object = image) <- assay
  object[[slice]] <- image
  return(object)
}

#' Load STARmap data
#'
#' @param data.dir location of data directory that contains the counts matrix,
#' gene name, qhull, and centroid files.
#' @param counts.file name of file containing the counts matrix (csv)
#' @param gene.file name of file containing the gene names (csv)
#' @param qhull.file name of file containing the hull coordinates (tsv)
#' @param centroid.file name of file containing the centroid positions (tsv)
#' @param assay Name of assay to associate spatial data to
#' @param image Name of "image" object storing spatial coordinates
#'
#' @return A \code{\link{Seurat}} object
#'
#' @importFrom methods new
#' @importFrom utils read.csv read.table
#'
#' @seealso \code{\link{STARmap}}
#'
#' @export
#' @concept preprocessing
#'
LoadSTARmap <- function(
  data.dir,
  counts.file = "cell_barcode_count.csv",
  gene.file = "genes.csv",
  qhull.file = "qhulls.tsv",
  centroid.file = "centroids.tsv",
  assay = "Spatial",
  image = "image"
) {
  if (!dir.exists(paths = data.dir)) {
    stop("Cannot find directory ", data.dir, call. = FALSE)
  }
  counts <- read.csv(
    file = file.path(data.dir, counts.file),
    as.is = TRUE,
    header = FALSE
  )
  gene.names <- read.csv(
    file = file.path(data.dir, gene.file),
    as.is = TRUE,
    header = FALSE
  )
  qhulls <- read.table(
    file = file.path(data.dir, qhull.file),
    sep = '\t',
    col.names = c('cell', 'y', 'x'),
    as.is = TRUE
  )
  centroids <- read.table(
    file = file.path(data.dir, centroid.file),
    sep = '\t',
    as.is = TRUE,
    col.names = c('y', 'x')
  )
  colnames(x = counts) <- gene.names[, 1]
  rownames(x = counts) <- paste0('starmap', seq(1:nrow(x = counts)))
  counts <- as.matrix(x = counts)
  rownames(x = centroids) <- rownames(x = counts)
  qhulls$cell <- paste0('starmap', qhulls$cell)
  centroids <- as.matrix(x = centroids)
  starmap <- CreateSeuratObject(counts = t(x = counts), assay = assay)
  starmap[[image]] <- new(
    Class = 'STARmap',
    assay = assay,
    coordinates = as.data.frame(x = centroids),
    qhulls = qhulls
  )
  return(starmap)
}

#' Demultiplex samples based on classification method from MULTI-seq (McGinnis et al., bioRxiv 2018)
#'
#' Identify singlets, doublets and negative cells from multiplexing experiments. Annotate singlets by tags.
#'
#' @param object Seurat object. Assumes that the specified assay data has been added
#' @param assay Name of the multiplexing assay (HTO by default)
#' @param quantile The quantile to use for classification
#' @param autoThresh Whether to perform automated threshold finding to define the best quantile. Default is FALSE
#' @param maxiter Maximum number of iterations if autoThresh = TRUE. Default is 5
#' @param qrange A range of possible quantile values to try if autoThresh = TRUE
#' @param verbose Prints the output
#'
#' @return A Seurat object with demultiplexing results stored at \code{object$MULTI_ID}
#'
#' @export
#' @concept preprocessing
#'
#' @references \url{https://www.biorxiv.org/content/10.1101/387241v1}
#'
#' @examples
#' \dontrun{
#' object <- MULTIseqDemux(object)
#' }
#'
MULTIseqDemux <- function(
  object,
  assay = "HTO",
  quantile = 0.7,
  autoThresh = FALSE,
  maxiter = 5,
  qrange = seq(from = 0.1, to = 0.9, by = 0.05),
  verbose = TRUE
) {
  assay <- assay %||% DefaultAssay(object = object)
  multi_data_norm <- t(x = GetAssayData(
    object = object,
    slot = "data",
    assay = assay
  ))
  if (autoThresh) {
    iter <- 1
    negatives <- c()
    neg.vector <- c()
    while (iter <= maxiter) {
      # Iterate over q values to find ideal barcode thresholding results by maximizing singlet classifications
      bar.table_sweep.list <- list()
      n <- 0
      for (q in qrange) {
        n <- n + 1
        # Generate list of singlet/doublet/negative classifications across q sweep
        bar.table_sweep.list[[n]] <- ClassifyCells(data = multi_data_norm, q = q)
        names(x = bar.table_sweep.list)[n] <- paste0("q=" , q)
      }
      # Determine which q values results in the highest pSinglet
      res_round <- FindThresh(call.list = bar.table_sweep.list)$res
      res.use <- res_round[res_round$Subset == "pSinglet", ]
      q.use <- res.use[which.max(res.use$Proportion),"q"]
      if (verbose) {
        message("Iteration ", iter)
        message("Using quantile ", q.use)
      }
      round.calls <- ClassifyCells(data = multi_data_norm, q = q.use)
      #remove negative cells
      neg.cells <- names(x = round.calls)[which(x = round.calls == "Negative")]
      neg.vector <- c(neg.vector, rep(x = "Negative", length(x = neg.cells)))
      negatives <- c(negatives, neg.cells)
      if (length(x = neg.cells) == 0) {
        break
      }
      multi_data_norm <- multi_data_norm[-which(x = rownames(x = multi_data_norm) %in% neg.cells), ]
      iter <- iter + 1
    }
    names(x = neg.vector) <- negatives
    demux_result <- c(round.calls,neg.vector)
    demux_result <- demux_result[rownames(x = object[[]])]
  } else{
    demux_result <- ClassifyCells(data = multi_data_norm, q = quantile)
  }
  demux_result <- demux_result[rownames(x = object[[]])]
  object[['MULTI_ID']] <- factor(x = demux_result)
  Idents(object = object) <- "MULTI_ID"
  bcs <- colnames(x = multi_data_norm)
  bc.max <- bcs[apply(X = multi_data_norm, MARGIN = 1, FUN = which.max)]
  bc.second <- bcs[unlist(x = apply(
    X = multi_data_norm,
    MARGIN = 1,
    FUN = function(x) {
      return(which(x == MaxN(x)))
    }
  ))]
  doublet.names <- unlist(x = lapply(
    X = 1:length(x = bc.max),
    FUN = function(x) {
      return(paste(sort(x = c(bc.max[x], bc.second[x])), collapse =  "_"))
    }
  ))
  doublet.id <- which(x = demux_result == "Doublet")
  MULTI_classification <- as.character(object$MULTI_ID)
  MULTI_classification[doublet.id] <- doublet.names[doublet.id]
  object$MULTI_classification <- factor(x = MULTI_classification)
  return(object)
}

#' Load in data from 10X
#'
#' Enables easy loading of sparse data matrices provided by 10X genomics.
#'
#' @param data.dir Directory containing the matrix.mtx, genes.tsv (or features.tsv), and barcodes.tsv
#' files provided by 10X. A vector or named vector can be given in order to load
#' several data directories. If a named vector is given, the cell barcode names
#' will be prefixed with the name.
#' @param gene.column Specify which column of genes.tsv or features.tsv to use for gene names; default is 2
#' @param cell.column Specify which column of barcodes.tsv to use for cell names; default is 1
#' @param unique.features Make feature names unique (default TRUE)
#' @param strip.suffix Remove trailing "-1" if present in all cell barcodes.
#'
#' @return If features.csv indicates the data has multiple data types, a list
#'   containing a sparse matrix of the data from each type will be returned.
#'   Otherwise a sparse matrix containing the expression data will be returned.
#'
#' @importFrom Matrix readMM
#' @importFrom utils read.delim
#'
#' @export
#' @concept preprocessing
#'
#' @examples
#' \dontrun{
#' # For output from CellRanger < 3.0
#' data_dir <- 'path/to/data/directory'
#' list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
#' expression_matrix <- Read10X(data.dir = data_dir)
#' seurat_object = CreateSeuratObject(counts = expression_matrix)
#'
#' # For output from CellRanger >= 3.0 with multiple data types
#' data_dir <- 'path/to/data/directory'
#' list.files(data_dir) # Should show barcodes.tsv.gz, features.tsv.gz, and matrix.mtx.gz
#' data <- Read10X(data.dir = data_dir)
#' seurat_object = CreateSeuratObject(counts = data$`Gene Expression`)
#' seurat_object[['Protein']] = CreateAssayObject(counts = data$`Antibody Capture`)
#' }
#'
Read10X <- function(
  data.dir,
  gene.column = 2,
  cell.column = 1,
  unique.features = TRUE,
  strip.suffix = FALSE
) {
  full.data <- list()
  has_dt <- requireNamespace("data.table", quietly = TRUE) && requireNamespace("R.utils", quietly = TRUE)
  for (i in seq_along(along.with = data.dir)) {
    run <- data.dir[i]
    if (!dir.exists(paths = run)) {
      stop("Directory provided does not exist")
    }
    barcode.loc <- file.path(run, 'barcodes.tsv')
    gene.loc <- file.path(run, 'genes.tsv')
    features.loc <- file.path(run, 'features.tsv.gz')
    matrix.loc <- file.path(run, 'matrix.mtx')
    # Flag to indicate if this data is from CellRanger >= 3.0
    pre_ver_3 <- file.exists(gene.loc)
    if (!pre_ver_3) {
      addgz <- function(s) {
        return(paste0(s, ".gz"))
      }
      barcode.loc <- addgz(s = barcode.loc)
      matrix.loc <- addgz(s = matrix.loc)
    }
    if (!file.exists(barcode.loc)) {
      stop("Barcode file missing. Expecting ", basename(path = barcode.loc))
    }
    if (!pre_ver_3 && !file.exists(features.loc) ) {
      stop("Gene name or features file missing. Expecting ", basename(path = features.loc))
    }
    if (!file.exists(matrix.loc)) {
      stop("Expression matrix file missing. Expecting ", basename(path = matrix.loc))
    }
    data <- readMM(file = matrix.loc)
    if (has_dt) {
      cell.barcodes <- as.data.frame(data.table::fread(barcode.loc, header = FALSE))
    } else {
      cell.barcodes <- read.table(file = barcode.loc, header = FALSE, sep = '\t', row.names = NULL)
    }

    if (ncol(x = cell.barcodes) > 1) {
      cell.names <- cell.barcodes[, cell.column]
    } else {
      cell.names <- readLines(con = barcode.loc)
    }
    if (all(grepl(pattern = "\\-1$", x = cell.names)) & strip.suffix) {
      cell.names <- as.vector(x = as.character(x = sapply(
        X = cell.names,
        FUN = ExtractField,
        field = 1,
        delim = "-"
      )))
    }
    if (is.null(x = names(x = data.dir))) {
      if (length(x = data.dir) < 2) {
        colnames(x = data) <- cell.names
      } else {
        colnames(x = data) <- paste0(i, "_", cell.names)
      }
    } else {
      colnames(x = data) <- paste0(names(x = data.dir)[i], "_", cell.names)
    }

    if (has_dt) {
      feature.names <- as.data.frame(data.table::fread(ifelse(test = pre_ver_3, yes = gene.loc, no = features.loc), header = FALSE))
    } else {
      feature.names <- read.delim(
        file = ifelse(test = pre_ver_3, yes = gene.loc, no = features.loc),
        header = FALSE,
        stringsAsFactors = FALSE
      )
    }

    if (any(is.na(x = feature.names[, gene.column]))) {
      warning(
        'Some features names are NA. Replacing NA names with ID from the opposite column requested',
        call. = FALSE,
        immediate. = TRUE
      )
      na.features <- which(x = is.na(x = feature.names[, gene.column]))
      replacement.column <- ifelse(test = gene.column == 2, yes = 1, no = 2)
      feature.names[na.features, gene.column] <- feature.names[na.features, replacement.column]
    }
    if (unique.features) {
      fcols = ncol(x = feature.names)
      if (fcols < gene.column) {
        stop(paste0("gene.column was set to ", gene.column,
                    " but feature.tsv.gz (or genes.tsv) only has ", fcols, " columns.",
                    " Try setting the gene.column argument to a value <= to ", fcols, "."))
      }
      rownames(x = data) <- make.unique(names = feature.names[, gene.column])
    }
    # In cell ranger 3.0, a third column specifying the type of data was added
    # and we will return each type of data as a separate matrix
    if (ncol(x = feature.names) > 2) {
      data_types <- factor(x = feature.names$V3)
      lvls <- levels(x = data_types)
      if (length(x = lvls) > 1 && length(x = full.data) == 0) {
        message("10X data contains more than one type and is being returned as a list containing matrices of each type.")
      }
      expr_name <- "Gene Expression"
      if (expr_name %in% lvls) { # Return Gene Expression first
        lvls <- c(expr_name, lvls[-which(x = lvls == expr_name)])
      }
      data <- lapply(
        X = lvls,
        FUN = function(l) {
          return(data[data_types == l, , drop = FALSE])
        }
      )
      names(x = data) <- lvls
    } else{
      data <- list(data)
    }
    full.data[[length(x = full.data) + 1]] <- data
  }
  # Combine all the data from different directories into one big matrix, note this
  # assumes that all data directories essentially have the same features files
  list_of_data <- list()
  for (j in 1:length(x = full.data[[1]])) {
    list_of_data[[j]] <- do.call(cbind, lapply(X = full.data, FUN = `[[`, j))
    # Fix for Issue #913
    list_of_data[[j]] <- as.sparse(x = list_of_data[[j]])
  }
  names(x = list_of_data) <- names(x = full.data[[1]])
  # If multiple features, will return a list, otherwise
  # a matrix.
  if (length(x = list_of_data) == 1) {
    return(list_of_data[[1]])
  } else {
    return(list_of_data)
  }
}

#' Read 10X hdf5 file
#'
#' Read count matrix from 10X CellRanger hdf5 file.
#' This can be used to read both scATAC-seq and scRNA-seq matrices.
#'
#' @param filename Path to h5 file
#' @param use.names Label row names with feature names rather than ID numbers.
#' @param unique.features Make feature names unique (default TRUE)
#'
#' @return Returns a sparse matrix with rows and columns labeled. If multiple
#' genomes are present, returns a list of sparse matrices (one per genome).
#'
#' @export
#' @concept preprocessing
#'
Read10X_h5 <- function(filename, use.names = TRUE, unique.features = TRUE) {
  if (!requireNamespace('hdf5r', quietly = TRUE)) {
    stop("Please install hdf5r to read HDF5 files")
  }
  if (!file.exists(filename)) {
    stop("File not found")
  }
  infile <- hdf5r::H5File$new(filename = filename, mode = 'r')
  genomes <- names(x = infile)
  output <- list()
  if (hdf5r::existsGroup(infile, 'matrix')) {
    # cellranger version 3
    if (use.names) {
      feature_slot <- 'features/name'
    } else {
      feature_slot <- 'features/id'
    }
  } else {
    if (use.names) {
      feature_slot <- 'gene_names'
    } else {
      feature_slot <- 'genes'
    }
  }
  for (genome in genomes) {
    counts <- infile[[paste0(genome, '/data')]]
    indices <- infile[[paste0(genome, '/indices')]]
    indptr <- infile[[paste0(genome, '/indptr')]]
    shp <- infile[[paste0(genome, '/shape')]]
    features <- infile[[paste0(genome, '/', feature_slot)]][]
    barcodes <- infile[[paste0(genome, '/barcodes')]]
    sparse.mat <- sparseMatrix(
      i = indices[] + 1,
      p = indptr[],
      x = as.numeric(x = counts[]),
      dims = shp[],
      repr = "T"
    )
    if (unique.features) {
      features <- make.unique(names = features)
    }
    rownames(x = sparse.mat) <- features
    colnames(x = sparse.mat) <- barcodes[]
    sparse.mat <- as.sparse(x = sparse.mat)
    # Split v3 multimodal
    if (infile$exists(name = paste0(genome, '/features'))) {
      types <- infile[[paste0(genome, '/features/feature_type')]][]
      types.unique <- unique(x = types)
      if (length(x = types.unique) > 1) {
        message(
          "Genome ",
          genome,
          " has multiple modalities, returning a list of matrices for this genome"
        )
        sparse.mat <- sapply(
          X = types.unique,
          FUN = function(x) {
            return(sparse.mat[which(x = types == x), ])
          },
          simplify = FALSE,
          USE.NAMES = TRUE
        )
      }
    }
    output[[genome]] <- sparse.mat
  }
  infile$close_all()
  if (length(x = output) == 1) {
    return(output[[genome]])
  } else{
    return(output)
  }
}

#' Load a 10X Genomics Visium Image
#'
#' @param image.dir Path to directory with 10X Genomics visium image data;
#' should include files \code{tissue_lowres_iamge.png},
#' \code{scalefactors_json.json} and \code{tissue_positions_list.csv}
#' @param filter.matrix Filter spot/feature matrix to only include spots that
#' have been determined to be over tissue.
#' @param ... Ignored for now
#'
#' @return A \code{\link{VisiumV1}} object
#'
#' @importFrom png readPNG
#' @importFrom jsonlite fromJSON
#'
#' @seealso \code{\link{VisiumV1}} \code{\link{Load10X_Spatial}}
#'
#' @export
#' @concept preprocessing
#'
Read10X_Image <- function(image.dir, filter.matrix = TRUE, ...) {
  image <- readPNG(source = file.path(image.dir, 'tissue_lowres_image.png'))
  scale.factors <- fromJSON(txt = file.path(image.dir, 'scalefactors_json.json'))
  tissue.positions.path <- Sys.glob(paths = file.path(image.dir, 'tissue_positions*'))
  tissue.positions <- read.csv(
    file = tissue.positions.path,
    col.names = c('barcodes', 'tissue', 'row', 'col', 'imagerow', 'imagecol'),
    header = ifelse(
      test = basename(tissue.positions.path) == "tissue_positions.csv",
      yes = TRUE,
      no = FALSE
    ),
    as.is = TRUE,
    row.names = 1
  )
  if (filter.matrix) {
    tissue.positions <- tissue.positions[which(x = tissue.positions$tissue == 1), , drop = FALSE]
  }
  unnormalized.radius <- scale.factors$fiducial_diameter_fullres * scale.factors$tissue_lowres_scalef
  spot.radius <-  unnormalized.radius / max(dim(x = image))
  return(new(
    Class = 'VisiumV1',
    image = image,
    scale.factors = scalefactors(
      spot = scale.factors$tissue_hires_scalef,
      fiducial = scale.factors$fiducial_diameter_fullres,
      hires = scale.factors$tissue_hires_scalef,
      scale.factors$tissue_lowres_scalef
    ),
    coordinates = tissue.positions,
    spot.radius = spot.radius
  ))
}

#' Read and Load Akoya CODEX data
#'
#' @param filename Path to matrix generated by upstream processing.
#' @param type Specify which type matrix is being provided.
#' \itemize{
#'  \item \dQuote{\code{processor}}: matrix generated by CODEX Processor
#'  \item \dQuote{\code{inform}}: matrix generated by inForm
#'  \item \dQuote{\code{qupath}}: matrix generated by QuPath
#' }
#' @param filter A pattern to filter features by; pass \code{NA} to
#' skip feature filtering
#' @param inform.quant When \code{type} is \dQuote{\code{inform}}, the
#' quantification level to read in
#'
#' @return \code{ReadAkoya}: A list with some combination of the following values
#' \itemize{
#'  \item \dQuote{\code{matrix}}: a
#'  \link[Matrix:dgCMatrix-class]{sparse matrix} with expression data; cells
#'   are columns and features are rows
#'  \item \dQuote{\code{centroids}}: a data frame with cell centroid
#'   coordinates in three columns: \dQuote{x}, \dQuote{y}, and \dQuote{cell}
#'  \item \dQuote{\code{metadata}}: a data frame with cell-level meta data;
#'   includes all columns in \code{filename} that aren't in
#'   \dQuote{\code{matrix}} or \dQuote{\code{centroids}}
#' }
#' When \code{type} is \dQuote{\code{inform}}, additional expression matrices
#' are returned and named using their segmentation type (eg.
#' \dQuote{nucleus}, \dQuote{membrane}). The \dQuote{Entire Cell} segmentation
#' type is returned in the \dQuote{\code{matrix}} entry of the list
#'
#' @export
#'
#' @order 1
#'
#' @concept preprocessing
#'
#' @template section-progressr
#'
#' @templateVar pkg data.table
#' @template note-reqdpkg
#'
ReadAkoya <- function(
  filename,
  type = c('inform', 'processor', 'qupath'),
  filter = 'DAPI|Blank|Empty',
  inform.quant = c('mean', 'total', 'min', 'max', 'std')
) {
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Please install 'data.table' for this function")
  }
  # Check arguments
  if (!file.exists(filename)) {
    stop(paste("Can't file file:", filename))
  }
  type <- tolower(x = type[1L])
  type <- match.arg(arg = type)
  # outs <- list(matrix = NULL, centroids = NULL)
  ratio <- getOption(x = 'Seurat.input.sparse_ratio', default = 0.4)
  p <- progressor()
  # Preload matrix
  p(message = "Preloading Akoya matrix", class = 'sticky', amount = 0)
  sep <- switch(EXPR = type, 'inform' = '\t', ',')
  mtx <- data.table::fread(
    file = filename,
    sep = sep,
    data.table = FALSE,
    verbose = FALSE
  )
  # Assemble outputs
  p(
    message = paste0("Parsing matrix in '", type, "' format"),
    class = 'sticky',
    amount = 0
  )
  outs <- switch(
    EXPR = type,
    'processor' = {
      # Create centroids data frame
      p(
        message = 'Creating centroids coordinates',
        class = 'sticky',
        amount = 0
      )
      centroids <- data.frame(
        x = mtx[['x:x']],
        y = mtx[['y:y']],
        cell = as.character(x = mtx[['cell_id:cell_id']]),
        stringsAsFactors = FALSE
      )
      rownames(x = mtx) <- as.character(x = mtx[['cell_id:cell_id']])
      # Create metadata data frame
      p(message = 'Creating meta data', class = 'sticky', amount = 0)
      md <- mtx[, !grepl(pattern = '^cyc', x = colnames(x = mtx)), drop = FALSE]
      colnames(x = md) <- vapply(
        X = strsplit(x = colnames(x = md), split = ':'),
        FUN = '[[',
        FUN.VALUE = character(length = 1L),
        2L
      )
      # Create expression matrix
      p(message = 'Creating expression matrix', class = 'sticky', amount = 0)
      mtx <- mtx[, grepl(pattern = '^cyc', x = colnames(x = mtx)), drop = FALSE]
      colnames(x = mtx) <- vapply(
        X = strsplit(x = colnames(x = mtx), split = ':'),
        FUN = '[[',
        FUN.VALUE = character(length = 1L),
        2L
      )
      if (!is.na(x = filter)) {
        p(
          message = paste0("Filtering features with pattern '", filter, "'"),
          class = 'sticky',
          amount = 0
        )
        mtx <- mtx[, !grepl(pattern = filter, x = colnames(x = mtx)), drop = FALSE]
      }
      mtx <- t(x = mtx)
      if ((sum(mtx == 0) / length(x = mtx)) > ratio) {
        p(
          message = 'Converting expression to sparse matrix',
          class = 'sticky',
          amount = 0
        )
        mtx <- as.sparse(x = mtx)
      }
      list(matrix = mtx, centroids = centroids, metadata = md)
    },
    'inform' = {
      inform.quant <- tolower(x = inform.quant[1L])
      inform.quant <- match.arg(arg = inform.quant)
      expr.key <- c(
        mean = 'Mean',
        total = 'Total',
        min = 'Min',
        max = 'Max',
        std = 'Std Dev'
      )[inform.quant]
      expr.pattern <- '\\(Normalized Counts, Total Weighting\\)'
      rownames(x = mtx) <- mtx[['Cell ID']]
      mtx <- mtx[, setdiff(x = colnames(x = mtx), y = 'Cell ID'), drop = FALSE]
      # Create centroids
      p(
        message = 'Creating centroids coordinates',
        class = 'sticky',
        amount = 0
      )
      centroids <- data.frame(
        x = mtx[['Cell X Position']],
        y = mtx[['Cell Y Position']],
        cell  = rownames(x = mtx),
        stringsAsFactors = FALSE
      )
      # Create metadata
      p(message = 'Creating meta data', class = 'sticky', amount = 0)
      cols <- setdiff(
        x = grep(
          pattern = expr.pattern,
          x = colnames(x = mtx),
          value = TRUE,
          invert = TRUE
        ),
        y = paste('Cell', c('X', 'Y'), 'Position')
      )
      md <- mtx[, cols, drop = FALSE]
      # Create expression matrices
      exprs <- data.frame(
        cols = grep(
          pattern = paste(expr.key, expr.pattern),
          x = colnames(x = mtx),
          value = TRUE
        )
      )
      exprs$feature <- vapply(
        X = trimws(x = gsub(
          pattern = paste(expr.key, expr.pattern),
          replacement = '',
          x = exprs$cols
        )),
        FUN = function(x) {
          x <- unlist(x = strsplit(x = x, split = ' '))
          x <- x[length(x = x)]
          return(gsub(pattern = '\\(|\\)', replacement = '', x = x))
        },
        FUN.VALUE = character(length = 1L)
      )
      exprs$class <- tolower(x = vapply(
        X = strsplit(x = exprs$cols, split = ' '),
        FUN = '[[',
        FUN.VALUE = character(length = 1L),
        1L
      ))
      classes <- unique(x = exprs$class)
      outs <- vector(
        mode = 'list',
        length = length(x = classes) + 2L
      )
      names(x = outs) <- c(
        'matrix',
        'centroids',
        'metadata',
        setdiff(x = classes, y = 'entire')
      )
      outs$centroids <- centroids
      outs$metadata <- md
      # browser()
      for (i in classes) {
        p(
          message = paste(
            'Creating',
            switch(EXPR = i, 'entire' = 'entire cell', i),
            'expression matrix'
          ),
          class = 'sticky',
          amount = 0
        )
        df <- exprs[exprs$class == i, , drop = FALSE]
        expr <- mtx[, df$cols]
        colnames(x = expr) <- df$feature
        if (!is.na(x = filter)) {
          p(
            message = paste0("Filtering features with pattern '", filter, "'"),
            class = 'sticky',
            amount = 0
          )
          expr <- expr[, !grepl(pattern = filter, x = colnames(x = expr)), drop = FALSE]
        }
        expr <- t(x = expr)
        if ((sum(expr == 0, na.rm = TRUE) / length(x = expr)) > ratio) {
          p(
            message = paste(
              'Converting',
              switch(EXPR = i, 'entire' = 'entire cell', i),
              'expression to sparse matrix'
            ),
            class = 'sticky',
            amount = 0
          )
          expr <- as.sparse(x = expr)
        }
        outs[[switch(EXPR = i, 'entire' = 'matrix', i)]] <- expr
      }
      outs
    },
    'qupath' = {
      rownames(x = mtx) <- as.character(x = seq_len(length.out = nrow(x = mtx)))
      # Create centroids
      p(
        message = 'Creating centroids coordinates',
        class = 'sticky',
        amount = 0
      )
      xpos <- sort(
        x = grep(pattern = 'Centroid X', x = colnames(x = mtx), value = TRUE),
        decreasing = TRUE
      )[1L]
      ypos <- sort(
        x = grep(pattern = 'Centroid Y', x = colnames(x = mtx), value = TRUE),
        decreasing = TRUE
      )[1L]
      centroids <- data.frame(
        x = mtx[[xpos]],
        y = mtx[[ypos]],
        cell = rownames(x = mtx),
        stringsAsFactors = FALSE
      )
      # Create metadata
      p(message = 'Creating meta data', class = 'sticky', amount = 0)
      cols <- setdiff(
        x = grep(
          pattern = 'Cell: Mean',
          x = colnames(x = mtx),
          ignore.case = TRUE,
          value = TRUE,
          invert = TRUE
        ),
        y = c(xpos, ypos)
      )
      md <- mtx[, cols, drop = FALSE]
      # Create expression matrix
      p(message = 'Creating expression matrix', class = 'sticky', amount = 0)
      idx <- which(x = grepl(
        pattern = 'Cell: Mean',
        x = colnames(x = mtx),
        ignore.case = TRUE
      ))
      mtx <- mtx[, idx, drop = FALSE]
      colnames(x = mtx) <- vapply(
        X = strsplit(x = colnames(x = mtx), split = ':'),
        FUN = '[[',
        FUN.VALUE = character(length = 1L),
        1L
      )
      if (!is.na(x = filter)) {
        p(
          message = paste0("Filtering features with pattern '", filter, "'"),
          class = 'sticky',
          amount = 0
        )
        mtx <- mtx[, !grepl(pattern = filter, x = colnames(x = mtx)), drop = FALSE]
      }
      mtx <- t(x = mtx)
      if ((sum(mtx == 0) / length(x = mtx)) > ratio) {
        p(
          message = 'Converting expression to sparse matrix',
          class = 'sticky',
          amount = 0
        )
        mtx <- as.sparse(x = mtx)
      }
      list(matrix = mtx, centroids = centroids, metadata = md)
    },
    stop("Unknown matrix type: ", type)
  )
  return(outs)
}

#' Load in data from remote or local mtx files
#'
#' Enables easy loading of sparse data matrices
#'
#' @param mtx Name or remote URL of the mtx file
#' @param cells Name or remote URL of the cells/barcodes file
#' @param features Name or remote URL of the features/genes file
#' @param cell.column Specify which column of cells file to use for cell names; default is 1
#' @param feature.column Specify which column of features files to use for feature/gene names; default is 2
#' @param skip.cell Number of lines to skip in the cells file before beginning to read cell names
#' @param skip.feature Number of lines to skip in the features file before beginning to gene names
#' @param unique.features Make feature names unique (default TRUE)
#' @param strip.suffix Remove trailing "-1" if present in all cell barcodes.
#'
#' @return A sparse matrix containing the expression data.
#'
#' @importFrom Matrix readMM
#' @importFrom utils read.delim
#' @importFrom httr build_url parse_url
#' @importFrom tools file_ext
#'
#'
#' @export
#' @concept preprocessing
#'
#' @examples
#' \dontrun{
#' # For local files:
#'
#' expression_matrix <- ReadMtx(
#'   mtx = "count_matrix.mtx.gz", features = "features.tsv.gz",
#'   cells = "barcodes.tsv.gz"
#' )
#' seurat_object <- CreateSeuratObject(counts = expression_matrix)
#'
#' # For remote files:
#'
#' expression_matrix <- ReadMtx(mtx = "http://localhost/matrix.mtx",
#' cells = "http://localhost/barcodes.tsv",
#' features = "http://localhost/genes.tsv")
#' seurat_object <- CreateSeuratObject(counts = data)
#' }
#'
ReadMtx <- function(
  mtx,
  cells,
  features,
  cell.column = 1,
  feature.column = 2,
  skip.cell = 0,
  skip.feature = 0,
  unique.features = TRUE,
  strip.suffix = FALSE
) {
  all.files <- list(
    "expression matrix" = mtx,
    "barcode list" = cells,
    "feature list" = features
  )
  for (i in seq_along(along.with = all.files)) {
    uri <- tryCatch(
      expr = {
        con <- url(description = all.files[[i]])
        close(con = con)
        all.files[[i]]
      },
      error = function(...) {
        return(normalizePath(path = all.files[[i]], winslash = '/'))
      }
    )
    err <- paste("Cannot find", names(x = all.files)[i], "at", uri)
    uri <- build_url(url = parse_url(url = uri))
    if (grepl(pattern = '^[A-Z]?:///', x = uri)) {
      uri <- gsub(pattern = '^://', replacement = '', x = uri)
      if (!file.exists(uri)) {
        stop(err, call. = FALSE)
      }
    } else {
      if (!Online(url = uri, seconds = 2L)) {
        stop(err, call. = FALSE)
      }
      if (file_ext(uri) == 'gz') {
        con <- url(description = uri)
        uri <- gzcon(con = con, text = TRUE)
      }
    }
    all.files[[i]] <- uri
  }
  cell.barcodes <- read.table(
    file = all.files[['barcode list']],
    header = FALSE,
    sep = '\t',
    row.names = NULL,
    skip = skip.cell
  )
  feature.names <- read.table(
    file = all.files[['feature list']],
    header = FALSE,
    sep = '\t',
    row.names = NULL,
    skip = skip.feature
  )
  # read barcodes
  bcols <- ncol(x = cell.barcodes)
  if (bcols < cell.column) {
    stop(
      "cell.column was set to ",
      cell.column,
      " but ",
      cells,
      " only has ",
      bcols,
      " columns.",
      " Try setting the cell.column argument to a value <= to ",
      bcols,
      "."
    )
  }
  cell.names <- cell.barcodes[, cell.column]
  if (all(grepl(pattern = "\\-1$", x = cell.names)) & strip.suffix) {
    cell.names <- as.vector(x = as.character(x = sapply(
      X = cell.names,
      FUN = ExtractField,
      field = 1,
      delim = "-"
    )))
  }
  # read features
  fcols <- ncol(x = feature.names)
  if (fcols < feature.column) {
    stop(
      "feature.column was set to ",
      feature.column,
      " but ",
      features,
      " only has ",
      fcols, " column(s).",
      " Try setting the feature.column argument to a value <= to ",
      fcols,
      "."
    )
  }
  if (any(is.na(x = feature.names[, feature.column]))) {
    na.features <- which(x = is.na(x = feature.names[, feature.column]))
    replacement.column <- ifelse(test = feature.column == 2, yes = 1, no = 2)
    if (replacement.column > fcols) {
      stop(
        "Some features names are NA in column ",
        feature.column,
        ". Try specifiying a different column.",
        call. = FALSE
        )
    } else {
      warning(
        "Some features names are NA in column ",
        feature.column,
        ". Replacing NA names with ID from column ",
        replacement.column,
        ".",
        call. = FALSE
        )
    }
    feature.names[na.features, feature.column] <- feature.names[na.features, replacement.column]
  }
  feature.names <- feature.names[, feature.column]
  if (unique.features) {
    feature.names <- make.unique(names = feature.names)
  }
  data <- readMM(file = all.files[['expression matrix']])
  if (length(x = cell.names) != ncol(x = data)) {
    stop(
      "Matrix has ",
      ncol(data),
      " columns but found ", length(cell.names),
      " barcodes. ",
      ifelse(
        test = length(x = cell.names) > ncol(x = data),
        yes = "Try increasing `skip.cell`. ",
        no = ""
      ),
      call. = FALSE
      )
  }
  if (length(x = feature.names) != nrow(x = data)) {
    stop(
      "Matrix has ",
      nrow(data),
      " rows but found ", length(feature.names),
      " features. ",
      ifelse(
        test = length(x = feature.names) > nrow(x = data),
        yes = "Try increasing `skip.feature`. ",
        no = ""
      ),
      call. = FALSE
      )
  }

  colnames(x = data) <- cell.names
  rownames(x = data) <- feature.names
  data <- as.sparse(x = data)
  return(data)
}

#' Read and Load Nanostring SMI data
#'
#' @param data.dir Directory containing all Nanostring SMI files with
#' default filenames
#' @param mtx.file Path to Nanostring cell x gene matrix CSV
#' @param metadata.file Contains metadata including cell center, area,
#' and stain intensities
#' @param molecules.file Path to molecules file
#' @param segmentations.file Path to segmentations CSV
#' @param type Type of cell spatial coordinate matrices to read; choose one
#' or more of:
#' \itemize{
#'  \item \dQuote{centroids}: cell centroids in pixel coordinate space
#'  \item \dQuote{segmentations}: cell segmentations in pixel coordinate space
#' }
#' @param mol.type Type of molecule spatial coordinate matrices to read;
#' choose one or more of:
#' \itemize{
#'  \item \dQuote{pixels}: molecule coordinates in pixel space
#' }
#' @param metadata Type of available metadata to read;
#' choose zero or more of:
#' \itemize{
#'  \item \dQuote{Area}: number of pixels in cell segmentation
#'  \item \dQuote{fov}: cell's fov
#'  \item \dQuote{Mean.MembraneStain}: mean membrane stain intensity
#'  \item \dQuote{Mean.DAPI}: mean DAPI stain intensity
#'  \item \dQuote{Mean.G}: mean green channel stain intensity
#'  \item \dQuote{Mean.Y}: mean yellow channel stain intensity
#'  \item \dQuote{Mean.R}: mean red channel stain intensity
#'  \item \dQuote{Max.MembraneStain}: max membrane stain intensity
#'  \item \dQuote{Max.DAPI}: max DAPI stain intensity
#'  \item \dQuote{Max.G}: max green channel stain intensity
#'  \item \dQuote{Max.Y}: max yellow stain intensity
#'  \item \dQuote{Max.R}: max red stain intensity
#' }
#' @param mols.filter Filter molecules that match provided string
#' @param genes.filter Filter genes from cell x gene matrix that match
#' provided string
#' @param fov.filter Only load in select FOVs. Nanostring SMI data contains
#' 30 total FOVs.
#' @param subset.counts.matrix If the counts matrix should be built from
#' molecule coordinates for a specific segmentation; One of:
#' \itemize{
#'  \item \dQuote{Nuclear}: nuclear segmentations
#'  \item \dQuote{Cytoplasm}: cell cytoplasm segmentations
#'  \item \dQuote{Membrane}: cell membrane segmentations
#' }
#' @param cell.mols.only If TRUE, only load molecules within a cell
#'
#' @return \code{ReadNanostring}: A list with some combination of the
#' following values:
#' \itemize{
#'  \item \dQuote{\code{matrix}}: a
#'  \link[Matrix:dgCMatrix-class]{sparse matrix} with expression data; cells
#'   are columns and features are rows
#'  \item \dQuote{\code{centroids}}: a data frame with cell centroid
#'   coordinates in three columns: \dQuote{x}, \dQuote{y}, and \dQuote{cell}
#'  \item \dQuote{\code{pixels}}: a data frame with molecule pixel coordinates
#'   in three columns: \dQuote{x}, \dQuote{y}, and \dQuote{gene}
#' }
#'
#' @importFrom future.apply future_lapply
#'
#' @export
#'
#' @order 1
#'
#' @concept preprocessing
#'
#' @template section-progressr
#' @template section-future
#'
#' @templateVar pkg data.table
#' @template note-reqdpkg
#'
ReadNanostring <- function(
  data.dir,
  mtx.file = NULL,
  metadata.file = NULL,
  molecules.file = NULL,
  segmentations.file = NULL,
  type = 'centroids',
  mol.type = 'pixels',
  metadata = NULL,
  mols.filter = NA_character_,
  genes.filter = NA_character_,
  fov.filter = NULL,
  subset.counts.matrix = NULL,
  cell.mols.only = TRUE
) {
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Please install 'data.table' for this function")
  }

  # Argument checking
  type <- match.arg(
    arg = type,
    choices = c('centroids', 'segmentations'),
    several.ok = TRUE
  )
  mol.type <- match.arg(
    arg = mol.type,
    choices = c('pixels'),
    several.ok = TRUE
  )
  if (!is.null(metadata)) {
    metadata <- match.arg(
      arg = metadata,
      choices = c(
        "Area", "fov", "Mean.MembraneStain", "Mean.DAPI", "Mean.G",
        "Mean.Y", "Mean.R", "Max.MembraneStain", "Max.DAPI", "Max.G",
        "Max.Y", "Max.R"
      ),
      several.ok = TRUE
    )
  }

  use.dir <- all(vapply(
    X = c(mtx.file, metadata.file, molecules.file),
    FUN = function(x) {
      return(is.null(x = x) || is.na(x = x))
    },
    FUN.VALUE = logical(length = 1L)
  ))

  if (use.dir && !dir.exists(paths = data.dir)) {
    stop("Cannot find Nanostring directory ", data.dir)
  }
  # Identify input files
  files <- c(
    matrix = mtx.file %||% '[_a-zA-Z0-9]*_exprMat_file.csv',
    metadata.file = metadata.file %||% '[_a-zA-Z0-9]*_metadata_file.csv',
    molecules.file = molecules.file %||% '[_a-zA-Z0-9]*_tx_file.csv',
    segmentations.file = segmentations.file %||% '[_a-zA-Z0-9]*-polygons.csv'
  )

  files <- vapply(
    X = files,
    FUN = function(x) {
      x <- as.character(x = x)
      if (isTRUE(x = dirname(path = x) == '.')) {
        fnames <- list.files(
          path = data.dir,
          pattern = x,
          recursive = FALSE,
          full.names = TRUE
        )
        return(sort(x = fnames, decreasing = TRUE)[1L])
      } else {
        return(x)
      }
    },
    FUN.VALUE = character(length = 1L),
    USE.NAMES = TRUE
  )
  files[!file.exists(files)] <- NA_character_

  if (all(is.na(x = files))) {
    stop("Cannot find Nanostring input files in ", data.dir)
  }
  # Checking for loading spatial coordinates
  if (!is.na(x = files[['metadata.file']])) {
    pprecoord <- progressor()
    pprecoord(
      message = "Preloading cell spatial coordinates",
      class = 'sticky',
      amount = 0
    )
    md <- data.table::fread(
      file = files[['metadata.file']],
      sep = ',',
      data.table = FALSE,
      verbose = FALSE
    )

    # filter metadata file by FOVs
    if (!is.null(x = fov.filter)) {
      md <- md[md$fov %in% fov.filter,]
    }
    pprecoord(type = 'finish')
  }
  if (!is.na(x = files[['segmentations.file']])) {
    ppresegs <- progressor()
    ppresegs(
      message = "Preloading cell segmentation vertices",
      class = 'sticky',
      amount = 0
    )
    segs <- data.table::fread(
      file = files[['segmentations.file']],
      sep = ',',
      data.table = FALSE,
      verbose = FALSE
    )

    # filter metadata file by FOVs
    if (!is.null(x = fov.filter)) {
      segs <- segs[segs$fov %in% fov.filter,]
    }
    ppresegs(type = 'finish')
  }
  # Check for loading of molecule coordinates
  if (!is.na(x = files[['molecules.file']])) {
    ppremol <- progressor()
    ppremol(
      message = "Preloading molecule coordinates",
      class = 'sticky',
      amount = 0
    )
    mx <- data.table::fread(
      file = files[['molecules.file']],
      sep = ',',
      verbose = FALSE
    )

    # filter molecules file by FOVs
    if (!is.null(x = fov.filter)) {
      mx <- mx[mx$fov %in% fov.filter,]
    }

    # Molecules outside of a cell have a cell_ID of 0
    if (cell.mols.only) {
      mx <- mx[mx$cell_ID != 0,]
    }

    if (!is.na(x = mols.filter)) {
      ppremol(
        message = paste("Filtering molecules with pattern", mols.filter),
        class = 'sticky',
        amount = 0
      )
      mx <- mx[!grepl(pattern = mols.filter, x = mx$target), , drop = FALSE]
    }
    ppremol(type = 'finish')
    mols <- rep_len(x = files[['molecules.file']], length.out = length(x = mol.type))
    names(x = mols) <- mol.type
    files <- c(files, mols)
    files <- files[setdiff(x = names(x = files), y = 'molecules.file')]
  }
  files <- files[!is.na(x = files)]

  outs <- list("matrix"=NULL, "pixels"=NULL, "centroids"=NULL)
  if (!is.null(metadata)) {
    outs <- append(outs, list("metadata" = NULL))
  }
  if ("segmentations" %in% type) {
    outs <- append(outs, list("segmentations" = NULL))
  }

  for (otype in names(x = outs)) {
    outs[[otype]] <- switch(
      EXPR = otype,
      'matrix' = {
        ptx <- progressor()
        ptx(message = 'Reading counts matrix', class = 'sticky', amount = 0)
        if (!is.null(subset.counts.matrix)) {
          tx <- build.cellcomp.matrix(mols.df=mx, class=subset.counts.matrix)
        } else {
          tx <- data.table::fread(
            file = files[[otype]],
            sep = ',',
            data.table = FALSE,
            verbose = FALSE
          )
          # Combination of Cell ID (for non-zero cell_IDs) and FOV are assumed to be unique. Used to create barcodes / rownames.
          bcs <- paste0(as.character(tx$cell_ID), "_", tx$fov)
          rownames(x = tx) <- bcs
          # remove all rows which represent counts of mols not assigned to a cell for each FOV
          tx <- tx[!tx$cell_ID == 0,]
          # filter fovs from counts matrix
          if (!is.null(x = fov.filter)) {
            tx <- tx[tx$fov %in% fov.filter,]
          }
          tx <- subset(tx, select = -c(fov, cell_ID))
        }

        tx <- as.data.frame(t(x = as.matrix(x = tx[, -1, drop = FALSE])))
        if (!is.na(x = genes.filter)) {
          ptx(
            message = paste("Filtering genes with pattern", genes.filter),
            class = 'sticky',
            amount = 0
          )
          tx <- tx[!grepl(pattern = genes.filter, x = rownames(x = tx)), , drop = FALSE]
        }
        # only keep cells with counts greater than 0
        tx <- tx[, which(colSums(tx) != 0)]
        ratio <- getOption(x = 'Seurat.input.sparse_ratio', default = 0.4)

        if ((sum(tx == 0) / length(x = tx)) > ratio) {
          ptx(
            message = 'Converting counts to sparse matrix',
            class = 'sticky',
            amount = 0
          )
          tx <- as.sparse(x = tx)
        }

        ptx(type = 'finish')

        tx
      },
      'centroids' = {
        pcents <- progressor()
        pcents(
          message = 'Creating centroid coordinates',
          class = 'sticky',
          amount = 0
        )
        pcents(type = 'finish')
        data.frame(
          x = md$CenterX_global_px,
          y = md$CenterY_global_px,
          cell = paste0(as.character(md$cell_ID), "_", md$fov),
          stringsAsFactors = FALSE
        )
      },
      'segmentations' = {
        pcents <- progressor()
        pcents(
          message = 'Creating segmentation coordinates',
          class = 'sticky',
          amount = 0
        )
        pcents(type = 'finish')
        data.frame(
          x = segs$x_global_px,
          y = segs$y_global_px,
          cell = paste0(as.character(segs$cellID), "_", segs$fov),  # cell_ID column in this file doesn't have an underscore
          stringsAsFactors = FALSE
        )
      },
      'metadata' = {
        pmeta <- progressor()
        pmeta(
          message = 'Loading metadata',
          class = 'sticky',
          amount = 0
        )
        pmeta(type = 'finish')
        df <- md[,metadata]
        df$cell <- paste0(as.character(md$cell_ID), "_", md$fov)
        df
      },
      'pixels' = {
        ppixels <- progressor()
        ppixels(
          message = 'Creating pixel-level molecule coordinates',
          class = 'sticky',
          amount = 0
        )
        df <- data.frame(
          x = mx$x_global_px,
          y = mx$y_global_px,
          gene = mx$target,
          stringsAsFactors = FALSE
        )
        ppixels(type = 'finish')
        df
      },
      # 'microns' = {
      #   pmicrons <- progressor()
      #   pmicrons(
      #     message = "Creating micron-level molecule coordinates",
      #     class = 'sticky',
      #     amount = 0
      #   )
      #   df <- data.frame(
      #     x = mx$global_x,
      #     y = mx$global_y,
      #     gene = mx$gene,
      #     stringsAsFactors = FALSE
      #   )
      #   pmicrons(type = 'finish')
      #   df
      # },
      stop("Unknown Nanostring input type: ", outs[[otype]])
    )
  }
  return(outs)
}

#' Read and Load 10x Genomics Xenium in-situ data
#'
#' @param data.dir Directory containing all Xenium output files with
#' default filenames
#' @param outs Types of molecular outputs to read; choose one or more of:
#' \itemize{
#'  \item \dQuote{matrix}: the counts matrix
#'  \item \dQuote{microns}: molecule coordinates
#' }
#' @param type Type of cell spatial coordinate matrices to read; choose one
#' or more of:
#' \itemize{
#'  \item \dQuote{centroids}: cell centroids in pixel coordinate space
#'  \item \dQuote{segmentations}: cell segmentations in pixel coordinate space
#' }
#' @param mols.qv.threshold Remove transcript molecules with
#' a QV less than this threshold. QV >= 20 is the standard threshold
#' used to construct the cell x gene count matrix.
#'
#' @return \code{ReadXenium}: A list with some combination of the
#' following values:
#' \itemize{
#'  \item \dQuote{\code{matrix}}: a
#'  \link[Matrix:dgCMatrix-class]{sparse matrix} with expression data; cells
#'   are columns and features are rows
#'  \item \dQuote{\code{centroids}}: a data frame with cell centroid
#'   coordinates in three columns: \dQuote{x}, \dQuote{y}, and \dQuote{cell}
#'  \item \dQuote{\code{pixels}}: a data frame with molecule pixel coordinates
#'   in three columns: \dQuote{x}, \dQuote{y}, and \dQuote{gene}
#' }
#'
#'
#' @export
#' @concept preprocessing
#'
ReadXenium <- function(
  data.dir,
  outs = c("matrix", "microns"),
  type = "centroids",
  mols.qv.threshold = 20
) {
  # Argument checking
  type <- match.arg(
    arg = type,
    choices = c("centroids", "segmentations"),
    several.ok = TRUE
  )

  outs <- match.arg(
    arg = outs,
    choices = c("matrix", "microns"),
    several.ok = TRUE
  )

  outs <- c(outs, type)

  has_dt <- requireNamespace("data.table", quietly = TRUE) && requireNamespace("R.utils", quietly = TRUE)

  data <- sapply(outs, function(otype) {
    switch(
      EXPR = otype,
      'matrix' = {
        pmtx <- progressor()
        pmtx(message = 'Reading counts matrix', class = 'sticky', amount = 0)
        matrix <- suppressWarnings(Read10X(data.dir = file.path(data.dir, "cell_feature_matrix/")))
        pmtx(type = "finish")
        matrix
      },
      'centroids' = {
        pcents <- progressor()
        pcents(
          message = 'Loading cell centroids',
          class = 'sticky',
          amount = 0
        )
        if (has_dt) {
          cell_info <- as.data.frame(data.table::fread(file.path(data.dir, "cells.csv.gz")))
        } else {
          cell_info <- read.csv(file.path(data.dir, "cells.csv.gz"))
        }
        cell_centroid_df <- data.frame(
          x = cell_info$x_centroid,
          y = cell_info$y_centroid,
          cell = cell_info$cell_id,
          stringsAsFactors = FALSE
        )
        pcents(type = 'finish')
        cell_centroid_df
      },
      'segmentations' = {
        psegs <- progressor()
        psegs(
          message = 'Loading cell segmentations',
          class = 'sticky',
          amount = 0
        )

        # load cell boundaries
        if (has_dt) {
          cell_boundaries_df <- as.data.frame(data.table::fread(file.path(data.dir, "cell_boundaries.csv.gz")))
        } else {
          cell_boundaries_df <- read.csv(file.path(data.dir, "cell_boundaries.csv.gz"), stringsAsFactors = FALSE)
        }
        names(cell_boundaries_df) <- c("cell", "x", "y")
        psegs(type = "finish")
        cell_boundaries_df
      },
      'microns' = {
        pmicrons <- progressor()
        pmicrons(
          message = "Loading molecule coordinates",
          class = 'sticky',
          amount = 0
        )

        # molecules
        if (has_dt) {
          tx_dt <- as.data.frame(data.table::fread(file.path(data.dir, "transcripts.csv.gz")))
          transcripts <- subset(tx_dt, qv >= mols.qv.threshold)
        } else {
          transcripts <- read.csv(file.path(data.dir, "transcripts.csv.gz"))
          transcripts <- subset(transcripts, qv >= mols.qv.threshold)
        }

        df <-
          data.frame(
            x = transcripts$x_location,
            y = transcripts$y_location,
            gene = transcripts$feature_name,
            stringsAsFactors = FALSE
          )
        pmicrons(type = 'finish')
        df
      },
      stop("Unknown Xenium input type: ", otype)
    )
  }, USE.NAMES = TRUE)
  return(data)
}

#' Load Slide-seq spatial data
#'
#' @param coord.file Path to csv file containing bead coordinate positions
#' @param assay Name of assay to associate image to
#'
#' @return A \code{\link{SlideSeq}} object
#'
#' @importFrom utils read.csv
#'
#' @seealso \code{\link{SlideSeq}}
#'
#' @export
#' @concept preprocessing
#'
ReadSlideSeq <- function(coord.file, assay = 'Spatial') {
  if (!file.exists(paths = coord.file)) {
    stop("Cannot find coord file ", coord.file, call. = FALSE)
  }
  slide.seq <- new(
    Class = 'SlideSeq',
    assay = assay,
    coordinates = read.csv(
      file = coord.file,
      header = TRUE,
      as.is = TRUE,
      row.names = 1
    )
  )
  return(slide.seq)
}

#' Read Data From Vitessce
#'
#' Read in data from Vitessce-formatted JSON files
#'
#' @param counts Path or URL to a Vitessce-formatted JSON file with
#' expression data; should end in \dQuote{\code{.genes.json}} or
#' \dQuote{\code{.clusters.json}}; pass \code{NULL} to skip
#' @param coords Path or URL to a Vitessce-formatted JSON file with cell/spot
#' spatial coordinates; should end in \dQuote{\code{.cells.json}};
#' pass \code{NULL} to skip
#' @param molecules Path or URL to a Vitessce-formatted JSON file with molecule
#' spatial coordinates; should end in \dQuote{\code{.molecules.json}};
#' pass \code{NULL} to skip
#' @param type Type of cell/spot spatial coordinates to return,
#' choose one or more from:
#' \itemize{
#'  \item \dQuote{segmentations} cell/spot segmentations
#'  \item \dQuote{centroids} cell/spot centroids
#' }
#' @param filter A character to filter molecules by, pass \code{NA} to skip
#' molecule filtering
#'
#' @return \code{ReadVitessce}: A list with some combination of the
#' following values:
#' \itemize{
#'  \item \dQuote{\code{counts}}: if \code{counts} is not \code{NULL}, an
#'   expression matrix with cells as columns and features as rows
#'  \item \dQuote{\code{centroids}}: if \code{coords} is not \code{NULL} and
#'   \code{type} is contains\dQuote{centroids}, a data frame with cell centroids
#'   in three columns: \dQuote{x}, \dQuote{y}, and \dQuote{cell}
#'  \item \dQuote{\code{segmentations}}: if \code{coords} is not \code{NULL} and
#'   \code{type} contains \dQuote{centroids}, a data frame with cell
#'   segmentations in three columns: \dQuote{x}, \dQuote{y} and \dQuote{cell}
#'  \item \dQuote{\code{molecules}}: if \code{molecules} is not \code{NULL}, a
#'   data frame with molecule spatial coordinates in three columns: \dQuote{x},
#'   \dQuote{y}, and \dQuote{gene}
#' }
#'
#' @importFrom jsonlite read_json
#' @importFrom tools file_ext file_path_sans_ext
#'
#' @export
#'
#' @order 1
#'
#' @concept preprocessing
#'
#' @template section-progressr
#'
#' @templateVar pkg jsonlite
#' @template note-reqdpkg
#'
#' @examples
#' \dontrun{
#' coords <- ReadVitessce(
#'   counts =
#'      "https://s3.amazonaws.com/vitessce-data/0.0.31/master_release/wang/wang.genes.json",
#'   coords =
#'      "https://s3.amazonaws.com/vitessce-data/0.0.31/master_release/wang/wang.cells.json",
#'   molecules =
#'      "https://s3.amazonaws.com/vitessce-data/0.0.31/master_release/wang/wang.molecules.json"
#' )
#' names(coords)
#' coords$counts[1:10, 1:10]
#' head(coords$centroids)
#' head(coords$segmentations)
#' head(coords$molecules)
#' }
#'
ReadVitessce <- function(
  counts = NULL,
  coords = NULL,
  molecules = NULL,
  type = c('segmentations', 'centroids'),
  filter = NA_character_
) {
  if (!requireNamespace('jsonlite', quietly = TRUE)) {
    stop("Please install 'jsonlite' for this function")
  }
  type <- match.arg(arg = type, several.ok = TRUE)
  nouts <- c(
    counts %iff% 'counts',
    coords %iff% type,
    molecules %iff% 'molecules'
  )
  outs <- vector(mode = 'list', length = length(x = nouts))
  names(x = outs) <- nouts
  if (!is.null(x = coords)) {
    ppreload <- progressor()
    ppreload(message = "Preloading coordinates", class = 'sticky', amount = 0)
    cells <- read_json(path = coords)
    ppreload(type = 'finish')
  }
  for (i in nouts) {
    outs[[i]] <- switch(
      EXPR = i,
      'counts' = {
        counts.type <- file_ext(x = basename(path = file_path_sans_ext(
          x = counts
        )))
        cts <- switch(
          EXPR = counts.type,
          'clusters' = .ReadVitessceClusters(counts = counts),
          'genes' = .ReadVitessceGenes(counts = counts),
          stop("Unknown Vitessce counts filetype: '", counts.type, "'")
        )
        pcts <- progressor()
        if (!is.na(x = filter)) {
          pcts(
            message = paste("Filtering genes with pattern", filter),
            class = 'sticky',
            amount = 0
          )
          cts <- cts[!grepl(pattern = filter, x = rownames(x = cts)), , drop = FALSE]
        }
        ratio <- getOption(x = 'Seurat.input.sparse_ratio', default = 0.4)
        if ((sum(cts == 0) / length(x = cts)) > ratio) {
          pcts(
            message = 'Converting counts to sparse matrix',
            class = 'sticky',
            amount = 0
          )
          cts <- as.sparse(x = cts)
        }
        pcts(type = 'finish')
        cts
      },
      'centroids' = {
        pcents <- progressor(steps = length(x = cells))
        pcents(message = "Reading centroids", class = 'sticky', amount = 0)
        centroids <- lapply(
          X = names(x = cells),
          FUN = function(x) {
            cents <- cells[[x]]$xy
            names(x = cents) <- c('x', 'y')
            cents <- as.data.frame(x = cents)
            cents$cell <- x
            pcents()
            return(cents)
          }
        )
        pcents(type = 'finish')
        do.call(what = 'rbind', args = centroids)
      },
      'segmentations' = {
        psegs <- progressor(steps = length(x = cells))
        psegs(message = "Reading segmentations", class = 'sticky', amount = 0)
        segmentations <- lapply(
          X = names(x = cells),
          FUN = function(x) {
            poly <- cells[[x]]$poly
            poly <- lapply(X = poly, FUN = unlist)
            poly <- as.data.frame(x = do.call(what = 'rbind', args = poly))
            colnames(x = poly) <- c('x', 'y')
            poly$cell <- x
            psegs()
            return(poly)
          }
        )
        psegs(type = 'finish')
        do.call(what = 'rbind', args = segmentations)
      },
      'molecules' = {
        pmols1 <- progressor()
        pmols1(message = "Reading molecules", class = 'sticky', amount = 0)
        pmols1(type = 'finish')
        mols <- read_json(path = molecules)
        pmols2 <- progressor(steps = length(x = mols))
        mols <- lapply(
          X = names(x = mols),
          FUN = function(m) {
            x <- mols[[m]]
            x <- lapply(X = x, FUN = unlist)
            x <- as.data.frame(x = do.call(what = 'rbind', args = x))
            colnames(x = x) <- c('x', 'y')
            x$gene <- m
            pmols2()
            return(x)
          }
        )
        mols <- do.call(what = 'rbind', args = mols)
        pmols2(type = 'finish')
        if (!is.na(x = filter)) {
          pmols3 <- progressor()
          pmols3(
            message = paste("Filtering molecules with pattern", filter),
            class = 'sticky',
            amount = 0
          )
          pmols3(type = 'finish')
          mols <- mols[!grepl(pattern = filter, x = mols$gene), , drop = FALSE]
        }
        mols
      },
      stop("Unknown data type: ", i)
    )
  }
  return(outs)
}

#' Read and Load MERFISH Input from Vizgen
#'
#' Read and load in MERFISH data from Vizgen-formatted files
#'
#' @inheritParams ReadVitessce
#' @param data.dir Path to the directory with Vizgen MERFISH files; requires at
#' least one of the following files present:
#' \itemize{
#'  \item \dQuote{\code{cell_by_gene.csv}}: used for reading count matrix
#'  \item \dQuote{\code{cell_metadata.csv}}: used for reading cell spatial
#'  coordinate matrices
#'  \item \dQuote{\code{detected_transcripts.csv}}: used for reading molecule
#'  spatial coordinate matrices
#' }
#' @param transcripts Optional file path for counts matrix; pass \code{NA} to
#' suppress reading counts matrix
#' @param spatial Optional file path for spatial metadata; pass \code{NA} to
#' suppress reading spatial coordinates. If \code{spatial} is provided and
#' \code{type} is \dQuote{segmentations}, uses \code{dirname(spatial)} instead of
#' \code{data.dir} to find HDF5 files
#' @param molecules Optional file path for molecule coordinates file; pass
#' \code{NA} to suppress reading spatial molecule information
#' @param type Type of cell spatial coordinate matrices to read; choose one
#' or more of:
#' \itemize{
#'  \item \dQuote{segmentations}: cell segmentation vertices; requires
#'  \href{https://cran.r-project.org/package=hdf5r}{\pkg{hdf5r}} to be
#'   installed and requires a directory \dQuote{\code{cell_boundaries}} within
#'   \code{data.dir}. Within \dQuote{\code{cell_boundaries}}, there must be
#'   one or more HDF5 file named \dQuote{\code{feature_data_##.hdf5}}
#'  \item \dQuote{centroids}: cell centroids in micron coordinate space
#'  \item \dQuote{boxes}: cell box outlines in micron coordinate space
#' }
#' @param mol.type Type of molecule spatial coordinate matrices to read;
#' choose one or more of:
#' \itemize{
#'  \item \dQuote{pixels}: molecule coordinates in pixel space
#'  \item \dQuote{microns}: molecule coordinates in micron space
#' }
#' @param metadata Type of available metadata to read;
#' choose zero or more of:
#' \itemize{
#'  \item \dQuote{volume}: estimated cell volume
#'  \item \dQuote{fov}: cell's fov
#' }
#' @param z Z-index to load; must be between 0 and 6, inclusive
#'
#' @return \code{ReadVizgen}: A list with some combination of the
#' following values:
#' \itemize{
#'  \item \dQuote{\code{transcripts}}: a
#'  \link[Matrix:dgCMatrix-class]{sparse matrix} with expression data; cells
#'   are columns and features are rows
#'  \item \dQuote{\code{segmentations}}: a data frame with cell polygon outlines in
#'   three columns: \dQuote{x}, \dQuote{y}, and \dQuote{cell}
#'  \item \dQuote{\code{centroids}}: a data frame with cell centroid
#'   coordinates in three columns: \dQuote{x}, \dQuote{y}, and \dQuote{cell}
#'  \item \dQuote{\code{boxes}}: a data frame with cell box outlines in three
#'   columns: \dQuote{x}, \dQuote{y}, and \dQuote{cell}
#'  \item \dQuote{\code{microns}}: a data frame with molecule micron
#'   coordinates in three columns: \dQuote{x}, \dQuote{y}, and \dQuote{gene}
#'  \item \dQuote{\code{pixels}}: a data frame with molecule pixel coordinates
#'   in three columns: \dQuote{x}, \dQuote{y}, and \dQuote{gene}
#'  \item \dQuote{\code{metadata}}: a data frame with the cell-level metadata
#'   requested by \code{metadata}
#' }
#'
#' @importFrom future.apply future_lapply
#'
#' @export
#'
#' @order 1
#'
#' @concept preprocessing
#'
#' @template section-progressr
#' @template section-future
#'
#' @templateVar pkg data.table
#' @template note-reqdpkg
#'
ReadVizgen <- function(
  data.dir,
  transcripts = NULL,
  spatial = NULL,
  molecules = NULL,
  type = 'segmentations',
  mol.type = 'microns',
  metadata = NULL,
  filter = NA_character_,
  z = 3L
) {
  # TODO: handle multiple segmentations per z-plane
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Please install 'data.table' for this function")
  }
  # hdf5r is only used for loading polygon boundaries
  # Not needed for all Vizgen input
  hdf5 <- requireNamespace("hdf5r", quietly = TRUE)
  # Argument checking
  type <- match.arg(
    arg = type,
    choices = c('segmentations', 'centroids', 'boxes'),
    several.ok = TRUE
  )
  mol.type <- match.arg(
    arg = mol.type,
    choices = c('pixels', 'microns'),
    several.ok = TRUE
  )
  if (!is.null(x = metadata)) {
    metadata <- match.arg(
      arg = metadata,
      choices = c("volume", "fov"),
      several.ok = TRUE
    )
  }
  if (!z %in% seq.int(from = 0L, to = 6L)) {
    stop("The z-index must be in the range [0, 6]")
  }
  use.dir <- all(vapply(
    X = c(transcripts, spatial, molecules),
    FUN = function(x) {
      return(is.null(x = x) || is.na(x = x))
    },
    FUN.VALUE = logical(length = 1L)
  ))
  if (use.dir && !dir.exists(paths = data.dir)) {
    stop("Cannot find Vizgen directory ", data.dir)
  }
  # Identify input files
  files <- c(
    transcripts = transcripts %||% 'cell_by_gene[_a-zA-Z0-9]*.csv',
    spatial = spatial %||% 'cell_metadata[_a-zA-Z0-9]*.csv',
    molecules = molecules %||% 'detected_transcripts[_a-zA-Z0-9]*.csv'
  )
  files[is.na(x = files)] <- NA_character_
  h5dir <- file.path(
    ifelse(
      test = dirname(path = files['spatial']) == '.',
      yes = data.dir,
      no = dirname(path = files['spatial'])
    ),
    'cell_boundaries'
  )
  zidx <- paste0('zIndex_', z)
  files <- vapply(
    X = files,
    FUN = function(x) {
      x <- as.character(x = x)
      if (isTRUE(x = dirname(path = x) == '.')) {
        fnames <- list.files(
          path = data.dir,
          pattern = x,
          recursive = FALSE,
          full.names = TRUE
        )
        return(sort(x = fnames, decreasing = TRUE)[1L])
      } else {
        return(x)
      }
    },
    FUN.VALUE = character(length = 1L),
    USE.NAMES = TRUE
  )
  files[!file.exists(files)] <- NA_character_
  if (all(is.na(x = files))) {
    stop("Cannot find Vizgen input files in ", data.dir)
  }
  # Checking for loading spatial coordinates
  if (!is.na(x = files[['spatial']])) {
    pprecoord <- progressor()
    pprecoord(
      message = "Preloading cell spatial coordinates",
      class = 'sticky',
      amount = 0
    )
    sp <- data.table::fread(
      file = files[['spatial']],
      sep = ',',
      data.table = FALSE,
      verbose = FALSE
      # showProgress = progressr:::progressr_in_globalenv(action = 'query')
      # showProgress = verbose
    )
    pprecoord(type = 'finish')
    rownames(x = sp) <- as.character(x = sp[, 1])
    sp <- sp[, -1, drop = FALSE]
    # Check to see if we should load segmentations
    if ('segmentations' %in% type) {
      poly <- if (isFALSE(x = hdf5)) {
        warning(
          "Cannot find hdf5r; unable to load segmentation vertices",
          immediate. = TRUE
        )
        FALSE
      } else if (!dir.exists(paths = h5dir)) {
        warning("Cannot find cell boundary H5 files", immediate. = TRUE)
        FALSE
      } else {
        TRUE
      }
      if (isFALSE(x = poly)) {
        type <- setdiff(x = type, y = 'segmentations')
      }
    }
    spatials <- rep_len(x = files[['spatial']], length.out = length(x = type))
    names(x = spatials) <- type
    files <- c(files, spatials)
    files <- files[setdiff(x = names(x = files), y = 'spatial')]
  } else if (!is.null(x = metadata)) {
    warning(
      "metadata can only be loaded when spatial coordinates are loaded",
      immediate. = TRUE
    )
    metadata <- NULL
  }
  # Check for loading of molecule coordinates
  if (!is.na(x = files[['molecules']])) {
    ppremol <- progressor()
    ppremol(
      message = "Preloading molecule coordinates",
      class = 'sticky',
      amount = 0
    )
    mx <- data.table::fread(
      file = files[['molecules']],
      sep = ',',
      verbose = FALSE
      # showProgress = verbose
    )
    mx <- mx[mx$global_z == z, , drop = FALSE]
    if (!is.na(x = filter)) {
      ppremol(
        message = paste("Filtering molecules with pattern", filter),
        class = 'sticky',
        amount = 0
      )
      mx <- mx[!grepl(pattern = filter, x = mx$gene), , drop = FALSE]
    }
    ppremol(type = 'finish')
    mols <- rep_len(x = files[['molecules']], length.out = length(x = mol.type))
    names(x = mols) <- mol.type
    files <- c(files, mols)
    files <- files[setdiff(x = names(x = files), y = 'molecules')]
  }
  files <- files[!is.na(x = files)]
  # Read input data
  outs <- vector(mode = 'list', length = length(x = files))
  names(x = outs) <- names(x = files)
  if (!is.null(metadata)) {
    outs <- c(outs, list(metadata = NULL))
  }
  for (otype in names(x = outs)) {
    outs[[otype]] <- switch(
      EXPR = otype,
      'transcripts' = {
        ptx <- progressor()
        ptx(message = 'Reading counts matrix', class = 'sticky', amount = 0)
        tx <- data.table::fread(
          file = files[[otype]],
          sep = ',',
          data.table = FALSE,
          verbose = FALSE
        )
        rownames(x = tx) <- as.character(x = tx[, 1])
        tx <- t(x = as.matrix(x = tx[, -1, drop = FALSE]))
        if (!is.na(x = filter)) {
          ptx(
            message = paste("Filtering genes with pattern", filter),
            class = 'sticky',
            amount = 0
          )
          tx <- tx[!grepl(pattern = filter, x = rownames(x = tx)), , drop = FALSE]
        }
        ratio <- getOption(x = 'Seurat.input.sparse_ratio', default = 0.4)
        if ((sum(tx == 0) / length(x = tx)) > ratio) {
          ptx(
            message = 'Converting counts to sparse matrix',
            class = 'sticky',
            amount = 0
          )
          tx <- as.sparse(x = tx)
        }
        ptx(type = 'finish')
        tx
      },
      'centroids' = {
        pcents <- progressor()
        pcents(
          message = 'Creating centroid coordinates',
          class = 'sticky',
          amount = 0
        )
        pcents(type = 'finish')
        data.frame(
          x = sp$center_x,
          y = sp$center_y,
          cell = rownames(x = sp),
          stringsAsFactors = FALSE
        )
      },
      'segmentations' = {
        ppoly <- progressor(steps = length(x = unique(x = sp$fov)))
        ppoly(
          message = "Creating polygon coordinates",
          class = 'sticky',
          amount = 0
        )
        pg <- future_lapply(
          X = unique(x = sp$fov),
          FUN = function(f, ...) {
            fname <- file.path(h5dir, paste0('feature_data_', f, '.hdf5'))
            if (!file.exists(fname)) {
              warning(
                "Cannot find HDF5 file for field of view ",
                f,
                immediate. = TRUE
              )
              return(NULL)
            }
            hfile <- hdf5r::H5File$new(filename = fname, mode = 'r')
            on.exit(expr = hfile$close_all())
            cells <- rownames(x = subset(x = sp, subset = fov == f))
            df <- lapply(
              X = cells,
              FUN = function(x) {
                return(tryCatch(
                  expr = {
                    cc <- hfile[['featuredata']][[x]][[zidx]][['p_0']][['coordinates']]$read()
                    cc <- as.data.frame(x = t(x = cc))
                    colnames(x = cc) <- c('x', 'y')
                    cc$cell <- x
                    cc
                  },
                  error = function(...) {
                    return(NULL)
                  }
                ))
              }
            )
            ppoly()
            return(do.call(what = 'rbind', args = df))
          }
        )
        ppoly(type = 'finish')
        pg <- do.call(what = 'rbind', args = pg)
        npg <- length(x = unique(x = pg$cell))
        if (npg < nrow(x = sp)) {
          warning(
            nrow(x = sp) - npg,
            " cells missing polygon information",
            immediate. = TRUE
          )
        }
        pg
      },
      'boxes' = {
        pbox <- progressor(steps = nrow(x = sp))
        pbox(message = "Creating box coordinates", class = 'sticky', amount = 0)
        bx <- future_lapply(
          X = rownames(x = sp),
          FUN = function(cell) {
            row <- sp[cell, ]
            df <- expand.grid(
              x = c(row$min_x, row$max_x),
              y = c(row$min_y, row$max_y),
              cell = cell,
              KEEP.OUT.ATTRS = FALSE,
              stringsAsFactors = FALSE
            )
            df <- df[c(1, 3, 4, 2), , drop = FALSE]
            pbox()
            return(df)
          }
        )
        pbox(type = 'finish')
        do.call(what = 'rbind', args = bx)
      },
      'metadata' = {
        pmeta <- progressor()
        pmeta(
          message = 'Loading metadata',
          class = 'sticky',
          amount = 0
        )
        pmeta(type = 'finish')
        sp[, metadata, drop = FALSE]
      },
      'pixels' = {
        ppixels <- progressor()
        ppixels(
          message = 'Creating pixel-level molecule coordinates',
          class = 'sticky',
          amount = 0
        )
        df <- data.frame(
          x = mx$x,
          y = mx$y,
          gene = mx$gene,
          stringsAsFactors = FALSE
        )
        # if (!is.na(x = filter)) {
        #   ppixels(
        #     message = paste("Filtering molecules with pattern", filter),
        #     class = 'sticky',
        #     amount = 0
        #   )
        #   df <- df[!grepl(pattern = filter, x = df$gene), , drop = FALSE]
        # }
        ppixels(type = 'finish')
        df
      },
      'microns' = {
        pmicrons <- progressor()
        pmicrons(
          message = "Creating micron-level molecule coordinates",
          class = 'sticky',
          amount = 0
        )
        df <- data.frame(
          x = mx$global_x,
          y = mx$global_y,
          gene = mx$gene,
          stringsAsFactors = FALSE
        )
        # if (!is.na(x = filter)) {
        #   pmicrons(
        #     message = paste("Filtering molecules with pattern", filter),
        #     class = 'sticky',
        #     amount = 0
        #   )
        #   df <- df[!grepl(pattern = filter, x = df$gene), , drop = FALSE]
        # }
        pmicrons(type = 'finish')
        df
      },
      stop("Unknown MERFISH input type: ", type)
    )
  }
  return(outs)
}

#' Normalize raw data to fractions
#'
#' Normalize count data to relative counts per cell by dividing by the total
#' per cell. Optionally use a scale factor, e.g. for counts per million (CPM)
#' use \code{scale.factor = 1e6}.
#'
#' @param data Matrix with the raw count data
#' @param scale.factor Scale the result. Default is 1
#' @param verbose Print progress
#' @return Returns a matrix with the relative counts
#'
#' @importFrom methods as
#' @importFrom Matrix colSums
#'
#' @export
#' @concept preprocessing
#'
#' @examples
#' mat <- matrix(data = rbinom(n = 25, size = 5, prob = 0.2), nrow = 5)
#' mat
#' mat_norm <- RelativeCounts(data = mat)
#' mat_norm
#'
RelativeCounts <- function(data, scale.factor = 1, verbose = TRUE) {
  if (is.data.frame(x = data)) {
    data <- as.matrix(x = data)
  }
  if (!inherits(x = data, what = 'dgCMatrix')) {
    data <- as.sparse(x = data)
  }
  if (verbose) {
    cat("Performing relative-counts-normalization\n", file = stderr())
  }
  norm.data <- data
  norm.data@x <- norm.data@x / rep.int(Matrix::colSums(norm.data), diff(norm.data@p)) * scale.factor
  return(norm.data)
}

#' Run the mark variogram computation on a given position matrix and expression
#' matrix.
#'
#' Wraps the functionality of markvario from the spatstat package.
#'
#' @param spatial.location A 2 column matrix giving the spatial locations of
#' each of the data points also in data
#' @param data Matrix containing the data used as "marks" (e.g. gene expression)
#' @param ... Arguments passed to markvario
#'
#' @importFrom spatstat.explore markvario
#' @importFrom spatstat.geom ppp
#'
#' @export
#' @concept preprocessing
#'
RunMarkVario <- function(
  spatial.location,
  data,
  ...
) {
  pp <- ppp(
    x = spatial.location[, 1],
    y = spatial.location[, 2],
    xrange = range(spatial.location[, 1]),
    yrange = range(spatial.location[, 2])
  )
  if (nbrOfWorkers() > 1) {
    chunks <- nbrOfWorkers()
    features <- rownames(x = data)
    features <- split(
      x = features,
      f = ceiling(x = seq_along(along.with = features) / (length(x = features) / chunks))
    )
    mv <- future_lapply(X = features, FUN = function(x) {
      pp[["marks"]] <- as.data.frame(x = t(x = data[x, ]))
      markvario(X = pp, normalise = TRUE, ...)
    })
    mv <- unlist(x = mv, recursive = FALSE)
    names(x = mv) <- rownames(x = data)
  } else {
    pp[["marks"]] <- as.data.frame(x = t(x = data))
    mv <- markvario(X = pp, normalise = TRUE, ...)
  }
  return(mv)
}

#' Compute Moran's I value.
#'
#' Wraps the functionality of the Moran.I function from the ape package.
#' Weights are computed as 1/distance.
#'
#' @param data Expression matrix
#' @param pos Position matrix
#' @param verbose Display messages/progress
#'
#' @importFrom stats dist
#'
#' @export
#' @concept preprocessing
#'
RunMoransI <- function(data, pos, verbose = TRUE) {
  mysapply <- sapply
  if (verbose) {
    message("Computing Moran's I")
    mysapply <- pbsapply
  }
  Rfast2.installed <- PackageCheck("Rfast2", error = FALSE)
  if (Rfast2.installed) {
    MyMoran <- Rfast2::moranI
  } else if (!PackageCheck('ape', error = FALSE)) {
    stop(
      "'RunMoransI' requires either Rfast2 or ape to be installed",
      call. = FALSE
    )
  } else {
    MyMoran <- ape::Moran.I
    if (getOption('Seurat.Rfast2.msg', TRUE)) {
      message(
        "For a more efficient implementation of the Morans I calculation,",
        "\n(selection.method = 'moransi') please install the Rfast2 package",
        "\n--------------------------------------------",
        "\ninstall.packages('Rfast2')",
        "\n--------------------------------------------",
        "\nAfter installation of Rfast2, Seurat will automatically use the more ",
        "\nefficient implementation (no further action necessary).",
        "\nThis message will be shown once per session"
      )
      options(Seurat.Rfast2.msg = FALSE)
    }
  }
  pos.dist <- dist(x = pos)
  pos.dist.mat <- as.matrix(x = pos.dist)
  # weights as 1/dist^2
  weights <- 1/pos.dist.mat^2
  diag(x = weights) <- 0
  results <- mysapply(X = 1:nrow(x = data), FUN = function(x) {
    tryCatch(
      expr = MyMoran(data[x, ], weights),
      error = function(x) c(1,1,1,1)
    )
  })
  pcol <- ifelse(test = Rfast2.installed, yes = 2, no = 4)
  results <- data.frame(
    observed = unlist(x = results[1, ]),
    p.value = unlist(x = results[pcol, ])
  )
  rownames(x = results) <- rownames(x = data)
  return(results)
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
#' @importFrom methods as
#'
#' @return Matrix with downsampled data
#'
#' @export
#' @concept preprocessing
#'
#' @examples
#' data("pbmc_small")
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
  data <- as.sparse(x = data)
  if (length(x = max.umi) == 1) {
    new_data <- RunUMISampling(
      data = data,
      sample_val = max.umi,
      upsample = upsample,
      display_progress = verbose
    )
  } else if (length(x = max.umi) != ncol(x = data)) {
    stop("max.umi vector not equal to number of cells")
  } else {
    new_data <- RunUMISamplingPerCell(
      data = data,
      sample_val = max.umi,
      upsample = upsample,
      display_progress = verbose
    )
  }
  dimnames(x = new_data) <- dimnames(x = data)
  return(new_data)
}

#' Use regularized negative binomial regression to normalize UMI count data
#'
#' This function calls sctransform::vst. The sctransform package is available at
#' https://github.com/satijalab/sctransform.
#' Use this function as an alternative to the NormalizeData,
#' FindVariableFeatures, ScaleData workflow. Results are saved in a new assay
#' (named SCT by default) with counts being (corrected) counts, data being log1p(counts),
#' scale.data being pearson residuals; sctransform::vst intermediate results are saved
#' in misc slot of new assay.
#'
#' @param object UMI counts matrix
#' @param cell.attr A metadata with cell attributes
#' @param reference.SCT.model If not NULL, compute residuals for the object
#' using the provided SCT model; supports only log_umi as the latent variable.
#' If residual.features are not specified, compute for the top variable.features.n
#' specified in the model which are also present in the object. If
#' residual.features are specified, the variable features of the resulting SCT
#' assay are set to the top variable.features.n in the model.
#' @param do.correct.umi Place corrected UMI matrix in assay counts slot; default is TRUE
#' @param ncells Number of subsampling cells used to build NB regression; default is 5000
#' @param residual.features Genes to calculate residual features for; default is NULL (all genes).
#' If specified, will be set to VariableFeatures of the returned object.
#' @param variable.features.n Use this many features as variable features after
#' ranking by residual variance; default is 3000. Only applied if residual.features is not set.
#' @param variable.features.rv.th Instead of setting a fixed number of variable features,
#' use this residual variance cutoff; this is only used when \code{variable.features.n}
#' is set to NULL; default is 1.3. Only applied if residual.features is not set.
#' @param vars.to.regress Variables to regress out in a second non-regularized linear
#' regression. For example, percent.mito. Default is NULL
#' @param do.scale Whether to scale residuals to have unit variance; default is FALSE
#' @param do.center Whether to center residuals to have mean zero; default is TRUE
#' @param clip.range Range to clip the residuals to; default is \code{c(-sqrt(n/30), sqrt(n/30))},
#' where n is the number of cells
#' @param vst.flavor When set to 'v2' sets method = glmGamPoi_offset, n_cells=2000,
#' and exclude_poisson = TRUE which causes the model to learn theta and intercept
#' only besides excluding poisson genes from learning and regularization
#' @param conserve.memory If set to TRUE the residual matrix for all genes is never
#' created in full; useful for large data sets, but will take longer to run;
#' this will also set return.only.var.genes to TRUE; default is FALSE
#' @param return.only.var.genes If set to TRUE the scale.data matrices in output assay are
#' subset to contain only the variable genes; default is TRUE
#' @param seed.use Set a random seed. By default, sets the seed to 1448145. Setting
#' NULL will not set a seed.
#' @param verbose Whether to print messages and progress bars
#' @param ... Additional parameters passed to \code{sctransform::vst}
#'
#' @return Returns a Seurat object with a new assay (named SCT by default) with
#' counts being (corrected) counts, data being log1p(counts), scale.data being
#' pearson residuals; sctransform::vst intermediate results are saved in misc
#' slot of the new assay.
#'
#' @importFrom stats setNames
#' @importFrom Matrix colSums
#' @importFrom sctransform vst get_residual_var get_residuals correct_counts
#'
#' @seealso \code{\link[sctransform]{correct_counts}} \code{\link[sctransform]{get_residuals}}
#'
#' @rdname SCTransform
#' @concept preprocessing
#' @export
#'
SCTransform.default <- function(
  object,
  cell.attr,
  reference.SCT.model = NULL,
  do.correct.umi = TRUE,
  ncells = 5000,
  residual.features = NULL,
  variable.features.n = 3000,
  variable.features.rv.th = 1.3,
  vars.to.regress = NULL,
  do.scale = FALSE,
  do.center = TRUE,
  clip.range = c(-sqrt(x = ncol(x = umi) / 30), sqrt(x = ncol(x = umi) / 30)),
  vst.flavor = 'v2',
  conserve.memory = FALSE,
  return.only.var.genes = TRUE,
  seed.use = 1448145,
  verbose = TRUE,
  ...
) {
  vst.args <- list(...)
  umi <- object
  # check for batch_var in meta data
  if ('batch_var' %in% names(x = vst.args)) {
    if (!(vst.args[['batch_var']] %in% colnames(x = cell.attr))) {
      stop('batch_var not found in seurat object meta data')
    }
  }
  # parameter checking when reference.SCT.model is set
  if (!is.null(x = reference.SCT.model) ) {
    if (inherits(x = reference.SCT.model, what = "SCTModel")) {
      reference.SCT.model <- SCTModel_to_vst(SCTModel = reference.SCT.model)
    }
    if (is.list(x = reference.SCT.model) & inherits(x = reference.SCT.model[[1]], what = "SCTModel")) {
      stop("reference.SCT.model must be one SCTModel rather than a list of SCTModel")
    }
    if ('latent_var' %in% names(x = vst.args)) {
      stop('custom latent variables are not supported when reference.SCT.model is given')
    }
    if (reference.SCT.model$model_str != 'y ~ log_umi') {
      stop('reference.SCT.model must be derived using default SCT regression formula, `y ~ log_umi`')
    }

  }
  # check for latent_var in meta data
  if ('latent_var' %in% names(x = vst.args)) {
    known.attr <- c('umi', 'gene', 'log_umi', 'log_gene', 'umi_per_gene', 'log_umi_per_gene')
    if (!all(vst.args[['latent_var']] %in% c(colnames(x = cell.attr), known.attr))) {
      stop('latent_var values are not from the set of cell attributes sctransform calculates by default and cannot be found in seurat object meta data')
    }
  }
  # check for vars.to.regress in meta data
  if (any(!vars.to.regress %in% colnames(x = cell.attr))) {
    stop('problem with second non-regularized linear regression; not all variables found in seurat object meta data; check vars.to.regress parameter')
  }
  if (any(c('cell_attr', 'verbosity', 'return_cell_attr', 'return_gene_attr', 'return_corrected_umi') %in% names(x = vst.args))) {
    warning(
      'the following arguments will be ignored because they are set within this function:',
      paste(
        c(
          'cell_attr',
          'verbosity',
          'return_cell_attr',
          'return_gene_attr',
          'return_corrected_umi'
        ),
        collapse = ', '
      ),
      call. = FALSE,
      immediate. = TRUE
    )
  }
  
  vst.args[['vst.flavor']] <- vst.flavor
  vst.args[['umi']] <- umi
  vst.args[['cell_attr']] <- cell.attr
  vst.args[['verbosity']] <- as.numeric(x = verbose) * 2
  vst.args[['return_cell_attr']] <- TRUE
  vst.args[['return_gene_attr']] <- TRUE
  vst.args[['return_corrected_umi']] <- do.correct.umi
  vst.args[['n_cells']] <- min(ncells, ncol(x = umi))
  residual.type <- vst.args[['residual_type']] %||% 'pearson'
  res.clip.range <- vst.args[['res_clip_range']] %||% c(-sqrt(x = ncol(x = umi)), sqrt(x = ncol(x = umi)))
 # set sct normalization method
  if (!is.null( reference.SCT.model)) {
    sct.method <- "reference.model"
  } else if (!is.null(x = residual.features)) {
    sct.method <- "residual.features"
  } else if (conserve.memory) {
    sct.method <- "conserve.memory"
  } else {
    sct.method <- "default"
  }
  # set vst model
  vst.out <- switch(
    EXPR = sct.method,
    'reference.model' = {
      if (verbose) {
        message("Using reference SCTModel to calculate pearson residuals")
      }
      do.center <- FALSE
      do.correct.umi <- FALSE
      vst.out <- reference.SCT.model
      clip.range <- vst.out$arguments$sct.clip.range
      cell_attr <-  data.frame(log_umi = log10(x = colSums(umi)))
      rownames(cell_attr) <- colnames(x = umi)
      vst.out$cell_attr <- cell_attr

      all.features  <- intersect(
        x =  rownames(x = vst.out$gene_attr),
        y = rownames(x = umi)
      )
      vst.out$gene_attr <- vst.out$gene_attr[all.features ,]
      vst.out$model_pars_fit <- vst.out$model_pars_fit[all.features,]
      vst.out
    },
    'residual.features' = {
      if (verbose) {
        message("Computing residuals for the ", length(x = residual.features), " specified features")
      }
      return.only.var.genes <- TRUE
      do.correct.umi <- FALSE
      vst.args[['return_corrected_umi']] <- FALSE
      vst.args[['residual_type']] <- 'none'
      vst.out <- do.call(what = 'vst', args = vst.args)
      vst.out$gene_attr$residual_variance <- NA_real_
      vst.out
    },
    'conserve.memory' = {
      return.only.var.genes <- TRUE
      vst.args[['residual_type']] <- 'none'
      vst.out <- do.call(what = 'vst', args = vst.args)
      feature.variance <- get_residual_var(
        vst_out = vst.out,
        umi = umi,
        residual_type = residual.type,
        res_clip_range = res.clip.range
      )
      vst.out$gene_attr$residual_variance <- NA_real_
      vst.out$gene_attr[names(x = feature.variance), 'residual_variance'] <- feature.variance
      vst.out
    },
    'default' = {
      vst.out <- do.call(what = 'vst', args = vst.args)
      vst.out
    })

  feature.variance <- vst.out$gene_attr[,"residual_variance"]
  names(x = feature.variance) <- rownames(x = vst.out$gene_attr)
  if (verbose) {
    message('Determine variable features')
  }
  feature.variance <- sort(x = feature.variance, decreasing = TRUE)
  if (!is.null(x = variable.features.n)) {
    top.features <- names(x = feature.variance)[1:min(variable.features.n, length(x = feature.variance))]
  } else {
    top.features <- names(x = feature.variance)[feature.variance >= variable.features.rv.th]
  }

  # get residuals
  vst.out <- switch(
    EXPR = sct.method,
    'reference.model' = {
      if (is.null(x = residual.features)) {
        residual.features <- top.features
      }
      residual.features <- Reduce(
        f = intersect,
        x = list(residual.features, rownames(x = umi), rownames(x = vst.out$model_pars_fit))
      )
      residual.feature.mat <- get_residuals(
        vst_out = vst.out,
        umi = umi[residual.features, , drop = FALSE],
        verbosity = as.numeric(x = verbose)*2
      )
      vst.out$gene_attr <- vst.out$gene_attr[residual.features ,]
      ref.residuals.mean <- vst.out$gene_attr[,"residual_mean"]
      vst.out$y <- sweep(
        x = residual.feature.mat,
        MARGIN = 1,
        STATS = ref.residuals.mean,
        FUN = "-"
      )
      vst.out
    },
    'residual.features' = {
      residual.features <- intersect(
        x = residual.features,
        y = rownames(x = vst.out$gene_attr)
      )
      residual.feature.mat <- get_residuals(
        vst_out = vst.out,
        umi = umi[residual.features, , drop = FALSE],
        verbosity = as.numeric(x = verbose)*2
      )
      vst.out$y <- residual.feature.mat
      vst.out$gene_attr$residual_mean <- NA_real_
      vst.out$gene_attr$residual_variance <- NA_real_
      vst.out$gene_attr[residual.features, "residual_mean"] <- rowMeans2(x = vst.out$y)
      vst.out$gene_attr[residual.features, "residual_variance"] <- RowVar(x = vst.out$y)
      vst.out
    },
    'conserve.memory' = {
      vst.out$y <- get_residuals(
        vst_out = vst.out,
        umi = umi[top.features, ],
        residual_type = residual.type,
        res_clip_range = res.clip.range,
        verbosity = as.numeric(x = verbose)*2
      )
      vst.out$gene_attr$residual_mean <- NA_real_
      vst.out$gene_attr[top.features, "residual_mean"] = rowMeans2(x =  vst.out$y)
      if (do.correct.umi & residual.type == 'pearson') {
        vst.out$umi_corrected <- correct_counts(
          x = vst.out,
          umi = umi,
          verbosity = as.numeric(x = verbose) * 2
        )
      }
      vst.out
    },
    'default' = {
      if (return.only.var.genes) {
        vst.out$y <- vst.out$y[top.features, ]
      }
      vst.out
    })

  scale.data <- vst.out$y
  # clip the residuals
  scale.data[scale.data < clip.range[1]] <- clip.range[1]
  scale.data[scale.data > clip.range[2]] <- clip.range[2]
  # 2nd regression
  scale.data <- ScaleData(
    scale.data,
    features = NULL,
    vars.to.regress = vars.to.regress,
    latent.data = cell.attr[, vars.to.regress, drop = FALSE],
    model.use = 'linear',
    use.umi = FALSE,
    do.scale = do.scale,
    do.center = do.center,
    scale.max = Inf,
    block.size = 750,
    min.cells.to.block = 3000,
    verbose = verbose
  )
  vst.out$y <- scale.data
  vst.out$variable_features <- residual.features %||% top.features

  return(vst.out)
}

#' @rdname SCTransform
#' @concept preprocessing
#' @export
#' @method SCTransform Assay
#'
SCTransform.Assay <- function(
    object,
    cell.attr,
    reference.SCT.model = NULL,
    do.correct.umi = TRUE,
    ncells = 5000,
    residual.features = NULL,
    variable.features.n = 3000,
    variable.features.rv.th = 1.3,
    vars.to.regress = NULL,
    do.scale = FALSE,
    do.center = TRUE,
    clip.range = c(-sqrt(x = ncol(x = object) / 30), sqrt(x = ncol(x = object) / 30)),
    vst.flavor = 'v2',
    conserve.memory = FALSE,
    return.only.var.genes = TRUE,
    seed.use = 1448145,
    verbose = TRUE,
    ...
) {
  if (!is.null(x = seed.use)) {
    set.seed(seed = seed.use)
  }
  if (!is.null(reference.SCT.model)){
    do.correct.umi <- FALSE
    do.center <- FALSE
  }
  umi <- GetAssayData(object = object, slot = 'counts')
  vst.out <- SCTransform(object = umi,
                         cell.attr = cell.attr,
                         reference.SCT.model = reference.SCT.model,
                         do.correct.umi = do.correct.umi,
                         ncells = ncells,
                         residual.features = residual.features,
                         variable.features.n = variable.features.n,
                         variable.features.rv.th = variable.features.rv.th,
                         vars.to.regress = vars.to.regress,
                         do.scale = do.scale,
                         do.center = do.center,
                         clip.range = clip.range,
                         vst.flavor = vst.flavor,
                         conserve.memory = conserve.memory,
                         return.only.var.genes = return.only.var.genes,
                         seed.use = seed.use,
                         verbose = verbose,
                         ...)
  residual.type <- vst.out[['residual_type']] %||% 'pearson'
  sct.method <- vst.out[["sct.method"]]
  # create output assay and put (corrected) umi counts in count slot
  if (do.correct.umi & residual.type == 'pearson') {
    if (verbose) {
      message('Place corrected count matrix in counts slot')
    }
    assay.out <- CreateAssayObject(counts = vst.out$umi_corrected)
    vst.out$umi_corrected <- NULL
  } else {
    # TODO: restore once check.matrix is in SeuratObject
    # assay.out <- CreateAssayObject(counts = umi, check.matrix = FALSE)
    assay.out <- CreateAssayObject(counts = umi)
  }
  # set the variable genes
  VariableFeatures(object = assay.out) <- vst.out$variable_features
  # put log1p transformed counts in data
  assay.out <- SetAssayData(
    object = assay.out,
    slot = 'data',
    new.data = log1p(x = GetAssayData(object = assay.out, slot = 'counts'))
  )
  scale.data <- vst.out$y
  assay.out <- SetAssayData(
    object = assay.out,
    slot = 'scale.data',
    new.data = scale.data
  )
  vst.out$y <- NULL
  # save clip.range into vst model
  vst.out$arguments$sct.clip.range <- clip.range
  vst.out$arguments$sct.method <- sct.method
  Misc(object = assay.out, slot = 'vst.out') <- vst.out
  assay.out <- as(object = assay.out, Class = "SCTAssay")
  return(assay.out)
}

#' @param assay Name of assay to pull the count data from; default is 'RNA'
#' @param new.assay.name Name for the new assay containing the normalized data; default is 'SCT'
#'
#' @rdname SCTransform
#' @concept preprocessing
#' @export
#' @method SCTransform Seurat
#'
SCTransform.Seurat <- function(
    object,
    assay = NULL,
    new.assay.name = 'SCT',
    reference.SCT.model = NULL,
    do.correct.umi = TRUE,
    ncells = 5000,
    residual.features = NULL,
    variable.features.n = 3000,
    variable.features.rv.th = 1.3,
    vars.to.regress = NULL,
    do.scale = FALSE,
    do.center = TRUE,
    clip.range = c(-sqrt(x = ncol(x = object[[assay]]) / 30), sqrt(x = ncol(x = object[[assay]]) / 30)),
    vst.flavor = "v2",
    conserve.memory = FALSE,
    return.only.var.genes = TRUE,
    seed.use = 1448145,
    verbose = TRUE,
    ...
) {
  assay <- assay %||% DefaultAssay(object = object)
  if (verbose){
    message("Running SCTransform on assay: ", assay)
  }
  cell.attr <- slot(object = object, name = 'meta.data')[colnames(object[[assay]]),]

  assay.data <- SCTransform(object = object[[assay]],
                            cell.attr = cell.attr,
                            reference.SCT.model = reference.SCT.model,
                            do.correct.umi = do.correct.umi,
                            ncells = ncells,
                            residual.features = residual.features,
                            variable.features.n = variable.features.n,
                            variable.features.rv.th = variable.features.rv.th,
                            vars.to.regress = vars.to.regress,
                            do.scale = do.scale,
                            do.center = do.center,
                            clip.range = clip.range,
                            vst.flavor = vst.flavor,
                            conserve.memory = conserve.memory,
                            return.only.var.genes = return.only.var.genes,
                            seed.use = seed.use,
                            verbose = verbose,
                            ...)
  assay.data <- SCTAssay(assay.data, assay.orig = assay)
  slot(object = slot(object = assay.data, name = "SCTModel.list")[[1]], name = "umi.assay") <- assay
  object[[new.assay.name]] <- assay.data

  if (verbose) {
    message(paste("Set default assay to", new.assay.name))
  }
  DefaultAssay(object = object) <- new.assay.name
  object <- LogSeuratCommand(object = object)
  return(object)
}

#' Subset a Seurat Object based on the Barcode Distribution Inflection Points
#'
#' This convenience function subsets a Seurat object based on calculated inflection points.
#'
#' See [CalculateBarcodeInflections()] to calculate inflection points and
#' [BarcodeInflectionsPlot()] to visualize and test inflection point calculations.
#'
#' @param object Seurat object
#'
#' @return Returns a subsetted Seurat object.
#'
#' @export
#' @concept preprocessing
#'
#' @author Robert A. Amezquita, \email{robert.amezquita@fredhutch.org}
#' @seealso \code{\link{CalculateBarcodeInflections}} \code{\link{BarcodeInflectionsPlot}}
#'
#' @examples
#' data("pbmc_small")
#' pbmc_small <- CalculateBarcodeInflections(
#'   object = pbmc_small,
#'   group.column = 'groups',
#'   threshold.low = 20,
#'   threshold.high = 30
#' )
#' SubsetByBarcodeInflections(object = pbmc_small)
#'
SubsetByBarcodeInflections <- function(object) {
  cbi.data <- Tool(object = object, slot = 'CalculateBarcodeInflections')
  if (is.null(x = cbi.data)) {
    stop("Barcode inflections not calculated, please run CalculateBarcodeInflections")
  }
  return(object[, cbi.data$cells_pass])
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Methods for Seurat-defined generics
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @param selection.method How to choose top variable features. Choose one of :
#' \itemize{
#'   \item{vst:}{ First, fits a line to the relationship of log(variance) and
#'   log(mean) using local polynomial regression (loess). Then standardizes the
#'   feature values using the observed mean and expected variance (given by the
#'   fitted line). Feature variance is then calculated on the standardized values
#'   after clipping to a maximum (see clip.max parameter).}
#'   \item{mean.var.plot (mvp):}{ First, uses a function to calculate average
#'   expression (mean.function) and dispersion (dispersion.function) for each
#'   feature. Next, divides features into num.bin (deafult 20) bins based on
#'   their average expression, and calculates z-scores for dispersion within
#'   each bin. The purpose of this is to identify variable features while
#'   controlling for the strong relationship between variability and average
#'   expression.}
#'   \item{dispersion (disp):}{ selects the genes with the highest dispersion values}
#' }
#' @param loess.span (vst method) Loess span parameter used when fitting the
#' variance-mean relationship
#' @param clip.max (vst method) After standardization values larger than
#' clip.max will be set to clip.max; default is 'auto' which sets this value to
#' the square root of the number of cells
#' @param mean.function Function to compute x-axis value (average expression).
#'  Default is to take the mean of the detected (i.e. non-zero) values
#' @param dispersion.function Function to compute y-axis value (dispersion).
#' Default is to take the standard deviation of all values
#' @param num.bin Total number of bins to use in the scaled analysis (default
#' is 20)
#' @param binning.method Specifies how the bins should be computed. Available
#' methods are:
#' \itemize{
#'   \item{equal_width:}{ each bin is of equal width along the x-axis [default]}
#'   \item{equal_frequency:}{ each bin contains an equal number of features (can
#'   increase statistical power to detect overdispersed features at high
#'   expression values, at the cost of reduced resolution along the x-axis)}
#' }
#' @param verbose show progress bar for calculations
#'
#' @rdname FindVariableFeatures
#' @concept preprocessing
#' @export
#'
FindVariableFeatures.V3Matrix <- function(
  object,
  selection.method = "vst",
  loess.span = 0.3,
  clip.max = 'auto',
  mean.function = FastExpMean,
  dispersion.function = FastLogVMR,
  num.bin = 20,
  binning.method = "equal_width",
  verbose = TRUE,
  ...
) {
  CheckDots(...)
  if (!inherits(x = object, 'Matrix')) {
    object <- as(object = as.matrix(x = object), Class = 'Matrix')
  }
  if (!inherits(x = object, what = 'dgCMatrix')) {
    object <- as.sparse(x = object)
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
    colnames(x = hvf.info) <- paste0('vst.', colnames(x = hvf.info))
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
    colnames(x = hvf.info) <- paste0('mvp.', c('mean', 'dispersion', 'dispersion.scaled'))
  }
  return(hvf.info)
}

#' @param nfeatures Number of features to select as top variable features;
#' only used when \code{selection.method} is set to \code{'dispersion'} or
#' \code{'vst'}
#' @param mean.cutoff A two-length numeric vector with low- and high-cutoffs for
#' feature means
#' @param dispersion.cutoff A two-length numeric vector with low- and high-cutoffs for
#' feature dispersions
#'
#' @rdname FindVariableFeatures
#' @concept preprocessing
#'
#' @importFrom utils head
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
  nfeatures = 2000,
  mean.cutoff = c(0.1, 8),
  dispersion.cutoff = c(1, Inf),
  verbose = TRUE,
  ...
) {
  if (length(x = mean.cutoff) != 2 || length(x = dispersion.cutoff) != 2) {
    stop("Both 'mean.cutoff' and 'dispersion.cutoff' must be two numbers")
  }
  if (selection.method == "vst") {
    data <- GetAssayData(object = object, slot = "counts")
    # if (ncol(x = data) < 1 || nrow(x = data) < 1) {
    if (IsMatrixEmpty(x = data)) {
      warning("selection.method set to 'vst' but count slot is empty; will use data slot instead")
      data <- GetAssayData(object = object, slot = "data")
    }
  } else {
    data <- GetAssayData(object = object, slot = "data")
  }
  hvf.info <- FindVariableFeatures(
    object = data,
    selection.method = selection.method,
    loess.span = loess.span,
    clip.max = clip.max,
    mean.function = mean.function,
    dispersion.function = dispersion.function,
    num.bin = num.bin,
    binning.method = binning.method,
    verbose = verbose,
    ...
  )
  object[names(x = hvf.info)] <- hvf.info
  hvf.info <- hvf.info[which(x = hvf.info[, 1, drop = TRUE] != 0), ]
  if (selection.method == "vst") {
    hvf.info <- hvf.info[order(hvf.info$vst.variance.standardized, decreasing = TRUE), , drop = FALSE]
  } else {
    hvf.info <- hvf.info[order(hvf.info$mvp.dispersion, decreasing = TRUE), , drop = FALSE]
  }
  selection.method <- switch(
    EXPR = selection.method,
    'mvp' = 'mean.var.plot',
    'disp' = 'dispersion',
    selection.method
  )
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
  vf.name <- ifelse(
    test = selection.method == 'vst',
    yes = 'vst',
    no = 'mvp'
  )
  vf.name <- paste0(vf.name, '.variable')
  object[vf.name] <- rownames(x = object[]) %in% top.features
  return(object)
}

#' @rdname FindVariableFeatures
#' @export
#' @method FindVariableFeatures SCTAssay
#'
FindVariableFeatures.SCTAssay <- function(
  object,
  nfeatures = 2000,
  ...
) {
  if (length(x = slot(object = object, name = "SCTModel.list")) > 1) {
    stop("SCT assay is comprised of multiple SCT models. To change the variable features, please set manually with VariableFeatures<-", call. = FALSE)
  }
  feature.attr <- SCTResults(object = object, slot = "feature.attributes")
  nfeatures <- min(nfeatures, nrow(x = feature.attr))
  top.features <- rownames(x = feature.attr)[order(feature.attr$residual_variance, decreasing = TRUE)[1:nfeatures]]
  VariableFeatures(object = object) <- top.features
  return(object)
}

#' @param assay Assay to use
#'
#' @rdname FindVariableFeatures
#' @concept preprocessing
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
  nfeatures = 2000,
  mean.cutoff = c(0.1, 8),
  dispersion.cutoff = c(1, Inf),
  verbose = TRUE,
  ...
) {
  assay <- assay[1L] %||% DefaultAssay(object = object)
  assay <- match.arg(arg = assay, Assays(object = object))
  assay.data <- FindVariableFeatures(
    object = object[[assay]],
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
    verbose = verbose,
    ...
  )
  object[[assay]] <- assay.data
  if (inherits(x = object[[assay]], what = "SCTAssay")) {
    object <- GetResidual(
      object = object,
      assay = assay,
      features = VariableFeatures(object = assay.data),
      verbose = FALSE
    )
  }
  object <- LogSeuratCommand(object = object)
  return(object)
}

#' @param object A Seurat object, assay, or expression matrix
#' @param spatial.location Coordinates for each cell/spot/bead
#' @param selection.method Method for selecting spatially variable features.
#'  \itemize{
#'   \item \code{markvariogram}: See \code{\link{RunMarkVario}} for details
#'   \item \code{moransi}: See \code{\link{RunMoransI}} for details.
#' }
#'
#' @param r.metric r value at which to report the "trans" value of the mark
#' variogram
#' @param x.cuts Number of divisions to make in the x direction, helps define
#' the grid over which binning is performed
#' @param y.cuts Number of divisions to make in the y direction, helps define
#' the grid over which binning is performed
#' @param verbose Print messages and progress
#'
#' @method FindSpatiallyVariableFeatures default
#' @rdname FindSpatiallyVariableFeatures
#' @concept preprocessing
#' @concept spatial
#' @export
#'
#'
FindSpatiallyVariableFeatures.default <- function(
  object,
  spatial.location,
  selection.method = c('markvariogram', 'moransi'),
  r.metric = 5,
  x.cuts = NULL,
  y.cuts = NULL,
  verbose = TRUE,
  ...
) {
  # error check dimensions
  if (ncol(x = object) != nrow(x = spatial.location)) {
    stop("Please provide the same number of observations as spatial locations.")
  }
  if (!is.null(x = x.cuts) & !is.null(x = y.cuts)) {
    binned.data <- BinData(
      data = object,
      pos = spatial.location,
      x.cuts = x.cuts,
      y.cuts = y.cuts,
      verbose = verbose
    )
    object <- binned.data$data
    spatial.location <- binned.data$pos
  }
  svf.info <- switch(
    EXPR = selection.method,
    'markvariogram' = RunMarkVario(
      spatial.location = spatial.location,
      data = object
    ),
    'moransi' = RunMoransI(
      data = object,
      pos = spatial.location,
      verbose = verbose
    ),
    stop("Invalid selection method. Please choose one of: markvariogram, moransi.")
  )
  return(svf.info)
}

#' @param slot Slot in the Assay to pull data from
#' @param features If provided, only compute on given features. Otherwise,
#' compute for all features.
#' @param nfeatures Number of features to mark as the top spatially variable.
#'
#' @method FindSpatiallyVariableFeatures Assay
#' @rdname FindSpatiallyVariableFeatures
#' @concept preprocessing
#' @concept spatial
#' @export
#'
FindSpatiallyVariableFeatures.Assay <- function(
  object,
  slot = "scale.data",
  spatial.location,
  selection.method = c('markvariogram', 'moransi'),
  features = NULL,
  r.metric = 5,
  x.cuts = NULL,
  y.cuts = NULL,
  nfeatures = nfeatures,
  verbose = TRUE,
  ...
) {
  features <- features %||% rownames(x = object)
  if (selection.method == "markvariogram" && "markvariogram" %in% names(x = Misc(object = object))) {
    features.computed <- names(x = Misc(object = object, slot = "markvariogram"))
    features <- features[! features %in% features.computed]
  }
  data <- GetAssayData(object = object, slot = slot)
  data <- as.matrix(x = data[features, ])
  data <- data[RowVar(x = data) > 0, ]
  if (nrow(x = data) != 0) {
    svf.info <- FindSpatiallyVariableFeatures(
      object = data,
      spatial.location = spatial.location,
      selection.method = selection.method,
      r.metric = r.metric,
      x.cuts = x.cuts,
      y.cuts = y.cuts,
      verbose = verbose,
      ...
    )
  } else {
    svf.info <- c()
  }
  if (selection.method == "markvariogram") {
    if ("markvariogram" %in% names(x = Misc(object = object))) {
      svf.info <- c(svf.info, Misc(object = object, slot = "markvariogram"))
    }
    suppressWarnings(expr = Misc(object = object, slot = "markvariogram") <- svf.info)
    svf.info <- ComputeRMetric(mv = svf.info, r.metric)
    svf.info <- svf.info[order(svf.info[, 1]), , drop = FALSE]
  }
  if (selection.method == "moransi") {
    colnames(x = svf.info) <- paste0("MoransI_", colnames(x = svf.info))
    svf.info <- svf.info[order(svf.info[, 2], -abs(svf.info[, 1])), , drop = FALSE]
  }
  var.name <- paste0(selection.method, ".spatially.variable")
  var.name.rank <- paste0(var.name, ".rank")
  svf.info[[var.name]] <- FALSE
  svf.info[[var.name]][1:(min(nrow(x = svf.info), nfeatures))] <- TRUE
  svf.info[[var.name.rank]] <- 1:nrow(x = svf.info)
  object[names(x = svf.info)] <- svf.info
  return(object)
}

#' @param assay Assay to pull the features (marks) from
#' @param image Name of image to pull the coordinates from
#'
#' @method FindSpatiallyVariableFeatures Seurat
#' @rdname FindSpatiallyVariableFeatures
#' @concept preprocessing
#' @concept spatial
#' @export
#'
FindSpatiallyVariableFeatures.Seurat <- function(
  object,
  assay = NULL,
  slot = "scale.data",
  features = NULL,
  image = NULL,
  selection.method = c('markvariogram', 'moransi'),
  r.metric = 5,
  x.cuts = NULL,
  y.cuts = NULL,
  nfeatures = 2000,
  verbose = TRUE,
  ...
) {
  assay <- assay %||% DefaultAssay(object = object)
  features <- features %||% rownames(x = object[[assay]])
  image <- image %||% DefaultImage(object = object)
  tc <- GetTissueCoordinates(object = object[[image]])
  # check if markvariogram has been run on necessary features
  # only run for new ones
  object[[assay]] <- FindSpatiallyVariableFeatures(
    object = object[[assay]],
    slot = slot,
    features = features,
    spatial.location = tc,
    selection.method = selection.method,
    r.metric = r.metric,
    x.cuts = x.cuts,
    y.cuts = y.cuts,
    nfeatures = nfeatures,
    verbose = verbose,
    ...
  )
  object <- LogSeuratCommand(object = object)
}

#' @rdname LogNormalize
#' @method LogNormalize data.frame
#' @export
#'
LogNormalize.data.frame <- function(
  data,
  scale.factor = 1e4,
  verbose = TRUE,
  ...
) {
  return(LogNormalize(
    data = as.matrix(x = data),
    scale.factor = scale.factor,
    verbose = verbose,
    ...
  ))
}

#' @rdname LogNormalize
#' @method LogNormalize V3Matrix
#' @export
#'
LogNormalize.V3Matrix <- function(
  data,
  scale.factor = 1e4,
  verbose = TRUE,
  ...
) {
  # if (is.data.frame(x = data)) {
  #   data <- as.matrix(x = data)
  # }
  if (!inherits(x = data, what = 'dgCMatrix')) {
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

#' @importFrom future.apply future_lapply
#' @importFrom future nbrOfWorkers
#'
#' @param normalization.method Method for normalization.
#'  \itemize{
#'   \item{LogNormalize: }{Feature counts for each cell are divided by the total
#'   counts for that cell and multiplied by the scale.factor. This is then
#'   natural-log transformed using log1p.}
#'   \item{CLR: }{Applies a centered log ratio transformation}
#'   \item{RC: }{Relative counts. Feature counts for each cell are divided by the total
#'   counts for that cell and multiplied by the scale.factor. No log-transformation is applied.
#'   For counts per million (CPM) set \code{scale.factor = 1e6}}
#' }
#' @param scale.factor Sets the scale factor for cell-level normalization
#' @param margin If performing CLR normalization, normalize across features (1) or cells (2)
# @param across If performing CLR normalization, normalize across either "features" or "cells".
#' @param block.size How many cells should be run in each chunk, will try to split evenly across threads
#' @param verbose display progress bar for normalization procedure
#'
#' @rdname NormalizeData
#' @concept preprocessing
#' @export
#'
NormalizeData.V3Matrix <- function(
  object,
  normalization.method = "LogNormalize",
  scale.factor = 1e4,
  margin = 1,
  block.size = NULL,
  verbose = TRUE,
  ...
) {
  CheckDots(...)
  if (is.null(x = normalization.method)) {
    return(object)
  }
  normalized.data <- if (nbrOfWorkers() > 1) {
    norm.function <- switch(
      EXPR = normalization.method,
      'LogNormalize' = LogNormalize,
      'CLR' = CustomNormalize,
      'RC' = RelativeCounts,
      stop("Unknown normalization method: ", normalization.method)
    )
    if (normalization.method != 'CLR') {
      margin <- 2
    }
    tryCatch(
      expr = Parenting(parent.find = 'Seurat', margin = margin),
      error = function(e) {
        invisible(x = NULL)
      }
    )
    dsize <- switch(
      EXPR = margin,
      '1' = nrow(x = object),
      '2' = ncol(x = object),
      stop("'margin' must be 1 or 2")
    )
    chunk.points <- ChunkPoints(
      dsize = dsize,
      csize = block.size %||% ceiling(x = dsize / nbrOfWorkers())
    )
    normalized.data <- future_lapply(
      X = 1:ncol(x = chunk.points),
      FUN = function(i) {
        block <- chunk.points[, i]
        data <- if (margin == 1) {
          object[block[1]:block[2], , drop = FALSE]
        } else {
          object[, block[1]:block[2], drop = FALSE]
        }
        clr_function <- function(x) {
          return(log1p(x = x / (exp(x = sum(log1p(x = x[x > 0]), na.rm = TRUE) / length(x = x)))))
        }
        args <- list(
          data = data,
          scale.factor = scale.factor,
          verbose = FALSE,
          custom_function = clr_function, margin = margin
        )
        args <- args[names(x = formals(fun = norm.function))]
        return(do.call(
          what = norm.function,
          args = args
        ))
      }
    )
    do.call(
      what = switch(
        EXPR = margin,
        '1' = 'rbind',
        '2' = 'cbind',
        stop("'margin' must be 1 or 2")
      ),
      args = normalized.data
    )
  } else {
    switch(
      EXPR = normalization.method,
      'LogNormalize' = LogNormalize(
        data = object,
        scale.factor = scale.factor,
        verbose = verbose
      ),
      'CLR' = CustomNormalize(
        data = object,
        custom_function = function(x) {
          return(log1p(x = x / (exp(x = sum(log1p(x = x[x > 0]), na.rm = TRUE) / length(x = x)))))
        },
        margin = margin,
        verbose = verbose
        # across = across
      ),
      'RC' = RelativeCounts(
        data = object,
        scale.factor = scale.factor,
        verbose = verbose
      ),
      stop("Unkown normalization method: ", normalization.method)
    )
  }
  return(normalized.data)
}

#' @rdname NormalizeData
#' @concept preprocessing
#' @export
#' @method NormalizeData Assay
#'
NormalizeData.Assay <- function(
  object,
  normalization.method = "LogNormalize",
  scale.factor = 1e4,
  margin = 1,
  verbose = TRUE,
  ...
) {
  object <- SetAssayData(
    object = object,
    slot = 'data',
    new.data = NormalizeData(
      object = GetAssayData(object = object, slot = 'counts'),
      normalization.method = normalization.method,
      scale.factor = scale.factor,
      verbose = verbose,
      margin = margin,
      ...
    )
  )
  return(object)
}

#' @param assay Name of assay to use
#'
#' @rdname NormalizeData
#' @concept preprocessing
#' @export
#' @method NormalizeData Seurat
#'
#' @examples
#' \dontrun{
#' data("pbmc_small")
#' pbmc_small
#' pmbc_small <- NormalizeData(object = pbmc_small)
#' }
#'
NormalizeData.Seurat <- function(
  object,
  assay = NULL,
  normalization.method = "LogNormalize",
  scale.factor = 1e4,
  margin = 1,
  verbose = TRUE,
  ...
) {
  assay <- assay %||% DefaultAssay(object = object)
  assay.data <- NormalizeData(
    object = object[[assay]],
    normalization.method = normalization.method,
    scale.factor = scale.factor,
    verbose = verbose,
    margin = margin,
    ...
  )
  object[[assay]] <- assay.data
  object <- LogSeuratCommand(object = object)
  return(object)
}

#' @importFrom future nbrOfWorkers
#'
#' @param features Vector of features names to scale/center. Default is variable features.
#' @param vars.to.regress Variables to regress out (previously latent.vars in
#' RegressOut). For example, nUMI, or percent.mito.
#' @param latent.data Extra data to regress out, should be cells x latent data
#' @param split.by Name of variable in object metadata or a vector or factor defining
#' grouping of cells. See argument \code{f} in \code{\link[base]{split}} for more details
#' @param model.use Use a linear model or generalized linear model
#' (poisson, negative binomial) for the regression. Options are 'linear'
#' (default), 'poisson', and 'negbinom'
#' @param use.umi Regress on UMI count data. Default is FALSE for linear
#' modeling, but automatically set to TRUE if model.use is 'negbinom' or 'poisson'
#' @param do.scale Whether to scale the data.
#' @param do.center Whether to center the data.
#' @param scale.max Max value to return for scaled data. The default is 10.
#' Setting this can help reduce the effects of features that are only expressed in
#' a very small number of cells. If regressing out latent variables and using a
#' non-linear model, the default is 50.
#' @param block.size Default size for number of features to scale at in a single
#' computation. Increasing block.size may speed up calculations but at an
#' additional memory cost.
#' @param min.cells.to.block If object contains fewer than this number of cells,
#' don't block for scaling calculations.
#' @param verbose Displays a progress bar for scaling procedure
#'
#' @importFrom future.apply future_lapply
#'
#' @rdname ScaleData
#' @concept preprocessing
#' @export
#'
ScaleData.default <- function(
  object,
  features = NULL,
  vars.to.regress = NULL,
  latent.data = NULL,
  split.by = NULL,
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
  CheckDots(...)
  features <- features %||% rownames(x = object)
  features <- as.vector(x = intersect(x = features, y = rownames(x = object)))
  object <- object[features, , drop = FALSE]
  object.names <- dimnames(x = object)
  min.cells.to.block <- min(min.cells.to.block, ncol(x = object))
  suppressWarnings(expr = Parenting(
    parent.find = "ScaleData.Assay",
    features = features,
    min.cells.to.block = min.cells.to.block
  ))
  split.by <- split.by %||% TRUE
  split.cells <- split(x = colnames(x = object), f = split.by)
  CheckGC()
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
        t(x = object[vars.to.regress[vars.to.regress %in% rownames(x = object)], , drop=FALSE])
      )
    }
    # Currently, RegressOutMatrix will do nothing if latent.data = NULL
    notfound <- setdiff(x = vars.to.regress, y = colnames(x = latent.data))
    if (length(x = notfound) == length(x = vars.to.regress)) {
      stop(
        "None of the requested variables to regress are present in the object.",
        call. = FALSE
      )
    } else if (length(x = notfound) > 0) {
      warning(
        "Requested variables to regress not in object: ",
        paste(notfound, collapse = ", "),
        call. = FALSE,
        immediate. = TRUE
      )
      vars.to.regress <- colnames(x = latent.data)
    }
    if (verbose) {
      message("Regressing out ", paste(vars.to.regress, collapse = ', '))
    }
    chunk.points <- ChunkPoints(dsize = nrow(x = object), csize = block.size)
    if (nbrOfWorkers() > 1) { # TODO: lapply
      chunks <- expand.grid(
        names(x = split.cells),
        1:ncol(x = chunk.points),
        stringsAsFactors = FALSE
      )
      object <- future_lapply(
        X = 1:nrow(x = chunks),
        FUN = function(i) {
          row <- chunks[i, ]
          group <- row[[1]]
          index <- as.numeric(x = row[[2]])
          return(RegressOutMatrix(
            data.expr = object[chunk.points[1, index]:chunk.points[2, index], split.cells[[group]], drop = FALSE],
            latent.data = latent.data[split.cells[[group]], , drop = FALSE],
            features.regress = features,
            model.use = model.use,
            use.umi = use.umi,
            verbose = FALSE
          ))
        }
      )
      if (length(x = split.cells) > 1) {
        merge.indices <- lapply(
          X = 1:length(x = split.cells),
          FUN = seq.int,
          to = length(x = object),
          by = length(x = split.cells)
        )
        object <- lapply(
          X = merge.indices,
          FUN = function(x) {
            return(do.call(what = 'rbind', args = object[x]))
          }
        )
        object <- do.call(what = 'cbind', args = object)
      } else {
        object <- do.call(what = 'rbind', args = object)
      }
    } else {
      object <- lapply(
        X = names(x = split.cells),
        FUN = function(x) {
          if (verbose && length(x = split.cells) > 1) {
            message("Regressing out variables from split ", x)
          }
          return(RegressOutMatrix(
            data.expr = object[, split.cells[[x]], drop = FALSE],
            latent.data = latent.data[split.cells[[x]], , drop = FALSE],
            features.regress = features,
            model.use = model.use,
            use.umi = use.umi,
            verbose = verbose
          ))
        }
      )
      object <- do.call(what = 'cbind', args = object)
    }
    dimnames(x = object) <- object.names
    CheckGC()
  }
  if (verbose && (do.scale || do.center)) {
    msg <- paste(
      na.omit(object = c(
        ifelse(test = do.center, yes = 'centering', no = NA_character_),
        ifelse(test = do.scale, yes = 'scaling', no = NA_character_)
      )),
      collapse = ' and '
    )
    msg <- paste0(
      toupper(x = substr(x = msg, start = 1, stop = 1)),
      substr(x = msg, start = 2, stop = nchar(x = msg)),
      ' data matrix'
    )
    message(msg)
  }
  if (inherits(x = object, what = c('dgCMatrix', 'dgTMatrix'))) {
    scale.function <- FastSparseRowScale
  } else {
    object <- as.matrix(x = object)
    scale.function <- FastRowScale
  }
  if (nbrOfWorkers() > 1) {
    blocks <- ChunkPoints(dsize = length(x = features), csize = block.size)
    chunks <- expand.grid(
      names(x = split.cells),
      1:ncol(x = blocks),
      stringsAsFactors = FALSE
    )
    scaled.data <- future_lapply(
      X = 1:nrow(x = chunks),
      FUN = function(index) {
        row <- chunks[index, ]
        group <- row[[1]]
        block <- as.vector(x = blocks[, as.numeric(x = row[[2]])])
        arg.list <- list(
          mat = object[features[block[1]:block[2]], split.cells[[group]], drop = FALSE],
          scale = do.scale,
          center = do.center,
          scale_max = scale.max,
          display_progress = FALSE
        )
        arg.list <- arg.list[intersect(x = names(x = arg.list), y = names(x = formals(fun = scale.function)))]
        data.scale <- do.call(what = scale.function, args = arg.list)
        dimnames(x = data.scale) <- dimnames(x = object[features[block[1]:block[2]], split.cells[[group]]])
        suppressWarnings(expr = data.scale[is.na(x = data.scale)] <- 0)
        CheckGC()
        return(data.scale)
      }
    )
    if (length(x = split.cells) > 1) {
      merge.indices <- lapply(
        X = 1:length(x = split.cells),
        FUN = seq.int,
        to = length(x = scaled.data),
        by = length(x = split.cells)
      )
      scaled.data <- lapply(
        X = merge.indices,
        FUN = function(x) {
          return(suppressWarnings(expr = do.call(what = 'rbind', args = scaled.data[x])))
        }
      )
      scaled.data <- suppressWarnings(expr = do.call(what = 'cbind', args = scaled.data))
    } else {
      suppressWarnings(expr = scaled.data <- do.call(what = 'rbind', args = scaled.data))
    }
  } else {
    scaled.data <- matrix(
      data = NA_real_,
      nrow = nrow(x = object),
      ncol = ncol(x = object),
      dimnames = object.names
    )
    max.block <- ceiling(x = length(x = features) / block.size)
    for (x in names(x = split.cells)) {
      if (verbose) {
        if (length(x = split.cells) > 1 && (do.scale || do.center)) {
          message(gsub(pattern = 'matrix', replacement = 'from split ', x = msg), x)
        }
        pb <- txtProgressBar(min = 0, max = max.block, style = 3, file = stderr())
      }
      for (i in 1:max.block) {
        my.inds <- ((block.size * (i - 1)):(block.size * i - 1)) + 1
        my.inds <- my.inds[my.inds <= length(x = features)]
        arg.list <- list(
          mat = object[features[my.inds], split.cells[[x]], drop = FALSE],
          scale = do.scale,
          center = do.center,
          scale_max = scale.max,
          display_progress = FALSE
        )
        arg.list <- arg.list[intersect(x = names(x = arg.list), y = names(x = formals(fun = scale.function)))]
        data.scale <- do.call(what = scale.function, args = arg.list)
        dimnames(x = data.scale) <- dimnames(x = object[features[my.inds], split.cells[[x]]])
        scaled.data[features[my.inds], split.cells[[x]]] <- data.scale
        rm(data.scale)
        CheckGC()
        if (verbose) {
          setTxtProgressBar(pb = pb, value = i)
        }
      }
      if (verbose) {
        close(con = pb)
      }
    }
  }
  dimnames(x = scaled.data) <- object.names
  scaled.data[is.na(x = scaled.data)] <- 0
  CheckGC()
  return(scaled.data)
}

#' @rdname ScaleData
#' @concept preprocessing
#' @export
#' @method ScaleData IterableMatrix
#'
ScaleData.IterableMatrix <- function(
    object,
    features = NULL,
    do.scale = TRUE,
    do.center = TRUE,
    scale.max = 10,
    ...
) {
  features <- features %||% rownames(x = object)
  features <- as.vector(x = intersect(x = features, y = rownames(x = object)))
  object <- object[features, , drop = FALSE]
  if (do.center) {
    features.mean <- BPCells::matrix_stats(
      matrix = object,
      row_stats = 'mean')$row_stats['mean',]
  } else {
    features.mean <- 0
  }
  if (do.scale) {
    features.sd <- sqrt(BPCells::matrix_stats(
      matrix = object,
      row_stats = 'variance')$row_stats['variance',])
  } else {
    features.sd <- 1
  }
  if (scale.max != Inf) {
    object <- BPCells::min_by_row(mat = object, vals = scale.max * features.sd + features.mean)
  }
  scaled.data <- (object - features.mean) / features.sd
return(scaled.data)
}


#' @rdname ScaleData
#' @concept preprocessing
#' @export
#' @method ScaleData Assay
#'
ScaleData.Assay <- function(
  object,
  features = NULL,
  vars.to.regress = NULL,
  latent.data = NULL,
  split.by = NULL,
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
  if (length(x = features) == 0) {
    features <- rownames(x = GetAssayData(object = object, slot = slot.use))
  }
  object <- SetAssayData(
    object = object,
    slot = 'scale.data',
    new.data = ScaleData(
      object = GetAssayData(object = object, slot = slot.use),
      features = features,
      vars.to.regress = vars.to.regress,
      latent.data = latent.data,
      split.by = split.by,
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
  suppressWarnings(expr = Parenting(
    parent.find = "ScaleData.Seurat",
    features = features,
    min.cells.to.block = min.cells.to.block,
    use.umi = use.umi
  ))
  return(object)
}

#' @param assay Name of Assay to scale
#'
#' @rdname ScaleData
#' @concept preprocessing
#' @export
#' @method ScaleData Seurat
#'
ScaleData.Seurat <- function(
  object,
  features = NULL,
  assay = NULL,
  vars.to.regress = NULL,
  split.by = NULL,
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
  assay <- assay[1L] %||% DefaultAssay(object = object)
  assay <- match.arg(arg = assay, choices = Assays(object = object))
  if (any(vars.to.regress %in% colnames(x = object[[]]))) {
    latent.data <- object[[vars.to.regress[vars.to.regress %in% colnames(x = object[[]])]]]
  } else {
    latent.data <- NULL
  }
  if (is.character(x = split.by) && length(x = split.by) == 1) {
    split.by <- object[[split.by]]
  }
  assay.data <- ScaleData(
    # object = assay.data,
    object = object[[assay]],
    features = features,
    vars.to.regress = vars.to.regress,
    latent.data = latent.data,
    split.by = split.by,
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
  return(object)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Internal
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Read Vitessce Expression Data
#'
#' @inheritParams ReadVitessce
#'
#' @return An expression matrix with cells as columns and features as rows
#'
#' @name vitessce-helpers
#' @rdname vitessce-helpers
#'
#' @importFrom jsonlite read_json
#'
#' @keywords internal
#'
#' @noRd
#'
.ReadVitessceGenes <- function(counts) {
  p1 <- progressor()
  p1(
    message = "Reading counts in Vitessce genes format",
    class = 'sticky',
    amount = 0
  )
  p1(type = 'finish')
  cts <- read_json(path = counts)
  p2 <- progressor(steps = length(x = cts))
  cts <- lapply(
    X = names(x = cts),
    FUN = function(x) {
      expr <- cts[[x]]$cells
      expr <- as.matrix(x = expr)
      colnames(x = expr) <- x
      p2()
      return(expr)
    }
  )
  p2(type = 'finish')
  cts <- Reduce(
    f = function(x, y) {
      a <- merge(x = x, y = y, by = 0, all = TRUE)
      rownames(x = a) <- a$Row.names
      a$Row.names <- NULL
      return(as.matrix(x = a))
    },
    x = cts
  )
  cts[is.na(x = cts)] <- 0
  return(t(x = cts))
}

#' @name vitessce-helpers
#' @rdname vitessce-helpers
#'
#' @importFrom jsonlite read_json
#'
#' @keywords internal
#'
#' @noRd
#'
.ReadVitessceClusters <- function(counts) {
  p1 <- progressor()
  p1(
    message = "Reading counts in Vitessce clusters format",
    class = 'sticky',
    amount = 0
  )
  p1(type = 'finish')
  cts <- read_json(path = counts)
  # p2 <- progressor(steps = length(x = cts))
  cells <- unlist(x = cts$cols)
  features <- unlist(x = cts$rows)
  cts <- lapply(X = cts[['matrix']], FUN = unlist)
  cts <- t(x = as.data.frame(x = cts))
  dimnames(x = cts) <- list(features, cells)
  return(cts)
}


#' @name nanostring-helpers
#' @rdname nanostring-helpers
#'
#' @return data frame containing counts for cells based on a single class of segmentation (eg Nuclear)
#'
#' @keywords internal
#'
#' @noRd
#'
build.cellcomp.matrix <- function(mols.df, class=NULL) {
  if (!is.null(class)) {
    if (!(class %in% c("Nuclear", "Membrane", "Cytoplasm"))) {
      stop(paste("Cannot subset matrix based on segmentation:", class))
    }
    mols.df <- mols.df[mols.df$CellComp == class,]  # subset based on cell class
  }
  mols.df$bc <- paste0(as.character(mols.df$cell_ID), "_", as.character(mols.df$fov))
  ncol <- length(unique(mols.df$target))
  nrow <- length(unique(mols.df$bc))  # will mols.df already have a cell barcode column at this point
  mtx <- matrix(data=rep(0, nrow*ncol), nrow=nrow, ncol=ncol)
  colnames(mtx) <- unique(mols.df$target)
  rownames(mtx) <- unique(mols.df$bc)
  for (row in 1:nrow(mols.df)) {
    mol <- mols.df[row, "target"]
    bc <- mols.df[row, "bc"]
    mtx[bc, mol] <- mtx[bc, mol] + 1
  }
  return(as.data.frame(mtx))
}

# Bin spatial regions into grid and average expression values
#
# @param dat Expression data
# @param pos Position information/coordinates for each sample
# @param x.cuts Number of cuts to make in the x direction (defines grid along
# with y.cuts)
# @param y.cuts Number of cuts to make in the y direction
#
# @return returns a list with positions as centers of the bins and average
# expression within the bins
#
#' @importFrom Matrix rowMeans
#
BinData <- function(data, pos, x.cuts = 10, y.cuts = x.cuts, verbose = TRUE) {
  if (verbose) {
    message("Binning spatial data")
  }
  pos$x.cuts <- cut(x = pos[, 1], breaks = x.cuts)
  pos$y.cuts <- cut(x = pos[, 2], breaks = y.cuts)
  pos$bin <- paste0(pos$x.cuts, "_", pos$y.cuts)
  all.bins <- unique(x = pos$bin)
  new.pos <- matrix(data = numeric(), nrow = length(x = all.bins), ncol = 2)
  new.dat <- matrix(data = numeric(), nrow = nrow(x = data), ncol = length(x = all.bins))
  for(i in 1:length(x = all.bins)) {
    samples <- rownames(x = pos)[which(x = pos$bin == all.bins[i])]
    dat <- data[, samples]
    if (is.null(x = dim(x = dat))) {
      new.dat[, i] <- dat
    } else {
      new.dat[, i] <- rowMeans(data[, samples])
    }
    new.pos[i, 1] <- mean(pos[samples, "x"])
    new.pos[i, 2] <- mean(pos[samples, "y"])
  }
  rownames(x = new.dat) <- rownames(x = data)
  colnames(x = new.dat) <- all.bins
  rownames(x = new.pos) <- all.bins
  colnames(x = new.pos) <- colnames(x = pos)[1:2]
  return(list(data = new.dat, pos = new.pos))
}

# Sample classification from MULTI-seq
#
# Identify singlets, doublets and negative cells from multiplexing experiments.
#
# @param data Data frame with the raw count data (cell x tags)
# @param q Scale the data. Default is 1e4
#
# @return Returns a named vector with demultiplexed identities
#
#' @importFrom KernSmooth bkde
#' @importFrom stats approxfun quantile
#
# @author Chris McGinnis, Gartner Lab, UCSF
#
# @examples
# demux_result <- ClassifyCells(data = counts_data, q = 0.7)
#
ClassifyCells <- function(data, q) {
  ## Generate Thresholds: Gaussian KDE with bad barcode detection, outlier trimming
  ## local maxima estimation with bad barcode detection, threshold definition and adjustment
  # n_BC <- ncol(x = data)
  n_cells <- nrow(x = data)
  bc_calls <- vector(mode = "list", length = n_cells)
  n_bc_calls <- numeric(length = n_cells)
  for (i in 1:ncol(x = data)) {
    model <- tryCatch(
      expr = approxfun(x = bkde(x = data[, i], kernel = "normal")),
      error = function(e) {
        message("No threshold found for ", colnames(x = data)[i], "...")
      }
    )
    if (is.null(x = model)) {
      next
    }
    x <- seq.int(
      from = quantile(x = data[, i], probs = 0.001),
      to = quantile(x = data[, i], probs = 0.999),
      length.out = 100
    )
    extrema <- LocalMaxima(x = model(x))
    if (length(x = extrema) <= 1) {
      message("No threshold found for ", colnames(x = data)[i], "...")
      next
    }
    low.extremum <- min(extrema)
    high.extremum <- max(extrema)
    thresh <- (x[high.extremum] + x[low.extremum])/2
    ## Account for GKDE noise by adjusting low threshold to most prominent peak
    low.extremae <- extrema[which(x = x[extrema] <= thresh)]
    new.low.extremum <- low.extremae[which.max(x = model(x)[low.extremae])]
    thresh <- quantile(x = c(x[high.extremum], x[new.low.extremum]), probs = q)
    ## Find which cells are above the ith threshold
    cell_i <- which(x = data[, i] >= thresh)
    n <- length(x = cell_i)
    if (n == 0) { ## Skips to next BC if no cells belong to the ith group
      next
    }
    bc <- colnames(x = data)[i]
    if (n == 1) {
      bc_calls[[cell_i]] <- c(bc_calls[[cell_i]], bc)
      n_bc_calls[cell_i] <- n_bc_calls[cell_i] + 1
    } else {
      # have to iterate, lame
      for (cell in cell_i) {
        bc_calls[[cell]] <- c(bc_calls[[cell]], bc)
        n_bc_calls[cell] <- n_bc_calls[cell] + 1
      }
    }
  }
  calls <- character(length = n_cells)
  for (i in 1:n_cells) {
    if (n_bc_calls[i] == 0) { calls[i] <- "Negative"; next }
    if (n_bc_calls[i] > 1) { calls[i] <- "Doublet"; next }
    if (n_bc_calls[i] == 1) { calls[i] <- bc_calls[[i]] }
  }
  names(x = calls) <- rownames(x = data)
  return(calls)
}

# Computes the metric at a given r (radius) value and stores in meta.features
#
# @param mv Results of running markvario
# @param r.metric r value at which to report the "trans" value of the mark
# variogram
#
# @return Returns a data.frame with r.metric values
#
#
ComputeRMetric <- function(mv, r.metric = 5) {
  r.metric.results <- unlist(x = lapply(
    X = mv,
    FUN = function(x) {
      x$trans[which.min(x = abs(x = x$r - r.metric))]
    }
  ))
  r.metric.results <- as.data.frame(x = r.metric.results)
  colnames(r.metric.results) <- paste0("r.metric.", r.metric)
  return(r.metric.results)
}

# Normalize a given data matrix
#
# Normalize a given matrix with a custom function. Essentially just a wrapper
# around apply. Used primarily in the context of CLR normalization.
#
# @param data Matrix with the raw count data
# @param custom_function A custom normalization function
# @param margin Which way to we normalize. Set 1 for rows (features) or 2 for columns (genes)
# @parm across Which way to we normalize? Choose form 'cells' or 'features'
# @param verbose Show progress bar
#
# @return Returns a matrix with the custom normalization
#
#' @importFrom Matrix t
#' @importFrom methods as
#' @importFrom pbapply pbapply
#
CustomNormalize <- function(data, custom_function, margin, verbose = TRUE) {
  if (is.data.frame(x = data)) {
    data <- as.matrix(x = data)
  }
  if (!inherits(x = data, what = 'dgCMatrix')) {
    data <- as.sparse(x = data)
  }
  myapply <- ifelse(test = verbose, yes = pbapply, no = apply)
  # margin <- switch(
  #   EXPR = across,
  #   'cells' = 2,
  #   'features' = 1,
  #   stop("'across' must be either 'cells' or 'features'")
  # )
  if (verbose) {
    message("Normalizing across ", c('features', 'cells')[margin])
  }
  norm.data <- myapply(
    X = data,
    MARGIN = margin,
    FUN = custom_function)
  if (margin == 1) {
    norm.data = Matrix::t(x = norm.data)
  }
  colnames(x = norm.data) <- colnames(x = data)
  rownames(x = norm.data) <- rownames(x = data)
  return(norm.data)
}

# Inter-maxima quantile sweep to find ideal barcode thresholds
#
# Finding ideal thresholds for positive-negative signal classification per multiplex barcode
#
# @param call.list A list of sample classification result from different quantiles using ClassifyCells
#
# @return A list with two values: \code{res} and \code{extrema}:
# \describe{
#   \item{res}{A data.frame named res_id documenting the quantile used, subset, number of cells and proportion}
#   \item{extrema}{...}
# }
#
# @author Chris McGinnis, Gartner Lab, UCSF
#
# @examples
# FindThresh(call.list = bar.table_sweep.list)
#
FindThresh <- function(call.list) {
  # require(reshape2)
  res <- as.data.frame(x = matrix(
    data = 0L,
    nrow = length(x = call.list),
    ncol = 4
  ))
  colnames(x = res) <- c("q","pDoublet","pNegative","pSinglet")
  q.range <- unlist(x = strsplit(x = names(x = call.list), split = "q="))
  res$q <- as.numeric(x = q.range[grep(pattern = "0", x = q.range)])
  nCell <- length(x = call.list[[1]])
  for (i in 1:nrow(x = res)) {
    temp <- table(call.list[[i]])
    if ("Doublet" %in% names(x = temp) == TRUE) {
      res$pDoublet[i] <- temp[which(x = names(x = temp) == "Doublet")]
    }
    if ( "Negative" %in% names(temp) == TRUE ) {
      res$pNegative[i] <- temp[which(x = names(x = temp) == "Negative")]
    }
    res$pSinglet[i] <- sum(temp[which(x = !names(x = temp) %in% c("Doublet", "Negative"))])
  }
  res.q <- res$q
  q.ind <- grep(pattern = 'q', x = colnames(x = res))
  res <- Melt(x = res[, -q.ind])
  res[, 1] <- rep.int(x = res.q, times = length(x = unique(res[, 2])))
  colnames(x = res) <- c('q', 'variable', 'value')
  res[, 4] <- res$value/nCell
  colnames(x = res)[2:4] <- c("Subset", "nCells", "Proportion")
  extrema <- res$q[LocalMaxima(x = res$Proportion[which(x = res$Subset == "pSinglet")])]
  return(list(res = res, extrema = extrema))
}

# Calculate pearson residuals of features not in the scale.data
# This function is the secondary function under GetResidual
#
# @param object A seurat object
# @param features Name of features to add into the scale.data
# @param assay Name of the assay of the seurat object generated by SCTransform
# @param vst_out The SCT parameter list
# @param clip.range Numeric of length two specifying the min and max values the Pearson residual
# will be clipped to
# Useful if you want to change the clip.range.
# @param verbose Whether to print messages and progress bars
#
# @return Returns a matrix containing not-centered pearson residuals of added features
#
#' @importFrom sctransform get_residuals
#
GetResidualSCTModel <- function(
  object,
  assay,
  SCTModel,
  new_features,
  clip.range,
  replace.value,
  verbose
) {
  clip.range <- clip.range %||% SCTResults(object = object[[assay]], slot = "clips", model = SCTModel)$sct
  model.features <- rownames(x = SCTResults(object = object[[assay]], slot = "feature.attributes", model = SCTModel))
  umi.assay <- SCTResults(object = object[[assay]], slot = "umi.assay", model = SCTModel)
  model.cells <- Cells(x = slot(object = object[[assay]], name = "SCTModel.list")[[SCTModel]])
  sct.method <-  SCTResults(object = object[[assay]], slot = "arguments", model = SCTModel)$sct.method %||% "default"
  scale.data.cells <- colnames(x = GetAssayData(object = object, assay = assay, slot = "scale.data"))
  if (length(x = setdiff(x = model.cells, y =  scale.data.cells)) == 0) {
  existing_features <- names(x = which(x = ! apply(
    X = GetAssayData(object = object, assay = assay, slot = "scale.data")[, model.cells],
    MARGIN = 1,
    FUN = anyNA)
  ))
 } else {
   existing_features <- character()
 }
  if (replace.value) {
    features_to_compute <- new_features
  } else {
    features_to_compute <- setdiff(x = new_features, y = existing_features)
  }
  if (sct.method == "reference.model") {
    if (verbose) {
      message("sct.model ", SCTModel, " is from reference, so no residuals will be recalculated")
    }
    features_to_compute <- character()
  }
  if (!umi.assay %in% Assays(object = object)) {
    warning("The umi assay (", umi.assay, ") is not present in the object. ",
             "Cannot compute additional residuals.", call. = FALSE, immediate. = TRUE)
    return(NULL)
  }
  diff_features <- setdiff(x = features_to_compute, y = model.features)
  intersect_features <- intersect(x = features_to_compute, y = model.features)
  if (length(x = diff_features) == 0) {
    umi <- GetAssayData(object = object, assay = umi.assay, slot = "counts" )[features_to_compute, model.cells, drop = FALSE]
  } else {
    warning(
      "In the SCTModel ", SCTModel, ", the following ", length(x = diff_features),
      " features do not exist in the counts slot: ", paste(diff_features, collapse = ", ")
    )
    if (length(x = intersect_features) == 0) {
      return(matrix(
        data = NA,
        nrow = length(x = features_to_compute),
        ncol = length(x = model.cells),
        dimnames = list(features_to_compute, model.cells)
      ))
    }
    umi <- GetAssayData(object = object, assay = umi.assay, slot = "counts")[intersect_features, model.cells, drop = FALSE]
  }
  clip.max <- max(clip.range)
  clip.min <- min(clip.range)
  if (nrow(x = umi) > 0) {
    vst_out <- SCTModel_to_vst(SCTModel = slot(object = object[[assay]], name = "SCTModel.list")[[SCTModel]])
    if (verbose) {
      message("sct.model: ", SCTModel)
    }
    new_residual <- get_residuals(
      vst_out = vst_out,
      umi = umi,
      residual_type = "pearson",
      res_clip_range = c(clip.min, clip.max),
      verbosity = as.numeric(x = verbose) * 2
    )
    new_residual <- as.matrix(x = new_residual)
    # centered data
    new_residual <- new_residual - rowMeans(x = new_residual)
  } else {
    new_residual <- matrix(data = NA, nrow = 0, ncol = length(x = model.cells), dimnames = list(c(), model.cells))
  }
  old.features <- setdiff(x = new_features, y = features_to_compute)
  if (length(x = old.features) > 0) {
    old_residuals <- GetAssayData(object = object[[assay]], slot = "scale.data")[old.features, model.cells, drop = FALSE]
    new_residual <- rbind(new_residual, old_residuals)[new_features, ]
  }
  return(new_residual)
}

# Convert SCTModel class to vst_out used in the sctransform
# @param SCTModel
# @return Return a list containing sct model
#
SCTModel_to_vst <- function(SCTModel) {
  feature.params <- c("theta", "(Intercept)",  "log_umi")
  feature.attrs <- c("residual_mean", "residual_variance" )
  vst_out <- list()
  vst_out$model_str <- slot(object = SCTModel, name = "model")
  vst_out$model_pars_fit <- as.matrix(x = slot(object = SCTModel, name = "feature.attributes")[, feature.params])
  vst_out$gene_attr <- slot(object = SCTModel, name = "feature.attributes")[, feature.attrs]
  vst_out$cell_attr <- slot(object = SCTModel, name = "cell.attributes")
  vst_out$arguments <- slot(object = SCTModel, name = "arguments")
  return(vst_out)
}

# Local maxima estimator
#
# Finding local maxima given a numeric vector
#
# @param x A continuous vector
#
# @return Returns a (named) vector showing positions of local maximas
#
# @author Tommy
# @references \url{https://stackoverflow.com/questions/6836409/finding-local-maxima-and-minima}
#
# @examples
# x <- c(1, 2, 9, 9, 2, 1, 1, 5, 5, 1)
# LocalMaxima(x = x)
#
LocalMaxima <- function(x) {
  # Use -Inf instead if x is numeric (non-integer)
  y <- diff(x = c(-.Machine$integer.max, x)) > 0L
  y <- cumsum(x = rle(x = y)$lengths)
  y <- y[seq.int(from = 1L, to = length(x = y), by = 2L)]
  if (x[[1]] == x[[2]]) {
    y <- y[-1]
  }
  return(y)
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
  if (is.numeric(x = fit)) {
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
# Remove unwanted effects from a matrix
#
# @parm data.expr An expression matrix to regress the effects of latent.data out
# of should be the complete expression matrix in genes x cells
# @param latent.data A matrix or data.frame of latent variables, should be cells
# x latent variables, the colnames should be the variables to regress
# @param features.regress An integer vector representing the indices of the
# genes to run regression on
# @param model.use Model to use, one of 'linear', 'poisson', or 'negbinom'; pass
# NULL to simply return data.expr
# @param use.umi Regress on UMI count data
# @param verbose Display a progress bar
#
#' @importFrom stats as.formula lm
#' @importFrom utils txtProgressBar setTxtProgressBar
#
RegressOutMatrix <- function(
  data.expr,
  latent.data = NULL,
  features.regress = NULL,
  model.use = NULL,
  use.umi = FALSE,
  verbose = TRUE
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
    pb <- txtProgressBar(char = '=', style = 3, file = stderr())
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
    data.resid <- log1p(x = Sweep(
      x = data.resid,
      MARGIN = 1,
      STATS = apply(X = data.resid, MARGIN = 1, FUN = min),
      FUN = '-'
    ))
  }
  dimnames(x = data.resid) <- dimnames(x = data.expr)
  return(data.resid)
}
