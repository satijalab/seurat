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
#' @param threshold.high Ignore barcodes of rank above this threshold in inflection calculation
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
    layer = 'counts'
  )
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
  average.expression <- suppressWarnings(
    AverageExpression(
      object = object,
      assays = assay,
      verbose = FALSE
    )[[assay]]
  )
  #checking for any cluster with all zero counts for any barcode
  if (sum(average.expression == 0) > 0) {
    stop("Cells with zero counts exist as a cluster.")
  }
  #create a matrix to store classification result
  discrete <- GetAssayData(object = object, assay = assay)
  discrete[discrete > 0] <- 0
  # for each HTO, we will use the minimum cluster for fitting
  for (iter in rownames(x = data)) {
    values <- counts[iter, ]
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
#' \dontrun{
#' data("pbmc_small")
#' pbmc_small <- SCTransform(object = pbmc_small, variable.features.n = 20)
#' pbmc_small <- GetResidual(object = pbmc_small, features = c('MS4A1', 'TCL1A'))
#' }
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
        model.cells <- Cells(x = slot(object = object[[assay]], name = "SCTModel.list")[[x]])
        FetchResidualSCTModel(object = object[[assay]],
                              umi.object = object[[umi.assay]],
                              layer.cells = model.cells,
                              SCTModel = x,
                              new_features = features,
                              replace.value = replace.value,
                              clip.range = clip.range,
                              verbose = verbose)
      }
    )
  }
  existing.data <- GetAssayData(object = object, layer = 'scale.data', assay = assay)
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
    layer = "scale.data",
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
#' @param bin.size Specifies the bin sizes to read in, can include "polygons" to load segmentations. Defaults to c(16, 8)
#' @param filter.matrix Only keep spots that have been determined to be over
#' tissue
#' @param to.upper Converts all feature names to upper case. Can be useful when
#' analyses require comparisons between human and mouse gene names for example.
#' @param image \code{VisiumV1}/\code{VisiumV2} instance(s) - if a vector is
#' passed in it should be co-indexed with \code{`bin.size`}
#' @param segmentation.type Which segmentations to load (cell or nucleus) when bin.size includes "polygons".
#' Defaults to "cell".
#' @param compact Whether to store segmentations in \emph{only} the \code{sf.data} slot
#' in the corresponding Segmentation object (default TRUE) to save memory and processing time.
#' If FALSE, segmentations are also stored in \code{\link[sp]{sp}} format in addition to the \code{sf.data} slot.
#' @param image.name Name of the tissue image to be plotted. Defaults to tissue_lowres_image.png
#' @param ... Arguments passed to \code{\link{Read10X_h5}}
#'
#' @return A \code{Seurat} object
#'
#' @importFrom png readPNG
#' @importFrom jsonlite fromJSON
#' @importFrom SeuratObject DefaultBoundary<-
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

Load10X_Spatial <- function (
  data.dir,
  filename = "filtered_feature_bc_matrix.h5",
  assay = "Spatial",
  slice = "slice1",
  bin.size = NULL,
  filter.matrix = TRUE,
  to.upper = FALSE,
  image = NULL,
  image.name = "tissue_lowres_image.png",
  segmentation.type = NULL,
  compact = TRUE,
  ...
) {
  # if more than one directory is passed in
  if (length(x = data.dir) > 1) {
    # party on with the first value
    data.dir <- data.dir[1]
    # but also raise a warning
    warning(
      paste0(
        "data.dir expects a single value but received multiple - ",
        "continuing using the first: '",
        data.dir,
        "'."
      ),
      immediate. = TRUE,
    )
  }
  # if the specified directory does not exist
  if (!file.exists(data.dir)) {
    # raise an error
    stop(paste0("No such file or directory: ", "'", data.dir, "'"))
  }

  # if bin.size is not set but data.dir points to a folder with binned data
  if (is.null(bin.size) & file.exists(paste0(data.dir, "/binned_outputs"))) {
    # point bin.size to the "standard" set - i.e. everything in the default
    # output except the 2 um binning because it's a memory hog
    bin.size <- c(16, 8)
  }

  # Seurat object to return
  object <- NULL

  bin.size.numeric <- bin.size

  if (!is.null(bin.size)) {
    # Store numeric bin sizes - these are used to load binned outputs
    bin.size.numeric <- as.numeric(bin.size[bin.size != "polygons"])
  }

  # Set flag to indicate if segmentations should be loaded
  load.segmentations <- length(bin.size.numeric) != length(bin.size)

  # read h5 files if bin.size is NULL (occurs when no /binned_outputs directory exists) or if bin.size contains numeric values
  load.binned.outputs <- is.null(bin.size) || (!is.null(bin.size.numeric) && length(bin.size.numeric) > 0)

  # If bin.size is specified and binned outputs need to be loaded
  if(!is.null(bin.size.numeric) && length(bin.size.numeric) > 0) {
    # convert bin.size to a character vector and pad values to three digits
    bin.size.pretty <- paste0(sprintf("%03d", bin.size.numeric), "um")
    # point data.dirs to the specified binnings
    data.dirs <- paste0(
      data.dir,
      "/binned_outputs/",
      "square_",
      bin.size.pretty
    )
    # suffix assay/slice names with each bin size
    assay.names <- paste0(assay, ".", bin.size.pretty)
    slice.names <- paste0(slice, ".", bin.size.pretty)
  } else {
    # otherwise just hold onto the top-level directory
    data.dirs <- data.dir
    # and keep the assay/slice names unchanged
    assay.names <- assay
    slice.names <- slice
  }

  # read the h5 files in the top-level / binned output directory
  if(load.binned.outputs) {
    # read in counts matrices from specified h5 files
    counts.paths <- lapply(data.dirs, file.path, filename)
    counts.list <- lapply(counts.paths, Read10X_h5, ...)
    # maybe convert Cell identifiers to uppercase
    if (to.upper) {
      rownames(counts) <- lapply(rownames(counts), toupper)
    }

    if (is.null(image)) {
      # read in the corresponding images and coordinate mappings
      image.list <- mapply(
        Read10X_Image,
        file.path(data.dirs, "spatial"),
        assay = assay.names,
        slice = slice.names,
        image.name = image.name,
        MoreArgs = list(filter.matrix = filter.matrix)
      )
    } else {
      # make sure any passed images are in a vector
      image.list <- c(image)
    }

    # check that for each counts matrix there is a corresponding image
    if (length(image.list) != length(counts.list)) {
      stop(
        paste0(
          "The number of images does not match the number of counts matrices. ",
          "Ensure each spatial dataset has a corresponding image."
        )
      )
    }

    # for each counts matrix, build a Seurat object
    object.list <- mapply(CreateSeuratObject, counts.list, assay = assay.names)
    # associate each counts matrix with its corresponding image
    object.list <- mapply(
      function(
        .object,
        .image,
        .assay,
        .slice
      ) {
        # align the image's identifiers with the object's
        .image <- .image[Cells(.object)]
        # add the image to the corresponding Seurat instance
        .object[[.slice]] <- .image
        return (.object)
      },
      object.list,
      image.list,
      assay.names,
      slice.names
    )
    # merge the Seurat instances - each assay should have unique Cell identifiers
    object <- merge(
      object.list[[1]],
      y = object.list[-1]
    )
  }

  # read segmentation data if requested
  if (load.segmentations) {
    # Check for required packages, stop with clear message if missing
    if (!requireNamespace("sf", quietly = TRUE)) {
      stop("The 'sf' package must be installed to load segmentation data.")
    }

    segmentation.assay.name <- paste0(assay, ".Polygons")
    seg.data.dir <- file.path(data.dir, "segmented_outputs")

    # Check for different possible file formats/names
    possible.files <- c(
      "filtered_feature_cell_matrix.h5",
      "filtered_feature_bc_matrix.h5",
      "raw_feature_bc_matrix.h5"
    )

    seg.counts.path <- NULL
    for (pf in possible.files) {
      test.path <- file.path(seg.data.dir, pf)
      if (file.exists(test.path)) {
        seg.counts.path <- test.path
        break
      }
    }

    if (is.null(seg.counts.path)) {
      stop("No cell segmentation matrix found. Looked for: ", paste(possible.files, collapse = ", "))
    }

    # Read raw counts matrix
    segmentation.counts <- Read10X_h5(seg.counts.path, ...)

    # Holds barcode names
    segmentation.counts.cell.ids <- colnames(segmentation.counts)

    # Check segmentation type
    if (is.null(segmentation.type)) {
      segmentation.type <- "cell"
    } else if (length(segmentation.type) != 1 || !(segmentation.type %in% c("cell", "nucleus"))) {
      stop("segmentation.type must be either 'cell' or 'nucleus'")
    }

    # Read the Visium (V2) object with segmentations loaded
    visium.segmentation <- Read10X_Segmentations(
      image.dir = file.path(seg.data.dir, "spatial"),
      data.dir = data.dir,
      image.name = image.name,
      segmentation.type = segmentation.type,
      compact = compact
    )

    # Create a new Seurat object with the raw counts
    segmentation.object <- CreateSeuratObject(
      segmentation.counts,
      assay = segmentation.assay.name
    )

    common_cells <- unique(Cells(visium.segmentation)[Cells(visium.segmentation) %in% Cells(segmentation.object)])
    visium.segmentation <- subset(
      x = visium.segmentation,
      cells = common_cells
    )

    # Set the default boundary type to centroids for plotting
    DefaultBoundary(object = visium.segmentation) <- "centroids"

    # Add the Visium object with segmentations to the Seurat object holding counts
    segmentation.object[[paste0(slice, ".polygons")]] <- visium.segmentation

    # Merge segmented outputs into the object containing binned outputs, if it exists
    if (!is.null(object)) {
      object <- merge(x = object, y = segmentation.object)
    } else {
      object <- segmentation.object
    }
    DefaultAssay(object = object) <- segmentation.assay.name
  }

  return(object)
}



#' Read10x Probe Metadata
#'
#' This function reads the probe metadata from a 10x Genomics probe barcode matrix file in HDF5 format.
#'
#' @param data.dir The directory where the file is located.
#' @param filename The name of the file containing the raw probe barcode matrix in HDF5 format. The default filename is 'raw_probe_bc_matrix.h5'.
#'
#' @return Returns a data.frame containing the probe metadata.
#'
#' @export
#' @concept preprocessing
#'
Read10X_probe_metadata <- function(
  data.dir,
  filename = 'raw_probe_bc_matrix.h5'
) {
  if (isFALSE(x = requireNamespace('hdf5r', quietly = TRUE))) {
    stop("Please install hdf5r to read HDF5 files")
  }
  file.path = paste0(data.dir,"/", filename)
  if (!file.exists(file.path)) {
    stop("File not found")
  }
  infile <- hdf5r::H5File$new(filename = file.path, mode = 'r')
  if("matrix/features/probe_region" %in% hdf5r::list.objects(infile)) {
    probe.name <- infile[['matrix/features/name']][]
    probe.region<- infile[['matrix/features/probe_region']][]
    meta.data <- data.frame(probe.name, probe.region)
    return(meta.data)
  }
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

#' Load Curio Seeker data
#'
#' @param data.dir location of data directory that contains the counts matrix,
#' gene names, barcodes/beads, and barcodes/bead location files.
#' @param assay Name of assay to associate spatial data to
#'
#' @return A \code{\link{Seurat}} object
#'
#' @importFrom Matrix readMM
#'
#' @export
#' @concept preprocessing
#'
LoadCurioSeeker <- function(data.dir, assay = "Spatial") {
  # check and find input files
  if (length(x = data.dir) > 1) {
    warning("'LoadCurioSeeker' accepts only one 'data.dir'",
            immediate. = TRUE)
    data.dir <- data.dir[1]
  }
  mtx.file <- list.files(
    data.dir,
    pattern = "*MoleculesPerMatchedBead.mtx",
    full.names = TRUE)
  if (length(x = mtx.file) > 1) {
    warning("Multiple files matched the pattern '*MoleculesPerMatchedBead.mtx'",
            immediate. = TRUE)
  } else if (length(x = mtx.file) == 0) {
    stop("No file matched the pattern '*MoleculesPerMatchedBead.mtx'", call. = FALSE)
  }
  mtx.file <- mtx.file[1]
  barcodes.file <- list.files(
    data.dir,
    pattern = "*barcodes.tsv",
    full.names = TRUE)
  if (length(x = barcodes.file) > 1) {
    warning("Multiple files matched the pattern '*barcodes.tsv'",
            immediate. = TRUE)
  } else if (length(x = barcodes.file) == 0) {
    stop("No file matched the pattern '*barcodes.tsv'", call. = FALSE)
  }
  barcodes.file <- barcodes.file[1]
  genes.file <- list.files(
    data.dir,
    pattern = "*genes.tsv",
    full.names = TRUE)
  if (length(x = genes.file) > 1) {
    warning("Multiple files matched the pattern '*genes.tsv'",
            immediate. = TRUE)
  } else if (length(x = genes.file) == 0) {
    stop("No file matched the pattern '*genes.tsv'", call. = FALSE)
  }
  genes.file <- genes.file[1]
  coordinates.file <- list.files(
    data.dir,
    pattern = "*MatchedBeadLocation.csv",
    full.names = TRUE)
  if (length(x = coordinates.file) > 1) {
    warning("Multiple files matched the pattern '*MatchedBeadLocation.csv'",
            immediate. = TRUE)
  } else if (length(x = coordinates.file) == 0) {
    stop("No file matched the pattern '*MatchedBeadLocation.csv'", call. = FALSE)
  }
  coordinates.file <- coordinates.file[1]

  # load counts matrix and create seurat object
  mtx <- readMM(mtx.file)
  mtx <- as.sparse(mtx)
  barcodes <- read.csv(barcodes.file, header = FALSE)
  genes <- read.csv(genes.file, header = FALSE)
  colnames(mtx) <- barcodes$V1
  rownames(mtx) <- genes$V1
  object <- CreateSeuratObject(counts = mtx, assay = assay)

  # load positions of each bead and store in a SlideSeq object in images slot
  coords <- read.csv(coordinates.file)
  colnames(coords) <- c("cell", "x", "y")
  coords$y <- -coords$y
  rownames(coords) <- coords$cell
  coords$cell <- NULL
  image <- new(Class = 'SlideSeq', assay = assay, coordinates = coords)
  object[["Slice"]] <- image
  return(object)
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
#' @references \doi{10.1038/s41592-019-0433-8}
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
    layer = "data",
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
      if ("Protein Expression" %in% lvls) {
        message("Xenium protein expression detected, but no scaling factor",
                " is supplied with the MEX matrices (vs HDF5). The loaded",
                " matrix is scaled by a constant from the original values.")
      }
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
  if (isFALSE(x = requireNamespace('hdf5r', quietly = TRUE))) {
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

        # Apply scaling factor that was used when serializing.
        if ("Protein Expression" %in% types.unique) {
          if ("protein_scaling_factor" %in% hdf5r::h5attr_names(infile)) {
            apply_scaling_factor <- 1.0 / hdf5r::h5attr(infile, "protein_scaling_factor")
            message("Scaling 'Protein Expression' by ", apply_scaling_factor)
            sparse.mat[["Protein Expression"]] <- sparse.mat[["Protein Expression"]] * apply_scaling_factor
          }
        }
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
#' should include files \code{tissue_lowres_image.png},
#' \code{scalefactors_json.json} and \code{tissue_positions_list.csv}
#' @param image.name PNG file to read in
#' @param assay Name of associated assay
#' @param slice Name for the image, used to populate the instance's key
#' @param filter.matrix Filter spot/feature matrix to only include spots that
#' have been determined to be over tissue
#' @param image.type Image type to return, one of: "VisiumV1" or "VisiumV2"
#'
#' @return A \code{\link{VisiumV2}} object
#'
#' @seealso \code{\link{VisiumV2}} \code{\link{Load10X_Spatial}}
#'
#' @export
#' @concept preprocessing
#'
Read10X_Image <- function(
  image.dir,
  image.name = "tissue_lowres_image.png",
  assay = "Spatial",
  slice = "slice1",
  filter.matrix = TRUE,
  image.type = "VisiumV2"
) {
  # Validate the `image.type` parameter.
  image.type <- match.arg(image.type, choices = c("VisiumV1", "VisiumV2"))

  # Read in the H&E stain image.
  primary.path <- file.path(image.dir, image.name)
  fallback.path <- file.path(dirname(dirname(dirname(image.dir))), "spatial", image.name)

  image <- tryCatch({
    png::readPNG(primary.path)
  }, error = function(e) {
    if (file.exists(fallback.path)) {
      png::readPNG(fallback.path)
    } else {
      stop("Neither primary nor fallback image could be read:\n", primary.path, "\n", fallback.path)
    }
  })

  # Read in the scale factors.
  scale.factors <- Read10X_ScaleFactors(
    filename = file.path(image.dir, "scalefactors_json.json")
  )

  # Read in the tissue coordinates as a data.frame.
  coordinates <- Read10X_Coordinates(
    filename = Sys.glob(file.path(image.dir, "*tissue_positions*")),
    filter.matrix
  )

  # Use the `slice` value to populate a Seurat-style identifier for the image.
  key <- Key(slice, quiet = TRUE)

  # Return the specified `image.type`.
  if (image.type == "VisiumV1") {
    visium.v1 <- new(
      Class = image.type,
      assay = assay,
      key = key,
      coordinates = coordinates,
      scale.factors = scale.factors,
      image = image
    )

    # As of v5.1.0 `Radius.VisiumV1` no longer returns the value of the
    # `spot.radius` slot and instead calculates the value on the fly, but we
    # can populate the static slot in case it's depended on.
    visium.v1@spot.radius <- Radius(visium.v1)

    return(visium.v1)
  }

  # If `image.type` is not "VisiumV1" then it must be "VisiumV2".
  stopifnot(image.type == "VisiumV2")

  # Create an `sp` compatible `FOV` instance.
  fov <- CreateFOV(
    coordinates[, c("imagecol", "imagerow")],
    type = "centroids",
    radius = scale.factors[["spot"]],
    assay = assay,
    key = key
  )

  #### NOTE ####
  # The Visium coordinate system takes the origin to be in the top-left corner,
  # where the x-axis is horizontal and associated with the image column.
  # We mark this with the coords_x_orientation flag.
  # Older Visium objects in Seurat have a different system (x-axis vertical, etc),
  # which is updated after checking whether the flag is set (SeuratObject::UpdateSeuratObject).
  ###############

  # Build the final `VisiumV2` instance
  visium.v2 <- new(
    Class = "VisiumV2",
    boundaries = fov@boundaries,
    molecules = fov@molecules,
    assay = fov@assay,
    key = fov@key,
    image = image,
    scale.factors = scale.factors,
    coords_x_orientation = "horizontal"
  )

  return(visium.v2)
}

#' Load 10X Genomics Visium Tissue Positions
#'
#' @param filename Path to a \code{tissue_positions_list.csv} file
#' @param filter.matrix Filter spot/feature matrix to only include spots that
#' have been determined to be over tissue
#'
#' @return A data.frame
#'
#' @export
#' @concept preprocessing
#'
Read10X_Coordinates <- function(filename, filter.matrix) {
  # output columns names
  col.names <- c("barcodes", "tissue", "row", "col", "imagerow", "imagecol")

  # if the coordinate mappings are in a parquet file
  if (tools::file_ext(filename) == "parquet") {
    # `arrow` must be installed to read parquet files
    if (isFALSE(x = requireNamespace('arrow', quietly = TRUE))) {
      stop("Please install arrow to read parquet files")
    }

    # read in coordinates and conver the resulting tibble into a data.frame
    coordinates <- as.data.frame(arrow::read_parquet(filename))
    # normalize column names for consistency with other datatypes
    input.col.names <- c(
      "barcode",
      "in_tissue",
      "array_row",
      "array_col",
      "pxl_row_in_fullres",
      "pxl_col_in_fullres"
    )
    col.map <- stats::setNames(col.names, input.col.names)
    colnames(coordinates) <- ifelse(
      colnames(coordinates) %in% names(col.map),
      col.map[colnames(coordinates)],
      colnames(coordinates)
    )

    # set rownames to "barcodes" then drop the column
    rownames(coordinates) <- coordinates[["barcodes"]]
    coordinates[["barcodes"]] <- NULL

  } else {
    # the coordinate mappings must be in a CSV - read it in
    coordinates <- read.csv(
        file = filename,
        col.names = col.names,
        header = ifelse(
          # assume files calles "tissue_positions.csv" have headers, otherwise
          # assume they do not (i.e. "tissue_positions_list.csv")
          test = basename(filename) == "tissue_positions.csv",
          yes = TRUE,
          no = FALSE
        ),
        as.is = TRUE,
        row.names = 1
      )
  }

  # the `tissue` column should contain a boolean indicating whether or not a
  # spot sits on top of the the tissue sample - maybe filter spots that do not
  if (filter.matrix) {
    coordinates <- coordinates[which(coordinates$tissue == 1), , drop = FALSE]
  }

  return (coordinates)
}

#' Load 10X Genomics Visium Scale Factors
#'
#' @param filename Path to a \code{scalefactors_json.json} file
#'
#' @return A scalefactors object
#'
#' @export
#' @concept preprocessing
#'
Read10X_ScaleFactors <- function(filename) {
  raw.data <- jsonlite::fromJSON(file.path(filename))

  scale.factors <- scalefactors(
    spot = raw.data$spot_diameter_fullres,
    fiducial = raw.data$fiducial_diameter_fullres,
    hires = raw.data$tissue_hires_scalef,
    lowres = raw.data$tissue_lowres_scalef
  )

  return (scale.factors)
}

#' Load 10X Genomics Visium Cell Segmentations
#'
#' @param image.dir Path to directory with 10X Genomics visium image data;
#' @param data.dir Directory of the base spaceranger outs
#' @param image.name Name of the tissue image to be plotted. tissue_lowres_image.png or tissue_hires_image.png
#' @param assay Name of assay to associate segmentations to
#' @param slice Name of the slice to associate the segmentations to
#' @param segmentation.type Which segmentations to load, cell or nucleus. If using nucleus the full matrix from cells is still used
#' @param compact Whether to store segmentations in only the \code{sf.data} slot; see \code{\link{Load10X_Spatial}} for details
#'
#'
#' @return A VisiumV2 object with segmentations
#'
#' @export
#' @concept preprocessing
#'
Read10X_Segmentations <- function (image.dir,
                                   data.dir,
                                   image.name = "tissue_lowres_image.png",
                                   assay = "Spatial.Polygons",
                                   slice = "slice1.polygons",
                                   segmentation.type = "cell",
                                   compact = TRUE) {


  sf.obj <- Read10X_HD_GeoJson(data.dir = data.dir,
                                segmentation.type = segmentation.type)

  # Create a Segmentation object; populate it based on the coordinates from the sf object
  segmentations <- CreateSegmentation(sf.obj, compact = compact)

  # Create a Centroids object; populate it based on the centroids from the sf object
  centroids <- CreateCentroids(sf.obj,
                              nsides = Inf,
                              radius = NULL,
                              theta = 0)

  # Named list with segmentations and centroids
  boundaries <- list(segmentations = segmentations, centroids = centroids)

  # Get image, scale factors, key
  image <- png::readPNG(source = file.path(image.dir, image.name))
  scale.factors <- Read10X_ScaleFactors(filename = file.path(image.dir,
                                                             "scalefactors_json.json"))
  key <- Key(slice, quiet = TRUE)

  # Build VisiumV2 object
  visium.v2 <- new(
    Class = "VisiumV2",
    boundaries = boundaries,
    assay = assay,
    key = key,
    image = image,
    scale.factors = scale.factors,
    coords_x_orientation = "horizontal"
  )

  return(visium.v2)
}

#' Format 10X Genomics GeoJson cell IDs
#'
#' @param ids Vector of cell IDs to format
#' @param prefix Optional prefix string
#' @param suffix Optional suffix string
#' @param digits Number of digits to zero-pad
#'
#' A helper function to format cell IDs from the segmentation GeoJson to the same type as in the matrix.h5
#' The GeoJson has cell IDs as integers (eg 1). They need to be in the format cellid_000000001-1
#'
#' @return Vector of formatted cell IDs
Format10X_GeoJson_CellID <- function(ids, prefix = "cellid_", suffix = "-1", digits = 9) {
  format_string <- paste0("%0", as.integer(digits), "d")

  formatted_ids <- sapply(ids, function(id) {
    numeric_part <- sprintf(format_string, as.integer(id))
    paste0(prefix, numeric_part, suffix)
  })

  return(formatted_ids)
}

#' Load 10X Genomics GeoJson
#'
#' @param data.dir Path to the directory containing matrix data
#' @param segmentation.type Which segmentations to load, cell or nucleus. If using nucleus the full matrix from cells is still used
#'
#' @return An \code{sf} object containing polygon segmentations from the GeoJSON provided by 10x, formatted for downstream coordinate retrieval
#'
#' @export
#' @concept preprocessing
#'
Read10X_HD_GeoJson <- function(data.dir, segmentation.type = "cell") {
  segmentation_polygons <- sf::read_sf(file.path(data.dir,"segmented_outputs", paste0(segmentation.type, "_segmentations.geojson")))

  # Restructure sf geometry for downstream compatibility
  segmentation_polygons$geometry <- sf::st_sfc(lapply(
    segmentation_polygons$geometry,
    function(geom) {
      coords <- geom[[1]]
      sf::st_polygon(list(coords))
    }
  ), crs = sf::st_crs(NA))

  segmentation_polygons$barcodes <- Format10X_GeoJson_CellID(segmentation_polygons$cell_id)
  segmentation_polygons
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
  if (isFALSE(x = requireNamespace('data.table', quietly = TRUE))) {
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
#' @param cell.sep Specify the delimiter in the cell name file
#' @param feature.sep Specify the delimiter in the feature name file
#' @param skip.cell Number of lines to skip in the cells file before beginning to read cell names
#' @param skip.feature Number of lines to skip in the features file before beginning to gene names
#' @param mtx.transpose Transpose the matrix after reading in
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
  cell.sep = "\t",
  feature.sep = "\t",
  skip.cell = 0,
  skip.feature = 0,
  mtx.transpose = FALSE,
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
    sep = cell.sep,
    row.names = NULL,
    skip = skip.cell
  )
  feature.names <- read.table(
    file = all.files[['feature list']],
    header = FALSE,
    sep = feature.sep,
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
  if (mtx.transpose) {
    data <- t(x = data)
  }
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
  if (isFALSE(x = requireNamespace('data.table', quietly = TRUE))) {
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

        tx <- as.data.frame(t(x = as.matrix(x = tx)))
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
#'  \item \dQuote{segmentation_method}: cell segmentation method (for runs which
#'  use multi-modal segmentation)
#' }
#' @param type Type of cell spatial coordinate matrices to read; choose one
#' or more of:
#' \itemize{
#'  \item \dQuote{centroids}: cell centroids in pixel coordinate space
#'  \item \dQuote{segmentations}: cell segmentations in pixel coordinate space
#'  \item \dQuote{nucleus_segmentations}: nucleus segmentations in pixel coordinate space
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
  outs = c("segmentation_method", "matrix", "microns"),
  type = "centroids",
  mols.qv.threshold = 20,
  flip.xy = F
) {
  # Argument checking
  type <- match.arg(
    arg = type,
    choices = c("centroids", "segmentations", "nucleus_segmentations"),
    several.ok = TRUE
  )

  outs <- match.arg(
    arg = outs,
    choices = c("segmentation_method", "matrix", "microns"),
    several.ok = TRUE
  )

  outs <- c(outs, type)

  has_dt <- requireNamespace("data.table", quietly = TRUE) && requireNamespace("R.utils", quietly = TRUE)
  has_arrow <- requireNamespace("arrow", quietly = TRUE)
  has_hdf5r <- requireNamespace("hdf5r", quietly = TRUE)

  binary_to_string <- function(arrow_binary) {
    if(typeof(arrow_binary) == 'list') {
      unlist(
        lapply(
          arrow_binary, function(x) rawToChar(as.raw(strtoi(x, 16L)))
        )
      )
    } else {
      arrow_binary
    }
  }

  data <- sapply(outs, function(otype) {
    switch(
      EXPR = otype,
      'matrix' = {
        pmtx <- progressor()
        pmtx(message = 'Reading counts matrix', class = 'sticky', amount = 0)

        for(option in Filter(function(x) x$req, list(
          list(filename = "cell_feature_matrix.h5", fn = Read10X_h5, req = has_hdf5r),
          list(filename = "cell_feature_matrix", fn = Read10X, req = TRUE)
        ))) {
          matrix <- try(suppressWarnings(option$fn(file.path(data.dir, option$filename))))
          if(!inherits(matrix, "try-error")) { break }
        }

        if(!exists('matrix') || inherits(matrix, "try-error")) {
          stop("Xenium outputs were incomplete: missing cell_feature_matrix")
        }

        pmtx(type = "finish")
        matrix
      },
      'segmentation_method' = {
        psegs <- progressor()
        psegs(
          message = 'Loading cell metadata',
          class = 'sticky',
          amount = 0
        )

        col.use <- c(
          cell_id = 'cell',
          segmentation_method = 'segmentation_method'
        )

        for(option in Filter(function(x) x$req, list(
          list(
            filename = "cells.parquet",
            fn = function(x) as.data.frame(arrow::read_parquet(x, col_select = names(col.use))),
            req = has_arrow
          ),
          list(
            filename = "cells.csv.gz",
            fn = function(x) data.table::fread(x, data.table = FALSE, stringsAsFactors = FALSE, select = names(col.use)),
            req = has_dt
          ),
          list(filename = "cells.csv.gz", fn = function(x) read.csv(x, stringsAsFactors = FALSE), req = TRUE)
        ))) {
          cell_seg <- try(suppressWarnings(option$fn(file.path(data.dir, option$filename))), silent = TRUE)
          if(!inherits(cell_seg, "try-error")) { break }
        }

        if (!exists("cell_seg") || inherits(cell_seg, "try-error")) {
          warning("cells did not contain a segmentation_method column. Skipping...", call. = FALSE, immediate. = TRUE)
          NULL
        } else {
          #Attempt to add default segmentation_method if nuclei/cell boundary files are provided
          message("Cell_seg columns: ", paste(colnames(cell_seg), collapse = ", "))
          if (!"segmentation_method" %in% colnames(cell_seg)) {
            message("Adding default segmentation_method = 'cell'")
            cell_seg$segmentation_method <- "cell"
          }

          #Try to detect unique cell identifier
          if (!"cell_id" %in% colnames(cell_seg)) {
            stop("Missing required column: cell_id")
          }

          cell_seg <- cell_seg[, c("cell_id", "segmentation_method")]
          colnames(cell_seg) <- c("cell", "segmentation_method")
          cell_seg$cell <- binary_to_string(cell_seg$cell)

          psegs(type = "finish")

          data.frame(segmentation_method = cell_seg$segmentation_method, row.names = cell_seg$cell)
        }
      },
      'centroids' = {
        pcents <- progressor()
        pcents(
          message = 'Loading cell centroids',
          class = 'sticky',
          amount = 0
        )

        col.use <- c(
          x_centroid = letters[24 + flip.xy],
          y_centroid = letters[25 - flip.xy],
          cell_id = 'cell'
        )

        for(option in Filter(function(x) x$req, list(
          list(
            filename = "cells.parquet",
            fn = function(x) as.data.frame(arrow::read_parquet(x, col_select = names(col.use))),
            req = has_arrow
          ),
          list(
            filename = "cells.csv.gz",
            fn = function(x) data.table::fread(x, data.table = FALSE, stringsAsFactors = FALSE, select = names(col.use)),
            req = has_dt
          ),
          list(filename = "cells.csv.gz", fn = function(x) read.csv(x, stringsAsFactors = FALSE), req = TRUE)
        ))) {
          cell_info <- try(suppressWarnings(option$fn(file.path(data.dir, option$filename))))
          if(!inherits(cell_info, "try-error")) { break }
        }

        if(!exists('cell_info') || inherits(cell_info, "try-error")) {
          stop("Xenium outputs were incomplete: missing cells")
        }

        cell_info$cell_id <- binary_to_string(cell_info$cell_id)

        cell_info <- cell_info[, names(col.use)]
        colnames(cell_info) <- col.use

        pcents(type = 'finish')

        cell_info
      },
      'segmentations' = {
        psegs <- progressor()
        psegs(
          message = 'Loading cell segmentations',
          class = 'sticky',
          amount = 0
        )

        for(option in Filter(function(x) x$req, list(
          list(
            filename = "cell_boundaries.parquet",
            fn = function(x) as.data.frame(arrow::read_parquet(x)),
            req = has_arrow
          ),
          list(
            filename = "cell_boundaries.csv.gz",
            fn = function(x) data.table::fread(x, data.table = FALSE, stringsAsFactors = FALSE),
            req = has_dt
          ),
          list(filename = "cell_boundaries.csv.gz", fn = function(x) read.csv(x, stringsAsFactors = FALSE), req = TRUE)
        ))) {
          cell_boundaries_df <- try(suppressWarnings(option$fn(file.path(data.dir, option$filename))))
          if(!inherits(cell_boundaries_df, "try-error")) { break }
        }

        if(!exists('cell_boundaries_df') || inherits(cell_boundaries_df, "try-error")) {
          stop("Xenium outputs were incomplete: missing cell_boundaries")
        }

        colnames(cell_boundaries_df) <- c(
          'cell',
          letters[24 + flip.xy],
          letters[25 - flip.xy]
        )

        cell_boundaries_df$cell <- binary_to_string(cell_boundaries_df$cell)

        psegs(type = "finish")

        cell_boundaries_df
      },
      'nucleus_segmentations' = {
        psegs <- progressor()
        psegs(
          message = 'Loading nucleus segmentations',
          class = 'sticky',
          amount = 0
        )

        for(option in Filter(function(x) x$req, list(
          list(
            filename = "nucleus_boundaries.parquet",
            fn = function(x) as.data.frame(arrow::read_parquet(x)),
            req = has_arrow
          ),
          list(
            filename = "nucleus_boundaries.csv.gz",
            fn = function(x) data.table::fread(x, data.table = FALSE, stringsAsFactors = FALSE),
            req = has_dt
          ),
          list(filename = "nucleus_boundaries.csv.gz", fn = function(x) read.csv(x, stringsAsFactors = FALSE), req = TRUE)
        ))) {
          nucleus_boundaries_df <- try(suppressWarnings(option$fn(file.path(data.dir, option$filename))))
          if(!inherits(nucleus_boundaries_df, "try-error")) { break }
        }

        if(!exists('nucleus_boundaries_df') || inherits(nucleus_boundaries_df, "try-error")) {
          stop("Xenium outputs were incomplete: missing nucleus_boundaries")
        }

        colnames(nucleus_boundaries_df) <- c(
          'cell',
          letters[24 + flip.xy],
          letters[25 - flip.xy]
        )
        nucleus_boundaries_df$cell <- binary_to_string(nucleus_boundaries_df$cell)

        psegs(type = "finish")

        nucleus_boundaries_df
      },
      'microns' = {
        pmicrons <- progressor()
        pmicrons(
          message = "Loading molecule coordinates",
          class = 'sticky',
          amount = 0
        )

        col.use = c(
          x_location = letters[24+flip.xy],
          y_location = letters[25-flip.xy],
          feature_name = 'gene'
        )

        for(option in Filter(function(x) x$req, list(
          list(
            filename = "transcripts.parquet",
            fn = function(x) as.data.frame(arrow::read_parquet(x, col_select = names(col.use))),
            req = has_arrow
          ),
          list(
            filename = "transcripts.csv.gz",
            fn = function(x) data.table::fread(x, data.table = FALSE, select = names(col.use), stringsAsFactors = FALSE),
            req = has_dt
          ),
          list(filename = "transcripts.csv.gz", fn = function(x) read.csv(x, stringsAsFactors = FALSE), req = TRUE)
        ))) {
          transcripts <- try(suppressWarnings(option$fn(file.path(data.dir, option$filename))))
          if(!inherits(transcripts, "try-error")) { break }
        }

        if(!exists('transcripts') || inherits(transcripts, "try-error")) {
          hint <- ""
          if(file.exists(file.path(data.dir, "transcripts.parquet"))) {
            hint <- ". Xenium outputs no longer include `transcripts.csv.gz`. Instead, please install `arrow` to read transcripts.parquet"
          }

          stop(paste0("Xenium outputs were incomplete: missing transcripts", hint))
        }

        transcripts <- transcripts[, names(col.use)]
        colnames(transcripts) <- col.use

        transcripts$gene <- binary_to_string(transcripts$gene)

        pmicrons(type = 'finish')

        transcripts
      },
      stop("Unknown Xenium input type: ", otype)
    )
  }, simplify = FALSE, USE.NAMES = TRUE)

  metadata <- file.path(data.dir, "experiment.xenium")
  if(file.exists(metadata) && requireNamespace("jsonlite", quietly = TRUE)) {
    meta <- jsonlite::read_json(metadata)
    data$metadata <- meta[
      intersect(
        names(meta),
        c(
          'run_start_time', 'preservation_method', 'panel_name',
          'panel_organism', 'panel_tissue_type',
          'instrument_sw_version', 'segmentation_stain'
        )
      )
    ]
  }
  return(data)
}

#' Check that the packages required to read Atera outputs are installed
#'
#' @keywords internal
#' @noRd
.AteraCheckDeps <- function() {
  pkgs <- c('blosc', 'jsonlite', 'data.table')
  have <- vapply(X = pkgs, FUN = requireNamespace, FUN.VALUE = logical(1), quietly = TRUE)
  if (!all(have)) {
    stop(
      "Reading Atera outputs requires the following package(s): ",
      paste(pkgs[!have], collapse = ', '),
      call. = FALSE
    )
  }
}

#' @keywords internal
#' @noRd
.AteraLeU16 <- function(raw, pos) {
  as.integer(raw[pos]) + as.integer(raw[pos + 1L]) * 256L
}

#' Widen a little-endian 4-byte field to a double, avoiding 32-bit signed
#' integer overflow (zip offsets/sizes routinely exceed 2^31)
#'
#' @keywords internal
#' @noRd
.AteraLeU32 <- function(raw, pos) {
  b <- as.integer(raw[pos:(pos + 3L)])
  b[1] + b[2] * 256 + b[3] * 65536 + b[4] * 16777216
}

#' @keywords internal
#' @noRd
.AteraLeU64 <- function(raw, pos) {
  b <- as.integer(raw[pos:(pos + 7L)])
  lo <- b[1] + b[2] * 256 + b[3] * 65536 + b[4] * 16777216
  hi <- b[5] + b[6] * 256 + b[7] * 65536 + b[8] * 16777216
  lo + hi * 4294967296
}

#' Locate the End Of Central Directory record of a zip file, following the
#' Zip64 EOCD locator/record when the classic EOCD's entry count/size/offset
#' fields are the \code{0xFFFF}/\code{0xFFFFFFFF} placeholder
#'
#' @keywords internal
#' @noRd
.AteraZipEOCD <- function(con, file.size) {
  eocd.sig <- as.raw(c(0x50, 0x4b, 0x05, 0x06))
  tail.size <- min(file.size, 22L + 65535L)
  seek(con, where = file.size - tail.size, origin = "start")
  tail.bytes <- readBin(con, what = "raw", n = tail.size)
  n <- length(tail.bytes)
  pos <- NA_integer_
  for (i in (n - 21L):1L) {
    if (i < 1L) break
    if (tail.bytes[i] == eocd.sig[1] && tail.bytes[i + 1L] == eocd.sig[2] &&
        tail.bytes[i + 2L] == eocd.sig[3] && tail.bytes[i + 3L] == eocd.sig[4]) {
      pos <- i
      break
    }
  }
  if (is.na(pos)) {
    stop("Could not locate End Of Central Directory record; not a valid zip file", call. = FALSE)
  }
  n.entries <- .AteraLeU16(tail.bytes, pos + 10L)
  cd.size <- .AteraLeU32(tail.bytes, pos + 12L)
  cd.offset <- .AteraLeU32(tail.bytes, pos + 16L)

  is.zip64 <- n.entries == 0xFFFF || cd.size >= 0xFFFFFFFF || cd.offset >= 0xFFFFFFFF
  if (is.zip64) {
    # the Zip64 EOCD locator is the fixed-size (20 byte) record immediately
    # preceding the EOCD record we just found; `pos` is a 1-based index into
    # tail.bytes, so the 0-based absolute file offset of the signature is
    # (file.size - tail.size) + pos - 1
    locator.abs.pos <- (file.size - tail.size) + pos - 21L
    seek(con, where = locator.abs.pos, origin = "start")
    locator <- readBin(con, what = "raw", n = 20L)
    zip64.eocd.offset <- .AteraLeU64(locator, 9L)
    seek(con, where = zip64.eocd.offset, origin = "start")
    zip64.eocd <- readBin(con, what = "raw", n = 56L)
    n.entries <- .AteraLeU64(zip64.eocd, 33L)
    cd.size <- .AteraLeU64(zip64.eocd, 41L)
    cd.offset <- .AteraLeU64(zip64.eocd, 49L)
  }
  list(n.entries = n.entries, cd.size = cd.size, cd.offset = cd.offset)
}

#' Parse the central directory of a zip file into a name -> (local header
#' offset, size) index, entirely in memory (no disk extraction). Atera
#' zarr.zip archives always store entries uncompressed (zip method
#' \dQuote{Stored}); the zarr/Blosc layer does all the compression, so this
#' index is all that's needed to seek directly to any chunk's raw bytes.
#'
#' @keywords internal
#' @noRd
.AteraZipIndex <- function(zip.file) {
  file.size <- file.info(zip.file)$size
  con <- file(zip.file, "rb")
  on.exit(close(con))
  eocd <- .AteraZipEOCD(con, file.size)

  seek(con, where = eocd$cd.offset, origin = "start")
  cd <- readBin(con, what = "raw", n = eocd$cd.size)

  cdfh.sig <- as.raw(c(0x50, 0x4b, 0x01, 0x02))
  names.out <- character(eocd$n.entries)
  offsets.out <- numeric(eocd$n.entries)
  sizes.out <- numeric(eocd$n.entries)
  methods.out <- integer(eocd$n.entries)

  pos <- 1L
  i <- 0L
  cd.len <- length(cd)
  while (pos <= cd.len - 45L) {
    if (!(cd[pos] == cdfh.sig[1] && cd[pos + 1L] == cdfh.sig[2] &&
          cd[pos + 2L] == cdfh.sig[3] && cd[pos + 3L] == cdfh.sig[4])) {
      break
    }
    i <- i + 1L
    method <- .AteraLeU16(cd, pos + 10L)
    csize <- .AteraLeU32(cd, pos + 20L)
    usize <- .AteraLeU32(cd, pos + 24L)
    fname.len <- .AteraLeU16(cd, pos + 28L)
    extra.len <- .AteraLeU16(cd, pos + 30L)
    comment.len <- .AteraLeU16(cd, pos + 32L)
    lho <- .AteraLeU32(cd, pos + 42L)

    name.start <- pos + 46L
    name <- rawToChar(cd[name.start:(name.start + fname.len - 1L)])

    if (extra.len > 0L) {
      extra.start <- name.start + fname.len
      extra <- cd[extra.start:(extra.start + extra.len - 1L)]
      # walk extra-field sub-records looking for the Zip64 tag (0x0001);
      # per the zip spec, only fields that were placeholder-valued in the
      # fixed header are present here, in order: usize, csize, offset
      ep <- 1L
      while (ep <= length(extra) - 3L) {
        tag <- .AteraLeU16(extra, ep)
        sz <- .AteraLeU16(extra, ep + 2L)
        if (tag == 1L) {
          dp <- ep + 4L
          if (usize >= 0xFFFFFFFF) { usize <- .AteraLeU64(extra, dp); dp <- dp + 8L }
          if (csize >= 0xFFFFFFFF) { csize <- .AteraLeU64(extra, dp); dp <- dp + 8L }
          if (lho >= 0xFFFFFFFF)   { lho   <- .AteraLeU64(extra, dp); dp <- dp + 8L }
          break
        }
        ep <- ep + 4L + sz
      }
    }

    names.out[i] <- name
    offsets.out[i] <- lho
    sizes.out[i] <- csize
    methods.out[i] <- method

    pos <- name.start + fname.len + extra.len + comment.len
  }

  if (i != eocd$n.entries) {
    names.out <- names.out[seq_len(i)]
    offsets.out <- offsets.out[seq_len(i)]
    sizes.out <- sizes.out[seq_len(i)]
    methods.out <- methods.out[seq_len(i)]
  }
  if (any(methods.out != 0L)) {
    stop("Atera zarr.zip entries are expected to be stored (uncompressed by zip); found a deflated entry", call. = FALSE)
  }

  list(
    name = names.out,
    offset = offsets.out,
    size = sizes.out,
    lookup = as.list(setNames(seq_along(names.out), names.out))
  )
}

#' Given a local file header offset (from the central directory), compute
#' the byte offset where the entry's raw data actually begins (skipping the
#' fixed 30-byte local header plus its variable name/extra fields)
#'
#' @keywords internal
#' @noRd
.AteraLocalDataOffset <- function(con, local.header.offset) {
  seek(con, where = local.header.offset, origin = "start")
  hdr <- readBin(con, what = "raw", n = 30L)
  fname.len <- .AteraLeU16(hdr, 27L)
  extra.len <- .AteraLeU16(hdr, 29L)
  local.header.offset + 30L + fname.len + extra.len
}

#' Read one zip entry's raw bytes directly into memory (no disk
#' extraction). \code{zidx} is a \code{.AteraZipIndex()} result; \code{con}
#' is an open \code{file(zip.file, "rb")} connection. Returns \code{NULL}
#' if the entry isn't present (eg a missing/all-fill-value zarr chunk).
#'
#' @keywords internal
#' @noRd
.AteraReadEntryRaw <- function(con, zidx, name) {
  i <- zidx$lookup[[name]]
  if (is.null(i)) {
    return(NULL)
  }
  data.offset <- .AteraLocalDataOffset(con, zidx$offset[i])
  seek(con, where = data.offset, origin = "start")
  readBin(con, what = "raw", n = zidx$size[i])
}

#' @keywords internal
#' @noRd
.AteraReadJSON <- function(con, zidx, name) {
  raw <- .AteraReadEntryRaw(con, zidx, name)
  if (is.null(raw)) {
    return(NULL)
  }
  jsonlite::fromJSON(rawToChar(raw), simplifyVector = TRUE)
}

#' R storage mode that a given zarr dtype string decodes to via
#' \code{blosc::blosc_decompress}, used to pre-allocate output vectors
#' without needing to decompress a chunk first. Note \code{<u4}/\code{<i8}/
#' \code{<u8} all decode to \code{double}, since R has no native type wide
#' enough to hold their full range losslessly.
#'
#' @keywords internal
#' @noRd
.AteraDtypeRType <- function(dtype) {
  switch(
    EXPR = dtype,
    "|i1" = , "|u1" = , "<i2" = , "<u2" = , "<i4" = "integer",
    "<u4" = , "<i8" = , "<u8" = , "<f2" = , "<f4" = , "<f8" = "double",
    "|b1" = "logical",
    stop("Unsupported zarr dtype: ", dtype, call. = FALSE)
  )
}

#' Decode one already-Blosc-compressed chunk's raw bytes into a typed R
#' vector, using the exact zarr dtype string from \code{.zarray} (eg
#' \code{"<i4"}, \code{"<u4"}, \code{"<f4"}, \code{"|b1"}) so that
#' \code{blosc}'s self-describing frame header drives decompression
#'
#' @keywords internal
#' @noRd
.AteraDecodeChunk <- function(raw, dtype, n) {
  vals <- blosc::blosc_decompress(raw, dtype = dtype)
  vals[seq_len(n)]
}

#' Number of worker processes to use for \code{.AteraReadArray}'s
#' Blosc-decompression step, honoring the package-wide \code{getThreads()}/
#' \code{setThreads()} option so it stays consistent with the rest of
#' Seurat's multithreading rather than defaulting to its own value.
#' \code{parallel::mclapply} is fork-based and unsupported on Windows, so
#' threading is disabled there instead of being requested and silently
#' downgraded.
#'
#' @keywords internal
#' @noRd
.AteraThreads <- function() {
  if (.Platform$OS.type == "windows") {
    return(1L)
  }
  min(getThreads(), parallel::detectCores(), na.rm = TRUE)
}

#' Read a zarr v2 array (1-D or 2-D), stored inside a zip archive, directly
#' into memory with no disk extraction. \code{row.range} (1-based,
#' inclusive \code{c(start, end)}) restricts which rows (first dimension)
#' are decoded and returned; only chunks overlapping that range are
#' read/decompressed, which is what makes gene-restricted transcript reads
#' cheap. Missing chunks (sparse zarr arrays) are filled with the array's
#' \code{fill_value}.
#'
#' @keywords internal
#' @noRd
.AteraReadArray <- function(con, zidx, array.path, row.range = NULL) {
  meta <- .AteraReadJSON(con, zidx, paste0(array.path, "/.zarray"))
  shape <- meta$shape
  # Coerced to double (chunk sizes parse as plain R integers, unlike `shape`,
  # which jsonlite already widens to double once it exceeds 32-bit range):
  # chunk-index * chunk-size products below can exceed 32-bit range for large
  # arrays (eg a multi-billion-element sparse matrix's X/data), and R's
  # native integer arithmetic silently overflows to NA rather than promoting
  # to double, so at least one operand of every such product must be double
  chunks <- as.double(meta$chunks)
  dtype <- meta$dtype
  order <- if (is.null(meta$order)) "C" else meta$order
  sep <- if (is.null(meta$dimension_separator)) "." else meta$dimension_separator
  fill.value <- if (is.null(meta$fill_value)) 0 else meta$fill_value
  ndim <- length(shape)

  if (ndim == 0L) {
    raw.chunk <- .AteraReadEntryRaw(con, zidx, paste0(array.path, "/0"))
    return(if (is.null(raw.chunk)) fill.value else .AteraDecodeChunk(raw.chunk, dtype, 1L))
  }

  # Chunk decoding is split into two phases so no I/O happens inside forked
  # workers: `con` is a single connection shared across this whole call, and
  # a `parallel::mclapply`/`mcmapply` fork inherits the *same* underlying
  # open file description, so concurrent seek()/readBin() calls from forked
  # children would race on that shared position and return corrupted bytes.
  # Phase 1 (below, sequential) reads each chunk's raw bytes through `con`;
  # phase 2 (`decode()`, possibly parallel) only touches those in-memory raw
  # vectors, which is safe to fork over.
  nthreads <- .AteraThreads()

  if (ndim == 1L) {
    n <- shape[1L]
    csize <- chunks[1L]
    if (is.null(row.range)) row.range <- c(1L, n)
    start <- row.range[1L]; end <- row.range[2L]
    c.first <- (start - 1L) %/% csize
    c.last <- (end - 1L) %/% csize
    chunk.idx <- c.first:c.last
    raw.chunks <- lapply(chunk.idx, function(ci) {
      .AteraReadEntryRaw(con, zidx, paste0(array.path, "/", ci))
    })
    actual.ns <- vapply(chunk.idx, function(ci) min(csize, n - ci * csize), numeric(1L))
    decode <- function(raw.chunk, actual.n) {
      if (is.null(raw.chunk)) rep(fill.value, actual.n) else .AteraDecodeChunk(raw.chunk, dtype, actual.n)
    }
    decoded <- if (length(raw.chunks) > 1L && nthreads > 1L) {
      parallel::mcmapply(decode, raw.chunks, actual.ns, SIMPLIFY = FALSE, mc.cores = nthreads)
    } else {
      Map(decode, raw.chunks, actual.ns)
    }
    out <- vector(mode = .AteraDtypeRType(dtype), length = end - start + 1L)
    for (k in seq_along(chunk.idx)) {
      ci <- chunk.idx[k]
      row0 <- ci * csize + 1L
      actual.n <- actual.ns[k]
      vals <- decoded[[k]]
      lo <- max(start, row0)
      hi <- min(end, row0 + actual.n - 1L)
      if (lo > hi) next
      out[(lo - start + 1L):(hi - start + 1L)] <- vals[(lo - row0 + 1L):(hi - row0 + 1L)]
    }
    return(out)
  }

  if (ndim == 2L) {
    nrow <- shape[1L]; ncol <- shape[2L]
    crow <- chunks[1L]; ccol <- chunks[2L]
    if (is.null(row.range)) row.range <- c(1L, nrow)
    start <- row.range[1L]; end <- row.range[2L]
    n.chunk.cols <- ceiling(ncol / ccol)
    c.first <- (start - 1L) %/% crow
    c.last <- (end - 1L) %/% crow

    jobs <- list()
    for (ri in c.first:c.last) {
      row0 <- ri * crow + 1L
      actual.nrow <- min(crow, nrow - ri * crow)
      lo <- max(start, row0)
      hi <- min(end, row0 + actual.nrow - 1L)
      if (lo > hi) next
      for (ci in seq_len(n.chunk.cols) - 1L) {
        col0 <- ci * ccol + 1L
        actual.ncol <- min(ccol, ncol - ci * ccol)
        key <- paste0(array.path, "/", ri, sep, ci)
        raw.chunk <- .AteraReadEntryRaw(con, zidx, key)
        jobs[[length(jobs) + 1L]] <- list(
          raw = raw.chunk,
          n = actual.nrow * actual.ncol,
          nrow = actual.nrow,
          ncol = actual.ncol,
          row0 = row0, col0 = col0, lo = lo, hi = hi
        )
      }
    }

    decode <- function(job) {
      if (is.null(job$raw)) rep(fill.value, job$n) else .AteraDecodeChunk(job$raw, dtype, job$n)
    }
    decoded <- if (length(jobs) > 1L && nthreads > 1L) {
      parallel::mclapply(jobs, decode, mc.cores = nthreads)
    } else {
      lapply(jobs, decode)
    }

    out <- matrix(vector(mode = .AteraDtypeRType(dtype), length = 1L), nrow = end - start + 1L, ncol = ncol)
    for (k in seq_along(jobs)) {
      job <- jobs[[k]]
      vals <- decoded[[k]]
      chunk.mat <- if (identical(order, "F")) {
        matrix(vals, nrow = job$nrow, ncol = job$ncol)
      } else {
        t(matrix(vals, nrow = job$ncol, ncol = job$nrow))
      }
      out[(job$lo - start + 1L):(job$hi - start + 1L), job$col0:(job$col0 + job$ncol - 1L)] <-
        chunk.mat[(job$lo - job$row0 + 1L):(job$hi - job$row0 + 1L), , drop = FALSE]
    }
    return(out)
  }

  stop("Unsupported array rank: ", ndim, call. = FALSE)
}

#' Decode a numcodecs VLenUTF8-filtered, already-Blosc-decompressed byte
#' buffer into a character vector: a \code{u32} count \code{N}, followed by
#' \code{N} \code{(u32 len, len raw utf8 bytes)} records, with no padding
#'
#' @keywords internal
#' @noRd
.AteraDecodeVlenUtf8 <- function(raw.bytes) {
  n <- as.integer(raw.bytes[1L]) + as.integer(raw.bytes[2L]) * 256L +
    as.integer(raw.bytes[3L]) * 65536L + as.integer(raw.bytes[4L]) * 16777216L
  strs <- character(n)
  pos <- 5L
  for (i in seq_len(n)) {
    len <- as.integer(raw.bytes[pos]) + as.integer(raw.bytes[pos + 1L]) * 256L +
      as.integer(raw.bytes[pos + 2L]) * 65536L + as.integer(raw.bytes[pos + 3L]) * 16777216L
    pos <- pos + 4L
    strs[i] <- rawToChar(raw.bytes[pos:(pos + len - 1L)])
    pos <- pos + len
  }
  strs
}

#' Read a \code{|O} (vlen-utf8) 1-D string array, chunk by chunk
#'
#' @keywords internal
#' @noRd
.AteraReadStringArray <- function(con, zidx, array.path) {
  meta <- .AteraReadJSON(con, zidx, paste0(array.path, "/.zarray"))
  n <- meta$shape[1L]
  csize <- meta$chunks[1L]
  n.chunks <- ceiling(n / csize)
  parts <- vector("list", n.chunks)
  for (ci in seq_len(n.chunks) - 1L) {
    key <- paste0(array.path, "/", ci)
    raw.chunk <- .AteraReadEntryRaw(con, zidx, key)
    dec <- blosc::blosc_decompress(raw.chunk)
    parts[[ci + 1L]] <- .AteraDecodeVlenUtf8(dec)
  }
  unlist(parts, use.names = FALSE)
}

#' Read an AnnData-style categorical column (a subgroup containing
#' \code{categories} and \code{codes} arrays) into an R factor
#'
#' @keywords internal
#' @noRd
.AteraReadCategorical <- function(con, zidx, group.path) {
  cat.meta <- .AteraReadJSON(con, zidx, paste0(group.path, "/categories/.zarray"))
  categories <- if (identical(cat.meta$dtype, "|O")) {
    .AteraReadStringArray(con, zidx, paste0(group.path, "/categories"))
  } else {
    .AteraReadArray(con, zidx, paste0(group.path, "/categories"))
  }
  codes <- .AteraReadArray(con, zidx, paste0(group.path, "/codes"))
  values <- rep_len(NA_character_, length(codes))
  keep <- codes >= 0
  values[keep] <- categories[codes[keep] + 1L]
  factor(values, levels = categories)
}

#' Decode a two-column (low, high) packed-uint32 id array (a matrix, as
#' returned by \code{.AteraReadArray}) into a single double id per row,
#' mapping the "all bits set" sentinel to NA (unassigned). Ids that were
#' stored as a plain (already-scalar) array are returned as-is.
#'
#' @keywords internal
#' @noRd
.AteraDecodePackedId <- function(ids) {
  if (!is.matrix(x = ids)) {
    return(ids)
  }
  sentinel <- 2^32 - 1
  unassigned <- ids[, 1] == sentinel & ids[, 2] == sentinel
  decoded <- ids[, 1] + ids[, 2] * 2^32
  decoded[unassigned] <- NA
  return(decoded)
}

#' Format a decoded numeric id as a string without falling back to
#' scientific notation, so ids remain usable as matrix/data frame join keys
#'
#' @keywords internal
#' @noRd
.AteraFormatId <- function(ids) {
  formatted <- sprintf(fmt = '%.0f', ids)
  formatted[is.na(x = ids)] <- NA
  return(formatted)
}

#' Enumerate the direct child entries (arrays or categorical subgroups) of
#' a zarr group from an already-parsed zip index, without any additional
#' I/O
#'
#' @keywords internal
#' @noRd
.AteraGroupEntries <- function(zidx, group.path) {
  prefix <- paste0(group.path, "/")
  matches <- zidx$name[startsWith(zidx$name, prefix)]
  rest <- substring(matches, nchar(prefix) + 1L)
  first.seg <- sub("/.*$", "", rest)
  setdiff(unique(first.seg), c(".zgroup", ".zattrs", ".zarray"))
}

#' Read a flat AnnData-style zarr group (eg \code{obs}/\code{var}) into a
#' data frame; plain columns are zarr arrays (string arrays use the
#' vlen-utf8 filter), categorical columns are subgroups containing
#' \code{categories} and \code{codes} arrays
#'
#' @keywords internal
#' @noRd
.AteraReadFlatGroup <- function(con, zidx, group.path, columns = NULL) {
  entries <- .AteraGroupEntries(zidx, group.path)
  if (!is.null(x = columns)) {
    entries <- intersect(x = entries, y = columns)
  }
  cols <- list()
  for (e in entries) {
    epath <- paste0(group.path, "/", e)
    if (!is.null(zidx$lookup[[paste0(epath, "/.zarray")]])) {
      meta <- .AteraReadJSON(con, zidx, paste0(epath, "/.zarray"))
      cols[[e]] <- if (identical(meta$dtype, "|O")) {
        .AteraReadStringArray(con, zidx, epath)
      } else {
        .AteraReadArray(con, zidx, epath)
      }
    } else if (!is.null(zidx$lookup[[paste0(epath, "/categories/.zarray")]])) {
      cols[[e]] <- .AteraReadCategorical(con, zidx, epath)
    }
  }
  return(as.data.frame(x = cols, stringsAsFactors = FALSE, check.names = FALSE))
}

#' Read one of Atera's flat (non-gridded) segmentation polygon sets
#' (\code{polygon_sets/0} = nucleus, \code{polygon_sets/1} = cell) into a
#' long-format data frame of \code{cell}/\code{x}/\code{y}, one row per
#' polygon vertex, joined against \code{cell.id} (a formatted, decoded
#' \code{cell_id} vector read from the same \code{cells.zarr.zip}) by
#' 0-based row index -- not assumed to be in the same row order as any
#' other file
#'
#' Every polygon's last stored vertex is a closing duplicate of its first
#' vertex (\code{x} always matches exactly), and in the majority of polygons
#' that duplicate's \code{y} is corrupted to exactly \code{0} in Atera's own
#' \code{vertices} array. Since consumers (eg \code{geom_polygon}, \code{sf})
#' already close rings back to the first vertex, that last vertex is dropped
#' here rather than propagating the corrupted coordinate.
#'
#' @keywords internal
#' @noRd
.AteraReadPolygonSet <- function(con, zidx, set.idx, cell.id) {
  base <- paste0("polygon_sets/", set.idx)
  cell_index <- .AteraReadArray(con, zidx, paste0(base, "/cell_index"))
  num_vertices <- .AteraReadArray(con, zidx, paste0(base, "/num_vertices")) - 1L
  vertices <- .AteraReadArray(con, zidx, paste0(base, "/vertices"))

  n <- nrow(vertices)
  poly.idx <- rep(seq_len(n), num_vertices)
  vert.idx <- sequence(num_vertices)
  col.x <- (vert.idx - 1L) * 2L + 1L
  col.y <- (vert.idx - 1L) * 2L + 2L

  data.frame(
    cell = cell.id[cell_index[poly.idx] + 1L],
    x = vertices[cbind(poly.idx, col.x)],
    y = vertices[cbind(poly.idx, col.y)]
  )
}

#' Build a cheap, reusable handle onto a \code{transcripts.zarr.zip}: parses
#' the zip's central directory and each grid tile's (tiny) \code{gene_offset}
#' index up front, but reads no \code{location}/\code{quality_score} data.
#' Used to defer the expensive part of transcript loading (decoding
#' \code{location}/\code{quality_score} chunks) until specific genes are
#' actually requested, so \code{molecule.coordinates = TRUE} doesn't have to
#' materialize the full transcript table just to be usable later.
#'
#' @keywords internal
#' @noRd
.AteraMoleculesHandle <- function(data.dir, mols.qv.threshold = 20) {
  zip.file <- file.path(data.dir, "transcripts.zarr.zip")
  zidx <- .AteraZipIndex(zip.file)
  con <- file(zip.file, "rb")
  on.exit(close(con))

  attrs <- .AteraReadJSON(con, zidx, ".zattrs")
  gene.names <- unlist(attrs$gene_names)

  grid.entries <- grep("^grid/[^/.]", zidx$name, value = TRUE)
  tile.dirs <- unique(vapply(
    strsplit(grid.entries, "/"),
    function(p) paste(p[1:2], collapse = "/"),
    character(1)
  ))
  tiles <- lapply(tile.dirs, function(tile.dir) {
    list(dir = tile.dir, gene_offset = .AteraReadArray(con, zidx, paste0(tile.dir, "/gene_offset")))
  })

  structure(
    list(
      zip.file = zip.file,
      zidx = zidx,
      gene.names = gene.names,
      tiles = tiles,
      mols.qv.threshold = mols.qv.threshold
    ),
    class = "AteraMoleculesHandle"
  )
}

#' Fetch transcript molecule coordinates from an \code{.AteraMoleculesHandle}
#' for \code{genes} (or all genes, if \code{NULL}), reading only the
#' (gene-sorted, contiguous) chunks needed for the requested genes out of
#' each grid tile
#'
#' A gene's rows are typically spread across most of a bundle's grid tiles
#' (eg ~85-88 of 102 for genes checked on a real whole-transcriptome bundle);
#' each tile's read is independent, so this is parallelized across tiles
#' rather than reading them one at a time. Each forked worker opens its own
#' connection, since a connection's file position can't safely be shared/
#' seeked concurrently across forked processes.
#'
#' @keywords internal
#' @noRd
.AteraFetchMolecules <- function(handle, genes = NULL) {
  nthreads <- .AteraThreads()

  read.tile <- function(tile) {
    con <- file(handle$zip.file, "rb")
    on.exit(close(con))
    gene_offset <- tile$gene_offset
    if (is.null(genes)) {
      location <- .AteraReadArray(con, handle$zidx, paste0(tile$dir, "/location"))
      quality_score <- .AteraReadArray(con, handle$zidx, paste0(tile$dir, "/quality_score"))
      gene.idx <- rep(seq_len(nrow(gene_offset)), gene_offset[, 2] - gene_offset[, 1])
      return(data.frame(
        x = location[, 1],
        y = location[, 2],
        gene = handle$gene.names[gene.idx],
        qv = as.vector(quality_score)
      ))
    }
    gene.rows <- match(genes, handle$gene.names)
    gene.rows <- gene.rows[!is.na(gene.rows) & (gene_offset[gene.rows, 2] - gene_offset[gene.rows, 1]) > 0]
    if (length(gene.rows) == 0) {
      return(data.frame(x = numeric(0), y = numeric(0), gene = character(0), qv = numeric(0)))
    }
    gene.dfs <- lapply(gene.rows, function(g) {
      rows <- c(gene_offset[g, 1] + 1L, gene_offset[g, 2])
      loc <- .AteraReadArray(con, handle$zidx, paste0(tile$dir, "/location"), row.range = rows)
      qv <- .AteraReadArray(con, handle$zidx, paste0(tile$dir, "/quality_score"), row.range = rows)
      data.frame(x = loc[, 1], y = loc[, 2], gene = handle$gene.names[g], qv = as.vector(qv))
    })
    data.table::rbindlist(gene.dfs)
  }

  tile.dfs <- if (nthreads > 1L) {
    parallel::mclapply(handle$tiles, read.tile, mc.cores = nthreads)
  } else {
    lapply(handle$tiles, read.tile)
  }

  df <- as.data.frame(data.table::rbindlist(tile.dfs))
  if (!is.null(handle$mols.qv.threshold)) {
    df <- df[!is.na(df$gene) & df$qv >= handle$mols.qv.threshold, , drop = FALSE]
  } else {
    df <- df[!is.na(df$gene), , drop = FALSE]
  }
  df$qv <- NULL
  df
}

#' Check that the \code{RBioFormats} package is installed; only called when
#' a morphology image is actually requested. The morphology OME-TIFFs are
#' pyramidal and JPEG2000-compressed, which the base \code{tiff} package
#' cannot reliably decode (it reports an "unknown" compression type and
#' cannot enumerate pyramid levels); \code{RBioFormats} wraps the Bio-Formats
#' Java library, which handles them correctly.
#'
#' @keywords internal
#' @noRd
.AteraCheckMorphologyDeps <- function() {
  if (!requireNamespace('RBioFormats', quietly = TRUE)) {
    stop(
      "Reading Atera morphology images requires the 'RBioFormats' package. ",
      "Install it with BiocManager::install('RBioFormats')",
      call. = FALSE
    )
  }
}

#' List the per-channel morphology OME-TIFF files in a
#' \dQuote{morphology_2d}/\dQuote{morphology_3d} directory, excluding macOS
#' resource-fork files (\code{._*})
#'
#' @keywords internal
#' @noRd
.AteraMorphologyFiles <- function(data.dir, three.d = FALSE) {
  dir <- file.path(data.dir, if (isTRUE(three.d)) "morphology_3d" else "morphology_2d")
  if (!dir.exists(dir)) {
    stop("No ", basename(dir), " directory found at ", dir, call. = FALSE)
  }
  files <- list.files(dir, pattern = "\\.ome\\.tif+$", full.names = TRUE)
  files <- files[!grepl("^\\._", basename(files))]
  if (!length(files)) {
    stop("No morphology OME-TIFF files found in ", dir, call. = FALSE)
  }
  files
}

#' Parse the 0-based channel index out of a \dQuote{chNNNN_<name>.ome.tif}
#' morphology image filename. Each file is itself a multi-page OME-TIFF
#' containing every channel of the panel, but only the page at this index
#' holds that channel's real image data (the rest are placeholders) -- see
#' \code{.AteraReadMorphologyImage}
#'
#' @keywords internal
#' @noRd
.AteraMorphologyChannelIndex <- function(filename) {
  m <- regmatches(basename(filename), regexec("^ch(\\d+)_.+\\.ome\\.tif+$", basename(filename)))[[1]]
  if (length(m) != 2) {
    stop(
      "Expected a morphology image filename of the form 'chNNNN_<name>.ome.tif', found ",
      basename(filename),
      call. = FALSE
    )
  }
  as.integer(m[2])
}

#' Parse the channel name out of a \dQuote{chNNNN_<name>.ome.tif} morphology
#' image filename
#'
#' @keywords internal
#' @noRd
.AteraMorphologyChannelName <- function(filename) {
  m <- regmatches(basename(filename), regexec("^ch\\d+_(.+)\\.ome\\.tif+$", basename(filename)))[[1]]
  if (length(m) != 2) {
    stop(
      "Expected a morphology image filename of the form 'chNNNN_<name>.ome.tif', found ",
      basename(filename),
      call. = FALSE
    )
  }
  m[2]
}

#' Build a cheap, reusable handle onto a channel's Atera morphology
#' OME-TIFF: resolves the matching file, its real-data channel index, and
#' per-resolution-level pixel dimensions/pixel size up front, but reads no
#' pixel data. Used to defer the expensive part (decoding image tiles) until
#' a specific region is actually requested via
#' \code{.AteraReadMorphologyRegion}, the same handle/fetch split already
#' used for transcripts by \code{.AteraMoleculesHandle}/
#' \code{.AteraFetchMolecules}.
#'
#' @keywords internal
#' @noRd
.AteraMorphologyHandle <- function(data.dir, channel = "dapi") {
  .AteraCheckMorphologyDeps()

  files <- .AteraMorphologyFiles(data.dir)
  names(files) <- vapply(files, .AteraMorphologyChannelName, character(1L))

  file <- if (is.numeric(channel)) {
    indices <- vapply(files, .AteraMorphologyChannelIndex, integer(1L))
    hits <- files[indices == as.integer(channel)]
    if (!length(hits)) {
      stop("No morphology image found for channel index ", channel, call. = FALSE)
    }
    hits[[1]]
  } else {
    hits <- names(files)[grepl(channel, names(files), ignore.case = TRUE)]
    if (!length(hits)) {
      stop(
        "No morphology image found matching channel '", channel, "'; available channels: ",
        paste(names(files), collapse = ", "),
        call. = FALSE
      )
    }
    files[[hits[1]]]
  }

  md <- RBioFormats::read.metadata(file)
  n.levels <- RBioFormats::seriesCount(md)
  cm <- RBioFormats::coreMetadata(md)
  # RBioFormats::coreMetadata() unwraps its usual per-series list and returns
  # a single series' metadata directly when there's only one series (eg a
  # non-pyramidal image with no resolution levels), so it must be re-wrapped
  # to keep the one-element-per-level shape `dims` below expects
  if (n.levels == 1L) {
    cm <- list(cm)
  }
  dims <- lapply(cm, function(x) c(x = x$sizeX, y = x$sizeY))

  specs.file <- file.path(data.dir, "experiment.spatial")
  pixel.size.full <- if (file.exists(specs.file) && requireNamespace('jsonlite', quietly = TRUE)) {
    jsonlite::read_json(specs.file)$pixel_size
  } else {
    NA_real_
  }

  structure(
    list(
      file = file,
      channel = .AteraMorphologyChannelName(file),
      channel.index = .AteraMorphologyChannelIndex(file),
      n.resolutions = n.levels,
      dim = dims,
      pixel.size.full = pixel.size.full
    ),
    class = "AteraMorphologyHandle"
  )
}

#' Convert a \code{region} (a list with optional \code{x}/\code{y} elements,
#' each a length-2 micron range) into 1-based pixel index ranges at a given
#' \code{.AteraMorphologyHandle}'s resolution \code{level}, clamped to the
#' plane's extent on that axis. A missing/\code{NULL} \code{region}, or a
#' missing \code{x}/\code{y} element, reads the whole plane along that axis.
#'
#' @keywords internal
#' @noRd
.AteraMorphologyRegionToPixels <- function(handle, level, region, pixel.size) {
  d <- handle$dim[[level]]
  to.pixels <- function(um.range, size) {
    if (is.null(um.range)) {
      return(NULL)
    }
    if (is.na(pixel.size)) {
      stop(
        "Cannot convert 'morphology.region' from microns to pixels: pixel size unavailable ",
        "(no 'experiment.spatial' file found)",
        call. = FALSE
      )
    }
    px <- round(um.range / pixel.size) + 1L
    px <- pmax(1L, pmin(size, px))
    px[1]:px[2]
  }
  list(
    x = to.pixels(region$x, d["x"]),
    y = to.pixels(region$y, d["y"])
  )
}

#' Read a (optionally cropped) single channel plane from an
#' \code{.AteraMorphologyHandle} at a given pyramid resolution level, by
#' default the lowest-resolution (smallest) level, which is normally all
#' that is needed for overview plots.
#'
#' \code{x.range}/\code{y.range} (1-based pixel indices at \code{resolution})
#' restrict the read to a rectangular window via
#' \code{RBioFormats::read.image}'s \code{subset} argument, which decodes
#' only the on-disk tiles overlapping that window rather than materializing
#' the whole plane -- needed because a full-resolution plane of these
#' whole-slide images can be tens of thousands of pixels per side (large
#' enough that reading it whole can exceed available Java heap memory).
#' Resolution levels are exposed by \code{RBioFormats} as separate
#' \dQuote{series} (\code{resolution} argument) of a single real image
#' series (\code{series = 1}), not as true multi-series data.
#'
#' @return A list with elements \code{image} (a numeric matrix, indexed
#' \verb{[x, y]} in the same top-left-origin pixel convention as the OME-TIFF
#' and, after scaling by \code{pixel.size} and offsetting by \code{origin},
#' the same convention as \code{ReadAtera}'s \code{centroids}/\code{microns}
#' micron coordinates -- no axis flip is needed), \code{channel} (the
#' matched channel name), \code{resolution} (the pyramid level read, 1 =
#' full resolution), \code{n.resolutions} (the number of pyramid levels
#' available), \code{pixel.size} (microns per pixel of \code{image}, i.e.
#' already scaled for the resolution level read; \code{NA} if
#' \code{experiment.spatial} could not be read), and \code{origin} (the
#' micron coordinates, \code{c(x=, y=)}, of \code{image}'s \verb{[1, 1]}
#' pixel -- \code{c(x=0, y=0)} unless \code{x.range}/\code{y.range} crop out
#' the top-left corner of the full plane)
#'
#' @keywords internal
#' @noRd
.AteraReadMorphologyRegion <- function(handle, resolution = NULL, x.range = NULL, y.range = NULL) {
  level <- resolution %||% handle$n.resolutions
  d <- handle$dim[[level]]
  downsample <- 2 ^ (level - 1L)
  pixel.size <- if (is.na(handle$pixel.size.full)) NA_real_ else handle$pixel.size.full * downsample

  x.range <- x.range %||% seq_len(d["x"])
  y.range <- y.range %||% seq_len(d["y"])

  img <- RBioFormats::read.image(
    handle$file, series = 1L, resolution = level, normalize = FALSE,
    subset = list(x = x.range, y = y.range, c = handle$channel.index + 1L)
  )
  ## requesting a single channel via 'subset' already drops the channel
  ## dimension (2D result); only index it away if it's still present (3D)
  if (length(dim(img)) == 3L) {
    img <- img[, , 1L]
  }

  list(
    image = img,
    channel = handle$channel,
    resolution = level,
    n.resolutions = handle$n.resolutions,
    pixel.size = pixel.size,
    origin = c(x = (x.range[1] - 1L), y = (y.range[1] - 1L)) * (if (is.na(pixel.size)) NA_real_ else pixel.size)
  )
}

#' Read a single channel plane of an Atera morphology image, optionally
#' cropped to \code{region} (a list with \code{x}/\code{y} elements, each a
#' length-2 micron range; \code{NULL}, the default, reads the whole plane).
#' A thin convenience wrapper around \code{.AteraMorphologyHandle} +
#' \code{.AteraReadMorphologyRegion} for the common one-shot case (see those
#' for details/return value).
#'
#' @keywords internal
#' @noRd
.AteraReadMorphologyImage <- function(data.dir, channel = "dapi", resolution = NULL, region = NULL) {
  handle <- .AteraMorphologyHandle(data.dir = data.dir, channel = channel)
  level <- resolution %||% handle$n.resolutions
  downsample <- 2 ^ (level - 1L)
  pixel.size <- if (is.na(handle$pixel.size.full)) NA_real_ else handle$pixel.size.full * downsample

  px.region <- .AteraMorphologyRegionToPixels(handle, level, region, pixel.size)
  .AteraReadMorphologyRegion(handle, resolution = level, x.range = px.region$x, y.range = px.region$y)
}

#' Turn a feature type name (eg \dQuote{Negative Control Probe}) into a
#' filesystem-safe directory name for use under \code{bpcells.dir}
#'
#' @keywords internal
#' @noRd
.AteraSanitizeName <- function(x) {
  gsub("[^A-Za-z0-9]+", "_", x)
}

#' Check that the \code{BPCells} package is installed; only called when a
#' \code{bpcells.dir} cache is actually requested
#'
#' @keywords internal
#' @noRd
.AteraCheckBPCellsDeps <- function() {
  if (!requireNamespace('BPCells', quietly = TRUE)) {
    stop("Reading/writing an Atera bpcells.dir cache requires the 'BPCells' package", call. = FALSE)
  }
}

#' Default, session-scoped \code{bpcells.dir} used when the caller hasn't
#' specified one: a subdirectory of \code{tempdir()} keyed off
#' \code{data.dir}, so repeated \code{LoadAtera()}/\code{ReadAtera()} calls
#' against the same bundle within a session reuse the same on-disk cache.
#'
#' @keywords internal
#' @noRd
.AteraDefaultBPCellsDir <- function(data.dir) {
  file.path(tempdir(), "atera_bpcells", .AteraSanitizeName(normalizePath(data.dir, mustWork = FALSE)))
}

#' Resolve the effective \code{bpcells.dir} to use, applying the
#' always-on-by-default policy: \code{NULL} (the parameter default) means
#' "cache on disk via BPCells if it's installed", falling back to in-memory
#' matrices with a one-time message if it isn't. Passing \code{FALSE}
#' explicitly opts out of BPCells entirely (plain in-memory matrices, no
#' dependency on the package). Any other value is used as-is, as an explicit
#' cache directory (and requires \code{BPCells} to be installed).
#'
#' @keywords internal
#' @noRd
.AteraResolveBPCellsDir <- function(bpcells.dir, data.dir) {
  if (isFALSE(bpcells.dir)) {
    return(NULL)
  }
  if (is.null(bpcells.dir)) {
    if (!requireNamespace('BPCells', quietly = TRUE)) {
      .AteraWarnOnceNoBPCells()
      return(NULL)
    }
    return(.AteraDefaultBPCellsDir(data.dir = data.dir))
  }
  return(bpcells.dir)
}

#' Emit a one-time (per session) message when falling back to in-memory
#' matrices because BPCells isn't installed
#'
#' @keywords internal
#' @noRd
.AteraWarnOnceNoBPCells <- function() {
  opt <- 'Seurat.atera.warned_no_bpcells'
  if (!isTRUE(getOption(opt))) {
    message(
      "Package 'BPCells' is not installed; loading Atera counts matrices ",
      "in-memory instead of on-disk. Install 'BPCells' for lazy, on-disk-",
      "backed loading (recommended for large bundles), or pass ",
      "`bpcells.dir = FALSE` to silence this message and always load in-",
      "memory."
    )
    options(structure(list(TRUE), names = opt))
  }
}

#' Split a contiguous feature row range \code{r1:r2} into smaller row
#' sub-ranges, each covering no more than \code{nnz.budget} non-zero matrix
#' entries (per the running element counts in \code{x.indptr}), so that
#' \code{.AteraReadCountsMatrix} can decode/build one bounded-size matrix
#' chunk at a time instead of materializing an entire (potentially
#' multi-billion-entry) feature type's worth of \code{X/data}/\code{X/indices}
#' in memory at once. A single row whose own nnz already exceeds
#' \code{nnz.budget} still gets its own (over-budget) batch, since rows can't
#' be split further.
#'
#' @return A list of length-2 integer vectors \code{c(start, end)}, each a
#' 1-based, inclusive row sub-range of \code{r1:r2}
#'
#' @keywords internal
#' @noRd
.AteraNnzRowBatches <- function(x.indptr, r1, r2, nnz.budget) {
  batches <- list()
  batch.start <- r1
  base <- x.indptr[r1]
  for (i in r1:r2) {
    if (x.indptr[i + 1L] - base >= nnz.budget) {
      batches[[length(batches) + 1L]] <- c(batch.start, i)
      batch.start <- i + 1L
      base <- x.indptr[i + 1L]
    }
  }
  if (batch.start <= r2) {
    batches[[length(batches) + 1L]] <- c(batch.start, r2)
  }
  batches
}

#' Decode a single feature (\code{X} row) sub-range \code{sr1:sr2} of the
#' counts matrix into an in-memory sparse matrix, named/dimensioned
#' consistently with the full matrix (all cells as columns)
#'
#' @keywords internal
#' @noRd
.AteraReadRowBlockMatrix <- function(con, zidx, x.indptr, sr1, sr2, var, cell.ids) {
  nnz.start <- x.indptr[sr1] + 1L
  nnz.end <- x.indptr[sr2 + 1L]
  if (nnz.end >= nnz.start) {
    d <- .AteraReadArray(con, zidx, "X/data", row.range = c(nnz.start, nnz.end))
    i <- .AteraReadArray(con, zidx, "X/indices", row.range = c(nnz.start, nnz.end))
  } else {
    d <- numeric(0)
    i <- numeric(0)
  }
  blk <- new(
    Class = "dgRMatrix",
    j = as.integer(i),
    p = as.integer(x.indptr[sr1:(sr2 + 1L)] - x.indptr[sr1]),
    x = as.double(d),
    Dim = c(sr2 - sr1 + 1L, length(cell.ids))
  )
  blk <- as(blk, "CsparseMatrix")
  rownames(blk) <- var$feature_name[sr1:sr2]
  colnames(blk) <- cell.ids
  blk
}

#' Build (and, when \code{bpcells.dir} is set, cache) one feature type's
#' counts matrix, in bounded-memory row batches (see
#' \code{.AteraNnzRowBatches}) rather than decoding the whole
#' \code{r1:r2}/\code{X/data} range at once -- the latter needs tens of GB of
#' RAM for a multi-billion-entry feature type (eg a whole-transcriptome
#' panel's "Gene Expression" rows), which doesn't fit on commodity machines.
#'
#' Each batch is decoded, (when caching) converted to a \code{BPCells}
#' \code{uint32_t} matrix and written to its own small temporary on-disk
#' directory, and then discarded from memory; batches are only ever
#' recombined (\code{rbind}) as cheap, lazy on-disk/\code{IterableMatrix}
#' references, and that combined reference is written to \code{cache.path}
#' (the real, permanent cache directory) \emph{before} the temporary
#' per-batch directories are cleaned up -- \code{BPCells::write_matrix_dir}'s
#' own chunked iteration streams through that final write without
#' re-materializing everything at once. Without \code{bpcells.dir}, batches
#' are instead recombined in-memory, so peak memory is still bounded by the
#' full matrix's size (an in-memory result is, by definition, not scalable
#' beyond available RAM).
#'
#' @return A \link[Matrix:dgCMatrix-class]{sparse matrix}, or (when
#' \code{bpcells.dir} is set) a \code{BPCells} \code{IterableMatrix} opened
#' from \code{cache.path}
#'
#' @keywords internal
#' @noRd
.AteraBuildFeatureTypeMatrix <- function(con, zidx, x.indptr, r1, r2, var, cell.ids, bpcells.dir, cache.path, nnz.budget) {
  row.batches <- .AteraNnzRowBatches(x.indptr, r1, r2, nnz.budget)
  use.bpcells <- !is.null(bpcells.dir)

  if (!use.bpcells) {
    parts <- lapply(row.batches, function(b) {
      .AteraReadRowBlockMatrix(con, zidx, x.indptr, b[1L], b[2L], var, cell.ids)
    })
    return(if (length(parts) == 1L) parts[[1L]] else do.call(rbind, parts))
  }

  .AteraCheckBPCellsDeps()
  staging.dir <- file.path(bpcells.dir, paste0(".tmp_", .AteraSanitizeName(basename(tempfile()))))
  dir.create(staging.dir, recursive = TRUE)
  on.exit(unlink(staging.dir, recursive = TRUE), add = TRUE)

  tmp.dirs <- vapply(seq_along(row.batches), function(k) {
    b <- row.batches[[k]]
    blk <- .AteraReadRowBlockMatrix(con, zidx, x.indptr, b[1L], b[2L], var, cell.ids)
    blk <- BPCells::convert_matrix_type(methods::as(blk, "IterableMatrix"), type = "uint32_t")
    tmp.dir <- file.path(staging.dir, k)
    BPCells::write_matrix_dir(mat = blk, dir = tmp.dir)
    tmp.dir
  }, character(1L))

  part.mats <- lapply(tmp.dirs, BPCells::open_matrix_dir)
  combined <- if (length(part.mats) == 1L) part.mats[[1L]] else do.call(rbind, part.mats)
  BPCells::write_matrix_dir(mat = combined, dir = cache.path)
  BPCells::open_matrix_dir(dir = cache.path)
}

#' Read the Atera counts matrix (\code{csc_cell_feature_matrix.zarr.zip}),
#' optionally restricted to a subset of feature types and optionally backed
#' by an on-disk \code{BPCells} cache.
#'
#' \code{var/feature_type} is row-contiguous in this file (each feature type
#' occupies one contiguous run of rows), so each requested feature type is
#' translated into a single contiguous row-range; that is in turn translated
#' (via a full, but trivially small, read of \code{X/indptr}) into nnz
#' element-ranges, and only those \code{X/data}/\code{X/indices} chunks are
#' decoded via \code{.AteraReadArray}'s \code{row.range} chunk-skipping. Each
#' feature type's row-range is further split into bounded-size row batches
#' (\code{.AteraBuildFeatureTypeMatrix}) so that building/caching even a
#' multi-billion-entry feature type stays within a few GB of peak memory.
#'
#' When \code{bpcells.dir} is set and every requested feature type already has
#' a cache directory under it, the zarr file is not read at all beyond the
#' (tiny) \code{var/feature_type} column needed to validate the requested
#' feature type names.
#'
#' @return A named list of matrices (one per feature type), each either a
#' \link[Matrix:dgCMatrix-class]{sparse matrix} or, when cached/written via
#' \code{bpcells.dir}, a \code{BPCells} \code{IterableMatrix}.
#'
#' @keywords internal
#' @noRd
.AteraReadCountsMatrix <- function(data.dir, feature.types = NULL, bpcells.dir = NULL, nnz.batch = 5e7) {
  zip.file <- file.path(data.dir, "csc_cell_feature_matrix.zarr.zip")
  zidx <- .AteraZipIndex(zip.file)
  con <- file(zip.file, "rb")
  on.exit(close(con), add = TRUE)

  ft.factor <- .AteraReadCategorical(con, zidx, "var/feature_type")
  ft.chr <- as.character(ft.factor)
  rl <- rle(ft.chr)
  ends <- cumsum(rl$lengths)
  starts <- ends - rl$lengths + 1L

  all.types <- rl$values
  if (is.null(feature.types)) {
    target.types <- all.types
  } else {
    unknown <- setdiff(feature.types, all.types)
    if (length(unknown) > 0L) {
      stop("Unknown Atera feature type(s): ", paste(unknown, collapse = ', '), call. = FALSE)
    }
    target.types <- feature.types
  }

  result <- vector(mode = 'list', length = length(target.types))
  names(result) <- target.types
  need.build <- character(0)
  for (ft in target.types) {
    cache.path <- if (!is.null(bpcells.dir)) file.path(bpcells.dir, .AteraSanitizeName(ft)) else NULL
    if (!is.null(cache.path) && dir.exists(cache.path)) {
      .AteraCheckBPCellsDeps()
      result[[ft]] <- BPCells::open_matrix_dir(dir = cache.path)
    } else {
      need.build <- c(need.build, ft)
    }
  }
  if (length(need.build) == 0L) {
    return(result)
  }

  var <- .AteraReadFlatGroup(con, zidx, "var", columns = c("feature_name", "feature_type"))
  obs <- .AteraReadFlatGroup(con, zidx, "obs", columns = "cell_id")
  cell.ids <- .AteraFormatId(.AteraDecodePackedId(obs$cell_id))

  x.indptr <- .AteraReadArray(con, zidx, "X/indptr")

  if (!is.null(bpcells.dir)) {
    dir.create(bpcells.dir, showWarnings = FALSE, recursive = TRUE)
  }
  for (ft in need.build) {
    ft.i <- match(ft, rl$values)
    cache.path <- if (!is.null(bpcells.dir)) file.path(bpcells.dir, .AteraSanitizeName(ft)) else NULL
    result[[ft]] <- .AteraBuildFeatureTypeMatrix(
      con = con, zidx = zidx, x.indptr = x.indptr,
      r1 = starts[ft.i], r2 = ends[ft.i],
      var = var, cell.ids = cell.ids,
      bpcells.dir = bpcells.dir, cache.path = cache.path, nnz.budget = nnz.batch
    )
  }
  return(result)
}

#' Load Atera spatial data
#'
#' Read the output of \href{https://www.10xgenomics.com}{10x Genomics} Atera,
#' a next-generation in situ sequencing platform. All Atera outputs (aside
#' from images) are stored in AnnData zarr (storage format v2) files, zipped
#' with \code{Stored} (uncompressed) zip entries; reading them requires the
#' \code{blosc} and \code{jsonlite} packages. Zarr chunks are read directly
#' out of the zip archive by byte-range seeks with no disk extraction.
#'
#' @param data.dir Directory containing all Atera output files with default
#' filenames
#' @param outs Types of outputs to read; choose one or more of:
#' \itemize{
#'  \item \dQuote{matrix}: the counts matrix
#'  \item \dQuote{centroids}: cell centroids in micron coordinate space
#'  \item \dQuote{segmentations}: cell segmentation boundary polygons in
#'  micron coordinate space
#'  \item \dQuote{nucleus_segmentations}: nucleus segmentation boundary
#'  polygons in micron coordinate space. Cells may have zero, one, or more
#'  than one nucleus polygon
#'  \item \dQuote{microns}: transcript molecule coordinates. This is
#'  optional as transcripts files can be very large.
#'  \item \dQuote{morphology}: a single channel plane of a
#'  \dQuote{morphology_2d} image (eg the DAPI stain), by default at the
#'  lowest-resolution pyramid level, for use as a background/overview image.
#'  Requires the \code{RBioFormats} package, since these images are
#'  pyramidal, JPEG2000-compressed OME-TIFFs that the base \code{tiff}
#'  package cannot reliably decode.
#' }
#' @param mols.qv.threshold Remove transcript molecules with a calibrated
#' Q-score less than this threshold when reading \dQuote{microns}. Set to
#' \code{NULL} to disable filtering.
#' @param morphology.channel Channel to read when \dQuote{morphology} is in
#' \code{outs}: either a channel name (matched case-insensitively as a
#' substring, eg \dQuote{dapi}) or a 0-based integer channel index (matching
#' the \dQuote{chNNNN} prefix of the image filename).
#' @param morphology.resolution Pyramid resolution level to read when
#' \dQuote{morphology} is in \code{outs}, where \code{1} is full resolution
#' and each subsequent level halves both dimensions. Defaults to \code{NULL},
#' which reads the lowest-resolution (smallest, fastest) level available --
#' normally all that's needed for an overview plot.
#' @param morphology.region Optional list with \code{x}/\code{y} elements
#' (each a length-2 micron range, in the same coordinate space as
#' \dQuote{centroids}/\dQuote{microns}/segmentations), used to crop
#' \dQuote{morphology} to a rectangular window instead of reading the whole
#' plane. Only the on-disk tiles overlapping that window are decoded, so
#' this is the recommended way to read a region at \code{morphology.resolution
#' = 1} (full resolution): reading a full-resolution plane in its entirety
#' can exceed available memory for these whole-slide images, since they can
#' be tens of thousands of pixels per side. Defaults to \code{NULL}, which
#' reads the whole plane.
#' @param genes Optional character vector of gene names to restrict
#' \dQuote{microns} to. When set, only the chunks spanning each requested
#' gene's (gene-sorted, contiguous) row range are read/decompressed instead
#' of the full transcript table; this avoids materializing all rows when
#' only a handful of genes are needed. Ignored when \code{"microns"} is not
#' in \code{outs}.
#' @param feature.types Optional character vector of feature types (eg
#' \dQuote{Gene Expression}) to restrict \dQuote{matrix} to. When set, only
#' the (contiguous) rows for the requested feature type(s) are
#' read/decompressed instead of the full counts matrix. Defaults to
#' \code{NULL}, which reads all feature types (current/default behavior).
#' Ignored when \code{"matrix"} is not in \code{outs}.
#' @param bpcells.dir Path to a directory used to cache the counts matrix on
#' disk, one subdirectory per feature type, via
#' \link[BPCells:write_matrix_dir]{BPCells}. When set, a feature type already
#' cached under this directory is opened directly with
#' \link[BPCells:open_matrix_dir]{BPCells::open_matrix_dir} (no zarr read at
#' all); a feature type not yet cached is read as usual and then written to
#' this directory for reuse by later calls. Defaults to \code{NULL}, which
#' auto-selects a session-scoped cache directory under \code{tempdir()} (so
#' matrices are BPCells-backed, on-disk, and lazy by default, even for small
#' bundles) if the \code{BPCells} package is installed, or falls back to
#' in-memory sparse matrices (with a one-time message) if it isn't. Pass
#' \code{FALSE} to explicitly opt out and always return in-memory sparse
#' matrices, without requiring \code{BPCells} at all. Ignored when
#' \code{"matrix"} is not in \code{outs}.
#'
#' @return \code{ReadAtera}: A list with some combination of the following
#' values:
#' \itemize{
#'  \item \dQuote{\code{matrix}}: a named list of
#'  \link[Matrix:dgCMatrix-class]{sparse matrices} with expression data, one
#'  per feature type (eg \dQuote{Gene Expression}); cells are columns and
#'  features are rows
#'  \item \dQuote{\code{centroids}}: a data frame with cell centroid
#'  coordinates in three columns: \dQuote{x}, \dQuote{y}, and \dQuote{cell}
#'  \item \dQuote{\code{segmentations}}/\dQuote{\code{nucleus_segmentations}}:
#'  a data frame with one row per polygon vertex, in three columns:
#'  \dQuote{cell}, \dQuote{x}, and \dQuote{y}
#'  \item \dQuote{\code{microns}}: a data frame with transcript coordinates
#'  in three columns: \dQuote{x}, \dQuote{y}, and \dQuote{gene}
#'  \item \dQuote{\code{morphology}}: a list with elements \code{image} (a
#'  numeric matrix, indexed \verb{[x, y]}, in the same top-left-origin pixel
#'  convention as \dQuote{centroids}/\dQuote{microns}' micron coordinates --
#'  no axis flip is needed), \code{channel}, \code{resolution},
#'  \code{n.resolutions}, \code{pixel.size} (microns per pixel of
#'  \code{image}), and \code{origin} (the micron coordinates,
#'  \code{c(x=, y=)}, of \code{image}'s \verb{[1, 1]} pixel -- \code{c(x=0,
#'  y=0)} unless \code{morphology.region} crops out the top-left corner of
#'  the full plane); see
#'  \code{morphology.channel}/\code{morphology.resolution}/\code{morphology.region}
#' }
#'
#' @export
#' @concept preprocessing
#'
ReadAtera <- function(
  data.dir,
  outs = c("matrix", "centroids"),
  mols.qv.threshold = 20,
  genes = NULL,
  feature.types = NULL,
  bpcells.dir = NULL,
  morphology.channel = "dapi",
  morphology.resolution = NULL,
  morphology.region = NULL
) {
  outs <- match.arg(
    arg = outs,
    choices = c("matrix", "centroids", "segmentations", "nucleus_segmentations", "microns", "morphology"),
    several.ok = TRUE
  )

  .AteraCheckDeps()

  data <- sapply(outs, function(otype) {
    switch(
      EXPR = otype,
      'matrix' = {
        pmtx <- progressor()
        pmtx(message = 'Reading counts matrix', class = 'sticky', amount = 0)

        bpcells.dir <- .AteraResolveBPCellsDir(bpcells.dir = bpcells.dir, data.dir = data.dir)
        if (!is.null(bpcells.dir)) {
          .AteraCheckBPCellsDeps()
        }
        result <- .AteraReadCountsMatrix(
          data.dir = data.dir,
          feature.types = feature.types,
          bpcells.dir = bpcells.dir
        )

        pmtx(type = "finish")

        result
      },
      'centroids' = {
        pcents <- progressor()
        pcents(message = 'Loading cell centroids', class = 'sticky', amount = 0)

        zip.file <- file.path(data.dir, "cells.zarr.zip")
        zidx <- .AteraZipIndex(zip.file)
        con <- file(zip.file, "rb")
        on.exit(close(con), add = TRUE)

        cell_id <- .AteraDecodePackedId(.AteraReadArray(con, zidx, "cell_id"))
        cell_summary <- .AteraReadArray(con, zidx, "cell_summary")

        pcents(type = 'finish')

        data.frame(
          x = cell_summary[, 1],
          y = cell_summary[, 2],
          cell = .AteraFormatId(cell_id)
        )
      },
      'segmentations' = ,
      'nucleus_segmentations' = {
        pseg <- progressor()
        pseg(
          message = if (otype == 'segmentations') 'Loading cell segmentations' else 'Loading nucleus segmentations',
          class = 'sticky',
          amount = 0
        )

        zip.file <- file.path(data.dir, "cells.zarr.zip")
        zidx <- .AteraZipIndex(zip.file)
        con <- file(zip.file, "rb")
        on.exit(close(con), add = TRUE)

        cell_id <- .AteraFormatId(.AteraDecodePackedId(.AteraReadArray(con, zidx, "cell_id")))
        set.idx <- if (otype == 'segmentations') 1L else 0L
        df <- .AteraReadPolygonSet(con, zidx, set.idx, cell_id)

        pseg(type = 'finish')

        df
      },
      'microns' = {
        pmicrons <- progressor()
        pmicrons(message = "Loading transcript coordinates", class = 'sticky', amount = 0)

        handle <- .AteraMoleculesHandle(data.dir = data.dir, mols.qv.threshold = mols.qv.threshold)
        df <- .AteraFetchMolecules(handle = handle, genes = genes)

        pmicrons(type = 'finish')

        df
      },
      'morphology' = {
        pmorph <- progressor()
        pmorph(message = 'Loading morphology image', class = 'sticky', amount = 0)

        result <- .AteraReadMorphologyImage(
          data.dir = data.dir,
          channel = morphology.channel,
          resolution = morphology.resolution,
          region = morphology.region
        )

        pmorph(type = 'finish')

        result
      },
      stop("Unknown Atera input type: ", otype)
    )
  }, simplify = FALSE, USE.NAMES = TRUE)

  metadata <- file.path(data.dir, "experiment.spatial")
  if (file.exists(metadata) && requireNamespace("jsonlite", quietly = TRUE)) {
    meta <- jsonlite::read_json(metadata)
    data$metadata <- meta[
      intersect(
        names(meta),
        c(
          'run_name', 'slide_name', 'region_name', 'chemistry_version',
          'num_cells', 'num_transcripts'
        )
      )
    ]
    if (!is.null(meta$panel_a)) {
      data$metadata$panel_name <- meta$panel_a$panel_name
    }
  }
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
  if (isFALSE(x = requireNamespace('jsonlite', quietly = TRUE))) {
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
  if (isFALSE(x = requireNamespace('data.table', quietly = TRUE))) {
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
      # drop = FALSE so that a chunk holding a single feature stays a matrix
      # and t() keeps one row per cell
      pp[["marks"]] <- as.data.frame(x = t(x = data[x, , drop = FALSE]))
      chunk <- markvario(X = pp, normalise = TRUE, ...)
      # markvario() hands back a bare 'fv' for a single mark and a named list
      # of them for several. Keep every chunk a list, so that the unlist()
      # below contributes one entry per feature rather than one per column of
      # an 'fv' whenever a chunk holds a single feature.
      if (inherits(x = chunk, what = 'fv')) {
        chunk <- list(chunk)
      }
      return(chunk)
    })
    mv <- unlist(x = mv, recursive = FALSE)
  } else {
    pp[["marks"]] <- as.data.frame(x = t(x = data))
    mv <- markvario(X = pp, normalise = TRUE, ...)
    if (inherits(x = mv, what = 'fv')) {
      mv <- list(mv)
    }
  }
  names(x = mv) <- rownames(x = data)
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
  Rfast2.installed <- requireNamespace('Rfast2', quietly = TRUE)
  if (isTRUE(x = Rfast2.installed)) {
    MyMoran <- Rfast2::moranI
  } else if (isFALSE(x = requireNamespace('ape', quietly = TRUE))) {
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
#' counts = as.matrix(x = GetAssayData(object = pbmc_small, assay = "RNA", layer = "counts"))
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

#' SCTransform: Regularized NB regression for UMI count normalization
#'
#' Perform a variance‐stabilizing transformation on UMI counts using
#' \code{sctransform::vst} (https://github.com/satijalab/sctransform). This
#' replaces the \code{NormalizeData} → \code{FindVariableFeatures} →
#' \code{ScaleData} workflow by fitting a regularized negative binomial model
#' per gene and returning:
#'
#' - A new assay (default name “SCT”), in which:
#'   - \code{counts}: depth‐corrected UMI counts (as if each cell had uniform
#'     sequencing depth; controlled by \code{do.correct.umi}).
#'   - \code{data}: \code{log1p} of corrected counts.
#'   - \code{scale.data}: Pearson residuals from the fitted NB model (optionally
#'     centered and/or scaled).
#'   - \code{misc}: intermediate outputs from \code{sctransform::vst}.
#'
#' When multiple \code{counts} layers exist (e.g. after \code{split()}),
#' each layer is modeled independently. A consensus variable‐feature set is
#' then defined by ranking features by how often they’re called “variable”
#' across different layers (ties broken by median rank).
#'
#' By default, \code{sctransform::vst} will drop features expressed in fewer
#' than five cells. In the multi-layer case, this can lead to consenus
#' variable-features being excluded from the output's \code{scale.data} when
#' a feature is "variable" across many layers but sparsely expressed in at
#' least one.
#'
#' @param object A Seurat object or UMI count matrix.
#' @param cell.attr Optional metadata frame (cells × attributes).
#' @param reference.SCT.model Pre‐fitted SCT model (supports only \code{log_umi}
#'   as latent variable). If provided, computes residuals via that model. When
#'   \code{residual.features} is NULL, uses the model’s top
#'   \code{variable.features.n}; otherwise, sets the assay’s variable features
#'   to \code{residual.features}.
#' @param do.correct.umi Logical; if TRUE (default), stores corrected UMIs in
#'   \code{counts}.
#' @param ncells Integer; number of cells to subsample when fitting NB
#'   regression (default: 5000).
#' @param residual.features Character vector of genes to compute residuals for.
#'   Default NULL (all genes). If set, these become the assay’s variable
#'   features.
#' @param variable.features.n Integer; when \code{residual.features} is NULL,
#'   select this many top features by residual variance (default: 3000).
#' @param variable.features.rv.th Numeric; if \code{variable.features.n} is NULL,
#'   select features exceeding this residual‐variance threshold (default: 1.3).
#' @param vars.to.regress Character vector of metadata columns (e.g.
#'   \code{percent.mito}) to regress out in a second, non‐regularized model.
#' @param latent.data Numeric matrix (cells × latent covariates) to regress out.
#' @param do.scale Logical; if TRUE, scale residuals to unit variance
#'   (default: FALSE).
#' @param do.center Logical; if TRUE, center residuals to mean zero
#'   (default: TRUE).
#' @param clip.range Numeric vector of length 2; range to clip residuals
#'   (default \code{c(-sqrt(n/30), sqrt(n/30))}, with n = number of cells).
#' @param vst.flavor Character; if \code{"v2"}, uses \code{method = "glmGamPoi_offset"},
#'   \code{n_cells = 2000}, and \code{exclude_poisson = TRUE} to fit \eqn{\theta} and
#'   intercept only.
#' @param conserve.memory Logical; if TRUE, never builds the full residual
#'   matrix (slower but memory‐efficient; forces \code{return.only.var.genes=TRUE};
#'   default: FALSE).
#' @param return.only.var.genes Logical; if TRUE (default), \code{scale.data}
#'   is subset to variable features only.
#' @param defer.residual.matrix Logical; if TRUE, skip materializing the
#'   Pearson residual matrix in this call. (For internal use by the
#'   v5 multi-layer \code{SCTransform} workflow, which computes the final
#'   merged \code{scale.data} matrix after processing all layers. Not used for
#'   the standard v3/single-matrix workflow.)
#' @param seed.use Integer; random seed for reproducibility (default: 1448145).
#'   Set to NULL to skip setting a seed.
#' @param verbose Logical; whether to print progress messages (default: TRUE).
#' @param ... Additional arguments passed to \code{sctransform::vst}.
#'
#' @return A Seurat object with a new \code{SCT} assay containing:
#' \code{counts} (corrected UMIs), \code{data} (log1p counts), and
#' \code{scale.data} (Pearson residuals), plus \code{misc} for intermediate
#' \code{vst} outputs.
#'
#' @importFrom stats setNames
#' @importFrom Matrix colSums
#' @importFrom SeuratObject as.sparse
#' @importFrom sctransform vst get_residual_var get_residuals correct_counts
#'
#' @seealso \code{\link[sctransform]{vst}},
#'   \code{\link[sctransform]{get_residuals}},
#'   \code{\link[sctransform]{correct_counts}}
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
  latent.data = NULL,
  do.scale = FALSE,
  do.center = TRUE,
  clip.range = c(-sqrt(x = ncol(x = umi) / 30), sqrt(x = ncol(x = umi) / 30)),
  vst.flavor = 'v2',
  conserve.memory = FALSE,
  return.only.var.genes = TRUE,
  defer.residual.matrix = FALSE,
  seed.use = 1448145,
  verbose = TRUE,
  ...
) {
  if (!is.null(x = seed.use)) {
    set.seed(seed = seed.use)
  }
  vst.args <- list(...)
  object <- as.sparse(x = object)
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

  if (!is.null(x = vst.flavor) && !vst.flavor %in% c("v1", "v2")){
    stop("vst.flavor can be 'v1' or 'v2'. Default is 'v2'")
  }
  if (!is.null(x = vst.flavor) && vst.flavor == "v1"){
    vst.flavor <- NULL
  }

  vst.args[['vst.flavor']] <- vst.flavor
  vst.args[['umi']] <- umi
  vst.args[['cell_attr']] <- cell.attr
  vst.args[['verbosity']] <- as.numeric(x = verbose) * 1
  vst.args[['return_cell_attr']] <- TRUE
  vst.args[['return_gene_attr']] <- TRUE
  vst.args[['return_corrected_umi']] <- do.correct.umi
  vst.args[['n_cells']] <- min(ncells, ncol(x = umi))
  residual.type <- vst.args[['residual_type']] %||% 'pearson'
  res.clip.range <- vst.args[['res_clip_range']] %||% c(-sqrt(x = ncol(x = umi)), sqrt(x = ncol(x = umi)))
  fast.sct.default <- (
    identical(x = residual.type, y = 'pearson') &&
    !any(c("batch_var", "latent_var", "latent_var_nonreg") %in% names(x = vst.args)) &&
    !any(vst.args[['min_variance']] %in% c('model_mean', 'model_median'))
  )
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
    'default' = {
      if (isTRUE(x = fast.sct.default)) {
        vst.args[['return_corrected_umi']] <- FALSE
        vst.args[['residual_type']] <- 'none'
      } else {
        vst.args[['return_corrected_umi']] <- do.correct.umi
      }
      vst.out <- do.call(what = 'vst', args = vst.args)
      vst.out
    },
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
      if (do.correct.umi && !is.null(x = vst.out$umi_corrected)) {
        dimnames(x = vst.out$umi_corrected) <- dimnames(x = umi)
      }
      vst.out
    })

  # get residuals
  vst.out <- switch(
    EXPR = sct.method,
    'default' = {
      if (!isTRUE(x = fast.sct.default)) {
        feature.variance <- vst.out$gene_attr[, "residual_variance"]
        names(x = feature.variance) <- rownames(x = vst.out$gene_attr)
        feature.variance <- sort(x = feature.variance, decreasing = TRUE)
        feature.idx <- if (is.null(x = variable.features.n)) {
          feature.variance >= variable.features.rv.th
        } else {
          seq_len(length.out = min(variable.features.n, length(x = feature.variance)))
        }
        top.features <- names(x = feature.variance)[feature.idx]
        if (isTRUE(x = return.only.var.genes)) {
          scale.data.features <- intersect(x = top.features, y = rownames(x = vst.out$y))
          vst.out$y <- vst.out$y[scale.data.features, , drop = FALSE]
        }
        vst.out
      } else {
        model.pars <- vst.out$model_pars_fit
        genes <- rownames(x = model.pars)
        if (!identical(x = genes, y = rownames(x = umi))) {
          umi <- umi[genes, , drop = FALSE]
        }
        min.variance <- vst.out$arguments$min_variance
        min.var <- if (identical(x = min.variance, y = "umi_median")) {
          (median(umi@x) / 5) ^ 2
        } else {
          min.variance
        }
        # Persist the resolved numeric min_var in the model. Downstream residual
        # recomputation (FetchResiduals / GetResidual) reads arguments$min_variance
        # and only recomputes (median(nonzeros)/5)^2 when it is the string "umi_median".
        # That recompute is order/subset dependent, so it can diverge from the
        # value used here for scale.data. Storing the resolved value makes later
        # residuals deterministic and consistent with scale.data.
        vst.out$arguments$min_variance <- min.var
        # should be set already by the vst call but just fixing in case its null
        res.clip.range <- vst.out$arguments$res_clip_range %||%
          c(-sqrt(x = ncol(x = umi)), sqrt(x = ncol(x = umi)))

        # Compute residual statistics and corrected UMI counts
        # Note: does not compute residual matrix yet (saves a lot of memory)
        # Just computes the residual variance for each gene and (if asked for) corrected UMI counts
        stats <- SCTResidualStatsAndCorrected(
          x = umi@x,
          i = umi@i,
          p = umi@p,
          rows = nrow(x = umi),
          cols = ncol(x = umi),
          theta = model.pars[, "theta"],
          intercept = model.pars[, "(Intercept)"],
          slope = model.pars[, "log_umi"],
          log_umi = vst.out$cell_attr[colnames(x = umi), "log_umi"],
          target_log_umi = median(vst.out$cell_attr[, "log_umi"]),
          min_var = min.var,
          residual_clip_min = min(res.clip.range),
          residual_clip_max = max(res.clip.range),
          n_threads = getThreads(verbose = FALSE),
          compute_corrected = do.correct.umi
        )
        vst.out$gene_attr[genes, "residual_mean"] <- stats$residual_mean
        vst.out$gene_attr[genes, "residual_variance"] <- stats$residual_variance
    
        # Determine variable features
        feature.variance <- vst.out$gene_attr[, "residual_variance"]
        names(x = feature.variance) <- rownames(x = vst.out$gene_attr)
        feature.variance <- sort(x = feature.variance, decreasing = TRUE)
        feature.idx <- if (is.null(x = variable.features.n)) {
          feature.variance >= variable.features.rv.th
        } else {
          seq_len(length.out = min(variable.features.n, length(x = feature.variance)))
        }
        top.features <- names(x = feature.variance)[feature.idx]

        # Store corrected UMI counts if requested, if not just restore original counts matrix
        if (do.correct.umi) {
          vst.out$umi_corrected <- stats$corrected
          dimnames(x = vst.out$umi_corrected) <- dimnames(x = umi)
        } else {
          vst.out$umi_corrected <- umi
        }

        # Compute matrix of Pearson residuals for features to be included in scale.data
        scale.data.features <- if (return.only.var.genes) {
          top.features
        } else {
          genes
        }
        if (isTRUE(x = defer.residual.matrix)) {
          vst.out$y <- matrix(
            data = numeric(length = 0L),
            nrow = 0L,
            ncol = ncol(x = umi),
            dimnames = list(character(length = 0L), colnames(x = umi))
          )
        } else {
          vst.out$y <- SCTPearsonResidualMatrix(
            x = umi@x,
            i = umi@i,
            p = umi@p,
            rows = nrow(x = umi),
            cols = ncol(x = umi),
            theta = model.pars[, "theta"],
            intercept = model.pars[, "(Intercept)"],
            slope = model.pars[, "log_umi"],
            log_umi = vst.out$cell_attr[colnames(x = umi), "log_umi"],
            feature_index = as.integer(x = match(x = scale.data.features, table = genes) - 1L),
            min_var = min.var,
            clip_min = min(clip.range),
            clip_max = max(clip.range),
            do_center = do.center,
            n_threads = getThreads(verbose = FALSE),
            display_progress = verbose
          )
          dimnames(x = vst.out$y) <- list(scale.data.features, colnames(x = umi))
        }
      
      vst.out
      }
    },
    'reference.model' = {
      feature.variance <- vst.out$gene_attr[, "residual_variance"]
      names(x = feature.variance) <- rownames(x = vst.out$gene_attr)

      feature.variance <- sort(x = feature.variance, decreasing = TRUE)

      feature.idx <- if (is.null(x = variable.features.n)) {
        feature.variance >= variable.features.rv.th
      } else {
        seq_len(length.out = min(variable.features.n, length(x = feature.variance)))
      }
      top.features <- names(x = feature.variance)[feature.idx]
      if (is.null(x = residual.features)) {
        residual.features <- top.features
      }

      residual.features <- Reduce(
        f = intersect,
        x = list(residual.features, rownames(x = umi), rownames(x = vst.out$model_pars_fit))
      )
      sub <- umi[residual.features, , drop = FALSE]
      min.variance <- vst.out$arguments$min_variance
      # Reproduce sctransform::get_residuals() -
      # res_clip_range = +/- sqrt(ncol(sub)), the scalar variance floor from the
      # model, and do_center = FALSE (reference centering by the reference
      # residual_mean is applied by the sweep below, not by the kernel).
      # "model_mean"/"model_median" use a per-gene variance floor -> fall back to
      # get_residuals() for exact behavior.
      if (min.variance %in% c("model_mean", "model_median")) {
        residual.feature.mat <- get_residuals(
          vst_out = vst.out,
          umi = sub,
          verbosity = as.numeric(x = verbose) * 2
        )
      } else {
        model.pars <- vst.out$model_pars_fit[residual.features, , drop = FALSE]
        min.var <- if (identical(x = min.variance, y = "umi_median")) {
          (median(x = sub@x) / 5) ^ 2
        } else {
          min.variance
        }
        res.clip.range <- c(-sqrt(x = ncol(x = sub)), sqrt(x = ncol(x = sub)))
        residual.feature.mat <- SCTPearsonResidualMatrix(
          x = sub@x,
          i = sub@i,
          p = sub@p,
          rows = nrow(x = sub),
          cols = ncol(x = sub),
          theta = model.pars[, "theta"],
          intercept = model.pars[, "(Intercept)"],
          slope = model.pars[, "log_umi"],
          log_umi = vst.out$cell_attr[colnames(x = sub), "log_umi"],
          feature_index = as.integer(x = seq_len(length.out = nrow(x = sub)) - 1L),
          min_var = min.var,
          clip_min = min(res.clip.range),
          clip_max = max(res.clip.range),
          do_center = FALSE,
          n_threads = getThreads(verbose = FALSE),
          display_progress = verbose
        )
        dimnames(x = residual.feature.mat) <- dimnames(x = sub)
      }
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
      feature.variance <- vst.out$gene_attr[, "residual_variance"]
      names(x = feature.variance) <- rownames(x = vst.out$gene_attr)
      feature.variance <- sort(x = feature.variance, decreasing = TRUE)
      feature.idx <- if (is.null(x = variable.features.n)) {
        feature.variance >= variable.features.rv.th
      } else {
        seq_len(length.out = min(variable.features.n, length(x = feature.variance)))
      }
      top.features <- names(x = feature.variance)[feature.idx]
      vst.out$y <- get_residuals(
        vst_out = vst.out,
        umi = umi[top.features, , drop = FALSE],
        residual_type = residual.type,
        res_clip_range = res.clip.range,
        verbosity = as.numeric(x = verbose) * 2
      )
      vst.out$gene_attr$residual_mean <- NA_real_
      vst.out$gene_attr[top.features, "residual_mean"] <- rowMeans2(x = vst.out$y)
      if (do.correct.umi && residual.type == 'pearson') {
        vst.out$umi_corrected <- correct_counts(
          x = vst.out,
          umi = umi,
          verbosity = as.numeric(x = verbose) * 1
        )
      }
      vst.out
    }
   )
  residuals.preprocessed <- identical(x = sct.method, y = "default") && isTRUE(x = fast.sct.default)
  if (!isTRUE(x = residuals.preprocessed)) {
    scale.data <- vst.out$y
    # clip the residuals
    scale.data[scale.data < clip.range[1]] <- clip.range[1]
    scale.data[scale.data > clip.range[2]] <- clip.range[2]
    # 2nd regression
    vst.out$y <- ScaleData(
      scale.data,
      features = NULL,
      vars.to.regress = vars.to.regress,
      latent.data = latent.data,
      model.use = 'linear',
      use.umi = FALSE,
      do.scale = do.scale,
      do.center = do.center,
      scale.max = Inf,
      block.size = 750,
      min.cells.to.block = 3000,
      verbose = verbose
    )
  }
  
  # User may (not common) want to regress out additional variables after SCTransform
  # when the residual matrix has already been centered and clipped.
  if (isTRUE(x = residuals.preprocessed) && (!is.null(x = vars.to.regress) || isTRUE(x = do.scale))) {
    if (is.null(x = rownames(x = vst.out$y)) && nrow(x = vst.out$y) == nrow(x = vst.out$gene_attr)) {
      rownames(x = vst.out$y) <- rownames(x = vst.out$gene_attr)
    }
    if (is.null(x = colnames(x = vst.out$y)) && ncol(x = vst.out$y) == nrow(x = vst.out$cell_attr)) {
      colnames(x = vst.out$y) <- rownames(x = vst.out$cell_attr)
    }
    vst.out$y <- ScaleData(
      vst.out$y,
      features = NULL,
      vars.to.regress = vars.to.regress,
      latent.data = latent.data,
      model.use = 'linear',
      use.umi = FALSE,
      do.scale = do.scale,
      do.center = do.center,
      scale.max = Inf,
      block.size = 750,
      min.cells.to.block = 3000,
      verbose = verbose
    )
  }

  vst.out$variable_features <- residual.features %||% top.features
  if (
    !is.null(x = vst.out$umi_corrected) &&
      identical(x = dim(x = vst.out$umi_corrected), y = dim(x = umi)) &&
      (is.null(x = rownames(x = vst.out$umi_corrected)) || is.null(x = colnames(x = vst.out$umi_corrected)))
  ) {
    dimnames(x = vst.out$umi_corrected) <- dimnames(x = umi)
  }
  if (do.correct.umi && is.null(x = vst.out$umi_corrected)) {
    vst.out$residual_type <- vst.out$residual_type %||% "none"
  }
  if (!do.correct.umi) {
    vst.out$umi_corrected <- umi
  }
  if (sct.method %in% c('default', 'conserve.memory')) {
    # Store the residual type used for output, including deferred residuals.
    vst.out$arguments$residual_type <- residual.type
  }
  min_var <- vst.out$arguments$min_variance
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
    latent.data = NULL,
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
    rlang::warn(
      "A reference SCT model was provided, therefore counts are not corrected (regardless of do.correct.umi)",
      .frequency = "once",
      .frequency_id = "SCTransform-reference-SCTmodel-correct-counts"
    )
  }

  umi <- GetAssayData(object = object, layer = 'counts')
  vst.out <- SCTransform(object = umi,
                         cell.attr = cell.attr,
                         reference.SCT.model = reference.SCT.model,
                         do.correct.umi = do.correct.umi,
                         ncells = ncells,
                         residual.features = residual.features,
                         variable.features.n = variable.features.n,
                         variable.features.rv.th = variable.features.rv.th,
                         vars.to.regress = vars.to.regress,
                         latent.data = latent.data,
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
    layer = 'data',
    new.data = log1p(x = GetAssayData(object = assay.out, layer = 'counts'))
  )
  scale.data <- vst.out$y
  assay.out <- SetAssayData(
    object = assay.out,
    layer = 'scale.data',
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
    assay = "RNA",
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
  if (!is.null(x = seed.use)) {
    set.seed(seed = seed.use)
  }
  if (any(vars.to.regress %in% colnames(x = object[[]]))) {
    vars.to.regress.subset <- vars.to.regress[vars.to.regress %in% colnames(x = object[[]])]
    latent.data <- object[[vars.to.regress.subset]]
  } else {
    latent.data <- NULL
  }
  assay <- assay %||% DefaultAssay(object = object)
  if (assay == "SCT") {
    # if re-running SCTransform, use the RNA assay
    assay <- "RNA"
    warning("Running SCTransform on the RNA assay while default assay is SCT.")
  }

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
                            latent.data = latent.data,
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
  
  # Extract all SCT models stored in assay
  sct_models <- slot(object = assay.data, name = "SCTModel.list")
  
  # Update umi.assay field for every SCT model 
  slot(object = assay.data, name = "SCTModel.list") <- lapply(sct_models, function(model) {
    slot(model, name = "umi.assay") <- assay
    model
  })

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
#'   \item \dQuote{\code{vst}}:  First, fits a line to the relationship of
#'     log(variance) and log(mean) using local polynomial regression (loess).
#'     Then standardizes the feature values using the observed mean and
#'     expected variance (given by the fitted line). Feature variance is then
#'     calculated on the standardized values
#'     after clipping to a maximum (see clip.max parameter).
#'   \item \dQuote{\code{mean.var.plot}} (mvp): First, uses a function to
#'     calculate average expression (mean.function) and dispersion
#'     (dispersion.function) for each feature. Next, divides features into
#'     \code{num.bin} (default 20) bins based on their average expression,
#'     and calculates z-scores for dispersion within each bin. The purpose of
#'     this is to identify variable features while controlling for the
#'     strong relationship between variability and average expression
#'   \item \dQuote{\code{dispersion}} (disp): selects the genes with the
#'     highest dispersion values
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
#'   \item \dQuote{\code{equal_width}}: each bin is of equal width along the
#'     x-axis (default)
#'   \item \dQuote{\code{equal_frequency}}: each bin contains an equal number
#'     of features (can increase statistical power to detect overdispersed
#'     features at high expression values, at the cost of reduced resolution
#'     along the x-axis)
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
  nfeatures = 2000,
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
    hvf.info <- .FindVariableFeaturesVSTInfo(
      object = object,
      loess.span = loess.span,
      clip.max = clip.max,
      nselect = nfeatures,
      verbose = verbose
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
    data <- GetAssayData(object = object, layer = "counts")
    # if (ncol(x = data) < 1 || nrow(x = data) < 1) {
    if (IsMatrixEmpty(x = data)) {
      warning("selection.method set to 'vst' but count slot is empty; will use data slot instead")
      data <- GetAssayData(object = object, layer = "data")
    }
  } else {
    data <- GetAssayData(object = object, layer = "data")
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
  object[[names(x = hvf.info)]] <- hvf.info
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
  object[[vf.name]] <- rownames(x = object[[]]) %in% top.features
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
  VariableFeatures(object) <- VariableFeatures(object, nfeatures = nfeatures)
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
  selection.method <- match.arg(arg = selection.method)
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

#' @param layer The layer in the specified assay to pull data from.
#' @param slot Deprecated, use `layer`.
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
  layer = "scale.data",
  slot = deprecated(),
  spatial.location,
  selection.method = c('markvariogram', 'moransi'),
  features = NULL,
  r.metric = 5,
  x.cuts = NULL,
  y.cuts = NULL,
  nfeatures = 2000,
  verbose = TRUE,
  ...
) {
  if (is_present(slot)) {
    deprecate_soft(
      when = '5.3.0',
      what = 'FindSpatiallyVariableFeatures(slot = )',
      with = 'FindSpatiallyVariableFeatures(layer = )'
    )
    layer <- slot %||% layer
  }
  features <- features %||% Features(object, layer = layer)
  selection.method <- match.arg(selection.method)
  if (selection.method == "markvariogram" && "markvariogram" %in% names(x = Misc(object = object))) {
    features.computed <- names(x = Misc(object = object, slot = "markvariogram"))
    features <- features[! features %in% features.computed]
  }
  cells <- rownames(spatial.location)
  cell.mismatch <- setdiff(cells, Cells(x = object, layer = layer))
  if (length(cell.mismatch) > 0L) {
    stop(
      "At least some of the row names in 'spatial.location' do not match cells in the '",
      layer, "' layer; check that the row names of 'spatial.location' are cell names.",
      call. = FALSE
    )
  }
  data <- LayerData(object, layer = layer, cells = cells, features = features)
  data <- as.matrix(x = data)
  # RowVar() is C++ and aborts the session on an empty matrix rather than
  # raising an R error, so catch any remaining path that produces one. One
  # column is caught here too: RowVar() divides by 'ncol - 1', so a single
  # cell yields NaN for every feature, and the NA subscript that follows
  # fails further down with a message that says nothing about the cause.
  if (ncol(x = data) < 2) {
    stop(
      "Fewer than two cells were returned from the '", layer, "' layer; ",
      "finding spatially variable features requires at least two cells.",
      call. = FALSE
    )
  }
  # Keep a single surviving feature as a one-row matrix; without drop = FALSE it
  # becomes a vector and `nrow()` returns NULL, triggering "argument is of length
  # zero" in the `if (nrow(x = data) != 0)` check below.
  data <- data[RowVar(x = data) > 0, , drop = FALSE]
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
  if (selection.method == "markvariogram" &&
      "markvariogram" %in% names(x = Misc(object = object))) {
    svf.info <- c(svf.info, Misc(object = object, slot = "markvariogram"))
  }
  # Nothing was computed and nothing was cached: every requested feature had
  # zero variance. The branches below index into 'svf.info', so return here
  # rather than letting them fail on a NULL.
  if (!length(x = svf.info)) {
    warning(
      "None of the requested features vary across the given cells in the '",
      layer, "' layer; returning the object unchanged.",
      call. = FALSE
    )
    return(object)
  }
  if (selection.method == "markvariogram") {
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
  object[[names(x = svf.info)]] <- svf.info
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
  layer = "scale.data",
  # Using `deprecated()` as the default for any arguments will break the
  # `LogSeuratCommand` call at the end of this method.
  slot = NULL,
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
  if (!is.null(slot)) {
    deprecate_soft(
      when = '5.3.0',
      what = 'FindSpatiallyVariableFeatures(slot = )',
      with = 'FindSpatiallyVariableFeatures(layer = )'
    )
    layer <- slot %||% layer
  }

  assay <- assay %||% DefaultAssay(object = object)
  selection.method <- match.arg(arg = selection.method)
  images <- Images(object = object, assay = assay)
  if (is.null(x = image)) {
    if (!length(x = images)) {
      stop(
        "No image is associated with assay ", sQuote(x = assay, q = FALSE),
        call. = FALSE
      )
    }
    image <- images[[1L]]
  }
  features <- features %||% Features(object, assay = assay, layer = layer)
  tc <- GetTissueCoordinates(object = object[[image]])
  if ('cell' %in% colnames(x = tc)) {
    cell.names <- as.character(x = tc[['cell']])
    if (anyNA(cell.names) || any(!nzchar(x = cell.names)) || anyDuplicated(x = cell.names)) {
      stop(
        "Spatial coordinates must contain exactly one row per cell. Use ",
        "'DefaultBoundary()' to select a centroid-based boundary or provide ",
        "'spatial.location'.",
        call. = FALSE
      )
    }
    rownames(x = tc) <- cell.names
    tc <- tc[, setdiff(x = colnames(x = tc), y = 'cell'), drop = FALSE]
  }

  object[[assay]] <- FindSpatiallyVariableFeatures(
    object = object[[assay]],
    layer = layer,
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

  object <- LogSeuratCommand(object)

  return(object)
}

#' @rdname LogNormalize
#' @method LogNormalize data.frame
#' @export
#'
LogNormalize.data.frame <- function(
  data,
  scale.factor = 1e4,
  margin = 2L,
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
  margin = 2L,
  verbose = TRUE,
  ...
) {
  if (!inherits(x = data, what = 'dgCMatrix')) {
    data <- as(object = data, Class = "dgCMatrix")
  }
  # call Rcpp function to normalize
  if (verbose) {
    cat("Performing log-normalization\n", file = stderr())
  }
  nthreads <- getThreads(verbose = FALSE)
  # LogNorm takes only x and p slots of the dgCMatrix
  # then replaces the x slot with normalized values - all other slots can be reused
  data@x <- LogNorm(x = data@x, p = data@p, scale_factor = scale.factor, nthreads = nthreads, display_progress = verbose)
  return(data)
}

#' @param normalization.method Method for normalization.
#'  \itemize{
#'   \item \dQuote{\code{LogNormalize}}: Feature counts for each cell are
#'    divided by the total counts for that cell and multiplied by the
#'    \code{scale.factor}. This is then natural-log transformed using \code{log1p}
#'   \item \dQuote{\code{CLR}}: Applies a centered log ratio transformation
#'   \item \dQuote{\code{RC}}: Relative counts. Feature counts for each cell
#'    are divided by the total counts for that cell and multiplied by the
#'    \code{scale.factor}. No log-transformation is applied. For counts per
#'    million (CPM) set \code{scale.factor = 1e6}
#' }
#' @param scale.factor Sets the scale factor for cell-level normalization
#' @param margin If performing CLR normalization, normalize across features (1) or cells (2)
# @param across If performing CLR normalization, normalize across either "features" or "cells".
#' @param block.size How many cells should be run in each chunk - will be split evenly across threads.
#' If supplied, temporarily sets the thread count to \code{ceiling(ncol(object) / block.size)},
#' capped at the number of available cores.
#' @param verbose Whether to display a progress bar
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
  if (!is.null(block.size)) {
    stopifnot("'block.size' must be a positive number" =
              length(x = block.size) == 1L && is.numeric(x = block.size) && block.size > 0)
    req_nthreads <- ceiling(length(Cells(x = object)) / block.size)
    available <- future::availableCores()
    if (is.na(x = available)) {
      available <- 1L
    }
    req_nthreads <- max(1L, min(as.integer(x = req_nthreads), as.integer(x = available)))
    old_options <- options(Seurat.nthreads = req_nthreads)
    on.exit(expr = options(old_options), add = TRUE)
  }
  normalized.data <- switch(EXPR = normalization.method,
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
                            ),
                            'RC' = RelativeCounts(
                              data = object,
                              scale.factor = scale.factor,
                              verbose = verbose
                            ),
                            stop("Unknown normalization method: ", normalization.method))
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
    layer = 'data',
    new.data = NormalizeData(
      object = GetAssayData(object = object, layer = 'counts'),
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
  if (!(assay %in% Assays(object = object))) {
    stop("Assay ", assay, " not found in object")
  }
  slot(object, "assays")[[assay]] <- NormalizeData(
    object = slot(object, "assays")[[assay]],
    normalization.method = normalization.method,
    scale.factor = scale.factor,
    verbose = verbose,
    margin = margin,
    ...
  )
  object <- LogSeuratCommand(object = object)
  return(object)
}

FastSparseRowScaleInternal <- FastSparseRowScale

FastSparseRowScale <- function(
  mat,
  features = as.integer(x = c()),
  scale = TRUE,
  center = TRUE,
  scale_max = 10,
  nthreads = 1L,
  display_progress = FALSE
) {
  if (inherits(x = mat, what = "dgTMatrix")) {
    mat <- as(object = mat, Class = "dgCMatrix")
  }
  if (!inherits(x = mat, what = "dgCMatrix")) {
    stop("FastSparseRowScale requires a dgCMatrix or dgTMatrix", call. = FALSE)
  }
  rows <- nrow(x = mat)
  cols <- ncol(x = mat)
  result <- FastSparseRowScaleInternal(
    x = mat@x,
    i = mat@i,
    p = mat@p,
    rows = rows,
    cols = cols,
    features,
    scale,
    center,
    scale_max,
    nthreads,
    display_progress
  )
  if (isTRUE(x = scale)) {
    selected <- if (length(x = features)) {
      features + 1L
    } else {
      seq_len(length.out = rows)
    }
    row.sum <- Matrix::rowSums(x = mat[selected, , drop = FALSE])
    row.sq.sum <- Matrix::rowSums(x = mat[selected, , drop = FALSE] ^ 2)
    variance.numerator <- if (isTRUE(x = center)) {
      row.sq.sum - (row.sum * row.sum / cols)
    } else {
      row.sq.sum
    }
    sigma <- sqrt(x = variance.numerator / (cols - 1))
    invalid <- which(x = !(sigma > 0) | is.na(x = sigma))
    if (length(x = invalid)) {
      result[invalid, ] <- NaN
    }
  }
  return(result)
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
  # pass the full matrix plus feature indices into C++ and select the requested rows there
  object.names <- list(features, colnames(x = object))
  min.cells.to.block <- min(min.cells.to.block, ncol(x = object))
  suppressWarnings(expr = Parenting(
    parent.find = "ScaleData.Assay",
    features = features,
    min.cells.to.block = min.cells.to.block
  ))
  split.by <- split.by %||% TRUE
  split.cells <- split(x = colnames(x = object), f = split.by)
  nthreads <- getThreads(verbose = FALSE)
  CheckGC()
  # When no regression is requested and the data are not split,
  # compute the row statistics and materialise the final dense matrix in a single C++ call.
  if (
    is.null(x = vars.to.regress) &&
    is.null(x = latent.data) &&
    length(x = split.cells) == 1 &&
    !anyNA(x = object)
  ) {
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
    # 0-based row indices of the requested features within the full matrix.
    feature.idx <- match(x = features, table = rownames(x = object)) - 1L
    if (inherits(x = object, what = 'dgTMatrix')) {
      object <- as(object = object, Class = 'dgCMatrix')
    }
    if (is(object = object, class2 = 'dgCMatrix')) {
      scaled.data <- FastSparseRowScaleInternal(
        x = object@x,
        i = object@i,
        p = object@p,
        rows = nrow(x = object),
        cols = ncol(x = object),
        features = feature.idx,
        scale = do.scale,
        center = do.center,
        scale_max = scale.max,
        nthreads = nthreads,
        display_progress = verbose
      )
    } else {
      scaled.data <- FastDenseRowScale(
        mat = as.matrix(x = object),
        features = feature.idx,
        scale = do.scale,
        center = do.center,
        scale_max = scale.max,
        nthreads = nthreads,
        display_progress = verbose
      )
    }
    dimnames(x = scaled.data) <- object.names
    return(scaled.data)
  }
  object <- object[features, , drop = FALSE]
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
            features.regress = NULL,
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
            features.regress = NULL,
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
    vars.to.regress = NULL,
    latent.data = NULL,
    do.scale = TRUE,
    do.center = TRUE,
    scale.max = 10,
    verbose = TRUE,
    ...
) {
  features <- features %||% rownames(x = object)
  features <- as.vector(x = intersect(x = features, y = rownames(x = object)))
  object <- object[features, , drop = FALSE]

  # Handle covariate regression using BPCells::regress_out
  if (!is.null(x = vars.to.regress)) {
    if (is.null(x = latent.data)) {
      latent.data <- data.frame(row.names = colnames(x = object))
    } else {
      latent.data <- latent.data[colnames(x = object), , drop = FALSE]
      rownames(x = latent.data) <- colnames(x = object)
    }
    # Check if any vars.to.regress are features in the matrix
    if (any(vars.to.regress %in% rownames(x = object))) {
      feature_vars <- vars.to.regress[vars.to.regress %in% rownames(x = object)]
      # For IterableMatrix, convert the subset to a regular matrix for latent.data
      feature_data <- t(as.matrix(object[feature_vars, , drop = FALSE]))
      latent.data <- cbind(latent.data, feature_data)
    }
    # Validate that we have the requested variables
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
    # Use BPCells regress_out function
    regress_data <- latent.data[, vars.to.regress, drop = FALSE]
    object <- BPCells::regress_out(mat = object, latent_data = regress_data, prediction_axis = "row")
  }

  # Proceed with scaling/centering
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

  if (do.center) {
    features.mean <- BPCells::matrix_stats(
      matrix = object,
      row_stats = 'mean')$row_stats['mean',]
  } else {
    features.mean <- 0
  }
  if (do.scale) {
    features.var <- BPCells::matrix_stats(
      matrix = object,
      row_stats = 'variance')$row_stats['variance',]
    if (do.center) {
      features.sd <- sqrt(features.var)
    } else {
      # When not centering, scale by sqrt(sum(x^2) / (n-1)) to match
      # FastSparseRowScale behavior (Bessel-corrected root mean square)
      n <- ncol(object)
      features.row.mean <- BPCells::matrix_stats(
        matrix = object,
        row_stats = 'mean')$row_stats['mean',]
      features.sd <- sqrt(features.var + n * features.row.mean^2 / (n - 1))
    }
    features.sd[features.sd == 0] <- 0.01
  } else {
    features.sd <- 1
  }
  if (scale.max != Inf && (do.scale || do.center)) {
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
    features <- rownames(x = GetAssayData(object = object, layer = slot.use))
  }
  object <- SetAssayData(
    object = object,
    layer = 'scale.data',
    new.data = ScaleData(
      object = GetAssayData(object = object, layer = slot.use),
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
  scale.data.cells <- colnames(x = GetAssayData(object = object, assay = assay, layer = "scale.data"))
  if (length(x = setdiff(x = model.cells, y =  scale.data.cells)) == 0) {
  existing_features <- names(x = which(x = ! apply(
    X = GetAssayData(object = object, assay = assay, layer = "scale.data")[, model.cells],
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
    umi <- GetAssayData(object = object, assay = umi.assay, layer = "counts" )[features_to_compute, model.cells, drop = FALSE]
  } else {
    warning(
      "In the SCTModel ", SCTModel, ", the following ", length(x = diff_features),
      " features do not exist in the counts slot: ", paste(diff_features, collapse = ", ")
    )
    if (length(x = intersect_features) == 0) {
      new_residual <- matrix(
        data = NA,
        nrow = length(x = features_to_compute),
        ncol = length(x = model.cells),
        dimnames = list(features_to_compute, model.cells)
      )
    } else {
      umi <- GetAssayData(object = object, assay = umi.assay, layer = "counts")[intersect_features, model.cells, drop = FALSE]
    }
  }
  clip.max <- max(clip.range)
  clip.min <- min(clip.range)
  if (exists(x = "umi", inherits = FALSE) && nrow(x = umi) > 0) {
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
    if (!identical(dim(x = new_residual), dim(x = umi))) {
      new_residual <- matrix(
        data = new_residual,
        nrow = nrow(x = umi),
        ncol = ncol(x = umi),
        dimnames = dimnames(x = umi)
      )
    } else {
      dimnames(x = new_residual) <- dimnames(x = umi)
    }
    # centered data
    new_residual <- sweep(
      x = new_residual,
      MARGIN = 1,
      STATS = rowMeans(x = new_residual),
      FUN = "-"
    )
  } else if (!exists(x = "new_residual", inherits = FALSE)) {
    new_residual <- matrix(data = NA, nrow = 0, ncol = length(x = model.cells), dimnames = list(c(), model.cells))
  }
  if (length(x = diff_features) > 0 && length(x = intersect_features) > 0) {
    padded_residual <- matrix(
      data = NA,
      nrow = length(x = features_to_compute),
      ncol = length(x = model.cells),
      dimnames = list(features_to_compute, model.cells)
    )
    padded_residual[rownames(x = new_residual), colnames(x = new_residual)] <- new_residual
    new_residual <- padded_residual
  }
  old.features <- setdiff(x = new_features, y = features_to_compute)
  if (length(x = old.features) > 0) {
    old_residuals <- GetAssayData(object = object[[assay]], layer = "scale.data")[old.features, model.cells, drop = FALSE]
    combined_residual <- matrix(
      data = NA,
      nrow = length(x = new_features),
      ncol = length(x = model.cells),
      dimnames = list(new_features, model.cells)
    )
    combined_residual[rownames(x = new_residual), colnames(x = new_residual)] <- new_residual
    combined_residual[rownames(x = old_residuals), colnames(x = old_residuals)] <- old_residuals
    new_residual <- combined_residual
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
  if (nrow(x = data.expr) == 0) {
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
    features.regress <- seq_len(length.out = nrow(x = data.expr))
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
