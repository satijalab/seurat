#' @include generics.R
#'
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
#'
#' @author Robert A. Amezquita, \email{robert.amezquita@fredhutch.org}
#' @seealso \code{\link{BarcodeInflectionsPlot}} \code{\link{SubsetByBarcodeInflections}}
#'
#' @examples
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

#' Convert a peak matrix to a gene activity matrix
#'
#' This function will take in a peak matrix and an annotation file (gtf) and collapse the peak
#' matrix to a gene activity matrix. It makes the simplifying assumption that all counts in the gene
#' body plus X kb up and or downstream should be attributed to that gene.
#'
#' @param peak.matrix Matrix of peak counts
#' @param annotation.file Path to GTF annotation file
#' @param seq.levels Which seqlevels to keep (corresponds to chromosomes usually)
#' @param include.body Include the gene body?
#' @param upstream Number of bases upstream to consider
#' @param downstream Number of bases downstream to consider
#' @param verbose Print progress/messages
#'
#' @importFrom future nbrOfWorkers
#' @export
#'
CreateGeneActivityMatrix <- function(
  peak.matrix,
  annotation.file,
  seq.levels = c(1:22, "X", "Y"),
  include.body = TRUE,
  upstream = 2000,
  downstream = 0,
  verbose = TRUE
) {
  if (!PackageCheck('GenomicRanges', error = FALSE)) {
    stop("Please install GenomicRanges from Bioconductor.")
  }
  if (!PackageCheck('rtracklayer', error = FALSE)) {
    stop("Please install rtracklayer from Bioconductor.")
  }
  # convert peak matrix to GRanges object
  peak.df <- rownames(x = peak.matrix)
  peak.df <- do.call(what = rbind, args = strsplit(x = gsub(peak.df, pattern = ":", replacement = "-"), split = "-"))
  peak.df <- as.data.frame(x = peak.df)
  colnames(x = peak.df) <- c("chromosome", 'start', 'end')
  peaks.gr <- GenomicRanges::makeGRangesFromDataFrame(df = peak.df)

  # if any peaks start at 0, change to 1
  # otherwise GenomicRanges::distanceToNearest will not work
  BiocGenerics::start(peaks.gr[BiocGenerics::start(peaks.gr) == 0, ]) <- 1

  # get annotation file, select genes
  gtf <- rtracklayer::import(con = annotation.file)
  gtf <- GenomeInfoDb::keepSeqlevels(x = gtf, value = seq.levels, pruning.mode = 'coarse')

  # change seqlevelsStyle if not the same
  if (!any(GenomeInfoDb::seqlevelsStyle(x = gtf) == GenomeInfoDb::seqlevelsStyle(x = peaks.gr))) {
    GenomeInfoDb::seqlevelsStyle(gtf) <- GenomeInfoDb::seqlevelsStyle(peaks.gr)
  }
  gtf.genes <- gtf[gtf$type == 'gene']

  # Extend definition up/downstream
  if (include.body) {
    gtf.body_prom <- Extend(x = gtf.genes, upstream = upstream, downstream = downstream)
  } else {
    gtf.body_prom <- SummarizedExperiment::promoters(x = gtf.genes, upstream = upstream, downstream = downstream)
  }
  gene.distances <- GenomicRanges::distanceToNearest(x = peaks.gr, subject = gtf.body_prom)
  keep.overlaps <- gene.distances[rtracklayer::mcols(x = gene.distances)$distance == 0]
  peak.ids <- peaks.gr[S4Vectors::queryHits(x = keep.overlaps)]
  gene.ids <- gtf.genes[S4Vectors::subjectHits(x = keep.overlaps)]

  # Some GTF rows will not have gene_name attribute
  # Replace it by gene_id attribute
  gene.ids$gene_name[is.na(gene.ids$gene_name)] <- gene.ids$gene_id[is.na(gene.ids$gene_name)]

  peak.ids$gene.name <- gene.ids$gene_name
  peak.ids <- as.data.frame(x = peak.ids)
  peak.ids$peak <- rownames(peak.matrix)[S4Vectors::queryHits(x = keep.overlaps)]
  annotations <- peak.ids[, c('peak', 'gene.name')]
  colnames(x = annotations) <- c('feature', 'new_feature')

  # collapse into expression matrix
  peak.matrix <- as(object = peak.matrix, Class = 'matrix')
  all.features <- unique(x = annotations$new_feature)

  if (nbrOfWorkers() > 1) {
    mysapply <- future_sapply
  } else {
    mysapply <- ifelse(test = verbose, yes = pbsapply, no = sapply)
  }
  newmat <- mysapply(X = 1:length(x = all.features), FUN = function(x){
    features.use <- annotations[annotations$new_feature == all.features[[x]], ]$feature
    submat <- peak.matrix[features.use, ]
    if (length(x = features.use) > 1) {
      return(Matrix::colSums(x = submat))
    } else {
      return(submat)
    }
  })
  newmat <- t(x = newmat)
  rownames(x = newmat) <- all.features
  colnames(x = newmat) <- colnames(x = peak.matrix)
  return(as(object = newmat, Class = 'dgCMatrix'))
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
#' @param umi.assay Name of the assay of the seurat object containing UMI matrix and the default is
#' RNA
#' @param clip.range Numeric of length two specifying the min and max values the Pearson residual
#' will be clipped to
#' @param replace.value Recalculate residuals for all features, even if they are already present.
#' Useful if you want to change the clip.range.
#' @param verbose Whether to print messages and progress bars
#'
#' @return Returns a Seurat object containing  pearson residuals of added features in its scale.data
#'
#' @importFrom sctransform get_residuals
#'
#' @export
#'
#' @seealso \code{\link[sctransform]{get_residuals}}
#'
#' @examples
#' pbmc_small <- SCTransform(object = pbmc_small, variable.features.n = 20)
#' pbmc_small <- GetResidual(object = pbmc_small, features = c('MS4A1', 'TCL1A'))
#'
GetResidual <- function(
  object,
  features,
  assay = "SCT",
  umi.assay = NULL,
  clip.range = NULL,
  replace.value = FALSE,
  verbose = TRUE
) {
  if (!IsSCT(assay = object[[assay]])) {
    stop(assay, " assay was not generated by SCTransform")
  }
  umi.assay <- umi.assay %||% Misc(object = object[[assay]], slot = "umi.assay")
  umi.assay <- umi.assay %||% "RNA" # for object created in 3.1.1 or earlier, default to RNA
  if (replace.value) {
    new_features <- features
  } else {
    new_features <- setdiff(
      x = features,
      y = rownames(x = GetAssayData(object = object, assay = assay, slot = "scale.data"))
    )
  }
  if (length(x = new_features) == 0) {
    if (verbose) {
      message("Pearson residuals of input features exist already")
    }
  } else {
    if (is.null(x = Misc(object = object[[assay]], slot = 'vst.set'))) {
      vst_out <- Misc(object = object[[assay]], slot = 'vst.out')
      # filter cells not in the object but in the SCT model
      vst_out$cell_attr <- vst_out$cell_attr[Cells(x = object), ]
      vst_out$cells_step1 <- intersect(x = vst_out$cells_step1, y = Cells(x = object))
      object <- GetResidualVstOut(
        object = object,
        assay = assay,
        umi.assay = umi.assay,
        new_features = new_features,
        vst_out = vst_out,
        clip.range = clip.range,
        verbose = verbose
      )
  } else {
    # Calculate Pearson Residual from integrated object SCT assay
    vst.set <- Misc(object = object[[assay]], slot = 'vst.set')
    scale.data <- GetAssayData(
      object = object,
      assay = assay,
      slot = "scale.data"
    )

    vst_set_genes <-  sapply(1:length(vst.set), function(x) rownames(vst.set[[x]]$model_pars_fit))
    vst_set_genes <- Reduce(intersect, vst_set_genes)
    diff_features <- setdiff(
      x = new_features,
      y = vst_set_genes
    )
    if (length(x = diff_features) !=0) {
      warning(
        "The following ", length(x = diff_features),
        " features do not exist in all SCT models: ",
        paste(diff_features, collapse = " ")
      )
    }
    new_features <- intersect(
      x = new_features,
      y = vst_set_genes
    )
    if (length(new_features) != 0){
    object <- SetAssayData(
      object = object,
      assay = assay,
      slot = "scale.data",
      new.data = scale.data[!rownames(x = scale.data) %in% new_features, , drop = FALSE]
    )
    new.scale.data <- matrix(nrow = length(new_features), ncol = 0)
    rownames(x = new.scale.data) <- new_features
    for (v in 1:length(x = vst.set)) {
      vst_out <- vst.set[[v]]
      # confirm that cells from SCT model also exist in the integrated object
      cells.v <- intersect(x = rownames(x = vst_out$cell_attr), y = Cells(x = object))
      vst_out$cell_attr <- vst_out$cell_attr[cells.v, ]
      vst_out$cells_step1 <- intersect(x = vst_out$cells_step1, y = cells.v)
      object.v <- subset(x = object, cells = cells.v)

      object.v <- GetResidualVstOut(
        object = object.v,
        assay = assay,
        umi.assay = umi.assay[[v]],
        new_features = new_features,
        vst_out = vst_out,
        clip.range = clip.range,
        verbose = verbose
      )
      new.scale.data <- cbind(
        new.scale.data,
        GetAssayData(object = object.v, assay = assay, slot ="scale.data" )[new_features, , drop = FALSE]
      )
    }
    object <- SetAssayData(
      object = object,
      assay = assay,
      slot = "scale.data",
      new.data = rbind(
        GetAssayData(object = object, slot = 'scale.data', assay = assay),
        new.scale.data
      )
    )
    }
  }
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
  if (is.data.frame(x = data)) {
    data <- as.matrix(x = data)
  }
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
#' @import Matrix
#'
#' @export
#'
#' @references \url{https://www.biorxiv.org/content/early/2018/08/08/387241}
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

#' Load in data from Alevin pipeline
#'
#' Enables easy loading of csv format matrix provided by Alevin
#' ran with `--dumpCsvCounts` flags.
#'
#' @param base.path Directory containing the alevin/quant_mat*
#' files provided by Alevin.
#'
#' @return Returns a matrix with rows and columns labeled
#'
#' @importFrom utils read.csv read.delim
#' @export
#'
#' @author Avi Srivastava
#'
#' @examples
#' \dontrun{
#' data_dir <- 'path/to/output/directory'
#' list.files(data_dir) # Should show alevin/quants_mat* files
#' expression_matrix <- ReadAlevinCsv(base.path = data_dir)
#' seurat_object = CreateSeuratObject(counts = expression_matrix)
#' }
#'
ReadAlevinCsv <- function(base.path) {
  if (!dir.exists(base.path)) {
    stop("Directory provided does not exist")
  }
  barcode.loc <- file.path(base.path, "alevin", "quants_mat_rows.txt")
  gene.loc <- file.path(base.path, "alevin", "quants_mat_cols.txt")
  matrix.loc <- file.path( base.path, "alevin", "quants_mat.csv" )
  if (!file.exists(barcode.loc)) {
    stop("Barcode file missing")
  }
  if (!file.exists(gene.loc)) {
    stop("Gene name file missing")
  }
  if (!file.exists(matrix.loc)) {
    stop("Expression matrix file missing")
  }
  matrix <- as.matrix(x = read.csv(file = matrix.loc, header = FALSE))
  matrix <- t(x = matrix[, 1:ncol(x = matrix) - 1])
  cell.names <- readLines(con = barcode.loc)
  gene.names <- readLines(con = gene.loc)
  colnames(x = matrix) <- cell.names
  rownames(x = matrix) <- gene.names
  matrix[is.na(x = matrix)] <- 0
  return(matrix)
}

#' Load in data from Alevin pipeline
#'
#' Enables easy loading of binary format matrix provided by Alevin
#'
#' @param base.path Directory containing the alevin/quant_mat*
#' files provided by Alevin.
#'
#' @return Returns a matrix with rows and columns labeled
#'
#' @export
#'
#' @author Avi Srivastava
#'
#' @examples
#' \dontrun{
#' data_dir <- 'path/to/output/directory'
#' list.files(data_dir) # Should show alevin/quants_mat* files
#' expression_matrix <- ReadAlevin(base.path = data_dir)
#' seurat_object = CreateSeuratObject(counts = expression_matrix)
#' }
#'
ReadAlevin <- function(base.path) {
  if (!dir.exists(base.path)) {
    stop("Directory provided does not exist")
  }
  barcode.loc <- file.path(base.path, "alevin", "quants_mat_rows.txt")
  gene.loc <- file.path(base.path, "alevin", "quants_mat_cols.txt")
  matrix.loc <- file.path(base.path, "alevin", "quants_mat.gz")
  if (!file.exists(barcode.loc)) {
    stop("Barcode file missing")
  }
  if (!file.exists(gene.loc)) {
    stop("Gene name file missing")
  }
  if (!file.exists(matrix.loc)) {
    stop("Expression matrix file missing")
  }
  cell.names <- readLines(con = barcode.loc)
  gene.names <- readLines(con = gene.loc)
  num.cells <- length(x = cell.names)
  num.genes <- length(x = gene.names)
  out.matrix <- matrix(data = NA, nrow = num.genes, ncol = num.cells)
  con <- gzcon(con = file(description = matrix.loc, open = "rb"))
  total.molecules <- 0.0
  for (n in seq_len(length.out = num.cells)) {
    out.matrix[, n] <- readBin(
      con = con,
      what = double(),
      endian = "little",
      n = num.genes
    )
    total.molecules <- total.molecules + sum(out.matrix[, n])
  }
  colnames(x = out.matrix) <- cell.names
  rownames(x = out.matrix) <- gene.names
  message("Found total ", total.molecules, " molecules")
  return(out.matrix)
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
#' @param unique.features Make feature names unique (default TRUE)
#'
#' @return If features.csv indicates the data has multiple data types, a list
#'   containing a sparse matrix of the data from each type will be returned.
#'   Otherwise a sparse matrix containing the expression data will be returned.
#'
#' @importFrom Matrix readMM
#'
#' @export
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
Read10X <- function(data.dir = NULL, gene.column = 2, unique.features = TRUE) {
  full.data <- list()
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
    cell.names <- readLines(barcode.loc)
    if (all(grepl(pattern = "\\-1$", x = cell.names))) {
      cell.names <- as.vector(x = as.character(x = sapply(
        X = cell.names,
        FUN = ExtractField,
        field = 1,
        delim = "-"
      )))
    }
    if (is.null(x = names(x = data.dir))) {
      if (i < 2) {
        colnames(x = data) <- cell.names
      } else {
        colnames(x = data) <- paste0(i, "_", cell.names)
      }
    } else {
      colnames(x = data) <- paste0(names(x = data.dir)[i], "_", cell.names)
    }
    feature.names <- read.delim(
      file = ifelse(test = pre_ver_3, yes = gene.loc, no = features.loc),
      header = FALSE,
      stringsAsFactors = FALSE
    )
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
          return(data[data_types == l, ])
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
    list_of_data[[j]] <- as(object = list_of_data[[j]], Class = "dgCMatrix")
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
  if (!infile$attr_exists("PYTABLES_FORMAT_VERSION")) {
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
      giveCsparse = FALSE
    )
    if (unique.features) {
      features <- make.unique(names = features)
    }
    rownames(x = sparse.mat) <- features
    colnames(x = sparse.mat) <- barcodes[]
    sparse.mat <- as(object = sparse.mat, Class = 'dgCMatrix')
    # Split v3 multimodal
    if (infile$exists(name = paste0(genome, '/features'))) {
      types <- infile[[paste0(genome, '/features/feature_type')]][]
      types.unique <- unique(x = types)
      if (length(x = types.unique) > 1) {
        message("Genome ", genome, " has multiple modalities, returning a list of matrices for this genome")
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
#' @import Matrix
#' @importFrom methods as
#'
#' @export
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
    data <- as(object = data, Class = "dgCMatrix")
  }
  if (verbose) {
    cat("Performing relative-counts-normalization\n", file = stderr())
  }
  norm.data <- data
  norm.data@x <- norm.data@x / rep.int(colSums(norm.data), diff(norm.data@p)) * scale.factor
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

#' Use regularized negative binomial regression to normalize UMI count data
#'
#' This function calls sctransform::vst. The sctransform package is available at
#' https://github.com/ChristophH/sctransform.
#' Use this function as an alternative to the NormalizeData,
#' FindVariableFeatures, ScaleData workflow. Results are saved in a new assay
#' (named SCT by default) with counts being (corrected) counts, data being log1p(counts),
#' scale.data being pearson residuals; sctransform::vst intermediate results are saved
#' in misc slot of new assay.
#'
#' @param object A seurat object
#' @param assay Name of assay to pull the count data from; default is 'RNA'
#' @param new.assay.name Name for the new assay containing the normalized data
#' @param do.correct.umi Place corrected UMI matrix in assay counts slot; default is TRUE
#' @param ncells Number of subsampling cells used to build NB regression; default is NULL
#' @param variable.features.n Use this many features as variable features after
#' ranking by residual variance; default is 3000
#' @param variable.features.rv.th Instead of setting a fixed number of variable features,
#' use this residual variance cutoff; this is only used when \code{variable.features.n}
#' is set to NULL; default is 1.3
#' @param vars.to.regress Variables to regress out in a second non-regularized linear
#' regression. For example, percent.mito. Default is NULL
#' @param do.scale Whether to scale residuals to have unit variance; default is FALSE
#' @param do.center Whether to center residuals to have mean zero; default is TRUE
#' @param clip.range Range to clip the residuals to; default is \code{c(-sqrt(n/30), sqrt(n/30))},
#' where n is the number of cells
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
#' @importFrom sctransform vst get_residual_var get_residuals correct_counts
#'
#' @seealso \code{\link[sctransform]{correct_counts}} \code{\link[sctransform]{get_residuals}}
#' @export
#'
#' @examples
#' SCTransform(object = pbmc_small)
#'
SCTransform <- function(
  object,
  assay = 'RNA',
  new.assay.name = 'SCT',
  do.correct.umi = TRUE,
  ncells = NULL,
  variable.features.n = 3000,
  variable.features.rv.th = 1.3,
  vars.to.regress = NULL,
  do.scale = FALSE,
  do.center = TRUE,
  clip.range = c(-sqrt(x = ncol(x = object[[assay]]) / 30), sqrt(x = ncol(x = object[[assay]]) / 30)),
  conserve.memory = FALSE,
  return.only.var.genes = TRUE,
  seed.use = 1448145,
  verbose = TRUE,
  ...
) {
  if (!is.null(x = seed.use)) {
    set.seed(seed = seed.use)
  }
  assay <- assay %||% DefaultAssay(object = object)
  assay.obj <- GetAssay(object = object, assay = assay)
  umi <- GetAssayData(object = assay.obj, slot = 'counts')
  cell.attr <- slot(object = object, name = 'meta.data')
  vst.args <- list(...)
  # check for batch_var in meta data
  if ('batch_var' %in% names(x = vst.args)) {
    if (!(vst.args[['batch_var']] %in% colnames(x = cell.attr))) {
      stop('batch_var not found in seurat object meta data')
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
  if (any(c('cell_attr', 'show_progress', 'return_cell_attr', 'return_gene_attr', 'return_corrected_umi') %in% names(x = vst.args))) {
    warning(
      'the following arguments will be ignored because they are set within this function:',
      paste(
        c(
          'cell_attr',
          'show_progress',
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
  vst.args[['umi']] <- umi
  vst.args[['cell_attr']] <- cell.attr
  vst.args[['show_progress']] <- verbose
  vst.args[['return_cell_attr']] <- TRUE
  vst.args[['return_gene_attr']] <- TRUE
  vst.args[['return_corrected_umi']] <- do.correct.umi
  vst.args[['n_cells']] <- ncells
  residual.type <- vst.args[['residual_type']] %||% 'pearson'
  res.clip.range <- vst.args[['res_clip_range']] %||% c(-sqrt(x = ncol(x = umi)), sqrt(x = ncol(x = umi)))
  if (conserve.memory) {
    return.only.var.genes <- TRUE
  }
  if (conserve.memory) {
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
  } else {
    vst.out <- do.call(what = 'vst', args = vst.args)
    feature.variance <- setNames(
      object = vst.out$gene_attr$residual_variance,
      nm = rownames(x = vst.out$gene_attr)
    )
  }
  if (verbose) {
    message('Determine variable features')
  }
  feature.variance <- sort(x = feature.variance, decreasing = TRUE)
  if (!is.null(x = variable.features.n)) {
    top.features <- names(x = feature.variance)[1:min(variable.features.n, length(x = feature.variance))]
  } else {
    top.features <- names(x = feature.variance)[feature.variance >= variable.features.rv.th]
  }
  if (verbose) {
    message('Set ', length(x = top.features), ' variable features')
  }
  if (conserve.memory) {
    # actually get the residuals this time
    if (verbose) {
      message("Return only variable features for scale.data slot of the output assay")
    }
    vst.out$y <- get_residuals(
      vst_out = vst.out,
      umi = umi[top.features, ],
      residual_type = residual.type,
      res_clip_range = res.clip.range
    )
    if (do.correct.umi & residual.type == 'pearson') {
      vst.out$umi_corrected <- correct_counts(
        x = vst.out,
        umi = umi,
        show_progress = verbose
      )
    }
  }
  # create output assay and put (corrected) umi counts in count slot
  if (do.correct.umi & residual.type == 'pearson') {
    if (verbose) {
      message('Place corrected count matrix in counts slot')
    }
    assay.out <- CreateAssayObject(counts = vst.out$umi_corrected)
    vst.out$umi_corrected <- NULL
  } else {
    assay.out <- CreateAssayObject(counts = umi)
  }
  # set the variable genes
  VariableFeatures(object = assay.out) <- top.features
  # put log1p transformed counts in data
  assay.out <- SetAssayData(
    object = assay.out,
    slot = 'data',
    new.data = log1p(x = GetAssayData(object = assay.out, slot = 'counts'))
  )
  if (return.only.var.genes & !conserve.memory) {
    scale.data <- vst.out$y[top.features, ]
  } else {
    scale.data <- vst.out$y
  }
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
  assay.out <- SetAssayData(
    object = assay.out,
    slot = 'scale.data',
    new.data = scale.data
  )
  # save vst output (except y) in @misc slot
  vst.out$y <- NULL
  # save clip.range into vst model
  vst.out$arguments$sct.clip.range <- clip.range
  Misc(object = assay.out, slot = 'vst.out') <- vst.out
  Misc(object = assay.out, slot = 'umi.assay') <- assay
  # also put gene attributes in meta.features
  assay.out[[paste0('sct.', names(x = vst.out$gene_attr))]] <- vst.out$gene_attr
  assay.out[['sct.variable']] <- rownames(x = assay.out[[]]) %in% top.features
  object[[new.assay.name]] <- assay.out
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
#'
#' @author Robert A. Amezquita, \email{robert.amezquita@fredhutch.org}
#' @seealso \code{\link{CalculateBarcodeInflections}} \code{\link{BarcodeInflectionsPlot}}
#'
#' @examples
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

#' Term frequency-inverse document frequency
#'
#' Normalize binary data per cell using the term frequency-inverse document frequency
#' normalization method (TF-IDF).
#' This is suitable for the normalization of binary ATAC peak datasets.
#'
#' @param data Matrix with the raw count data
#' @param verbose Print progress
#'
#' @return Returns a matrix with the normalized data
#'
#' @importFrom Matrix colSums rowSums
#'
#' @export
#'
#' @examples
#' mat <- matrix(data = rbinom(n = 25, size = 5, prob = 0.2), nrow = 5)
#' mat_norm <- TF.IDF(data = mat)
#'
TF.IDF <- function(data, verbose = TRUE) {
  if (is.data.frame(x = data)) {
    data <- as.matrix(x = data)
  }
  if (!inherits(x = data, what = 'dgCMatrix')) {
    data <- as(object = data, Class = "dgCMatrix")
  }
  if (verbose) {
    message("Performing TF-IDF normalization")
  }
  npeaks <- colSums(x = data)
  tf <- t(x = t(x = data) / npeaks)
  idf <- ncol(x = data) / rowSums(x = data)
  norm.data <- Diagonal(n = length(x = idf), x = idf) %*% tf
  norm.data[which(x = is.na(x = norm.data))] <- 0
  return(norm.data)
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
  verbose = TRUE,
  ...
) {
  CheckDots(...)
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

#' @inheritParams FindVariableFeatures.Assay
#' @param assay Assay to use
#'
#' @rdname FindVariableFeatures
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
    verbose = verbose,
    ...
  )
  object[[assay]] <- assay.data
  object <- LogSeuratCommand(object = object)
  return(object)
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
#' @export
#'
NormalizeData.default <- function(
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
  margin = 1,
  verbose = TRUE,
  ...
) {
  assay <- assay %||% DefaultAssay(object = object)
  assay.data <- GetAssay(object = object, assay = assay)
  assay.data <- NormalizeData(
    object = assay.data,
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

#' @param k  The rank of the rank-k approximation. Set to NULL for automated choice of k.
#' @param q  The number of additional power iterations in randomized SVD when
#' computing rank k approximation. By default, q=10.
#' @param quantile.prob The quantile probability to use when calculating threshold.
#' By default, quantile.prob = 0.001.
#' @param use.mkl Use the Intel MKL based implementation of SVD. Needs to be
#' installed from https://github.com/KlugerLab/rpca-mkl. \strong{Note}: this requires
#' the \href{https://github.com/satijalab/seurat-wrappers}{SeuratWrappers} implementation
#' of \code{RunALRA}
#' @param mkl.seed Only relevant if \code{use.mkl = TRUE}. Set the seed for the random
#' generator for the Intel MKL implementation of SVD. Any number <0 will
#' use the current timestamp. If \code{use.mkl = FALSE}, set the seed using
#' \code{\link{set.seed}()} function as usual.
#'
#' @rdname RunALRA
#' @export
#'
RunALRA.default <- function(
  object,
  k = NULL,
  q = 10,
  quantile.prob = 0.001,
  use.mkl = FALSE,
  mkl.seed = -1,
  ...
) {
  CheckDots(...)
  A.norm <- t(x = as.matrix(x = object))
  message("Identifying non-zero values")
  originally.nonzero <- A.norm > 0
  message("Computing Randomized SVD")
  if (use.mkl) {
    warning(
      "Using the Intel MKL-based implementation of SVD requires RunALRA from SeuratWrappers\n",
      "For more details, see https://github.com/satijalab/seurat-wrappers\n",
      "Continuing with standard SVD implementation",
      call. = FALSE,
      immediate. = TRUE
    )
  }
  fastDecomp.noc <- rsvd(A = A.norm, k = k, q = q)
  A.norm.rank.k <- fastDecomp.noc$u[, 1:k] %*%
    diag(x = fastDecomp.noc$d[1:k]) %*%
    t(x = fastDecomp.noc$v[,1:k])
  message(sprintf("Find the %f quantile of each gene", quantile.prob))
  A.norm.rank.k.mins <- abs(x = apply(
    X = A.norm.rank.k,
    MARGIN = 2,
    FUN = function(x) {
      return(quantile(x = x, probs = quantile.prob))
    }
  ))
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
  toscale <- !is.na(sigma.1) & !is.na(sigma.2) & !(sigma.1 == 0 & sigma.2 == 0) & !(sigma.1 == 0)
  message(sprintf(fmt = "Scaling all except for %d columns", sum(!toscale)))
  sigma.1.2 <- sigma.2 / sigma.1
  toadd <- -1 * mu.1 * sigma.2 / sigma.1 + mu.2
  A.norm.rank.k.temp <- A.norm.rank.k.cor[, toscale]
  A.norm.rank.k.temp <- Sweep(
    x = A.norm.rank.k.temp,
    MARGIN = 2,
    STATS = sigma.1.2[toscale],
    FUN = "*"
  )
  A.norm.rank.k.temp <- Sweep(
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
#' than the smallest dimension of the matrix. Default 100 or smallest dimension.
#' @param thresh The threshold for ''significance'' when choosing k. Default 1e-10.
#' @param noise.start Index for which all smaller singular values are considered noise.
#' Default K - 20.
#' @param q.k Number of additional power iterations when choosing k. Default 2.
#' @param k.only If TRUE, only computes optimal k WITHOUT performing ALRA
#'
#' @importFrom rsvd rsvd
#' @importFrom Matrix Matrix
#' @importFrom stats sd setNames quantile
#'
#' @rdname RunALRA
#' @export
#' @method RunALRA Seurat
#'
RunALRA.Seurat <- function(
  object,
  k = NULL,
  q = 10,
  quantile.prob = 0.001,
  use.mkl = FALSE,
  mkl.seed=-1,
  assay = NULL,
  slot = "data",
  setDefaultAssay = TRUE,
  genes.use = NULL,
  K = NULL,
  thresh = 6,
  noise.start = NULL,
  q.k = 2,
  k.only = FALSE,
  ...
) {
  if (!is.null(x = k) && k.only) {
    warning("Stop: k is already given, set k.only = FALSE or k = NULL")
  }
  genes.use <- genes.use %||% rownames(x = object)
  assay <- assay %||% DefaultAssay(object = object)
  alra.previous <- Tool(object = object, slot = 'RunALRA')
  alra.info <- list()
  # Check if k is already stored
  if (is.null(x = k) & !is.null(alra.previous[["k"]])) {
    k <- alra.previous[["k"]]
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
    if (use.mkl) {
      warning(
        "Using the Intel MKL-based implementation of SVD requires RunALRA from SeuratWrappers\n",
        "For more details, see https://github.com/satijalab/seurat-wrappers\n",
        "Continuing with standard SVD implementation",
        call. = FALSE,
        immediate. = TRUE
      )
    }
    rsvd.out <- rsvd(A = t(x = as.matrix(x = data.used)), k = K, q = q.k)
    diffs <- rsvd.out$d[1:(length(x = rsvd.out$d)-1)] - rsvd.out$d[2:length(x = rsvd.out$d)]
    mu <- mean(x = diffs[noise.svals - 1])
    sigma <- sd(x = diffs[noise.svals - 1])
    num_of_sds <- (diffs - mu) / sigma
    k <- max(which(x = num_of_sds > thresh))
    alra.info[["d"]] <- rsvd.out$d
    alra.info[["k"]] <- k
    alra.info[["diffs"]] <- diffs
    Tool(object = object) <- alra.info
  }
  if (k.only) {
    message("Chose rank k = ", k, ", WITHOUT performing ALRA")
    return(object)
  }
  message("Rank k = ", k)
  # Perform ALRA on data.used
  output.alra <- RunALRA(
    object = data.used,
    k = k,
    q = q,
    quantile.prob = quantile.prob,
    use.mkl = use.mkl,
    mkl.seed = mkl.seed
  )
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
#' Setting this can help reduce the effects of feautres that are only expressed in
#' a very small number of cells. If regressing out latent variables and using a
#' non-linear model, the default is 50.
#' @param block.size Default size for number of feautres to scale at in a single
#' computation. Increasing block.size may speed up calculations but at an
#' additional memory cost.
#' @param min.cells.to.block If object contains fewer than this number of cells,
#' don't block for scaling calculations.
#' @param verbose Displays a progress bar for scaling procedure
#'
#' @importFrom future.apply future_lapply
#'
#' @rdname ScaleData
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
        t(x = object[vars.to.regress[vars.to.regress %in% rownames(x = object)], ])
      )
    }
    # Currently, RegressOutMatrix will do nothing if latent.data = NULL
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
        data.scale <- scale.function(
          mat = object[features[block[1]:block[2]], split.cells[[group]], drop = FALSE],
          scale = do.scale,
          center = do.center,
          scale_max = scale.max,
          display_progress = FALSE
        )
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
        data.scale <- scale.function(
          mat = object[features[my.inds], split.cells[[x]], drop = FALSE],
          scale = do.scale,
          center = do.center,
          scale_max = scale.max,
          display_progress = FALSE
        )
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
  assay <- assay %||% DefaultAssay(object = object)
  assay.data <- GetAssay(object = object, assay = assay)
  if (any(vars.to.regress %in% colnames(x = object[[]]))) {
    latent.data <- object[[vars.to.regress[vars.to.regress %in% colnames(x = object[[]])]]]
  } else {
    latent.data <- NULL
  }
  if (is.character(x = split.by) && length(x = split.by) == 1) {
    split.by <- object[[split.by]]
  }
  assay.data <- ScaleData(
    object = assay.data,
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
    if (is.character(x = model)) {
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
#' @importFrom methods as
#' @importFrom pbapply pbapply
# @import Matrix
#
CustomNormalize <- function(data, custom_function, margin, verbose = TRUE) {
  if (is.data.frame(x = data)) {
    data <- as.matrix(x = data)
  }
  if (!inherits(x = data, what = 'dgCMatrix')) {
    data <- as(object = data, Class = "dgCMatrix")
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
    norm.data = t(x = norm.data)
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
# @return Returns a Seurat object containing  pearson residuals of added features in its scale.data
#
#' @importFrom sctransform get_residuals
#
GetResidualVstOut <- function(
  object,
  assay,
  umi.assay,
  new_features,
  vst_out,
  clip.range,
  verbose
) {
  diff_features <- setdiff(
    x = new_features,
    y = rownames(x = vst_out$model_pars_fit)
  )
  intersect_feature <- intersect(
    x = new_features,
    y = rownames(x = vst_out$model_pars_fit)
  )
  if (length(x = diff_features) == 0) {
    umi <- GetAssayData(object = object, assay = umi.assay, slot = "counts" )[new_features, , drop = FALSE]
  } else {

    warning(
      "The following ", length(x = diff_features),
      " features do not exist in the counts slot: ",
      paste(diff_features, collapse = " ")
    )

    if (length(x = intersect_feature) == 0) {
      return(object)
    }
    umi <- GetAssayData(object = object, assay = umi.assay, slot = "counts" )[intersect_feature, , drop = FALSE]
  }
  if (is.null(x = clip.range)) {
    if(length(vst_out$arguments$sct.clip.range)!=0 ){
    clip.max <- max(vst_out$arguments$sct.clip.range)
    clip.min <- min(vst_out$arguments$sct.clip.range)
    } else{
      clip.max <- max(vst_out$arguments$res_clip_range)
      clip.min <- min(vst_out$arguments$res_clip_range)
    }
  } else {
    clip.max <- max(clip.range)
    clip.min <- min(clip.range)
  }
  new_residual <- get_residuals(
    vst_out = vst_out,
    umi = umi,
    residual_type = "pearson",
    res_clip_range = c(clip.min, clip.max),
    show_progress = verbose
  )
  new_residual <- as.matrix(x = new_residual)
  # centered data
  new_residual <- new_residual - rowMeans(new_residual)
  # remove genes from the scale.data if genes are part of new_features
  scale.data <- GetAssayData(object = object, assay = assay, slot = "scale.data")
  object <- SetAssayData(
    object = object,
    assay = assay,
    slot = "scale.data",
    new.data = scale.data[!rownames(x = scale.data) %in% new_features, , drop = FALSE]
  )
  if (nrow(x = GetAssayData(object = object, slot = 'scale.data', assay = assay)) == 0 ) {
    object <- SetAssayData(
      object = object,
      slot = 'scale.data',
      new.data = new_residual,
      assay = assay
    )
  } else {
    object <- SetAssayData(
      object = object,
      slot = 'scale.data',
      new.data = rbind(
        GetAssayData(object = object, slot = 'scale.data', assay = assay),
        new_residual
      ),
      assay = assay
    )
  }
  return(object)
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
