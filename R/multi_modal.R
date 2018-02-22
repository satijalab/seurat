#' @include seurat.R
NULL

# Set up assay class to hold multimodal data sets

assay <- setClass(
  Class = "assay",
  slots = list(
    raw.data = "ANY",
    data = "ANY",
    scale.data = "ANY",
    key = "character",
    misc = "ANY",
    var.genes="vector",
    mean.var="data.frame"
  )
)

#' Accessor function for multimodal data
#'
#' Pull information for specified stored dimensional reduction analysis
#'
#' @param object Seurat object
#' @param assay.type Type of assay to fetch data for (default is RNA)
#' @param slot Specific information to pull (i.e. raw.data, data, scale.data,...). Default is data
#'
#' @return Returns assay data
#'
#' @export
#'
#' @examples
#' # Simulate CITE-Seq results
#' df <- t(x = data.frame(
#'   x = round(x = rnorm(n = 80, mean = 20, sd = 2)),
#'   y = round(x = rbinom(n = 80, size = 100, prob = 0.2))
#' ))
#' pbmc_small <- SetAssayData(
#'   object = pbmc_small,
#'   assay.type = 'CITE',
#'   new.data = df,
#'   slot = 'raw.data'
#' )
#' GetAssayData(object = pbmc_small, assay.type = 'CITE', slot = 'raw.data')
#'
GetAssayData <- function(object, assay.type = "RNA", slot = "data") {
  if (assay.type == "RNA") {
    if (slot == "raw.data") {
      to.return <- object@raw.data
    } else if (slot == "data") {
      to.return <- object@data
    } else if (slot == "scale.data") {
      if (length(x = object@scale.data) == 0) {
        stop("Object@scale.data has not been set. Run ScaleData() and then retry.")
      }
      to.return <- object@scale.data
    }
    if(is.null(to.return)) {
      return(to.return)
    }
    #note that we check for this to avoid a long subset for large matrices if it can be avoided
    if (length(x = object@cell.names) == ncol(to.return)) {
      return(to.return)
    }
    return(to.return[, object@cell.names])
  }
  if (! (assay.type %in% names(object@assay))) {
    stop(paste(assay.type, "data has not been added"))
  }
  if (! (slot %in% slotNames(eval(expr = parse(text = paste0("object@assay$", assay.type)))))) {
    stop(paste(slot, "slot doesn't exist"))
  }
  to.return <- (eval(expr = parse(text = paste0("object@assay$", assay.type, "@", slot))))
  if (length(x = object@cell.names) == ncol(x = to.return)) {
    return(to.return)
  }
  return(to.return[, object@cell.names])
}

#' Assay Data Mutator Function
#'
#' Store information for specified assay, for multimodal analysis. new.data
#' needs to have cells as the columns and measurement features (e.g. genes,
#' proteins, etc ...) as rows. Additionally, all the cell names in the new.data
#' must match the cell names in the object (object@@cell.names).
#'
#' @inheritParams GetAssayData
#' @param new.data New data to insert
#'
#' @return Seurat object with updated slot
#'
#' @export
#'
#' @examples
#' # Simulate CITE-Seq results
#' df <- t(x = data.frame(
#'   x = round(x = rnorm(n = 80, mean = 20, sd = 2)),
#'   y = round(x = rbinom(n = 80, size = 100, prob = 0.2))
#' ))
#' pbmc_small = SetAssayData(
#'   object = pbmc_small,
#'   assay.type = 'CITE',
#'   new.data = df,
#'   slot = 'raw.data'
#' )
#'
SetAssayData <- function(object, assay.type, slot, new.data) {
  if(! all(colnames(new.data) %in% object@cell.names)) {
    stop("Cell names in new assay data matrix don't exactly match the cell names of the object")
  }
  if (assay.type == "RNA") {
    if (slot == "raw.data") {
      (object@raw.data <- new.data)
    } else if (slot == "data") {
      (object@data <- new.data)
    } else if (slot == "scale.data") {
      (object@scale.data <- new.data)
    }
    return(object)
  }
  if (assay.type %in% names(object@assay)) {
    eval(expr = parse(text = paste0("object@assay$", assay.type, "@", slot, "<- new.data")))
  } else {
    new.assay <- new(Class = "assay")
    eval(expr = parse(text = paste0("new.assay@", slot, "<- new.data")))
    eval(expr = parse(text = paste0("object@assay$", assay.type, "<- new.assay")))
  }
  return(object)
}


#' Slim down a multi-species expression matrix, when only one species is primarily of interenst.
#'
#' Valuable for CITE-seq analyses, where we typically spike in rare populations of 'negative control' cells from a different species.
#'
#' @param data.matrix A UMI count matrix. Should contain rownames that start with the ensuing arguments prefix.1 or prefix.2
#' @param prefix.1 The prefix denoting rownames for the species of interest. Default is "HUMAN_". These rownames will have this prefix removed in the returned matrix.
#' @param prefix.controls The prefix denoting rownames for the species of 'negative control' cells. Default is "MOUSE_".
#' @param features.controls.toKeep How many of the most highly expressed (average) negative control features (by default, 100 mouse genes), should be kept? All other rownames starting with prefix.2 are discarded.
#' @return A UMI count matrix. Rownames that started with prefix.1 have this prefix discarded. For rownames starting with prefix.2, only the most highly expressed features are kept, and the prefix is kept. All other rows are retained.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' cbmc.rna.collapsed <- CollapseSpeciesExpressionMatrix(cbmc.rna)
#' }
CollapseSpeciesExpressionMatrix <- function(
  data.matrix,
  prefix.1 = "HUMAN_",
  prefix.controls = "MOUSE_",
  features.controls.toKeep = 100
) {
  data.matrix.1 <- SubsetRow(data = data.matrix, code = prefix.1)
  data.matrix.2 <- SubsetRow(data = data.matrix, code = prefix.controls)
  data.matrix.3 <- data.matrix[setdiff(
    x = rownames(x = data.matrix),
    y = c(rownames(x = data.matrix.1), rownames(x = data.matrix.2))
  ),]
  rownames(x = data.matrix.1) <- make.unique(names = sapply(
    X = rownames(x = data.matrix.1),
    FUN = function(x) {
      return(gsub(pattern = prefix.1, replacement = "", x = x))
    }
  ))
  control.sums <- rowSums(data.matrix.2)
  control.keep <- names(x = head(
    x = sort(x = control.sums, decreasing = TRUE),
    n = features.controls.toKeep
  ))
  control.matrix.keep <- data.matrix.2[control.keep, ]
  final.matrix <- rbind(data.matrix.1, control.matrix.keep, data.matrix.3)
  return(final.matrix)
}


#' Demultiplex samples based on data from cell 'hashing'
#'
#' Assign sample-of-origin for each cell, annotate doublets.
#'
#' @param object Seurat object. Assumes that the hash tag oligo (HTO) data has been added and normalized in the HTO slot.
#' @param percent_cutoff The quantile of inferred 'negative' distribution for each HTO - over which the cell is considered 'positive'. Default is 0.999
#' @param init_centers Initial number of clusters for kmeans of the HTO oligos. Default is the # of samples + 1 (to account for negatives)
#' @param cluster_nstarts nstarts value for the initial k-means clustering
#' @param k_function Clustering function for initial HTO grouping. Default is "kmeans", also support "clara" for fast k-medoids clustering on large applications
#' @param nsamples Number of samples to be drawn from the dataset used for clustering, for k_function = "clara"
#' @param print.output Prints the output
#' @param assay.type Naming of HTO assay
#' @param confidence_threshold The quantile of the inferred 'positive' distribution for each HTO. Cells that have lower counts than this threshold are labeled as uncertain in the confidence field. Default is 0.05
#' @return Seurat object. Demultiplexed information is stored in the object meta data.
#'
#' @importFrom stats pnbinom
#' @importFrom cluster clara
#' @importFrom fitdistrplus fitdist
#' @export
#'
#' @examples
#' \dontrun{
#' object <- HTODemux(object)
#' }
HTODemux <- function(
  object,
  percent_cutoff = 0.999,
  init_centers = NULL,
  cluster_nstarts = 100,
  k_function = "kmeans",
  nsamples = 100,
  print.output = TRUE,
  assay.type = "HTO",
  confidence_threshold = 0.05
) {
  #hashing
  hash_data <- GetAssayData(object = object, assay.type = assay.type)
  hash_raw_data <- GetAssayData(
    object = object,
    assay.type = assay.type,
    slot = "raw.data"
  )[, object@cell.names]
  ncenters <- SetIfNull(x = init_centers, default = nrow(hash_data) + 1)
  if (k_function == "kmeans") {
    hto_init_clusters <- kmeans(
      x = t(x = GetAssayData(object = object, assay.type = assay.type)),
      centers = ncenters,
      nstart = cluster_nstarts
    )
    #identify positive and negative signals for all HTO
    object <- SetIdent(
      object = object,
      cells.use = names(hto_init_clusters$cluster),
      ident.use = hto_init_clusters$cluster
    )
  } else {
    #use fast k-medoid clustering
    hto_init_clusters <- clara(
      x = t(x = GetAssayData(object = object, assay.type = assay.type)),
      k = ncenters, samples = nsamples
    )
    #identify positive and negative signals for all HTO
    object <- SetIdent(
      object = object,
      cells.use = names(x = hto_init_clusters$clustering),
      ident.use = hto_init_clusters$clustering
    )
  }
  #work around so we don't average all the RNA levels which takes time
  object2 <- object
  object2@data <- object2@data[1:10, ]
  hto_averages <- AverageExpression(
    object = object2,
    return.seurat = TRUE,
    show.progress = FALSE
  )
  average_hto <- GetAssayData(
    object = hto_averages,
    assay.type = assay.type,
    slot = "raw.data"
  )
  hto_discrete <- GetAssayData(object = object, assay.type = assay.type)
  hto_discrete[hto_discrete > 0] <- 0
  hto_prob_pos <- hto_discrete
  rownames(x = hto_prob_pos) <- paste(rownames(x = hto_prob_pos), "pos", sep = "_")
  hto_prob_neg <- hto_discrete
  rownames(x = hto_prob_neg) <- paste(rownames(x = hto_prob_neg), "neg", sep = "_")
  hash_data <- GetAssayData(object = object, assay.type = assay.type)
  hash_raw_data <- GetAssayData(
    object = object,
    assay.type = assay.type,
    slot = "raw.data"
  )
  # for each HTO, we will use the minimum cluster for fitting
  # we will also store positive and negative distributions for all cells after determining the cutoff
  hto_dist_pos = list()
  hto_dist_neg = list()
  for (hto_iter in rownames(x = hash_data)) {
    hto_values <- hash_raw_data[hto_iter, object@cell.names]
    #commented out if we take all but the top cluster as background
    #hto_values_negative=hto_values[setdiff(object@cell.names,WhichCells(object,which.max(average_hto[hto_iter,])))]
    hto_values_use <- hto_values[WhichCells(object = object, ident = which.min(x = average_hto[hto_iter, ])
    )]
    #throw off the top 0% of values, probably should be excluded for the min
    #hto_values_use <- hto_values_negative[hto_values_negative<quantile(hto_values_negative,1)]
    hto_fit <- fitdist(hto_values_use, "nbinom")
    hto_cutoff <- as.numeric(x = quantile(x = hto_fit, probs = percent_cutoff)$quantiles[1])
    hto_discrete[hto_iter, names(x = which(x = hto_values > hto_cutoff))] <- 1
    if (print.output) {
      print(paste0("Cutoff for ", hto_iter, " : ", hto_cutoff, " reads"))
    }
    hto_values_pos <- hto_values[hto_values > hto_cutoff]
    hto_values_neg <- hto_values[hto_values <= hto_cutoff]
    hto_dist_neg[[hto_iter]] <- fitdist(hto_values_neg, "nbinom")
    hto_dist_pos[[hto_iter]] <- fitdist(hto_values_pos, "nbinom")
    hto_prob_pos[paste(hto_iter, "pos", sep = "_"),] <- sapply(
      X = hto_values,
      FUN = function(x) {
        return(pnbinom(
          q = x,
          size = hto_dist_pos[[hto_iter]]$estimate["size"],
          mu = hto_dist_pos[[hto_iter]]$estimate["mu"],
          log.p = TRUE
        ))
      }
    )
    hto_prob_neg[paste(hto_iter, "neg", sep = "_"),] <- -1 * sapply(
      X = hto_values,
      FUN = function(x) {
        return(pnbinom(
          q = x,
          size = hto_dist_neg[[hto_iter]]$estimate["size"],
          mu = hto_dist_neg[[hto_iter]]$estimate["mu"],
          log.p = TRUE,
          lower.tail = FALSE
        ))
      }
    )
    #hto_prob_neg=sapply(hto_values,function(x) pnbinom(x,size = hto_dist_neg[[hto_iter]]$estimate["size"],mu = hto_dist_neg[[hto_iter]]$estimate["mu"]))
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
  # hash_second <- apply(X = hash_data, MARGIN = 2, FUN = function(x) MaxN(x, 2))
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
  doublet_id <- sapply(X = 1:length(x = hash_maxID), function(x) paste(sort(c(hash_maxID[x], hash_secondID[x])), collapse = "_"))
  # doublet_names <- names(x = table(doublet_id))[-1] # Not used
  hto_classification <- hto_classification_global
  hto_classification[hto_classification_global == "Negative"] <- "Negative"
  hto_classification[hto_classification_global == "Singlet"] <- hash_maxID[which(hto_classification_global == "Singlet")]
  hto_classification[hto_classification_global == "Doublet"] <- doublet_id[which(hto_classification_global == "Doublet")]
  classification_metadata <- data.frame(
    hash_max,
    hash_maxID,
    hash_second,
    hash_secondID, hash_margin,
    hto_classification,
    hto_classification_global
  )
  object <- AddMetaData(object = object, metadata = classification_metadata)
  hto_shortID <- paste(object@meta.data$hash_maxID, object@meta.data$hto_classification_global, sep = "_")
  names(x = hto_shortID) <- object@cell.names
  object <- AddMetaData(object = object, metadata = hto_shortID, col.name = "hto_shortID")
  object <- SetAllIdent(object = object, id = "hto_shortID")
  object_pn <- data.frame(t(x = rbind(hto_prob_neg, hto_prob_pos)))
  object <- AddMetaData(object = object, metadata = object_pn)
  confidence_threshold <- log(x = confidence_threshold)
  object@meta.data$confidence <- "Confident"
  singlet_data <- subset(x = object@meta.data, subset = hto_classification_global == "Singlet")
  singlet_pos <- sapply(
    X = 1:nrow(x = singlet_data),
    FUN = function(x) {
      return(singlet_data[
        x,
        paste(as.character(x = singlet_data[x, "hash_maxID"]), "pos", sep = "_")
      ])
    }
  )
  names(x = singlet_pos) <- rownames(x = singlet_data)
  object@meta.data[names(x = which(x = singlet_pos < confidence_threshold)), "confidence"] <- "Uncertain"
  doublet_data <- subset(x = object@meta.data, hto_classification_global == "Doublet")
  doublet_pos <- sapply(
    X = 1:nrow(x = doublet_data),
    FUN = function(x) {
      return(doublet_data[
        x,
        paste(as.character(x = doublet_data[x, "hash_maxID"]), "pos", sep = "_")
      ])
    }
  )
  names(x = doublet_pos) <- rownames(x = doublet_data)
  object@meta.data[names(x = which(x = doublet_pos < confidence_threshold)), "confidence"] <- "Uncertain"
  doublet_pos <- sapply(
    X = 1:nrow(doublet_data),
    FUN = function(x) {
      return(doublet_data[
        x,
        paste(as.character(x = doublet_data[x, "hash_secondID"]), "pos", sep = "_")
      ])
    }
  )
  names(x = doublet_pos) <- rownames(x = doublet_data)
  object@meta.data[names(x = which(x = doublet_pos < confidence_threshold)), "confidence"] <- "Uncertain"
  if (print.output) {
    print(x = table(object@meta.data$hto_classification_global, object@meta.data$confidence))
  }
  return(object)
}
