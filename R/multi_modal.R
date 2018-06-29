#' @include seurat.R
#' @importFrom methods setClass setMethod
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
    var.genes = "vector",
    mean.var = "data.frame"
  )
)

setMethod(
  f = 'show',
  signature = 'assay',
  definition = function(object) {
    cat(
      'Seurat assay data with',
      nrow(x = object@data),
      'measurements for',
      ncol(x = object@data), 'cells\n'
    )
    if (length(x = object@var.genes) > 0) {
      cat(
        "Top 10 variable measurements:\n",
        strwrap(x = paste(head(x = object@var.genes, n = 10L), collapse = ', '))
      )
    }
  }
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
#' @importFrom methods slotNames
#'
#' @export
#'
#' @examples
#' # Simulate CITE-Seq results
#' df <- t(x = data.frame(
#'   x = round(x = rnorm(n = 80, mean = 20, sd = 2)),
#'   y = round(x = rbinom(n = 80, size = 100, prob = 0.2)),
#'   row.names = pbmc_small@cell.names
#' ))
#' pbmc_small <- SetAssayData(
#'   object = pbmc_small,
#'   assay.type = 'CITE',
#'   new.data = df,
#'   slot = 'data'
#' )
#' GetAssayData(object = pbmc_small, assay.type = 'CITE', slot = 'data')
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
    if (is.null(x = to.return)) {
      return(to.return)
    }
    #note that we check for this to avoid a long subset for large matrices if it can be avoided
    if (length(x = object@cell.names) == ncol(to.return)) {
      return(to.return)
    }
    return(to.return[, object@cell.names])
  }
  if (!(assay.type %in% names(object@assay))) {
    stop(paste(assay.type, "data has not been added"))
  }
  if (!(slot %in% slotNames(eval(expr = parse(text = paste0("object@assay$", assay.type)))))) {
    stop(paste(slot, "slot doesn't exist"))
  }
  to.return <- slot(object = object@assay[[assay.type]], name = slot)
  if(is.null(to.return)){
    stop(paste0("The ", slot, " for the ", assay.type, " assay is empty"))
  }
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
#' @importFrom methods new slot<-
#'
#' @export
#'
#' @examples
#' # Simulate CITE-Seq results
#' df <- t(x = data.frame(
#'   x = round(x = rnorm(n = 80, mean = 20, sd = 2)),
#'   y = round(x = rbinom(n = 80, size = 100, prob = 0.2)),
#'   row.names = pbmc_small@cell.names
#' ))
#' pbmc_small <- SetAssayData(
#'   object = pbmc_small,
#'   assay.type = 'CITE',
#'   new.data = df,
#'   slot = 'data'
#' )
#' pbmc_small@assay
#'
SetAssayData <- function(object, assay.type, slot, new.data) {
  if (ncol(x = new.data) != length(x = object@cell.names)) {
    stop("Wrong number of cells being added: should be ", length(x = object@cell.names))
  } else if (!all(colnames(x = new.data) %in% object@cell.names)) {
    stop("Cell names in new assay data matrix don't exactly match the cell names of the object")
  }
  new.data <- new.data[, match(x = object@cell.names, table = colnames(x = new.data))]
  if (any(colnames(x = new.data) != object@cell.names)) {
    stop("Please ensure cell names for assay data match exactly (including in order) the cell names in object@cell.names")
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
    slot(object = new.assay, name = slot) <- new.data
    object@assay <- c(object@assay, new.assay)
    names(x = object@assay) <- c(
      names(x = object@assay)[1:(length(x = object@assay) - 1)],
      assay.type
    )
    # eval(expr = parse(text = paste0("new.assay@", slot, "<- new.data")))
    # eval(expr = parse(text = paste0("object@assay$", assay.type, "<- new.assay")))
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
#' @param object Seurat object. Assumes that the hash tag oligo (HTO) data has been added and normalized.
#' @param assay.type Name of the Hashtag assay (HTO by default)
#' @param positive_quantile The quantile of inferred 'negative' distribution for each hashtag - over which the cell is considered 'positive'. Default is 0.99
#' @param init_centers Initial number of clusters for hashtags. Default is the # of hashtag oligo names + 1 (to account for negatives)
#' @param k_function Clustering function for initial hashtag grouping. Default is "clara" for fast k-medoids clustering on large applications, also support "kmeans" for kmeans clustering
#' @param nsamples Number of samples to be drawn from the dataset used for clustering, for k_function = "clara"
#' @param cluster_nstarts nstarts value for k-means clustering (for k_function = "kmeans"). 100 by default
#' @param print.output Prints the output
#' 
#' @return The Seurat object with the following demultiplexed information stored in the meta data:
#' \item{hash_maxID}{Name of hashtag with the highest signal}
#' \item{hash_secondID}{Name of hashtag with the second highest signal}
#' \item{hash_margin}{The difference between signals for hash_maxID and hash_secondID}
#' \item{hto_classification}{Classification result, with doublets/multiplets named by the top two highest hashtags}
#' \item{hto_classification_global}{Global classification result (singlet, doublet or negative)}
#' \item{hash_ID}{Classification result where doublet IDs are collapsed}
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
  assay.type = "HTO",
  positive_quantile = 0.99,
  init_centers = NULL,
  cluster_nstarts = 100,
  k_function = "clara",
  nsamples = 100,
  print.output = TRUE
) {
  #initial clustering
  hash_data <- GetAssayData(object = object, assay.type = assay.type)
  hash_raw_data <- GetAssayData(
    object = object,
    assay.type = assay.type,
    slot = "raw.data"
  )[, object@cell.names]
  hash_raw_data <- as.matrix(hash_raw_data)
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
  
  #average hto signals per cluster
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
  
  #create a matrix to store classification result
  hto_discrete <- GetAssayData(object = object, assay.type = assay.type)
  hto_discrete[hto_discrete > 0] <- 0
  
  # for each HTO, we will use the minimum cluster for fitting
  for (hto_iter in rownames(x = hash_data)) {
    hto_values <- hash_raw_data[hto_iter, object@cell.names]
    #commented out if we take all but the top cluster as background
    #hto_values_negative=hto_values[setdiff(object@cell.names,WhichCells(object,which.max(average_hto[hto_iter,])))]
    hto_values_use <- hto_values[WhichCells(object = object, ident = which.min(x = average_hto[hto_iter, ])
    )]
    
    hto_fit <- fitdist(hto_values_use, "nbinom")
    hto_cutoff <- as.numeric(x = quantile(x = hto_fit, probs = positive_quantile)$quantiles[1])
    hto_discrete[hto_iter, names(x = which(x = hto_values > hto_cutoff))] <- 1
    if (print.output) {
      print(paste0("Cutoff for ", hto_iter, " : ", hto_cutoff, " reads"))
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
  doublet_id <- sapply(X = 1:length(x = hash_maxID), function(x) paste(sort(c(hash_maxID[x], hash_secondID[x])), collapse = "_"))
  # doublet_names <- names(x = table(doublet_id))[-1] # Not used
  hto_classification <- hto_classification_global
  hto_classification[hto_classification_global == "Negative"] <- "Negative"
  hto_classification[hto_classification_global == "Singlet"] <- hash_maxID[which(hto_classification_global == "Singlet")]
  hto_classification[hto_classification_global == "Doublet"] <- doublet_id[which(hto_classification_global == "Doublet")]
  classification_metadata <- data.frame(
    hash_maxID,
    hash_secondID, hash_margin,
    hto_classification,
    hto_classification_global
  )
  object <- AddMetaData(object = object, metadata = classification_metadata)
  
  if (print.output) {
    print(x = table(object@meta.data$hto_classification_global))
  }
  object=SetAllIdent(object = object,id = "hto_classification")
  object=SetIdent(object,cells.use = WhichCells(object,subset.name = "hto_classification_global",accept.value = "Doublet"),ident.use = "Doublet")
  object@meta.data$hash_ID=object@ident[rownames(object@meta.data)]
  
  return(object)
}



#' Hashtag oligo heatmap
#'
#' Draws a heatmap of hashtag oligo signals across singlets/doublets/negative cells. Allows for the visualization of HTO demultiplexing results.
#'
#' @param object Seurat object. Assumes that the hash tag oligo (HTO) data has been added and normalized, and demultiplexing has been run with HTODemux().
#' @param hto.classification The naming for object@meta.data slot with classification result from HTODemux().
#' @param global.classification The slot for object@meta.data slot specifying a cell as singlet/doublet/negative.
#' @param assay.type Hashtag assay name.
#' @param num.cells Number of cells to plot. Default is to choose 5000 cells by random subsampling, to avoid having to draw exceptionally large heatmaps.
#' @param singlet.names Namings for the singlets. Default is to use the same names as HTOs.
#' @param ... Additional arguments for DoHeatmap().
#' 
#' @return Returns a ggplot2 plot object.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' object <- HTODemux(object)
#' HTOHeatmap(object)
#' }

HTOHeatmap <- function(
  object,
  hto.classification = "hto_classification",
  global.classification = "hto_classification_global",
  assay.type = "HTO", 
  num.cells = 5000, 
  singlet.names = NULL,
  ...
){
  
  object <- SetAllIdent(object,id = hto.classification)
  objmini <- SubsetData(object,cells.use = sample(object@cell.names,num.cells))
  
  metadata_use <- objmini@meta.data
  singlet_id <- sort(unique(as.character(metadata_use[metadata_use[,global.classification]=="Singlet",hto.classification])))
  doublet_id <- sort(unique(as.character(metadata_use[metadata_use[,global.classification]=="Doublet",hto.classification])))
  
  heatmap_levels <- c(singlet_id,doublet_id,"Negative")
  objmini <- SetIdent(objmini,FastWhichCells(objmini,group.by = global.classification,"Doublet"),"Multiplet")
  
  objmini@ident <- factor(objmini@ident, c(singlet_id,"Multiplet","Negative"))
  cells.ordered=as.character(unlist(sapply(heatmap_levels,function(x) sample(FastWhichCells(objmini,group.by = hto.classification,x)))))
  objmini <- ScaleData(objmini,assay.type = assay.type,display.progress = FALSE)
  
  if (!is.null(singlet.names)){
    levels(objmini@ident) <- c(singlet.names, "Multiplet", "Negative")
  } 
  DoHeatmap(objmini,slim.col.label = T,genes.use = singlet_id,assay.type = assay.type,cells.use = cells.ordered,group.label.rot = T)
  
}