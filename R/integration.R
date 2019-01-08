#' @include generics.R
#'
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Find integration anchors
#'
#' Finds the integration anchors
#'
#' @param object.list A list of objects between which to find anchors for downstream integration.
#' @param assay A vector of assay names specifying which assay to use when constructing anchors. If
#' NULL, the current default assay for each object is used.
#' @param anchor.features Can be either:
#' \itemize{
#'   \item{A numeric value. This will call \code{\link{SelectIntegrationFeatures}} to select the
#'   provided number of features to be used in anchor finding}
#'   \item{A vector of features to be used as input to the anchor finding process}
#' }
#' @param scale Whether or not to scale the features provided. Only set to FALSE if you have
#' previously scaled the features you want to use for each object in the object.list
#' @param l2.norm Perform L2 normalization on the CCA cell embeddings after dimensional reduction
#' @param dims Which dimensions to use from the CCA to specify the neighbor search space
#' @param k.anchor How many neighbors (k) to use when picking anchors
#' @param k.filter How many neighbors (k) to use when filtering anchors
#' @param k.score How many neighbors (k) to use when scoring anchors
#' @param max.features The maximum number of features to use when specifying the neighborhood search
#' space in the anchor filtering
#' @param eps Error bound on the neighbor finding algorithm (from RANN)
#' @param verbose Print progress bars and output
#'
#' @return Returns an AnchorSet object
#'
#' @importFrom pbapply pblapply
#' @importFrom future.apply future_lapply
#'
#' @export
#'
FindIntegrationAnchors <- function(
  object.list = NULL,
  assay = NULL,
  anchor.features = 2000,
  scale = TRUE,
  l2.norm = TRUE,
  dims = 1:30,
  k.anchor = 5,
  k.filter = 200,
  k.score = 30,
  max.features = 200,
  eps = 0,
  verbose = TRUE
) {
  my.lapply <- ifelse(
    test = verbose && PlanThreads() == 1,
    yes = pblapply,
    no = future_lapply
  )
  if (!is.null(x = assay)) {
    if (length(x = assay) != length(x = object.list)) {
      stop("If specifying the assay, please specify one assay per object in the object.list")
    }
    for (ii in length(x = object.list)) {
      DefaultAssay(object = object.list[[ii]]) <- assay[ii]
    }
  } else {
    assay <- sapply(X = object.list, FUN = DefaultAssay)
  }
  object.list <- CheckDuplicateCellNames(object.list = object.list)
  if (is.numeric(x = anchor.features)) {
    if (verbose) {
      message(paste("Computing", anchor.features, "integration features"))
    }
    anchor.features <- SelectIntegrationFeatures(
      object.list = object.list,
      nfeatures = anchor.features,
      assay = assay
    )
  }
  if (scale) {
    if (verbose) {
      message("Scaling features for provided objects")
    }
    object.list <- my.lapply(
      X = object.list,
      FUN = function(object) {
        ScaleData(object = object, features = anchor.features, verbose = FALSE)
      }
    )
  }
  # determine pairwise combinations
  combinations <- expand.grid(1:length(x = object.list), 1:length(x = object.list))
  combinations <- combinations[combinations$Var1 < combinations$Var2, , drop = FALSE]
  # determine the proper offsets for indexing anchors
  objects.ncell <- sapply(X = object.list, FUN = ncol)
  offsets <- as.vector(x = cumsum(x = c(0, objects.ncell)))[1:length(x = object.list)]
  if (verbose) {
    message("Finding all pairwise anchors")
  }
  reduction <- "cca"
  # determine all pairwise anchors
  all.anchors <- my.lapply(
    X = 1:nrow(x = combinations),
    FUN = function(row) {
      i <- combinations[row, 1]
      j <- combinations[row, 2]
      object.1 <- object.list[[i]]
      object.2 <- object.list[[j]]
      object.pair <- RunCCA(
        object1 = object.1,
        object2 = object.2,
        features = anchor.features,
        num.cc = max(dims),
        renormalize = FALSE,
        rescale = FALSE,
        verbose = verbose
      )
      if (l2.norm){
        object.pair <- L2Dim(object = object.pair, reduction = reduction)
        reduction <- paste0(reduction, ".l2")
      }
      anchors <- FindAnchors(
        object.pair = object.pair,
        assay = c(assay[i], assay[j]),
        cells1 = colnames(x = object.1),
        cells2 = colnames(x = object.2),
        reduction = reduction,
        dims = dims,
        k.anchor = k.anchor,
        k.filter = k.filter,
        k.score = k.score,
        max.features = max.features,
        eps = eps,
        verbose = verbose
      )
      anchors[, 1] <- anchors[, 1] + offsets[i]
      anchors[, 2] <- anchors[, 2] + offsets[j]
      return(anchors)
    }
  )
  all.anchors <- do.call(what = 'rbind', args = all.anchors)
  all.anchors <- rbind(all.anchors, all.anchors[, c(2, 1, 3)])
  all.anchors <- AddDatasetID(anchor.df = all.anchors, offsets = offsets, obj.lengths = objects.ncell)
  anchor.set <- new(Class = "AnchorSet",
                    object.list = object.list,
                    anchors = all.anchors,
                    offsets = offsets,
                    anchor.features = anchor.features
  )
  return(anchor.set)
}

#' Find transfer anchors
#'
#' Finds the transfer anchors
#'
#' @param reference Seurat object to use as the reference
#' @param query Seurat object to use as the query
#' @param reference.assay Assay to use from reference
#' @param query.assay Assay to use from query
#' @param reduction Dimensional reduction to perform when finding anchors. Options are:
#' \itemize{
#'    \item{pcaproject: Project the PCA from the reference onto the query. We recommend using PCA
#'    when reference and query datasets are from scRNA-seq}
#'    \item{cca: Run a CCA on the reference and query }
#' }
#' @param project.query Project the PCA from the query dataset onto the reference. Use only in rare
#' cases where the query dataset has a much larger cell number, but the reference dataset has a
#' unique assay for transfer.
#' @param features Features to use for dimensional reduction
#' @param npcs Number of PCs to compute on reference. If null, then use an existing PCA structure in
#' the reference object
#' @param l2.norm Perform L2 normalization on the cell embeddings after dimensional reduction
#' @param dims Which dimensions to use from the reduction to specify the neighbor search space
#' @param k.anchor How many neighbors (k) to use when picking anchors
#' @param k.filter How many neighbors (k) to use when filtering anchors
#' @param k.score How many neighbors (k) to use when scoring anchors
#' @param max.features The maximum number of features to use when specifying the neighborhood search
#' space in the anchor filtering
#' @param eps Error bound on the neighbor finding algorithm (from RANN)
#' @param approx.pca Use truncated singular value decomposition to approximate PCA
#' @param verbose Print progress bars and output
#'
#' @return Returns an AnchorSet object
#'
#'
#' @export
#'
FindTransferAnchors <- function(
  reference,
  query,
  reference.assay = NULL,
  query.assay = NULL,
  reduction = "pcaproject",
  project.query = FALSE,
  features = NULL,
  npcs = 30,
  l2.norm = TRUE,
  dims = 1:30,
  k.anchor = 5,
  k.filter = 200,
  k.score = 30,
  max.features = 200,
  eps = 0,
  approx.pca = TRUE,
  verbose = TRUE
) {
  if (length(x = reference) > 1 | length(x = query) > 1) {
    stop("We currently only support transfer between a single query and reference")
  }
  if (!reduction %in% c("pcaproject", "cca", "pcaqueryproject")) {
    stop("Please select either pcaproject, cca, or pcaqueryproject for the reduction parameter.")
  }
  query <- RenameCells(object = query, new.names = paste0(Cells(object = query), "_", "query"))
  reference <- RenameCells(object = reference, new.names = paste0(Cells(object = reference), "_", "reference"))
  features <- features %||% VariableFeatures(object = reference)
  reference.assay <- reference.assay %||% DefaultAssay(object = reference)
  query.assay <- query.assay %||% DefaultAssay(object = query)
  DefaultAssay(object = reference) <- reference.assay
  DefaultAssay(object = query) <- query.assay

  ## find anchors using PCA projection
  if (reduction == 'pcaproject') {
    if (project.query){
      if (!is.null(x = npcs)) {
        if (verbose) {
          message("Performing PCA on the provided query using ", length(x = features), " features as input.")
        }
        query <- ScaleData(object = query, features = features, verbose = FALSE)
        query <- RunPCA(object = query, npcs = npcs, verbose = FALSE, features = features, approx = approx.pca)
      }
      projected.pca <- ProjectCellEmbeddings(
        reference = query,
        query = reference,
        dims = dims,
        verbose = verbose
      )
      query.pca <- Embeddings(object = query[["pca"]])[, dims]
      combined.pca <- CreateDimReducObject(
        embeddings = as.matrix(x = rbind(projected.pca,query.pca))[, dims],
        key = "ProjectPC_",
        assay = reference.assay
      )
      combined.ob <- merge(x = reference, y = query)
      combined.ob[["pcaproject"]] <- combined.pca
      old.loadings <- Loadings(object = query[["pca"]])
      colnames(x = old.loadings) <- paste0("ProjectPC_", 1:ncol(x = old.loadings))
      Loadings(object = combined.ob[["pcaproject"]]) <- old.loadings[, dims]
    } else {
      if (!is.null(x = npcs)) {
        if (verbose) {
          message("Performing PCA on the provided reference using ", length(x = features), " features as input.")
        }
        reference <- ScaleData(object = reference, features = features, verbose = FALSE)
        reference <- RunPCA(object = reference, npcs = npcs, verbose = FALSE,features = features, approx = approx.pca)
      }
      projected.pca <- ProjectCellEmbeddings(
        reference = reference,
        query = query,
        dims = dims,
        verbose = verbose
      )
      ref.pca <- Embeddings(object = reference[["pca"]])[, dims]
      combined.pca <- CreateDimReducObject(
        embeddings = as.matrix(x = rbind(ref.pca, projected.pca))[, dims],
        key = "ProjectPC_",
        assay = reference.assay
      )
      combined.ob <- merge(x = reference, y = query)
      combined.ob[["pcaproject"]] <- combined.pca
      old.loadings <- Loadings(object = reference[["pca"]])
      colnames(x = old.loadings) <- paste0("ProjectPC_", 1:ncol(x = old.loadings))
      Loadings(object = combined.ob[["pcaproject"]]) <- old.loadings[, dims]
    }
  }

  ## find anchors using CCA
  if (reduction == 'cca') {
    reference <- ScaleData(object = reference, features = features, verbose = FALSE)
    query <- ScaleData(object = query, features = features, verbose = FALSE)
    combined.ob <- RunCCA(
      object1 = reference,
      object2 = query,
      features = features,
      num.cc = max(dims),
      renormalize = FALSE,
      rescale = FALSE,
      verbose = verbose
    )
  }

  if (l2.norm){
    combined.ob <- L2Dim(object = combined.ob, reduction = reduction)
    reduction <- paste0(reduction, ".l2")
  }
  anchors <- FindAnchors(
    object.pair = combined.ob,
    assay = c(reference.assay, query.assay),
    cells1 = colnames(x = reference),
    cells2 = colnames(x = query),
    reduction = reduction,
    dims = dims,
    k.anchor = k.anchor,
    k.filter = k.filter,
    k.score = k.score,
    max.features = max.features,
    eps = eps,
    verbose = verbose
  )
  anchor.set <- new(
    Class = "AnchorSet",
    object.list = list(combined.ob),
    reference.cells = colnames(x = reference),
    query.cells = colnames(x = query),
    anchors = anchors,
    anchor.features = features
  )
  return(anchor.set)
}


#' Integrate data
#'
#' Integrates the data
#'
#' @param anchorset Results from FindIntegrationAnchors
#' @param new.assay.name Name for the new assay containing the integrated data
#' @param features Vector of features to use when computing the PCA to determine the weights. Only set
#' if you want a different set from those used in the anchor finding process
#' @param features.to.integrate Vector of features to integrate. By default, will use the features
#' used in anchor finding.
#' @param dims Number of PCs to use in the weighting procedure
#' @param k.weight Number of neighbors to consider when weighting
#' @param sd.weight Controls the bandwidth of the Gaussian kernel for weighting
#' @param sample.tree Specify the order of integration. If null, will compute automatically.
#' @param preserve.order Integrate the objects in the order they were provided to FindAnchors
#' @param do.cpp Run cpp code where applicable
#' @param eps Error bound on the neighbor finding algorithm (from RANN)
#' @param verbose Print progress bars and output
#'
#' @return Returns a Seurat object with a new integrated Assay
#'
#' @export
#'
IntegrateData <- function(
  anchorset,
  new.assay.name = "integrated",
  features = NULL,
  features.to.integrate = NULL,
  dims = 1:30,
  k.weight = 100,
  sd.weight = 1,
  sample.tree = NULL,
  preserve.order = FALSE,
  do.cpp = TRUE,
  eps = 0,
  verbose = TRUE
) {
  object.list <- slot(object = anchorset, name = "object.list")
  anchors <- slot(object = anchorset, name = "anchors")
  offsets <- slot(object = anchorset, name = "offsets")
  features <- features %||% slot(object = anchorset, name = "anchor.features")
  features.to.integrate <- features.to.integrate %||% features
  objects.ncell <- sapply(X = object.list, FUN = ncol) #function(x) ncol(x = x))
  if (is.null(x = sample.tree)) {
    similarity.matrix <- CountAnchors(
      anchor.df = anchors,
      offsets = offsets,
      obj.lengths = objects.ncell
    )
    sample.tree <- BuildSampleTree(similarity.matrix = similarity.matrix)
  }
  cellnames.list <- list()
  for (ii in 1:length(x = object.list)) {
    cellnames.list[[ii]] <- colnames(x = object.list[[ii]])
  }
  unintegrated <- merge(
    x = object.list[[1]],
    y = object.list[2:length(x = object.list)]
  )
  names(x = object.list) <- as.character(-(1:length(x = object.list)))
  for (ii in 1:nrow(x = sample.tree)) {
    merge.pair <- as.character(x = sample.tree[ii, ])
    length1 <- ncol(x = object.list[[merge.pair[1]]])
    length2 <- ncol(x = object.list[[merge.pair[2]]])
    if (!(preserve.order) & (length2 > length1)) {
      merge.pair <- rev(merge.pair)
      sample.tree[ii, ] <- as.numeric(merge.pair)
    }
    object.1 <- object.list[[merge.pair[1]]]
    object.2 <- object.list[[merge.pair[2]]]
    datasets <- ParseMergePair(sample.tree, ii)
    if (verbose) {
      message(
        "Merging dataset ",
        paste(datasets$object2, collapse = " "),
        " into ",
        paste(datasets$object1, collapse = " ")
      )
    }
    cells1 <- colnames(x = object.1)
    cells2 <- colnames(x = object.2)
    merged.obj <- merge(x = object.1, y = object.2, merge.data = TRUE)
    if (verbose) {
      message("Extracting anchors for merged samples")
    }
    filtered.anchors <- anchors[anchors$dataset1 %in% datasets$object1 & anchors$dataset2 %in% datasets$object2, ]
    cell1.name <- sapply(X = 1:nrow(x = filtered.anchors), FUN = function(x) {
      cellnames.list[[filtered.anchors$dataset1[x]]][filtered.anchors$cell1[x]]
    })
    cell2.name <- sapply(X = 1:nrow(x = filtered.anchors), FUN = function(x) {
      cellnames.list[[filtered.anchors$dataset2[x]]][filtered.anchors$cell2[x]]
    })
    cell1.offset <- sapply(
      X = 1:length(x = cell1.name),
      FUN = function(x) {
        return(which(cells1 == cell1.name[x]))
      }
    )
    cell2.offset <- sapply(
      X = 1:length(x = cell2.name),
      FUN = function(x) {
        return(which(cells2 == cell2.name[x]))
      }
    )
    filtered.anchors.new <- filtered.anchors
    filtered.anchors.new[, 1] <- cell1.offset
    filtered.anchors.new[, 2] <- cell2.offset
    integration.name <- "integrated"
    merged.obj <- SetIntegrationData(
      object = merged.obj,
      integration.name = integration.name,
      slot = 'anchors',
      new.data = filtered.anchors.new
    )
    merged.obj <- SetIntegrationData(
      object = merged.obj,
      integration.name = integration.name,
      slot = 'neighbors',
      new.data = list('cells1' = cells1, 'cells2' = cells2)
    )
    merged.obj <- FindIntegrationMatrix(
      object = merged.obj,
      integration.name = integration.name,
      features.integrate = features.to.integrate,
      verbose = verbose
    )
    assay <- DefaultAssay(object = merged.obj)
    merged.obj <- ScaleData(
      object = merged.obj,
      features = features,
      verbose = FALSE
    )
    merged.obj <- RunPCA(
      object = merged.obj,
      npcs = max(dims),
      verbose = FALSE,
      features = features
    )
    merged.obj <- FindWeights(
      object = merged.obj,
      integration.name = integration.name,
      reduction = 'pca',
      nn.reduction = 'pca',
      do.cpp = do.cpp,
      nn.dims = dims,
      k.weight = k.weight,
      sd.weight = sd.weight,
      eps = eps,
      verbose = verbose
    )
    merged.obj <- TransformDataMatrix(
      object = merged.obj,
      new.assay.name = new.assay.name,
      features.to.integrate = features.to.integrate,
      integration.name = integration.name,
      do.cpp = do.cpp,
      verbose = verbose
    )
    integrated.matrix <- GetAssayData(
      object = merged.obj,
      assay = new.assay.name,
      slot = 'data'
    )
    merged.obj <- SetAssayData(
      object = merged.obj,
      assay = assay,
      slot = 'data',
      new.data = integrated.matrix
    )
    object.list[[as.character(x = ii)]] <- merged.obj
    object.list[[merge.pair[[1]]]] <- NULL
    object.list[[merge.pair[[2]]]] <- NULL
    invisible(x = gc(verbose = FALSE))
  }
  integrated.data <- GetAssayData(
    object = object.list[[as.character(x = ii)]],
    assay = assay,
    slot = 'data'
  )
  integrated.data <- integrated.data[, colnames(x = unintegrated)]
  new.assay <- new(
    Class = 'Assay',
    counts =  new(Class = "dgCMatrix"),
    data = integrated.data,
    scale.data = matrix(),
    var.features = vector(),
    meta.features = data.frame(row.names = rownames(x = integrated.data)),
    misc = NULL
  )
  unintegrated[[new.assay.name]] <- new.assay
  # "unintegrated" now contains the integrated assay
  DefaultAssay(object = unintegrated) <- new.assay.name
  VariableFeatures(object = unintegrated) <- features
  unintegrated <- SetIntegrationData(
    object = unintegrated,
    integration.name = "Integration",
    slot = "anchors",
    new.data = anchors
  )
  unintegrated <- SetIntegrationData(
    object = unintegrated,
    integration.name = "Integration",
    slot = "sample.tree",
    new.data = sample.tree
  )
  return(unintegrated)
}

#' Calculate the local structure preservation metric
#'
#' Calculates a metric that describes how well the local structure of each group
#' prior to integration is preserved after integration. This procedure works as
#' follows: For each group, compute a PCA, compute the top num.neighbors in pca
#' space, compute the top num.neighbors in corrected pca space, compute the
#' size of the intersection of those two sets of neighbors.
#' Return the average over all groups.
#'
#' @param object Seurat object
#' @param grouping.var Grouping variable
#' @param idents Optionally specify a set of idents to compute metric for
#' @param neighbors Number of neighbors to compute in pca/corrected pca space
#' @param reduction Dimensional reduction to use for corrected space
#' @param reduced.dims Number of reduced dimensions to use
#' @param orig.dims Number of PCs to use in original space
#' @param verbose Display progress bar
#'
#' @return Returns the average preservation metric
#'
#' @importFrom RANN nn2
#' @importFrom utils txtProgressBar setTxtProgressBar
#'
#' @export
#'
LocalStruct <- function(
  object,
  grouping.var,
  idents = NULL,
  neighbors = 100,
  reduction = "pca",
  reduced.dims = 1:10,
  orig.dims = 1:10,
  verbose = TRUE
) {
  if (is.null(x = idents)) {
    cells.use <- colnames(x = object)
  } else {
    cells.use <- WhichCells(object = object, idents = idents)
  }
  Idents(object = object) <- grouping.var
  local.struct <- list()
  ob.list <- SplitObject(object = object, split.by = grouping.var)
  if (verbose) {
    pb <- txtProgressBar(
      min = 1,
      max = length(x = ob.list),
      style = 3,
      file = stderr()
    )
  }
  embeddings <- Embeddings(object = object[[reduction]])[, reduced.dims]

  for (i in 1:length(x = ob.list)) {
    ob <- ob.list[[i]]
    ob <- FindVariableFeatures(
      object = ob,
      verbose = FALSE,
      selection.method = "dispersion",
      nfeatures = 2000
    )
    ob <- ScaleData(
      object = ob,
      features = VariableFeatures(object = ob),
      verbose = FALSE
    )
    ob <- RunPCA(
      object = ob,
      features = VariableFeatures(object = ob),
      verbose = FALSE,
      npcs = max(orig.dims)
    )
    ob.cells <- intersect(x = cells.use, y = colnames(x = ob))
    if (length(x = ob.cells) == 0) next
    nn.corrected <- nn2(
      data = embeddings[colnames(x = ob), ],
      query = embeddings[ob.cells, ],
      k = neighbors
    )$nn.idx
    nn.orig <- nn2(
      data = Embeddings(object = ob[["pca"]])[, orig.dims],
      query = Embeddings(object = ob[["pca"]])[ob.cells, orig.dims],
      k = neighbors
    )$nn.idx
    local.struct[[i]] <- sapply(X = 1:nrow(x = nn.orig), FUN = function(x) {
      length(x = intersect(x = nn.orig[x, ], y = nn.corrected[x, ])) / neighbors
    })
    if (verbose) {
      setTxtProgressBar(pb = pb, value = i)
    }
  }
  names(x = local.struct) <- names(x = ob.list)
  return(local.struct)
}

#' Calculates a mixing metric
#'
#' Here we compute a measure of how well mixed a composite dataset is. To compute, we first examine
#' the local neighborhood for each cell (looking at max.k neighbors) and determine for each group
#' (could be the dataset after integration) the k nearest neighbor and what rank that neighbor was
#' in the overall neighborhood. We then take the median across all groups as the mixing metric per
#' cell.
#'
#' @param object Seurat object
#' @param grouping.var Grouping variable for dataset
#' @param reduction Which dimensionally reduced space to use
#' @param dims Dimensions to use
#' @param k Neighbor number to examine per group
#' @param max.k Maximum size of local neighborhood to compute
#' @param eps Error bound on the neighbor finding algorithm (from RANN)
#' @param verbose Displays progress bar
#'
#' @return Returns a vector of values representing the entropy metric from each
#' bootstrapped iteration.
#'
#' @importFrom RANN nn2
#' @importFrom pbapply pbsapply
#' @importFrom future.apply future_sapply
#' @export
#'
MixingMetric <- function(
  object,
  grouping.var,
  reduction = "pca",
  dims = 1:2,
  k = 5,
  max.k = 300,
  eps = 0,
  verbose = TRUE
) {
  my.sapply <- ifelse(
    test = verbose && PlanThreads() == 1,
    yes = pbsapply,
    no = future_sapply
  )
  embeddings <- Embeddings(object = object[[reduction]])[, dims]
  nn <- nn2(
    data = embeddings,
    k = max.k,
    eps = eps
  )
  group.info <- object[[grouping.var, drop = TRUE]]
  groups <- unique(x = group.info)
  mixing <- my.sapply(
    X = 1:ncol(x = object),
    FUN = function(x) {
      sapply(X = groups, FUN = function(y) {
        which(x = group.info[nn$nn.idx[x, ]] == y)[k]
      })
    }
  )
  mixing[is.na(x = mixing)] <- max.k
  mixing <- apply(
    X = mixing,
    MARGIN = 2,
    FUN = median
  )
  return(mixing)
}

#' Select integration features
#'
#' Choose the features to use when integrating multiple datasets. This function ranks features by
#' the number of datasets they appear in, breaking ties by the median rank across datasets. It
#' returns the highest features by this ranking.
#'
#' @param object.list List of seurat objects
#' @param nfeatures Number of features to return
#' @param assay Name of assay from which to pull the variable features.
#' @param verbose Print messages
#' @param fvf.nfeatures nfeatures for FindVariableFeatures. Used if VariableFeatures have not been
#' set for any object in object.list.
#' @param ... Additional parameters to \code{\link{FindVariableFeatures}}
#'
#' @export
#'
SelectIntegrationFeatures <- function(
  object.list,
  nfeatures = 2000,
  assay = NULL,
  verbose = TRUE,
  fvf.nfeatures = 2000,
  ...
) {
  if (!is.null(x = assay)) {
    if (length(x = assay) != length(x = object.list)) {
      stop("If specifying the assay, please specify one assay per object in the object.list")
    }
    for(ii in length(x = object.list)) {
      DefaultAssay(object = object.list[[ii]]) <- assay[ii]
    }
  } else {
    assay <- sapply(X = object.list, FUN = DefaultAssay)
  }
  for(ii in 1:length(x = object.list)) {
    if (length(x = VariableFeatures(object = object.list[[ii]])) == 0) {
      if (verbose) {
        message(paste0("No variable features found for object", ii, " in the object.list. Running FindVariableFeatures ..."))
      }
      object.list[[ii]] <- FindVariableFeatures(object = object.list[[ii]], nfeatures = fvf.nfeatures, verbose = verbose, ...)
    }
  }
  var.features <- unname(obj = unlist(x = lapply(
    X = 1:length(x = object.list),
    FUN = function(x) VariableFeatures(object = object.list[[x]], assay = assay[x]))
  ))
  var.features <- sort(x = table(var.features), decreasing = TRUE)
  for(i in 1:length(x = object.list)) {
    var.features <- var.features[names(x = var.features) %in% rownames(x = object.list[[i]][[assay[i]]])]
  }
  tie.val <- var.features[min(nfeatures, length(x = var.features))]
  features <- names(x = var.features[which(x = var.features > tie.val)])
  if (length(x = features) > 0) {
    feature.ranks <- sapply(X = features, FUN = function(x) {
      ranks <- sapply(X = object.list, FUN = function(y) {
        vf <- VariableFeatures(object = y)
        if (x %in% vf){
          return(which(x = x == vf))
        }
        return(NULL)
      })
      median(x = unlist(x = ranks))
    })
    features <- names(x = sort(x = feature.ranks))
  }
  features.tie <- var.features[which(x = var.features == tie.val)]
  tie.ranks <- sapply(X = names(x = features.tie), FUN = function(x) {
    ranks <- sapply(X = object.list, FUN = function(y) {
      vf <- VariableFeatures(object = y)
      if (x %in% vf){
        return(which(x = x == vf))
      }
      return(NULL)
    })
    median(x = unlist(x = ranks))
  })
  features <- c(features, names(x = head(x = sort(x = tie.ranks), nfeatures - length(x = features))))
}

#' Transfer Labels
#'
#' Transfers the labels
#'
#' @param anchorset Results from FindTransferAnchors
#' @param refdata Data to transfer. Should be either a vector where the names correspond to
#' reference cells, or a matrix, where the column names correspond to the reference cells.
#' @param reduction Dimensional reduction to use for the weighting. Options are:
#' \itemize{
#'    \item{pcaproject: Use the projected PCA used for anchor building}
#'    \item{pca: Use an internal PCA on the query only}
#'    \item{cca: Use the CCA used for anchor building}
#' }
#' @param l2.norm Perform L2 normalization on the cell embeddings after dimensional reduction
#' @param dims Number of PCs to use in the weighting procedure
#' @param k.weight Number of neighbors to consider when weighting
#' @param sd.weight Controls the bandwidth of the Gaussian kernel for weighting
#' @param eps Error bound on the neighbor finding algorithm (from RANN)
#' @param do.cpp Run cpp code where applicable
#' @param verbose Print progress bars and output
#' @param slot Slot to store the imputed data
#'
#' @return If refdata is a vector, returns a dataframe with label predictions. If refdata is a
#' matrix, returns an Assay object where the imputed data has been stored in the provided slot.
#'
#' @export
#'
TransferData <- function(
  anchorset,
  refdata,
  reduction = 'pcaproject',
  l2.norm = FALSE,
  dims = 1:30,
  k.weight = 50,
  sd.weight = 1,
  eps = 0,
  do.cpp = TRUE,
  verbose = TRUE,
  slot = "data"
) {

  combined.ob <- slot(object = anchorset, name = "object.list")[[1]]
  anchors <- slot(object = anchorset, name = "anchors")
  reference.cells <- slot(object = anchorset, name = "reference.cells")
  query.cells <- slot(object = anchorset, name = "query.cells")

  if (inherits(x = refdata, what = c("character", "factor"))) {
    if (length(x = refdata) != length(x = reference.cells)) {
      stop(paste0("Please provide a vector that is the same length as the number of reference cells",
                  " used in anchor finding.\n",
                  "Length of vector provided: ", length(x = refdata), "\n",
                  "Length of vector required: ", length(x = reference.cells)))
    }
    label.transfer <- TRUE
  } else if (inherits(x = refdata, what = c("dgCMatrix", "matrix"))) {
    if (ncol(x = refdata) != length(x = reference.cells)) {
      stop(paste0("Please provide a matrix that has the same number of columns as the number of reference cells",
                  " used in anchor finding.\n",
                  "Number of columns in provided matrix : ", ncol(x = refdata), "\n",
                  "Number of columns required           : ", length(x = reference.cells)))
    }
    if (any(!colnames(x = refdata) == reference.cells)) {
      if (any(!colnames(x = refdata) %in% reference.cells) | any(!reference.cells %in% colnames(x = refdata))) {
        stop("Some (or all) of the column names of the provided refdata don't match the reference cells used in anchor finding.")
      }
      refdata <- refdata[, reference.cells]
    }
    if (!slot %in% c("counts", "data")) {
      stop("Please specify slot as either 'counts' or 'data'.")
    }
    label.transfer <- FALSE
  } else {
    stop(paste0("Please provide either a vector (character or factor) for label transfer or a matrix",
         "for feature transfer.\n", "Type provided: ", class(x = refdata)))
  }

  if (reduction == 'pca') {
    message("Running PCA on query dataset")
    features <- slot(object = anchorset, name = "anchor.features")
    query <- combined.ob[features, query.cells]
    query <- ScaleData(object = query, features = features, verbose = FALSE)
    query <- RunPCA(object = query, npcs = max(dims), features = features, verbose = FALSE)
    query.pca <- Embeddings(query[['pca']])

    #fill with 0s
    ref.pca <- matrix(
      data = 0,
      nrow = length(x = reference.cells),
      ncol = ncol(x = query.pca),
      dimnames = list(reference.cells, colnames(x = query.pca))
    )
    combined.pca.embeddings <- rbind(ref.pca, query.pca)[colnames(x = combined.ob), ]
    combined.pca <- CreateDimReducObject(
      embeddings = combined.pca.embeddings,
      key = "PC_",
      assay = DefaultAssay(object = combined.ob)
    )
    combined.ob[["pca"]] <- combined.pca
    if (l2.norm) {
      combined.ob <- L2Dim(object = combined.ob, reduction = 'pca')
    }
  }
  if (l2.norm) {
    reduction <- paste0(reduction, ".l2")
  }

  combined.ob <- SetIntegrationData(
    object = combined.ob,
    integration.name = "integrated",
    slot = 'anchors',
    new.data = anchors
  )
  combined.ob <- SetIntegrationData(
    object = combined.ob,
    integration.name = "integrated",
    slot = 'neighbors',
    new.data = list('cells1' = reference.cells, 'cells2' = query.cells)
  )
  combined.ob <- FindIntegrationMatrix(
    object = combined.ob,
    verbose = verbose
  )
  combined.ob <- FindWeights(
    object = combined.ob,
    reduction = reduction,
    nn.reduction = reduction,
    nn.dims = dims,
    k.weight = k.weight,
    sd.weight = sd.weight,
    eps = eps,
    do.cpp = do.cpp,
    verbose = verbose
  )
  weights <- GetIntegrationData(
    object = combined.ob,
    integration.name = "integrated",
    slot = 'weights'
  )
  anchors <- as.data.frame(x = anchors)
  query.cells <- unname(obj = sapply(
    X = query.cells,
    FUN = function(x) gsub(pattern = "_query", replacement = "", x = x)
  ))
  # case for projection
  if (label.transfer) {
    anchors$id1 <- refdata[anchors[, "cell1"]]
    reference.ids <- factor(x = anchors$id1, levels = unique(x = refdata))
    possible.ids <- levels(x = reference.ids)
    prediction.mat <- matrix(nrow = nrow(x = anchors), ncol = length(x = possible.ids), data = 0)
    for(i in 1:length(x = possible.ids)) {
      prediction.mat[which(reference.ids == possible.ids[i]), i] = 1
    }
    if (verbose) {
      message("Predicting cell labels")
    }
    prediction.scores <- t(x = weights) %*% prediction.mat
    colnames(x = prediction.scores) <- possible.ids
    rownames(x = prediction.scores) <- query.cells
    prediction.ids <- possible.ids[apply(X = prediction.scores, MARGIN = 1, FUN = which.max)]
    prediction.ids <- as.character(prediction.ids)
    prediction.scores <- cbind(prediction.scores, max = apply(X = prediction.scores, MARGIN = 1, FUN = max))
    predictions <- (data.frame(
      predicted.id = prediction.ids,
      prediction.score = as.matrix(prediction.scores),
      row.names = query.cells,
      stringsAsFactors = FALSE)
    )
    return(predictions)
  } else {  # case for transferring features
    reference.cell.indices <- reference.cells[anchors$cell1]
    refdata.anchors <- refdata[, reference.cell.indices]
    nfeatures <- nrow(x = refdata)
    if (verbose) {
      message(paste0("Transfering ", nfeatures, " features onto reference data"))
    }
    new.data <- refdata.anchors %*% weights
    rownames(new.data) <- rownames(refdata)
    colnames(new.data) <- query.cells
    if (slot == "counts") {
      new.assay <- CreateAssayObject(counts = new.data)
    } else if (slot == "data") {
      new.assay <- CreateAssayObject(data = new.data)
    }
    return(new.assay)
  }
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Methods for Seurat-defined generics
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Internal
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Add dataset number and remove cell offset
#
# Record which dataset number in the original list of Seurat objects
# each anchor cell came from, and correct the cell index so it corresponds to
# the position of the anchor cell in its own dataset
#
# @param anchor.df Dataframe of anchors
# @param offsets size of each dataset in anchor dataframe
# @param obj.length Vector of object lengths
#
# @return Anchor dataframe with additional columns corresponding to the dataset of each cell
#
AddDatasetID <- function(
  anchor.df,
  offsets,
  obj.lengths
) {
  ndataset <- length(x = offsets)
  total.cells <- sum(obj.lengths)
  offsets <- c(offsets, total.cells)
  row.offset <- rep.int(x = offsets[1:ndataset], times = obj.lengths)
  dataset <- rep.int(x = 1:ndataset, times = obj.lengths)
  anchor.df <- data.frame(
    'cell1' = anchor.df[, 1] - row.offset[anchor.df[, 1]],
    'cell2' = anchor.df[, 2] - row.offset[anchor.df[, 2]],
    'score' = anchor.df[, 3],
    'dataset1' = dataset[anchor.df[, 1]],
    'dataset2' = dataset[anchor.df[, 2]]
  )
  return(anchor.df)
}

# Add info to anchor matrix
#
# @param object Seurat object
# @param toolname Name in tool slot to pull from
# @param annotation Name in metadata to annotate anchors with
# @param object.list List of objects using in FindIntegrationAnchors call
#
# @return Returns the anchor dataframe with additional columns for annotation metadata
#
AnnotateAnchors <- function(
  object,
  toolname = "integrated",
  annotation = NULL,
  object.list = NULL
) {
  anchors <- GetIntegrationData(
    object = object,
    integration.name = toolname,
    slot = 'anchors'
  )
  for(i in annotation) {
    if (! i %in% colnames(x = object[[]])) {
      warning(i, " not in object metadata")
      next
    }
    if(!is.null(x = object.list)) {
      anchors[, paste0("cell1.", i)] <- apply(X = anchors, MARGIN = 1, function(x){
        as.character(object.list[[as.numeric(x[["dataset1"]])]][[]][as.numeric(x[["cell1"]]), i])
      })
      anchors[, paste0("cell2.", i)] <- apply(X = anchors, MARGIN = 1, function(x){
        as.character(object.list[[as.numeric(x[["dataset2"]])]][[]][as.numeric(x[["cell2"]]), i])
      })
    } else {
      cells1 <- GetIntegrationData(
        object = object,
        integration.name = toolname,
        slot = 'neighbors'
      )$cells1
      cells2 <- GetIntegrationData(
        object = object,
        integration.name = toolname,
        slot = 'neighbors'
      )$cells2
      anchors[, paste0("cell1.", i)] <- object[[i]][cells1[anchors$cell1], , drop = TRUE]
      anchors[, paste0("cell2.", i)] <- object[[i]][cells2[anchors$cell2], , drop = TRUE]
      anchors[, paste0(i, ".match")] <- anchors[, paste0("cell1.", i)] == anchors[, paste0("cell2.", i)]
    }
  }
  return(anchors)
}

# Build tree of datasets based on cell similarity
#
# @param similarity.matrix Dataset similarity matrix
#
# @return Returns a heirarchical clustering of datasets
#
#' @importFrom stats hclust
#
BuildSampleTree <- function(similarity.matrix) {
  dist.mat <- as.dist(m = 1 / similarity.matrix)
  clusters <- hclust(d = dist.mat)
  return(clusters$merge)
}

# Construct nearest neighbor matrix from nn.idx
#
# @param nn.idx Nearest neighbor index matrix (nn.idx from RANN)
# @param offset1 Offsets for the first neighbor
# @param offset2 Offsets for the second neighbor
#
# @return returns a sparse matrix representing the NN matrix
#
ConstructNNMat <- function(nn.idx, offset1, offset2, dims) {
  k <- ncol(x = nn.idx)
  j <- as.numeric(x = t(x = nn.idx)) + offset2
  i <- ((1:length(x = j)) - 1) %/% k + 1 + offset1
  nn.mat <- sparseMatrix(i = i, j = j, x = 1, dims = dims)
  return(nn.mat)
}

# Count anchors between all datasets
#
# Counts anchors between each dataset and scales based on total number of cells in the datasets
#
# @param anchor.df Matrix of anchors
# @param offsets Dataset sizes in anchor matrix. Used to identify boundaries of each dataset in
# matrix, so that total pairwise anchors between all datasets can be counted
#
# @return Returns a similarity matrix
#
CountAnchors <- function(
  anchor.df,
  offsets,
  obj.lengths
) {
  similarity.matrix <- matrix(data = 0, ncol = length(x = offsets), nrow = length(x = offsets))
  similarity.matrix[upper.tri(x = similarity.matrix, diag = TRUE)] <- NA
  total.cells <- sum(obj.lengths)
  offsets <- c(offsets, total.cells)
  for (i in 1:nrow(x = similarity.matrix)){
    for (j in 1:ncol(x = similarity.matrix)){
      if (!is.na(x = similarity.matrix[i, j])){
        relevant.rows <- anchor.df[(anchor.df$dataset1 %in% c(i, j)) & (anchor.df$dataset2 %in% c(i, j)), ]
        score <- nrow(x = relevant.rows)
        ncell <- min(obj.lengths[[i]], obj.lengths[[j]])
        similarity.matrix[i, j] <- score / ncell
      }
    }
  }
  return(similarity.matrix)
}

FilterAnchors <- function(
  object,
  assay = NULL,
  integration.name = 'integrated',
  features = NULL,
  k.filter = 200,
  eps = 0,
  verbose = TRUE
) {
  if (verbose) {
    message("Filtering Anchors")
  }
  assay <- assay %||% DefaultAssay(object = object)
  features <- features %||% VariableFeatures(object = object)
  if (length(x = features) == 0) {
    stop("No features provided and no VariableFeatures computed.")
  }
  features <- unique(x = features)
  neighbors <- GetIntegrationData(object = object, integration.name = integration.name, slot = 'neighbors')
  nn.cells1 <- neighbors$cells1
  nn.cells2 <- neighbors$cells2
  cn.data1 <- L2Norm(
    mat = as.matrix(x = t(x = GetAssayData(
      object = object[[assay[1]]],
      slot = "data")[features, nn.cells1])),
    MARGIN = 1)
  cn.data2 <- L2Norm(
    mat = as.matrix(x = t(x = GetAssayData(
      object = object[[assay[2]]],
      slot = "data")[features, nn.cells2])),
    MARGIN = 1)
  nn <- nn2(
    data = cn.data2[nn.cells2, ],
    query = cn.data1[nn.cells1, ],
    k = k.filter,
    eps = eps
  )

  anchors <- GetIntegrationData(object = object, integration.name = integration.name, slot = "anchors")
  position <- sapply(X = 1:nrow(x = anchors), FUN = function(x) {
    which(x = anchors[x, "cell2"] == nn$nn.idx[anchors[x, "cell1"], ])[1]
  })
  anchors <- anchors[!is.na(x = position), ]
  if (verbose) {
    message("\tRetained ", nrow(x = anchors), " anchors")
  }
  object <- SetIntegrationData(
    object = object,
    integration.name = integration.name,
    slot = "anchors",
    new.data = anchors
  )
  return(object)
}

FindAnchors <- function(
  object.pair,
  assay,
  cells1,
  cells2,
  reduction,
  dims = 1:10,
  k.anchor = 5,
  k.filter = 200,
  k.score = 30,
  max.features = 200,
  eps = 0,
  verbose = TRUE
) {
  # compute local neighborhoods, use max of k.anchor and k.score if also scoring to avoid
  # recomputing neighborhoods
  k.neighbor <- k.anchor
  if (!is.na(x = k.score)) {
    k.neighbor <- max(k.anchor, k.score)
  }
  object.pair <- FindNN(
    object = object.pair,
    cells1 = cells1,
    cells2 = cells2,
    dims = dims,
    reduction = reduction,
    k = k.neighbor,
    eps = eps,
    verbose = verbose
  )
  object.pair <- FindAnchorPairs(
    object = object.pair,
    integration.name = "integrated",
    k.anchor = k.anchor,
    verbose = verbose
  )
  if (!is.na(x = k.filter)) {
    top.features <- TopDimFeatures(
      object = object.pair,
      reduction = reduction,
      dims = dims,
      features.per.dim = 100,
      max.features = max.features
    )
    object.pair <- FilterAnchors(
      object = object.pair,
      assay = assay,
      integration.name = 'integrated',
      features = top.features,
      k.filter = k.filter,
      eps = eps,
      verbose = verbose
    )
  }
  if (!is.na(x = k.score)) {
    object.pair = ScoreAnchors(
      object = object.pair,
      assay = DefaultAssay(object = object.pair),
      integration.name = "integrated",
      k.score = k.score
    )
  }
  anchors <- GetIntegrationData(
    object = object.pair,
    integration.name = 'integrated',
    slot = 'anchors'
  )
  return(anchors)
}

# Find Anchor pairs
#
FindAnchorPairs <- function(
  object,
  integration.name = 'integrated',
  cells1 = NULL,
  cells2 = NULL,
  k.anchor = 5,
  verbose = TRUE
) {
  neighbors <- GetIntegrationData(object = object, integration.name = integration.name, slot = 'neighbors')
  max.nn <- c(ncol(x = neighbors$nnab$nn.idx), ncol(x = neighbors$nnba$nn.idx))
  if (any(k.anchor > max.nn)) {
    message(paste0('warning: requested k.anchor = ', k.anchor, ', only ', min(max.nn), ' in dataset'))
    k.anchor <- min(max.nn)
  }
  if (verbose) {
    message("Finding mutual nearest neighborhoods")
  }
  if (is.null(x = cells1)) {
    cells1 <- colnames(x = object)
  }
  if (is.null(x = cells2)) {
    cells2 <- colnames(x = object)
  }
  if (!(cells1 %in% colnames(object)) || !(cells2 %in% colnames(object))) {
    warning("Requested cells not contained in Seurat object. Subsetting list of cells.")
    cells1 <- intersect(x = cells1, y = colnames(x = object))
    cells2 <- intersect(x = cells2, y = colnames(x = object))
  }
  # convert cell name to neighbor index
  nn.cells1 <- neighbors$cells1
  nn.cells2 <- neighbors$cells2
  cell1.index <- sapply(X = cells1, FUN = function(x) return(which(x == nn.cells1)))
  cell2.index <- sapply(X = cells2, FUN = function(x) return(which(x == nn.cells2)))

  ncell <- 1:nrow(x = neighbors$nnab$nn.idx)
  ncell <- ncell[ncell %in% cell1.index]
  anchors <- list()
  # pre allocate vector
  anchors$cell1 <- rep(x = 0, length(x = ncell) * 5)
  anchors$cell2 <- anchors$cell1
  anchors$score <- anchors$cell1 + 1
  idx <- 0
  for (cell in ncell) {
    neighbors.ab <- neighbors$nnab$nn.idx[cell, 1:k.anchor]
    mutual.neighbors <- which(
      x = neighbors$nnba$nn.idx[neighbors.ab, 1:k.anchor, drop = FALSE] == cell,
      arr.ind = TRUE
    )[, 1]
    for (i in neighbors.ab[mutual.neighbors]){
      idx <- idx + 1
      anchors$cell1[idx] <- cell
      anchors$cell2[idx] <- i
      anchors$score[idx] <- 1
    }
  }
  anchors$cell1 <- anchors$cell1[1:idx]
  anchors$cell2 <- anchors$cell2[1:idx]
  anchors$score <- anchors$score[1:idx]
  anchors <- t(x = do.call(what = rbind, args = anchors))
  anchors <- as.matrix(x = anchors)
  object <- SetIntegrationData(
    object = object,
    integration.name = integration.name,
    slot = 'anchors',
    new.data = anchors
  )
  if (verbose) {
    message(paste0("\tFound ", nrow(x = anchors), " anchors"))
  }
  object <- LogSeuratCommand(object = object)
  command.name <- LogSeuratCommand(object = object, return.command = TRUE)
  return(object)
}

FindIntegrationMatrix <- function(
  object,
  assay = NULL,
  integration.name = 'integrated',
  features.integrate = NULL,
  verbose = TRUE
) {
  assay <- assay %||% DefaultAssay(object = object)
  neighbors <- GetIntegrationData(object = object, integration.name = integration.name, slot = 'neighbors')
  nn.cells1 <- neighbors$cells1
  nn.cells2 <- neighbors$cells2
  anchors <- GetIntegrationData(
    object = object,
    integration.name = integration.name,
    slot = 'anchors'
  )
  if (verbose) {
    message("Finding integration vectors")
  }
  features.integrate <- features.integrate %||% rownames(
    x = GetAssayData(object = object, assay = assay, slot = "data")
  )
  data.use1 <- t(x = GetAssayData(
    object = object,
    assay = assay,
    slot = "data")[features.integrate, nn.cells1]
  )
  data.use2 <- t(x = GetAssayData(
    object = object,
    assay = assay,
    slot = "data")[features.integrate, nn.cells2]
  )
  anchors1 <- nn.cells1[anchors[, "cell1"]]
  anchors2 <- nn.cells2[anchors[, "cell2"]]
  data.use1 <- data.use1[anchors1, ]
  data.use2 <- data.use2[anchors2, ]
  integration.matrix <- data.use2 - data.use1
  object <- SetIntegrationData(
    object = object,
    integration.name = integration.name,
    slot = 'integration.matrix',
    new.data = integration.matrix
  )
  object <- LogSeuratCommand(object = object)
  command.name <- LogSeuratCommand(object = object, return.command = TRUE)
  return(object)
}

# Find nearest neighbors
#
FindNN <- function(
  object,
  cells1 = NULL,
  cells2 = NULL,
  grouping.var = NULL,
  dims = 1:10,
  reduction = "cca.l2",
  nn.dims = dims,
  nn.reduction = reduction,
  k = 300,
  eps = 0,
  integration.name = 'integrated',
  verbose = TRUE
) {
  if (xor(x = is.null(x = cells1), y = is.null(x = cells2))) {
    stop("cells1 and cells2 must both be specified")
  }
  if (!is.null(x = cells1) && !is.null(x = cells2) && !is.null(x = grouping.var)) {
    stop("Specify EITHER grouping.var or cells1/2.")
  }
  if (is.null(x = cells1) && is.null(x = cells2) && is.null(x = grouping.var)) {
    stop("Please set either cells1/2 or grouping.var")
  }
  if (!is.null(x = grouping.var)) {
    if (nrow(x = unique(x = object[[grouping.var]])) != 2) {
      stop("Number of groups in grouping.var not equal to 2.")
    }
    groups <- names(x = sort(x = table(... = object[[grouping.var]]), decreasing = TRUE))
    cells1 <- colnames(x = object)[object[[grouping.var]] == groups[[1]]]
    cells2 <- colnames(x = object)[object[[grouping.var]] == groups[[2]]]
  }
  if (verbose) {
    message("Finding neighborhoods")
  }
  dim.data.self <- Embeddings(object = object[[nn.reduction]])[ ,nn.dims]
  dim.data.opposite <- Embeddings(object = object[[reduction]])[ ,dims]
  dims.cells1.self <- dim.data.self[cells1, ]
  dims.cells1.opposite <- dim.data.opposite[cells1, ]
  dims.cells2.self <- dim.data.self[cells2, ]
  dims.cells2.opposite <- dim.data.opposite[cells2, ]

  nnaa <- nn2(
    data = dims.cells1.self,
    k = k + 1,
    eps = eps
  )
  nnab <- nn2(
    data = dims.cells2.opposite,
    query = dims.cells1.opposite,
    k = k,
    eps = eps
  )
  nnbb <- nn2(
    data = dims.cells2.self,
    k = k + 1,
    eps = eps
  )
  nnba <- nn2(
    data = dims.cells1.opposite,
    query = dims.cells2.opposite,
    k = k,
    eps = eps
  )
  object <- SetIntegrationData(
    object = object,
    integration.name = integration.name,
    slot = 'neighbors',
    new.data = list('nnaa' = nnaa, 'nnab' = nnab, 'nnba' = nnba, 'nnbb' = nnbb, 'cells1' = cells1, 'cells2' = cells2)
  )
  object <- LogSeuratCommand(object = object)
  command.name <- LogSeuratCommand(object = object, return.command = TRUE)
  return(object)
}

FindWeights <- function(
  object,
  assay = NULL,
  integration.name = 'integrated',
  reduction = "cca",
  nn.reduction = "cca",
  nn.dims = 1:10,
  features = NULL,
  k.weight = 300,
  sd.weight = 1,
  eps = 0,
  min.dist = 0,
  verbose = TRUE,
  do.cpp = FALSE
) {
  if (verbose) {
    message("Finding integration vector weights")
  }
  assay <- assay %||% DefaultAssay(object = object)
  neighbors <- GetIntegrationData(object = object, integration.name = integration.name, slot = 'neighbors')
  nn.cells1 <- neighbors$cells1
  nn.cells2 <- neighbors$cells2
  anchors <- GetIntegrationData(
    object = object,
    integration.name = integration.name,
    slot = 'anchors'
  )
  anchors.cells2 <- nn.cells2[anchors[, "cell2"]]
  if (is.null(x = features)) {
    data.use <- Embeddings(object[[nn.reduction]])[nn.cells2, nn.dims]
  } else {
    data.use <- t(x = GetAssayData(object = object, slot = 'data', assay = assay)[features, nn.cells2])
  }
  knn_2_2 <- nn2(
    data = data.use[anchors.cells2, ],
    query = data.use,
    k = k.weight + 1,
    eps = eps
  )
  distances <- knn_2_2$nn.dists[, -1]
  distances <- 1 - (distances / distances[, ncol(x = distances)])
  cell.index <- knn_2_2$nn.idx[, -1]
  integration.matrix <- GetIntegrationData(
    object = object,
    integration.name = integration.name,
    slot = "integration.matrix"
  )
  if (do.cpp) {
    weights <- FindWeightsC(
      integration_matrix = as(integration.matrix, "dgCMatrix"),
      cells2 = 0:(length(x = nn.cells2) - 1),
      distances = as.matrix(x = distances),
      anchor_cells2 = anchors.cells2,
      integration_matrix_rownames = rownames(x = integration.matrix),
      cell_index = cell.index,
      anchor_score = anchors[, "score"],
      min_dist = min.dist,
      sd = sd.weight,
      display_progress = verbose
    )
  } else {
    if (verbose) {
      pb <- txtProgressBar(min = 1, max = length(x = nn.cells2), initial = 1, style = 3, file = stderr())
    }
    dist.weights <- matrix(
      data = min.dist,
      nrow = nrow(x = integration.matrix),
      ncol = length(x = nn.cells2)
    )
    for (cell in 1:length(x = nn.cells2)) {
      wt <- distances[cell, ]
      cellnames <- anchors.cells2[cell.index[cell, ]]
      names(x = wt) <- cellnames
      for (i in cellnames){
        anchor.index <- which(rownames(integration.matrix) == i)
        dist.weights[anchor.index, cell] <- wt[[i]]
      }
      if(verbose) setTxtProgressBar(pb, cell)
    }
    if(verbose) message("")
    dist.anchor.weight <- dist.weights * anchors[, "score"]
    weights <- 1 - exp(-1 * dist.anchor.weight / (2 * (1 / sd.weight)) ^ 2)
    weights <- sweep(weights, 2, Matrix::colSums(weights), "/")
  }
  object <- SetIntegrationData(
    object = object,
    integration.name = integration.name,
    slot = 'weights',
    new.data = weights
  )
  object <- LogSeuratCommand(object = object)
  command.name <- LogSeuratCommand(object = object, return.command = TRUE)
  return(object)
}

# Convert nearest neighbor information to a sparse matrix
#
# @param idx Nearest neighbor index
# @param distance Nearest neighbor distance
# @param k Number of nearest neighbors
#
NNtoMatrix <- function(idx, distance, k) {
  nn <- list()
  x <- 1
  for (i in 1:nrow(x = idx)) {
    for (j in 2:k) {
      nn.idx <- idx[i, j]
      nn.dist <- distance[i, j]
      nn[[x]] <- c('i' = i, 'j' = nn.idx, 'x' = 1/nn.dist)
      x <- x + 1
    }
  }
  nn <- do.call(what = rbind, args = nn)
  nn.matrix <- new(
    Class = 'dgTMatrix',
    i = as.integer(x = nn[, 1] - 1),
    j = as.integer(x = nn[, 2] - 1),
    x = as.numeric(x = nn[, 3]),
    Dim = as.integer(x = c(nrow(idx), nrow(x = idx)))
  )
  nn.matrix <- as(object = nn.matrix, Class = 'dgCMatrix')
  return(nn.matrix)
}

# Parse merge information from dataset clustering
#
# @param clustering clustering dataframe from hclust ($merge).
#  Gives the order of merging datasets to get to the root of the tree.
# @param i current row in clustering dataframe
#
ParseMergePair <- function(clustering, i){
  # return 2-element list of datasets in first and second object
  datasets <- list('object1' = clustering[i, 1], 'object2' = clustering[i, 2])
  if (datasets$object1 > 0) {
    datasets$object1 <- ParseRow(clustering, datasets$object1)
  }
  if (datasets$object2 > 0) {
    datasets$object2 <- ParseRow(clustering, datasets$object2)
  }
  datasets$object1 <- abs(x = datasets$object1)
  datasets$object2 <- abs(x = datasets$object2)
  return(datasets)
}

# Parse row of clustering order
#
# Used recursively to work out the dataset composition of a merged object
#
# @param clustering clustering dataframe from hclust ($merge).
#  Gives the order of merging datasets to get to the root of the tree.
# @param i current row in clustering dataframe
#
ParseRow <- function(clustering, i){
  # returns vector of datasets
  datasets <- as.list(x = clustering[i, ])
  if (datasets[[1]] > 0) {
    datasets[[1]] <- ParseRow(clustering = clustering, i = datasets[[1]])
  }
  if (datasets[[2]] > 0) {
    datasets[[2]] <- ParseRow(clustering = clustering, i = datasets[[2]])
  }
  return(unlist(datasets))
}

ProjectCellEmbeddings <- function(
  reference,
  query,
  reference.assay = NULL,
  query.assay = NULL,
  dims = 1:50,
  verbose = TRUE,
  feature.mean = NULL,
  feature.sd = NULL
) {
  if (verbose) {
    message("Projecting PCA")
  }
  reduction <- "pca"
  reference.assay <- reference.assay %||% DefaultAssay(object = reference)
  query.assay <- query.assay %||% DefaultAssay(object = query)
  features <- rownames(x = Loadings(object = reference[[reduction]]))
  features <- intersect(x = features, y = rownames(x = query[[query.assay]]))

  reference.data <-  GetAssayData(
    object = reference,
    assay.use = reference.assay,
    slot = "data")[features, ]
  query.data <- GetAssayData(
    object = query,
    assay.use = query.assay,
    slot = "data")[features, ]

  if (is.null(x = feature.mean)){
    feature.mean <- rowMeans(x = reference.data)
    feature.sd <- sqrt(SparseRowVar2(mat = reference.data, mu = feature.mean, display_progress = FALSE))
    feature.sd[is.na(x = feature.sd)] <- 1
    feature.mean[is.na(x = feature.mean)] <- 1
  }
  query.feature.mean <- rowMeans(x = query.data)
  proj.data <- GetAssayData(
    object = query,
    assay = query.assay,
    slot = "data"
  )[features, ]
  store.names <- dimnames(x = proj.data)
  proj.data <- FastSparseRowScaleWithKnownStats(
    mat = proj.data,
    mu = feature.mean,
    sigma = feature.sd,
    display_progress = FALSE
  )
  dimnames(x = proj.data) <- store.names
  ref.feature.loadings <- Loadings(object = reference[[reduction]])[features, dims]
  proj.pca <- t(crossprod(x = ref.feature.loadings, y = proj.data))
  return(proj.pca)
}

# Calculate position along a defined reference range for a given vector of
# numerics. Will range from 0 to 1.
#
# @param x      Vector of numeric type
# @param lower  Lower end of reference range
# @param upper  Upper end of reference range
#
#' @importFrom stats quantile
#
# @return       Returns a vector that describes the position of each element in
#               x along the defined reference range
#
ReferenceRange <- function(x, lower = 0.025, upper = 0.975) {
  return((x - quantile(x = x, probs = lower)) /
           (quantile(x = x, probs = upper) - quantile(x = x, probs = lower)))
}

ScoreAnchors <- function(
  object,
  assay = NULL,
  integration.name = 'integrated',
  verbose = TRUE,
  k.score = 30,
  do.cpp = TRUE
) {
  assay <- assay %||% DefaultAssay(object = object)
  anchor.df <- as.data.frame(x = GetIntegrationData(object = object, integration.name = integration.name, slot = 'anchors'))
  neighbors <- GetIntegrationData(object = object, integration.name = integration.name, slot = "neighbors")
  offset <- length(x = neighbors$cells1)
  anchor.df$cell2 <- anchor.df$cell2 + offset
  # make within dataset df
  if (verbose) {
    message("Extracting within-dataset neighbors!")
  }
  total.cells <- offset + length(neighbors$cells2)
  nn.m1 <- ConstructNNMat(nn.idx = neighbors$nnaa$nn.idx[,1:k.score], offset1 = 0, offset2 = 0, dims = c(total.cells, total.cells))
  nn.m2 <- ConstructNNMat(nn.idx = neighbors$nnab$nn.idx[,1:k.score], offset1 = 0, offset2 = offset, dims = c(total.cells, total.cells))
  nn.m3 <- ConstructNNMat(nn.idx = neighbors$nnba$nn.idx[,1:k.score], offset1 = offset, offset2 = 0, dims = c(total.cells, total.cells))
  nn.m4 <- ConstructNNMat(nn.idx = neighbors$nnbb$nn.idx[,1:k.score], offset1 = offset, offset2 = offset, dims = c(total.cells, total.cells))
  k.matrix <- nn.m1 + nn.m2 + nn.m3 + nn.m4
  anchor.only <- sparseMatrix(i = anchor.df[, 1], j = anchor.df[, 2], x = 1, dims = c(total.cells, total.cells))

  if (do.cpp){
    anchor.matrix <- SNNAnchor(k_matrix = k.matrix, anchor_only = anchor.only)
  } else {
    jaccard.dist <- tcrossprod(x = k.matrix)
    anchor.matrix <- jaccard.dist * anchor.only
  }

  anchor.matrix <- as(object = anchor.matrix, Class = "dgTMatrix")
  anchor.new <- data.frame(
    'cell1' = anchor.matrix@i + 1,
    'cell2' = anchor.matrix@j + 1,
    'score' = anchor.matrix@x
  )
  anchor.new$cell2 <- anchor.new$cell2 - offset
  max.score <- quantile(anchor.new$score, 0.9)
  min.score <- quantile(anchor.new$score, 0.01)
  anchor.new$score <- anchor.new$score - min.score
  anchor.new$score <- anchor.new$score / (max.score - min.score)
  anchor.new$score[anchor.new$score > 1] <-  1
  anchor.new$score[anchor.new$score < 0] <- 0
  anchor.new <- as.matrix(x = anchor.new)
  object <- SetIntegrationData(
    object = object,
    integration.name = integration.name,
    slot = 'anchors',
    new.data = anchor.new
  )
  object <- LogSeuratCommand(object = object)
  command.name <- LogSeuratCommand(object = object, return.command = TRUE)
  return(object)
}

# Get top n features across given set of dimensions
#
# @param object Seurat object
# @param reduction Which dimension reduction to use
# @param dims Which dimensions to use
# @param features.per.dim How many features to consider per dimension
# @param max.features Number of features to return at most
#
TopDimFeatures <- function(
  object,
  reduction,
  dims = 1:10,
  features.per.dim = 100,
  max.features = 200
) {
  dim.reduction <- object[[reduction]]
  max.features <- max(length(x = dims) * 2, max.features)
  num.features <- sapply(X = 1:features.per.dim, FUN = function(y) {
    length(x = unique(x = as.vector(x = sapply(X = dims, FUN = function(x) {
        unlist(x = TopFeatures(object = dim.reduction, dim = x, nfeatures = y, balanced = TRUE))
      }))))
    })
  max.per.pc <- which.max(x = num.features[num.features < max.features])
  features <- unique(x = as.vector(x = sapply(X = dims, FUN = function(x) {
    unlist(x = TopFeatures(object = dim.reduction, dim = x, nfeatures = max.per.pc, balanced = TRUE))
    })))
  features <- unique(x = features)
  return(features)
}

TransformDataMatrix <- function(
  object,
  assay = NULL,
  new.assay.name = 'integrated',
  integration.name = 'integrated',
  features.to.integrate = NULL,
  reduction = "cca",
  do.cpp = TRUE,
  verbose = TRUE
) {
  if(verbose) {
    message("Integrating data")
  }
  assay <- assay %||% DefaultAssay(object = object)
  weights <- GetIntegrationData(
    object = object,
    integration.name = integration.name,
    slot = 'weights'
  )
  integration.matrix <- GetIntegrationData(
    object = object,
    integration.name = integration.name,
    slot = 'integration.matrix'
  )
  neighbors <- GetIntegrationData(object = object, integration.name = integration.name, slot = 'neighbors')
  nn.cells1 <- neighbors$cells1
  nn.cells2 <- neighbors$cells2

  data.use1 <- t(x = GetAssayData(
    object = object,
    assay = assay,
    slot = "data")[features.to.integrate, nn.cells1]
  )
  data.use2 <- t(x = GetAssayData(
    object = object,
    assay = assay,
    slot = "data")[features.to.integrate, nn.cells2]
  )
  if (do.cpp) {
    integrated <- IntegrateDataC(integration_matrix = as(integration.matrix, "dgCMatrix"),
                                 weights = as(weights, "dgCMatrix"),
                                 expression_cells2 = as(data.use2, "dgCMatrix"))
    dimnames(integrated) <- dimnames(data.use2)
  } else {
    bv <-  t(weights) %*% integration.matrix
    integrated <- data.use2 - bv
  }

  new.expression <- t(rbind(data.use1, integrated))
  new.expression <- new.expression[, colnames(object)]
  new.assay <- new(
    Class = 'Assay',
    counts = new(Class = "dgCMatrix"),
    data = new.expression,
    scale.data = matrix(),
    var.features = vector(),
    meta.features = data.frame(row.names = rownames(x = new.expression)),
    misc = NULL
  )
  object[[new.assay.name]] <- new.assay
  object <- LogSeuratCommand(object = object)
  command.name <- LogSeuratCommand(object = object, return.command = TRUE)
  return(object)
}
