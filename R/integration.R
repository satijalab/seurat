#' @include generics.R
#'
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Find integration anchors
#'
#' Find a set of anchors between a list of \code{\link{Seurat}} objects.
#' These  anchors can later be used to integrate the objects using the
#' \code{\link{IntegrateData}} function.
#'
#' The main steps of this procedure are outlined below. For a more detailed
#' description of the methodology, please see Stuart, Butler, et al Cell 2019:
#' \doi{10.1016/j.cell.2019.05.031}; \doi{10.1101/460147}
#'
#' First, determine anchor.features if not explicitly specified using
#' \code{\link{SelectIntegrationFeatures}}. Then for all pairwise combinations
#' of reference and query datasets:
#'
#' \itemize{
#'   \item{Perform dimensional reduction on the dataset pair as specified via
#'   the \code{reduction} parameter. If \code{l2.norm} is set to \code{TRUE},
#'   perform L2 normalization of the embedding vectors.}
#'   \item{Identify anchors - pairs of cells from each dataset
#'   that are contained within each other's neighborhoods (also known as mutual
#'   nearest neighbors).}
#'   \item{Filter low confidence anchors to ensure anchors in the low dimension
#'   space are in broad agreement with the high dimensional measurements. This
#'   is done by looking at the neighbors of each query cell in the reference
#'   dataset using \code{max.features} to define this space. If the reference
#'   cell isn't found within the first \code{k.filter} neighbors, remove the
#'   anchor.}
#'   \item{Assign each remaining anchor a score. For each anchor cell, determine
#'   the nearest \code{k.score} anchors within its own dataset and within its
#'   pair's dataset. Based on these neighborhoods, construct an overall neighbor
#'   graph and then compute the shared neighbor overlap between anchor and query
#'   cells (analogous to an SNN graph). We use the 0.01 and 0.90 quantiles on
#'   these scores to dampen outlier effects and rescale to range between 0-1.}
#' }
#'
#' @param object.list A list of \code{\link{Seurat}} objects between which to
#' find anchors for downstream integration.
#' @param assay A vector of assay names specifying which assay to use when
#' constructing anchors. If NULL, the current default assay for each object is
#' used.
#' @param reference A vector specifying the object/s to be used as a reference
#' during integration. If NULL (default), all pairwise anchors are found (no
#' reference/s). If not NULL, the corresponding objects in \code{object.list}
#' will be used as references. When using a set of specified references, anchors
#' are first found between each query and each reference. The references are
#' then integrated through pairwise integration. Each query is then mapped to
#' the integrated reference.
#' @param anchor.features Can be either:
#' \itemize{
#'   \item{A numeric value. This will call \code{\link{SelectIntegrationFeatures}}
#'   to select the provided number of features to be used in anchor finding}
#'   \item{A vector of features to be used as input to the anchor finding process}
#' }
#' @param scale Whether or not to scale the features provided. Only set to FALSE
#' if you have previously scaled the features you want to use for each object in
#' the object.list
#' @param normalization.method Name of normalization method used: LogNormalize
#' or SCT
#' @param sct.clip.range Numeric of length two specifying the min and max values
#' the Pearson residual will be clipped to
#' @param reduction Dimensional reduction to perform when finding anchors. Can
#' be one of:
#' \itemize{
#'   \item{cca: Canonical correlation analysis}
#'   \item{rpca: Reciprocal PCA}
#'   \item{jpca: Joint PCA}
#'   \item{rlsi: Reciprocal LSI}
#' }
#' @param l2.norm Perform L2 normalization on the CCA cell embeddings after
#' dimensional reduction
#' @param dims Which dimensions to use from the CCA to specify the neighbor
#' search space
#' @param k.anchor How many neighbors (k) to use when picking anchors
#' @param k.filter How many neighbors (k) to use when filtering anchors
#' @param k.score How many neighbors (k) to use when scoring anchors
#' @param max.features The maximum number of features to use when specifying the
#' neighborhood search space in the anchor filtering
#' @param nn.method Method for nearest neighbor finding. Options include: rann,
#' annoy
#' @param n.trees More trees gives higher precision when using annoy approximate
#' nearest neighbor search
#' @param eps Error bound on the neighbor finding algorithm (from RANN/Annoy)
#' @param verbose Print progress bars and output
#'
#' @return Returns an \code{\link{AnchorSet}} object that can be used as input to
#' \code{\link{IntegrateData}}.
#'
#' @references Stuart T, Butler A, et al. Comprehensive Integration of
#' Single-Cell Data. Cell. 2019;177:1888-1902 \doi{10.1016/j.cell.2019.05.031}
#'
#' @importFrom pbapply pblapply
#' @importFrom future.apply future_lapply
#' @importFrom future nbrOfWorkers
#'
#' @export
#' @concept integration
#'
#' @examples
#' \dontrun{
#' # to install the SeuratData package see https://github.com/satijalab/seurat-data
#' library(SeuratData)
#' data("panc8")
#'
#' # panc8 is a merged Seurat object containing 8 separate pancreas datasets
#' # split the object by dataset
#' pancreas.list <- SplitObject(panc8, split.by = "tech")
#'
#' # perform standard preprocessing on each object
#' for (i in 1:length(pancreas.list)) {
#'   pancreas.list[[i]] <- NormalizeData(pancreas.list[[i]], verbose = FALSE)
#'   pancreas.list[[i]] <- FindVariableFeatures(
#'     pancreas.list[[i]], selection.method = "vst",
#'     nfeatures = 2000, verbose = FALSE
#'   )
#' }
#'
#' # find anchors
#' anchors <- FindIntegrationAnchors(object.list = pancreas.list)
#'
#' # integrate data
#' integrated <- IntegrateData(anchorset = anchors)
#' }
#'
FindIntegrationAnchors <- function(
  object.list = NULL,
  assay = NULL,
  reference = NULL,
  anchor.features = 2000,
  scale = TRUE,
  normalization.method = c("LogNormalize", "SCT"),
  sct.clip.range = NULL,
  reduction = c("cca", "rpca", "jpca", "rlsi"),
  l2.norm = TRUE,
  dims = 1:30,
  k.anchor = 5,
  k.filter = 200,
  k.score = 30,
  max.features = 200,
  nn.method = "annoy",
  n.trees = 50,
  eps = 0,
  verbose = TRUE
) {
  normalization.method <- match.arg(arg = normalization.method)
  reduction <- match.arg(arg = reduction)
  if (reduction == "rpca") {
    reduction <- "pca"
  }
  if (reduction == "rlsi") {
    reduction <- "lsi"
    if (normalization.method == "SCT") {
      warning("Requested normalization method 'SCT' is not applicable for LSI")
      normalization.method <- "LogNormalize"
    }
    scale <- FALSE
    k.filter <- NA
  }
  my.lapply <- ifelse(
    test = verbose && nbrOfWorkers() == 1,
    yes = pblapply,
    no = future_lapply
  )
  object.ncells <- sapply(X = object.list, FUN = function(x) dim(x = x)[2])
  if (any(object.ncells <= max(dims))) {
    bad.obs <- which(x = object.ncells <= max(dims))
    stop("Max dimension too large: objects ", paste(bad.obs, collapse = ", "),
         " contain fewer than ", max(dims), " cells. \n Please specify a",
         " maximum dimensions that is less than the number of cells in any ",
         "object (", min(object.ncells), ").")
  }
  if (!is.null(x = assay)) {
    if (length(x = assay) != length(x = object.list)) {
      stop("If specifying the assay, please specify one assay per object in the object.list")
    }
    object.list <- sapply(
      X = 1:length(x = object.list),
      FUN = function(x) {
        DefaultAssay(object = object.list[[x]]) <- assay[x]
        return(object.list[[x]])
      }
    )
  } else {
    assay <- sapply(X = object.list, FUN = DefaultAssay)
  }
  # check tool
 object.list <- lapply(
   X = object.list,
   FUN = function (obj) {
      slot(object = obj, name = "tools")$Integration <- NULL
      return(obj)
  })
  object.list <- CheckDuplicateCellNames(object.list = object.list)
  slot <- "data"
  if (reduction == "lsi") {
    all.rownames <- lapply(X = object.list, FUN = rownames)
    anchor.features <- Reduce(f = intersect, x = all.rownames)
  }
  if (normalization.method == "SCT") {
    slot <- "scale.data"
    scale <- FALSE
    if (is.numeric(x = anchor.features)) {
      stop("Please specify the anchor.features to be used. The expected ",
      "workflow for integratinge assays produced by SCTransform is ",
      "SelectIntegrationFeatures -> PrepSCTIntegration -> ",
      "FindIntegrationAnchors.")
    }
    sct.check <- sapply(
      X = 1:length(x = object.list),
      FUN = function(x) {
        sct.cmd <- grep(
          pattern = 'PrepSCTIntegration',
          x = Command(object = object.list[[x]]),
          value = TRUE
        )
        # check assay has gone through PrepSCTIntegration
        if (!any(grepl(pattern = "PrepSCTIntegration", x = Command(object = object.list[[x]]))) ||
            Command(object = object.list[[x]], command = sct.cmd, value = "assay") != assay[x]) {
          stop("Object ", x, " assay - ", assay[x], " has not been processed ",
          "by PrepSCTIntegration. Please run PrepSCTIntegration prior to ",
          "FindIntegrationAnchors if using assays generated by SCTransform.", call. = FALSE)
        }
        # check that the correct features are being used
        if (all(Command(object = object.list[[x]], command = sct.cmd, value = "anchor.features") != anchor.features)) {
          stop("Object ", x, " assay - ", assay[x], " was processed using a ",
          "different feature set than in PrepSCTIntegration. Please rerun ",
          "PrepSCTIntegration with the same anchor.features for all objects in ",
          "the object.list.", call. = FALSE)
        }
      }
    )
  }
  if (is.numeric(x = anchor.features) && normalization.method != "SCT") {
    if (verbose) {
      message("Computing ", anchor.features, " integration features")
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
  nn.reduction <- reduction
  # if using pca or lsi, only need to compute the internal neighborhood structure once
  # for each dataset
  internal.neighbors <- list()
  if (nn.reduction %in% c("pca", "lsi","jpca")) {
    if (nn.reduction == 'jpca') {
      nn.reduction <- 'joint.pca'
      reduction <- 'joint.pca'
    }
    k.filter <- NA
    if (verbose) {
      message("Computing within dataset neighborhoods")
    }
    k.neighbor <- max(k.anchor, k.score)
    internal.neighbors <- my.lapply(
      X = 1:length(x = object.list),
      FUN = function(x) {
        NNHelper(
          data = Embeddings(object = object.list[[x]][[nn.reduction]])[, dims],
          k = k.neighbor + 1,
          method = nn.method,
          n.trees = n.trees,
          eps = eps
        )
      }
    )
  }
  # determine pairwise combinations
  combinations <- expand.grid(1:length(x = object.list), 1:length(x = object.list))
  combinations <- combinations[combinations$Var1 < combinations$Var2, , drop = FALSE]
  # determine the proper offsets for indexing anchors
  objects.ncell <- sapply(X = object.list, FUN = ncol)
  offsets <- as.vector(x = cumsum(x = c(0, objects.ncell)))[1:length(x = object.list)]
  if (is.null(x = reference)) {
    # case for all pairwise, leave the combinations matrix the same
    if (verbose) {
      message("Finding all pairwise anchors")
    }
  } else {
    reference <- unique(x = sort(x = reference))
    if (max(reference) > length(x = object.list)) {
      stop('Error: requested reference object ', max(reference), " but only ",
           length(x = object.list), " objects provided")
    }
    # modify the combinations matrix to retain only R-R and R-Q comparisons
    if (verbose) {
      message("Finding anchors between all query and reference datasets")
      ok.rows <- (combinations$Var1 %in% reference) | (combinations$Var2 %in% reference)
      combinations <- combinations[ok.rows, ]
    }
  }
  # determine all anchors
  anchoring.fxn <- function(row) {
    i <- combinations[row, 1]
    j <- combinations[row, 2]
    object.1 <- DietSeurat(
      object = object.list[[i]],
      assays = assay[i],
      features = anchor.features,
      counts = FALSE,
      scale.data = TRUE,
      dimreducs = reduction
    )
    object.2 <- DietSeurat(
      object = object.list[[j]],
      assays = assay[j],
      features = anchor.features,
      counts = FALSE,
      scale.data = TRUE,
      dimreducs = reduction
    )
    # suppress key duplication warning
    suppressWarnings(object.1[["ToIntegrate"]] <- object.1[[assay[i]]])
    DefaultAssay(object = object.1) <- "ToIntegrate"
    if (reduction %in% Reductions(object = object.1)) {
      slot(object = object.1[[reduction]], name = "assay.used") <- "ToIntegrate"
    }
    object.1 <- DietSeurat(object = object.1, assays = "ToIntegrate", scale.data = TRUE, dimreducs = reduction)
    suppressWarnings(object.2[["ToIntegrate"]] <- object.2[[assay[j]]])
    DefaultAssay(object = object.2) <- "ToIntegrate"
    if (reduction %in% Reductions(object = object.2)) {
      slot(object = object.2[[reduction]], name = "assay.used") <- "ToIntegrate"
    }
    object.2 <- DietSeurat(object = object.2, assays = "ToIntegrate", scale.data = TRUE, dimreducs = reduction)
    object.pair <- switch(
      EXPR = reduction,
      'cca' = {
        object.pair <- RunCCA(
          object1 = object.1,
          object2 = object.2,
          assay1 = "ToIntegrate",
          assay2 = "ToIntegrate",
          features = anchor.features,
          num.cc = max(dims),
          renormalize = FALSE,
          rescale = FALSE,
          verbose = verbose
        )
        if (l2.norm){
          object.pair <- L2Dim(object = object.pair, reduction = reduction)
          reduction <- paste0(reduction, ".l2")
          nn.reduction <- reduction
        }
        reduction.2 <- character()
        object.pair
      },
      'pca' = {
        object.pair <- ReciprocalProject(
          object.1 = object.1,
          object.2 = object.2,
          reduction = 'pca',
          projected.name = 'projectedpca',
          features = anchor.features,
          do.scale = FALSE,
          do.center = FALSE,
          slot = 'scale.data',
          l2.norm = l2.norm,
          verbose = verbose
        )
        reduction <- "projectedpca.ref"
        reduction.2 <- "projectedpca.query"
        if (l2.norm) {
          reduction <- paste0(reduction, ".l2")
          reduction.2 <- paste0(reduction.2, ".l2")
        }
        object.pair
      },
      'lsi' = {
        object.pair <- ReciprocalProject(
          object.1 = object.1,
          object.2 = object.2,
          reduction = 'lsi',
          projected.name = 'projectedlsi',
          features = anchor.features,
          do.center = TRUE,
          do.scale = FALSE,
          slot = 'data',
          l2.norm = l2.norm,
          verbose = verbose
        )
        reduction <- "projectedlsi.ref"
        reduction.2 <- "projectedlsi.query"
        if (l2.norm) {
          reduction <- paste0(reduction, ".l2")
          reduction.2 <- paste0(reduction.2, ".l2")
        }
        object.pair
      },
      'joint.pca' = {
        object.pair <- merge(x = object.1, y = object.2)
        reduction.2 <- "joint.pca"
        object.pair[['joint.pca']] <- CreateDimReducObject(
          embeddings = rbind(Embeddings(object.1[['joint.pca']]),
                             Embeddings(object.2[['joint.pca']])),
          loadings = Loadings(object.1[['joint.pca']]),
            key = 'Joint_',
          assay = 'ToIntegrate')
        if (l2.norm) {
          object.pair <- L2Dim(object = object.pair,
                               reduction = 'joint.pca',
                               new.dr = 'joint.pca.l2',
                               new.key = 'Jl2_'
                               )
          reduction <- paste0(reduction, ".l2")
          reduction.2 <- paste0(reduction.2, ".l2")
        }
        object.pair
      },
      stop("Invalid reduction parameter. Please choose either cca, rpca, or rlsi")
    )
    internal.neighbors <- internal.neighbors[c(i, j)]
    anchors <- FindAnchors(
      object.pair = object.pair,
      assay = c("ToIntegrate", "ToIntegrate"),
      slot = slot,
      cells1 = colnames(x = object.1),
      cells2 = colnames(x = object.2),
      internal.neighbors = internal.neighbors,
      reduction = reduction,
      reduction.2 = reduction.2,
      nn.reduction = nn.reduction,
      dims = dims,
      k.anchor = k.anchor,
      k.filter = k.filter,
      k.score = k.score,
      max.features = max.features,
      nn.method = nn.method,
      n.trees = n.trees,
      eps = eps,
      verbose = verbose
    )
    anchors[, 1] <- anchors[, 1] + offsets[i]
    anchors[, 2] <- anchors[, 2] + offsets[j]
    return(anchors)
  }
  if (nbrOfWorkers() == 1) {
    all.anchors <- pblapply(
      X = 1:nrow(x = combinations),
      FUN = anchoring.fxn
    )
  } else {
    all.anchors <- future_lapply(
      X = 1:nrow(x = combinations),
      FUN = anchoring.fxn,
      future.seed = TRUE
    )
  }
  all.anchors <- do.call(what = 'rbind', args = all.anchors)
  all.anchors <- rbind(all.anchors, all.anchors[, c(2, 1, 3)])
  all.anchors <- AddDatasetID(anchor.df = all.anchors, offsets = offsets, obj.lengths = objects.ncell)
  command <- LogSeuratCommand(object = object.list[[1]], return.command = TRUE)
  anchor.set <- new(Class = "IntegrationAnchorSet",
                    object.list = object.list,
                    reference.objects = reference %||% seq_along(object.list),
                    anchors = all.anchors,
                    offsets = offsets,
                    anchor.features = anchor.features,
                    command = command
  )
  return(anchor.set)
}

# Merge dataset and perform reciprocal SVD projection, adding new dimreducs
# for each projection and the merged original SVDs.
#
# @param object.1 First Seurat object to merge
# @param object.2 Second Seurat object to merge
# @param reduction Name of DimReduc to use. Must be an SVD-based DimReduc (eg, PCA or LSI)
# so that the loadings can be used to project new embeddings. Must be present
# in both input objects, with a substantial overlap in the features use to construct
# the SVDs.
# @param dims dimensions used for rpca
# @param projected.name Name to store projected SVDs under (eg, "projectedpca")
# @param features Features to use. Will subset the SVD loadings to use these features
# before performing projection. Typically uses the anchor.features for integration.
# @param do.center Center projected values (subtract mean)
# @param do.scale Scale projected values (divide by SD)
# @param slot Name of slot to pull data from. Should be scale.data for PCA and data for LSI
# @param verbose Display messages
# @return Returns a merged Seurat object with two projected SVDs (object.1 -> object.2, object.2 -> object.1)
# and a merged SVD (needed for within-dataset neighbors)
ReciprocalProject <- function(
  object.1,
  object.2,
  reduction,
  dims,
  projected.name,
  features,
  do.scale,
  do.center,
  slot,
  l2.norm,
  verbose = TRUE
) {
  common.features <- intersect(
    x = rownames(x = Loadings(object = object.1[[reduction]])),
    y = rownames(x = Loadings(object = object.2[[reduction]]))
  )
  common.features <- intersect(
    x = common.features,
    y = features
  )
  object.pair <- merge(x = object.1, y = object.2, merge.data = TRUE)
  data.1 <- GetAssayData(
    object = object.1,
    slot = slot
  )
  data.2 <- GetAssayData(
    object = object.2,
    slot = slot
  )

  proj.1 <- ProjectSVD(
    reduction = object.2[[reduction]],
    data = data.1,
    mode = reduction,
    features = common.features,
    do.scale = do.scale,
    do.center = do.center,
    use.original.stats = FALSE,
    verbose = verbose
  )
  proj.2 <- ProjectSVD(
    reduction = object.1[[reduction]],
    data = data.2,
    mode = reduction,
    features = common.features,
    do.scale = do.scale,
    do.center = do.center,
    use.original.stats = FALSE,
    verbose = verbose
  )
  # object.1 is ref, and object.2 is query
  reduction.dr.name.1 <- paste0(projected.name, ".ref")
  reduction.dr.name.2 <- paste0(projected.name, ".query")
  object.pair[[reduction.dr.name.1]] <- CreateDimReducObject(
    embeddings = rbind(Embeddings(object = object.1[[reduction]]), proj.2)[,dims],
    loadings = Loadings(object = object.1[[reduction]])[,dims],
    assay = DefaultAssay(object = object.1),
    key = paste0(projected.name, "ref_")
  )
  object.pair[[reduction.dr.name.2]] <- CreateDimReducObject(
    embeddings = rbind(proj.1, Embeddings(object = object.2[[reduction]]))[,dims],
    loadings = Loadings(object = object.2[[reduction]])[,dims],
    assay = DefaultAssay(object = object.2),
    key = paste0(projected.name, "query_")
  )
  object.pair[[reduction]] <- CreateDimReducObject(
    embeddings = rbind(
      Embeddings(object = object.1[[reduction]]),
      Embeddings(object = object.2[[reduction]]))[,dims],
    loadings = Loadings(object = object.1[[reduction]])[,dims],
    assay = DefaultAssay(object = object.1),
    key = paste0(projected.name, "_")
  )
  if (l2.norm) {
    slot(object = object.pair[[reduction.dr.name.1]], name = "cell.embeddings") <- Sweep(
      x = Embeddings(object = object.pair[[reduction.dr.name.1]]),
      MARGIN = 2,
      STATS = apply(X = Embeddings(object = object.pair[[reduction.dr.name.1]]), MARGIN = 2, FUN = sd),
      FUN = "/"
    )
    slot(object = object.pair[[reduction.dr.name.2]], name = "cell.embeddings") <- Sweep(
      x = Embeddings(object = object.pair[[reduction.dr.name.2]]),
      MARGIN = 2,
      STATS = apply(X = Embeddings(object = object.pair[[reduction.dr.name.2]]), MARGIN = 2, FUN = sd),
      FUN = "/"
    )

    object.pair <- L2Dim(object = object.pair, reduction = reduction.dr.name.1)
    object.pair <- L2Dim(object = object.pair, reduction = reduction.dr.name.2)
  }
  return(object.pair)
}

#' Find transfer anchors
#'
#' Find a set of anchors between a reference and query object. These
#' anchors can later be used to transfer data from the reference to
#' query object using the \code{\link{TransferData}} object.
#'
#' The main steps of this procedure are outlined below. For a more detailed
#' description of the methodology, please see Stuart, Butler, et al Cell 2019.
#' \doi{10.1016/j.cell.2019.05.031}; \doi{10.1101/460147}
#'
#' \itemize{
#'
#'   \item{Perform dimensional reduction. Exactly what is done here depends on
#'   the values set for the \code{reduction} and \code{project.query}
#'   parameters. If \code{reduction = "pcaproject"}, a PCA is performed on
#'   either the reference (if \code{project.query = FALSE}) or the query (if
#'   \code{project.query = TRUE}), using the \code{features} specified. The data
#'   from the other dataset is then projected onto this learned PCA structure.
#'   If \code{reduction = "cca"}, then CCA is performed on the reference and
#'   query for this dimensional reduction step. If
#'   \code{reduction = "lsiproject"}, the stored LSI dimension reduction in the
#'   reference object is used to project the query dataset onto the reference.
#'   If \code{l2.norm} is set to \code{TRUE}, perform L2 normalization of the
#'   embedding vectors.}
#'   \item{Identify anchors between the reference and query - pairs of cells
#'   from each dataset that are contained within each other's neighborhoods
#'   (also known as mutual nearest neighbors).}
#'   \item{Filter low confidence anchors to ensure anchors in the low dimension
#'   space are in broad agreement with the high dimensional measurements. This
#'   is done by looking at the neighbors of each query cell in the reference
#'   dataset using \code{max.features} to define this space. If the reference
#'   cell isn't found within the first \code{k.filter} neighbors, remove the
#'   anchor.}
#'   \item{Assign each remaining anchor a score. For each anchor cell, determine
#'   the nearest \code{k.score} anchors within its own dataset and within its
#'   pair's dataset. Based on these neighborhoods, construct an overall neighbor
#'   graph and then compute the shared neighbor overlap between anchor and query
#'   cells (analogous to an SNN graph). We use the 0.01 and 0.90 quantiles on
#'   these scores to dampen outlier effects and rescale to range between 0-1.}
#' }
#'
#' @param reference \code{\link{Seurat}} object to use as the reference
#' @param query \code{\link{Seurat}} object to use as the query
#' @param reference.assay Name of the Assay to use from reference
#' @param reference.neighbors Name of the Neighbor to use from the reference.
#' Optionally enables reuse of precomputed neighbors.
#' @param query.assay Name of the Assay to use from query
#' @param reduction Dimensional reduction to perform when finding anchors.
#' Options are:
#' \itemize{
#'    \item{pcaproject: Project the PCA from the reference onto the query. We
#'    recommend using PCA when reference and query datasets are from scRNA-seq}
#'    \item{lsiproject: Project the LSI from the reference onto the query. We
#'    recommend using LSI when reference and query datasets are from scATAC-seq.
#'    This requires that LSI has been computed for the reference dataset, and the
#'    same features (eg, peaks or genome bins) are present in both the reference
#'    and query. See \code{\link[Signac]{RunTFIDF}} and
#'    \code{\link[Signac]{RunSVD}}}
#'    \item{rpca: Project the PCA from the reference onto the query, and the PCA
#'    from the query onto the reference (reciprocal PCA projection).}
#'    \item{cca: Run a CCA on the reference and query }
#' }
#' @param reference.reduction Name of dimensional reduction to use from the
#' reference if running the pcaproject workflow. Optionally enables reuse of
#' precomputed reference dimensional reduction. If NULL (default), use a PCA
#' computed on the reference object.
#' @param project.query Project the PCA from the query dataset onto the
#' reference. Use only in rare cases where the query dataset has a much larger
#' cell number, but the reference dataset has a unique assay for transfer. In
#' this case, the default features will be set to the variable features of the
#' query object that are alos present in the reference.
#' @param features Features to use for dimensional reduction. If not specified,
#' set as variable features of the reference object which are also present in
#' the query.
#' @param scale Scale query data.
#' @param normalization.method Name of normalization method used: LogNormalize
#' or SCT.
#' @param recompute.residuals If using SCT as a normalization method, compute
#' query Pearson residuals using the reference SCT model parameters.
#' @param npcs Number of PCs to compute on reference if reference.reduction is
#' not provided.
#' @param l2.norm Perform L2 normalization on the cell embeddings after
#' dimensional reduction
#' @param dims Which dimensions to use from the reduction to specify the
#' neighbor search space
#' @param k.anchor How many neighbors (k) to use when finding anchors
#' @param k.filter How many neighbors (k) to use when filtering anchors. Set to
#' NA to turn off filtering.
#' @param k.score How many neighbors (k) to use when scoring anchors
#' @param max.features The maximum number of features to use when specifying the
#' neighborhood search space in the anchor filtering
#' @param nn.method Method for nearest neighbor finding. Options include: rann,
#' annoy
#' @param n.trees More trees gives higher precision when using annoy approximate
#' nearest neighbor search
#' @param eps Error bound on the neighbor finding algorithm (from
#' \code{\link{RANN}} or \code{\link{RcppAnnoy}})
#' @param approx.pca Use truncated singular value decomposition to approximate
#' PCA
#' @param mapping.score.k Compute and store nearest k query neighbors in the
#' AnchorSet object that is returned. You can optionally set this if you plan
#' on computing the mapping score and want to enable reuse of some downstream
#' neighbor calculations to make the mapping score function more efficient.
#' @param verbose Print progress bars and output
#'
#' @return Returns an \code{AnchorSet} object that can be used as input to
#' \code{\link{TransferData}}, \code{\link{IntegrateEmbeddings}} and
#' \code{\link{MapQuery}}. The dimension reduction used for finding anchors is
#' stored in the \code{AnchorSet} object and can be used for computing anchor
#' weights in downstream functions. Note that only the requested dimensions are
#' stored in the dimension reduction object in the \code{AnchorSet}. This means
#' that if \code{dims=2:20} is used, for example, the dimension of the stored
#' reduction is \code{1:19}.
#'
#' @references Stuart T, Butler A, et al. Comprehensive Integration of
#' Single-Cell Data. Cell. 2019;177:1888-1902 \doi{10.1016/j.cell.2019.05.031};
#'
#' @export
#' @importFrom methods slot slot<-
#' @concept integration
#' @examples
#' \dontrun{
#' # to install the SeuratData package see https://github.com/satijalab/seurat-data
#' library(SeuratData)
#' data("pbmc3k")
#'
#' # for demonstration, split the object into reference and query
#' pbmc.reference <- pbmc3k[, 1:1350]
#' pbmc.query <- pbmc3k[, 1351:2700]
#'
#' # perform standard preprocessing on each object
#' pbmc.reference <- NormalizeData(pbmc.reference)
#' pbmc.reference <- FindVariableFeatures(pbmc.reference)
#' pbmc.reference <- ScaleData(pbmc.reference)
#'
#' pbmc.query <- NormalizeData(pbmc.query)
#' pbmc.query <- FindVariableFeatures(pbmc.query)
#' pbmc.query <- ScaleData(pbmc.query)
#'
#' # find anchors
#' anchors <- FindTransferAnchors(reference = pbmc.reference, query = pbmc.query)
#'
#' # transfer labels
#' predictions <- TransferData(
#'   anchorset = anchors,
#'   refdata = pbmc.reference$seurat_annotations
#' )
#' pbmc.query <- AddMetaData(object = pbmc.query, metadata = predictions)
#' }
#'
FindTransferAnchors <- function(
  reference,
  query,
  normalization.method = "LogNormalize",
  recompute.residuals = TRUE,
  reference.assay = NULL,
  reference.neighbors = NULL,
  query.assay = NULL,
  reduction = "pcaproject",
  reference.reduction = NULL,
  project.query = FALSE,
  features = NULL,
  scale = TRUE,
  npcs = 30,
  l2.norm = TRUE,
  dims = 1:30,
  k.anchor = 5,
  k.filter = 200,
  k.score = 30,
  max.features = 200,
  nn.method = "annoy",
  n.trees = 50,
  eps = 0,
  approx.pca = TRUE,
  mapping.score.k = NULL,
  verbose = TRUE
) {
  op <- options(Seurat.object.assay.calcn = FALSE)
  on.exit(expr = options(op), add = TRUE)
  # input validation
  ValidateParams_FindTransferAnchors(
    reference = reference,
    query = query,
    normalization.method = normalization.method,
    recompute.residuals = recompute.residuals,
    reference.assay = reference.assay,
    reference.neighbors = reference.neighbors,
    query.assay = query.assay,
    reduction = reduction,
    reference.reduction = reference.reduction,
    project.query = project.query,
    features = features,
    scale = scale,
    npcs = npcs,
    l2.norm = l2.norm,
    dims = dims,
    k.anchor = k.anchor,
    k.filter = k.filter,
    k.score = k.score,
    max.features = max.features,
    nn.method = nn.method,
    n.trees = n.trees,
    eps = eps,
    approx.pca = approx.pca,
    mapping.score.k = mapping.score.k,
    verbose = verbose
  )
  projected <- ifelse(test = reduction == "pcaproject", yes = TRUE, no = FALSE)
  reduction.2 <- character()
  feature.mean <- NULL
  reference.reduction.init <- reference.reduction
  if (inherits(x = reference[[reference.assay]], what = 'Assay5')) {
    if (length(Layers(reference, search = "data")) > 1) {
      reference[[reference.assay]] <- JoinLayers(reference[[reference.assay]], 
                                                 layers = "data", new = "data") 
    }
  }
    if (normalization.method == "SCT") {
      if (is.null(x = reference.reduction)) {
        reference <- suppressWarnings(expr = GetResidual(
          object = reference,
          assay = reference.assay,
          features = features,
          verbose = FALSE
          ))
        features <- intersect(
          x = features,
          y = rownames(reference[[reference.assay]]$scale.data)
        )
        VariableFeatures(reference) <- features
      }
      if (IsSCT(assay = query[[query.assay]])) {
        query <- suppressWarnings(expr = GetResidual(
          object = query,
          assay = query.assay,
          features = features,
          verbose = FALSE
          ))
      }
  }
  # Rename query assay w same name as reference assay
  if (query.assay != reference.assay) {
    suppressWarnings(expr = query <- RenameAssays(
      object = query,
      assay.name = query.assay,
      new.assay.name = reference.assay,
      verbose = FALSE
      ))
    DefaultAssay(query) <- reference.assay
  }
  # only keep necessary info from objects
  suppressWarnings(
  query <- DietSeurat(
    object = query,
    assays = reference.assay,
    dimreducs = reference.reduction,
    features = features,
    scale.data = TRUE
  )
  )
  # check assay in the reference.reduction
  if (!is.null(reference.reduction) &&
    slot(object = reference[[reference.reduction]], name = "assay.used") != reference.assay) {
    warnings("reference assay is diffrent from the assay.used in", reference.reduction)
    slot(object = reference[[reference.reduction]], name = "assay.used") <- reference.assay
  }
  suppressWarnings(
    reference <- DietSeurat(
      object = reference,
      assays = reference.assay,
      dimreducs = reference.reduction,
      features = features,
      scale.data = TRUE
    )
  )
  # append query and reference to cell names - mainly to avoid name conflicts
  query <- RenameCells(
    object = query,
    new.names = paste0(Cells(x = query), "_", "query")
  )
  reference <- RenameCells(
    object = reference,
    new.names = paste0(Cells(x = reference), "_", "reference")
  )
  # Perform PCA projection
  if (reduction == 'pcaproject') {
    if (project.query) {
      if (is.null(x = reference.reduction)) {
        reference.reduction <- "pca"
        if (verbose) {
          message("Performing PCA on the provided query using ", length(x = features), " features as input.")
        }
        if (normalization.method == "LogNormalize") {
          query <- ScaleData(
            object = query,
            features = features,
            do.scale = scale,
            verbose = FALSE
          )
        }
        query <- RunPCA(
          object = query,
          npcs = npcs,
          reduction.name = reference.reduction,
          verbose = FALSE,
          features = features,
          approx = approx.pca
        )
      }
      projected.pca <- ProjectCellEmbeddings(
        reference = query,
        reduction = reference.reduction,
        query = reference,
        scale = scale,
        dims = dims,
        verbose = verbose
      )
      orig.embeddings <- Embeddings(object = query[[reference.reduction]])[, dims]
      orig.loadings <- Loadings(object = query[[reference.reduction]])
    } else {
      if (is.null(x = reference.reduction)) {
        reference.reduction <- "pca"
        if (verbose) {
          message("Performing PCA on the provided reference using ", length(x = features), " features as input.")
        }
        if (normalization.method == "LogNormalize") {
          reference <- ScaleData(object = reference, features = features, do.scale = scale, verbose = FALSE)
        }
        reference <- RunPCA(
          object = reference,
          npcs = npcs,
          verbose = FALSE,
          features = features,
          approx = approx.pca
        )
      }
      query_nCount_UMI <- query[[]][, paste0("nCount_", query.assay)]
      names(x = query_nCount_UMI) <- colnames(x = query)
      projected.pca <- ProjectCellEmbeddings(
         reference = reference,
         reduction = reference.reduction,
         normalization.method = normalization.method,
         query = query,
         scale = scale,
         dims = dims,
         nCount_UMI = query_nCount_UMI,
         feature.mean = feature.mean,
         verbose = verbose
       )
      orig.embeddings <- Embeddings(object = reference[[reference.reduction]])[, dims]
      orig.loadings <- Loadings(object = reference[[reference.reduction]])
    }
    combined.pca <- CreateDimReducObject(
      embeddings = as.matrix(x = rbind(orig.embeddings, projected.pca)),
      key = "ProjectPC_",
      assay = reference.assay
    )
    combined.ob <- suppressWarnings(expr = merge(
      x = DietSeurat(object = reference, counts = FALSE),
      y = DietSeurat(object = query, counts = FALSE),
    ))
    combined.ob[["pcaproject"]] <- combined.pca
    colnames(x = orig.loadings) <- paste0("ProjectPC_", 1:ncol(x = orig.loadings))
    Loadings(object = combined.ob[["pcaproject"]]) <- orig.loadings[, dims]
  }
  # Use reciprocal PCA projection in anchor finding
  if (reduction == "rpca") {
    # Run PCA on reference and query
    if (is.null(x = reference.reduction)) {
      reference.reduction <- "pca"
      if (verbose) {
        message("Performing PCA on the provided reference using ", length(x = features), " features as input.")
      }
      if (normalization.method == "LogNormalize") {
        reference <- ScaleData(
          object = reference,
          features = features,
          do.scale = scale,
          verbose = verbose
        )
      }
      reference <- RunPCA(
        object = reference,
        npcs = npcs,
        verbose = FALSE,
        features = features,
        approx = approx.pca
      )
    }
    if (verbose) {
      message("Performing PCA on the provided query using ", length(x = features), " features as input.")
    }
    if (normalization.method == "LogNormalize") {
      query <- ScaleData(
        object = query,
        features = features,
        do.scale = scale,
        verbose = verbose
      )
    }
    query <- RunPCA(
      object = query,
      npcs =  ncol(x = reference[[reference.reduction]]),
      reduction.name = reference.reduction,
      verbose = FALSE,
      features = features,
      approx = approx.pca
    )
    combined.ob <- ReciprocalProject(
      object.1 = reference,
      object.2 = query,
      reduction = reference.reduction,
      dims = dims,
      projected.name = reduction,
      features = features,
      do.scale = FALSE,
      do.center = FALSE,
      slot = 'scale.data',
      l2.norm = l2.norm,
      verbose = verbose
    )
    # pcaproject is used as the weight.matrix in MapQuery
    projected.pca <- ProjectCellEmbeddings(
      reference = reference,
      reduction = reference.reduction,
      query = query,
      scale = scale,
      normalization.method = normalization.method,
      dims = dims,
      feature.mean = feature.mean,
      verbose = verbose
    )
    orig.embeddings <- Embeddings(object = reference[[reference.reduction]])[, dims]
    orig.loadings <- Loadings(object = reference[[reference.reduction]])
  combined.pca <- CreateDimReducObject(
    embeddings = as.matrix(x = rbind(orig.embeddings, projected.pca)),
    key = "ProjectPC_",
    assay = reference.assay
  )
  combined.ob[["pcaproject"]] <- combined.pca
  colnames(x = orig.loadings) <- paste0("ProjectPC_", 1:ncol(x = orig.loadings))
  Loadings(object = combined.ob[["pcaproject"]]) <- orig.loadings[, dims]

    if (l2.norm) {
      # L2 norm is done on each projected PCA in ReciprocalProject, so turn it off here
      # avoids later error as we now have two reductions (rpca.ref and rpca.query)
      l2.norm <- FALSE
      reduction <- "rpca.ref.l2"
      reduction.2 <- "rpca.query.l2"
    } else {
      reduction <- "rpca.ref"
      reduction.2 <- "rpca.query"
    }
  if (project.query) {
    reduction <- gsub(".ref", ".query", reduction)
    reduction.2 <- gsub(".query", ".ref", reduction.2)
  }
  }
  # Run CCA as dimension reduction to be used in anchor finding
  if (reduction == 'cca') {
    if (normalization.method == "LogNormalize") {
      reference <- ScaleData(object = reference, features = features, do.scale = scale, verbose = FALSE)
      query <- ScaleData(object = query, features = features, do.scale = scale, verbose = FALSE)
    }
    combined.ob <- RunCCA(
      object1 = reference,
      object2 = query,
      features = features,
      num.cc = max(dims),
      renormalize = FALSE,
      rescale = FALSE,
      verbose = verbose
    )
    slot(object = combined.ob[["cca"]], name = "cell.embeddings") <- Embeddings(combined.ob[["cca"]])[, dims]
    slot(object = combined.ob[["cca"]], name = "feature.loadings") <- Loadings(combined.ob[["cca"]])[, dims]
    slot(object = combined.ob[["cca"]], name = "feature.loadings.projected") <- Loadings(object = combined.ob[["cca"]], projected = TRUE)[, dims]
  }
  if (reduction == "lsiproject") {
    if (project.query) {
      projected.lsi <- ProjectSVD(
        reduction = query[[reference.reduction]],
        data = GetAssayData(object = reference, assay = reference.assay, slot = "data"),
        mode = "lsi",
        do.center = FALSE,
        do.scale = FALSE,
        use.original.stats = FALSE,
        verbose = verbose
      )
      orig.embeddings <- Embeddings(object = query[[reference.reduction]])
      orig.loadings <- Loadings(object = query[[reference.reduction]])
    } else {
      projected.lsi <- ProjectSVD(
        reduction = reference[[reference.reduction]],
        data = GetAssayData(object = query, assay = reference.assay, slot = "data"),
        mode = "lsi",
        do.center = FALSE,
        do.scale = FALSE,
        use.original.stats = FALSE,
        verbose = verbose
      )
      orig.embeddings <- Embeddings(object = reference[[reference.reduction]])
      orig.loadings <- Loadings(object = reference[[reference.reduction]])
    }
    combined.lsi <- CreateDimReducObject(
      embeddings = as.matrix(x = rbind(orig.embeddings, projected.lsi))[,dims],
      key = "ProjectLSI_",
      assay = reference.assay
    )
    combined.ob <- merge(
      x = DietSeurat(object = reference),
      y = DietSeurat(object = query)
    )
    combined.ob[["lsiproject"]] <- combined.lsi
    colnames(x = orig.loadings) <- paste0("ProjectLSI_", 1:ncol(x = orig.loadings))
    Loadings(object = combined.ob[["lsiproject"]]) <- orig.loadings[,dims]
  }
  if (l2.norm) {
    combined.ob <- L2Dim(object = combined.ob, reduction = reduction)
    reduction <- paste0(reduction, ".l2")
  }
  precomputed.neighbors <- list(ref.neighbors = NULL, query.neighbors = NULL)
  nn.idx1 <- NULL
  nn.idx2 <- NULL
  # if computing the mapping score later, compute large enough query
  # neighborhood here to reuse
  if (!is.null(x = mapping.score.k)) {
    if (verbose) {
      message("Finding query neighbors")
    }
    k.nn <- max(k.score, k.anchor)
    query.neighbors <- NNHelper(
      data = Embeddings(object = combined.ob[[reduction]])[colnames(x = query), ],
      k = max(mapping.score.k, k.nn + 1),
      method = nn.method,
      n.trees = n.trees,
      cache.index = TRUE
    )
    query.neighbors.sub <- query.neighbors
    slot(object = query.neighbors.sub, name = "nn.idx") <- slot(
      object = query.neighbors.sub,
      name = "nn.idx")[, 1:(k.nn + 1)]
    slot(object = query.neighbors.sub, name = "nn.dist") <- slot(
      object = query.neighbors.sub,
      name = "nn.dist")[, 1:(k.nn + 1)]
    precomputed.neighbors[["query.neighbors"]] <- query.neighbors.sub
    nn.idx2 <- Index(object = query.neighbors.sub)
  }
  if (!is.null(x = reference.neighbors)) {
    precomputed.neighbors[["ref.neighbors"]] <- reference[[reference.neighbors]]
  } else {
    precomputed.neighbors[["ref.neighbors"]] <- NNHelper(
      data = Embeddings(combined.ob[[reduction]])[
        colnames(x = reference),
        1:length(x = dims)
        ],
      k = max(k.score, k.anchor) + 1,
      method = nn.method,
      cache.index = TRUE
      )
  }
  nn.idx1 <- Index(object = precomputed.neighbors[["ref.neighbors"]])
  anchors <- FindAnchors(
    object.pair = combined.ob,
    assay = c(reference.assay, reference.assay),
    slot = "data",
    cells1 = colnames(x = reference),
    cells2 = colnames(x = query),
    reduction = reduction,
    reduction.2 = reduction.2,
    internal.neighbors = precomputed.neighbors,
    dims = 1:length(x = dims),
    k.anchor = k.anchor,
    k.filter = k.filter,
    k.score = k.score,
    max.features = max.features,
    nn.method = nn.method,
    n.trees = n.trees,
    nn.idx1 = nn.idx1,
    nn.idx2 = nn.idx2,
    eps = eps,
    projected = projected,
    verbose = verbose
  )
  reductions <- slot(object = combined.ob, name = "reductions")
  for (i in unique(x = c(reference.assay))) {
    dummy.assay <- paste0(i, "DUMMY")
    suppressWarnings(
      expr = combined.ob[[dummy.assay]] <- CreateDummyAssay(assay = combined.ob[[i]])
    )
    DefaultAssay(combined.ob) <- dummy.assay
    combined.ob[[i]] <- NULL
    suppressWarnings(
      expr = combined.ob[[i]] <- combined.ob[[dummy.assay]]
    )
    DefaultAssay(combined.ob) <- i
    combined.ob[[dummy.assay]] <- NULL
  }
  slot(object = combined.ob, name = "reductions") <- reductions
  command <- LogSeuratCommand(object = combined.ob, return.command = TRUE)
  slot(command, name = 'params')$reference.reduction <- reference.reduction.init
  anchor.set <- new(
    Class = "TransferAnchorSet",
    object.list = list(combined.ob),
    reference.cells = colnames(x = reference),
    query.cells = colnames(x = query),
    anchors = anchors,
    anchor.features = features,
    command = command
  )
  if (!is.null(x = precomputed.neighbors[["query.neighbors"]])) {
    slot(object = anchor.set, name = "neighbors") <- list(
      query.neighbors = query.neighbors)
  }
  return(anchor.set)
}

#' Get the predicted identity
#'
#' Utility function to easily pull out the name of the class with the maximum
#' prediction. This is useful if you've set \code{prediction.assay = TRUE} in
#' \code{\link{TransferData}} and want to have a vector with the predicted class.
#'
#' @param object Seurat object
#' @param assay Name of the assay holding the predictions
#' @param slot Slot of the assay in which the prediction scores are stored
#' @param score.filter Return "Unassigned" for any cell with a score less than
#' this value
#'
#' @return Returns a vector of predicted class names
#'
#' @examples
#' \dontrun{
#'   prediction.assay <- TransferData(anchorset = anchors, refdata = reference$class)
#'   query[["predictions"]] <- prediction.assay
#'   query$predicted.id <- GetTransferPredictions(query)
#' }
#' @export
#' @concept integration
#'
GetTransferPredictions <- function(object, assay = "predictions", slot = "data", score.filter = 0.75) {
  dat <- GetAssayData(object[[assay]], slot = slot)
  predictions <- apply(
    X = dat,
    MARGIN = 2,
    FUN = function(x){
      if (x['max'] < score.filter) {
        "Unassigned"
      } else {
        x <- x[-which(x = names(x = x) == "max")]
        names(x = which.max(x = x))
      }
    }
  )
  return(predictions)
}

#' Integrate data
#'
#' Perform dataset integration using a pre-computed \code{\link{AnchorSet}}.
#'
#' The main steps of this procedure are outlined below. For a more detailed
#' description of the methodology, please see Stuart, Butler, et al Cell 2019.
#' \doi{10.1016/j.cell.2019.05.031}; \doi{10.1101/460147}
#'
#' For pairwise integration:
#'
#' \itemize{
#'   \item{Construct a weights matrix that defines the association between each
#'   query cell and each anchor. These weights are computed as 1 - the distance
#'   between the query cell and the anchor divided by the distance of the query
#'   cell to the \code{k.weight}th anchor multiplied by the anchor score
#'   computed in \code{\link{FindIntegrationAnchors}}. We then apply a Gaussian
#'   kernel width a bandwidth defined by \code{sd.weight} and normalize across
#'   all \code{k.weight} anchors.}
#'   \item{Compute the anchor integration matrix as the difference between the
#'   two expression matrices for every pair of anchor cells}
#'   \item{Compute the transformation matrix as the product of the integration
#'   matrix and the weights matrix.}
#'   \item{Subtract the transformation matrix from the original expression
#'   matrix.}
#' }
#'
#' For multiple dataset integration, we perform iterative pairwise integration.
#' To determine the order of integration (if not specified via
#' \code{sample.tree}), we
#' \itemize{
#'   \item{Define a distance between datasets as the total number of cells in
#'   the smaller dataset divided by the total number of anchors between the two
#'   datasets.}
#'   \item{Compute all pairwise distances between datasets}
#'   \item{Cluster this distance matrix to determine a guide tree}
#' }
#'
#'
#' @param anchorset An \code{\link{AnchorSet}} object generated by
#' \code{\link{FindIntegrationAnchors}}
#' @param new.assay.name Name for the new assay containing the integrated data
#' @param normalization.method Name of normalization method used: LogNormalize
#' or SCT
#' @param features Vector of features to use when computing the PCA to determine
#' the weights. Only set if you want a different set from those used in the
#' anchor finding process
#' @param features.to.integrate Vector of features to integrate. By default,
#' will use the features used in anchor finding.
#' @param dims Number of dimensions to use in the anchor weighting procedure
#' @param k.weight Number of neighbors to consider when weighting anchors
#' @param weight.reduction Dimension reduction to use when calculating anchor
#' weights. This can be one of:
#' \itemize{
#'    \item{A string, specifying the name of a dimension reduction present in
#'    all objects to be integrated}
#'    \item{A vector of strings, specifying the name of a dimension reduction to
#'    use for each object to be integrated}
#'    \item{A vector of \code{\link{DimReduc}} objects, specifying the object to
#'    use for each object in the integration}
#'    \item{NULL, in which case a new PCA will be calculated and used to
#'    calculate anchor weights}
#' }
#' Note that, if specified, the requested dimension reduction will only be used
#' for calculating anchor weights in the first merge between reference and
#' query, as the merged object will subsequently contain more cells than was in
#' query, and weights will need to be calculated for all cells in the object.
#' @param sd.weight Controls the bandwidth of the Gaussian kernel for weighting
#' @param sample.tree Specify the order of integration. Order of integration
#' should be encoded in a matrix, where each row represents one of the pairwise
#' integration steps. Negative numbers specify a dataset, positive numbers
#' specify the integration results from a given row (the format of the merge
#' matrix included in the \code{\link{hclust}} function output). For example:
#' `matrix(c(-2, 1, -3, -1), ncol = 2)` gives:
#'
#' ```
#'             [,1]  [,2]
#'        [1,]   -2   -3
#'        [2,]    1   -1
#' ```
#'
#' Which would cause dataset 2 and 3 to be integrated first, then the resulting
#' object integrated with dataset 1.
#'
#'  If NULL, the sample tree will be computed automatically.
#' @param preserve.order Do not reorder objects based on size for each pairwise
#' integration.
#' @param eps Error bound on the neighbor finding algorithm (from
#' \code{\link{RANN}})
#' @param verbose Print progress bars and output
#'
#' @return Returns a \code{\link{Seurat}} object with a new integrated
#' \code{\link{Assay}}. If \code{normalization.method = "LogNormalize"}, the
#' integrated data is returned to the \code{data} slot and can be treated as
#' log-normalized, corrected data. If \code{normalization.method = "SCT"}, the
#' integrated data is returned to the \code{scale.data} slot and can be treated
#' as centered, corrected Pearson residuals.
#'
#' @references Stuart T, Butler A, et al. Comprehensive Integration of
#' Single-Cell Data. Cell. 2019;177:1888-1902 \doi{10.1016/j.cell.2019.05.031}
#'
#' @export
#' @concept integration
#' @md
#' @examples
#' \dontrun{
#' # to install the SeuratData package see https://github.com/satijalab/seurat-data
#' library(SeuratData)
#' data("panc8")
#'
#' # panc8 is a merged Seurat object containing 8 separate pancreas datasets
#' # split the object by dataset
#' pancreas.list <- SplitObject(panc8, split.by = "tech")
#'
#' # perform standard preprocessing on each object
#' for (i in 1:length(pancreas.list)) {
#'   pancreas.list[[i]] <- NormalizeData(pancreas.list[[i]], verbose = FALSE)
#'   pancreas.list[[i]] <- FindVariableFeatures(
#'     pancreas.list[[i]], selection.method = "vst",
#'     nfeatures = 2000, verbose = FALSE
#'   )
#' }
#'
#' # find anchors
#' anchors <- FindIntegrationAnchors(object.list = pancreas.list)
#'
#' # integrate data
#' integrated <- IntegrateData(anchorset = anchors)
#' }
#'
IntegrateData <- function(
  anchorset,
  new.assay.name = "integrated",
  normalization.method = c("LogNormalize", "SCT"),
  features = NULL,
  features.to.integrate = NULL,
  dims = 1:30,
  k.weight = 100,
  weight.reduction = NULL,
  sd.weight = 1,
  sample.tree = NULL,
  preserve.order = FALSE,
  eps = 0,
  verbose = TRUE
) {
  normalization.method <- match.arg(arg = normalization.method)
  reference.datasets <- slot(object = anchorset, name = 'reference.objects')
  object.list <- slot(object = anchorset, name = 'object.list')
  anchors <- slot(object = anchorset, name = 'anchors')
  ref <- object.list[reference.datasets]
  features <- features %||% slot(object = anchorset, name = "anchor.features")

  unintegrated <- suppressWarnings(expr = merge(
    x = object.list[[1]],
    y = object.list[2:length(x = object.list)]
  ))
  if (!is.null(x = features.to.integrate)) {
    features.to.integrate <- intersect(
      x = features.to.integrate,
      y = Reduce(
        f = intersect,
        x = lapply(
          X = object.list,
          FUN = rownames
        )
      )
    )
  }
  if (normalization.method == "SCT") {
    model.list <- list()
    for (i in 1:length(x = object.list)) {
      assay <- DefaultAssay(object = object.list[[i]])
      if (length(x = setdiff(x = features.to.integrate, y = features)) != 0) {
        object.list[[i]] <- GetResidual(
          object = object.list[[i]],
          features = setdiff(x = features.to.integrate, y = features),
          verbose = verbose
        )
      }
      print(i)
      model.list[[i]] <- slot(object = object.list[[i]][[assay]], name = "SCTModel.list")
      object.list[[i]][[assay]] <- suppressWarnings(expr = CreateSCTAssayObject(
        data = GetAssayData(
          object = object.list[[i]],
          assay = assay,
          slot = "scale.data")
        )
      )
    }
    model.list <- unlist(x = model.list)
    slot(object = anchorset, name = "object.list") <- object.list
  }
  # perform pairwise integration of reference objects
  reference.integrated <- PairwiseIntegrateReference(
    anchorset = anchorset,
    new.assay.name = new.assay.name,
    normalization.method = normalization.method,
    features = features,
    features.to.integrate = features.to.integrate,
    dims = dims,
    k.weight = k.weight,
    weight.reduction = weight.reduction,
    sd.weight = sd.weight,
    sample.tree = sample.tree,
    preserve.order = preserve.order,
    eps = eps,
    verbose = verbose
  )

  # set SCT model
  if (normalization.method == "SCT") {
    if (is.null(x = Tool(object = reference.integrated, slot = "Integration"))) {
      reference.sample <- slot(object = anchorset, name = "reference.objects")
    } else {
      reference.sample <- SampleIntegrationOrder(
        tree = slot(
          object = reference.integrated,
          name = "tools"
        )$Integration@sample.tree
      )[1]
    }
    reference.cells <- Cells(x = object.list[[reference.sample]])
    reference.model <- NULL
    if (length(x = model.list) > 0) {
      reference.model <- sapply(X = model.list, FUN = function(model) {
        reference.check <- FALSE
        model.cells <- Cells(x = model)
        if (length(x = model.cells) > 0 &
            length(x = setdiff(x = model.cells, y = reference.cells)) == 0) {
          reference.check <- TRUE
        }
        return(reference.check)
        }
      )
      reference.model <- model.list[[which(reference.model)]]
    }
  }
  if (length(x = reference.datasets) == length(x = object.list)) {
    if (normalization.method == "SCT") {
      reference.integrated[[new.assay.name]] <- CreateSCTAssayObject(
        data = GetAssayData(object = reference.integrated, assay = new.assay.name, slot = "data"),
        scale.data = ScaleData(
          object = GetAssayData(object = reference.integrated, assay = new.assay.name, slot = "scale.data"),
          do.scale = FALSE,
          do.center = TRUE,
          verbose = FALSE),
        SCTModel.list = reference.model
      )
      levels(x =  reference.integrated[[new.assay.name]]) <- "refmodel"
      reference.integrated[[assay]] <- unintegrated[[assay]]
    }
    DefaultAssay(object = reference.integrated) <- new.assay.name
    VariableFeatures(object = reference.integrated) <- features
    reference.integrated[["FindIntegrationAnchors"]] <- slot(object = anchorset, name = "command")
    reference.integrated <- suppressWarnings(LogSeuratCommand(object = reference.integrated))
    return(reference.integrated)
  } else {
    active.assay <- DefaultAssay(object = ref[[1]])
    reference.integrated[[active.assay]] <- NULL
    reference.integrated[[active.assay]] <- CreateAssayObject(
      data = GetAssayData(
        object = reference.integrated[[new.assay.name]],
        slot = 'data'
      ),
      check.matrix = FALSE
    )
    DefaultAssay(object = reference.integrated) <- active.assay
    reference.integrated[[new.assay.name]] <- NULL
    VariableFeatures(object = reference.integrated) <- features
    # Extract the query objects (if any) and map to reference
    integrated.data <- MapQueryData(
      anchorset = anchorset,
      reference = reference.integrated,
      new.assay.name = new.assay.name,
      normalization.method = normalization.method,
      features = features,
      features.to.integrate = features.to.integrate,
      dims = dims,
      k.weight = k.weight,
      weight.reduction = weight.reduction,
      sd.weight = sd.weight,
      preserve.order = preserve.order,
      eps = eps,
      verbose = verbose
    )

    # Construct final assay object
    integrated.assay <- CreateAssayObject(
      data = integrated.data,
      check.matrix = FALSE
    )
    if (normalization.method == "SCT") {
      integrated.assay <- CreateSCTAssayObject(
        data =  integrated.data,
        scale.data = ScaleData(
          object = integrated.data,
          do.scale = FALSE,
          do.center = TRUE,
          verbose = FALSE),
        SCTModel.list = reference.model
      )
      levels(x = integrated.assay) <- "refmodel"
    }
    unintegrated[[new.assay.name]] <- integrated.assay
    unintegrated <- SetIntegrationData(
      object = unintegrated,
      integration.name = "Integration",
      slot = "anchors",
      new.data = anchors
    )
    if (!is.null(x = Tool(object = reference.integrated, slot = "Integration"))) {
      sample.tree <- GetIntegrationData(
        object = reference.integrated,
        integration.name = "Integration",
        slot = "sample.tree"
      )
    }
    unintegrated <- SetIntegrationData(
      object = unintegrated,
      integration.name = "Integration",
      slot = "sample.tree",
      new.data = sample.tree
    )
    DefaultAssay(object = unintegrated) <- new.assay.name
    VariableFeatures(object = unintegrated) <- features
    unintegrated[["FindIntegrationAnchors"]] <- slot(object = anchorset, name = "command")
    unintegrated <- suppressWarnings(LogSeuratCommand(object = unintegrated))
    return(unintegrated)
  }
}

#' @inheritParams IntegrateData
#'
#' @rdname IntegrateEmbeddings
#' @concept integration
#' @export
#' @method IntegrateEmbeddings IntegrationAnchorSet
#'
IntegrateEmbeddings.IntegrationAnchorSet <- function(
  anchorset,
  new.reduction.name = "integrated_dr",
  reductions = NULL,
  dims.to.integrate = NULL,
  k.weight = 100,
  weight.reduction = NULL,
  sd.weight = 1,
  sample.tree = NULL,
  preserve.order = FALSE,
  verbose = TRUE,
  ...
) {
  CheckDots(...)
  reference.datasets <- slot(object = anchorset, name = 'reference.objects')
  object.list <- slot(object = anchorset, name = 'object.list')
  anchors <- slot(object = anchorset, name = 'anchors')
  reductions <- reductions %||%  slot(
    object = anchorset,
    name = 'weight.reduction'
  )
  ValidateParams_IntegrateEmbeddings_IntegrationAnchors(
    anchorset = anchorset,
    object.list = object.list,
    reductions = reductions,
    dims.to.integrate = dims.to.integrate,
    k.weight = k.weight,
    weight.reduction = weight.reduction,
    sample.tree = sample.tree
  )
  unintegrated <- merge(
    x = object.list[[1]],
    y = object.list[2:length(x = object.list)]
  )
  # make DimReducs into Assays temporarily
  intdr.assay <- DefaultAssay(object = reductions)
  int.assay <- DefaultAssay(object = object.list[[1]])
  dims.names <- paste0("drtointegrate-", dims.to.integrate)
  # cell.names.map <- Cells(x = unintegrated)
  cell.names.map <- colnames(x = unintegrated)
  names(x = cell.names.map) <- make.unique(names = unname(obj = do.call(
    what = c,
    args = lapply(X = object.list, FUN = colnames)))
  )
  for (i in 1:length(x = object.list)) {
    embeddings <- t(x = Embeddings(object = reductions)[cell.names.map[Cells(x = object.list[[i]])], dims.to.integrate])
    rownames(x = embeddings) <- dims.names
    fake.assay <- suppressWarnings(
      expr = CreateAssayObject(
        data = embeddings,
        check.matrix = FALSE
      )
    )
    object.list[[i]][['drtointegrate']] <- fake.assay
    DefaultAssay(object = object.list[[i]]) <- "drtointegrate"
  }
  slot(object = anchorset, name = "object.list") <- object.list
  new.reduction.name.safe <- gsub(pattern = "_", replacement = "", x = new.reduction.name)
  new.reduction.name.safe <- gsub(pattern = "[.]", replacement = "", x = new.reduction.name.safe)

  reference.integrated <- PairwiseIntegrateReference(
    anchorset = anchorset,
    new.assay.name = new.reduction.name.safe,
    normalization.method = "LogNormalize",
    features = dims.names,
    features.to.integrate = dims.names,
    dims = NULL,
    k.weight = k.weight,
    weight.reduction = weight.reduction,
    sd.weight = sd.weight,
    sample.tree = sample.tree,
    preserve.order = preserve.order,
    verbose = verbose
  )
  if (length(x = reference.datasets) == length(x = object.list)) {
    reference.dr <- CreateDimReducObject(
      embeddings = as.matrix(x = t(GetAssayData(
        object = reference.integrated[[new.reduction.name.safe]]
      ))),
      assay = intdr.assay,
      loadings = Loadings(object = reductions),
      key = paste0(new.reduction.name.safe, "_")
    )
    DefaultAssay(object = reference.integrated) <- int.assay
    reference.integrated[["drtointegrate"]] <- NULL
    reference.integrated[[new.reduction.name.safe]] <- NULL
    reference.integrated[[new.reduction.name]] <- reference.dr
    return(reference.integrated)
  }
  active.assay <- DefaultAssay(object = object.list[reference.datasets][[1]])
  reference.integrated[[active.assay]] <- NULL
  reference.integrated[[active.assay]] <- CreateAssayObject(
    data = GetAssayData(
      object = reference.integrated[[new.reduction.name.safe]],
      slot = 'data'
    )
  )
  DefaultAssay(object = reference.integrated) <- active.assay
  reference.integrated[[new.reduction.name.safe]] <- NULL
  VariableFeatures(object = reference.integrated) <- dims.names
  # Extract the query objects (if any) and map to reference
  integrated.data <- MapQueryData(
    anchorset = anchorset,
    reference = reference.integrated,
    new.assay.name = new.reduction.name.safe,
    normalization.method = "LogNormalize",
    features = dims.names,
    features.to.integrate = dims.names,
    dims = NULL,
    k.weight = k.weight,
    weight.reduction = weight.reduction,
    sd.weight = sd.weight,
    preserve.order = preserve.order,
    verbose = verbose
  )
  suppressWarnings(expr = unintegrated[[new.reduction.name]] <- CreateDimReducObject(
    embeddings = as.matrix(x = t(x = integrated.data)),
    assay = intdr.assay,
    loadings = Loadings(object = reductions),
    key = paste0(new.reduction.name.safe, "_")
  ))
  unintegrated <- SetIntegrationData(
    object = unintegrated,
    integration.name = "Integration",
    slot = "anchors",
    new.data = anchors
  )
  if (!is.null(x = Tool(object = reference.integrated, slot = "Integration"))) {
    sample.tree <- GetIntegrationData(
      object = reference.integrated,
      integration.name = "Integration",
      slot = "sample.tree"
    )
  }
  unintegrated <- SetIntegrationData(
    object = unintegrated,
    integration.name = "Integration",
    slot = "sample.tree",
    new.data = sample.tree
  )
  unintegrated[["FindIntegrationAnchors"]] <- slot(object = anchorset, name = "command")
  suppressWarnings(unintegrated <- LogSeuratCommand(object = unintegrated))
  return(unintegrated)
}
#' @param reference Reference object used in anchorset construction
#' @param query Query object used in anchorset construction
#' @param reuse.weights.matrix Can be used in conjunction with the store.weights
#' parameter in TransferData to reuse a precomputed weights matrix.
#'
#' @rdname IntegrateEmbeddings
#' @concept integration
#' @export
#' @method IntegrateEmbeddings TransferAnchorSet
#'
IntegrateEmbeddings.TransferAnchorSet <- function(
  anchorset,
  reference,
  query,
  new.reduction.name = "integrated_dr",
  reductions = 'pcaproject',
  dims.to.integrate = NULL,
  k.weight = 100,
  weight.reduction = NULL,
  reuse.weights.matrix = TRUE,
  sd.weight = 1,
  preserve.order = FALSE,
  verbose = TRUE,
  ...
) {
  CheckDots(...)
  combined.object <- slot(object = anchorset, name = 'object.list')[[1]]
  anchors <- slot(object = anchorset, name = 'anchors')
  weights.matrix <- NULL
  ValidateParams_IntegrateEmbeddings_TransferAnchors(
    anchorset = anchorset,
    combined.object = combined.object,
    reference = reference,
    query = query,
    reductions = reductions,
    dims.to.integrate = dims.to.integrate,
    k.weight = k.weight,
    weight.reduction = weight.reduction,
    reuse.weights.matrix = reuse.weights.matrix
  )
  object.list <- list(reference, query)
  # make DimReducs into Assays temporarily
  intdr.assay <- DefaultAssay(object = object.list[[1]][[reductions[[1]]]])
  int.assay <- DefaultAssay(object = object.list[[1]])
  dims.names <- paste0("drtointegrate-", dims.to.integrate)
  for (i in 1:length(x = object.list)) {
    embeddings <- t(x = Embeddings(
      object = object.list[[i]], reduction = reductions[[i]]
    )[ , dims.to.integrate])
    rownames(x = embeddings) <- dims.names
    fake.assay <- suppressWarnings(
      expr = CreateAssayObject(
        data = embeddings,
        check.matrix = FALSE
      )
    )
    object.list[[i]][['drtointegrate']] <- fake.assay
    DefaultAssay(object = object.list[[i]]) <- "drtointegrate"
    object.list[[i]] <- DietSeurat(object = object.list[[i]], assays = "drtointegrate")
  }
  slot(object = anchorset, name = "object.list") <- object.list
  new.reduction.name.safe <- gsub(pattern = "_", replacement = "", x = new.reduction.name)
  new.reduction.name.safe <- gsub(pattern = "[.]", replacement = "", x = new.reduction.name)
  slot(object = anchorset, name = "reference.objects") <- 1
  anchors <- as.data.frame(x = anchors)
  anchors$dataset1 <- 1
  anchors$dataset2 <- 2
  slot(object = anchorset, name = "anchors") <- anchors
  integrated.embeddings <- MapQueryData(
    anchorset = anchorset,
    reference = object.list[[1]],
    new.assay.name = new.reduction.name.safe,
    normalization.method = "LogNormalize",
    features = dims.names,
    features.to.integrate = dims.names,
    dims = NULL,
    k.weight = k.weight,
    weight.reduction = weight.reduction,
    weights.matrix = weights.matrix,
    no.offset = TRUE,
    sd.weight = sd.weight,
    preserve.order = preserve.order,
    verbose = verbose
  )
  integrated.embeddings <- as.matrix(x = integrated.embeddings)
  query[[new.reduction.name]]  <- CreateDimReducObject(
    embeddings = t(x = integrated.embeddings[, Cells(x = query)]),
    assay = DefaultAssay(object = query[[reductions[1]]]),
    key = paste0(new.reduction.name.safe, "_")
  )
  query <- RenameCells(
    object = query,
    new.names = gsub(pattern = "_query$", replacement = "", x = Cells(x = query))
  )
  query[[reductions[[1]]]] <- NULL
  return(query)
}


#' Integrate embeddings from the integrated sketched.assay
#'
#' The main steps of this procedure are outlined below. For a more detailed
#' description of the methodology, please see Hao,  et al Biorxiv 2022:
#' \doi{10.1101/2022.02.24.481684}
#'
#' First learn a atom dictionary representation to reconstruct each cell.
#' Then, using this dictionary representation,
#' reconstruct the embeddings of each cell from the integrated atoms.
#'
#' @param object A Seurat object with all cells for one dataset
#' @param sketched.assay Assay name for sketched-cell expression (default is 'sketch')
#' @param assay Assay name for original expression (default is 'RNA')
#' @param features Features used for atomic sketch integration
#' @param reduction Dimensional reduction name for batch-corrected embeddings
#' in the sketched object (default is 'integrated_dr')
#' @param method Methods to construct sketch-cell representation
#' for all cells (default is 'sketch'). Can be one of:
#' \itemize{
#'  \item \dQuote{\code{sketch}}: Use random sketched data slot
#'  \item \dQuote{\code{data}}: Use data slot
#' }
#' @param ratio Sketch ratio of data slot when \code{dictionary.method} is set
#' to \dQuote{\code{sketch}}; defaults to 0.8
#' @param reduction.name Name to save new reduction as; defaults to
#' \code{paste0(reduction, '.orig')}
#' @param reduction.key Key for new dimensional reduction; defaults to creating
#' one from \code{reduction.name}
#' @param layers Names of layers for correction.
#' @param verbose Print progress and message
#'
#' @return Returns a Seurat object with an integrated dimensional reduction
#'
#' @importFrom MASS ginv
#' @importFrom Matrix t
#'
#' @export
#'
ProjectIntegration <- function(
  object,
  sketched.assay = 'sketch', # DefaultAssay(object)
  assay = 'RNA',
  reduction = 'integrated_dr', # harmony; rerun UMAP on this
  features = NULL, # VF from object[[atom.assay]]
  layers = 'data',
  reduction.name = NULL,
  reduction.key = NULL,
  method = c('sketch', 'data'),
  ratio = 0.8,
  sketched.layers = NULL,
  seed = 123,
  verbose = TRUE
) {
  
  layers <- Layers(object = object[[assay]], search = layers)
  # Check input and output dimensional reductions
  sketched.layers <- sketched.layers %||% layers
  reduction <- match.arg(arg = reduction, choices = Reductions(object = object))
  reduction.name <- reduction.name %||% paste0(reduction, '.full')
  reduction.key <- reduction.key %||% Key(object = reduction.name, quiet = TRUE)
  if (reduction.name %in% Reductions(object = object)) {
    warning(
      "'",
      reduction.name,
      "' already exists, overwriting",
      call. = FALSE,
      immediate. = TRUE
    )
  }
  # Check the method being used
  method <- method[1L]
  method <- match.arg(arg = method)
  # Check our layers
  sketched.assay <- match.arg(arg = sketched.assay, choices = Assays(object = object))
  assay <- match.arg(arg = assay, choices = Assays(object = object))
  layer.full <- layers
  layers <- layers %||% intersect(
    x = DefaultLayer(object[[sketched.assay]]),
    y = Layers(object[[assay]])
  )
  if (is.null(x = layer.full)) {
    sketched.assay.missing <- setdiff(x = layers, DefaultLayer(object = object[[sketched.assay]]))
    if (length(x = sketched.assay.missing) == length(x = layers)) {
      stop("None of the requested layers are present in the sketched.assay")
    } else if (length(x = sketched.assay.missing)) {
      warning(
        length(x = sketched.assay.missing),
        " layers missing from the sketched.assay",
        call. = FALSE,
        immediate. = TRUE
      )
      layers <- intersect(x = layers, y = DefaultLayer(object = object[[sketched.assay]]))
    }
  }
  # check layers
  layers.missing <- setdiff(layers, Layers(object = object[[assay]]))
  if (length(x = layers.missing)) {
    stop('layer ', layers.missing[1L], ' are not present in ', assay, " assay")
  }
  # check features
  features <- features %||% VariableFeatures(object = object[[sketched.assay]])
  # TODO: see if we can handle missing features with `union`
  features.atom <- Reduce(
    f = intersect,
    x = lapply(
      X = sketched.layers,
      FUN = function(lyr) {
        return(Features(x = object[[sketched.assay]], layer = lyr))
      }
    )
  )
  features <- intersect(x = features, y = features.atom)
  ncells <- c(
    0,
    sapply(
      X = layers,
      FUN = function(lyr) {
        return(length(x = Cells(x = object[[assay]], layer = lyr)))
      }
    )
  )
  if (length(sketched.layers) == 1) {
    sketched.layers <- rep(sketched.layers, length(layers))
  }
  sketch.matrix <- switch(
    EXPR = method,
    data = {
      R = as.sparse(
        x = diag(
          x = length(
            x = features)
          )
        )
      R
      },
    sketch = {
      R <- FeatureSketch(features = features,
                         ratio = ratio,
                         seed = seed
                         )
      R
    }
  )
  emb.list <- list()
  cells.list <- list()
  for (i in seq_along(along.with = layers)) {
    if (length(unique(sketched.layers)) == length(layers)) {
      cells.sketch <- Cells(x = object[[sketched.assay]], layer = sketched.layers[i])
    } else if (length(unique(sketched.layers)) == 1) {
      cells.sketch <- intersect(Cells(x = object[[sketched.assay]][[sketched.layers[[1]]]]),
                                Cells(object[[assay]][[layers[i] ]] ))
    }
    if (isTRUE(x = verbose)) {
      message(
        length(x = cells.sketch),
        ' atomic cells identified in the sketched.assay'
      )
      message("Correcting embeddings")
    }
    emb <- UnSketchEmbeddings(
      atom.data = LayerData(
      object = object[[sketched.assay]],
      layer = layers[i],
      features = features
    ),
    atom.cells = cells.sketch,
    orig.data = LayerData(
      object = object[[assay]],
      layer = layers[i],
      features = features
    ),
    embeddings = Embeddings(object = object[[reduction]]),
    sketch.matrix = sketch.matrix
    )
    emb.list[[i]] <- emb
    cells.list[[i]] <- colnames(x = emb)
  }
   emb.all <- t(matrix(data = unlist(emb.list),
                     nrow = ncol(x = object[[reduction]]),
                     ncol = length(unlist(cells.list))
                     ))
   rownames(emb.all) <- unlist(cells.list)
   emb.all <- emb.all[colnames(object[[assay]]), ]
  object[[reduction.name]] <- CreateDimReducObject(
    embeddings = emb.all,
    loadings = Loadings(object = object[[reduction]]),
    key = reduction.key,
    assay = assay
  )
  CheckGC()
  return(object)
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
#' @concept integration
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

#' Map query cells to a reference
#'
#' This is a convenience wrapper function around the following three functions
#' that are often run together when mapping query data to a reference:
#' \code{\link{TransferData}}, \code{\link{IntegrateEmbeddings}},
#' \code{\link{ProjectUMAP}}. Note that by default, the \code{weight.reduction}
#' parameter for all functions will be set to the dimension reduction method
#' used in the \code{\link{FindTransferAnchors}} function call used to construct
#' the anchor object, and the \code{dims} parameter will be the same dimensions
#' used to find anchors.
#'
#' @inheritParams IntegrateEmbeddings
#' @inheritParams TransferData
#' @inheritParams ProjectUMAP
#' @param store.weights Determine if the weight and anchor matrices are stored.
#' @param transferdata.args A named list of additional arguments to
#' \code{\link{TransferData}}
#' @param integrateembeddings.args A named list of additional arguments to
#' \code{\link{IntegrateEmbeddings}}
#' @param projectumap.args A named list of additional arguments to
#' \code{\link{ProjectUMAP}}
#'
#' @return Returns a modified query Seurat object containing:#'
#' \itemize{
#'   \item{New Assays corresponding to the features transferred and/or their
#'   corresponding prediction scores from \code{\link{TransferData}}}
#'   \item{An integrated reduction from \code{\link{IntegrateEmbeddings}}}
#'   \item{A projected UMAP reduction of the query cells projected into the
#'   reference UMAP using \code{\link{ProjectUMAP}}}
#' }
#'
#' @importFrom rlang invoke
#'
#' @export
#' @concept integration
#'
MapQuery <- function(
  anchorset,
  query,
  reference,
  refdata = NULL,
  new.reduction.name = NULL,
  reference.reduction = NULL,
  reference.dims = NULL,
  query.dims = NULL,
  store.weights = FALSE,
  reduction.model = NULL,
  transferdata.args = list(),
  integrateembeddings.args = list(),
  projectumap.args = list(),
  verbose = TRUE
) {
  transfer.reduction <- slot(object = anchorset, name = "command")$reduction
  if (DefaultAssay(anchorset@object.list[[1]]) %in% Assays(reference)) {
    DefaultAssay(reference) <- DefaultAssay(anchorset@object.list[[1]])
  } else {
    stop('The assay used to create the anchorset does not match any', 
         'of the assays in the reference object.')
  }
  # determine anchor type
  if (grepl(pattern = "pca", x = transfer.reduction)) {
    anchor.reduction <- "pcaproject"
    # check if the anchorset can be used for mapping
    if (is.null(x = slot(object = anchorset, name = "command")$reference.reduction)) {
      stop('The reference.reduction parameter was not set when running ',
      'FindTransferAnchors, so the resulting AnchorSet object cannot be used ',
      'in the MapQuery function.')
    }
  } else if (grepl(pattern = "cca", x = transfer.reduction)) {
    anchor.reduction <- "cca"
    ref.cca.embedding <- Embeddings(
      slot(object = anchorset, name = "object.list")[[1]][["cca"]]
      )[slot(object = anchorset, name = "reference.cells"), ]
    rownames(x = ref.cca.embedding) <- gsub(
      pattern =  "_reference",
      replacement = "",
      x = rownames(x = ref.cca.embedding)
      )
    query.cca.embedding <- Embeddings(
      slot(object = anchorset, name = "object.list")[[1]][["cca"]]
      )[slot(object = anchorset, name = "query.cells"), ]
    rownames(x = query.cca.embedding) <- gsub(
      pattern = "_query",
      replacement = "",
      x = rownames(x = query.cca.embedding)
      )
    reference[["cca"]] <- CreateDimReducObject(
      embeddings = ref.cca.embedding,
      key = "CCA_",
      assay = DefaultAssay(reference)
      )
    query[["cca"]] <- CreateDimReducObject(
      embeddings = query.cca.embedding,
      key = "CCA_",
      assay = DefaultAssay(query)
      )
    reference.reduction <- new.reduction.name <- "cca"
    reference.dims <- query.dims <- 1:ncol(x = ref.cca.embedding)
  } else if (grepl(pattern = "lsi", x = transfer.reduction)) {
    anchor.reduction <- "lsiproject"
  }  else if (grepl(pattern = "direct", x = transfer.reduction)) {
    anchor.reduction <- paste0(
      slot(object = anchorset,
           name = "command")$bridge.assay.name,
      ".reduc"
      )
    ref.reduction.emb <- Embeddings(
      object =
        slot(
          object = anchorset,
          name = "object.list"
          )[[1]][[anchor.reduction]])[
            slot(object = anchorset, name = "reference.cells"),]
    rownames(ref.reduction.emb) <- gsub(
      pattern = "_reference",
      replacement = "",
      x = rownames(ref.reduction.emb)
      )
     reference[[anchor.reduction]] <- CreateDimReducObject(
       embeddings = ref.reduction.emb,
       key = "L_",
       assay = DefaultAssay(reference)
       )
  }
  else {
    stop("unkown type of anchors")
  }
  reference.reduction <- reference.reduction %||%
    slot(object = anchorset, name = "command")$reference.reduction %||%
    anchor.reduction
  new.reduction.name <- new.reduction.name %||%
    paste0("ref.", reference.reduction)
  # checking TransferData parameters
  td.badargs <- names(x = transferdata.args)[!names(x = transferdata.args) %in% names(x = formals(fun = TransferData))]
  if (length(x = td.badargs) > 0) {
    warning("The following arguments in transferdata.args are not valid: ",
            paste(td.badargs, collapse = ", "), immediate. = TRUE, call. = FALSE)
  }
  transferdata.args <- transferdata.args[names(x = transferdata.args) %in% names(x = formals(fun = TransferData))]
  transferdata.args$weight.reduction <- transferdata.args$weight.reduction %||% anchor.reduction
  # checking IntegrateEmbeddings parameters
  ie.badargs <- names(x = integrateembeddings.args)[!names(x = integrateembeddings.args) %in% names(x = formals(fun = IntegrateEmbeddings.TransferAnchorSet))]
  if (length(x = ie.badargs) > 0) {
    warning("The following arguments in integrateembeddings.args are not valid: ",
            paste(ie.badargs, collapse = ", "), immediate. = TRUE, call. = FALSE)
  }
  integrateembeddings.args <- integrateembeddings.args[names(x = integrateembeddings.args) %in% names(x = formals(fun = IntegrateEmbeddings.TransferAnchorSet))]
  integrateembeddings.args$reductions <- integrateembeddings.args$reductions %||% anchor.reduction
  integrateembeddings.args$weight.reduction <- integrateembeddings.args$weight.reduction %||% anchor.reduction
  slot(object = query, name = "tools")$TransferData <- NULL
  reuse.weights.matrix <- FALSE
  query <- invoke(
    .fn = TransferData,
    .args  = c(list(
      anchorset = anchorset,
      reference = reference,
      query = query,
      refdata = refdata,
      store.weights = TRUE,
      only.weights =  is.null(x = refdata),
      verbose = verbose
    ), transferdata.args
    )
  )
  if (inherits(x = transferdata.args$weight.reduction , "character") &&
      transferdata.args$weight.reduction == integrateembeddings.args$weight.reduction) {
    reuse.weights.matrix <- TRUE
  }
  if (anchor.reduction != "cca") {
    query <- invoke(
      .fn = IntegrateEmbeddings,
      .args  = c(list(
        anchorset = anchorset,
        reference = reference,
        query = query,
        new.reduction.name = new.reduction.name,
        reuse.weights.matrix = reuse.weights.matrix,
        verbose = verbose
      ), integrateembeddings.args
      )
    )
    Misc(
      object = query[[new.reduction.name]],
      slot = 'ref.dims'
      ) <-  slot(object = anchorset, name = "command")$dims
  }
  slot(object = query, name = "tools")$MapQuery <- NULL
  if (store.weights) {
    slot(object = query, name = "tools")$MapQuery <- slot(
      object = query,
      name = "tools"
      )$TransferData
    slot(object = query, name = "tools")$MapQuery$anchor <- slot(
      object = anchorset,
      name = "anchors"
      )
  }
  slot(object = query, name = "tools")$TransferData <- NULL
  if (!is.null(x = reduction.model)) {
    reference.dims <- reference.dims %||% slot(object = anchorset, name = "command")$dims
    query.dims <- query.dims %||% 1:ncol(x = query[[new.reduction.name]])
    if (length(x = query.dims) != length(x = reference.dims)) {
      message("Query and reference dimensions are not equal, proceeding with reference dimensions.")
      query.dims <- reference.dims
    }
    ref_nn.num <- Misc(object = reference[[reduction.model]], slot = "model")$n_neighbors
    query <- invoke(
      .fn = ProjectUMAP,
      .args  = c(list(
        query = query,
        query.reduction = new.reduction.name,
        query.dims = query.dims,
        reference = reference,
        reference.dims = reference.dims,
        reference.reduction = reference.reduction,
        reduction.model = reduction.model,
        k.param = ref_nn.num
        ), projectumap.args
      )
    )
  }
  return(query)
}

#' @param anchors AnchorSet object or just anchor matrix from the
#' Anchorset object returned from FindTransferAnchors
#' @param combined.object  Combined object (ref + query) from the
#' Anchorset object returned
#' @param query.neighbors Neighbors object computed on query cells
#' @param ref.embeddings Reference embeddings matrix
#' @param query.embeddings Query embeddings matrix
#' @param kanchors Number of anchors to use in projection steps when computing
#' weights
#' @param ndim Number of dimensions to use when working with low dimensional
#' projections of the data
#' @param ksmooth Number of cells to average over when computing transition
#' probabilities
#' @param ksnn Number of cells to average over when determining the kernel
#' bandwidth from the SNN graph
#' @param snn.prune Amount of pruning to apply to edges in SNN graph
#' @param subtract.first.nn Option to the scoring function when computing
#' distances to subtract the distance to the first nearest neighbor
#' @param nn.method Nearest neighbor method to use (annoy or RANN)
#' @param n.trees More trees gives higher precision when using annoy approximate
#' nearest neighbor search
#' @param query.weights Query weights matrix for reuse
#' @param verbose Display messages/progress
#' @param ... Reserved for internal use
#'
#' @return Returns a vector of cell scores
#'
#' @importClassesFrom SeuratObject Neighbor
#'
#' @rdname MappingScore
#' @concept integration
#' @export
#'
MappingScore.default <- function(
  anchors,
  combined.object,
  query.neighbors,
  ref.embeddings,
  query.embeddings,
  kanchors = 50,
  ndim = 50,
  ksmooth = 100,
  ksnn = 20,
  snn.prune = 0,
  subtract.first.nn = TRUE,
  nn.method = "annoy",
  n.trees = 50,
  query.weights = NULL,
  verbose = TRUE,
  ...
) {
  CheckDots(...)
  # Input checks
  start.time <- Sys.time()
  if (is.null(x = query.neighbors) || ncol(x = query.neighbors) < ksmooth) {
    message("Recomputing query neighborhoods.\nSetting mapping.score.k in ",
            "FindTransferAnchors to the ksmooth \nvalue here (",
            ksmooth, "), can bypass this calculation in future runs.")
    query.neighbors <- FindNeighbors(
      object = query.embeddings,
      k.param = ksmooth,
      nn.method = nn.method,
      n.trees = n.trees,
      cache.index = TRUE,
      return.neighbor = TRUE,
      verbose = FALSE
    )
  }
  ref.cells <- rownames(x = ref.embeddings)
  query.cells <- rownames(query.embeddings)
  # Project reference values onto query
  if (verbose) {
    message("Projecting reference PCA onto query")
  }
  ## Need to set up an IntegrationData object to use FindWeights here
  int.mat <- matrix(data = NA, nrow = nrow(x = anchors), ncol = 0)
  rownames(x = int.mat) <- query.cells[anchors[, "cell2"]]
  slot(object = combined.object, name = 'tools')[["IT1"]] <- new(
    Class = "IntegrationData",
    anchors = anchors,
    neighbors = list(cells1 = ref.cells, cells2 = query.cells),
    integration.matrix = int.mat
  )
  ## Finding weights of anchors in query pca space
  ref.pca.orig <- ref.embeddings[, 1:ndim]
  query.pca.orig <- query.embeddings[, 1:ndim]
  dr.weights <- suppressWarnings(expr = CreateDimReducObject(
    embeddings = rbind(query.pca.orig, ref.pca.orig)
  ))
  if (!is.null(x = query.weights)) {
    weights.matrix <- query.weights
  } else {
    combined.object <- FindWeights(
      object = combined.object,
      integration.name = "IT1",
      reduction = dr.weights,
      dims = 1:ncol(x = dr.weights),
      k = kanchors,
      sd.weight = 1,
      eps = 0,
      nn.method = nn.method,
      n.trees = n.trees,
      verbose = verbose
    )
    weights.matrix <- GetIntegrationData(
      object = combined.object,
      integration.name = "IT1",
      slot = "weights"
    )
  }
  ## Perform projection of ref pca values using weights matrix
  ref.pca <- ref.embeddings[ref.cells[anchors[, 1]], 1:ndim]
  rownames(x = ref.pca) <- paste0(rownames(x = ref.pca), "_reference")
  query.cells.projected <- Matrix::crossprod(
    x = as.sparse(x = ref.pca),
    y = weights.matrix
  )
  colnames(x = query.cells.projected) <- query.cells
  rownames(x = query.cells.projected) <- colnames(x = ref.pca)

  # Re-project the query cells back onto query
  if (verbose) {
    message("Projecting back the query cells into original PCA space")
  }
  ## Compute new weights
  dr.weights <- suppressWarnings(CreateDimReducObject(
    embeddings = rbind(
      t(x = as.matrix(x = query.cells.projected)),
      ref.pca.orig[ref.cells, ]
    ),
  ))
  combined.object <- FindWeights(
    object = combined.object,
    integration.name = "IT1",
    reduction = dr.weights,
    dims = 1:ndim,
    k = kanchors,
    sd.weight = 1,
    eps = 0,
    nn.method = nn.method,
    n.trees = n.trees,
    reverse = TRUE,
    verbose = verbose
  )
  weights.matrix <- GetIntegrationData(
    object = combined.object,
    integration.name = "IT1",
    slot = "weights"
  )
  ## Project back onto query
  orig.pca <- query.embeddings[query.cells[anchors[, 2]], ]
  query.cells.back.corrected <- Matrix::t(
    x = Matrix::crossprod(
      x = as.sparse(x = orig.pca),
      y = weights.matrix)[1:ndim, ]
  )
  query.cells.back.corrected <- as.matrix(x = query.cells.back.corrected)
  rownames(x = query.cells.back.corrected) <- query.cells
  query.cells.pca <- query.embeddings[query.cells, 1:ndim]
  if (verbose) {
    message("Computing scores:")
    message("    Finding neighbors of original query cells")
  }
  ## Compute original neighborhood of query cells
  if (is.null(x = query.neighbors)) {
    query.neighbors <- NNHelper(
      data = query.cells.pca,
      query = query.cells.pca,
      k = max(ksmooth, ksnn),
      method = nn.method,
      n.trees = n.trees,
      cache.index = TRUE
    )
  }
  if (verbose) {
    message("    Finding neighbors of transformed query cells")
  }
  ## Compute new neighborhood of query cells after projections
  if (nn.method == "annoy") {
    if (is.null(x = Index(object = query.neighbors))) {
      corrected.neighbors <- NNHelper(
        data = query.cells.pca,
        query = query.cells.back.corrected,
        k = max(ksmooth, ksnn),
        method = nn.method,
        n.treees = n.trees,
        cache.index = TRUE
      )
    } else {
      corrected.neighbors <- AnnoySearch(
        index = Index(object = query.neighbors),
        query = query.cells.back.corrected,
        k = max(ksmooth, ksnn)
      )
      corrected.neighbors <- new(
        Class = 'Neighbor',
        nn.idx = corrected.neighbors$nn.idx,
        nn.dist = corrected.neighbors$nn.dists
      )
    }
  }
  if (verbose) {
    message("    Computing query SNN")
  }
  snn <- ComputeSNN(
    nn_ranked = Indices(query.neighbors)[, 1:ksnn],
    prune = snn.prune
  )
  query.cells.pca <- t(x = query.cells.pca)
  if (verbose) {
    message("    Determining bandwidth and computing transition probabilities")
  }
  scores <- ScoreHelper(
    snn = snn,
    query_pca = query.cells.pca,
    query_dists = Distances(object = query.neighbors),
    corrected_nns = Indices(object = corrected.neighbors),
    k_snn = ksnn,
    subtract_first_nn = subtract.first.nn,
    display_progress = verbose
  )
  scores[scores > 1] <- 1
  names(x = scores) <- query.cells
  end.time <- Sys.time()
  if (verbose) {
    message("Total elapsed time: ", end.time - start.time)
  }
  return(scores)
}

#' @rdname MappingScore
#' @export
#' @concept integration
#' @method MappingScore AnchorSet
#'
MappingScore.AnchorSet <- function(
  anchors,
  kanchors = 50,
  ndim = 50,
  ksmooth = 100,
  ksnn = 20,
  snn.prune = 0,
  subtract.first.nn = TRUE,
  nn.method = "annoy",
  n.trees = 50,
  query.weights = NULL,
  verbose = TRUE,
  ...
) {
  CheckDots(...)
  combined.object <- slot(object = anchors, name = "object.list")[[1]]
  combined.object <- RenameCells(
    object = combined.object,
    new.names = unname(obj = sapply(
      X = Cells(x = combined.object),
      FUN = RemoveLastField
    ))
  )
  query.cells <- sapply(
    X = slot(object = anchors, name = "query.cells"),
    FUN = RemoveLastField
  )
  ref.cells <- sapply(
    X = slot(object = anchors, name = "reference.cells"),
    FUN = RemoveLastField
  )
  query.embeddings <- Embeddings(object = subset(
    x = combined.object[["pcaproject.l2"]],
    cells = query.cells
  ))
  ref.embeddings <- Embeddings(object = subset(
    x = combined.object[["pcaproject.l2"]],
    cells = ref.cells
  ))
  query.neighbors <- slot(object = anchors, name = "neighbors")[["query.neighbors"]]
  # reduce size of anchorset combined object
  combined.object <- DietSeurat(object = combined.object)
  combined.object <- subset(
    x = combined.object,
    features = c(rownames(x = combined.object)[1])
  )
  for (i in colnames(x = combined.object[[]])) {
    combined.object[[i]] <- NULL
  }
  return(MappingScore(
    anchors = slot(object = anchors, name = "anchors"),
    combined.object = combined.object,
    query.neighbors = query.neighbors,
    ref.embeddings = ref.embeddings,
    query.embeddings = query.embeddings,
    kanchors = kanchors,
    ndim = ndim,
    ksmooth = ksmooth,
    ksnn = ksnn,
    snn.prune = snn.prune,
    subtract.first.nn = subtract.first.nn,
    nn.method = nn.method,
    n.trees = n.trees,
    query.weights = query.weights,
    verbose = verbose
  ))
}

#' Calculates a mixing metric
#'
#' Here we compute a measure of how well mixed a composite dataset is. To
#' compute, we first examine the local neighborhood for each cell (looking at
#' max.k neighbors) and determine for each group (could be the dataset after
#' integration) the k nearest neighbor and what rank that neighbor was in the
#' overall neighborhood. We then take the median across all groups as the mixing
#' metric per cell.
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
#' @return Returns a vector of values of the mixing metric for each cell
#'
#' @importFrom RANN nn2
#' @importFrom pbapply pbsapply
#' @importFrom future.apply future_sapply
#' @importFrom future nbrOfWorkers
#' @export
#' @concept integration
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
    test = verbose && nbrOfWorkers() == 1,
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

#' Prepare an object list normalized with sctransform for integration.
#'
#' This function takes in a list of objects that have been normalized with the
#' \code{\link{SCTransform}} method and performs the following steps:
#' \itemize{
#'   \item{If anchor.features is a numeric value, calls \code{\link{SelectIntegrationFeatures}}
#'   to determine the features to use in the downstream integration procedure.}
#'   \item{Ensures that the sctransform residuals for the features specified
#'   to anchor.features are present in each object in the list. This is
#'   necessary because the default behavior of \code{\link{SCTransform}} is to
#'   only store the residuals for the features determined to be variable.
#'   Residuals are recomputed for missing features using the stored model
#'   parameters via the \code{\link{GetResidual}} function.}
#'   \item{Subsets the \code{scale.data} slot to only contain the residuals for
#'   anchor.features for efficiency in downstream processing. }
#' }
#'
#' @param object.list A list of \code{\link{Seurat}} objects to prepare for integration
#' @param assay The name of the \code{\link{Assay}} to use for integration. This can be a
#' single name if all the assays to be integrated have the same name, or a character vector
#' containing the name of each \code{\link{Assay}} in each object to be integrated. The
#' specified assays must have been normalized using \code{\link{SCTransform}}.
#' If NULL (default), the current default assay for each object is used.
#' @param anchor.features Can be either:
#' \itemize{
#'   \item{A numeric value. This will call \code{\link{SelectIntegrationFeatures}}
#'   to select the provided number of features to be used in anchor finding}
#'   \item{A vector of features to be used as input to the anchor finding
#'   process}
#' }
#' @param sct.clip.range Numeric of length two specifying the min and max values
#' the Pearson residual will be clipped to
#' @param verbose Display output/messages
#'
#' @return A list of \code{\link{Seurat}} objects with the appropriate \code{scale.data} slots
#' containing only the required \code{anchor.features}.
#'
#' @importFrom pbapply pblapply
#' @importFrom methods slot slot<-
#' @importFrom future nbrOfWorkers
#' @importFrom future.apply future_lapply
#'
#' @export
#' @concept integration
#' @examples
#' \dontrun{
#' # to install the SeuratData package see https://github.com/satijalab/seurat-data
#' library(SeuratData)
#' data("panc8")
#'
#' # panc8 is a merged Seurat object containing 8 separate pancreas datasets
#' # split the object by dataset and take the first 2 to integrate
#' pancreas.list <- SplitObject(panc8, split.by = "tech")[1:2]
#'
#' # perform SCTransform normalization
#' pancreas.list <- lapply(X = pancreas.list, FUN = SCTransform)
#'
#' # select integration features and prep step
#' features <- SelectIntegrationFeatures(pancreas.list)
#' pancreas.list <- PrepSCTIntegration(
#'   pancreas.list,
#'   anchor.features = features
#' )
#'
#' # downstream integration steps
#' anchors <- FindIntegrationAnchors(
#'   pancreas.list,
#'   normalization.method = "SCT",
#'   anchor.features = features
#' )
#' pancreas.integrated <- IntegrateData(anchors, normalization.method = "SCT")
#' }
#'
PrepSCTIntegration <- function(
  object.list,
  assay = NULL,
  anchor.features = 2000,
  sct.clip.range = NULL,
  verbose = TRUE
) {
  my.lapply <- ifelse(
    test = verbose && nbrOfWorkers() == 1,
    yes = pblapply,
    no = future_lapply
  )
  assay <- assay %||% sapply(X = object.list, FUN = DefaultAssay)
  assay <- rep_len(x = assay, length.out = length(x = object.list))
  objects.names <- names(x = object.list)
  object.list <- lapply(
    X = 1:length(x = object.list),
    FUN = function(i) {
      DefaultAssay(object = object.list[[i]]) <- assay[i]
      object.list[[i]][[assay[i]]] <- as(object = object.list[[i]][[assay[i]]], Class = "SCTAssay")
      return(object.list[[i]])
    }
  )
  sct.check <- vapply(
    X = 1:length(x = object.list),
    FUN = function(i) {
      sct.check <- IsSCT(assay = object.list[[i]][[assay[i]]])
      if (!sct.check) {
        if ("FindIntegrationAnchors" %in% Command(object = object.list[[i]]) &&
            Command(object = object.list[[i]], command = "FindIntegrationAnchors", value = "normalization.method") == "SCT") {
          sct.check <- TRUE
        }
      }
      return(sct.check)
    },
    FUN.VALUE = logical(length = 1L),
    USE.NAMES = FALSE
  )
  if (!all(sct.check)) {
    stop(
      "The following assays have not been processed with SCTransform:\n",
      paste(
        ' object:',
        which(x = !sct.check, useNames = FALSE),
        '- assay:',
        assay[!sct.check],
        collapse = '\n'
      ),
      call. = FALSE
    )
  }
  if (is.numeric(x = anchor.features)) {
    anchor.features <- SelectIntegrationFeatures(
      object.list = object.list,
      nfeatures = anchor.features,
      verbose = verbose
    )
  }
  object.list <- my.lapply(
    X = 1:length(x = object.list),
    FUN = function(i) {
      obj <- GetResidual(
        object = object.list[[i]],
        assay = assay[i],
        features = anchor.features,
        replace.value = ifelse(test = is.null(x = sct.clip.range), yes = FALSE, no = TRUE),
        clip.range = sct.clip.range,
        verbose = FALSE
      )
      scale.data <- GetAssayData(
        object = obj,
        assay = assay[i],
        slot = 'scale.data'
      )
      obj <- SetAssayData(
        object = obj,
        slot = 'scale.data',
        new.data = scale.data[anchor.features, ],
        assay = assay[i]
      )
      return(obj)
    }
  )
  assays.used <- assay
  for (i in 1:length(x = object.list)) {
    assay <- as.character(x = assays.used[i])
    object.list[[i]] <- LogSeuratCommand(object = object.list[[i]])
  }
  names(x = object.list) <- objects.names
  return(object.list)
}

#' Select integration features
#'
#' Choose the features to use when integrating multiple datasets. This function
#' ranks features by the number of datasets they are deemed variable in,
#' breaking ties by the median variable feature rank across datasets. It returns
#' the top scoring features by this ranking.
#'
#' If for any assay in the list, \code{\link{FindVariableFeatures}} hasn't been
#' run, this method will try to run it using the \code{fvf.nfeatures} parameter
#' and any additional ones specified through the \dots.
#'
#' @param object.list List of seurat objects
#' @param nfeatures Number of features to return
#' @param assay Name or vector of assay names (one for each object) from which
#' to pull the variable features.
#' @param verbose Print messages
#' @param fvf.nfeatures nfeatures for \code{\link{FindVariableFeatures}}. Used
#' if \code{VariableFeatures} have not been set for any object in
#' \code{object.list}.
#' @param ... Additional parameters to \code{\link{FindVariableFeatures}}
#'
#' @return A vector of selected features
#'
#' @importFrom utils head
#'
#' @export
#' @concept integration
#'
#' @examples
#' \dontrun{
#' # to install the SeuratData package see https://github.com/satijalab/seurat-data
#' library(SeuratData)
#' data("panc8")
#'
#' # panc8 is a merged Seurat object containing 8 separate pancreas datasets
#' # split the object by dataset and take the first 2
#' pancreas.list <- SplitObject(panc8, split.by = "tech")[1:2]
#'
#' # perform SCTransform normalization
#' pancreas.list <- lapply(X = pancreas.list, FUN = SCTransform)
#'
#' # select integration features
#' features <- SelectIntegrationFeatures(pancreas.list)
#' }
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
    for (ii in length(x = object.list)) {
      DefaultAssay(object = object.list[[ii]]) <- assay[ii]
    }
  } else {
    assay <- sapply(X = object.list, FUN = DefaultAssay)
  }
  for (ii in 1:length(x = object.list)) {
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
  for (i in 1:length(x = object.list)) {
    var.features <- var.features[names(x = var.features) %in% rownames(x = object.list[[i]][[assay[i]]])]
  }
  tie.val <- var.features[min(nfeatures, length(x = var.features))]
  features <- names(x = var.features[which(x = var.features > tie.val)])
  vf.list <- lapply(X = object.list, FUN = VariableFeatures)
  if (length(x = features) > 0) {
    feature.ranks <- sapply(X = features, FUN = function(x) {
      ranks <- sapply(X = vf.list, FUN = function(vf) {
        if (x %in% vf) {
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
    ranks <- sapply(X = vf.list, FUN = function(vf) {
      if (x %in% vf) {
        return(which(x = x == vf))
      }
      return(NULL)
    })
    median(x = unlist(x = ranks))
  })
  features <- c(
    features,
    names(x = head(x = sort(x = tie.ranks), nfeatures - length(x = features)))
  )
  return(features)
}

.FeatureRank <- function(features, flist, ranks = FALSE) {
  franks <- vapply(
    X = features,
    FUN = function(x) {
      return(median(x = unlist(x = lapply(
        X = flist,
        FUN = function(fl) {
          if (x %in% fl) {
            return(which(x = x == fl))
          }
          return(NULL)
        }
      ))))
    },
    FUN.VALUE = numeric(length = 1L)
  )
  franks <- sort(x = franks)
  if (!isTRUE(x = ranks)) {
    franks <- names(x = franks)
  }
  return(franks)
}

#' @export
#'
SelectIntegrationFeatures5 <- function(
  object,
  nfeatures = 2000,
  assay = NULL,
  method = NULL,
  layers = NULL,
  verbose = TRUE,
  ...
) {
  assay <- assay %||% DefaultAssay(object = object)
  layers <- Layers(object = object[[assay]], search = layers)
  var.features <- VariableFeatures(
    object = object,
    assay = assay,
    nfeatures = nfeatures,
    method = method,
    layer = layers,
    simplify = TRUE
  )
  return(var.features)
}

#' @export
#'
SelectSCTIntegrationFeatures <- function(
  object,
  nfeatures = 3000,
  assay = NULL,
  verbose = TRUE,
  ...
) {
  assay <- assay %||% DefaultAssay(object = object)
  if (!inherits(x = object[[assay]], what = 'SCTAssay')) {
    abort(message = "'assay' must be an SCTAssay")
  }
  models <- levels(x = object[[assay]])
  vf.list <- VariableFeatures(
    object = object[[assay]],
    layer = models,
    nfeatures = nfeatures,
    simplify = FALSE
  )
  var.features <- sort(
    x = table(unlist(x = vf.list, use.names = FALSE)),
    decreasing = TRUE
  )
  for (i in 1:length(x = models)) {
    vst_out <- SCTModel_to_vst(SCTModel = slot(object = object[[assay]], name = "SCTModel.list")[[models[[i]]]])
    var.features <- var.features[names(x = var.features) %in% rownames(x = vst_out$gene_attr)]
  }
  tie.val <- var.features[min(nfeatures, length(x = var.features))]
  features <- names(x = var.features[which(x = var.features > tie.val)])
  if (length(x = features) > 0) {
    feature.ranks <- sapply(X = features, FUN = function(x) {
      ranks <- sapply(X = vf.list, FUN = function(vf) {
        if (x %in% vf) {
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
    ranks <- sapply(X = vf.list, FUN = function(vf) {
      if (x %in% vf) {
        return(which(x = x == vf))
      }
      return(NULL)
    })
    median(x = unlist(x = ranks))
  })
  features <- c(
    features,
    names(x = head(x = sort(x = tie.ranks), nfeatures - length(x = features)))
  )
  return(features)
}

#' Transfer data
#'
#' Transfer categorical or continuous data across single-cell datasets. For
#' transferring categorical information, pass a vector from the reference
#' dataset (e.g. \code{refdata = reference$celltype}). For transferring
#' continuous information, pass a matrix from the reference dataset (e.g.
#' \code{refdata = GetAssayData(reference[['RNA']])}).
#'
#' The main steps of this procedure are outlined below. For a more detailed
#' description of the methodology, please see Stuart, Butler, et al Cell 2019.
#' \doi{10.1016/j.cell.2019.05.031}; \doi{10.1101/460147}
#'
#' For both transferring discrete labels and also feature imputation, we first
#' compute the weights matrix.
#'
#' \itemize{
#'   \item{Construct a weights matrix that defines the association between each
#'   query cell and each anchor. These weights are computed as 1 - the distance
#'   between the query cell and the anchor divided by the distance of the query
#'   cell to the \code{k.weight}th anchor multiplied by the anchor score
#'   computed in \code{\link{FindIntegrationAnchors}}. We then apply a Gaussian
#'   kernel width a bandwidth defined by \code{sd.weight} and normalize across
#'   all \code{k.weight} anchors.}
#' }
#'
#' The main difference between label transfer (classification) and feature
#' imputation is what gets multiplied by the weights matrix. For label transfer,
#' we perform the following steps:
#'
#' \itemize{
#'   \item{Create a binary classification matrix, the rows corresponding to each
#'   possible class and the columns corresponding to the anchors. If the
#'   reference cell in the anchor pair is a member of a certain class, that
#'   matrix entry is filled with a 1, otherwise 0.}
#'   \item{Multiply this classification matrix by the transpose of weights
#'   matrix to compute a prediction score for each class for each cell in the
#'   query dataset.}
#' }
#'
#' For feature imputation, we perform the following step:
#' \itemize{
#'   \item{Multiply the expression matrix for the reference anchor cells by the
#'   weights matrix. This returns a predicted expression matrix for the
#'   specified features for each cell in the query dataset.}
#' }
#'
#'
#' @param anchorset An \code{\link{AnchorSet}} object generated by
#' \code{\link{FindTransferAnchors}}
#' @param refdata Data to transfer. This can be specified in one of two ways:
#' \itemize{
#'   \item{The reference data itself as either a vector where the names
#'   correspond to the reference cells, or a matrix, where the column names
#'   correspond to the reference cells.}
#'   \item{The name of the metadata field or assay from the reference object
#'   provided. This requires the reference parameter to be specified. If pulling
#'   assay data in this manner, it will pull the data from the data slot. To
#'   transfer data from other slots, please pull the data explicitly with
#'   \code{\link{GetAssayData}} and provide that matrix here.}
#' }
#' @param reference Reference object from which to pull data to transfer
#' @param query Query object into which the data will be transferred.
#' @param weight.reduction Dimensional reduction to use for the weighting
#' anchors. Options are:
#' \itemize{
#'    \item{pcaproject: Use the projected PCA used for anchor building}
#'    \item{lsiproject: Use the projected LSI used for anchor building}
#'    \item{pca: Use an internal PCA on the query only}
#'    \item{cca: Use the CCA used for anchor building}
#'    \item{custom DimReduc: User provided \code{\link{DimReduc}} object
#'    computed on the query cells}
#' }
#' @param l2.norm Perform L2 normalization on the cell embeddings after
#' dimensional reduction
#' @param dims Set of dimensions to use in the anchor weighting procedure. If
#' NULL, the same dimensions that were used to find anchors will be used for
#' weighting.
#' @param k.weight Number of neighbors to consider when weighting anchors
#' @param sd.weight Controls the bandwidth of the Gaussian kernel for weighting
#' @param eps Error bound on the neighbor finding algorithm (from
#' \code{\link{RANN}})
#' @param n.trees More trees gives higher precision when using annoy approximate
#' nearest neighbor search
#' @param verbose Print progress bars and output
#' @param slot Slot to store the imputed data. Must be either "data" (default)
#' or "counts"
#' @param prediction.assay Return an \code{Assay} object with the prediction
#' scores for each class stored in the \code{data} slot.
#' @param only.weights Only return weights matrix
#' @param store.weights Optionally store the weights matrix used for predictions
#' in the returned query object.
#'
#' @return
#' If \code{query} is not provided, for the categorical data in \code{refdata},
#' returns a data.frame with label predictions. If \code{refdata} is a matrix,
#' returns an Assay object where the imputed data has been stored in the
#' provided slot.
#'
#' If \code{query} is provided, a modified query object is returned. For
#' the categorical data in refdata, prediction scores are stored as Assays
#' (prediction.score.NAME) and two additional metadata fields: predicted.NAME
#' and predicted.NAME.score which contain the class prediction and the score for
#' that predicted class. For continuous data, an Assay called NAME is returned.
#' NAME here corresponds to the name of the element in the refdata list.
#'
#'
#' @references Stuart T, Butler A, et al. Comprehensive Integration of
#' Single-Cell Data. Cell. 2019;177:1888-1902 \doi{10.1016/j.cell.2019.05.031}
#'
#' @importFrom Matrix t
#'
#' @export
#'
#' @concept integration
#'
#' @examples
#' \dontrun{
#' # to install the SeuratData package see https://github.com/satijalab/seurat-data
#' library(SeuratData)
#' data("pbmc3k")
#'
#' # for demonstration, split the object into reference and query
#' pbmc.reference <- pbmc3k[, 1:1350]
#' pbmc.query <- pbmc3k[, 1351:2700]
#'
#' # perform standard preprocessing on each object
#' pbmc.reference <- NormalizeData(pbmc.reference)
#' pbmc.reference <- FindVariableFeatures(pbmc.reference)
#' pbmc.reference <- ScaleData(pbmc.reference)
#'
#' pbmc.query <- NormalizeData(pbmc.query)
#' pbmc.query <- FindVariableFeatures(pbmc.query)
#' pbmc.query <- ScaleData(pbmc.query)
#'
#' # find anchors
#' anchors <- FindTransferAnchors(reference = pbmc.reference, query = pbmc.query)
#'
#' # transfer labels
#' predictions <- TransferData(anchorset = anchors, refdata = pbmc.reference$seurat_annotations)
#' pbmc.query <- AddMetaData(object = pbmc.query, metadata = predictions)
#' }
#'
TransferData <- function(
  anchorset,
  refdata,
  reference = NULL,
  query = NULL,
  weight.reduction = 'pcaproject',
  l2.norm = FALSE,
  dims = NULL,
  k.weight = 50,
  sd.weight = 1,
  eps = 0,
  n.trees = 50,
  verbose = TRUE,
  slot = "data",
  prediction.assay = FALSE,
  only.weights = FALSE,
  store.weights = TRUE
) {
  combined.ob <- slot(object = anchorset, name = "object.list")[[1]]
  anchors <- slot(object = anchorset, name = "anchors")
  reference.cells <- slot(object = anchorset, name = "reference.cells")
  query.cells <- slot(object = anchorset, name = "query.cells")
  label.transfer <- list()
  ValidateParams_TransferData(
    anchorset = anchorset,
    combined.ob = combined.ob,
    anchors = anchors,
    reference.cells = reference.cells,
    query.cells = query.cells,
    refdata = refdata,
    reference = reference,
    query = query,
    weight.reduction = weight.reduction,
    l2.norm = l2.norm,
    dims = dims,
    k.weight = k.weight,
    sd.weight = sd.weight,
    eps = eps,
    n.trees = n.trees,
    verbose = verbose,
    only.weights = only.weights,
    slot = slot,
    prediction.assay = prediction.assay,
    label.transfer = label.transfer
  )
  if (!inherits(x = weight.reduction, what = "DimReduc") && weight.reduction == 'pca') {
    if (verbose) {
      message("Running PCA on query dataset")
    }

    features <- slot(object = anchorset, name = "anchor.features")
    query.ob <- query
    query.ob <- ScaleData(object = query.ob, features = features, verbose = FALSE)
    query.ob <- RunPCA(object = query.ob, npcs = max(dims), features = features, verbose = FALSE)
    query.pca <- Embeddings(query.ob[['pca']])
    rownames(x = query.pca) <- paste0(rownames(x = query.pca), "_query")
    #fill with 0s
    ref.pca <- matrix(
      data = 0,
      nrow = length(x = reference.cells),
      ncol = ncol(x = query.pca),
      dimnames = list(reference.cells, colnames(x = query.pca))
    )
    rm(query.ob)
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
  if (!inherits(x = weight.reduction, what = "DimReduc") && weight.reduction == "lsi") {
    if (!("lsi" %in% Reductions(object = query))) {
      stop("Requested lsi for weight.reduction, but lsi not stored in query object.")
    } else {
      weight.reduction <- query[["lsi"]]
    }
  }
  if (inherits(x = weight.reduction, what = "DimReduc")) {
    weight.reduction <- RenameCells(
      object = weight.reduction,
      new.names = paste0(Cells(x = weight.reduction), "_query")
    )
  } else {
    if (l2.norm) {
      weight.reduction.l2 <- paste0(weight.reduction, ".l2")
      if (weight.reduction.l2 %in% Reductions(object = combined.ob)) {
        combined.ob <- L2Dim(object = combined.ob, reduction = weight.reduction)
      }
      weight.reduction <- weight.reduction.l2
    }
    weight.reduction <- combined.ob[[weight.reduction]]
  }
  dims <- dims %||% seq_len(length.out = ncol(x = weight.reduction))
  if (max(dims) > ncol(x = weight.reduction)) {
    stop("dims is larger than the number of available dimensions in ",
         "weight.reduction (", ncol(x = weight.reduction), ").", call. = FALSE)
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
    reduction = weight.reduction,
    dims = dims,
    k = k.weight,
    sd.weight = sd.weight,
    eps = eps,
    n.trees = n.trees,
    verbose = verbose
  )
  weights <- GetIntegrationData(
    object = combined.ob,
    integration.name = "integrated",
    slot = 'weights'
  )
  if (only.weights) {
    if (is.null(x = query)) {
      return(weights)
    } else {
      slot(object = query, name = "tools")[["TransferData"]] <- list(weights.matrix = weights)
      return(query)
    }
  }
  anchors <- as.data.frame(x = anchors)
  query.cells <- unname(obj = sapply(
    X = query.cells,
    FUN = function(x) gsub(pattern = "_query", replacement = "", x = x)
  ))
  transfer.results <- list()
  for (rd in 1:length(x = refdata)) {
    if (isFALSE(x = refdata[[rd]])) {
      transfer.results[[rd]] <- NULL
      next
    }
    rd.name <- names(x = refdata)[rd]
    # case for projection
    if (label.transfer[[rd]]) {
      anchors$id1 <- refdata[[rd]][anchors[, "cell1"]]
      reference.ids <- factor(x = anchors$id1, levels = unique(x = refdata[[rd]]))
      possible.ids <- levels(x = reference.ids)
      prediction.mat <- matrix(
        nrow = nrow(x = anchors),
        ncol = length(x = possible.ids),
        data = 0
      )
      for (i in 1:length(x = possible.ids)) {
        prediction.mat[which(reference.ids == possible.ids[i]), i] = 1
      }
      if (verbose) {
        message("Predicting cell labels")
      }
      prediction.scores <- t(x = weights) %*% prediction.mat
      colnames(x = prediction.scores) <- possible.ids
      rownames(x = prediction.scores) <- query.cells
      if ("bridge.sets" %in% names(anchorset@weight.reduction@misc)) {
        bridge.weight <- anchorset@weight.reduction@misc$bridge.sets
        bridge.prediction.matrix <- as.sparse(
          x = dummy_cols(
            refdata[[rd]][ bridge.weight$bridge.ref_anchor ]
          )[, -1]
        )
        colnames(bridge.prediction.matrix) <- gsub(
          pattern = ".data_",
          replacement = "",
          x = colnames(bridge.prediction.matrix)
        )
        extra.id <- setdiff(possible.ids, colnames(bridge.prediction.matrix))
        if (length(extra.id) > 0) {
          extra.prediction <- as.sparse(x = matrix(
            data = 0,
            nrow = nrow(bridge.prediction.matrix),
            ncol = length(extra.id)
          ))
          colnames(extra.prediction) <- extra.id
          bridge.prediction.matrix <- cbind(
            bridge.prediction.matrix,
            extra.prediction
            )
        }
        bridge.prediction.matrix <- bridge.prediction.matrix[,possible.ids, drop = FALSE]
         bridge.prediction.scores <- t(bridge.weight$query.weights) %*%
        (t(bridge.weight$bridge.weights) %*%
           bridge.prediction.matrix)[bridge.weight$query.ref_anchor,]
         prediction.scores <- (prediction.scores + bridge.prediction.scores)/2
         prediction.scores <- as.matrix(x = prediction.scores)
      }
      prediction.ids <- possible.ids[apply(X = prediction.scores, MARGIN = 1, FUN = which.max)]
      prediction.ids <- as.character(prediction.ids)
      prediction.max <- apply(X = prediction.scores, MARGIN = 1, FUN = max)
      if (is.null(x = query)) {
        prediction.scores <- cbind(prediction.scores, max = prediction.max)
      }
      predictions <- data.frame(
        predicted.id = prediction.ids,
        prediction.score = as.matrix(prediction.scores),
        row.names = query.cells,
        stringsAsFactors = FALSE
      )
      if (prediction.assay || !is.null(x = query)) {
        predictions <- CreateAssayObject(
          data = t(x = as.matrix(x = prediction.scores)), check.matrix = FALSE
        )
        Key(object = predictions) <- Key(paste0("predictionscore", rd.name), quiet = TRUE)
      }
      if (is.null(x = query)) {
        transfer.results[[rd]] <- predictions
      } else {
        query <- AddMetaData(object = query, metadata = prediction.max, col.name = paste0("predicted.", rd.name, ".score"))
        query <- AddMetaData(object = query, metadata = prediction.ids, col.name = paste0("predicted.", rd.name))
        query[[paste0("prediction.score.", rd.name)]] <- predictions
      }
    } else {  # case for transferring features
      reference.cell.indices <- reference.cells[anchors$cell1]
      refdata.anchors <- refdata[[rd]][, reference.cell.indices]
      nfeatures <- nrow(x = refdata[[rd]])
      if (verbose) {
        message(paste0("Transfering ", nfeatures, " features onto reference data"))
      }
      new.data <- refdata.anchors %*% weights
      rownames(x = new.data) <- rownames(x = refdata[[rd]])
      colnames(x = new.data) <- query.cells
      if (inherits(x = new.data, what = "Matrix")) {
        new.data <- as.sparse(x = new.data)
      }
      if (slot == "counts") {
        new.assay <- CreateAssayObject(counts = new.data, check.matrix = FALSE)
      } else if (slot == "data") {
        new.assay <- CreateAssayObject(data = new.data, check.matrix = FALSE)
      }
      Key(object = new.assay) <- Key(rd.name, quiet = TRUE)
      if (is.null(x = query)) {
        transfer.results[[rd]] <- new.assay
      } else {
        if (rd.name %in% Assays(object = query)) {
          message(
            rd.name,
            " already present in query. ",
            "Storing as ",
            paste0("predicted_", rd.name)
          )
          rd.name <- paste0("predicted_", rd.name)
        }
        query[[rd.name]] <- new.assay
      }
    }
  }
  if (is.null(x = query)) {
    names(x = transfer.results) <- names(x = refdata)
    if (length(x = transfer.results) == 1) {
      transfer.results <- transfer.results[[1]]
    }
    return(transfer.results)
  } else {
    if (store.weights) {
      slot(object = query, name = "tools")[["TransferData"]] <- list(weights.matrix = weights)
    }
    return(query)
  }
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Methods for Seurat-defined generics
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @param object.list List of Seurat objects
#' @rdname AnnotateAnchors
#' @export
#' @method AnnotateAnchors default
#' @concept integration
#'
AnnotateAnchors.default <- function(
  anchors,
  vars = NULL,
  slot = NULL,
  object.list,
  assay = NULL,
  ...
) {
  # reorder columns
  anchors <- anchors[, c("cell1", "dataset1", "cell2", "dataset2", "score")]
  colnames(x = anchors)[5] <- "anchor.score"
  cell.names <- lapply(X = object.list, FUN = Cells)
  cell1.names <- character(length = nrow(x = anchors))
  for (dataset in unique(x = anchors$dataset1)) {
    dataset.cells <- which(x = anchors$dataset1 == dataset)
    cell1.names[dataset.cells] <- cell.names[[dataset]][anchors[dataset.cells, "cell1"]]
  }
  anchors$cell1 <- cell1.names
  cell2.names <- character(length(x = nrow(x = anchors)))
  for (dataset in unique(x = anchors$dataset2)) {
    dataset.cells <- which(x = anchors$dataset2 == dataset)
    cell2.names[dataset.cells] <- cell.names[[dataset]][anchors[dataset.cells, "cell2"]]
  }
  anchors$cell2 <- cell2.names
  slot <- slot %||% "data"
  assay <- assay %||% sapply(X = object.list, FUN = DefaultAssay)
  if (length(x = assay) == 1) {
    assay <- rep(x = assay, times = length(x = object.list))
  }
  if (length(x = assay) != length(x = object.list)) {
    stop("Number of assays provided should either be one or the length of object.list")
  }
  for (ob in 1:length(x = object.list)) {
    DefaultAssay(object = object.list[[ob]]) <- assay[ob]
  }
  if (length(x = slot) == 1) {
    slot <- rep(x = slot, times = length(x = vars))
  }
  if (length(x = vars) > 0) {
    for(v in 1:length(x = vars)) {
      var <- vars[v]
      var.list <- lapply(X = object.list, FUN = function(x) {
        tryCatch(
          expr = FetchData(object = x, vars = var, slot = slot[v]),
          error = function(e) {
            data.fetched <- as.data.frame(
              x = rep(x = NA, times = ncol(x = x)),
              row.names = Cells(x = x),
              stringsAsFactors = FALSE
            )
            colnames(x = data.fetched) <- var
            return(data.fetched)
          }
        )
      })
      if (all(unlist(x = lapply(X = var.list, FUN = isFALSE)))) {
        warning(
          var, " not found in all objects",
          call. = FALSE,
          immediate. = TRUE
        )
        next
      }
      if (any(unlist(x = lapply(X = var.list, FUN = isFALSE)))) {
        warning(
          var, " not in all objects. Filling missing objects with NA",
          call. = FALSE,
          immediate. = TRUE
        )
      }
      if (is.null(x = names(x = object.list))) {
        names(x = var.list) <- 1:length(x = object.list)
      } else {
        names(x = var.list) <- names(x = object.list)
      }
      for(i in c(1, 2)) {
        cell <- paste0("cell", i)
        if (is.factor(x = anchors[, cell])) {
          anchors[, cell] <- as.character(x = anchors[, cell])
        }
        for (j in unique(x = anchors[, paste0("dataset", i)])) {
          var.df <- var.list[[j]]
          dataset.cells <- which(x = anchors[, paste0("dataset", i)] == j)
          anchors[dataset.cells, paste0(cell, ".", var)] <- var.df[anchors[, cell][dataset.cells], ]
        }
      }
      # column specifying whether the annotation matches across pair of datasets
      anchors[, paste0(var, ".match")] <- anchors[, paste0("cell1.", var)] ==
        anchors[, paste0("cell2.", var)]
    }
  }
  return(anchors)
}

#' @rdname AnnotateAnchors
#' @export
#' @method AnnotateAnchors IntegrationAnchorSet
#'
AnnotateAnchors.IntegrationAnchorSet <- function(
  anchors,
  vars = NULL,
  slot = NULL,
  object.list = NULL,
  assay = NULL,
  ...
) {
  anchor.df <- slot(object = anchors, name = 'anchors')
  object.list <- object.list %||% slot(object = anchors, name = 'object.list')
  anchor.df <- as.data.frame(x = anchor.df)
  anchor.df <- AnnotateAnchors(
    anchors = anchor.df,
    vars = vars,
    slot = slot,
    object.list = object.list,
    assay = assay
  )
  return(anchor.df)
}

#' @param reference Reference object used in \code{\link{FindTransferAnchors}}
#' @param query Query object used in \code{\link{FindTransferAnchors}}
#' @rdname AnnotateAnchors
#' @export
#' @method AnnotateAnchors TransferAnchorSet
#'
AnnotateAnchors.TransferAnchorSet <- function(
  anchors,
  vars = NULL,
  slot = NULL,
  reference = NULL,
  query = NULL,
  assay = NULL,
  ...
) {
  anchor.df <- slot(object = anchors, name = 'anchors')
  if (class(x = reference) != class(x = query)) {
    stop("If setting reference/query, please set both parameters.")
  }
  if (is.null(x = reference)) {
    object.list <- slot(object = anchors, name = 'object.list')[[1]]
    reference.cells <- slot(object = anchors, name = "reference.cells")
    reference <- subset(x = object.list, cells = reference.cells, recompute = FALSE)
    reference <- RenameCells(
      object = reference,
      new.names = gsub(pattern = "_reference$", replacement = "", x = reference.cells)
    )
    query.cells <- slot(object = anchors, name = "query.cells")
    query <- subset(x = object.list, cells = query.cells, recompute = FALSE)
    query <- RenameCells(
      object = query,
      new.names = gsub(pattern = "_query$", replacement = "", x = query.cells)
    )
  }
  object.list <- list(reference = reference, query = query)
  anchor.df <- as.data.frame(x = anchor.df)
  anchor.df$dataset1 <- "reference"
  anchor.df$dataset2 <- "query"
  anchor.df <- AnnotateAnchors(
    anchors = anchor.df,
    vars = vars,
    slot = slot,
    object.list = object.list,
    assay = assay
  )
  return(anchor.df)
}

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
# @return Anchor dataframe with additional columns corresponding to the dataset
# of each cell

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

# Adjust sample tree to only include given reference objects
#
# @param x A sample tree
# @param reference.objects a sorted list of reference object IDs
#
AdjustSampleTree <- function(x, reference.objects) {
  for (i in 1:nrow(x = x)) {
    obj.id <- -(x[i, ])
    if (obj.id[[1]] > 0) {
      x[i, 1] <- -(reference.objects[[obj.id[[1]]]])
    }
    if (obj.id[[2]] > 0) {
      x[i, 2] <- -(reference.objects[[obj.id[[2]]]])
    }
  }
  return(x)
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
# Counts anchors between each dataset and scales based on total number of cells
# in the datasets
#
# @param anchor.df Matrix of anchors
# @param offsets Dataset sizes in anchor matrix. Used to identify boundaries of
# each dataset in matrix, so that total pairwise anchors between all datasets
# can be counted
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
  slot = "data",
  integration.name = 'integrated',
  features = NULL,
  k.filter = 200,
  nn.method = "annoy",
  n.trees = 50,
  eps = 0,
  verbose = TRUE
) {
  if (verbose) {
    message("Filtering anchors")
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
  if (min(length(x = nn.cells1), length(x = nn.cells2)) < k.filter) {
    warning("Number of anchor cells is less than k.filter. Retaining all anchors.")
    k.filter <- min(length(x = nn.cells1), length(x = nn.cells2))
    anchors <- GetIntegrationData(object = object, integration.name = integration.name, slot = "anchors")
  } else {
    cn.data1 <- L2Norm(
      mat = as.matrix(x = t(x = GetAssayData(
        object = object[[assay[1]]],
        slot = slot)[features, nn.cells1])),
      MARGIN = 1)
    cn.data2 <- L2Norm(
      mat = as.matrix(x = t(x = GetAssayData(
        object = object[[assay[2]]],
        slot = slot)[features, nn.cells2])),
      MARGIN = 1)
    nn <- NNHelper(
      data = cn.data2[nn.cells2, ],
      query = cn.data1[nn.cells1, ],
      k = k.filter,
      method = nn.method,
      n.trees = n.trees,
      eps = eps
    )

    anchors <- GetIntegrationData(object = object, integration.name = integration.name, slot = "anchors")
    position <- sapply(X = 1:nrow(x = anchors), FUN = function(x) {
      which(x = anchors[x, "cell2"] == Indices(object = nn)[anchors[x, "cell1"], ])[1]
    })
    anchors <- anchors[!is.na(x = position), ]
    if (verbose) {
      message("\tRetained ", nrow(x = anchors), " anchors")
    }
  }
  object <- SetIntegrationData(
    object = object,
    integration.name = integration.name,
    slot = "anchors",
    new.data = anchors
  )
  return(object)
}



FindAnchors_v3 <- function(
  object.pair,
  assay,
  slot,
  cells1,
  cells2,
  internal.neighbors,
  reduction,
  reduction.2 = character(),
  nn.reduction = reduction,
  dims = 1:10,
  k.anchor = 5,
  k.filter = 200,
  k.score = 30,
  max.features = 200,
  nn.method = "annoy",
  n.trees = 50,
  nn.idx1 = NULL,
  nn.idx2 = NULL,
  eps = 0,
  projected = FALSE,
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
    internal.neighbors = internal.neighbors,
    dims = dims,
    reduction = reduction,
    reduction.2 = reduction.2,
    nn.reduction = nn.reduction,
    k = k.neighbor,
    nn.method = nn.method,
    n.trees = n.trees,
    nn.idx1 = nn.idx1,
    nn.idx2 = nn.idx2,
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
      max.features = max.features,
      projected = projected
    )

    if(length(top.features) == 2){
      top.features <- intersect(top.features[[1]], top.features[[2]])
    } else{
      top.features <- as.vector(top.features)
    }
    top.features <- top.features[top.features %in% rownames(x = object.pair)]
    object.pair <- FilterAnchors(
      object = object.pair,
      assay = assay,
      slot = slot,
      integration.name = 'integrated',
      features = top.features,
      k.filter = k.filter,
      nn.method = nn.method,
      n.trees = n.trees,
      eps = eps,
      verbose = verbose
    )
  }
  if (!is.na(x = k.score)) {
    object.pair = ScoreAnchors(
      object = object.pair,
      assay = DefaultAssay(object = object.pair),
      integration.name = "integrated",
      verbose = verbose,
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


FindAnchors_v5 <- function(
  object.pair,
  assay,
  slot,
  cells1,
  cells2,
  internal.neighbors,
  reduction,
  reduction.2 = character(),
  nn.reduction = reduction,
  dims = 1:10,
  k.anchor = 5,
  k.filter = 200,
  k.score = 30,
  max.features = 200,
  nn.method = "annoy",
  n.trees = 50,
  nn.idx1 = NULL,
  nn.idx2 = NULL,
  eps = 0,
  projected = FALSE,
  verbose = TRUE
) {
  ref.assay <- assay[1]
  query.assay <- assay[2]
  reference.layers <- Layers(object.pair[[ref.assay]], search = 'data')[1]
  query.layers <- setdiff(Layers(object.pair[[query.assay]], search = 'data'), reference.layers)
  anchor.list <- list()
  for (i in seq_along(query.layers)) {
    cells2.i <- Cells(
      x = object.pair[[query.assay]],
      layer = query.layers[i]
    )
    object.pair.i <- subset(
      x = object.pair,
      cells = c(cells1, cells2.i)
    )
    anchor.list[[i]] <- FindAnchors_v3(
      object.pair = object.pair.i,
      assay = assay,
      slot = slot,
      cells1 = cells1,
      cells2 = cells2.i,
      internal.neighbors = internal.neighbors,
      reduction = reduction,
      reduction.2 = reduction.2,
      nn.reduction = nn.reduction,
      dims = dims,
      k.anchor = k.anchor,
      k.filter = k.filter,
      k.score = k.score,
      max.features = max.features,
      nn.method = nn.method,
      n.trees = n.trees,
      nn.idx1 = nn.idx1,
      nn.idx2 = nn.idx2,
      eps = eps,
      projected = projected,
      verbose = verbose
    )
    anchor.list[[i]][,2] <- match(x = cells2.i, table = cells2)[anchor.list[[i]][,2]]
    anchor.list[[i]] <- t(anchor.list[[i]])
  }
  anchors <- t(x = matrix(
    data = unlist(x = anchor.list),
    nrow = 3,
    ncol = sum(
      sapply(X = anchor.list, FUN = function(x) ncol(x))
    )
  )
  )
  colnames(anchors) <- c('cell1', 'cell2', 'score')
  return(anchors)
}

FindAnchors <- function(
  object.pair,
  assay,
  slot,
  cells1,
  cells2,
  internal.neighbors,
  reduction,
  reduction.2 = character(),
  nn.reduction = reduction,
  dims = 1:10,
  k.anchor = 5,
  k.filter = 200,
  k.score = 30,
  max.features = 200,
  nn.method = "annoy",
  n.trees = 50,
  nn.idx1 = NULL,
  nn.idx2 = NULL,
  eps = 0,
  projected = FALSE,
  verbose = TRUE
) {
  if (inherits(x = object.pair[[assay[1]]], what = 'Assay')) {
    FindAnchors.function <- FindAnchors_v3
  } else if (inherits(x = object.pair[[assay[1]]], what = 'Assay5')) {
    FindAnchors.function <- FindAnchors_v5
  }
  anchors <- FindAnchors.function(
    object.pair = object.pair,
    assay = assay,
    slot = slot,
    cells1 = cells1,
    cells2 = cells2,
    internal.neighbors = internal.neighbors,
    reduction = reduction,
    reduction.2 = reduction.2,
    nn.reduction = nn.reduction,
    dims = dims,
    k.anchor = k.anchor,
    k.filter = k.filter,
    k.score = k.score,
    max.features = max.features,
    nn.method = nn.method,
    n.trees = n.trees,
    nn.idx1 = nn.idx1,
    nn.idx2 = nn.idx2,
    eps = eps,
    projected = projected,
    verbose = verbose
  )
  return(anchors)
}

# Find Anchor pairs
#
FindAnchorPairs <- function(
  object,
  integration.name = 'integrated',
  k.anchor = 5,
  verbose = TRUE
) {
  neighbors <- GetIntegrationData(object = object, integration.name = integration.name, slot = 'neighbors')
  max.nn <- c(ncol(x = neighbors$nnab), ncol(x = neighbors$nnba))
  if (any(k.anchor > max.nn)) {
    message(paste0('warning: requested k.anchor = ', k.anchor, ', only ', min(max.nn), ' in dataset'))
    k.anchor <- min(max.nn)
  }
  if (verbose) {
    message("Finding anchors")
  }
  # convert cell name to neighbor index
  nn.cells1 <- neighbors$cells1
  nn.cells2 <- neighbors$cells2
  cell1.index <-  suppressWarnings(which(colnames(x = object) == nn.cells1, arr.ind = TRUE))
  ncell <- 1:nrow(x = neighbors$nnab)
  ncell <- ncell[ncell %in% cell1.index]
  anchors <- list()
  # pre allocate vector
  anchors$cell1 <- rep(x = 0, length(x = ncell) * 5)
  anchors$cell2 <- anchors$cell1
  anchors$score <- anchors$cell1 + 1
  idx <- 0
  indices.ab <- Indices(object = neighbors$nnab)
  indices.ba <- Indices(object = neighbors$nnba)
  for (cell in ncell) {
    neighbors.ab <- indices.ab[cell, 1:k.anchor]
    mutual.neighbors <- which(
      x = indices.ba[neighbors.ab, 1:k.anchor, drop = FALSE] == cell,
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
  return(object)
}

# Find nearest neighbors
#
FindNN <- function(
  object,
  cells1 = NULL,
  cells2 = NULL,
  internal.neighbors,
  grouping.var = NULL,
  dims = 1:10,
  reduction = "cca.l2",
  reduction.2 = character(),
  nn.dims = dims,
  nn.reduction = reduction,
  k = 300,
  nn.method = "annoy",
  n.trees = 50,
  nn.idx1 = NULL,
  nn.idx2 = NULL,
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
    groups <- names(x = sort(x = table(object[[grouping.var]]), decreasing = TRUE))
    cells1 <- colnames(x = object)[object[[grouping.var]] == groups[[1]]]
    cells2 <- colnames(x = object)[object[[grouping.var]] == groups[[2]]]
  }
  if (verbose) {
    message("Finding neighborhoods")
  }
  dim.data.self <- Embeddings(object = object[[nn.reduction]])[, nn.dims]
  if (!is.null(x = internal.neighbors[[1]])) {
    nnaa <- internal.neighbors[[1]]
  } else {
    dims.cells1.self <- dim.data.self[cells1, ]
    nnaa <- NNHelper(
      data = dims.cells1.self,
      k = k + 1,
      method = nn.method,
      n.trees = n.trees,
      eps = eps,
      cache.index = TRUE,
      index = nn.idx1
    )
    nn.idx1 <- Index(object = nnaa)
  }
  if (!is.null(x = internal.neighbors[[2]])) {
    nnbb <- internal.neighbors[[2]]
  } else {
    dims.cells2.self <- dim.data.self[cells2, ]
    nnbb <- NNHelper(
      data = dims.cells2.self,
      k = k + 1,
      method = nn.method,
      n.trees = n.trees,
      eps = eps,
      cache.index = TRUE
    )
    nn.idx2 <- Index(object = nnbb)
  }
  if (length(x = reduction.2) > 0) {
    nnab <- NNHelper(
      data = Embeddings(object = object[[reduction.2]])[cells2, nn.dims],
      query = Embeddings(object = object[[reduction.2]])[cells1, nn.dims],
      k = k,
      method = nn.method,
      n.trees = n.trees,
      eps = eps,
      index = if (reduction.2 == nn.reduction) nn.idx2 else NULL
    )
    
    nnba <- NNHelper(
      data = Embeddings(object = object[[reduction]])[cells1, nn.dims],
      query = Embeddings(object = object[[reduction]])[cells2, nn.dims],
      k = k,
      method = nn.method,
      n.trees = n.trees,
      eps = eps,
      index = if (reduction == nn.reduction) nn.idx1 else NULL
    )
  } else {
    dim.data.opposite <- Embeddings(object = object[[reduction]])[ ,dims]
    dims.cells1.opposite <- dim.data.opposite[cells1, ]
    dims.cells2.opposite <- dim.data.opposite[cells2, ]
    nnab <- NNHelper(
      data = dims.cells2.opposite,
      query = dims.cells1.opposite,
      k = k,
      method = nn.method,
      n.trees = n.trees,
      eps = eps,
      index =  if (reduction == nn.reduction) nn.idx2 else NULL
    )
    nnba <- NNHelper(
      data = dims.cells1.opposite,
      query = dims.cells2.opposite,
      k = k,
      method = nn.method,
      n.trees = n.trees,
      eps = eps,
      index =  if (reduction == nn.reduction) nn.idx1 else NULL
    )
  }
  object <- SetIntegrationData(
    object = object,
    integration.name = integration.name,
    slot = 'neighbors',
    new.data = list('nnaa' = nnaa, 'nnab' = nnab, 'nnba' = nnba, 'nnbb' = nnbb, 'cells1' = cells1, 'cells2' = cells2)
  )
  return(object)
}

# @param reduction a DimReduc object containing cells in the query object
# @param reverse Compute weights matrix for reference anchors that are nearest
# to query cells. Used in mapping metric to perform projection of query cells
# back from reference space.

FindWeights <- function(
  object,
  reduction = NULL,
  assay = NULL,
  integration.name = 'integrated',
  dims = 1:10,
  features = NULL,
  k = 300,
  sd.weight = 1,
  nn.method = "annoy",
  n.trees = 50,
  eps = 0,
  reverse = FALSE,
  verbose = TRUE
) {
  if (verbose) {
    message("Finding integration vector weights")
  }
  if (is.null(x = reduction) & is.null(x = features)) {
    stop("Need to specify either dimension reduction object or a set of features")
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
  if (reverse) {
    anchors.cells2 <- nn.cells2[anchors[, "cell2"]]
    anchors.cells1 <- nn.cells1[anchors[, "cell1"]]
    to.keep <- !duplicated(x = anchors.cells1)
    anchors.cells1 <- anchors.cells1[to.keep]
    anchors.cells2 <- anchors.cells2[to.keep]
    if (is.null(x = features)) {
      data.use <- Embeddings(object = reduction)[nn.cells1, dims]
      data.use.query <- Embeddings(object = reduction)[nn.cells2, dims]
    } else {
      data.use <- t(x = GetAssayData(
        object = object,
        slot = 'data',
        assay = assay)[features, nn.cells1]
      )
      data.use.query <- t(x = GetAssayData(
        object = object,
        slot = 'data',
        assay = assay)[features, nn.cells2]
      )
    }
    knn_2_2 <- NNHelper(
      data = data.use[anchors.cells1, ],
      query = data.use.query,
      k = k,
      method = nn.method,
      n.trees = n.trees,
      eps = eps
    )
  } else {
    anchors.cells2 <- unique(x = nn.cells2[anchors[, "cell2"]])
    if (is.null(x = features)) {
      data.use <- Embeddings(reduction)[nn.cells2, dims]
    } else {
      data.use <- t(x = GetAssayData(object = object, slot = 'data', assay = assay)[features, nn.cells2])
    }
    knn_2_2 <- NNHelper(
      data = data.use[anchors.cells2, ],
      query = data.use,
      k = k,
      method = nn.method,
      n.trees = n.trees,
      eps = eps
    )
  }
  distances <- Distances(object = knn_2_2)
  distances <- 1 - (distances / distances[, ncol(x = distances)])
  cell.index <- Indices(object = knn_2_2)
  integration.matrix <- GetIntegrationData(
    object = object,
    integration.name = integration.name,
    slot = "integration.matrix"
  )

  weights <- FindWeightsC(
    cells2 = 0:(length(x = nn.cells2) - 1),
    distances = as.matrix(x = distances),
    anchor_cells2 = anchors.cells2,
    integration_matrix_rownames = rownames(x = integration.matrix),
    cell_index = cell.index,
    anchor_score = anchors[, "score"],
    min_dist = 0,
    sd = sd.weight,
    display_progress = verbose
  )
  object <- SetIntegrationData(
    object = object,
    integration.name = integration.name,
    slot = 'weights',
    new.data = weights
  )
  return(object)
}


# Find weight matrix between query and reference cells from a neighbor object
#
#
FindWeightsNN <- function(
  nn.obj,
  query.cells,
  reference.cells,
  verbose = TRUE
) {
  distances <- Distances(object = nn.obj)
  distances <- 1 - (distances / distances[, ncol(x = distances)])
  cell.index <- Indices(object = nn.obj)
  weights <- Seurat:::FindWeightsC(
    cells2 = 0:(length(query.cells) - 1),
    distances = as.matrix(x = distances),
    anchor_cells2 = reference.cells,
    integration_matrix_rownames = reference.cells,
    cell_index = cell.index,
    anchor_score = rep(1, length(reference.cells)),
    min_dist = 0,
    sd = 1,
    display_progress = verbose
  )
  colnames(weights) <- query.cells
  return(weights)
}


# Work out the anchor cell offsets for given set of cells in anchor list
#
# @param anchors A dataframe of anchors, from AnchorSet object
# @param dataset Dataset number (1 or 2)
# @param cell Cell number (1 or 2)
# @param cellnames.list List of cell names in all objects
# @param cellnames list of cell names for only the object in question
#
# @return Returns a list of offsets
#
GetCellOffsets <- function(anchors, dataset, cell, cellnames.list, cellnames) {
  cell.id <- sapply(X = 1:nrow(x = anchors), FUN = function(x) {
    cellnames.list[[anchors[, dataset+3][x]]][anchors[, cell][x]]
  })
  cell.offset <- sapply(
    X = 1:length(x = cell.id),
    FUN = function(x) {
      return(which(x = cellnames == cell.id[x]))
    }
  )
  return(cell.offset)
}

# Map queries to reference
#
# Map query objects onto assembled reference dataset
#
# @param anchorset Anchorset found by FindIntegrationAnchors
# @param reference Pre-integrated reference dataset to map query datasets to
# @param new.assay.name Name for the new assay containing the integrated data
# @param normalization.method Name of normalization method used: LogNormalize
# or SCT
# @param features Vector of features to use when computing the PCA to determine the weights. Only set
# if you want a different set from those used in the anchor finding process
# @param features.to.integrate Vector of features to integrate. By default, will use the features
# used in anchor finding.
# @param dims Number of PCs to use in the weighting procedure
# @param k.weight Number of neighbors to consider when weighting
# @param weight.reduction Dimension reduction to use when calculating anchor weights.
# This can be either:
# \itemize{
#    \item{A string, specifying the name of a dimension reduction present in all objects to be integrated}
#    \item{A vector of strings, specifying the name of a dimension reduction to use for each object to be integrated}
#    \item{NULL, in which case a new PCA will be calculated and used to calculate anchor weights}
# }
# Note that, if specified, the requested dimension reduction will only be used for calculating anchor weights in the
# first merge between reference and query, as the merged object will subsequently contain more cells than was in
# query, and weights will need to be calculated for all cells in the object.
# @param sd.weight Controls the bandwidth of the Gaussian kernel for weighting
# @param preserve.order Do not reorder objects based on size for each pairwise integration.
# @param eps Error bound on the neighbor finding algorithm (from \code{\link{RANN}})
# @param verbose Print progress bars and output
#
# @return Returns an integrated matrix
#
MapQueryData <- function(
  anchorset,
  reference,
  new.assay.name = "integrated",
  normalization.method = c("LogNormalize", "SCT"),
  features = NULL,
  features.to.integrate = NULL,
  dims = 1:30,
  k.weight = 100,
  weight.reduction = NULL,
  weights.matrix = NULL,
  no.offset = FALSE,
  sd.weight = 1,
  preserve.order = FALSE,
  eps = 0,
  verbose = TRUE
) {
  normalization.method <- match.arg(arg = normalization.method)
  reference.datasets <- slot(object = anchorset, name = 'reference.objects')
  object.list <- slot(object = anchorset, name = 'object.list')
  anchors <- slot(object = anchorset, name = 'anchors')
  features <- features %||% slot(object = anchorset, name = "anchor.features")
  features.to.integrate <- features.to.integrate %||% features
  cellnames.list <- list()
  for (ii in 1:length(x = object.list)) {
    cellnames.list[[ii]] <- colnames(x = object.list[[ii]])
  }
  if (length(x = reference.datasets) == length(x = object.list)) {
    query.datasets <- NULL
  } else {
    query.datasets <- setdiff(x = seq_along(along.with = object.list), y = reference.datasets)
  }
  my.lapply <- ifelse(
    test = verbose && nbrOfWorkers() == 1,
    yes = pblapply,
    no = future_lapply
  )
  query.corrected <- my.lapply(
    X = query.datasets,
    FUN = function(dataset1) {
      if (verbose) {
        message("\nIntegrating dataset ", dataset1, " with reference dataset")
      }
      filtered.anchors <- anchors[anchors$dataset1 %in% reference.datasets & anchors$dataset2 == dataset1, ]
      integrated <- RunIntegration(
        filtered.anchors = filtered.anchors,
        reference = reference,
        query = object.list[[dataset1]],
        new.assay.name = new.assay.name,
        normalization.method = normalization.method,
        cellnames.list = cellnames.list,
        features.to.integrate = features.to.integrate,
        weight.reduction = weight.reduction,
        weights.matrix = weights.matrix,
        no.offset = no.offset,
        features = features,
        dims = dims,
        k.weight = k.weight,
        sd.weight = sd.weight,
        eps = eps,
        verbose = verbose
      )
      return(integrated)
    }
  )
  reference.integrated <- GetAssayData(
    object = reference,
    slot = 'data'
  )[features.to.integrate, ]
  query.corrected[[length(x = query.corrected) + 1]] <- reference.integrated
  all.integrated <- do.call(cbind, query.corrected)
  return(all.integrated)
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
  nn.matrix <- as.sparse(x = nn.matrix)
  return(nn.matrix)
}

# Pairwise dataset integration
#
# Used for reference construction
#
# @param anchorset Results from FindIntegrationAnchors
# @param new.assay.name Name for the new assay containing the integrated data
# @param normalization.method Name of normalization method used: LogNormalize
# or SCT
# @param features Vector of features to use when computing the PCA to determine
# the weights. Only set if you want a different set from those used in the
# anchor finding process
# @param features.to.integrate Vector of features to integrate. By default,
# will use the features used in anchor finding.
# @param dims Number of PCs to use in the weighting procedure
# @param k.weight Number of neighbors to consider when weighting
# @param weight.reduction Dimension reduction to use when calculating anchor
# weights. This can be either:
# \itemize{
#    \item{A string, specifying the name of a dimension reduction present in
#    all objects to be integrated}
#    \item{A vector of strings, specifying the name of a dimension reduction to
#    use for each object to be integrated}
#    \item{NULL, in which case a new PCA will be calculated and used to
#    calculate anchor weights}
# }
# Note that, if specified, the requested dimension reduction will only be used
# for calculating anchor weights in the first merge between reference and
# query, as the merged object will subsequently contain more cells than was in
# query, and weights will need to be calculated for all cells in the object.
# @param sd.weight Controls the bandwidth of the Gaussian kernel for weighting
# @param sample.tree Specify the order of integration. If NULL, will compute
# automatically.
# @param preserve.order Do not reorder objects based on size for each pairwise
# integration.
# @param eps Error bound on the neighbor finding algorithm (from
# \code{\link{RANN}})
# @param verbose Print progress bars and output
#
# @return Returns a Seurat object with a new integrated Assay
#
PairwiseIntegrateReference <- function(
  anchorset,
  new.assay.name = "integrated",
  normalization.method = c("LogNormalize", "SCT"),
  features = NULL,
  features.to.integrate = NULL,
  dims = 1:30,
  k.weight = 100,
  weight.reduction = NULL,
  sd.weight = 1,
  sample.tree = NULL,
  preserve.order = FALSE,
  eps = 0,
  verbose = TRUE
) {
  object.list <- slot(object = anchorset, name = "object.list")
  reference.objects <- slot(object = anchorset, name = "reference.objects")
  features <- features %||% slot(object = anchorset, name = "anchor.features")
  features.to.integrate <- features.to.integrate %||% features
  if (length(x = reference.objects) == 1) {
    ref.obj <- object.list[[reference.objects]]
    ref.obj[[new.assay.name]] <- CreateAssayObject(
      data = GetAssayData(ref.obj, slot = 'data')[features.to.integrate, ],
      check.matrix = FALSE
    )
    DefaultAssay(object = ref.obj) <- new.assay.name
    return(ref.obj)
  }
  anchors <- slot(object = anchorset, name = "anchors")
  offsets <- slot(object = anchorset, name = "offsets")
  objects.ncell <- sapply(X = object.list, FUN = ncol)
  if (!is.null(x = weight.reduction)) {
    if (length(x = weight.reduction) == 1 | inherits(x = weight.reduction, what = "DimReduc")) {
      if (length(x = object.list) == 2) {
        weight.reduction <- list(NULL, weight.reduction)
      } else if (inherits(x = weight.reduction, what = "character")) {
        weight.reduction <- as.list(x = rep(x = weight.reduction, times = length(x = object.list)))
      } else {
        stop("Invalid input for weight.reduction. Please specify either the names of the dimension",
             "reduction for each object in the list or provide DimReduc objects.")
      }
    }
    if (length(x = weight.reduction) != length(x = object.list)) {
      stop("Please specify a dimension reduction for each object, or one dimension reduction to be used for all objects")
    }
    if (inherits(x = weight.reduction, what = "character")) {
      weight.reduction <- as.list(x = weight.reduction)
    }
    available.reductions <- lapply(X = object.list, FUN = FilterObjects, classes.keep = 'DimReduc')
    for (ii in 1:length(x = weight.reduction)) {
      if (ii == 1 & is.null(x = weight.reduction[[ii]])) next
      if (!inherits(x = weight.reduction[[ii]], what = "DimReduc")) {
        if (!weight.reduction[[ii]] %in% available.reductions[[ii]]) {
          stop("Requested dimension reduction (", weight.reduction[[ii]], ") is not present in object ", ii)
        }
        weight.reduction[[ii]] <- object.list[[ii]][[weight.reduction[[ii]]]]
      }
    }
  }
  if (is.null(x = sample.tree)) {
    similarity.matrix <- CountAnchors(
      anchor.df = anchors,
      offsets = offsets,
      obj.lengths = objects.ncell
    )
    similarity.matrix <- similarity.matrix[reference.objects, reference.objects]
    sample.tree <- BuildSampleTree(similarity.matrix = similarity.matrix)
    sample.tree <- AdjustSampleTree(x = sample.tree, reference.objects = reference.objects)
  }
  cellnames.list <- list()
  for (ii in 1:length(x = object.list)) {
    cellnames.list[[ii]] <- colnames(x = object.list[[ii]])
  }
  unintegrated <- suppressWarnings(expr = merge(
    x = object.list[[reference.objects[[1]]]],
    y = object.list[reference.objects[2:length(x = reference.objects)]]
  ))
  names(x = object.list) <- as.character(-(1:length(x = object.list)))
  if (!is.null(x = weight.reduction)) {
    names(x = weight.reduction) <- names(x = object.list)
  }
  if (verbose & (length(x = reference.objects) != length(x = object.list))) {
    message("Building integrated reference")
  }
  for (ii in 1:nrow(x = sample.tree)) {
    merge.pair <- as.character(x = sample.tree[ii, ])
    length1 <- ncol(x = object.list[[merge.pair[1]]])
    length2 <- ncol(x = object.list[[merge.pair[2]]])
    if (!(preserve.order) & (length2 > length1)) {
      merge.pair <- rev(x = merge.pair)
      sample.tree[ii, ] <- as.numeric(merge.pair)
    }
    if (!is.null(x = weight.reduction)) {
      # extract the correct dimreduc objects, in the correct order
      weight.pair <- weight.reduction[merge.pair]
    } else {
      weight.pair <- NULL
    }
    object.1 <- DietSeurat(
      object = object.list[[merge.pair[1]]],
      assays = DefaultAssay(object =  object.list[[merge.pair[1]]]),
      counts = FALSE
    )
    object.2 <- DietSeurat(
      object = object.list[[merge.pair[2]]],
      assays = DefaultAssay(object =  object.list[[merge.pair[2]]]),
      counts = FALSE
    )
    # suppress key duplication warning
    suppressWarnings(object.1[["ToIntegrate"]] <- object.1[[DefaultAssay(object = object.1)]])
    DefaultAssay(object = object.1) <- "ToIntegrate"
    object.1 <- DietSeurat(object = object.1, assays = "ToIntegrate")
    suppressWarnings(object.2[["ToIntegrate"]] <- object.2[[DefaultAssay(object = object.2)]])
    DefaultAssay(object = object.2) <- "ToIntegrate"
    object.2 <- DietSeurat(object = object.2, assays = "ToIntegrate")
    datasets <- ParseMergePair(sample.tree, ii)
    if (verbose) {
      message(
        "Merging dataset ",
        paste(datasets$object2, collapse = " "),
        " into ",
        paste(datasets$object1, collapse = " ")
      )
    }
    merged.obj <- merge(x = object.1, y = object.2, merge.data = TRUE)
    if (verbose) {
      message("Extracting anchors for merged samples")
    }
    filtered.anchors <- anchors[anchors$dataset1 %in% datasets$object1 & anchors$dataset2 %in% datasets$object2, ]
    integrated.matrix <- RunIntegration(
      filtered.anchors = filtered.anchors,
      normalization.method = normalization.method,
      reference = object.1,
      query = object.2,
      cellnames.list = cellnames.list,
      new.assay.name = new.assay.name,
      features.to.integrate = features.to.integrate,
      features = features,
      dims = dims,
      weight.reduction = weight.reduction,
      k.weight = k.weight,
      sd.weight = sd.weight,
      eps = eps,
      verbose = verbose
    )
    integrated.matrix <- cbind(integrated.matrix, GetAssayData(object = object.1, slot = 'data')[features.to.integrate, ])
    merged.obj[[new.assay.name]] <- CreateAssayObject(data = integrated.matrix, check.matrix = FALSE)
    DefaultAssay(object = merged.obj) <- new.assay.name
    object.list[[as.character(x = ii)]] <- merged.obj
    object.list[[merge.pair[[1]]]] <- NULL
    object.list[[merge.pair[[2]]]] <- NULL
    invisible(x = CheckGC())
  }
  integrated.data <- GetAssayData(
    object = object.list[[as.character(x = ii)]],
    assay = new.assay.name,
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
    misc = NULL,
    key = paste0(new.assay.name, "_")
  )
  unintegrated[[new.assay.name]] <- new.assay
  # "unintegrated" now contains the integrated assay
  DefaultAssay(object = unintegrated) <- new.assay.name
  VariableFeatures(object = unintegrated) <- features
  if (normalization.method == "SCT"){
    unintegrated[[new.assay.name]] <- SetAssayData(
      object = unintegrated[[new.assay.name]],
      slot = "scale.data",
      new.data = as.matrix(x = GetAssayData(object = unintegrated[[new.assay.name]], slot = "data"))
    )
  }
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
  unintegrated[["FindIntegrationAnchors"]] <- slot(object = anchorset, name = "command")
  suppressWarnings(expr = unintegrated <- LogSeuratCommand(object = unintegrated))
  return(unintegrated)
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


#' @rdname ProjectCellEmbeddings
#' @method ProjectCellEmbeddings Seurat
#' @export
#'
#'
ProjectCellEmbeddings.Seurat <- function(
  query,
  reference,
  query.assay = NULL,
  reference.assay = NULL,
  reduction = "pca",
  dims = 1:50,
  normalization.method = c("LogNormalize", "SCT"),
  scale = TRUE,
  verbose = TRUE,
  nCount_UMI = NULL,
  feature.mean = NULL,
  feature.sd = NULL
) {
  if (verbose) {
    message("Projecting cell embeddings")
  }
  normalization.method <- match.arg(arg = normalization.method)
  query.assay <- query.assay %||% DefaultAssay(object = query)
  reference.assay <- reference.assay %||% DefaultAssay(object = reference)
  if (normalization.method == 'SCT') {
    if (!IsSCT(assay = reference[[reference.assay]])) {
      stop('reference in ', reference.assay, ' assay does not have a SCT model' )
    }
    reference.model.num <- length(slot(object = reference[[reference.assay]], name = "SCTModel.list"))
    if (reference.model.num > 1) {
      stop("Given reference assay (", reference.assay, ") has ", reference.model.num ,
           " reference sct models. Please provide a reference assay with a ",
           " single reference sct model.", call. = FALSE)
    } else if (reference.model.num == 0) {
        stop("Given reference assay (", reference.assay,
             ") doesn't contain a reference SCT model.")
    }
  }
  proj.pca <- ProjectCellEmbeddings(
    query = query[[query.assay]],
    reference = reference,
    reference.assay = reference.assay,
    reduction = reduction,
    dims = dims,
    scale = scale,
    normalization.method = normalization.method,
    verbose = verbose,
    nCount_UMI = nCount_UMI,
    feature.mean = feature.mean,
    feature.sd = feature.sd
  )
  return(proj.pca)
}

#' @rdname ProjectCellEmbeddings
#' @method ProjectCellEmbeddings Assay
#' @export
#'
ProjectCellEmbeddings.Assay <- function(
  query,
  reference,
  reference.assay = NULL,
  reduction = "pca",
  dims = 1:50,
  scale = TRUE,
  normalization.method = NULL,
  verbose = TRUE,
  nCount_UMI = NULL,
  feature.mean = NULL,
  feature.sd = NULL
) {
  features <- Reduce(
    f = intersect,
    x = list(
      rownames(x = Loadings(object = reference[[reduction]])),
    rownames(x = reference[[reference.assay]]),
    rownames(x = query)
  )
  )
 if (normalization.method == 'SCT') {
   slot <- 'counts'
 } else {
   slot <- 'data'
 }
  proj.pca <- ProjectCellEmbeddings(
    query = GetAssayData(
      object = query,
      slot = slot),
    reference = reference,
    reference.assay = reference.assay,
    reduction = reduction,
    dims = dims,
    scale = scale,
    normalization.method = normalization.method,
    verbose = verbose,
    features = features,
    nCount_UMI = nCount_UMI,
    feature.mean = feature.mean,
    feature.sd = feature.sd
  )
  return(proj.pca)
}

#' @rdname ProjectCellEmbeddings
#' @method ProjectCellEmbeddings SCTAssay
#' @export
#'
ProjectCellEmbeddings.SCTAssay <- function(
  query,
  reference,
  reference.assay = NULL,
  reduction = "pca",
  dims = 1:50,
  scale = TRUE,
  normalization.method = NULL,
  verbose = TRUE,
  nCount_UMI = NULL,
  feature.mean = NULL,
  feature.sd = NULL
) {
  if (normalization.method != 'SCT') {
    warning('Query data is SCT normalized, but normalization.method is set to LogNormalize')
  }
  features <- Reduce(
    f = intersect,
    x = list(
      rownames(x = Loadings(object = reference[[reduction]])),
      rownames(x = reference[[reference.assay]]),
      rownames(x = query$scale.data)
    )
  )
  query.data <- GetAssayData(
    object = query,
    slot = "scale.data")[features,]
  ref.feature.loadings <- Loadings(object = reference[[reduction]])[features, dims]
  proj.pca <- t(crossprod(x = ref.feature.loadings, y = query.data))
  return(proj.pca)
}

#' @rdname ProjectCellEmbeddings
#' @method ProjectCellEmbeddings StdAssay
#' @export
#'
ProjectCellEmbeddings.StdAssay <- function(
  query,
  reference,
  reference.assay = NULL,
  reduction = "pca",
  dims = 1:50,
  scale = TRUE,
  normalization.method = NULL,
  verbose = TRUE,
  nCount_UMI = NULL,
  feature.mean = NULL,
  feature.sd = NULL
) {
  reference.assay <- reference.assay %||% DefaultAssay(object = reference)
  features <- Reduce(
    f = intersect,
    x = list(
      rownames(x = Loadings(object = reference[[reduction]])),
      rownames(x = reference[[reference.assay]])
    )
  )
  if (normalization.method == 'SCT') {
    layers.set <- Layers(object = query, search = 'counts')
  } else {
    layers.set <- Layers(object = query, search = 'data')
  }
  proj.pca.list <- list()
  cell.list <- list()
  for (i in seq_along(layers.set)) {
    proj.pca.list[[i]] <- t(ProjectCellEmbeddings(
      query = LayerData(object = query, layer = layers.set[i]),
      reference = reference,
      reference.assay = reference.assay,
      reduction = reduction,
      dims = dims,
      scale = scale,
      normalization.method = normalization.method,
      verbose = verbose,
      features = features,
      nCount_UMI = nCount_UMI[Cells(x = query, layer = layers.set[i])],
      feature.mean = feature.mean,
      feature.sd = feature.sd
    ))
    cell.list[[i]] <- colnames(proj.pca.list[[i]])
  }
  proj.pca <- matrix(
    data = unlist(proj.pca.list),
    nrow = nrow(proj.pca.list[[1]]),
    ncol = ncol(query)
    )
  rownames(proj.pca) <- rownames(proj.pca.list[[1]])
  colnames(proj.pca) <- unlist(cell.list)
  proj.pca <- t(proj.pca)
  proj.pca <- proj.pca[colnames(query),]
  return(proj.pca)
}

#' @rdname ProjectCellEmbeddings
#' @method ProjectCellEmbeddings default
#' @export
#'
ProjectCellEmbeddings.default <- function(
  query,
  reference,
  reference.assay = NULL,
  reduction = "pca",
  dims = 1:50,
  scale = TRUE,
  normalization.method = NULL,
  verbose = TRUE,
  features = NULL,
  nCount_UMI = NULL,
  feature.mean = NULL,
  feature.sd = NULL
){
  features <- features %||% rownames(x = Loadings(object = reference[[reduction]]))
if (normalization.method == 'SCT') {
  reference.SCT.model <- slot(object = reference[[reference.assay]], name = "SCTModel.list")[[1]]
  query <- FetchResiduals_reference(
    object = query,
    reference.SCT.model = reference.SCT.model,
    features = features,
    nCount_UMI = nCount_UMI)
} else {
  query <- query[features,]
  reference.data <-  GetAssayData(
    object = reference,
    assay = reference.assay,
    slot = "data")[features, ]
  if (is.null(x = feature.mean)) {
    if (inherits(x = reference.data, what = 'dgCMatrix')) {
      feature.mean <- RowMeanSparse(mat = reference.data)
    } else if (inherits(x = reference.data, what = "IterableMatrix")) {
      bp.stats <- BPCells::matrix_stats(matrix = reference.data, 
                                        row_stats = "variance")
      feature.mean <- bp.stats$row_stats["mean",]
    } else {
      feature.mean <- rowMeans2(x = reference.data)
    }
    if (scale) {
      if (inherits(x = reference.data, what = "IterableMatrix")) {
        feature.sd <- sqrt(bp.stats$row_stats["variance",])
      } else {
        feature.sd <- sqrt(x = RowVarSparse(mat = as.sparse(reference.data)))
      }
      feature.sd[is.na(x = feature.sd)] <- 1
    } else {
      feature.sd <- rep(x = 1, nrow(x = reference.data))
    }
    feature.mean[is.na(x = feature.mean)] <- 1
  }
  store.names <- dimnames(x = query)
  if (is.numeric(x = feature.mean)) {
    query <- FastSparseRowScaleWithKnownStats(
      mat = as.sparse(x = query),
      mu = feature.mean,
      sigma = feature.sd,
      display_progress = FALSE
    )
  }
  dimnames(x = query) <- store.names
}
  ref.feature.loadings <- Loadings(object = reference[[reduction]])[features, dims]
  proj.pca <- t(crossprod(x = ref.feature.loadings, y = query))
  return(proj.pca)
}

#' @rdname ProjectCellEmbeddings
#' @method ProjectCellEmbeddings IterableMatrix
#' @export
#'
#'
ProjectCellEmbeddings.IterableMatrix <- function(
  query,
  reference,
  reference.assay = NULL,
  reduction = "pca",
  dims = 1:50,
  scale = TRUE,
  normalization.method = NULL,
  verbose = TRUE,
  features = features,
  nCount_UMI = NULL,
  feature.mean = NULL,
  feature.sd = NULL,
  block.size = 10000
) {
  features <- features %||% rownames(x = Loadings(object = reference[[reduction]]))
  features <- intersect(x = features, y = rownames(query))
  if (normalization.method == 'SCT') {
    reference.SCT.model <- slot(object = reference[[reference.assay]], name = "SCTModel.list")[[1]]
    cells.grid <- split(
      x = 1:ncol(query),
      f = ceiling(seq_along(along.with = 1:ncol(query))/block.size)
      )
    proj.list <- list()
    for (i in seq_along(along.with = cells.grid)) {
      query.i <- FetchResiduals_reference(
        object = as.sparse(query[,cells.grid[[i]]]),
        reference.SCT.model = reference.SCT.model,
        features = features,
        nCount_UMI = nCount_UMI[colnames(query)[cells.grid[[i]]]])
      proj.list[[i]] <- t(Loadings(object = reference[[reduction]])[features,dims]) %*% query.i
    }
    proj.pca <- t(matrix(
      data = unlist(x = proj.list),
      nrow = length(x = dims),
      ncol = ncol(x = query),
      dimnames = list(
        colnames(x = Embeddings(object = reference[[reduction]]))[dims],
        colnames(x = query))
      ))
  } else {
    query <- query[features,]
    reference.data <- LayerData(object = reference[[reference.assay]], layer = 'data')[features,]
    if (is.null(x = feature.mean)) {
      if (inherits(x = reference.data, what = 'dgCMatrix')) {
        feature.mean <- RowMeanSparse(mat = reference.data)
      } else if (inherits(x = reference.data, what = "IterableMatrix")) {
        bp.stats <- BPCells::matrix_stats(matrix = reference.data, 
                                          row_stats = "variance")
        feature.mean <- bp.stats$row_stats["mean",]
      } else {
        feature.mean <- rowMeans(mat = reference.data)
      }
      if (scale) {
        if (inherits(x = reference.data, what = "IterableMatrix")) {
          feature.sd <- sqrt(x = bp.stats$row_stats["variance",])
        } else {
          feature.sd <- sqrt(
            x = RowVarSparse(
              mat = as.sparse(x = reference.data)
            )
          )
        }
        feature.sd[is.na(x = feature.sd)] <- 1
      } else {
        feature.sd <- rep(x = 1, nrow(x = reference.data))
      }
      feature.mean[is.na(x = feature.mean)] <- 1
    }
    query.scale <- BPCells::min_by_row(mat = query, vals = 10*feature.sd + feature.mean)
    query.scale <- (query.scale - feature.mean)/feature.sd
    proj.pca <- t(query.scale) %*% Loadings(object = reference[[reduction]])[features,dims]
    rownames(x = proj.pca) <- colnames(x = query)
    colnames(x = proj.pca) <- colnames(x = Embeddings(object = reference[[reduction]]))[dims]
  }
  return(proj.pca)
}

#' @rdname ProjectCellEmbeddings
#' @method ProjectCellEmbeddings DelayedMatrix
#' @export
#'
ProjectCellEmbeddings.DelayedMatrix <- function(
  query.data,
  block.size = 1e9,
  reference,
  assay = NULL,
  reduction,
  normalization.method = NULL,
  dims = NULL,
  feature.mean = NULL,
  feature.sd = NULL
) {
  dims <- dims %||% 1:ncol(reference[[reduction]])
  assay <- assay %||% DefaultAssay(reference)
  features <- intersect(rownames(query.data),
                        rownames(reference[[reduction]]@feature.loadings))
  query.data <- query.data[features,]
  if (IsSCT(object[[assay]])) {
    # TODO: SCT reiduals projection
  } else {
    feature.mean <- feature.mean[features] %||%
      RowMeanSparse(mat =  LayerData(object = reference[[assay]], layer = 'data')[features,])
    feature.sd <- feature.sd[features] %||%
      sqrt(RowVarSparse(mat = LayerData(object = reference[[assay]], layer = 'data')[features,]))
    feature.sd <- MinMax(feature.sd, max = max(feature.sd), min = 0.1)
    suppressMessages(setAutoBlockSize(size = block.size))
    cells.grid <- DelayedArray::colAutoGrid(x = query.data)
    emb.list <- list()
    for (i in seq_len(length.out = length(x = cells.grid))) {
      vp <- cells.grid[[i]]
      data.block <- DelayedArray::read_block(x = query.data,
                                             viewport = vp,
                                             as.sparse = TRUE)
      data.block <- apply(data.block, MARGIN = 2, function(x) {
        x <- (x - feature.mean)/feature.sd
        return(x)
      })
      emb.block <- t(reference[[reduction]]@feature.loadings[features,dims]) %*%  data.block
      emb.list[[i]] <- emb.block
    }
    # list to matrix, column has to be cells
    emb.mat <- t(matrix(data = unlist(emb.list), nrow = length(dims) , ncol = ncol(query.data)))
    rownames(emb.mat) <- colnames(query.data)
    colnames(emb.mat) <- colnames(reference[[reduction]]@cell.embeddings)[dims]
  }
  return(emb.mat)
}




# Project new data onto SVD (LSI or PCA)
#
# A = U∑V         SVD
# U' = VA'/∑      LSI projection
#
# Note that because in LSI we don't multiply by ∑ to get the embeddings (it's just U),
# we need to divide by ∑ in the projection to get the equivalent. Therefore need
# the singular values, which (in Signac RunLSI) we store in the DimReduc misc slot.
#
# @param reduction A \code{DimReduc} object containing the SVD dimension
# reduction. Assumes original irlba output is stored in the misc slot of the dimreduc.
# @param data A data matrix to project onto the SVD. Must contain the same
# features used to construct the original SVD.
# @param mode "pca" or "lsi". Determines if we divide projected values by singular values.
# @param features Features to use. If NULL, use all common features between
# the dimreduc and the data matrix.
# @param do.center Center the projected cell embeddings (subtract mean across cells)
# @param do.scale Scale the projected cell embeddings (divide by standard deviation across cells)
# @param use.original.stats When standardizing the vectors, use the mean and standard deviation
# of the original vectors from the SVD, rather than the mean and standard deviation of the
# projected vectors.
# @param dims A vector containing the dimensions to use in the projection. If NULL (default),
# project to all dimensions in the input SVD.
# @param verbose Display messages
#
# @return Returns a matrix
#' @importFrom Matrix crossprod
# @export
ProjectSVD <- function(
  reduction,
  data,
  mode = "pca",
  features = NULL,
  do.center = FALSE,
  do.scale = FALSE,
  use.original.stats = FALSE,
  dims = NULL,
  verbose = TRUE
) {
  vt <- Loadings(object = reduction)
  dims <- dims %||% seq_len(length.out = ncol(x = vt))
  features <- features %||% rownames(x = vt)
  features <- intersect(x = features, y = rownames(x = data))
  vt <- vt[features, dims]
  data <- data[features, ]
  if (verbose) {
    message("Projecting new data onto SVD")
  }
  projected.u <- as.matrix(x = crossprod(x = vt, y = data))
  if (mode == "lsi") {
    components <- slot(object = reduction, name = 'misc')
    sigma <- components$d
    projected.u <- projected.u / sigma[dims]
  }
  if (do.center) {
    if (use.original.stats) {
      components <- slot(object = reduction, name = 'misc')
      if ("u" %in% names(x = components)) {
        # preferentially use original irlba output stored in misc
        # signac scales and centers embeddings by default
        embed.mean <- apply(X = components$u, MARGIN = 2, FUN = mean)
      } else {
        # raw irlba output not stored, fall back to the reference embeddings
        ref.emb <- Embeddings(object = reduction)
        embed.mean <- apply(X = ref.emb, MARGIN = 2, FUN = mean)
      }
    } else {
      # projected.u is transposed so use MARGIN = 1
      embed.mean <- apply(X = projected.u, MARGIN = 1, FUN = mean)
    }
    projected.u <- projected.u - embed.mean
  }
  if (do.scale) {
    if (use.original.stats) {
      components <- slot(object = reduction, name = 'misc')
      if ("u" %in% names(x = components)) {
        embed.sd <- apply(X = components$u, MARGIN = 2, FUN = sd)
      } else {
        ref.emb <- Embeddings(object = reduction)
        embed.sd <- apply(X = ref.emb, MARGIN = 2, FUN = sd)
      }
    } else {
      embed.sd <- apply(X = projected.u, MARGIN = 1, FUN = sd)
    }
    projected.u <- projected.u / embed.sd
  }
  return(t(x = projected.u))
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

# Run integration between a reference and query object
#
# Should only be called from within another function
#
# @param filtered.anchors A dataframe containing only anchors between reference and query
# @param reference A reference object
# @param query A query object
# @param cellnames.list List of all cell names in all objects to be integrated
# @param new.assay.name Name for the new assay containing the integrated data
# @param features Vector of features to use when computing the PCA to determine the weights. Only set
# if you want a different set from those used in the anchor finding process
# @param features.to.integrate Vector of features to integrate. By default, will use the features
# used in anchor finding.
# @param dims Number of PCs to use in the weighting procedure
# @param k.weight Number of neighbors to consider when weighting
# @param weight.reduction Dimension reduction to use when calculating anchor weights.
# This can be either:
# \itemize{
#    \item{A string, specifying the name of a dimension reduction present in all objects to be integrated}
#    \item{A vector of strings, specifying the name of a dimension reduction to use for each object to be integrated}
#    \item{NULL, in which case a new PCA will be calculated and used to calculate anchor weights}
# }
# Note that, if specified, the requested dimension reduction will only be used for calculating anchor weights in the
# first merge between reference and query, as the merged object will subsequently contain more cells than was in
# query, and weights will need to be calculated for all cells in the object.
# @param sd.weight Controls the bandwidth of the Gaussian kernel for weighting
# @param sample.tree Specify the order of integration. If NULL, will compute automatically.
# @param eps Error bound on the neighbor finding algorithm (from \code{\link{RANN}})
# @param verbose Print progress bars and output
#
RunIntegration <- function(
  filtered.anchors,
  normalization.method,
  reference,
  query,
  cellnames.list,
  new.assay.name,
  features.to.integrate,
  weight.reduction,
  weights.matrix = NULL,
  no.offset = FALSE,
  features,
  dims,
  k.weight,
  sd.weight,
  eps,
  verbose
) {
  cells1 <- colnames(x = reference)
  cells2 <- colnames(x = query)
  if (nrow(x = filtered.anchors) < k.weight) {
    warning("Number of anchors is less than k.weight. Lowering k.weight for sample pair.")
    k.weight <- nrow(x = filtered.anchors)
  }
  merged.obj <- merge(x = reference, y = query, merge.data = TRUE)
  if (no.offset) {
    cell1.offset <- filtered.anchors[, 1]
    cell2.offset <- filtered.anchors[, 2]
  } else {
    cell1.offset <- GetCellOffsets(
      anchors = filtered.anchors,
      dataset = 1,
      cell = 1,
      cellnames.list = cellnames.list,
      cellnames = cells1
    )
    cell2.offset <- GetCellOffsets(
      anchors = filtered.anchors,
      dataset = 2,
      cell = 2,
      cellnames.list = cellnames.list,
      cellnames = cells2
    )
  }
  filtered.anchors[, 1] <- cell1.offset
  filtered.anchors[, 2] <- cell2.offset
  integration.name <- "integrated"
  merged.obj <- SetIntegrationData(
    object = merged.obj,
    integration.name = integration.name,
    slot = 'anchors',
    new.data = filtered.anchors
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
  if (is.null(x = weights.matrix)) {
    if (is.null(x = weight.reduction) && !is.null(x = dims)) {
      if (normalization.method == "SCT"){
        # recenter residuals
        centered.resids <- ScaleData(
          object = GetAssayData(object = merged.obj, assay = assay, slot = "data"),
          do.scale = FALSE,
          do.center = TRUE,
          verbose = FALSE
        )
        merged.obj[["pca"]] <- RunPCA(
          object = centered.resids[features, ],
          assay = assay,
          npcs = max(dims),
          verbose = FALSE,
          features = features
        )
      } else {
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
      }
      dr.weights <- merged.obj[['pca']]
    } else if(is.null(x = weight.reduction) && is.null(x = dims)) {
      dr.weights <- CreateDimReducObject(
        embeddings = as.matrix(x = t(x = GetAssayData(object = merged.obj))),
        key = "int_",
        assay = "ToIntegrate"
      )
      dims <- 1:ncol(x = dr.weights)
    } else {
      # need to match order of objects
      dr <- weight.reduction[[2]]
      if (!all(cells2 %in% rownames(x = dr))) {
        stop("Query cells not present in supplied DimReduc object. Set weight.reduction to a DimReduc object containing the query cells.")
      }
      if (inherits(x = dr, what = "DimReduc")) {
        dr.weights <- dr
      } else {
        dr.weights <- query[[dr]]
      }
      dims <- 1:ncol(x = dr.weights)
    }
    merged.obj <- FindWeights(
      object = merged.obj,
      integration.name = integration.name,
      reduction = dr.weights,
      dims = dims,
      k = k.weight,
      sd.weight = sd.weight,
      eps = eps,
      verbose = verbose
    )
  } else {
    merged.obj <- SetIntegrationData(
      object = merged.obj,
      integration.name = "integrated",
      slot = "weights",
      new.data = weights.matrix
    )
  }
  merged.obj <- TransformDataMatrix(
    object = merged.obj,
    new.assay.name = new.assay.name,
    features.to.integrate = features.to.integrate,
    integration.name = integration.name,
    verbose = verbose
  )
  integrated.matrix <- GetAssayData(
    object = merged.obj,
    assay = new.assay.name,
    slot = 'data'
  )
  return(integrated.matrix[, cells2])
}

# order samples based on sample tree
# the first sample is reference sample
SampleIntegrationOrder <- function(tree) {
  order <- tree[nrow(x = tree), ]
  while (sum(order > 0) != 0) {
    replace.idx <- which(x = order > 0)[1]
    replace <- tree[order[replace.idx], ]
    if (replace.idx == 1) {
      left <- vector()
      right <- order[(replace.idx + 1):length(x = order)]
      replace <- tree[order[replace.idx], ]
      order <- c(left, replace, right)
    } else if (replace.idx == length(x = order)) {
      left <- order[1:(replace.idx - 1)]
      right <- vector()
    } else {
      left <- order[1:(replace.idx - 1)]
      right <- order[(replace.idx + 1):length(x = order)]
    }
    order <- c(left, replace, right)
  }
  order <- order * (-1)
  return(order)
}

ScoreAnchors <- function(
  object,
  assay = NULL,
  integration.name = 'integrated',
  verbose = TRUE,
  k.score = 30
) {
  assay <- assay %||% DefaultAssay(object = object)
  anchor.df <- as.data.frame(x = GetIntegrationData(object = object, integration.name = integration.name, slot = 'anchors'))
  neighbors <- GetIntegrationData(object = object, integration.name = integration.name, slot = "neighbors")
  offset <- length(x = neighbors$cells1)
  indices.aa <- Indices(object = neighbors$nnaa)
  indices.bb <- Indices(object = neighbors$nnbb)
  indices.ab <- Indices(object = neighbors$nnab)
  indices.ba <- Indices(object = neighbors$nnba)
  nbrsetA <- function(x) c(indices.aa[x, 1:k.score], indices.ab[x, 1:k.score] + offset)
  nbrsetB <- function(x) c(indices.ba[x, 1:k.score], indices.bb[x, 1:k.score] + offset)
  # score = number of shared neighbors
  anchor.new <- data.frame(
    'cell1' = anchor.df[, 1],
    'cell2' = anchor.df[, 2],
    'score' = mapply(
      FUN = function(x, y) {
        length(x = intersect(x = nbrsetA(x = x), nbrsetB(x = y)))},
      anchor.df[, 1],
      anchor.df[, 2]
    )
  )
  # normalize the score
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
  return(object)
}

# Get top n features across given set of dimensions
#
# @param object Seurat object
# @param reduction Which dimension reduction to use
# @param dims Which dimensions to use
# @param features.per.dim How many features to consider per dimension
# @param max.features Number of features to return at most
# @param projected Use projected loadings
#
TopDimFeatures <- function(
  object,
  reduction,
  dims = 1:10,
  features.per.dim = 100,
  max.features = 200,
  projected = FALSE
) {
  dim.reduction <- object[[reduction]]
  max.features <- max(length(x = dims) * 2, max.features)
  num.features <- sapply(X = 1:features.per.dim, FUN = function(y) {
    length(x = unique(x = as.vector(x = sapply(X = dims, FUN = function(x) {
      unlist(x = TopFeatures(object = dim.reduction, dim = x, nfeatures = y, balanced = TRUE, projected = projected))
    }))))
  })
  max.per.pc <- which.max(x = num.features[num.features < max.features])
  features <- unique(x = as.vector(x = sapply(X = dims, FUN = function(x) {
    unlist(x = TopFeatures(object = dim.reduction, dim = x, nfeatures = max.per.pc, balanced = TRUE, projected = projected))
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

  integrated <- IntegrateDataC(integration_matrix = as.sparse(x = integration.matrix),
                               weights = as.sparse(x = weights),
                               expression_cells2 = as.sparse(x = data.use2))
  dimnames(integrated) <- dimnames(data.use2)

  new.expression <- t(rbind(data.use1, integrated))
  new.expression <- new.expression[, colnames(object)]
  new.assay <- new(
    Class = 'Assay',
    # key = paste0(new.assay.name,"_"),
    counts = new(Class = "dgCMatrix"),
    data = new.expression,
    scale.data = matrix(),
    var.features = vector(),
    meta.features = data.frame(row.names = rownames(x = new.expression)),
    misc = NULL,
    # key = paste0(new.assay.name, "_")
    key = Key(object = new.assay.name, quiet = TRUE)
  )
  object[[new.assay.name]] <- new.assay
  return(object)
}

# Helper function to validate parameters for FindTransferAnchors
#
ValidateParams_FindTransferAnchors <- function(
  reference,
  query,
  normalization.method,
  recompute.residuals,
  reference.assay,
  reference.neighbors,
  query.assay,
  reduction,
  reference.reduction,
  project.query,
  features,
  scale,
  npcs,
  l2.norm,
  dims,
  k.anchor,
  k.filter,
  k.score,
  max.features,
  nn.method,
  n.trees,
  eps,
  approx.pca,
  mapping.score.k,
  verbose
) {
  reference.assay <- reference.assay %||% DefaultAssay(object = reference)
  ModifyParam(param = "reference.assay", value = reference.assay)
  query.assay <- query.assay %||% DefaultAssay(object = query)
  ModifyParam(param = "query.assay", value = query.assay)
  DefaultAssay(object = reference) <- reference.assay
  ModifyParam(param = "reference", value = reference)
  DefaultAssay(object = query) <- query.assay
  ModifyParam(param = "query", value = query)
  if (!is.logical(x = scale)) {
    stop("Scale should be TRUE or FALSE")
  }
  if (length(x = reference) > 1 | length(x = query) > 1) {
    stop("We currently only support transfer between a single query and reference",
         call. = FALSE)
  }
  if (!reduction %in% c("pcaproject", "cca", "lsiproject", "rpca")) {
    stop("Please select either pcaproject, rpca, cca, or lsiproject for the reduction parameter.",
         call. = FALSE)
  }
  if (reduction == "cca" && !is.null(x = reference.reduction)) {
    stop("Specifying a reference reduction is only compatible with reduction = 'pcaproject'",
         call. = FALSE)
  }
  if (!normalization.method %in% c("LogNormalize", "SCT")) {
    stop("Please select either LogNormalize or SCT, for the normalization.method parameter.",
         call. = FALSE)
  }
  if (normalization.method == "SCT") {
    ModifyParam(param = "k.filter", value = NA)
  }
  if (reduction == "lsiproject") {
    ModifyParam(param = "k.filter", value = NA)
  }
  if (inherits(x = reference[[reference.assay]], what = 'Assay5') ||
      inherits(x = query[[query.assay]], what = 'Assay5')) {
    # current filter anchors not support for v5 assay
    ModifyParam(param = "k.filter", value = NA)
  }
  if (!is.na(x = k.filter) && k.filter > ncol(x = query)) {
    warning("k.filter is larger than the number of cells present in the query.\n",
            "Continuing without anchor filtering.",
            immediate. = TRUE, call. = FALSE)
    ModifyParam(param = "k.filter", value = NA)
  }
  if ((k.anchor + 1) > min(ncol(x = query), ncol(x = reference))) {
    stop("Please set k.anchor to be smaller than the number of cells in query (",
         ncol(x = query), ") and reference (", ncol(x = reference), ") objects.",
         call. = FALSE)
  }
  if ((k.score + 1) > min(ncol(x = query), ncol(x = reference))) {
    stop("Please set k.score to be smaller than the number of cells in query (",
         ncol(x = query), ") and reference (", ncol(x = reference), ") objects.",
         call. = FALSE)
  }
  if (reduction == "cca" && isTRUE(x = project.query)) {
    stop("The project.query workflow is not compatible with reduction = 'cca'",
         call. = FALSE)
  }
  if (normalization.method == "SCT" && isTRUE(x = project.query) && !IsSCT(query[[query.assay]])) {
    stop("In the project.query workflow, normalization is SCT, but query is not SCT normalized",
         call. = FALSE)
  }
  if (IsSCT(assay = query[[query.assay]]) && IsSCT(assay = reference[[reference.assay]]) &&
      normalization.method != "SCT") {
    warning("Both reference and query assays have been processed with SCTransform.",
            "Setting normalization.method = 'SCT' and continuing.")
    normalization.method <- "SCT"
    ModifyParam(param = "normalization.method", value = "SCT")
  }
  if (IsSCT(assay = query[[query.assay]]) && normalization.method == "LogNormalize") {
    stop("An SCT assay (", query.assay, ") was provided for query.assay but ",
         "normalization.method was set as LogNormalize", call. = FALSE)
  }
  if (IsSCT(assay = query[[query.assay]]) && !inherits(x = query[[query.assay]], what = "SCTAssay")) {
    query[[query.assay]] <- as(object = query[[query.assay]], Class = "SCTAssay")
    ModifyParam(param = "query", value = query)
  }
  if (IsSCT(assay = reference[[reference.assay]]) && !inherits(x = reference[[reference.assay]], what = "SCTAssay")) {
    reference[[reference.assay]] <- as(object = reference[[reference.assay]], Class = "SCTAssay")
    ModifyParam(param = "reference", value = reference)
  }
  if (normalization.method != "SCT") {
    recompute.residuals <- FALSE
    ModifyParam(param = "recompute.residuals", value = recompute.residuals)
  }
  if (recompute.residuals) {
    # recompute.residuals only happens in ProjectCellEmbeddings, so k.filter set to NA.
    ModifyParam(param = "k.filter", value = NA)
    reference.model.num <- length(x = slot(object = reference[[reference.assay]], name = "SCTModel.list"))
    if (reference.model.num > 1) {
      stop("Given reference assay (", reference.assay, ") has ", reference.model.num ,
           " reference sct models. Please provide a reference assay with a ",
           " single reference sct model.", call. = FALSE)
    } else if (reference.model.num == 0) {
      if (IsSCT(query[[query.assay]])) {
        stop("Given reference assay (", reference.assay,
             ") doesn't contain a reference SCT model.\n",
             "Query assay is a SCTAssay. ",
             "You can set recompute.residuals to FALSE ",
             "to use Query residuals to continue the analysis",
             call. = FALSE)
      }
      stop("Given reference assay (", reference.assay,
           ") doesn't contain a reference SCT model. ",
           call. = FALSE)
    } else if (reference.model.num == 1) {
      new.sct.assay <- reference.assay
      if (verbose) {
        message("Normalizing query using reference SCT model")
      }
    }
    query.umi.assay <- query.assay
    if (IsSCT(assay = query[[query.assay]])) {
      query.sct.models <- slot(object = query[[query.assay]], name = "SCTModel.list")
      query.umi.assay <- unique(x = unname(obj = unlist(x = lapply(X = query.sct.models, FUN = slot, name = "umi.assay"))))
      if (length(x = query.umi.assay) > 1) {
        stop("Query assay provided is an SCTAssay with multiple different original umi assays", call = FALSE)
      }
      if (!query.umi.assay %in% Assays(object = query)) {
        stop("Query assay provided is an SCTAssay based on an orignal UMI assay",
             " that is no longer present in the query Seurat object. Unable to",
             " recompute residuals based on the reference SCT model.\n",
             "If you want to use Query SCTAssay residuals to continue the analysis, ",
             "you can set recompute.residuals to FALSE", call. = FALSE)
      }
    }
    if (reduction %in% c('cca', 'rpca')) {
      query <- SCTransform(
        object = query,
        reference.SCT.model = slot(object = reference[[reference.assay]], name = "SCTModel.list")[[1]],
        residual.features = features,
        assay = query.umi.assay,
        new.assay.name = new.sct.assay,
        verbose = FALSE
      )
    } else {
      new.sct.assay <- query.umi.assay
    }
    
    
    DefaultAssay(query) <- new.sct.assay
    ModifyParam(param = "query.assay", value = new.sct.assay)
    ModifyParam(param = "query", value = query)
    ModifyParam(param = "reference", value = reference)
  }
  if (IsSCT(assay = reference[[reference.assay]]) && normalization.method == "LogNormalize") {
    stop("An SCT assay (", reference.assay, ") was provided for reference.assay but ",
         "normalization.method was set as LogNormalize.", call. = FALSE)
  }
  if (!IsSCT(assay = reference[[reference.assay]]) && normalization.method == "SCT") {
    stop("Given reference.assay (", reference.assay, ") has not been processed with ",
         "SCTransform. Please either run SCTransform or set normalization.method = 'LogNormalize'.",
         call. = FALSE)
  }
  # Make data slot if DNE
  if (inherits(x = query[[query.assay]], what = "Assay5")){
    if (is.null(
      tryCatch(expr = Layers(object = query[[query.assay]], search = 'data'),
               error = function (e) return(NULL))
        )
    ) {
      LayerData(
        object = query[[query.assay]],
        layer = "data") <- sparseMatrix(
          i = 1,
          j = 1,
          x = 1,
          dims = c(nrow(x = query[[query.assay]]),
                   ncol(x = query[[query.assay]])
                   )
          )
      ModifyParam(param = "query", value = query)
    }
  }
  # features must be in both reference and query
  query.assay.check <- query.assay
  reference.assay.check <- reference.assay
  ref.features <- rownames(x = reference[[reference.assay.check]])
  query.features <- rownames(x = query[[query.assay.check]])
  if (normalization.method == "SCT") {
    if (IsSCT(query[[query.assay.check]])) {
      query.features <- rownames(x = query[[query.assay.check]]$scale.data)
    }
    query.model.features <- rownames(x = Misc(object = query[[query.assay]])$vst.out$gene_attr)
    query.features <- unique(c(query.features, query.model.features))
    ref.model.features <- rownames(x = Misc(object = reference[[reference.assay]])$vst.out$gene_attr)
    ref.features <- unique(c(ref.features, ref.model.features))
  }
  if (!is.null(x = features)) {
    if (project.query) {
      features.new <- intersect(x = features, y = ref.features)
    } else {
      features.new <- intersect(x = features, y = query.features)
    }
    if (length(x = features.new) != length(x = features)) {
      warning(length(x = features) - length(x = features.new), " features of ",
              "the features specified were not present in both the reference ",
              "query assays. \nContinuing with remaining ", length(x = features.new),
              " features.", immediate. = TRUE, call. = FALSE)
      features <- features.new
    }
  } else {
    if (project.query) {
      features <- intersect(
        x = VariableFeatures(object = query[[query.assay]]),
        y = ref.features
      )
    } else {
      features <- intersect(
        x = VariableFeatures(object = reference[[reference.assay]]),
        y = query.features
      )
    }
  }
  if (length(x = features) == 0) {
    stop("No features to use in finding transfer anchors. To troubleshoot, try ",
         "explicitly providing features to the features parameter and ensure that ",
         "they are present in both reference and query assays.", call. = FALSE)
  }
  ModifyParam(param = "features", value = features)
  if (!is.null(x = reference.reduction)) {
    if (project.query) {
      if (!reference.reduction %in% Reductions(object = query)){
        stop("reference.reduction (", reference.reduction, ") is not present in ",
             "the provided query object (Note: project.query was set to TRUE).",
             call. = FALSE)
      }
      if (ncol(x = reference[[reference.reduction]]) < max(dims)) {
        stop("reference.reduction (", reference.reduction, ") does not contain ",
             "all the dimensions required by the dims parameter (Note: ",
             "project.query was set to TRUE).", call. = FALSE)
      }
    } else {
      if (!reference.reduction %in% Reductions(object = reference)){
        stop("reference.reduction (", reference.reduction, ") is not present in ",
             "the provided reference object.", call. = FALSE)
      }
      if (ncol(x = reference[[reference.reduction]]) < max(dims)) {
        stop("reference.reduction (", reference.reduction, ") does not contain ",
             "all the dimensions required by the dims parameter.", call. = FALSE)
      }
    }
  } else {
    if (reduction == "lsiproject") {
      stop("Must supply a reference reduction if reduction='lsiproject'")
    }
    mdim <- max(dims)
    if (npcs < mdim) {
      warning("npcs is smaller than the largest value requested by the dims ",
              "parameter.\nSetting npcs to ", mdim, " and continuing.",
              immediate. = TRUE, call. = FALSE)
      ModifyParam(param = "npcs", value = mdim)
      if (mdim >= length(x = features)) {
        stop("npcs (", npcs, ") must be smaller than the number of features (",
             length(x = features), "). Please either lower the npcs and/or dims ",
             "parameter settings or increase the size of the feature set.",
             call. = FALSE)
      }
    }
  }
  if (!is.null(x = reference.neighbors)) {
    if (!reference.neighbors %in% Neighbors(object = reference)) {
      stop("Specified reference.neighbors (", reference.neighbors, ") is not ",
           "available in the provided reference object.", call. = FALSE)
    }
    k.nn <- max(k.score, k.anchor)
    if (ncol(x = Indices(reference[[reference.neighbors]])) < (k.nn + 1)){
      stop("k.score or k.anchor is larger than the number of neighbors ",
           "contained in reference.nn. Recompute reference.nn using ",
           "FindNeighbors with k > k.score and k > k.anchor", call. = FALSE)
    }
  }
}

# Helper function to validate parameters for TransferData
#
ValidateParams_TransferData <- function(
  anchorset,
  combined.ob,
  anchors,
  reference.cells,
  query.cells,
  reference,
  query,
  refdata,
  weight.reduction,
  l2.norm,
  dims,
  k.weight,
  sd.weight,
  eps,
  n.trees,
  verbose,
  slot,
  only.weights,
  prediction.assay,
  label.transfer
) {
  ## check refdata
  if (is.null(refdata)) {
    if (!only.weights) {
      stop("refdata is NULL and only.weights is FALSE")
    }
  } else {
    if (!inherits(x = refdata, what = "list")) {
      refdata <- list(id = refdata)
    }
    for (i in 1:length(x = refdata)) {
      if (inherits(x = refdata[[i]], what = c("character", "factor"))) {
        # check is it's in the reference object
        if (length(x = refdata[[i]]) == 1) {
          if (is.null(x = reference)) {
            warning("If providing a single string to refdata element number ", i,
                    ", please provide the reference object. Skipping element ", i,
                    ".", call. = FALSE, immediate. = TRUE)
            refdata[[i]] <- FALSE
            next
          }
          if (refdata[[i]] %in% Assays(object = reference)) {
            refdata[[i]] <- GetAssayData(object = reference, assay = refdata[[i]])
            colnames(x = refdata[[i]]) <- paste0(colnames(x = refdata[[i]]), "_reference")
            label.transfer[[i]] <- FALSE
            next
          } else if (refdata[[i]] %in% colnames(x = reference[[]])) {
            refdata[[i]] <- reference[[refdata[[i]]]][, 1]
          } else {
            warning("Element number ", i, " provided to refdata does not exist in ",
                    "the provided reference object.", call. = FALSE, immediate. = TRUE)
            refdata[[i]] <- FALSE
            next
          }
        } else if (length(x = refdata[[i]]) != length(x = reference.cells)) {
          warning("Please provide a vector that is the same length as the number ",
                  "of reference cells used in anchor finding.\n",
                  "Length of vector provided: ", length(x = refdata[[i]]), "\n",
                  "Length of vector required: ", length(x = reference.cells),
                  "\nSkipping element ", i, ".", call. = FALSE, immediate. = TRUE)
          refdata[[i]] <- FALSE
        }
        label.transfer[[i]] <- TRUE
      } else if (inherits(x = refdata[[i]], what = c("dgCMatrix", "matrix"))) {
        if (ncol(x = refdata[[i]]) != length(x = reference.cells)) {
          warning("Please provide a matrix that has the same number of columns as ",
                  "the number of reference cells used in anchor finding.\n",
                  "Number of columns in provided matrix : ", ncol(x = refdata[[i]]), "\n",
                  "Number of columns required           : ", length(x = reference.cells),
                  "\nSkipping element ", i, ".", call. = FALSE, immediate. = TRUE)
          refdata[[i]] <- FALSE
        } else {
          colnames(x = refdata[[i]]) <- paste0(colnames(x = refdata[[i]]), "_reference")
          if (any(!colnames(x = refdata[[i]]) == reference.cells)) {
            if (any(!colnames(x = refdata[[i]]) %in% reference.cells) || any(!reference.cells %in% colnames(x = refdata[[i]]))) {
              warning("Some (or all) of the column names of the provided refdata ",
                      "don't match the reference cells used in anchor finding ",
                      "\nSkipping element", i, ".", call. = FALSE, immediate. = TRUE)
              refdata[[i]] <- FALSE
            } else {
              refdata[[i]] <- refdata[[i]][, reference.cells]
            }
          }
        }
        if (!slot %in% c("counts", "data")) {
          stop("Please specify slot as either 'counts' or 'data'.")
        }
        label.transfer[[i]] <- FALSE
      } else {
        warning("Please provide either a vector (character or factor) for label ",
                "transfer or a matrix for feature transfer. \nType provided: ",
                class(x = refdata[[i]]))
        refdata[[i]] <- FALSE
      }
      if (names(x = refdata)[i] == "") {
        possible.names <- make.unique(names = c(names(x = refdata), paste0("e", i)))
        names(x = refdata)[i] <- possible.names[length(x = possible.names)]
        if (verbose) {
          message("refdata element ", i, " is not named. Setting name as ", names(x = refdata)[i])
        }
      }
    }
    ModifyParam(param = "label.transfer", value = label.transfer)
    if (all(unlist(x = lapply(X = refdata, FUN = isFALSE)))) {
      stop("None of the provided refdata elements are valid.", call. = FALSE)
    }
    ModifyParam(param = "refdata", value = refdata)
  }




  object.reduction <- Reductions(object = slot(object = anchorset, name = "object.list")[[1]])
  valid.weight.reduction <- c("pcaproject", "pca", "cca", "rpca.ref","lsiproject", "lsi", object.reduction)
  if (!inherits(x = weight.reduction, "DimReduc")) {
    if (!weight.reduction %in% valid.weight.reduction) {
      stop("Please provide one of ", paste(valid.weight.reduction, collapse = ", "), " or a custom DimReduc to ",
           "the weight.reduction parameter.", call. = FALSE)
    }
    if (weight.reduction %in% c("pcaproject", "cca", "rpca.ref", "lsiproject") &&
        !weight.reduction %in% Reductions(object = combined.ob)) {
      stop("Specified weight.reduction (", weight.reduction, ") is not present ",
           "in the provided anchorset.", call. = FALSE)
    }
    if (weight.reduction %in% c("pca", "lsi") && is.null(x = query)) {
      stop("To use an internal PCA on the query only for weight.reduction, ",
           "please provide the query object.", call. = FALSE)
    }
  }
  if (inherits(x = weight.reduction, "DimReduc")) {
    if (is.null(x = dims)) {
      stop("Please specify dims", call. = FALSE)
    }
    if (max(dims) > ncol(x = weight.reduction)) {
      stop("The max of dims specified (", max(dims), ") is greater than the ",
           "number of dimensions in the given DimReduc (",
           ncol(x = weight.reduction), ").", call. = FALSE)
    }
  } else {
    if (is.null(x = dims) && !is.null(x = slot(object = anchorset, name = "command")$dims)) {
     ModifyParam(param = "dims", value = 1:length(x = slot(object = anchorset, name = "command")$dims))
    }
  }

  if (!is.null(x = query)) {
    if (!isTRUE(x = all.equal(
      target = gsub(pattern = "_query", replacement = "", x = query.cells),
      current = colnames(x = query),
      check.attributes = FALSE)
      )) {
      stop("Query object provided contains a different set of cells from the ",
           "query used to construct the AnchorSet provided.", call. = FALSE)
    }
  }
  if(k.weight > nrow(x = anchors)) {
    stop("Please set k.weight to be smaller than the number of anchors (",
         nrow(x = anchors), ").", call. = FALSE)
  }
}

# Internal function to validate the parameters for IntegrateEmbeddings run on
# an IntegrationAnchorSet object
#
ValidateParams_IntegrateEmbeddings_IntegrationAnchors <- function(
  anchorset,
  object.list,
  reductions,
  dims.to.integrate,
  k.weight,
  weight.reduction,
  sample.tree
) {
  nobs <- length(x = object.list)
  if (is.null(x = reductions)) {
    stop("Must supply reductions to integrate")
  }
  if (!inherits(x = reductions, what = "DimReduc")) {
    stop("Please provide a single pre-computed DimReduc object to the ",
         "reductions parameter", call. = FALSE)
  } else {
    all.cells <- make.unique(names = unname(obj = do.call(
      what = c,
      args = lapply(X = object.list, FUN = Cells)))
    )
    if (nrow(x = reductions) != length(x = all.cells)) {
      stop("The number of cells in the reduction provided (", nrow(x = reductions),
           ") doesn't match the number of cells in the objects used to build the ",
           "AnchorSet (", length(x = all.cells), ").", call. = FALSE)
    }
    if (!all(Cells(x = reductions) %in% all.cells)) {
      stop("The cell names in the reduction provided don't match the cell names ",
           "present in the objects used to build the AnchorSet", call. = FALSE)
    }
    dims.to.integrate <- dims.to.integrate %||% 1:ncol(x = reductions)
    if (max(dims.to.integrate) > ncol(x = reductions)) {
      warning("Max dims.to.integrate is larger than the number of dimensions in ",
              "the provided reduction. Setting dims.to.integrate to 1:",
              ncol(x = reductions), " and continuing.", immediate. = TRUE, call. = FALSE)
    }
    ModifyParam(param = 'dims.to.integrate', value = 1:ncol(x = reductions))
  }
  if (!is.null(x = weight.reduction)) {
    if (inherits(x = weight.reduction, what = "character")) {
      if (length(x = weight.reduction) == 1) {
        weight.reduction <- as.list(x = rep(x = weight.reduction, times = nobs))
      }
      ModifyParam(param = 'weight.reduction', value = weight.reduction)
      for (i in 1:nobs) {
        if (!weight.reduction[[i]] %in% Reductions(object = object.list[[i]])) {
          stop("weight.reduction (", weight.reduction[[i]], ") is not present ",
               "in object number ", i, ".", call. = FALSE)
        }
      }
    }
    if (inherits(x = weight.reduction[[1]], what = "DimReduc")) {
      if (length(x = weight.reduction) != nobs) {
        stop("Please provide one weight.reduction for each object. ",
             length(x = weight.reduction), " provided, ", nobs, " required.",
             call. = FALSE)
      }
      for (i in 1:nobs) {
        if (!isTRUE(all.equal(
          target = Cells(x = weight.reduction[[i]]),
          current = Cells(x = object.list[[i]]),
          check.attributes = FALSE
          ))
        ) {
          stop("Cell names in the provided weight.reduction ", i, " don't ",
               "match with the cell names in object ", i, ".", call. = FALSE)
        }
      }
    }
  }
  min.object.size <- min(sapply(X = object.list, FUN = ncol))
  if (k.weight > min.object.size) {
    stop("k.weight (", k.weight, ") is set larger than the number of cells in ",
         "the smallest object (", min.object.size, "). Please choose a smaller ",
         "k.weight.", call. = FALSE)
  }
  if (!is.null(x = sample.tree)) {
    if (ncol(x = sample.tree) != 2) {
      stop("Invalid sample tree. Please provide a two column matrix specifying
           the order of integration.")
    }
    if (min(sample.tree) < (-1 * nobs)) {
      stop("Invalid sample tree. Dataset index greater than the number of ",
           "objects was provided.")
    }
  }
}

# Internal function to validate the parameters for IntegrateEmbeddings run on
# a TransferAnchorSet object
#
ValidateParams_IntegrateEmbeddings_TransferAnchors <- function(
  anchorset,
  combined.object ,
  reference,
  query,
  reductions,
  dims.to.integrate,
  k.weight,
  weight.reduction,
  reuse.weights.matrix
) {
  if (missing(x = reference)) {
    stop("Please provide the reference object.", call. = FALSE)
  }
  if (missing(x = query)) {
    stop("Please provide the query object.", call. = FALSE)
  }
  reference.cells <- slot(object = anchorset, name = "reference.cells")
  reference.cells <- gsub(pattern = "_reference", replacement = "", x = reference.cells)
  if (!isTRUE(x = all.equal(target = reference.cells, current = Cells(x = reference), check.attributes = FALSE))) {
    stop("The set of cells used as a reference in the AnchorSet does not match ",
         "the set of cells provided in the reference object.")
  }
  query.cells <- slot(object = anchorset, name = "query.cells")
  query.cells <- gsub(pattern = "_query", replacement = "", x = query.cells)
  if (!isTRUE(x = all.equal(target = query.cells, current = colnames(x = query), check.attributes = FALSE))) {
    stop("The set of cells used as a query in the AnchorSet does not match ",
         "the set of cells provided in the query object.")
  }
  if (length(x = reductions) != 1) {
    stop("Please provide a single reduction name to reductions that is present ",
         "in the anchorset.", call. = FALSE)
  }
  if (!reductions %in% Reductions(object = combined.object)) {
    stop("Please specify a reduction that is present in the anchorset: ",
         paste(Reductions(object = combined.object), collapse = ", "), call. = FALSE)
  }
  reference <- RenameCells(object = reference, new.names = paste0(Cells(x = reference), "_reference"))
  reference.embeddings <- Embeddings(object = combined.object[[reductions]])[Cells(x = reference), ]
  reference[[reductions]] <- CreateDimReducObject(embeddings = reference.embeddings, assay = DefaultAssay(object = reference))
  ModifyParam(param = "reference", value = reference)
  query <- RenameCells(object = query, new.names = paste0(Cells(x = query), "_query"))
  query.embeddings <- Embeddings(object = combined.object[[reductions]])[Cells(x = query), ]
  query[[reductions]] <- CreateDimReducObject(embeddings = query.embeddings, assay = DefaultAssay(object = query))
  ModifyParam(param = "query", value = query)
  ModifyParam(param = "reductions", value = c(reductions, reductions))
  min.ndim <- min(ncol(x = query[[reductions[2]]]), ncol(x = reference[[reductions[1]]]))
  if (is.null(x = dims.to.integrate)) {
    dims.to.integrate <- 1:min.ndim
  } else {
    if (max(dims.to.integrate) > min.ndim) {
      dims.to.integrate <- dims.to.integrate[dims.to.integrate <= min.ndim]
      warning("Max dims.to.integrate is larger than the max dims for at least ",
              "one of the reductions specified. Setting dims.to.integrate to ",
              paste(dims.to.integrate, collapse = ","), " and continuing.",
              immediate. = TRUE, call. = FALSE)
    }
  }
  ModifyParam(param = "dims.to.integrate", value = dims.to.integrate)
  if (isTRUE(x = reuse.weights.matrix)) {
    weights.matrix <- Tool(object = query, slot = "TransferData")$weights.matrix
    if (is.null(x = weights.matrix)) {
      message("Requested to reuse weights matrix, but no weights found. Computing new weights.")
      reuse.weights.matrix <- FALSE
    } else if (nrow(x = weights.matrix) != nrow(x = slot(object = anchorset, name = "anchors"))) {
      stop("The number of anchors in the weights matrix stored in the query (",
           nrow(x = weights.matrix), ") doesn't match the number of anchors ",
           "in the anchorset (", nrow(x = slot(object = anchorset, name = "anchors")),
           ").", call. = FALSE)
    } else {
      ModifyParam(param = 'weights.matrix', value = weights.matrix)
    }
  }
  # check T/F again due to possible modification in above
  if (isFALSE(x = reuse.weights.matrix)) {
    if (k.weight > ncol(x = query)) {
      stop("k.weight (", k.weight, ") is set larger than the number of cells in ",
           "the query object (", ncol(x = query), "). Please choose a smaller ",
           "k.weight.", call. = FALSE)
    }
    if (inherits(x = weight.reduction, what = "list")) {
      if (length(x = weight.reduction) > 2) {
        stop("Supplied too many dimension reduction objects for weight.reduction. ",
             "Should supply a single DimReduc object.")
      }
      if (length(x = weight.reduction) == 2) {
        # take the second element as the dimreduc to use for query
        weight.reduction <- weight.reduction[[2]]
      }
    }
    if (inherits(x = weight.reduction, what = "character")) {
      if (length(x = weight.reduction) > 2) {
        stop("Supplied too many dimension reduction names for weight.reduction. ",
             "Should supply the name of a single DimReduc present in the query.")
      }
      if (length(x = weight.reduction) == 2) {
        # take the second element as the dimreduc to use for query
        weight.reduction <- weight.reduction[[2]]
      }
      if (!weight.reduction %in% Reductions(object = query)) {
        stop("The weight.reduction ", weight.reduction, " is not present in the ",
             "query object.", call. = FALSE)
      }
      ModifyParam(param = 'weight.reduction', value = list(NULL, query[[weight.reduction]]))
    }
    if (inherits(x = weight.reduction, what = "DimReduc")) {
      weight.reduction <- RenameCells(object = weight.reduction, new.names = paste0(Cells(x = weight.reduction), "_query"))
      if (!isTRUE(all.equal(
        target = Cells(x = weight.reduction),
        current = Cells(x = query),
        check.attributes = FALSE
      ))) {
        stop("Cell names in the provided weight.reduction  don't ",
             "match with the cell names in the query object.", call. = FALSE)
      }
      ModifyParam(param = 'weight.reduction', value = list(NULL, weight.reduction))
    }
  }
}


#' Convert Neighbor class to an asymmetrical Graph class
#'
#'
#' @param nn.object A neighbor class object
#' @param col.cells Cells names of the neighbors, cell names in nn.object is used by default
#' @param weighted Determine if use distance in the Graph
#'
#' @return Returns a Graph object
#'
#' @importFrom Matrix sparseMatrix
#'
#' @export
#'
NNtoGraph <- function(
  nn.object,
  col.cells = NULL,
  weighted = FALSE
) {
  select_nn <- Indices(object = nn.object)
  col.cells <- col.cells %||% Cells(x = nn.object)
  ncol.nn <- length(x = col.cells)
  k.nn <- ncol(x = select_nn)
  j <- as.numeric(x = t(x = select_nn))
  i <- ((1:length(x = j)) - 1) %/% k.nn + 1
  if (weighted) {
    select_nn_dist <- Distances(object = nn.object)
    dist.element <- as.numeric(x = t(x = select_nn_dist))
    nn.matrix <- sparseMatrix(
      i = i,
      j = j,
      x = dist.element,
      dims = c(nrow(x = select_nn), ncol.nn)
    )
  } else {
    nn.matrix <- sparseMatrix(
      i = i,
      j = j,
      x = 1,
      dims = c(nrow(x = select_nn), ncol.nn)
    )
  }
  rownames(x = nn.matrix) <- Cells(x = nn.object)
  colnames(x = nn.matrix) <- col.cells
  nn.matrix <- as.Graph(x = nn.matrix)
  return(nn.matrix)
}


# Find Anchor directly from assay
#
#
# @return Returns a TranserAnchor or Integration set
FindAssayAnchor <- function(
  object.list,
  reference = NULL,
  anchor.type = c("Transfer", "Integration"),
  assay = "Bridge",
  slot = "data",
  reduction =  NULL,
  k.anchor = 20,
  k.score = 50,
  verbose = TRUE
) {
  anchor.type <- match.arg(arg = anchor.type)
  reduction.name <- reduction %||% paste0(assay, ".reduc")
  if ( is.null(x = reduction) || !reduction %in% Reductions(object.list[[1]])) {
    object.list <- lapply(object.list, function(x) {
      if (is.null(reduction)) {
        x[[reduction.name]] <- CreateDimReducObject(
          embeddings = t(GetAssayData(
            object = x,
            slot = slot,
            assay = assay
          )),
          key = "L_",
          assay = assay
        )
      }
    DefaultAssay(x) <- assay
    x <- DietSeurat(x, assays = assay, dimreducs = reduction.name)
    return(x)
    }
  )
}
    object.both <- merge(object.list[[1]], object.list[[2]], merge.dr = reduction.name)
    objects.ncell <- sapply(X = object.list, FUN = function(x) dim(x = x)[2])
    offsets <- as.vector(x = cumsum(x = c(0, objects.ncell)))[1:length(x = object.list)]
    if (verbose) {
      message("Finding ",  anchor.type," anchors from assay ", assay)
    }
    anchors <- FindAnchors(object.pair = object.both,
                           assay = c(DefaultAssay(object.both), DefaultAssay(object.both)),
                           slot = 'data',
                           cells1 = colnames(object.list[[1]]),
                           cells2 = colnames(object.list[[2]]),
                           internal.neighbors = NULL,
                           reduction = reduction.name,
                           k.anchor = k.anchor,
                           k.score = k.score,
                           dims = 1:ncol(object.both[[reduction.name]]),
                           k.filter = NA,
                           verbose = verbose
    )
    inte.anchors <- anchors
    inte.anchors[, 1] <- inte.anchors[, 1] + offsets[1]
    inte.anchors[, 2] <- inte.anchors[, 2] + offsets[2]
    # determine all anchors
    inte.anchors <- rbind(inte.anchors, inte.anchors[, c(2, 1, 3)])
    inte.anchors <- AddDatasetID(
      anchor.df = inte.anchors,
      offsets = offsets,
      obj.lengths = objects.ncell
      )
    command <- LogSeuratCommand(object = object.list[[1]], return.command = TRUE)
    anchor.features <- rownames(object.both)
    if (anchor.type == "Integration") {
      anchor.set <- new(Class = "IntegrationAnchorSet",
                        object.list = object.list,
                        reference.objects = reference %||% seq_along(object.list),
                        anchors = inte.anchors,
                        weight.reduction = object.both[[reduction.name]],
                        offsets = offsets,
                        anchor.features = anchor.features,
                        command = command
      )
    } else if (anchor.type == "Transfer") {
      reference.index <- reference
      reference <- object.list[[reference.index]]
      query  <- object.list[[setdiff(c(1,2), reference.index)]]
      query <- RenameCells(
        object = query,
        new.names = paste0(Cells(x = query), "_", "query")
      )
      reference <- RenameCells(
        object = reference,
        new.names = paste0(Cells(x = reference), "_", "reference")
      )
      combined.ob <- suppressWarnings(expr = merge(
        x = reference,
        y = query,
        merge.dr = reduction.name
      ))
      anchor.set <- new(
        Class = "TransferAnchorSet",
        object.list = list(combined.ob),
        reference.cells = colnames(x = reference),
        query.cells = colnames(x = query),
        anchors = anchors,
        anchor.features = anchor.features,
        command = command
      )
    }
    return(anchor.set)
}


#' Construct a dictionary representation for each unimodal dataset
#'
#'
#' @param object.list A list of Seurat objects
#' @param bridge.object A multi-omic bridge Seurat which is used as the basis to
#' represent unimodal datasets
#' @param object.reduction A list of dimensional reductions from object.list used
#' to be reconstructed by bridge.object
#' @param bridge.reduction A list of dimensional reductions from bridge.object used
#' to reconstruct object.reduction
#' @param laplacian.reduction Name of bridge graph laplacian dimensional reduction
#' @param laplacian.dims Dimensions used for bridge graph laplacian dimensional reduction
#' @param bridge.assay.name Assay name used for bridge object reconstruction value (default is 'Bridge')
#' @param return.all.assays Whether to return all assays in the object.list.
#' Only bridge assay is returned by default.
#' @param l2.norm Whether to l2 normalize the dictionary representation
#' @param verbose Print messages and progress
#'
#' @importFrom MASS ginv
#' @return Returns a object list in which each object has a bridge cell derived assay
#' @export
#'
BridgeCellsRepresentation <- function(object.list,
                                      bridge.object,
                                      object.reduction,
                                      bridge.reduction,
                                      laplacian.reduction = 'lap',
                                      laplacian.dims = 1:50,
                                      bridge.assay.name = "Bridge",
                                      return.all.assays = FALSE,
                                      l2.norm = TRUE,
                                      verbose = TRUE
) {
  my.lapply <- ifelse(
    test = verbose && nbrOfWorkers() == 1,
    yes = pblapply,
    no = future_lapply
  )
  if (verbose) {
    message("Constructing Bridge-cells representation")
  }
  single.object = FALSE
  if (length(x = object.list) == 1 &
      inherits(x = object.list, what = 'Seurat')
  ) {
    object.list <- list(object.list)
    single.object = TRUE
  }
  dims.list <- list()
  for (i in 1:length(object.reduction)) {
   ref.dims <- list(
    object= Misc(object.list[[i]][[object.reduction[[i]]]], slot = 'ref.dims'),
    bridge = Misc( bridge.object[[bridge.reduction[[i]]]], slot = 'ref.dims')
   )
   all.dims <- list(
     object = 1:ncol(object.list[[i]][[object.reduction[[i]]]]),
     bridge = 1:ncol( bridge.object[[bridge.reduction[[i]] ]])
     )
   projected.dims.index <- which(sapply(ref.dims, function(x) !is.null(x)))
   if (length(projected.dims.index) == 0) {
     warning('No reference dims found in the dimensional reduction,',
             ' all dims in the dimensional reduction will be used.')
     if (all.dims[[1]] == all.dims[[2]]) {
       dims.list[[i]]  <- all.dims
     } else {
       stop( 'The number of dimensions in the object.list ',
             object.reduction[[i]],
             ' (', length(all.dims[[1]]), ') ',
       ' and the number of dimensions in the bridge object ',
       bridge.reduction[[i]],
       ' (', length(all.dims[[2]]), ') ',
       ' is different.')
     }
   } else {
     reference.dims.index <- setdiff(c(1:2), projected.dims.index)
     dims.list[[i]] <- list()
     dims.list[[i]][[reference.dims.index]] <- ref.dims[[projected.dims.index ]]
     dims.list[[i]][[projected.dims.index]] <- all.dims[[projected.dims.index]]
     names(dims.list[[i]]) <- c('object', 'bridge')
   }
    }
  object.list <- my.lapply(
    X = 1:length(x = object.list),
    FUN = function(x) {
      SA.inv <- ginv(
        X = Embeddings(
          object = bridge.object,
          reduction = bridge.reduction[[x]]
        )[ ,dims.list[[x]]$bridge]
      )
        if (!is.null(laplacian.reduction)) {
          lap.vector <- Embeddings(bridge.object[[laplacian.reduction]])[,laplacian.dims]
          X <- Embeddings(
            object = object.list[[x]],
            reduction = object.reduction[[x]]
          )[, dims.list[[x]]$object] %*% (SA.inv %*% lap.vector)
        } else {
          X <- Embeddings(
            object = object.list[[x]],
            reduction = object.reduction[[x]]
          )[,  dims.list[[x]]$object] %*% SA.inv
          colnames(X) <- Cells(bridge.object)
        }
      if (l2.norm) {
        X <- L2Norm(mat = X, MARGIN = 1)
      }
      colnames(x = X) <- paste0('bridge_',  colnames(x = X))
      suppressWarnings(
        object.list[[x]][[bridge.assay.name]] <- CreateAssayObject(data = t(X))
        )
      object.list[[x]][[bridge.assay.name]]@misc$SA.inv <- SA.inv
      DefaultAssay(object.list[[x]]) <- bridge.assay.name
      VariableFeatures(object = object.list[[x]]) <- rownames(object.list[[x]])
      return (object.list[[x]])
    }
  )
  if (!return.all.assays) {
    object.list <- my.lapply(
      X = object.list,
      FUN = function(x) {
        x <- DietSeurat(object = x, assays = bridge.assay.name, scale.data = TRUE)
        return(x)
      }
    )
  }
  if (single.object) {
    object.list <- object.list[[1]]
  }
  return(object.list)
}

#' Find bridge anchors between two unimodal datasets
#'
#' First, bridge object is used to reconstruct two single-modality profiles and
#' then project those cells into bridage graph laplacian space.
#' Next, find a set of anchors between two single-modality objects. These
#' anchors can later be used to integrate embeddings or transfer data from the reference to
#' query object using the \code{\link{MapQuery}} object.
#'
#'  \itemize{
#'  \item{ Bridge cells reconstruction
#'  }
#'   \item{ Find anchors between objects. It can be either IntegrationAnchors or TransferAnchor.
#'  }
#' }
#'
#' @inheritParams BridgeCellsRepresentation
#' @param anchor.type The type of anchors. Can
#' be one of:
#' \itemize{
#'   \item{Integration: Generate IntegrationAnchors for integration}
#'   \item{Transfer: Generate TransferAnchors for transfering data}
#' }
#' @param reference A vector specifying the object/s to be used as a reference
#' during integration or transfer data.
#' @param reduction Dimensional reduction to perform when finding anchors. Can
#' be one of:
#' \itemize{
#'   \item{cca: Canonical correlation analysis}
#'   \item{direct: Use assay data as a dimensional reduction}
#' }
#' @param reference.bridge.stored If refernece has stored the bridge dictionary representation
#' @param k.anchor How many neighbors (k) to use when picking anchors
#' @param k.score How many neighbors (k) to use when scoring anchors
#' @param verbose Print messages and progress
#' @param ... Additional parameters passed to \code{FindIntegrationAnchors} or
#' \code{FindTransferAnchors}
#'
#'
#' @return Returns an \code{\link{AnchorSet}} object that can be used as input to
#' \code{\link{IntegrateEmbeddings}}.or \code{\link{MapQuery}}
#'

FindBridgeAnchor <- function(object.list,
                             bridge.object,
                             object.reduction,
                             bridge.reduction,
                             anchor.type = c("Transfer", "Integration"),
                             reference = NULL,
                             laplacian.reduction = "lap",
                             laplacian.dims = 1:50,
                             reduction = c("direct", "cca"),
                             bridge.assay.name = "Bridge",
                             reference.bridge.stored = FALSE,
                             k.anchor = 20,
                             k.score = 50,
                             verbose = TRUE,
                             ...
                             ) {
  anchor.type <- match.arg(arg = anchor.type)
  reduction <- match.arg(arg = reduction)
  if (!is.null(laplacian.reduction)) {
    bridge.method <- "bridge graph"
  } else {
    bridge.method <- "bridge cells"
  }
  if (verbose) {
    message("Finding ", anchor.type," anchors")
    switch(
      EXPR = bridge.method,
      "bridge graph" = {
        message('Transform cells to bridge graph laplacian space')
      },
      "bridge cells" = {
        message('Transform cells to bridge cells space')
      }
    )
  }
  reference <- reference %||% c(1)
  query <- setdiff(c(1,2), reference)
  if (anchor.type == "Transfer") {
    stored.bridge.weights <- FALSE
    # check weight matrix
    if (is.null(bridge.object@tools$MapQuery)) {
      warning("No weights stored between reference and bridge obejcts.",
           "Please set store.weights to TRUE in MapQuery")
    } else if (is.null(object.list[[query]]@tools$MapQuery)) {
      warning("No weights stored between query and bridge obejcts.",
           "Please set store.weights to TRUE in MapQuery")
    } else {
      stored.bridge.weights <- TRUE
    }
  }
  if (reference.bridge.stored) {
    object.list[[query]] <- BridgeCellsRepresentation(
      object.list = object.list[[query]] ,
      bridge.object = bridge.object,
      object.reduction = object.reduction[[query]] ,
      bridge.reduction = bridge.reduction[[query]] ,
      bridge.assay.name = bridge.assay.name,
      laplacian.reduction = laplacian.reduction,
      laplacian.dims = laplacian.dims,
      verbose = verbose
    )
  } else {
    object.list <- BridgeCellsRepresentation(
      object.list = object.list ,
      bridge.object = bridge.object,
      object.reduction = object.reduction,
      bridge.reduction = bridge.reduction,
      bridge.assay.name = bridge.assay.name,
      laplacian.reduction = laplacian.reduction,
      laplacian.dims = laplacian.dims,
      verbose = verbose
    )
  }
  if (reduction == "direct") {
    anchor <- FindAssayAnchor(
      object.list = object.list ,
      reference = reference,
      slot = "data",
      anchor.type = anchor.type,
      assay = bridge.assay.name,
      k.anchor = k.anchor,
      k.score = k.score,
      verbose = verbose
    )
  } else if (reduction == "cca") {
    # set data slot to scale.data slot
    object.list <- lapply(
      X = object.list,
      FUN = function(x) {
     x <- SetAssayData(
       object = x,
       slot = "scale.data",
       new.data = as.matrix(
         x = GetAssayData(object = x, slot = "data")
         ))
     return(x)
     }
     )
    anchor <- switch(EXPR = anchor.type,
                     "Integration" = {
                       anchor <- FindIntegrationAnchors(
                         object.list = object.list,
                         k.filter = NA,
                         reference = reference,
                         reduction = "cca",
                         scale = FALSE,
                         k.anchor = k.anchor,
                         k.score = k.score,
                         verbose = verbose,
                         ...)
                       object.merge <- merge(x = object.list[[1]],
                                             y = object.list[2:length(object.list)]
                                             )
                       slot(
                         object = anchor,
                         name = "weight.reduction"
                         ) <- CreateDimReducObject(
                           embeddings = t(GetAssayData(
                             object = object.merge,
                             slot = 'data'
                           )),
                           key = "L_",
                           assay = bridge.assay.name
                         )
                       anchor
                     },
                     "Transfer" = {
                       anchor <-  FindTransferAnchors(
                         reference = object.list[[reference]],
                         query = object.list[[query]],
                         reduction = "cca",
                         scale = FALSE,
                         k.filter = NA,
                         k.anchor = k.anchor,
                         k.score = k.score,
                         verbose = verbose,
                         ...
                       )
                     }
    )
  }
  if (anchor.type == "Transfer") {
    if (stored.bridge.weights) {
      slot( object = anchor,name = "weight.reduction"
      )@misc$bridge.sets <- list(
        bridge.weights =   slot(object = bridge.object,
                                name = "tools"
        )$MapQuery_PrepareBridgeReference$weights.matrix,
        bridge.ref_anchor =  slot(object = bridge.object,
                                  name = "tools"
        )$MapQuery_PrepareBridgeReference$anchor[,1],
        query.weights =  slot(object = object.list[[query]],
                              name = "tools"
        )$MapQuery$weights.matrix,
        query.ref_anchor =  slot(object = object.list[[query]],
                                 name = "tools"
        )$MapQuery$anchor[,1]
      )
    }
  }
  slot(object = anchor, name = "command") <- LogSeuratCommand(
    object = object.list[[1]],
    return.command = TRUE
    )
  return(anchor)
}


# Helper function to transfer labels based on neighbors object
# @param nn.object  the query neighbors object
# @param reference.object the reference seurat object
# @param group.by  A vector of variables to group cells by
# @param weight.matrix A reference x query cell weight matrix
# @return Returns a list for predicted labels, prediction score and matrix
#' @importFrom Matrix sparseMatrix
#' @importFrom fastDummies dummy_cols
#' @importFrom Matrix rowMeans t
#'
TransferLablesNN <- function(
  nn.object = NULL,
  weight.matrix = NULL,
  reference.labels
){
  reference.labels.matrix <- CreateCategoryMatrix(labels = as.character(reference.labels))
  if (!is.null(x = weight.matrix) & !is.null(x = nn.object)) {
    warning('both nn.object and weight matrix are set. Only weight matrix is used for label transfer')
  }
  if (is.null(x = weight.matrix)) {
    select_nn <- Indices(nn.object)
    k.nn <- ncol(select_nn)
    j <- as.numeric(x = t(x = select_nn ))
    i <- ((1:length(x = j)) - 1) %/% k.nn + 1
    nn.matrix <- sparseMatrix(
      i = i,
      j = j,
      x = 1,
      dims = c(nrow(select_nn), nrow(reference.labels.matrix))
    )
    rownames(nn.matrix) <- Cells(nn.object)
  } else if (nrow(weight.matrix) == nrow(reference.labels.matrix)) {
    nn.matrix <- t(weight.matrix)
    k.nn <- 1
  } else if (ncol(weight.matrix) == nrow(reference.labels.matrix)) {
    nn.matrix <- weight.matrix
    k.nn <- 1
  } else {
    stop('wrong weights matrix input')
  }
  query.label.mat <- nn.matrix %*% reference.labels.matrix
  query.label.mat <- query.label.mat/k.nn
  prediction.max <- apply(X = query.label.mat, MARGIN = 1, FUN = which.max)

  query.label <- colnames(x = query.label.mat)[prediction.max]
  query.label.score <- apply(X = query.label.mat, MARGIN = 1, FUN = max)
  names(query.label) <- names(query.label.score) <- rownames(query.label.mat)
  if (is.factor(reference.labels)) {
    levels(query.label) <- levels(reference.labels)
  }
  output.list <- list(labels = query.label,
                      scores = query.label.score,
                      prediction.mat = query.label.mat
                      )
  return(output.list)
}

# transfer continuous value based on neighbors
#
TransferExpressionNN<- function(
  nn.object,
  reference.object,
  var.name = NULL
) {
  nn.matrix <- NNtoGraph(nn.object = nn.object,
                         col.cells = Cells(reference.object)
                         )
  reference.exp.matrix <- FetchData(object = reference.object, vars = var.name)
  # remove NA
  reference.exp.matrix <- reference.exp.matrix[complete.cases(reference.exp.matrix), ,drop= F]
  nn.matrix <- nn.matrix[, rownames(reference.exp.matrix)]

  # remove NO neighbor query
  nn.sum <- RowSumSparse(mat = nn.matrix)
  nn.matrix <- nn.matrix[nn.sum > 2, ]
  nn.sum <- nn.sum[nn.sum>2]

  # transfer data
  reference.exp.matrix <- as.matrix(reference.exp.matrix)
  query.exp.mat <- nn.matrix %*% reference.exp.matrix
  query.exp.mat <- sweep(x = query.exp.mat, MARGIN = 1, STATS = nn.sum, FUN = "/")

  # set output for all query cells
  query.exp.all <- data.frame(row.names = Cells(nn.object))
  query.exp.all[rownames(query.exp.mat),1] <- query.exp.mat[,1]
  colnames(query.exp.all) <- var.name
  return(query.exp.all)
}


#' @param reduction.name dimensional reduction name, lap by default
#' @param graph The name of graph
#' @rdname RunGraphLaplacian
#' @concept dimensional_reduction
#' @export
#' @method RunGraphLaplacian Seurat
#'
RunGraphLaplacian.Seurat <- function(
  object,
  graph,
  reduction.name = "lap",
  reduction.key ="LAP_",
  n = 50,
  verbose = TRUE,
  ...
) {
  lap_dir <- RunGraphLaplacian(object = object[[graph]],
                               n = n,
                               reduction.key = reduction.key ,
                               verbose = verbose,
                               ...
                               )
  object[[reduction.name]] <- lap_dir
  return(object)
}



#' @param n Total Number of Eigenvectors to compute and store (50 by default)
#' @param reduction.key dimensional reduction key, specifies the string before
#' the number for the dimension names. LAP by default
#' @param verbose Print message and process
#' @param ... Arguments passed to eigs_sym
#'
#'
#' @concept dimensional_reduction
#' @rdname RunGraphLaplacian
#' @export
#'
#' @importFrom Matrix diag t rowSums
#' @importFrom RSpectra eigs_sym
RunGraphLaplacian.default <- function(object,
                                      n = 50,
                                      reduction.key ="LAP_",
                                      verbose = TRUE,
                                      ...
) {
 if (!all(
   slot(object = t(x = object), name = "x") == slot(object = object, name = "x")
   )) {
   stop("Input graph is not symmetric")
 }
  if (verbose) {
    message("Generating normalized laplacian graph")
  }
  D_half <- sqrt(x = rowSums(x = object))
  L <- -1 * (t(object / D_half) / D_half)
  diag(L) <- 1 + diag(L)
  if (verbose) {
    message("Performing eigendecomposition of the normalized laplacian graph")
  }
  L_eigen <- eigs_sym(L, k = n + 1, which = "SM", ...)
  #delete the first eigen vector
  new_order <- n:1
  lap_output <- list(eigen_vector = Re(L_eigen$vectors[, new_order]),
                         eigen_value = L_eigen$values[new_order]
  )
  rownames(lap_output$eigen_vector) <- colnames(object)
  colnames(lap_output$eigen_vector) <- paste0(reduction.key, 1:n )
  lap_dir <- CreateDimReducObject(embeddings = lap_output$eigen_vector,
                       key = reduction.key,
                       assay = DefaultAssay(object),
                       stdev = lap_output$eigen_value
  )
  return(lap_dir)
}


# Check if the var.name already existed in the meta.data
#
CheckMetaVarName <- function(object, var.name) {
  if (var.name %in% colnames(x = object[[]])) {
    var.name.exist <- var.name
    var.name <- rev(
      x = make.unique(
        names = c(colnames(object[[]]), var.name.exist)
        )
      )[1]
    warning(var.name.exist, " is already existed in the meta.data. ",
            var.name, " will store leverage score value")
  }
  return(var.name)
}



# Run hnsw to find neighbors
#
# @param data Data to build the index with
# @param query A set of data to be queried against data
# @param metric Distance metric; can be one of "euclidean", "cosine", "manhattan",
# "hamming"
# @param k Number of neighbors
# @param ef_construction  A larger value means a better quality index, but increases build time.
# @param ef Higher values lead to improved recall at the expense of longer search time.
# @param n_threads Maximum number of threads to use.
# @param index optional index object, will be recomputed if not provided
#' @importFrom RcppHNSW hnsw_build hnsw_search
#
HnswNN <- function(data,
                    query = data,
                    metric = "euclidean",
                    k,
                    ef_construction = 200,
                    ef = 10,
                    index = NULL,
                    n_threads = 0
) {
  idx <- index %||% hnsw_build(
    X = data,
    distance = metric,
    ef = ef_construction,
    n_threads = n_threads
    )
  nn <- hnsw_search(
    X = query,
    ann = idx,
    k = k,
    ef = ef,
    n_threads = n_threads
    )
  names(nn) <- c("nn.idx", "nn.dists")
  nn$idx <- idx
  nn$alg.info <- list(metric = metric, ndim = ncol(x = data))
  return(nn)
}


# Calculate reference index from the integrated object
#
IntegrationReferenceIndex <- function(object) {
  if (is.null(object@tools$Integration@sample.tree)) {
    reference.index <- object@commands$FindIntegrationAnchors$reference
    if (length(x = reference.index) > 1) {
      stop('the number of the reference is bigger than 1')
    }
  } else {
    reference.index <- SampleIntegrationOrder(tree = object@tools$Integration@sample.tree)[1]
  }
  return(reference.index)
}


# Calculate mean and sd
#
SparseMeanSd <- function(object,
                         assay = NULL,
                         slot = 'data',
                         features = NULL,
                         eps = 1e-8
){
  assay <- assay%||% DefaultAssay(object)
  features <- features %||% rownames(object[[assay]])
  assay <- assay %||% DefaultAssay(object = object)
  mat <- GetAssayData(object = object[[assay]], slot = slot)[features,]
  if (class(mat)[1] !='dgCMatrix'){
    stop('Matrix is not sparse')
  }
  mat.mean <-  RowMeanSparse(mat)
  mat.sd <-  sqrt(RowVarSparse(mat))
  names(mat.mean) <- names(mat.sd) <- rownames(mat)
  mat.sd <- MinMax(data = mat.sd, min = eps, max = max(mat.sd))
  output <- list(mean = mat.mean, sd = mat.sd)
  return(output)
}



# Run PCA on sparse matrix
#
#' @importFrom Matrix t
#' @importFrom rlang exec
#' @importFrom irlba irlba
#
#
RunPCA_Sparse <- function(
  object,
  features = NULL,
  reduction.key = "PCsp_",
  reduction.name = "pca.sparse",
  npcs = 50,
  do.scale = TRUE,
  verbose = TRUE
) {
  features <- features %||% VariableFeatures(object)
  data <- GetAssayData(object = object, slot = "data")[features,]
  n <- npcs
  args <- list(A = t(data), nv = n)
  args$center <- RowMeanSparse(data)
  feature.var <- RowVarSparse(data)
  args$totalvar <- sum(feature.var)
  if (do.scale) {
    args$scale <- sqrt(feature.var)
    args$scale <- MinMax(args$scale, min = 1e-8, max = max(args$scale))
  } else {
    args$scale <- FALSE
  }
  if (verbose) {
    message("Running PCA")
  }
  pca.irlba <- exec(.fn = irlba, !!!args)
  sdev <- pca.irlba$d/sqrt(max(1, ncol(data) - 1))
  feture.loadings <- pca.irlba$v
  rownames(feture.loadings) <- rownames(data)
  embeddings <- sweep(x = pca.irlba$u, MARGIN = 2, STATS = pca.irlba$d, FUN = "*")
  rownames(embeddings) <- colnames(data)
  colnames(feture.loadings) <- colnames(embeddings) <- paste0(reduction.key, 1:npcs)
  object[[reduction.name]] <- CreateDimReducObject(
    embeddings = embeddings,
    loadings = feture.loadings,
    stdev = sdev,
    key = reduction.key,
    assay = DefaultAssay(object),
    misc = list(d = pca.irlba$d)
  )
  return(object)
}

# Smoothing labels based on the clusters
# @param labels the original labels
# @param clusters the clusters that are used to smooth labels
#
SmoothLabels <- function(labels, clusters) {
  cluster.set <- unique(clusters)
  smooth.labels <- labels
  for (c in cluster.set) {
    cell.c <- which(clusters == c)
    smooth.labels[cell.c] <- names(sort(table(labels[cell.c]), decreasing = T)[1])
  }
  return(smooth.labels)
}



#' Project query data to reference dimensional reduction
#'
#' @param query Query object
#' @param reference Reference object
#' @param mode Projection mode name for projection
#'  \itemize{
#' \item{pcaproject: PCA projection}
#' \item{lsiproject: LSI projection}
#' }
#' @param reference.reduction Name of dimensional reduction in the reference object
#' @param combine Determine if query and reference objects are combined
#' @param query.assay Assay used for query object
#' @param reference.assay Assay used for reference object
#' @param features Features used for projection
#' @param do.scale Determine if scale expression matrix in the pcaproject mode
#' @param reduction.name dimensional reduction name, reference.reduction is used by default
#' @param reduction.key dimensional reduction key, the key in reference.reduction
#' is used by default
#' @param verbose Print progress and message
#'
#' @return Returns a query-only or query-reference combined seurat object
#' @export
ProjectDimReduc <- function(query,
                            reference,
                            mode = c('pcaproject', 'lsiproject'),
                            reference.reduction,
                            combine = FALSE,
                            query.assay = NULL,
                            reference.assay = NULL,
                            features = NULL,
                            do.scale = TRUE,
                            reduction.name = NULL,
                            reduction.key= NULL,
                            verbose = TRUE
) {
  query.assay <- query.assay %||% DefaultAssay(object = query)
  reference.assay <- reference.assay %||% DefaultAssay(object = reference)
  DefaultAssay(object = query) <- query.assay
  DefaultAssay(object = reference) <- reference.assay
  reduction.name <- reduction.name %||% reference.reduction
  reduction.key <- reduction.key %||% Key(object = reference[[reference.reduction]])
  if (reduction.name %in% Reductions(object = query)) {
    warning(reduction.name,
            ' already exists in the query object. It will be overwritten.'
    )
  }
  features <- features %||% rownames(x = Loadings(object = reference[[reference.reduction]]))
  features <- intersect(x = features, y = rownames(x = query))
  if (mode == 'lsiproject') {
    if (verbose) {
      message('LSI projection to ', reference.reduction)
    }
    projected.embeddings <- ProjectSVD(
      reduction = reference[[reference.reduction]],
      data = GetAssayData(object = query, assay = query.assay, slot = "data"),
      mode = "lsi",
      do.center = FALSE,
      do.scale = FALSE,
      features = features,
      use.original.stats = FALSE,
      verbose = verbose
    )
  } else if (mode == 'pcaproject') {
    if (inherits(query[[query.assay]], what = 'SCTAssay')) {
      if (verbose) {
        message('PCA projection to ', reference.reduction, ' in SCT assay')
      }
      query <- suppressWarnings(
        expr = GetResidual(object = query,
                           assay = query.assay,
                           features = features,
                           verbose = FALSE)
      )
      query.mat <- GetAssayData(object = query, slot = 'scale.data')[features,]

      projected.embeddings <- t(
        crossprod(x = Loadings(
          object = reference[[reference.reduction]])[features, ],
          y = query.mat
        )
      )
    } else {
      if (verbose) {
        message('PCA projection to ', reference.reduction)
      }
      projected.embeddings <- ProjectCellEmbeddings(
        reference = reference,
        reduction = reference.reduction,
        query = query,
        scale = do.scale,
        dims = 1:ncol(reference[[reference.reduction]]),
        feature.mean = NULL,
        verbose = verbose
      )
    }
  }
  query[[reduction.name]] <- CreateDimReducObject(
    embeddings = projected.embeddings,
    loadings = Loadings(reference[[reference.reduction]])[features,],
    assay = query.assay,
    key = reduction.key,
    misc = Misc(reference[[reference.reduction]])
  )
  if (combine) {
    query <- DietSeurat(object = query,
                        dimreducs = reduction.name,
                        features = features,
                        assays = query.assay
    )
    reference <- DietSeurat(object = reference,
                            dimreducs = reference.reduction,
                            features = features,
                            assays = reference.assay)
    suppressWarnings(
      combine.obj <- merge(query, reference,
                           merge.dr = c(reduction.name, reference.reduction)
      )
    )
    Idents(combine.obj) <- c(rep(x = 'query', times = ncol(query)),
                            rep(x = 'reference', times = ncol(reference))
                            )
    return(combine.obj)
  } else {
    return(query)
  }
}


#' Prepare the bridge and reference datasets
#'
#' Preprocess the multi-omic bridge and unimodal reference datasets into
#' an extended reference.
#' This function performs the following three steps:
#' 1. Performs within-modality harmonization between bridge and reference
#' 2. Performs dimensional reduction on the SNN graph of bridge datasets via
#' Laplacian Eigendecomposition
#' 3. Constructs a bridge dictionary representation for unimodal reference cells
#'
#' @param reference A reference Seurat object
#' @param bridge A multi-omic bridge Seurat object
#' @param reference.reduction Name of dimensional reduction of the reference object (default is 'pca')
#' @param reference.dims Number of dimensions used for the reference.reduction (default is 50)
#' @param normalization.method Name of normalization method used: LogNormalize
#' or SCT
#' @param reference.assay Assay name for reference (default is \code{\link{DefaultAssay}})
#' @param bridge.ref.assay Assay name for bridge used for reference mapping. RNA by default
#' @param bridge.query.assay Assay name for bridge used for query mapping. ATAC by default
#' @param supervised.reduction Type of supervised dimensional reduction to be performed
#' for integrating the bridge and query.
#' #' Options are:
#' \itemize{
#'    \item{slsi: Perform supervised LSI as the dimensional reduction for
#'    the bridge-query integration}
#'    \item{spca: Perform supervised PCA as the dimensional reduction for
#'    the bridge-query integration}
#'    \item{NULL: no supervised dimensional reduction will be calculated.
#'    bridge.query.reduction is used for the bridge-query integration}
#' }
#' @param bridge.query.reduction Name of dimensions used for the bridge-query harmonization.
#' 'bridge.query.reduction' and 'supervised.reduction' cannot be NULL together.
#' @param bridge.query.features Features used for bridge query dimensional reduction
#' (default is NULL which uses VariableFeatures from the bridge object)
#' @param laplacian.reduction.name Name of dimensional reduction name of graph laplacian eigenspace (default is 'lap')
#' @param laplacian.reduction.key Dimensional reduction key (default is 'lap_')
#' @param laplacian.reduction.dims Number of dimensions used for graph laplacian eigenspace (default is 50)
#' @param verbose Print progress and message (default is TRUE)
#'
#' @export
#' @return Returns a \code{BridgeReferenceSet} that can be used as input to
#'  \code{\link{FindBridgeTransferAnchors}}.
#' The parameters used are stored in the \code{BridgeReferenceSet} as well
#'
PrepareBridgeReference <- function (
  reference,
  bridge,
  reference.reduction = 'pca',
  reference.dims = 1:50,
  normalization.method = c('SCT', 'LogNormalize'),
  reference.assay = NULL,
  bridge.ref.assay = 'RNA',
  bridge.query.assay = 'ATAC',
  supervised.reduction = c('slsi', 'spca', NULL),
  bridge.query.reduction = NULL,
  bridge.query.features = NULL,
  laplacian.reduction.name = 'lap',
  laplacian.reduction.key = 'lap_',
  laplacian.reduction.dims = 1:50,
  verbose = TRUE
) {
  ## checking
  if (!is.null(supervised.reduction)) {
  supervised.reduction <- match.arg(arg = supervised.reduction)
  }
  if (!is.null(x = bridge.query.reduction) & !is.null(x = supervised.reduction)) {
    stop('bridge.query.reduction and supervised.reduction can only set one.',
         'If you want to set bridge.query.reduction, supervised.reduction should set to NULL')
  }
  if (is.null(x = bridge.query.reduction) & is.null(x = supervised.reduction)) {
    stop('Both bridge.query.reduction and supervised.reduction are NULL. One of them needs to be set')
  }
  bridge.query.features <- bridge.query.features %||%
    VariableFeatures(object = bridge[[bridge.query.assay]])
  if (length(x = bridge.query.features) == 0) {
    stop('bridge object ', bridge.query.assay,
         ' assay has no variable genes and bridge.query.features has no input')
  }
  # modality harmonization
  reference.assay <- reference.assay %||% DefaultAssay(reference)
  DefaultAssay(reference) <- reference.assay
  DefaultAssay(bridge) <- bridge.ref.assay
  ref.anchor  <- FindTransferAnchors(
    reference =  reference,
    reference.reduction = reference.reduction,
    normalization.method = normalization.method,
    dims = reference.dims,
    query = bridge,
    recompute.residuals = TRUE,
    features = rownames(reference[[reference.reduction]]@feature.loadings),
    k.filter = NA,
    verbose = verbose
  )
  bridge <- MapQuery(anchorset = ref.anchor,
                     reference = reference,
                     query = bridge,
                     store.weights = TRUE,
                     verbose = verbose
  )
  bridge.ref.reduction <- paste0('ref.', reference.reduction)
  bridge <- FindNeighbors(object = bridge,
                          reduction = bridge.ref.reduction,
                          dims = 1:ncol(x = bridge[[bridge.ref.reduction]]),
                          return.neighbor = FALSE,
                          graph.name = c('bridge.ref.nn', 'bridge.ref.snn'),
                          prune.SNN = 0)
  bridge <- RunGraphLaplacian(object = bridge,
                              graph = "bridge.ref.snn",
                              reduction.name = laplacian.reduction.name,
                              reduction.key = laplacian.reduction.key,
                              verbose = verbose)
  DefaultAssay(object = bridge) <- bridge.query.assay
  if (!is.null(supervised.reduction)) {
    bridge <- switch(EXPR = supervised.reduction,
                     'slsi' = {
                       bridge.reduc <- RunSLSI(object = bridge,
                                               features = VariableFeatures(bridge),
                                               graph = 'bridge.ref.nn',
                                               assay = bridge.query.assay
                       )
                       bridge.reduc
                     },
                     'spca' = {
                       bridge.reduc <- RunSPCA(object = bridge,
                                               features = VariableFeatures(bridge),
                                               graph = 'bridge.ref.snn',
                                               assay = bridge.query.assay
                       )
                       bridge.reduc
                     }
    )
  }
  # bridge representation
  reference.bridge <- BridgeCellsRepresentation(
    object.list =  reference,
    bridge.object = bridge,
    object.reduction = c(reference.reduction),
    bridge.reduction =  c(bridge.ref.reduction),
    laplacian.reduction = laplacian.reduction.name,
    laplacian.dims = laplacian.reduction.dims
  )
  reference[['Bridge']] <- reference.bridge[['Bridge']]
  reference <- merge(x = reference, y = bridge, merge.dr = NA)
  reference@tools$MapQuery_PrepareBridgeReference <- bridge@tools$MapQuery
  command <- LogSeuratCommand(object = reference, return.command = TRUE)
  slot(object = command, name = "params")$bridge.query.features <- NULL
  command.name <- slot(object = command, name = "name")
  reference[[command.name]] <- command
  return(reference)
}


#' Find bridge anchors between query and extended bridge-reference
#'
#' Find a set of anchors between unimodal query and the other unimodal reference
#' using a pre-computed \code{\link{BridgeReferenceSet}}.
#' This function performs three steps:
#' 1. Harmonize the bridge and query cells in the bridge query reduction space
#' 2. Construct the bridge dictionary representations for query cells
#' 3. Find a set of anchors between query and reference in the bridge graph laplacian eigenspace
#' These anchors can later be used to integrate embeddings or transfer data from the reference to
#' query object using the \code{\link{MapQuery}} object.

#' @param extended.reference BridgeReferenceSet object generated from
#'  \code{\link{PrepareBridgeReference}}
#' @param query A query Seurat object
#' @param query.assay Assay name for query-bridge integration
#' @param scale Determine if scale the query data for projection
#' @param dims Number of dimensions for query-bridge integration
#' @param reduction Dimensional reduction to perform when finding anchors.
#' Options are:
#' \itemize{
#'    \item{pcaproject: Project the PCA from the bridge onto the query. We
#'    recommend using PCA when bridge and query datasets are from scRNA-seq}
#'    \item{lsiproject: Project the LSI from the bridge onto the query. We
#'    recommend using LSI when bridge and query datasets are from scATAC-seq or scCUT&TAG data.
#'    This requires that LSI or supervised LSI has been computed for the bridge dataset, and the
#'    same features (eg, peaks or genome bins) are present in both the bridge
#'    and query.
#' }
#' }
#' @param verbose Print messages and progress
#'
#' @export
#' @return Returns an \code{AnchorSet} object that can be used as input to
#' \code{\link{TransferData}}, \code{\link{IntegrateEmbeddings}} and
#' \code{\link{MapQuery}}.
#'
FindBridgeTransferAnchors <- function(
  extended.reference,
  query,
  query.assay = NULL,
  dims = 1:30,
  scale = FALSE,
  reduction = c('lsiproject', 'pcaproject'),
  bridge.reduction = c('direct', 'cca'),
  verbose = TRUE
) {
  bridge.reduction <- match.arg(arg = bridge.reduction)
  reduction <-  match.arg(arg = reduction)
  query.assay <- query.assay %||% DefaultAssay(query)
  DefaultAssay(query) <- query.assay
  params <- Command(object = extended.reference, command = 'PrepareBridgeReference')
  bridge.query.assay <- params$bridge.query.assay
  bridge.query.reduction <- params$bridge.query.reduction %||% params$supervised.reduction
  reference.reduction <- params$reference.reduction
  bridge.ref.reduction <-  paste0('ref.', reference.reduction)
  DefaultAssay(extended.reference) <- bridge.query.assay
  extended.reference.bridge <- DietSeurat(
    object = extended.reference,
    assays = bridge.query.assay,
    dimreducs = c(bridge.ref.reduction, bridge.query.reduction, params$laplacian.reduction.name)
    )
    query.anchor <- FindTransferAnchors(
      reference = extended.reference.bridge,
      reference.reduction = bridge.query.reduction,
      dims = dims,
      query = query,
      reduction = reduction,
      scale = scale,
      features = rownames(Loadings(extended.reference[[bridge.query.reduction]])),
      k.filter = NA,
      verbose = verbose
    )

  query <- MapQuery(anchorset =  query.anchor,
                    reference = extended.reference.bridge,
                    query = query,
                    store.weights = TRUE
  )
  DefaultAssay(extended.reference) <- 'Bridge'
  bridge_anchor  <- FindBridgeAnchor(
    object.list = list(DietSeurat(object = extended.reference, assays = 'Bridge'), query),
    bridge.object = extended.reference.bridge,
    object.reduction = c(reference.reduction, paste0('ref.', bridge.query.reduction)),
    bridge.reduction = c(bridge.ref.reduction, bridge.query.reduction),
    anchor.type = "Transfer",
    reduction = bridge.reduction,
    reference.bridge.stored = TRUE,
    verbose = verbose
  )
  return(bridge_anchor)
}



#' Find integration bridge anchors between query and extended bridge-reference
#'
#' Find a set of anchors between unimodal query and the other unimodal reference
#' using a pre-computed \code{\link{BridgeReferenceSet}}.
#' These integration anchors can later be used to integrate query and reference
#' using the \code{\link{IntegrateEmbeddings}} object.
#'
#' @inheritParams FindBridgeTransferAnchors
#' @param integration.reduction Dimensional reduction to perform when finding anchors
#' between query and reference.
#' Options are:
#' \itemize{
#'    \item{direct: find anchors directly on the bridge representation space}
#'    \item{cca: perform cca on the on the bridge representation space and then find anchors
#' }
#' }
#'
#' @export
#' @return Returns an \code{AnchorSet} object that can be used as input to
#' \code{\link{IntegrateEmbeddings}}.
#'
FindBridgeIntegrationAnchors <- function(
  extended.reference,
  query,
  query.assay = NULL,
  dims = 1:30,
  scale = FALSE,
  reduction = c('lsiproject', 'pcaproject'),
  integration.reduction = c('direct', 'cca'),
  verbose = TRUE
) {
  reduction <-  match.arg(arg = reduction)
  integration.reduction <-  match.arg(arg = integration.reduction)
  query.assay <- query.assay %||% DefaultAssay(query)
  DefaultAssay(query) <- query.assay

  params <- Command(object = extended.reference, command = 'PrepareBridgeReference')
  bridge.query.assay <- params$bridge.query.assay
  bridge.query.reduction <- params$bridge.query.reduction %||% params$supervised.reduction
  reference.reduction <- params$reference.reduction
  bridge.ref.reduction <- paste0( 'ref.', params$bridge.ref.reduction)
  DefaultAssay(extended.reference) <- bridge.query.assay

  extended.reference.bridge <- DietSeurat(
    object = extended.reference,
    assays = bridge.query.assay,
    dimreducs = c(bridge.query.reduction, bridge.ref.reduction, params$laplacian.reduction.name)
    )

  query.anchor <- FindTransferAnchors(
    reference = extended.reference.bridge,
    reference.reduction = bridge.query.reduction,
    dims = dims,
    query = query,
    reduction = reduction,
    scale = scale,
    features = rownames(Loadings(extended.reference.bridge[[bridge.query.reduction]])),
    k.filter = NA,
    verbose = verbose
  )
  query <- MapQuery(anchorset =  query.anchor,
                    reference = extended.reference.bridge,
                    query = query,
                    store.weights = TRUE
  )
  DefaultAssay(extended.reference) <- 'Bridge'
  bridge_anchor  <- FindBridgeAnchor(
    object.list = list(DietSeurat(object = extended.reference, assays = 'Bridge'), query),
    bridge.object = extended.reference.bridge,
    reduction = integration.reduction,
    object.reduction = c(reference.reduction, paste0('ref.', bridge.query.reduction)),
    bridge.reduction = c(bridge.ref.reduction, bridge.query.reduction),
    anchor.type = "Integration",
    reference.bridge.stored = TRUE,
    verbose = verbose
  )
  return(bridge_anchor)
}


#' Perform integration on the joint PCA cell embeddings.
#'
#' This is a convenience wrapper function around the following three functions
#' that are often run together when perform integration.
#' #' \code{\link{FindIntegrationAnchors}}, \code{\link{RunPCA}},
#' \code{\link{IntegrateEmbeddings}}.
#'
#' @inheritParams FindIntegrationAnchors
#' @param new.reduction.name Name of integrated dimensional reduction
#' @param npcs Total Number of PCs to compute and store (50 by default)
#' @param findintegrationanchors.args A named list of additional arguments to
#' \code{\link{FindIntegrationAnchors}}
#' @param verbose Print messages and progress
#'
#' @importFrom rlang invoke
#' @return Returns a Seurat object with integrated dimensional reduction
#' @export
#'
FastRPCAIntegration <- function(
  object.list,
  reference = NULL,
  anchor.features = 2000,
  k.anchor = 20,
  dims = 1:30,
  scale = TRUE,
  normalization.method = c("LogNormalize", "SCT"),
  new.reduction.name = 'integrated_dr',
  npcs = 50,
  findintegrationanchors.args = list(),
  verbose = TRUE
) {
  npcs <- max(npcs, dims)
  my.lapply <- ifelse(
    test = verbose && nbrOfWorkers() == 1,
    yes = pblapply,
    no = future_lapply
  )
  reduction <- 'rpca'
  if (is.numeric(x = anchor.features)) {
    anchor.features <- SelectIntegrationFeatures(
      object.list = object.list,
      nfeatures = anchor.features,
      verbose = FALSE
    )
  }
  if (normalization.method == 'SCT') {
    scale <- FALSE
    object.list <- PrepSCTIntegration(object.list = object.list,
                                      anchor.features = anchor.features
    )
  }

  if (verbose) {
    message('Performing PCA for each object')
  }
  object.list <- my.lapply(X = object.list,
                           FUN = function(x) {
                             if (normalization.method != 'SCT') {
                               x <- ScaleData(x, features = anchor.features, do.scale = scale, verbose = FALSE)
                             }
                             x <- RunPCA(x, features = anchor.features, verbose = FALSE, npcs = npcs)
                             return(x)
                           }
  )

  anchor <- invoke(
    .fn = FindIntegrationAnchors,
    .args = c(list(
      object.list = object.list,
      reference = reference,
      anchor.features = anchor.features,
      reduction = reduction,
      normalization.method = normalization.method,
      scale = scale,
      k.anchor = k.anchor,
      dims = dims,
      verbose = verbose
    ), findintegrationanchors.args
    )
  )
  object_merged <- merge(x = object.list[[1]],
                         y = object.list[2:length(object.list)]

  )

  anchor.feature <- slot(object = anchor, name = 'anchor.features')
  if (normalization.method != 'SCT') {
    object_merged <- ScaleData(object = object_merged,
                               features = anchor.feature,
                               do.scale = scale,
                               verbose = FALSE
    )
  }
  object_merged <- RunPCA(object_merged,
                          features = anchor.feature,
                          verbose = FALSE,
                          npcs = npcs

  )

  temp <- object_merged[["pca"]]
  object_merged <- IntegrateEmbeddings(
    anchorset = anchor,
    reductions = object_merged[['pca']],
    new.reduction.name = new.reduction.name,
    verbose = verbose)
  object_merged[['pca']] <- temp
  VariableFeatures(object = object_merged) <- anchor.feature
  return(object_merged)

}


#' Transfer embeddings from sketched cells to the full data
#'
#' @importFrom MASS ginv
#' @importFrom Matrix t
#'
#' @export
#'
UnSketchEmbeddings <- function(atom.data,
                               atom.cells = NULL,
                               orig.data,
                               embeddings,
                               sketch.matrix = NULL
) {
  if(!all(rownames(atom.data) == rownames(orig.data))) {
    stop('fetures in atom.data and orig.data are not identical')
  } else {
    features = rownames(atom.data)
  }
  atom.cells <- atom.cells %||% colnames(x = atom.data)
  if (inherits(x = orig.data, what = 'DelayedMatrix') ) {
    matrix.prod.function <- crossprod_DelayedAssay
  } else if(inherits(x = orig.data, what = 'IterableMatrix')) {
    matrix.prod.function <- crossprod_BPCells
  } else {
    matrix.prod.function <- crossprod
  }
  sketch.matrix <- sketch.matrix %||% as.sparse(diag(length(features)))
  atom.data <- atom.data[, atom.cells]
  embeddings <- embeddings[atom.cells,]
  exp.mat <- as.matrix(x = t(x = atom.data) %*% sketch.matrix)
  sketch.transform <- ginv(X = exp.mat) %*% embeddings
  emb <- matrix.prod.function(
    x = as.matrix(sketch.matrix %*% sketch.transform),
    y = orig.data
  )
  emb <- as.matrix(x = emb)
  return(emb)
}


FeatureSketch <- function(features, ratio = 0.8, seed = 123) {
  sketch.R <- t(x = CountSketch(
    nsketch = round(x = ratio *  length(x = features)),
    ncells = length(x = features),
    seed = seed)
  )
  return(sketch.R)
}


