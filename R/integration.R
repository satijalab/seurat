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
  reduction = c("cca", "rpca", "rlsi"),
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
  if (nn.reduction %in% c("pca", "lsi")) {
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
  if (normalization.method == "SCT") {
      # ensure all residuals required are computed
      query <- suppressWarnings(expr = GetResidual(object = query, assay = query.assay, features = features, verbose = FALSE))
      if (is.null(x = reference.reduction)) {
        reference <- suppressWarnings(expr = GetResidual(object = reference, assay = reference.assay, features = features, verbose = FALSE))
        features <- intersect(
          x = features,
          y = intersect(
            x = rownames(x = GetAssayData(object = query[[query.assay]], slot = "scale.data")),
            y = rownames(x = GetAssayData(object = reference[[reference.assay]], slot = "scale.data"))
          )
        )
        reference[[reference.assay]] <- as(
          object = CreateAssayObject(
            data = GetAssayData(object = reference[[reference.assay]], slot = "scale.data")[features, ]),
          Class = "SCTAssay"
        )
        reference <- SetAssayData(
          object = reference,
          slot = "scale.data",
          assay = reference.assay,
          new.data =  as.matrix(x = GetAssayData(object = reference[[reference.assay]], slot = "data"))
      )
    }
    query[[query.assay]] <- as(
      object = CreateAssayObject(
        data = GetAssayData(object = query[[query.assay]], slot = "scale.data")[features, ]),
      Class = "SCTAssay"
    )
    query <- SetAssayData(
      object = query,
      slot = "scale.data",
      assay = query.assay,
      new.data = as.matrix(x = GetAssayData(object = query[[query.assay]], slot = "data"))
    )
    feature.mean <- "SCT"
  }
  # only keep necessary info from objects
  query <- DietSeurat(
    object = query,
    assays = query.assay,
    dimreducs = reference.reduction,
    features = features,
    scale.data = TRUE
  )
  # check assay in the reference.reduction
  if (!is.null(reference.reduction) &&
    slot(object = reference[[reference.reduction]], name = "assay.used") != reference.assay) {
    warnings("reference assay is diffrent from the assay.used in", reference.reduction)
    slot(object = reference[[reference.reduction]], name = "assay.used") <- reference.assay
  }

  reference <- DietSeurat(
    object = reference,
    assays = reference.assay,
    dimreducs = reference.reduction,
    features = features,
    scale.data = TRUE
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
      projected.pca <- ProjectCellEmbeddings(
        reference = reference,
        reduction = reference.reduction,
        query = query,
        scale = scale,
        dims = dims,
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
        data = GetAssayData(object = query, assay = query.assay, slot = "data"),
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
      data = Embeddings(object = combined.ob[[reduction]])[Cells(x = query), ],
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
    nn.idx1 <- Index(object = reference[[reference.neighbors]])
  }
  anchors <- FindAnchors(
    object.pair = combined.ob,
    assay = c(reference.assay, query.assay),
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
  for (i in unique(x = c(reference.assay, query.assay))) {
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
  cell.names.map <- Cells(x = unintegrated)
  names(x = cell.names.map) <- make.unique(names = unname(obj = do.call(
    what = c,
    args = lapply(X = object.list, FUN = Cells)))
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
      embeddings = as.matrix(x = t(GetAssayData(reference.integrated[[new.reduction.name.safe]]))),
      assay = intdr.assay,
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
      slot = 'data',
      check.matrix = FALSE
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
  unintegrated[[new.reduction.name]] <- CreateDimReducObject(
    embeddings = as.matrix(x = t(x = integrated.data)),
    assay = intdr.assay,
    key = paste0(new.reduction.name.safe, "_")
  )
  unintegrated <- SetIntegrationData(
    object = unintegrated,
    integration.name = "Integration",
    slot = "anchors",
    new.data = anchors
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
#' @param transferdata.args A named list of additional arguments to
#' \code{\link{TransferData}}
#' @param integrateembeddings.args A named list of additional arguments to
#' \code{\link{IntegrateEmbeddings}}
#' @param projectumap.args A named list of additional arguments to
#' \code{\link{ProjectUMAP}}
#' @return Returns a modified query Seurat object containing:
#'
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
  reduction.model = NULL,
  transferdata.args = list(),
  integrateembeddings.args = list(),
  projectumap.args = list(),
  verbose = TRUE
) {

  # determine anchor type
  if (grepl(pattern = "pca", x = slot(object = anchorset, name = "command")$reduction)) {
    anchor.reduction <- "pcaproject"
    # check if the anchorset can be used for mapping
    if (is.null(x = slot(object = anchorset, name = "command")$reference.reduction)) {
      stop('The reference.reduction parameter was not set when running ',
      'FindTransferAnchors, so the resulting AnchorSet object cannot be used ',
      'in the MapQuery function.')
    }
  } else if (grepl(pattern = "cca", x = slot(object = anchorset, name = "command")$reduction)) {
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
  } else if (grepl(pattern = "lsi", x = slot(object = anchorset, name = "command")$reduction))  {
    anchor.reduction <- "lsiproject"
  } else {
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
  if (!is.null(x = refdata)) {
    query <- invoke(
      .fn = TransferData,
      .args  = c(list(
        anchorset = anchorset,
        reference = reference,
        query = query,
        refdata = refdata,
        store.weights = TRUE,
        verbose = verbose
        ), transferdata.args
      )
    )
    if (transferdata.args$weight.reduction == integrateembeddings.args$weight.reduction) {
      reuse.weights.matrix <- TRUE
    }
  }
  if (anchor.reduction != "cca"){
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
#' @export
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
      prediction.max <- apply(X = prediction.scores, MARGIN = 1, FUN = max)
      if (is.null(x = query)){
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
        Key(object = predictions) <- paste0("predictionscore", rd.name, "_")
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
      Key(object = new.assay) <- paste0(rd.name, "_")
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
      index = nn.idx1
    )
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
      index = nn.idx1
    )
  }
  if (length(x = reduction.2) > 0) {
    nnab <- NNHelper(
      data = Embeddings(object = object[[reduction.2]])[cells2, nn.dims],
      query = Embeddings(object = object[[reduction.2]])[cells1, nn.dims],
      k = k,
      method = nn.method,
      n.trees = n.trees,
      eps = eps,
      index = nn.idx2
    )
    nnba <- NNHelper(
      data = Embeddings(object = object[[reduction]])[cells1, nn.dims],
      query = Embeddings(object = object[[reduction]])[cells2, nn.dims],
      k = k,
      method = nn.method,
      n.trees = n.trees,
      eps = eps,
      index = nn.idx1
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
      index = nn.idx2
    )
    nnba <- NNHelper(
      data = dims.cells1.opposite,
      query = dims.cells2.opposite,
      k = k,
      method = nn.method,
      n.trees = n.trees,
      eps = eps,
      index = nn.idx1
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
    misc = NULL
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

# Rescale query with mean and sd from reference, or known mean and SD
#
# @param reference A reference object
# @param query A query object
# @param features Features to scale
# @param scale Scale data (divide by SD)
# @return Returns a matrix containing the scaled query data
RescaleQuery <- function(
  reference,
  query,
  reference.assay = NULL,
  query.assay = NULL,
  features = NULL,
  feature.mean = NULL,
  feature.sd = NULL,
  scale = TRUE
) {
  reference.assay <- reference.assay %||% DefaultAssay(object = reference)
  query.assay <- query.assay %||% DefaultAssay(object = query)
  features <- features %||% intersect(
    rownames(x = reference[[reference.assay]]),
    rownames(x = query[[query.assay]])
  )
  reference.data <-  GetAssayData(
    object = reference,
    assay = reference.assay,
    slot = "data")[features, ]
  query.data <- GetAssayData(
    object = query,
    assay = query.assay,
    slot = "data")[features, ]
  if (is.null(x = feature.mean)) {
    feature.mean <- rowMeans(x = reference.data)
    if (scale) {
      feature.sd <- sqrt(
        x = SparseRowVar2(
          mat = as.sparse(x = reference.data),
          mu = feature.mean,
          display_progress = FALSE
        )
      )
      feature.sd[is.na(x = feature.sd)] <- 1
    } else {
      feature.sd <- rep(x = 1, nrow( reference.data))
    }
    feature.mean[is.na(x = feature.mean)] <- 1
  }
  proj.data <- GetAssayData(
    object = query,
    assay = query.assay,
    slot = "data"
  )[features, ]
  store.names <- dimnames(x = proj.data)
  if (is.numeric(x = feature.mean) && feature.mean[[1]] != "SCT") {
    proj.data <- FastSparseRowScaleWithKnownStats(
      mat = as.sparse(x = proj.data),
      mu = feature.mean,
      sigma = feature.sd,
      display_progress = FALSE
    )
  }
  dimnames(x = proj.data) <- store.names
  return(proj.data)
}

ProjectCellEmbeddings <- function(
  reference,
  query,
  reduction = "pca",
  reference.assay = NULL,
  query.assay = NULL,
  dims = 1:50,
  scale = TRUE,
  verbose = TRUE,
  feature.mean = NULL,
  feature.sd = NULL
) {
  if (verbose) {
    message("Projecting cell embeddings")
  }
  reference.assay <- reference.assay %||% DefaultAssay(object = reference)
  query.assay <- query.assay %||% DefaultAssay(object = query)
  features <- rownames(x = Loadings(object = reference[[reduction]]))
  features <- intersect(x = features, y = rownames(x = query[[query.assay]]))
  proj.data <- RescaleQuery(
    reference = reference,
    query = query,
    features = features,
    scale = scale,
    feature.mean = feature.mean,
    feature.sd = feature.sd
  )
  ref.feature.loadings <- Loadings(object = reference[[reduction]])[features, dims]
  proj.pca <- t(crossprod(x = ref.feature.loadings, y = proj.data))
  return(proj.pca)
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
    counts = new(Class = "dgCMatrix"),
    data = new.expression,
    scale.data = matrix(),
    var.features = vector(),
    meta.features = data.frame(row.names = rownames(x = new.expression)),
    misc = NULL
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
    query <- SCTransform(
      object = query,
      reference.SCT.model = slot(object = reference[[reference.assay]], name = "SCTModel.list")[[1]],
      residual.features = features,
      assay = query.umi.assay,
      new.assay.name = new.sct.assay,
      verbose = FALSE
    )
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
  # features must be in both reference and query
  feature.slot <- ifelse(test = normalization.method == "SCT", yes = "scale.data", no = "data")
  query.assay.check <- query.assay
  reference.assay.check <- reference.assay
  ref.features <- rownames(x = GetAssayData(object = reference[[reference.assay.check]], slot = feature.slot))
  query.features <- rownames(x = GetAssayData(object = query[[query.assay.check]], slot = feature.slot))
  if (normalization.method == "SCT") {
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
  prediction.assay,
  label.transfer
) {
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
  valid.weight.reduction <- c("pcaproject", "pca", "cca", "rpca.ref","lsiproject", "lsi")
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
    if (is.null(x = dims)) {
      ModifyParam(param = "dims", value = 1:length(x = slot(object = anchorset, name = "command")$dims))
    }
  }

  if (!is.null(x = query)) {
    if (!isTRUE(x = all.equal(
      target = gsub(pattern = "_query", replacement = "", x = query.cells),
      current = Cells(x = query),
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
          current = Cells(x = object.list[[i]])))
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
  if (!isTRUE(x = all.equal(target = reference.cells, current = as.character(x = Cells(x = reference))))) {
    stop("The set of cells used as a reference in the AnchorSet does not match ",
         "the set of cells provided in the reference object.")
  }
  query.cells <- slot(object = anchorset, name = "query.cells")
  query.cells <- gsub(pattern = "_query", replacement = "", x = query.cells)
  if (!isTRUE(x = all.equal(target = query.cells, current = Cells(x = query), check.attributes = FALSE))) {
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
        current = Cells(x = query)
      ))) {
        stop("Cell names in the provided weight.reduction  don't ",
             "match with the cell names in the query object.", call. = FALSE)
      }
      ModifyParam(param = 'weight.reduction', value = list(NULL, weight.reduction))
    }
  }
}
