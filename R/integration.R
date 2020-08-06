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
#' \url{https://doi.org/10.1016/j.cell.2019.05.031};
#' \url{https://doi.org/10.1101/460147}
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
#' @param eps Error bound on the neighbor finding algorithm (from RANN)
#' @param verbose Print progress bars and output
#'
#' @return Returns an \code{\link{AnchorSet}} object that can be used as input to 
#' \code{\link{IntegrateData}}. 
#' 
#' @references Stuart T, Butler A, et al. Comprehensive Integration of 
#' Single-Cell Data. Cell. 2019;177:1888-1902 \url{https://doi.org/10.1016/
#' j.cell.2019.05.031}
#' 
#' @importFrom pbapply pblapply
#' @importFrom future.apply future_lapply
#' @importFrom future nbrOfWorkers
#'
#' @export
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
  reduction = c("cca", "rpca"),
  l2.norm = TRUE,
  dims = 1:30,
  k.anchor = 5,
  k.filter = 200,
  k.score = 30,
  max.features = 200,
  nn.method = "rann",
  eps = 0,
  verbose = TRUE
) {
  normalization.method <- match.arg(arg = normalization.method)
  reduction <- match.arg(arg = reduction)
  if (reduction == "rpca") {
    reduction <- "pca"
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
  object.list <- CheckDuplicateCellNames(object.list = object.list)

  slot <- "data"
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
  # if using pca, only need to compute the internal neighborhood structure once
  # for each dataset
  internal.neighbors <- list()
  if (nn.reduction == "pca") {
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
  all.anchors <- my.lapply(
    X = 1:nrow(x = combinations),
    FUN = function(row) {
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
          common.features <- intersect(
            x = rownames(x = Loadings(object = object.1[["pca"]])),
            y = rownames(x = Loadings(object = object.2[["pca"]]))
          )
          object.pair <- merge(x = object.1, y = object.2, merge.data = TRUE)
          projected.embeddings.1<- t(x = GetAssayData(object = object.1, slot = "scale.data")[common.features, ]) %*%
            Loadings(object = object.2[["pca"]])[common.features, ]
          object.pair[['projectedpca.1']] <- CreateDimReducObject(
            embeddings = rbind(projected.embeddings.1, Embeddings(object = object.2[["pca"]])),
            assay = DefaultAssay(object = object.1),
            key = "projectedpca1_"
          )
          projected.embeddings.2 <- t(x = GetAssayData(object = object.2, slot = "scale.data")[common.features, ]) %*%
            Loadings(object = object.1[["pca"]])[common.features, ]
          object.pair[['projectedpca.2']] <- CreateDimReducObject(
            embeddings = rbind(projected.embeddings.2, Embeddings(object = object.1[["pca"]])),
            assay = DefaultAssay(object = object.2),
            key = "projectedpca2_"
          )
          object.pair[["pca"]] <- CreateDimReducObject(
            embeddings = rbind(
              Embeddings(object = object.1[["pca"]]),
              Embeddings(object = object.2[["pca"]])),
            assay = DefaultAssay(object = object.1),
            key = "pca_"
          )
          reduction <- "projectedpca.1"
          reduction.2 <- "projectedpca.2"
          if (l2.norm){
            slot(object = object.pair[["projectedpca.1"]], name = "cell.embeddings") <- Sweep(
              x = Embeddings(object = object.pair[["projectedpca.1"]]),
              MARGIN = 2,
              STATS = apply(X = Embeddings(object = object.pair[["projectedpca.1"]]), MARGIN = 2, FUN = sd),
              FUN = "/"
            )
            slot(object = object.pair[["projectedpca.2"]], name = "cell.embeddings") <- Sweep(
              x = Embeddings(object = object.pair[["projectedpca.2"]]),
              MARGIN = 2,
              STATS = apply(X = Embeddings(object = object.pair[["projectedpca.2"]]), MARGIN = 2, FUN = sd),
              FUN = "/"
            )
            object.pair <- L2Dim(object = object.pair, reduction = "projectedpca.1")
            object.pair <- L2Dim(object = object.pair, reduction = "projectedpca.2")
            reduction <- paste0(reduction, ".l2")
            reduction.2 <- paste0(reduction.2, ".l2")
          }
          object.pair
        },
        stop("Invalid reduction parameter. Please choose either cca or rpca")
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
  command <- LogSeuratCommand(object = object.list[[1]], return.command = TRUE)
  anchor.set <- new(Class = "AnchorSet",
                    object.list = object.list,
                    reference.objects = reference %||% seq_along(object.list),
                    anchors = all.anchors,
                    offsets = offsets,
                    anchor.features = anchor.features,
                    command = command
  )
  return(anchor.set)
}

#' Find transfer anchors
#'
#' Find a set of anchors between a reference and query object. These 
#' anchors can later be used to transfer data from the reference to 
#' query object using the \code{\link{TransferData}} object.
#'
#' The main steps of this procedure are outlined below. For a more detailed 
#' description of the methodology, please see Stuart, Butler, et al Cell 2019.
#' \url{https://doi.org/10.1016/j.cell.2019.05.031};
#' \url{https://doi.org/10.1101/460147}
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
#'   query for this dimensional reduction step. If \code{l2.norm} is set to 
#'   \code{TRUE}, perform L2 normalization of the embedding vectors.}
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
#' @param query.assay Name of the Assay to use from query
#' @param reduction Dimensional reduction to perform when finding anchors. 
#' Options are:
#' \itemize{
#'    \item{pcaproject: Project the PCA from the reference onto the query. We 
#'    recommend using PCA when reference and query datasets are from scRNA-seq}
#'    \item{cca: Run a CCA on the reference and query }
#' }
#' @param project.query Project the PCA from the query dataset onto the 
#' reference. Use only in rare cases where the query dataset has a much larger 
#' cell number, but the reference dataset has a unique assay for transfer.
#' @param features Features to use for dimensional reduction. If not specified,
#' set as variable features of the reference object which are also present in 
#' the query. 
#' @param normalization.method Name of normalization method used: LogNormalize 
#' or SCT
#' @param npcs Number of PCs to compute on reference. If null, then use an 
#' existing PCA structure in the reference object
#' @param l2.norm Perform L2 normalization on the cell embeddings after 
#' dimensional reduction
#' @param dims Which dimensions to use from the reduction to specify the 
#' neighbor search space
#' @param k.anchor How many neighbors (k) to use when finding anchors
#' @param k.filter How many neighbors (k) to use when filtering anchors
#' @param k.score How many neighbors (k) to use when scoring anchors
#' @param max.features The maximum number of features to use when specifying the 
#' neighborhood search space in the anchor filtering
#'@param nn.method Method for nearest neighbor finding. Options include: rann,
#' annoy
#' @param eps Error bound on the neighbor finding algorithm (from 
#' \code{\link{RANN}})
#' @param approx.pca Use truncated singular value decomposition to approximate 
#' PCA
#' @param verbose Print progress bars and output
#'
#' @return Returns an \code{AnchorSet} object that can be used as input to 
#' \code{\link{TransferData}}
#' 
#' @references Stuart T, Butler A, et al. Comprehensive Integration of 
#' Single-Cell Data. Cell. 2019;177:1888-1902 \url{https://doi.org/10.1016/
#' j.cell.2019.05.031};
#' 
#' @export
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
  normalization.method = c("LogNormalize", "SCT"),
  reference.assay = NULL,
  query.assay = NULL,
  reference.reduction = "pca", 
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
  nn.method = "rann",
  eps = 0,
  scale = TRUE, 
  approx.pca = TRUE,
  verbose = TRUE
) {
  if (length(x = reference) > 1 | length(x = query) > 1) {
    stop("We currently only support transfer between a single query and reference")
  }
  if (!reduction %in% c("pcaproject", "cca", "pcaqueryproject", "lda")) {
    stop("Please select either pcaproject, cca, or pcaqueryproject for the reduction parameter.")
  }
  if (reduction %in% c('pcaproject', 'pcaqueryproject')) {
    projected = TRUE
  } else {
    projected = FALSE
  }
  normalization.method <- match.arg(arg = normalization.method)
  query <- RenameCells(
    object = query,
    new.names = paste0(Cells(x = query), "_", "query")
  )
  reference <- RenameCells(
    object = reference,
    new.names = paste0(Cells(x = reference), "_", "reference")
  )
  features <- features %||% VariableFeatures(object = reference)
  reference.assay <- reference.assay %||% DefaultAssay(object = reference)
  query.assay <- query.assay %||% DefaultAssay(object = query)
  DefaultAssay(object = reference) <- reference.assay
  DefaultAssay(object = query) <- query.assay
  feature.mean <- NULL
  slot <- "data"
  reference[[reference.reduction  ]]@assay.used <- reference.assay

  if (normalization.method == "SCT") {
    features <- intersect(x = features, y = rownames(x = query))
    query <- GetResidual(object = query, features = features, verbose = FALSE)
    query[[query.assay]] <- CreateAssayObject(
      counts =  as.sparse(x = GetAssayData(object = query[[query.assay]], slot = "scale.data")[features, ])
    )
    query <- SetAssayData(
      object = query,
      slot = "data",
      assay = query.assay,
      new.data = GetAssayData(object = query[[query.assay]], slot = "counts")
    )
    query <- SetAssayData(
      object = query,
      slot = "scale.data",
      assay = query.assay,
      new.data = as.matrix(x = GetAssayData(object = query[[query.assay]], slot = "counts"))
    )
    if (IsSCT(assay = reference[[reference.assay]])) {
      reference <- GetResidual(object = reference, features = features, verbose = FALSE)
    }
    reference[[reference.assay]] <- CreateAssayObject(
      counts =  as.sparse(x = GetAssayData(object = reference[[reference.assay]], slot = "scale.data")[features, ])
    )
    reference <- SetAssayData(
      object = reference,
      slot = "data",
      assay = reference.assay,
      new.data = GetAssayData(object = reference[[reference.assay]], slot = "counts")
    )
    reference <- SetAssayData(
      object = reference,
      slot = "scale.data",
      assay = reference.assay,
      new.data =  as.matrix(x = GetAssayData(object = reference[[reference.assay]], slot = "counts"))
    )
    feature.mean <- "SCT"
    slot <- "scale.data"
  }
  reference <- DietSeurat(object = reference,
                          assays = reference.assay,
                          dimreducs = reference.reduction, 
                          features = features, 
                          scale.data = TRUE )
  query <- DietSeurat(object = query,
                      assays = query.assay, 
                      features = features,
                      scale.data = TRUE  )
  
  ## find anchors using PCA projection
  if (reduction == 'pcaproject') {
    if (project.query) {
      if (!is.null(x = npcs)) {
        if (verbose) {
          message("Performing PCA on the provided query using ", length(x = features), " features as input.")
        }
        if (normalization.method == "LogNormalize") {
          query <- ScaleData(object = query, features = features, verbose = FALSE)
        }
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
        if (normalization.method == "LogNormalize") {
          reference <- ScaleData(object = reference, features = features, verbose = FALSE)
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
        dims = dims,
        scale = scale,
        feature.mean = feature.mean,
        verbose = verbose
      )
      ref.pca <- Embeddings(object = reference[[reference.reduction]])[, dims]
      combined.pca <- CreateDimReducObject(
        embeddings = as.matrix(x = rbind(ref.pca, projected.pca))[, dims],
        key = "ProjectPC_",
        assay = reference.assay
      )
      combined.ob <- merge(
        x = DietSeurat(object = reference, features = features), 
        y = DietSeurat(object = query, features = features)
      )
      combined.ob[["pcaproject"]] <- combined.pca
      old.loadings <- Loadings(object = reference[[reference.reduction]])
      colnames(x = old.loadings) <- paste0("ProjectPC_", 1:ncol(x = old.loadings))
      Loadings(object = combined.ob[["pcaproject"]]) <- old.loadings[, dims]
    }
  }
  ## find anchors using CCA
  if (reduction == 'cca') {
    if (normalization.method == "LogNormalize") {
      reference <- ScaleData(object = reference, features = features, do.scale = scale,verbose = FALSE)
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
  }
  
  if (reduction == "lda") {
    query.lda <- PredictLDA(reference = reference, 
                            query = query,
                            model = reference[[reference.reduction]]@misc$model,
                            reference.assay = reference.assay,
                            query.assay = query.assay, 
                            features = features, 
                            normalization.method = normalization.method )
    combined.ob <- merge(x = reference, y = query)
    combine.lda.embeddings <- rbind(Embeddings(reference[[reference.reduction]])[, dims], 
                                    query.lda$x[, dims])
    combined.ob[["ldaproject"]] <- CreateDimReducObject(embeddings = combine.lda.embeddings,
                                                        key =  "projLDA_", 
                                                        assay = reference.assay)
    reduction <- "ldaproject"
    k.filter <- NA
  }
  if (l2.norm) {
    combined.ob <- L2Dim(object = combined.ob, reduction = reduction)
    reduction <- paste0(reduction, ".l2")
  }
  slot <- "data"
  anchors <- FindAnchors(
    object.pair = combined.ob,
    assay = c(reference.assay, query.assay),
    slot = slot,
    cells1 = colnames(x = reference),
    cells2 = colnames(x = query),
    reduction = reduction,
    internal.neighbors = list(NULL, NULL),
    dims = dims,
    k.anchor = k.anchor,
    k.filter = k.filter,
    k.score = k.score,
    max.features = max.features,
    nn.method = nn.method,
    eps = eps,
    projected = projected,
    verbose = verbose
  )
  command <- LogSeuratCommand(object = combined.ob, return.command = TRUE)
  anchor.set <- new(
    Class = "AnchorSet",
    object.list = list(combined.ob),
    reference.cells = colnames(x = reference),
    query.cells = colnames(x = query),
    anchors = anchors,
    anchor.features = features,
    command = command
  )
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
#' \url{https://doi.org/10.1016/j.cell.2019.05.031};
#' \url{https://doi.org/10.1101/460147}
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
#' @param sample.tree Specify the order of integration. If NULL, will compute 
#' automatically.
#' @param preserve.order Do not reorder objects based on size for each pairwise 
#' integration.
#' @param do.cpp Run cpp code where applicable
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
#' Single-Cell Data. Cell. 2019;177:1888-1902 \url{https://doi.org/10.1016/
#' j.cell.2019.05.031}
#'
#' @export
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
  do.cpp = TRUE,
  eps = 0,
  verbose = TRUE
) {
  normalization.method <- match.arg(arg = normalization.method)
  reference.datasets <- slot(object = anchorset, name = 'reference.objects')
  object.list <- slot(object = anchorset, name = 'object.list')
  anchors <- slot(object = anchorset, name = 'anchors')
  ref <- object.list[reference.datasets]
  features <- features %||% slot(object = anchorset, name = "anchor.features")
  unintegrated <- merge(
    x = object.list[[1]],
    y = object.list[2:length(x = object.list)]
  )
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
    for (i in 1:length(x = object.list)) {
      assay <- DefaultAssay(object = object.list[[i]])
      if (length(x = setdiff(x = features.to.integrate, y = features)) != 0) {
        object.list[[i]] <- GetResidual(
          object = object.list[[i]], 
          features = setdiff(x = features.to.integrate, y = features), 
          verbose = verbose
        )
      }
      object.list[[i]][[assay]] <- CreateAssayObject(
        data = GetAssayData(object = object.list[[i]], assay = assay, slot = "scale.data")
      )
    }
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
    do.cpp = do.cpp,
    eps = eps,
    verbose = verbose
  )
  if (length(x = reference.datasets) == length(x = object.list)) {
    if (normalization.method == "SCT") {
      reference.integrated <- SetAssayData(
        object = reference.integrated,
        assay = new.assay.name,
        slot = "scale.data",
        new.data = ScaleData(
          object = GetAssayData(object = reference.integrated, assay = new.assay.name, slot = "scale.data"),
          do.scale = FALSE,
          do.center = TRUE,
          verbose = FALSE
        )
      )
      reference.integrated[[assay]] <- unintegrated[[assay]]
    }
    return(reference.integrated)
  } else {
    active.assay <- DefaultAssay(object = ref[[1]])
    reference.integrated[[active.assay]] <- NULL
    reference.integrated[[active.assay]] <- CreateAssayObject(
      data = GetAssayData(
        object = reference.integrated[[new.assay.name]],
        slot = 'data'
      )
    )
    DefaultAssay(object = reference.integrated) <- active.assay
    reference.integrated[[new.assay.name]] <- NULL
    VariableFeatures(object = reference.integrated) <- features
    # Extract the query objects (if any) and map to reference
    integrated.data <- MapQuery(
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
      sample.tree = sample.tree,
      preserve.order = preserve.order,
      do.cpp = do.cpp,
      eps = eps,
      verbose = verbose
    )

    # Construct final assay object
    integrated.assay <- CreateAssayObject(
      data = integrated.data
    )
    if (normalization.method == "SCT") {
      integrated.assay <- SetAssayData(
        object = integrated.assay,
        slot = "scale.data",
        new.data =  ScaleData(
          object = GetAssayData(object = integrated.assay, slot = "data"),
          do.scale = FALSE,
          do.center = TRUE,
          verbose = FALSE
        )
      )
      integrated.assay <- SetAssayData(
        object = integrated.assay,
        slot = "data",
        new.data = GetAssayData(object = integrated.assay, slot = "scale.data")
      )
    }
    unintegrated[[new.assay.name]] <- integrated.assay
    unintegrated <- SetIntegrationData(
      object = unintegrated,
      integration.name = "Integration",
      slot = "anchors",
      new.data = anchors
    )
    DefaultAssay(object = unintegrated) <- new.assay.name
    VariableFeatures(object = unintegrated) <- features
    unintegrated[["FindIntegrationAnchors"]] <- slot(object = anchorset, name = "command")
    unintegrated <- LogSeuratCommand(object = unintegrated)
    return(unintegrated)
  }
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
#' pancreas.integrated <- IntegrateData(anchors)
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

  object.list <- lapply(
    X = 1:length(x = object.list),
    FUN = function(i) {
      vst_out <- Misc(object = object.list[[i]][[assay[i]]], slot = "vst.out")
      vst_out$cell_attr <- vst_out$cell_attr[Cells(x = object.list[[i]]), ]
      vst_out$cells_step1 <- intersect(x = vst_out$cells_step1, y = Cells(x = object.list[[i]]))
      suppressWarnings(expr = Misc(object = object.list[[i]][[assay[i]]], slot = "vst.out") <- vst_out)
      return(object.list[[i]])
    }
  )

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
      if (!IsSCT(assay = object.list[[i]][[assay[i]]])) {
        return(object.list[[i]])
      }
      obj <- if (is.null(x = sct.clip.range)) {
        GetResidual(
          object = object.list[[i]],
          features = anchor.features,
          assay = assay[i],
          verbose = FALSE
        )
      } else {
        GetResidual(
          object = object.list[[i]],
          assay = assay[i],
          features = anchor.features,
          replace.value = TRUE,
          clip.range = sct.clip.range,
          verbose = FALSE
        )
      }
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
#' @export
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
  if (length(x = features) > 0) {
    feature.ranks <- sapply(X = features, FUN = function(x) {
      ranks <- sapply(X = object.list, FUN = function(y) {
        vf <- VariableFeatures(object = y)
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
    ranks <- sapply(X = object.list, FUN = function(y) {
      vf <- VariableFeatures(object = y)
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
#' \url{https://doi.org/10.1016/j.cell.2019.05.031};
#' \url{https://doi.org/10.1101/460147}
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
#' @param refdata Data to transfer. Should be either a vector where the names
#' correspond to reference cells, or a matrix, where the column names correspond
#' to the reference cells.
#' @param weight.reduction Dimensional reduction to use for the weighting 
#' anchors. Options are:
#' \itemize{
#'    \item{pcaproject: Use the projected PCA used for anchor building}
#'    \item{pca: Use an internal PCA on the query only}
#'    \item{cca: Use the CCA used for anchor building}
#'    \item{custom DimReduc: User provided \code{\link{DimReduc}} object 
#'    computed on the query cells}
#' }
#' @param l2.norm Perform L2 normalization on the cell embeddings after
#' dimensional reduction
#' @param dims Number of dimensions to use in the anchor weighting procedure
#' @param k.weight Number of neighbors to consider when weighting anchors
#' @param sd.weight Controls the bandwidth of the Gaussian kernel for weighting
#' @param eps Error bound on the neighbor finding algorithm (from 
#' \code{\link{RANN}})
#' @param do.cpp Run cpp code where applicable
#' @param verbose Print progress bars and output
#' @param slot Slot to store the imputed data. Must be either "data" (default) 
#' or "counts"
#' @param prediction.assay Return an \code{Assay} object with the prediction
#' scores for each class stored in the \code{data} slot. 
#'
#' @return If \code{refdata} is a vector, returns a data.frame with label 
#' predictions. If \code{refdata} is a matrix, returns an Assay object where the 
#' imputed data has been stored in the provided slot.
#'
#' @references Stuart T, Butler A, et al. Comprehensive Integration of 
#' Single-Cell Data. Cell. 2019;177:1888-1902 \url{https://doi.org/10.1016/
#' j.cell.2019.05.031}
#' 
#' @export
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
  weight.reduction = 'pcaproject',
  multi.modal = FALSE, 
  l2.norm = FALSE,
  dims = 1:30,
  k.weight = 50,
  sd.weight = 1,
  eps = 0,
  do.cpp = TRUE,
  verbose = TRUE,
  slot = "data",
  prediction.assay = FALSE
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
    if(   grepl( pattern = "_reference",x = reference.cells[1])){
      colnames(x = refdata) <- paste0(colnames(x = refdata), "_reference")
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
  if (!inherits(x = weight.reduction, what = "DimReduc") && weight.reduction == 'pca') {
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
    weight.reduction <- paste0(weight.reduction, ".l2")
  }
  if (inherits(x = weight.reduction, what = "DimReduc")) {
    weight.reduction <- RenameCells(
      object = weight.reduction,
      new.names = paste0(Cells(x = weight.reduction), "_query")
    )
  } else {
    if(multi.modal){
      weight.reduction <- lapply( X = combined.ob@misc$proj.reduction,
                                  FUN = function(r) combined.ob[[r]])
      names(weight.reduction) <- unlist(combined.ob@misc$proj.reduction)
    } else{
      weight.reduction <- combined.ob[[weight.reduction]]
    }
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
    cpp = do.cpp,
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
    if (prediction.assay) {
      pa <- CreateAssayObject(data = t(x = as.matrix(x = prediction.scores)))
      return(pa)
    } else {
      return(predictions)
    }
  } else {  # case for transferring features
    reference.cell.indices <- reference.cells[anchors$cell1]
    refdata.anchors <- refdata[, reference.cell.indices]
    nfeatures <- nrow(x = refdata)
    if (verbose) {
      message(paste0("Transfering ", nfeatures, " features onto reference data"))
    }
    new.data <- refdata.anchors %*% weights
    rownames(x = new.data) <- rownames(x = refdata)
    colnames(x = new.data) <- query.cells
    if (inherits(x = new.data, what = "Matrix")) {
      new.data <- as(object = new.data, Class = "dgCMatrix")
    }
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
  nn.method = "rann",
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



FilterMultiModalAnchors <- function(
  object, 
  reference.assay.list, 
  query.assay.list, 
  slot = "data",
  features.list,
  k.nn = 20, 
  s.nn = 20, 
  k.filter = 200,
  integration.name = 'integrated',
  verbose = TRUE
){
  features.list <- lapply(X = 1:length(features.list), function(f){
    feature = unique(features.list[[f]])
    feature = intersect(feature,  rownames(object[[query.assay.list[[f]] ]] ))
    return(feature)
  } )
  neighbors <- GetIntegrationData(object = object, integration.name = integration.name, slot = 'neighbors')
  nn.cells1 <- neighbors$cells1
  nn.cells2 <- neighbors$cells2
  
  cn.data1 <- lapply(1:length(reference.assay.list), 
                     FUN = function(r){   
                       L2Norm(
                         mat = as.matrix(x = t(x = GetAssayData(
                           object = object[[reference.assay.list[[r]] ]],
                           slot = slot)[ features.list[[r]] , nn.cells1])),
                         MARGIN = 1)
                       }) 
  
  cn.data2 <- lapply(1:length(reference.assay.list), 
                     FUN = function(r){   
                       L2Norm(
                         mat = as.matrix(x = t(x = GetAssayData(
                           object = object[[query.assay.list[[r]] ]],
                           slot = slot)[ features.list[[r]] , nn.cells2])),
                         MARGIN = 1)
                     }) 
  reference.top <-  subset(object, cells = nn.cells1)
  nfeature.list <- list()
  DimR.list<- list()
  for (r in 1:length( reference.assay.list )){
    nfeature.list[[r]] <-  1:ncol( cn.data1[[r]])
    colnames(cn.data1[[r]]) <- paste0(reference.assay.list[[r]], "_", nfeature.list[[r]] )
    DimR.list[[r]] <-  paste0(reference.assay.list[[r]],".DimR" )
    reference.top[[ DimR.list[[r]]  ]] <- CreateDimReducObject( embeddings = cn.data1[[r]], 
                                                             key = paste0(reference.assay.list[[r]], "_"),
                                                             assay =  reference.assay.list[[r]] )  
    }
  
  refernce.top.weight <- FindModalityWeights(object = reference.top,
                                           reduction.list = DimR.list,
                                           dims.list = nfeature.list, 
                                           l2.norm = F,
                                           snn.far.nn = T, 
                                           k.nn = k.nn , 
                                           s.nn =  s.nn, 
                                           verbose = F
   )
  
  nn <- MultiModalNN(object = cn.data2,
               query = cn.data1,
               k.nn = k.filter, 
               knn.range = 2*k.filter,
               modality.weight = refernce.top.weight$first.modality.weight, 
               l2.norm =  FALSE, 
               nearest.dist = refernce.top.weight$params$nearest.dist,
               sigma.list = refernce.top.weight$params$sigma.list)

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
  nn.method = "rann",
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
    top.features <- sapply(X = unique(assay),
                            FUN = function(a){
                              assay.feature <- rownames(x = GetAssayData(object = object.pair,
                                                                         slot = slot, 
                                                                         assay = a))
                             features <- intersect( assay.feature, top.features)
                             return( features )
                            })
    if(length(top.features) == 2){
      top.features <- intersect(top.features[[1]], top.features[[2]])
    } else{ 
      top.features <- as.vector(top.features)
      }
    object.pair <- FilterAnchors(
      object = object.pair,
      assay = assay,
      slot = slot,
      integration.name = 'integrated',
      features = top.features,
      k.filter = k.filter,
      nn.method = nn.method,
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
  max.nn <- c(ncol(x = neighbors$nnab$nn.idx), ncol(x = neighbors$nnba$nn.idx))
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
  nn.method = "rann",
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
  if (!is.null(x = internal.neighbors[[1]])) {
    nnaa <- internal.neighbors[[1]]
    nnbb <- internal.neighbors[[2]]
  } else {
    dim.data.self <- Embeddings(object = object[[nn.reduction]])[ ,nn.dims]
    dims.cells1.self <- dim.data.self[cells1, ]
    dims.cells2.self <- dim.data.self[cells2, ]
    nnaa <- NNHelper(
      data = dims.cells1.self,
      k = k + 1,
      method = nn.method,
      eps = eps
    )
    nnbb <- NNHelper(
      data = dims.cells2.self,
      k = k + 1,
      method = nn.method,
      eps = eps
    )
  }
  if (length(x = reduction.2) > 0) {
    nnab <- NNHelper(
      data = Embeddings(object = object[[reduction.2]])[cells2, ],
      query = Embeddings(object = object[[reduction.2]])[cells1, ],
      k = k,
      method = nn.method,
      eps = eps
    )
    nnba <- NNHelper(
      data = Embeddings(object = object[[reduction]])[cells1, ],
      query = Embeddings(object = object[[reduction]])[cells2, ],
      k = k,
      method = nn.method,
      eps = eps
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
      eps = eps
    )
    nnba <- NNHelper(
      data = dims.cells1.opposite,
      query = dims.cells2.opposite,
      k = k,
      method = nn.method,
      eps = eps
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
FindWeights <- function(
  object,
  reduction = NULL,
  assay = NULL,
  integration.name = 'integrated',
  dims = 1:10,
  features = NULL,
  k = 300,
  sd.weight = 1,
  nn.method = "rann",
  eps = 0,
  verbose = TRUE,
  cpp = FALSE
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
  k <- min(k, nrow(anchors) - 1)
  anchors.cells2 <- nn.cells2[anchors[, "cell2"]]
  if (is.null(x = features)) {
    if(is.list(reduction)){
      message("Multi-modal anchors transfer")
      query.data.use <- lapply(X = reduction, 
                               FUN = function(reduc){ 
                                 Embeddings(reduc)[nn.cells2, ]
                               })  
      
      data.use <- lapply(X = reduction, 
                         FUN = function(reduc){ 
                           Embeddings(reduc)[anchors.cells2, ]
                         })  
      
    } else{
      data.use <- Embeddings(reduction)[nn.cells2, dims]
    }
  } else {
    data.use <- t(x = GetAssayData(object = object, slot = 'data', assay = assay)[features, nn.cells2])
  }
  
  if( is.list( reduction )){
    knn_2_2 <- MultiModalNN(object = data.use, 
                            query = query.data.use,
                            k.nn = k, 
                            knn.range = max(200, 5*k),
                            l2.norm =  FALSE, 
                            modality.weight = object@misc$query.modality.weight, 
                            sigma.list  = object@misc$query.sigma.list, 
                            nearest.dist = object@misc$nearest.dist)
    
    distances <- knn_2_2$nn.dists + 1e-10
    distances <-  1 - (distances / distances[, ncol(x = distances)])+ 1e-10
    cell.index <- knn_2_2$nn.idx
  } else{
    knn_2_2 <- NNHelper(
      data = data.use[anchors.cells2, ],
      query = data.use,
      k = k + 1,
      method = nn.method,
      eps = eps
    )
    distances <- knn_2_2$nn.dists[, -1]
    distances <- 1 - (distances / distances[, ncol(x = distances)])
    cell.index <- knn_2_2$nn.idx[, -1]
  }
  integration.matrix <- GetIntegrationData(
    object = object,
    integration.name = integration.name,
    slot = "integration.matrix"
  )
  if (cpp) {
    weights <- FindWeightsC(
      integration_matrix = as(integration.matrix, "dgCMatrix"),
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
  } else {
    if (verbose) {
      pb <- txtProgressBar(min = 1, max = length(x = nn.cells2), initial = 1, style = 3, file = stderr())
    }
    dist.weights <- matrix(
      data = 0,
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
      if (verbose) setTxtProgressBar(pb, cell)
    }
    if (verbose) message("")
    dist.anchor.weight <- dist.weights * anchors[, "score"]
    weights <- 1 - exp(-1 * dist.anchor.weight / (2 * (1 / sd.weight)) ^ 2)
    weights <- Sweep(x = weights, MARGIN = 2, STATS = Matrix::colSums(weights), FUN = "/")
  }
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
# @param sample.tree Specify the order of integration. If NULL, will compute automatically.
# @param preserve.order Do not reorder objects based on size for each pairwise integration.
# @param do.cpp Run cpp code where applicable
# @param eps Error bound on the neighbor finding algorithm (from \code{\link{RANN}})
# @param verbose Print progress bars and output
#
# @return Returns an integrated matrix
#
MapQuery <- function(
  anchorset,
  reference,
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
  do.cpp = TRUE,
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
        message("Integrating dataset ", dataset1, " with reference dataset")
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
        features = features,
        dims = dims,
        do.cpp = do.cpp,
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
  nn.matrix <- as(object = nn.matrix, Class = 'dgCMatrix')
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
# @param do.cpp Run cpp code where applicable
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
  do.cpp = TRUE,
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
      data = GetAssayData(ref.obj, slot = 'data')[features.to.integrate, ]
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
        weight.reduction <- rep(x = weight.reduction, times = length(x = object.list))
      } else {
        stop("Invalid input for weight.reduction. Please specify either the names of the dimension",
             "reduction for each object in the list or provide DimReduc objects.")
      }
    }
    if (length(x = weight.reduction) != length(x = object.list)) {
      stop("Please specify a dimension reduction for each object, or one dimension reduction to be used for all objects")
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
  unintegrated <- merge(
    x = object.list[[reference.objects[[1]]]],
    y = object.list[reference.objects[2:length(x = reference.objects)]]
  )
  names(x = object.list) <- as.character(-(1:length(x = object.list)))
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
      do.cpp = do.cpp,
      k.weight = k.weight,
      sd.weight = sd.weight,
      eps = eps,
      verbose = verbose
    )
    integrated.matrix <- cbind(integrated.matrix, GetAssayData(object = object.1, slot = 'data')[features.to.integrate, ])
    merged.obj[[new.assay.name]] <- CreateAssayObject(data = integrated.matrix)
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
  unintegrated <- LogSeuratCommand(object = unintegrated)
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
    message("Projecting PCA")
  }
  reference.assay <- reference.assay %||% DefaultAssay(object = reference)
  query.assay <- query.assay %||% DefaultAssay(object = query)
  features <- rownames(x = Loadings(object = reference[[reduction]]))
  features <- intersect(x = features, y = rownames(x = query[[query.assay]]))

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
    if(scale){
    feature.sd <- sqrt(
      x = SparseRowVar2(
        mat = as(object = reference.data, Class = "dgCMatrix"), 
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
  if (is.numeric(x = feature.mean) && feature.mean != "SCT") {
    proj.data <- FastSparseRowScaleWithKnownStats(
      mat = as(object = proj.data, Class = "dgCMatrix"),
      mu = feature.mean,
      sigma = feature.sd,
      display_progress = FALSE
    )
  }
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
# @param do.cpp Run cpp code where applicable
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
  features,
  dims,
  do.cpp,
  k.weight,
  sd.weight,
  eps,
  verbose
) {
  cells1 <- colnames(x = reference)
  cells2 <- colnames(x = query)
  merged.obj <- merge(x = reference, y = query, merge.data = TRUE)
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
  if (is.null(x = weight.reduction)) {
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
  } else {
    dr <- weight.reduction[[2]]
    if (!all(cells2 %in% rownames(x = dr))) {
      stop("Query cells not present in supplied DimReduc object. Set weight.reduction to a DimReduc object containing the query cells.")
    }
    if (inherits(x = dr, what = "DimReduc")) {
      dr.weights <- dr
    } else {
      dr.weights <- query[[dr]]
    }
  }
  merged.obj <- FindWeights(
    object = merged.obj,
    integration.name = integration.name,
    reduction = dr.weights,
    cpp = do.cpp,
    dims = dims,
    k = k.weight,
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
  return(integrated.matrix[, cells2])
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
  nbrsetA <- function(x) c(neighbors$nnaa$nn.idx[x, 1:k.score], neighbors$nnab$nn.idx[x, 1:k.score] + offset)
  nbrsetB <- function(x) c(neighbors$nnba$nn.idx[x, 1:k.score], neighbors$nnbb$nn.idx[x, 1:k.score] + offset)
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
  return(object)
}
