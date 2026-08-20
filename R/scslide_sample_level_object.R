#' Perform preprocessing procedures for a Seurat object to prepare for the sample-level analysis
#'
#' This is a wrapper function to apply landmark sketching, training data sketching (optional), PLS learning,
#' and a nearest neighbor workflow to a single-cell Seurat object to get the necessary components for
#' sample-level analyses later. The output object can be directly used by \code{\link{GenerateSampleObject}} to obtain
#' a sample-level matrix.
#'
#' Two modalities are supported (see the \code{modality} argument). The default multi-modal path combines an
#' existing reduction (e.g. PCA, via \code{name.reduction.1}) with the PLS embedding using the weighted nearest
#' neighbor (WNN) workflow (\code{\link{FindmmNN}}). The single-modal path builds the neighbor graph from the PLS
#' embedding alone (\code{\link{FindsmNN}}), which is useful when no complementary reduction is available or desired.
#'
#' @param object a Seurat object.
#' @param assay the assay to perform the analyses on
#' @param add.hvg a boolen to specify whether the function runs additional correlation test between features and response (Y) to
#'  obtain more response-relevant features to be included in VariableFeatures and used for PLS learning.
#'  See \code{\link{QuickCorTest}} for details. Default is TRUE.
#' @param group.by.CorTest A metadata column name to group cells by before QuickCorTest. If NULL,
#'   falls back to QuickCorTest without grouping. Default is NULL.
#' @param Y A metadata column name of responses. Will also be used as the responses for PLS learning.
#' @param sketch.training a boolen to specify whether the function performs sketching for the Seurat object to get a training subset for PLS learning
#'  (the PLS learned from subset will be projected to the full data). For efficiency, it is recommended to do so when the Seurat object is very
#'  large. If the 'assay' being specified for this Seurat object is an on-disk assay, this parameter must be set to TRUE. Default is TRUE.
#' @param group.by.Sketch A metadata column name to group cells by before sketching. If NULL,
#'   falls back to standard sketching without grouping. Default is NULL.
#' @param ncells.per.group A positive integer, named vector, or list specifying the number of cells to sample.
#'   Default is 1000 per category in 'group.by.Sketch'. See \code{\link{SketchDataByGroup}} for detailed usage.
#' @param training.assay.name The assay name of the training data.
#' @param training.sketch.method Sketching method to use for the training subset. Can be 'LeverageScore' or 'Uniform'.
#'   Default is 'Uniform' because for training we want to obtain an unbiased subset of the original data.
#'   See \code{\link{LeverageScore}} for details.
#' @param ncells.landmark A positive integer specifying the number of 'landmark cells'. The function will generate a new
#'   assay for these 'landmark cells' (please do not confuse it with the training subset). The landmark cells should be a
#'   subset that covers the diverse cell types and cell states of the data. Default is 2000.
#' @param landmark.assay.name The assay name of the landmark cells. Default is "LANDMARK".
#' @param landmark.sketch.method Sketching method to use for the landmark cells. Can be 'LeverageScore' or 'Uniform'.
#'   Default is 'LeverageScore' because we want the landmark cells to be as diverse/representative as possible.
#'   See \code{\link{LeverageScore}} for details.
#' @param ncomp Number of components to compute
#' @param pls.function PLS function from pls package to run (options: plsr, spls, cppls)
#' @param pls.reduction.name PLS dimensional reduction name
#' @param modality Whether to build the neighbor graph from two modalities or one. \code{"multi"} (default)
#'   combines \code{name.reduction.1} with the PLS embedding via WNN (\code{\link{FindmmNN}}); \code{"single"}
#'   uses the PLS embedding alone (\code{\link{FindsmNN}}). In \code{"single"} mode, \code{name.reduction.1},
#'   \code{dims.reduction.1}, \code{knn.range}, and \code{fix.wnn.weights} are ignored.
#'
#' @param k.nn The number of nearest neighbors to compute for each modality
#' @param knn.range Range parameter for nearest neighbor search. Passed to \code{\link{FindmmNN}}, which
#'   clamps it to the number of landmark cells if it exceeds that. Default is 200. Ignored when \code{modality = "single"}.
#' @param l2.norm Perform L2 normalization on the cell embeddings during the nearest neighbor process. Default is TRUE.
#' @param nn.name Name of the output neighbor object. Defaults to \code{"weighted.nn"} for
#'   \code{modality = "multi"}, and to \code{"single.nn"} for \code{modality = "single"} (unless set explicitly).
#' @param weighted.nn.name Deprecated. Use \code{nn.name} instead.
#' @param fix.wnn.weights Pre-specified modality weights. If provided, skips the calculation and uses these weights directly.
#'   Should be a list with the same length as reduction.list. Ignored when \code{modality = "single"}.
#' @param name.reduction.1 The name of the DimReduc to use as the 1st embedding (the 2nd is the PLS embedding) in the WNN process.
#'   Ignored when \code{modality = "single"}.
#' @param dims.reduction.1 The dimensions for reduction.1 to use during the WNN process. Ignored when \code{modality = "single"}.
#' @param rm.training.assay Whether to remove the training assay after running PrepareSampleObject(). This is used to reduce the
#'   memory usage when the seurat object is large (e.g., a large on-disk object). Default is FALSE.
#' @param verbose Print progress and diagnostic messages
#' @export
#' @concept scSLIDE
#' @return return a Seurat object that contains a weighted.nn Neighbor object between the landmark cells and all the other cells.
#'   the Seurat object will also contain a landmark assay (and a training assay if sketch.training == TRUE).
#'
#' @importFrom SeuratObject DefaultAssay DefaultAssay<- Assays VariableFeatures VariableFeatures<- CreateDimReducObject
#'
PrepareSampleObject <- function(
    object = NULL,
    assay = NULL,
    add.hvg = TRUE,
    group.by.CorTest = NULL,
    Y = NULL,
    sketch.training = TRUE,
    group.by.Sketch = NULL,
    ncells.per.group = 1000,
    training.assay.name = "TRAINING",
    training.sketch.method = "Uniform",
    ncells.landmark = 2000,
    landmark.assay.name = "LANDMARK",
    landmark.sketch.method = "LeverageScore",
    ncomp = 10,
    pls.function = c("plsr", "spls", "cppls"),
    pls.reduction.name = "pls",
    modality = c("multi", "single"),
    k.nn = 5,
    knn.range = 200,
    l2.norm = TRUE,
    name.reduction.1 = "pca",
    dims.reduction.1 = 1:30,
    nn.name = "weighted.nn",
    weighted.nn.name = NULL,
    fix.wnn.weights = c(0.5, 0.5),
    rm.training.assay = FALSE,
    verbose = TRUE
){
  # house-keeping checks
  assay <- assay[1L] %||% DefaultAssay(object = object)
  assay <- match.arg(arg = assay, choices = Assays(object = object))
  modality <- match.arg(arg = modality)
  # reduction.1 is only used for the multi-modal (WNN) path
  if (identical(modality, "multi")) {
    name.reduction.1 <- match.arg(arg = name.reduction.1, choices = names(object@reductions))
  }
  pls.function <- match.arg(arg = pls.function)

  # deprecation shim: 'weighted.nn.name' -> 'nn.name'
  nn.name.supplied <- !missing(nn.name)
  if (!is.null(weighted.nn.name)) {
    warning("'weighted.nn.name' is deprecated; please use 'nn.name' instead.")
    nn.name <- weighted.nn.name
    nn.name.supplied <- TRUE
  }

  # for the single-modality path, avoid the 'weighted' misnomer in the default NN name
  if (identical(modality, "single") && !nn.name.supplied) {
    nn.name <- "single.nn"
  }

  # whether to perform sketching to get a subset of the full data for PLS learning
  if(isTRUE(x = sketch.training)){
    object <- SketchDataByGroup(object,
                                assay = assay,
                                ncells = ncells.per.group,
                                sketched.assay = training.assay.name,
                                method = training.sketch.method,
                                group.by = group.by.Sketch,
                                verbose = verbose)
  } else {
    training.assay.name <- assay
  }
  
  # an extra step to make sure variableFeatures is set 
  var_feature <- VariableFeatures(object)
  if(length(var_feature) <= 0) {
      stop("Please run FindVariableFeatures() or manually set VariableFeatures for this object.")
  }
  
  # whether to perform additional cor.test to get more relevant genes for PLS training
  if(isTRUE(x = add.hvg)){
    deg_list <- QuickCorTest(object,
                             assay = training.assay.name,
                             layer = "data",
                             group.by = group.by.CorTest,
                             Y = Y,
                             verbose = verbose)
    all_HVG <- union(VariableFeatures(object), deg_list)
    VariableFeatures(object) <- all_HVG
    #
    if(verbose){
      message("The number of HVGs used for PLS learning is ", length(all_HVG))
    }
  } else {
    all_HVG <- VariableFeatures(object)
    if(verbose){
      message("Using the default VariableFeatures for PLS learning, n = ", length(all_HVG))
    }
  }
  if(length(all_HVG) < 1) stop("The number of highly variable genes (HVGs) for running PLS is 0, something is wrong.")
  
  # to get a assay for the landmark cells
  object <- SketchData(object,
                       assay = training.assay.name,
                       ncells = ncells.landmark,
                       sketched.assay = landmark.assay.name,
                       method = landmark.sketch.method,
                       features = var_feature,
                       verbose = verbose)

  # perform PLS learning for the training.assay
  DefaultAssay(object) <- training.assay.name
  object <- ScaleData(object, assay = training.assay.name, features = all_HVG)
  message("\nRunning PLS learning for the data.")
  object <- RunPLS(object = object,
                   assay = training.assay.name,
                   features = all_HVG,
                   ncomp = ncomp,
                   Y = Y,
                   pls.function = pls.function,
                   reduction.name = pls.reduction.name,
                   verbose = verbose)

  # to project the PLS back to the full data (if we are using a subset of training data)
  if(isTRUE(x = sketch.training)){
    proj.pls <- ProjectCellEmbeddings(
      query = object,
      query.assay = assay,
      reference = object,
      reference.assay = training.assay.name,
      reduction = pls.reduction.name,
      dims = 1:ncomp,
      scale = TRUE,
      normalization.method = "LogNormalize",
      verbose = verbose
    )
    # append the projected PLS embeddings to the full data
    object[[paste0("proj.", pls.reduction.name)]] <- CreateDimReducObject(embeddings = proj.pls, key = "pPLS_", assay = assay)
  }

  if(rm.training.assay && sketch.training){
    if(verbose) message("removing sketch.training.assay to reduce memory usage.")
    DefaultAssay(object) <- assay
    object[[training.assay.name]] <- NULL
    gc()
  }

  # to perform nearest neighbors between the landmarks and the other cells
  DefaultAssay(object) <- assay
  # resolve the PLS reduction to use (the projected embedding if we sketched a training subset)
  pls.reduction.use <- if (isTRUE(x = sketch.training)) paste0("proj.", pls.reduction.name) else pls.reduction.name

  if (identical(modality, "multi")) {
    message("Constructing WNN graph between the landmark cells and all other cells.")
    # get the reduction.list ready
    reduction.list <- list(name.reduction.1, pls.reduction.use)
    # get the dims.list ready
    if(dims.reduction.1[length(dims.reduction.1)] > dim(object[[name.reduction.1]])[2]){
      dims.reduction.1 <- dims.reduction.1[1]:dim(object[[name.reduction.1]])[2]
      warning("The dims.reduction.1 specified exceeds the number of components in ", name.reduction.1,
              ", so switching to ", dims.reduction.1[1], ":", dims.reduction.1[length(dims.reduction.1)])
    }
    dims.list <- list(dims.reduction.1, 1:ncomp)
  } else {
    message("Constructing single-modality NN graph ('", pls.reduction.use,
            "') between the landmark cells and all other cells.")
  }
  
  # build the NN object between landmarks and all other cells
  if (identical(modality, "multi")) {
    object = FindmmNN(object,
                      sketch.assay = landmark.assay.name,
                      reduction.list = reduction.list,
                      k.nn = k.nn,
                      knn.range = knn.range,
                      l2.norm = l2.norm,
                      nn.name = nn.name,
                      dims.list = dims.list,
                      fix.wnn.weights = fix.wnn.weights,
                      verbose = verbose
                     )
  } else {
    object = FindsmNN(object,
                      sketch.assay = landmark.assay.name,
                      reduction = pls.reduction.use,
                      dims = 1:ncomp,
                      k.nn = k.nn,
                      l2.norm = l2.norm,
                      nn.name = nn.name,
                      verbose = verbose
                     )
  }
  # return the object
  return(object)
}


#' Generate a sample-level count matrix based on landmark assay and its weighted.nn object
#'
#' Returns summed counts of weighted.nn for each landmark cell within each sample.
#'
#' @param object Seurat object
#' @param nn.name Name of the Neighbor object to use for the calculation
#'   (e.g. \code{"weighted.nn"} as produced by \code{\link{FindmmNN}} /
#'   \code{\link{PrepareSampleObject}}, or \code{"single.nn"} as produced by
#'   \code{\link{FindsmNN}}). If \code{NULL} (default), the function auto-detects
#'   a Neighbor object named \code{"weighted.nn"} or \code{"single.nn"}; if both
#'   are present it errors and asks you to specify this parameter.
#' @param return.seurat Whether to return the data as a Seurat object. Default is TRUE
#' @param k.nn the number of nearest neighbors to perform the summing
#' @param sketch.assay the name of the sketch.assay you used to perform the FindmmNN()
#' @param group.by A single meta-data column to group cells by (e.g. Donor ID); 'ident' by default.
#'  To group on a combination of columns, add a merged column to the object beforehand and pass
#'  that, e.g. \code{object$batch_rep <- paste(object$batch, object$replicate, sep = "_")}.
#' @param normalization.method Method for normalization. Supports LogNormalize and ChiSquared. see details at
#' \code{\link{NormalizeData}} and \code{\link{NormalizeChiSquared}}
#' @param scale.factor Scale factor for Log-Normalization, see \code{\link{NormalizeData}}
#' @param rename.group.by if rename.group.by is NULL, the rownames of the landmark matrix will used the original cell IDs.
#'  But, user can indicate which meta-data column to use to rename the rows.
#'  A suffix of "_LM" + order number will be added automatically.
#' @param add.meta.data if TRUE, the function will automatically detect sample-level meta-data (based on 'group.by') and append it
#'  to the sample-level object; if FALSE, it will not do so.
#' @param new_assay_name Name for the new assay containing landmark counts
#' @param cells.use An optional character vector of cell names to subset before building
#'   the sample-level matrix. When non-NULL, only cells in this vector (that are also present
#'   in the NN object) are retained. Default is NULL (use all cells).
#' @param verbose Print progress and diagnostic messages
#' @param ... Arguments to be passed to methods such as \code{\link{CreateSeuratObject}}
#'
#' @export
#' @concept scSLIDE
#' @return return a Seurat object that contains a count matrix with number of landmark as rows, sample as columns
#'
#' @importFrom Matrix sparseMatrix
#' @importFrom SeuratObject FetchData CreateSeuratObject AddMetaData
#'
GenerateSampleObject <- function(
    object,
    nn.name = NULL,
    k.nn = 5,
    sketch.assay = "LANDMARK",
    return.seurat = TRUE,
    new_assay_name = "LMC",
    group.by = 'ident',
    normalization.method = "ChiSquared",
    scale.factor = 10000,
    rename.group.by = NULL,
    add.meta.data = TRUE,
    cells.use = NULL,
    verbose = TRUE,
    ...
) {
  # Get internal function from Seurat
  
  normalization.method <- match.arg(arg = normalization.method,
                                    choices = c("ChiSquared", "LogNormalize"))

  # quick check on the group.by
  data <- FetchData(object = object, vars = rev(x = group.by))
  group.by <- intersect(group.by, colnames(data))
  if (length(group.by) < 1) stop("Please specify the correct meta-data column names.")
  if (length(group.by) > 1) {
    stop("'group.by' must be a single meta-data column, but got: ",
         paste(group.by, collapse = ", "), ". To group on a combination of columns, add a ",
         "merged column to the object beforehand and pass that, e.g. ",
         "object$batch_rep <- paste(object$batch, object$replicate, sep = \"_\").")
  }

  # get the cells in the nn object
  if (is.null(x = nn.name)) {
    # auto-detect among the default NN object names produced by FindmmNN / FindsmNN
    default.nn.names <- c("weighted.nn", "single.nn")
    present.nn <- intersect(default.nn.names, names(object@neighbors))
    if (length(present.nn) == 0) {
      stop("No default NN object ('weighted.nn' or 'single.nn') was found. ",
           "Please specify 'nn.name'.")
    }
    if (length(present.nn) > 1) {
      stop("Multiple default NN objects were found (",
           paste(present.nn, collapse = ", "),
           "). Please specify 'nn.name' to disambiguate.")
    }
    nn.name <- present.nn
  } else if (!nn.name %in% names(object@neighbors)) {
    stop("Please specify the correct name of the NN object.")
  }
  slct_cells <- Cells(object[[nn.name]])
  # the NN object must carry at least k.nn neighbours per cell
  nn.k.available <- ncol(x = object[[nn.name]]@nn.idx)
  if (k.nn > nn.k.available) {
    stop("k.nn (", k.nn, ") exceeds the number of neighbors stored in '", nn.name,
         "' (", nn.k.available, ").")
  }

  # now extract the cell identity based on group.by
  data <- data[slct_cells, , drop = F]
  data <- data[which(rowSums(x = is.na(x = data)) == 0), , drop = F]
  for (i in 1:ncol(x = data)) {
    data[, i] <- as.factor(x = data[, i])
  }
  num.levels <- sapply(
    X = 1:ncol(x = data),
    FUN = function(i) {
      length(x = levels(x = data[, i]))
    }
  )
  # 'group.by' is a single column, so a single level leaves nothing to group on
  if (num.levels < 2) {
    stop("The grouping variable '", group.by, "' has only one value ('",
         levels(x = data[, 1]), "') among the cells in '", nn.name,
         "', so no sample-level grouping is possible. Please choose a meta-data ",
         "column with at least two values.")
  }
  category.matrix <- CreateCategoryMatrix(labels = data, method = "aggregate")

  # now starts to generate cell-level count matrix
  rows <- as.vector(t(object[[nn.name]]@nn.idx[, 1:k.nn, drop = FALSE]))  # No copy if nn.idx is integer
  cols <- rep(1:length(object[[nn.name]]@cell.names), each = k.nn)
  raw_ct_mat <- sparseMatrix(i = rows, j = cols,
                             dims = c(ncol(object[[sketch.assay]]), length(object[[nn.name]]@cell.names)),
                             x = 1)
  rownames(raw_ct_mat) <- colnames(object[[sketch.assay]])
  colnames(raw_ct_mat) <- object[[nn.name]]@cell.names

  # optionally subset to user-specified cells
  if (!is.null(cells.use)) {
    keep <- which(colnames(raw_ct_mat) %in% cells.use)
    if (length(keep) == 0) stop("None of the cells in 'cells.use' found in the NN object.")
    raw_ct_mat <- raw_ct_mat[, keep, drop = FALSE]
    if (verbose) message("Subsetting to ", length(keep), " cells from 'cells.use'.")
  }

  # landmark cells must never be counted towards their own sample
  rm_col <- which(colnames(raw_ct_mat) %in% rownames(raw_ct_mat))
  if (length(x = rm_col) > 0){
    raw_ct_mat <- raw_ct_mat[, -rm_col, drop = FALSE]
  }

  # collapse the cell-level matrix into sample-level matrix, based on the category.matrix
  landmark_ct_mat <- matrix(nrow = nrow(raw_ct_mat),
                            ncol = ncol(category.matrix))
  i <- 1
  valid_cols <- colnames(raw_ct_mat)
  for(SAMPLE_ID in colnames(category.matrix)){
    cell_idx <- rownames(category.matrix)[which(category.matrix[, SAMPLE_ID] == 1)]
    cell_idx <- cell_idx[! cell_idx %in% rownames(raw_ct_mat)]
    cell_idx <- cell_idx[cell_idx %in% valid_cols]
    if (length(cell_idx) == 0) {
      landmark_ct_mat[, i] <- 0
      i <- i + 1
      next
    }
    landmark_ct_mat[, i] <- SeuratObject::rowSums(raw_ct_mat[, cell_idx, drop = FALSE])
    i <- i + 1
  }

  # assign row and col names to the count matrix
  landmark.ids <- colnames(object[[sketch.assay]])
  if (!is.null(rename.group.by)) {
    mdata <- FetchData(object = object, vars = rev(x = rename.group.by))
    rename.group.by <- intersect(rename.group.by, colnames(mdata))
    if(length(x = rename.group.by) > 1) {
      rename.group.by <- rename.group.by[1]
    }
    #
    if (length(x = rename.group.by) < 1) {
      # FetchData can return columns under names other than the ones requested;
      # fall back to the landmark cell IDs rather than pasting a zero-length vector
      warning("Cannot find the correct meta-data columns to rename the landmark matrix.")
      row.labels <- landmark.ids
    } else {
      row.labels <- mdata[landmark.ids, rename.group.by]
    }
  } else {
    row.labels <- landmark.ids
  }
  rownames(landmark_ct_mat) <- paste0(row.labels, "_LM", 1:nrow(landmark_ct_mat))
  rownames(landmark_ct_mat) <- gsub(pattern = "\\s", replacement = "_", rownames(landmark_ct_mat))
  colnames(landmark_ct_mat) <- colnames(category.matrix)

  if(return.seurat == FALSE) {
    return(landmark_ct_mat)
  }

  # generate a new Seurat object for the sample-level matrix
  landmark_obj <- CreateSeuratObject(counts = landmark_ct_mat, assay = new_assay_name)

  # we can now add donor-level meta data to this sample-level object and perform some basic QC and processings
  if(isTRUE(x = add.meta.data)){
    meta_data <- object@meta.data
    find_donor_specific_cols <- function(meta_data, donor_col = group.by) {
      # Check each column for within-donor consistency
      donor_specific <- sapply(names(meta_data), function(col_name) {
        if (col_name == donor_col) return(TRUE)  # Always include donor column

        # Split by donor and check if each donor has only 1 unique value
        by_donor <- split(meta_data[[col_name]], meta_data[[donor_col]])
        all(sapply(by_donor, function(x) length(unique(x)) == 1))
      })
      return(names(donor_specific)[donor_specific])
    }

    meta_data <- meta_data[, find_donor_specific_cols(meta_data, group.by)]
    meta_data <- meta_data[!duplicated(meta_data), ]
    # rownames(meta_data) <- meta_data[[group.by]]
    rownames(meta_data) <- gsub(pattern = "\\_", replacement = "-", meta_data[[group.by]])
    landmark_obj <- AddMetaData(landmark_obj, metadata = meta_data)
  }

  # normalization
  if(normalization.method == "LogNormalize"){
    landmark_obj <- NormalizeData(landmark_obj, normalization.method = normalization.method, scale.factor = scale.factor)
  } else if(normalization.method == "ChiSquared"){
    landmark_obj <- NormalizeChiSquared(object = landmark_obj, assay = new_assay_name)
  }
  #
  return(landmark_obj)
}

