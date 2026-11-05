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
#'   No cell in 'assay' can have a missing (\code{NA}) response, including cells that would not be
#'   sketched into the training subset.
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
#' @return return a Seurat object that contains a weighted.nn Neighbor object between the landmark cells and all the other cells.
#'   the Seurat object will also contain a landmark assay (and a training assay if sketch.training == TRUE).
#'
#' @importFrom SeuratObject DefaultAssay DefaultAssay<- Assays VariableFeatures VariableFeatures<- CreateDimReducObject
#'
#' @concept sample_level_analysis
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
  # missing responses would only fail once PLS is reached, and only if the
  # affected cells happen to be sketched into the training subset
  .check_Y_missing(
    Y = if (is.character(x = Y)) {
      FetchData(object = object, vars = Y, cells = Cells(x = object[[assay]]), clean = FALSE)
    } else {
      Y
    },
    label = "Y"
  )

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
#' @return return a Seurat object that contains a count matrix with number of landmark as rows, sample as columns
#'
#' @importFrom Matrix sparseMatrix
#' @importFrom SeuratObject FetchData CreateSeuratObject AddMetaData
#'
#' @concept sample_level_analysis
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
    # only the cells that were grouped above; 'group.by' may not be a meta-data
    # column, as FetchData also resolves e.g. 'ident'
    meta_data <- object@meta.data[rownames(x = data), , drop = FALSE]
    if (!group.by %in% colnames(x = meta_data)) {
      meta_data[[group.by]] <- data[[group.by]]
    }
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

    meta_data <- meta_data[, find_donor_specific_cols(meta_data, group.by), drop = FALSE]
    meta_data <- meta_data[!duplicated(meta_data), , drop = FALSE]
    # rownames(meta_data) <- meta_data[[group.by]]
    rownames(meta_data) <- gsub(pattern = "\\_", replacement = "-", meta_data[[group.by]])
    landmark_obj <- AddMetaData(landmark_obj, metadata = meta_data)
  }

  # normalization
  if(normalization.method == "LogNormalize"){
    landmark_obj <- NormalizeData(landmark_obj, normalization.method = normalization.method, scale.factor = scale.factor)
  } else if(normalization.method == "ChiSquared"){
    landmark_obj <- NormalizeChiSquared(object = landmark_obj, assay = new_assay_name, verbose = verbose)
  }
  #
  return(landmark_obj)
}

#' Sketch Data by Group
#'
#' This function uses sketching methods to downsample high-dimensional single-cell RNA expression data
#' within specified groups/categories, which can help with scalability for large datasets while
#' maintaining representation across different cell types or conditions.
#'
#' @param object A Seurat object.
#' @param group.by A metadata column name to group cells by before sketching. If NULL,
#'   falls back to standard sketching without grouping. Default is NULL.
#' @param assay Assay name. Default is NULL, in which case the default assay of the object is used.
#' @param ncells A positive integer, named vector, or named list specifying the number of cells to sample.
#'   \itemize{
#'     \item If a single integer: same number of cells sampled from each group and layer
#'     \item If a named vector with group names: specific number per group
#'     \item If a named list with group names containing layer vectors: specific number per group per layer
#'   }
#'   Without \code{group.by}, this is the number of cells per layer, and names refer to layers.
#'   Default is 1000.
#' @param sketched.assay Sketched assay name. A sketch assay is created or overwritten with the sketch data.
#'   Default is 'sketch'.
#' @param method Sketching method to use. Can be 'LeverageScore' or 'Uniform'.
#'   Default is 'Uniform'.
#' @param var.name A metadata column name to store the leverage scores. Default is 'leverage.score'.
#' @param over.write Whether to overwrite existing column in the metadata. Default is FALSE.
#' @param seed A positive integer for the seed of the random number generator. Default is 123.
#' @param cast The type to cast the resulting assay to. Default is NULL (no
#'   casting), which preserves the original matrix format (e.g. on-disk
#'   BPCells \code{IterableMatrix} objects stay on disk). Set to
#'   \code{"dgCMatrix"} to force the sketched assay into an in-memory sparse
#'   matrix.
#' @param verbose Print progress and diagnostic messages. Default is TRUE.
#' @param features A character vector of feature names to include in the sketched assay.
#' @param min.cells.per.group Minimum number of cells required per group to perform sketching.
#'   Groups with fewer cells will be included entirely. Default is 10.
#' @param ... Arguments passed to other methods
#'
#' @return A Seurat object with the sketched data added as a new assay. The metadata will contain
#'   information about which groups each sketched cell belongs to.
#'
#' @details
#' Normalized or transformed expression must be stored in a layer named
#' \code{data}, or split layers such as \code{data.batch1}
#'
#' When \code{group.by} is specified, the function performs the following steps: \cr
#' 1. Splits cells into groups based on the metadata column \cr
#' 2. Calculates leverage scores (if method = 'LeverageScore') within each group \cr
#' 3. Samples the specified number of cells from each group \cr
#' 4. Combines all sampled cells into the sketched assay
#'
#' This approach ensures that rare cell types or conditions are not underrepresented
#' in the final sketched dataset.
#'
#' @examples
#' \dontrun{
#' # Basic usage with grouping by cell type
#' sketched_obj <- SketchDataByGroup(
#'   object = seurat_obj,
#'   group.by = "cell_type",
#'   ncells = 500
#' )
#'
#' # Different number of cells per group
#' sketched_obj <- SketchDataByGroup(
#'   object = seurat_obj,
#'   group.by = "condition",
#'   ncells = c("control" = 1000, "treatment" = 800)
#' )
#'
#' # Without grouping (falls back to original behavior)
#' sketched_obj <- SketchDataByGroup(
#'   object = seurat_obj,
#'   ncells = 5000
#' )
#' }
#'
#' @importFrom SeuratObject CastAssay Key Key<- Layers DefaultAssay Assays Cells
#' @importFrom rlang abort
#' @importFrom stats setNames
#'
#' @export
#' @concept sketching
#'
SketchDataByGroup <- function(
    object,
    group.by = NULL,
    assay = NULL,
    ncells = 1000L,
    sketched.assay = 'sketch',
    method = c('Uniform', 'LeverageScore'),
    var.name = "leverage.score",
    over.write = FALSE,
    seed = 123L,
    cast = NULL,
    verbose = TRUE,
    features = NULL,
    min.cells.per.group = 10L,
    ...
) {
  # Get internal functions from Seurat
  
  # Input validation and setup
  assay <- assay[1L] %||% DefaultAssay(object = object)
  assay <- match.arg(arg = assay, choices = Assays(object = object))
  method <- match.arg(arg = method)
  layer.names <- suppressWarnings(Layers(object = object[[assay]], search = 'data'))

  if (length(x = layer.names) == 0L) {
    abort(message = paste0(
      "No 'data' layer found. Store normalized or transformed expression in 'data' ",
      "(or split layers such as 'data.batch1') before running SketchDataByGroup()"
    ))
  }

  if (sketched.assay == assay) {
    abort(message = "Cannot overwrite existing assays")
  }
  if (sketched.assay %in% Assays(object = object)) {
    if (sketched.assay == DefaultAssay(object = object)) {
      DefaultAssay(object = object) <- assay
    }
    object[[sketched.assay]] <- NULL
  }
  if (!over.write) {
    var.name <- CheckMetaVarName(object = object, var.name = var.name)
  }

  # If no grouping specified, use original function logic
  if (is.null(x = group.by)) {
    return(SketchData(
      object = object,
      assay = assay,
      ncells = ncells,
      sketched.assay = sketched.assay,
      method = method,
      var.name = var.name,
      over.write = over.write,
      seed = seed,
      cast = cast,
      verbose = verbose,
      features = features,
      ...
    ))
  }

  # Validate group.by parameter
  if (!group.by %in% colnames(object[[]])) {
    abort(message = paste0("Column '", group.by, "' not found in object metadata"))
  }

  # Get grouping information
  group_data <- object[[group.by]]
  groups <- as.character(unique(group_data[, 1]))
  groups <- groups[!is.na(groups)]  # Remove NA values

  if (verbose) {
    message(paste0("Found ", length(groups), " groups in '", group.by, "': ",
                   paste(groups, collapse = ", ")))
  }

  # Process ncells parameter for groups
  if (is.list(x = ncells)) {
    if (is.null(x = names(x = ncells)) || !all(groups %in% names(x = ncells))) {
      abort(message = "If specified as a list, ncells must contain a named entry for every group")
    }
    ncells_per_group <- ncells[groups]
  } else if (length(ncells) == 1) {
    ncells_per_group <- setNames(rep(ncells, length(groups)), groups)
  } else if (is.null(names(ncells))) {
    if (length(ncells) != length(groups)) {
      abort(message = "When ncells is a vector without names, length must match number of groups")
    }
    ncells_per_group <- setNames(ncells, groups)
  } else {
    # Named vector - check that all groups are covered
    missing_groups <- setdiff(groups, names(ncells))
    if (length(missing_groups) > 0) {
      if (verbose) {
        message(paste0("Using default ncells (", ncells[1], ") for groups: ",
                       paste(missing_groups, collapse = ", ")))
      }
      ncells_per_group <- setNames(rep(ncells[1], length(groups)), groups)
      ncells_per_group[names(ncells)] <- ncells
    } else {
      ncells_per_group <- ncells[groups]
    }
  }

  # Calculate leverage scores if needed
  if (method == 'LeverageScore') {
    if (verbose) {
      message("Calculating Leverage Scores")
    }
    object <- LeverageScore(
      object = object,
      assay = assay,
      var.name = var.name,
      over.write = over.write,
      seed = seed,
      verbose = FALSE,
      features = features,
      ...
    )
  } else if (method == 'Uniform') {
    if (verbose) {
      message("Using uniform sampling")
    }
    object[[var.name]] <- 1
  }

  leverage.score <- object[[var.name]]

  # Sample cells within each group
  all_sampled_cells <- c()

  for (group in groups) {
    if (verbose) {
      message(paste0("Processing group: ", group))
    }

    # Get cells in this group
    group_cells <- rownames(object[[]])[group_data[, 1] == group & !is.na(group_data[, 1])]

    if (length(group_cells) == 0) {
      if (verbose) {
        message(paste0("No cells found in group '", group, "', skipping"))
      }
      next
    }

    if (length(group_cells) < min.cells.per.group) {
      if (verbose) {
        message(paste0("Group '", group, "' has ", length(group_cells),
                       " cells (< min.cells.per.group), including all cells"))
      }
      all_sampled_cells <- c(all_sampled_cells, group_cells)
      next
    }

    # Sample cells from each layer within this group
    group_sampled_cells <- c()
    ncells_group <- ncells_per_group[[group]]

    for (layer_name in layer.names) {
      layer_cells <- Cells(object[[assay]], layer = layer_name)
      group_layer_cells <- intersect(group_cells, layer_cells)
      if (length(group_layer_cells) == 0) {
        next
      }

      # Calculate how many cells to sample from this layer for this group
      ncells_layer <- if (is.null(names(ncells_group))) ncells_group else ncells_group[layer_name]
      if (anyNA(ncells_layer)) {
        abort(message = paste0("Please specify ncells for group '", group, "' in layer '", layer_name, "'"))
      }

      ncells_layer <- min(ncells_layer, length(group_layer_cells))

      if (length(group_layer_cells) <= ncells_layer) {
        cells_to_keep <- group_layer_cells
      } else {
        if (!is.null(seed)) {
          set.seed(seed + which(groups == group) + which(layer.names == layer_name))
        }
        cells_to_keep <- sample(
          x = group_layer_cells,
          size = ncells_layer,
          prob = leverage.score[group_layer_cells, ]
        )
      }

      group_sampled_cells <- c(group_sampled_cells, cells_to_keep)
    }

    if (verbose) {
      message(paste0("Sampled ", length(group_sampled_cells), " cells from group '", group, "'"))
    }

    all_sampled_cells <- c(all_sampled_cells, group_sampled_cells)
  }

  if (length(all_sampled_cells) == 0) {
    abort(message = "No cells were sampled. Check your grouping variable and ncells parameters.")
  }

  if (verbose) {
    message(paste0("Total sampled cells: ", length(all_sampled_cells)))
  }

  # Create sketched assay
  sketched <- suppressWarnings(expr = subset(
    x = object[[assay]],
    cells = all_sampled_cells,
    layers = Layers(object = object[[assay]], search = c('counts', 'data'))
  ))

  if (!is.null(x = cast) && inherits(x = sketched, what = 'Assay5')) {
    sketched <- CastAssay(object = sketched, to = cast, ...)
  }

  Key(object = sketched) <- Key(object = sketched.assay, quiet = TRUE)
  object[[sketched.assay]] <- sketched
  DefaultAssay(object = object) <- sketched.assay

  if (verbose) {
    message(paste0("Created sketched assay '", sketched.assay, "' with ",
                   ncol(sketched), " cells"))
  }

  return(object)
}

#' Find Multi-Modal Nearest Neighbors
#'
#' A wrapper function that finds Multi-Modal weights and generate a Multi-Modal NN object.
#'
#' This function is essentially a wrapper function that perform the following 2 procedures:
#' \itemize{
#'   \item Seurat:::FindModalityWeights
#'   \item Seurat:::MultiModalNN
#' }
#'
#' @param object A Seurat object
#' @param sketch.assay The name of the sketched assay of the landmark cells.
#' @param reduction.list A list of dimensional reductions, one for each modality
#' @param dims.list A list containing the dimensions for each reduction to use
#' @param k.nn The number of nearest neighbors to compute for each modality.
#'   Internally clamped to a minimum of 20 (i.e. values below 20 are silently
#'   raised to 20).
#' @param knn.range Range parameter for nearest neighbor search
#' @param l2.norm Perform L2 normalization on the cell embeddings
#' @param fix.wnn.weights Pre-specified modality weights. If provided, skips the calculation
#'   and uses these weights directly. Should be a list with the same length as reduction.list.
#' @param nn.name Name under which the resulting Neighbor object is stored.
#' @param weighted.nn.name Deprecated. Use \code{nn.name} instead.
#' @param verbose Print progress bars and output
#' @param ...	Additional arguments for \code{Seurat:::FindModalityWeights} (for
#'   example \code{sigma.idx}, \code{sd.scale} or \code{prune.SNN}) or
#'   \code{Seurat:::MultiModalNN} (for example \code{kernel.power}). Each is
#'   routed to whichever of the two declares it; arguments matching neither are
#'   reported.
#' 
#' @return return a Seurat object that contains a weighted.nn Neighbor object between the landmark cells and all the other cells.
#' @export
#'
#' @concept sample_level_analysis
FindmmNN <- function(
    object,
    sketch.assay = NULL,
    reduction.list,
    dims.list,
    k.nn = 20,
    knn.range = 200,
    l2.norm = TRUE,
    fix.wnn.weights = NULL,
    nn.name = "weighted.nn",
    weighted.nn.name = NULL,
    verbose = TRUE,
    ...
){
  # deprecation shim: 'weighted.nn.name' -> 'nn.name'
  if (!is.null(weighted.nn.name)) {
    warning("'weighted.nn.name' is deprecated; please use 'nn.name' instead.")
    nn.name <- weighted.nn.name
  }

  # Get internal functions from Seurat

  if(k.nn <= 2) {
    stop("k.nn must be a integer larger than 1.")
  }
  # FindModalityWeights() and MultiModalNN() each ignore the dots meant for the
  # other, so both can be handed the dots directly. CheckDots reports an
  # argument belonging to neither, which would otherwise be absorbed silently.
  #
  # Do not route these through do.call(). LogSeuratCommand() deparses the
  # calling frame, and under do.call() that frame carries the function body and
  # every argument by value, so deparsing walks the whole Seurat object --
  # 6-7x slower here at 30k cells and worse as the object grows. See
  # satijalab/seurat#10153 for the same failure in SCTransform.
  # CheckDots() is not usable here: it only reports unused arguments when none
  # of the functions named accepts dots, and both of these now do.
  unused <- setdiff(
    x = names(x = list(...)),
    y = c(names(x = formals(fun = FindModalityWeights)),
          names(x = formals(fun = MultiModalNN)))
  )
  if (length(x = unused) > 0) {
    warning(
      "The following arguments are not used: ", paste(unused, collapse = ', '),
      call. = FALSE,
      immediate. = TRUE
    )
  }
  # calculate WNN weights
  wnn.weight <- FindModalityWeights(object = object,
                                    reduction.list = reduction.list,
                                    dims.list = dims.list,
                                    fix.wnn.weights = fix.wnn.weights,
                                    snn.far.nn = FALSE,
                                    l2.norm = l2.norm,
                                    k.nn = ifelse(20 >= k.nn, 20, k.nn),
                                    verbose = verbose,
                                    ...)

  # prepare for the objects used to search for WNN
  if (is.null(sketch.assay) | !sketch.assay %in% names(object@assays)){
    stop("sketch.assay is not correctly defined. please check.")
  }
  sketch_cells <- Cells(object[[sketch.assay]])
  
  # check if knn.range is too large and needs to be fixed
  if(knn.range > length(sketch_cells)){
      knn.range <- length(sketch_cells)
      warning("The number of knn.range is exceeding the number of available landmarks. Resetting it to length(landmarks) = ", 
              length(sketch_cells))
  }
  
  sketch_obj <- subset(object, cells = sketch_cells)
  # query_obj <- subset(object, cells = sketch_cells, invert = T)
  # run WNN
  weighted.nn <- MultiModalNN(object = sketch_obj,
                              query = object,
                              k.nn = ifelse(20 >= k.nn, 20, k.nn),
                              modality.weight = wnn.weight,
                              knn.range = knn.range,
                              verbose = verbose,
                              l2.norm = l2.norm,
                              ...)
  #
  object[[nn.name]] <- weighted.nn
  return (object)
}

#' Find Single-Modality Nearest Neighbors
#'
#' A wrapper function that builds a Neighbor object from a single reduction
#' (e.g. the PLS embedding), mapping every cell to its nearest landmark cells.
#'
#' This is the single-modality counterpart of \code{\link{FindmmNN}}: instead of
#' combining several reductions with WNN weights, it runs a single
#' \code{Seurat::FindNeighbors} call on one reduction, with the landmark (sketched)
#' cells as the reference and all cells as the query. The resulting Neighbor object
#' is a drop-in replacement for the \code{weighted.nn} produced by \code{FindmmNN}
#' and is consumed by \code{\link{GenerateSampleObject}} in the same way.
#'
#' @param object A Seurat object
#' @param sketch.assay The name of the sketched assay of the landmark cells.
#' @param reduction The name of a single dimensional reduction to use (e.g. \code{"proj.pls"}).
#' @param dims The dimensions of the reduction to use. If \code{NULL} (default), all
#'   dimensions are used.
#' @param k.nn The number of nearest neighbors to compute for each cell. Honored as
#'   given (unlike \code{\link{FindmmNN}}, which clamps to a minimum of 20).
#' @param l2.norm Perform L2 normalization on the cell embeddings
#' @param nn.name Name under which the resulting Neighbor object is stored.
#' @param verbose Print progress bars and output
#' @param ... Arguments passed to \code{Seurat::FindNeighbors}
#'
#' @return return a Seurat object that contains a Neighbor object between the landmark cells and all the other cells.
#' @export
#'
#' @importFrom SeuratObject Cells Embeddings
#'
#' @concept sample_level_analysis
FindsmNN <- function(
    object,
    sketch.assay = NULL,
    reduction,
    dims = NULL,
    k.nn = 20,
    l2.norm = TRUE,
    nn.name = "single.nn",
    verbose = TRUE,
    ...
){
  if (k.nn <= 2) {
    stop("k.nn must be a integer larger than 1.")
  }
  # prepare for the objects used to search for NN
  if (is.null(sketch.assay) || !sketch.assay %in% names(object@assays)) {
    stop("sketch.assay is not correctly defined. please check.")
  }
  if (!reduction %in% names(object@reductions)){
    stop("reduction '", reduction, "' is not found in the object. please check.")
  }

  sketch_cells <- Cells(object[[sketch.assay]])

  # pull the embedding directly (avoids subsetting/copying the assays)
  emb <- Embeddings(object, reduction = reduction)
  if (!is.null(dims)) {
    emb <- emb[, dims, drop = FALSE]
  }

  # reference = landmarks, query = all cells (mirrors FindmmNN's query = object)
  single.nn <- FindNeighbors(object = emb[sketch_cells, , drop = FALSE],
                             query = emb,
                             k.param = k.nn,
                             return.neighbor = TRUE,
                             l2.norm = l2.norm,
                             verbose = verbose,
                             ...)

  object[[nn.name]] <- single.nn
  return (object)
}

#' Apply Chi-squared normalization to a Seurat object
#'
#' This function applies Chi-squared normalization to a (sample-level) Seurat object. It calculates the expected density (count of cells) for each
#' landmark-sample pair assuming random distribution, then normalizes observed counts using:
#' (observed - expected) / sqrt(expected). The resulting values show how much each landmark's density deviates
#' from the expected baseline density values. Positive values indicate higher-than-expected density,
#' negative values indicate lower-than-expected density.
#'
#' @param object a count matrix (a landmark-sample density matrix), Assay, StdAssay, or Seurat object.
#'   Values must be raw counts (non-negative whole numbers)
#' @param assay Name of assay to use (for Seurat objects)
#' @param layer Layer to use (for StdAssay objects, and passed on by the Seurat method). It must
#'   resolve to a single layer whose name starts with "counts"; for split layers, run
#'   \code{JoinLayers} first or name one layer
#' @param save Layer to save results to (for StdAssay objects, and passed on by the Seurat method).
#'   Its name cannot start with "counts"
#' @param verbose display progress bar for normalization procedure
#' @param ... Additional arguments
#'
#' @return Normalized object of the same class as input
#' @export
#'
#' @examples
#' \dontrun{
#' # Normalize a sample-level Seurat object
#' sample_obj <- NormalizeChiSquared(sample_obj)
#' }
#'
#' @importFrom SeuratObject LayerData LayerData<-
#'
#' @concept sample_level_analysis
NormalizeChiSquared <- function(object, ...) {
  UseMethod(generic = 'NormalizeChiSquared', object = object)
}

#' @rdname NormalizeChiSquared
#' @method NormalizeChiSquared default
#' @export
NormalizeChiSquared.default <- function(
    object = NULL,
    verbose = TRUE,
    ...
){
  if (verbose) {
    cat("Performing Chi-squared-normalization\n")
  }
  #
  ct_matrix <- object
  # other values would not fail below, they would give plausible but wrong results
  if (any(ct_matrix < 0) || any(ct_matrix != round(x = ct_matrix))) {
    stop("NormalizeChiSquared() requires raw counts (non-negative whole numbers). ",
         "Please check the 'layer' being used.")
  }
  # calculate the row sums and column sums
  col_Sums_ct <- SeuratObject::colSums(ct_matrix)
  row_Sums_ct <- SeuratObject::rowSums(ct_matrix)
  total_ct <- sum(ct_matrix)
  # the expected count matrix
  expect_ct_matrix <- t(outer(col_Sums_ct, row_Sums_ct)/total_ct)
  # calculate the Pearson residual
  relative_matrix <- (ct_matrix - expect_ct_matrix) / sqrt(expect_ct_matrix)
  # there are some 0(NA) in the relative_matrix, need to convert them
  relative_matrix[is.na(relative_matrix)] <- 0

  # a standard matrix object should be fine since the matrix will always be dense
  relative_matrix <- as.matrix(relative_matrix)
  #
  return(relative_matrix)
}

#' @rdname NormalizeChiSquared
#' @method NormalizeChiSquared Assay
#' @export
NormalizeChiSquared.Assay <- function(
    object,
    verbose = TRUE,
    ...
) {
  LayerData(object, layer = "data") <- NormalizeChiSquared(
    object = LayerData(object = object, layer = "counts"),
    verbose = verbose
  )
  return(object)
}

#' @rdname NormalizeChiSquared
#' @method NormalizeChiSquared StdAssay
#' @export
NormalizeChiSquared.StdAssay <- function(
    object,
    layer = 'counts',
    save = 'data',
    verbose = TRUE,
    ...
) {
  if (startsWith(x = save, prefix = "counts")) {
    stop("'save' cannot start with 'counts', as this would overwrite raw counts ",
         "or be mistaken for them, but got: ", save, ".")
  }
  olayer <- unique(x = layer)
  layer <- Layers(object = object, search = olayer)
  if (length(x = layer) == 0) {
    stop("No layer matching '", paste(olayer, collapse = ", "), "' was found.")
  }
  if (!all(startsWith(x = layer, prefix = "counts"))) {
    stop("NormalizeChiSquared() requires a layer whose name starts with 'counts', but got: ",
         paste(layer, collapse = ", "), ".")
  }
  # the expected counts depend on the matrix being normalized, so normalizing
  # layers separately would not match normalizing them together
  if (length(x = layer) > 1) {
    stop("NormalizeChiSquared() works on a single layer, but found: ",
         paste(layer, collapse = ", "), ". Please run JoinLayers() first, or name ",
         "one layer through 'layer'.")
  }
  LayerData(object, layer = save) <- NormalizeChiSquared(
    object = LayerData(object = object, layer = layer),
    verbose = verbose
  )
  return(object)
}

#' @rdname NormalizeChiSquared
#' @method NormalizeChiSquared Seurat
#' @export
NormalizeChiSquared.Seurat <- function(
    object,
    assay = NULL,
    verbose = TRUE,
    ...
) {
  assay <- assay %||% DefaultAssay(object = object)
  assay <- match.arg(arg = assay, choices = Assays(object = object))

  assay.data <- NormalizeChiSquared(
    object = object[[assay]],
    verbose = verbose,
    ...
  )
  object[[assay]] <- assay.data
  object <- LogSeuratCommand(object = object)
  return(object)
}

#' Perform cor.test between each gene and a response variable given a Seurat object
#'
#' This function take a Seurat object's gene expression matrix and a meta-data column (a response variable), and
#' run cor.test between each gene and the response variable and keep track of the top-correlated genes and return
#' that gene list to user. It is intended to give users a more response-informative gene list for downstream supervised
#' dimension reduction analyses like partial least squares.
#'
#' @param object A Seurat object
#' @param assay The name of the assay to retrieve the expression matrix
#' @param layer The name of the layer to use for cor.test
#' @param group.by A metadata column name to group cells by before cor.test. If NULL,
#'   falls back to cor.test without grouping. Default is NULL.
#' @param Y A metadata column name of responses. No cell in 'assay' can have a missing (\code{NA}) response.
#' @param cor.cutoff The correlation coefficient cutoff to select the top genes. Default is NULL.
#' @param top.prop The proportion of genes to be selected as the top genes. Default is 0.01, i.e. the most
#'   correlated 1\% of genes will be selected. When used together with cor.cutoff, genes must pass both filters;
#'   set top.prop to NULL to select every gene that passes cor.cutoff. They cannot be NULL at the same time.
#' @param min.cells the minimum number of cells per group to perform the cor.test. Default is 100. groups with < 100 will be skipped.
#' @param min.pct only test genes that are detected in a minimum fraction of min.pct cells in each group. Default is 0.1
#' @param verbose Print progress bars and output
#' @param unlist whether to return a unlist object (vector) or a named list (each element is a vector of DEGs for a category in group.by)
#' @param ... Additional arguments
#'
#' @return return a list of gene names for the genes that highly correlated with the responses
#' @export
#'
#' @examples
#' \dontrun{
#' # Find genes correlated with pseudotime
#' cor_genes <- QuickCorTest(seurat_obj, Y = "pseudotime")
#' }
#'
#' @importFrom stats model.matrix sd
#' @importFrom SeuratObject FetchData DefaultAssay Assays LayerData Layers CellsByIdentities Idents<- Cells
#'
#' @concept sample_level_analysis
QuickCorTest <- function(
    object = NULL,
    assay = NULL,
    layer = "data",
    group.by = NULL,
    Y = NULL,
    cor.cutoff = NULL,
    top.prop = 0.01,
    min.pct = 0.1,
    min.cells = 100,
    verbose = TRUE,
    unlist = TRUE,
    ...
){
  assay <- assay %||% DefaultAssay(object = object)
  assay <- match.arg(arg = assay, choices = Assays(object = object))
  # fetch the target variable
  if (is.null(x = Y) | length(x = Y) != 1) {
    stop("When setting 'Y', one and only one variable is allowed.")
  } else {
    Y <- FetchData(
      object = object,
      vars = Y,
      cells = Cells(x = object[[assay]]),
      clean = FALSE
    )
    .check_Y_missing(Y = Y, label = "Y")
  }
  for(i in 1:ncol(Y)){
    if(is.character(Y[,i])) Y[,i] <- as.factor(Y[,i])
  }
  # fetch the expression assay
  data.use <- object[[assay]]
  if (length(x = Layers(object = data.use, search = layer)) > 1) {
    stop(layer, " layers are not joined. Please run JoinLayers")
  }
  # get cell IDs within each category of the "group.by"
  if (!is.null(x = group.by) && group.by %in% colnames(object@meta.data)) {
    Idents(object) <- group.by
    cell_celltype_list <- CellsByIdentities(object, cells = Cells(data.use))
  } else {
    cell_celltype_list <- list(all = Cells(data.use))
  }
  # check
  if(is.null(x = cor.cutoff) & is.null(x = top.prop)) stop("Must provide at least one of 'top.prop' and 'cor.cutoff'.")
  #### build each group's response matrix first, so that groups which are too small,
  #### or whose responses turn out to be constant, never reach the expression data
  responses <- lapply(names(cell_celltype_list), function(CELLTYPE) {
    idx <- cell_celltype_list[[CELLTYPE]]
    # Skip if the number of cells is less than min.cells
    if(length(idx) < min.cells){
      if(verbose) message(paste0("Skipping ", CELLTYPE, " since the number of cells is less than ", min.cells))
      return(NULL)  # Return NULL for skipped cell types
    } else {
      if(verbose) message(paste0("Including ", CELLTYPE, " (", length(idx), " cells)"))
    }
    # get Y
    Y2 <- Y[idx, , drop = F]
    Y2 <- model.matrix(~. + 0, data = Y2)
    rm_idx = vector()
    for(i in 1:ncol(Y2)){
      if(length(unique(Y2[, i])) <= 1) rm_idx <- c(rm_idx, i)
    }
    if(length(rm_idx) > 0) Y2 <- Y2[, -rm_idx, drop = F]
    if(ncol(Y2) == 0) return(NULL)
    return(Y2)
  })
  names(responses) <- names(cell_celltype_list)
  tested <- names(responses)[!sapply(responses, is.null)]
  slct_genes_list <- vector(mode = "list", length = length(x = cell_celltype_list))
  names(slct_genes_list) <- names(cell_celltype_list)

  #### read the expression matrix once and derive every group's statistics and
  #### cross-products from that, instead of subsetting cells group by group
  if(length(tested) > 0){
    responses <- responses[tested]
    df_data <- LayerData(object = data.use, assay = assay, layer = layer)
    index <- lapply(cell_celltype_list[tested], match, table = colnames(df_data))
    groups <- factor(rep_len(NA_character_, ncol(df_data)), levels = tested)
    for(CELLTYPE in tested) groups[index[[CELLTYPE]]] <- CELLTYPE
    # cells outside every tested group still need a level to be counted under
    if(anyNA(groups)){
      levels(groups) <- c(tested, ".unused")
      groups[is.na(groups)] <- ".unused"
    }
    # Detection rate, mean and SD of every gene in every group, in a single pass
    if(verbose) message(paste0("Pass 1/2: gene statistics for ", length(tested), " groups over ",
                               ncol(x = df_data), " cells..."))
    group.stats <- GroupStats(x = df_data, groups = groups)
    # Cross-products of every gene against every group's responses, in a second pass
    if(verbose) message("Pass 2/2: gene-response cross-products...")
    xy <- GroupCrossprod(x = df_data, responses = responses, index = index)
    if(verbose) message("Selecting the top-correlated genes...")

    ####  loop through the cell types to get the top DEGs (separately for each column in Y)
    slct_genes_list[tested] <- lapply(tested, function(CELLTYPE) {
      Y2 <- responses[[CELLTYPE]]
      n <- nrow(x = Y2)
      idx_rowCount <- group.stats$nonzero[, CELLTYPE] >= min.pct*n
      features <- rownames(df_data)[idx_rowCount]
      # Calculate correlation more efficiently
      cor_res <- CorFromCrossprod(
        xy = xy[[CELLTYPE]][idx_rowCount, , drop = FALSE],
        n = n,
        x.mean = group.stats$mean[idx_rowCount, CELLTYPE],
        x.sd = group.stats$sd[idx_rowCount, CELLTYPE],
        y = Y2
      )
      # a single response collapses to the named vector that apply() used to return
      if (nrow(x = cor_res) == 1) {
        cor_res <- cor_res[1, ]
      }

      # Handle both vector (continuous Y) and matrix (categorical Y) cases
      if (is.vector(cor_res)) {
        # For continuous variables, cor_res is a vector
        abs_cor <- abs(cor_res)
        if(is.null(x = top.prop)){
          idx_slct <- seq_along(along.with = abs_cor)
        } else {
          n_top <- ceiling(top.prop * length(cor_res))
          idx_slct <- order(abs_cor, decreasing = TRUE)[seq_len(length.out = n_top)]
        }
        if(!is.null(x = cor.cutoff)){
          cutoff_idx <- which(abs_cor >= cor.cutoff)
          idx_slct <- intersect(idx_slct, cutoff_idx)
        }
      } else {
        # For categorical variables, cor_res is a matrix
        if(!is.null(x = top.prop)) n_top <- ceiling(top.prop * ncol(cor_res))
        idx_slct <- apply(abs(cor_res), 1, function(x) {
          if(is.null(x = top.prop)) return(which(x >= cor.cutoff))
          a <- order(x, decreasing = TRUE)[seq_len(length.out = n_top)]
          if(!is.null(x = cor.cutoff)){
            b <- which(x >= cor.cutoff)
            return(intersect(a, b))
          } else {
            return(a)
          }
        })
        idx_slct <- sort(unique(as.vector(unlist(idx_slct))))
      }

      # Return the selected genes
      features[idx_slct]
    })
  }
  # Remove NULL entries (skipped cell types)
  if(isTRUE(x = unlist)){
    slct_genes_list <- slct_genes_list[!sapply(slct_genes_list, is.null)]
    all_genes <- Reduce(union, slct_genes_list)
    return(all_genes)
  } else {
    return(slct_genes_list)
  }
}

# Helper function: nonzero count, mean and standard deviation of every gene within
# every group, as features x groups matrices. BPCells::pseudobulk_matrix() returns
# all three for all groups from a single streamed pass, so an on-disk matrix is read
# once no matter how many groups there are -- subsetting its cells group by group
# would instead cost one random-access pass each.
GroupStats <- function(x, groups) {
  if (inherits(x = x, what = 'IterableMatrix')) {
    pb <- BPCells::pseudobulk_matrix(
      mat = x,
      cell_groups = groups,
      method = c('nonzeros', 'mean', 'variance')
    )
    return(list(nonzero = pb$nonzeros, mean = pb$mean, sd = sqrt(x = pb$variance)))
  }
  n <- tabulate(bin = groups, nbins = nlevels(x = groups))
  indicator <- Matrix::sparseMatrix(
    i = seq_along(along.with = groups),
    j = as.integer(x = groups),
    x = 1,
    dims = c(length(x = groups), nlevels(x = groups)),
    dimnames = list(colnames(x = x), levels(x = groups))
  )
  # `x > 0` builds a logical copy of the whole matrix that the product below then
  # coerces straight back to double. Reusing the sparsity pattern instead skips both
  # steps; it counts stored entries rather than positive ones, which is the same
  # thing for the counts and data layers this runs on
  if (inherits(x = x, what = 'dgCMatrix')) {
    pattern <- x
    pattern@x <- rep.int(x = 1, times = length(x = x@x))
  } else {
    pattern <- x > 0
  }
  group.mean <- t(x = t(x = as.matrix(x = x %*% indicator)) / n)
  group.sumsq <- as.matrix(x = (x ^ 2) %*% indicator)
  group.var <- t(x = (t(x = group.sumsq) - n * t(x = group.mean) ^ 2) / (n - 1))
  return(list(
    nonzero = as.matrix(x = pattern %*% indicator),
    mean = group.mean,
    sd = sqrt(x = pmax(group.var, 0))
  ))
}

# Helper function: sum_i x_gi * y_ik over each group's cells, for every gene g,
# response k and group. Every group's responses are packed side by side into one
# right-hand side, zero outside that group's cells, so a single x %*% y covers all
# of them in one pass. Returns one features x responses matrix per group.
#
# Since each cell carries responses for its own group only, that packed right-hand
# side is mostly structural zeros once there is more than one group. Against a
# sparse x a dense right-hand side therefore costs nnz * total responses, growing
# linearly with the number of groups, while a sparse one stays flat (11.6x quicker
# at 40 groups). BPCells is the other way round -- it answers sparse products with
# a lazy MatrixMultiply but has an eager dense kernel -- so on-disk matrices, and
# dense in-memory ones, keep the dense right-hand side.
GroupCrossprod <- function(x, responses, index, max.bytes = 2.5e8) {
  n <- ncol(x = x)
  widths <- vapply(X = responses, FUN = ncol, FUN.VALUE = integer(length = 1L))
  sparse <- length(x = responses) > 1L && inherits(x = x, what = 'sparseMatrix')
  # a sparse right-hand side holds only the responses themselves, so every group
  # fits in one batch; a dense one is batched to keep its size bounded
  batch.size <- if (sparse) sum(widths) else max(1, floor(x = max.bytes / (8 * n)))
  res <- vector(mode = 'list', length = length(x = responses))
  names(x = res) <- names(x = responses)
  for (batch in split(x = seq_along(along.with = responses), f = ceiling(x = cumsum(x = widths) / batch.size))) {
    offsets <- cumsum(x = c(0, widths[batch]))
    if (isTRUE(x = sparse)) {
      rhs <- Matrix::sparseMatrix(
        i = unlist(x = lapply(
          X = seq_along(along.with = batch),
          FUN = function(i) rep(x = index[[batch[i]]], times = widths[batch[i]])
        )),
        j = unlist(x = lapply(
          X = seq_along(along.with = batch),
          FUN = function(i) rep(
            x = offsets[i] + seq_len(length.out = widths[batch[i]]),
            each = length(x = index[[batch[i]]])
          )
        )),
        x = unlist(x = lapply(X = responses[batch], FUN = as.vector)),
        dims = c(n, sum(widths[batch]))
      )
    } else {
      rhs <- matrix(data = 0, nrow = n, ncol = sum(widths[batch]))
      for (i in seq_along(along.with = batch)) {
        rhs[index[[batch[i]]], offsets[i] + seq_len(length.out = widths[batch[i]])] <- responses[[batch[i]]]
      }
    }
    product <- as.matrix(x = x %*% rhs)
    for (i in seq_along(along.with = batch)) {
      res[[batch[i]]] <- product[, offsets[i] + seq_len(length.out = widths[batch[i]]), drop = FALSE]
    }
  }
  return(res)
}

# Helper function: Pearson correlations from precomputed cross-products. Given
# xy[g, k] = sum_i x_gi * y_ik over the group's n cells, the correlation follows in
# closed form from the gene and response means/SDs, with no further access to the
# expression data. Returns a responses x features matrix, matching the orientation
# apply(x, 1, cor, y) produced.
CorFromCrossprod <- function(xy, n, x.mean, x.sd, y) {
  y.mean <- colMeans(x = y)
  y.sd <- apply(X = y, MARGIN = 2, FUN = sd)
  # cov(x_i, y_j) = (sum_k x_ik * y_kj - n * mean(x_i) * mean(y_j)) / (n - 1)
  covariance <- (xy - n * outer(X = x.mean, Y = y.mean)) / (n - 1)
  res <- t(x = covariance / outer(X = x.sd, Y = y.sd))
  # constant genes have no correlation, as cor() would also report
  res[!is.finite(x = res)] <- NA
  return(res)
}
