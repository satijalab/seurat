#' @include generics.R
#' @importFrom Rcpp evalCpp
#' @importFrom Matrix colSums rowSums colMeans rowMeans
#' @importFrom methods setClass setOldClass setClassUnion slot
#' slot<- setMethod new signature
#' @importClassesFrom Matrix dgCMatrix
#' @useDynLib Seurat
#'
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Class definitions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

setOldClass(Classes = 'package_version')
setClassUnion(name = 'AnyMatrix', c("matrix", "dgCMatrix"))

#' The Assay Class
#'
#' The Assay object is the basic unit of Seurat; each Assay stores raw, normalized, and scaled data
#' as well as cluster information, variable features, and any other assay-specific metadata.
#' Assays should contain single cell expression data such as RNA-seq, protein, or imputed expression
#' data.
#'
#' @slot counts Unnormalized data such as raw counts or TPMs
#' @slot data Normalized expression data
#' @slot scale.data Scaled expression data
#' @slot key Key for the Assay
#' @slot var.features Vector of features exhibiting high variance across single cells
#' @slot meta.features Feature-level metadata
#' @slot misc Utility slot for storing additional data associated with the assay
#'
#' @name Assay
#' @exportClass Assay
#'
Assay <- setClass(
  Class = 'Assay',
  slots = c(
    counts = 'AnyMatrix',
    data = 'AnyMatrix',
    scale.data = 'matrix',
    key = 'character',
    var.features = 'vector',
    meta.features = 'data.frame',
    misc = 'ANY'
  )
)

#' The JackStrawData Class
#'
#' The JackStrawData is used to store the results of a JackStraw computation.
#'
#' @slot empirical.p.values Empirical p-values
#' @slot fake.reduction.scores Fake reduction scores
#' @slot empirical.p.values.full Empirical p-values on full
#' @slot overall.p.values Overall p-values from ScoreJackStraw
#'
#' @name JackStrawData
#' @exportClass JackStrawData
#'
JackStrawData <- setClass(
  Class = "JackStrawData",
  slots = list(
    empirical.p.values = "matrix",
    fake.reduction.scores = "matrix",
    empirical.p.values.full = "matrix",
    overall.p.values = "matrix"
  )
)

#' The Dimmensional Reduction Class
#'
#' The DimReduc object is ...
#'
#' @slot cell.embeddings ...
#' @slot feature.loadings ...
#' @slot feature.loadings.projected ...
#' @slot assay.used ...
#' @slot stdev ...
#' @slot key ...
#' @slot jackstraw ...
#' @slot misc ...
#'
#' @name DimReduc
#' @exportClass DimReduc
#'
DimReduc <- setClass(
  Class = 'DimReduc',
  slots = c(
    cell.embeddings = 'matrix',
    feature.loadings = 'matrix',
    feature.loadings.projected = 'matrix',
    assay.used = 'character',
    stdev = 'numeric',
    key = 'character',
    jackstraw = 'JackStrawData',
    misc = 'list'
  )
)

#' The Graph Class
#'
#' The Graph class simply inherits from dgCMatrix
#'
#' @name Graph
#' @exportClass Graph
#'
Graph <- setClass(
  Class = 'Graph',
  contains = "dgCMatrix"
)

#' The SeuratCommand Class
#'
#' The SeuratCommand is used for logging commands that are run on a SeuratObject. It stores parameters and timestamps
#'
#' @slot name Command name
#' @slot timestamp Timestamp of when command was tun
#' @slot call_string String of the command call
#' @slot params List of parameters used in the command call
#' @name SeuratCommand
#' @exportClass SeuratCommand
#'
SeuratCommand <- setClass(
  Class = 'SeuratCommand',
  slots = c(
    name = 'character',
    time.stamp = 'POSIXct',
    call.string = 'character',
    params = 'ANY'
  )
)

#' The SeuratWorkflow class
#'
#' The SeuratWorkflow class aims to create a Makefile-like approach for Seurat
#' analysis. A full set of parameters can be defined with a single command, and
#' the workflow encodes the pipeline command structure (for example,
#' NormalizeData -> FindVariableFeatures -> ScaleData -> RunPCA -> FindClusters).
#' This leads to very simple workflows being encoded in just a couple commands,
#' without losing the flexibility to run each step for clarity if desired.
#'
#' @slot name Workflow name
#' @slot depends Dependency graph encoding the relationship between commands
#' (for example, RunPCA depends on ScaleData)
#' @slot update Vector specifying for each command, if it needs to be
#' updated/rerun
#' @slot params List of parameters used across the workflow
#' @slot mostRecent Vector specifying for each command, the most recent
#' timestamp. If the current timestamp matches this, should not need to be rerun.
#' @name SeuratWorkflow
#' @exportClass SeuratWorkflow
#'
SeuratWorkflow <- setClass(
  Class = 'SeuratWorkflow',
  slots = c(
    name = 'character',
    depends = 'ANY',
    params = 'ANY',
    mostRecent = 'ANY'
  )
)

#' The Seurat Class
#'
#' The Seurat object is ...
#'
#' @slot assays A list of assays for this project
#' @slot meta.data Contains meta-information about each cell, starting with number of genes detected (nGene)
#' and the original identity class (orig.ident); more information is added using \code{AddMetaData}
#' @slot active.assay ...
#' @slot active.ident The active cluster identity for this Seurat object
#' @slot graphs ...
#' @slot neighbors ...
#' @slot reductions A list of dimmensional reduction objects for this Assay
#' @slot project.name ...
#' @slot calc.params A list of calculation parameters performed on this Seurat object
#' @slot misc A list of miscellaneous information
#' @slot version Version of Seurat this object was built under
#'
#' @name Seurat
#' @exportClass Seurat
#'
Seurat <- setClass(
  Class = 'Seurat',
  slots = c(
    assays = 'list',
    meta.data = 'data.frame',
    active.assay = 'character',
    active.ident = 'factor',
    graphs = 'list',
    neighbors = 'list',
    reductions = 'list',
    project.name = 'character',
    calc.params = 'list',
    misc = 'list',
    version = 'package_version',
    commands = 'list',
    workflows = 'list',
    tools = 'list'
  )
)

# The Seurat Class
#
# The Seurat object is the center of each single cell analysis. It stores all information
# associated with the dataset, including data, annotations, analyes, etc. All that is needed
# to construct a Seurat object is an expression matrix (rows are genes, columns are cells), which
# should be log-scale
#
# Each Seurat object has a number of slots which store information. Key slots to access
# are listed below.
#
# @slot raw.data The raw project data
# @slot data The normalized expression matrix (log-scale)
# @slot scale.data scaled (default is z-scoring each gene) expression matrix; used for dimmensional reduction and heatmap visualization
# @slot var.genes Vector of genes exhibiting high variance across single cells
# @slot is.expr Expression threshold to determine if a gene is expressed (0 by default)
# @slot ident THe 'identity class' for each cell
# @slot meta.data Contains meta-information about each cell, starting with number of genes detected (nGene)
# and the original identity class (orig.ident); more information is added using \code{AddMetaData}
# @slot project.name Name of hte project (for record keeping)
# @slot dr List of stored dimmensional reductions; named by technique
# @slot assay List of additional assays for multimodal analysis; named by technique
# @slot hvg.info The output of the mean/variability analysis for all genes
# @slot imputed Matrix of imputed gene scores
# @slot cell.names Names of all single cells (column names of the expression matrix)
# @slot cluster.tree List where the first element is a phylo object containing the phylogenetic tree relating different identity classes
# @slot snn Spare matrix object representation of the SNN graph
# @slot calc.params Named list to store all calculation-related parameter choices
# @slot kmeans Stores output of gene-based clustering from \code{DoKMeans}
# @slot spatial Stores internal data and calculations for spatial mapping of single cells
# @slot misc Miscellaneous spot to store any data alongisde the object (for example, gene lists)
# @slot version Version of package used in object creation
#
# @name seurat
# @rdname seurat
# @aliases seurat-class
#' @exportClass seurat
#
seurat <- setClass(
  "seurat",
  slots = c(
    raw.data = "ANY",
    data = "ANY",
    scale.data = "ANY",
    var.genes = "vector",
    is.expr = "numeric",
    ident = "factor",
    meta.data = "data.frame",
    project.name = "character",
    dr = "list",
    assay = "list",
    hvg.info = "data.frame",
    imputed = "data.frame",
    cell.names = "vector",
    cluster.tree = "list",
    snn = "dgCMatrix",
    calc.params = "list",
    kmeans = "ANY",
    spatial = "ANY",
    misc = "ANY",
    version = "ANY"
  )
)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Add Metadata
#'
#' Adds additional data for single cells to the Seurat object. Can be any piece
#' of information associated with a cell (examples include read depth,
#' alignment rate, experimental batch, or subpopulation identity). The
#' advantage of adding it to the Seurat object is so that it can be
#' analyzed/visualized using FetchData, VlnPlot, GenePlot, SubsetData, etc.
#'
#' @param object Seurat object
#' @param metadata Data frame where the row names are cell names (note : these
#' must correspond exactly to the items in object@@cell.names), and the columns
#' are additional metadata items.
#' @param col.name Name for metadata if passing in single vector of information
#'
#' @return Seurat object where the additional metadata has been added as
#' columns in object@@meta.data
#'
#' @export
#'
#' @examples
#' cluster_letters <- LETTERS[Idents(object = pbmc_small)]
#' names(cluster_letters) <- colnames(x = pbmc_small)
#' pbmc_small <- AddMetaData(
#'   object = pbmc_small,
#'   metadata = cluster_letters,
#'   col.name = 'letter.idents'
#' )
#' head(x = pbmc_small[[]])
#'
AddMetaData <- function(object, metadata, col.name = NULL) {
  if (typeof(x = metadata) != "list") {
    metadata <- as.data.frame(x = metadata)
    if (is.null(x = col.name)) {
      stop("Please provide a name for provided metadata")
    }
    colnames(x = metadata) <- col.name
  }
  cols.add <- colnames(x = metadata)
  #meta.add <- metadata[rownames(x = object@meta.data), cols.add]
  meta.order <- match(x = rownames(x = object[[]]), rownames(x = metadata))
  meta.add <- metadata[meta.order, ]
  if (all(is.null(x = meta.add))) {
    stop("Metadata provided doesn't match the cells in this object")
  }
  slot(object = object, name = "meta.data")[, cols.add] <- meta.add
  return(object)
}

#' Checks if a workflow is defined for a Seurat object
#'
#' Checks if a workflow is defined for a Seurat object
#'
#' @param object Seurat object
#' @param workflow.name Workflow name, should already be initialized using InitializeWorkflow#'
#' @return TRUE if workflow is defined. STOP otherwise
#'
#' @export
#'
#' @examples
#' \dontrun{
#' CheckWorkflow(pbmc_small, "cluster")
#' }
#'
CheckWorkflow <- function(object, workflow.name) {
  # Check if workflow is there
  workflow.present <- FALSE
  if (workflow.name %in% names(x = object)) {
    if (class(x = object[[workflow.name]])[[1]] == "SeuratWorkflow") {
      workflow.present <- TRUE
    }
  }
  if (!workflow.present) {
    stop("Workflow not present, initialize first.")
  }
  return(TRUE)
}

#' Check if workflow command needs update
#'
#' Compares the stored timestamp with the most recently recorded timestamp to see if a dependency has been updated
#'
#' @param object Seurat object
#' @param workflow.name Workflow name, should already be initialized using InitializeWorkflow
#' @param command.name Name of the command to check
#'
#' @return Returns TRUE if the dependency has changed (or has not been run), and an update is needed. FALSE otherwise
#'
#' @export
#'
#' @examples
#' \dontrun{
#' CheckWorkflowUpdate(object = pbmc_small,workflow.name = "cluster", command.name = "ScaleData")
#' }
#'
CheckWorkflowUpdate <- function(object, workflow.name, command.name) {
  CheckWorkflow(object = object, workflow.name = workflow.name)

  # According to the workflow, the most recent update
  mostRecent <- slot(object[[workflow.name]],"mostRecent")
  workflow.timestamp <- mostRecent[command.name]
  seurat.timestamp <- as.POSIXct("1900-01-01");

  #means Seurat command has never been run in the workflow
  if (workflow.timestamp==seurat.timestamp) {
    return(TRUE)
  }
  # According to SeuratCommand, the most recent update
  # go to workflow to look up assay and DR
  params <- slot(object = object[[workflow.name]], name = "params")
  assay <- params$global$assay
  assay <- assay %iff% params[[command.name]]$assay
  assay <- assay %||% DefaultAssay(object)
  reduction <- params$global$reduction
  reduction <- reduction %iff% params[[command.name]]$reduction
  reduction <- reduction %||% formals(fun = paste0(command.name, ".Seurat"))$reduction

  command.name <- paste0(command.name, ".", assay, ".", reduction)
  command.name <- sub(pattern = "[\\.]+$", replacement = "", x = command.name, perl = TRUE)
  command.name <- sub(pattern = "\\.\\.", replacement = "\\.", x = command.name, perl = TRUE)

  if (length(x = command.name)==1) {
    seurat.timestamp <- slot(object = object[[command.name]], name = "time.stamp")
  }
  if (seurat.timestamp == workflow.timestamp) {
    return(FALSE)
  }
  return(TRUE)
}

#' Create an Assay object
#'
#' Create an Assay object from a feature (e.g. gene) expression matrix. The expected format of the
#' input matrix is features x cells.
#'
#' Non-unique cell or feature names are not allowed. Please make unique before calling this function.
#'
#' @param counts Unnormalized data such as raw counts or TPMs
#' @param data Prenormalized data; if provided, donot pass \code{counts}
#' @param min.cells Include features detected in at least this many cells. Will subset the counts
#' matrix as well. To reintroduce excluded features, create a new object with a lower cutoff.
#' @param min.features Include cells where at least this many features are detected.
#'
#' @importFrom methods as
#' @importFrom Matrix colSums
#'
#' @export
#'
#' @examples
#' pbmc_raw <- read.table(
#'   file = system.file('extdata', 'pbmc_raw.txt', package = 'Seurat'),
#'   as.is = TRUE
#' )
#' pbmc_rna <- CreateAssayObject(counts = pbmc_raw)
#' pbmc_rna
#'
CreateAssayObject <- function(
  counts,
  data,
  min.cells = 0,
  min.features = 0
) {
  if (missing(x = counts) && missing(x = data)) {
    stop("Must provide either 'counts' or 'data'")
  } else if (!missing(x = counts) && !missing(x = data)) {
    stop("Either 'counts' or 'data' must be missing; both cannot be provided")
  } else if (!missing(x = counts)) {
    # check that dimnames of input counts are unique
    if (anyDuplicated(rownames(x = counts))) {
      stop("Non-unique features (rownames) present in the input matrix")
    }
    if (anyDuplicated(colnames(x = counts))) {
      stop("Non-unique cell names (colnames) present in the input matrix")
    }
    if (!inherits(x = counts, what = 'dgCMatrix')) {
      counts <- as(object = as.matrix(x = counts), Class = 'dgCMatrix')
    }
    # Filter based on min.features
    nfeatures <- Matrix::colSums(x = counts > 0)
    counts <- counts[, which(x = nfeatures > min.features)]
    # filter genes on the number of cells expressing
    if (min.cells > 0) {
      num.cells <- rowSums(x = counts > 0)
      counts <- counts[which(x = num.cells >= min.cells), ]
    }
    data <- counts
  } else if (!missing(x = data)) {
    # check that dimnames of input counts are unique
    if (anyDuplicated(rownames(x = data))) {
      stop("Non-unique features (rownames) present in the input matrix")
    }
    if (anyDuplicated(colnames(x = data))) {
      stop("Non-unique cell names (colnames) present in the input matrix")
    }
    counts <- new(Class = 'matrix')
  }
  # Initialize meta.features
  init.meta.features <- data.frame(row.names = rownames(x = data))
  assay <- new(
    Class = 'Assay',
    counts = counts,
    data = data,
    scale.data = new(Class = 'matrix'),
    meta.features = init.meta.features
  )
  return(assay)
}

#' Create a DimReduc object
#'
#' @param embeddings ...
#' @param loadings ...
#' @param projected ...
#' @param assay ...
#' @param stdev ...
#' @param key ...
#' @param jackstraw ...
#' @param misc ...
#' @param ... Ignored for now
#'
#' @export
#'
CreateDimReducObject <- function(
  embeddings = new(Class = 'matrix'),
  loadings = new(Class = 'matrix'),
  projected = new(Class = 'matrix'),
  assay = NULL,
  stdev = numeric(),
  key = NULL,
  jackstraw = NULL,
  misc = list(),
  ...
) {
  # if (is.null(assay.used)) {
  #   stop("Please specify the assay that was used to construct the reduction")
  # }
  if (is.null(key)) {
    stop("Please specify a key for the DimReduc object")
  }
  jackstraw <- jackstraw %||% new(Class = 'JackStrawData')
  dim.reduc <- new(
    Class = 'DimReduc',
    cell.embeddings = embeddings,
    feature.loadings = loadings,
    feature.loadings.projected = projected,
    assay.used = assay,
    stdev = stdev,
    key = key,
    jackstraw = jackstraw,
    misc = misc
  )
  return(dim.reduc)
}

#' Create a Seurat object
#'
#' Create a Seurat object from a feature (e.g. gene) expression matrix. The expected format of the
#' input matrix is features x cells.
#'
#'
#' Note: In previous versions (<3.0), this function also accepted a parameter to set the expression
#' threshold for a 'detected' feature (gene). This functionality has been removed to simplify the
#' initialization process/assumptions. If you would still like to impose this threshold for your
#' particular dataset, simply filter the input expression matrix before calling this function.
#'
#' @inheritParams CreateAssayObject
#' @param project Sets the project name for the Seurat object.
#' @param assay Name of the assay corresponding to the initial input data.
#' @param names.field For the initial identity class for each cell, choose this field from the
#' cell's name. E.g. If your cells are named as BARCODE_CLUSTER_CELLTYPE in the input matrix, set
#' names.field to 3 to set the initial identities to CELLTYPE.
#' @param names.delim For the initial identity class for each cell, choose this delimiter from the
#' cell's column name. E.g. If you cells are named as BARCODE-CLUSTER-CELLTYPE, set this to "-" to
#' separate the cell name into it's component parts for picking the relevant field.
#' @param meta.data Additional metadata to add to the Seurat object. Should be a data frame where
#' the rows are cell names, and the columns are additional metadata fields.
#'
#' @importFrom utils packageVersion
#' @importFrom Matrix colSums
#' @export
#'
#' @examples
#' pbmc_raw <- read.table(
#'   file = system.file('extdata', 'pbmc_raw.txt', package = 'Seurat'),
#'   as.is = TRUE
#' )
#' pbmc_small <- CreateSeuratObject(counts = pbmc_raw)
#' pbmc_small
#'
CreateSeuratObject <- function(
  counts,
  project = 'SeuratProject',
  assay = 'RNA',
  min.cells = 0,
  min.features = 0,
  names.field = 1,
  names.delim = "_",
  meta.data = NULL
) {
  assay.data <- CreateAssayObject(
    counts = counts,
    min.cells = min.cells,
    min.features = min.features
  )
  Key(object = assay.data) <- paste0(tolower(x = assay), '_')
  assay.list <- list(assay.data)
  names(x = assay.list) <- assay
  init.meta.data <- data.frame(row.names = colnames(x = assay.list[[assay]]))
  # Set idents
  idents <- factor(x = unlist(x = lapply(
    X = colnames(x = counts),
    FUN = ExtractField,
    field = names.field,
    delim = names.delim
  )))
  # if there are more than 100 idents, set all idents to ... name
  ident.levels <- length(x = unique(x = idents))
  if (ident.levels > 100 || ident.levels == 0 || ident.levels == length(x = idents)) {
    idents <- rep.int(x = factor(x = project), times = ncol(x = counts))
  }
  names(x = idents) <- colnames(x = counts)
  object <- new(
    Class = 'Seurat',
    assays = assay.list,
    meta.data = init.meta.data,
    active.assay = assay,
    active.ident = idents,
    project.name = project,
    version = packageVersion(pkg = 'Seurat')
  )
  object[['orig.ident']] <- idents
  # Calculate nCount and nFeature
  filtered.counts <- GetAssayData(object = object, assay = assay, slot = 'counts')
  object[[paste('nCount', assay, sep = '_')]] <- Matrix::colSums(filtered.counts)
  object[[paste('nFeature', assay, sep = '_')]] <- Matrix::colSums(filtered.counts > 0)
  if (!is.null(meta.data)) {
    object <- AddMetaData(object = object, metadata = meta.data)
  }
  return(object)
}

#' Access cellular data
#'
#' Retreives data (feature expression, PCA scores, metrics, etc.) for a set
#' of cells in a Seurat object
#'
#' @param object Seurat object
#' @param vars List of all variables to fetch
#' @param cells Cells to collect data for (default is all cells)
#' @param slot Slot to pull feature data for
#'
#' @return A data frame with cells as rows and cellular data as columns
#'
#' @export
#'
#' @examples
#' pc1 <- FetchData(object = pbmc_small, vars = 'PC1')
#' head(x = pc1)
#'
FetchData <- function(object, vars, cells = NULL, slot = 'data') {
  cells <- cells %||% colnames(x = object)
  if (is.numeric(x = cells)) {
    cells <- colnames(x = object)[cells]
  }
  objects.use <- FilterObjects(object = object)
  object.keys <- sapply(X = objects.use, FUN = function(i) {return(Key(object[[i]]))})
  keyed.vars <- lapply(
    X = object.keys,
    FUN = function(key) {
      if (length(x = key) == 0) {
        return(integer(length = 0L))
      }
      return(grep(pattern = paste0('^', key), x = vars))
    }
  )
  keyed.vars <- Filter(f = length, x = keyed.vars)
  data.fetched <- lapply(
    X = names(x = keyed.vars),
    FUN = function(x) {
      vars.use <- vars[keyed.vars[[x]]]
      key.use <- object.keys[x]
      data.return <- switch(
        EXPR = class(x = object[[x]]),
        'DimReduc' = {
          vars.use <- grep(
            pattern = paste0('^', key.use, '[[:digit:]]+$'),
            x = vars.use,
            value = TRUE
          )
          if (length(x = vars.use) > 0) {
            object[[x]][[cells, vars.use, drop = FALSE]]
          } else {
            NULL
          }
        },
        'Assay' = {
          vars.use <- gsub(pattern = paste0('^', key.use), replacement = '', x = vars.use)
          data.vars <- t(x = as.matrix(x = GetAssayData(
            object = object,
            slot = slot,
            assay = x
          )[vars.use, cells, drop = FALSE]))
          colnames(x = data.vars) <- paste0(key.use, vars.use)
          data.vars
        }
      )
      data.return <- as.list(x = as.data.frame(x = data.return))
      return(data.return)
    }
  )
  data.fetched <- unlist(x = data.fetched, recursive = FALSE)
  meta.vars <- vars[vars %in% colnames(x = object[[]])]
  data.fetched <- c(data.fetched, object[[meta.vars]][cells, , drop = FALSE])
  default.vars <- vars[vars %in% rownames(x = object)]
  data.fetched <- c(
    data.fetched,
    as.data.frame(x = t(x = as.matrix(x = GetAssayData(
      object = object,
      slot = slot
    )[default.vars, cells, drop = FALSE])))
  )
  fetched <- names(x = data.fetched)
  vars.missing <- setdiff(x = vars, y = fetched)
  m2 <- if (length(x = vars.missing) > 10) {
    paste0(' (10 out of ', length(x = vars.missing), ' shown)')
  } else {
    ''
  }
  if (length(x = vars.missing) == length(x = vars)) {
    stop(
      "None of the requested variables were found",
      m2,
      ': ',
      paste(head(x = vars.missing, n = 10L), collapse = ', ')
    )
  } else if (length(x = vars.missing) > 0) {
    warning(
      "The following requested variables were not found",
      m2,
      ': ',
      paste(head(x = vars.missing, n = 10L), collapse = ', ')
    )
  }
  data.fetched <- as.data.frame(
    x = data.fetched,
    row.names = cells,
    stringsAsFactors = FALSE
  )
  data.order <- na.omit(object = pmatch(
    x = vars,
    table = fetched
  ))
  if (length(x = data.order) > 1) {
    data.fetched <- data.fetched[, data.order]
  }
  colnames(x = data.fetched) <- vars[vars %in% fetched]
  return(data.fetched)
}

#' Initialize a Seurat workflow
#'
#' Reads dependencies from a file and initializes a Seurat workflow
#'
#' @param object Seurat object
#' @param file Ini configuration file. See cluster.workflow.ini for an example.
#'
#' @return Object with modified workflows
#'
#' @importFrom ini read.ini
#' @export
#'
#' @examples
#' \dontrun{
#' pbmc_small <- InitializeWorkflow(object = pbmc_small, file = 'workflows/cluster.workflow.txt')
#' }
#'
InitializeWorkflow <- function(object, file) {
  if (!file.exists(... = file)) {
    stop("Provided workflow file does not exist.")
  }
  config <- read.ini(filepath = file)
  ValidateWorkflowFile(config = config)
  workflow.name <- gsub(
    pattern = ".workflow.ini",
    replacement = "",
    x = basename(path = file)
  )
  depend.fxns <- unlist(x = strsplit(
    x = unname(obj = unlist(x = config$dependencies)),
    split = ","
  ))
  fxns <- union(x = depend.fxns, y = names(x = config$dependencies))
  depends <- matrix(nrow = length(x = fxns), ncol = length(x = fxns))
  rownames(x = depends) <- colnames(x = depends) <- fxns
  for(cmd in 1:length(x = config$dependencies)) {
    cmd.name <- names(x = config$dependencies[cmd])
    cmd.vals <- unlist(x = strsplit(x = config$dependencies[[cmd]], split = ","))
    for(cv in cmd.vals) {
      depends[cmd.name, cv] <- 1
    }
  }
  mostRecent <- rep(x = as.POSIXct(x = "1900-01-01"), length(x = fxns))
  names(x = mostRecent) <- fxns
  for (mr in names(x = mostRecent)) {
    assay <- config$global$assay
    assay <- assay %iff% config[mr]$assay
    assay <- assay %||% DefaultAssay(object = object)
    reduction <- config$global$reduction
    reduction <- config[mr]$reduction
    reduction <- reduction %||% formals(fun = paste0(mr, ".Seurat"))$reduction
    command.name <- paste0(mr, ".", assay, ".", reduction)
    command.name <- sub(
      pattern = "[\\.]+$",
      replacement = "",
      x = command.name,
      perl = TRUE
    )
    command.name <- sub(
      pattern = "\\.\\.",
      replacement = "\\.",
      x = command.name,
      perl = TRUE
    )
    if (command.name %in% names(x = object)) {
      seurat.timestamp <- slot(object = object[[command.name]], name = "time.stamp")
      mostRecent[mr] <- seurat.timestamp
    }
  }
  params <- list()
  if (!is.null(x = config$global)) {
    params[["global"]] <- config$global
    for (p in 1:length(x = params$global)){
      params$global[names(x = params$global[p])] <- ToNumeric(x = params$global[[p]])
    }
  }
  # set fxn specific params
  fxn.param.names <- setdiff(
    x = names(x = config),
    y = c("dependencies", "global")
  )
  if (length(x = fxn.param.names) > 0) {
    for(i in 1:length(x = fxn.param.names)) {
      params[fxn.param.names[i]] <- config[fxn.param.names[i]]
      for (p in 1:length(x = params[[fxn.param.names[i]]])) {
        params[[fxn.param.names[i]]][[p]] <- ToNumeric(x = params[[fxn.param.names[i]]][[p]])
      }
    }
  }
  seurat.workflow <- new(
    Class = 'SeuratWorkflow',
    name = workflow.name,
    depends = depends,
    params = params,
    mostRecent = mostRecent
  )
  object[[workflow.name]] <- seurat.workflow
  return(object)
}

#' Output individual function calls to recreate workflow
#'
#' Output all commands to reproduce your analysis without shortcuts. Should enhance reproducibility, but can be confused by custom modifcations, usage of SubsetData, etc.
#' We hope this will be very useful, but use with care and verify that it does indeed reproduce your work.
#'
#' @param object Seurat object
#' @param workflow.name Workflow name, should already be initialized using InitializeWorkflow
#' @param command.name Name of the command at the end of the workflow
#' @param depth depth of the recursive call. Only depth 1 outputs the parameters
#'
#' @return Seurat object with updated workflow
#'
#' @export
#'
#' @examples
#' \dontrun{
#' RecreateWorkflows(object = pbmc_small,workflow.name = "cluster", command.name = "FindClusters")
#' }
#'
RecreateWorkflows <- function(object, workflow.name, command.name,depth=1) {
  CheckWorkflow(object = object, workflow.name = workflow.name)
  depends <- slot(object = object[[workflow.name]],name = "depends")
  prereq.commands <- colnames(depends)[which(depends[command.name,]==1)]
  if (depth == 1) {
    message(paste0("\tNeed to output SetParams"))
  }
  for(i in prereq.commands) {
    RecreateWorkflows(object = object,workflow.name = workflow.name,command.name = i,depth = depth + 1)
  }
  #TODO deal with Assay better
  command.name <- intersect(c(command.name, paste0(command.name,".",DefaultAssay(object))), names(object))
  if (length(x = command.name)==1) {
    call.string <- slot(object[[command.name]],"call.string")
    #browser()
    message(paste0("\t",call.string))
  }
}

#' Set Seurat workflow parameters
#'
#' Sets parameters for a workflow
#'
#' @param object Seurat object
#' @param workflow.name Workflow name, should already be initialized using
#' InitializeWorkflow
#' @param fxn.param.names Name of the function and parameter to set (formatted as
#' FXNNAME_PARAM). Can take a vector to set multiple functions
#' @param fxn.param.values Value of the parameter to set. Should be of equal length
#' to fxn.param.names
#' @param \dots Global parameters to set, will be fed into workflow Seurat functions.
#' Parameters ending with a "." will populate all similar variable names (i.e.
#'  setting dims. will set both dims.compute, and dims.cluster)
#'
#' @return Object with modified workflows
#'
#' @export
#'
#' @examples
#' \dontrun{
#' pbmc_small <- SetWorkflowParams(object = pbmc_small, seed.use = 31, dims. = 20)
#' }
#'
SetWorkflowParams <- function(
  object,
  workflow.name = NULL,
  ...,
  fxn.param.names = NULL,
  fxn.param.values = NULL
) {
  CheckWorkflow(object = object, workflow.name = workflow.name)
  if (length(x = fxn.param.names) != length(x = fxn.param.values)) {
    stop("length of fxn.parameter.names needs to equal length of fxn.parameter.values")
  }
  params <- slot(object = object[[workflow.name]], name = "params")
  # set global params
  global.params <- list(...)
  for (gp in 1:length(global.params)) {
    params[["global"]][names(global.params[gp])] <- global.params[gp]
  }
  # set fxn specific params
  if (!is.null(fxn.param.names)) {
    for (i in 1:length(x = fxn.param.names)) {
      fxn <- unlist(strsplit(x = fxn.param.names[i], split = "_"))
      params[[fxn[1]]][fxn[2]] <- fxn.param.values[i]
    }
  }
  slot(object = object[[workflow.name]], name = "params") <- params
  return(object)
}

#' Splits object into a list of subsetted objects.
#'
#' Splits object based on a single attribute into a list of subsetted objects,
#' one for each level of the attribute. For example, useful for taking an object
#' that contains cells from many patients, and subdividing it into
#' patient-specific objects.
#'
#' @param object Seurat object
#' @param split.by Attribute for splitting. Default is "ident". Currently
#' only supported for class-level (i.e. non-quantitative) attributes.
#' @param ... Ignored
#'
#' @return A named list of Seurat objects, each containing a subset of cells
#' from the original object.
#'
#' @export
#'
#' @examples
#' # Assign the test object a three level attribute
#' groups <- sample(c("group1", "group2", "group3"), size = 80, replace = TRUE)
#' names(groups) <- colnames(pbmc_small)
#' pbmc_small <- AddMetaData(object = pbmc_small, metadata = groups, col.name = "group")
#' obj.list <- SplitObject(pbmc_small, split.by = "group")
#'
SplitObject <- function(object, split.by = "ident", ...) {
  if (split.by == 'ident') {
    groupings <- Idents(object = object)
  } else {
    groupings <- FetchData(object = object, vars = split.by)[, 1]
  }
  groupings <- unique(x = as.character(x = groupings))
  obj.list <- list()
  for (i in groupings) {
    if (split.by == "ident") {
      obj.list[[i]] <- subset(x = object, select = i)
    }
    else {
      cells <- which(x = object[[split.by, drop = TRUE]] == i)
      cells <- colnames(x = object)[cells]
      obj.list[[i]] <- subset(x = object, select = cells)
    }
  }
  return(obj.list)
}

#' Find features with highest scores for a given dimensional reduction technique
#'
#' Return a list of features with the strongest contribution to a set of components
#'
#' @param object DimReduc object
#' @param dim.use Dimension to use
#' @param nfeatures Number of features to return
#' @param projected Use the full PCA (projected PCA). Default i s FALSE
#' @param balanced Return an equal number of features with both + and - scores.
#'
#' @return Returns a vector of features
#'
#' @export
#'
#' @examples
#' pbmc_small
#' TopFeatures(object = pbmc_small[["pca"]], dim = 1)
#' # After projection:
#' TopFeatures(object = pbmc_small[["pca"]], dim = 1,  projected = TRUE)
#'
TopFeatures <- function(
  object,
  dim = 1,
  nfeatures = 20,
  projected = FALSE,
  balanced = FALSE
) {
  loadings <- Loadings(object = object, projected = projected)[, dim, drop = FALSE]
  return(Top(
    data = loadings,
    num = nfeatures,
    balanced = balanced
  ))
}

#' Find cells with highest scores for a given dimensional reduction technique
#'
#' Return a list of genes with the strongest contribution to a set of components
#'
#' @param object DimReduc object
#' @param dim.use Dimension to use
#' @param num.cells Number of cells to return
#' @param balanced Return an equal number of cells with both + and - scores.
#'
#' @return Returns a vector of cells
#'
#' @export
#'
#' @examples
#' pbmc_small
#' head(TopCells(object = pbmc_small[["pca"]]))
#' # Can specify which dimension and how many cells to return
#' TopCells(object = pbmc_small[["pca"]], dim = 2, ncells = 5)
#'
TopCells <- function(
  object,
  dim = 1,
  ncells = 20,
  balanced = FALSE
) {
  embeddings <- Embeddings(object = object)[, dim, drop = FALSE]
  return(Top(
    data = embeddings,
    num = ncells,
    balanced = balanced
  ))
}

#' Updates workflow timestamps
#'
#' Like the touch command in linux. Updates a workflow command's timestamp, and its dependencies
#'
#' @param object Seurat object
#' @param workflow.name Workflow name, should already be initialized using InitializeWorkflow
#' @param command.name Name of the command to touch
#' @param time.stamp Timestamp to assign
#'
#' @return Seurat object with updated workflow
#'
#' @export
#'
#' @examples
#' \dontrun{
#' TouchWorkflow(object = pbmc_small,workflow.name = "cluster", command.name = "ScaleData")
#' }
#'
TouchWorkflow <- function(object, workflow.name, command.name, time.stamp = Sys.time()) {
  CheckWorkflow(object = object, workflow.name = workflow.name)
  #Now update all dependencies, recursively
  depends <- slot(object = object[[workflow.name]],name = "depends")
  depend.commands <- colnames(depends)[which(depends[,command.name]==1)]
  mostRecent <- slot(object[[workflow.name]],"mostRecent")
  mostRecent[command.name] <- time.stamp
  slot(object[[workflow.name]],"mostRecent") <- mostRecent
  for(i in depend.commands) {
    object <- TouchWorkflow(object,workflow.name = workflow.name,command.name = i,time.stamp = time.stamp)
  }
  return(object)
}

#' Update old Seurat object to accomodate new features
#'
#' Updates Seurat objects to new structure for storing data/calculations.
#'
#' @param object Seurat object
#'
#' @return Returns a Seurat object compatible with latest changes
#'
#' @importFrom utils packageVersion
#' @importFrom methods .hasSlot new slotNames as
#'
#' @export
#'
#' @examples
#' \dontrun{
#' updated_seurat_object = UpdateSeuratObject(object = old_seurat_object)
#' }
#'
UpdateSeuratObject <- function(object) {
  if (.hasSlot(object, "version")) {
    if (package_version(x = object@version) >= package_version(x = "3.0.0")) {
      message("Object representation is consistent with the most current Seurat version")
      return(object)
    } else if (package_version(x = object@version) >= package_version(x = "2.0.0")) {
      seurat.version <- packageVersion(pkg = "Seurat")
      new.assay <- UpdateAssay(old.assay = object, assay = "RNA")
      assay.list <- list(new.assay)
      names(assay.list) <- "RNA"
      for(i in names(object@assay)) {
        assay.list[[i]] <- UpdateAssay(object@assay[[i]], assay = i)
      }
      new.dr <- UpdateDimReduction(old.dr = object@dr, assay = "RNA")
      new.object <- new(
        Class = "Seurat",
        version = seurat.version,
        assays = assay.list,
        active.assay = "RNA",
        project.name = object@project.name,
        calc.params = object@calc.params,
        misc = object@misc %||% list(),
        active.ident = object@ident,
        reductions = new.dr,
        meta.data = object@meta.data,
        tools = list()
      )
      return(new.object)
    }
  }
  stop("Cannot convert version <2")
}

#' Output status of each command in the workflow
#'
#' For each command in the workflow, indicate whether it is up-to-date.
#'
#' @param object Seurat object
#' @param workflow.name Workflow name, should already be initialized using InitializeWorkflow
#' @param command.name Name of the command at the end of the workflow
#'
#' @return Seurat object with updated workflow
#'
#' @export
#'
#' @examples
#' \dontrun{
#' WorkflowStatus(object = pbmc_small, workflow.name = "cluster")
#' }
#'
WorkflowStatus <- function(object, workflow.name, command.name) {
  CheckWorkflow(object = object, workflow.name = workflow.name)
  message(paste("Status  for", workflow.name, "workflow"))
  depends <- slot(object = object[[workflow.name]], name = "depends")
  all.cmds <- rownames(x = depends)
  for (i in all.cmds) {
    is.updated <- !CheckWorkflowUpdate(
      object = object,
      workflow.name = workflow.name,
      command.name = i
    )
    if (is.updated) {
      message(paste0("\t", i, " up to date"))
    } else {
      message(paste0("\t\t", i, " is out of date"))
    }
  }
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Methods for Seurat-defined generics
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @rdname Convert
#' @export
#' @method as.SingleCellExperiment seurat
#'
as.SingleCellExperiment.seurat <- function(from) {
  return(Convert(from = from, to = 'sce'))
}

#' @rdname Convert
#' @export
#' @method as.seurat SingleCellExperiment
#'
as.seurat.SingleCellExperiment <- function(from) {
  return(Convert(from = from, to = 'seurat'))
}

#' @describeIn Command Get the SeuratCommands
#' @export
#' @method Command Seurat
#'
Command.Seurat <- function(object, command, value = NULL) {
  commands <- slot(object = object, name = "commands")
  if (is.null(x = commands[[command]])) {
    stop(command, " has not been run or is not a valid command.")
  }
  command <- commands[[command]]
  if (is.null(x = value)) {
    return(command)
  }
  params <- slot(object = command, name = "params")
  if (!value %in% names(x = params)) {
    stop(value, " is not a valid parameter for ", slot(object = command, name = "name"))
  }
  return(params[[value]])
}

#' @param X.slot Seurat slot to transfer anndata X into. Default is scale.data
#' @param raw.slot Seurat slot to transfer anndata raw into. Default is data
#' @describeIn Convert from Anndata file to a Seurat object
#' @importFrom reticulate py_to_r
#' @export
#' @method Convert anndata.base.AnnData
#'
Convert.anndata.base.AnnData <- function(
  from,
  to,
  X.slot = "scale.data",
  raw.slot = "data",
  ...
) {
  object.to <- switch(
    EXPR = to,
    'seurat' = {
      raw.data.matrix <- sparseMatrix(
        i = as.numeric(x = from$raw$X$indices),
        p = as.numeric(x = from$raw$X$indptr),
        x = as.numeric(x = from$raw$X$data),
        index1 = FALSE
      )
      rownames(x = raw.data.matrix) <- rownames(x = py_to_r(from$raw$var))
      colnames(x = raw.data.matrix) <- rownames(x = py_to_r(from$obs))
      data.matrix <- t(x = py_to_r(from$X))
      rownames(x = data.matrix) <- rownames(x = py_to_r(from$var))
      colnames(x = data.matrix) <- rownames(x = py_to_r(from$obs))
      meta.data <- py_to_r(from$obs)
      if ("n_counts" %in% colnames(x = meta.data)) {
        colnames(x = meta.data) <- gsub(
          pattern = "n_counts",
          replacement = "nUMI",
          x = colnames(x = meta.data)
        )
      }
      if ("n_gene" %in% colnames(x = meta.data)) {
        colnames(x = meta.data) <- gsub(
          pattern = "n_gene",
          replacement = "nGene",
          x = colnames(x = meta.data)
        )
      }
      seurat.object <- CreateSeuratObject(
        raw.data = raw.data.matrix,
        meta.data = meta.data
      )
      seurat.object <- SetAssayData(
        object = seurat.object,
        assay.type = "RNA",
        slot = X.slot,
        new.data = data.matrix
      )
      #todo, deal with obsm fields that are not dimensional reductions, or have different name structures
      drs <- unlist(x = py_to_r(from$obsm$keys()))
      for (dr in drs) {
        dr.embed <- py_to_r(from$obsm[[eval(dr)]])
        dr.name <- ExtractField(string = dr, field = 2)
        if (is.na(dr.name)) {
          dr.name <- dr
        }
        dr.dict <- list(tSNE_ = "tsne", PC = "pca")
        if (dr.name %in% dr.dict) {
          dr.key <- names(x = which(x = dr.dict == dr.name))
        } else {
          dr.key <- toupper(x = dr.name)
        }
        colnames(x = dr.embed) <- paste0(dr.key, 1:ncol(x = dr.embed))
        rownames(x = dr.embed) <- seurat.object@cell.names
        seurat.object <- SetDimReduction(
          object = seurat.object,
          reduction.type = dr.name,
          slot = "cell.embeddings",
          new.data = dr.embed
        )
        seurat.object <- SetDimReduction(
          object = seurat.object,
          reduction.type = dr.name,
          slot = "key",
          new.data = dr.key
        )
      }
      seurat.object
    },
    stop(paste0("Cannot convert AnnData objects to class '", to, "'"))
  )
  return(object.to)
}

#' @param raw.data.slot name of the SingleCellExperiment assay to slot into @@raw.data
#' @param data.slot name of the SingleCellExperiment assay to slot into @@data
#'
#' @describeIn Convert Convert from SingleCellExperiment to a Seurat object
#' @export
#' @method Convert SingleCellExperiment
#'
Convert.SingleCellExperiment <- function(
  from,
  to,
  raw.data.slot = "counts",
  data.slot = "logcounts",
  ...
) {
  object.to <- switch(
    EXPR = to,
    'seurat' = {
      raw.data <- tryCatch(
        expr = SummarizedExperiment::assay(from, raw.data.slot),
        error = function(e) {
          stop(paste0("No data in provided assay - ", raw.data.slot))
        }
      )
      data <- tryCatch(
        expr = SummarizedExperiment::assay(from, data.slot),
        error = function(e) {
          stop(paste0("No data in provided assay - ", data.slot))
        }
      )
      meta.data <- as.data.frame(SummarizedExperiment::colData(from))
      seurat.object <- CreateSeuratObject(raw.data = raw.data, meta.data = meta.data)
      seurat.object@data <- data
      if (length(x = SingleCellExperiment::reducedDimNames(from)) > 0) {
        for (dr in SingleCellExperiment::reducedDimNames(from)) {
          seurat.object <- SetDimReduction(
            object = seurat.object,
            reduction.type = dr,
            slot = "cell.embeddings",
            new.data = SingleCellExperiment::reducedDim(x = from, type = dr)
          )
          key <- gsub(
            pattern = "[[:digit:]]",
            replacement = "",
            x = colnames(x = SingleCellExperiment::reducedDim(x = from, type = dr)
            )[1])
          seurat.object <- SetDimReduction(
            object = seurat.object,
            reduction.type = dr,
            slot = "key",
            new.data = key
          )
        }
      }
      seurat.object
    },
    stop(paste0("Cannot convert SingleCellExperiment objects to class '", to, "'"))
  )
  return(object.to)
}

#' @param filename Filename for writing files
#' @param chunk.dims Internal HDF5 chunk size
#' @param chunk.size Number of cells to stream to loom file at a time
#' @param overwrite Overwrite existing file at \code{filename}?
#' @param display.progress Display a progress bar
#' @param anndata.raw Name of matrix (raw.data, data) to put in the anndata raw slot
#' @param anndata.X Name of matrix (data, scale.data) to put in the anndata X slot
#'
#' @describeIn Convert Convert a Seurat object
#'
#' @importFrom utils installed.packages
#' @importFrom methods as slot
#' @importFrom reticulate import np_array tuple dict r_to_py
#'
#' @export
#' @method Convert seurat
#'
Convert.seurat <- function(
  from,
  to,
  filename,
  chunk.dims = 'auto',
  chunk.size = 1000,
  overwrite = FALSE,
  display.progress = TRUE,
  anndata.raw = "raw.data",
  anndata.X = "data",
  ...
) {
  object.to <- switch(
    EXPR = to,
    'loom' = {
      if (!'loomR' %in% rownames(x = installed.packages())) {
        stop("Please install loomR from GitHub before converting to a loom object")
      }
      cell.order <- from@cell.names
      gene.order <- rownames(x = from@raw.data)
      loomfile <- loomR::create(
        filename = filename,
        data = from@raw.data[, cell.order],
        cell.attrs = from@meta.data[cell.order, ],
        layers = list('norm_data' = t(x = from@data[, cell.order])),
        chunk.dims = chunk.dims,
        chunk.size = chunk.size,
        overwrite = overwrite,
        display.progress = display.progress
      )
      if (nrow(x = from@hvg.info) > 0) {
        hvg.info <- from@hvg.info
        colnames(x = hvg.info) <- gsub(
          pattern = '.',
          replacement = '_',
          x = colnames(x = hvg.info),
          fixed = TRUE
        )
        loomfile$add.row.attribute(hvg.info[gene.order, ])
      }
      if (length(x = from@var.genes) > 0) {
        loomfile$add.row.attribute(list('var_genes' = gene.order %in% from@var.genes))
      }
      if (!is.null(x = from@scale.data) && dim(x = from@scale.data) != c(1, 1)) {
        loomfile$add.layer(list(
          'scale_data' = as.matrix(x = t(x = as.data.frame(x = from@scale.data)[gene.order, cell.order]))
        ))
      }
      for (dim.reduc in names(x = from@dr)) {
        cell.embeddings <- from@dr[[dim.reduc]]@cell.embeddings
        ce.dims <- unique(x = dim(x = cell.embeddings))
        if (length(x = ce.dims) != 1 || ce.dims != 0) {
          if (nrow(x = cell.embeddings) < ncol(x = from@raw.data)) {
            cell.embeddings.padded <- matrix(
              nrow = ncol(x = from@raw.data),
              ncol = ncol(x = cell.embeddings)
            )
            if (is.null(x = rownames(x = cell.embeddings)) || is.null(x = from@cell.names)) {
              pad.order <- 1:nrow(x = cell.embeddings)
            } else {
              pad.order <- match(
                x = rownames(x = cell.embeddings),
                table = from@cell.names
              )
            }
            cell.embeddings.padded[pad.order, ] <- cell.embeddings
          } else if (nrow(x = cell.embeddings) > ncol(x = from@raw.data)) {
            stop("Cannot have more cells in the dimmensional reduction than in the dataset")
          } else {
            cell.embeddings.padded <- cell.embeddings
          }
          cell.embeddings.padded <- list(cell.embeddings.padded)
          names(x = cell.embeddings.padded) <- paste0(dim.reduc, '_cell_embeddings')
          loomfile$add.col.attribute(cell.embeddings.padded)
        }
        gene.loadings <- from@dr[[dim.reduc]]@gene.loadings
        gl.dims <- unique(x = dim(x = gene.loadings))
        if (length(x = gl.dims) == 1 && gl.dims == 0) {
          gene.loadings <- from@dr[[dim.reduc]]@gene.loadings.full
        }
        gl.dims <- unique(x = dim(x = gene.loadings))
        if (length(x = gl.dims) != 1 || gl.dims != 0) {
          if (nrow(x = gene.loadings) < nrow(x = from@raw.data)) {
            gene.loadings.padded <- matrix(
              nrow = nrow(x = from@raw.data),
              ncol = ncol(x = gene.loadings)
            )
            if (is.null(x = rownames(x = gene.loadings)) || is.null(x = rownames(x = from@raw.data))) {
              pad.order <- seq_len(nrow(x = gene.loadings))
            } else {
              pad.order <- match(
                x = rownames(x = gene.loadings),
                table = rownames(x = from@raw.data)
              )
            }
            gene.loadings.padded[pad.order, ] <- gene.loadings
          } else if (nrow(x = gene.loadings) > nrow(x = from@raw.data)) {
            stop("Cannot have more genes in the dimmensional reduction than in the dataset")
          } else {
            gene.loadings.padded <- gene.loadings
          }
          gene.loadings.padded <- list(gene.loadings.padded)
          names(x = gene.loadings.padded) <- paste0(dim.reduc, '_gene_loadings')
          loomfile$add.row.attribute(gene.loadings.padded)
        }
      }
      loomfile
    },
    'sce' = {
      if (!'SingleCellExperiment' %in% rownames(x = installed.packages())) {
        stop("Please install SingleCellExperiment from Bioconductor before converting to a SingeCellExperiment object")
      }
      if (inherits(x = from@raw.data, what = "data.frame")) {
        from@raw.data <- as.matrix(from@raw.data)
      }
      if (inherits(x = from@data, what = "data.frame")) {
        from@data <- as.matrix(from@data)
      }
      sce <- if (class(from@raw.data) %in% c("matrix", "dgTMatrix")) {
        SingleCellExperiment::SingleCellExperiment(assays = list(counts = as(from@raw.data[rownames(from@data), from@cell.names], "dgCMatrix")))
      } else if (inherits(x = from@raw.data, what = "dgCMatrix")) {
        SingleCellExperiment::SingleCellExperiment(assays = list(counts = from@raw.data[rownames(from@data), from@cell.names]))
      } else {
        stop("Invalid class stored in seurat object's raw.data slot")
      }
      if (class(from@data) %in% c("matrix", "dgTMatrix")) {
        SummarizedExperiment::assay(sce, "logcounts") <- as(from@data, "dgCMatrix")
      } else if (inherits(x = from@data, what = "dgCMatrix")) {
        SummarizedExperiment::assay(sce, "logcounts") <- from@data
      } else {
        stop("Invalid class stored in seurat object's data slot")
      }
      meta.data <- from@meta.data
      meta.data$ident <- from@ident
      SummarizedExperiment::colData(sce) <- S4Vectors::DataFrame(meta.data)
      row.data <- from@hvg.info[rownames(from@data), ]
      row.data <- cbind(gene = rownames(x = from@data), row.data, stringsAsFactors = FALSE)
      SummarizedExperiment::rowData(sce) <- S4Vectors::DataFrame(row.data)
      for (dr in names(from@dr)) {
        SingleCellExperiment::reducedDim(sce, toupper(x = dr)) <- slot(
          object = slot(object = from, name = "dr")[[dr]],
          name = "cell.embeddings"
        )
      }
      sce
    },
    'anndata' = {
      if (!py_module_available("anndata")) {
        stop("Please install the anndata python module")
      }
      ad <- import("anndata")
      raw <- switch(
        EXPR = anndata.raw,
        "raw.data" = from@raw.data,
        "data" = from@data,
        stop("Invalid Seurat data slot. Please choose one of: raw.data, data")
      )
      raw <- raw[,from@cell.names]
      X <- switch(
        EXPR = anndata.X,
        "data" = from@data,
        "scale.data" = from@scale.data,
        stop("Invalid Seurat data slot. Please choose one of: data, scale.data")
      )
      cell_names <- colnames(x = X)
      gene_names <- rownames(x = X)
      if (inherits(x = raw, what = c('matrix', 'Matrix'))) {
        raw <- as(object = raw, Class = "dgCMatrix")
      } else {
        raw <- as(object = as.matrix(x = raw), Class = "dgCMatrix")
      }
      scipy <- import(module = 'scipy.sparse', convert = FALSE)
      sp_sparse_csc <- scipy$csc_matrix
      raw.rownames <- rownames(x = raw)
      raw <- sp_sparse_csc(
        tuple(np_array(raw@x), np_array(raw@i), np_array(raw@p)),
        shape = tuple(raw@Dim[1], raw@Dim[2])
      )
      if (inherits(x = raw, what = c('matrix', 'Matrix', 'data.frame'))) {
        raw <- r_to_py(x = raw)
      }
      raw <- raw$T
      raw <- dict(X = raw, var = dict(var_names = raw.rownames))
      if (anndata.X == 'data') {
        X <- sp_sparse_csc(
          tuple(np_array(X@x), np_array(X@i), np_array(X@p)),
          shape = tuple(X@Dim[1], X@Dim[2])
        )
        X <- X$T
      } else {
        X <- np_array(t(x = X))
      }
      obsm <- list()
      for (dr in names(from@dr)) {
        obsm[[paste0("X_",dr)]] <- np_array(GetCellEmbeddings(
          object = from,
          reduction.type = dr
        ))
      }
      obsm <- if (!identical(obsm, list())) dict(obsm) else NULL
      meta_data <- from@meta.data
      if ("nUMI" %in% colnames(x = meta_data)) {
        colnames(x = meta_data) <- gsub(
          pattern = "nUMI",
          replacement = "n_counts",
          x = colnames(x = meta_data)
        )
      }
      if ("nGene" %in% colnames(x = meta_data)) {
        colnames(x = meta_data) <- gsub(
          pattern = "nGene",
          replacement = "n_genes",
          x = colnames(x = meta_data)
        )
      }
      colnames(x = meta_data) <- gsub(
        pattern = "\\.",
        replacement = "_",
        x = colnames(x = meta_data)
      )
      anndata.object <- ad$AnnData(
        raw = raw,
        X = X,
        obs = meta_data,
        var = from@hvg.info,
        obsm = obsm
      )
      anndata.object$var_names <- gene_names
      anndata.object$obs_names <- cell_names
      if (!missing(x = filename)) {
        anndata.object$write(filename)
      }
      anndata.object
    },
    stop(paste0("Cannot convert Seurat objects to class '", to, "'"))
  )
  return(object.to)
}

#' @describeIn DefaultAssay Get the name of the assay.used to calculate this DimReduc
#' @export
#' @method DefaultAssay DimReduc
#'
DefaultAssay.DimReduc <- function(object, ...) {
  return(slot(object = object, name = 'assay.used'))
}

#' @describeIn DefaultAssay Get the default assay of a Seurat object
#' @export
#' @method DefaultAssay Seurat
#'
DefaultAssay.Seurat <- function(object, ...) {
  return(slot(object = object, name = 'active.assay'))
}

#' @export
#' @method DefaultAssay<- DimReduc
#'
"DefaultAssay<-.DimReduc" <- function(object, ..., value) {
  slot(object = object, name = 'assay.used') <- value
  return(object)
}

#' @describeIn DefaultAssay Set the default assay of a Seurat object
#' @export
#' @method DefaultAssay<- Seurat
#'
"DefaultAssay<-.Seurat" <- function(object, ..., value) {
  if (!value %in% names(x = slot(object = object, name = 'assays'))) {
    stop("Cannot find assay ", value)
  }
  slot(object = object, name = 'active.assay') <- value
  return(object)
}

#' @describeIn Embeddings Get the cell embeddings from a DimReduc object
#' @export
#' @method Embeddings DimReduc
#'
Embeddings.DimReduc <- function(object, ...) {
  return(slot(object = object, name = 'cell.embeddings'))
}

#' @describeIn GetAssay Get an assay from a Seurat object
#' @export
#' @method GetAssay Seurat
#'
GetAssay.Seurat <- function(object, assay = NULL) {
  assay <- assay %||% DefaultAssay(object = object)
  object.assays <- FilterObjects(object = object, classes.keep = 'Assay')
  if (!assay %in% object.assays) {
    stop(paste0(
      assay,
      " is not an assay present in the given object. Available assays are: ",
      paste(object.assays, collapse = ", ")
    ))
  }
  return(slot(object = object, name = 'assays')[[assay]])
}

#' @describeIn GetAssayData Get assay data for an Assay object
#' @export
#' @method GetAssayData Assay
#'
GetAssayData.Assay <- function(object, slot = 'data') {
  return(slot(object = object, name = slot))
}

#' @param assay Name of assay to pull data from
#'
#' @describeIn GetAssayData Get assay data from a Seurat object
#' @export
#' @method GetAssayData Seurat
#'
GetAssayData.Seurat <- function(object, slot = 'data', assay = NULL, ...) {
  assay <- assay %||% DefaultAssay(object = object)
  return(GetAssayData(
    object = GetAssay(object = object, assay = assay),
    slot = slot
  ))
}

#' @describeIn HVFInfo Get highly variable feature information from an Assay
#' @export
#' @method HVFInfo Assay
#'
HVFInfo.Assay <- function(object, ...) {
  vars <- c(
    'mean',
    if ('variance.standardized' %in% colnames(x = object[[]])) {
      c('variance', 'variance.standardized')
    } else {
      c('dispersion', 'dispersion.scaled')
    }
  )
  hvf.info <- object[[vars]]
  return(hvf.info)
}

#' @param assay Name of assay to pull highly variable feature information for
#'
#' @describeIn HVFInfo Get highly variable feature information from a Seurat object
#' @export
#' @method HVFInfo Seurat
#'
HVFInfo.Seurat <- function(object, assay = NULL, ...) {
  assay <- assay %||% DefaultAssay(object = object)
  return(HVFInfo(object = GetAssay(object = object, assay = assay)))
}

#' @rdname Idents
#' @export
#' @method Idents Seurat
#'
Idents.Seurat <- function(object, ...) {
  return(slot(object = object, name = 'active.ident'))
}

#' @param cells Set cell identities for specific cells
#'
#' @rdname Idents
#' @export
#' @method Idents<- Seurat
#'
"Idents<-.Seurat" <- function(object, cells = NULL, ..., value) {
  cells <- cells %||% colnames(x = object)
  if (is.numeric(x = cells)) {
    cells <- colnames(x = object)[cells]
  }
  cells <- intersect(x = cells, y = colnames(x = object))
  cells <- match(x = cells, table = colnames(x = object))
  if (length(x = cells) == 0) {
    stop("Cannot find cells provided")
  }
  idents.new <- if (length(x = value) == 1 && value %in% colnames(x = object[[]])) {
    unlist(x = object[[value]], use.names = FALSE)[cells]
  } else {
    if (is.list(x = value)) {
      value <- unlist(x = value, use.names = FALSE)
    }
    rep_len(x = value, length.out = length(x = cells))
  }
  idents.new <- as.vector(x = idents.new)
  idents <- as.vector(x = Idents(object = object))
  idents[cells] <- idents.new
  idents[is.na(x = idents)] <- 'NA'
  names(x = idents) <- colnames(x = object)
  missing.cells <- which(x = is.na(x = names(x = idents)))
  if (length(x = missing.cells) > 0) {
    idents <- idents[-missing.cells]
  }
  idents <- factor(x = idents)
  slot(object = object, name = 'active.ident') <- idents
  return(object)
}

#' @describeIn JS Get JackStraw information from a DimReduc object
#' @export
#' @method JS DimReduc
#'
JS.DimReduc <- function(object, slot = NULL, ...) {
  jackstraw <- slot(object = object, name = 'jackstraw')
  if (!is.null(x = slot)) {
    jackstraw <- JS(object = jackstraw, slot = slot)
  }
  return(jackstraw)
}

#' @describeIn JS Get Jackstraw information from a JackStrawData object
#' @export
#' @method JS JackStrawData
#'
JS.JackStrawData <- function(object, slot, ...) {
  slot <- switch(
    EXPR = slot,
    'empirical' = 'empirical.p.values',
    'fake' = 'fake.reduction.scores',
    'full' = 'empirical.p.values.full',
    'overall' = 'overall.p.values',
    slot
  )
  return(slot(object = object, name = slot))
}

#' @describeIn JS Set JackStraw information for a DimReduc object
#' @export
#' @method JS<- DimReduc
#'
"JS<-.DimReduc" <- function(object, slot = NULL, ..., value) {
  if (inherits(x = value, what = 'JackStrawData')) {
    slot(object = object, name = 'jackstraw') <- value
  } else if (is.null(x = NULL)) {
    stop("A slot must be specified")
  } else {
    JS(object = JS(object = object), slot = slot) <- value
  }
  return(object)
}

#' @describeIn JS Get Jackstraw information from a JackStrawData object
#' @export
#' @method JS<- JackStrawData
#'
"JS<-.JackStrawData" <- function(object, slot, ..., value) {
  slot <- switch(
    EXPR = slot,
    'empirical' = 'empirical.p.values',
    'fake' = 'fake.reduction.scores',
    'full' = 'empirical.p.values.full',
    'overall' = 'overall.p.values',
    slot
  )
  slot(object = object, name = slot) <- value
  return(object)
}

#' @describeIn Key Get the key for an Assay object
#' @export
#' @method Key Assay
#'
Key.Assay <- function(object, ...) {
  return(slot(object = object, name = 'key'))
}

#' @describeIn Key Get the key for a DimReduc object
#' @export
#' @method Key DimReduc
#'
Key.DimReduc <- function(object, ...) {
  return(slot(object = object, name = 'key'))
}

#' @describeIn Key Set the key for an Assay object
#' @export
#' @method Key<- Assay
#'
"Key<-.Assay" <- function(object, ..., value) {
  slot(object = object, name = 'key') <- value
  return(object)
}

#' @describeIn Loadings Get the feature loadings from a DimReduc object
#' @export
#' @method Loadings DimReduc
#'
Loadings.DimReduc <- function(object, projected = NULL, ...) {
  projected <- projected %||% Projected(object = object)
  slot <- ifelse(
    test = projected,
    yes = 'feature.loadings.projected',
    no = 'feature.loadings'
  )
  return(slot(object = object, name = slot))
}

#' @describeIn Loadings Add projected feature loadings to a DimReduc object
#' @export
#' @method Loadings<- DimReduc
#'
"Loadings<-.DimReduc" <- function(object, projected = TRUE, ..., value) {
  slot.use <- ifelse(
    test = projected,
    yes = 'feature.loadings.projected',
    no = 'feature.loadings'
  )
  if (ncol(x = value) != length(x = object)) {
    stop("New feature loadings must have the dimensions as currently calculated")
  }
  slot(object = object, name = slot.use) <- value
  return(object)
}

#' @describeIn Misc Get miscellaneous data from a Seurat object
#' @export
#' @method Misc Seurat
#'
Misc.Seurat <- function(object, slot = NULL, ...) {
  if (is.null(x = slot)) {
    return(slot(object = object, name = 'misc'))
  }
  return(slot(object = object, name = 'misc')[[slot]])
}

#' @describeIn Misc Set miscellaneous data for a Seurat object
#' @export
#' @method Misc<- Seurat
#'
"Misc<-.Seurat" <- function(object, slot, ..., value) {
  if (slot %in% names(x = Misc(object = object))) {
    warning("Overwriting miscellanous data for ", slot)
  }
  slot(object = object, name = 'misc')[[slot]] <- value
  return(object)
}

#' @describeIn WhichCells Get cells from an Assay object
#' @export
#' @method OldWhichCells Assay
#'
OldWhichCells.Assay <- function(
  object,
  cells,
  subset.name = NULL,
  low.threshold = -Inf,
  high.threshold = Inf,
  accept.value = NULL,
  ...
) {
  cells <- cells %||% colnames(x = object)
  # input checking
  if (length(x = subset.name) > 1) {
    stop("subset.name must be a single parameter")
  }
  if (length(x = low.threshold) > 1 | length(x = high.threshold) > 1) {
    stop("Multiple values passed to low.threshold or high.threshold")
  }
  if (low.threshold >= high.threshold) {
    stop("low.threshold is greater than or equal to high.threshold")
  }
  if (!is.null(x = subset.name)) {
    subset.name <- as.character(x = subset.name)
    data.use <- GetAssayData(
      object = object,
      ... = ...
    )
    data.use <- t(x = data.use[subset.name, cells, drop = FALSE])
    if (!is.null(x = accept.value)) {
      if (!all(accept.value %in% unique(x = data.use[, 1]))) {
        bad.vals <- accept.value[!(accept.value %in% unique(x = data.use[, 1]))]
        stop("Identity: ", bad.vals, " not found.")
      }
      pass.inds <- which(x = apply(data.use, MARGIN = 1, function(x) x %in% accept.value))
    } else {
      pass.inds <- which(x = (data.use > low.threshold) & (data.use < high.threshold))
    }
    cells <- rownames(x = data.use)[pass.inds]
  }
  return(cells)
}

#' @param ident.keep Create a cell subset based on the provided identity classes
#' @param ident.remove Subtract out cells from these identity classes (used for
#' filtration)
#' @param max.cells.per.ident Can be used to downsample the data to a certain
#' max per cell ident. Default is INF.
#' @param random.seed Random seed for downsampling
#' @param assay Which assay to filter on
#' @param ... Extra parameters passed to \code{FetchData}
#'
#' @seealso \code{\link{FetchData}}
#'
#' @describeIn WhichCells Get cells from a Seurat object
#' @export
#' @method OldWhichCells Seurat
#'
OldWhichCells.Seurat <- function(
  object,
  cells = NULL,
  subset.name = NULL,
  low.threshold = -Inf,
  high.threshold = Inf,
  accept.value = NULL,
  ident.keep = NULL,
  ident.remove = NULL,
  max.cells.per.ident = Inf,
  random.seed = 1,
  assay = NULL,
  ...
) {
  # input checking
  if (length(x = subset.name) > 1) {
    stop("subset.name must be a single parameter")
  }
  if (length(x = low.threshold) > 1 | length(x = high.threshold) > 1) {
    stop("Multiple values passed to low.threshold or high.threshold")
  }
  if (low.threshold >= high.threshold) {
    stop("low.threshold is greater than or equal to high.threshold")
  }
  if (!is.na(x = random.seed)) {
    set.seed(seed = random.seed)
  }
  cells <- cells %||% colnames(x = object)
  assay <- assay %||% DefaultAssay(object = object)
  ident.keep <- ident.keep %||% unique(x = Idents(object = object))
  bad.remove.idents <- ident.remove[!ident.remove %in% unique(x = Idents(object = object))]
  if (length(bad.remove.idents) > 0) {
    stop(paste("Identity :", bad.remove.idents, "not found.   "))
  }
  ident.keep <- setdiff(x = ident.keep, y = ident.remove)
  if (!all(ident.keep %in% unique(Idents(object = object)))) {
    bad.idents <- ident.keep[!(ident.keep %in% unique(x = Idents(object = object)))]
    stop("Identity: ", bad.idents, " not found.")
  }
  cells.to.use <- character()
  for (id in ident.keep) {
    cells.in.ident <- Idents(object = object)[cells]
    cells.in.ident <- names(x = cells.in.ident[cells.in.ident == id])
    cells.in.ident <- cells.in.ident[!is.na(x = cells.in.ident)]
    if (length(x = cells.in.ident) > max.cells.per.ident) {
      cells.in.ident <- sample(x = cells.in.ident, size = max.cells.per.ident)
    }
    cells.to.use <- c(cells.to.use, cells.in.ident)
  }
  cells <- cells.to.use
  if (!is.null(x = subset.name)) {
    subset.name <- as.character(subset.name)
    data.use <- FetchData(
      object = object,
      vars = subset.name,
      cells = cells,
      ...
    )
    if (!is.null(x = accept.value)) {
      if (!all(accept.value %in% unique(x = data.use[, 1]))) {
        bad.vals <- accept.value[!accept.value %in% unique(x = data.use[, 1])]
        stop("Identity: ", bad.vals, " not found.")
      }
      pass.inds <- which(x = apply(X = data.use, MARGIN = 1, FUN = function(x) x %in% accept.value))
    } else {
      pass.inds <- which(x = (data.use > low.threshold) & (data.use < high.threshold))
    }
    cells <- rownames(x = data.use)[pass.inds]
  }
  return(cells)
}

#' @param dims Number of dimensions to display
#' @param nfeatures Number of genes to display
#' @param projected Use projected slot
#'
#' @export
#'
#' @describeIn Print Print top features for DimReduc
#' @method Print DimReduc
#'
Print.DimReduc <- function(
  object,
  dims = 1:5,
  nfeatures = 20,
  projected = FALSE
) {
  loadings <- Loadings(object = object, projected = projected)
  nfeatures <- min(nfeatures, nrow(x = loadings))
  if (ncol(x = loadings) == 0) {
    warning("Dimensions have not been projected. Setting projected = FALSE")
    projected <- FALSE
    loadings <- Loadings(object, projected = projected)
  }
  if (max(dims) > ncol(x = loadings)) {
    stop(paste0("Only ", ncol(x = loadings), " dimensions have been computed."))
  }
  for (dim in dims) {
    features <- TopFeatures(
      object = object,
      dim = dim,
      nfeatures = nfeatures * 2,
      projected = projected,
      balanced = TRUE
    )
   message(paste0(Key(object = object), dim))
   pos.features <- split(x = features$positive, f = ceiling(x = seq_along(along.with = features$positive) / 10))
   message(paste0("Positive: "), paste(pos.features[[1]], collapse = ", "))
   pos.features[[1]] <- NULL
   if (length(x = pos.features) > 0) {
     for (i in pos.features) {
       message(paste0("\t  ", paste(i, collapse = ", ")))
     }
   }
   neg.features <- split(x = features$negative, f = ceiling(x = seq_along(along.with = features$negative) / 10))
   message(paste0("Negative: "), paste(neg.features[[1]], collapse = ", "))
   neg.features[[1]] <- NULL
   if (length(x = neg.features) > 0) {
     for (i in neg.features) {
       message(paste0("\t  ", paste(i, collapse = ", ")))
     }
   }
   message("")
  }
}

#' @describeIn RenameCells Rename cells in an Assay object
#' @export
#' @method RenameCells Assay
#'
RenameCells.Assay <- function(
  object,
  new.names = NULL
) {
  for (data.slot in c("counts", "data", "scale.data")) {
    old.data <- GetAssayData(object = object, slot = data.slot)
    if (ncol(x = old.data) == 0) {
      next
    }
    colnames(x = slot(object = object, name = data.slot)) <- new.names
  }
  return(object)
}

#' @describeIn RenameCells Rename cells in a DimReduc object
#' @export
#' @method RenameCells DimReduc
#'
RenameCells.DimReduc <- function(
  object,
  new.names = NULL
) {
  old.data <- Embeddings(object = object)
  rownames(old.data) <- new.names
  slot(object = object, name = "cell.embeddings") <- old.data
  return(object)
}

#' @param for.merge Only rename slots needed for merging Seurat objects.
#' Currently only renames the raw.data and meta.data slots.
#' @param add.cell.id prefix to add cell names
#'
#' @describeIn RenameCells Rename cells in a Seurat object
#' @export
#' @method RenameCells Seurat
#'
RenameCells.Seurat <- function(
  object,
  add.cell.id = NULL,
  new.names = NULL,
  for.merge = FALSE
) {
  if (missing(x = add.cell.id) && missing(x = new.names)) {
    stop("One of 'add.cell.id' and 'new.names' must be set")
  }
  if (!missing(x = add.cell.id) && !missing(x = new.names)) {
    stop("Only one of 'add.cell.id' and 'new.names' may be set")
  }
  if (!missing(x = add.cell.id)) {
    new.cell.names <- paste(add.cell.id, colnames(x = object), sep = "_")
  } else {
    if (length(x = new.names) == ncol(x = object)) {
      new.cell.names <- new.names
    } else {
      stop(
        "the length of 'new.names' (",
        length(x = new.names),
        ") must be the same as the number of cells (",
        ncol(x = object),
        ")"
      )
    }
  }
  # rename in the assay objects
  assays <- FilterObjects(object = object, classes.keep = 'Assay')
  for (assay in assays) {
    slot(object = object, name = "assays")[[assay]] <- RenameCells(
      object = object[[assay]],
      new.names = new.cell.names
    )
  }
  # rename in the DimReduc objects
  dimreducs <- FilterObjects(object = object, classes.keep = 'DimReduc')
  for (dr in dimreducs) {
    object[[dr]] <- RenameCells(
      object = object[[dr]],
      new.names = new.cell.names
    )
  }
  # rename the active.idents
  old.ids <- Idents(object = object)
  names(x = old.ids) <- new.cell.names
  Idents(object = object, cells = rownames(x = object[[]])) <- old.ids
  # rename the cell-level metadata
  old.meta.data <- object[[]]
  rownames(x = old.meta.data) <- new.cell.names
  slot(object = object, name = "meta.data") <- old.meta.data
  # rename the graphs
  graphs <- FilterObjects(object = object, classes.keep = "Graph")
  for (g in graphs) {
    colnames(x = object[[g]]) <- new.cell.names
    rownames(x = object[[g]]) <- new.cell.names
  }
  return(object)
}

#' @rdname Idents
#' @export
#' @method RenameIdents Seurat
#'
RenameIdents.Seurat <- function(object, ...) {
  ident.pairs <- list(...)
  if (is.null(x = names(x = ident.pairs))) {
    stop("All arguments must be named with the old identity class")
  }
  if (!all(sapply(X = ident.pairs, FUN = length) == 1)) {
    stop("Can only rename identity classes to one value")
  }
  if (!any(names(x = ident.pairs) %in% levels(x = object))) {
    stop("Cannot find any of the provided identities")
  }
  cells.idents <- CellsByIdentities(object = object)
  for (i in names(x = ident.pairs)) {
    if (!i %in% names(x = cells.idents)) {
      warning("Cannot find identity ", i, call. = FALSE, immediate. = TRUE)
      next
    }
    Idents(object = object, cells = cells.idents[[i]]) <- ident.pairs[[i]]
  }
  return(object)
}

#' @describeIn SetAssayData Set assay data for an Assay object
#' @export
#' @method SetAssayData Assay
#'
SetAssayData.Assay <- function(object, slot, new.data) {
  slots.use <- c('counts', 'data', 'scale.data')
  if (!slot %in% slots.use) {
    stop("'slot' must be one of ", paste(slots.use, collapse = ', '))
  }
  if (!IsMatrixEmpty(x = new.data) && !slot %in% c('counts', 'scale.data')) {
    if (ncol(x = new.data) != ncol(x = object)) {
      stop("The new data doesn't have the same number of cells as the current data")
    }
    num.counts <- nrow(x = GetAssayData(object = object, slot = 'counts'))
    counts.names <- rownames(x = GetAssayData(object = object, slot = 'counts'))
    if (num.counts <= 1) {
      num.counts <- nrow(x = object)
      counts.names <- rownames(x = object)
    }
    if (slot == 'counts' && nrow(x = new.data) != num.counts) {
      warning("The new data doesn't have the same number of features as the current data")
    } else if (slot %in% c('data', 'scale.data') && nrow(x = new.data) > num.counts) {
      warning("Adding more features than present in current data")
    }
    if (!all(rownames(x = new.data) %in% counts.names)) {
      warning("Adding features not currently present in the object")
    }
    new.features <- na.omit(object = match(
      x = counts.names,
      table = rownames(x = new.data)
    ))
    #if (slot == 'scale.data' && nrow(x = new.data) > nrow(x = object)) {
    #  stop("Cannot add more features than present in current data")
    #} else if (slot != 'scale.data' && nrow(x = new.data) != nrow(x = object)) {
    #  stop("The new data doesn't have the same number of features as the current data")
    #}
    new.cells <- colnames(x = new.data)
    if (!all(new.cells %in% colnames(x = object))) {
      stop("All cell names must match current cell names")
    }
    new.data <- new.data[new.features, colnames(x = object)]
  }
  slot(object = object, name = slot) <- new.data
  return(object)
}

#' @param assay Name of assay whose data should be set
#'
#' @describeIn SetAssayData Set assay data for an Assay object in a Seurat object
#' @export
#' @method SetAssayData Seurat
#'
SetAssayData.Seurat <- function(
  object,
  slot = 'data',
  new.data,
  assay = NULL,
  ...
) {
  assay <- assay %||% DefaultAssay(object = object)
  object[[assay]] <- SetAssayData(object = object[[assay]], slot = slot, new.data = new.data)
  return(object)
}

#' @inheritParams Idents
#'
#' @rdname Idents
#' @export
#' @method SetIdent Seurat
#'
SetIdent.Seurat <- function(object, cells = NULL, value) {
  Idents(object = object, cells = cells) <- value
  return(object)
}

#' @rdname Idents
#' @export
#' @method StashIdent Seurat
#'
StashIdent.Seurat <- function(object, save.name = 'orig.ident', ...) {
  object[[save.name]] <- Idents(object = object)
  return(object)
}

#' @describeIn Stdev Get the standard deviations from a DimReduc object
#' @export
#' @method Stdev DimReduc
#'
Stdev.DimReduc <- function(object) {
  return(slot(object = object, name = 'stdev'))
}

#' @param reduction Name of reduction to use
#'
#' @describeIn Stdev Get the standard deviations of a dimensional reduction from a Seurat object
#' @export
#' @method Stdev Seurat
#'
Stdev.Seurat <- function(object, reduction, ...) {
  return(Stdev(object = object[[reduction]]))
}

#' @describeIn SubsetData Subset an Assay object
#' @export
#' @method SubsetData Assay
#'
SubsetData.Assay <- function(
  object,
  cells = NULL,
  subset.name = NULL,
  low.threshold = -Inf,
  high.threshold = Inf,
  accept.value = NULL,
  do.clean = FALSE,
  ...
) {
  cells <- cells %||% colnames(x = object)
  cells <- OldWhichCells(
    object = object,
    cells = cells,
    subset.name = subset.name,
    low.threshold = low.threshold,
    high.threshold = high.threshold,
    accept.value = accept.value,
    ...
  )
  if (ncol(x = GetAssayData(object = object, slot = 'counts')) == ncol(x = object)) {
    slot(object = object, name = "counts") <- GetAssayData(object = object, slot = "counts")[, cells]
  }
  slot(object = object, name = "data") <- GetAssayData(object = object, slot = "data")[, cells]
  cells.scaled <- colnames(x = GetAssayData(object = object, slot = "scale.data"))
  cells.scaled <- cells.scaled[cells.scaled %in% cells]
  if (length(x = cells.scaled) > 0) {
    slot(object = object, name = "scale.data") <- GetAssayData(object = object, slot = "scale.data")[, cells]
  }
  return(object)
}

#' @param assay Assay to subset on
#' @param ident.use Create a cell subset based on the provided identity classes
#' @param ident.remove Subtract out cells from these identity classes (used for
#' filtration)
#' @param max.cells.per.ident Can be used to downsample the data to a certain
#' max per cell ident. Default is INF.
#' @param random.seed Random seed for downsampling
#'
#' @describeIn SubsetData Subset an Seurat object
#' @export
#' @method SubsetData Seurat
#'
SubsetData.Seurat <- function(
  object,
  assay = NULL,
  cells = NULL,
  subset.name = NULL,
  ident.use = NULL,
  ident.remove = NULL,
  low.threshold = -Inf,
  high.threshold = Inf,
  accept.value = NULL,
  max.cells.per.ident = Inf,
  random.seed = 1,
  do.clean = FALSE,
  ...
) {
  assay <- assay %||% DefaultAssay(object = object)
  cells <- OldWhichCells(
    object = object,
    assay = assay,
    ident = ident.use,
    ident.remove = ident.remove,
    subset.name = subset.name,
    cells = cells,
    max.cells.per.ident = max.cells.per.ident,
    random.seed = random.seed,
    low.threshold = low.threshold,
    high.threshold = high.threshold,
    accept.value = accept.value,
    ...
  )
  # Subset all the Assays
  assays <- FilterObjects(object = object, classes.keep = 'Assay')
  for (assay in assays) {
    slot(object = object, name = "assays")[[assay]] <- SubsetData(
      object = object[[assay]],
      cells = cells
    )
  }
  # Subset all the DimReducs
  drs <- FilterObjects(object = object, classes.keep = 'DimReduc')
  for (dr in drs) {
    object[[dr]] <- CreateDimReducObject(
      embeddings = Embeddings(object = object[[dr]])[cells, ],
      loadings = Loadings(object = object[[dr]], projected = FALSE),
      projected = Loadings(object = object[[dr]], projected = TRUE),
      assay = DefaultAssay(object = object[[dr]]),
      stdev = Stdev(object = object[[dr]]),
      key = Key(object = object[[dr]]),
      jackstraw = slot(object = object[[dr]], name = "jackstraw"),
      misc = slot(object[[dr]], name = "misc")
    )
  }
  slot(object = object, name = "active.ident") <- Idents(object = object)[cells]
  slot(object = object, name = "meta.data") <- slot(object = object, name = "meta.data")[cells, ]
  return(object)
}

#' @describeIn VariableFeatures Get the variable features of an assay object
#' @export
#' @method VariableFeatures Assay
#'
VariableFeatures.Assay <- function(object, ...) {
  return(slot(object = object, name = 'var.features'))
}

#' @param assay Name of assay to pull variable features for
#'
#' @describeIn VariableFeatures Get the variable features of a Seurat object
#' @export
#' @method VariableFeatures Seurat
#'
VariableFeatures.Seurat <- function(object, assay = NULL, ...) {
  assay <- assay %||% DefaultAssay(object = object)
  return(VariableFeatures(object = object[[assay]]))
}

#' @describeIn VariableFeatures Set the variable features of an assay object
#' @export
#' @method VariableFeatures<- Assay
#'
"VariableFeatures<-.Assay" <- function(object, ..., value) {
  slot(object = object, name = 'var.features') <- value
  return(object)
}

#' @inheritParams VariableFeatures.Seurat
#'
#' @describeIn VariableFeatures Set variable features for a Seurat object
#' @export
#' @method VariableFeatures<- Seurat
#'
"VariableFeatures<-.Seurat" <- function(object, assay = NULL, ..., value) {
  assay <- assay %||% DefaultAssay(object = object)
  VariableFeatures(object = object[[assay]]) <- value
  return(object)
}

#' @describeIn WhichCells Get cells from an Assay object
#' @export
#' @method WhichCells Assay
#'
WhichCells.Assay <- function(
  object,
  cells = NULL,
  expression,
  invert = FALSE,
  ...
) {
  cells <- cells %||% colnames(x = object)
  if (!missing(x = expression)) {
    key.pattern <- paste0('^', Key(object = object))
    expr <- substitute(expr = expression)
    expr.char <- as.character(x = expr)
    expr.char <- unlist(x = lapply(X = expr.char, FUN = strsplit, split = ' '))
    expr.char <- gsub(
      pattern = key.pattern,
      replacement = '',
      x = expr.char,
      perl = TRUE
    )
    expr.char <- gsub(
      pattern = '(',
      replacement = '',
      x = expr.char,
      fixed = TRUE
    )
    vars.use <- which(x = expr.char %in% rownames(x = object))
    expr.char <- expr.char[vars.use]
    data.subset <- as.data.frame(x = t(x = as.matrix(x = object[expr.char, ])))
    colnames(x = data.subset) <- expr.char
    data.subset <- subset.data.frame(x = data.subset, subset = eval(expr = expr))
    cells <- rownames(x = data.subset)
  }
  if (invert) {
    cells <- colnames(x = object)[!colnames(x = object) %in% cells]
  }
  return(cells)
}

#' @param idents A vector of identity classes to keep
#' @param downsample Maximum number of cells per identity class, default is \code{Inf};
#' downsampling will happen after all other operations, including inverting the
#' cell selection
#' @param seed Random seed for downsampling
#'
#' @describeIn WhichCells Get cells from a Seurat object
#' @export
#' @method WhichCells Seurat
#'
WhichCells.Seurat <- function(
  object,
  cells = NULL,
  idents = NULL,
  expression,
  invert = FALSE,
  downsample = Inf,
  seed = 1,
  ...
) {
  cells <- cells %||% colnames(x = object)
  if (is.numeric(x = cells)) {
    cells <- colnames(x = object)[cells]
  }
  if (!is.null(x = idents)) {
    set.seed(seed = seed)
    if (any(!idents %in% levels(x = Idents(object = object)))) {
      stop(
        "Cannot find the following identities in the object: ",
        paste(
          idents[!idents %in% levels(x = Idents(object = object))],
          sep = ', '
        )
      )
    }
    cells.idents <- unlist(x = lapply(
      X = idents,
      FUN = function(i) {
        cells.use <- which(x = as.vector(x = Idents(object = object)) == i)
        cells.use <- names(x = Idents(object = object)[cells.use])
        return(cells.use)
      }
    ))
    cells <- intersect(x = cells, y = cells.idents)
  }
  if (!missing(x = expression)) {
    objects.use <- FilterObjects(object = object)
    object.keys <- sapply(
      X = objects.use,
      FUN = function(i) {
        return(Key(object = object[[i]]))
      }
    )
    key.pattern <- paste0('^', object.keys, collapse = '|')
    expr <- if (is.call(x = substitute(expr = expression))) {
      substitute(expr = expression)
    } else {
      parse(text = expression)
    }
    expr.char <- as.character(x = expr)
    expr.char <- unlist(x = lapply(X = expr.char, FUN = strsplit, split = ' '))
    expr.char <- gsub(
      pattern = '(',
      replacement = '',
      x = expr.char,
      fixed = TRUE
    )
    vars.use <- which(
      x = expr.char %in% rownames(x = object) |
        expr.char %in% colnames(x = object[[]]) |
        grepl(pattern = key.pattern, x = expr.char, perl = TRUE)
    )
    data.subset <- FetchData(
      object = object,
      vars = expr.char[vars.use],
      cells = cells
    )
    data.subset <- subset.data.frame(x = data.subset, subset = eval(expr = expr))
    cells <- rownames(x = data.subset)
  }
  if (invert) {
    cells <- colnames(x = object)[!colnames(x = object) %in% cells]
  }
  cells <- CellsByIdentities(object = object, cells = cells)
  cells <- lapply(
    X = cells,
    FUN = function(x) {
      if (length(x = x) > downsample) {
        x <- sample(x = x, size = downsample, replace = FALSE)
      }
      return(x)
    }
  )
  return(unlist(x = cells, use.names = FALSE))
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Methods for R-defined generics
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @export
#'
".DollarNames.JackStrawData" <- function(x, pattern = '') {
  utils:::findMatches(pattern = pattern, values = slotNames(x = x))
}

#' @importFrom utils .DollarNames
#' @export
#'
".DollarNames.Seurat" <- function(x, pattern = '') {
  utils:::findMatches(pattern = pattern, values = colnames(x = x[[]]))
}

#' @importFrom utils .DollarNames
#' @export
#'
".DollarNames.SeuratCommand" <- function(x, pattern = '') {
  params <- slot(object = x, name = "params")
  utils:::findMatches(pattern = pattern, names(x = params))
}

#' @export
#'
"$.JackStrawData" <- function(x, i, ...) {
  return(slot(object = x, name = i))
}

#' @export
#'
"$.Seurat" <- function(x, i, ...) {
  return(x[[i, drop = TRUE]])
}

#' @export
#'
"$.SeuratCommand" <- function(x, i, ...) {
  params <- slot(object = x, name = "params")
  return(params[[i]])
}

#' @export
#'
"$<-.Seurat" <- function(x, i, ..., value) {
  x[[i]] <- value
  return(x)
}

#' @export
#'
"[.Assay" <- function(x, i, j, ...) {
  if (missing(x = i)) {
    i <- 1:nrow(x = x)
  }
  if (missing(x = j)) {
    j <- 1:ncol(x = x)
  }
  return(GetAssayData(object = x)[i, j, ..., drop = FALSE])
}

#' @export
#'
"[.DimReduc" <- function(x, i, j, ...) {
  if (missing(x = i)) {
    i <- 1:nrow(x = x)
  }
  if (missing(x = j)) {
    j <- 1:length(x = x)
  }
  return(Loadings(object = x)[i, j, ...])
}

#' @inheritParams subset.Seurat
#' @param i A vector of features to keep
#' @param j A vector of cells to keep
#'
#' @rdname subset.Seurat
#' @export
#'
#' @examples
#' pbmc_small[VariableFeatures(object = pbmc_small), ]
#' pbmc_small[, 1:10]
#'
"[.Seurat" <- function(x, i, j, ...) {
  if (missing(x = i) && missing(x = j)) {
    return(x)
  }
  if (missing(x = i)) {
    i <- NULL
  } else if (missing(x = j)) {
    j <- colnames(x = x)
  }
  if (is.numeric(x = i)) {
    i <- rownames(x = x)[i]
  }
  if (is.numeric(x = j)) {
    j <- colnames(x = x)[j]
  }
  return(subset.Seurat(x = x, select = c(i, j), ...))
}

#' @export
#'
"[.SeuratCommand" <- function(x, i, ...) {
  slot.use <- c("name", "timestamp", "call_string", "params")
  if (!i %in% slot.use) {
    stop("Invalid slot")
  }
  return(slot(object = x, name = i))
}

#' @export
#'
"[[.Assay" <- function(x, i, ..., drop = FALSE) {
  if (missing(x = i)) {
    i <- colnames(x = slot(object = x, name = 'meta.features'))
  }
  data.return <- slot(object = x, name = 'meta.features')[, i, drop = FALSE, ...]
  if (drop) {
    data.return <- unlist(x = data.return, use.names = FALSE)
    names(x = data.return) <- rep.int(x = rownames(x = x), times = length(x = i))
  }
  return(data.return)
}

#' @export
#'
"[[.DimReduc" <- function(x, i, j, ...) {
  if (missing(x = i)) {
    i <- 1:ncol(x = x)
  }
  if (missing(x = j)) {
    j <- 1:length(x = x)
  }
  return(Embeddings(object = x)[i, j, ...])
}

#' @export
#'
"[[.Seurat" <- function(x, i, ..., drop = FALSE) {
  if (missing(x = i)) {
    i <- colnames(x = slot(object = x, name = 'meta.data'))
  }
  if (length(x = i) == 0) {
    return(data.frame(row.names = colnames(x = x)))
  } else if (length(x = i) > 1 || any(i %in% colnames(x = slot(object = x, name = 'meta.data')))) {
    if (any(!i %in% colnames(x = slot(object = x, name = 'meta.data')))) {
      warning(
        "Cannot find the following bits of meta data: ",
        paste0(
          i[!i %in% colnames(x = slot(object = x, name = 'meta.data'))],
          collapse = ', '
        )
      )
    }
    i <- i[i %in% colnames(x = slot(object = x, name = 'meta.data'))]
    data.return <- slot(object = x, name = 'meta.data')[, i, drop = FALSE, ...]
    if (drop) {
      data.return <- unlist(x = data.return, use.names = FALSE)
      names(x = data.return) <- rep.int(x = colnames(x = x), times = length(x = i))
    }
  } else {
    slot.use <- unlist(x = lapply(
      X = c('assays', 'reductions', 'graphs', 'neighbors', 'commands', 'workflows'),
      FUN = function(s) {
        if (any(i %in% names(x = slot(object = x, name = s)))) {
          return(s)
        }
        return(NULL)
      }
    ))
    if (is.null(x = slot.use)) {
      stop("Cannot find '", i, "' in this Seurat object")
    }
    data.return <- slot(object = x, name = slot.use)[[i]]
  }
  return(data.return)
}

#' @export
#' @method as.logical JackStrawData
#'
as.logical.JackStrawData <- function(x, ...) {
  empP <- JS(object = x, slot = 'empirical')
  return(!(all(dim(x = empP) == 0) || all(is.na(x = empP))))
}

#' @export
#' @method dim Assay
#'
dim.Assay <- function(x) {
  return(dim(x = GetAssayData(object = x)))
}

#' @export
#' @method dim DimReduc
#'
dim.DimReduc <- function(x) {
  return(c(
    nrow(x = Loadings(object = x)),
    nrow(x = Embeddings(object = x))
  ))
}

#' @export
#' @method dim Seurat
#'
dim.Seurat <- function(x) {
  return(dim(x = GetAssay(object = x)))
}

#' @export
#' @method dimnames Assay
#'
dimnames.Assay <- function(x) {
  return(dimnames(x = GetAssayData(object = x)))
}

#' @export
#' @method dimnames DimReduc
#'
dimnames.DimReduc <- function(x) {
  return(list(
    rownames(x = Loadings(object = x)),
    rownames(x = Embeddings(object = x)))
  )
}

#' @export
#' @method dimnames Seurat
#'
dimnames.Seurat <- function(x) {
  return(dimnames(x = GetAssay(object = x)))
}

#' @importFrom ggplot2 ggplot aes
#' @export
#'
ggplot.DimReduc <- function(
  data = NULL,
  type = 'embeddings',
  colors = NULL,
  projected = NULL,
  rows.use = NULL,
  pt.size = NULL,
  pt.shape = NULL,
  mapping = aes(), ...,
  environment = parent.frame()
) {
  data.plot <- if (type == 'embeddings') {
    Embeddings(object = data)
  } else if (type == 'loadings') {
    Loadings(object = data, projected = projected)
  } else {
    stop("'type' must be either 'embeddings' or 'loadings'")
  }
  data.plot <- as.data.frame(x = data.plot)
  if (!is.null(x = colors)) {
    if (!is.null(x = names(x = colors))) {
      colors <- colors[colnames(x = data)]
    }
    data.plot$color <- colors
  }
  if (!is.null(x = pt.size)) {
    data.plot$size <- pt.size
  }
  if (!is.null(x = pt.shape)) {
    data.plot$shape <- pt.shape
  }
  if (!is.null(x = rows.use)) {
    data.plot <- data.plot[rows.use, , drop = FALSE]
  }
  return(ggplot(
    data = data.plot,
    mapping = mapping,
    ...,
    environment = environment
  ))
}

#' @export
#' @method length DimReduc
#'
length.DimReduc <- function(x) {
  return(ncol(x = Embeddings(object = x)))
}

#' @rdname Idents
#' @export
#' @method levels Seurat
#'
levels.Seurat <- function(x) {
  return(levels(x = Idents(object = x)))
}

#' @rdname merge.Seurat
#' @export
#' @method merge Assay
#'
merge.Assay <- function(
  x = NULL,
  y = NULL,
  add.cell.ids = NULL,
  min.cells = 0,
  min.features = 0,
  merge.data = TRUE
) {
  assays <- c(x, y)
  if (!is.null(x = add.cell.ids)) {
    for (i in 1:length(assays)) {
      assays[[i]] <- RenameCells(object = assays[[i]], new.names = add.cell.ids[i])
    }
  }
  # Merge the counts
  merged.counts <- GetAssayData(object = assays[[1]], slot = "counts")
  for (i in 2:length(x = assays)) {
    merged.counts <- RowMergeSparseMatrices(
      mat1 = merged.counts,
      mat2 = GetAssayData(object = assays[[i]], slot = "counts")
    )
  }
  combined.assay <- CreateAssayObject(
    counts = merged.counts,
    min.cells = min.cells,
    min.features = min.features
  )
  if (merge.data) {
    merged.data <- GetAssayData(object = assays[[1]], slot = "data")
    for (i in 2:length(x = assays)) {
      merged.data <- RowMergeSparseMatrices(
        mat1 = merged.data,
        mat2 = GetAssayData(object = assays[[i]], slot = "data")
      )
    }
    # only keep cells that made it through counts filtering params
    merged.data <- merged.data[, colnames(combined.assay)]
    combined.assay <- SetAssayData(
      object = combined.assay,
      slot = "data",
      new.data = merged.data
    )
  }
  return(combined.assay)
}

#' Merge Seurat Objects
#'
#' Merge two or more objects.
#'
#' When merging Seurat objects, the merge procedure will merge the Assay level
#' counts and potentially the data slots (depending on the merge.data parameter).
#' It will also merge the cell-level meta data that was stored with each object
#' and preserve the cell identities that were active in the objects pre-merge.
#' The merge will not preserve reductions, graphs or logged commands that were
#' present in the original objects.
#'
#' @param x Object
#' @param y Object (or a list of multiple objects)
#' @param add.cell.ids A character vector of length(x = c(x, y)). Appends the
#' corresponding values to the start of each objects' cell names.
#' @param merge.data Merge the data slots instead of just merging the counts
#' (which requires renormalization). This is recommended if the same normalization
#' approach was applied to all objects.
#' @inheritParams CreateSeuratObject
#'
#' @return Merged object
#'
#' @rdname merge.Seurat
#' @aliases merge
#'
#' @export
#' @method merge Seurat
#'
merge.Seurat <- function(
  x = NULL,
  y = NULL,
  add.cell.ids = NULL,
  merge.data = TRUE,
  project = "SeuratProject",
  min.cells = 0,
  min.features = 0
) {
  objects <- c(x, y)
  if (!is.null(add.cell.ids)) {
    if (length(x = add.cell.ids) != length(x = objects)) {
      stop("Please provide a cell identifier for each object provided to merge")
    }
    for(i in 1:length(objects)) {
      objects[[i]] <- RenameCells(object = objects[[i]], add.cell.id = add.cell.ids[i])
    }
  }
  assays.to.merge <- c()
  for (i in 1:length(objects)) {
    assays.to.merge <- c(assays.to.merge, FilterObjects(object = objects[[i]], classes.keep = "Assay"))
  }
  assays.to.merge <- names(which(x = table(... = assays.to.merge) == length(x = objects)))
  combined.assays <- list()
  for (assay in assays.to.merge) {
    assay1 <- objects[[1]][[assay]]
    assay2 <- list()
    for(i in 2:length(objects)) {
      assay2[[i-1]] <- objects[[i]][[assay]]
    }
    combined.assays[[assay]] <- merge(
      x = assay1,
      y = assay2,
      merge.data = merge.data,
      min.cells = min.cells,
      min.features = min.features
    )
  }
  # Merge the meta.data
  # get rid of nCount_ and nFeature_*
  combined.meta.data <- data.frame(row.names = colnames(combined.assays[[1]]))
  new.idents <- c()
  for (object in objects) {
    old.meta.data <- object[[]]
    old.meta.data[, which(grepl(pattern = "nCount_", x = colnames(old.meta.data)))] <- NULL
    old.meta.data[, which(grepl(pattern = "nFeature_", x = colnames(old.meta.data)))] <- NULL
    if (any(!colnames(x = old.meta.data) %in% colnames(combined.meta.data))) {
      cols.to.add <- colnames(x = old.meta.data)[!colnames(x = old.meta.data) %in% colnames(combined.meta.data)]
      combined.meta.data[, cols.to.add] <- NA
    }
    # unfactorize any factor columns
    i <- sapply(X = old.meta.data, FUN = is.factor)
    old.meta.data[i] <- lapply(old.meta.data[i], as.vector)
    combined.meta.data[rownames(old.meta.data), colnames(old.meta.data)] <- old.meta.data
    new.idents <- c(new.idents, as.vector(Idents(object = object)))
  }
  names(new.idents) <- rownames(combined.meta.data)
  new.idents <- factor(new.idents)
  if (DefaultAssay(object = x) %in% assays.to.merge){
    new.default.assay <- DefaultAssay(object = x)
  } else if (DefaultAssay(object = y) %in% assays.to.merge) {
    new.default.assay <- DefaultAssay(object = y)
  } else {
    new.default.assay <- assays.to.merge[1]
  }
  merged.object <- new(
    Class = 'Seurat',
    assays = combined.assays,
    meta.data = combined.meta.data,
    active.assay = new.default.assay,
    active.ident = new.idents,
    project.name = project,
    version = packageVersion(pkg = 'Seurat')
  )
  for (assay in assays.to.merge) {
    merged.object[[paste('nFeature', assay, sep = '_')]] <-
      Matrix::colSums(x = GetAssayData(
        object = merged.object,
        assay = assay, slot = "counts") > 0)
    merged.object[[paste('nCount', assay, sep = '_')]] <-
      Matrix::colSums(x = GetAssayData(
        object = merged.object,
        assay = assay,
        slot = 'counts'))
  }
  return(merged.object)
}

#' @export
#' @method names Seurat
#'
names.Seurat <- function(x) {
  return(unlist(
    x = lapply(
      X = c('assays', 'reductions', 'graphs'),
      FUN = function(n) {
        return(names(x = slot(object = x, name = n)))
      }
    ),
    use.names = FALSE
  ))
}

#' @export
#' @method subset Assay
#'
subset.Assay <- function(x, cells = NULL, features = NULL) {
  cells <- cells %||% colnames(x = x)
  features <- features %||% rownames(x = x)
  if (all(sapply(X = list(features, cells), FUN = length) == dim(x = x))) {
    return(x)
  }
  if (is.numeric(x = features)) {
    features <- rownames(x = x)[features]
  }
  features <- gsub(
    pattern = paste0('^', Key(object = x)),
    replacement = '',
    x = features
  )
  features <- intersect(x = rownames(x = x), y = features)
  if (length(x = features) == 0) {
    stop("Cannot find features provided")
  }
  if (ncol(x = GetAssayData(object = x, slot = 'counts')) == ncol(x = x)) {
    slot(object = x, name = "counts") <- GetAssayData(object = x, slot = "counts")[features, cells, drop = FALSE]
  }
  slot(object = x, name = "data") <- GetAssayData(object = x, slot = "data")[features, cells, drop = FALSE]
  cells.scaled <- colnames(x = GetAssayData(object = x, slot = "scale.data"))
  cells.scaled <- cells.scaled[cells.scaled %in% cells]
  features.scaled <- rownames(x = GetAssayData(object = x, slot = 'scale.data'))
  features.scaled <- features.scaled[features.scaled %in% features]
  slot(object = x, name = "scale.data") <- if (length(x = cells.scaled) > 0 && length(x = features.scaled) > 0) {
    GetAssayData(object = x, slot = "scale.data")[features.scaled, cells.scaled, drop = FALSE]
  } else {
    new(Class = 'matrix')
  }
  VariableFeatures(object = x) <- VariableFeatures(object = x)[VariableFeatures(object = x) %in% features]
  slot(object = x, name = 'meta.features') <- x[[]][features, , drop = FALSE]
  return(x)
}

#' @export
#' @method subset DimReduc
#'
subset.DimReduc <- function(x, cells = NULL, features = NULL) {
  cells <- colnames(x = x) %iff% cells %||% colnames(x = x)
  features <- rownames(x = x) %iff% features %||% rownames(x = x)
  if (all(sapply(X = list(features, cells), FUN = length) == dim(x = x))) {
    return(x)
  }
  slot(object = x, name = 'cell.embeddings') <- if (is.null(x = cells)) {
    new(Class = 'matrix')
  } else {
    if (is.numeric(x = cells)) {
      cells <- colnames(x = x)[cells]
    }
    cells <- intersect(x = cells, y = colnames(x = x))
    if (length(x = cells) == 0) {
      stop("Cannot find cell provided", call. = FALSE)
    }
    x[[cells, , drop = FALSE]]
  }
  slot(object = x, name = 'feature.loadings') <- if (is.null(x = features)) {
    new(Class = 'matrix')
  } else {
    if (is.numeric(x = features)) {
      features <- rownames(x = x)[features]
    }
    features.loadings <- intersect(
      x = rownames(x = Loadings(object = x, projected = FALSE)),
      y = features
    )
    if (length(x = features.loadings) == 0) {
      stop("Cannot find features provided", call. = FALSE)
    }
    Loadings(object = x, projected = FALSE)[features.loadings, , drop = FALSE]
  }
  slot(object = x, name = 'feature.loadings.projected') <- if (is.null(x = features) && !Projected(object = x)) {
    new(Class = 'matrix')
  } else {
    features.projected <- intersect(
      x = rownames(x = Loadings(object = x, projected = TRUE)),
      y = features
    )
    if (length(x = features.projected) == 0) {
      stop("Cannot find features provided", call. = FALSE)
    }
    Loadings(object = x, projected = TRUE)[features.projected, , drop = FALSE]
  }
  slot(object = x, name = 'jackstraw') <- new(Class = 'JackStrawData')
  return(x)
}

#' Subset a Seurat object
#'
#' @param x Seurat object to be subsetted
#' @param subset Logical expression indicating features/variables to keep
#' @param select A vector of cells, identity classes, or features to keep
#' @param ... Arguments passed to \code{WhichCells}
#'
#' @return A subsetted Seurat object
#'
#' @rdname subset.Seurat
#' @aliases subset SubsetData
#' @seealso \code{\link{base::subset}} \code{\link{WhichCells}}
#'
#' @export
#' @method subset Seurat
#'
#' @examples
#' subset(x = pbmc_small, subset = MS4A1 > 7)
#' subset(x = pbmc_small, select = VariableFeatures(object = pbmc_small))
#'
subset.Seurat <- function(x, subset, select = NULL, ...) {
  if (!missing(x = subset)) {
    subset <- deparse(expr = substitute(expr = subset))
  }
  cells <- select %iff% if (any(select %in% colnames(x = x))) {
    select[select %in% colnames(x = x)]
  } else {
    NULL
  }
  idents <- select %iff% if (any(select %in% levels(x = Idents(object = x)))) {
    select[select %in% levels(x = Idents(object = x))]
  } else {
    NULL
  }
  cells <- WhichCells(
    object = x,
    cells = cells,
    idents = idents,
    expression = subset,
    ...
  )
  if (length(x = cells) == 0) {
    stop("No cells found", call. = FALSE)
  }
  assays <- FilterObjects(object = x, classes.keep = 'Assay')
  assay.keys <- sapply(
    X = assays,
    FUN = function(i) {
      return(Key(object = x[[i]]))
    }
  )
  assay.features <- unlist(
    x = lapply(
      X = assays,
      FUN = function(i) {
        return(rownames(x = x[[i]]))
      }
    ),
    use.names = FALSE
  )
  assay.features <- unique(x = assay.features)
  key.patterns <- paste0('^', assay.keys, collapse = '|')
  features <- select %iff% if (any(grepl(pattern = key.patterns, x = select, perl = TRUE) | select %in% assay.features)) {
    select[grepl(pattern = key.patterns, x = select, perl = TRUE) | select %in% assay.features]
  } else {
    NULL
  }
  for (assay in assays) {
    assay.features <- features %||% rownames(x = x[[assay]])
    slot(object = x, name = 'assays')[[assay]] <- tryCatch(
      expr = subset.Assay(x = x[[assay]], cells = cells, features = assay.features),
      error = function(e) {
        return(NULL)
      }
    )
  }
  slot(object = x, name = 'assays') <- Filter(
    f = Negate(f = is.null),
    x = slot(object = x, name = 'assays')
  )
  if (length(x = FilterObjects(object = x, classes.keep = 'Assay')) == 0 || is.null(x = x[[DefaultAssay(object = x)]])) {
    stop("Cannot delete the default assay", call. = FALSE)
  }
  for (dimreduc in FilterObjects(object = x, classes.keep = 'DimReduc')) {
    x[[dimreduc]] <- tryCatch(
      expr = subset.DimReduc(x = x[[dimreduc]], cells = cells, features = features),
      error = function(e) {
        return(NULL)
      }
    )
  }
  slot(object = x, name = 'meta.data') <- slot(object = x, name = 'meta.data')[cells, , drop = FALSE]
  slot(object = x, name = 'graphs') <- list()
  Idents(object = x) <- Idents(object = x)[cells]
  return(x)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# S4 methods
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

setMethod(
  f = '[<-',
  signature = c('x' = 'Seurat'),
  definition = function(x, i, ..., value) {
    .Defunct(
      new = '[[<-',
      package = 'Seurat',
      msg = "When setting meta data, please use '[[<-' instead"
    )
  }
)

setMethod(
  f = '[[<-',
  signature = c('x' = 'Assay'),
  definition = function(x, i, ..., value) {
    meta.data <- x[[]]
    feature.names <- rownames(x = meta.data)
    if (length(x = i) > 1) {
      value <- rep_len(x = value, length.out = length(x = i))
      for (index in 1:length(x = i)) {
        meta.data[i[index]] <- value[index]
      }
    } else {
      if (length(x = intersect(x = names(x = value), y = feature.names)) > 0) {
        meta.data[, i] <- value[feature.names]
      } else if (length(x = value) %in% c(nrow(x = meta.data), 1)) {
        meta.data[, i] <- value
      } else {
        stop("Cannot add more or fewer meta.features information without values being named with feature names")
      }
    }
    slot(object = x, name = 'meta.features') <- meta.data
    return(x)
  }
)

setMethod( # because R doesn't allow S3-style [[<- for S4 classes
  f = '[[<-',
  signature = c('x' = 'Seurat'),
  definition = function(x, i, ..., value) {
    if (!is.character(x = i)) {
      stop("'i' must be a character", call. = FALSE)
    }
    if (is.null(x = value)) {
      slot.use <- FindObject(object = x, name = i)
      if (is.null(x = slot.use)) {
        stop("Cannot find object ", i, call. = FALSE)
      }
      if (i == DefaultAssay(object = x)) {
        stop("Cannot delete the default assay", call. = FALSE)
      }
    }
    slot.use <- switch(
      EXPR = as.character(x = class(x = value)),
      'Assay' = {
        if (all(colnames(x = value) %in% colnames(x = x)) && !all(colnames(x = value) == colnames(x = x))) {
          for (slot in c('counts', 'data', 'scale.data')) {
            assay.data <- GetAssayData(object = value, slot = slot)
            # if (nrow(x = assay.data) > 0) {
            if (!IsMatrixEmpty(x = assay.data)) {
              assay.data <- assay.data[, colnames(x = x), drop = FALSE]
            }
            value <- SetAssayData(object = value, slot = slot, new.data = assay.data)
          }
        }
        'assays'
      },
      'Graph' = 'graphs',
      'DimReduc' = {
        if (is.null(x = DefaultAssay(object = value))) {
          stop("Cannot add a DimReduc without an assay associated with it", call. = FALSE)
        }
        'reductions'
      },
      'SeuratCommand' = 'commands',
      'SeuratWorkflow' = 'workflows',
      'NULL' = slot.use,
      # stop("Unknown object type: ", class(x = value))
      'meta.data'
    )
    if (slot.use == 'meta.data') {
      meta.data <- x[[]]
      cell.names <- rownames(x = meta.data)
      if (is.data.frame(x = value) && !is.null(x = rownames(x = value))) {
        value <- value[cell.names, , drop = FALSE]
      }
      if (length(x = i) > 1) {
        value <- rep_len(x = value, length.out = length(x = i))
        for (index in 1:length(x = i)) {
          meta.data[i[index]] <- value[index]
        }
      } else {
        if (length(x = intersect(x = names(x = value), y = cell.names)) > 0) {
          meta.data[, i] <- value[cell.names]
        } else if (length(x = value) %in% c(nrow(x = meta.data), 1)) {
          meta.data[, i] <- value
        } else {
          stop("Cannot add more or fewer cell meta.data information without values being named with cell names")
        }
      }
      slot(object = x, name = 'meta.data') <- meta.data
    } else {
      if (!(class(x = value) %in% c('SeuratCommand', 'NULL')) && !all(colnames(x = value) == colnames(x = x))) {
        stop("All cells in the object being added must match the cells in this object", call. = FALSE)
      }
      if (!is.null(x = FindObject(object = x, name = i)) && !(class(x = value) %in% c(class(x = x[[i]]), 'NULL'))) {
        stop(
          "This object already contains ",
          i,
          " as a",
          ifelse(
            test = tolower(x = substring(text = class(x = x[[i]]), first = 1, last = 1)) %in% c('a', 'e', 'i', 'o', 'u'),
            yes = 'n ',
            no = ' '
          ),
          class(x = x[[i]]),
          "; duplicate names are not allowed",
          call. = FALSE
        )
      }
      if (class(x = value) %in% c('Assay', 'DimReduc') && length(x = Key(object = value)) == 0) {
        Key(object = value) <- paste0(tolower(x = i), '_')
      }
      slot(object = x, name = slot.use)[[i]] <- value
      slot(object = x, name = slot.use) <- Filter(
        f = Negate(f = is.null),
        x = slot(object = x, name = slot.use)
      )
    }
    gc(verbose = FALSE)
    return(x)
  }
)

setMethod(
  f = 'colMeans',
  signature = c('x' = 'Assay'),
  definition = function(x, na.rm = FALSE, dims = 1, ...) {
    return(Matrix::colMeans(
      x = GetAssayData(object = x),
      na.rm = na.rm,
      dims = dims,
      ...
    ))
  }
)

setMethod(
  f = 'colMeans',
  signature = c('x' = 'Seurat'),
  definition = function(x, na.rm = FALSE, dims = 1, ...) {
    return(colMeans(
      x = GetAssayData(object = x),
      na.rm = na.rm,
      dims = dims,
      ...
    ))
  }
)

setMethod(
  f = 'colSums',
  signature = c('x' = 'Assay'),
  definition = function(x, na.rm = FALSE, dims = 1, ...) {
    return(Matrix::colSums(
      x = GetAssayData(object = x),
      na.rm = na.rm,
      dims = dims,
      ...
    ))
  }
)

setMethod(
  f = 'colSums',
  signature = c('x' = 'Seurat'),
  definition = function(x, na.rm = FALSE, dims = 1, ...) {
    return(Matrix::colSums(
      x = GetAssayData(object = x),
      na.rm = na.rm,
      dims = dims,
      ...
    ))
  }
)

setMethod(
  f = 'rowMeans',
  signature = c('x' = 'Assay'),
  definition = function(x, na.rm = FALSE, dims = 1, ...) {
    return(Matrix::rowMeans(
      x = GetAssayData(object = x),
      na.rm = na.rm,
      dims = dims,
      ...
    ))
  }
)

setMethod(
  f = 'rowMeans',
  signature = c('x' = 'Seurat'),
  definition = function(x, na.rm = FALSE, dims = 1, ...) {
    return(Matrix::rowMeans(
      x = GetAssayData(object = x),
      na.rm = na.rm,
      dims = dims,
      ...
    ))
  }
)

setMethod(
  f = 'rowSums',
  signature = c('x' = 'Assay'),
  definition = function(x, na.rm = FALSE, dims = 1, ...) {
    return(Matrix::rowSums(
      x = GetAssayData(object = x),
      na.rm = na.rm,
      dims = dims,
      ...
    ))
  }
)

setMethod(
  f = 'rowSums',
  signature = c('x' = 'Seurat'),
  definition = function(x, na.rm = FALSE, dims = 1, ...) {
    return(Matrix::rowSums(
      x = GetAssayData(object = x),
      na.rm = na.rm,
      dims = dims,
      ...
    ))
  }
)

setMethod(
  f = 'show',
  signature = 'Assay',
  definition = function(object) {
    cat('Assay data with', nrow(x = object), 'features for', ncol(x = object), 'cells\n')
    if (length(x = VariableFeatures(object = object)) > 0) {
      top.ten <- head(x = VariableFeatures(object = object), n = 10L)
      cat(
        "Top",
        length(x = top.ten),
        paste0("variable feature", if (length(x = top.ten) > 1) {'s'}, ":\n"),
        paste(strwrap(x = paste(top.ten, collapse = ', ')), collapse = '\n '),
        '\n'
      )
    }
  }
)

setMethod(
  f = 'show',
  signature = 'DimReduc',
  definition = function(object) {
    cat(
      "A dimensional reduction object with key", Key(object = object), '\n',
      'Number of dimensions:', length(x = object), '\n',
      'Projected dimensional reduction calculated: ', Projected(object = object), '\n',
      'Jackstraw run:', as.logical(x = JS(object = object)), '\n',
      'Computed using assay:', DefaultAssay(object), '\n'
    )
  }
)

setMethod(
  f = 'show',
  signature = 'JackStrawData',
  definition = function(object) {
    # empp <- GetJS(object = object, slot = "empirical.p.values")
    empp <- object$empirical.p.values
    # scored <- GetJS(object = object, slot = "overall.p.values")
    scored <- object$overall.p.values
    cat(
      "A JackStrawData object simulated on",
      nrow(x = empp),
      "features for",
      ncol(x = empp),
      "dimensions.\n",
      "Scored for:",
      nrow(x = scored),
      "dimensions.\n"
    )
  }
)

setMethod(
  f = "show",
  signature = "Seurat",
  definition = function(object) {
    assays <- FilterObjects(object = object, classes.keep = 'Assay')
    nfeatures <- sum(vapply(
      X = assays,
      FUN = function(x) {
        return(nrow(x = object[[x]]))
      },
      FUN.VALUE = integer(length = 1L)
    ))
    num.assays <- length(x = assays)
    cat("An object of class", class(x = object), "\n")
    cat(
      nfeatures,
      'features across',
      ncol(x = object),
      'samples within',
      num.assays,
      ifelse(test = num.assays == 1, yes = 'assay', no = 'assays'),
      "\n"
    )
    cat("Active assay:", DefaultAssay(object = object))
    reductions <- FilterObjects(object = object, classes.keep = 'DimReduc')
    if (length(x = reductions) > 0) {
      cat(
        '\n',
        length(x = reductions),
        'dimensional',
        ifelse(test = length(x = reductions) == 1, yes = 'reduction', no = 'reductions'),
        'calculated:',
        strwrap(x = paste(reductions, collapse = ', '))
      )
    }
    cat('\n')
  }
)

setMethod(
  f = 'show',
  signature = 'SeuratCommand',
  definition = function(object) {
    params <- slot(object = object, name = "params")
    cat(
      "Command: ", slot(object = object, name = "call.string"), '\n',
      "Time: ", as.character(slot(object = object, name = "time.stamp")), '\n',
      sep = ""
    )
    for(p in 1:length(params)){
      cat(
       names(params[p]), ":", params[[p]], "\n"
      )
    }
  }
)

setMethod(
  f = 'show',
  signature = 'seurat',
  definition = function(object) {
    cat(
      "An old seurat object\n",
      nrow(x = object@data),
      'genes across',
      ncol(x = object@data),
      'samples\n'
    )
  }
)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Internal
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Get cell names grouped by identity class
#
# @param object A Seurat object
# @param cells A vector of cells to grouping to
#
# @return A named list where names are identity classes and values are vectors
# of cells beloning to that class
#
CellsByIdentities <- function(object, cells = NULL) {
  cells <- cells %||% colnames(x = object)
  cells <- intersect(x = cells, y = colnames(x = object))
  if (length(x = cells) == 0) {
    stop("Cannot find cells provided")
  }
  idents <- levels(x = object)
  cells.idents <- sapply(
    X = idents,
    FUN = function(i) {
      return(cells[as.vector(x = Idents(object = object)[cells]) == i])
    },
    simplify = FALSE,
    USE.NAMES = TRUE
  )
  return(cells.idents)
}

# Get the names of objects within a Seurat object that are of a certain class
#
# @param object A Seurat object
# @param classes.keep A vector of names of classes to get
#
# @return A vector with the names of objects within the Seurat object that are of class \code{classes.keep}
#
FilterObjects <- function(object, classes.keep = c('Assay', 'DimReduc')) {
  object.classes <- sapply(
    X = names(x = object),
    FUN = function(i) {
      return(class(x = object[[i]]))
    }
  )
  object.classes <- object.classes[object.classes %in% classes.keep]
  return(names(x = object.classes))
}

# Find the collection of an object within a Seurat object
#
# @param object A Seurat object
# @param name Name of object to find
#
# @return The collection (slot) of the object
#
FindObject <- function(object, name) {
  collections <- c('assays', 'graphs', 'neighbors', 'reductions', 'commands', 'workflows')
  object.names <- lapply(
    X = collections,
    FUN = function(x) {
      return(names(x = slot(object = object, name = x)))
    }
  )
  names(x = object.names) <- collections
  object.names <- Filter(f = Negate(f = is.null), x = object.names)
  for (i in names(x = object.names)) {
    if (name %in% names(x = slot(object = object, name = i))) {
      return(i)
    }
  }
  return(NULL)
}

# Prepare a Seurat function for a workflow run
#
# Checks dependencies (and runs them if necessary), and then reads in parameter values
#
PrepareWorkflow <- function(object, workflow.name) {
  CheckWorkflow(object = object, workflow.name = workflow.name)
  command.name <- as.character(deparse(sys.calls()[[sys.nframe()-1]]))
  command.name <- gsub(pattern = ".Seurat",replacement = "",x = command.name)
  command.name <- ExtractField(string = command.name, field = 1, delim = "\\(")
  depends <- slot(object = object[[workflow.name]], name = "depends")
  prereq.commands <- colnames(depends)[which(depends[command.name, ] == 1)]
  for(i in prereq.commands) {
    check.prereqs <- CheckWorkflowUpdate(
      object = object,
      workflow.name = workflow.name,
      command.name = i)
    if (check.prereqs) {
      # run the dependency
      workflow.name.quotes <- paste0('\"', workflow.name, '\"')
      new.cmd <- paste0("object <- ", i, "(object, workflow = ", workflow.name.quotes, ")")
      message(paste0("Updating ", i))
      message(new.cmd)
      eval(expr = parse(text = new.cmd))
    }
  }
  ReadWorkflowParams( object = object, workflow.name = workflow.name, depth = 2)
  return(object)
}

# Check to see if projected loadings have been set
#
# @param object a DimReduc object
#
# @return TRUE if proejcted loadings have been set, else FALSE
#
Projected <- function(object) {
  projected.dims <- dim(x = slot(object = object, name = 'feature.loadings.projected'))
  if (all(projected.dims == 1)) {
    return(!all(is.na(x = slot(object = object, name = 'feature.loadings.projected'))))
  }
  return(!all(projected.dims == 0))
}

# Read Seurat workflow parameters
#
# Reads parameters for a workflow and assigns them to the parent function call.
#
# @param depth is the depth of the call in the function stack. If called
# directly from, for example, ScaleData, then depth=1. If called indirectly
# from PrepareWorkflow, then depth=2
#
ReadWorkflowParams <- function(object, workflow.name, depth = 2) {
  CheckWorkflow(object = object, workflow.name = workflow.name)
  param.list <- names(formals(fun = sys.function(sys.parent(depth))))
  workflow.params <- slot(object = object[[workflow.name]], name = "params")

  # global variables
  to.set <- intersect(x = param.list, y = names(workflow.params$global))
  p.env <- parent.frame(depth)
  for(i in to.set) {
    assign(x = names(workflow.params$global[i]),
           value = workflow.params$global[[i]],
           envir = p.env)
  }

  # parameter-specific variables
  command.name <- as.character(deparse(sys.calls()[[sys.nframe()-depth]]))
  command.name <- gsub(pattern = ".Seurat",replacement = "",x = command.name)
  command.name <- ExtractField(string = command.name, field = 1, delim = "\\(")
  to.set <- intersect(x = param.list, y = names(workflow.params[[command.name]]))
  for(i in to.set) {
    assign(x = names(workflow.params[[command.name]][i]),
           value = workflow.params[[command.name]][[i]],
           envir = p.env)
  }

  # overwrite any arguments passed in on the command line
  argnames <- sys.call(which = depth)
  argList <- as.list(argnames[-1])
  args_ignore <- c("", "object", "workflow.name")
  args_use <- setdiff(x = names(argList), y = args_ignore)
  for(i in args_use) {
    if(as.character(unlist(argList[i])[[1]]) == "F") {
      arg.val <- FALSE
    } else if(as.character(unlist(argList[i])[[1]]) == "T") {
      arg.val <- TRUE
    } else {
      arg.val <- unlist(argList[i])[[1]]
    }
    assign(x = i, value = arg.val, envir = p.env)
  }
}

# Get the top
#
# @param data Data to pull the top from
# @param num Pull top \code{num}
# @param balanced Pull even amounts of from positive and negative values
#
# @return The top \code{num}
# @seealso \{code{\link{TopCells}}} \{code{\link{TopFeatures}}}
#
Top <- function(
  data,
  num,
  balanced
) {
  top <- if (balanced) {
    num <- round(x = num / 2)
    data <- data[order(data), , drop = FALSE]
    positive <- head(x = rownames(x = data), n = num)
    negative <- rev(x = tail(x = rownames(x = data), n = num))
    list(positive = positive, negative = negative)
  } else {
    data <- data[rev(x = order(abs(x = data))), , drop = FALSE]
    top <- head(x = rownames(x = data), n = num)
    top[order(data[top, ])]
  }
  return(top)
}

# Update Seurat assay
#
# @param old.assay Seurat2 assay
# @param assay Name to store for assay in new object
#
UpdateAssay <- function(old.assay, assay){
  cells <- colnames(old.assay@data)
  counts <- old.assay@raw.data
  data <- old.assay@data
  if (!inherits(x = counts, what = 'dgCMatrix')) {
    counts <- as(object = as.matrix(x = counts), Class = 'dgCMatrix')
  }
  if (!inherits(x = data, what = 'dgCMatrix')) {
    data <- as(object = as.matrix(x = data), Class = 'dgCMatrix')
  }
  new.assay <- new(
    Class = 'Assay',
    counts = counts[, cells],
    data = data,
    scale.data = old.assay@scale.data %||% matrix(),
    meta.features = data.frame(row.names = rownames(x = counts)),
    var.features = old.assay@var.genes,
    key = paste0(assay, "_")
  )
  return(new.assay)
}

# Update dimension reduction
#
# @param old.dr Seurat2 dimension reduction slot
#
UpdateDimReduction <- function(old.dr, assay){
  new.dr <- list()
  for(i in names(old.dr)){
    cell.embeddings <- old.dr[[i]]@cell.embeddings %||% matrix()
    feature.loadings <- old.dr[[i]]@gene.loadings %||% matrix()
    stdev <- old.dr[[i]]@sdev %||% numeric()
    misc <- old.dr[[i]]@misc %||% list()
    new.jackstraw <- UpdateJackstraw(old.dr[[i]]@jackstraw)

    new.dr[[i]] <- new(
      Class = 'DimReduc',
      cell.embeddings = as(cell.embeddings, 'matrix'),
      feature.loadings = as(feature.loadings, 'matrix'),
      assay.used = assay,
      stdev = as(stdev, 'numeric'),
      key = as(old.dr[[i]]@key, 'character'),
      jackstraw = new.jackstraw,
      misc = as(misc, 'list')
    )
  }
  return(new.dr)
}

# Update jackstraw
#
# @param old.jackstraw
#
UpdateJackstraw <- function(old.jackstraw) {
  if (is.null(old.jackstraw)) {
    new.jackstraw <- new(
      Class = 'JackStrawData',
      empirical.p.values = matrix(),
      fake.reduction.scores = matrix(),
      empirical.p.values.full = matrix(),
      overall.p.values = matrix()
    )
  } else {
    if(.hasSlot(old.jackstraw, 'overall.p.values')) {
      overall.p <- old.jackstraw@overall.p.values %||% matrix()
    } else {
      overall.p <- matrix()
    }
    new.jackstraw <- new(
      Class = 'JackStrawData',
      empirical.p.values = old.jackstraw@emperical.p.value %||% matrix(),
      fake.reduction.scores = old.jackstraw@fake.pc.scores %||% matrix(),
      empirical.p.values.full = old.jackstraw@emperical.p.value.full %||% matrix(),
      overall.p.values = overall.p
    )
  }
  return(new.jackstraw)
}

# UpdateWorkflow
#
# Updates a workflow object after a command is run, and makes sure timestamps are properly set.
#
UpdateWorkflow <- function(object, workflow.name, command.name = NULL) {
  CheckWorkflow(object = object,workflow.name = workflow.name)
  if (is.null(x = command.name)) {
    command.name <- as.character(x = deparse(expr = sys.calls()[[sys.nframe() - 1]]))
    command.name <- gsub(pattern = ".Seurat", replacement = "", x = command.name)
    command.name <- ExtractField(string = command.name, field = 1, delim = "\\(")
    command.name.seurat <- intersect(
      x = c(command.name, paste0(command.name, ".", DefaultAssay(object = object))),
      names(x = object)
    )
  } else {
    command.name.seurat <- command.name
    command.name <- ExtractField(string = command.name, field = 1, delim = "\\.")
  }
  #TODO - Deal with Assay better
  seurat.timestamp <- Sys.time()
  if (length(x = command.name) == 1) {
    seurat.timestamp <- slot(
      object = object[[command.name.seurat]],
      name = "time.stamp"
    )
  }
  object <- TouchWorkflow(
    object = object,
    workflow.name = workflow.name,
    command.name = command.name,
    time.stamp = seurat.timestamp
  )
  return(object)
}

# Validates the workflow file that was provided.
#
# @param config results of reading in the ini file
#
# @return No return
#
ValidateWorkflowFile <- function(config) {
  sections <- names(x = config)
  if (!'dependencies' %in% sections) {
    stop("Workflow file is missing the dependencies section")
  }
}
