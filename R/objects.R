#' @include generics.R
#' @importFrom Rcpp evalCpp
#' @importFrom Matrix colSums rowSums colMeans rowMeans
#' @importFrom methods setClass setOldClass setClassUnion slot slot<- setMethod new
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
#' @importClassesFrom Matrix dgCMatrix
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
#' @importClassesFrom Matrix dgCMatrix
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
#' cluster_letters <- LETTERS[pbmc_small@ident]
#' pbmc_small <- AddMetaData(
#'   object = pbmc_small,
#'   metadata = cluster_letters,
#'   col.name = 'letter.idents'
#' )
#' head(x = pbmc_small@meta.data)
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
  meta.order <- match(rownames(object[]), rownames(metadata))
  meta.add <- metadata[meta.order, ]
  if (all(is.null(x = meta.add))) {
    stop("Metadata provided doesn't match the cells in this object")
  }
  slot(object = object, name = "meta.data")[, cols.add] <- meta.add
  return(object)
}

#' Create an Assay object
#'
#' Create an Assay object from a feature (e.g. gene) expression matrix. The expected format of the
#' input matrix is features x cells.
#'
#' Non-unique cell or feature names are not allowed. Please make unique before calling this function.
#'
#' @param counts Unnormalized data such as raw counts or TPMs
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
  min.cells = 0,
  min.features = 0
) {
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
  num.features <- Matrix::colSums(x = counts > 0)
  counts <- counts[, which(x = num.features > min.features)]
  # filter genes on the number of cells expressing
  if (min.cells > 0) {
    num.cells <- rowSums(x = counts > 0)
    counts <- counts[which(x = num.cells >= min.cells), ]
  }
  # Initialize meta.features
  init.meta.features <- data.frame(row.names = rownames(x = counts))
  assay <- new(
    Class = 'Assay',
    counts = counts,
    data = counts,
    meta.features = init.meta.features
  )
  return(assay)
}

#' Create a DimReduc object
#'
#' @param cell.embeddings ...
#' @param feature.loadings ...
#' @param feature.loadings.projected ...
#' @param assay.used ...
#' @param stdev ...
#' @param key ...
#' @param jackstraw ...
#' @param misc ...
#' @param ... Ignored for now
#'
#' @export
#'
CreateDimReducObject <- function(
  cell.embeddings = matrix(),
  feature.loadings = matrix(),
  feature.loadings.projected = matrix(),
  assay.used = NULL,
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
    cell.embeddings = cell.embeddings,
    feature.loadings = feature.loadings,
    feature.loadings.projected = feature.loadings.projected,
    assay.used = assay.used,
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
#' @param assay.use Name of the assay corresponding to the initial input data.
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
  assay.use = 'RNA',
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
  Key(object = assay.data) <- paste0(tolower(x = assay.use), '_')
  assay.list <- list(assay.data)
  names(x = assay.list) <- assay.use
  init.meta.data <- data.frame(row.names = colnames(x = assay.list[[assay.use]]))
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
    active.assay = assay.use,
    active.ident = idents,
    project.name = project,
    version = packageVersion(pkg = 'Seurat')
  )
  object['orig.ident'] <- idents
  # Calculate nUMI and nFeature
  object['nUMI'] <- Matrix::colSums(x = object)
  object[paste('nFeature', assay.use, sep = '_')] <- Matrix::colSums(counts > 0)
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
#' @param vars.fetch List of all variables to fetch
#' @param cells.use Cells to collect data for (default is all cells)
#' @param slot Slot to pull feature data for
#'
#' @return A data frame with cells as rows and cellular data as columns
#'
#' @export
#'
#' @examples
#' pc1 <- FetchData(object = pbmc_small, vars.all = 'PC1')
#' head(x = pc1)
#'
FetchData <- function(object, vars.fetch, cells.use = NULL, slot = 'data') {
  cells.use <- cells.use %||% colnames(x = object)
  if (is.numeric(x = cells.use)) {
    cells.use <- colnames(x = object)[cells.use]
  }
  objects.use <- FilterObjects(object = object)
  object.keys <- sapply(X = objects.use, FUN = function(i) {return(Key(object[[i]]))})
  keyed.vars <- lapply(
    X = object.keys,
    FUN = function(key) {
      if (length(x = key) == 0) {
        return(integer(length = 0L))
      }
      return(grep(pattern = paste0('^', key), x = vars.fetch))
    }
  )
  keyed.vars <- Filter(f = length, x = keyed.vars)
  data.fetched <- lapply(
    X = names(x = keyed.vars),
    FUN = function(x) {
      vars.use <- vars.fetch[keyed.vars[[x]]]
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
            object[[x]][[cells.use, vars.use, drop = FALSE]]
          } else {
            NULL
          }
        },
        'Assay' = {
          vars.use <- gsub(pattern = paste0('^', key.use), replacement = '', x = vars.use)
          data.vars <- t(x = as.matrix(x = GetAssayData(
            object = object,
            slot = slot,
            assay.use = x
          )[vars.use, cells.use, drop = FALSE]))
          colnames(x = data.vars) <- paste0(key.use, vars.use)
          data.vars
        }
      )
      data.return <- as.list(x = as.data.frame(x = data.return))
      return(data.return)
    }
  )
  data.fetched <- unlist(x = data.fetched, recursive = FALSE)
  meta.vars <- vars.fetch[vars.fetch %in% colnames(x = object[])]
  data.fetched <- c(data.fetched, object[meta.vars][cells.use, , drop = FALSE])
  default.vars <- vars.fetch[vars.fetch %in% rownames(x = object)]
  data.fetched <- c(
    data.fetched,
    as.data.frame(x = t(x = as.matrix(x = GetAssayData(
      object = object,
      slot = slot
    )[default.vars, cells.use, drop = FALSE])))
  )
  vars.fetched <- names(x = data.fetched)
  vars.missing <- setdiff(x = vars.fetch, y = vars.fetched)
  m2 <- if (length(x = vars.missing) > 10) {
    paste0(' (10 out of ', length(x = vars.missing), ' shown)')
  } else {
    ''
  }
  if (length(x = vars.missing) == length(x = vars.fetch)) {
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
    # warning('T', msg, immediate. = TRUE)
  }
  data.fetched <- as.data.frame(
    x = data.fetched,
    row.names = cells.use,
    stringsAsFactors = FALSE
  )
  data.order <- na.omit(object = pmatch(
    x = vars.fetch,
    table = vars.fetched
  ))
  if (length(x = data.order) > 1) {
    data.fetched <- data.fetched[, data.order]
  }
  colnames(x = data.fetched) <- vars.fetch[vars.fetch %in% vars.fetched]
  return(data.fetched)
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
#' @param \dots Additional parameters to pass to SubsetData
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
    groupings <- FetchData(object = object, vars.fetch = split.by)[, 1]
  }
  groupings <- unique(x = as.character(x = groupings))
  obj.list <- list()
  for (i in groupings) {
    if (split.by == "ident") {
      obj.list[[i]] <- SubsetData(object = object, ident.use = i, ...)
    }
    else {
      obj.list[[i]] <- SubsetData(
        object = object,
        subset.name = split.by,
        accept.value = i, ...
      )
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
#' @param num.features Number of features to return
#' @param projected Use the full PCA (projected PCA). Default i s FALSE
#' @param do.balanced Return an equal number of features with both + and - scores.
#'
#' @return Returns a vector of features
#'
#' @export
#'
#' @examples
#' pbmc_small
#' TopFeatures(object = pbmc_small, dim.use = 1, reduction.type = "pca")
#' # After projection:
#' TopFeatures(object = pbmc_small, dim.use = 1, reduction.type = "pca", use.full = TRUE)
#'
TopFeatures <- function(
  object,
  dim.use = 1,
  num.features = 20,
  projected = FALSE,
  do.balanced = FALSE
) {
  loadings <- Loadings(object = object, projected = projected)[, dim.use, drop = FALSE]
  return(Top(
    data.use = loadings,
    dim.use = dim.use,
    num.use = num.features,
    do.balanced = do.balanced
  ))
}

#' Find cells with highest scores for a given dimensional reduction technique
#'
#' Return a list of genes with the strongest contribution to a set of components
#'
#' @param object DimReduc object
#' @param dim.use Dimension to use
#' @param num.cells Number of cells to return
#' @param do.balanced Return an equal number of cells with both + and - scores.
#'
#' @return Returns a vector of cells
#'
#' @export
#'
#' @examples
#' pbmc_small
#' head(TopCells(object = pbmc_small, reduction.type = "pca"))
#' # Can specify which dimension and how many cells to return
#' TopCells(object = pbmc_small, reduction.type = "pca", dim.use = 2, num.cells = 5)
#'
TopCells <- function(
  object,
  dim.use = 1,
  num.cells = 20,
  do.balanced = FALSE
) {
  embeddings <- Embeddings(object = object)[, dim.use, drop = FALSE]
  return(Top(
    data.use = embeddings,
    dim.use = dim.use,
    num.use = num.cells,
    do.balanced = do.balanced
  ))
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Methods for Seurat-defined generics
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

#' @describeIn DefaultAssay Get the name of the assay used to calculate this DimReduc
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
GetAssay.Seurat <- function(object, assay.use = NULL) {
  assay.use <- assay.use %||% DefaultAssay(object = object)
  if (!assay.use %in% names(x = slot(object = object, name = 'assays'))) {
    stop(paste0(assay.use, " is not an assay present in the given object. Available assays are: ",
                paste(names(x = slot(object = object, name = 'assays')), collapse = ", ")))
  }
  return(slot(object = object, name = 'assays')[[assay.use]])
}

#' @describeIn GetAssayData Get assay data for an Assay object
#' @export
#' @method GetAssayData Assay
#'
GetAssayData.Assay <- function(object, slot = 'data') {
  return(slot(object = object, name = slot))
}

#' @param assay.use Name of assay to pull data from
#'
#' @describeIn GetAssayData Get assay data from a Seurat object
#' @export
#' @method GetAssayData Seurat
#'
GetAssayData.Seurat <- function(object, slot = 'data', assay.use = NULL, ...) {
  assay.use <- assay.use %||% DefaultAssay(object = object)
  return(GetAssayData(
    object = GetAssay(object = object, assay.use = assay.use),
    slot = slot
  ))
}

#' @describeIn HVFInfo Get highly variable feature information from an Assay
#' @export
#' @method HVFInfo Assay
#'
HVFInfo.Assay <- function(object, ...) {
  return(object[[c('mean', 'dispersion', 'dispersion.scaled')]])
}

#' @param assay.use Name of assay to pull highly variable feature information for
#'
#' @describeIn HVFInfo Get highly variable feature information from a Seurat object
#' @export
#' @method HVFInfo Seurat
#'
HVFInfo.Seurat <- function(object, assay.use = NULL, ...) {
  assay.use <- assay.use %||% DefaultAssay(object = object)
  assay.data <- GetAssay(object = object, assay.use = assay.use)
  return(HVFInfo(object = assay.data))
}

#' @describeIn Idents Get the active identities of a Seurat object
#' @export
#' @method Idents Seurat
#'
Idents.Seurat <- function(object, ...) {
  return(slot(object = object, name = 'active.ident'))
}

#' @param cells.use Set cell identities for specific cells
#'
#' @describeIn Idents Set the active identities of a Seurat object
#' @export
#' @method Idents<- Seurat
#'
"Idents<-.Seurat" <- function(object, cells.use = NULL, ..., value) {
  cells.use <- cells.use %||% colnames(x = object)
  if (is.numeric(x = cells.use)) {
    cells.use <- colnames(x = object)[cells.use]
  }
  cells.use <- intersect(x = cells.use, y = colnames(x = object))
  cells.use <- match(x = cells.use, table = colnames(x = object))
  idents.new <- if (length(x = value) == 1 && value %in% colnames(x = object[])) {
    unlist(x = object[value], use.names = FALSE)[cells.use]
  } else {
    if (is.list(x = value)) {
      value <- unlist(x = value, use.names = FALSE)
    }
    rep_len(x = value, length.out = length(x = cells.use))
  }
  idents.new <- as.vector(x = idents.new)
  idents <- as.vector(x = Idents(object = object))
  idents[cells.use] <- idents.new
  idents <- factor(x = idents)
  names(x = idents) <- colnames(x = object)
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
  slot.use <- if (is.null(x = projected)) {
    projected.data <- slot(object = object, name = 'feature.loadings.projected')
    ifelse(
      test = all(is.na(x = projected.data)) && unique(x = dim(x = projected.data)) == 1,
      yes = 'feature.loadings',
      no = 'feature.loadings.projected'
    )
  } else if (projected) {
    'feature.loadings.projected'
  } else {
    'feature.loadings'
  }
  return(slot(object = object, name = slot.use))
}

#' @describeIn Loadings Add projected feature loadings to a DimReduc object
#' @export
#' @method Loadings<- DimReduc
#'
"Loadings<-.DimReduc" <- function(object, ..., value) {
  if (nrow(x = value) != length(x = object)) {
    stop("New feature loadings must have the dimensions as currently calculated")
  }
  slot(object = object, name = 'feature.loadings.projected') <- value
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

#' @param dims Number of dimensions to display
#' @param num.features Number of genes to display
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
  num.features = 20,
  projected = FALSE
) {
  loadings <- Loadings(object = object, projected = projected)
  num.features <- min(num.features, nrow(x = loadings))
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
      dim.use = dim,
      num.features = num.features * 2,
      projected = projected,
      do.balanced = TRUE
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
  Idents(object = object, cells.use = rownames(x = object[])) <- old.ids
  # rename the cell-level metadata
  old.meta.data <- object[]
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

#' @describeIn SetAssayData Set assay data for an Assay object
#' @export
#' @method SetAssayData Assay
#'
SetAssayData.Assay <- function(object, slot, new.data) {
  slots.use <- c('counts', 'data', 'scale.data')
  if (!slot %in% slots.use) {
    stop("'slot' must be one of ", paste(slots.use, collapse = ', '))
  }
  if (ncol(x = new.data) != ncol(x = object)) {
    stop("The new data doesn't have the same number of cells as the current data")
  }
  num.counts <- nrow(x = GetAssayData(object = object, slot = 'counts'))
  counts.names <- rownames(x = GetAssayData(object = object, slot = 'counts'))
  if (num.counts == 0) {
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
  new.features <- na.omit(object = match(x = counts.names, table = rownames(x = new.data)))
  #if (slot == 'scale.data' && nrow(x = new.data) > nrow(x = object)) {
  #  stop("Cannot add more features than present in current data")
  #} else if (slot != 'scale.data' && nrow(x = new.data) != nrow(x = object)) {
  #  stop("The new data doesn't have the same number of features as the current data")
  #}
  new.cells <- colnames(x = new.data)
  if (!all(new.cells %in% colnames(x = object))) {
    stop("All cell names must match current cell names")
  }
  slot(object = object, name = slot) <- new.data[new.features, colnames(x = object)]
  return(object)
}

#' @param assay.use Name of assay whose data should be set
#'
#' @describeIn SetAssayData Set assay data for an Assay object in a Seurat object
#' @export
#' @method SetAssayData Seurat
#'
SetAssayData.Seurat <- function(
  object,
  slot = 'data',
  new.data,
  assay.use = NULL,
  ...
) {
  assay.use <- assay.use %||% DefaultAssay(object = object)
  assay.data <- GetAssay(object = object, assay.use = assay.use)
  assay.data <- SetAssayData(object = assay.data, slot = slot, new.data = new.data)
  object[[assay.use]] <- assay.data
  return(object)
}

#' @describeIn Idents Stash cell identities of a Seurat object
#' @export
#' @method StashIdent Seurat
#'
StashIdent.Seurat <- function(object, save.name = 'orig.ident', ...) {
  object[save.name] <- Idents(object = object)
  return(object)
}

#' @describeIn Stdev Get the standard deviations from a DimReduc object
#' @export
#' @method Stdev DimReduc
#'
Stdev.DimReduc <- function(object) {
  return(slot(object = object, name = 'stdev'))
}

#' @param reduction.use Name of reduction to use
#'
#' @describeIn Stdev Get the standard deviations of a dimensional reduction from a Seurat object
#' @export
#' @method Stdev Seurat
#'
Stdev.Seurat <- function(object, reduction.use, ...) {
  return(Stdev(object = object[[reduction.use]]))
}

#' @describeIn SubsetData Subset an Assay object
#' @export
#' @method SubsetData Assay
#'
SubsetData.Assay <- function(
  object,
  cells.use = NULL,
  subset.name = NULL,
  low.threshold = -Inf,
  high.threshold = Inf,
  accept.value = NULL,
  do.clean = FALSE,
  ...
) {
  cells.use <- cells.use %||% colnames(x = object)
  cells.use <- WhichCells(
    object = object,
    cells.use = cells.use,
    subset.name = subset.name,
    low.threshold = low.threshold,
    high.threshold = high.threshold,
    accept.value = accept.value,
    ...
  )
  slot(object = object, name = "counts") <- GetAssayData(object = object, slot = "counts")[, cells.use]
  slot(object = object, name = "data") <- GetAssayData(object = object, slot = "data")[, cells.use]
  cells.scaled <- colnames(x = GetAssayData(object = object, slot = "scale.data"))
  cells.scaled <- cells.scaled[cells.scaled %in% cells.use]
  if (length(x = cells.scaled) > 0) {
    slot(object = object, name = "scale.data") <- GetAssayData(object = object, slot = "scale.data")[, cells.use]
  }
  return(object)
}

#' @param assay.use Assay to subset on
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
  assay.use = NULL,
  cells.use = NULL,
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
  assay.use <- assay.use %||% DefaultAssay(object = object)
  cells.use <- WhichCells(
    object = object,
    assay.use = assay.use,
    ident = ident.use,
    ident.remove = ident.remove,
    subset.name = subset.name,
    cells.use = cells.use,
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
      cells.use = cells.use
    )
  }
  # Subset all the DimReducs
  drs <- FilterObjects(object = object, classes.keep = 'DimReduc')
  for (dr in drs) {
    object[[dr]] <- CreateDimReducObject(
      cell.embeddings = Embeddings(object = object[[dr]])[cells.use, ],
      feature.loadings = Loadings(object = object[[dr]], projected = FALSE),
      feature.loadings.projected = Loadings(object = object[[dr]], projected = TRUE),
      assay.used = DefaultAssay(object = object[[dr]]),
      # assay.used = slot(object = object[[dr]], name = "assay.used"),
      # stdev = slot(object = object[[dr]], name = "stdev"),
      stdev = Stdev(object = object[[dr]]),
      key = Key(object = object[[dr]]),
      jackstraw = slot(object = object[[dr]], name = "jackstraw"),
      misc = slot(object[[dr]], name = "misc")
    )
  }
  slot(object = object, name = "active.ident") <- Idents(object = object)[cells.use]
  slot(object = object, name = "meta.data") <- slot(object = object, name = "meta.data")[cells.use, ]
  return(object)
}

#' @describeIn VariableFeatures Get the variable features of an assay object
#' @export
#' @method VariableFeatures Assay
#'
VariableFeatures.Assay <- function(object, ...) {
  return(slot(object = object, name = 'var.features'))
}

#' @param assay.use Name of assay to pull variable features for
#'
#' @describeIn VariableFeatures Get the variable features of a Seurat object
#' @export
#' @method VariableFeatures Seurat
#'
VariableFeatures.Seurat <- function(object, assay.use = NULL, ...) {
  assay.use <- assay.use %||% DefaultAssay(object = object)
  assay.data <- GetAssay(object = object, assay.use = assay.use)
  return(VariableFeatures(object = assay.data))
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
"VariableFeatures<-.Seurat" <- function(object, assay.use = NULL, ..., value) {
  assay.use <- assay.use %||% DefaultAssay(object = object)
  assay.data <- object[[assay.use]]
  VariableFeatures(object = assay.data) <- value
  object[[assay.use]] <- assay.data
  return(object)
}

#' @describeIn WhichCells Get cells from an Assay object
#' @export
#' @method WhichCells Assay
#'
WhichCells.Assay <- function(
  object,
  cells.use,
  subset.name = NULL,
  low.threshold = -Inf,
  high.threshold = Inf,
  accept.value = NULL,
  ...
) {
  cells.use <- cells.use %||% colnames(x = object)
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
    data.use <- t(x = data.use[subset.name, cells.use, drop = FALSE])
    if (!is.null(x = accept.value)) {
      if (!all(accept.value %in% unique(x = data.use[, 1]))) {
        bad.vals <- accept.value[!(accept.value %in% unique(x = data.use[, 1]))]
        stop("Identity: ", bad.vals, " not found.")
      }
      pass.inds <- which(x = apply(data.use, MARGIN = 1, function(x) x %in% accept.value))
    } else {
      pass.inds <- which(x = (data.use > low.threshold) & (data.use < high.threshold))
    }
    cells.use <- rownames(x = data.use)[pass.inds]
  }
  return(cells.use)
}

#' @param ident.keep Create a cell subset based on the provided identity classes
#' @param ident.remove Subtract out cells from these identity classes (used for
#' filtration)
#' @param max.cells.per.ident Can be used to downsample the data to a certain
#' max per cell ident. Default is INF.
#' @param random.seed Random seed for downsampling
#' @param assay.use Which assay to filter on
#' @param ... Extra parameters passed to \code{FetchData}
#'
#' @seealso \code{\link{FetchData}}
#'
#' @describeIn WhichCells Get cells from a Seurat object
#' @export
#' @method WhichCells Seurat
#'
WhichCells.Seurat <- function(
  object,
  cells.use = NULL,
  subset.name = NULL,
  low.threshold = -Inf,
  high.threshold = Inf,
  accept.value = NULL,
  ident.keep = NULL,
  ident.remove = NULL,
  max.cells.per.ident = Inf,
  random.seed = 1,
  assay.use = NULL,
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
  cells.use <- cells.use %||% colnames(x = object)
  assay.use <- assay.use %||% DefaultAssay(object = object)
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
    cells.in.ident <- Idents(object = object)[cells.use]
    cells.in.ident <- names(x = cells.in.ident[cells.in.ident == id])
    cells.in.ident <- cells.in.ident[!is.na(x = cells.in.ident)]
    if (length(x = cells.in.ident) > max.cells.per.ident) {
      cells.in.ident <- sample(x = cells.in.ident, size = max.cells.per.ident)
    }
    cells.to.use <- c(cells.to.use, cells.in.ident)
  }
  cells.use <- cells.to.use
  if (!is.null(x = subset.name)) {
    subset.name <- as.character(subset.name)
    data.use <- FetchData(
      object = object,
      vars.fetch = subset.name,
      cells.use = cells.use,
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
    cells.use <- rownames(x = data.use)[pass.inds]
  }
  return(cells.use)
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
  utils:::findMatches(pattern, colnames(x = x[]))
}

#' @export
#'
".DollarNames.SeuratCommand" <- function(x, pattern = ''){
  params <- slot(object = x, name = "params")
  utils:::findMatches(pattern, names(x = params))
}

#' @export
#'
"$.JackStrawData" <- function(x, i, ...) {
  return(slot(object = x, name = i))
}

#' @export
#'
"$.Seurat" <- function(x, i, ...) {
  return(x[i])
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
  x[i] <- value
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
  return(GetAssayData(object = x)[i, j, ...])
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

#' @export
#'
"[.Seurat" <- function(x, i, j, ..., drop = FALSE) {
  if (missing(x = i)) {
    i <- if (missing(x = j)) {
      colnames(x = slot(object = x, name = 'meta.data'))
    } else {
      rownames(x = x)
    }
  }
  meta.count <- sum(i %in% colnames(x = slot(object = x, name = 'meta.data')))
  feat.count <- sum(
    i %in% rownames(x = x),
    vapply(X = i, FUN = is.numeric, FUN.VALUE = logical(length = 1L))
  )
  if (feat.count == 0 || meta.count > feat.count) {
    data.return <- slot(object = x, name = 'meta.data')[, i, drop = FALSE, ...]
    if (drop) {
      data.return <- unlist(x = data.return, use.names = FALSE)
      names(x = data.return) <- rep.int(x = colnames(x = x), times = length(x = i))
    }
  } else {
    if (missing(x = j)) {
      j <- colnames(x = x)
    }
    data.return <- GetAssayData(object = x)[i, j, drop = drop]
  }
  return(data.return)
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
"[[.Seurat" <- function(x, i, ...) {
  slot.use <- unlist(x = lapply(
    X = c('assays', 'reductions', 'graphs', 'neighbors', 'commands', 'workflows'),
    FUN = function(s) {
      if (i %in% names(x = slot(object = x, name = s))) {
        return(s)
      }
      return(NULL)
    }
  ))
  if (is.null(x = slot.use)) {
    stop("Cannot find '", i, "' in this Seurat object")
  }
  return(slot(object = x, name = slot.use)[[i]])
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
length.DimReduc <- function(x) {
  return(ncol(x = Embeddings(object = x)))
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
  for(i in 1:length(objects)) {
    assays.to.merge <- c(assays.to.merge, FilterObjects(object = objects[[i]], classes.keep = "Assay"))
  }
  assays.to.merge <- names(which(x = table(... = assays.to.merge) == length(x = objects)))
  combined.assays <- list()
  for(assay in assays.to.merge) {
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
  # get rid of nUMI and nFeature_*
  combined.meta.data <- data.frame(row.names = colnames(combined.assays[[1]]))
  new.idents <- c()
  for(object in objects) {
    old.meta.data <- object[]
    old.meta.data$nUMI <- NULL
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
  merged.object <- new(
    Class = 'Seurat',
    assays = combined.assays,
    meta.data = combined.meta.data,
    active.assay = assays.to.merge[1],
    active.ident = new.idents,
    project.name = project,
    version = packageVersion(pkg = 'Seurat')
  )
  merged.object['nUMI'] <- Matrix::colSums(x = merged.object)
  for(assay in assays.to.merge) {
    merged.object[paste('nFeature', assay, sep = '_')] <-
      Matrix::colSums(x = GetAssayData(
        object = merged.object,
        assay = assay, slot = "counts") > 0)
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

#' Subset a Seurat object
#'
#' @param x Seurat object to be subsetted
#' @param subset Logical expression indicating features/variables to keep
#' @param ... Arguments passed to other methods
#'
#' @return A subsetted Seurat object
#'
#' @rdname subset.Seurat
#' @aliases subset
#' @seealso \code{\link{base::subset}}
#'
#' @export
#' @method subset Seurat
#'
#' @examples
#' subset(x = pbmc, subset = MS4A1 > 7)
#'
subset.Seurat <- function(x, subset, ...) {
  objects.use <- FilterObjects(object = x)
  object.keys <- sapply(X = objects.use, FUN = function(i) {return(Key(x[[i]]))})
  key.pattern <- paste0('^', object.keys, collapse = '|')
  expr <- substitute(expr = subset)
  expr.char <- as.character(x = expr)
  expr.char <- unlist(x = lapply(X = expr.char, FUN = strsplit, split = ' '))
  vars.use <- which(
    x = expr.char %in% rownames(x = x) | expr.char %in% colnames(x = x[]) | grepl(pattern = key.pattern, x = expr.char, perl = TRUE)
  )
  data.subset <- FetchData(object = x, vars.fetch = expr.char[vars.use])
  data.subset <- subset.data.frame(x = data.subset, subset = eval(expr = expr))
  return(SubsetData(object = x, cells.use = rownames(x = data.subset)))
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# S4 methods
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

setMethod(
  f = '[<-',
  signature = c('x' = 'Seurat'),
  definition = function(x, i, ..., value) {
    meta.data <- x[]
    cell.names <- rownames(x = meta.data)
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
    return(x)
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
      stop("'i' must be a character")
    }
    slot.use <- switch(
      EXPR = as.character(x = class(x = value)),
      'Assay' = 'assays',
      'Graph' = 'graphs',
      'DimReduc' = {
        if (is.null(x = DefaultAssay(object = value))) {
          stop("Cannot add a DimReduc without an assay associated with it")
        }
        'reductions'
      },
      'SeuratCommand' = 'commands',
      'SeuratWorkflow' = 'workflows',
      stop("Unknown object type: ", class(x = value))
    )
    if (class(x = value) != 'SeuratCommand' && !all(colnames(x = value) == colnames(x = x))) {
      stop("All cells in the object being added must match the cells in this object")
    }
    if (i %in% names(x = x) && class(x = value) != class(x = x[[i]])) {
      stop("This object already contains ", i, " as a ", class(x = x[[i]]), "; duplicate names are not allowed", call. = FALSE)
    }
    if (class(x = value) %in% c('Assay', 'DimReduc') && length(x = Key(object = value)) == 0) {
      Key(object = value) <- paste0(tolower(x = i), '_')
    }
    slot(object = x, name = slot.use)[[i]] <- value
    return(x)
  }
)

setMethod(
  f = 'colMeans',
  signature = c('x' = 'Assay'),
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
    return(colSums(
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
    return(rowMeans(
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
    return(rowSums(
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
      top.ten <- VariableFeatures(object = object)[1:10]
      cat(
        "Top 10 variable features:\n",
        paste(strwrap(x = paste(top.ten, collapse = ', ')), collapse = '\n ')
      )
    }
  }
)

setMethod(
  f = 'show',
  signature = 'DimReduc',
  definition = function(object) {
    projected.data <- slot(object = object, name = 'feature.loadings.projected')
    projected <- !(all(is.na(x = projected.data)) && unique(x = dim(x = projected.data)) %in% c(0, 1))
    cat(
      "A dimensional reduction object with key", Key(object = object), '\n',
      'Number of dimensions:', length(x = object), '\n',
      'Projected dimensional reduction calculated:', projected, '\n',
    #   'Jackstraw run:', !is.null(x = object@jackstraw), '\n'
      'Jackstraw run:', as.logical(x = JS(object = object)), '\n'
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
      "dimensions."
    )
  }
)

setMethod(
  f = "show",
  signature = "Seurat",
  definition = function(object) {
    assays <- FilterObjects(object = object, classes.keep = 'Assay')
    num.features <- sum(vapply(
      X = assays,
      FUN = function(x) {
        return(nrow(x = object[[x]]))
      },
      FUN.VALUE = integer(length = 1L)
    ))
    num.assays <- length(x = assays)
    cat("An object of class", class(x = object))
    cat(
      '\n',
      num.features,
      'features across',
      ncol(x = object),
      'samples within',
      num.assays,
      ifelse(test = num.assays == 1, yes = 'assay', no = 'assays')
    )
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
    invisible(x = NULL)
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
    cat("An old seurat object\n", nrow(x = object@data), 'genes across', ncol(x = object@data), 'samples')
  }
)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Internal
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

# ...
#
# @param data.use ...
# @param dim.use ...
# @param num.use ...
# @param do.balanced ...
#
# @return ...
#
Top <- function(
  data.use,
  dim.use,
  num.use,
  do.balanced
) {
  top <- if (do.balanced) {
    num.use <- round(x = num.use / 2)
    data.use <- data.use[order(data.use), , drop = FALSE]
    positive <- head(x = rownames(x = data.use), n = num.use)
    negative <- rev(x = tail(x = rownames(x = data.use), n = num.use))
    list(positive = positive, negative = negative)
  } else {
    data.use <- data.use[rev(x = order(abs(x = data.use))), , drop = FALSE]
    top <- head(x = rownames(x = data.use), n = num.use)
    top[order(data.use[top, ])]
  }
  return(top)
}
