#' @include generics.R
#' @importFrom Rcpp evalCpp
#' @importFrom Matrix colSums rowSums colMeans rowMeans
#' @importFrom methods setClass setOldClass setClassUnion slot
#' slot<- setMethod new signature slotNames is
#' @importClassesFrom Matrix dgCMatrix
#' @useDynLib Seurat
#'
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Class definitions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

setOldClass(Classes = 'package_version')
setClassUnion(name = 'AnyMatrix', members = c("matrix", "dgCMatrix"))
setClassUnion(name = 'OptionalCharacter', members = c('NULL', 'character'))

#' The AnchorSet Class
#'
#' The AnchorSet class is an intermediate data storage class that stores the anchors and other
#' related information needed for performing downstream analyses - namely data integration
#' (\code{\link{IntegrateData}}) and data transfer (\code{\link{TransferData}}).
#'
#' @slot object.list List of objects used to create anchors
#' @slot reference.cells List of cell names in the reference dataset - needed when performing data
#' transfer.
#' @slot reference.objects Position of reference object/s in object.list
#' @slot query.cells List of cell names in the query dataset - needed when performing data transfer
#' @slot anchors The anchor matrix. This contains the cell indices of both anchor pair cells, the
#' anchor score, and the index of the original dataset in the object.list for cell1 and cell2 of
#' the anchor.
#' @slot offsets The offsets used to enable cell look up in downstream functions
#' @slot anchor.features The features used when performing anchor finding.
#' @slot command Store log of parameters that were used
#'
#' @name AnchorSet-class
#' @rdname AnchorSet-class
#' @exportClass AnchorSet
#'
AnchorSet <- setClass(
  Class = "AnchorSet",
  slots = list(
    object.list = "list",
    reference.cells = "vector",
    reference.objects = "vector",
    query.cells = "vector",
    anchors = "ANY",
    offsets = "ANY",
    anchor.features = "ANY",
    command = "ANY"
  )
)

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
#' @slot assay.orig Original assay that this assay is based off of. Used to track
#' assay provenence
#' @slot var.features Vector of features exhibiting high variance across single cells
#' @slot meta.features Feature-level metadata
#' @slot misc Utility slot for storing additional data associated with the assay
#'
#' @name Assay-class
#' @rdname Assay-class
#' @exportClass Assay
#'
Assay <- setClass(
  Class = 'Assay',
  slots = c(
    counts = 'AnyMatrix',
    data = 'AnyMatrix',
    scale.data = 'matrix',
    key = 'character',
    assay.orig = 'OptionalCharacter',
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
#' @name JackStrawData-class
#' @rdname JackStrawData-class
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
#' The DimReduc object stores a dimensionality reduction taken out in Seurat; each DimReduc
#' consists of a cell embeddings matrix, a feature loadings matrix, and a projected feature
#' loadings matrix.
#'
#' @slot cell.embeddings Cell embeddings matrix (required)
#' @slot feature.loadings Feature loadings matrix (optional)
#' @slot feature.loadings.projected Projected feature loadings matrix (optional)
#' @slot assay.used Name of assay used to generate \code{DimReduc} object
#' @slot global Is this \code{DimReduc} global/persistent? If so, it will not be
#' removed when removing its associated assay
#' @slot stdev A vector of standard deviations
#' @slot key Key for the \code{DimReduc}, must be alphanumerics followed by an underscore
#' @slot jackstraw A \code{\link{JackStrawData-class}} object associated with
#' this \code{DimReduc}
#' @slot misc Utility slot for storing additional data associated with the
#' \code{DimReduc} (e.g. the total variance of the PCA)
#'
#' @name DimReduc-class
#' @rdname DimReduc-class
#' @exportClass DimReduc
#'
DimReduc <- setClass(
  Class = 'DimReduc',
  slots = c(
    cell.embeddings = 'matrix',
    feature.loadings = 'matrix',
    feature.loadings.projected = 'matrix',
    assay.used = 'character',
    global = 'logical',
    stdev = 'numeric',
    key = 'character',
    jackstraw = 'JackStrawData',
    misc = 'list'
  )
)

#' The Graph Class
#'
#' The Graph class inherits from dgCMatrix. We do this to enable future expandability of graphs.
#'
#' @slot assay.used Optional name of assay used to generate \code{Graph} object
#'
#' @name Graph-class
#' @rdname Graph-class
#' @exportClass Graph
#'
#' @seealso \code{\link[Matrix]{dgCMatrix-class}}
#'
Graph <- setClass(
  Class = 'Graph',
  contains = "dgCMatrix",
  slots = list(
    assay.used = 'OptionalCharacter'
  )
)

#' The IntegrationData Class
#'
#' The IntegrationData object is an intermediate storage container used internally throughout the
#' integration procedure to hold bits of data that are useful downstream.
#'
#' @slot neighbors List of neighborhood information for cells (outputs of \code{RANN::nn2})
#' @slot weights Anchor weight matrix
#' @slot integration.matrix Integration matrix
#' @slot anchors Anchor matrix
#' @slot offsets The offsets used to enable cell look up in downstream functions
#' @slot objects.ncell Number of cells in each object in the object.list
#' @slot sample.tree Sample tree used for ordering multi-dataset integration
#'
#' @name IntegrationData-class
#' @rdname IntegrationData-class
#' @exportClass IntegrationData
#'
IntegrationData <- setClass(
  Class = "IntegrationData",
  slots = list(
    neighbors = "ANY",
    weights = "ANY",
    integration.matrix = "ANY",
    anchors = "ANY",
    offsets = "ANY",
    objects.ncell = "ANY",
    sample.tree = "ANY"
  )
)

#' The SeuratCommand Class
#'
#' The SeuratCommand is used for logging commands that are run on a SeuratObject. It stores parameters and timestamps
#'
#' @slot name Command name
#' @slot time.stamp Timestamp of when command was tun
#' @slot assay.used Optional name of assay used to generate \code{SeuratCommand} object
#' @slot call.string String of the command call
#' @slot params List of parameters used in the command call
#'
#' @name SeuratCommand-class
#' @rdname SeuratCommand-class
#' @exportClass SeuratCommand
#'
SeuratCommand <- setClass(
  Class = 'SeuratCommand',
  slots = c(
    name = 'character',
    time.stamp = 'POSIXct',
    assay.used = 'OptionalCharacter',
    call.string = 'character',
    params = 'ANY'
  )
)

#' The Seurat Class
#'
#' The Seurat object is a representation of single-cell expression data for R; each Seurat
#' object revolves around a set of cells and consists of one or more \code{\link{Assay-class}}
#' objects, or individual representations of expression data (eg. RNA-seq, ATAC-seq, etc).
#' These assays can be reduced from their high-dimensional state to a lower-dimension state
#' and stored as \code{\link{DimReduc-class}} objects. Seurat objects also store additional
#' meta data, both at the cell and feature level (contained within individual assays). The
#' object was designed to be as self-contained as possible, and easily extendible to new methods.
#'
#' @slot assays A list of assays for this project
#' @slot meta.data Contains meta-information about each cell, starting with number of genes detected (nGene)
#' and the original identity class (orig.ident); more information is added using \code{AddMetaData}
#' @slot active.assay Name of the active, or default, assay; settable using \code{\link{DefaultAssay}}
#' @slot active.ident The active cluster identity for this Seurat object; settable using \code{\link{Idents}}
#' @slot graphs A list of \code{\link{Graph-class}} objects
#' @slot neighbors ...
#' @slot reductions A list of dimmensional reduction objects for this object
#' @slot project.name Name of the project
#' @slot misc A list of miscellaneous information
#' @slot version Version of Seurat this object was built under
#' @slot commands A list of logged commands run on this \code{Seurat} object
#' @slot tools A list of miscellaneous data generated by other tools, should be filled by developers only using \code{\link{Tool}<-}
#'
#' @name Seurat-class
#' @rdname Seurat-class
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
    misc = 'list',
    version = 'package_version',
    commands = 'list',
    tools = 'list'
  )
)

#' The Seurat Class
#'
#' The Seurat object is the center of each single cell analysis. It stores all information
#' associated with the dataset, including data, annotations, analyes, etc. All that is needed
#' to construct a Seurat object is an expression matrix (rows are genes, columns are cells), which
#' should be log-scale
#'
#' Each Seurat object has a number of slots which store information. Key slots to access
#' are listed below.
#'
#' @slot raw.data The raw project data
#' @slot data The normalized expression matrix (log-scale)
#' @slot scale.data scaled (default is z-scoring each gene) expression matrix; used for dimmensional reduction and heatmap visualization
#' @slot var.genes Vector of genes exhibiting high variance across single cells
#' @slot is.expr Expression threshold to determine if a gene is expressed (0 by default)
#' @slot ident THe 'identity class' for each cell
#' @slot meta.data Contains meta-information about each cell, starting with number of genes detected (nGene)
#' and the original identity class (orig.ident); more information is added using \code{AddMetaData}
#' @slot project.name Name of the project (for record keeping)
#' @slot dr List of stored dimmensional reductions; named by technique
#' @slot assay List of additional assays for multimodal analysis; named by technique
#' @slot hvg.info The output of the mean/variability analysis for all genes
#' @slot imputed Matrix of imputed gene scores
#' @slot cell.names Names of all single cells (column names of the expression matrix)
#' @slot cluster.tree List where the first element is a phylo object containing the phylogenetic tree relating different identity classes
#' @slot snn Spare matrix object representation of the SNN graph
#' @slot calc.params Named list to store all calculation-related parameter choices
#' @slot kmeans Stores output of gene-based clustering from \code{DoKMeans}
#' @slot spatial Stores internal data and calculations for spatial mapping of single cells
#' @slot misc Miscellaneous spot to store any data alongisde the object (for example, gene lists)
#' @slot version Version of package used in object creation
#'
#' @name seurat-class
#' @rdname oldseurat-class
#' @aliases seurat-class
#'
seurat <- setClass(
  Class = "seurat",
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

#' Pull Assays or assay names
#'
#' Lists the names of \code{\link{Assay}} objects present in
#' a Seurat object. If slot is provided, pulls specified Assay object.
#'
#' @param object A Seurat object
#' @param slot Name of Assay to return
#'
#' @return If \code{slot} is \code{NULL}, the names of all \code{Assay} objects
#' in this Seurat object. Otherwise, the \code{Assay} object specified
#'
#' @export
#'
#' @examples
#' Assays(object = pbmc_small)
#'
Assays <- function(object, slot = NULL) {
  assays <- FilterObjects(object = object, classes.keep = 'Assay')
  if (is.null(x = slot)) {
    return(assays)
  }
  if (!slot %in% assays) {
    warning(
      "Cannot find an assay of name ",
      slot,
      " in this Seurat object",
      call. = FALSE,
      immediate. = TRUE
    )
  }
  return(slot(object = object, name = 'assays')[[slot]])
}

#' Get cell names grouped by identity class
#'
#' @param object A Seurat object
#' @param idents A vector of identity class levels to limit resulting list to;
#' defaults to all identity class levels
#' @param cells A vector of cells to grouping to
#'
#' @return A named list where names are identity classes and values are vectors
#' of cells beloning to that class
#'
#' @export
#'
#' @examples
#' CellsByIdentities(object = pbmc_small)
#'
CellsByIdentities <- function(object, idents = NULL, cells = NULL) {
  cells <- cells %||% colnames(x = object)
  cells <- intersect(x = cells, y = colnames(x = object))
  if (length(x = cells) == 0) {
    stop("Cannot find cells provided")
  }
  idents <- idents %||% levels(x = object)
  idents <- intersect(x = idents, y = levels(x = object))
  if (length(x = idents) == 0) {
    stop("None of the provided identity class levels were found", call. = FALSE)
  }
  cells.idents <- sapply(
    X = idents,
    FUN = function(i) {
      return(cells[as.vector(x = Idents(object = object)[cells]) == i])
    },
    simplify = FALSE,
    USE.NAMES = TRUE
  )
  if (any(is.na(x = Idents(object = object)[cells]))) {
    cells.idents["NA"] <- names(x = which(x = is.na(x = Idents(object = object)[cells])))
  }
  return(cells.idents)
}

#' Create an Assay object
#'
#' Create an Assay object from a feature (e.g. gene) expression matrix. The
#' expected format of the input matrix is features x cells.
#'
#' Non-unique cell or feature names are not allowed. Please make unique before
#' calling this function.
#'
#' @param counts Unnormalized data such as raw counts or TPMs
#' @param data Prenormalized data; if provided, do not pass \code{counts}
#' @param min.cells Include features detected in at least this many cells. Will
#' subset the counts matrix as well. To reintroduce excluded features, create a
#' new object with a lower cutoff.
#' @param min.features Include cells where at least this many features are
#' detected.
#'
#' @importFrom methods as
#' @importFrom Matrix colSums rowSums
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
      warning(
        "Non-unique features (rownames) present in the input matrix, making unique",
        call. = FALSE,
        immediate. = TRUE
      )
      rownames(x = counts) <- make.unique(names = rownames(x = counts))
    }
    if (anyDuplicated(colnames(x = counts))) {
      warning(
        "Non-unique cell names (colnames) present in the input matrix, making unique",
        call. = FALSE,
        immediate. = TRUE
      )
      colnames(x = counts) <- make.unique(names = colnames(x = counts))
    }
    if (is.null(x = colnames(x = counts))) {
      stop("No cell names (colnames) names present in the input matrix")
    }
    if (any(rownames(x = counts) == '')) {
      stop("Feature names of counts matrix cannot be empty", call. = FALSE)
    }
    if (nrow(x = counts) > 0 && is.null(x = rownames(x = counts))) {
      stop("No feature names (rownames) names present in the input matrix")
    }
    if (!inherits(x = counts, what = 'dgCMatrix')) {
      counts <- as(object = as.matrix(x = counts), Class = 'dgCMatrix')
    }
    # Filter based on min.features
    if (min.features > 0) {
      nfeatures <- Matrix::colSums(x = counts > 0)
      counts <- counts[, which(x = nfeatures >= min.features)]
    }
    # filter genes on the number of cells expressing
    if (min.cells > 0) {
      num.cells <- Matrix::rowSums(x = counts > 0)
      counts <- counts[which(x = num.cells >= min.cells), ]
    }
    data <- counts
  } else if (!missing(x = data)) {
    # check that dimnames of input data are unique
    if (anyDuplicated(rownames(x = data))) {
      warning(
        "Non-unique features (rownames) present in the input matrix, making unique",
        call. = FALSE,
        immediate. = TRUE
      )
      rownames(x = data) <- make.unique(names = rownames(x = data))
    }
    if (anyDuplicated(colnames(x = data))) {
      warning(
        "Non-unique cell names (colnames) present in the input matrix, making unique",
        call. = FALSE,
        immediate. = TRUE
      )
      colnames(x = data) <- make.unique(names = colnames(x = data))
    }
    if (is.null(x = colnames(x = data))) {
      stop("No cell names (colnames) names present in the input matrix")
    }
    if (any(rownames(x = data) == '')) {
      stop("Feature names of data matrix cannot be empty", call. = FALSE)
    }
    if (nrow(x = data) > 0 && is.null(x = rownames(x = data))) {
      stop("No feature names (rownames) names present in the input matrix")
    }
    if (min.cells != 0 | min.features != 0) {
      warning(
        "No filtering performed if passing to data rather than counts",
        call. = FALSE,
        immediate. = TRUE
      )
    }
    counts <- new(Class = 'matrix')
  }
  # Ensure row- and column-names are vectors, not arrays
  if (!is.vector(x = rownames(x = counts))) {
    rownames(x = counts) <- as.vector(x = rownames(x = counts))
  }
  if (!is.vector(x = colnames(x = counts))) {
    colnames(x = counts) <- as.vector(x = colnames(x = counts))
  }
  if (!is.vector(x = rownames(x = data))) {
    rownames(x = data) <- as.vector(x = rownames(x = data))
  }
  if (!is.vector(x = colnames(x = data))) {
    colnames(x = data) <- as.vector(x = colnames(x = data))
  }
  if (any(grepl(pattern = '_', x = rownames(x = counts))) || any(grepl(pattern = '_', x = rownames(x = data)))) {
    warning(
      "Feature names cannot have underscores ('_'), replacing with dashes ('-')",
      call. = FALSE,
      immediate. = TRUE
    )
    rownames(x = counts) <- gsub(
      pattern = '_',
      replacement = '-',
      x = rownames(x = counts)
    )
    rownames(x = data) <- gsub(
      pattern = '_',
      replacement = '-',
      x = rownames(x = data)
    )
  }
  if (any(grepl(pattern = '|', x = rownames(x = counts), fixed = TRUE)) || any(grepl(pattern = '|', x = rownames(x = data), fixed = TRUE))) {
    warning(
      "Feature names cannot have pipe characters ('|'), replacing with dashes ('-')",
      call. = FALSE,
      immediate. = TRUE
    )
    rownames(x = counts) <- gsub(
      pattern = '|',
      replacement = '-',
      x = rownames(x = counts),
      fixed = TRUE
    )
    rownames(x = data) <- gsub(
      pattern = '|',
      replacement = '-',
      x = rownames(x = data),
      fixed = TRUE
    )
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
#' @param embeddings A matrix with the cell embeddings
#' @param loadings A matrix with the feature loadings
#' @param projected A matrix with the projected feature loadings
#' @param assay Assay used to calculate this dimensional reduction
#' @param stdev Standard deviation (if applicable) for the dimensional reduction
#' @param key A character string to facilitate looking up features from a
#' specific DimReduc
#' @param global Specify this as a global reduction (useful for visualizations)
#' @param jackstraw Results from the JackStraw function
#' @param misc list for the user to store any additional information associated
#' with the dimensional reduction
#'
#' @aliases SetDimReduction
#'
#' @export
#'
#' @examples
#' data <- GetAssayData(pbmc_small[["RNA"]], slot = "scale.data")
#' pcs <- prcomp(x = data)
#' pca.dr <- CreateDimReducObject(
#'   embeddings = pcs$rotation,
#'   loadings = pcs$x,
#'   stdev = pcs$sdev,
#'   key = "PC",
#'   assay = "RNA"
#' )
#'
CreateDimReducObject <- function(
  embeddings = new(Class = 'matrix'),
  loadings = new(Class = 'matrix'),
  projected = new(Class = 'matrix'),
  assay = NULL,
  stdev = numeric(),
  key = NULL,
  global = FALSE,
  jackstraw = NULL,
  misc = list()
) {
  if (is.null(x = assay)) {
    warning(
      "No assay specified, setting assay as RNA by default.",
      call. = FALSE,
      immediate. = TRUE
    )
    assay <- "RNA"
  }
  # Try to infer key from column names
  if (is.null(x = key) && is.null(x = colnames(x = embeddings))) {
    stop("Please specify a key for the DimReduc object")
  } else if (is.null(x = key)) {
    key <- regmatches(
      x = colnames(x = embeddings),
      m = regexec(pattern = '^[[:alnum:]]+_', text = colnames(x = embeddings))
    )
    key <- unique(x = unlist(x = key, use.names = FALSE))
  }
  if (length(x = key) != 1) {
    stop("Please specify a key for the DimReduc object")
  } else if (!grepl(pattern = '^[[:alnum:]]+_$', x = key)) {
    old.key  <- key
    key <- UpdateKey(key = old.key)
    colnames(x = embeddings) <- gsub(
      x = colnames(x = embeddings),
      pattern = old.key,
      replacement = key
    )
    warning(
      "All keys should be one or more alphanumeric characters followed by an underscore '_', setting key to ",
      key,
      call. = FALSE,
      immediate. = TRUE
    )
  }
  # ensure colnames of the embeddings are the key followed by a numeric
  if (is.null(x = colnames(x = embeddings))) {
    warning(
      "No columnames present in cell embeddings, setting to '",
      key,
      "1:",
      ncol(x = embeddings),
      "'",
      call. = FALSE,
      immediate. = TRUE
    )
    colnames(x = embeddings) <- paste0(key, 1:ncol(x = embeddings))
  } else if (!all(grepl(pattern = paste0('^', key, "[[:digit:]]+$"), x = colnames(x = embeddings)))) {
    digits <- unlist(x = regmatches(
      x = colnames(x = embeddings),
      m = regexec(pattern = '[[:digit:]]+$', text = colnames(x = embeddings))
    ))
    if (length(x = digits) != ncol(x = embeddings)) {
      stop("Please ensure all column names in the embeddings matrix are the key plus a digit representing a dimension number")
    }
    colnames(x = embeddings) <- paste0(key, digits)
  }
  if (!IsMatrixEmpty(x = loadings)) {
    if (any(rownames(x = loadings) == '')) {
      stop("Feature names of loadings matrix cannot be empty", call. = FALSE)
    }
    colnames(x = loadings) <- colnames(x = embeddings)
  }
  if (!IsMatrixEmpty(x = projected)) {
    if (any(rownames(x = loadings) == '')) {
      stop("Feature names of projected loadings matrix cannot be empty", call. = FALSE)
    }
    colnames(x = projected) <- colnames(x = embeddings)
  }
  jackstraw <- jackstraw %||% new(Class = 'JackStrawData')
  dim.reduc <- new(
    Class = 'DimReduc',
    cell.embeddings = embeddings,
    feature.loadings = loadings,
    feature.loadings.projected = projected,
    assay.used = assay,
    global = global,
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
#' cell's column name. E.g. If your cells are named as BARCODE-CLUSTER-CELLTYPE, set this to "-" to
#' separate the cell name into its component parts for picking the relevant field.
#' @param meta.data Additional cell-level metadata to add to the Seurat object. Should be a data
#' frame where the rows are cell names and the columns are additional metadata fields.
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
  if (!is.null(x = meta.data)) {
    if (is.null(x = rownames(x = meta.data))) {
      stop("Row names not set in metadata. Please ensure that rownames of metadata match column names of data matrix")
    }
    if (length(x = setdiff(x = rownames(x = meta.data), y = colnames(x = counts)))) {
      warning("Some cells in meta.data not present in provided counts matrix.")
      meta.data <- meta.data[intersect(x = rownames(x = meta.data), y = colnames(x = counts)), ]
    }
    if (is.data.frame(x = meta.data)) {
      new.meta.data <- data.frame(row.names = colnames(x = counts))
      for (ii in 1:ncol(x = meta.data)) {
        new.meta.data[rownames(x = meta.data), colnames(x = meta.data)[ii]] <- meta.data[, ii, drop = FALSE]
      }
      meta.data <- new.meta.data
    }
  }
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
    X = colnames(x = assay.data),
    FUN = ExtractField,
    field = names.field,
    delim = names.delim
  )))
  if (any(is.na(x = idents))) {
    warning("Input parameters result in NA values for initial cell identities. Setting all initial idents to the project name")
  }
  # if there are more than 100 idents, set all idents to ... name
  ident.levels <- length(x = unique(x = idents))
  if (ident.levels > 100 || ident.levels == 0 || ident.levels == length(x = idents)) {
    idents <- rep.int(x = factor(x = project), times = ncol(x = assay.data))
  }
  names(x = idents) <- colnames(x = assay.data)
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
  n.calc <- CalcN(object = assay.data)
  if (!is.null(x = n.calc)) {
    names(x = n.calc) <- paste(names(x = n.calc), assay, sep = '_')
    object[[names(x = n.calc)]] <- n.calc
  }
  if (!is.null(x = meta.data)) {
    object <- AddMetaData(object = object, metadata = meta.data)
  }
  return(object)
}

#' Slim down a Seurat object
#'
#' Keep only certain aspects of the Seurat object. Can be useful in functions that utilize merge as
#' it reduces the amount of data in the merge.
#'
#' @param object Seurat object
#' @param counts Preserve the count matrices for the assays specified
#' @param data Preserve the data slot for the assays specified
#' @param scale.data Preserve the scale.data slot for the assays specified
#' @param features Only keep a subset of features, defaults to all features
#' @param assays Only keep a subset of assays specified here
#' @param dimreducs Only keep a subset of DimReducs specified here (if NULL,
#' remove all DimReducs)
#' @param graphs Only keep a subset of Graphs specified here (if NULL, remove
#' all Graphs)
#'
#' @export
#'
DietSeurat <- function(
  object,
  counts = TRUE,
  data = TRUE,
  scale.data = FALSE,
  features = NULL,
  assays = NULL,
  dimreducs = NULL,
  graphs = NULL
) {
  assays <- assays %||% FilterObjects(object = object, classes.keep = "Assay")
  assays <- assays[assays %in% FilterObjects(object = object, classes.keep = 'Assay')]
  if (length(x = assays) == 0) {
    stop("No assays provided were found in the Seurat object")
  }
  if (!DefaultAssay(object = object) %in% assays) {
    stop("The default assay is slated to be removed, please change the default assay")
  }
  if (!counts && !data) {
    stop("Either one or both of 'counts' and 'data' must be kept")
  }
  for (assay in FilterObjects(object = object, classes.keep = 'Assay')) {
    if (!(assay %in% assays)) {
      object[[assay]] <- NULL
    } else {
      features.assay <- features %||% rownames(x = object[[assay]])
      features.assay <- intersect(x = features.assay, y = rownames(x = object[[assay]]))
      if (length(x = features.assay) == 0) {
        if (assay == DefaultAssay(object = object)) {
          stop("The default assay is slated to be removed, please change the default assay")
        } else {
          warning("No features found in assay '", assay, "', removing...")
          object[[assay]] <- NULL
        }
      } else {
        if (counts) {
          if (!is.null(x = features)) {
            slot(object = object[[assay]], name = 'counts') <- slot(object = object[[assay]], name = 'counts')[features.assay, ]
          }
        } else {
          slot(object = object[[assay]], name = 'counts') <- new(Class = 'matrix')
        }
        if (data) {
          if (!is.null(x = features)) {
            slot(object = object[[assay]], name = 'data') <- slot(object = object[[assay]], name = 'data')[features.assay, ]
          }
        } else {
          stop('data = FALSE currently not supported')
          slot(object = object[[assay]], name = 'data') <- new(Class = 'matrix')
        }
        features.scaled <- features.assay[features.assay %in% rownames(x = slot(object = object[[assay]], name = 'scale.data'))]
        if (scale.data && length(x = features.scaled) > 0) {
          if (! all(rownames(x = slot(object = object[[assay]], name = 'scale.data')) %in% features.scaled)) {
            slot(object = object[[assay]], name = 'scale.data') <-  slot(object = object[[assay]], name = 'scale.data')[features.scaled, ]
          }
        } else {
          slot(object = object[[assay]], name = 'scale.data') <- new(Class = 'matrix')
        }
      }
    }
  }
  # remove unspecified DimReducs and Graphs
  all.objects <- FilterObjects(object = object, classes.keep = c('DimReduc', 'Graph'))
  objects.to.remove <- all.objects[!all.objects %in% c(dimreducs, graphs)]
  for (ob in objects.to.remove) {
    object[[ob]] <- NULL
  }
  return(object)
}

#' Access cellular data
#'
#' Retreives data (feature expression, PCA scores, metrics, etc.) for a set
#' of cells in a Seurat object
#'
#' @param object Seurat object
#' @param vars List of all variables to fetch, use keyword 'ident' to pull identity classes
#' @param cells Cells to collect data for (default is all cells)
#' @param slot Slot to pull feature data for
#'
#' @return A data frame with cells as rows and cellular data as columns
#'
#' @export
#'
#' @examples
#' pc1 <- FetchData(object = pbmc_small, vars = 'PC_1')
#' head(x = pc1)
#' head(x = FetchData(object = pbmc_small, vars = c('groups', 'ident')))
#'
FetchData <- function(object, vars, cells = NULL, slot = 'data') {
  cells <- cells %||% colnames(x = object)
  if (is.numeric(x = cells)) {
    cells <- colnames(x = object)[cells]
  }
  # Get a list of all objects to search through and their keys
  objects.use <- FilterObjects(object = object)
  object.keys <- sapply(X = objects.use, FUN = function(i) {return(Key(object[[i]]))})
  # Find all vars that are keyed
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
      data.return <- if (inherits(x = object[[x]], what = 'DimReduc')) {
        vars.use <- grep(
          pattern = paste0('^', key.use, '[[:digit:]]+$'),
          x = vars.use,
          value = TRUE
        )
        if (length(x = vars.use) > 0) {
          tryCatch(
            expr = object[[x]][[cells, vars.use, drop = FALSE]],
            error = function(...) {
              return(NULL)
            }
          )
        } else {
          NULL
        }
      } else if (inherits(x = object[[x]], what = 'Assay')) {
        vars.use <- gsub(pattern = paste0('^', key.use), replacement = '', x = vars.use)
        data.assay <- GetAssayData(
          object = object,
          slot = slot,
          assay = x
        )
        vars.use <- vars.use[vars.use %in% rownames(x = data.assay)]
        data.vars <- t(x = as.matrix(data.assay[vars.use, cells, drop = FALSE]))
        if (ncol(data.vars) > 0) {
          colnames(x = data.vars) <- paste0(key.use, vars.use)
        }
        data.vars
      }
      data.return <- as.list(x = as.data.frame(x = data.return))
      return(data.return)
    }
  )
  data.fetched <- unlist(x = data.fetched, recursive = FALSE)
  # Pull vars from object metadata
  meta.vars <- vars[vars %in% colnames(x = object[[]]) & ! vars %in% names(x = data.fetched)]
  data.fetched <- c(data.fetched, object[[meta.vars]][cells, , drop = FALSE])
  # Pull vars from the default assay
  default.vars <- vars[vars %in% rownames(x = GetAssayData(object = object, slot = slot)) & ! vars %in% names(x = data.fetched)]
  data.fetched <- c(
    data.fetched,
    tryCatch(
      expr = as.data.frame(x = t(x = as.matrix(x = GetAssayData(
        object = object,
        slot = slot
      )[default.vars, cells, drop = FALSE]))),
      error = function(...) {
        return(NULL)
      }
    )
  )
  # Pull identities
  if ('ident' %in% vars && !'ident' %in% colnames(x = object[[]])) {
    data.fetched[['ident']] <- Idents(object = object)[cells]
  }
  # Try to find ambiguous vars
  fetched <- names(x = data.fetched)
  vars.missing <- setdiff(x = vars, y = fetched)
  if (length(x = vars.missing) > 0) {
    # Search for vars in alternative assays
    vars.alt <- vector(mode = 'list', length = length(x = vars.missing))
    names(x = vars.alt) <- vars.missing
    for (assay in FilterObjects(object = object, classes.keep = 'Assay')) {
      vars.assay <- Filter(
        f = function(x) {
          features.assay <- rownames(x = GetAssayData(
            object = object,
            assay = assay,
            slot = slot
          ))
          return(x %in% features.assay)
        },
        x = vars.missing
      )
      for (var in vars.assay) {
        vars.alt[[var]] <- append(x = vars.alt[[var]], values = assay)
      }
    }
    # Vars found in multiple alternative assays are truly ambiguous, will not pull
    vars.many <- names(x = Filter(
      f = function(x) {
        return(length(x = x) > 1)
      },
      x = vars.alt
    ))
    if (length(x = vars.many) > 0) {
      warning(
        "Found the following features in more than one assay, excluding the default. We will not include these in the final dataframe: ",
        paste(vars.many, collapse = ', '),
        call. = FALSE,
        immediate. = TRUE
      )
    }
    vars.missing <- names(x = Filter(
      f = function(x) {
        return(length(x = x) != 1)
      },
      x = vars.alt
    ))
    # Pull vars found in only one alternative assay
    # Key this var to highlight that it was found in an alternate assay
    vars.alt <- Filter(
      f = function(x) {
        return(length(x = x) == 1)
      },
      x = vars.alt
    )
    for (var in names(x = vars.alt)) {
      assay <- vars.alt[[var]]
      warning(
        'Could not find ',
        var,
        ' in the default search locations, found in ',
        assay,
        ' assay instead',
        immediate. = TRUE,
        call. = FALSE
      )
      keyed.var <- paste0(Key(object = object[[assay]]), var)
      data.fetched[[keyed.var]] <- as.vector(
        x = GetAssayData(object = object, assay = assay, slot = slot)[var, cells]
      )
      vars <- sub(
        pattern = paste0('^', var, '$'),
        replacement = keyed.var,
        x = vars
      )
    }
    fetched <- names(x = data.fetched)
  }
  # Name the vars not found in a warning (or error if no vars found)
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
  # Assembled fetched vars in a dataframe
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

#' Get integation data
#'
#' @param object Seurat object
#' @param integration.name Name of integration object
#' @param slot Which slot in integration object to get
#'
#' @return Returns data from the requested slot within the integrated object
#'
#' @export
#'
GetIntegrationData <- function(object, integration.name, slot) {
  tools <- slot(object = object, name = 'tools')
  if (!(integration.name %in% names(tools))) {
    stop('Requested integration key does not exist')
  }
  int.data <- tools[[integration.name]]
  return(slot(object = int.data, name = slot))
}

#' Log a command
#'
#' Logs command run, storing the name, timestamp, and argument list. Stores in
#' the Seurat object
#'
#' @param object Name of Seurat object
#' @param return.command Return a \link{SeuratCommand} object instead
#'
#' @return If \code{return.command}, returns a SeuratCommand object. Otherwise,
#' returns the Seurat object with command stored
#'
#' @export
#'
#' @seealso \code{\link{Command}}
#'
LogSeuratCommand <- function(object, return.command = FALSE) {
  time.stamp <- Sys.time()
  #capture function name
  which.frame <- sys.nframe() - 1
  if (which.frame < 1) {
    stop("'LogSeuratCommand' cannot be called at the top level", call. = FALSE)
  }
  command.name <- as.character(x = deparse(expr = sys.calls()[[which.frame]]))
  command.name <- gsub(pattern = "\\.Seurat", replacement = "", x = command.name)
  call.string <- command.name
  command.name <- ExtractField(string = command.name, field = 1, delim = "\\(")
  #capture function arguments
  argnames <- names(x = formals(fun = sys.function(which = sys.parent(n = 1))))
  argnames <- grep(pattern = "object", x = argnames, invert = TRUE, value = TRUE)
  argnames <- grep(pattern = "anchorset", x = argnames, invert = TRUE, value = TRUE)
  argnames <- grep(pattern = "\\.\\.\\.", x = argnames, invert = TRUE, value = TRUE)
  params <- list()
  p.env <- parent.frame(n = 1)
  argnames <- intersect(x = argnames, y = ls(name = p.env))
  # fill in params list
  for (arg in argnames) {
    param_value <- get(x = arg, envir = p.env)
    if (inherits(x = param_value, what = 'Seurat')) {
      next
    }
    #TODO Institute some check of object size?
    params[[arg]] <- param_value
  }
  # check if function works on the Assay and/or the DimReduc Level
  assay <- params[["assay"]]
  reduction <- params[["reduction"]]
  # Get assay used for command
  cmd.assay <- assay %||% (reduction %iff% if (inherits(x = reduction, what = 'DimReduc')) {
    DefaultAssay(object = reduction)
  } else if (reduction %in% Reductions(object = object)) {
    DefaultAssay(object = object[[reduction]])
  })
  if (inherits(x = reduction, what = 'DimReduc')) {
    reduction <- 'DimReduc'
  }
  # rename function name to include Assay/DimReduc info
  if (length(x = assay) == 1) {
    command.name <- paste(command.name, assay, reduction, sep = '.')
  }
  command.name <- sub(pattern = "[\\.]+$", replacement = "", x = command.name, perl = TRUE)
  command.name <- sub(pattern = "\\.\\.", replacement = "\\.", x = command.name, perl = TRUE)
  # store results
  seurat.command <- new(
    Class = 'SeuratCommand',
    name = command.name,
    params = params,
    time.stamp = time.stamp,
    call.string = call.string,
    assay.used = cmd.assay
  )
  if (return.command) {
    return(seurat.command)
  }
  object[[command.name]] <- seurat.command
  return(object)
}

#' Pull DimReducs or DimReduc names
#'
#' Lists the names of \code{\link{DimReduc}} objects present in
#' a Seurat object. If slot is provided, pulls specified DimReduc object.
#'
#' @param object A Seurat object
#' @param slot Name of DimReduc
#'
#' @return If \code{slot} is \code{NULL}, the names of all \code{DimReduc} objects
#' in this Seurat object. Otherwise, the \code{DimReduc} object requested
#'
#' @export
#'
#' @examples
#' Reductions(object = pbmc_small)
#'
Reductions <- function(object, slot = NULL) {
  reductions <- FilterObjects(object = object, classes.keep = 'DimReduc')
  if (is.null(x = slot)) {
    return(reductions)
  }
  if (!slot %in% reductions) {
    warning(
      "Cannot find a DimReduc of name ",
      slot,
      " in this Seurat object",
      call. = FALSE,
      immediate. = TRUE
    )
  }
  return(slot(object = object, name = 'reductions')[[slot]])
}

#' Rename assays in a \code{Seurat} object
#'
#' @param object A \code{Seurat} object
#' @param ... Named arguments as \code{old.assay = new.assay}
#'
#' @return \code{object} with assays renamed
#'
#' @export
#' @examples
#' RenameAssays(object = pbmc_small, RNA = 'rna')
#'
RenameAssays <- function(object, ...) {
  assay.pairs <- tryCatch(
    expr = as.list(x = ...),
    error = function(e) {
      return(list(...))
    }
  )
  old.assays <- names(x = assay.pairs)
  # Handle missing assays
  missing.assays <- setdiff(x = old.assays, y = Assays(object = object))
  if (length(x = missing.assays) == length(x = old.assays)) {
    stop("None of the assays provided are present in this object", call. = FALSE)
  } else if (length(x = missing.assays)) {
    warning(
      "The following assays could not be found: ",
      paste(missing.assays, collapse = ', '),
      call. = FALSE,
      immediate. = TRUE
    )
  }
  old.assays <- setdiff(x = old.assays, missing.assays)
  assay.pairs <- assay.pairs[old.assays]
  # Check to see that all old assays are named
  if (is.null(x = names(x = assay.pairs)) || any(sapply(X = old.assays, FUN = nchar) < 1)) {
    stop("All arguments must be named with the old assay name", call. = FALSE)
  }
  # Ensure each old assay is going to one new assay
  if (!all(sapply(X = assay.pairs, FUN = length) == 1) || length(x = old.assays) != length(x = unique(x = old.assays))) {
    stop("Can only rename assays to one new name", call. = FALSE)
  }
  # Ensure each new assay is coming from one old assay
  if (length(x = assay.pairs) != length(x = unique(x = assay.pairs))) {
    stop(
      "One or more assays are set to be lost due to duplicate new assay names",
      call. = FALSE
    )
  }
  # Rename assays
  for (old in names(x = assay.pairs)) {
    new <- assay.pairs[[old]]
    # If we aren't actually renaming any
    if (old == new) {
      next
    }
    old.key <- Key(object = object[[old]])
    suppressWarnings(expr = object[[new]] <- object[[old]])
    if (old == DefaultAssay(object = object)) {
      message("Renaming default assay from ", old, " to ", new)
      DefaultAssay(object = object) <- new
    }
    Key(object = object[[new]]) <- old.key
    object[[old]] <- NULL
  }
  return(object)
}

#' Set integation data
#'
#' @param object Seurat object
#' @param integration.name Name of integration object
#' @param slot Which slot in integration object to set
#' @param new.data New data to insert
#'
#' @return Returns a \code{\link{Seurat}} object
#'
#' @export
#'
SetIntegrationData <- function(object, integration.name, slot, new.data) {
  tools <- slot(object = object, name = 'tools')
  if (!(integration.name %in% names(tools))) {
    new.integrated <- new(Class = 'IntegrationData')
    slot(object = new.integrated, name = slot) <- new.data
    tools[[integration.name]] <- new.integrated
    slot(object = object, name = 'tools') <- tools
    return(object)
  }
  int.data <- tools[[integration.name]]
  slot(object = int.data, name = slot) <- new.data
  tools[[integration.name]] <- int.data
  slot(object = object, name = 'tools') <- tools
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
SplitObject <- function(object, split.by = "ident") {
  if (split.by == 'ident') {
    groupings <- Idents(object = object)
  } else {
    groupings <- FetchData(object = object, vars = split.by)[, 1]
  }
  groupings <- unique(x = as.character(x = groupings))
  obj.list <- list()
  for (i in groupings) {
    if (split.by == "ident") {
      obj.list[[i]] <- subset(x = object, idents = i)
    }
    else {
      cells <- which(x = object[[split.by, drop = TRUE]] == i)
      cells <- colnames(x = object)[cells]
      obj.list[[i]] <- subset(x = object, cells = cells)
    }
  }
  return(obj.list)
}

#' Find features with highest scores for a given dimensional reduction technique
#'
#' Return a list of features with the strongest contribution to a set of components
#'
#' @param object DimReduc object
#' @param dim Dimension to use
#' @param nfeatures Number of features to return
#' @param projected Use the projected feature loadings
#' @param balanced Return an equal number of features with both + and - scores.
#' @param ... Extra parameters passed to \code{\link{Loadings}}
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
  balanced = FALSE,
  ...
) {
  loadings <- Loadings(object = object, projected = projected, ...)[, dim, drop = FALSE]
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
#' @param dim Dimension to use
#' @param ncells Number of cells to return
#' @param balanced Return an equal number of cells with both + and - scores.
#' @param ... Extra parameters passed to \code{\link{Embeddings}}
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
TopCells <- function(object, dim = 1, ncells = 20, balanced = FALSE, ...) {
  embeddings <- Embeddings(object = object, ...)[, dim, drop = FALSE]
  return(Top(
    data = embeddings,
    num = ncells,
    balanced = balanced
  ))
}

#' Update old Seurat object to accomodate new features
#'
#' Updates Seurat objects to new structure for storing data/calculations.
#' For Seurat v3 objects, will validate object structure ensuring all keys and feature
#' names are formed properly.
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
    if (slot(object = object, name = 'version') >= package_version(x = "2.0.0") && slot(object = object, name = 'version') < package_version(x = '3.0.0')) {
      # Run update
      message("Updating from v2.X to v3.X")
      seurat.version <- packageVersion(pkg = "Seurat")
      new.assay <- UpdateAssay(old.assay = object, assay = "RNA")
      assay.list <- list(new.assay)
      names(x = assay.list) <- "RNA"
      for (i in names(x = object@assay)) {
        assay.list[[i]] <- UpdateAssay(old.assay = object@assay[[i]], assay = i)
      }
      new.dr <- UpdateDimReduction(old.dr = object@dr, assay = "RNA")
      object <- new(
        Class = "Seurat",
        version = seurat.version,
        assays = assay.list,
        active.assay = "RNA",
        project.name = object@project.name,
        misc = object@misc %||% list(),
        active.ident = object@ident,
        reductions = new.dr,
        meta.data = object@meta.data,
        tools = list()
      )
      # Run CalcN
      for (assay in Assays(object = object)) {
        n.calc <- CalcN(object = object[[assay]])
        if (!is.null(x = n.calc)) {
          names(x = n.calc) <- paste(names(x = n.calc), assay, sep = '_')
          object[[names(x = n.calc)]] <- n.calc
        }
        to.remove <- c("nGene", "nUMI")
        for (i in to.remove) {
          if (i %in% colnames(x = object[[]])) {
            object[[i]] <- NULL
          }
        }
      }
    }
    if (package_version(x = slot(object = object, name = 'version')) >= package_version(x = "3.0.0")) {
      # Run validation
      message("Validating object structure")
      # Update object slots
      message("Updating object slots")
      object <- UpdateSlots(object = object)
      # Rename assays
      assays <- make.names(names = Assays(object = object))
      names(x = assays) <- Assays(object = object)
      object <- do.call(what = RenameAssays, args = c('object' = object, assays))
      for (obj in FilterObjects(object = object, classes.keep = c('Assay', 'DimReduc', 'Graph'))) {
        suppressWarnings(expr = object[[obj]] <- UpdateSlots(object = object[[obj]]))
      }
      for (cmd in Command(object = object)) {
        slot(object = object, name = 'commands')[[cmd]] <- UpdateSlots(
          object = Command(object = object, command = cmd)
        )
      }
      # Validate object keys
      message("Ensuring keys are in the proper strucutre")
      for (ko in FilterObjects(object = object)) {
        Key(object = object[[ko]]) <- UpdateKey(key = Key(object = object[[ko]]))
      }
      # Check feature names
      message("Ensuring feature names don't have underscores or pipes")
      for (assay.name in FilterObjects(object = object, classes.keep = 'Assay')) {
        assay <- object[[assay.name]]
        for (slot in c('counts', 'data', 'scale.data')) {
          if (!IsMatrixEmpty(x = slot(object = assay, name = slot))) {
            rownames(x = slot(object = assay, name = slot)) <- gsub(
              pattern = '_',
              replacement = '-',
              x = rownames(x = slot(object = assay, name = slot))
            )
            rownames(x = slot(object = assay, name = slot)) <- gsub(
              pattern = '|',
              replacement = '-',
              x = rownames(x = slot(object = assay, name = slot)),
              fixed = TRUE
            )
          }
        }
        VariableFeatures(object = assay) <- gsub(
          pattern = '_',
          replacement = '-',
          x = VariableFeatures(object = assay)
        )
        VariableFeatures(object = assay) <- gsub(
          pattern = '|',
          replacement = '-',
          x = VariableFeatures(object = assay),
          fixed = TRUE
        )
        rownames(x = slot(object = assay, name = "meta.features")) <-  gsub(
          pattern = '_',
          replacement = '-',
          x = rownames(x = assay[[]])
        )
        rownames(x = slot(object = assay, name = "meta.features")) <-  gsub(
          pattern = '|',
          replacement = '-',
          x = rownames(x = assay[[]]),
          fixed = TRUE
        )
        object[[assay.name]] <- assay
      }
      for (reduc.name in FilterObjects(object = object, classes.keep = 'DimReduc')) {
        reduc <- object[[reduc.name]]
        for (slot in c('feature.loadings', 'feature.loadings.projected')) {
          if (!IsMatrixEmpty(x = slot(object = reduc, name = slot))) {
            rownames(x = slot(object = reduc, name = slot)) <- gsub(
              pattern = '_',
              replacement = '-',
              x = rownames(x = slot(object = reduc, name = slot))
            )
            rownames(x = slot(object = reduc, name = slot)) <- gsub(
              pattern = '_',
              replacement = '-',
              x = rownames(x = slot(object = reduc, name = slot)),
              fixed = TRUE
            )
          }
        }
        object[[reduc.name]] <- reduc
      }
    }
    if (package_version(x = slot(object = object, name = 'version')) <= package_version(x = '3.1.1')) {
      # Update Assays, DimReducs, and Graphs
      for (x in names(x = object)) {
        message("Updating slots in ", x)
        xobj <- object[[x]]
        xobj <- UpdateSlots(object = xobj)
        if (inherits(x = xobj, what = 'DimReduc')) {
          if (any(sapply(X = c('tsne', 'umap'), FUN = grepl, x = tolower(x = x)))) {
            message("Setting ", x, " DimReduc to global")
            slot(object = xobj, name = 'global') <- TRUE
          }
        } else if (inherits(x = xobj, what = 'Graph')) {
          graph.assay <- unlist(x = strsplit(x = x, split = '_'))[1]
          if (graph.assay %in% Assays(object = object)) {
            message("Setting default assay of ", x, " to ", graph.assay)
            DefaultAssay(object = xobj) <- graph.assay
          }
        }
        object[[x]] <- xobj
      }
      # Update SeuratCommands
      for (cmd in Command(object = object)) {
        cobj <- Command(object = object, command = cmd)
        cobj <- UpdateSlots(object = cobj)
        cmd.assay <- unlist(x = strsplit(x = cmd, split = '\\.'))
        cmd.assay <- cmd.assay[length(x = cmd.assay)]
        cmd.assay <- if (cmd.assay %in% Assays(object = object)) {
          cmd.assay
        } else if (cmd.assay %in% Reductions(object = object)) {
          DefaultAssay(object = object[[cmd.assay]])
        } else {
          NULL
        }
        if (is.null(x = cmd.assay)) {
          message("No assay information could be found for ", cmd)
        } else {
          message("Setting assay used for ", cmd, " to ", cmd.assay)
        }
        slot(object = cobj, name = 'assay.used') <- cmd.assay
        object[[cmd]] <- cobj
      }
      # Update object version
      slot(object = object, name = 'version') <- packageVersion(pkg = 'Seurat')
    }
    message("Object representation is consistent with the most current Seurat version")
    return(object)
  }
  stop(
    "We are unable to convert Seurat objects less than version 2.X to version 3.X\n",
    'Please use devtools::install_version to install Seurat v2.3.4 and update your object to a 2.X object',
    call. = FALSE
  )
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Methods for Seurat-defined generics
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @rdname AddMetaData
#' @export
#' @method AddMetaData Assay
#'
AddMetaData.Assay <- function(object, metadata, col.name = NULL) {
  return(.AddMetaData(object = object, metadata = metadata, col.name = col.name))
}

#' @rdname AddMetaData
#' @export
#' @method AddMetaData Seurat
#'
AddMetaData.Seurat <- function(object, metadata, col.name = NULL) {
  return(.AddMetaData(object = object, metadata = metadata, col.name = col.name))
}

#' @param assay Assay to convert
#' @param reduction Name of DimReduc to set to main reducedDim in cds
#'
#' @rdname as.CellDataSet
#' @export
#' @method as.CellDataSet Seurat
#'
as.CellDataSet.Seurat <- function(x, assay = NULL, reduction = NULL, ...) {
  CheckDots(...)
  if (!PackageCheck('monocle', error = FALSE)) {
    stop("Please install monocle from Bioconductor before converting to a CellDataSet object")
  } else if (packageVersion(pkg = 'monocle') >= package_version(x = '2.99.0')) {
    stop("Seurat can only convert to/from Monocle v2.X objects")
  }
  assay <- assay %||% DefaultAssay(object = x)
  # make variables, then run `newCellDataSet`
  # create cellData counts
  counts <- GetAssayData(object = x, assay = assay, slot = "counts")
  # metadata
  cell.metadata <- x[[]]
  feature.metadata <- x[[assay]][[]]
  if (!"gene_short_name" %in% colnames(x = feature.metadata)) {
    feature.metadata$gene_short_name <- rownames(x = feature.metadata)
  }
  pd <- new(Class = "AnnotatedDataFrame", data = cell.metadata)
  fd <- new(Class = "AnnotatedDataFrame", data = feature.metadata)
  # Now, determine the expressionFamily
  if ("monocle" %in% names(x = Misc(object = x))) {
    expressionFamily <- Misc(object = x, slot = "monocle")[["expressionFamily"]]
  } else {
    if (all(counts == floor(x = counts))) {
      expressionFamily <- VGAM::negbinomial.size()
    } else if (any(counts < 0)) {
      expressionFamily <- VGAM::uninormal()
    } else {
      expressionFamily <- VGAM::tobit()
    }
  }
  cds <- monocle::newCellDataSet(
    cellData = counts,
    phenoData = pd,
    featureData = fd,
    expressionFamily = expressionFamily
  )
  if ("monocle" %in% names(x = Misc(object = x))) {
    monocle::cellPairwiseDistances(cds = cds) <- Misc(object = x, slot = "monocle")[["cellPairwiseDistances"]]
    monocle::minSpanningTree(cds = cds) <- Misc(object = x, slot = "monocle")[["minSpanningTree"]]
    Biobase::experimentData(cds = cds) <- Misc(object = x, slot = "monocle")[["experimentData"]]
    Biobase::protocolData(cds = cds) <- Misc(object = x, slot = "monocle")[["protocolData"]]
    Biobase::classVersion(cds = cds) <- Misc(object = x, slot = "monocle")[["classVersion"]]
    # no setter methods found for following slots
    slot(object = cds, name = "lowerDetectionLimit") <- Misc(object = x, slot = "monocle")[["lowerDetectionLimit"]]
    slot(object = cds, name = "dispFitInfo") <- Misc(object = x, slot = "monocle")[["dispFitInfo"]]
    slot(object = cds, name = "auxOrderingData") <- Misc(object = x, slot = "monocle")[["auxOrderingData"]]
    slot(object = cds, name = "auxClusteringData") <- Misc(object = x, slot = "monocle")[["auxClusteringData"]]
  }
  # adding dimensionality reduction data to the CDS
  dr.slots <- c("reducedDimS", "reducedDimK", "reducedDimW", "reducedDimA")
  reduction <- reduction %||% DefaultDimReduc(object = x, assay = assay)
  if (!is.null(x = reduction)) {
    if (grepl(pattern = 'tsne', x = tolower(x = reduction))) {
      slot(object = cds, name = "dim_reduce_type") <- "tSNE"
      monocle::reducedDimA(cds = cds) <- t(x = Embeddings(object = x[[reduction]]))
    } else {
      slot(object = cds, name = "dim_reduce_type") <- reduction
      monocle::reducedDimA(cds = cds) <- Loadings(object = x[[reduction]])
      slot(object = cds, name = "reducedDimS") <- Embeddings(object = x[[reduction]])
    }
    for (ii in dr.slots) {
      if (ii %in% names(x = slot(object = x[[reduction]], name = "misc"))) {
        slot(object = cds, name = ii) <- slot(object = x[[reduction]], name = "misc")[[ii]]
      }
    }
  }
  return(cds)
}

#' @rdname as.Graph
#' @export
#' @method as.Graph Matrix
#'
#' @examples
#' # converting sparse matrix
#' mat <- Matrix::rsparsematrix(nrow = 10, ncol = 10, density = 0.1)
#' rownames(x = mat) <- paste0("feature_", 1:10)
#' colnames(x = mat) <- paste0("cell_", 1:10)
#' g <- as.Graph(x = mat)
#'
as.Graph.Matrix <- function(x, ...) {
  CheckDots(...)
  x <- as.sparse(x = x)
  if (is.null(x = rownames(x = x))) {
    stop("Please provide rownames to the matrix before converting to a Graph.")
  }
  if (is.null(x = colnames(x = x))) {
    stop("Please provide colnames to the matrix before converting to a Graph.")
  }
  return(as(object = x, Class = "Graph"))
}

#' @rdname as.Graph
#' @export
#' @method as.Graph matrix
#'
#' @examples
#' # converting dense matrix
#' mat <- matrix(data = 1:16, nrow = 4)
#' rownames(x = mat) <- paste0("feature_", 1:4)
#' colnames(x = mat) <- paste0("cell_", 1:4)
#' g <- as.Graph(x = mat)
#'
as.Graph.matrix <- function(x, ...) {
  CheckDots(...)
  return(as.Graph.Matrix(x = as(object = x, Class = 'Matrix')))
}

#' @details
#' The Seurat method for \code{as.loom} will try to automatically fill in datasets based on data presence.
#' For example, if an assay's scaled data slot isn't filled, then dimensional reduction and graph information
#' will not be filled, since those depend on scaled data. The following is a list of how datasets will be filled
#' \itemize{
#'   \item \code{counts} will be stored in \code{matrix}
#'   \item Cell names will be stored in \code{col_attrs/CellID}; feature names will be stored in \code{row_attrs/Gene}
#'   \item \code{data} will be stored in \code{layers/norm_data}
#'   \item \code{scale.data} will be stored in \code{layers/scale_data}
#'   \item Cell-level metadata will be stored in \code{col_attrs}; all periods '.' in metadata will be replaced with underscores '_'
#'   \item Clustering information from \code{Idents(object = x)} will be stored in \code{col_attrs/ClusterID} and \code{col_attrs/ClusterName}
#'   for the numeric and string representation of the factor, respectively
#'   \item Feature-level metadata will be stored in \code{Feature_attrs}; all periods '.' in metadata will be replaced with underscores '_'
#'   \item Variable features, if set, will be stored in \code{row_attrs/Selected}; features declared as variable will be stored as '1',
#'   others will be stored as '0'
#'   \item Dimensional reduction information for the assay provided will be stored in \code{col_attrs} for cell embeddings and \code{row_attrs}
#'    for feature loadings; datasets will be named as \code{name_type} where \code{name} is the name within the Seurat object
#'    and \code{type} is \code{cell_embeddings} or \code{feature_loadings}; if feature loadings have been projected for all features,
#'    then projected loadings will be stored instead and \code{type} will be \code{feature_loadings_projected}
#'   \item Nearest-neighbor graphs that start with the name of the assay will be stored in \code{col_graphs}
#'   \item Assay information will be stored as an HDF5 attribute called \code{assay} at the root level
#' }
#'
#' @inheritParams loomR::create
#' @param assay Assay to store in loom file
#'
#' @rdname as.loom
#' @export
#' @method as.loom Seurat
#'
#' @examples
#' \dontrun{
#' lfile <- as.loom(x = pbmc_small)
#' }
#'
as.loom.Seurat <- function(
  x,
  assay = NULL,
  filename = file.path(getwd(), paste0(Project(object = x), '.loom')),
  max.size = '400mb',
  chunk.dims = NULL,
  chunk.size = NULL,
  overwrite = FALSE,
  verbose = TRUE,
  ...
) {
  if (!PackageCheck('loomR', error = FALSE)) {
    stop("Please install loomR from GitHub before converting to a loom object")
  }
  CheckDots(..., fxns = 'loomR::create')
  # Set the default assay to make life easy
  assay <- assay %||% DefaultAssay(object = x)
  DefaultAssay(object = x) <- assay
  # Pull ordering information
  cell.order <- colnames(x = x)
  feature.order <- rownames(x = x)
  # Get cell- and feature-level metadata
  meta.data <- x[[]][cell.order, ]
  colnames(x = meta.data) <- gsub(
    pattern = '\\.',
    replacement = '_',
    x = colnames(x = meta.data)
  )
  meta.data$ClusterID <- as.integer(x = Idents(object = x)[rownames(x = meta.data)])
  meta.data$ClusterName <- as.character(x = Idents(object = x)[rownames(x = meta.data)])
  meta.feature <- x[[assay]][[]][feature.order, ]
  colnames(x = meta.feature) <- gsub(
    pattern = '\\.',
    replacement = '_',
    x = colnames(x = meta.feature)
  )
  if (length(x = VariableFeatures(object = x)) > 0) {
    meta.feature[VariableFeatures(object = x), 'Selected'] <- 1
    meta.feature[is.na(x = meta.feature$Selected), 'Selected'] <- 0
  }
  if (IsMatrixEmpty(x = GetAssayData(object = x, slot = 'counts'))) {
    data <- GetAssayData(object = x, slot = 'data')
    layers <- NULL
  } else {
    data <- GetAssayData(object = x, slot = 'counts') # Raw counts matrix
    layers = list('norm_data' = GetAssayData(object = x, slot = 'data')) # Add data slot as norm_data
  }
  # Make the initial loom object
  lfile <- loomR::create(
    filename = filename,
    data = data[feature.order, cell.order],
    feature.attrs = as.list(x = meta.feature), # Feature-level metadata
    cell.attrs = as.list(x = meta.data), # Cell-level metadata
    layers = layers,
    transpose = TRUE,
    calc.count = FALSE,
    max.size = max.size,
    chunk.size = chunk.size,
    chunk.dims = chunk.dims,
    overwrite = overwrite,
    verbose = verbose,
    ...
  )
  # Add scale.data
  if (!IsMatrixEmpty(x = GetAssayData(object = x, slot = 'scale.data'))) {
    if (verbose) {
      message("Adding scaled data matrix to /layers/scale_data")
    }
    lfile$add.layer(
      layers = list(
        'scale_data' = as.matrix(
          x = t(
            x = as.data.frame(
              x = GetAssayData(object = x, slot = 'scale.data')
            )[feature.order, cell.order]
          )
        )
      ),
      verbose = verbose
    )
    dim.reducs <- FilterObjects(object = x, classes.keep = 'DimReduc')
    dim.reducs <- Filter(
      f = function(d) {
        return(DefaultAssay(object = x[[d]]) == assay)
      },
      x = dim.reducs
    )
    # Add dimensional reduction information
    for (dr in dim.reducs) {
      if (verbose) {
        message("Adding dimensional reduction information for ", dr)
      }
      embeddings <- Embeddings(object = x, reduction = dr)[cell.order, ]
      embeddings <- list(embeddings)
      names(x = embeddings) <- paste0(dr, '_cell_embeddings')
      if (verbose) {
        message("Adding cell embedding information for ", dr)
      }
      lfile$add.col.attribute(attributes = embeddings)
      loadings <- Loadings(
        object = x,
        reduction = dr,
        projected = Projected(object = x[[dr]])
      )
      # Add feature loading information
      if (!IsMatrixEmpty(x = loadings)) {
        if (verbose) {
          message("Adding feature loading information for ", dr)
        }
        loadings <- as.matrix(x = as.data.frame(x = loadings)[feature.order, ])
        loadings <- list(loadings)
        names(x = loadings) <- paste0(dr, '_feature_loadings')
        if (Projected(object = x[[dr]])) {
          names(x = loadings) <- paste0(names(x = loadings), '_projected')
        }
        lfile$add.row.attribute(attributes = loadings)
      } else if (verbose) {
        message("No feature loading information for ", dr)
      }
    }
    # Add graph information
    graphs <- FilterObjects(object = x, classes.keep = 'Graph')
    graphs <- grep(pattern = paste0('^', assay), x = graphs, value = TRUE)
    for (gr in graphs) {
      if (verbose) {
        message("Adding graph ", gr)
      }
      lfile$add.graph.matrix(mat = x[[gr]], name = gr, MARGIN = 2)
    }
  } else if (verbose) {
    message("No scaled data present, not adding scaled data, dimensional reduction information, or neighbor graphs")
  }
  # Store assay
  hdf5r::h5attr(x = lfile, which = 'assay') <- assay
  return(lfile)
}

#' @param slot Slot to store expression data as
#'
#' @importFrom utils packageVersion
#'
#' @rdname as.Seurat
#' @export
#' @method as.Seurat CellDataSet
#'
as.Seurat.CellDataSet <- function(
  x,
  slot = 'counts',
  assay = 'RNA',
  verbose = TRUE,
  ...
) {
  CheckDots(...)
  if (!PackageCheck('monocle', error = FALSE)) {
    stop("Please install monocle from Bioconductor before converting to a CellDataSet object")
  } else if (packageVersion(pkg = 'monocle') >= package_version(x = '2.99.0')) {
    stop("Seurat can only convert to/from Monocle v2.X objects")
  }
  slot <- match.arg(arg = slot, choices = c('counts', 'data'))
  if (verbose) {
    message("Pulling expression data")
  }
  expr <- Biobase::exprs(object = x)
  if (IsMatrixEmpty(x = expr)) {
    stop("No data provided in this CellDataSet object", call. = FALSE)
  }
  meta.data <- as.data.frame(x = Biobase::pData(object = x))
  # if cell names are NULL, fill with cell_X
  if (is.null(x = colnames(x = expr))) {
    warning(
      "The column names of the 'counts' and 'data' matrices are NULL. Setting cell names to cell_columnidx (e.g 'cell_1').",
      call. = FALSE,
      immediate. = TRUE
    )
    rownames(x = meta.data) <- colnames(x = expr) <- paste0("cell_", 1:ncol(x = expr))
  }
  # Creating the object
  if (verbose) {
    message("Building Seurat object")
  }
  if (slot == 'data') {
    assays <- list(CreateAssayObject(data = expr))
    names(x = assays) <- assay
    Key(object = assays[[assay]]) <- suppressWarnings(expr = UpdateKey(key = assay))
    object <- new(
      Class = 'Seurat',
      assays = assays,
      meta.data = meta.data,
      version = packageVersion(pkg = 'Seurat'),
      project.name = 'SeuratProject'
    )
    DefaultAssay(object = object) <- assay
  } else {
    object <- CreateSeuratObject(
      counts = expr,
      meta.data = meta.data,
      assay = assay
    )
  }
  # feature metadata
  if (verbose) {
    message("Adding feature-level metadata")
  }
  feature.metadata <- Biobase::fData(object = x)
  object[[assay]][[names(x = feature.metadata)]] <- feature.metadata
  # mean/dispersion values
  disp.table <- tryCatch(
    expr = suppressWarnings(expr = monocle::dispersionTable(cds = x)),
    error = function(...) {
      return(NULL)
    }
  )
  if (!is.null(x = disp.table)) {
    if (verbose) {
      message("Adding dispersion information")
    }
    rownames(x = disp.table) <- disp.table[, 1]
    disp.table[, 1] <- NULL
    colnames(x = disp.table) <- paste0('monocle_', colnames(x = disp.table))
    object[[assay]][[names(x = disp.table)]] <- disp.table
  } else if (verbose) {
    message("No dispersion information in CellDataSet object")
  }
  # variable features
  if ("use_for_ordering" %in% colnames(x = feature.metadata)) {
    if (verbose) {
      message("Setting variable features")
    }
    VariableFeatures(object = object, assay = assay) <- rownames(x = feature.metadata)[which(x = feature.metadata[, "use_for_ordering"])]
  } else if (verbose) {
    message("No variable features present")
  }
  # add dim reduction
  dr.name <- slot(object = x, name = "dim_reduce_type")
  if (length(x = dr.name) > 0) {
    if (verbose) {
      message("Adding ", dr.name, " dimensional reduction")
    }
    reduced.A <- t(x = slot(object = x, name = 'reducedDimA'))
    reduced.S <- t(x = slot(object = x, name = 'reducedDimS'))
    if (IsMatrixEmpty(x = reduced.S)) {
      embeddings <- reduced.A
      loadings <- new(Class = 'matrix')
    } else {
      embeddings <- reduced.S
      loadings <- t(x = reduced.A)
    }
    rownames(x = embeddings) <- colnames(x = object)
    misc.dr <- list(
      reducedDimS = slot(object = x, name = "reducedDimS"),
      reducedDimK = slot(object = x, name = "reducedDimK"),
      reducedDimW = slot(object = x, name = "reducedDimW"),
      reducedDimA = slot(object = x, name = "reducedDimA")
    )
    dr <- suppressWarnings(expr = CreateDimReducObject(
      embeddings = embeddings,
      loadings = loadings,
      assay = assay,
      key = UpdateKey(key = tolower(x = dr.name)),
      misc = misc.dr
    ))
    object[[dr.name]] <- dr
  } else if (verbose) {
    message("No dimensional reduction information found")
  }
  monocle.specific.info <- list(
    expressionFamily = slot(object = x, name = "expressionFamily"),
    lowerDetectionLimit = slot(object = x, name = "lowerDetectionLimit"),
    dispFitInfo = slot(object = x, name = "dispFitInfo"),
    cellPairwiseDistances = slot(object = x, name = "cellPairwiseDistances"),
    minSpanningTree = slot(object = x, name = "minSpanningTree"),
    auxOrderingData = slot(object = x, name = "auxOrderingData"),
    auxClusteringData = slot(object = x, name = "auxClusteringData"),
    experimentData = slot(object = x, name = "experimentData"),
    protocolData = slot(object = x, name = "protocolData"),
    classVersion = slot(object = x, name = ".__classVersion__")
  )
  Misc(object = object, slot = "monocle")  <- monocle.specific.info
  return(object)
}

#' @details
#' The \code{loom} method for \code{as.Seurat} will try to automatically fill in a Seurat object based on data presence.
#' For example, if no normalized data is present, then scaled data, dimensional reduction informan, and neighbor graphs
#' will not be pulled as these depend on normalized data. The following is a list of how the Seurat object will be constructed
#' \itemize{
#'   \item If no assay information is provided, will default to an assay name in a root-level HDF5 attribute called \code{assay};
#'   if no attribute is present, will default to "RNA"
#'   \item Cell-level metadata will consist of all one-dimensional datasets in \code{col_attrs} \strong{except} datasets named "ClusterID", "ClusterName",
#'   and whatever is passed to \code{cells}
#'   \item Identity classes will be set if either \code{col_attrs/ClusterID} or \code{col_attrs/ClusterName} are present; if both are present, then
#'   the values in \code{col_attrs/ClusterID} will set the order (numeric value of a factor) for values in \code{col_attrs/ClusterName}
#'   (charater value of a factor)
#'   \item Feature-level metadata will consist of all one-dimensional datasets in \code{row_attrs} \strong{except} datasets named "Selected" and whatever
#'   is passed to \code{features}; any feature-level metadata named "variance_standardized", "variance_expected", or "dispersion_scaled" will have
#'   underscores "_" replaced with a period "."
#'   \item Variable features will be set if \code{row_attrs/Selected} exists and it is a numeric type
#'   \item If a dataset is passed to \code{normalized}, stored as a sparse matrix in \code{data};
#'   if no dataset provided, \code{scaled} will be set to \code{NULL}
#'   \item If a dataset is passed to \code{scaled}, stored as a dense matrix in \code{scale.data}; all rows entirely consisting of \code{NA}s
#'   will be removed
#'   \item If a dataset is passed to \code{scaled}, dimensional reduction information will assembled from cell embedding information
#'   stored in \code{col_attrs}; cell embeddings will be pulled from two-dimensional datasets ending with "_cell_embeddings"; priority will
#'   be given to cell embeddings that have the name of \code{assay} in their name; feature loadings will be added from two-dimensional
#'   datasets in \code{row_attrs} that start with the name of the dimensional reduction and end with either "feature_loadings" or
#'   "feature_loadings_projected" (priority given to the latter)
#'   \item If a dataset is passed to \code{scaled}, neighbor graphs will be pulled from \code{col_graphs}, provided the name starts
#'   with the value of \code{assay}
#' }
#'
#' @param cells The name of the dataset within \code{col_attrs} containing cell names
#' @param features The name of the dataset within \code{row_attrs} containing feature names
#' @param normalized The name of the dataset within \code{layers} containing the
#' normalized expression matrix; pass \code{/matrix} (with preceeding forward slash) to store
#' \code{/matrix} as normalized data
#' @param scaled The name of the dataset within \code{layers} containing the scaled expression matrix
#' @param verbose Display progress updates
#'
#' @importFrom Matrix sparseMatrix
#'
#' @rdname as.Seurat
#' @export
#' @method as.Seurat loom
#'
#' @examples
#' \dontrun{
#' lfile <- as.loom(x = pbmc_small)
#' pbmc <- as.Seurat(x = lfile)
#' }
#'
as.Seurat.loom <- function(
  x,
  cells = 'CellID',
  features = 'Gene',
  normalized = NULL,
  scaled = NULL,
  assay = NULL,
  verbose = TRUE,
  ...
) {
  CheckDots(...)
  # Shouldn't be necessary
  if (!PackageCheck('loomR', error = FALSE)) {
    stop("Please install loomR")
  }
  # Check prerequisite datasets
  if (!x[['col_attrs']]$exists(name = cells)) {
    stop("Cannot find provided cell name attribute in the loom file")
  }
  if (!x[['row_attrs']]$exists(name = features)) {
    stop("Cannot find provided feature name attribute in the loom file")
  }
  assay <- assay %||% hdf5r::h5attributes(x = x)$assay %||% 'RNA'
  # Read in the counts matrix
  if (verbose) {
    message(
      "Pulling ",
      ifelse(
        test = !is.null(x = normalized) && normalized == '/matrix',
        yes = 'normalized data',
        no = 'counts'
      )
      ," matrix"
    )
  }
  counts <- x$get.sparse(
    dataset = 'matrix',
    feature.names = features,
    cell.names = cells,
    verbose = verbose
  )
  if (!is.null(x = normalized) && normalized == '/matrix') {
    assays <- list(CreateAssayObject(data = counts))
    names(x = assays) <- assay
    object <- new(
      Class = 'Seurat',
      assays = assays,
      meta.data = data.frame(row.names = colnames(x = assays[[assay]])),
      version = packageVersion(pkg = 'Seurat'),
      project.name = 'SeuratProject'
    )
    DefaultAssay(object = object) <- assay
  } else {
    object <- CreateSeuratObject(
      counts = counts,
      assay = assay
    )
  }
  # Read in normalized and scaled data
  if (!is.null(x = normalized) && normalized != '/matrix') {
    normalized <- basename(path = normalized)
    if (!x[['layers']]$exists(name = normalized)) {
      warning(
        "Cannot find provided normalized data in the loom file",
        call. = FALSE,
        immediate. = TRUE
      )
      scaled <- NULL
    } else {
      if (verbose) {
        message("Adding normalized data")
      }
      norm.data <- x$get.sparse(
        dataset = paste0('layers/', normalized),
        feature.names = features,
        cell.names = cells
      )
      object <- SetAssayData(object = object, slot = 'data', new.data = norm.data)
    }
  } else if (is.null(x = normalized) || normalized != '/matrix') {
    if (verbose) {
      message("No normalized data provided, not adding scaled data")
    }
    scaled <- NULL
  }
  if (!is.null(x = scaled)) {
    scaled <- basename(path = scaled)
    if (!x[['layers']]$exists(name = scaled)) {
      warning(
        "Cannot find provided scaled data in the loom file",
        call. = FALSE,
        immediate. = TRUE
      )
      scaled <- NULL
    } else {
      if (verbose) {
        message("Adding scaled data")
      }
      scale.data <- t(x = x[['layers']][[scaled]][, ])
      rownames(x = scale.data) <- x[['row_attrs']][[features]][]
      colnames(x = scale.data) <- x[['col_attrs']][[cells]][]
      row.drop <- apply(
        X = scale.data,
        MARGIN = 1,
        FUN = function(row) {
          return(all(is.na(x = row)))
        }
      )
      scale.data <- scale.data[!row.drop, , drop = FALSE]
      object <- SetAssayData(
        object = object,
        slot = 'scale.data',
        new.data = scale.data
      )
    }
  } else if (verbose) {
    message("No scaled data provided")
  }
  # Read in cell-level metadata
  meta.data <- hdf5r::list.datasets(
    object = x,
    path = 'col_attrs',
    full.names = FALSE,
    recursive = FALSE
  )
  meta.data <- meta.data[-which(x = meta.data %in% c(cells, 'ClusterID', 'ClusterName'))]
  meta.data <- Filter(
    f = function(m) {
      return(length(x = x[['col_attrs']][[m]]$dims) == 1)
    },
    x = meta.data
  )
  if (length(x = meta.data) > 0) {
    meta.data <- sapply(
      X = meta.data,
      FUN = function(m) {
        return(x[['col_attrs']][[m]][])
      },
      simplify = FALSE,
      USE.NAMES = TRUE
    )
    meta.data <- as.data.frame(x = meta.data)
    rownames(x = meta.data) <- make.unique(names = x[['col_attrs']][[cells]][])
    colnames(x = meta.data) <- gsub(
      pattern = 'orig_ident',
      replacement = 'orig.ident',
      x = colnames(x = meta.data)
    )
    object[[colnames(x = meta.data)]] <- meta.data
  }
  # Set clustering information
  idents <- if (x[['col_attrs']]$exists(name = 'ClusterID')) {
    if (length(x = x[['col_attrs/ClusterID']]$dims) == 1) {
      x[['col_attrs/ClusterID']][]
    } else {
      NULL
    }
  } else {
    NULL
  }
  if (x[['col_attrs']]$exists(name = 'ClusterName')) {
    if (length(x = x[['col_attrs/ClusterName']]$dims) == 1) {
      ident.order <- idents
      idents <- x[['col_attrs/ClusterName']][]
    } else {
      ident.order <- NULL
    }
  } else {
    ident.order <- NULL
  }
  if (!is.null(x = idents)) {
    if (verbose) {
      message("Setting cluster IDs")
    }
    names(x = idents) <- x[['col_attrs']][[cells]][]
    levels <- if (is.null(x = ident.order)) {
      idents
    } else {
      idents[order(ident.order)]
    }
    levels <- unique(x = levels)
    idents <- factor(x = idents, levels = levels)
    Idents(object = object) <- idents
  } else if (verbose) {
    message("No clustering information present")
  }
  # Read in feature-level metadata
  meta.features <- hdf5r::list.datasets(
    object = x,
    path = 'row_attrs',
    full.names = FALSE,
    recursive = FALSE
  )
  meta.features <- meta.features[-which(x = meta.features %in% c(features, 'Selected'))]
  meta.features <- Filter(
    f = function(m) {
      return(length(x = x[['row_attrs']][[m]]$dims) == 1)
    },
    x = meta.features
  )
  if (length(x = meta.features) > 0) {
    meta.features <- sapply(
      X = meta.features,
      FUN = function(m) {
        return(x[['row_attrs']][[m]][])
      },
      simplify = FALSE,
      USE.NAMES = TRUE
    )
    meta.features <- as.data.frame(x = meta.features)
    rownames(x = meta.features) <- make.unique(names = x[['row_attrs']][[features]][])
    colnames(x = meta.features) <- gsub(
      pattern = 'variance_standardized',
      replacement = 'variance.standardized',
      x = colnames(x = meta.features)
    )
    colnames(x = meta.features) <- gsub(
      pattern = 'variance_expected',
      replacement = 'variance.expected',
      x = colnames(x = meta.features)
    )
    colnames(x = meta.features) <- gsub(
      pattern = 'dispersion_scaled',
      replacement = 'dispersion.scaled',
      x = colnames(x = meta.features)
    )
    object[[assay]][[colnames(x = meta.features)]] <- meta.features
  }
  # Look for variable features
  if (x[['row_attrs']]$exists(name = 'Selected')) {
    if (inherits(x = x[['row_attrs/Selected']]$get_type(), what = c('H5T_FLOAT', 'H5T_INTEGER'))) {
      var.features <- which(x = x[['row_attrs/Selected']][] == 1)
      VariableFeatures(object = object) <- x[['row_attrs']][[features]][var.features]
    } else if (verbose) {
      message("'Selected' must be a dataset of floats or integers, with '1' signifiying variable")
    }
  }
  # If scaled, look for dimensional reduction information
  if (!is.null(x = scaled)) {
    reductions <- hdf5r::list.datasets(
      object = x,
      path = 'col_attrs',
      full.names = FALSE,
      recursive = FALSE
    )
    reductions <- grep(
      pattern = '_cell_embeddings$',
      x = reductions,
      value = TRUE
    )
    reductions <- Filter(
      f = function(r) {
        return(length(x = x[['col_attrs']][[r]]$dims) == 2)
      },
      x = reductions
    )
    reduc.names <- sapply(
      X = strsplit(x = reductions, split = '_'),
      FUN = '[',
      1
    )
    reductions <- sapply(
      X = reduc.names,
      FUN = function(r) {
        reducs <- grep(pattern = paste0('^', r), x = reductions, value = TRUE)
        if (sum(grepl(pattern = assay, x = reducs)) == 1) {
          return(grep(pattern = assay, x = reducs, value = TRUE))
        }
        return(reducs[which.min(x = nchar(x = reducs))])
      },
      USE.NAMES = FALSE
    )
    all.loadings <- grep(
      pattern = '_feature_loadings[_projected]',
      x = names(x = x[['row_attrs']]),
      value = TRUE,
      perl = TRUE
    )
    for (reduc in reductions) {
      dim.name <- gsub(pattern = '_cell_embeddings', replacement = '', x = reduc)
      if (verbose) {
        message("Adding ", dim.name, " dimensional reduction information")
      }
      key <- switch(
        EXPR = dim.name,
        'pca' = 'PC',
        'tsne' = 'tSNE',
        toupper(x = dim.name)
      )
      key <- paste0(key, '_')
      embeddings <- t(x = x[['col_attrs']][[reduc]][, ])
      rownames(x = embeddings) <- x[['col_attrs']][[cells]][]
      dr <- CreateDimReducObject(
        embeddings = embeddings,
        assay = assay,
        key = key
      )
      loadings <- grep(pattern = dim.name, x = all.loadings, value = TRUE)
      if (length(x = loadings) == 1) {
        if (verbose) {
          message("Pulling feature loadings for ", dim.name)
        }
        projected <- grepl(pattern = '_projected$', x = loadings)
        loadings <- t(x = x[['row_attrs']][[loadings]][, ])
        rownames(x = loadings) <- if (projected) {
          x[['row_attrs']][[features]][]
        } else {
          rownames(x = GetAssayData(object = object, slot = 'scale.data'))
        }
        Loadings(object = dr, projected = projected) <- loadings
      } else if (verbose) {
        message("No loadings present for ", dim.name)
      }
      object[[dim.name]] <- dr
    }
  } else if (verbose) {
    message("No scaled data, not searching for dimensional reduction information")
  }
  # If scaled, look for graphs
  if (!is.null(x = scaled)) {
    for (gname in names(x = x[['col_graphs']])) {
      if (!grepl(pattern = paste0('^', assay), x = gname)) {
        next
      }
      if (verbose) {
        message("Loading graph ", gname)
      }
      graph <- sparseMatrix(
        i = x[['col_graphs']][[gname]][['a']][] + 1,
        j = x[['col_graphs']][[gname]][['b']][],
        x = x[['col_graphs']][[gname]][['w']][]
      )
      rownames(x = graph) <- colnames(x = graph) <- x[['col_attrs']][[cells]][]
      object[[gname]] <- as.Graph(x = graph)
    }
  } else if (verbose) {
    message("No scaled data, not searching for nearest neighbor graphs")
  }
  return(object)
}

#' @param counts name of the SingleCellExperiment assay to store as \code{counts};
#' set to \code{NULL} if only normalized data are present
#' @param data name of the SingleCellExperiment assay to slot as \code{data}.
#' Set to NULL if only counts are present
#' @param assay Name to store expression matrices as
#' @param project Project name for new Seurat object
#'
#' @rdname as.Seurat
#' @export
#' @method as.Seurat SingleCellExperiment
#'
as.Seurat.SingleCellExperiment <- function(
  x,
  counts = 'counts',
  data = 'logcounts',
  assay = 'RNA',
  project = 'SingleCellExperiment',
  ...
) {
  CheckDots(...)
  if (!PackageCheck('SingleCellExperiment', error = FALSE)) {
    stop(
      "Please install SingleCellExperiment from Bioconductor before converting to a SingeCellExperiment object",
      call. = FALSE
    )
  }
  meta.data <- as.data.frame(x = SummarizedExperiment::colData(x = x))
  # Pull expression matrices
  mats <- list(counts = counts, data = data)
  mats <- Filter(f = Negate(f = is.null), x = mats)
  if (length(x = mats) == 0) {
    stop("Cannot pass 'NULL' to both 'counts' and 'data'")
  }
  for (m in 1:length(x = mats)) {
    # if (is.null(x = mats[[m]])) next
    mats[[m]] <- tryCatch(
      expr = SummarizedExperiment::assay(x = x, i = mats[[m]]),
      error = function(e) {
        stop("No data in provided assay - ", mats[[m]], call. = FALSE)
      }
    )
    # if cell names are NULL, fill with cell_X
    if (is.null(x = colnames(x = mats[[m]]))) {
      warning(
        "The column names of the",
        names(x = mats)[m],
        " matrix is NULL. Setting cell names to cell_columnidx (e.g 'cell_1').",
        call. = FALSE,
        immediate. = TRUE
      )
      cell.names <- paste0("cell_", 1:ncol(x = mats[[m]]))
      colnames(x = mats[[m]]) <- cell.names
      rownames(x = meta.data) <- cell.names
    }
  }
  assays <- if (is.null(x = mats$counts)) {
    list(CreateAssayObject(data = mats$data))
  } else if (is.null(x = mats$data)) {
    list(CreateAssayObject(counts = mats$counts))
  } else {
    a <- CreateAssayObject(counts = mats$counts)
    a <- SetAssayData(object = a, slot = 'data', new.data = mats$data)
    list(a)
  }
  names(x = assays) <- assay
  Key(object = assays[[assay]]) <- paste0(tolower(x = assay), '_')
  # Create the Seurat object
  object <- new(
    Class = 'Seurat',
    assays = assays,
    meta.data = meta.data,
    version = packageVersion(pkg = 'Seurat'),
    project.name = project
  )
  DefaultAssay(object = object) <- assay
  Idents(object = object) <- project
  # Get DimReduc information
  if (length(x = SingleCellExperiment::reducedDimNames(x = x)) > 0) {
    for (dr in SingleCellExperiment::reducedDimNames(x = x)) {
      embeddings <- SingleCellExperiment::reducedDim(x = x, type = dr)
      if (is.null(x = rownames(x = embeddings))) {
        rownames(x = embeddings)  <- cell.names
      }
      key <- gsub(
        pattern = "[[:digit:]]",
        replacement = "_",
        x = colnames(x = SingleCellExperiment::reducedDim(x = x, type = dr))[1]
      )
      if (length(x = key) == 0) {
        key <- paste0(dr, "_")
      }
      colnames(x = embeddings) <- paste0(key, 1:ncol(x = embeddings))
      object[[dr]] <- CreateDimReducObject(
        embeddings = embeddings,
        key = key,
        assay = DefaultAssay(object = object)
      )
    }
  }
  return(object)
}

#' @param assay Assay to convert
#'
#' @rdname as.SingleCellExperiment
#' @export
#' @method as.SingleCellExperiment Seurat
#'
as.SingleCellExperiment.Seurat <- function(x, assay = NULL, ...) {
  CheckDots(...)
  if (!PackageCheck('SingleCellExperiment', error = FALSE)) {
    stop("Please install SingleCellExperiment from Bioconductor before converting to a SingeCellExperiment object")
  }
  assay <- assay %||% DefaultAssay(object = x)
  assays = list(
    counts = GetAssayData(object = x, assay = assay, slot = "counts"),
    logcounts = GetAssayData(object = x, assay = assay, slot = "data")
  )
  assays <- assays[sapply(X = assays, FUN = nrow) != 0]
  sce <- SingleCellExperiment::SingleCellExperiment(assays = assays)
  metadata <- x[[]]
  metadata$ident <- Idents(object = x)
  SummarizedExperiment::colData(sce) <- S4Vectors::DataFrame(metadata)
  SummarizedExperiment::rowData(sce) <- S4Vectors::DataFrame(x[[assay]][[]])
  for (dr in FilterObjects(object = x, classes.keep = "DimReduc")) {
    SingleCellExperiment::reducedDim(sce, toupper(x = dr)) <- Embeddings(object = x[[dr]])
  }
  return(sce)
}

#' @importFrom methods as
#' @importClassesFrom Matrix dgCMatrix
#'
#' @rdname as.sparse
#' @export
#' @method as.sparse data.frame
#'
as.sparse.data.frame <- function(x, ...) {
  CheckDots(...)
  return(as(object = as.matrix(x = x), Class = 'dgCMatrix'))
}

#' @importFrom methods is
#' @importFrom Matrix sparseMatrix
#'
#' @rdname as.sparse
#' @export
#' @method as.sparse H5Group
#'
as.sparse.H5Group <- function(x, ...) {
  CheckDots(...)
  for (i in c('data', 'indices', 'indptr')) {
    if (!x$exists(name = i) || !is(object = x[[i]], class2 = 'H5D')) {
      stop("Invalid H5Group specification for a sparse matrix, missing dataset ", i)
    }
  }
  if ('h5sparse_shape' %in% hdf5r::h5attr_names(x = x)) {
    return(sparseMatrix(
      i = x[['indices']][] + 1,
      p = x[['indptr']][],
      x = x[['data']][],
      dims = rev(x = hdf5r::h5attr(x = x, which = 'h5sparse_shape'))
    ))
  }
  return(sparseMatrix(
    i = x[['indices']][] + 1,
    p = x[['indptr']][],
    x = x[['data']][]
  ))
}

#' @importFrom methods as
#' @importClassesFrom Matrix dgCMatrix
#'
#' @rdname as.sparse
#' @export
#' @method as.sparse Matrix
#'
as.sparse.Matrix <- function(x, ...) {
  CheckDots(...)
  return(as(object = x, Class = 'dgCMatrix'))
}

#' @rdname as.sparse
#' @export
#' @method as.sparse matrix
#'
as.sparse.matrix <- function(x, ...) {
  return(as.sparse.Matrix(x = x, ...))
}

#' @rdname Cells
#' @export
#'
Cells.default <- function(x) {
  return(colnames(x = x))
}

#' @rdname Cells
#' @export
#' @method Cells DimReduc
#'
Cells.DimReduc <- function(x) {
  return(rownames(x = x))
}

#' @param command Name of the command to pull, pass \code{NULL} to get the names of all commands run
#' @param value Name of the parameter to pull the value for
#'
#' @rdname Command
#' @export
#' @method Command Seurat
#'
Command.Seurat <- function(object, command = NULL, value = NULL, ...) {
  CheckDots(...)
  commands <- slot(object = object, name = "commands")
  if (is.null(x = command)) {
    return(names(x = commands))
  }
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

#' @rdname DefaultAssay
#' @export
#' @method DefaultAssay Assay
#'
DefaultAssay.Assay <- function(object, ...) {
  object <- UpdateSlots(object = object)
  return(slot(object = object, name = 'assay.orig'))
}

#' @rdname DefaultAssay
#' @export
#' @method DefaultAssay DimReduc
#'
DefaultAssay.DimReduc <- function(object, ...) {
  CheckDots(...)
  return(slot(object = object, name = 'assay.used'))
}

#' @rdname DefaultAssay
#' @export
#' @method DefaultAssay Graph
#'
DefaultAssay.Graph <- function(object, ...) {
  object <- UpdateSlots(object = object)
  return(slot(object = object, name = 'assay.used'))
}

#' @rdname DefaultAssay
#' @export
#' @method DefaultAssay Seurat
#'
#' @examples
#' # Get current default assay
#' DefaultAssay(object = pbmc_small)
#'
DefaultAssay.Seurat <- function(object, ...) {
  CheckDots(...)
  return(slot(object = object, name = 'active.assay'))
}

#' @rdname DefaultAssay
#' @export
#' @method DefaultAssay SeuratCommand
#'
DefaultAssay.SeuratCommand <- function(object, ...) {
  object <- UpdateSlots(object = object)
  return(slot(object = object, name = 'assay.used'))
}

#' @export
#' @method DefaultAssay<- Assay
#'
"DefaultAssay<-.Assay" <- function(object, ..., value) {
  object <- UpdateSlots(object = object)
  return(slot(object = object, name = 'assay.used'))
  object <- UpdateSlots(object = object)
  slot(object = object, name = 'assay.orig') <- value
  return(object)
}

#' @export
#' @method DefaultAssay<- DimReduc
#'
"DefaultAssay<-.DimReduc" <- function(object, ..., value) {
  CheckDots(...)
  slot(object = object, name = 'assay.used') <- value
  return(object)
}

#' @export
#' @method DefaultAssay<- Graph
#'
"DefaultAssay<-.Graph" <- function(object, ..., value) {
  object <- UpdateSlots(object = object)
  slot(object = object, name = 'assay.used') <- value
  return(object)
}

#' @rdname DefaultAssay
#' @export
#' @method DefaultAssay<- Seurat
#'
#' @examples
#' # Create dummy new assay to demo switching default assays
#' new.assay <- pbmc_small[["RNA"]]
#' Key(object = new.assay) <- "RNA2_"
#' pbmc_small[["RNA2"]] <- new.assay
#' # switch default assay to RNA2
#' DefaultAssay(object = pbmc_small) <- "RNA2"
#' DefaultAssay(object = pbmc_small)
#'
"DefaultAssay<-.Seurat" <- function(object, ..., value) {
  CheckDots(...)
  if (!value %in% names(x = slot(object = object, name = 'assays'))) {
    stop("Cannot find assay ", value)
  }
  slot(object = object, name = 'active.assay') <- value
  return(object)
}

#' @rdname Embeddings
#' @export
#' @method Embeddings DimReduc
#'
#' @examples
#' # Get the embeddings directly from a DimReduc object
#' Embeddings(object = pbmc_small[["pca"]])[1:5, 1:5]
#'
Embeddings.DimReduc <- function(object, ...) {
  CheckDots(...)
  return(slot(object = object, name = 'cell.embeddings'))
}

#' @param reduction Name of reduction to pull cell embeddings for
#'
#' @rdname Embeddings
#' @export
#' @method Embeddings Seurat
#'
#' @examples
#' # Get the embeddings from a specific DimReduc in a Seurat object
#' Embeddings(object = pbmc_small, reduction = "pca")[1:5, 1:5]
#'
Embeddings.Seurat <- function(object, reduction = 'pca', ...) {
  return(Embeddings(object = object[[reduction]], ...))
}

#' @param assay Assay to get
#'
#' @rdname GetAssay
#' @export
#' @method GetAssay Seurat
#'
#' @examples
#' GetAssay(object = pbmc_small, assay = "RNA")
#'
GetAssay.Seurat <- function(object, assay = NULL, ...) {
  CheckDots(...)
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

#' @param slot Specific information to pull (i.e. counts, data, scale.data, ...)
#'
#' @rdname GetAssayData
#' @export
#' @method GetAssayData Assay
#'
#' @examples
#' # Get the data directly from an Assay object
#' GetAssayData(object = pbmc_small[["RNA"]], slot = "data")[1:5,1:5]
#'
GetAssayData.Assay <- function(object, slot = 'data', ...) {
  CheckDots(...)
  return(slot(object = object, name = slot))
}

#' @param assay Name of assay to pull data from
#'
#' @rdname GetAssayData
#' @export
#' @method GetAssayData Seurat
#'
#' @examples
#' # Get the data from a specific Assay in a Seurat object
#' GetAssayData(object = pbmc_small, assay = "RNA", slot = "data")[1:5,1:5]
#'
GetAssayData.Seurat <- function(object, slot = 'data', assay = NULL, ...) {
  CheckDots(...)
  assay <- assay %||% DefaultAssay(object = object)
  return(GetAssayData(
    object = GetAssay(object = object, assay = assay),
    slot = slot
  ))
}

#' @param selection.method Which method to pull; choose one from \code{c('sctransform', 'sct')}
#' or \code{c('mean.var.plot', 'dispersion', 'mvp', 'disp')}
#' @param status Add variable status to the resulting data.frame
#'
#' @rdname HVFInfo
#' @export
#' @method HVFInfo Assay
#'
#' @examples
#' # Get the HVF info directly from an Assay object
#' HVFInfo(object = pbmc_small[["RNA"]], selection.method = 'vst')[1:5, ]
#'
HVFInfo.Assay <- function(object, selection.method, status = FALSE, ...) {
  CheckDots(...)
  disp.methods <- c('mean.var.plot', 'dispersion', 'disp')
  if (tolower(x = selection.method) %in% disp.methods) {
    selection.method <- 'mvp'
  }
  selection.method <- switch(
    EXPR = tolower(x = selection.method),
    'sctransform' = 'sct',
    selection.method
  )
  vars <- switch(
    EXPR = selection.method,
    'vst' = c('mean', 'variance', 'variance.standardized'),
    'mvp' = c('mean', 'dispersion', 'dispersion.scaled'),
    'sct' = c('gmean', 'variance', 'residual_variance'),
    stop("Unknown method: '", selection.method, "'", call. = FALSE)
  )
  tryCatch(
    expr = hvf.info <- object[[paste(selection.method, vars, sep = '.')]],
    error = function(e) {
      stop(
        "Unable to find highly variable feature information for method '",
        selection.method,
        "'",
        call. = FALSE
      )
    }
  )
  colnames(x = hvf.info) <- vars
  if (status) {
    hvf.info$variable <- object[[paste0(selection.method, '.variable')]]
  }
  return(hvf.info)
}

#' @param assay Name of assay to pull highly variable feature information for
#'
#' @importFrom tools file_path_sans_ext
#'
#' @rdname HVFInfo
#' @export
#' @method HVFInfo Seurat
#'
#' @examples
#' # Get the HVF info from a specific Assay in a Seurat object
#' HVFInfo(object = pbmc_small, assay = "RNA")[1:5, ]
#'
HVFInfo.Seurat <- function(
  object,
  selection.method = NULL,
  assay = NULL,
  status = FALSE,
  ...
) {
  CheckDots(...)
  assay <- assay %||% DefaultAssay(object = object)
  if (is.null(x = selection.method)) {
    cmds <- apply(
      X = expand.grid(
        c('FindVariableFeatures', 'SCTransform'),
        FilterObjects(object = object, classes.keep = 'Assay')
      ),
      MARGIN = 1,
      FUN = paste,
      collapse = '.'
    )
    find.command <- Command(object = object)[Command(object = object) %in% cmds]
    if (length(x = find.command) < 1) {
      stop(
        "Please run either 'FindVariableFeatures' or 'SCTransform'",
        call. = FALSE
      )
    }
    find.command <- find.command[length(x = find.command)]
    test.command <- paste(file_path_sans_ext(x = find.command), assay, sep = '.')
    find.command <- ifelse(
      test = test.command %in% Command(object = object),
      yes = test.command,
      no = find.command
    )
    selection.method <- switch(
      EXPR = file_path_sans_ext(x = find.command),
      'FindVariableFeatures' = Command(
        object = object,
        command = find.command,
        value = 'selection.method'
      ),
      'SCTransform' = 'sct',
      stop("Unknown command for finding variable features: '", find.command, "'", call. = FALSE)
    )
  }
  return(HVFInfo(
    object = GetAssay(object = object, assay = assay),
    selection.method = selection.method,
    status = status
  ))
}

#' @rdname Idents
#' @export
#' @method Idents Seurat
#'
Idents.Seurat <- function(object, ...) {
  CheckDots(...)
  return(slot(object = object, name = 'active.ident'))
}

#' @param cells Set cell identities for specific cells
#' @param drop Drop unused levels
#'
#' @rdname Idents
#' @export
#' @method Idents<- Seurat
#'
"Idents<-.Seurat" <- function(object, cells = NULL, drop = FALSE, ..., value) {
  CheckDots(...)
  cells <- cells %||% colnames(x = object)
  if (is.numeric(x = cells)) {
    cells <- colnames(x = object)[cells]
  }
  cells <- intersect(x = cells, y = colnames(x = object))
  cells <- match(x = cells, table = colnames(x = object))
  if (length(x = cells) == 0) {
    warning("Cannot find cells provided")
    return(object)
  }
  idents.new <- if (length(x = value) == 1 && value %in% colnames(x = object[[]])) {
    unlist(x = object[[value]], use.names = FALSE)[cells]
  } else {
    if (is.list(x = value)) {
      value <- unlist(x = value, use.names = FALSE)
    }
    rep_len(x = value, length.out = length(x = cells))
  }
  new.levels <- if (is.factor(x = idents.new)) {
    levels(x = idents.new)
  } else {
    unique(x = idents.new)
  }
  old.levels <- levels(x = object)
  levels <- c(new.levels, old.levels)
  idents.new <- as.vector(x = idents.new)
  idents <- as.vector(x = Idents(object = object))
  idents[cells] <- idents.new
  idents[is.na(x = idents)] <- 'NA'
  levels <- intersect(x = levels, y = unique(x = idents))
  names(x = idents) <- colnames(x = object)
  missing.cells <- which(x = is.na(x = names(x = idents)))
  if (length(x = missing.cells) > 0) {
    idents <- idents[-missing.cells]
  }
  idents <- factor(x = idents, levels = levels)
  slot(object = object, name = 'active.ident') <- idents
  if (drop) {
    object <- droplevels(x = object)
  }
  return(object)
}

#' @rdname IsGlobal
#' @export
#' @method IsGlobal default
#'
IsGlobal.default <- function(object, ...) {
  return(FALSE)
}

#' @rdname IsGlobal
#' @export
#' @method IsGlobal DimReduc
#'
IsGlobal.DimReduc <- function(object, ...) {
  object <- UpdateSlots(object = object)
  return(slot(object = object, name = 'global'))
}

#' @param slot Name of slot to store JackStraw scores to
#' Can shorten to 'empirical', 'fake', 'full', or 'overall'
#'
#' @rdname JS
#' @export
#' @method JS DimReduc
#'
JS.DimReduc <- function(object, slot = NULL, ...) {
  CheckDots(...)
  jackstraw <- slot(object = object, name = 'jackstraw')
  if (!is.null(x = slot)) {
    jackstraw <- JS(object = jackstraw, slot = slot)
  }
  return(jackstraw)
}

#' @rdname JS
#' @export
#' @method JS JackStrawData
#'
JS.JackStrawData <- function(object, slot, ...) {
  CheckDots(...)
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

#' @rdname JS
#' @export
#' @method JS<- DimReduc
#'
"JS<-.DimReduc" <- function(object, slot = NULL, ..., value) {
  CheckDots(...)
  if (inherits(x = value, what = 'JackStrawData')) {
    slot(object = object, name = 'jackstraw') <- value
  } else if (is.null(x = NULL)) {
    stop("A slot must be specified")
  } else {
    JS(object = JS(object = object), slot = slot) <- value
  }
  return(object)
}

#' @rdname JS
#' @export
#' @method JS<- JackStrawData
#'
"JS<-.JackStrawData" <- function(object, slot, ..., value) {
  CheckDots(...)
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

#' @rdname Key
#' @export
#' @method Key Assay
#'
#' @examples
#' # Get an Assay key
#' Key(object = pbmc_small[["RNA"]])
#'
Key.Assay <- function(object, ...) {
  CheckDots(...)
  return(slot(object = object, name = 'key'))
}

#' @rdname Key
#' @export
#' @method Key DimReduc
#'
#' @examples
#' # Get a DimReduc key
#' Key(object = pbmc_small[["pca"]])
#'
Key.DimReduc <- function(object, ...) {
  CheckDots(...)
  return(slot(object = object, name = 'key'))
}

#' @rdname Key
#' @export
#' @method Key Seurat
#'
#' @examples
#' # Show all keys associated with a Seurat object
#' Key(object = pbmc_small)
#'
Key.Seurat <- function(object, ...) {
  CheckDots(...)
  keyed.objects <- FilterObjects(object = object)
  return(sapply(
    X = keyed.objects,
    FUN = function(x) {
    return(Key(object = object[[x]]))
    }
  ))
}

#' @rdname Key
#' @export
#' @method Key<- Assay
#'
#' @examples
#' # Set the key for an Assay
#' Key(object = pbmc_small[["RNA"]]) <- "newkey_"
#' Key(object = pbmc_small[["RNA"]])
#'
"Key<-.Assay" <- function(object, ..., value) {
  CheckDots(...)
  slot(object = object, name = 'key') <- value
  return(object)
}

#' @rdname Key
#' @export
#' @method Key<- DimReduc
#'
#' @examples
#' # Set the key for DimReduc
#' Key(object = pbmc_small[["pca"]]) <- "newkey2_"
#' Key(object = pbmc_small[["pca"]])
#'
"Key<-.DimReduc" <- function(object, ..., value) {
  CheckDots(...)
  object <- UpdateSlots(object = object)
  old.key <- Key(object = object)
  slots <- Filter(
    f = function(x) {
      return(class(x = slot(object = object, name = x)) == 'matrix')
    },
    x = slotNames(x = object)
  )
  for (s in slots) {
    mat <- slot(object = object, name = s)
    if (!IsMatrixEmpty(x = mat)) {
      colnames(x = mat) <- sub(
        pattern = paste0('^', old.key),
        replacement = value,
        x = colnames(x = mat)
      )
    }
    slot(object = object, name = s) <- mat
  }
  slot(object = object, name = 'key') <- value
  return(object)
}

#' @param projected Pull the projected feature loadings?
#'
#' @rdname Loadings
#' @export
#' @method Loadings DimReduc
#'
#' @examples
#' # Get the feature loadings for a given DimReduc
#' Loadings(object = pbmc_small[["pca"]])[1:5,1:5]
#'
Loadings.DimReduc <- function(object, projected = FALSE, ...) {
  CheckDots(...)
  projected <- projected %||% Projected(object = object)
  slot <- ifelse(
    test = projected,
    yes = 'feature.loadings.projected',
    no = 'feature.loadings'
  )
  return(slot(object = object, name = slot))
}

#' @param reduction Name of reduction to pull feature loadings for
#'
#' @rdname Loadings
#' @export
#' @method Loadings Seurat
#'
#' @examples
#' # Get the feature loadings for a specified DimReduc in a Seurat object
#' Loadings(object = pbmc_small, reduction = "pca")[1:5,1:5]
#'
Loadings.Seurat <- function(object, reduction = 'pca', projected = FALSE, ...) {
  return(Loadings(object = object[[reduction]], projected = projected, ...))
}

#' @rdname Loadings
#' @export
#' @method Loadings<- DimReduc
#'
#' @examples
#' # Set the feature loadings for a given DimReduc
#' new.loadings <- Loadings(object = pbmc_small[["pca"]])
#' new.loadings <- new.loadings + 0.01
#' Loadings(object = pbmc_small[["pca"]]) <- new.loadings
#'
"Loadings<-.DimReduc" <- function(object, projected = TRUE, ..., value) {
  CheckDots(...)
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

#' @param slot Name of specific bit of meta data to pull
#'
#' @rdname Misc
#' @export
#' @method Misc Assay
#'
Misc.Assay <- function(object, slot = NULL, ...) {
  CheckDots(...)
  if (is.null(x = slot)) {
    return(slot(object = object, name = 'misc'))
  }
  return(slot(object = object, name = 'misc')[[slot]])
}

#' @rdname Misc
#' @export
#' @method Misc Seurat
#'
#' @examples
#' # Get the misc info
#' Misc(object = pbmc_small, slot = "example")
#'
Misc.Seurat <- function(object, slot = NULL, ...) {
  CheckDots(...)
  if (is.null(x = slot)) {
    return(slot(object = object, name = 'misc'))
  }
  return(slot(object = object, name = 'misc')[[slot]])
}

#' @rdname Misc
#' @export
#' @method Misc<- Assay
#'
"Misc<-.Assay" <- function(object, slot, ..., value) {
  CheckDots(...)
  if (slot %in% names(x = Misc(object = object))) {
    warning("Overwriting miscellanous data for ", slot)
  }
  if (is.list(x = value)) {
    slot(object = object, name = 'misc')[[slot]] <- c(value)
  } else {
    slot(object = object, name = 'misc')[[slot]] <- value
  }
  return(object)
}

#' @rdname Misc
#' @export
#' @method Misc<- Seurat
#'
#' @examples
#'# Add misc info
#' Misc(object = pbmc_small, slot = "example") <- "testing_misc"
#'
"Misc<-.Seurat" <- function(object, slot, ..., value) {
  CheckDots(...)
  if (slot %in% names(x = Misc(object = object))) {
    warning("Overwriting miscellanous data for ", slot)
  }
  if (is.list(x = value)) {
    slot(object = object, name = 'misc')[[slot]] <- c(value)
  } else {
    slot(object = object, name = 'misc')[[slot]] <- value
  }
  return(object)
}

#' @param cells Subset of cell names
#' @param subset.name Parameter to subset on. Eg, the name of a gene, PC_1, a
#' column name in object@@meta.data, etc. Any argument that can be retreived
#' using FetchData
#' @param low.threshold Low cutoff for the parameter (default is -Inf)
#' @param high.threshold High cutoff for the parameter (default is Inf)
#' @param accept.value Returns all cells with the subset name equal to this value
#'
#' @rdname OldWhichCells
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
#'
#' @seealso \code{\link{FetchData}}
#'
#' @rdname OldWhichCells
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
  .Deprecated(new = "WhichCells", old = "OldWhichCells")
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
  expression <- character(length = 0L)
  if (!is.null(x = subset.name)) {
    sub <- gsub(
      pattern = '"',
      replacement = '',
      x = deparse(expr = substitute(expr = subset.name))
    )
    if (!is.infinite(x = low.threshold)) {
      expression <- c(
        expression,
        paste(sub, '>', deparse(expr = substitute(expr = low.threshold)))
      )
    }
    if (!is.infinite(x = high.threshold)) {
      expression <- c(
        expression,
        paste(sub, '<', deparse(expr = substitute(expr = high.threshold)))
      )
    }
    if (!is.null(x = accept.value)) {
      expression <- c(
        expression,
        paste(sub, '==', deparse(expr = substitute(expr = accept.value)))
      )
    }
  }
  #message(
  #  'With Seurat 3.X, identifying cells can now be done with:\n',
  #  'WhichCells(object = ',
  #  deparse(expr = substitute(expr = object)),
  #  if (length(x = expression) > 0) {
  #    paste0(', subset = ', paste(expression, collapse = ' & '))
  #  },
  #  if (!is.null(x = cells)) {
  #    paste(', cells =', deparse(expr = substitute(expr = cells)))
  #  },
  #  if (!is.null(x = ident.keep)) {
  #    paste(', idents =', deparse(expr = substitute(expr = ident.keep)))
  #  },
  #  if (!is.infinite(x = max.cells.per.ident)) {
  #    paste0(', downsample = ', max.cells.per.ident, ', seed = ', random.seed)
  #  },
  #  ')'
  #)
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

#' @rdname Project
#' @export
#' @method Project Seurat
#'
Project.Seurat <- function(object, ...) {
  CheckDots(...)
  return(slot(object = object, name = 'project.name'))
}

#' @rdname Project
#' @export
#' @method Project<- Seurat
#'
"Project<-.Seurat" <- function(object, ..., value) {
  CheckDots(...)
  slot(object = object, name = 'project.name') <- as.character(x = value)
  return(object)
}

#' @param assay Name of assay to store
#' @param layers Slot to store layers as; choose from 'counts' or 'data'; pass
#' \code{FALSE} to not pull layers; may pass one value of 'counts' or 'data' for
#' each layer in the H5AD file, must be in order
#' @param verbose Show progress updates
#'
#' @rdname h5ad
#' @export
#' @method ReadH5AD character
#'
ReadH5AD.character <- function(
  file,
  assay = 'RNA',
  layers = 'data',
  verbose = TRUE,
  ...
) {
  CheckDots(...)
  if (!PackageCheck('hdf5r', error = FALSE)) {
    stop("Please install hdf5r' for h5ad capabilities")
  }
  if (!file.exists(file)) {
    stop("Unable to find input H5AD file ", file)
  }
  hfile <- hdf5r::h5file(filename = file, mode = 'r')
  object <- ReadH5AD(
    file = hfile,
    assay = assay,
    layers = layers,
    verbose = verbose,
    ...
  )
  hfile$close_all()
  return(object)
}

#' @importFrom methods is
#' @importFrom Matrix sparseMatrix
#' @importFrom utils packageVersion
#'
#' @rdname h5ad
#' @export
#' @method ReadH5AD H5File
#'
ReadH5AD.H5File <- function(
  file,
  assay = 'RNA',
  layers = 'data',
  verbose = TRUE,
  ...
) {
  CheckDots(...)
  # Pull assay data
  # If X is an H5D, assume scaled
  # Otherwise, if file$exists(name = 'raw'), assume X is normalized
  # Otherwise, assume file[['X']] is raw counts
  if (verbose) {
    message("Pulling expression matrices and metadata")
  }
  if (is(object = file[['X']], class2 = 'H5Group')) {
    x <- as.sparse(x = file[['X']])
  } else {
    x <- file[['X']][, ]
  }
  # x will be an S3 matrix if X was scaled, otherwise will be a dgCMatrix
  scaled <- is.matrix(x = x)
  if (verbose) {
    message("Data is ", ifelse(test = scaled, yes = 'scaled', no = 'unscaled'))
  }
  # Pull cell- and feature-level metadata
  obs <- file[['obs']][]
  x.var <- file[['var']][]
  rownames(x = x) <- rownames(x = x.var) <- x.var$index
  colnames(x = x) <- rownames(x = obs) <- obs$index
  # Pull raw expression matrix and feature-level metadata
  if (file$exists(name = 'raw.X')) {
    raw <- as.sparse(x = file[['raw.X']])
    raw.var <- file[['raw.var']][]
    slot(object = raw, name = 'Dim') <- c(nrow(x = raw.var), nrow(x = obs))
    rownames(x = raw) <- rownames(x = raw.var) <- raw.var$index
    colnames(x = raw) <- obs$index
    raw.var <- raw.var[, -which(x = colnames(x = raw.var) == 'index'), drop = FALSE]
    x.slot <- ifelse(test = scaled, yes = 'scale.data', no = 'data')
  } else {
    # If X is scaled, we required normalized data present in raw
    if (scaled) {
      stop("Seurat requires normalized data present in the raw slot when X is scaled")
    } else {
      x.slot <- 'raw'
    }
  }
  obs <- obs[, -which(x = colnames(x = obs) == 'index'), drop = FALSE]
  x.var <- x.var[, -which(x = colnames(x = x.var) == 'index'), drop = FALSE]
  # Merge raw.var and x.var
  # Only happens when we have a raw.X and raw.var in the h5ad file
  if (x.slot != 'raw') {
    if (verbose) {
      message("Merging feature-level metadata dataframes")
    }
    x.var <- x.var[, -which(x = colnames(x = x.var) %in% colnames(x = raw.var))]
    meta.features <- merge(x = raw.var, y = x.var, by = 0, all = TRUE)
    rownames(x = meta.features) <- meta.features$Row.names
    meta.features <- meta.features[, -which(x = colnames(x = meta.features) == 'Row.names'), drop = FALSE]
    rm(raw.var)
  } else {
    meta.features <- x.var
  }
  # Fix meta feature colnames
  colnames(x = meta.features) <- gsub(
    pattern = 'dispersions_norm',
    replacement = 'mvp.dispersion.scaled',
    x = colnames(x = meta.features)
  )
  colnames(x = meta.features) <- gsub(
    pattern = 'dispersions',
    replacement = 'mvp.dispersion',
    x = colnames(x = meta.features)
  )
  colnames(x = meta.features) <- gsub(
    pattern = 'means',
    replacement = 'mvp.mean',
    x = colnames(x = meta.features)
  )
  colnames(x = meta.features) <- gsub(
    pattern = '_',
    replacement = '.',
    x = colnames(x = meta.features)
  )
  if ('highly.variable' %in% colnames(x = meta.features)) {
    meta.features$highly.variable[is.na(x = meta.features$highly.variable)] <- FALSE
  }
  rm(x.var)
  CheckGC()
  # Fix metadata colnames
  colnames(x = obs) <- gsub(
    pattern = '_',
    replacement = '.',
    x = colnames(x = obs)
  )
  colnames(x = obs) <- gsub(
    pattern = 'n.genes',
    replacement = paste0('nFeatures_', assay),
    x = colnames(x = obs)
  )
  colnames(x = obs) <- gsub(
    pattern = 'n.counts',
    replacement = paste0('nCount_', assay),
    x = colnames(x = obs)
  )
  # Assemble assay object
  if (verbose) {
    message("Creating assay object")
    message(
      "Storing X as ",
      x.slot,
      ifelse(
        test = x.slot != 'counts',
        yes = paste(" and raw as", ifelse(test = scaled, yes = 'data', no = 'counts')),
        no = ''
      )
    )
  }
  if (scaled) {
    assays <- list(CreateAssayObject(data = raw))
    assays[[1]] <- SetAssayData(
      object = assays[[1]],
      slot = 'scale.data',
      new.data = x
    )
    rm(raw)
  } else if (x.slot == 'data') {
    assays <- list(CreateAssayObject(counts = raw))
    assays[[1]] <- SetAssayData(
      object = assays[[1]],
      slot = 'data',
      new.data = x
    )
    rm(raw)
  } else {
    assays <- list(CreateAssayObject(counts = x))
  }
  names(x = assays) <- assay
  # Add meta feature information
  if (ncol(x = meta.features) > 0) {
    assays[[assay]][[names(x = meta.features)]] <- meta.features
  }
  # Add highly variable feature information
  if ('highly.variable' %in% colnames(x = assays[[assay]][[]])) {
    if (verbose) {
      message("Setting highly variable features")
    }
    hvf.info <- HVFInfo(object = assays[[assay]], selection.method = 'mvp')
    hvf.info <- hvf.info[order(hvf.info$dispersion, decreasing = TRUE), , drop = FALSE]
    means.use <- (hvf.info$mean > 0.1) & (hvf.info$mean < 8)
    dispersions.use <- (hvf.info$dispersion.scaled > 1) & (hvf.info$dispersion.scaled < Inf)
    top.features <- rownames(x = hvf.info)[which(x = means.use & dispersions.use)]
    VariableFeatures(object = assays[[assay]]) <- top.features
  } else if (verbose) {
    message("No variable feature expression found in h5ad file")
  }
  Key(object = assays[[assay]]) <- paste0(tolower(x = assay), '_')
  rm(x)
  CheckGC()
  # Get dimensional reduction information
  # If data isn't scaled, don't bother
  if (scaled && file$exists(name = 'obsm')) {
    if (verbose) {
      message("Pulling dimensional reduction information")
      message("Pulling cell embeddings")
    }
    # Pull cell embeddings
    if (inherits(x = file[['obsm']], what = 'H5Group')) {
      embed.reduc <- names(x = file[['obsm']])
      embeddings <- sapply(
        X = embed.reduc,
        FUN = function(x) {
          return(t(x = file[['obsm']][[x]][, ]))
        },
        simplify = FALSE,
        USE.NAMES = TRUE
      )
    } else {
      embed.reduc <- file[['obsm']]$get_type()$get_cpd_labels()
      embed.n <- sapply(
        X = file[['obsm']]$get_type()$describe()$cpd_types,
        FUN = '[[',
        'array_dims'
      )
      names(x = embed.n) <- embed.reduc
      ncells <- file[['obsm']]$dims
      embeddings <- lapply(
        X = embed.reduc,
        FUN = function(r) {
          return(t(x = vapply(
            X = 1:ncells,
            FUN = function(i) {
              return(file[['obsm']][i][[r]])
            },
            FUN.VALUE = numeric(length = embed.n[[r]])
          )))
        }
      )
      names(x = embeddings) <- embed.reduc
    }
    # Set cell names for embeddings matrices
    for (i in 1:length(x = embeddings)) {
      rownames(x = embeddings[[i]]) <- colnames(x = assays[[assay]])
    }
    # Pull feature loadings
    if (file$exists(name = 'varm')) {
      if (verbose) {
        message("Pulling feature loadings")
      }
      if (inherits(x = file[['varm']], what = 'H5Group')) {
        load.reduc <- names(x = file[['varm']])
        loadings <- sapply(
          X = load.reduc,
          FUN = function(x) {
            return(t(x = file[['varm']][[x]][, ]))
          },
          simplify = FALSE,
          USE.NAMES = TRUE
        )
      } else {
        load.reduc <- file[['varm']]$get_type()$get_cpd_labels()
        load.n <- sapply(
          X = file[['varm']]$get_type()$describe()$cpd_types,
          FUN = '[[',
          'array_dims'
        )
        names(x = load.n) <- load.reduc
        nfeatures <- file[['varm']]$dims
        loadings <- lapply(
          X = load.reduc,
          FUN = function(r) {
            return(t(x = vapply(
              X = 1:nfeatures,
              FUN = function(i) {
                return(file[['varm']][i][[r]])
              },
              FUN.VALUE = numeric(length = load.n[[load.reduc]])
            )))
          }
        )
      }
      match.ind <- lapply(
        X = gsub(pattern = 's$', replacement = '', x = tolower(x = load.reduc)),
        FUN = grep,
        x = embed.reduc
      )
      no.match <- which(x = sapply(X = match.ind, FUN = length) != 1)
      if (length(x = no.match) >= 1) {
        warning(
          "Unable to determine where the following feature loadings belong: ",
          paste(load.reduc[no.match], collapse = ', '),
          call. = FALSE,
          immediate. = TRUE
        )
        loadings <- loadings[-no.match]
        load.reduc <- load.reduc[-no.match]
        match.ind <- match.ind[-no.match]
      }
      names(x = loadings) <- embed.reduc[unlist(x = match.ind)]
      for (i in 1:length(x = loadings)) {
        rownames(x = loadings[[i]]) <- rownames(x = GetAssayData(
          object = assays[[assay]],
          slot = 'scale.data'
        ))
      }
    } else {
      if (verbose) {
        message("No feature loadings found")
      }
      loadings <- list()
    }
    # Create DimReduc objects
    dim.reducs <- vector(mode = 'list', length = length(x = embed.reduc))
    for (i in 1:length(x = embed.reduc)) {
      r <- embed.reduc[i]
      key <- tolower(x = gsub(pattern = 'X_', replacement = '', x = r))
      key <- switch(
        EXPR = key,
        'pca' = 'PC',
        'tsne' = 'tSNE',
        toupper(x = key)
      )
      key <- paste0(key, '_')
      stdev <- if (r == 'X_pca' && file$exists(name = 'uns') && file$exists(name = 'uns/pca/variance')) {
        sqrt(x = file[['uns/pca/variance']][])
      } else {
        numeric(length = 0L)
      }
      dim.reducs[[i]] <- CreateDimReducObject(
        embeddings = embeddings[[r]],
        loadings = loadings[[r]] %||% new(Class = 'matrix'),
        assay = assay,
        stdev = stdev,
        key = key
      )
    }
    # Properly name dimensional reductions
    names(x = dim.reducs) <- gsub(
      pattern = 'X_',
      replacement = '',
      x = embed.reduc
    )
    # Clean up
    rm(embeddings, loadings)
    CheckGC()
  } else {
    if (verbose) {
      message("No dimensional reduction information found")
    }
    dim.reducs <- list()
  }
  # Create the Seurat object
  if (verbose) {
    message("Assembling Seurat object")
  }
  # Create a project name, will be used as identity classes
  project <- gsub(
    pattern = '\\.h5ad',
    replacement = '',
    x = basename(path = file$filename)
  )
  object <- new(
    Class = 'Seurat',
    assays = assays,
    meta.data = obs,
    version = packageVersion(pkg = 'Seurat'),
    project.name = project
  )
  # Set default assay and identity information
  DefaultAssay(object = object) <- assay
  Idents(object = object) <- project
  # Add dimensional reduction infrom
  if (scaled && length(x = dim.reducs) >= 1) {
    for (r in names(x = dim.reducs)) {
      object[[r]] <- dim.reducs[[r]]
    }
  }
  # Get graph information
  if (scaled && file$exists(name = 'uns') && file$exists(name = 'uns/neighbors')) {
    if (verbose) {
      message("Finding nearest neighbor graph")
    }
    graph <- as.sparse(x = file[['uns/neighbors/distances']])
    colnames(x = graph) <- rownames(x = graph) <- colnames(x = object)
    method <- ifelse(
      test = file[['uns/neighbors/params']]$exists(name = 'method'),
      yes = file[['uns/neighbors/params/method']][],
      no = 'adata'
    )
    object[[paste(assay, method, sep = '_')]] <- as.Graph(x = graph)
  } else if (verbose) {
    message("No nearest-neighbor graph")
  }
  # Add layers
  if (isFALSE(x = layers)) {
    if (verbose) {
      message("Not pulling layers")
    }
  } else if (file$exists(name = 'layers')) {
    file.layers <- names(x = file[['layers']])
    layers <- rep_len(
      x = tolower(x = layers),
      length.out = length(x = file.layers)
    )
    if (!all(layers %in% c('counts', 'data'))) {
      stop("'layers' must be either 'counts' or 'data'", call. = FALSE)
    }
    names(x = layers) <- file.layers
    for (layer in file.layers) {
      layer.dest <- layers[[layer]]
      if (verbose) {
        message(
          "Reading ",
          layer,
          " into new assay, putting data into ",
          layer.dest
        )
      }
      layer.data <- if (inherits(x = file[['layers']][[layer]], what = 'H5Group')) {
        as.sparse(x = file[['layers']][[layer]])
      } else {
        file[['layers']][[layer]][, ]
      }
      dimnames(x = layer.data) <- dimnames(x = object)
      layer.assay <- switch(
        EXPR = layer.dest,
        'counts' = CreateAssayObject(
          counts = layer.data,
          min.cells = -1,
          min.features = -1
        ),
        'data' = CreateAssayObject(data = layer.data),
        stop("Unknown layer destination: ", layer.data, call. = FALSE)
      )
      object[[layer]] <- layer.assay
    }
  } else if (verbose) {
    message("No additional layers found")
  }
  return(object)
}

#' @param reverse Reverse ordering
#' @param afxn Function to evaluate each identity class based on; default is
#' \code{\link[base]{mean}}
#' @param reorder.numeric Rename all identity classes to be increasing numbers
#' starting from 1 (default is FALSE)
#'
#' @rdname Idents
#' @export
#' @method ReorderIdent Seurat
#'
ReorderIdent.Seurat <- function(
  object,
  var,
  reverse = FALSE,
  afxn = mean,
  reorder.numeric = FALSE,
  ...
) {
  data.use <- FetchData(object = object, vars = var, ...)[, 1]
  rfxn <- ifelse(
    test = reverse,
    yes = function(x) {
      return(max(x) + 1 - x)
    },
    no = Same
  )
  new.levels <- names(x = rfxn(x = sort(x = tapply(
    X = data.use,
    INDEX = Idents(object = object),
    FUN = afxn
  ))))
  new.idents <- factor(
    x = Idents(object = object),
    levels = new.levels,
    ordered = TRUE
  )
  if (reorder.numeric) {
    new.idents <- rfxn(x = rank(x = tapply(
      X = data.use,
      INDEX = as.numeric(x = new.idents),
      FUN = mean
    )))[as.numeric(x = new.idents)]
    new.idents <- factor(
      x = new.idents,
      levels = 1:length(x = new.idents),
      ordered = TRUE
    )
  }
  Idents(object = object) <- new.idents
  return(object)
}

#' @param new.names vector of new cell names
#'
#' @rdname RenameCells
#' @export
#' @method RenameCells Assay
#'
#' @examples
#' # Rename cells in an Assay
#' head(x = colnames(x = pbmc_small[["RNA"]]))
#' renamed.assay <- RenameCells(
#'     object = pbmc_small[["RNA"]],
#'     new.names = paste0("A_", colnames(x = pbmc_small[["RNA"]]))
#' )
#' head(x = colnames(x = renamed.assay))
#'
RenameCells.Assay <- function(object, new.names = NULL, ...) {
  CheckDots(...)
  if (IsSCT(assay = object)) {
    if (is.null(x = Misc(object = object, slot = 'vst.set'))) {
      suppressWarnings(Misc(object = object, slot = "vst.out")$cells_step1 <- new.names)
      suppressWarnings(rownames(x = Misc(object = object, slot = "vst.out")$cell_attr) <- new.names)
    } else{
      suppressWarnings(
        Misc(object, slot = "vst.set") <- lapply(
          X = Misc(object = object, slot = "vst.set"),
          FUN = function(x) {
            new.names.vst <- new.names[which(x = x$cells_step1 %in% Cells(x = object))]
            x$cells_step1 <- new.names.vst
            rownames(x = x$cell_attr) <- new.names.vst
            return(x)
          }
        )
      )
    }
  }
  for (data.slot in c("counts", "data", "scale.data")) {
    old.data <- GetAssayData(object = object, slot = data.slot)
    if (ncol(x = old.data) <= 1) {
      next
    }
    colnames(x = slot(object = object, name = data.slot)) <- new.names
  }
  return(object)
}

#' @rdname RenameCells
#' @export
#' @method RenameCells DimReduc
#'
#' @examples
#' # Rename cells in a DimReduc
#' head(x = Cells(x = pbmc_small[["pca"]]))
#' renamed.dimreduc <- RenameCells(
#'     object = pbmc_small[["pca"]],
#'     new.names = paste0("A_", Cells(x = pbmc_small[["pca"]]))
#' )
#' head(x = Cells(x = renamed.dimreduc))
#'
RenameCells.DimReduc <- function(object, new.names = NULL, ...) {
  CheckDots(...)
  old.data <- Embeddings(object = object)
  rownames(x = old.data) <- new.names
  slot(object = object, name = "cell.embeddings") <- old.data
  return(object)
}

#' @param for.merge Only rename slots needed for merging Seurat objects.
#' Currently only renames the raw.data and meta.data slots.
#' @param add.cell.id prefix to add cell names
#'
#' @details
#' If \code{add.cell.id} is set a prefix is added to existing cell names. If
#' \code{new.names} is set these will be used to replace existing names.
#'
#' @rdname RenameCells
#' @export
#' @method RenameCells Seurat
#'
#' @examples
#' # Rename cells in a Seurat object
#' head(x = colnames(x = pbmc_small))
#' pbmc_small <- RenameCells(object = pbmc_small, add.cell.id = "A")
#' head(x = colnames(x = pbmc_small))
#'
RenameCells.Seurat <- function(
  object,
  add.cell.id = NULL,
  new.names = NULL,
  for.merge = FALSE,
  ...
) {
  CheckDots(...)
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
  Idents(object = object) <- old.ids
  # rename the cell-level metadata
  old.meta.data <- object[[]]
  rownames(x = old.meta.data) <- new.cell.names
  slot(object = object, name = "meta.data") <- old.meta.data
  # rename the graphs
  graphs <- FilterObjects(object = object, classes.keep = "Graph")
  for (g in graphs) {
    rownames(x = object[[g]]) <- colnames(x = object[[g]]) <- new.cell.names
  }
  return(object)
}

#' @rdname Idents
#' @export
#' @method RenameIdents Seurat
#'
RenameIdents.Seurat <- function(object, ...) {
  ident.pairs <- tryCatch(
    expr = as.list(x = ...),
    error = function(e) {
      return(list(...))
    }
  )
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
  for (i in rev(x = names(x = ident.pairs))) {
    if (!i %in% names(x = cells.idents)) {
      warning("Cannot find identity ", i, call. = FALSE, immediate. = TRUE)
      next
    }
    Idents(object = object, cells = cells.idents[[i]]) <- ident.pairs[[i]]
  }
  return(object)
}

#' @param slot Where to store the new data
#' @param new.data New data to insert
#'
#'
#' @importFrom stats na.omit
#'
#' @rdname SetAssayData
#' @export
#' @method SetAssayData Assay
#'
#' @examples
#' # Set an Assay slot directly
#' count.data <- GetAssayData(object = pbmc_small[["RNA"]], slot = "counts")
#' count.data <- as.matrix(x = count.data + 1)
#' new.assay <- SetAssayData(object = pbmc_small[["RNA"]], slot = "counts", new.data = count.data)
#'
SetAssayData.Assay <- function(object, slot, new.data, ...) {
  CheckDots(...)
  slots.use <- c('counts', 'data', 'scale.data')
  if (!slot %in% slots.use) {
    stop(
      "'slot' must be one of ",
      paste(slots.use, collapse = ', '),
      call. = FALSE
    )
  }
  if (!IsMatrixEmpty(x = new.data)) {
    if (any(grepl(pattern = '_', x = rownames(x = new.data)))) {
      warning(
        "Feature names cannot have underscores ('_'), replacing with dashes ('-')",
        call. = FALSE,
        immediate. = TRUE
      )
      rownames(x = new.data) <- gsub(
        pattern = '_',
        replacement = '-',
        x = rownames(x = new.data)
      )
    }
    if (ncol(x = new.data) != ncol(x = object)) {
      stop(
        "The new data doesn't have the same number of cells as the current data",
        call. = FALSE
      )
    }
    num.counts <- nrow(x = object)
    counts.names <- rownames(x = object)
    if (slot == 'scale.data' && nrow(x = new.data) > num.counts) {
      warning(
        "Adding more features than present in current data",
        call. = FALSE,
        immediate. = TRUE
      )
    } else if (slot %in% c('counts', 'data') && nrow(x = new.data) != num.counts) {
      warning(
        "The new data doesn't have the same number of features as the current data",
        call. = FALSE,
        immediate. = TRUE
      )
    }
    if (!all(rownames(x = new.data) %in% counts.names)) {
      warning(
        "Adding features not currently present in the object",
        call. = FALSE,
        immediate. = TRUE
      )
    }
    new.features <- na.omit(object = match(
      x = counts.names,
      table = rownames(x = new.data)
    ))
    new.cells <- colnames(x = new.data)
    if (!all(new.cells %in% colnames(x = object))) {
      stop(
        "All cell names must match current cell names",
        call. = FALSE
      )
    }
    new.data <- new.data[new.features, colnames(x = object), drop = FALSE]
    if (slot %in% c('counts', 'data') && !all(dim(x = new.data) == dim(x = object))) {
      stop(
        "Attempting to add a different number of cells and/or features",
        call. = FALSE
      )
    }
  }
  if (!is.vector(x = rownames(x = new.data))) {
    rownames(x = new.data) <- as.vector(x = rownames(x = new.data))
  }
  if (!is.vector(x = colnames(x = new.data))) {
    colnames(x = new.data) <- as.vector(x = colnames(x = new.data))
  }
  slot(object = object, name = slot) <- new.data
  return(object)
}

#' @param assay Name of assay whose data should be set
#'
#' @rdname SetAssayData
#' @export
#' @method SetAssayData Seurat
#'
#' @examples
#' # Set an Assay slot through the Seurat object
#' count.data <- GetAssayData(object = pbmc_small[["RNA"]], slot = "counts")
#' count.data <- as.matrix(x = count.data + 1)
#' new.seurat.object <- SetAssayData(
#'     object = pbmc_small,
#'     slot = "counts",
#'     new.data = count.data,
#'     assay = "RNA"
#' )
#'
SetAssayData.Seurat <- function(
  object,
  slot = 'data',
  new.data,
  assay = NULL,
  ...
) {
  CheckDots(...)
  assay <- assay %||% DefaultAssay(object = object)
  object[[assay]] <- SetAssayData(object = object[[assay]], slot = slot, new.data = new.data, ...)
  return(object)
}

#' @inheritParams Idents
#'
#' @rdname Idents
#' @export
#' @method SetIdent Seurat
#'
SetIdent.Seurat <- function(object, cells = NULL, value, ...) {
  #message(
  #  'With Seurat 3.X, setting identity classes can be done as follows:\n',
  #  'Idents(object = ',
  #  deparse(expr = substitute(expr = object)),
  #  if (!is.null(x = cells)) {
  #    paste0(', cells = ', deparse(expr = substitute(expr = cells)))
  #  },
  #  ') <- ',
  #  deparse(expr = substitute(expr = value))
  #)
  CheckDots(...)
  Idents(object = object, cells = cells) <- value
  return(object)
}

#' @inheritParams Idents
#' @param save.name Store current identity information under this name
#'
#' @rdname Idents
#' @export
#' @method StashIdent Seurat
#'
StashIdent.Seurat <- function(object, save.name = 'orig.ident', ...) {
  message(
    'With Seurat 3.X, stashing identity classes can be accomplished with the following:\n',
    deparse(expr = substitute(expr = object)),
    '[[',
    deparse(expr = substitute(expr = save.name)),
    ']] <- Idents(object = ',
    deparse(expr = substitute(expr = object)),
    ')'
  )
  CheckDots(...)
  object[[save.name]] <- Idents(object = object)
  return(object)
}

#' @rdname Stdev
#' @export
#' @method Stdev DimReduc
#'
#' @examples
#' # Get the standard deviations for each PC from the DimReduc object
#' Stdev(object = pbmc_small[["pca"]])
#'
Stdev.DimReduc <- function(object, ...) {
  CheckDots(...)
  return(slot(object = object, name = 'stdev'))
}

#' @param reduction Name of reduction to use
#'
#' @rdname Stdev
#' @export
#' @method Stdev Seurat
#'
#' @examples
#' # Get the standard deviations for each PC from the Seurat object
#' Stdev(object = pbmc_small, reduction = "pca")
#'
Stdev.Seurat <- function(object, reduction = 'pca', ...) {
  CheckDots(...)
  return(Stdev(object = object[[reduction]]))
}

#' @param cells A vector of cell names to use as a subset. If NULL
#' (default), then this list will be computed based on the next three
#' arguments. Otherwise, will return an object consissting only of these cells
#' @param subset.name Parameter to subset on. Eg, the name of a gene, PC_1, a
#' column name in object@@meta.data, etc. Any argument that can be retreived
#' using FetchData
#' @param low.threshold Low cutoff for the parameter (default is -Inf)
#' @param high.threshold High cutoff for the parameter (default is Inf)
#' @param accept.value Returns cells with the subset name equal to this value
#'
#' @rdname SubsetData
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
#' @rdname SubsetData
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
  ...
) {
  .Deprecated(old = "SubsetData", new = "subset")
  expression <- character(length = 0L)
  if (!is.null(x = subset.name)) {
    sub <- gsub(
      pattern = '"',
      replacement = '',
      x = deparse(expr = substitute(expr = subset.name))
    )
    if (!is.infinite(x = low.threshold)) {
      expression <- c(
        expression,
        paste(sub, '>', deparse(expr = substitute(expr = low.threshold)))
      )
    }
    if (!is.infinite(x = high.threshold)) {
      expression <- c(
        expression,
        paste(sub, '<', deparse(expr = substitute(expr = high.threshold)))
      )
    }
    if (!is.null(x = accept.value)) {
      expression <- c(
        expression,
        paste(sub, '==', deparse(expr = substitute(expr = accept.value)))
      )
    }
  }
  #message(
  #  'With Seurat 3.X, subsetting Seurat objects can now be done with:\n',
  #  'subset(x = ',
  #  deparse(expr = substitute(expr = object)),
  #  if (length(x = expression) > 0) {
  #    paste0(', subset = ', paste(expression, collapse = ' & '))
  #  },
  #  if (length(x = c(cells, ident.use) > 0)) {
  #    paste0(', select = c("', paste0(c(cells, ident.use), collapse = '", '), '")')
  #  },
  #  if (!is.infinite(x = max.cells.per.ident)) {
  #    paste0(', downsample = ', max.cells.per.ident, ', seed = ', random.seed)
  #  },
  #  ')'
  #)
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

#' @param slot Name of tool to pull
#'
#' @rdname Tool
#' @export
#' @method Tool Seurat
#'
#' @examples
#' Tool(object = pbmc_small)
#'
Tool.Seurat <- function(object, slot = NULL, ...) {
  CheckDots(...)
  if (is.null(x = slot)) {
    return(names(x = slot(object = object, name = 'tools')))
  }
  return(slot(object = object, name = 'tools')[[slot]])
}

#' @rdname Tool
#' @export
#' @method Tool<- Seurat
#'
#' @examples
#' \dontrun{
#' sample.tool.output <- matrix(data = rnorm(n = 16), nrow = 4)
#' # must be run from within a function
#' Tool(object = pbmc_small) <- sample.tool.output
#' }
"Tool<-.Seurat" <- function(object, ..., value) {
  CheckDots(...)
  calls <- as.character(x = sys.calls())
  calls <- lapply(
    X = strsplit(x = calls, split = '(', fixed = TRUE),
    FUN = '[',
    1
  )
  tool.call <- min(grep(pattern = 'Tool<-', x = calls))
  if (tool.call <= 1) {
    stop("'Tool<-' cannot be called at the top level", call. = FALSE)
  }
  tool.call <- calls[[tool.call - 1]]
  class.call <- unlist(x = strsplit(
    x = as.character(x = sys.call())[1],
    split = '.',
    fixed = TRUE
  ))
  class.call <- class.call[length(x = class.call)]
  tool.call <- sub(
    pattern = paste0('\\.', class.call, '$'),
    replacement = '',
    x = tool.call,
    perl = TRUE
  )
  slot(object = object, name = 'tools')[[tool.call]] <- value
  return(object)
}

#' @rdname VariableFeatures
#' @export
#' @method VariableFeatures Assay
#'
VariableFeatures.Assay <- function(object, selection.method = NULL, ...) {
  CheckDots(...)
  if (!is.null(x = selection.method)) {
    vf <- HVFInfo(object = object, selection.method = selection.method, status = TRUE)
    return(rownames(x = vf)[which(x = vf[, "variable"][, 1])])
  }
  return(slot(object = object, name = 'var.features'))
}

#' @param assay Name of assay to pull variable features for
#'
#' @rdname VariableFeatures
#' @export
#' @method VariableFeatures Seurat
#'
VariableFeatures.Seurat <- function(object, assay = NULL, selection.method = NULL, ...) {
  CheckDots(...)
  assay <- assay %||% DefaultAssay(object = object)
  return(VariableFeatures(object = object[[assay]], selection.method = selection.method))
}

#' @rdname VariableFeatures
#' @export
#' @method VariableFeatures<- Assay
#'
"VariableFeatures<-.Assay" <- function(object, ..., value) {
  CheckDots(...)
  if (length(x = value) == 0) {
    slot(object = object, name = 'var.features') <- character(length = 0)
    return(object)
  }
  if (any(grepl(pattern = '_', x = value))) {
    warning(
      "Feature names cannot have underscores '_', replacing with dashes '-'",
      call. = FALSE,
      immediate = TRUE
    )
    value <- gsub(pattern = '_', replacement = '-', x = value)
  }
  value <- split(x = value, f = value %in% rownames(x = object))
  if (length(x = value[['FALSE']]) > 0) {
    if (length(x = value[['TRUE']]) == 0) {
      stop("None of the features provided are in this Assay object", call. = FALSE)
    } else {
      warning(
        "Not all features provided are in this Assay object, removing the following feature(s): ",
        paste(value[['FALSE']], collapse = ', '),
        call. = FALSE,
        immediate. = TRUE
      )
    }
  }
  slot(object = object, name = 'var.features') <- value[['TRUE']]
  return(object)
}

#' @inheritParams VariableFeatures.Seurat
#'
#' @rdname VariableFeatures
#' @export
#' @method VariableFeatures<- Seurat
#'
"VariableFeatures<-.Seurat" <- function(object, assay = NULL, ..., value) {
  CheckDots(...)
  assay <- assay %||% DefaultAssay(object = object)
  VariableFeatures(object = object[[assay]]) <- value
  return(object)
}

#' @param cells Subset of cell names
#' @param expression A predicate expression for feature/variable expression, can
#' evalue anything that can be pulled by \code{FetchData}; please note, you may
#' need to wrap feature names in backticks (\code{``}) if dashes between numbers
#' are present in the feature name
#' @param invert Invert the selection of cells
#'
#' @importFrom stats na.omit
#'
#' @rdname WhichCells
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
  CheckDots(...)
  cells <- cells %||% colnames(x = object)
  if (!missing(x = expression) && !is.null(x = substitute(expr = expression))) {
    key.pattern <- paste0('^', Key(object = object))
    expr <- if (is.call(x = substitute(expr = expression))) {
      substitute(expr = expression)
    } else {
      parse(text = expression)
    }
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
    expr.char <- gsub(
      pattern = '`',
      replacement = '',
      x = expr.char
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
  cells <- na.omit(object = unlist(x = cells, use.names = FALSE))
  return(as.character(x = cells))
}

#' @param idents A vector of identity classes to keep
#' @param slot Slot to pull feature data for
#' @param downsample Maximum number of cells per identity class, default is \code{Inf};
#' downsampling will happen after all other operations, including inverting the
#' cell selection
#' @param seed Random seed for downsampling. If NULL, does not set a seed
#'
#' @importFrom stats na.omit
#'
#' @rdname WhichCells
#' @export
#' @method WhichCells Seurat
#'
WhichCells.Seurat <- function(
  object,
  cells = NULL,
  idents = NULL,
  expression,
  slot = 'data',
  invert = FALSE,
  downsample = Inf,
  seed = 1,
  ...
) {
  CheckDots(...)
  cells <- cells %||% colnames(x = object)
  if (is.numeric(x = cells)) {
    cells <- colnames(x = object)[cells]
  }
  cell.order <- cells
  if (!is.null(x = idents)) {
    if (!is.null(x = seed)) {
      set.seed(seed = seed)
    }
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
    expr.char <- gsub(
      pattern = '`',
      replacement = '',
      x = expr.char
    )
    vars.use <- which(
      x = expr.char %in% rownames(x = object) |
        expr.char %in% colnames(x = object[[]]) |
        grepl(pattern = key.pattern, x = expr.char, perl = TRUE)
    )
    data.subset <- FetchData(
      object = object,
      vars = expr.char[vars.use],
      cells = cells,
      slot = slot
    )
    data.subset <- subset.data.frame(x = data.subset, subset = eval(expr = expr))
    cells <- rownames(x = data.subset)
  }
  if (invert) {
    cell.order <- colnames(x = object)
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
  cells <- as.character(x = na.omit(object = unlist(x = cells, use.names = FALSE)))
  cells <- cells[na.omit(object = match(x = cell.order, table = cells))]
  return(cells)
}

#' @note
#' \code{WriteH5AD} is not currently functional, please use \code{\link{as.loom}} instead
#'
#' @seealso \code{\link{as.loom}}
#'
#' @param graph Name of graph to write out, defaults to \code{paste0(assay, '_snn')}
#' @param overwrite Overwrite existing file
#'
#' @importFrom methods slot
#' @importFrom reticulate py_module_available import tuple np_array dict
#'
#' @rdname h5ad
#' @export
#' @method WriteH5AD Seurat
#'
WriteH5AD.Seurat <- function(
  object,
  file,
  assay = NULL,
  graph = NULL,
  verbose = TRUE,
  overwrite = FALSE,
  ...
) {
  message("WriteH5AD is not currently operational, please use as.loom")
  .NotYetImplemented()
  if (!PackageCheck('hdf5r', error = FALSE)) {
    stop("Please install hdf5r to enable h5ad functionality")
  }
  CheckDots(...)
  if (file.exists(file) && !overwrite) {
    stop("Output file exists, not overwriting")
  }
  assay <- assay %||% DefaultAssay(object = object)
  graph <- graph %||% paste0(assay, '_snn')
  DefaultAssay(object = object) <- assay
  object[['active_assay']] <- Idents(object = object)
  # Figure out which slot to store as X
  x.slot <- if (!IsMatrixEmpty(x = GetAssayData(object = object, slot = 'scale.data'))) {
    'scale.data'
  } else if (identical(x = GetAssayData(object = object, slot = 'counts'), y = GetAssayData(object = object, slot = 'data'))) {
    'counts'
  } else {
    'data'
  }
  if (verbose) {
    message("Storing '", x.slot, "' into 'X'")
  }
  # Figure out which slot to store as raw
  raw.slot <- switch(
    EXPR = x.slot,
    'scale.data' = 'data',
    'data' = 'counts',
    NULL
  )
  if (verbose) {
    message("Storing '", raw.slot, "' into 'raw.X'")
  }
  # Handle cases where we have data but no counts
  if (x.slot == 'counts' && IsMatrixEmpty(x = GetAssayData(object = object, slot = x.slot))) {
    if (verbose) {
      warning("Counts slot empty, storing data slot into 'X', not storing a raw.X")
    }
    x.slot <- 'data'
    raw.slot <- NULL
  }
  # Fix meta feature column names
  meta.features <- object[[assay]][[]]
  colnames(x = meta.features) <- gsub(
    pattern = 'dispersion.scaled',
    replacement = 'dispersions_norm',
    x = colnames(x = meta.features)
  )
  colnames(x = meta.features) <- gsub(
    pattern = 'dispersion',
    replacement = 'dispersions',
    x = colnames(x = meta.features)
  )
  colnames(x = meta.features) <- gsub(
    pattern = 'mean',
    replacement = 'means',
    x = colnames(x = meta.features)
  )
  colnames(x = meta.features) <- gsub(
    pattern = '\\.',
    replacement = '_',
    x = colnames(x = meta.features)
  )
  # Add variable feature information
  meta.features$highly_variable <- FALSE
  meta.features[VariableFeatures(object = object), 'highly_variable'] <- TRUE
  # Reorder feature-level metadata
  meta.features$index <- rownames(x = meta.features)
  mf.order <- c(
    'index',
    grep(
      pattern = 'index',
      x = colnames(x = meta.features),
      invert = TRUE,
      value = TRUE
    )
  )
  meta.features <- meta.features[, mf.order, drop = FALSE]
  # Fix metadata column names
  meta.data <- object[[]]
  assays.remove <- grep(
    pattern = assay,
    x = FilterObjects(object = object, classes.keep = 'Assay'),
    invert = TRUE,
    value = TRUE
  )
  if (length(x = assays.remove)) {
    assays.remove <- grep(
      pattern = assays.remove,
      x = colnames(x = meta.data)
    )
    meta.data <- meta.data[, -assays.remove, drop = FALSE]
  }
  colnames(x = meta.data) <- gsub(
    pattern = paste0('nCount_', assay),
    replacement = 'n_counts',
    x = colnames(x = meta.data)
  )
  colnames(x = meta.data) <- gsub(
    pattern = paste0('nFeatures_', assay),
    replacement = 'n_umis',
    x = colnames(x = meta.data)
  )
  colnames(x = meta.data) <- gsub(
    pattern = '\\.',
    replacement = '_',
    x = colnames(x = meta.data)
  )
  # Reorder cell-level metadata
  meta.data$index <- rownames(x = meta.data)
  md.order <- c(
    'index',
    grep(
      pattern = 'index',
      x = colnames(x = meta.data),
      invert = TRUE,
      value = TRUE
    )
  )
  meta.data <- meta.data[, md.order, drop = FALSE]
  # Write X
  hfile <- hdf5r::h5file(filename = file, mode = 'w')
  if (verbose) {
    message("Writing 'X' matrix")
  }
  x.data <- GetAssayData(object = object, slot = x.slot, assay = assay)
  switch(
    EXPR = x.slot,
    'scale.data' = hfile[['X']] <- x.data,
    {
      x.data <- as.sparse(x = x.data)
      hfile[['X/indices']] <- slot(object = x.data, 'i') - 1
      hfile[['X/indptr']] <- slot(object = x.data, 'p')
      hfile[['X/data']] <- slot(object = x.data, 'x')
    }
  )
  # Write var (feature-level metadata)
  # var only has the features that are present in X
  if (verbose) {
    message("Writing 'var' metadata")
  }
  hfile[['var']] <- meta.features[rownames(x = x.data), , drop = FALSE]
  # Write raw.X and raw.var
  if (!is.null(x = raw.slot)) {
    if (verbose) {
      message("Writing 'raw.X' sparse matrix")
    }
    # Write raw.X
    raw.data <- GetAssayData(object = object, slot = raw.slot, assay = assay)
    hfile[['raw.X/indices']] <- slot(object = raw.data, 'i') - 1
    hfile[['raw.X/indptr']] <- slot(object = raw.data, 'p')
    hfile[['raw.X/data']] <- slot(object = raw.data, 'x')
    # Write raw.var
    if (verbose) {
      message("Writing 'raw.var' metadata")
    }
    hfile[['raw.var']] <- meta.features
  }
  # Write obs (cell-level metadata)
  if (verbose) {
    message("Writing 'obs' metadata")
  }
  hfile[['obs']] <- meta.data
  # Write out dimensional reduction information
  if (x.slot == 'scale.data') {
    # Find dimensional reduction objects for this assay
    dim.reducs <- FilterObjects(object = object, classes.keep = 'DimReduc')
    dim.reducs <- Filter(
      f = function(x) {
        return(DefaultAssay(object = object[[x]]) == assay)
      },
      x = dim.reducs
    )
    # If we have any dimensional reduction objects, write them out
    if (length(x = dim.reducs) >= 1) {
      embedding.names <- paste0('X_', dim.reducs)
      names(x = embedding.names) <- dim.reducs
      loading.names <- gsub(
        pattern = '_$',
        replacement = 's',
        x = vapply(
          X = dim.reducs,
          FUN = function(x) {
            return(Key(object[[x]]))
          },
          FUN.VALUE = character(length = 1L)
        )
      )
      # TODO: Write obsm (cell embeddings)
      embeddings <- sapply(
        X = dim.reducs,
        FUN = function(x) {
          return(t(x = Embeddings(object = object, reduction = x)))
        },
        USE.NAMES = TRUE,
        simplify = FALSE
      )
      names(x = embeddings) <- embedding.names[names(x = embeddings)]
      hfile[['obsm']] <- embeddings
      # TODO: Write varm (feature loadings)
      # TODO: Write Stdev information to uns
    } else if (verbose) {
      warning("No dimensional reduction objects for assay '", assay, "' found")
    }
  } else if (verbose) {
    warning("Intial object unscaled, not storing dimensional reduction information")
  }
  # TODO: Write neighbors
  if (x.slot == 'scale.data') {
    # Find a Graph with the name provided
    graphs <- FilterObjects(object = object, classes.keep = 'Graph')
    graphs <- grep(pattern = graph, x = graphs, value = TRUE)
    # If we have a grpah, write it out
    if (length(x = graphs) == 1) {
      ''
    } else if (verbose) {
      warning("Could not find a graph named '", graph, "'")
    }
  } else if (verbose) {
    warning("Initial object unscaled, not storing graph information")
  }
  # Flush our h5ad file and return it invisibly
  hfile$flush()
  invisible(x = hfile)
  # adata <- anndata$AnnData()
  # if (!is.null(x = raw.slot)) {
  #   raw <- GetAssayData(object = object, slot = raw.slot)
  #   raw.mf <- object[[assay]][[]]
  #   raw.mf <- raw.mf[rownames(x = raw), , drop = FALSE]
  #   adata$X <- as.scipy(x = raw)
  #   adata$var <- raw.mf
  #   adata$raw <- adata
  # }
  # if (inherits(x = raw, what = c('matrix', 'Matrix'))) {
  #   raw <- as(object = raw, Class = "dgCMatrix")
  # } else {
  #   raw <- as(object = as.matrix(x = raw), Class = "dgCMatrix")
  # }
  # sp_sparse_csc <- scipy$csc_matrix
  # raw.rownames <- rownames(x = raw)
  # raw <- sp_sparse_csc(
  #   tuple(np_array(raw@x), np_array(raw@i), np_array(raw@p)),
  #   shape = tuple(raw@Dim[1], raw@Dim[2])
  # )
  # if (inherits(x = raw, what = c('matrix', 'Matrix', 'data.frame'))) {
  #   raw <- r_to_py(x = raw)
  # }
  # raw <- raw$T
  # raw <- dict(X = raw, var = dict(var_names = raw.rownames))
  # if (anndata.X == 'data') {
  #   X <- sp_sparse_csc(
  #     tuple(np_array(X@x), np_array(X@i), np_array(X@p)),
  #     shape = tuple(X@Dim[1], X@Dim[2])
  #   )
  #   X <- X$T
  # } else {
  #   X <- np_array(t(x = X))
  # }
  # obsm <- list()
  # for (dr in names(object@dr)) {
  #   obsm[[paste0("X_",dr)]] <- np_array(Embeddings(object = object[[dr]]))
  # }
  # obsm <- if (!identical(obsm, list())) dict(obsm) else NULL
  # meta_data <- object@meta.data
  # if ("nUMI" %in% colnames(x = meta_data)) {
  #   colnames(x = meta_data) <- gsub(
  #     pattern = "nUMI",
  #     replacement = "n_counts",
  #     x = colnames(x = meta_data)
  #   )
  # }
  # if ("nGene" %in% colnames(x = meta_data)) {
  #   colnames(x = meta_data) <- gsub(
  #     pattern = "nGene",
  #     replacement = "n_genes",
  #     x = colnames(x = meta_data)
  #   )
  # }
  # colnames(x = meta_data) <- gsub(
  #   pattern = "\\.",
  #   replacement = "_",
  #   x = colnames(x = meta_data)
  # )
  # anndata.object <- ad$AnnData(
  #   raw = raw,
  #   X = X,
  #   obs = meta_data,
  #   var = object@hvg.info,
  #   obsm = obsm
  # )
  # anndata.object$var_names <- gene_names
  # anndata.object$obs_names <- cell_names
  # if (!missing(x = file)) {
  #   anndata.object$write(file)
  # }
  # invisible(x = NULL)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Methods for R-defined generics
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @importFrom utils .DollarNames
#' @export
#' @method .DollarNames JackStrawData
#'
".DollarNames.JackStrawData" <- function(x, pattern = '') {
  slotnames <- as.list(x = slotNames(x = x))
  names(x = slotnames) <- unlist(x = slotnames)
  return(.DollarNames(x = slotnames, pattern = pattern))
}

#' @importFrom utils .DollarNames
#' @export
#' @method .DollarNames Seurat
#'
".DollarNames.Seurat" <- function(x, pattern = '') {
  meta.data <- as.list(x = colnames(x = x[[]]))
  names(x = meta.data) <- unlist(x = meta.data)
  return(.DollarNames(x = meta.data, pattern = pattern))
}

#' @importFrom utils .DollarNames
#' @export
#' @method .DollarNames SeuratCommand
#'
".DollarNames.SeuratCommand" <- function(x, pattern = '') {
  return(.DollarNames(x = slot(object = x, name = "params"), pattern = pattern))
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
#' @method [ Assay
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
#' @method [ DimReduc
#'
"[.DimReduc" <- function(x, i, j, drop = FALSE, ...) {
  loadings <- Loadings(object = x)
  if (missing(x = i)) {
    i <- 1:nrow(x = loadings)
  }
  if (missing(x = j)) {
    j <- names(x = x)
  } else if (is.numeric(x = j)) {
    j <- names(x = x)[j]
  }
  bad.j <- j[!j %in% colnames(x = loadings)]
  j <- j[!j %in% bad.j]
  if (length(x = j) == 0) {
    stop("None of the requested loadings are present.")
  }
  if (length(x = bad.j) > 0) {
    warning(paste0("The following loadings are not present: ", paste(bad.j, collapse = ", ")))
  }
  return(Loadings(object = x)[i, j, drop = drop, ...])
}

#' @inheritParams subset.Seurat
#'
#' @rdname subset.Seurat
#' @export
#' @method [ Seurat
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
  if (is.logical(x = i)) {
    if (length(i) != nrow(x = x)) {
      stop("Incorrect number of logical values provided to subset features")
    }
    i <- rownames(x = x)[i]
  }
  if (is.logical(x = j)) {
    if (length(j) != ncol(x = x)) {
      stop("Incorrect number of logical values provided to subset cells")
    }
    j <- colnames(x = x)[j]
  }
  if (is.numeric(x = i)) {
    i <- rownames(x = x)[i]
  }
  if (is.numeric(x = j)) {
    j <- colnames(x = x)[j]
  }
  return(subset.Seurat(x = x, features = i, cells = j, ...))
}

#' @export
#' @method [ SeuratCommand
#'
"[.SeuratCommand" <- function(x, i, ...) {
  slot.use <- c("name", "timestamp", "call_string", "params")
  if (!i %in% slot.use) {
    stop("Invalid slot")
  }
  return(slot(object = x, name = i))
}

#' @export
#' @method [[ Assay
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
#' @method [[ DimReduc
#'
"[[.DimReduc" <- function(x, i, j, drop = FALSE, ...) {
  if (missing(x = i)) {
    i <- 1:nrow(x = x)
  }
  if (missing(x = j)) {
    j <- names(x = x)
  } else if (is.numeric(x = j)) {
    j <- names(x = x)[j]
  }
  embeddings <- Embeddings(object = x)
  bad.j <- j[!j %in% colnames(x = embeddings)]
  j <- j[!j %in% bad.j]
  if (length(x = j) == 0) {
    stop("None of the requested embeddings are present.")
  }
  if (length(x = bad.j) > 0) {
    warning(paste0("The following embeddings are not present: ", paste(bad.j, collapse = ", ")))
  }
  return(embeddings[i, j, drop = drop, ...])
}

#' @export
#' @method [[ Seurat
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
      X = c('assays', 'reductions', 'graphs', 'neighbors', 'commands'),
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

#' Coerce a SeuratCommand to a list
#'
#' @inheritParams base::as.list
#' @param complete Include slots besides just parameters (eg. call string, name, timestamp)
#'
#' @return A list with the parameters and, if \code{complete = TRUE}, the call string, name, and timestamp
#'
#' @export
#' @method as.list SeuratCommand
#'
as.list.SeuratCommand <- function(x, complete = FALSE, ...) {
  CheckDots(...)
  cmd <- slot(object = x, name = 'params')
  if (complete) {
    cmd <- append(
      x = cmd,
      values = sapply(
        X = grep(
          pattern = 'params',
          x = slotNames(x = x),
          invert = TRUE,
          value = TRUE
        ),
        FUN = slot,
        object = x,
        simplify = FALSE,
        USE.NAMES = TRUE
      ),
      after = 0
    )
  }
  for (i in 1:length(x = cmd)) {
    if (is.character(x = cmd[[i]])) {
      cmd[[i]] <- paste(trimws(x = cmd[[i]]), collapse = ' ')
    }
  }
  return(cmd)
}

#' @export
#' @method as.logical JackStrawData
#'
as.logical.JackStrawData <- function(x, ...) {
  CheckDots(...)
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
  return(dim(x = Embeddings(object = x)))
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
  return(dimnames(x = Embeddings(object = x)))
}

#' @export
#' @method dimnames Seurat
#'
dimnames.Seurat <- function(x) {
  return(dimnames(x = GetAssay(object = x)))
}

#' @export
#' @method droplevels Seurat
#'
droplevels.Seurat <- function(x, ...) {
  slot(object = x, name = 'active.ident') <- droplevels(x = Idents(object = x), ...)
  return(x)
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
#' @examples
#' # Get the levels of identity classes of a Seurat object
#' levels(x = pbmc_small)
#'
levels.Seurat <- function(x) {
  return(levels(x = Idents(object = x)))
}

#' @rdname Idents
#' @export
#' @method levels<- Seurat
#'
#' @examples
#' # Reorder identity classes
#' levels(x = pbmc_small)
#' levels(x = pbmc_small) <- c('C', 'A', 'B')
#' levels(x = pbmc_small)
#'
"levels<-.Seurat" <- function(x, value) {
  idents <- Idents(object = x)
  if (!all(levels(x = idents) %in% value)) {
    stop("NA's generated by missing levels", call. = FALSE)
  }
  idents <- factor(x = idents, levels = value)
  Idents(object = x) <- idents
  return(x)
}

#' @rdname merge.Seurat
#' @export
#' @method merge Assay
#'
merge.Assay <- function(
  x = NULL,
  y = NULL,
  add.cell.ids = NULL,
  merge.data = TRUE,
  ...
) {
  CheckDots(...)
  assays <- c(x, y)
  if (!is.null(x = add.cell.ids)) {
    for (i in 1:length(assays)) {
      assays[[i]] <- RenameCells(object = assays[[i]], new.names = add.cell.ids[i])
    }
  }
  # Merge the counts (if present)
  merged.counts <- ValidateDataForMerge(assay = assays[[1]], slot = "counts")
  keys <- Key(object = assays[[1]])
  for (i in 2:length(x = assays)) {
    merged.counts <- RowMergeSparseMatrices(
      mat1 = merged.counts,
      mat2 = ValidateDataForMerge(assay = assays[[i]], slot = "counts")
    )
    if (length(Key(object = assays[[i]]) > 0)) {
      keys[i] <- Key(object = assays[[i]])
    }
  }
  combined.assay <- CreateAssayObject(
    counts = merged.counts,
    min.cells = -1,
    min.features = -1
  )
  if (length(x = unique(x = keys)) == 1) {
    Key(object = combined.assay) <- keys[1]
  }
  if (merge.data) {
    merged.data <- ValidateDataForMerge(assay = assays[[1]], slot = "data")
    for (i in 2:length(x = assays)) {
      merged.data <- RowMergeSparseMatrices(
        mat1 = merged.data,
        mat2 = ValidateDataForMerge(assay = assays[[i]], slot = "data")
      )
    }
    # only keep cells that made it through counts filtering params
    merged.data <- merged.data[, colnames(x = combined.assay)]
    combined.assay <- SetAssayData(
      object = combined.assay,
      slot = "data",
      new.data = merged.data
    )
  }
  # merge SCT assay misc vst info and scale.data
  if (all(IsSCT(assay = assays))) {
    vst.set.new <- list()
    idx <- 1
    umi.assay.new <- list()
    for (i in 1:length(x = assays)) {
      vst.set.old <- Misc(object = assays[[i]], slot = "vst.set")
      umi.assay.old <- Misc(object = assays[[i]], slot = "umi.assay")
      if (!is.null(x = vst.set.old)) {
        for (j in 1:length(x = vst.set.old)) {
          vst.set.new[[idx]] <- vst.set.old[[j]]
          umi.assay.new[[idx]] <- umi.assay.old[[j]]
          idx <- idx + 1
        }
      } else if (!is.null(x = Misc(object = assays[[i]], slot = "vst.out"))) {
        vst.set.new[[idx]] <- Misc(object = assays[[i]], slot = "vst.out")
        umi.assay.new[[idx]] <- Misc(object = assays[[i]], slot = "umi.assay")
        idx <- idx + 1
      }
    }
    Misc(object = combined.assay, slot = "vst.set") <- vst.set.new
    Misc(object = combined.assay, slot = "umi.assay") <- umi.assay.new
    scale.data <- do.call(
      what = cbind,
      args = lapply(X = assays, FUN = function(x) GetAssayData(object = x, slot = "scale.data"))
    )
    combined.assay <- SetAssayData(
      object = combined.assay,
      slot = "scale.data",
      new.data = scale.data
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
#' The merge will not preserve reductions, graphs, logged commands, or feature-level metadata
#' that were present in the original objects. If add.cell.ids isn't specified
#' and any cell names are duplicated, cell names will be appended with _X, where
#' X is the numeric index of the object in c(x, y).
#'
#' @inheritParams CreateSeuratObject
#' @param x Object
#' @param y Object (or a list of multiple objects)
#' @param add.cell.ids A character vector of length(x = c(x, y)). Appends the
#' corresponding values to the start of each objects' cell names.
#' @param merge.data Merge the data slots instead of just merging the counts
#' (which requires renormalization). This is recommended if the same normalization
#' approach was applied to all objects.
#' @param ... Arguments passed to other methods
#'
#' @return Merged object
#'
#' @rdname merge.Seurat
#' @aliases merge MergeSeurat AddSamples
#'
#' @export
#' @method merge Seurat
#'
#' @examples
#' # merge two objects
#' merge(x = pbmc_small, y = pbmc_small)
#' # to merge more than two objects, pass one to x and a list of objects to y
#' merge(x = pbmc_small, y = c(pbmc_small, pbmc_small))
#'
merge.Seurat <- function(
  x = NULL,
  y = NULL,
  add.cell.ids = NULL,
  merge.data = TRUE,
  project = "SeuratProject",
  ...
) {
  CheckDots(...)
  objects <- c(x, y)
  if (!is.null(x = add.cell.ids)) {
    if (length(x = add.cell.ids) != length(x = objects)) {
      stop("Please provide a cell identifier for each object provided to merge")
    }
    for (i in 1:length(x = objects)) {
      objects[[i]] <- RenameCells(object = objects[[i]], add.cell.id = add.cell.ids[i])
    }
  }
  # ensure unique cell names
  objects <- CheckDuplicateCellNames(object.list = objects)
  assays <- lapply(
    X = objects,
    FUN = FilterObjects,
    classes.keep = 'Assay'
  )
  fake.feature <- RandomName(length = 17)
  assays <- unique(x = unlist(x = assays, use.names = FALSE))
  combined.assays <- vector(mode = 'list', length = length(x = assays))
  names(x = combined.assays) <- assays
  for (assay in assays) {
    assays.merge <- lapply(
      X = objects,
      FUN = function(object) {
        return(tryCatch(
          expr = object[[assay]],
          error = function(e) {
            return(CreateAssayObject(counts = Matrix(
              data = 0,
              ncol = ncol(x = object),
              dimnames = list(fake.feature, colnames(x = object)),
              sparse = TRUE
            )))
          }
        ))
      }
    )
    if (all(IsSCT(assay = assays.merge))) {
      scaled.features <- unique(x = unlist(x = lapply(
        X = assays.merge,
        FUN = function(x) rownames(x = GetAssayData(object = x, slot = "scale.data")))
      ))
      for (ob in 1:length(x = objects)) {
        if (assay %in% FilterObjects(object = objects[[ob]], classes.keep = "Assay")) {
          objects[[ob]] <- suppressWarnings(GetResidual(object = objects[[ob]], features = scaled.features, assay = assay, verbose = FALSE))
          assays.merge[[ob]] <- objects[[ob]][[assay]]
        }
      }
      # handle case where some features aren't in counts and can't be retrieved with
      # GetResidual - take intersection
      scaled.features <- names(x = which(x = table(x = unlist(x = lapply(
        X = assays.merge,
        FUN = function(x) rownames(x = GetAssayData(object = x, slot = "scale.data")))
      )) == length(x = assays.merge)))
      for (a in 1:length(x = assays.merge)) {
        assays.merge[[a]] <- SetAssayData(
          object = assays.merge[[a]],
          slot = "scale.data",
          new.data = GetAssayData(object = assays.merge[[a]], slot = "scale.data")[scaled.features, ])
      }
    }
    merged.assay <- merge(
      x = assays.merge[[1]],
      y = assays.merge[2:length(x = assays.merge)],
      merge.data = merge.data
    )
    merged.assay <- subset(
      x = merged.assay,
      features = rownames(x = merged.assay)[rownames(x = merged.assay) != fake.feature]
    )
    if (length(x = Key(object = merged.assay)) == 0) {
      Key(object = merged.assay) <- paste0(assay, '_')
    }
    combined.assays[[assay]] <- merged.assay
  }
  # Merge the meta.data
  combined.meta.data <- data.frame(row.names = colnames(x = combined.assays[[1]]))
  new.idents <- c()
  for (object in objects) {
    old.meta.data <- object[[]]
    if (any(!colnames(x = old.meta.data) %in% colnames(x = combined.meta.data))) {
      cols.to.add <- colnames(x = old.meta.data)[!colnames(x = old.meta.data) %in% colnames(x = combined.meta.data)]
      combined.meta.data[, cols.to.add] <- NA
    }
    # unfactorize any factor columns
    i <- sapply(X = old.meta.data, FUN = is.factor)
    old.meta.data[i] <- lapply(X = old.meta.data[i], FUN = as.vector)
    combined.meta.data[rownames(x = old.meta.data), colnames(x = old.meta.data)] <- old.meta.data
    new.idents <- c(new.idents, as.vector(Idents(object = object)))
  }
  names(x = new.idents) <- rownames(x = combined.meta.data)
  new.idents <- factor(x = new.idents)
  if (DefaultAssay(object = x) %in% assays) {
    new.default.assay <- DefaultAssay(object = x)
  } else if (DefaultAssay(object = y) %in% assays) {
    new.default.assay <- DefaultAssay(object = y)
  } else {
    new.default.assay <- assays[1]
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
  return(merged.object)
}

#' @export
#' @method names DimReduc
#'
names.DimReduc <- function(x) {
  return(colnames(x = Embeddings(object = x)))
}

#' @export
#' @method names Seurat
#'
names.Seurat <- function(x) {
  return(FilterObjects(object = x, classes.keep = c('Assay', 'DimReduc', 'Graph')))
}

#' Print the results of a dimensional reduction analysis
#'
#' Prints a set of features that most strongly define a set of components
#'
#' @param x An object
#' @param dims Number of dimensions to display
#' @param nfeatures Number of genes to display
#' @param projected Use projected slot
#' @param ... Arguments passed to other methods
#'
#' @return Set of features defining the components
#'
#' @aliases print
#' @seealso \code{\link[base]{cat}}
#'
#' @export
#' @method print DimReduc
#'
print.DimReduc <- function(x, dims = 1:5, nfeatures = 20, projected = FALSE, ...) {
  CheckDots(...)
  loadings <- Loadings(object = x, projected = projected)
  nfeatures <- min(nfeatures, nrow(x = loadings))
  if (ncol(x = loadings) == 0) {
    warning("Dimensions have not been projected. Setting projected = FALSE")
    projected <- FALSE
    loadings <- Loadings(object = x, projected = projected)
  }
  if (min(dims) > ncol(x = loadings)) {
    stop("Cannot print dimensions greater than computed")
  }
  if (max(dims) > ncol(x = loadings)) {
    warning(paste0("Only ", ncol(x = loadings), " dimensions have been computed."))
    dims <- min(dims):ncol(x = loadings)
  }
  for (dim in dims) {
    features <- TopFeatures(
      object = x,
      dim = dim,
      nfeatures = nfeatures * 2,
      projected = projected,
      balanced = TRUE
    )
    cat(Key(object = x), dim, '\n')
    pos.features <- split(x = features$positive, f = ceiling(x = seq_along(along.with = features$positive) / 10))
    cat("Positive: ", paste(pos.features[[1]], collapse = ", "), '\n')
    pos.features[[1]] <- NULL
    if (length(x = pos.features) > 0) {
      for (i in pos.features) {
        cat("\t  ", paste(i, collapse = ", "), '\n')
      }
    }
    neg.features <- split(x = features$negative, f = ceiling(x = seq_along(along.with = features$negative) / 10))
    cat("Negative: ", paste(neg.features[[1]], collapse = ", "), '\n')
    neg.features[[1]] <- NULL
    if (length(x = neg.features) > 0) {
      for (i in neg.features) {
        cat("\t  ", paste(i, collapse = ", "), '\n')
      }
    }
  }
}

#' @importFrom stats na.omit
#'
#' @export
#' @method subset Assay
#'
subset.Assay <- function(x, cells = NULL, features = NULL, ...) {
  CheckDots(...)
  cells <- cells %||% colnames(x = x)
  if (all(is.na(x = cells))) {
    cells <- colnames(x = x)
  } else if (any(is.na(x = cells))) {
    warning("NAs passed in cells vector, removing NAs")
    cells <- na.omit(object = cells)
  }
  features <- features %||% rownames(x = x)
  if (all(is.na(x = features))) {
    features <- rownames(x = x)
  } else if (any(is.na(x = features))) {
    warning("NAs passed in the features vector, removing NAs")
    features <- na.omit(object = features)
  }
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
  cells.scaled <- cells.scaled[na.omit(object = match(x = colnames(x = x), table = cells.scaled))]
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
subset.DimReduc <- function(x, cells = NULL, features = NULL, ...) {
  CheckDots(...)
  cells <- Cells(x = x) %iff% cells %||% Cells(x = x)
  if (all(is.na(x = cells))) {
    cells <- Cells(x = x)
  } else if (any(is.na(x = cells))) {
    warning("NAs passed in cells vector, removing NAs")
    cells <- na.omit(object = cells)
  }
  # features <- rownames(x = x) %iff% features %||% rownames(x = x)
  features <- rownames(x = Loadings(object = x)) %iff% features %||% rownames(x = Loadings(object = x))
  if (all(sapply(X = list(features, cells), FUN = length) == dim(x = x))) {
    return(x)
  }
  slot(object = x, name = 'cell.embeddings') <- if (is.null(x = cells)) {
    new(Class = 'matrix')
  } else {
    if (is.numeric(x = cells)) {
      cells <- Cells(x = x)[cells]
    }
    cells <- intersect(x = cells, y = Cells(x = x))
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
  slot(object = x, name = 'feature.loadings.projected') <- if (is.null(x = features) || !Projected(object = x)) {
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
#' @param i,features A vector of features to keep
#' @param j,cells A vector of cells to keep
#' @param idents A vector of identity classes to keep
#' @param ... Extra parameters passed to \code{\link{WhichCells}},
#' such as \code{slot}, \code{invert}, or \code{downsample}
#'
#' @return A subsetted Seurat object
#'
#' @rdname subset.Seurat
#' @aliases subset
#' @seealso \code{\link[base]{subset}} \code{\link{WhichCells}}
#'
#' @export
#' @method subset Seurat
#'
#' @examples
#' subset(x = pbmc_small, subset = MS4A1 > 4)
#' subset(x = pbmc_small, subset = `DLGAP1-AS1` > 2)
#' subset(x = pbmc_small, idents = '0', invert = TRUE)
#' subset(x = pbmc_small, subset = MS4A1 > 3, slot = 'counts')
#' subset(x = pbmc_small, features = VariableFeatures(object = pbmc_small))
#'
subset.Seurat <- function(x, subset, cells = NULL, features = NULL, idents = NULL, ...) {
  if (!missing(x = subset)) {
    subset <- deparse(expr = substitute(expr = subset))
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
  if (all(cells %in% Cells(x = x)) && length(x = cells) == length(x = Cells(x = x)) && is.null(x = features)) {
    return(x)
  }
  assays <- FilterObjects(object = x, classes.keep = 'Assay')
  # Filter Assay objects
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
    stop("Under current subsetting parameters, the default assay will be removed. Please adjust subsetting parameters or change default assay.", call. = FALSE)
  }
  # Filter DimReduc objects
  for (dimreduc in FilterObjects(object = x, classes.keep = 'DimReduc')) {
    x[[dimreduc]] <- tryCatch(
      expr = subset.DimReduc(x = x[[dimreduc]], cells = cells, features = features),
      error = function(e) {
        return(NULL)
      }
    )
  }
  # Remove metadata for cells not present
  slot(object = x, name = 'meta.data') <- slot(object = x, name = 'meta.data')[cells, , drop = FALSE]
  # Recalculate nCount and nFeature
  for (assay in FilterObjects(object = x, classes.keep = 'Assay')) {
    n.calc <- CalcN(object = x[[assay]])
    if (!is.null(x = n.calc)) {
      names(x = n.calc) <- paste(names(x = n.calc), assay, sep = '_')
      x[[names(x = n.calc)]] <- n.calc
    }
  }
  slot(object = x, name = 'graphs') <- list()
  Idents(object = x, drop = TRUE) <- Idents(object = x)[cells]

  return(x)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# S4 methods
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @rdname AddMetaData
#'
setMethod(
  f = '[[<-',
  signature = c('x' = 'Assay'),
  definition = function(x, i, ..., value) {
    meta.data <- x[[]]
    feature.names <- rownames(x = meta.data)
    if (is.data.frame(x = value)) {
      value <- lapply(
        X = 1:ncol(x = value),
        FUN = function(index) {
          v <- value[[index]]
          names(x = v) <- rownames(x = value)
          return(v)
        }
      )
    }
    err.msg <- "Cannot add more or fewer meta.features information without values being named with feature names"
    if (length(x = i) > 1) {
      # Add multiple bits of feature-level metadata
      value <- rep_len(x = value, length.out = length(x = i))
      for (index in 1:length(x = i)) {
        names.intersect <- intersect(x = names(x = value[[index]]), feature.names)
        if (length(x = names.intersect) > 0) {
          meta.data[names.intersect, i[index]] <- value[[index]][names.intersect]
        } else if (length(x = value) %in% c(nrow(x = meta.data), 1) %||% is.null(x = value)) {
          meta.data[i[index]] <- value[index]
        } else {
          stop(err.msg, call. = FALSE)
        }
      }
    } else {
      # Add a single column to feature-level metadata
      value <- unlist(x = value)
      if (length(x = intersect(x = names(x = value), y = feature.names)) > 0) {
        meta.data[, i] <- value[feature.names]
      } else if (length(x = value) %in% c(nrow(x = meta.data), 1) || is.null(x = value)) {
        meta.data[, i] <- value
      } else {
        stop(err.msg, call. = FALSE)
      }
    }
    slot(object = x, name = 'meta.features') <- meta.data
    return(x)
  }
)

#' @rdname AddMetaData
#'
setMethod( # because R doesn't allow S3-style [[<- for S4 classes
  f = '[[<-',
  signature = c('x' = 'Seurat'),
  definition = function(x, i, ..., value) {
    # Require names, no index setting
    if (!is.character(x = i)) {
      stop("'i' must be a character", call. = FALSE)
    }
    # Allow removing of other object
    if (is.null(x = value)) {
      slot.use <- if (i %in% colnames(x = x[[]])) {
        'meta.data'
      } else {
        FindObject(object = x, name = i)
      }
      if (is.null(x = slot.use)) {
        stop("Cannot find object ", i, call. = FALSE)
      }
      if (i == DefaultAssay(object = x)) {
        stop("Cannot delete the default assay", call. = FALSE)
      }
    }
    # remove disallowed characters from object name
    newi <- if (is.null(x = value)) {
      i
    } else {
      make.names(names = i)
    }
    if (any(i != newi)) {
      warning(
        "Invalid name supplied, making object name syntactically valid. New object name is ",
         newi,
        "; see ?make.names for more details on syntax validity",
        call. = FALSE,
        immediate. = TRUE
      )
      i <- newi
    }
    # Figure out where to store data
    slot.use <- if (inherits(x = value, what = 'Assay')) {
      # Ensure we have the same number of cells
      if (ncol(x = value) != ncol(x = x)) {
        stop(
          "Cannot add a different number of cells than already present",
          call. = FALSE
        )
      }
      # Ensure cell order stays the same
      if (all(Cells(x = value) %in% Cells(x = x)) && !all(Cells(x = value) == Cells(x = x))) {
        for (slot in c('counts', 'data', 'scale.data')) {
          assay.data <- GetAssayData(object = value, slot = slot)
          if (!IsMatrixEmpty(x = assay.data)) {
            assay.data <- assay.data[, Cells(x = x), drop = FALSE]
          }
          # Use slot because SetAssayData is being weird
          slot(object = value, name = slot) <- assay.data
        }
      }
      'assays'
    } else if (inherits(x = value, what = 'Graph')) {
      # Ensure Assay that Graph is associated with is present in the Seurat object
      if (is.null(x = DefaultAssay(object = value))) {
        warning(
          "Adding a Graph without an assay associated with it",
          call. = FALSE,
          immediate. = TRUE
        )
      } else if (!any(DefaultAssay(object = value) %in% Assays(object = x))) {
        stop("Cannot find assay '", DefaultAssay(object = value), "' in this Seurat object", call. = FALSE)
      }
      # Ensure Graph object is in order
      if (all(Cells(x = value) %in% Cells(x = x)) && !all(Cells(x = value) == Cells(x = x))) {
        value <- value[Cells(x = x), Cells(x = x)]
      }
      'graphs'
    } else if (inherits(x = value, what = 'DimReduc')) {
      # All DimReducs must be associated with an Assay
      if (is.null(x = DefaultAssay(object = value))) {
        stop("Cannot add a DimReduc without an assay associated with it", call. = FALSE)
      }
      # Ensure Assay that DimReduc is associated with is present in the Seurat object
      if (!IsGlobal(object = value) && !any(DefaultAssay(object = value) %in% Assays(object = x))) {
        stop("Cannot find assay '", DefaultAssay(object = value), "' in this Seurat object", call. = FALSE)
      }
      # Ensure DimReduc object is in order
      if (all(Cells(x = value) %in% Cells(x = x)) && !all(Cells(x = value) == Cells(x = x))) {
        slot(object = value, name = 'cell.embeddings') <- value[[Cells(x = x), ]]
      }
      'reductions'
    } else if (inherits(x = value, what = 'SeuratCommand')) {
      # Ensure Assay that SeuratCommand is associated with is present in the Seurat object
      if (is.null(x = DefaultAssay(object = value))) {
        warning(
          "Adding a command log without an assay associated with it",
          call. = FALSE,
          immediate. = TRUE
        )
      } else if (!any(DefaultAssay(object = value) %in% Assays(object = x))) {
        stop("Cannot find assay '", DefaultAssay(object = value), "' in this Seurat object", call. = FALSE)
      }
      'commands'
    } else if (is.null(x = value)) {
      slot.use
    } else {
      'meta.data'
    }
    if (slot.use == 'meta.data') {
      # Add data to object metadata
      meta.data <- x[[]]
      cell.names <- rownames(x = meta.data)
      # If we have metadata with names, ensure they match our order
      if (is.data.frame(x = value) && !is.null(x = rownames(x = value))) {
        meta.order <- match(x = rownames(x = meta.data), table = rownames(x = value))
        value <- value[meta.order, , drop = FALSE]
      }
      if (length(x = i) > 1) {
        # Add multiple pieces of metadata
        value <- rep_len(x = value, length.out = length(x = i))
        for (index in 1:length(x = i)) {
          meta.data[i[index]] <- value[index]
        }
      } else {
        # Add a single column to metadata
        if (length(x = intersect(x = names(x = value), y = cell.names)) > 0) {
          meta.data[, i] <- value[cell.names]
        } else if (length(x = value) %in% c(nrow(x = meta.data), 1) || is.null(x = value)) {
          meta.data[, i] <- value
        } else {
          stop("Cannot add more or fewer cell meta.data information without values being named with cell names", call. = FALSE)
        }
      }
      # Check to ensure that we aren't adding duplicate names
      if (any(colnames(x = meta.data) %in% FilterObjects(object = x))) {
        bad.cols <- colnames(x = meta.data)[which(colnames(x = meta.data) %in% FilterObjects(object = x))]
        stop(
          paste0(
            "Cannot add a metadata column with the same name as an Assay or DimReduc - ",
            paste(bad.cols, collapse = ", ")),
          call. = FALSE
        )
      }
      # Store the revised metadata
      slot(object = x, name = 'meta.data') <- meta.data
    } else {
      # Add other object to Seurat object
      # Ensure cells match in value and order
      if (!(class(x = value) %in% c('SeuratCommand', 'NULL')) && !all(Cells(x = value) == Cells(x = x))) {
        stop("All cells in the object being added must match the cells in this object", call. = FALSE)
      }
      # Ensure we're not duplicating object names
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
      # Check keyed objects
      if (inherits(x = value, what = c('Assay', 'DimReduc'))) {
        if (length(x = Key(object = value)) == 0) {
          Key(object = value) <- paste0(tolower(x = i), '_')
        }
        Key(object = value) <- UpdateKey(key = Key(object = value))
        # Check for duplicate keys
        object.keys <- sapply(
          X = FilterObjects(object = x),
          FUN = function(i) {
            return(Key(object = x[[i]]))
          }
        )
        if (Key(object = value) %in% object.keys && is.null(x = FindObject(object = x, name = i))) {
          # Attempt to create a duplicate key based off the name of the object being added
          new.keys <- c(paste0(tolower(x = i), c('_', paste0(RandomName(length = 2L), '_'))))
          # Select new key to use
          key.use <- min(which(x = !new.keys %in% object.keys))
          new.key <- if (is.infinite(x = key.use)) {
            RandomName(length = 17L)
          } else {
            new.keys[key.use]
          }
          warning(
            "Cannot add objects with duplicate keys (offending key: ",
            Key(object = value),
            "), setting key to '",
            new.key,
            "'",
            call. = FALSE
          )
          # Set new key
          Key(object = value) <- new.key
        }
      }
      # For Assays, run CalcN
      if (inherits(x = value, what = 'Assay')) {
        if ((!i %in% Assays(object = x)) |
            (i %in% Assays(object = x) && ! identical(
              x = GetAssayData(object = x, assay = i, slot = "counts"),
              y = GetAssayData(object = value, slot = "counts"))
            )) {
          n.calc <- CalcN(object = value)
          if (!is.null(x = n.calc)) {
            names(x = n.calc) <- paste(names(x = n.calc), i, sep = '_')
            x[[names(x = n.calc)]] <- n.calc
          }
        }
      }
      # When removing an Assay, clear out associated DimReducs, Graphs, and SeuratCommands
      if (is.null(x = value) && inherits(x = x[[i]], what = 'Assay')) {
        objs.assay <- FilterObjects(
          object = x,
          classes.keep = c('DimReduc', 'SeuratCommand', 'Graph')
        )
        objs.assay <- Filter(
          f = function(o) {
            return(all(DefaultAssay(object = x[[o]]) == i) && !IsGlobal(object = x[[o]]))
          },
          x = objs.assay
        )
        for (o in objs.assay) {
          x[[o]] <- NULL
        }
      }
      # If adding a command, ensure it gets put at the end of the command list
      if (inherits(x = value, what = 'SeuratCommand')) {
        slot(object = x, name = slot.use)[[i]] <- NULL
        slot(object = x, name = slot.use) <- Filter(
          f = Negate(f = is.null),
          x = slot(object = x, name = slot.use)
        )
      }
      slot(object = x, name = slot.use)[[i]] <- value
      slot(object = x, name = slot.use) <- Filter(
        f = Negate(f = is.null),
        x = slot(object = x, name = slot.use)
      )
    }
    CheckGC()
    return(x)
  }
)

setMethod(
  f = 'colMeans',
  signature = c('x' = 'Assay'),
  definition = function(x, na.rm = FALSE, dims = 1, ..., slot = 'data') {
    return(Matrix::colMeans(
      x = GetAssayData(object = x, slot = slot),
      na.rm = na.rm,
      dims = dims,
      ...
    ))
  }
)

setMethod(
  f = 'colMeans',
  signature = c('x' = 'Seurat'),
  definition = function(x, na.rm = FALSE, dims = 1, ..., slot = 'data') {
    return(colMeans(
      x = GetAssayData(object = x, slot = slot),
      na.rm = na.rm,
      dims = dims,
      ...
    ))
  }
)

setMethod(
  f = 'colSums',
  signature = c('x' = 'Assay'),
  definition = function(x, na.rm = FALSE, dims = 1, ..., slot = 'data') {
    return(Matrix::colSums(
      x = GetAssayData(object = x, slot = slot),
      na.rm = na.rm,
      dims = dims,
      ...
    ))
  }
)

setMethod(
  f = 'colSums',
  signature = c('x' = 'Seurat'),
  definition = function(x, na.rm = FALSE, dims = 1, ..., slot = 'data') {
    return(Matrix::colSums(
      x = GetAssayData(object = x, slot = slot),
      na.rm = na.rm,
      dims = dims,
      ...
    ))
  }
)

setMethod(
  f = 'rowMeans',
  signature = c('x' = 'Assay'),
  definition = function(x, na.rm = FALSE, dims = 1, ..., slot = 'data') {
    return(Matrix::rowMeans(
      x = GetAssayData(object = x, slot = slot),
      na.rm = na.rm,
      dims = dims,
      ...
    ))
  }
)

setMethod(
  f = 'rowMeans',
  signature = c('x' = 'Seurat'),
  definition = function(x, na.rm = FALSE, dims = 1, ..., slot = 'data') {
    return(Matrix::rowMeans(
      x = GetAssayData(object = x, slot = slot),
      na.rm = na.rm,
      dims = dims,
      ...
    ))
  }
)

setMethod(
  f = 'rowSums',
  signature = c('x' = 'Assay'),
  definition = function(x, na.rm = FALSE, dims = 1, ..., slot = 'data') {
    return(Matrix::rowSums(
      x = GetAssayData(object = x, slot = slot),
      na.rm = na.rm,
      dims = dims,
      ...
    ))
  }
)

setMethod(
  f = 'rowSums',
  signature = c('x' = 'Seurat'),
  definition = function(x, na.rm = FALSE, dims = 1, ..., slot = 'data') {
    return(Matrix::rowSums(
      x = GetAssayData(object = x, slot = slot),
      na.rm = na.rm,
      dims = dims,
      ...
    ))
  }
)

setMethod(
  f = 'show',
  signature = 'AnchorSet',
  definition = function(object) {
    cat('An AnchorSet object containing', nrow(x = slot(object = object, name = "anchors")),
        "anchors between", length(x = slot(object = object, name = "object.list")), "Seurat objects \n",
        "This can be used as input to IntegrateData or TransferData.")
  }
)

setMethod(
  f = 'show',
  signature = 'Assay',
  definition = function(object) {
    cat('Assay data with', nrow(x = object), 'features for', ncol(x = object), 'cells\n')
    if (length(x = VariableFeatures(object = object)) > 0) {
      top.ten <- head(x = VariableFeatures(object = object), n = 10L)
      top <- 'Top'
      variable <- 'variable'
    } else {
      top.ten <- head(x = rownames(x = object), n = 10L)
      top <- 'First'
      variable <- ''
    }
    features <- paste0(
      variable,
      ' feature',
      if (length(x = top.ten) != 1) {'s'}, ":\n"
    )
    features <- gsub(pattern = '^\\s+', replacement = '', x = features)
    cat(
      top,
      length(x = top.ten),
      features,
      paste(strwrap(x = paste(top.ten, collapse = ', ')), collapse = '\n'),
      '\n'
    )
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
      'Computed using assay:', DefaultAssay(object = object), '\n'
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
    cat(
      "Active assay:",
      DefaultAssay(object = object),
      paste0('(', nrow(x = object), ' features)')
    )
    other.assays <- assays[assays != DefaultAssay(object = object)]
    if (length(x = other.assays) > 0) {
      cat(
        '\n',
        length(x = other.assays),
        'other',
        ifelse(test = length(x = other.assays) == 1, yes = 'assay', no = 'assays'),
        'present:',
        strwrap(x = paste(other.assays, collapse = ', '))
      )
    }
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
    params <- params[sapply(X = params, FUN = class) != "function"]
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

# Internal AddMetaData defintion
#
# @param object An object
# @param metadata A vector, list, or data.frame with metadata to add
# @param col.name A name for meta data if not a named list or data.frame
#
# @return object with metadata added
#
.AddMetaData <- function(object, metadata, col.name = NULL) {
  if (is.null(x = col.name) && is.atomic(x = metadata)) {
    stop("'col.name' must be provided for atomic metadata types (eg. vectors)")
  }
  if (inherits(x = metadata, what = c('matrix', 'Matrix'))) {
    metadata <- as.data.frame(x = metadata)
  }
  col.name <- col.name %||% names(x = metadata) %||% colnames(x = metadata)
  if (is.null(x = col.name)) {
    stop("No metadata name provided and could not infer it from metadata object")
  }
  object[[col.name]] <- metadata
  # if (class(x = metadata) == "data.frame") {
  #   for (ii in 1:ncol(x = metadata)) {
  #     object[[colnames(x = metadata)[ii]]] <- metadata[, ii, drop = FALSE]
  #   }
  # } else {
  #   object[[col.name]] <- metadata
  # }
  return(object)
}

# Find the names of collections in an object
#
# @return A vector with the names of slots that are a list
#
Collections <- function(object) {
  collections <- vapply(
    X = slotNames(x = object),
    FUN = function(x) {
      return(any(grepl(pattern = 'list', x = class(x = slot(object = object, name = x)))))
    },
    FUN.VALUE = logical(length = 1L)
  )
  collections <- Filter(f = isTRUE, x = collections)
  return(names(x = collections))
}

# Calculate nCount and nFeature
#
# @param object An Assay object
#
# @return A named list with nCount and nFeature
#
#' @importFrom Matrix colSums
#
CalcN <- function(object) {
  if (IsMatrixEmpty(x = GetAssayData(object = object, slot = "counts"))) {
    return(NULL)
  }
  return(list(
    nCount = colSums(x = object, slot = 'counts'),
    nFeature = colSums(x = GetAssayData(object = object, slot = 'counts') > 0)
  ))
}

# Get the names of objects within a Seurat object that are of a certain class
#
# @param object A Seurat object
# @param classes.keep A vector of names of classes to get
#
# @return A vector with the names of objects within the Seurat object that are of class \code{classes.keep}
#
FilterObjects <- function(object, classes.keep = c('Assay', 'DimReduc')) {
  slots <- na.omit(object = Filter(
    f = function(x) {
      sobj <- slot(object = object, name = x)
      return(is.list(x = sobj) && !is.data.frame(x = sobj) && !is.package_version(x = sobj))
    },
    x = slotNames(x = object)
  ))
  slots <- grep(pattern = 'tools', x = slots, value = TRUE, invert = TRUE)
  slots <- grep(pattern = 'misc', x = slots, value = TRUE, invert = TRUE)
  slots.objects <- unlist(
    x = lapply(
      X = slots,
      FUN = function(x) {
        return(names(x = slot(object = object, name = x)))
      }
    ),
    use.names = FALSE
  )
  object.classes <- sapply(
    X = slots.objects,
    FUN = function(i) {
      return(inherits(x = object[[i]], what = classes.keep))
    }
  )
  object.classes <- which(x = object.classes, useNames = TRUE)
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
  collections <- c('assays', 'graphs', 'neighbors', 'reductions', 'commands')
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

# Get the top
#
# @param data Data to pull the top from
# @param num Pull top \code{num}
# @param balanced Pull even amounts of from positive and negative values
#
# @return The top \code{num}
# @seealso \{code{\link{TopCells}}} \{code{\link{TopFeatures}}}
#
Top <- function(data, num, balanced) {
  top <- if (balanced) {
    num <- round(x = num / 2)
    data <- data[order(data, decreasing = TRUE), , drop = FALSE]
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
  cells <- colnames(x = old.assay@data)
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
    scale.data = old.assay@scale.data %||% new(Class = 'matrix'),
    meta.features = data.frame(row.names = rownames(x = counts)),
    var.features = old.assay@var.genes,
    key = paste0(assay, "_")
  )
  return(new.assay)
}

# Update dimension reduction
#
# @param old.dr Seurat2 dimension reduction slot
# @param assay.used Name of assay used to compute dimension reduction
#
UpdateDimReduction <- function(old.dr, assay) {
  new.dr <- list()
  for (i in names(x = old.dr)) {
    cell.embeddings <- old.dr[[i]]@cell.embeddings %||% new(Class = 'matrix')
    feature.loadings <- old.dr[[i]]@gene.loadings %||% new(Class = 'matrix')
    stdev <- old.dr[[i]]@sdev %||% numeric()
    misc <- old.dr[[i]]@misc %||% list()
    new.jackstraw <- UpdateJackstraw(old.jackstraw = old.dr[[i]]@jackstraw)
    old.key <- old.dr[[i]]@key
    if (length(x = old.key) == 0) {
      old.key <- gsub(pattern = "(.+?)(([0-9]+).*)", replacement = "\\1",  x = colnames(cell.embeddings)[[1]])
      if (length(x = old.key) == 0) {
        old.key <- i
      }
    }
    new.key <- suppressWarnings(expr = UpdateKey(key = old.key))
    colnames(x = cell.embeddings) <- gsub(
      pattern = old.key,
      replacement = new.key,
      x = colnames(x = cell.embeddings)
    )
    colnames(x = feature.loadings) <- gsub(
      pattern = old.key,
      replacement = new.key,
      x = colnames(x = feature.loadings)
    )
    new.dr[[i]] <- new(
      Class = 'DimReduc',
      cell.embeddings = as(object = cell.embeddings, Class = 'matrix'),
      feature.loadings = as(object = feature.loadings, Class = 'matrix'),
      assay.used = assay,
      stdev = as(object = stdev, Class = 'numeric'),
      key = as(object = new.key, Class = 'character'),
      jackstraw = new.jackstraw,
      misc = as(object = misc, Class = 'list')
    )
  }
  return(new.dr)
}

# Update jackstraw
#
# @param old.jackstraw
#
UpdateJackstraw <- function(old.jackstraw) {
  if (is.null(x = old.jackstraw)) {
    new.jackstraw <- new(
      Class = 'JackStrawData',
      empirical.p.values = new(Class = 'matrix'),
      fake.reduction.scores = new(Class = 'matrix'),
      empirical.p.values.full = new(Class = 'matrix'),
      overall.p.values = new(Class = 'matrix')
    )
  } else {
    if (.hasSlot(object = old.jackstraw, name = 'overall.p.values')) {
      overall.p <- old.jackstraw@overall.p.values %||% new(Class = 'matrix')
    } else {
      overall.p <- new(Class = 'matrix')
    }
    new.jackstraw <- new(
      Class = 'JackStrawData',
      empirical.p.values = old.jackstraw@emperical.p.value %||% new(Class = 'matrix'),
      fake.reduction.scores = old.jackstraw@fake.pc.scores %||% new(Class = 'matrix'),
      empirical.p.values.full = old.jackstraw@emperical.p.value.full %||% new(Class = 'matrix'),
      overall.p.values = overall.p
    )
  }
  return(new.jackstraw)
}

# Update a Key
#
# @param key A character to become a Seurat Key
#
# @return An updated Key that's valid for Seurat
#
UpdateKey <- function(key) {
  if (grepl(pattern = '^[[:alnum:]]+_$', x = key)) {
    return(key)
  } else {
    new.key <- regmatches(
      x = key,
      m = gregexpr(pattern = '[[:alnum:]]+', text = key)
    )
    new.key <- paste0(paste(unlist(x = new.key), collapse = ''), '_')
    if (new.key == '_') {
      new.key <- paste0(RandomName(length = 3), '_')
    }
    warning(
      "Keys should be one or more alphanumeric characters followed by an underscore, setting key from ",
      key,
      " to ",
      new.key,
       call. = FALSE,
      immediate. = TRUE
    )
    return(new.key)
  }
}

# Update slots in an object
#
# @param object An object to update
#
# @return \code{object} with the latest slot definitions
#
UpdateSlots <- function(object) {
  object.list <- sapply(
    X = slotNames(x = object),
    FUN = function(x) {
      return(tryCatch(
        expr = slot(object = object, name = x),
        error = function(...) {
          return(NULL)
        }
      ))
    },
    simplify = FALSE,
    USE.NAMES = TRUE
  )
  object.list <- Filter(f = Negate(f = is.null), x = object.list)
  object.list <- c('Class' = class(x = object)[1], object.list)
  object <- do.call(what = 'new', args = object.list)
  for (x in setdiff(x = slotNames(x = object), y = names(x = object.list))) {
    xobj <- slot(object = object, name = x)
    if (is.vector(x = xobj) && !is.list(x = xobj) && length(x = xobj) == 0) {
      slot(object = object, name = x) <- vector(mode = class(x = xobj), length = 1L)
    }
  }
  return(object)
}

# Pulls the proper data matrix for merging assay data. If the slot is empty, will return an empty
# matrix with the proper dimensions from one of the remaining data slots.
#
# @param assay Assay to pull data from
# @param slot Slot to pull from
#
# @return Returns the data matrix if present (i.e.) not 0x0. Otherwise, returns an
# appropriately sized empty sparse matrix
#
ValidateDataForMerge <- function(assay, slot) {
  mat <- GetAssayData(object = assay, slot = slot)
  if (any(dim(x = mat) == c(0, 0))) {
    slots.to.check <- setdiff(x = c("counts", "data", "scale.data"), y = slot)
    for (ss in slots.to.check) {
      data.dims <- dim(x = GetAssayData(object = assay, slot = ss))
      data.slot <- ss
      if (!any(data.dims == c(0, 0))) {
        break
      }
    }
    if (any(data.dims == c(0, 0))) {
      stop("The counts, data, and scale.data slots are all empty for the provided assay.")
    }
    mat <- Matrix(
      data = 0,
      nrow = data.dims[1],
      ncol = data.dims[2],
      dimnames = dimnames(x = GetAssayData(object = assay, slot = data.slot))
    )
  }
  return(mat)
}
