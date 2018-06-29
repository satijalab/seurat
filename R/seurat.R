#' @importFrom utils globalVariables
#' @importFrom methods setClass setMethod setOldClass
#' @importFrom Rcpp evalCpp
#' @useDynLib Seurat
NULL

setOldClass(Classes = 'package_version')

#' The Assay Class
#'
#' The Assay object is the basic unit of Seurat; each Assay stores raw, normalized, and scaled data
#' as well as cluster information, variable features, and any other assay-specific metadata.
#' Assays should contain single cell expression data such as RNA-seq, protein, or imputed expression data.
#'
#' @slot raw.data Raw expression data
#' @slot data Normalized expression data
#' @slot scale.data Scaled expression data
#' @slot hvf.info Output of mean/variability analysis for all features
#' @slot var.features Vector of features exhibiting high variance across single cells
#' @slot dim.reduc A list of dimmensional reduction objects for this Assay
#' @slot meta.data Contains meta-information about each cell, starting with number of genes detected (nGene)
#' and the original identity class (orig.ident); more information is added using \code{AddMetaData}
#' @slot meta.features Feature-level meteadata
#' @slot cluster.tree ...
#' @slot kmeans ...
#' @slot spatial ...
#'
#' @name Assay
#' @exportClass Assay
#' @importClassesFrom Matrix dgCMatrix
#'
Assay <- setClass(
  Class = 'Assay',
  slots = c(
    raw.data = 'dgCMatrix',
    data = 'dgCMatrix',
    scale.data = 'matrix',
    ident = 'factor',
    hvf.info = 'data.frame',
    var.features = 'vector',
    dim.reduc = 'list',
    meta.data = 'data.frame',
    meta.features = 'data.frame',
    cluster.tree = 'ANY',
    kmeans = 'ANY',
    spatial = 'ANY'
  )
)

# The Dimmensional Reduction Class
#
# The DimReduc object is ...
#
# @slot cell.embeddings ...
# @slot gene.loadings ...
# @slot gene.loadings.projected ...
# @slot stdev ...
# @slot key ...
# @slot jackstraw ...
# @slot misc ...
#
DimReduc <- setClass(
  Class = 'DimReduc',
  slots = c(
    cell.embeddings = 'matrix',
    gene.loadings = 'matrix',
    gene.loadings.projected = 'matrix',
    stdev = 'numeric',
    key = 'character',
    jackstraw = 'ANY',
    misc = 'list'
  )
)

#' The Seurat Class
#'
#' The Seurat object is ...
#'
#' @slot assays A list of assays for this project
#' @slot meta.data Contains meta-information about each cell, starting with number of genes detected (nGene)
#' and the original identity class (orig.ident); more information is added using \code{AddMetaData}
#' @slot active.ident The active cluster identity for this Seurat object
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
    snn = 'list',
    project.name = 'character',
    calc.params = 'list',
    misc = 'list',
    version = 'package_version'
  )
)

#' @export
#'
MakeAssayObject <- function(
  raw.data,
  min.cells = 0,
  min.genes = 0,
  is.expr = 0,
  names.field = 1,
  names.delim = '_',
  ...
) {
  if (!inherits(x = raw.data, what = 'dgCMatrix')) {
    raw.data <- as(object = as.matrix(x = raw.data), Class = 'dgCMatrix')
  }
  if (is.expr > 0) {
    # suppress Matrix package note:
    # Note: method with signature 'CsparseMatrix#Matrix#missing#replValue' chosen for function '[<-',
    # target signature 'dgCMatrix#lgeMatrix#missing#numeric'.
    # "Matrix#ldenseMatrix#missing#replValue" would also be valid
    suppressMessages(expr = raw.data[raw.data < is.expr] <- 0)
  }
  # TODO: Add ngene and nUMI here
  # Filter based on min.genes
  num.genes <- colSums(x = raw.data > is.expr)
  raw.data <- raw.data[, which(x = num.genes > min.genes)]
  # filter genes on the number of cells expressing
  if (min.cells > 0) {
    num.cells <- rowSums(x = raw.data > 0)
    raw.data <- raw.data[which(x = num.cells >= min.cells), ]
  }
  # Set idents
  idents <- factor(x = unlist(x = lapply(
    X = colnames(x = raw.data),
    FUN = ExtractField,
    field = names.field,
    delim = names.delim
  )))
  names(x = idents) <- colnames(x = raw.data)
  # if there are more than 100 idents, set all idents to ... name
  ident.levels <- length(x = unique(x = idents))
  if (ident.levels > 100 || ident.levels == 0 || ident.levels == length(x = idents)) {
    idents <- factor(x = RandomName()) # TODO: change this
  }
  assay <- new(
    Class = 'Assay',
    raw.data = raw.data,
    data = raw.data,
    ident = idents
  )
  return(assay)
}

#' @importFrom utils packageVersion
#' @export
#'
MakeSeuratObject <- function(
  raw.data,
  project = 'SeuratProject',
  assay.name = 'RNA',
  min.cells = 0,
  min.genes = 0,
  is.expr = 0,
  normalization.method = NULL,
  scale.factor = 1e4,
  do.scale = FALSE,
  do.center = FALSE,
  names.field = 1,
  names.delim = "_",
  meta.data = NULL,
  display.progress = TRUE,
  ...
) {
  assay.list <- list()
  assay.list[assay.name] <- MakeAssayObject(
    raw.data = raw.data,
    min.cells = min.cells,
    min.genes = min.genes,
    is.expr = is.expr
  )
  object <- new(
    Class = 'Seurat',
    assays = assay.list,
    active.assay = assay.name,
    active.ident = assay.list[[assay.name]]@ident, # TODO: replace this
    project.name = project,
    version = packageVersion(pkg = 'Seurat')
  )
  # TODO: Add normalization routine
  # TODO: Add scaling routine
  # TODO: Add MetaData
  return(object)
}

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
# @exportClass seurat
# @importFrom Rcpp evalCpp
# @useDynLib Seurat
#
# seurat <- setClass(
#   "seurat",
#   slots = c(
#     raw.data = "ANY",
#     data = "ANY",
#     scale.data = "ANY",
#     var.genes = "vector",
#     is.expr = "numeric",
#     ident = "factor",
#     meta.data = "data.frame",
#     project.name = "character",
#     dr = "list",
#     assay = "list",
#     hvg.info = "data.frame",
#     imputed = "data.frame",
#     cell.names = "vector",
#     cluster.tree = "list",
#     snn = "dgCMatrix",
#     calc.params = "list",
#     kmeans = "ANY",
#     spatial = "ANY",
#     misc = "ANY",
#     version = "ANY"
#   )
# )

# show method for seurat
#
# @param object A Seurat object
# @name show
# @aliases show,seurat-method
# @docType methods
# @rdname show-methods
#
setMethod(
  f = "show",
  signature = "Seurat",
  definition = function(object) {
    cat(
      "An object of class",
      class(x = object)
    #   "in project",
    #   object@project.name,
    #   "\n",
    #   nrow(x = object@data),
    #   "genes across",
    #   ncol(x = object@data),
    #   "samples.\n"
    )
    invisible(x = NULL)
  }
)
