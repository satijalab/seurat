#' @importFrom methods setClass setOldClass
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
#' @slot key ...
#' @slot hvf.info Output of mean/variability analysis for all features
#' @slot var.features Vector of features exhibiting high variance across single cells
#' @slot meta.features Feature-level meteadata
#' @slot cluster.tree ...
#' @slot kmeans ...
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
    key = 'character',
    var.features = 'vector',
    meta.features = 'data.frame',
    cluster.tree = 'ANY',
    kmeans = 'ANY'
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
    jackstraw = 'ANY',
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
    commands='list'
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
    time.stamp = 'ANY',
    call.string = 'character',
    params = 'ANY'
  )
)
