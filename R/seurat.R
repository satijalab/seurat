#' @importFrom utils globalVariables
#' @importFrom methods setClass setMethod
NULL

################################################################################
### Seurat

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
#' @slot project.name Name of hte project (for record keeping)
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
#' @name seurat
#' @rdname seurat
#' @aliases seurat-class
#' @exportClass seurat
#' @importFrom Rcpp evalCpp
#' @useDynLib Seurat

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
  signature = "seurat",
  definition = function(object) {
    cat(
      "An object of class",
      class(object),
      "in project",
      object@project.name,
      "\n",
      nrow(x = object@data),
      "genes across",
      ncol(x = object@data),
      "samples.\n"
    )
    invisible(x = NULL)
  }
)
