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
#'
#'@section Slots:
#'  \describe{
#'    \item{\code{raw.data}:}{\code{"ANY"}, The raw project data }
#'    \item{\code{data}:}{\code{"ANY"}, The expression matrix (log-scale) }
#'    \item{\code{scale.data}:}{\code{"ANY"}, The scaled (after z-scoring
#'    each gene) expression matrix. Used for PCA, ICA, and heatmap plotting}
#'    \item{\code{var.genes}:}{\code{"vector"},  Variable genes across single cells }
#'    \item{\code{is.expr}:}{\code{"numeric"}, Expression threshold to determine if a gene is expressed }
#'    \item{\code{ident}:}{\code{"factor"},  The 'identity class' for each single cell }
#'    \item{\code{data.info}:}{\code{"data.frame"}, Contains information about each cell, starting with # of genes detected (nGene)
#'    the original identity class (orig.ident), user-provided information (through AddMetaData), etc.  }
#'    \item{\code{project.name}:}{\code{"character"}, Name of the project (for record keeping) }
#'    \item{\code{dr}:}{\code{"list"}, List of stored dimensional reductions. Named by technique }
#'    \item{\code{assay}:}{\code{"list"}, List of additional assays for multimodal analysis. Named by technique }
#'    \item{\code{hvg.info}:}{\code{"data.frame"}, The output of the mean/variability analysis for all genes }
#'    \item{\code{imputed}:}{\code{"data.frame"}, Matrix of imputed gene scores }
#'    \item{\code{cell.names}:}{\code{"vector"},  Names of all single cells (column names of the expression matrix) }
#'    \item{\code{cluster.tree}:}{\code{"list"},  List where the first element is a phylo object containing the
#'    phylogenetic tree relating different identity classes }
#'    \item{\code{snn}:}{\code{"dgCMatrix"}, Sparse matrix object representation of the SNN graph }
#'    \item{\code{calc.params}:}{\code{"list"}, Named list to store all calculation related parameters choices}
#'    \item{\code{version}:}{\code{"ANY"}, Version of package used in object creation}
#'}
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
    data.info = "data.frame",
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




# Documentation
###############
#' @export
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


