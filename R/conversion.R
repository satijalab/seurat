#' @include seurat.R
#' @import loomR
#' @importFrom methods signature
NULL

#' Convert Seurat objects to loom and vice versa
#'
#' @param from ...
#' @param to ...
#' @param filename Filename for writing loom files
#'
#' @return An object of class \code{to}
#'
#' @rdname Convert
#' @export Convert
#'
setGeneric(
  name = 'Convert',
  def = function(from, to, ...) {
    return(standardGeneric(f = 'Convert'))
  }
)

#' @rdname Convert
#' @exportMethod Convert
#'
setMethod(
  # Conversion method from Seurat to...
  f = 'Convert',
  signature = signature('from' = 'seurat', 'to' = 'character'),
  definition = function(from, to, filename, chunk.dims = 'auto') {
    object.to <- switch(
      EXPR = to,
      'loom' = {
        cell.order <- from@cell.names
        gene.order <- rownames(x = from@raw.data)
        loomfile <- create(
          filename = filename,
          data = t(x = as.matrix(x = from@raw.data[, cell.order])),
          cell.attrs = from@meta.data[cell.order, ],
          layers = list('norm_data' = t(x = from@data[, cell.order])),
          chunk.dims = chunk.dims
        )
        if (nrow(x = from@hvg.info) > 0) {
          loomfile$add.row.attribute(from@hvg.info[gene.order, ])
        }
        if (length(x = from@var.genes) > 0) {
          loomfile$add.row.attribute(list('var_genes' = gene.order %in% from@var.genes))
        }
        if (!is.null(x = from@scale.data)) {
          loomfile$add.layer(list(
            'scale_data' = as.matrix(x = t(x = as.data.frame(x = from@scale.data)[gene.order, cell.order]))
          ))
        }
        loomfile
      },
      stop(paste0("Cannot convert Seurat objects to class '", to, "'"))
    )
    return(object.to)
  }
)

#' @param raw.data Path to raw data matrix (required)
#' @param gene.names Path to gene names array (required)
#' @param cell.names Path to cell names array (required)
#' @param norm.data Path to normalized data matrix (Optional, pass NULL to skip)
#' @param scale.data Path to scaled data matrix (Optional, pass NULL to skip)
#' @param gene.means Path to gene means array (Optional, pass NULL to skip)
#' @param gene.dispersion Path to gene disperion array (Optional, pass NULL to skip)
#' @param gene.scaled Path to scaled gene dispersion array (Optional, pass NULL to skip)
#' @param var.gens Path to variable genes array (Optional, pass NULL to skip)
#' @rdname Convert
#' @exportMethod Convert
#'
setMethod(
  # Conversion method from loom to...
  f = 'Convert',
  signature = signature('from' = 'loom', 'to' = 'character'),
  definition = function(
    from,
    to,
    raw.data = 'matrix',
    gene.names = 'row_attrs/gene_names',
    cell.names = 'col_attrs/cell_names',
    norm.data = 'layers/norm_data',
    scale.data = 'layers/scale_data',
    gene.means = 'row_attrs/gene_means',
    gene.dispersion = 'row_attrs/gene_dispersion',
    gene.scaled = 'row_attrs/gene_dispersion_scaled',
    var.genes = 'row_attrs/var_genes'
  ) {
    show(object = from)
    object.to <- switch(
      EXPR = to,
      'seurat' = {
        'looooooooooooooooooom'
        # object <- CreateSeuratObject(raw.data = from[[raw.data]][, ])
        # rownames(x = object@raw.data)
      },
      stop(paste0("Cannot convert loom objects to class '", to, "'"))
    )
    return(object.to)
  }
)
