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
#' @importFrom Matrix Matrix
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
    gene.means = 'row_attrs/gene.mean',
    gene.dispersion = 'row_attrs/gene.dispersion',
    gene.scaled = 'row_attrs/gene.dispersion.scaled',
    var.genes = 'row_attrs/var_genes'
  ) {
    object.to <- switch(
      EXPR = to,
      'seurat' = {
        raw.matrix <- t(x = from[[raw.data]][, ])
        rownames(x = raw.matrix) <- from[[gene.names]][]
        colnames(x = raw.matrix) <- from[[cell.names]][]
        raw.matrix <- Matrix(data = raw.matrix, sparse = TRUE)
        object <- CreateSeuratObject(raw.data = raw.matrix)
        if (!is.null(x = norm.data)) { # Add normalized data
          norm.matrix <- t(x = from[[norm.data]][, ])
          rownames(x = norm.matrix) <- rownames(x = raw.matrix)
          colnames(x = norm.matrix) <- colnames(x = raw.matrix)
          norm.matrix <- Matrix(data = norm.matrix, sparse = TRUE)
          object@data <- norm.matrix
          gc()
        }
        if (!is.null(x = scale.data)) { # Add scaled data
          scale.matrix <- t(x = from[[scale.data]][, ])
          rownames(x = scale.matrix) <- rownames(x = raw.matrix)
          colnames(x = scale.matrix) <- colnames(x = scale.matrix)
          object@scale.data <- scale.matrix
          gc()
        }
        # Add variable genes
        if (!is.null(x = gene.means) && !is.null(x = gene.dispersion) && !(is.null(gene.scaled))) {
          object@hvg.info <- data.frame(
            gene.mean = from[[gene.means]][],
            gene.dispersion = from[[gene.dispersion]][],
            gene.dispersion.scaled = from[[gene.scaled]][],
            row.names = rownames(x = raw.matrix)
          )
          if (!is.null(x = var.genes)) {
            object@var.genes <- rownames(x = raw.matrix)[from[[var.genes]][]]
          }
          gc()
        }
        # Add meta data
        meta.data <- names(x = from[['col_attrs']])
        meta.data <- meta.data[!(meta.data %in% basename(path = cell.names))]
        if (length(x = meta.data) > 0) {
          row.attrs <- from[['col_attrs']]
          meta.df <- data.frame('NA' = rep.int(x = NA, times = from$shape[1]))
          rownames(x = meta.df) <- colnames(x = raw.matrix)
          for (i in meta.data) {
            meta.df[, i] <- row.attrs[[i]][]
          }
          meta.df <- meta.df[, meta.data]
          object@meta.data <- meta.df
          gc()
        }
        object
      },
      stop(paste0("Cannot convert loom objects to class '", to, "'"))
    )
    return(object.to)
  }
)
