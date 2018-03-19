#' @include seurat.R
#' @import loomR
#' @importFrom methods signature
NULL

#' Convert Seurat objects to loom and vice versa
#'
#' @param from Object to convert from
#' @param to Class of object to convert to
#'
#' @return An object of class \code{to}
#'
#' @rdname Convert
#' @export Convert
#'
Convert <- function(from, ...) {
  UseMethod(generic = 'Convert', object = from)
}

#' @param filename Filename for writing loom files
#' @param chunk.dims Internal HDF5 chunk size
#' @param chunk.size Number of cells to stream to loom file at a time
#' @param overwrite Overwrite existing file at \code{filename}?
#'
#' @describeIn Convert Convert a Seurat object
#' @export Convert.seurat
#' @method Convert seurat
#'
Convert.seurat <- function(
  from,
  to,
  filename,
  chunk.dims = 'auto',
  chunk.size = 1000,
  overwrite = FALSE
) {
  object.to <- switch(
    EXPR = to,
    'loom' = {
      cell.order <- from@cell.names
      gene.order <- rownames(x = from@raw.data)
      loomfile <- create(
        filename = filename,
        data = from@raw.data[, cell.order],
        cell.attrs = from@meta.data[cell.order, ],
        layers = list('norm_data' = t(x = from@data[, cell.order])),
        chunk.dims = chunk.dims,
        chunk.size = chunk.size,
        overwrite = overwrite
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

#' @param raw.data Path to raw data matrix (required)
#' @param gene.names Path to gene names array (required)
#' @param cell.names Path to cell names array (required)
#' @param norm.data Path to normalized data matrix (Optional, pass NULL to skip)
#' @param scale.data Path to scaled data matrix (Optional, pass NULL to skip)
#' @param gene.means Path to gene means array (Optional, pass NULL to skip)
#' @param gene.dispersion Path to gene disperion array (Optional, pass NULL to skip)
#' @param gene.scaled Path to scaled gene dispersion array (Optional, pass NULL to skip)
#' @param var.gens Path to variable genes array (Optional, pass NULL to skip)
#'
#' @importFrom Matrix Matrix
#'
#' @describeIn Convert Convert a loom object
#' @export Convert.loom
#' @method Convert loom
#'
Convert.loom <- function(
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
        if (norm.data %in% list.datasets(object = from)) {
          norm.matrix <- t(x = from[[norm.data]][, ])
          rownames(x = norm.matrix) <- rownames(x = raw.matrix)
          colnames(x = norm.matrix) <- colnames(x = raw.matrix)
          norm.matrix <- Matrix(data = norm.matrix, sparse = TRUE)
          object@data <- norm.matrix
          gc()
        } else {
          warning(paste0("Cannot find normalized dataset '", norm.data, "'"))
        }
      }
      if (!is.null(x = scale.data)) { # Add scaled data
        if (scale.data %in% list.datasets(object = from)) {
          scale.matrix <- t(x = from[[scale.data]][, ])
          rownames(x = scale.matrix) <- rownames(x = raw.matrix)
          colnames(x = scale.matrix) <- colnames(x = raw.matrix)
          object@scale.data <- scale.matrix
          gc()
        } else {
          warning(paste0("Cannot find scaled dataset '", scale.data, "'"))
        }
      }
      # Add variable genes
      if (!is.null(x = gene.means) && !is.null(x = gene.dispersion) && !(is.null(gene.scaled))) {
        var.genes.check <- vapply(
          X = c(gene.means, gene.dispersion, gene.scaled),
          FUN = function(x) {return(x %in% list.datasets(object = from))},
          FUN.VALUE = logical(length = 1L)
        )
        if (all(var.genes.check)) {
          object@hvg.info <- data.frame(
            gene.mean = from[[gene.means]][],
            gene.dispersion = from[[gene.dispersion]][],
            gene.dispersion.scaled = from[[gene.scaled]][],
            row.names = rownames(x = raw.matrix)
          )
          if (!is.null(x = var.genes)) {
            if (var.genes %in% list.datasets(object = from)) {
              object@var.genes <- rownames(x = raw.matrix)[from[[var.genes]][]]
            } else {
              warning(paste0("Cannot find variable genes list '", var.genes, "'"))
            }
          }
          gc()
        } else {
          warning("Cannot find hvg.info data")
        }
      }
      # Add DRs
      drs <- c("pca_cell_embeddings")
      if ("pca_cell_embeddings" %in% names(x = from[['col_attrs']])) {
        pca.cell.embeddings <- t(from[['col_attrs/pca_cell_embeddings']][,])
        rownames(x = pca.cell.embeddings) <- from[[cell.names]][]
        colnames(x = pca.cell.embeddings) <- paste0("PC", 1:ncol(x = pca.cell.embeddings))
        pca.gene.loadings <- t(from[['row_attrs/pca_gene_loadings']][,])
        rownames(x = pca.gene.loadings) <- from[[gene.names]][]
        colnames(x = pca.gene.loadings) <- paste0("PC", 1:ncol(x = pca.gene.loadings))

        pca.dr <- new(
          Class = "dim.reduction",
          cell.embeddings = pca.cell.embeddings,
          gene.loadings = pca.gene.loadings,
          key = "PC"
        )
        object@dr$pca <- pca.dr
      }
      # Add meta data
      meta.data <- names(x = from[['col_attrs']])
      meta.data <- meta.data[!(meta.data %in% c(drs, basename(path = cell.names)))]
      if (length(x = meta.data) > 0) {
        row.attrs <- from[['col_attrs']]
        meta.df <- data.frame('NA' = rep.int(x = NA, times = from$shape[2]))
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
