#' @include seurat.R
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
#' @param display.progress Display a progress bar
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
  overwrite = FALSE,
  display.progress = TRUE
) {
  object.to <- switch(
    EXPR = to,
    'loom' = {
      tryCatch(
          expr = library(loomR),
          error = function(e) {
            stop("Please install loomR from GitHub before converting to a loom object")
          }
      )
      cell.order <- from@cell.names
      gene.order <- rownames(x = from@raw.data)
      loomfile <- create(
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
      if (!is.null(x = from@scale.data)) {
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
              pad.order <- 1:nrow(x = gene.loadings)
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
      tryCatch(
        expr = library(SingleCellExperiment),
        error = function(e) {
          stop("Please install SingleCellExperiment from Bioconductor before converting to a SingeCellExperiment object")
        }
      )
      if (class(from@raw.data) %in% c("matrix", "dgTMatrix")) {
        sce <- SingleCellExperiment(assays = list(counts = as(from@raw.data[, from@cell.names], "dgCMatrix")))
      } else {
        sce <- SingleCellExperiment(assays = list(counts = from@raw.data[, from@cell.names]))
      }

      if (class(from@data) %in% c("matrix", "dgTMatrix")) {
        assay(sce, "logcounts") <- as(from@data, "dgCMatrix")
      } else {
        assay(sce, "logcounts") <- from@data
      }
      meta.data <- from@meta.data
      meta.data$ident <- from@ident
      colData(sce) <- DataFrame(meta.data)
      for (dr in names(from@dr)){
        reducedDim(sce, dr) <-
          slot(slot(from, "dr")[[dr]], "cell.embeddings")
      }
      if (!all(dim(from@hvg.info) != c(0, 0))) {
        gene_data <- from@hvg.info
        rowData(sce) <- gene_data
      }
      sce
    },
    stop(paste0("Cannot convert Seurat objects to class '", to, "'"))
  )
  return(object.to)
}

#' @param raw.data.slot name of the SingleCellExperiment assay to slot into @@raw.data
#' @param data.slot name of the SingleCellExperiment assay to slot into @@data
#'
#' @describeIn Convert Convert from SingleCellExperiment to a Seurat object
#' @export Convert.SingleCellExperiment
#' @method Convert SingleCellExperiment
#'
Convert.SingleCellExperiment <- function(
  from,
  to,
  raw.data.slot = "counts",
  data.slot = "logcounts"
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
      meta.data <- as.data.frame(colData(from))
      seurat.object <- CreateSeuratObject(raw.data = raw.data, meta.data = meta.data)
      seurat.object@data <- data
      if(length(reducedDimNames(from)) > 0) {
        for(dr in reducedDimNames(from)) {
          seurat.object <- SetDimReduction(object = seurat.object,
                                           reduction.type = dr,
                                           slot = "cell.embeddings",
                                           new.data = reducedDim(x = from, type = dr))
        }
      }
      seurat.object
    },
    stop(paste0("Cannot convert SingleCellExperiment objects to class '", to, "'"))
  )
  return(object.to)
}
