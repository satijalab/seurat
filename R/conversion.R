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

#' @param filename Filename for writing files
#' @param chunk.dims Internal HDF5 chunk size
#' @param chunk.size Number of cells to stream to loom file at a time
#' @param overwrite Overwrite existing file at \code{filename}?
#' @param display.progress Display a progress bar
#' @param anndata.raw Name of matrix (raw.data, data) to put in the anndata raw slot
#' @param anndata.X Name of matrix (data, scale.data) to put in the anndata X slot
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
  display.progress = TRUE,
  anndata.raw = "raw.data",
  anndata.X = "data"
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
      row.data <- from@hvg.info[rownames(from@data), ]
      row.data <- cbind(gene = rownames(from@data), row.data)
      rowData(sce) <- DataFrame(row.data)
      for (dr in names(from@dr)){
        reducedDim(sce, toupper(dr)) <-
          slot(slot(from, "dr")[[dr]], "cell.embeddings")
      }
      sce
    },
    'anndata' = {
      if (! py_module_available("anndata")) {
        stop("Please install the anndata python module")
      }
      ad <- import("anndata")

      raw <- switch(
        EXPR = anndata.raw,
        "raw.data" = from@raw.data,
        "data" = from@data,
        stop("Invalid Seurat data slot. Please choose one of: raw.data, data")
      )
      raw <- raw[,from@cell.names]
      X <- switch(
        EXPR = anndata.X,
        "data" = from@data,
        "scale.data" = from@scale.data,
        stop("Invalid Seurat data slot. Please choose one of: data, scale.data")
      )

      cell_names <- colnames(X)
      gene_names <- rownames(X)

      if (class(raw) %in% c("matrix", "dgTMatrix")) {
        raw <- as(raw, "dgCMatrix")
      } else {
        raw <- as(as.matrix(raw), "dgCMatrix")
      }
      scipy <- import('scipy.sparse')
      sp_sparse_csc <- scipy$csc_matrix
      raw.rownames <- rownames(raw)
      raw <- sp_sparse_csc(tuple(np_array(raw@x), np_array(raw@i), np_array(raw@p)), shape = tuple(raw@Dim[1], raw@Dim[2]))
      raw <-  raw$T
      raw <- dict(X = raw, var = dict(var_names = raw.rownames))
      X <- np_array(t(X))
      obsm <- list()
      for (dr in names(from@dr)){
        obsm[[paste0("X_",dr)]] <- np_array(GetCellEmbeddings(from, reduction.type = dr))
      }
      obsm <- dict(obsm)

      meta_data <- from@meta.data
      colnames(meta_data) <- gsub("nUMI","n_counts",colnames(meta_data))
      colnames(meta_data) <- gsub("\\.","_",colnames(meta_data))
      colnames(meta_data) <- gsub("nGene","n_genes",colnames(meta_data))

      anndata.object <- ad$AnnData(raw = raw,
                                   X = X,
                                   obs = meta_data,
                                   var = from@hvg.info,
                                   obsm = obsm)

      anndata.object$var_names <- gene_names
      anndata.object$obs_names <- cell_names
      if (! missing(filename)) {
        anndata.object$write(filename)
      }
      anndata.object
    },
    stop(paste0("Cannot convert Seurat objects to class '", to, "'"))
  )
  return(object.to)
}

setAs(
  from = "seurat",
  to = "SingleCellExperiment",
  def = function(from) {
    return(Convert(from = from, to = "sce"))
  }
)

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
          seurat.object <- SetDimReduction(object = seurat.object,
                                           reduction.type = dr,
                                           slot = "key",
                                           new.data = dr)
        }
      }
      seurat.object
    },
    stop(paste0("Cannot convert SingleCellExperiment objects to class '", to, "'"))
  )
  return(object.to)
}

setAs(
  from = "SingleCellExperiment",
  to = "seurat",
  def = function(from) {
    return(Convert(from = from, to = "seurat"))
  }
)

#' @param filename.from filename holding the AnnData object, should end in h5ad
#' @param X.slot Which Seurat slot should AnnData.X be transferred into? Default is data
#' @describeIn Convert from Anndata file to a Seurat object
#' @export Convert.AnnData
#' @method Convert AnnData
#'
#'
Convert.AnnData <- function(
  filename.from,
  to,
  X.slot = "data"
) {
  object.to <- switch(
    EXPR = to,
    'seurat' = {
      #to do - make sure we load in hdf5r functions
      h5 <- h5file(filename)
      raw_data <- h5[["raw.X"]]

      #to do, make sure we properly import Matrix and sparseMatrix
      raw_data_matrix <- Matrix(sparseMatrix(i=raw_data[["indices"]][],p=raw_data[["indptr"]][],x=raw_data[["data"]][],index1 = F),sparse=T)
      rownames(raw_data_matrix) <- h5[["raw.var"]][][,1]
      colnames(raw_data_matrix) <- h5[["obs"]][][,1]

      data_matrix <- h5[["X"]][,]
      rownames(data_matrix) <- h5[["var"]][][,1]
      colnames(data_matrix) <- h5[["obs"]][][,1]

      meta_data <- h5[["obs"]][]
      rownames(meta_data) <- meta_data$index
      meta_data <- meta_data[,-1]
      colnames(meta_data) <- gsub("n_counts", "nUMI",colnames(meta_data))
      colnames(meta_data) <- gsub("n_gene","nGene",colnames(meta_data))
      object.to <- CreateSeuratObject(raw.data = raw_data_matrix, meta.data = meta_data)

      #todo, check for sparse matrix possibilities
      object.to <- SetAssayData(object, "RNA", X.slot, data_matrix)

      #todo, deal with obsm fields that are not dimensional reductions, or have different name structures
      drs <- names(h5[["obsm"]][])
      dr_names <- sapply(drs,function(x)ExtractField(x,2))
      for(i in 1:length(dr_names)) {
        dr_embed <- matrix(h5[["obsm"]][][i][,1],nrow = length(object.to@cell.names))
        dr <- dr_names[i]
        #todo: translate DR codes
        dr_key <- toupper(dr)
        colnames(dr_embed) <- paste0(dr_key,1:ncol(dr_embed))
        rownames(dr_embed) <- object.to@cell.names
        object.to <- SetDimReduction(object.to,dr,"cell.embeddings",dr_embed)
        object.to <- SetDimReduction(object.to,dr,"key",dr_key)
      }
      return(object.to)
    },
    stop(paste0("Cannot convert AnnData objects to class '", to, "'"))
  )
  return(object.to)
}
