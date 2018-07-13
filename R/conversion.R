#' @include seurat.R
#' @importFrom methods signature
NULL

#' Convert Seurat objects to other classes and vice versa
#'
#' Currently, we support direct conversion to/from loom (\url{http://loompy.org/}),
#' SingleCellExperiment (\url{https://bioconductor.org/packages/release/bioc/html/SingleCellExperiment.html}),
#' and Anndata(\url{https://anndata.readthedocs.io/en/latest/}) objects.
#'
#' @param from Object to convert from
#' @param to Class of object to convert to
#' @param ... Arguments passed to and from other methods
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
#'
#' @importFrom utils installed.packages
#' @importFrom methods as slot
#' @importFrom reticulate import np_array tuple dict r_to_py
#'
#' @export
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
  anndata.X = "data",
  ...
) {
  object.to <- switch(
    EXPR = to,
    'loom' = {
      if (!'loomR' %in% rownames(x = installed.packages())) {
        stop("Please install loomR from GitHub before converting to a loom object")
      }
      cell.order <- from@cell.names
      gene.order <- rownames(x = from@raw.data)
      loomfile <- loomR::create(
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
      if (!is.null(x = from@scale.data) && dim(x = from@scale.data) != c(1, 1)) {
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
              pad.order <- seq_len(nrow(x = gene.loadings))
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
      if (!'SingleCellExperiment' %in% rownames(x = installed.packages())) {
        stop("Please install SingleCellExperiment from Bioconductor before converting to a SingeCellExperiment object")
      }
      if (inherits(x = from@raw.data, what = "data.frame")) {
        from@raw.data <- as.matrix(from@raw.data)
      }
      if (inherits(x = from@data, what = "data.frame")) {
        from@data <- as.matrix(from@data)
      }
      sce <- if (class(from@raw.data) %in% c("matrix", "dgTMatrix")) {
        SingleCellExperiment::SingleCellExperiment(assays = list(counts = as(from@raw.data[rownames(from@data), from@cell.names], "dgCMatrix")))
      } else if (inherits(x = from@raw.data, what = "dgCMatrix")) {
        SingleCellExperiment::SingleCellExperiment(assays = list(counts = from@raw.data[rownames(from@data), from@cell.names]))
      } else {
        stop("Invalid class stored in seurat object's raw.data slot")
      }
      if (class(from@data) %in% c("matrix", "dgTMatrix")) {
        SummarizedExperiment::assay(sce, "logcounts") <- as(from@data, "dgCMatrix")
      } else if (inherits(x = from@data, what = "dgCMatrix")) {
        SummarizedExperiment::assay(sce, "logcounts") <- from@data
      } else {
        stop("Invalid class stored in seurat object's data slot")
      }
      meta.data <- from@meta.data
      meta.data$ident <- from@ident
      SummarizedExperiment::colData(sce) <- S4Vectors::DataFrame(meta.data)
      row.data <- from@hvg.info[rownames(from@data), ]
      row.data <- cbind(gene = rownames(x = from@data), row.data, stringsAsFactors = FALSE)
      SummarizedExperiment::rowData(sce) <- S4Vectors::DataFrame(row.data)
      for (dr in names(from@dr)) {
        SingleCellExperiment::reducedDim(sce, toupper(x = dr)) <- slot(
          object = slot(object = from, name = "dr")[[dr]],
          name = "cell.embeddings"
        )
      }
      sce
    },
    'anndata' = {
      if (!py_module_available("anndata")) {
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
      cell_names <- colnames(x = X)
      gene_names <- rownames(x = X)
      if (inherits(x = raw, what = c('matrix', 'Matrix'))) {
        raw <- as(object = raw, Class = "dgCMatrix")
      } else {
        raw <- as(object = as.matrix(x = raw), Class = "dgCMatrix")
      }
      scipy <- import(module = 'scipy.sparse', convert = FALSE)
      sp_sparse_csc <- scipy$csc_matrix
      raw.rownames <- rownames(x = raw)
      raw <- sp_sparse_csc(
        tuple(np_array(raw@x), np_array(raw@i), np_array(raw@p)),
        shape = tuple(raw@Dim[1], raw@Dim[2])
      )
      if (inherits(x = raw, what = c('matrix', 'Matrix', 'data.frame'))) {
        raw <- r_to_py(x = raw)
      }
      raw <- raw$T
      raw <- dict(X = raw, var = dict(var_names = raw.rownames))
      if (anndata.X == 'data') {
        X <- sp_sparse_csc(
          tuple(np_array(X@x), np_array(X@i), np_array(X@p)),
          shape = tuple(X@Dim[1], X@Dim[2])
        )
        X <- X$T
      } else {
        X <- np_array(t(x = X))
      }
      obsm <- list()
      for (dr in names(from@dr)) {
        obsm[[paste0("X_",dr)]] <- np_array(GetCellEmbeddings(
          object = from,
          reduction.type = dr
        ))
      }
      obsm <- if (!identical(obsm, list())) dict(obsm) else NULL
      meta_data <- from@meta.data
      if ("nUMI" %in% colnames(x = meta_data)) {
        colnames(x = meta_data) <- gsub(
          pattern = "nUMI",
          replacement = "n_counts",
          x = colnames(x = meta_data)
        )
      }
      if ("nGene" %in% colnames(x = meta_data)) {
        colnames(x = meta_data) <- gsub(
          pattern = "nGene",
          replacement = "n_genes",
          x = colnames(x = meta_data)
        )
      }
      colnames(x = meta_data) <- gsub(
        pattern = "\\.",
        replacement = "_",
        x = colnames(x = meta_data)
      )
      anndata.object <- ad$AnnData(
        raw = raw,
        X = X,
        obs = meta_data,
        var = from@hvg.info,
        obsm = obsm
      )
      anndata.object$var_names <- gene_names
      anndata.object$obs_names <- cell_names
      if (!missing(x = filename)) {
        anndata.object$write(filename)
      }
      anndata.object
    },
    stop(paste0("Cannot convert Seurat objects to class '", to, "'"))
  )
  return(object.to)
}

#' @param raw.data.slot name of the SingleCellExperiment assay to slot into @@raw.data
#' @param data.slot name of the SingleCellExperiment assay to slot into @@data
#'
#' @describeIn Convert Convert from SingleCellExperiment to a Seurat object
#' @export
#' @method Convert SingleCellExperiment
#'
Convert.SingleCellExperiment <- function(
  from,
  to,
  raw.data.slot = "counts",
  data.slot = "logcounts",
  ...
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
      meta.data <- as.data.frame(SummarizedExperiment::colData(from))
      seurat.object <- CreateSeuratObject(raw.data = raw.data, meta.data = meta.data)
      seurat.object@data <- data
      if (length(x = SingleCellExperiment::reducedDimNames(from)) > 0) {
        for (dr in SingleCellExperiment::reducedDimNames(from)) {
          seurat.object <- SetDimReduction(
            object = seurat.object,
            reduction.type = dr,
            slot = "cell.embeddings",
            new.data = SingleCellExperiment::reducedDim(x = from, type = dr)
          )
          key <- gsub(
            pattern = "[[:digit:]]",
            replacement = "",
            x = colnames(x = SingleCellExperiment::reducedDim(x = from, type = dr)
          )[1])
          seurat.object <- SetDimReduction(
            object = seurat.object,
            reduction.type = dr,
            slot = "key",
            new.data = key
          )
        }
      }
      seurat.object
    },
    stop(paste0("Cannot convert SingleCellExperiment objects to class '", to, "'"))
  )
  return(object.to)
}

#' @param X.slot Seurat slot to transfer anndata X into. Default is scale.data
#' @param raw.slot Seurat slot to transfer anndata raw into. Default is data
#' @describeIn Convert from Anndata file to a Seurat object
#' @importFrom reticulate py_to_r
#' @export
#' @method Convert anndata.base.AnnData
#'
#'
Convert.anndata.base.AnnData <- function(
  from,
  to,
  X.slot = "scale.data",
  raw.slot = "data",
  ...
) {
  object.to <- switch(
    EXPR = to,
    'seurat' = {
      raw.data.matrix <- sparseMatrix(
        i = as.numeric(x = from$raw$X$indices),
        p = as.numeric(x = from$raw$X$indptr),
        x = as.numeric(x = from$raw$X$data),
        index1 = FALSE
      )
      rownames(x = raw.data.matrix) <- rownames(x = py_to_r(from$raw$var))
      colnames(x = raw.data.matrix) <- rownames(x = py_to_r(from$obs))
      data.matrix <- t(x = py_to_r(from$X))
      rownames(x = data.matrix) <- rownames(x = py_to_r(from$var))
      colnames(x = data.matrix) <- rownames(x = py_to_r(from$obs))
      meta.data <- py_to_r(from$obs)
      if ("n_counts" %in% colnames(x = meta.data)) {
        colnames(x = meta.data) <- gsub(
          pattern = "n_counts",
          replacement = "nUMI",
          x = colnames(x = meta.data)
        )
      }
      if ("n_gene" %in% colnames(x = meta.data)) {
        colnames(x = meta.data) <- gsub(
          pattern = "n_gene",
          replacement = "nGene",
          x = colnames(x = meta.data)
        )
      }
      seurat.object <- CreateSeuratObject(
        raw.data = raw.data.matrix,
        meta.data = meta.data
      )
      seurat.object <- SetAssayData(
        object = seurat.object,
        assay.type = "RNA",
        slot = X.slot,
        new.data = data.matrix
      )
      #todo, deal with obsm fields that are not dimensional reductions, or have different name structures
      drs <- unlist(x = py_to_r(from$obsm$keys()))
      for (dr in drs) {
        dr.embed <- py_to_r(from$obsm[[eval(dr)]])
        dr.name <- ExtractField(string = dr, field = 2)
        if (is.na(dr.name)) {
          dr.name <- dr
        }
        dr.dict <- list(tSNE_ = "tsne", PC = "pca")
        if (dr.name %in% dr.dict) {
          dr.key <- names(x = which(x = dr.dict == dr.name))
        } else {
          dr.key <- toupper(x = dr.name)
        }
        colnames(x = dr.embed) <- paste0(dr.key, 1:ncol(x = dr.embed))
        rownames(x = dr.embed) <- seurat.object@cell.names
        seurat.object <- SetDimReduction(
          object = seurat.object,
          reduction.type = dr.name,
          slot = "cell.embeddings",
          new.data = dr.embed
        )
        seurat.object <- SetDimReduction(
          object = seurat.object,
          reduction.type = dr.name,
          slot = "key",
          new.data = dr.key
        )
      }
      seurat.object
    },
    stop(paste0("Cannot convert AnnData objects to class '", to, "'"))
  )
  return(object.to)
}
