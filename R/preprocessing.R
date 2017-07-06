#' Setup Seurat object
#'
#' Setup and initialize basic parameters of the Seurat object
#'
#' @param object Seurat object
#' @param project Project name (string)
#' @param min.cells Include genes with detected expression in at least this
#' many cells
#' @param min.genes Include cells where at least this many genes are detected
#' @param is.expr Expression threshold for 'detected' gene
#' @param do.logNormalize whether to normalize the expression data per cell and transform to log space.
#' @param total.expr scale factor in the log normalization
#' @param do.scale In object@@scale.data, perform row-scaling (gene-based
#' z-score)
#' @param do.center In object@@scale.data, perform row-centering (gene-based
#' centering)
#' @param names.field For the initial identity class for each cell, choose this
#' field from the cell's column name
#' @param names.delim For the initial identity class for each cell, choose this
#' delimiter from the cell's column name
#' @param meta.data Additional metadata to add to the Seurat object. Should be
#' a data frame where the rows are cell names, and the columns are additional
#' metadata fields
#' @param save.raw TRUE by default. If FALSE, do not save the unmodified data in object@@raw.data
#' which will save memory downstream for large datasets
#'
#' @return Seurat object. Fields modified include object@@data,
#' object@@scale.data, object@@data.info, object@@ident
#'
#' @import stringr
#' @import pbapply
#' @importFrom Matrix colSums rowSums
#'
#' @export
#'
Setup <- function(
  object,
  project = "default",
  min.cells = 3,
  min.genes = 1000,
  is.expr = 0,
  do.logNormalize = TRUE,
  total.expr = 1e4,
  do.scale = TRUE,
  do.center = TRUE,
  names.field = 1,
  names.delim = "_",
  meta.data = NULL,
  save.raw = TRUE
) {
  object@project.name <- project
  object@is.expr <- is.expr
  num.genes <- colSums(object@raw.data > is.expr)
  num.mol=colSums(object@raw.data)
  cells.use <- names(num.genes[which(num.genes > min.genes)])
  object@data <- object@raw.data[, cells.use]
  #to save memory downstream, especially for large object
  if (!(save.raw)) {
    object@raw.data <- matrix()
  }
  genes.use <- rownames(object@data)
  if (min.cells > 0) {
    num.cells <- rowSums(object@data > is.expr)
    genes.use <- names(num.cells[which(num.cells >= min.cells)])
    object@data <- object@data[genes.use, ]
  }
  if (do.logNormalize) {
    object@data=LogNormalize(object@data,scale.factor = total.expr)
  }
  object@ident <- factor(
    x = unlist(
      x = lapply(
        X = colnames(x = object@data),
        FUN = extract_field,
        field = names.field,
        delim = names.delim
      )
    )
  )
  names(x = object@ident) <- colnames(x = object@data)
  object@cell.names <- names(x = object@ident)
  # if there are more than 100 idents, set all ident to project name
  ident.levels <- length(x = unique(x = object@ident))
  if ((ident.levels > 100 || ident.levels == 0) || ident.levels == length(x = object@ident)) {
    object <- SetIdent(object, ident.use = project)
  }
  data.ngene <- num.genes[cells.use]
  data.nmol <- num.mol[cells.use]
  object@gene.scores <- data.frame(data.ngene)
  colnames(object@gene.scores)[1] <- "nGene"
  nGene <- data.ngene
  nUMI <- data.nmol
  object@data.info <- data.frame(nGene, nUMI)
  if (! is.null(x = meta.data)) {
    object <- AddMetaData(object = object, metadata = meta.data)
  }
  object@mix.probs <- data.frame(data.ngene)
  colnames(x = object@mix.probs)[1] <- "nGene"
  rownames(x = object@gene.scores) <- colnames(x = object@data)
  object@data.info[names(object@ident),"orig.ident"] <- object@ident
  object@scale.data <- matrix()
  object@assay <- list()
  if(do.scale | do.center) {
    object <- ScaleData(object = object, do.scale = do.scale, do.center = do.center)
  }
  #if(calc.noise) {
  #  object=CalcNoiseModels(object,...)
  #  object=GetWeightMatrix(object)
  #}
  return(object)
}

#' Load in data from 10X
#'
#' Enables easy loading of sparse data matrices provided by 10X genomics.
#'
#' @param data.dir Directory containing the matrix.mtx, genes.tsv, and barcodes.tsv
#' files provided by 10X. A vector or named vector can be given in order to load
#' several data directories. If a named vector is given, the cell barcode names
#' will be prefixed with the name.
#'
#' @return Returns a sparse matrix with rows and columns labeled
#'
#' @importFrom Matrix readMM
#'
#' @export
#'
Read10X <- function(data.dir = NULL){
  full.data <- list()
  for (i in seq_along(data.dir)) {
    run <- data.dir[i]
    if (! dir.exists(run)){
      stop("Directory provided does not exist")
    }
    if(!grepl("\\/$", run)){
      run <- paste(run, "/", sep = "")
    }
    barcode.loc <- paste0(run, "barcodes.tsv")
    gene.loc <- paste0(run, "genes.tsv")
    matrix.loc <- paste0(run, "matrix.mtx")
    if (!file.exists(barcode.loc)){
      stop("Barcode file missing")
    }
    if (! file.exists(gene.loc)){
      stop("Gene name file missing")
    }
    if (! file.exists(matrix.loc)){
      stop("Expression matrix file missing")
    }
    data <- readMM(file = matrix.loc)
    cell.names <- readLines(barcode.loc)
    gene.names <- readLines(gene.loc)
    if (all(grepl(pattern = "\\-1$", x = cell.names))) {
      cell.names <- as.vector(
        x = as.character(
          x = sapply(
            X = cell.names,
            FUN = extract_field, field = 1,
            delim = "-"
          )
        )
      )
    }
    rownames(x = data) <- make.unique(
      names = as.character(
        x = sapply(
          X = gene.names,
          FUN = extract_field,
          field = 2,
          delim = "\\t"
        )
      )
    )
    if (is.null(x = names(x = data.dir))) {
      if(i < 2){
        colnames(x = data) <- cell.names
      }
      else {
        colnames(x = data) <- paste0(i, "_", cell.names)
      }
    } else {
      colnames(x = data) <- paste0(names(x = data.dir)[i], "_", cell.names)
    }
    full.data <- append(x = full.data, values = data)
  }
  full.data <- do.call(cbind, full.data)
  return(full.data)
}

#' Scale and center the data
#'
#' @param object Seurat object
#' @param genes.use Vector of gene names to scale/center. Default is all genes in object@@data.
#' @param data.use Can optionally pass a matrix of data to scale, default is object@data[genes.use,]
#' @param do.scale Whether to scale the data.
#' @param do.center Whether to center the data.
#' @param scale.max Max value to accept for scaled data. The default is 10. Setting this can help
#' reduce the effects of genes that are only expressed in a very small number of cells.
#'
#' @return Returns a seurat object with object@@scale.data updated with scaled and/or centered data.
#'
#' @export
#'
ScaleData <- function(
  object,
  genes.use = NULL,
  data.use = NULL,
  do.scale = TRUE,
  do.center = TRUE,
  scale.max = 10
) {
  genes.use <- set.ifnull(genes.use,rownames(object@data))
  genes.use <- ainb(genes.use,rownames(object@data))
  data.use <- set.ifnull(data.use,object@data[genes.use, ])
  object@scale.data <- matrix(
    data = NA,
    nrow = length(x = genes.use),
    ncol = ncol(x = object@data)
  )
  #rownames(object@scale.data) <- genes.use
  #colnames(object@scale.data) <- colnames(object@data)
  dimnames(x = object@scale.data) <- dimnames(x = data.use)
  if (do.scale | do.center) {
    bin.size <- 1000
    max.bin <- floor(length(genes.use)/bin.size) + 1
    print("Scaling data matrix")
    pb <- txtProgressBar(min = 0, max = max.bin, style = 3)
    for (i in 1:max.bin) {
      my.inds <- ((bin.size * (i - 1)):(bin.size * i - 1)) + 1
      my.inds <- my.inds[my.inds <= length(x = genes.use)]
      #print(my.inds)
      new.data <- t(
        x = scale(
          x = t(x = as.matrix(x = data.use[genes.use[my.inds], ])),
          center = do.center,
          scale = do.scale
        )
      )
      new.data[new.data > scale.max] <- scale.max
      object@scale.data[genes.use[my.inds], ] <- new.data
      setTxtProgressBar(pb, i)
    }
    close(pb)
  }
  return(object)
}


#' Scale and center the data using C++
#'
#' Note: will give slightly different results than R due to differences in numerical precision
#' between R and C++. This could cause the results of some stochastic processes like tSNE to change.
#'

#' @param object Seurat object
#' @param genes.use Vector of gene names to scale/center. Default is all genes in object@@data.
#' @param data.use Can optionally pass a matrix of data to scale, default is object@data[genes.use,]
#' @param do.scale Whether to scale the data.
#' @param do.center Whether to center the data.
#' @param scale.max Max value to return for scaled data. The default is 10. Setting this can help
#' reduce the effects of genes that are only expressed in a very small number of cells.
#' @param assay.type Assay to scale data for. Default is RNA. Can be changed for multimodal analyses.
#'
#' @return Returns a seurat object with object@@scale.data updated with scaled and/or centered data.
#'
#' @export
#'
FastScaleData <- function(
  object,
  genes.use = NULL,
  data.use = NULL,
  do.scale = TRUE,
  do.center = TRUE,
  scale.max = 10,
  display.progress = TRUE,
  assay.type = "RNA"
) {
  orig.data <- GetAssayData(object = object, assay.type,"data")
  genes.use <- set.ifnull(genes.use, rownames(x = orig.data))
  genes.use <- as.vector(
    x = intersect(
      x = genes.use,
      y = rownames(x = orig.data)
    )
  )
  data.use <- set.ifnull(data.use, orig.data[genes.use, ])
  if (ncol(x = data.use) > 50000) {
    return(ScaleData(
      object = object,
      genes.use = genes.use,
      data.use = data.use,
      do.scale = do.scale,
      do.center = do.center,
      scale.max = scale.max
    ))
  }
  if (class(data.use) == "dgCMatrix" || class(data.use) == "dgTMatrix") {
    data.scale <- FastSparseRowScale(
      mat = data.use,
      scale = do.scale,
      center = do.center,
      scale_max = scale.max,
      display_progress = display.progress
    )
  } else {
    data.use <- as.matrix(x = data.use)
    data.scale <- FastRowScale(
      mat = data.use,
      scale = do.scale,
      center = do.center,
      scale_max = scale.max,
      display_progress = display.progress
    )
  }
  dimnames(x = data.scale) <- dimnames(x = data.use)
  object <- SetAssayData(
    object = object,
    assay.type = assay.type,
    slot = 'scale.data',
    new.data = data.scale
  )
  return(object)
}

#' Normalize raw data
#'
#' Normalize count data per cell and transform to log scale
#'
#' @param data Matrix with the raw count data
#' @param scale.factor Scale the data. Default is 1e4
#' @param display.progress Print progress
#'
#' @return Returns a matrix with the normalize and log transformed data
#'
#' @import Matrix
#'
#' @export
#'
LogNormalize <- function(data, scale.factor = 1e4, display.progress = TRUE) {
  if (class(x = data) == "data.frame") {
    data <- as.matrix(x = data)
  }
  if (class(x = data) != "dgCMatrix") {
    data <- as(object = data, Class = "dgCMatrix")
  }
  # call Rcpp function to normalize
  if (display.progress) {
    print("Performing log-normalization")
  }
  norm.data <- LogNorm(data, scale_factor = scale.factor, display_progress = display.progress)
  colnames(x = norm.data) <- colnames(x = data)
  rownames(x = norm.data) <- rownames(x = data)
  return(norm.data)
}

#' Sample UMI
#'
#' Downsample each cell to a specified number of UMIs. Includes
#' an option to upsample cells below specified UMI as well.
#'
#' @param data Matrix with the raw count data
#' @param max.umi Number of UMIs to sample to
#' @param upsample Upsamples all cells with fewer than max.umi
#' @param progress.bar Display the progress bar
#'
#' @import Matrix
#'
#' @return Matrix with downsampled data
#'
#' @export
#'
SampleUMI <- function(
  data,
  max.umi = 1000,
  upsample = FALSE,
  progress.bar = TRUE
) {
  data <- as(data, "dgCMatrix")
  if (length(x = max.umi) == 1) {
    return(
      RunUMISampling(
        data = data,
        sample_val = max.umi,
        upsample = upsample,
        display_progress = progress.bar
      )
    )
  } else if (length(x = max.umi) != ncol(x = data)) {
    stop("max.umi vector not equal to number of cells")
  }
  return(
    RunUMISamplingPerCell(
      data = data,
      sample_val = max.umi,
      upsample = upsample,
      display_progress = progress.bar
    )
  )
}
