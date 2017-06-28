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
#'    \item{\code{ident}:}{\code{"vector"},  The 'identity class' for each single cell }
#'    \item{\code{data.info}:}{\code{"data.frame"}, Contains information about each cell, starting with # of genes detected (nGene)
#'    the original identity class (orig.ident), user-provided information (through AddMetaData), etc.  }
#'    \item{\code{project.name}:}{\code{"character"}, Name of the project (for record keeping) }
#'    \item{\code{dr}:}{\code{"list"}, List of stored dimensional reductions. Named by technique }
#'    \item{\code{assay}:}{\code{"list"}, List of additional assays for multimodal analysis. Named by technique }
#'    \item{\code{tsne.rot}:}{\code{"data.frame"}, Cell coordinates on the t-SNE map }
#'    \item{\code{mean.var}:}{\code{"data.frame"}, The output of the mean/variability analysis for all genes }
#'    \item{\code{imputed}:}{\code{"data.frame"}, Matrix of imputed gene scores }
#'    \item{\code{final.prob}:}{\code{"data.frame"}, For spatial inference, posterior probability of each cell mapping to each bin }
#'    \item{\code{insitu.matrix}:}{\code{"data.frame"}, For spatial inference, the discretized spatial reference map }
#'    \item{\code{cell.names}:}{\code{"vector"},  Names of all single cells (column names of the expression matrix) }
#'    \item{\code{cluster.tree}:}{\code{"list"},  List where the first element is a phylo object containing the
#'    phylogenetic tree relating different identity classes }
#'    \item{\code{snn.sparse}:}{\code{"dgCMatrix"}, Sparse matrix object representation of the SNN graph }
#'    \item{\code{snn.dense}:}{\code{"matrix"}, Dense matrix object representation of the SNN graph }
#'    \item{\code{snn.k}:}{\code{"numeric"}, k used in the construction of the SNN graph }
#'
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
    ident = "vector",
    dr = "list",
    assay = "list",
    emp.pval = "data.frame",
    kmeans.obj = "list",
    gene.scores = "data.frame",
    drop.coefs = "data.frame",
    wt.matrix = "data.frame",
    drop.wt.matrix = "data.frame",
    trusted.genes = "vector",
    drop.expr = "numeric",
    data.info = "data.frame",
    project.name = "character",
    kmeans.gene = "list",
    kmeans.cell = "list",
    jackStraw.empP = "data.frame",
    jackStraw.fakePC = "data.frame",
    jackStraw.empP.full = "data.frame",
    kmeans.col = "list",
    mean.var = "data.frame",
    imputed = "data.frame",
    mix.probs = "data.frame",
    mix.param = "data.frame",
    final.prob = "data.frame",
    insitu.matrix = "data.frame",
    tsne.rot = "data.frame",
    cell.names = "vector",
    cluster.tree = "list",
    snn.sparse = "dgCMatrix",
    snn.dense = "matrix",
    snn.k="numeric"
  )
)


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
#'
#' @param object Seurat object
#' @param genes.use Vector of gene names to scale/center. Default is all genes in object@@data.
#' @param data.use Can optionally pass a matrix of data to scale, default is object@data[genes.use,]
#' @param do.scale Whether to scale the data.
#' @param do.center Whether to center the data.
#' @param scale.max Max value to return for scaled data. The default is 10. Setting this can help
#' reduce the effects of genes that are only expressed in a very small number of cells.
#' @param assay.type Assay to scale data for. Default is RNA. Can be changed for multimodal analyses.
#' @return Returns a seurat object with object@@scale.data updated with scaled and/or centered data.
#' @export
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


#' Make object sparse
#'
#' Converts stored data matrices to sparse matrices to save space. Converts object@@raw.data and object@@data to sparse matrices.
#' If the snn has been stored as a dense matrix, this will convert it to a sparse matrix, store it in object@@snn.sparse and
#' remove object@@snn.dense.
#'
#' @param object Seurat object
#'
#' @return Returns a seurat object with data converted to sparse matrices.
#'
#' @import Matrix
#'
#' @export
#'
MakeSparse <- function(object) {
  if (class(object@raw.data) == "data.frame") {
    object@raw.data <- as.matrix(x = object@raw.data)
  }
  if (class(object@data) == "data.frame") {
    object@data <- as.matrix(x = object@data)
  }
  if (length(x = object@snn.sparse) == 1 && length(x = object@snn.dense) > 1) {
    if (class(object@snn.dense) == "data.frame") {
      object@snn.dense <- as.matrix(x = object@snn.dense)
    }
    object@snn.sparse <- as(object = object@snn.dense, Class = "dgCMatrix")
    object@snn.dense <- matrix()
  }
  object@raw.data <- as(object = object@raw.data, Class = "dgCMatrix")
  object@data <- as(object = object@data, Class = "dgCMatrix")
  return(object)
}

#' Convert old Seurat object to accomodate new features
#'
#' Adds the object@@dr slot to older objects and moves the stored PCA/ICA analyses to new slot
#'
#' @param object Seurat object
#'
#' @return Returns a Seurat object compatible with latest changes
#'
#' @export
#'
ConvertSeurat <- function(object) {
    if (! (.hasSlot(object = object, name = "assay"))) {
        object@assay=list()
    }
    if ((.hasSlot(object,"dr"))) {
        return(object)
    }
    object@dr <- list()
    pca.x <- matrix()
    pca.x.full <- matrix()
    pca.rotation <- matrix()
    pca.sdev <- numeric()
    pca.misc <- NULL
    if (length(x = object@pca.x) > 0) {
        pca.x <- as.matrix(object@pca.x)
    }
    if (length(x = object@pca.x.full) > 0) {
        pca.x.full <- as.matrix(object@pca.x.full)
    }
    if (length(x = object@pca.rot) > 0) {
        pca.rotation <- as.matrix(object@pca.rot)
    }
    if (length(x = object@pca.obj) > 0) {
        pca.sdev <- object@pca.obj[[1]]$sdev
        if (is.null(x = pca.sdev)) {
            pca.sdev <- object@pca.obj[[1]]$d
        }
        pca.misc <- object@pca.obj[[1]]
    }
    pca.checks <- c(
        vapply(
            X = list(pca.x, pca.x.full, pca.rotation),
            FUN = function(pca) {
                return(length(x = pca) > 1)
            },
            FUN.VALUE = logical(1)
        ),
        length(x = pca.sdev) > 0,
        ! is.null(x = pca.misc)
    )
    # if(length(pca.x) > 1 || length(pca.x.full) > 1 || length(pca.rotation) > 1 || length(pca.sdev) > 0 || !is.null(pca.misc)) {
    if (all(pca.checks)) {
        pca.obj <- new(
            Class = "dim.reduction",
            x = pca.x,
            x.full = pca.x.full,
            rotation = pca.rotation,
            sdev = pca.sdev,
            key = "PC",
            misc = pca.misc
        )
        object@dr$pca <- pca.obj
    }
    ica.x <- matrix()
    ica.rotation <- matrix()
    ica.sdev <- numeric()
    ica.misc <- NULL
    if (length(x = object@ica.x) > 0) {
        ica.x <- as.matrix(object@ica.x)
    }
    if (length(x = object@ica.rot) > 0) {
        ica.rotation <- as.matrix(object@ica.rot)
    }
    if (length(x = object@ica.obj) > 0) {
        ica.sdev <- sqrt(object@ica.obj[[1]]$vafs)
        ica.misc <- object@ica.obj[[1]]
    }
    ica.checks <- c(
        vapply(
            X = list(ica.x, ica.rotation),
            FUN = function(ica) {
                return(length(x = ica) > 1)
            },
            FUN.VALUE = logical(1)
        ),
        length(x = ica.sdev) > 0,
        ! is.null(x = ica.misc)
    )
  # if(length(ica.x) > 1 || length(ica.rotation) > 1 || length(ica.sdev) > 0  || !is.null(ica.misc)) {
    if (all(ica.checks)) {
        ica.obj <- new(
            "dim.reduction",
            x = ica.x,
            rotation = ica.rotation,
            sdev = ica.sdev,
            key = "IC",
            misc = "ica.misc"
        )
        object@dr$ica <- ica.obj
    }
    tsne.rotation <- matrix()
    if (length(x = object@tsne.rot) > 0) {
        tsne.rotation <- as.matrix(object@tsne.rot)
    }
    if (length(x = tsne.rotation) > 1) {
        tsne.obj <- new(
            "dim.reduction",
            rotation = tsne.rotation,
            key = "tSNE_"
        )
        object@dr$tsne <- tsne.obj
    }
    return(object)
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

#' Add samples into existing Seurat object.
#'
#' @param object Seurat object
#' @param project Project name (string)
#' @param new.data Data matrix for samples to be added
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
#' @param add.cell.id String to be appended to the names of all cells in new.data. E.g. if add.cell.id = "rep1",
#' "cell1" becomes "cell1.rep1"
#'
#' @import Matrix
#' @importFrom dplyr full_join
#'
#' @export
#'
AddSamples <- function(
  object,
  new.data,
  project = NULL,
  min.cells = 3,
  min.genes = 1000,
  is.expr = 0,
  do.logNormalize = TRUE,
  total.expr = 1e4,
  do.scale=TRUE,
  do.center = TRUE,
  names.field = 1,
  names.delim = "_",
  meta.data = NULL,
  save.raw = TRUE,
  add.cell.id = NULL
) {
  if (length(x = object@raw.data) < 2) {
    stop("Object provided has an empty raw.data slot. Adding/Merging performed on raw count data.")
  }
  if (! missing(x = add.cell.id)) {
    colnames(x= new.data) <- paste(colnames(x = new.data), add.cell.id, sep = ".")
  }
  if (any(colnames(x = new.data) %in% object@cell.names)) {
    warning("Duplicate cell names, enforcing uniqueness via make.unique()")
    new.data.names <- as.list(
      x = make.unique(
        names = c(
          colnames(x = object@raw.data),
          colnames(x = new.data)
        )
      )[(ncol(x = object@raw.data) + 1):(ncol(x = object@raw.data) + ncol(x = new.data))]
    )
    names(x = new.data.names) <- colnames(x = new.data)
    colnames(x = new.data) <- new.data.names
    if (! is.null(x = meta.data)){
      rownames(x = meta.data) <- unlist(
        x = unname(
          obj = new.data.names[rownames(x = meta.data)]
        )
      )
    }
  }
  combined.data <- RowMergeSparseMatrices(
    mat1 = object@raw.data[, object@cell.names],
    mat2 = new.data
  )
  new.object <- new(Class = "seurat", raw.data = combined.data)
  if (is.null(x = meta.data)) {
    filler <- matrix(NA, nrow = ncol(new.data), ncol = ncol(object@data.info))
    rownames(filler) <- colnames(new.data)
    colnames(filler) <- colnames(object@data.info)
    filler <- as.data.frame(filler)
    combined.meta.data <- rbind(object@data.info, filler)
  } else {
    combined.meta.data <- suppressMessages(
      suppressWarnings(
        full_join(x = object@data.info, y = meta.data)
      )
    )
  }
  project <- set.ifnull(project, object@project.name)
  new.object <- Setup(
    new.object,
    project,
    min.cells = min.cells,
    min.genes = min.genes,
    is.expr = is.expr,
    do.logNormalize = do.logNormalize,
    total.expr = total.expr,
    do.scale = do.scale,
    do.center = do.center,
    names.field = names.field,
    names.delim = names.delim,
    save.raw = save.raw
  )
  new.object@data.info <- combined.meta.data[new.object@cell.names,]
  return(new.object)
}


#' Merge Seurat Objects
#'
#' Merge two Seurat objects
#'
#' @param object1 First Seurat object to merge
#' @param object2 Second Seurat object to merge
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
#' @param add.cell.id1 String to be appended to the names of all cells in object1
#' @param add.cell.id2 String to be appended to the names of all cells in object2
#'
#' @return Merged Seurat object
#'
#' @import Matrix
#' @importFrom dplyr full_join filter
#'
#' @export
#'
MergeSeurat <- function(
    object1,
    object2,
    project = NULL,
    min.cells = 0,
    min.genes = 0,
    is.expr = 0,
    do.logNormalize = TRUE,
    total.expr = 1e4,
    do.scale = TRUE,
    do.center = TRUE,
    names.field = 1,
    names.delim = "_",
    save.raw = TRUE,
    add.cell.id1 = NULL,
    add.cell.id2 = NULL
) {
  if (length(x = object1@raw.data) < 2) {
    stop("First object provided has an empty raw.data slot. Adding/Merging performed on raw count data.")
  }
  if (length(x = object2@raw.data) < 2) {
    stop("Second object provided has an empty raw.data slot. Adding/Merging performed on raw count data.")
  }
  if (! missing(add.cell.id1)) {
    object1@cell.names <- paste(object1@cell.names, add.cell.id1, sep = ".")
    colnames(x = object1@raw.data) <- paste(
      colnames(x = object1@raw.data),
      add.cell.id1,
      sep = "."
    )
    rownames(x = object1@data.info) <- paste(
      rownames(x = object1@data.info),
      add.cell.id1,
      sep = "."
    )
  }
  if (! missing(add.cell.id2)) {
    object2@cell.names <- paste(object2@cell.names, add.cell.id2, sep = ".")
    colnames(x = object2@raw.data) <- paste(
      colnames(x = object2@raw.data),
      add.cell.id2,
      sep = "."
    )
    rownames(x = object2@data.info) <- paste(
      rownames(x = object2@data.info),
      add.cell.id2,
      sep = "."
    )
  }
  if (any(object1@cell.names %in% object2@cell.names)) {
    warning("Duplicate cell names, enforcing uniqueness via make.unique()")
    object2.names <- as.list(
      x = make.unique(
        names = c(
          colnames(x = object1@raw.data),
          colnames(x = object2@raw.data)
        )
      )[(ncol(x = object1@raw.data) + 1):(ncol(x = object1@raw.data) + ncol(x = object2@raw.data))]
    )
    names(x = object2.names) <- colnames(x = object2@raw.data)
    colnames(x = object2@raw.data) <- object2.names
    object2@cell.names <- unlist(
      x = unname(
        obj = object2.names[object2@cell.names]
      )
    )
    rownames(x = object2@data.info) <- unlist(
      x = unname(
        obj = object2.names[rownames(x = object2@data.info)]
      )
    )
  }
  merged.raw.data <- RowMergeSparseMatrices(
    mat1 = object1@raw.data[,object1@cell.names],
    mat2 = object2@raw.data[,object2@cell.names]
  )
  object1@data.info <- object1@data.info[object1@cell.names, ]
  object2@data.info <- object2@data.info[object2@cell.names, ]
  new.object <- new(Class = "seurat", raw.data = merged.raw.data)
  project <- set.ifnull(project, object1@project.name)
  object1@data.info$cell.name <- rownames(x = object1@data.info)
  object2@data.info$cell.name <- rownames(x = object2@data.info)
  merged.meta.data <- suppressMessages(
    suppressWarnings(
      full_join(x = object1@data.info, y = object2@data.info)
    )
  )
  merged.object <- Setup(
    new.object,
    project = project,
    min.cells = min.cells,
    min.genes = min.genes,
    is.expr = is.expr,
    do.logNormalize = do.logNormalize,
    total.expr = total.expr,
    do.scale = do.scale,
    do.center = do.center,
    names.field = names.field,
    names.delim = names.delim,
    save.raw = save.raw
  )
  merged.meta.data %>% filter(
    cell.name %in% merged.object@cell.names
  ) -> merged.meta.data
  rownames(x= merged.meta.data) <- merged.object@cell.names
  merged.meta.data$cell.name <- NULL
  merged.object@data.info <- merged.meta.data
  return(merged.object)
}

# Internal function for merging two matrices by rowname
RowMergeSparseMatrices <- function(mat1, mat2){
  if (inherits(mat1, "data.frame")) {
    mat1 = as.matrix(mat1)
  }
  if (inherits(mat2, "data.frame")) {
    mat2 = as.matrix(mat2)
  }
  mat1 <- as(mat1, "RsparseMatrix")
  mat2 <- as(mat2, "RsparseMatrix")
  mat1.names <- rownames(mat1)
  mat2.names <- rownames(mat2)
  all.names <- union(mat1.names, mat2.names)
  new.mat <- RowMergeMatrices(mat1 = mat1, mat2 = mat2, mat1_rownames = mat1.names, mat2_rownames = mat2.names, all_rownames = all.names)
  rownames(new.mat) <- make.unique(all.names)
  #colnames(mat2) <- sprintf('%s_2', colnames(mat2))
  colnames(new.mat) <- make.unique(c(colnames(mat1), colnames(mat2)))
  return(new.mat)
}


#' Classify New Data
#'
#' Classify new data based on the cluster information of the provided object.
#' Random Forests are used as the basis of the classification.
#'
#' @param object Seurat object on which to train the classifier
#' @param classifier Random Forest classifier from BuildRFClassifier. If not provided,
#' it will be built from the training data provided.
#' @param training.genes Vector of genes to build the classifier on
#' @param training.classes Vector of classes to build the classifier on
#' @param new.data New data to classify
#' @param ... additional parameters passed to ranger
#'
#' @return Vector of cluster ids
#'
#' @import Matrix
#' @importFrom ranger ranger
#' @importFrom plyr mapvalues
#'
#' @export
#'
ClassifyCells <- function(
  object,
  classifier,
  training.genes = NULL,
  training.classes = NULL,
  new.data = NULL,
  ...
) {
  # build the classifier
  if (missing(classifier)){
    classifier <- BuildRFClassifier(
      object = object,
      training.genes = training.genes,
      training.classes = training.classes,
      ...
    )
  }
  # run the classifier on the new data
  features <- classifier$forest$independent.variable.names
  genes.to.add <- setdiff(x = features, y = rownames(x = new.data))
  data.to.add <- matrix(
    data = 0,
    nrow = length(x = genes.to.add),
    ncol = ncol(x = new.data)
  )
  rownames(x = data.to.add) <- genes.to.add
  new.data <- rbind(new.data, data.to.add)
  new.data <- new.data[features, ]
  new.data <- as.matrix(x = t(x = new.data))
  print("Running Classifier ...")
  prediction <- predict(classifier, new.data)
  new.classes <- prediction$predictions
  return(new.classes)
}

#' Build Random Forest Classifier
#'
#' Train the random forest classifier
#'
#'
#' @param object Seurat object on which to train the classifier
#' @param training.genes Vector of genes to build the classifier on
#' @param training.classes Vector of classes to build the classifier on
#' @param verbose Additional progress print statements
#' @param ... additional parameters passed to ranger
#'
#' @return Returns the random forest classifier
#'
#' @import Matrix
#' @importFrom ranger ranger
#' @importFrom plyr mapvalues
#'
#' @export
#'
BuildRFClassifier <- function(
  object,
  training.genes = NULL,
  training.classes = NULL,
  verbose = TRUE,
  ...
) {
  training.classes <- as.vector(x = training.classes)
  training.genes <- set.ifnull(training.genes, rownames(x = object@data))
  training.data <- as.data.frame(
    x = as.matrix(
      x = t(
        x = object@data[training.genes, ]
      )
    )
  )
  training.data$class <- factor(x = training.classes)
  if (verbose) {
    print("Training Classifier ...")
  }
  classifier <- ranger(
    data = training.data,
    dependent.variable.name = "class",
    classification = TRUE,
    write.forest = TRUE,
    ...
  )
  return(classifier)
}

#' Highlight classification results
#'
#' This function is useful to view where proportionally the clusters returned from
#' classification map to the clusters present in the given object. Utilizes the FeaturePlot()
#' function to color clusters in object.
#'
#' @param object Seurat object on which the classifier was trained and
#' onto which the classification results will be highlighted
#' @param clusters vector of cluster ids (output of ClassifyCells)
#' @param ... additional parameters to pass to FeaturePlot()
#'
#' @return Returns a feature plot with clusters highlighted by proportion of cells
#' mapping to that cluster
#'
#' @export
#'
VizClassification <- function(object, clusters, ...) {
  cluster.dist <- prop.table(x = table(out)) # What is out?
  object@data.info$Classification <- numeric(nrow(x = object@data.info))
  for (cluster in 1:length(x = cluster.dist)) {
    cells.to.highlight <- WhichCells(object, names(cluster.dist[cluster]))
    if (length(x = cells.to.highlight) > 0) {
      object@data.info[cells.to.highlight, ]$Classification <- cluster.dist[cluster]
    }
  }
  if (any(grepl(pattern = "cols.use", x = deparse(match.call())))) {
    return(FeaturePlot(object, "Classification", ...))
  }
  cols.use = c("#f6f6f6", "black")
  return(FeaturePlot(object, "Classification", cols.use = cols.use, ...))
}

#   Documentation...
####################
calc.drop.prob <- function(x, a, b) {
  return(exp(x = a + b * x) / (1 + exp(x = a + b * x)))
}

#' Find all markers for a node
#'
#' This function finds markers for all splits at or below the specified node
#'
#' @param object Seurat object. Must have object@@cluster.tree slot filled. Use BuildClusterTree() if not.
#' @param node Node from which to start identifying split markers, default is top node.
#' @param genes.use Genes to test. Default is to use all genes
#' @param thresh.use Limit testing to genes which show, on average, at least
#' X-fold difference (log-scale) between the two groups of cells.
#' @param test.use Denotes which test to use. Seurat currently implements
#' "bimod" (likelihood-ratio test for single cell gene expression, McDavid et
#' al., Bioinformatics, 2011, default), "roc" (standard AUC classifier), "t"
#' (Students t-test), and "tobit" (Tobit-test for differential gene expression,
#' as in Trapnell et al., Nature Biotech, 2014), 'poisson', and 'negbinom'.
#' The latter two options should only be used on UMI datasets, and assume an underlying
#' poisson or negative-binomial distribution.
#' @param min.pct - only test genes that are detected in a minimum fraction of min.pct cells
#' in either of the two populations. Meant to speed up the function by not testing genes that are very infrequently expression
#' @param min.diff.pct - only test genes that show a minimum difference in the fraction of detection between the two groups. Set to -Inf by default
#' @param only.pos Only return positive markers (FALSE by default)
#' @param print.bar Print a progress bar once expression testing begins (uses pbapply to do this)
#' @param max.cells.per.ident Down sample each identity class to a max number. Default is no downsampling.
#' @param random.seed Random seed for downsampling
#' @param return.thresh Only return markers that have a p-value < return.thresh, or a power > return.thresh (if the test is ROC)
#' @param min.cells Minimum number of cells expressing the gene in at least one of the two groups
#'
#' @return Returns a dataframe with a ranked list of putative markers for each node and associated statistics
#'
#' @importFrom ape drop.tip
#'
#' @export
#'
FindAllMarkersNode <- function(
  object,
  node = NULL,
  genes.use = NULL,
  thresh.use = 0.25,
  test.use = "bimod",
  min.pct = 0.1,
  min.diff.pct = 0.05,
  print.bar = TRUE,
  only.pos = FALSE,
  max.cells.per.ident = Inf,
  return.thresh = 1e-2,
  do.print = FALSE,
  random.seed = 1,
  min.cells = 3
) {
  genes.use <- set.ifnull(genes.use, rownames(object@data))
  node <- set.ifnull(node, tree$edge[1,1])
  ident.use <- object@ident
  tree.use <- object@cluster.tree[[1]]
  descendants <- DFT(tree = tree.use, node = node, path = NULL, include.children = TRUE)
  all.children <- sort(x = tree.use$edge[,2][!tree.use$edge[,2] %in% tree.use$edge[,1]])
  descendants <- suppressMessages(mapvalues(descendants, from = all.children, to = tree.use$tip.label))
  drop.children <- setdiff(tree.use$tip.label, descendants)
  keep.children <- setdiff(tree.use$tip.label, drop.children)
  orig.nodes <- c(node, as.numeric(setdiff(descendants, keep.children)))
  tree.use <- drop.tip(tree.use, drop.children)
  new.nodes <- unique(tree.use$edge[,1])
  if ((test.use == 'roc') && (return.thresh==1e-2)) {
      return.thresh <- 0.7
  }
  genes.de <- list()
  for (i in ((tree.use$Nnode+2):max(tree.use$edge))) {
    genes.de[[i]] <- FindMarkersNode(
      object = object,
      node = i,
      tree.use = tree.use,
      genes.use = genes.use,
      thresh.use = thresh.use,
      test.use = test.use,
      min.pct = min.pct,
      min.diff.pct = min.diff.pct,
      print.bar = print.bar,
      only.pos = only.pos,
      max.cells.per.ident = max.cells.per.ident,
      random.seed = random.seed,
      min.cells = min.cells
    )
    if (do.print) {
      print(paste("Calculating node", i))
    }
  }
  gde.all <- data.frame()
  for (i in ((tree.use$Nnode+2):max(tree.use$edge))) {
    gde <- genes.de[[i]]
    if (is.null(x = gde)) {
      next
    }
    if (nrow(x = gde) > 0) {
      if (test.use == 'roc') {
        gde <- subset(
          x = gde,
          subset = (myAUC > return.thresh | myAUC < (1 - return.thresh))
        )
      }
      if ( (test.use == 'bimod') || (test.use == 't')) {
        gde <- gde[order(gde$p_val,-gde$avg_diff), ]
        gde <- subset(data = gde, subset = p_val < return.thresh)
      }
      if (nrow(x = gde) > 0) {
        gde$cluster <- i
        gde$gene <- rownames(x = gde)
      }
      if (nrow(x = gde) > 0) {
        gde.all <- rbind(gde.all,gde)
      }
    }
  }
  gde.all$cluster <- mapvalues(
    x = gde.all$cluster,
    from = new.nodes,
    to = orig.nodes
  )
  return(gde.all)
}

#   Documentation...
####################
weighted.euclidean <- function(x, y, w) {
  v.dist <- sum(sqrt(x = w * (x - y) ^ 2))
  return(v.dist)
}

#   Documentation...
####################
#from Jean Fan - thanks!!
custom.dist <- function(my.mat, my.function, ...) {
  n <- ncol(x = my.mat)
  mat <- matrix(data = 0, ncol = n, nrow = n)
  colnames(x = mat) <- rownames(x = mat) <- colnames(x = my.mat)
  for (i in 1:nrow(x = mat)) {
    for(j in 1:ncol(x = mat)) {
      mat[i,j] <- my.function(my.mat[, i], my.mat[, j], ...)
    }
  }
  return(as.dist(mat))
}

#' Phylogenetic Analysis of Identity Classes
#'
#' Constructs a phylogenetic tree relating the 'average' cell from each
#' identity class. Tree is estimated based on a distance matrix constructed in
#' either gene expression space or PCA space.
#'
#' Note that the tree is calculated for an 'average' cell, so gene expression
#' or PC scores are averaged across all cells in an identity class before the
#' tree is constructed.
#'
#' @param object Seurat object
#' @param genes.use Genes to use for the analysis. Default is the set of
#' variable genes (object@@var.genes). Assumes pcs.use=NULL (tree calculated in
#' gene expression space)
#' @param pcs.use If set, tree is calculated in PCA space, using the
#' eigenvalue-weighted euclidean distance across these PC scores.
#' @param SNN.use If SNN is passed, build tree based on SNN graph connectivity between clusters
#' @param do.plot Plot the resulting phylogenetic tree
#' @param do.reorder Re-order identity classes (factor ordering), according to
#' position on the tree. This groups similar classes together which can be
#' helpful, for example, when drawing violin plots.
#' @param reorder.numeric Re-order identity classes according to position on
#' the tree, assigning a numeric value ('1' is the leftmost node)
#'
#' @return A Seurat object where the cluster tree is stored in
#' object@@cluster.tree[[1]]
#'
#' @importFrom ape as.phylo
#'
#' @export
#'
BuildClusterTree <- function(
    object,
    genes.use = NULL,
    pcs.use = NULL,
    SNN.use = NULL,
    do.plot = TRUE,
    do.reorder = FALSE,
    reorder.numeric = FALSE
) {
  genes.use <- set.ifnull(genes.use,object@var.genes)
  ident.names <- as.character(x = unique(x = object@ident))
  if (! is.null(x = genes.use)) {
    genes.use <- ainb(genes.use, rownames(x = object@data))
    data.avg <- AverageExpression(object = object, genes.use = genes.use)
    data.dist <- dist(t(x = data.avg[genes.use, ]))
  }
  if (! is.null(x = pcs.use)) {
    data.pca <- AveragePCA(object = object)
    data.use <- data.pca[pcs.use,]
    if (is.null(x = object@pca.obj[[1]]$sdev)) {
      data.eigenval <- (object@pca.obj[[1]]$d) ^ 2
    } else {
      data.eigenval <- (object@pca.obj[[1]]$sdev) ^ 2
    }
    data.weights <- (data.eigenval / sum(data.eigenval))[pcs.use]
    data.weights <- data.weights / sum(data.weights)
    data.dist <- custom.dist(
      data.pca[pcs.use, ],
      weighted.euclidean,
      data.weights
    )
  }
  if (! is.null(x = SNN.use)) {
      num.clusters <- length(x = ident.names)
      data.dist = matrix(data = 0, nrow = num.clusters, ncol = num.clusters)
      for (i in 1:(num.clusters - 1)) {
          for (j in (i + 1):num.clusters) {
              subSNN <- SNN.use[
                  match(
                    x = WhichCells(object = object, ident = i),
                    table = colnames(x = SNN.use)
                  ), # Row
                  match(
                    x = WhichCells(object = object, ident = j),
                    table = rownames(x = SNN.use)
                  ) # Column
              ]
              d <- mean(subSNN)
              if (is.na(x = d)) {
                  data.dist[i, j] <- 0
              } else {
                  data.dist[i, j] = d
              }
          }
      }
      diag(x = data.dist) <- 1
      data.dist <- dist(data.dist)
  }
  data.tree <- as.phylo(x = hclust(d = data.dist))
  object@cluster.tree[[1]] <- data.tree
  if (do.reorder) {
    old.ident.order <- sort(x = unique(x = object@ident))
    data.tree <- object@cluster.tree[[1]]
    all.desc <- getDescendants(tree = data.tree, node = (data.tree$Nnode + 2))
    all.desc <- old.ident.order[all.desc[all.desc <= (data.tree$Nnode + 1)]]
    object@ident <- factor(x = object@ident, levels = all.desc, ordered = TRUE)
    if (reorder.numeric) {
        object <- SetIdent(
          object = object,
          cells.use = object@cell.names,
          ident.use = as.integer(x = object@ident)
        )
        object@data.info[object@cell.names, "tree.ident"] <- as.integer(x = object@ident)
    }
    object <- BuildClusterTree(
        object = object,
        genes.use = genes.use,
        pcs.use = pcs.use,
        do.plot = FALSE,
        do.reorder = FALSE
    )
  }
  if (do.plot) {
    PlotClusterTree(object)
  }
  return(object)
}

#' Plot phylogenetic tree
#'
#' Plots previously computed phylogenetic tree (from BuildClusterTree)
#'
#' @param object Seurat object
#' @param \dots Additional arguments for plotting the phylogeny
#'
#' @return Plots dendogram (must be precomputed using BuildClusterTree), returns no value
#'
#' @importFrom ape plot.phylo
#' @importFrom ape nodelabels
#'
#' @export
#'
PlotClusterTree <- function(object, ...) {
  if (length(x = object@cluster.tree) == 0) {
    stop("Phylogenetic tree does not exist, build using BuildClusterTree")
  }
  data.tree <- object@cluster.tree[[1]]
  plot.phylo(x = data.tree, direction = "downwards", ...)
  nodelabels()
}

#' Visualize expression/dropout curve
#'
#' Plot the probability of detection vs average expression of a gene.
#'
#' Assumes that this 'noise' model has been precomputed with CalcNoiseModels
#'
#' @param object Seurat object
#' @param cell.ids Cells to use
#' @param col.use Color code or name
#' @param lwd.use Line width for curve
#' @param do.new Create a new plot (default) or add to existing
#' @param x.lim Maximum value for X axis
#' @param \dots Additional arguments to pass to lines function
#' @return Returns no value, displays a plot
#'
#' @export
#'
PlotNoiseModel <- function(
  object,
  cell.ids = c(1, 2),
  col.use = 'black',
  lwd.use = 2,
  do.new = TRUE,
  x.lim = 10,
  ...
) {
  cell.coefs <- object@drop.coefs[cell.ids,]
  if (do.new) {
    plot(
      x = 1,
      y = 1,
      pch = 16,
      type = 'n',
      xlab = 'Average expression',
      ylab = 'Probability of detection',
      xlim = c(0, x.lim),
      ylim =c(0, 1)
    )
  }
  unlist(
    x = lapply(
      X = 1:length(x = cell.ids),
      FUN = function(y) {
        x.vals <- seq(from = 0, to = x.lim, by = 0.05)
        y.vals <- unlist(x = lapply(
          X = x.vals,
          FUN = calc.drop.prob,
          a = cell.coefs[y, 1],
          b = cell.coefs[y, 2])
        )
        lines(
          x = x.vals,
          y = y.vals,
          lwd = lwd.use,
          col = col.use,
          ...
        )
      }
    )
  )
}

#' Regress out technical effects and cell cycle
#'
#' Remove unwanted effects from scale.data
#'
#' @param object Seurat object
#' @param latent.vars effects to regress out
#' @param genes.regress gene to run regression for (default is all genes)
#' @param model.use Use a linear model or generalized linear model (poisson, negative binomial) for the regression. Options are 'linear' (default), 'poisson', and 'negbinom'
#' @param use.umi Regress on UMI count data. Default is FALSE for linear modeling, but automatically set to TRUE if model.use is 'negbinom' or 'poisson'
#' @param do.scale Whether to scale the data.
#' @param do.center Whether to center the data.
#' @param scale.max Max value to accept for scaled data. The default is 10. Setting this can help
#' reduce the effects of genes that are only expressed in a very small number of cells.
#'
#' @return Returns Seurat object with the scale.data (object@scale.data) genes returning the residuals from the regression model
#'
#' @import Matrix
#'
#' @export
#'
RegressOut <- function(
  object,
  latent.vars,
  genes.regress = NULL,
  model.use = 'linear',
  use.umi = FALSE,
  do.scale = TRUE,
  do.center = TRUE,
  scale.max = 10
) {
  possible.models <- c("linear", "poisson", "negbinom")
  if (! model.use %in% possible.models){
    stop(
      paste0(
        model.use,
        " is not a valid model. Please use one the following: ",
        paste0(possible.models, collapse = ", "),
        "."
      )
    )
  }
  genes.regress <- set.ifnull(genes.regress,rownames(object@data))
  genes.regress <- ainb(genes.regress,rownames(object@data))
  latent.data <- FetchData(object = object, vars.all = latent.vars)
  bin.size <- 100
  if (model.use == 'negbinom') {
    bin.size <- 5
  }
  bin.ind <- ceiling(x = 1:length(x = genes.regress) / bin.size)
  max.bin <- max(bin.ind)
  print(paste("Regressing out", latent.vars))
  pb <- txtProgressBar(min = 0, max = max.bin, style = 3)
  data.resid <- c()
  data.use <- object@data[genes.regress, , drop = FALSE];
  if (model.use != "linear") {
    use.umi <- TRUE
  }
  if (use.umi) {
    data.use <- object@raw.data[genes.regress, object@cell.names, drop = FALSE]
  }
  for (i in 1:max.bin) {
    genes.bin.regress <- rownames(x = data.use)[bin.ind == i]
    gene.expr <- as.matrix(x = data.use[genes.bin.regress, , drop = FALSE])
    new.data <- do.call(
      rbind,
      lapply(
        X = genes.bin.regress,
        FUN = function(x) {
          regression.mat <- cbind(latent.data, gene.expr[x,])
          colnames(x = regression.mat) <- c(colnames(x = latent.data), "GENE")
          fmla <- as.formula(
            object = paste0(
              "GENE ",
              " ~ ",
              paste(latent.vars, collapse = "+")
            )
          )
          if (model.use == 'linear') {
            return(lm(formula = fmla, data = regression.mat)$residuals)
          }
          if (model.use == 'poisson') {
            return(residuals(
              object = glm(
                formula = fmla,
                data = regression.mat,
                family = "poisson"
              ),
              type='pearson'
            ))
          }
          if (model.use == 'negbinom') {
            return(nb.residuals(
              formula = fmla,
              regression.mat = regression.mat,
              gene = x
            ))
          }
        }
      )
    )
    if (i == 1) {
      data.resid=new.data
    }
    if (i > 1) {
      data.resid=rbind(data.resid,new.data)
    }
    setTxtProgressBar(pb, i)
  }
  close(pb)
  rownames(x = data.resid) <- genes.regress
  if (use.umi) {
    data.resid <- log1p(
      x = sweep(
        x = data.resid,
        MARGIN = 1,
        STATS = apply(X = data.resid, MARGIN = 1, FUN = min),
        FUN = "-"
      )
    )
  }
  object@scale.data <- data.resid
  if (do.scale) {
    if (use.umi && base::missing(x = scale.max)) {
      scale.max <- 50
    }
    object <- ScaleData(
      object = object,
      genes.use = rownames(x = data.resid),
      data.use = data.resid,
      do.center = do.center,
      do.scale = do.scale,
      scale.max = scale.max
    )
  }
  object@scale.data[is.na(x = object@scale.data)] <- 0
  return(object)
}

#' Regress out technical effects and cell cycle using regularized Negative Binomial regression
#'
#' Remove unwanted effects from umi data and set scale.data to Pearson residuals
#' Uses mclapply; you can set the number of cores it will use to n with command options(mc.cores = n)
#'
#' @param object Seurat object
#' @param latent.vars effects to regress out
#' @param genes.regress gene to run regression for (default is all genes)
#' @param pr.clip.range numeric of length two specifying the min and max values the results will be clipped to
#'
#' @return Returns Seurat object with the scale.data (object@scale.data) genes returning the residuals from the regression model
#'
#' @import Matrix
#' @importFrom MASS theta.ml negative.binomial
#' @import parallel
#'
#' @export
#'
RegressOutNBreg <- function(
  object,
  latent.vars,
  genes.regress = NULL,
  pr.clip.range = c(-30, 30),
  min.theta = 0.01
) {
  genes.regress <- set.ifnull(genes.regress, rownames(x = object@data))
  genes.regress <- ainb(genes.regress,rownames(x = object@data))
  cm <- object@raw.data[genes.regress, colnames(x = object@data), drop = FALSE]
  latent.data <- FetchData(object = object, vars.all = latent.vars)
  cat(sprintf('Regressing out %s for %d genes\n', paste(latent.vars), length(x = genes.regress)))
  theta.fit <- theta.reg(cm = cm, latent.data = latent.data, min.theta = 0.01, bin.size = 128)
  print('Second run NB regression with fixed theta')
  bin.size <- 128
  bin.ind <- ceiling(1:length(genes.regress)/bin.size)
  max.bin <- max(bin.ind)
  pb <- txtProgressBar(min = 0, max = max.bin, style = 3)
  pr <- c()
  for (i in 1:max.bin) {
    genes.bin.regress <- genes.regress[bin.ind == i]
    bin.pr.lst <- parallel::mclapply(
      X = genes.bin.regress,
      FUN = function(j) {
        fit <- 0
        try(
          expr = fit <- glm(
            cm[j, ] ~ .,
            data = latent.data,
            family = MASS::negative.binomial(theta = theta.fit[j])
          ),
          silent=TRUE
        )
        if (class(fit)[1] == 'numeric') {
          message(
            sprintf(
              'glm and family=negative.binomial(theta=%f) failed for gene %s; falling back to scale(log10(y+1))',
              theta.fit[j],
              j
            )
          )
            res <- scale(log10(cm[j, ] + 1))[, 1]
        } else {
            res <- residuals(fit, type = 'pearson')
        }
        return(res)
      }
    )
    pr <- rbind(pr, do.call(rbind, bin.pr.lst))
    setTxtProgressBar(pb, i)
  }
  close(pb)
  dimnames(x = pr) <- dimnames(x = cm)
  pr[pr < pr.clip.range[1]] <- pr.clip.range[1]
  pr[pr > pr.clip.range[2]] <- pr.clip.range[2]
  object@scale.data <- pr
  return(object)
}

#' Regress out technical effects and cell cycle using regularized Negative Binomial regression
#'
#' Remove unwanted effects from umi data and set scale.data to Pearson residuals
#' Uses mclapply; you can set the number of cores it will use to n with command options(mc.cores = n)
#'
#' @param object Seurat object
#' @param latent.vars effects to regress out
#' @param genes.regress gene to run regression for (default is all genes)
#' @param pr.clip.range numeric of length two specifying the min and max values the results will be clipped to
#'
#' @return Returns Seurat object with the scale.data (object@scale.data) genes returning the residuals from the regression model
#'
#' @import Matrix
#' @importFrom MASS theta.ml negative.binomial
#' @import parallel
#'
#' @export

RegressOutNBreg <- function(
  object,
  latent.vars,
  genes.regress = NULL,
  pr.clip.range = c(-30, 30),
  min.theta = 0.01
) {
  genes.regress <- set.ifnull(genes.regress, rownames(x = object@data))
  genes.regress <- ainb(genes.regress, rownames(x = object@data))
  cm <- object@raw.data[genes.regress, colnames(x = object@data), drop=FALSE]
  latent.data <- FetchData(boject = object, vars.all = latent.vars)
  bin.size <- 128
  bin.ind <- ceiling(x = 1:length(x = genes.regress) / bin.size)
  max.bin <- max(bin.ind)
  print(paste("Regressing out", latent.vars))
  print('First run Poisson regression (to get initial mean), and estimate theta per gene')
  pb <- txtProgressBar(min = 0, max = max.bin, style = 3)
  theta.estimate <- c()
  for (i in 1:max.bin) {
    genes.bin.regress <- genes.regress[bin.ind == i]
    bin.theta.estimate <- unlist(
      parallel::mclapply(
        X = genes.bin.regress,
        FUN = function(j) {
          as.numeric(
            x = MASS::theta.ml(
              cm[j, ],
              glm(cm[j, ] ~ ., data = latent.data, family=poisson)$fitted
            )
          )
        }
      ),
      use.names = FALSE
    )
    theta.estimate <- c(theta.estimate, bin.theta.estimate)
    setTxtProgressBar(pb, i)
  }
  close(pb)
  UMI.mean <- apply(X = cm, MARGIN = 1, FUN = mean)
  var.estimate <- UMI.mean + (UMI.mean ^ 2) / theta.estimate
  fit <- loess(log10(var.estimate) ~ log10(UMI.mean), span = 0.33)
  theta.fit <- (UMI.mean ^ 2) / (10 ^ fit$fitted - UMI.mean)
  names(x = theta.fit) <- genes.regress
  to.fix <- theta.fit <= min.theta | is.infinite(x = theta.fit)
  if (any(to.fix)) {
    cat(
      'Fitted theta below',
      min.theta,
      'for',
      sum(to.fix),
      'genes, setting them to',
      min.theta,
      '\n'
    )
    theta.fit[to.fix] <- min.theta
  }
  print('Second run NB regression with fixed theta')
  pb <- txtProgressBar(min = 0, max = max.bin, style = 3)
  pr <- c()
  for(i in 1:max.bin) {
    genes.bin.regress <- genes.regress[bin.ind == i]
    bin.pr.lst <- parallel::mclapply(
      X = genes.bin.regress,
      FUN = function(j) {
        fit <- 0
        try(
          fit <- glm(
            cm[j, ] ~ .,
            data = latent.data,
            family=MASS::negative.binomial(theta = theta.fit[j])
          ),
          silent=TRUE
        )
        if (class(fit)[1] == 'numeric') {
          message(
            sprintf(
              'glm and family=negative.binomial(theta=%f) failed for gene %s; falling back to scale(log10(y+1))',
              theta.fit[j],
              j
            )
          )
          res <- scale(x = log10(cm[j, ] + 1))[, 1]
        } else {
          res <- residuals(object = fit, type='pearson')
        }
        return(res)
      }
    )
    pr <- rbind(pr, do.call(rbind, bin.pr.lst))
    setTxtProgressBar(pb, i)
  }
  close(pb)
  dimnames(pr) <- dimnames(cm)
  pr[pr < pr.clip.range[1]] <- pr.clip.range[1]
  pr[pr > pr.clip.range[2]] <- pr.clip.range[2]
  object@scale.data <- r
  return(object)
}

#' Return a subset of the Seurat object
#'
#' Creates a Seurat object containing only a subset of the cells in the
#' original object. Takes either a list of cells to use as a subset, or a
#' parameter (for example, a gene), to subset on.
#'
#' @param object Seurat object
#' @param cells.use A vector of cell names to use as a subset. If NULL
#' (default), then this list will be computed based on the next three
#' arguments. Otherwise, will return an object consissting only of these cells
#' @param subset.name Parameter to subset on. Eg, the name of a gene, PC1, a
#' column name in object@@data.info, etc. Any argument that can be retreived
#' using FetchData
#' @param ident.use Create a cell subset based on the provided identity classes
#' @param ident.remove Subtract out cells from these identity classes (used for filtration)
#' @param accept.low Low cutoff for the parameter (default is -Inf)
#' @param accept.high High cutoff for the parameter (default is Inf)
#' @param do.center Recenter the new object@@scale.data
#' @param do.scale Rescale the new object@@scale.data. FALSE by default
#' @param max.cells.per.ident Can be used to downsample the data to a certain max per cell ident. Default is inf.
#' @param random.seed Random seed for downsampling
#' @param \dots Additional arguments to be passed to FetchData (for example,
#' use.imputed=TRUE)
#'
#' @return Returns a Seurat object containing only the relevant subset of cells
#'
#' @export
#'
SubsetData <- function(
  object,
  cells.use = NULL,
  subset.name = NULL,
  ident.use = NULL,
  ident.remove = NULL,
  accept.low = -Inf,
  accept.high = Inf,
  do.center = FALSE,
  do.scale = FALSE,
  max.cells.per.ident = Inf,
  random.seed = 1,
  ...
) {
  data.use <- NULL
  cells.use <- set.ifnull(cells.use,object@cell.names)
  if (!is.null(x = ident.use)) {
    ident.use <- anotinb(ident.use, ident.remove)
    cells.use <- WhichCells(object, ident.use)
  }
  if ((is.null(x = ident.use)) && ! is.null(x = ident.remove)) {
    ident.use <- anotinb(unique(object@ident), ident.remove)
    cells.use <- WhichCells(object, ident.use)
  }
  if (! is.null(x = subset.name)) {
    data.use <- FetchData(object, subset.name, ...)
    if (length(x = data.use) == 0) {
      return(object)
    }
    subset.data <- data.use[, subset.name]
    pass.inds <- which(x = (subset.data > accept.low) & (subset.data < accept.high))
    cells.use <- rownames(data.use)[pass.inds]
  }
  cells.use <- ainb(cells.use, object@cell.names)
  cells.use <-  WhichCells(
    object = object,
    cells.use = cells.use,
    max.cells.per.ident = max.cells.per.ident,
    random.seed = random.seed
  )
  object@data <- object@data[, cells.use]
  if(! is.null(x = object@scale.data)) {
    if (length(x = colnames(x = object@scale.data) > 0)) {
      object@scale.data[, cells.use]
      object@scale.data <- object@scale.data[
        complete.cases(object@scale.data), # Row
        cells.use # Columns
      ]
    }
  }
  if (do.scale) {
    object <- ScaleData(
      object = object,
      do.scale = do.scale,
      do.center = do.center
    )
    object@scale.data <- object@scale.data[
      complete.cases(object@scale.data), # Row
      cells.use # Column
    ]
  }
  object@ident <- drop.levels(x = object@ident[cells.use])
  if (length(x = object@dr) > 0) {
    for (i in 1:length(object@dr)) {
      object@dr[[i]]@rotation <- object@dr[[i]]@rotation[cells.use, ,drop = FALSE]
    }
  }
  # handle multimodal casess
  if (! .hasSlot(object = object, name = "assay")) {
    object@assay <- list()
  }
  if (length(object@assay) > 0) {
    for(i in 1:length(object@assay)) {
      if ((! is.null(x = object@assay[[i]]@raw.data)) && (ncol(x = object@assay[[i]]@raw.data) > 1)) {
        object@assay[[i]]@raw.data <- object@assay[[i]]@raw.data[, cells.use]
      }
      if ((! is.null(x = object@assay[[i]]@data)) && (ncol(x = object@assay[[i]]@data) > 1)) {
        object@assay[[i]]@data <- object@assay[[i]]@data[, cells.use]
      }
      if ((! is.null(x = object@assay[[i]]@scale.data)) && (ncol(x = object@assay[[i]]@scale.data) > 1)) {
        object@assay[[i]]@scale.data <- object@assay[[i]]@scale.data[, cells.use]
      }
    }
  }
  #object@tsne.rot=object@tsne.rot[cells.use, ]
  object@cell.names <- cells.use
  object@gene.scores <- data.frame(object@gene.scores[cells.use,])
  colnames(x = object@gene.scores)[1] <- "nGene"
  rownames(x = object@gene.scores) <- colnames(x = object@data)
  object@data.info <- data.frame(object@data.info[cells.use,])
  #object@mix.probs=data.frame(object@mix.probs[cells.use,]); colnames(object@mix.probs)[1]="nGene"; rownames(object@mix.probs)=colnames(object@data)
  return(object)
}

#' Return a subset of the Seurat object
#'
#' Creates a Seurat object containing only a subset of the cells in the
#' original object. Takes either a list of cells to use as a subset, or a
#' parameter (for example, a gene), to subset on.
#'
#' @param object Seurat object
#' @param cells.use A vector of cell names to use as a subset. If NULL
#' (default), then this list will be computed based on the next three
#' arguments. Otherwise, will return an object consissting only of these cells
#' @param subset.name Parameter to subset on. Eg, the name of a gene, PC1, a
#' column name in object@@data.info, etc. Any argument that can be retreived
#' using FetchData
#' @param accept.low Low cutoff for the parameter (default is -Inf)
#' @param accept.high High cutoff for the parameter (default is Inf)
#' @param do.center Recenter the new object@@scale.data
#' @param do.scale Rescale the new object@@scale.data
#' @param \dots Additional arguments to be passed to FetchData (for example,
#' use.imputed=TRUE)
#'
#' @return Returns a Seurat object containing only the relevant subset of cells
#'
#' @export
#'
SubsetCells <- function(
  object,
  cells.use = NULL,
  subset.name = NULL,
  accept.low = -Inf,
  accept.high = Inf,
  do.center = TRUE,
  do.scale = TRUE,
  ...
) {
  data.use <- FetchData(
    object = object,
    vars.all = subset.name,
    cells.use = cells.use,
    ...
  )
  if (length(x = data.use) == 0) {
    return(object)
  }
  subset.data <- data.use[,subset.name]
  pass.inds <- which(x = (subset.data>accept.low) & (subset.data<accept.high))
  cells.use <- rownames(x = data.use)[pass.inds]
  return(cells.use)
}

#   Documentation
#################
ProjectSamples <- function(object, new.samples) {
  genes.use <- rownames(x = object@data)
  genes.pca <- rownames(x = object@pca.x)
  data.project <- object@scale.data[genes.pca, ]
  data.project[is.na(x = data.project)] <- 0
  new.rot = t(x = data.project) %*% as.matrix(x = object@pca.x)
  scale.vec = apply(
    X = new.rot,
    MARGIN = 2,
    FUN = function(x) {
      return(sqrt(x = sum(x ^ 2)))
    }
  )
  new.rot.scale <- scale(x = new.rot, center = FALSE, scale = scale.vec)
  object@pca.rot <- as.data.frame(x = new.rot.scale)
return(object)
}

#' Run t-distributed Stochastic Neighbor Embedding
#'
#' Run t-SNE dimensionality reduction on selected features. Has the option of running in a reduced
#' dimensional space (i.e. spectral tSNE, recommended), or running based on a set of genes
#'
#' @param object Seurat object
#' @param cells.use Which cells to analyze (default, all cells)
#' @param dims.use Which dimensions to use as input features
#' @param k.seed Random seed for the t-SNE
#' @param do.fast If TRUE, uses the Barnes-hut implementation, which runs
#' faster, but is less flexible
#' @param add.iter If an existing tSNE has already been computed, uses the
#' current tSNE to seed the algorithm and then adds additional iterations on top of this
#' @param genes.use If set, run the tSNE on this subset of genes
#' (instead of running on a set of reduced dimensions). Not set (NULL) by default
#' @param reduction.use Which dimensional reduction (PCA or ICA) to use for the tSNE. Default is PCA
#' @param dim_embed The dimensional space of the resulting tSNE embedding (default is 2).
#' For example, set to 3 for a 3d tSNE
#' @param q.use Quantile to use
#' @param max.dim Max dimension to keep from diffusion calculation
#' @param scale.clip Max/min value for scaled data. Default is 3
#' @param ... Additional arguments to the tSNE call. Most commonly used is
#' perplexity (expected number of neighbors default is 30)
#'
#' @return Returns a Seurat object with a tSNE embedding in object@@tsne_rot
#'
#' @import Rtsne
#' @import tsne
#'
#' @export
#'
RunDiffusion <- function(
  object,
  cells.use = NULL,
  dims.use = 1:5,
  k.seed = 1,
  do.fast = FALSE,
  add.iter = 0,
  genes.use = NULL,
  reduction.use = 'pca',
  dim_embed = 2,
  q.use = 0.01,
  max.dim = 2,
  scale.clip = 10,
  ...
) {
  cells.use <- set.ifnull(cells.use, colnames(x = object@data))
  if (is.null(x = genes.use)) {
    dim.code <- translate.dim.code(object = object, object = reduction.use)
    dim.codes <- paste0(dim.code, dims.use)
    data.use <- FetchData(object = object, vars.all = dim.codes)
  }
  if (! is.null(x = genes.use)) {
    genes.use <- ainb(genes.use, rownames(x = object@scale.data))
    data.use <- minmax(
      data = t(x = object@data[genes.use, cells.use]),
      min = -1 * scale.clip,
      max = scale.clip
    )
  }
  data.dist <- dist(data.use)
  data.diffusion <- data.frame(
    diffuse( # Where is diffuse?
      data.dist,
      neigen = max.dim,
      maxdim = max.dim,
      ...
    )$X
  )
  colnames(x = data.diffusion) <- paste0("DM", 1:ncol(x = data.diffusion))
  rownames(x = data.diffusion) <- cells.use
  for (i in 1:max.dim) {
    x <- data.diffusion[,i]
    x <- minmax(
      data = x,
      min = quantile(x = x, probs = q.use),
      quantile(x = x, probs = 1-q.use)
    )
    data.diffusion[, i] <- x
  }
  object <- SetDimReduction(
    object = object,
    reduction.type = "dm",
    slot = "rotation",
    new.data = as.matrix(x = data.diffusion)
  )
  object <- SetDimReduction(
    object = object,
    reduction.type = "dm",
    slot = "key",
    new.data = "DM"
  )
  return(object)
}

#   Documentation
#################
#Not currently supported
add_tsne <- function(
  object,
  cells.use = NULL,
  pcs.use = 1:10,
  do.plot = TRUE,
  k.seed = 1,
  add.iter = 1000,
  ...
) {
  cells.use <- set.ifnull(cells.use, colnames(x = object@data))
  data.use <- object@pca.rot[cells.use, pcs.use]
  #data.dist=as.dist(mahalanobis.dist(data.use))
  set.seed(seed = k.seed)
  data.tsne <- data.frame(
    tsne(
      X = data.use,
      initial_config = as.matrix(x = object@tsne.rot[cells.use,]),
      max_iter = add.iter,
      ...
    )
  )
  colnames(x = data.tsne) <- paste0("tSNE_", 1:ncol(data.tsne))
  rownames(x = data.tsne) <- cells.use
  object@tsne.rot <- data.tsne
  return(object)
}

#' Probability of detection by identity class
#'
#' For each gene, calculates the probability of detection for each identity
#' class.
#'
#' @param object Seurat object
#' @param thresh.min Minimum threshold to define 'detected' (log-scale)
#'
#' @return Returns a matrix with genes as rows, identity classes as columns.
#'
#' @export
#'
ClusterAlpha <- function(object, thresh.min = 0) {
  ident.use <- object@ident
  data.all <- data.frame(row.names = rownames(x = object@data))
  for (i in sort(x = unique(x = ident.use))) {
    temp.cells <- WhichCells(object = object, ident = i)
    data.temp <- apply(
      X = object@data[, temp.cells],
      MARGIN = 1,
      FUN = function(x) {
        return(sum(x > thresh.min)/length(x = x))
      }
    )
    data.all <- cbind(data.all, data.temp)
    colnames(x = data.all)[ncol(x = data.all)] <- i
  }
  colnames(x = data.all) <- sort(x = unique(x = ident.use))
  return(data.all)
}

#' Reorder identity classes
#'
#' Re-assigns the identity classes according to the average expression of a particular feature (i.e, gene expression, or PC score)
#' Very useful after clustering, to re-order cells, for example, based on PC scores
#'
#' @param object Seurat object
#' @param feature Feature to reorder on. Default is PC1
#' @param rev Reverse ordering (default is FALSE)
#' @param aggregate.fxn Function to evaluate each identity class based on (default is mean)
#' @param reorder.numeric Rename all identity classes to be increasing numbers starting from 1 (default is FALSE)
#' @param \dots additional arguemnts (i.e. use.imputed=TRUE)
#'
#' @return A seurat object where the identity have been re-oredered based on the average.
#'
#' @export
#'
ReorderIdent <- function(
  object,
  feature = "PC1",
  rev = FALSE,
  aggregate.fxn = mean,
  reorder.numeric = FALSE,
  ...
) {
  ident.use <- object@ident
  data.use <- FetchData(object = object, vars.all = feature, ...)[, 1]
  revFxn <- same
  if (rev) {
    revFxn <- function(x) {
      return(max(x) + 1 - x)
    }
  }
  names.sort <- names(
    x = revFxn(
      sort(
        x = tapply(
          X = data.use,
          INDEX = (ident.use),
          FUN = aggregate.fxn
        )
      )
    )
  )
  ident.new <- factor(x = ident.use, levels = names.sort, ordered = TRUE)
  if (reorder.numeric) {
    ident.new <- factor(
      x = revFxn(
        rank(
          tapply(
            X = data.use,
            INDEX = as.numeric(x = ident.new),
            FUN = mean
          )
        )
      )[as.numeric(ident.new)],
      levels = 1:length(x = levels(x = ident.new)),
      ordered = TRUE
    )
  }
  names(x = ident.new) <- names(x = ident.use)
  object@ident <- ident.new
  return(object)
}

#' Average PCA scores by identity class
#'
#' Returns the PCA scores for an 'average' single cell in each identity class
#'
#' @param object Seurat object
#'
#' @return Returns a matrix with genes as rows, identity classes as columns
#'
#' @export
#'
AveragePCA <- function(object) {
  ident.use <- object@ident
  data.all <- data.frame(row.names = colnames(x = object@pca.rot))
  for (i in levels(x = ident.use)) {
    temp.cells <- WhichCells(object = object, ident = i)
    if (length(x = temp.cells) == 1) {
      data.temp <- apply(
        X = data.frame((object@pca.rot[c(temp.cells, temp.cells), ])),
        MARGIN = 2,
        FUN = mean
      )
    }
    if (length(x = temp.cells) > 1) {
      data.temp <- apply(
        X = object@pca.rot[temp.cells, ],
        MARGIN = 2,
        FUN = mean
      )
    }
    data.all <- cbind(data.all, data.temp)
    colnames(x = data.all)[ncol(x = data.all)] <- i
  }
  #colnames(data.all)=levels(ident.use)
  return((data.all))
}

#' Averaged gene expression by identity class
#'
#' Returns gene expression for an 'average' single cell in each identity class
#'
#' Output is in log-space, but averaging is done in non-log space.
#'
#' @param object Seurat object
#' @param genes.use Genes to analyze. Default is all genes.
#' @param return.seurat Whether to return the data as a Seurat object. Default is false.
#' @param add.ident Place an additional label on each cell prior to averaging (very useful if you want to observe cluster averages, separated by replicate, for example)
#' @inheritParams FetchData
#' @param show.progress Show progress bar (default is T)
#' @param ... Arguments to be passed to methods such as \code{\link{Setup}}
#' @return Returns a matrix with genes as rows, identity classes as columns.
#' @export
AverageExpression <- function(
  object,
  genes.use = NULL,
  return.seurat = FALSE,
  add.ident = NULL,
  use.scale = FALSE,
  use.raw = FALSE,
  show.progress = TRUE,
  ...
) {

  ident.orig <- object@ident
  orig.levels <- levels(x = object@ident)
  ident.new <- c()
  if (! is.null(x = add.ident)) {
    new.data <- FetchData(object = object, vars.all = add.ident)
    new.ident <- paste(
      object@ident[rownames(x = new.data)],
      new.data[, 1],
      sep = '_'
    )
    object <- SetIdent(
      object = object,
      cells.use = rownames(x = new.data),
      ident.use = new.ident
    )
  }
  if (return.seurat) {
    assays.use <- c("RNA", names(x = object@assay))
  } else {
    assays.use <- "RNA"
  }
  slot.use <- "data"
  fxn.average <- expMean
  if (show.progress) {
    fxn.loop <- pbsapply
  } else {
    fxn.loop <- sapply
  }
  if (use.scale) {
    slot.use <- "scale.data"
    fxn.average <- mean
  }
  if (use.raw) {
    slot.use <- "raw.data"
    fxn.average <- mean
  }
  data.return <- list()
  for (i in 1:length(x = assays.use)) {
    data.use <- GetAssayData(
      object = object,
      assay.type = assays.use[i],
      slot = slot.use
    )
    genes.assay <- genes.use
    if (length(x = intersect(x = genes.use,y = rownames(x = data.use))) <1 ) {
      genes.assay <- rownames(x = data.use)
    }
    data.all <- data.frame(row.names = genes.assay)
    for (j in levels(x = object@ident)) {
      temp.cells <- WhichCells(object = object, ident = j)
      genes.assay <- unique(x = intersect(x = genes.assay, y = rownames(x = data.use)))
      if (length(x = temp.cells) == 1) {
        data.temp <- (data.use[genes.assay, temp.cells])
      }
      if (length(x = temp.cells) >1 ) {
        data.temp <- apply(
          X = data.use[genes.assay, temp.cells],
          MARGIN = 1,
          FUN = fxn.average
        )
      }
      data.all <- cbind(data.all, data.temp)
      colnames(x = data.all)[ncol(x = data.all)] <- j
      if (show.progress) {
        print(paste0("Finished averaging ", assays.use[i], " for cluster ", j))
      }
      if(i == 1) {
        ident.new <- c(ident.new, as.character(x = ident.orig[temp.cells[1]]))
      }
    }
    names(x = ident.new) <- levels(x = object@ident)
    data.return[[i]] <- data.all
    names(x = data.return)[i] <- assays.use[[i]]
  }
  if (return.seurat) {
    toRet <- new(Class = "seurat", raw.data = data.return[[1]])
    toRet <- Setup(
      object = toRet,
      project = "Average",
      min.cells = 0,
      min.genes = 0,
      is.expr = 0,
      ...
    )
    #for multimodal data
    if (length(x = data.return) > 1) {
      for (i in 2:length(x = data.return)) {
        toRet <- SetAssayData(
          object = toRet,
          assay.type = names(x = data.return)[i],
          slot = "raw.data",
          new.data = data.return[[i]]
        )
      }
    }
    toRet <- SetIdent(
      object = toRet,
      cells.use = toRet@cell.names,
      ident.use = ident.new[toRet@cell.names]
    )
    toRet@ident <- factor(
      x = toRet@ident,
      levels = as.character(x = orig.levels),
      ordered = TRUE
    )
    return(toRet)
  } else {
    return(data.return[[1]])
  }
}

#' Access cellular data
#'
#' Retreives data (gene expression, PCA scores, etc, metrics, etc.) for a set
#' of cells in a Seurat object
#'
#' @param object Seurat object
#' @param vars.all List of all variables to fetch
#' @param cells.use Cells to collect data for (default is all cells)
#' @param use.imputed For gene expression, use imputed values. Default is FALSE
#' @param use.scaled For gene expression, use scaled values. Default is FALSE
#' @param use.raw For gene expression, use raw values. Default is FALSE
#'
#' @return A data frame with cells as rows and cellular data as columns
#'
#' @export
#'
FetchData <- function(
  object,
  vars.all = NULL,
  cells.use = NULL,
  use.imputed = FALSE,
  use.scaled = FALSE,
  use.raw = FALSE
) {
  cells.use <- set.ifnull(x = cells.use, y = object@cell.names)
  data.return <- data.frame(row.names = cells.use)
  data.expression <- as.matrix(x = data.frame(row.names = cells.use))
  # if any vars passed are genes, subset expression data
  gene.check <- vars.all %in% rownames(object@data)
  #data.expression <- matrix()
  if (all(gene.check)){
    if (use.imputed) {
      data.expression <- object@imputed[vars.all, cells.use,drop = FALSE]
    }
    if (use.scaled) {
      data.expression <- object@scale.data[vars.all, cells.use, drop = FALSE]
    }
    if (use.raw) {
      data.expression <-  object@raw.data[vars.all, cells.use, drop = FALSE]
    } else {
      data.expression <- object@data[vars.all, cells.use, drop = FALSE ]
    }
    return(t(x = as.matrix(x = data.expression)))
  } else if (any(gene.check)) {
    if (use.imputed) {
      data.expression <- object@imputed[vars.all[gene.check], cells.use, drop = FALSE]
    }
    if(use.scaled) {
      data.expression <-  object@scale.data[vars.all[gene.check], cells.use, drop = FALSE]
    }
    if (use.raw) {
      data.expression <- object@raw.data[vars.all[gene.check], cells.use, drop = FALSE]
    } else {
      data.expression <- object@data[vars.all[gene.check], cells.use, drop = FALSE]
    }
    data.expression <- t(x = data.expression)
  }
  #now check for multimodal data
  if (length(x = object@assay) > 0) {
    data.types <- names(x = object@assay)
    slot.use <- "data"
    if (use.scaled) {
      slot.use <- "scale.data"
    }
    if (use.raw) {
      slot.use <- "raw.data"
    }
    for (data.type in data.types) {
      all_data <- (GetAssayData(
        object = object,
        assay.type = data.type,
        slot = slot.use
      ))
      genes.include <- intersect(x = vars.all, y = rownames(x = all_data))
      data.expression <- cbind(
        data.expression,
        t(x = all_data[genes.include, , drop = FALSE])
      )
    }
  }
  var.options <- c("data.info", "mix.probs", "gene.scores")
  if (length(x = names(x = object@dr)) > 0) {
    dr.options <- names(x = object@dr)
    dr.names <- paste0("dr$", names(x = object@dr), "@key")
    dr.names <- sapply(
      X = dr.names,
      FUN = function(x) {
        return(eval(expr = parse(text = paste0("object@", x))))
      }
    )
    names(x = dr.names) <- dr.options
    var.options <- c(var.options, dr.names)
  }
  object@data.info[,"ident"] <- object@ident[rownames(x = object@data.info)]
  for (my.var in vars.all) {
    data.use=data.frame()
    if (my.var %in% colnames(data.expression)) {
      data.use <- data.expression
    } else {
      for(i in var.options) {
        if (all(unlist(x = strsplit(x = my.var, split = "[0-9]+")) == i)) {
          eval(
            expr = parse(
              text = paste0(
                "data.use <- object@dr$",
                names(x = var.options[which(i == var.options)]),
                "@rotation"
              )
            )
          )
          colnames(x = data.use) <- paste0(i, 1:ncol(x = data.use))
          break
        }
      }
    }
    if (my.var %in% colnames(object@data.info)) {
      data.use <- object@data.info[, my.var, drop = FALSE]
    }
    if (ncol(x = data.use) == 0) {
      stop(paste("Error:", my.var, "not found"))
    }
    cells.use <- ainb(a = cells.use, b = rownames(x = data.use))
    if (! my.var %in% colnames(x = data.use)) {
      stop(paste("Error:", my.var, "not found"))
    }
    data.add <- data.use[cells.use, my.var]
    if (is.null(x = data.add)) {
      stop(paste("Error:", my.var, "not found"))
    }
    data.return <- cbind(data.return, data.add)
  }
  colnames(x = data.return) <- vars.all
  rownames(x = data.return) <- cells.use
  return(data.return)
}

#   Documentation...
####################
GetWeightMatrix <- function(object) {
  data <- object@data
  data.humpAvg <- apply(
    X = data,
    MARGIN = 1,
    FUN = humpMean,
    min = object@drop.expr
  )
  wt.matrix <- data.frame(
    t(
      x = pbsapply(
        X = data.humpAvg,
        FUN = expAlpha,
        coefs = object@drop.coefs
      )
    )
  )
  colnames(x = wt.matrix) <- colnames(x = data)
  rownames(x = wt.matrix) <- rownames(x = data)
  wt.matrix[is.na(x = wt.matrix)] <- 0
  object@wt.matrix <- wt.matrix
  wt1.matrix <- data.frame(
    pbsapply(
      X = 1:ncol(x = data),
      FUN = function(x) {
        setWt1(
          x = data[, x],
          wts = wt.matrix[, x],
          min = object@drop.expr
        )
      }
    )
  )
  colnames(x = wt1.matrix) <- colnames(x = data)
  rownames(x = wt1.matrix) <- rownames(x = data)
  wt1.matrix[is.na(x = wt1.matrix)] <- 0
  object@drop.wt.matrix <- wt1.matrix
  return(object)
}

#   Documentation...
####################
RegulatorScore <- function(
  object,
  candidate.reg,
  score.name,
  cells.use = NULL
) {
  cells.use <- set.ifnull(x = cells.use, y = colnames(x = object@data))
  candidate.reg <- candidate.reg[candidate.reg %in% rownames(x = object@data)]
  my.score <- retreiveScore(object, score.name)[cells.use]
  my.data <- object@data[, cells.use]
  my.ident <- object@ident[cells.use]
  reg.score <- unlist(
    lapply(
      X = candidate.reg,
      FUN = regressionSig,
      score = my.score,
      data = my.data,
      latent = my.ident,
      code = "rsem"
    )
  )
  names(x = reg.score) <- candidate.reg
  return(reg.score)
}

#' Gene expression markers of identity classes defined by a phylogenetic clade
#'
#' Finds markers (differentially expressed genes) based on a branching point (node) in
#' the phylogenetic tree. Markers that define clusters in the left branch are positive markers.
#' Markers that define the right branch are negative markers.
#'
#' @inheritParams FindMarkers
#' @param node The node in the phylogenetic tree to use as a branch point
#' @param tree.use Can optionally pass the tree to be used. Default uses the tree in object@@cluster.tree
#' @param ... Additional arguments passed to FindMarkers
#'
#' @return Matrix containing a ranked list of putative markers, and associated
#' statistics (p-values, ROC score, etc.)
#'
#' @export
#'
FindMarkersNode <- function(
  object,
  node,
  tree.use = NULL,
  genes.use = NULL,
  thresh.use = 0.25,
  test.use = "bimod",
  ...
) {
  genes.use <- set.ifnull(x = genes.use, y = rownames(x = object@data))
  tree <- set.ifnull(x = tree.use, y = object@cluster.tree[[1]])
  ident.order <- tree$tip.label
  nodes.1 <- ident.order[getLeftDecendants(tree = tree, node = node)]
  nodes.2 <- ident.order[getRightDecendants(tree = tree, node = node)]
  #print(nodes.1)
  #print(nodes.2)
  to.return <- FindMarkers(
    object = object,
    ident.1 = nodes.1,
    ident.2 = nodes.2,
    genes.use = genes.use,
    thresh.use = thresh.use,
    test.use = test.use,
    ...
  )
  return(to.return)
}

#' Gene expression markers of identity classes
#'
#' Finds markers (differentially expressed genes) for identity classes
#'
#' @param object Seurat object
#' @param ident.1 Identity class to define markers for
#' @param ident.2 A second identity class for comparison. If NULL (default) -
#' use all other cells for comparison.
#' @param genes.use Genes to test. Default is to use all genes
#' @param thresh.use Limit testing to genes which show, on average, at least
#' X-fold difference (log-scale) between the two groups of cells. Default is 0.25
#' Increasing thresh.use speeds up the function, but can miss weaker signals.
#' @param test.use Denotes which test to use. Seurat currently implements
#' "bimod" (likelihood-ratio test for single cell gene expression, McDavid et
#' al., Bioinformatics, 2011, default), "roc" (standard AUC classifier), "t"
#' (Students t-test), and "tobit" (Tobit-test for differential gene expression,
#' as in Trapnell et al., Nature Biotech, 2014), 'poisson', and 'negbinom'.
#' The latter two options should only be used on UMI datasets, and assume an underlying
#' poisson or negative-binomial distribution
#' @param min.pct - only test genes that are detected in a minimum fraction of min.pct cells
#' in either of the two populations. Meant to speed up the function by not testing genes that are very infrequently expressed. Default is 0.1
#' @param min.diff.pct - only test genes that show a minimum difference in the fraction of detection between the two groups. Set to -Inf by default
#' @param only.pos Only return positive markers (FALSE by default)
#' @param print.bar Print a progress bar once expression testing begins (uses pbapply to do this)
#' @param max.cells.per.ident Down sample each identity class to a max number. Default is no downsampling. Not activated by default (set to Inf)
#' @param random.seed Random seed for downsampling
#' @param min.cells Minimum number of cells expressing the gene in at least one of the two groups
#'
#' @return Matrix containing a ranked list of putative markers, and associated statistics (p-values, ROC score, etc.)
#'
#' @import VGAM
#' @import pbapply
#'
#' @export
#'
FindMarkers <- function(
  object,
  ident.1,
  ident.2 = NULL,
  genes.use = NULL,
  thresh.use = 0.25,
  test.use = "bimod",
  min.pct = 0.1,
  min.diff.pct = -Inf,
  print.bar = TRUE,
  only.pos = FALSE,
  max.cells.per.ident = Inf,
  random.seed = 1,
  latent.vars = "nUMI",
  min.cells = 3
) {
  genes.use <- set.ifnull(x = genes.use, y = rownames(x = object@data))
  if (max.cells.per.ident < Inf) {
    object <- SubsetData(
      object = object,
      max.cells.per.ident = max.cells.per.ident,
      random.seed = random.seed
    )
  }
  # in case the user passed in cells instead of identity classes
  # Umm... confused about the first part of the conditional...
  if (length(x = as.vector(x = ident.1) > 1) && any(as.character(x = ident.1) %in% object@cell.names)) {
    cells.1 <- ainb(a = ident.1, b = object@cell.names)
  } else {
    cells.1 <- WhichCells(object = object, ident = ident.1)
  }
  # if NULL for ident.2, use all other cells
  if (length(x = as.vector(x = ident.2) > 1) && any(as.character(x = ident.2) %in% object@cell.names)) {
    cells.2 <- ainb(a = ident.2, b = object@cell.names)
  } else {
    if (is.null(x = ident.2)) {
      cells.2 <- object@cell.names
    } else {
      cells.2 <- WhichCells(object = object, ident = ident.2)
    }
  }
  cells.2 <- anotinb(x = cells.2, y = cells.1)
  #error checking
  if (length(x = cells.1) == 0) {
    print(paste("Cell group 1 is empty - no cells with identity class", ident.1))
    return(NULL)
  }
  if (length(x = cells.2) == 0) {
    print(paste("Cell group 2 is empty - no cells with identity class", ident.2))
    return(NULL)
  }
  #gene selection (based on percent expressed)
  thresh.min <- object@is.expr
  data.temp1 <- round(
    x = apply(
      X = object@data[genes.use, cells.1],
      MARGIN = 1,
      FUN = function(x) {
        return(sum(x > thresh.min) / length(x = x))
        # return(length(x = x[x>thresh.min]) / length(x = x))
      }
    ),
    digits = 3
  )
  data.temp2 <- round(
    x = apply(
      X = object@data[genes.use, cells.2],
      MARGIN = 1,
      FUN = function(x) {
        return(sum(x > thresh.min) / length(x = x))
        # return(length(x = x[x > thresh.min]) / length(x = x))
      }
    ),
    digits = 3
  )
  data.alpha <- cbind(data.temp1, data.temp2)
  colnames(x = data.alpha) <- c("pct.1","pct.2")
  alpha.min <- apply(X = data.alpha, MARGIN = 1, FUN = max)
  names(x = alpha.min) <- rownames(x = data.alpha)
  genes.use <- names(x = which(x = alpha.min > min.pct))
  alpha.diff <- alpha.min - apply(X = data.alpha, MARGIN = 1, FUN = min)
  genes.use <- names(
    x = which(x = alpha.min > min.pct & alpha.diff > min.diff.pct)
  )
  #gene selection (based on average difference)
  data.1 <- apply(X = object@data[genes.use, cells.1], MARGIN = 1, FUN = expMean)
  data.2 <- apply(X = object@data[genes.use, cells.2], MARGIN = 1, FUN = expMean)
  total.diff <- (data.1 - data.2)
  genes.diff <- names(x = which(x = abs(x = total.diff) > thresh.use))
  genes.use <- ainb(a = genes.use, b = genes.diff)
  #perform DR
  if (test.use == "bimod") {
    to.return <- DiffExpTest(
      object = object,
      cells.1 = cells.1,
      cells.2 = cells.2,
      genes.use = genes.use,
      print.bar = print.bar
    )
  }
  if (test.use == "roc") {
    to.return <- MarkerTest(
      object = object,
      cells.1 = cells.1,
      cells.2 = cells.2,
      genes.use = genes.use,
      print.bar = print.bar
    )
  }
  if (test.use == "t") {
    to.return <- DiffTTest(
      object = object,
      cells.1 = cells.1,
      cells.2 = cells.2,
      genes.use = genes.use,
      print.bar = print.bar
    )
  }
  if (test.use == "tobit") {
    to.return <- TobitTest(
      object = object,
      cells.1 = cells.1,
      cells.2 = cells.2,
      genes.use = genes.use,
      print.bar = print.bar
    )
  }
  if (test.use == "negbinom") {
    to.return <- NegBinomDETest(
      object = object,
      cells.1 = cells.1,
      cells.2 = cells.2,
      genes.use = genes.use,
      latent.vars = latent.vars,
      print.bar = print.bar,
      min.cells = min.cells
    )
  }
  if (test.use == "poisson") {
    to.return <- PoissonDETest(
      object = object,
      cells.1 = cells.1,
      cells.2 = cells.2,
      genes.use = genes.use,
      latent.vars = latent.vars,
      print.bar = print.bar,
      min.cells # PoissonDETest doesn't have something for min.cells
    )
  }
  #return results
  to.return[, "avg_diff"] <- total.diff[rownames(x = to.return)]
  to.return <- cbind(to.return, data.alpha[rownames(x = to.return), ])
  if (test.use == "roc") {
    to.return <- to.return[order(-to.return$power, -to.return$avg_diff), ]
  } else {
    to.return <- to.return[order(to.return$p_val, -to.return$avg_diff), ]
  }
  if (only.pos) {
    to.return <- subset(x = to.return, subset = avg_diff > 0)
  }
  return(to.return)
}

#' Finds markers that are conserved between the two groups
#'
#' @param object Seurat object
#' @param ident.1 Identity class to define markers for
#' @param ident.2 A second identity class for comparison. If NULL (default) - use all other cells
#' for comparison.
#' @param grouping.var grouping variable
#' @param \dots parameters to pass to FindMarkers
#'
#' @return Matrix containing a ranked list of putative conserved markers, and associated statistics
#' (p-values within each group and a combined p-value (fisher_pval), percentage of cells expressing
#' the marker, average differences)
#'
#' @export
#'
FindConservedMarkers <- function(
  object,
  ident.1,
  ident.2 = NULL,
  grouping.var,
  ...
) {
  object.var <- FetchData(object = object, vars.all = grouping.var)
  object <- SetIdent(
    object = object,
    cells.use = object@cell.names,
    ident.use = paste(object@ident, object.var[, 1], sep = "_")
  )
  levels.split <- names(x = sort(x = table(object.var[, 1])))
  if (length(x = levels.split) != 2) {
    stop(
      paste0(
        "There are not two options for ",
        grouping.var,
        ". \n Current groups include: ",
        paste(levels.split, collapse = ", ")
      )
    )
  }
  cells <- list()
  for (i in 1:2) {
    cells[[i]] <- rownames(
      x = object.var[object.var[, 1] == levels.split[i], , drop = FALSE]
    )
  }
  marker.test <- list()
  # do marker tests
  for (i in 1:2) {
    level.use <- levels.split[i]
    ident.use.1 <- paste(ident.1, level.use, sep = "_")
    cells.1 <- WhichCells(object = object, ident = ident.use.1)
    if (is.null(x = ident.2)) {
      cells.2 <- anotinb(x = cells[[i]], y = cells.1)
      ident.use.2 <- names(x = which(x = table(object@ident[cells.2]) > 0))
      if (length(x = ident.use.2) == 0) {
        stop(paste("Only one identity class present:", ident.1))
      }
    }
    if (! is.null(x = ident.2)) {
      ident.use.2 <- paste(ident.2, level.use, sep = "_")
    }
    cat(
      paste0(
        "Testing ",
        ident.use.1,
        " vs ",
        paste(ident.use.2, collapse = ", "), "\n"
      ),
      file = stderr()
    )
    marker.test[[i]] <- FindMarkers(
      object = object,
      ident.1 = ident.use.1,
      ident.2 = ident.use.2,
      ...
    )
  }
  genes.conserved <- intersect(
    x = rownames(x = marker.test[[1]]),
    y = rownames(x = marker.test[[2]])
  )
  markers.conserved <- list()
  for (i in 1:2) {
    markers.conserved[[i]] <- marker.test[[i]][genes.conserved, ]
    colnames(x = markers.conserved[[i]]) <- paste(
      levels.split[i],
      colnames(x = markers.conserved[[i]]),
      sep="_"
    )
  }
  markers.combined <- cbind(markers.conserved[[1]], markers.conserved[[2]])
  pval.codes <- paste(levels.split, "p_val", sep = "_")
  markers.combined$max_pval <- apply(
    X = markers.combined[, pval.codes],
    MARGIN = 1,
    FUN = max
  )
  markers.combined$fisher_pval <- apply(
    X = markers.combined[, pval.codes],
    MARGIN = 1,
    FUN = fisher.integrate
  )
  markers.combined <- markers.combined[order(markers.combined$fisher_pval), ]
  return(markers.combined)
}

#   Documentation...
####################
fisher.integrate <- function(pvals) {
  return(1 - pchisq(q = -2 * sum(log(x = pvals)), df = 2 * length(x = pvals)))
}

#' Gene expression markers for all identity classes
#'
#' Finds markers (differentially expressed genes) for each of the identity classes in a dataset
#'
#' @param object Seurat object
#' @param ident.1 Identity class to define markers for
#' @param ident.2 A second identity class for comparison. If NULL (default) -
#' use all other cells for comparison.
#' @param genes.use Genes to test. Default is to all genes
#' @param thresh.use Limit testing to genes which show, on average, at least
#' X-fold difference (log-scale) between the two groups of cells.
#' Increasing thresh.use speeds up the function, but can miss weaker signals.
#' @param test.use Denotes which test to use. Seurat currently implements
#' "bimod" (likelihood-ratio test for single cell gene expression, McDavid et
#' al., Bioinformatics, 2011, default), "roc" (standard AUC classifier), "t"
#' (Students t-test), and "tobit" (Tobit-test for differential gene expression,
#' as in Trapnell et al., Nature Biotech, 2014), 'poisson', and 'negbinom'.
#' The latter two options should only be used on UMI datasets, and assume an underlying
#' poisson or negative-binomial distribution
#' @param min.pct - only test genes that are detected in a minimum fraction of min.pct cells
#' in either of the two populations. Meant to speed up the function by not testing genes that are very infrequently expression
#' @param min.diff.pct - only test genes that show a minimum difference in the fraction of detection between the two groups. Set to -Inf by default
#' @param only.pos Only return positive markers (FALSE by default)
#' @param print.bar Print a progress bar once expression testing begins (uses pbapply to do this)
#' @param max.cells.per.ident Down sample each identity class to a max number. Default is no downsampling.
#' @param random.seed Random seed for downsampling
#' @param return.thresh Only return markers that have a p-value < return.thresh, or a power > return.thresh (if the test is ROC)
#' @param do.print FALSE by default. If TRUE, outputs updates on progress.
#' @param min.cells Minimum number of cells expressing the gene in at least one of the two groups
#' @param latent.vars remove the effects of these variables
#'
#' @return Matrix containing a ranked list of putative markers, and associated
#' statistics (p-values, ROC score, etc.)
#'
#' @export
#'
FindAllMarkers <- function(
  object,
  ident.1,
  ident.2 = NULL,
  genes.use = NULL,
  thresh.use = 0.25,
  test.use = "bimod",
  min.pct = 0.1,
  min.diff.pct = 0.05,
  print.bar = TRUE,
  only.pos = FALSE,
  max.cells.per.ident = Inf,
  return.thresh = 1e-2,
  do.print = FALSE,
  random.seed = 1,
  min.cells = 3,
  latent.vars = "nUMI"
) {
  genes.use <- set.ifnull(x = genes.use, y = rownames(x = object@data))
  ident.use <- object@ident
  if ((test.use == "roc") && (return.thresh == 1e-2)) {
    return.thresh = 0.7
  }
  idents.all <- sort(x = unique(x = object@ident))
  genes.de <- list()
  if (max.cells.per.ident < Inf) {
    object <- SubsetData(
      object = object,
      max.cells.per.ident = max.cells.per.ident,
      random.seed = random.seed
    )
  }
  for (i in 1:length(x = idents.all)) {
    genes.de[[i]] <- FindMarkers(
      object = object,
      ident.1 = idents.all[i],
      ident.2 = NULL,
      genes.use = genes.use,
      thresh.use = thresh.use,
      test.use = test.use,
      min.pct = min.pct,
      min.diff.pct = min.diff.pct,
      print.bar = print.bar,
      min.cells = min.cells,
      latent.vars = latent.vars
    )
    if (do.print) {
      print(paste("Calculating cluster", idents.all[i]))
    }
  }
  gde.all <- data.frame()
  for (i in 1:length(x = idents.all)) {
    gde <- genes.de[[i]]
    if (nrow(x = gde) > 0) {
      if (test.use == "roc") {
        gde <- subset(
          x = gde,
          subset = (myAUC > return.thresh | myAUC < (1 - return.thresh))
        )
      } else {
        gde <- gde[order(gde$p_val, -gde$avg_diff), ]
        gde <- subset(x = gde, subset = p_val < return.thresh)
      }
      if (nrow(x = gde) > 0) {
        gde$cluster <- idents.all[i]
        gde$gene <- rownames(x = gde)
      }
      if (nrow(x = gde) > 0) {
        gde.all <- rbind(gde.all, gde)
      }
    }
  }
  if (only.pos) {
    return(subset(x = gde.all, subset = avg_diff > 0))
  }
  return(gde.all)
}

#' Likelihood ratio test for zero-inflated data
#'
#' Identifies differentially expressed genes between two groups of cells using
#' the LRT model proposed in Mcdavid et al, Bioinformatics, 2011
#'
#' @inheritParams FindMarkers
#' @param object Seurat object
#' @param cells.1 Group 1 cells
#' @param cells.2 Group 2 cells
#' @return Returns a p-value ranked matrix of putative differentially expressed
#' genes.
#'
#' @export
#'
DiffExpTest <- function(
  object,
  cells.1,
  cells.2,
  genes.use = NULL,
  print.bar = TRUE
) {
  genes.use <- set.ifnull(x = genes.use, y = rownames(x = object@data))
  if (print.bar) {
    iterate.fxn <- pblapply
  } else {
    iterate.fxn <- lapply
  }
  p_val <- unlist(
    x = iterate.fxn(
      X = genes.use,
      FUN = function(x) {
        return(
          diffLRT(
            x = as.numeric(x = object@data[x, cells.1]),
            y = as.numeric(x = object@data[x, cells.2])
          )
        )
      }
    )
  )
  to.return <- data.frame(p_val, row.names = genes.use)
  return(to.return)
}

#' Negative binomial test for UMI-count based data
#'
#' Identifies differentially expressed genes between two groups of cells using
#' a negative binomial generalized linear model
#
#'
#' @inheritParams FindMarkers
#' @param object Seurat object
#' @param cells.1 Group 1 cells
#' @param cells.2 Group 2 cells
#'
#' @return Returns a p-value ranked matrix of putative differentially expressed
#' genes.
#'
#' @importFrom MASS glm.nb
#' @importFrom pbapply pbapply
#'
#' @export
#'
NegBinomDETest <- function(
  object,
  cells.1,
  cells.2,
  genes.use = NULL,
  latent.vars = NULL,
  print.bar = TRUE,
  min.cells = 3
) {
  genes.use <- set.ifnull(x = genes.use, y = rownames(x = object@data))
  # check that the gene made it through the any filtering that was done
  genes.use <- genes.use[genes.use %in% rownames(x = object@data)]
  my.latent <- FetchData(
    object = object,
    vars.all = latent.vars,
    cells.use = c(cells.1, cells.2),
    use.raw = TRUE
  )
  to.test.data <- object@raw.data[genes.use, c(cells.1, cells.2)]
  to.test <- data.frame(my.latent, row.names = c(cells.1, cells.2))
  to.test[cells.1, "group"] <- "A"
  to.test[cells.2, "group"] <- "B"
  to.test$group <- factor(x = to.test$group)
  latent.vars <- c("group", latent.vars)
  if (print.bar) {
    iterate.fxn <- pblapply
  } else {
    iterate.fxn <- lapply
  }
  p_val <- unlist(
    x = iterate.fxn(
      X = genes.use,
      FUN = function(x) {
        to.test[, "GENE"] <- as.numeric(x = to.test.data[x, ])
        # check that gene is expressed in specified number of cells in one group
        if (sum(to.test$GENE[to.test$group == "A"]) < min.cells ||
            sum(to.test$GENE[to.test$group == "B"]) < min.cells) {
          warning(paste0(
            "Skipping gene --- ",
            x,
            ". Fewer than ",
            min.cells,
            " in at least one of the two clusters."
          ))
          return(2)
        }
        # check that variance between groups is not 0
        if (var(x = to.test$GENE) == 0) {
          warning(paste0(
            "Skipping gene -- ",
            x,
            ". No variance in expression between the two clusters."
          ))
          return(2)
        }
        fmla <- as.formula(paste0("GENE ", " ~ ", paste(latent.vars, collapse = "+")))
        p.estimate <- 2
        try(
          expr = p.estimate <- summary(
            object = glm.nb(formula = fmla, data = to.test)
          )$coef[2, 4],
          silent = TRUE
        )
        return(p.estimate)
      }
    )
  )
  if (length(x = which(x = p_val == 2)) > 0){
    genes.use <- genes.use[-which(x = p_val == 2)]
    p_val <- p_val[! p_val == 2]
  }
  to.return <- data.frame(p_val, row.names = genes.use)
  return(to.return)
}

#' Negative binomial test for UMI-count based data (regularized version)
#'
#' Identifies differentially expressed genes between two groups of cells using
#' a likelihood ratio test of negative binomial generalized linear models where
#' the overdispersion parameter theta is determined by pooling information
#' across genes.
#'
#' @inheritParams FindMarkers
#' @param object Seurat object
#' @param cells.1 Group 1 cells
#' @param cells.2 Group 2 cells
#'
#' @return Returns a p-value ranked data frame of test results.
#'
#' @export
#'
NegBinomRegDETest <- function(
  object,
  cells.1,
  cells.2,
  genes.use = NULL,
  latent.vars = NULL,
  print.bar = TRUE,
  min.cells = 3
) {
  genes.use <- set.ifnull(x = genes.use, y = rownames(x = object@data))
  # check that the gene made it through the any filtering that was done
  genes.use <- genes.use[genes.use %in% rownames(x = object@data)]
  print(
    sprintf(
      'NegBinomRegDETest for %d genes and %d and %d cells',
      length(x = genes.use),
      length(x = cells.1),
      length(x = cells.2)
    )
  )
  grp.fac <- factor(
    x = c(
      rep.int(x = 'A', times = length(x = cells.1)),
      rep.int(x = 'B', times = length(x = cells.2))
    )
  )
  to.test.data <- object@raw.data[genes.use, c(cells.1, cells.2), drop = FALSE]
  print('Calculating mean per gene per group')
  above.threshold <- pmax(
    apply(X = to.test.data[, cells.1] > 0, MARGIN = 1, FUN = mean),
    apply(X = to.test.data[, cells.2] > 0, MARGIN = 1, FUN = mean)
  ) >= 0.02
  print(
    sprintf(
      '%d genes are detected in at least 2%% of the cells in at least one of the groups and will be tested',
      sum(above.threshold)
    )
  )
  genes.use <- genes.use[above.threshold]
  to.test.data <- to.test.data[genes.use, , drop = FALSE]
  my.latent <- FetchData(
    object = object,
    vars.all = latent.vars,
    cells.use = c(cells.1, cells.2),
    use.raw = TRUE
  )
  to.test <- data.frame(my.latent, row.names = c(cells.1, cells.2))
  print(paste('Latent variables are', latent.vars))
  # get regularized theta (ignoring group factor)
  theta.fit <- theta.reg(
    cm = to.test.data,
    latent.data = to.test,
    min.theta = 0.01,
    bin.size = 128
  )
  print('Running NB regression model comparison')
  to.test$NegBinomRegDETest.group <- grp.fac
  bin.size <- 128
  bin.ind <- ceiling(1:length(x = genes.use) / bin.size)
  max.bin <- max(bin.ind)
  pb <- txtProgressBar(min = 0, max = max.bin, style = 3)
  res <- c()
  for (i in 1:max.bin) {
    genes.bin.use <- genes.use[bin.ind == i]
    bin.out.lst <- parallel::mclapply(
      X = genes.bin.use,
      FUN = function(j) {
        return(de.nb.reg(
          y = to.test.data[j, ],
          theta = theta.fit[j],
          latent.data = to.test,
          com.frac = latent.vars,
          grp.frac = 'NegBinomRegDETest.group'
        ))
      }
    )
    res <- rbind(res, do.call(rbind, bin.out.lst))
    setTxtProgressBar(pb = pb, value = i)
  }
  close(pb)
  rownames(res) <- genes.use
  res <- as.data.frame(x = res)
  res$adj.pval <- p.adjust(p = res$pval, method='fdr')
  res <- res[order(res$pval, -abs(x = res$log.fc)), ]
  return(res)
}

#' Poisson test for UMI-count based data
#'
#' Identifies differentially expressed genes between two groups of cells using
#' a poisson generalized linear model
#
#' @inheritParams FindMarkers
#' @param object Seurat object
#' @param cells.1 Group 1 cells
#' @param cells.2 Group 2 cells
#'
#' @return Returns a p-value ranked matrix of putative differentially expressed
#' genes.
#'
#' @importFrom pbapply pbapply
#'
#' @export
#'
PoissonDETest <- function(
  object,
  cells.1,
  cells.2,
  genes.use = NULL,
  latent.vars = NULL,
  print.bar = TRUE
) {
  genes.use <- set.ifnull(x = genes.use, y = rownames(x = object@data))
  # check that the gene made it through the any filtering that was done
  genes.use <- genes.use[genes.use %in% rownames(x = object@data)]
  my.latent <- FetchData(
    object = object,
    vars.all = latent.vars,
    cells.use = c(cells.1, cells.2),
    use.raw = TRUE
  )
  to.test.data <- object@raw.data[genes.use, c(cells.1, cells.2)]
  to.test <- data.frame(my.latent, row.names = c(cells.1, cells.2))
  to.test[cells.1,"group"] <- "A"
  to.test[cells.2,"group"] <- "B"
  to.test$group <- factor(x = to.test$group)
  latent.vars <- c("group", latent.vars)
  if (print.bar) {
    iterate.fxn <- pblapply
  } else {
    iterate.fxn <- lapply
  }
  p_val <- unlist(
    x = iterate.fxn(
      X = genes.use,
      FUN = function(x) {
        to.test[,"GENE"] <- as.numeric(x = to.test.data[x, ])
        # check that gene is expressed in specified number of cells in one group
        if (sum(to.test$GENE[to.test$group == "A"]) < min.cells ||
            sum(to.test$GENE[to.test$group == "B"]) < min.cells) {
          warning(paste0(
            "Skipping gene --- ",
            x,
            ". Fewer than",
            min.cells,
            " in at least one of the two clusters."
          ))
          return(2)
        }
        # check that variance between groups is not 0
        if (var(to.test$GENE) == 0) {
          print("what") # what?
          warning(paste0(
            "Skipping gene -- ",
            x,
            ". No variance in expression between the two clusters."
          ))
          return(2)
        }
        fmla <- as.formula(
          object = paste0("GENE ", " ~ ", paste(latent.vars, collapse="+"))
        )
        return(
          summary(
            object = glm(
              formula = fmla,
              data = to.test,
              family = "poisson"
            )
          )$coef[2,4]
        )
      }
    )
  )
  if (length(x = which(x = p_val == 2)) > 0) {
    genes.use <- genes.use[-which(x = p_val == 2)]
    p_val <- p_val[! p_val == 2]
  }
  to.return <- data.frame(p_val, row.names = genes.use)
  return(to.return)
}

#' Differential expression testing using Tobit models
#'
#' Identifies differentially expressed genes between two groups of cells using
#' Tobit models, as proposed in Trapnell et al., Nature Biotechnology, 2014
#'
#' @inheritParams FindMarkers
#' @inheritParams DiffExpTest
#'
#' @return Returns a p-value ranked matrix of putative differentially expressed
#' genes.
#'
#' @export
#'
TobitTest <- function(
  object,
  cells.1,
  cells.2,
  genes.use = NULL,
  print.bar = TRUE
) {
  genes.use <- set.ifnull(x = genes.use, y = rownames(x = object@data))
  #print(genes.diff)
  to.return <- TobitDiffExpTest(
    data1 = object@data[, cells.1],
    data2 = object@data[, cells.2],
    mygenes = genes.use,
    print.bar = print.bar
  )
  return(to.return)
}

#' Identify potential genes associated with batch effects
#'
#' Test for genes whose expression value is strongly predictive of batch (based
#' on ROC classification). Important note: Assumes that the 'batch' of each
#' cell is assigned to the cell's identity class (will be improved in a future
#' release)
#'
#' @param object Seurat object
#' @param idents.use Batch names to test
#' @param genes.use Gene list to test
#' @param auc.cutoff Minimum AUC needed to qualify as a 'batch gene'
#' @param thresh.use Limit testing to genes which show, on average, at least
#' X-fold difference (log-scale) in any one batch
#'
#' @return Returns a list of genes that are strongly correlated with batch.
#'
#' @export
#'
BatchGene <- function(
  object,
  idents.use,
  genes.use = NULL,
  auc.cutoff = 0.6,
  thresh.use = 0
) {
  batch.genes <- c()
  genes.use <- set.ifnull(x = genes.use, y = rownames(x = object@data))
  for (ident in idents.use) {
    cells.1 <- names(x = object@ident)[object@ident == ident]
    cells.2 <- names(x = object@ident)[object@ident != ident]
    if ((length(x = cells.1) < 5) | (length(x = cells.2) < 5)) {
      break
    }
    markers.ident <- MarkerTest(
      object = object,
      cells.1 = cells.1,
      cells.2 = cells.2,
      genes.use = genes.use,
      print.bar = thresh.use
    )
    batch.genes <- unique(
      x = c(
        batch.genes,
        rownames(x = subset(x = markers.ident, subset = myAUC > auc.cutoff))
      )
    )
  }
  return(batch.genes)
}

#' ROC-based marker discovery
#'
#' Identifies 'markers' of gene expression using ROC analysis. For each gene,
#' evaluates (using AUC) a classifier built on that gene alone, to classify
#' between two groups of cells.
#'
#' An AUC value of 1 means that expression values for this gene alone can
#' perfectly classify the two groupings (i.e. Each of the cells in cells.1
#' exhibit a higher level than each of the cells in cells.2). An AUC value of 0
#' also means there is perfect classification, but in the other direction. A
#' value of 0.5 implies that the gene has no predictive power to classify the
#' two groups.
#'
#' @inheritParams FindMarkers
#' @inheritParams DiffExpTest
#' @param object Seurat object
#'
#' @return Returns a 'predictive power' (abs(AUC-0.5)) ranked matrix of
#' putative differentially expressed genes.
#'
#' @import ROCR
#'
#' @export
#'
MarkerTest <- function(
  object,
  cells.1,
  cells.2,
  genes.use = NULL,
  print.bar = TRUE
) {
  genes.use <- set.ifnull(x = genes.use, y = rownames(x = object@data))
  to.return <- marker.auc.test(
    data1 = object@data[, cells.1],
    data2 = object@data[, cells.2],
    mygenes = genes.use,
    print.bar = print.bar
  )
  to.return$power <- abs(x = to.return$myAUC - 0.5) * 2
  #print(head(to.return))
  return(to.return)
}

#' Differential expression testing using Student's t-test
#'
#' Identify differentially expressed genes between two groups of cells using
#' the Student's t-test
#'
#' @inheritParams FindMarkers
#' @inheritParams DiffExpTest
#'
#' @return Returns a p-value ranked matrix of putative differentially expressed
#' genes.
#'
#' @importFrom pbapply pblapply
#'
#' @export
#'
DiffTTest <- function(
  object,
  cells.1,
  cells.2,
  genes.use = NULL,
  print.bar = TRUE
) {
  genes.use <- set.ifnull(x = genes.use, y = rownames(x = object@data))
  data.use <- object@data
  if (print.bar) {
    iterate.fxn=pblapply
  } else {
    iterate.fxn <- lapply
  }
  p_val <- unlist(
    x = iterate.fxn(
      X = genes.use,
      FUN = function(x) {
        t.test(x = object@data[x, cells.1], y = object@data[x, cells.2])$p.value
      }
    )
  )
  to.return <- data.frame(p_val,row.names = genes.use)
  return(to.return)
}

#' FastWhichCells
#' Identify cells matching certain criteria (limited to character values)
#' @param object Seurat object
#' @param group.by Group cells in different ways (for example, orig.ident). Should be a column name in object@data.info
#' @param subset.value  Return cells matching this value
#' @param invert invert cells to return.FALSE by default
#'
#' @export
#'
FastWhichCells <- function(object, group.by, subset.value, invert = FALSE) {
  object <- SetAllIdent(object = object, id = group.by)
  cells.return <- WhichCells(object = object, ident = subset.value)
  if (invert) {
    cells.return <- setdiff(x = object@cell.names, y = cells.return)
  }
  return(cells.return)
}



#' Identify cells matching certain criteria
#'
#' Returns a list of cells that match a particular set of criteria such as
#' identity class, high/low values for particular PCs, ect..
#'
#' @param object Seurat object
#' @param ident Identity classes to subset. Default is all identities.
#' @param ident.remove Indentity classes to remove. Default is NULL.
#' @param cells.use Subset of cell names
#' @param subset.name Parameter to subset on. Eg, the name of a gene, PC1, a
#' column name in object@@data.info, etc. Any argument that can be retreived
#' using FetchData
#' @param accept.low Low cutoff for the parameter (default is -Inf)
#' @param accept.high High cutoff for the parameter (default is Inf)
#' @param accept.value Returns all cells with the subset name equal to this value
#' @param max.cells.per.ident Can be used to downsample the data to a certain max per cell ident. Default is inf.
#' @param random.seed Random seed for downsampling
#'
#' @return A vector of cell names
#'
#' @export
#'
WhichCells <- function(
  object,
  ident = NULL,
  ident.remove = NULL,
  cells.use = NULL,
  subset.name = NULL,
  accept.low = -Inf,
  accept.high = Inf,
  accept.value = NULL,
  max.cells.per.ident = Inf,
  random.seed = 1
) {
  set.seed(seed = random.seed)
  cells.use <- set.ifnull(x = cells.use, y = object@cell.names)
  ident <- set.ifnull(x = ident, y = unique(x = object@ident))
  ident <- anotinb(x = ident, y = ident.remove)
  if (! all(ident %in% unique(x = object@ident))) {
    bad.idents <- ident[! (ident %in% unique(x = object@ident))]
    stop(paste("Identity :", bad.idents, "not found.   "))
  }
  cells.to.use <- character()
  for (id in ident) {
    cells.in.ident <- object@ident[cells.use]
    cells.in.ident <- names(x = cells.in.ident[cells.in.ident == id])
    cells.in.ident <- cells.in.ident[! is.na(x = cells.in.ident)]
    if (length(x = cells.in.ident) > max.cells.per.ident) {
      cells.in.ident <- sample(x = cells.in.ident, size = max.cells.per.ident)
    }
    cells.to.use <- c(cells.to.use, cells.in.ident)
  }
  cells.use <- cells.to.use
  if (! is.null(x = subset.name)){
    data.use <- FetchData(
      object = object,
      vars.all = subset.name,
      cells.use = cells.use
    )
    if (length(x = data.use) == 0) {
      stop(paste("Error : ", id, " not found"))
    }
    subset.data <- data.use[, subset.name]
    if(! is.null(x = accept.value)) {
      pass.inds <- which(x = subset.data == accept.value)
    } else {
      pass.inds <- which(x = (subset.data > accept.low) & (subset.data < accept.high))
    }
    cells.use <- rownames(x = data.use)[pass.inds]
  }
  return(cells.use)
}

#' Switch identity class definition to another variable
#'
#' @param object Seurat object
#' @param id Variable to switch identity class to (for example, 'DBclust.ident', the output
#' of density clustering) Default is orig.ident - the original annotation pulled from the cell name.
#'
#' @return A Seurat object where object@@ident has been appropriately modified
#'
#' @export
#'
SetAllIdent <- function(object, id = NULL) {
  id <- set.ifnull(x = id, y = "orig.ident")
  if (id %in% colnames(x = object@data.info)) {
    cells.use <- rownames(x = object@data.info)
    ident.use <- object@data.info[, id]
    object <- SetIdent(
      object = object,
      cells.use = cells.use,
      ident.use = ident.use
    )
  }
  return(object)
}

#' Rename one identity class to another
#'
#' Can also be used to join identity classes together (for example, to merge clusters).
#'
#' @param object Seurat object
#' @param old.ident.name The old identity class (to be renamed)
#' @param new.ident.name The new name to apply
#'
#' @return A Seurat object where object@@ident has been appropriately modified
#'
#' @export
#'
RenameIdent <- function(object, old.ident.name = NULL, new.ident.name = NULL) {
  if (! old.ident.name %in% object@ident) {
    stop(paste("Error : ", old.ident.name, " is not a current identity class"))
  }
  new.levels <- old.levels <- levels(x = object@ident)
  # new.levels <- old.levels
  if (new.ident.name %in% old.levels) {
    new.levels <- new.levels[new.levels != old.ident.name]
  }
  if(! (new.ident.name %in% old.levels)) {
    new.levels[new.levels == old.ident.name] <- new.ident.name
  }
  ident.vector <- as.character(x = object@ident)
  names(x = ident.vector) <- names(object@ident)
  ident.vector[WhichCells(object = object, ident = old.ident.name)] <- new.ident.name
  object@ident <- factor(x = ident.vector, levels = new.levels)
  return(object)
}

#' Set identity class information
#'
#' Stashes the identity in data.info to be retrieved later. Useful if, for example, testing multiple clustering parameters
#'
#' @param object Seurat object
#' @param save.name Store current object@@ident under this column name in object@@data.info. Can be easily retrived with SetAllIdent
#'
#' @return A Seurat object where object@@ident has been appropriately modified
#'
#' @export
#'
StashIdent <- function(object, save.name = "oldIdent") {
  object@data.info[, save.name] <- as.character(x = object@ident)
  return(object)
}

#' Set identity class information
#'
#' Sets the identity class value for a subset (or all) cells
#'
#' @param object Seurat object
#' @param cells.use Vector of cells to set identity class info for (default is
#' all cells)
#' @param ident.use Vector of identity class values to assign (character
#' vector)
#'
#' @return A Seurat object where object@@ident has been appropriately modified
#'
#' @importFrom gdata drop.levels
#'
#' @export
#'
SetIdent <- function(object, cells.use = NULL, ident.use = NULL) {
  cells.use <- set.ifnull(x = cells.use, y = object@cell.names)
  if (length(x = anotinb(x = cells.use, y = object@cell.names) > 0)) {
    stop(paste(
      "ERROR : Cannot find cells ",
      anotinb(x = cells.use, y = object@cell.names)
    ))
  }
  ident.new <- anotinb(x = ident.use, y = levels(x = object@ident))
  object@ident <- factor(
    x = object@ident,
    levels = unique(
      x = c(
        as.character(x = object@ident),
        as.character(x = ident.new)
      )
    )
  )
  object@ident[cells.use] <- ident.use
  object@ident <- drop.levels(x = object@ident)
  return(object)
}

#' Merge subchilden of a node
#'
#' Merge the subchilden of a node into a single identity class
#'
#' @param object Seurat object
#' @param node.use Merge subchildren of this node
#'
#' @export
#'
MergeNode <- function(object, node.use = NULL) {
  object.tree <- object@cluster.tree[[1]]
  node.children <- DFT(
    tree = object.tree,
    node = node.use,
    include.children = TRUE
  )
  node.children <- ainb(a = node.children, b = levels(x = object@ident))
  children.cells <- WhichCells(object = object, ident = node.children)
  if (length(x = children.cells > 0)) {
    object <- SetIdent(
      object = object,
      cells.use = children.cells,
      ident.use = min(node.children)
    )
  }
  return(object)
}

# Documentation
###############
#Not documented for now
#' @export
PosteriorPlot <- function(object, name) {
  post.names <- colnames(x = subc(data = object@mix.probs, code = name))
  VlnPlot(
    object = object,
    features.plot = post.names,
    inc.first=TRUE,
    inc.final=TRUE,
    by.k=TRUE
  )
}

# Documentation
###############
#Internal, not documented for now
map.cell.score <- function(gene, gene.value, insitu.bin, mu, sigma, alpha) {
  code.1 <- paste(gene, insitu.bin, sep=".")
  mu.use <- mu[paste(code.1, "mu", sep="."), 1]
  sigma.use <- sigma[paste(code.1, "sigma", sep="."), 1]
  alpha.use <- alpha[paste(code.1, "alpha", sep="."), 1]
  bin.prob <- unlist(
    x = lapply(
      X = 1:length(x = insitu.bin),
      FUN = function(x) {
        return(dnorm(
          x = gene.value,
          mean = mu.use[x],
          sd = sigma.use[x],
          log = TRUE) + log(x = alpha.use[x])
        )
      }
    )
  )
  return(bin.prob)
}

#Internal, not documented for now
#' @export
#'
MapCell <- function(
  object,
  cell.name,
  do.plot = FALSE,
  safe.use = TRUE,
  text.val = NULL,
  do.rev = FALSE
) {
  insitu.matrix <- object@insitu.matrix
  insitu.genes <- colnames(x = insitu.matrix)
  insitu.genes <- insitu.genes[insitu.genes %in% rownames(x = object@imputed)]
  insitu.use <- insitu.matrix[, insitu.genes]
  imputed.use <- object@imputed[insitu.genes, ]
  if (safe.use) {
    safe_fxn <- log_add
  } else {
    safe_fxn <- sum
  }
  all.needed.cols <- unique(
    x = unlist(
      x = lapply(
        X = insitu.genes,
        FUN = function(x) {
          return(paste(x, insitu.use[, x], "post", sep="."))
        }
      )
    )
  )
  missing.cols <- which(! (all.needed.cols %in% colnames(x = object@mix.probs)))
  if (length(x = missing.cols) > 0) {
    stop(paste(
      "Error: ",
      all.needed.cols[missing.cols],
      "is missing from the mixture fits"
    ))
  }
  all.probs <- data.frame(
    sapply(
      X = insitu.genes,
      FUN = function(x) {
        return(
          log(x = as.numeric(x = object@mix.probs[
            cell.name, # Row
            paste(x, insitu.use[, x], "post", sep=".") # Column
          ])))
      }
    )
  )
  scale.probs <- t(
    x = t(x = all.probs) - apply(X = t(x = all.probs), MARGIN = 1, FUN = log_add)
  )
  scale.probs[scale.probs < (-9.2)] <- (-9.2)
  #head(scale.probs)
  total.prob <- exp(x = apply(X = scale.probs, MARGIN = 1, FUN = safe_fxn))
  total.prob <- total.prob / sum(total.prob)
  if (do.plot) {
    #plot(total.prob,main=cell.name)
    par(mfrow = c(1, 2))
    txt.matrix <- matrix(data = rep(x = "", 64), nrow=8, ncol=8)
    if (! is.null(x = text.val)) {
      txt.matrix[text.val] <- "X"
    }
    if (do.rev) {
      scale.probs <- scale.probs[unlist(
        x = lapply(
          X = 0:7,
          FUN = function(x) {
            return(seq(from = 1, to = 57, by = 8) + x)
          }
        )
      ), ]
    }
    aheatmap(
      x = matrix(data = total.prob, nrow=8, ncol=8),
      Rowv = NA,
      Colv = NA,
      txt = txt.matrix,
      col = bwCols
    )
    aheatmap(x = scale.probs, Rowv = NA, Colv = NA)
    rp()
  }
  return(total.prob)
}

#' Get cell centroids
#'
#' Calculate the spatial mapping centroids for each cell, based on previously
#' calculated mapping probabilities for each bin.
#'
#' Currently, Seurat assumes that the tissue of interest has an 8x8 bin
#' structure. This will be broadened in a future release.
#'
#' @param object Seurat object
#' @param cells.use Cells to calculate centroids for (default is all cells)
#' @param get.exact Get exact centroid (Default is TRUE). If FALSE, identify
#' the single closest bin.
#'
#' @return Data frame containing the x and y coordinates for each cell
#' centroid.
#'
#' @export
#'
GetCentroids <- function(object, cells.use = NULL, get.exact = TRUE) {
  cells.use <- set.ifnull(x = cells.use, y = colnames(x = object@final.prob))
  #Error checking
  cell.names <- ainb(a = cells.use, b = colnames(x = object@final.prob))
  if (length(x = cell.names) != length(x = cells.use)) {
    print(paste(
      "Error",
      anotinb(x = cells.use, y = colnames(x = object@final.prob)),
      " have not been mapped"
    ))
    return(0)
  }
  if (get.exact) {
    my.centroids <- data.frame(t(x = sapply(
      X = colnames(x = object@data),
      FUN = function(x) {
        return(exact.cell.centroid(cell.probs = object@final.prob[, x]))
      }
    )))
  } else {
    my.centroids <- data.frame(t(x = sapply(
      X = colnames(x = object@data),
      FUN = function(x) {
        return(cell.centroid(cell.probs = object@final.prob[, x]))
      }
    )))
  }
  colnames(x = my.centroids) <- c("bin.x","bin.y")
  return(my.centroids)
}

#' Quantitative refinement of spatial inferences
#'
#' Refines the initial mapping with more complex models that allow gene
#' expression to vary quantitatively across bins (instead of 'on' or 'off'),
#' and that also considers the covariance structure between genes.
#'
#' Full details given in spatial mapping manuscript.
#'
#' @param object Seurat object
#' @param genes.use Genes to use to drive the refinement procedure.
#'
#' @return Seurat object, where mapping probabilities for each bin are stored
#' in object@@final.prob
#'
#' @import fpc
#'
#' @export
#'
RefinedMapping <- function(object, genes.use) {
  genes.use <- ainb(a = genes.use, b = rownames(x = object@imputed))
  cells.max <- t(x = sapply(
    X = colnames(object@data),
    FUN = function(x) {
      return(exact.cell.centroid(object@final.prob[, x]))
    }
  ))
  all.mu <- sapply(
    X = genes.use,
    FUN = function(gene) {
      return(sapply(X = 1:64, FUN = function(bin) {
        mean(x = as.numeric(x = object@imputed[
          gene, # Row
          fetch.closest(
            bin = bin,
            all.centroids = cells.max,
            num.cell = 2*length(x = genes.use)
          ) # Column
        ]))
      }))
    }
  )
  all.cov <- list()
  for (x in 1:64) {
    all.cov[[x]] <- cov(
      x = t(
        x = object@imputed[
          genes.use, # Row
          fetch.closest(
            bin = x,
            all.centroids = cells.max,
            num.cell = 2*length(x = genes.use)
          ) # Columns
        ]
      )
    )
  }
  mv.probs <- sapply(
    X = colnames(x = object@data),
    FUN = function(my.cell) {
      return(sapply(X = 1:64, FUN = function(bin) {
        return(slimdmvnorm(
          x = as.numeric(x = object@imputed[genes.use, my.cell]),
          mean = as.numeric(x = all.mu[bin, genes.use]),
          sigma = all.cov[[bin]])
        )
      }))
    }
  )
  mv.final <- exp(
    x = sweep(
      x = mv.probs,
      MARGIN = 2,
      STATS = apply(X = mv.probs, MARGIN = 2, FUN = log_add)
    )
  )
  object@final.prob <- data.frame(mv.final)
  return(object)
}

#' Infer spatial origins for single cells
#'
#' Probabilistically maps single cells based on (imputed) gene expression
#' estimates, a set of mixture models, and an in situ spatial reference map.
#'
#' @param object Seurat object
#' @param cells.use Which cells to map
#'
#' @return Seurat object, where mapping probabilities for each bin are stored
#' in object@@final.prob
#'
#' @export
#'
InitialMapping <- function(object, cells.use = NULL) {
  cells.use <- set.ifnull(x = cells.use, y = colnames(x = object@data))
  every.prob <- sapply(
    X = cells.use,
    FUN = function(x) {
      return(MapCell(
        object = object,
        cell.name = x,
        do.plot = FALSE,
        safe.use = FALSE
      ))
    }
  )
  object@final.prob <- data.frame(every.prob)
  rownames(x = object@final.prob) <- paste0("bin.", rownames(x = object@final.prob))
  return(object)
}

# Documentation
###############
#Internal, not documented for now
CalcInsitu <- function(
  object,
  gene,
  do.plot = TRUE,
  do.write = FALSE,
  write.dir = "~/window/insitu/",
  col.use = bwCols,
  do.norm = FALSE,
  cells.use = NULL,
  do.return = FALSE,
  probs.min = 0,
  do.log = FALSE,
  use.imputed = FALSE,
  bleach.use = 0
) {
  cells.use <- set.ifnull(x = cells.use, y = colnames(x = object@final.prob))
  probs.use <- object@final.prob
  if (use.imputed) {
    data.use <- exp(x = object@imputed) - 1
  } else {
    data.use <- exp(object@data) - 1
  }
  cells.use <- cells.use[cells.use %in% colnames(x = probs.use)]
  cells.use <- cells.use[cells.use %in% colnames(x = data.use)]
  #insilico.stain=matrix(unlist(lapply(1:64,function(x) sum(probs.use[x,]*data.use[gene,]))),nrow=8,ncol=8)
  insilico.vector <- unlist(
    x = lapply(
      X = 1:64,
      FUN = function(x) {
        return(sum(
          as.numeric(x = probs.use[x, cells.use]) *
            as.numeric(x = data.use[gene, cells.use])
        ))
      }
    )
  )
  probs.total <- apply(X = probs.use, MARGIN = 1, FUN = sum)
  probs.total[probs.total < probs.min] <- probs.min
  insilico.stain <- (matrix(data = insilico.vector / probs.total, nrow=8, ncol=8))
  if (do.log) {
    insilico.stain <- log(x = insilico.stain + 1)
  }
  if (bleach.use > 0) {
    insilico.stain <- insilico.stain - bleach.use
    insilico.stain <- minmax(data = insilico.stain, min=0, max=1e6)
  }
  if (do.norm) {
    insilico.stain <- (insilico.stain - min(insilico.stain)) /
      (max(insilico.stain) - min(insilico.stain))
  }
  title.use <- gene
  if (gene %in% colnames(x = object@insitu.matrix)) {
    pred.use <- prediction(
      predictions = insilico.vector / probs.total,
      labels = object@insitu.matrix[, gene],
      label.ordering = 0:1
    )
    perf.use <- performance(prediction.obj = pred.use, measure = "auc")
    auc.use <- round(x = perf.use@y.values[[1]], digits = 3)
    title.use <- paste(gene, sep=" ")
  }
  if (do.write) {
    write.table(
      x = insilico.stain,
      file = paste0(write.dir, gene, ".txt"),
      quote=FALSE,
      row.names=FALSE,
      col.names=FALSE
    )
  }
  if (do.plot) {
    aheatmap(
      x = insilico.stain,
      Rowv = NA,
      Colv = NA,
      col = col.use,
      main=title.use
    )
  }
  if (do.return) {
    return(as.vector(x = insilico.stain))
  }
  return(object)
}

#' Build mixture models of gene expression
#'
#' Models the imputed gene expression values as a mixture of gaussian
#' distributions. For a two-state model, estimates the probability that a given
#' cell is in the 'on' or 'off' state for any gene. Followed by a greedy
#' k-means step where cells are allowed to flip states based on the overall
#' structure of the data (see Manuscript for details)
#'
#' @param object Seurat object
#' @param gene Gene to fit
#' @param do.k Number of modes for the mixture model (default is 2)
#' @param num.iter Number of 'greedy k-means' iterations (default is 1)
#' @param do.plot Plot mixture model results
#' @param genes.use Genes to use in the greedy k-means step (See manuscript for details)
#' @param start.pct Initial estimates of the percentage of cells in the 'on'
#' state (usually estimated from the in situ map)
#'
#' @return A Seurat object, where the posterior of each cell being in the 'on'
#' or 'off' state for each gene is stored in object@@mix.probs
#'
#' @importFrom mixtools normalmixEM
#'
#' @export
#'
FitGeneK <- function(
  object,
  gene,
  do.k = 2,
  num.iter = 1,
  do.plot = FALSE,
  genes.use = NULL,
  start.pct = NULL
) {
  data <- object@imputed
  data.use <- data[gene, ]
  names(x = data.use) <- colnames(x = data.use)
  scale.data <- t(x = scale(x = t(x = object@imputed)))
  genes.use <- set.ifnull(x = genes.use,y = rownames(x = scale.data))
  genes.use <- genes.use[genes.use %in% rownames(x = scale.data)]
  scale.data <- scale.data[genes.use, ]
  data.cut <- as.numeric(x = data.use[gene, ])
  cell.ident <- as.numeric(x = cut(x = data.cut, breaks = do.k))
  if (! (is.null(x = start.pct))) {
    cell.ident <- rep.int(x = 1, times = length(x = data.cut))
    cell.ident[data.cut > quantile(x = data.cut, probs = 1 - start.pct)] <- 2
  }
  cell.ident <- order(tapply(
    X = as.numeric(x = data.use),
    INDEX = cell.ident,
    FUN = mean
  ))[cell.ident]
  ident.table <- table(cell.ident)
  if (num.iter > 0) {
    for (i2 in 1:num.iter) {
      cell.ident <- iter.k.fit(
        scale.data = scale.data,
        cell.ident = cell.ident,
        data.use = data.use
      )
      ident.table <- table(cell.ident)
    }
  }
  ident.table <- table(cell.ident)
  raw.probs <- t(
    x = sapply(
      X = data.use,
      FUN = function(y) {
        return(unlist(
          x = lapply(
            X = 1:do.k,
            FUN = function(x) {
              return(
                (ident.table[x] / sum(ident.table)) * dnorm(
                  x = y,
                  mean = mean(x = as.numeric(x = data.use[cell.ident == x])),
                  sd = sd(x = as.numeric(x = data.use[cell.ident == x]))
                )
              )
            }
          )
        ))
      }
    )
  )
  norm.probs <- raw.probs / apply(X = raw.probs, MARGIN = 1, FUN = sum)
  colnames(x = norm.probs) <- unlist(
    x = lapply(
      X = 1:do.k,
      FUN = function(x) {
        paste(gene, x - 1, "post", sep=".")
      }
    )
  )
  norm.probs <- cbind(norm.probs, cell.ident)
  colnames(x = norm.probs)[ncol(x = norm.probs)] <- paste0(gene, ".ident")
  new.mix.probs <- data.frame(
    minusc(data = object@mix.probs, code = paste0(gene, ".")),
    row.names = rownames(x = object@mix.probs)
  )
  colnames(x = new.mix.probs)[1] <- "nGene"
  object@mix.probs <- cbind(new.mix.probs, norm.probs)
  if (do.plot) {
    nCol <- 2
    num.row <- floor(x = (do.k + 1) / nCol - (1e-5)) + 1
    hist(
      x = as.numeric(x = data.use),
      probability = TRUE,
      ylim = c(0, 1),
      xlab = gene,
      main = gene
    )
    for (i in 1:do.k) {
      lines(
        x = seq(from = -10, to = 10, by = 0.01),
        y = (ident.table[i] / sum(ident.table)) * dnorm(
          x = seq(from = -10, to = 10, by = 0.01),
          mean = mean(x = as.numeric(x = data.use[cell.ident == i])),
          sd = sd(x = as.numeric(x = data.use[cell.ident == i]))
        ),
        col=i,
        lwd=2
      )
    }
  }
  return(object)
}

#Internal, not documented for now
iter.k.fit <- function(scale.data, cell.ident, data.use) {
  means.all <- sapply(
    X = sort(x = unique(x = cell.ident)),
    FUN = function(x) {
      return(apply(X = scale.data[, cell.ident == x], MARGIN = 1, FUN = mean))
    }
  )
  all.dist <- data.frame(
    t(x = sapply(
      X = 1:ncol(x = scale.data),
      FUN = function(x) {
        return(unlist(x = lapply(
          X = sort(x = unique(x = cell.ident)),
          FUN = function(y) {
            return(dist(x = rbind(scale.data[, x], means.all[, y])))
          }
        )))
      }
    ))
  )
  cell.ident <- apply(X = all.dist, MARGIN = 1, FUN = which.min)
  cell.ident <- order(tapply(
    X = as.numeric(x = data.use),
    INDEX = cell.ident,
    FUN = mean
  ))[cell.ident]
  return(cell.ident)
}

# Documentation
###############
#Internal, not documented for now
#' @export
FitGeneMix <- function(
  object,
  gene,
  do.k = 3,
  use.mixtools = TRUE,
  do.plot = FALSE,
  plot.with.imputed = TRUE,
  min.bin.size = 10
) {
  data.fit <- as.numeric(x = object@imputed[gene, ])
  mixtools.fit <- normalmixEM(x = data.fit, k = do.k)
  comp.order <- order(mixtools.fit$mu)
  mixtools.posterior <- data.frame(mixtools.fit$posterior[, comp.order])
  colnames(x = mixtools.posterior) <- unlist(
    x = lapply(
      X = 1:do.k,
      FUN = function(x) {
        return(paste(gene, x - 1, "post", sep="."))
      }
    )
  )
  #mixtools.mu=data.frame(mixtools.fit$mu[comp.order])
  #mixtools.sigma=data.frame(mixtools.fit$sigma[comp.order])
  #mixtools.alpha=data.frame(mixtools.fit$lambda[comp.order])
  #rownames(mixtools.mu)=unlist(lapply(1:do.k,function(x)paste(gene,x-1,"mu",sep=".")))
  #rownames(mixtools.sigma)=unlist(lapply(1:do.k,function(x)paste(gene,x-1,"sigma",sep=".")))
  #rownames(mixtools.alpha)=unlist(lapply(1:do.k,function(x)paste(gene,x-1,"alpha",sep=".")))
  #object@mix.mu = rbind(minusr(object@mix.mu,gene), mixtools.mu);
  #object@mix.sigma = rbind(minusr(object@mix.sigma,gene), mixtools.sigma);
  #o#bject@mu.alpha =rbind(minusr(object@mu.alpha,gene), mixtools.alpha);
  if (do.plot) {
    nCol <- 2
    num.row <- floor(x = (do.k + 1) / nCol - (1e-5)) + 1
    par(mfrow = c(num.row, nCol))
    plot.mixEM(x = mixtools.fit, which = 2)
    plot.data <- as.numeric(x = object@imputed[gene, ])
    if (! plot.with.imputed) {
      plot.data <- as.numeric(x = object@data[gene, ])
    }
    unlist(
      x = lapply(
        X = 1:do.k,
        FUN = function(x) {
          plot(
            x = plot.data,
            y = mixtools.posterior[, x],
            ylab = paste0("Posterior for Component ", x - 1),
            xlab = gene,
            main = gene
          )
        }
      )
    )
  }
  new.mix.probs <- data.frame(
    minusc(data = object@mix.probs, code = paste0(gene, ".")),
    row.names = rownames(x = object@mix.probs)
  )
  colnames(x = new.mix.probs)[1] <- "nGene"
  object@mix.probs <- cbind(new.mix.probs, mixtools.posterior)
  return(object)
}

# Documentation
###############
#Internal, not documented for now
lasso.fxn <- function(
  lasso.input,
  genes.obs,
  s.use = 20,
  gene.name = NULL,
  do.print = FALSE,
  gram = TRUE
) {
  lasso.model <- lars(
    x = lasso.input,
    y = as.numeric(x = genes.obs),
    type = "lasso",
    max.steps = s.use * 2,
    use.Gram = gram
  )
  #lasso.fits=predict.lars(lasso.model,lasso.input,type="fit",s=min(s.use,max(lasso.model$df)))$fit
  lasso.fits <- predict.lars(
    object = lasso.model,
    newx = lasso.input,
    type = "fit",
    s = s.use
  )$fit
  if (do.print) {
    print(gene.name)
  }
  return(lasso.fits)
}

#' Calculate smoothed expression values
#'
#'
#' Smooths expression values across the k-nearest neighbors based on dimensional reduction
#'
#' @inheritParams FeaturePlot
#' @inheritParams AddImputedScore
#' @param genes.fit Genes to calculate smoothed values for
#' @param k k-param for k-nearest neighbor calculation
#' @param do.log Whether to perform smoothing in log space. Default is false.
#'
#' @importFrom FNN get.knn
#'
#' @export
#'
AddSmoothedScore <- function(
  object,
  genes.fit = NULL,
  dim.1 = 1,
  dim.2 = 2,
  reduction.use = "tSNE",
  k = 30,
  do.log = FALSE,
  do.print = FALSE
) {
  genes.fit <- set.ifnull(x = genes.fit, y = object@var.genes)
  genes.fit <- genes.fit[genes.fit %in% rownames(x = object@data)]
  dim.code <- translate.dim.code(object = object, reduction.use = reduction.use)
  dim.codes <- paste0(dim.code, c(dim.1, dim.2))
  data.plot <- FetchData(object = object, vars.all = dim.codes)
  knn.smooth <- get.knn(data = data.plot, k = k)$nn.index
  avg.fxn <- mean
  if (! do.log) {
    avg.fxn <- expMean
  }
  lasso.fits <- data.frame(
    t(
      x = sapply(
        X = genes.fit,
        FUN = function(g) {
          return(unlist(
            x = lapply(
              X = 1:nrow(x = data.plot),
              FUN = function(y) {
                avg.fxn(as.numeric(x = object@data[g, knn.smooth[y, ]]))
              }
            )
          ))
        }
      )
    )
  )
  colnames(x = lasso.fits) <- rownames(x = data.plot)
  genes.old <- genes.fit[genes.fit %in% rownames(x = object@imputed)]
  genes.new <- genes.fit[! (genes.fit %in% rownames(x = object@imputed))]
  if (length(x = genes.old) > 0) {
    object@imputed[genes.old, ] <- lasso.fits[genes.old, ]
  }
  object@imputed <- rbind(object@imputed, lasso.fits[genes.new, ])
  return(object)
}

#' Calculate imputed expression values
#'
#' Uses L1-constrained linear models (LASSO) to impute single cell gene
#' expression values.
#'
#'
#' @param object Seurat object
#' @param genes.use A vector of genes (predictors) that can be used for
#' building the LASSO models.
#' @param genes.fit A vector of genes to impute values for
#' @param s.use Maximum number of steps taken by the algorithm (lower values
#' indicate a greater degree of smoothing)
#' @param do.print Print progress (output the name of each gene after it has
#' been imputed).
#' @param gram The use.gram argument passed to lars
#'
#' @return Returns a Seurat object where the imputed values have been added to
#' object@@data
#'
#' @import lars
#'
#' @export
#'
AddImputedScore <- function(
  object,
  genes.use = NULL,
  genes.fit = NULL,
  s.use = 20,
  do.print = FALSE,
  gram = TRUE
) {
  genes.use <- set.ifnull(x = genes.use, y = object@var.genes)
  genes.fit <- set.ifnull(x = genes.fit, y = object@var.genes)
  genes.use <- genes.use[genes.use %in% rownames(x = object@data)]
  genes.fit <- genes.fit[genes.fit %in% rownames(x = object@data)]
  lasso.input <- t(x = object@data[genes.use, ])
  lasso.fits <- data.frame(t(
    x = sapply(
      X = genes.fit,
      FUN = function(x) {
        return(
          lasso.fxn(
            lasso.input = t(x = object@data[genes.use[genes.use != x], ]),
            genes.obs = object@data[x, ],
            s.use = s.use,
            gene.name = x,
            do.print = do.print,
            gram = gram
          )
        )
      }
    )
  ))
  genes.old <- genes.fit[genes.fit %in% rownames(x = object@imputed)]
  genes.new <- genes.fit[! (genes.fit %in% rownames(x = object@imputed))]
  if (length(x = genes.old) > 0) {
    object@imputed[genes.old, ] <- lasso.fits[genes.old, ]
  }
  object@imputed <- rbind(object@imputed, lasso.fits[genes.new, ])
  return(object)
}

# Documentation
###############
# Not currently supported, but a cool scoring function
#' @export
GetNewScore <- function(
  object,
  score.name,
  score.genes,
  cell.ids = NULL,
  score.func = weighted.mean,
  scramble = FALSE,
  no.tech.wt = FALSE,
  biol.wts = NULL,
  use.scaled = FALSE
) {
  data.use <- object@data
  if (use.scaled) {
    data.use <- minmax(data = object@scale.data, min = -2, max=2)
  }
  if (no.tech.wt) {
    score.genes <- score.genes[score.genes %in% rownames(x = data.use)]
  } else {
    score.genes <- score.genes[score.genes %in% rownames(x = object@wt.matrix)]
  }
  cell.ids <- set.ifnull(x = cell.ids, y = colnames(x = data.use))
  wt.matrix <- data.frame(
    matrix(
      data = 1,
      nrow = length(x = score.genes),
      ncol = length(x = cell.ids),
      dimnames = list(score.genes, cell.ids)
    )
  )
  if (! (no.tech.wt)) {
    wt.matrix <- object@drop.wt.matrix[score.genes, cell.ids]
  }
  score.data <- data.use[score.genes, cell.ids]
  if (scramble) {
    score.data <- score.data[, sample(x = ncol(x = score.data))]
  }
  wt.matrix <- wt.matrix * (wt.matrix)
  if (no.tech.wt) {
    wt.matrix[wt.matrix < 1] = 1
  }
  biol.wts <- set.ifnull(x = biol.wts, y = rep(x = 1, nrow(x = wt.matrix)))
  if (mean(x = biol.wts) == 1) {
    names(x = biol.wts) <- score.genes
  }
  my.scores <- unlist(
    x = lapply(
      X = colnames(x = score.data),
      FUN = function(x) {
        return(score.func(score.data[, x], wt.matrix[, x] * biol.wts[score.genes]))
      }
    )
  )
  names(x = my.scores) <- colnames(x = score.data)
  object@gene.scores[cell.ids, score.name] <- my.scores
  return(object)
}

# Documentation
###############
# Not currently supported, but a cool function for QC
#' @export
CalcNoiseModels <- function(
  object,
  cell.ids = NULL,
  trusted.genes = NULL,
  n.bin = 20,
  drop.expr = 1
) {
  object@drop.expr <- drop.expr
  cell.ids <- set.ifnull(x = cell.ids, y = 1:ncol(x = object@data))
  trusted.genes <- set.ifnull(x = trusted.genes, y = rownames(x = object@data))
  trusted.genes <- trusted.genes[trusted.genes %in% rownames(x = object@data)]
  object@trusted.genes <- trusted.genes
  data <- object@data[trusted.genes, ]
  idents <- data.frame(data[, 1])
  code_humpAvg <- apply(X = data, MARGIN = 1, FUN = humpMean, min = object@drop.expr)
  code_humpAvg[code_humpAvg > 9] <- 9
  code_humpAvg[is.na(x = code_humpAvg)] <- 0
  idents$code_humpAvg <- code_humpAvg
  data[data > object@drop.expr] <- 1
  data[data < object@drop.expr] <- 0
  data$bin <- cut(x = code_humpAvg, breaks = n.bin)
  data$avg <- code_humpAvg
  rownames(x = idents) <- rownames(x = data)
  my.coefs <- data.frame(t(x = pbsapply(
    X = colnames(x = data[1:(ncol(x = data) - 2)]),
    FUN = getAB,
    data = data,
    data2 = idents,
    status = "code",
    code2 = "humpAvg",
    hasBin = TRUE,
    doPlot = FALSE
  )))
  colnames(x = my.coefs) <- c("a", "b")
  object@drop.coefs <- my.coefs
  return(object)
}

#' Visualize 'features' on a dimensional reduction plot
#'
#' Colors single cells on a dimensional reduction plot according to a 'feature'
#' (i.e. gene expression, PC scores, number of genes detected, etc.)
#'
#'
#' @param object Seurat object
#' @param features.plot Vector of features to plot
#' @param min.cutoff Vector of minimum cutoff values for each feature, may specify quantile in the form of 'q##' where '##' is the quantile (eg, 1, 10)
#' @param max.cutoff Vector of maximum cutoff values for each feature, may specify quantile in the form of 'q##' where '##' is the quantile (eg, 1, 10)
#' @param dim.1 Dimension for x-axis (default 1)
#' @param dim.2 Dimension for y-axis (default 2)
#' @param cells.use Vector of cells to plot (default is all cells)
#' @param pt.size Adjust point size for plotting
#' @param cols.use The two colors to form the gradient over. Provide as string vector with
#' the first color corresponding to low values, the second to high. Also accepts a Brewer
#' color scale or vector of colors. Note: this will bin the data into number of colors provided.
#' @param pch.use Pch for plotting
#' @param overlay Plot two features overlayed one on top of the other
#' @param do.hover Enable hovering over points to view information
#' @param data.hover Data to add to the hover, pass a character vector of features to add. Defaults to cell name
#' @param do.identify Opens a locator session to identify clusters of cells
#' @param reduction.use Which dimensionality reduction to use. Default is
#' "tsne", can also be "pca", or "ica", assuming these are precomputed.
#' @param use.imputed Use imputed values for gene expression (default is FALSE)
#' @param nCol Number of columns to use when plotting multiple features.
#' @param no.axes Remove axis labels
#' @param no.legend Remove legend from the graph. Default is TRUE.
#' @param dark.theme Plot in a dark theme
#' @param do.return return the ggplot2 object
#'
#' @importFrom RColorBrewer brewer.pal.info
#'
#' @return No return value, only a graphical output
#'
#' @export
#'
FeaturePlot <- function(
    object,
    features.plot,
    min.cutoff = NA,
    max.cutoff = NA,
    dim.1 = 1,
    dim.2 = 2,
    cells.use = NULL,
    pt.size = 1,
    cols.use = c("yellow", "red"),
    pch.use = 16,
    overlay = FALSE,
    do.hover = FALSE,
    data.hover = NULL,
    do.identify = FALSE,
    reduction.use = "tsne",
    use.imputed = FALSE,
    nCol = NULL,
    no.axes = FALSE,
    no.legend = TRUE,
    dark.theme = FALSE,
    do.return = FALSE
) {
    cells.use <- set.ifnull(cells.use, colnames(object@data))
    if (is.null(nCol)) {
      nCol <- 2
      if (length(features.plot) == 1) nCol <- 1
      if (length(features.plot) > 6) nCol <- 3
      if (length(features.plot) > 9) nCol <- 4
    }
    num.row <- floor(length(features.plot) / nCol - 1e-5) + 1
    if (overlay | do.hover) {
        num.row <- 1
        nCol <- 1
    }
    par(mfrow = c(num.row, nCol))
    dim.code <- translate.dim.code(object, reduction.use)
    dim.codes <- paste(dim.code, c(dim.1, dim.2), sep = "")
    data.plot <- FetchData(object, dim.codes, cells.use = cells.use)

    x1 <- paste(dim.code, dim.1, sep = "")
    x2 <- paste(dim.code, dim.2, sep = "")

    data.plot$x <- data.plot[, x1]
    data.plot$y <- data.plot[, x2]
    data.plot$pt.size <- pt.size
    names(x = data.plot) <- c('x', 'y')
    # data.plot$pt.size <- pt.size
    data.use <- t(FetchData(object, features.plot, cells.use = cells.use,
                                       use.imputed = use.imputed))
    #   Check mins and maxes
    min.cutoff <- mapply(
        FUN = function(cutoff, feature) {
            ifelse(
                test = is.na(x = cutoff),
                yes = min(data.use[feature, ]),
                no = cutoff
            )
        },
        cutoff = min.cutoff,
        feature = features.plot
    )
    max.cutoff <- mapply(
        FUN = function(cutoff, feature) {
            ifelse(
                test = is.na(x = cutoff),
                yes = max(data.use[feature, ]),
                no = cutoff
            )
        },
        cutoff = max.cutoff,
        feature = features.plot
    )
    check_lengths = unique(x = vapply(X = list(features.plot, min.cutoff, max.cutoff), FUN = length, FUN.VALUE = 1))
    if (length(x = check_lengths) != 1) {
        stop('There must be the same number of minimum and maximum cuttoffs as there are features')
    }
    if (overlay) {
        #   Wrap as a list for MutiPlotList
        pList <- list(
            BlendPlot(
                data.use = data.use,
                features.plot = features.plot,
                data.plot = data.plot,
                pt.size = pt.size,
                pch.use = pch.use,
                cols.use = cols.use,
                dim.codes = dim.codes,
                min.cutoff = min.cutoff,
                max.cutoff = max.cutoff,
                no.axes = no.axes,
                no.legend = no.legend,
                dark.theme = dark.theme
            )
        )
    } else {
        #   Use mapply instead of lapply for multiple iterative variables.
        pList <- mapply(
            FUN = SingleFeaturePlot,
            feature = features.plot,
            min.cutoff = min.cutoff,
            max.cutoff = max.cutoff,
            MoreArgs = list( # Arguments that are not being repeated
                data.use = data.use,
                data.plot = data.plot,
                pt.size = pt.size,
                pch.use = pch.use,
                cols.use = cols.use,
                dim.codes = dim.codes,
                no.axes = no.axes,
                no.legend = no.legend,
                dark.theme = dark.theme
            ),
            SIMPLIFY = FALSE # Get list, not matrix
        )
    }
    if (do.hover) {
        if (length(x = pList) != 1) {
            stop("'do.hover' only works on a single feature or an overlayed FeaturePlot")
        }
        if (is.null(x = data.hover)) {
            features.info <- NULL
        } else {
            features.info <- FetchData(object = object, vars.all = data.hover)
        }
        #   Use pList[[1]] to properly extract the ggplot out of the plot list
        return(HoverLocator(plot = pList[[1]], data.plot = data.plot, features.info = features.info, dark.theme = dark.theme, title = features.plot))
        # invisible(readline(prompt = 'Press <Enter> to continue\n'))
    } else if (do.identify) {
        if (length(x = pList) != 1) {
            stop("'do.identify' only works on a single feature or an overlayed FeaturePlot")
        }
        #   Use pList[[1]] to properly extract the ggplot out of the plot list
        return(FeatureLocator(plot = pList[[1]], data.plot = data.plot, dark.theme = dark.theme))
    } else {
        MultiPlotList(pList, cols = nCol)
    }
    rp()
    if(do.return){
      return(pList)
    }
}

# Documentation
###############
SingleFeaturePlot <- function(
    data.use,
    feature,
    data.plot,
    pt.size,
    pch.use,
    cols.use,
    dim.codes,
    min.cutoff,
    max.cutoff,
    no.axes,
    no.legend,
    dark.theme
) {
    data.gene <- na.omit(object = data.frame(data.use[feature, ]))
    #   Check for quantiles
    min.cutoff <- setQuantile(cutoff = min.cutoff, data = data.gene)
    max.cutoff <- setQuantile(cutoff = max.cutoff, data = data.gene)
    #   Mask any values below the minimum and above the maximum values
    data.gene <- sapply(X = data.gene, FUN = function(x) ifelse(test = x < min.cutoff, yes = min.cutoff, no = x))
    data.gene <- sapply(X = data.gene, FUN = function(x) ifelse(test = x > max.cutoff, yes = max.cutoff, no = x))
    data.plot$gene <- data.gene
    #   Stuff for break points
    if(length(cols.use) == 1){
        brewer.gran <- brewer.pal.info[cols.use, ]$maxcolors
    }
    else{
        brewer.gran <- length(cols.use)
    }
    #   Cut points
    if(all(data.gene == 0)){
        data.cut <- 0
    }
    else{
        data.cut <- as.numeric(x = as.factor(x = cut(x = as.numeric(x = data.gene), breaks = brewer.gran)))
    }
    data.plot$col <- as.factor(x = data.cut)
    #   Start plotting
    p <- ggplot(data.plot, aes(x, y))
    if(brewer.gran != 2){
        if(length(cols.use) == 1){
            p <- p + geom_point(aes(color=col), size=pt.size, shape=pch.use) +
                scale_color_brewer(palette=cols.use)
        }
        else{
            p <- p + geom_point(aes(color=col), size=pt.size, shape=pch.use) +
                scale_color_manual(values=cols.use)
        }
    }
    else{
        if(all(data.plot$gene == data.plot$gene[1])){
            warning(paste0("All cells have the same value of ", feature, "."))
            p <- p + geom_point(color=cols.use[1], size=pt.size, shape=pch.use)
        }
        else{
            p <- p + geom_point(aes(color=gene), size=pt.size, shape=pch.use) +
                scale_color_gradientn(colors=cols.use, guide = guide_colorbar(title = feature))
        }
    }
    if(no.axes){
        p <- p + labs(title = feature, x ="", y="") + theme(
            axis.line = element_blank(),
            axis.text.x = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks = element_blank(),
            axis.title.x = element_blank(),
            axis.title.y = element_blank()
        )
    }
    else{
        p <- p + labs(title = feature, x = dim.codes[1], y = dim.codes[2])
    }
    if(no.legend){
        p <- p + theme(legend.position = 'none')
    }
    if (dark.theme) {
        p <- p + DarkTheme()
    }
    return(p)
}

# Documentation
###############
BlendPlot <- function(
    data.use,
    features.plot,
    data.plot,
    pt.size,
    pch.use,
    cols.use,
    dim.codes,
    min.cutoff,
    max.cutoff,
    no.axes,
    no.legend
) {
    num.cols <- length(x = cols.use)
    #   Create a vector of colors that weren't provided
    cols.not.provided <- colors(distinct = TRUE)
    cols.not.provided <- cols.not.provided[!(grepl(pattern = paste(cols.use, collapse = '|'), x = cols.not.provided, ignore.case = TRUE))]
    if (num.cols > 4) {
        #   If provided more than four colors, take only the first four
        cols.use <- cols.use[c(1:4)]
    } else if ((num.cols == 2) || (num.cols == 3)) {
        #   If two or three colors, use the last two as high values for blending
        #   and add to our vector of colors
        blend <- blendColors(cols.use[c(num.cols - 1, num.cols)])
        cols.use <- c(cols.use, blend)
        if (num.cols == 2) {
            #   If two colors, provided,
            #   we still need a low color
            cols.use <- c(sample(x = cols.not.provided, size = 1), cols.use)
        }
    } else if ((num.cols == 1)) {
        #   If only one color provided
        if (cols.use %in% rownames(x = brewer.pal.info)) {
            #   Was it a palette from RColorBrewer? If so, create
            #   our colors based on the palette
            palette <- brewer.pal(n = 3, name = cols.use)
            cols.use <- c(palette, blendColors(palette[c(2, 3)]))
        } else {
            #   If not, randomly create our colors
            cols.high <- sample(x = cols.not.provided, size = 2, replace = FALSE)
            cols.use <- c(cols.use, cols.high, blendColors(cols.high))
        }
    } else if (num.cols <= 0) {
        cols.use <- c('yellow','red', 'blue', blendColors('red', 'blue'))
    }
    names(x = cols.use) <- c('low', 'high1', 'high2', 'highboth')
    if((length(x = features.plot) != 2) || (length(x = min.cutoff) != 2) || (length(x = max.cutoff) != 2)) {
        stop("An overlayed FeaturePlot only works with two features and requires two minimum and maximum cutoffs")
    }
    #   Check for quantiles
    min.cutoff <- c(
        setQuantile(cutoff = min.cutoff[1], data = data.gene[features.plot[1], ]),
        setQuantile(cutoff = min.cutoff[2], data = data.gene[features.plot[2], ])
    )
    max.cutoff <- c(
        setQuantile(cutoff = max.cutoff[1], data = data.gene[features.plot[1], ]),
        setQuantile(cutoff = max.cutoff[2], data = data.gene[features.plot[2], ])
    )
    data.gene <- na.omit(object = data.frame(data.use[features.plot, ]))
    cell.names <- colnames(x = data.gene)
    #   Minimum and maximum masking
    data.gene <- matrix(
        data = vapply(
            X = data.gene,
            FUN = function(x) ifelse(test = x < min.cutoff, yes = min.cutoff, no = x),
            FUN.VALUE = c(1, 1)
        ),
        nrow = 2
    )
    data.gene <- matrix(
        data = vapply(
            X = as.data.frame(x = data.gene),
            FUN = function(x) ifelse(test = x > max.cutoff, yes = max.cutoff, no = x),
            FUN.VALUE = c(1, 1)
        ),
        nrow = 2
    )
    data.gene <- as.data.frame(x = data.gene)
    rownames(x = data.gene) <- features.plot
    colnames(x = data.gene) <- cell.names
    #   Stuff for break points
    if(all(data.gene ==0)) {
        data.cut <- 0
    } else {
        #   Cut the expression of both features
        cuts <- apply(X = data.gene, MARGIN = 1, FUN = cut, breaks = 2, labels = FALSE)
        cuts.dim <- dim(cuts)
        if (cuts.dim[1] > cuts.dim[2]){
            cuts <- t(x = cuts)
        }
        #   Apply colors dependent on if the cell expresses
        #   none, one, or both features
        data.cut = apply(
            X = cuts,
            MARGIN = 2,
            FUN = function(x) {
                if((x[1] == 1) && (x[2] == 2)) { # Expression in 2
                    'high2'
                } else if((x[1] == 2) && (x[2] == 1)) { # Expression in 1
                    'high1'
                } else if ((x[1] == 2) && (x[2] == 2)) { # Expression in both
                    'highboth'
                } else { # Expression in neither
                    'low'
                }
            }
        )
        data.cut <- as.factor(x = data.cut)
    }
    data.plot$colors <- data.cut
    #   Start plotting
    legend.names <- c('high1' = paste('High', features.plot[1]), 'high2' = paste('High', features.plot[2]), 'highboth' = 'High both')
    title <- paste0(features.plot, collapse = ' x ')
    p <- ggplot(data = data.plot, mapping = aes(x = x, y = y))
    p <- p + geom_point(mapping = aes(color = colors), size = pt.size, shape = pch.use)
    p <- p + scale_color_manual(
        values = cols.use,
        limits = c('high1', 'high2', 'highboth'),
        labels = legend.names,
        guide = guide_legend(title = NULL, override.aes = list(size = 2))
    )
    #   Deal with axes and legends
    if(no.axes) {
        p <- p + labs(title = title, x ="", y="") + theme(
            axis.line = element_blank(),
            axis.text.x = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks = element_blank(),
            axis.title.x = element_blank(),
            axis.title.y = element_blank()
        )
    } else {
        p <- p + labs(title = title, x = dim.codes[1], y = dim.codes[2])
    }
    if(no.legend){
        p <- p + theme(legend.position = 'none')
    }
    if (dark.theme) {
        p <- p + DarkTheme()
    }
    return(p)
}

# Documentation
###############
blendColors <- function(..., as.rgb = FALSE) {
    #   Assemble the arguments passed into a character vector
    colors <- as.character(x = c(...))
    if (length(x = colors) < 2) {
        stop("Please provide two or more colors to blend")
    }
    #   Check for hexadecimal values for automatic alpha blending
    alpha.value <- 255
    if (sum(sapply(X = colors, FUN = grepl, pattern = '^#')) != 0) {
        hex <- colors[which(x = grepl(pattern = '^#', x = colors))]
        hex.length <- sapply(X = hex, FUN = nchar)
        #   9-character hexadecimal values specify alpha levels
        if (9 %in% hex.length) {
            hex.alpha <- hex[which(x = hex.length == 9)]
            hex.vals <- sapply(X = hex.alpha, FUN = substr, start = 8, stop = 9)
            dec.vals <- sapply(X = hex.vals, FUN = strtoi, base = 16)
            dec.vals <- dec.vals / 255 # Convert to 0:1 scale for calculations
            alpha.value <- dec.vals[1]
            #   Blend alpha levels, going top-down
            for (val in dec.vals[-1]) {
                alpha.value <- alpha.value + (val * (1 - alpha.value))
            }
            alpha.value <- alpha.value * 255 # Convert back to 0:255 scale
        }
    }
    #   Convert to a 3 by `length(colors)` matrix of RGB values
    rgb.vals <- sapply(X = colors, FUN = col2rgb)
    if (nrow(x = rgb.vals) != 3) {
        rgb.vals <- t(x = rgb.vals)
    }
    #   Blend together using the additive method
    #   Basically, resulting colors are the mean of the component colors
    blend <- apply(X = rgb.vals, MARGIN = 1, FUN = function(component) mean(x = component))
    #   If we're returning RGB values, convert to matrix, just like col2rgb
    #   Otherwise, return as hexadecimal; can be used directly for plotting
    if (as.rgb) {
        result <- matrix(
            data = blend,
            nrow = 3,
            dimnames = list(c('red', 'green', 'blue'), 'blend')
        )
    } else {
        result <- rgb(matrix(data = blend, ncol = 3), alpha = alpha.value, maxColorValue = 255)
    }
    return(result)
}

# Documentation
###############
setQuantile <- function(cutoff, data) {
    regex.quantile <- '^q[0-9]{1,2}$'
    if (grepl(pattern = regex.quantile, x = as.character(x = cutoff), perl = TRUE)) {
        this.quantile <- as.numeric(x = sub(pattern = 'q', replacement = '', x = as.character(x = cutoff))) / 100
        cutoff <- quantile(x = unlist(x = data), probs = this.quantile)
    }
    return(as.numeric(x = cutoff))
}

#' Vizualization of multiple features
#'
#' Similar to FeaturePlot, however, also splits the plot by visualizing each
#' identity class separately.
#'
#' Particularly useful for seeing if the same groups of cells co-exhibit a
#' common feature (i.e. co-express a gene), even within an identity class. Best
#' understood by example.
#'
#'
#' @param object Seurat object
#' @param features.plot Vector of features to plot
#' @param dim.1 Dimension for x-axis (default 1)
#' @param dim.2 Dimension for y-axis (default 2)
#' @param idents.use Which identity classes to display (default is all identity
#' classes)
#' @param pt.size Adjust point size for plotting
#' @param cols.use Ordered vector of colors to use for plotting. Default is
#' heat.colors(10).
#' @param pch.use Pch for plotting
#' @param reduction.use Which dimensionality reduction to use. Default is
#' "tsne", can also be "pca", or "ica", assuming these are precomputed.
#' @param group.by Group cells in different ways (for example, orig.ident)
#' @param sep.scale Scale each group separately. Default is FALSE.
#' @param max.exp Max cutoff for scaled expression value
#' @param min.exp Min cutoff for scaled expression value
#' @param rotate.key rotate the legend
#' @param plot.horiz rotate the plot such that the features are columns, groups are the rows
#' @param key.position position of the legend ("top", "right", "bottom", "left")
#' @param do.return Return the ggplot2 object
#'
#' @return No return value, only a graphical output
#'
#' @importFrom dplyr %>% mutate_each group_by select ungroup
#'
#' @export
#'
FeatureHeatmap <- function(
  object, features.plot, dim.1 = 1, dim.2 = 2, idents.use = NULL, pt.size = 2,
  cols.use = c("grey", "red"), pch.use = 16, reduction.use = "tsne",
  group.by = NULL, sep.scale = FALSE, do.return = FALSE, min.exp = -Inf,
  max.exp = Inf, rotate.key = F, plot.horiz = FALSE, key.position = "right"
) {
  if (! is.null(x = group.by)) {
    object <- SetAllIdent(object = object, id = group.by)
  }
  idents.use <- set.ifnull(x = idents.use, y = sort(x = unique(x = object@ident)))
  par(mfrow = c(length(x = features.plot), length(x = idents.use)))
  dim.code <- translate.dim.code(object = object, reduction.use = reduction.use)
  dim.codes <- paste0(dim.code, c(dim.1, dim.2))
  data.plot <- data.frame(FetchData(
    object = object,
    vars.all = c(dim.codes, features.plot)
  ))
  colnames(x = data.plot)[1:2] <- c("dim1", "dim2")
  data.plot$ident <- as.character(x = object@ident)
  data.plot$cell <- rownames(x = data.plot)
  features.plot <- gsub('-', '\\.', features.plot)
  data.plot  %>% gather(gene, expression, features.plot, -dim1, -dim2, -ident, -cell) -> data.plot
  if (sep.scale) {
    data.plot %>% group_by(ident, gene) %>% mutate(scaled.expression = scale(expression)) -> data.plot
  } else {
    data.plot %>%  group_by(gene) %>% mutate(scaled.expression = scale(expression)) -> data.plot
  }
  data.plot$gene <- factor(x = data.plot$gene, levels = features.plot)
  data.plot$scaled.expression <- minmax(
    data = data.plot$scaled.expression,
    min = min.exp,
    max = max.exp
  )
  if (rotate.key) {
    key.direction <- "horizontal"
    key.title.pos <- "top"
  } else {
    key.direction <- "vertical"
    key.title.pos <- "left"
  }
  p <- ggplot(data = data.plot, mapping = aes(x = dim1, y = dim2)) +
    geom_point(mapping = aes(colour = scaled.expression), size = pt.size)
  if (rotate.key) {
    p <- p + scale_colour_gradient(
      low = cols.use[1],
      high = cols.use[2],
      guide = guide_colorbar(
        direction = key.direction,
        title.position = key.title.pos,
        title = "Scaled Expression"
      )
    )
  } else {
    p <- p + scale_colour_gradient(
      low = cols.use[1],
      high = cols.use[2],
      guide = guide_colorbar(title = "Scaled Expression")
    )
  }
  if(plot.horiz){
    p <- p + facet_grid(ident ~ gene)
  }
  else{
    p <- p + facet_grid(gene ~ ident)
  }
  p2 <- p +
    theme_bw() +
    nogrid +
    ylab(label = dim.codes[2]) +
    xlab(label = dim.codes[1])
  p2 <- p2 + theme(legend.position = key.position)
  if (do.return) {
    return(p2)
  }
  print(p2)
}

# Documentation
###############
translate.dim.code <- function(object, reduction.use) {
  if (is.null(x = reduction.use)) {
    return.code <- "PC"
  }
  else {
    return.code <- object@dr[[reduction.use]]@key
  }
  return(return.code)
}

# Documentation
###############
#Cool, but not supported right now
SpatialDe <- function(object, marker.cells, genes.use = NULL) {
  object <- p15
  embed.map <- object@tsne.rot
  mult.use <- 2
  mult.use.far <- 10
  if ((mult.use.far * length(x = marker.cells)) > nrow(x = embed.map)) {
    mult.use.far <- 1
    mult.use <- 1
  }
  genes.use <- set.ifnull(x = genes.use, y = object@var.genes)
  marker.pos <- apply(X = embed.map[marker.cells, ], MARGIN = 2, FUN = mean)
  embed.map <- rbind(embed.map, marker.pos)
  rownames(x = embed.map)[nrow(x = embed.map)] <- "marker"
  embed.dist <- sort(x = as.matrix(x = dist(x = (embed.map)))["marker", ])
  embed.diff <- names(x = embed.dist[! (names(x = embed.dist) %in% marker.cells)][1:(mult.use * length(x = marker.cells))][-1])
  embed.diff.far <- names(x = embed.dist[! (names(x = embed.dist) %in% marker.cells)][1:(mult.use.far * length(x = marker.cells))][-1])
  diff.genes <- rownames(
    x = subset(
      x = DiffExpTest(
        object = p15,
        cells.1 = marker.cells,
        cells.2 = embed.diff,
        genes.use = genes.use
      ),
      subset = p_val < (1e-5)
    )
  )
  diff.genes <- subset(
    x = DiffExpTest(
      object = p15,
      cells.1 = marker.cells,
      cells.2 = embed.diff,
      genes.use = diff.genes
    ),
    subset = p_val<(1e-10)
  )
  return(diff.genes)
}

#' Perform spectral density clustering on single cells
#'
#' Find point clounds single cells in a two-dimensional space using density clustering (DBSCAN).
#'
#' @param object Seurat object
#' @param dim.1 First dimension to use
#' @param dim.2 second dimension to use
#' @param reduction.use Which dimensional reduction to use (either 'pca' or 'ica')
#' @param G.use Parameter for the density clustering. Lower value to get more fine-scale clustering
#' @param set.ident TRUE by default. Set identity class to the results of the density clustering.
#' Unassigned cells (cells that cannot be assigned a cluster) are placed in cluster 1, if there are any.
#' @param seed.use Random seed for the dbscan function
#' @param ... Additional arguments to be passed to the dbscan function
#'
#' @export
#'
DBClustDimension <- function(
  object,
  dim.1 = 1,
  dim.2 = 2,
  reduction.use = "tsne",
  G.use = NULL,
  set.ident = TRUE,
  seed.use = 1,
  ...
) {
  dim.code <- translate.dim.code(object = object, reduction.use = reduction.use)
  dim.codes <- paste0(dim.code, c(dim.1, dim.2))
  data.plot <- FetchData(object = object, vars.all = dim.codes)
  x1 <- paste0(dim.code, dim.1)
  x2 <- paste0(dim.code, dim.2)
  data.plot$x <- data.plot[, x1]
  data.plot$y <- data.plot[, x2]
  set.seed(seed = seed.use)
  data.mclust <- ds <- dbscan(data = data.plot[, c("x", "y")], eps = G.use, ...)
  to.set <- as.numeric(x = data.mclust$cluster + 1)
  data.names <- names(x = object@ident)
  object@data.info[data.names, "DBclust.ident"] <- to.set
  if (set.ident) {
    object@ident <- factor(x = to.set)
    names(x = object@ident) <- data.names
  }
  return(object)
}

#' Perform spectral k-means clustering on single cells
#'
#' Find point clounds single cells in a low-dimensional space using k-means clustering.
#'
#' CAn be useful for smaller datasets, not documented here yet
#' @export
#'
KClustDimension <- function(
  object,
  dim.1 = 1,
  dim.2 = 2,
  cells.use = NULL,
  pt.size = 4,
  reduction.use = "tsne",
  k.use = 5,
  set.ident = FALSE,
  seed.use = 1,
  ...
) {
  cells.use <- set.ifnull(x = cells.use, y = colnames(x = object@data))
  dim.code <- "PC"
  if (reduction.use == "pca") {
    data.plot <- object@pca.rot[cells.use, ]
  }
  if (reduction.use == "tsne") {
    data.plot <- object@tsne.rot[cells.use, ]
    dim.code <- "Seurat_"
  }
  if (reduction.use == "ica") {
    data.plot <- object@ica.rot[cells.use, ]
    dim.code <- "IC"
  }
  x1 <- paste0(dim.code, dim.1)
  x2 <- paste0(dim.code, dim.2)
  data.plot$x <- data.plot[, x1]
  data.plot$y <- data.plot[, x2]
  if (reduction.use != "pca") {
    set.seed(seed = seed.use)
    data.mclust <- ds <- kmeans(x = data.plot[, c("x", "y")], centers = k.use)
  }
  if (reduction.use == "pca") {
    set.seed(seed = seed.use)
    data.mclust <- ds <- kmeans(x = object@pca.rot[cells.use, dim.1], centers = k.use)
  }
  to.set <- as.numeric(x = data.mclust$cluster)
  data.names <- names(x = object@ident)
  object@data.info[data.names, "kdimension.ident"] <- to.set
  if (set.ident) {
    object@ident <- factor(x = to.set)
    names(x = object@ident) <- data.names
  }
  return(object)
}

#' Significant genes from a PCA
#'
#' Returns a set of genes, based on the JackStraw analysis, that have
#' statistically significant associations with a set of PCs.
#'
#' @param object Seurat object
#' @param pcs.use PCS to use.
#' @param pval.cut P-value cutoff
#' @param use.full Use the full list of genes (from the projected PCA). Assumes
#' that ProjectPCA has been run. Default is TRUE
#' @param max.per.pc Maximum number of genes to return per PC. Used to avoid genes from one PC dominating the entire analysis.
#'
#' @return A vector of genes whose p-values are statistically significant for
#' at least one of the given PCs.
#'
#' @export
#'
PCASigGenes <- function(
  object,
  pcs.use,
  pval.cut = 0.1,
  use.full = TRUE,
  max.per.pc = NULL
) {
  pvals.use <- object@jackStraw.empP
  pcx.use <- object@pca.x
  if (use.full) {
    pvals.use <- object@jackStraw.empP.full
    pcx.use <- object@pca.x.full
  }
  if (length(x = pcs.use) == 1) {
    pvals.min <- pvals.use[, pcs.use]
  }
  if (length(x = pcs.use) > 1) {
    pvals.min <- apply(X = pvals.use[, pcs.use], MARGIN = 1, FUN = min)
  }
  names(x = pvals.min) <- rownames(x = pvals.use)
  genes.use <- names(x = pvals.min)[pvals.min < pval.cut]
  if (! is.null(x = max.per.pc)) {
    pc.top.genes <- PCTopGenes(
      object = object,
      pc.use = pcs.use,
      num.genes = max.per.pc,
      use.full = use.full,
      do.balanced = FALSE
    )
    genes.use <- ainb(a = pc.top.genes, b = genes.use)
  }
  return(genes.use)
}

#' Gene expression heatmap
#'
#' Draws a heatmap of single cell gene expression using the heatmap.2 function.
#'
#' @param object Seurat object
#' @param cells.use Cells to include in the heatmap (default is all cells)
#' @param genes.use Genes to include in the heatmap (ordered)
#' @param disp.min Minimum display value (all values below are clipped)
#' @param disp.max Maximum display value (all values above are clipped)
#' @param draw.line Draw vertical lines delineating cells in different identity
#' classes.
#' @param do.return Default is FALSE. If TRUE, return a matrix of scaled values
#' which would be passed to heatmap.2
#' @param order.by.ident Order cells in the heatmap by identity class (default
#' is TRUE). If FALSE, cells are ordered based on their order in cells.use
#' @param col.use Color palette to use
#' @param slim.col.label if (order.by.ident==TRUE) then instead of displaying
#' every cell name on the heatmap, display only the identity class name once
#' for each group
#' @param group.by If (order.by.ident==TRUE) default,  you can group cells in
#' different ways (for example, orig.ident)
#' @param remove.key Removes the color key from the plot.
#' @param cex.col positive numbers, used as cex.axis in for the column axis labeling.
#' The defaults currently only use number of columns
#' @param do.scale whether to use the data or scaled data
#' @param ... Additional parameters to heatmap.2. Common examples are cexRow
#' and cexCol, which set row and column text sizes
#'
#' @return If do.return==TRUE, a matrix of scaled values which would be passed
#' to heatmap.2. Otherwise, no return value, only a graphical output
#'
#' @importFrom gplots heatmap.2
#'
#' @export
#'
DoHeatmap <- function(
  object,
  cells.use = NULL,
  genes.use = NULL,
  disp.min = NULL,
  disp.max = NULL,
  draw.line = TRUE,
  do.return = FALSE,
  order.by.ident = TRUE,
  col.use = pyCols,
  slim.col.label = FALSE,
  group.by = NULL,
  remove.key = FALSE,
  cex.col = NULL,
  do.scale = TRUE,
  ...
) {
  cells.use <- set.ifnull(x = cells.use, y = object@cell.names)
  cells.use <- ainb(a = cells.use, b = object@cell.names)
  cells.ident <- object@ident[cells.use]
  if (! is.null(x = group.by)) {
    cells.ident <- factor(x = FetchData(object = object, vars.all = group.by)[, 1])
  }
  cells.ident <- factor(
    x = cells.ident,
    labels = ainb(a = levels(x = cells.ident), b = cells.ident)
  )
  if (order.by.ident) {
    cells.use <- cells.use[order(cells.ident)]
  } else {
    cells.ident <- factor(
      x = cells.ident,
      levels = as.vector(x = unique(x = cells.ident))
    )
  }
  #determine assay type
  data.use <- NULL
  assays.use <- c("RNA", names(x = object@assay))
  slot.use <- "scale.data"
  if (do.scale) {
    if ((is.null(x = disp.min) || is.null(x = disp.max))) {
      disp.min <- -2.5
      disp.max <- 2.5
    }
  } else {
    slot.use <- "data"
    if ((is.null(x = disp.min) || is.null(x = disp.max))) {
      disp.min <- -Inf
      disp.max <- Inf
    }
  }
  for (assay.check in assays.use) {
    data.assay <- GetAssayData(object = object, assay.type = assay.check, slot = slot.use)
    genes.intersect <- intersect(x = genes.use, y = rownames(x = data.assay))
    new.data <- data.assay[genes.intersect, cells.use, drop = FALSE]
    if (! (is.matrix(x = new.data))) {
      new.data <- as.matrix(x = new.data)
    }
    data.use <- rbind(data.use, new.data)
  }
  data.use <- minmax(data = data.use, min = disp.min, max = disp.max)
  vline.use <- NULL
  colsep.use <- NULL
  if (remove.key) {
    hmFunction <- heatmap2NoKey
  } else {
    hmFunction <- heatmap.2
  }
  if (draw.line) {
    colsep.use <- cumsum(x = table(cells.ident))
  }
  if (slim.col.label && order.by.ident) {
    col.lab <- rep("", length(x = cells.use))
    col.lab[round(x = cumsum(x = table(cells.ident)) - table(cells.ident) / 2) + 1] <- levels(x = cells.ident)
    cex.col <- set.ifnull(x = cex.col,y = 0.2 + 1 / log10(x = length(x = unique(x = cells.ident))))
    hmFunction(
      data.use,
      Rowv = NA,
      Colv = NA,
      trace = "none",
      col = col.use,
      colsep = colsep.use,
      labCol = col.lab,
      cexCol = cex.col,
      ...
    )
  } else if (slim.col.label) {
    col.lab = rep("", length(x = cells.use))
    cex.col <- set.ifnull(x = cex.col, y = 0.2 + 1 / log10(x = length(x = unique(x = cells.ident))))
    hmFunction(
      data.use,
      Rowv = NA,
      Colv = NA,
      trace = "none",
      col = col.use,
      colsep = colsep.use,
      labCol = col.lab,
      cexCol = cex.col,
      ...
    )
  } else {
    hmFunction(
      data.use,
      Rowv = NA,
      Colv = NA,
      trace = "none",
      col = col.use,
      colsep = colsep.use,
      ...
    )
  }
  if (do.return) {
    return(data.use)
  }
}

#' Node Heatmap
#'
#' Takes an object, a marker list (output of FindAllMarkers), and a node
#' and plots a heatmap where genes are ordered vertically by the splits present
#' in the object@@cluster.tree slot.
#'
#' @param object Seurat object. Must have the cluster.tree slot filled (use BuildClusterTree)
#' @param marker.list List of marker genes given from the FindAllMarkersNode function
#' @param node Node in the cluster tree from which to start the plot, defaults to highest node in marker list
#' @param max.genes Maximum number of genes to keep for each division
#' @param ... Additional parameters to pass to DoHeatmap
#'
#' @importFrom dplyr %>% group_by filter top_n select
#'
#' @return Plots heatmap. No return value.
#'
#' @export
#'
HeatmapNode <- function(object, marker.list, node = NULL, max.genes = 10, ...) {
  tree <- object@cluster.tree[[1]]
  node <- set.ifnull(node, min(marker.list$cluster))
  node.order <- c(node, DFT(tree = tree, node = node))
  marker.list$rank <- seq(1:nrow(x = marker.list))
  marker.list %>% group_by(cluster) %>% filter(avg_diff > 0) %>%
    top_n(max.genes, -rank) %>% select(gene, cluster) -> pos.genes
  marker.list %>% group_by(cluster) %>% filter(avg_diff < 0) %>%
    top_n(max.genes, -rank) %>% select(gene, cluster) -> neg.genes
  gene.list <- vector()
  node.stack <- vector()
  for (n in node.order) {
    if (NodeHasChild(tree = tree, node = n)) {
      gene.list <- c(
        gene.list,
        c(
          subset(x = pos.genes, subset = cluster == n)$gene,
          subset(x = neg.genes, subset = cluster == n)$gene
        )
      )
      if (NodeHasOnlyChildren(tree = tree, node = n)) {
        gene.list <- c(
          gene.list,
          subset(x = neg.genes, subset = cluster == node.stack[length(node.stack)])$gene
        )
        node.stack <- node.stack[-length(x = node.stack)]
      }
    }
    else {
      gene.list <- c(gene.list, subset(x = pos.genes, subset = cluster == n)$gene)
      node.stack <- append(x = node.stack, values = n)
    }
  }
  #gene.list <- rev(unique(rev(gene.list)))
  descendants <- getDescendants(tree = tree, node = node)
  children <- descendants[!descendants %in% tree$edge[, 1]]
  all.children <- tree$edge[,2][!tree$edge[,2] %in% tree$edge[, 1]]
  DoHeatmap(
    object = object,
    cells.use = WhichCells(object = object, ident = children),
    genes.use = gene.list,
    slim.col.label = TRUE,
    remove.key = TRUE,
    ...
  )
}

#' K-Means Clustering
#'
#' Perform k=means clustering on both genes and single cells
#'
#' K-means and heatmap are calculated on object@@scale.data
#'
#' @param object Seurat object
#' @param genes.use Genes to use for clustering
#' @param k.genes K value to use for clustering genes
#' @param k.cells K value to use for clustering cells (default is NULL, cells
#' are not clustered)
#' @param k.seed Random seed
#' @param do.plot Draw heatmap of clustered genes/cells (default is TRUE)
#' @param data.cut Clip all z-scores to have an absolute value below this.
#' Reduces the effect of huge outliers in the data.
#' @param k.cols Color palette for heatmap
#' @param pc.row.order Order gene clusters based on the average PC score within
#' a cluster. Can be useful if you want to visualize clusters, for example,
#' based on their average score for PC1.
#' @param pc.col.order Order cell clusters based on the average PC score within
#' a cluster
#' @param rev.pc.order Use the reverse PC ordering for gene and cell clusters
#' (since the sign of a PC is arbitrary)
#' @param use.imputed Cluster imputed values (default is FALSE)
#' @param set.ident If clustering cells (so k.cells>0), set the cell identity
#' class to its K-means cluster (default is TRUE)
#' @param do.constrained FALSE by default. If TRUE, use the constrained K-means function implemented in the tclust package.
#' @param \dots Additional parameters passed to kmeans (or tkmeans)
#'
#' @importFrom tclust tkmeans
#'
#' @return Seurat object where the k-means results for genes is stored in
#'
#' object@@kmeans.obj[[1]], and the k-means results for cells is stored in
#' object@@kmeans.col[[1]]. The cluster for each cell is stored in object@@data.info[,"kmeans.ident"]
#' and also object@@ident (if set.ident=TRUE)
#'
#' @export
DoKMeans <- function(
  object,
  genes.use = NULL,
  k.genes = NULL,
  k.cells = 0,
  k.seed = 1,
  do.plot = TRUE,
  data.cut = 2.5,
  k.cols = pyCols,
  pc.row.order = NULL,
  pc.col.order = NULL,
  rev.pc.order = FALSE,
  use.imputed = FALSE,
  set.ident = TRUE,
  do.constrained = FALSE,
  ...
) {
  if (use.imputed) {
    data.use.orig <- data.frame(t(x = scale(x = t(x = object@imputed))))
  } else {
    data.use.orig <- object@scale.data
  }
  data.use <- minmax(data = data.use.orig, min = data.cut * (-1), max = data.cut)
  if (rev.pc.order) {
    revFxn <- function(x) {
      max(x) + 1 - x
    }
  } else {
    revFxn <- same
  }
  kmeans.col <- NULL
  genes.use <- set.ifnull(x = genes.use, y = object@var.genes)
  genes.use <- genes.use[genes.use %in% rownames(x = data.use)]
  cells.use <- object@cell.names
  kmeans.data <- data.use[genes.use, cells.use]
  if (do.constrained) {
    set.seed(seed = k.seed)
    kmeans.obj <- tkmeans(objectkmeans.data, k = k.genes, ...)
  } else {
    set.seed(seed = k.seed)
    kmeans.obj <- kmeans(x = kmeans.data, centers = k.genes, ...)
  }
  if (! is.null(x = pc.row.order)) {
    if (nrow(x = object@pca.x.full > 0)) {
      pcx.use <- object@pca.x.full
      pc.genes <- ainb(a = genes.use, b = rownames(x = pcx.use))
    } else {
      pcx.use <- object@pca.x
      pc.genes <- ainb(a = genes.use, b = rownames(x = pcx.use))
    }
    kmeans.obj$cluster <- as.numeric(
      x = revFxn(
        rank(
          x = tapply(
            X = pcx.use[genes.use, pc.row.order],
            INDEX = as.numeric(x = kmeans.obj$cluster),
            FUN = mean
          )
        )
      )[as.numeric(x = kmeans.obj$cluster)]
    )
  }
  names(x = kmeans.obj$cluster) <- genes.use
  if (k.cells > 0) {
    kmeans.col <- kmeans(x = t(x = kmeans.data), centers = k.cells)
    if (! is.null(x = pc.col.order)) {
      kmeans.col$cluster <- as.numeric(
        x = revFxn(
          rank(
            x = tapply(
              X = object@pca.rot[cells.use, pc.col.order],
              INDEX = kmeans.col$cluster,
              FUN = mean
            )
          )
        )[as.numeric(x = kmeans.col$cluster)]
      )
    }
    names(x = kmeans.col$cluster) <- cells.use
  }
  object@kmeans.obj <- list(kmeans.obj)
  object@kmeans.col <- list(kmeans.col)
  kmeans.obj <- object@kmeans.obj[[1]]
  kmeans.col <- object@kmeans.col[[1]]
  if (k.cells > 0) {
    object@data.info[names(x = kmeans.col$cluster), "kmeans.ident"] <- kmeans.col$cluster
  }
  if (set.ident && (k.cells > 0)) {
    object <- SetIdent(
      object = object,
      cells.use = names(x = kmeans.col$cluster),
      ident.use = kmeans.col$cluster
    )
  }
  if (do.plot) {
    disp.data <- minmax(
      data = kmeans.data[order(kmeans.obj$cluster[genes.use]), ],
      min = data.cut * (-1),
      max = data.cut
    )
    #DoHeatmap(object,object@cell.names,names(sort(kmeans.obj$cluster)),data.cut*(-1),data.cut,col.use = k.cols,...)
    KMeansHeatmap(object = object)
  }
  return(object)
}

# Documentation
###############
#' @export
#'
GenesInCluster <- function(object, cluster.num, max.genes = 1e6) {
  toReturn <- unlist(
    x = lapply(
      X = cluster.num,
      FUN = function(x) {
        return(head(
          x = sort(x = names(x = which(x = object@kmeans.obj[[1]]$cluster==x))),
          n = max.genes
        ))
      }
    )
  )
  return(toReturn)
}

# Documentation
###############
#Needs comments
#' @export
#'
KMeansHeatmap <- function(
  object,
  cells.use = object@cell.names,
  genes.cluster = NULL,
  max.genes = 1e6,
  slim.col.label = TRUE,
  remove.key = TRUE,
  row.lines = TRUE,
  ...
) {
  genes.cluster <- set.ifnull(
    x = genes.cluster,
    y = unique(x = object@kmeans.obj[[1]]$cluster)
  )
  genes.use <- GenesInCluster(
    object = object,
    cluster.num = genes.cluster,
    max.genes = max.genes
  )
  cluster.lengths <- sapply(
    X = genes.cluster,
    FUN = function(x) {
      return(length(x = GenesInCluster(object = object, cluster.num = x)))
    }
  )
  print(cluster.lengths)
  if (row.lines) {
    rowsep.use <- cumsum(x = cluster.lengths)
  } else {
    rowsep.use <- NA
  }
  DoHeatmap(
    object = object,
    cells.use = cells.use,
    genes.use = genes.use,
    slim.col.label = slim.col.label,
    remove.key = remove.key,
    rowsep = rowsep.use,
    ...
  )
}

# Documentation
###############
#' @export
#'
CellCorMatrix <- function(
  object,
  cor.genes = NULL,
  cell.inds = NULL,
  do.k = FALSE,
  k.seed = 1,
  k.num = 4,
  vis.low = (-1),
  vis.high = 1,
  vis.one = 0.8,
  pcs.use = 1:3,
  col.use = pyCols
) {
  cor.genes <- set.ifnull(x = cor.genes, y = object@var.genes)
  cell.inds <- set.ifnull(x = cell.inds, y = colnames(x = object@data))
  cor.genes <- cor.genes[cor.genes %in% rownames(x = object@data)]
  data.cor <- object@scale.data[cor.genes, cell.inds]
  cor.matrix <- cor(x = data.cor)
  set.seed(seed = k.seed)
  kmeans.cor <- kmeans(x = cor.matrix, centers = k.num)
  if (do.k) {
    cor.matrix <- cor.matrix[order(kmeans.cor$cluster), order(kmeans.cor$cluster)]
  }
  kmeans.names <- rownames(x = cor.matrix)
  row.annot <- data.frame(
    cbind(
      kmeans.cor$cluster[kmeans.names],
      object@pca.rot[kmeans.names, pcs.use]
    )
  )
  colnames(x = row.annot) <- c("K", paste0("PC", pcs.use))
  cor.matrix[cor.matrix == 1] <- vis.one
  cor.matrix <- minmax(data = cor.matrix, min = vis.low, max = vis.high)
  object@kmeans.cell <- list(kmeans.cor)
  if (do.k) {
    aheatmap(
      x = cor.matrix,
      col = col.use,
      Rowv = NA,
      Colv = NA,
      annRow = row.annot
    )
  } else {
    heatmap.2(
      x = cor.matrix,
      trace = "none",
      Rowv = NA,
      Colv = NA,
      col = pyCols
    )
  }
  return(object)
}

# Documentation
###############
#' @export
#'
GeneCorMatrix <- function(
  object,
  cor.genes = NULL,
  cell.inds = NULL,
  do.k = FALSE,
  k.seed = 1,
  k.num = 4,
  vis.low = (-1),
  vis.high = 1,
  vis.one = 0.8,
  pcs.use = 1:3,
  col.use = pyCols
) {
  cor.genes <- set.ifnull(x = cor.genes, y = object@var.genes)
  cell.inds <- set.ifnull(x = cell.inds, y = colnames(x = object@data))
  cor.genes <- cor.genes[cor.genes %in% rownames(x = object@data)]
  data.cor <- object@data[cor.genes, cell.inds]
  cor.matrix <- cor(x = t(x = data.cor))
  set.seed(seed = k.seed)
  kmeans.cor <- kmeans(x = cor.matrix, centers = k.num)
  cor.matrix <- cor.matrix[order(kmeans.cor$cluster), order(kmeans.cor$cluster)]
  kmeans.names <- rownames(x = cor.matrix)
  row.annot <- data.frame(
    cbind(
      kmeans.cor$cluster[kmeans.names],
      object@pca.x[kmeans.names, pcs.use]
    )
  )
  colnames(x = row.annot) <- c("K", paste0("PC", pcs.use))
  cor.matrix[cor.matrix == 1] <- vis.one
  cor.matrix <- minmax(data = cor.matrix, min = vis.low, max = vis.high)
  object@kmeans.gene <- list(kmeans.cor)
  if (do.k) {
    aheatmap(
      x = cor.matrix,
      col = col.use,
      Rowv = NA,
      Colv = NA,
      annRow = row.annot
    )
  } else {
    aheatmap(
      x = cor.matrix,
      col = col.use,
      annRow = row.annot
    )
  }
  return(object)
}

# Documentation
###############
#' @export
#'
CalinskiPlot <- function(
  object,
  pcs.use,
  pval.cut = 0.1,
  gene.max = 15,
  col.max = 25,
  use.full = TRUE
) {
  if (length(x = pcs.use) == 1) {
    pvals.min <- object@jackStraw.empP.full[, pcs.use]
  } else if (length(x = pcs.use) > 1) {
    pvals.min <- apply(
      X = object@jackStraw.empP.full[, pcs.use],
      MARGIN = 1,
      FUN = min
    )
  }
  names(x = pvals.min) <- rownames(x = object@jackStraw.empP.full)
  genes.use <- names(x = pvals.min)[pvals.min < pval.cut]
  genes.use <- genes.use[genes.use %in% rownames(do.NULL = object@scale.data)]
  par(mfrow = c(1, 2))
  mydata <- object@scale.data[genes.use, ]
  wss <- (nrow(x = mydata) - 1) * sum(apply(X = mydata, MARGIN = 2, FUN = var))
  for (i in 1:gene.max) {
    wss[i] <- sum(kmeans(x = mydata, centers=i)$withinss)
  }
  plot(
    x = 1:gene.max,
    y = wss,
    type = "b",
    xlab= "Number of Clusters for Genes",
    ylab = "Within groups sum of squares"
  )
  mydata <- t(x = object@scale.data[genes.use, ])
  wss <- (nrow(x = mydata) - 1) * sum(apply(X = mydata, MARGIN = 2, FUN = var))
  for (i in 1:col.max) {
    wss[i] <- sum(kmeans(x = mydata, centers = i)$withinss)
  }
  plot(
    x = 1:col.max,
    y = wss,
    type = "b",
    xlab = "Number of Clusters for Cells",
    ylab = "Within groups sum of squares"
  )
  rp()
  return(object)
}

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

#' Add Metadata
#'
#' Adds additional data for single cells to the Seurat object. Can be any piece
#' of information associated with a cell (examples include read depth,
#' alignment rate, experimental batch, or subpopulation identity). The
#' advantage of adding it to the Seurat object is so that it can be
#' analyzed/visualized using FetchData, VlnPlot, GenePlot, SubsetData, etc.
#'
#' @param object Seurat object
#' @param metadata Data frame where the row names are cell names (note : these
#' must correspond exactly to the items in object@@cell.names), and the columns
#' are additional metadata items.
#' @param col.name Name for metadata if passing in single vector of information
#'
#' @return Seurat object where the additional metadata has been added as
#' columns in object@@data.info
#'
#' @export
#'
AddMetaData <- function(object, metadata, col.name = NULL) {
  if (typeof(x = metadata) != "list") {
    metadata <- as.data.frame(x = metadata)
    if (is.null(x = col.name)) {
      stop("Please provide a name for provided metadata")
    }
    colnames(x = metadata) <- col.name
  }
  cols.add <- colnames(x = metadata)
  object@data.info[, cols.add] <- metadata[rownames(x = object@data.info), ]
  return(object)
}

#' JackStraw Plot
#'
#' Plots the results of the JackStraw analysis for PCA significance. For each
#' PC, plots a QQ-plot comparing the distribution of p-values for all genes
#' across each PC, compared with a uniform distribution. Also determines a
#' p-value for the overall significance of each PC (see Details).
#'
#' Significant PCs should show a p-value distribution (black curve) that is
#' strongly skewed to the left compared to the null distribution (dashed line)
#' The p-value for each PC is based on a proportion test comparing the number
#' of genes with a p-value below a particular threshold (score.thresh), compared with the
#' proportion of genes expected under a uniform distribution of p-values.
#'
#' @param object Seurat plot
#' @param plot.x.lim X-axis maximum on each QQ plot.
#' @param plot.y.lim Y-axis maximum on each QQ plot.
#' @param PCs Which PCs to examine
#' @param nCol Number of columns
#' @param score.thresh Threshold to use for the proportion test of PC
#' significance (see Details)
#'
#' @return No value returned, just the PC plots.
#'
#' @author Thanks to Omri Wurtzel for integrating with ggplot
#'
#' @import gridExtra
#'
#' @export
#'
JackStrawPlot <- function(
  object,
  PCs = 1:5,
  nCol = 3,
  score.thresh = 1e-5,
  plot.x.lim = 0.1,
  plot.y.lim = 0.3
) {
  pAll <- object@jackStraw.empP
  pAll <- pAll[, PCs, drop = FALSE]
  pAll$Contig <- rownames(x = pAll)
  pAll.l <- melt(data = pAll, id.vars = "Contig")
  colnames(x = pAll.l) <- c("Contig", "PC", "Value")
  qq.df <- NULL
  score.df <- NULL
  for (i in PCs) {
    q <- qqplot(x = pAll[, i], y = runif(n = 1000), plot.it = FALSE)
    #pc.score=mean(q$y[which(q$x <=score.thresh)])
    pc.score <- prop.test(
      x = c(
        length(x = which(x = pAll[, i] <= score.thresh)),
        floor(x = nrow(x = pAll) * score.thresh)
      ),
      n = c(nrow(pAll), nrow(pAll))
    )$p.val
    if (length(x = which(x = pAll[, i] <= score.thresh)) == 0) {
      pc.score=1
    }
    if (is.null(x = score.df)) {
      score.df <- data.frame(PC=paste("PC",i, sep=""), Score=pc.score)
    } else {
      score.df <- rbind(score.df, data.frame(PC=paste("PC",i, sep=""), Score=pc.score))
    }
    if (is.null(x = qq.df)) {
      qq.df <- data.frame(x = q$x, y = q$y, PC = paste0("PC", i))
    } else {
      qq.df <- rbind(qq.df, data.frame(x = q$x, y = q$y, PC = paste0("PC", i)))
    }
  }
  # create new dataframe column to wrap on that includes the PC number and score
  pAll.l$PC.Score <- rep(
    x = paste0(score.df$PC, " ", sprintf("%1.3g", score.df$Score)),
    each = length(x = unique(x = pAll.l$Contig))
  )
  pAll.l$PC.Score <- factor(
    x = pAll.l$PC.Score,
    levels = paste0(score.df$PC, " ", sprintf("%1.3g", score.df$Score))
  )
  gp <- ggplot(data = pAll.l, mapping = aes(sample=Value)) +
    stat_qq(distribution = qunif) +
    facet_wrap("PC.Score", ncol = nCol) +
    labs(x = "Theoretical [runif(1000)]", y = "Empirical") +
    xlim(0, plot.y.lim) +
    ylim(0, plot.x.lim) +
    coord_flip() +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", na.rm = TRUE) +
    theme_bw()
  return(gp)
}

#' Scatter plot of single cell data
#'
#' Creates a scatter plot of two features (typically gene expression), across a
#' set of single cells. Cells are colored by their identity class.
#'
#' @param object Seurat object
#' @inheritParams FetchData
#' @param gene1 First feature to plot. Typically gene expression but can also
#' be metrics, PC scores, etc. - anything that can be retreived with FetchData
#' @param gene2 Second feature to plot.
#' @param cell.ids Cells to include on the scatter plot.
#' @param col.use Colors to use for identity class plotting.
#' @param pch.use Pch argument for plotting
#' @param cex.use Cex argument for plotting
#' @param use.imputed Use imputed values for gene expression (Default is FALSE)
#' @param use.scaled Use scaled data
#' @param use.raw Use raw data
#' @param do.hover Enable hovering over points to view information
#' @param data.hover Data to add to the hover, pass a character vector of features to add. Defaults to cell name
#' @param do.identify Opens a locator session to identify clusters of cells.
#' @param dark.theme Use a dark theme for the plot
#' @param do.spline Add a spline (currently hardwired to df=4, to be improved)
#' @param spline.span spline span in loess function call
#' @param \dots Additional arguments to be passed to plot.
#'
#' @return No return, only graphical output
#'
#' @export
#'
GenePlot <- function(
  object,
  gene1,
  gene2,
  cell.ids = NULL,
  col.use = NULL,
  pch.use = 16,
  cex.use = 1.5,
  use.imputed = FALSE,
  use.scaled = FALSE,
  use.raw = FALSE,
  do.hover = FALSE,
  data.hover = NULL,
  do.identify = FALSE,
  dark.theme = FALSE,
  do.spline = FALSE,
  spline.span = 0.75,
  ...

) {
  cell.ids <- set.ifnull(x = cell.ids, y = object@cell.names)
  #   Don't transpose the data.frame for better compatability with FeatureLocator and the rest of Seurat
  data.use <- as.data.frame(
    x = FetchData(
      object = object,
      vars.all = c(gene1, gene2),
      cells.use = cell.ids,
      use.imputed = use.imputed,
      use.scaled = use.scaled,
      use.raw = use.raw
    )
  )
  #   Ensure that our data is only the cells we're working with and
  #   the genes we want. This step seems kind of redundant though...
  data.plot <- data.use[cell.ids, c(gene1, gene2)]
  #   Set names to 'x' and 'y' for easy calling later on
  names(x = data.plot) <- c('x', 'y')
  ident.use <- as.factor(x = object@ident[cell.ids])
  if (length(x = col.use) > 1) {
    col.use <- col.use[as.numeric(x = ident.use)]
  } else {
    col.use <- set.ifnull(x = col.use,y = as.numeric(x = ident.use))
  }
  gene.cor <- round(x = cor(x = data.plot$x, y = data.plot$y), digits = 2)
  if (dark.theme) {
    par(bg = 'black')
    col.use <- sapply(
      X = col.use,
      FUN = function(color) ifelse(
        test = all(col2rgb(color) == 0),
        yes = 'white',
        no = color
      )
    )
    axes = FALSE
    col.lab = 'white'
  } else {
    axes = TRUE
    col.lab = 'black'
  }
  #   Plot the data
  plot(
    x = data.plot$x,
    y = data.plot$y,
    xlab = gene1,
    ylab = gene2,
    col = col.use,
    cex = cex.use,
    main = gene.cor,
    pch = pch.use,
    axes = axes,
    col.lab = col.lab,
    col.main = col.lab,
    ...
  )
  if (dark.theme) {
    axis(
      side = 1,
      at = NULL,
      labels = TRUE,
      col.axis = col.lab,
      col = col.lab
    )
    axis(
      side = 2,
      at = NULL,
      labels = TRUE,
      col.axis = col.lab,
      col = col.lab
    )
  }
  if (do.spline) {
    # spline.fit <- smooth.spline(x = g1, y = g2, df = 4)
    spline.fit <- smooth.spline(x = data.plot$x, y = data.plot$y, df = 4)
    #lines(spline.fit$x,spline.fit$y,lwd=3)
    #spline.fit=smooth.spline(g1,g2,df = 4)
    # loess.fit <- loess(formula = g2 ~ g1, span=spline.span)
    loess.fit <- loess(formula = y ~ x, data = data.plot, span = spline.span)
    #lines(spline.fit$x,spline.fit$y,lwd=3)
    # points(x = g1, y = loess.fit$fitted, col="darkblue")
    points(x = data.plot$x, y = loess.fit$fitted, col = 'darkblue')
  }
  if (do.identify | do.hover) {
    #   This is where that untransposed renamed data.frame comes in handy
    p <- ggplot2::ggplot(data = data.plot, mapping = aes(x = x, y = y))
    p <- p + geom_point(
      mapping = aes(color = colors),
      size = cex.use,
      shape = pch.use,
      color = col.use
    )
    p <- p + labs(title = gene.cor, x = gene1, y = gene2)
    if (do.hover) {
      names(x = data.plot) <- c(gene1, gene2)
      if (is.null(x = data.hover)) {
        features.info <- NULL
      } else {
        features.info <- FetchData(object = object, vars.all = data.hover)
      }
      return(HoverLocator(
        plot = p,
        data.plot = data.plot,
        features.info = features.info,
        dark.theme = dark.theme,
        title = gene.cor
      ))
    } else if (do.identify) {
      return(FeatureLocator(
        plot = p,
        data.plot = data.plot,
        dark.theme = dark.theme
      ))
    }
  }
}

# Documentation
###############
#' @export
#'
RemovePC <- function(object, pcs.remove, use.full = FALSE, ...) {
  data.old <- object@data
  pcs.use <- anotinb(x = 1:ncol(x = object@pca.obj[[1]]$rotation), y = pcs.remove)
  if (use.full) {
    data.x <- as.matrix(
      x = object@pca.x.full[, ainb(
        a = pcs.use,
        b = 1:ncol(x = object@pca.x.full)
      )]
    )
  } else {
    data.x <- as.matrix(x = object@pca.obj[[1]]$x[, pcs.use])
  }
  data.1 <- data.x %*% t(x = as.matrix(x = object@pca.obj[[1]]$rotation[, pcs.use]))
  data.2 <- sweep(
    x = data.1,
    MARGIN = 2,
    STATS = colMeans(x = object@scale.data),
    FUN = "+"
  )
  data.3 <- sweep(
    x = data.2,
    MARGIN = 1,
    STATS = apply(X = object@data[rownames(x = data.2), ], MARGIN = 1, FUN = sd),
    FUN = "*"
  )
  data.3 <- sweep(
    X = data.3,
    MARGIN = 1,
    STATS = apply(X = object@data[rownames(x = data.2), ], MARGIN = 1, FUN = mean),
    FUN = "+"
  )
  object@scale.data <- (data.2)
  data.old <- data.old[rownames(x = data.3), ]
  data.4 <- data.3
  data.4[data.old == 0] <- 0
  data.4[data.4 < 0] <- 0
  object@data[rownames(x = data.4), ] <- data.frame(data.4)
  return(object)
}

# Documentation
###############
GeneScorePlot <- function(
  object,
  gene1,
  score.name,
  cell.ids = NULL,
  col.use = NULL,
  nrpoints.use = Inf,
  pch.use = 16,
  cex.use = 2,
  ...
) {
  cell.ids <- set.ifnull(x = cell.ids, y = colnames(x = object@data))
  g1 <- as.numeric(x = object@data[gene1, cell.ids])
  my.score <- retreiveScore(object, score.name)
  s1 <- as.numeric(x = my.score[cell.ids])
  col.use <- set.ifnull(as.numeric(x = as.factor(x = object@ident[cell.ids])))
  gene.cor <- round(x = cor(x = g1, y = s1), digits = 2)
  smoothScatter(
    x = g1,
    y = s1,
    xlab = gene1,
    ylab = score.name,
    col = col.use,
    nrpoints = nrpoints.use,
    cex = cex.use,
    main = gene.cor,
    pch = pch.use
  )
}

#' Cell-cell scatter plot
#'
#' Creates a plot of scatter plot of genes across two single cells
#'
#' @param object Seurat object
#' @param cell1 Cell 1 name (can also be a number, representing the position in
#' object@@cell.names)
#' @param cell2 Cell 2 name (can also be a number, representing the position in
#' object@@cell.names)
#' @param gene.ids Genes to plot (default, all genes)
#' @param col.use Colors to use for the points
#' @param nrpoints.use Parameter for smoothScatter
#' @param pch.use Point symbol to use
#' @param cex.use Point size
#' @param do.hover Enable hovering over points to view information
#' @param do.identify Opens a locator session to identify clusters of cells.
#' points to reveal gene names (hit ESC to stop)
#' @param \dots Additional arguments to pass to smoothScatter
#'
#' @return No return value (plots a scatter plot)
#'
#' @export
#'
CellPlot <- function(
  object,
  cell1,
  cell2,
  gene.ids = NULL,
  col.use = "black",
  nrpoints.use = Inf,
  pch.use = 16,
  cex.use = 0.5,
  do.hover = FALSE,
  do.identify = FALSE,
  ...
) {
  gene.ids <- set.ifnull(x = gene.ids, y = rownames(x = object@data))
  #   Transpose this data.frame so that the genes are in the row for
  #   easy selecting with do.identify
  data.plot <- as.data.frame(
    x = t(
      x = FetchData(
        object = object,
        vars.all = gene.ids,
        cells.use = c(cell1, cell2)
      )
    )
  )
  #   Set names for easy calling with ggplot
  names(x = data.plot) <- c('x', 'y')
  gene.cor <- round(x = cor(x = data.plot$x, y = data.plot$y), digits = 2)
  smoothScatter(
    x = data.plot$x,
    y = data.plot$y,
    xlab = cell1,
    ylab = cell2,
    col = col.use,
    nrpoints = nrpoints.use,
    pch = pch.use,
    cex = cex.use,
    main = gene.cor
  )
  if (do.identify | do.hover) {
    #   This is where that untransposed renamed data.frame comes in handy
    p <- ggplot2::ggplot(data = data.plot, mapping = aes(x = x, y = y))
    p <- p + geom_point(
      mapping = aes(color = colors),
      size = cex.use,
      shape = pch.use,
      color = col.use
    )
    p <- p + labs(title = gene.cor, x = cell1, y = cell2)
    if (do.hover) {
      names(x = data.plot) <- c(cell1, cell2)
      return(HoverLocator(plot = p, data.plot = data.plot, title = gene.cor))
    } else if (do.identify) {
      return(FeatureLocator(plot = p, data.plot = data.plot, ...))
    }
  }
}

# Documentation
###############
#' @export
#'
JackStrawPermutationTest <- function(
  object,
  genes.use = NULL,
  num.iter = 100,
  thresh.use = 0.05,
  do.print = TRUE,
  k.seed = 1
) {
  genes.use <- set.ifnull(x = genes.use, y = rownames(x = object@pca.x))
  genes.use <- ainb(a = genes.use, b = rownames(x = object@scale.data))
  data.use <- t(x = as.matrix(x = object@scale.data[genes.use, ]))
  if (do.print) {
    print(paste("Running", num.iter, "iterations"))
  }
  pa.object <- permutationPA(
    data.use,
    B = num.iter,
    threshold = thresh.use,
    verbose = do.print,
    seed = k.seed
  )
  if (do.print) {
    cat("\n\n")
  }
  if (do.print) {
    print(paste("JackStraw returns", pa.object$r, "significant components"))
  }
  return(pa.object)
}

# Documentation
###############
#multicore version of jackstraw
#DOES NOT WORK WITH WINDOWS
#' @export
#'
JackStrawMC <- function(
  object,
  num.pc = 30,
  num.replicate = 100,
  prop.freq = 0.01,
  do.print = FALSE,
  num.cores = 8
) {
  pc.genes <- rownames(x = object@pca.x)
  if (length(x = pc.genes) < 200) {
    prop.freq <- max(prop.freq, 0.015)
  }
  md.x <- as.matrix(x = object@pca.x)
  md.rot <- as.matrix(x = object@pca.rot)
  if (do.print) {
    fake.pcVals.raw <- mclapply(
      X = 1:num.replicate,
      FUN = function(x) {
        print(x)
        return(jackRandom(
          scaled.data = object@scale.data[pc.genes, ],
          prop = prop.freq,
          r1.use = 1,
          r2.use = num.pc,
          seed.use = x
        ))
      },
      mc.cores = num.cores
    )
  } else {
    fake.pcVals.raw <- mclapply(
      X = 1:num.replicate,
      FUN = function(x) {
        return(jackRandom(
          scaled.data = object@scale.data[pc.genes, ],
          prop = prop.freq,
          r1.use = 1,
          r2.use = num.pc,
          seed.use=x
        ))
      }, mc.cores = num.cores
    )
  }
  fake.pcVals <- simplify2array(
    x = mclapply(
      X = 1:num.pc,
      FUN = function(x) {
        return(as.numeric(x = unlist(x = lapply(
          X = 1:num.replicate,
          FUN = function(y) {
            return(fake.pcVals.raw[[y]][, x])
          }
        ))))
      },
      mc.cores = num.cores
    )
  )
  object@jackStraw.fakePC <- data.frame(fake.pcVals)
  object@jackStraw.empP <- data.frame(
    simplify2array(
      x = mclapply(
        X = 1:num.pc,
        FUN = function(x) {
          return(unlist(x = lapply(
            X = abs(md.x[, x]),
            FUN = empP,
            nullval = abs(x = fake.pcVals[, x])
          )))
        },
        mc.cores = num.cores
      )
    )
  )
  colnames(x = object@jackStraw.empP) <- paste0("PC", 1:ncol(x = object@jackStraw.empP))
  return(object)
}

#' Determine statistical significance of PCA scores.
#'
#' Randomly permutes a subset of data, and calculates projected PCA scores for
#' these 'random' genes. Then compares the PCA scores for the 'random' genes
#' with the observed PCA scores to determine statistical signifance. End result
#' is a p-value for each gene's association with each principal component.
#'
#' @param object Seurat object
#' @param num.pc Number of PCs to compute significance for
#' @param num.replicate Number of replicate samplings to perform
#' @param prop.freq Proportion of the data to randomly permute for each
#' replicate
#' @param do.print Print the number of replicates that have been processed.
#' @param rev.pca By default computes the PCA on the cell x gene matrix. Setting to
#' true will compute it on gene x cell matrix. This should match what was set when the intial PCA was run.
#' @param do.fast Compute the PCA using the fast approximate calculation from the IRLBA package. Values stored with object
#' must also have been computed using the PCAFast() function.
#'
#' @return Returns a Seurat object where object@@jackStraw.empP represents
#' p-values for each gene in the PCA analysis. If ProjectPCA is subsequently
#' run, object@@jackStraw.empP.full then represents p-values for all genes.
#'
#' @importFrom pbapply pbsapply
#'
#' @references Inspired by Chung et al, Bioinformatics (2014)
#'
#' @export
#'
JackStraw <- function(
  object,
  num.pc = 30,
  num.replicate = 100,
  prop.freq = 0.01,
  do.print = FALSE,
  rev.pca = FALSE,
  do.fast = FALSE
) {
  # check that PCA calculation method matches
  if (do.fast) {
    if (is.null(x = object@pca.obj[[1]]$d)) {
      stop("For fast JackStraw, store PCA values computed with PCAFast()")
    }
  } else {
    if (is.null(x = object@pca.obj[[1]]$sdev)) {
      stop("For regular JackStraw, store PCA values computed with PCA()")
    }
  }
  # error checking for number of PCs
  if (num.pc > ncol(x = object@pca.rot)) {
    num.pc <- ncol(x = object@pca.rot)
    warning("Number of PCs specified is greater than PCs available. Setting num.pc to ", num.pc, " and continuing.")
  }
  if (num.pc > length(x = object@cell.names)) {
    num.pc <- length(x = object@cell.names)
    warning("Number of PCs specified is greater than number of cells. Setting num.pc to ", num.pc, " and continuing.")
  }
  pc.genes <- rownames(x = object@pca.x)
  if (length(x = pc.genes) < 3) {
    stop("Too few variable genes")
  }
  if (length(x = pc.genes) * prop.freq < 3) {
    warning(
      "Number of variable genes given ",
      prop.freq,
      " as the prop.freq is low. Consider including more variable genes and/or increasing prop.freq. ",
      "Continuing with 3 genes in every random sampling."
    )
  }
  md.x <- as.matrix(x = object@pca.x)
  md.rot <- as.matrix(x = object@pca.rot)
  # if (do.print) {
  #   fake.pcVals.raw <- pbsapply(
  #     X = 1:num.replicate,
  #     FUN = function(x){
  #       return(jackRandom(
  #         scaled.data = object@scale.data[pc.genes, ],
  #         prop = prop.freq,
  #         r1.use = 1,
  #         r2.use = num.pc,
  #         seed.use = x,
  #         rev.pca = rev.pca,
  #         do.fast = do.fast
  #       ))
  #     },
  #     simplify = FALSE
  #   )
  # } else {
  #   fake.pcVals.raw <- sapply(
  #     X = 1:num.replicate,
  #     FUN = function(x)
  #       return(jackRandom(
  #         scaled.data = object@scale.data[pc.genes, ],
  #         prop = prop.freq,
  #         r1.use = 1,
  #         r2.use = num.pc,
  #         seed.use = x,
  #         rev.pca = rev.pca,
  #         do.fast = do.fast
  #       )),
  #     simplify = FALSE
  #   )
  # }
  if (do.print) {
    applyFunction <- pbsapply
  } else {
    applyFunction <- sapply
  }
  fake.pcVals.raw <- applyFunction(
    X = 1:num.replicate,
    FUN = function(x)
      return(jackRandom(
        scaled.data = object@scale.data[pc.genes, ],
        prop = prop.freq,
        r1.use = 1,
        r2.use = num.pc,
        seed.use = x,
        rev.pca = rev.pca,
        do.fast = do.fast
      )),
    simplify = FALSE
  )
  fake.pcVals <- sapply(
    X = 1:num.pc,
    FUN = function(x) {
      return(as.numeric(x = unlist(x = lapply(
        X = 1:num.replicate,
        FUN = function(y) {
          return(fake.pcVals.raw[[y]][, x])
        }
      ))))
    }
  )
  object@jackStraw.fakePC = data.frame(fake.pcVals)
  object@jackStraw.empP <- data.frame(
    sapply(
      X = 1:num.pc,
      FUN = function(x) {
        return(unlist(x = lapply(
          X = abs(md.x[, x]),
          FUN = empP,
          nullval = abs(fake.pcVals[,x])
        )))
      }
    )
  )
  colnames(x = object@jackStraw.empP) <- paste0("PC", 1:ncol(x = object@jackStraw.empP))
  return(object)
}

# Documentatin
##############
#' @export
#'
jackRandom <- function(
  scaled.data,
  prop.use = 0.01,
  r1.use = 1,
  r2.use = 5,
  seed.use = 1,
  rev.pca = FALSE,
  do.fast = FALSE
) {
  set.seed(seed = seed.use)
  rand.genes <- sample(
    x = rownames(x = scaled.data),
    size = nrow(x = scaled.data) * prop.use
  )
  # make sure that rand.genes is at least 3
  if (length(x = rand.genes) < 3){
    rand.genes <- sample(x = rownames(x = scaled.data), size = 3)
  }
  data.mod <- scaled.data
  data.mod[rand.genes, ] <- shuffleMatRow(x = scaled.data[rand.genes, ])
  if (rev.pca){
    if (do.fast) {
      fake.pca <- irlba(A = data.mod, nv = r2.use)
      fake.rot <- fake.pca$v[, 1:r2.use]
      rownames(x = fake.rot) <- colnames(x = data.mod)
      colnames(x = fake.rot) <- paste0("PC", 1:r2.use)
      fake.x <- fake.pca$u[, 1:r2.use]
      rownames(x = fake.x) <- rownames(x = data.mod)
      colnames(x = fake.x) <- colnames(x = fake.rot)
    } else {
      fake.pca <- prcomp(x = data.mod)
      fake.x <- fake.pca$x
      fake.rot <- fake.pca$rotation
    }
  } else {
    data.mod <- t(x = data.mod)
    if (do.fast) {
      fake.pca <- irlba(A = data.mod, nv = r2.use)
      fake.rot <- fake.pca$u[, 1:r2.use]
      rownames(x = fake.rot) <- rownames(x = data.mod)
      colnames(x = fake.rot) <- paste0("PC", 1:r2.use)
      fake.x <- fake.pca$v[, 1:r2.use]
      rownames(x = fake.x) <- colnames(x = data.mod)
      colnames(x = fake.x) <- colnames(x = fake.rot)
    } else {
      fake.pca <- prcomp(x = data.mod)
      fake.x <- fake.pca$rotation
      fake.rot <- fake.pca$x
    }
  }
  return(fake.x[rand.genes, r1.use:r2.use])
}

# Documentation
###############
JackStrawFull <- function(
  object,
  num.pc = 5,
  num.replicate = 100,
  prop.freq = 0.01
) {
  pc.genes <- rownames(x = object@pca.x)
  if (length(x = pc.genes) < 200) {
    prop.freq <- max(prop.freq, 0.015)
  }
  md.x <- as.matrix(x = object@pca.x)
  md.rot <- as.matrix(x = object@pca.rot)
  real.fval <- sapply(
    X = 1:num.pc,
    FUN = function(x) {
      return(unlist(x = lapply(
        X = pc.genes,
        FUN = jackF,
        r1 = x,
        r2 = x,
        x = md.x,
        rot = md.rot
      )))
    }
  )
  rownames(x = real.fval) <- pc.genes
  object@real.fval <- data.frame(real.fval)
  fake.fval <- sapply(
    X = 1:num.pc,
    FUN = function(x) {
      return(unlist(x = replicate(
        n = num.replicate,
        expr = jackStrawF(
          prop = prop.freq,
          data = object@scale.data[pc.genes, ],
          myR1 = x,
          myR2 = x
        ),
        simplify = FALSE
      )))
    }
  )
  rownames(x = fake.fval) <- 1:nrow(x = fake.fval)
  object@fake.fval <- data.frame(fake.fval)
  object@emp.pval <- data.frame(
    sapply(
      X = 1:num.pc,
      FUN = function(x) {
        return(unlist(x = lapply(
          X = object@real.fval[, x],
          FUN = empP,
          nullval = object@fake.fval[, x]
        )))
      }
    )
  )
  rownames(x = object@emp.pval) <- pc.genes
  colnames(x = object@emp.pval) <- paste0("PC", 1:ncol(x = object@emp.pval),)
  return(object)
}

#' Identify variable genes
#'
#' Identifies genes that are outliers on a 'mean variability plot'. First, uses
#' a function to calculate average expression (fxn.x) and dispersion (fxn.y)
#' for each gene. Next, divides genes into num.bin (deafult 20) bins based on
#' their average expression, and calculates z-scores for dispersion within each
#' bin. The purpose of this is to identify variable genes while controlling for
#' the strong relationship between variability and average expression.
#'
#' Exact parameter settings may vary empirically from dataset to dataset, and
#' based on visual inspection of the plot.
#' Setting the y.cutoff parameter to 2 identifies genes that are more than two standard
#' deviations away from the average dispersion within a bin. The default X-axis function
#' is the mean expression level, and for Y-axis it is the log(Variance/mean). All mean/variance
#' calculations are not performed in log-space, but the results are reported in log-space -
#' see relevant functions for exact details.
#'
#' @param object Seurat object
#' @param fxn.x Function to compute x-axis value (average expression). Default
#' is to take the mean of the detected (i.e. non-zero) values
#' @param fxn.y Function to compute y-axis value (dispersion). Default is to
#' take the standard deviation of all values/
#' @param do.plot Plot the average/dispersion relationship
#' @param set.var.genes Set object@@var.genes to the identified variable genes
#' (default is TRUE)
#' @param do.text Add text names of variable genes to plot (default is TRUE)
#' @param x.low.cutoff Bottom cutoff on x-axis for identifying variable genes
#' @param x.high.cutoff Top cutoff on x-axis for identifying variable genes
#' @param y.cutoff Bottom cutoff on y-axis for identifying variable genes
#' @param y.high.cutoff Top cutoff on y-axis for identifying variable genes
#' @param cex.use Point size
#' @param cex.text.use Text size
#' @param do.spike FALSE by default. If TRUE, color all genes starting with ^ERCC a different color
#' @param pch.use Pch value for points
#' @param col.use Color to use
#' @param spike.col.use if do.spike, color for spike-in genes
#' @param plot.both Plot both the scaled and non-scaled graphs.
#' @param do.contour Draw contour lines calculated based on all genes
#' @param contour.lwd Contour line width
#' @param contour.col Contour line color
#' @param contour.lty Contour line type
#' @param num.bin Total number of bins to use in the scaled analysis (default
#' is 20)
#' @param do.recalc TRUE by default. If FALSE, plots and selects variable genes without recalculating statistics for each gene.
#' @param sort.results If TRUE (by default), sort results in object@mean.var in decreasing order of dispersion
#'
#' @importFrom MASS kde2d
#'
#' @return Returns a Seurat object, placing variable genes in object@@var.genes.
#' The result of all analysis is stored in object@@mean.var
#'
#' @export
#'
MeanVarPlot <- function(
  object,
  fxn.x = expMean,
  fxn.y = logVarDivMean,
  do.plot = TRUE,
  set.var.genes = TRUE,
  do.text = TRUE,
  x.low.cutoff = 0.1,
  x.high.cutoff = 8,
  y.cutoff = 1,
  y.high.cutoff = Inf,
  cex.use = 0.5,
  cex.text.use = 0.5,
  do.spike = FALSE,
  pch.use=16,
  col.use = "black",
  spike.col.use = "red",
  plot.both = FALSE,
  do.contour = TRUE,
  contour.lwd = 3,
  contour.col = "white",
  contour.lty = 2,
  num.bin = 20,
  do.recalc = TRUE,
  sort.results=TRUE
) {
  data <- object@data
  if (do.recalc) {
    genes.use <- rownames(x = object@data)
    data.x <- rep(x = 0, length(x = genes.use))
    names(x = data.x) <- genes.use
    data.y <- data.x
    data.norm.y <- data.x
    bin.size <- 1000
    max.bin <- floor(x = length(x = genes.use) / bin.size) + 1
    print("Calculating gene dispersion")
    pb <- txtProgressBar(min = 0, max = max.bin, style = 3)
    for (i in 1:max.bin) {
      my.inds <- ((bin.size * (i - 1)):(bin.size * i - 1)) + 1
      my.inds <- my.inds[my.inds <= length(x = genes.use)]
      genes.iter <- genes.use[my.inds]
      data.iter <- data[genes.iter, ]
      data.x[genes.iter] <- apply(X = data.iter, MARGIN = 1, FUN = fxn.x)
      data.y[genes.iter] <- apply(X = data.iter, MARGIN = 1, FUN = fxn.y)
      setTxtProgressBar(pb = pb, value = i)
    }
    close(con = pb)
    data.y[is.na(x = data.y)] <- 0
    data.x[is.na(x = data.x)] <- 0
    data_x_bin <- cut(x = data.x, breaks = num.bin)
    names(x = data_x_bin) <- names(x = data.x)
    mean_y <- tapply(X = data.y, INDEX = data_x_bin, FUN = mean)
    sd_y <- tapply(X = data.y, INDEX = data_x_bin, FUN = sd)
    data.norm.y <- (data.y - mean_y[as.numeric(x = data_x_bin)]) /
      sd_y[as.numeric(x = data_x_bin)]
    #data.x=apply(data,1,fxn.x); data.y=apply(data,1,fxn.y); data.x[is.na(data.x)]=0
    #data.norm.y=meanNormFunction(data,fxn.x,fxn.y,num.bin)
    data.norm.y[is.na(x = data.norm.y)] <- 0
    names(x = data.norm.y) <- names(x = data.x)
    mv.df <- data.frame(data.x, data.y, data.norm.y)
    rownames(x = mv.df) <- rownames(x = data)
    object@mean.var <- mv.df
  }
  data.x <- object@mean.var[, 1]
  data.y <- object@mean.var[, 2]
  data.norm.y <- object@mean.var[, 3]
  names(x = data.x) <- names(x = data.y) <- names(x = data.norm.y) <- rownames(x = object@data)
  pass.cutoff <- names(x = data.x)[which(
    x = (
      (data.x > x.low.cutoff) & (data.x < x.high.cutoff)
    ) &
      (data.norm.y > y.cutoff) &
      (data.norm.y < y.high.cutoff)
  )]
  if (do.spike) {
    spike.genes <- rownames(x = subr(data = data, code = "^ERCC"))
  }
  if (do.plot) {
    if (plot.both) {
      par(mfrow = c(1, 2))
      smoothScatter(
        x = data.x,
        y = data.y,
        pch = pch.use,
        cex = cex.use,
        col = col.use,
        xlab = "Average expression",
        ylab = "Dispersion",
        nrpoints = Inf
      )
      if (do.contour) {
        data.kde <- kde2d(x = data.x, y = data.y)
        contour(
          x = data.kde,
          add = TRUE,
          lwd = contour.lwd,
          col = contour.col,
          lty = contour.lty
        )
      }
      if (do.spike) {
        points(
          x = data.x[spike.genes],
          y = data.y[spike.genes],
          pch = 16,
          cex = cex.use,
          col = spike.col.use
        )
      }
      if (do.text) {
        text(
          x = data.x[pass.cutoff],
          y = data.y[pass.cutoff],
          labels = pass.cutoff,
          cex = cex.text.use
        )
      }
    }
    smoothScatter(
      x = data.x,
      y = data.norm.y,
      pch = pch.use,
      cex = cex.use,
      col = col.use,
      xlab = "Average expression",
      ylab = "Dispersion",
      nrpoints = Inf
    )
    if (do.contour) {
      data.kde <- kde2d(x = data.x, y = data.norm.y)
      contour(
        x = data.kde,
        add = TRUE,
        lwd = contour.lwd,
        col = contour.col,
        lty = contour.lty
      )
    }
    if (do.spike) {
      points(
        x = data.x[spike.genes],
        y = data.norm.y[spike.genes],
        pch = 16,
        cex = cex.use,
        col = spike.col.use,
        nrpoints = Inf
      )
    }
    if (do.text) {
      text(
        x = data.x[pass.cutoff],
        y = data.norm.y[pass.cutoff],
        labels = pass.cutoff,
        cex = cex.text.use
      )
    }
  }
  if (set.var.genes) {
    object@var.genes <- pass.cutoff
    if (sort.results) {
      object@mean.var=object@mean.var[order(object@mean.var$data.y,decreasing = T),]
    }
    return(object)
  } else {
    return(pass.cutoff)
  }
}

#' Calculate enrichment scores for gene expression programs in single cells
#'
#' Calculate the average expression levels of each program (cluster) on single cell level,
#' subtracted by the aggregated expression of control gene sets.
#' All analyzed genes are binned based on averaged expression, and the control genes are
#' randomly selected from each bin.
#'
#' @param object Seurat object
#' @param genes.list Gene expression programs in list
#' @param n.bin Number of bins of aggregate expression levels for all analyzed genes
#' @param seed.use Random seed for sampling
#' @param ctrl.size Number of control genes selected from the same bin per analyzed gene
#' @param use.k Use gene clusters returned from DoKMeans()
#' @param enrich.name Name for the expression programs
#'
#' @return Returns a Seurat object with enrichment scores added to object@data.info
#'
#' @importFrom Hmisc cut2
#'
#' @references Tirosh et al, Science (2016)
#'
#' @export
#'
AddEnrichScore <- function(
  object,
  genes.list = NULL,
  n.bin = 20,
  seed.use = 1,
  ctrl.size = 20,
  use.k = FALSE,
  enrich.name = "Cluster"
) {
  if (use.k) {
    genes.k <- names(x = object@kmeans.obj[[1]]$cluster)
    genes.list <- list()
    for (i in as.numeric(x = names(x = table(object@kmeans.obj[[1]]$cluster)))) {
      genes.list[[i]] <- names(x = which(x = object@kmeans.obj[[1]]$cluster == i))
    }
    cluster.length <- length(x = genes.list)
  } else {
    if (is.null(x = genes.list)) {
      stop("Missing input gene list")
    }
    genes.k <- ainb(
      a = unique(x = unlist(x = genes.list)),
      b = rownames(x = object@data)
    )
    cluster.length <- length(x = genes.list)
  }
  data.avg <- apply(X = object@data[genes.k, ], MARGIN = 1, FUN = mean)
  data.avg <- data.avg[order(data.avg)]
  data.cut <- as.numeric(x = cut2(
    x = data.avg,
    m = round(x = length(x = data.avg) / n.bin)
  ))
  names(x = data.cut) <- names(x = data.avg)
  ctrl.use <- vector("list", cluster.length)
  for (i in 1:cluster.length) {
    genes.use <- genes.list[[i]]
    for (j in 1:length(x = genes.use)) {
      ctrl.use[[i]] <- c(
        ctrl.use[[i]],
        names(x = sample(
          x = data.cut[which(x = data.cut == data.cut[genes.use[j]])],
          size = ctrl.size,
          replace = FALSE
        ))
      )
    }
  }
  ctrl.use <- lapply(X = ctrl.use, FUN = unique)
  ctrl.scores <- c()
  for (i in 1:length(ctrl.use)) {
    genes.use <- ctrl.use[[i]]
    ctrl.scores <- rbind(
      ctrl.scores,
      apply(X = object@data[genes.use, ], MARGIN = 2, FUN = mean)
    )
  }
  genes.scores <- c()
  for (i in 1:cluster.length) {
    genes.use <- genes.list[[i]]
    genes.scores <- rbind(
      genes.scores,
      apply(X = object@data[genes.use, ], MARGIN = 2, FUN = mean)
    )
  }

  genes.scores.use <- genes.scores - ctrl.scores
  rownames(x = genes.scores.use) <- paste0(enrich.name, 1:cluster.length)
  genes.scores.use <- t(x = as.data.frame(x = genes.scores.use))
  object <- AddMetaData(
    object = object,
    metadata = genes.scores.use,
    col.name = colnames(x = genes.scores.use)
  )
  return (object)
}
