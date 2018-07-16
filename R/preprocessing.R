#' Initialize and setup the Seurat object
#'
#' Initializes the Seurat object and some optional filtering
#' @param raw.data Raw input data
#' @param project Project name (string)
#' @param min.cells Include genes with detected expression in at least this
#' many cells. Will subset the raw.data matrix as well. To reintroduce excluded
#' genes, create a new object with a lower cutoff.
#' @param min.genes Include cells where at least this many genes are detected.
#' @param is.expr Expression threshold for 'detected' gene. For most datasets, particularly UMI
#' datasets, will be set to 0 (default). If not, when initializing, this should be set to a level
#' based on pre-normalized counts (i.e. require at least 5 counts to be treated as expresesd) All
#' values less than this will be set to 0 (though maintained in object@raw.data).
#' @param normalization.method Method for cell normalization. Default is no normalization.
#' In this case, run NormalizeData later in the workflow. As a shortcut, you can specify a
#' normalization method (i.e. LogNormalize) here directly.
#' @param scale.factor If normalizing on the cell level, this sets the scale factor.
#' @param do.scale In object@@scale.data, perform row-scaling (gene-based
#' z-score). FALSE by default. In this case, run ScaleData later in the workflow. As a shortcut, you
#' can specify do.scale = TRUE (and do.center = TRUE) here.
#' @param do.center In object@@scale.data, perform row-centering (gene-based centering)
#' @param names.field For the initial identity class for each cell, choose this field from the
#' cell's column name
#' @param names.delim For the initial identity class for each cell, choose this delimiter from the
#' cell's column name
#' @param meta.data Additional metadata to add to the Seurat object. Should be a data frame where
#' the rows are cell names, and the columns are additional metadata fields
#' @param verbose display progress bar for normalization and/or scaling procedure.
#' @param ... Ignored
#'
#' @return Returns a Seurat object with the raw data stored in object@@raw.data.
#' object@@data, object@@meta.data, object@@ident, also initialized.
#'
#' @importFrom methods new
#' @importFrom utils packageVersion
#' @importFrom Matrix colSums rowSums
#'
#' @export
#'
#' @examples
#' pbmc_raw <- read.table(
#'   file = system.file('extdata', 'pbmc_raw.txt', package = 'Seurat'),
#'   as.is = TRUE
#' )
#' pbmc_small <- CreateSeuratObject(raw.data = pbmc_raw)
#' pbmc_small
#'
CreateSeuratObject <- function(
  raw.data,
  project = "SeuratProject",
  min.cells = 0,
  min.genes = 0,
  is.expr = 0,
  normalization.method = NULL,
  scale.factor = 1e4,
  do.scale = FALSE,
  do.center = FALSE,
  names.field = 1,
  names.delim = "_",
  meta.data = NULL,
  verbose = TRUE,
  ...
) {
  seurat.version <- packageVersion("Seurat")
  object <- new(
    Class = "seurat",
    raw.data = raw.data,
    is.expr = is.expr,
    project.name = project,
    version = seurat.version
  )
  # filter cells on number of genes detected
  # modifies the raw.data slot as well now
  object.raw.data <- object@raw.data
  if (is.expr > 0) {
    # suppress Matrix package note:
    # Note: method with signature 'CsparseMatrix#Matrix#missing#replValue' chosen for function '[<-',
    # target signature 'dgCMatrix#lgeMatrix#missing#numeric'.
    # "Matrix#ldenseMatrix#missing#replValue" would also be valid
    suppressMessages(expr = object.raw.data[object.raw.data < is.expr] <- 0)
  }
  num.genes <- colSums(object.raw.data > is.expr)
  num.mol <- colSums(object.raw.data)
  cells.use <- names(x = num.genes[which(x = num.genes > min.genes)])
  object@raw.data <- object@raw.data[, cells.use]
  object@data <- object.raw.data[, cells.use]
  # filter genes on the number of cells expressing
  # modifies the raw.data slot as well now
  genes.use <- rownames(object@data)
  if (min.cells > 0) {
    num.cells <- rowSums(object@data > 0)
    genes.use <- names(x = num.cells[which(x = num.cells >= min.cells)])
    object@raw.data <- object@raw.data[genes.use, ]
    object@data <- object@data[genes.use, ]
  }
  object@ident <- factor(x = unlist(x = lapply(
    X = colnames(x = object@data),
    FUN = ExtractField,
    field = names.field,
    delim = names.delim
  )))
  names(x = object@ident) <- colnames(x = object@data)
  object@cell.names <- names(x = object@ident)
  # if there are more than 100 idents, set all idents to project name
  ident.levels <- length(x = unique(x = object@ident))
  if ((ident.levels > 100 || ident.levels == 0) || ident.levels == length(x = object@ident)) {
    object <- SetIdent(object, ident.use = project)
  }
  nGene <- num.genes[cells.use]
  nUMI <- num.mol[cells.use]
  object@meta.data <- data.frame(nGene, nUMI)
  if (!is.null(x = meta.data)) {
    object <- AddMetaData(object = object, metadata = meta.data)
  }
  object@meta.data[["orig.ident"]] <- NULL
  object@meta.data[names(object@ident), "orig.ident"] <- object@ident
  if (!is.null(normalization.method)) {
    object <- NormalizeData(
      object = object,
      assay.type = "RNA",
      normalization.method = normalization.method,
      scale.factor = scale.factor,
      verbose = verbose
    )
  }
  if (do.scale | do.center) {
    object <- ScaleData(
      object = object,
      do.scale = do.scale,
      do.center = do.center,
      verbose = verbose
    )
  }
  spatial.obj <- new(
    Class = "spatial.info",
    mix.probs = data.frame(nGene)
  )
  object@spatial <- spatial.obj
  parameters.to.store <- as.list(x = environment(), all = TRUE)[names(formals("CreateSeuratObject"))]
  parameters.to.store$raw.data <- NULL
  parameters.to.store$meta.data <- NULL
  object <- SetCalcParams(
    object = object,
    calculation = "CreateSeuratObject",
    ... = parameters.to.store
  )

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
#' @examples
#' \dontrun{
#' data_dir <- 'path/to/data/directory'
#' list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
#' expression_matrix <- Read10X(data.dir = data_dir)
#' seurat_object = CreateSeuratObject(raw.data = expression_matrix)
#' }
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
            FUN = ExtractField, field = 1,
            delim = "-"
          )
        )
      )
    }
    rownames(x = data) <- make.unique(
      names = as.character(
        x = sapply(
          X = gene.names,
          FUN = ExtractField,
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

#' Read 10X hdf5 file
#'
#' Read gene expression matrix from 10X CellRanger hdf5 file
#'
#' @param filename Path to h5 file
#' @param ensg.names Label row names with ENSG names rather than unique gene names
#'
#' @return Returns a sparse matrix with rows and columns labeled. If multiple genomes are present,
#' returns a list of sparse matrices (one per genome).
#'
#' @importFrom hdf5r H5File
#'
#' @export
#'
Read10X_h5 <- function(filename, ensg.names = FALSE){
  if(!file.exists(filename)){
    stop("File not found")
  }
  infile <- H5File$new(filename)
  genomes <- names(infile)
  output <- list()
  for(genome in genomes){
    raw.data <- infile[[paste0(genome, '/data')]]
    indices <- infile[[paste0(genome, '/indices')]]
    indptr <- infile[[paste0(genome, '/indptr')]]
    shp <- infile[[paste0(genome, '/shape')]]
    if(ensg.names){
      gene_names <- infile[[paste0(genome, '/genes')]][]
    } else {
      gene_names <- make.unique(infile[[paste0(genome, '/gene_names')]][])
    }
    barcodes <- infile[[paste0(genome, '/barcodes')]]
    sparse.mat <- sparseMatrix(i = indices[] + 1, p = indptr[],
                               x = as.numeric(raw.data[]),
                               dims = shp[], giveCsparse = FALSE)
    rownames(sparse.mat) <- gene_names
    colnames(sparse.mat) <- barcodes[]
    output[[genome]] <- sparse.mat
  }
  infile$close_all()
  if(length(output) == 1) {
    return(output[[genome]])
  } else{
    return(output)
  }
}

#' @export
#'
NormalizeData.default <- function(
  object,
  normalization.method = "LogNormalize",
  scale.factor = 1e4,
  verbose = TRUE,
  ...
) {
  if (is.null(x = normalization.method)) {
    return(object)
  }
  normalized.data <- switch(
    EXPR = normalization.method,
    'LogNormalize' = LogNormalize(
      data = object,
      scale.factor = scale.factor,
      verbose = verbose
    ),
    'genesCLR' = CustomNormalize(
      data = object,
      custom_function = function(x) {
        return(log1p(x = x / (exp(x = sum(log1p(x = x[x > 0]), na.rm = TRUE) / length(x = x + 1)))))
      },
      across = 'genes'
    ),
    stop("Unkown normalization method: ", normalization.method)
  )
  return(normalized.data)
}

#' @describeIn NormalizeData Normalize data in an Assay object
#' @export
#' @method NormalizeData Assay
#'
NormalizeData.Assay <- function(
  object,
  normalization.method = "LogNormalize",
  scale.factor = 1e4,
  verbose = TRUE,
  ...
) {
  object <- SetAssayData(
    object = object,
    slot = 'data',
    new.data = NormalizeData(
      object = GetAssayData(object = object, slot = 'raw.data'),
      normalization.method = normalization.method,
      scale.factor = scale.factor,
      verbose = verbose,
      ...
    )
  )
  return(object)
}

#' @param assay.use Name of assay to use
#'
#' @describeIn NormalizeData Normalize data in a Seurat object
#' @export
#' @method NormalizeData Seurat
#'
#' @examples
#' \dontrun{
#' pbmc_small
#' pmbc_small <- NormalizeData(object = pbmc_small)
#' }
#'
NormalizeData.Seurat <- function(
  object,
  assay.use = NULL,
  normalization.method = "LogNormalize",
  scale.factor = 1e4,
  verbose = TRUE,
  ...
) {
  # parameters.to.store <- as.list(environment(), all = TRUE)[names(formals("NormalizeData"))]
  assay.use <- assay.use %||% DefaultAssay(object = object)
  # object <- SetCalcParams(
  #   object = object,
  #   calculation = "NormalizeData",
  #   ... = parameters.to.store
  # )
  assay.data <- GetAssay(object = object, assay.use = assay.use)
  assay.data <- NormalizeData(
    object = assay.data,
    normalization.method = normalization.method,
    scale.factor = scale.factor,
    verbose = verbose,
    ...
  )
  object[[assay.use]] <- assay.data
  object <- LogSeuratCommand(object = object)
  return(object)
}

#' Old R based implementation of ScaleData. Scales and centers the data
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
#' @importFrom utils txtProgressBar setTxtProgressBar
#'
#' @export
#'
#' @examples
#' \dontrun{
#' pbmc_small <- ScaleDataR(object = pbmc_small)
#' }
#'
ScaleDataR <- function(
  object,
  genes.use = NULL,
  data.use = NULL,
  do.scale = TRUE,
  do.center = TRUE,
  scale.max = 10
) {
  genes.use <- SetIfNull(x = genes.use, default = rownames(x = object@data))
  genes.use <- intersect(x = genes.use, y = rownames(x = object@data))
  data.use <- SetIfNull(x = data.use, default = object@data[genes.use, ])
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
    message("Scaling data matrix")
    pb <- txtProgressBar(min = 0, max = max.bin, style = 3, file = stderr())
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
  object@scale.data[is.na(object@scale.data)] <- 0
  return(object)
}

#' Scale and center the data.
#'
#' Scales and centers genes in the dataset. If variables are provided in vars.to.regress,
#' they are individually regressed against each gene, and the resulting residuals are
#' then scaled and centered.
#'
#' ScaleData now incorporates the functionality of the function formerly known
#' as RegressOut (which regressed out given the effects of provided variables
#' and then scaled the residuals). To make use of the regression functionality,
#' simply pass the variables you want to remove to the vars.to.regress parameter.
#'
#' Setting center to TRUE will center the expression for each gene by subtracting
#' the average expression for that gene. Setting scale to TRUE will scale the
#' expression level for each gene by dividing the centered gene expression
#' levels by their standard deviations if center is TRUE and by their root mean
#' square otherwise.
#'
#'
#' @param object Seurat object
#' @param genes.use Vector of gene names to scale/center. Default is all genes
#' in object@@data.
#' @param data.use Can optionally pass a matrix of data to scale, default is
#' object@data[genes.use, ]
#' @param vars.to.regress Variables to regress out (previously latent.vars in
#' RegressOut). For example, nUMI, or percent.mito.
#' @param model.use Use a linear model or generalized linear model
#' (poisson, negative binomial) for the regression. Options are 'linear'
#' (default), 'poisson', and 'negbinom'
#' @param use.umi Regress on UMI count data. Default is FALSE for linear
#' modeling, but automatically set to TRUE if model.use is 'negbinom' or 'poisson'
#' @param do.scale Whether to scale the data.
#' @param do.center Whether to center the data.
#' @param scale.max Max value to return for scaled data. The default is 10.
#' Setting this can help reduce the effects of genes that are only expressed in
#' a very small number of cells. If regressing out latent variables and using a
#' non-linear model, the default is 50.
#' @param block.size Default size for number of genes to scale at in a single
#' computation. Increasing block.size may speed up calculations but at an
#' additional memory cost.
#' @param min.cells.to.block If object contains fewer than this number of cells,
#' don't block for scaling calculations.
#' @param verbose Displays a progress bar for scaling procedure
#' @param assay.type Assay to scale data for. Default is RNA. Can be changed for
#' multimodal analyses.
#' @param do.cpp By default (TRUE), most of the heavy lifting is done in c++.
#' We've maintained support for our previous implementation in R for
#' reproducibility (set this to FALSE) as results can change slightly due to
#' differences in numerical precision which could affect downstream calculations.
#' @param check.for.norm Check to see if data has been normalized, if not,
#' output a warning (TRUE by default)
#' @param do.par use parallel processing for regressing out variables faster.
#' If set to TRUE, will use half of the machines available cores (FALSE by default)
#' @param num.cores If do.par = TRUE, specify the number of cores to use.
#'
#' @return Returns a seurat object with object@@scale.data updated with scaled
#' and/or centered data.
#'
#' @importFrom utils txtProgressBar setTxtProgressBar
#'
#' @export
#'
#' @examples
#' \dontrun{
#' pbmc_small <- ScaleData(object = pbmc_small)
#' # To regress out certain effects
#' pbmc_small = ScaleData(object = pbmc_small, vars.to.regress = effects_list)
#' }
#'
ScaleDataOld <- function(
  object,
  genes.use = NULL,
  data.use = NULL,
  vars.to.regress,
  model.use = 'linear',
  use.umi = FALSE,
  do.scale = TRUE,
  do.center = TRUE,
  scale.max = 10,
  block.size = 1000,
  min.cells.to.block = 3000,
  verbose = TRUE,
  assay.type = "RNA",
  do.cpp = TRUE,
  check.for.norm = TRUE,
  do.par = FALSE,
  num.cores = 1
) {
  data.use <- SetIfNull(
    x = data.use,
    default = GetAssayData(
      object = object,
      assay.type = assay.type,
      slot = "data"
    )
  )
  if (check.for.norm) {
    if (!("NormalizeData" %in% names(object@calc.params))) {
      cat("NormalizeData has not been run, therefore ScaleData is running on non-normalized values. Recommended workflow is to run NormalizeData first.\n")
    }
    if (is.null(object@calc.params$NormalizeData$normalization.method)) {
      cat("ScaleData is running on non-normalized values. Recommended workflow is to run NormalizeData first.\n")
    }
  }
  genes.use <- SetIfNull(x = genes.use, default = rownames(x = data.use))
  genes.use <- as.vector(
    x = intersect(
      x = genes.use,
      y = rownames(x = data.use)
    )
  )
  data.use <- data.use[genes.use, ]
  if (!missing(x = vars.to.regress) && !is.null(x = vars.to.regress)) {
    data.use <- RegressOutResid(
      object = object,
      vars.to.regress = vars.to.regress,
      genes.regress = genes.use,
      use.umi = use.umi,
      model.use = model.use,
      verbose = verbose,
      do.par = do.par,
      num.cores = num.cores
    )
    if (model.use != "linear") {
      use.umi <- TRUE
    }
    if (use.umi && missing(scale.max)) {
      scale.max <- 50
    }
  }
  parameters.to.store <- as.list(environment(), all = TRUE)[names(formals("ScaleData"))]
  parameters.to.store$data.use <- NULL

  object <- SetCalcParams(
    object = object,
    calculation = "ScaleData",
    ... = parameters.to.store
  )
  if (!do.cpp) {
    return(ScaleDataR(
      object = object,
      data.use = data.use,
      do.scale = do.scale,
      do.center = do.center,
      scale.max = scale.max,
      genes.use = genes.use
    ))
  }
  scaled.data <- matrix(
    data = NA,
    nrow = length(x = genes.use),
    ncol = ncol(x = object@data
    )
  )
  rownames(scaled.data) <- genes.use
  if (length(object@cell.names) <= min.cells.to.block) {
    block.size <- length(genes.use)
  }
  gc()
  colnames(scaled.data) <- colnames(object@data)
  max.block <- ceiling(x = length(x = genes.use) / block.size)
  gc()
  if (verbose) {
    message("Scaling data matrix")
    pb <- txtProgressBar(min = 0, max = max.block, style = 3, file = stderr())
  }
  for (i in 1:max.block) {
    my.inds <- ((block.size * (i - 1)):(block.size * i - 1)) + 1
    my.inds <- my.inds[my.inds <= length(x = genes.use)]
    if (class(x = data.use) == "dgCMatrix" | class(x = data.use) == "dgTMatrix") {
      data.scale <- FastSparseRowScale(
        mat = data.use[genes.use[my.inds], , drop = F],
        scale = do.scale,
        center = do.center,
        scale_max = scale.max,
        display_progress = FALSE
      )
    } else {
      data.scale <- FastRowScale(
        mat = as.matrix(x = data.use[genes.use[my.inds], , drop = F]),
        scale = do.scale,
        center = do.center,
        scale_max = scale.max,
        display_progress = FALSE
      )
    }
    dimnames(x = data.scale) <- dimnames(x = data.use[genes.use[my.inds], ])
    scaled.data[genes.use[my.inds], ] <- data.scale
    rm(data.scale)
    gc()
    if (verbose) {
      setTxtProgressBar(pb, i)
    }
  }
  if (verbose) {
    close(pb)
  }
  object <- SetAssayData(
    object = object,
    assay.type = assay.type,
    slot = 'scale.data',
    new.data = scaled.data
  )
  gc()
  object@scale.data[is.na(object@scale.data)] <- 0
  return(object)
}

#' @export
#'
ScaleData.default <- function(
  object,
  features.use = NULL,
  vars.to.regress = NULL,
  latent.data = NULL,
  model.use = 'linear',
  use.umi = FALSE,
  do.scale = TRUE,
  do.center = TRUE,
  scale.max = 10,
  block.size = 1000,
  min.cells.to.block = 3000,
  verbose = TRUE,
  ...
) {
  features.use <- features.use %||% rownames(x = object)
  features.use <- as.vector(x = intersect(x = features.use, y = rownames(x = object)))
  object <- object[features.use, , drop = FALSE]
  scaled.data <- matrix(data = NA, nrow = nrow(x = object), ncol = ncol(x = object))
  dimnames(x = scaled.data) <- dimnames(x = object)
  min.cells.to.block <- min(min.cells.to.block, ncol(x = object))
  gc(verbose = FALSE)
  if (!is.null(x = vars.to.regress)) {
    if (is.null(x = latent.data)) {
      latent.data <- data.frame(row.names = colnames(x = object))
    } else {
      latent.data <- latent.data[colnames(x = object), , drop = FALSE]
      rownames(x = latent.data) <- colnames(x = object)
    }
    if (any(vars.to.regress %in% rownames(x = object))) {
      latent.data <- cbind(
        latent.data,
        t(x = object[vars.to.regress[vars.to.regress %in% rownames(x = object)], ])
      )
    }
    # Currently, RegressOutMatrix will do nothing if latent.data = NULL
    # TODO: implement parallelization for RegressOutMatrix
    if (verbose) {
      message("Regressing out ", paste(vars.to.regress, collapse = ', '))
    }
    object <- RegressOutMatrix(
      data.expr = object,
      latent.data = latent.data,
      features.regress = features.use,
      model.use = model.use,
      use.umi = use.umi,
      verbose = verbose
    )
    gc(verbose = FALSE)
  }
  max.block <- ceiling(x = length(x = features.use) / block.size)
  if (verbose) {
    message("Scaling data matrix")
    pb <- txtProgressBar(min = 0, max = max.block, style = 3, file = stderr())
  }
  for (i in 1:max.block) {
    my.inds <- ((block.size * (i - 1)):(block.size * i - 1)) + 1
    my.inds <- my.inds[my.inds <= length(x = features.use)]
    if (inherits(x = object, what = c('dgCMatrix', 'dgTMatrix'))) {
      scale.function <- FastSparseRowScale
    } else {
      object <- as.matrix(x = object)
      scale.function <- FastRowScale
    }
    data.scale <- scale.function(
      mat = object[features.use[my.inds], , drop = FALSE],
      scale = do.scale,
      center = do.center,
      scale_max = scale.max,
      display_progress = FALSE
    )
    dimnames(x = data.scale) <- dimnames(x = object[features.use[my.inds], ])
    scaled.data[features.use[my.inds], ] <- data.scale
    rm(data.scale)
    gc(verbose = FALSE)
    if (verbose) {
      setTxtProgressBar(pb = pb, value = i)
    }
  }
  if (verbose) {
    close(con = pb)
  }
  scaled.data[is.na(x = scaled.data)] <- 0
  gc(verbose = FALSE)
  return(scaled.data)
}

#' @param latent.data Extra data to regress out, should be cells x latent data
#'
#' @describeIn ScaleData Scale an Assay object
#' @export
#' @method ScaleData Assay
#'
ScaleData.Assay <- function(
  object,
  features.use = NULL,
  vars.to.regress = NULL,
  latent.data = NULL,
  model.use = 'linear',
  use.umi = FALSE,
  do.scale = TRUE,
  do.center = TRUE,
  scale.max = 10,
  block.size = 1000,
  min.cells.to.block = 3000,
  verbose = TRUE,
  ...
) {
  use.umi <- ifelse(test = model.use != 'linear', yes = use.umi, no = use.umi)
  slot.use <- ifelse(test = use.umi, yes = 'matrix', no = 'data')
  features.use <- features.use %||% VariableFeatures(object)
  if (length(features.use) == 0) features.use <- rownames(GetAssayData(object = object, slot = slot.use))
  object <- SetAssayData(
    object = object,
    slot = 'scale.data',
    new.data = ScaleData(
      object = GetAssayData(object = object, slot = slot.use),
      features.use = features.use,
      vars.to.regress = vars.to.regress,
      latent.data = latent.data,
      model.use = model.use,
      use.umi = use.umi,
      do.scale = do.scale,
      do.center = do.center,
      scale.max = scale.max,
      block.size = block.size,
      min.cells.to.block = min.cells.to.block,
      verbose = verbose,
      ...
    )
  )
  return(object)
}

#' @param assay.use Name of Assay to scale
#'
#' @describeIn ScaleData Scale a Seurat object
#' @export
#' @method ScaleData Seurat
#'
ScaleData.Seurat <- function(
  object,
  features.use = NULL,
  assay.use = NULL,
  vars.to.regress = NULL,
  model.use = 'linear',
  use.umi = FALSE,
  do.scale = TRUE,
  do.center = TRUE,
  scale.max = 10,
  block.size = 1000,
  min.cells.to.block = 3000,
  verbose = TRUE,
  ...
) {
  assay.use <- assay.use %||% DefaultAssay(object = object)
  assay.data <- GetAssay(object = object, assay.use = assay.use)
  if (any(vars.to.regress %in% colnames(x = object[]))) {
    latent.data <- object[vars.to.regress[vars.to.regress %in% colnames(x = object[])]]
  } else {
    latent.data <- NULL
  }
  assay.data <- ScaleData(
    object = assay.data,
    features.use = features.use,
    vars.to.regress = vars.to.regress,
    latent.data = latent.data,
    model.use = model.use,
    use.umi = use.umi,
    do.scale = do.scale,
    do.center = do.center,
    scale.max = scale.max,
    block.size = block.size,
    min.cells.to.block = min.cells.to.block,
    verbose = verbose,
    ...
  )
  object[[assay.use]] <- assay.data
  object <- LogSeuratCommand(object = object)
  return(object)
}

#' Normalize raw data
#'
#' Normalize count data per cell and transform to log scale
#'
#' @param data Matrix with the raw count data
#' @param scale.factor Scale the data. Default is 1e4
#' @param verbose Print progress
#'
#' @return Returns a matrix with the normalize and log transformed data
#'
#' @import Matrix
#' @importFrom methods as
#'
#' @export
#'
#' @examples
#' mat <- matrix(data = rbinom(n = 25, size = 5, prob = 0.2), nrow = 5)
#' mat
#' mat_norm <- LogNormalize(data = mat)
#' mat_norm
#'
LogNormalize <- function(data, scale.factor = 1e4, verbose = TRUE) {
  if (class(x = data) == "data.frame") {
    data <- as.matrix(x = data)
  }
  if (class(x = data) != "dgCMatrix") {
    data <- as(object = data, Class = "dgCMatrix")
  }
  # call Rcpp function to normalize
  if (verbose) {
    cat("Performing log-normalization\n", file = stderr())
  }
  norm.data <- LogNorm(data, scale_factor = scale.factor, display_progress = verbose)
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
#' @importFrom methods as
#'
#' @return Matrix with downsampled data
#'
#' @export
#'
#' @examples
#' raw_data = as.matrix(x = pbmc_small@raw.data)
#' downsampled = SampleUMI(data = raw_data)
#' head(x = downsampled)
#'
SampleUMI <- function(
  data,
  max.umi = 1000,
  upsample = FALSE,
  progress.bar = FALSE
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
  new_data = RunUMISamplingPerCell(
    data = data,
    sample_val = max.umi,
    upsample = upsample,
    display_progress = progress.bar
  )
  dimnames(new_data) <- dimnames(data)
  return(new_data)
}

#' @export
#'
FindVariableFeatures.default <- function(
  object,
  mean.function = FastExpMean,
  dispersion.function = FastLogVMR,
  num.bin = 20,
  binning.method = "equal_width",
  verbose = TRUE,
  ...
) {
  if (!inherits(x = mean.function, what = 'function')) {
    stop("'mean.function' must be a function")
  }
  if (!inherits(x = dispersion.function, what = 'function')) {
    stop("'dispersion.function' must be a function")
  }
  if (!inherits(x = object, 'Matrix')) {
    object <- as(object = as.matrix(x = object), Class = 'Matrix')
  }
  if (!inherits(x = object, what = 'dgCMatrix')) {
    object <- as(object = object, Class = 'dgCMatrix')
  }
  feature.mean <- mean.function(object, verbose)
  feature.dispersion <- dispersion.function(object, verbose)
  names(x = feature.mean) <- names(x = feature.dispersion) <- rownames(x = object)
  feature.dispersion[is.na(x = feature.dispersion)] <- 0
  feature.mean[is.na(x = feature.mean)] <- 0
  data.x.breaks <- switch(
    EXPR = binning.method,
    'equal_width' = num.bin,
    'equal_frequency' = c(
      -1,
      quantile(
        x = feature.mean[feature.mean > 0],
        probs = seq.int(from = 0, to = 1, length.out = num.bin)
      )
    ),
    stop("Unknown binning method: ", binning.method)
  )
  data.x.bin <- cut(x = feature.mean, breaks = data.x.breaks)
  names(x = data.x.bin) <- names(x = feature.mean)
  mean.y <- tapply(X = feature.dispersion, INDEX = data.x.bin, FUN = mean)
  sd.y <- tapply(X = feature.dispersion, INDEX = data.x.bin, FUN = sd)
  feature.dispersion.scaled <- (feature.dispersion - mean.y[as.numeric(x = data.x.bin)]) /
    sd.y[as.numeric(x = data.x.bin)]
  names(x = feature.dispersion.scaled) <- names(x = feature.mean)
  hvf.info <- data.frame(feature.mean, feature.dispersion, feature.dispersion.scaled)
  rownames(x = hvf.info) <- rownames(x = object)
  colnames(x = hvf.info) <- c('mean', 'dispersion', 'dispersion.scaled')
  return(hvf.info)
}

#' @param num.features Number of features to select as top variable features;
#' only used when \code{selection.method = 'dispersion'}
#' @param mean.cutoff A two-length numeric vector with low- and high-cutoffs for
#' feature means
#' @param dispersion.cutoff A two-length numeric vector with low- and high-cutoffs for
#' feature dispersions
#' @param selection.method How to choose top variable features, either 'mean.var.plot' or 'dispersion
#'
#' @describeIn FindVariableFeatures Find variable features in an Assay object
#' @export
#' @method FindVariableFeatures Assay
#'
FindVariableFeatures.Assay <- function(
  object,
  mean.function = FastExpMean,
  dispersion.function = FastLogVMR,
  num.bin = 20,
  binning.method = "equal_width",
  num.features = 1000,
  mean.cutoff = c(0.1, 8),
  dispersion.cutoff = c(1, Inf),
  selection.method = "mean.var.plot",
  verbose = TRUE,
  ...
) {
  if (length(x = mean.cutoff) != 2 || length(x = dispersion.cutoff) != 2) {
    stop("Both 'mean.cutoff' and 'dispersion.cutoff' must be two numbers")
  }
  hvf.info <- FindVariableFeatures(
    object = GetAssayData(object = object),
    mean.function = mean.function,
    dispersion.function = dispersion.function,
    num.bin = num.bin,
    binning.method = binning.method,
    verbose = verbose,
    ...
  )
  object[[names(x = hvf.info)]] <- list(hvf.info)
  hvf.info <- hvf.info[order(hvf.info$dispersion, decreasing = TRUE), , drop = FALSE]
  top.features <- switch(
    EXPR = selection.method,
    'mean.var.plot' = {
      means.use <- (hvf.info[, 1] > mean.cutoff[1]) & (hvf.info[, 1] < mean.cutoff[2])
      dispersions.use <- (hvf.info[, 3] > dispersion.cutoff[1]) & (hvf.info[, 3] < dispersion.cutoff[2])
      rownames(x = hvf.info)[which(x = means.use & dispersions.use)]
    },
    'dispersion' = head(x = rownames(x = hvf.info), n = num.features),
    stop("Unkown selection method: ", selection.method)
  )
  VariableFeatures(object = object) <- top.features
  return(object)
}

#' @inheritParams FindVariableFeatures.Assay
#' @param assay.use ...
#'
#' @describeIn FindVariableFeatures Find variable features in a Seurat object
#' @export
#' @method FindVariableFeatures Seurat
#'
FindVariableFeatures.Seurat <- function(
  object,
  assay.use = NULL,
  mean.function = FastExpMean,
  dispersion.function = FastLogVMR,
  num.bin = 20,
  binning.method = "equal_width",
  num.features = 1000,
  mean.cutoff = c(0.1, 8),
  dispersion.cutoff = c(1, Inf),
  selection.method = "dispersion",
  verbose = TRUE,
  ...
) {
  assay.use <- assay.use %||% DefaultAssay(object = object)
  assay.data <- GetAssay(object = object, assay.use = assay.use)
  assay.data <- FindVariableFeatures(
    object = assay.data,
    mean.function = mean.function,
    dispersion.function = dispersion.function,
    num.bin = num.bin,
    binning.method = binning.method,
    mean.cutoff = mean.cutoff,
    dispersion.cutoff = dispersion.cutoff,
    verbose = verbose,
    selection.method = selection.method,
    ...
  )
  object[[assay.use]] <- assay.data
  object <- LogSeuratCommand(object = object)
  return(object)
}

#' Identify variable genes
#'
#' Identifies genes that are outliers on a 'mean variability plot'. First, uses
#' a function to calculate average expression (mean.function) and dispersion (dispersion.function)
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
#' @param mean.function Function to compute x-axis value (average expression). Default
#' is to take the mean of the detected (i.e. non-zero) values
#' @param dispersion.function Function to compute y-axis value (dispersion). Default is to
#' take the standard deviation of all values/
#' @param do.plot Plot the average/dispersion relationship
#' @param set.var.genes Set object@@var.genes to the identified variable genes
#' (default is TRUE)
#' @param x.low.cutoff Bottom cutoff on x-axis for identifying variable genes
#' @param x.high.cutoff Top cutoff on x-axis for identifying variable genes
#' @param y.cutoff Bottom cutoff on y-axis for identifying variable genes
#' @param y.high.cutoff Top cutoff on y-axis for identifying variable genes
#' @param num.bin Total number of bins to use in the scaled analysis (default
#' is 20)
#' @param binning.method Specifies how the bins should be computed. Available methods are:
#' \itemize{
#' \item{equal_width:}{ each bin is of equal width along the x-axis [default]}
#' \item{equal_frequency:}{ each bin contains an equal number of genes (can increase
#' statistical power to detect overdispersed genes at high expression values, at
#' the cost of reduced resolution along the x-axis)}
#' }
#' @param selection.method Specifies how to select the genes to store in @@var.genes.
#' \itemize{
#' \item{mean.var.plot: }{Default method, placing cutoffs on the mean variablility plot}
#' \item{dispersion: }{Choose the top.genes with the highest dispersion}
#' }
#' @param top.genes Selects the genes with the highest value according to the
#' selection method.
#' @param do.recalc TRUE by default. If FALSE, plots and selects variable genes without recalculating statistics for each gene.
#' @param sort.results If TRUE (by default), sort results in object@hvg.info in decreasing order of dispersion
#' @param do.cpp Run c++ version of mean.function and dispersion.function if they
#' exist.
#' @param verbose show progress bar for calculations
#' @param ... Extra parameters to VariableGenePlot
#' @inheritParams VariableGenePlot
#'
#' @importFrom MASS kde2d
#' @importFrom methods as
#' @importFrom utils txtProgressBar setTxtProgressBar
#'
#' @return Returns a Seurat object, placing variable genes in object@@var.genes.
#' The result of all analysis is stored in object@@hvg.info
#'
#' @seealso \code{VariableGenePlot}
#'
#' @export
#'
#' @examples
#' pbmc_small <- FindVariableGenes(object = pbmc_small, do.plot = FALSE)
#' pbmc_small@var.genes
#'
FindVariableGenes <- function(
  object,
  mean.function = ExpMean,
  dispersion.function = LogVMR,
  do.plot = TRUE,
  set.var.genes = TRUE,
  x.low.cutoff = 0.1,
  x.high.cutoff = 8,
  y.cutoff = 1,
  y.high.cutoff = Inf,
  num.bin = 20,
  binning.method = "equal_width",
  selection.method = "mean.var.plot",
  top.genes = 1000,
  do.recalc = TRUE,
  sort.results = TRUE,
  do.cpp = TRUE,
  verbose = TRUE,
  ...
) {
  parameters.to.store <- as.list(environment(), all = TRUE)[names(formals("FindVariableGenes"))]
  parameters.to.store$mean.function <- as.character(substitute(mean.function))
  parameters.to.store$dispersion.function <- as.character(substitute(dispersion.function))
  object <- SetCalcParams(
    object = object,
    calculation = "FindVariableGenes",
    ... = parameters.to.store
  )
  data <- object@data
  genes.use <- rownames(x = object@data)
  if (do.recalc) {
    if (do.cpp) {
      if (!identical(mean.function, ExpMean)) {
        warning("No equivalent mean.function implemented in c++ yet, falling back to R version")
        do.cpp <- FALSE
      }
      if (!identical(dispersion.function, LogVMR)) {
        warning("No equivalent dispersion.function implemented in c++ yet, falling back to R version")
        do.cpp <- FALSE
      }
    }
    if (do.cpp) {
      if (class(x = data) != "dgCMatrix") {
        data <- as(as.matrix(data), "dgCMatrix")
      }
      gene.mean <- FastExpMean(data, verbose)
      names(gene.mean) <- genes.use
      gene.dispersion <- FastLogVMR(data, verbose)
      names(gene.dispersion) <- genes.use
    }
    if (!do.cpp) {
      gene.mean <- rep(x = 0, length(x = genes.use))
      names(x = gene.mean) <- genes.use
      gene.dispersion <- gene.mean
      gene.dispersion.scaled <- gene.mean
      bin.size <- 1000
      max.bin <- floor(x = length(x = genes.use) / bin.size) + 1
      if (verbose) {
        message("Calculating gene dispersion")
        pb <- txtProgressBar(min = 0, max = max.bin, style = 3, file = stderr())
      }
      for (i in 1:max.bin) {
        my.inds <- ((bin.size * (i - 1)):(bin.size * i - 1)) + 1
        my.inds <- my.inds[my.inds <= length(x = genes.use)]
        genes.iter <- genes.use[my.inds]
        data.iter <- data[genes.iter, , drop = F]
        gene.mean[genes.iter] <- apply(X = data.iter, MARGIN = 1, FUN = mean.function)
        gene.dispersion[genes.iter] <- apply(X = data.iter, MARGIN = 1, FUN = dispersion.function)
        if (verbose) {
          setTxtProgressBar(pb = pb, value = i)
        }
      }
      if (verbose) {
        close(con = pb)
      }
    }
    gene.dispersion[is.na(x = gene.dispersion)] <- 0
    gene.mean[is.na(x = gene.mean)] <- 0
    if (binning.method == "equal_width") {
      data_x_bin <- cut(x = gene.mean, breaks = num.bin)
    }
    else if (binning.method == "equal_frequency") {
      data_x_bin <- cut(x = gene.mean, breaks = c(-1, quantile(gene.mean[gene.mean > 0], probs = seq(from = 0, to = 1,length.out = num.bin))))
    }
    else {
      stop(paste0("Invalid selection: '",binning.method,"' for 'binning.method'."))
    }
    names(x = data_x_bin) <- names(x = gene.mean)
    mean_y <- tapply(X = gene.dispersion, INDEX = data_x_bin, FUN = mean)
    sd_y <- tapply(X = gene.dispersion, INDEX = data_x_bin, FUN = sd)
    gene.dispersion.scaled <- (gene.dispersion - mean_y[as.numeric(x = data_x_bin)]) /
      sd_y[as.numeric(x = data_x_bin)]
    gene.dispersion.scaled[is.na(x = gene.dispersion.scaled)] <- 0
    names(x = gene.dispersion.scaled) <- names(x = gene.mean)
    mv.df <- data.frame(gene.mean, gene.dispersion, gene.dispersion.scaled)
    rownames(x = mv.df) <- rownames(x = data)
    object@hvg.info <- mv.df
  }
  gene.mean <- object@hvg.info[, 1]
  gene.dispersion <- object@hvg.info[, 2]
  gene.dispersion.scaled <- object@hvg.info[, 3]
  names(x = gene.mean) <- names(x = gene.dispersion) <- names(x = gene.dispersion.scaled) <- rownames(x = object@hvg.info)
  pass.cutoff <- names(x = gene.mean)[which(
    x = (
      (gene.mean > x.low.cutoff) & (gene.mean < x.high.cutoff)
    ) &
      (gene.dispersion.scaled > y.cutoff) &
      (gene.dispersion.scaled < y.high.cutoff)
  )]
  if (do.plot) {
    VariableGenePlot(
      object = object,
      x.low.cutoff = x.low.cutoff,
      x.high.cutoff = x.high.cutoff,
      y.cutoff = y.cutoff,
      y.high.cutoff = y.high.cutoff,
      ...
    )
  }
  if (set.var.genes) {
    if(selection.method == "dispersion") {
      sort.results <- TRUE
    }
    if (sort.results) {
      object@hvg.info <- object@hvg.info[order(
        object@hvg.info$gene.dispersion,
        decreasing = TRUE
      ),]
    }
    if(selection.method == "dispersion"){
      object@var.genes <- head(rownames(object@hvg.info), top.genes)
    } else if(selection.method == "mean.var.plot") {
      object@var.genes <- pass.cutoff
    } else {
      stop("Invalid selection.method")
    }
    return(object)
  } else {
    return(pass.cutoff)
  }
}

#' Return a subset of the Seurat object
#'
#' Creates a Seurat object containing only a subset of the cells in the
#' original object. Takes either a list of cells to use as a subset, or a
#' parameter (for example, a gene), to subset on.
#'
#' @param object Seurat object
#' @param subset.names Parameters to subset on. Eg, the name of a gene, PC1, a
#' column name in object@@meta.data, etc. Any argument that can be retreived
#' using FetchData
#' @param low.thresholds Low cutoffs for the parameters (default is -Inf)
#' @param high.thresholds High cutoffs for the parameters (default is Inf)
#' @param cells.use A vector of cell names to use as a subset
#'
#' @return Returns a Seurat object containing only the relevant subset of cells
#'
#' @export
#'
#' @examples
#' head(x = FetchData(object = pbmc_small, vars.all = 'LTB'))
#' pbmc_filtered <- FilterCells(
#'   object = pbmc_small,
#'   subset.names = 'LTB',
#'   high.thresholds = 6
#' )
#' head(x = FetchData(object = pbmc_filtered, vars.all = 'LTB'))
#'
FilterCells <- function(
  object,
  subset.names,
  low.thresholds,
  high.thresholds,
  cells.use = NULL
) {
  if (missing(x = low.thresholds)) {
    low.thresholds <- replicate(n = length(x = subset.names), expr = -Inf)
  }
  if (missing(x = high.thresholds)) {
    high.thresholds <- replicate(n = length(x = subset.names), expr = Inf)
  }
  length.check <- sapply(
    X = list(subset.names, low.thresholds, high.thresholds),
    FUN = length
  )
  if (length(x = unique(x = length.check)) != 1) {
    stop("'subset.names', 'low.thresholds', and 'high.thresholds' must all have the same length")
  }
  data.subsets <- data.frame(subset.names, low.thresholds, high.thresholds)
  cells.use <- cells.use %||% colnames(x = object)
  for (i in seq(nrow(data.subsets))) {
    cells.use <- WhichCells(
      object = object,
      cells.use = cells.use,
      subset.name = data.subsets[i, 1],
      low.threshold = data.subsets[i, 2],
      high.threshold = data.subsets[i, 3]
    )
  }
  object <- SubsetData(object, cells.use = cells.use)
  object <- LogSeuratCommand(object)
  return(object)
}
