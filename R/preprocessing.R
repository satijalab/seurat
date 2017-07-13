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

  
  nGene <- data.ngene
  nUMI <- data.nmol
  object@data.info <- data.frame(nGene, nUMI)
  if (! is.null(x = meta.data)) {
    object <- AddMetaData(object = object, metadata = meta.data)
  }
  object@data.info[names(object@ident),"orig.ident"] <- object@ident
  object@scale.data <- matrix()
  object@assay <- list()
  if(do.scale | do.center) {
    object <- ScaleData(object = object, do.scale = do.scale, do.center = do.center)
  }

  spatial.obj <- new(
    Class = "spatial.info"
  )
  object@spatial <- spatial.obj
  object@spatial@mix.probs <- data.frame(data.ngene)
  colnames(x = object@spatial@mix.probs)[1] <- "nGene"
  
  
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
#' @export
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
  genes.use <- ainb(genes.use,rownames(object@data))
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


#' Scale and center the data.
#'
#' Scales and centers the data. If latent variables are provided (latent.vars), their effects are 
#' removed through regression and the resulting residuals are then scaled and centered.
#'
#' 
#' @param object Seurat object
#' @param genes.use Vector of gene names to scale/center. Default is all genes in object@@data.
#' @param data.use Can optionally pass a matrix of data to scale, default is object@data[genes.use,]
#' @param latent.vars effects to regress out
#' @param genes.regress gene to run regression for (default is all genes)
#' @param model.use Use a linear model or generalized linear model (poisson, negative binomial) for the regression. Options are 'linear' (default), 'poisson', and 'negbinom'
#' @param use.umi Regress on UMI count data. Default is FALSE for linear modeling, but automatically set to TRUE if model.use is 'negbinom' or 'poisson'
#' @param do.scale Whether to scale the data.
#' @param do.center Whether to center the data.
#' @param scale.max Max value to return for scaled data. The default is 10. Setting this can help
#' reduce the effects of genes that are only expressed in a very small number of cells.
#' @param assay.type Assay to scale data for. Default is RNA. Can be changed for multimodal analyses.
#' @param do.cpp By default (TRUE), most of the heavy lifting is done in c++. We've maintained 
#' support for our previous implementation in R for reproducibility (set this to FALSE) as results 
#' can change slightly due to differences in numerical precision which could affect downstream 
#' calculations.  
#'
#' @return Returns a seurat object with object@@scale.data updated with scaled and/or centered data.
#'
#' @export
#'
ScaleData <- function(
  object,
  genes.use = NULL,
  data.use = NULL,
  latent.vars,
  model.use = 'linear',
  use.umi = FALSE,
  do.scale = TRUE,
  do.center = TRUE,
  scale.max = 10,
  block.size = 1000,
  display.progress = TRUE,
  assay.type = "RNA",
  do.cpp = TRUE
) {
  data.use <- SetIfNull(x = data.use, default = GetAssayData(object = object, assay.type,"data"))
  genes.use <- SetIfNull(x = genes.use, default = rownames(x = data.use))
  genes.use <- as.vector(
    x = intersect(
      x = genes.use,
      y = rownames(x = data.use)
    )
  )
  data.use <- data.use[genes.use, ]
  
  if(!missing(latent.vars)){
    data.use <- RegressOut(object = object, 
                           latent.vars = latent.vars, 
                           use.umi = use.umi,
                           model.use = model.use)
  }
  if(!do.cpp){
    return(ScaleDataR(object = object, 
                     data.use = data.use, 
                     do.scale = do.scale,
                     do.center = do.center,
                     scale.max = scale.max))
  }
  scaled.data <- matrix(data = NA,
                        nrow = length(x = genes.use),
                        ncol = ncol(x = object@data)
  )
  rownames(scaled.data) <- genes.use
  gc()
  colnames(scaled.data) <- colnames(object@data)
  max.block <- ceiling(length(genes.use)/block.size)
  gc()
  print("Scaling data matrix")
  pb <- txtProgressBar(min = 0, max = max.block, style = 3)
  for (i in 1:max.block) {
    my.inds <- ((block.size * (i - 1)):(block.size * i - 1)) + 1
    my.inds <- my.inds[my.inds <= length(x = genes.use)]
    if (class(data.use) == "dgCMatrix" | class(data.use) == "dgTMatrix") {
      data.scale <- FastSparseRowScale(
        mat = data.use[genes.use[my.inds], , drop = F],
        scale = do.scale,
        center = do.center,
        scale_max = scale.max,
        display_progress = FALSE
      )
    } else {
      data.scale <- FastRowScale(
        mat = as.matrix(data.use[genes.use[my.inds], , drop = F]),
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
    setTxtProgressBar(pb, i)
  }
  close(pb)
  object <- SetAssayData(
    object = object,
    assay.type = assay.type,
    slot = 'scale.data',
    new.data = scaled.data
  )
  gc()
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
#' @param x.low.cutoff Bottom cutoff on x-axis for identifying variable genes
#' @param x.high.cutoff Top cutoff on x-axis for identifying variable genes
#' @param y.cutoff Bottom cutoff on y-axis for identifying variable genes
#' @param y.high.cutoff Top cutoff on y-axis for identifying variable genes
#' @param num.bin Total number of bins to use in the scaled analysis (default
#' is 20)
#' @param do.recalc TRUE by default. If FALSE, plots and selects variable genes without recalculating statistics for each gene.
#' @param sort.results If TRUE (by default), sort results in object@mean.var in decreasing order of dispersion
#' @param ... Extra parameters to VariableGenePlot
#' @inheritParams VariableGenePlot
#'
#' @importFrom MASS kde2d
#'
#' @return Returns a Seurat object, placing variable genes in object@@var.genes.
#' The result of all analysis is stored in object@@mean.var
#'
#' @seealso \code{\link{VariableGenePlot}}
#'
#' @export
#'
FindVariableGenes <- function(
  object,
  fxn.x = expMean,
  fxn.y = logVarDivMean,
  do.plot = TRUE,
  set.var.genes = TRUE,
  x.low.cutoff = 0.1,
  x.high.cutoff = 8,
  y.cutoff = 1,
  y.high.cutoff = Inf,
  num.bin = 20,
  do.recalc = TRUE,
  sort.results = TRUE,
  ...
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
  if (do.plot) {
    VariableGenePlot(
      object = object,
      ...
    )
  }
  # if (do.plot) {
  #   if (plot.both) {
  #     par(mfrow = c(1, 2))
  #     smoothScatter(
  #       x = data.x,
  #       y = data.y,
  #       pch = pch.use,
  #       cex = cex.use,
  #       col = col.use,
  #       xlab = "Average expression",
  #       ylab = "Dispersion",
  #       nrpoints = Inf
  #     )
  #     if (do.contour) {
  #       data.kde <- kde2d(x = data.x, y = data.y)
  #       contour(
  #         x = data.kde,
  #         add = TRUE,
  #         lwd = contour.lwd,
  #         col = contour.col,
  #         lty = contour.lty
  #       )
  #     }
  #     if (do.spike) {
  #       points(
  #         x = data.x[spike.genes],
  #         y = data.y[spike.genes],
  #         pch = 16,
  #         cex = cex.use,
  #         col = spike.col.use
  #       )
  #     }
  #     if (do.text) {
  #       text(
  #         x = data.x[pass.cutoff],
  #         y = data.y[pass.cutoff],
  #         labels = pass.cutoff,
  #         cex = cex.text.use
  #       )
  #     }
  #   }
  #   smoothScatter(
  #     x = data.x,
  #     y = data.norm.y,
  #     pch = pch.use,
  #     cex = cex.use,
  #     col = col.use,
  #     xlab = "Average expression",
  #     ylab = "Dispersion",
  #     nrpoints = Inf
  #   )
  #   if (do.contour) {
  #     data.kde <- kde2d(x = data.x, y = data.norm.y)
  #     contour(
  #       x = data.kde,
  #       add = TRUE,
  #       lwd = contour.lwd,
  #       col = contour.col,
  #       lty = contour.lty
  #     )
  #   }
  #   if (do.spike) {
  #     points(
  #       x = data.x[spike.genes],
  #       y = data.norm.y[spike.genes],
  #       pch = 16,
  #       cex = cex.use,
  #       col = spike.col.use,
  #       nrpoints = Inf
  #     )
  #   }
  #   if (do.text) {
  #     text(
  #       x = data.x[pass.cutoff],
  #       y = data.norm.y[pass.cutoff],
  #       labels = pass.cutoff,
  #       cex = cex.text.use
  #     )
  #   }
  # }
  if (set.var.genes) {
    object@var.genes <- pass.cutoff
    if (sort.results) {
      object@mean.var <- object@mean.var[order(
        object@mean.var$data.y,
        decreasing = TRUE
      ),]
    }
    return(object)
  } else {
    return(pass.cutoff)
  }
}
