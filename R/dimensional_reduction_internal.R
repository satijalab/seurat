#' @include seurat.R
NULL
# Set up dim.reduction class

dim.reduction <- setClass(
  Class = "dim.reduction",
  slots = list(
    cell.embeddings = "matrix",
    gene.loadings = "matrix",
    gene.loadings.full = "matrix",
    sdev = "numeric",
    key = "character",
    jackstraw="ANY",
    misc = "ANY"
  )
)

# Prep data for dimensional reduction
#
# Common checks and preparatory steps before running certain dimensional
# reduction techniques
#
# @param object        Seurat object
# @param genes.use     Genes to use as input for the dimensional reduction technique.
#                      Default is object@@var.genes
# @param dims.compute  Number of dimensions to compute
# @param use.imputed   Whether to run the dimensional reduction on imputed values
# @param assay.type Assay to scale data for. Default is RNA. Can be changed for multimodal analysis

PrepDR <- function(
  object,
  genes.use = NULL,
  use.imputed = FALSE,
  assay.type="RNA"
) {

  if (length(object@var.genes) == 0 && is.null(x = genes.use)) {
    stop("Variable genes haven't been set. Run MeanVarPlot() or provide a vector
          of genes names in genes.use and retry.")
  }
  if (use.imputed) {
    data.use <- t(x = scale(x = t(x = object@imputed)))
  } else {
    data.use <- GetAssayData(object, assay.type = assay.type,slot = "scale.data")
  }
  genes.use <- SetIfNull(x = genes.use, default = object@var.genes)
  genes.use <- unique(x = genes.use[genes.use %in% rownames(x = data.use)])
  genes.var <- apply(X = data.use[genes.use, ], MARGIN = 1, FUN = var)
  genes.use <- genes.use[genes.var > 0]
  genes.use <- genes.use[! is.na(x = genes.use)]
  data.use <- data.use[genes.use, ]
  return(data.use)
}

# Get the top genes associated with given dimensional reduction scores
#
# @param i            Dimension for which to pull genes
# @param dim.scores   Matrix containing the dimensional reduction scores
# @param do.balanced  Whether to pull genes associated with both large and small
#                     scores (+/-)
# @param num.genes    Number of genes to return

GetTopGenes <- function(
  i,
  dim.scores,
  do.balanced = FALSE,
  num.genes = 30
) {
  if (do.balanced) {
    num.genes <- round(x = num.genes / 2)
    sx <- dim.scores[order(dim.scores[, i]), , drop = FALSE]
    genes.1 <- (rownames(x = sx[1:num.genes, , drop = FALSE]))
    genes.2 <- (rownames(x = sx[(nrow(x = sx) - num.genes + 1):nrow(x = sx), , drop = FALSE]))
    return(c(genes.1, genes.2))
  } else {
    sx <- dim.scores[rev(x = order(abs(x = dim.scores[, i]))), ,drop = FALSE]
    genes.1 <- (rownames(x = sx[1:num.genes, , drop = FALSE]))
    genes.1 <- genes.1[order(dim.scores[genes.1, i])]
    return(genes.1)
  }
}

# Check group exists either as an ident or that all cells passed as vector are
# present
#
# @param object    Seurat object
# @param group     Identity or vector of cell names
# @param group.id  Corresponds to the the either group1 or group2 parameter from
#                  RunCCA

CheckGroup <- function(object, group, group.id) {
  if (all(group %in% unique(x = object@ident))) {
    cells.use <- WhichCells(object = object, ident = group)
  } else {
    if (all(group %in% object@cell.names)) {
      cells.use <- group
    } else {
      stop(paste(
        group.id,
        "must be either a vector of valid cell names or idents"
      ))
    }
  }
  return(cells.use)
}

# Check that genes have non-zero variance
#
# @param data.use   Gene expression matrix (genes are rows)
# @param genes.use  Genes in expression matrix to check
#
# @return           Returns a vector of genes that is the subset of genes.use
#                   that have non-zero variance
#
CheckGenes <- function(data.use, genes.use) {
  genes.var <- apply(X = data.use[genes.use, ], MARGIN = 1, FUN = var)
  genes.use <- genes.use[genes.var > 0]
  genes.use <- genes.use[! is.na(x = genes.use)]
  return(genes.use)
}

# Run the diagonal canonical correlation procedure
#
# @param mat1         First matrix
# @param mat2         Second matrix
# @param standardize  Standardize matrices - scales columns to have unit
#                     variance and mean 0
# @param k            Number of canonical correlation vectors (CCs) to calculate
#
# @return             Returns the canonical correlation vectors - corresponding
#                     to the left and right singular vectors after SVD - as well
#                     as the singular values.
#
CanonCor <- function(mat1, mat2, standardize = TRUE, k = 20) {
  set.seed(seed = 42)
  if (standardize) {
    mat1 <- Standardize(mat = mat1, display_progress = FALSE)
    mat2 <- Standardize(mat = mat2, display_progress = FALSE)
  }
  mat3 <- FastMatMult(m1 = t(x = mat1), m2 = mat2)
  cca.svd <- irlba(A = mat3, nv = k)
  return(list(u = cca.svd$u, v = cca.svd$v, d = cca.svd$d))
}

# Calculate percent variance explained
#
# Projects dataset onto the orthonormal space defined by some dimensional
# reduction technique (e.g. PCA, CCA) and calculates the percent of the
# variance in gene expression explained by each cell in that lower dimensional
# space.
#
# @param object          Seurat object
# @param reduction.type  Name of the reduction to use for the projection
# @param dims.use        Vector of dimensions to project onto (default is the
#                        1:number stored for given technique)
# @param genes.use       vector of genes to use in calculation
#
# @return                Returns a Seurat object wih the variance in gene
#                        expression explained by each cell in a low dimensional
#                        space stored as metadata.
#
CalcProjectedVar <- function(
  object,
  low.dim.data,
  reduction.type = "pca",
  dims.use,
  genes.use
) {
  if (missing(x = low.dim.data)) {
    low.dim.data <- CalcLDProj(
      object = object,
      reduction.type = reduction.type,
      dims.use = dims.use,
      genes.use = genes.use
    )
  }
  projected.var <- apply(X = low.dim.data, MARGIN = 2, FUN = var)
  calc.name <- paste0(reduction.type, ".var")
  object <- AddMetaData(
    object = object,
    metadata = projected.var,
    col.name = calc.name
  )
  return(object)
}

# Calculate a low dimensional projection of the data. First forms an orthonormal
# basis of the gene loadings via QR decomposition, projects the data onto that
# basis, and reconstructs the data using on the dimensions specified.
#
# @param object          Seurat object
# @param reduction.type  Type of dimensional reduction to use
# @param dims.use        Dimensions to use in calculation
# @param genes.use       Genes to consider when calculating
#
# @return                Returns a matrix with the low dimensional reconstruction
#
CalcLDProj <- function(object, reduction.type, dims.use, genes.use) {
  if (missing(x = dims.use)){
    dims.use <- 1:ncol(x = GetCellEmbeddings(
      object = object,
      reduction.type = reduction.type
    ))
  }
  x.vec <- GetGeneLoadings(
    object = object,
    reduction.type = reduction.type,
    dims.use = dims.use,
    genes.use = genes.use
  )
  # form orthonormal basis via QR
  x.norm <- qr.Q(qr = qr(x = x.vec))
  data.use <- object@scale.data[rownames(x.vec), ]
  # project data onto othronormal basis
  projected.data <- t(x = data.use) %*% x.norm
  # reconstruct data using only dims specified
  low.dim.data <- x.norm %*% t(x = projected.data)
  return(low.dim.data)
}

# MultiCCA helper function - calculates critical value (when to stop iterating
# in the while loop)
#
# Modified from PMA package
# @references Witten, Tibshirani, and Hastie, Biostatistics 2009
# @references \url{https://github.com/cran/PMA/blob/master/R/MultiCCA.R}
#
# @param mat.list list of matrices
# @param ws vector of projection vectors
# @param num.sets number of datasets
#
# @return returns updated critical value
#
GetCrit <- function(mat.list, ws, num.sets){
  crit <- 0
  for(i in 2:num.sets){
    for(j in 1:(i-1)){
      crit <- crit + t(ws[[i]])%*%t(mat.list[[i]])%*%mat.list[[j]]%*%ws[[j]]
    }
  }
  return(crit)
}

# MultiCCA helper function - updates W
#
# Modified from PMA package
# @references Witten, Tibshirani, and Hastie, Biostatistics 2009
# @references \url{https://github.com/cran/PMA/blob/master/R/MultiCCA.R}
#
# @param mat.list list of matrices
# @param i index of current matrix
# @param num.sets number of datasets
# @param ws initial vector of projection vectors
# @param ws.final final vector of projection vectors
#
# @return returns updated w value
#
UpdateW <- function(mat.list, i, num.sets, ws, ws.final){
  tots <- 0
  for(j in (1:num.sets)[-i]){
    diagmat <- (t(ws.final[[i]])%*%t(mat.list[[i]]))%*%(mat.list[[j]]%*%ws.final[[j]])
    diagmat[row(diagmat)!=col(diagmat)] <- 0
    tots <- tots + t(mat.list[[i]])%*%(mat.list[[j]]%*%ws[[j]]) - ws.final[[i]]%*%(diagmat%*%(t(ws.final[[j]])%*%ws[[j]]))
  }
  w <- tots/l2n(tots)
  return(w)
}

# Calculates the l2-norm of a vector
#
# Modified from PMA package
# @references Witten, Tibshirani, and Hastie, Biostatistics 2009
# @references \url{https://github.com/cran/PMA/blob/master/R/PMD.R}
#
# @param vec numeric vector
#
# @return returns the l2-norm.
#
l2n <- function(vec){
  a <- sqrt(sum(vec^2))
  if(a==0){
    a <- .05
  }
  return(a)
}

# MultiCCA helper function - calculates correlation
#
# Modified from PMA package
# @references Witten, Tibshirani, and Hastie, Biostatistics 2009
# @references \url{https://github.com/cran/PMA/blob/master/R/MultiCCA.R}
#
# @param mat.list list of matrices to calculate correlation
# @param ws vector of projection vectors
# @param num.sets number of datasets
#
# @return total correlation
#
GetCors <- function(mat.list, ws, num.sets){
  cors <- 0
  for(i in 2:num.sets){
    for(j in 1:(i-1)){
      thiscor  <-  cor(mat.list[[i]]%*%ws[[i]], mat.list[[j]]%*%ws[[j]])
      if(is.na(thiscor)) thiscor <- 0
      cors <- cors + thiscor
    }
  }
  return(cors)
}
