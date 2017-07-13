#' Shuffle a vector
#' @param x A vector
#' @return A vector with the same values of x, just in random order
#' @export
#'
Shuffle <- function(x) {
  return(x[base::sample.int(
    n = base::length(x = x),
    size = base::length(x = x),
    replace = FALSE
  )])
}

#' Remove data from a table
#'
#' This function will remove any rows from a data frame or matrix
#' that contain certain values
#'
#' @param to.remove A vector of values that indicate removal
#' @param data A data frame or matrix
#'
#' @return A data frame or matrix with values removed by row
#'
#' @export
#'
RemoveFromTable <- function(to.remove, data) {
  remove.indecies <- apply(
    X = data,
    MARGIN = 2,
    FUN = function(col) {
      return(which(x = col %in% to.remove))
    }
  )
  remove.indecies <- unlist(x = remove.indecies)
  remove.indecies <- as.numeric(x = remove.indecies)
  if (length(x = remove.indecies) == 0) {
    return(data)
  } else {
    return(data[-remove.indecies, ])
  }
}

#' Make object sparse
#'
#' Converts stored data matrices to sparse matrices to save space. Converts
#' object@@raw.data and object@@data to sparse matrices.
#' @param object Seurat object
#' @return Returns a seurat object with data converted to sparse matrices.
#' @import Matrix
#' @export
#'
MakeSparse <- function(object) {
  if (class(object@raw.data) == "data.frame") {
    object@raw.data <- as.matrix(x = object@raw.data)
  }
  if (class(object@data) == "data.frame") {
    object@data <- as.matrix(x = object@data)
  }
  object@raw.data <- as(object = object@raw.data, Class = "dgCMatrix")
  object@data <- as(object = object@data, Class = "dgCMatrix")
  return(object)
}

#' Convert old Seurat object to accomodate new features
#'
#' Adds the object@@dr slot to older objects and moves the stored PCA/ICA
#' analyses to new slot. Moves snn to new slot.
#'
#' @param object Seurat object
#'
#' @return Returns a Seurat object compatible with latest changes
#'
#' @export
#'
ConvertSeurat <- function(object) {
  if (.hasSlot(object, "version")) {
    if(packageVersion("Seurat") >= package_version("2.0.0")){
      return(object)
    }
  }

  new.object <- new("seurat",
                    raw.data = object@raw.data)
  new.slots <- slotNames(new.object)
  for(s in new.slots){
    new.object <- FillSlot(slot.name = s, old.object = object,
                              new.object = new.object)
  }
  # Conversion from development versions prior to 2.0.0
  # Slots to replace: pca.x, pca.rot, pca.x.full, tsne.rot, ica.rot, ica.x,
  #                   tsne.rot
  if ((.hasSlot(object, "dr"))) {
    if(!is.null(object@dr$pca)){
      pca.obj <- new(
        Class = "dim.reduction",
        gene.loadings = object@dr$pca@x,
        gene.loadings.full = object@dr$pca@x.full,
        cell.embeddings = object@dr$pca@rotation,
        sdev = object@dr$pca@sdev,
        key = "PC",
        misc = object@dr$pca@misc
      )
      new.object@dr$pca <- pca.obj
    }
   if(!is.null(object@dr$ica)){
     ica.obj <- new(
       Class = "dim.reduction",
       gene.loadings = object@dr$ica@x,
       cell.embeddings = object@dr$ica@rotation,
       sdev = object@dr$ica@sdev,
       key = "IC",
       misc = object@dr$ica@misc
     )
     new.object@dr$ica <- ica.obj
   }
  if(!is.null(object@dr$tsne)){
    tsne.obj <- new(
      Class = "dim.reduction",
      cell.embeddings = object@dr$tsne@rotation,
      key = "tSNE_"
    )
    new.object@dr$tsne <- tsne.obj
  }
  }
  # Conversion from release versions prior to 2.0.0
  # Slots to replace: pca.x, pca.rot, pca.x.full, tsne.rot, ica.rot, ica.x,
  #                   tsne.rot
  else{
    pca.sdev <- object@pca.obj[[1]]$sdev
    if (is.null(x = pca.sdev)) {
      pca.sdev <- object@pca.obj[[1]]$d
    }
    pca.obj <- new(
      Class = "dim.reduction",
      gene.loadings = as.matrix(object@pca.x),
      gene.loadings.full = as.matrix(object@pca.x.full),
      cell.embeddings = as.matrix(object@pca.rot),
      sdev = pca.sdev,
      key = "PC"
    )
    new.object@dr$pca <- pca.obj
    ica.obj <- new(
      Class = "dim.reduction",
      gene.loadings = as.matrix(object@ica.x),
      cell.embeddings = as.matrix(object@ica.rot),
      key = "IC"
    )
    new.object@dr$ica <- ica.obj
    tsne.obj <- new(
      Class = "dim.reduction",
      cell.embeddings = as.matrix(object@tsne.rot),
      key = "tSNE_"
    )
    new.object@dr$tsne <- tsne.obj
  }

  if (length(x = object@snn.sparse) == 1 && length(x = object@snn.dense) > 1) {
    if (class(object@snn.dense) == "data.frame") {
      object@snn.dense <- as.matrix(x = object@snn.dense)
    }
    new.object@snn <- as(object = object@snn.dense, Class = "dgCMatrix")
  }
  else{
    new.object@snn <- object@snn.sparse
  }

  return(new.object)
}


#' Return a subset of rows for a matrix or data frame
#'
#' @param data Matrix or data frame with row names
#' @param code Pattern for matching within row names
#' @return Returns a subset of data, using only rownames that did not yield a match to the pattern
#' @export
minusr <- function(data, code) {
  matchCode <- rownames(x = data)[grep(pattern = code, x = rownames(x = data))]
  toIgnore <- which(x = rownames(x = data) %in% matchCode)
  if (length(x = toIgnore) == 0) {
    return(data)
  }
  toRet <- data.frame(data[-toIgnore, ])
  rownames(x = toRet) <- rownames(x = data)[-toIgnore]
  colnames(x = toRet) <- colnames(x = data)
  return(toRet)
}

#' Return a subset of columns for a matrix or data frame
#'
#' @param data Matrix or data frame with column names
#' @param code Pattern for matching within column names
#' @return Returns a subset of data, using only column names that did not yield a match to the pattern
#' @export
minusc <- function(data, code) {
  matchCode <- colnames(x = data)[grep(pattern = code, colnames(x = data))]
  toIgnore <- which(x = colnames(x = data) %in% matchCode)
  if (length(x = toIgnore) == 0) {
    return(data)
  }
  return(data[, -toIgnore])
}

#' Return a subset of rows for a matrix or data frame
#'
#' @param data Matrix or data frame with row names
#' @param code Pattern for matching within row names
#' @return Returns a subset of data, using only rownames that yielded a match to the pattern
#' @export
subr <- function(data, code) {
  return(data[grep(pattern = code, x = rownames(x = data)), ])
}

#' Independently shuffle values within each row of a matrix
#'
#' Creates a matrix where correlation structure has been removed, but overall values are the same
#'
#' @param x Matrix to shuffle
#' @return Returns a scrambled matrix, where each row is shuffled independently
#' @export
shuffleMatRow <- function(x) {
  x2 <- x
  x2 <- t(x = x)
  ind <- order(c(col(x = x2)), runif(n = length(x = x2)))
  x2 <- matrix(
    data = x2[ind],
    nrow = nrow(x = x),
    ncol = ncol(x = x),
    byrow = TRUE
  )
  return(x2)
}


#' #returns the setdiff of two vectors, i.e all elements of a that are not in b
#' #we should be switching to use the setdiff function instead
#' #' @export
#' setdiff <- function(x, y) {
#'   return(x[! x %in% y])
#' }

#' Return a subset of columns for a matrix or data frame
#'
#' @param data Matrix or data frame with column names
#' @param code Pattern for matching within column names
#' @return Returns a subset of data, using only column names that yield a match to the pattern
#' @export
subc <- function(data, code) {
  return(data[, grep(pattern = code, x = colnames(x = data))])
}

#' Apply a ceiling and floor to all values in a matrix
#'
#' @param data Matrix or data frame
#' @param min all values below this min value will be replaced with min
#' @param max all values above this max value will be replaced with max
#' @return Returns matrix after performing these floor and ceil operations
#' @export
minmax <- function(data, min, max) {
  data2 <- data
  data2[data2 > max] <- max
  data2[data2 < min] <- min
  return(data2)
}



#' Extract delimiter information from a string.
#'
#' Parses a string (usually a cell name) and extracts fields based on a delimiter
#'
#' @param string String to parse.
#' @param field Integer(s) indicating which field(s) to extract. Can be a vector multiple numbers.
#' @param delim Delimiter to use, set to underscore by default.
#'
#' @return A new string, that parses out the requested fields, and (if multiple), rejoins them with the same delimiter
#' @export
extract_field <- function(string, field = 1, delim = "_") {
  fields <- as.numeric(x = unlist(x = strsplit(x = as.character(x = field), split = ",")))
  if (length(x = fields) == 1) {
    return(strsplit(x = string, split = delim)[[1]][field])
  }
  return(paste(strsplit(x = string, split = delim)[[1]][fields], collapse = delim))
}

#' Calculate the variance of logged values
#'
#' Calculate variance of logged values in non-log space (return answer in
#' log-space)
#'
#' @param x value or vector of values
#'
#' @return Returns the variance in log-space
#' @export
expVar <- function(x) {
  return(log1p(var(expm1(x))))
}

#' Calculate the standard deviation of logged values
#'
#' Calculate SD of logged values in non-log space (return answer in log-space)
#'
#' @param x value or vector of values
#'
#' @return Returns the standard deviation in log-space
#' @export
expSD <- function(x) {
  return(log1p(sd(expm1(x))))
}

#' Calculate the mean of logged values
#'
#' Calculate mean of logged values in non-log space (return answer in log-space)
#'
#' @param x value or vector of values
#'
#' @return Returns the mean in log-space
#' @export
expMean <- function(x) {
  return(log(x = mean(x = exp(x = x) - 1) + 1))
}

#' Calculate the dispersion of logged values
#'
#' Calculate the dispersion (variance divided by mean) in non-log space
#' (return answer in log-space)
#'
#' @param x value or vector of values
#'
#' @return Returns the dispersion in log-space
#' @export
logVarDivMean <- function(x) {
  return(log(x = var(x = exp(x = x) - 1) / mean(x = exp(x = x) - 1)))
}

#   Run a custom distance function on an input data matrix
#from Jean Fan - thanks!!
#' @export
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
AverageDetectionRate <- function(object, thresh.min = 0) {
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
  node.children <- intersect(x = node.children, y = levels(x = object@ident))
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
  genes.fit <- SetIfNull(x = genes.fit, default = object@var.genes)
  genes.fit <- genes.fit[genes.fit %in% rownames(x = object@data)]
  dim.code <- GetDimReduction(
    object = object,
    reduction.type = reduction.use,
    slot = 'key'
  )
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
#' object@@imputed
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
  genes.use <- SetIfNull(x = genes.use, default = object@var.genes)
  genes.fit <- SetIfNull(x = genes.fit, default = object@var.genes)
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

#' GenesInCluster
#'
#' After k-means analysis, previously run with DoKMeans, returns a set of genes associated with each cluster
#'
#' @param object Seurat object. Assumes DoKMeans has already been run
#' @param cluster.num K-means cluster(s) to return genes for
#' @param max.genes max number of genes to return
#' @return A vector of genes who are members in the cluster.num k-means cluster(s)
#'
#' @export
GenesInCluster <- function(object, cluster.num, max.genes = 1e6) {
  toReturn <- unlist(
    x = lapply(
      X = cluster.num,
      FUN = function(x) {
        return(head(
          x = sort(x = names(x = which(x = object@kmeans@gene.kmeans.obj$cluster==x))),
          n = max.genes
        ))
      }
    )
  )
  return(toReturn)
}


