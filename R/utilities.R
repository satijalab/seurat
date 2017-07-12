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


#returns the setdiff of two vectors, i.e all elements of a that are not in b
#we should be switching to use the setdiff function instead
#' @export
anotinb <- function(x, y) {
  return(x[! x %in% y])
}

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
