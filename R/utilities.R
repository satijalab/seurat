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
