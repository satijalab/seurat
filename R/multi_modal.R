#' @include seurat.R
NULL

# Set up assay class to hold multimodal data sets

assay <- setClass(
  Class = "assay",
  slots = list(
    raw.data = "ANY",
    data = "ANY",
    scale.data = "ANY",
    key = "character",
    misc = "ANY",
    var.genes="vector",
    mean.var="data.frame"
  )
)

#' Accessor function for multimodal data
#'
#' Pull information for specified stored dimensional reduction analysis
#'
#' @param object Seurat object
#' @param assay.type Type of assay to fetch data for (default is RNA)
#' @param slot Specific information to pull (i.e. raw.data, data, scale.data,...). Default is data
#'
#' @return Returns assay data
#'
#' @export
#'
#' @examples
#' # Simulate CITE-Seq results
#' df <- t(x = data.frame(
#'   x = round(x = rnorm(n = 80, mean = 20, sd = 2)),
#'   y = round(x = rbinom(n = 80, size = 100, prob = 0.2))
#' ))
#' pbmc_small <- SetAssayData(
#'   object = pbmc_small,
#'   assay.type = 'CITE',
#'   new.data = df,
#'   slot = 'raw.data'
#' )
#' GetAssayData(object = pbmc_small, assay.type = 'CITE', slot = 'raw.data')
#'
GetAssayData <- function(object, assay.type = "RNA", slot = "data") {
  if (assay.type == "RNA") {
    if (slot == "raw.data") {
      to.return <- object@raw.data
    } else if (slot == "data") {
      to.return <- object@data
    } else if (slot == "scale.data") {
      if (length(x = object@scale.data) == 0) {
        stop("Object@scale.data has not been set. Run ScaleData() and then retry.")
      }
      to.return <- object@scale.data
    }
    if(is.null(to.return)) {
      return(to.return)
    }
    #note that we check for this to avoid a long subset for large matrices if it can be avoided
    if (length(x = object@cell.names) == ncol(to.return)) {
      return(to.return)
    }
    return(to.return[, object@cell.names])
  }
  if (! (assay.type %in% names(object@assay))) {
    stop(paste(assay.type, "data has not been added"))
  }
  if (! (slot %in% slotNames(eval(expr = parse(text = paste0("object@assay$", assay.type)))))) {
    stop(paste(slot, "slot doesn't exist"))
  }
  to.return <- (eval(expr = parse(text = paste0("object@assay$", assay.type, "@", slot))))
  if (length(x = object@cell.names) == ncol(x = to.return)) {
    return(to.return)
  }
  return(to.return[, object@cell.names])
}

#' Assay Data Mutator Function
#'
#' Store information for specified assay, for multimodal analysis
#'
#' @inheritParams GetAssayData
#' @param new.data New data to insert
#'
#' @return Seurat object with updated slot
#'
#' @export
#'
#' @examples
#' # Simulate CITE-Seq results
#' df <- t(x = data.frame(
#'   x = round(x = rnorm(n = 80, mean = 20, sd = 2)),
#'   y = round(x = rbinom(n = 80, size = 100, prob = 0.2))
#' ))
#' pbmc_small = SetAssayData(
#'   object = pbmc_small,
#'   assay.type = 'CITE',
#'   new.data = df,
#'   slot = 'raw.data'
#' )
#'
SetAssayData <- function(object, assay.type, slot, new.data) {
  if (assay.type == "RNA") {
    if (slot == "raw.data") {
      (object@raw.data <- new.data)
    } else if (slot == "data") {
      (object@data <- new.data)
    } else if (slot == "scale.data") {
      (object@scale.data <- new.data)
    }
    return(object)
  }
  if (assay.type %in% names(object@assay)) {
    eval(expr = parse(text = paste0("object@assay$", assay.type, "@", slot, "<- new.data")))
  } else {
    new.assay <- new(Class = "assay")
    eval(expr = parse(text = paste0("new.assay@", slot, "<- new.data")))
    eval(expr = parse(text = paste0("object@assay$", assay.type, "<- new.assay")))
  }
  return(object)
}
