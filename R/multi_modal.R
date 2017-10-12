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
#' Store information for specified assay, for multimodal analysis. new.data 
#' needs to have cells as the columns and measurement features (e.g. genes, 
#' proteins, etc ...) as rows. Additionally, all the cell names in the new.data
#' must match the cell names in the object (object@@cell.names).
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
  if(! all(colnames(new.data) %in% object@cell.names)) {
    stop("Cell names in new assay data matrix don't exactly match the cell names of the object")
  }
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


#' Slim down a multi-species expression matrix, when only one species is primarily of interenst.
#'
#' Valuable for CITE-seq analyses, where we typically spike in rare populations of 'negative control' cells from a different species.
#'
#' @param data.matrix A UMI count matrix. Should contain rownames that start with the ensuing arguments prefix.1 or prefix.2
#' @param prefix.1 The prefix denoting rownames for the species of interest. Default is "HUMAN_". These rownames will have this prefix removed in the returned matrix.
#' @param prefix.controls The prefix denoting rownames for the species of 'negative control' cells. Default is "MOUSE_".
#' @param features.controls.toKeep How many of the most highly expressed (average) negative control features (by default, 100 mouse genes), should be kept? All other rownames starting with prefix.2 are discarded.
#' @return A UMI count matrix. Rownames that started with prefix.1 have this prefix discarded. For rownames starting with prefix.2, only the most highly expressed features are kept, and the prefix is kept. All other rows are retained.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' cbmc.rna.collapsed <- CollapseSpeciesExpressionMatrix(cbmc.rna)
#' }
CollapseSpeciesExpressionMatrix <- function(
  data.matrix,
  prefix.1 = "HUMAN_",
  prefix.controls = "MOUSE_",
  features.controls.toKeep = 100
) {
  data.matrix.1 <- SubsetRow(data = data.matrix, code = prefix.1)
  data.matrix.2 <- SubsetRow(data = data.matrix, code = prefix.controls)
  data.matrix.3 <- data.matrix[setdiff(
    x = rownames(x = data.matrix),
    y = c(rownames(x = data.matrix.1), rownames(x = data.matrix.2))
  ),]
  rownames(x = data.matrix.1) <- make.unique(names = sapply(
    X = rownames(x = data.matrix.1),
    FUN = function(x) {
      return(gsub(pattern = prefix.1, replacement = "", x = x))
    }
  ))
  control.sums <- rowSums(data.matrix.2)
  control.keep <- names(x = head(
    x = sort(x = control.sums, decreasing = TRUE),
    n = features.controls.toKeep
  ))
  control.matrix.keep <- data.matrix.2[control.keep, ]
  final.matrix <- rbind(data.matrix.1, control.matrix.keep, data.matrix.3)
  return(final.matrix)
}
