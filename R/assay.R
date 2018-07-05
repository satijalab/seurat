# Methods for Assay objects
#' @include objects.R
#' @include assay_generics.R
#' @importFrom methods slot slot<- setMethod
NULL

#' @export
#'
MakeAssayObject <- function(
  raw.data,
  min.cells = 0,
  min.genes = 0,
  is.expr = 0,
  names.field = 1,
  names.delim = '_',
  ...
) {
  if (!inherits(x = raw.data, what = 'dgCMatrix')) {
    raw.data <- as(object = as.matrix(x = raw.data), Class = 'dgCMatrix')
  }
  if (is.expr > 0) {
    # suppress Matrix package note:
    # Note: method with signature 'CsparseMatrix#Matrix#missing#replValue' chosen for function '[<-',
    # target signature 'dgCMatrix#lgeMatrix#missing#numeric'.
    # "Matrix#ldenseMatrix#missing#replValue" would also be valid
    suppressMessages(expr = raw.data[raw.data < is.expr] <- 0)
  }
  # TODO: Add ngene and nUMI here
  # Filter based on min.genes
  num.genes <- colSums(x = raw.data > is.expr)
  raw.data <- raw.data[, which(x = num.genes > min.genes)]
  # filter genes on the number of cells expressing
  if (min.cells > 0) {
    num.cells <- rowSums(x = raw.data > 0)
    raw.data <- raw.data[which(x = num.cells >= min.cells), ]
  }
  # Set idents
  idents <- factor(x = unlist(x = lapply(
    X = colnames(x = raw.data),
    FUN = ExtractField,
    field = names.field,
    delim = names.delim
  )))
  names(x = idents) <- colnames(x = raw.data)
  # if there are more than 100 idents, set all idents to ... name
  ident.levels <- length(x = unique(x = idents))
  if (ident.levels > 100 || ident.levels == 0 || ident.levels == length(x = idents)) {
    idents <- factor(x = RandomName()) # TODO: change this
  }
  assay <- new(
    Class = 'Assay',
    raw.data = raw.data,
    data = raw.data,
    ident = idents
  )
  return(assay)
}

#' @describeIn GetAssayData Get assay data for an Assay object
#' @export
#' @method GetAssayData Assay
#'
GetAssayData.Assay <- function(object, slot = 'data') {
  return(slot(object = object, name = slot))
}

#' @describeIn SetAssayData Set assay data for an Assay object
#' @export
#' @method SetAssayData Assay
#'
SetAssayData.Assay <- function(object, slot, new.data) {
  slots.use <- c('raw.data', 'data', 'scale.data')
  if (!slot %in% slots.use) {
    stop("'slot' must be one of ", paste(slots.use, collapse = ', '))
  }
  slot(object = object, name = slot) <- new.data
  return(object)
}

#' @export
#' @method dim Assay
#'
dim.Assay <- function(x) {
  return(c(
    nrow(x = GetAssayData(object = x, slot = 'data')),
    ncol(x = GetAssayData(object = x, slot = 'data'))
  ))
}

setMethod(
  f = 'show',
  signature = 'Assay',
  definition = function(object) {
    cat(
      'Seurat assay data with',
      nrow(x = object),
      'features for',
      ncol(x = object), 'cells\n'
    )
    if (length(x = object@var.features) > 0) {
      cat(
        "Top 10 variable features:\n",
        strwrap(x = paste(head(x = object@var.features, n = 10L), collapse = ', '))
      )
    }
  }
)
