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
GetAssayData <- function(object, assay.type = "RNA", slot = "data") {
  if (assay.type == "RNA") {
    if (slot == "raw.data") {
      to.return <- object@raw.data
    } else if (slot == "data") {
      to.return <- object@data
    } else if (slot == "scale.data") {
      to.return <- object@scale.data
    }
    #note that we check for this to avoid a long subset for large matrices if it can be avoided
    if (length(x = object@cell.names) == ncol(to.return)) {
      return(to.return)
    }
    return(to.return[, object@cell.names])
  }
  if (! (assay.type %in% names(objectobject@assay))) {
    stop(paste(assay.type, "data has not been added"))
  }
  if (! (slot %in% slotNames(eval(envir = parse(text = paste0("object@assay$", assay.type)))))) {
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
  if (assay.type %in% names(objectobject@assay)) {
    eval(expr = parse(text = paste0("object@assay$", assay.type, "@", slot, "<- new.data")))
  } else {
    new.assay <- new(Class = "assay")
    eval(expr = parse(text = paste0("new.assay@", slot, "<- new.data")))
    eval(expr = parse(text = paste0("object@assay$", assay.type, "<- new.assay")))
  }
  return(object)
}

#' Normalize Assay Data
#'
#' Normalize data for a given assay
#'
#' @param object Seurat object
#' @param assay.type Type of assay to normalize for (default is RNA), but can be changed for multimodal analyses.
#' @param normalization.method Method for normalization. Default is log-normalization (LogNormalize). Other options include CLR (CLR), regularized NB normalization (NBReg; RNA only)
#'
#' @importFrom compositions clr
#'
#' @return Returns object after normalization. Normalized data is stored in data or scale.data slot, depending on the method
#'
#' @export
#'
NormalizeData <- function(
  object,
  assay.type = "RNA",
  normalization.method = "LogNormalize",
  slot = "data",
  scale.factor = 1e4,
  display.progress = TRUE,
  ...
) {
  if (normalization.method == "LogNormalize") {
    raw.data <- GetAssayData(
      object = object,
      assay.type = assay.type,
      slot = "raw.data"
    )
    if (is.null(x = raw.data)) {
      stop(paste("Raw data for", assay.type, "has not been set"))
    }
    normalized.data <- LogNormalize(
      data = raw.data,
      scale.factor = scale.factor,
      display.progress = display.progress
    )
    object <- SetAssayData(
      object = object,
      assay.type = assay.type,
      slot = "data",
      new.data = normalized.data
    )
  }
  if (normalization.method == "CLR") {
    raw.data <- GetAssayData(
      object = object,
      assay.type = assay.type,
      slot = "raw.data"
    )
    normalized.data <- t(x = as.matrix(x = t(x = clr(x = raw.data))))
    object <- SetAssayData(
      object = object,
      assay.type = assay.type,
      slot = "data",
      new.data = normalized.data
    )
    object <- SetAssayData(
      object = object,
      assay.type = assay.type,
      slot = "scale.data",
      new.data = normalized.data
    )
  }
  return(object)
}

#' Run Canonical Correlation Analysis (CCA) on multimodal data
#'
#' CCA finds a shared correlation structure betwen two different datasets, enabling integrated downstream analysis
#'
#' @param object Seurat object
#' @param assay.1 First assay for multimodal analysis. Default is RNA
#' @param assay.2 Second assay for multimodal analysis. Default is CITE for CITE-Seq analysis.
#' @param features.1 Features of assay 1 to consider (default is variable genes)
#' @param features.2 Features of assay 2 to consider (default is all features, i.e. for CITE-Seq, all antibodies)
#' @param num.cc Number of canonical correlations to compute and store. Default is 20, but will calculate less if either assay has <20 features.
#' @param normalize.variance Z-score the embedding of each CC to 1, so each CC contributes equally in downstream analysis (default is T)
#'
#' @importFrom PMA CCA
#'
#' @return Returns object after CCA, with results stored in dimensional reduction cca.assay1 (ie. cca.RNA) and cca.assay2. For example, results can be visualized using DimPlot(object,reduction.use="cca.RNA")
#'
#' @export
#'
MultiModal_CCA <- function(
  object,
  assay.1 = "RNA",
  assay.2 = "CITE",
  features.1 = NULL,
  features.2 = NULL,
  num.cc = 20,
  normalize.variance = TRUE
) {
  #first pull out data, define features
  data.1 <- GetAssayData(
    object = object,
    assay.type = assay.1,
    slot = "scale.data"
  )
  data.2 <- GetAssayData(
    object = object,
    assay.type = assay.2,
    slot = "scale.data"
  )
  if (is.null(x = features.1)) {
    if ((assay.1 == "RNA") && length(x = object@var.genes) > 0) {
      features.1 <- object@var.genes
    } else {
      features.1 <- rownames(x = data.1)
    }
  }
  features.2 <- set.ifnull(x = features.2, y = rownames(x = data.2))
  data.1 <- t(x = data.1[features.1, ])
  data.2 <- t(x = data.2[features.2, ])
  num.cc <- min(20, min(length(x = features.1), length(x = features.2)))
  cca.data <- list(data.1, data.2)
  names(x = cca.data) <- c(assay.1, assay.2)
  # now run CCA
  out <- CCA(
    x = cca.data[[1]],
    z = cca.data[[2]],
    typex = "standard",
    typez = "standard",
    K = num.cc,
    penaltyz = 1,
    penaltyx = 1
  )
  cca.output <- list(out$u, out$v)
  embeddings.cca <- list()
  for (i in 1:length(x = cca.data)) {
    assay.use <- names(x = cca.data)[i]
    rownames(x = cca.output[[i]]) <- colnames(x = cca.data[[i]])
    embeddings.cca[[i]] <- cca.data[[i]] %*% cca.output[[i]]
    colnames(x = embeddings.cca[[i]]) <- paste0(
      assay.use,
      "CC",
      1:ncol(x = embeddings.cca[[i]])
    )
    colnames(x = cca.output[[i]]) <- colnames(x = embeddings.cca[[i]])
    if (normalize.variance) {
      embeddings.cca[[i]] <- scale(x = embeddings.cca[[i]])
    }
    object <- SetDimReduction(
      object = object,
      reduction.type = paste0(assay.use, "CCA"),
      slot = "rotation",
      new.data = embeddings.cca[[i]]
    )
    object <- SetDimReduction(
      object = object,
      reduction.type = paste0(assay.use, "CCA"),
      slot = "key",
      new.data = paste0(assay.use, "CC")
    )
    object <- SetDimReduction(
      object = object,
      reduction.type = paste0(assay.use, "CCA"),
      slot = "x",
      new.data =  cca.output[[i]]
    )
  }
  return(object)
}

#' Run coinertia analysis on multimodal data
#'
#' CIA finds a shared correlation structure betwen two different datasets, enabling integrated downstream analysis
#'
#' @param object Seurat object
#' @param assay.1 First assay for multimodal analysis. Default is RNA
#' @param assay.2 Second assay for multimodal analysis. Default is CITE for CITE-Seq analysis.
#' @param features.1 Features of assay 1 to consider (default is variable genes)
#' @param features.2 Features of assay 2 to consider (default is all features, i.e. for CITE-Seq, all antibodies)
#' @param num.axes Number of principal axes to compute and store. Default is 20, but will calculate less if either assay has <20 features.
#' @param normalize.variance Return the normalized row scares, so each aexis contributes equally in downstream analysis (default is T)
#'
#' @importFrom made4 cia
#'
#' @return Returns object after CIA, with results stored in dimensional reduction cia.assay1 (ie. cia.RNA) and cia.assay2. For example, results can be visualized using DimPlot(object,reduction.use="cia.RNA")
#'
#' @export
#'
MultiModal_CIA <- function(
  object,
  assay.1 = "RNA",
  assay.2 = "CITE",
  features.1 = NULL,
  features.2 = NULL,
  num.axes = 20,
  normalize.variance = T
) {
  #first pull out data, define features
  data.1 <- GetAssayData(
    object = object,
    assay.type = assay.1,
    slot = "scale.data"
  )
  data.2 <- GetAssayData(
    object = object,
    assay.type = assay.2,
    slot = "scale.data"
  )
  if (is.null(x = features.1)) {
    if ((assay.1 == "RNA") && length(x = object@var.genes) > 0) {
      features.1 <- object@var.genes
    } else {
      features.1 <- rownames(x = data.1)
    }
  }
  features.2 <- set.ifnull(x = features.2, y = rownames(x = data.2))
  data.1 <- t(x = data.1[features.1, ])
  data.2 <- t(x = data.2[features.2, ])
  num.axes <- min(20, min(length(x = features.1), length(x = features.2)))
  cia.data <- list(data.1, data.2)
  names(x = cia.data) <- c(assay.1, assay.2)
  # now run cia
  out <- cia(
    df1 = t(x = cia.data[[1]]),
    df2 = t(x = cia.data[[2]]),
    cia.nf = num.axes
  )
  out <- out$coinertia
  cia.output <- list(as.matrix(x = out$c1), as.matrix(x = out$l1))
  embeddings.cia.norm <- list(as.matrix(x = out$mX), as.matrix(x = out$mY))
  embeddings.cia <- list(as.matrix(x = out$lX), as.matrix(x = out$lY))
  for (i in 1:length(x = cia.data)) {
    assay.use <- names(x = cia.data)[i]
    #rownames(cia.output[[i]])=colnames(cia.data[[i]])
    if (normalize.variance) {
      embeddings.cia[[i]] <- (embeddings.cia.norm[[i]])
    }
    colnames(x = embeddings.cia[[i]]) <- paste0(
      assay.use,
      "CI",
      1:ncol(x = embeddings.cia[[i]])
    )
    colnames(x = cia.output[[i]]) <- colnames(x = embeddings.cia[[i]])
    object <- SetDimReduction(
      object = object,
      reduction.type = paste("cia", assay.use, sep="_"),
      slot = "rotation",
      new.data = embeddings.cia[[i]]
    )
    object <- SetDimReduction(
      object = object,
      reduction.type = paste("cia", assay.use, sep="_"),
      slot = "key",
      new.data = paste0(assay.use," CI")
    )
    object <- SetDimReduction(
      object = object,
      reduction.type = paste("cia", assay.use, sep="_"),
      slot = "x",
      new.data =  cia.output[[i]]
    )
  }
  return(object)
}
