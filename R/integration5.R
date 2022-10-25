#' @include zzz.R
#' @include generics.R
#'
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Generics
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Integrate Layers
#'
#' @param object A \code{\link[SeuratObject]{Seurat}} object
#' @param method Integration method function; can choose from:
#' \Sexpr[stage=render,results=rd]{Seurat::.rd_methods("integration")}
#' @param ... Arguments passed on to \code{method}
#'
#' @return \code{object} with integration data added to it
#'
#' @section Integration Methods Functions:
#' Integration method functions can be written by anyone to implement any
#' integration method in Seurat. These methods should expect to take a
#' \link[SeuratObject:Assay5]{v5 assay} as input and return a named list of
#' objects that can be added back to a \code{Seurat} object (eg. a
#' \link[SeuratObject:DimReduc]{dimensional reduction} or cell-level meta data)
#'
#' Every integration method function should expect the following arguments:
#' \itemize{
#'  \item \dQuote{\code{object}}: an \code{\link[SeuratObject]{Assay5}} object
#'  \item \dQuote{\code{assay}}: name of \code{object} in the original
#'  \code{\link[SeuratObject]{Seurat}} object
#'  \item \dQuote{\code{layers}}: names of normalized layers in \code{object}
#'  \item \dQuote{\code{scale.layer}}: name(s) of scaled layer(s) in
#'  \code{object}
#'  \item \dQuote{\code{features}}: a vector of features for integration
#'  \item \dQuote{\code{groups}}: a one-column data frame with the groups for
#'  each cell in \code{object}; the column name will be \dQuote{group}
#' }
#'
#' @export
#'
IntegrateLayers <- function(
  object,
  method,
  group.by = NULL,
  assay = NULL,
  features = NULL,
  layers = NULL,
  scale.layer = 'scale.data',
  ...
) {
  # Get the integration method
  if (is_quosure(x = method)) {
    method <- eval(
      expr = quo_get_expr(quo = method),
      envir = quo_get_env(quo = method)
    )
  }
  if (is.character(x = method)) {
    method <- get(x = method)
  }
  if (!is.function(x = method)) {
    abort(message = "'method' must be a function for integrating layers")
  }
  # Check our assay
  assay <- assay %||% DefaultAssay(object = object)
  if (!inherits(x = object[[assay]], what = 'StdAssay')) {
    abort(message = "'assay' must be a v5 assay")
  }
  layers <- Layers(object = object, assay = assay, search = layers)
  features <- features %||% SelectIntegrationFeatures5(object = object, assay = assay)
  scale.layer <- Layers(object = object, search = scale.layer)
  features <- intersect(
    x = features,
    y = Features(x = object[[assay]], layer = scale.layer)
  )
  if (!length(x = features)) {
    abort(message = "None of the features provided are found in this assay")
  }
  # Check our groups
  groups <- if (is.null(x = group.by) && length(x = layers) > 1L) {
    cmap <- slot(object = object[[assay]], name = 'cells')[, layers]
    as.data.frame(x = labels(
      object = cmap,
      values = Cells(x = object[[assay]], layer = scale.layer)
    ))
  } else if (is_scalar_character(x = group.by)  && group.by %in% names(x = object[[]])) {
    FetchData(
      object = object,
      vars = group.by,
      cells = colnames(x = object[[assay]])
    )
  } else {
    abort(message = "'group.by' must correspond to a column of cell-level meta data")
  }
  names(x = groups) <- "group"
  # Run the integration method
  value <- method(
    object = object[[assay]],
    assay = assay,
    layers = layers,
    scale.layer = scale.layer,
    features = features,
    groups = groups,
    ...
  )
  for (i in names(x = value)) {
    object[[i]] <- value[[i]]
  }
  return(object)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Harmony Integration
#'
#' @export
#'
HarmonyIntegration <- function(
  object,
  assay,
  groups,
  features = NULL,
  scale.layer = 'scale.data',
  layers = NULL,
  npcs = 50L,
  key = 'harmony_',
  theta = NULL,
  lambda = NULL,
  sigma = 0.1,
  nclust = NULL,
  tau = 0,
  block.size = 0.05,
  max.iter.harmony = 10L,
  max.iter.cluster = 20L,
  epsilon.cluster = 1e-05,
  epsilon.harmony = 1e-04,
  # project.dim = TRUE,
  verbose = TRUE,
  ...
) {
  check_installed(
    pkg = "harmony",
    reason = "for running integration with Harmony"
  )
  if (!inherits(x = object, what = 'StdAssay')) {
    abort(message = "'object' must be a v5 assay object")
  }
  # Run joint PCA
  features <- features %||% Features(x = object, layer = scale.layer)
  pca <- RunPCA(
    object = object,
    assay = assay,
    features = features,
    layer = scale.layer,
    npcs = npcs,
    verbose = verbose
  )
  # Run Harmony
  harmony.embed <- harmony::HarmonyMatrix(
    data_mat = Embeddings(object = pca),
    meta_data = groups,
    vars_use = 'group',
    do_pca = FALSE,
    npcs = 0L,
    theta = theta,
    lambda = lambda,
    sigma = sigma,
    nclust = nclust,
    tau = tau,
    block.size = block.size,
    max.iter.harmony = max.iter.harmony,
    max.iter.cluster = max.iter.cluster,
    epsilon.cluster = epsilon.cluster,
    epsilon.harmony = epsilon.harmony,
    return_object = FALSE,
    verbose = verbose
  )
  rownames(x = harmony.embed) <- Cells(x = pca)
  # TODO add feature loadings from PCA
  dr <- suppressWarnings(expr = CreateDimReducObject(
    embeddings = harmony.embed,
    key = key,
    assay = assay
  ))
  return(list(pca = pca, harmony = dr))
}

attr(x = HarmonyIntegration, which = 'Seurat.method') <- 'integration'

#' Seurat-CCA Integration
#' 
#' @inheritParams FindIntegrationAnchors
#' @export
#'

CCAIntegration <- function(
    object = NULL,
    assay = NULL,
    layers = NULL,
    orig.reduction = 'pca.rna',
    new.reduction = 'integrated.dr',
    reference = NULL,
    anchor.features = NULL,
    normalization.method = c("LogNormalize", "SCT"),
    dims = 1:30,
    k.filter = NA,
    scale.data.layer = 'scale.data',
    verbose = TRUE,
    ...) {
  anchor.features <- anchor.features %||% SelectIntegrationFeatures5(object = object)
  assay <- assay %||% DefaultAssay(object = object)
  layers <- layers %||% Layers(object, search = 'data')
  object <- RunPCA(object = object,
                   assay = assay,
                   features = anchor.features,
                   reduction.name = orig.reduction,
                   reduction.key = paste0(orig.reduction,"_"),
                   verbose = verbose
  )
  
  object.list <- list()
  for (i in seq_along(along.with = layers)) {
    object.list[[i]] <- CreateSeuratObject(counts = object[[assay]][[layers[i]]] )
    object.list[[i]][['RNA']][[scale.data.layer]] <- object[[assay]]$scale.data[,Cells(object.list[[i]])]
    object.list[[i]][['RNA']]$counts <- NULL
  }
  anchor <- FindIntegrationAnchors(object.list = object.list, 
                                   anchor.features = anchor.features, 
                                   scale = FALSE, 
                                   reduction = 'cca', 
                                   normalization.method = normalization.method,
                                   dims = dims,
                                   k.filter = k.filter,
                                   reference = reference,
                                   verbose = verbose,
                                   ...
  )

  ## diet Seurat object 
  ###

  ###
  object_merged <- IntegrateEmbeddings(anchorset = anchor,
                                       reductions = object[[orig.reduction]],
                                       new.reduction.name = new.reduction,
                                       verbose = verbose)
  object[[new.reduction]] <- object_merged[[new.reduction]]
  return(object)
}




#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Methods for Seurat-defined generics
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# @method IntegrateLayers StdAssay
# @export
#'
IntegrateLayers.StdAssay <- function(
  object,
  method,
  assay,
  features, # TODO: allow selectintegrationfeatures to run on v5 assays
  layers = NULL,
  group.by = NULL,
  ...
) {
  if (is_quosure(x = method)) {
    method <- eval(
      expr = quo_get_expr(quo = method),
      envir = quo_get_env(quo = method)
    )
  }
  if (is.character(x = method)) {
    method <- get(x = method)
  }
  if (!is.function(x = method)) {
    abort(message = "'method' must be a function for integrating layers")
  }
  layers <- Layers(object = object, search = layers)
  if (!is.null(x = group.by)) {
    object <- split(x = object, f = group.by, drop = TRUE)
  }
  return(method(object = object, features = features, assay = assay, layers = layers, ...))
}

# @method IntegrateLayers Seurat
# @export
#'
IntegrateLayers.Seurat <- function(
  object,
  method,
  group.by = NULL,
  assay = NULL,
  features = NULL,
  layers = NULL,
  ...
) {
  method <- enquo(arg = method)
  assay <- assay %||% DefaultAssay(object = object)
  if (!inherits(x = object[[assay]], what = 'StdAssay')) {
    abort(message = "'assay' must be a v5 assay")
  }
  features <- features %||% SelectIntegrationFeatures5(object = object, assay = assay)
  group.by <- if (is.null(x = group.by)) {
    group.by
  } else if (rlang::is_na(x = group.by)) {
    Idents(object = object)
  } else if (is_scalar_character(x = group.by) && group.by %in% names(x = object[[]])) {
    object[[group.by, drop = TRUE]]
  } else {
    abort(message = "'group.by' must be the name of a column in cell-level meta data")
  }
  value <- IntegrateLayers(
    object = object[[assay]],
    method = method,
    assay = assay,
    features = features,
    layers = layers,
    group.by = group.by,
    ...
  )
  for (i in names(x = value)) {
    object[[i]] <- value[[i]]
  }
  return(object)
}



#' Seurat-RPCA Integration
#' 
#' @inheritParams FindIntegrationAnchors
#' @export
#'

RPCAIntegration <- function(
    object = NULL,
    assay = NULL,
    layers = NULL,
    orig.reduction = 'pca.rna',
    new.reduction = 'integrated.dr',
    reference = NULL,
    anchor.features = NULL,
    normalization.method = c("LogNormalize", "SCT"),
    dims = 1:30,
    k.filter = NA,
    scale.data.layer = 'scale.data',
    verbose = TRUE,
    ...) {
  anchor.features <- anchor.features %||% SelectIntegrationFeatures5(object = object)
  assay <- assay %||% DefaultAssay(object = object)
  layers <- layers %||% Layers(object, search = 'data')
  object <- RunPCA(object = object,
                   assay = assay,
                   features = anchor.features,
                   reduction.name = orig.reduction,
                   reduction.key = paste0(orig.reduction,"_"),
                   verbose = verbose
  )
  
  object.list <- list()
  for (i in seq_along(along.with = layers)) {
    object.list[[i]] <- CreateSeuratObject(counts = object[[assay]][[layers[i]]] )
    VariableFeatures(object =  object.list[[i]]) <- anchor.features
    object.list[[i]] <- ScaleData( object.list[[i]], verbose = FALSE)
    object.list[[i]] <- RunPCA( object.list[[i]], verbose = FALSE)
    object.list[[i]][['RNA']]$counts <- NULL
  }

  anchor <- FindIntegrationAnchors(object.list = object.list, 
                                   anchor.features = anchor.features, 
                                   scale = FALSE, 
                                   reduction = 'rpca', 
                                   normalization.method = normalization.method,
                                   dims = dims,
                                   k.filter = k.filter,
                                   reference = reference,
                                   verbose = verbose,
                                   ...
  )
  ## diet Seurat object 
  ###
  
  ###
  object_merged <- IntegrateEmbeddings(anchorset = anchor,
                                       reductions = object[[orig.reduction]],
                                       new.reduction.name = new.reduction,
                                       verbose = verbose)
  object[[new.reduction]] <- object_merged[[new.reduction]]
  return(object)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Internal
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
