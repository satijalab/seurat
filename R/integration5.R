#' @include zzz.R
#' @include generics.R
#'
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Generics
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @export
#'
IntegrateLayers <- function(object, ...) {
  UseMethod(generic = 'IntegrateLayers', object = object)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @export
#'
HarmonyIntegration <- function(
  object,
  features,
  assay,
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
  project.dim = TRUE,
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
  layers <- Layers(object = object, search = layers)
  object <- ScaleData(
    object = object,
    features = features,
    layer = layers,
    verbose = verbose
  )
  pca <- RunPCA(object = object, npcs = npcs, verbose = verbose)
  cmap <- slot(object = object, name = 'cells')[, layers]
  md <- as.data.frame(x = labels(
    object = cmap,
    values = Cells(x = object, layer = 'scale.data'),
    select = 'first'
  ))
  names(x = md) <- 'layer'
  md <- md[Cells(x = object, layer = 'scale.data'), , drop = FALSE]
  harmony.embed <- harmony::HarmonyMatrix(
    data_mat = Embeddings(object = pca),
    meta_data = md,
    vars_use = 'layer',
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
  dr <- suppressWarnings(expr = CreateDimReducObject(
    embeddings = harmony.embed,
    key = key,
    assay = assay
  ))
  if (isTRUE(x = project.dim)) {
    warn("projection")
  }
  return(list(harmony = dr))
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Methods for Seurat-defined generics
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @method IntegrateLayers StdAssay
#' @export
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
  return(method(object = object, features = features, assay = assay, ...))
}

#' @method IntegrateLayers Seurat
#' @export
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
  features <- SelectIntegrationFeatures5(object = object, assay = assay)
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

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Internal
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
