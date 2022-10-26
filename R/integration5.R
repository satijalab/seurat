#' @include zzz.R
#' @include generics.R
#'
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Generics
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Harmony Integration
#'
#' @inheritParams harmony::HarmonyMatrix
#' @param object An \code{\link[SeuratObject]{Assay5}} object
# @param assay Name of \code{object} in the containing \code{Seurat} object
#' @param orig \link[SeuratObject:DimReduc]{Dimensional reduction} to correct
#' @param groups A data frame ...
#' @param features ...
#' @param scale.layer ...
#' @param layers ...
#' @param ... Ignored
#'
#' @return ...
#'
#' @note This function requires the
#' \href{https://cran.r-project.org/package=harmony}{\pkg{harmony}} package
#' to be installed
#'
# @templateVar pkg harmony
# @template note-reqdpkg
#'
#' @export
#'
#' @concept integration
#'
#' @seealso \code{\link[harmony:HarmonyMatrix]{harmony::HarmonyMatrix}()}
#'
HarmonyIntegration <- function(
  object,
  orig,
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
  verbose = TRUE,
  ...
) {
  check_installed(
    pkg = "harmony",
    reason = "for running integration with Harmony"
  )
  if (!inherits(x = object, what = 'StdAssay')) {
    abort(message = "'object' must be a v5 assay object")
  } else if (!inherits(x = orig, what = 'DimReduc')) {
    abort(message = "'orig' must be a dimensional reduction")
  }
  # # Run joint PCA
  # features <- features %||% Features(x = object, layer = scale.layer)
  # pca <- RunPCA(
  #   object = object,
  #   assay = assay,
  #   features = features,
  #   layer = scale.layer,
  #   npcs = npcs,
  #   verbose = verbose
  # )
  # Run Harmony
  harmony.embed <- harmony::HarmonyMatrix(
    data_mat = Embeddings(object = orig),
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
  rownames(x = harmony.embed) <- Cells(x = orig)
  # TODO add feature loadings from PCA
  dr <- suppressWarnings(expr = CreateDimReducObject(
    embeddings = harmony.embed,
    key = key,
    # assay = assay
    assay = DefaultAssay(object = orig)
  ))
  return(list(harmony = dr))
}

attr(x = HarmonyIntegration, which = 'Seurat.method') <- 'integration'

#' Integrate Layers
#'
#' @param object A \code{\link[SeuratObject]{Seurat}} object
#' @param method Integration method function
#' @param orig Name of dimensional reduction for correction
#' @param group.by Name of meta data to group cells by; defaults to splits
#' assay layers
#' @param assay Name of assay for integration
#' @param features A vector of features to use for integration
#' @param layers Names of normalized layers in \code{assay}
#' @param scale.layer Name(s) of scaled layer(s) in \code{assay}
#' @param ... Arguments passed on to \code{method}
#'
#' @return \code{object} with integration data added to it
#'
#' @section Integration Method Functions:
#' The following integration method functions are available:
#' \Sexpr[stage=render,results=rd]{Seurat::.rd_methods("integration")}
#'
#' @export
#'
#' @concept integration
#'
#' @seealso \link[Seurat:writing-integration]{Writing integration method functions}
#'
IntegrateLayers <- function(
  object,
  method,
  orig = NULL,
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
  # Check our dimensional reduction
  orig <- orig %||% DefaultDimReduc(object = object, assay = assay)
  if (!orig %in% Reductions(object = object)) {
    abort(message = paste(sQuote(x = orig), 'is not a dimensional reduction'))
  }
  obj.orig <- object[[orig]]
  if (is.null(x = DefaultAssay(object = obj.orig))) {
    DefaultAssay(object = obj.orig) <- assay
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
    orig = obj.orig,
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
# Methods for Seurat-defined generics
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Internal
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Writing Integration Method Functions
#'
#' Integration method functions can be written by anyone to implement any
#' integration method in Seurat. These methods should expect to take a
#' \link[SeuratObject:Assay5]{v5 assay} as input and return a named list of
#' objects that can be added back to a \code{Seurat} object (eg. a
#' \link[SeuratObject:DimReduc]{dimensional reduction} or cell-level meta data)
#'
#' @section Provided Parameters:
#' Every integration method function should expect the following arguments:
#' \itemize{
#'  \item \dQuote{\code{object}}: an \code{\link[SeuratObject]{Assay5}} object
#  \item \dQuote{\code{assay}}: name of \code{object} in the original
#  \code{\link[SeuratObject]{Seurat}} object
#'  \item \dQuote{\code{orig}}: \link[SeuratObject:DimReduc]{dimensional
#'  reduction} to correct
#'  \item \dQuote{\code{layers}}: names of normalized layers in \code{object}
#'  \item \dQuote{\code{scale.layer}}: name(s) of scaled layer(s) in
#'  \code{object}
#'  \item \dQuote{\code{features}}: a vector of features for integration
#'  \item \dQuote{\code{groups}}: a one-column data frame with the groups for
#'  each cell in \code{object}; the column name will be \dQuote{group}
#' }
#'
#' @section Method Discovery:
#' The documentation for \code{\link{IntegrateLayers}()} will automatically
#' link to integration method functions provided by packages in the
#' \code{\link[base]{search}()} space. To make an integration method function
#' discoverable by the documentation, simply add an attribute named
#' \dQuote{\code{Seurat.method}} to the function with a value of
#' \dQuote{\code{integration}}
#' \preformatted{
#' attr(MyIntegrationFunction, which = "Seurat.method") <- "integration"
#' }
#'
#' @keywords internal
#'
#' @concept integration
#'
#' @name writing-integration
#' @rdname writing-integration
#'
#' @seealso \code{\link{IntegrateLayers}()}
#'
NULL
