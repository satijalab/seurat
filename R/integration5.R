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
#' @param object An \code{\link[SeuratObject]{Assay5}} object
# @param assay Name of \code{object} in the containing \code{Seurat} object
#' @param orig A \link[SeuratObject:DimReduc]{dimensional reduction} to correct
#' @param features Ignored
#' @param scale.layer Ignored
#' @param new.reduction Name of new integrated dimensional reduction
#' @param layers Ignored
#' @param key Key for Harmony dimensional reduction
#' @param npcs If doing PCA on input matrix, number of PCs to compute
#' @param theta Diversity clustering penalty parameter
#' @param lambda Ridge regression penalty parameter
#' @param sigma Width of soft kmeans clusters
#' @param nclust Number of clusters in model
#' @param tau Protection against overclustering small datasets with large ones
#' @param block.size What proportion of cells to update during clustering
#' @param max.iter.harmony Maximum number of rounds to run Harmony
#' @param max.iter.cluster Maximum number of rounds to run clustering at each round of Harmony
#' @param epsilon.cluster Convergence tolerance for clustering round of Harmony
#' @param epsilon.harmony Convergence tolerance for Harmony
#' @param verbose Whether to print progress messages. TRUE to print, FALSE to suppress
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
#' @examples
#' \dontrun{
#' # Preprocessing
#' obj <- SeuratData::LoadData("pbmcsca")
#' obj[["RNA"]] <- split(obj[["RNA"]], f = obj$Method)
#' obj <- NormalizeData(obj)
#' obj <- FindVariableFeatures(obj)
#' obj <- ScaleData(obj)
#' obj <- RunPCA(obj)
#'
#' # After preprocessing, we integrate layers with added parameters specific to Harmony:
#' obj <- IntegrateLayers(object = obj, method = HarmonyIntegration, orig.reduction = "pca",
#'   new.reduction = 'harmony', verbose = FALSE)
#'
#' # Modifying Parameters
#' # We can also add arguments specific to Harmony such as theta, to give more diverse clusters
#' obj <- IntegrateLayers(object = obj, method = HarmonyIntegration, orig.reduction = "pca",
#'   new.reduction = 'harmony', verbose = FALSE, theta = 3)
#' # Integrating SCTransformed data
#' obj <- SCTransform(object = obj)
#' obj <- IntegrateLayers(object = obj, method = HarmonyIntegration,
#'   orig.reduction = "pca", new.reduction = 'harmony',
#'   assay = "SCT", verbose = FALSE)
#' }
#'
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
  features = NULL,
  scale.layer = 'scale.data',
  new.reduction = 'harmony',
  layers = NULL,
  npcs = NULL,
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
  epsilon.harmony = 0.01,
  verbose = TRUE,
  ...
) {  

  # If the `features` param has been set, inform the user that is ignored.
  if (!is.null(features)) {
    inform(
      "The `features` argument is ignored by `HarmonyIntegration`.",
      .frequency = "once",
      .frequency_id = "HarmonyIntegration@features"
    )
  }
  # If the `npcs` param has been set, inform the user that is ignored.
  if (!is.null(npcs)) {
    inform(
      "The `npcs` argument is ignored by `HarmonyIntegration`.",
      .frequency = "once",
      .frequency_id = "HarmonyIntegration@npcs"
    )
  }
  check_installed(
    pkg = "harmony",
    reason = "for running integration with Harmony"
  )
  if (!inherits(x = object, what = c('StdAssay', 'SCTAssay'))) {
    abort(message = "'object' must be a v5 or SCT assay")
  } else if (!inherits(x = orig, what = 'DimReduc')) {
    abort(message = "'orig' must be a dimensional reduction")
  }

  #create grouping variables
  groups <- CreateIntegrationGroups(
    object, 
    layers = layers, 
    scale.layer = scale.layer
  )

  # Set up advanced Harmony options.
  advanced_options <- harmony::harmony_options(
    tau = tau,
    block.size = block.size,
    max.iter.cluster = max.iter.cluster,
    epsilon.cluster = epsilon.cluster,
    epsilon.harmony = epsilon.harmony
  )

  # Run Harmony
  harmony.embed <- harmony::RunHarmony(
    data_mat = Embeddings(object = orig),
    meta_data = groups,
    vars_use = 'group',
    theta = theta,
    lambda = lambda,
    sigma = sigma,
    nclust = nclust,
    max_iter = max.iter.harmony,
    return_object = FALSE,
    verbose = verbose,
    .options = advanced_options
  )
  rownames(harmony.embed) <- Cells(orig)
  # TODO add feature loadings from PCA
  dr <- suppressWarnings(expr = CreateDimReducObject(
    embeddings = harmony.embed,
    key = key,
    assay = DefaultAssay(object = orig)
  ))
  output.list <- list(dr)
  names(output.list) <- c(new.reduction)
  return(output.list)
}

attr(x = HarmonyIntegration, which = 'Seurat.method') <- 'integration'

#' Seurat-CCA Integration
#'
#' @inheritParams RPCAIntegration
#' @export
#' @concept integration
#'
#' @examples
#' \dontrun{
#' # Preprocessing
#' obj <- SeuratData::LoadData("pbmcsca")
#' obj[["RNA"]] <- split(obj[["RNA"]], f = obj$Method)
#' obj <- NormalizeData(obj)
#' obj <- FindVariableFeatures(obj)
#' obj <- ScaleData(obj)
#' obj <- RunPCA(obj)
#'
#' # After preprocessing, we integrate layers.
#' obj <- IntegrateLayers(object = obj, method = CCAIntegration,
#'   orig.reduction = "pca", new.reduction = "integrated.cca",
#'   verbose = FALSE)
#'
#' # Modifying parameters
#' # We can also specify parameters such as `k.anchor` to increase the strength of integration
#' obj <- IntegrateLayers(object = obj, method = CCAIntegration,
#'   orig.reduction = "pca", new.reduction = "integrated.cca",
#'   k.anchor = 20, verbose = FALSE)
#'
#' # Integrating SCTransformed data
#' obj <- SCTransform(object = obj)
#' obj <- IntegrateLayers(object = obj, method = CCAIntegration,
#'   orig.reduction = "pca", new.reduction = "integrated.cca",
#'   assay = "SCT", verbose = FALSE)
#' }
#'
CCAIntegration <- function(
    object = NULL,
    assay = NULL,
    layers = NULL,
    orig = NULL,
    new.reduction = 'integrated.dr',
    reference = NULL,
    features = NULL,
    normalization.method = c("LogNormalize", "SCT"),
    dims = 1:30,
    k.filter = NA,
    scale.layer = 'scale.data',
    dims.to.integrate = NULL,
    k.weight = 100,
    weight.reduction = NULL,
    sd.weight = 1,
    sample.tree = NULL,
    preserve.order = FALSE,
    verbose = TRUE,
    ...
) {
  op <- options(Seurat.object.assay.version = "v3", Seurat.object.assay.calcn = FALSE)
  on.exit(expr = options(op), add = TRUE)
  normalization.method <- match.arg(arg = normalization.method)
  features <- features %||% SelectIntegrationFeatures5(object = object)
  assay <- assay %||% 'RNA'
  layers <- layers %||% Layers(object, search = 'data')
  if (normalization.method == 'SCT') {
    #create grouping variables
    groups <- CreateIntegrationGroups(object, layers = layers, scale.layer = scale.layer)
    object.sct <- CreateSeuratObject(counts = object, assay = 'SCT')
    object.sct$split <- groups[,1]
    object.list <- SplitObject(object = object.sct,split.by = 'split')
    object.list  <- PrepSCTIntegration(object.list, anchor.features = features)
  } else {
  object.list <- list()
  for (i in seq_along(along.with = layers)) {
    if (inherits(x = object[layers[i]], what = "IterableMatrix")) {
      warning("Converting BPCells matrix to dgCMatrix for integration ",
        "as on-disk CCA Integration is not currently supported", call. = FALSE, immediate. = TRUE)
      counts <- as(object = object[layers[i]][features, ],
                   Class = "dgCMatrix")
    }
    else {
      counts <- object[layers[i]][features, ]
    }
    object.list[[i]] <- CreateSeuratObject(counts = counts)
    if (inherits(x = object[scale.layer], what = "IterableMatrix")) {
      scale.data.layer <- as.matrix(object[scale.layer][features,
                                                          Cells(object.list[[i]])])
      object.list[[i]][["RNA"]]$scale.data <- scale.data.layer
    }
    else {
      object.list[[i]][["RNA"]]$scale.data <- object[scale.layer][features,
                                                                    Cells(object.list[[i]])]
    }
    object.list[[i]][['RNA']]$counts <- NULL
  }
  }

  anchor <- FindIntegrationAnchors(object.list = object.list,
                                   anchor.features = features,
                                   scale = FALSE,
                                   reduction = 'cca',
                                   normalization.method = normalization.method,
                                   dims = dims,
                                   k.filter = k.filter,
                                   reference = reference,
                                   verbose = verbose,
                                   ...
  )
  suppressWarnings({
    anchor@object.list <- lapply(anchor@object.list, function(x) {
      x <- DietSeurat(x, features = features[1:2])
      return(x)
    })
  }, classes = "dimWarning")
  object_merged <- IntegrateEmbeddings(anchorset = anchor,
                                       reductions = orig,
                                       new.reduction.name = new.reduction,
                                       dims.to.integrate = dims.to.integrate,
                                       k.weight = k.weight,
                                       weight.reduction = weight.reduction,
                                       sd.weight = sd.weight,
                                       sample.tree = sample.tree,
                                       preserve.order = preserve.order,
                                       verbose = verbose
  )
  output.list <- list(object_merged[[new.reduction]])
  names(output.list) <- c(new.reduction)
  return(output.list)
}

attr(x = CCAIntegration, which = 'Seurat.method') <- 'integration'

#' Seurat-RPCA Integration
#'
#' @param object A \code{Seurat} object
#' @param assay Name of \code{Assay} in the \code{Seurat} object
#' @param layers Names of layers in \code{assay}
#' @param orig A \link[SeuratObject]{DimReduc} to correct
#' @param new.reduction Name of new integrated dimensional reduction
#' @param reference A reference \code{Seurat} object
#' @param features A vector of features to use for integration
#' @param normalization.method Name of normalization method used: LogNormalize
#' or SCT
#' @param dims Dimensions of dimensional reduction to use for integration
#' @param k.filter Number of anchors to filter
#' @param scale.layer Name of scaled layer in \code{Assay}
#' @param verbose Print progress
#' @param ... Additional arguments passed to \code{FindIntegrationAnchors}
#'
#' @examples
#' \dontrun{
#' # Preprocessing
#' obj <- SeuratData::LoadData("pbmcsca")
#' obj[["RNA"]] <- split(obj[["RNA"]], f = obj$Method)
#' obj <- NormalizeData(obj)
#' obj <- FindVariableFeatures(obj)
#' obj <- ScaleData(obj)
#' obj <- RunPCA(obj)
#'
#' # After preprocessing, we run integration
#' obj <- IntegrateLayers(object = obj, method = RPCAIntegration,
#'   orig.reduction = "pca", new.reduction = 'integrated.rpca',
#'   verbose = FALSE)
#'
#' # Reference-based Integration
#' # Here, we use the first layer as a reference for integraion
#' # Thus, we only identify anchors between the reference and the rest of the datasets,
#' # saving computational resources
#' obj <- IntegrateLayers(object = obj, method = RPCAIntegration,
#'   orig.reduction = "pca", new.reduction = 'integrated.rpca',
#'   reference = 1, verbose = FALSE)
#'
#' # Modifying parameters
#' # We can also specify parameters such as `k.anchor` to increase the strength of
#' # integration
#' obj <- IntegrateLayers(object = obj, method = RPCAIntegration,
#'   orig.reduction = "pca", new.reduction = 'integrated.rpca',
#'   k.anchor = 20, verbose = FALSE)
#'
#' # Integrating SCTransformed data
#' obj <- SCTransform(object = obj)
#' obj <- IntegrateLayers(object = obj, method = RPCAIntegration,
#'   orig.reduction = "pca", new.reduction = 'integrated.rpca',
#'   assay = "SCT", verbose = FALSE)
#' }
#'

#' @inheritParams FindIntegrationAnchors
#' @inheritParams IntegrateEmbeddings
#' @param ... Arguments passed on to \code{FindIntegrationAnchors}
#' @export
#'
#' @concept utilities
#'
RPCAIntegration <- function(
    object = NULL,
    assay = NULL,
    layers = NULL,
    orig = NULL,
    new.reduction = 'integrated.dr',
    reference = NULL,
    features = NULL,
    normalization.method = c("LogNormalize", "SCT"),
    dims = 1:30,
    k.filter = NA,
    scale.layer = 'scale.data',
    dims.to.integrate = NULL,
    k.weight = 100,
    weight.reduction = NULL,
    sd.weight = 1,
    sample.tree = NULL,
    preserve.order = FALSE,
    verbose = TRUE,
    ...
) {
  op <- options(Seurat.object.assay.version = "v3", Seurat.object.assay.calcn = FALSE)
  on.exit(expr = options(op), add = TRUE)
  normalization.method <- match.arg(arg = normalization.method)
  features <- features %||% SelectIntegrationFeatures5(object = object)
  assay <- assay %||% 'RNA'
  layers <- layers %||% Layers(object = object, search = 'data')
  #check that there enough cells present
  ncells <- sapply(X = layers, FUN = function(x) {ncell <-  dim(object[x])[2]
  return(ncell) })
  if (min(ncells) < max(dims))  {
    abort(message = "At least one layer has fewer cells than dimensions specified, please lower 'dims' accordingly.")
  }
  if (normalization.method == 'SCT') {
    #create grouping variables
    groups <- CreateIntegrationGroups(object, layers = layers, scale.layer = scale.layer)
    object.sct <- CreateSeuratObject(counts = object, assay = 'SCT')
    object.sct$split <- groups[,1]
    object.list <- SplitObject(object = object.sct, split.by = 'split')
    object.list <- PrepSCTIntegration(object.list = object.list, anchor.features = features)
    object.list <- lapply(X = object.list, FUN = function(x) {
      x <- RunPCA(object = x, features = features, verbose = FALSE, npcs = max(dims))
      return(x)
    }
    )
  } else {
    object.list <- list()
    for (i in seq_along(along.with = layers)) {
      object.list[[i]] <- suppressMessages(suppressWarnings(
        CreateSeuratObject(counts = NULL, data = object[layers[i]][features,])
      ))
      VariableFeatures(object =  object.list[[i]]) <- features
      object.list[[i]] <- suppressWarnings(ScaleData(object = object.list[[i]], verbose = FALSE))
      object.list[[i]] <- RunPCA(object = object.list[[i]], verbose = FALSE, npcs=max(dims))
      suppressWarnings(object.list[[i]][['RNA']]$counts <- NULL)
    }
  }
  anchor <- FindIntegrationAnchors(object.list = object.list,
                                   anchor.features = features,
                                   scale = FALSE,
                                   reduction = 'rpca',
                                   normalization.method = normalization.method,
                                   dims = dims,
                                   k.filter = k.filter,
                                   reference = reference,
                                   verbose = verbose,
                                   ...
  )
  slot(object = anchor, name = "object.list") <- lapply(
      X = slot(
        object = anchor,
        name = "object.list"),
      FUN = function(x) {
      suppressWarnings(expr = x <- DietSeurat(x, features = features[1:2]))
    return(x)
  })
  object_merged <- IntegrateEmbeddings(anchorset = anchor,
                                       reductions = orig,
                                       new.reduction.name = new.reduction,
                                       dims.to.integrate = dims.to.integrate,
                                       k.weight = k.weight,
                                       weight.reduction = weight.reduction,
                                       sd.weight = sd.weight,
                                       sample.tree = sample.tree,
                                       preserve.order = preserve.order,
                                       verbose = verbose
  )

  output.list <- list(object_merged[[new.reduction]])
  names(output.list) <- c(new.reduction)
  return(output.list)
}

attr(x = RPCAIntegration, which = 'Seurat.method') <- 'integration'

#' Seurat-Joint PCA Integration
#'
#' @inheritParams RPCAIntegration
#' @inheritParams FindIntegrationAnchors
#' @inheritParams IntegrateEmbeddings
#' @param ... Arguments passed on to \code{FindIntegrationAnchors}
#' @export
#' @concept integration
#'
JointPCAIntegration <- function(
    object = NULL,
    assay = NULL,
    layers = NULL,
    orig = NULL,
    new.reduction = 'integrated.dr',
    reference = NULL,
    features = NULL,
    normalization.method = c("LogNormalize", "SCT"),
    dims = 1:30,
    k.anchor = 20,
    scale.layer = 'scale.data',
    dims.to.integrate = NULL,
    k.weight = 100,
    weight.reduction = NULL,
    sd.weight = 1,
    sample.tree = NULL,
    preserve.order = FALSE,
    verbose = TRUE,
    ...
) {
  op <- options(Seurat.object.assay.version = "v3", Seurat.object.assay.calcn = FALSE)
  on.exit(expr = options(op), add = TRUE)
  normalization.method <- match.arg(arg = normalization.method)
  features <- features %||% SelectIntegrationFeatures5(object = object)
  features.diet <- features[1:2]
  assay <- assay %||%  DefaultAssay(object)
  layers <- layers %||% Layers(object, search = 'data')

  if (normalization.method == 'SCT') {
    #create grouping variables
    groups <- CreateIntegrationGroups(object, layers = layers, scale.layer = scale.layer)
    object.sct <- CreateSeuratObject(counts = object, assay = 'SCT')
    object.sct <- DietSeurat(object = object.sct, features = features.diet)
    object.sct[['joint.pca']] <- CreateDimReducObject(
      embeddings = Embeddings(object = orig),
      assay = 'SCT',
      loadings = Loadings(orig),
      key = 'J_'
    )
    object.sct$split <- groups[,1]
    object.list <- SplitObject(object = object.sct,split.by = 'split')
    object.list  <- PrepSCTIntegration(object.list, anchor.features = features.diet)
    object.list <- lapply(object.list, function(x) {
      x[['SCT']]@SCTModel.list <- list()
      return(x)
    })
    } else {
    object.list <- list()
    for (i in seq_along(along.with = layers)) {
      object.list[[i]] <- CreateSeuratObject(counts = object[layers[i]][features.diet, ] )
      object.list[[i]][['RNA']]$counts <- NULL
      object.list[[i]][['joint.pca']] <- CreateDimReducObject(
        embeddings = Embeddings(object = orig)[Cells(object.list[[i]]),],
        assay = 'RNA',
        loadings = Loadings(orig),
        key = 'J_'
      )
    }
  }

  anchor <- FindIntegrationAnchors(object.list = object.list,
                                   anchor.features = features.diet,
                                   scale = FALSE,
                                   reduction = 'jpca',
                                   normalization.method = normalization.method,
                                   dims = dims,
                                   k.anchor = k.anchor,
                                   k.filter = NA,
                                   reference = reference,
                                   verbose = verbose,
                                   ...
  )
  object_merged <- IntegrateEmbeddings(anchorset = anchor,
                                       reductions = orig,
                                       new.reduction.name = new.reduction,
                                       dims.to.integrate = dims.to.integrate,
                                       k.weight = k.weight,
                                       weight.reduction = weight.reduction,
                                       sd.weight = sd.weight,
                                       sample.tree = sample.tree,
                                       preserve.order = preserve.order,
                                       verbose = verbose
  )
  output.list <- list(object_merged[[new.reduction]])
  names(output.list) <- c(new.reduction)
  return(output.list)
}

attr(x = JointPCAIntegration, which = 'Seurat.method') <- 'integration'

#' Integrate Layers
#'
#' @param object A \code{\link[SeuratObject]{Seurat}} object
#' @param method Integration method function
#' @param orig.reduction Name of dimensional reduction for correction
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
#' \Sexpr[stage=render,results=rd]{Seurat:::.rd_methods("integration")}
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
  orig.reduction = 'pca',
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
  if (inherits(x = object[[assay]], what = 'SCTAssay')) {
    layers <- 'data'
    scale.layer <- 'scale.data'
    features <- features %||% SelectSCTIntegrationFeatures(
      object = object,
      assay = assay
    )
  } else if (inherits(x = object[[assay]], what = 'StdAssay')) {
    layers <- Layers(object = object, assay = assay, search = layers %||% 'data')
    scale.layer <- Layers(object = object, search = scale.layer)
    features <- features %||% VariableFeatures(
      object = object,
      assay = assay,
      nfeatures = 2000L
    )
  } else {
    abort(message = "'assay' must be a v5 or SCT assay")
  }
  if (!is.null(scale.layer)) {
    features <- intersect(
      x = features,
      y = Features(x = object, assay = assay, layer = scale.layer)
    )
  }
  if (!length(x = features)) {
    abort(message = "None of the features provided are found in this assay")
  }
  if (!is.null(orig.reduction)) {
    # Check our dimensional reduction
    orig.reduction <- orig.reduction %||% DefaultDimReduc(object = object, assay = assay)
    if (!orig.reduction %in% Reductions(object = object)) {
      abort(message = paste(sQuote(x = orig.reduction), 'is not a dimensional reduction'))
    }
    obj.orig <- object[[orig.reduction]]
    if (is.null(x = DefaultAssay(object = obj.orig))) {
      DefaultAssay(object = obj.orig) <- assay
    }
  }
  # Run the integration method
  value <- method(
    object = object[[assay]],
    assay = assay,
    orig = obj.orig,
    layers = layers,
    scale.layer = scale.layer,
    features = features,
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
# Methods for Seurat-defined generics
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Internal
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Creates data.frame with cell group assignments for integration
# uses SCT models if SCTAssay and layers otherwise
CreateIntegrationGroups <- function(object, layers, scale.layer) {
    groups <- if (inherits(x = object, what = 'SCTAssay')) {
        df <- SeuratObject::EmptyDF(n = ncol(x = object))
        row.names(x = df) <- colnames(x = object)
        for (model in levels(x = object)) {
            cc <- Cells(x = object, layer = model)
            df[cc, "group"] <- model
        }
       df
      } else if (length(x = layers) > 1L) {
          cmap <- slot(object = object, name = 'cells')[, layers]
          as.data.frame(x = labels(
              object = cmap,
              values = Cells(x = object, layer = scale.layer)
            ))
      }
    names(x = groups) <- 'group'
    return(groups)
}

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
