#' @include reexports.R
#' @include generics.R
#' @importFrom Rcpp evalCpp
#' @importFrom Matrix colSums rowSums colMeans rowMeans
#' @importFrom methods setClass setOldClass setClassUnion slot
#' slot<- setMethod new signature slotNames is setAs setValidity .hasSlot
#' @importClassesFrom Matrix dgCMatrix
#' @useDynLib Seurat
#'
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Class definitions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

setOldClass(Classes = 'package_version')

#' The AnchorSet Class
#'
#' The AnchorSet class is an intermediate data storage class that stores the anchors and other
#' related information needed for performing downstream analyses - namely data integration
#' (\code{\link{IntegrateData}}) and data transfer (\code{\link{TransferData}}).
#'
#' @slot object.list List of objects used to create anchors
#' @slot reference.cells List of cell names in the reference dataset - needed when performing data
#' transfer.
#' @slot reference.objects Position of reference object/s in object.list
#' @slot query.cells List of cell names in the query dataset - needed when performing data transfer
#' @slot anchors The anchor matrix. This contains the cell indices of both anchor pair cells, the
#' anchor score, and the index of the original dataset in the object.list for cell1 and cell2 of
#' the anchor.
#' @slot offsets The offsets used to enable cell look up in downstream functions
#' @slot weight.reduction The weight dimensional reduction used to calculate weight matrix
#' @slot anchor.features The features used when performing anchor finding.
#' @slot neighbors List containing Neighbor objects for reuse later (e.g. mapping)
#' @slot command Store log of parameters that were used
#'
#' @name AnchorSet-class
#' @rdname AnchorSet-class
#' @concept objects
#' @exportClass AnchorSet
#'
AnchorSet <- setClass(
  Class = "AnchorSet",
  contains = 'VIRTUAL',
  slots = list(
    object.list = "list",
    reference.cells = "vector",
    reference.objects = "vector",
    query.cells = "vector",
    anchors = "ANY",
    offsets = "ANY",
    weight.reduction = "DimReduc",
    anchor.features = "ANY",
    neighbors = "list",
    command = "ANY"
  )
)

#' The TransferAnchorSet Class
#'
#' Inherits from the Anchorset class. Implemented mainly for method dispatch
#' purposes.  See \code{\link{AnchorSet}} for slot details.
#'
#' @name TransferAnchorSet-class
#' @rdname TransferAnchorSet-class
#' @concept objects
#' @exportClass TransferAnchorSet
#'
TransferAnchorSet <- setClass(
  Class = "TransferAnchorSet",
  contains = "AnchorSet"
)

#' The IntegrationAnchorSet Class
#'
#' Inherits from the Anchorset class. Implemented mainly for method dispatch
#' purposes.  See \code{\link{AnchorSet}} for slot details.
#'
#' @name IntegrationAnchorSet-class
#' @rdname IntegrationAnchorSet-class
#' @concept objects
#' @exportClass IntegrationAnchorSet
#'
IntegrationAnchorSet <- setClass(
  Class = "IntegrationAnchorSet",
  contains = "AnchorSet"
)

#' The ModalityWeights Class
#'
#' The ModalityWeights class is an intermediate data storage class that stores the modality weight and other
#' related information needed for performing downstream analyses - namely data integration
#' (\code{FindModalityWeights}) and data transfer (\code{\link{FindMultiModalNeighbors}}).
#'
#' @slot modality.weight.list A list of modality weights value from all modalities
#' @slot modality.assay Names of assays for the list of dimensional reductions
#' @slot params A list of parameters used in the FindModalityWeights
#' @slot score.matrix a list of score matrices representing cross and within-modality prediction
#' score, and kernel value
#' @slot command Store log of parameters that were used
#'
#' @name ModalityWeights-class
#' @rdname ModalityWeights-class
#' @concept objects
#' @exportClass ModalityWeights
#'
ModalityWeights <- setClass(
  Class = "ModalityWeights",
  slots = list(
    modality.weight.list = "list",
    modality.assay = "vector",
    params = "list",
    score.matrix = "list",
    command = "ANY"
  )
)




#' The BridgeReferenceSet Class
#' The BridgeReferenceSet is an output from PrepareBridgeReference
#' @slot bridge The multi-omic object
#' @slot reference The Reference object only containing bridge representation assay
#' @slot params A list of parameters used in the PrepareBridgeReference
#' @slot command Store log of parameters that were used
#'
#' @name BridgeReferenceSet-class
#' @rdname BridgeReferenceSet-class
#' @concept objects
#' @exportClass BridgeReferenceSet
#'
BridgeReferenceSet <- setClass(
  Class = "BridgeReferenceSet",
  slots = list(
    bridge = "ANY",
    reference = "ANY",
    params = "list",
    command = "ANY"
  )
)

#' The IntegrationData Class
#'
#' The IntegrationData object is an intermediate storage container used internally throughout the
#' integration procedure to hold bits of data that are useful downstream.
#'
#' @slot neighbors List of neighborhood information for cells (outputs of \code{RANN::nn2})
#' @slot weights Anchor weight matrix
#' @slot integration.matrix Integration matrix
#' @slot anchors Anchor matrix
#' @slot offsets The offsets used to enable cell look up in downstream functions
#' @slot objects.ncell Number of cells in each object in the object.list
#' @slot sample.tree Sample tree used for ordering multi-dataset integration
#'
#' @name IntegrationData-class
#' @rdname IntegrationData-class
#' @concept objects
#' @exportClass IntegrationData
#'
IntegrationData <- setClass(
  Class = "IntegrationData",
  slots = list(
    neighbors = "ANY",
    weights = "ANY",
    integration.matrix = "ANY",
    anchors = "ANY",
    offsets = "ANY",
    objects.ncell = "ANY",
    sample.tree = "ANY"
  )
)

#' The SCTModel Class
#'
#' The SCTModel object is a model and parameters storage from SCTransform.
#' It can be used to calculate Pearson residuals for new genes.
#'
#' @slot feature.attributes A data.frame with feature attributes in SCTransform
#' @slot cell.attributes A data.frame with cell attributes in SCTransform
#' @slot clips A list of two numeric of length two specifying the min and max
#' values the Pearson residual will be clipped to. One for vst and one for
#' SCTransform
#' @slot umi.assay Name of the assay of the seurat object containing UMI matrix
#' and the default is RNA
#' @slot model A formula used in SCTransform
#' @slot arguments other information used in SCTransform
#' @slot median_umi Median UMI (or scale factor) used to calculate corrected counts
#'
#' @seealso \code{\link{Assay}}
#'
#' @name SCTAssay-class
#' @rdname SCTAssay-class
#' @concept objects
#'
#' @examples
#' \dontrun{
#' # SCTAssay objects are generated from SCTransform
#' pbmc_small <- SCTransform(pbmc_small)
#' }
#'
SCTModel <- setClass(
  Class = 'SCTModel',
  slots = c(
    feature.attributes = 'data.frame',
    cell.attributes = 'data.frame',
    clips = 'list',
    umi.assay = 'character',
    model = 'character',
    arguments = 'list',
    median_umi = 'numeric'
  )
)

#' The SCTAssay Class
#'
#' The SCTAssay object contains all the information found in an \code{\link{Assay}}
#' object, with extra information from the results of \code{\link{SCTransform}}
#'
#' @slot SCTModel.list A list containing SCT models
#'
#' @seealso \code{\link{Assay}}
#'
#' @name SCTAssay-class
#' @rdname SCTAssay-class
#' @concept objects
#'
#' @examples
#' # SCTAssay objects are generated from SCTransform
#' pbmc_small <- SCTransform(pbmc_small)
#' pbmc_small[["SCT"]]
#'
SCTAssay <- setClass(
  Class = 'SCTAssay',
  contains = 'Assay',
  slots = c(
    SCTModel.list = 'list'
  )
)


#' @note \code{scalefactors} objects can be created with \code{scalefactors()}
#'
#' @param spot Spot full resolution scale factor
#' @param fiducial Fiducial full resolution scale factor
#' @param hires High resolutoin scale factor
#' @param lowres Low resolution scale factor
#'
#' @rdname ScaleFactors
#' @concept objects
#' @concept spatial
#' @export
#'
scalefactors <- function(spot, fiducial, hires, lowres) {
  object <- list(
    spot = spot,
    fiducial = fiducial,
    hires = hires,
    lowres = lowres
  )
  object <- sapply(X = object, FUN = as.numeric, simplify = FALSE, USE.NAMES = TRUE)
  return(structure(.Data = object, class = 'scalefactors'))
}

setOldClass(Classes = c('scalefactors'))

#' The SlideSeq class
#'
#' The SlideSeq class represents spatial information from the Slide-seq platform
#'
#' @inheritSection SeuratObject::SpatialImage Slots
#' @slot coordinates ...
#' @concept spatial
#'
SlideSeq <- setClass(
  Class = 'SlideSeq',
  contains = 'SpatialImage',
  slots = list(
    'coordinates' = 'data.frame'
  )
)

#' The STARmap class
#'
#'
#' @inheritSection SeuratObject::SpatialImage Slots
#' @concept objects
#' @concept spatial
#'
STARmap <- setClass(
  Class = 'STARmap',
  contains = 'SpatialImage',
  slots = list(
    'coordinates' = 'data.frame',
    'qhulls' = 'data.frame'
  )
)

#' The VisiumV1 class
#'
#' The VisiumV1 class represents spatial information from the 10X Genomics Visium
#' platform
#'
#' @slot image A three-dimensional array with PNG image data, see
#' \code{\link[png]{readPNG}} for more details
#' @slot scale.factors An object of class \code{\link{scalefactors}}; see
#' \code{\link{scalefactors}} for more information
#' @slot coordinates A data frame with tissue coordinate information
#' @slot spot.radius Single numeric value giving the radius of the spots
#'
#' @name VisiumV1-class
#' @rdname VisiumV1-class
#' @concept objects
#' @concept spatial
#' @exportClass VisiumV1
#'
VisiumV1 <- setClass(
  Class = 'VisiumV1',
  contains = 'SpatialImage',
  slots = list(
    'image' = 'array',
    'scale.factors' = 'scalefactors',
    'coordinates' = 'data.frame',
    'spot.radius' = 'numeric'
  )
)

setClass(Class = 'SliceImage', contains = 'VisiumV1')

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Get a vector of cell names associated with an image (or set of images)
#'
#' @param object Seurat object
#' @param images Vector of image names
#' @param unlist Return as a single vector of cell names as opposed to a list,
#' named by image name.
#'
#' @return A vector of cell names
#'
#' @examples
#' \dontrun{
#' CellsByImage(object = object, images = "slice1")
#' }
#'
CellsByImage <- function(object, images = NULL, unlist = FALSE) {
  images <- images %||% Images(object = object)
  cells <- sapply(
    X = images,
    FUN = function(x) {
      Cells(x = object[[x]])
    },
    simplify = FALSE,
    USE.NAMES = TRUE
  )
  if (unlist) {
    cells <- unname(obj = unlist(x = cells))
  }
  return(cells)
}

#' Create a SCT Assay object
#'
#' Create a SCT object from a feature (e.g. gene) expression matrix and a list of SCTModels.
#' The expected format of the input matrix is features x cells.
#'
#' Non-unique cell or feature names are not allowed. Please make unique before
#' calling this function.
#' @param scale.data a residual matrix
#' @param SCTModel.list list of SCTModels
#' @param umi.assay The UMI assay name. Default is RNA
#' @inheritParams SeuratObject::CreateAssayObject
#'
#' @importFrom methods as
#' @importFrom Matrix colSums rowSums
#'
#' @export
#' @concept objects
#'
CreateSCTAssayObject <- function(
  counts,
  data,
  scale.data = NULL,
  umi.assay = "RNA",
  min.cells = 0,
  min.features = 0,
  SCTModel.list = NULL
) {
  assay <- CreateAssayObject(
    counts = counts,
    data = data,
    min.cells = min.cells,
    min.features = min.features
  )
  if (!is.null(scale.data)) {
    assay <- SetAssayData(object = assay, slot = "scale.data", new.data = scale.data)
  }

  slot(object = assay, name = "assay.orig") <- umi.assay

  #checking SCTModel.list format
  if (is.null(x = SCTModel.list)) {
    SCTModel.type <- "none"
    warning("An empty SCTModel will be generated due to no SCTModel input")
  } else {
    if (inherits(x = SCTModel.list, what = "SCTModel")) {
      SCTModel.list <- list(model1 = SCTModel.list)
      SCTModel.type <- "SCTModel.list"
    } else if (inherits(x = SCTModel.list, what = "list")) {
      if (inherits(x = SCTModel.list[[1]], what = "SCTModel")){
        SCTModel.type <-  "SCTModel.list"
      } else if (IsVSTout(vst.out = SCTModel.list)){
        SCTModel.type <- "vst.out"
      } else if (IsVSTout(SCTModel.list[[1]])) {
        SCTModel.type <- "vst.set"
      } else {
        stop("SCTModel input is not a correct format")
      }
    }
  }
  model.list <- switch(
    EXPR = SCTModel.type,
    "none" = {
      list()
    },
    "SCTModel.list" = {
      SCTModel.list <- lapply(X = SCTModel.list, FUN = function(model) {
        select.cell <- intersect(x = Cells(x = model), Cells(x = assay))
        if (length(x = select.cell) == 0) {
          stop("Cells in SCTModel.list don't match Cells in assay")
        } else {
          model@cell.attributes <- model@cell.attributes[select.cell, , drop = FALSE]
        }
        return(model)
      })
      SCTModel.list
    },
    "vst.out" = {
      SCTModel.list$umi.assay <- umi.assay
      SCTModel.list <- PrepVSTResults(
        vst.res = SCTModel.list,
        cell.names = Cells(x = assay)
      )
      list(model1 = SCTModel.list)
    },
    "vst.set" = {
      new.model <- lapply(
        X = SCTModel.list,
        FUN = function(vst.res) {
          vst.res$umi.assay <- umi.assay
          return(PrepVSTResults(vst.res = vst.res, cell.names = colnames(x = assay)))
        }
      )
      names(x = new.model) <- paste0("model", 1:length(x = new.model))
      new.model
    }
  )
  assay <- new(
    Class = "SCTAssay",
    assay,
    SCTModel.list = model.list
  )
  return(assay)
}

#' Slim down a Seurat object
#'
#' Keep only certain aspects of the Seurat object. Can be useful in functions
#' that utilize merge as it reduces the amount of data in the merge
#'
#' @param object A \code{\link[SeuratObject]{Seurat}} object
#' @param layers A vector or named list of layers to keep
#' @param features Only keep a subset of features, defaults to all features
#' @param assays Only keep a subset of assays specified here
#' @param dimreducs Only keep a subset of DimReducs specified here (if
#' \code{NULL}, remove all DimReducs)
#' @param graphs Only keep a subset of Graphs specified here (if \code{NULL},
#' remove all Graphs)
#' @param misc Preserve the \code{misc} slot; default is \code{TRUE}
#' @param ... Ignored
#'
#' @return \code{object} with only the sub-object specified retained
#'
#' @importFrom SeuratObject .FilterObjects .PropagateList Assays
#' Layers UpdateSlots
#'
#' @export
#'
#' @concept objects
#'
DietSeurat <- function(
  object,
  layers = NULL,
  features = NULL,
  assays = NULL,
  dimreducs = NULL,
  graphs = NULL,
  misc = TRUE,
  counts = deprecated(),
  data = deprecated(),
  scale.data = deprecated(),
  ...
) {
  CheckDots(...)
  dep.args <- c(counts = counts, data = data, scale.data = scale.data)
  for (lyr in names(x = dep.args)) {
    if (is_present(arg = dep.args[[lyr]])) {
      if (is.null(x = layers)) {
        layers <- unique(x = unlist(x = lapply(
          X = Assays(object = object),
          FUN = function(x) {
            return(Layers(object = object[[x]]))
          }
        )))
      }
      deprecate_soft(
        when = '5.0.0',
        what = paste0('DietSeurat(', lyr, ' = )'),
        with = 'DietSeurat(layers = )'
      )
      layers <- if (isTRUE(x = dep.args[[lyr]])) {
        c(layers, lyr)
      } else {
        Filter(f = \(x) x != lyr, x = layers)
      }
    }
  }
  object <- UpdateSlots(object = object)
  assays <- assays %||% Assays(object = object)
  assays <- intersect(x = assays, y = Assays(object = object))
  if (!length(x = assays)) {
    abort(message = "No assays provided were found in the Seurat object")
  }
  if (!DefaultAssay(object = object) %in% assays) {
    abort(
      message = "The default assay is slated to be removed, please change the default assay"
    )
  }
  layers <- layers %||% assays
  layers <- .PropagateList(x = layers, names = assays)
  for (assay in names(x = layers)) {
    layers[[assay]] <- tryCatch(
      expr = Layers(object = object[[assay]], search = layers[[assay]]),
      error = function(...) {
        return(character(length = 0L))
      }
    )
  }
  layers <- Filter(f = length, x = layers)
  if (!length(x = layers)) {
    abort(message = "None of the requested layers found")
  }
  for (assay in Assays(object = object)) {
    if (!(assay %in% assays)) {
      object[[assay]] <- NULL
      next
    }
    layers.rm <- setdiff(
      x = Layers(object = object[[assay]]),
      y = layers[[assay]]
    )
    if (length(x = layers.rm)) {
      if (inherits(x = object[[assay]], what = 'Assay') && all(c('counts', 'data') %in% layers.rm)) {
        abort(message = "Cannot remove both 'counts' and 'data' from v3 Assays")
      }
      for (lyr in layers.rm) {
        object <- tryCatch(expr = {
          object[[assay]][[lyr]] <- NULL
          object
        }, error = function(e) {
          if (lyr == "data"){
            object[[assay]][[lyr]] <- sparseMatrix(i = 1, j = 1, x = 1,
                         dims = dim(object[[assay]][[lyr]]), 
                         dimnames = dimnames(object[[assay]][[lyr]]))
          } else{
            slot(object = object[[assay]], name = lyr) <- new(Class = "dgCMatrix")
          }
          message("Converting layer ", lyr, " in assay ", 
                  assay, " to empty dgCMatrix")
          object
        })
      }
    }
    if (!is.null(x = features)) {
      features.assay <- intersect(
        x = features,
        y = rownames(x = object[[assay]])
      )
      if (!length(x = features.assay)) {
        warn(message = paste0(
          'No features found in assay ',
          sQuote(x = assay),
          ', removing...'
        ))
        object[[assay]] <- NULL
        next
      }
      object[[assay]] <- subset(x = object[[assay]], features = features.assay)
    }
  }
  # remove misc when desired
  if (!isTRUE(x = misc)) {
    slot(object = object, name = "misc") <- list()
  }
  # remove unspecified DimReducs and Graphs
  all.objects <- .FilterObjects(
    object = object,
    classes.keep = c('DimReduc', 'Graph')
  )
  objects.to.remove <- all.objects[!all.objects %in% c(dimreducs, graphs)]
  for (ob in objects.to.remove) {
    object[[ob]] <- NULL
  }
  cells.keep <- list()
  for (assay in  Assays(object = object)) {
    cells.keep[[assay]] <- colnames(x = object[[assay]] )
  }
  cells.keep <- intersect(colnames(x = object), unlist(cells.keep))
  if (length(cells.keep) <- ncol(x = object)) {
    object <- subset(object, cells = cells.keep)
  }
  return(object)
}

#' Filter stray beads from Slide-seq puck
#'
#' This function is useful for removing stray beads that fall outside the main
#' Slide-seq puck area. Essentially, it's a circular filter where you set a
#' center and radius defining a circle of beads to keep. If the center is not
#' set, it will be estimated from the bead coordinates (removing the 1st and
#' 99th quantile to avoid skewing the center by the stray beads). By default,
#' this function will display a \code{\link{SpatialDimPlot}} showing which cells
#' were removed for easy adjustment of the center and/or radius.
#'
#' @param object Seurat object with slide-seq data
#' @param image Name of the image where the coordinates are stored
#' @param center Vector specifying the x and y coordinates for the center of the
#' inclusion circle
#' @param radius Radius of the circle of inclusion
#' @param do.plot Display a \code{\link{SpatialDimPlot}} with the cells being
#' removed labeled.
#'
#' @return Returns a Seurat object with only the subset of cells that pass the
#' circular filter
#'
#' @concept objects
#' @concept spatial
#' @examples
#' \dontrun{
#' # This example uses the ssHippo dataset which you can download
#' # using the SeuratData package.
#' library(SeuratData)
#' data('ssHippo')
#' # perform filtering of beads
#' ssHippo.filtered <- FilterSlideSeq(ssHippo, radius = 2300)
#' # This radius looks to small so increase and repeat until satisfied
#' }
#' @export
#'
FilterSlideSeq <- function(
  object,
  image = "image",
  center = NULL,
  radius = NULL,
  do.plot = TRUE
) {
  if (!inherits(x = object[[image]], what = "SlideSeq")) {
    warning(
      "This fxn is intended for filtering SlideSeq data and is untested ",
      "outside of that context."
    )
  }
  dat <- GetTissueCoordinates(object[[image]])
  if (is.null(x = center)) {
    # heuristic for determining center of puck
    center <- c()
    x.vals <- dat[, 1]
    center[1] <- mean(
      x = x.vals[x.vals < quantile(x = x.vals, probs = 0.99) &
                   x.vals > quantile(x = x.vals, probs = 0.01)]
    )
    y.vals <- dat[, 2]
    center[2] <- mean(
      x = y.vals[y.vals < quantile(x = y.vals, probs = 0.99) &
                   y.vals > quantile(x = y.vals, probs = 0.01)]
    )
  }
  if (is.null(x = radius)) {
    stop("Please provide a radius.")
  }
  dists <- apply(X = dat, MARGIN = 1, FUN = function(x) {
    as.numeric(dist(rbind(x[c(1, 2)], center)))
  })
  cells.to.remove <- names(x = which(x = (dists > radius)))
  if (do.plot) {
    Idents(object) <- "keep"
    object <- SetIdent(object = object, cells = cells.to.remove, value = "remove")
    print(SpatialDimPlot(object = object))
  }
  return(subset(x = object, cells = cells.to.remove, invert = TRUE))
}

#' Get integration data
#'
#' @param object Seurat object
#' @param integration.name Name of integration object
#' @param slot Which slot in integration object to get
#'
#' @return Returns data from the requested slot within the integrated object
#'
#' @export
#' @concept objects
#'
GetIntegrationData <- function(object, integration.name, slot) {
  tools <- slot(object = object, name = 'tools')
  if (!(integration.name %in% names(tools))) {
    stop('Requested integration key does not exist')
  }
  int.data <- tools[[integration.name]]
  return(slot(object = int.data, name = slot))
}

#' Set integration data
#'
#' @param object Seurat object
#' @param integration.name Name of integration object
#' @param slot Which slot in integration object to set
#' @param new.data New data to insert
#'
#' @return Returns a \code{\link{Seurat}} object
#'
#' @export
#' @concept objects
#'
SetIntegrationData <- function(object, integration.name, slot, new.data) {
  tools <- slot(object = object, name = 'tools')
  if (!(integration.name %in% names(tools))) {
    new.integrated <- new(Class = 'IntegrationData')
    slot(object = new.integrated, name = slot) <- new.data
    tools[[integration.name]] <- new.integrated
    slot(object = object, name = 'tools') <- tools
    return(object)
  }
  int.data <- tools[[integration.name]]
  slot(object = int.data, name = slot) <- new.data
  tools[[integration.name]] <- int.data
  slot(object = object, name = 'tools') <- tools
  return(object)
}

#' Splits object into a list of subsetted objects.
#'
#' Splits object based on a single attribute into a list of subsetted objects,
#' one for each level of the attribute. For example, useful for taking an object
#' that contains cells from many patients, and subdividing it into
#' patient-specific objects.
#'
#' @param object Seurat object
#' @param split.by Attribute for splitting. Default is "ident". Currently
#' only supported for class-level (i.e. non-quantitative) attributes.
#'
#' @return A named list of Seurat objects, each containing a subset of cells
#' from the original object.
#'
#' @export
#' @concept objects
#'
#' @examples
#' data("pbmc_small")
#' # Assign the test object a three level attribute
#' groups <- sample(c("group1", "group2", "group3"), size = 80, replace = TRUE)
#' names(groups) <- colnames(pbmc_small)
#' pbmc_small <- AddMetaData(object = pbmc_small, metadata = groups, col.name = "group")
#' obj.list <- SplitObject(pbmc_small, split.by = "group")
#'
SplitObject <- function(object, split.by = "ident") {
  if (split.by == 'ident') {
    groupings <- Idents(object = object)
  } else {
    groupings <- FetchData(object = object, vars = split.by)[, 1]
  }
  groupings <- unique(x = as.character(x = groupings))
  obj.list <- list()
  for (i in groupings) {
    if (split.by == "ident") {
      obj.list[[i]] <- subset(x = object, idents = i)
    }
    else {
      cells <- which(x = object[[split.by, drop = TRUE]] == i)
      cells <- colnames(x = object)[cells]
      obj.list[[i]] <- subset(x = object, cells = cells)
    }
  }
  return(obj.list)
}

#' Find features with highest scores for a given dimensional reduction technique
#'
#' Return a list of features with the strongest contribution to a set of components
#'
#' @param object DimReduc object
#' @param dim Dimension to use
#' @param nfeatures Number of features to return
#' @param projected Use the projected feature loadings
#' @param balanced Return an equal number of features with both + and - scores.
#' @param ... Extra parameters passed to \code{\link{Loadings}}
#'
#' @return Returns a vector of features
#'
#' @export
#' @concept objects
#'
#' @examples
#' data("pbmc_small")
#' pbmc_small
#' TopFeatures(object = pbmc_small[["pca"]], dim = 1)
#' # After projection:
#' TopFeatures(object = pbmc_small[["pca"]], dim = 1,  projected = TRUE)
#'
TopFeatures <- function(
  object,
  dim = 1,
  nfeatures = 20,
  projected = FALSE,
  balanced = FALSE,
  ...
) {
  loadings <- Loadings(object = object, projected = projected, ...)[, dim, drop = FALSE]
  return(Top(
    data = loadings,
    num = nfeatures,
    balanced = balanced
  ))
}

#' Find cells with highest scores for a given dimensional reduction technique
#'
#' Return a list of genes with the strongest contribution to a set of components
#'
#' @param object DimReduc object
#' @param dim Dimension to use
#' @param ncells Number of cells to return
#' @param balanced Return an equal number of cells with both + and - scores.
#' @param ... Extra parameters passed to \code{\link{Embeddings}}
#'
#' @return Returns a vector of cells
#'
#' @export
#' @concept objects
#'
#' @examples
#' data("pbmc_small")
#' pbmc_small
#' head(TopCells(object = pbmc_small[["pca"]]))
#' # Can specify which dimension and how many cells to return
#' TopCells(object = pbmc_small[["pca"]], dim = 2, ncells = 5)
#'
TopCells <- function(object, dim = 1, ncells = 20, balanced = FALSE, ...) {
  embeddings <- Embeddings(object = object, ...)[, dim, drop = FALSE]
  return(Top(
    data = embeddings,
    num = ncells,
    balanced = balanced
  ))
}

#' Get nearest neighbors for given cell
#'
#' Return a vector of cell names of the nearest n cells.
#'
#' @param object \code{\link{Neighbor}} object
#' @param cell Cell of interest
#' @param n Number of neighbors to return
#'
#' @return Returns a vector of cell names
#'
#' @export
#' @concept objects
#'
TopNeighbors <- function(object, cell, n = 5) {
  indices <- Indices(object = object)[cell, 1:n]
  return(Cells(x = object)[indices])
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Methods for Seurat-defined generics
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @param assay Assay to convert
#' @param reduction Name of DimReduc to set to main reducedDim in cds
#'
#' @rdname as.CellDataSet
#' @concept objects
#' @export
#' @method as.CellDataSet Seurat
#'
as.CellDataSet.Seurat <- function(x, assay = NULL, reduction = NULL, ...) {
  CheckDots(...)
  if (!PackageCheck('monocle', error = FALSE)) {
    stop("Please install monocle from Bioconductor before converting to a CellDataSet object")
  } else if (packageVersion(pkg = 'monocle') >= package_version(x = '2.99.0')) {
    stop("Seurat can only convert to/from Monocle v2.X objects")
  }
  assay <- assay %||% DefaultAssay(object = x)
  # make variables, then run `newCellDataSet`
  # create cellData counts
  counts <- GetAssayData(object = x, assay = assay, slot = "counts")
  # metadata
  cell.metadata <- x[[]]
  feature.metadata <- x[[assay]][[]]
  if (!"gene_short_name" %in% colnames(x = feature.metadata)) {
    feature.metadata$gene_short_name <- rownames(x = feature.metadata)
  }
  pd <- new(Class = "AnnotatedDataFrame", data = cell.metadata)
  fd <- new(Class = "AnnotatedDataFrame", data = feature.metadata)
  # Now, determine the expressionFamily
  if ("monocle" %in% names(x = Misc(object = x))) {
    expressionFamily <- Misc(object = x, slot = "monocle")[["expressionFamily"]]
  } else {
    if (all(counts == floor(x = counts))) {
      expressionFamily <- VGAM::negbinomial.size()
    } else if (any(counts < 0)) {
      expressionFamily <- VGAM::uninormal()
    } else {
      expressionFamily <- VGAM::tobit()
    }
  }
  cds <- monocle::newCellDataSet(
    cellData = counts,
    phenoData = pd,
    featureData = fd,
    expressionFamily = expressionFamily
  )
  if ("monocle" %in% names(x = Misc(object = x))) {
    monocle::cellPairwiseDistances(cds = cds) <- Misc(object = x, slot = "monocle")[["cellPairwiseDistances"]]
    monocle::minSpanningTree(cds = cds) <- Misc(object = x, slot = "monocle")[["minSpanningTree"]]
    Biobase::experimentData(cds = cds) <- Misc(object = x, slot = "monocle")[["experimentData"]]
    Biobase::protocolData(cds = cds) <- Misc(object = x, slot = "monocle")[["protocolData"]]
    Biobase::classVersion(cds = cds) <- Misc(object = x, slot = "monocle")[["classVersion"]]
    # no setter methods found for following slots
    slot(object = cds, name = "lowerDetectionLimit") <- Misc(object = x, slot = "monocle")[["lowerDetectionLimit"]]
    slot(object = cds, name = "dispFitInfo") <- Misc(object = x, slot = "monocle")[["dispFitInfo"]]
    slot(object = cds, name = "auxOrderingData") <- Misc(object = x, slot = "monocle")[["auxOrderingData"]]
    slot(object = cds, name = "auxClusteringData") <- Misc(object = x, slot = "monocle")[["auxClusteringData"]]
  }
  # adding dimensionality reduction data to the CDS
  dr.slots <- c("reducedDimS", "reducedDimK", "reducedDimW", "reducedDimA")
  reduction <- reduction %||% DefaultDimReduc(object = x, assay = assay)
  if (!is.null(x = reduction)) {
    if (grepl(pattern = 'tsne', x = tolower(x = reduction))) {
      slot(object = cds, name = "dim_reduce_type") <- "tSNE"
      monocle::reducedDimA(cds = cds) <- t(x = Embeddings(object = x[[reduction]]))
    } else {
      slot(object = cds, name = "dim_reduce_type") <- reduction
      monocle::reducedDimA(cds = cds) <- Loadings(object = x[[reduction]])
      slot(object = cds, name = "reducedDimS") <- Embeddings(object = x[[reduction]])
    }
    for (ii in dr.slots) {
      if (ii %in% names(x = slot(object = x[[reduction]], name = "misc"))) {
        slot(object = cds, name = ii) <- slot(object = x[[reduction]], name = "misc")[[ii]]
      }
    }
  }
  return(cds)
}

#' Convert objects to \code{Seurat} objects
#'
#' @inheritParams SeuratObject::as.Seurat
#' @param slot Slot to store expression data as
#' @param verbose Show progress updates
#'
#' @return A \code{Seurat} object generated from \code{x}
#'
#' @importFrom utils packageVersion
#'
#' @rdname as.Seurat
#' @concept objects
#' @export
#' @method as.Seurat CellDataSet
#'
#' @seealso \code{\link[SeuratObject:as.Seurat]{SeuratObject::as.Seurat}}
#'
as.Seurat.CellDataSet <- function(
  x,
  slot = 'counts',
  assay = 'RNA',
  verbose = TRUE,
  ...
) {
  CheckDots(...)
  if (!PackageCheck('monocle', error = FALSE)) {
    stop("Please install monocle from Bioconductor before converting to a CellDataSet object")
  } else if (packageVersion(pkg = 'monocle') >= package_version(x = '2.99.0')) {
    stop("Seurat can only convert to/from Monocle v2.X objects")
  }
  slot <- match.arg(arg = slot, choices = c('counts', 'data'))
  if (verbose) {
    message("Pulling expression data")
  }
  expr <- Biobase::exprs(object = x)
  if (IsMatrixEmpty(x = expr)) {
    stop("No data provided in this CellDataSet object", call. = FALSE)
  }
  meta.data <- as.data.frame(x = Biobase::pData(object = x))
  # if cell names are NULL, fill with cell_X
  if (is.null(x = colnames(x = expr))) {
    warning(
      "The column names of the 'counts' and 'data' matrices are NULL. Setting cell names to cell_columnidx (e.g 'cell_1').",
      call. = FALSE,
      immediate. = TRUE
    )
    rownames(x = meta.data) <- colnames(x = expr) <- paste0("cell_", 1:ncol(x = expr))
  }
  # Creating the object
  if (verbose) {
    message("Building Seurat object")
  }
  if (slot == 'data') {
    assays <- list(CreateAssayObject(data = expr))
    names(x = assays) <- assay
    Key(object = assays[[assay]]) <- suppressWarnings(expr = UpdateKey(key = assay))
    object <- new(
      Class = 'Seurat',
      assays = assays,
      meta.data = meta.data,
      version = packageVersion(pkg = 'Seurat'),
      project.name = 'SeuratProject'
    )
    DefaultAssay(object = object) <- assay
  } else {
    object <- CreateSeuratObject(
      counts = expr,
      meta.data = meta.data,
      assay = assay
    )
  }
  # feature metadata
  if (verbose) {
    message("Adding feature-level metadata")
  }
  feature.metadata <- Biobase::fData(object = x)
  object[[assay]][[names(x = feature.metadata)]] <- feature.metadata
  # mean/dispersion values
  disp.table <- tryCatch(
    expr = suppressWarnings(expr = monocle::dispersionTable(cds = x)),
    error = function(...) {
      return(NULL)
    }
  )
  if (!is.null(x = disp.table)) {
    if (verbose) {
      message("Adding dispersion information")
    }
    rownames(x = disp.table) <- disp.table[, 1]
    disp.table[, 1] <- NULL
    colnames(x = disp.table) <- paste0('monocle_', colnames(x = disp.table))
    object[[assay]][[names(x = disp.table)]] <- disp.table
  } else if (verbose) {
    message("No dispersion information in CellDataSet object")
  }
  # variable features
  if ("use_for_ordering" %in% colnames(x = feature.metadata)) {
    if (verbose) {
      message("Setting variable features")
    }
    VariableFeatures(object = object, assay = assay) <- rownames(x = feature.metadata)[which(x = feature.metadata[, "use_for_ordering"])]
  } else if (verbose) {
    message("No variable features present")
  }
  # add dim reduction
  dr.name <- slot(object = x, name = "dim_reduce_type")
  if (length(x = dr.name) > 0) {
    if (verbose) {
      message("Adding ", dr.name, " dimensional reduction")
    }
    reduced.A <- t(x = slot(object = x, name = 'reducedDimA'))
    reduced.S <- t(x = slot(object = x, name = 'reducedDimS'))
    if (IsMatrixEmpty(x = reduced.S)) {
      embeddings <- reduced.A
      loadings <- new(Class = 'matrix')
    } else {
      embeddings <- reduced.S
      loadings <- t(x = reduced.A)
    }
    rownames(x = embeddings) <- colnames(x = object)
    misc.dr <- list(
      reducedDimS = slot(object = x, name = "reducedDimS"),
      reducedDimK = slot(object = x, name = "reducedDimK"),
      reducedDimW = slot(object = x, name = "reducedDimW"),
      reducedDimA = slot(object = x, name = "reducedDimA")
    )
    dr <- suppressWarnings(expr = CreateDimReducObject(
      embeddings = embeddings,
      loadings = loadings,
      assay = assay,
      key = UpdateKey(key = tolower(x = dr.name)),
      misc = misc.dr
    ))
    object[[dr.name]] <- dr
  } else if (verbose) {
    message("No dimensional reduction information found")
  }
  monocle.specific.info <- list(
    expressionFamily = slot(object = x, name = "expressionFamily"),
    lowerDetectionLimit = slot(object = x, name = "lowerDetectionLimit"),
    dispFitInfo = slot(object = x, name = "dispFitInfo"),
    cellPairwiseDistances = slot(object = x, name = "cellPairwiseDistances"),
    minSpanningTree = slot(object = x, name = "minSpanningTree"),
    auxOrderingData = slot(object = x, name = "auxOrderingData"),
    auxClusteringData = slot(object = x, name = "auxClusteringData"),
    experimentData = slot(object = x, name = "experimentData"),
    protocolData = slot(object = x, name = "protocolData"),
    classVersion = slot(object = x, name = ".__classVersion__")
  )
  Misc(object = object, slot = "monocle")  <- monocle.specific.info
  return(object)
}

#' @param counts name of the SingleCellExperiment assay to store as \code{counts};
#' set to \code{NULL} if only normalized data are present
#' @param data name of the SingleCellExperiment assay to slot as \code{data}.
#' Set to NULL if only counts are present
#' @param assay Name of assays to convert; set to \code{NULL} for all assays to be converted
#' @param project Project name for new Seurat object
#'
#' @rdname as.Seurat
#' @concept objects
#' @export
#' @method as.Seurat SingleCellExperiment
#'
as.Seurat.SingleCellExperiment <- function(
  x,
  counts = 'counts',
  data = 'logcounts',
  assay = NULL,
  project = 'SingleCellExperiment',
  ...
) {
  CheckDots(...)
  if (!PackageCheck('SingleCellExperiment', error = FALSE)) {
    stop(
      "Please install SingleCellExperiment from Bioconductor before converting to a SingeCellExperiment object.",
      "\nhttps://bioconductor.org/packages/SingleCellExperiment/",
      call. = FALSE
    )
  }
  meta.data <- as.data.frame(x = SummarizedExperiment::colData(x = x))
  if (packageVersion(pkg = "SingleCellExperiment") >= "1.14.0") {
    orig.exp <- SingleCellExperiment::mainExpName(x = x) %||% "originalexp"
  } else {
    orig.exp <- "originalexp"
  }
  if (!is.null(SingleCellExperiment::altExpNames(x = x))) {
    assayn <- assay %||% SingleCellExperiment::altExpNames(x = x)
    if (!all(assay %in% SingleCellExperiment::altExpNames(x = x))) {
      stop("One or more of the assays you are trying to convert is not in the SingleCellExperiment object")
    }
    assayn <- c(orig.exp, assayn)
  } else {
    assayn <- orig.exp
  }
  for (assay in assayn) {
    if (assay != orig.exp) {
      x <- SingleCellExperiment::swapAltExp(x = x, name = assay, saved = NULL)
    }
    # Pull matrices
    mats <- list(counts = counts, data = data)
    mats <- Filter(f = Negate(f = is.null), x = mats)
    if (length(x = mats) == 0) {
      stop("Cannot pass 'NULL' to both 'counts' and 'data'")
    }
    for (m in 1:length(x = mats)) {
      mats[[m]] <- tryCatch(
        expr = SummarizedExperiment::assay(x = x, i = mats[[m]]),
        error = function(e) {
          stop("No data in provided assay - ", mats[[m]], call. = FALSE)
        }
      )
      # if cell names are NULL, fill with cell_X
      if (is.null(x = colnames(x = mats[[m]]))) {
        warning(
          "The column names of the ",
          names(x = mats)[m],
          " matrix is NULL. Setting cell names to cell_columnidx (e.g 'cell_1').",
          call. = FALSE,
          immediate. = TRUE
        )
        cell.names <- paste0("cell_", 1:ncol(x = mats[[m]]))
        colnames(x = mats[[m]]) <- cell.names
        rownames(x = meta.data) <- cell.names
      }
    }
    assays <- if (is.null(x = mats$counts)) {
      list(CreateAssayObject(data = mats$data))
    } else if (is.null(x = mats$data)) {
      list(CreateAssayObject(counts = mats$counts))
    } else {
      a <- CreateAssayObject(counts = mats$counts)
      a <- SetAssayData(object = a, slot = 'data', new.data = mats$data)
      list(a)
    }
    names(x = assays) <- assay
    Key(object = assays[[assay]]) <- paste0(tolower(x = assay), '_')
    # Create the Seurat object
    if (!exists(x = "object")) {
      object <- CreateSeuratObject(
        counts = assays[[assay]],
        Class = 'Seurat',
        assay = assay,
        meta.data = meta.data,
        version = packageVersion(pkg = 'Seurat'),
        project.name = project
      )
    } else {
      object[[assay]] <- assays[[assay]]
    }
    DefaultAssay(object = object) <- assay
    # add feature level meta data
    md <- SingleCellExperiment::rowData(x = x)
    if (ncol(x = md) > 0) {
      # replace underscores
      rownames(x = md) <- gsub(pattern = "_", replacement = "-", x = rownames(x = md))
      md <- as.data.frame(x = md)
      # ensure order same as data
      md <- md[rownames(x = object[[assay]]), , drop = FALSE]
      object[[assay]] <- AddMetaData(
        object = object[[assay]],
        metadata = md
      )
    }
    Idents(object = object) <- project
    # Get DimReduc information, add underscores if needed and pull from different alt EXP
    if (length(x = SingleCellExperiment::reducedDimNames(x = x)) > 0) {
      for (dr in SingleCellExperiment::reducedDimNames(x = x)) {
        embeddings <- as.matrix(x = SingleCellExperiment::reducedDim(x = x, type = dr))
        if (is.null(x = rownames(x = embeddings))) {
          rownames(x = embeddings)  <- cell.names
        }
        if (isTRUE(x = !grepl('_$',
        gsub(pattern = "[[:digit:]]",
          replacement = "_",
          x = colnames(x = SingleCellExperiment::reducedDim(x = x, type = dr))[1]
        )))) {
        key <- gsub(
          pattern = "[[:digit:]]",
          replacement = "_",
          x = colnames(x = SingleCellExperiment::reducedDim(x = x, type = dr))[1]
        )
        } else
        {
          key <- gsub(
            pattern = "[[:digit:]]",
            replacement = "",
            x = colnames(x = SingleCellExperiment::reducedDim(x = x, type = dr))[1]
          )
       }
        if (length(x = key) == 0) {
          key <- paste0(dr, "_")
        }
        colnames(x = embeddings) <- paste0(key, 1:ncol(x = embeddings))
        object[[dr]] <- CreateDimReducObject(
          embeddings = embeddings,
          key = key,
          assay = DefaultAssay(object = object)
        )
      }
    }
  }
  return(object)
}

#' @param assay Assays to convert
#'
#' @rdname as.SingleCellExperiment
#' @concept objects
#' @export
#' @method as.SingleCellExperiment Seurat
#'
as.SingleCellExperiment.Seurat <- function(x, assay = NULL, ...) {
  CheckDots(...)
  if (!PackageCheck('SingleCellExperiment', error = FALSE)) {
    stop("Please install SingleCellExperiment from Bioconductor before converting to a SingeCellExperiment object")
  }
  assay <- assay %||% Assays(object = x)
  if (!all(assay %in% Assays(object = x))) {
    stop("One or more of the assays you are trying to convert is not in the Seurat object")
  }
  if (DefaultAssay(object = x) %in% assay) {
    assay <- union(DefaultAssay(object = x), assay)
  }
  experiments <- list()
  for (assayn in assay) {
    assays <- list(
      counts = GetAssayData(object = x, assay = assayn, slot = "counts"),
      logcounts = GetAssayData(object = x, assay = assayn, slot = "data")
    )
    scaledata_a <- GetAssayData(object = x, assay = assayn, slot = "scale.data")
    if (isTRUE(x = all.equal(
      target = dim(x = assays[["counts"]]),
      current = dim(x = scaledata_a))
    )) {
      assays[["scaledata"]] <- scaledata_a
    }
    assays <- assays[sapply(X = assays, FUN = nrow) != 0]
    sume <- SummarizedExperiment::SummarizedExperiment(assays = assays)
    experiments[[assayn]] <- sume
  }
  # create one single cell experiment
  sce <- as(object = experiments[[1]], Class = "SingleCellExperiment")
  orig.exp.name <- names(x = experiments[1])
  if (packageVersion(pkg = "SingleCellExperiment") >= "1.14.0") {
    SingleCellExperiment::mainExpName(sce) <- names(x = experiments[1])
  }
  if (length(x = experiments) > 1) {
    sce <- SingleCellExperiment::SingleCellExperiment(sce, altExps = experiments)
    sce <- SingleCellExperiment::swapAltExp(
      x = sce,
      name = orig.exp.name,
      saved = NULL
    )
  }
  metadata <- x[[]]
  metadata$ident <- Idents(object = x)
  SummarizedExperiment::colData(x = sce) <- S4Vectors::DataFrame(metadata)
  for (assayn in assay) {
    if (assayn != orig.exp.name) {
      sce <- SingleCellExperiment::swapAltExp(
        x = sce,
        name = assayn,
        saved = orig.exp.name
      )
      SummarizedExperiment::rowData(x = sce) <- S4Vectors::DataFrame(x[[assayn]][[]])
      sce <- SingleCellExperiment::swapAltExp(
        x = sce,
        name = orig.exp.name,
        saved = assayn
      )
    }
  }
  for (dr in FilterObjects(object = x, classes.keep = "DimReduc")) {
    assay.used <- DefaultAssay(object = x[[dr]])
    swap.exp <- assay.used %in% SingleCellExperiment::altExpNames(x = sce) & assay.used != orig.exp.name
    if (swap.exp) {
      sce <- SingleCellExperiment::swapAltExp(
        x = sce,
        name = assay.used,
        saved = orig.exp.name
      )
    }
    SingleCellExperiment::reducedDim(x = sce, type = toupper(x = dr)) <- Embeddings(object = x[[dr]])
    if (swap.exp) {
      sce <- SingleCellExperiment::swapAltExp(
        x = sce,
        name = orig.exp.name,
        saved = assay.used
      )
    }
  }
  return(sce)
}

#' Cast to Sparse
#'
#' @inheritParams SeuratObject::as.sparse
#'
#' @importFrom methods is
#' @importFrom Matrix sparseMatrix
#'
#' @rdname as.sparse
#' @concept objects
#' @export
#' @method as.sparse H5Group
#'
#'
#' @seealso \code{\link[SeuratObject:as.sparse]{SeuratObject::as.sparse}}
#'
as.sparse.H5Group <- function(x, ...) {
  CheckDots(...)
  for (i in c('data', 'indices', 'indptr')) {
    if (!x$exists(name = i) || !is(object = x[[i]], class2 = 'H5D')) {
      stop("Invalid H5Group specification for a sparse matrix, missing dataset ", i)
    }
  }
  if ('h5sparse_shape' %in% hdf5r::h5attr_names(x = x)) {
    return(sparseMatrix(
      i = x[['indices']][] + 1,
      p = x[['indptr']][],
      x = x[['data']][],
      dims = rev(x = hdf5r::h5attr(x = x, which = 'h5sparse_shape'))
    ))
  }
  return(sparseMatrix(
    i = x[['indices']][] + 1,
    p = x[['indptr']][],
    x = x[['data']][]
  ))
}


#' @method as.sparse IterableMatrix
#' @export
#'
as.sparse.IterableMatrix <- function(x, ...) {
  return(as(object = x, Class = 'dgCMatrix'))
}


#' Get Cell Names
#'
#' @inheritParams SeuratObject::Cells
#'
#' @rdname Cells
#' @concept objects
#' @method Cells SCTModel
#' @export
#'
Cells.SCTModel <- function(x, ...) {
  return(rownames(x = slot(object = x, name = "cell.attributes")))
}

#' @method Cells SCTAssay
#' @export
#'
Cells.SCTAssay <- function(x, layer = NA) {
  layer <- layer %||% levels(x = x)[1L]
  if (rlang::is_na(x = layer)) {
    return(colnames(x = x))
  }
  return(Cells(x = components(object = x, model = layer)))
}

#' @rdname Cells
#' @concept objects
#' @concept spatial
#' @method Cells SlideSeq
#' @export
#'
#' @seealso \code{\link[SeuratObject:Cells]{SeuratObject::Cells}}
#'
Cells.SlideSeq <- function(x, ...) {
  return(rownames(x = GetTissueCoordinates(object = x)))
}

#' @rdname Cells
#' @concept objects
#' @concept spatial
#' @method Cells STARmap
#' @export
#'
Cells.STARmap <- function(x, ...) {
  return(rownames(x = GetTissueCoordinates(object = x)))
}

#' @rdname Cells
#' @concept objects
#' @method Cells VisiumV1
#' @export
#'
Cells.VisiumV1 <- function(x, ...) {
  return(rownames(x = GetTissueCoordinates(object = x, scale = NULL)))
}

#' @importFrom SeuratObject DefaultLayer Layers
#'
#' @method Features SCTAssay
#' @export
#'
Features.SCTAssay <- function(x, layer = NA) {
  layer <- layer %||% DefaultLayer(object = x)
  if (rlang::is_na(x = layer)) {
    return(rownames(x = x))
  }
  layer <- rlang::arg_match(arg = layer, values = c(Layers(object = x), levels(x = x)))
  if (layer %in% levels(x = x)) {
    return(Features(x = components(object = x, model = layer)))
  }
  return(NextMethod())
}

#' @method Features SCTModel
#' @export
#'
Features.SCTModel <- function(x, ...) {
  return(rownames(x = SCTResults(object = x, slot = 'feature.attributes')))
}

#' @param assay Assay to get
#'
#' @rdname GetAssay
#' @concept objects
#' @export
#' @method GetAssay Seurat
#'
#' @examples
#' data("pbmc_small")
#' GetAssay(object = pbmc_small, assay = "RNA")
#'
GetAssay.Seurat <- function(object, assay = NULL, ...) {
  CheckDots(...)
  assay <- assay %||% DefaultAssay(object = object)
  object.assays <- FilterObjects(object = object, classes.keep = 'Assay')
  if (!assay %in% object.assays) {
    stop(paste0(
      assay,
      " is not an assay present in the given object. Available assays are: ",
      paste(object.assays, collapse = ", ")
    ))
  }
  return(slot(object = object, name = 'assays')[[assay]])
}

#' Get Image Data
#'
#' @inheritParams SeuratObject::GetImage
#'
#' @rdname GetImage
#' @method GetImage SlideSeq
#' @concept objects
#' @concept spatial
#' @export
#'
#' @seealso \code{\link[SeuratObject:GetImage]{SeuratObject::GetImage}}
#'
GetImage.SlideSeq <- function(
  object,
  mode = c('grob', 'raster', 'plotly', 'raw'),
  ...
) {
  mode <- match.arg(arg = mode)
  return(NullImage(mode = mode))
}

#' @rdname GetImage
#' @method GetImage STARmap
#' @concept objects
#' @concept spatial
#' @export
#'
GetImage.STARmap <- function(
  object,
  mode = c('grob', 'raster', 'plotly', 'raw'),
  ...
) {
  mode <- match.arg(arg = mode)
  return(NullImage(mode = mode))
}

#' @importFrom plotly raster2uri
#' @importFrom grDevices as.raster
#' @importFrom grid rasterGrob unit
#'
#' @rdname GetImage
#' @concept objects
#' @concept spatial
#' @method GetImage VisiumV1
#' @export
#'
GetImage.VisiumV1 <- function(
  object,
  mode = c('grob', 'raster', 'plotly', 'raw'),
  ...
) {
  mode <- match.arg(arg = mode)
  image <- slot(object = object, name = 'image')
  image <- switch(
    EXPR = mode,
    'grob' = rasterGrob(
      image = image,
      width = unit(x = 1, units = 'npc'),
      height = unit(x = 1, units = 'npc')
    ),
    'raster' = as.raster(x = image),
    'plotly' = list(
      source = raster2uri(r = GetImage(object = object, mode = 'raster')),
      xref = 'x',
      yref = 'y',
      # x = -7,
      # y = -7,
      sizex = ncol(x = object),
      sizey = nrow(x = object),
      sizing = 'stretch',
      opacity = 1,
      layer = 'below'
    ),
    'raw' = image,
    stop("Unknown image mode: ", mode, call. = FALSE)
  )
  return(image)
}

#' Get Tissue Coordinates
#'
#' @inheritParams SeuratObject::GetTissueCoordinates
#'
#' @rdname GetTissueCoordinates
#' @method GetTissueCoordinates SlideSeq
#' @concept objects
#' @concept spatial
#' @export
#'
#' @seealso \code{\link[SeuratObject:GetTissueCoordinates]{SeuratObject::GetTissueCoordinates}}
#'
GetTissueCoordinates.SlideSeq <- function(object, ...) {
  coords <- slot(object = object, name = 'coordinates')
  colnames(x = coords) <- c('x', 'y')
  # coords$y <- -rev(x = coords$y) + 1
  # coords$y <- FlipCoords(x = coords$y)
  coords$cells <- rownames(x = coords)
  return(coords)
}

#' @param qhulls return qhulls instead of centroids
#'
#' @rdname GetTissueCoordinates
#' @method GetTissueCoordinates STARmap
#' @concept objects
#' @concept spatial
#' @export
#'
GetTissueCoordinates.STARmap <- function(object, qhulls = FALSE, ...) {
  if (qhulls) {
    return(slot(object = object, name = 'qhulls'))
  }
  return(slot(object = object, name = 'coordinates'))
}

#' @param scale A factor to scale the coordinates by; choose from: 'tissue',
#' 'fiducial', 'hires', 'lowres', or \code{NULL} for no scaling
#' @param cols Columns of tissue coordinates data.frame to pull
#'
#' @rdname GetTissueCoordinates
#' @method GetTissueCoordinates VisiumV1
#' @concept objects
#' @concept spatial
#' @export
#'
GetTissueCoordinates.VisiumV1 <- function(
  object,
  scale = 'lowres',
  cols = c('imagerow', 'imagecol'),
  ...
) {
  cols <- cols %||% colnames(x = slot(object = object, name = 'coordinates'))
  if (!is.null(x = scale)) {
    coordinates <- slot(object = object, name = 'coordinates')[, c('imagerow', 'imagecol')]
    scale <- match.arg(arg = scale, choices = c('spot', 'fiducial', 'hires', 'lowres'))
    scale.use <- ScaleFactors(object = object)[[scale]]
    coordinates <- coordinates * scale.use
  } else {
    coordinates <- slot(object = object, name = 'coordinates')[, cols]
  }
  return(coordinates)
}

#' Get Variable Feature Information
#'
#' Get variable feature information from \code{\link{SCTAssay}} objects
#'
#' @inheritParams SeuratObject::HVFInfo
#'
#' @export
#' @method HVFInfo SCTAssay
#'
#' @seealso \code{\link[SeuratObject]{HVFInfo}}
#'
#' @examples
#' # Get the HVF info directly from an SCTAssay object
#' pbmc_small <- SCTransform(pbmc_small)
#' HVFInfo(pbmc_small[["SCT"]], selection.method = 'sct')[1:5, ]
#'
HVFInfo.SCTAssay <- function(object, selection.method, status = FALSE, ...) {
  CheckDots(...)
  disp.methods <- c('mean.var.plot', 'dispersion', 'disp')
  if (tolower(x = selection.method) %in% disp.methods) {
    selection.method <- 'mvp'
  }
  selection.method <- switch(
    EXPR = tolower(x = selection.method),
    'sctransform' = 'sct',
    selection.method
  )
  vars <- c('gmean', 'variance', 'residual_variance')
  hvf.info <- SCTResults(object = object, slot = "feature.attributes")[,vars]
  if (status) {
    hvf.info$variable <- FALSE
    hvf.info[VariableFeatures(object = object), "variable"] <-  TRUE
  }
  return(hvf.info)
}

#' Get Spot Radius
#'
#' @inheritParams SeuratObject::Radius
#'
#' @rdname Radius
#' @concept objects
#' @concept spatial
#' @method Radius SlideSeq
#' @export
#'
#' @seealso \code{\link[SeuratObject:Radius]{SeuratObject::Radius}}
#'
Radius.SlideSeq <- function(object) {
  return(0.005)
}

#' @rdname Radius
#' @concept objects
#' @concept spatial
#' @method Radius STARmap
#' @export
#'
Radius.STARmap <- function(object) {
  return(NULL)
}

#' @rdname Radius
#' @concept objects
#' @concept spatial
#' @method Radius VisiumV1
#' @export
#'
Radius.VisiumV1 <- function(object) {
  return(slot(object = object, name = 'spot.radius'))
}

#' @rdname RenameCells
#' @export
#' @concept objects
#' @method RenameCells SCTAssay
#'
RenameCells.SCTAssay <- function(object, new.names = NULL, ...) {
  CheckDots(...)
  old.names <- Cells(x = object)
  names(x = new.names) <- old.names
  cell.attributes <- SCTResults(object = object, slot = "cell.attributes")
  if (length(x = cell.attributes) > 0) {
    if (is.data.frame(x = cell.attributes)) {
      old.names <- rownames(x = cell.attributes)
      rownames(x = cell.attributes) <- unname(obj = new.names[old.names])
    } else {
      cell.attributes <- lapply(
        X = cell.attributes,
        FUN = function(x) {
          old.names <- rownames(x = x)
          rownames(x = x) <- unname(obj = new.names[old.names])
          return(x)
        }
      )
    }
    SCTResults(object = object, slot = "cell.attributes") <- cell.attributes
  }
  new.names <- unname(obj = new.names)
  object <- NextMethod()
  return(object)
}

#' Rename Cells in an Object
#'
#' @inheritParams SeuratObject::RenameCells
#'
#' @rdname RenameCells
#' @concept objects
#' @method RenameCells SlideSeq
#' @export
#'
#' @seealso \code{\link[SeuratObject:RenameCells]{SeuratObject::RenameCells}}
#'
RenameCells.SlideSeq <- function(object, new.names = NULL, ...) {
  return(RenameCells.VisiumV1(object = object, new.names = new.names))
}

#' @rdname RenameCells
#' @concept objects
#' @method RenameCells STARmap
#' @export
#'
RenameCells.STARmap <- function(object, new.names = NULL, ...) {
  names(x = new.names) <- Cells(x = object)
  object <- RenameCells.VisiumV1(object = object, new.names = new.names)
  qhulls <- GetTissueCoordinates(object = object, qhull = TRUE)
  qhulls$cell <- new.names[qhulls$cell]
  slot(object = object, name = "qhulls") <- qhulls
  return(object)
}

#' @rdname RenameCells
#' @concept objects
#' @method RenameCells VisiumV1
#' @export
#'
RenameCells.VisiumV1 <- function(object, new.names = NULL, ...) {
  if (is.null(x = new.names)) {
    return(object)
  } else if (length(x = new.names) != length(x = Cells(x = object))) {
    stop("Wrong number of cell/spot names", call. = FALSE)
  }
  names(x = new.names) <- Cells(x = object)
  coordinates <- GetTissueCoordinates(object = object, scale = NULL, cols = NULL)
  rownames(x = coordinates) <- new.names[rownames(x = coordinates)]
  slot(object = object, name = 'coordinates') <- coordinates
  return(object)
}

#' @rdname SCTResults
#' @export
#' @method SCTResults SCTModel
#'
SCTResults.SCTModel <- function(object, slot, ...) {
  CheckDots(...)
  slots.use <- c('feature.attributes', 'cell.attributes', 'clips','umi.assay',  'model', 'arguments', 'median_umi')
  if (!slot %in% slots.use) {
    stop(
      "'slot' must be one of ",
      paste(slots.use, collapse = ', '),
      call. = FALSE
    )
  }
  return(slot(object = object, name = slot))
}

#' @rdname SCTResults
#' @concept objects
#' @export
#' @method SCTResults<- SCTModel
#'
"SCTResults<-.SCTModel" <- function(object, slot, ..., value) {
  slots.use <- c('feature.attributes', 'cell.attributes', 'clips','umi.assay', 'model', 'arguments', 'median_umi')
  if (!slot %in% slots.use) {
    stop(
      "'slot' must be one of ",
      paste(slots.use, collapse = ', '),
      call. = FALSE
    )
  }
  slot(object = object, name = slot) <- value
  return(object)
}

#' @param slot Which slot to pull the SCT results from
#' @param model Name of SCModel to pull result from. Available names can be
#' retrieved with \code{levels}.
#'
#' @return Returns the value present in the requested slot for the requested
#' group. If group is not specified, returns a list of slot results for each
#' group unless there is only one group present (in which case it just returns
#' the slot directly).
#'
#' @rdname SCTResults
#' @concept objects
#' @export
#' @method SCTResults SCTAssay
#'
SCTResults.SCTAssay <- function(object, slot, model = NULL, ...) {
  CheckDots(...)
  slots.use <- c('feature.attributes', 'cell.attributes', 'clips', 'umi.assay',  'model', 'arguments', 'median_umi')
  if (!slot %in% slots.use) {
    stop(
      "'slot' must be one of ",
      paste(slots.use, collapse = ', '),
      call. = FALSE
    )
  }
  model <- model %||% levels(x = object)
  model.list <- slot(object = object, name = "SCTModel.list")[model]
  results.list <- lapply(X = model.list, FUN = function(x) SCTResults(object = x, slot = slot))
  if (length(x = results.list) == 1) {
    results.list <- results.list[[1]]
  }
  return(results.list)
}

#' @rdname SCTResults
#' @concept objects
#' @export
#' @method SCTResults<- SCTAssay
#'
"SCTResults<-.SCTAssay" <- function(object, slot, model = NULL, ..., value) {
  slots.use <- c('feature.attributes', 'cell.attributes', 'clips','umi.assay', 'model', 'arguments', 'median_umi')
  if (!slot %in% slots.use) {
    stop(
      "'slot' must be one of ",
      paste(slots.use, collapse = ', '),
      call. = FALSE
    )
  }
  model <- model %||% levels(x = object)
  model.list <- slot(object = object, name = "SCTModel.list")[model]
  if (!is.list(x = value) | is.data.frame(x = value)) {
    value <- list(value)
  }
  model.names <- names(x = model.list)
  model.list <- lapply(
    X = 1:length(x = model.list),
    FUN = function(x) {
      SCTResults(object = model.list[[x]], slot = slot) <- value[[x]]
      return(model.list[[x]])
    }
  )
  names(x = model.list) <- model.names
  slot(object = object, name = "SCTModel.list")[model.names] <- model.list
  return(object)
}

#' @param assay Assay in the Seurat object to pull from
#'
#' @rdname SCTResults
#' @export
#' @concept objects
#' @method SCTResults Seurat
#'
SCTResults.Seurat <- function(object, assay = "SCT", slot, model = NULL, ...) {
  CheckDots(...)
  return(SCTResults(object = object[[assay]], slot = slot, model = model, ...))
}

#' @importFrom utils head
#' @method VariableFeatures SCTModel
#' @export
#'
VariableFeatures.SCTModel <- function(object, nfeatures = 3000, ...) {
  if (!is_scalar_integerish(x = nfeatures) || (!is_na(x = nfeatures < 1L) && nfeatures < 1L)) {
    abort(message = "'nfeatures' must be a single positive integer")
  }
  feature.attr <- SCTResults(object = object, slot = 'feature.attributes')
  feature.variance <- feature.attr[, 'residual_variance']
  names(x = feature.variance) <- row.names(x = feature.attr)
  feature.variance <- sort(x = feature.variance, decreasing = TRUE)
  if (is_na(x = nfeatures)) {
    return(names(x = feature.variance))
  }
  return(head(x = names(x = feature.variance), n = nfeatures))
}

#' @importFrom utils head
#' @method VariableFeatures SCTAssay
#' @export
#'
VariableFeatures.SCTAssay <- function(
  object,
  layer = NULL,
  nfeatures = NULL,
  simplify = TRUE,
  use.var.features = TRUE,
  ...
) {
  # Is the information already in var.features?
  var.features.existing <- slot(object = object, name = "var.features")
  nfeatures <- nfeatures %||% length(x = var.features.existing) %||% 3000
  if (is.null(x = layer)) {
    layer <- levels(x = object)
  }
  if (simplify == TRUE & use.var.features == TRUE & length(var.features.existing)>=nfeatures){
     return (head(x = var.features.existing, n = nfeatures))
  }

  layer <- match.arg(arg = layer, choices = levels(x = object), several.ok = TRUE)
  # run variable features on each model

  vf.list <- sapply(
    X = layer,
    FUN = function(lyr) {
      return(VariableFeatures(
        object = components(object = object, model = lyr),
        nfeatures = nfeatures,
        ...
      ))
    },
    simplify = FALSE,
    USE.NAMES = TRUE
  )
  if (isFALSE(x = simplify)){
    return (vf.list)
  }
  var.features <- sort(
    x = table(unlist(x = vf.list, use.names = FALSE)),
    decreasing = TRUE
  )
  if (length(x = var.features) == 0) {
    return(NULL)
  }
  for (i in 1:length(x = layer)) {
    vst_out <- SCTModel_to_vst(SCTModel = slot(object = object, name = "SCTModel.list")[[layer[[i]]]])
    var.features <- var.features[names(x = var.features) %in% rownames(x = vst_out$gene_attr)]
  }
  tie.val <- var.features[min(nfeatures, length(x = var.features))]
  features <- names(x = var.features[which(x = var.features > tie.val)])
  if (length(x = features) > 0) {
    feature.ranks <- sapply(X = features, FUN = function(x) {
      ranks <- sapply(X = vf.list, FUN = function(vf) {
        if (x %in% vf) {
          return(which(x = x == vf))
        }
        return(NULL)
      })
      median(x = unlist(x = ranks))
    })
    features <- names(x = sort(x = feature.ranks))
  }
  features.tie <- var.features[which(x = var.features == tie.val)]
  tie.ranks <- sapply(X = names(x = features.tie), FUN = function(x) {
    ranks <- sapply(X = vf.list, FUN = function(vf) {
      if (x %in% vf) {
        return(which(x = x == vf))
      }
      return(NULL)
    })
    median(x = unlist(x = ranks))
  })
  features <- c(
    features,
    names(x = head(x = sort(x = tie.ranks), nfeatures - length(x = features)))
  )
  return(features)
}

#' @rdname ScaleFactors
#' @method ScaleFactors VisiumV1
#' @export
#' @concept spatial
#'
ScaleFactors.VisiumV1 <- function(object, ...) {
  return(slot(object = object, name = 'scale.factors'))
}

#' @rdname ScaleFactors
#' @method ScaleFactors VisiumV1
#' @export
#' @concept spatial
#'
ScaleFactors.VisiumV1 <- function(object, ...) {
  return(slot(object = object, name = 'scale.factors'))
}

#' @rdname FetchData
#' @method FetchData VisiumV1
#' @export
#' @concept spatial
#'
FetchData.VisiumV1 <- function(
  object,
  vars,
  cells = NULL,
  ...
) {
  if (is.numeric(x = cells)) {
    cells <- Cells(x = object)[cells]
  } else if (is.null(x = cells)) {
    cells <- Cells(x = object)
  }
  vars.unkeyed <- gsub(pattern = paste0('^', Key(object)), replacement = '', x = vars)
  coords <- GetTissueCoordinates(object = object)[cells, vars.unkeyed, drop = FALSE]
  colnames(x = coords) <- vars
  return(coords)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Methods for R-defined generics
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @method [ SlideSeq
#' @concept objects
#' @export
#'
"[.SlideSeq" <- function(x, i, ...) {
  return(subset(x = x, cells = i, ...))
}

#' @method [ VisiumV1
#' @export
#'
"[.VisiumV1" <- function(x, i, ...) {
  return(subset(x = x, cells = i))
}

#' @method components SCTAssay
#' @export
#'
components.SCTAssay <- function(object, model) {
  model <- rlang::arg_match(arg = model, values = levels(x = object))
  return(slot(object = object, name = 'SCTModel.list')[[model]])
}

#' @method dim SlideSeq
#' @concept objects
#' @export
#'
dim.SlideSeq <- function(x) {
  # return(dim(x = GetImage(object = x, mode = 'raw')))
  return(c(599, 600))
}

#' @method dim STARmap
#' @concept objects
#' @export
#'
dim.STARmap <- function(x) {
  coords <- GetTissueCoordinates(object = x)
  return(c(
    max(coords[, 1]) - min(coords[, 1]),
    max(coords[, 2]) - min(coords[, 2])
  ))
}

#' @method dim VisiumV1
#' @concept objects
#' @export
#'
dim.VisiumV1 <- function(x) {
  return(dim(x = GetImage(object = x)$raster))
}

#' @rdname SCTAssay-class
#' @name SCTAssay-class
#'
#' @section Get and set SCT model names:
#' SCT results are named by initial run of \code{\link{SCTransform}} in order
#' to keep SCT parameters straight between runs. When working with merged
#' \code{SCTAssay} objects, these model names are important. \code{levels}
#' allows querying the models present. \code{levels<-} allows the changing of
#' the names of the models present, useful when merging \code{SCTAssay} objects.
#' Note: unlike normal \code{\link[base]{levels<-}}, \code{levels<-.SCTAssay}
#' allows complete changing of model names, not reordering.
#'
#' @param x An \code{SCTAssay} object
#'
#' @return \code{levels}: SCT model names
#'
#' @export
#' @concept objects
#' @method levels SCTAssay
#'
#' @examples
#' \dontrun{
#' # Query and change SCT model names
#' levels(pbmc_small[['SCT']])
#' levels(pbmc_small[['SCT']]) <- '3'
#' levels(pbmc_small[['SCT']])
#' }
#'
levels.SCTAssay <- function(x) {
  return(names(x = slot(object = x, name = "SCTModel.list")))
}

#' @rdname SCTAssay-class
#' @name SCTAssay-class
#'
#' @param value New levels, must be in the same order as the levels present
#'
#' @return \code{levels<-}: \code{x} with updated SCT model names
#'
#' @export
#' @concept objects
#' @method levels<- SCTAssay
#'
"levels<-.SCTAssay" <- function(x, value) {
  value <- sapply(X = value, FUN = function(v) {
    if (suppressWarnings(expr = !is.na(x = as.numeric(x = v)))) {
      warning("SCTModel groups cannot be number, group is added in front of ", v)
      v <- paste0("group", v)
    }
    return (v)
  })
  # Get current levels
  levels <- levels(x = x)
  if (length(x = value) != length(x = levels)) {
    stop("Must provide a vector of length ", length(x = levels), " as new levels.", call. = FALSE)
  }
  names(x = slot(object = x, name = "SCTModel.list")) <- value
  return(x)
}

#' Merge SCTAssay objects
#'
#' @inheritParams SeuratObject::merge
#' @param x A \code{\link[SeuratObject]{Seurat}} object
#' @param na.rm If na.rm = TRUE, this will only preserve residuals that are
#' present in all SCTAssays being merged. Otherwise, missing residuals will be
#' populated with NAs.
#' @export
#' @method merge SCTAssay
#' @concept objects
#'
merge.SCTAssay <- function(
  x = NULL,
  y = NULL,
  add.cell.ids = NULL,
  merge.data = TRUE,
  na.rm = TRUE,
  ...
) {
  assays <- c(x, y)
  if (any(sapply(
    X = assays,
    FUN = function(assay.i) inherits(x = assay.i, what = "Assay5")
  ))) {
    return(merge(x = as(x, "Assay5"), y, ...))
  }
  parent.call <- grep(pattern = "merge.Seurat", x = sys.calls())
  if (length(x = parent.call) > 0) {
    # Try and fill in missing residuals if called in the context of merge.Seurat
    all.features <- unique(
      x = unlist(
        x = lapply(
          X = assays,
          FUN = function(assay) {
      if (inherits(x = assay, what = "SCTAssay")) {
        return(rownames(x = GetAssayData(object = assay, slot = "scale.data")))
      }
    })
    )
    )
    if (!is.null(all.features)) {
      assays <- lapply(X = 1:length(x = assays), FUN = function(assay) {
        if (inherits(x = assays[[assay]], what = "SCTAssay")) {
          parent.environ <- sys.frame(which = parent.call[1])
          seurat.object <- parent.environ$objects[[assay]]
          seurat.object <- suppressWarnings(expr = GetResidual(object = seurat.object, features = all.features,
                                                               assay = parent.environ$assay, verbose = FALSE))
          return(seurat.object[[parent.environ$assay]])
        }
        return(assays[[assay]])
      })
    }
  }
  sct.check <- sapply(X = assays, FUN = function(x) inherits(x = x, what = "SCTAssay"))
  if (any(!sct.check)) {
    warning("Attempting to merge an SCTAssay with another Assay type \n",
            "Converting all to standard Assay objects.", call. = FALSE)
    assays <- lapply(1:length(x = assays), FUN = function(x) {
      if (sct.check[x]) {
        assays[[x]] <- as(object = assays[[x]], Class = "Assay")
      }
      return(assays[[x]])
    })
    combined.assay <- merge(
        x = assays[[1]],
        y = assays[2:length(x = assays)],
        add.cell.ids = add.cell.ids,
        merge.data = merge.data
    )
    return(combined.assay)
  }
  combined.assay <- NextMethod()
  all.levels <- unlist(x = lapply(X = assays, FUN = levels))
  while (anyDuplicated(x = all.levels)) {
    levels.duplicate <- which(x = duplicated(x = all.levels))
    all.levels <- sapply(X = 1:length(x = all.levels), FUN = function(l) {
      if (l %in% levels.duplicate) {
        return(tryCatch(
          expr = as.numeric(x = all.levels[l]) + 1,
          warning = function(...) {
            make.unique(names = all.levels)[l]
          },
          error = function(...){
            make.unique(names = all.levels)[l]
          }
        ))
      } else {
        return(all.levels[l])
      }
    })
  }
  scale.data <- lapply(X = assays, FUN = function(x) {
    dat <- GetAssayData(object = x, slot = "scale.data")
    if (ncol(x = dat) == 0) {
      dat <- matrix(ncol = ncol(x = x))
    }
    return(dat)
  })
  all.features <- lapply(X = scale.data, FUN = rownames)
  if (na.rm) {
    # merge intersection of possible residuals
    scaled.features <- names(x = which(x = table(x = unlist(x = all.features)) == length(x = assays)))
    if (length(x = scaled.features) == 0) {
      scale.data <- list(new(Class = "matrix"))
    } else {
      scale.data <- lapply(X = scale.data, FUN = function(x) x[scaled.features, ])
    }
  } else {
    scaled.features <- unique(x = unlist(x = all.features))
    scale.data <- lapply(X = 1:length(x = scale.data), FUN = function(x) {
      na.features <- setdiff(x = scaled.features, y = rownames(x = scale.data[[x]]))
      na.mat <- matrix(
        data = NA,
        nrow = length(x = na.features),
        ncol = ncol(x = assays[[x]]),
        dimnames = list(na.features, colnames(x = assays[[x]]))
      )
      return(rbind(scale.data[[x]], na.mat)[scaled.features, ])
    })
  }
  scale.data <- do.call(what = cbind, args = scale.data)
  combined.assay <- SetAssayData(object = combined.assay, slot = "scale.data", new.data = scale.data)
  model.list <- unlist(x = lapply(
    X = assays,
    FUN = slot,
    name = "SCTModel.list"
  ))
  names(x = model.list) <- all.levels
  model.list <- model.list %||% list()
  combined.assay <- new(
    Class = "SCTAssay",
    combined.assay,
    SCTModel.list = model.list
  )
  features <- VariableFeatures(object = combined.assay)
  VariableFeatures(object = combined.assay) <- features
  return(combined.assay)
}

#' Subset an AnchorSet object
#'
#' @inheritParams base::subset
#' @param score.threshold Only anchor pairs with scores greater than this value
#' are retained.
#' @param disallowed.dataset.pairs Remove any anchors formed between the
#' provided pairs. E.g. \code{list(c(1, 5), c(1, 2))} filters out any anchors between
#' datasets 1 and 5 and datasets 1 and 2.
#' @param dataset.matrix Provide a binary matrix specifying whether a dataset
#' pair is allowable (1) or not (0). Should be a dataset x dataset matrix.
#' @param group.by Grouping variable to determine allowable ident pairs
#' @param disallowed.ident.pairs Remove any anchors formed between provided
#' ident pairs. E.g. \code{list(c("CD4", "CD8"), c("B-cell", "T-cell"))}
#' @param ident.matrix Provide a binary matrix specifying whether an ident pair
#' is allowable (1) or not (0). Should be an ident x ident symmetric matrix
#'
#' @return Returns an \code{\link{AnchorSet}} object with specified anchors
#' filtered out
#'
#' @export
#' @method subset AnchorSet
#' @concept objects
#'
subset.AnchorSet <- function(
  x,
  score.threshold = NULL,
  disallowed.dataset.pairs = NULL,
  dataset.matrix = NULL,
  group.by = NULL,
  disallowed.ident.pairs = NULL,
  ident.matrix = NULL,
  ...
) {
  if (!is.null(x = disallowed.dataset.pairs) && !is.null(x = dataset.matrix)) {
    stop("Please use either disallowed.dataset.pairs OR dataset.matrix, not both.")
  }
  # Filter based on scores
  if (!is.null(x = score.threshold)) {
    if (score.threshold > 1 | score.threshold < 0) {
      stop(
        "Anchors are scored on a scale between 0 and 1. Please provide a value",
        " in that range to score.threshold."
      )
    }
    anchors <- slot(object = x, name = "anchors")
    anchors <- anchors[anchors[, 'score'] > score.threshold, , drop = FALSE]
    slot(object = x, name = "anchors") <- anchors
  }
  object.names <- names(x = slot(object = x, name = "object.list"))
  num.obs <- length(x = object.names)
  # Filter based on dataset pairings
  if (!is.null(x = disallowed.dataset.pairs)) {
    dataset.matrix <- matrix(data = 1, nrow = num.obs, ncol = num.obs)
    for(i in 1:length(x = disallowed.dataset.pairs)) {
      pair <- disallowed.dataset.pairs[[i]]
      if (length(x = pair) != 2) {
        stop("Please ensure all list items in disallowed.dataset.pairs are of length 2.")
      }
      if (any(pair %in% object.names)) {
        pair[which(pair %in% object.names)] <- sapply(
          X = pair[which(pair %in% object.names)],
          FUN = function(x) {
            which(object.names == x)
          })
      }
      pair <- as.numeric(x = pair)
      dataset.matrix[pair[1], pair[2]] <- 0
    }
  }
  if (!is.null(x = dataset.matrix)) {
    if (any(dim(x = dataset.matrix) != c(num.obs, num.obs))){
      stop("Please provide a dataset.matrix that is ", num.obs, " x ", num.obs, ".")
    }
    anchors <- slot(object = x, name = "anchors")
    pairs <- which(dataset.matrix == 0, arr.ind = TRUE)
    for (i in 1:nrow(x = pairs)) {
      anchors <- anchors[-which(x = anchors$dataset1 == pairs[i, 1] & anchors$dataset2 == pairs[i, 2]), ]
      anchors <- anchors[-which(x = anchors$dataset1 == pairs[i, 2] & anchors$dataset2 == pairs[i, 1]), ]
    }
    slot(object = x, name = "anchors") <- anchors
  }
  # Filter based on ident pairings
  if (!is.null(x = group.by)) {
    anchors <- AnnotateAnchors(anchors = x, vars = group.by)
    if (!is.null(x = disallowed.ident.pairs) && !is.null(x = ident.matrix)) {
      stop("Please use either disallowed.ident.pairs OR ident.matrix, not both.")
    }
    unique.ids <- unique(x = c(
      as.character(x = anchors[, paste0("cell1.", group.by)]),
      as.character(x = anchors[, paste0("cell2.", group.by)]))
    )
    unique.ids <- unique.ids[!is.na(x = unique.ids)]
    num.ids <- length(x = unique.ids)
    if (!is.null(x = disallowed.ident.pairs)) {
      ident.matrix <- matrix(data = 1, nrow = num.ids, ncol = num.ids)
      rownames(x = ident.matrix) <- unique.ids
      colnames(x = ident.matrix) <- unique.ids
      for(i in 1:length(x = disallowed.ident.pairs)) {
        pair <- disallowed.ident.pairs[[i]]
        if (length(x = pair) != 2) {
          stop("Please ensure all list items in disallowed.dataset.pairs are of length 2.")
        }
        ident.matrix[pair[1], pair[2]] <- 0
      }
    }
    if (!is.null(x = ident.matrix)) {
      if (any(dim(x = ident.matrix) != c(num.ids, num.ids))){
        stop("Please provide a dataset.matrix that is ", num.ids, " x ", num.ids, ".")
      }
      to.remove <- c()
      pairs <- which(ident.matrix == 0, arr.ind = TRUE)
      for (i in 1:nrow(x = pairs)) {
        id1 <- rownames(x = ident.matrix)[pairs[i, 1]]
        id2 <- colnames(x = ident.matrix)[pairs[i, 2]]
        to.remove <- c(to.remove, which(x = anchors[, paste0("cell1.", group.by)] == id1 & anchors[, paste0("cell2.", group.by)] == id2))
        to.remove <- c(to.remove, which(x = anchors[, paste0("cell1.", group.by)] == id2 & anchors[, paste0("cell2.", group.by)] == id1))
      }
      anchors <- slot(object = x, name = "anchors")
      anchors <- anchors[-to.remove, ]
      slot(object = x, name = "anchors") <- anchors
    }
  }
  return(x)
}

#' @export
#' @method subset SCTAssay
#' @concept objects
#'
subset.SCTAssay <- function(x, cells = NULL, features = NULL, ...) {
  x <- NextMethod()
  models <- levels(x = x)
  for (m in models) {
    attr <- SCTResults(object = x, slot = "cell.attributes", model = m)
    attr <- attr[intersect(x = rownames(x = attr), y = Cells(x = x)), , drop = FALSE]
    SCTResults(object = x, slot = "cell.attributes", model = m) <- attr
   if (nrow(x = attr) == 0) {
     slot(object = x, name = 'SCTModel.list')[[m]] <- NULL
   }
    }
  return(x)
}

#' @method subset SlideSeq
#' @concept objects
#' @export
#'
subset.SlideSeq <- function(x, cells, ...) {
  x <- subset.VisiumV1(x = x, cells = cells, ...)
  return(x)
}

#' @method subset STARmap
#' @concept objects
#' @export
#'
subset.STARmap <- function(x, cells, ...) {
  x <- subset.VisiumV1(x = x, cells = cells, ...)
  qhulls <- GetTissueCoordinates(object = x, qhulls = TRUE)
  qhulls <- qhulls[qhulls$cell %in% cells, ]
  slot(object = x, name = 'qhulls') <- qhulls
  return(x)
}

#' @method subset VisiumV1
#' @concept objects
#' @export
#'
subset.VisiumV1 <- function(x, cells, ...) {
  coordinates <- GetTissueCoordinates(object = x, scale = NULL, cols = NULL)
  cells <- cells[cells %in% rownames(x = coordinates)]
  coordinates <- coordinates[cells, ]
  slot(object = x, name = 'coordinates') <- coordinates
  return(x)
}

#' Update pre-V4 Assays generated with SCTransform in the Seurat to the new
#' SCTAssay class
#
#' @param object A Seurat object
#' @export
#' @concept objects
#' @return A Seurat object with updated SCTAssays
#'
UpdateSCTAssays <- function(object) {
  assays <- Assays(object = object)
  for (assay in assays) {
    if (IsSCT(assay = object[[assay]]) && !inherits(x = object[[assay]], what = "SCTAssay")) {
      object[[assay]] <- as(object =  object[[assay]], Class = "SCTAssay")
    }
  }
  return(object)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# S4 methods
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @rdname SCTAssay-class
#' @name SCTAssay-class
#'
#' @section Creating an \code{SCTAssay} from an \code{Assay}:
#' Conversion from an \code{Assay} object to an \code{SCTAssay} object by
#' is done by adding the additional slots to the object. If \code{from} has
#' results generated by \code{\link{SCTransform}} from Seurat v3.0.0 to v3.1.1,
#' the conversion will automagically fill the new slots with the data
#'
setAs(
  from = 'Assay',
  to = 'SCTAssay',
  def = function(from) {
    object.list <- sapply(
      X = slotNames(x = from),
      FUN = slot,
      object = from,
      simplify = FALSE,
      USE.NAMES = TRUE
    )
    object.list <- c(
      list('Class' = 'SCTAssay'),
      object.list
    )
    if (IsSCT(assay = from)) {
      vst.slots <- c('vst.set', 'vst.out')
      vst.use <- vst.slots[vst.slots %in% names(x = Misc(object = from))][1]
      vst.res <- Misc(object = from, slot = vst.use)
      umi.assay <- Misc(object = from, slot = "umi.assay")
      if (vst.use == 'vst.out') {
        vst.res <- list(vst.res)
        umi.assay <- list(umi.assay)
      }
      if (length(x = vst.res) == 0) {
        vst.res <- list()
      } else if (length(x = vst.res) > 0) {
        vst.res <- lapply(
          X = 1:length(x = vst.res),
          FUN = function(i) {
            vst.res[[i]]$umi.assay <- umi.assay[[i]]
            return(PrepVSTResults(
              vst.res = vst.res[[i]],
              cell.names = colnames(x = from)
            ))
          }
        )
        names(x = vst.res) <- paste0("model", 1:length(x = vst.res))
      }
      object.list$misc[[vst.use]] <- NULL
      object.list$SCTModel.list <- vst.res
    }
    return(do.call(what = 'new', args = object.list))
  }
)

setMethod(
  f = 'show',
  signature = 'TransferAnchorSet',
  definition = function(object) {
    cat('An AnchorSet object containing', nrow(x = slot(object = object, name = "anchors")),
        "anchors between the reference and query Seurat objects. \n",
        "This can be used as input to TransferData.\n")
  }
)

setMethod(
  f = 'show',
  signature = 'IntegrationAnchorSet',
  definition = function(object) {
    cat('An AnchorSet object containing', nrow(x = slot(object = object, name = "anchors")),
        "anchors between", length(x = slot(object = object, name = "object.list")), "Seurat objects \n",
        "This can be used as input to IntegrateData.\n")
  }
)

setMethod(
  f = 'show',
  signature = 'ModalityWeights',
  definition = function(object) {
    cat(
      'A ModalityWeights object containing modality weights between',
      paste(slot(object = object, name = "modality.assay"), collapse = " and "),
      "assays \n", "This can be used as input to FindMultiModelNeighbors.\n")
  }
)

setMethod(
  f = 'show',
  signature = 'BridgeReferenceSet',
  definition = function(object) {
    cat(
      'A BridgeReferenceSet object has a bridge object with ',
      ncol(slot(object = object, name = 'bridge')),
      'cells and a reference object with ',
      ncol(slot(object = object, name = 'reference')),
      'cells. \n','The bridge query reduction is ',
      slot(object = object, name = 'params')$bridge.query.reduction %||%
        slot(object = object, name = 'params')$supervised.reduction,
   "\n This can be used as input to FindBridgeTransferAnchors and FindBridgeIntegrationAnchors")
  }
)

setMethod(
  f = 'show',
  signature = 'SCTModel',
  definition = function(object) {
    cat(
      "An sctransform model.\n",
      " Model formula: ", slot(object = object, name = "model"),
      "\n  Parameters stored for", nrow(x = SCTResults(object = object, slot = "feature.attributes")), "features,",
      nrow(x = SCTResults(object = object, slot = "cell.attributes")), "cells.\n")
  }
)

#' @importFrom utils head
#
setMethod(
  f = 'show',
  signature = 'SCTAssay',
  definition = function(object) {
    cat('SCTAssay data with', nrow(x = object), 'features for', ncol(x = object),
        'cells, and', length(x = levels(x = object)) , 'SCTModel(s) \n')

    if (length(x = VariableFeatures(object = object)) > 0) {
      top.ten <- head(x = VariableFeatures(object = object), n = 10L)
      top <- 'Top'
      variable <- 'variable'
    } else {
      top.ten <- head(x = rownames(x = object), n = 10L)
      top <- 'First'
      variable <- ''
    }
    features <- paste0(
      variable,
      ' feature',
      if (length(x = top.ten) != 1) {'s'}, ":\n"
    )
    features <- gsub(pattern = '^\\s+', replacement = '', x = features)
    cat(
      top,
      length(x = top.ten),
      features,
      paste(strwrap(x = paste(top.ten, collapse = ', ')), collapse = '\n'),
      '\n'
    )
  }
)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Internal
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Internal AddMetaData defintion
#
# @param object An object
# @param metadata A vector, list, or data.frame with metadata to add
# @param col.name A name for meta data if not a named list or data.frame
#
# @return object with metadata added
#
.AddMetaData <- function(object, metadata, col.name = NULL) {
  object <- UpdateSlots(object = object)
  if (is.null(x = col.name) && is.atomic(x = metadata)) {
    stop("'col.name' must be provided for atomic metadata types (eg. vectors)")
  }
  if (inherits(x = metadata, what = c('matrix', 'Matrix'))) {
    metadata <- as.data.frame(x = metadata)
  }
  col.name <- col.name %||% names(x = metadata) %||% colnames(x = metadata)
  if (is.null(x = col.name)) {
    stop("No metadata name provided and could not infer it from metadata object")
  }
  object[[col.name]] <- metadata
  # if (class(x = metadata) == "data.frame") {
  #   for (ii in 1:ncol(x = metadata)) {
  #     object[[colnames(x = metadata)[ii]]] <- metadata[, ii, drop = FALSE]
  #   }
  # } else {
  #   object[[col.name]] <- metadata
  # }
  return(object)
}

# Find the names of collections in an object
#
# @return A vector with the names of slots that are a list
#
Collections <- function(object) {
  collections <- vapply(
    X = slotNames(x = object),
    FUN = function(x) {
      return(any(grepl(pattern = 'list', x = class(x = slot(object = object, name = x)))))
    },
    FUN.VALUE = logical(length = 1L)
  )
  collections <- Filter(f = isTRUE, x = collections)
  return(names(x = collections))
}

# Get the default image of an object
#
# Attempts to find all images associated with the default assay of the object.
# If none present, finds all images present in the object. Returns the name of
# the first image
#
# @param object A Seurat object
#
# @return The name of the default image
#
DefaultImage <- function(object) {
  object <- UpdateSlots(object = object)
  images <- Images(object = object, assay = DefaultAssay(object = object))
  if (length(x = images) < 1) {
    images <- Images(object = object)
  }
  return(images[[1]])
}

# Get the names of objects within a Seurat object that are of a certain class
#
# @param object A Seurat object
# @param classes.keep A vector of names of classes to get
#
# @return A vector with the names of objects within the Seurat object that are of class \code{classes.keep}
#
#' @importFrom stats na.omit
#
FilterObjects <- function(object, classes.keep = c('Assay', 'DimReduc')) {
  object <- UpdateSlots(object = object)
  slots <- na.omit(object = Filter(
    f = function(x) {
      sobj <- slot(object = object, name = x)
      return(is.list(x = sobj) && !is.data.frame(x = sobj) && !is.package_version(x = sobj))
    },
    x = slotNames(x = object)
  ))
  slots <- grep(pattern = 'tools', x = slots, value = TRUE, invert = TRUE)
  slots <- grep(pattern = 'misc', x = slots, value = TRUE, invert = TRUE)
  slots.objects <- unlist(
    x = lapply(
      X = slots,
      FUN = function(x) {
        return(names(x = slot(object = object, name = x)))
      }
    ),
    use.names = FALSE
  )
  object.classes <- sapply(
    X = slots.objects,
    FUN = function(i) {
      return(inherits(x = object[[i]], what = classes.keep))
    }
  )
  object.classes <- which(x = object.classes, useNames = TRUE)
  return(names(x = object.classes))
}

# Find the collection of an object within a Seurat object
#
# @param object A Seurat object
# @param name Name of object to find
#
# @return The collection (slot) of the object
#
FindObject <- function(object, name) {
  collections <- c(
    'assays',
    'graphs',
    'neighbors',
    'reductions',
    'commands',
    'images'
  )
  object.names <- lapply(
    X = collections,
    FUN = function(x) {
      return(names(x = slot(object = object, name = x)))
    }
  )
  names(x = object.names) <- collections
  object.names <- Filter(f = Negate(f = is.null), x = object.names)
  for (i in names(x = object.names)) {
    if (name %in% names(x = slot(object = object, name = i))) {
      return(i)
    }
  }
  return(NULL)
}

# Prepare VST results for use with SCTAssay objects
#
# @param vst.res Results from sctransform::vst
# @param cell.names Vector of valid cell names still in object
#
# @return An SCTModel object.
#
#
PrepVSTResults <- function(vst.res, cell.names) {
  # Prepare cell attribute information
  cell.attrs <- vst.res$cell_attr
  cell.names <- intersect(x = cell.names, y = rownames(x = cell.attrs))
  cell.cols <- c(
    'umi',
    'gene',
    'log_umi',
    'log_gene',
    'umi_per_gene',
    'log_umi_per_gene'
  )
  cell.cols <- intersect(x = cell.cols, y = colnames(x = cell.attrs))
  cell.attrs <- cell.attrs[cell.names, cell.cols, drop = FALSE]
  colnames(x = cell.attrs) <- gsub(
    pattern = 'gene',
    replacement = 'feature',
    x = colnames(x = cell.attrs)
  )
  if (!is.null(x = vst.res$cells_step1)) {
    cell.attrs[, "cells_step1"] <- FALSE
    cells_step1 <- intersect(x = vst.res$cells_step1,
                             y = rownames(x = cell.attrs))
    cell.attrs[cells_step1, "cells_step1"] <- TRUE
  }
  # Prepare feature attribute information
  feature.attrs <- vst.res$gene_attr
  feature.cols <- c(
    'detection_rate',
    'gmean',
    'variance',
    'residual_mean',
    'residual_variance'
  )
  feature.cols <- intersect(x = feature.cols, y = colnames(x = feature.attrs))
  feature.attrs <- feature.attrs[, feature.cols, drop = FALSE]
  feature.attrs <- cbind(feature.attrs, vst.res$model_pars_fit[rownames(feature.attrs), , drop = FALSE])

  if (!is.null(x = vst.res$genes_log_gmean_step1)) {
    feature.attrs[,"genes_log_gmean_step1"] <- FALSE
    genes_step1 <- intersect(
      x = names(vst.res$genes_log_gmean_step1),
      y = rownames(feature.attrs)
    )
    feature.attrs[genes_step1,"genes_log_gmean_step1"] <- TRUE

    # add parameters from step1
    feature.attrs[, paste0("step1_", colnames(vst.res$model_pars))] <- NA
    feature.attrs[genes_step1, paste0("step1_", colnames(vst.res$model_pars))] <- vst.res$model_pars[genes_step1,]

  }
  # Prepare clipping information
  clips <- list(
    'vst' = vst.res$arguments$res_clip_range,
    'sct' = vst.res$arguments$sct.clip.range
  )
  median_umi <- NA
  # check if a custom scale_factor was provided to vst()
  if ("scale_factor" %in% names(vst.res$arguments)){
    median_umi <- vst.res$arguments$scale_factor
  }
  if (is.na(median_umi)) {
    if ("umi" %in% colnames(x = cell.attrs)) {
      median_umi <- median(cell.attrs$umi)
    } else if ("log_umi" %in% colnames(x = cell.attrs)) {
      median_umi <- median(10 ^ cell.attrs$log_umi)
    }
  }
  vst.res.SCTModel  <- SCTModel(
    feature.attributes = feature.attrs,
    cell.attributes = cell.attrs,
    clips = clips,
    umi.assay = vst.res$umi.assay %||% "RNA",
    model =  vst.res$model_str,
    arguments =  vst.res$arguments,
    median_umi = median_umi
  )
  return(vst.res.SCTModel)
}

# Return a null image
#
# @param mode Image representation to return
# see \code{\link{GetImage}} for more details
#
#' @importFrom grid nullGrob
#' @importFrom grDevices as.raster
#
NullImage <- function(mode) {
  image <- switch(
    EXPR = mode,
    'grob' = nullGrob(),
    'raster' = as.raster(x = new(Class = 'matrix')),
    'plotly' = list('visible' = FALSE),
    'raw' = NULL,
    stop("Unknown image mode: ", mode, call. = FALSE)
  )
  return(image)
}

# Check to see if projected loadings have been set
#
# @param object a DimReduc object
#
# @return TRUE if proejcted loadings have been set, else FALSE
#
Projected <- function(object) {
  projected.dims <- dim(x = slot(object = object, name = 'feature.loadings.projected'))
  if (all(projected.dims == 1)) {
    return(!all(is.na(x = slot(object = object, name = 'feature.loadings.projected'))))
  }
  return(!all(projected.dims == 0))
}

# Subset cells in vst data
# @param sct.info A vst.out list
# @param cells vector of cells to retain
# @param features vector of features to retain
SubsetVST <- function(sct.info, cells, features) {
  cells.keep <- intersect(x = cells, y = rownames(x = sct.info$cell_attr))
  sct.info$cell_attr <- sct.info$cell_attr[cells.keep, ]
  # find which subset of features are in the SCT assay
  feat.keep <- intersect(x = features, y = rownames(x = sct.info$gene_attr))
  sct.info$gene_attr <- sct.info$gene_attr[feat.keep, ]
  return(sct.info)
}

# Get the top
#
# @param data Data to pull the top from
# @param num Pull top \code{num}
# @param balanced Pull even amounts of from positive and negative values
#
# @return The top \code{num}
# @seealso \{code{\link{TopCells}}} \{code{\link{TopFeatures}}}
#
#' @importFrom utils head tail
#
Top <- function(data, num, balanced) {
  nr <- nrow(x = data)
  if (num > nr) {
    warning("Requested number is larger than the number of available items (",
            nr, "). Setting to ", nr , ".", call. = FALSE)
    num <- nr
  }
  if (num == 1) {
    balanced <- FALSE
  }
  top <- if (balanced) {
    num <- round(x = num / 2)
    data <- data[order(data, decreasing = TRUE), , drop = FALSE]
    positive <- head(x = rownames(x = data), n = num)
    negative <- rev(x = tail(x = rownames(x = data), n = num))

    # remove duplicates
    if (positive[num] == negative[num]) {
      negative <- negative[-num]
    }
    list(positive = positive, negative = negative)
  } else {
    data <- data[rev(x = order(abs(x = data))), , drop = FALSE]
    top <- head(x = rownames(x = data), n = num)
    top[order(data[top, ])]
  }
  return(top)
}

# Update Seurat assay
#
# @param old.assay Seurat2 assay
# @param assay Name to store for assay in new object
#
UpdateAssay <- function(old.assay, assay){
  cells <- colnames(x = old.assay@data)
  counts <- old.assay@raw.data
  data <- old.assay@data
  if (!inherits(x = counts, what = 'dgCMatrix')) {
    counts <- as.sparse(x = as.matrix(x = counts))
  }
  if (!inherits(x = data, what = 'dgCMatrix')) {
    data <- as.sparse(x = as.matrix(x = data))
  }
  new.assay <- new(
    Class = 'Assay',
    counts = counts[, cells],
    data = data,
    scale.data = old.assay@scale.data %||% new(Class = 'matrix'),
    meta.features = data.frame(row.names = rownames(x = counts)),
    var.features = old.assay@var.genes,
    key = paste0(assay, "_")
  )
  return(new.assay)
}

# Update dimension reduction
#
# @param old.dr Seurat2 dimension reduction slot
# @param assay.used Name of assay used to compute dimension reduction
#
UpdateDimReduction <- function(old.dr, assay) {
  new.dr <- list()
  for (i in names(x = old.dr)) {
    cell.embeddings <- old.dr[[i]]@cell.embeddings %||% new(Class = 'matrix')
    feature.loadings <- old.dr[[i]]@gene.loadings %||% new(Class = 'matrix')
    stdev <- old.dr[[i]]@sdev %||% numeric()
    misc <- old.dr[[i]]@misc %||% list()
    new.jackstraw <- UpdateJackstraw(old.jackstraw = old.dr[[i]]@jackstraw)
    old.key <- old.dr[[i]]@key
    if (length(x = old.key) == 0) {
      old.key <- gsub(pattern = "(.+?)(([0-9]+).*)", replacement = "\\1",  x = colnames(cell.embeddings)[[1]])
      if (length(x = old.key) == 0) {
        old.key <- i
      }
    }
    new.key <- suppressWarnings(expr = UpdateKey(key = old.key))
    colnames(x = cell.embeddings) <- gsub(
      pattern = old.key,
      replacement = new.key,
      x = colnames(x = cell.embeddings)
    )
    colnames(x = feature.loadings) <- gsub(
      pattern = old.key,
      replacement = new.key,
      x = colnames(x = feature.loadings)
    )
    new.dr[[i]] <- new(
      Class = 'DimReduc',
      cell.embeddings = as(object = cell.embeddings, Class = 'matrix'),
      feature.loadings = as(object = feature.loadings, Class = 'matrix'),
      assay.used = assay,
      stdev = as(object = stdev, Class = 'numeric'),
      key = as(object = new.key, Class = 'character'),
      jackstraw = new.jackstraw,
      misc = as(object = misc, Class = 'list')
    )
  }
  return(new.dr)
}

# Update jackstraw
#
# @param old.jackstraw
#
UpdateJackstraw <- function(old.jackstraw) {
  if (is.null(x = old.jackstraw)) {
    new.jackstraw <- new(
      Class = 'JackStrawData',
      empirical.p.values = new(Class = 'matrix'),
      fake.reduction.scores = new(Class = 'matrix'),
      empirical.p.values.full = new(Class = 'matrix'),
      overall.p.values = new(Class = 'matrix')
    )
  } else {
    if (.hasSlot(object = old.jackstraw, name = 'overall.p.values')) {
      overall.p <- old.jackstraw@overall.p.values %||% new(Class = 'matrix')
    } else {
      overall.p <- new(Class = 'matrix')
    }
    new.jackstraw <- new(
      Class = 'JackStrawData',
      empirical.p.values = old.jackstraw@emperical.p.value %||% new(Class = 'matrix'),
      fake.reduction.scores = old.jackstraw@fake.pc.scores %||% new(Class = 'matrix'),
      empirical.p.values.full = old.jackstraw@emperical.p.value.full %||% new(Class = 'matrix'),
      overall.p.values = overall.p
    )
  }
  return(new.jackstraw)
}

# Update a Key
#
# @param key A character to become a Seurat Key
#
# @return An updated Key that's valid for Seurat
#
UpdateKey <- function(key) {
  if (grepl(pattern = '^[[:alnum:]]+_$', x = key)) {
    return(key)
  } else {
    new.key <- regmatches(
      x = key,
      m = gregexpr(pattern = '[[:alnum:]]+', text = key)
    )
    new.key <- paste0(paste(unlist(x = new.key), collapse = ''), '_')
    if (new.key == '_') {
      new.key <- paste0(RandomName(length = 3), '_')
    }
    warning(
      "Keys should be one or more alphanumeric characters followed by an underscore, setting key from ",
      key,
      " to ",
      new.key,
      call. = FALSE,
      immediate. = TRUE
    )
    return(new.key)
  }
}

# Update slots in an object
#
# @param object An object to update
#
# @return \code{object} with the latest slot definitions
#
UpdateSlots <- function(object) {
  object.list <- sapply(
    X = slotNames(x = object),
    FUN = function(x) {
      return(tryCatch(
        expr = slot(object = object, name = x),
        error = function(...) {
          return(NULL)
        }
      ))
    },
    simplify = FALSE,
    USE.NAMES = TRUE
  )
  object.list <- Filter(f = Negate(f = is.null), x = object.list)
  object.list <- c('Class' = class(x = object)[1], object.list)
  object <- rlang::invoke(
     .fn = new,
     .args  = object.list
   )
  for (x in setdiff(x = slotNames(x = object), y = names(x = object.list))) {
    xobj <- slot(object = object, name = x)
    if (is.vector(x = xobj) && !is.list(x = xobj) && length(x = xobj) == 0) {
      slot(object = object, name = x) <- vector(mode = class(x = xobj), length = 1L)
    }
  }
  return(object)
}

# Pulls the proper data matrix for merging assay data. If the slot is empty, will return an empty
# matrix with the proper dimensions from one of the remaining data slots.
#
# @param assay Assay to pull data from
# @param slot Slot to pull from
#
# @return Returns the data matrix if present (i.e.) not 0x0. Otherwise, returns an
# appropriately sized empty sparse matrix
#
#' @importFrom Matrix Matrix
#
ValidateDataForMerge <- function(assay, slot) {
  mat <- GetAssayData(object = assay, slot = slot)
  if (any(dim(x = mat) == c(0, 0))) {
    slots.to.check <- setdiff(x = c("counts", "data", "scale.data"), y = slot)
    for (ss in slots.to.check) {
      data.dims <- dim(x = GetAssayData(object = assay, slot = ss))
      data.slot <- ss
      if (!any(data.dims == c(0, 0))) {
        break
      }
    }
    if (any(data.dims == c(0, 0))) {
      stop("The counts, data, and scale.data slots are all empty for the provided assay.")
    }
    mat <- Matrix(
      data = 0,
      nrow = data.dims[1],
      ncol = data.dims[2],
      dimnames = dimnames(x = GetAssayData(object = assay, slot = data.slot))
    )
    mat <- as.sparse(x = mat)
  }
  return(mat)
}
