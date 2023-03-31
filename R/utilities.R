#' @include generics.R
#' @importFrom SeuratObject PackageCheck
#'
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Add Azimuth Results
#'
#' Add mapping and prediction scores, UMAP embeddings, and imputed assay (if
#' available)
#' from Azimuth to an existing or new \code{\link[SeuratObject]{Seurat}} object
#'
#' @param object A \code{\link[SeuratObject]{Seurat}} object
#' @param filename Path to Azimuth mapping scores file
#'
#' @return \code{object} with Azimuth results added
#'
#' @examples
#' \dontrun{
#' object <- AddAzimuthResults(object, filename = "azimuth_results.Rds")
#' }
#'
#' @export
AddAzimuthResults <- function(object = NULL, filename) {
  if (is.null(x = filename)) {
    stop("No Azimuth results provided.")
  }
  azimuth_results <- readRDS(file = filename)
  if (!is.list(x = azimuth_results) || any(!(c('umap', 'pred.df') %in% names(x = azimuth_results)))) {
    stop("Expected following format for azimuth_results:
           `list(umap = <DimReduc>, pred.df = <data.frame>[, impADT = <Assay>])`")
  }

  if (is.null(x = object)) {
    message("No existing Seurat object provided. Creating new one.")
    object <- CreateSeuratObject(
      counts = matrix(
        nrow = 1,
        ncol = nrow(x = azimuth_results$umap),
        dimnames = list(
          row.names = 'Dummy.feature',
          col.names = rownames(x = azimuth_results$umap))
      ),
      assay = 'Dummy'
    )
  } else {
    overlap.cells <- intersect(
      x = Cells(x = object),
      y = rownames(x = azimuth_results$umap)
    )
    if (!(all(overlap.cells %in% Cells(x = object)))) {
      stop("Cells in object do not match cells in download")
    } else if (length(x = overlap.cells) < length(x = Cells(x = object))) {
      warning(paste0("Subsetting out ", length(x = Cells(x = object)) - length(x = overlap.cells),
                     " cells that are absent in downloaded results (perhaps filtered by Azimuth)"))
      object <- subset(x = object, cells = overlap.cells)
    }
  }

  azimuth_results$pred.df$cell <- NULL
  object <- AddMetaData(object = object, metadata = azimuth_results$pred.df)
  object[['umap.proj']] <- azimuth_results$umap
  if ('impADT' %in% names(x = azimuth_results)) {
    object[['impADT']] <- azimuth_results$impADT
    if ('Dummy' %in% Assays(object = object)) {
      DefaultAssay(object = object) <- 'impADT'
      object[['Dummy']] <- NULL
    }
  }
  return(object)
}

#' Add Azimuth Scores
#'
#' Add mapping and prediction scores from Azimuth to a
#' \code{\link[SeuratObject]{Seurat}} object
#'
#' @param object A \code{\link[SeuratObject]{Seurat}} object
#' @param filename Path to Azimuth mapping scores file
#'
#' @return \code{object} with the mapping scores added
#'
#' @examples
#' \dontrun{
#' object <- AddAzimuthScores(object, filename = "azimuth_pred.tsv")
#' }
#'
AddAzimuthScores <- function(object, filename) {
  if (!file.exists(filename)) {
    stop("Cannot find Azimuth scores file ", filename, call. = FALSE)
  }
  object <- AddMetaData(
    object = object,
    metadata = read.delim(file = filename, row.names = 1)
  )
  return(object)
}

#' Calculate module scores for feature expression programs in single cells
#'
#' Calculate the average expression levels of each program (cluster) on single
#' cell level, subtracted by the aggregated expression of control feature sets.
#' All analyzed features are binned based on averaged expression, and the
#' control features are randomly selected from each bin.
#'
#' @param object Seurat object
#' @param features A list of vectors of features for expression programs; each
#' entry should be a vector of feature names
#' @param pool List of features to check expression levels against, defaults to
#' \code{rownames(x = object)}
#' @param nbin Number of bins of aggregate expression levels for all
#' analyzed features
#' @param ctrl Number of control features selected from the same bin per
#' analyzed feature
#' @param k Use feature clusters returned from DoKMeans
#' @param assay Name of assay to use
#' @param name Name for the expression programs; will append a number to the
#' end for each entry in \code{features} (eg. if \code{features} has three
#' programs, the results will be stored as \code{name1}, \code{name2},
#' \code{name3}, respectively)
#' @param seed Set a random seed. If NULL, seed is not set.
#' @param search Search for symbol synonyms for features in \code{features} that
#' don't match features in \code{object}? Searches the HGNC's gene names
#' database; see \code{\link{UpdateSymbolList}} for more details
#' @param ... Extra parameters passed to \code{\link{UpdateSymbolList}}
#'
#' @return Returns a Seurat object with module scores added to object meta data;
#' each module is stored as \code{name#} for each module program present in
#' \code{features}
#'
#' @importFrom ggplot2 cut_number
#' @importFrom Matrix rowMeans colMeans
#'
#' @references Tirosh et al, Science (2016)
#'
#' @export
#' @concept utilities
#'
#' @examples
#' \dontrun{
#' data("pbmc_small")
#' cd_features <- list(c(
#'   'CD79B',
#'   'CD79A',
#'   'CD19',
#'   'CD180',
#'   'CD200',
#'   'CD3D',
#'   'CD2',
#'   'CD3E',
#'   'CD7',
#'   'CD8A',
#'   'CD14',
#'   'CD1C',
#'   'CD68',
#'   'CD9',
#'   'CD247'
#' ))
#' pbmc_small <- AddModuleScore(
#'   object = pbmc_small,
#'   features = cd_features,
#'   ctrl = 5,
#'   name = 'CD_Features'
#' )
#' head(x = pbmc_small[])
#' }
#'
AddModuleScore <- function(
  object,
  features,
  pool = NULL,
  nbin = 24,
  ctrl = 100,
  k = FALSE,
  assay = NULL,
  name = 'Cluster',
  seed = 1,
  search = FALSE,
  ...
) {
  if (!is.null(x = seed)) {
    set.seed(seed = seed)
  }
  assay.old <- DefaultAssay(object = object)
  assay <- assay %||% assay.old
  DefaultAssay(object = object) <- assay
  assay.data <- GetAssayData(object = object)
  features.old <- features
  if (k) {
    .NotYetUsed(arg = 'k')
    features <- list()
    for (i in as.numeric(x = names(x = table(object@kmeans.obj[[1]]$cluster)))) {
      features[[i]] <- names(x = which(x = object@kmeans.obj[[1]]$cluster == i))
    }
    cluster.length <- length(x = features)
  } else {
    if (is.null(x = features)) {
      stop("Missing input feature list")
    }
    features <- lapply(
      X = features,
      FUN = function(x) {
        missing.features <- setdiff(x = x, y = rownames(x = object))
        if (length(x = missing.features) > 0) {
          warning(
            "The following features are not present in the object: ",
            paste(missing.features, collapse = ", "),
            ifelse(
              test = search,
              yes = ", attempting to find updated synonyms",
              no = ", not searching for symbol synonyms"
            ),
            call. = FALSE,
            immediate. = TRUE
          )
          if (search) {
            tryCatch(
              expr = {
                updated.features <- UpdateSymbolList(symbols = missing.features, ...)
                names(x = updated.features) <- missing.features
                for (miss in names(x = updated.features)) {
                  index <- which(x == miss)
                  x[index] <- updated.features[miss]
                }
              },
              error = function(...) {
                warning(
                  "Could not reach HGNC's gene names database",
                  call. = FALSE,
                  immediate. = TRUE
                )
              }
            )
            missing.features <- setdiff(x = x, y = rownames(x = object))
            if (length(x = missing.features) > 0) {
              warning(
                "The following features are still not present in the object: ",
                paste(missing.features, collapse = ", "),
                call. = FALSE,
                immediate. = TRUE
              )
            }
          }
        }
        return(intersect(x = x, y = rownames(x = object)))
      }
    )
    cluster.length <- length(x = features)
  }
  if (!all(LengthCheck(values = features))) {
    warning(paste(
      'Could not find enough features in the object from the following feature lists:',
      paste(names(x = which(x = !LengthCheck(values = features)))),
      'Attempting to match case...'
    ))
    features <- lapply(
      X = features.old,
      FUN = CaseMatch,
      match = rownames(x = object)
    )
  }
  if (!all(LengthCheck(values = features))) {
    stop(paste(
      'The following feature lists do not have enough features present in the object:',
      paste(names(x = which(x = !LengthCheck(values = features)))),
      'exiting...'
    ))
  }
  pool <- pool %||% rownames(x = object)
  data.avg <- Matrix::rowMeans(x = assay.data[pool, ])
  data.avg <- data.avg[order(data.avg)]
  data.cut <- cut_number(x = data.avg + rnorm(n = length(data.avg))/1e30, n = nbin, labels = FALSE, right = FALSE)
  #data.cut <- as.numeric(x = Hmisc::cut2(x = data.avg, m = round(x = length(x = data.avg) / (nbin + 1))))
  names(x = data.cut) <- names(x = data.avg)
  ctrl.use <- vector(mode = "list", length = cluster.length)
  for (i in 1:cluster.length) {
    features.use <- features[[i]]
    for (j in 1:length(x = features.use)) {
      ctrl.use[[i]] <- c(
        ctrl.use[[i]],
        names(x = sample(
          x = data.cut[which(x = data.cut == data.cut[features.use[j]])],
          size = ctrl,
          replace = FALSE
        ))
      )
    }
  }
  ctrl.use <- lapply(X = ctrl.use, FUN = unique)
  ctrl.scores <- matrix(
    data = numeric(length = 1L),
    nrow = length(x = ctrl.use),
    ncol = ncol(x = object)
  )
  for (i in 1:length(ctrl.use)) {
    features.use <- ctrl.use[[i]]
    ctrl.scores[i, ] <- Matrix::colMeans(x = assay.data[features.use, ])
  }
  features.scores <- matrix(
    data = numeric(length = 1L),
    nrow = cluster.length,
    ncol = ncol(x = object)
  )
  for (i in 1:cluster.length) {
    features.use <- features[[i]]
    data.use <- assay.data[features.use, , drop = FALSE]
    features.scores[i, ] <- Matrix::colMeans(x = data.use)
  }
  features.scores.use <- features.scores - ctrl.scores
  rownames(x = features.scores.use) <- paste0(name, 1:cluster.length)
  features.scores.use <- as.data.frame(x = t(x = features.scores.use))
  rownames(x = features.scores.use) <- colnames(x = object)
  object[[colnames(x = features.scores.use)]] <- features.scores.use
  CheckGC()
  DefaultAssay(object = object) <- assay.old
  return(object)
}

#' Aggregated feature expression by identity class
#'
#' Returns aggregated (summed) expression values for each identity class
#'
#' If slot is set to 'data', this function assumes that the data has been log
#' normalized and therefore feature values are exponentiated prior to aggregating
#' so that sum is done in non-log space. Otherwise, if slot is set to
#' either 'counts' or 'scale.data', no exponentiation is performed prior to
#' aggregating
#' If \code{return.seurat = TRUE} and slot is not 'scale.data', aggregated values
#' are placed in the 'counts' slot of the returned object and the log of aggregated values
#' are placed in the 'data' slot. For the \code{\link{ScaleData}} is then run on the default assay
#' before returning the object.
#' If \code{return.seurat = TRUE} and slot is 'scale.data', the 'counts' slot is left empty,
#' the 'data' slot is filled with NA, and 'scale.data' is set to the aggregated values.
#'
#' @param object Seurat object
#' @param assays Which assays to use. Default is all assays
#' @param features Features to analyze. Default is all features in the assay
#' @param return.seurat Whether to return the data as a Seurat object. Default is FALSE
#' @param group.by Categories for grouping (e.g, ident, replicate, celltype); 'ident' by default
#' @param add.ident (Deprecated) Place an additional label on each cell prior to pseudobulking
#' (very useful if you want to observe cluster pseudobulk values, separated by replicate, for example)
#' @param slot Slot(s) to use; if multiple slots are given, assumed to follow
#' the order of 'assays' (if specified) or object's assays
#' @param verbose Print messages and show progress bar
#' @param ... Arguments to be passed to methods such as \code{\link{CreateSeuratObject}}#'
#' @return Returns a matrix with genes as rows, identity classes as columns.
#' If return.seurat is TRUE, returns an object of class \code{\link{Seurat}}.
#' @export
#' @concept utilities
#'
#' @examples
#' data("pbmc_small")
#' head(AggregateExpression(object = pbmc_small))
#'
AggregateExpression <- function(
  object,
  assays = NULL,
  features = NULL,
  return.seurat = FALSE,
  group.by = 'ident',
  add.ident = NULL,
  slot = 'data',
  verbose = TRUE,
  ...
) {
  return(
    AverageExpression(
      object = object,
      method = 'aggregate',
      assays = assays,
      features = features,
      return.seurat = return.seurat,
      group.by = group.by,
      add.ident = add.ident,
      slot = slot,
      verbose = verbose,
      ...
    )
  )
}

#' Averaged feature expression by identity class
#'
#' Returns averaged expression values for each identity class
#'
#' If slot is set to 'data', this function assumes that the data has been log
#' normalized and therefore feature values are exponentiated prior to averaging
#' so that averaging is done in non-log space. Otherwise, if slot is set to
#' either 'counts' or 'scale.data', no exponentiation is performed prior to
#' averaging
#' If \code{return.seurat = TRUE} and slot is not 'scale.data', averaged values
#' are placed in the 'counts' slot of the returned object and the log of averaged values
#' are placed in the 'data' slot. \code{\link{ScaleData}} is then run on the default assay
#' before returning the object.
#' If \code{return.seurat = TRUE} and slot is 'scale.data', the 'counts' slot is left empty,
#' the 'data' slot is filled with NA, and 'scale.data' is set to the aggregated values.
#'
#' @param object Seurat object
#' @param assays Which assays to use. Default is all assays
#' @param features Features to analyze. Default is all features in the assay
#' @param return.seurat Whether to return the data as a Seurat object. Default is FALSE
#' @param group.by Categories for grouping (e.g, ident, replicate, celltype); 'ident' by default
#' @param add.ident (Deprecated) Place an additional label on each cell prior to pseudobulking
#' (very useful if you want to observe cluster pseudobulk values, separated by replicate, for example)
#' @param slot Slot(s) to use; if multiple slots are given, assumed to follow
#' the order of 'assays' (if specified) or object's assays
#' @param verbose Print messages and show progress bar
#' @param ... Arguments to be passed to methods such as \code{\link{CreateSeuratObject}}
#'
#' @return Returns a matrix with genes as rows, identity classes as columns.
#' If return.seurat is TRUE, returns an object of class \code{\link{Seurat}}.
#' @export
#' @concept utilities
#'
#' @examples
#' data("pbmc_small")
#' head(AverageExpression(object = pbmc_small))
#'
AverageExpression <- function(
  object,
  assays = NULL,
  features = NULL,
  return.seurat = FALSE,
  group.by = 'ident',
  add.ident = NULL,
  slot = 'counts',
  method = 'average',
  verbose = TRUE,
  ...
) {
  CheckDots(..., fxns = 'CreateSeuratObject')
  if (!is.null(x = add.ident)) {
    .Deprecated(msg = "'add.ident' is a deprecated argument, please use the 'group.by' argument instead")
    group.by <- c('ident', add.ident)
  }
  if (!(method %in% c('average', 'aggregate'))) {
    stop("'method' must be either 'average' or 'aggregate'")
  }
  object.assays <- FilterObjects(object = object, classes.keep = c('Assay', 'Assay5'))
  assays <- assays %||% object.assays
  if (!all(assays %in% object.assays)) {
    assays <- assays[assays %in% object.assays]
    if (length(x = assays) == 0) {
      stop("None of the requested assays are present in the object")
    } else {
      warning("Requested assays that do not exist in object. Proceeding with existing assays only.")
    }
  }
  if (length(x = slot) == 1) {
    slot <- rep_len(x = slot, length.out = length(x = assays))
  } else if (length(x = slot) != length(x = assays)) {
    stop("Number of slots provided does not match number of assays")
  }
  data <- FetchData(object = object, vars = rev(x = group.by))
  data <- data[which(rowSums(x = is.na(x = data)) == 0), , drop = F]
  if (nrow(x = data) < ncol(x = object)) {
    message("Removing cells with NA for 1 or more grouping variables")
    object <- subset(x = object, cells = rownames(x = data))
  }
  for (i in 1:ncol(x = data)) {
    data[, i] <- as.factor(x = data[, i])
  }
  num.levels <- sapply(
    X = 1:ncol(x = data),
    FUN = function(i) {
      length(x = levels(x = data[, i]))
    }
  )
  if (any(num.levels == 1)) {
    message(paste0("The following grouping variables have 1 value and will be ignored: ",
                   paste0(colnames(x = data)[which(num.levels <= 1)], collapse = ", ")))
    group.by <- colnames(x = data)[which(num.levels > 1)]
    data <- data[, which(num.levels > 1), drop = F]
  }
  category.matrix <- CreateCategoryMatrix(labels = data, method = method)
  data.return <- list()
  for (i in 1:length(x = assays)) {
    if (inherits(x = features, what = "list")) {
      features.i <- features[[i]]
    } else {
      features.i <- features
    }
    data.return[[assays[i]]] <- PseudobulkExpression(
      object = object[[assays[i]]],
      assay = assays[i],
      category.matrix = category.matrix,
      features = features.i,
      slot = slot[i],
      verbose = verbose,
      ...
    )
  }
  if (return.seurat) {
    op <- options(Seurat.object.assay.version = "v3", Seurat.object.assay.calcn = FALSE)
    on.exit(expr = options(op), add = TRUE)
    if (slot[1] == 'scale.data') {
      na.matrix <- as.matrix(x = data.return[[1]])
      na.matrix[1:length(x = na.matrix)] <- NA
      toRet <- CreateSeuratObject(
        counts = na.matrix,
        project = if (method == "average") "Average" else "Aggregate",
        assay = names(x = data.return)[1],
        check.matrix = FALSE,
        ...
      )
      toRet <- SetAssayData(
        object = toRet,
        assay = names(x = data.return)[1],
        slot = "counts",
        new.data = matrix()
      )
      toRet <- SetAssayData(
        object = toRet,
        assay = names(x = data.return)[1],
        slot = "data",
        new.data = na.matrix
      )
      toRet <- SetAssayData(
        object = toRet,
        assay = names(x = data.return)[1],
        slot = "scale.data",
        new.data = data.return[[1]]
      )
    } else {
      toRet <- CreateSeuratObject(
        counts = data.return[[1]],
        project = if (method == "average") "Average" else "Aggregate",
        assay = names(x = data.return)[1],
        check.matrix = FALSE,
        ...
      )
      toRet <- SetAssayData(
        object = toRet,
        assay = names(x = data.return)[1],
        slot = "data",
        new.data = log1p(x = as.matrix(x = data.return[[1]]))
      )
    }
    #for multimodal data
    if (length(x = data.return) > 1) {
      for (i in 2:length(x = data.return)) {
        if (slot[i] == 'scale.data') {
          na.matrix <- as.matrix(x = data.return[[i]])
          na.matrix[1:length(x = na.matrix)] <- NA
          toRet[[names(x = data.return)[i]]] <- CreateAssayObject(counts = na.matrix, check.matrix = FALSE)
          toRet <- SetAssayData(
            object = toRet,
            assay = names(x = data.return)[i],
            slot = "counts",
            new.data = matrix()
          )
          toRet <- SetAssayData(
            object = toRet,
            assay = names(x = data.return)[i],
            slot = "data",
            new.data = na.matrix
          )
          toRet <- SetAssayData(
            object = toRet,
            assay = names(x = data.return)[i],
            slot = "scale.data",
            new.data = as.matrix(x = data.return[[i]])
          )
        } else {
          toRet[[names(x = data.return)[i]]] <- CreateAssayObject(counts = data.return[[i]], check.matrix = FALSE)
          toRet <- SetAssayData(
            object = toRet,
            assay = names(x = data.return)[i],
            slot = "data",
            new.data = log1p(x = as.matrix(x = data.return[[i]]))
          )
        }
        
      }
    }
    if (DefaultAssay(object = object) %in% names(x = data.return)) {
      DefaultAssay(object = toRet) <- DefaultAssay(object = object)
      if (slot[which(DefaultAssay(object = object) %in% names(x = data.return))[1]] != 'scale.data') {
        toRet <- ScaleData(object = toRet, verbose = verbose)
      }
    }
    if ('ident' %in% group.by) {
      first.cells <- sapply(X = 1:ncol(x = category.matrix),
                            FUN = function(x) {
                              return(category.matrix[,x, drop = FALSE ]@i[1] + 1)
                            }
      )
      Idents(object = toRet,
             cells = colnames(x = toRet)
             ) <- Idents(object = object)[first.cells]
    }
    return(toRet)
  } else {
    return(data.return)
  }
}

#' Match the case of character vectors
#'
#' @param search A vector of search terms
#' @param match A vector of characters whose case should be matched
#'
#' @return Values from search present in match with the case of match
#'
#' @export
#' @concept utilities
#'
#' @examples
#' data("pbmc_small")
#' cd_genes <- c('Cd79b', 'Cd19', 'Cd200')
#' CaseMatch(search = cd_genes, match = rownames(x = pbmc_small))
#'
CaseMatch <- function(search, match) {
  search.match <- sapply(
    X = search,
    FUN = function(s) {
      return(grep(
        pattern = paste0('^', s, '$'),
        x = match,
        ignore.case = TRUE,
        perl = TRUE,
        value = TRUE
      ))
    }
  )
  return(unlist(x = search.match))
}

#' Score cell cycle phases
#'
#' @param object A Seurat object
#' @param s.features A vector of features associated with S phase
#' @param g2m.features A vector of features associated with G2M phase
#' @param ctrl Number of control features selected from the same bin per
#' analyzed feature supplied to \code{\link{AddModuleScore}}.
#' Defaults to value equivalent to minimum number of features
#' present in 's.features' and 'g2m.features'.
#' @param set.ident If true, sets identity to phase assignments
#' Stashes old identities in 'old.ident'
#' @param ... Arguments to be passed to \code{\link{AddModuleScore}}
#'
#' @return A Seurat object with the following columns added to object meta data: S.Score, G2M.Score, and Phase
#'
#' @seealso \code{AddModuleScore}
#'
#' @export
#' @concept utilities
#'
#' @examples
#' \dontrun{
#' data("pbmc_small")
#' # pbmc_small doesn't have any cell-cycle genes
#' # To run CellCycleScoring, please use a dataset with cell-cycle genes
#' # An example is available at http://satijalab.org/seurat/cell_cycle_vignette.html
#' pbmc_small <- CellCycleScoring(
#'   object = pbmc_small,
#'   g2m.features = cc.genes$g2m.genes,
#'   s.features = cc.genes$s.genes
#' )
#' head(x = pbmc_small@meta.data)
#' }
#'
CellCycleScoring <- function(
  object,
  s.features,
  g2m.features,
  ctrl = NULL,
  set.ident = FALSE,
  ...
) {
  name <- 'Cell.Cycle'
  features <- list('S.Score' = s.features, 'G2M.Score' = g2m.features)
  if (is.null(x = ctrl)) {
    ctrl <- min(vapply(X = features, FUN = length, FUN.VALUE = numeric(length = 1)))
  }
  object.cc <- AddModuleScore(
    object = object,
    features = features,
    name = name,
    ctrl = ctrl,
    ...
  )
  cc.columns <- grep(pattern = name, x = colnames(x = object.cc[[]]), value = TRUE)
  cc.scores <- object.cc[[cc.columns]]
  rm(object.cc)
  CheckGC()
  assignments <- apply(
    X = cc.scores,
    MARGIN = 1,
    FUN = function(scores, first = 'S', second = 'G2M', null = 'G1') {
      if (all(scores < 0)) {
        return(null)
      } else {
        if (length(which(x = scores == max(scores))) > 1) {
          return('Undecided')
        } else {
          return(c(first, second)[which(x = scores == max(scores))])
        }
      }
    }
  )
  cc.scores <- merge(x = cc.scores, y = data.frame(assignments), by = 0)
  colnames(x = cc.scores) <- c('rownames', 'S.Score', 'G2M.Score', 'Phase')
  rownames(x = cc.scores) <- cc.scores$rownames
  cc.scores <- cc.scores[, c('S.Score', 'G2M.Score', 'Phase')]
  object[[colnames(x = cc.scores)]] <- cc.scores
  if (set.ident) {
    object[['old.ident']] <- Idents(object = object)
    Idents(object = object) <- 'Phase'
  }
  return(object)
}

#' Slim down a multi-species expression matrix, when only one species is primarily of interenst.
#'
#' Valuable for CITE-seq analyses, where we typically spike in rare populations of 'negative control' cells from a different species.
#'
#' @param object A UMI count matrix. Should contain rownames that start with
#' the ensuing arguments prefix.1 or prefix.2
#' @param prefix The prefix denoting rownames for the species of interest.
#' Default is "HUMAN_". These rownames will have this prefix removed in the returned matrix.
#' @param controls The prefix denoting rownames for the species of 'negative
#' control' cells. Default is "MOUSE_".
#' @param ncontrols How many of the most highly expressed (average) negative
#' control features (by default, 100 mouse genes), should be kept? All other
#' rownames starting with prefix.2 are discarded.
#'
#' @return A UMI count matrix. Rownames that started with \code{prefix} have this
#' prefix discarded. For rownames starting with \code{controls}, only the
#' \code{ncontrols} most highly expressed features are kept, and the
#' prefix is kept. All other rows are retained.
#'
#' @importFrom utils head
#' @importFrom Matrix rowSums
#'
#' @export
#' @concept utilities
#'
#' @examples
#' \dontrun{
#' cbmc.rna.collapsed <- CollapseSpeciesExpressionMatrix(cbmc.rna)
#' }
#'
CollapseSpeciesExpressionMatrix <- function(
  object,
  prefix = "HUMAN_",
  controls = "MOUSE_",
  ncontrols = 100
) {
  features <- grep(pattern = prefix, x = rownames(x = object), value = TRUE)
  controls <- grep(pattern = controls, x = rownames(x = object), value = TRUE)
  others <- setdiff(x = rownames(x = object), y = c(features, controls))
  controls <- rowSums(x = object[controls, ])
  controls <- names(x = head(
    x = sort(x = controls, decreasing = TRUE),
    n = ncontrols
  ))
  object <- object[c(features, controls, others), ]
  rownames(x = object) <- gsub(
    pattern = prefix,
    replacement = '',
    x = rownames(x = object)
  )
  return(object)
}

# Create an Annoy index
#
# @note Function exists because it's not exported from \pkg{uwot}
#
# @param name Distance metric name
# @param ndim Number of dimensions
#
# @return An nn index object
#
#' @importFrom methods new
#' @importFrom RcppAnnoy AnnoyAngular AnnoyManhattan AnnoyEuclidean AnnoyHamming
#
CreateAnn <- function(name, ndim) {
  return(switch(
    EXPR = name,
    cosine = new(Class = AnnoyAngular, ndim),
    manhattan = new(Class = AnnoyManhattan, ndim),
    euclidean = new(Class = AnnoyEuclidean, ndim),
    hamming = new(Class = AnnoyHamming, ndim),
    stop("BUG: unknown Annoy metric '", name, "'")
  ))
}

#' Run a custom distance function on an input data matrix
#'
#' @author Jean Fan
#'
#' @param my.mat A matrix to calculate distance on
#' @param my.function A function to calculate distance
#' @param ... Extra parameters to my.function
#'
#' @return A distance matrix
#'
#' @importFrom stats as.dist
#'
#' @export
#' @concept utilities
#'
#' @examples
#' data("pbmc_small")
#' # Define custom distance matrix
#' manhattan.distance <- function(x, y) return(sum(abs(x-y)))
#'
#' input.data <- GetAssayData(pbmc_small, assay.type = "RNA", slot = "scale.data")
#' cell.manhattan.dist <- CustomDistance(input.data, manhattan.distance)
#'
CustomDistance <- function(my.mat, my.function, ...) {
  CheckDots(..., fxns = my.function)
  n <- ncol(x = my.mat)
  mat <- matrix(data = 0, ncol = n, nrow = n)
  colnames(x = mat) <- rownames(x = mat) <- colnames(x = my.mat)
  for (i in 1:nrow(x = mat)) {
    for (j in 1:ncol(x = mat)) {
      mat[i,j] <- my.function(my.mat[, i], my.mat[, j], ...)
    }
  }
  return(as.dist(m = mat))
}

#' Calculate the mean of logged values
#'
#' Calculate mean of logged values in non-log space (return answer in log-space)
#'
#' @param x A vector of values
#' @param ... Other arguments (not used)
#'
#' @return Returns the mean in log-space
#'
#' @export
#' @concept utilities
#'
#' @examples
#' ExpMean(x = c(1, 2, 3))
#'
ExpMean <- function(x, ...) {
  if (inherits(x = x, what = 'AnyMatrix')) {
    return(apply(X = x, FUN = function(i) {log(x = mean(x = exp(x = i) - 1) + 1)}, MARGIN = 1))
  } else {
    return(log(x = mean(x = exp(x = x) - 1) + 1))
  }
}

#' Calculate the standard deviation of logged values
#'
#' Calculate SD of logged values in non-log space (return answer in log-space)
#'
#' @param x A vector of values
#'
#' @return Returns the standard deviation in log-space
#'
#' @importFrom stats sd
#'
#' @export
#' @concept utilities
#'
#' @examples
#' ExpSD(x = c(1, 2, 3))
#'
ExpSD <- function(x) {
  return(log1p(x = sd(x = expm1(x = x))))
}

#' Calculate the variance of logged values
#'
#' Calculate variance of logged values in non-log space (return answer in
#' log-space)
#'
#' @param x A vector of values
#'
#' @return Returns the variance in log-space
#'
#' @importFrom stats var
#'
#' @export
#' @concept utilities
#'
#' @examples
#' ExpVar(x = c(1, 2, 3))
#'
ExpVar <- function(x) {
  return(log1p(x = var(x = expm1(x = x))))
}

#' Scale and/or center matrix rowwise
#'
#' Performs row scaling and/or centering. Equivalent to using t(scale(t(mat)))
#' in R except in the case of NA values.
#'
#' @param mat A matrix
#' @param center a logical value indicating whether to center the rows
#' @param scale a logical value indicating whether to scale the rows
#' @param scale_max clip all values greater than scale_max to scale_max. Don't
#' clip if Inf.
#' @return Returns the center/scaled matrix
#'
#' @importFrom matrixStats rowMeans2 rowSds rowSums2
#'
#' @export
#' @concept utilities
#'
FastRowScale <- function(
  mat,
  center = TRUE,
  scale = TRUE,
  scale_max = 10
) {
  # inspired by https://www.r-bloggers.com/a-faster-scale-function/
  if (center) {
    rm <- rowMeans2(x = mat, na.rm = TRUE)
  }
  if (scale) {
    if (center) {
      rsd <- rowSds(mat, center = rm)
    } else {
      rsd <- sqrt(x = rowSums2(x = mat^2)/(ncol(x = mat) - 1))
    }
  }
  if (center) {
    mat <- mat - rm
  }
  if (scale) {
    mat <- mat / rsd
  }
  if (scale_max != Inf) {
    mat[mat > scale_max] <- scale_max
  }
  return(mat)
}


#' Get updated synonyms for gene symbols
#'
#' Find current gene symbols based on old or alias symbols using the gene
#' names database from the HUGO Gene Nomenclature Committee (HGNC)
#'
#' @details For each symbol passed, we query the HGNC gene names database for
#' current symbols that have the provided symbol as either an alias
#' (\code{alias_symbol}) or old (\code{prev_symbol}) symbol. All other queries
#' are \strong{not} supported.
#'
#' @note This function requires internet access
#'
#' @param symbols A vector of gene symbols
#' @param timeout Time to wait before canceling query in seconds
#' @param several.ok Allow several current gene symbols for each
#' provided symbol
#' @param search.types Type of query to perform:
#' \describe{
#'  \item{\dQuote{\code{alias_symbol}}}{Find alternate symbols for the genes
#'  described by \code{symbols}}
#'  \item{\dQuote{\code{prev_symbol}}}{Find new new symbols for the genes
#'  described by \code{symbols}}
#' }
#' This parameter accepts multiple options and short-hand options
#' (eg. \dQuote{\code{prev}} for \dQuote{\code{prev_symbol}})
#' @param verbose Show a progress bar depicting search progress
#' @param ... Extra parameters passed to \code{\link[httr]{GET}}
#'
#' @return \code{GeneSymbolThesarus}:, if \code{several.ok}, a named list
#' where each entry is the current symbol found for each symbol provided and
#' the names are the provided symbols. Otherwise, a named vector with the
#' same information.
#'
#' @source \url{https://www.genenames.org/} \url{https://www.genenames.org/help/rest/}
#'
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @importFrom httr GET accept_json timeout status_code content
#'
#' @rdname UpdateSymbolList
#' @name UpdateSymbolList
#'
#' @export
#' @concept utilities
#'
#' @seealso \code{\link[httr]{GET}}
#'
#' @examples
#' \dontrun{
#' GeneSybmolThesarus(symbols = c("FAM64A"))
#' }
#'
GeneSymbolThesarus <- function(
  symbols,
  timeout = 10,
  several.ok = FALSE,
  search.types = c('alias_symbol', 'prev_symbol'),
  verbose = TRUE,
  ...
) {
  db.url <- 'http://rest.genenames.org/fetch'
  # search.types <- c('alias_symbol', 'prev_symbol')
  search.types <- match.arg(arg = search.types, several.ok = TRUE)
  synonyms <- vector(mode = 'list', length = length(x = symbols))
  not.found <- vector(mode = 'logical', length = length(x = symbols))
  multiple.found <- vector(mode = 'logical', length = length(x = symbols))
  names(x = multiple.found) <- names(x = not.found) <- names(x = synonyms) <- symbols
  if (verbose) {
    pb <- txtProgressBar(max = length(x = symbols), style = 3, file = stderr())
  }
  for (symbol in symbols) {
    sym.syn <- character()
    for (type in search.types) {
      response <- GET(
        url = paste(db.url, type, symbol, sep = '/'),
        config = c(accept_json(), timeout(seconds = timeout)),
        ...
      )
      if (!identical(x = status_code(x = response), y = 200L)) {
        next
      }
      response <- content(x = response)
      if (response$response$numFound != 1) {
        if (response$response$numFound > 1) {
          warning(
            "Multiple hits found for ",
            symbol,
            " as ",
            type,
            ", skipping",
            call. = FALSE,
            immediate. = TRUE
          )
        }
        next
      }
      sym.syn <- c(sym.syn, response$response$docs[[1]]$symbol)
    }
    not.found[symbol] <- length(x = sym.syn) < 1
    multiple.found[symbol] <- length(x = sym.syn) > 1
    if (length(x = sym.syn) == 1 || (length(x = sym.syn) > 1 && several.ok)) {
      synonyms[[symbol]] <- sym.syn
    }
    if (verbose) {
      setTxtProgressBar(pb = pb, value = pb$getVal() + 1)
    }
  }
  if (verbose) {
    close(con = pb)
  }
  if (sum(not.found) > 0) {
    warning(
      "The following symbols had no synonyms: ",
      paste(names(x = which(x = not.found)), collapse = ', '),
      call. = FALSE,
      immediate. = TRUE
    )
  }
  if (sum(multiple.found) > 0) {
    msg <- paste(
      "The following symbols had multiple synonyms:",
      paste(names(x = which(x = multiple.found)), sep = ', ')
    )
    if (several.ok) {
      message(msg)
      message("Including anyways")
    } else {
      warning(msg, call. = FALSE, immediate. = TRUE)
    }
  }
  synonyms <- Filter(f = Negate(f = is.null), x = synonyms)
  if (!several.ok) {
    synonyms <- unlist(x = synonyms)
  }
  return(synonyms)
}

#' Compute the correlation of features broken down by groups with another
#' covariate
#'
#' @param object Seurat object
#' @param assay Assay to pull the data from
#' @param slot Slot in the assay to pull feature expression data from (counts,
#' data, or scale.data)
#' @param var Variable with which to correlate the features
#' @param group.assay Compute the gene groups based off the data in this assay.
#' @param min.cells Only compute for genes in at least this many cells
#' @param ngroups Number of groups to split into
#' @param do.plot Display the group correlation boxplot (via
#' \code{GroupCorrelationPlot})
#'
#' @return A Seurat object with the correlation stored in metafeatures
#'
#' @export
#' @concept utilities
#'
GroupCorrelation <- function(
  object,
  assay = NULL,
  slot = "scale.data",
  var = NULL,
  group.assay = NULL,
  min.cells = 5,
  ngroups = 6,
  do.plot = TRUE
) {
  assay <- assay %||% DefaultAssay(object = object)
  group.assay <- group.assay %||% assay
  var <- var %||% paste0("nCount_", group.assay)
  gene.grp <- GetFeatureGroups(
    object = object,
    assay = group.assay,
    min.cells = min.cells,
    ngroups = ngroups
  )
  data <- as.matrix(x = GetAssayData(object = object[[assay]], slot = slot))
  data <- data[rowMeans(x = data) != 0, ]
  grp.cors <- apply(
    X = data,
    MARGIN = 1,
    FUN = function(x) {
      cor(x = x, y = object[[var]])
    }
  )
  grp.cors <- grp.cors[names(x = gene.grp)]
  grp.cors <- as.data.frame(x = grp.cors[which(x = !is.na(x = grp.cors))])
  grp.cors$gene_grp <- gene.grp[rownames(x = grp.cors)]
  colnames(x = grp.cors) <- c(paste0(var, "_cor"), "feature.grp")
  object[[assay]][] <- grp.cors
  if (isTRUE(x = do.plot)) {
    print(GroupCorrelationPlot(
      object = object,
      assay = assay,
      feature.group = "feature.grp",
      cor = paste0(var, "_cor")
    ))
  }
  return(object)
}

#' Load the Annoy index file
#'
#' @param object Neighbor object
#' @param file Path to file with annoy index
#'
#' @return Returns the Neighbor object with the index stored
#' @export
#' @concept utilities
#'
LoadAnnoyIndex <- function(object, file){
  metric <- slot(object = object, name = "alg.info")$metric
  ndim <- slot(object = object, name = "alg.info")$ndim
  if (is.null(x = metric)) {
    stop("Provided Neighbor object wasn't generated with annoy")
  }
  annoy.idx <- CreateAnn(name = metric, ndim = ndim)
  annoy.idx$load(path.expand(path = file))
  Index(object = object) <- annoy.idx
  return(object)
}

#' Calculate the variance to mean ratio of logged values
#'
#' Calculate the variance to mean ratio (VMR) in non-logspace (return answer in
#' log-space)
#'
#' @param x A vector of values
#' @param ... Other arguments (not used)
#'
#' @return Returns the VMR in log-space
#'
#' @importFrom stats var
#'
#' @export
#' @concept utilities
#'
#' @examples
#' LogVMR(x = c(1, 2, 3))
#'
LogVMR <- function(x, ...) {
  if (inherits(x = x, what = 'AnyMatrix')) {
    return(apply(X = x, FUN = function(i) {log(x = var(x = exp(x = i) - 1) / mean(x = exp(x = i) - 1))}, MARGIN = 1))
  } else {
    return(log(x = var(x = exp(x = x) - 1) / mean(x = exp(x = x) - 1)))
  }
}

#' Aggregate expression of multiple features into a single feature
#'
#' Calculates relative contribution of each feature to each cell
#' for given set of features.
#'
#' @param object A Seurat object
#' @param features List of features to aggregate
#' @param meta.name Name of column in metadata to store metafeature
#' @param cells List of cells to use (default all cells)
#' @param assay Which assay to use
#' @param slot Which slot to take data from (default data)
#'
#' @return Returns a \code{Seurat} object with metafeature stored in objct metadata
#'
#' @importFrom Matrix rowSums colMeans
#'
#' @export
#' @concept utilities
#'
#' @examples
#' data("pbmc_small")
#' pbmc_small <- MetaFeature(
#'   object = pbmc_small,
#'   features = c("LTB", "EAF2"),
#'   meta.name = 'var.aggregate'
#' )
#' head(pbmc_small[[]])
#'
MetaFeature <- function(
  object,
  features,
  meta.name = 'metafeature',
  cells = NULL,
  assay = NULL,
  slot = 'data'
) {
  cells <- cells %||% colnames(x = object)
  assay <- assay %||% DefaultAssay(object = object)
  newmat <- GetAssayData(object = object, assay = assay, slot = slot)
  newmat <- newmat[features, cells]
  if (slot == 'scale.data') {
    newdata <- Matrix::colMeans(newmat)
  } else {
    rowtotals <- Matrix::rowSums(newmat)
    newmat <- newmat / rowtotals
    newdata <- Matrix::colMeans(newmat)
  }
  object[[meta.name]] <- newdata
  return(object)
}

#' Apply a ceiling and floor to all values in a matrix
#'
#' @param data Matrix or data frame
#' @param min all values below this min value will be replaced with min
#' @param max all values above this max value will be replaced with max
#' @return Returns matrix after performing these floor and ceil operations
#' @export
#' @concept utilities
#'
#' @examples
#' mat <- matrix(data = rbinom(n = 25, size = 20, prob = 0.2 ), nrow = 5)
#' mat
#' MinMax(data = mat, min = 4, max = 5)
#'
MinMax <- function(data, min, max) {
  data2 <- data
  data2[data2 > max] <- max
  data2[data2 < min] <- min
  return(data2)
}

#' Calculate the percentage of a vector above some threshold
#'
#' @param x Vector of values
#' @param threshold Threshold to use when calculating percentage
#'
#' @return Returns the percentage of \code{x} values above the given threshold
#'
#' @export
#' @concept utilities
#'
#' @examples
#' set.seed(42)
#' PercentAbove(sample(1:100, 10), 75)
#'
PercentAbove <- function(x, threshold) {
  return(length(x = x[x > threshold]) / length(x = x))
}

#' Calculate the percentage of all counts that belong to a given set of features
#'
#' This function enables you to easily calculate the percentage of all the counts belonging to a
#' subset of the possible features for each cell. This is useful when trying to compute the percentage
#' of transcripts that map to mitochondrial genes for example. The calculation here is simply the
#' column sum of the matrix present in the counts slot for features belonging to the set divided by
#' the column sum for all features times 100.
#'
#' @param object A Seurat object
#' @param pattern A regex pattern to match features against
#' @param features A defined feature set. If features provided, will ignore the pattern matching
#' @param col.name Name in meta.data column to assign. If this is not null, returns a Seurat object
#' with the proportion of the feature set stored in metadata.
#' @param assay Assay to use
#'
#' @return Returns a vector with the proportion of the feature set or if md.name is set, returns a
#' Seurat object with the proportion of the feature set stored in metadata.
#' @importFrom Matrix colSums
#' @export
#' @concept utilities
#'
#' @examples
#' data("pbmc_small")
#' # Calculate the proportion of transcripts mapping to mitochondrial genes
#' # NOTE: The pattern provided works for human gene names. You may need to adjust depending on your
#' # system of interest
#' pbmc_small[["percent.mt"]] <- PercentageFeatureSet(object = pbmc_small, pattern = "^MT-")
#'
PercentageFeatureSet <- function(
  object,
  pattern = NULL,
  features = NULL,
  col.name = NULL,
  assay = NULL
) {
  assay <- assay %||% DefaultAssay(object = object)
  if (!is.null(x = features) && !is.null(x = pattern)) {
    warn(message = "Both pattern and features provided. Pattern is being ignored.")
  }
  percent.featureset <- list()
  layers <- Layers(object = object, pattern = "counts")
  for (i in seq_along(along.with = layers)) {
    layer <- layers[i]
    features.layer <- features %||% grep(
      pattern = pattern, 
      x = rownames(x = object[[assay]][[layer]]), 
      value = TRUE)
    layer.data <- LayerData(object = object, 
                            assay = assay, 
                            layer = layer)
    layer.sums <- colSums(x = layer.data[features.layer, , drop = FALSE])
    layer.perc <- layer.sums / object[[]][colnames(layer.data), paste0("nCount_", assay)] * 100
    percent.featureset[[i]] <- layer.perc
  }
  percent.featureset <- unlist(percent.featureset)
  if (!is.null(x = col.name)) {
    object <- AddMetaData(object = object, metadata = percent.featureset, col.name = col.name)
    return(object)
  }
  return(percent.featureset)
}

# Pseudobulk feature expression by identity class
#
# Returns a representative expression value for each identity class
#
# @param object Seurat object
# @param method Whether to 'average' (default) or 'aggregate' expression levels
# @param assays Which assays to use. Default is all assays
# @param features Features to analyze. Default is all features in the assay
# @param return.seurat Whether to return the data as a Seurat object. Default is FALSE
# @param group.by Categories for grouping (e.g, ident, replicate, celltype); 'ident' by default
# @param add.ident (Deprecated) Place an additional label on each cell prior to pseudobulking
# (very useful if you want to observe cluster pseudobulk values, separated by replicate, for example)
# @param slot Slot(s) to use; if multiple slots are given, assumed to follow
# the order of 'assays' (if specified) or object's assays
# @param verbose Print messages and show progress bar
# @param ... Arguments to be passed to methods such as \code{\link{CreateSeuratObject}}
#
# @return Returns a matrix with genes as rows, identity classes as columns.
# If return.seurat is TRUE, returns an object of class \code{\link{Seurat}}.
#' @method PseudobulkExpression Assay
#' @export
#
#
PseudobulkExpression.Assay <- function(
  object,
  assay,
  category.matrix,
  features = NULL,
  slot = 'data',
  verbose = TRUE,
  ...
) {
    data.use <- GetAssayData(
      object = object, 
      slot = slot
    )
    features.to.avg <- features %||% rownames(x = data.use)
    if (IsMatrixEmpty(x = data.use)) {
      warning(
        "The ", slot, " slot for the ", assay,
        " assay is empty. Skipping assay.", immediate. = TRUE, call. = FALSE)
      return(NULL)
    }
    bad.features <- setdiff(x = features.to.avg, y = rownames(x = data.use))
    if (length(x = bad.features) > 0) {
      warning(
        "The following ", length(x = bad.features),
        " features were not found in the ", assay, " assay: ",
        paste(bad.features, collapse = ", "), call. = FALSE, immediate. = TRUE)
    }
    features.assay <- intersect(x = features.to.avg, y = rownames(x = data.use))
    if (length(x = features.assay) > 0) {
      data.use <- data.use[features.assay, ]
    } else {
      warning("None of the features specified were found in the ", assay,
              " assay.", call. = FALSE, immediate. = TRUE)
      return(NULL)
    }
    if (slot == 'data') {
      data.use <- expm1(x = data.use)
      if (any(data.use == Inf)) {
        warning("Exponentiation yielded infinite values. `data` may not be log-normed.")
      }
    }
    data.return <- data.use %*% category.matrix
   return(data.return)
 
 
}

#' @method PseudobulkExpression StdAssay
#' @export
#
#
PseudobulkExpression.StdAssay <- function(
  object,
  assay,
  category.matrix,
  features = NULL,
  slot = 'data',
  verbose = TRUE,
  ...
) {
  if (slot == 'data') {
    message("Assay5 will use arithmetic mean for data slot.")
  }
  layers.set <- Layers(object = object, search = slot)
  features.to.avg <- features %||% rownames(x = object)
  bad.features <- setdiff(x = features.to.avg, y = rownames(x = object))
  if (length(x = bad.features) > 0) {
    warning(
      "The following ", length(x = bad.features),
      " features were not found in the ", assay, " assay: ",
      paste(bad.features, collapse = ", "), call. = FALSE, immediate. = TRUE)
  }
  features.assay <- Reduce(
    f = intersect,
    x =  c(list(features.to.avg),
           lapply(X = layers.set, FUN = function(l) rownames(object[[l]]))
                )
         )
  if (length(x = features.assay) == 0) {
    warning("None of the features specified were found in the ", assay,
            " assay.", call. = FALSE, immediate. = TRUE)
    return(NULL)
  }
  data.return <- as.sparse(
    x = matrix(
      data = 0,
      nrow = length(x = features.assay),
      ncol = ncol(x = category.matrix)
      )
    )
  for (i in seq_along(layers.set)) {
    data.i <- LayerData(object = object,
                        layer = layers.set[i],
                        features = features.assay
                        )
    category.matrix.i <- category.matrix[colnames(x = data.i),]
    if (inherits(x = data.i, what = 'DelayedArray')) {
      data.return.i<- tcrossprod_DelayedAssay(x = data.i, y = t(category.matrix.i))
    } else {
      data.return.i <- as.sparse(x = data.i %*% category.matrix.i)
    }
    data.return <- data.return + data.return.i
  }
  if (slot == 'data') {
    data.return <- expm1(x = data.return)
  }
  return(data.return)
}

#' Regroup idents based on meta.data info
#'
#' For cells in each ident, set a new identity based on the most common value
#' of a specified metadata column.
#'
#' @param object Seurat object
#' @param metadata Name of metadata column
#' @return A Seurat object with the active idents regrouped
#'
#' @export
#' @concept utilities
#'
#' @examples
#' data("pbmc_small")
#' pbmc_small <- RegroupIdents(pbmc_small, metadata = "groups")
#'
RegroupIdents <- function(object, metadata) {
  for (ii in levels(x = object)) {
    ident.cells <- WhichCells(object = object, idents = ii)
    if (length(x = ident.cells) == 0) {
      next
    }
    new.ident <- names(x = which.max(x = table(object[[metadata]][ident.cells, ])))
    if (is.null(x = new.ident)) {
      stop("Cluster ", ii, " contains only cells with NA values in the '", metadata, "' metadata column.")
    }
    Idents(object = object, cells = ident.cells) <- new.ident
  }
  return(object)
}

#' Save the Annoy index
#'
#' @param object A Neighbor object with the annoy index stored
#' @param file Path to file to write index to
#'
#' @export
#' @concept utilities
#'
SaveAnnoyIndex <- function(
  object,
  file
) {
  index <- Index(object = object)
  if (is.null(x = index)) {
    stop("Index for provided Neighbor object is NULL")
  }
  index$save(path.expand(path = file))
}

#' Find the Quantile of Data
#'
#' Converts a quantile in character form to a number regarding some data.
#' String form for a quantile is represented as a number prefixed with
#' \dQuote{q}; for example, 10th quantile is \dQuote{q10} while 2nd quantile is
#' \dQuote{q2}. Will only take a quantile of non-zero data values
#'
#' @param cutoff The cutoff to turn into a quantile
#' @param data The data to turn find the quantile of
#'
#' @return The numerical representation of the quantile
#'
#' @importFrom stats quantile
#'
#' @export
#' @concept utilities
#'
#' @examples
#' set.seed(42)
#' SetQuantile('q10', sample(1:100, 10))
#'
SetQuantile <- function(cutoff, data) {
  if (grepl(pattern = '^q[0-9]{1,2}$', x = as.character(x = cutoff), perl = TRUE)) {
    this.quantile <- as.numeric(x = sub(
      pattern = 'q',
      replacement = '',
      x = as.character(x = cutoff)
    )) / 100
    data <- unlist(x = data)
    data <- data[data > 0]
    cutoff <- quantile(x = data, probs = this.quantile)
  }
  return(as.numeric(x = cutoff))
}

#' @rdname UpdateSymbolList
#'
#' @return \code{UpdateSymbolList}: \code{symbols} with updated symbols from
#' HGNC's gene names database
#'
#' @export
#' @concept utilities
#'
#' @examples
#' \dontrun{
#' UpdateSymbolList(symbols = cc.genes$s.genes)
#' }
#'
UpdateSymbolList <- function(
  symbols,
  timeout = 10,
  several.ok = FALSE,
  verbose = TRUE,
  ...
) {
  new.symbols <- suppressWarnings(expr = GeneSymbolThesarus(
    symbols = symbols,
    timeout = timeout,
    several.ok = several.ok,
    search.types = 'prev_symbol',
    verbose = verbose,
    ...
  ))
  if (length(x = new.symbols) < 1) {
    warning("No updated symbols found", call. = FALSE, immediate. = TRUE)
  } else {
    if (verbose) {
      message("Found updated symbols for ", length(x = new.symbols), " symbols")
      x <- sapply(X = new.symbols, FUN = paste, collapse = ', ')
      message(paste(names(x = x), x, sep = ' -> ', collapse = '\n'))
    }
    for (sym in names(x = new.symbols)) {
      index <- which(x = symbols == sym)
      symbols <- append(
        x = symbols[-index],
        values = new.symbols[[sym]],
        after = index - 1
      )
    }
  }
  return(symbols)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Methods for Seurat-defined generics
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Methods for R-defined generics
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @inheritParams base::as.data.frame
#'
#' @return \code{as.data.frame.Matrix}: A data frame representation of the S4 Matrix
#'
#' @importFrom Matrix as.matrix
#'
#' @rdname as.sparse
#' @concept utilities
#' @export
#' @method as.data.frame Matrix
#'
as.data.frame.Matrix <- function(
  x,
  row.names = NULL,
  optional = FALSE,
  ...,
  stringsAsFactors = getOption(x = "stringsAsFactors", default = FALSE)
) {
  return(as.data.frame(
    x = as.matrix(x = x),
    row.names = row.names,
    optional = optional,
    stringsAsFactors = stringsAsFactors,
    ...
  ))
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Internal
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Create Abbreviations
#'
#' @param x A character vector
#' @param digits Include digits in the abbreviation
#'
#' @return Abbreviated versions of \code{x}
#'
#' @keywords internal
#'
#' @examples
#' .Abbrv(c('HelloWorld, 'LetsGo3', 'tomato'))
#' .Abbrv(c('HelloWorld, 'LetsGo3', 'tomato'), digits = FALSE)
#' .Abbrv('Wow3', digits = FALSE)
#'
#' @noRd
#'
.Abbrv <- function(x, digits = TRUE) {
  pattern <- ifelse(test = isTRUE(x = digits), yes = '[A-Z0-9]+', no = '[A-Z]+')
  y <- vapply(
    X = regmatches(x = x, m = gregexec(pattern = pattern, text = x)),
    FUN = paste,
    FUN.VALUE = character(length = 1L),
    collapse = ''
  )
  na <- nchar(x = y) <= 1L
  y[na] <- x[na]
  return(tolower(x = y))
}

.AsList <- function(x) {
  x <- as.list(x = x)
  return(sapply(
    X = unique(x = names(x = x)),
    FUN = function(i) {
      return(unlist(
        x = x[which(x = names(x = x) == i)],
        recursive = FALSE,
        use.names = FALSE
      ))
    },
    simplify = FALSE,
    USE.NAMES = TRUE
  ))
}

#' @importFrom ggplot2 cut_number
#'
.Cut <- function(min, max, n) {
  breaks <- levels(x = cut_number(x = c(min, max), n = n))
  breaks <- gsub(pattern = '.*,', replacement = '', x = breaks)
  breaks <- gsub(pattern = ']$', replacement = '', x = breaks)
  as.numeric(x = breaks)
}

.FindE <- function(x) {
  x <- as.character(x = x)
  if (grepl(pattern = 'e', x = x)) {
    return(as.integer(x = gsub(pattern = '.*e', replacement = '', x = x)))
  } else if (grepl(pattern = '^0\\.', x = x)) {
    x <- unlist(x = strsplit(
      x = gsub(pattern = '.*\\.', replacement = '', x = x),
      split = ''
    ))
    idx <- which(x = x != '0')
    return(-idx)
  }
  stop("Invalid format")
}

#' @importFrom SeuratObject Boundaries
#'
.BoundariesByImage <- function(object, fov, boundaries) {
  if (!is.list(x = boundaries)) {
    if (is.null(x = names(x = boundaries))) {
      boundaries <- rep_len(x = list(boundaries), length.out = length(x = fov))
      names(x = boundaries) <- fov
    } else {
      boundaries <- .AsList(x = boundaries)
    }
  }
  if (any(!nchar(x = names(x = boundaries)))) {
    missing <- setdiff(x = fov, y = names(x = boundaries))
    idx <- which(x = !nchar(x = names(x = boundaries)))
    boundaries <- c(
      boundaries[intersect(x = names(x = boundaries), y = fov)],
      rep_len(x = boundaries[idx], length.out = length(x = missing))
    )
    names(x = boundaries)[!nchar(x = names(x = boundaries))] <- missing
  }
  if (any(!fov %in% names(x = boundaries))) {
    for (i in setdiff(x = fov, y = names(x = boundaries))) {
      boundaries[[i]] <- Boundaries(object = object[[i]])[1L]
    }
  }
  fov <- union(x = fov, y = names(x = boundaries))
  if (length(x = boundaries) != length(x = fov)) {
    fov <- intersect(x = fov, y = names(x = boundaries))
  }
  boundaries <- boundaries[fov]
  for (i in fov) {
    boundaries[[i]] <- Filter(
      f = function(x) {
        return(x %in% Boundaries(object = object[[i]]) || is_na(x = x))
      },
      x = boundaries[[i]]
    )
  }
  boundaries <- Filter(f = length, x = boundaries)
  return(boundaries)
}

# Generate chunk points
#
# @param dsize How big is the data being chunked
# @param csize How big should each chunk be
#
# @return A matrix where each column is a chunk, row 1 is start points, row 2 is end points
#
ChunkPoints <- function(dsize, csize) {
  return(vapply(
    X = 1L:ceiling(x = dsize / csize),
    FUN = function(i) {
      return(c(
        start = (csize * (i - 1L)) + 1L,
        end = min(csize * i, dsize)
      ))
    },
    FUN.VALUE = numeric(length = 2L)
  ))
}

# L2 normalize the columns (or rows) of a given matrix
# @param mat Matrix to cosine normalize
# @param MARGIN Perform normalization over rows (1) or columns (2)
#
#
# @return returns l2-normalized matrix
#
#
L2Norm <- function(mat, MARGIN = 1){
  normalized <- Sweep(
    x = mat,
    MARGIN = MARGIN,
    STATS = apply(
      X = mat,
      MARGIN = MARGIN,
      FUN = function(x){
        sqrt(x = sum(x ^ 2))
      }
    ),
    FUN = "/"
  )
  normalized[!is.finite(x = normalized)] <- 0
  return(normalized)
}

# Check the use of ...
#
# @param ... Arguments passed to a function that fall under ...
# @param fxns A list/vector of functions or function names
#
# @return ...
#
# @importFrom utils argsAnywhere getAnywhere
#' @importFrom utils isS3stdGeneric methods argsAnywhere isS3method
#
# @examples
#
CheckDots <- function(..., fxns = NULL) {
  args.names <- names(x = list(...))
  if (length(x = list(...)) == 0) {
    return(invisible(x = NULL))
  }
  if (is.null(x = args.names)) {
    stop("No named arguments passed")
  }
  if (length(x = fxns) == 1) {
    fxns <- list(fxns)
  }
  for (f in fxns) {
    if (!(is.character(x = f) || is.function(x = f))) {
      stop("CheckDots only works on characters or functions, not ", class(x = f))
    }
  }
  fxn.args <- suppressWarnings(expr = sapply(
    X = fxns,
    FUN = function(x) {
      x <- tryCatch(
        expr = if (isS3stdGeneric(f = x)) {
          as.character(x = methods(generic.function = x))
        } else {
          x
        },
        error = function(...) {
          return(x)
        }
      )
      x <- if (is.character(x = x)) {
        sapply(X = x, FUN = argsAnywhere, simplify = FALSE, USE.NAMES = TRUE)
      } else if (length(x = x) <= 1) {
        list(x)
      }
      return(sapply(
        X = x,
        FUN = function(f) {
          return(names(x = formals(fun = f)))
        },
        simplify = FALSE,
        USE.NAMES = TRUE
      ))
    },
    simplify = FALSE,
    USE.NAMES = TRUE
  ))
  fxn.args <- unlist(x = fxn.args, recursive = FALSE)
  fxn.null <- vapply(X = fxn.args, FUN = is.null, FUN.VALUE = logical(length = 1L))
  if (all(fxn.null) && !is.null(x = fxns)) {
    stop("None of the functions passed could be found")
  } else if (any(fxn.null)) {
    warning(
      "The following functions passed could not be found: ",
      paste(names(x = which(x = fxn.null)), collapse = ', '),
      call. = FALSE,
      immediate. = TRUE
    )
    fxn.args <- Filter(f = Negate(f = is.null), x = fxn.args)
  }
  dfxns <- vector(mode = 'logical', length = length(x = fxn.args))
  names(x = dfxns) <- names(x = fxn.args)
  for (i in 1:length(x = fxn.args)) {
    dfxns[i] <- any(grepl(pattern = '...', x = fxn.args[[i]], fixed = TRUE))
  }
  if (any(dfxns)) {
    dfxns <- names(x = which(x = dfxns))
    if (any(nchar(x = dfxns) > 0)) {
      fx <- vapply(
        X = Filter(f = nchar, x = dfxns),
        FUN = function(x) {
          if (isS3method(method = x)) {
            x <- unlist(x = strsplit(x = x, split = '\\.'))
            x <- x[length(x = x) - 1L]
          }
          return(x)
        },
        FUN.VALUE = character(length = 1L)
      )
      message(
        "The following functions and any applicable methods accept the dots: ",
        paste(unique(x = fx), collapse = ', ')
      )
      if (any(nchar(x = dfxns) < 1)) {
        message(
          "In addition, there is/are ",
          length(x = Filter(f = Negate(f = nchar), x = dfxns)),
          " other function(s) that accept(s) the dots"
        )
      }
    } else {
      message("There is/are ", length(x = dfxns), 'function(s) that accept(s) the dots')
    }
  } else {
    unused <- Filter(
      f = function(x) {
        return(!x %in% unlist(x = fxn.args))
      },
      x = args.names
    )
    if (length(x = unused) > 0) {
      msg <- paste0(
        "The following arguments are not used: ",
        paste(unused, collapse = ', ')
      )
      switch(
        EXPR = getOption(x = "Seurat.checkdots"),
        "warn" = warning(msg, call. = FALSE, immediate. = TRUE),
        "stop" = stop(msg),
        "silent" = NULL,
        stop("Invalid Seurat.checkdots option. Please choose one of warn, stop, silent")
      )
      unused.hints <- sapply(X = unused, FUN = OldParamHints)
      names(x = unused.hints) <- unused
      unused.hints <- na.omit(object = unused.hints)
      if (length(x = unused.hints) > 0) {
        message(
          "Suggested parameter: ",
          paste(unused.hints, "instead of", names(x = unused.hints), collapse = '; '),
          "\n"
        )
      }
    }
  }
}

# Call gc() to perform garbage collection
#
CheckGC <- function() {
  if (getOption(x = "Seurat.memsafe")) {
    gc(verbose = FALSE)
  }
}

# Check a list of objects for duplicate cell names
#
# @param object.list List of Seurat objects
# @param verbose Print message about renaming
# @param stop Error out if any duplicate names exist
#
# @return Returns list of objects with duplicate cells renamed to be unique
#
# @keywords internal
#
# @noRd
#
CheckDuplicateCellNames <- function(object.list, verbose = TRUE, stop = FALSE) {
  cell.names <- unlist(x = lapply(X = object.list, FUN = colnames))
  if (any(duplicated(x = cell.names))) {
    if (stop) {
      stop("Duplicate cell names present across objects provided.")
    }
    if (verbose) {
      warning("Some cell names are duplicated across objects provided. Renaming to enforce unique cell names.")
    }
    object.list <- lapply(
      X = 1:length(x = object.list),
      FUN = function(x) {
        return(RenameCells(
          object = object.list[[x]],
          new.names = paste0(Cells(x = object.list[[x]]), "_", x)
        ))
      }
    )
  }
  return(object.list)
}


# Create an empty dummy assay to replace existing assay
#' @importFrom Matrix sparseMatrix
CreateDummyAssay <- function(assay) {
  cm <- sparseMatrix(
    i = {},
    j = {},
    dims = c(nrow(x = assay), ncol(x = assay))
  )
  cm <- as.sparse(x = cm)
  rownames(x = cm) <- rownames(x = assay)
  colnames(x = cm) <- colnames(x = assay)
  return(CreateAssayObject(
    counts = cm,
    check.matrix = FALSE
  ))
}

# Extract delimiter information from a string.
#
# Parses a string (usually a cell name) and extracts fields based on a delimiter
#
# @param string String to parse.
# @param field Integer(s) indicating which field(s) to extract. Can be a vector multiple numbers.
# @param delim Delimiter to use, set to underscore by default.
#
# @return A new string, that parses out the requested fields, and (if multiple), rejoins them with the same delimiter
#
# @export
#
# @examples
# ExtractField(string = 'Hello World', field = 1, delim = '_')
#
ExtractField <- function(string, field = 1, delim = "_") {
  fields <- as.numeric(x = unlist(x = strsplit(x = as.character(x = field), split = ",")))
  if (length(x = fields) == 1) {
    return(strsplit(x = string, split = delim)[[1]][field])
  }
  return(paste(strsplit(x = string, split = delim)[[1]][fields], collapse = delim))
}

# Resize GenomicRanges upstream and or downstream
# from https://support.bioconductor.org/p/78652/
#
Extend <- function(x, upstream = 0, downstream = 0) {
  if (any(GenomicRanges::strand(x = x) == "*")) {
    warning("'*' ranges were treated as '+'")
  }
  on_plus <- GenomicRanges::strand(x = x) == "+" | GenomicRanges::strand(x = x) == "*"
  new_start <- GenomicRanges::start(x = x) - ifelse(test = on_plus, yes = upstream, no = downstream)
  new_end <- GenomicRanges::end(x = x) + ifelse(test = on_plus, yes = downstream, no = upstream)
  IRanges::ranges(x = x) <- IRanges::IRanges(start = new_start, end = new_end)
  x <- GenomicRanges::trim(x = x)
  return(x)
}

# Interleave vectors together
#
# @param ... Vectors to be interleaved
#
# @return A vector with the values from each vector in ... interleaved
#
Interleave <- function(...) {
  return(as.vector(x = t(x = as.data.frame(x = list(...)))))
}

# Check if a matrix is empty
#
# Takes a matrix and asks if it's empty (either 0x0 or 1x1 with a value of NA)
#
# @param x A matrix
#
# @return Whether or not \code{x} is empty
#
IsMatrixEmpty <- function(x) {
  matrix.dims <- dim(x = x)
  matrix.na <- all(matrix.dims == 1) && all(is.na(x = x))
  return(all(matrix.dims == 0) || matrix.na)
}

# Check if externalptr is null
# From https://stackoverflow.com/questions/26666614/how-do-i-check-if-an-externalptr-is-null-from-within-r
#
is.null.externalptr <- function(pointer) {
  stopifnot(is(pointer, "externalptr"))
  .Call("isnull", pointer)
}

# Check whether an assay has been processed by sctransform
#
# @param assay assay to check
#
# @return Boolean
#
IsSCT <- function(assay) {
  if (is.list(x = assay)) {
    sct.check <- lapply(X = assay, FUN = function(x) {
      return(!is.null(x = Misc(object = x, slot = 'vst.out')) | !is.null(x = Misc(object = x, slot = 'vst.set')) | inherits(x = x, what = "SCTAssay"))
    })
    return(unlist(x = sct.check))
  }
  return(!is.null(x = Misc(object = assay, slot = 'vst.out')) | !is.null(x = Misc(object = assay, slot = 'vst.set')) | inherits(x = assay, what = "SCTAssay"))
}

# Check whether a vst.out is from sctransform
#
# @param vst.out a sct model from sctransform
#
# @return Boolean
#
IsVSTout <- function(vst.out) {
  vst.element <- c("model_str", "model_pars_fit", "cell_attr" )
   vst.check <- setdiff(x = vst.element, y = names(x = vst.out))
   if (length(x = setdiff(x = vst.element, y = names(x = vst.out))) == 0) {
     vst.check <- TRUE
   } else {
     vst.check <- FALSE
   }
  return(vst.check)
}

# Calculate euclidean distance the x and y,
# and subtract the nearest neighbors of x distance to keep local connectivity
# It is used in FindModalityWeights to calculate the with and cross modality distance
impute_dist <- function(x, y, nearest.dist) {
  dist <- sqrt(x = rowSums(x = (x - y)**2)) - nearest.dist
  dist <- ReLu(x = dist)
  return(dist)
}

# Check the length of components of a list
#
# @param values A list whose components should be checked
# @param cutoff A minimum value to check for
#
# @return a vector of logicals
#
LengthCheck <- function(values, cutoff = 0) {
  return(vapply(
    X = values,
    FUN = function(x) {
      return(length(x = x) > cutoff)
    },
    FUN.VALUE = logical(1)
  ))
}

# Function to map values in a vector `v` as defined in `from`` to the values
# defined in `to`.
#
# @param v     vector of values to map
# @param from  vector of original values
# @param to    vector of values to map original values to (should be of equal
#              length as from)
# @return      returns vector of mapped values
#
MapVals <- function(v, from, to) {
  if (length(x = from) != length(x = to)) {
    stop("from and to vectors are not the equal length.")
  }
  vals.to.match <- match(x = v, table = from)
  vals.to.match.idx  <- !is.na(x = vals.to.match)
  v[vals.to.match.idx] <- to[vals.to.match[vals.to.match.idx]]
  return(v)
}

# Independently shuffle values within each row of a matrix
#
# Creates a matrix where correlation structure has been removed, but overall values are the same
#
# @param x Matrix to shuffle
#
# @return Returns a scrambled matrix, where each row is shuffled independently
#
#' @importFrom stats runif
#
# @export
#
# @examples
# mat <- matrix(data = rbinom(n = 25, size = 20, prob = 0.2 ), nrow = 5)
# mat
# MatrixRowShuffle(x = mat)
#
MatrixRowShuffle <- function(x) {
  x2 <- x
  x2 <- t(x = x)
  ind <- order(c(col(x = x2)), runif(n = length(x = x2)))
  x2 <- matrix(
    data = x2[ind],
    nrow = nrow(x = x),
    ncol = ncol(x = x),
    byrow = TRUE
  )
  return(x2)
}

# Reverse the vector x and return the value at the Nth index. If N is larger
# than the length of the vector, return the last value in the reversed vector.
#
# @param x vector of interest
# @param N index in reversed vector
#
# @return returns element at given index
#
MaxN <- function(x, N = 2){
  len <- length(x)
  if (N > len) {
    warning('N greater than length(x).  Setting N=length(x)')
    N <- length(x)
  }
  sort(x, partial = len - N + 1)[len - N + 1]
}

# Given a range from cut, compute the mean
#
# @x range from cut as a string (e.g. (10, 20] )
# @return returns a numeric with the mean of the range
#
MeanRange <- function(x) {
  left <- gsub(pattern = "\\]", replacement = "", x = sub(pattern = "\\([[:digit:]\\.e+]*,", x = x, replacement = ""))
  right <- gsub(pattern = "\\(", replacement = "", x = sub(pattern = ",[[:digit:]\\.e+]*]", x = x, replacement = ""))
  return(mean(c(as.numeric(x = left), as.numeric(x = right))))
}

# Melt a data frame
#
# @param x A data frame
#
# @return A molten data frame
#
Melt <- function(x) {
  if (!is.data.frame(x = x)) {
    x <- as.data.frame(x = x)
  }
  return(data.frame(
    rows = rep.int(x = rownames(x = x), times = ncol(x = x)),
    cols = unlist(x = lapply(X = colnames(x = x), FUN = rep.int, times = nrow(x = x))),
    vals = unlist(x = x, use.names = FALSE)
  ))
}

# Modify parameters in calling environment
#
# Used exclusively for helper parameter validation functions
#
# @param param name of parameter to change
# @param value new value for parameter
#
ModifyParam <- function(param, value) {
  # modify in original function environment
  env1 <- sys.frame(which = length(x = sys.frames()) - 2)
  env1[[param]] <- value
  # also modify in validation function environment
  env2 <- sys.frame(which = length(x = sys.frames()) - 1)
  env2[[param]] <- value
}

# Give hints for old parameters and their newer counterparts
#
# This is a non-exhaustive list. If your function isn't working properly based
# on the parameters you give it, please read the documentation for your function
#
# @param param A vector of parameters to get hints for
#
# @return Parameter hints for the specified parameters
#
OldParamHints <- function(param) {
  param.conversion <- c(
    'raw.data' = 'counts',
    'min.genes' = 'min.features',
    'features.plot' = 'features',
    'pc.genes' = 'features',
    'do.print' = 'verbose',
    'genes.print' = 'nfeatures.print',
    'pcs.print' = 'ndims.print',
    'pcs.use' = 'dims',
    'reduction.use' = 'reduction',
    'cells.use' = 'cells',
    'do.balanced' = 'balanced',
    'display.progress' = 'verbose',
    'print.output' = 'verbose',
    'dims.use' = 'dims',
    'reduction.type' = 'reduction',
    'y.log' = 'log',
    'cols.use' = 'cols',
    'assay.use' = 'assay'
  )
  return(param.conversion[param])
}

# Check if a web resource is available
#
# @param url A URL
# @param strict Perform a strict web availability test
# @param seconds Timeout in seconds
#
# @return \code{TRUE} if \url{is available} otherwise \code{FALSE}
#
#' @importFrom httr GET status_code timeout
#
# @keywords internal
#
Online <- function(url, strict = FALSE, seconds = 5L) {
  if (isTRUE(x = strict)) {
    code <- 200L
    comp <- identical
  } else {
    code <- 404L
    comp <- Negate(f = identical)
  }
  request <- tryCatch(
    expr = GET(url = url, timeout(seconds = seconds)),
    error = function(err) {
      code <- if (grepl(pattern = '^Timeout was reached', x = err$message)) {
        408L
      } else {
        404L
      }
      return(code)
    }
  )
  return(comp(x = status_code(x = request), y = code))
}

# Parenting parameters from one environment to the next
#
# This function allows one to modify a parameter in a parent environment
# The primary use of this is to ensure logging functions store correct parameters
# if they've been modified by a child function or method
#
# @param parent.find Regex pattern of name of parent function call to modify.
# For example, this can be the class name for a method that was dispatched previously
# @param ... Parameter names and values to parent; both name and value must be supplied
# in the standard \code{name = value} format; any number of name/value pairs can be specified
#
# @return No return, modifies parent environment directly
#
# @examples
# Parenting(parent.find = 'Seurat', features = features[features > 7])
#
Parenting <- function(parent.find = 'Seurat', ...) {
  calls <- as.character(x = sys.calls())
  calls <- lapply(
    X = strsplit(x = calls, split = '(', fixed = TRUE),
    FUN = '[',
    1
  )
  parent.index <- grep(pattern = parent.find, x = calls)
  if (length(x = parent.index) != 1) {
    warning(
      "Cannot find a parent environment called ",
      parent.find,
      immediate. = TRUE,
      call. = FALSE
    )
  } else {
    to.parent <- list(...)
    if (length(x = to.parent) == 0) {
      warning("Nothing to parent", immediate. = TRUE, call. = FALSE)
    } else if (is.null(x = names(x = to.parent))) {
      stop("All input must be in a key = value pair")
    } else if (length(x = Filter(f = nchar, x = names(x = to.parent))) != length(x = to.parent)) {
      stop("All inputs must be named")
    } else {
      parent.environ <- sys.frame(which = parent.index)
      for (i in 1:length(x = to.parent)) {
        parent.environ[[names(x = to.parent)[i]]] <- to.parent[[i]]
      }
    }
  }
}

# Generate a random name
#
# Make a name from randomly sampled lowercase letters,
# pasted together with no spaces or other characters
#
# @param length How long should the name be
# @param ... Extra parameters passed to sample
#
# @return A character with nchar == length of randomly sampled letters
#
# @seealso \code{\link{sample}}
#
RandomName <- function(length = 5L, ...) {
  CheckDots(..., fxns = 'sample')
  return(paste(sample(x = letters, size = length, ...), collapse = ''))
}


# Rectified linear units function. Calculate positive part of its argument
# The input can be a vector and a matrix
ReLu <- function(x) {
  x[x < 0] <- 0
  return(x)
}

# Remove the last field from a string
#
# Parses a string (usually a cell name) and removes the last field based on a delimter
#
# @param string String to parse
# @param delim Delimiter to use, set to underscore by default.
#
# @return A new string sans the last field
#
RemoveLastField <- function(string, delim = "_") {
  ss <- strsplit(x = string, split = delim)[[1]]
  if (length(x = ss) == 1) {
    return(string)
  } else {
    return(paste(ss[1:(length(x = ss)-1)], collapse = delim))
  }
}

# Calculate row mean of a sparse matrix
# @param mat sparse matrix
# @return A vector of row mean
#
RowMeanSparse <- function(mat) {
  mat <- RowSparseCheck(mat = mat)
  output <- row_mean_dgcmatrix(
    x = slot(object = mat, name = "x"),
    i = slot(object = mat, name = "i"),
    rows = nrow(x = mat),
    cols = ncol(x = mat)
    )
  names(x = output) <- rownames(x = mat)
  return(output)
}

# Calculate row sum of a sparse matrix
#
# @param mat sparse matrix
# @return A vector of row sum
#
RowSumSparse <- function(mat) {
  mat <- RowSparseCheck(mat = mat)
  output <- row_sum_dgcmatrix(
    x = slot(object = mat, name = "x"),
    i = slot(object = mat, name = "i"),
    rows = nrow(x = mat),
    cols = ncol(x = mat)
  )
  names(x = output) <- rownames(x = mat)
  return(output)
}

# Calculate row variance of a sparse matrix
#
# @param mat sparse matrix
# @return A vector of row variance
#
RowVarSparse <- function(mat) {
  mat <- RowSparseCheck(mat = mat)
  output <- row_var_dgcmatrix(
    x = slot(object = mat, name = "x"),
    i = slot(object = mat, name = "i"),
    rows = nrow(x = mat),
    cols = ncol(x = mat)
  )
  names(x = output) <- rownames(x = mat)
  return(output)
}

# Check if the input matrix is dgCMatrix
#
# @param mat sparse matrix
# @return A dgCMatrix
#
RowSparseCheck <- function(mat) {
  if (!inherits(x = mat, what = "sparseMatrix")) {
    stop("Input should be sparse matrix")
  } else if (!is(object = mat, class2 = "dgCMatrix")) {
    warning("Input matrix is converted to dgCMatrix.")
    mat <- as.sparse(x = mat)
  }
  return(mat)
}

# Sweep out array summaries
#
# Reimplmentation of \code{\link[base]{sweep}} to maintain compatability with
# both R 3.X and 4.X
#
# @inheritParams base::sweep
# @param x an array.
#
# @seealso \code{\link[base]{sweep}}
#
Sweep <- function(x, MARGIN, STATS, FUN = '-', check.margin = TRUE, ...) {
  if (any(grepl(pattern = 'X', x = names(x = formals(fun = sweep))))) {
    return(sweep(
      X = x,
      MARGIN = MARGIN,
      STATS = STATS,
      FUN = FUN,
      check.margin = check.margin,
      ...
    ))
  } else {
    return(sweep(
      x = x,
      MARGIN = MARGIN,
      STATS = STATS,
      FUN = FUN,
      check.margin = check.margin,
      ...
    ))
  }
}

# Get program paths in a system-agnostic way
#
# @param progs A vector of program names
# @param error Throw an error if any programs are not found
# @param add.exe Add '.exe' extension to program names that don't have it
#
# @return A named vector of program paths; missing programs are returned as
# \code{NA} if \code{error = FALSE}
#
#' @importFrom tools file_ext
#
SysExec <- function(
  progs,
  error = ifelse(test = length(x = progs) == 1, yes = TRUE, no = FALSE),
  add.exe = .Platform$OS.type == 'windows'
) {
  cmd <- ifelse(
    test = .Platform$OS.type == 'windows',
    yes = 'where.exe',
    no = 'which'
  )
  if (add.exe) {
    missing.exe <- file_ext(x = progs) != 'exe'
    progs[missing.exe] <- paste0(progs[missing.exe], '.exe')
  }
  paths <- sapply(
    X = progs,
    FUN = function(x) {
      return(tryCatch(
        expr = system2(command = cmd, args = x, stdout = TRUE)[1],
        warning = function(...) {
          return(NA_character_)
        }
      ))
    }
  )
  if (error && any(is.na(x = paths))) {
    stop(
      "Could not find the following programs: ",
      paste(names(x = paths[is.na(x = paths)]), collapse = ', '),
      call. = FALSE
    )
  }
  return(paths)
}

# Try to convert x to numeric, if NA's introduced return x as is
#
ToNumeric <- function(x){
  # check for x:y range
  if (is.numeric(x = x)) {
    return(x)
  }
  if (length(x = unlist(x = strsplit(x = x, split = ":"))) == 2) {
    num <- unlist(x = strsplit(x = x, split = ":"))
    return(num[1]:num[2])
  }
  num <- suppressWarnings(expr = as.numeric(x = x))
  if (!is.na(x = num)) {
    return(num)
  }
  return(x)
}



# cross product from delayed array
#
crossprod_DelayedAssay <- function(x, y, block.size = 1e8) {
  # perform t(x) %*% y in blocks for y
  if (!inherits(x = y, 'DelayedMatrix')) {
    stop('y should a DelayedMatrix')
  }
  if (nrow(x) != nrow(y)) {
    stop('row of x and y should be the same')
  }
  sparse <- DelayedArray::is_sparse(x = y)
  suppressMessages(setAutoBlockSize(size = block.size))
  cells.grid <- DelayedArray::colAutoGrid(x = y)
  product.list <- list()
  for (i in seq_len(length.out = length(x = cells.grid))) {
    vp <- cells.grid[[i]]
    block <- DelayedArray::read_block(x = y, viewport = vp, as.sparse = sparse)
    if (sparse) {
      block <- as(object = block, Class = 'dgCMatrix')
    } else {
      block <- as(object = block, Class = 'Matrix')
    }
    product.list[[i]] <- as.matrix(t(x) %*% block)
  }
  product.mat <- matrix(data = unlist(product.list), nrow = ncol(x) , ncol = ncol(y))
  colnames(product.mat) <- colnames(y)
  rownames(product.mat) <- colnames(x)
  return(product.mat)
}


# cross product from BPCells
#
crossprod_BPCells <- function(x, y) {
  # perform t(x) %*% y, y is from BPCells
  product.mat <- t(x) %*% y
  colnames(product.mat) <- colnames(y)
  rownames(product.mat) <- colnames(x)
  return(product.mat)
}

# transpose cross product from delayed array
#
tcrossprod_DelayedAssay <- function(x, y, block.size = 1e8) {
  # perform  x  %*% t(y) in blocks for x
  if (!inherits(x = x, 'DelayedMatrix')) {
    stop('y should a DelayedMatrix')
  }
  if (ncol(x) != ncol(y)) {
    stop('column of x and y should be the same')
  }
  sparse <- DelayedArray::is_sparse(x = x)
  suppressMessages(setAutoBlockSize(size = block.size))
  cells.grid <- DelayedArray::colAutoGrid(x = x)
  product.list <- list()
  for (i in seq_len(length.out = length(x = cells.grid))) {
    vp <- cells.grid[[i]]
    vp.range <- vp@ranges[2]@start : (vp@ranges[2]@start + vp@ranges[2]@width - 1)
    block <- DelayedArray::read_block(x = x, viewport = vp, as.sparse = sparse)
    if (sparse) {
      block <- as(object = block, Class = 'dgCMatrix')
    } else {
      block <- as(object = block, Class = 'Matrix')
    }
    product.list[[i]] <- as.matrix( block %*% t(y[,vp.range]))
  }
  product.mat <-  Reduce(f = '+', product.list)
  colnames(product.mat) <- rownames(y)
  rownames(product.mat) <- rownames(x)
  return(product.mat)
}

# cross product row norm from delayed array
#
crossprodNorm_DelayedAssay <- function(x, y, block.size = 1e8) {
  # perform t(x) %*% y in blocks for y
  if (!inherits(x = y, 'DelayedMatrix')) {
    stop('y should a DelayedMatrix')
  }
  if (nrow(x) != nrow(y)) {
    stop('row of x and y should be the same')
  }
  sparse <- DelayedArray::is_sparse(x = y)
  suppressMessages(setAutoBlockSize(size = block.size))
  cells.grid <- DelayedArray::colAutoGrid(x = y)
  norm.list <- list()
  for (i in seq_len(length.out = length(x = cells.grid))) {
    vp <- cells.grid[[i]]
    block <- DelayedArray::read_block(x = y, viewport = vp, as.sparse = sparse)
    if (sparse) {
      block <- as(object = block, Class = 'dgCMatrix')
    } else {
      block <- as(object = block, Class = 'Matrix')
    }
    norm.list[[i]] <- colSums(x = as.matrix(t(x) %*% block) ^ 2)
  }
  norm.vector <- unlist(norm.list)
  return(norm.vector)

}

# row mean from delayed array
#
RowMeanDelayedAssay <- function(x, block.size = 1e8) {
  if (!inherits(x = x, 'DelayedMatrix')) {
    stop('input x should a DelayedMatrix')
  }
  sparse <- DelayedArray::is_sparse(x = x)
  if (sparse ) {
    row.sum.function <- RowSumSparse
  } else {
    row.sum.function <- rowSums2
  }
  suppressMessages(setAutoBlockSize(size = block.size))
  cells.grid <- DelayedArray::colAutoGrid(x = x)
  sum.list <- list()
  for (i in seq_len(length.out = length(x = cells.grid))) {
    vp <- cells.grid[[i]]
    block <- DelayedArray::read_block(x = x, viewport = vp, as.sparse = sparse)
    if (sparse) {
      block <- as(object = block, Class = 'dgCMatrix')
    } else {
      block <- as(object = block, Class = 'Matrix')
    }
    sum.list[[i]] <- row.sum.function(mat = block)
  }
  mean.mat <- Reduce('+', sum.list)
  mean.mat <- mean.mat/ncol(x)
  return(mean.mat)
}

# row variance from delayed array
#
RowVarDelayedAssay <- function(x, block.size = 1e8) {
  if (!inherits(x = x, 'DelayedMatrix')) {
    stop('input x should a DelayedMatrix')
  }
  sparse <- DelayedArray::is_sparse(x = x)
  if (sparse ) {
    row.sum.function <- RowSumSparse
  } else {
    row.sum.function <- rowSums2
  }

  suppressMessages(setAutoBlockSize(size = block.size))
  cells.grid <- DelayedArray::colAutoGrid(x = x)
  sum2.list <- list()
  sum.list <- list()

  for (i in seq_len(length.out = length(x = cells.grid))) {
    vp <- cells.grid[[i]]
    block <- DelayedArray::read_block(x = x, viewport = vp, as.sparse = sparse)
    if (sparse) {
      block <- as(object = block, Class = 'dgCMatrix')
    } else {
      block <- as(object = block, Class = 'Matrix')
    }
    sum2.list[[i]] <- row.sum.function(mat = block**2)
    sum.list[[i]] <- row.sum.function(mat = block)
  }
  sum.mat <- Reduce('+', sum.list)
  sum2.mat <- Reduce('+', sum2.list)
  var.mat <- sum2.mat/ncol(x) - (sum.mat/ncol(x))**2
  var.mat <- var.mat * ncol(counts) / (ncol(counts) - 1)
  return(var.mat)
}



# nonzero element version of sweep
#
SweepNonzero <- function(
    x,
    MARGIN,
    STATS,
    FUN = "/"
) {
  if (!inherits(x = x, what = 'dgCMatrix')) {
    stop('input should be dgCMatrix. eg: x <- as(x, "CsparseMatrix")')
  }
  if (dim(x = x)[MARGIN] != length(STATS)){
    warning("Length of STATS is not equal to dim(x)[MARGIN]")
  }
  fun <- match.fun(FUN)
  if (MARGIN == 1) {
    idx <- x@i + 1
    x@x <- fun(x@x, STATS[idx])
  } else if (MARGIN == 2) {
    x <- as(x, "RsparseMatrix")
    idx <- x@j + 1
    x@x <- fun(x@x, STATS[idx])
    x <- as(x, "CsparseMatrix")
  }
  return(x)
}


#' Create one hot matrix for a given label
#' @importFrom Matrix colSums sparse.model.matrix
#' @importFrom stats as.formula
#' @export

CreateCategoryMatrix <- function(
  labels,
  method = c('aggregate', 'average'),
  cells.name = NULL
  ) {
  method <- match.arg(arg = method)
  if (is.null(dim(labels))) {
    if (length(x = unique(x = labels)) == 1) {
      data <- matrix(nrow = length(x = labels), ncol = 0)
    } else {
      data <- cbind(labels = labels)
    }
  } else {
    data <- labels 
  }
  cells.name <- cells.name %||% rownames(data)
  if (!is.null(cells.name) & length(cells.name) != nrow(data)) {
    stop('length of cells name should be equal to the length of input labels')
  }
  if (ncol(x = data) == 0) {
    message("All grouping variables have 1 value only. Computing across all cells.")
    category.matrix <- matrix(
      data = 1,
      nrow = nrow(x = data),
      dimnames = list(cells.name, 'all')
    )
    if (method == 'average') {
      category.matrix <- category.matrix / sum(category.matrix)
    }
    return(category.matrix)
  }
  group.by <- colnames(x = data)
  category.matrix <- sparse.model.matrix(object = as.formula(
    object = paste0(
      '~0+',
      paste0(
        "data[,",
        1:length(x = group.by),
        "]",
        collapse = ":"
      )
    )
  ))
  colsums <- colSums(x = category.matrix)
  category.matrix <- category.matrix[, colsums > 0]
  colsums <- colsums[colsums > 0]
  
  if (method =='average') {
    category.matrix <- SweepNonzero(
      x = category.matrix,
      MARGIN = 2,
      STATS = colsums,
      FUN = "/")
  }
  colnames(x = category.matrix) <- gsub(pattern = '_',
                                        replacement = '-',
                                        x = colnames(x = category.matrix)
                                        )
  colnames(x = category.matrix) <- sapply(
    X = colnames(x = category.matrix),
    FUN = function(name) {
      name <- gsub(pattern = "data\\[, [1-9]*\\]", replacement = "", x = name)
      return(paste0(rev(x = unlist(x = strsplit(x = name, split = ":"))), collapse = "_"))
    })
  rownames(category.matrix) <- cells.name
  return(category.matrix)
}

#' Construct an assay for spatial niche analysis
#'
#' This function will construct a new assay where each feature is a
#' cell label The values represents the sum of a particular cell label
#' neighboring a given cell.
#'
#' @param object A Seurat object
#' @param fov FOV object to gather cell positions from
#' @param group.by Cell classifications to count in spatial neighborhood
#' @param assay Name for spatial neighborhoods assay
#' @param neighbors.k Number of neighbors to consider for each cell
#' @param niches.k Number of clusters to return based on the niche assay
#' 
#' @importFrom stats kmeans
#' @return Seurat object containing a new assay
#' @concept clustering
#' @export
#'
BuildNicheAssay <- function(
  object,
  fov,
  group.by,
  assay = "niche",
  neighbors.k = 20,
  niches.k = 4
) {
  # find neighbors based on tissue position
  coords <- GetTissueCoordinates(object[[fov]], which = "centroids")
  cells <- coords$cell
  rownames(coords) <- cells
  coords <- as.matrix(coords[ , c("x", "y")])
  neighbors <- FindNeighbors(coords, k.param = neighbors.k)
  neighbors$nn <- neighbors$nn[Cells(object), Cells(object)]
  
  # build cell x cell type matrix
  ct.mtx <- matrix(
    data = 0,
    nrow = length(cells),
    ncol = length(unlist(unique(object[[group.by]])))
  )
  rownames(ct.mtx) <- cells
  colnames(ct.mtx) <- unique(unlist(object[[group.by]]))
  cts <- object[[group.by]]
  for (i in 1:length(cells)) {
    ct <- as.character(cts[cells[[i]], ])
    ct.mtx[cells[[i]], ct] <- 1
  }
  
  # create niche assay
  sum.mtx <- as.matrix(neighbors$nn %*% ct.mtx)
  niche.assay <- CreateAssayObject(counts = t(sum.mtx))
  object[[assay]] <- niche.assay
  DefaultAssay(object) <- assay
  
  # cluster niches assay
  object <- ScaleData(object)
  results <- kmeans(
    x = t(object[[assay]]@scale.data),
    centers = niches.k,
    nstart = 30
  )
  object$niches <- results[["cluster"]]
  
  return(object)  
}
