#' @include generics.R
#'
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

globalVariables(
  names = c('myAUC', 'p_val', 'avg_logFC'),
  package = 'Seurat',
  add = TRUE
)
#' Gene expression markers for all identity classes
#'
#' Finds markers (differentially expressed genes) for each of the identity classes in a dataset
#'
#' @inheritParams FindMarkers
#' @param node A node to find markers for and all its children; requires
#' \code{\link{BuildClusterTree}} to have been run previously; replaces \code{FindAllMarkersNode}
#' @param return.thresh Only return markers that have a p-value < return.thresh, or a power > return.thresh (if the test is ROC)
#'
#' @return Matrix containing a ranked list of putative markers, and associated
#' statistics (p-values, ROC score, etc.)
#'
#' @importFrom stats setNames
#'
#' @export
#'
#' @aliases FindAllMarkersNode
#' @concept differential_expression
#'
#' @examples
#' data("pbmc_small")
#' # Find markers for all clusters
#' all.markers <- FindAllMarkers(object = pbmc_small)
#' head(x = all.markers)
#' \dontrun{
#' # Pass a value to node as a replacement for FindAllMarkersNode
#' pbmc_small <- BuildClusterTree(object = pbmc_small)
#' all.markers <- FindAllMarkers(object = pbmc_small, node = 4)
#' head(x = all.markers)
#' }
#'
FindAllMarkers <- function(
  object,
  assay = NULL,
  features = NULL,
  group.by = NULL,
  logfc.threshold = 0.1,
  test.use = 'wilcox',
  slot = 'data',
  min.pct = 0.01,
  min.diff.pct = -Inf,
  node = NULL,
  verbose = TRUE,
  only.pos = FALSE,
  max.cells.per.ident = Inf,
  random.seed = 1,
  latent.vars = NULL,
  min.cells.feature = 3,
  min.cells.group = 3,
  mean.fxn = NULL,
  fc.name = NULL,
  base = 2,
  return.thresh = 1e-2,
  densify = FALSE,
  ...
) {
  MapVals <- function(vec, from, to) {
    vec2 <- setNames(object = to, nm = from)[as.character(x = vec)]
    vec2[is.na(x = vec2)] <- vec[is.na(x = vec2)]
    return(unname(obj = vec2))
  }
  GetNormMethod <- function(object, assay) {
    norm.command <- paste0("NormalizeData.", assay)
    if (norm.command %in% Command(object = object)) {
      return(Command(
        object = object,
        command = norm.command,
        value = "normalization.method"
      ))
    }
    transfer.command <- intersect(
      x = c("FindIntegrationAnchors", "FindTransferAnchors"),
      y = Command(object = object)
    )
    if (length(x = transfer.command)) {
      return(Command(
        object = object,
        command = transfer.command[1],
        value = "normalization.method"
      ))
    }
    return(NULL)
  }
  if ((test.use == "roc") && (return.thresh == 1e-2)) {
    return.thresh <- 0.7
  }
  if (is.null(x = node)) {
    if (!is.null(x = group.by) && !identical(x = group.by, y = "ident")) {
      if (length(x = group.by) == 1 && ! group.by %in% colnames(x = object@meta.data)) {
        stop("'", group.by, "' not found in object metadata")
      }
      Idents(object = object) <- group.by
    }
    idents.all <- sort(x = unique(x = Idents(object = object)))
  } else {
    if (!requireNamespace("ape", quietly = TRUE)) {
      stop(cluster.ape, call. = FALSE)
    }
    if (!is.null(group.by)) {
      warning(
        paste0(
          "The `group.by` parameter for `FindAllMarkers` ",
          "is ignored when `node` is set."
        )
      )
    }
    tree <- Tool(object = object, slot = 'BuildClusterTree')
    if (is.null(x = tree)) {
      stop("Please run 'BuildClusterTree' before finding markers on nodes")
    }
    descendants <- DFT(tree = tree, node = node, include.children = TRUE)
    all.children <- sort(x = tree$edge[, 2][!tree$edge[, 2] %in% tree$edge[, 1]])
    descendants <- MapVals(
      vec = descendants,
      from = all.children,
      to = tree$tip.label
    )
    drop.children <- setdiff(x = tree$tip.label, y = descendants)
    keep.children <- setdiff(x = tree$tip.label, y = drop.children)
    orig.nodes <- c(
      node,
      as.numeric(x = setdiff(x = descendants, y = keep.children))
    )
    tree <- ape::drop.tip(phy = tree, tip = drop.children)
    new.nodes <- unique(x = tree$edge[, 1, drop = TRUE])
    idents.all <- (tree$Nnode + 2):max(tree$edge)
  }
  genes.de <- vector(mode = "list", length = length(x = idents.all))
  messages <- vector(mode = "list", length = length(x = idents.all))
  assay <- assay %||% DefaultAssay(object = object)
  data.use <- NULL
  cellnames.use <- NULL
  cells.by.ident <- NULL
  latent.vars.data <- NULL
  norm.method <- NULL
  if (is.null(x = node)) {
    data.use <- object[[assay]]
    cellnames.use <- colnames(x = data.use)
    idents.use <- Idents(object = object)[cellnames.use]
    ident.values <- as.character(x = idents.use)
    cells.by.ident <- split(x = cellnames.use, f = idents.use)
    latent.vars.data <- if (is.null(x = latent.vars)) {
      NULL
    } else {
      FetchData(
        object = object,
        vars = latent.vars,
        cells = cellnames.use
      )
    }
    norm.method <- GetNormMethod(object = object, assay = assay)
  }
  # For the common Wilcoxon one-vs-rest case, presto can test all groups in
  # one matrix pass; unsupported cases return NULL and use the original loop.
  genes.de.presto <- FindAllMarkersWilcoxPresto(
    object = data.use,
    cells.by.ident = cells.by.ident,
    cellnames.use = cellnames.use,
    idents.use = idents.use,
    idents.all = idents.all,
    features = features,
    slot = slot,
    test.use = test.use,
    min.pct = min.pct,
    min.diff.pct = min.diff.pct,
    min.cells.group = min.cells.group,
    logfc.threshold = logfc.threshold,
    latent.vars = latent.vars,
    max.cells.per.ident = max.cells.per.ident,
    only.pos = only.pos,
    densify = densify,
    mean.fxn = mean.fxn,
    fc.name = fc.name,
    base = base,
    norm.method = norm.method,
    ...
  )
  if (!is.null(x = genes.de.presto)) {
    return(FindAllMarkersFinalize(
      genes.de = genes.de.presto,
      idents.all = idents.all,
      test.use = test.use,
      return.thresh = return.thresh,
      only.pos = only.pos
    ))
  }
  old_opt <- options(Seurat.warn.findmarkers.bpcells.colmajor = TRUE)
  on.exit(options(old_opt), add = TRUE)
  withCallingHandlers({
    for (i in seq_along(along.with = idents.all)) {
      if (verbose) {
        message("Calculating cluster ", idents.all[i])
      }
      genes.de[[i]] <- tryCatch(
       expr = {
          if (is.null(x = node)) {
            cells.1 <- cells.by.ident[[as.character(x = idents.all[i])]]
            cells.2 <- cellnames.use[ident.values != as.character(x = idents.all[i])]
            latent.vars.i <- if (is.null(x = latent.vars.data)) {
              NULL
            } else {
              latent.vars.data[c(cells.1, cells.2), , drop = FALSE]
            }
            FindMarkers(
              object = data.use,
              cells.1 = cells.1,
              cells.2 = cells.2,
              features = features,
              logfc.threshold = logfc.threshold,
              test.use = test.use,
              slot = slot,
              min.pct = min.pct,
              min.diff.pct = min.diff.pct,
              verbose = verbose,
              only.pos = only.pos,
              max.cells.per.ident = max.cells.per.ident,
              random.seed = random.seed,
              latent.vars = latent.vars.i,
              min.cells.feature = min.cells.feature,
              min.cells.group = min.cells.group,
              mean.fxn = mean.fxn,
              fc.name = fc.name,
              base = base,
              densify = densify,
              norm.method = norm.method,
              ...
            )
          } else {
            FindMarkers(
              object = object,
              assay = assay,
              ident.1 = tree,
              ident.2 = idents.all[i],
              features = features,
              logfc.threshold = logfc.threshold,
              test.use = test.use,
              slot = slot,
              min.pct = min.pct,
              min.diff.pct = min.diff.pct,
              verbose = verbose,
              only.pos = only.pos,
              max.cells.per.ident = max.cells.per.ident,
              random.seed = random.seed,
              latent.vars = latent.vars,
              min.cells.feature = min.cells.feature,
              min.cells.group = min.cells.group,
              mean.fxn = mean.fxn,
              fc.name = fc.name,
              base = base,
              densify = densify,
              ...
            )
          }
        },
        error = function(cond) {
          return(cond$message)
        }
      )
      if (is.character(x = genes.de[[i]])) {
        messages[[i]] <- genes.de[[i]]
        genes.de[[i]] <- NULL
      }
    }
  }, warning = function(w) {
    if (grepl("Column-major order detected", conditionMessage(w))) {
      msg <- conditionMessage(w)
      msg <- sub("\\.\\nThis message will be shown once per session\\.$", " to avoid repeated transpositions.", msg)
      warning(msg, immediate. = TRUE, call. = FALSE)
      invokeRestart("muffleWarning")
    }
  })
  return(FindAllMarkersFinalize(
    genes.de = genes.de,
    idents.all = idents.all,
    test.use = test.use,
    return.thresh = return.thresh,
    only.pos = only.pos,
    messages = messages,
    node = node,
    new.nodes = if (exists(x = "new.nodes")) new.nodes else NULL,
    orig.nodes = if (exists(x = "orig.nodes")) orig.nodes else NULL
  ))
}

#' Finalize FindAllMarkers Results
#'
#' Combine per-ident \code{FindMarkers}-style data.frames into the public
#' \code{FindAllMarkers} table. Both the legacy loop and fast Presto path call
#' this so row ordering, \code{return.thresh}, cluster/gene columns, and
#' warnings stay aligned.
#'
#' @param genes.de List of per-ident differential expression result data.frames.
#' @param idents.all Identity classes tested by \code{FindAllMarkers}.
#' @param test.use Differential expression test used to generate \code{genes.de}.
#' @param return.thresh Threshold for retaining markers in the final table.
#' @param only.pos Only retain markers with positive fold-change values.
#' @param messages Optional list of per-ident messages for tests that were not
#' run.
#' @param node Optional tree node tested by \code{FindAllMarkers}.
#' @param new.nodes,orig.nodes Optional node label mappings used when
#' \code{node} is set.
#'
#' @return A data.frame in \code{FindAllMarkers} output format.
#' @noRd
FindAllMarkersFinalize <- function(
  genes.de,
  idents.all,
  test.use,
  return.thresh,
  only.pos,
  messages = NULL,
  node = NULL,
  new.nodes = NULL,
  orig.nodes = NULL
) {
  gde.all <- data.frame()
  for (i in seq_along(along.with = idents.all)) {
    gde <- genes.de[[i]]
    if (is.null(x = gde) || nrow(x = gde) == 0) {
      next
    }
    if (test.use == "roc") {
      gde <- subset(
        x = gde,
        subset = (myAUC > return.thresh | myAUC < (1 - return.thresh))
      )
    } else if (is.null(x = node) || test.use %in% c('bimod', 't')) {
      gde <- gde[order(gde$p_val, -abs(gde$pct.1-gde$pct.2)), , drop = FALSE]
      gde <- subset(x = gde, subset = p_val < return.thresh)
    }
    if (nrow(x = gde) == 0) {
      next
    }
    gde$cluster <- idents.all[i]
    gde$gene <- rownames(x = gde)
    gde.all <- rbind(gde.all, gde)
  }
  if ((only.pos) && nrow(x = gde.all) > 0) {
    return(subset(x = gde.all, subset = gde.all[, 2] > 0))
  }
  rownames(x = gde.all) <- make.unique(names = as.character(x = gde.all$gene))
  if (nrow(x = gde.all) == 0) {
    warning("No DE genes identified", call. = FALSE, immediate. = TRUE)
  }
  if (any(vapply(
    X = messages,
    FUN = function(x) !is.null(x),
    FUN.VALUE = logical(length = 1L)
  ))) {
    warning("The following tests were not performed: ", call. = FALSE, immediate. = TRUE)
    for (i in seq_along(along.with = messages)) {
      if (!is.null(x = messages[[i]])) {
        warning("When testing ", idents.all[i], " versus all:\n\t", messages[[i]], call. = FALSE, immediate. = TRUE)
      }
    }
  }
  if (!is.null(x = node)) {
    MapVals <- function(vec, from, to) {
      vec2 <- setNames(object = to, nm = from)[as.character(x = vec)]
      vec2[is.na(x = vec2)] <- vec[is.na(x = vec2)]
      return(unname(obj = vec2))
    }
    gde.all$cluster <- MapVals(
      vec = gde.all$cluster,
      from = new.nodes,
      to = orig.nodes
    )
  }
  return(gde.all)
}

#' Get a full per-ident DE list for the fast Presto Wilcoxon path
#'
#' @inheritParams FindAllMarkers
#' @param cells.by.ident Named list of cell names split by identity.
#' @param cellnames.use Cells included in the one-vs-rest tests.
#' @param idents.use Identity vector named by cell.
#' @param idents.all Identity classes to test.
#' @param norm.method Normalization method recorded for the assay, when
#' available.
#'
#' @return A list of per-ident DE data.frames (\code{NULL} when fallback is required)
#' @noRd
FindAllMarkersWilcoxPresto <- function(
  object,
  cells.by.ident,
  cellnames.use,
  idents.use,
  idents.all,
  features,
  slot,
  test.use,
  min.pct,
  min.diff.pct,
  min.cells.group,
  logfc.threshold,
  latent.vars,
  max.cells.per.ident,
  only.pos,
  densify,
  mean.fxn,
  fc.name,
  base,
  norm.method,
  ...
) {
  if (is.null(x = object) || test.use != "wilcox" || max.cells.per.ident < Inf) {
    return(NULL)
  }
  ident.sizes <- vapply(X = cells.by.ident, FUN = length, FUN.VALUE = integer(length = 1L))
  layer.data <- tryCatch(
    expr = LayerData(object = object, layer = slot),
    error = function(...) NULL
  )
  data.is.iterable <- is.null(x = layer.data) || inherits(x = layer.data, what = "IterableMatrix")
  # presto::wilcoxauc computes all one-vs-rest Wilcoxon tests in one pass.
  # Keep BPCells and sampled comparisons on the old path because they have
  # backend-specific behavior or per-cluster sampling that is not equivalent.
  if (
    length(x = ident.sizes) <= 1 ||
    any(ident.sizes < min.cells.group) ||
    any((length(x = cellnames.use) - ident.sizes) < min.cells.group) ||
    data.is.iterable ||
    !requireNamespace("presto", quietly = TRUE)
  ) {
    return(NULL)
  }
  if (!is.null(x = latent.vars)) {
    warning(
      "'latent.vars' is only used for the following tests: ",
      paste(DEmethods_latent(), collapse=", "),
      call. = FALSE,
      immediate. = TRUE
    )
  }

  tryCatch(
    expr = FindAllMarkersWilcoxPrestoDE(
      object = object,
      cells.by.ident = cells.by.ident,
      cellnames.use = cellnames.use,
      idents.use = idents.use,
      idents.all = idents.all,
      features = features,
      slot = slot,
      min.pct = min.pct,
      min.diff.pct = min.diff.pct,
      logfc.threshold = logfc.threshold,
      only.pos = only.pos,
      densify = densify,
      mean.fxn = mean.fxn,
      fc.name = fc.name,
      base = base,
      norm.method = norm.method,
      layer.data = layer.data,
      ...
    ),
    error = function(...) NULL
  )
}

#' Filter features By FoldChange statistics
#'
#' Apply the same pct and fold-change feature filters used by
#' \code{FindMarkers.default} to one per-ident fold-change table
#'
#' @param fc.results Fold-change result data.frame for one identity.
#' @param slot Assay data slot used for differential expression testing.
#' @param min.pct Minimum fraction of cells expressing a feature in either
#' group.
#' @param min.diff.pct Minimum difference in detection fraction between groups.
#' @param logfc.threshold Minimum fold-change threshold.
#' @param only.pos Only retain features with positive fold-change values.
#'
#' @return A character vector of features that pass filters
#' @noRd
FindAllMarkersFilterFeatures <- function(
  fc.results,
  slot,
  min.pct,
  min.diff.pct,
  logfc.threshold,
  only.pos
) {
  # Match FindMarkers.default prefiltering: first require expression in either
  # group, then require a minimum pct difference, then apply the FC threshold
  alpha.min <- pmax(fc.results$pct.1, fc.results$pct.2)
  names(x = alpha.min) <- rownames(x = fc.results)
  features.use <- names(x = which(x = alpha.min >= min.pct))
  if (!length(x = features.use)) {
    return(features.use)
  }
  alpha.diff <- alpha.min - pmin(fc.results$pct.1, fc.results$pct.2)
  features.use <- names(
    x = which(x = alpha.min >= min.pct & alpha.diff >= min.diff.pct)
  )
  if (!length(x = features.use) || slot == "scale.data") {
    return(features.use)
  }
  total.diff <- fc.results[, 1]
  names(x = total.diff) <- rownames(x = fc.results)
  features.diff <- if (only.pos) {
    names(x = which(x = total.diff >= logfc.threshold))
  } else {
    names(x = which(x = abs(x = total.diff) >= logfc.threshold))
  }
  intersect(x = features.use, y = features.diff)
}

#' Bulk FoldChange For FindAllMarkers
#'
#' Compute fold-change and detection percentages for all identities in one pass when 
#' the requested mean function is available as a default
#'
#' @inheritParams FindAllMarkers
#' @param cellnames.use Cells included in the one-vs-rest tests.
#' @param idents.use Identity vector named by cell.
#' @param idents.all Identity classes to test.
#' @param fc.slot Assay data slot used for fold-change calculations.
#' @param pseudocount.use Pseudocount added to average expression values.
#' @param norm.method Normalization method recorded for the assay, when
#' available.
#' @param data Optional matrix to reuse for fold-change calculations.
#'
#' @return A named list of fold-change data.frames (\code{NULL} when fallback is required)
#' @noRd
FoldChangeFindAllMarkers <- function(
  object,
  features,
  cellnames.use,
  idents.use,
  idents.all,
  fc.slot,
  pseudocount.use,
  fc.name,
  mean.fxn,
  base,
  norm.method,
  slot,
  min.pct,
  min.diff.pct,
  logfc.threshold,
  only.pos,
  data = NULL
) {
  # Bulk FC is exact only for Seurat's built-in mean functions, not custom functions
  if (!is.null(x = mean.fxn) || !fc.slot %in% c("counts", "data", "scale.data")) {
    return(NULL)
  }
  data <- data %||% tryCatch(
    expr = GetAssayData(object = object, layer = fc.slot),
    error = function(...) NULL
  )
  if (is.null(x = data)) {
    return(NULL)
  }
  features <- intersect(x = features, y = rownames(x = data))
  if (!length(x = features)) {
    return(vector(mode = "list", length = length(idents.all)))
  }
  pseudocount.use <- pseudocount.use %||% 1
  groups <- factor(
    x = as.character(x = idents.use[cellnames.use]),
    levels = as.character(x = idents.all)
  )
  group.n <- tabulate(bin = as.integer(x = groups), nbins = length(x = idents.all))
  rest.n <- length(x = cellnames.use) - group.n
  data <- if (
    identical(x = features, y = rownames(x = data)) &&
    identical(x = cellnames.use, y = colnames(x = data))
  ) {
    data
  } else {
    data[features, cellnames.use, drop = FALSE]
  }
  feature.names <- rownames(x = data)
  sparse.stats <- FindAllMarkersSparseStats(
    data = data,
    groups = groups,
    n.groups = length(x = idents.all),
    fc.slot = fc.slot,
    norm.method = norm.method
  )
  if (!is.null(x = sparse.stats)) {
    dimnames(x = sparse.stats$detected) <- list(feature.names, as.character(x = idents.all))
    dimnames(x = sparse.stats$group_sum) <- list(feature.names, as.character(x = idents.all))
    dimnames(x = sparse.stats$rest_sum) <- list(feature.names, as.character(x = idents.all))
    names(x = sparse.stats$total_detected) <- feature.names
  }
  # Detection percentages use the same >0 threshold and 3-digit rounding as
  # FoldChange.default, but derive pct.2 from total-minus-group counts
  if (is.null(x = sparse.stats)) {
    # Cell x group membership matrix, multiplying feature x cell data by this
    # matrix gives feature x group sums/counts for all identities at once
    group.mat <- Matrix::sparseMatrix(
      i = seq_along(along.with = groups),
      j = as.integer(x = groups),
      x = 1,
      dims = c(length(x = groups), length(x = idents.all))
    )
    detected.1 <- as.matrix(x = (data > 0) %*% group.mat)
    detected.total <- Matrix::rowSums(x = data > 0)
  } else {
    detected.1 <- sparse.stats$detected
    detected.total <- sparse.stats$total_detected
  }
  pct.1 <- round(x = DivideColumnsByDenomList(x = detected.1, denominator = group.n), digits = 3)
  pct.2 <- round(x = DivideColumnsByDenomList(x = detected.total - detected.1, denominator = rest.n), digits = 3)
  base.text <- ifelse(test = base == exp(1), yes = "", no = base)
  fc.name <- fc.name %||% ifelse(
    test = fc.slot == "scale.data",
    yes = "avg_diff",
    no = paste0("avg_log", base.text, "FC")
  )
  features.by.pct <- vector(mode = "list", length = length(idents.all))
  names(x = features.by.pct) <- as.character(x = idents.all)
  for (i in seq_along(along.with = idents.all)) {
    alpha.min <- pmax(pct.1[, i], pct.2[, i])
    features.use <- feature.names[which(x = alpha.min >= min.pct)]
    if (length(x = features.use)) {
      alpha.diff <- alpha.min - pmin(pct.1[, i], pct.2[, i])
      features.use <- feature.names[which(x = alpha.min >= min.pct & alpha.diff >= min.diff.pct)]
    }
    features.by.pct[[i]] <- features.use
  }
  features.fc <- unique(x = unlist(x = features.by.pct, use.names = FALSE))
  fc.by.ident <- vector(mode = "list", length = length(idents.all))
  names(x = fc.by.ident) <- as.character(x = idents.all)
  if (!length(x = features.fc)) {
    for (i in seq_along(along.with = idents.all)) {
      fc.by.ident[[i]] <- data.frame(
        fc = numeric(length = 0L),
        pct.1 = numeric(length = 0L),
        pct.2 = numeric(length = 0L)
      )
      colnames(x = fc.by.ident[[i]])[1] <- fc.name
    }
    attr(x = fc.by.ident, which = "features.filtered") <- TRUE
    return(fc.by.ident)
  }
  use.full.fc <- length(x = features.fc) > (0.5 * nrow(x = data))
  data.fc <- if (use.full.fc) {
    data
  } else {
    data[features.fc, , drop = FALSE]
  }
  if (is.null(x = sparse.stats)) {
    data.mean <- switch(
      EXPR = fc.slot,
      "scale.data" = data.fc,
      "data" = if (is.null(x = norm.method) || norm.method == "LogNormalize") {
        expm1(x = data.fc)
      } else {
        data.fc
      },
      data.fc
    )
    group.sum <- as.matrix(x = data.mean %*% group.mat)
    rest.sum <- matrix(
      data = 0,
      nrow = nrow(x = group.sum),
      ncol = ncol(x = group.sum),
      dimnames = dimnames(x = group.sum)
    )
    for (i in seq_len(length.out = ncol(x = group.sum))) {
      rest.sum[, i] <- if (ncol(x = group.sum) == 2L) {
        group.sum[, 3L - i]
      } else {
        Matrix::rowSums(x = group.sum[, -i, drop = FALSE])
      }
    }
  } else {
    group.sum <- sparse.stats$group_sum[features.fc, , drop = FALSE]
    rest.sum <- sparse.stats$rest_sum[features.fc, , drop = FALSE]
  }
  if (fc.slot == "scale.data") {
    mean.1 <- DivideColumnsByDenomList(x = group.sum, denominator = group.n)
    mean.2 <- DivideColumnsByDenomList(x = rest.sum, denominator = rest.n)
  } else {
    mean.1 <- log(
      x = DivideColumnsByDenomList(x = group.sum + pseudocount.use, denominator = group.n),
      base = base
    )
    mean.2 <- log(
      x = DivideColumnsByDenomList(
        x = rest.sum + pseudocount.use,
        denominator = rest.n
      ),
      base = base
    )
  }
  fc <- mean.1 - mean.2
  for (i in seq_along(along.with = idents.all)) {
    features.use <- features.by.pct[[i]]
    fc.use <- numeric(length = 0L)
    if (length(x = features.use)) {
      fc.use <- fc[features.use, i]
      names(x = fc.use) <- features.use
    }
    if (length(x = features.use) && slot != "scale.data") {
      total.diff <- fc.use
      features.diff <- if (only.pos) {
        names(x = which(x = total.diff >= logfc.threshold))
      } else {
        names(x = which(x = abs(x = total.diff) >= logfc.threshold))
      }
      features.use <- intersect(x = features.use, y = features.diff)
    }
    fc.by.ident[[i]] <- data.frame(
      fc = fc.use[features.use],
      pct.1 = pct.1[features.use, i],
      pct.2 = pct.2[features.use, i],
      row.names = features.use
    )
    colnames(x = fc.by.ident[[i]])[1] <- fc.name
  }
  attr(x = fc.by.ident, which = "features.filtered") <- TRUE
  return(fc.by.ident)
}

#' Sparse FoldChange Summary Statistics
#'
#' Compute detection counts and group/rest sums for supported sparse matrices
#' using the C++ FindAllMarkers fold-change helper.
#'
#' @param data Expression matrix.
#' @param groups Factor of identity assignments for the columns of \code{data}.
#' @param n.groups Number of identity groups.
#' @param fc.slot Assay data slot used for fold-change calculations.
#' @param norm.method Normalization method recorded for the assay, when
#' available.
#'
#' @return A list of sparse summary matrices, or \code{NULL} when unsupported.
#' @noRd
FindAllMarkersSparseStats <- function(
  data,
  groups,
  n.groups,
  fc.slot,
  norm.method
) {
  if (!inherits(x = data, what = "dgCMatrix") || !fc.slot %in% c("counts", "data")) {
    return(NULL)
  }
  FindAllMarkersSparseFoldChangeStats(
    x = data@x,
    i = data@i,
    p = data@p,
    rows = nrow(x = data),
    cols = ncol(x = data),
    groups = as.integer(x = groups),
    n_groups = n.groups,
    log_normalize = fc.slot == "data" && (is.null(x = norm.method) || norm.method == "LogNormalize"),
    nthreads = getThreads()
  )
}

#' FindAllMarkers Wilcoxon w/ presto
#'
#' Compute all one-vs-rest Wilcoxon p-values with one \code{presto::wilcoxauc} call,
#' and attach \code{Seurat::FoldChange} results for each identity.
#'
#' @inheritParams FindAllMarkers
#' @param cells.by.ident Named list of cell names split by identity.
#' @param cellnames.use Cells included in the one-vs-rest tests.
#' @param idents.use Identity vector named by cell.
#' @param idents.all Identity classes to test.
#' @param norm.method Normalization method recorded for the assay, when
#' available.
#' @param layer.data Optional preloaded assay data for \code{slot}.
#'
#' @return A list of per-ident DE data.frames
#' @noRd
FindAllMarkersWilcoxPrestoDE <- function(
  object,
  cells.by.ident,
  cellnames.use,
  idents.use,
  idents.all,
  features,
  slot,
  min.pct,
  min.diff.pct,
  logfc.threshold,
  only.pos,
  densify,
  mean.fxn,
  fc.name,
  base,
  norm.method,
  layer.data = NULL,
  ...
) {
  dots <- list(...)
  fc.slot <- dots$fc.slot %||% "data"
  dots$fc.slot <- NULL
  if (inherits(x = object, what = "SCTAssay")) {
    ValidateSCTFindMarkers(
      object = object,
      recorrect_umi = dots$recorrect_umi %||% TRUE
    )
    layer.data <- NULL
  }
  dots$recorrect_umi <- NULL
  data.use <- layer.data %||% LayerData(object = object, layer = slot)
  features <- features %||% rownames(x = data.use)
  features <- intersect(x = features, y = rownames(x = data.use))
  if (!length(x = features)) {
    return(vector(mode = "list", length = length(idents.all)))
  }
  fc.by.ident <- FoldChangeFindAllMarkers(
    object = object,
    features = features,
    cellnames.use = cellnames.use,
    idents.use = idents.use,
    idents.all = idents.all,
    fc.slot = fc.slot,
    pseudocount.use = dots$pseudocount.use,
    fc.name = fc.name,
    mean.fxn = mean.fxn,
    base = base,
    norm.method = norm.method,
    slot = slot,
    min.pct = min.pct,
    min.diff.pct = min.diff.pct,
    logfc.threshold = logfc.threshold,
    only.pos = only.pos,
    data = if (identical(x = fc.slot, y = slot)) data.use else NULL
  )
  if (is.null(x = fc.by.ident)) {
    # Fall back to the original FoldChange loop for custom or unsupported FC
    # semantics, while still using the single Presto call for p-values
    fc.by.ident <- vector(mode = "list", length = length(idents.all))
    names(x = fc.by.ident) <- as.character(x = idents.all)
    for (i in seq_along(along.with = idents.all)) {
      ident <- as.character(x = idents.all[i])
      cells.1 <- cells.by.ident[[ident]]
      cells.2 <- setdiff(x = cellnames.use, y = cells.1)
      fc.args <- list(
        object = object,
        cells.1 = cells.1,
        cells.2 = cells.2,
        features = features,
        slot = fc.slot,
        fc.name = fc.name,
        mean.fxn = mean.fxn,
        base = base,
        norm.method = norm.method
      )
      if (!is.null(x = dots$pseudocount.use)) {
        fc.args$pseudocount.use <- dots$pseudocount.use
      }
      fc.args <- c(fc.args, dots[setdiff(x = names(x = dots), y = "pseudocount.use")])
      fc.by.ident[[i]] <- do.call(what = FoldChange, args = fc.args)
    }
  }
  features.by.ident <- vector(mode = "list", length = length(idents.all))
  names(x = features.by.ident) <- as.character(x = idents.all)
  if (isTRUE(x = attr(x = fc.by.ident, which = "features.filtered"))) {
    for (i in seq_along(along.with = idents.all)) {
      features.by.ident[[i]] <- rownames(x = fc.by.ident[[i]])
    }
  } else {
    for (i in seq_along(along.with = idents.all)) {
      fc.results <- fc.by.ident[[i]]
      features.use <- FindAllMarkersFilterFeatures(
        fc.results = fc.results,
        slot = slot,
        min.pct = min.pct,
        min.diff.pct = min.diff.pct,
        logfc.threshold = logfc.threshold,
        only.pos = only.pos
      )
      fc.by.ident[[i]] <- fc.results[features.use, , drop = FALSE]
      features.by.ident[[i]] <- features.use
    }
  }
  features.test <- unique(x = unlist(x = features.by.ident, use.names = FALSE))
  if (!length(x = features.test)) {
    return(fc.by.ident)
  }
  use.full.test <- identical(x = cellnames.use, y = colnames(x = data.use))
  if (use.full.test) {
    data.test <- data.use
    features.test <- rownames(x = data.use)
  } else {
    data.test <- data.use[features.test, cellnames.use, drop = FALSE]
  }
  if (densify) {
    data.test <- as.matrix(x = data.test)
  }
  groups <- factor(
    x = as.character(x = idents.use[cellnames.use]),
    levels = as.character(x = idents.all)
  )
  presto.results <- presto::wilcoxauc(X = data.test, y = groups)
  pval.mat <- matrix(
    data = presto.results$pval,
    nrow = length(x = features.test),
    ncol = length(x = idents.all),
    dimnames = list(features.test, as.character(x = idents.all))
  )
  genes.de <- vector(mode = "list", length = length(idents.all))
  for (i in seq_along(along.with = idents.all)) {
    ident <- as.character(x = idents.all[i])
    fc.results <- fc.by.ident[[i]]
    if (is.null(x = fc.results) || !nrow(x = fc.results)) {
      genes.de[[i]] <- fc.results
      next
    }
    de.results <- data.frame(
      p_val = pval.mat[rownames(x = fc.results), ident],
      row.names = rownames(x = fc.results)
    )
    de.results <- cbind(de.results, fc.results[rownames(x = de.results), , drop = FALSE])
    de.results$p_val_adj <- pmin(1, de.results$p_val * nrow(x = data.use))
    genes.de[[i]] <- de.results
  }
  return(genes.de)
}

#' Finds markers that are conserved between the groups
#'
#' @inheritParams FindMarkers
#' @param ident.1 Identity class to define markers for
#' @param ident.2 A second identity class for comparison. If NULL (default) -
#' use all other cells for comparison.
#' @param grouping.var grouping variable
#' @param assay of assay to fetch data for (default is RNA)
#' @param meta.method method for combining p-values. Should be a function from
#' the metap package (NOTE: pass the function, not a string)
#' @param \dots parameters to pass to FindMarkers
#'
#' @return data.frame containing a ranked list of putative conserved markers, and
#' associated statistics (p-values within each group and a combined p-value
#' (such as Fishers combined p-value or others from the metap package),
#' percentage of cells expressing the marker, average differences). Name of group is appended to each
#' associated output column (e.g. CTRL_p_val). If only one group is tested in the grouping.var, max
#' and combined p-values are not returned.
#'
#' @export
#' @concept differential_expression
#'
#' @examples
#' \dontrun{
#' data("pbmc_small")
#' pbmc_small
#' # Create a simulated grouping variable
#' pbmc_small[['groups']] <- sample(x = c('g1', 'g2'), size = ncol(x = pbmc_small), replace = TRUE)
#' FindConservedMarkers(pbmc_small, ident.1 = 0, ident.2 = 1, grouping.var = "groups")
#' }
#'
FindConservedMarkers <- function(
  object,
  ident.1,
  ident.2 = NULL,
  grouping.var,
  assay = 'RNA',
  slot = 'data',
  min.cells.group = 3,
  meta.method = metap::minimump,
  verbose = TRUE,
  ...
) {
  metap.installed <- requireNamespace("metap", quietly = TRUE)
  if (!metap.installed[1]) {
    stop(
      "Please install the metap package to use FindConservedMarkers.",
      "\nThis can be accomplished with the following commands: ",
      "\n----------------------------------------",
      "\ninstall.packages('BiocManager')",
      "\nBiocManager::install('multtest')",
      "\ninstall.packages('metap')",
      "\n----------------------------------------",
      call. = FALSE
    )
  }
  if (!is.function(x = meta.method)) {
    stop("meta.method should be a function from the metap package. Please see https://cran.r-project.org/web/packages/metap/metap.pdf for a detailed description of the available functions.")
  }
  object.var <- FetchData(object = object, vars = grouping.var)
  object <- SetIdent(
    object = object,
    cells = colnames(x = object),
    value = paste(Idents(object = object), object.var[, 1], sep = "_")
  )
  levels.split <- names(x = sort(x = table(object.var[, 1])))
  num.groups <- length(levels.split)
  cells <- list()
  for (i in 1:num.groups) {
    cells[[i]] <- rownames(
      x = object.var[object.var[, 1] == levels.split[i], , drop = FALSE]
    )
  }
  marker.test <- list()
  # do marker tests
  ident.2.save <- ident.2
  for (i in 1:num.groups) {
    level.use <- levels.split[i]
    ident.use.1 <- paste(ident.1, level.use, sep = "_")
    ident.use.1.exists <- ident.use.1 %in% Idents(object = object)
    if (!all(ident.use.1.exists)) {
      bad.ids <- ident.1[!ident.use.1.exists]
      warning(
        "Identity: ",
        paste(bad.ids, collapse = ", "),
        " not present in group ",
        level.use,
        ". Skipping ",
        level.use,
        call. = FALSE,
        immediate. = TRUE
      )
      next
    }
    ident.2 <- ident.2.save
    cells.1 <- WhichCells(object = object, idents = ident.use.1)
    if (length(cells.1) < min.cells.group) {
      warning(
        level.use,
        " has fewer than ",
        min.cells.group,
        " cells in Identity: ",
        paste(ident.1, collapse = ", "),
        ". Skipping ",
        level.use,
        call. = FALSE,
        immediate. = TRUE
      )
      next
    }
    if (is.null(x = ident.2)) {
      cells.2 <- setdiff(x = cells[[i]], y = cells.1)
      ident.use.2 <- names(x = which(x = table(Idents(object = object)[cells.2]) > 0))
      ident.2 <- gsub(pattern = paste0("_", level.use), replacement = "", x = ident.use.2)
      if (length(x = ident.use.2) == 0) {
        stop(paste("Only one identity class present:", ident.1))
      }
    } else {
      ident.use.2 <- paste(ident.2, level.use, sep = "_")
    }
    if (verbose) {
      message(
        "Testing group ",
        level.use,
        ": (",
        paste(ident.1, collapse = ", "),
        ") vs (",
        paste(ident.2, collapse = ", "),
        ")"
      )
    }
    ident.use.2.exists <- ident.use.2 %in% Idents(object = object)
    if (!all(ident.use.2.exists)) {
      bad.ids <- ident.2[!ident.use.2.exists]
      warning(
        "Identity: ",
        paste(bad.ids, collapse = ", "),
        " not present in group ",
        level.use,
        ". Skipping ",
        level.use,
        call. = FALSE,
        immediate. = TRUE
      )
      next
    }
    marker.test[[i]] <- FindMarkers(
      object = object,
      assay = assay,
      slot = slot,
      ident.1 = ident.use.1,
      ident.2 = ident.use.2,
      verbose = verbose,
      min.cells.group =   min.cells.group,
      ...
    )
    names(x = marker.test)[i] <- levels.split[i]
  }
  marker.test <- Filter(f = Negate(f = is.null), x = marker.test)
  genes.conserved <- Reduce(
    f = intersect,
    x = lapply(
      X = marker.test,
      FUN = function(x) {
        return(rownames(x = x))
      }
    )
  )
  markers.conserved <- list()
  for (i in 1:length(x = marker.test)) {
    markers.conserved[[i]] <- marker.test[[i]][genes.conserved, ]
    colnames(x = markers.conserved[[i]]) <- paste(
      names(x = marker.test)[i],
      colnames(x = markers.conserved[[i]]),
      sep = "_"
    )
  }
  markers.combined <- Reduce(cbind, markers.conserved)
  pval.codes <- colnames(x = markers.combined)[grepl(pattern = "*_p_val$", x = colnames(x = markers.combined))]
  if (length(x = pval.codes) > 1) {
    markers.combined$max_pval <- apply(
      X = markers.combined[, pval.codes, drop = FALSE],
      MARGIN = 1,
      FUN = max
    )
    combined.pval <- data.frame(cp = apply(
      X = markers.combined[, pval.codes, drop = FALSE],
      MARGIN = 1,
      FUN = function(x) {
        clean_x <- x[!is.na(x) & !is.nan(x)]
        return(meta.method(clean_x)$p)
      }
    ))
    meta.method.name <- as.character(x = formals()$meta.method)
    if (length(x = meta.method.name) == 3) {
      meta.method.name <- meta.method.name[3]
    }
    colnames(x = combined.pval) <- paste0(meta.method.name, "_p_val")
    markers.combined <- cbind(markers.combined, combined.pval)
    markers.combined <- markers.combined[order(markers.combined[, paste0(meta.method.name, "_p_val")]), ]
  } else {
    warning("Only a single group was tested", call. = FALSE, immediate. = TRUE)
  }
  return(markers.combined)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Methods for Seurat-defined generics
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @param slot Slot to pull data from; note that if \code{test.use} is
#' "negbinom", "poisson", or "DESeq2", \code{slot} will be set to "counts"
#' @param cells.1 Vector of cell names belonging to group 1
#' @param cells.2 Vector of cell names belonging to group 2
#' @param features Genes to test. Default is to use all genes
#' @param slot Slot to pull data from; note that if \code{test.use} is
#' "negbinom", "poisson", or "DESeq2", \code{slot} will be set to "counts"
#' @param logfc.threshold Limit testing to genes which show, on average, at least
#' X-fold difference (log-scale) between the two groups of cells. Default is 0.1
#' Increasing logfc.threshold speeds up the function, but can miss weaker signals.
#' If the \code{slot} parameter is "scale.data" no filtering is performed.
#' @param test.use Denotes which test to use. Available options are:
#' \itemize{
#'  \item{"wilcox"} : Identifies differentially expressed genes between two
#'  groups of cells using a Wilcoxon Rank Sum test (default); will use a fast
#'  implementation by Presto if installed
#'  \item{"wilcox_limma"} : Identifies differentially expressed genes between two
#'  groups of cells using the limma implementation of the Wilcoxon Rank Sum test;
#'  set this option to reproduce results from Seurat v4
#'  \item{"bimod"} : Likelihood-ratio test for single cell gene expression,
#'  (McDavid et al., Bioinformatics, 2013)
#'  \item{"roc"} : Identifies 'markers' of gene expression using ROC analysis.
#'  For each gene, evaluates (using AUC) a classifier built on that gene alone,
#'  to classify between two groups of cells. An AUC value of 1 means that
#'  expression values for this gene alone can perfectly classify the two
#'  groupings (i.e. Each of the cells in cells.1 exhibit a higher level than
#'  each of the cells in cells.2). An AUC value of 0 also means there is perfect
#'  classification, but in the other direction. A value of 0.5 implies that
#'  the gene has no predictive power to classify the two groups. Returns a
#'  'predictive power' (abs(AUC-0.5) * 2) ranked matrix of putative differentially
#'  expressed genes.
#'  \item{"t"} : Identify differentially expressed genes between two groups of
#'  cells using the Student's t-test.
#'  \item{"negbinom"} : Identifies differentially expressed genes between two
#'   groups of cells using a negative binomial generalized linear model.
#'   Use only for UMI-based datasets
#'  \item{"poisson"} : Identifies differentially expressed genes between two
#'   groups of cells using a poisson generalized linear model.
#'   Use only for UMI-based datasets
#'  \item{"LR"} : Uses a logistic regression framework to determine differentially
#'  expressed genes. Constructs a logistic regression model predicting group
#'  membership based on each feature individually and compares this to a null
#'  model with a likelihood ratio test.
#'  \item{"MAST"} : Identifies differentially expressed genes between two groups
#'  of cells using a hurdle model tailored to scRNA-seq data. Utilizes the MAST
#'  package to run the DE testing.
#'  \item{"DESeq2"} : Identifies differentially expressed genes between two groups
#'  of cells based on a model using DESeq2 which uses a negative binomial
#'  distribution (Love et al, Genome Biology, 2014).This test does not support
#'  pre-filtering of genes based on average difference (or percent detection rate)
#'  between cell groups. However, genes may be pre-filtered based on their
#'  minimum detection rate (min.pct) across both cell groups. To use this method,
#'  please install DESeq2, using the instructions at
#'  https://bioconductor.org/packages/release/bioc/html/DESeq2.html
#' }
#' @param min.pct  only test genes that are detected in a minimum fraction of
#' min.pct cells in either of the two populations. Meant to speed up the function
#' by not testing genes that are very infrequently expressed. Default is 0.01
#' @param min.diff.pct  only test genes that show a minimum difference in the
#' fraction of detection between the two groups. Set to -Inf by default
#' @param verbose Print a progress bar once expression testing begins
#' @param only.pos Only return positive markers (FALSE by default)
#' @param max.cells.per.ident Down sample each identity class to a max number.
#' Default is no downsampling. Not activated by default (set to Inf)
#' @param random.seed Random seed for downsampling
#' @param latent.vars Variables to test, used only when \code{test.use} is one of
#' 'LR', 'negbinom', 'poisson', or 'MAST'
#' @param min.cells.feature Minimum number of cells expressing the feature in at least one
#' of the two groups, currently only used for poisson and negative binomial tests
#' @param min.cells.group Minimum number of cells in one of the groups
#' @param fc.results data.frame from FoldChange
#' @param densify Convert the sparse matrix to a dense form before running the 
#' DE test. This can provide speedups but might require higher memory; default is FALSE
#'
#'
#' @importFrom Matrix rowMeans
#' @importFrom stats p.adjust
#'
#' @rdname FindMarkers
#' @concept differential_expression
#' @export
#' @method FindMarkers default
#'
FindMarkers.default <- function(
  object,
  slot = "data",
  cells.1 = NULL,
  cells.2 = NULL,
  features = NULL,
  logfc.threshold = 0.1,
  test.use = "wilcox",
  min.pct = 0.01,
  min.diff.pct = -Inf,
  verbose = TRUE,
  only.pos = FALSE,
  max.cells.per.ident = Inf,
  random.seed = 1,
  latent.vars = NULL,
  min.cells.feature = 3,
  min.cells.group = 3,
  fc.results = NULL,
  densify = FALSE,
  ...
) {
  ValidateCellGroups(
    object = object,
    cells.1 = cells.1,
    cells.2 = cells.2,
    min.cells.group = min.cells.group
  )
  features <- features %||% rownames(x = object)
  # reset parameters so no feature filtering is performed
  if (test.use %in% DEmethods_noprefilter()) {
    features <- rownames(x = object)
    min.diff.pct <- -Inf
    logfc.threshold <- 0
  }
  # feature selection (based on percentages)
  alpha.min <- pmax(fc.results$pct.1, fc.results$pct.2)
  names(x = alpha.min) <- rownames(x = fc.results)
  features <- names(x = which(x = alpha.min >= min.pct))
  if (length(x = features) == 0) {
    warning("No features pass min.pct threshold; returning empty data.frame")
    return(fc.results[features, ])
  }
  alpha.diff <- alpha.min - pmin(fc.results$pct.1, fc.results$pct.2)
  features <- names(
    x = which(x = alpha.min >= min.pct & alpha.diff >= min.diff.pct)
  )
  if (length(x = features) == 0) {
    warning("No features pass min.diff.pct threshold; returning empty data.frame")
    return(fc.results[features, ])
  }

  # feature selection (based on logFC)
  if (slot != "scale.data") {
    total.diff <- fc.results[, 1] #first column is logFC
    names(total.diff) <- rownames(fc.results)
    features.diff <- if (only.pos) {
      names(x = which(x = total.diff >= logfc.threshold))
    } else {
      names(x = which(x = abs(x = total.diff) >= logfc.threshold))
    }
    features <- intersect(x = features, y = features.diff)
    if (length(x = features) == 0) {
      warning("No features pass logfc.threshold threshold; returning empty data.frame")
      return(fc.results[features, ])
    }
  }

  # subsample cell groups if they are too large
  if (max.cells.per.ident < Inf) {
    set.seed(seed = random.seed)
    if (length(x = cells.1) > max.cells.per.ident) {
      cells.1 <- sample(x = cells.1, size = max.cells.per.ident)
    }
    if (length(x = cells.2) > max.cells.per.ident) {
      cells.2 <- sample(x = cells.2, size = max.cells.per.ident)
    }
    if (!is.null(x = latent.vars)) {
      latent.vars <- latent.vars[c(cells.1, cells.2), , drop = FALSE]
    }
  }
  if (inherits(x = object, what = "IterableMatrix")){
    if(test.use != "wilcox"){
      stop("Differential expression with BPCells currently only supports the 'wilcox' method.",
           " Please rerun with test.use = 'wilcox'")
    }
    if (BPCells::storage_order(object) == "col" && isTRUE(getOption('Seurat.warn.findmarkers.bpcells.colmajor'))) {
      warning(paste("Column-major order detected; FindMarkers requires row-major order.", 
              "Consider first running BPCells::transpose_storage_order().",
              "This message will be shown once per session.", sep = "\n"),
              immediate. = TRUE, call. = FALSE)
      options(Seurat.warn.findmarkers.bpcells.colmajor = FALSE)
    }
    data.use <- object[features, c(cells.1, cells.2), drop = FALSE]
    groups <- c(rep("foreground", length(cells.1)), rep("background", length(cells.2)))
    de.results <- suppressMessages(
      BPCells::marker_features(data.use, group = groups, method = "wilcoxon")
    )
    de.results <- subset(de.results, de.results$foreground == "foreground")
    de.results <- data.frame(feature = de.results$feature,
                             p_val = de.results$p_val_raw)
    rownames(de.results) <- de.results$feature
    de.results$feature <- NULL
  } else {
    de.results <- PerformDE(
      object = object,
      cells.1 = cells.1,
      cells.2 = cells.2,
      features = features,
      test.use = test.use,
      verbose = verbose,
      min.cells.feature = min.cells.feature,
      latent.vars = latent.vars,
      densify = densify,
      ...
    )
  }
  de.results <- cbind(de.results, fc.results[rownames(x = de.results), , drop = FALSE])
  if (only.pos) {
    de.results <- de.results[de.results[, 2] > 0, , drop = FALSE]
  }
  if (test.use %in% DEmethods_nocorrect()) {
    de.results <- de.results[order(-de.results$power, -de.results[, 1]), ]
  } else {
    de.results <- de.results[order(de.results$p_val, -abs(de.results$pct.1-de.results$pct.2)), ]
    de.results$p_val_adj = p.adjust(
      p = de.results$p_val,
      method = "bonferroni",
      n = nrow(x = object)
    )
  }
  return(de.results)
}


#' @param fc.slot Slot used to calculate fold-change - will also affect the 
#' default for \code{mean.fxn}, see below for more details. 
#' @param pseudocount.use Pseudocount to add to averaged expression values when
#' calculating logFC. 1 by default.
#' @param norm.method Normalization method for fold change calculation when
#' \code{slot} is \dQuote{\code{data}}
#' @param mean.fxn Function to use for fold change or average difference calculation.
#' The default depends on the the value of \code{fc.slot}:
#' \itemize{
#'  \item{"counts"} : difference in the log of the mean counts, with pseudocount.
#'  \item{"data"} : difference in the log of the average exponentiated data, with pseudocount.
#'  This adjusts for differences in sequencing depth between cells, and assumes that "data"
#'  has been log-normalized.
#'  \item{"scale.data"} : difference in the means of scale.data.
#' }
#' @param fc.name Name of the fold change, average difference, or custom function column
#' in the output data.frame. If NULL, the fold change column will be named
#' according to the logarithm base (eg, "avg_log2FC"), or if using the scale.data
#' slot "avg_diff".
#' @param base The base with respect to which logarithms are computed.
#'
#' @rdname FindMarkers
#' @concept differential_expression
#' @export
#' @method FindMarkers Assay
#'
FindMarkers.Assay <- function(
  object,
  slot = "data",
  cells.1 = NULL,
  cells.2 = NULL,
  features = NULL,
  test.use = "wilcox",
  fc.slot = "data",
  pseudocount.use = 1,
  norm.method = NULL,
  mean.fxn = NULL,
  fc.name = NULL,
  base = 2,
  ...
) {
  data.slot <- ifelse(
    test = test.use %in% DEmethods_counts(),
    yes = "counts",
    no = slot
  )
  if (length(x = Layers(object = object, search = slot)) > 1) {
    stop(slot, " layers are not joined. Please run JoinLayers")
  }
  data.use <- LayerData(object = object, layer = data.slot)
  fc.results <- FoldChange(
    object = object,
    slot = fc.slot,
    cells.1 = cells.1,
    cells.2 = cells.2,
    features = features,
    pseudocount.use = pseudocount.use,
    mean.fxn = mean.fxn,
    fc.name = fc.name,
    base = base,
    norm.method = norm.method
  )
  de.results <- FindMarkers(
    object = data.use,
    cells.1 = cells.1,
    cells.2 = cells.2,
    features = features,
    test.use = test.use,
    fc.results = fc.results,
    ...
  )
  return(de.results)
}

#' @method FindMarkers StdAssay
#' @export
#'
FindMarkers.StdAssay <- FindMarkers.Assay

#' @noRd
#'
ValidateSCTFindMarkers <- function(object, recorrect_umi = TRUE) {
  if (!recorrect_umi || length(x = levels(x = object)) <= 1) {
    return(invisible(x = NULL))
  }
  cell_attributes <- SCTResults(object = object, slot = "cell.attributes")
  observed_median_umis <- lapply(
    X = cell_attributes,
    FUN = function(x) median(x[, "umi"])
  )
  model.list <- slot(object = object, "SCTModel.list")
  median_umi.status <- lapply(
    X = model.list,
    FUN = function(x) {
      return(tryCatch(
        expr = slot(object = x, name = 'median_umi'),
        error = function(...) {
          return(NULL)
        }
      ))
    }
  )
  if (any(is.null(unlist(median_umi.status)))) {
    stop(
      "SCT assay does not contain median UMI information.",
      "Run `PrepSCTFindMarkers()` before running `FindMarkers()` or invoke `FindMarkers(recorrect_umi=FALSE)`."
    )
  }
  model_median_umis <- SCTResults(object = object, slot = "median_umi")
  min_median_umi <- min(unlist(x = observed_median_umis))
  if (any(unlist(model_median_umis) != min_median_umi)) {
    stop("Object contains multiple models with unequal library sizes. Run `PrepSCTFindMarkers()` before running `FindMarkers()`.")
  }
  return(invisible(x = NULL))
}

#' @param recorrect_umi Recalculate corrected UMI counts using minimum of the 
#' median UMIs when performing DE using multiple SCT objects; default is TRUE
#'
#' @rdname FindMarkers
#' @concept differential_expression
#' @export
#' @method FindMarkers SCTAssay
#'
FindMarkers.SCTAssay <- function(
  object,
  cells.1 = NULL,
  cells.2 = NULL,
  features = NULL,
  test.use = "wilcox",
  pseudocount.use = 1,
  slot = "data",
  fc.slot = "data",
  mean.fxn = NULL,
  fc.name = NULL,
  base = 2,
  recorrect_umi = TRUE,
  ...
) {
  data.slot <- ifelse(
    test = test.use %in% DEmethods_counts(),
    yes = "counts",
    no = slot
  )
  if (test.use %in% DEmethods_counts()){
    # set slot to counts
    if (slot != "counts") {
      message(paste0("Setting slot to counts for ", test.use, " (counts based test: "))
      slot <- "counts"
    }
  }
  ValidateSCTFindMarkers(object = object, recorrect_umi = recorrect_umi)

  data.use <-  GetAssayData(object = object, layer = data.slot)
  # Default assumes the input is log1p(corrected counts)
  default.mean.fxn <- function(x) {
    return(log(x = (rowSums(x = expm1(x = x)) + pseudocount.use)/NCOL(x), base = base))
  }
  mean.fxn <- mean.fxn %||% switch(
    EXPR = fc.slot,
    "counts" = function(x) {
      return(log(x = (rowSums(x = x) + pseudocount.use)/NCOL(x), base = base))
    },
    "scale.data" = rowMeans,
    default.mean.fxn
  )
  fc.results <- FoldChange(
    object = object,
    slot = fc.slot,
    cells.1 = cells.1,
    cells.2 = cells.2,
    features = features,
    pseudocount.use = pseudocount.use,
    mean.fxn = mean.fxn,
    fc.name = fc.name,
    base = base
  )
  de.results <- FindMarkers(
    object = data.use,
    cells.1 = cells.1,
    cells.2 = cells.2,
    features = features,
    test.use = test.use,
    fc.results = fc.results,
    ...
  )
  return(de.results)
}


#' @importFrom Matrix rowMeans
#' @rdname FindMarkers
#' @concept differential_expression
#' @export
#' @method FindMarkers DimReduc
#'
FindMarkers.DimReduc <- function(
  object,
  cells.1 = NULL,
  cells.2 = NULL, 
  features = NULL,
  logfc.threshold = 0.1,
  test.use = "wilcox",
  min.pct = 0.01,
  min.diff.pct = -Inf,
  verbose = TRUE,
  only.pos = FALSE,
  max.cells.per.ident = Inf,
  random.seed = 1,
  latent.vars = NULL,
  min.cells.feature = 3,
  min.cells.group = 3,
  densify = FALSE,
  mean.fxn = rowMeans,
  fc.name = NULL,
  ...

) {
  if (test.use %in% DEmethods_counts()) {
    stop("The following tests cannot be used for differential expression on a reduction as they assume a count model: ",
         paste(DEmethods_counts(), collapse=", "))
  }
  data <- t(x = Embeddings(object = object))
  ValidateCellGroups(
    object = data,
    cells.1 = cells.1,
    cells.2 = cells.2,
    min.cells.group = min.cells.group
  )
  features <- features %||% rownames(x = data)
  # reset parameters so no feature filtering is performed
  if (test.use %in% DEmethods_noprefilter()) {
    features <- rownames(x = data)
    min.diff.pct <- -Inf
    logfc.threshold <- 0
  }
  fc.results <- FoldChange(
    object = object,
    cells.1 = cells.1,
    cells.2 = cells.2,
    features = features,
    mean.fxn = mean.fxn,
    fc.name = fc.name
  )
  # subsample cell groups if they are too large
  if (max.cells.per.ident < Inf) {
    set.seed(seed = random.seed)
    if (length(x = cells.1) > max.cells.per.ident) {
      cells.1 <- sample(x = cells.1, size = max.cells.per.ident)
    }
    if (length(x = cells.2) > max.cells.per.ident) {
      cells.2 <- sample(x = cells.2, size = max.cells.per.ident)
    }
    if (!is.null(x = latent.vars)) {
      latent.vars <- latent.vars[c(cells.1, cells.2), , drop = FALSE]
    }
  }
  de.results <- PerformDE(
    object = data,
    cells.1 = cells.1,
    cells.2 = cells.2,
    features = features,
    test.use = test.use,
    verbose = verbose,
    min.cells.feature = min.cells.feature,
    latent.vars = latent.vars,
    densify = densify,
    ...
  )
  de.results <- cbind(de.results, fc.results)
  if (only.pos) {
    de.results <- de.results[de.results$avg_diff > 0, , drop = FALSE]
  }
  if (test.use %in% DEmethods_nocorrect()) {
    de.results <- de.results[order(-de.results$power, -de.results$avg_diff), ]
  } else {
    de.results <- de.results[order(de.results$p_val, -de.results$avg_diff), ]
    de.results$p_val_adj = p.adjust(
      p = de.results$p_val,
      method = "bonferroni",
      n = ncol(x = object)
    )
  }
  return(de.results)
}

#' @param ident.1 Identity class to define markers for; pass an object of class
#' \code{phylo} or 'clustertree' to find markers for a node in a cluster tree;
#' passing 'clustertree' requires \code{\link{BuildClusterTree}} to have been run
#' @param ident.2 A second identity class for comparison; if \code{NULL},
#' use all other cells for comparison; if an object of class \code{phylo} or
#' 'clustertree' is passed to \code{ident.1}, must pass a node to find markers for
#' @param group.by Regroup cells into a different identity class prior to 
#' performing differential expression (see example); \code{"ident"} to use Idents
#' @param subset.ident Subset a particular identity class prior to regrouping. 
#' Only relevant if group.by is set (see example)
#' @param assay Assay to use in differential expression testing
#' @param reduction Reduction to use in differential expression testing - will 
#' test for DE on cell embeddings
#'
#' @rdname FindMarkers
#' @concept differential_expression
#' @export
#' @method FindMarkers Seurat
#'
FindMarkers.Seurat <- function(
  object,
  ident.1 = NULL,
  ident.2 = NULL,
  latent.vars = NULL,
  group.by = NULL,
  subset.ident = NULL,
  assay = NULL,
  reduction = NULL,
  ...
) {
  if (!is.null(x = group.by) && !identical(x = group.by, y = "ident")) {
    if (!is.null(x = subset.ident)) {
      object <- subset(x = object, idents = subset.ident)
    }
    if (length(x = group.by) == 1 && ! group.by %in% colnames(x = object@meta.data)) {
      stop("'", group.by, "' not found in object metadata")
    }
    Idents(object = object) <- group.by
  }
  if (!is.null(x = assay) && !is.null(x = reduction)) {
    stop("Please only specify either assay or reduction.")
  }
  if (length(x = ident.1) == 0) {
    stop("At least 1 ident must be specified in `ident.1`")
  }

  # select which data to use
  if (is.null(x = reduction)) {
    assay <- assay %||% DefaultAssay(object = object)
    data.use <- object[[assay]]
    cellnames.use <-  colnames(x = data.use)
  } else {
    data.use <- object[[reduction]]
    cellnames.use <- rownames(x = data.use)
  }

  cells <- IdentsToCells(
    object = object,
    ident.1 = ident.1,
    ident.2 = ident.2,
    cellnames.use = cellnames.use
  )
  cells <- sapply(
    X = cells,
    FUN = intersect,
    y = cellnames.use,
    simplify = FALSE,
    USE.NAMES = TRUE
  )
  if (!all(vapply(X = cells, FUN = length, FUN.VALUE = integer(length = 1L)))) {
    abort(
      message = "Cells in one or both identity groups are not present in the data requested"
    )
  }

  # fetch latent.vars
  if (!is.null(x = latent.vars)) {
    latent.vars <- FetchData(
      object = object,
      vars = latent.vars,
      cells = c(cells$cells.1, cells$cells.2)
    )
  }

  # check normalization method
  norm.command <- paste0("NormalizeData.", assay)
  norm.method <- if (norm.command %in% Command(object = object) && is.null(x = reduction)) {
    Command(
      object = object,
      command = norm.command,
      value = "normalization.method"
    )
  } else if (length(x = intersect(x = c("FindIntegrationAnchors", "FindTransferAnchors"), y = Command(object = object)))) {
    command <- intersect(x = c("FindIntegrationAnchors", "FindTransferAnchors"), y = Command(object = object))[1]
    Command(
      object = object,
      command = command,
      value = "normalization.method"
      )
    } else {
    NULL
  }

  de.results <- FindMarkers(
    object = data.use,
    latent.vars = latent.vars,
    cells.1 = cells$cells.1,
    cells.2 = cells$cells.2,
    norm.method = norm.method,
    ...
  )

  return(de.results)
}

#' @param cells.1 Vector of cell names belonging to group 1
#' @param cells.2 Vector of cell names belonging to group 2
#' @param features Features to calculate fold change for.
#' If NULL, use all features
#' @importFrom Matrix rowSums
#' @rdname FoldChange
#' @concept differential_expression
#' @export
#' @method FoldChange default
FoldChange.default <- function(
  object,
  cells.1,
  cells.2,
  mean.fxn,
  fc.name,
  features = NULL,
  ...
) {
  features <- features %||% rownames(x = object)
  thresh.min <- 0
  # Subset each group's slice a single time and reuse it for both the
  # percent-expressed and mean-expression computations.
  mat.1 <- object[features, cells.1, drop = FALSE]
  mat.2 <- object[features, cells.2, drop = FALSE]
  # Calculate percent expressed
  pct.1 <- round(
    x = rowSums(x = mat.1 > thresh.min) / length(x = cells.1),
    digits = 3
  )
  pct.2 <- round(
    x = rowSums(x = mat.2 > thresh.min) / length(x = cells.2),
    digits = 3
  )
  # Calculate fold change
  data.1 <- mean.fxn(mat.1)
  data.2 <- mean.fxn(mat.2)
  fc <- (data.1 - data.2)
  fc.results <- as.data.frame(x = cbind(fc, pct.1, pct.2))
  colnames(fc.results) <- c(fc.name, "pct.1", "pct.2")
  return(fc.results)
}

#' @param norm.method Normalization method for mean function selection
#' when \code{slot} is \dQuote{\code{data}}
#'
#' @importFrom Matrix rowMeans
#' @importFrom Matrix rowSums
#' @rdname FoldChange
#' @concept differential_expression
#' @export
#' @method FoldChange Assay
FoldChange.Assay <- function(
  object,
  cells.1,
  cells.2,
  features = NULL,
  slot = "data",
  pseudocount.use = 1,
  fc.name = NULL,
  mean.fxn = NULL,
  base = 2,
  norm.method = NULL,
  ...
) {
  data <- GetAssayData(object = object, layer = slot)
  # By default run as if LogNormalize is done
  log1pdata.mean.fxn <- function(x) {
    return(log(x = (rowSums(x = expm1(x = x)) + pseudocount.use)/NCOL(x), base = base))
  }
  scaledata.mean.fxn <- rowMeans
  counts.mean.fxn <- function(x) {
    return(log(x = (rowSums(x = x) + pseudocount.use)/NCOL(x), base = base))
  }
  if (!is.null(x = norm.method)) {
    # For anything apart from log normalization set to rowMeans
    if (norm.method != "LogNormalize") {
      new.mean.fxn <- counts.mean.fxn
    } else {
      new.mean.fxn <- counts.mean.fxn
      if (slot == "data") {
        new.mean.fxn <- log1pdata.mean.fxn
      }  else if (slot == "scale.data") {
        new.mean.fxn <- scaledata.mean.fxn
      }
    }
  } else {
    # If no normalization method is passed use slots to decide mean function
    new.mean.fxn <- switch(
      EXPR = slot,
      'data' = log1pdata.mean.fxn,
      'scale.data' = scaledata.mean.fxn,
      'counts' = counts.mean.fxn,
      log1pdata.mean.fxn
    )
  }
  mean.fxn <- mean.fxn %||% new.mean.fxn
  # Omit the decimal value of e from the column name if base == exp(1)
  base.text <- ifelse(
    test = base == exp(1),
    yes = "",
    no = base
  )
  fc.name <- fc.name %||% ifelse(
    test = slot == "scale.data",
    yes = "avg_diff",
    no = paste0("avg_log", base.text, "FC")
  )
  FoldChange(
    object = data,
    cells.1 = cells.1,
    cells.2 = cells.2,
    features = features,
    mean.fxn = mean.fxn,
    fc.name = fc.name
  )
}

#' @method FoldChange StdAssay
#' @export
#'
FoldChange.StdAssay <- FoldChange.Assay

#' @importFrom Matrix rowMeans
#' @importFrom Matrix rowSums
#' @rdname FoldChange
#' @concept differential_expression
#' @export
#' @method FoldChange SCTAssay
FoldChange.SCTAssay <- function(
    object,
    cells.1,
    cells.2,
    features = NULL,
    slot = "data",
    pseudocount.use = 1,
    fc.name = NULL,
    mean.fxn = NULL,
    base = 2,
    ...
) {
  pseudocount.use <- pseudocount.use %||% 1
  data <- GetAssayData(object = object, layer = slot)
  default.mean.fxn <- function(x) {
    return(log(x = (rowSums(x = expm1(x = x)) + pseudocount.use)/NCOL(x), base = base))
  }
  mean.fxn <- mean.fxn %||% switch(
    EXPR = slot,
    'data' = default.mean.fxn,
    'scale.data' = rowMeans,
    'counts' = function(x) {
      return(log(x = (rowSums(x = x) + pseudocount.use)/NCOL(x), base = base))
    },
    default.mean.fxn
  )
  # Omit the decimal value of e from the column name if base == exp(1)
  base.text <- ifelse(
    test = base == exp(1),
    yes = "",
    no = base
  )
  fc.name <- fc.name %||% ifelse(
    test = slot == "scale.data",
    yes = "avg_diff",
    no = paste0("avg_log", base.text, "FC")
  )
  FoldChange(
    object = data,
    cells.1 = cells.1,
    cells.2 = cells.2,
    features = features,
    mean.fxn = mean.fxn,
    fc.name = fc.name
  )
}


#' @importFrom Matrix rowMeans
#' @rdname FoldChange
#' @concept differential_expression
#' @export
#' @method FoldChange DimReduc
FoldChange.DimReduc <- function(
  object,
  cells.1,
  cells.2,
  features = NULL,
  slot = NULL,
  pseudocount.use = 1,
  fc.name = NULL,
  mean.fxn = NULL,
  ...
) {
  mean.fxn <- mean.fxn %||% rowMeans
  fc.name <- fc.name %||% "avg_diff"
  data <- t(x = Embeddings(object = object))
  features <- features %||% rownames(x = data)
  # Calculate avg difference
  data.1 <- mean.fxn(data[features, cells.1, drop = FALSE])
  data.2 <- mean.fxn(data[features, cells.2, drop = FALSE])
  fc <- (data.1 - data.2)
  fc.results <- data.frame(fc)
  colnames(fc.results) <- fc.name
  return(fc.results)
}

#' @param ident.1 Identity class to calculate fold change for; pass an object of class
#' \code{phylo} or 'clustertree' to calculate fold change for a node in a cluster tree;
#' passing 'clustertree' requires \code{\link{BuildClusterTree}} to have been run
#' @param ident.2 A second identity class for comparison; if \code{NULL},
#' use all other cells for comparison; if an object of class \code{phylo} or
#' 'clustertree' is passed to \code{ident.1}, must pass a node to calculate fold change for
#' @param reduction Reduction to use - will calculate average difference on cell embeddings
#' @param group.by Regroup cells into a different identity class prior to
#' calculating fold change (see example in \code{\link{FindMarkers}})
#' @param subset.ident Subset a particular identity class prior to regrouping.
#' Only relevant if group.by is set (see example in \code{\link{FindMarkers}})
#' @param assay Assay to use in fold change calculation
#' @param slot Slot to pull data from
#' @param pseudocount.use Pseudocount to add to averaged expression values when
#' calculating logFC.
#' @param mean.fxn Function to use for fold change or average difference calculation
#' @param base The base with respect to which logarithms are computed.
#' @param fc.name Name of the fold change, average difference, or custom function column
#' in the output data.frame
#'
#' @rdname FoldChange
#' @concept differential_expression
#' @export
#' @method FoldChange Seurat
FoldChange.Seurat <- function(
  object,
  ident.1 = NULL,
  ident.2 = NULL,
  group.by = NULL,
  subset.ident = NULL,
  assay = NULL,
  slot = 'data',
  reduction = NULL,
  features = NULL,
  pseudocount.use = 1,
  mean.fxn = NULL,
  base = 2,
  fc.name = NULL,
  ...
) {
  if (!is.null(x = group.by)) {
    if (!is.null(x = subset.ident)) {
      object <- subset(x = object, idents = subset.ident)
    }
    Idents(object = object) <- group.by
  }
  if (!is.null(x = assay) && !is.null(x = reduction)) {
    stop("Please only specify either assay or reduction.")
  }
  # select which data to use
  if (is.null(x = reduction)) {
    assay <- assay %||% DefaultAssay(object = object)
    data.use <- object[[assay]]
    cellnames.use <-  colnames(x = data.use)
  } else {
    data.use <- object[[reduction]]
    cellnames.use <- rownames(data.use)
  }
  cells <- IdentsToCells(
    object = object,
    ident.1 = ident.1,
    ident.2 = ident.2,
    cellnames.use = cellnames.use
  )
  # check normalization method
  norm.command <- paste0("NormalizeData.", assay)
  norm.method <- if (norm.command %in% Command(object = object) && is.null(x = reduction)) {
    Command(
      object = object,
      command = norm.command,
      value = "normalization.method"
    )
  } else if (length(x = intersect(x = c("FindIntegrationAnchors", "FindTransferAnchors"), y = Command(object = object)))) {
    command <- intersect(x = c("FindIntegrationAnchors", "FindTransferAnchors"), y = Command(object = object))[1]
    Command(
      object = object,
      command = command,
      value = "normalization.method"
    )
  } else {
    NULL
  }
  fc.results <- FoldChange(
    object = data.use,
    cells.1 = cells$cells.1,
    cells.2 = cells$cells.2,
    features = features,
    slot = slot,
    pseudocount.use = pseudocount.use,
    mean.fxn = mean.fxn,
    base = base,
    fc.name = fc.name,
    norm.method = norm.method
  )
  return(fc.results)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Internal
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# internal function to calculate AUC values
#' @importFrom pbapply pblapply
#
AUCMarkerTest <- function(data1, data2, mygenes, print.bar = TRUE) {
  myAUC <- unlist(x = lapply(
    X = mygenes,
    FUN = function(x) {
      return(DifferentialAUC(
        x = as.numeric(x = data1[x, ]),
        y = as.numeric(x = data2[x, ])
      ))
    }
  ))
  myAUC[is.na(x = myAUC)] <- 0
  iterate.fxn <- ifelse(test = print.bar, yes = pblapply, no = lapply)
  avg_diff <- unlist(x = iterate.fxn(
    X = mygenes,
    FUN = function(x) {
      return(
        ExpMean(
          x = as.numeric(x = data1[x, ])
        ) - ExpMean(
          x = as.numeric(x = data2[x, ])
        )
      )
    }
  ))
  toRet <- data.frame(cbind(myAUC, avg_diff), row.names = mygenes)
  toRet <- toRet[rev(x = order(toRet$myAUC)), ]
  return(toRet)
}

#internal function to run mcdavid et al. DE test
#
#' @importFrom stats sd dnorm
#
bimodLikData <- function(x, xmin = 0) {
  x1 <- x[x <= xmin]
  x2 <- x[x > xmin]
  xal <- MinMax(
    data = length(x = x2) / length(x = x),
    min = 1e-5,
    max = (1 - 1e-5)
  )
  likA <- length(x = x1) * log(x = 1 - xal)
  if (length(x = x2) < 2) {
    mysd <- 1
  } else {
    mysd <- sd(x = x2)
  }
  likB <- length(x = x2) *
    log(x = xal) +
    sum(dnorm(x = x2, mean = mean(x = x2), sd = mysd, log = TRUE))
  return(likA + likB)
}

# returns tests that do not support feature pre-filtering
DEmethods_noprefilter <- function() {
  c("DESeq2")
}

# returns tests that support latent variables (latent.vars)
DEmethods_latent <- function() {
  c('negbinom', 'poisson', 'MAST', "LR")
}

# returns tests that require CheckDots
DEmethods_checkdots <- function() {
  c('wilcox', 'wilcox_limma', 'MAST', 'DESeq2')
}

# returns tests that do not use Bonferroni correction on the DE results
DEmethods_nocorrect <- function() {
  c('roc')
}

# returns tests that require count data
DEmethods_counts <- function() {
  c("negbinom", "poisson", "DESeq2")
}

# Differential expression using DESeq2
#
# Identifies differentially expressed genes between two groups of cells using
# DESeq2
#
# @references Love MI, Huber W and Anders S (2014). "Moderated estimation of
# fold change and dispersion for RNA-seq data with DESeq2." Genome Biology.
# https://bioconductor.org/packages/release/bioc/html/DESeq2.html
# @param data.use Data matrix to test
# @param cells.1 Group 1 cells
# @param cells.2 Group 2 cells
# @param verbose Print a progress bar
# @param ... Extra parameters to pass to DESeq2::results
# @return Returns a p-value ranked matrix of putative differentially expressed
# genes.
#
# @details
# This test does not support pre-filtering of genes based on average difference
# (or percent detection rate) between cell groups. However, genes may be
# pre-filtered based on their minimum detection rate (min.pct) across both cell
# groups. To use this method, please install DESeq2, using the instructions at
#  https://bioconductor.org/packages/release/bioc/html/DESeq2.html
#
# @export
#
# @examples
# \dontrun{
#   data("pbmc_small")
#   pbmc_small
#   DESeq2DETest(pbmc_small, cells.1 = WhichCells(object = pbmc_small, idents = 1),
#               cells.2 = WhichCells(object = pbmc_small, idents = 2))
# }
#
DESeq2DETest <- function(
  data.use,
  cells.1,
  cells.2,
  verbose = TRUE,
  ...
) {
  if (!requireNamespace("DESeq2", quietly=TRUE)) {
    stop("Please install DESeq2 - learn more at https://bioconductor.org/packages/release/bioc/html/DESeq2.html")
  }
  CheckDots(..., fxns = 'DESeq2::results')
  group.info <- data.frame(row.names = c(cells.1, cells.2))
  group.info[cells.1, "group"] <- "Group1"
  group.info[cells.2, "group"] <- "Group2"
  group.info[, "group"] <- factor(x = group.info[, "group"])
  group.info$wellKey <- rownames(x = group.info)
  dds1 <- DESeq2::DESeqDataSetFromMatrix(
    countData = data.use,
    colData = group.info,
    design = ~ group
  )
  dds1 <- DESeq2::estimateSizeFactors(object = dds1)
  dds1 <- DESeq2::estimateDispersions(object = dds1, fitType = "local")
  dds1 <- DESeq2::nbinomWaldTest(object = dds1)
  res <- DESeq2::results(
    object = dds1,
    contrast = c("group", "Group1", "Group2"),
    alpha = 0.05,
    ...
  )
  to.return <- data.frame(p_val = res$pvalue, row.names = rownames(res))
  return(to.return)
}

# internal function to calculate AUC values
#' @importFrom ROCR prediction performance
#'
DifferentialAUC <- function(x, y) {
  prediction.use <- prediction(
    predictions = c(x, y),
    labels = c(rep(x = 1, length(x = x)), rep(x = 0, length(x = y))),
    label.ordering = 0:1
  )
  perf.use <- performance(prediction.obj = prediction.use, measure = "auc")
  auc.use <- round(x = perf.use@y.values[[1]], digits = 3)
  return(auc.use)
}

#internal function to run mcdavid et al. DE test
#
#' @importFrom stats pchisq
#
DifferentialLRT <- function(x, y, xmin = 0) {
  lrtX <- bimodLikData(x = x)
  lrtY <- bimodLikData(x = y)
  lrtZ <- bimodLikData(x = c(x, y))
  lrt_diff <- 2 * (lrtX + lrtY - lrtZ)
  return(pchisq(q = lrt_diff, df = 3, lower.tail = F))
}

# Likelihood ratio test for zero-inflated data
#
# Identifies differentially expressed genes between two groups of cells using
# the LRT model proposed in McDavid et al, Bioinformatics, 2013
#
# @inheritParams FindMarkers
# @param object Seurat object
# @param cells.1 Group 1 cells
# @param cells.2 Group 2 cells
# @param assay.type Type of assay to fetch data for (default is RNA)
# @return Returns a p-value ranked matrix of putative differentially expressed
# genes.
#
#' @importFrom pbapply pbsapply
#' @importFrom future.apply future_sapply
#' @importFrom future nbrOfWorkers
#
# @export
# @examples
# data("pbmc_small")
# pbmc_small
# DiffExpTest(pbmc_small, cells.1 = WhichCells(object = pbmc_small, idents = 1),
#             cells.2 = WhichCells(object = pbmc_small, idents = 2))
#
DiffExpTest <- function(
  data.use,
  cells.1,
  cells.2,
  verbose = TRUE
) {
  my.sapply <- ifelse(
    test = verbose && nbrOfWorkers() == 1,
    yes = pbsapply,
    no = future_sapply
  )
  p_val <- unlist(
    x = my.sapply(
      X = 1:nrow(x = data.use),
      FUN = function(x) {
        return(DifferentialLRT(
          x = as.numeric(x = data.use[x, cells.1]),
          y = as.numeric(x = data.use[x, cells.2])
        ))
      }
    )
  )
  to.return <- data.frame(p_val, row.names = rownames(x = data.use))
  return(to.return)
}

# Differential expression testing using Student's t-test
#
# Identify differentially expressed genes between two groups of cells using
# the Student's t-test
#
# @return Returns a p-value ranked matrix of putative differentially expressed
# genes.
#
#' @importFrom stats t.test
#' @importFrom pbapply pbsapply
#' @importFrom future.apply future_sapply
#' @importFrom future nbrOfWorkers
#
# @export
#
# @examples
# data("pbmc_small")
# pbmc_small
# DiffTTest(pbmc_small, cells.1 = WhichCells(object = pbmc_small, idents = 1),
#             cells.2 = WhichCells(object = pbmc_small, idents = 2))
DiffTTest <- function(
  data.use,
  cells.1,
  cells.2,
  verbose = TRUE
) {
  my.sapply <- ifelse(
    test = verbose && nbrOfWorkers() == 1,
    yes = pbsapply,
    no = future_sapply
  )
  p_val <- unlist(
    x = my.sapply(
      X = 1:nrow(data.use),
      FUN = function(x) {
        t.test(x = data.use[x, cells.1], y = data.use[x, cells.2])$p.value
      }
    )
  )
  to.return <- data.frame(p_val,row.names = rownames(x = data.use))
  return(to.return)
}

# Tests for UMI-count based data
#
# Identifies differentially expressed genes between two groups of cells using
# either a negative binomial or poisson generalized linear model
#
# @param data.use Data to test
# @param cells.1 Group 1 cells
# @param cells.2 Group 2 cells
# @param min.cells Minimum number of cells threshold
# @param latent.vars Latent variables to test
# @param test.use parameterizes the glm
# @param verbose Print progress bar
#
# @return Returns a p-value ranked matrix of putative differentially expressed
# genes.
#
#' @importFrom MASS glm.nb
#' @importFrom pbapply pbsapply
#' @importFrom stats var as.formula
#' @importFrom future.apply future_sapply
#' @importFrom future nbrOfWorkers
#'
# @export
#
# @examples
# data("pbmc_small")
# pbmc_small
# # Note, not recommended for particularly small datasets - expect warnings
# NegBinomDETest(pbmc_small, cells.1 = WhichCells(object = pbmc_small, idents = 1),
#             cells.2 = WhichCells(object = pbmc_small, idents = 2))
#
GLMDETest <- function(
  data.use,
  cells.1,
  cells.2,
  min.cells = 3,
  latent.vars = NULL,
  test.use = NULL,
  verbose = TRUE
) {
  group.info <- data.frame(
    group = rep(
      x = c('Group1', 'Group2'),
      times = c(length(x = cells.1), length(x = cells.2))
    )
  )
  rownames(group.info) <- c(cells.1, cells.2)
  group.info[, "group"] <- factor(x = group.info[, "group"])
  latent.vars <- if (is.null(x = latent.vars)) {
    group.info
  } else {
    cbind(x = group.info, latent.vars)
  }
  latent.var.names <- colnames(x = latent.vars)
  my.sapply <- ifelse(
    test = verbose && nbrOfWorkers() == 1,
    yes = pbsapply,
    no = future_sapply
  )
  p_val <- unlist(
    x = my.sapply(
      X = 1:nrow(data.use),
      FUN = function(x) {
        latent.vars[, "GENE"] <- as.numeric(x = data.use[x, ])
        # check that gene is expressed in specified number of cells in one group
        if (sum(latent.vars$GENE[latent.vars$group == "Group1"] > 0) < min.cells &&
            sum(latent.vars$GENE[latent.vars$group == "Group2"] > 0) < min.cells) {
          warning(paste0(
            "Skipping gene --- ",
            x,
            ". Fewer than ",
            min.cells,
            " cells in both clusters."
          ))
          return(2)
        }
        # check that variance between groups is not 0
        if (var(x = latent.vars$GENE) == 0) {
          warning(paste0(
            "Skipping gene -- ",
            x,
            ". No variance in expression between the two clusters."
          ))
          return(2)
        }
        fmla <- as.formula(object = paste(
          "GENE ~",
          paste(latent.var.names, collapse = "+")
        ))
        p.estimate <- 2
        if (test.use == "negbinom") {
          try(
            expr = p.estimate <- summary(
              object = glm.nb(formula = fmla, data = latent.vars)
            )$coef[2, 4],
            silent = TRUE
          )
          return(p.estimate)
        } else if (test.use == "poisson") {
          return(summary(object = glm(
            formula = fmla,
            data = latent.vars,
            family = "poisson"
          ))$coef[2,4])
        }
      }
    )
  )
  features.keep <- rownames(data.use)
  if (length(x = which(x = p_val == 2)) > 0) {
    features.keep <- features.keep[-which(x = p_val == 2)]
    p_val <- p_val[!p_val == 2]
  }
  to.return <- data.frame(p_val, row.names = features.keep)
  return(to.return)
}

# Helper function for FindMarkers.Seurat and FoldChange.Seurat
# Convert idents to cells
#
#' @importFrom methods is
#
IdentsToCells <- function(
  object,
  ident.1,
  ident.2,
  cellnames.use
) {
  #
  if (is.null(x = ident.1)) {
    stop("Please provide ident.1")
  } else if ((length(x = ident.1) == 1 && ident.1[1] == 'clustertree') || is(object = ident.1, class2 = 'phylo')) {
    if (is.null(x = ident.2)) {
      stop("Please pass a node to 'ident.2' to run FindMarkers on a tree")
    }
    tree <- if (is(object = ident.1, class2 = 'phylo')) {
      ident.1
    } else {
      Tool(object = object, slot = 'BuildClusterTree')
    }
    if (is.null(x = tree)) {
      stop("Please run 'BuildClusterTree' or pass an object of class 'phylo' as 'ident.1'")
    }
    ident.1 <- tree$tip.label[GetLeftDescendants(tree = tree, node = ident.2)]
    ident.2 <- tree$tip.label[GetRightDescendants(tree = tree, node = ident.2)]
  }
  if (length(x = as.vector(x = ident.1)) > 1 &&
      any(as.character(x = ident.1) %in% cellnames.use)) {
    bad.cells <- cellnames.use[which(x = !as.character(x = ident.1) %in% cellnames.use)]
    if (length(x = bad.cells) > 0) {
      stop(paste0("The following cell names provided to ident.1 are not present in the object: ", paste(bad.cells, collapse = ", ")))
    }
  } else {
    ident.1 <- WhichCells(object = object, idents = ident.1)
  }
  # if NULL for ident.2, use all other cells
  if (length(x = as.vector(x = ident.2)) > 1 &&
      any(as.character(x = ident.2) %in% cellnames.use)) {
    bad.cells <- cellnames.use[which(!as.character(x = ident.2) %in% cellnames.use)]
    if (length(x = bad.cells) > 0) {
      stop(paste0("The following cell names provided to ident.2 are not present in the object: ", paste(bad.cells, collapse = ", ")))
    }
  } else {
    if (is.null(x = ident.2)) {
      ident.2 <- setdiff(x = cellnames.use, y = ident.1)
    } else {
      ident.2 <- WhichCells(object = object, idents = ident.2)
    }
  }
  return(list(cells.1 = ident.1, cells.2 = ident.2))
}

# Perform differential expression testing using a logistic regression framework
#
# Constructs a logistic regression model predicting group membership based on a
# given feature and compares this to a null model with a likelihood ratio test.
#
# @param data.use expression matrix
# @param cells.1 Vector of cells in group 1
# @param cells2. Vector of cells in group 2
# @param latent.vars Latent variables to include in model
# @param verbose Print messages
#
#' @importFrom lmtest lrtest
#' @importFrom pbapply pbsapply
#' @importFrom stats as.formula glm
#' @importFrom future.apply future_sapply
#' @importFrom future nbrOfWorkers
#
LRDETest <- function(
  data.use,
  cells.1,
  cells.2,
  latent.vars = NULL,
  verbose = TRUE
) {
  group.info <- data.frame(row.names = c(cells.1, cells.2))
  group.info[cells.1, "group"] <- "Group1"
  group.info[cells.2, "group"] <- "Group2"
  group.info[, "group"] <- factor(x = group.info[, "group"])
  data.use <- data.use[, rownames(group.info), drop = FALSE]
  latent.vars <- latent.vars[rownames(group.info), , drop = FALSE]
  my.sapply <- ifelse(
    test = verbose && nbrOfWorkers() == 1,
    yes = pbsapply,
    no = future_sapply
  )
  p_val <- my.sapply(
    X = 1:nrow(x = data.use),
    FUN = function(x) {
      if (is.null(x = latent.vars)) {
        model.data <- cbind(GENE = data.use[x, ], group.info)
        fmla <- as.formula(object = "group ~ GENE")
        fmla2 <- as.formula(object = "group ~ 1")
      } else {
        model.data <- cbind(GENE = data.use[x, ], group.info, latent.vars)
        fmla <- as.formula(object = paste(
          "group ~ GENE +",
          paste(colnames(x = latent.vars), collapse = "+")
        ))
        fmla2 <- as.formula(object = paste(
          "group ~",
          paste(colnames(x = latent.vars), collapse = "+")
        ))
      }
      model1 <- glm(formula = fmla, data = model.data, family = "binomial")
      model2 <- glm(formula = fmla2, data = model.data, family = "binomial")
      lrtest <- lrtest(model1, model2)
      return(lrtest$Pr[2])
    }
  )
  to.return <- data.frame(p_val, row.names = rownames(data.use))
  return(to.return)
}

# ROC-based marker discovery
#
# Identifies 'markers' of gene expression using ROC analysis. For each gene,
# evaluates (using AUC) a classifier built on that gene alone, to classify
# between two groups of cells.
#
# An AUC value of 1 means that expression values for this gene alone can
# perfectly classify the two groupings (i.e. Each of the cells in cells.1
# exhibit a higher level than each of the cells in cells.2). An AUC value of 0
# also means there is perfect classification, but in the other direction. A
# value of 0.5 implies that the gene has no predictive power to classify the
# two groups.
#
# @return Returns a 'predictive power' (abs(AUC-0.5) * 2) ranked matrix of
# putative differentially expressed genes.
#
# @export
#
# @examples
# data("pbmc_small")
# pbmc_small
# MarkerTest(pbmc_small, cells.1 = WhichCells(object = pbmc_small, idents = 1),
#             cells.2 = WhichCells(object = pbmc_small, idents = 2))
#
MarkerTest <- function(
  data.use,
  cells.1,
  cells.2,
  verbose = TRUE
) {
  to.return <- AUCMarkerTest(
    data1 = data.use[, cells.1, drop = FALSE],
    data2 = data.use[, cells.2, drop = FALSE],
    mygenes = rownames(x = data.use),
    print.bar = verbose
  )
  to.return$power <- abs(x = to.return$myAUC - 0.5) * 2
  return(to.return)
}

# Differential expression using MAST
#
# Identifies differentially expressed genes between two groups of cells using
# a hurdle model tailored to scRNA-seq data. Utilizes the MAST package to run
# the DE testing.
#
# @references Andrew McDavid, Greg Finak and Masanao Yajima (2017). MAST: Model-based
# Analysis of Single Cell Transcriptomics. R package version 1.2.1.
# https://github.com/RGLab/MAST/
#
# @param data.use Data to test
# @param cells.1 Group 1 cells
# @param cells.2 Group 2 cells
# @param latent.vars Confounding variables to adjust for in DE test
# @param verbose print output
# @param \dots Additional parameters to zero-inflated regression (zlm) function
# in MAST
# @details
# To use this method, please install MAST, using instructions at https://github.com/RGLab/MAST/
#
# @return Returns a p-value ranked matrix of putative differentially expressed
# genes.
#
#' @importFrom stats relevel
MASTDETest <- function(
  data.use,
  cells.1,
  cells.2,
  latent.vars = NULL,
  verbose = TRUE,
  ...
) {
  # Check for MAST
  if (!requireNamespace("MAST", quietly = TRUE)) {
    stop("Please install MAST - learn more at https://github.com/RGLab/MAST")
  }
  group.info <- data.frame(row.names = c(cells.1, cells.2))
  latent.vars <- latent.vars %||% group.info
  group.info[cells.1, "group"] <- "Group1"
  group.info[cells.2, "group"] <- "Group2"
  group.info[, "group"] <- factor(x = group.info[, "group"])
  latent.vars.names <- c("condition", colnames(x = latent.vars))
  latent.vars <- cbind(latent.vars, group.info)
  latent.vars$wellKey <- rownames(x = latent.vars)
  fdat <- data.frame(rownames(x = data.use))
  colnames(x = fdat)[1] <- "primerid"
  rownames(x = fdat) <- fdat[, 1]
  sca <- MAST::FromMatrix(
    exprsArray = as.matrix(x = data.use),
    check_sanity = FALSE,
    cData = latent.vars,
    fData = fdat
  )
  cond <- factor(x = SummarizedExperiment::colData(sca)$group)
  cond <- relevel(x = cond, ref = "Group1")
  SummarizedExperiment::colData(sca)$condition <- cond
  fmla <- as.formula(
    object = paste0(" ~ ", paste(latent.vars.names, collapse = "+"))
  )
  zlmCond <- MAST::zlm(formula = fmla, sca = sca, ...)
  summaryCond <- MAST::summary(object = zlmCond, doLRT = 'conditionGroup2')
  summaryDt <- summaryCond$datatable
  # fcHurdle <- merge(
  #   summaryDt[contrast=='conditionGroup2' & component=='H', .(primerid, `Pr(>Chisq)`)], #hurdle P values
  #   summaryDt[contrast=='conditionGroup2' & component=='logFC', .(primerid, coef, ci.hi, ci.lo)], by='primerid'
  # ) #logFC coefficients
  # fcHurdle[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
  p_val <- summaryDt[summaryDt[, "component"] == "H", 4]
  genes.return <- summaryDt[summaryDt[, "component"] == "H", 1]
  # p_val <- subset(summaryDt, component == "H")[, 4]
  # genes.return <- subset(summaryDt, component == "H")[, 1]
  to.return <- data.frame(p_val, row.names = genes.return)
  return(to.return)
}

# compare two negative binomial regression models
# model one uses only common factors (com.fac)
# model two additionally uses group factor (grp.fac)
#
#' @importFrom stats glm anova coef
#
NBModelComparison <- function(y, theta, latent.data, com.fac, grp.fac) {
  tab <- as.matrix(x = table(y > 0, latent.data[, grp.fac]))
  freqs <- tab['TRUE', ] / apply(X = tab, MARGIN = 2, FUN = sum)
  fit2 <- 0
  fit4 <- 0
  try(
    expr = fit2 <- glm(
      formula = y ~ .,
      data = latent.data[, com.fac, drop = FALSE],
      family = MASS::negative.binomial(theta = theta)
    ),
    silent=TRUE
  )
  try(
    fit4 <- glm(
      formula = y ~ .,
      data = latent.data[, c(com.fac, grp.fac)],
      family = MASS::negative.binomial(theta = theta)
    ),
    silent = TRUE
  )
  if (is.numeric(x = fit2) || is.numeric(x = fit4)) {
    message('One of the glm.nb calls failed')
    return(c(rep(x = NA, 5), freqs))
  }
  pval <- anova(fit2, fit4, test = 'Chisq')$'Pr(>Chi)'[2]
  foi <- 2 + length(x = com.fac)
  log2.fc <- log2(x = 1 / exp(x = coef(object = fit4)[foi]))
  ret <- c(
    fit2$deviance,
    fit4$deviance,
    pval,
    coef(object = fit4)[foi],
    log2.fc,
    freqs
  )
  names(x = ret) <- c(
    'dev1',
    'dev2',
    'pval',
    'coef',
    'log2.fc',
    'freq1',
    'freq2'
  )
  return(ret)
}

# Run a per-feature DE test on a sparse matrix by densifying in row-blocks.
#
# To keep peak memory bounded, the subset of `data.use` is densified at most `block.elems`
# elements at a time (a block of ceiling(block.elems / ncol) features). Normal-
# sized comparisons fit in a single block, so are run without splitting.
RunDEBlocked <- function(data.use, testfun, ..., block.elems = 2e7) {
  # Non-sparse input (e.g. densify = TRUE) already has fast row access.
  if (!inherits(x = data.use, what = "sparseMatrix")) {
    return(testfun(data.use = data.use, ...))
  }
  n <- nrow(x = data.use)
  block <- max(1L, floor(x = block.elems / ncol(x = data.use)))
  if (block >= n) {
    return(testfun(data.use = as.matrix(x = data.use), ...))
  }
  idx <- split(x = seq_len(length.out = n), f = ceiling(x = seq_len(length.out = n) / block))
  do.call(
    what = rbind,
    args = unname(obj = lapply(
      X = idx,
      FUN = function(rows) {
        testfun(data.use = as.matrix(x = data.use[rows, , drop = FALSE]), ...)
      }
    ))
  )
}

PerformDE <- function(
  object,
  cells.1,
  cells.2,
  features,
  test.use,
  verbose,
  min.cells.feature,
  latent.vars,
  densify,
  ...
) {
  if (!(test.use %in% DEmethods_latent()) && !is.null(x = latent.vars)) {
    warning(
      "'latent.vars' is only used for the following tests: ",
      paste(DEmethods_latent(), collapse=", "),
      call. = FALSE,
      immediate. = TRUE
    )
  }
  if (!test.use %in% DEmethods_checkdots()) {
    CheckDots(...)
  }
  data.use <- object[features, c(cells.1, cells.2), drop = FALSE]
  if (densify){
    data.use <- as.matrix(x = data.use)
  }
  de.results <- switch(
    EXPR = test.use,
    'wilcox' = WilcoxDETest(
      data.use = data.use,
      cells.1 = cells.1,
      cells.2 = cells.2,
      verbose = verbose,
      ...
    ),
    'wilcox_limma' = WilcoxDETest(
      data.use = data.use,
      cells.1 = cells.1,
      cells.2 = cells.2,
      verbose = verbose,
      limma = TRUE,
      ...
    ),
    'bimod' = RunDEBlocked(
      data.use = data.use,
      testfun = DiffExpTest,
      cells.1 = cells.1,
      cells.2 = cells.2,
      verbose = verbose
    ),
    'roc' = RunDEBlocked(
      data.use = data.use,
      testfun = MarkerTest,
      cells.1 = cells.1,
      cells.2 = cells.2,
      verbose = verbose
    ),
    't' = RunDEBlocked(
      data.use = data.use,
      testfun = DiffTTest,
      cells.1 = cells.1,
      cells.2 = cells.2,
      verbose = verbose
    ),
    'negbinom' = GLMDETest(
      data.use = data.use,
      cells.1 = cells.1,
      cells.2 = cells.2,
      min.cells = min.cells.feature,
      latent.vars = latent.vars,
      test.use = test.use,
      verbose = verbose
    ),
    'poisson' = GLMDETest(
      data.use = data.use,
      cells.1 = cells.1,
      cells.2 = cells.2,
      min.cells = min.cells.feature,
      latent.vars = latent.vars,
      test.use = test.use,
      verbose = verbose
    ),
    'MAST' = MASTDETest(
      data.use = data.use,
      cells.1 = cells.1,
      cells.2 = cells.2,
      latent.vars = latent.vars,
      verbose = verbose,
      ...
    ),
    "DESeq2" = DESeq2DETest(
      data.use = data.use,
      cells.1 = cells.1,
      cells.2 = cells.2,
      verbose = verbose,
      ...
    ),
    "LR" = LRDETest(
      data.use = data.use,
      cells.1 = cells.1,
      cells.2 = cells.2,
      latent.vars = latent.vars,
      verbose = verbose
    ),
    stop("Unknown test: ", test.use)
  )
  return(de.results)
}

#' Prepare object to run differential expression on SCT assay with multiple models
#'
#' Given a merged object with multiple SCT models, this function uses minimum
#' of the median UMI (calculated using the raw UMI counts) of individual objects
#' to reverse the individual SCT regression model using minimum of median UMI
#' as the sequencing depth covariate.
#' The counts slot of the SCT assay is replaced with recorrected counts and
#' the data slot is replaced with log1p of recorrected counts.
#' @param object Seurat object with SCT assays
#' @param assay Assay name where for SCT objects are stored; Default is 'SCT'
#' @param verbose Print messages and progress
#' @importFrom Matrix Matrix
#' @importFrom SeuratObject SparseEmptyMatrix
#' @importFrom pbapply pblapply
#' @importFrom future.apply future_lapply
#' @importFrom future nbrOfWorkers
#' @importFrom sctransform correct_counts
#' @importFrom SeuratObject JoinLayers
#'
#' @return Returns a Seurat object with recorrected counts and data in the SCT assay.
#' @export
#'
#' @concept differential_expression
#' @template section-progressr
#' @template section-future
#' @examples
#' data("pbmc_small")
#' pbmc_small1 <- SCTransform(object = pbmc_small, variable.features.n = 20, vst.flavor="v1")
#' pbmc_small2 <- SCTransform(object = pbmc_small, variable.features.n = 20, vst.flavor="v1")
#' pbmc_merged <- merge(x = pbmc_small1, y = pbmc_small2)
#' pbmc_merged <- PrepSCTFindMarkers(object = pbmc_merged)
#' markers <- FindMarkers(
#'   object = pbmc_merged,
#'   ident.1 = "0",
#'   ident.2 = "1",
#'   assay = "SCT"
#' )
#' pbmc_subset <- subset(pbmc_merged, idents = c("0", "1"))
#' markers_subset <- FindMarkers(
#'   object = pbmc_subset,
#'   ident.1 = "0",
#'   ident.2 = "1",
#'   assay = "SCT",
#'   recorrect_umi = FALSE
#' )
#'
PrepSCTFindMarkers <- function(object, assay = "SCT", verbose = TRUE) {
  if (verbose && nbrOfWorkers() == 1) {
    my.lapply <- pblapply
  } else {
    my.lapply <- future_lapply
  }
  if (length(x = levels(x = object[[assay]])) == 1) {
    if (verbose) {
      message("Only one SCT model is stored - skipping recalculating corrected counts")
    }
    return(object)
  }
  observed_median_umis <- lapply(
    X = SCTResults(object = object[[assay]], slot = "cell.attributes"),
    FUN = function(x) median(x[, "umi"])
  )
  model.list <- slot(object = object[[assay]], name = "SCTModel.list")
  median_umi.status <- lapply(X = model.list,
                              FUN = function(x) { return(tryCatch(
                                expr = slot(object = x, name = 'median_umi'),
                                error = function(...) {return(NULL)})
                              )})
  if (any(is.null(x = unlist(x = median_umi.status)))){
    # For old SCT objects  median_umi is set to median umi as calculated from obserbed UMIs
    slot(object = object[[assay]], name = "SCTModel.list") <- lapply(X = model.list,
                                                                     FUN = UpdateSlots)
    SCTResults(object = object[[assay]], slot = "median_umi") <- observed_median_umis

  }
  model_median_umis <- SCTResults(object = object[[assay]], slot = "median_umi")
  min_median_umi <- min(unlist(x = observed_median_umis), na.rm = TRUE)
  if (all(unlist(x = model_median_umis) > min_median_umi)){
    if (verbose){
      message("Minimum UMI unchanged. Skipping re-correction.")
    }
    return(object)
  }
  if (verbose) {
    message(paste0("Found ",
                   length(x = levels(x = object[[assay]])),
                   " SCT models.",
                   " Recorrecting SCT counts using minimum median counts: ",
                   min_median_umi))
  }
  umi.assay <- unique(
    x = unlist(
      x = SCTResults(object = object[[assay]], slot = "umi.assay")
    )
  )
  if (length(x = umi.assay) > 1) {
    stop("Multiple UMI assays are used for SCTransform: ",
         paste(umi.assay, collapse = ", ")
    )
  }
  umi.layers <- Layers(object = object, assay = umi.assay, search = 'counts')
  if (length(x = umi.layers) > 1) {
    object[[umi.assay]] <- JoinLayers(
      object = object[[umi.assay]],
      layers = "counts", new = "counts")
  }
  raw_umi <- GetAssayData(object = object, assay = umi.assay, layer = "counts")
  corrected_counts <- Matrix(
    nrow = nrow(x = raw_umi),
    ncol = ncol(x = raw_umi),
    data = 0,
    dimnames = dimnames(x = raw_umi),
    sparse = TRUE
  )
  cell_attr <- SCTResults(object = object[[assay]], slot = "cell.attributes")
  model_pars_fit <- lapply(
    X = SCTResults(object = object[[assay]], slot = "feature.attributes"),
    FUN = function(x) x[, c("theta", "(Intercept)", "log_umi")]
  )
  arguments <- SCTResults(object = object[[assay]], slot = "arguments")
  model_str <- SCTResults(object = object[[assay]], slot = "model")
  set_median_umi <- rep(min_median_umi, length(levels(x = object[[assay]])))
  names(set_median_umi) <- levels(x = object[[assay]])
  set_median_umi <- as.list(set_median_umi)
  all_genes <- rownames(x = object[[assay]])
  # correct counts
  # Recorrect raw UMI counts for one SCT model at a shared target depth
  # We only recorrect genes whose fitted SCT parameters are all finite
  # Genes with non finite fitted parameters can produce invalid values during recorrection,
  # so we exclude them here and add them back later as zero rows to preserve the full feature set
  my.correct_counts <- function(model_name){
    # Retrieve model parameters SCT model per-gene
    pars <- model_pars_fit[[model_name]]

    # Keep only genes with finite fitted parameters
    # This is because in split-layer SCT workflows, some sparse genes may have NaN values for
    # fitted parameters due to low expression levels in some layers
    valid_genes <- rownames(pars)[
      !is.na(pars[, "theta"]) &
      is.finite(pars[, "(Intercept)"]) &
      is.finite(pars[, "log_umi"])
    ]

    cells <- rownames(cell_attr[[model_name]])
    # Prepare input list for correct_counts function
    x <- list(
      model_str = model_str[[model_name]],
      arguments = arguments[[model_name]],
      model_pars_fit = as.matrix(pars[valid_genes, , drop = FALSE]),
      cell_attr = cell_attr[[model_name]]
    )
    # Retrieve raw UMI counts for valid genes and cells
    umi <- raw_umi[valid_genes, cells, drop = FALSE]

    # Recorrect counts to the common minimum median UMI depth so that counts from different 
    #SCT models are on the same scale
    umi_corrected <- correct_counts(
      x = x,
      umi = umi,
      verbosity = 0,
      scale_factor = min_median_umi
    )

    # Remove invalid rows to prevent downstream issues with non-finite values in corrected counts
    # Want to avoid any NaN propagation into SCT assay
    bad_idx <- which(!is.finite(umi_corrected@x))
    # correct_counts produces a sparse matrix so we can use the x and i slots for efficiency
    if (length(bad_idx) > 0) {
      # For sparse matrices, @x stores nonzero values and @i stores their zero based row indices
      bad_rows <- unique(umi_corrected@i[bad_idx] + 1) 
      umi_corrected <- umi_corrected[-bad_rows, , drop = FALSE]
    }
    # Add back any genes that were not recorrected for this model (removed above)
    # This keeps the corrected matrix aligned to the full assay wide gene set
    missing_features <- setdiff(x = all_genes, y = rownames(x = umi_corrected))
    gc(verbose = FALSE)
    empty <- SparseEmptyMatrix(nrow = length(x = missing_features), ncol = ncol(x = umi_corrected))
    rownames(x = empty) <- missing_features
    colnames(x = empty) <- colnames(x = umi_corrected)

    # Restore full gene set, original row order
    umi_corrected <- rbind(umi_corrected, empty)[all_genes,]

    return(umi_corrected)
  }
  corrected_counts.list <- my.lapply(X = levels(x = object[[assay]]),
                                     FUN = my.correct_counts)
  names(x = corrected_counts.list) <- levels(x = object[[assay]])

  corrected_counts <- do.call(what = MergeSparseMatrices, args = corrected_counts.list)
  corrected_counts <- as.sparse(x = corrected_counts)
  corrected_data <- log1p(x = corrected_counts)
  suppressWarnings({object <- SetAssayData(object = object,
                                           assay = assay,
                                           layer = "counts",
                                           new.data = corrected_counts)})
  suppressWarnings({object <- SetAssayData(object = object,
                                           assay = assay,
                                           layer = "data",
                                           new.data = corrected_data)})
  SCTResults(object = object[[assay]], slot = "median_umi") <- set_median_umi
  return(object)
}

PrepSCTFindMarkers.V5 <- function(object, assay = "SCT", umi.assay = "RNA", layer = "counts", verbose = TRUE) {
  layers <- Layers(object = object[[umi.assay]], search = layer)
  dataset.names <- gsub(pattern = paste0(layer, "."), replacement = "", x = layers)
  for (i in seq_along(along.with = layers)) {
    l <- layers[i]
    counts <- LayerData(
      object = object[[umi.assay]],
      layer = l
      )
  }
  cells.grid <- DelayedArray::colAutoGrid(x = counts, ncol = min(length(Cells(object)), ncol(counts)))
}

# given a UMI count matrix, estimate NB theta parameter for each gene
# and use fit of relationship with mean to assign regularized theta to each gene
#
#' @importFrom stats glm loess poisson
#' @importFrom utils txtProgressBar setTxtProgressBar
#
RegularizedTheta <- function(cm, latent.data, min.theta = 0.01, bin.size = 128) {
  genes.regress <- rownames(x = cm)
  bin.ind <- ceiling(x = 1:length(x = genes.regress) / bin.size)
  max.bin <- max(bin.ind)
  message('Running Poisson regression (to get initial mean), and theta estimation per gene')
  pb <- txtProgressBar(min = 0, max = max.bin, style = 3, file = stderr())
  theta.estimate <- c()
  for (i in 1:max.bin) {
    genes.bin.regress <- genes.regress[bin.ind == i]
    bin.theta.estimate <- unlist(
      x = parallel::mclapply(
        X = genes.bin.regress,
        FUN = function(j) {
          return(as.numeric(x = MASS::theta.ml(
            y = cm[j, ],
            mu = glm(
              formula = cm[j, ] ~ .,
              data = latent.data,
              family = poisson
            )$fitted
          )))
        }
      ),
      use.names = FALSE
    )
    theta.estimate <- c(theta.estimate, bin.theta.estimate)
    setTxtProgressBar(pb = pb, value = i)
  }
  close(con = pb)
  UMI.mean <- apply(X = cm, MARGIN = 1, FUN = mean)
  var.estimate <- UMI.mean + (UMI.mean ^ 2) / theta.estimate
  for (span in c(1/3, 1/2, 3/4, 1)) {
    fit <- loess(
      formula = log10(x = var.estimate) ~ log10(x = UMI.mean),
      span = span
    )
    if (! any(is.na(x = fit$fitted))) {
      message(sprintf(
        'Used loess with span %1.2f to fit mean-variance relationship\n',
        span
      ))
      break
    }
  }
  if (any(is.na(x = fit$fitted))) {
    stop('Problem when fitting NB gene variance in RegularizedTheta - NA values were fitted.')
  }
  theta.fit <- (UMI.mean ^ 2) / ((10 ^ fit$fitted) - UMI.mean)
  names(x = theta.fit) <- genes.regress
  to.fix <- theta.fit <= min.theta | is.infinite(x = theta.fit)
  if (any(to.fix)) {
    message(
      'Fitted theta below ',
      min.theta,
      ' for ',
      sum(to.fix),
      ' genes, setting them to ',
      min.theta
    )
    theta.fit[to.fix] <- min.theta
  }
  return(theta.fit)
}

# FindMarkers helper function for cell grouping error checking
ValidateCellGroups <- function(
  object,
  cells.1,
  cells.2,
  min.cells.group
) {
  if (length(x = cells.1) == 0) {
    stop("Cell group 1 is empty - no cells with identity class ", cells.1)
  } else if (length(x = cells.2) == 0) {
    stop("Cell group 2 is empty - no cells with identity class ", cells.2)
    return(NULL)
  } else if (length(x = cells.1) < min.cells.group) {
    stop("Cell group 1 has fewer than ", min.cells.group, " cells")
  } else if (length(x = cells.2) < min.cells.group) {
    stop("Cell group 2 has fewer than ", min.cells.group, " cells")
  } else if (any(!cells.1 %in% colnames(x = object))) {
    bad.cells <- colnames(x = object)[which(x = !as.character(x = cells.1) %in% colnames(x = object))]
    stop(
      "The following cell names provided to cells.1 are not present: ",
      paste(bad.cells, collapse = ", ")
    )
  } else if (any(!cells.2 %in% colnames(x = object))) {
    bad.cells <- colnames(x = object)[which(x = !as.character(x = cells.2) %in% colnames(x = object))]
    stop(
      "The following cell names provided to cells.2 are not present: ",
      paste(bad.cells, collapse = ", ")
    )
  }
}

# Differential expression using Wilcoxon Rank Sum
#
# Identifies differentially expressed genes between two groups of cells using
# a Wilcoxon Rank Sum test. Makes use of presto::wilcoxauc for a more efficient
# implementation of the wilcoxon test. If presto is not installed, or if limma
# is requested, makes use of limma::rankSumTestWithCorrelation for a
# more efficient implementation of the wilcoxon test. Thanks to Yunshun Chen and
# Gordon Smyth for suggesting the limma implementation. If limma is also not installed,
# uses wilcox.test.
#
# @param data.use Data matrix to test
# @param cells.1 Group 1 cells
# @param cells.2 Group 2 cells
# @param verbose Print a progress bar
# @param limma If limma should be used for testing; default is FALSE
# @param ... Extra parameters passed to wilcox.test
#
# @return Returns a p-value ranked matrix of putative differentially expressed
# features
#
#' @importFrom pbapply pbsapply
#' @importFrom stats wilcox.test
#' @importFrom future.apply future_sapply
#' @importFrom future nbrOfWorkers
#
# @export
#
# @examples
# data("pbmc_small")
# pbmc_small
# WilcoxDETest(pbmc_small, cells.1 = WhichCells(object = pbmc_small, idents = 1),
#             cells.2 = WhichCells(object = pbmc_small, idents = 2))
#
WilcoxDETest <- function(
  data.use,
  cells.1,
  cells.2,
  verbose = TRUE,
  limma = FALSE,
  ...
) {
  j <- seq_len(length.out = length(x = cells.1))
  my.sapply <- ifelse(
    test = verbose && nbrOfWorkers() == 1,
    yes = pbsapply,
    no = future_sapply
  )
  overflow.check <- ifelse(
    test = is.na(x = suppressWarnings(length(x = data.use[1, ]) * length(x = data.use[1, ]))),
    yes = FALSE,
    no = TRUE
  )
  presto.check <- requireNamespace("presto", quietly = TRUE)
  limma.check <- requireNamespace("limma", quietly = TRUE)
  group.info <- data.frame(row.names = c(cells.1, cells.2))
  group.info[cells.1, "group"] <- "Group1"
  group.info[cells.2, "group"] <- "Group2"
  group.info[, "group"] <- factor(x = group.info[, "group"])
  if (presto.check[1] && (!limma)) {
    res <- presto::wilcoxauc(X = data.use, y = group.info[, "group"])
    res <- res[1:(nrow(x = res)/2),]
    p_val <- res$pval
  } else {
    if (getOption('Seurat.presto.wilcox.msg', TRUE) && (!limma)) {
      message(
        "For a (much!) faster implementation of the Wilcoxon Rank Sum Test,",
        "\n(default method for FindMarkers) please install the presto package",
        "\n--------------------------------------------",
        "\ninstall.packages('devtools')",
        "\ndevtools::install_github('immunogenomics/presto')",
        "\n--------------------------------------------",
        "\nAfter installation of presto, Seurat will automatically use the more ",
        "\nefficient implementation (no further action necessary).",
        "\nThis message will be shown once per session"
      )
      options(Seurat.presto.wilcox.msg = FALSE)
    }
    # These fallback paths extract one gene row per iteration; convert the
    # column-major subset to row-major once so per-gene access is efficient.
    if (inherits(x = data.use, what = "CsparseMatrix")) {
      data.use <- as(object = data.use, Class = "RsparseMatrix")
    }
    if (limma.check[1] && overflow.check) {
      p_val <- my.sapply(
        X = 1:nrow(x = data.use),
        FUN = function(x) {
          return(min(2 * min(limma::rankSumTestWithCorrelation(index = j, statistics = data.use[x, ])), 1))
        }
      )
    } else {
      if (limma && overflow.check) {
        stop(
          "To use the limma implementation of the Wilcoxon Rank Sum Test,
        please install the limma package:
        --------------------------------------------
        install.packages('BiocManager')
        BiocManager::install('limma')
        --------------------------------------------"
        )
      } else {
        # Use the vector interface of wilcox.test rather than the formula
        # interface, which re-parses a model.frame for every gene. Columns of
        # data.use are ordered c(cells.1, cells.2), so `j` selects group 1 and
        # the remaining columns are group 2 (matching the factor level order
        # used by the formula method: Group1 -> x, Group2 -> y).
        p_val <- my.sapply(
          X = 1:nrow(x = data.use),
          FUN = function(x) {
            vals <- data.use[x, ]
            return(wilcox.test(x = vals[j], y = vals[-j], ...)$p.value)
          }
        )
      }
    }
  }
  return(data.frame(p_val, row.names = rownames(x = data.use)))
}
