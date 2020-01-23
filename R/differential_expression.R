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
#' @importFrom ape drop.tip
#' @importFrom stats setNames
#'
#' @export
#'
#' @aliases FindAllMarkersNode
#'
#' @examples
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
  logfc.threshold = 0.25,
  test.use = 'wilcox',
  slot = 'data',
  min.pct = 0.1,
  min.diff.pct = -Inf,
  node = NULL,
  verbose = TRUE,
  only.pos = FALSE,
  max.cells.per.ident = Inf,
  random.seed = 1,
  latent.vars = NULL,
  min.cells.feature = 3,
  min.cells.group = 3,
  pseudocount.use = 1,
  return.thresh = 1e-2,
  ...
) {
  MapVals <- function(vec, from, to) {
    vec2 <- setNames(object = to, nm = from)[as.character(x = vec)]
    vec2[is.na(x = vec2)] <- vec[is.na(x = vec2)]
    return(unname(obj = vec2))
  }
  if ((test.use == "roc") && (return.thresh == 1e-2)) {
    return.thresh <- 0.7
  }
  if (is.null(x = node)) {
    idents.all <- sort(x = unique(x = Idents(object = object)))
  } else {
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
    tree <- drop.tip(phy = tree, tip = drop.children)
    new.nodes <- unique(x = tree$edge[, 1, drop = TRUE])
    idents.all <- (tree$Nnode + 2):max(tree$edge)
  }
  genes.de <- list()
  messages <- list()
  for (i in 1:length(x = idents.all)) {
    if (verbose) {
      message("Calculating cluster ", idents.all[i])
    }
    genes.de[[i]] <- tryCatch(
      expr = {
        FindMarkers(
          object = object,
          assay = assay,
          ident.1 = if (is.null(x = node)) {
            idents.all[i]
          } else {
            tree
          },
          ident.2 = if (is.null(x = node)) {
            NULL
          } else {
            idents.all[i]
          },
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
          pseudocount.use = pseudocount.use,
          ...
        )
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
  gde.all <- data.frame()
  for (i in 1:length(x = idents.all)) {
    if (is.null(x = unlist(x = genes.de[i]))) {
      next
    }
    gde <- genes.de[[i]]
    if (nrow(x = gde) > 0) {
      if (test.use == "roc") {
        gde <- subset(
          x = gde,
          subset = (myAUC > return.thresh | myAUC < (1 - return.thresh))
        )
      } else if (is.null(x = node) || test.use %in% c('bimod', 't')) {
        gde <- gde[order(gde$p_val, -gde[, 2]), ]
        gde <- subset(x = gde, subset = p_val < return.thresh)
      }
      if (nrow(x = gde) > 0) {
        gde$cluster <- idents.all[i]
        gde$gene <- rownames(x = gde)
      }
      if (nrow(x = gde) > 0) {
        gde.all <- rbind(gde.all, gde)
      }
    }
  }
  if ((only.pos) && nrow(x = gde.all) > 0) {
    return(subset(x = gde.all, subset = gde.all[, 2] > 0))
  }
  rownames(x = gde.all) <- make.unique(names = as.character(x = gde.all$gene))
  if (nrow(x = gde.all) == 0) {
    warning("No DE genes identified", call. = FALSE, immediate. = TRUE)
  }
  if (length(x = messages) > 0) {
    warning("The following tests were not performed: ", call. = FALSE, immediate. = TRUE)
    for (i in 1:length(x = messages)) {
      if (!is.null(x = messages[[i]])) {
        warning("When testing ", idents.all[i], " versus all:\n\t", messages[[i]], call. = FALSE, immediate. = TRUE)
      }
    }
  }
  if (!is.null(x = node)) {
    gde.all$cluster <- MapVals(
      vec = gde.all$cluster,
      from = new.nodes,
      to = orig.nodes
    )
  }
  return(gde.all)
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
#' @importFrom metap minimump
#'
#' @export
#'
#' @examples
#' \dontrun{
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
  meta.method = minimump,
  verbose = TRUE,
  ...
) {
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
        return(meta.method(x)$p)
      }
    ))
    colnames(x = combined.pval) <- paste0(
      as.character(x = formals()$meta.method),
      "_p_val"
    )
    markers.combined <- cbind(markers.combined, combined.pval)
    markers.combined <- markers.combined[order(markers.combined[, paste0(as.character(x = formals()$meta.method), "_p_val")]), ]
  } else {
    warning("Only a single group was tested", call. = FALSE, immediate. = TRUE)
  }
  return(markers.combined)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Methods for Seurat-defined generics
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @param cells.1 Vector of cell names belonging to group 1
#' @param cells.2 Vector of cell names belonging to group 2
#' @param counts Count matrix if using scale.data for DE tests. This is used for
#' computing pct.1 and pct.2 and for filtering features based on fraction
#' expressing
#' @param features Genes to test. Default is to use all genes
#' @param logfc.threshold Limit testing to genes which show, on average, at least
#' X-fold difference (log-scale) between the two groups of cells. Default is 0.25
#' Increasing logfc.threshold speeds up the function, but can miss weaker signals.
#' @param test.use Denotes which test to use. Available options are:
#' \itemize{
#'  \item{"wilcox"} : Identifies differentially expressed genes between two
#'  groups of cells using a Wilcoxon Rank Sum test (default)
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
#' by not testing genes that are very infrequently expressed. Default is 0.1
#' @param min.diff.pct  only test genes that show a minimum difference in the
#' fraction of detection between the two groups. Set to -Inf by default
#' @param only.pos Only return positive markers (FALSE by default)
#' @param verbose Print a progress bar once expression testing begins
#' @param max.cells.per.ident Down sample each identity class to a max number.
#' Default is no downsampling. Not activated by default (set to Inf)
#' @param random.seed Random seed for downsampling
#' @param latent.vars Variables to test, used only when \code{test.use} is one of
#' 'LR', 'negbinom', 'poisson', or 'MAST'
#' @param min.cells.feature Minimum number of cells expressing the feature in at least one
#' of the two groups, currently only used for poisson and negative binomial tests
#' @param min.cells.group Minimum number of cells in one of the groups
#' @param pseudocount.use Pseudocount to add to averaged expression values when
#' calculating logFC. 1 by default.
#'
#' @importFrom Matrix rowSums rowMeans
#' @importFrom stats p.adjust
#'
#' @rdname FindMarkers
#' @export
#' @method FindMarkers default
#'
FindMarkers.default <- function(
  object,
  slot = "data",
  counts = numeric(),
  cells.1 = NULL,
  cells.2 = NULL,
  features = NULL,
  reduction = NULL,
  logfc.threshold = 0.25,
  test.use = 'wilcox',
  min.pct = 0.1,
  min.diff.pct = -Inf,
  verbose = TRUE,
  only.pos = FALSE,
  max.cells.per.ident = Inf,
  random.seed = 1,
  latent.vars = NULL,
  min.cells.feature = 3,
  min.cells.group = 3,
  pseudocount.use = 1,
  ...
) {
  features <- features %||% rownames(x = object)
  methods.noprefiliter <- c("DESeq2")
  if (test.use %in% methods.noprefiliter) {
    features <- rownames(x = object)
    min.diff.pct <- -Inf
    logfc.threshold <- 0
  }
  # error checking
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
  # feature selection (based on percentages)
  data <- switch(
    EXPR = slot,
    'scale.data' = counts,
    object
  )
  if (is.null(x = reduction)) {
    thresh.min <- 0
    pct.1 <- round(
      x = rowSums(x = data[features, cells.1, drop = FALSE] > thresh.min) /
        length(x = cells.1),
      digits = 3
    )
    pct.2 <- round(
      x = rowSums(x = data[features, cells.2, drop = FALSE] > thresh.min) /
        length(x = cells.2),
      digits = 3
    )
    data.alpha <- cbind(pct.1, pct.2)
    colnames(x = data.alpha) <- c("pct.1", "pct.2")
    alpha.min <- apply(X = data.alpha, MARGIN = 1, FUN = max)
    names(x = alpha.min) <- rownames(x = data.alpha)
    features <- names(x = which(x = alpha.min > min.pct))
    if (length(x = features) == 0) {
      stop("No features pass min.pct threshold")
    }
    alpha.diff <- alpha.min - apply(X = data.alpha, MARGIN = 1, FUN = min)
    features <- names(
      x = which(x = alpha.min > min.pct & alpha.diff > min.diff.pct)
    )
    if (length(x = features) == 0) {
      stop("No features pass min.diff.pct threshold")
    }
  } else {
    data.alpha <- data.frame(
      pct.1 = rep(x = NA, times = length(x = features)),
      pct.2 = rep(x = NA, times = length(x = features))
    )
  }
  # feature selection (based on average difference)
  mean.fxn <- if (is.null(x = reduction) && slot != "scale.data") {
    switch(
      EXPR = slot,
      'data' = function(x) {
        return(log(x = rowMeans(x = expm1(x = x)) + pseudocount.use))
      },
      function(x) {
        return(log(x = rowMeans(x = x) + pseudocount.use))
      }
    )
  } else {
    rowMeans
  }
  data.1 <- mean.fxn(data[features, cells.1, drop = FALSE])
  data.2 <- mean.fxn(data[features, cells.2, drop = FALSE])

  total.diff <- (data.1 - data.2)
  if (is.null(x = reduction) && slot != "scale.data") {
    features.diff <- if (only.pos) {
      names(x = which(x = total.diff > logfc.threshold))
    } else {
      names(x = which(x = abs(x = total.diff) > logfc.threshold))
    }
    features <- intersect(x = features, y = features.diff)
    if (length(x = features) == 0) {
      stop("No features pass logfc.threshold threshold")
    }
  }
  if (max.cells.per.ident < Inf) {
    set.seed(seed = random.seed)
    # Should be cells.1 and cells.2?
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
  # perform DE
  if (!(test.use %in% c('negbinom', 'poisson', 'MAST', "LR")) && !is.null(x = latent.vars)) {
    warning(
      "'latent.vars' is only used for 'negbinom', 'poisson', 'LR', and 'MAST' tests",
      call. = FALSE,
      immediate. = TRUE
    )
  }
  if (!test.use %in% c('wilcox', 'MAST', 'DESeq2')) {
    CheckDots(...)
  }
  de.results <- switch(
    EXPR = test.use,
    'wilcox' = WilcoxDETest(
      data.use = object[features, c(cells.1, cells.2), drop = FALSE],
      cells.1 = cells.1,
      cells.2 = cells.2,
      verbose = verbose,
      ...
    ),
    'bimod' = DiffExpTest(
      data.use = object[features, c(cells.1, cells.2), drop = FALSE],
      cells.1 = cells.1,
      cells.2 = cells.2,
      verbose = verbose
    ),
    'roc' = MarkerTest(
      data.use = object[features, c(cells.1, cells.2), drop = FALSE],
      cells.1 = cells.1,
      cells.2 = cells.2,
      verbose = verbose
    ),
    't' = DiffTTest(
      data.use = object[features, c(cells.1, cells.2), drop = FALSE],
      cells.1 = cells.1,
      cells.2 = cells.2,
      verbose = verbose
    ),
    'negbinom' = GLMDETest(
      data.use = object[features, c(cells.1, cells.2), drop = FALSE],
      cells.1 = cells.1,
      cells.2 = cells.2,
      min.cells = min.cells.feature,
      latent.vars = latent.vars,
      test.use = test.use,
      verbose = verbose
    ),
    'poisson' = GLMDETest(
      data.use = object[features, c(cells.1, cells.2), drop = FALSE],
      cells.1 = cells.1,
      cells.2 = cells.2,
      min.cells = min.cells.feature,
      latent.vars = latent.vars,
      test.use = test.use,
      verbose = verbose
    ),
    'MAST' = MASTDETest(
      data.use = object[features, c(cells.1, cells.2), drop = FALSE],
      cells.1 = cells.1,
      cells.2 = cells.2,
      latent.vars = latent.vars,
      verbose = verbose,
      ...
    ),
    "DESeq2" = DESeq2DETest(
      data.use = object[features, c(cells.1, cells.2), drop = FALSE],
      cells.1 = cells.1,
      cells.2 = cells.2,
      verbose = verbose,
      ...
    ),
    "LR" = LRDETest(
      data.use = object[features, c(cells.1, cells.2), drop = FALSE],
      cells.1 = cells.1,
      cells.2 = cells.2,
      latent.vars = latent.vars,
      verbose = verbose
    ),
    stop("Unknown test: ", test.use)
  )
  if (is.null(x = reduction)) {
    diff.col <- ifelse(
      test = slot == "scale.data" || test.use == 'roc',
      yes = "avg_diff",
      no = "avg_logFC"
    )
    de.results[, diff.col] <- total.diff[rownames(x = de.results)]
    de.results <- cbind(de.results, data.alpha[rownames(x = de.results), , drop = FALSE])
  } else {
    diff.col <- "avg_diff"
    de.results[, diff.col] <- total.diff[rownames(x = de.results)]
  }
  if (only.pos) {
    de.results <- de.results[de.results[, diff.col] > 0, , drop = FALSE]
  }
  if (test.use == "roc") {
    de.results <- de.results[order(-de.results$power, -de.results[, diff.col]), ]
  } else {
    de.results <- de.results[order(de.results$p_val, -de.results[, diff.col]), ]
    de.results$p_val_adj = p.adjust(
      p = de.results$p_val,
      method = "bonferroni",
      n = nrow(x = object)
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
#' @param reduction Reduction to use in differential expression testing - will test for DE on cell embeddings
#' @param group.by Regroup cells into a different identity class prior to performing differential expression (see example)
#' @param subset.ident Subset a particular identity class prior to regrouping. Only relevant if group.by is set (see example)
#' @param assay Assay to use in differential expression testing
#' @param slot Slot to pull data from; note that if \code{test.use} is "negbinom", "poisson", or "DESeq2",
#' \code{slot} will be set to "counts"
#'
#' @importFrom methods is
#'
#' @rdname FindMarkers
#' @export
#' @method FindMarkers Seurat
#'
FindMarkers.Seurat <- function(
  object,
  ident.1 = NULL,
  ident.2 = NULL,
  group.by = NULL,
  subset.ident = NULL,
  assay = NULL,
  slot = 'data',
  reduction = NULL,
  features = NULL,
  logfc.threshold = 0.25,
  test.use = "wilcox",
  min.pct = 0.1,
  min.diff.pct = -Inf,
  verbose = TRUE,
  only.pos = FALSE,
  max.cells.per.ident = Inf,
  random.seed = 1,
  latent.vars = NULL,
  min.cells.feature = 3,
  min.cells.group = 3,
  pseudocount.use = 1,
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
  data.slot <- ifelse(
    test = test.use %in% c("negbinom", "poisson", "DESeq2"),
    yes = 'counts',
    no = slot
  )
  if (is.null(x = reduction)) {
    assay <- assay %||% DefaultAssay(object = object)
    data.use <-  GetAssayData(object = object[[assay]], slot = data.slot)
  } else {
    if (data.slot == "counts") {
      stop("The following tests cannot be used when specifying a reduction as they assume a count model: negbinom, poisson, DESeq2")
    }
    data.use <- t(x = Embeddings(object = object, reduction = reduction))
  }
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
      any(as.character(x = ident.1) %in% colnames(x = data.use))) {
    bad.cells <- colnames(x = data.use)[which(x = !as.character(x = ident.1) %in% colnames(x = data.use))]
    if (length(x = bad.cells) > 0) {
      stop(paste0("The following cell names provided to ident.1 are not present in the object: ", paste(bad.cells, collapse = ", ")))
    }
  } else {
    ident.1 <- WhichCells(object = object, idents = ident.1)
  }
  # if NULL for ident.2, use all other cells
  if (length(x = as.vector(x = ident.2)) > 1 &&
      any(as.character(x = ident.2) %in% colnames(x = data.use))) {
    bad.cells <- colnames(x = data.use)[which(!as.character(x = ident.2) %in% colnames(x = data.use))]
    if (length(x = bad.cells) > 0) {
      stop(paste0("The following cell names provided to ident.2 are not present in the object: ", paste(bad.cells, collapse = ", ")))
    }
  } else {
    if (is.null(x = ident.2)) {
      ident.2 <- setdiff(x = colnames(x = data.use), y = ident.1)
    } else {
      ident.2 <- WhichCells(object = object, idents = ident.2)
    }
  }
  if (!is.null(x = latent.vars)) {
    latent.vars <- FetchData(
      object = object,
      vars = latent.vars,
      cells = c(ident.1, ident.2)
    )
  }
  counts <- switch(
    EXPR = data.slot,
    'scale.data' = GetAssayData(object = object[[assay]], slot = "counts"),
    numeric()
  )
  de.results <- FindMarkers(
    object = data.use,
    slot = data.slot,
    counts = counts,
    cells.1 = ident.1,
    cells.2 = ident.2,
    features = features,
    reduction = reduction,
    logfc.threshold = logfc.threshold,
    test.use = test.use,
    min.pct = min.pct,
    min.diff.pct = min.diff.pct,
    verbose = verbose,
    only.pos = only.pos,
    max.cells.per.ident = max.cells.per.ident,
    random.seed = random.seed,
    latent.vars = latent.vars,
    min.cells.feature = min.cells.feature,
    min.cells.group = min.cells.group,
    pseudocount.use = pseudocount.use,
    ...
  )
  return(de.results)
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
  if (!PackageCheck('DESeq2', error = FALSE)) {
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
# @param latent.vars Confounding variables to adjust for in DE test. Default is
# "nUMI", which adjusts for cellular depth (i.e. cellular detection rate). For
# non-UMI based data, set to nGene instead.
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
#
# @export
#
# @examples
# \dontrun{
#   pbmc_small
#   MASTDETest(pbmc_small, cells.1 = WhichCells(object = pbmc_small, idents = 1),
#               cells.2 = WhichCells(object = pbmc_small, idents = 2))
# }
#
MASTDETest <- function(
  data.use,
  cells.1,
  cells.2,
  latent.vars = NULL,
  verbose = TRUE,
  ...
) {
  # Check for MAST
  if (!PackageCheck('MAST', error = FALSE)) {
    stop("Please install MAST - learn more at https://github.com/RGLab/MAST")
  }
  if (length(x = latent.vars) > 0) {
    latent.vars <- scale(x = latent.vars)
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
  summaryCond <- summary(object = zlmCond, doLRT = 'conditionGroup2')
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

# Differential expression using Wilcoxon Rank Sum
#
# Identifies differentially expressed genes between two groups of cells using
# a Wilcoxon Rank Sum test
#
# @param data.use Data matrix to test
# @param cells.1 Group 1 cells
# @param cells.2 Group 2 cells
# @param verbose Print a progress bar
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
# pbmc_small
# WilcoxDETest(pbmc_small, cells.1 = WhichCells(object = pbmc_small, idents = 1),
#             cells.2 = WhichCells(object = pbmc_small, idents = 2))
#
WilcoxDETest <- function(
  data.use,
  cells.1,
  cells.2,
  verbose = TRUE,
  ...
) {
  group.info <- data.frame(row.names = c(cells.1, cells.2))
  group.info[cells.1, "group"] <- "Group1"
  group.info[cells.2, "group"] <- "Group2"
  group.info[, "group"] <- factor(x = group.info[, "group"])
  data.use <- data.use[, rownames(x = group.info), drop = FALSE]
  my.sapply <- ifelse(
    test = verbose && nbrOfWorkers() == 1,
    yes = pbsapply,
    no = future_sapply
  )
  p_val <- my.sapply(
    X = 1:nrow(x = data.use),
    FUN = function(x) {
      return(wilcox.test(data.use[x, ] ~ group.info[, "group"], ...)$p.value)
    }
  )
  return(data.frame(p_val, row.names = rownames(x = data.use)))
}
