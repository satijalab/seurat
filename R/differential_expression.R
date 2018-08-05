#' @param cells.1 Vector of cell names belonging to group 1
#' @param cells.2 Vector of cell names belonging to group 2
#'
#' @describeIn FindMarkers Run differential expression test on matrix
#' @export
#' @method FindMarkers default
#'
FindMarkers.default <- function(
  object,
  cells.1 = NULL,
  cells.2 = NULL,
  features.use = NULL,
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
  features.use <- features.use %||% rownames(object)
  methods.noprefiliter <- c("DESeq2", "zingeR")
  if (test.use %in% methods.noprefiliter) {
    features.use <- rownames(object)
    min.diff.pct <- -Inf
    logfc.threshold <- 0
  }
  # error checking
  if (length(x = cells.1) == 0) {
    message(paste("Cell group 1 is empty - no cells with identity class", cells.1))
    return(NULL)
  }
  if (length(x = cells.2) == 0) {
    message(paste("Cell group 2 is empty - no cells with identity class", cells.2))
    return(NULL)
  }
  if (length(cells.1) < min.cells.group) {
    stop(paste("Cell group 1 has fewer than", as.character(min.cells.group), "cells"))
  }
  if (length(cells.2) < min.cells.group) {
    stop(paste("Cell group 2 has fewer than", as.character(min.cells.group), " cells"))
  }
  if(any(!cells.1 %in% colnames(object))) {
    bad.cells <- colnames(object)[which(!as.character(x = cells.1) %in% colnames(object))]
    stop(paste0("The following cell names provided to cells.1 are not present: ", paste(bad.cells, collapse = ", ")))
  }
  if(any(!cells.2 %in% colnames(object))) {
    bad.cells <- colnames(object)[which(!as.character(x = cells.2) %in% colnames(object))]
    stop(paste0("The following cell names provided to cells.2 are not present: ", paste(bad.cells, collapse = ", ")))
  }

  # feature selection (based on percentages)
  thresh.min <- 0
  pct.1 <- round(
    x = apply(
      X = object[features.use, cells.1, drop = F],
      MARGIN = 1,
      FUN = function(x) {
        return(sum(x > thresh.min) / length(x = x))
      }
    ),
    digits = 3
  )
  pct.2 <- round(
    x = apply(
      X = object[features.use, cells.2, drop = F],
      MARGIN = 1,
      FUN = function(x) {
        return(sum(x > thresh.min) / length(x = x))
      }
    ),
    digits = 3
  )
  data.alpha <- cbind(pct.1, pct.2)
  colnames(x = data.alpha) <- c("pct.1", "pct.2")
  alpha.min <- apply(X = data.alpha, MARGIN = 1, FUN = max)
  names(x = alpha.min) <- rownames(x = data.alpha)
  features.use <- names(x = which(x = alpha.min > min.pct))
  if (length(x = features.use) == 0) {
    stop("No features pass min.pct threshold")
  }
  alpha.diff <- alpha.min - apply(X = data.alpha, MARGIN = 1, FUN = min)
  features.use <- names(
    x = which(x = alpha.min > min.pct & alpha.diff > min.diff.pct)
  )
  if (length(x = features.use) == 0) {
    stop("No features pass min.diff.pct threshold")
  }
  # gene selection (based on average difference)
  data.1 <- apply(X = object[features.use, cells.1, drop = F],
                  MARGIN = 1,
                  FUN = function(x) log(x = mean(x = expm1(x = x)) + pseudocount.use))
  data.2 <- apply(X = object[features.use, cells.2, drop = F],
                  MARGIN = 1,
                  FUN = function(x) log(x = mean(x = expm1(x = x)) + pseudocount.use))
  total.diff <- (data.1 - data.2)
  if (!only.pos) features.diff <- names(x = which(x = abs(x = total.diff) > logfc.threshold))
  if (only.pos) features.diff <- names(x = which(x = total.diff > logfc.threshold))
  features.use <- intersect(x = features.use, y = features.diff)
  if (length(x = features.use) == 0) {
    stop("No features pass logfc.threshold threshold")
  }
  if (max.cells.per.ident < Inf) {
    set.seed(seed = random.seed)
    if (length(ident.1) > max.cells.per.ident) {
      cells.1 <- sample(x = cells.1, size = max.cells.per.ident)
    }
    if (length(ident.2) > max.cells.per.ident) {
      cells.2 = sample(x = cells.2, size = max.cells.per.ident)
    }
  }
  # perform DE
  if (!(test.use %in% c('negbinom', 'poisson', 'MAST')) && !is.null(x = latent.vars)) {
    warning("'latent.vars' is only used for 'negbinom', 'poisson', and 'MAST' tests")
  }
  de.results <- switch(
    EXPR = test.use,
    'wilcox' = WilcoxDETest(
      data.use = object[features.use, c(cells.1, cells.2)],
      cells.1 = cells.1,
      cells.2 = cells.2,
      verbose = verbose,
      ...
    ),
    'bimod' = DiffExpTest(
      data.use = object[features.use, c(cells.1, cells.2)],
      cells.1 = cells.1,
      cells.2 = cells.2,
      verbose = verbose
    ),
    'roc' = MarkerTest(
      data.use = object[features.use, c(cells.1, cells.2)],
      cells.1 = cells.1,
      cells.2 = cells.2,
      verbose = verbose
    ),
    't' = DiffTTest(
      data.use = object[features.use, c(cells.1, cells.2)],
      cells.1 = cells.1,
      cells.2 = cells.2,
      verbose = verbose
    ),
    'tobit' = TobitTest(
      data.use = object[features.use, c(cells.1, cells.2)],
      cells.1 = cells.1,
      cells.2 = cells.2,
      verbose = verbose
    ),
    'negbinom' = GLMDETest(
      data.use = object[features.use, c(cells.1, cells.2)],
      cells.1 = cells.1,
      cells.2 = cells.2,
      min.cells = min.cells.feature,
      latent.vars = latent.vars,
      test.use = test.use,
      verbose = verbose
    ),
    'poisson' = GLMDETest(
      data.use = object[features.use, c(cells.1, cells.2)],
      cells.1 = cells.1,
      cells.2 = cells.2,
      min.cells = min.cells.feature,
      latent.vars = latent.vars,
      test.use = test.use,
      verbose = verbose
    ),
    'MAST' = MASTDETest(
      data.use = object[features.use, c(cells.1, cells.2)],
      cells.1 = cells.1,
      cells.2 = cells.2,
      latent.vars = latent.vars,
      verbose = verbose
    ),
    "DESeq2" = DESeq2DETest(
      data.use = object[features.use, c(cells.1, cells.2)],
      cells.1 = cells.1,
      cells.2 = cells.2,
      verbose = verbose
    ),
    "LR" = LRDETest(
      data.use = object[features.use, c(cells.1, cells.2)],
      cells.1 = cells.1,
      cells.2 = cells.2,
      verbose = verbose
    ),
    stop("Unknown test: ", test.use)
  )
  de.results[, "avg_logFC"] <- total.diff[rownames(x = de.results)]
  de.results <- cbind(de.results, data.alpha[rownames(x = de.results), , drop = FALSE])
  de.results$p_val_adj = p.adjust(
    p = de.results$p_val,method = "bonferroni",
    n = nrow(object)
  )
  if (test.use == "roc") {
    de.results <- de.results[order(-de.results$power, -de.results$avg_logFC), ]
  } else {
    de.results <- de.results[order(de.results$p_val, -de.results$avg_logFC), ]
  }
  if (only.pos) {
    de.results <- subset(x = de.results, subset = avg_logFC > 0)
  }
  return(de.results)
}


#' @param ident.1 Identity class to define markers for
#' @param ident.2 A second identity class for comparison. If NULL (default) -
#' use all other cells for comparison.
#' @param assay.use Assay to use in differential expression testing
#'
#' @describeIn FindMarkers Run differential expression test on a Seurat object
#' @export
#' @method FindMarkers Seurat
#'
FindMarkers.Seurat <- function(
  object,
  ident.1 = NULL,
  ident.2 = NULL,
  assay.use = NULL,
  features.use = NULL,
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
  assay.use <- assay.use %||% DefaultAssay(object = object)
  data.slot <- "data"
  if (test.use %in% c("negbinom", "poisson", "DESeq2")) {
    data.slot <- "raw.data"
  }
  data.use <- GetAssayData(object = object[[assay.use]], slot = data.slot)
  if (is.null(ident.1)) {
    stop("Please provide ident.1")
  }
  if (length(x = as.vector(x = ident.1)) > 1 &&
      any(as.character(x = ident.1) %in% colnames(data.use))) {
    bad.cells <- colnames(data.use)[which(!as.character(x = ident.1) %in% colnames(data.use))]
    if(length(bad.cells) > 0) {
      stop(paste0("The following cell names provided to ident.1 are not present in the object: ", paste(bad.cells, collapse = ", ")))
    }
  } else {
    ident.1 <- WhichCells(object = object, ident.keep = ident.1)
  }
  # if NULL for ident.2, use all other cells
  if (length(x = as.vector(x = ident.2)) > 1 &&
      any(as.character(x = ident.2) %in% colnames(data.use))) {
    bad.cells <- colnames(data.use)[which(!as.character(x = ident.2) %in% colnames(data.use))]
    if(length(bad.cells) > 0) {
      stop(paste0("The following cell names provided to ident.2 are not present in the object: ", paste(bad.cells, collapse = ", ")))
    }
  } else {
    if (is.null(x = ident.2)) {
      ident.2 <- setdiff(colnames(data.use), ident.1)
    } else {
      ident.2 <- WhichCells(object = object, ident.keep = ident.2)
    }
  }
  if (!is.null(latent.vars)) {
    latent.vars <- FetchData(
      object = object,
      vars.fetch = latent.vars,
      cells.use = c(ident.1, ident.2)
    )
  }
  de.results <- FindMarkers(
    object = data.use,
    cells.1 = ident.1,
    cells.2 = ident.2,
    features.use = features.use,
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
#' @param return.thresh Only return markers that have a p-value < return.thresh, or a power > return.thresh (if the test is ROC)
#'
#' @return Matrix containing a ranked list of putative markers, and associated
#' statistics (p-values, ROC score, etc.)
#'
#' @export
#' @examples
#' all_markers <- FindAllMarkers(object = pbmc_small)
#' head(x = all_markers)
#'
FindAllMarkers <- function(
  object,
  assay.use = NULL,
  features.use = NULL,
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
  return.thresh = 1e-2,
  ...
) {
  if ((test.use == "roc") && (return.thresh == 1e-2)) {
    return.thresh <- 0.7
  }
  idents.all <- sort(x = unique(x = Idents(object = object)))
  genes.de <- list()
  for (i in 1:length(x = idents.all)) {
    genes.de[[i]] <- tryCatch(
      {
        FindMarkers(
          object = object,
          assay.use = assay.use,
          ident.1 = idents.all[i],
          ident.2 = NULL,
          features.use = features.use,
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
      },
      error = function(cond){
        return(NULL)
      }
    )
    if (verbose) {
      message(paste("Calculating cluster", idents.all[i]))
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
      } else {
        gde <- gde[order(gde$p_val, -gde$avg_logFC), ]
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
  if ((only.pos) && nrow(gde.all) > 0) {
    return(subset(x = gde.all, subset = avg_logFC > 0))
  }
  rownames(x = gde.all) <- make.unique(names = as.character(x = gde.all$gene))
  if (nrow(gde.all) == 0) {
    warning("No DE genes identified.")
  }
  return(gde.all)
}

#' Gene expression markers of identity classes defined by a phylogenetic clade
#'
#' Finds markers (differentially expressed genes) based on a branching point (node) in
#' the phylogenetic tree. Markers that define clusters in the left branch are positive markers.
#' Markers that define the right branch are negative markers.
#'
#' @inheritParams FindMarkers
#' @param node The node in the phylogenetic tree to use as a branch point
#' @param tree.use Can optionally pass the tree to be used. Default uses the tree in object@@cluster.tree
#' @param assay.type Type of assay to fetch data for (default is RNA)
#' @param ... Additional arguments passed to FindMarkers
#'
#' @return Matrix containing a ranked list of putative markers, and associated
#' statistics (p-values, ROC score, etc.)
#'
#' @export
#'
#' @examples
#' FindMarkersNode(pbmc_small, 5)
#'
FindMarkersNode <- function(
  object,
  node,
  tree.use = NULL,
  genes.use = NULL,
  logfc.threshold = 0.25,
  test.use = "wilcox",
  assay.type = "RNA",
  ...
) {
  data.use <- GetAssayData(
    object = object,
    assay.type = assay.type
  )
  genes.use <- SetIfNull(x = genes.use, default = rownames(x = data.use))
  tree <- SetIfNull(x = tree.use, default = object@cluster.tree[[1]])
  ident.order <- tree$tip.label
  nodes.1 <- ident.order[GetLeftDescendants(tree = tree, node = node)]
  nodes.2 <- ident.order[GetRightDescendants(tree = tree, node = node)]
  #print(nodes.1)
  #print(nodes.2)
  to.return <- FindMarkers(
    object = object,
    assay.type = assay.type,
    ident.1 = nodes.1,
    ident.2 = nodes.2,
    genes.use = genes.use,
    logfc.threshold = logfc.threshold,
    test.use = test.use,
    ...
  )
  return(to.return)
}

globalVariables(names = c('myAUC', 'p_val'), package = 'Seurat', add = TRUE)
#' Find all markers for a node
#'
#' This function finds markers for all splits at or below the specified node
#'
#' @param object Seurat object. Must have object@@cluster.tree slot filled. Use BuildClusterTree() if not.
#' @param node Node from which to start identifying split markers, default is top node.
#' @param genes.use Genes to test. Default is to use all genes
#' @param logfc.threshold Limit testing to genes which show, on average, at least
#' X-fold difference (log-scale) between the two groups of cells.
#' @param test.use Denotes which test to use. Seurat currently implements
#' "bimod" (likelihood-ratio test for single cell gene expression, McDavid et
#' al., Bioinformatics, 2013, default), "roc" (standard AUC classifier), "t"
#' (Students t-test), and "tobit" (Tobit-test for differential gene expression,
#' as in Trapnell et al., Nature Biotech, 2014), 'poisson', and 'negbinom'.
#' The latter two options should only be used on UMI datasets, and assume an underlying
#' poisson or negative-binomial distribution.
#' @param min.pct - only test genes that are detected in a minimum fraction of min.pct cells
#' in either of the two populations. Meant to speed up the function by not testing genes that are very infrequently expression
#' @param min.diff.pct - only test genes that show a minimum difference in the fraction of detection between the two groups. Set to -Inf by default
#' @param only.pos Only return positive markers (FALSE by default)
#' @param print.bar Print a progress bar once expression testing begins (uses pbapply to do this)
#' @param max.cells.per.ident Down sample each identity class to a max number. Default is no downsampling.
#' @param random.seed Random seed for downsampling
#' @param return.thresh Only return markers that have a p-value < return.thresh, or a power > return.thresh (if the test is ROC)
#' @param do.print Print status updates
#' @param min.cells.gene Minimum number of cells expressing the gene in at least one
#' of the two groups, currently only used for poisson and negative binomial tests
#' @param min.cells.group Minimum number of cells in one of the groups
#' @param assay.type Type of assay to fetch data for (default is RNA)
#' @param \dots Additional parameters to pass to specific DE functions
#'
#' @return Returns a dataframe with a ranked list of putative markers for each node and associated statistics
#'
#' @importFrom ape drop.tip
#'
#' @export
#'
#' @examples
#' pbmc_small
#'
#' FindAllMarkersNode(pbmc_small)
#'
FindAllMarkersNode <- function(
  object,
  node = NULL,
  genes.use = NULL,
  logfc.threshold = 0.25,
  test.use = "wilcox",
  min.pct = 0.1,
  min.diff.pct = 0.05,
  print.bar = TRUE,
  only.pos = FALSE,
  max.cells.per.ident = Inf,
  return.thresh = 1e-2,
  do.print = FALSE,
  random.seed = 1,
  min.cells.gene = 3,
  min.cells.group = 3,
  assay.type = "RNA",
  ...
) {
  if (length(object@cluster.tree) == 0) {
    stop("Tree hasn't been built yet. Run BuildClusterTree to build.")
  }
  data.use <- GetAssayData(object = object,assay.type = assay.type,slot = "data")
  genes.use <- SetIfNull(x = genes.use, default = rownames(x = data.use))
  node <- SetIfNull(x = node, default = object@cluster.tree[[1]]$edge[1, 1])
  tree.use <- object@cluster.tree[[1]]
  descendants <- DFT(tree = tree.use, node = node, path = NULL, include.children = TRUE)
  all.children <- sort(x = tree.use$edge[,2][!tree.use$edge[,2] %in% tree.use$edge[,1]])
  descendants <- MapVals(v = descendants, from = all.children, to = tree.use$tip.label)
  drop.children <- setdiff(tree.use$tip.label, descendants)
  keep.children <- setdiff(tree.use$tip.label, drop.children)
  orig.nodes <- c(node, as.numeric(setdiff(descendants, keep.children)))
  tree.use <- drop.tip(tree.use, drop.children)
  new.nodes <- unique(tree.use$edge[,1])
  if ((test.use == 'roc') && (return.thresh == 1e-2)) {
    return.thresh <- 0.7
  }
  genes.de <- list()
  for (i in ((tree.use$Nnode + 2):max(tree.use$edge))) {
    genes.de[[i]] <- FindMarkersNode(
      object = object,
      assay.type = assay.type,
      node = i,
      tree.use = tree.use,
      genes.use = genes.use,
      logfc.threshold = logfc.threshold,
      test.use = test.use,
      min.pct = min.pct,
      min.diff.pct = min.diff.pct,
      print.bar = print.bar,
      only.pos = only.pos,
      max.cells.per.ident = max.cells.per.ident,
      random.seed = random.seed,
      min.cells.gene = min.cells.gene,
      min.cells.group = min.cells.group
    )
    if (do.print) {
      message(paste("Calculating node", i))
    }
  }
  gde.all <- data.frame()
  for (i in ((tree.use$Nnode + 2):max(tree.use$edge))) {
    if (is.null(x = unlist(x = genes.de[i]))) {
      next
    }
    gde <- genes.de[[i]]
    if (nrow(x = gde) > 0) {
      if (test.use == 'roc') {
        gde <- subset(
          x = gde,
          subset = (myAUC > return.thresh | myAUC < (1 - return.thresh))
        )
      }
      if ( (test.use == 'bimod') || (test.use == 't')) {
        gde <- gde[order(gde$p_val,-gde$avg_logFC), ]
        gde <- subset(x = gde, subset = p_val < return.thresh)
      }
      if (nrow(x = gde) > 0) {
        gde$cluster <- i
        gde$gene <- rownames(x = gde)
      }
      if (nrow(x = gde) > 0) {
        gde.all <- rbind(gde.all,gde)
      }
    }
  }
  gde.all$cluster <- MapVals(
    v = gde.all$cluster,
    from = new.nodes,
    to = orig.nodes
  )
  return(gde.all)
}

#' Finds markers that are conserved between the two groups
#'
#' @param object Seurat object
#' @param ident.1 Identity class to define markers for
#' @param ident.2 A second identity class for comparison. If NULL (default) -
#' use all other cells for comparison.
#' @param grouping.var grouping variable
#' @param assay.type Type of assay to fetch data for (default is RNA)
#' @param meta.method method for combining p-values. Should be a function from
#' the metap package (NOTE: pass the function, not a string)
#' @param \dots parameters to pass to FindMarkers
#'
#' @return Matrix containing a ranked list of putative conserved markers, and
#' associated statistics (p-values within each group and a combined p-value
#' (such as Fishers combined p-value or others from the MetaDE package),
#' percentage of cells expressing the marker, average differences)
#'
#' @import metap
#' @export
#'
#' @examples
#' \dontrun{
#' pbmc_small
#' # Create a simulated grouping variable
#' pbmc_small@meta.data$groups <- sample(
#'   x = c("g1", "g2"),
#'   size = length(x = pbmc_small@cell.names),
#'   replace = TRUE
#' )
#' FindConservedMarkers(pbmc_small, ident.1 = 0, ident.2 = 1, grouping.var = "groups")
#' }
#'
FindConservedMarkers <- function(
  object,
  ident.1,
  ident.2 = NULL,
  grouping.var,
  assay.type = "RNA",
  meta.method = minimump,
  ...
) {
  if(class(meta.method) != "function") {
    stop("meta.method should be a function from the metap package. Please see https://cran.r-project.org/web/packages/metap/metap.pdf for a detail description of the available functions.")
  }
  object.var <- FetchData(object = object, vars.all = grouping.var)
  object <- SetIdent(
    object = object,
    cells.use = object@cell.names,
    ident.use = paste(object@ident, object.var[, 1], sep = "_")
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
  for (i in 1:num.groups) {
    level.use <- levels.split[i]
    ident.use.1 <- paste(ident.1, level.use, sep = "_")
    if(!ident.use.1 %in% object@ident) {
      stop(paste0("Identity: ", ident.1, " not present in group ", level.use))
    }
    cells.1 <- WhichCells(object = object, ident = ident.use.1)
    if (is.null(x = ident.2)) {
      cells.2 <- setdiff(x = cells[[i]], y = cells.1)
      ident.use.2 <- names(x = which(x = table(object@ident[cells.2]) > 0))
      if (length(x = ident.use.2) == 0) {
        stop(paste("Only one identity class present:", ident.1))
      }
    }
    if (! is.null(x = ident.2)) {
      ident.use.2 <- paste(ident.2, level.use, sep = "_")
    }
    cat(
      paste0(
        "Testing ",
        ident.use.1,
        " vs ",
        paste(ident.use.2, collapse = ", "), "\n"
      ),
      file = stderr()
    )
    if(!ident.use.2 %in% object@ident) {
      stop(paste0("Identity: ", ident.2, " not present in group ", level.use))
    }
    marker.test[[i]] <- FindMarkers(
      object = object,
      assay.type = assay.type,
      ident.1 = ident.use.1,
      ident.2 = ident.use.2,
      ...
    )
  }
  genes.conserved <- Reduce(intersect, lapply(marker.test, FUN = function(x) rownames(x)))
  markers.conserved <- list()
  for (i in 1:num.groups) {
    markers.conserved[[i]] <- marker.test[[i]][genes.conserved, ]
    colnames(x = markers.conserved[[i]]) <- paste(
      levels.split[i],
      colnames(x = markers.conserved[[i]]),
      sep="_"
    )
  }
  markers.combined <- Reduce(cbind, markers.conserved)
  pval.codes <- paste(levels.split, "p_val", sep = "_")
  markers.combined$max_pval <- apply(
    X = markers.combined[, pval.codes],
    MARGIN = 1,
    FUN = max
  )
  combined.pval <- data.frame(cp = apply(X = markers.combined[, pval.codes], MARGIN = 1, FUN = function(x) meta.method(x)$p))
  colnames(combined.pval) <- paste0(as.character(formals()$meta.method), "_p_val")
  markers.combined <- cbind(markers.combined, combined.pval)
  markers.combined <- markers.combined[order(markers.combined[,paste0(as.character(formals()$meta.method), "_p_val")]), ]
  return(markers.combined)
}


#' Negative binomial test for UMI-count based data (regularized version)
#'
#' Identifies differentially expressed genes between two groups of cells using
#' a likelihood ratio test of negative binomial generalized linear models where
#' the overdispersion parameter theta is determined by pooling information
#' across genes.
#'
#' @inheritParams FindMarkers
#' @param object Seurat object
#' @param cells.1 Group 1 cells
#' @param cells.2 Group 2 cells
#' @param genes.use Genes to use for test
#' @param latent.vars Latent variables to test
#' @param print.bar Print progress bar
#' @param min.cells Minimum number of cells threshold
#' @param assay.type Type of assay to fetch data for (default is RNA)
#'
#' @return Returns a p-value ranked data frame of test results.
#'
#' @importFrom stats p.adjust
#' @importFrom utils txtProgressBar setTxtProgressBar
#'
#' @export
#'
#' @examples
#' # Note, not recommended for particularly small datasets - expect warnings
#' NegBinomDETest(
#'   object = pbmc_small,
#'   cells.1 = WhichCells(object = pbmc_small, ident = 1),
#'   cells.2 = WhichCells(object = pbmc_small, ident = 2)
#' )
#'
NegBinomRegDETest <- function(
  object,
  cells.1,
  cells.2,
  genes.use = NULL,
  latent.vars = NULL,
  print.bar = TRUE,
  min.cells = 3,
  assay.type = "RNA"
) {
  if (!is.null(genes.use)) {
    message('Make sure that genes.use contains mostly genes that are not expected to be
             differentially expressed to allow unbiased theta estimation')
  }
  genes.use <- SetIfNull(x = genes.use, default = rownames(x = GetAssayData(object = object,assay.type = assay.type,slot = "data")))
  # check that the gene made it through the any filtering that was done
  genes.use <- genes.use[genes.use %in% rownames(x = GetAssayData(object = object,assay.type = assay.type,slot = "data"))]
  message(
    sprintf(
      'NegBinomRegDETest for %d genes and %d and %d cells',
      length(x = genes.use),
      length(x = cells.1),
      length(x = cells.2)
    )
  )
  grp.fac <- factor(
    x = c(
      rep.int(x = 'A', times = length(x = cells.1)),
      rep.int(x = 'B', times = length(x = cells.2))
    )
  )
  to.test.data <- GetAssayData(object = object,assay.type = assay.type,slot = "raw.data")[genes.use, c(cells.1, cells.2), drop = FALSE]
  message('Calculating mean per gene per group')
  above.threshold <- pmax(
    apply(X = to.test.data[, cells.1] > 0, MARGIN = 1, FUN = mean),
    apply(X = to.test.data[, cells.2] > 0, MARGIN = 1, FUN = mean)
  ) >= 0.02
  message(
    sprintf(
      '%d genes are detected in at least 2%% of the cells in at least one of the groups and will be tested',
      sum(above.threshold)
    )
  )
  genes.use <- genes.use[above.threshold]
  to.test.data <- to.test.data[genes.use, , drop = FALSE]
  my.latent <- FetchData(
    object = object,
    vars.all = latent.vars,
    cells.use = c(cells.1, cells.2),
    use.raw = TRUE
  )
  to.test <- data.frame(my.latent, row.names = c(cells.1, cells.2))
  message(paste('Latent variables are', paste(latent.vars, collapse = " ")))
  # get regularized theta (ignoring group factor)
  theta.fit <- RegularizedTheta(
    cm = to.test.data,
    latent.data = to.test,
    min.theta = 0.01,
    bin.size = 128
  )
  message('Running NB regression model comparison')
  to.test$NegBinomRegDETest.group <- grp.fac
  bin.size <- 128
  bin.ind <- ceiling(1:length(x = genes.use) / bin.size)
  max.bin <- max(bin.ind)
  pb <- txtProgressBar(min = 0, max = max.bin, style = 3, file = stderr())
  res <- c()
  for (i in 1:max.bin) {
    genes.bin.use <- genes.use[bin.ind == i]
    bin.out.lst <- parallel::mclapply(
      X = genes.bin.use,
      FUN = function(j) {
        return(NBModelComparison(
          y = to.test.data[j, ],
          theta = theta.fit[j],
          latent.data = to.test,
          com.fac = latent.vars,
          grp.fac = 'NegBinomRegDETest.group'
        ))
      }
    )
    res <- rbind(res, do.call(rbind, bin.out.lst))
    setTxtProgressBar(pb = pb, value = i)
  }
  close(pb)
  rownames(res) <- genes.use
  res <- as.data.frame(x = res)
  res$adj.pval <- p.adjust(p = res$pval, method='fdr')
  res <- res[order(res$pval, -abs(x = res$log2.fc)), ]
  return(res)
}
