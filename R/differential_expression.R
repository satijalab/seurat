#' @include seurat.R
NULL

globalVariables(names = 'avg_logFC', package = 'Seurat', add = TRUE)
#' Gene expression markers of identity classes
#'
#' Finds markers (differentially expressed genes) for identity classes
#'
#' @param object Seurat object
#' @param ident.1 Identity class to define markers for
#' @param ident.2 A second identity class for comparison. If NULL (default) -
#' use all other cells for comparison.
#' @param genes.use Genes to test. Default is to use all genes
#' @param logfc.threshold Limit testing to genes which show, on average, at least
#' X-fold difference (log-scale) between the two groups of cells. Default is 0.25
#' Increasing logfc.threshold speeds up the function, but can miss weaker signals.
#' @param test.use Denotes which test to use. Available options are:
##' \itemize{
##'  \item{"wilcox"} : Wilcoxon rank sum test (default)
##'  \item{"bimod"} : Likelihood-ratio test for single cell gene expression,
##'  (McDavid et al., Bioinformatics, 2013)
##'  \item{"roc"} : Standard AUC classifier
##'  \item{"t"} : Student's t-test
##'  \item{"tobit"} : Tobit-test for differential gene expression (Trapnell et
##'  al., Nature Biotech, 2014)
##'  \item{"poisson"} : Likelihood ratio test assuming an underlying poisson
##'   distribution. Use only for UMI-based datasets
##'  \item{"negbinom"} :  Likelihood ratio test assuming an underlying negative
##'  binomial distribution. Use only for UMI-based datasets
##'  \item{"MAST} : GLM-framework that treates cellular detection rate as a
##'  covariate (Finak et al, Genome Biology, 2015)
##'  \item{"DESeq2} : DE based on a model using the negative binomial
##'  distribution (Love et al, Genome Biology, 2014)
##' }
#' @param min.pct  only test genes that are detected in a minimum fraction of
#' min.pct cells in either of the two populations. Meant to speed up the function
#' by not testing genes that are very infrequently expressed. Default is 0.1
#' @param min.diff.pct  only test genes that show a minimum difference in the
#' fraction of detection between the two groups. Set to -Inf by default
#' @param only.pos Only return positive markers (FALSE by default)
#' @param print.bar Print a progress bar once expression testing begins (uses
#' pbapply to do this)
#' @param max.cells.per.ident Down sample each identity class to a max number.
#' Default is no downsampling. Not activated by default (set to Inf)
#' @param random.seed Random seed for downsampling
#' @param latent.vars Variables to test
#' @param min.cells.gene Minimum number of cells expressing the gene in at least one
#' of the two groups, currently only used for poisson and negative binomial tests
#' @param min.cells.group Minimum number of cells in one of the groups
#' @param pseudocount.use Pseudocount to add to averaged expression values when
#' calculating logFC. 1 by default.
#' @param assay.type Type of assay to fetch data for (default is RNA)
#' @param \dots Additional parameters to pass to specific DE functions
#' @seealso \code{\link{MASTDETest}}, and \code{\link{DESeq2DETest}} for more information on these methods
#' @return Matrix containing a ranked list of putative markers, and associated
#' statistics (p-values, ROC score, etc.)
#' @details p-value adjustment is performed using bonferroni correction based on
#' the total number of genes in the dataset. Other correction methods are not
#' recommended, as Seurat pre-filters genes using the arguments above, reducing
#' the number of tests performed. Lastly, as Aaron Lun has pointed out, p-values
#' should be interpreted cautiously, as the genes used for clustering are the
#' same genes tested for differential expression.
#' @import pbapply
#' @importFrom lmtest lrtest
#'
#' @seealso \code{\link{NegBinomDETest}}
#'
#' @export
#'
#' @examples
#' markers <- FindMarkers(object = pbmc_small, ident.1 = 3)
#' head(markers)
#'
FindMarkers <- function(
  object,
  ident.1,
  ident.2 = NULL,
  genes.use = NULL,
  logfc.threshold = 0.25,
  test.use = "wilcox",
  min.pct = 0.1,
  min.diff.pct = -Inf,
  print.bar = TRUE,
  only.pos = FALSE,
  max.cells.per.ident = Inf,
  random.seed = 1,
  latent.vars = "nUMI",
  min.cells.gene = 3,
  min.cells.group = 3,
  pseudocount.use = 1,
  assay.type = "RNA",
  ...
) {
  data.use <- GetAssayData(object = object,assay.type = assay.type,slot = "data")
  genes.use <- SetIfNull(x = genes.use, default = rownames(x = data.use))
  methods.noprefiliter <- c("DESeq2", "zingeR")
  if (test.use %in% methods.noprefiliter) {
    genes.use <- rownames(x = data.use)
    min.diff.pct <- -Inf
    logfc.threshold <- 0
  }
  # in case the user passed in cells instead of identity classes
  if (length(x = as.vector(x = ident.1) > 1) && any(as.character(x = ident.1) %in% object@cell.names)) {
    cells.1 <- intersect(x = ident.1, y = object@cell.names)
  } else {
    cells.1 <- WhichCells(object = object, ident = ident.1)
  }
  # if NULL for ident.2, use all other cells
  if (length(x = as.vector(x = ident.2) > 1) && any(as.character(x = ident.2) %in% object@cell.names)) {
    cells.2 <- intersect(x = ident.2, y = object@cell.names)
  } else {
    if (is.null(x = ident.2)) {
      # cells.2 <- object@cell.names
      cells.2 <- WhichCells(object = object,cells.use = setdiff(object@cell.names,cells.1))
    } else {
      cells.2 <- WhichCells(object = object, ident = ident.2)
    }
  }
  cells.2 <- setdiff(x = cells.2, y = cells.1)
  # error checking
  if (length(x = cells.1) == 0) {
    print(paste("Cell group 1 is empty - no cells with identity class", ident.1))
    return(NULL)
  }
  if (length(x = cells.2) == 0) {
    print(paste("Cell group 2 is empty - no cells with identity class", ident.2))
    return(NULL)
  }
  if (length(cells.1) < min.cells.group) {
    stop(paste("Cell group 1 has fewer than", as.character(min.cells.group), "cells in identity class", ident.1))
  }
  if (length(cells.2) < min.cells.group) {
    stop(paste("Cell group 2 has fewer than", as.character(min.cells.group), " cells in identity class", ident.2))
  }
  # gene selection (based on percent expressed)
  thresh.min <- 0
  data.temp1 <- round(
    x = apply(
      X = data.use[genes.use, cells.1, drop = F],
      MARGIN = 1,
      FUN = function(x) {
        return(sum(x > thresh.min) / length(x = x))
        # return(length(x = x[x>thresh.min]) / length(x = x))
      }
    ),
    digits = 3
  )
  data.temp2 <- round(
    x = apply(
      X = data.use[genes.use, cells.2, drop = F],
      MARGIN = 1,
      FUN = function(x) {
        return(sum(x > thresh.min) / length(x = x))
        # return(length(x = x[x > thresh.min]) / length(x = x))
      }
    ),
    digits = 3
  )
  data.alpha <- cbind(data.temp1, data.temp2)
  colnames(x = data.alpha) <- c("pct.1","pct.2")
  alpha.min <- apply(X = data.alpha, MARGIN = 1, FUN = max)
  names(x = alpha.min) <- rownames(x = data.alpha)
  genes.use <- names(x = which(x = alpha.min > min.pct))
  if (length(x = genes.use) == 0) {
    stop("No genes pass min.pct threshold")
  }
  alpha.diff <- alpha.min - apply(X = data.alpha, MARGIN = 1, FUN = min)
  genes.use <- names(
    x = which(x = alpha.min > min.pct & alpha.diff > min.diff.pct)
  )
  if (length(x = genes.use) == 0) {
    stop("No genes pass min.diff.pct threshold")
  }

  #gene selection (based on average difference)
  data.1 <- apply(X = data.use[genes.use, cells.1, drop = F], MARGIN = 1, FUN = function(x) log(x = mean(x = expm1(x = x)) + pseudocount.use))
  data.2 <- apply(X = data.use[genes.use, cells.2, drop = F], MARGIN = 1, FUN = function(x) log(x = mean(x = expm1(x = x)) + pseudocount.use))
  total.diff <- (data.1 - data.2)
  if (!only.pos) genes.diff <- names(x = which(x = abs(x = total.diff) > logfc.threshold))
  if (only.pos) genes.diff <- names(x = which(x = total.diff > logfc.threshold))
  genes.use <- intersect(x = genes.use, y = genes.diff)
  if (length(x = genes.use) == 0) {
    stop("No genes pass logfc.threshold threshold")
  }
  if (max.cells.per.ident < Inf) {
    set.seed(seed = random.seed)
    if (length(cells.1) > max.cells.per.ident) cells.1 = sample(x = cells.1, size = max.cells.per.ident)
    if (length(cells.2) > max.cells.per.ident) cells.2 = sample(x = cells.2, size = max.cells.per.ident)
  }
  #perform DR
  if (test.use == "bimod") {
    to.return <- DiffExpTest(
      object = object,
      assay.type = assay.type,
      cells.1 = cells.1,
      cells.2 = cells.2,
      genes.use = genes.use,
      print.bar = print.bar
    )
  }
  if (test.use == "roc") {
    to.return <- MarkerTest(
      object = object,
      assay.type = assay.type,
      cells.1 = cells.1,
      cells.2 = cells.2,
      genes.use = genes.use,
      print.bar = print.bar
    )
  }
  if (test.use == "t") {
    to.return <- DiffTTest(
      object = object,
      assay.type = assay.type,
      cells.1 = cells.1,
      cells.2 = cells.2,
      genes.use = genes.use,
      print.bar = print.bar
    )
  }
  if (test.use == "tobit") {
    to.return <- TobitTest(
      object = object,
      assay.type = assay.type,
      cells.1 = cells.1,
      cells.2 = cells.2,
      genes.use = genes.use,
      print.bar = print.bar
    )
  }
  if (test.use == "negbinom") {
    to.return <- NegBinomDETest(
      object = object,
      assay.type = assay.type,
      cells.1 = cells.1,
      cells.2 = cells.2,
      genes.use = genes.use,
      latent.vars = latent.vars,
      print.bar = print.bar,
      min.cells = min.cells.gene
    )
  }
  if (test.use == "poisson") {
    to.return <- PoissonDETest(
      object = object,
      assay.type = assay.type,
      cells.1 = cells.1,
      cells.2 = cells.2,
      genes.use = genes.use,
      latent.vars = latent.vars,
      print.bar = print.bar,
      min.cells = min.cells.gene
    )
  }
  if (test.use == "MAST") {
    to.return <- MASTDETest(
      object = object,
      assay.type = assay.type,
      cells.1 = cells.1,
      cells.2 = cells.2,
      genes.use = genes.use,
      latent.vars = latent.vars,
      ...
    )
  }

  if (test.use == "wilcox") {
    to.return <- WilcoxDETest(
      object = object,
      assay.type = assay.type,
      cells.1 = cells.1,
      cells.2 = cells.2,
      genes.use = genes.use,
      print.bar = print.bar,
      ...
    )
  }
  if (test.use == "LR") {
    to.return <- LRDETest(
      object = object,
      assay.type = assay.type,
      cells.1 = cells.1,
      cells.2 = cells.2,
      genes.use = genes.use,
      print.bar = print.bar,
      ...
    )
  }
    if (test.use == "DESeq2") {
      to.return <- DESeq2DETest(
        object = object,
        assay.type = assay.type,
        cells.1 = cells.1,
        cells.2 = cells.2,
        genes.use = genes.use,
        ...
      )
  }
  #return results
  to.return[, "avg_logFC"] <- total.diff[rownames(x = to.return)]
  to.return <- cbind(to.return, data.alpha[rownames(x = to.return), ])
  to.return$p_val_adj = p.adjust(p = to.return$p_val,method = "bonferroni",
                                 n =   nrow(GetAssayData(object = object,
                                                         assay.type = assay.type,
                                                         slot = "data")))
  if (test.use == "roc") {
    to.return <- to.return[order(-to.return$power, -to.return$avg_logFC), ]
  } else {
    to.return <- to.return[order(to.return$p_val, -to.return$avg_logFC), ]
  }
  if (only.pos) {
    to.return <- subset(x = to.return, subset = avg_logFC > 0)
  }
  return(to.return)
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
#' @param print.bar Print a progress bar once expression testing begins (uses pbapply to do this)
#' @param max.cells.per.ident Down sample each identity class to a max number. Default is no downsampling.
#' @param random.seed Random seed for downsampling
#' @param return.thresh Only return markers that have a p-value < return.thresh, or a power > return.thresh (if the test is ROC)
#' @param do.print FALSE by default. If TRUE, outputs updates on progress.
#' @param min.cells.gene Minimum number of cells expressing the gene in at least one
#' of the two groups, currently only used for poisson and negative binomial tests
#' @param min.cells.group Minimum number of cells in one of the groups
#' @param latent.vars remove the effects of these variables
#' @param assay.type Type of assay to perform DE for (default is RNA)
#' @param \dots Additional parameters to pass to specific DE functions
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
  genes.use = NULL,
  logfc.threshold = 0.25,
  test.use = "wilcox",
  min.pct = 0.1,
  min.diff.pct = -Inf,
  print.bar = TRUE,
  only.pos = FALSE,
  max.cells.per.ident = Inf,
  return.thresh = 1e-2,
  do.print = FALSE,
  random.seed = 1,
  min.cells.gene = 3,
  min.cells.group = 3,
  latent.vars = "nUMI",
  assay.type = "RNA",
  ...
) {
  data.1 <- GetAssayData(object = object,assay.type = assay.type,slot = "data")
  genes.use <- SetIfNull(x = genes.use, default = rownames(x = data.1))
  if ((test.use == "roc") && (return.thresh == 1e-2)) {
    return.thresh = 0.7
  }
  idents.all <- sort(x = unique(x = object@ident))
  genes.de <- list()
  #if (max.cells.per.ident < Inf) {
  #  object <- SubsetData(
  #    object = object,
  #    max.cells.per.ident = max.cells.per.ident,
  #    random.seed = random.seed
  #  )
  #}
  for (i in 1:length(x = idents.all)) {
    genes.de[[i]] <- tryCatch(
      {
        FindMarkers(
          object = object,
          assay.type = assay.type,
          ident.1 = idents.all[i],
          ident.2 = NULL,
          genes.use = genes.use,
          logfc.threshold = logfc.threshold,
          test.use = test.use,
          min.pct = min.pct,
          min.diff.pct = min.diff.pct,
          print.bar = print.bar,
          min.cells.gene = min.cells.gene,
          min.cells.group = min.cells.group,
          latent.vars = latent.vars,
          max.cells.per.ident = max.cells.per.ident,
          ...
        )
      },
      error = function(cond){
        return(NULL)
      }
    )
    if (do.print) {
      print(paste("Calculating cluster", idents.all[i]))
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
  if(nrow(gde.all) == 0) {
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
      print(paste("Calculating node", i))
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

#' Likelihood ratio test for zero-inflated data
#'
#' Identifies differentially expressed genes between two groups of cells using
#' the LRT model proposed in McDavid et al, Bioinformatics, 2013
#'
#' @inheritParams FindMarkers
#' @param object Seurat object
#' @param cells.1 Group 1 cells
#' @param cells.2 Group 2 cells
#' @param assay.type Type of assay to fetch data for (default is RNA)
#' @return Returns a p-value ranked matrix of putative differentially expressed
#' genes.
#'
#' @export
#' @examples
#' pbmc_small
#' DiffExpTest(pbmc_small, cells.1 = WhichCells(object = pbmc_small, ident = 1),
#'             cells.2 = WhichCells(object = pbmc_small, ident = 2))
#'
DiffExpTest <- function(
  object,
  cells.1,
  cells.2,
  assay.type = "RNA",
  genes.use = NULL,
  print.bar = TRUE
) {
  data.test <- GetAssayData(object = object,assay.type = assay.type,slot = "data")
  genes.use <- SetIfNull(x = genes.use, default = rownames(data.test))
  if (print.bar) {
    iterate.fxn <- pblapply
  } else {
    iterate.fxn <- lapply
  }
  p_val <- unlist(
    x = iterate.fxn(
      X = genes.use,
      FUN = function(x) {
        return(
          DifferentialLRT(
            x = as.numeric(x = data.test[x, cells.1]),
            y = as.numeric(x = data.test[x, cells.2])
          )
        )
      }
    )
  )
  to.return <- data.frame(p_val, row.names = genes.use)
  return(to.return)
}

#' Negative binomial test for UMI-count based data
#'
#' Identifies differentially expressed genes between two groups of cells using
#' a negative binomial generalized linear model
#'
#' @param object Seurat object
#' @param cells.1 Group 1 cells
#' @param cells.2 Group 2 cells
#' @param genes.use Genes to use for test
#' @param latent.vars Latent variables to test
#' @param print.bar Print progress bar
#' @param min.cells Minimum number of cells threshold
#' @param assay.type Type of assay to fetch data for (default is RNA)
#'
#' @return Returns a p-value ranked matrix of putative differentially expressed
#' genes.
#'
#' @importFrom MASS glm.nb
#' @importFrom pbapply pbapply
#' @importFrom stats var as.formula
#'
#' @export
#'
#'@examples
#' pbmc_small
#' # Note, not recommended for particularly small datasets - expect warnings
#' NegBinomDETest(pbmc_small, cells.1 = WhichCells(object = pbmc_small, ident = 1),
#'             cells.2 = WhichCells(object = pbmc_small, ident = 2))
#'
NegBinomDETest <- function(
  object,
  cells.1,
  cells.2,
  genes.use = NULL,
  latent.vars = NULL,
  print.bar = TRUE,
  min.cells = 3,
  assay.type = "RNA"
) {
  genes.use <- SetIfNull(x = genes.use, default = rownames(x = GetAssayData(object = object,assay.type = assay.type,slot = "data")))
  # check that the gene made it through the any filtering that was done
  genes.use <- genes.use[genes.use %in% rownames(x = GetAssayData(object = object,assay.type = assay.type,slot = "data"))]
  my.latent <- FetchData(
    object = object,
    vars.all = latent.vars,
    cells.use = c(cells.1, cells.2),
    use.raw = TRUE
  )
  to.test.data <- GetAssayData(object = object,assay.type = assay.type,slot = "raw.data")[genes.use, c(cells.1, cells.2)]
  to.test <- data.frame(my.latent, row.names = c(cells.1, cells.2))
  to.test[cells.1, "group"] <- "A"
  to.test[cells.2, "group"] <- "B"
  to.test$group <- factor(x = to.test$group)
  latent.vars <- c("group", latent.vars)
  if (print.bar) {
    iterate.fxn <- pblapply
  } else {
    iterate.fxn <- lapply
  }
  p_val <- unlist(
    x = iterate.fxn(
      X = genes.use,
      FUN = function(x) {
        to.test[, "GENE"] <- as.numeric(x = to.test.data[x, ])
        # check that gene is expressed in specified number of cells in one group
        if (sum(to.test$GENE[to.test$group == "A"]) < min.cells ||
            sum(to.test$GENE[to.test$group == "B"]) < min.cells) {
          warning(paste0(
            "Skipping gene --- ",
            x,
            ". Fewer than ",
            min.cells,
            " in at least one of the two clusters."
          ))
          return(2)
        }
        # check that variance between groups is not 0
        if (var(x = to.test$GENE) == 0) {
          warning(paste0(
            "Skipping gene -- ",
            x,
            ". No variance in expression between the two clusters."
          ))
          return(2)
        }
        fmla <- as.formula(paste0("GENE ", " ~ ", paste(latent.vars, collapse = "+")))
        p.estimate <- 2
        try(
          expr = p.estimate <- summary(
            object = glm.nb(formula = fmla, data = to.test)
          )$coef[2, 4],
          silent = TRUE
        )
        return(p.estimate)
      }
    )
  )
  if (length(x = which(x = p_val == 2)) > 0){
    genes.use <- genes.use[-which(x = p_val == 2)]
    p_val <- p_val[! p_val == 2]
  }
  to.return <- data.frame(p_val, row.names = genes.use)
  return(to.return)
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
  genes.use <- SetIfNull(x = genes.use, default = rownames(x = GetAssayData(object = object,assay.type = assay.type,slot = "data")))
  # check that the gene made it through the any filtering that was done
  genes.use <- genes.use[genes.use %in% rownames(x = GetAssayData(object = object,assay.type = assay.type,slot = "data"))]
  print(
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
  print('Calculating mean per gene per group')
  above.threshold <- pmax(
    apply(X = to.test.data[, cells.1] > 0, MARGIN = 1, FUN = mean),
    apply(X = to.test.data[, cells.2] > 0, MARGIN = 1, FUN = mean)
  ) >= 0.02
  print(
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
  print(paste('Latent variables are', latent.vars))
  # get regularized theta (ignoring group factor)
  theta.fit <- RegularizedTheta(
    cm = to.test.data,
    latent.data = to.test,
    min.theta = 0.01,
    bin.size = 128
  )
  print('Running NB regression model comparison')
  to.test$NegBinomRegDETest.group <- grp.fac
  bin.size <- 128
  bin.ind <- ceiling(1:length(x = genes.use) / bin.size)
  max.bin <- max(bin.ind)
  pb <- txtProgressBar(min = 0, max = max.bin, style = 3)
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
  res <- res[order(res$pval, -abs(x = res$log.fc)), ]
  return(res)
}

globalVariables(names = 'min.cells', package = 'Seurat', add = TRUE)
#' Poisson test for UMI-count based data
#'
#' Identifies differentially expressed genes between two groups of cells using
#' a poisson generalized linear model
#
#' @param object Seurat object
#' @param cells.1 Group 1 cells
#' @param cells.2 Group 2 cells
#' @param min.cells Minimum number of cells expressing the gene in at least one of the two groups
#' @param genes.use Genes to use for test
#' @param latent.vars Latent variables to test
#' @param print.bar Print progress bar
#' @param assay.type Type of assay to fetch data for (default is RNA)
#'
#' @return Returns a p-value ranked matrix of putative differentially expressed
#' genes.
#'
#' @importFrom pbapply pbapply
#' @importFrom stats var as.formula glm
#'
#' @export
#'
#'@examples
#' pbmc_small
#' # Note, expect warnings with example dataset due to min.cells threshold.
#' PoissonDETest(pbmc_small, cells.1 = WhichCells(object = pbmc_small, ident = 1),
#'             cells.2 = WhichCells(object = pbmc_small, ident = 2))
#'
PoissonDETest <- function(
  object,
  cells.1,
  cells.2,
  min.cells = 3,
  genes.use = NULL,
  latent.vars = NULL,
  print.bar = TRUE,
  assay.type = "RNA"
) {
  genes.use <- SetIfNull(x = genes.use, default = rownames(x = GetAssayData(object = object,assay.type = assay.type,slot = "data")))
  # check that the gene made it through the any filtering that was done
  genes.use <- genes.use[genes.use %in% rownames(x = GetAssayData(object = object,assay.type = assay.type,slot = "data"))]
  my.latent <- FetchData(
    object = object,
    vars.all = latent.vars,
    cells.use = c(cells.1, cells.2),
    use.raw = TRUE
  )
  to.test.data <- GetAssayData(object = object,assay.type = assay.type,slot = "raw.data")[genes.use, c(cells.1, cells.2)]
  to.test <- data.frame(my.latent, row.names = c(cells.1, cells.2))
  to.test[cells.1,"group"] <- "A"
  to.test[cells.2,"group"] <- "B"
  to.test$group <- factor(x = to.test$group)
  latent.vars <- c("group", latent.vars)
  if (print.bar) {
    iterate.fxn <- pblapply
  } else {
    iterate.fxn <- lapply
  }
  p_val <- unlist(
    x = iterate.fxn(
      X = genes.use,
      FUN = function(x) {
        to.test[,"GENE"] <- as.numeric(x = to.test.data[x, ])
        # check that gene is expressed in specified number of cells in one group
        if (sum(to.test$GENE[to.test$group == "A"]) < min.cells ||
            sum(to.test$GENE[to.test$group == "B"]) < min.cells) {
          warning(paste0(
            "Skipping gene --- ",
            x,
            ". Fewer than",
            min.cells,
            " in at least one of the two clusters."
          ))
          return(2)
        }
        # check that variance between groups is not 0
        if (var(to.test$GENE) == 0) {
          print("what") # what?
          warning(paste0(
            "Skipping gene -- ",
            x,
            ". No variance in expression between the two clusters."
          ))
          return(2)
        }
        fmla <- as.formula(
          object = paste0("GENE ", " ~ ", paste(latent.vars, collapse="+"))
        )
        return(
          summary(
            object = glm(
              formula = fmla,
              data = to.test,
              family = "poisson"
            )
          )$coef[2,4]
        )
      }
    )
  )
  if (length(x = which(x = p_val == 2)) > 0) {
    genes.use <- genes.use[-which(x = p_val == 2)]
    p_val <- p_val[! p_val == 2]
  }
  to.return <- data.frame(p_val, row.names = genes.use)
  return(to.return)
}

#' Differential expression using MAST
#'
#' Identifies differentially expressed genes between two groups of cells using
#' a hurdle model tailored to scRNA-seq data. Utilizes the MAST package to run
#' the DE testing.
#'
#' @references Andrew McDavid, Greg Finak and Masanao Yajima (2017). MAST: Model-based
#' Analysis of Single Cell Transcriptomics. R package version 1.2.1.
#' https://github.com/RGLab/MAST/
#' @param object Seurat object
#' @param cells.1 Group 1 cells
#' @param cells.2 Group 2 cells
#' @param genes.use Genes to use for test
#' @param latent.vars Confounding variables to adjust for in DE test. Default is
#' "nUMI", which adjusts for cellular depth (i.e. cellular detection rate). For
#' non-UMI based data, set to nGene instead.
#' @param assay.type Type of assay to fetch data for (default is RNA)
#' @param \dots Additional parameters to zero-inflated regression (zlm) function
#' in MAST
#' @details
#' To use this method, please install MAST, using instructions at https://github.com/RGLab/MAST/
#'
#' @return Returns a p-value ranked matrix of putative differentially expressed
#' genes.
#'
#' @importFrom stats relevel
#' @importFrom utils installed.packages
#'
#' @export
#'
#' @examples
#' \dontrun{
#'   pbmc_small
#'   MASTDETest(pbmc_small, cells.1 = WhichCells(object = pbmc_small, ident = 1),
#'               cells.2 = WhichCells(object = pbmc_small, ident = 2))
#' }
#'
MASTDETest <- function(
  object,
  cells.1,
  cells.2,
  genes.use = NULL,
  latent.vars = NULL,
  assay.type = "RNA",
  ...
) {
  # Check for MAST
  if (!'MAST' %in% rownames(x = installed.packages())) {
    stop("Please install MAST - learn more at https://github.com/RGLab/MAST")
  }
  data.test <- GetAssayData(object = object,assay.type = assay.type,slot = "data")
  genes.use <- SetIfNull(x = genes.use, default = rownames(x = data.test))
  # check that the gene made it through the any filtering that was done
  genes.use <- genes.use[genes.use %in% rownames(x = data.test)]
  my.latent <- FetchData(
    object = object,
    vars.all = latent.vars,
    cells.use = c(cells.1, cells.2),
    use.raw = TRUE
  )
  if (length(x = latent.vars) > 0) {
    my.latent <- scale(x = my.latent)
  }
  coldata <- data.frame(my.latent, row.names = c(cells.1, cells.2))
  coldata[cells.1, "group"] <- "Group1"
  coldata[cells.2, "group"] <- "Group2"
  coldata$group <- factor(x = coldata$group)
  coldata$wellKey <- rownames(x = coldata)
  latent.vars <- c("condition", latent.vars)
  countdata.test <- data.test[genes.use, rownames(x = coldata)]
  fdat <- data.frame(rownames(x = countdata.test))
  colnames(x = fdat)[1] <- "primerid"
  rownames(x = fdat) <- fdat[, 1]
  sca <- MAST::FromMatrix(
    exprsArray = as.matrix(x = countdata.test),
    cData = coldata,
    fData = fdat
  )
  cond <- factor(x = SummarizedExperiment::colData(sca)$group)
  cond <- relevel(x = cond, ref = "Group1")
  SummarizedExperiment::colData(sca)$condition <- cond
  fmla <- as.formula(
    object = paste0(" ~ ", paste(latent.vars, collapse = "+"))
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

#' Differential expression using DESeq2
#'
#' Identifies differentially expressed genes between two groups of cells using
#' DESeq2
#'
#' @references Love MI, Huber W and Anders S (2014). "Moderated estimation of
#' fold change and dispersion for RNA-seq data with DESeq2." Genome Biology.
#' https://bioconductor.org/packages/release/bioc/html/DESeq2.html
#' @param object Seurat object
#' @param cells.1 Group 1 cells
#' @param cells.2 Group 2 cells
#' @param genes.use Genes to use for test
#' @param assay.type Type of assay to fetch data for (default is RNA)
#' @param ... Extra parameters to pass to DESeq2::results
#' @return Returns a p-value ranked matrix of putative differentially expressed
#' genes.
#'
#' @details
#' This test does not support pre-filtering of genes based on average difference
#' (or percent detection rate) between cell groups. However, genes may be
#' pre-filtered based on their minimum detection rate (min.pct) across both cell
#' groups. To use this method, please install DESeq2, using the instructions at
#'  https://bioconductor.org/packages/release/bioc/html/DESeq2.html
#'
#' @importFrom utils installed.packages
#'
#' @export
#'
#' @examples
#' \dontrun{
#'   pbmc_small
#'   DESeq2DETest(pbmc_small, cells.1 = WhichCells(object = pbmc_small, ident = 1),
#'               cells.2 = WhichCells(object = pbmc_small, ident = 2))
#' }
#'
DESeq2DETest <- function(
  object,
  cells.1,
  cells.2,
  genes.use = NULL,
  assay.type = "RNA",
  ...
) {
  if (!'DESeq2' %in% rownames(x = installed.packages())) {
    stop("Please install DESeq2 - learn more at https://bioconductor.org/packages/release/bioc/html/DESeq2.html")
  }
  genes.use <- SetIfNull(x = genes.use, default = rownames(x = GetAssayData(object = object,assay.type = assay.type,slot = "data")))
  # check that the gene made it through the any filtering that was done
  genes.use <- genes.use[genes.use %in% rownames(x = GetAssayData(object = object,assay.type = assay.type,slot = "data"))]
  coldata <- object@meta.data[c(cells.1, cells.2), ]
  coldata[cells.1, "group"] <- "Group1"
  coldata[cells.2, "group"] <- "Group2"
  coldata$group <- factor(x = coldata$group)
  coldata$wellKey <- rownames(x = coldata)
  countdata.test <- GetAssayData(object = object,assay.type = assay.type,slot = "raw.data")[genes.use, rownames(x = coldata)]
  fdat <- data.frame(rownames(x = countdata.test))
  colnames(x = fdat)[1] <- "primerid"
  rownames(x = fdat) <- fdat[, 1]
  dds1 <- DESeq2::DESeqDataSetFromMatrix(
    countData = countdata.test,
    colData = coldata,
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
  p_val <- res$pvalue
  genes.return <- rownames(x = res)
  to.return <- data.frame(p_val, row.names = genes.return)
  return(to.return)
}

#' Differential expression using Wilcoxon Rank Sum
#'
#' Identifies differentially expressed genes between two groups of cells using
#' a Wilcoxon Rank Sum test
#' @param object Seurat object
#' @param cells.1 Group 1 cells
#' @param cells.2 Group 2 cells
#' @param genes.use Genes to use for test
#' @param print.bar Print a progress bar
#' @param assay.type Type of assay to perform DE for (default is RNA)
#' @param ... Extra parameters passed to wilcox.test
#'
#' @return Returns a p-value ranked matrix of putative differentially expressed
#' genes.
#'
#' @importFrom pbapply pbsapply
#' @importFrom stats wilcox.test
#'
#' @export
#'
#' @examples
#' pbmc_small
#' WilcoxDETest(pbmc_small, cells.1 = WhichCells(object = pbmc_small, ident = 1),
#'             cells.2 = WhichCells(object = pbmc_small, ident = 2))
#'
WilcoxDETest <- function(
  object,
  cells.1,
  cells.2,
  genes.use = NULL,
  print.bar = TRUE,
  assay.type = "RNA",
  ...
) {
  data.test <- GetAssayData(object = object,assay.type = assay.type,slot = "data")
  genes.use <- SetIfNull(x = genes.use, default = rownames(x = data.test))
  # check that the gene made it through the any filtering that was done
  genes.use <- genes.use[genes.use %in% rownames(x = data.test)]
  coldata <- object@meta.data[c(cells.1, cells.2), ]
  coldata[cells.1, "group"] <- "Group1"
  coldata[cells.2, "group"] <- "Group2"
  coldata$group <- factor(x = coldata$group)
  coldata$wellKey <- rownames(x = coldata)
  countdata.test <- data.test[genes.use, rownames(x = coldata)]
  mysapply <- if (print.bar) {pbsapply} else {sapply}
  p_val <- mysapply(
    X = 1:nrow(x = countdata.test),
    FUN = function(x) {
      return(wilcox.test(countdata.test[x, ] ~ coldata$group, ...)$p.value)
    }
  )
  genes.return <- rownames(x = countdata.test)
  to.return <- data.frame(p_val, row.names = genes.return)
  return(to.return)
}


LRDETest <- function(
  object,
  cells.1,
  cells.2,
  min.cells = 3,
  genes.use = NULL,
  print.bar = TRUE,
  assay.type = "RNA",
  ...
) {
  data.test <- GetAssayData(object = object, assay.type = assay.type, slot = "data")
  genes.use <- SetIfNull(x = genes.use, default = rownames(x = data.test))
  # check that the gene made it through the any filtering that was done
  genes.use <- genes.use[genes.use %in% rownames(x = data.test)]
  coldata <- object@meta.data[c(cells.1, cells.2), ]
  coldata[cells.1, "group"] <- "Group1"
  coldata[cells.2, "group"] <- "Group2"
  coldata$group <- factor(x = coldata$group)
  coldata$wellKey <- rownames(x = coldata)
  countdata.test <- data.test[genes.use, rownames(x = coldata)]
  mysapply <- if (print.bar) {pbsapply} else {sapply}
  p_val <- mysapply(
    X = 1:nrow(x = countdata.test),
    FUN = function(x) {
      model1 <- glm(coldata$group ~ countdata.test[x,],family = "binomial")
      model2 <- glm(coldata$group ~ 1, family = "binomial")
      lrtest <- lrtest(model1, model2)
      return(lrtest$Pr[2])
    }
  )
  genes.return <- rownames(x = countdata.test)
  to.return <- data.frame(p_val, row.names = genes.return)
  return(to.return)
}

#' Differential expression testing using Tobit models
#'
#' Identifies differentially expressed genes between two groups of cells using
#' Tobit models, as proposed in Trapnell et al., Nature Biotechnology, 2014
#'
#' @inheritParams FindMarkers
#' @inheritParams DiffExpTest
#'
#' @return Returns a p-value ranked matrix of putative differentially expressed
#' genes.
#'
#' @export
#'
#'@examples
#' pbmc_small
#' \dontrun{
#' TobitTest(pbmc_small, cells.1 = WhichCells(object = pbmc_small, ident = 1),
#'             cells.2 = WhichCells(object = pbmc_small, ident = 2))
#' }
#'
TobitTest <- function(
  object,
  cells.1,
  cells.2,
  genes.use = NULL,
  print.bar = TRUE,
  assay.type = "RNA"
) {
  data.test <- GetAssayData(object = object,assay.type = assay.type,slot = "data")
  genes.use <- SetIfNull(x = genes.use, default = rownames(x = data.test))
  #print(genes.diff)
  to.return <- TobitDiffExpTest(
    data1 = data.test[, cells.1],
    data2 = data.test[, cells.2],
    mygenes = genes.use,
    print.bar = print.bar
  )
  return(to.return)
}

#' ROC-based marker discovery
#'
#' Identifies 'markers' of gene expression using ROC analysis. For each gene,
#' evaluates (using AUC) a classifier built on that gene alone, to classify
#' between two groups of cells.
#'
#' An AUC value of 1 means that expression values for this gene alone can
#' perfectly classify the two groupings (i.e. Each of the cells in cells.1
#' exhibit a higher level than each of the cells in cells.2). An AUC value of 0
#' also means there is perfect classification, but in the other direction. A
#' value of 0.5 implies that the gene has no predictive power to classify the
#' two groups.
#'
#' @inheritParams FindMarkers
#' @inheritParams DiffExpTest
#' @param object Seurat object
#'
#' @return Returns a 'predictive power' (abs(AUC-0.5)) ranked matrix of
#' putative differentially expressed genes.
#'
#' @export
#'
#' @examples
#' pbmc_small
#' MarkerTest(pbmc_small, cells.1 = WhichCells(object = pbmc_small, ident = 1),
#'             cells.2 = WhichCells(object = pbmc_small, ident = 2))
#'
MarkerTest <- function(
  object,
  cells.1,
  cells.2,
  genes.use = NULL,
  print.bar = TRUE,
  assay.type = "RNA"
) {
  data.test <- GetAssayData(object = object,assay.type = assay.type,slot = "data")
  genes.use <- SetIfNull(x = genes.use, default = rownames(x = data.test))
  to.return <- AUCMarkerTest(
    data1 = data.test[, cells.1],
    data2 = data.test[, cells.2],
    mygenes = genes.use,
    print.bar = print.bar
  )
  to.return$power <- abs(x = to.return$myAUC - 0.5) * 2
  #print(head(to.return))
  return(to.return)
}

#' Differential expression testing using Student's t-test
#'
#' Identify differentially expressed genes between two groups of cells using
#' the Student's t-test
#'
#' @inheritParams FindMarkers
#' @inheritParams DiffExpTest
#'
#' @return Returns a p-value ranked matrix of putative differentially expressed
#' genes.
#'
#' @importFrom stats t.test
#' @importFrom pbapply pblapply
#'
#' @export
#'
#' @examples
#' pbmc_small
#' DiffTTest(pbmc_small, cells.1 = WhichCells(object = pbmc_small, ident = 1),
#'             cells.2 = WhichCells(object = pbmc_small, ident = 2))
DiffTTest <- function(
  object,
  cells.1,
  cells.2,
  genes.use = NULL,
  print.bar = TRUE,
  assay.type = "RNA"
) {
  data.test <- GetAssayData(object = object,assay.type = assay.type,slot = "data")
  genes.use <- SetIfNull(x = genes.use, default = rownames(x = data.test))
  # data.use <- object@data
  if (print.bar) {
    iterate.fxn=pblapply
  } else {
    iterate.fxn <- lapply
  }
  p_val <- unlist(
    x = iterate.fxn(
      X = genes.use,
      FUN = function(x) {
        t.test(x = data.test[x, cells.1], y = data.test[x, cells.2])$p.value
      }
    )
  )
  to.return <- data.frame(p_val,row.names = genes.use)
  return(to.return)
}
