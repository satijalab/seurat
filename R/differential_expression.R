#' @include seurat.R
NULL

globalVariables(names = 'avg_diff', package = 'Seurat', add = TRUE)
#' Gene expression markers of identity classes
#'
#' Finds markers (differentially expressed genes) for identity classes
#'
#' @param object Seurat object
#' @param ident.1 Identity class to define markers for
#' @param ident.2 A second identity class for comparison. If NULL (default) -
#' use all other cells for comparison.
#' @param genes.use Genes to test. Default is to use all genes
#' @param thresh.use Limit testing to genes which show, on average, at least
#' X-fold difference (log-scale) between the two groups of cells. Default is 0.25
#' Increasing thresh.use speeds up the function, but can miss weaker signals.
#' @param test.use Denotes which test to use. Seurat currently implements
#' "bimod" (likelihood-ratio test for single cell gene expression, McDavid et
#' al., Bioinformatics, 2013, default), "roc" (standard AUC classifier), "t"
#' (Students t-test), and "tobit" (Tobit-test for differential gene expression,
#' as in Trapnell et al., Nature Biotech, 2014), 'poisson', and 'negbinom'.
#' The latter two options should only be used on UMI datasets, and assume an underlying
#' poisson or negative-binomial distribution
#' @param min.pct - only test genes that are detected in a minimum fraction of min.pct cells
#' in either of the two populations. Meant to speed up the function by not testing genes that are very infrequently expressed. Default is 0.1
#' @param min.diff.pct - only test genes that show a minimum difference in the fraction of detection between the two groups. Set to -Inf by default
#' @param only.pos Only return positive markers (FALSE by default)
#' @param print.bar Print a progress bar once expression testing begins (uses pbapply to do this)
#' @param max.cells.per.ident Down sample each identity class to a max number. Default is no downsampling. Not activated by default (set to Inf)
#' @param random.seed Random seed for downsampling
#' @param latent.vars Variables to test
#' @param min.cells Minimum number of cells expressing the gene in at least one of the two groups
#'
#' @return Matrix containing a ranked list of putative markers, and associated statistics (p-values, ROC score, etc.)
#'
#' @import pbapply
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
  thresh.use = 0.25,
  test.use = "bimod",
  min.pct = 0.1,
  min.diff.pct = -Inf,
  print.bar = TRUE,
  only.pos = FALSE,
  max.cells.per.ident = Inf,
  random.seed = 1,
  latent.vars = "nUMI",
  min.cells = 3
) {
  genes.use <- SetIfNull(x = genes.use, default = rownames(x = object@data))
  if (max.cells.per.ident < Inf) {
    object <- SubsetData(
      object = object,
      max.cells.per.ident = max.cells.per.ident,
      random.seed = random.seed
    )
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
      cells.2 <- object@cell.names
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
  if (length(cells.1) < min.cells) {
    stop(paste("Cell group 1 has fewer than", as.character(min.cells), "cells in identity class", ident.1))
  }
  if (length(cells.2) < min.cells) {
    stop(paste("Cell group 2 has fewer than", as.character(min.cells), " cells in identity class", ident.2))
  }

  # gene selection (based on percent expressed)
  thresh.min <- 0
  data.temp1 <- round(
    x = apply(
      X = object@data[genes.use, cells.1, drop = F],
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
      X = object@data[genes.use, cells.2, drop = F],
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
  if(length(genes.use) == 0) {
    stop("No genes pass min.pct threshold")
  }
  alpha.diff <- alpha.min - apply(X = data.alpha, MARGIN = 1, FUN = min)
  genes.use <- names(
    x = which(x = alpha.min > min.pct & alpha.diff > min.diff.pct)
  )
  if(length(genes.use) == 0) {
    stop("No genes pass min.diff.pct threshold")
  }

  #gene selection (based on average difference)
  data.1 <- apply(X = object@data[genes.use, cells.1, drop = F], MARGIN = 1, FUN = ExpMean)
  data.2 <- apply(X = object@data[genes.use, cells.2, drop = F], MARGIN = 1, FUN = ExpMean)
  total.diff <- (data.1 - data.2)
  genes.diff <- names(x = which(x = abs(x = total.diff) > thresh.use))
  genes.use <- intersect(x = genes.use, y = genes.diff)
  if(length(genes.use) == 0) {
    stop("No genes pass thresh.use threshold")
  }
  #perform DR
  if (test.use == "bimod") {
    to.return <- DiffExpTest(
      object = object,
      cells.1 = cells.1,
      cells.2 = cells.2,
      genes.use = genes.use,
      print.bar = print.bar
    )
  }
  if (test.use == "roc") {
    to.return <- MarkerTest(
      object = object,
      cells.1 = cells.1,
      cells.2 = cells.2,
      genes.use = genes.use,
      print.bar = print.bar
    )
  }
  if (test.use == "t") {
    to.return <- DiffTTest(
      object = object,
      cells.1 = cells.1,
      cells.2 = cells.2,
      genes.use = genes.use,
      print.bar = print.bar
    )
  }
  if (test.use == "tobit") {
    to.return <- TobitTest(
      object = object,
      cells.1 = cells.1,
      cells.2 = cells.2,
      genes.use = genes.use,
      print.bar = print.bar
    )
  }
  if (test.use == "negbinom") {
    to.return <- NegBinomDETest(
      object = object,
      cells.1 = cells.1,
      cells.2 = cells.2,
      genes.use = genes.use,
      latent.vars = latent.vars,
      print.bar = print.bar,
      min.cells = min.cells
    )
  }
  if (test.use == "poisson") {
    to.return <- PoissonDETest(
      object = object,
      cells.1 = cells.1,
      cells.2 = cells.2,
      genes.use = genes.use,
      latent.vars = latent.vars,
      print.bar = print.bar
      # min.cells # PoissonDETest doesn't have something for min.cells
    )
  }
  #return results
  to.return[, "avg_diff"] <- total.diff[rownames(x = to.return)]
  to.return <- cbind(to.return, data.alpha[rownames(x = to.return), ])
  if (test.use == "roc") {
    to.return <- to.return[order(-to.return$power, -to.return$avg_diff), ]
  } else {
    to.return <- to.return[order(to.return$p_val, -to.return$avg_diff), ]
  }
  if (only.pos) {
    to.return <- subset(x = to.return, subset = avg_diff > 0)
  }
  return(to.return)
}

globalVariables(
  names = c('myAUC', 'p_val', 'avg_diff'),
  package = 'Seurat',
  add = TRUE
)
#' Gene expression markers for all identity classes
#'
#' Finds markers (differentially expressed genes) for each of the identity classes in a dataset
#'
#' @param object Seurat object
#' @param genes.use Genes to test. Default is to all genes
#' @param thresh.use Limit testing to genes which show, on average, at least
#' X-fold difference (log-scale) between the two groups of cells.
#' Increasing thresh.use speeds up the function, but can miss weaker signals.
#' @param test.use Denotes which test to use. Seurat currently implements
#' "bimod" (likelihood-ratio test for single cell gene expression, McDavid et
#' al., Bioinformatics, 2013, default), "roc" (standard AUC classifier), "t"
#' (Students t-test), and "tobit" (Tobit-test for differential gene expression,
#' as in Trapnell et al., Nature Biotech, 2014), 'poisson', and 'negbinom'.
#' The latter two options should only be used on UMI datasets, and assume an underlying
#' poisson or negative-binomial distribution
#' @param min.pct - only test genes that are detected in a minimum fraction of min.pct cells
#' in either of the two populations. Meant to speed up the function by not testing genes that are very infrequently expression
#' @param min.diff.pct - only test genes that show a minimum difference in the fraction of detection between the two groups. Set to -Inf by default
#' @param only.pos Only return positive markers (FALSE by default)
#' @param print.bar Print a progress bar once expression testing begins (uses pbapply to do this)
#' @param max.cells.per.ident Down sample each identity class to a max number. Default is no downsampling.
#' @param random.seed Random seed for downsampling
#' @param return.thresh Only return markers that have a p-value < return.thresh, or a power > return.thresh (if the test is ROC)
#' @param do.print FALSE by default. If TRUE, outputs updates on progress.
#' @param min.cells Minimum number of cells expressing the gene in at least one of the two groups
#' @param latent.vars remove the effects of these variables
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
  thresh.use = 0.25,
  test.use = "bimod",
  min.pct = 0.1,
  min.diff.pct = 0.05,
  print.bar = TRUE,
  only.pos = FALSE,
  max.cells.per.ident = Inf,
  return.thresh = 1e-2,
  do.print = FALSE,
  random.seed = 1,
  min.cells = 3,
  latent.vars = "nUMI"
) {
  genes.use <- SetIfNull(x = genes.use, default = rownames(x = object@data))
  ident.use <- object@ident
  if ((test.use == "roc") && (return.thresh == 1e-2)) {
    return.thresh = 0.7
  }
  idents.all <- sort(x = unique(x = object@ident))
  genes.de <- list()
  if (max.cells.per.ident < Inf) {
    object <- SubsetData(
      object = object,
      max.cells.per.ident = max.cells.per.ident,
      random.seed = random.seed
    )
  }
  for (i in 1:length(x = idents.all)) {
    genes.de[[i]] <- tryCatch(
      {
        FindMarkers(object = object,
                    ident.1 = idents.all[i],
                    ident.2 = NULL,
                    genes.use = genes.use,
                    thresh.use = thresh.use,
                    test.use = test.use,
                    min.pct = min.pct,
                    min.diff.pct = min.diff.pct,
                    print.bar = print.bar,
                    min.cells = min.cells,
                    latent.vars = latent.vars
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
    gde <- genes.de[[i]]
    if (is.null(unlist(gde))) next
    if (nrow(x = gde) > 0) {
      if (test.use == "roc") {
        gde <- subset(
          x = gde,
          subset = (myAUC > return.thresh | myAUC < (1 - return.thresh))
        )
      } else {
        gde <- gde[order(gde$p_val, -gde$avg_diff), ]
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
  if (only.pos) {
    return(subset(x = gde.all, subset = avg_diff > 0))
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
  thresh.use = 0.25,
  test.use = "bimod",
  ...
) {
  genes.use <- SetIfNull(x = genes.use, default = rownames(x = object@data))
  tree <- SetIfNull(x = tree.use, default = object@cluster.tree[[1]])
  ident.order <- tree$tip.label
  nodes.1 <- ident.order[GetLeftDescendants(tree = tree, node = node)]
  nodes.2 <- ident.order[GetRightDescendants(tree = tree, node = node)]
  #print(nodes.1)
  #print(nodes.2)
  to.return <- FindMarkers(
    object = object,
    ident.1 = nodes.1,
    ident.2 = nodes.2,
    genes.use = genes.use,
    thresh.use = thresh.use,
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
#' @param thresh.use Limit testing to genes which show, on average, at least
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
#' @param min.cells Minimum number of cells expressing the gene in at least one of the two groups
#'
#' @return Returns a dataframe with a ranked list of putative markers for each node and associated statistics
#'
#' @importFrom ape drop.tip
#'
#' @export
#'
#' @examples
#' FindAllMarkersNode(pbmc_small)
#'
FindAllMarkersNode <- function(
  object,
  node = NULL,
  genes.use = NULL,
  thresh.use = 0.25,
  test.use = "bimod",
  min.pct = 0.1,
  min.diff.pct = 0.05,
  print.bar = TRUE,
  only.pos = FALSE,
  max.cells.per.ident = Inf,
  return.thresh = 1e-2,
  do.print = FALSE,
  random.seed = 1,
  min.cells = 3
) {
  if(length(object@cluster.tree) == 0){
    stop("Tree hasn't been built yet. Run BuildClusterTree to build.")
  }
  genes.use <- SetIfNull(x = genes.use, default = rownames(object@data))
  node <- SetIfNull(x = node, default = object@cluster.tree[[1]]$edge[1, 1])
  ident.use <- object@ident
  tree.use <- object@cluster.tree[[1]]
  descendants <- DFT(tree = tree.use, node = node, path = NULL, include.children = TRUE)
  all.children <- sort(x = tree.use$edge[,2][!tree.use$edge[,2] %in% tree.use$edge[,1]])
  descendants <- MapVals(v = descendants, from = all.children, to = tree.use$tip.label)
  drop.children <- setdiff(tree.use$tip.label, descendants)
  keep.children <- setdiff(tree.use$tip.label, drop.children)
  orig.nodes <- c(node, as.numeric(setdiff(descendants, keep.children)))
  tree.use <- drop.tip(tree.use, drop.children)
  new.nodes <- unique(tree.use$edge[,1])
  if ((test.use == 'roc') && (return.thresh==1e-2)) {
    return.thresh <- 0.7
  }
  genes.de <- list()
  for (i in ((tree.use$Nnode+2):max(tree.use$edge))) {
    genes.de[[i]] <- FindMarkersNode(
      object = object,
      node = i,
      tree.use = tree.use,
      genes.use = genes.use,
      thresh.use = thresh.use,
      test.use = test.use,
      min.pct = min.pct,
      min.diff.pct = min.diff.pct,
      print.bar = print.bar,
      only.pos = only.pos,
      max.cells.per.ident = max.cells.per.ident,
      random.seed = random.seed,
      min.cells = min.cells
    )
    if (do.print) {
      print(paste("Calculating node", i))
    }
  }
  gde.all <- data.frame()
  for (i in ((tree.use$Nnode+2):max(tree.use$edge))) {
    gde <- genes.de[[i]]
    if (is.null(x = unlist(gde))) {
      next
    }
    if (nrow(x = gde) > 0) {
      if (test.use == 'roc') {
        gde <- subset(
          x = gde,
          subset = (myAUC > return.thresh | myAUC < (1 - return.thresh))
        )
      }
      if ( (test.use == 'bimod') || (test.use == 't')) {
        gde <- gde[order(gde$p_val,-gde$avg_diff), ]
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
  gde.all$cluster <- MapVals(v = gde.all$cluster,
                             from = new.nodes,
                             to = orig.nodes
  )
  return(gde.all)
}

#' Finds markers that are conserved between the two groups
#'
#' @param object Seurat object
#' @param ident.1 Identity class to define markers for
#' @param ident.2 A second identity class for comparison. If NULL (default) - use all other cells
#' for comparison.
#' @param grouping.var grouping variable
#' @param \dots parameters to pass to FindMarkers
#'
#' @return Matrix containing a ranked list of putative conserved markers, and associated statistics
#' (p-values within each group and a combined p-value (fisher_pval), percentage of cells expressing
#' the marker, average differences)
#'
#' @export
#'
#' @examples
#' pbmc_small
#' # Create a simulated grouping variable
#' pbmc_small@meta.data$groups <- sample(
#'   x = c("g1", "g2"),
#'   size = length(x = pbmc_small@cell.names),
#'   replace = TRUE
#' )
#' FindConservedMarkers(pbmc_small, ident.1 = 1, ident.2 = 2, grouping.var = "groups")
#'
FindConservedMarkers <- function(
  object,
  ident.1,
  ident.2 = NULL,
  grouping.var,
  ...
) {
  object.var <- FetchData(object = object, vars.all = grouping.var)
  object <- SetIdent(
    object = object,
    cells.use = object@cell.names,
    ident.use = paste(object@ident, object.var[, 1], sep = "_")
  )
  levels.split <- names(x = sort(x = table(object.var[, 1])))
  if (length(x = levels.split) != 2) {
    stop(
      paste0(
        "There are not two options for ",
        grouping.var,
        ". \n Current groups include: ",
        paste(levels.split, collapse = ", ")
      )
    )
  }
  cells <- list()
  for (i in 1:2) {
    cells[[i]] <- rownames(
      x = object.var[object.var[, 1] == levels.split[i], , drop = FALSE]
    )
  }
  marker.test <- list()
  # do marker tests
  for (i in 1:2) {
    level.use <- levels.split[i]
    ident.use.1 <- paste(ident.1, level.use, sep = "_")
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
    marker.test[[i]] <- FindMarkers(
      object = object,
      ident.1 = ident.use.1,
      ident.2 = ident.use.2,
      ...
    )
  }
  genes.conserved <- intersect(
    x = rownames(x = marker.test[[1]]),
    y = rownames(x = marker.test[[2]])
  )
  markers.conserved <- list()
  for (i in 1:2) {
    markers.conserved[[i]] <- marker.test[[i]][genes.conserved, ]
    colnames(x = markers.conserved[[i]]) <- paste(
      levels.split[i],
      colnames(x = markers.conserved[[i]]),
      sep="_"
    )
  }
  markers.combined <- cbind(markers.conserved[[1]], markers.conserved[[2]])
  pval.codes <- paste(levels.split, "p_val", sep = "_")
  markers.combined$max_pval <- apply(
    X = markers.combined[, pval.codes],
    MARGIN = 1,
    FUN = max
  )
  markers.combined$fisher_pval <- apply(
    X = markers.combined[, pval.codes],
    MARGIN = 1,
    FUN = FisherIntegrate
  )
  markers.combined <- markers.combined[order(markers.combined$fisher_pval), ]
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
  genes.use = NULL,
  print.bar = TRUE
) {
  genes.use <- SetIfNull(x = genes.use, default = rownames(x = object@data))
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
            x = as.numeric(x = object@data[x, cells.1]),
            y = as.numeric(x = object@data[x, cells.2])
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
  min.cells = 3
) {
  genes.use <- SetIfNull(x = genes.use, default = rownames(x = object@data))
  # check that the gene made it through the any filtering that was done
  genes.use <- genes.use[genes.use %in% rownames(x = object@data)]
  my.latent <- FetchData(
    object = object,
    vars.all = latent.vars,
    cells.use = c(cells.1, cells.2),
    use.raw = TRUE
  )
  to.test.data <- object@raw.data[genes.use, c(cells.1, cells.2)]
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
  min.cells = 3
) {
  genes.use <- SetIfNull(x = genes.use, default = rownames(x = object@data))
  # check that the gene made it through the any filtering that was done
  genes.use <- genes.use[genes.use %in% rownames(x = object@data)]
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
  to.test.data <- object@raw.data[genes.use, c(cells.1, cells.2), drop = FALSE]
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
  print.bar = TRUE
) {
  genes.use <- SetIfNull(x = genes.use, default = rownames(x = object@data))
  # check that the gene made it through the any filtering that was done
  genes.use <- genes.use[genes.use %in% rownames(x = object@data)]
  my.latent <- FetchData(
    object = object,
    vars.all = latent.vars,
    cells.use = c(cells.1, cells.2),
    use.raw = TRUE
  )
  to.test.data <- object@raw.data[genes.use, c(cells.1, cells.2)]
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
  print.bar = TRUE
) {
  genes.use <- SetIfNull(x = genes.use, default = rownames(x = object@data))
  #print(genes.diff)
  to.return <- TobitDiffExpTest(
    data1 = object@data[, cells.1],
    data2 = object@data[, cells.2],
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
  print.bar = TRUE
) {
  genes.use <- SetIfNull(x = genes.use, default = rownames(x = object@data))
  to.return <- AUCMarkerTest(
    data1 = object@data[, cells.1],
    data2 = object@data[, cells.2],
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
  print.bar = TRUE
) {
  genes.use <- SetIfNull(x = genes.use, default = rownames(x = object@data))
  data.use <- object@data
  if (print.bar) {
    iterate.fxn=pblapply
  } else {
    iterate.fxn <- lapply
  }
  p_val <- unlist(
    x = iterate.fxn(
      X = genes.use,
      FUN = function(x) {
        t.test(x = object@data[x, cells.1], y = object@data[x, cells.2])$p.value
      }
    )
  )
  to.return <- data.frame(p_val,row.names = genes.use)
  return(to.return)
}
