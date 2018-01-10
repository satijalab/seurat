globalVariables(
  names = c('alignment.index1', 'dups', 'cc_data2'),
  package = 'Seurat',
  add = TRUE
)
#' Align subspaces using dynamic time warping (DTW)
#'
#' Aligns subspaces across a given grouping variable.
#'
#' Following is a description for the two group case but this can be extended to
#' arbitrarily many groups which works by performing pairwise alignment to a
#' reference group (the largest group). First, we identify genes that are driving
#' variation in both datasets by looking at the correlation of gene expression
#' with each projection vector (e.g. CC1) in both datasets. For this we use the
#' biweight midcorrelation (bicor) and choose the top num.genes with the strongest
#' bicor to construct a 'metagene' for each dataset. We then scale each metagene
#' to match its 95\% reference range and linearly shift them by the minimum
#' difference between the two metagenes over the 10-90 quantile range. We then
#' map each cell in the smaller dataset to a cell in the larger dataset using
#' dynamic time warping (DTW) and apply the same map to the projection vectors (
#' CC vectors) to place both datasets on a common aligned scale. We apply this
#' procedue to each pair (group) of vectors individually for all specified in
#' dims.align. For a full description of the method, see Butler et al 2017.
#'
#' @param object Seurat object
#' @param reduction.type Reduction to align scores for. Default is "cca".
#' @param grouping.var Name of the grouping variable for which to align the scores
#' @param dims.align Dims to align, default is all
#' @param num.possible.genes Number of possible genes to search when choosing
#' genes for the metagene. Set to 2000 by default. Lowering will decrease runtime
#' but may result in metagenes constructed on fewer than num.genes genes.
#' @param num.genes Number of genes to use in construction of "metagene" (default
#' is 30).
#' @param show.plots Show debugging plots
#' @param verbose Displays progress and other output
#' @param ... Additional parameters to ScaleData
#'
#' @return Returns Seurat object with the dims aligned, stored in
#'  object@@dr$reduction.type.aligned
#'
#' @importFrom dtw dtw
#' @importFrom stats density
#' @importFrom pbapply pbapply
#' @importFrom graphics lines plot
#'
#' @export
#'
#' @examples
#' \dontrun{
#' pbmc_small
#' # Requires CCA to have previously been run
#' # As CCA requires two datasets, we will split our test object into two just for this example
#' pbmc1 <- SubsetData(pbmc_small,cells.use = pbmc_small@cell.names[1:40])
#' pbmc2 <- SubsetData(pbmc_small,cells.use = pbmc_small@cell.names[41:80])
#' pbmc1@meta.data$group <- "group1"
#' pbmc2@meta.data$group <- "group2"
#' pbmc_cca <- RunCCA(pbmc1,pbmc2)
#' pbmc_cca <- AlignSubspace(pbmc_cca,reduction.type = "cca", grouping.var = "group", dims.align = 1:2)
#' }
#'
AlignSubspace <- function(
  object,
  reduction.type = "cca",
  grouping.var,
  dims.align,
  num.possible.genes = 2000,
  num.genes = 30,
  show.plots = FALSE,
  verbose = TRUE,
  ...
) {
  parameters.to.store <- as.list(environment(), all = TRUE)[names(formals("AlignSubspace"))]
  object <- SetCalcParams(object = object,
                          calculation = paste0("AlignSubspace.", reduction.type),
                          ... = parameters.to.store)
  ident.orig <- object@ident
  object <- SetAllIdent(object = object, id = grouping.var)
  levels.split <- names(x = sort(x = table(object@ident), decreasing = T))
  num.groups <- length(levels.split)
  objects <- list()
  for (i in 1:num.groups){
    objects[[i]] <- SubsetData(object = object, ident.use = levels.split[i])
  }
  object@ident <- ident.orig
  cc.loadings <- list()
  scaled.data <- list()
  cc.embeds <- list()
  for (i in 1:num.groups) {
    if (verbose){
      cat(paste0("Rescaling group ", i, "\n"), file = stderr())
    }
    objects[[i]] <- ScaleData(object = objects[[i]], display.progress = verbose, ...)
    objects[[i]]@scale.data[is.na(x = objects[[i]]@scale.data)] <- 0
    objects[[i]] <- ProjectDim(
      object = objects[[i]],
      reduction.type = reduction.type,
      do.print = FALSE
    )
    cc.loadings[[i]] <- GetGeneLoadings(
      object = objects[[i]],
      reduction.type = reduction.type,
      use.full = TRUE
    )
    cc.embeds[[i]] <- GetCellEmbeddings(
      object = objects[[i]],
      reduction.type = reduction.type
    )
    scaled.data[[i]] <- objects[[i]]@scale.data
  }
  cc.embeds.all <- GetCellEmbeddings(object = object,
                                     reduction.type = reduction.type,
                                     dims.use = dims.align)
  colnames(cc.embeds.all) <- paste0("A", colnames(x = cc.embeds.all))
  cc.embeds.orig <- cc.embeds.all
  for (cc.use in dims.align) {
    for (g in 2:num.groups){
      if (verbose) {
        cat(paste0("Aligning dimension ", cc.use, "\n"), file = stderr())
      }
      genes.rank <- data.frame(
        rank(x = abs(x = cc.loadings[[1]][, cc.use])),
        rank(x = abs(x = cc.loadings[[g]][, cc.use])),
        cc.loadings[[1]][, cc.use],
        cc.loadings[[g]][, cc.use]
      )
      genes.rank$min <- apply(X = genes.rank[,1:2], MARGIN = 1, FUN = min)
      genes.rank <- genes.rank[order(genes.rank$min, decreasing = TRUE), ]
      genes.top <- rownames(x = genes.rank)[1:min(num.possible.genes, nrow(genes.rank))]
      bicors <- list()
      for (i in c(1, g)) {
        cc.vals <- cc.embeds[[i]][, cc.use]
        if(verbose) {
          bicors[[i]] <- pbsapply(
            X = genes.top,
            FUN = function(x) {
              return(BiweightMidcor(x = cc.vals, y = scaled.data[[i]][x, ]))
            }
          )
        } else {
          bicors[[i]] <- sapply(
            X = genes.top,
            FUN = function(x) {
              return(BiweightMidcor(x = cc.vals, y = scaled.data[[i]][x, ]))
            }
          )
        }
      }
      genes.rank <- data.frame(
        rank(x = abs(x = bicors[[1]])),
        rank(x = abs(x = bicors[[g]])),
        bicors[[1]],
        bicors[[g]]
      )
      genes.rank$min <- apply(X = abs(x = genes.rank[, 1:2]), MARGIN = 1, FUN = min)
      # genes must be correlated in same direction in both datasets
      genes.rank <- genes.rank[sign(genes.rank[,3]) == sign(genes.rank[,4]), ]
      genes.rank <- genes.rank[order(genes.rank$min, decreasing = TRUE), ]
      genes.use <- rownames(x = genes.rank)[1:min(num.genes, nrow(genes.rank))]
      if(length(genes.use) == 0) {
        stop("Can't align group ", g, " for dimension ", cc.use)
      }
      metagenes <- list()
      multvar.data <- list()
      for (i in c(1, g)) {
        scaled.use <- sweep(
          x = scaled.data[[i]][genes.use, ],
          MARGIN = 1,
          STATS = sign(x = genes.rank[genes.use, which(c(1, g) == i) + 2]),
          FUN = "*"
        )
        scaled.use <- scaled.use[, names(x = sort(x = cc.embeds[[i]][, cc.use]))]
        metagenes[[i]] <- (
          cc.loadings[[i]][genes.use, cc.use] %*% scaled.data[[i]][genes.use, ]
        )[1, colnames(x = scaled.use)]
      }
      mean.difference <- mean(x = ReferenceRange(x = metagenes[[g]])) -
        mean(x = ReferenceRange(x = metagenes[[1]]))
      metric.use <- "Euclidean"
      align.1 <- ReferenceRange(x = metagenes[[g]])
      align.2 <- ReferenceRange(x = metagenes[[1]])
      a1q <- sapply(
        X = seq(from = 0, to = 1, by = 0.001),
        FUN = function(x) {
          return(quantile(x = align.1, probs = x))
        }
      )
      a2q <- sapply(
        X = seq(from = 0, to = 1, by = 0.001),
        FUN = function(x) {
          quantile(x = align.2, probs = x)
        }
      )
      iqr <- (a1q - a2q)[100:900]
      iqr.x <- which.min(x = abs(x = iqr))
      iqrmin <- iqr[iqr.x]
      if (show.plots) {
        print(iqrmin)
      }
      align.2 <- align.2 + iqrmin
      alignment <- dtw(x = align.1,
                       y = align.2,
                       keep.internals = TRUE,
                       dist.method = metric.use)
      alignment.map <- data.frame(alignment$index1, alignment$index2)
      alignment.map$cc_data1 <- sort(cc.embeds[[g]][, cc.use])[alignment$index1]
      alignment.map$cc_data2 <- sort(cc.embeds[[1]][, cc.use])[alignment$index2]
      alignment.map.orig <- alignment.map
      alignment.map$dups <- duplicated(x = alignment.map$alignment.index1) |
        duplicated(x = alignment.map$alignment.index1, fromLast = TRUE)
      alignment.map %>% group_by(alignment.index1) %>% mutate(cc_data1_mapped = ifelse(dups, mean(cc_data2), cc_data2)) -> alignment.map
      alignment.map <- alignment.map[! duplicated(x = alignment.map$alignment.index1), ]
      cc.embeds.all[names(x = sort(x = cc.embeds[[g]][, cc.use])), cc.use] <- alignment.map$cc_data1_mapped
      if (show.plots) {
        par(mfrow = c(3, 2))
        plot(x = ReferenceRange(x = metagenes[[1]]), main = cc.use)
        plot(x = ReferenceRange(x = metagenes[[g]]))
        plot(
          x = ReferenceRange(x = metagenes[[1]])[(alignment.map.orig$alignment.index2)],
          pch = 16
        )
        points(
          x = ReferenceRange(metagenes[[g]])[(alignment.map.orig$alignment.index1)],
          col = "red",
          pch = 16,
          cex = 0.4
        )
        plot(x = density(x = alignment.map$cc_data1_mapped))
        lines(x = density(x = sort(x = cc.embeds[[1]][, cc.use])), col = "red")
        plot(x = alignment.map.orig$cc_data1)
        points(x = alignment.map.orig$cc_data2, col = "red")
      }
    }
  }
  new.type <- paste0(reduction.type, ".aligned")
  new.key <- paste0(
    "A",
    GetDimReduction(
      object = object,
      reduction.type = reduction.type,
      slot = "key"
    )
  )
  object <- SetDimReduction(
    object = object,
    reduction.type = new.type,
    slot = "cell.embeddings",
    new.data = scale(x = cc.embeds.all)
  )
  object <- SetDimReduction(
    object = object,
    reduction.type = new.type,
    slot = "key",
    new.data = new.key
  )
  return(object)
}

#' Calculate an alignment score
#'
#' Calculates an alignment score to determine how well aligned two (or more)
#' groups have been aligned. We first split the data into groups based on the
#' grouping.var provided and randomly downsample all groups to have as many cells
#' as in the smallest group. We then construct a nearest neighbor graph and ask
#' for each cell, how many of its neighbors have the same group identity as it
#' does. We then take the average over all cells, compare it to the expected
#' value for perfectly mixed neighborhoods, and scale it to range from 0 to 1.
#'
#' xbar is the average number of neighbors belonging to any cells' same group,
#' N is the number of groups in the given grouping.var, k is the number of
#' neighbors in the KNN graph.
#' \deqn{1 - \frac{\bar{x} - \frac{k}{N}}{k - \frac{k}{N}}}{1 - (xbar - k/N)/(k - k/N)}
#'
#' @param object Seurat object
#' @param reduction.use Stored dimensional reduction on which to build NN graph.
#' Usually going to be cca.aligned.
#' @param dims.use Dimensions to use in building the NN graph
#' @param grouping.var Grouping variable used in the alignment.
#' @param nn Number of neighbors to calculate in the NN graph construction
#'
#' @importFrom FNN get.knn
#' @export
#'
#' @examples
#' \dontrun{
#' pbmc_small
#' # As CCA requires two datasets, we will split our test object into two just for this example
#' pbmc1 <- SubsetData(pbmc_small, cells.use = pbmc_small@cell.names[1:40])
#' pbmc2 <- SubsetData(pbmc_small, cells.use = pbmc_small@cell.names[41:80])
#' pbmc1@meta.data$group <- "group1"
#' pbmc2@meta.data$group <- "group2"
#' pbmc_cca <- RunCCA(pbmc1,pbmc2)
#' pbmc_cca <- AlignSubspace(pbmc_cca, reduction.type = "cca",
#'                           grouping.var = "group", dims.align = 1:5)
#' CalcAlignmentMetric(pbmc_cca, reduction.use = "cca.aligned",
#'                     dims.use = 1:5, grouping.var =  "group")
#' }
#'
CalcAlignmentMetric <- function(object, reduction.use = "cca.aligned", dims.use,
                                grouping.var, nn){
  object <- SetAllIdent(object, grouping.var)
  object <- SubsetData(object, max.cells.per.ident = min(table(object@ident)))
  num.groups <- length(unique(object@ident))
  if(missing(nn)){
    nn <- ceiling(table(object@ident)[1] * 0.01 * num.groups)
  }
  dist.mat <- GetCellEmbeddings(object, reduction.type = reduction.use, dims.use = dims.use)
  object.fnn <- get.knn(dist.mat, k = nn)
  alignment.score <- sapply(1:length(object@cell.names), function(x) {
    cell.id <- object@ident[x]
    num.same.id <- length(which(object@ident[object.fnn$nn.index[x, ]] == cell.id))
  })
  alignment.score <- 1 - ((mean(alignment.score) - nn /num.groups) / (nn - nn/num.groups))
  return(unname(alignment.score))
}

#' Calculate the ratio of variance explained by ICA or PCA to CCA
#'
#' @param object Seurat object
#' @param reduction.type type of dimensional reduction to compare to CCA (pca,
#' pcafast, ica)
#' @param grouping.var variable to group by
#' @param dims.use Vector of dimensions to project onto (default is the 1:number
#' stored for cca)
#' @param verbose Display progress and other output
#'
#' @return Returns Seurat object with ratio of variance explained stored in
#' object@@meta.data$var.ratio
#' @export
#'
#' @examples
#' pbmc_small
#' # Requires CCA to have previously been run
#' # As CCA requires two datasets, we will split our test object into two just for this example
#' pbmc1 <- SubsetData(pbmc_small,cells.use = pbmc_small@cell.names[1:40])
#' pbmc2 <- SubsetData(pbmc_small,cells.use = pbmc_small@cell.names[41:80])
#' pbmc1@meta.data$group <- "group1"
#' pbmc2@meta.data$group <- "group2"
#' pbmc_cca <- RunCCA(pbmc1,pbmc2)
#' pbmc_cca <- CalcVarExpRatio(pbmc_cca,reduction.type = "pca", grouping.var = "group", dims.use = 1:5)
#'
CalcVarExpRatio <- function(
  object,
  reduction.type = "pca",
  grouping.var,
  dims.use,
  verbose = TRUE
) {
  if (missing(x = grouping.var)) {
    stop("Need to provide grouping variable")
  }
  if (missing(x = dims.use)) {
    dims.use <- 1:ncol(x = GetCellEmbeddings(object = object, reduction.type = "cca"))
  }
  parameters.to.store <- as.list(environment(), all = TRUE)[names(formals("CalcVarExpRatio"))]
  object <- SetCalcParams(object = object,
                          calculation = "CalcVarExpRatio",
                          ... = parameters.to.store)
  groups <- as.vector(x = unique(x = FetchData(
    object = object,
    vars.all = grouping.var
  )[, 1]))
  genes.use <- rownames(x = GetGeneLoadings(object = object, reduction.type = "cca"))
  var.ratio <- data.frame()
  for (group in groups) {
    if (verbose) {
      cat(paste("Calculating for", group, "\n"), file = stderr())
    }
    group.cells <- WhichCells(
      object = object,
      subset.name = grouping.var,
      accept.value = group
    )
    if (verbose) {
      cat(paste("\t Separating", group, "cells\n"), file = stderr())
    }
    group.object <- SubsetData(object = object, cells.use = group.cells)
    if (verbose) {
      cat("\t Running Dimensional Reduction \n", file = stderr())
    }
    ldp.cca <- CalcLDProj(
      object = group.object,
      reduction.type = "cca",
      dims.use = dims.use,
      genes.use = genes.use
    )
    group.object <- CalcProjectedVar(
      object = group.object,
      low.dim.data = ldp.cca,
      reduction.type = "cca",
      dims.use = dims.use,
      genes.use = genes.use
    )
    if (reduction.type == "pca") {
      temp.matrix <- PrepDR(group.object, genes.use = genes.use)
      group.object <- RunPCA(
        object = group.object,
        pc.genes = genes.use,
        do.print = FALSE,
        center = rowMeans(temp.matrix),
        pcs.compute = max(dims.use)
      )
      ldp.pca <- CalcLDProj(
        object = group.object,
        reduction.type = "pca",
        dims.use = dims.use,
        genes.use = genes.use
      )
      group.object <- CalcProjectedVar(
        object = group.object,
        low.dim.data = ldp.pca,
        reduction.type = "pca",
        dims.use = dims.use,
        genes.use = genes.use
      )
      group.var.ratio <- group.object@meta.data[, "cca.var", drop = FALSE] /
        group.object@meta.data[, "pca.var", drop = FALSE]
    } else if (reduction.type == "ica") {
      group.object <- RunICA(
        object = group.object,
        ic.genes = genes.use,
        print.results = FALSE,
        ics.compute = max(dims.use)
      )
      ldp.ica <- CalcLDProj(
        object = group.object,
        reduction.type = "ica",
        dims.use = dims.use,
        genes.use = genes.use
      )
      group.object <- CalcProjectedVar(
        object = group.object,
        low.dim.data = ldp.ica,
        reduction.type = "ica",
        dims.use = dims.use,
        genes.use = genes.use
      )
      group.var.ratio <- group.object@meta.data[, "cca.var", drop = FALSE] /
        group.object@meta.data[, "ica.var", drop = FALSE]
    } else {
      stop(paste("reduction.type", reduction.type, "not supported"))
    }
    var.ratio <- rbind(var.ratio, group.var.ratio)
  }
  var.ratio$cell.name <- rownames(x = var.ratio)
  eval(expr = parse(text = paste0(
    "object@meta.data$var.ratio.",
    reduction.type,
    "<- NULL"
  )))
  colnames(x = var.ratio) <- c(
    paste0("var.ratio.", reduction.type),
    "cell.name"
  )
  object <- AddMetaData(object, metadata = var.ratio)
  object@meta.data$cell.name <- NULL
  return(object)
}
