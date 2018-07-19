#' Find gene terms from Enrichr
#'
#' Fing gene terms from Enrichr and return the XML information
#'
#' @param QueryGene A gene to query on Enrichr
#'
#' @export
#' @importFrom httr GET status_code content
#'
#' @return An XML document with information on the queried gene
#'
FindGeneTerms <- function(QueryGene = NULL) {
  if (is.null(x = QueryGene)) {
    stop("Missing query gene")
  }
  path.use <- "Enrichr/genemap"
  api.get <- GET(
    url = "http://amp.pharm.mssm.edu/",
    path = path.use,
    query = list(gene = QueryGene)
  )
  api.status <- status_code(x = api.get)
  if (api.status != 200) {
    stop("Error searching for terms")
  }
  api.data <- content(x = api.get)
  return (api.data)
}

# Run t-distributed Stochastic Neighbor Embedding
#
# Run t-SNE dimensionality reduction on selected features. Has the option of running in a reduced
# dimensional space (i.e. spectral tSNE, recommended), or running based on a set of genes
#
# @param object Seurat object
# @param cells.use Which cells to analyze (default, all cells)
# @param dims.use Which dimensions to use as input features
# @param k.seed Random seed for the t-SNE
# @param do.fast If TRUE, uses the Barnes-hut implementation, which runs
# faster, but is less flexible
# @param add.iter If an existing tSNE has already been computed, uses the
# current tSNE to seed the algorithm and then adds additional iterations on top of this
# @param genes.use If set, run the tSNE on this subset of genes
# (instead of running on a set of reduced dimensions). Not set (NULL) by default
# @param reduction.use Which dimensional reduction (PCA or ICA) to use for the tSNE. Default is PCA
# @param dim_embed The dimensional space of the resulting tSNE embedding (default is 2).
# For example, set to 3 for a 3d tSNE
# @param q.use Quantile to use
# @param max.dim Max dimension to keep from diffusion calculation
# @param scale.clip Max/min value for scaled data. Default is 3
# @param ... Additional arguments to the tSNE call. Most commonly used is
# perplexity (expected number of neighbors default is 30)
#
# @return Returns a Seurat object with a tSNE embedding in object@@tsne_rot
#
# @import Rtsne
# @import tsne
#
# Not currently supported
#
AddTSNE <- function(
  object,
  cells.use = NULL,
  pcs.use = 1:10,
  do.plot = TRUE,
  k.seed = 1,
  add.iter = 1000,
  ...
) {
  cells.use <- SetIfNull(x = cells.use, default = colnames(x = object@data))
  data.use <- object@pca.rot[cells.use, pcs.use]
  #data.dist=as.dist(mahalanobis.dist(data.use))
  set.seed(seed = k.seed)
  data.tsne <- data.frame(
    tsne(
      X = data.use,
      initial_config = as.matrix(x = object@tsne.rot[cells.use,]),
      max_iter = add.iter,
      ...
    )
  )
  colnames(x = data.tsne) <- paste0("tSNE_", 1:ncol(data.tsne))
  rownames(x = data.tsne) <- cells.use
  object@tsne.rot <- data.tsne
  return(object)
}

#calculate true positives, used in AUC
calcTP <- function(cutoff, data, score, real, nTP) {
  return(length(x = which(x = (data[, score] > cutoff) & (data[, real] > 0))) / nTP)
}

#calculate false positives, used in AUC
calcFP <- function(cutoff, data, score, real, nFP) {
  return(length(x = which(x = (data[, score] > cutoff) & (data[, real] == 0))) / nFP)
}

#i do not believe we use this function, but leaving it in to be safe
auc <- function(data, score, real, n = 20) {
  totalPos <- length(x = which(x = data[, real] == 1))
  totalNeg <- length(x = which(x = data[, real] == 0))
  scores <- data[, score]
  data$myScore <- (scores + min(scores)) / (max(scores) + min(scores))
  tp <- unlist(x = lapply(
    X = seq(from = -0.0001, to = 0.9999, by = 1 / n),
    FUN = calcTP,
    data = data,
    score = "myScore",
    real = real,
    nTP = totalPos
  ))
  fp <- unlist(x = lapply(
    X = seq(from = -0.0001, to = 0.9999, by = 1 / n),
    FUN = calcFP,
    data = data,
    score = "myScore",
    real = real,
    nFP = totalNeg
  ))
  plot(x = c(fp, 1), y = c(tp, 1), xlim = c(0, 1), ylim = c(0, 1))
  x1 <- c(1, fp)
  x2 <- c(1, tp)
  print(
    sum(diff(x = rev(x = x2)) * diff(x = rev(x = x1))) /
      2 + sum(diff(x = rev(x = x1)) * rev(x = x2[-1]))
  )
  return(list(c(1, fp), c(1, tp)))
}

#useful with FindAllMarkers, not ready for main package yet
returnTopX <- function(data, group.by, n.return, col.return = NA) {
  to.ret <- c()
  levels.use=unique(group.by); if (is.factor(group.by)) levels.use=levels(group.by)
  if (!is.na(col.return)) return(unlist(lapply(levels.use, function(x) head(data[group.by==x,col.return],n.return)))) else {
    return(unlist(lapply(levels.use, function(x) head(rownames(data[group.by==x,])))))
  }
}

#i like this, but not used too much yet
genes.ca.range <- function(object, my.min, my.max) {
  ca <- AverageDetectionRate(object = object)
  ca.min <- apply(X = ca, MARGIN = 1, FUN = min)
  ca.max <- apply(X = ca, MARGIN = 1, FUN = max)
  genes.1 <- names(x = ca.min[ca.min < my.max])
  genes.2 <- names(x = ca.max[ca.max > my.min])
  return(intersect(x = genes.1, y = genes.2))
}

# # Not currently supported, but a cool function for QC
# CalcNoiseModels <- function(
#   object,
#   cell.ids = NULL,
#   trusted.genes = NULL,
#   n.bin = 20,
#   drop.expr = 1
# ) {
#   object@drop.expr <- drop.expr
#   cell.ids <- SetIfNull(x = cell.ids, default = 1:ncol(x = object@data))
#   trusted.genes <- SetIfNull(x = trusted.genes, default = rownames(x = object@data))
#   trusted.genes <- trusted.genes[trusted.genes %in% rownames(x = object@data)]
#   object@trusted.genes <- trusted.genes
#   data <- object@data[trusted.genes, ]
#   idents <- data.frame(data[, 1])
#   code_humpAvg <- apply(X = data, MARGIN = 1, FUN = MeanGreaterThan, min = object@drop.expr)
#   code_humpAvg[code_humpAvg > 9] <- 9
#   code_humpAvg[is.na(x = code_humpAvg)] <- 0
#   idents$code_humpAvg <- code_humpAvg
#   data[data > object@drop.expr] <- 1
#   data[data < object@drop.expr] <- 0
#   data$bin <- cut(x = code_humpAvg, breaks = n.bin)
#   data$avg <- code_humpAvg
#   rownames(x = idents) <- rownames(x = data)
#   my.coefs <- data.frame(t(x = pbsapply(
#     X = colnames(x = data[1:(ncol(x = data) - 2)]),
#     FUN = getAB,
#     data = data,
#     data2 = idents,
#     status = "code",
#     code2 = "humpAvg",
#     hasBin = TRUE,
#     doPlot = FALSE
#   )))
#   colnames(x = my.coefs) <- c("a", "b")
#   object@drop.coefs <- my.coefs
#   return(object)
# }

# Visualize expression/dropout curve
#
# Plot the probability of detection vs average expression of a gene.
#
# Assumes that this 'noise' model has been precomputed with CalcNoiseModels
#
# @param object Seurat object
# @param cell.ids Cells to use
# @param col.use Color code or name
# @param lwd.use Line width for curve
# @param do.new Create a new plot (default) or add to existing
# @param x.lim Maximum value for X axis
# @param \dots Additional arguments to pass to lines function
# @return Returns no value, displays a plot
#
# PlotNoiseModel <- function(
#   object,
#   cell.ids = c(1, 2),
#   col.use = 'black',
#   lwd.use = 2,
#   do.new = TRUE,
#   x.lim = 10,
#   ...
# ) {
#   cell.coefs <- object@drop.coefs[cell.ids,]
#   if (do.new) {
#     plot(
#       x = 1,
#       y = 1,
#       pch = 16,
#       type = 'n',
#       xlab = 'Average expression',
#       ylab = 'Probability of detection',
#       xlim = c(0, x.lim),
#       ylim =c(0, 1)
#     )
#   }
#   unlist(
#     x = lapply(
#       X = 1:length(x = cell.ids),
#       FUN = function(y) {
#         x.vals <- seq(from = 0, to = x.lim, by = 0.05)
#         y.vals <- unlist(x = lapply(
#           X = x.vals,
#           FUN = calc.drop.prob,
#           a = cell.coefs[y, 1],
#           b = cell.coefs[y, 2])
#         )
#         lines(
#           x = x.vals,
#           y = y.vals,
#           lwd = lwd.use,
#           col = col.use,
#           ...
#         )
#       }
#     )
#   )
# }

# Deleted, but used in regression.sig (which was kept)...
vsubc <- function(data,code) {
  return(data[grep(pattern = code, x = names(x = data))])
}

regression.sig <- function(x, score, data, latent, code = "rsem") {
  if (var(x = as.numeric(x = SubsetColumn(data = data, code = code)[x, ])) == 0) {
    return(0)
  }
  latent <- latent[grep(pattern = code, x = names(x = data))]
  data <- rbind(SubsetColumn(data = data, code = code), vsubc(data = score, code = code))
  rownames(x = data)[nrow(x = data)] <- "score"
  data2 <- data[c(x, "score"), ]
  rownames(x = data2)[1] <- "fac"
  if (length(x = unique(x = latent)) > 1) {
    mylm <- lm(formula = score ~ fac + latent, data = data.frame(t(x = data2)))
  } else {
    mylm <- lm(formula = score ~ fac, data = data.frame(t(x = data2)))
  }
  return(coef(object = summary(object = mylm))["fac", 3])
}

# Documentation
# @export
#
RemovePC <- function(object, pcs.remove, use.full = FALSE, ...) {
  data.old <- object@data
  pcs.use <- setdiff(x = 1:ncol(x = object@pca.obj[[1]]$rotation), y = pcs.remove)
  if (use.full) {
    data.x <- as.matrix(
      x = object@pca.x.full[, intersect(
        x = pcs.use,
        y = 1:ncol(x = object@pca.x.full)
      )]
    )
  } else {
    data.x <- as.matrix(x = object@pca.obj[[1]]$x[, pcs.use])
  }
  data.1 <- data.x %*% t(x = as.matrix(x = object@pca.obj[[1]]$rotation[, pcs.use]))
  data.2 <- sweep(
    x = data.1,
    MARGIN = 2,
    STATS = colMeans(x = object@scale.data),
    FUN = "+"
  )
  data.3 <- sweep(
    x = data.2,
    MARGIN = 1,
    STATS = apply(X = object@data[rownames(x = data.2), ], MARGIN = 1, FUN = sd),
    FUN = "*"
  )
  data.3 <- sweep(
    X = data.3,
    MARGIN = 1,
    STATS = apply(X = object@data[rownames(x = data.2), ], MARGIN = 1, FUN = mean),
    FUN = "+"
  )
  object@scale.data <- (data.2)
  data.old <- data.old[rownames(x = data.3), ]
  data.4 <- data.3
  data.4[data.old == 0] <- 0
  data.4[data.4 < 0] <- 0
  object@data[rownames(x = data.4), ] <- data.frame(data.4)
  return(object)
}

# Documentation
# @export
#
CalinskiPlot <- function(
  object,
  pcs.use,
  pval.cut = 0.1,
  gene.max = 15,
  col.max = 25,
  use.full = TRUE
) {
  if (length(x = pcs.use) == 1) {
    pvals.min <- object@jackStraw.empP.full[, pcs.use]
  } else if (length(x = pcs.use) > 1) {
    pvals.min <- apply(
      X = object@jackStraw.empP.full[, pcs.use],
      MARGIN = 1,
      FUN = min
    )
  }
  names(x = pvals.min) <- rownames(x = object@jackStraw.empP.full)
  genes.use <- names(x = pvals.min)[pvals.min < pval.cut]
  genes.use <- genes.use[genes.use %in% rownames(do.NULL = object@scale.data)]
  par(mfrow = c(1, 2))
  mydata <- object@scale.data[genes.use, ]
  wss <- (nrow(x = mydata) - 1) * sum(apply(X = mydata, MARGIN = 2, FUN = var))
  for (i in 1:gene.max) {
    wss[i] <- sum(kmeans(x = mydata, centers=i)$withinss)
  }
  plot(
    x = 1:gene.max,
    y = wss,
    type = "b",
    xlab= "Number of Clusters for Genes",
    ylab = "Within groups sum of squares"
  )
  mydata <- t(x = object@scale.data[genes.use, ])
  wss <- (nrow(x = mydata) - 1) * sum(apply(X = mydata, MARGIN = 2, FUN = var))
  for (i in 1:col.max) {
    wss[i] <- sum(kmeans(x = mydata, centers = i)$withinss)
  }
  plot(
    x = 1:col.max,
    y = wss,
    type = "b",
    xlab = "Number of Clusters for Cells",
    ylab = "Within groups sum of squares"
  )
  ResetPar()
  return(object)
}

# Documentation
# #' @export
# #' @importFrom stats cor kmeans
# #' @importFrom NMF aheatmap
#
# CellCorMatrix <- function(
#   object,
#   cor.genes = NULL,
#   cell.inds = NULL,
#   do.k = FALSE,
#   k.seed = 1,
#   k.num = 4,
#   vis.low = (-1),
#   vis.high = 1,
#   vis.one = 0.8,
#   pcs.use = 1:3,
#   col.use = pyCols
# ) {
#   cor.genes <- SetIfNull(x = cor.genes, default = object@var.genes)
#   cell.inds <- SetIfNull(x = cell.inds, default = colnames(x = object@data))
#   cor.genes <- cor.genes[cor.genes %in% rownames(x = object@data)]
#   data.cor <- object@scale.data[cor.genes, cell.inds]
#   cor.matrix <- cor(x = data.cor)
#   set.seed(seed = k.seed)
#   kmeans.cor <- kmeans(x = cor.matrix, centers = k.num)
#   if (do.k) {
#     cor.matrix <- cor.matrix[order(kmeans.cor$cluster), order(kmeans.cor$cluster)]
#   }
#   kmeans.names <- rownames(x = cor.matrix)
#   row.annot <- data.frame(
#     cbind(
#       kmeans.cor$cluster[kmeans.names],
#       object@pca.rot[kmeans.names, pcs.use]
#     )
#   )
#   colnames(x = row.annot) <- c("K", paste0("PC", pcs.use))
#   cor.matrix[cor.matrix == 1] <- vis.one
#   cor.matrix <- MinMax(data = cor.matrix, min = vis.low, max = vis.high)
#   object@kmeans.cell <- list(kmeans.cor)
#   if (do.k) {
#     aheatmap(
#       x = cor.matrix,
#       col = col.use,
#       Rowv = NA,
#       Colv = NA,
#       annRow = row.annot
#     )
#   } else {
#     heatmap.2(
#       x = cor.matrix,
#       trace = "none",
#       Rowv = NA,
#       Colv = NA,
#       col = pyCols
#     )
#   }
#   return(object)
# }

# Documentation
# #' @export
#
# GeneCorMatrix <- function(
#   object,
#   cor.genes = NULL,
#   cell.inds = NULL,
#   do.k = FALSE,
#   k.seed = 1,
#   k.num = 4,
#   vis.low = (-1),
#   vis.high = 1,
#   vis.one = 0.8,
#   pcs.use = 1:3,
#   col.use = pyCols
# ) {
#   cor.genes <- SetIfNull(x = cor.genes, default = object@var.genes)
#   cell.inds <- SetIfNull(x = cell.inds, default = colnames(x = object@data))
#   cor.genes <- cor.genes[cor.genes %in% rownames(x = object@data)]
#   data.cor <- object@data[cor.genes, cell.inds]
#   cor.matrix <- cor(x = t(x = data.cor))
#   set.seed(seed = k.seed)
#   kmeans.cor <- kmeans(x = cor.matrix, centers = k.num)
#   cor.matrix <- cor.matrix[order(kmeans.cor$cluster), order(kmeans.cor$cluster)]
#   kmeans.names <- rownames(x = cor.matrix)
#   row.annot <- data.frame(
#     cbind(
#       kmeans.cor$cluster[kmeans.names],
#       object@pca.x[kmeans.names, pcs.use]
#     )
#   )
#   colnames(x = row.annot) <- c("K", paste0("PC", pcs.use))
#   cor.matrix[cor.matrix == 1] <- vis.one
#   cor.matrix <- MinMax(data = cor.matrix, min = vis.low, max = vis.high)
#   object@kmeans.gene <- list(kmeans.cor)
#   if (do.k) {
#     aheatmap(
#       x = cor.matrix,
#       col = col.use,
#       Rowv = NA,
#       Colv = NA,
#       annRow = row.annot
#     )
#   } else {
#     aheatmap(
#       x = cor.matrix,
#       col = col.use,
#       annRow = row.annot
#     )
#   }
#   return(object)
# }


#   Documentation
ProjectSamples <- function(object, new.samples) {
  genes.use <- rownames(x = object@data)
  genes.pca <- rownames(x = object@pca.x)
  data.project <- object@scale.data[genes.pca, ]
  data.project[is.na(x = data.project)] <- 0
  new.rot = t(x = data.project) %*% as.matrix(x = object@pca.x)
  scale.vec = apply(
    X = new.rot,
    MARGIN = 2,
    FUN = function(x) {
      return(sqrt(x = sum(x ^ 2)))
    }
  )
  new.rot.scale <- scale(x = new.rot, center = FALSE, scale = scale.vec)
  object@pca.rot <- as.data.frame(x = new.rot.scale)
  return(object)
}

# Documentation
#Cool, but not supported right now
SpatialDe <- function(object, marker.cells, genes.use = NULL) {
  embed.map <- object@tsne.rot
  mult.use <- 2
  mult.use.far <- 10
  if ((mult.use.far * length(x = marker.cells)) > nrow(x = embed.map)) {
    mult.use.far <- 1
    mult.use <- 1
  }
  genes.use <- SetIfNull(x = genes.use, default = object@var.genes)
  marker.pos <- apply(X = embed.map[marker.cells, ], MARGIN = 2, FUN = mean)
  embed.map <- rbind(embed.map, marker.pos)
  rownames(x = embed.map)[nrow(x = embed.map)] <- "marker"
  embed.dist <- sort(x = as.matrix(x = dist(x = (embed.map)))["marker", ])
  embed.diff <- names(x = embed.dist[! (names(x = embed.dist) %in% marker.cells)][1:(mult.use * length(x = marker.cells))][-1])
  embed.diff.far <- names(x = embed.dist[! (names(x = embed.dist) %in% marker.cells)][1:(mult.use.far * length(x = marker.cells))][-1])
  diff.genes <- rownames(
    x = subset(
      x = DiffExpTest(
        object = object,
        cells.1 = marker.cells,
        cells.2 = embed.diff,
        genes.use = genes.use
      ),
      subset = p_val < (1e-5)
    )
  )
  diff.genes <- subset(
    x = DiffExpTest(
      object = object,
      cells.1 = marker.cells,
      cells.2 = embed.diff,
      genes.use = diff.genes
    ),
    subset = p_val<(1e-10)
  )
  return(diff.genes)
}

#' Identify potential genes associated with batch effects
#'
#' Test for genes whose expression value is strongly predictive of batch (based
#' on ROC classification). Important note: Assumes that the 'batch' of each
#' cell is assigned to the cell's identity class (will be improved in a future
#' release)
#'
#' @param object Seurat object
#' @param idents.use Batch names to test
#' @param genes.use Gene list to test
#' @param auc.cutoff Minimum AUC needed to qualify as a 'batch gene'
#' @param thresh.use Limit testing to genes which show, on average, at least
#' X-fold difference (log-scale) in any one batch
#'
#' @return Returns a list of genes that are strongly correlated with batch.
#'
#' @export
#'
BatchGene <- function(
  object,
  idents.use,
  genes.use = NULL,
  auc.cutoff = 0.6,
  thresh.use = 0
) {
  batch.genes <- c()
  genes.use <- SetIfNull(x = genes.use, default = rownames(x = object@data))
  for (ident in idents.use) {
    cells.1 <- names(x = object@ident)[object@ident == ident]
    cells.2 <- names(x = object@ident)[object@ident != ident]
    if ((length(x = cells.1) < 5) | (length(x = cells.2) < 5)) {
      break
    }
    markers.ident <- MarkerTest(
      object = object,
      cells.1 = cells.1,
      cells.2 = cells.2,
      genes.use = genes.use,
      print.bar = thresh.use
    )
    batch.genes <- unique(
      x = c(
        batch.genes,
        rownames(x = subset(x = markers.ident, subset = myAUC > auc.cutoff))
      )
    )
  }
  return(batch.genes)
}



# Documentation
#multicore version of jackstraw
#DOES NOT WORK WITH WINDOWS
# @export
#
JackStrawMC <- function(
  object,
  num.pc = 30,
  num.replicate = 100,
  prop.freq = 0.01,
  do.print = FALSE,
  num.cores = 8
) {
  pc.genes <- rownames(x = object@pca.x)
  if (length(x = pc.genes) < 200) {
    prop.freq <- max(prop.freq, 0.015)
  }
  md.x <- as.matrix(x = object@pca.x)
  md.rot <- as.matrix(x = object@pca.rot)
  if (do.print) {
    fake.pcVals.raw <- mclapply(
      X = 1:num.replicate,
      FUN = function(x) {
        print(x)
        return(JackRandom(
          scaled.data = object@scale.data[pc.genes, ],
          prop.use = prop.freq,
          r1.use = 1,
          r2.use = num.pc,
          seed.use = x
        ))
      },
      mc.cores = num.cores
    )
  } else {
    fake.pcVals.raw <- mclapply(
      X = 1:num.replicate,
      FUN = function(x) {
        return(JackRandom(
          scaled.data = object@scale.data[pc.genes, ],
          prop.use = prop.freq,
          r1.use = 1,
          r2.use = num.pc,
          seed.use=x
        ))
      }, mc.cores = num.cores
    )
  }
  fake.pcVals <- simplify2array(
    x = mclapply(
      X = 1:num.pc,
      FUN = function(x) {
        return(as.numeric(x = unlist(x = lapply(
          X = 1:num.replicate,
          FUN = function(y) {
            return(fake.pcVals.raw[[y]][, x])
          }
        ))))
      },
      mc.cores = num.cores
    )
  )
  object@jackStraw.fakePC <- data.frame(fake.pcVals)
  object@jackStraw.empP <- data.frame(
    simplify2array(
      x = mclapply(
        X = 1:num.pc,
        FUN = function(x) {
          return(unlist(x = lapply(
            X = abs(md.x[, x]),
            FUN = EmpiricalP,
            nullval = abs(x = fake.pcVals[, x])
          )))
        },
        mc.cores = num.cores
      )
    )
  )
  colnames(x = object@jackStraw.empP) <- paste0("PC", 1:ncol(x = object@jackStraw.empP))
  return(object)
}


# # Documentation
# JackStrawFull <- function(
#   object,
#   num.pc = 5,
#   num.replicate = 100,
#   prop.freq = 0.01
# ) {
#   pc.genes <- rownames(x = object@pca.x)
#   if (length(x = pc.genes) < 200) {
#     prop.freq <- max(prop.freq, 0.015)
#   }
#   md.x <- as.matrix(x = object@pca.x)
#   md.rot <- as.matrix(x = object@pca.rot)
#   real.fval <- sapply(
#     X = 1:num.pc,
#     FUN = function(x) {
#       return(unlist(x = lapply(
#         X = pc.genes,
#         FUN = JackF,
#         r1 = x,
#         r2 = x,
#         x = md.x,
#         rot = md.rot
#       )))
#     }
#   )
#   rownames(x = real.fval) <- pc.genes
#   object@real.fval <- data.frame(real.fval)
#   fake.fval <- sapply(
#     X = 1:num.pc,
#     FUN = function(x) {
#       return(unlist(x = replicate(
#         n = num.replicate,
#         expr = JackstrawF(
#           prop = prop.freq,
#           data = object@scale.data[pc.genes, ],
#           myR1 = x,
#           myR2 = x
#         ),
#         simplify = FALSE
#       )))
#     }
#   )
#   rownames(x = fake.fval) <- 1:nrow(x = fake.fval)
#   object@fake.fval <- data.frame(fake.fval)
#   object@emp.pval <- data.frame(
#     sapply(
#       X = 1:num.pc,
#       FUN = function(x) {
#         return(unlist(x = lapply(
#           X = object@real.fval[, x],
#           FUN = EmpiricalP,
#           nullval = object@fake.fval[, x]
#         )))
#       }
#     )
#   )
#   rownames(x = object@emp.pval) <- pc.genes
#   colnames(x = object@emp.pval) <- paste0("PC", 1:ncol(x = object@emp.pval),)
#   return(object)
# }

# Documentation
# #' @export
#
# JackStrawPermutationTest <- function(
#   object,
#   genes.use = NULL,
#   num.iter = 100,
#   thresh.use = 0.05,
#   do.print = TRUE,
#   k.seed = 1
# ) {
#   genes.use <- SetIfNull(x = genes.use, default = rownames(x = object@pca.x))
#   genes.use <- intersect(x = genes.use, y = rownames(x = object@scale.data))
#   data.use <- t(x = as.matrix(x = object@scale.data[genes.use, ]))
#   if (do.print) {
#     print(paste("Running", num.iter, "iterations"))
#   }
#   pa.object <- permutationPA(
#     data.use,
#     B = num.iter,
#     threshold = thresh.use,
#     verbose = do.print,
#     seed = k.seed
#   )
#   if (do.print) {
#     cat("\n\n")
#   }
#   if (do.print) {
#     print(paste("JackStraw returns", pa.object$r, "significant components"))
#   }
#   return(pa.object)
# }

#' Run Canonical Correlation Analysis (CCA) on multimodal data
#'
#' CCA finds a shared correlation structure betwen two different datasets, enabling integrated downstream analysis
#'
#' @param object Seurat object
#' @param assay.1 First assay for multimodal analysis. Default is RNA
#' @param assay.2 Second assay for multimodal analysis. Default is CITE for CITE-Seq analysis.
#' @param features.1 Features of assay 1 to consider (default is variable genes)
#' @param features.2 Features of assay 2 to consider (default is all features, i.e. for CITE-Seq, all antibodies)
#' @param num.cc Number of canonical correlations to compute and store. Default is 20, but will calculate less if either assay has <20 features.
#' @param normalize.variance Z-score the embedding of each CC to 1, so each CC contributes equally in downstream analysis (default is T)
#'
#' @return Returns object after CCA, with results stored in dimensional reduction cca.assay1 (ie. cca.RNA) and cca.assay2. For example, results can be visualized using DimPlot(object,reduction.use="cca.RNA")
#'
#' @export
#'
MultiModal_CCA <- function(
  object,
  assay.1 = "RNA",
  assay.2 = "CITE",
  features.1 = NULL,
  features.2 = NULL,
  num.cc = 20,
  normalize.variance = TRUE
) {
  #first pull out data, define features
  data.1 <- GetAssayData(
    object = object,
    assay.type = assay.1,
    slot = "scale.data"
  )
  data.2 <- GetAssayData(
    object = object,
    assay.type = assay.2,
    slot = "scale.data"
  )
  if (is.null(x = features.1)) {
    if ((assay.1 == "RNA") && length(x = object@var.genes) > 0) {
      features.1 <- object@var.genes
    } else {
      features.1 <- rownames(x = data.1)
    }
  }
  features.2 <- SetIfNull(x = features.2, default = rownames(x = data.2))
  data.1 <- t(x = data.1[features.1, ])
  data.2 <- t(x = data.2[features.2, ])
  num.cc <- min(20, min(length(x = features.1), length(x = features.2)))
  cca.data <- list(data.1, data.2)
  names(x = cca.data) <- c(assay.1, assay.2)
  # now run CCA
  out <- CanonCor(mat1 = cca.data[[1]],
                  mat2 = cca.data[[2]],
                  standardize = TRUE,
                  k = num.cc)
  cca.output <- list(out$u, out$v)
  embeddings.cca <- list()
  for (i in 1:length(x = cca.data)) {
    assay.use <- names(x = cca.data)[i]
    rownames(x = cca.output[[i]]) <- colnames(x = cca.data[[i]])
    embeddings.cca[[i]] <- cca.data[[i]] %*% cca.output[[i]]
    colnames(x = embeddings.cca[[i]]) <- paste0(
      assay.use,
      "CC",
      1:ncol(x = embeddings.cca[[i]])
    )
    colnames(x = cca.output[[i]]) <- colnames(x = embeddings.cca[[i]])
    if (normalize.variance) {
      embeddings.cca[[i]] <- scale(x = embeddings.cca[[i]])
    }
    object <- SetDimReduction(
      object = object,
      reduction.type = paste0(assay.use, "CCA"),
      slot = "cell.embeddings",
      new.data = embeddings.cca[[i]]
    )
    object <- SetDimReduction(
      object = object,
      reduction.type = paste0(assay.use, "CCA"),
      slot = "key",
      new.data = paste0(assay.use, "CC")
    )
    object <- SetDimReduction(
      object = object,
      reduction.type = paste0(assay.use, "CCA"),
      slot = "x",
      new.data =  cca.output[[i]]
    )
  }
  return(object)
}

#' Run coinertia analysis on multimodal data
#'
#' CIA finds a shared correlation structure betwen two different datasets, enabling integrated downstream analysis
#'
#' @param object Seurat object
#' @param assay.1 First assay for multimodal analysis. Default is RNA
#' @param assay.2 Second assay for multimodal analysis. Default is CITE for CITE-Seq analysis.
#' @param features.1 Features of assay 1 to consider (default is variable genes)
#' @param features.2 Features of assay 2 to consider (default is all features, i.e. for CITE-Seq, all antibodies)
#' @param num.axes Number of principal axes to compute and store. Default is 20, but will calculate less if either assay has <20 features.
#' @param normalize.variance Return the normalized row scares, so each aexis contributes equally in downstream analysis (default is T)
#'
#' @importFrom utils installed.packages
#'
#' @return Returns object after CIA, with results stored in dimensional reduction cia.assay1 (ie. cia.RNA) and cia.assay2. For example, results can be visualized using DimPlot(object,reduction.use="cia.RNA")
#'
#' @export
#'
MultiModal_CIA <- function(
  object,
  assay.1 = "RNA",
  assay.2 = "CITE",
  features.1 = NULL,
  features.2 = NULL,
  num.axes = 20,
  normalize.variance = TRUE
) {
  if (!'made4' %in% rownames(x = installed.packages())) {
      stop("Please install made4")
  }
  #first pull out data, define features
  data.1 <- GetAssayData(
    object = object,
    assay.type = assay.1,
    slot = "scale.data"
  )
  data.2 <- GetAssayData(
    object = object,
    assay.type = assay.2,
    slot = "scale.data"
  )
  if (is.null(x = features.1)) {
    if ((assay.1 == "RNA") && length(x = object@var.genes) > 0) {
      features.1 <- object@var.genes
    } else {
      features.1 <- rownames(x = data.1)
    }
  }
  features.2 <- SetIfNull(x = features.2, default = rownames(x = data.2))
  data.1 <- t(x = data.1[features.1, ])
  data.2 <- t(x = data.2[features.2, ])
  num.axes <- min(20, min(length(x = features.1), length(x = features.2)))
  cia.data <- list(data.1, data.2)
  names(x = cia.data) <- c(assay.1, assay.2)
  # now run cia
  out <- made4::cia(
    df1 = t(x = cia.data[[1]]),
    df2 = t(x = cia.data[[2]]),
    cia.nf = num.axes
  )
  out <- out$coinertia
  cia.output <- list(as.matrix(x = out$c1), as.matrix(x = out$l1))
  embeddings.cia.norm <- list(as.matrix(x = out$mX), as.matrix(x = out$mY))
  embeddings.cia <- list(as.matrix(x = out$lX), as.matrix(x = out$lY))
  for (i in 1:length(x = cia.data)) {
    assay.use <- names(x = cia.data)[i]
    #rownames(cia.output[[i]])=colnames(cia.data[[i]])
    if (normalize.variance) {
      embeddings.cia[[i]] <- (embeddings.cia.norm[[i]])
    }
    colnames(x = embeddings.cia[[i]]) <- paste0(
      assay.use,
      "CI",
      1:ncol(x = embeddings.cia[[i]])
    )
    colnames(x = cia.output[[i]]) <- colnames(x = embeddings.cia[[i]])
    object <- SetDimReduction(
      object = object,
      reduction.type = paste("cia", assay.use, sep="_"),
      slot = "cell.embeddings",
      new.data = embeddings.cia[[i]]
    )
    object <- SetDimReduction(
      object = object,
      reduction.type = paste("cia", assay.use, sep="_"),
      slot = "key",
      new.data = paste0(assay.use," CI")
    )
    object <- SetDimReduction(
      object = object,
      reduction.type = paste("cia", assay.use, sep="_"),
      slot = "x",
      new.data =  cia.output[[i]]
    )
  }
  return(object)
}
