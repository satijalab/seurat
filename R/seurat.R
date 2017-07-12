################################################################################
### Seurat

#' The Seurat Class
#'
#' The Seurat object is the center of each single cell analysis. It stores all information
#' associated with the dataset, including data, annotations, analyes, etc. All that is needed
#' to construct a Seurat object is an expression matrix (rows are genes, columns are cells), which
#' should be log-scale
#'
#' Each Seurat object has a number of slots which store information. Key slots to access
#' are listed below.
#'
#'
#'@section Slots:
#'  \describe{
#'    \item{\code{raw.data}:}{\code{"ANY"}, The raw project data }
#'    \item{\code{data}:}{\code{"ANY"}, The expression matrix (log-scale) }
#'    \item{\code{scale.data}:}{\code{"ANY"}, The scaled (after z-scoring
#'    each gene) expression matrix. Used for PCA, ICA, and heatmap plotting}
#'    \item{\code{var.genes}:}{\code{"vector"},  Variable genes across single cells }
#'    \item{\code{is.expr}:}{\code{"numeric"}, Expression threshold to determine if a gene is expressed }
#'    \item{\code{ident}:}{\code{"factor"},  The 'identity class' for each single cell }
#'    \item{\code{data.info}:}{\code{"data.frame"}, Contains information about each cell, starting with # of genes detected (nGene)
#'    the original identity class (orig.ident), user-provided information (through AddMetaData), etc.  }
#'    \item{\code{project.name}:}{\code{"character"}, Name of the project (for record keeping) }
#'    \item{\code{dr}:}{\code{"list"}, List of stored dimensional reductions. Named by technique }
#'    \item{\code{assay}:}{\code{"list"}, List of additional assays for multimodal analysis. Named by technique }
#'    \item{\code{tsne.rot}:}{\code{"data.frame"}, Cell coordinates on the t-SNE map }
#'    \item{\code{mean.var}:}{\code{"data.frame"}, The output of the mean/variability analysis for all genes }
#'    \item{\code{imputed}:}{\code{"data.frame"}, Matrix of imputed gene scores }
#'    \item{\code{final.prob}:}{\code{"data.frame"}, For spatial inference, posterior probability of each cell mapping to each bin }
#'    \item{\code{insitu.matrix}:}{\code{"data.frame"}, For spatial inference, the discretized spatial reference map }
#'    \item{\code{cell.names}:}{\code{"vector"},  Names of all single cells (column names of the expression matrix) }
#'    \item{\code{cluster.tree}:}{\code{"list"},  List where the first element is a phylo object containing the
#'    phylogenetic tree relating different identity classes }
#'    \item{\code{snn}:}{\code{"dgCMatrix"}, Sparse matrix object representation of the SNN graph }
#'    \item{\code{snn.k}:}{\code{"numeric"}, k used in the construction of the SNN graph }
#'    \item{\code{calc.params}:}{\code{"list"}, Named list to store all calculation related parameters choices}
#'}
#' @name seurat
#' @rdname seurat
#' @aliases seurat-class
#' @exportClass seurat
#' @importFrom Rcpp evalCpp
#' @useDynLib Seurat

seurat <- setClass(
  "seurat",
  slots = c(
    raw.data = "ANY",
    data = "ANY",
    scale.data = "ANY",
    var.genes = "vector",
    is.expr = "numeric",
    ident = "factor",
    dr = "list",
    assay = "list",
    emp.pval = "data.frame",
    gene.scores = "data.frame",
    drop.coefs = "data.frame",
    wt.matrix = "data.frame",
    drop.wt.matrix = "data.frame",
    trusted.genes = "vector",
    drop.expr = "numeric",
    data.info = "data.frame",
    project.name = "character",
    jackStraw.empP = "data.frame",
    jackStraw.fakePC = "data.frame",
    jackStraw.empP.full = "data.frame",
    kmeans.col = "list",
    mean.var = "data.frame",
    imputed = "data.frame",
    mix.probs = "data.frame",
    mix.param = "data.frame",
    final.prob = "data.frame",
    insitu.matrix = "data.frame",
    tsne.rot = "data.frame",
    cell.names = "vector",
    cluster.tree = "list",
    snn = "dgCMatrix",
    snn.k = "numeric",
    calc.params = "list",
    kmeans="ANY"
  )
)

#   Documentation...
####################
calc.drop.prob <- function(x, a, b) {
  return(exp(x = a + b * x) / (1 + exp(x = a + b * x)))
}

#   Documentation...
####################
#from Jean Fan - thanks!!
custom.dist <- function(my.mat, my.function, ...) {
  n <- ncol(x = my.mat)
  mat <- matrix(data = 0, ncol = n, nrow = n)
  colnames(x = mat) <- rownames(x = mat) <- colnames(x = my.mat)
  for (i in 1:nrow(x = mat)) {
    for(j in 1:ncol(x = mat)) {
      mat[i,j] <- my.function(my.mat[, i], my.mat[, j], ...)
    }
  }
  return(as.dist(mat))
}

#   Documentation
#################
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

#' Probability of detection by identity class
#'
#' For each gene, calculates the probability of detection for each identity
#' class.
#'
#' @param object Seurat object
#' @param thresh.min Minimum threshold to define 'detected' (log-scale)
#'
#' @return Returns a matrix with genes as rows, identity classes as columns.
#'
#' @export
#'
ClusterAlpha <- function(object, thresh.min = 0) {
  ident.use <- object@ident
  data.all <- data.frame(row.names = rownames(x = object@data))
  for (i in sort(x = unique(x = ident.use))) {
    temp.cells <- WhichCells(object = object, ident = i)
    data.temp <- apply(
      X = object@data[, temp.cells],
      MARGIN = 1,
      FUN = function(x) {
        return(sum(x > thresh.min)/length(x = x))
      }
    )
    data.all <- cbind(data.all, data.temp)
    colnames(x = data.all)[ncol(x = data.all)] <- i
  }
  colnames(x = data.all) <- sort(x = unique(x = ident.use))
  return(data.all)
}

#' Average PCA scores by identity class
#'
#' Returns the PCA scores for an 'average' single cell in each identity class
#'
#' @param object Seurat object
#'
#' @return Returns a matrix with genes as rows, identity classes as columns
#'
#' @export
#'
AveragePCA <- function(object) {
  ident.use <- object@ident
  data.all <- data.frame(row.names = colnames(x = object@pca.rot))
  for (i in levels(x = ident.use)) {
    temp.cells <- WhichCells(object = object, ident = i)
    if (length(x = temp.cells) == 1) {
      data.temp <- apply(
        X = data.frame((object@pca.rot[c(temp.cells, temp.cells), ])),
        MARGIN = 2,
        FUN = mean
      )
    }
    if (length(x = temp.cells) > 1) {
      data.temp <- apply(
        X = object@pca.rot[temp.cells, ],
        MARGIN = 2,
        FUN = mean
      )
    }
    data.all <- cbind(data.all, data.temp)
    colnames(x = data.all)[ncol(x = data.all)] <- i
  }
  #colnames(data.all)=levels(ident.use)
  return((data.all))
}

#' Averaged gene expression by identity class
#'
#' Returns gene expression for an 'average' single cell in each identity class
#'
#' Output is in log-space, but averaging is done in non-log space.
#'
#' @param object Seurat object
#' @param genes.use Genes to analyze. Default is all genes.
#' @param return.seurat Whether to return the data as a Seurat object. Default is false.
#' @param add.ident Place an additional label on each cell prior to averaging (very useful if you want to observe cluster averages, separated by replicate, for example)
#' @inheritParams FetchData
#' @param show.progress Show progress bar (default is T)
#' @param ... Arguments to be passed to methods such as \code{\link{Setup}}
#' @return Returns a matrix with genes as rows, identity classes as columns.
#' @export
AverageExpression <- function(
  object,
  genes.use = NULL,
  return.seurat = FALSE,
  add.ident = NULL,
  use.scale = FALSE,
  use.raw = FALSE,
  show.progress = TRUE,
  ...
) {

  ident.orig <- object@ident
  orig.levels <- levels(x = object@ident)
  ident.new <- c()
  if (! is.null(x = add.ident)) {
    new.data <- FetchData(object = object, vars.all = add.ident)
    new.ident <- paste(
      object@ident[rownames(x = new.data)],
      new.data[, 1],
      sep = '_'
    )
    object <- SetIdent(
      object = object,
      cells.use = rownames(x = new.data),
      ident.use = new.ident
    )
  }
  if (return.seurat) {
    assays.use <- c("RNA", names(x = object@assay))
  } else {
    assays.use <- "RNA"
  }
  slot.use <- "data"
  fxn.average <- expMean
  if (show.progress) {
    fxn.loop <- pbsapply
  } else {
    fxn.loop <- sapply
  }
  if (use.scale) {
    slot.use <- "scale.data"
    fxn.average <- mean
  }
  if (use.raw) {
    slot.use <- "raw.data"
    fxn.average <- mean
  }
  data.return <- list()
  for (i in 1:length(x = assays.use)) {
    data.use <- GetAssayData(
      object = object,
      assay.type = assays.use[i],
      slot = slot.use
    )
    genes.assay <- genes.use
    if (length(x = intersect(x = genes.use,y = rownames(x = data.use))) <1 ) {
      genes.assay <- rownames(x = data.use)
    }
    data.all <- data.frame(row.names = genes.assay)
    for (j in levels(x = object@ident)) {
      temp.cells <- WhichCells(object = object, ident = j)
      genes.assay <- unique(x = intersect(x = genes.assay, y = rownames(x = data.use)))
      if (length(x = temp.cells) == 1) {
        data.temp <- (data.use[genes.assay, temp.cells])
      }
      if (length(x = temp.cells) >1 ) {
        data.temp <- apply(
          X = data.use[genes.assay, temp.cells],
          MARGIN = 1,
          FUN = fxn.average
        )
      }
      data.all <- cbind(data.all, data.temp)
      colnames(x = data.all)[ncol(x = data.all)] <- j
      if (show.progress) {
        print(paste0("Finished averaging ", assays.use[i], " for cluster ", j))
      }
      if(i == 1) {
        ident.new <- c(ident.new, as.character(x = ident.orig[temp.cells[1]]))
      }
    }
    names(x = ident.new) <- levels(x = object@ident)
    data.return[[i]] <- data.all
    names(x = data.return)[i] <- assays.use[[i]]
  }
  if (return.seurat) {
    toRet <- new(Class = "seurat", raw.data = data.return[[1]])
    toRet <- Setup(
      object = toRet,
      project = "Average",
      min.cells = 0,
      min.genes = 0,
      is.expr = 0,
      ...
    )
    #for multimodal data
    if (length(x = data.return) > 1) {
      for (i in 2:length(x = data.return)) {
        toRet <- SetAssayData(
          object = toRet,
          assay.type = names(x = data.return)[i],
          slot = "raw.data",
          new.data = data.return[[i]]
        )
      }
    }
    toRet <- SetIdent(
      object = toRet,
      cells.use = toRet@cell.names,
      ident.use = ident.new[toRet@cell.names]
    )
    toRet@ident <- factor(
      x = toRet@ident,
      levels = as.character(x = orig.levels),
      ordered = TRUE
    )
    return(toRet)
  } else {
    return(data.return[[1]])
  }
}

#   Documentation...
####################
RegulatorScore <- function(
  object,
  candidate.reg,
  score.name,
  cells.use = NULL
) {
  cells.use <- SetIfNull(x = cells.use, default = colnames(x = object@data))
  candidate.reg <- candidate.reg[candidate.reg %in% rownames(x = object@data)]
  my.score <- retreiveScore(object, score.name)[cells.use]
  my.data <- object@data[, cells.use]
  my.ident <- object@ident[cells.use]
  reg.score <- unlist(
    lapply(
      X = candidate.reg,
      FUN = regressionSig,
      score = my.score,
      data = my.data,
      latent = my.ident,
      code = "rsem"
    )
  )
  names(x = reg.score) <- candidate.reg
  return(reg.score)
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
  genes.use <- SetIfNull(x = genes.use, y = rownames(x = object@data))
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

#' Merge subchilden of a node
#'
#' Merge the subchilden of a node into a single identity class
#'
#' @param object Seurat object
#' @param node.use Merge subchildren of this node
#'
#' @export
#'
MergeNode <- function(object, node.use = NULL) {
  object.tree <- object@cluster.tree[[1]]
  node.children <- DFT(
    tree = object.tree,
    node = node.use,
    include.children = TRUE
  )
  node.children <- ainb(a = node.children, b = levels(x = object@ident))
  children.cells <- WhichCells(object = object, ident = node.children)
  if (length(x = children.cells > 0)) {
    object <- SetIdent(
      object = object,
      cells.use = children.cells,
      ident.use = min(node.children)
    )
  }
  return(object)
}

# Documentation
###############
#Not documented for now
#' @export
PosteriorPlot <- function(object, name) {
  post.names <- colnames(x = subc(data = object@mix.probs, code = name))
  VlnPlot(
    object = object,
    features.plot = post.names,
    inc.first=TRUE,
    inc.final=TRUE,
    by.k=TRUE
  )
}

# Documentation
###############
#Internal, not documented for now
map.cell.score <- function(gene, gene.value, insitu.bin, mu, sigma, alpha) {
  code.1 <- paste(gene, insitu.bin, sep=".")
  mu.use <- mu[paste(code.1, "mu", sep="."), 1]
  sigma.use <- sigma[paste(code.1, "sigma", sep="."), 1]
  alpha.use <- alpha[paste(code.1, "alpha", sep="."), 1]
  bin.prob <- unlist(
    x = lapply(
      X = 1:length(x = insitu.bin),
      FUN = function(x) {
        return(dnorm(
          x = gene.value,
          mean = mu.use[x],
          sd = sigma.use[x],
          log = TRUE) + log(x = alpha.use[x])
        )
      }
    )
  )
  return(bin.prob)
}

#Internal, not documented for now
#' @export
#'
MapCell <- function(
  object,
  cell.name,
  do.plot = FALSE,
  safe.use = TRUE,
  text.val = NULL,
  do.rev = FALSE
) {
  insitu.matrix <- object@insitu.matrix
  insitu.genes <- colnames(x = insitu.matrix)
  insitu.genes <- insitu.genes[insitu.genes %in% rownames(x = object@imputed)]
  insitu.use <- insitu.matrix[, insitu.genes]
  imputed.use <- object@imputed[insitu.genes, ]
  if (safe.use) {
    safe_fxn <- log_add
  } else {
    safe_fxn <- sum
  }
  all.needed.cols <- unique(
    x = unlist(
      x = lapply(
        X = insitu.genes,
        FUN = function(x) {
          return(paste(x, insitu.use[, x], "post", sep="."))
        }
      )
    )
  )
  missing.cols <- which(! (all.needed.cols %in% colnames(x = object@mix.probs)))
  if (length(x = missing.cols) > 0) {
    stop(paste(
      "Error: ",
      all.needed.cols[missing.cols],
      "is missing from the mixture fits"
    ))
  }
  all.probs <- data.frame(
    sapply(
      X = insitu.genes,
      FUN = function(x) {
        return(
          log(x = as.numeric(x = object@mix.probs[
            cell.name, # Row
            paste(x, insitu.use[, x], "post", sep=".") # Column
          ])))
      }
    )
  )
  scale.probs <- t(
    x = t(x = all.probs) - apply(X = t(x = all.probs), MARGIN = 1, FUN = log_add)
  )
  scale.probs[scale.probs < (-9.2)] <- (-9.2)
  #head(scale.probs)
  total.prob <- exp(x = apply(X = scale.probs, MARGIN = 1, FUN = safe_fxn))
  total.prob <- total.prob / sum(total.prob)
  if (do.plot) {
    #plot(total.prob,main=cell.name)
    par(mfrow = c(1, 2))
    txt.matrix <- matrix(data = rep(x = "", 64), nrow=8, ncol=8)
    if (! is.null(x = text.val)) {
      txt.matrix[text.val] <- "X"
    }
    if (do.rev) {
      scale.probs <- scale.probs[unlist(
        x = lapply(
          X = 0:7,
          FUN = function(x) {
            return(seq(from = 1, to = 57, by = 8) + x)
          }
        )
      ), ]
    }
    aheatmap(
      x = matrix(data = total.prob, nrow=8, ncol=8),
      Rowv = NA,
      Colv = NA,
      txt = txt.matrix,
      col = bwCols
    )
    aheatmap(x = scale.probs, Rowv = NA, Colv = NA)
    ResetPar()
  }
  return(total.prob)
}

# Documentation
###############
#Internal, not documented for now
CalcInsitu <- function(
  object,
  gene,
  do.plot = TRUE,
  do.write = FALSE,
  write.dir = "~/window/insitu/",
  col.use = bwCols,
  do.norm = FALSE,
  cells.use = NULL,
  do.return = FALSE,
  probs.min = 0,
  do.log = FALSE,
  use.imputed = FALSE,
  bleach.use = 0
) {
  cells.use <- SetIfNull(x = cells.use, default = colnames(x = object@final.prob))
  probs.use <- object@final.prob
  if (use.imputed) {
    data.use <- exp(x = object@imputed) - 1
  } else {
    data.use <- exp(object@data) - 1
  }
  cells.use <- cells.use[cells.use %in% colnames(x = probs.use)]
  cells.use <- cells.use[cells.use %in% colnames(x = data.use)]
  #insilico.stain=matrix(unlist(lapply(1:64,function(x) sum(probs.use[x,]*data.use[gene,]))),nrow=8,ncol=8)
  insilico.vector <- unlist(
    x = lapply(
      X = 1:64,
      FUN = function(x) {
        return(sum(
          as.numeric(x = probs.use[x, cells.use]) *
            as.numeric(x = data.use[gene, cells.use])
        ))
      }
    )
  )
  probs.total <- apply(X = probs.use, MARGIN = 1, FUN = sum)
  probs.total[probs.total < probs.min] <- probs.min
  insilico.stain <- (matrix(data = insilico.vector / probs.total, nrow=8, ncol=8))
  if (do.log) {
    insilico.stain <- log(x = insilico.stain + 1)
  }
  if (bleach.use > 0) {
    insilico.stain <- insilico.stain - bleach.use
    insilico.stain <- minmax(data = insilico.stain, min=0, max=1e6)
  }
  if (do.norm) {
    insilico.stain <- (insilico.stain - min(insilico.stain)) /
      (max(insilico.stain) - min(insilico.stain))
  }
  title.use <- gene
  if (gene %in% colnames(x = object@insitu.matrix)) {
    pred.use <- prediction(
      predictions = insilico.vector / probs.total,
      labels = object@insitu.matrix[, gene],
      label.ordering = 0:1
    )
    perf.use <- performance(prediction.obj = pred.use, measure = "auc")
    auc.use <- round(x = perf.use@y.values[[1]], digits = 3)
    title.use <- paste(gene, sep=" ")
  }
  if (do.write) {
    write.table(
      x = insilico.stain,
      file = paste0(write.dir, gene, ".txt"),
      quote=FALSE,
      row.names=FALSE,
      col.names=FALSE
    )
  }
  if (do.plot) {
    aheatmap(
      x = insilico.stain,
      Rowv = NA,
      Colv = NA,
      col = col.use,
      main=title.use
    )
  }
  if (do.return) {
    return(as.vector(x = insilico.stain))
  }
  return(object)
}

#' Build mixture models of gene expression
#'
#' Models the imputed gene expression values as a mixture of gaussian
#' distributions. For a two-state model, estimates the probability that a given
#' cell is in the 'on' or 'off' state for any gene. Followed by a greedy
#' k-means step where cells are allowed to flip states based on the overall
#' structure of the data (see Manuscript for details)
#'
#' @param object Seurat object
#' @param gene Gene to fit
#' @param do.k Number of modes for the mixture model (default is 2)
#' @param num.iter Number of 'greedy k-means' iterations (default is 1)
#' @param do.plot Plot mixture model results
#' @param genes.use Genes to use in the greedy k-means step (See manuscript for details)
#' @param start.pct Initial estimates of the percentage of cells in the 'on'
#' state (usually estimated from the in situ map)
#'
#' @return A Seurat object, where the posterior of each cell being in the 'on'
#' or 'off' state for each gene is stored in object@@mix.probs
#'
#' @importFrom mixtools normalmixEM
#'
#' @export
#'
FitGeneK <- function(
  object,
  gene,
  do.k = 2,
  num.iter = 1,
  do.plot = FALSE,
  genes.use = NULL,
  start.pct = NULL
) {
  data <- object@imputed
  data.use <- data[gene, ]
  names(x = data.use) <- colnames(x = data.use)
  scale.data <- t(x = scale(x = t(x = object@imputed)))
  genes.use <- SetIfNull(x = genes.use, default = rownames(x = scale.data))
  genes.use <- genes.use[genes.use %in% rownames(x = scale.data)]
  scale.data <- scale.data[genes.use, ]
  data.cut <- as.numeric(x = data.use[gene, ])
  cell.ident <- as.numeric(x = cut(x = data.cut, breaks = do.k))
  if (! (is.null(x = start.pct))) {
    cell.ident <- rep.int(x = 1, times = length(x = data.cut))
    cell.ident[data.cut > quantile(x = data.cut, probs = 1 - start.pct)] <- 2
  }
  cell.ident <- order(tapply(
    X = as.numeric(x = data.use),
    INDEX = cell.ident,
    FUN = mean
  ))[cell.ident]
  ident.table <- table(cell.ident)
  if (num.iter > 0) {
    for (i2 in 1:num.iter) {
      cell.ident <- iter.k.fit(
        scale.data = scale.data,
        cell.ident = cell.ident,
        data.use = data.use
      )
      ident.table <- table(cell.ident)
    }
  }
  ident.table <- table(cell.ident)
  raw.probs <- t(
    x = sapply(
      X = data.use,
      FUN = function(y) {
        return(unlist(
          x = lapply(
            X = 1:do.k,
            FUN = function(x) {
              return(
                (ident.table[x] / sum(ident.table)) * dnorm(
                  x = y,
                  mean = mean(x = as.numeric(x = data.use[cell.ident == x])),
                  sd = sd(x = as.numeric(x = data.use[cell.ident == x]))
                )
              )
            }
          )
        ))
      }
    )
  )
  norm.probs <- raw.probs / apply(X = raw.probs, MARGIN = 1, FUN = sum)
  colnames(x = norm.probs) <- unlist(
    x = lapply(
      X = 1:do.k,
      FUN = function(x) {
        paste(gene, x - 1, "post", sep=".")
      }
    )
  )
  norm.probs <- cbind(norm.probs, cell.ident)
  colnames(x = norm.probs)[ncol(x = norm.probs)] <- paste0(gene, ".ident")
  new.mix.probs <- data.frame(
    minusc(data = object@mix.probs, code = paste0(gene, ".")),
    row.names = rownames(x = object@mix.probs)
  )
  colnames(x = new.mix.probs)[1] <- "nGene"
  object@mix.probs <- cbind(new.mix.probs, norm.probs)
  if (do.plot) {
    nCol <- 2
    num.row <- floor(x = (do.k + 1) / nCol - (1e-5)) + 1
    hist(
      x = as.numeric(x = data.use),
      probability = TRUE,
      ylim = c(0, 1),
      xlab = gene,
      main = gene
    )
    for (i in 1:do.k) {
      lines(
        x = seq(from = -10, to = 10, by = 0.01),
        y = (ident.table[i] / sum(ident.table)) * dnorm(
          x = seq(from = -10, to = 10, by = 0.01),
          mean = mean(x = as.numeric(x = data.use[cell.ident == i])),
          sd = sd(x = as.numeric(x = data.use[cell.ident == i]))
        ),
        col=i,
        lwd=2
      )
    }
  }
  return(object)
}

#Internal, not documented for now
iter.k.fit <- function(scale.data, cell.ident, data.use) {
  means.all <- sapply(
    X = sort(x = unique(x = cell.ident)),
    FUN = function(x) {
      return(apply(X = scale.data[, cell.ident == x], MARGIN = 1, FUN = mean))
    }
  )
  all.dist <- data.frame(
    t(x = sapply(
      X = 1:ncol(x = scale.data),
      FUN = function(x) {
        return(unlist(x = lapply(
          X = sort(x = unique(x = cell.ident)),
          FUN = function(y) {
            return(dist(x = rbind(scale.data[, x], means.all[, y])))
          }
        )))
      }
    ))
  )
  cell.ident <- apply(X = all.dist, MARGIN = 1, FUN = which.min)
  cell.ident <- order(tapply(
    X = as.numeric(x = data.use),
    INDEX = cell.ident,
    FUN = mean
  ))[cell.ident]
  return(cell.ident)
}

# Documentation
###############
#Internal, not documented for now
#' @export
FitGeneMix <- function(
  object,
  gene,
  do.k = 3,
  use.mixtools = TRUE,
  do.plot = FALSE,
  plot.with.imputed = TRUE,
  min.bin.size = 10
) {
  data.fit <- as.numeric(x = object@imputed[gene, ])
  mixtools.fit <- normalmixEM(x = data.fit, k = do.k)
  comp.order <- order(mixtools.fit$mu)
  mixtools.posterior <- data.frame(mixtools.fit$posterior[, comp.order])
  colnames(x = mixtools.posterior) <- unlist(
    x = lapply(
      X = 1:do.k,
      FUN = function(x) {
        return(paste(gene, x - 1, "post", sep="."))
      }
    )
  )
  #mixtools.mu=data.frame(mixtools.fit$mu[comp.order])
  #mixtools.sigma=data.frame(mixtools.fit$sigma[comp.order])
  #mixtools.alpha=data.frame(mixtools.fit$lambda[comp.order])
  #rownames(mixtools.mu)=unlist(lapply(1:do.k,function(x)paste(gene,x-1,"mu",sep=".")))
  #rownames(mixtools.sigma)=unlist(lapply(1:do.k,function(x)paste(gene,x-1,"sigma",sep=".")))
  #rownames(mixtools.alpha)=unlist(lapply(1:do.k,function(x)paste(gene,x-1,"alpha",sep=".")))
  #object@mix.mu = rbind(minusr(object@mix.mu,gene), mixtools.mu);
  #object@mix.sigma = rbind(minusr(object@mix.sigma,gene), mixtools.sigma);
  #o#bject@mu.alpha =rbind(minusr(object@mu.alpha,gene), mixtools.alpha);
  if (do.plot) {
    nCol <- 2
    num.row <- floor(x = (do.k + 1) / nCol - (1e-5)) + 1
    par(mfrow = c(num.row, nCol))
    plot.mixEM(x = mixtools.fit, which = 2)
    plot.data <- as.numeric(x = object@imputed[gene, ])
    if (! plot.with.imputed) {
      plot.data <- as.numeric(x = object@data[gene, ])
    }
    unlist(
      x = lapply(
        X = 1:do.k,
        FUN = function(x) {
          plot(
            x = plot.data,
            y = mixtools.posterior[, x],
            ylab = paste0("Posterior for Component ", x - 1),
            xlab = gene,
            main = gene
          )
        }
      )
    )
  }
  new.mix.probs <- data.frame(
    minusc(data = object@mix.probs, code = paste0(gene, ".")),
    row.names = rownames(x = object@mix.probs)
  )
  colnames(x = new.mix.probs)[1] <- "nGene"
  object@mix.probs <- cbind(new.mix.probs, mixtools.posterior)
  return(object)
}

# Documentation
###############
#Internal, not documented for now
lasso.fxn <- function(
  lasso.input,
  genes.obs,
  s.use = 20,
  gene.name = NULL,
  do.print = FALSE,
  gram = TRUE
) {
  lasso.model <- lars(
    x = lasso.input,
    y = as.numeric(x = genes.obs),
    type = "lasso",
    max.steps = s.use * 2,
    use.Gram = gram
  )
  #lasso.fits=predict.lars(lasso.model,lasso.input,type="fit",s=min(s.use,max(lasso.model$df)))$fit
  lasso.fits <- predict.lars(
    object = lasso.model,
    newx = lasso.input,
    type = "fit",
    s = s.use
  )$fit
  if (do.print) {
    print(gene.name)
  }
  return(lasso.fits)
}

#' Calculate smoothed expression values
#'
#'
#' Smooths expression values across the k-nearest neighbors based on dimensional reduction
#'
#' @inheritParams FeaturePlot
#' @inheritParams AddImputedScore
#' @param genes.fit Genes to calculate smoothed values for
#' @param k k-param for k-nearest neighbor calculation
#' @param do.log Whether to perform smoothing in log space. Default is false.
#'
#' @importFrom FNN get.knn
#'
#' @export
#'
AddSmoothedScore <- function(
  object,
  genes.fit = NULL,
  dim.1 = 1,
  dim.2 = 2,
  reduction.use = "tSNE",
  k = 30,
  do.log = FALSE,
  do.print = FALSE
) {
  genes.fit <- SetIfNull(x = genes.fit, default = object@var.genes)
  genes.fit <- genes.fit[genes.fit %in% rownames(x = object@data)]
  dim.code <- GetDimReduction(
    object = object,
    reduction.type = reduction.use,
    slot = 'key'
  )
  dim.codes <- paste0(dim.code, c(dim.1, dim.2))
  data.plot <- FetchData(object = object, vars.all = dim.codes)
  knn.smooth <- get.knn(data = data.plot, k = k)$nn.index
  avg.fxn <- mean
  if (! do.log) {
    avg.fxn <- expMean
  }
  lasso.fits <- data.frame(
    t(
      x = sapply(
        X = genes.fit,
        FUN = function(g) {
          return(unlist(
            x = lapply(
              X = 1:nrow(x = data.plot),
              FUN = function(y) {
                avg.fxn(as.numeric(x = object@data[g, knn.smooth[y, ]]))
              }
            )
          ))
        }
      )
    )
  )
  colnames(x = lasso.fits) <- rownames(x = data.plot)
  genes.old <- genes.fit[genes.fit %in% rownames(x = object@imputed)]
  genes.new <- genes.fit[! (genes.fit %in% rownames(x = object@imputed))]
  if (length(x = genes.old) > 0) {
    object@imputed[genes.old, ] <- lasso.fits[genes.old, ]
  }
  object@imputed <- rbind(object@imputed, lasso.fits[genes.new, ])
  return(object)
}

#' Calculate imputed expression values
#'
#' Uses L1-constrained linear models (LASSO) to impute single cell gene
#' expression values.
#'
#'
#' @param object Seurat object
#' @param genes.use A vector of genes (predictors) that can be used for
#' building the LASSO models.
#' @param genes.fit A vector of genes to impute values for
#' @param s.use Maximum number of steps taken by the algorithm (lower values
#' indicate a greater degree of smoothing)
#' @param do.print Print progress (output the name of each gene after it has
#' been imputed).
#' @param gram The use.gram argument passed to lars
#'
#' @return Returns a Seurat object where the imputed values have been added to
#' object@@data
#'
#' @import lars
#'
#' @export
#'
AddImputedScore <- function(
  object,
  genes.use = NULL,
  genes.fit = NULL,
  s.use = 20,
  do.print = FALSE,
  gram = TRUE
) {
  genes.use <- SetIfNull(x = genes.use, default = object@var.genes)
  genes.fit <- SetIfNull(x = genes.fit, default = object@var.genes)
  genes.use <- genes.use[genes.use %in% rownames(x = object@data)]
  genes.fit <- genes.fit[genes.fit %in% rownames(x = object@data)]
  lasso.input <- t(x = object@data[genes.use, ])
  lasso.fits <- data.frame(t(
    x = sapply(
      X = genes.fit,
      FUN = function(x) {
        return(
          lasso.fxn(
            lasso.input = t(x = object@data[genes.use[genes.use != x], ]),
            genes.obs = object@data[x, ],
            s.use = s.use,
            gene.name = x,
            do.print = do.print,
            gram = gram
          )
        )
      }
    )
  ))
  genes.old <- genes.fit[genes.fit %in% rownames(x = object@imputed)]
  genes.new <- genes.fit[! (genes.fit %in% rownames(x = object@imputed))]
  if (length(x = genes.old) > 0) {
    object@imputed[genes.old, ] <- lasso.fits[genes.old, ]
  }
  object@imputed <- rbind(object@imputed, lasso.fits[genes.new, ])
  return(object)
}

# Documentation
###############
#Cool, but not supported right now
SpatialDe <- function(object, marker.cells, genes.use = NULL) {
  object <- p15
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
        object = p15,
        cells.1 = marker.cells,
        cells.2 = embed.diff,
        genes.use = genes.use
      ),
      subset = p_val < (1e-5)
    )
  )
  diff.genes <- subset(
    x = DiffExpTest(
      object = p15,
      cells.1 = marker.cells,
      cells.2 = embed.diff,
      genes.use = diff.genes
    ),
    subset = p_val<(1e-10)
  )
  return(diff.genes)
}

#' Perform spectral density clustering on single cells
#'
#' Find point clounds single cells in a two-dimensional space using density clustering (DBSCAN).
#'
#' @param object Seurat object
#' @param dim.1 First dimension to use
#' @param dim.2 second dimension to use
#' @param reduction.use Which dimensional reduction to use (either 'pca' or 'ica')
#' @param G.use Parameter for the density clustering. Lower value to get more fine-scale clustering
#' @param set.ident TRUE by default. Set identity class to the results of the density clustering.
#' Unassigned cells (cells that cannot be assigned a cluster) are placed in cluster 1, if there are any.
#' @param seed.use Random seed for the dbscan function
#' @param ... Additional arguments to be passed to the dbscan function
#'
#' @export
#'
DBClustDimension <- function(
  object,
  dim.1 = 1,
  dim.2 = 2,
  reduction.use = "tsne",
  G.use = NULL,
  set.ident = TRUE,
  seed.use = 1,
  ...
) {
  dim.code <- GetDimReduction(
    object = object,
    reduction.type = reduction.use,
    slot = 'key'
  )
  dim.codes <- paste0(dim.code, c(dim.1, dim.2))
  data.plot <- FetchData(object = object, vars.all = dim.codes)
  x1 <- paste0(dim.code, dim.1)
  x2 <- paste0(dim.code, dim.2)
  data.plot$x <- data.plot[, x1]
  data.plot$y <- data.plot[, x2]
  set.seed(seed = seed.use)
  data.mclust <- ds <- dbscan(data = data.plot[, c("x", "y")], eps = G.use, ...)
  to.set <- as.numeric(x = data.mclust$cluster + 1)
  data.names <- names(x = object@ident)
  object@data.info[data.names, "DBclust.ident"] <- to.set
  if (set.ident) {
    object@ident <- factor(x = to.set)
    names(x = object@ident) <- data.names
  }
  return(object)
}

#' Perform spectral k-means clustering on single cells
#'
#' Find point clounds single cells in a low-dimensional space using k-means clustering.
#'
#' CAn be useful for smaller datasets, not documented here yet
#' @export
#'
KClustDimension <- function(
  object,
  dim.1 = 1,
  dim.2 = 2,
  cells.use = NULL,
  pt.size = 4,
  reduction.use = "tsne",
  k.use = 5,
  set.ident = FALSE,
  seed.use = 1,
  ...
) {
  cells.use <- SetIfNull(x = cells.use, default = colnames(x = object@data))
  dim.code <- "PC"
  if (reduction.use == "pca") {
    data.plot <- object@pca.rot[cells.use, ]
  }
  if (reduction.use == "tsne") {
    data.plot <- object@tsne.rot[cells.use, ]
    dim.code <- "Seurat_"
  }
  if (reduction.use == "ica") {
    data.plot <- object@ica.rot[cells.use, ]
    dim.code <- "IC"
  }
  x1 <- paste0(dim.code, dim.1)
  x2 <- paste0(dim.code, dim.2)
  data.plot$x <- data.plot[, x1]
  data.plot$y <- data.plot[, x2]
  if (reduction.use != "pca") {
    set.seed(seed = seed.use)
    data.mclust <- ds <- kmeans(x = data.plot[, c("x", "y")], centers = k.use)
  }
  if (reduction.use == "pca") {
    set.seed(seed = seed.use)
    data.mclust <- ds <- kmeans(x = object@pca.rot[cells.use, dim.1], centers = k.use)
  }
  to.set <- as.numeric(x = data.mclust$cluster)
  data.names <- names(x = object@ident)
  object@data.info[data.names, "kdimension.ident"] <- to.set
  if (set.ident) {
    object@ident <- factor(x = to.set)
    names(x = object@ident) <- data.names
  }
  return(object)
}

#' Significant genes from a PCA
#'
#' Returns a set of genes, based on the JackStraw analysis, that have
#' statistically significant associations with a set of PCs.
#'
#' @param object Seurat object
#' @param pcs.use PCS to use.
#' @param pval.cut P-value cutoff
#' @param use.full Use the full list of genes (from the projected PCA). Assumes
#' that ProjectPCA has been run. Default is TRUE
#' @param max.per.pc Maximum number of genes to return per PC. Used to avoid genes from one PC dominating the entire analysis.
#'
#' @return A vector of genes whose p-values are statistically significant for
#' at least one of the given PCs.
#'
#' @export
#'
PCASigGenes <- function(
  object,
  pcs.use,
  pval.cut = 0.1,
  use.full = TRUE,
  max.per.pc = NULL
) {
  pvals.use <- object@jackStraw.empP
  pcx.use <- object@pca.x
  if (use.full) {
    pvals.use <- object@jackStraw.empP.full
    pcx.use <- object@pca.x.full
  }
  if (length(x = pcs.use) == 1) {
    pvals.min <- pvals.use[, pcs.use]
  }
  if (length(x = pcs.use) > 1) {
    pvals.min <- apply(X = pvals.use[, pcs.use], MARGIN = 1, FUN = min)
  }
  names(x = pvals.min) <- rownames(x = pvals.use)
  genes.use <- names(x = pvals.min)[pvals.min < pval.cut]
  if (! is.null(x = max.per.pc)) {
    pc.top.genes <- PCTopGenes(
      object = object,
      pc.use = pcs.use,
      num.genes = max.per.pc,
      use.full = use.full,
      do.balanced = FALSE
    )
    genes.use <- ainb(a = pc.top.genes, b = genes.use)
  }
  return(genes.use)
}


#' GenesInCluster
#'
#' After k-means analysis, previously run with DoKMeans, returns a set of genes associated with each cluster
#'
#' @param object Seurat object. Assumes DoKMeans has already been run
#' @param cluster.num K-means cluster(s) to return genes for
#' @param max.genes max number of genes to return 
#' @return A vector of genes who are members in the cluster.num k-means cluster(s)
#'
#' @export
GenesInCluster <- function(object, cluster.num, max.genes = 1e6) {
  toReturn <- unlist(
    x = lapply(
      X = cluster.num,
      FUN = function(x) {
        return(head(
          x = sort(x = names(x = which(x = object@kmeans@gene.kmeans.obj$cluster==x))),
          n = max.genes
        ))
      }
    )
  )
  return(toReturn)
}

# Documentation
###############
#' @export
#'
CellCorMatrix <- function(
  object,
  cor.genes = NULL,
  cell.inds = NULL,
  do.k = FALSE,
  k.seed = 1,
  k.num = 4,
  vis.low = (-1),
  vis.high = 1,
  vis.one = 0.8,
  pcs.use = 1:3,
  col.use = pyCols
) {
  cor.genes <- SetIfNull(x = cor.genes, default = object@var.genes)
  cell.inds <- SetIfNull(x = cell.inds, default = colnames(x = object@data))
  cor.genes <- cor.genes[cor.genes %in% rownames(x = object@data)]
  data.cor <- object@scale.data[cor.genes, cell.inds]
  cor.matrix <- cor(x = data.cor)
  set.seed(seed = k.seed)
  kmeans.cor <- kmeans(x = cor.matrix, centers = k.num)
  if (do.k) {
    cor.matrix <- cor.matrix[order(kmeans.cor$cluster), order(kmeans.cor$cluster)]
  }
  kmeans.names <- rownames(x = cor.matrix)
  row.annot <- data.frame(
    cbind(
      kmeans.cor$cluster[kmeans.names],
      object@pca.rot[kmeans.names, pcs.use]
    )
  )
  colnames(x = row.annot) <- c("K", paste0("PC", pcs.use))
  cor.matrix[cor.matrix == 1] <- vis.one
  cor.matrix <- minmax(data = cor.matrix, min = vis.low, max = vis.high)
  object@kmeans.cell <- list(kmeans.cor)
  if (do.k) {
    aheatmap(
      x = cor.matrix,
      col = col.use,
      Rowv = NA,
      Colv = NA,
      annRow = row.annot
    )
  } else {
    heatmap.2(
      x = cor.matrix,
      trace = "none",
      Rowv = NA,
      Colv = NA,
      col = pyCols
    )
  }
  return(object)
}

# Documentation
###############
#' @export
#'
GeneCorMatrix <- function(
  object,
  cor.genes = NULL,
  cell.inds = NULL,
  do.k = FALSE,
  k.seed = 1,
  k.num = 4,
  vis.low = (-1),
  vis.high = 1,
  vis.one = 0.8,
  pcs.use = 1:3,
  col.use = pyCols
) {
  cor.genes <- SetIfNull(x = cor.genes, default = object@var.genes)
  cell.inds <- SetIfNull(x = cell.inds, default = colnames(x = object@data))
  cor.genes <- cor.genes[cor.genes %in% rownames(x = object@data)]
  data.cor <- object@data[cor.genes, cell.inds]
  cor.matrix <- cor(x = t(x = data.cor))
  set.seed(seed = k.seed)
  kmeans.cor <- kmeans(x = cor.matrix, centers = k.num)
  cor.matrix <- cor.matrix[order(kmeans.cor$cluster), order(kmeans.cor$cluster)]
  kmeans.names <- rownames(x = cor.matrix)
  row.annot <- data.frame(
    cbind(
      kmeans.cor$cluster[kmeans.names],
      object@pca.x[kmeans.names, pcs.use]
    )
  )
  colnames(x = row.annot) <- c("K", paste0("PC", pcs.use))
  cor.matrix[cor.matrix == 1] <- vis.one
  cor.matrix <- minmax(data = cor.matrix, min = vis.low, max = vis.high)
  object@kmeans.gene <- list(kmeans.cor)
  if (do.k) {
    aheatmap(
      x = cor.matrix,
      col = col.use,
      Rowv = NA,
      Colv = NA,
      annRow = row.annot
    )
  } else {
    aheatmap(
      x = cor.matrix,
      col = col.use,
      annRow = row.annot
    )
  }
  return(object)
}

# Documentation
###############
#' @export
#'
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
###############
#' @export
setMethod(
  f = "show",
  signature = "seurat",
  definition = function(object) {
    cat(
      "An object of class",
      class(object),
      "in project",
      object@project.name,
      "\n",
      nrow(x = object@data),
      "genes across",
      ncol(x = object@data),
      "samples.\n"
    )
    invisible(x = NULL)
  }
)

# Documentation
###############
#' @export
#'
RemovePC <- function(object, pcs.remove, use.full = FALSE, ...) {
  data.old <- object@data
  pcs.use <- anotinb(x = 1:ncol(x = object@pca.obj[[1]]$rotation), y = pcs.remove)
  if (use.full) {
    data.x <- as.matrix(
      x = object@pca.x.full[, ainb(
        a = pcs.use,
        b = 1:ncol(x = object@pca.x.full)
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
###############
GeneScorePlot <- function(
  object,
  gene1,
  score.name,
  cell.ids = NULL,
  col.use = NULL,
  nrpoints.use = Inf,
  pch.use = 16,
  cex.use = 2,
  ...
) {
  cell.ids <- SetIfNull(x = cell.ids, default = colnames(x = object@data))
  g1 <- as.numeric(x = object@data[gene1, cell.ids])
  my.score <- retreiveScore(object, score.name)
  s1 <- as.numeric(x = my.score[cell.ids])
  col.use <- SetIfNull(as.numeric(x = as.factor(x = object@ident[cell.ids])))
  gene.cor <- round(x = cor(x = g1, y = s1), digits = 2)
  smoothScatter(
    x = g1,
    y = s1,
    xlab = gene1,
    ylab = score.name,
    col = col.use,
    nrpoints = nrpoints.use,
    cex = cex.use,
    main = gene.cor,
    pch = pch.use
  )
}
