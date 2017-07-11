# Find gene terms from Enrichr
#
# @param QueryGene A gene to query on Enrichr
#
# @return An XML document with information on the queried gene
#
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

# Regress out technical effects and cell cycle using regularized Negative Binomial regression
#
# Remove unwanted effects from umi data and set scale.data to Pearson residuals
# Uses mclapply; you can set the number of cores it will use to n with command options(mc.cores = n)
#
# @param object Seurat object
# @param latent.vars effects to regress out
# @param genes.regress gene to run regression for (default is all genes)
# @param pr.clip.range numeric of length two specifying the min and max values the results will be clipped to
#
# @return Returns Seurat object with the scale.data (object@scale.data) genes returning the residuals fromthe regression model
#
# @import Matrix
# @importFrom MASS theta.ml negative.binomial
# @import parallel
#
# @export
#
RegressOutNBreg <- function(
  object,
  latent.vars,
  genes.regress = NULL,
  pr.clip.range = c(-30, 30),
  min.theta = 0.01
) {
  genes.regress <- SetIfNull(x = genes.regress, default = rownames(x = object@data))
  genes.regress <- ainb(genes.regress,rownames(x = object@data))
  cm <- object@raw.data[genes.regress, colnames(x = object@data), drop = FALSE]
  latent.data <- FetchData(object = object, vars.all = latent.vars)
  cat(sprintf('Regressing out %s for %d genes\n', paste(latent.vars), length(x = genes.regress)))
  theta.fit <- theta.reg(cm = cm, latent.data = latent.data, min.theta = 0.01, bin.size = 128)
  print('Second run NB regression with fixed theta')
  bin.size <- 128
  bin.ind <- ceiling(1:length(genes.regress)/bin.size)
  max.bin <- max(bin.ind)
  pb <- txtProgressBar(min = 0, max = max.bin, style = 3)
  pr <- c()
  for (i in 1:max.bin) {
    genes.bin.regress <- genes.regress[bin.ind == i]
    bin.pr.lst <- parallel::mclapply(
      X = genes.bin.regress,
      FUN = function(j) {
        fit <- 0
        try(
          expr = fit <- glm(
            cm[j, ] ~ .,
            data = latent.data,
            family = MASS::negative.binomial(theta = theta.fit[j])
          ),
          silent=TRUE
        )
        if (class(fit)[1] == 'numeric') {
          message(
            sprintf(
              'glm and family=negative.binomial(theta=%f) failed for gene %s; falling back to scale(log10(y+1))',
              theta.fit[j],
              j
            )
          )
          res <- scale(log10(cm[j, ] + 1))[, 1]
        } else {
          res <- residuals(fit, type = 'pearson')
        }
        return(res)
      }
    )
    pr <- rbind(pr, do.call(rbind, bin.pr.lst))
    setTxtProgressBar(pb, i)
  }
  close(pb)
  dimnames(x = pr) <- dimnames(x = cm)
  pr[pr < pr.clip.range[1]] <- pr.clip.range[1]
  pr[pr > pr.clip.range[2]] <- pr.clip.range[2]
  object@scale.data <- pr
  return(object)
}

# Regress out technical effects and cell cycle using regularized Negative Binomial regression
#
# Remove unwanted effects from umi data and set scale.data to Pearson residuals
# Uses mclapply; you can set the number of cores it will use to n with command options(mc.cores = n)
#
# @param object Seurat object
# @param latent.vars effects to regress out
# @param genes.regress gene to run regression for (default is all genes)
# @param pr.clip.range numeric of length two specifying the min and max values the results will be clipped to
#
# @return Returns Seurat object with the scale.data (object@scale.data) genes returning the residuals from the regression model
#
# @import Matrix
# @importFrom MASS theta.ml negative.binomial
# @import parallel
#
# @export

RegressOutNBreg <- function(
  object,
  latent.vars,
  genes.regress = NULL,
  pr.clip.range = c(-30, 30),
  min.theta = 0.01
) {
  genes.regress <- SetIfNull(x = genes.regress, default = rownames(x = object@data))
  genes.regress <- ainb(genes.regress, rownames(x = object@data))
  cm <- object@raw.data[genes.regress, colnames(x = object@data), drop=FALSE]
  latent.data <- FetchData(boject = object, vars.all = latent.vars)
  bin.size <- 128
  bin.ind <- ceiling(x = 1:length(x = genes.regress) / bin.size)
  max.bin <- max(bin.ind)
  print(paste("Regressing out", latent.vars))
  print('First run Poisson regression (to get initial mean), and estimate theta per gene')
  pb <- txtProgressBar(min = 0, max = max.bin, style = 3)
  theta.estimate <- c()
  for (i in 1:max.bin) {
    genes.bin.regress <- genes.regress[bin.ind == i]
    bin.theta.estimate <- unlist(
      parallel::mclapply(
        X = genes.bin.regress,
        FUN = function(j) {
          as.numeric(
            x = MASS::theta.ml(
              cm[j, ],
              glm(cm[j, ] ~ ., data = latent.data, family=poisson)$fitted
            )
          )
        }
      ),
      use.names = FALSE
    )
    theta.estimate <- c(theta.estimate, bin.theta.estimate)
    setTxtProgressBar(pb, i)
  }
  close(pb)
  UMI.mean <- apply(X = cm, MARGIN = 1, FUN = mean)
  var.estimate <- UMI.mean + (UMI.mean ^ 2) / theta.estimate
  fit <- loess(log10(var.estimate) ~ log10(UMI.mean), span = 0.33)
  theta.fit <- (UMI.mean ^ 2) / (10 ^ fit$fitted - UMI.mean)
  names(x = theta.fit) <- genes.regress
  to.fix <- theta.fit <= min.theta | is.infinite(x = theta.fit)
  if (any(to.fix)) {
    cat(
      'Fitted theta below',
      min.theta,
      'for',
      sum(to.fix),
      'genes, setting them to',
      min.theta,
      '\n'
    )
    theta.fit[to.fix] <- min.theta
  }
  print('Second run NB regression with fixed theta')
  pb <- txtProgressBar(min = 0, max = max.bin, style = 3)
  pr <- c()
  for(i in 1:max.bin) {
    genes.bin.regress <- genes.regress[bin.ind == i]
    bin.pr.lst <- parallel::mclapply(
      X = genes.bin.regress,
      FUN = function(j) {
        fit <- 0
        try(
          fit <- glm(
            cm[j, ] ~ .,
            data = latent.data,
            family=MASS::negative.binomial(theta = theta.fit[j])
          ),
          silent=TRUE
        )
        if (class(fit)[1] == 'numeric') {
          message(
            sprintf(
              'glm and family=negative.binomial(theta=%f) failed for gene %s; falling back to scale(log10(y+1))',
              theta.fit[j],
              j
            )
          )
          res <- scale(x = log10(cm[j, ] + 1))[, 1]
        } else {
          res <- residuals(object = fit, type='pearson')
        }
        return(res)
      }
    )
    pr <- rbind(pr, do.call(rbind, bin.pr.lst))
    setTxtProgressBar(pb, i)
  }
  close(pb)
  dimnames(pr) <- dimnames(cm)
  pr[pr < pr.clip.range[1]] <- pr.clip.range[1]
  pr[pr > pr.clip.range[2]] <- pr.clip.range[2]
  object@scale.data <- r
  return(object)
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
