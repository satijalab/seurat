globalVariables(names = 'i', package = 'Seurat', add = TRUE)
# Regress out technical effects and cell cycle
#
# Remove unwanted effects from scale.data
#
# @keywords internal
# @param object Seurat object
# @param vars.to.regress effects to regress out
# @param genes.regress gene to run regression for (default is all genes)
# @param model.use Use a linear model or generalized linear model (poisson, negative binomial) for the regression. Options are 'linear' (default), 'poisson', and 'negbinom'
# @param use.umi Regress on UMI count data. Default is FALSE for linear modeling, but automatically set to TRUE if model.use is 'negbinom' or 'poisson'
# @param display.progress display progress bar for regression procedure.
# @param do.par use parallel processing for regressing out variables faster.
# If set to TRUE, will use half of the machines available cores (FALSE by default)
# @param num.cores If do.par = TRUE, specify the number of cores to use.
#
# @return Returns the residuals from the regression model
#
#' @import Matrix
#' @import doSNOW
#' @importFrom stats as.formula lm residuals glm
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @importFrom foreach foreach %dopar%
#
RegressOutResid <- function(
  object,
  vars.to.regress,
  genes.regress = NULL,
  model.use = 'linear',
  use.umi = FALSE,
  display.progress = TRUE,
  do.par = FALSE,
  num.cores = 1
) {
  possible.models <- c("linear", "poisson", "negbinom")
  if (!model.use %in% possible.models) {
    stop(
      paste0(
        model.use,
        " is not a valid model. Please use one the following: ",
        paste0(possible.models, collapse = ", "),
        "."
      )
    )
  }
  genes.regress <- SetIfNull(x = genes.regress, default = rownames(x = object@data))
  genes.regress <- intersect(x = genes.regress, y = rownames(x = object@data))
  latent.data <- FetchData(object = object, vars.all = vars.to.regress)
  bin.size <- ifelse(test = model.use == 'negbinom', yes = 5, no = 100)
  bin.ind <- ceiling(x = 1:length(x = genes.regress) / bin.size)
  max.bin <- max(bin.ind)
  if (display.progress) {
    message(paste("Regressing out:", paste(vars.to.regress, collapse = ", ")))
    pb <- txtProgressBar(min = 0, max = max.bin, style = 3, file = stderr())
  }
  data.resid <- c()
  data.use <- object@data[genes.regress, , drop = FALSE];
  if (model.use != "linear") {
    use.umi <- TRUE
  }
  if (use.umi) {
    data.use <- object@raw.data[genes.regress, object@cell.names, drop = FALSE]
  }
  # input checking for parallel options
  if (do.par) {
    if (num.cores == 1) {
      num.cores <- detectCores() / 2
    } else if (num.cores > detectCores()) {
      num.cores <- detectCores() - 1
      warning(paste0("num.cores set greater than number of available cores(", detectCores(), "). Setting num.cores to ", num.cores, "."))
    }
  } else if (num.cores != 1) {
    num.cores <- 1
    warning("For parallel processing, please set do.par to TRUE.")
  }
  cl <- parallel::makeCluster(num.cores)#, outfile = "")
  # using doSNOW library because it supports progress bar update
  registerDoSNOW(cl)
  opts <- list()
  if (display.progress) {
    # define progress bar function
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress = progress)
    time_elapsed <- Sys.time()
  }
  data.resid <- foreach(i = 1:max.bin, .combine = "c", .options.snow = opts) %dopar% {
    genes.bin.regress <- rownames(x = data.use)[bin.ind == i]
    gene.expr <- as.matrix(x = data.use[genes.bin.regress, , drop = FALSE])
    new.data <- sapply(
      X = genes.bin.regress,
      FUN = function(x) {
        regression.mat <- cbind(latent.data, gene.expr[x,])
        colnames(x = regression.mat) <- c(colnames(x = latent.data), "GENE")
        fmla <- as.formula(
          object = paste0(
            "GENE ",
            " ~ ",
            paste(vars.to.regress, collapse = "+")
          )
        )
        resid <- switch(
          EXPR = model.use,
          'linear' = lm(formula = fmla, data = regression.mat)$residuals,
          'poisson' = residuals(
            object = glm(
              formula = fmla,
              data = regression.mat,
              family = "poisson"
            ),
            type = 'pearson'
          ),
          'negbinom' = NBResiduals(
            fmla = fmla,
            regression.mat = regression.mat,
            gene = x,
            return.mode = TRUE
          )
        )
        if (!is.list(x = resid)) {
          resid <- list('resid' = resid, 'mode' = character(length = length(x = resid)))
        }
        return(resid)
      }
    )
    new.data.resid <- new.data[seq.int(from = 1, to = length(x = new.data), by = 2)]
    new.data.resid <- as.data.frame(x = new.data.resid)
    colnames(x = new.data.resid) <- genes.bin.regress
    new.data.resid <- as.matrix(x = new.data.resid)
    new.data.mode <- unlist(x = new.data[seq.int(from = 2, to = length(x = new.data), by = 2)])
    names(x = new.data.mode) <- genes.bin.regress
    new.data <- list('resid' = new.data.resid, 'mode' = new.data.mode)
    return(new.data)
  }
  if (display.progress) {
    time_elapsed <- Sys.time() - time_elapsed
    cat(paste("\nTime Elapsed: ",time_elapsed, units(time_elapsed)))
    close(pb)
  }
  stopCluster(cl)
  modes <- unlist(x = data.resid[seq.int(from = 2, to = length(x = data.resid), by = 2)])
  modes <- modes[modes == 'scale']
  names(x = modes) <- gsub(
    pattern = 'mode.',
    replacement = '',
    x = names(x = modes),
    fixed = TRUE
  )
  data.resid <- data.resid[seq.int(from = 1, to = length(x = data.resid), by = 2)]
  data.resid <- as.matrix(x = as.data.frame(x = data.resid))
  data.resid <- t(x = data.resid)
  if (length(x = modes)) {
    message(
      "The following genes failed with glm.nb, and fell back to scale(log(y+1))\n\t",
      paste(names(x = modes), collapse = ', ')
    )
  }
  rownames(x = data.resid) <- genes.regress
  suppressWarnings(expr = gc(verbose = FALSE))
  if (use.umi) {
    data.resid <- log1p(
      x = sweep(
        x = data.resid,
        MARGIN = 1,
        STATS = apply(X = data.resid, MARGIN = 1, FUN = min),
        FUN = "-"
      )
    )
  }
  return(data.resid)
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
#' @import Matrix
#' @import parallel
#' @importFrom stats glm residuals
#' @importFrom MASS theta.ml negative.binomial
#' @importFrom utils txtProgressBar setTxtProgressBar
#
RegressOutNB <- function(
  object,
  latent.vars,
  genes.regress = NULL,
  pr.clip.range = c(-30, 30),
  min.theta = 0.01
) {
  genes.regress <- SetIfNull(x = genes.regress, default = rownames(x = object@data))
  genes.regress <- intersect(x = genes.regress, y = rownames(x = object@data))
  cm <- object@raw.data[genes.regress, colnames(x = object@data), drop = FALSE]
  latent.data <- FetchData(object = object, vars.all = latent.vars)
  message(sprintf('Regressing out %s for %d genes\n', paste(latent.vars), length(x = genes.regress)))
  theta.fit <- RegularizedTheta(cm = cm, latent.data = latent.data, min.theta = 0.01, bin.size = 128)
  message('Second run NB regression with fixed theta')
  bin.size <- 128
  bin.ind <- ceiling(1:length(genes.regress)/bin.size)
  max.bin <- max(bin.ind)
  pb <- txtProgressBar(min = 0, max = max.bin, style = 3, file = stderr())
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
# If n.genes.step1 is set, only a (somewhat-random) subset of genes is used for estimating theta.
#
# @param object Seurat object
# @param latent.vars effects to regress out (character vector)
# @param n.genes.step1 number of genes to use when estimating theta (default uses all genes)
# @param genes.regress gene to run regression for (default uses all genes)
# @param res.clip.range numeric of length two specifying the min and max values the results will be clipped to
# @param min.theta minimum theta to use in NB regression
# @param residual.type string specifying the type of residual used (default is pearson)
# @param bin.size number of genes to put in each bin (to show progress)
# @param use.stored.theta skip the first step and use the fitted thetas from a previous run
# @param min.cells only use genes that have been observed in at least this many cells
#
# @return Returns Seurat object with the scale.data (object@scale.data) genes returning the residuals from the regression model
#
#' @import Matrix
#' @import parallel
#' @importFrom MASS theta.ml negative.binomial
#' @importFrom stats glm loess residuals approx
#' @importFrom utils txtProgressBar setTxtProgressBar
#
RegressOutNBreg <- function(
  object,
  latent.vars,
  n.genes.step1 = NULL,
  genes.regress = NULL,
  res.clip.range = c(-30, 30),
  min.theta = 0.01,
  residual.type = 'pearson',
  bin.size = 128,
  use.stored.theta = FALSE,
  min.cells = 3
) {
  # in the first step we use all genes, except if n.genes.step1 has been set
  cm <- Matrix(object@raw.data[, colnames(x = object@data), drop = FALSE])
  gene.observed <- rowSums(cm > 0)
  genes.regress <- SetIfNull(x = genes.regress, default = rownames(x = cm))
  genes.regress <- intersect(x = genes.regress, y = rownames(x = cm)[gene.observed >= min.cells])
  genes.step1 <- rownames(cm)[gene.observed >= min.cells]
  if (!is.null(n.genes.step1)) {
    # density-sample genes to speed up the first step
    raw.mean <- log10(rowMeans(cm[genes.step1, ]))
    raw.det.rate <- rowMeans(cm[genes.step1, ] > 0)
    dens <- apply(
      X = cbind(raw.mean, raw.det.rate),
      MARGIN = 2,
      FUN = function(y) {
        y.dens <- density(x = y, bw = 'nrd', adjust = 1)
        ret <- approx(x = y.dens$x, y = y.dens$y, xout = y)$y
        return(ret / sum(ret))
      }
    )
    sampling.prob <- 1 / apply(X = dens, MARGIN = 1, FUN = min)
    genes.step1 <- sample(x = genes.step1, size = n.genes.step1, prob = sampling.prob)
  }
  latent.data <- FetchData(object = object, vars.all = latent.vars)
  bin.ind <- ceiling(x = 1:length(x = genes.step1) / bin.size)
  max.bin <- max(bin.ind)
  message(paste("Regressing out", paste(latent.vars, collapse = ' ')))
  if (use.stored.theta) {
    message('Using previously fitted theta values for NB regression')
    theta.fit <- object@misc[['NBreg.theta.fit']]
  } else {
    message('First step: Poisson regression (to get initial mean), and estimate theta per gene')
    message('Using ', length(x = genes.step1), ' genes')
    pb <- txtProgressBar(min = 0, max = max.bin, style = 3)
    theta.estimate <- c()
    for (i in 1:max.bin) {
      genes.bin.regress <- genes.step1[bin.ind == i]
      bin.theta.estimate <- unlist(
        parallel::mclapply(
          X = genes.bin.regress,
          FUN = function(j) {
            as.numeric(
              x = MASS::theta.ml(
                as.numeric(x = unlist(x = cm[j, ])),
                glm(as.numeric(x = unlist(x = cm[j, ])) ~ ., data = latent.data, family = poisson)$fitted
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
    names(theta.estimate) <- genes.step1
    UMI.mean <- rowMeans(x = cm[genes.step1, ])
    var.estimate <- UMI.mean + (UMI.mean ^ 2) / theta.estimate
    fit <- loess(log10(var.estimate) ~ log10(UMI.mean), span = 0.33)
    genes.regress.mean <- rowMeans(x = cm[genes.regress, ])
    var.fit.log10 <- predict(fit, log10(genes.regress.mean))
    theta.fit <- (genes.regress.mean ^ 2) / (10 ^ var.fit.log10 - genes.regress.mean)
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
    # save theta estimate and fitted theta in object
    object@misc <- as.list(object@misc)
    object@misc[['NBreg.theta.estimate']] <- theta.estimate
    object@misc[['NBreg.theta.fit']] <- theta.fit
  }
  message('Second step: NB regression with fixed theta for ', length(x = genes.regress), ' genes')
  bin.ind <- ceiling(x = 1:length(x = genes.regress) / bin.size)
  max.bin <- max(bin.ind)
  pb <- txtProgressBar(min = 0, max = max.bin, style = 3, file = stderr())
  res <- c()
  for (i in 1:max.bin) {
    genes.bin.regress <- genes.regress[bin.ind == i]
    names(genes.bin.regress) <- genes.bin.regress
    bin.res.lst <- parallel::mclapply(
      X = genes.bin.regress,
      FUN = function(j) {
        fit <- 0
        try(
          fit <- glm(
            as.numeric(x = unlist(x = cm[j, ])) ~ .,
            data = latent.data,
            family = MASS::negative.binomial(theta = theta.fit[j])
          ),
          silent = TRUE
        )
        if (class(fit)[1] == 'numeric') {
          message <-
            sprintf(
              'glm and family=negative.binomial(theta=%f) failed for gene %s; falling back to scale(log10(y+1))',
              theta.fit[j],
              j
            )
          res <- scale(x = log10(as.numeric(x = unlist(x = cm[j, ])) + 1))[, 1]
        } else {
          message <- NULL
          res <- residuals(object = fit, type = residual.type)
        }
        return(list(res = res, message = message))
      }
    )
    # Print message to keep track of the genes for which glm failed to converge
    message <- unlist(x = lapply(X = bin.res.lst, FUN = function(l) { return(l$message) }), use.names = FALSE)
    if (!is.null(x = message)) {
      message(paste(message, collapse = "\n"))
    }
    bin.res.lst <- lapply(X = bin.res.lst, FUN = function(l) { return(l$res) })
    res <- rbind(res, do.call(rbind, bin.res.lst))
    setTxtProgressBar(pb, i)
  }
  close(pb)
  dimnames(x = res) <- list(genes.regress, colnames(x = cm))
  res[res < res.clip.range[1]] <- res.clip.range[1]
  res[res > res.clip.range[2]] <- res.clip.range[2]
  object@scale.data <- res
  return(object)
}

# Normalize raw data
#
# Normalize count data per cell and transform to centered log ratio
#
# @param data Matrix with the raw count data
# @param custom_function A custom normalization function
# @parm across Which way to we normalize? Choose form 'cells' or 'genes'
#
# @return Returns a matrix with the custom normalization
#
# @import Matrix
#
CustomNormalize <- function(data, custom_function, across) {
  if (class(x = data) == "data.frame") {
    data <- as.matrix(x = data)
  }
  if (class(x = data) != "dgCMatrix") {
    data <- as(object = data, Class = "dgCMatrix")
  }
  margin <- switch(
    EXPR = across,
    'cells' = 2,
    'genes' = 1,
    stop("'across' must be either 'cells' or 'genes'")
  )
  norm.data <- apply(
    X = data,
    MARGIN = margin,
    FUN = custom_function)
  if (margin == 1) {
    norm.data = t(x = norm.data)
  }
  colnames(x = norm.data) <- colnames(x = data)
  rownames(x = norm.data) <- rownames(x = data)
  return(norm.data)
}
