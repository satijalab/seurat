#' Likelihood ratio test for zero-inflated data
#'
#' Identifies differentially expressed genes between two groups of cells using
#' the LRT model proposed in Mcdavid et al, Bioinformatics, 2011
#'
#' @inheritParams FindMarkers
#' @param object Seurat object
#' @param cells.1 Group 1 cells
#' @param cells.2 Group 2 cells
#' @return Returns a p-value ranked matrix of putative differentially expressed
#' genes.
#'
#' @export
#'
DiffExpTest <- function(
  object,
  cells.1,
  cells.2,
  genes.use = NULL,
  print.bar = TRUE
) {
  genes.use <- set.ifnull(x = genes.use, y = rownames(x = object@data))
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
          diffLRT(
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
#
#'
#' @inheritParams FindMarkers
#' @param object Seurat object
#' @param cells.1 Group 1 cells
#' @param cells.2 Group 2 cells
#'
#' @return Returns a p-value ranked matrix of putative differentially expressed
#' genes.
#'
#' @importFrom MASS glm.nb
#' @importFrom pbapply pbapply
#'
#' @export
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
  genes.use <- set.ifnull(x = genes.use, y = rownames(x = object@data))
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
#'
#' @return Returns a p-value ranked data frame of test results.
#'
#' @export
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
  genes.use <- set.ifnull(x = genes.use, y = rownames(x = object@data))
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
  theta.fit <- theta.reg(
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
        return(de.nb.reg(
          y = to.test.data[j, ],
          theta = theta.fit[j],
          latent.data = to.test,
          com.frac = latent.vars,
          grp.frac = 'NegBinomRegDETest.group'
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

#' Poisson test for UMI-count based data
#'
#' Identifies differentially expressed genes between two groups of cells using
#' a poisson generalized linear model
#
#' @inheritParams FindMarkers
#' @param object Seurat object
#' @param cells.1 Group 1 cells
#' @param cells.2 Group 2 cells
#'
#' @return Returns a p-value ranked matrix of putative differentially expressed
#' genes.
#'
#' @importFrom pbapply pbapply
#'
#' @export
#'
PoissonDETest <- function(
  object,
  cells.1,
  cells.2,
  genes.use = NULL,
  latent.vars = NULL,
  print.bar = TRUE
) {
  genes.use <- set.ifnull(x = genes.use, y = rownames(x = object@data))
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
TobitTest <- function(
  object,
  cells.1,
  cells.2,
  genes.use = NULL,
  print.bar = TRUE
) {
  genes.use <- set.ifnull(x = genes.use, y = rownames(x = object@data))
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
#' @import ROCR
#'
#' @export
#'
MarkerTest <- function(
  object,
  cells.1,
  cells.2,
  genes.use = NULL,
  print.bar = TRUE
) {
  genes.use <- set.ifnull(x = genes.use, y = rownames(x = object@data))
  to.return <- marker.auc.test(
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
#' @importFrom pbapply pblapply
#'
#' @export
#'
DiffTTest <- function(
  object,
  cells.1,
  cells.2,
  genes.use = NULL,
  print.bar = TRUE
) {
  genes.use <- set.ifnull(x = genes.use, y = rownames(x = object@data))
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
