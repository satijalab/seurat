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

#internal function to run Tobit DE test
TobitDiffExpTest <- function(data1, data2, mygenes, print.bar) {
  p_val <- unlist(x = lapply(
    X = mygenes,
    FUN = function(x) {
      return(DifferentialTobit(
        x1 = as.numeric(x = data1[x, ]),
        x2 = as.numeric(x = data2[x, ])
      ))}
  ))
  p_val[is.na(x = p_val)] <- 1
  if (print.bar) {
    iterate.fxn <- pblapply
  } else {
    iterate.fxn <- lapply
  }
  toRet <- data.frame(p_val, row.names = mygenes)
  return(toRet)
}

#internal function to run Tobit DE test
#
#' @importFrom stats pchisq logLik
#
DifferentialTobit <- function(x1, x2, lower = 1, upper = Inf) {
  my.df <- data.frame(
    c(x1, x2),
    c(rep(x = 0, length(x = x1)), rep(x = 1, length(x = x2)))
  )
  colnames(x = my.df) <- c("Expression", "Stat")
  #model.v1=vgam(Expression~1,family = tobit(Lower = lower,Upper = upper),data = my.df)
  model.v1 <- TobitFitter(
    x = my.df,
    modelFormulaStr = "Expression~1",
    lower = lower,
    upper = upper
  )
  #model.v2=vgam(Expression~Stat+1,family = tobit(Lower = lower,Upper = upper),data = my.df)
  model.v2 <- TobitFitter(
    x = my.df,
    modelFormulaStr = "Expression~Stat+1",
    lower = lower,
    upper = upper
  )
  # if (is.null(x = model.v1) == FALSE && is.null(x = model.v2) == FALSE) {
  if (! is.null(x = model.v1) && ! is.null(x = model.v2)) {
    p <- pchisq(
      q = 2 * (logLik(object = model.v2) - logLik(object = model.v1)),
      df = 1,
      lower.tail = FALSE
    )
  } else {
    p <- 1
  }
  return(p)
}

#internal function to run Tobit DE test
#credit to Cole Trapnell for this
#
#' @importFrom VGAM vgam tobit
#' @importFrom stats as.formula
#
TobitFitter <- function(x, modelFormulaStr, lower = 1, upper = Inf){
  tryCatch(
    expr = return(suppressWarnings(expr = vgam(
      formula = as.formula(object = modelFormulaStr),
      family = tobit(Lower = lower, Upper = upper),
      data = x
    ))),
    #warning = function(w) { FM_fit },
    error = function(e) { NULL }
  )
}

# internal function to calculate AUC values
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
  if (print.bar) {
    iterate.fxn <- pblapply
  } else {
    iterate.fxn <- lapply
  }
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

# internal function to calculate AUC values
#' @importFrom ROCR prediction performance
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
  if (class(x = fit2)[1] == 'numeric' | class(x = fit4)[1] == 'numeric') {
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
