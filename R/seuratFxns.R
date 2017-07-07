
situ3d <- function(data, label = NULL, ...) {
  # Call Seurat function to get the in situ values out.
  exp.1 <- data
  exp.1 <- (exp.1 - min(exp.1)) / (max(exp.1) - min(exp.1))
  # Reformat them into an expression matrix as expected by the plotting function
  expression.matrix <- data.frame(matrix(data = exp.1, nrow = 8, ncol = 8))
  rownames(x = expression.matrix) <- c(
    "24-30",
    "17-23",
    "13-16",
    "9-12",
    "7-8",
    "5-6",
    "3-4",
    "1-2"
  )
  names(x = expression.matrix) <- c(
    "1-4",
    "5-8",
    "9-12",
    "13-16",
    "17-20",
    "21-24",
    "25-28",
    "29-32"
  )
  # Call the plotting function.
  zf.insitu.side(expression.matrix)
  par3d(windowRect = c(0, 0, 800, 800))
  # Label or not and then set the view.
  if (! is.null(x = label)) {
    text3d(x = 0, y = 0, z = 1.5, text = label, cex = 3)
  }
  view3d(zoom = .75, theta = 0, phi = -90, fov = 0)
}

aucFxn <- function(preds, truth, do.plot = FALSE, lab.main = "", ...) {
  pred.use <- ROCR::prediction(
    predictions = preds,
    labels = truth,
    label.ordering = 0:1
  )
  perf.use1 <- ROCR::performance(
    prediction.obj = pred.use,
    measure = "tpr",
    x.measure = "fpr"
  )
  perf.use <- ROCR::performance(prediction.obj = pred.use, measure = "auc")
  auc.use <- round(x = perf.use@y.values[[1]], digits = 3)
  if (do.plot) {
    plot(
      x = perf.use1@x.values[[1]],
      y = perf.use1@y.values[[1]],
      type = "l",
      xlab = "FPR",
      ylab = "TPR",
      main = paste(lab.main, auc.use)
    )
  }
  return(auc.use)
}

#' @export
humpMean <- function(x, min = 0) {
  return(mean(x = x[x > min]))
}

#' @export
humpVar <- function(x, min = 0) {
  return(var(x = x[x > min]))
}

debugdmvnorm <- function (x, mean = rep(0, p), sigma = diag(p), log = FALSE) {
  if (is.vector(x = x)) {
    x <- matrix(x, ncol = length(x))
  }
  p <- ncol(x = x)
  dec <- tryCatch(chol(x = sigma), error = function(e) e)
  tmp <- backsolve(r = dec, x = t(x) - mean, transpose = TRUE)
  rss <- colSums(x = tmp ^ 2)
  logretval <- -(log(x = diag(x = dec))) - 0.5 * p * log(x = 2 * pi) - 0.5 * rss
  names(x = logretval) <- names(x = x)
  return(logretval)
}

slimdmvnorm <- function (x, mean = rep(0, p), sigma = diag(p), log = FALSE) {
  x <- matrix(data = x, ncol = length(x = x))
  p <- ncol(x = x)
  dec <- tryCatch(chol(x = sigma), error = function(e) e)
  tmp <- backsolve(r = dec, t(x = x) - mean, transpose = TRUE)
  rss <- colSums(tmp ^ 2)
  logretval <- -sum(log(x = diag(x = dec))) - 0.5 * p * log(x = 2 * pi) - 0.5 * rss
  names(x = logretval) <- rownames(x = x)
  return(logretval)
}

slimdmvnorm_nosum <- function (
  x,
  mean = rep(x = 0, p),
  sigma = diag(x = p),
  log = FALSE
) {
  x <- matrix(data = x, ncol = length(x = x))
  p <- ncol(x = x)
  dec <- tryCatch(chol(x = sigma), error = function(e) e)
  tmp <- backsolve(r = dec, t(x = x) - mean, transpose = TRUE)
  rss <- colSums(x = tmp ^ 2)
  logretval <- -(log(x = diag(x = dec))) - 0.5 * p * log(x = 2 * pi) - 0.5 * rss
  names(x = logretval) <- rownames(x = x)
  return(logretval)
}

compareBins <- function(object, cell.use, bin.1, bin.2, bins.mu, bins.cov) {
  num.genes <- ncol(x = bins.mu)
  genes.use <- colnames(x = bins.mu)
  to.par <- floor(x = sqrt(x = num.genes)) + 1
  par(mfrow = c(to.par, to.par))
  data.use <- object@imputed[genes.use, cell.use]
  names(x = data.use) <- genes.use
  lik.diff <- sort(x = unlist(x = lapply(
    X = genes.use,
    FUN = function(g) {
      return(dnorm(
        x = data.use[g],
        mean = bins.mu[bin.1, g],
        sd = sqrt(x = bins.cov[[bin.1]][g, g])
      ) / dnorm(
        x = data.use[g],
        mean = bins.mu[bin.2, g],
        sd = sqrt(x = bins.cov[[bin.2]][g, g])
      ))
    }))
  )
  for (g in names(x = lik.diff)) {
    plot(
      x = 0,
      y = 0,
      type = "n",
      xlim = c(0, 8),
      ylim = c(0, 2),
      main = paste(g, round(x = lik.diff[g], digits = 2))
    )
    lines(
      x = density(
        x = rnorm(
          n = 10000,
          mean = bins.mu[bin.1, g],
          sd = sqrt(x = bins.cov[[bin.1]][g, g])
        )
      ),
      lwd = 2,
      col = "black"
    )
    lines(
      x = density(
        x = rnorm(
          n = 10000,
          mean = bins.mu[bin.2, g],sd = sqrt(x = bins.cov[[bin.2]][g, g])
        )
      ),
      lwd = 2,
      col = "red"
    )
    points(x = data.use[g], y = 0, cex = 1.6, col = "darkgreen", pch = 16)
  }
  #rp()
}

#' @export
subr <- function(data, code) {
  return(data[grep(pattern = code, x = rownames(x = data)), ])
}

#' @export
cv <- function(x) {
  return(sd(x = x) / mean(x = x))
}

#' @export
humpCt <- function(x, min = 0) {
  return(length(x = x[x > min]))
}

#' @export
log_add <- function(x) {
  mpi <- max(x)
  return(mpi + log(x = sum(exp(x = x - mpi))))
}

#' @export
minusr <- function(data, code) {
  matchCode <- rownames(x = data)[grep(pattern = code, x = rownames(x = data))]
  toIgnore <- which(x = rownames(x = data) %in% matchCode)
  if (length(x = toIgnore) == 0) {
    return(data)
  }
  toRet <- data.frame(data[-toIgnore, ])
  rownames(x = toRet) <- rownames(x = data)[-toIgnore]
  colnames(x = toRet) <- colnames(x = data)
  return(toRet)
}

#' @export
minusc <- function(data, code) {
  matchCode <- colnames(x = data)[grep(pattern = code, colnames(x = data))]
  toIgnore <- which(x = colnames(x = data) %in% matchCode)
  if (length(x = toIgnore) == 0) {
    return(data)
  }
  return(data[, -toIgnore])
}


#' @export
ainb <- function(a, b) {
  return(a[a %in% b])
}

#' @export
meanNormFunction <- function(data, myfuncX, myfuncY, nBin = 20) {
  data_x <- apply(X = data, MARGIN = 1, FUN = myfuncX)
  data_y <- apply(X = data, MARGIN = 1, FUN = myfuncY)
  data_x_bin <- cut(x = data_x, breaks = nBin)
  names(x = data_x_bin) <- names(x = data_x)
  mean_y <- tapply(X = data_y, INDEX = data_x_bin, FUN = mean)
  sd_y <- tapply(X = data_y, INDEX = data_x_bin, FUN = sd)
  return((data_y-mean_y[as.numeric(x = data_x_bin)]) / sd_y[as.numeric(x = data_x_bin)])
}

shift.cell <- function(bin, x, y) {
  bin.y <- (bin - 1) %/% 8 + 1
  bin.x <- (bin - 1) %% 8 + 1
  new.x <- minmax(data = bin.x + x, min = 1, max = 8)
  new.y <- minmax(data = bin.y + y, min = 1, max = 8)
  new.bin <- 8 * (new.y - 1) + new.x
  return(new.bin)
}

empP <- function(x, nullval) {
  return(sum(nullval > x) / length(x = nullval))
}

neighbor.cells <- function(bin) {
  return(unique(x = c(
    bin,
    shift.cell(bin = bin, x = 0, y = 1),
    shift.cell(bin = bin, x = 1, y = 0),
    shift.cell(bin = bin, x = -1, y = 0),
    shift.cell(bin = bin, x = 0, y = -1)
  )))
}

all.neighbor.cells <- function(bin, dist = 1) {
  all.comb <- expand.grid(rep(x = list(-dist:dist), 2))
  return(unique(x = unlist(x = lapply(
    X = 1:nrow(x = all.comb),
    FUN = function(x) {
    return(shift.cell(bin = bin, x = all.comb[x, 1], y = all.comb[x, 2]))
  }))))
}

#' @export
no.legend.title <- theme(legend.title = element_blank())
#' @export
gg.legend.text <- function(x = 12, y = "bold") {
  return(theme(legend.text = element_text(size = x, face = y)))
}
#' @export
gg.legend.pts <- function(x = 6) {
  return(guides(colour = guide_legend(override.aes = list(size = x))))
}
#' @export
gg.xax <- function(x = 16, y = "#990000", z = "bold", x2 = 12) {
  return(theme(
    axis.title.x = element_text(face = z, colour = y, size = x),
    axis.text.x = element_text(angle = 90, vjust = 0.5, size = x2)
  ))
}

#' @export
gg.yax <- function(x = 16, y = "#990000", z = "bold", x2 = 12) {
  return(theme(
    axis.title.y = element_text(face = z, colour = y, size = x),
    axis.text.y = element_text(angle = 90, vjust = 0.5, size = x2)
  ))
}

#' @export
sub.string <- function(x, s1, s2) {
  return(unlist(x = lapply(
    X = x,
    FUN = function(y) {
      return(gsub(pattern = s1, replacement = s2, x = y))
    }
  )))
}

#' @export
jackStrawF <- function(prop = 0.1, myR1, myR2 = 3, data = smD) {
  randGenes <- sample(x = rownames(x = data), size = nrow(x = data) * prop)
  smD.mod <- data
  smD.mod[randGenes, ] <- shuffleMatRow(x = data[randGenes, ])
  fmd.pca <- prcomp(x = smD.mod)
  fmd.x <- fmd.pca$x
  fmd.rot <- fmd.pca$rotation
  fakeF <- unlist(x = lapply(
    X = randGenes,
    FUN = jackF,
    r1 = myR1,
    r2 = myR2,
    x = fmd.x,
    rot = fmd.rot
  ))
}

#' @export
jackF <- function(gene, r1 = 1,r2 = 2, x = md.x, rot = md.rot) {
  if (r2 == 1) { #assuming r1, r2=1
    mod.x <- x[, r1]
    mod.x[gene] <- 0
    return(var.test(
      x = (x[, r1] %*% t(x = rot[, r1])),
      y = (mod.x %*% t(x = rot[, r1]))
    )$statistic)
  }
  mod.x <- x[, 1:r2]
  mod.x[gene, r1:r2] <- rep(x = 0, r2 - r1 + 1)
  return(var.test(
    x = (x[, 1:r2] %*% t(x = rot[, 1:r2])),
    y = (mod.x[, 1:r2] %*% t(x = rot[, 1:r2]))
  )$statistic)
}

#' @export
shuffleMatRow <- function(x) {
  x2 <- x
  x2 <- t(x = x)
  ind <- order(c(col(x = x2)), runif(n = length(x = x2)))
  x2 <- matrix(
    data = x2[ind],
    nrow = nrow(x = x),
    ncol = ncol(x = x),
    byrow = TRUE
  )
  return(x2)
}

#' @export
logMeanMinus <- function(x) {
  return(log(x = mean(x = exp(x = as.numeric(x = x)) - 1) + 1))
}
#' @export
logVarMinus <- function(x) {
  return(mean(x = var(x = as.numeric(x = x)) - 1) + 1)
}

#' @export
logVarMinus2 <- function(x) {
  return(var(x = exp(x = as.numeric(x = x)) - 1) + 1)
}

quickRNAHuman <- function(x) {
  dataFile <- paste0("~/big/", x, "/summary/", x,".expMatrix.txt")
  data <- log(x = read.table(file = dataFile, header = TRUE, sep = "\t")[, -1] + 1)
  data <- subc(data = data, code = "rsem")
  return(data)
}

#' @export
expVar <- function(x) {
  return(log(x = var(x = exp(x = x) - 1) + 1))
}

#' @export
expSD <- function(x) {
  return(log(x = sd(x = exp(x = x) - 1) + 1))
}

#' @export
expMean <- function(x) {
  return(log(x = mean(x = exp(x = x) - 1) + 1))
}

#' @export
normal.sample <- function(x) {
  return(rnorm(n = 10000, mean = mean(x = x), sd = sd(x = x)))
}

fetch.closest <- function(bin, all.centroids, num.cell) {
  bin.y <- (bin - 1) %/% 8 + 1
  bin.x <- (bin - 1) %% 8 + 1
  all.centroids <- rbind(all.centroids, c(bin.x, bin.y))
  all.dist <- as.matrix(x = dist(x = all.centroids))
  return(names(x = sort(x = all.dist[nrow(x = all.dist), ]))[2:(num.cell + 2)])
}

fetch.mincells <- function(bin, cells.max, min.cells) {
  for (i in 1:5) {
    my.names <- names(x = ainb(
      a = cells.max,
      b = all.neighbor.cells(bin = bin, dist = i)
    ))
    if (length(x = my.names) > min.cells) {
      break
    }
  }
  return(my.names)
}

#' @export
cell.centroid <- function(cell.probs) {
  centroid.x <- round(x = sum(sapply(
    X = 1:64,
    FUN = function(x) {
      return((x - 1) %% 8 + 1)
    }
  ) * cell.probs))
  centroid.y <- round(x = sum(sapply(
    X = 1:64,
    FUN = function(x) {
      return((x - 1) %/% 8 + 1)
    }
  ) * cell.probs))
  centroid.bin <- 8 * (centroid.y - 1) + centroid.x
  return(centroid.bin)
}

#' @export
cell.centroid.x <- function(cell.probs) {
  centroid.x <- round(x = sum(sapply(
    X = 1:64,
    FUN = function(x) {
      return((x - 1) %% 8 + 1)
    }
  ) * cell.probs))
  return(centroid.x)
}

#' @export
cell.centroid.y <- function(cell.probs) {
  centroid.y <- round(x = sum(sapply(
    X = 1:64,
    FUN = function(x) {
      return((x - 1) %/% 8 + 1)
    }
  ) * cell.probs))
  return(centroid.y)
}

#' @export
exact.cell.centroid <- function(cell.probs) {
  # centroid.x=(sum(sapply(1:64,function(x)(x-1)%%8+1)*cell.probs))
  centroid.x <- cell.centroid.x(cell.probs = cell.probs)
  # centroid.y=(sum(sapply(1:64,function(x)(x-1)%/%8+1)*cell.probs))
  centroid.y <- cell.centroid.y(cell.probs = cell.probs)
  return(c(centroid.x, centroid.y))
}

#' @export
marker.auc.test <- function(data1, data2, mygenes, print.bar = TRUE) {
  myAUC <- unlist(x = lapply(
    X = mygenes,
    FUN = function(x) {
      return(diffAUC(
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
        expMean(
          x = as.numeric(x = data1[x, ])
        ) - expMean(
          x = as.numeric(x = data2[x, ])
        )
      )
    }
  ))
  toRet <- data.frame(cbind(myAUC, avg_diff), row.names <- mygenes)
  toRet <- toRet[rev(x = order(toRet$myAUC)), ]
  return(toRet)
}

#credit to Cole Trapnell for this
tobit_fitter <- function(x, modelFormulaStr, lower = 1, upper = Inf){
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

#' @export
anotinb <- function(x, y) {
  return(x[! x %in% y])
}

diffTobit <- function(x1, x2, lower = 1, upper = Inf) {
  my.df <- data.frame(
    c(x1, x2),
    c(rep(x = 0, length(x = x1)), rep(x = 1, length(x = x2)))
  )
  colnames(x = my.df) <- c("Expression", "Stat")
  #model.v1=vgam(Expression~1,family = tobit(Lower = lower,Upper = upper),data = my.df)
  model.v1 <- tobit_fitter(
    x = my.df,
    modelFormulaStr = "Expression~1",
    lower = lower,
    upper = upper
  )
  #model.v2=vgam(Expression~Stat+1,family = tobit(Lower = lower,Upper = upper),data = my.df)
  model.v2 <- tobit_fitter(
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

TobitDiffExpTest <- function(data1, data2, mygenes, print.bar) {
  p_val <- unlist(x = lapply(
    X = mygenes,
    FUN = function(x) {
      return(diffTobit(
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

NegBinomDiffExpTest <- function(data1, data2, mygenes, print.bar) {
  p_val <- unlist(x = lapply(
    X = mygenes,
    FUN = function(x) {
      return(diffLRT(
        x = as.numeric(x = data1[x, ]),
        y = as.numeric(x = data2[x, ])
      ))
    }
  ))
  p_val[is.na(x = p_val)] <- 1
  if (print.bar) {
    iterate.fxn <- pblapply
  } else {
    iterate.fxn <- lapply
  }
  avg_diff <- unlist(x = iterate.fxn(
    X = mygenes,
    FUN = function(x) {
      return(
        expMean(
          x = as.numeric(x = data1[x, ])
        ) - expMean(
          x = as.numeric(x = data2[x, ])
        )
      )
    }
  ))
  toRet <- data.frame(cbind(p_val, avg_diff), row.names = mygenes)
  toRet <- toRet[order(toRet$p_val), ]
  return(toRet)
}

diffNegBinom <- function(x, data, vars) {

}

diffAUC <- function(x, y) {
  prediction.use <- prediction(
    predictions = c(x, y),
    labels = c(rep(x = 1, length(x = x)), rep(x = 0, length(x = y))),
    label.ordering = 0:1
  )
  perf.use <- performance(prediction.obj = prediction.use, measure = "auc")
  auc.use <- round(x = perf.use@y.values[[1]], digits = 3)
  return(auc.use)
}

diffLRT <- function(x, y, xmin = 1) {
  lrtX <- bimodLikData(x = x)
  lrtY <- bimodLikData(x = y)
  lrtZ <- bimodLikData(x = c(x, y))
  lrt_diff <- 2 * (lrtX + lrtY - lrtZ)
  return(pchisq(q = lrt_diff, df = 3, lower.tail = F))
}

bimodLikData <- function(x, xmin = 0) {
  x1 <- x[x <= xmin]
  x2 <- x[x > xmin]
  xal <- minmax(
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

makeAlnPlot <- function(proj) {
  alnFile <- paste0("~/big/", proj, "/summary/", proj, ".all.aln.metrics.txt")
  alnData <- read.table(file = alnFile)
  par(mar = c(10, 5, 4, 1))
  mymax <- max(100, max(apply(X = alnData[, 4:7], MARGIN = 1, FUN = sum)))
  x <- barplot(
    height = as.matrix(x = t(x = alnData[, 4:7])),
    names.arg = alnData$V1,
    las = 2,
    col = 1:4,
    ylim = c(0, mymax)
  )
  text(x = x, y = mymax - 10, labels = alnData$V2)
  rp()
}

meanVarPlot <- function(x, ...) {
  myMean <- apply(X = x, MARGIN = 1, FUN = logMeanMinus)
  myVar <- apply(X = x, MARGIN = 1, FUN = logVarMinus)
  plot(x = myMean, y = myVar)
}

#' @export
vsubc <- function(data, code) {
  return(data[grep(pattern = code, x = names(x = data))])
}

#' @export
vminusc <- function(data, code) {
  matchCode <- names(x = data)[grep(pattern = code, x = names(x = data))]
  toIgnore <- which(x = names(x = data) %in% matchCode)
  if (length(x = toIgnore) == 0) {
    return(data)
  }
  return(data[-toIgnore])
}

#' @export
plotVln <- function(
  gene,
  data = dc2,
  code = "rsem",
  mmax = 12,
  getStaty = 1,
  doRet = FALSE,
  doSort = FALSE
) {
  data$GENE <- as.character(x = rownames(x = data))
  a1 <- data[
    gene,
    c(colnames(x = data)[grep(pattern = code, x = colnames(x = data))], "GENE")
  ]
  a2 <- melt(a1, id = "GENE")
  a2$stat <- unlist(x = lapply(
    X = as.character(a2$variable),
    FUN = getStat,
    y = getStaty
  ))
  noise <- rnorm(n = length(x = a2$value)) / 100000
  a2$value <- a2$value + noise
  if (doSort) {
    a2$stat <- factor(
      x = a2$stat,
      levels = names(x = rev(x = sort(x = tapply(
        X = a2$value,
        INDEX = a2$stat,
        FUN = mean
      ))))
    )
  }
  p <- ggplot(data = a2, mapping = aes(x = factor(x = stat), y = value))
  p2 <- p +
    geom_violin(scale = "width", adjust = 0.75, aes(fill = factor(x = stat))) +
    ylab(label = "Expression level (log TPM)")
  p3 <- p2 +
    theme(legend.position = "top") +
    guides(fill = guide_legend(title = NULL)) +
    geom_jitter(height = 0) +
    blackbg +
    xlab(label = "Cell Type")
  p4 <- p3 +
    theme(
      axis.title.x = element_text(face = "bold", colour = "#990000", size = 16),
      axis.text.x  = element_text(angle = 90, vjust = 0.5, size = 12)
    )
  p5 <- p4 +
    theme(
      axis.title.y = element_text(face = "bold", colour = "#990000", size = 16),
      axis.text.y  = element_text(angle = 90, vjust = 0.5, size = 12)
    ) +
    ggtitle(label = gene) +
    theme(plot.title = element_text(size = 20, face = "bold"))
  if (doRet == TRUE) {
    return(p5)
  } else {
    print(p5)
  }
}

#' @export
extract_field <- function(string, field = 1, delim = "_") {
  fields <- as.numeric(x = unlist(x = strsplit(x = as.character(x = field), split = ",")))
  if (length(x = fields) == 1) {
    return(strsplit(x = string, split = delim)[[1]][field])
  }
  return(paste(strsplit(x = string, split = delim)[[1]][fields], collapse = delim))
}

#' @export
genes.ca.range <- function(object, my.min, my.max) {
  ca <- ClusterAlpha(object = object)
  ca.min <- apply(X = ca, MARGIN = 1, FUN = min)
  ca.max <- apply(X = ca, MARGIN = 1, FUN = max)
  genes.1 <- names(x = ca.min[ca.min < my.max])
  genes.2 <- names(x = ca.max[ca.max > my.min])
  return(ainb(a = genes.1, b = genes.2))
}

#' @export
getStat <- function(x, y = 1) {
  return(strsplit(x = x, split = "_")[[1]][y])
}

#' @importFrom grid grid.newpage pushViewport viewport grid.layout
#' @export
MultiPlotList <- function(plots, file, cols = 1, layout = NULL) {
  # Make a list from the ... arguments and plotlist
  #plots <- c(list(...), plotlist)
  numPlots = length(x = plots)
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(x = layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(
      data = seq(from = 1, to = cols * ceiling(x = numPlots / cols)),
      ncol = cols,
      nrow = ceiling(x = numPlots / cols)
    )
  }
  if (numPlots == 1) {
    print(plots[[1]])
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(
      nrow = nrow(x = layout),
      ncol = ncol(x = layout)
    )))
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(x = which(layout == i, arr.ind = TRUE))
      print(
        plots[[i]],
        vp = viewport(
          layout.pos.row = matchidx$row,
          layout.pos.col = matchidx$col
        )
      )
    }
  }
}

#' @export
subc <- function(data, code) {
  return(data[, grep(pattern = code, x = colnames(x = data))])
}

#' @export
minmax <- function(data, min, max) {
  data2 <- data
  data2[data2 > max] <- max
  data2[data2 < min] <- min
  return(data2)
}

#' @export
vp.layout <- function(x, y) {
  viewport(layout.pos.row = x, layout.pos.col = y)
}

#' @export
arrange <- function(..., nrow = NULL, ncol = NULL, as.table = FALSE) {
  dots <- list(...)
  n <- length(x = dots)
  if (is.null(x = nrow) & is.null(x = ncol)) {
    nrow <- floor(x = n / 2)
    ncol <- ceiling(x = n / nrow)
  }
  if (is.null(x = nrow)) {
    nrow <- ceiling(x = n / ncol)
  }
  if (is.null(x = ncol)) {
    ncol <- ceiling(n/nrow)
  }
  ## NOTE see n2mfrow in grDevices for possible alternative
  grid.newpage()
  pushViewport(viewport(layout=grid.layout(nrow,ncol) ) )
  ii.p <- 1
  for (ii.row in seq(from = 1, to = nrow)) {
    if (as.table) {
      ii.table.row <- nrow - ii.table.row + 1
    } else {
      ii.table.row <- ii.row
    }
    for (ii.col in seq(from = 1, to = ncol)) {
      ii.table <- ii.p
      if (ii.p > n) {
        break
      }
      print(dots[[ii.table]], vp = vp.layout(x = ii.table.row, y = ii.col))
      ii.p <- ii.p + 1
    }
  }
}

subSort <- function(vdat, my, fNum, sortBy) {
  vNames <- colnames(x = vdat)
  mycol <- which(x = vNames == sortBy)
  v2 <- vdat[order(vdat[, mycol]), ]
  vsort <- v2
  return(v2)
}

calcMedians <- function(p, start, end) {
  medians <- c()
  for (i in start:end) {
    scores <- p[which(x = p$factorNum == i), i + 6]
    medians <- c(medians, mean(x = scores))
  }
  return(medians)
}

#' @export
myPalette <- function(
  low = "white",
  high = c("green", "red"),
  mid = NULL,
  k = 50
) {
  low <- col2rgb(col = low) / 255
  high <- col2rgb(col = high) / 255
  if (is.null(x = mid)) {
    r <- seq(from = low[1], to = high[1], len = k)
    g <- seq(from = low[2], to = high[2], len = k)
    b <- seq(from = low[3], to = high[3], len = k)
  } else {
    k2 <- round(x = k / 2)
    mid <- col2rgb(col = mid) / 255
    r <- c(
      seq(from = low[1], to = mid[1], len = k2),
      seq(from = mid[1], to = high[1], len = k2)
    )
    g <- c(
      seq(from = low[2], to = mid[2], len = k2),
      seq(from = mid[2], to = high[2],len = k2)
    )
    b <- c(
      seq(from = low[3], to = mid[3], len = k2),
      seq(from = mid[3], to = high[3], len = k2)
    )
  }
  return(rgb(red = r, green = g, blue = b))
}

#' @export
bwCols <- myPalette(low = "white", high="black", k = 50)

comparePCA <- function(a, b) {
  inds <- c(6:26)
  pa <- prcomp(x = a[, inds], scale = TRUE, center = TRUE)
  pb <- prcomp(x = b[, inds], scale = TRUE, center = TRUE)
  print(summary(object = pa))
  print(summary(object = pb))
}

getSmooth <- function(
  vsort,
  myBin,
  n = 0,
  smooth = 0,
  overlap = 0.9,
  type = 0
) {
  if (smooth == 0) {
    delta <- 1 / (1 - overlap)
    nbin <- round(x = (length(x = vsort[, 1]) - delta) / (n - 1 + delta))
    smooth <- (nbin + 1) * delta
  }
  return(smooth)
}

calcTP <- function(cutoff, data, score, real, nTP) {
  return(length(x = which(x = (data[, score] > cutoff) & (data[, real] > 0))) / nTP)
}

calcFP <- function(cutoff, data, score, real, nFP) {
  return(length(x = which(x = (data[, score] > cutoff) & (data[, real] == 0))) / nFP)
}

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

#' @export
pyCols <- myPalette(low = "magenta", high = "yellow", mid = "black")

#' @export
rp <- function() {
  par(mfrow = c(1, 1))
}

calcResidLog <- function(x1, y1, mcut = 30, toAdd = 1) {
  touse <- which(x = (x1 > mcut) & (y1 > mcut))
  x <- log(x = x1 + toAdd, base = 2)
  y <- log(x = y1 + toAdd, base = 2)
  a <- x[touse]
  b <- y[touse]
  myLM <- lm(formula = b ~ a)
  myResid <- y - predict(object = myLM, data.frame(a = x))
  return(myResid)
}

#' @export
logVarDivMean <- function(x) {
  return(log(x = var(x = exp(x = x) - 1) / mean(x = exp(x = x) - 1)))
}

calcResid <- function(x1, y1, mcut = 30, toAdd = 1) {
  touse <- which(x = (x1 > mcut) & (y1 > mcut))
  x <- x1
  y <- y1
  a <- x[touse]
  b <- y[touse]
  myLM <- lm(formula = b ~ a)
  myResid <- y - predict(object = myLM, data.frame(a = x))
  return(myResid)
}

#' @export
findNGene <- function(data, is.expr = 1) {
  toRet <- unlist(x = lapply(
    X = 1:ncol(x = data),
    FUN = function(x) {
      return(gtCut(x = data[, x], cutoff = is.expr))
    }
  ))
  names(x = toRet) <- colnames(x = data)
  return(toRet)
}

#' @export
returnTopX <- function(data, group.by, n.return, col.return = NA) {
  to.ret <- c()
  levels.use=unique(group.by); if (is.factor(group.by)) levels.use=levels(group.by)
  if (!is.na(col.return)) return(unlist(lapply(levels.use, function(x) head(data[group.by==x,col.return],n.return)))) else {
    return(unlist(lapply(levels.use, function(x) head(rownames(data[group.by==x,])))))
  }
}

getCoefs <- function(data, nbin = 20, mycut = 1) {
  my_stats <- data.frame(data[, 1])
  code_humpAvg <- apply(X = data, MARGIN = 1, FUN = humpMean, min = mycut)
  code_humpAvg[code_humpAvg > 9] <- 9
  code_humpAvg[is.na(x = code_humpAvg)] <- 0
  my_stats$code_humpAvg <- code_humpAvg
  data[data > mycut] <- 1
  data[data < mycut] <- 0
  data$bin <- cut(x = code_humpAvg, breaks = nbin)
  data$avg <- code_humpAvg
  rownames(x = my_stats) <- rownames(x = data)
  my_coefs <- data.frame(t(x = sapply(
    X = colnames(x = data[1:(ncol(x = data) - 2)]),
    FUN = getAB,
    data = data,
    data2 = my_stats,
    status = "code",
    code2 = "humpAvg",
    hasBin = TRUE,
    doPlot = FALSE
  )))
  colnames(x = my_coefs) <- c("a", "b")
  return(my_coefs)
}

makeScorePlot2 <- function(
  allscores,
  getStaty = 2,
  mytitle = "Title"
) {
  alls <- data.frame(allscores)
  alls$stat <- unlist(x = lapply(
    X = names(x = allscores),
    FUN = getStat,
    y = getStaty
  ))
  p <- ggplot(data = alls, mapping = aes(x = factor(x = stat), y = allscores))
  p2 <- p +
    geom_violin(
      scale = "width",
      adjust = 0.75,
      mapping = aes(fill = factor(x = stat))
    ) + ylab(label = "Expression level (log TPM)")
  p3 <- p2 +
    theme(legend.position = "top") +
    guides(fill = guide_legend(title = NULL)) +
    geom_jitter(height = 0) +
    blackbg +
    xlab(label = "Cell Type")
  p4 <- p3 +
    theme(
      axis.title.x = element_text(face = "bold", colour = "#990000", size = 16),
      axis.text.x  = element_text(angle = 90, vjust = 0.5, size = 12)
    )
  print(
    p4 +
      theme(
        axis.title.y = element_text(face = "bold", colour = "#990000", size = 16),
        axis.text.y  = element_text(angle = 90, vjust = 0.5, size = 12)
      ) +
      ggtitle(label = mytitle) +
      theme(plot.title = element_text(size = 20, face = "bold"))
  )
}

corCellWeightFast <- function(
  cell1,
  cell2,
  wt_matrix,
  data = subc(data = cell, code = "lps_t2"),
  spear = FALSE
) {
  cell1Drop <- wt_matrix[rownames(x = data), cell1]
  cell2Drop <- wt_matrix[rownames(x = data), cell2]
  my_weights <- cell1Drop * cell2Drop
  cell1Data <- an(data[, cell1])
  cell2Data <- an(data[, cell2])
  # print(my_weights)
  my_weights[is.na(x = my_weights)] <- 0
  if (spear) {
    ret <- corr(
      d = as.matrix(x = cbind(rank(x = cell1Data), rank(x = cell2Data))),
      w = my_weights
    )
  } else {
    ret <- corr(d = as.matrix(x = cbind(cell1Data, cell2Data)), w = my_weights)
  }
  return(ret)
}

covCellWeightFast <- function(
  cell1,
  cell2,
  wt_matrix,
  data = subc(data = cell, data = "lps_t2"),
  spear = FALSE
) {
  cell1Drop <- wt_matrix[rownames(x = data), cell1]
  cell2Drop <- wt_matrix[rownames(x = data), cell2]
  my_weights <- cell1Drop * cell2Drop
  cell1Data <- an(data[, cell1])
  cell2Data <- an(data[, cell2])
  # print(my_weights)
  if (spear) {
    ret <- wtCov(x = rank(x = cell1Data), y = rank(x = cell2Data), w = my_weights)
  } else {
    ret <- wtCov(x = cell1Data, y = cell2Data, w = my_weights)
  }
  return(ret)
}

scaleSCMatrix2 <- function(data, wt_matrix, code = "rsem") {
  wtX <- unlist(x = lapply(
    X = rownames(x = data),
    FUN = function(x) {
      return(sum(data[x, ] * wt_matrix[x, ]) / sum(wt_matrix[x, ]))
    }
  ))
  sdX <- unlist(x = lapply(
    X = rownames(x = data),
    FUN = function(x) {
      return(wtCov(x = data[x, ], y = data[x,], y = wt_matrix[x, ]))
    }
  ))
  sData <- (data - wtX) / sqrt(x = sdX)
  return(sData)
}

wtCov <- function(x, y ,w) {
  w <- w / sum(w)
  wtX <- sum(x * w)
  wtY <- sum(y * w)
  wt_cov <- sum(w * (x - wtX) * (y - wtY))
  return(wt_cov)
}

expAlpha <- function(mu, coefs) {
  logA <- coefs$a
  logB <- coefs$b
  return(exp(x = logA + logB * mu) / (1 + (exp(x = logA + logB * mu))))
}

setWt1 <- function(x, wts, min = 1) {
  wts[x > min] <- 1
  return(wts)
}

#' @export
gtCut <- function(x, cutoff = 1) {
  return(length(x = which(x = x > cutoff)))
}

setWtMatrix1 <- function(data, wt_matrix, mycut = 1) {
  wt1_matrix <- sapply(
    X = 1:ncol(x = data),
    FUN = function(x) {
      return(setWt1(x = data[, x], wts = wt_matrix[, x], min = mycut))
    }
  )
  colnames(x = wt1_matrix) <- colnames(x = data)
  return(wt1_matrix)
}

getAB <- function(
  cn = "lps_t1_S1_rsem",
  code = "lps_t1",
  data = cell,
  data2 = cs,
  code2 = "avg",
  status = "",
  ncut = 25,
  hasBin = FALSE,
  doPlot = FALSE,
  myfunc = plot,
  func2 = lines,
  ...
) {
  if (status == "") {
    status <- getStatus(cn)
  }
  if (! hasBin) {
    data[data > 1] <- 1
    data[data < 1] <- 0
    data$avg <- data2[rownames(x = data), paste0(status, "_", code2)]
    data[data > 9] <- 9
    data$bin <- cut(x = data$avg, breaks = 20)
  }
  data$val <- data[, cn]
  x1 <- (tapply(X = data[, cn], INDEX = data$bin, FUN = mean))
  x2 <- (tapply(X = data[, "avg"], INDEX = data$bin, FUN = mean))
  #glm.out = glm(val~avg,data=data,family=binomial)
  glm.out <- glm(formula = x1 ~ x2, data = data, family = binomial)
  if (doPlot) {
    pred <- (predict(
      object = glm.out,
      data.frame(avg = as.numeric(x = x2)),
      type = "response"
    ))
    myfunc(
      x2,
      x1,
      pch = 16,
      xlab = "Log TPM",
      ylab = "Detection Rate",
      main = cn,
      ylim = c(0, 1)
    )
    func2(x2, pred, ...)
  }
  return(glm.out$coefficients)
}

sensitivityCurve <- function(
  cellName,
  scData,
  bulkData,
  mycut = 1,
  mycex = 1,
  mynew = TRUE,
  ...
) {
  cutLocs <- cut2(x = bulkData, g = 100, onlycuts = TRUE)
  bulkBin <- cut2(x = bulkData, g = 100)
  binaryData <- scData[, cellName]
  binaryData[binaryData > mycut] <- 1
  binaryData[binaryData < mycut] <- 0
  yBin <- tapply(X = binaryData, INDEX = bulkBin, FUN = mean)
  xBin <- tapply(X = bulkData, INDEX = bulkBin, FUN = mean)
  #glm.out = glm(val~avg,data=data,family=binomial)
  options(warn = -1) #otherwise glm throws an unnecessary error
  glm.out = glm(formula = binaryData ~ bulkData, family = binomial)
  options(warn = 0)
  x_vals <- seq(from = 0, to = 10, by = 0.1)
  y_vals <- predict(
    object = glm.out,
    data.frame(bulkData = x_vals),
    type = "response"
  )
  if (mynew) {
    plot(
      x = xBin,
      y = yBin,
      pch = 16,
      xlab = "Average expression",
      ylab = "Probability of detection",
      ...
    )
  }
  lines(x = x_vals, y = y_vals, lwd = 2, ...)
}

#' @export
heatmap2NoKey <- function (
  x,
  Rowv = TRUE,
  Colv = if (symm) "Rowv" else TRUE,
  distfun = dist,
  hclustfun = hclust,
  dendrogram = c("both", "row", "column", "none"),
  symm = FALSE,
  scale = c("none", "row", "column"),
  na.rm = TRUE,
  revC = identical(x = Colv, y = "Rowv"),
  add.expr,
  breaks,
  symbreaks = min(x < 0, na.rm = TRUE) || scale != "none",
  col = "heat.colors",
  colsep,
  rowsep,
  sepcolor = "white",
  sepwidth = c(0.05, 0.05),
  cellnote,
  notecex = 1,
  notecol = "cyan",
  na.color = par("bg"),
  trace = c("column", "row", "both", "none"),
  tracecol = "cyan",
  hline = median(breaks),
  vline = median(breaks),
  linecol = tracecol,
  margins = c(5, 5),
  ColSideColors,
  RowSideColors,
  cexRow = 0.2 + 1 / log10(x = nr),
  cexCol = 0.2 + 1 / log10(x = nc),
  labRow = NULL,
  labCol = NULL,
  key = TRUE,
  keysize = 1.5,
  density.info = c("histogram", "density", "none"),
  denscol = tracecol,
  symkey = min(x < 0, na.rm = TRUE) || symbreaks,
  densadj = 0.25,
  main = NULL,
  xlab = NULL,
  ylab = NULL,
  lmat = NULL,
  lhei = NULL,
  axRowCol="black",
  lwid = NULL,
  dimTitle = NULL,
  ...
) {
  scale01 <- function(x, low = min(x), high = max(x)) {
    return((x - low) / (high - low))
  }
  retval <- list()
  scale <- if (symm && missing(x = scale)) {
    "none"
  } else {
    match.arg(arg = scale)
  }
  dendrogram <- match.arg(arg = dendrogram)
  trace <- match.arg(arg = trace)
  density.info <- match.arg(density.info)
  if (length(x = col) == 1 && is.character(x = col)) {
    col <- get(col, mode = "function")
  }
  if (! missing(x = breaks) && (scale != "none")) {
    warning(
      "Using scale=\"row\" or scale=\"column\" when breaks are",
      "specified can produce unpredictable results.",
      "Please consider using only one or the other."
    )
  }
  if (is.null(x = Rowv) || is.na(x = Rowv)) {
    Rowv <- FALSE
  }
  if (is.null(x = Colv) || is.na(x = Colv)) {
    Colv <- FALSE
  } else if (Colv == "Rowv" && !isTRUE(x = Rowv)) {
    Colv <- FALSE
  }
  if (length(x = di <- dim(x = x)) != 2 || !is.numeric(x = x)) {
    stop("`x' must be a numeric matrix")
  }
  nr <- di[1]
  nc <- di[2]
  if (nr <= 1 || nc <= 1) {
    stop("`x' must have at least 2 rows and 2 columns")
  }
  if (! is.numeric(x = margins) || length(x = margins) != 2) {
    stop("`margins' must be a numeric vector of length 2")
  }
  if (missing(x = cellnote)) {
    cellnote <- matrix(data = "", ncol = ncol(x = x), nrow = nrow(x = x))
  }
  if (! inherits(x = Rowv, what = "dendrogram")) {
    if (((! isTRUE(x = Rowv)) || (is.null(x = Rowv))) &&
        (dendrogram %in% c("both", "row"))) {
      if (is.logical(x = Colv) && (Colv)) {
        dendrogram <- "column"
      } else {
        dedrogram <- "none"
      }
    }
  }
  if (! inherits(x = Colv, what = "dendrogram")) {
    if (((!isTRUE(x = Colv)) || (is.null(x = Colv))) &&
        (dendrogram %in% c("both", "column"))) {
      if (is.logical(x = Rowv) && (Rowv)) {
        dendrogram <- "row"
      } else {
        dendrogram <- "none"
      }
    }
  }
  if (inherits(x = Rowv, what = "dendrogram")) {
    ddr <- Rowv
    rowInd <- order.dendrogram(x = ddr)
  } else if (is.integer(x = Rowv)) {
    hcr <- hclustfun(distfun(x))
    ddr <- as.dendrogram(object = hcr)
    ddr <- reorder(x = ddr, X = Rowv)
    rowInd <- order.dendrogram(x = ddr)
    if (nr != length(x = rowInd)) {
      stop("row dendrogram ordering gave index of wrong length")
    }
  } else if (isTRUE(x = Rowv)) {
    Rowv <- rowMeans(x = x, na.rm = na.rm)
    hcr <- hclustfun(distfun(x))
    ddr <- as.dendrogram(object = hcr)
    ddr <- reorder(x = ddr, X = Rowv)
    rowInd <- order.dendrogram(x = ddr)
    if (nr != length(x = rowInd)) {
      stop("row dendrogram ordering gave index of wrong length")
    }
  } else {
    rowInd <- nr:1
  }
  if (inherits(x = Colv, what = "dendrogram")) {
    ddc <- Colv
    colInd <- order.dendrogram(x = ddc)
  } else if (identical(x = Colv, y = "Rowv")) {
    if (nr != nc) {
      stop("Colv = \"Rowv\" but nrow(x) != ncol(x)")
    }
    if (exists(x = "ddr")) {
      ddc <- ddr
      colInd <- order.dendrogram(x = ddc)
    } else {
      colInd <- rowInd
    }
  } else if (is.integer(x = Colv)) {
    hcc <- hclustfun(distfun(if (symm) {
      x
    } else {
      t(x = x)
    }))
    ddc <- as.dendrogram(object = hcc)
    ddc <- reorder(x = ddc, X = Colv)
    colInd <- order.dendrogram(x = ddc)
    if (nc != length(x = colInd)) {
      stop("column dendrogram ordering gave index of wrong length")
    }
  } else if (isTRUE(x = Colv)) {
    Colv <- colMeans(x = x, na.rm = na.rm)
    hcc <- hclustfun(distfun(if (symm) {
      x
    } else {
      t(x = x)
    }))
    ddc <- as.dendrogram(object = hcc)
    ddc <- reorder(x = ddc, X = Colv)
    colInd <- order.dendrogram(x = ddc)
    if (nc != length(x = colInd)) {
      stop("column dendrogram ordering gave index of wrong length")
    }
  } else {
    colInd <- 1:nc
  }
  retval$rowInd <- rowInd
  retval$colInd <- colInd
  retval$call <- match.call()
  x <- x[rowInd, colInd]
  x.unscaled <- x
  cellnote <- cellnote[rowInd, colInd]
  if (is.null(x = labRow)) {
    labRow <- if (is.null(rownames(x = x))) {
      (1:nr)[rowInd]
    } else {
      rownames(x = x)
    }
  } else {
    labRow <- labRow[rowInd]
  }
  if (is.null(x = labCol)) {
    labCol <- if (is.null(x = colnames(x = x))) {
      (1:nc)[colInd]
    } else {
      colnames(x = x)
    }
  } else {
    labCol <- labCol[colInd]
  }
  if (scale == "row") {
    retval$rowMeans <- rm <- rowMeans(x = x, na.rm = na.rm)
    x <- sweep(x = x, MARGIN = 1, STATS = rm)
    retval$rowSDs <- sx <- apply(X = x, MARGIN = 1, FUN = sd, na.rm = na.rm)
    x <- sweep(x = x, MARGIN = 1, STATS = sx, FUN = "/")
  } else if (scale == "column") {
    retval$colMeans <- rm <- colMeans(x = x, na.rm = na.rm)
    x <- sweep(x = x, MARGIN = 2, STATS = rm)
    retval$colSDs <- sx <- apply(X = x, MARGIN = 2, FUN = sd, na.rm = na.rm)
    x <- sweep(x = x, MARGIN = 2, STATS = sx, FUN = "/")
  }
  if (missing(x = breaks) || is.null(x = breaks) || length(x = breaks) <1) {
    if (missing(x = col) || is.function(x = col)) {
      breaks <- 16
    } else {
      breaks <- length(x = col) + 1
    }
  }
  if (length(x = breaks) == 1) {
    if (! symbreaks) {
      breaks <- seq(
        from = min(x, na.rm = na.rm),
        to = max(x, na.rm = na.rm),
        length = breaks
      )
    } else {
      extreme <- max(abs(x = x), na.rm = TRUE)
      breaks <- seq(from = -extreme, to = extreme, length = breaks)
    }
  }
  nbr <- length(x = breaks)
  ncol <- length(x = breaks) - 1
  if (class(x = col) == "function") {
    col <- col(x = ncol)
  }
  min.breaks <- min(breaks)
  max.breaks <- max(breaks)
  x[x < min.breaks] <- min.breaks
  x[x > max.breaks] <- max.breaks
  #  if (missing(lhei) || is.null(lhei))
  #    lhei <- c(keysize, 4)
  #  if (missing(lwid) || is.null(lwid))
  #    lwid <- c(keysize, 4)
  #  if (missing(lmat) || is.null(lmat)) {
  #    lmat <- rbind(4:3, 2:1)
  #    if (!missing(ColSideColors)) {
  #      if (!is.character(ColSideColors) || length(ColSideColors) !=
  #        nc)
  #        stop("'ColSideColors' must be a character vector of length ncol(x)")
  #      lmat <- rbind(lmat[1, ] + 1, c(NA, 1), lmat[2, ] +
  #        1)
  #      lhei <- c(lhei[1], 0.2, lhei[2])
  #    }
  #    if (!missing(RowSideColors)) {
  #      if (!is.character(RowSideColors) || length(RowSideColors) !=
  #        nr)
  #        stop("'RowSideColors' must be a character vector of length nrow(x)")
  #      lmat <- cbind(lmat[, 1] + 1, c(rep(NA, nrow(lmat) -
  #        1), 1), lmat[, 2] + 1)
  #      lwid <- c(lwid[1], 0.2, lwid[2])
  #    }
  #    lmat[is.na(lmat)] <- 0
  #  }
  #  if (length(lhei) != nrow(lmat))
  #    stop("lhei must have length = nrow(lmat) = ", nrow(lmat))
  #  if (length(lwid) != ncol(lmat))
  #    stop("lwid must have length = ncol(lmat) =", ncol(lmat))
  #  op <- par(no.readonly = TRUE)
  #  on.exit(par(op))
  #  layout(lmat, widths = lwid, heights = lhei, respect = FALSE)

  if (! missing(x = RowSideColors)) {
    par(mar = c(margins[1], 0, 0, 0.5))
    image(x = rbind(1:nr), col = RowSideColors[rowInd], axes = FALSE)
  }
  if (! missing(x = ColSideColors)) {
    par(mar = c(0.5, 0, 0, margins[2]))
    image(x = cbind(1:nc), col = ColSideColors[colInd], axes = FALSE)
  }
  oldMar <- par()$mar
  if (labCol[1] == "") {
    par(mar = c(margins[1]-3, margins[2]-2, margins[1]-3, margins[2]))
  } else {
    par(mar = c(margins[1], margins[2], margins[1], margins[2]))
  }
  x <- t(x = x)
  cellnote <- t(x = cellnote)
  if (revC) {
    iy <- nr:1
    if (exists(x = "ddr")) {
      ddr <- rev(x = ddr)
    }
    x <- x[, iy]
    cellnote <- cellnote[, iy]
  } else {
    iy <- 1:nr
  }
  # add pc number as title if plotting pc heatmaps
  if(is.null(x = dimTitle)) {
    dimTitle <- ""
  }
  #print(dimTitle)
  image(
    x = 1:nc,
    y = 1:nr,
    z = x,
    xlim = 0.5 + c(0, nc),
    ylim = 0.5 + c(0, nr),
    axes = FALSE,
    xlab = "",
    ylab = "",
    main = dimTitle,
    col = col,
    breaks = breaks,
    ...
  )
  retval$carpet <- x
  if (exists(x = "ddr")) {
    retval$rowDendrogram <- ddr
  }
  if (exists(x = "ddc")) {
    retval$colDendrogram <- ddc
  }
  retval$breaks <- breaks
  retval$col <- col
  if (any(is.na(x = x))) {
    mmat <- ifelse(test = is.na(x = x), yes = 1, no = NA)
    image(
      x = 1:nc,
      y = 1:nr,
      z = mmat,
      axes = FALSE,
      xlab = "",
      ylab = "",
      main = pc_title,
      col = na.color,
      add = TRUE
    )
  }
  axis(
    side = 1,
    at = 1:nc,
    labels = labCol,
    las = 2,
    line = -0.5,
    tick = 0,
    cex.axis = cexCol
  )
  if (! is.null(x = xlab)) {
    mtext(text = xlab, side = 1, line = margins[1] - 1.25)
  }
  axis(
    side = 4,
    at = iy,
    labels = labRow,
    las = 2,
    line = -0.5,
    tick = 0,
    cex.axis = cexRow,
    col = axRowCol
  )
  if (! is.null(x = ylab)) {
    mtext(text = ylab, side = 4, line = margins[2] - 1.25)
  }
  if (! missing(x = add.expr)) {
    eval(expr = substitute(expr = add.expr))
  }
  if (! missing(x = colsep)) {
    for (csep in colsep) {
      rect(
        xleft = csep + 0.5,
        ybottom = rep(x = 0, length(x = csep)),
        xright = csep + 0.5 + sepwidth[1],
        ytop = rep(x = ncol(x = x) + 1, csep),
        lty = 1,
        lwd = 1,
        col = sepcolor,
        border = sepcolor
      )
    }
  }
  if (! missing(x = rowsep)) {
    for (rsep in rowsep) {
      rect(
        xleft = 0,
        ybottom = (ncol(x = x) + 1 - rsep) - 0.5,
        xright = nrow(x = x) + 1,
        ytop = (ncol(x) + 1 - rsep) - 0.5 - sepwidth[2],
        lty = 1,
        lwd = 1,
        col = sepcolor,
        border = sepcolor
      )
    }
  }
  min.scale <- min(breaks)
  max.scale <- max(breaks)
  x.scaled <- scale01(x = t(x = x), low = min.scale, high = max.scale)
  if (trace %in% c("both", "column")) {
    retval$vline <- vline
    vline.vals <- scale01(x = vline, low = min.scale, high = max.scale)
    for (i in colInd) {
      if (! is.null(x = vline)) {
        abline(
          v = i - 0.5 + vline.vals,
          col = linecol,
          lty = 2
        )
      }
      xv <- rep(x = i, nrow(x = x.scaled)) + x.scaled[, i] - 0.5
      xv <- c(xv[1], xv)
      yv <- 1:length(x = xv) - 0.5
      ##lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
    }
  }
  if (trace %in% c("both", "row")) {
    retval$hline <- hline
    hline.vals <- scale01(x = hline, low = min.scale, high = max.scale)
    for (i in rowInd) {
      if (! is.null(x = hline)) {
        abline(h = i + hline, col = linecol, lty = 2)
      }
      yv <- rep(x = i, ncol(x = x.scaled)) + x.scaled[i, ] - 0.5
      yv <- rev(x = c(yv[1], yv))
      xv <- length(x = yv):1 - 0.5
      ##lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
    }
  }
  if (! missing(x = cellnote)) {
    text(
      x = c(row(x = cellnote)),
      y = c(col(x = cellnote)),
      labels = c(cellnote),
      col = notecol,
      cex = notecex
    )
  }
  #par(mar = c(margins[1], 0, 0, 0))
  if (dendrogram %in% c("both", "row")) {
    ##plot(ddr, horiz = TRUE, axes = FALSE, yaxs = "i", leaflab = "none")
  }
  ##else plot.new()
  #par(mar = c(0, 0, if (!is.null(main)) 5 else 0, margins[2]))
  if (dendrogram %in% c("both", "column")) {
    ##plot(ddc, axes = FALSE, xaxs = "i", leaflab = "none")
  }
  ##else plot.new()
  key <- FALSE
  if (! is.null(x = main))
    title(main = main, cex.main = 1.5 * op[["cex.main"]])
  if (key) {
    par(mar = c(5, 4, 2, 1), cex = 0.75)
    tmpbreaks <- breaks
    if (symkey) {
      max.raw <- max(abs(x = c(x, breaks)), na.rm = TRUE)
      min.raw <- -max.raw
      tmpbreaks[1] <- -max(abs(x = x), na.rm = TRUE)
      tmpbreaks[length(x = tmpbreaks)] <- max(abs(x = x), na.rm = TRUE)
    } else {
      min.raw <- min(x, na.rm = TRUE)
      max.raw <- max(x, na.rm = TRUE)
    }
    z <- seq(from = min.raw, to = max.raw, length = length(x = col))
    #image(z = matrix(z, ncol = 1), col = col, breaks = tmpbreaks,
    #      xaxt = "n", yaxt = "n")
    par(usr = c(0, 1, 0, 1))
    lv <- pretty(x = breaks)
    xv <- scale01(x = as.numeric(x = lv), low = min.raw, high = max.raw)
    axis(side = 1, at = xv, labels = lv)
    if (scale == "row") {
      mtext(side = 1, "Row Z-Score", line = 2)
    }
    else if (scale == "column") {
      mtext(side = 1, "Column Z-Score", line = 2)
    }
    else {
      mtext(side = 1, "Value", line = 2)
    }
    if (density.info == "density") {
      dens <- density(x = x, adjust = densadj, na.rm = TRUE)
      omit <- dens$x < min(breaks) | dens$x > max(breaks)
      dens$x <- dens$x[-omit]
      dens$y <- dens$y[-omit]
      dens$x <- scale01(x = dens$x, low = min.raw, high = max.raw)
      lines(
        x = dens$x,
        y = dens$y / max(dens$y) * 0.95,
        col = denscol,
        lwd = 1
      )
      axis(
        side = 2,
        at = pretty(x = dens$y) / max(dens$y) * 0.95,
        pretty(dens$y)
      )
      title(main = "Color Key\nand Density Plot")
      par(cex = 0.5)
      mtext(side = 2, "Density", line = 2)
    } else if (density.info == "histogram") {
      h <- hist(x = x, plot = FALSE, breaks = breaks)
      hx <- scale01(x = breaks, low = min.raw, high = max.raw)
      hy <- c(h$counts, h$counts[length(h$counts)])
      lines(
        x = hx,
        y = hy / max(hy) * 0.95,
        lwd = 1,
        type = "s",
        col = denscol
      )
      axis(
        side = 2,
        at = pretty(x = hy) / max(hy) * 0.95,
        pretty(x = hy)
      )
      title(main = "Color Key\nand Histogram")
      par(cex = 0.5)
      mtext(side = 2, "Count", line = 2)
    } else {
      title(main = "Color Key")
    }
  }
  ##else plot.new()
  retval$colorTable <- data.frame(
    low = retval$breaks[-length(x = retval$breaks)],
    high = retval$breaks[-1], color = retval$col
  )
  invisible(x = retval)
  par(mar = oldMar)
}

mergeDescendents <- function(object, tree, node, pcs, top.genes, acc.cutoff) {
  # find daughter cells of given node in given tree
  daughters <- tree$edge[which(x = tree$edge[, 1] == node), 2]
  # get the children of both daughters
  childNodes <- 1:(tree$Nnode + 1)
  if (length(x = ainb(a = c(daughters[1], daughters[2]), b = childNodes)) == 2) {
    d1 = WhichCells(object = object, ident = daughters[1])
    d2 = WhichCells(object = object, ident = daughters[2])
    y  = as.numeric(x = object@ident[c(d1, d2)]) - 1
    x  = data.frame(t(x = object@data[
      PCTopGenes(object = object, pc.use = pcs,num.genes = top.genes),
      c(d1, d2)
    ]))
    xv = apply(X = x, MARGIN = 2, FUN = var)
    x  = x[, names(x = xv > 0)]
    # run k-fold cross validation
    ctrl = trainControl(method = "repeatedcv", repeats = 5)
    set.seed(seed = 1500)
    model = train(
      formula = as.factor(x = y) ~ .,
      data = x,
      method = "svmLinear",
      trControl = ctrl
    )
    acc = model$results[, 2]
    # if classifier can't classify them well enough, merge clusters
    if (acc < acc.cutoff) {
      object <- SetIdent(
        object = object,
        cells.use = WhichCells(object = object, ident =daughters[1]),
        ident.use = daughters[2]
      )
    }
    return(object)
  }
  # recursion to traverse tree
  if(daughters[1] %in% childNodes){
    object = mergeDescendents(
      object = object,
      tree = tree,
      node = daughters[2],
      pcs = pcs,
      top.genes = top.genes,
      acc.cutoff = acc.cutoff
    )
    return(object)
  }
  if (daughters[2] %in% childNodes) {
    object = mergeDescendents(
      object = object,
      tree = tree,
      node = daughters[1],
      pcs = pcs,
      top.genes = top.genes,
      acc.cutoff = acc.cutoff
    )
    return(object)
  }
  object <- mergeDescendents(
    object = object,
    tree = tree,
    node = daughters[1],
    pcs = pcs,
    top.genes = top.genes,
    acc.cutoff = acc.cutoff
  )
  object <- mergeDescendents(
    object = object,
    tree = tree,
    node = daughters[2],
    pcs = pcs,
    top.genes = top.genes,
    acc.cutoff = acc.cutoff
  )
  return(object)
}

set.ifnull <- function(x, y) {
  if (is.null(x = x)) {
    return(y)
  }
  return(x)
}

expAlpha <- function(mu, coefs) {
  logA <- coefs$a
  logB <- coefs$b
  return(exp(x = logA + logB * mu) / (1 + (exp(x = logA + logB * mu))))
}

regression.sig <- function(x, score, data, latent, code = "rsem") {
  if (var(x = as.numeric(x = subc(data = data, code = code)[x, ])) == 0) {
    return(0)
  }
  latent <- latent[grep(pattern = code, x = names(x = data))]
  data <- rbind(subc(data = data, code = code), vsubc(data = score, code = code))
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

same <- function(x) {
  return(x)
}

nb.residuals <- function(fmla, regression.mat, gene) {
  fit <- 0
  try(expr = fit <- glm.nb(formula = fmla, data = regression.mat), silent=TRUE)
  if (class(fit)[1] == 'numeric') {
    message(sprintf(
      'glm.nb failed for gene %s; trying again with glm and family=negative.binomial(theta=0.1)',
      gene
    ))
    try(
      expr = fit <- glm(
        formula = fmla,
        data = regression.mat,
        family = negative.binomial(theta = 0.1)
      ),
      silent=TRUE
    )
    if (class(fit)[1] == 'numeric') {
      message('glm and family=negative.binomial(theta=0.1) failed; falling back to scale(log10(y+1))')
      return(scale(x = log10(x = regression.mat[, 'GENE'] + 1))[, 1])
    }
  }
  return(residuals(object = fit, type='pearson'))
}

# given a UMI count matrix, estimate NB theta parameter for each gene
# and use fit of relationship with mean to assign regularized theta to each gene
theta.reg <- function(cm, latent.data, min.theta = 0.01, bin.size = 128) {
  genes.regress <- rownames(x = cm)
  bin.ind <- ceiling(x = 1:length(x = genes.regress) / bin.size)
  max.bin <- max(bin.ind)
  print('Running Poisson regression (to get initial mean), and theta estimation per gene')
  pb <- txtProgressBar(min = 0, max = max.bin, style = 3)
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
      cat(sprintf(
        'Used loess with span %1.2f to fit mean-variance relationship\n',
        span
      ))
      break
    }
  }
  if (any(is.na(x = fit$fitted))) {
    stop('Problem when fitting NB gene variance in theta.reg - NA values were fitted.')
  }
  theta.fit <- (UMI.mean ^ 2) / ((10 ^ fit$fitted) - UMI.mean)
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
  return(theta.fit)
}

# compare two negative binomial regression models
# model one uses only common factors (com.fac)
# model two additionally uses group factor (grp.fac)
de.nb.reg <- function(y, theta, latent.data, com.fac, grp.fac) {
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
  log.fc <- log2(x = exp(x = coef(object = fit4)[foi])) #log.fc <- log2(1/exp(coef(fit4)[foi]))
  ret <- c(
    fit2$deviance,
    fit4$deviance,
    pval,
    coef(object = fit4)[foi],
    log.fc,
    freqs
  )
  names(x = ret) <- c(
    'dev1',
    'dev2',
    'pval',
    'coef',
    'log.fc',
    'freq1',
    'freq2'
  )
  return(ret)
}
