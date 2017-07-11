
#' Draw 3D in situ predictions from Zebrafish dataset
#'
#' From Jeff Farrell
#'
#' @param data Predicted expression levels across Zebrafish bins
#' @param label Plot label
#'
#' @export
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

#if we don't use this in the package, we should delete it
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

#return average of all values greater than a threshold
humpMean <- function(x, min = 0) {
  return(mean(x = x[x > min]))
}

#return variance of all values greater than a threshold
humpVar <- function(x, min = 0) {
  return(var(x = x[x > min]))
}


#calculate refined mapping probabilites based on multivariate distribution
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

#' Return a subset of rows for a matrix or data frame
#'
#' @param data Matrix or data frame with row names
#' @param code Pattern for matching within row names 
#' @return Returns a subset of data, using only rownames that yielded a match to the pattern
#' @export
subr <- function(data, code) {
  return(data[grep(pattern = code, x = rownames(x = data)), ])
}

#calculate the coefficient of variation 
cv <- function(x) {
  return(sd(x = x) / mean(x = x))
}

#return la count of all values greater than a threshold
humpCt <- function(x, min = 0) {
  return(length(x = x[x > min]))
}

#add values in log-space
log_add <- function(x) {
  mpi <- max(x)
  return(mpi + log(x = sum(exp(x = x - mpi))))
}

#' Return a subset of rows for a matrix or data frame
#'
#' @param data Matrix or data frame with row names
#' @param code Pattern for matching within row names 
#' @return Returns a subset of data, using only rownames that did not yield a match to the pattern
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

#' Return a subset of columns for a matrix or data frame
#'
#' @param data Matrix or data frame with column names
#' @param code Pattern for matching within column names 
#' @return Returns a subset of data, using only column names that did not yield a match to the pattern
#' @export
minusc <- function(data, code) {
  matchCode <- colnames(x = data)[grep(pattern = code, colnames(x = data))]
  toIgnore <- which(x = colnames(x = data) %in% matchCode)
  if (length(x = toIgnore) == 0) {
    return(data)
  }
  return(data[, -toIgnore])
}

#returns the intersection of two vectors, i.e all elements of a that are also in b
#we should be switching to use the intersect function instead
ainb <- function(a, b) {
  return(a[a %in% b])
}


#calculates binned scores for finding highly variable genes
meanNormFunction <- function(data, myfuncX, myfuncY, nBin = 20) {
  data_x <- apply(X = data, MARGIN = 1, FUN = myfuncX)
  data_y <- apply(X = data, MARGIN = 1, FUN = myfuncY)
  data_x_bin <- cut(x = data_x, breaks = nBin)
  names(x = data_x_bin) <- names(x = data_x)
  mean_y <- tapply(X = data_y, INDEX = data_x_bin, FUN = mean)
  sd_y <- tapply(X = data_y, INDEX = data_x_bin, FUN = sd)
  return((data_y-mean_y[as.numeric(x = data_x_bin)]) / sd_y[as.numeric(x = data_x_bin)])
}

#internal function for spatial mapping
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

#remove legend title
no.legend.title <- theme(legend.title = element_blank())

#set legend text
gg.legend.text <- function(x = 12, y = "bold") {
  return(theme(legend.text = element_text(size = x, face = y)))
}

#set legend point size
gg.legend.pts <- function(x = 6) {
  return(guides(colour = guide_legend(override.aes = list(size = x))))
}

#set x axis features
gg.xax <- function(x = 16, y = "#990000", z = "bold", x2 = 12) {
  return(theme(
    axis.title.x = element_text(face = z, colour = y, size = x),
    axis.text.x = element_text(angle = 90, vjust = 0.5, size = x2)
  ))
}

#set y axis features
gg.yax <- function(x = 16, y = "#990000", z = "bold", x2 = 12) {
  return(theme(
    axis.title.y = element_text(face = z, colour = y, size = x),
    axis.text.y = element_text(angle = 90, vjust = 0.5, size = x2)
  ))
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


#' Independently shuffle values within each row of a matrix
#' 
#' Creates a matrix where correlation structure has been removed, but overall values are the same
#'
#' @param x Matrix to shuffle
#' @return Returns a scrambled matrix, where each row is shuffled independently
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


#Calculate variance of logged values in non-log space (return answer in log-space)
#' @export
expVar <- function(x) {
  return(log1p(var(expm1(x))))
}

#Calculate SD of logged values in non-log space (return answer in log-space)
#' @export
expSD <- function(x) {
  return(log1p(sd(expm1(x))))
}

#Calculate mean of logged values in non-log space (return answer in log-space)
expMean <- function(x) {
  return(log(x = mean(x = exp(x = x) - 1) + 1))
}


#fetch closest bin, used internally in spatial mapping
fetch.closest <- function(bin, all.centroids, num.cell) {
  bin.y <- (bin - 1) %/% 8 + 1
  bin.x <- (bin - 1) %% 8 + 1
  all.centroids <- rbind(all.centroids, c(bin.x, bin.y))
  all.dist <- as.matrix(x = dist(x = all.centroids))
  return(names(x = sort(x = all.dist[nrow(x = all.dist), ]))[2:(num.cell + 2)])
}


#return cell centroid after spatial mappings (both X and Y)
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

#return x-coordinate cell centroid 
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

#return y-coordinate cell centroid 
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

#return x and y-coordinate cell centroid 
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

#returns the setdiff of two vectors, i.e all elements of a that are not in b
#we should be switching to use the setdiff function instead
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

#' Return a subset of columns for a matrix or data frame
#'
#' @param data Matrix or data frame with column names
#' @param code Pattern for matching within column names 
#' @return Returns a subset of data, using only column names that yield a match to the pattern
#' @export
subc <- function(data, code) {
  return(data[, grep(pattern = code, x = colnames(x = data))])
}

#' Apply a ceiling and floor to all values in a matrix
#'
#' @param data Matrix or data frame 
#' @param min all values below this min value will be replaced with min
#' @param max all values above this max value will be replaced with max
#' @return Returns matrix after performing these floor and ceil operations
#' @export
minmax <- function(data, min, max) {
  data2 <- data
  data2[data2 > max] <- max
  data2[data2 < min] <- min
  return(data2)
}

#used for zebrafish plotting
vp.layout <- function(x, y) {
  viewport(layout.pos.row = x, layout.pos.col = y)
}




#' Create a custom color palette
#'
#' Creates a custom color palette based on low, middle, and high color values
#'
#' @param low low color 
#' @param high high color
#' @param mid middle color. Optional.
#' @param k number of steps (colors levels) to include between low and high values 
#' @export
myPalette <- function(
  low = "white",
  high = "red",
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

#shortcut to make black-white palette
#' @export
bwCols <- myPalette(low = "white", high="black", k = 50)

#' @export
#shortcut to make purple-yellow palette, which is default in most Seurat heatmaps
pyCols <- myPalette(low = "magenta", high = "yellow", mid = "black")


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


#heatmap.2, but does not draw a key.
#unclear if this is necessary, but valuable to have the function coded in for modifications
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
