# Downsamples a matrix, returns vector of row (cell) names
#
# @param mat An embedding with cells as rows
# @param size A non-negative integer giving the number of rows (cells) to choose
# @param method Sampling strategy, one of 'random', 'max', 'prod', 'mean'
# @param eps A small value added to the density estimate
# @param bw.adjust Factor applied to default bandwidth

DownsampleMatrix <- function(
  mat,
  size,
  method = 'max',
  eps = .Machine$double.eps,
  bw.adjust = 0.5,
  do.pca = FALSE,
  norm.dens = FALSE,
  q = 0.99,
  prob.exp = 1,
  add.attr = FALSE
) {
  if (method == 'random') {
    return(sample(x = rownames(x = mat), size = size))
  }
  if (size == nrow(x = mat)) {
    return(rownames(x = mat))
  }
  if (do.pca) {
    mat <- prcomp(mat, center = TRUE, scale. = TRUE)$x
  }
  # get the density of the individual dimensions per cell
  dens <- apply(
    X = mat,
    MARGIN = 2,
    FUN = function(y) {
      y.dens <- density(x = y, bw = 'nrd', adjust = bw.adjust)
      return(approx(x = y.dens$x, y = y.dens$y, xout = y)$y + eps)
    }
  )
  if (norm.dens) {
    dens <- sweep(
      x = dens,
      MARGIN = 2,
      STATS = apply(X = dens, MARGIN = 2, FUN = sum),
      FUN = '/'
    )
  }
  # probability of sampling a cell should be inversely related
  # to the density around that cell
  weights <- 1 / dens
  # 
  #weights <- apply(X = weights, MARGIN = 2, FUN = function(x) x * rank(x))
  # for a given dimension all weights should add up to one
  # weights <- sweep(
  #  x = weights,
  #  MARGIN = 2,
  #  STATS = apply(X = weights, MARGIN = 2, FUN = sum),
  #  FUN = '/'
  # )
  #norm.factor <- apply(dens, 2, sum)
  norm.factor <- apply(weights, 2, quantile, q)
  weights <- sweep(
    x = weights,
    MARGIN = 2,
    STATS = norm.factor,
    FUN = '/'
  )
  if (method == 'dev') {
    dens.norm <- dens
    for (pc in 1:ncol(dens)) {
      o <- order(mat[, pc])
      n <- 5001
      dens.rm <- roll_mean(dens[o, pc], n=n, align='center')
      dens.rm <- c(rep(dens.rm[1], floor(n/2)), dens.rm, rep(dens.rm[length(dens.rm)], floor(n/2)))
      dens.norm[o, pc] <- dens.norm[o, pc] / dens.rm
    }
    #dens.norm <- sweep(x = dens, MARGIN = 2, STATS = apply(X = dens, MARGIN = 2, FUN = sum), FUN = '/')
    #dens.norm <- sweep(x = dens.norm, MARGIN = 2, STATS = apply(dens.norm, 2, function(x) sum(sort(x)[1:1000])), FUN = '*')
    #dens.norm <- dens / apply(dens, 1, sum)
    #dens.norm <- apply(dens.norm, 2, function(x) cumsum(sort(x))[rownames(dens.norm)])
    p.mat <- 1/(dens.norm)
    #p.mat <- sweep(p.mat, 2, apply(p.mat, 2, sum), FUN = '/')
    #p.mat <- sweep(p.mat, 2, apply(p.mat, 2, function(x) -sum(sort(-x)[500])), FUN = '/')
    p.mat <- sweep(p.mat, 2, apply(p.mat, 2, quantile, 0.99), FUN = '/')
    #p.mat <- sweep(p.mat, 2, apply(p.mat, 2, quantile, 0.75) - apply(p.mat, 2, quantile, 0.25), FUN = '/')
    #p.mat <- sweep(p.mat, 2, apply(p.mat, 2, sd), FUN = '/')
    #p.mat <- p.mat / apply(p.mat, 1, sum)
    p <- apply(p.mat, 1, max)
    #dens.min <- apply(dens.norm, 1, min)
    #p <- 1/(dens.min^prob.exp)
    o <- sample(1:nrow(dens), size, prob = p^prob.exp)
    return(rownames(mat)[o])
    th <- 0.05
    idx.union <- c()
    for (pc in 1:ncol(mat)) {
      o <- sample(1:nrow(dens), size, prob = 1/(dens[, pc]^prob.exp))
      cum.dens <- cumsum(dens[o, pc])
      cum.dens.norm <- cum.dens/sum(dens[, pc])
      cum.weight <- cumsum(1/dens[o, pc])
      cum.weight.norm <- cum.weight/sum(1/dens[, pc])
      #idx.union <- union(idx.union, o[cum.dens.norm <= cdn.th])
      idx.union <- union(idx.union, o[cum.weight.norm <= th])
      message(pc, ' ', sum(cum.weight.norm <= th))
    }
    return(rownames(mat)[idx.union])
  }
  if (method == 'nnd') {
    k <- floor(nrow(mat) * 0.01)
    nn.res <- RANN::nn2(mat, k = k + 1, searchtype = 'standard', eps = 1/2)
    sampling.prob <- nn.res$nn.dists[, k + 1]^prob.exp
    cells <- sample(x = rownames(mat), size = size, prob = sampling.prob)
  }
  if (method == 'indep') {
    # determine how to distribute the cells
    #size.factor <- 1 / apply(dens, 2, sd)
    #size.factor <- sqrt(apply(dens, 2, sum))
    size.factor <- rep(1, ncol(dens))
    size.pc <- ceiling((size.factor * size) / sum(size.factor))
    #print(size.pc)
    sampled <- rep(FALSE, nrow(dens))
    for (i in 1:ncol(dens)) {
      p <- 1 / dens[, i]
      p[sampled] <- 0
      sampled[sample(1:nrow(dens), size = size.pc[i], prob = p^prob.exp)] <- TRUE
    }
    cells <- rownames(x = mat)[sampled]
  }
  if (method %in% c('max', 'mean', 'gm.mean', 'prod', 'min')) {
    # combine the weights across dimensions using the function specified in method
    sampling.prob <- apply(X = weights, MARGIN = 1, FUN = method)^prob.exp
    cells <- sample(x = rownames(x = mat), size = size, prob = sampling.prob)
  }
  if (add.attr) {
    attr(cells, 'dens') <- dens
    attr(cells, 'weights') <- weights
    attr(cells, 'sampling.prob') <- sampling.prob
  }
  return(cells)
}

