# Downsamples a matrix, returns vector of row (cell) names
#
# @param matrix An embedding with cells as rows
# @param size A non-negative integer giving the number of rows (cells) to choose
# @param method Sampling strategy, one of 'random', 'max', 'prod', 'mean'
# @param eps A small value added to the density estimate
# @param bw.adjust Factor applied to default bandwidth

DownsampleMatrix <- function(
  mat,
  size,
  method = 'max',
  eps = 0.001,
  bw.adjust = 2,
  ...
) {
  if (!(method %in% c('random', 'max', 'prod', 'mean'))) {
    stop('method in DownsampleMatrix needs to be one of "random", "max", "prod", "mean"')
  }
  if (method == 'random') {
    return(sample(x = rownames(x = mat), size = size))
  }
  # get the density of the individual dimensions per cell
  dens <- apply(
    X = mat,
    MARGIN = 2,
    FUN = function(y) {
      y.dens <- density(x = y, bw = 'nrd', adjust = bw.adjust)
      return(approx(x = y.dens$x, y = y.dens$y, xout = y)$y)
    }
  )
  # probability of sampling a cell should be inversely related
  # to the density around that cell
  weights <- 1 / (dens + eps)
  # for a given dimension all weights should add up to one
  weights <- sweep(
    x = weights,
    MARGIN = 2,
    STATS = apply(X = weights, MARGIN = 2, FUN = sum),
    FUN = '/'
  )
  # combine the weights across dimensions using the function specified in method
  sampling.prob <- apply(X = weights, MARGIN = 1, FUN = method)
  cells <- sample(x = rownames(x = mat), size = size, prob = sampling.prob)
  return(cells)
}
