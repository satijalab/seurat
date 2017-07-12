#internal function for spatial mapping
shift.cell <- function(bin, x, y) {
  bin.y <- (bin - 1) %/% 8 + 1
  bin.x <- (bin - 1) %% 8 + 1
  new.x <- minmax(data = bin.x + x, min = 1, max = 8)
  new.y <- minmax(data = bin.y + y, min = 1, max = 8)
  new.bin <- 8 * (new.y - 1) + new.x
  return(new.bin)
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

#fetch closest bin, used internally in spatial mapping
fetch.closest <- function(bin, all.centroids, num.cell) {
  bin.y <- (bin - 1) %/% 8 + 1
  bin.x <- (bin - 1) %% 8 + 1
  all.centroids <- rbind(all.centroids, c(bin.x, bin.y))
  all.dist <- as.matrix(x = dist(x = all.centroids))
  return(names(x = sort(x = all.dist[nrow(x = all.dist), ]))[2:(num.cell + 2)])
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
