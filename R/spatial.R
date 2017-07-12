
#' Get cell centroids
#'
#' Calculate the spatial mapping centroids for each cell, based on previously
#' calculated mapping probabilities for each bin.
#'
#' Currently, Seurat assumes that the tissue of interest has an 8x8 bin
#' structure. This will be broadened in a future release.
#'
#' @param object Seurat object
#' @param cells.use Cells to calculate centroids for (default is all cells)
#' @param get.exact Get exact centroid (Default is TRUE). If FALSE, identify
#' the single closest bin.
#'
#' @return Data frame containing the x and y coordinates for each cell
#' centroid.
#'
#' @export
#'
GetCentroids <- function(object, cells.use = NULL, get.exact = TRUE) {
  cells.use <- SetIfNull(x = cells.use, default = colnames(x = object@final.prob))
  #Error checking
  cell.names <- ainb(a = cells.use, b = colnames(x = object@final.prob))
  if (length(x = cell.names) != length(x = cells.use)) {
    print(paste(
      "Error",
      anotinb(x = cells.use, y = colnames(x = object@final.prob)),
      " have not been mapped"
    ))
    return(0)
  }
  if (get.exact) {
    my.centroids <- data.frame(t(x = sapply(
      X = colnames(x = object@data),
      FUN = function(x) {
        return(exact.cell.centroid(cell.probs = object@final.prob[, x]))
      }
    )))
  } else {
    my.centroids <- data.frame(t(x = sapply(
      X = colnames(x = object@data),
      FUN = function(x) {
        return(cell.centroid(cell.probs = object@final.prob[, x]))
      }
    )))
  }
  colnames(x = my.centroids) <- c("bin.x","bin.y")
  return(my.centroids)
}

#' Quantitative refinement of spatial inferences
#'
#' Refines the initial mapping with more complex models that allow gene
#' expression to vary quantitatively across bins (instead of 'on' or 'off'),
#' and that also considers the covariance structure between genes.
#'
#' Full details given in spatial mapping manuscript.
#'
#' @param object Seurat object
#' @param genes.use Genes to use to drive the refinement procedure.
#'
#' @return Seurat object, where mapping probabilities for each bin are stored
#' in object@@final.prob
#'
#' @import fpc
#'
#' @export
#'
RefinedMapping <- function(object, genes.use) {
  genes.use <- ainb(a = genes.use, b = rownames(x = object@imputed))
  cells.max <- t(x = sapply(
    X = colnames(object@data),
    FUN = function(x) {
      return(exact.cell.centroid(object@final.prob[, x]))
    }
  ))
  all.mu <- sapply(
    X = genes.use,
    FUN = function(gene) {
      return(sapply(X = 1:64, FUN = function(bin) {
        mean(x = as.numeric(x = object@imputed[
          gene, # Row
          fetch.closest(
            bin = bin,
            all.centroids = cells.max,
            num.cell = 2*length(x = genes.use)
          ) # Column
          ]))
      }))
    }
  )
  all.cov <- list()
  for (x in 1:64) {
    all.cov[[x]] <- cov(
      x = t(
        x = object@imputed[
          genes.use, # Row
          fetch.closest(
            bin = x,
            all.centroids = cells.max,
            num.cell = 2*length(x = genes.use)
          ) # Columns
          ]
      )
    )
  }
  mv.probs <- sapply(
    X = colnames(x = object@data),
    FUN = function(my.cell) {
      return(sapply(X = 1:64, FUN = function(bin) {
        return(slimdmvnorm(
          x = as.numeric(x = object@imputed[genes.use, my.cell]),
          mean = as.numeric(x = all.mu[bin, genes.use]),
          sigma = all.cov[[bin]])
        )
      }))
    }
  )
  mv.final <- exp(
    x = sweep(
      x = mv.probs,
      MARGIN = 2,
      STATS = apply(X = mv.probs, MARGIN = 2, FUN = log_add)
    )
  )
  object@final.prob <- data.frame(mv.final)
  return(object)
}

#' Infer spatial origins for single cells
#'
#' Probabilistically maps single cells based on (imputed) gene expression
#' estimates, a set of mixture models, and an in situ spatial reference map.
#'
#' @param object Seurat object
#' @param cells.use Which cells to map
#'
#' @return Seurat object, where mapping probabilities for each bin are stored
#' in object@@final.prob
#'
#' @export
#'
InitialMapping <- function(object, cells.use = NULL) {
  cells.use <- SetIfNull(x = cells.use, default = colnames(x = object@data))
  every.prob <- sapply(
    X = cells.use,
    FUN = function(x) {
      return(MapCell(
        object = object,
        cell.name = x,
        do.plot = FALSE,
        safe.use = FALSE
      ))
    }
  )
  object@final.prob <- data.frame(every.prob)
  rownames(x = object@final.prob) <- paste0("bin.", rownames(x = object@final.prob))
  return(object)
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
