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
