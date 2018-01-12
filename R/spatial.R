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
#' @examples
#' \dontrun{
#' # Note that the PBMC test example object does not contain spatially restricted
#' # examples below are only demonstrate code
#' pmbc_small <- GetCentroids(pbmc_small, cells.use=pbmc_small@cell.names)
#' }
#'
GetCentroids <- function(object, cells.use = NULL, get.exact = TRUE) {
  cells.use <- SetIfNull(x = cells.use, default = colnames(x = object@spatial@finalprob))
  #Error checking
  cell.names <- intersect(x = cells.use, y = colnames(x = object@spatial@finalprob))
  if (length(x = cell.names) != length(x = cells.use)) {
    print(paste(
      "Error",
      setdiff(x = cells.use, y = colnames(x = object@spatial@finalprob)),
      " have not been mapped"
    ))
    return(0)
  }
  if (get.exact) {
    my.centroids <- data.frame(t(x = sapply(
      X = colnames(x = object@data),
      FUN = function(x) {
        return(ExactCellCentroid(cell.probs = object@spatial@finalprob[, x]))
      }
    )))
  } else {
    my.centroids <- data.frame(t(x = sapply(
      X = colnames(x = object@data),
      FUN = function(x) {
        return(CellCentroid(cell.probs = object@spatial@finalprob[, x]))
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
#' @importFrom stats cov
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Note that the PBMC test example object does not contain spatially restricted
#' # examples below are only demonstrate code
#' pmbc_small <- RefinedMapping(pbmc_small, genes.use=pbmc_small@var.genes)
#' }
#'
RefinedMapping <- function(object, genes.use) {
  genes.use <- intersect(x = genes.use, y = rownames(x = object@imputed))
  cells.max <- t(x = sapply(
    X = colnames(object@data),
    FUN = function(x) {
      return(ExactCellCentroid(object@spatial@finalprob[, x]))
    }
  ))
  all.mu <- sapply(
    X = genes.use,
    FUN = function(gene) {
      return(sapply(X = 1:64, FUN = function(bin) {
        mean(x = as.numeric(x = object@imputed[
          gene, # Row
          FetchClosest(
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
          FetchClosest(
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
      STATS = apply(X = mv.probs, MARGIN = 2, FUN = LogAdd)
    )
  )
  object@spatial@finalprob <- data.frame(mv.final)
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
#' @examples
#' \dontrun{
#' # Note that the PBMC test example object does not contain spatially restricted
#' # examples below are only demonstrate code
#' pmbc_small <- InitialMapping(pbmc_small)
#' }
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
  object@spatial@finalprob <- data.frame(every.prob)
  rownames(x = object@spatial@finalprob) <- paste0("bin.", rownames(x = object@spatial@finalprob))
  return(object)
}

#' Build mixture models of gene expression
#'
#' Models the imputed gene expression values as a mixture of gaussian
#' distributions. For a two-state model, estimates the probability that a given
#' cell is in the 'on' or 'off' state for any gene. Followed by a greedy
#' k-means step where cells are allowed to flip states based on the overall
#' structure of the data (see Manuscript for details)
#'
#' @param object Seurat object
#' @param gene Gene to fit
#' @param do.k Number of modes for the mixture model (default is 2)
#' @param num.iter Number of 'greedy k-means' iterations (default is 1)
#' @param do.plot Plot mixture model results
#' @param genes.use Genes to use in the greedy k-means step (See manuscript for details)
#' @param start.pct Initial estimates of the percentage of cells in the 'on'
#' state (usually estimated from the in situ map)
#'
#' @return A Seurat object, where the posterior of each cell being in the 'on'
#' or 'off' state for each gene is stored in object@@spatial@@mix.probs
#'
#' @importFrom graphics hist lines
#' @importFrom mixtools normalmixEM
#' @importFrom stats dnorm sd quantile
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Note that the PBMC test example object does not contain spatially restricted
#' # examples below are only demonstrate code
#' pmbc_small <- FitGeneK(object = pbmc_small, gene = "MS4A1")
#' }
#'
FitGeneK <- function(
  object,
  gene,
  do.k = 2,
  num.iter = 1,
  do.plot = FALSE,
  genes.use = NULL,
  start.pct = NULL
) {
  data <- object@imputed
  data.use <- data[gene, ]
  names(x = data.use) <- colnames(x = data.use)
  scale.data <- t(x = scale(x = t(x = object@imputed)))
  genes.use <- SetIfNull(x = genes.use, default = rownames(x = scale.data))
  genes.use <- genes.use[genes.use %in% rownames(x = scale.data)]
  scale.data <- scale.data[genes.use, ]
  data.cut <- as.numeric(x = data.use[gene, ])
  cell.ident <- as.numeric(x = cut(x = data.cut, breaks = do.k))
  if (! (is.null(x = start.pct))) {
    cell.ident <- rep.int(x = 1, times = length(x = data.cut))
    cell.ident[data.cut > quantile(x = data.cut, probs = 1 - start.pct)] <- 2
  }
  cell.ident <- order(tapply(
    X = as.numeric(x = data.use),
    INDEX = cell.ident,
    FUN = mean
  ))[cell.ident]
  ident.table <- table(cell.ident)
  if (num.iter > 0) {
    for (i2 in 1:num.iter) {
      cell.ident <- iter.k.fit(
        scale.data = scale.data,
        cell.ident = cell.ident,
        data.use = data.use
      )
      ident.table <- table(cell.ident)
    }
  }
  ident.table <- table(cell.ident)
  raw.probs <- t(
    x = sapply(
      X = data.use,
      FUN = function(y) {
        return(unlist(
          x = lapply(
            X = 1:do.k,
            FUN = function(x) {
              return(
                (ident.table[x] / sum(ident.table)) * dnorm(
                  x = y,
                  mean = mean(x = as.numeric(x = data.use[cell.ident == x])),
                  sd = sd(x = as.numeric(x = data.use[cell.ident == x]))
                )
              )
            }
          )
        ))
      }
    )
  )
  norm.probs <- raw.probs / apply(X = raw.probs, MARGIN = 1, FUN = sum)
  colnames(x = norm.probs) <- unlist(
    x = lapply(
      X = 1:do.k,
      FUN = function(x) {
        paste(gene, x - 1, "post", sep=".")
      }
    )
  )
  norm.probs <- cbind(norm.probs, cell.ident)
  colnames(x = norm.probs)[ncol(x = norm.probs)] <- paste0(gene, ".ident")
  new.mix.probs <- data.frame(
    SubsetColumn(
      data = object@spatial@mix.probs,
      code = paste0(gene, "."),
      invert = TRUE
    ),
    row.names = rownames(x = object@spatial@mix.probs)
  )
  colnames(x = new.mix.probs)[1] <- "nGene"
  object@spatial@mix.probs <- cbind(new.mix.probs, norm.probs)
  if (do.plot) {
    nCol <- 2
    num.row <- floor(x = (do.k + 1) / nCol - (1e-5)) + 1
    hist(
      x = as.numeric(x = data.use),
      probability = TRUE,
      ylim = c(0, 1),
      xlab = gene,
      main = gene
    )
    for (i in 1:do.k) {
      lines(
        x = seq(from = -10, to = 10, by = 0.01),
        y = (ident.table[i] / sum(ident.table)) * dnorm(
          x = seq(from = -10, to = 10, by = 0.01),
          mean = mean(x = as.numeric(x = data.use[cell.ident == i])),
          sd = sd(x = as.numeric(x = data.use[cell.ident == i]))
        ),
        col=i,
        lwd=2
      )
    }
  }
  return(object)
}
