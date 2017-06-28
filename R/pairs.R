#' Find marker pairs for cell-cycle phase prediction
#'
#' @param object A Seurat object
#' @param g1.cells A vector of cells known to be in G1 phase
#' @param g2m.cells A vector of cells known to be in G2M phase
#' @param s.cells A vector of cells known to be in S phase
#' @param genes.list A vector of genes to train the model on
#' @param threshold.fraction Fraction of differences needed to be significant
#'
#' @return A list with the following fields
#' @return \item{g1.pairs}{A data frame of gene pairs marking G1 phase}
#' @return \item{g2m.pairs}{A data frame of gene pairs marking G2M phase}
#' @return \item{s.pairs}{A data frame of gene pairs marking S phase}
#'
#' @references Scialdone A, Natarajan KN, Saraiva LR, Proserpio V, Teichmann SA, Stegle O, Marioni JC, Buettner F. (2015) Computational assignment of cell-cycle stage from single-cell transcriptome data. Methods
#'
#' @export
#'
FindGeneMarkerPairs <- function(
  object,
  g1.cells,
  g2m.cells,
  s.cells,
  genes.list,
  threshold.fraction = 0.5
) {
  data.use <- object@scale.data
  #   Quality control
  cells.data <- colnames(x = data.use)
  g1.cells <- base::intersect(x = g1.cells, y = cells.data)
  g2m.cells <- base::intersect(x = g2m.cells, y = cells.data)
  s.cells <- base::intersect(x = s.cells, y = cells.data)
  genes.list <- base::intersect(x = genes.list, y = rownames(x = data.use))
  #   Find marker pairs
  g1.pairs <- FindMarkerPairs(
    id1 = g1.cells,
    id2 = g2m.cells,
    id3 = s.cells,
    genes.list = genes.list,
    training.data = data.use,
    thr.frac = threshold.fraction
  )
  g2m.pairs <- FindMarkerPairs(
    id1 = g2m.cells,
    id2 = s.cells,
    id3 = g1.cells,
    genes.list = genes.list,
    training.data = data.use,
    thr.frac = threshold.fraction
  )
  s.pairs <- FindMarkerPairs(
    id1 = s.cells,
    id2 = g1.cells,
    id3 = g2m.cells,
    genes.list = genes.list,
    training.data = data.use,
    thr.frac = threshold.fraction
  )
  return(list(
    'g1.pairs' = g1.pairs,
    'g2m.pairs' = g2m.pairs,
    's.pairs' = s.pairs
  ))
}

#' Predict cell cycle phases
#'
#' @param object A Seurat object
#' @param g1.pairs A data frame of gene pairs marking G1 phase
#' @param g2m.pairs A data frame of gene pairs marking G2M phase
#' @param s.pairs A data frame of gene pairs marking S phase
#' @param genes.training A list of genes that the model was trained on,
#' will assemble from pairs data frames if not provided
#' @param null.phase Which phase is considered the 'null' phase when allocating scores
#' ie. which phase is not compared for score allocation. See Scialdone et al (2015) for more details
#' @param score.threshold Threshold for comparing scores. If neither comapred score is above this
#' threshold, then the 'null' phase is assigned.
#' @param num.replicates The number of times we should repeat our testing
#' @param random.threshold The minimum number of times a randomly selected
#' pair of genes provides a significant hit
#' @param couples.threshold The minimum number of gene pairs that provide a significant hit
#' @param num.cores The number of cores to run on, set as 0 to use all available cores
#' @param set.ident Set the identity of the new object to cell-cycle assignment,
#' will stash old identity as 'old.ident'
#'
#' @return A Seurat object with cell-cycle scores and allocations added to object@data.info
#'
#' @import snowfall
#' @references Scialdone A, Natarajan KN, Saraiva LR, Proserpio V, Teichmann SA, Stegle O, Marioni JC, Buettner F. (2015) Computational assignment of cell-cycle stage from single-cell transcriptome data. Methods
#'
#' @export
#'
PredictCellCyclePhases <- function(
  object,
  g1.pairs,
  g2m.pairs,
  s.pairs,
  genes.training = NULL,
  null.phase = 'G1',
  score.threshold = 0.5,
  num.replicates = 1000,
  random.threshold = 100,
  couples.threshold = 10,
  num.cores = 1,
  set.ident = FALSE
) {
  expression.data <- object@scale.data
  #   Quality control
  #   Figure out what the null phase is and
  #   what scores columns we use for score allocation
  null.phase <- base::toupper(x = null.phase)
  scores.allocation <- switch(
    EXPR = null.phase,
    'G1' = list(first = 'G2M', second = 'S', scores.cols = c(2, 3), norm.cols = c(5, 6)),
    'G2M' = list(first = 'G1', second = 'S', scores.cols = c(1, 3), norm.cols = c(4, 6)),
    'S' = list(first = 'G1', second = 'G2M', scores.cols = c(1, 2), norm.cols = c(4, 5)),
    stop("'null' must be one of 'G1', 'G2M' or 'S'")
  )
  scores.allocation$null <- null.phase
  #   If no genes.training was provided,
  #   make it out of the unique genes from
  #   the pairs data frames
  if (is.null(x = genes.training)) {
    genes.training <- unique(
      x = as.character(
        x = unlist(
          x = c(
            g1.pairs[, 1:2],
            g2m.pairs[, 1:2],
            s.pairs[, 1:2]
          )
        )
      )
    )
  }
  #   Get scores
  tryCatch(
    expr = scores <- AssignScore(
      data = expression.data,
      g1.marker.pairs = g1.pairs,
      g2m.marker.pairs = g2m.pairs,
      s.marker.pairs = s.pairs,
      genes.training = genes.training,
      num.replicates = num.replicates,
      random.threshold = random.threshold,
      couples.threshold = couples.threshold,
      num.cores = num.cores
    ),
    interrupt = function(condition) {
      sfStop()
      stop('killed')
    }
  )
  #   Assemble the scores
  scores <- AssembleScores(results = scores)
  #   Allocate the scores to get cell-cycle phases
  scores$assignment <- apply(
    X = scores[, scores.allocation$scores.cols],
    MARGIN = 1,
    FUN = AllocateScores,
    first = scores.allocation$first,
    second = scores.allocation$second,
    null = scores.allocation$null,
    score.threshold = score.threshold
  )
  scores$norm.assignment <- apply(
    X = scores[, scores.allocation$norm.cols],
    MARGIN = 1,
    FUN = AllocateScores,
    first = scores.allocation$first,
    second = scores.allocation$second,
    null = scores.allocation$null,
    score.threshold = score.threshold
  )
  #   Add the cell-cycle scores and assignment information to object@data.info
  object <- Seurat::AddMetaData(object = object, metadata = scores)
  #   If we're setting the identity to cell-cycle phase
  #   first, stash the current identity as 'old.ident'
  #   then, set the identity to the normalized assignment
  if (set.ident) {
    object <- Seurat::StashIdent(object = object, save.name = 'old.ident')
    object <- Seurat::SetAllIdent(object = object, id = 'norm.assignment')
  }
  return(object)
}


#' Shuffle a vector
#' @param x A vector
#' @return A vector with the same values of x, just in random order
#' @export
#'
Shuffle <- function(x) {
  return(x[base::sample.int(
    n = base::length(x = x),
    size = base::length(x = x),
    replace = FALSE
  )])
}

#' Remove data from a table
#'
#' This function will remove any rows from a data frame or matrix
#' that contain certain values
#'
#' @param to.remove A vector of values that indicate removal
#' @param data A data frame or matrix
#'
#' @return A data frame or matrix with values removed by row
#'
#' @export
#'
RemoveFromTable <- function(to.remove, data) {
  remove.indecies <- apply(
    X = data,
    MARGIN = 2,
    FUN = function(col) {
      return(which(x = col %in% to.remove))
    }
  )
  remove.indecies <- unlist(x = remove.indecies)
  remove.indecies <- as.numeric(x = remove.indecies)
  if (length(x = remove.indecies) == 0) {
    return(data)
  } else {
    return(data[-remove.indecies, ])
  }
}
