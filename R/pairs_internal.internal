# Find marker pairs
#
# @param id1 Cells for the first grouping (the ones finding the pairs for)
# @param id2 Cells for the second grouping
# @param id3 Cells for the third grouping
# @param genes.list A vector of genes to be considered
# @param training.data A dataframe to train the model on
# @param thr.frac That is the threshold fraction to consider significant
#
# @return A dataframe with the genes pairs that are significantly different
#
FindMarkerPairs <- function(
  id1,
  id2,
  id3,
  genes.list,
  training.data,
  thr.frac
) {
  #   Quality control:
  #   ensure all genes in our list are in our training dataset
  genes.list <- genes.list[genes.list %in% rownames(x = training.data)]
  print(paste('Using', length(x = genes.list), 'genes to create the model'))
  #   A function for testing
  significant.difference <- function(diffs, threshold) {
    return(sum(diffs) >= threshold)
  }
  #   Separate our data by cell
  data1 <- training.data[genes.list, id1]
  data2 <- training.data[genes.list, id2]
  data3 <- training.data[genes.list, id3]
  #   Set the thresholds
  Ngenes <- length(x = genes.list)
  Nthr1 <- ceiling(x = length(x = id1) * thr.frac)
  Nthr2 <- ceiling(x = length(x = id2) * thr.frac)
  Nthr3 <- ceiling(x = length(x = id3) * thr.frac)
  Nthrs <- c(Nthr1, Nthr2, Nthr3)
  #   Find differences
  couple.markers <- matrix(data = -1, nrow = 1, ncol = 2)
  for (gene1 in 1:(Ngenes-1) ) {
    for (gene2 in gene1:Ngenes) {
      diff1 <- (data1[gene1, ] - data1[gene2, ])
      diff2 <- (data2[gene1, ] - data2[gene2, ])
      diff3 <- (data3[gene1, ] - data3[gene2, ])
      #   Check to see if the difference between our two genes is significant
      #   and in which order it the significance occurs
      sig.diffs <- mapply( # Gene1 is greater than Gene2
        FUN = significant.difference,
        diffs = list(
          diff1 > 0,
          diff2 < 0,
          diff3 < 0
        ),
        threshold = Nthrs
      )
      sig.diffs2 <- mapply( # Gene2 is greater than Gene1
        FUN = significant.difference,
        diffs = list(
          diff1 < 0,
          diff2 > 0,
          diff3 > 0
        ),
        threshold = Nthrs
      )
      if (all(sig.diffs)) {
        couple.markers <- rbind(couple.markers, c(gene1, gene2))
      } else if (all(sig.diffs2)) {
        couple.markers <- rbind(couple.markers, c(gene2, gene1))
      }
    }
  }
  couple.markers <- couple.markers[-1, ]
  return(data.frame(
    Gene1 = genes.list[couple.markers[, 1]],
    Gene2 = genes.list[couple.markers[, 2]]
  ))
}

# Classify a whether or not a single cell has enough differences for a set of marker pairs
#
# @param cell Expression data for a single cell
# @param markers Marker pairs to consider
# @param couples.threshold How many pairs of markers must there be a difference
#
# @return The number of significant differences divided by all differences, or NA
#
SingleClassification <- function(cell, markers, couples.threshold) {
  #   Find the difference between every marker in this cell
  differences <- unlist(x = cell[markers[, 1]] - cell[markers[, 2]])
  #   Note: summing a vector of booleans is the same as
  #   subsetting a vector and taking the length of it based
  #   on a condition eg. sum(x > 0) == length(x[x > 0])
  #   except it's about twice as fast
  #   Get the 'hits', where gene1 is greater than gene2
  hits <- sum(differences > 0)
  #   Find all differences
  total.diffs <- sum(differences != 0)
  #   Test to see if we've met our minimum threshold
  if(total.diffs >= couples.threshold) { #min number of couples that can be used
    return(hits / total.diffs)
  } else {
    #   If not, return NA
    print("Not enough pairs where there is a difference in expression, returning 'NA'")
    return(NA)
  }
}

# Randomly pair genes together and ask if the true marker pair list is
# significantly different than the random model
#
# @param cell Expression data for a cell
# @param markers Marker pairs to consider
# @param num.replicates How many replicates for the random model
# @param random.threshold How many cases in the random model must there be to be considered
# @param couples.threshold How many pairs of markers must there be a difference
#
# @return A p-value for this set of marker pairs, or NA
#
RandomSuccess <- function(
  cell,
  markers,
  num.replicates,
  random.threshold,
  couples.threshold
) {
  #   For times 1-N, randomize marker pairs
  #   and ask how many hits we get
  #   (null model)
  null.model <- replicate(
    n = num.replicates,
    expr = Seurat:::SingleClassification(
      cell = Shuffle(x = cell),
      markers = markers,
      couples.threshold = couples.threshold
    )
  )
  #   Remove the NAs
  null.model <- as.numeric(x = na.omit(object = null.model))
  #   If we have too few successes, return NA
  if(length(x = null.model) < random.threshold) {
    print("Too few cases under null model, returning 'NA'")
    return(NA)
  }
  #   Run once without randomly sampling marker pairs
  #   (alternative hypothesis)
  alt.model <- Seurat:::SingleClassification(
    cell = cell,
    markers = markers,
    couples.threshold = couples.threshold
  )
  #   If we get an NA, return an NA
  if(is.na(x = alt.model)) {
    print("Alternative model failed, returning 'NA'")
    return(NA)
  } else {
    #   Return the proportion of times our random pairing
    #   successes were less than the real run,
    #   the smaller the better (p-value)
    return(sum(null.model < alt.model) / length(x = null.model))
  }
}

# Assign a cell-cycle score
#
# @param data A matrix of expression data in cell by gene format
# @param g1.marker.pairs Marker pairs for G1 phase
# @param g2m.marker.pairs Marker pairs for G2M phase
# @param s.marker.pairs Marker pairs for S phase
# @param genes.training Genes that were used for training
# @param num.replicates How many replicates for the null model
# @param random.threshold How many cases in the random model must there be to be considered
# @param couples.threshold How many pairs of markers must there be a difference
# @param num.cores How many cores should we use?
#
# @return A list of two dataframes
# @return \item{scores} Regular cell-cyle scores
# @return \tiem{scores.norm} Normalized cell-cycle scores
#
AssignScore <- function(
  data,
  g1.marker.pairs,
  g2m.marker.pairs,
  s.marker.pairs,
  genes.training,
  num.replicates = 1000,
  random.threshold = 100,
  couples.threshold = 10,
  num.cores = 1
) {
  #   Quality control
  #   ensure that 'couples.threshold' is less than the smallest number of pairs
  #   otherwise, this thing won't work
  num.pairs <- vapply(
    X = list(g1.marker.pairs, g2m.marker.pairs, s.marker.pairs),
    FUN = nrow,
    FUN.VALUE = numeric(1)
  )
  if (any(num.pairs < couples.threshold)) {
    couples.threshold <- min(num.pairs)
    warning(
      paste0(
        "There are fewer marker pairs than 'couples.threshold', using ",
        couples.threshold,
        ", the minimum number of marker pairs as 'couples.threshold', instead"
      )
    )
  }
  #select genes
  data <- data[intersect(x = genes.training, y = rownames(x = data)), ]
  genes.data <- rownames(x = data)
  markers <- unique(x = as.character(x = unlist(x = c(
    g1.marker.pairs[, 1:2],
    g2m.marker.pairs[, 1:2],
    s.marker.pairs[, 1:2]
  ))))
  if (all(markers %in% genes.data)) {
    #checking if markers are in the list of genes of the testing dataset
    g1.markers <- g1.marker.pairs
    g2m.markers <- g2m.marker.pairs
    s.markers <- s.marker.pairs
  } else {
    #remove markers that are not in the dataset
    genes.to.remove <- markers[!(markers %in% genes.data)]
    g1.markers <- RemoveFromTable(
      to.remove = genes.to.remove,
      data = g1.marker.pairs
    )
    g2m.markers <- RemoveFromTable(
      to.remove = genes.to.remove,
      data = g2m.marker.pairs
    )
    s.markers <- RemoveFromTable(
      to.remove = genes.to.remove,
      data = s.marker.pairs
    )
  }
  #   Info messages
  print(sprintf("Number of G1 pairs: %d", nrow(x = g1.markers)))
  print(sprintf("Number of G2M pairs: %d", nrow(x = g2m.markers)))
  print(sprintf("Number of S pairs: %d", nrow(x = s.markers)))
  #run the allocation algorithm
  #   Applied by cell
  if (num.cores == 0) {
    num.cores <- parallel::detectCores()
  }
  sfInit(parallel = TRUE, cpus = num.cores, type = 'SOCK')
  sfExportAll()
  print("Calculcating G1-phase scores")
  scores.g1 <- sfApply(
    x = data,
    margin = 2,
    fun = RandomSuccess,
    markers = g1.markers,
    num.replicates = num.replicates,
    random.threshold = random.threshold,
    couples.threshold = couples.threshold
  )
  print("Calculcating G2M-phase scores")
  scores.g2m <- sfApply(
    x = data,
    margin = 2,
    fun = RandomSuccess,
    markers = g2m.markers,
    num.replicates = num.replicates,
    random.threshold = random.threshold,
    couples.threshold = couples.threshold
  )
  print("Calculcating S-phase scores")
  scores.s <- sfApply(
    x = data,
    margin = 2,
    fun = RandomSuccess,
    markers = s.markers,
    num.replicates = num.replicates,
    random.threshold = random.threshold,
    couples.threshold = couples.threshold
  )
  sfStop()
  #   Assemble the scores and normalize
  scores <- data.frame(
    g1.score = scores.g1,
    g2m.score = scores.g2m,
    s.score = scores.s
  )
  scores.normalised <- t(
    x = apply(
      X = scores,
      MARGIN = 1,
      FUN = function(x) (x)/sum(x)
    )
  )
  #   Assemble our results and return
  results <- list(
    "scores" = scores,
    "scores.norm" = scores.normalised
  )
  return(results)
}

# Determine which score is most significant
#
# @param cell.scores A two-length vector of cell-cycle scores
# @param first Which phase represents cell.scores[1]
# @param second Which phase represents cell.scores[2]
# @param null Which phase if neither score is significant
# @param score.threshold The threshold in which all of cell.scores must be greater than
# to count as signifcant. If either score is below score.threshold, null is assigned
#
# @return first, second, or null, depending on which phase is deemed most significant
#
AllocateScores <- function(
  cell.scores,
  first = 'G1',
  second = 'G2M',
  null = 'S',
  score.threshold = 0.5
) {
  #   Quality control
  num.scores <- length(x = cell.scores)
  if (num.scores > 2) {
    warning("'cell.scores' is longer than 2, only considering first two scores")
    cell.scores <- cell.scores[1:2]
  } else if (num.scores < 2) {
    stop("Too few scores to compare")
  }
  #   Check to see neither score is above our threshold
  #   If so, assign the 'null' phase
  #   Otherwise, check to see which phase has the greater score
  null.phase <- vapply(
    X = cell.scores,
    FUN = function(score) {
      return(score < score.threshold)
    },
    FUN.VALUE = logical(1)
  )
  if (all(null.phase)) {
    return(null)
  } else if (cell.scores[1] >= cell.scores[2]) {
    if (cell.scores[1] == cell.scores[2]) {
      warning("Scores are equal, returning first allocation")
    }
    return(first)
  } else {
    return(second)
  }
}

# Assemble the list of score dataframes into a single dataframe
#
# @param results The results from AssignScore
#
# @return A dataframe with regular and normalized scores
#
AssembleScores <- function(results) {
  scores <- data.frame(
    G1 = results$scores[, 'g1.score'],
    G2M = results$scores[, 'g2m.score'],
    S = results$scores[, 's.score'],
    G1.norm = results$scores.norm[, 'g1.score'],
    G2M.norm = results$scores.norm[, 'g2m.score'],
    S.norm = results$scores.norm[, 's.score']
  )
  return(scores)
}
