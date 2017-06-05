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
FindGeneMarkerPairs <- function(object, g1.cells, g2m.cells, s.cells, genes.list, threshold.fraction = 0.5) {
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
#' @param num.cores The number of cores to run on
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
    num.cores = NULL
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
    object <- Seurat::AddMetaData(object = object, metadata = scores)
    return(object)
}

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

#' Shuffle a vector
#' @param x A vector
#' @return A vector with the same values of x, just in random order
#' @export
#'
Shuffle <- function(x) {
    return(x[base::sample.int(n = base::length(x = x), size = base::length(x = x), replace = FALSE)])
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

RandomSuccess <- function(cell, markers, num.replicates, random.threshold, couples.threshold) {
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

AssignScore <- function(
    data,
    g1.marker.pairs,
    g2m.marker.pairs,
    s.marker.pairs,
    genes.training,
    num.replicates = 1000,
    random.threshold = 100,
    couples.threshold = 10,
    num.cores = NULL
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
    markers <- unique(
        x = as.character(
            x = unlist(
                x = c(
                    g1.marker.pairs[, 1:2],
                    g2m.marker.pairs[, 1:2],
                    s.marker.pairs[, 1:2]
                )
            )
        )
    )
    if (all(markers %in% genes.data)) {
        #checking if markers are in the list of genes of the testing dataset
        g1.markers <- g1.marker.pairs
        g2m.markers <- g2m.marker.pairs
        s.markers <- s.marker.pairs
    } else {
        #remove markers that are not in the dataset
        genes.to.remove <- markers[!(markers %in% genes.data)]
        g1.markers <- RemoveFromTable(to.remove = genes.to.remove, data = g1.marker.pairs)
        g2m.markers <- RemoveFromTable(to.remove = genes.to.remove, data = g2m.marker.pairs)
        s.markers <- RemoveFromTable(to.remove = genes.to.remove, data = s.marker.pairs)
    }
    #   Info messages
    print(sprintf("Number of G1 pairs: %d", nrow(x = g1.markers)))
    print(sprintf("Number of G2M pairs: %d", nrow(x = g2m.markers)))
    print(sprintf("Number of S pairs: %d", nrow(x = s.markers)))
    #run the allocation algorithm
    #   Applied by cell
    if (is.null(x = num.cores)) {
        num.cores <- parallel::detectCores()
    }
    sfInit(parallel = TRUE, cpus = num.cores, type = 'SOCK')
    sfExportAll()
    scores.g1 <- sfApply(
        x = data,
        margin = 2,
        fun = RandomSuccess,
        markers = g1.markers,
        num.replicates = num.replicates,
        random.threshold = random.threshold,
        couples.threshold = couples.threshold
    )
    scores.g2m <- sfApply(
        x = data,
        margin = 2,
        fun = RandomSuccess,
        markers = g2m.markers,
        num.replicates = num.replicates,
        random.threshold = random.threshold,
        couples.threshold = couples.threshold
    )
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

AllocateScores <- function(cell.scores, first = 'G1', second = 'G2M', null = 'S', score.threshold = 0.5) {
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
        FUN = function(score) score < score.threshold,
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

plot.scores <- function(scores, truths = FALSE, do.return = FALSE) {
    if (truths) {
        stitle <- 'Score Assignment Truths'
        ntitle <- 'Normalized Score Truths'
    } else {
        stitle <- 'G1 vs G2M scores'
        ntitle <- 'G1 vs G2M normalized scores'
    }
    pscores <- ggplot2::ggplot(
        data = scores,
        mapping = ggplot2::aes(x = G1, y = G2M)
    )
    if (truths) {
        pscores <- pscores + ggplot2::geom_point(
            size = 1,
            mapping = ggplot2::aes(color = Match)
        )
    } else {
        pscores <- pscores + ggplot2::geom_point(
            size = 1,
            mapping = ggplot2::aes(color = Assignment)
        )
    }
    pscores <- pscores + ggplot2::scale_color_discrete(name = '')
    pscores <- pscores + ggplot2::labs(title = stitle)
    pscores <- pscores + ggplot2::theme(legend.position = 'bottom')
    pnorm <- ggplot2::ggplot(
        data = scores,
        mapping = ggplot2::aes(x = G1.norm, y = G2M.norm)
    )
    if (truths) {
        pnorm <- pnorm + ggplot2::geom_point(
            size = 1,
            mapping = ggplot2::aes(color = Norm.match)
        )
    } else {
        pnorm <- pnorm + ggplot2::geom_point(
            size = 1,
            mapping = ggplot2::aes(color = Norm.assign)
        )
    }
    pnorm <- pnorm + ggplot2::scale_color_discrete(name = '')
    pnorm <- pnorm + ggplot2::labs(
        title = ntitle,
        x = 'G1 Normalized',
        y = 'G2M Normalized'
    )
    pnorm <- pnorm + ggplot2::theme(legend.position = 'bottom')
    if (do.return) {
        return(list('scores' = pscores, 'norm' = pnorm))
    } else {
        gridExtra::grid.arrange(pscores, pnorm, ncol = 2)
    }
}

pie.chart <- function(slices, ...) {
    percents <- round(x = slices / sum(slices) * 100, 2)
    if (is.null(x = names(x = slices))) {
        labels <- as.character(x = percents)
    } else {
        labels <- paste0(names(x = slices), ': ', percents, '%')
    }
    pie(x = slices, labels = labels, ...)
}
