#' Calculate enrichment scores for gene expression programs in single cells
#'
#' Calculate the average expression levels of each program (cluster) on single cell level,
#' subtracted by the aggregated expression of control gene sets.
#' All analyzed genes are binned based on averaged expression, and the control genes are
#' randomly selected from each bin.
#'
#' @param object Seurat object
#' @param genes.list Gene expression programs in list
#' @param n.bin Number of bins of aggregate expression levels for all analyzed genes
#' @param seed.use Random seed for sampling
#' @param ctrl.size Number of control genes selected from the same bin per analyzed gene
#' @param use.k Use gene clusters returned from DoKMeans()
#' @param enrich.name Name for the expression programs
#'
#' @return Returns a Seurat object with enrichment scores added to object@meta.data
#'
#' @importFrom Hmisc cut2
#'
#' @references Tirosh et al, Science (2016)
#'
#' @export
#'
AddEnrichScore <- function(
  object,
  genes.list = NULL,
  genes.pool = NULL,
  n.bin = 25,
  seed.use = 1,
  ctrl.size = 100,
  use.k = FALSE,
  enrich.name = "Cluster"
) {
  if (use.k) {
    genes.list <- list()
    for (i in as.numeric(x = names(x = table(object@kmeans.obj[[1]]$cluster)))) {
      genes.list[[i]] <- names(x = which(x = object@kmeans.obj[[1]]$cluster == i))
    }
    cluster.length <- length(x = genes.list)
  } else {
    if (is.null(x = genes.list)) {
      stop("Missing input gene list")
    }
    genes.list <- lapply(
      X = genes.list,
      FUN = function(x) {
        return(intersect(x = x, y = rownames(x = object@data)))
      }
    )
    cluster.length <- length(x = genes.list)
  }
  if (is.null(x = genes.pool)) genes.pool = rownames(object@data)
  data.avg <- apply(X = object@data[genes.pool,], MARGIN = 1, FUN = mean)
  data.avg <- data.avg[order(data.avg)]
  data.cut <- as.numeric(x = cut2(
    x = data.avg,
    m = round(x = length(x = data.avg) / n.bin)
  ))
  names(x = data.cut) <- names(x = data.avg)
  ctrl.use <- vector("list", cluster.length)
  for (i in 1:cluster.length) {
    genes.use <- genes.list[[i]]
    for (j in 1:length(x = genes.use)) {
      ctrl.use[[i]] <- c(
        ctrl.use[[i]],
        names(x = sample(
          x = data.cut[which(x = data.cut == data.cut[genes.use[j]])],
          size = ctrl.size,
          replace = FALSE
        ))
      )
    }
  }
  ctrl.use <- lapply(X = ctrl.use, FUN = unique)
  ctrl.scores <- c()
  for (i in 1:length(ctrl.use)) {
    genes.use <- ctrl.use[[i]]
    ctrl.scores <- rbind(
      ctrl.scores,
      apply(X = object@data[genes.use, ], MARGIN = 2, FUN = mean)
    )
  }
  genes.scores <- c()
  for (i in 1:cluster.length) {
    genes.use <- genes.list[[i]]
    genes.scores <- rbind(
      genes.scores,
      apply(X = object@data[genes.use, ], MARGIN = 2, FUN = mean)
    )
  }
  
  genes.scores.use <- genes.scores - ctrl.scores
  rownames(x = genes.scores.use) <- paste0(enrich.name, 1:cluster.length)
  genes.scores.use <- t(x = as.data.frame(x = genes.scores.use))
  object <- AddMetaData(
    object = object,
    metadata = genes.scores.use,
    col.name = colnames(x = genes.scores.use)
  )
  return (object)
}

#' Run Enrichr
#'
#' Run Enrichr in R: enrichment against available libraries in Enrichr
#'
#' @param GeneList Input gene list
#' @param EnrichrID Enrichr identifier for input gene list, if input genes are already uploaded
#' @param GeneSetLib Library to be enriched against
#' @param PrintList Print gene list
#' @param PrintID Print Enrichr identifier
#' @param Download Download enrichment results as a txt file
#' @param FileName Name for the txt file
#'
#' @return If Download=TRUE, a tab-delimited file is written in the directory specified in FileName
#'
#' @import httr enrichR
#'
#' @export
#'

RunEnrichr <- function(
  GeneList = NULL,
  EnrichrID = NULL,
  GeneSetLib = "GO_Biological_Process_2017",
  PrintList = FALSE,
  PrintID = FALSE,
  Download = FALSE,
  FileName = "test.txt",
  ...
) {
  if (is.null(x = EnrichrID)) {
    EnrichrID <- AddGeneList(GeneList = GeneList)$id
  }
  if (! GeneSetLib %in% EnrichrLibs()) {
    stop("Error getting Enrichr libraries")
  }
  path.use <- "Enrichr/export"
  api.get <- GET(
    url = "http://amp.pharm.mssm.edu/",
    path = path.use,
    query = list(
      userListId = EnrichrID,
      filename = FileName,
      backgroundType = GeneSetLib
    ),
    stream = TRUE
  )
  api.status = status_code(x = api.get)
  bin <- content(x = api.get, as = "raw")
  writeBin(object = bin, con = FileName)
  if (api.status != 200) {
    stop("Error fetching enrichment results")
  }
  if (PrintList) {
    print(GetGeneList(EnrichrID = EnrichrID))
  }
  if (PrintID) {
    print(EnrichrID)
  }
  api.data <- read.delim(file = FileName, header = TRUE, sep = "\t")
  if (! Download) {
    file.remove(FileName)
    return (api.data)
  }
}

#' Enrich for cell cycle phases
#'
#' @param object A Seurat object
#' @param g2m.genes A vector of genes associated with G2M phase
#' @param s.genes A vector of genes associated with S phases
#' @param set.ident If true, sets identity to phase assignments
#' Stashes old identities in 'old.ident'
#'
#' @return A Seurat object with the following columns added to object@meta.data
#' @return \list{\code{S.Score}} S phase score
#' @return \list{\code{G2M.Score}} G2M phase score
#' @return \list{\code{Phase}} Cell-cycle phase assignment
#'
#' @seealso \link{\code{AddEnrichScore}}
#'
#' @export
#'
CellCycleEnrichment <- function(
  object,
  g2m.genes,
  s.genes,
  set.ident = FALSE
) {
  enrich.name <- 'Cell Cycle'
  genes.list <- list('S.Score' = s.genes, 'G2M.Score' = g2m.genes)
  object.cc <- AddEnrichScore(
    object = object,
    genes.list = genes.list,
    enrich.name = enrich.name,
    ctrl.size = min(vapply(X = genes.list, FUN = length, FUN.VALUE = numeric(1)))
  )
  cc.columns <- grep(pattern = enrich.name, x = colnames(x = object.cc@meta.data))
  cc.scores <- object.cc@meta.data[, cc.columns]
  assignments <- apply(
    X = cc.scores,
    MARGIN = 1,
    FUN = function(scores, first = 'S', second = 'G2M', null = 'G1') {
      if (all(scores < 0)) {
        return(null)
      } else {
        return(c(first, second)[which(x = scores == max(scores))])
      }
    }
  )
  cc.scores <- merge(x = cc.scores, y = data.frame(assignments), by = 0)
  colnames(x = cc.scores) <- c('rownames', 'S.Score', 'G2M.Score', 'Phase')
  rownames(x = cc.scores) <- cc.scores$rownames
  cc.scores <- cc.scores[, c('S.Score', 'G2M.Score', 'Phase')]
  object <- AddMetaData(object = object, metadata = cc.scores)
  if (set.ident) {
    object <- StashIdent(object = object, save.name = 'old.ident')
    object <- SetAllIdent(object = object, id = 'Phase')
  }
  return(object)
}

#' View Enrichr libraries
#'
#' @return A vector of Enrichr Libraries
#'
EnrichrLibs <- function() {
  return (listEnrichrDbs()$libraryName)
}
