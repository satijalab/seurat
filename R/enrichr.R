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
#' @return Returns a Seurat object with enrichment scores added to object@data.info
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
  n.bin = 20,
  seed.use = 1,
  ctrl.size = 20,
  use.k = FALSE,
  enrich.name = "Cluster"
) {
  if (use.k) {
    genes.k <- names(x = object@kmeans.obj[[1]]$cluster)
    genes.list <- list()
    for (i in as.numeric(x = names(x = table(object@kmeans.obj[[1]]$cluster)))) {
      genes.list[[i]] <- names(x = which(x = object@kmeans.obj[[1]]$cluster == i))
    }
    cluster.length <- length(x = genes.list)
  } else {
    if (is.null(x = genes.list)) {
      stop("Missing input gene list")
    }
    genes.k <- ainb(
      a = unique(x = unlist(x = genes.list)),
      b = rownames(x = object@data)
    )
    cluster.length <- length(x = genes.list)
  }
  data.avg <- apply(X = object@data[genes.k, ], MARGIN = 1, FUN = mean)
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
    EnrichrID <- AddGeneList(GeneList = GeneList)[[2]]
  }
  Enrichrlibs <- listEnrichrDbs()$libraryName
  if (! GeneSetLib %in% Enrichrlibs) {
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
    ViewGeneList(EnrichrID = EnrichrID)
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

PrintLibs <- function() {
  return (listEnrichrDbs()$libraryName)
}

AddGeneList <- function(GeneList = NULL, Description = "Example gene list") {
  path.use <- "Enrichr/addList"
  if (is.null(x = GeneList)) {
    stop("Missing gene list")
  }
  genes.use <- paste(GeneList, collapse = "\n")
  query.use <- list(list = genes.use, description = Description)
  api.post <- POST(
    url = "http://amp.pharm.mssm.edu/",
    path = path.use,
    body = query.use
  )
  api.status <- status_code(x = api.post)
  if (api.status != 200) {
    stop("Error analyzing gene list")
  }
  api.data <- content(x = api.post, as = "text")
  enrichr.id <- extract_field(extract_field(api.data, 3, "\n"), 2, ": ")
  return(list(api.data, enrichr.id))
}

ViewGeneList <- function(EnrichrID) {
  path.use <- "Enrichr/view"
  api.get <- GET(
    url = "http://amp.pharm.mssm.edu/",
    path = path.use,
    query = list(userListId = EnrichrID)
  )
  api.status <- status_code(x = api.get)
  if (api.status != 200) {
    stop("Error getting gene list")
  }
  api.data <- content(x = api.get)$genes
  print(sort(x = unlist(x = api.data)))
}

FindGeneTerms <- function(QueryGene = NULL) {
  if (is.null(x = QueryGene)) {
    stop("Missing query gene")
  }
  path.use <- "Enrichr/genemap"
  api.get <- GET(
    url = "http://amp.pharm.mssm.edu/",
    path = path.use,
    query = list(gene = QueryGene)
  )
  api.status <- status_code(x = api.get)
  if (api.status != 200) {
    stop("Error searching for terms")
  }
  api.data <- content(x = api.get)
  return (api.data)
}
