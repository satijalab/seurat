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

#' View Enrichr libraries
#'
#' @return A vector of Enrichr Libraries
#'
EnrichrLibs <- function() {
  return (listEnrichrDbs()$libraryName)
}
