# Add genes to Enrichr to analyze
#
# @param GeneList Input gene list
# @param Description A description of the gene list
#
# @return A list with the following values
# @return \item{data} The data from analyzing the gene list
# @return \item{id} The enrichr ID
#
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
  enrichr.id <- ExtractField(ExtractField(api.data, 3, "\n"), 2, ": ")
  return(list(data = api.data, id = enrichr.id))
}

# View a list of genes given an Enrichr ID
#
# @param EnrichrID The Enrichr ID to view
#
# @return A vector of genes for this Enrichr ID
#
GetGeneList <- function(EnrichrID) {
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
  return(sort(x = unlist(x = api.data)))
}
