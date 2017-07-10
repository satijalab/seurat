# Find gene terms from Enrichr
#
# @param QueryGene A gene to query on Enrichr
#
# @return An XML document with information on the queried gene
#
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
