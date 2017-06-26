
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

PrintLibs=function(){
  return (listEnrichrDbs()$libraryName)
}

RunEnrichr=function(GeneList=NULL,EnrichrID=NULL,GeneSetLib="GO_Biological_Process_2017",PrintList=FALSE,PrintID=FALSE,Download=FALSE,FileName="test.txt",...){
  
  if (is.null(EnrichrID)) EnrichrID=AddGeneList(GeneList)[[2]]
  
  Enrichrlibs=listEnrichrDbs()$libraryName
  if (! GeneSetLib %in% Enrichrlibs) stop("Error getting Enrichr libraries")
  
  path.use="Enrichr/export"
  api.get=GET(url = "http://amp.pharm.mssm.edu/",path = path.use,query = list(userListId = EnrichrID, filename = FileName, backgroundType = GeneSetLib),stream=TRUE)
  api.status=status_code(api.get)
  bin=content(api.get,"raw")
  writeBin(bin,FileName)
  
  if (api.status!=200) stop("Error fetching enrichment results")
  
  if (PrintList) ViewGeneList(EnrichrID)
  if (PrintID) print (EnrichrID)
  
  api.data=read.delim(FileName,header = T,sep = "\t")
  if (!Download) file.remove(FileName);return (api.data)

}

AddGeneList=function(GeneList=NULL, Description = "Example gene list"){
  path.use="Enrichr/addList"
  if (is.null(GeneList)) stop("Missing gene list")
  genes.use=paste(GeneList,collapse = "\n")
  query.use=list(list = genes.use, description = Description)
  
  api.post=POST(url = "http://amp.pharm.mssm.edu/", path = path.use, body = query.use)
  
  api.status=status_code(api.post)
  if (api.status!=200) stop("Error analyzing gene list")
  api.data=content(api.post,"text")
  enrichr.id=extract_field(extract_field(api.data,3,"\n"),2,": ")
  
  return(list(api.data,enrichr.id))
}

ViewGeneList=function(EnrichrID){
  path.use="Enrichr/view"
  api.get=GET(url = "http://amp.pharm.mssm.edu/",path = path.use,query = list(userListId = EnrichrID))
  
  api.status=status_code(api.get)
  if (api.status!=200) stop("Error getting gene list")
  api.data=content(api.get)$genes
  print (sort(unlist(api.data)))
}

FindGeneTerms=function(QueryGene=NULL){
  
  if (is.null(QueryGene)) stop("Missing query gene")
  
  path.use="Enrichr/genemap"
  api.get=GET(url = "http://amp.pharm.mssm.edu/",path = path.use,query = list(gene = QueryGene))
  api.status=status_code(api.get)
  if (api.status!=200) stop("Error searching for terms")
  
  api.data=content(api.get)
  return (api.data)
}


