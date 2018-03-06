#' @include seurat.R
NULL

tcc <- setClass(
  Class = "tcc",
  slots = list(
    tcc.raw = "dgCMatrix",
    tcc.norm = "dgCMatrix",
    ec.to.tid.map = "environment",
    tid.to.ec.map = "environment",
    gene.map = "data.frame"
  )
)

# Function to perform mapping functions for TCC data. To specify to or from
# types as EC (equivalence class), TIDX (transcript id index), TID (transcript
# id name), GENE (gene name), MAT (matrix with full mapping info)
#
# @param from        Value to map from - should be a string.
# @param from.type   Type to map from - one of the following: EC, TIDX, TID, GENE
# @param to          Type to map to - one of the following: EC, TIDX, TID, GENE,
#                    MAT
#
# @return returns either a vector of mapped values or matrix with mapping info
#
TCCMap <- function(object, from, from.type, to) {
  if(length(from) > 1){
    stop("from must be of length 1")
  }
  if(from.type == "EC") {
    ec <- from
    transcript.ids <- as.character(unname(sapply(X = from, FUN = function(x) HashTableLookup(x, object@tcc@ec.to.tid.map))))
    if(to == "TIDX") {
      return(transcript.ids)
    }
    ec <- rep(ec, length(transcript.ids))
    from <- transcript.ids
    from.type <- "TIDX"
  }
  if(from.type == "TIDX") {
    transcript.ids <- from
    if(to == "EC"){
      ec <- as.vector(sapply(X = transcript.ids, FUN = function(x) HashTableLookup(x, object@tcc@tid.to.ec.map)))
    }
    transcript.names <- object@tcc@gene.map[as.numeric(transcript.ids) + 1, 1]
    if (to == "TID") {
      return(transcript.names)
    }
    from <- transcript.names
    from.type <- "TID"
  }
  if(from.type == "TID") {
    transcript.names <- from
    if (! exists("transcript.ids")){
      transcript.ids <- which(object@tcc@gene.map[, 1] %in% transcript.names) - 1
    }
    if(to == "TIDX") {
      return(transcript.ids)
    }
    if(to == "EC" || to == "MAT"){
      if(! exists("ec")){
        ec <- as.vector(sapply(X = as.character(transcript.ids), FUN = function(x) HashTableLookup(x, object@tcc@tid.to.ec.map)))
      }
      if(to == "EC") {
        return(ec)
      }
    }
    gene.names <- object@tcc@gene.map[as.numeric(transcript.ids) + 1, 2]
    if(to == "GENE") {
      return(gene.names)
    }
  }
  if(from.type == "GENE") {
    gene.names <- from
    transcript.ids <- which(object@tcc@gene.map[, 2] %in% gene.names)
    transcript.names <- object@tcc@gene.map[transcript.ids, 1]
    if(to == "TID") {
      return(transcript.names)
    }
    if(to == "TIDX") {
      return(transcript.ids)
    }
    ec <- as.vector(sapply(X = transcript.ids, FUN = function(x) HashTableLookup(x, object@tcc@tid.to.ec.map)))
    if(to == "EC") {
      return(ec)
    }
  }
  map.mat <- data.frame(EC = ec,
                        TIDX = transcript.ids,
                        TID = transcript.names,
                        GENE = gene.names,
                        stringsAsFactors = FALSE)
  return(map.mat)
}
