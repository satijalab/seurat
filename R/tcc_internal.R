#' @include seurat.R
NULL

tcc <- setClass(
  Class = "tcc",
  slots = list(
    tcc.raw = "dgCMatrix",
    tcc.norm = "dgCMatrix",
    tx.raw = "dgCMatrix",
    ec.to.tid.map = "environment",
    tid.to.ec.map = "environment",
    gene.map = "matrix"
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
      return(unname(transcript.names))
    }
    from <- transcript.names
    from.type <- "TID"
  }
  if(from.type == "TID") {
    transcript.names <- from
    if (! exists("transcript.ids")){
      transcript.ids <- as.character(which(object@tcc@gene.map[, 1] %in% transcript.names) - 1)
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
      return(unique(unname(gene.names)))
    }
  }
  if(from.type == "GENE") {
    gene.names <- from
    transcript.ids <- as.character(which(object@tcc@gene.map[, 2] %in% gene.names) - 1)
    transcript.names <- object@tcc@gene.map[as.numeric(transcript.ids) + 1, 1]
    if(to == "TID") {
      return(unname(transcript.names))
    }
    if(to == "TIDX") {
      return(transcript.ids)
    }
    ec <- sapply(X = as.character(transcript.ids), FUN = function(x) HashTableLookup(x, object@tcc@tid.to.ec.map))
    if(to == "EC") {
      return(as.character(unique(unlist(unname(ec)))))
    }
    gene.names <- rep(gene.names, time = length(unname(unlist(ec))))
    reps <- unname(unlist(lapply(ec, length)))
    transcript.ids <- rep(transcript.ids, times = reps)
    transcript.names <- rep(transcript.names, times = reps)
    ec <- unname(unlist(ec))
  }
  map.mat <- data.frame(EC = ec,
                        TIDX = transcript.ids,
                        TID = transcript.names,
                        GENE = gene.names,
                        stringsAsFactors = FALSE)
  return(map.mat)
}

# Mapping function from gene to EC. Option to only return unambiguously mapping
# ECs or not
#
# @param object      Seurat object
# @param gene        Gene name
# @param ambig       whether to return ambiguously mapping ECs. Default is TRUE
#
# @return returns a list of ECs that map to the gene
#
GeneToECMap <- function(object, gene, ambig = TRUE) {
  ecs <- TCCMap(object = object, from = gene, from.type = "GENE", to = "EC")
  if(ambig) {
    return(ecs)
  }
  num.mapped.genes <- sapply(ecs, FUN = function(x){
    length(TCCMap(object = object, from = x, from.type = "EC", to = "GENE"))
  })
  return(names(which(num.mapped.genes == 1)))
}

# Mapping function from gene to transcript names.Option to only return
# unambiguously mapping transcripts or not
#
# @param object      Seurat object
# @param gene        Gene name
# @param ambig       whether to return ambiguously mapping transcripts.
#                    Default is TRUE
#
# @return returns a list of ECs that map to the gene
#
GeneToTIDMap <- function(object, gene, ambig = TRUE) {
  tids <- TCCMap(object = object, from = gene, from.type = "GENE", to = "TID")
  if(ambig) {
    return(tids)
  }
  num.mapped.genes <- sapply(tids, FUN = function(x){
    length(TCCMap(object = object, from = x, from.type = "TID", to = "GENE"))
  })
  return(names(which(num.mapped.genes == 1)))
}

