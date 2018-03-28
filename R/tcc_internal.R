#' @include seurat.R
NULL

tcc <- setClass(
  Class = "tcc",
  slots = list(
    tcc.raw = "dgCMatrix",
    tcc.norm = "dgCMatrix",
    tx.raw = "dgCMatrix",
    ec.to.tx.map = "environment",
    tx.to.ec.map = "environment",
    tx.to.gene.map = "environment",
    gene.to.tx.map = "environment"
  )
)


# Function to perform mapping functions for TCC data. To specify to or from
# types as EC (equivalence class), TX (transcript name), GENE (gene name),
# MAT (matrix with full mapping info)
#
# @param from        Value to map from - should be a string.
# @param from.type   Type to map from - one of the following: EC, TX, GENE
# @param to          Type to map to - one of the following: EC, TX, GENE, MAT
#
# @return returns either a vector of mapped values or matrix with mapping info
#
TCCMap <- function(object, from, from.type, to) {
  if(length(from) > 1){
    stop("from must be of length 1")
  }
  if (!from.type %in% c("EC", "TX", "GENE")) {
    stop("Invalid from.type argument")
  }
  if (!to %in% c("EC", "TX", "GENE", "MAT")) {
    stop("Invalid to argument")
  }

  if(from.type == "EC") {
    ec <- from
    tx <- HashTableLookup(ec, object@tcc@ec.to.tx.map)
    if(to == "TX") {
      return(tx)
    }
    ec <- rep(ec, length(tx))
    from <- tx
    from.type <- "TX"
  }
  if(from.type == "TX") {
    tx <- from
    if (to == "EC" | (to == "MAT" & !exists("ec"))) {
      ec <- as.vector(sapply(X = tx, FUN = function(x) HashTableLookup(x, object@tcc@tx.to.ec.map)))
      if(to == "EC") {
        return(unique(ec))
      }
    }
    if (to == "GENE" | to == "MAT") {
      gene <- as.vector(sapply(X = tx, FUN = function(x) HashTableLookup(x, object@tcc@tx.to.gene.map)))
      if (to == "GENE") {
        return(unique(gene))
      }
    }
  }
  if(from.type == "GENE") {
    gene <- from
    tx <- HashTableLookup(gene, object@tcc@gene.to.tx.map)
    if(to == "TX") {
      return(unique(tx))
    }
    ec <- sapply(X = tx, FUN = function(x) HashTableLookup(x, object@tcc@tx.to.ec.map))
    if(to == "EC") {
      return(unique(unname(unlist(ec))))
    }
    gene <- rep(gene, time = length(unname(unlist(ec))))
    reps <- unname(unlist(lapply(ec, length)))
    tx <- rep(tx, times = reps)
    ec <- unname(unlist(ec))
  }
  map.mat <- data.frame(EC = ec,
                        TX = tx,
                        GENE = gene,
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
  tids <- TCCMap(object = object, from = gene, from.type = "GENE", to = "TX")
  if(ambig) {
    return(tids)
  }
  num.mapped.genes <- sapply(tids, FUN = function(x){
    length(TCCMap(object = object, from = x, from.type = "TX", to = "GENE"))
  })
  return(names(which(num.mapped.genes == 1)))
}

