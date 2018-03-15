#' Add in raw TCC data to Seurat object
#'
#' @param object Seurat object
#' @param raw.counts Path to TCC matrix
#' @param ec.map  Path to file mapping between equivalence classes and
#' transcript IDs
#' @param gene.map Path to file mapping transcript IDs to gene IDs
#' @param min.ec.filter Only keep ECs with more than this many counts across all
#' cells
#' @param display.progress prints output/progress bars
#'
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @importFrom Matrix rowSums
#' @export
#'
AddTCC <- function(object, raw.counts, ec.map, gene.map, min.ec.filter = 0,
                   display.progress = TRUE){
  if (display.progress) {
    cat("Reading files\n", file = stderr())
  }
  tcc.mat <- read.table(file = raw.counts, sep = "\t", header = TRUE,
                        row.names = 1)
  tcc.mat <- as(as.matrix(tcc.mat), "dgCMatrix")
  ec.map <- as.matrix(read.table(file = ec.map, stringsAsFactors = FALSE,
                                 sep = "\t", row.names = 1))
  gene.map <- read.table(file = gene.map, stringsAsFactors = FALSE)
  ecs.to.keep <- which(Matrix::rowSums(tcc.mat) > min.ec.filter)
  tcc.mat <- tcc.mat[ecs.to.keep, ]
  ec.map <- ec.map[ecs.to.keep, ,drop = FALSE]

  if (display.progress) {
    cat("Building EC/transcript/gene maps. Note: this may take a few minutes\n", file = stderr())
    pb <- txtProgressBar(min = 0, max = nrow(ec.map), style = 3)
  }

  ht <- HashTable()
  ht2 <- HashTable()
  for(i in 1:nrow(ec.map)){
    new.tids <- as.numeric(unlist(strsplit(x = ec.map[i, ], split = ",")))
    HashTableInsert(key = rownames(ec.map)[i], value = new.tids, ht = ht)
    for(j in new.tids){
      HashTableAdd(key = as.character(j), value = as.numeric(rownames(ec.map)[i]), ht = ht2)
    }
    if(i %% 10000 == 0){
      if (display.progress) {
        setTxtProgressBar(pb, i)
      }
    }
  }
  if (display.progress) {
    setTxtProgressBar(pb, nrow(ec.map))
    close(pb)
  }
  tids.to.keep <- sort(as.numeric(ls(ht2)))
  gene.map <- gene.map[tids.to.keep, ]
  tcc <- new(
    Class = "tcc",
    tcc.raw = tcc.mat,
    ec.to.tid.map = ht,
    tid.to.ec.map = ht2,
    gene.map = gene.map
  )
  object@tcc <- tcc
  return(object)
}


