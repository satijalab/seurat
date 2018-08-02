#' Add in raw TCC data to Seurat object
#'
#' @param object Seurat object
#' @param tcc.counts Path to TCC matrix
#' @param tx.counts Path to transcript counts file
#' @param ec.map  Path to file mapping between equivalence classes and
#' transcript IDs
#' @param gene.map Path to file mapping transcript IDs to gene IDs
#' @param min.ec.filter Only keep ECs with at least this many counts across all
#' cells
#' @param display.progress prints output/progress bars
#'
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @importFrom Matrix rowSums
#' @export
#'
AddTCC <- function(object, tcc.counts, tx.counts, ec.map, gene.map,
                   min.ec.filter = 0, display.progress = TRUE){
  if (display.progress) {
    cat("Reading files\n", file = stderr())
  }
  tcc.mat <- read.table(file = tcc.counts)
  tcc.mat <- sparseMatrix(i = tcc.mat$V2 + 1, j = tcc.mat$V1 +1, x = tcc.mat$V3)
  tcc.mat <- as(tcc.mat, "dgCMatrix")
  tx.mat <- readRDS(file = tx.counts)
  tx.mat <- as(as.matrix(tx.mat$counts), "dgCMatrix")
  ec.map <- as.matrix(read.table(file = ec.map, stringsAsFactors = FALSE,
                                 sep = "\t", row.names = 1))
  gene.map <- read.table(file = gene.map, stringsAsFactors = FALSE)
  gene.map <- as.matrix(gene.map)
  rownames(gene.map) <- gene.map[, 1]
  ecs.to.keep <- which(Matrix::rowSums(tcc.mat) >= min.ec.filter)
  tcc.mat <- tcc.mat[ecs.to.keep, ]
  ec.map <- ec.map[ecs.to.keep, ,drop = FALSE]

  if (display.progress) {
    cat("Building EC/transcript maps. Note: this may take a few minutes\n", file = stderr())
    pb <- txtProgressBar(min = 0, max = nrow(ec.map), style = 3)
  }

  ec.to.tx.ht <- HashTable()
  tx.to.ec.ht <- HashTable()

  for(i in 1:nrow(ec.map)){
    new.tids <- as.numeric(unlist(strsplit(x = ec.map[i, ], split = ",")))
    new.tx <- unname(gene.map[new.tids + 1, 1])
    HashTableInsert(key = rownames(ec.map)[i], value = new.tx, ht = ec.to.tx.ht)
    for(j in new.tx){
      HashTableAdd(key = j, value = rownames(ec.map)[i], ht = tx.to.ec.ht)
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
  tx.to.keep <- ls(tx.to.ec.ht)
  gene.map <- gene.map[tx.to.keep, ]

  if (display.progress) {
    cat("Building gene/transcript maps. Note: this may take a few minutes\n", file = stderr())
    pb <- txtProgressBar(min = 0, max = nrow(gene.map), style = 3)
  }

  tx.to.gene.ht <- HashTable()
  gene.to.tx.ht <- HashTable()
  for(i in 1:nrow(gene.map)) {
    HashTableInsert(key = gene.map[i, 1], value = gene.map[i, 2], ht = tx.to.gene.ht)
    HashTableAdd(key = gene.map[i, 2], value = gene.map[i, 1], ht = gene.to.tx.ht)
    if(i %% 1000 == 0){
      if (display.progress) {
        setTxtProgressBar(pb, i)
      }
    }
  }
  if (display.progress) {
    close(pb)
  }

  tcc <- new(
    Class = "tcc",
    tcc.raw = tcc.mat,
    tx.raw = tx.mat,
    ec.to.tx.map = ec.to.tx.ht,
    tx.to.ec.map = tx.to.ec.ht,
    tx.to.gene.map = tx.to.gene.ht,
    gene.to.tx.map = gene.to.tx.ht
  )
  object@tcc <- tcc
  return(object)
}
