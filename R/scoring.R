#' Calculate module scores for gene expression programs in single cells
#'
#' Calculate the average expression levels of each program (cluster) on single cell level,
#' subtracted by the aggregated expression of control gene sets.
#' All analyzed genes are binned based on averaged expression, and the control genes are
#' randomly selected from each bin.
#'
#' @param object Seurat object
#' @param genes.list Gene expression programs in list
#' @param genes.pool List of genes to check expression levels agains, defaults to rownames(x = object@data)
#' @param n.bin Number of bins of aggregate expression levels for all analyzed genes
#' @param seed.use Random seed for sampling
#' @param ctrl.size Number of control genes selected from the same bin per analyzed gene
#' @param use.k Use gene clusters returned from DoKMeans()
#' @param enrich.name Name for the expression programs
#' @param random.seed Set a random seed
#'
#' @return Returns a Seurat object with module scores added to object@meta.data
#'
#' @importFrom Hmisc cut2
#' @importFrom Matrix rowMeans colMeans
#'
#' @references Tirosh et al, Science (2016)
#'
#' @export
#'
#' @examples
#' cd_genes <- list(c(
#'   'CD79B',
#'   'CD79A',
#'   'CD19',
#'   'CD180',
#'   'CD200',
#'   'CD3D',
#'   'CD2',
#'   'CD3E',
#'   'CD7',
#'   'CD8A',
#'   'CD14',
#'   'CD1C',
#'   'CD68',
#'   'CD9',
#'   'CD247'
#' ))
#' pbmc_small <- AddModuleScore(
#'   object = pbmc_small,
#'   genes.list = cd_genes,
#'   ctrl.size = 5,
#'   enrich.name = 'CD_Genes'
#' )
#' head(x = pbmc_small@meta.data)
#'
AddModuleScore <- function(
  object,
  genes.list = NULL,
  genes.pool = NULL,
  n.bin = 25,
  seed.use = 1,
  ctrl.size = 100,
  use.k = FALSE,
  enrich.name = "Cluster",
  random.seed = 1
) {
  set.seed(seed = random.seed)
  genes.old <- genes.list
  if (use.k) {
    genes.list <- list()
    for (i in as.numeric(x = names(x = table(object@kmeans.obj[[1]]$cluster)))) {
      genes.list[[i]] <- names(x = which(x = object@kmeans.obj[[1]]$cluster == i))
    }
    cluster.length <- length(x = genes.list)
  } else {
    if (is.null(x = genes.list)) {
      stop("Missing input gene list")
    }
    genes.list <- lapply(
      X = genes.list,
      FUN = function(x) {
        return(intersect(x = x, y = rownames(x = object@data)))
      }
    )
    cluster.length <- length(x = genes.list)
  }
  if (!all(LengthCheck(values = genes.list))) {
    warning(paste(
      'Could not find enough genes in the object from the following gene lists:',
      paste(names(x = which(x = ! LengthCheck(values = genes.list)))),
      'Attempting to match case...'
    ))
    genes.list <- lapply(
      X = genes.old,
      FUN = CaseMatch, match = rownames(x = object@data)
    )
  }
  if (!all(LengthCheck(values = genes.list))) {
    stop(paste(
      'The following gene lists do not have enough genes present in the object:',
      paste(names(x = which(x = ! LengthCheck(values = genes.list)))),
      'exiting...'
    ))
  }
  if (is.null(x = genes.pool)) {
    genes.pool = rownames(x = object@data)
  }
  data.avg <- Matrix::rowMeans(x = object@data[genes.pool, ])
  data.avg <- data.avg[order(data.avg)]
  data.cut <- as.numeric(x = Hmisc::cut2(
    x = data.avg,
    m = round(x = length(x = data.avg) / n.bin)
  ))
  names(x = data.cut) <- names(x = data.avg)
  ctrl.use <- vector(mode = "list", length = cluster.length)
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
  ctrl.scores <- matrix(
    data = numeric(length = 1L),
    nrow = length(x = ctrl.use),
    ncol = ncol(x = object@data)
  )
  for (i in 1:length(ctrl.use)) {
    genes.use <- ctrl.use[[i]]
    ctrl.scores[i, ] <- Matrix::colMeans(x = object@data[genes.use, ])
  }
  genes.scores <- matrix(
    data = numeric(length = 1L),
    nrow = cluster.length,
    ncol = ncol(x = object@data)
  )
  for (i in 1:cluster.length) {
    genes.use <- genes.list[[i]]
    data.use <- object@data[genes.use, , drop = FALSE]
    genes.scores[i, ] <- Matrix::colMeans(x = data.use)
  }
  genes.scores.use <- genes.scores - ctrl.scores
  rownames(x = genes.scores.use) <- paste0(enrich.name, 1:cluster.length)
  genes.scores.use <- as.data.frame(x = t(x = genes.scores.use))
  rownames(x = genes.scores.use) <- colnames(x = object@data)
  object <- AddMetaData(
    object = object,
    metadata = genes.scores.use,
    col.name = colnames(x = genes.scores.use)
  )
  gc(verbose = FALSE)
  return(object)
}

#' Score cell cycle phases
#'
#' @param object A Seurat object
#' @param g2m.genes A vector of genes associated with G2M phase
#' @param s.genes A vector of genes associated with S phases
#' @param set.ident If true, sets identity to phase assignments
#' Stashes old identities in 'old.ident'
#'
#' @return A Seurat object with the following columns added to object@meta.data: S.Score, G2M.Score, and Phase
#'
#' @seealso \code{AddModuleScore}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # pbmc_small doesn't have any cell-cycle genes
#' # To run CellCycleScoring, please use a dataset with cell-cycle genes
#' # An example is available at http://satijalab.org/seurat/cell_cycle_vignette.html
#' pbmc_small <- CellCycleScoring(
#'   object = pbmc_small,
#'   g2m.genes = cc.genes$g2m.genes,
#'   s.genes = cc.genes$s.genes
#' )
#' head(x = pbmc_small@meta.data)
#' }
#'
CellCycleScoring <- function(
  object,
  g2m.genes,
  s.genes,
  set.ident = FALSE
) {
  enrich.name <- 'Cell Cycle'
  genes.list <- list('S.Score' = s.genes, 'G2M.Score' = g2m.genes)
  object.cc <- AddModuleScore(
    object = object,
    genes.list = genes.list,
    enrich.name = enrich.name,
    ctrl.size = min(vapply(X = genes.list, FUN = length, FUN.VALUE = numeric(1)))
  )
  cc.columns <- grep(pattern = enrich.name, x = colnames(x = object.cc@meta.data))
  cc.scores <- object.cc@meta.data[, cc.columns]
  rm(object.cc)
  gc(verbose = FALSE)
  assignments <- apply(
    X = cc.scores,
    MARGIN = 1,
    FUN = function(scores, first = 'S', second = 'G2M', null = 'G1') {
      if (all(scores < 0)) {
        return(null)
      } else {
        return(c(first, second)[which(x = scores == max(scores))])
      }
    }
  )
  cc.scores <- merge(x = cc.scores, y = data.frame(assignments), by = 0)
  colnames(x = cc.scores) <- c('rownames', 'S.Score', 'G2M.Score', 'Phase')
  rownames(x = cc.scores) <- cc.scores$rownames
  cc.scores <- cc.scores[, c('S.Score', 'G2M.Score', 'Phase')]
  object <- AddMetaData(object = object, metadata = cc.scores)
  if (set.ident) {
    object <- StashIdent(object = object, save.name = 'old.ident')
    object <- SetAllIdent(object = object, id = 'Phase')
  }
  return(object)
}
