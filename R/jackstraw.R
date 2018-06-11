#' Determine statistical significance of PCA scores.
#'
#' Randomly permutes a subset of data, and calculates projected PCA scores for
#' these 'random' genes. Then compares the PCA scores for the 'random' genes
#' with the observed PCA scores to determine statistical signifance. End result
#' is a p-value for each gene's association with each principal component.
#'
#' @param object Seurat object
#' @param num.pc Number of PCs to compute significance for
#' @param num.replicate Number of replicate samplings to perform
#' @param prop.freq Proportion of the data to randomly permute for each
#' replicate
#' @param display.progress Print progress bar showing the number of replicates
#' that have been processed.
#' @param do.par use parallel processing for regressing out variables faster.
#' If set to TRUE, will use half of the machines available cores (FALSE by default)
#' @param num.cores If do.par = TRUE, specify the number of cores to use.
#' Note that for higher number of cores, larger free memory is needed.
#' If \code{num.cores = 1} and \code{do.par = TRUE}, \code{num.cores} will be set to half
#' of all available cores on the machine.
#' @param maxit maximum number of iterations to be performed by the irlba function of RunPCA
#'
#' @return Returns a Seurat object where object@@dr$pca@@jackstraw@@emperical.p.value
#' represents p-values for each gene in the PCA analysis. If ProjectPCA is
#' subsequently run, object@dr$pca@jackstraw@emperical.p.value.full then
#' represents p-values for all genes.
#'
#' @import doSNOW
#' @importFrom methods new
#' @importFrom utils txtProgressBar setTxtProgressBar
#'
#' @references Inspired by Chung et al, Bioinformatics (2014)
#'
#' @export
#'
#' @examples
#' \dontrun{
#' pbmc_small = suppressWarnings(JackStraw(pbmc_small))
#' head(pbmc_small@dr$pca@jackstraw@emperical.p.value)
#' }
#'
JackStraw <- function(
  object,
  num.pc = 20,
  num.replicate = 100,
  prop.freq = 0.01,
  display.progress = TRUE,
  do.par = FALSE,
  num.cores = 1,
  maxit = 1000
) {
  if (is.null(object@dr$pca)) {
    stop("PCA has not been computed yet. Please run RunPCA().")
  }
  # error checking for number of PCs
  if (num.pc > ncol(x = GetDimReduction(object,"pca","cell.embeddings"))) {
    num.pc <- ncol(x = GetDimReduction(object,"pca","cell.embeddings"))
    warning("Number of PCs specified is greater than PCs available. Setting num.pc to ", num.pc, " and continuing.")
  }
  if (num.pc > length(x = object@cell.names)) {
    num.pc <- length(x = object@cell.names)
    warning("Number of PCs specified is greater than number of cells. Setting num.pc to ", num.pc, " and continuing.")
  }
  pc.genes <- rownames(x = GetDimReduction(object,"pca","gene.loadings"))
  if (length(x = pc.genes) < 3) {
    stop("Too few variable genes")
  }
  if (length(x = pc.genes) * prop.freq < 3) {
    warning(
      "Number of variable genes given ",
      prop.freq,
      " as the prop.freq is low. Consider including more variable genes and/or increasing prop.freq. ",
      "Continuing with 3 genes in every random sampling."
    )
  }
  md.x <- as.matrix(x = GetDimReduction(object,"pca","gene.loadings"))
  md.rot <- as.matrix(x = GetDimReduction(object,"pca","cell.embeddings"))

  rev.pca <- GetCalcParam(
    object = object,
    calculation = "RunPCA",
    parameter = "rev.pca"
  )
  weight.by.var <- GetCalcParam(
    object = object,
    calculation = "RunPCA",
    parameter = "weight.by.var"
  )
  data.use.scaled <- GetAssayData(
    object = object,
    assay.type = "RNA",
    slot = "scale.data"
  )[pc.genes,]

  # input checking for parallel options
  if (do.par) {
    if (num.cores == 1) {
      num.cores <- detectCores() / 2
      warning(paste0("do.par set to TRUE but num.cores set to 1. Setting num.cores to ", num.cores, "."))
    } else if (num.cores > detectCores()) {
      num.cores <- detectCores() - 1
      warning(paste0("num.cores set greater than number of available cores(", detectCores(), "). Setting num.cores to ", num.cores, "."))
    }
  } else if (num.cores != 1) {
      num.cores <- 1
      warning("For parallel processing, please set do.par to TRUE.")
  }

  cl <- parallel::makeCluster(num.cores)

  registerDoSNOW(cl)

  if (display.progress) {
    time_elapsed <- Sys.time()
  }

  opts <- list()
  if (display.progress) {
    # define progress bar function
    pb <- txtProgressBar(min = 0, max = num.replicate, style = 3)
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress = progress)
    time_elapsed <- Sys.time()
  }

  fake.pcVals.raw <- foreach(
    x = 1:num.replicate,
    .options.snow = opts,
    .export = c('JackRandom')
  ) %dopar% {
    JackRandom(
      scaled.data = data.use.scaled,
      prop.use = prop.freq,
      r1.use = 1,
      r2.use = num.pc,
      seed.use = x,
      rev.pca = rev.pca,
      weight.by.var = weight.by.var,
      maxit = maxit
    )
  }

  if (display.progress) {
    time_elapsed <- Sys.time() - time_elapsed
    cat(paste("\nTime Elapsed: ",time_elapsed, units(time_elapsed), "\n"))
    close(pb)
  }

  stopCluster(cl)

  fake.pcVals <- sapply(
    X = 1:num.pc,
    FUN = function(x) {
      return(as.numeric(x = unlist(x = lapply(
        X = 1:num.replicate,
        FUN = function(y) {
          return(fake.pcVals.raw[[y]][, x])
        }
      ))))
    }
  )
  jackStraw.fakePC <- as.matrix(x = fake.pcVals)
  jackStraw.empP <- as.matrix(
    sapply(
      X = 1:num.pc,
      FUN = function(x) {
        return(unlist(x = lapply(
          X = abs(md.x[, x]),
          FUN = EmpiricalP,
          nullval = abs(fake.pcVals[,x])
        )))
      }
    )
  )
  colnames(x = jackStraw.empP) <- paste0("PC", 1:ncol(x = jackStraw.empP))

  jackstraw.obj <- new(
    Class = "jackstraw.data",
    emperical.p.value  = jackStraw.empP,
    fake.pc.scores = fake.pcVals,
    emperical.p.value.full = matrix()
  )
  object <- SetDimReduction(
    object = object,
    reduction.type = "pca",
    slot = "jackstraw",
    new.data = jackstraw.obj
  )

  return(object)
}

#' Significant genes from a PCA
#'
#' Returns a set of genes, based on the JackStraw analysis, that have
#' statistically significant associations with a set of PCs.
#'
#' @param object Seurat object
#' @param pcs.use PCS to use.
#' @param pval.cut P-value cutoff
#' @param use.full Use the full list of genes (from the projected PCA). Assumes
#' that ProjectPCA has been run. Currently, must be set to FALSE.
#' @param max.per.pc Maximum number of genes to return per PC. Used to avoid genes from one PC dominating the entire analysis.
#'
#' @return A vector of genes whose p-values are statistically significant for
#' at least one of the given PCs.
#'
#' @export
#'
#' @examples
#' PCASigGenes(pbmc_small, pcs.use = 1:2)
#'
PCASigGenes <- function(
  object,
  pcs.use,
  pval.cut = 0.1,
  use.full = FALSE,
  max.per.pc = NULL
) {
  pvals.use <- GetDimReduction(object,reduction.type = "pca",slot = "jackstraw")@emperical.p.value
  pcx.use <- GetDimReduction(object,reduction.type = "pca",slot = "gene.loadings")
  if (use.full) {
    pvals.use <- GetDimReduction(object,reduction.type = "pca",slot = "jackstraw")@emperical.p.value.full
    pcx.use <- GetDimReduction(object,reduction.type = "pca",slot = "gene.loadings.full")
  }
  if (length(x = pcs.use) == 1) {
    pvals.min <- pvals.use[, pcs.use]
  }
  if (length(x = pcs.use) > 1) {
    pvals.min <- apply(X = pvals.use[, pcs.use], MARGIN = 1, FUN = min)
  }
  names(x = pvals.min) <- rownames(x = pvals.use)
  genes.use <- names(x = pvals.min)[pvals.min < pval.cut]
  if (! is.null(x = max.per.pc)) {
    pc.top.genes <- PCTopGenes(
      object = object,
      pc.use = pcs.use,
      num.genes = max.per.pc,
      use.full = use.full,
      do.balanced = FALSE
    )
    genes.use <- intersect(x = pc.top.genes, y = genes.use)
  }
  return(genes.use)
}
