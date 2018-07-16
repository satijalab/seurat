#' Determine statistical significance of PCA scores.
#'
#' Randomly permutes a subset of data, and calculates projected PCA scores for
#' these 'random' genes. Then compares the PCA scores for the 'random' genes
#' with the observed PCA scores to determine statistical signifance. End result
#' is a p-value for each gene's association with each principal component.
#'
#' @param object Seurat object
#' @param reduction.use DimReduc to use. ONLY PCA CURRENTLY SUPPORTED.
#' @param assay.use Assay used to calculate reduction.
#' @param dims Number of PCs to compute significance for
#' @param num.replicate Number of replicate samplings to perform
#' @param prop.freq Proportion of the data to randomly permute for each
#' replicate
#' @param verbose Print progress bar showing the number of replicates
#' that have been processed.
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
  reduction.use = "pca",
  assay.use = NULL,
  dims = 20,
  num.replicate = 100,
  prop.freq = 0.01,
  verbose = TRUE,
  maxit = 1000
) {
  if (reduction.use != "pca") {
    stop("Only pca for reduction.use is currently supported")
  }
  embeddings <- Embeddings(object = object[[reduction.use]])
  if (dims > ncol(embeddings)) {
    dims <- ncol(embeddings)
    warning("Number of dimensions specified is greater than those available. Setting dims to ", dims, " and continuing", immediate. = TRUE)
  }
  if (dims > ncol(object)) {
    dims <- ncol(object)
    warning("Number of dimensions specified is greater than the number of cells. Setting dims to ", dims, " and continuing", immediate. = TRUE)
  }
  reduc.features <- rownames(Loadings(object = object[[reduction.use]]))
  if (length(x = reduc.features) < 3) {
    stop("Too few features")
  }
  if (length(x = reduc.features) * prop.freq < 3) {
    warning(
      "Number of variable genes given ",
      prop.freq,
      " as the prop.freq is low. Consider including more variable genes and/or increasing prop.freq. ",
      "Continuing with 3 genes in every random sampling."
    )
  }
  assay.use <- assay.use %||% DefaultAssay(object = object)
  loadings <- Loadings(object = object[[reduction.use]])
  data.use <- GetAssayData(object = object, assay.use = assay.use, slot = "scale.data")[reduc.features, ]
  ## TODO: Proper accessors for command params
  rev.pca <- object[[paste0("RunPCA.", assay.use)]]@params$rev.pca
  weight.by.var <- object[[paste0("RunPCA.", assay.use)]]@params$weight.by.var

  ## TODO: Parallelization
  fake.vals.raw <- lapply(X = 1:num.replicate, FUN = function(x) {
    JackRandom(
      scaled.data = data.use,
      prop.use = prop.freq,
      r1.use = 1,
      r2.use = dims,
      seed.use = x,
      rev.pca = rev.pca,
      weight.by.var = weight.by.var,
      maxit = maxit
    )}
  )
  fake.vals <- sapply(
    X = 1:dims,
    FUN = function(x) {
      return(as.numeric(x = unlist(x = lapply(
        X = 1:num.replicate,
        FUN = function(y) {
          return(fake.vals.raw[[y]][, x])
        }
      ))))
    }
  )
  jackStraw.fakePC <- as.matrix(x = fake.vals)
  jackStraw.empP <- as.matrix(
    sapply(
      X = 1:dims,
      FUN = function(x) {
        return(unlist(x = lapply(
          X = abs(loadings[, x]),
          FUN = EmpiricalP,
          nullval = abs(fake.vals[,x])
        )))
      }
    )
  )
  colnames(x = jackStraw.empP) <- paste0("PC", 1:ncol(x = jackStraw.empP))
  jackstraw.obj <- new(
    Class = "jackstraw.data",
    emperical.p.value  = jackStraw.empP,
    fake.pc.scores = fake.vals,
    emperical.p.value.full = matrix()
  )
  object[[reduction.use]] <- SetDimReduc(
    object = object[[reduction.use]],
    slot = "jackstraw",
    new.data = jackstraw.obj
  )
  object <- LogSeuratCommand(object = object)
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
