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
  rev.pca <- slot(object = object[[paste0("RunPCA.", assay.use)]], name = "params")$rev.pca
  weight.by.var <- slot(object = object[[paste0("RunPCA.", assay.use)]], name = "params")$weight.by.var

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
    Class = "JackStrawData",
    empirical.p.values  = jackStraw.empP,
    fake.reduction.scores = fake.vals,
    empirical.p.values.full = matrix()
  )
  object[[reduction.use]] <- SetDimReduc(
    object = object[[reduction.use]],
    slot = "jackstraw",
    new.data = jackstraw.obj
  )
  object <- LogSeuratCommand(object = object)
  return(object)
}

#' @describeIn ScoreJackStraw Score JackStraw results given a JackStrawData
#' @export
#' @method ScoreJackStraw JackStrawData
#'
ScoreJackStraw.JackStrawData <- function(
  object,
  dims = 1:5,
  score.thresh = 1e-5,
  do.plot = FALSE,
  ...
) {
  pAll <- GetJS(object = object, slot = "empirical.p.values")
  pAll <- pAll[, dims, drop = FALSE]
  pAll <- as.data.frame(pAll)
  pAll$Contig <- rownames(x = pAll)
  pAll.l <- reshape2::melt(data = pAll, id.vars = "Contig")
  colnames(x = pAll.l) <- c("Contig", "PC", "Value")
  qq.df <- NULL
  score.df <- NULL
  for (i in dims) {
    q <- qqplot(x = pAll[, i], y = runif(n = 1000), plot.it = FALSE)
    pc.score <- suppressWarnings(prop.test(
      x = c(
        length(x = which(x = pAll[, i] <= score.thresh)),
        floor(x = nrow(x = pAll) * score.thresh)
      ),
      n = c(nrow(pAll), nrow(pAll))
    )$p.val)
    if (length(x = which(x = pAll[, i] <= score.thresh)) == 0) {
      pc.score <- 1
    }
    if (is.null(x = score.df)) {
      score.df <- data.frame(PC = paste0("PC", i), Score = pc.score)
    } else {
      score.df <- rbind(score.df, data.frame(PC = paste0("PC",i), Score = pc.score))
    }
    if (is.null(x = qq.df)) {
      qq.df <- data.frame(x = q$x, y = q$y, PC = paste0("PC", i))
    } else {
      qq.df <- rbind(qq.df, data.frame(x = q$x, y = q$y, PC = paste0("PC", i)))
    }
  }
  score.df$PC <- dims
  score.df <- as.matrix(score.df)
  object <- SetJS(object = object, slot = "overall.p.values", new.data = score.df)
  return(object)
}


#' @describeIn ScoreJackStraw Score JackStraw results given a DimReduc
#' @export
#' @method ScoreJackStraw DimReduc
#'
ScoreJackStraw.DimReduc <- function(
  object,
  dims = 1:5,
  score.thresh = 1e-5,
  do.plot = FALSE,
  ...
){
  jackstraw.data <- ScoreJackStraw(
    object = GetDimReduc(object = object, slot = "jackstraw"),
    dims = dims,
    score.thresh = 1e-5,
    do.plot = FALSE,
    ...
  )
  object <- SetDimReduc(object = object, slot = "jackstraw", new.data = jackstraw.data)
  return(object)
}

#' @describeIn ScoreJackStraw Score JackStraw results given a Seurat object
#' @param reduction.use Reduction associated with JackStraw to score
#' @export
#' @method ScoreJackStraw Seurat
#'
ScoreJackStraw.Seurat <- function(
  object,
  reduction.use = "pca",
  dims = 1:5,
  score.thresh = 1e-5,
  do.plot = FALSE,
  ...
){
  object[[reduction.use]] <- ScoreJackStraw(
    object = object[[reduction.use]],
    dims = dims,
    ...
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
  pvals.use <- GetDimReduction(object,reduction.type = "pca",slot = "jackstraw")@empirical.p.values
  pcx.use <- GetDimReduction(object,reduction.type = "pca",slot = "gene.loadings")
  if (use.full) {
    pvals.use <- GetDimReduction(object,reduction.type = "pca",slot = "jackstraw")@empirical.p.values.full
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


#' @describeIn GetJS Get a slot for a given JackStrawData
#' @export
#' @method GetJS JackStrawData
#'
GetJS.JackStrawData <- function(object, slot) {
  return(slot(object = object, name = slot))
}

#' @describeIn SetJS Set a slot for a given JackStrawData
#' @export
#' @method SetJS JackStrawData
#'
SetJS.JackStrawData <- function(object, slot, new.data) {
  slots.use <- c("empirical.p.values", "fake.reduction.scores",
                 "empirical.p.values.full", "overall.p.values")
  if (!slot %in% slots.use) {
    stop("'slot' must be one of ", paste(slots.use, collapse = ', '))
  }
  slot(object = object, name = slot) <- new.data
  return(object)
}

setMethod(
  f = 'show',
  signature = 'JackStrawData',
  definition = function(object) {
    empp <- GetJS(object = object, slot = "empirical.p.values")
    scored <- GetJS(object = object, slot = "overall.p.values")
    cat(
      "A JackStrawData object simulated on", nrow(empp), "features for", ncol(empp), "dimensions.\n",
      "Scored for:", nrow(scored), "dimensions."
    )
  }
)
