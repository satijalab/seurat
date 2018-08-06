#' Run Independent Component Analysis on gene expression
#'
#' Run fastica algorithm from the ica package for ICA dimensionality reduction.
#' For details about stored ICA calculation parameters, see
#' \code{PrintICAParams}.
#'
#' @param object Seurat object
#' @param ic.genes Genes to use as input for ICA. Default is object@@var.genes
#' @param ics.compute Number of ICs to compute
#' @param use.imputed Run ICA on imputed values (FALSE by default)
#' @param rev.ica By default, computes the dimensional reduction on the cell x
#' gene matrix. Setting to true will compute it on the transpose (gene x cell
#' matrix).
#' @param print.results Print the top genes associated with each dimension
#' @param ics.print ICs to print genes for
#' @param genes.print Number of genes to print for each IC
#' @param ica.function ICA function from ica package to run (options: icafast,
#' icaimax, icajade)
#' @param seed.use Random seed to use for fastica
#' @param reduction.name dimensional reduction name, specifies the position in the object$dr list. ica by default
#' @param reduction.key dimensional reduction key, specifies the string before the number for the dimension names. IC by default
#' @param \dots Additional arguments to be passed to fastica
#'
#' @importFrom methods new
#' @importFrom ica icafast icaimax icajade
#'
#' @return Returns Seurat object with an ICA calculation stored in
#' object@@dr$ica
#'
#' @export
#'
#' @examples
#' pbmc_small
#' # Run ICA on variable genes (default)
#' pbmc_small <- RunICA(pbmc_small, ics.compute=5)
#' # Run ICA on different gene set (in this case all genes)
#' pbmc_small <- RunICA(pbmc_small, ic.genes = rownames(pbmc_small@data))
#' # Plot results
#' ICAPlot(pbmc_small)
#'
RunICA <- function(
  object,
  ic.genes = NULL,
  ics.compute = 50,
  use.imputed = FALSE,
  rev.ica = FALSE,
  print.results = TRUE,
  ics.print = 1:5,
  genes.print = 50,
  ica.function = "icafast",
  seed.use = 1,
  reduction.name = "ica",
  reduction.key = "IC",
  ...
) {
  data.use <- PrepDR(
    object = object,
    genes.use = ic.genes,
    use.imputed = use.imputed)
  set.seed(seed = seed.use)
  ics.compute <- min(ics.compute, ncol(x = data.use))
  ica.fxn <- eval(parse(text = ica.function))
  if (rev.ica) {
    ica.results <- ica.fxn(data.use, nc = ics.compute,...)
    cell.embeddings <- ica.results$M
  } else {
    ica.results <- ica.fxn(t(x = data.use), nc = ics.compute,...)
    cell.embeddings <- ica.results$S
  }
  gene.loadings <- (as.matrix(x = data.use ) %*% as.matrix(x = cell.embeddings))
  colnames(x = gene.loadings) <- paste0(reduction.key, 1:ncol(x = gene.loadings))
  colnames(x = cell.embeddings) <- paste0(reduction.key, 1:ncol(x = cell.embeddings))
  ica.obj <- new(
    Class = "dim.reduction",
    gene.loadings = gene.loadings,
    cell.embeddings = cell.embeddings,
    sdev = sqrt(x = ica.results$vafs),
    key = "IC"
  )

  eval(expr = parse(text = paste0("object@dr$", reduction.name, "<- ica.obj")))
  parameters.to.store <- as.list(environment(), all = TRUE)[names(formals("ICA"))]
  object <- SetCalcParams(object = object, calculation = "ICA", ... = parameters.to.store)
  if(is.null(object@calc.params$ICA$ic.genes)){
    object@calc.params$ICA$ic.genes <- rownames(data.use)
  }
  if(print.results){
    PrintDim(object = object, dims.print = ics.print, genes.print = genes.print,reduction.type = reduction.name)
  }
  return(object)
}

#' Run diffusion map
#'
#' NOTE: Prior to v2.3.4, this function used the R package diffusionMap to compute
#' the diffusion map components. This package was being archived and thus
#' RunDiffusion now uses the destiny package for the diffusion computations.
#' Please be aware that this will result in different default values as the two
#' underlying package implementations are different.
#'
#' @param object Seurat object
#' @param cells.use Which cells to analyze (default, all cells)
#' @param dims.use Which dimensions to use as input features
#' @param genes.use If set, run the diffusion map procedure on this subset of
#' genes (instead of running on a set of reduced dimensions). Not set (NULL) by
#' default
#' @param reduction.use Which dimensional reduction (PCA or ICA) to use for the
#' diffusion map input. Default is PCA
#' @param q.use Quantile to clip diffusion map components at. This addresses an
#' issue where 1-2 cells will have extreme values that obscure all other points.
#' 0.01 by default
#' @param max.dim Max dimension to keep from diffusion calculation
#' @param scale.clip Max/min value for scaled data. Default is 3
#' @param reduction.name dimensional reduction name, specifies the position in
#' the object$dr list. dm by default
#' @param reduction.key dimensional reduction key, specifies the string before
#' the number for the dimension names. DM by default
#' @param ... Additional arguments to the DiffusionMap call
#'
#' @return Returns a Seurat object with a diffusion map
#'
#' @importFrom utils installed.packages
#' @importFrom stats dist quantile
#'
#' @export
#'
#' @examples
#' \dontrun{
#' pbmc_small
#' # Run Diffusion on variable genes
#' pbmc_small <- RunDiffusion(pbmc_small,genes.use = pbmc_small@var.genes)
#' # Run Diffusion map on first 10 PCs
#' pbmc_small <- RunDiffusion(pbmc_small,genes.use = pbmc_small@var.genes)
#' # Plot results
#' DMPlot(pbmc_small)
#' }
#'
RunDiffusion <- function(
  object,
  cells.use = NULL,
  dims.use = 1:5,
  genes.use = NULL,
  reduction.use = 'pca',
  q.use = 0.01,
  max.dim = 2,
  scale.clip = 10,
  reduction.name = "dm",
  reduction.key = "DM",
  ...
) {
  # Check for destiny
  if (!'destiny' %in% rownames(x = installed.packages())) {
    stop("Please install destiny - learn more at https://bioconductor.org/packages/release/bioc/html/destiny.html")
  }
  cells.use <- SetIfNull(x = cells.use, default = colnames(x = object@data))
  if (is.null(x = genes.use)) {
    dim.code <- GetDimReduction(
      object = object,
      reduction.type = reduction.use,
      slot = 'key'
    )
    dim.codes <- paste0(dim.code, dims.use)
    data.use <- FetchData(object = object, vars.all = dim.codes)
  }
  if (! is.null(x = genes.use)) {
    genes.use <- intersect(x = genes.use, y = rownames(x = object@scale.data))
    data.use <- MinMax(
      data = t(x = object@data[genes.use, cells.use]),
      min = -1 * scale.clip,
      max = scale.clip
    )
  }
  parameters.to.store <- as.list(environment(), all = TRUE)[names(formals("RunDiffusion"))]
  object <- SetCalcParams(object = object,
                          calculation = "RunDiffusion",
                          ... = parameters.to.store)
  data.dist <- dist(data.use)
  data.diffusion <- data.frame(
    destiny::DiffusionMap(data = as.matrix(data.dist),
                          n_eigs = max.dim, ...)@eigenvectors
  )
  colnames(x = data.diffusion) <- paste0(reduction.key, 1:ncol(x = data.diffusion))
  rownames(x = data.diffusion) <- cells.use
  for (i in 1:max.dim) {
    x <- data.diffusion[,i]
    x <- MinMax(
      data = x,
      min = quantile(x = x, probs = q.use),
      quantile(x = x, probs = 1-q.use)
    )
    data.diffusion[, i] <- x
  }
  object <- SetDimReduction(
    object = object,
    reduction.type = reduction.name,
    slot = "cell.embeddings",
    new.data = as.matrix(x = data.diffusion)
  )
  object <- SetDimReduction(
    object = object,
    reduction.type = reduction.name,
    slot = "key",
    new.data = "DM"
  )
  return(object)
}

#' Run PHATE
#'
#' PHATE is a data reduction method specifically designed for visualizing
#' **high** dimensional data in **low** dimensional spaces.
#' To run, you must first install the `phate` python
#' package (e.g. via pip install phate). Details on this package can be
#' found here: \url{https://github.com/KrishnaswamyLab/PHATE}. For a more in depth
#' discussion of the mathematics underlying PHATE, see the bioRxiv paper here:
#' \url{https://www.biorxiv.org/content/early/2017/12/01/120378}.
#'
#' @param object Seurat object
#' @param cells.use Which cells to analyze (default, all cells)
#' @param genes.use If set, run PHATE on this subset of genes.
#' Not set (NULL) by default
#' @param assay.type Assay to pull data for (default: 'RNA')
#' @param max.dim Total number of dimensions to embed in PHATE.
#' @param k int, optional, default: 15
#' number of nearest neighbors on which to build kernel
#' @param alpha int, optional, default: 10
#' sets decay rate of kernel tails.
#' If NA, alpha decaying kernel is not used
#' @param use.alpha boolean, default: NA
#' forces the use of alpha decaying kernel
#' If NA, alpha decaying kernel is used for small inputs
#' (n_samples < n_landmark) and not used otherwise
#' @param n.landmark int, optional, default: 2000
#' number of landmarks to use in fast PHATE
#' @param potential.method string, optional, default: 'log'
#' choose from 'log' and 'sqrt'
#' which transformation of the diffusional operator is used
#' to compute the diffusion potential
#' @param t int, optional, default: 'auto'
#' power to which the diffusion operator is powered
#' sets the level of diffusion
#' @param knn.dist.method string, optional, default: 'euclidean'.
#' The desired distance function for calculating pairwise distances on the data.
#' If 'precomputed', `data` is treated as a
#' (n_samples, n_samples) distance or affinity matrix
#' @param mds.method string, optional, default: 'metric'
#' choose from 'classic', 'metric', and 'nonmetric'
#' which MDS algorithm is used for dimensionality reduction
#' @param mds.dist.method string, optional, default: 'euclidean'
#' recommended values: 'euclidean' and 'cosine'
#' @param t.max int, optional, default: 100.
#' Maximum value of t to test for automatic t selection.
#' @param npca int, optional, default: 100
#' Number of principal components to use for calculating
#' neighborhoods. For extremely large datasets, using
#' n_pca < 20 allows neighborhoods to be calculated in
#' log(n_samples) time.
#' @param plot.optimal.t boolean, optional, default: FALSE
#' If TRUE, produce a plot showing the Von Neumann Entropy
#' curve for automatic t selection.
#' @param verbose `int` or `boolean`, optional (default : 1)
#' If `TRUE` or `> 0`, print verbose updates.
#' @param n.jobs `int`, optional (default: 1)
#' The number of jobs to use for the computation.
#' If -1 all CPUs are used. If 1 is given, no parallel computing code is
#' used at all, which is useful for debugging.
#' For n_jobs below -1, (n.cpus + 1 + n.jobs) are used. Thus for
#' n_jobs = -2, all CPUs but one are used
#' @param seed.use int or `NA`, random state (default: `NA`)
#' @param reduction.name dimensional reduction name, specifies the position in
#' the object$dr list. phate by default
#' @param reduction.key dimensional reduction key, specifies the string before
#' the number for the dimension names. PHATE by default
#' @param ... Additional arguments for `phateR::phate`
#'
#' @return Returns a Seurat object containing a PHATE representation
#'
#' @importFrom utils installed.packages
#' @export
#'
#' @references Moon K, van Dijk D, Wang Z, Burkhardt D, Chen W, van den Elzen A,
#' Hirn M, Coifman R, Ivanova N, Wolf G and Krishnaswamy S (2017).
#' "Visualizing Transitions and Structure for High Dimensional Data
#' Exploration." _bioRxiv_, pp. 120378. doi: 10.1101/120378
#' (URL: http://doi.org/10.1101/120378),
#' <URL: https://www.biorxiv.org/content/early/2017/12/01/120378>.
#' @examples
#' if (reticulate::py_module_available("phate")) {
#'
#' # Load data
#' pbmc_small
#'
#' # Run PHATE with default parameters
#' pbmc_small <- RunPHATE(object = pbmc_small)
#' # Plot results
#' DimPlot(object = pbmc_small, reduction.use = 'phate')
#'
#' # Try smaller `k` for a small dataset, and larger `t` for a noisy embedding
#' pbmc_small <- RunPHATE(object = pbmc_small, k = 4, t = 12)
#' # Plot results
#' DimPlot(object = pbmc_small, reduction.use = 'phate')
#'1
#' # For increased emphasis on local structure, use sqrt potential
#' pbmc_small <- RunPHATE(object = pbmc_small, potential.method='sqrt')
#' # Plot results
#' DimPlot(object = pbmc_small, reduction.use = 'phate')
#' }
#'
RunPHATE <- function(
  object,
  cells.use = NULL,
  genes.use = NULL,
  assay.type = 'RNA',
  max.dim = 2L,
  k = 15,
  alpha = 10,
  use.alpha = NA,
  n.landmark = 2000,
  potential.method = "log",
  t = "auto",
  knn.dist.method = "euclidean",
  mds.method = "metric",
  mds.dist.method = "euclidean",
  t.max = 100,
  npca = 100,
  plot.optimal.t = FALSE,
  verbose = 1,
  n.jobs = 1,
  seed.use = NA,
  reduction.name = "phate",
  reduction.key = "PHATE",
  ...
) {
  if (!'phateR' %in% rownames(x = installed.packages())) {
    stop("Please install phateR")
  }
  data.use <- GetAssayData(object, assay.type = assay.type, slot = "scale.data")
  if (!is.null(x = cells.use)) {
    data.use <- data.use[, cells.use]
  }
  if (!is.null(x = genes.use)) {
    data.use <- data.use[genes.use, ]
  }
  data.use <- t(x = data.use)
  parameters.to.store <- as.list(x = environment(), all = TRUE)[names(x = formals(fun = "RunPHATE"))]
  object <- SetCalcParams(
    object = object,
    calculation = "RunPHATE",
    ... = parameters.to.store
  )
  phate_output <- phateR::phate(
    data.use,
    ndim = max.dim,
    k = k,
    alpha = alpha,
    use.alpha = alpha,
    n.landmark = n.landmark,
    potential.method = potential.method,
    t = t,
    knn.dist.method = knn.dist.method,
    init = NULL,
    mds.method = mds.method,
    mds.dist.method = mds.dist.method,
    t.max = t.max,
    npca = npca,
    plot.optimal.t = plot.optimal.t,
    verbose = verbose,
    n.jobs = n.jobs,
    seed = seed.use,
    ...
  )
  phate_output <- as.matrix(x = phate_output)
  colnames(x = phate_output) <- paste0(reduction.key, 1:ncol(x = phate_output))
  object <- SetDimReduction(
    object = object,
    reduction.type = reduction.name,
    slot = "cell.embeddings",
    new.data = phate_output
  )
  object <- SetDimReduction(
    object = object,
    reduction.type = reduction.name,
    slot = "key",
    new.data = reduction.key
  )
  return(object)
}

# Check group exists either as an ident or that all cells passed as vector are
# present
#
# @param object    Seurat object
# @param group     Identity or vector of cell names
# @param group.id  Corresponds to the the either group1 or group2 parameter from RunCCA
#
CheckGroup <- function(object, group, group.id) {
  if (all(group %in% unique(x = object@ident))) {
    cells.use <- WhichCells(object = object, ident = group)
  } else {
    if (all(group %in% object@cell.names)) {
      cells.use <- group
    } else {
      stop(paste(
        group.id,
        "must be either a vector of valid cell names or idents"
      ))
    }
  }
  if (length(cells.use) == 0) {
    stop(paste0("No cells present in group: ", group.id))
  }
  return(cells.use)
}

#' Gene expression markers of identity classes defined by a phylogenetic clade
#'
#' Finds markers (differentially expressed genes) based on a branching point (node) in
#' the phylogenetic tree. Markers that define clusters in the left branch are positive markers.
#' Markers that define the right branch are negative markers.
#'
#' @inheritParams FindMarkers
#' @param node The node in the phylogenetic tree to use as a branch point
#' @param tree.use Can optionally pass the tree to be used. Default uses the tree in object@@cluster.tree
#' @param assay.type Type of assay to fetch data for (default is RNA)
#' @param ... Additional arguments passed to FindMarkers
#'
#' @return Matrix containing a ranked list of putative markers, and associated
#' statistics (p-values, ROC score, etc.)
#'
#' @export
#'
#' @examples
#' FindMarkersNode(pbmc_small, 5)
#'
FindMarkersNode <- function(
  object,
  node,
  tree.use = NULL,
  genes.use = NULL,
  logfc.threshold = 0.25,
  test.use = "wilcox",
  assay.type = "RNA",
  ...
) {
  data.use <- GetAssayData(
    object = object,
    assay.type = assay.type
  )
  genes.use <- SetIfNull(x = genes.use, default = rownames(x = data.use))
  tree <- SetIfNull(x = tree.use, default = object@cluster.tree[[1]])
  ident.order <- tree$tip.label
  nodes.1 <- ident.order[GetLeftDescendants(tree = tree, node = node)]
  nodes.2 <- ident.order[GetRightDescendants(tree = tree, node = node)]
  #print(nodes.1)
  #print(nodes.2)
  to.return <- FindMarkers(
    object = object,
    assay.type = assay.type,
    ident.1 = nodes.1,
    ident.2 = nodes.2,
    genes.use = genes.use,
    logfc.threshold = logfc.threshold,
    test.use = test.use,
    ...
  )
  return(to.return)
}

globalVariables(names = c('myAUC', 'p_val'), package = 'Seurat', add = TRUE)
#' Find all markers for a node
#'
#' This function finds markers for all splits at or below the specified node
#'
#' @param object Seurat object. Must have object@@cluster.tree slot filled. Use BuildClusterTree() if not.
#' @param node Node from which to start identifying split markers, default is top node.
#' @param genes.use Genes to test. Default is to use all genes
#' @param logfc.threshold Limit testing to genes which show, on average, at least
#' X-fold difference (log-scale) between the two groups of cells.
#' @param test.use Denotes which test to use. Seurat currently implements
#' "bimod" (likelihood-ratio test for single cell gene expression, McDavid et
#' al., Bioinformatics, 2013, default), "roc" (standard AUC classifier), "t"
#' (Students t-test), and "tobit" (Tobit-test for differential gene expression,
#' as in Trapnell et al., Nature Biotech, 2014), 'poisson', and 'negbinom'.
#' The latter two options should only be used on UMI datasets, and assume an underlying
#' poisson or negative-binomial distribution.
#' @param min.pct - only test genes that are detected in a minimum fraction of min.pct cells
#' in either of the two populations. Meant to speed up the function by not testing genes that are very infrequently expression
#' @param min.diff.pct - only test genes that show a minimum difference in the fraction of detection between the two groups. Set to -Inf by default
#' @param only.pos Only return positive markers (FALSE by default)
#' @param print.bar Print a progress bar once expression testing begins (uses pbapply to do this)
#' @param max.cells.per.ident Down sample each identity class to a max number. Default is no downsampling.
#' @param random.seed Random seed for downsampling
#' @param return.thresh Only return markers that have a p-value < return.thresh, or a power > return.thresh (if the test is ROC)
#' @param do.print Print status updates
#' @param min.cells.gene Minimum number of cells expressing the gene in at least one
#' of the two groups, currently only used for poisson and negative binomial tests
#' @param min.cells.group Minimum number of cells in one of the groups
#' @param assay.type Type of assay to fetch data for (default is RNA)
#' @param \dots Additional parameters to pass to specific DE functions
#'
#' @return Returns a dataframe with a ranked list of putative markers for each node and associated statistics
#'
#' @importFrom ape drop.tip
#'
#' @export
#'
#' @examples
#' pbmc_small
#'
#' FindAllMarkersNode(pbmc_small)
#'
FindAllMarkersNode <- function(
  object,
  node = NULL,
  genes.use = NULL,
  logfc.threshold = 0.25,
  test.use = "wilcox",
  min.pct = 0.1,
  min.diff.pct = 0.05,
  print.bar = TRUE,
  only.pos = FALSE,
  max.cells.per.ident = Inf,
  return.thresh = 1e-2,
  do.print = FALSE,
  random.seed = 1,
  min.cells.gene = 3,
  min.cells.group = 3,
  assay.type = "RNA",
  ...
) {
  if (length(object@cluster.tree) == 0) {
    stop("Tree hasn't been built yet. Run BuildClusterTree to build.")
  }
  data.use <- GetAssayData(object = object,assay.type = assay.type,slot = "data")
  genes.use <- SetIfNull(x = genes.use, default = rownames(x = data.use))
  node <- SetIfNull(x = node, default = object@cluster.tree[[1]]$edge[1, 1])
  tree.use <- object@cluster.tree[[1]]
  descendants <- DFT(tree = tree.use, node = node, path = NULL, include.children = TRUE)
  all.children <- sort(x = tree.use$edge[,2][!tree.use$edge[,2] %in% tree.use$edge[,1]])
  descendants <- MapVals(v = descendants, from = all.children, to = tree.use$tip.label)
  drop.children <- setdiff(tree.use$tip.label, descendants)
  keep.children <- setdiff(tree.use$tip.label, drop.children)
  orig.nodes <- c(node, as.numeric(setdiff(descendants, keep.children)))
  tree.use <- drop.tip(tree.use, drop.children)
  new.nodes <- unique(tree.use$edge[,1])
  if ((test.use == 'roc') && (return.thresh == 1e-2)) {
    return.thresh <- 0.7
  }
  genes.de <- list()
  for (i in ((tree.use$Nnode + 2):max(tree.use$edge))) {
    genes.de[[i]] <- FindMarkersNode(
      object = object,
      assay.type = assay.type,
      node = i,
      tree.use = tree.use,
      genes.use = genes.use,
      logfc.threshold = logfc.threshold,
      test.use = test.use,
      min.pct = min.pct,
      min.diff.pct = min.diff.pct,
      print.bar = print.bar,
      only.pos = only.pos,
      max.cells.per.ident = max.cells.per.ident,
      random.seed = random.seed,
      min.cells.gene = min.cells.gene,
      min.cells.group = min.cells.group
    )
    if (do.print) {
      message(paste("Calculating node", i))
    }
  }
  gde.all <- data.frame()
  for (i in ((tree.use$Nnode + 2):max(tree.use$edge))) {
    if (is.null(x = unlist(x = genes.de[i]))) {
      next
    }
    gde <- genes.de[[i]]
    if (nrow(x = gde) > 0) {
      if (test.use == 'roc') {
        gde <- subset(
          x = gde,
          subset = (myAUC > return.thresh | myAUC < (1 - return.thresh))
        )
      }
      if ( (test.use == 'bimod') || (test.use == 't')) {
        gde <- gde[order(gde$p_val,-gde$avg_logFC), ]
        gde <- subset(x = gde, subset = p_val < return.thresh)
      }
      if (nrow(x = gde) > 0) {
        gde$cluster <- i
        gde$gene <- rownames(x = gde)
      }
      if (nrow(x = gde) > 0) {
        gde.all <- rbind(gde.all,gde)
      }
    }
  }
  gde.all$cluster <- MapVals(
    v = gde.all$cluster,
    from = new.nodes,
    to = orig.nodes
  )
  return(gde.all)
}

#' Negative binomial test for UMI-count based data (regularized version)
#'
#' Identifies differentially expressed genes between two groups of cells using
#' a likelihood ratio test of negative binomial generalized linear models where
#' the overdispersion parameter theta is determined by pooling information
#' across genes.
#'
#' @inheritParams FindMarkers
#' @param object Seurat object
#' @param cells.1 Group 1 cells
#' @param cells.2 Group 2 cells
#' @param genes.use Genes to use for test
#' @param latent.vars Latent variables to test
#' @param print.bar Print progress bar
#' @param min.cells Minimum number of cells threshold
#' @param assay.type Type of assay to fetch data for (default is RNA)
#'
#' @return Returns a p-value ranked data frame of test results.
#'
#' @importFrom stats p.adjust
#' @importFrom utils txtProgressBar setTxtProgressBar
#'
#' @export
#'
#' @examples
#' # Note, not recommended for particularly small datasets - expect warnings
#' NegBinomDETest(
#'   object = pbmc_small,
#'   cells.1 = WhichCells(object = pbmc_small, ident = 1),
#'   cells.2 = WhichCells(object = pbmc_small, ident = 2)
#' )
#'
NegBinomRegDETest <- function(
  object,
  cells.1,
  cells.2,
  genes.use = NULL,
  latent.vars = NULL,
  print.bar = TRUE,
  min.cells = 3,
  assay.type = "RNA"
) {
  if (!is.null(genes.use)) {
    message('Make sure that genes.use contains mostly genes that are not expected to be
            differentially expressed to allow unbiased theta estimation')
  }
  genes.use <- SetIfNull(x = genes.use, default = rownames(x = GetAssayData(object = object,assay.type = assay.type,slot = "data")))
  # check that the gene made it through the any filtering that was done
  genes.use <- genes.use[genes.use %in% rownames(x = GetAssayData(object = object,assay.type = assay.type,slot = "data"))]
  message(
    sprintf(
      'NegBinomRegDETest for %d genes and %d and %d cells',
      length(x = genes.use),
      length(x = cells.1),
      length(x = cells.2)
    )
  )
  grp.fac <- factor(
    x = c(
      rep.int(x = 'A', times = length(x = cells.1)),
      rep.int(x = 'B', times = length(x = cells.2))
    )
  )
  to.test.data <- GetAssayData(object = object,assay.type = assay.type,slot = "raw.data")[genes.use, c(cells.1, cells.2), drop = FALSE]
  message('Calculating mean per gene per group')
  above.threshold <- pmax(
    apply(X = to.test.data[, cells.1] > 0, MARGIN = 1, FUN = mean),
    apply(X = to.test.data[, cells.2] > 0, MARGIN = 1, FUN = mean)
  ) >= 0.02
  message(
    sprintf(
      '%d genes are detected in at least 2%% of the cells in at least one of the groups and will be tested',
      sum(above.threshold)
    )
  )
  genes.use <- genes.use[above.threshold]
  to.test.data <- to.test.data[genes.use, , drop = FALSE]
  my.latent <- FetchData(
    object = object,
    vars.all = latent.vars,
    cells.use = c(cells.1, cells.2),
    use.raw = TRUE
  )
  to.test <- data.frame(my.latent, row.names = c(cells.1, cells.2))
  message(paste('Latent variables are', paste(latent.vars, collapse = " ")))
  # get regularized theta (ignoring group factor)
  theta.fit <- RegularizedTheta(
    cm = to.test.data,
    latent.data = to.test,
    min.theta = 0.01,
    bin.size = 128
  )
  message('Running NB regression model comparison')
  to.test$NegBinomRegDETest.group <- grp.fac
  bin.size <- 128
  bin.ind <- ceiling(1:length(x = genes.use) / bin.size)
  max.bin <- max(bin.ind)
  pb <- txtProgressBar(min = 0, max = max.bin, style = 3, file = stderr())
  res <- c()
  for (i in 1:max.bin) {
    genes.bin.use <- genes.use[bin.ind == i]
    bin.out.lst <- parallel::mclapply(
      X = genes.bin.use,
      FUN = function(j) {
        return(NBModelComparison(
          y = to.test.data[j, ],
          theta = theta.fit[j],
          latent.data = to.test,
          com.fac = latent.vars,
          grp.fac = 'NegBinomRegDETest.group'
        ))
      }
    )
    res <- rbind(res, do.call(rbind, bin.out.lst))
    setTxtProgressBar(pb = pb, value = i)
  }
  close(pb)
  rownames(res) <- genes.use
  res <- as.data.frame(x = res)
  res$adj.pval <- p.adjust(p = res$pval, method = 'fdr')
  res <- res[order(res$pval, -abs(x = res$log2.fc)), ]
  return(res)
}

globalVariables(names = 'i', package = 'Seurat', add = TRUE)
# Regress out technical effects and cell cycle
#
# Remove unwanted effects from scale.data
#
# @keywords internal
# @param object Seurat object
# @param vars.to.regress effects to regress out
# @param genes.regress gene to run regression for (default is all genes)
# @param model.use Use a linear model or generalized linear model (poisson, negative binomial) for the regression. Options are 'linear' (default), 'poisson', and 'negbinom'
# @param use.umi Regress on UMI count data. Default is FALSE for linear modeling, but automatically set to TRUE if model.use is 'negbinom' or 'poisson'
# @param verbose display progress bar for regression procedure.
# @param do.par use parallel processing for regressing out variables faster.
# If set to TRUE, will use half of the machines available cores (FALSE by default)
# @param num.cores If do.par = TRUE, specify the number of cores to use.
#
# @return Returns the residuals from the regression model
#
#' @import Matrix
#' @import doSNOW
#' @importFrom stats as.formula lm residuals glm
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @importFrom foreach foreach %dopar%
#
RegressOutResid <- function(
  object,
  vars.to.regress,
  genes.regress = NULL,
  model.use = 'linear',
  use.umi = FALSE,
  verbose = TRUE,
  do.par = FALSE,
  num.cores = 1
) {
  possible.models <- c("linear", "poisson", "negbinom")
  if (!model.use %in% possible.models) {
    stop(
      paste0(
        model.use,
        " is not a valid model. Please use one the following: ",
        paste0(possible.models, collapse = ", "),
        "."
      )
    )
  }
  genes.regress <- SetIfNull(x = genes.regress, default = rownames(x = object@data))
  genes.regress <- intersect(x = genes.regress, y = rownames(x = object@data))
  latent.data <- FetchData(object = object, vars.all = vars.to.regress)
  bin.size <- ifelse(test = model.use == 'negbinom', yes = 5, no = 100)
  bin.ind <- ceiling(x = 1:length(x = genes.regress) / bin.size)
  max.bin <- max(bin.ind)
  if (verbose) {
    message(paste("Regressing out:", paste(vars.to.regress, collapse = ", ")))
    pb <- txtProgressBar(min = 0, max = max.bin, style = 3, file = stderr())
  }
  data.resid <- c()
  data.use <- object@data[genes.regress, , drop = FALSE];
  if (model.use != "linear") {
    use.umi <- TRUE
  }
  if (use.umi) {
    data.use <- object@raw.data[genes.regress, object@cell.names, drop = FALSE]
  }
  # input checking for parallel options
  if (do.par) {
    if (num.cores == 1) {
      num.cores <- detectCores() / 2
    } else if (num.cores > detectCores()) {
      num.cores <- detectCores() - 1
      warning(paste0("num.cores set greater than number of available cores(", detectCores(), "). Setting num.cores to ", num.cores, "."))
    }
  } else if (num.cores != 1) {
    num.cores <- 1
    warning("For parallel processing, please set do.par to TRUE.")
  }
  cl <- parallel::makeCluster(num.cores)#, outfile = "")
  # using doSNOW library because it supports progress bar update
  registerDoSNOW(cl)
  opts <- list()
  if (verbose) {
    # define progress bar function
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress = progress)
    time_elapsed <- Sys.time()
  }
  reg.mat.colnames <- c(colnames(x = latent.data), "GENE")
  fmla_str = paste0("GENE ", " ~ ", paste(vars.to.regress, collapse = "+"))
  if (model.use == "linear") {
    # In this code, we'll repeatedly regress different Y against the same X
    # (latent.data) in order to calculate residuals.  Rather that repeatedly
    # call lm to do this, we'll avoid recalculating the QR decomposition for the
    # latent.data matrix each time by reusing it after calculating it once and
    # storing it in a fastResiduals function.
    regression.mat <- cbind(latent.data, data.use[1,])
    colnames(regression.mat) <- reg.mat.colnames
    qr = lm(as.formula(fmla_str), data = regression.mat, qr = TRUE)$qr
    rm(regression.mat)
  }
  data.resid <- foreach(i = 1:max.bin, .combine = "c", .options.snow = opts) %dopar% {
    genes.bin.regress <- rownames(x = data.use)[bin.ind == i]
    gene.expr <- as.matrix(x = data.use[genes.bin.regress, , drop = FALSE])
    empty_char = character(length = dim(gene.expr)[1]) # Empty vector to reuse
    new.data <- sapply(
      X = genes.bin.regress,
      FUN = function(x) {
        # Fast path for std. linear models
        if(model.use=="linear") {
          resid <- qr.resid(qr, gene.expr[x,])
        } else {
          regression.mat <- cbind(latent.data, gene.expr[x,])
          colnames(x = regression.mat) <- reg.mat.colnames
          fmla = as.formula(fmla_str)
          resid <- switch(
            EXPR = model.use,
            'poisson' = residuals(
              object = glm(
                formula = fmla,
                data = regression.mat,
                family = "poisson"
              ),
              type = 'pearson'
            ),
            'negbinom' = NBResiduals(
              fmla = fmla,
              regression.mat = regression.mat,
              gene = x,
              return.mode = TRUE
            )
          )
        }
        if (!is.list(x = resid)) {
          resid <- list('resid' = resid, 'mode' = empty_char)
        }
        return(resid)
      }
    )
    new.data.resid <- new.data[seq.int(from = 1, to = length(x = new.data), by = 2)]
    new.data.resid = matrix(unlist(new.data.resid), nrow = length(new.data.resid[[1]]))
    colnames(x = new.data.resid) <- genes.bin.regress
    new.data.mode <- unlist(x = new.data[seq.int(from = 2, to = length(x = new.data), by = 2)])
    names(x = new.data.mode) <- genes.bin.regress
    new.data <- list('resid' = new.data.resid, 'mode' = new.data.mode)
    return(new.data)
  }
  if (verbose) {
    time_elapsed <- Sys.time() - time_elapsed
    cat(paste("\nTime Elapsed: ",time_elapsed, units(time_elapsed)))
    close(pb)
  }
  stopCluster(cl)
  modes <- unlist(x = data.resid[seq.int(from = 2, to = length(x = data.resid), by = 2)])
  modes <- modes[modes == 'scale']
  names(x = modes) <- gsub(
    pattern = 'mode.',
    replacement = '',
    x = names(x = modes),
    fixed = TRUE
  )
  data.resid <- data.resid[seq.int(from = 1, to = length(x = data.resid), by = 2)]
  data.resid <- as.matrix(x = as.data.frame(x = data.resid))
  data.resid <- t(x = data.resid)
  if (length(x = modes)) {
    message(
      "The following genes failed with glm.nb, and fell back to scale(log(y+1))\n\t",
      paste(names(x = modes), collapse = ', ')
    )
  }
  rownames(x = data.resid) <- genes.regress
  suppressWarnings(expr = gc(verbose = FALSE))
  if (use.umi) {
    data.resid <- log1p(
      x = sweep(
        x = data.resid,
        MARGIN = 1,
        STATS = apply(X = data.resid, MARGIN = 1, FUN = min),
        FUN = "-"
      )
    )
  }
  return(data.resid)
}

# Regress out technical effects and cell cycle using regularized Negative Binomial regression
#
# Remove unwanted effects from umi data and set scale.data to Pearson residuals
# Uses mclapply; you can set the number of cores it will use to n with command options(mc.cores = n)
# If n.genes.step1 is set, only a (somewhat-random) subset of genes is used for estimating theta.
#
# @param object Seurat object
# @param latent.vars effects to regress out (character vector)
# @param n.genes.step1 number of genes to use when estimating theta (default uses all genes)
# @param genes.regress gene to run regression for (default uses all genes)
# @param res.clip.range numeric of length two specifying the min and max values the results will be clipped to
# @param min.theta minimum theta to use in NB regression
# @param residual.type string specifying the type of residual used (default is pearson)
# @param bin.size number of genes to put in each bin (to show progress)
# @param use.stored.theta skip the first step and use the fitted thetas from a previous run
# @param min.cells only use genes that have been observed in at least this many cells
#
# @return Returns Seurat object with the scale.data (object@scale.data) genes returning the residuals from the regression model
#
#' @import Matrix
#' @import parallel
#' @importFrom MASS theta.ml negative.binomial
#' @importFrom stats glm loess residuals approx
#' @importFrom utils txtProgressBar setTxtProgressBar
#
RegressOutNBregOld <- function(
  object,
  latent.vars,
  n.genes.step1 = NULL,
  genes.regress = NULL,
  res.clip.range = c(-30, 30),
  min.theta = 0.01,
  residual.type = 'pearson',
  bin.size = 128,
  use.stored.theta = FALSE,
  min.cells = 3
) {
  # in the first step we use all genes, except if n.genes.step1 has been set
  cm <- Matrix(object@raw.data[, colnames(x = object@data), drop = FALSE])
  gene.observed <- Matrix::rowSums(cm > 0)
  genes.regress <- SetIfNull(x = genes.regress, default = rownames(x = cm))
  genes.regress <- intersect(x = genes.regress, y = rownames(x = cm)[gene.observed >= min.cells])
  genes.step1 <- rownames(cm)[gene.observed >= min.cells]
  if (!is.null(n.genes.step1)) {
    # density-sample genes to speed up the first step
    raw.mean <- log10(rowMeans(cm[genes.step1, ]))
    raw.det.rate <- rowMeans(cm[genes.step1, ] > 0)
    dens <- apply(
      X = cbind(raw.mean, raw.det.rate),
      MARGIN = 2,
      FUN = function(y) {
        y.dens <- density(x = y, bw = 'nrd', adjust = 1)
        ret <- approx(x = y.dens$x, y = y.dens$y, xout = y)$y
        return(ret / sum(ret))
      }
    )
    sampling.prob <- 1 / apply(X = dens, MARGIN = 1, FUN = min)
    genes.step1 <- sample(x = genes.step1, size = n.genes.step1, prob = sampling.prob)
  }
  latent.data <- FetchData(object = object, vars.all = latent.vars)
  bin.ind <- ceiling(x = 1:length(x = genes.step1) / bin.size)
  max.bin <- max(bin.ind)
  message(paste("Regressing out", paste(latent.vars, collapse = ' ')))
  if (use.stored.theta) {
    message('Using previously fitted theta values for NB regression')
    theta.fit <- object@misc[['NBreg.theta.fit']]
  } else {
    message('First step: Poisson regression (to get initial mean), and estimate theta per gene')
    message('Using ', length(x = genes.step1), ' genes')
    pb <- txtProgressBar(min = 0, max = max.bin, style = 3)
    theta.estimate <- c()
    for (i in 1:max.bin) {
      genes.bin.regress <- genes.step1[bin.ind == i]
      bin.theta.estimate <- unlist(
        parallel::mclapply(
          X = genes.bin.regress,
          FUN = function(j) {
            as.numeric(
              x = MASS::theta.ml(
                as.numeric(x = unlist(x = cm[j, ])),
                glm(as.numeric(x = unlist(x = cm[j, ])) ~ ., data = latent.data, family = poisson)$fitted
              )
            )
          }
        ),
        use.names = FALSE
      )
      theta.estimate <- c(theta.estimate, bin.theta.estimate)
      setTxtProgressBar(pb, i)
    }
    close(pb)
    names(theta.estimate) <- genes.step1
    UMI.mean <- rowMeans(x = cm[genes.step1, ])
    var.estimate <- UMI.mean + (UMI.mean ^ 2) / theta.estimate
    fit <- loess(log10(var.estimate) ~ log10(UMI.mean), span = 0.33)
    genes.regress.mean <- rowMeans(x = cm[genes.regress, ])
    var.fit.log10 <- predict(fit, log10(genes.regress.mean))
    theta.fit <- (genes.regress.mean ^ 2) / (10 ^ var.fit.log10 - genes.regress.mean)
    to.fix <- theta.fit <= min.theta | is.infinite(x = theta.fit)
    if (any(to.fix)) {
      message(
        'Fitted theta below ',
        min.theta,
        ' for ',
        sum(to.fix),
        ' genes, setting them to ',
        min.theta
      )
      theta.fit[to.fix] <- min.theta
    }
    # save theta estimate and fitted theta in object
    object@misc <- as.list(object@misc)
    object@misc[['NBreg.theta.estimate']] <- theta.estimate
    object@misc[['NBreg.theta.fit']] <- theta.fit
  }
  message('Second step: NB regression with fixed theta for ', length(x = genes.regress), ' genes')
  bin.ind <- ceiling(x = 1:length(x = genes.regress) / bin.size)
  max.bin <- max(bin.ind)
  pb <- txtProgressBar(min = 0, max = max.bin, style = 3, file = stderr())
  res <- c()
  for (i in 1:max.bin) {
    genes.bin.regress <- genes.regress[bin.ind == i]
    names(genes.bin.regress) <- genes.bin.regress
    bin.res.lst <- parallel::mclapply(
      X = genes.bin.regress,
      FUN = function(j) {
        fit <- 0
        try(
          fit <- glm(
            as.numeric(x = unlist(x = cm[j, ])) ~ .,
            data = latent.data,
            family = MASS::negative.binomial(theta = theta.fit[j])
          ),
          silent = TRUE
        )
        if (class(fit)[1] == 'numeric') {
          message <-
            sprintf(
              'glm and family=negative.binomial(theta=%f) failed for gene %s; falling back to scale(log10(y+1))',
              theta.fit[j],
              j
            )
          res <- scale(x = log10(as.numeric(x = unlist(x = cm[j, ])) + 1))[, 1]
        } else {
          message <- NULL
          res <- residuals(object = fit, type = residual.type)
        }
        return(list(res = res, message = message))
      }
    )
    # Print message to keep track of the genes for which glm failed to converge
    message <- unlist(x = lapply(X = bin.res.lst, FUN = function(l) { return(l$message) }), use.names = FALSE)
    if (!is.null(x = message)) {
      message(paste(message, collapse = "\n"))
    }
    bin.res.lst <- lapply(X = bin.res.lst, FUN = function(l) { return(l$res) })
    res <- rbind(res, do.call(rbind, bin.res.lst))
    setTxtProgressBar(pb, i)
  }
  close(pb)
  dimnames(x = res) <- list(genes.regress, colnames(x = cm))
  res[res < res.clip.range[1]] <- res.clip.range[1]
  res[res > res.clip.range[2]] <- res.clip.range[2]
  object@scale.data <- res
  return(object)
}
