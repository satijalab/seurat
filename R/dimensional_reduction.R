#' @include seurat.R
NULL

# Set up dim.reduction class

dim.reduction <- setClass("dim.reduction", slots = list(
  rotation = "matrix", x = "matrix", x.full = "matrix", sdev = "numeric", key = "character", 
  misc = "ANY"
))

#' Dimensional Reduction Accessor Function
#' 
#' Pull information for specified stored dimensional reduction analysis
#' 
#' @param object Seurat object
#' @param reduction.type Type of dimensional reduction to fetch (default is PCA)
#' @param slot Specific information to pull (i.e. rotation, x, x.full, ...)
#' @return Returns results from reduction technique
#' @export
GetDimReduction <- function(object, reduction.type = "pca", slot = "x") {
  if (!(reduction.type %in% names(object@dr))) {
    stop(paste(reduction.type, " dimensional reduction has not been computed"))
  }
  if (!(slot %in% slotNames(eval(parse(text = paste0("object@dr$", reduction.type)))))) {
    stop(paste0(slot, " slot doesn't exist"))
  }
  return(eval(parse(text = paste0("object@dr$", reduction.type, "@", slot))))
}

#' Dimensional Reduction Rotation Accessor Function
#' 
#' Pull rotation matrix for specified stored dimensional reduction analysis
#' 
#' @param object Seurat object
#' @param reduction.type Type of dimensional reduction to fetch (default is PCA)
#' @param dims.use Dimensions to include (default is all stored dims)
#' @param cells.use Cells to include (default is all cells)
#' @return Rotation matrix for given reduction, cells, and dimensions
#' @export
DimRot <- function(object, reduction.type = "pca", dims.use = NULL, cells.use = NULL) {
  object.rot <- GetDimReduction(object, reduction.type = reduction.type, slot = "rotation")
  if(length(object.rot) == 0){
    stop(paste0("Rotation slot for ", reduction.type, " is empty."))
  }
  cells.use <- set.ifnull(cells.use, rownames(object.rot))
  if (any(!cells.use %in% rownames(object.rot))){
    missing.cells <- paste0(cells.use[which(!cells.use %in% rownames(object.rot))], collapse = ", ")
    warning(paste0("Could not find the following cell names: ", missing.cells))
    cells.use <- ainb(cells.use, rownames(object.rot))
  }
  dims.use <- set.ifnull(dims.use, 1:ncol(object.rot))
  if(any(!dims.use %in% 1:ncol(object.rot))){
    missing.dims <- paste0(dims.use[which(!dims.use %in% 1:ncol(object.rot))], collapse = ", ")
    stop(paste0("Could not find the following dimensions: ", missing.dims))
  }
  object.rot <- object.rot[cells.use, dims.use, drop = F]
  object.key <- GetDimReduction(object, reduction.type = reduction.type, slot = "key")
  if(length(object.key) == 0){
    colnames(object.rot) <- NULL
  }
  else{
    colnames(object.rot) <- paste0(object.key, dims.use)
  }
  return(object.rot)
}

#' Dimensional Reduction Gene Loadings Accessor Function
#' 
#' Pull "x" matrix for specified stored dimensional reduction analysis. 
#' 
#' @param object Seurat object
#' @param reduction.type Type of dimensional reduction to fetch (default is PCA)
#' @param dims.use Dimensions to include (default is all stored dims)
#' @param genes.use Genes to include (default is all genes)
#' @param use.full Return projected gene loadings (default is FALSE)
#' @return X matrix for given reduction, cells, and genes
#' @export
DimX <- function(object, reduction.type = "pca", dims.use = NULL, genes.use = NULL, use.full=FALSE) {
  if (use.full){
    object.x <- GetDimReduction(object, reduction.type = reduction.type, slot = "x.full")
  }
  else{
    object.x <- GetDimReduction(object, reduction.type = reduction.type, slot = "x")
  }
  if(length(object.x) == 0){
    stop(paste0("X slot for ", reduction.type, " is empty."))
  }
  genes.use <- set.ifnull(genes.use, rownames(object.x))
  if (any(!genes.use %in% rownames(object.x))){
    missing.genes <- paste0(genes.use[which(!genes.use %in% rownames(object.x))], collapse = ", ")
    warning(paste0("Could not find the following gene names: ", missing.genes))
    genes.use <- ainb(genes.use, rownames(object.x))
  }
  dims.use <- set.ifnull(dims.use, 1:ncol(object.x))
  if(any(!dims.use %in% 1:ncol(object.x))){
    missing.dims <- paste0(dims.use[which(!dims.use %in% 1:ncol(object.x))], collapse = ", ")
    stop(paste0("Could not find the following dimensions: ", missing.dims))
  }
  object.x <- object.x[genes.use, dims.use, drop = F]
  object.key <- GetDimReduction(object, reduction.type = reduction.type, slot = "key")
  if(length(object.key) == 0){
    colnames(object.x) <- NULL
  }
  else{
    colnames(object.x) <- paste0(object.key, dims.use)
  }
  return(object.x)
}


#' Diffusion Maps Rotation Accessor Function
#' 
#' Pull Diffusion maps rotation matrix
#' 
#' @param object Seurat object
#' @param dims.use Dimensions to include (default is all stored dims)
#' @param cells.use Cells to include (default is all cells)
#' @return Diffusion maps rotation matrix for given cells and DMs
#' @export
DMRot <- function(object, dims.use = NULL, cells.use = NULL) {
  return(DimRot(object = object, reduction.type = "dm", dims.use = dims.use, cells.use = cells.use))
}

#' PCA Rotation Accessor Function
#' 
#' Pull PCA rotation matrix
#' 
#' @param object Seurat object
#' @param dims.use, Dimensions to include (default is all stored dims)
#' @param cells.use, Cells to include (default is all cells)
#' @return PCA rotation matrix for given cells and PCs
#' @export
PCARot <- function(object, dims.use = NULL, cells.use = NULL) {
  return(DimRot(object = object, reduction.type = "pca", dims.use = dims.use, cells.use = cells.use))
}

#' ICA Rotation Accessor Function
#' 
#' Pull ICA rotation matrix
#' 
#' @param object Seurat object
#' @param dims.use, Dimensions to include (default is all stored dims)
#' @param cells.use, Cells to include (default is all cells)
#' @return ICA rotation matrix for given cells and ICs
#' @export
ICARot <- function(object, dims.use = NULL, cells.use = NULL) {
  return(DimRot(object = object, reduction.type = "ica", dims.use = dims.use, cells.use = cells.use))
}

#' PCA Gene Loadings Accessor Function
#' 
#' Pull PCA Gene Loadings
#' 
#' @param object Seurat object
#' @param dims.use, Dimensions to include (default is all stored dims)
#' @param genes.use, Genes to include (default is all genes)
#' @return PCA loading matrix for given genes and PCs
#' @export
PCAX <- function(object, dims.use = NULL, genes.use = NULL, use.full = FALSE) {
  return(DimX(object = object, reduction.type = "pca", dims.use = dims.use, genes.use = genes.use,
              use.full = use.full))
}

#' ICA Gene Loadings Accessor Function
#' 
#' Pull ICA Gene Loadings
#' 
#' @param object Seurat object
#' @param dims.use, Dimensions to include (default is all stored dims)
#' @param genes.use, Genes to include (default is all)
#' @return ICA loading matrix for given genes and ICs
#' @export
ICAX <- function(object, dims.use = NULL, genes.use = NULL, use.full = F) {
  return(DimX(object = object, reduction.type = "ica", dims.use = dims.use, genes.use = genes.use,
              use.full = use.full))
}

#' Dimensional Reduction Mutator Function
#' 
#' Set information for specified stored dimensional reduction analysis
#' 
#' @param object Seurat object
#' @param reduction.type Type of dimensional reduction to set
#' @param slot Specific information to set (i.e. rotation, x, x.full, ...)
#' @param new.data New data to insert
#' @return Seurat object with updated slot
#' @export
SetDimReduction <- function(object, reduction.type, slot, new.data) {
  if (reduction.type %in% names(object@dr)) {
    eval(parse(text = paste0("object@dr$", reduction.type, "@", slot, "<- new.data")))
  }
  else{
    new.dr <- new("dim.reduction")
    eval(parse(text = paste0("new.dr@", slot, "<- new.data")))
    eval(parse(text = paste0("object@dr$", reduction.type, "<- new.dr")))
  }
  return(object)
}


#' Dimensional Reduction
#' 
#' Various methods for dimensional reductions
#' 
#' @param object Seurat object
#' @param reduction.type Type of dimensional reduction to run. Options include "pca", "pcafast", "ica" 
#' @param genes.use Genes to use as input for the dimensional reduction technique. Default is 
#' object@@var.genes.
#' @param dims.store Number of dimensions to store
#' @param dims.compute Number of dimensions to compute (for fast approximations)
#' @param use.imputed Whether to run the dimensional reduction on imputed values.
#' @param rev.reduction By default, computes the dimensional reduction on the cell x gene matrix.
#' Setting to true will compute it on the transpose (gene x cell matrix).
#' @param print.results Print the top genes associated with each dimension
#' @param dims.print Number of dimensions to print genes for
#' @param genes.print Number of genes to print for each dimension
#' @param ica.fxn Function to use if calculating ica
#' @param seed.use Random seed
#' @param ... Additional arguments to be passed to specific reduction technique
#' @return Returns a Seurat object with the dimensional reduction information stored 
#' @importFrom ica icafast icaimax icajade
#' @export
DimReduction <- function(object, reduction.type = NULL, genes.use = NULL, dims.store = 40, 
                         dims.compute = 40, use.imputed = FALSE, rev.reduction = FALSE, 
                         print.results = TRUE, dims.print = 1:5, genes.print = 30, ica.fxn = icafast,
                         seed.use = 1, ...){
  
  if (length(object@scale.data) == 0){
    stop("Object@scale.data has not been set. Run ScaleData() and then retry.")
  }
  if (length(object@var.genes) == 0 && is.null(genes.use)) {
    stop("Variable genes haven't been set. Run MeanVarPlot() or provide a vector of genes names in 
         genes.use and retry.")
  }
  
  dims.store <- min(dims.store, dims.compute)
  if (use.imputed) {
    data.use <- t(scale(t(object@imputed)))
  }
  else{
    data.use <- object@scale.data
  }
  
  genes.use <- set.ifnull(genes.use, object@var.genes)
  genes.use <- unique(genes.use[genes.use %in% rownames(data.use)])
  genes.var <- apply(data.use[genes.use, ], 1, var)
  genes.use <- genes.use[genes.var > 0]
  genes.use <- genes.use[!is.na(genes.use)]
  
  data.use <- data.use[genes.use, ]
  
  # call reduction technique
  reduction.type <- tolower(reduction.type)
  if(reduction.type == "pca"){
    pcaobj <- RunPCA(data.use = data.use, rev.pca = rev.reduction, pcs.store = dims.store, ...)
    object@dr$pca <- pcaobj
  }
  if(reduction.type == "pcafast") {
    pcafastobj <- RunPCAFast(data.use = data.use, rev.pca = rev.reduction, pcs.store = dims.store, 
                             pcs.compute = dims.compute, ...)
    object@dr$pca <- pcafastobj
    reduction.type <- "pca"
  }
  if(reduction.type == "ica") {
    icaobj=RunICA(data.use = data.use, ics.compute = dims.store, rev.ica = rev.reduction, 
                  ics.store = dims.store, ica.fxn = ica.fxn, seed.use = seed.use, ...)
    object@dr$ica=icaobj
  }
  
  if(print.results){
    PrintDim(object, reduction.type = reduction.type, dims.print = dims.print,
             genes.print = genes.print)
  }
  return(object)
} 

RunPCA <- function(data.use, rev.pca, pcs.store, ...){
  pca.results <- NULL
  if(rev.pca){
    pcs.store <- min(pcs.store, ncol(data.use))
    pca.results <- prcomp(data.use, ...)
    x <- pca.results$x[, 1:pcs.store]
    rotation <- pca.results$rotation[, 1:pcs.store]
  }
  else{
    pcs.store <- min(pcs.store, nrow(data.use))
    pca.results <- prcomp(t(data.use), ...)
    x <- pca.results$rotation[, 1:pcs.store]
    rotation <- pca.results$x[, 1:pcs.store]
  }
  pca.obj <- new("dim.reduction", x = x, rotation = rotation, sdev = pca.results$sdev, key = "PC")
  return(pca.obj)  
}


RunICA <- function(data.use, ics.compute, rev.ica, ica.fxn = icafast, ics.store, seed.use, ...) {
  set.seed(seed.use)
  ics.store <- min(ics.store, ncol(data.use))
  ica.results <- NULL
  if(rev.ica){
    ica.results <- ica.fxn(data.use, nc = ics.compute,...)
    rotation <- ica.results$M[, 1:ics.store]
  }
  else{
    ica.results <- ica.fxn(t(data.use), nc = ics.compute,...)
    rotation <- ica.results$S[, 1:ics.store]
  }
  
  x <- (as.matrix(data.use )%*% as.matrix(rotation))
  colnames(x) <- paste0("IC", 1:ncol(x))
  colnames(rotation) <- paste0("IC", 1:ncol(x))
  
  ica.obj <- new("dim.reduction", x = x, rotation = rotation, sdev = sqrt(ica.results$vafs), key = "IC")
  return(ica.obj)  
}

RunPCAFast <- function(data.use, rev.pca, pcs.store, pcs.compute, ...){
  pcs.compute <- min(pcs.compute, ncol(data.use))
  pcs.store <- min(pcs.store, pcs.compute)
  pca.results <- NULL
  if(rev.pca){
    pca.results <- irlba(data.use, nv = pcs.compute, ...)
    x <- pca.results$u[, 1:pcs.store]
    rotation <- pca.results$v[, 1:pcs.store]
  }
  else{
    pca.results <- irlba(t(data.use), nv = pcs.compute, ...)
    x <- pca.results$v[, 1:pcs.store]
    rotation <- pca.results$u[, 1:pcs.store]
    
  }
  rownames(x) <- rownames(data.use)
  colnames(x) <- paste0("PC", 1:pcs.compute)
  rownames(rotation) <- colnames(data.use)
  colnames(rotation) <- colnames(x)
  pca.obj <- new("dim.reduction", x = x, rotation = rotation, sdev = pca.results$d, key = "PC")
  return(pca.obj)
}

#' Run Principal Component Analysis on gene expression
#'
#' Run prcomp for PCA dimensionality reduction
#'
#' @param object Seurat object
#' @param pc.genes Genes to use as input for PCA. Default is object@@var.genes
#' @param do.print Print the top genes associated with high/low loadings for
#' the PCs
#' @param pcs.print Number of PCs to print genes for
#' @param pcs.store Number of PCs to store
#' @param genes.print Number of genes to print for each PC
#' @param use.imputed Run PCA on imputed values (FALSE by default)
#' @param rev.pca By default computes the PCA on the cell x gene matrix. Setting to true will 
#' compute it on gene x cell matrix. 
#' @param \dots Additional arguments to be passed to prcomp
#' @return Returns Seurat object with an PCA embedding (object@@pca.rot) and
#' gene projection matrix (object@@pca.x). The PCA object itself is stored in
#' object@@pca.obj[[1]]
#' @export
PCA <- function(object, pc.genes = NULL, do.print = TRUE, pcs.print = 5, pcs.store = 40, 
                genes.print = 30, use.imputed = FALSE, rev.pca = FALSE, ...) {
  return(DimReduction(object, reduction.type = "pca", genes.use = pc.genes, print.results = do.print,
                      dims.store = pcs.store, genes.print = genes.print, use.imputed = use.imputed, 
                      rev.reduction = rev.pca, ...))
}


#' Run Principal Component Analysis on gene expression using IRLBA
#'
#' Run Fast PCA dimensionality reduction
#'
#' @param object Seurat object
#' @param pc.genes Genes to use as input for PCA. Default is object@@var.genes
#' @param do.print Print the top genes associated with high/low loadings for
#' the PCs
#' @param pcs.print PCs to print genes for
#' @param pcs.store Number of PCs to store
#' @param pcs.compute Total Number of PCs to compute and store
#' @param genes.print Number of genes to print for each PC
#' @param \dots Additional arguments to be passed to prcomp
#' @return Returns Seurat object with an PCA embedding (object@@pca.rot) and
#' gene projection matrix (object@@pca.x). The PCA object itself is stored in
#' object@@pca.obj[[1]]
#' @importFrom irlba irlba
#' @export
PCAFast <- function(object, pc.genes = NULL, do.print = TRUE, pcs.print = 1:5, pcs.store = 40, 
                    pcs.compute = 20, genes.print = 30, ...) {
  return(DimReduction(object, reduction.type = "pcafast", genes.use = pc.genes, 
                      print.results = do.print, dims.store = pcs.store, dims.compute = pcs.compute, 
                      genes.print = genes.print,dims.print = pcs.print))
}


#' Run Independent Component Analysis on gene expression
#'
#' Run fastica algorithm from the ica package for ICA dimensionality reduction
#'
#'
#' @param object Seurat object
#' @param ic.genes Genes to use as input for ICA. Default is object@@var.genes
#' @param ics.store Number of ICs to store
#' @param ics.compute Number of ICs to compute
#' @param use.imputed Run ICA on imputed values (FALSE by default)
#' @param rev.ica By default, computes the dimensional reduction on the cell x gene matrix. 
#' Setting to true will compute it on the transpose (gene x cell matrix).
#' @param print.results Print the top genes associated with each dimension
#' @param ics.print ICs to print genes for
#' @param genes.print Number of genes to print for each IC
#' @param seed.use Random seed to use for fastica
#' @param \dots Additional arguments to be passed to fastica
#' @return Returns Seurat object with an ICA embedding (object@@ica.rot) and
#' gene projection matrix (object@@ica.x). The ICA object itself is stored in
#' object@@ica.obj[[1]]
#' @export
ICA <- function(object, ic.genes = NULL, ics.store = 40, ics.compute = 50, use.imputed = FALSE,
                rev.ica = FALSE, print.results = TRUE, ics.print = 1:5, genes.print = 50, 
                seed.use = 1, ...) {
  return(DimReduction(object, reduction.type = "ica", genes.use = ic.genes, dims.store = ics.compute, 
                      dims.compute = ics.compute, use.imputed = use.imputed, rev.reduction = rev.ica, 
                      print.results = print.results, dims.print = ics.print, genes.print = genes.print, 
                      seed.use = seed.use, ... ))
}


#' Run t-distributed Stochastic Neighbor Embedding
#'
#' Run t-SNE dimensionality reduction on selected features. Has the option of running in a reduced
#' dimensional space (i.e. spectral tSNE, recommended), or running based on a set of genes
#'
#'
#' @param object Seurat object
#' @param reduction.use Which dimensional reduction (e.g. PCA, ICA) to use for the tSNE. Default is PCA
#' @param cells.use Which cells to analyze (default, all cells)
#' @param dims.use Which dimensions to use as input features
#' @param genes.use If set, run the tSNE on this subset of genes
#' (instead of running on a set of reduced dimensions). Not set (NULL) by default
#' @param seed.use Random seed for the t-SNE
#' @param do.fast If TRUE, uses the Barnes-hut implementation, which runs
#' faster, but is less flexible
#' @param add.iter If an existing tSNE has already been computed, uses the
#' current tSNE to seed the algorithm and then adds additional iterations on top of this
#' @param dim.embed The dimensional space of the resulting tSNE embedding (default is 2).
#' For example, set to 3 for a 3d tSNE
#' @param \dots Additional arguments to the tSNE call. Most commonly used is
#' perplexity (expected number of neighbors default is 30)
#' @return Returns a Seurat object with a tSNE embedding in object@@dr$tsne@rotation
#' @importFrom Rtsne Rtsne
#' @importFrom tsne tsne
#' @export
RunTSNE <- function(object, reduction.use = "pca", cells.use = NULL, dims.use = 1:5, genes.use = NULL,
                    seed.use = 1, do.fast = FALSE, add.iter = 0, dim.embed = 2, ...) {
  if (is.null(genes.use)) {
    data.use <- GetDimReduction(object, reduction.type = reduction.use, slot = "rotation")[, dims.use]
    
  }
  if (!is.null(genes.use)) {
    if (length(object@scale.data) == 0){
      stop("Object@scale.data has not been set. Run ScaleData() and then retry.")
    }
    cells.use <- set.ifnull(cells.use, colnames(object@scale.data))
    genes.use <- ainb(genes.use, rownames(object@scale.data))
    data.use <- t(object@scale.data[genes.use, cells.use])
  }
  set.seed(seed.use)
  if (do.fast) {
    data.tsne <- Rtsne(as.matrix(data.use), dims = dim.embed, ...)
    data.tsne <- data.tsne$Y
  }
  else{
    data.tsne <- tsne(data.use, k = dim.embed, ...)
  }
  if (add.iter > 0) {
    data.tsne <- tsne(data.use, initial_config = as.matrix(data.tsne), max_iter = add.iter, ...)
  }
  colnames(data.tsne) <- paste0("tSNE_",1:ncol(data.tsne))
  rownames(data.tsne) <- rownames(data.use)
  object <- SetDimReduction(object, reduction.type = "tsne", slot = "rotation", new.data = data.tsne)
  object <- SetDimReduction(object, reduction.type = "tsne", slot = "key", new.data = "tSNE_")
  return(object)
}

#' Project Dimensional reduction onto full dataset
#'
#' Takes a pre-computed dimensional reduction (typically calculated on a subset of genes) and
#' projects this onto the entire dataset (all genes). Note that the cell
#' loadings (rotation matrices) remains unchanged, but now there are gene
#' scores for all genes.
#'
#'
#' @param object Seurat object
#' @param dims.print Number of dims to print genes for
#' @param dims.store Number of dims to store (default is 30)
#' @param genes.print Number of genes with highest/lowest loadings to print for
#' each PC
#' @param replace.dim Replace the existing data (overwite object@@dr$XXX@x), not done
#' by default.
#' @param do.center Center the dataset prior to projection (should be set to TRUE)
#' @param do.print Print top genes associated with the projected dimensions
#' @return Returns Seurat object with the projected values in object@@dr$XXX@x.full
#' @export
ProjectDim <- function(object, reduction.type = "pca", dims.print = 1:5, dims.store = 30, 
                       genes.print = 30, replace.dim = FALSE, do.center = FALSE, do.print = TRUE) {
  if (!(reduction.type%in%names(object@dr))) {
    stop(paste(reduction.type, " dimensional reduction has not been computed"))
  }
  data.use <- object@scale.data
  if (do.center) {
    data.use <- scale(as.matrix(object@scale.data), center = TRUE, scale = FALSE)
  }
  rotation <- GetDimReduction(object, reduction.type = reduction.type, slot = "rotation")
  new.full.x <- FastMatMult(data.use, rotation)
  rownames(new.full.x) <- rownames(object@scale.data)
  colnames(new.full.x) <- colnames(rotation)
  object <- SetDimReduction(object, reduction.type = reduction.type, slot = "x.full", 
                            new.data = new.full.x)
  if(replace.dim == TRUE){
    object <- SetDimReduction(object, reduction.type = reduction.type, slot = "x", 
                              new.data = new.full.x)
  }
  if(do.print) PrintDim(object, reduction.type = reduction.type, genes.print = genes.print, 
                        use.full = TRUE, dims.print = dims.print)
  return(object)
}


#' Project Principal Components Analysis onto full dataset
#'
#' Takes a pre-computed PCA (typically calculated on a subset of genes) and
#' projects this onto the entire dataset (all genes). Note that the cell
#' loadings (PCA rotation matrices) remains unchanged, but now there are gene
#' scores for all genes.
#'
#'
#' @param object Seurat object
#' @param do.print Print top genes associated with the projected PCs
#' @param pcs.print Number of PCs to print genes for
#' @param pcs.store Number of PCs to store (default is 30)
#' @param genes.print Number of genes with highest/lowest loadings to print for
#' each PC
#' @param replace.pc Replace the existing PCA (overwite object@@pca.x), not done
#' by default.
#' @param do.center Center the dataset prior to projection (should be set to TRUE)
#' @return Returns Seurat object with the projected PCA values in
#' object@@pca.x.full
#' @export
ProjectPCA <- function(object, do.print = TRUE, pcs.print = 1:5, pcs.store = 30, genes.print = 30,
                       replace.pc = FALSE, do.center = FALSE) {
  return(ProjectDim(object, reduction.type = "pca", dims.print = pcs.print, genes.print = 30, 
                    replace.dim = replace.pc, do.center = do.center, do.print = do.print,
                    dims.store = pcs.store))
}


#Internal, not documented for now
topGenesForDim=function(i, dim.scores, do.balanced=FALSE, num.genes=30, reduction.type="pca", 
                        key = "") {
  if (do.balanced) {
    num.genes=round(num.genes/2)
    sx=dim.scores[order(dim.scores[, i]),,drop=F]
    genes.1=(rownames(sx[1:num.genes,,drop=F]))
    genes.2=(rownames(sx[(nrow(sx)-num.genes+1):nrow(sx),,drop=F]))
    return(c(genes.1,genes.2))
  }
  if (!(do.balanced)) {
    sx=dim.scores[rev(order(abs(dim.scores[, i]))),,drop=F]
    genes.1=(rownames(sx[1:num.genes,,drop=F]))
    genes.1=genes.1[order(dim.scores[genes.1, i])]
    return(genes.1)
  }
}


#' Find genes with highest ICA scores
#'
#' Return a list of genes with the strongest contribution to a set of indepdendent components
#'
#' @param object Seurat object
#' @param ic.use Independent components to use
#' @param num.genes Number of genes to return
#' @param do.balanced Return an equal number of genes with both + and - IC scores.
#' @return Returns a vector of genes
#' @export
ICTopGenes <- function(object, ic.use = 1, num.genes = 30, use.full=F,do.balanced = FALSE) {
  return(DimTopGenes(object, dim.use = ic.use, reduction.type = "ica", use.full = use.full,num.genes = num.genes,
                     do.balanced = do.balanced))
}

#' Find genes with highest PCA scores
#'
#' Return a list of genes with the strongest contribution to a set of principal components
#'
#' @param object Seurat object
#' @param pc.use Principal components to use
#' @param num.genes Number of genes to return
#' @param use.full Use the full PCA (projected PCA). Default i s FALSE
#' @param do.balanced Return an equal number of genes with both + and - PC scores.
#' @return Returns a vector of genes
#' @export
PCTopGenes <- function(object, pc.use = 1, num.genes = 30, use.full = FALSE, do.balanced = FALSE) {
  return(DimTopGenes(object, dim.use = pc.use, reduction.type = "pca", num.genes = num.genes,
                     use.full = use.full, do.balanced = do.balanced))
}

#' Find cells with highest PCA scores
#'
#' Return a list of genes with the strongest contribution to a set of principal components
#'
#' @param object Seurat object
#' @param pc.use Principal component to use
#' @param num.cells Number of cells to return
#' @param do.balanced Return an equal number of cells with both + and - PC scores.
#' @return Returns a vector of cells
#' @export
PCTopCells <- function(object, pc.use = 1, num.cells = NULL, do.balanced = FALSE) {
  return(DimTopCells(object, dim.use = pc.use, reduction.type = "pca", num.cells = num.cells, 
                     do.balanced = do.balanced))
}

#' Find cells with highest ICA scores
#'
#' Return a list of genes with the strongest contribution to a set of principal components
#'
#' @param object Seurat object
#' @param ic.use Independent component to use
#' @param num.cells Number of cells to return
#' @param do.balanced Return an equal number of cells with both + and - PC scores.
#' @return Returns a vector of cells
#' @export
ICTopCells <- function(object, ic.use = 1, num.cells = NULL, do.balanced = FALSE) {
  return(DimTopCells(object, dim.use = ic.use, reduction.type = "ica", num.cells = num.cells, 
                     do.balanced = do.balanced))
}

#' Find cells with highest scores for a given dimensional reduction technique
#'
#' Return a list of genes with the strongest contribution to a set of components
#'
#' @param object Seurat object
#' @param reduction.type Dimensional reduction to find the highest score for
#' @param dim.use Components to use
#' @param num.cells Number of cells to return
#' @param do.balanced Return an equal number of cells with both + and - scores.
#' @return Returns a vector of cells
#' @export
DimTopCells <- function(object, dim.use = 1, reduction.type = "pca", num.cells = NULL, 
                        do.balanced=FALSE) {
  #note that we use topGenesForDim, but it still works
  #error checking
  if (!(reduction.type%in%names(object@dr))) {
    stop(paste(reduction.type, " dimensional reduction has not been computed"))
  }
  num.cells <- set.ifnull(num.cells, length(object@cell.names))
  dim.scores <- GetDimReduction(object, reduction.type = reduction.type, slot = "rotation")
  key <- GetDimReduction(object, reduction.type = reduction.type, slot = "key")
  i <- dim.use
  dim.top.cells <- unique(unlist(lapply(i, topGenesForDim, dim.scores, do.balanced, num.cells, 
                                        reduction.type, key)))
  return(dim.top.cells)
}

#' Find genes with highest scores for a given dimensional reduction technique
#'
#' Return a list of genes with the strongest contribution to a set of components
#'
#' @param object Seurat object
#' @param reduction.type Dimensional reduction to find the highest score for
#' @param pc.use Components to use
#' @param num.genes Number of genes to return
#' @param use.full Use the full PCA (projected PCA). Default i s FALSE
#' @param do.balanced Return an equal number of genes with both + and - scores.
#' @return Returns a vector of genes
#' @export
DimTopGenes <- function(object, dim.use = 1, reduction.type = "pca", num.genes = 30, use.full = F,
                        do.balanced = FALSE) {
  #note that we use topGenesForDim, but it still works
  #error checking
  if (!(reduction.type%in%names(object@dr))) {
    stop(paste(reduction.type, " dimensional reduction has not been computed"))
  }
  dim.scores <- GetDimReduction(object, reduction.type = reduction.type, slot = "x")
  if (use.full) dim.scores <- GetDimReduction(object, reduction.type = reduction.type, slot = "x.full")
  i <- dim.use
  key <- GetDimReduction(object, reduction.type = reduction.type, slot = "key")
  dim.top.genes <- unique(unlist(lapply(i, topGenesForDim, dim.scores, do.balanced, num.genes,
                                        reduction.type, key)))
  return(dim.top.genes)
}



#' Dimensional reduction heatmap
#'
#' Draws a heatmap focusing on a principal component. Both cells and genes are sorted by their 
#' principal component scores. Allows for nice visualization of sources of heterogeneity in the dataset.
#'
#' @inheritParams DoHeatmap
#' @inheritParams PCTopGenes
#' @inheritParams VizPCA
#' @param cells.use A list of cells to plot. If numeric, just plots the top cells.
#' @param use.scale Default is TRUE: plot scaled data. If FALSE, plot raw data on the heatmap.
#' @param label.columns Whether to label the columns. Default is TRUE for 1 PC, FALSE for > 1 PC
#' @return If do.return==TRUE, a matrix of scaled values which would be passed
#' to heatmap.2. Otherwise, no return value, only a graphical output
#' @export
DimHeatmap <- function(object, reduction.type = "pca", dim.use = 1, cells.use = NULL, 
                       num.genes = 30, use.full = FALSE, disp.min = -2.5, disp.max = 2.5,
                       do.return = FALSE, col.use = pyCols, use.scale = TRUE,
                       do.balanced = FALSE, remove.key = FALSE, label.columns=NULL, ...){
  
  num.row <- floor(length(dim.use) / 3.01) + 1
  orig_par <- par()$mfrow
  par(mfrow=c(num.row, min(length(dim.use), 3)))
  cells <- cells.use
  plots <- c()
  
  if (is.null(label.columns)){
    if (length(dim.use) > 1){
      label.columns <- FALSE
    }
    else{
      label.columns <- TRUE
    }
  }
  
  for(ndim in dim.use){
    if (is.numeric((cells))) {
      cells.use <- DimTopCells(object = object,dim.use = ndim, reduction.type = reduction.type, 
                               num.cells = cells,do.balanced = do.balanced)
    }
    else {
      cells.use <- set.ifnull(cells, object@cell.names)
    }
    genes.use <- rev(DimTopGenes(object = object,dim.use = ndim, reduction.type = reduction.type,
                                 num.genes = num.genes, use.full = use.full, 
                                 do.balanced = do.balanced))
    dim.scores <- GetDimReduction(object, reduction.type = reduction.type, slot = "rotation")
    dim.key <- GetDimReduction(object, reduction.type = reduction.type, slot = "key")
    cells.ordered <- cells.use[order(dim.scores[cells.use, paste(dim.key, ndim, sep = "")])]
    data.use <- object@scale.data[genes.use, cells.ordered]
    data.use <- minmax(data.use, min = disp.min, max = disp.max)
    if (!(use.scale)) data.use <- as.matrix(object@data[genes.use, cells.ordered])
    vline.use <- NULL
    hmTitle <- paste(dim.key, ndim)
    if (remove.key || length(dim.use) > 1){
      hmFunction <- "heatmap2NoKey(data.use, Rowv = NA, Colv = NA, trace = \"none\", col = col.use, dimTitle = hmTitle, "
    }
    else{
      hmFunction <- "heatmap.2(data.use,Rowv=NA,Colv=NA,trace = \"none\",col=col.use, dimTitle = hmTitle, "
    }
    
    if (!label.columns){
      
      hmFunction <- paste(hmFunction, "labCol=\"\", ", sep="")
    }
    hmFunction <- paste(hmFunction, "...)", sep="")
    #print(hmFunction)
    eval(parse(text = hmFunction))
  }
  if (do.return) {
    return(data.use)
  }
  # reset graphics parameters
  par(mfrow = orig_par)
}

PlotDim <- function(ndim, object, reduction.type, use.scaled, use.full, cells.use, num.genes, 
                    group.by, disp.min, disp.max, col.low, col.mid, col.high, slim.col.label, 
                    do.balanced, remove.key, cex.col, cex.row, group.label.loc, group.label.rot, 
                    group.cex, group.spacing){
  if (is.numeric((cells.use))) {
    cells.use <- DimTopCells(object = object, dim.use = ndim, reduction.type = reduction.type, 
                             num.cells = cells.use ,do.balanced = do.balanced)
  }
  else {
    cells.use <- set.ifnull(cells.use, object@cell.names)
  }
  genes.use <- rev(DimTopGenes(object = object,dim.use = ndim, reduction.type = reduction.type,
                               num.genes = num.genes, use.full = use.full, 
                               do.balanced = do.balanced))
  dim.scores <- GetDimReduction(object, reduction.type = reduction.type, slot = "rotation")
  dim.key <- GetDimReduction(object, reduction.type = reduction.type, slot = "key")
  cells.ordered <- cells.use[order(dim.scores[cells.use, paste0(dim.key, ndim)])]
  data.use <- object@scale.data[genes.use, cells.ordered]
  data.use <- minmax(data.use, min = disp.min, max = disp.max)
  if (!(use.scaled)) data.use <- as.matrix(object@data[genes.use, cells.ordered])
  return(DoHeatmapGG(object, data.use = data.use, cells.use = cells.use, 
              genes.use = genes.use, group.by = group.by, disp.min = disp.min, disp.max = disp.max,
              col.low = col.low, col.mid = col.mid, col.high = col.high, 
              slim.col.label = slim.col.label, remove.key = remove.key, cex.col = cex.col, 
              cex.row = cex.row, group.label.loc = group.label.loc, group.label.rot = group.label.rot, 
              group.cex = group.cex, group.spacing = group.spacing, title = paste0(dim.key, ndim), 
              do.plot = FALSE))
}


#' Principal component heatmap
#'
#' Draws a heatmap focusing on a principal component. Both cells and genes are sorted by their principal component scores.
#' Allows for nice visualization of sources of heterogeneity in the dataset.
#'
#' @inheritParams DoHeatmap
#' @inheritParams PCTopGenes
#' @inheritParams VizPCA
#' @param cells.use A list of cells to plot. If numeric, just plots the top cells.
#' @param use.scale Default is TRUE: plot scaled data. If FALSE, plot raw data on the heatmap.
#' @param label.columns Whether to label the columns. Default is TRUE for 1 PC, FALSE for > 1 PC
#' @return If do.return==TRUE, a matrix of scaled values which would be passed
#' to heatmap.2. Otherwise, no return value, only a graphical output
#' @export
PCHeatmap <- function(object, pc.use = 1, cells.use = NULL, num.genes = 30, use.full = FALSE, 
                      disp.min = -2.5, disp.max = 2.5, do.return = FALSE, col.use=pyCols,
                      use.scale = TRUE, do.balanced = FALSE, remove.key = FALSE, 
                      label.columns = NULL, ...) {
  
  return(DimHeatmap(object, reduction.type = "pca", dim.use = pc.use, cells.use = cells.use, 
                    num.genes = num.genes, use.full = use.full, disp.min = disp.min, disp.max = disp.max,
                    do.return = do.return, col.use = col.use, use.scale = use.scale, 
                    do.balanced = do.balanced, remove.key = remove.key, label.columns = label.columns, ...))
  
}


#' Independent component heatmap
#'
#' Draws a heatmap focusing on a principal component. Both cells and genes are sorted by their 
#' principal component scores. Allows for nice visualization of sources of heterogeneity 
#' in the dataset."()
#'
#' @inheritParams DoHeatmap
#' @inheritParams ICTopGenes
#' @inheritParams VizICA
#' @param use.scale Default is TRUE: plot scaled data. If FALSE, plot raw data on the heatmap.
#' @return If do.return==TRUE, a matrix of scaled values which would be passed
#' to heatmap.2. Otherwise, no return value, only a graphical output
#' @export
ICHeatmap <- function(object, ic.use = 1, cells.use = NULL, num.genes = 30, disp.min = -2.5, 
                      disp.max = 2.5, do.return = FALSE, col.use = pyCols, use.scale = TRUE, 
                      do.balanced = FALSE, remove.key = FALSE, label.columns = NULL, ...) {
  
  return(DimHeatmap(object, reduction.type = "ica", dim.use = ic.use, cells.use = cells.use, 
                    num.genes = num.genes, disp.min = disp.min, disp.max = disp.max, 
                    do.return = do.return, col.use = col.use, use.scale = use.scale, 
                    do.balanced = do.balanced, remove.key = remove.key, 
                    label.columns = label.columns, ...))
}

#' Print the results of a dimensional reduction analysis
#' 
#' Prints a set of genes that most strongly define a set of components
#' 
#' @param object Seurat object
#' @param reduction.type Reduction technique to print results for
#' @param dims.print Number of dimensions to display
#' @param genes.print Number of genes to display
#' @param use.full Use full PCA (i.e. the projected PCA, by default FALSE)
#' @return Set of genes defining the components 
#' @export 
PrintDim <- function(object, reduction.type = "pca", dims.print = 1:5, genes.print = 30, 
                     use.full = FALSE){
  
  if(length(GetDimReduction(object, reduction.type = reduction.type, slot = "x.full")) == 0 && use.full){
    warning("Dimensions have not been projected. Setting use.full = FALSE")
    use.full <- FALSE
  }
  for(i in dims.print) {
    code <- paste0(GetDimReduction(object, reduction.type = reduction.type, slot = "key"), i)
    sx <- DimTopGenes(object, dim.use = i, reduction.type = reduction.type, 
                      num.genes = genes.print * 2, use.full = use.full, do.balanced = TRUE)
    print(code)
    print((sx[1:genes.print]))
    print ("")
    
    print(rev((sx[(length(sx)-genes.print+1):length(sx)])))
    print ("")
    print ("")
  }
}

#' Print the results of a PCA analysis
#'
#' Prints a set of genes that most strongly define a set of principal components
#'
#' @inheritParams VizPCA
#' @param pcs.print Set of PCs to print genes for
#' @param genes.print Number of genes to print for each PC
#' @return Only text output
#' @export
PrintICA <- function(object, ics.print = 1:5, genes.print = 30, use.full = FALSE) {
  PrintDim(object, reduction.type = "ica", dims.print = ics.print, genes.print = genes.print, 
           use.full = use.full)
}

#' Print the results of a PCA analysis
#'
#' Prints a set of genes that most strongly define a set of principal components
#'
#' @inheritParams VizPCA
#' @param pcs.print Set of PCs to print genes for
#' @param genes.print Number of genes to print for each PC
#' @return Only text output
#' @export
PrintPCA <- function(object, pcs.print = 1:5, genes.print = 30, use.full = FALSE) {
  PrintDim(object, reduction.type = "pca", dims.print = pcs.print, genes.print = genes.print, 
           use.full = use.full)
}

#' Visualize Dimensional Reduction genes
#'
#' Visualize top genes associated with reduction components
#' 
#' @param object Seurat object
#' @param reduction.type Reduction technique to visualize results for
#' @param dims.use Number of dimensions to display
#' @param num.genes Number of genes to display
#' @param use.full Use reduction values for full dataset (i.e. projected dimensional reduction values)
#' @param font.size Font size
#' @param nCol Number of columns to display
#' @param do.balanced Return an equal number of genes with + and - scores. If FALSE (default), returns
#' the top genes ranked by the scores absolute values
#' @return Graphical, no return value
#' @export 
VizDimReduction <- function(object, reduction.type = "pca", dims.use = 1:5, num.genes = 30, 
                            use.full = FALSE, font.size = 0.5, nCol = NULL, do.balanced = FALSE) {
  dim.scores <- GetDimReduction(object, reduction.type = reduction.type, slot = "x")
  if(use.full == TRUE) dim.scores <- GetDimReduction(object, reduction.type = reduction.type, "x.full")
  
  if(is.null(nCol)) {
    nCol <- 2
    if (length(dims.use) > 6) nCol <- 3
    if (length(dims.use) > 9) nCol <- 4
  }
  num.row <- floor(length(dims.use) / nCol - 1e-5) + 1
  par(mfrow = c(num.row, nCol))
  
  for(i in dims.use){
    subset.use <- dim.scores[DimTopGenes(object, i, reduction.type, num.genes, use.full, do.balanced), ]
    plot(subset.use[, i], 1:nrow(subset.use) ,pch = 16, col = "blue", xlab = paste0("PC", i), 
         yaxt="n", ylab="")
    axis(2, at = 1:nrow(subset.use), labels = rownames(subset.use), las = 1, cex.axis = font.size)
  }
  rp()
}


#' Visualize PCA genes
#'
#' Visualize top genes associated with principal components
#'
#'
#' @param object Seurat object
#' @param pcs.use Number of PCs to display
#' @param num.genes Number of genes to display
#' @param use.full Use full PCA (i.e. the projected PCA, by default FALSE)
#' @param font.size Font size
#' @param nCol Number of columns to display
#' @param do.balanced Return an equal number of genes with both + and - PC scores.
#' If FALSE (by default), returns the top genes ranked by the score's absolute values
#' @return Graphical, no return value
#' @export
VizPCA <- function(object, pcs.use = 1:5, num.genes = 30, use.full = FALSE, font.size = 0.5, 
                   nCol = NULL, do.balanced = FALSE){
  VizDimReduction(object, reduction.type = "pca", dims.use = pcs.use, num.genes = num.genes, 
                  use.full = use.full, font.size = font.size, nCol = nCol, do.balanced = do.balanced)
}


#' Visualize ICA genes
#'
#' Visualize top genes associated with principal components
#'
#'
#' @param object Seurat object
#' @param ics.use Number of ICs to display
#' @param num.genes Number of genes to display
#' @param use.full Use full ICA (i.e. the projected ICA, by default FALSE)
#' @param font.size Font size
#' @param nCol Number of columns to display
#' @param do.balanced Return an equal number of genes with both + and - IC scores.
#' If FALSE (by default), returns the top genes ranked by the score's absolute values
#' @return Graphical, no return value
#' @export
VizICA <- function(object, ics.use = 1:5, num.genes = 30, use.full = FALSE, font.size = 0.5, 
                   nCol = NULL, do.balanced = FALSE) {
  VizDimReduction(object, reduction.type = "ica", dims.use = pcs.use, num.genes = num.genes, 
                  use.full = use.full, font.size = font.size, nCol = nCol, do.balanced = do.balanced)
}


#' Dimensional reduction plot
#'
#' Graphs the output of a dimensional reduction technique (PCA by default).
#' Cells are colored by their identity class.
#'
#'
#' @param object Seurat object
#' @param reduction.use Which dimensionality reduction to use. Default is
#' "pca", can also be "tsne", or "ica", assuming these are precomputed.
#' @param dim.1 Dimension for x-axis (default 1)
#' @param dim.2 Dimension for y-axis (default 2)
#' @param cells.use Vector of cells to plot (default is all cells)
#' @param pt.size Adjust point size for plotting
#' @param do.return Return a ggplot2 object (default : FALSE)
#' @param do.bare Do only minimal formatting (default : FALSE)
#' @param cols.use Vector of colors, each color corresponds to an identity
#' class. By default, ggplot assigns colors.
#' @param group.by Group (color) cells in different ways (for example, orig.ident)
#' @param pt.shape If NULL, all points are circles (default). You can specify any
#' cell attribute (that can be pulled with FetchData) allowing for both different colors and 
#' different shapes on cells.
#' @param do.label Whether to label the clusters
#' @param label.size Sets size of labels
#' @param no.legend Setting to TRUE will remove the legend
#' @return If do.return==TRUE, returns a ggplot2 object. Otherwise, only
#' graphical output.
#' @importFrom dplyr summarize group_by
#' @export
DimPlot <- function(object, reduction.use = "pca", dim.1 = 1, dim.2 = 2, cells.use = NULL, 
                    pt.size = 3, do.return = FALSE, do.bare = FALSE, cols.use = NULL, 
                    group.by = "ident", pt.shape = NULL, do.label = FALSE, label.size = 1, 
                    no.legend = FALSE) {
  if(length(GetDimReduction(object, reduction.type = reduction.use, slot = "rotation")) == 0) {
    stop(paste0(reduction.use, "has not been run for this object yet."))
  }
  cells.use <- set.ifnull(cells.use, colnames(object@data))
  dim.code <- GetDimReduction(object, reduction.type = reduction.use, slot = "key")
  dim.codes <- paste0(dim.code, c(dim.1, dim.2))
  data.plot <- as.data.frame(GetDimReduction(object, reduction.type = reduction.use, slot = "rotation"))
  cells.use=intersect(cells.use,rownames(data.plot))
  data.plot=data.plot[cells.use, dim.codes]
  ident.use <- as.factor(object@ident[cells.use])
  if (group.by != "ident") ident.use <- as.factor(FetchData(object, group.by)[cells.use, 1])
  data.plot$ident <- ident.use
  data.plot$x=data.plot[, dim.codes[1]]
  data.plot$y=data.plot[, dim.codes[2]]
  data.plot$pt.size <- pt.size
  p <- ggplot(data.plot, aes(x = x,y = y)) + geom_point(aes(colour = factor(ident)), size = pt.size)
  if (!is.null(pt.shape)) {
    shape.val <- FetchData(object, pt.shape)[cells.use, 1]
    if (is.numeric(shape.val)) {
      shape.val <- cut(shape.val, breaks = 5)
    }
    data.plot[,"pt.shape"] <- shape.val
    p <- ggplot(data.plot, aes(x = x, y = y)) + geom_point(aes(colour = factor(ident), 
                                                               shape = factor(pt.shape)), size = pt.size)
  }
  if (!is.null(cols.use)) {
    p <- p + scale_colour_manual(values = cols.use)
  }
  p2 <- p + xlab(dim.codes[[1]]) + ylab(dim.codes[[2]]) + scale_size(range = c(pt.size, pt.size))
  p3 <- p2 + gg.xax() + gg.yax() + gg.legend.pts(6) + gg.legend.text(12) + no.legend.title +
    theme_bw() + nogrid
  p3 <- p3 + theme(legend.title = element_blank())
  if (do.label) {
    data.plot %>% dplyr::group_by(ident) %>% summarize(x = median(x), y = median(y)) -> centers
    p3 <- p3 + geom_point(data = centers, aes(x = x, y = y), size = 0, alpha = 0) + 
      geom_text(data = centers, aes(label = ident), size = label.size)
  }
  if (no.legend) {
    p3 <- p3 + theme(legend.position = "none")
  }
  if (do.return) {
    if (do.bare) return(p)
    return(p3)
  }
  if (do.bare) print(p)
  else print(p3)
}


#' Plot PCA map
#'
#' Graphs the output of a PCA analysis
#' Cells are colored by their identity class.
#'
#' This function is a wrapper for DimPlot. See ?DimPlot for a full list of possible
#' arguments which can be passed in here.
#'
#' @param object Seurat object
#' @param \dots Additional parameters to DimPlot, for example, which dimensions to plot.
#' @export
PCAPlot <- function(object,...) {
  return(DimPlot(object, reduction.use = "pca", label.size = 6, ...))
}

#' Plot Diffusion map
#'
#' Graphs the output of a Diffusion map analysis
#' Cells are colored by their identity class.
#'
#' This function is a wrapper for DimPlot. See ?DimPlot for a full list of possible
#' arguments which can be passed in here.
#'
#' @param object Seurat object
#' @param \dots Additional parameters to DimPlot, for example, which dimensions to plot.
#' @export
DMPlot <- function(object,...) {
  return(DimPlot(object, reduction.use = "dm", label.size = 6, ...))
}

#' Plot ICA map
#'
#' Graphs the output of a ICA analysis
#' Cells are colored by their identity class.
#'
#' This function is a wrapper for DimPlot. See ?DimPlot for a full list of possible
#' arguments which can be passed in here.
#'
#' @param object Seurat object
#' @param \dots Additional parameters to DimPlot, for example, which dimensions to plot.
#' @export
ICAPlot <- function(object,...) {
  return(DimPlot(object,reduction.use = "ica",...))
}


#' Plot tSNE map
#'
#' Graphs the output of a tSNE analysis
#' Cells are colored by their identity class.
#'
#' This function is a wrapper for DimPlot. See ?DimPlot for a full list of possible
#' arguments which can be passed in here.
#'
#' @param object Seurat object
#' @param do.label FALSE by default. If TRUE, plots an alternate view where the center of each
#' cluster is labeled
#' @param pt.size Set the point size
#' @param label.size Set the size of the text labels
#' @param cells.use Vector of cell names to use in the plot.
#' @param colors.use Manually set the color palette to use for the points
#' @param \dots Additional parameters to DimPlot, for example, which dimensions to plot.
#' @seealso DimPlot
#' @export
TSNEPlot <- function(object, do.label = FALSE, pt.size=1, label.size=4, cells.use = NULL, colors.use = NULL,...) {
  return(DimPlot(object, reduction.use = "tsne", cells.use = cells.use, pt.size = pt.size, 
                 do.label = do.label, label.size = label.size, cols.use = colors.use, ...))
}



#' Quickly Pick Relevant Dimensions
#'
#' Plots the standard deviations (or approximate singular values if running PCAFast)
#' of the principle components for easy identification of an elbow in the graph. 
#' This elbow often corresponds well with the significant dims and is much faster to run than 
#' Jackstraw
#'
#'
#' @param object Seurat object
#' @param reduction.type  Type of dimensional reduction to plot data for
#' @param dims.plot Number of dimensions to plot sd for 
#' @param xlab X axis label
#' @param ylab Y axis label
#' @param title Plot title
#' @return Returns ggplot object
#' @export
DimElbowPlot <- function(object, reduction.type = "pca", dims.plot = 20, xlab = "", ylab = "", 
                         title = "") {
  if(length(GetDimReduction(object, reduction.type = reduction.type, slot = "sdev")) == 0) {
    stop(paste0("No standard deviation info stored for", reduction.use))
  }
  data.use <- GetDimReduction(object, reduction.type = reduction.type, slot = "sdev")
  if (length(data.use) < dims.plot) {
    warning(paste("The object only has information for", length(data.use), "PCs." ))
    dims.plot <- length(data.use)
  }
  data.use <- data.use[1:dims.plot]
  dims <- 1:length(data.use)
  data.plot <- data.frame(dims, data.use)
  if(reduction.type == "pca"){
    plot <- ggplot(data.plot, aes(dims, data.use)) + geom_point() + labs(y = "Standard Deviation of PC", 
                                                                         x = "PC", title = title)
  }
  else if(reduction.type == "pcafast"){
    plot <- ggplot(data.plot, aes(dims, data.use)) + geom_point() + labs(y = "Eigen values of PC", 
                                                                         x = "Eigen value", 
                                                                         title = title)
  }
  else if(reduction.type == "ica"){
    plot <- ggplot(data.plot, aes(dims, data.use)) + geom_point() + labs(y = "Standard Deviation of IC", 
                                                                         x = "IC", title = title)
  }
  else{
    plot <- ggplot(data.plot, aes(dims, data.use)) + geom_point() + labs(y = ylab, x = xlab, 
                                                                         title = title)
  }
  return (plot)
}

#' Quickly Pick Relevant PCs
#'
#' Plots the standard deviations (or approximate singular values if running PCAFast)
#' of the principle components for easy identification of an elbow in the graph. 
#' This elbow often corresponds well with the significant PCs and is much faster to run.
#'
#'
#' @param object Seurat object
#' @param num.pc Number of PCs to plot
#' @return Returns ggplot object
#' @export
PCElbowPlot <- function(object, num.pc = 20) {
  return(DimElbowPlot(object, reduction.type = "pca", dims.plot = num.pc))
}

#' Perform Canonical Correlation Analysis
#'
#' Runs a canonical correlation analysis using the sparse implementation of CCA from the PMA package.
#'
#'
#' @param object Seurat object
#' @param object2 Optional second object. If object2 is passed, object with be considered as group1
#' and object2 as group2.
#' @param group1 First set of cells (or IDs) for CCA
#' @param group2 Second set of cells (or IDs) for CCA
#' @param group.by Factor to group by (column vector stored in object@@data.info)
#' @param num.cc Number of canonical vectors to calculate
#' @param genes.use Set of genes to use in CCA. Default is object@@var.genes. If two objects are given,
#' the default is the union of both variable gene sets that are also present in both objects.
#' @param scale.data Use the scaled data from the object
#' @param PMA.version Use PMA version of CCA
#' @param rescale.groups Rescale each set of cells independently
#' @importFrom PMA CCA
#' @return Returns Seurat object with the CCA stored in the @@dr$cca slot. If one object is passed,
#' the same object is returned. If two are passed, a combined object is returned. 
#' @export
RunCCA <- function(object, object2, group1, group2, group.by, num.cc = 20, genes.use, 
                   scale.data = TRUE, rescale.groups = FALSE, PMA.version = FALSE) {
  if(!missing(object2) && (!missing(group1) || !missing(group2))){
    warning("Both object2 and group set. Continuing with objects defining the groups")
  }
  if(!missing(object2)){
    if(missing(genes.use)){
      genes.use <- union(object@var.genes, object2@var.genes)
      if(length(genes.use) == 0) stop("No variable genes present. Run MeanVarPlot and retry")
    }
    if(scale.data){
      possible.genes <- intersect(rownames(object@scale.data), rownames(object2@scale.data))
      genes.use <- genes.use[genes.use %in% possible.genes]
      data.use1 <- object@scale.data[genes.use, ]
      data.use2 <- object2@scale.data[genes.use, ]
    }
    else{
      possible.genes <- intersect(rownames(object@data), rownames(object2@data))
      genes.use <- genes.use[genes.use %in% possible.genes]
      data.use1 <- object@data[genes.use, ]
      data.use2 <- object2@data[genes.use, ]
    }
    if(length(genes.use) == 0) stop("0 valid genes in genes.use")
  }
  else{
    if(missing(group1)) stop("group1 not set")
    if(missing(group2)) stop("group2 not set")
    if(!missing(group.by)){
      if(!group.by %in% colnames(object@data.info)) stop("invalid group.by parameter")
    }
    if(missing(genes.use)){
      genes.use <- object@var.genes
      if(length(genes.use) == 0) stop("No variable genes present. Run MeanVarPlot and retry")
    }
    if(missing(group.by)){
      cells.1 <- CheckGroup(object, group = group1, group.id = "group1")
      cells.2 <- CheckGroup(object, group = group2, group.id = "group2")
    }
    else{
      object.current.ids <- object@ident
      object <- SetAllIdent(object, id = group.by)
      cells.1 <- CheckGroup(object, group = group1, group.id = "group1")
      cells.2 <- CheckGroup(object, group = group2, group.id = "group2")
      object <- SetIdent(object, cells.use = object@cell.names, ident.use = object.current.ids)
    }
    if(scale.data){
      if(rescale.groups){
        data.use1 <- ScaleData(object, data.use = object@data[genes.use, cells.1])
        data.use1 <- data.use1@scale.data
        data.use2 <- ScaleData(object, data.use = object@data[genes.use, cells.2])
        data.use2 <- data.use2@scale.data
      }
      else{
        data.use1 <- object@scale.data[genes.use, cells.1]
        data.use2 <- object@scale.data[genes.use, cells.2]
      }
    }
    else{
      data.use1 <- object@data[genes.use, cells.1]
      data.use2 <- object@data[genes.use, cells.2]
    }
  }
  if(PMA.version){
    cca.results <- CCA(data.use1, data.use2, typex = "standard", typez = "standard", K = num.cc, 
                       penaltyz = 1, penaltyx = 1, trace = F)
  }
  else{
    cca.results <- SparseCanonCor(data.use1, data.use2, standardize = TRUE, k = num.cc)
  }
  cca.data <- rbind(cca.results$u, cca.results$v)
  rownames(cca.data) <- c(colnames(data.use1), colnames(data.use2))
  colnames(cca.data) <- paste0("CC", 1:num.cc)

  if(!missing(object2)){
    cat("Merging objects\n", file = stderr())
    combined.object <- MergeSeurat(object, object2, do.scale = F, do.center = F)
    combined.object@var.genes=genes.use
    combined.object <- FastScaleData(combined.object)
    combined.object <- SetDimReduction(combined.object, reduction.type = "cca", slot = "rotation",
                                       new.data = cca.data)
    combined.object <- SetDimReduction(combined.object, reduction.type = "cca", slot = "key",
                                       new.data = "CC")
    combined.object=ProjectDim(combined.object,reduction.type = "cca",do.print = F)
    combined.object=SetDimReduction(combined.object,reduction.type = "cca","x",new.data = DimX(combined.object,reduction.type = "cca",use.full = T,genes.use = genes.use))
    return(combined.object)
  }
  else{
    object <- SetDimReduction(object, reduction.type = "cca", slot = "rotation", new.data = cca.data)
    object <- SetDimReduction(object, reduction.type = "cca", slot = "key", new.data = "CC")
    
    object=ProjectDim(object,reduction.type = "cca")
    return(object)
  }
}

CheckGroup <- function(object, group, group.id){
  if(all(group %in% unique(object@ident))) {
    cells.use <- WhichCells(object, ident = group)
  }
  else{
    if(all(group %in% object@cell.names)){
      cells.use <- group
    }
    else{
      stop(paste0(group.id, " must be either a vector of valid cell names or idents"))
    }
  }
  return(cells.use)
}

# Note: not really "sparse" yet -- could add in the penalties discussed in Witten 

SparseCanonCor <- function(mat1, mat2, standardize = TRUE, k = 20){
  if(standardize){
    mat1 <- Standardize(mat1, FALSE)
    mat2 <- Standardize(mat2, FALSE)
  }
  u.mat <- matrix(nrow = ncol(mat1), ncol = k)
  v.mat <- matrix(nrow = ncol(mat2), ncol = k)
  d <- numeric(k)
  v.init <- svd(t(mat1) %*% mat2)$v
  for (i in 1:k){
    v <- v.init[,i]
    u <- t(t(mat2 %*% v) %*% mat1)
    u <- u / norm(as.vector(u), "2")
    v = t(t(mat1 %*% u) %*% mat2)
    v = v / norm(as.vector(v), "2")
    u.mat[, i] <- u
    v.mat[, i] <- v
    d[i] <- sum((mat1 %*% u) * (mat2 %*% v))
    mat1 <- rbind(mat1, sqrt(d[i]) * t(u))
    mat2 <- rbind(mat2, -sqrt(d[i]) * t(v))
  }
  return(list(u = u.mat, v = v.mat, d = d))
}



#' Shift dim scores
#'
#' Shifts dim scores so that they line up across grouping variable (only implemented for case with
#' 2 categories in grouping.var)
#'
#'
#' @param object Seurat object
#' @param group.var Name of the grouping variable for which to shift the scores 
#' @param reduction.type reduction to shift scores for 
#' @param dims.shift Dims to shift, default is all
#' @param ds.amt percent of cells in each group for downsampling
#' @return Returns Seurat object with the dims shifted, stored in object@@dr$reduction.type.shifted
#' @export
ShiftDim <- function(object, grouping.var, reduction.type, dims.shift, ds.amt = 0.4,constrain.var=T) {
  groups <- as.vector(unique(FetchData(object, grouping.var)[, 1]))
  cells.1.all <- WhichCells(object, subset.name = grouping.var, accept.value = groups[1])
  cells.2.all <- WhichCells(object, subset.name = grouping.var, accept.value = groups[2])
  if(missing(dims.shift)){
    dims.shift <- 1:ncol(DimRot(object, reduction.type = reduction.type))
  }
  shifted.rot <- DimRot(object, reduction.type = reduction.type)
  if (constrain.var) {
    prior.proj=FastMatMult(t(shifted.rot),t(object@scale.data[rownames(DimX(object,reduction.type = reduction.type)),]))
    prior.var=apply(prior.proj,1,var)
  }
  object.all=object
  object=SetAllIdent(object,grouping.var)
  object=SubsetData(object,max.cells.per.ident = min(table(object@ident)))
  cells.1 <- WhichCells(object, subset.name = grouping.var, accept.value = groups[1])
  cells.2 <- WhichCells(object, subset.name = grouping.var, accept.value = groups[2])
  
  #grrr, we had discussed on the phone that:
  #1. optim does not work here (need a grid search)
  #2. need to downsample first to have equal representation from both
  
  for (i in dims.shift){
    dim.1 <- DDDS(DimRot(object, reduction.type = reduction.type, cells.use = cells.1, dims.use = i), 
                  amt = ds.amt * length(cells.1))
    dim.2 <- DDDS(DimRot(object, reduction.type = reduction.type, cells.use = cells.2, dims.use = i),
                  amt = ds.amt * length(cells.2))
    cells.use = c(names(dim.1), names(dim.2))
    #shifting.params <- optim(c(0, 1), fn = ShiftObj, 
    #                         drrot = DimRot(object, reduction.type = reduction.type, dims.use = i, 
    #                                        cells.use = cells.use), 
    #                         ids = object@ident[cells.use], cells.1 = names(dim.1), fun.opt = 2)$par
    
    dimrot.use=DimRot(object, reduction.type = reduction.type, dims.use = i, cells.use = cells.use)
    shifting.params=gridSearch(ShiftObj, drrot = dimrot.use, ids = object@ident[cells.use], cells.1 = names(dim.1), fun.opt = 2,lower = c(-0.05,0.05),upper = c(0.5,2),npar = 2,n = 50)$minlevels
    shifted.rot[cells.1.all, i] <- shifting.params[1] + shifted.rot[cells.1.all, i] * shifting.params[2]
    colnames(shifted.rot)[i] <- paste0("S", GetDimReduction(object, reduction.type = reduction.type, 
                                                            slot = "key"), i)
    shifted.rot=scale(shifted.rot)
    if (constrain.var) shifted.rot=sweep(shifted.rot,MARGIN = 2,STATS = sqrt(prior.var),FUN = "*")
    object.all <- SetDimReduction(object.all, reduction.type = paste0(reduction.type, ".shifted"), 
                                  slot = "rotation", new.data = (shifted.rot))
  }
  
  object.all <- SetDimReduction(object.all, reduction.type = paste0(reduction.type, ".shifted"), 
                                slot = "key", new.data = 
                                  paste0("S", GetDimReduction(object, reduction.type = reduction.type, 
                                                              slot = "key")))
  return(object.all)
}
# Density dependent downsampling

DDDS <- function(data, pct = 0.5, amt = 300){
  data.approx <- approxfun(density(data))(data)
  names(data.approx) <- rownames(data)
  pct <- quantile(data.approx, pct)
  data.downsampled <- sample(data.approx, size = amt, prob = 1/(data.approx + pct))
  return(data.downsampled)
}

# Shifting objective function 

ShiftObj <- function(par, drrot, ids, cells.1, fun.opt = 2) {
  if(fun.opt == 1) drrot[cells.1, ] <- par[1] * drrot[cells.1, ]
  else if(fun.opt == 2) drrot[cells.1, ] <- par[2]*drrot[cells.1, ] + par[1]
  else stop("invalid objective function option")
  drrot <- drrot[order(drrot),,drop=F]
  drident <- as.numeric(ids[rownames(drrot)])
  drmax <- max(rle(drident)$lengths,0.975)
  return(drmax)
}
