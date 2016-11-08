#' @include seurat.R
NULL

#' Dimensional Reduction
#' 
#' Various methods for dimensional reductions
#' 
#' @param object Seurat object
#' @param reduction.type Type of dimensional reduction to run. Options include "pca", "pcafast", 
#' "ica", "icafast" 
#' @param genes.use Genes to use as input for the dimensional reduction technique. Default is 
#' object@@var.genes.
#' @param dims.store Number of dimensions to store
#' @param dims.compute Number of dimensions to compute (for fast approximations)
#' @param use.imputed Whether to run the dimensional reduction on imputed values.
#' @param rev.reduction By default, computes the dimensional reduction on the cell x gene matrix.
#' Setting to true will compute it on the transpose (gene x cell matrix).
#' @param print.results Print the top genes associated with each dimension
#' @param dims.print Number of dimensions to print genes for
#' @param genes.print Number of genes to print for each PC
#' @param ... Additional arguments to be passed to specific reduction technique
#' @return Returns a Seurat object with the dimensional reduction information stored 
#' @export
DimReduction <- function(object, reduction.type = NULL, genes.use = NULL, dims.store = 40, 
                         dims.compute = 20, use.imputed = FALSE, rev.reduction = FALSE, print.results = TRUE, 
                         dims.print = 5, genes.print = 30, ...){
  
  if (length(object@scale.data) == 0){
    stop("Object@scale.data has not been set. Run ScaleData() and then retry.")
  }
  if (length(object@var.genes) == 0 && is.null(pc.genes)) {
    stop("Variable genes haven't been set. Run MeanVarPlot() or provide a vector of genes names in 
         pc.genes and retry.")
  }
  browser()
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
    pca <- RunPCA(data.use = data.use, rev.pca = rev.reduction, pcs.store = dims.store, ...)
    object@dim.reduction$pca <- pca
  }
  if(reduction.type == "pcafast") {
    pcafast <- RunPCAFast(data.use = data.use, rev.pca = rev.reduction, pcs.store = dims.store, 
                      pcs.compute = dims.compute, ...)
    object@dim.reduction$pcafast <- pcafast
  }
  if(reduction.type == "ica") print("doit")
  if(reduction.type == "icafast") print("doit")
  
  # print results
  if(print.results){
    results <- eval(parse(text = paste0("object@dim.reduction$", reduction.type, "$x")))
    for(i in 1:dims.print) {
      genes.ordered <- results[order(results[, i]), ]
      top.genes <- genes.ordered[1:genes.print, ]
      print(colnames(top.genes)[i])
      print(rownames(top.genes))
      print("")
      
      top.genes <- genes.ordered[(nrow(genes.ordered) - genes.print + 1) : nrow(genes.ordered), ]
      print(rev(rownames(top.genes)))
      print("")
      print("")
    }
  }
  return(object)
} 

RunPCA <- function(data.use, rev.pca, pcs.store, ...){
  pcs.store <- min(pcs.store, ncol(data.use))
  if(rev.pca){
    pca.obj <- prcomp(data.use, ...)
    pca.obj$x <- pca.obj$x[, 1:pcs.store]
    pca.obj$rotation <- pca.obj$rotation[, 1:pcs.store]
  }
  else{
    pca.obj = prcomp(t(data.use), ...)
    pca.obj$x <- pca.obj$rotation[, 1:pcs.store]
    pca.obj$rotation <- pca.obj$x[, 1:pcs.store]
  }
  pca.obj$center <- NULL
  pca.obj$scale <- NULL
  return(pca.obj)  
}

RunPCAFast <- function(data.use, rev.pca, pcs.store, pcs.compute, ...){
  pcs.compute <- min(pcs.compute, ncol(data.use))
  pcs.store <- min(pcs.store, pcs.compute)
  if(rev.pca){
    pca.obj <- irlba(data.use, nv = pcs.compute, ...)
    pca.obj$x <- pca.obj$u[, 1:pcs.store]
    pca.obj$rotation <- pca.obj$v[, 1:pcs.store]
  }
  else{
    pca.obj <- irlba(t(data.use), nv = pcs.compute, ...)
    pca.obj$x <- pca.obj$v[, 1:pcs.store]
    pca.obj$rotation <- pca.obj$u[, 1:pcs.store]

  }
  rownames(pca.obj$x) <- rownames(data.use)
  colnames(pca.obj$x) <- paste0("PC", 1:pcs.compute)
  rownames(pca.obj$rotation) <- colnames(data.use)
  colnames(pca.obj$rotation) <- colnames(pca.obj$x)
  return(pca.obj)
}
