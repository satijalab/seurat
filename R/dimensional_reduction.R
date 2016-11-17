#' @include seurat.R
NULL

# Set up dim.reduction class

dim.reduction <- setClass("dim.reduction", slots = list(
  rotation = "matrix", x = "matrix", x.full = "matrix", sdev = "numeric", key = "character", 
  misc = "ANY"
))

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
                         dims.compute = 40, use.imputed = FALSE, rev.reduction = FALSE, 
                         print.results = TRUE, dims.print = 1:5, genes.print = 30, ica.fxn = icafast,
                         ...){
  
  if (length(object@scale.data) == 0){
    stop("Object@scale.data has not been set. Run ScaleData() and then retry.")
  }
  if (length(object@var.genes) == 0 && is.null(genes.use)) {
    stop("Variable genes haven't been set. Run MeanVarPlot() or provide a vector of genes names in 
         genes.use and retry.")
  }
  
  dims.store=min(dims.store,dims.compute)
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
  }
  if(reduction.type == "ica") {
    icaobj=RunICA(data.use = data.use, ics.compute = dims.store, rev.ica = rev.reduction, 
                  ics.store = dims.store, ica.fxn = ica.fxn,...)
    object@dr$ica=icaobj
  }
  
  #if(reduction.type == "icafast") print("doit")
  
  if(print.results){
    PrintDim(object, reduction.type = reduction.type, dims.print = 1:dims.print,
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
    pca.results = prcomp(t(data.use), ...)
    x <- pca.results$rotation[, 1:pcs.store]
    rotation <- pca.results$x[, 1:pcs.store]
  }
  pca.obj <- new("dim.reduction", x = x, rotation = rotation, sdev = pca.results$sdev, key = "PC")
  return(pca.obj)  
}

RunICA <- function(data.use, ics.compute, rev.ica, ica.fxn=icafast, ics.store,...) {
  ics.store <- min(ics.store, ncol(data.use))
  ica.results <- NULL
  if(rev.ica){
    ica.results <- ica.fxn(data.use, ics.compute,...)
    rotation <- ica.results$M[,1:ics.store]
  }
  else{
    ica.results = ica.fxn(t(data.use),ics.compute,...)
    rotation <- ica.results$S[,1:ics.store]
  }
  
  x=(as.matrix(data.use)%*%as.matrix(rotation))
  colnames(x)=paste("IC",1:ncol(x),sep="")
  colnames(rotation)=paste("IC",1:ncol(x),sep="")
  
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
#' @param rev.pca By default computes the PCA on the cell x gene matrix. Setting to true will compute it on gene x cell matrix. 
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
#' @param pcs.print Number of PCs to print genes for
#' @param pcs.store Number of PCs to store
#' @param pcs.compute Total Number of PCs to compute and store
#' @param genes.print Number of genes to print for each PC
#' @param \dots Additional arguments to be passed to prcomp
#' @return Returns Seurat object with an PCA embedding (object@@pca.rot) and
#' gene projection matrix (object@@pca.x). The PCA object itself is stored in
#' object@@pca.obj[[1]]
#' @importFrom irlba irlba
#' @export
PCAFast <- function(object, pc.genes = NULL, do.print = TRUE, pcs.print = 5, pcs.store = 40, 
                    pcs.compute = 20, genes.print = 30, ...) {
  return(DimReduction(object, reduction.type = "pcafast", genes.use = pc.genes, 
                      print.results = do.print, dims.store = pcs.store, dims.compute = pcs.compute, 
                      genes.print = genes.print))
}


#' Run Independent Component Analysis on gene expression
#'
#' Run fastICA algorithm for ICA dimensionality reduction
#'
#'
#' @param object Seurat object
#' @param ic.genes Genes to use as input for ICA. Default is object@@var.genes
#' @param do.print Print the top genes associated with high/low loadings for
#' the ICs
#' @param ics.print Number of ICs to print genes for
#' @param ics.store Number of ICs to store
#' @param genes.print Number of genes to print for each IC
#' @param use.imputed Run ICA on imputed values (FALSE by default)
#' @param seed.use Random seed to use for fastICA
#' @param \dots Additional arguments to be passed to fastICA
#' @return Returns Seurat object with an ICA embedding (object@@ica.rot) and
#' gene projection matrix (object@@ica.x). The ICA object itself is stored in
#' object@@ica.obj[[1]]
#' @import fastICA
#' @export
ICA <- function(object, ic.genes = NULL, do.print = TRUE, ics.print = 5, ics.store = 50, 
                genes.print = 50, use.imputed = FALSE, seed.use = 1, ...) {
  return(DimReduction(object, reduction.type = "ica", genes.use = ic.genes, dims.store = ics.store, 
                      use.imputed = use.imputed, print.results = do.print, dims.print = ics.print, 
                      genes.print = genes.print, ... ))
}


#Internal, not documented for now
topGenesForDim=function(i,dim_scores,do.balanced=FALSE,num.genes=30,reduction.use="pca") {
  code=paste(translate.dim.code(reduction.use),i,sep="")
  if (do.balanced) {
    num.genes=round(num.genes/2)
    sx=dim_scores[order(dim_scores[,code]),]
    genes.1=(rownames(sx[1:num.genes,]))
    genes.2=(rownames(sx[(nrow(sx)-num.genes+1):nrow(sx),]))
    return(c(genes.1,genes.2))
  }
  if (!(do.balanced)) {
    sx=dim_scores[rev(order(abs(dim_scores[,code]))),]
    genes.1=(rownames(sx[1:num.genes,]))
    genes.1=genes.1[order(dim_scores[genes.1,code])]
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
  num.cells=set.ifnull(num.cells,length(object@cell.names))
  dim_scores=eval(parse(text=paste("object@dr$",reduction.type,"@rotation",sep="")))
  i=dim.use
  dim.top.cells=unique(unlist(lapply(i,topGenesForDim,dim_scores,do.balanced,num.cells,reduction.type)))
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
  dim_scores=eval(parse(text=paste("object@dr$",reduction.type,"@x",sep="")))
  if (use.full) dim_scores=eval(parse(text=paste("object@dr$",reduction.type,"@x.full",sep="")))
  i=dim.use
  dim.top.genes=unique(unlist(lapply(i,topGenesForDim,dim_scores,do.balanced,num.genes,reduction.type)))
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
    dim_scores <- eval(parse(text = paste("object@dr$", reduction.type, "@rotation", sep="")))
    dim_key <- eval(parse(text = paste("object@dr$", reduction.type, "@key", sep="")))
    cells.ordered <- cells.use[order(dim_scores[cells.use, paste(dim_key, ndim, sep="")])]
    data.use <- object@scale.data[genes.use, cells.ordered]
    data.use <- minmax(data.use, min = disp.min, max = disp.max)
    if (!(use.scale)) data.use <- as.matrix(object@data[genes.use, cells.ordered])
    vline.use <- NULL
    hmTitle <- paste(dim_key,ndim)
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
    eval(parse(text=hmFunction))
  }
  if (do.return) {
    return(data.use)
  }
  # reset graphics parameters
  par(mfrow=orig_par)
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
  
  if(length(eval(parse(text = paste0("object@dr$", reduction.type, "@x.full")))) == 0 && use.full){
    warning("Dimensions have not been projected. Setting use.full = FALSE")
    use.full <- FALSE
  }
  for(i in dims.print) {
    code <- paste0(eval(parse(text = paste0("object@dr$", reduction.type, "@key"))), i)
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

