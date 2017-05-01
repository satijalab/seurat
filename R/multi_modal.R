#' @include seurat.R
NULL

# Set up assay class to hold multimodal data sets

assay <- setClass("assay", slots = list(
  raw.data = "ANY", data = "ANY", scale.data = "ANY", key = "character", 
  misc = "ANY",var.genes="vector",mean.var="data.frame"
))

#' Accessor function for multimodal data
#' 
#' Pull information for specified stored dimensional reduction analysis
#' 
#' @param object Seurat object
#' @param assay.type Type of assay to fetch data for (default is RNA)
#' @param slot Specific information to pull (i.e. raw.data, data, scale.data,...). Default is data
#' @return Returns assay data
#' @export
GetAssayData <- function(object, assay.type = "RNA", slot = "data") {
  if (assay.type=="RNA") {
    if (slot=="raw.data") {
      to.return=object@raw.data
    }
    if (slot=="data") {
      to.return=object@data
    }    
    if (slot=="scale.data") {
      to.return=object@scale.data
    }
    
    #note that we check for this to avoid a long subset for large matrices if it can be avoided
    if (length(object@cell.names)==ncol(to.return)) return(to.return)
    return(to.return[,object@cell.names])
    
  }
  if (!(assay.type %in% names(object@assay))) {
    stop(paste(assay.type, " data has not been added"))
  }
  if (!(slot %in% slotNames(eval(parse(text = paste0("object@assay$", assay.type)))))) {
    stop(paste0(slot, " slot doesn't exist"))
  }
  to.return=(eval(parse(text = paste0("object@assay$", assay.type, "@", slot))))
  if (length(object@cell.names)==ncol(to.return)) return(to.return)
  return(to.return[,object@cell.names])
}

#' Assay Data Mutator Function
#' 
#' Store information for specified assay, for multimodal analysis
#' 
#' @inheritParams GetAssayData
#' @param new.data New data to insert
#' @return Seurat object with updated slot
#' @export
SetAssayData <- function(object, assay.type, slot, new.data) {
  if (assay.type=="RNA") {
    if (slot=="raw.data") {
      (object@raw.data=new.data)
    }
    if (slot=="data") {
      (object@data=new.data)
    }    
    if (slot=="scale.data") {
      (object@scale.data=new.data)
    }
    return(object)
  }
  if (assay.type %in% names(object@assay)) {
    eval(parse(text = paste0("object@assay$", assay.type, "@", slot, "<- new.data")))
  }
  else{
    new.assay <- new("assay")
    eval(parse(text = paste0("new.assay@", slot, "<- new.data")))
    eval(parse(text = paste0("object@assay$", assay.type, "<- new.assay")))
  }
  return(object)
}



#' Normalize Assay Data
#' 
#' Normalize data for a given assay
#' 
#' @param object Seurat object
#' @param assay.type Type of assay to normalize for (default is RNA), but can be changed for multimodal analyses.
#' @param normalization.method Method for normalization. Default is log-normalization (LogNormalize). Other options include CLR (CLR), regularized NB normalization (NBReg; RNA only)
#' @importFrom compositions clr
#' @return Returns object after normalization. Normalized data is stored in data or scale.data slot, depending on the method
#' @export
NormalizeData <- function(object, assay.type = "RNA", normalization.method = "LogNormalize", slot = "data", scale.factor=1e4, display.progress=T, ...) {

  if (normalization.method == "LogNormalize") {
    raw.data=GetAssayData(object,assay.type,"raw.data")
    if (is.null(raw.data))  stop(paste0("Raw data for ", assay.type, " has not been set"))
    normalized.data=LogNormalize(raw.data,scale.factor,display.progress)
    object=SetAssayData(object,assay.type,"data",normalized.data)
  }
  if (normalization.method == "CLR") {
    raw.data=GetAssayData(object,assay.type,"raw.data")
    normalized.data=t(as.matrix(t(clr(raw.data))))
    object=SetAssayData(object,assay.type,"data",normalized.data)
    object=SetAssayData(object,assay.type,"scale.data",normalized.data)
  }
  return(object)
}



#' Run Canonical Correlation Analysis (CCA) on multimodal data
#' 
#' CCA finds a shared correlation structure betwen two different datasets, enabling integrated downstream analysis
#' 
#' @param object Seurat object
#' @param assay.1 First assay for multimodal analysis. Default is RNA
#' @param assay.2 Second assay for multimodal analysis. Default is CITE for CITE-Seq analysis.
#' @param features.1 Features of assay 1 to consider (default is variable genes)
#' @param features.2 Features of assay 2 to consider (default is all features, i.e. for CITE-Seq, all antibodies)
#' @param num.cc Number of canonical correlations to compute and store. Default is 20, but will calculate less if either assay has <20 features.
#' @param normalize.variance Z-score the embedding of each CC to 1, so each CC contributes equally in downstream analysis (default is T)
#' @importFrom PMA CCA
#' @return Returns object after CCA, with results stored in dimensional reduction cca.assay1 (ie. cca.RNA) and cca.assay2. For example, results can be visualized using DimPlot(object,reduction.use="cca.RNA")
#' @export
MultiModal_CCA <- function(object,assay.1="RNA",assay.2="CITE",features.1=NULL,features.2=NULL,num.cc=20,normalize.variance=T) {
  
  #first pull out data, define features
  data.1=GetAssayData(object,assay.1,"scale.data")
  data.2=GetAssayData(object,assay.2,"scale.data")
  
  if (is.null(features.1)) {
    if ((assay.1=="RNA") && length(object@var.genes)>0) {
      features.1=object@var.genes
    } else {
      features.1=rownames(data.1)
    }
  }
  
  features.2=set.ifnull(features.2,rownames(data.2))
  
  data.1=t(data.1[features.1,])
  data.2=t(data.2[features.2,])
  num.cc=min(20,min(length(features.1),length(features.2)))
  cca.data=list(data.1,data.2)
  names(cca.data)=c(assay.1,assay.2)
  # now run CCA
  out <- CCA(cca.data[[1]],cca.data[[2]],typex="standard",typez="standard",K=num.cc,penaltyz = 1,penaltyx = 1)
  cca.output=list(out$u,out$v); 
  embeddings.cca=list()
  
  for(i in 1:length(cca.data)) {
    assay.use=names(cca.data)[i]
    rownames(cca.output[[i]])=colnames(cca.data[[i]])
    embeddings.cca[[i]]=cca.data[[i]]%*%cca.output[[i]]
    colnames(embeddings.cca[[i]])=paste0(assay.use,"CC",1:ncol( embeddings.cca[[i]]),sep="")
    colnames(cca.output[[i]])=colnames(embeddings.cca[[i]])
    if (normalize.variance) {
      embeddings.cca[[i]]=scale(embeddings.cca[[i]])
    }
    
    object=SetDimReduction(object,reduction.type = paste0(assay.use,"CCA",sep=""),slot = "rotation",new.data = embeddings.cca[[i]])
    object=SetDimReduction(object,reduction.type = paste0(assay.use,"CCA",sep=""),slot = "key",new.data =  paste0(assay.use,"CC",sep=""))
    object=SetDimReduction(object,reduction.type = paste0(assay.use,"CCA",sep=""),slot = "x",new.data =  cca.output[[i]])
  }
  return(object)
}


#' Run coinertia analysis on multimodal data
#' 
#' CIA finds a shared correlation structure betwen two different datasets, enabling integrated downstream analysis
#' 
#' @param object Seurat object
#' @param assay.1 First assay for multimodal analysis. Default is RNA
#' @param assay.2 Second assay for multimodal analysis. Default is CITE for CITE-Seq analysis.
#' @param features.1 Features of assay 1 to consider (default is variable genes)
#' @param features.2 Features of assay 2 to consider (default is all features, i.e. for CITE-Seq, all antibodies)
#' @param num.axes Number of principal axes to compute and store. Default is 20, but will calculate less if either assay has <20 features.
#' @param normalize.variance Return the normalized row scares, so each aexis contributes equally in downstream analysis (default is T)
#' @importFrom made4 cia
#' @return Returns object after CIA, with results stored in dimensional reduction cia.assay1 (ie. cia.RNA) and cia.assay2. For example, results can be visualized using DimPlot(object,reduction.use="cia.RNA")
#' @export
MultiModal_CIA <- function(object,assay.1="RNA",assay.2="CITE",features.1=NULL,features.2=NULL,num.axes=20,normalize.variance=T) {
  
  #first pull out data, define features
  data.1=GetAssayData(object,assay.1,"scale.data")
  data.2=GetAssayData(object,assay.2,"scale.data")
  
  if (is.null(features.1)) {
    if ((assay.1=="RNA") && length(object@var.genes)>0) {
      features.1=object@var.genes
    } else {
      features.1=rownames(data.1)
    }
  }
  
  features.2=set.ifnull(features.2,rownames(data.2))
  
  data.1=t(data.1[features.1,])
  data.2=t(data.2[features.2,])
  num.axes=min(20,min(length(features.1),length(features.2)))
  cia.data=list(data.1,data.2)
  names(cia.data)=c(assay.1,assay.2)
  # now run cia
  out=cia(t(cia.data[[1]]),t(cia.data[[2]]),cia.nf = num.axes)
  out=out$coinertia
  cia.output=list(as.matrix(out$c1),as.matrix(out$l1))
  embeddings.cia.norm=list(as.matrix(out$mX),as.matrix(out$mY))
  embeddings.cia=list(as.matrix(out$lX),as.matrix(out$lY))
  
  
  for(i in 1:length(cia.data)) {
    assay.use=names(cia.data)[i]
    #rownames(cia.output[[i]])=colnames(cia.data[[i]])
    if (normalize.variance) {
      embeddings.cia[[i]]=(embeddings.cia.norm[[i]])
    }
    colnames(embeddings.cia[[i]])=paste0(assay.use,"CI",1:ncol( embeddings.cia[[i]]),sep="")
    colnames(cia.output[[i]])=colnames(embeddings.cia[[i]])

    
    object=SetDimReduction(object,reduction.type = paste("cia",assay.use, sep="_"),slot = "rotation",new.data = embeddings.cia[[i]])
    object=SetDimReduction(object,reduction.type = paste("cia",assay.use, sep="_"),slot = "key",new.data =  paste0(assay.use,"CI",sep=""))
    object=SetDimReduction(object,reduction.type = paste("cia",assay.use, sep="_"),slot = "x",new.data =  cia.output[[i]])
  }
  return(object)
}

