

################################################################################
### Seurat

#' The Seurat Class
#'
#' The Seurat object is the center of each single cell analysis. It stores all information 
#' associated with the dataset, including data, annotations, analyes, etc. All that is needed 
#' to construct a Seurat object is an expression matrix (rows are genes, columns are cells), which
#' should be log-scale
#' 
#' Each Seurat object has a number of slots which store information. Key slots to access
#' are listed below.
#'
#'
#'@section Slots:
#'  \describe{
#'    \item{\code{data}:}{\code{"data.frame"}, The expression matrix (log-scale) }
#'    \item{\code{scale.data}:}{\code{"data.frame"}, The scaled (after z-scoring 
#'    each gene) expression matrix. Used for PCA, ICA, and heatmap plotting}
#'    \item{\code{var.genes}:}{\code{"vector"},  Variable genes across single cells }
#'    \item{\code{is.expr}:}{\code{"numeric"}, Expression threshold to determine if a gene is expressed }
#'    \item{\code{ident}:}{\code{"vector"},  The 'identity class' for each single cell }
#'    \item{\code{data.info}:}{\code{"data.frame"}, Contains information about each cell, starting with # of genes detected (nGene)
#'    the original identity class (orig.ident), user-provided information (through addMetaData), etc.  }
#'    \item{\code{project.name}:}{\code{"character"}, Name of the project (for record keeping) }
#'    \item{\code{pca.x}:}{\code{"data.frame"}, Gene projection scores for the PCA analysis  }
#'    \item{\code{pca.x.full}:}{\code{"data.frame"}, Gene projection scores for the projected PCA  (contains all genes)  }
#'    \item{\code{pca.rot}:}{\code{"data.frame"}, The rotation matrix (eigenvectors) of the PCA   }
#'    \item{\code{ica.x}:}{\code{"data.frame"}, Gene projection scores for ICA   }
#'    \item{\code{ica.rot}:}{\code{"data.frame"}, The estimated source matrix from ICA }
#'    \item{\code{tsne.rot}:}{\code{"data.frame"}, Cell coordinates on the t-SNE map }
#'    \item{\code{mean.var}:}{\code{"data.frame"}, The output of the mean/variability analysis for all genes }
#'    \item{\code{imputed}:}{\code{"data.frame"}, Matrix of imputed gene scores }
#'    \item{\code{final.prob}:}{\code{"data.frame"}, For spatial inference, posterior probability of each cell mapping to each bin }
#'    \item{\code{insitu.matrix}:}{\code{"data.frame"}, For spatial inference, the discretized spatial reference map }
#'    \item{\code{cell.names}:}{\code{"vector"},  Names of all single cells (column names of the expression matrix) }
#'    \item{\code{cluster.tree}:}{\code{"list"},  List where the first element is a phylo object containing the 
#'    phylogenetic tree relating different identity classes }
#'      
#'}
#' @name seurat
#' @rdname seurat
#' @aliases seurat-class
#' @exportClass seurat
#' @importFrom useful corner

seurat <- setClass("seurat", slots = 
                     c(raw.data = "data.frame", data="data.frame",scale.data="matrix",var.genes="vector",is.expr="numeric",
                       ident="vector",pca.x="data.frame",pca.rot="data.frame",
                       emp.pval="data.frame",kmeans.obj="list",pca.obj="list",
                       gene.scores="data.frame", drop.coefs="data.frame",
                       wt.matrix="data.frame", drop.wt.matrix="data.frame",trusted.genes="vector",drop.expr="numeric",data.info="data.frame",
                       project.name="character", kmeans.gene="list", kmeans.cell="list",jackStraw.empP="data.frame", 
                       jackStraw.fakePC = "data.frame",jackStraw.empP.full="data.frame",pca.x.full="data.frame", kmeans.col="list",mean.var="data.frame", imputed="data.frame",mix.probs="data.frame",
                       mix.param="data.frame",final.prob="data.frame",insitu.matrix="data.frame",
                       tsne.rot="data.frame", ica.rot="data.frame", ica.x="data.frame", ica.obj="list",cell.names="vector",cluster.tree="list"))

calc.drop.prob=function(x,a,b) {
  return(exp(a+b*x)/(1+exp(a+b*x)))
}


setGeneric("find_all_markers_node", function(object, thresh.test=1,test.use="bimod",return.thresh=1e-2,do.print=FALSE) standardGeneric("find_all_markers_node"))
setMethod("find_all_markers_node","seurat",
          function(object, thresh.test=1,test.use="bimod",return.thresh=1e-2,do.print=FALSE) {
            ident.use=object@ident
            tree.use=object@cluster.tree[[1]]
            
            if ((test.use=="roc") && (return.thresh==1e-2)) return.thresh=0.8
            genes.de=list()
            for(i in ((tree.use$Nnode+2):max(tree.use$edge))) {
              genes.de[[i]]=find.markers.node(object,i,genes.use=rownames(object@data),thresh.use = thresh.test,test.use = test.use)
              if (do.print) print(paste("Calculating node", i))
            }
            gde.all=data.frame()
            for(i in ((tree.use$Nnode+2):max(tree.use$edge))) {
              gde=genes.de[[i]]
              if (is.null(gde)) next;
              if (nrow(gde)>0) {
                if (test.use=="roc") gde=subset(gde,(myAUC>return.thresh|myAUC<(1-return.thresh)))
                if (test.use=="bimod") {
                  gde=gde[order(gde$p_val,-gde$avg_diff),]
                  gde=subset(gde,p_val<return.thresh)
                }
                if (nrow(gde)>0) gde$cluster=i; gde$gene=rownames(gde)
                if (nrow(gde)>0) gde.all=rbind(gde.all,gde)
              }
            }
            return(gde.all)
          }
)

weighted.euclidean=function(x,y,w) {
  v.dist=sum(sqrt(w*(x-y)^2))
  return(v.dist)
}

#from Jean Fan - thanks!!
custom.dist <- function(my.mat, my.function,...) {
  n <- ncol(my.mat)
  mat <- matrix(0, ncol = n, nrow = n)
  colnames(mat) <- rownames(mat) <- colnames(my.mat)
  for(i in 1:nrow(mat)) {
    for(j in 1:ncol(mat)) {
      mat[i,j] <- my.function(my.mat[,i],my.mat[,j],...)
    }}
  return(as.dist(mat))
}


#' Phylogenetic Analysis of Identity Classes
#'
#' Constructs a phylogenetic tree relating the 'average' cell from each
#' identity class. Tree is estimated based on a distance matrix constructed in
#' either gene expression space or PCA space.
#'
#' Note that the tree is calculated for an 'average' cell, so gene expression
#' or PC scores are averaged across all cells in an identity class before the
#' tree is constructed.
#' 
#' @param object Seurat object
#' @param genes.use Genes to use for the analysis. Default is the set of
#' variable genes (object@@var.genes). Assumes pcs.use=NULL (tree calculated in
#' gene expression space)
#' @param pcs.use If set, tree is calculated in PCA space, using the
#' eigenvalue-weighted euclidean distance across these PC scores.
#' @param do.plot Plot the resulting phylogenetic tree
#' @param do.reorder Re-order identity classes (factor ordering), according to
#' position on the tree. This groups similar classes together which can be
#' helpful, for example, when drawing violin plots.
#' @param reorder.numeric Re-order identity classes according to position on
#' the tree, assigning a numeric value ('1' is the leftmost node)
#' @return A Seurat object where the cluster tree is stored in
#' object@@cluster.tree[[1]]
#' @importFrom ape as.phylo
#' @export
setGeneric("buildClusterTree", function(object, genes.use=NULL,pcs.use=NULL,do.plot=TRUE,do.reorder=FALSE,reorder.numeric=FALSE) standardGeneric("buildClusterTree"))
#' @export
setMethod("buildClusterTree","seurat",
          function(object,genes.use=NULL,pcs.use=NULL,do.plot=TRUE,do.reorder=FALSE,reorder.numeric=FALSE) {
            genes.use=set.ifnull(genes.use,object@var.genes)
            ident.names=as.character(unique(object@ident))
            if (is.null(pcs.use)) {
              genes.use=ainb(genes.use,rownames(object@data))
              data.avg=average.expression(object,genes.use = genes.use)
              data.dist=dist(t(data.avg[genes.use,]))
            }
            if (!is.null(pcs.use)) {
              data.pca=average.pca(object)
              data.use=data.pca[pcs.use,]
              data.eigenval=(object@pca.obj[[1]]$sdev)^2
              data.weights=(data.eigenval/sum(data.eigenval))[pcs.use]; data.weights=data.weights/sum(data.weights)
              data.dist=custom.dist(data.pca[pcs.use,],weighted.euclidean,data.weights)
            }
            data.tree=as.phylo(hclust(data.dist))
            object@cluster.tree[[1]]=data.tree
            
            if (do.reorder) {
              old.ident.order=sort(unique(object@ident))
              data.tree=object@cluster.tree[[1]]
              all.desc=getDescendants(data.tree,data.tree$Nnode+2); all.desc=old.ident.order[all.desc[all.desc<=(data.tree$Nnode+1)]]
              object@ident=factor(object@ident,levels = all.desc,ordered = TRUE) 
              if (reorder.numeric) {
                object=set.ident(object,object@cell.names,as.integer(object@ident))
                object@data.info[object@cell.names,"tree.ident"]=as.integer(object@ident)
              }
              object=buildClusterTree(object,genes.use,pcs.use,do.plot=FALSE,do.reorder=FALSE)
            }
            if (do.plot) plotClusterTree(object)
            return(object)
          }
)

#' Plot phylogenetic tree
#'
#' Plots previously computed phylogenetic tree (from buildClusterTree)
#'
#' @param object Seurat object
#' @return Plots dendogram (must be precomputed using buildClusterTree), returns no value
#' @importFrom ape plot.phylo
#' @importFrom ape nodelabels
#' @export
setGeneric("plotClusterTree", function(object) standardGeneric("plotClusterTree"))
#' @export
setMethod("plotClusterTree","seurat",
          function(object) {
            if (length(object@cluster.tree)==0) stop("Phylogenetic tree does not exist, build using buildClusterTree")
            data.tree=object@cluster.tree[[1]]
            plot.phylo(data.tree,direction="downwards")
            nodelabels()
          }
)


#' Visualize expression/dropout curve
#'
#' Plot the probability of detection vs average expression of a gene.
#'
#' Assumes that this 'noise' model has been precomputed with calcNoiseModels
#'
#'
#' @param object Seurat object
#' @param cell.ids Cells to use
#' @param col.use Color code or name
#' @param lwd.use Line width for curve
#' @param do.new Create a new plot (default) or add to existing
#' @param x.lim Maximum value for X axis
#' @param \dots Additional arguments to pass to lines function
#' @return Returns no value, displays a plot
#' @export
setGeneric("plotNoiseModel", function(object, cell.ids=c(1,2), col.use="black",lwd.use=2,do.new=TRUE,x.lim=10,...) standardGeneric("plotNoiseModel"))
#' @export
setMethod("plotNoiseModel","seurat",
          function(object, cell.ids=c(1,2), col.use="black",lwd.use=2,do.new=TRUE,x.lim=10,...) {
            cell.coefs=object@drop.coefs[cell.ids,]
            if (do.new) plot(1,1,pch=16,type="n",xlab="Average expression",ylab="Probability of detection",xlim=c(0,x.lim),ylim=c(0,1))
            unlist(lapply(1:length(cell.ids), function(y) {
              x.vals=seq(0,x.lim,0.05)
              y.vals=unlist(lapply(x.vals,calc.drop.prob,cell.coefs[y,1],cell.coefs[y,2]))
              lines(x.vals,y.vals,lwd=lwd.use,col=col.use,...)
            }))
          }
)

#' Setup Seurat object
#'
#' Setup and initialize basic parameters of the Seurat object
#'
#'
#' @param object Seurat object
#' @param project Project name (string)
#' @param min.cells Include genes with detected expression in at least this
#' many cells
#' @param min.genes Include cells where at least this many genes are detected
#' @param is.expr Expression threshold for 'detected' gene
#' @param do.scale In object@@scale.data, perform row-scaling (gene-based
#' z-score)
#' @param do.center In object@@scale.data, perform row-centering (gene-based
#' centering)
#' @param names.field For the initial identity class for each cell, choose this
#' field from the cell's column name
#' @param names.delim For the initial identity class for each cell, choose this
#' delimiter from the cell's column name
#' @param meta.data Additional metadata to add to the Seurat object. Should be
#' a data frame where the rows are cell names, and the columns are additional
#' metadata fields
#' @param save.raw TRUE by default. If FALSE, do not save the unmodified data in object@@raw.data 
#' which will save memory downstream for large datasets
#' @param \dots Additional arguments, not currently used.
#' @return Seurat object. Fields modified include object@@data,
#' object@@scale.data, object@@data.info, object@@ident
#' @import stringr
#' @export
setGeneric("setup", function(object, project, min.cells=3, min.genes=2500, is.expr=1, do.scale=TRUE, do.center=TRUE,names.field=1,names.delim="_",meta.data=NULL,save.raw=TRUE,...) standardGeneric("setup"))
#' @export
setMethod("setup","seurat",
          function(object, project, min.cells=3, min.genes=2500, is.expr=1, do.scale=TRUE, do.center=TRUE,names.field=1,names.delim="_",meta.data=NULL,save.raw=TRUE,...) {
            object@is.expr = is.expr
            num.genes=findNGene(object@raw.data,object@is.expr)
            cells.use=names(num.genes[which(num.genes>min.genes)]) 
            
            object@data=object@raw.data[,cells.use]
            
            #to save memory downstream, especially for large object
            if (!(save.raw)) object@raw.data=data.frame();
            
            num.cells=apply(object@data,1,humpCt,min=object@is.expr)
            genes.use=names(num.cells[which(num.cells>min.cells)])
            object@data=object@data[genes.use,]
            
            object@ident=factor(unlist(lapply(colnames(object@data),extract.field,names.field,names.delim)))
            names(object@ident)=colnames(object@data)
            object@cell.names=names(object@ident)
            object@scale.data=t(scale(t(object@data),center=do.center,scale=do.scale))
            data.ngene=num.genes[cells.use]
            object@gene.scores=data.frame(data.ngene); colnames(object@gene.scores)[1]="nGene"
            object@data.info=data.frame(data.ngene); colnames(object@data.info)[1]="nGene"
            if (!is.null(meta.data)) {
              object=addMetaData(object,metadata = meta.data)
            }
            object@mix.probs=data.frame(data.ngene); colnames(object@mix.probs)[1]="nGene"
            rownames(object@gene.scores)=colnames(object@data)
            
            object@data.info[names(object@ident),"orig.ident"]=object@ident
            
            object@project.name=project
            #if(calc.noise) {
            #  object=calcNoiseModels(object,...)
            #  object=getWeightMatrix(object)
            #}
            return(object)
          }         
)

#' Return a subset of the Seurat object
#'
#' Creates a Seurat object containing only a subset of the cells in the
#' original object. Takes either a list of cells to use as a subset, or a
#' parameter (for example, a gene), to subset on.
#'
#' @param object Seurat object
#' @param cells.use A vector of cell names to use as a subset. If NULL
#' (default), then this list will be computed based on the next three
#' arguments. Otherwise, will return an object consissting only of these cells
#' @param subset.name Parameter to subset on. Eg, the name of a gene, PC1, a
#' column name in object@@data.info, etc. Any argument that can be retreived
#' using fetch.data
#' @param accept.low Low cutoff for the parameter (default is -Inf)
#' @param accept.high High cutoff for the parameter (default is Inf)
#' @param do.center Recenter the new object@@scale.data
#' @param do.scale Rescale the new object@@scale.data
#' @param \dots Additional arguments to be passed to fetch.data (for example,
#' use.imputed=TRUE)
#' @return Returns a Seurat object containing only the relevant subset of cells
#' @export
setGeneric("subsetData",  function(object,cells.use=NULL,subset.name=NULL,accept.low=-Inf, accept.high=Inf,do.center=TRUE,do.scale=TRUE,...) standardGeneric("subsetData"))
#' @export
setMethod("subsetData","seurat",
          function(object,cells.use=NULL,subset.name=NULL,accept.low=-Inf, accept.high=Inf,do.center=TRUE,do.scale=TRUE,...) {
            data.use=NULL
            if (is.null(cells.use)) {
              data.use=fetch.data(object,subset.name,...)
              if (length(data.use)==0) return(object)
              subset.data=data.use[,subset.name]
              pass.inds=which((subset.data>accept.low) & (subset.data<accept.high))
              cells.use=rownames(data.use)[pass.inds]
            }
            object@data=object@data[,cells.use]
            object@scale.data=t(scale(t(object@data),center=do.center,scale=do.scale))
            object@scale.data=object@scale.data[complete.cases(object@scale.data),]
            object@ident=drop.levels(object@ident[cells.use])
            object@tsne.rot=object@tsne.rot[cells.use,]
            object@pca.rot=object@pca.rot[cells.use,]
            object@cell.names=cells.use
            
            object@gene.scores=data.frame(object@gene.scores[cells.use,]); colnames(object@gene.scores)[1]="nGene"; rownames(object@gene.scores)=colnames(object@data)
            object@data.info=data.frame(object@data.info[cells.use,])
            object@mix.probs=data.frame(object@mix.probs[cells.use,]); colnames(object@mix.probs)[1]="nGene"; rownames(object@mix.probs)=colnames(object@data)
            
            return(object)
          }         
)

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
setGeneric("project.pca", function(object,do.print=TRUE,pcs.print=5,pcs.store=30,genes.print=30,replace.pc=FALSE,do.center=TRUE) standardGeneric("project.pca"))
#' @export
setMethod("project.pca", "seurat", 
          function(object,do.print=TRUE,pcs.print=5,pcs.store=30,genes.print=30,replace.pc=FALSE,do.center=TRUE) {
            if (!(do.center)) object@pca.x.full=data.frame(as.matrix(object@scale.data)%*%as.matrix(object@pca.rot))
            if (do.center) object@pca.x.full=data.frame(scale(as.matrix(object@scale.data),center = TRUE,scale = FALSE)%*%as.matrix(object@pca.rot))
            if (ncol(object@jackStraw.fakePC)>0) {
              object@jackStraw.empP.full=data.frame(sapply(1:ncol(object@jackStraw.fakePC),function(x)unlist(lapply(abs(object@pca.x.full[,x]),empP,abs(object@jackStraw.fakePC[,x])))))
              colnames(object@jackStraw.empP.full)=paste("PC",1:ncol(object@jackStraw.empP),sep="")
              rownames(object@jackStraw.empP.full)=rownames(object@scale.data)
            }
            
            if (replace.pc==TRUE) {
              object@jackStraw.empP=object@jackStraw.empP.full
              object@pca.x=object@pca.x.full
            }
            
            if (do.print) {
                print.pca(object,1:pcs.print,genes.print,TRUE)  
            }
            return(object)
          }
)

#' Run t-distributed Stochastic Neighbor Embedding
#'
#' Run t-SNE dimensionality reduction on selected features. Has the option of running in a reduced 
#' dimensional space (i.e. spectral tSNE, recommended), or running based on a set of genes
#' 
#'
#' @param object Seurat object
#' @param cells.use Which cells to analyze (default, all cells)
#' @param dims.use Which dimensions to use as input features 
#' @param k.seed Random seed for the t-SNE
#' @param do.fast If TRUE, uses the Barnes-hut implementation, which runs
#' faster, but is less flexible
#' @param add.iter If an existing tSNE has already been computed, uses the 
#' current tSNE to seed the algorithm and then adds additional iterations on top of this
#' @param genes.use If set, run the tSNE on this subset of genes 
#' (instead of running on a set of reduced dimensions). Not set (NULL) by default
#' @param reduction.use Which dimensional reduction (PCA or ICA) to use for the tSNE. Default is PCA
#' @param dim_embed The dimensional space of the resulting tSNE embedding (default is 2).
#' For example, set to 3 for a 3d tSNE
#' @param \dots Additional arguments to the tSNE call. Most commonly used is
#' perplexity (expected number of neighbors default is 30)
#' @return Returns a Seurat object with a tSNE embedding in object@@tsne_rot
#' @import Rtsne
#' @import tsne
#' @export
setGeneric("run_tsne", function(object,cells.use=NULL,dims.use=1:5,k.seed=1,do.fast=FALSE,add.iter=0,genes.use=NULL,reduction.use="pca",dim_embed=2,...) standardGeneric("run_tsne"))
#' @export
setMethod("run_tsne", "seurat", 
          function(object,cells.use=NULL,dims.use=1:5,k.seed=1,do.fast=FALSE,add.iter=0,genes.use=NULL,reduction.use="pca",dim_embed=2,...) {
            cells.use=set.ifnull(cells.use,colnames(object@data))
            if (is.null(genes.use)) {
              dim.code=translate.dim.code(reduction.use); dim.codes=paste(dim.code,dims.use,sep="")
              data.use=fetch.data(object,dim.codes)
            }
            if (!is.null(genes.use)) {
              genes.use=ainb(genes.use,rownames(object@scale.data))
              data.use=t(object@scale.data[genes.use,cells.use])
            }
            
            
            #data.dist=as.dist(mahalanobis.dist(data.use))
            if (do.fast) {
              set.seed(k.seed); data.tsne=Rtsne(as.matrix(data.use),dims=dim_embed,...)
              data.tsne=data.frame(data.tsne$Y)
            }
            if (!(do.fast)) {
              set.seed(k.seed); data.tsne=data.frame(tsne(data.use,k=dim_embed,...))
            }
            if (add.iter > 0) {
              data.tsne=data.frame(tsne(data.use,initial_config = as.matrix(data.tsne),max_iter = add.iter,...))
            }
            colnames(data.tsne)=paste("tSNE_",1:ncol(data.tsne),sep="")
            rownames(data.tsne)=cells.use
            object@tsne.rot=data.tsne
            return(object)
          }
)

#Not currently supported
setGeneric("add_tsne", function(object,cells.use=NULL,pcs.use=1:10,do.plot=TRUE,k.seed=1,add.iter=1000,...) standardGeneric("add_tsne"))
setMethod("add_tsne", "seurat", 
          function(object,cells.use=NULL,pcs.use=1:10,do.plot=TRUE,k.seed=1,add.iter=1000,...) {
            cells.use=set.ifnull(cells.use,colnames(object@data))
            data.use=object@pca.rot[cells.use,pcs.use]
            #data.dist=as.dist(mahalanobis.dist(data.use))
            set.seed(k.seed); data.tsne=data.frame(tsne(data.use,initial_config = as.matrix(object@tsne.rot[cells.use,]),max_iter = add.iter,...))
            colnames(data.tsne)=paste("tSNE_",1:ncol(data.tsne),sep="")
            rownames(data.tsne)=cells.use
            object@tsne.rot=data.tsne
            return(object)
          }
)

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
setGeneric("ica", function(object,ic.genes=NULL,do.print=TRUE,ics.print=5,ics.store=30,genes.print=30,use.imputed=FALSE,seed.use=1,...) standardGeneric("ica"))
#' @export
setMethod("ica", "seurat", 
          function(object,ic.genes=NULL,do.print=TRUE,ics.print=5,ics.store=30,genes.print=30,use.imputed=FALSE,seed.use=1,...) {            
            data.use=object@scale.data
            if (use.imputed) data.use=data.frame(t(scale(t(object@imputed))))
            ic.genes=set.ifnull(ic.genes,object@var.genes)
            ic.genes = unique(ic.genes[ic.genes%in%rownames(data.use)])
            ic.genes.var = apply(data.use[ic.genes,],1,var)
            ic.data = data.use[ic.genes[ic.genes.var>0],]
            set.seed(seed.use); ica.obj = fastICA(t(ic.data),n.comp=ics.store,...)
            object@ica.obj=list(ica.obj)
            ics.store=min(ics.store,ncol(ic.data))
            ics.print=min(ics.print,ncol(ic.data))
            ic_scores=data.frame(ic.data%*%ica.obj$S)
            colnames(ic_scores)=paste("IC",1:ncol(ic_scores),sep="")
            object@ica.x=ic_scores
            object@ica.rot=data.frame(ica.obj$S[,1:ics.store])
            colnames(object@ica.rot)=paste("IC",1:ncol(object@ica.rot),sep="")
            rownames(object@ica.rot)=colnames(ic.data)
            if (do.print) {
              for(i in 1:ics.print) {
                code=paste("IC",i,sep="")
                sx=ic_scores[order(ic_scores[,code]),]
                print(code)
                print(rownames(sx[1:genes.print,]))
                print ("")
                
                print(rev(rownames(sx[(nrow(sx)-genes.print):nrow(sx),])))
                print ("")
                print ("")      
              }
            }
            return(object)
          }
)

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
#' @param \dots Additional arguments to be passed to prcomp
#' @return Returns Seurat object with an PCA embedding (object@@pca.rot) and
#' gene projection matrix (object@@pca.x). The PCA object itself is stored in
#' object@@pca.obj[[1]]
#' @export
setGeneric("pca", function(object,pc.genes=NULL,do.print=TRUE,pcs.print=5,pcs.store=40,genes.print=30,use.imputed=FALSE,...) standardGeneric("pca"))
#' @export
setMethod("pca", "seurat", 
          function(object,pc.genes=NULL,do.print=TRUE,pcs.print=5,pcs.store=40,genes.print=30,use.imputed=FALSE,...) {
            data.use=object@scale.data
            if (use.imputed) data.use=data.frame(t(scale(t(object@imputed))))
            pc.genes=set.ifnull(pc.genes,object@var.genes)
            pc.genes = unique(pc.genes[pc.genes%in%rownames(data.use)])
            pc.genes.var = apply(data.use[pc.genes,],1,var)
            pc.data = data.use[pc.genes[pc.genes.var>0],]
            pca.obj = prcomp(pc.data,...)
            object@pca.obj=list(pca.obj)
            
              
            pcs.store=min(pcs.store,ncol(pc.data))
            pcs.print=min(pcs.print,ncol(pc.data))
            object@pca.x=data.frame(pca.obj$x[,1:pcs.store])
            object@pca.rot=data.frame(pca.obj$rotation[,1:pcs.store])

            
            if (do.print) {
              pc_scores=object@pca.x
              for(i in 1:pcs.print) {
                code=paste("PC",i,sep="")
                sx=pc_scores[order(pc_scores[,code]),]
                print(code)
                print(rownames(sx[1:genes.print,]))
                print ("")
                
                print(rev(rownames(sx[(nrow(sx)-genes.print):nrow(sx),])))
                print ("")
                print ("")      
              }
            }
            return(object)
          }
)

#' Probability of detection by identity class
#'
#' For each gene, calculates the probability of detection for each identity
#' class.
#'
#'
#' @param object Seurat object
#' @param thresh.min Minimum threshold to define 'detected' (log-scale)
#' @return Returns a matrix with genes as rows, identity classes as columns.
#' @export
setGeneric("cluster.alpha", function(object,thresh.min=0) standardGeneric("cluster.alpha"))
#' @export
setMethod("cluster.alpha", "seurat", 
          function(object,thresh.min=0) {
            ident.use=object@ident
            data.all=data.frame(row.names = rownames(object@data))
            for(i in sort(unique(ident.use))) {
              temp.cells=which.cells(object,i)
              data.temp=apply(object@data[,temp.cells],1,function(x)return(length(x[x>thresh.min])/length(x)))
              data.all=cbind(data.all,data.temp)
              colnames(data.all)[ncol(data.all)]=i
            }
            colnames(data.all)=sort(unique(ident.use))
            return(data.all)
          }
)

#' Average PCA scores by identity class
#' 
#' Returns the PCA scores for an 'average' single cell in each identity class
#' 
#' @param object Seurat object
#' @return Returns a matrix with genes as rows, identity classes as columns 
#' @export
setGeneric("average.pca", function(object) standardGeneric("average.pca"))
#' @export
setMethod("average.pca", "seurat", 
          function(object) {
            ident.use=object@ident
            data.all=data.frame(row.names = colnames(object@pca.rot))
            for(i in levels(ident.use)) {
              temp.cells=which.cells(object,i)
              if (length(temp.cells)==1) {
                data.temp=apply(data.frame((object@pca.rot[c(temp.cells,temp.cells),])),2,mean)
              }
              if (length(temp.cells)>1) data.temp=apply(object@pca.rot[temp.cells,],2,mean)
              data.all=cbind(data.all,data.temp)
              colnames(data.all)[ncol(data.all)]=i
            }
            #colnames(data.all)=levels(ident.use)
            return((data.all))
          }
)

#' Averaged gene expression by identity class
#'
#' Returns gene expression for an 'average' single cell in each identity class
#'
#' Output is in log-space, but averaging is done in non-log space.
#'
#' @param object Seurat object
#' @param genes.use Genes to analyze. Default is all genes.
#' @return Returns a matrix with genes as rows, identity classes as columns.
#' @export
setGeneric("average.expression", function(object,genes.use=NULL) standardGeneric("average.expression"))
#' @export
setMethod("average.expression", "seurat", 
          function(object,genes.use=NULL) {
            genes.use=set.ifnull(genes.use,rownames(object@data))
            genes.use=ainb(genes.use,rownames(object@data))
            ident.use=object@ident
            data.all=data.frame(row.names = genes.use)
            for(i in levels(ident.use)) {
              temp.cells=which.cells(object,i)
              if (length(temp.cells)==1) data.temp=(object@data[genes.use,temp.cells])
              if (length(temp.cells)>1) data.temp=apply(object@data[genes.use,temp.cells],1,expMean)
              data.all=cbind(data.all,data.temp)
              colnames(data.all)[ncol(data.all)]=i
            }
            colnames(data.all)=levels(ident.use)
            return(data.all)
          }
)

#Internal, not documented for now
topGenesForDim=function(i,dim_scores,do.balanced=FALSE,num.genes=30,reduction.use="pca") {
  code=paste(translate.dim.code(reduction.use),i,sep="")
  if (do.balanced) {
    num.genes=round(num.genes/2)
    sx=dim_scores[order(dim_scores[,code]),]
    genes.1=(rownames(sx[1:num.genes,]))
    genes.2=(rownames(sx[(nrow(sx)-num.genes):nrow(sx),]))
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
setGeneric("icTopGenes", function(object,ic.use=1,num.genes=30,do.balanced=FALSE) standardGeneric("icTopGenes"))
#' @export
setMethod("icTopGenes", "seurat", 
          function(object,ic.use=1,num.genes=30,do.balanced=FALSE) {
            ic_scores=object@ica.x
            i=ic.use
            ic.top.genes=unique(unlist(lapply(i,topGenesForDim,ic_scores,do.balanced,num.genes,"ica")))
            return(ic.top.genes)
          }
)

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
setGeneric("pcTopGenes", function(object,pc.use=1,num.genes=30,use.full=FALSE,do.balanced=FALSE) standardGeneric("pcTopGenes"))
#' @export
setMethod("pcTopGenes", "seurat", 
          function(object,pc.use=1,num.genes=30,use.full=FALSE,do.balanced=FALSE) {
            pc_scores=object@pca.x
            i=pc.use
            if (use.full==TRUE) pc_scores = object@pca.x.full
            pc.top.genes=unique(unlist(lapply(i,topGenesForDim,pc_scores,do.balanced,num.genes,"pca")))
            return(pc.top.genes)
          }
)

#' Print the results of a PCA analysis
#'
#' Prints a set of genes that most strongly define a set of principal components
#'
#' @inheritParams viz.pca
#' @param pcs.print Set of PCs to print genes for
#' @param genes.print Number of genes to print for each PC
#' @return Only text output
#' @export
setGeneric("print.pca", function(object,pcs.print=1:5,genes.print=30,use.full=FALSE) standardGeneric("print.pca"))
#' @export
setMethod("print.pca", "seurat", 
          function(object,pcs.print=1:5,genes.print=30,use.full=FALSE) {
            genes.print.use=round(genes.print/2)
            for(i in pcs.print) {
              code=paste("PC",i,sep="")
              sx=pcTopGenes(object,i,genes.print.use*2,use.full,TRUE)
              print(code)
              print((sx[1:genes.print.use]))
              print ("")
              
              print(rev((sx[(length(sx)-genes.print.use):length(sx)])))
              print ("")
              print ("")      
            }
          }
)

#' Access cellular data
#'
#' Retreives data (gene expression, PCA scores, etc, metrics, etc.) for a set
#' of cells in a Seurat object
#'
#' 
#' @param object Seurat object
#' @param vars.all List of all variables to fetch
#' @param cells.use Cells to collect data for (default is all cells)
#' @param use.imputed For gene expression, use imputed values
#' @param use.scaled For gene expression, use scaled values
#' @return A data frame with cells as rows and cellular data as columns
#' @export
setGeneric("fetch.data",  function(object, vars.all=NULL,cells.use=NULL,use.imputed=FALSE, use.scaled=FALSE) standardGeneric("fetch.data"))
#' @export
setMethod("fetch.data","seurat",
          function(object, vars.all=NULL,cells.use=NULL,use.imputed=FALSE, use.scaled=FALSE) {
            cells.use=set.ifnull(cells.use,object@cell.names)
            data.return=data.frame(row.names = cells.use)
            # if any vars passed are genes, subset expression data
            gene_check <- vars.all %in% rownames(object@data)
            data.expression <- matrix()
            if (any(gene_check)){
              if (all(gene_check)){
                if(use.imputed) data.expression = object@imputed[vars.all, ]
                if(use.scaled) data.expression = object@scale.data[vars.all, ]
                else data.expression = object@data[vars.all, ]
                return(t(data.expression))
              }
              else{
                if(use.imputed) data.expression = object@imputed[vars.all[gene_check], ]
                if(use.scaled) data.expression = object@scale.data[vars.all[gene_check], ]
                else data.expression = object@data[vars.all[gene_check], ]
                data.expression = t(data.expression)
              }
            }
            var.options=c("data.info","pca.rot","ica.rot","tsne.rot","mix.probs","gene.scores")
            object@data.info[,"ident"]=object@ident[rownames(object@data.info)]
            for (my.var in vars.all) {
              data.use=data.frame()
              if (my.var %in% colnames(data.expression)) {
                data.use=data.expression
              } else {
                for(i in var.options) {
                  eval(parse(text=paste("data.use = object@",i,sep="")))
                  if (my.var %in% colnames(data.use)) {
                    break;
                  }
                }
              }              
              if (ncol(data.use)==0) {
                print(paste("Error : ", my.var, " not found", sep=""))
                return(0);
              }
              cells.use=ainb(cells.use,rownames(data.use))
              data.add=data.use[cells.use,my.var]
              if (is.null(data.add)) {
                print(paste("Error : ", my.var, " not found", sep=""))
                return(0);
              }
              data.return=cbind(data.return,data.add)  
            }
            colnames(data.return)=vars.all
            rownames(data.return)=cells.use
            return(data.return)
          }
)

#' Visualize ICA genes
#' 
#' Visualize top genes associated with principal components
#' 
#' 
#' @param object Seurat object
#' @param ics.use Number of ICs to display
#' @param num.genes Number of genes to display
#' @param font.size Font size
#' @param nCol Number of columns to display
#' @param do.balanced Return an equal number of genes with both + and - IC scores. 
#' If FALSE (by default), returns the top genes ranked by the score's absolute values 
#' @return Graphical, no return value
#' @export
setGeneric("viz.ica", function(object,ics.use=1:5,num.genes=30,font.size=0.5,nCol=NULL,do.balanced=FALSE) standardGeneric("viz.ica"))
#' @export
setMethod("viz.ica", "seurat", 
          function(object,ics.use=1:5,num.genes=30,font.size=0.5,nCol=NULL,do.balanced=FALSE) {
            ic_scores=object@ica.x            
            if (is.null(nCol)) {
              nCol=2
              if (length(ics.use)>6) nCol=3
              if (length(ics.use)>9) nCol=4
            }         
            num.row=floor(length(ics.use)/nCol-1e-5)+1
            par(mfrow=c(num.row,nCol))
            
            for(i in ics.use) {
              subset.use=ic_scores[icTopGenes(object,i,num.genes,do.balanced),]
              print(head(subset.use))
              plot(subset.use[,i],1:nrow(subset.use),pch=16,col="blue",xlab=paste("ic",i,sep=""),yaxt="n",ylab="")
              axis(2,at=1:nrow(subset.use),labels = rownames(subset.use),las=1,cex.axis=font.size)   
            }
            rp()
          }
)


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
setGeneric("viz.pca", function(object,pcs.use=1:5,num.genes=30,use.full=FALSE,font.size=0.5,nCol=NULL,do.balanced=FALSE) standardGeneric("viz.pca"))
#' @export
setMethod("viz.pca", "seurat", 
          function(object,pcs.use=1:5,num.genes=30,use.full=FALSE,font.size=0.5,nCol=NULL,do.balanced=FALSE) {
            pc_scores=object@pca.x
            if (use.full==TRUE) pc_scores = object@pca.x.full
            
            if (is.null(nCol)) {
              nCol=2
              if (length(pcs.use)>6) nCol=3
              if (length(pcs.use)>9) nCol=4
            }         
            num.row=floor(length(pcs.use)/nCol-1e-5)+1
            par(mfrow=c(num.row,nCol))
            
            for(i in pcs.use) {
              subset.use=pc_scores[pcTopGenes(object,i,num.genes,use.full,do.balanced),]
              plot(subset.use[,i],1:nrow(subset.use),pch=16,col="blue",xlab=paste("PC",i,sep=""),yaxt="n",ylab="")
              axis(2,at=1:nrow(subset.use),labels = rownames(subset.use),las=1,cex.axis=font.size)   
            }
            rp()
          }
)


set.ifnull=function(x,y) {
  if(is.null(x)) x=y
  return(x)
}

kill.ifnull=function(x,message="Error:Execution Halted") {
  if(is.null(x)) {
    stop(message)
  }
}

expAlpha=function(mu,coefs) {
  logA=coefs$a
  logB=coefs$b
  return(exp(logA+logB*mu)/(1+(exp(logA+logB*mu))))
}

setGeneric("getWeightMatrix", function(object) standardGeneric("getWeightMatrix"))
setMethod("getWeightMatrix", "seurat",
          function(object) {
            data=object@data
            data.humpAvg=apply(data,1,humpMean,min=object@drop.expr)
            wt.matrix=data.frame(t(sapply(data.humpAvg,expAlpha,object@drop.coefs)))
            colnames(wt.matrix)=colnames(data); rownames(wt.matrix)=rownames(data)
            wt.matrix[is.na(wt.matrix)]=0
            object@wt.matrix=wt.matrix
            wt1.matrix=data.frame(sapply(1:ncol(data),function(x)setWt1(data[,x],wt.matrix[,x],min=object@drop.expr)))
            colnames(wt1.matrix)=colnames(data); rownames(wt1.matrix)=rownames(data)
            wt1.matrix[is.na(wt1.matrix)]=0
            object@drop.wt.matrix=wt1.matrix
            return(object)
          }      
)

regression.sig=function(x,score,data,latent,code="rsem") {
  if(var(as.numeric(subc(data,code)[x,]))==0) {
    return(0)
  }
  latent=latent[grep(code,names(data))]
  data=rbind(subc(data,code),vsubc(score,code))
  rownames(data)[nrow(data)]="score"
  data2=data[c(x,"score"),]
  rownames(data2)[1]="fac"
  if (length(unique(latent))>1) {
    mylm=lm(score ~ fac+latent, data = data.frame(t(data2)))
  }
  else {
    mylm=lm(score ~ fac, data = data.frame(t(data2)))
  }
  return(coef(summary(mylm))["fac",3])
}

setGeneric("regulatorScore", function(object, candidate.reg, score.name, cells.use=NULL) standardGeneric("regulatorScore"))
setMethod("regulatorScore", "seurat",
          function(object, candidate.reg, score.name, cells.use=NULL) {
            cells.use=set.ifnull(cells.use, colnames(object@data))
            candidate.reg=candidate.reg[candidate.reg%in%rownames(object@data)]
            my.score=retreiveScore(object,score.name)[cells.use]
            my.data=object@data[,cells.use]
            my.ident=object@ident[cells.use]
            reg.score=unlist(lapply(candidate.reg,regressionSig,score = my.score,data = my.data,latent = my.ident,code = "rsem"))
            names(reg.score)=candidate.reg
            return(reg.score)
          }
)

#' Gene expression markers of identity classes defined by a phylogenetic clade
#' 
#' Finds markers (differentially expressed genes) based on a branching point (node) in
#' the phylogenetic tree. Markers that define clusters in the left branch are positive markers.
#' Markers that define the right branch are negative markers.
#' 
#' @inheritParams find.markers
#' @param node The node in the phylogenetic tree to use as a branch point 
#' @return Matrix containing a ranked list of putative markers, and associated
#' identistics (p-values, ROC score, etc.)
#' @export
setGeneric("find.markers.node", function(object,node,genes.use=NULL,thresh.use=log(2),test.use="bimod") standardGeneric("find.markers.node"))
#' @export
setMethod("find.markers.node", "seurat",
          function(object,node,genes.use=NULL,thresh.use=log(2),test.use="bimod") {
            genes.use=set.ifnull(genes.use,rownames(object@data))
            tree=object@cluster.tree[[1]]
            ident.order=tree$tip.label
            nodes.1=ident.order[getLeftDecendants(tree,node)]
            nodes.2=ident.order[getRightDecendants(tree,node)]
            #print(nodes.1)
            #print(nodes.2)
            to.return=find.markers(object,nodes.1,nodes.2,genes.use,thresh.use,test.use)
            return(to.return)
          } 
)

#' Gene expression markers of identity classes
#' 
#' Finds markers (differentially expressed genes) for identity classes
#' 
#' 
#' @param object Seurat object
#' @param ident.1 Identity class to define markers for
#' @param ident.2 A second identity class for comparison. If NULL (default) -
#' use all other cells for comparison.
#' @param genes.use Genes to test. Default is to use all genes.
#' @param thresh.use Limit testing to genes which show, on average, at least
#' X-fold difference (log-scale) between the two groups of cells.
#' 
#' Increasing thresh.use speeds up the function, but can miss weaker signals.
#' @param test.use Denotes which test to use. Seurat currently implements
#' "bimod" (likelihood-ratio test for single cell gene expression, McDavid et
#' al., Bioinformatics, 2011, default), "roc" (standard AUC classifier), "t"
#' (Students t-test), and "tobit" (Tobit-test for differential gene expression,
#' as in Trapnell et al., Nature Biotech, 2014)
#' @return Matrix containing a ranked list of putative markers, and associated
#' identistics (p-values, ROC score, etc.)
#' @import VGAM
#' @export
setGeneric("find.markers", function(object, ident.1,ident.2=NULL,genes.use=NULL,thresh.use=log(2),test.use="bimod") standardGeneric("find.markers"))
#' @export
setMethod("find.markers", "seurat",
          function(object, ident.1,ident.2=NULL,genes.use=NULL,thresh.use=log(2), test.use="bimod") {
            genes.use=set.ifnull(genes.use,rownames(object@data))
            
            cells.1=which.cells(object,ident.1)
            # in case the user passed in cells instead of identity classes
            if (length(ident.1>1)&&any(ident.1%in%object@cell.names)) {
              cells.1=ainb(ident.1,object@cell.names)
            }
            
            # if NULL for ident.2, use all other cells
            if (is.null(ident.2)) {
              cells.2=object@cell.names
            }
            else {
              cells.2=which.cells(object,ident.2)
            }
            if (length(ident.2>1)&&any(ident.2%in%object@cell.names)) {
              cells.2=ainb(ident.2,object@cell.names)
            }
            cells.2=anotinb(cells.2,cells.1)
            
            #error checking
            if (length(cells.1)==0) {
              print(paste("Cell group 1 is empty - no cells with identity class", ident.1))
              return(NULL)
            }
            if (length(cells.2)==0) {
              print(paste("Cell group 2 is empty - no cells with identity class", ident.2))
              return(NULL)
            }
            
            if (test.use=="bimod") to.return=diffExp.test(object,cells.1,cells.2,genes.use,thresh.use) 
            if (test.use=="roc") to.return=marker.test(object,cells.1,cells.2,genes.use,thresh.use) 
            if (test.use=="t") to.return=diff.t.test(object,cells.1,cells.2,genes.use,thresh.use) 
            if (test.use=="tobit") to.return=tobit.test(object,cells.1,cells.2,genes.use,thresh.use) 
            
            return(to.return)
          } 
)


#' Gene expression markers for all identity classes
#' 
#' Finds markers (differentially expressed genes) for each of the identity classes in a dataset
#' 
#' 
#' @inheritParams find.markers
#' @param thresh.test Limit testing to genes which show, on average, at least X-fold difference (log-scale) between cells in an identity class, and all other cells.
#' @param return.thresh Only return markers that have a p-value < return.thresh, or a power > return.thresh (if the test is ROC)
#' @param do.print FALSE by default. If TRUE, outputs updates on progress.
#' @return Matrix containing a ranked list of putative markers, and associated
#' identistics (p-values, ROC score, etc.)
#' @export
setGeneric("find_all_markers", function(object, thresh.test=1,test.use="bimod",return.thresh=1e-2,do.print=FALSE) standardGeneric("find_all_markers"))
#' @export
setMethod("find_all_markers","seurat",
          function(object, thresh.test=1,test.use="bimod",return.thresh=1e-2,do.print=FALSE) {
            ident.use=object@ident
            if ((test.use=="roc") && (return.thresh==1e-2)) return.thresh=0.7
            idents.all=sort(unique(object@ident))
            genes.de=list()
            for(i in 1:length(idents.all)) {
              genes.de[[i]]=find.markers(object,idents.all[i],genes.use=rownames(object@data),thresh.use = thresh.test,test.use = test.use)
              if (do.print) print(paste("Calculating cluster", idents.all[i]))
            }
            gde.all=data.frame()
            for(i in 1:length(idents.all)) {
              gde=genes.de[[i]]
              if (nrow(gde)>0) {
                if (test.use=="roc") gde=subset(gde,(myAUC>return.thresh|myAUC<(1-return.thresh)))
                if (test.use=="bimod") {
                  gde=gde[order(gde$p_val,-gde$avg_diff),]
                  gde=subset(gde,p_val<return.thresh)
                }
                if (nrow(gde)>0) gde$cluster=idents.all[i]; gde$gene=rownames(gde)
                if (nrow(gde)>0) gde.all=rbind(gde.all,gde)
              }
            }
            return(gde.all)
          }
)


#' Likelihood ratio test for zero-inflated data
#'
#' Identifies differentially expressed genes between two groups of cells using
#' the LRT model proposed in Mcdavid et al, Bioinformatics, 2011
#' 
#' @param object Seurat object
#' @param cells.1 Group 1 cells
#' @param cells.2 Group 2 cells
#' @param genes.use Genes to test. Default is to use all genes.
#' @param thresh.use Limit testing to genes which show, on average, at least
#' X-fold difference (log-scale) between the two groups of cells.
#' 
#' Increasing thresh.use speeds up the function, but can miss weaker signals.
#' @return Returns a p-value ranked matrix of putative differentially expressed
#' genes.
#' @export
setGeneric("diffExp.test", function(object, cells.1,cells.2,genes.use=NULL,thresh.use=log(2)) standardGeneric("diffExp.test"))
#' @export
setMethod("diffExp.test", "seurat",
          function(object, cells.1,cells.2,genes.use=NULL,thresh.use=log(2)) {
            genes.use=set.ifnull(genes.use,object@var.genes)
            data.1=apply(object@data[genes.use,cells.1],1,expMean)
            data.2=apply(object@data[genes.use,cells.2],1,expMean)
            total.diff=abs(data.1-data.2)
            genes.diff = names(which(total.diff>thresh.use))
            #print(genes.diff)
            to.return=bimod.diffExp.test(object@data[,cells.1],object@data[,cells.2],genes.diff)
            to.return=to.return[order(to.return$p_val,-abs(to.return$avg_diff)),]
            return(to.return)
          } 
)

#' Differential expression testing using Tobit models
#'
#' Identifies differentially expressed genes between two groups of cells using
#' Tobit models, as proposed in Trapnell et al., Nature Biotechnology, 2014
#'
#' @param object Seurat object
#' @param cells.1 Group 1 cells
#' @param cells.2 Group 2 cells
#' @param genes.use Genes to test. Default is to use all genes.
#' @param thresh.use Limit testing to genes which show, on average, at least
#' X-fold difference (log-scale) between the two groups of cells.
#'
#' Increasing thresh.use speeds up the function, but can miss weaker signals.
#' @return Returns a p-value ranked matrix of putative differentially expressed
#' genes.
#' @export
setGeneric("tobit.test", function(object, cells.1,cells.2,genes.use=NULL,thresh.use=log(2)) standardGeneric("tobit.test"))
#' @export
setMethod("tobit.test", "seurat",
          function(object, cells.1,cells.2,genes.use=NULL,thresh.use=log(2)) {
            genes.use=set.ifnull(genes.use,object@var.genes)
            data.1=apply(object@data[genes.use,cells.1],1,expMean)
            data.2=apply(object@data[genes.use,cells.2],1,expMean)
            total.diff=abs(data.1-data.2)
            genes.diff = names(which(total.diff>thresh.use))
            #print(genes.diff)
            to.return=tobit.diffExp.test(object@data[,cells.1],object@data[,cells.2],genes.diff)
            to.return=to.return[order(to.return$p_val,-abs(to.return$avg_diff)),]
            return(to.return)
          } 
)

#' Identify potential genes associated with batch effects
#' 
#' Test for genes whose expression value is strongly predictive of batch (based
#' on ROC classification). Important note: Assumes that the 'batch' of each
#' cell is assigned to the cell's identity class (will be improved in a future
#' release)
#' 
#' @param object Seurat object
#' @param idents.use Batch names to test
#' @param genes.use Gene list to test
#' @param auc.cutoff Minimum AUC needed to qualify as a 'batch gene'
#' @param thresh.use Limit testing to genes which show, on average, at least
#' X-fold difference (log-scale) in any one batch
#' @return Returns a list of genes that are strongly correlated with batch.
#' @export
setGeneric("batch.gene", function(object, idents.use,genes.use=NULL,auc.cutoff=0.6,thresh.use=0) standardGeneric("batch.gene"))
#' @export
setMethod("batch.gene", "seurat",
          function(object, idents.use,genes.use=NULL,auc.cutoff=0.6,thresh.use=0) {
            batch.genes=c()
            genes.use=set.ifnull(genes.use,rownames(object@data))
            for(ident in idents.use ) {
              cells.1=names(object@ident)[object@ident==ident]
              cells.2=names(object@ident)[object@ident!=ident]
              if ((length(cells.1)<5)|(length(cells.2)<5)) {
                break;
              }
              markers.ident=marker.test(object,cells.1,cells.2,genes.use,thresh.use)
              batch.genes=unique(c(batch.genes,rownames(subset(markers.ident,myAUC>auc.cutoff))))
            }
            return(batch.genes)
          } 
)

#' ROC-based marker discovery
#'
#' Identifies 'markers' of gene expression using ROC analysis. For each gene,
#' evaluates (using AUC) a classifier built on that gene alone, to classify
#' between two groups of cells.
#'
#' An AUC value of 1 means that expression values for this gene alone can
#' perfectly classify the two groupings (i.e. Each of the cells in cells.1
#' exhibit a higher level than each of the cells in cells.2). An AUC value of 0
#' also means there is perfect classification, but in the other direction. A
#' value of 0.5 implies that the gene has no predictive power to classify the
#' two groups.
#'
#' @param object Seurat object
#' @param cells.1 Group 1 cells
#' @param cells.2 Group 2 cells
#' @param genes.use Genes to test. Default is to use all genes.
#' @param thresh.use Limit testing to genes which show, on average, at least
#' X-fold difference (log-scale) between the two groups of cells.
#'
#' Increasing thresh.use speeds up the function, but can miss weaker signals.
#' @return Returns a 'predictive power' (abs(AUC-0.5)) ranked matrix of
#' putative differentially expressed genes.
#' @import ROCR
#' @export
setGeneric("marker.test", function(object, cells.1,cells.2,genes.use=NULL,thresh.use=log(2)) standardGeneric("marker.test"))
#' @export
setMethod("marker.test", "seurat",
          function(object, cells.1,cells.2,genes.use=NULL,thresh.use=log(2)) {
            genes.use=set.ifnull(genes.use,object@var.genes)
            data.use=object@data
            data.1=apply(object@data[genes.use,cells.1],1,expMean)
            data.2=apply(object@data[genes.use,cells.2],1,expMean)
            total.diff=abs(data.1-data.2)
            genes.diff = names(which(total.diff>thresh.use))
            genes.use=ainb(genes.diff,rownames(data.use))
            to.return=marker.auc.test(object@data[,cells.1],object@data[,cells.2],genes.use)
            to.return=to.return[rev(order(abs(to.return$myAUC-0.5))),]
            to.return$power=abs(to.return$myAUC-0.5)*2
            return(to.return)
          } 
)

#' Differential expression testing using Student's t-test
#'
#' Identify differentially expressed genes between two groups of cells using
#' the Student's t-test
#' 
#' @param object Seurat object
#' @param cells.1 Group 1 cells
#' @param cells.2 Group 2 cells
#' @param genes.use Genes to test. Default is to use all genes.
#' @param thresh.use Limit testing to genes which show, on average, at least
#' X-fold difference (log-scale) between the two groups of cells.
#' 
#' Increasing thresh.use speeds up the function, but can miss weaker signals.
#' @return Returns a p-value ranked matrix of putative differentially expressed
#' genes.
#' @export
setGeneric("diff.t.test", function(object, cells.1,cells.2,genes.use=NULL,thresh.use=log(2)) standardGeneric("diff.t.test"))
#' @export
setMethod("diff.t.test", "seurat",
          function(object, cells.1,cells.2,genes.use=NULL,thresh.use=log(2)) {
            genes.use=set.ifnull(genes.use,object@var.genes)
            data.use=object@data
            data.1=apply(object@data[genes.use,cells.1],1,expMean)
            data.2=apply(object@data[genes.use,cells.2],1,expMean)
            total.diff=abs(data.1-data.2)
            genes.diff = names(which(total.diff>thresh.use))
            genes.use=ainb(genes.diff,rownames(data.use))
            p_val=unlist(lapply(genes.use,function(x)t.test(object@data[x,cells.1],object@data[x,cells.2])$p.value))
            avg_diff=(data.1-data.2)[genes.use]
            to.return=data.frame(p_val,avg_diff,row.names = genes.use)
            to.return=to.return[with(to.return, order(p_val, -abs(avg_diff))), ]
            return(to.return)
          } 
)

#' Identify matching cells
#'
#' Returns a list of cells that match a particular query (usually, query is
#' based on identity class). For example, to find the names of all cells in cluster 1.
#'
#' @param object Seurat object
#' @param value Query value to match
#' @param id Variable to query (by default, identity class)
#' @return A vector of cell names
#' @export
setGeneric("which.cells", function(object,value=1, id=NULL) standardGeneric("which.cells"))
#' @export
setMethod("which.cells", "seurat",
          function(object, value=1,id=NULL) {
            id=set.ifnull(id,"ident")
            data.use=NULL;
            if (id=="ident") {
              data.use=object@ident
            } else {
              if (id %in% colnames(object@data.info)) {
                data.use=object@data.info[,id]; names(data.use)=rownames(object@data.info)
              }
            }
            if (is.null(data.use)) {
              stop(paste("Error : ", id, " not found"))
            }
            return(names(data.use[which(data.use%in%value)]))
          } 
)


#' Switch identity class definition to another variable
#'
#'
#' @param object Seurat object
#' @param id Variable to switch identity class to (for example, 'DBclust.ident', the output
#' of density clustering) Default is orig.ident - the original annotation pulled from the cell name.
#' @return A Seurat object where object@@ident has been appropriately modified
#' @export
setGeneric("set.all.ident", function(object,id=NULL) standardGeneric("set.all.ident"))
#' @export
setMethod("set.all.ident", "seurat",
          function(object, id=NULL) {
            id=set.ifnull(id,"orig.ident")
            if (id %in% colnames(object@data.info)) {
              cells.use=rownames(object@data.info)
              ident.use=object@data.info[,id]
              object=set.ident(object,cells.use,ident.use)
            }
            return(object)
          } 
)

#' Rename one identity class to another
#'
#' Can also be used to join identity classes together (for example, to merge clusters).
#'
#' @param object Seurat object
#' @param old.ident.name The old identity class (to be renamed)
#' @param new.ident.name The new name to apply
#' @return A Seurat object where object@@ident has been appropriately modified
#' @export
setGeneric("rename.ident", function(object,old.ident.name=NULL,new.ident.name=NULL) standardGeneric("rename.ident"))
#' @export
setMethod("rename.ident", "seurat",
          function(object,old.ident.name=NULL,new.ident.name=NULL) {
            if (!old.ident.name%in%object@ident) {
              stop(paste("Error : ", old.ident.name, " is not a current identity class"))
            }
            old.levels=levels(object@ident); new.levels=old.levels
            if(new.ident.name%in%old.levels) {
              new.levels=new.levels[new.levels!=old.ident.name]
            }
            if(!(new.ident.name%in%old.levels)) {
              new.levels[new.levels==old.ident.name]=new.ident.name
            }
            ident.vector=as.character(object@ident); names(ident.vector)=names(object@ident)
            ident.vector[which.cells(object,old.ident.name)]=new.ident.name            
            object@ident=factor(ident.vector,levels = new.levels)
            return(object)
          } 
)

#' Set identity class information
#'
#' Sets the identity class value for a subset (or all) cells
#'
#' @param object Seurat object
#' @param cells.use Vector of cells to set identity class info for (default is
#' all cells)
#' @param ident.use Vector of identity class values to assign (character
#' vector)
#' @return A Seurat object where object@@ident has been appropriately modified
#' @importFrom gdata drop.levels
#' @export
setGeneric("set.ident", function(object,cells.use=NULL,ident.use=NULL) standardGeneric("set.ident"))
#' @export
setMethod("set.ident", "seurat",
          function(object, cells.use=NULL,ident.use=NULL) {
            cells.use=set.ifnull(cells.use,object@cell.names)
            if (length(anotinb(cells.use,object@cell.names)>0)) {
              stop(paste("ERROR : Cannot find cells ",anotinb(cells.use,object@cell.names)))
            }
            ident.new=anotinb(ident.use,levels(object@ident))
            object@ident=factor(object@ident,levels = unique(c(as.character(object@ident),as.character(ident.new))))
            object@ident[cells.use]=ident.use
            object@ident=drop.levels(object@ident)
            return(object)
          } 
)

#Not documented for now
#' @export
setGeneric("posterior.plot", function(object, name) standardGeneric("posterior.plot"))
#' @export
setMethod("posterior.plot", "seurat",
          function(object, name) {
            post.names=colnames(subc(object@mix.probs,name))
            vlnPlot(object,post.names,inc.first=TRUE,inc.final=TRUE,by.k=TRUE)
            
            
          } 
)

#Internal, not documented for now
map.cell.score=function(gene,gene.value,insitu.bin,mu,sigma,alpha) {
  code.1=paste(gene,insitu.bin,sep=".")
  mu.use=mu[paste(code.1,"mu",sep="."),1]
  sigma.use=sigma[paste(code.1,"sigma",sep="."),1]
  alpha.use=alpha[paste(code.1,"alpha",sep="."),1]
  bin.prob=unlist(lapply(1:length(insitu.bin),function(x) dnorm(gene.value,mean = mu.use[x],sd = sigma.use[x],log = TRUE) + log(alpha.use[x])))
  return(bin.prob)
}

#Internal, not documented for now
#' @export
setGeneric("map.cell",  function(object,cell.name,do.plot=FALSE,safe.use=TRUE,text.val=NULL,do.rev=FALSE) standardGeneric("map.cell"))
#' @export
setMethod("map.cell", "seurat",
          function(object,cell.name,do.plot=FALSE,safe.use=TRUE,text.val=NULL,do.rev=FALSE) {
            insitu.matrix=object@insitu.matrix
            insitu.genes=colnames(insitu.matrix)
            insitu.genes=insitu.genes[insitu.genes%in%rownames(object@imputed)]
            insitu.use=insitu.matrix[,insitu.genes]
            imputed.use=object@imputed[insitu.genes,]
            safe_fxn=sum
            if (safe.use) safe_fxn=log_add
            
            all.needed.cols=unique(unlist(lapply(insitu.genes,function(x) paste(x,insitu.use[,x],"post",sep="."))))
            missing.cols=which(!(all.needed.cols%in%colnames(object@mix.probs)))
            if (length(missing.cols)>0) stop(paste("Error : ", all.needed.cols[missing.cols], " is missing from the mixture fits",sep=""))
            all.probs=data.frame(sapply(insitu.genes,function(x) log(as.numeric(object@mix.probs[cell.name,paste(x,insitu.use[,x],"post",sep=".")]))))
            scale.probs=t(t(all.probs)-apply(t(all.probs),1,log_add))
            scale.probs[scale.probs<(-9.2)]=(-9.2)
            #head(scale.probs)
            total.prob=exp(apply(scale.probs,1,safe_fxn))
            total.prob=total.prob/sum(total.prob)
            if (do.plot) {
              #plot(total.prob,main=cell.name)
              par(mfrow=c(1,2))
              txt.matrix=matrix(rep("",64),nrow=8,ncol=8)
              if (!is.null(text.val)) txt.matrix[text.val]="X"
              if (do.rev) scale.probs=scale.probs[unlist(lapply(0:7,function(x)seq(1,57,8)+x)),] 
              aheatmap(matrix(total.prob,nrow=8,ncol=8),Rowv=NA,Colv=NA,txt=txt.matrix,col=bwCols)   
              aheatmap(scale.probs,Rowv=NA,Colv=NA)
              rp()
            }
            return(total.prob)
          } 
)

#' Get cell centroids
#'
#' Calculate the spatial mapping centroids for each cell, based on previously
#' calculated mapping probabilities for each bin.
#'
#' Currently, Seurat assumes that the tissue of interest has an 8x8 bin
#' structure. This will be broadened in a future release.
#' 
#' @param object Seurat object
#' @param cells.use Cells to calculate centroids for (default is all cells)
#' @param get.exact Get exact centroid (Default is TRUE). If FALSE, identify
#' the single closest bin.
#' @return Data frame containing the x and y coordinates for each cell
#' centroid.
#' @export
setGeneric("get.centroids", function(object, cells.use=NULL,get.exact=TRUE) standardGeneric("get.centroids")) 
#' @export
setMethod("get.centroids", "seurat",
          function(object, cells.use=NULL,get.exact=TRUE) {            
            cells.use=set.ifnull(cells.use,colnames(object@final.prob))
            
            #Error checking
            cell.names=ainb(cells.use, colnames(object@final.prob))
            if (length(cell.names)!=length(cells.use)) {
              print(paste("Error", anotinb(cells.use,colnames(object@final.prob)), " have not been mapped"))
              return(0);
            }
            
            if (get.exact) my.centroids=data.frame(t(sapply(colnames(object@data),function(x) exact.cell.centroid(object@final.prob[,x])))); colnames(my.centroids)=c("bin.x","bin.y")
            if (!(get.exact)) my.centroids=data.frame(t(sapply(colnames(object@data),function(x) cell.centroid(object@final.prob[,x])))); colnames(my.centroids)=c("bin.x","bin.y")
            
            return(my.centroids)
          }
)

#' Quantitative refinement of spatial inferences
#'
#' Refines the initial mapping with more complex models that allow gene
#' expression to vary quantitatively across bins (instead of 'on' or 'off'),
#' and that also considers the covariance structure between genes.
#'
#' Full details given in spatial mapping manuscript.
#'
#' @param object Seurat object
#' @param genes.use Genes to use to drive the refinement procedure.
#' @return Seurat object, where mapping probabilities for each bin are stored
#' in object@@final.prob
#' @import fpc
#' @export
setGeneric("refined.mapping",  function(object,genes.use) standardGeneric("refined.mapping"))
#' @export
setMethod("refined.mapping", "seurat",
          function(object,genes.use) {
            
            genes.use=ainb(genes.use, rownames(object@imputed))            
            cells.max=t(sapply(colnames(object@data),function(x) exact.cell.centroid(object@final.prob[,x])))
            all.mu=sapply(genes.use,function(gene) sapply(1:64, function(bin) mean(as.numeric(object@imputed[gene,fetch.closest(bin,cells.max,2*length(genes.use))]))))
            all.cov=list(); for(x in 1:64) all.cov[[x]]=cov(t(object@imputed[genes.use,fetch.closest(x,cells.max,2*length(genes.use))]))
            
            mv.probs=sapply(colnames(object@data),function(my.cell) sapply(1:64,function(bin) slimdmvnorm(as.numeric(object@imputed[genes.use,my.cell]),as.numeric(all.mu[bin,genes.use]),all.cov[[bin]])))
            mv.final=exp(sweep(mv.probs,2,apply(mv.probs,2,log_add)))
            object@final.prob=data.frame(mv.final)            
            return(object)
          } 
)

#' Infer spatial origins for single cells
#'
#' Probabilistically maps single cells based on (imputed) gene expression
#' estimates, a set of mixture models, and an in situ spatial reference map.
#'
#'
#' @param object Seurat object
#' @param cells.use Which cells to map
#' @return Seurat object, where mapping probabilities for each bin are stored
#' in object@@final.prob
#' @export
setGeneric("initial.mapping", function(object,cells.use=NULL) standardGeneric("initial.mapping"))
#' @export
setMethod("initial.mapping", "seurat",
          function(object,cells.use=NULL) {
            cells.use=set.ifnull(cells.use,colnames(object@data))
            every.prob=sapply(cells.use,function(x)map.cell(object,x,FALSE,FALSE))
            object@final.prob=data.frame(every.prob)
            rownames(object@final.prob)=paste("bin.",rownames(object@final.prob),sep="")
            return(object)
          } 
)

#Internal, not documented for now
setGeneric("calc.insitu", function(object,gene,do.plot=TRUE,do.write=FALSE,write.dir="~/window/insitu/",col.use=bwCols,do.norm=FALSE,cells.use=NULL,do.return=FALSE,probs.min=0,do.log=FALSE, use.imputed=FALSE, bleach.use=0) standardGeneric("calc.insitu"))
setMethod("calc.insitu", "seurat",
          function(object,gene,do.plot=TRUE,do.write=FALSE,write.dir="~/window/insitu/",col.use=bwCols,do.norm=FALSE,cells.use=NULL,do.return=FALSE,probs.min=0,do.log=FALSE,use.imputed=FALSE,bleach.use=0) {
            cells.use=set.ifnull(cells.use,colnames(object@final.prob))
            probs.use=object@final.prob
            data.use=exp(object@data)-1  
            if (use.imputed) data.use=exp(object@imputed)-1
            cells.use=cells.use[cells.use%in%colnames(probs.use)]; cells.use=cells.use[cells.use%in%colnames(data.use)]
            #insilico.stain=matrix(unlist(lapply(1:64,function(x) sum(probs.use[x,]*data.use[gene,]))),nrow=8,ncol=8)
            insilico.vector=unlist(lapply(1:64,function(x) sum(as.numeric(probs.use[x,cells.use])*as.numeric(data.use[gene,cells.use]))))
            probs.total=apply(probs.use,1,sum)
            probs.total[probs.total<probs.min]=probs.min
            insilico.stain=(matrix(insilico.vector/probs.total,nrow=8,ncol=8))
            if (do.log) insilico.stain=log(insilico.stain+1)
            if (bleach.use > 0) {
              insilico.stain=insilico.stain-bleach.use
              insilico.stain=minmax(insilico.stain,min=0,max=1e6)
            }
            if (do.norm) insilico.stain=(insilico.stain-min(insilico.stain))/(max(insilico.stain)-min(insilico.stain))
            title.use=gene
            if (gene %in% colnames(object@insitu.matrix)) {
              pred.use=prediction(insilico.vector/probs.total,object@insitu.matrix[,gene],0:1)
              perf.use=performance(pred.use,"auc")
              auc.use=round(perf.use@y.values[[1]],3)
              title.use=paste(gene,sep=" ")
            }
            if (do.write) {
              write.table(insilico.stain,paste(write.dir,gene,".txt",sep=""),quote=FALSE,row.names=FALSE,col.names=FALSE)
            }
            if (do.plot) {
              aheatmap(insilico.stain,Rowv=NA,Colv=NA,col=col.use, main=title.use)
            }
            if (do.return) {
              return(as.vector(insilico.stain))
            }
            return(object)
          } 
)

#' Build mixture models of gene expression
#'
#' Models the imputed gene expression values as a mixture of gaussian
#' distributions. For a two-state model, estimates the probability that a given
#' cell is in the 'on' or 'off' state for any gene. Followed by a greedy
#' k-means step where cells are allowed to flip states based on the overall
#' structure of the data (see Manuscript for details)
#' 
#' 
#' @param object Seurat object
#' @param gene Gene to fit
#' @param do.k Number of modes for the mixture model (default is 2)
#' @param num.iter Number of 'greedy k-means' iterations (default is 1)
#' @param do.plot Plot mixture model results
#' @param genes.use Genes to use in the greedy k-means step (See manuscript for details)
#' @param start.pct Initial estimates of the percentage of cells in the 'on'
#' state (usually estimated from the in situ map)
#' @return A Seurat object, where the posterior of each cell being in the 'on'
#' or 'off' state for each gene is stored in object@@mix.probs
#' @importFrom mixtools normalmixEM
#' @export
setGeneric("fit.gene.k", function(object, gene, do.k=2,num.iter=1,do.plot=FALSE,genes.use=NULL,start.pct=NULL) standardGeneric("fit.gene.k"))
#' @export
setMethod("fit.gene.k", "seurat",
          function(object, gene, do.k=2,num.iter=1,do.plot=FALSE,genes.use=NULL,start.pct=NULL) {
            data=object@imputed            
            data.use=data[gene,]
            names(data.use)=colnames(data.use)
            scale.data=t(scale(t(object@imputed)))
            genes.use=set.ifnull(genes.use,rownames(scale.data))
            genes.use=genes.use[genes.use%in%rownames(scale.data)]
            scale.data=scale.data[genes.use,]

            data.cut=as.numeric(data.use[gene,])
            cell.ident=as.numeric(cut(data.cut,do.k))
            if (!(is.null(start.pct))) {
              cell.ident=rep(1,length(data.cut))
              cell.ident[data.cut>quantile(data.cut,1-start.pct)]=2
            }
            cell.ident=order(tapply(as.numeric(data.use),cell.ident,mean))[cell.ident]
            ident.table=table(cell.ident)
            if (num.iter > 0) {
              for(i2 in 1:num.iter) {
                cell.ident=iter.k.fit(scale.data,cell.ident,data.use)
                ident.table=table(cell.ident)
              }
            }
            ident.table=table(cell.ident)
            raw.probs=t(sapply(data.use,function(y) unlist(lapply(1:do.k,function(x) ((ident.table[x]/sum(ident.table))*dnorm(y,mean(as.numeric(data.use[cell.ident==x])),sd(as.numeric(data.use[cell.ident==x]))))))))
            norm.probs=raw.probs/apply(raw.probs,1,sum)
            colnames(norm.probs)=unlist(lapply(1:do.k,function(x)paste(gene,x-1,"post",sep=".")))
            norm.probs=cbind(norm.probs,cell.ident); colnames(norm.probs)[ncol(norm.probs)]=paste(gene,".ident",sep="")
            new.mix.probs=data.frame(minusc(object@mix.probs,paste(gene,".",sep="")),row.names = rownames(object@mix.probs)); colnames(new.mix.probs)[1]="nGene"
            object@mix.probs=cbind(new.mix.probs,norm.probs)
            
            if (do.plot) {
              nCol=2
              num.row=floor((do.k+1)/nCol-1e-5)+1
              hist(as.numeric(data.use),probability = TRUE,ylim=c(0,1),xlab=gene,main=gene);
              for(i in 1:do.k) lines(seq(-10,10,0.01),(ident.table[i]/sum(ident.table)) * dnorm(seq(-10,10,0.01),mean(as.numeric(data.use[cell.ident==i])),sd(as.numeric(data.use[cell.ident==i]))),col=i,lwd=2); 
            }
            return(object)
          }
)

#Internal, not documented for now
iter.k.fit=function(scale.data,cell.ident,data.use) {
  means.all=sapply(sort(unique(cell.ident)),function(x)apply(scale.data[,cell.ident==x],1,mean))
  all.dist=data.frame(t(sapply(1:ncol(scale.data),function(x) unlist(lapply(sort(unique(cell.ident)),function(y)dist(rbind(scale.data[,x],means.all[,y])))))))
  cell.ident=apply(all.dist,1,which.min)
  cell.ident=order(tapply(as.numeric(data.use),cell.ident,mean))[cell.ident]
  return(cell.ident)
}

#Internal, not documented for now
#' @export
setGeneric("fit.gene.mix", function(object, gene, do.k=3,use.mixtools=TRUE,do.plot=FALSE,plot.with.imputed=TRUE,min.bin.size=10) standardGeneric("fit.gene.mix"))
#' @export
setMethod("fit.gene.mix", "seurat",
          function(object, gene, do.k=3,use.mixtools=TRUE,do.plot=FALSE,plot.with.imputed=TRUE,min.bin.size=10) {
            data.fit=as.numeric(object@imputed[gene,])
            mixtools.fit=normalmixEM(data.fit,k=do.k)
            comp.order=order(mixtools.fit$mu)
            mixtools.posterior=data.frame(mixtools.fit$posterior[,comp.order])
            colnames(mixtools.posterior)=unlist(lapply(1:do.k,function(x)paste(gene,x-1,"post",sep=".")))
            
            #mixtools.mu=data.frame(mixtools.fit$mu[comp.order])
            #mixtools.sigma=data.frame(mixtools.fit$sigma[comp.order])
            #mixtools.alpha=data.frame(mixtools.fit$lambda[comp.order])
            #rownames(mixtools.mu)=unlist(lapply(1:do.k,function(x)paste(gene,x-1,"mu",sep=".")))
            #rownames(mixtools.sigma)=unlist(lapply(1:do.k,function(x)paste(gene,x-1,"sigma",sep=".")))
            #rownames(mixtools.alpha)=unlist(lapply(1:do.k,function(x)paste(gene,x-1,"alpha",sep=".")))
            #object@mix.mu = rbind(minusr(object@mix.mu,gene), mixtools.mu); 
            #object@mix.sigma = rbind(minusr(object@mix.sigma,gene), mixtools.sigma); 
            #o#bject@mu.alpha =rbind(minusr(object@mu.alpha,gene), mixtools.alpha); 
            
            if (do.plot) {
              nCol=2
              num.row=floor((do.k+1)/nCol-1e-5)+1
              par(mfrow=c(num.row,nCol))
              plot.mixEM(mixtools.fit,which=2)
              plot.data=as.numeric(object@imputed[gene,])
              if (!plot.with.imputed) plot.data=as.numeric(object@data[gene,])
              unlist(lapply(1:do.k,function(x) plot(plot.data,mixtools.posterior[,x],ylab=paste("Posterior for Component ",x-1,sep=""),xlab=gene,main=gene)))
            }
            new.mix.probs=data.frame(minusc(object@mix.probs,paste(gene,".",sep="")),row.names = rownames(object@mix.probs)); colnames(new.mix.probs)[1]="nGene"
            object@mix.probs=cbind(new.mix.probs,mixtools.posterior)
            return(object)
          } 
)

#Internal, not documented for now
lasso.fxn = function(lasso.input,genes.obs,s.use=20,gene.name=NULL,do.print=FALSE,gram=TRUE) {
  lasso.model=lars(lasso.input,as.numeric(genes.obs),type="lasso",max.steps = s.use*2,use.Gram=gram)
  #lasso.fits=predict.lars(lasso.model,lasso.input,type="fit",s=min(s.use,max(lasso.model$df)))$fit
  lasso.fits=predict.lars(lasso.model,lasso.input,type="fit",s=s.use)$fit
  if (do.print) print(gene.name)
  return(lasso.fits)  
}

#' Calculate imputed expression values
#'
#' Uses L1-constrained linear models (LASSO) to impute single cell gene
#' expression values.
#'
#'
#' @param object Seurat object
#' @param genes.use A vector of genes (predictors) that can be used for
#' building the LASSO models.
#' @param genes.fit A vector of genes to impute values for
#' @param s.use Maximum number of steps taken by the algorithm (lower values
#' indicate a greater degree of smoothing)
#' @param do.print Print progress (output the name of each gene after it has
#' been imputed).
#' @param gram The use.gram argument passed to lars
#' @return Returns a Seurat object where the imputed values have been added to
#' object@@data
#' @import lars
#' @export
setGeneric("addImputedScore", function(object, genes.use=NULL,genes.fit=NULL,s.use=20,do.print=FALSE,gram=TRUE) standardGeneric("addImputedScore"))
#' @export
setMethod("addImputedScore", "seurat",
          function(object, genes.use=NULL,genes.fit=NULL,s.use=20,do.print=FALSE,gram=TRUE) {
            genes.use=set.ifnull(genes.use,object@var.genes)
            genes.fit=set.ifnull(genes.fit,object@var.genes)
            genes.use=genes.use[genes.use%in%rownames(object@data)]
            genes.fit=genes.fit[genes.fit%in%rownames(object@data)]
            
            lasso.input=t(object@data[genes.use,])
            lasso.fits=data.frame(t(sapply(genes.fit,function(x)lasso.fxn(t(object@data[genes.use[genes.use!=x],]),object@data[x,],s.use=s.use,x,do.print,gram))))
            genes.old=genes.fit[genes.fit%in%rownames(object@imputed)]
            genes.new=genes.fit[!(genes.fit%in%rownames(object@imputed))]
            
            if (length(genes.old)>0) object@imputed[genes.old,]=lasso.fits[genes.old,]
            object@imputed=rbind(object@imputed,lasso.fits[genes.new,])
            return(object)
          }
)    


# Not currently supported, but a cool scoring function
setGeneric("getNewScore", function(object, score.name,score.genes, cell.ids=NULL, score.func=weighted.mean,scramble=FALSE, no.tech.wt=FALSE, biol.wts=NULL,use.scaled=FALSE) standardGeneric("getNewScore"))
setMethod("getNewScore", "seurat",
          function(object, score.name,score.genes, cell.ids=NULL, score.func=weighted.mean,scramble=FALSE, no.tech.wt=FALSE, biol.wts=NULL,use.scaled=FALSE) {
            data.use=object@data; if (use.scaled) data.use=minmax(object@scale.data,min = -2,max=2)
            if (!(no.tech.wt)) score.genes=score.genes[score.genes%in%rownames(object@wt.matrix)]
            if (no.tech.wt) score.genes=score.genes[score.genes%in%rownames(data.use)]
            cell.ids=set.ifnull(cell.ids,colnames(data.use))
            wt.matrix=data.frame(matrix(1,nrow=length(score.genes),ncol = length(cell.ids),dimnames = list(score.genes,cell.ids))); 
            if (!(no.tech.wt)) wt.matrix=object@drop.wt.matrix[score.genes,cell.ids]
            score.data=data.use[score.genes,cell.ids]
            if(scramble) {
              score.data=score.data[,sample(ncol(score.data))]
            }
            wt.matrix=wt.matrix*(wt.matrix)
            if (no.tech.wt) wt.matrix[wt.matrix<1]=1
            biol.wts=set.ifnull(biol.wts,rep(1,nrow(wt.matrix)))
            if (mean(biol.wts)==1) names(biol.wts)=score.genes
            my.scores=unlist(lapply(colnames(score.data),function(x)score.func(score.data[,x],wt.matrix[,x]*biol.wts[score.genes])))
            names(my.scores)=colnames(score.data)
            object@gene.scores[cell.ids,score.name]=my.scores
            return(object)
          }
)

# Not currently supported, but a cool function for QC
#' @export
setGeneric("calcNoiseModels", function(object, cell.ids=NULL, trusted.genes=NULL,n.bin=20,drop.expr=1) standardGeneric("calcNoiseModels"))
#' @export
setMethod("calcNoiseModels","seurat",
          function(object, cell.ids=NULL, trusted.genes=NULL,n.bin=20,drop.expr=1) {
            object@drop.expr=drop.expr
            cell.ids=set.ifnull(cell.ids,1:ncol(object@data))
            trusted.genes=set.ifnull(trusted.genes,rownames(object@data))
            trusted.genes=trusted.genes[trusted.genes%in%rownames(object@data)]
            object@trusted.genes=trusted.genes
            data=object@data[trusted.genes,]
            idents=data.frame(data[,1])
            code_humpAvg=apply(data,1,humpMean,min=object@drop.expr)
            code_humpAvg[code_humpAvg>9]=9
            code_humpAvg[is.na(code_humpAvg)]=0
            idents$code_humpAvg=code_humpAvg
            data[data>object@drop.expr]=1
            data[data<object@drop.expr]=0
            data$bin=cut(code_humpAvg,n.bin)
            data$avg=code_humpAvg
            rownames(idents)=rownames(data)
            my.coefs=data.frame(t(sapply(colnames(data[1:(ncol(data)-2)]),
                                         getAB,data=data,data2=idents,status="code",code2="humpAvg",hasBin=TRUE,doPlot=FALSE)))
            colnames(my.coefs)=c("a","b")
            object@drop.coefs = my.coefs
            return(object)
          }
)  

#' Visualize 'features' on a dimensional reduction plot
#' 
#' Colors single cells on a dimensional reduction plot according to a 'feature'
#' (i.e. gene expression, PC scores, number of genes detected, etc.)
#' 
#' To determine the color, the feature values across all cells are placed into
#' discrete bins, and then assigned a color based on cols.use. The number of
#' bins is determined by the number of colors in cols.use
#' 
#' @param object Seurat object
#' @param features.plot Vector of features to plot
#' @param dim.1 Dimension for x-axis (default 1)
#' @param dim.2 Dimension for y-axis (default 2)
#' @param cells.use Vector of cells to plot (default is all cells)
#' @param pt.size Adjust point size for plotting
#' @param cols.use Ordered vector of colors to use for plotting. Default is
#' heat.colors(10).
#' @param pch.use Pch for plotting
#' @param reduction.use Which dimensionality reduction to use. Default is
#' "tsne", can also be "pca", or "ica", assuming these are precomputed.
#' @param nCol Number of columns to use when plotting multiple features.
#' @return No return value, only a graphical output
#' @export
setGeneric("feature.plot", function(object,features.plot,dim.1=1,dim.2=2,cells.use=NULL,pt.size=1,cols.use=heat.colors(10),pch.use=16,reduction.use="tsne",nCol=NULL) standardGeneric("feature.plot"))
#' @export
setMethod("feature.plot", "seurat", 
          function(object,features.plot,dim.1=1,dim.2=2,cells.use=NULL,pt.size=1,cols.use=heat.colors(10),pch.use=16,reduction.use="tsne",nCol=NULL) {
            cells.use=set.ifnull(cells.use,colnames(object@data))
            dim.code="PC"
            if (is.null(nCol)) {
              nCol=2
              if (length(features.plot)>6) nCol=3
              if (length(features.plot)>9) nCol=4
            }         
            num.row=floor(length(features.plot)/nCol-1e-5)+1
            par(mfrow=c(num.row,nCol))
            dim.code=translate.dim.code(reduction.use); dim.codes=paste(dim.code,c(dim.1,dim.2),sep="")
            data.plot=fetch.data(object,dim.codes)
      
            x1=paste(dim.code,dim.1,sep=""); x2=paste(dim.code,dim.2,sep="")
            data.plot$x=data.plot[,x1]; data.plot$y=data.plot[,x2]
            data.plot$pt.size=pt.size
            data.use=data.frame(t(fetch.data(object,features.plot,cells.use = cells.use)))
            for(i in features.plot) {
              data.gene=na.omit(data.frame(data.use[i,]))
              data.cut=as.numeric(as.factor(cut(as.numeric(data.gene),breaks = length(cols.use))))
              data.col=rev(cols.use)[data.cut]
              plot(data.plot$x,data.plot$y,col=data.col,cex=pt.size,pch=pch.use,main=i,xlab=x1,ylab=x2)
            }
            rp()
          }
)

#' Vizualization of multiple features
#'
#' Similar to feature.plot, however, also splits the plot by visualizing each
#' identity class separately.
#'
#' Particularly useful for seeing if the same groups of cells co-exhibit a
#' common feature (i.e. co-express a gene), even within an identity class. Best
#' understood by example.
#'
#'
#' @param object Seurat object
#' @param features.plot Vector of features to plot
#' @param dim.1 Dimension for x-axis (default 1)
#' @param dim.2 Dimension for y-axis (default 2)
#' @param idents.use Which identity classes to display (default is all identity
#' classes)
#' @param pt.size Adjust point size for plotting
#' @param cols.use Ordered vector of colors to use for plotting. Default is
#' heat.colors(10).
#' @param pch.use Pch for plotting
#' @param reduction.use Which dimensionality reduction to use. Default is
#' "tsne", can also be "pca", or "ica", assuming these are precomputed.
#' @return No return value, only a graphical output
#' @export
setGeneric("feature.heatmap", function(object,features.plot,dim.1=1,dim.2=2,idents.use=NULL,pt.size=2,cols.use=rev(heat.colors(10)),pch.use=16,reduction.use="tsne") standardGeneric("feature.heatmap"))
#' @export
setMethod("feature.heatmap", "seurat", 
          function(object,features.plot,dim.1=1,dim.2=2,idents.use=NULL,pt.size=2,cols.use=rev(heat.colors(10)),pch.use=16,reduction.use="tsne") {
            idents.use=set.ifnull(idents.use,sort(unique(object@ident)))
            dim.code="PC"
            par(mfrow=c(length(features.plot),length(idents.use)))
            dim.code=translate.dim.code(reduction.use); dim.codes=paste(dim.code,c(dim.1,dim.2),sep="")
            data.plot=data.frame(fetch.data(object,dim.codes))
            
            ident.use=as.factor(object@ident)
            data.plot$ident=ident.use
            x1=paste(dim.code,dim.1,sep=""); x2=paste(dim.code,dim.2,sep="")
            data.plot$x=data.plot[,x1]; data.plot$y=data.plot[,x2]
            data.plot$pt.size=pt.size
            data.use=data.frame(t(data.frame(fetch.data(object,features.plot))))
            data.scale=apply(t(data.use),2,function(x)(factor(cut(x,breaks=length(cols.use),labels = FALSE),levels=1:length(cols.use),ordered=TRUE)))
            data.plot.all=data.frame(cbind(data.plot,data.scale))
            data.reshape=melt((data.plot.all),id = colnames(data.plot))
            data.reshape=data.reshape[data.reshape$ident%in%idents.use,]
            data.reshape$value=factor(data.reshape$value,levels=1:length(cols.use),ordered=TRUE)
            #p <- ggplot(data.reshape, aes(x,y)) + geom_point(aes(colour=reorder(value,1:length(cols.use)),size=pt.size)) + scale_colour_manual(values=cols.use)
            p <- ggplot(data.reshape, aes(x,y)) + geom_point(aes(colour=value,size=pt.size)) + scale_colour_manual(values=cols.use)
            
            p=p + facet_grid(variable~ident) + scale_size(range = c(pt.size, pt.size))
            p2=p+gg.xax()+gg.yax()+gg.legend.pts(6)+gg.legend.text(12)+no.legend.title+theme_bw()+nogrid+theme(legend.title=element_blank())
            print(p2)
          }
)

#' Plot tSNE map
#'
#' Graphs the output of a tSNE analysis
#' Cells are colored by their identity class.
#'
#' This function is a wrapper for dim.plot. See ?dim.plot for a full list of possible
#' arguments which can be passed in here.
#' 
#' @param object Seurat object
#' @param do.label FALSE by default. If TRUE, plots an alternate view where the center of each
#' cluster is lebeled
#' @param label.pt.size If do.label is set, the point size
#' @param label.cex.text If label.cex.text is set, the size of the text labels
#' @param label.cols.use If do.label is set, the color palette to use for the points
#' @param \dots Additional parameters to dim.plot, for example, which dimensions to plot. 
#' @seealso dim.plot
#' @export
setGeneric("tsne.plot", function(object,do.label=FALSE,label.pt.size=1,label.cex.text=1,label.cols.use=NULL,...) standardGeneric("tsne.plot"))
#' @export
setMethod("tsne.plot", "seurat", 
          function(object,do.label=FALSE,label.pt.size=1,label.cex.text=1,label.cols.use=NULL,...) {
            if (do.label==TRUE) {
              label.cols.use=set.ifnull(label.cols.use,rainbow(length(levels(object@ident)))); 
              plot(object@tsne.rot[,1],object@tsne.rot[,2],col=label.cols.use[as.integer(object@ident)],pch=16,xlab="tSNE_1",ylab="tSNE_2",cex=label.pt.size)
              k.centers=t(sapply(levels(object@ident),function(x) apply(object@tsne.rot[which.cells(object,x),],2,mean)))
              points(k.centers[,1],k.centers[,2],cex=1.3,col="white",pch=16); text(k.centers[,1],k.centers[,2],levels(object@ident),cex=label.cex.text)
            }
            else {
              return(dim.plot(object,reduction.use = "tsne",...))
            }
          }
)

#' Plot ICA map
#'
#' Graphs the output of a ICA analysis
#' Cells are colored by their identity class.
#'
#' This function is a wrapper for dim.plot. See ?dim.plot for a full list of possible
#' arguments which can be passed in here.
#' 
#' @param object Seurat object
#' @param \dots Additional parameters to dim.plot, for example, which dimensions to plot. 
#' @export
setGeneric("ica.plot", function(object,...) standardGeneric("ica.plot"))
#' @export
setMethod("ica.plot", "seurat", 
          function(object,...) {
            return(dim.plot(object,reduction.use = "ica",...))
          }
)

#' Plot PCA map
#'
#' Graphs the output of a PCA analysis
#' Cells are colored by their identity class.
#'
#' This function is a wrapper for dim.plot. See ?dim.plot for a full list of possible
#' arguments which can be passed in here.
#' 
#' @param object Seurat object
#' @param \dots Additional parameters to dim.plot, for example, which dimensions to plot. 
#' @export
setGeneric("pca.plot", function(object,...) standardGeneric("pca.plot"))
#' @export
setMethod("pca.plot", "seurat", 
          function(object,...) {
              return(dim.plot(object,reduction.use = "pca",...))
          }
)

translate.dim.code=function(reduction.use) {
  return.code="PC"
  if (reduction.use=="ica") return.code="IC"
  if (reduction.use=="tsne") return.code="tSNE_"
  if (reduction.use=="mds") return.code="MDS"
  return(return.code)
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
#' cell attribute (that can be pulled with fetch.data) allowing for both different colors and different shapes on cells.
#' @return If do.return==TRUE, returns a ggplot2 object. Otherwise, only
#' graphical output.
#' @export
setGeneric("dim.plot", function(object,reduction.use="pca",dim.1=1,dim.2=2,cells.use=NULL,pt.size=3,do.return=FALSE,do.bare=FALSE,cols.use=NULL,group.by="ident",pt.shape=NULL) standardGeneric("dim.plot"))
#' @export
setMethod("dim.plot", "seurat", 
          function(object,reduction.use="pca",dim.1=1,dim.2=2,cells.use=NULL,pt.size=3,do.return=FALSE,do.bare=FALSE,cols.use=NULL,group.by="ident",pt.shape=NULL) {
            cells.use=set.ifnull(cells.use,colnames(object@data))
            dim.code=translate.dim.code(reduction.use); dim.codes=paste(dim.code,c(dim.1,dim.2),sep="")
            data.plot=fetch.data(object,dim.codes,cells.use)
            
            ident.use=as.factor(object@ident[cells.use])
            if (group.by != "ident") ident.use=as.factor(fetch.data(object,group.by)[,1])
            data.plot$ident=ident.use
            x1=paste(dim.code,dim.1,sep=""); x2=paste(dim.code,dim.2,sep="")
            data.plot$x=data.plot[,x1]; data.plot$y=data.plot[,x2]
            data.plot$pt.size=pt.size
            p=ggplot(data.plot,aes(x=x,y=y))+geom_point(aes(colour=factor(ident)),size=pt.size)
            if (!is.null(pt.shape)) {
              shape.val=fetch.data(object,pt.shape)[cells.use,1]
              if (is.numeric(shape.val)) {
                shape.val=cut(shape.val,breaks = 5)
              }
              data.plot[,"pt.shape"]=shape.val
              p=ggplot(data.plot,aes(x=x,y=y))+geom_point(aes(colour=factor(ident),shape=factor(pt.shape)),size=pt.size)
              
            }
            if (!is.null(cols.use)) {
              p=p+scale_colour_manual(values=cols.use)
            }
            p2=p+xlab(x1)+ylab(x2)+scale_size(range = c(pt.size, pt.size))
            p3=p2+gg.xax()+gg.yax()+gg.legend.pts(6)+gg.legend.text(12)+no.legend.title+theme_bw()+nogrid
            p3=p3+theme(legend.title=element_blank())
            if (do.return) {
              if (do.bare) return(p)
              return(p3)
            }
            if (do.bare) print(p)
            else print(p3)
          }
)

#Cool, but not supported right now
setGeneric("spatial.de", function(object,marker.cells,genes.use=NULL,...) standardGeneric("spatial.de"))
setMethod("spatial.de", "seurat", 
          function(object,marker.cells,genes.use=NULL) {
            object=p15
            embed.map=object@tsne.rot
            mult.use=2
            mult.use.far=10
            if ((mult.use.far*length(marker.cells))>nrow(embed.map)) {
              mult.use.far=1
              mult.use=1
            }
            genes.use=set.ifnull(genes.use,object@var.genes)
            marker.pos=apply(embed.map[marker.cells,],2,mean)
            embed.map=rbind(embed.map,marker.pos)
            rownames(embed.map)[nrow(embed.map)]="marker"
            embed.dist=sort(as.matrix(dist((embed.map)))["marker",])
            embed.diff=names(embed.dist[!(names(embed.dist)%in%marker.cells)][1:(mult.use*length(marker.cells))][-1])
            embed.diff.far=names(embed.dist[!(names(embed.dist)%in%marker.cells)][1:(mult.use.far*length(marker.cells))][-1])
            
            diff.genes=rownames(subset(diffExp.test(p15,marker.cells,embed.diff,genes.use=genes.use),p_val<(1e-5)))
            diff.genes=subset(diffExp.test(p15,marker.cells,embed.diff,genes.use = diff.genes),p_val<(1e-10))
            return(diff.genes)
          }
)

#' Perform spectral density clustering on single cells
#'
#' Find point clounds single cells in a two-dimensional space using density clustering (DBSCAN). 
#'
#'
#' @param object Seurat object
#' @param dim.1 First dimension to use
#' @param dim.2 second dimension to use
#' @param reduction.use Which dimensional reduction to use (either 'pca' or 'ica')
#' @param G.use Parameter for the density clustering. Lower value to get more fine-scale clustering
#' @param set.ident TRUE by default. Set identity class to the results of the density clustering. 
#' Unassigned cells (cells that cannot be assigned a cluster) are placed in cluster 1, if there are any.
#' @param seed.use Random seed for the dbscan function
#' @param \dots Additional arguments to be passed to the dbscan function
#' @export
setGeneric("DBclust_dimension", function(object,dim.1=1,dim.2=2,reduction.use="tsne",G.use=NULL,set.ident=TRUE,seed.use=1,...) standardGeneric("DBclust_dimension"))
#' @export
setMethod("DBclust_dimension", "seurat", 
          function(object,dim.1=1,dim.2=2,reduction.use="tsne",G.use=NULL,set.ident=TRUE,seed.use=1,...) {
            dim.code=translate.dim.code(reduction.use); dim.codes=paste(dim.code,c(dim.1,dim.2),sep="")
            data.plot=fetch.data(object,dim.codes)
            x1=paste(dim.code,dim.1,sep=""); x2=paste(dim.code,dim.2,sep="")
            data.plot$x=data.plot[,x1]; data.plot$y=data.plot[,x2]
            set.seed(seed.use); data.mclust=ds <- dbscan(data.plot[,c("x","y")],eps = G.use,...)
            
            to.set=as.numeric(data.mclust$cluster+1)
            data.names=names(object@ident)
            object@data.info[data.names,"DBclust.ident"]=to.set
            if (set.ident) {
              object@ident=factor(to.set); names(object@ident)=data.names;               
            }
            
            return(object)
          }
)

#' Perform spectral k-means clustering on single cells
#'
#' Find point clounds single cells in a low-dimensional space using k-means clustering. 
#'
#' CAn be useful for smaller datasets, not documented here yet
#' @export
setGeneric("Kclust_dimension", function(object,dim.1=1,dim.2=2,cells.use=NULL,pt.size=4,reduction.use="tsne",k.use=5,set.ident=FALSE,seed.use=1,...) standardGeneric("Kclust_dimension"))
#' @export
setMethod("Kclust_dimension", "seurat", 
          function(object,dim.1=1,dim.2=2,cells.use=NULL,pt.size=4,reduction.use="tsne",k.use=5,set.ident=FALSE,seed.use=1,...) {
            cells.use=set.ifnull(cells.use,colnames(object@data))
            dim.code="PC"
            if (reduction.use=="pca") data.plot=object@pca.rot[cells.use,]
            if (reduction.use=="tsne") {
              data.plot=object@tsne.rot[cells.use,]
              dim.code="Seurat_"
            }
            if (reduction.use=="ica") {
              data.plot=object@ica.rot[cells.use,]
              dim.code="IC"
            }
            x1=paste(dim.code,dim.1,sep=""); x2=paste(dim.code,dim.2,sep="")
            data.plot$x=data.plot[,x1]; data.plot$y=data.plot[,x2]
            if (reduction.use!="pca") {
              set.seed(seed.use); data.mclust=ds <- kmeans(data.plot[,c("x","y")], k.use)   
            }
            if (reduction.use=="pca") {
              set.seed(seed.use); data.mclust=ds <- kmeans(object@pca.rot[cells.use,dim.1], k.use)   
            }
            to.set=as.numeric(data.mclust$cluster)
            data.names=names(object@ident)
            object@data.info[data.names,"kdimension.ident"]=to.set
            if (set.ident) {
              object@ident=factor(to.set); names(object@ident)=data.names;               
            }
            
            return(object)
          }
)

#' Significant genes from a PCA
#'
#' Returns a set of genes, based on the jackStraw analysis, that have
#' statistically significant associations with a set of PCs.
#'
#'
#' @param object Seurat object
#' @param pcs.use PCS to use.
#' @param pval.cut P-value cutoff
#' @param use.full Use the full list of genes (from the projected PCA). Assumes
#' that project.pca has been run. Default is TRUE
#' @param max.per.pc Maximum number of genes to return per PC. Used to avoid genes from one PC dominating the entire analysis.
#' @return A vector of genes whose p-values are statistically significant for
#' at least one of the given PCs.
#' @export
setGeneric("pca.sig.genes", function(object,pcs.use,pval.cut=0.1,use.full=TRUE,max.per.pc=NULL) standardGeneric("pca.sig.genes"))
#' @export
setMethod("pca.sig.genes", "seurat", 
          function(object,pcs.use,pval.cut=0.1,use.full=TRUE,max.per.pc=NULL) {
            pvals.use=object@jackStraw.empP
            pcx.use=object@pca.x
            if (use.full)  {
              pvals.use=object@jackStraw.empP.full
              pcx.use=object@pca.x.full
            }
            if (length(pcs.use)==1) pvals.min=pvals.use[,pcs.use]
            if (length(pcs.use)>1) pvals.min=apply(pvals.use[,pcs.use],1,min)
            names(pvals.min)=rownames(pvals.use)
            genes.use=names(pvals.min)[pvals.min<pval.cut]
            if (!is.null(max.per.pc)) {
              pc.top.genes=pcTopGenes(object,pcs.use,max.per.pc,use.full,FALSE)
              genes.use=ainb(pc.top.genes,genes.use)
            }
            return(genes.use)       
          }
)



same=function(x) return(x)

#' Gene expression heatmap
#'
#' Draws a heatmap of single cell gene expression using the heatmap.2 function.
#' 
#' @param object Seurat object
#' @param cells.use Cells to include in the heatmap (default is all cells)
#' @param genes.use Genes to include in the heatmap (ordered)
#' @param disp.min Minimum display value (all values below are clipped)
#' @param disp.max Maximum display value (all values above are clipped)
#' @param draw.line Draw vertical lines delineating cells in different identity
#' classes.
#' @param do.return Default is FALSE. If TRUE, return a matrix of scaled values
#' which would be passed to heatmap.2
#' @param order.by.ident Order cells in the heatmap by identity class (default
#' is TRUE). If FALSE, cells are ordered based on their order in cells.use
#' @param col.use Color palette to use
#' @param slim.col.label if (order.by.ident==TRUE) then instead of displaying 
#' every cell name on the heatmap, display only the identity class name once 
#' for each group
#' @param group.by If (order.by.ident==TRUE) default,  you can group cells in 
#' different ways (for example, orig.ident)
#' @param remove.key Removes the color key from the plot.
#' @param \dots Additional parameters to heatmap.2. Common examples are cexRow
#' and cexCol, which set row and column text sizes
#' @return If do.return==TRUE, a matrix of scaled values which would be passed
#' to heatmap.2. Otherwise, no return value, only a graphical output
#' @importFrom gplots heatmap.2
#' @export
setGeneric("doHeatMap", function(object,cells.use=NULL,genes.use=NULL,disp.min=-2.5,disp.max=2.5,draw.line=TRUE,do.return=FALSE,order.by.ident=TRUE,col.use=pyCols,slim.col.label=FALSE,group.by=NULL,remove.key=FALSE,...) standardGeneric("doHeatMap"))
#' @export
setMethod("doHeatMap","seurat",
          function(object,cells.use=NULL,genes.use=NULL,disp.min=-2.5,disp.max=2.5,draw.line=TRUE,do.return=FALSE,order.by.ident=TRUE,col.use=pyCols,slim.col.label=FALSE,group.by=NULL,remove.key=FALSE,...) {
            cells.use=set.ifnull(cells.use,object@cell.names)
            genes.use=ainb(genes.use,rownames(object@scale.data))
            cells.use=ainb(cells.use,object@cell.names)
            cells.ident=object@ident[cells.use]
            if (!is.null(group.by)) cells.ident=factor(fetch.data(object,group.by)[,1])
            cells.ident=factor(cells.ident,labels = ainb(levels(cells.ident),cells.ident))
            if (order.by.ident) {
              cells.use=cells.use[order(cells.ident)]
            }
            data.use=object@scale.data[genes.use,cells.use]
            data.use=minmax(data.use,min=disp.min,max=disp.max)
            vline.use=NULL;
            colsep.use=NULL
            hmFunction=heatmap.2
            if (remove.key) hmFunction=heatmap2NoKey
            if (draw.line) {
              colsep.use=cumsum(table(cells.ident))
            }
            if(slim.col.label && order.by.ident) {
              col.lab=rep("",length(cells.use))
              col.lab[round(cumsum(table(cells.ident))-table(cells.ident)/2)+1]=levels(cells.ident)
              hmFunction(data.use,Rowv=NA,Colv=NA,trace = "none",col=col.use,colsep = colsep.use,labCol=col.lab,cexCol=0.2+1/log10(length(unique(cells.ident))),...)
            }
            else {
              hmFunction(data.use,Rowv=NA,Colv=NA,trace = "none",col=col.use,colsep = colsep.use,...)
            }
            if (do.return) {
              return(data.use)
            }
          }
)

#' Independent component heatmap
#' 
#' Draws a heatmap focusing on a principal component. Both cells and genes are sorted by their principal component scores. Allows for nice visualization of sources of heterogeneity in the dataset.
#' 
#' @inheritParams doHeatMap
#' @inheritParams icTopGenes
#' @inheritParams viz.ica
#' @param use.scale Default is TRUE: plot scaled data. If FALSE, plot raw data on the heatmap.
#' @return If do.return==TRUE, a matrix of scaled values which would be passed
#' to heatmap.2. Otherwise, no return value, only a graphical output
#' @export
setGeneric("icHeatmap", function(object,ic.use=1,cells.use=NULL,num.genes=30, disp.min=-2.5,disp.max=2.5,do.return=FALSE,col.use=pyCols,use.scale=TRUE,do.balanced=FALSE,...) standardGeneric("icHeatmap"))
#' @export
setMethod("icHeatmap","seurat",
          function(object,ic.use=1,cells.use=NULL,num.genes=30,disp.min=-2.5,disp.max=2.5,do.return=FALSE,col.use=pyCols,use.scale=TRUE,do.balanced=FALSE,...) {
            cells.use=set.ifnull(cells.use,object@cell.names)
            genes.use=icTopGenes(object,ic.use,num.genes,do.balanced)
            cells.ordered=cells.use[order(object@ica.rot[cells.use,ic.use])]
            data.use=object@scale.data[genes.use,cells.ordered]
            data.use=minmax(data.use,min=disp.min,max=disp.max)
            if (!(use.scale)) data.use=as.matrix(object@data[genes.use,cells.ordered])
            vline.use=NULL;
            heatmap.2(data.use,Rowv=NA,Colv=NA,trace = "none",col=col.use)
            if (do.return) {
              return(data.use)
            }
          }
)

#' Principal component heatmap
#' 
#' Draws a heatmap focusing on a principal component. Both cells and genes are sorted by their principal component scores. Allows for nice visualization of sources of heterogeneity in the dataset.
#' 
#' @inheritParams doHeatMap
#' @inheritParams pcTopGenes
#' @inheritParams viz.pca
#' @param use.scale Default is TRUE: plot scaled data. If FALSE, plot raw data on the heatmap.
#' @return If do.return==TRUE, a matrix of scaled values which would be passed
#' to heatmap.2. Otherwise, no return value, only a graphical output
#' @export
setGeneric("pcHeatmap", function(object,pc.use=1,cells.use=NULL,num.genes=30,use.full=FALSE, disp.min=-2.5,disp.max=2.5,do.return=FALSE,col.use=pyCols,use.scale=TRUE,do.balanced=FALSE,remove.key=FALSE,...) standardGeneric("pcHeatmap"))
#' @export
setMethod("pcHeatmap","seurat",
          function(object,pc.use=1,cells.use=NULL,num.genes=30,use.full=FALSE, disp.min=-2.5,disp.max=2.5,do.return=FALSE,col.use=pyCols,use.scale=TRUE,do.balanced=FALSE,remove.key=FALSE,...) {
            cells.use=set.ifnull(cells.use,object@cell.names)
            genes.use=rev(pcTopGenes(object,pc.use,num.genes,use.full,do.balanced))
            cells.ordered=cells.use[order(object@pca.rot[cells.use,pc.use])]
            data.use=object@scale.data[genes.use,cells.ordered]
            data.use=minmax(data.use,min=disp.min,max=disp.max)
            if (!(use.scale)) data.use=as.matrix(object@data[genes.use,cells.ordered])
            vline.use=NULL;
            hmFunction=heatmap.2; if (remove.key) hmFunction=heatmap2NoKey
            hmFunction(data.use,Rowv=NA,Colv=NA,trace = "none",col=col.use,...)
            if (do.return) {
              return(data.use)
            }
          }
)


#' K-Means Clustering
#' 
#' Perform k=means clustering on both genes and single cells
#' 
#' K-means and heatmap are calculated on object@@scale.data
#' 
#' @param object Seurat object
#' @param genes.use Genes to use for clustering
#' @param k.genes K value to use for clustering genes
#' @param k.cells K value to use for clustering cells (default is NULL, cells
#' are not clustered)
#' @param k.seed Random seed
#' @param do.plot Draw heatmap of clustered genes/cells (default is TRUE)
#' @param data.cut Clip all z-scores to have an absolute value below this.
#' Reduces the effect of huge outliers in the data.
#' @param k.cols Color palette for heatmap
#' @param pc.row.order Order gene clusters based on the average PC score within
#' a cluster. Can be useful if you want to visualize clusters, for example,
#' based on their average score for PC1.
#' @param pc.col.order Order cell clusters based on the average PC score within
#' a cluster
#' @param rev.pc.order Use the reverse PC ordering for gene and cell clusters
#' (since the sign of a PC is arbitrary)
#' @param use.imputed Cluster imputed values (default is FALSE)
#' @param set.ident If clustering cells (so k.cells>0), set the cell identity
#' class to its K-means cluster (default is TRUE)
#' @param \dots Additional parameters passed to doHeatMap for plotting
#' @return Seurat object where the k-means results for genes is stored in
#' object@@kmeans.obj[[1]], and the k-means results for cells is stored in
#' object@@kmeans.col[[1]]. The cluster for each cell is stored in object@@data.info[,"kmeans.ident"]
#' and also object@@ident (if set.ident=TRUE)
#' @export
setGeneric("doKMeans", function(object,genes.use=NULL,k.genes=NULL,k.cells=NULL,k.seed=1,do.plot=TRUE,data.cut=2.5,k.cols=pyCols,
                                pc.row.order=NULL,pc.col.order=NULL, rev.pc.order=FALSE, use.imputed=FALSE,set.ident=TRUE,...) standardGeneric("doKMeans"))
#' @export
setMethod("doKMeans","seurat",
          function(object,genes.use=NULL,k.genes=NULL,k.cells=0,k.seed=1,do.plot=TRUE,data.cut=2.5,k.cols=pyCols,
                   pc.row.order=NULL,pc.col.order=NULL, rev.pc.order=FALSE, use.imputed=FALSE,set.ident=TRUE,...) {
            
            data.use.orig=object@scale.data
            if (use.imputed) data.use.orig=data.frame(t(scale(t(object@imputed))))
            data.use=minmax(data.use.orig,min=data.cut*(-1),max=data.cut)
            revFxn=same; if (rev.pc.order) revFxn=function(x)max(x)+1-x;
            kmeans.col=NULL

            genes.use=set.ifnull(genes.use,object@var.genes)          
            genes.use=genes.use[genes.use%in%rownames(data.use)]
            cells.use=object@cell.names
            
            kmeans.data=data.use[genes.use,cells.use]      
            set.seed(k.seed); kmeans.obj=kmeans(kmeans.data,k.genes); 
            if (!(is.null(pc.row.order))) {
              pcx.use=object@pca.x; pc.genes=ainb(genes.use,rownames(pcx.use))
              if(nrow(object@pca.x.full>0)) {
                pcx.use=object@pca.x.full; pc.genes=ainb(genes.use,rownames(pcx.use))
              }
              kmeans.obj$cluster=as.numeric(revFxn(rank(tapply(pcx.use[genes.use,pc.row.order],as.numeric(kmeans.obj$cluster),mean)))[as.numeric(kmeans.obj$cluster)])
            }
            names(kmeans.obj$cluster)=genes.use
            
            if (k.cells>0) {
              kmeans.col=kmeans(t(kmeans.data),k.cells)
              if (!(is.null(pc.col.order))) {
                    kmeans.col$cluster=as.numeric(revFxn(rank(tapply(object@pca.rot[cells.use,pc.col.order],kmeans.col$cluster,mean)))[as.numeric(kmeans.col$cluster)])
              }
              names(kmeans.col$cluster)=cells.use
            }
            
            object@kmeans.obj=list(kmeans.obj)
            object@kmeans.col=list(kmeans.col)
          
            kmeans.obj=object@kmeans.obj[[1]]
            kmeans.col=object@kmeans.col[[1]]
            object@data.info[names(kmeans.col$cluster),"kmeans.ident"]=kmeans.col$cluster
            
            if ((set.ident) && (k.cells > 0)) {
              object=set.ident(object,cells.use=names(kmeans.col$cluster),ident.use = kmeans.col$cluster)
            }
            if (do.plot) {       
              disp.data=minmax(kmeans.data[order(kmeans.obj$cluster[genes.use]),],min=data.cut*(-1),max=data.cut)
              doHeatMap(object,object@cell.names,names(sort(kmeans.obj$cluster)),data.cut*(-1),data.cut,col.use = k.cols,...)
            }
            return(object)
          }
)

#' @export
setGeneric("genes.in.cluster", function(object, cluster.num)  standardGeneric("genes.in.cluster"))
#' @export
setMethod("genes.in.cluster", signature = "seurat",
          function(object, cluster.num) {
            print(unlist(lapply(cluster.num,function(x)sort(names(which(object@kmeans.obj[[1]]$cluster==x))))))
          }    
)



#' @export
setGeneric("cell.cor.matrix", function(object, cor.genes=NULL,cell.inds=NULL, do.k=FALSE,k.seed=1,k.num=4,vis.low=(-1),vis.high=1,vis.one=0.8,pcs.use=1:3,col.use=pyCols)  standardGeneric("cell.cor.matrix"))
#' @export
setMethod("cell.cor.matrix", signature = "seurat",
          function(object, cor.genes=NULL,cell.inds=NULL, do.k=FALSE,k.seed=1,k.num=4,vis.low=(-1),vis.high=1,vis.one=0.8,pcs.use=1:3,col.use=pyCols) {
            cor.genes=set.ifnull(cor.genes,object@var.genes)
            cell.inds=set.ifnull(cell.inds,colnames(object@data))
            cor.genes=cor.genes[cor.genes%in%rownames(object@data)]
            data.cor=object@data[cor.genes,cell.inds]
            cor.matrix=cor((data.cor))
            set.seed(k.seed); kmeans.cor=kmeans(cor.matrix,k.num)
            cor.matrix=cor.matrix[order(kmeans.cor$cluster),order(kmeans.cor$cluster)]
            kmeans.names=rownames(cor.matrix)
            row.annot=data.frame(cbind(kmeans.cor$cluster[kmeans.names],object@pca.rot[kmeans.names,pcs.use]))
            colnames(row.annot)=c("K",paste("PC",pcs.use,sep=""))
            cor.matrix[cor.matrix==1]=vis.one
            cor.matrix=minmax(cor.matrix,min = vis.low,max=vis.high)
            object@kmeans.cell=list(kmeans.cor)
            if (do.k) aheatmap(cor.matrix,col=col.use,Rowv=NA,Colv=NA,annRow=row.annot)
            if (!(do.k)) aheatmap(cor.matrix,col=col.use,annRow=row.annot)
            return(object)
          }    
)

#' @export
setGeneric("gene.cor.matrix", function(object, cor.genes=NULL,cell.inds=NULL, do.k=FALSE,k.seed=1,k.num=4,vis.low=(-1),vis.high=1,vis.one=0.8,pcs.use=1:3,col.use=pyCols)  standardGeneric("gene.cor.matrix"))
#' @export
setMethod("gene.cor.matrix", signature = "seurat",
          function(object, cor.genes=NULL,cell.inds=NULL, do.k=FALSE,k.seed=1,k.num=4,vis.low=(-1),vis.high=1,vis.one=0.8,pcs.use=1:3,col.use=pyCols) {
            cor.genes=set.ifnull(cor.genes,object@var.genes)
            cell.inds=set.ifnull(cell.inds,colnames(object@data))
            cor.genes=cor.genes[cor.genes%in%rownames(object@data)]
            data.cor=object@data[cor.genes,cell.inds]
            cor.matrix=cor(t(data.cor))
            set.seed(k.seed); kmeans.cor=kmeans(cor.matrix,k.num)
            cor.matrix=cor.matrix[order(kmeans.cor$cluster),order(kmeans.cor$cluster)]
            kmeans.names=rownames(cor.matrix)
            row.annot=data.frame(cbind(kmeans.cor$cluster[kmeans.names],object@pca.x[kmeans.names,pcs.use]))
            colnames(row.annot)=c("K",paste("PC",pcs.use,sep=""))
            cor.matrix[cor.matrix==1]=vis.one
            cor.matrix=minmax(cor.matrix,min = vis.low,max=vis.high)
            object@kmeans.gene=list(kmeans.cor)
            if (do.k) aheatmap(cor.matrix,col=col.use,Rowv=NA,Colv=NA,annRow=row.annot)
            if (!(do.k)) aheatmap(cor.matrix,col=col.use,annRow=row.annot)
            return(object)
          }    
)

#' @export
setGeneric("calinskiPlot", function(object,pcs.use,pval.cut=0.1,gene.max=15,col.max=25,use.full=TRUE)  standardGeneric("calinskiPlot"))
#' @export
setMethod("calinskiPlot","seurat",
          function(object,pcs.use,pval.cut=0.1,gene.max=15,col.max=25,use.full=TRUE) {
            if (length(pcs.use)==1) pvals.min=object@jackStraw.empP.full[,pcs.use]
            if (length(pcs.use)>1) pvals.min=apply(object@jackStraw.empP.full[,pcs.use],1,min)
            names(pvals.min)=rownames(object@jackStraw.empP.full)
            genes.use=names(pvals.min)[pvals.min<pval.cut]
            genes.use=genes.use[genes.use%in%rownames(object@scale.data)]
            
            par(mfrow=c(1,2))
            
            mydata <- object@scale.data[genes.use,]
            wss <- (nrow(mydata)-1)*sum(apply(mydata,2,var))
            for (i in 1:gene.max) wss[i] <- sum(kmeans(mydata,
                                                       centers=i)$withinss)
            plot(1:gene.max, wss, type="b", xlab="Number of Clusters for Genes",
                 ylab="Within groups sum of squares")
            
            
            mydata <- t(object@scale.data[genes.use,])
            wss <- (nrow(mydata)-1)*sum(apply(mydata,2,var))
            for (i in 1:col.max) wss[i] <- sum(kmeans(mydata,
                                                      centers=i)$withinss)
            plot(1:col.max, wss, type="b", xlab="Number of Clusters for Cells",
                 ylab="Within groups sum of squares")
            rp()
            return(object)
          }
)

#' @export
setMethod("show", "seurat",
          function(object) {
            cat("An object of class ", class(object), " in project ", object@project.name, "\n", sep = "")
            cat(" ", nrow(object@data), " genes across ",
                ncol(object@data), " samples.\n", sep = "")
            invisible(NULL)
          }
)

plot.Vln=function(gene,data,cell.ident,ylab.max=12,do.ret=FALSE,do.sort=FALSE,size.x.use=16,size.y.use=16,size.title.use=20,adjust.use=1,size.use=1,cols.use=NULL) {
  data$gene=as.character(rownames(data))
  data.use=data.frame(data[gene,])
  if (length(gene)==1) {
    data.melt=data.frame(rep(gene,length(cell.ident))); colnames(data.melt)[1]="gene"
    data.melt$value=as.numeric(data[1,1:length(cell.ident)])
    data.melt$id=names(data)[1:length(cell.ident)]
  }
  #print(head(data.melt))
  
  if (length(gene)>1) data.melt=melt(data.use,id="gene")
  data.melt$ident=cell.ident
  
  noise <- rnorm(length(data.melt$value))/100000
  data.melt$value=as.numeric(as.character(data.melt$value))+noise
  if(do.sort) {
    data.melt$ident=factor(data.melt$ident,levels=names(rev(sort(tapply(data.melt$value,data.melt$ident,mean)))))
  }
  p=ggplot(data.melt,aes(factor(ident),value))
  p2=p + geom_violin(scale="width",adjust=adjust.use,trim=TRUE,aes(fill=factor(ident))) + ylab("Expression level (log TPM)")
  if (!is.null(cols.use)) {
    p2=p2+scale_fill_manual(values=cols.use)
  }
  p3=p2+theme(legend.position="top")+guides(fill=guide_legend(title=NULL))+geom_jitter(height=0,size=size.use)+xlab("Cell Type")
  p4=p3+ theme(axis.title.x = element_text(face="bold", colour="#990000", size=size.x.use), axis.text.x  = element_text(angle=90, vjust=0.5, size=12))+theme_bw()+nogrid
  p5=(p4+theme(axis.title.y = element_text(face="bold", colour="#990000", size=size.y.use), axis.text.y  = element_text(angle=90, vjust=0.5, size=12))+ggtitle(gene)+theme(plot.title = element_text(size=size.title.use, face="bold")))
  if(do.ret==TRUE) {
    return(p5)
  }
  else {
    print(p5)
  }
}

#' Dot plot visualization
#' 
#' Intuitive way of visualizing how gene expression changes across different identity classes (clusters). 
#' The size of the dot encodes the percentage of cells within a class, while the color encodes the 
#' average expression level of 'expressing' cells (green is high).
#' 
#' @param genes.plot Input vector of genes
#' @param cex.use Scaling factor for the dots (scales all dot sizes)
#' @param thresh.col The raw data value which corresponds to a red dot (lowest expression)
#' @param dot.min The fraction of cells at which to draw the smallest dot (default is 0.05)
#' @inheritParams vlnPlot
#' @return Only graphical output
#' @export
setGeneric("dot.plot", function(object,genes.plot,cex.use=2,cols.use=NULL,thresh.col=2.5,dot.min=0.05,group.by=NULL)  standardGeneric("dot.plot"))
#' @export
setMethod("dot.plot","seurat",
          function(object,genes.plot,cex.use=2,cols.use=NULL,thresh.col=2.5,dot.min=0.05,group.by=NULL) {
            genes.plot=ainb(genes.plot,rownames(object@data))
            if (!(is.null(group.by))) object=set.all.ident(object,id = group.by)
            object@data=object@data[genes.plot,]
            avg.exp=average.expression(object)
            avg.alpha=cluster.alpha(object)
            cols.use=set.ifnull(cols.use,myPalette(low = "red",high="green"))
            exp.scale=t(scale(t(avg.exp)))
            exp.scale=minmax(exp.scale,max=thresh.col,min=(-1)*thresh.col)
            n.col=length(cols.use)
            data.y=rep(1:ncol(avg.exp),nrow(avg.exp))
            data.x=unlist(lapply(1:nrow(avg.exp),rep,ncol(avg.exp)))
            data.avg=unlist(lapply(1:length(data.y),function(x) exp.scale[data.x[x],data.y[x]]))
            exp.col=cols.use[floor(n.col*(data.avg+thresh.col)/(2*thresh.col)+.5)]
            data.cex=unlist(lapply(1:length(data.y),function(x) avg.alpha[data.x[x],data.y[x]]))*cex.use+dot.min
            plot(data.x,data.y,cex=data.cex,pch=16,col=exp.col,xaxt="n",xlab="",ylab="",yaxt="n")
            axis(1,at = 1:length(genes.plot),genes.plot)
            axis(2,at=1:ncol(avg.alpha),colnames(avg.alpha),las=1)
          }
          
) 

#' Single cell violin plot
#' 
#' Draws a violin plot of single cell data (gene expression, metrics, PC
#' scores, etc.)
#' 
#' @param object Seurat object
#' @param features.plot Features to plot (gene expression, metrics, PC scores,
#' anything that can be retreived by fetch.data)
#' @param nCol Number of columns if multiple plots are displayed
#' @param ylab.max Maximum ylab value
#' @param do.ret FALSE by default. If TRUE, returns a list of ggplot objects.
#' @param do.sort Sort identity classes (on the x-axis) by the average
#' expression of the attribute being potted
#' @param size.x.use X axis font size
#' @param size.y.use Y axis font size
#' @param size.title.use Title font size
#' @param use.imputed Use imputed values for gene expression (default is FALSE)
#' @param adjust.use Adjust parameter for geom_violin
#' @param size.use Point size for geom_violin
#' @param cols.use Colors to use for plotting
#' @param group.by Group (color) cells in different ways (for example, orig.ident)
#' @import grid
#' @import gridExtra
#' @import ggplot2
#' @import reshape2
#' @return By default, no return, only graphical output. If do.return=TRUE,
#' returns a list of ggplot objects.
#' @export
setGeneric("vlnPlot", function(object,features.plot,nCol=NULL,ylab.max=12,do.ret=TRUE,do.sort=FALSE,
                               size.x.use=16,size.y.use=16,size.title.use=20, use.imputed=FALSE,adjust.use=1,size.use=1,cols.use=NULL,group.by="ident")  standardGeneric("vlnPlot"))
#' @export
setMethod("vlnPlot","seurat",
          function(object,features.plot,nCol=NULL,ylab.max=12,do.ret=FALSE,do.sort=FALSE,size.x.use=16,size.y.use=16,size.title.use=20,use.imputed=FALSE,adjust.use=1,
                   size.use=1,cols.use=NULL,group.by=NULL) {
            if (is.null(nCol)) {
              nCol=2
              if (length(features.plot)>6) nCol=3
              if (length(features.plot)>9) nCol=4
            }
            data.use=data.frame(t(fetch.data(object,features.plot,use.imputed=use.imputed)))
            #print(head(data.use))
            ident.use=object@ident
            if (!is.null(group.by)) ident.use=as.factor(fetch.data(object,group.by)[,1])
            pList=lapply(features.plot,function(x) plot.Vln(x,data.use[x,],ident.use,ylab.max,TRUE,do.sort,size.x.use,size.y.use,size.title.use,adjust.use,size.use,cols.use))
            
            if(do.ret) {
              return(pList)
            }
            else {
              multiplotList(pList,cols = nCol)
              rp()
            }
          }
)  

#' Add Metadata
#'
#' Adds additional data for single cells to the Seurat object. Can be any piece
#' of information associated with a cell (examples include read depth,
#' alignment rate, experimental batch, or subpopulation identity). The
#' advantage of adding it to the Seurat object is so that it can be
#' analyzed/visualized using fetch.data, vlnPlot, genePlot, subsetData, etc.
#'
#'
#' @param object Seurat object
#' @param metadata Data frame where the row names are cell names (note : these
#' must correspond exactly to the items in object@@cell.names), and the columns
#' are additional metadata items.
#' @return Seurat object where the additional metadata has been added as
#' columns in object@@data.info
#' @export
setGeneric("addMetaData", function(object,metadata)  standardGeneric("addMetaData"))
#' @export
setMethod("addMetaData","seurat",
          function(object,metadata) {
            cols.add=colnames(metadata)
            object@data.info[,cols.add]=metadata[rownames(object@data.info),]
            return(object)
          }
)

facet_wrap_labeller <- function(gg.plot,labels=NULL) {  
  # code from stackoverflow: http://stackoverflow.com/questions/19282897/how-to-add-expressions-to-labels-in-facet-wrap
  #works with R 3.0.1 and ggplot2 0.9.3.1
  
  g <- ggplotGrob(gg.plot)
  gg <- g$grobs      
  strips <- grep("strip_t", names(gg))
  
  for(ii in seq_along(labels))  {
    modgrob <- getGrob(gg[[strips[ii]]], "strip.text", 
                       grep=TRUE, global=TRUE)
    gg[[strips[ii]]]$children[[modgrob$name]] <- editGrob(modgrob,label=labels[ii])
  }
  
  g$grobs <- gg
  class(g) = c("arrange", "ggplot",class(g)) 
  g
}



#' JackStraw Plot
#'
#' Plots the results of the JackStraw analysis for PCA significance. For each
#' PC, plots a QQ-plot comparing the distribution of p-values for all genes
#' across each PC, compared with a uniform distribution. Also determines a
#' p-value for the overall significance of each PC (see Details).
#'
#' Significant PCs should show a p-value distribution (black curve) that is 
#' strongly skewed to the left compared to the null distribution (dashed line)
#' The p-value for each PC is based on a proportion test comparing the number
#' of genes with a p-value below a particular threshold (score.thresh), compared with the
#' proportion of genes expected under a uniform distribution of p-values. 
#'
#' @param object Seurat plot
#' @param plot.x.lim X-axis maximum on each QQ plot.
#' @param plot.y.lim Y-axis maximum on each QQ plot.
#' @param PCs Which PCs to examine
#' @param nCol Number of columns
#' @param score.thresh Threshold to use for the proportion test of PC
#' significance (see Details)
#' @return No value returned, just the PC plots.
#' @author Thanks to Omri Wurtzel for integrating with ggplot
#' @import gridExtra
#' @export
setGeneric("jackStrawPlot", function(object,PCs=1:5, nCol=3, score.thresh=1e-5,plot.x.lim=0.1,plot.y.lim=0.3)  standardGeneric("jackStrawPlot"))
#' @export
setMethod("jackStrawPlot","seurat",
          function(object,PCs=1:5, nCol=3, score.thresh=1e-5,plot.x.lim=0.1,plot.y.lim=0.3) {
            pAll=object@jackStraw.empP
            pAll <- pAll[,PCs, drop=F]
            pAll$Contig <- rownames(pAll)
            pAll.l <- melt(pAll, id.vars = "Contig")
            colnames(pAll.l) <- c("Contig", "PC", "Value")
            
            qq.df <- NULL
            score.df <- NULL
            for (i in PCs){
              q <- qqplot(pAll[, i],runif(1000),plot.it=FALSE)
              
              #pc.score=mean(q$y[which(q$x <=score.thresh)])
              pc.score=prop.test(c(length(which(pAll[, i] <= score.thresh)), floor(nrow(pAll) * score.thresh)), c(nrow(pAll), nrow(pAll)))$p.val
              if (length(which(pAll[, i] <= score.thresh))==0) pc.score=1
              
              if(is.null(score.df))
                score.df <- data.frame(PC=paste("PC",i, sep=""), Score=pc.score)
              else
                score.df <- rbind(score.df, data.frame(PC=paste("PC",i, sep=""), Score=pc.score))
              
              
              if (is.null(qq.df))
                qq.df <- data.frame(x=q$x, y=q$y, PC=paste("PC",i, sep=""))
              else
                qq.df <- rbind(qq.df, data.frame(x=q$x, y=q$y, PC=paste("PC",i, sep="")))
            }

            # create new dataframe column to wrap on that includes the PC number and score
            pAll.l$PC.Score <- rep(paste0(score.df$PC, " ", sprintf("%1.3g", score.df$Score)), each=length(unique(pAll.l$Contig)))
            pAll.l$PC.Score <- factor(pAll.l$PC.Score, levels = paste0(score.df$PC, " ", sprintf("%1.3g", score.df$Score)))
            gp <- ggplot(pAll.l, aes(sample=Value)) + stat_qq(distribution=qunif) + facet_wrap("PC.Score", ncol = nCol) + labs(x="Theoretical [runif(1000)]", y = "Empirical") +  xlim(0,plot.y.lim) + ylim(0,plot.x.lim) + coord_flip() + geom_abline(intercept=0, slope=1, linetype="dashed",na.rm=T) + theme_bw()
            return(gp)

          })

#' Scatter plot of single cell data
#'
#' Creates a scatter plot of two features (typically gene expression), across a
#' set of single cells. Cells are colored by their identity class.
#' 
#' @param object Seurat object
#' @param gene1 First feature to plot. Typically gene expression but can also
#' be metrics, PC scores, etc. - anything that can be retreived with fetch.data
#' @param gene2 Second feature to plot.
#' @param cell.ids Cells to include on the scatter plot.
#' @param col.use Colors to use for identity class plotting.
#' @param pch.use Pch argument for plotting
#' @param cex.use Cex argument for plotting
#' @param use.imputed Use imputed values for gene expression (Default is FALSE)
#' @param do.ident False by default. If TRUE,
#' @param do.spline Add a spline (currently hardwired to df=4, to be improved)
#' @param \dots Additional arguments to be passed to plot.
#' @return No return, only graphical output
#' @export
setGeneric("genePlot", function(object, gene1, gene2, cell.ids=NULL,col.use=NULL,
                                pch.use=16,cex.use=1.5,use.imputed=FALSE,do.ident=FALSE,do.spline=FALSE,...)  standardGeneric("genePlot"))
#' @export
setMethod("genePlot","seurat",
          function(object, gene1, gene2, cell.ids=NULL,col.use=NULL,
                   pch.use=16,cex.use=1.5,use.imputed=FALSE,do.ident=FALSE,do.spline=FALSE,...) {
            cell.ids=set.ifnull(cell.ids,object@cell.names)
            data.use=data.frame(t(fetch.data(object,c(gene1,gene2),cells.use = cell.ids,use.imputed=use.imputed)))
            corner(data.use)
            g1=as.numeric(data.use[gene1,cell.ids])
            g2=as.numeric(data.use[gene2,cell.ids])
            ident.use=as.factor(object@ident[cell.ids])
            if (length(col.use)>1) {
              col.use=col.use[as.numeric(ident.use)]
            }
            else {
              col.use=set.ifnull(col.use,as.numeric(ident.use))
            }
            gene.cor=round(cor(g1,g2),2)
            plot(g1,g2,xlab=gene1,ylab=gene2,col=col.use,cex=cex.use,main=gene.cor,pch=pch.use,...)
            if (do.spline) {
              spline.fit=smooth.spline(g1,g2,df = 4)
              lines(spline.fit$x,spline.fit$y,lwd=3)
            }
            if (do.ident) {
              return(identify(g1,g2,labels = cell.ids))
            }
          }
)

#' @export
setGeneric("removePC", function(object, pcs.remove,...)  standardGeneric("removePC"))
#' @export
setMethod("removePC","seurat",
          function(object, pcs.remove,...) {
            data.old=object@data
            pcs.use=anotinb(1:ncol(object@pca.obj[[1]]$rotation),pcs.remove)
            data.1=as.matrix(object@pca.obj[[1]]$x[,pcs.use])%*%t(as.matrix(object@pca.obj[[1]]$rotation[,pcs.use]))
            data.2=sweep(data.1,2,colMeans(object@scale.data),"+")
            data.3=sweep(data.2,MARGIN = 1,apply(object@data,1,sd),"*")         
            data.3=sweep(data.3,MARGIN = 1,apply(object@data,1,mean),"+")         
            object@scale.data=(data.2); object@data=data.frame(data.3)
            object@data[data.old==0]=0; object@data[object@data<0]=0
            return(object)
          }
)

setGeneric("geneScorePlot", function(object, gene1, score.name, cell.ids=NULL,col.use=NULL,nrpoints.use=Inf,pch.use=16,cex.use=2,...)  standardGeneric("geneScorePlot"))
setMethod("geneScorePlot","seurat",
          function(object, gene1, score.name, cell.ids=NULL,col.use=NULL,nrpoints.use=Inf,pch.use=16,cex.use=2,...) {
            cell.ids=set.ifnull(cell.ids,colnames(object@data))
            g1=as.numeric(object@data[gene1,cell.ids])
            my.score=retreiveScore(object,score.name)
            s1=as.numeric(my.score[cell.ids])
            col.use=set.ifnull(as.numeric(as.factor(object@ident[cell.ids])))
            gene.cor=round(cor(g1,s1),2)
            smoothScatter(g1,s1,xlab=gene1,ylab=score.name,col=col.use,nrpoints=nrpoints.use,cex=cex.use,main=gene.cor,pch=pch.use)
          }
)

#' Cell-cell scatter plot
#'
#' Creates a plot of scatter plot of genes across two single cells
#' 
#' @param object Seurat object
#' @param cell1 Cell 1 name (can also be a number, representing the position in
#' object@@cell.names)
#' @param cell2 Cell 2 name (can also be a number, representing the position in
#' object@@cell.names)
#' @param gene.ids Genes to plot (default, all genes)
#' @param col.use Colors to use for the points
#' @param nrpoints.use Parameter for smoothScatter
#' @param pch.use Point symbol to use
#' @param cex.use Point size
#' @param do.ident FALSE by default. If TRUE, allows you to click on individual
#' points to reveal gene names (hit ESC to stop)
#' @param \dots Additional arguments to pass to smoothScatter
#' @return No return value (plots a scatter plot)
#' @export
setGeneric("cellPlot", function(object, cell1, cell2, gene.ids=NULL,col.use="black",nrpoints.use=Inf,pch.use=16,cex.use=0.5,do.ident=FALSE,...)  standardGeneric("cellPlot"))
#' @export
setMethod("cellPlot","seurat",
          function(object, cell1, cell2, gene.ids=NULL,col.use="black",nrpoints.use=Inf,pch.use=16,cex.use=0.5,do.ident=FALSE,...) {
            gene.ids=set.ifnull(gene.ids,rownames(object@data))
            c1=as.numeric(object@data[gene.ids,cell1])
            c2=as.numeric(object@data[gene.ids,cell2])
            gene.cor=round(cor(c1,c2),2)
            smoothScatter(c1,c2,xlab=cell1,ylab=cell2,col=col.use,nrpoints=nrpoints.use,pch=pch.use,cex=cex.use,main=gene.cor)
            if (do.ident) {
              identify(c1,c2,labels = gene.ids)
            }
          }
)

#' @export
#' @import jackstraw
setGeneric("jackStraw.permutation.test", function(object,genes.use=NULL,num.iter=100, thresh.use=0.05,do.print=TRUE,k.seed=1)  standardGeneric("jackStraw.permutation.test"))
#' @export
setMethod("jackStraw.permutation.test","seurat",
          function(object,genes.use=NULL,num.iter=100, thresh.use=0.05,do.print=TRUE,k.seed=1) {
            genes.use=set.ifnull(genes.use,rownames(object@pca.x))
            genes.use=ainb(genes.use,rownames(object@scale.data))
            data.use=t(as.matrix(object@scale.data[genes.use,]))
            if (do.print) print(paste("Running ", num.iter, " iterations",sep=""))
            pa.object=permutationPA(data.use,B = num.iter,threshold = thresh.use,verbose = do.print,seed = k.seed)
            if (do.print) cat("\n\n")
            if (do.print) print(paste("JackStraw returns ", pa.object$r, " significant components", sep=""))
            return(pa.object)
          }
)


#multicore version of jackstraw
#DOES NOT WORK WITH WINDOWS
#' @export
setGeneric("jackStrawMC", function(object,num.pc=30,num.replicate=100,prop.freq=0.01,do.print=FALSE, num.cores=8)  standardGeneric("jackStrawMC"))
#' @export
setMethod("jackStrawMC","seurat",
          function(object,num.pc=30,num.replicate=100,prop.freq=0.01,do.print=FALSE, num.cores=8) {
            pc.genes=rownames(object@pca.x)
            if (length(pc.genes)<200) prop.freq=max(prop.freq,0.015)
            md.x=as.matrix(object@pca.x)
            md.rot=as.matrix(object@pca.rot)
            if (!(do.print)) fake.pcVals.raw=mclapply(1:num.replicate, function(x)jackRandom(scaled.data=object@scale.data[pc.genes,],prop=prop.freq,r1.use = 1,r2.use = num.pc,seed.use=x), mc.cores = num.cores)
            if ((do.print)) fake.pcVals.raw=mclapply(1:num.replicate,function(x){ print(x); jackRandom(scaled.data=object@scale.data[pc.genes,],prop=prop.freq,r1.use = 1,r2.use = num.pc,seed.use=x)}, mc.cores=num.cores)
            
            fake.pcVals=simplify2array(mclapply(1:num.pc,function(x)as.numeric(unlist(lapply(1:num.replicate,function(y)fake.pcVals.raw[[y]][,x]))), mc.cores=num.cores))
            object@jackStraw.fakePC = data.frame(fake.pcVals)
            object@jackStraw.empP=data.frame(simplify2array(mclapply(1:num.pc,function(x)unlist(lapply(abs(md.x[,x]),empP,abs(fake.pcVals[,x]))), mc.cores=num.cores)))
            colnames(object@jackStraw.empP)=paste("PC",1:ncol(object@jackStraw.empP),sep="")
            return(object)
          }
)

#' Determine statistical significance of PCA scores.
#'
#' Randomly permutes a subset of data, and calculates projected PCA scores for
#' these 'random' genes. Then compares the PCA scores for the 'random' genes
#' with the observed PCA scores to determine statistical signifance. End result
#' is a p-value for each gene's association with each principal component.
#'
#'
#' @param object Seurat object
#' @param num.pc Number of PCs to compute significance for
#' @param num.replicate Number of replicate samplings to perform
#' @param prop.freq Proportion of the data to randomly permute for each
#' replicate
#' @param do.print Print the number of replicates that have been processed.
#' @return Returns a Seurat object where object@@jackStraw.empP represents
#' p-values for each gene in the PCA analysis. If project.pca is subsequently
#' run, object@@jackStraw.empP.full then represents p-values for all genes.
#' @references Inspired by Chung et al, Bioinformatics (2014)
#' @export
setGeneric("jackStraw", function(object,num.pc=30,num.replicate=100,prop.freq=0.01,do.print=FALSE)  standardGeneric("jackStraw"))
#' @export
setMethod("jackStraw","seurat",
          function(object,num.pc=30,num.replicate=100,prop.freq=0.01,do.print=FALSE) {
            
            # error checking for number of PCs
            if (num.pc >= ncol(object@pca.rot)){
              num.pc <- ncol(object@pca.rot)
              warning("Number of PCs specified is greater than PCs available. Setting num.pc to ", num.pc, " and continuing.")
            }
            if (num.pc > length(object@cell.names)){
              num.pc <- length(object@cell.names)
              warning("Number of PCs specified is greater than number of cells. Setting num.pc to ", num.pc, " and continuing.")
            }
            
            pc.genes=rownames(object@pca.x)
            
            if (length(pc.genes) < 3){
              stop("Too few variable genes")
            }
            if (length(pc.genes) * prop.freq < 3){
              warning("Number of variable genes given ", prop.freq, " as the prop.freq is low. Consider including more variable genes and/or increasing prop.freq. ", 
                      "Continuing with 3 genes in every random sampling.")
            }
            
            md.x=as.matrix(object@pca.x)
            md.rot=as.matrix(object@pca.rot)
            if (!(do.print)) fake.pcVals.raw=sapply(1:num.replicate,function(x)jackRandom(scaled.data=object@scale.data[pc.genes,],prop=prop.freq,r1.use = 1,r2.use = num.pc,seed.use=x),simplify = FALSE)
            if ((do.print)) fake.pcVals.raw=sapply(1:num.replicate,function(x){ print(x); jackRandom(scaled.data=object@scale.data[pc.genes,],prop=prop.freq,r1.use = 1,r2.use = num.pc,seed.use=x)},simplify = FALSE)
            
            fake.pcVals=sapply(1:num.pc,function(x)as.numeric(unlist(lapply(1:num.replicate,function(y)fake.pcVals.raw[[y]][,x]))))
            object@jackStraw.fakePC = data.frame(fake.pcVals)
            object@jackStraw.empP=data.frame(sapply(1:num.pc,function(x)unlist(lapply(abs(md.x[,x]),empP,abs(fake.pcVals[,x])))))
            colnames(object@jackStraw.empP)=paste("PC",1:ncol(object@jackStraw.empP),sep="")
            return(object)
          }
)

#' @export
jackRandom=function(scaled.data,prop.use=0.01,r1.use=1,r2.use=5, seed.use=1) {
  set.seed(seed.use); rand.genes=sample(rownames(scaled.data),nrow(scaled.data)*prop.use)
  
  # make sure that rand.genes is at least 3
  if (length(rand.genes) < 3){
    rand.genes <- sample(rownames(scaled.data), 3)
  }

  data.mod=scaled.data
  data.mod[rand.genes,]=shuffleMatRow(scaled.data[rand.genes,])
  fake.pca=prcomp(data.mod)
  fake.x=fake.pca$x
  fake.rot=fake.pca$rotation
  return(fake.x[rand.genes,r1.use:r2.use])
}

setGeneric("jackStrawFull", function(object,num.pc=5,num.replicate=100,prop.freq=0.01)  standardGeneric("jackStrawFull"))
setMethod("jackStrawFull","seurat",
          function(object,num.pc=5,num.replicate=100,prop.freq=0.01) {
            pc.genes=rownames(object@pca.x)
            if (length(pc.genes)<200) prop.freq=max(prop.freq,0.015)
            md.x=as.matrix(object@pca.x)
            md.rot=as.matrix(object@pca.rot)
            real.fval=sapply(1:num.pc,function(x)unlist(lapply(pc.genes,jackF,r1=x,r2=x,md.x,md.rot)))
            rownames(real.fval)=pc.genes
            object@real.fval=data.frame(real.fval)
            
            fake.fval=sapply(1:num.pc,function(x)unlist(replicate(num.replicate,
                                                                  jackStrawF(prop=prop.freq,data=object@scale.data[pc.genes,],myR1 = x,myR2 = x),simplify=FALSE)))
            rownames(fake.fval)=1:nrow(fake.fval)
            object@fake.fval=data.frame(fake.fval)
            
            object@emp.pval=data.frame(sapply(1:num.pc,function(x)unlist(lapply(object@real.fval[,x],empP,object@fake.fval[,x]))))
            
            rownames(object@emp.pval)=pc.genes
            colnames(object@emp.pval)=paste("PC",1:ncol(object@emp.pval),sep="")
            return(object)
          }
)


#' Identify variable genes
#'
#' Identifies genes that are outliers on a 'mean variability plot'. First, uses
#' a function to calculate average expression (fxn.x) and dispersion (fxn.y)
#' for each gene. Next, divides genes into num.bin (deafult 20) bins based on
#' their average expression, and calculates z-scores for dispersion within each
#' bin. The purpose of this is to identify variable genes while controlling for
#' the strong relationship between variability and average expression.
#'
#' Exact parameter settings may vary empirically from dataset to dataset, and
#' based on visual inspection of the plot. 
#' Setting the y.cutoff parameter to 2 identifies genes that are more than two standard
#' deviations away from the average dispersion within a bin. The default X-axis function
#' is the mean expression level, and for Y-axis it is the log(Variance/mean). All mean/variance
#' calculations are not performed in log-space, but the results are reported in log-space - 
#' see relevant functions for exact details.
#'
#' @param object Seurat object
#' @param fxn.x Function to compute x-axis value (average expression). Default
#' is to take the mean of the detected (i.e. non-zero) values
#' @param fxn.y Function to compute y-axis value (dispersion). Default is to
#' take the standard deviation of all values/
#' @param do.plot Plot the average/dispersion relationship
#' @param set.var.genes Set object@@var.genes to the identified variable genes
#' (default is TRUE)
#' @param do.text Add text names of variable genes to plot (default is TRUE)
#' @param x.low.cutoff Bottom cutoff on x-axis for identifying variable genes
#' @param x.high.cutoff Top cutoff on x-axis for identifying variable genes
#' @param y.cutoff Bottom cutoff on y-axis for identifying variable genes
#' @param y.high.cutoff Top cutoff on y-axis for identifying variable genes
#' @param cex.use Point size
#' @param cex.text.use Text size
#' @param do.spike FALSE by default. If TRUE, color all genes starting with ^ERCC a different color
#' @param pch.use Pch value for points
#' @param col.use Color to use
#' @param spike.col.use if do.spike, color for spike-in genes
#' @param plot.both Plot both the scaled and non-scaled graphs.
#' @param do.contour Draw contour lines calculated based on all genes
#' @param contour.lwd Contour line width
#' @param contour.col Contour line color
#' @param contour.lty Contour line type
#' @param num.bin Total number of bins to use in the scaled analysis (default
#' is 20)
#' @importFrom MASS kde2d
#' @return Returns a Seurat object, placing variable genes in object@@var.genes.
#' The result of all analysis is stored in object@@mean.var
#' @export
setGeneric("mean.var.plot", function(object, fxn.x=expMean, fxn.y=logVarDivMean,do.plot=TRUE,set.var.genes=TRUE,do.text=TRUE,
                                     x.low.cutoff=4,x.high.cutoff=8,y.cutoff=2,y.high.cutoff=12,cex.use=0.5,cex.text.use=0.5,do.spike=FALSE, 
                                     pch.use=16, col.use="black", spike.col.use="red",plot.both=FALSE,do.contour=TRUE,
                                     contour.lwd=3, contour.col="white", contour.lty=2,num.bin=20) standardGeneric("mean.var.plot"))
#' @export
setMethod("mean.var.plot", signature = "seurat",
          function(object, fxn.x=humpMean, fxn.y=sd,do.plot=TRUE,set.var.genes=TRUE,do.text=TRUE,
                   x.low.cutoff=4,x.high.cutoff=8,y.cutoff=1,y.high.cutoff=12,cex.use=0.5,cex.text.use=0.5,do.spike=FALSE, 
                   pch.use=16, col.use="black", spike.col.use="red",plot.both=FALSE,do.contour=TRUE,
                   contour.lwd=3, contour.col="white", contour.lty=2,num.bin=20) {
            data=object@data
            data.x=apply(data,1,fxn.x); data.y=apply(data,1,fxn.y); data.x[is.na(data.x)]=0
            data.norm.y=meanNormFunction(data,fxn.x,fxn.y,num.bin)
            data.norm.y[is.na(data.norm.y)]=0
            names(data.norm.y)=names(data.x)
            pass.cutoff=names(data.x)[which(((data.x>x.low.cutoff) & (data.x<x.high.cutoff)) & (data.norm.y>y.cutoff) & (data.norm.y < y.high.cutoff))]
            mv.df=data.frame(data.x,data.y,data.norm.y)
            rownames(mv.df)=rownames(data)
            object@mean.var=mv.df
            if (do.spike) spike.genes=rownames(subr(data,"^ERCC"))
            if (do.plot) {
              if (plot.both) {
                par(mfrow=c(1,2)) 
                smoothScatter(data.x,data.y,pch=pch.use,cex=cex.use,col=col.use,xlab="Average expression",ylab="Dispersion",nrpoints=Inf)
                
                if (do.contour) {
                  data.kde=kde2d(data.x,data.y)
                  contour(data.kde,add=TRUE,lwd=contour.lwd,col=contour.col,lty=contour.lty)
                }
                if (do.spike) points(data.x[spike.genes],data.y[spike.genes],pch=16,cex=cex.use,col=spike.col.use)
                if(do.text) text(data.x[pass.cutoff],data.y[pass.cutoff],pass.cutoff,cex=cex.text.use)
              }
              smoothScatter(data.x,data.norm.y,pch=pch.use,cex=cex.use,col=col.use,xlab="Average expression",ylab="Dispersion",nrpoints=Inf)
              if (do.contour) {
                data.kde=kde2d(data.x,data.norm.y)
                contour(data.kde,add=TRUE,lwd=contour.lwd,col=contour.col,lty=contour.lty)
              }
              if (do.spike) points(data.x[spike.genes],data.norm.y[spike.genes],pch=16,cex=cex.use,col=spike.col.use,nrpoints=Inf)
              if(do.text) text(data.x[pass.cutoff],data.norm.y[pass.cutoff],pass.cutoff,cex=cex.text.use)
            }
            if (set.var.genes) { 
              object@var.genes=pass.cutoff
              return(object)
              if (!set.var.genes) return(pass.cutoff)
            }
          }
          
)

