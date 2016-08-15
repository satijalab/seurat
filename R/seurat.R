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
#'    \item{\code{raw.data}:}{\code{"ANY"}, The raw project data }
#'    \item{\code{data}:}{\code{"ANY"}, The expression matrix (log-scale) }
#'    \item{\code{scale.data}:}{\code{"ANY"}, The scaled (after z-scoring
#'    each gene) expression matrix. Used for PCA, ICA, and heatmap plotting}
#'    \item{\code{var.genes}:}{\code{"vector"},  Variable genes across single cells }
#'    \item{\code{is.expr}:}{\code{"numeric"}, Expression threshold to determine if a gene is expressed }
#'    \item{\code{ident}:}{\code{"vector"},  The 'identity class' for each single cell }
#'    \item{\code{data.info}:}{\code{"data.frame"}, Contains information about each cell, starting with # of genes detected (nGene)
#'    the original identity class (orig.ident), user-provided information (through AddMetaData), etc.  }
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
#'    \item{\code{snn.sparse}:}{\code{"dgCMatrix"}, Sparse matrix object representation of the SNN graph }
#'    \item{\code{snn.dense}:}{\code{"matrix"}, Dense matrix object representation of the SNN graph }
#'    \item{\code{snn.k}:}{\code{"numeric"}, k used in the construction of the SNN graph }
#'
#'}
#' @name seurat
#' @rdname seurat
#' @aliases seurat-class
#' @exportClass seurat

seurat <- setClass("seurat", slots =
                     c(raw.data = "ANY", data="ANY",scale.data="ANY",var.genes="vector",is.expr="numeric",
                       ident="vector",pca.x="data.frame",pca.rot="data.frame",
                       emp.pval="data.frame",kmeans.obj="list",pca.obj="list",
                       gene.scores="data.frame", drop.coefs="data.frame",
                       wt.matrix="data.frame", drop.wt.matrix="data.frame",trusted.genes="vector",drop.expr="numeric",data.info="data.frame",
                       project.name="character", kmeans.gene="list", kmeans.cell="list",jackStraw.empP="data.frame",
                       jackStraw.fakePC = "data.frame",jackStraw.empP.full="data.frame",pca.x.full="data.frame", kmeans.col="list",mean.var="data.frame", imputed="data.frame",mix.probs="data.frame",
                       mix.param="data.frame",final.prob="data.frame",insitu.matrix="data.frame",
                       tsne.rot="data.frame", ica.rot="data.frame", ica.x="data.frame", ica.obj="list",cell.names="vector",cluster.tree="list",
                       snn.sparse="dgCMatrix", snn.dense="matrix", snn.k="numeric"))


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
#' @import pbapply
#' @importFrom Matrix colSums rowSums
#' @export
setGeneric("Setup", function(object, project, min.cells=3, min.genes=1000, is.expr=0, do.scale=TRUE, do.center=TRUE,names.field=1,names.delim="_",meta.data=NULL,save.raw=TRUE,...) standardGeneric("Setup"))
#' @export
setMethod("Setup","seurat",
          function(object, project, min.cells=3, min.genes=1000, is.expr=0, do.scale=TRUE, do.center=TRUE,names.field=1,names.delim="_",meta.data=NULL,save.raw=TRUE,...) {
            object@is.expr <- is.expr
            num.genes <- colSums(object@raw.data > is.expr)
            cells.use <- names(num.genes[which(num.genes > min.genes)])
            object@data <- object@raw.data[, cells.use]
            
            #to save memory downstream, especially for large object
            if (!(save.raw)) object@raw.data <- matrix();
            genes.use <- rownames(object@data)
            if (min.cells > 0) {
              num.cells <- rowSums(object@data > is.expr)
              genes.use <- names(num.cells[which(num.cells >= min.cells)])
              object@data <- object@data[genes.use, ]
            }
            
            object@ident <- factor(unlist(lapply(colnames(object@data), extract_field, names.field, names.delim)))
            names(object@ident) <- colnames(object@data)
            object@cell.names <- names(object@ident)
            
            # if there are more than 100 idents, set all ident to project name
            if(length(unique(object@ident)) > 100 || length(unique(object@ident)) == 0) {
              object <- SetIdent(object, ident.use = project)
            }
            object@scale.data <- matrix()
            if(do.scale | do.center) {
              object=ScaleData(object,do.scale = do.scale,do.center = do.center)
            }
            
            data.ngene <- num.genes[cells.use]
            object@gene.scores <- data.frame(data.ngene)
            colnames(object@gene.scores)[1] <- "nGene"
            
            object@data.info <- data.frame(data.ngene)
            colnames(object@data.info)[1] <- "nGene"
            
            if (!is.null(meta.data)) {
              object <- AddMetaData(object ,metadata = meta.data)
            }
            object@mix.probs <- data.frame(data.ngene)
            colnames(object@mix.probs)[1] <- "nGene"
            rownames(object@gene.scores) <- colnames(object@data)
            
            object@data.info[names(object@ident),"orig.ident"] <- object@ident
            
            object@project.name <- project
            #if(calc.noise) {
            #  object=CalcNoiseModels(object,...)
            #  object=GetWeightMatrix(object)
            #}
            return(object)
          }
)



#' Scale and center the data
#'
#'
#' @param object Seurat object
#' @param genes.use Vector of gene names to scale/center. Default is all genes in object@@data.
#' @param data.use Can optionally pass a matrix of data to scale 
#' @param do.scale Whether to scale the data. 
#' @param do.center Whether to center the data.
#' @param scale.max Max value to accept for scaled data. The default is 10. Setting this can help 
#' reduce the effects of genes that are only expressed in a very small number of cells. 
#' @return Returns a seurat object with object@@scale.data updated with scaled and/or centered data.
setGeneric("ScaleData", function(object, genes.use=NULL, data.use=NULL, do.scale=TRUE, do.center=TRUE, scale.max=10) standardGeneric("ScaleData"))
#' @export
setMethod("ScaleData", "seurat",
          function(object, genes.use=NULL, data.use=NULL, do.scale=TRUE, do.center=TRUE, scale.max=10) {
            genes.use <- set.ifnull(genes.use,rownames(object@data))
            data.use <- set.ifnull(data.use,object@data[genes.use, ])
            object@scale.data <- matrix(NA, nrow = length(genes.use), ncol = ncol(object@data))
            #rownames(object@scale.data) <- genes.use 
            #colnames(object@scale.data) <- colnames(object@data)
            dimnames(object@scale.data) <- dimnames(data.use)
            if(do.scale | do.center) {
              bin.size <- 1000
              max.bin <- floor(length(genes.use)/bin.size) + 1
              pb <- txtProgressBar(min = 0, max = max.bin, style = 3)
              for(i in 1:max.bin) {
                my.inds <- ((bin.size * (i - 1)):(bin.size * i - 1))+1
                my.inds <- my.inds[my.inds <= length(genes.use)]
                #print(my.inds)
                new.data <- t(scale(t(as.matrix(data.use[genes.use[my.inds], ])), center = do.center, scale = do.scale))
                new.data[new.data>scale.max] <- scale.max
                object@scale.data[genes.use[my.inds], ] <- new.data
                setTxtProgressBar(pb, i)  
              }
              close(pb)
            }
            return(object)
          }
)

#' Normalize raw data
#'
#' Normalize count data per cell and transform to log scale
#'
#'
#' @param data Matrix with the raw count data
#' @return Returns a matrix with the normalize and log transformed data
#' @importFrom Matrix colSums 
#' @export
setGeneric("LogNormalize", function(data, scale.factor = 1e4) standardGeneric("LogNormalize"))
#' @export
setMethod("LogNormalize", "ANY",
          function(data, scale.factor = 1e4) {
            if(is.matrix(data) || class(data) == "data.frame") {
              return(log(sweep(data, 2, colSums(data), FUN = "/") * scale.factor + 1))
            }
            else {
              cells.use <- colnames(data)
              bin.size <- 1000
              max.bin <- floor(length(cells.use)/bin.size) + 1
              pb <- txtProgressBar(min = 0, max = max.bin, style = 3)
              for(i in 1:max.bin) {
                my.inds <- ((bin.size * (i - 1)):(bin.size * i - 1))+1
                my.inds <- my.inds[my.inds <= length(cells.use)]
                cells.iter <- cells.use[my.inds]
                data.iter <- data[, my.inds]
                data.new <- sweep(data.iter, 2, colSums(data.iter), FUN = "/")
                data.new <- log1p(data.new * scale.factor)
                if (i == 1) all.norm <- data.new
                if (i > 1) all.norm <- cbind(all.norm, data.new)
                setTxtProgressBar(pb, i)  
              }
              close(pb)
              return(all.norm)
            }
          }
)


calc.drop.prob=function(x,a,b) {
  return(exp(a+b*x)/(1+exp(a+b*x)))
}


setGeneric("FindAllMarkersNode", function(object, thresh.test=1,test.use="bimod",return.thresh=1e-2,do.print=FALSE,...) standardGeneric("FindAllMarkersNode"))
setMethod("FindAllMarkersNode","seurat",
          function(object, thresh.test=1,test.use="bimod",return.thresh=1e-2,do.print=FALSE,...) {
            ident.use=object@ident
            tree.use=object@cluster.tree[[1]]

            if ((test.use=="roc") && (return.thresh==1e-2)) return.thresh=0.8
            genes.de=list()
            for(i in ((tree.use$Nnode+2):max(tree.use$edge))) {
              genes.de[[i]]=FindMarkersNode(object,i,genes.use=rownames(object@data),thresh.use = thresh.test,test.use = test.use,...)
              if (do.print) print(paste("Calculating node", i))
            }
            gde.all=data.frame()
            for(i in ((tree.use$Nnode+2):max(tree.use$edge))) {
              gde=genes.de[[i]]
              if (is.null(gde)) next;
              if (nrow(gde)>0) {
                if (test.use=="roc") gde=subset(gde,(myAUC>return.thresh|myAUC<(1-return.thresh)))
                if ((test.use=="bimod")||(test.use=="t")) {
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
#' @param SNN.use If SNN is passed, build tree based on SNN graph connectivity between clusters
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
setGeneric("BuildClusterTree", function(object, genes.use=NULL,pcs.use=NULL, SNN.use=NULL, do.plot=TRUE,do.reorder=FALSE,reorder.numeric=FALSE) standardGeneric("BuildClusterTree"))
#' @export
setMethod("BuildClusterTree","seurat",
          function(object,genes.use=NULL,pcs.use=NULL, SNN.use=NULL, do.plot=TRUE,do.reorder=FALSE,reorder.numeric=FALSE) {
            genes.use=set.ifnull(genes.use,object@var.genes)
            ident.names=as.character(unique(object@ident))
            if (!is.null(genes.use)) {
              genes.use=ainb(genes.use,rownames(object@data))
              data.avg=AverageExpression(object,genes.use = genes.use)
              data.dist=dist(t(data.avg[genes.use,]))
            }
            if (!is.null(pcs.use)) {
              data.pca=AveragePCA(object)
              data.use=data.pca[pcs.use,]
              data.eigenval=(object@pca.obj[[1]]$sdev)^2
              data.weights=(data.eigenval/sum(data.eigenval))[pcs.use]; data.weights=data.weights/sum(data.weights)
              data.dist=custom.dist(data.pca[pcs.use,],weighted.euclidean,data.weights)
            }
            if(!is.null(SNN.use)){
              num_clusters = length(ident.names)
              data.dist = matrix(0, nrow=num_clusters, ncol= num_clusters)
              for (i in 1:num_clusters-1){
                for (j in (i+1):num_clusters){
                  subSNN = SNN.use[match(WhichCells(object, i), colnames(SNN.use)), match(WhichCells(object, j), rownames(SNN.use))]
                  d = mean(subSNN)
                  if(is.na(d)) data.dist[i,j] = 0
                  else data.dist[i,j] = d
                }
              }
              diag(data.dist)=1
              data.dist=dist(data.dist)
            }
            data.tree=as.phylo(hclust(data.dist))
            object@cluster.tree[[1]]=data.tree

            if (do.reorder) {
              old.ident.order=sort(unique(object@ident))
              data.tree=object@cluster.tree[[1]]
              all.desc=getDescendants(data.tree,data.tree$Nnode+2); all.desc=old.ident.order[all.desc[all.desc<=(data.tree$Nnode+1)]]
              object@ident=factor(object@ident,levels = all.desc,ordered = TRUE)
              if (reorder.numeric) {
                object=SetIdent(object,object@cell.names,as.integer(object@ident))
                object@data.info[object@cell.names,"tree.ident"]=as.integer(object@ident)
              }
              object=BuildClusterTree(object,genes.use,pcs.use,do.plot=FALSE,do.reorder=FALSE)
            }
            if (do.plot) PlotClusterTree(object)
            return(object)
          }
)

#' Plot phylogenetic tree
#'
#' Plots previously computed phylogenetic tree (from BuildClusterTree)
#'
#' @param object Seurat object
#' @return Plots dendogram (must be precomputed using BuildClusterTree), returns no value
#' @importFrom ape plot.phylo
#' @importFrom ape nodelabels
#' @export
setGeneric("PlotClusterTree", function(object) standardGeneric("PlotClusterTree"))
#' @export
setMethod("PlotClusterTree","seurat",
          function(object) {
            if (length(object@cluster.tree)==0) stop("Phylogenetic tree does not exist, build using BuildClusterTree")
            data.tree=object@cluster.tree[[1]]
            plot.phylo(data.tree,direction="downwards")
            nodelabels()
          }
)


#' Visualize expression/dropout curve
#'
#' Plot the probability of detection vs average expression of a gene.
#'
#' Assumes that this 'noise' model has been precomputed with CalcNoiseModels
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
setGeneric("PlotNoiseModel", function(object, cell.ids=c(1,2), col.use="black",lwd.use=2,do.new=TRUE,x.lim=10,...) standardGeneric("PlotNoiseModel"))
#' @export
setMethod("PlotNoiseModel","seurat",
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

#' Regress out technical effects and cell cycle
#'
#' Remove unwanted effects from scale.data
#'
#'
#' @param object Seurat object
#' @param latent.vars effects to regress out
#' @param genes.regress gene to run regression for (default is all genes)
#' @param do.scale Z-normalize the residual values (default is TRUE)
#' @return Returns Seurat object with the scale.data (object@scale.data) genes returning the residuals from the regression model
#' @import pbapply
#' @import Matrix
#' @export
setGeneric("RegressOut", function(object,latent.vars,genes.regress=NULL,do.scale=TRUE,...) standardGeneric("RegressOut"))
#' @export
setMethod("RegressOut", "seurat",
          function(object,latent.vars,genes.regress=NULL,do.scale=TRUE,...) {
            require(Matrix)
            genes.regress=set.ifnull(genes.regress,rownames(object@data))
            genes.regress=ainb(genes.regress,rownames(object@data))
            latent.data=FetchData(object,latent.vars)
            
            #NEEDS TO BE FIXED!!!!
            exp.data=t(as.matrix(object@data[genes.regress,]))
            regression.mat=cbind(latent.data,exp.data)
            new.data=t(pbsapply(genes.regress, function(x) {
              regression.mat.2=latent.data
              regression.mat.2[,"GENE"] = regression.mat[,x];
              fmla=as.formula(paste("GENE ", " ~ ", paste(latent.vars,collapse="+"),sep=""));
              return(lm(fmla,data = regression.mat.2)$residuals)
            }))
            if (do.scale==TRUE) {
              object=ScaleData(object,genes.use = rownames(new.data),data.use = new.data,...)
            }
            object@scale.data[is.na(object@scale.data)]=0
            return(object)
          }
)

#' @export
setGeneric("AddSamples", function(object,new.data, min.genes=2500 ,names.field=1,names.delim="_") standardGeneric("AddSamples"))
#' @export
setMethod("AddSamples","seurat",
          function(object,new.data, min.genes=2500 ,names.field=1,names.delim="_") {
            geneMeans.old = apply(object@data, 1, mean)
            geneSd.old = apply(object@data, 1, sd)
            levels.old=levels(object@ident)
            num.genes=findNGene(new.data,object@is.expr)
            cells.use=names(num.genes[which(num.genes>min.genes)]); new.data=new.data[,cells.use]
            genes.old=rownames(object@data); genes.new=rownames(new.data)
            genes.same.1=which(genes.old%in%genes.new); genes.same.2=which(genes.new%in%genes.old)
            genes.diff.1=which(!(genes.old%in%genes.new)); genes.diff.2=which(!(genes.new%in%genes.old))
            new.data.1=new.data
            if (length(genes.diff.1)!=0) {
              data.bind.new.1=matrix(rep(0,length(genes.diff.1)*ncol(new.data)),nrow = length(genes.diff.1)); colnames(data.bind.new.1)=colnames(new.data); new.data.1=rbind(new.data,data.bind.new.1)
              rownames(new.data.1)[(nrow(new.data)+1):nrow(new.data.1)]=genes.old[genes.diff.1]; new.data.1=new.data.1[sort(rownames(new.data.1)),]
            }
            data.bind.new.2=matrix(rep(0,length(genes.diff.2)*ncol(object@data)),nrow = length(genes.diff.2)); colnames(data.bind.new.2)=colnames(object@data); new.data.2=rbind(object@data,data.bind.new.2)
            rownames(new.data.2)[(nrow(object@data)+1):nrow(new.data.2)]=genes.new[genes.diff.2]; new.data.2=new.data.2[sort(rownames(new.data.2)),]

            #genes.new=setdiff(rownames(new.data),rownames(object@data))
            #new.data=new.data[genes.use,cells.use]; new.data[is.na(new.data)]=0; rownames(new.data)=genes.use
            object@data=data.frame(cbind(new.data.2,new.data.1))

            new.means=apply(object@data[genes.new[genes.diff.2],],1,mean)
            new.sd=apply(object@data[genes.new[genes.diff.2],],1,sd)

            new.means=c(geneMeans.old,new.means); new.means=new.means[rownames(object@data)]
            new.sd=c(geneSd.old,new.sd); new.sd=new.sd[rownames(object@data)]

            data.new.scale = t(scale(t(object@data),center=new.means,scale=new.sd))
            #data.new.scale=data.new.scale[rownames(object@scale.data),]; data.new.scale[is.na(data.new.scale)]=-10
            #genes.diff=anotinb(rownames(object@scale.data),rownames(data.new.scale))
            #print(genes.diff)
            #object@scale.data = cbind(object@scale.data, data.new.scale)

            object@scale.data=data.new.scale
            new.ident=(unlist(lapply(colnames(new.data),extract_field,names.field,names.delim)))
            names(new.ident)=colnames(new.data)
            object@ident=factor(c(as.character(object@ident),as.character(new.ident)))
            names(object@ident)=colnames(object@data)
            object@cell.names=names(object@ident)

            info.rows=ncol(new.data); info.cols=ncol(object@data.info)
            new.names=colnames(new.data)
            new.data.info=data.frame(matrix(rep(0,info.rows*info.cols),nrow = info.rows),row.names = new.names)
            colnames(new.data.info)=colnames(object@data.info)
            new.data.info[new.names,"nGene"]=num.genes[new.names]
            new.data.info[new.names,"orig.ident"]=as.character(object@ident[new.names])
            object@data.info=rbind(object@data.info,new.data.info)
            object=project.samples(object,new.data)
            #NOT SUPPORTED - adding mix.probs, gene.scores, etc.
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
#' using FetchData
#' @param ident.use Create a cell subset based on the provided identity classes
#' @param accept.low Low cutoff for the parameter (default is -Inf)
#' @param accept.high High cutoff for the parameter (default is Inf)
#' @param do.center Recenter the new object@@scale.data
#' @param do.scale Rescale the new object@@scale.data
#' @param \dots Additional arguments to be passed to FetchData (for example,
#' use.imputed=TRUE)
#' @return Returns a Seurat object containing only the relevant subset of cells
#' @export
setGeneric("SubsetData",  function(object,cells.use=NULL,subset.name=NULL,ident.use=NULL,accept.low=-Inf, accept.high=Inf,do.center=TRUE,do.scale=TRUE,...) standardGeneric("SubsetData"))
#' @export
setMethod("SubsetData","seurat",
          function(object,cells.use=NULL,subset.name=NULL,ident.use=NULL,accept.low=-Inf, accept.high=Inf,do.center=TRUE,do.scale=TRUE,...) {
            data.use=NULL
            if (!is.null(ident.use)) {
              cells.use=WhichCells(object,ident.use)
            }
            if (is.null(cells.use)) {
              data.use=FetchData(object,subset.name,...)
              if (length(data.use)==0) return(object)
              subset.data=data.use[,subset.name]
              pass.inds=which((subset.data>accept.low) & (subset.data<accept.high))
              cells.use=rownames(data.use)[pass.inds]
            }
            cells.use=ainb(cells.use,object@cell.names)
            object@data=object@data[,cells.use]
            if(!(is.null(object@scale.data))) {
              if (length(colnames(object@scale.data)>0)) {
                object@scale.data[,cells.use]
                object@scale.data=object@scale.data[complete.cases(object@scale.data),cells.use]
              }
            }
            if (do.scale) {
              object=ScaleData(object,do.scale = do.scale,do.center = do.center)
              object@scale.data=object@scale.data[complete.cases(object@scale.data),cells.use]
            }
            object@ident=drop.levels(object@ident[cells.use])
            object@tsne.rot=object@tsne.rot[cells.use,]
            object@pca.rot=object@pca.rot[cells.use,]
            object@cell.names=cells.use

            object@gene.scores=data.frame(object@gene.scores[cells.use,]); colnames(object@gene.scores)[1]="nGene"; rownames(object@gene.scores)=colnames(object@data)
            object@data.info=data.frame(object@data.info[cells.use,])
            #object@mix.probs=data.frame(object@mix.probs[cells.use,]); colnames(object@mix.probs)[1]="nGene"; rownames(object@mix.probs)=colnames(object@data)

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
#' using FetchData
#' @param accept.low Low cutoff for the parameter (default is -Inf)
#' @param accept.high High cutoff for the parameter (default is Inf)
#' @param do.center Recenter the new object@@scale.data
#' @param do.scale Rescale the new object@@scale.data
#' @param \dots Additional arguments to be passed to FetchData (for example,
#' use.imputed=TRUE)
#' @return Returns a Seurat object containing only the relevant subset of cells
#' @export
setGeneric("SubsetCells",  function(object,cells.use=NULL,subset.name=NULL,accept.low=-Inf, accept.high=Inf,do.center=TRUE,do.scale=TRUE,...) standardGeneric("SubsetCells"))
#' @export
setMethod("SubsetCells","seurat",
          function(object,cells.use=NULL,subset.name=NULL,accept.low=-Inf, accept.high=Inf,do.center=TRUE,do.scale=TRUE,...) {
            data.use=NULL
              data.use=FetchData(object,subset.name,cells.use,...)
              if (length(data.use)==0) return(object)
              subset.data=data.use[,subset.name]
              pass.inds=which((subset.data>accept.low) & (subset.data<accept.high))
              cells.use=rownames(data.use)[pass.inds]
            return(cells.use)
          }
)


#' @export
setGeneric("ProjectSamples", function(object,new.samples) standardGeneric("ProjectSamples"))
#' @export
setMethod("ProjectSamples", "seurat",
          function(object,new.samples) {
            genes.use = rownames(object@data)
            genes.pca = rownames(object@pca.x)
            data.project = object@scale.data[genes.pca,]; data.project[is.na(data.project)]=0;
            new.rot = t(data.project) %*% as.matrix(object@pca.x)
            scale.vec = apply(new.rot,2, function(x) sqrt(sum(x^2)))
            new.rot.scale = scale(new.rot, center=FALSE, scale=scale.vec)
            object@pca.rot = as.data.frame(new.rot.scale)
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
setGeneric("ProjectPCA", function(object,do.print=TRUE,pcs.print=5,pcs.store=30,genes.print=30,replace.pc=FALSE,do.center=FALSE) standardGeneric("ProjectPCA"))
#' @export
setMethod("ProjectPCA", "seurat",
          function(object,do.print=TRUE,pcs.print=5,pcs.store=30,genes.print=30,replace.pc=FALSE,do.center=FALSE) {
            if (!(do.center)) {
              genes.use=rownames(object@scale.data)
              #object@pca.x.full=data.frame(as.matrix(object@scale.data)%*%as.matrix(object@pca.rot))
              object.rot=as.matrix(object@pca.rot)
              object@pca.x.full <- data.frame(matrix(NA, nrow = length(genes.use), ncol = ncol(object@pca.rot)))
              rownames(object@pca.x.full) <- genes.use 
              colnames(object@pca.x.full) <- colnames(object@pca.rot)
             # dimnames(object@scale.data)=dimnames(data.use)

                bin.size <- 1000
                max.bin <- floor(length(genes.use)/bin.size) + 1
                pb <- txtProgressBar(min = 0, max = max.bin, style = 3)
                for(i in 1:max.bin) {
                  my.inds <- ((bin.size * (i - 1)):(bin.size * i - 1))+1
                  my.inds <- my.inds[my.inds <= length(genes.use)]
                  #print(my.inds)
                  object@pca.x.full[genes.use[my.inds], ] <- (object@scale.data[genes.use[my.inds], ])%*%object.rot
                  setTxtProgressBar(pb, i)  
                }
                close(pb)
              }
              
            if (do.center) object@pca.x.full=data.frame(scale(as.matrix(object@scale.data),center = TRUE,scale = FALSE)%*%as.matrix(object@pca.rot))
            if (ncol(object@jackStraw.fakePC)>0) {
              object@jackStraw.empP.full=data.frame(sapply(1:ncol(object@jackStraw.fakePC),function(x)unlist(lapply(abs(object@pca.x.full[,x]),empP,abs(object@jackStraw.fakePC[,x])))))
              colnames(object@jackStraw.empP.full)=paste("PC",1:ncol(object@jackStraw.empP),sep="")
              rownames(object@jackStraw.empP.full)=rownames(object@scale.data)
            }
            object@pca.x.full[is.na(object@pca.x.full)]=0

            if (replace.pc==TRUE) {
              object@jackStraw.empP=object@jackStraw.empP.full
              object@pca.x=object@pca.x.full
            }

            if (do.print) {
                PrintPCA(object,1:pcs.print,genes.print,TRUE)
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
setGeneric("RunTSNE", function(object,cells.use=NULL,dims.use=1:5,k.seed=1,do.fast=FALSE,add.iter=0,genes.use=NULL,reduction.use="pca",dim_embed=2,...) standardGeneric("RunTSNE"))
#' @export
setMethod("RunTSNE", "seurat",
          function(object,cells.use=NULL,dims.use=1:5,k.seed=1,do.fast=FALSE,add.iter=0,genes.use=NULL,reduction.use="pca",dim_embed=2,...) {
            cells.use=set.ifnull(cells.use,colnames(object@data))
            if (is.null(genes.use)) {
              dim.code=translate.dim.code(reduction.use); dim.codes=paste(dim.code,dims.use,sep="")
              data.use=FetchData(object,dim.codes)
            }
            if (!is.null(genes.use)) {
              genes.use=ainb(genes.use,rownames(object@data))
              data.use=t(object@data[genes.use,cells.use])
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
setGeneric("RunDiffusion", function(object,cells.use=NULL,dims.use=1:5,k.seed=1,do.fast=FALSE,add.iter=0,genes.use=NULL,reduction.use="pca",dim_embed=2,q.use=0.05,max.dim=2,scale.clip=3,...) standardGeneric("RunDiffusion"))
#' @export
setMethod("RunDiffusion", "seurat",
          function(object,cells.use=NULL,dims.use=1:5,k.seed=1,do.fast=FALSE,add.iter=0,genes.use=NULL,reduction.use="pca",dim_embed=2,q.use=0.05,max.dim=2,scale.clip=3,...) {
            cells.use=set.ifnull(cells.use,colnames(object@data))
            if (is.null(genes.use)) {
              dim.code=translate.dim.code(reduction.use); dim.codes=paste(dim.code,dims.use,sep="")
              data.use=FetchData(object,dim.codes)
            }
            if (!is.null(genes.use)) {
              genes.use=ainb(genes.use,rownames(object@scale.data))
              data.use=minmax(t(object@scale.data[genes.use,cells.use]),-1*scale.clip,scale.clip)
            }
            data.dist=dist(data.use)
            data.diffusion=data.frame(diffuse(data.dist,neigen = max.dim,maxdim = max.dim,...)$X)
            colnames(data.diffusion)=paste("tSNE_",1:ncol(data.diffusion),sep="")
            rownames(data.diffusion)=cells.use
            for(i in 1:max.dim) {
              x=data.diffusion[,i]; x=minmax(x,min = quantile(x,q.use),quantile(x,1-q.use)); data.diffusion[,i]=x
            }
            object@tsne.rot=data.diffusion
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
setGeneric("ICA", function(object,ic.genes=NULL,do.print=TRUE,ics.print=5,ics.store=30,genes.print=30,use.imputed=FALSE,seed.use=1,...) standardGeneric("ICA"))
#' @export
setMethod("ICA", "seurat",
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
#' @param rev.pca By default computes the PCA on the cell x gene matrix. Setting to true will compute it on gene x cell matrix. 
#' @param \dots Additional arguments to be passed to prcomp
#' @return Returns Seurat object with an PCA embedding (object@@pca.rot) and
#' gene projection matrix (object@@pca.x). The PCA object itself is stored in
#' object@@pca.obj[[1]]
#' @export
setGeneric("PCA", function(object,pc.genes=NULL,do.print=TRUE,pcs.print=5,pcs.store=40,genes.print=30,use.imputed=FALSE,rev.pca=FALSE,...) standardGeneric("PCA"))
#' @export
setMethod("PCA", "seurat",
          function(object,pc.genes=NULL,do.print=TRUE,pcs.print=5,pcs.store=40,genes.print=30,use.imputed=FALSE,rev.pca=FALSE,...) {
            data.use=object@scale.data
            if (use.imputed) data.use=data.frame(t(scale(t(object@imputed))))
            pc.genes=set.ifnull(pc.genes,object@var.genes)
            pc.genes = unique(pc.genes[pc.genes%in%rownames(data.use)])
            pc.genes.var = apply(data.use[pc.genes,],1,var)

            if (rev.pca) {
              pc.genes.use=pc.genes[pc.genes.var>0]; pc.genes.use=pc.genes.use[!is.na(pc.genes.use)]
              pc.data = data.use[pc.genes.use,]
              pca.obj = prcomp(pc.data,...)
              object@pca.obj=list(pca.obj)

              pcs.store=min(pcs.store,ncol(pc.data))
              pcs.print=min(pcs.print,ncol(pc.data))
              object@pca.x=data.frame(pca.obj$x[,1:pcs.store])
              object@pca.rot=data.frame(pca.obj$rotation[,1:pcs.store])
            }
            if (!rev.pca) {
              pc.genes.use=pc.genes[pc.genes.var>0]; pc.genes.use=pc.genes.use[!is.na(pc.genes.use)]
              pc.data = data.use[pc.genes.use,]
              pca.obj = prcomp(t(pc.data),...)
              object@pca.obj=list(pca.obj)
  
              pcs.store=min(pcs.store,ncol(pc.data))
              pcs.print=min(pcs.print,ncol(pc.data))
              object@pca.x=data.frame(pca.obj$rotation[,1:pcs.store])
              object@pca.rot=data.frame(pca.obj$x[,1:pcs.store])
              
            }

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

#' Run Principal Component Analysis on gene expression using IRLBA
#'
#' Run Fast PCA dimensionality reduction
#'
#' @param object Seurat object
#' @param pc.genes Genes to use as input for PCA. Default is object@@var.genes
#' @param do.print Print the top genes associated with high/low loadings for
#' the PCs
#' @param pcs.print Number of PCs to print genes for
#' @param pcs.compute Total Number of PCs to compute and store
#' @param genes.print Number of genes to print for each PC
#' @param \dots Additional arguments to be passed to prcomp
#' @return Returns Seurat object with an PCA embedding (object@@pca.rot) and
#' gene projection matrix (object@@pca.x). The PCA object itself is stored in
#' object@@pca.obj[[1]]
#' @export
setGeneric("PCAFast", function(object,pc.genes=NULL,do.print=TRUE,pcs.print=5,pcs.compute=20,genes.print=30,...) standardGeneric("PCAFast"))
#' @export
setMethod("PCAFast", "seurat",
          function(object,pc.genes=NULL,do.print=TRUE,pcs.print=5,pcs.compute=20,genes.print=30,...) {
            data.use=object@scale.data
            pc.genes=set.ifnull(pc.genes,object@var.genes)
            pc.genes = unique(pc.genes[pc.genes%in%rownames(data.use)])
            pc.genes.var = apply(data.use[pc.genes,],1,var)

              pc.genes.use=pc.genes[pc.genes.var>0]; pc.genes.use=pc.genes.use[!is.na(pc.genes.use)]
              pc.data = data.use[pc.genes.use,]
              pca.obj = irlba(t(pc.data),nv = pcs.compute,...)
              object@pca.obj=list(pca.obj)
              
              pcs.store=min(pcs.compute,ncol(pc.data))
              pcs.print=min(pcs.print,ncol(pc.data))
              object@pca.x=data.frame(pca.obj$v[,1:pcs.store],row.names = rownames(pc.data)); colnames(object@pca.x)=paste("PC",1:pcs.compute,sep="")
              object@pca.rot=data.frame(pca.obj$u[,1:pcs.store],row.names = colnames(pc.data)); colnames(object@pca.rot)=colnames(object@pca.x)
              
            
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
setGeneric("ClusterAlpha", function(object,thresh.min=0) standardGeneric("ClusterAlpha"))
#' @export
setMethod("ClusterAlpha", "seurat",
          function(object,thresh.min=0) {
            ident.use=object@ident
            data.all=data.frame(row.names = rownames(object@data))
            for(i in sort(unique(ident.use))) {
              temp.cells=WhichCells(object,i)
              data.temp=apply(object@data[,temp.cells],1,function(x)return(length(x[x>thresh.min])/length(x)))
              data.all=cbind(data.all,data.temp)
              colnames(data.all)[ncol(data.all)]=i
            }
            colnames(data.all)=sort(unique(ident.use))
            return(data.all)
          }
)


#' Reorder identity classes
#'
#' Re-assigns the identity classes according to the average expression of a particular feature (i.e, gene expression, or PC score)
#' Very useful after clustering, to re-order cells, for example, based on PC scores
#'
#' @param object Seurat object
#' @param feature Feature to reorder on. Default is PC1
#' @param rev Reverse ordering (default is FALSE)
#' @param aggregate.fxn Function to evaluate each identity class based on (default is mean)
#' @param reorder.numeric Rename all identity classes to be increasing numbers starting from 1 (default is FALSE)
#' @param \dots additional arguemnts (i.e. use.imputed=TRUE)
#' @return A seurat object where the identity have been re-oredered based on the average.
#' @export
setGeneric("ReorderIdent", function(object,feature="PC1",rev=FALSE,aggregate.fxn=mean,reorder.numeric=FALSE,...) standardGeneric("ReorderIdent"))
#' @export
setMethod("ReorderIdent", "seurat",
          function(object,feature="PC1",rev=FALSE,aggregate.fxn=mean,reorder.numeric=FALSE,...) {
            ident.use=object@ident
            data.use=FetchData(object,feature,...)[,1]
            revFxn=same; if (rev) revFxn=function(x)max(x)+1-x;
            names.sort=names(revFxn(sort(tapply(data.use,(ident.use),aggregate.fxn))))
            ident.new=factor(ident.use,levels = names.sort,ordered = TRUE)
            if (reorder.numeric)  ident.new=factor(revFxn(rank(tapply(data.use,as.numeric(ident.new),mean)))[as.numeric(ident.new)],levels = 1:length(levels(ident.new)),ordered = TRUE)
            names(ident.new)=names(ident.use)
            object@ident=ident.new
            return(object)
          }
)

#' Average PCA scores by identity class
#'
#' Returns the PCA scores for an 'average' single cell in each identity class
#'
#' @param object Seurat object
#' @return Returns a matrix with genes as rows, identity classes as columns
#' @export
setGeneric("AveragePCA", function(object) standardGeneric("AveragePCA"))
#' @export
setMethod("AveragePCA", "seurat",
          function(object) {
            ident.use=object@ident
            data.all=data.frame(row.names = colnames(object@pca.rot))
            for(i in levels(ident.use)) {
              temp.cells=WhichCells(object,i)
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
setGeneric("AverageExpression", function(object,genes.use=NULL,return.seurat=F,add.ident=NULL,...) standardGeneric("AverageExpression"))
#' @export
setMethod("AverageExpression", "seurat",
          function(object,genes.use=NULL,return.seurat=F,add.ident=NULL,...) {
            genes.use=set.ifnull(genes.use,rownames(object@data))
            genes.use=ainb(genes.use,rownames(object@data))
            ident.use=object@ident
            if (!is.null(add.ident)) {
              new.data=FetchData(object,add.ident)
              new.ident=paste(object@ident[rownames(new.data)],new.data[,1],sep="_")
              object=SetIdent(object,rownames(new.data),new.ident)
            }
            data.all=data.frame(row.names = genes.use)
            for(i in levels(object@ident)) {
              temp.cells=WhichCells(object,i)
              if (length(temp.cells)==1) data.temp=(object@data[genes.use,temp.cells])
              if (length(temp.cells)>1) data.temp=apply(object@data[genes.use,temp.cells],1,expMean)
              data.all=cbind(data.all,data.temp)
              colnames(data.all)[ncol(data.all)]=i
            }
            colnames(data.all)=levels(object@ident)
            if (return.seurat) {
              toRet=new("seurat",raw.data=data.all)
              toRet=Setup(toRet,project = "Average",min.cells = 0,min.genes = 0,is.expr = 0,...)
              return(toRet)
            }
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
setGeneric("ICTopGenes", function(object,ic.use=1,num.genes=30,do.balanced=FALSE) standardGeneric("ICTopGenes"))
#' @export
setMethod("ICTopGenes", "seurat",
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
setGeneric("PCTopGenes", function(object,pc.use=1,num.genes=30,use.full=FALSE,do.balanced=FALSE) standardGeneric("PCTopGenes"))
#' @export
setMethod("PCTopGenes", "seurat",
          function(object,pc.use=1,num.genes=30,use.full=FALSE,do.balanced=FALSE) {
            pc_scores=object@pca.x
            i=pc.use
            if (use.full==TRUE) pc_scores = object@pca.x.full
            pc.top.genes=unique(unlist(lapply(i,topGenesForDim,pc_scores,do.balanced,num.genes,"pca")))
            return(pc.top.genes)
          }
)

#' Find cells with highest PCA scores
#'
#' Return a list of genes with the strongest contribution to a set of principal components
#'
#' @param object Seurat object
#' @param pc.use Principal component to use
#' @param num.genes Number of cells to return
#' @param do.balanced Return an equal number of cells with both + and - PC scores.
#' @return Returns a vector of cells
#' @export
setGeneric("PCTopCells",  function(object,pc.use=1,num.cells=NULL,do.balanced=FALSE) standardGeneric("PCTopCells"))
#' @export
setMethod("PCTopCells", "seurat",
          function(object,pc.use=1,num.cells=NULL,do.balanced=FALSE) {

            #note that we use topGenesForDim, but it still works
            num.cells=set.ifnull(num.cells,length(object@cell.names))
            pc_scores=object@pca.rot
            i=pc.use
            pc.top.cells=unique(unlist(lapply(i,topGenesForDim,pc_scores,do.balanced,num.cells,"pca")))
            return(pc.top.cells)
          }
)

#' Print the results of a PCA analysis
#'
#' Prints a set of genes that most strongly define a set of principal components
#'
#' @inheritParams VizPCA
#' @param pcs.print Set of PCs to print genes for
#' @param genes.print Number of genes to print for each PC
#' @return Only text output
#' @export
setGeneric("PrintPCA", function(object,pcs.print=1:5,genes.print=30,use.full=FALSE) standardGeneric("PrintPCA"))
#' @export
setMethod("PrintPCA", "seurat",
          function(object,pcs.print=1:5,genes.print=30,use.full=FALSE) {
            genes.print.use=round(genes.print/2)
            for(i in pcs.print) {
              code=paste("PC",i,sep="")
              sx=PCTopGenes(object,i,genes.print.use*2,use.full,TRUE)
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
setGeneric("FetchData",  function(object, vars.all=NULL,cells.use=NULL,use.imputed=FALSE, use.scaled=FALSE) standardGeneric("FetchData"))
#' @export
setMethod("FetchData","seurat",
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
                else data.expression = object@data[vars.all, , drop = FALSE ]
                return(t(as.matrix(data.expression)))
              }
              else{
                if(use.imputed) data.expression = object@imputed[vars.all[gene_check], ]
                if(use.scaled) data.expression = object@scale.data[vars.all[gene_check], ]
                else data.expression = object@data[vars.all[gene_check], , drop = FALSE]
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
                stop(paste("Error : ", my.var, " not found", sep=""))
              }
              cells.use=ainb(cells.use,rownames(data.use))
              data.add=data.use[cells.use,my.var]
              if (is.null(data.add)) {
                stop(paste("Error : ", my.var, " not found", sep=""))
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
setGeneric("VizICA", function(object,ics.use=1:5,num.genes=30,font.size=0.5,nCol=NULL,do.balanced=FALSE) standardGeneric("VizICA"))
#' @export
setMethod("VizICA", "seurat",
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
              subset.use=ic_scores[ICTopGenes(object,i,num.genes,do.balanced),]
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
setGeneric("VizPCA", function(object,pcs.use=1:5,num.genes=30,use.full=FALSE,font.size=0.5,nCol=NULL,do.balanced=FALSE) standardGeneric("VizPCA"))
#' @export
setMethod("VizPCA", "seurat",
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
              subset.use=pc_scores[PCTopGenes(object,i,num.genes,use.full,do.balanced),]
              plot(subset.use[,i],1:nrow(subset.use),pch=16,col="blue",xlab=paste("PC",i,sep=""),yaxt="n",ylab="")
              axis(2,at=1:nrow(subset.use),labels = rownames(subset.use),las=1,cex.axis=font.size)
            }
            rp()
          }
)


setGeneric("GetWeightMatrix", function(object) standardGeneric("GetWeightMatrix"))
setMethod("GetWeightMatrix", "seurat",
          function(object) {
            data=object@data
            data.humpAvg=apply(data,1,humpMean,min=object@drop.expr)
            wt.matrix=data.frame(t(pbsapply(data.humpAvg,expAlpha,object@drop.coefs)))
            colnames(wt.matrix)=colnames(data); rownames(wt.matrix)=rownames(data)
            wt.matrix[is.na(wt.matrix)]=0
            object@wt.matrix=wt.matrix
            wt1.matrix=data.frame(pbsapply(1:ncol(data),function(x)setWt1(data[,x],wt.matrix[,x],min=object@drop.expr)))
            colnames(wt1.matrix)=colnames(data); rownames(wt1.matrix)=rownames(data)
            wt1.matrix[is.na(wt1.matrix)]=0
            object@drop.wt.matrix=wt1.matrix
            return(object)
          }
)


setGeneric("RegulatorScore", function(object, candidate.reg, score.name, cells.use=NULL) standardGeneric("RegulatorScore"))
setMethod("RegulatorScore", "seurat",
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
#' @inheritParams FindMarkers
#' @param node The node in the phylogenetic tree to use as a branch point
#' @return Matrix containing a ranked list of putative markers, and associated
#' statistics (p-values, ROC score, etc.)
#' @export
setGeneric("FindMarkersNode", function(object,node,genes.use=NULL,thresh.use=log(2),test.use="bimod",...) standardGeneric("FindMarkersNode"))
#' @export
setMethod("FindMarkersNode", "seurat",
          function(object,node,genes.use=NULL,thresh.use=log(2),test.use="bimod",...) {
            genes.use=set.ifnull(genes.use,rownames(object@data))
            tree=object@cluster.tree[[1]]
            ident.order=tree$tip.label
            nodes.1=ident.order[getLeftDecendants(tree,node)]
            nodes.2=ident.order[getRightDecendants(tree,node)]
            #print(nodes.1)
            #print(nodes.2)
            to.return=FindMarkers(object,nodes.1,nodes.2,genes.use,thresh.use,test.use,...)
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
#' @param min.pct - only test genes that are detected in a minimum fraction of min.pct cells
#' in either of the two populations. Meant to speed up the function by not testing genes that are very infrequently expression
#' @param only.pos Only return positive markers (FALSE by default)
#' @return Matrix containing a ranked list of putative markers, and associated statistics (p-values, ROC score, etc.)
#' @param print.bar Print a progress bar once expression testing begins (uses pbapply to do this)
#' @import VGAM
#' @import pbapply
#' @export
setGeneric("FindMarkers", function(object, ident.1,ident.2=NULL,genes.use=NULL,thresh.use=log(2),test.use="bimod",min.pct=0,print.bar=TRUE,only.pos=FALSE) standardGeneric("FindMarkers"))
#' @export
setMethod("FindMarkers", "seurat",
          function(object, ident.1,ident.2=NULL,genes.use=NULL,thresh.use=log(2), test.use="bimod",min.pct=0,print.bar=TRUE,only.pos=FALSE) {
            genes.use=set.ifnull(genes.use,rownames(object@data))

            cells.1=WhichCells(object,ident.1)
            # in case the user passed in cells instead of identity classes
            if (length(ident.1>1)&&any(ident.1%in%object@cell.names)) {
              cells.1=ainb(ident.1,object@cell.names)
            }

            # if NULL for ident.2, use all other cells
            if (is.null(ident.2)) {
              cells.2=object@cell.names
            }
            else {
              cells.2=WhichCells(object,ident.2)
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
            thresh.min=object@is.expr
            data.temp1=round(apply(object@data[genes.use,cells.1],1,function(x)return(length(x[x>thresh.min])/length(x))),3)
            data.temp2=round(apply(object@data[genes.use,cells.2],1,function(x)return(length(x[x>thresh.min])/length(x))),3)
            data.alpha=cbind(data.temp1,data.temp2); colnames(data.alpha)=c("pct.1","pct.2")
            alpha.min=apply(data.alpha,1,max); names(alpha.min)=rownames(data.alpha); genes.use=names(which(alpha.min>min.pct))

            if (test.use=="bimod") to.return=DiffExpTest(object,cells.1,cells.2,genes.use,thresh.use,print.bar)
            if (test.use=="roc") to.return=MarkerTest(object,cells.1,cells.2,genes.use,thresh.use,print.bar)
            if (test.use=="t") to.return=DiffTTest(object,cells.1,cells.2,genes.use,thresh.use,print.bar)
            if (test.use=="tobit") to.return=TobitTest(object,cells.1,cells.2,genes.use,thresh.use,print.bar)
            to.return=cbind(to.return,data.alpha[rownames(to.return),])
            if(only.pos) to.return=subset(to.return,avg_diff>0)
            return(to.return)
          }
)


#' Gene expression markers for all identity classes
#'
#' Finds markers (differentially expressed genes) for each of the identity classes in a dataset
#'
#'
#' @inheritParams FindMarkers
#' @param thresh.test Limit testing to genes which show, on average, at least X-fold difference (log-scale) between cells in an identity class, and all other cells.
#' @param return.thresh Only return markers that have a p-value < return.thresh, or a power > return.thresh (if the test is ROC)
#' @param do.print FALSE by default. If TRUE, outputs updates on progress.
#' @return Matrix containing a ranked list of putative markers, and associated
#' statistics (p-values, ROC score, etc.)
#' @export
setGeneric("FindAllMarkers", function(object, thresh.test=1,test.use="bimod",return.thresh=1e-2,do.print=FALSE,min.pct=0,print.bar=TRUE,only.pos=FALSE) standardGeneric("FindAllMarkers"))
#' @export
setMethod("FindAllMarkers","seurat",
      function(object, thresh.test=1,test.use="bimod",return.thresh=1e-2,do.print=FALSE,min.pct=0,print.bar=TRUE,only.pos=FALSE) {
            ident.use=object@ident
            if ((test.use=="roc") && (return.thresh==1e-2)) return.thresh=0.7
            idents.all=sort(unique(object@ident))
            genes.de=list()
            for(i in 1:length(idents.all)) {
              genes.de[[i]]=FindMarkers(object,ident.1 = idents.all[i],ident.2 = NULL,genes.use=rownames(object@data),thresh.use = thresh.test,test.use = test.use,min.pct,print.bar)
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
            if(only.pos) return(subset(gde.all,avg_diff>0))
            return(gde.all)
          }
)


#' Likelihood ratio test for zero-inflated data
#'
#' Identifies differentially expressed genes between two groups of cells using
#' the LRT model proposed in Mcdavid et al, Bioinformatics, 2011
#'
#'
#' @inheritParams FindMarkers
#' @param object Seurat object
#' @param cells.1 Group 1 cells
#' @param cells.2 Group 2 cells
#' @return Returns a p-value ranked matrix of putative differentially expressed
#' genes.
#' @export
setGeneric("DiffExpTest", function(object, cells.1,cells.2,genes.use=NULL,thresh.use=log(2),print.bar=TRUE) standardGeneric("DiffExpTest"))
#' @export
setMethod("DiffExpTest", "seurat",
          function(object, cells.1,cells.2,genes.use=NULL,thresh.use=log(2),print.bar=TRUE) {
            genes.use=set.ifnull(genes.use,object@var.genes)
            data.1=apply(object@data[genes.use,cells.1],1,expMean)
            data.2=apply(object@data[genes.use,cells.2],1,expMean)
            total.diff=abs(data.1-data.2)
            genes.diff = names(which(total.diff>thresh.use))
            #print(genes.diff)
            to.return=BimodDiffExpTest(object@data[,cells.1],object@data[,cells.2],genes.diff,print.bar)
            to.return=to.return[order(to.return$p_val,-abs(to.return$avg_diff)),]
            return(to.return)
          }
)

#' Differential expression testing using Tobit models
#'
#' Identifies differentially expressed genes between two groups of cells using
#' Tobit models, as proposed in Trapnell et al., Nature Biotechnology, 2014
#'
#' @inheritParams FindMarkers
#' @inheritParams DiffExpTest
#' @return Returns a p-value ranked matrix of putative differentially expressed
#' genes.
#' @export
setGeneric("TobitTest", function(object, cells.1,cells.2,genes.use=NULL,thresh.use=log(2),print.bar=TRUE) standardGeneric("TobitTest"))
#' @export
setMethod("TobitTest", "seurat",
          function(object, cells.1,cells.2,genes.use=NULL,thresh.use=log(2),print.bar=TRUE) {
            genes.use=set.ifnull(genes.use,object@var.genes)
            data.1=apply(object@data[genes.use,cells.1],1,expMean)
            data.2=apply(object@data[genes.use,cells.2],1,expMean)
            total.diff=abs(data.1-data.2)
            genes.diff = names(which(total.diff>thresh.use))
            #print(genes.diff)
            to.return=TobitDiffExpTest(object@data[,cells.1],object@data[,cells.2],genes.diff,print.bar)
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
setGeneric("BatchGene", function(object, idents.use,genes.use=NULL,auc.cutoff=0.6,thresh.use=0) standardGeneric("BatchGene"))
#' @export
setMethod("BatchGene", "seurat",
          function(object, idents.use,genes.use=NULL,auc.cutoff=0.6,thresh.use=0) {
            batch.genes=c()
            genes.use=set.ifnull(genes.use,rownames(object@data))
            for(ident in idents.use ) {
              cells.1=names(object@ident)[object@ident==ident]
              cells.2=names(object@ident)[object@ident!=ident]
              if ((length(cells.1)<5)|(length(cells.2)<5)) {
                break;
              }
              markers.ident=MarkerTest(object,cells.1,cells.2,genes.use,thresh.use)
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
#' @inheritParams FindMarkers
#' @inheritParams DiffExpTest
#' @param object Seurat object
#' @return Returns a 'predictive power' (abs(AUC-0.5)) ranked matrix of
#' putative differentially expressed genes.
#' @import ROCR
#' @export
setGeneric("MarkerTest", function(object, cells.1,cells.2,genes.use=NULL,thresh.use=log(2),print.bar=TRUE) standardGeneric("MarkerTest"))
#' @export
setMethod("MarkerTest", "seurat",
          function(object, cells.1,cells.2,genes.use=NULL,thresh.use=log(2),print.bar=TRUE) {
            genes.use=set.ifnull(genes.use,object@var.genes)
            data.use=object@data
            data.1=apply(object@data[genes.use,cells.1],1,expMean)
            data.2=apply(object@data[genes.use,cells.2],1,expMean)
            total.diff=abs(data.1-data.2)
            genes.diff = names(which(total.diff>thresh.use))
            genes.use=ainb(genes.diff,rownames(data.use))
            to.return=marker.auc.test(object@data[,cells.1],object@data[,cells.2],genes.use,print.bar=TRUE)
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
#' @inheritParams FindMarkers
#' @inheritParams DiffExpTest
#' @return Returns a p-value ranked matrix of putative differentially expressed
#' genes.
#' @export
setGeneric("DiffTTest", function(object, cells.1,cells.2,genes.use=NULL,thresh.use=log(2),print.bar=TRUE) standardGeneric("DiffTTest"))
#' @export
setMethod("DiffTTest", "seurat",
          function(object, cells.1,cells.2,genes.use=NULL,thresh.use=log(2),print.bar=TRUE) {
            genes.use=set.ifnull(genes.use,object@var.genes)
            data.use=object@data
            data.1=apply(object@data[genes.use,cells.1],1,expMean)
            data.2=apply(object@data[genes.use,cells.2],1,expMean)
            total.diff=abs(data.1-data.2)
            genes.diff = names(which(total.diff>thresh.use))
            genes.use=ainb(genes.diff,rownames(data.use))
            iterate.fxn=lapply; if (print.bar) iterate.fxn=pblapply
            p_val=unlist(iterate.fxn(genes.use,function(x)t.test(object@data[x,cells.1],object@data[x,cells.2])$p.value))
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
setGeneric("WhichCells", function(object,value=1, id=NULL) standardGeneric("WhichCells"))
#' @export
setMethod("WhichCells", "seurat",
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
setGeneric("SetAllIdent", function(object,id=NULL) standardGeneric("SetAllIdent"))
#' @export
setMethod("SetAllIdent", "seurat",
          function(object, id=NULL) {
            id=set.ifnull(id,"orig.ident")
            if (id %in% colnames(object@data.info)) {
              cells.use=rownames(object@data.info)
              ident.use=object@data.info[,id]
              object=SetIdent(object,cells.use,ident.use)
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
setGeneric("RenameIdent", function(object,old.ident.name=NULL,new.ident.name=NULL) standardGeneric("RenameIdent"))
#' @export
setMethod("RenameIdent", "seurat",
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
            ident.vector[WhichCells(object,old.ident.name)]=new.ident.name
            object@ident=factor(ident.vector,levels = new.levels)
            return(object)
          }
)

#' Set identity class information
#'
#' Stashes the identity in data.info to be retrieved later. Useful if, for example, testing multiple clustering parameters
#'
#' @param object Seurat object
#' @param save.name Store current object@@ident under this column name in object@@data.info. Can be easily retrived with SetAllIdent
#' @return A Seurat object where object@@ident has been appropriately modified
#' @export
setGeneric("StashIdent", function(object,save.name="oldIdent") standardGeneric("StashIdent"))
#' @export
setMethod("StashIdent", "seurat",
          function(object,save.name="oldIdent") {
            object@data.info[,save.name]=as.character(object@ident)
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
setGeneric("SetIdent", function(object,cells.use=NULL,ident.use=NULL) standardGeneric("SetIdent"))
#' @export
setMethod("SetIdent", "seurat",
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
setGeneric("PosteriorPlot", function(object, name) standardGeneric("PosteriorPlot"))
#' @export
setMethod("PosteriorPlot", "seurat",
          function(object, name) {
            post.names=colnames(subc(object@mix.probs,name))
            VlnPlot(object,post.names,inc.first=TRUE,inc.final=TRUE,by.k=TRUE)


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
setGeneric("MapCell",  function(object,cell.name,do.plot=FALSE,safe.use=TRUE,text.val=NULL,do.rev=FALSE) standardGeneric("MapCell"))
#' @export
setMethod("MapCell", "seurat",
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
setGeneric("GetCentroids", function(object, cells.use=NULL,get.exact=TRUE) standardGeneric("GetCentroids"))
#' @export
setMethod("GetCentroids", "seurat",
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
setGeneric("RefinedMapping",  function(object,genes.use) standardGeneric("RefinedMapping"))
#' @export
setMethod("RefinedMapping", "seurat",
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
setGeneric("InitialMapping", function(object,cells.use=NULL) standardGeneric("InitialMapping"))
#' @export
setMethod("InitialMapping", "seurat",
          function(object,cells.use=NULL) {
            cells.use=set.ifnull(cells.use,colnames(object@data))
            every.prob=sapply(cells.use,function(x)MapCell(object,x,FALSE,FALSE))
            object@final.prob=data.frame(every.prob)
            rownames(object@final.prob)=paste("bin.",rownames(object@final.prob),sep="")
            return(object)
          }
)

#Internal, not documented for now
setGeneric("CalcInsitu", function(object,gene,do.plot=TRUE,do.write=FALSE,write.dir="~/window/insitu/",col.use=bwCols,do.norm=FALSE,cells.use=NULL,do.return=FALSE,probs.min=0,do.log=FALSE, use.imputed=FALSE, bleach.use=0) standardGeneric("CalcInsitu"))
setMethod("CalcInsitu", "seurat",
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
setGeneric("FitGeneK", function(object, gene, do.k=2,num.iter=1,do.plot=FALSE,genes.use=NULL,start.pct=NULL) standardGeneric("FitGeneK"))
#' @export
setMethod("FitGeneK", "seurat",
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
setGeneric("FitGeneMix", function(object, gene, do.k=3,use.mixtools=TRUE,do.plot=FALSE,plot.with.imputed=TRUE,min.bin.size=10) standardGeneric("FitGeneMix"))
#' @export
setMethod("FitGeneMix", "seurat",
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

#' Calculate smoothed expression values
#'
#'
#' Smooths expression values across the k-nearest neighbors based on dimensional reduction
#'
#' @inheritParams FeaturePlot
#' @inheritParams AddImputedScore
#' @param inheritParams Seurat object
#' @param genes.fit Genes to calculate smoothed values for
#' @importFrom FNN get.knn
#' @export
setGeneric("AddSmoothedScore", function(object,genes.fit=NULL,dim.1=1,dim.2=2,reduction.use="tsne",k=30,do.log=FALSE,do.print=FALSE) standardGeneric("AddSmoothedScore"))
#' @export
setMethod("AddSmoothedScore", "seurat",
          function(object,genes.fit=NULL,dim.1=1,dim.2=2,reduction.use="tSNE",k=30,do.log=FALSE,do.print=FALSE) {
            genes.fit=set.ifnull(genes.fit,object@var.genes)
            genes.fit=genes.fit[genes.fit%in%rownames(object@data)]

            dim.code=translate.dim.code(reduction.use); dim.codes=paste(dim.code,c(dim.1,dim.2),sep="")
            data.plot=FetchData(object,dim.codes)
            knn.smooth=get.knn(data.plot,k)$nn.index
            avg.fxn=mean;
            if (do.log==FALSE) avg.fxn=expMean;
            lasso.fits=data.frame(t(sapply(genes.fit,function(g) unlist(lapply(1:nrow(data.plot),function(y) avg.fxn(as.numeric(object@data[g,knn.smooth[y,]])))))))
            colnames(lasso.fits)=rownames(data.plot)
            genes.old=genes.fit[genes.fit%in%rownames(object@imputed)]
            genes.new=genes.fit[!(genes.fit%in%rownames(object@imputed))]

            if (length(genes.old)>0) object@imputed[genes.old,]=lasso.fits[genes.old,]
            object@imputed=rbind(object@imputed,lasso.fits[genes.new,])
            return(object)
          }
)
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
setGeneric("AddImputedScore", function(object, genes.use=NULL,genes.fit=NULL,s.use=20,do.print=FALSE,gram=TRUE) standardGeneric("AddImputedScore"))
#' @export
setMethod("AddImputedScore", "seurat",
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
#' @export
setGeneric("GetNewScore", function(object, score.name,score.genes, cell.ids=NULL, score.func=weighted.mean,scramble=FALSE, no.tech.wt=FALSE, biol.wts=NULL,use.scaled=FALSE) standardGeneric("GetNewScore"))
#' @export
setMethod("GetNewScore", "seurat",
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
setGeneric("CalcNoiseModels", function(object, cell.ids=NULL, trusted.genes=NULL,n.bin=20,drop.expr=1) standardGeneric("CalcNoiseModels"))
#' @export
setMethod("CalcNoiseModels","seurat",
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
            my.coefs=data.frame(t(pbsapply(colnames(data[1:(ncol(data)-2)]),
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
#'
#' @param object Seurat object
#' @param features.plot Vector of features to plot
#' @param dim.1 Dimension for x-axis (default 1)
#' @param dim.2 Dimension for y-axis (default 2)
#' @param cells.use Vector of cells to plot (default is all cells)
#' @param pt.size Adjust point size for plotting
#' @param cols.use The two colors to form the gradient over. Provide as string vector with
#' the first color corresponding to low values, the second to high. Also accepts a Brewer
#' color scale or vector of colors. Note: this will bin the data into number of colors provided.
#' @param pch.use Pch for plotting
#' @param reduction.use Which dimensionality reduction to use. Default is
#' "tsne", can also be "pca", or "ica", assuming these are precomputed.
#' @param use.imputed Use imputed values for gene expression (default is FALSE)
#' @param nCol Number of columns to use when plotting multiple features.
#' @param no.axes Remove axis labels
#' @importFrom RColorBrewer brewer.pal.info
#' @return No return value, only a graphical output
#' @export
setGeneric("FeaturePlot", function(object,features.plot,dim.1=1,dim.2=2,cells.use=NULL,pt.size=1,cols.use=c("yellow", "red"), pch.use=16,reduction.use="tsne",use.imputed=FALSE,nCol=NULL,no.axes=FALSE,...) standardGeneric("FeaturePlot"))
#' @export
setMethod("FeaturePlot", "seurat",
          function(object,features.plot,dim.1=1,dim.2=2,cells.use=NULL,pt.size=1,cols.use=c("yellow", "red"), pch.use=16,reduction.use="tsne",use.imputed=FALSE,nCol=NULL,no.axes=FALSE,...) {
            cells.use=set.ifnull(cells.use,colnames(object@data))
            dim.code="PC"
            if (is.null(nCol)) {
              nCol=2
              if (length(features.plot) == 1) ncol=1
              if (length(features.plot)>6) nCol=3
              if (length(features.plot)>9) nCol=4
            }
            num.row=floor(length(features.plot)/nCol-1e-5)+1
            par(mfrow=c(num.row,nCol))
            dim.code=translate.dim.code(reduction.use); dim.codes=paste(dim.code,c(dim.1,dim.2),sep="")
            data.plot=FetchData(object,dim.codes)

            x1=paste(dim.code,dim.1,sep=""); x2=paste(dim.code,dim.2,sep="")
            
            data.plot$x=data.plot[,x1]; data.plot$y=data.plot[,x2]
            data.plot$pt.size=pt.size
            data.use=data.frame(t(FetchData(object,features.plot,cells.use = cells.use,use.imputed = use.imputed)))
            pList=lapply(features.plot,function(x) SingleFeaturePlot(data.use, x, data.plot, pt.size, pch.use, cols.use, x1, x2, no.axes))
            
            multiplotList(pList, cols = nCol)
            rp()
          }
)

SingleFeaturePlot <- function(data.use, feature, data.plot, pt.size, pch.use, cols.use, x1, x2, no.axes){
  data.gene=na.omit(data.frame(data.use[feature,]))
  data.plot$gene = t(data.gene)
  brewer.gran <- 1
  if(length(cols.use) == 1){
    brewer.gran <- brewer.pal.info[cols.use,]$maxcolors
  }
  else{
    brewer.gran <- length(cols.use)
  }
  data.cut=as.numeric(as.factor(cut(as.numeric(data.gene),breaks = brewer.gran)))
  data.plot$col=as.factor(data.cut)
  p <- ggplot(data.plot, aes(x,y))
  if(brewer.gran != 2){
    if(length(cols.use) == 1){
      p <- p + geom_point(aes(color=col), size=pt.size, shape=pch.use) + theme(legend.position='none') + 
        scale_color_brewer(palette=cols.use)
    }
    else{
      p <- p + geom_point(aes(color=col), size=pt.size, shape=pch.use) + theme(legend.position='none') + 
        scale_color_manual(values=cols.use)
    }
  }
  else{
    p <- p + geom_point(aes(color=gene), size=pt.size, shape=pch.use) + theme(legend.position='none') +
      scale_color_gradientn(colors=cols.use) 
  }
  if(no.axes){
    p <- p + labs(title = feature, x ="", y="") +  theme(axis.line=element_blank(),axis.text.x=element_blank(),
                                                         axis.text.y=element_blank(),axis.ticks=element_blank(),
                                                         axis.title.x=element_blank(),
                                                         axis.title.y=element_blank())
  }
  else{
    p <- p + labs(title = feature, x = x1, y = x2)
  }
  return(p)
}

#' Vizualization of multiple features
#'
#' Similar to FeaturePlot, however, also splits the plot by visualizing each
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
setGeneric("FeatureHeatmap", function(object,features.plot,dim.1=1,dim.2=2,idents.use=NULL,pt.size=2,cols.use=rev(heat.colors(10)),pch.use=16,reduction.use="tsne") standardGeneric("FeatureHeatmap"))
#' @export
setMethod("FeatureHeatmap", "seurat",
          function(object,features.plot,dim.1=1,dim.2=2,idents.use=NULL,pt.size=2,cols.use=rev(heat.colors(10)),pch.use=16,reduction.use="tsne") {
            idents.use=set.ifnull(idents.use,sort(unique(object@ident)))
            dim.code="PC"
            par(mfrow=c(length(features.plot),length(idents.use)))
            dim.code=translate.dim.code(reduction.use); dim.codes=paste(dim.code,c(dim.1,dim.2),sep="")
            data.plot=data.frame(FetchData(object,dim.codes))

            ident.use=as.factor(object@ident)
            data.plot$ident=ident.use
            x1=paste(dim.code,dim.1,sep=""); x2=paste(dim.code,dim.2,sep="")
            data.plot$x=data.plot[,x1]; data.plot$y=data.plot[,x2]
            data.plot$pt.size=pt.size
            data.use=data.frame(t(data.frame(FetchData(object,features.plot))))
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
#' This function is a wrapper for DimPlot. See ?DimPlot for a full list of possible
#' arguments which can be passed in here.
#'
#' @param object Seurat object
#' @param do.label FALSE by default. If TRUE, plots an alternate view where the center of each
#' cluster is labeled
#' @param pt.size Set the point size
#' @param label.size Set the size of the text labels
#' @param colors.use Manually set the color palette to use for the points
#' @param \dots Additional parameters to DimPlot, for example, which dimensions to plot.
#' @seealso DimPlot
#' @export
setGeneric("TSNEPlot", function(object,do.label=FALSE, pt.size=1, label.size = 4, cells.use = NULL, colors.use = NULL, ...) standardGeneric("TSNEPlot"))
#' @export
setMethod("TSNEPlot", "seurat",
          function(object,do.label=FALSE, pt.size=1, label.size=4, cells.use = NULL, colors.use = NULL,...) {
            cells.use=set.ifnull(cells.use,object@cell.names)
            #print(head(cells.use))
            cells.use=ainb(cells.use,object@cell.names)
            return(DimPlot(object,reduction.use = "tsne",cells.use = cells.use, pt.size = pt.size, do.label = do.label, label.size = label.size, cols.use = colors.use, ...))
          }
)

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
setGeneric("ICAPlot", function(object,...) standardGeneric("ICAPlot"))
#' @export
setMethod("ICAPlot", "seurat",
          function(object,...) {
            return(DimPlot(object,reduction.use = "ica",...))
          }
)

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
setGeneric("PCAPlot", function(object,...) standardGeneric("PCAPlot"))
#' @export
setMethod("PCAPlot", "seurat",
          function(object,...) {
              return(DimPlot(object,reduction.use = "pca",...))
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
#' cell attribute (that can be pulled with FetchData) allowing for both different colors and different shapes on cells.
#' @param do.label Whether to label the clusters
#' @param label.size Sets size of labels
#' @param no.legend Setting to TRUE will remove the legend
#' @return If do.return==TRUE, returns a ggplot2 object. Otherwise, only
#' graphical output.
#' @importFrom dplyr summarize group_by
#' @export
setGeneric("DimPlot", function(object,reduction.use="pca",dim.1=1,dim.2=2,cells.use=NULL,pt.size=3,do.return=FALSE,do.bare=FALSE,cols.use=NULL,group.by="ident",pt.shape=NULL, do.label = FALSE, label.size = 1, no.legend = FALSE) standardGeneric("DimPlot"))
#' @export
setMethod("DimPlot", "seurat",
          function(object,reduction.use="pca",dim.1=1,dim.2=2,cells.use=NULL,pt.size=3,do.return=FALSE,do.bare=FALSE,cols.use=NULL,group.by="ident",pt.shape=NULL, do.label = FALSE, label.size = 1, no.legend = FALSE) {
            cells.use=set.ifnull(cells.use,colnames(object@data))
            dim.code=translate.dim.code(reduction.use); dim.codes=paste(dim.code,c(dim.1,dim.2),sep="")
            data.plot=FetchData(object,dim.codes,cells.use)

            ident.use=as.factor(object@ident[cells.use])
            if (group.by != "ident") ident.use=as.factor(FetchData(object,group.by)[,1])
            data.plot$ident=ident.use
            x1=paste(dim.code,dim.1,sep=""); x2=paste(dim.code,dim.2,sep="")
            data.plot$x=data.plot[,x1]; data.plot$y=data.plot[,x2]
            data.plot$pt.size=pt.size
            p=ggplot(data.plot,aes(x=x,y=y))+geom_point(aes(colour=factor(ident)),size=pt.size)
            if (!is.null(pt.shape)) {
              shape.val=FetchData(object,pt.shape)[cells.use,1]
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
            if (do.label) {
              data.plot %>% dplyr::group_by(ident) %>% summarize(x = median(x), y = median(y)) -> centers
              p3 <- p3 + geom_point(data = centers, aes(x=x, y=y), size=0, alpha=0) + geom_text(data=centers, aes(label=ident), size = label.size)
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
)

#Cool, but not supported right now
setGeneric("SpatialDe", function(object,marker.cells,genes.use=NULL,...) standardGeneric("SpatialDe"))
setMethod("SpatialDe", "seurat",
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

            diff.genes=rownames(subset(DiffExpTest(p15,marker.cells,embed.diff,genes.use=genes.use),p_val<(1e-5)))
            diff.genes=subset(DiffExpTest(p15,marker.cells,embed.diff,genes.use = diff.genes),p_val<(1e-10))
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
#' @param SetIdent TRUE by default. Set identity class to the results of the density clustering.
#' Unassigned cells (cells that cannot be assigned a cluster) are placed in cluster 1, if there are any.
#' @param seed.use Random seed for the dbscan function
#' @param \dots Additional arguments to be passed to the dbscan function
#' @export
setGeneric("DBClustDimension", function(object,dim.1=1,dim.2=2,reduction.use="tsne",G.use=NULL,set.ident=TRUE,seed.use=1,...) standardGeneric("DBClustDimension"))
#' @export
setMethod("DBClustDimension", "seurat",
          function(object,dim.1=1,dim.2=2,reduction.use="tsne",G.use=NULL,set.ident=TRUE,seed.use=1,...) {
            dim.code=translate.dim.code(reduction.use); dim.codes=paste(dim.code,c(dim.1,dim.2),sep="")
            data.plot=FetchData(object,dim.codes)
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
setGeneric("KClustDimension", function(object,dim.1=1,dim.2=2,cells.use=NULL,pt.size=4,reduction.use="tsne",k.use=5,set.ident=FALSE,seed.use=1,...) standardGeneric("KClustDimension"))
#' @export
setMethod("KClustDimension", "seurat",
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
#' Returns a set of genes, based on the JackStraw analysis, that have
#' statistically significant associations with a set of PCs.
#'
#'
#' @param object Seurat object
#' @param pcs.use PCS to use.
#' @param pval.cut P-value cutoff
#' @param use.full Use the full list of genes (from the projected PCA). Assumes
#' that ProjectPCA has been run. Default is TRUE
#' @param max.per.pc Maximum number of genes to return per PC. Used to avoid genes from one PC dominating the entire analysis.
#' @return A vector of genes whose p-values are statistically significant for
#' at least one of the given PCs.
#' @export
setGeneric("PCASigGenes", function(object,pcs.use,pval.cut=0.1,use.full=TRUE,max.per.pc=NULL) standardGeneric("PCASigGenes"))
#' @export
setMethod("PCASigGenes", "seurat",
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
              pc.top.genes=PCTopGenes(object,pcs.use,max.per.pc,use.full,FALSE)
              genes.use=ainb(pc.top.genes,genes.use)
            }
            return(genes.use)
          }
)

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
setGeneric("DoHeatmap", function(object,cells.use=NULL,genes.use=NULL,disp.min=-2.5,disp.max=2.5,draw.line=TRUE,do.return=FALSE,order.by.ident=TRUE,col.use=pyCols,slim.col.label=FALSE,group.by=NULL,remove.key=FALSE,cex.col=NULL,do.scale=TRUE,...) standardGeneric("DoHeatmap"))
#' @export
setMethod("DoHeatmap","seurat",
          function(object,cells.use=NULL,genes.use=NULL,disp.min=-2.5,disp.max=2.5,draw.line=TRUE,do.return=FALSE,order.by.ident=TRUE,col.use=pyCols,slim.col.label=FALSE,group.by=NULL,remove.key=FALSE,cex.col=NULL,do.scale=TRUE,...) {
            cells.use=set.ifnull(cells.use,object@cell.names)
            genes.use=ainb(genes.use,rownames(object@scale.data))
            cells.use=ainb(cells.use,object@cell.names)
            cells.ident=object@ident[cells.use]
            if (!is.null(group.by)) cells.ident=factor(FetchData(object,group.by)[,1])
            cells.ident=factor(cells.ident,labels = ainb(levels(cells.ident),cells.ident))
            if (order.by.ident) {
              cells.use=cells.use[order(cells.ident)]
            }
            data.use=object@scale.data[genes.use,cells.use]
            if (!do.scale) data.use=as.matrix(object@data[genes.use,cells.use])

            if (do.scale) data.use=minmax(data.use,min=disp.min,max=disp.max)
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
              cex.col=set.ifnull(cex.col,0.2+1/log10(length(unique(cells.ident))))
              hmFunction(data.use,Rowv=NA,Colv=NA,trace = "none",col=col.use,colsep = colsep.use,labCol=col.lab,cexCol=cex.col,...)
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
#' @inheritParams DoHeatmap
#' @inheritParams ICTopGenes
#' @inheritParams VizICA
#' @param use.scale Default is TRUE: plot scaled data. If FALSE, plot raw data on the heatmap.
#' @return If do.return==TRUE, a matrix of scaled values which would be passed
#' to heatmap.2. Otherwise, no return value, only a graphical output
#' @export
setGeneric("ICHeatmap", function(object,ic.use=1,cells.use=NULL,num.genes=30, disp.min=-2.5,disp.max=2.5,do.return=FALSE,col.use=pyCols,use.scale=TRUE,do.balanced=FALSE,...) standardGeneric("ICHeatmap"))
#' @export
setMethod("ICHeatmap","seurat",
          function(object,ic.use=1,cells.use=NULL,num.genes=30,disp.min=-2.5,disp.max=2.5,do.return=FALSE,col.use=pyCols,use.scale=TRUE,do.balanced=FALSE,...) {
            cells.use=set.ifnull(cells.use,object@cell.names)
            genes.use=ICTopGenes(object,ic.use,num.genes,do.balanced)
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
setGeneric("PCHeatmap", function(object,pc.use=1,cells.use=NULL,num.genes=30,use.full=FALSE, disp.min=-2.5,disp.max=2.5,do.return=FALSE,col.use=pyCols,use.scale=TRUE,do.balanced=FALSE,remove.key=FALSE, label.columns=NULL, ...) standardGeneric("PCHeatmap"))
#' @export
setMethod("PCHeatmap","seurat",
          function(object,pc.use=1,cells.use=NULL,num.genes=30,use.full=FALSE, disp.min=-2.5,disp.max=2.5,do.return=FALSE,col.use=pyCols,use.scale=TRUE,do.balanced=FALSE,remove.key=FALSE, label.columns=NULL, ...) {
            num.row=floor(length(pc.use)/3.01)+1
            orig_par <- par()$mfrow
            par(mfrow=c(num.row, min(length(pc.use),3)))
            cells <- cells.use
            plots <- c()

            if (is.null(label.columns)){
              if (length(pc.use) > 1){
                label.columns = FALSE
              }
              else{
                label.columns = TRUE
              }
            }

            for(pc in pc.use){
              if (is.numeric((cells))) {
                cells.use=PCTopCells(object,pc,cells,do.balanced)
              }
              else {
                cells.use=set.ifnull(cells,object@cell.names)
              }
              genes.use=rev(PCTopGenes(object,pc,num.genes,use.full,do.balanced))
              cells.ordered=cells.use[order(object@pca.rot[cells.use,pc])]
              data.use=object@scale.data[genes.use,cells.ordered]
              data.use=minmax(data.use,min=disp.min,max=disp.max)
              if (!(use.scale)) data.use=as.matrix(object@data[genes.use,cells.ordered])
              vline.use=NULL;

              if (remove.key || length(pc.use) > 1){
                hmFunction <- "heatmap2NoKey(data.use,Rowv=NA,Colv=NA,trace = \"none\",col=col.use, pc = pc, "
              }
              else{
                hmFunction <- "heatmap.2(data.use,Rowv=NA,Colv=NA,trace = \"none\",col=col.use, main = paste(\"PC\",pc) , "
              }

              if (!label.columns){

                hmFunction <- paste(hmFunction, "labCol=\"\", ", sep="")
              }
              hmFunction <- paste(hmFunction, "...)", sep="")
              eval(parse(text=hmFunction))
            }
            if (do.return) {
              return(data.use)
            }
            # reset graphics parameters
            par(mfrow=orig_par)
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
#' @param \dots Additional parameters passed to DoHeatmap for plotting
#' @return Seurat object where the k-means results for genes is stored in
#' object@@kmeans.obj[[1]], and the k-means results for cells is stored in
#' object@@kmeans.col[[1]]. The cluster for each cell is stored in object@@data.info[,"kmeans.ident"]
#' and also object@@ident (if set.ident=TRUE)
#' @export
setGeneric("DoKMeans", function(object,genes.use=NULL,k.genes=NULL,k.cells=NULL,k.seed=1,do.plot=TRUE,data.cut=2.5,k.cols=pyCols,
                                pc.row.order=NULL,pc.col.order=NULL, rev.pc.order=FALSE, use.imputed=FALSE,set.ident=TRUE,...) standardGeneric("DoKMeans"))
#' @export
setMethod("DoKMeans","seurat",
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
            if (k.cells>0) object@data.info[names(kmeans.col$cluster),"kmeans.ident"]=kmeans.col$cluster

            if ((set.ident) && (k.cells > 0)) {
              object=SetIdent(object,cells.use=names(kmeans.col$cluster),ident.use = kmeans.col$cluster)
            }
            if (do.plot) {
              disp.data=minmax(kmeans.data[order(kmeans.obj$cluster[genes.use]),],min=data.cut*(-1),max=data.cut)
              #DoHeatmap(object,object@cell.names,names(sort(kmeans.obj$cluster)),data.cut*(-1),data.cut,col.use = k.cols,...)
              KMeansHeatmap(object)
            }
            return(object)
          }
)

#' @export
setGeneric("GenesInCluster", function(object, cluster.num,max.genes=1e6)  standardGeneric("GenesInCluster"))
#' @export
setMethod("GenesInCluster", signature = "seurat",
          function(object, cluster.num,max.genes=1e6) {
            toReturn=(unlist(lapply(cluster.num,function(x)head(sort(names(which(object@kmeans.obj[[1]]$cluster==x))),max.genes))))
            return(toReturn)
          }
)

#Needs comments
#' @export
setGeneric("KMeansHeatmap", function(object, cells.use=object@cell.names,genes.cluster=NULL,max.genes=1e6,slim.col.label=TRUE,remove.key=TRUE,row.lines=T,...)  standardGeneric("KMeansHeatmap"))
#' @export
setMethod("KMeansHeatmap", signature = "seurat",
          function(object, cells.use=object@cell.names,genes.cluster=NULL,max.genes=1e6,slim.col.label=TRUE,remove.key=TRUE,row.lines=T,...) {
            genes.cluster=set.ifnull(genes.cluster,unique(object@kmeans.obj[[1]]$cluster))
            genes.use=GenesInCluster(object,genes.cluster,max.genes)
            cluster.lengths=sapply(genes.cluster,function(x)length(GenesInCluster(object,x)))
            print(cluster.lengths)
            rowsep.use=NA; if (row.lines) rowsep.use=cumsum(cluster.lengths)
            DoHeatmap(object,cells.use = cells.use,genes.use = genes.use,slim.col.label=slim.col.label,remove.key=remove.key,rowsep=rowsep.use,...)
          }
)

#' @export
setGeneric("CellCorMatrix", function(object, cor.genes=NULL,cell.inds=NULL, do.k=FALSE,k.seed=1,k.num=4,vis.low=(-1),vis.high=1,vis.one=0.8,pcs.use=1:3,col.use=pyCols)  standardGeneric("CellCorMatrix"))
#' @export
setMethod("CellCorMatrix", signature = "seurat",
          function(object, cor.genes=NULL,cell.inds=NULL, do.k=FALSE,k.seed=1,k.num=4,vis.low=(-1),vis.high=1,vis.one=0.8,pcs.use=1:3,col.use=pyCols) {
            cor.genes=set.ifnull(cor.genes,object@var.genes)
            cell.inds=set.ifnull(cell.inds,colnames(object@data))
            cor.genes=cor.genes[cor.genes%in%rownames(object@data)]
            data.cor=object@scale.data[cor.genes,cell.inds]
            cor.matrix=cor((data.cor))
            set.seed(k.seed); kmeans.cor=kmeans(cor.matrix,k.num)
            if (do.k) cor.matrix=cor.matrix[order(kmeans.cor$cluster),order(kmeans.cor$cluster)]
            kmeans.names=rownames(cor.matrix)
            row.annot=data.frame(cbind(kmeans.cor$cluster[kmeans.names],object@pca.rot[kmeans.names,pcs.use]))
            colnames(row.annot)=c("K",paste("PC",pcs.use,sep=""))
            cor.matrix[cor.matrix==1]=vis.one
            cor.matrix=minmax(cor.matrix,min = vis.low,max=vis.high)
            object@kmeans.cell=list(kmeans.cor)
            if (do.k) aheatmap(cor.matrix,col=col.use,Rowv=NA,Colv=NA,annRow=row.annot)
            if (!(do.k)) heatmap.2(cor.matrix,trace="none",Rowv=NA,Colv=NA,col=pyCols)
            return(object)
          }
)

#' @export
setGeneric("GeneCorMatrix", function(object, cor.genes=NULL,cell.inds=NULL, do.k=FALSE,k.seed=1,k.num=4,vis.low=(-1),vis.high=1,vis.one=0.8,pcs.use=1:3,col.use=pyCols)  standardGeneric("GeneCorMatrix"))
#' @export
setMethod("GeneCorMatrix", signature = "seurat",
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
setGeneric("CalinskiPlot", function(object,pcs.use,pval.cut=0.1,gene.max=15,col.max=25,use.full=TRUE)  standardGeneric("CalinskiPlot"))
#' @export
setMethod("CalinskiPlot","seurat",
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


#' Dot plot visualization
#'
#' Intuitive way of visualizing how gene expression changes across different identity classes (clusters).
#' The size of the dot encodes the percentage of cells within a class, while the color encodes the
#' AverageExpression level of 'expressing' cells (green is high).
#'
#' @param genes.plot Input vector of genes
#' @param cex.use Scaling factor for the dots (scales all dot sizes)
#' @param thresh.col The raw data value which corresponds to a red dot (lowest expression)
#' @param dot.min The fraction of cells at which to draw the smallest dot (default is 0.05)
#' @inheritParams VlnPlot
#' @return Only graphical output
#' @export
setGeneric("DotPlot", function(object,genes.plot,cex.use=2,cols.use=NULL,thresh.col=2.5,dot.min=0.05,group.by=NULL,...)  standardGeneric("DotPlot"))
#' @export
setMethod("DotPlot","seurat",
          function(object,genes.plot,cex.use=2,cols.use=NULL,thresh.col=2.5,dot.min=0.05,group.by=NULL,...) {
            if (!(is.null(group.by))) object=SetAllIdent(object,id = group.by)
            #object@data=object@data[genes.plot,]
            object@data=data.frame(t(FetchData(object,genes.plot)))
            avg.exp=AverageExpression(object)
            avg.alpha=ClusterAlpha(object)
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
#' anything that can be retreived by FetchData)
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
setGeneric("VlnPlot", function(object,features.plot,nCol=NULL,ylab.max=12,do.ret=TRUE,do.sort=FALSE,
                               size.x.use=16,size.y.use=16,size.title.use=20, use.imputed=FALSE,adjust.use=1,size.use=1,cols.use=NULL,group.by="ident")  standardGeneric("VlnPlot"))
#' @export
setMethod("VlnPlot","seurat",
          function(object,features.plot,nCol=NULL,ylab.max=12,do.ret=FALSE,do.sort=FALSE,size.x.use=16,size.y.use=16,size.title.use=20,use.imputed=FALSE,adjust.use=1,
                   size.use=1,cols.use=NULL,group.by=NULL) {
            if (is.null(nCol)) {
              nCol=2
              if (length(features.plot)>6) nCol=3
              if (length(features.plot)>9) nCol=4
            }
            
            if(length(features.plot == 1)) {
              data.use=data.frame(FetchData(object,features.plot,use.imputed=use.imputed))
              if(nrow(data.use) > 1) {
                data.use = t(data.use)
              }
              rownames(data.use) <- features.plot
            }
            
            else {
              data.use=data.frame(t(FetchData(object,features.plot,use.imputed=use.imputed)))
            }
            ident.use=object@ident
            if (!is.null(group.by)) ident.use=as.factor(FetchData(object,group.by)[,1])
            gene.names <- rownames(data.use)[rownames(data.use) %in% rownames(object@data)]
            pList=lapply(features.plot,function(x) plot.Vln(x,data.use[x,,drop=FALSE],ident.use,ylab.max,TRUE,do.sort,size.x.use,size.y.use,size.title.use,adjust.use,size.use,cols.use,gene.names))

            if(do.ret) {
              return(pList)
            }
            
            else {
              multiplotList(pList,cols = nCol)
              rp()
            }
          }
)

plot.Vln=function(gene,data,cell.ident,ylab.max=12,do.ret=FALSE,do.sort=FALSE,size.x.use=16,size.y.use=16,size.title.use=20,adjust.use=1,size.use=1,cols.use=NULL,gene.names) {
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
  p2=p + geom_violin(scale="width",adjust=adjust.use,trim=TRUE,aes(fill=factor(ident)))
  if(gene %in% gene.names){
    p2=p2 + ylab("Expression level (log TPM)")
  }
  else{
    p2 = p2 + ylab("")
  }
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



#' Add Metadata
#'
#' Adds additional data for single cells to the Seurat object. Can be any piece
#' of information associated with a cell (examples include read depth,
#' alignment rate, experimental batch, or subpopulation identity). The
#' advantage of adding it to the Seurat object is so that it can be
#' analyzed/visualized using FetchData, VlnPlot, GenePlot, SubsetData, etc.
#'
#'
#' @param object Seurat object
#' @param metadata Data frame where the row names are cell names (note : these
#' must correspond exactly to the items in object@@cell.names), and the columns
#' are additional metadata items.
#' @return Seurat object where the additional metadata has been added as
#' columns in object@@data.info
#' @export
setGeneric("AddMetaData", function(object,metadata)  standardGeneric("AddMetaData"))
#' @export
setMethod("AddMetaData","seurat",
          function(object,metadata) {
            cols.add=colnames(metadata)
            object@data.info[,cols.add]=metadata[rownames(object@data.info),]
            return(object)
          }
)

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
setGeneric("JackStrawPlot", function(object,PCs=1:5, nCol=3, score.thresh=1e-5,plot.x.lim=0.1,plot.y.lim=0.3)  standardGeneric("JackStrawPlot"))
#' @export
setMethod("JackStrawPlot","seurat",
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
#' be metrics, PC scores, etc. - anything that can be retreived with FetchData
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
setGeneric("GenePlot", function(object, gene1, gene2, cell.ids=NULL,col.use=NULL,
                                pch.use=16,cex.use=1.5,use.imputed=FALSE,do.ident=FALSE,do.spline=FALSE,spline.span=0.75,...)  standardGeneric("GenePlot"))
#' @export
setMethod("GenePlot","seurat",
          function(object, gene1, gene2, cell.ids=NULL,col.use=NULL,
                   pch.use=16,cex.use=1.5,use.imputed=FALSE,do.ident=FALSE,do.spline=FALSE,spline.span=0.75,...) {
            cell.ids=set.ifnull(cell.ids,object@cell.names)
            data.use=data.frame(t(FetchData(object,c(gene1,gene2),cells.use = cell.ids,use.imputed=use.imputed)))
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
              #lines(spline.fit$x,spline.fit$y,lwd=3)
              #spline.fit=smooth.spline(g1,g2,df = 4)
              loess.fit=loess(g2~g1,span=spline.span)
              #lines(spline.fit$x,spline.fit$y,lwd=3)
              points(g1,loess.fit$fitted,col="darkblue")
            }
            if (do.ident) {
              return(identify(g1,g2,labels = cell.ids))
            }
          }
)

#' @export
setGeneric("RemovePC", function(object, pcs.remove,use.full=FALSE,...)  standardGeneric("RemovePC"))
#' @export
setMethod("RemovePC","seurat",
          function(object, pcs.remove,use.full=FALSE,...) {
            data.old=object@data
            pcs.use=anotinb(1:ncol(object@pca.obj[[1]]$rotation),pcs.remove)
            data.x=as.matrix(object@pca.obj[[1]]$x[,pcs.use])
            if (use.full) data.x=as.matrix(object@pca.x.full[,ainb(pcs.use,1:ncol(object@pca.x.full))])
            data.1=data.x%*%t(as.matrix(object@pca.obj[[1]]$rotation[,pcs.use]))
            data.2=sweep(data.1,2,colMeans(object@scale.data),"+")
            data.3=sweep(data.2,MARGIN = 1,apply(object@data[rownames(data.2),],1,sd),"*")
            data.3=sweep(data.3,MARGIN = 1,apply(object@data[rownames(data.2),],1,mean),"+")
            object@scale.data=(data.2);
            data.old=data.old[rownames(data.3),]
            data.4=data.3; data.4[data.old==0]=0; data.4[data.4<0]=0
            object@data[rownames(data.4),]=data.frame(data.4)
            return(object)
          }
)

setGeneric("GeneScorePlot", function(object, gene1, score.name, cell.ids=NULL,col.use=NULL,nrpoints.use=Inf,pch.use=16,cex.use=2,...)  standardGeneric("GeneScorePlot"))
setMethod("GeneScorePlot","seurat",
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
setGeneric("CellPlot", function(object, cell1, cell2, gene.ids=NULL,col.use="black",nrpoints.use=Inf,pch.use=16,cex.use=0.5,do.ident=FALSE,...)  standardGeneric("CellPlot"))
#' @export
setMethod("CellPlot","seurat",
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
setGeneric("JackStrawPermutationTest", function(object,genes.use=NULL,num.iter=100, thresh.use=0.05,do.print=TRUE,k.seed=1)  standardGeneric("JackStrawPermutationTest"))
#' @export
setMethod("JackStrawPermutationTest","seurat",
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
setGeneric("JackStrawMC", function(object,num.pc=30,num.replicate=100,prop.freq=0.01,do.print=FALSE, num.cores=8)  standardGeneric("JackStrawMC"))
#' @export
setMethod("JackStrawMC","seurat",
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
#' @param rev.pca By default computes the PCA on the cell x gene matrix. Setting to 
#' true will compute it on gene x cell matrix. This should match what was set when the intial PCA was run.
#' @return Returns a Seurat object where object@@jackStraw.empP represents
#' p-values for each gene in the PCA analysis. If ProjectPCA is subsequently
#' run, object@@jackStraw.empP.full then represents p-values for all genes.
#' @importFrom pbapply pbsapply
#' @references Inspired by Chung et al, Bioinformatics (2014)
#' @export
setGeneric("JackStraw", function(object,num.pc=30,num.replicate=100,prop.freq=0.01,do.print=FALSE, rev.pca=FALSE)  standardGeneric("JackStraw"))
#' @export
setMethod("JackStraw","seurat",
          function(object,num.pc=30,num.replicate=100,prop.freq=0.01,do.print=FALSE, rev.pca=FALSE) {

            # error checking for number of PCs
            if (num.pc > ncol(object@pca.rot)){
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
            
            if (!(do.print)) fake.pcVals.raw=sapply(1:num.replicate,function(x)jackRandom(scaled.data=object@scale.data[pc.genes,],prop=prop.freq,r1.use = 1,r2.use = num.pc,seed.use=x, rev.pca=rev.pca),simplify = FALSE)
            if ((do.print)) fake.pcVals.raw=pbsapply(1:num.replicate,function(x){jackRandom(scaled.data=object@scale.data[pc.genes,],prop=prop.freq,r1.use = 1,r2.use = num.pc,seed.use=x,rev.pca=rev.pca)},simplify = FALSE)
            
            fake.pcVals=sapply(1:num.pc,function(x)as.numeric(unlist(lapply(1:num.replicate,function(y)fake.pcVals.raw[[y]][,x]))))
            object@jackStraw.fakePC = data.frame(fake.pcVals)
            object@jackStraw.empP=data.frame(sapply(1:num.pc,function(x)unlist(lapply(abs(md.x[,x]),empP,abs(fake.pcVals[,x])))))
            colnames(object@jackStraw.empP)=paste("PC",1:ncol(object@jackStraw.empP),sep="")
            return(object)
          }
)

#' @export
jackRandom=function(scaled.data,prop.use=0.01,r1.use=1,r2.use=5, seed.use=1,rev.pca=FALSE) {
  set.seed(seed.use)
  rand.genes <- sample(rownames(scaled.data), nrow(scaled.data) * prop.use)

  # make sure that rand.genes is at least 3
  if (length(rand.genes) < 3){
    rand.genes <- sample(rownames(scaled.data), 3)
  }

  data.mod <- scaled.data
  data.mod[rand.genes, ] <- shuffleMatRow(scaled.data[rand.genes, ])
  
  if(rev.pca){
    fake.pca <- prcomp(data.mod)
    fake.x <- fake.pca$x
    fake.rot <- fake.pca$rotation
  }
  else {
    fake.pca <- prcomp(t(data.mod))
    fake.x <- fake.pca$rotation
    fake.rot <- fake.pca$x
  }

  return(fake.x[rand.genes, r1.use:r2.use])
}

setGeneric("JackStrawFull", function(object,num.pc=5,num.replicate=100,prop.freq=0.01)  standardGeneric("JackStrawFull"))
setMethod("JackStrawFull","seurat",
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

#' Quickly Pick Relevant PCs
#'
#' Plots the standard deviations of the principle components for easy
#' identification of an elbow in the graph. This often corresponds well with the
#' significant PCs.
#'
#'
#' @param object Seurat object
#' @param num.pc Number of PCs to plot
#' @return Returns ggplot object
#' @export
setGeneric("PCElbowPlot", function(object,num.pc = 20)  standardGeneric("PCElbowPlot"))
#' @export
setMethod("PCElbowPlot","seurat",
          function(object,num.pc = 20) {
            if (length(object@pca.obj) == 0) {
              stop("This object has no PCA associated with it. Please run PCA() and then retry.")
            }
            if (length(object@pca.obj[[1]]$sdev) < num.pc) {
              num.pc <- length(object@pca.obj[[1]]$sdev)
              warning(paste("The object only has information for", num.pc, "PCs." ))
            }
            sdev <- object@pca.obj[[1]]$sdev[1:num.pc]
            pc <- 1:length(sdev)
            data <- data.frame(pc, sdev)
            plot <- ggplot(data, aes(pc, sdev)) + geom_point() + labs(y = "Standard Deviation of PC", x = "PC")
            return (plot)
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
#' @param do.recalc TRUE by default. If FALSE, plots and selects variable genes without recalculating statistics for each gene.
#' @importFrom MASS kde2d
#' @return Returns a Seurat object, placing variable genes in object@@var.genes.
#' The result of all analysis is stored in object@@mean.var
#' @export
setGeneric("MeanVarPlot", function(object, fxn.x=expMean, fxn.y=logVarDivMean,do.plot=TRUE,set.var.genes=TRUE,do.text=TRUE,
                                     x.low.cutoff=4,x.high.cutoff=8,y.cutoff=2,y.high.cutoff=Inf,cex.use=0.5,cex.text.use=0.5,do.spike=FALSE,
                                     pch.use=16, col.use="black", spike.col.use="red",plot.both=FALSE,do.contour=TRUE,
                                     contour.lwd=3, contour.col="white", contour.lty=2,num.bin=20,do.recalc=TRUE) standardGeneric("MeanVarPlot"))
#' @export
setMethod("MeanVarPlot", signature = "seurat",
          function(object, fxn.x=expMean, fxn.y=logVarDivMean, do.plot=TRUE,set.var.genes=TRUE,do.text=TRUE,
                   x.low.cutoff=4,x.high.cutoff=8,y.cutoff=1,y.high.cutoff=Inf,cex.use=0.5,cex.text.use=0.5,do.spike=FALSE,
                   pch.use=16, col.use="black", spike.col.use="red",plot.both=FALSE,do.contour=TRUE,
                   contour.lwd=3, contour.col="white", contour.lty=2,num.bin=20,do.recalc=TRUE) {
            data=object@data
            
            if (do.recalc) {    
                genes.use <- rownames(object@data)
                data.x=rep(0,length(genes.use)); names(data.x)=genes.use; data.y=data.x; data.norm.y=data.x;
    
                bin.size <- 1000
                max.bin <- floor(length(genes.use)/bin.size) + 1
                pb <- txtProgressBar(min = 0, max = max.bin, style = 3)
                for(i in 1:max.bin) {
                  my.inds <- ((bin.size * (i - 1)):(bin.size * i - 1))+1
                  my.inds <- my.inds[my.inds <= length(genes.use)]
                  genes.iter=genes.use[my.inds]; data.iter=data[genes.iter,]
                  data.x[genes.iter]=apply(data.iter,1,fxn.x); data.y[genes.iter]=apply(data.iter,1,fxn.y)
                  setTxtProgressBar(pb, i)  
                }
                close(pb)
                data.y[is.na(data.y)]=0
                data.x[is.na(data.x)]=0
                
                data_x_bin=cut(data.x,num.bin)
                names(data_x_bin)=names(data.x)
                mean_y=tapply(data.y,data_x_bin,mean)
                sd_y=tapply(data.y,data_x_bin,sd)
                data.norm.y=(data.y-mean_y[as.numeric(data_x_bin)])/sd_y[as.numeric(data_x_bin)]
                #data.x=apply(data,1,fxn.x); data.y=apply(data,1,fxn.y); data.x[is.na(data.x)]=0
                #data.norm.y=meanNormFunction(data,fxn.x,fxn.y,num.bin)
                
                data.norm.y[is.na(data.norm.y)]=0
                names(data.norm.y)=names(data.x)
                
                mv.df=data.frame(data.x,data.y,data.norm.y)
                rownames(mv.df)=rownames(data)
                object@mean.var=mv.df
            }
            data.x=object@mean.var[,1]; data.y=object@mean.var[,2]; data.norm.y=object@mean.var[,3]; 
            names(data.x)=names(data.y)=names(data.norm.y)=rownames(object@data)
            
            pass.cutoff=names(data.x)[which(((data.x>x.low.cutoff) & (data.x<x.high.cutoff)) & (data.norm.y>y.cutoff) & (data.norm.y < y.high.cutoff))]
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
