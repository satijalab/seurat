#install.packages(c("ROCR","ggplot2","Hmisc","reshape","gplots","stringr","NMF","mixtools","lars","XLConnect","reshape2","vioplot","fastICA","tsne","Rtsne","fpc","ape"))

seurat <- setClass("seurat", slots = 
                     c(raw.data = "data.frame", data="data.frame",scale.data="matrix",var.genes="vector",is.expr="numeric",
                       stat.fxn="function",data.stat="vector",data.ngene="vector",pca.x="data.frame",pca.rot="data.frame",
                       real.fval="data.frame",fake.fval="data.frame",emp.pval="data.frame",kmeans.obj="list",pca.obj="list",
                       calinski.best="numeric", gene.scores="data.frame", k.num = "numeric", drop.coefs="data.frame",
                       wt.matrix="data.frame", drop.wt.matrix="data.frame",trusted.genes="vector",drop.expr="numeric",data.info="data.frame",
                       project.name="character",project.dir="character", kmeans.gene="list", kmeans.cell="list",jackStraw.empP="data.frame", pc.x.full="data.frame",
                       jackStraw.fakePC = "data.frame",jackStraw.empP.full="data.frame",pca.x.full="data.frame", kmeans.col="list",mean.var="data.frame", imputed="data.frame",mix.probs="data.frame",
                       mix.mu="data.frame",mix.sigma="data.frame",mu.alpha="data.frame",mix.param="data.frame",final.prob="data.frame",insitu.matrix="data.frame",
                       tsne.rot="data.frame", ica.rot="data.frame", ica.x="data.frame",ica.obj="list"))

calc.drop.prob=function(x,a,b) {
  return(exp(a+b*x)/(1+exp(a+b*x)))
}

setGeneric("find_all_markers", function(object, thresh.test=1,test.use="bimod",by.k=FALSE,return.thresh=1e-2,do.print=FALSE) standardGeneric("find_all_markers"))
setMethod("find_all_markers","seurat",
    function(object, thresh.test=1,test.use="bimod",by.k=FALSE,return.thresh=1e-2,do.print=FALSE) {
      stat.use=object@data.stat
      if ((test.use=="roc") && (return.thresh==1e-2)) return.thresh=0.8
      if (by.k) stat.use=as.factor(retreiveCluster(object))
      stats.all=sort(unique(object@data.stat))
      genes.de=list()
      for(i in 1:length(stats.all)) {
        genes.de[[i]]=find.markers(object,stats.all[i],genes.use=rownames(object@data),thresh.use = thresh.test,test.use = test.use)
        if (do.print) print(paste("Calculating cluster", stats.all[i]))
      }
      gde.all=data.frame()
      for(i in 1:length(stats.all)) {
        gde=genes.de[[i]]
        if (nrow(gde)>0) {
          if (test.use=="roc") gde=subset(gde,(myAUC>return.thresh|myAUC<(1-return.thresh)))
          if (test.use=="bimod") {
            gde=gde[order(gde$myP,-gde$myDiff),]
            gde=subset(gde,myP<return.thresh)
          }
          if (nrow(gde)>0) gde$cluster=stats.all[i]; gde$gene=rownames(gde)
          if (nrow(gde)>0) gde.all=rbind(gde.all,gde)
        }
      }
      return(gde.all)
    }
)
    
    
setGeneric("plotNoiseModel", function(object, cell.ids=c(1,2), col.use="black",lwd.use=2,do.new=TRUE,x.lim=10,...) standardGeneric("plotNoiseModel"))
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

setGeneric("setup", function(object, project, min.cells=3, min.genes=2500, is.expr=1, do.scale=TRUE, do.center=TRUE,proj.dir=NULL,calc.noise=FALSE,...) standardGeneric("setup"))
setMethod("setup","seurat",
          function(object, project, min.cells=3, min.genes=2500, is.expr=1, do.scale=TRUE, do.center=TRUE,proj.dir=NULL,calc.noise=TRUE,...) {
            object@is.expr = is.expr
            num.genes=findNGene(object@raw.data,object@is.expr)
            cells.use=names(num.genes[which(num.genes>min.genes)]) 
            
            object@data=object@raw.data[,cells.use]
            num.cells=apply(object@data,1,humpCt,min=object@is.expr)
            genes.use=names(num.cells[which(num.cells>min.cells)])
            object@data=object@data[genes.use,]
            
            object@data.stat=unlist(lapply(colnames(object@data),object@stat.fxn))
            names(object@data.stat)=colnames(object@data)
            object@scale.data=t(scale(t(object@data),center=do.center,scale=do.scale))
            object@data.ngene=num.genes[cells.use]
            object@gene.scores=data.frame(object@data.ngene); colnames(object@gene.scores)[1]="nGene"
            object@data.info=data.frame(object@data.ngene); colnames(object@data.info)[1]="nGene"
            object@mix.probs=data.frame(object@data.ngene); colnames(object@mix.probs)[1]="nGene"
            rownames(object@gene.scores)=colnames(object@data)
            
            object@data.info[names(object@data.stat),"stat"]=object@data.stat
            
            object@project.name=project
            object@project.dir=set.ifnull(proj.dir,paste("~/big/",project,"/",sep=""))
            if(calc.noise) {
              object=calcNoiseModels(object,...)
              object=getWeightMatrix(object)
            }
            return(object)
          }         
)

setGeneric("subsetData",  function(object, subset.name=NULL, cells.use=NULL,accept.low=0, accept.high=Inf,do.center=TRUE,do.scale=TRUE) standardGeneric("subsetData"))
setMethod("subsetData","seurat",
          function(object, subset.name=NULL, cells.use=NULL,accept.low=0, accept.high=Inf,do.center=TRUE,do.scale=TRUE) {
            data.use=NULL
            if (is.null(cells.use)) {
              if (subset.name%in%rownames(object@data)) data.use=object@data
              if (subset.name%in%colnames(object@gene.scores)) data.use=data.frame(t(object@gene.scores))
              if (subset.name%in%colnames(object@data.info)) data.use=data.frame(t(object@data.info))
              if (is.null(data.use)) {
                print(paste("ERROR ", subset.name, " not found", sep=""))
                return(object)
              }
              subset.data=data.use[subset.name,]
              pass.inds=which((subset.data>accept.low) & (subset.data<accept.high))
              cells.use=colnames(data.use)[pass.inds]
            }
            object@data=object@data[,cells.use]
            object@scale.data=t(scale(t(object@data),center=do.center,scale=do.scale))
            object@scale.data=object@scale.data[complete.cases(object@scale.data),]
            object@data.ngene=object@data.ngene[cells.use]
            object@data.stat=object@data.stat[cells.use]
            object@tsne.rot=object@tsne.rot[cells.use,]
            object@pca.rot=object@pca.rot[cells.use,]
            
            object@gene.scores=data.frame(object@gene.scores[cells.use,]); colnames(object@gene.scores)[1]="nGene"; rownames(object@gene.scores)=colnames(object@data)
            object@data.info=data.frame(object@data.info[cells.use,])
            object@mix.probs=data.frame(object@mix.probs[cells.use,]); colnames(object@mix.probs)[1]="nGene"; rownames(object@mix.probs)=colnames(object@data)
            
            return(object)
          }         
)

setGeneric("loadMetrics", function(object, metrics.file=NULL, col.names=NULL,sep.use="\t",row.add="_rsem",...) standardGeneric("loadMetrics"))
setMethod("loadMetrics","seurat",
          function(object, metrics.file=NULL, col.names=NULL,sep.use="\t",row.add="_rsem",...) {
            metrics.file=set.ifnull(metrics.file,paste(object@project.dir,"summary/",object@project.name,".all.aln.metrics.txt",sep=""))
            col.names=set.ifnull(col.names,c("n.read","n.aln.read","pct.aln.rsem","pct.aln.ribo","pct.aln.spike","pct.aln.coli"))
            metrics.data=read.table(metrics.file,sep=sep.use,row.names="V1",...)
            colnames(metrics.data)=col.names
            rownames(metrics.data)=paste(sub.string(rownames(metrics.data),"-","_"),row.add,sep="")
            object@data.info=cbind(object@data.info,metrics.data[rownames(object@data.info),])
            return(object)
          }
)

setGeneric("project.pca", function(object,pc.genes=NULL,do.print=TRUE,pcs.print=5,pcs.store=20,genes.print=30,replace.pc=FALSE) standardGeneric("project.pca"))
setMethod("project.pca", "seurat", 
          function(object,pc.genes=NULL,do.print=TRUE,pcs.print=5,pcs.store=20,genes.print=30,replace.pc=FALSE) {
            object@pca.x.full=data.frame(as.matrix(object@scale.data)%*%as.matrix(object@pca.rot))
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
              pc_scores=object@pca.x.full
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

setGeneric("run_tsne", function(object,cells.use=NULL,pcs.use=1:10,do.plot=TRUE,by.k=FALSE,col.use=NULL,k.seed=1,do.fast=FALSE,...) standardGeneric("run_tsne"))
setMethod("run_tsne", "seurat", 
          function(object,cells.use=NULL,pcs.use=1:10,do.plot=TRUE,by.k=FALSE,col.use=NULL,k.seed=1,do.fast=FALSE,...) {
            cells.use=set.ifnull(cells.use,colnames(object@data))
            data.use=object@pca.rot[cells.use,pcs.use]
            if (do.fast) {
              set.seed(k.seed); data.tsne=Rtsne(as.matrix(data.use),...)
              data.tsne=data.frame(data.tsne$Y)
            }
            if (!(do.fast)) {
              set.seed(k.seed); data.tsne=data.frame(tsne(data.use,...))
            }
            print(head(data.tsne))
            colnames(data.tsne)=paste("TSNE_",1:ncol(data.tsne),sep="")
            rownames(data.tsne)=cells.use
            object@tsne.rot=data.tsne
            return(object)
          }
)

setGeneric("ica", function(object,ic.genes=NULL,do.print=TRUE,ics.print=5,ics.store=10,genes.print=30,use.imputed=FALSE) standardGeneric("ica"))
setMethod("ica", "seurat", 
          function(object,ic.genes=NULL,do.print=TRUE,ics.print=5,ics.store=10,genes.print=30,use.imputed=FALSE) {
            data.use=object@scale.data
            if (use.imputed) data.use=data.frame(t(scale(t(object@imputed))))
            ic.genes=set.ifnull(ic.genes,object@var.genes)
            ic.genes = ic.genes[ic.genes%in%rownames(data.use)]
            ic.genes.var = apply(data.use[ic.genes,],1,var)
            ic.data = data.use[ic.genes[ic.genes.var>0],]
            ica.obj = (fastICA(t(ic.data),n.comp=ics.store))
            object@ica.obj=list(ica.obj)
            ics.store=min(ics.store,ncol(ic.data))
            ics.print=min(ics.print,ncol(ic.data))
            ic_scores=data.frame(ic.data%*%ica.obj$S)
            colnames(ic_scores)=paste("IC",1:ncol(ic_scores),sep="")
            object@ica.x=ic_scores
            object@ica.rot=data.frame(ica.obj$S[,1:ics.store])
            colnames(object@ica.rot)=paste("IC",1:ncol(object@ica.rot),sep="")
            print(head(object@ica.rot))
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

setGeneric("pca", function(object,pc.genes=NULL,do.print=TRUE,pcs.print=5,pcs.store=40,genes.print=30,use.imputed=FALSE) standardGeneric("pca"))
setMethod("pca", "seurat", 
          function(object,pc.genes=NULL,do.print=TRUE,pcs.print=5,pcs.store=40,genes.print=30,use.imputed=FALSE) {
            data.use=object@scale.data
            if (use.imputed) data.use=data.frame(t(scale(t(object@imputed))))
            pc.genes=set.ifnull(pc.genes,object@var.genes)
            pc.genes = pc.genes[pc.genes%in%rownames(data.use)]
            pc.genes.var = apply(data.use[pc.genes,],1,var)
            pc.data = data.use[pc.genes[pc.genes.var>0],]
            pca.obj = (prcomp(pc.data))
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

setGeneric("cluster.alpha", function(object,by.k=FALSE,thresh.min=0) standardGeneric("cluster.alpha"))
setMethod("cluster.alpha", "seurat", 
          function(object,by.k=FALSE,thresh.min=0) {
            stat.use=object@data.stat
            if (by.k) stat.use=retreiveCluster(object)
            data.all=data.frame(row.names = rownames(object@data))
            for(i in sort(unique(stat.use))) {
              temp.cells=names(which(stat.use==i))
              data.temp=apply(object@data[,temp.cells],1,function(x)return(length(x[x>thresh.min])/length(x)))
              data.all=cbind(data.all,data.temp)
              colnames(data.all)[ncol(data.all)]=i
            }
            colnames(data.all)=sort(unique(stat.use))
            return(data.all)
          }
)



setGeneric("average.expression", function(object,by.k=FALSE) standardGeneric("average.expression"))
setMethod("average.expression", "seurat", 
          function(object,by.k=FALSE) {
            stat.use=object@data.stat
            if (by.k) stat.use=retreiveCluster(object)
            data.all=data.frame(row.names = rownames(object@data))
            for(i in sort(unique(stat.use))) {
              temp.cells=names(which(stat.use==i))
              data.temp=apply(object@data[,temp.cells],1,expMean)
              data.all=cbind(data.all,data.temp)
              colnames(data.all)[ncol(data.all)]=i
            }
            colnames(data.all)=sort(unique(stat.use))
            return(data.all)
          }
)

setGeneric("print.pca", function(object,pcs.print=1:5,genes.print=30,use.full=FALSE) standardGeneric("print.pca"))
setMethod("print.pca", "seurat", 
          function(object,pcs.print=1:5,genes.print=30,use.full=FALSE) {
            pc_scores=object@pca.x
            if (use.full==TRUE) pc_scores = object@pca.x.full
            for(i in pcs.print) {
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
)


set.ifnull=function(x,y) {
  if(is.null(x)) x=y
  return(x)
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
            my.stat=object@data.stat[cells.use]
            reg.score=unlist(lapply(candidate.reg,regressionSig,score = my.score,data = my.data,latent = my.stat,code = "rsem"))
            names(reg.score)=candidate.reg
            return(reg.score)
          }
)

setGeneric("find.markers.node", function(object,tree,node,genes.use=NULL,thresh.use=log(2),by.k=FALSE,test.use="bimod") standardGeneric("find.markers.node"))
setMethod("find.markers.node", "seurat",
          function(object,tree,node,genes.use=NULL,thresh.use=log(2),by.k=FALSE,test.use="bimod") {
            nodes.1=getLeftDecendants(tree,node)
            nodes.2=getRightDecendants(tree,node)
            to.return=find.markers(object,nodes.1,nodes.2,genes.use,thresh.use,by.k,test.use)
            return(to.return)
          } 
)

setGeneric("find.markers", function(object, stat.1,stat.2=NULL,genes.use=NULL,thresh.use=log(2),by.k=FALSE,test.use="bimod") standardGeneric("find.markers"))
setMethod("find.markers", "seurat",
          function(object, stat.1,stat.2=NULL,genes.use=NULL,thresh.use=log(2),by.k=FALSE,test.use="bimod") {
            genes.use=set.ifnull(genes.use,rownames(object@data))
            stat.use=object@data.stat
            if (by.k) stat.use=retreiveCluster(object)
            cells.1=names(stat.use[which(stat.use%in%stat.1)])
            if (is.null(stat.2)) {
              cells.2=names(stat.use)
              cells.2=cells.2[!(cells.2%in%cells.1)]
            }
            else {
              cells.2=names(stat.use[which(stat.use%in%stat.2)])
            }
            if (test.use=="bimod") to.return=diffExp.test(object,cells.1,cells.2,genes.use,thresh.use) 
            if (test.use=="roc") to.return=marker.test(object,cells.1,cells.2,genes.use,thresh.use) 
            return(to.return)
          } 
)

setGeneric("diffExp.test", function(object, cells.1,cells.2,genes.use=NULL,thresh.use=log(2)) standardGeneric("diffExp.test"))
setMethod("diffExp.test", "seurat",
          function(object, cells.1,cells.2,genes.use=NULL,thresh.use=log(2)) {
            genes.use=set.ifnull(genes.use,object@var.genes)
            data.1=apply(object@data[genes.use,cells.1],1,expMean)
            data.2=apply(object@data[genes.use,cells.2],1,expMean)
            total.diff=abs(data.1-data.2)
            genes.diff = names(which(total.diff>thresh.use))
            #print(genes.diff)
            to.return=bimod.diffExp.test(object@data[,cells.1],object@data[,cells.2],genes.diff)
            to.return=to.return[order(to.return$myP,-to.return$myDiff),]
            return(to.return)
          } 
)

setGeneric("batch.gene", function(object, stats.use,genes.use=NULL,use.imputed=FALSE,auc.cutoff=0.6) standardGeneric("batch.gene"))
setMethod("batch.gene", "seurat",
          function(object, stats.use,genes.use=NULL,use.imputed=FALSE,auc.cutoff=0.6) {
            batch.genes=c()
            genes.use=set.ifnull(genes.use,rownames(object@data))
            for(stat in stats.use ) {
              cells.1=names(object@data.stat)[object@data.stat==stat]
              cells.2=names(object@data.stat)[object@data.stat!=stat]
              if ((length(cells.1)<5)|(length(cells.2)<5)) {
                break;
              }
              markers.stat=marker.test(object,cells.1,cells.2,genes.use)
              batch.genes=unique(c(batch.genes,rownames(subset(markers.stat,myAUC>auc.cutoff))))
            }
            return(batch.genes)
          } 
)

setGeneric("marker.test", function(object, cells.1,cells.2,genes.use=NULL,thresh.use=log(2)) standardGeneric("marker.test"))
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
            return(to.return)
          } 
)

setGeneric("retreiveScore", function(object, score.name) standardGeneric("retreiveScore"))
setMethod("retreiveScore", "seurat",
          function(object, score.name) {
            my.score=object@gene.scores[,score.name]
            names(my.score)=rownames(object@gene.scores)
            return(my.score)
          } 
)

setGeneric("which.cells", function(object,value=1, id=NULL) standardGeneric("which.cells"))
setMethod("which.cells", "seurat",
          function(object, value=1,id=NULL) {
            id=set.ifnull(id,"data.stat")
            if (id=="data.stat") {
              data.use=object@data.stat
            }
            if (id %in% colnames(object@data.info)) {
              data.use=object@data.info[,id]; names(data.use)=rownames(object@data.info)
            }
            return(names(which(data.use==value)))
          } 
)

setGeneric("set.all.stat", function(object,id=NULL) standardGeneric("set.all.stat"))
setMethod("set.all.stat", "seurat",
          function(object, id=NULL) {
            id=set.ifnull(id,"stat")
            if (id %in% colnames(object@data.info)) {
              object@data.stat[rownames(object@data.info)]=object@data.info[,id]
            }
            return(object)
          } 
)

setGeneric("set.stat", function(object,cells.use=NULL,stat.use=NULL) standardGeneric("set.stat"))
setMethod("set.stat", "seurat",
          function(object, cells.use=NULL,stat.use=NULL) {
            id=set.ifnull(id,"stat")
            if (id %in% colnames(object@data.info)) {
              object@data.stat[rownames(object@data.info)]=object@data.info[,id]
            }
            return(object)
          } 
)


setGeneric("retreiveCellInfo", function(object, cells.use=NULL,id=NULL) standardGeneric("retreiveCellInfo"))
setMethod("retreiveCellInfo", "seurat",
          function(object, cells.use=NULL,id=NULL) {
            cells.use=set.ifnull(cells.use,colnames(object@data))
            id=set.ifnull(id,"data.stat")
            if (id=="data.stat") {
              data.use=object@data.stat
            }
            if (id %in% colnames(object@data.info)) {
              data.use=object@data.info[,id]; names(data.use)=rownames(object@data.info)
            }
            return(data.use[cells.use])
          } 
)


setGeneric("retreiveCluster", function(object, names=NULL) standardGeneric("retreiveCluster"))
setMethod("retreiveCluster", "seurat",
          function(object, names=NULL) {
            names=set.ifnull(names,colnames(object@data))
            if (names[1] %in% rownames(object@data)) return(object@kmeans.obj[[1]]$cluster[names])
            if (names[1] %in% colnames(object@data)) return(object@kmeans.col[[1]]$cluster[names])
          } 
)

setGeneric("posterior.plot", function(object, name) standardGeneric("posterior.plot"))
setMethod("posterior.plot", "seurat",
          function(object, name) {
            post.names=colnames(subc(object@mix.probs,name))
            vlnPlot(object,post.names,inc.first=TRUE,inc.final=TRUE,by.k=TRUE)
            
            
          } 
)

map.cell.score=function(gene,gene.value,insitu.bin,mu,sigma,alpha) {
  code.1=paste(gene,insitu.bin,sep=".")
  mu.use=mu[paste(code.1,"mu",sep="."),1]
  sigma.use=sigma[paste(code.1,"sigma",sep="."),1]
  alpha.use=alpha[paste(code.1,"alpha",sep="."),1]
  bin.prob=unlist(lapply(1:length(insitu.bin),function(x) dnorm(gene.value,mean = mu.use[x],sd = sigma.use[x],log = TRUE) + log(alpha.use[x])))
  return(bin.prob)
}

setGeneric("map.cell",  function(object,cell.name,do.plot=FALSE,safe.use=TRUE,text.val=NULL,do.rev=FALSE) standardGeneric("map.cell"))
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
            if (length(missing.cols)>0) print(paste("Error : ", all.needed.cols[missing.cols], " is missing from the mixture fits",sep=""))
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


setGeneric("refined.mapping",  function(object,genes.use) standardGeneric("refined.mapping"))
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


setGeneric("initial.mapping", function(object,cells.use=NULL,do.plot=FALSE,safe=FALSE) standardGeneric("initial.mapping"))
setMethod("initial.mapping", "seurat",
          function(object,cells.use=NULL,do.plot=FALSE,safe=FALSE) {
            cells.use=set.ifnull(cells.use,colnames(object@data))
            every.prob=sapply(cells.use,function(x)map.cell(object,x,do.plot=FALSE,safe.use=safe))
            object@final.prob=data.frame(every.prob)
            rownames(object@final.prob)=paste("bin.",rownames(object@final.prob),sep="")
            return(object)
          } 
)

bwCols=myPalette(low = "white",high="black",k = 50)

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

setGeneric("fit.gene.k", function(object, gene, do.k=3,num.iter=2,do.plot=FALSE,genes.use=NULL,start.pct=NULL) standardGeneric("fit.gene.k"))
setMethod("fit.gene.k", "seurat",
          function(object, gene, do.k=3,num.iter=2,do.plot=FALSE,genes.use=NULL,start.pct=NULL) {
            data=object@imputed            
            data.use=data[gene,]
            names(data.use)=colnames(data.use)
            scale.data=t(scale(t(object@imputed)))
            genes.use=set.ifnull(genes.use,rownames(scale.data))
            genes.use=genes.use[genes.use%in%rownames(scale.data)]
            scale.data=scale.data[genes.use,]
            #print(genes.use)
            #seed the k-means based on the 0'th component
            data.cut=as.numeric(data.use[gene,])
            cell.ident=as.numeric(cut(data.cut,do.k))
            if (!(is.null(start.pct))) {
              cell.ident=rep(1,length(data.cut))
              cell.ident[data.cut>quantile(data.cut,1-start.pct)]=2
            }
            cell.ident=order(tapply(as.numeric(data.use),cell.ident,mean))[cell.ident]
            #cell.ident=sample(cell.ident)
            ident.table=table(cell.ident)
#            if (do.plot) {
#              par(mfrow=c(2,2))
#              hist(as.numeric(data.use),probability = TRUE,ylim=c(0,1),xlab=gene,main=gene);
#              for(i in 1:do.k) lines(seq(-10,10,0.01),(ident.table[i]/sum(ident.table)) * dnorm(seq(-10,10,0.01),mean(as.numeric(data.use[cell.ident==i])),sd(as.numeric(data.use[cell.ident==i]))),col=i,lwd=2); 
#            }
            if (num.iter > 0) {
              for(i2 in 1:num.iter) {
                cell.ident=iter.k.fit(scale.data,cell.ident,data.use)
                ident.table=table(cell.ident)
#                if (do.plot) {
#                  hist(as.numeric(data.use),probability = TRUE,ylim=c(0,1),xlab=gene,main=gene);
#                  for(i in 1:do.k) lines(seq(-10,10,0.01),(ident.table[i]/sum(ident.table)) * dnorm(seq(-10,10,0.01),mean(as.numeric(data.use[cell.ident==i])),sd(as.numeric(data.use[cell.ident==i]))),col=i,lwd=2); 
#                }
              }
            }
            ident.table=table(cell.ident)
            #for(i in 2:do.k) ident.table[i]=ident.table[1]
            raw.probs=t(sapply(data.use,function(y) unlist(lapply(1:do.k,function(x) ((ident.table[x]/sum(ident.table))*dnorm(y,mean(as.numeric(data.use[cell.ident==x])),sd(as.numeric(data.use[cell.ident==x]))))))))
            norm.probs=raw.probs/apply(raw.probs,1,sum)
            colnames(norm.probs)=unlist(lapply(1:do.k,function(x)paste(gene,x-1,"post",sep=".")))
            norm.probs=cbind(norm.probs,cell.ident); colnames(norm.probs)[ncol(norm.probs)]=paste(gene,".ident",sep="")
            new.mix.probs=data.frame(minusc(object@mix.probs,paste(gene,".",sep="")),row.names = rownames(object@mix.probs)); colnames(new.mix.probs)[1]="nGene"
            object@mix.probs=cbind(new.mix.probs,norm.probs)
            
            if (do.plot) {
              nCol=2
              #vlnPlot(object,gene,by.k=TRUE,use.imputed=TRUE)
              num.row=floor((do.k+1)/nCol-1e-5)+1
              #par(mfrow=c(num.row,nCol))         
              #par(mfrow=c(1,1))
              hist(as.numeric(data.use),probability = TRUE,ylim=c(0,1),xlab=gene,main=gene);
              for(i in 1:do.k) lines(seq(-10,10,0.01),(ident.table[i]/sum(ident.table)) * dnorm(seq(-10,10,0.01),mean(as.numeric(data.use[cell.ident==i])),sd(as.numeric(data.use[cell.ident==i]))),col=i,lwd=2); 
              #unlist(lapply(1:do.k,function(x) plot(as.numeric(data.use),norm.probs[,x],ylab=paste("Posterior for Component ",x-1,sep=""),main=gene)))     
              #barplot(ident.table,main=round(ident.table[2]/sum(ident.table),3))
            }
            return(object)
          }
)

iter.k.fit=function(scale.data,cell.ident,data.use) {
  means.all=sapply(sort(unique(cell.ident)),function(x)apply(scale.data[,cell.ident==x],1,mean))
  all.dist=data.frame(t(sapply(1:ncol(scale.data),function(x) unlist(lapply(sort(unique(cell.ident)),function(y)dist(rbind(scale.data[,x],means.all[,y])))))))
  cell.ident=apply(all.dist,1,which.min)
  cell.ident=order(tapply(as.numeric(data.use),cell.ident,mean))[cell.ident]
  return(cell.ident)
}


setGeneric("fit.gene.mix", function(object, gene, do.k=3,use.mixtools=TRUE,do.plot=FALSE,plot.with.imputed=TRUE,min.bin.size=10) standardGeneric("fit.gene.mix"))
setMethod("fit.gene.mix", "seurat",
          function(object, gene, do.k=3,use.mixtools=TRUE,do.plot=FALSE,plot.with.imputed=TRUE,min.bin.size=10) {
            require(mixtools)
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

lasso.fxn = function(lasso.input,genes.obs,s.use=20,gene.name=NULL,do.print=FALSE,gram=TRUE) {
  lasso.model=lars(lasso.input,as.numeric(genes.obs),type="lasso",max.steps = s.use*2,use.Gram=gram)
  #lasso.fits=predict.lars(lasso.model,lasso.input,type="fit",s=min(s.use,max(lasso.model$df)))$fit
  lasso.fits=predict.lars(lasso.model,lasso.input,type="fit",s=s.use)$fit
  if (do.print) print(gene.name)
  return(lasso.fits)  
}

setGeneric("addImputedScore", function(object, genes.use=NULL,genes.fit=NULL,s.use=20,do.print=FALSE,gram=TRUE) standardGeneric("addImputedScore"))
setMethod("addImputedScore", "seurat",
          function(object, genes.use=NULL,genes.fit=NULL,s.use=20,do.print=FALSE,gram=TRUE) {
            genes.use=set.ifnull(genes.use,object@var.genes)
            genes.fit=set.ifnull(genes.fit,object@var.genes)
            genes.use=genes.use[genes.use%in%rownames(object@data)]
            genes.fit=genes.fit[genes.fit%in%rownames(object@data)]
            
            lasso.input=t(object@data[genes.use,])
            print(genes.fit)
            lasso.fits=data.frame(t(sapply(genes.fit,function(x)lasso.fxn(t(object@data[genes.use[genes.use!=x],]),object@data[x,],s.use=s.use,x,do.print,gram))))
            genes.old=genes.fit[genes.fit%in%rownames(object@imputed)]
            genes.new=genes.fit[!(genes.fit%in%rownames(object@imputed))]
            
            if (length(genes.old)>0) object@imputed[genes.old,]=lasso.fits[genes.old,]
            object@imputed=rbind(object@imputed,lasso.fits[genes.new,])
            return(object)
          }
)    



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

setGeneric("calcNoiseModels", function(object, cell.ids=NULL, trusted.genes=NULL,n.bin=20,drop.expr=1) standardGeneric("calcNoiseModels"))
setMethod("calcNoiseModels","seurat",
          function(object, cell.ids=NULL, trusted.genes=NULL,n.bin=20,drop.expr=1) {
            object@drop.expr=drop.expr
            cell.ids=set.ifnull(cell.ids,1:ncol(object@data))
            trusted.genes=set.ifnull(trusted.genes,rownames(object@data))
            trusted.genes=trusted.genes[trusted.genes%in%rownames(object@data)]
            object@trusted.genes=trusted.genes
            data=object@data[trusted.genes,]
            data.stats=data.frame(data[,1])
            code_humpAvg=apply(data,1,humpMean,min=object@drop.expr)
            code_humpAvg[code_humpAvg>9]=9
            code_humpAvg[is.na(code_humpAvg)]=0
            data.stats$code_humpAvg=code_humpAvg
            data[data>object@drop.expr]=1
            data[data<object@drop.expr]=0
            data$bin=cut(code_humpAvg,n.bin)
            data$avg=code_humpAvg
            rownames(data.stats)=rownames(data)
            my.coefs=data.frame(t(sapply(colnames(data[1:(ncol(data)-2)]),
                                         getAB,data=data,data2=data.stats,status="code",code2="humpAvg",hasBin=TRUE,doPlot=FALSE)))
            colnames(my.coefs)=c("a","b")
            object@drop.coefs = my.coefs
            return(object)
          }
)  

setGeneric("feature.plot", function(object,genes.plot,pc.1=1,pc.2=2,cells.use=NULL,pt.size=1,do.return=FALSE,do.bare=FALSE,by.k=FALSE,cols.use=NULL,pch.use=16,reduction.use="pca",nCol=NULL) standardGeneric("feature.plot"))
setMethod("feature.plot", "seurat", 
          function(object,genes.plot,pc.1=1,pc.2=2,cells.use=NULL,pt.size=1,do.return=FALSE,do.bare=FALSE,by.k=FALSE,cols.use=NULL,pch.use=16,reduction.use="tsne",nCol=NULL) {
            cells.use=set.ifnull(cells.use,colnames(object@data))
            dim.code="PC"
            if (is.null(nCol)) {
              nCol=2
              if (length(genes.plot)>6) nCol=3
              if (length(genes.plot)>9) nCol=4
            }         
            num.row=floor(length(genes.plot)/nCol-1e-5)+1
            par(mfrow=c(num.row,nCol))
            if (reduction.use=="pca") {
              data.plot=object@pca.rot[cells.use,]
              dim.code="PC"
            }
            if (reduction.use=="tsne") {
              data.plot=object@tsne.rot[cells.use,]
              dim.code="TSNE_"
            }
            if (reduction.use=="ica") {
              data.plot=object@ica.rot[cells.use,]
              dim.code="IC"
            }
            
            stat.use=as.factor(object@data.stat[cells.use])
            if (by.k) stat.use=as.factor(retreiveCluster(object))
            data.plot$stat=stat.use
            x1=paste(dim.code,pc.1,sep=""); x2=paste(dim.code,pc.2,sep="")
            data.plot$x=data.plot[,x1]; data.plot$y=data.plot[,x2]
            data.plot$pt.size=pt.size
            print(head(data.plot))
            data.use=object@data[,cells.use]
            if (length(which(!(genes.plot%in%rownames(data.use))))>0) {
              if (ncol(object@pca.rot)>=2) {
                data.use=data.frame(t(object@pca.rot))
              }
              if (length(ainb(genes.plot,"nGene")>0)) {
                data.use=rbind(data.use,object@data.ngene)
                rownames(data.use)[nrow(data.use)]="nGene"
              }
            }
            
            for(i in genes.plot) {
              data.gene=na.omit(data.frame(data.use[i,]))
              data.cut=as.numeric(as.factor(cut(as.numeric(data.gene),breaks = 10)))
              data.col=rev(heat.colors(10))[data.cut]
              plot(data.plot$x,data.plot$y,col=data.col,cex=pt.size,pch=pch.use,main=i,xlab=x1,ylab=x2)
            }
            rp()
          }
)

setGeneric("pca.plot", function(object,pc.1=1,pc.2=2,cells.use=NULL,pt.size=4,do.return=FALSE,do.bare=FALSE,by.k=FALSE,cols.use=NULL,reduction.use="pca") standardGeneric("pca.plot"))
setMethod("pca.plot", "seurat", 
          function(object,pc.1=1,pc.2=2,cells.use=NULL,pt.size=3,do.return=FALSE,do.bare=FALSE,by.k=FALSE,cols.use=NULL,reduction.use="pca") {
            cells.use=set.ifnull(cells.use,colnames(object@data))
            dim.code="PC"
            if (reduction.use=="pca") data.plot=object@pca.rot[cells.use,]
            if (reduction.use=="tsne") {
              data.plot=object@tsne.rot[cells.use,]
              dim.code="TSNE_"
            }
            if (reduction.use=="ica") {
              data.plot=object@ica.rot[cells.use,]
              dim.code="IC"
            }
            
            stat.use=as.factor(object@data.stat[cells.use])
            if (by.k) stat.use=as.factor(retreiveCluster(object))
            data.plot$stat=stat.use
            x1=paste(dim.code,pc.1,sep=""); x2=paste(dim.code,pc.2,sep="")
            data.plot$x=data.plot[,x1]; data.plot$y=data.plot[,x2]
            data.plot$pt.size=pt.size
            p=ggplot(data.plot,aes(x=x,y=y))+geom_point(aes(colour=factor(stat),size=pt.size))
            if (!is.null(cols.use)) {
              p=p+scale_colour_manual(values=cols.use)
            }
            p2=p+xlab(x1)+ylab(x2)+scale_size(range = c(pt.size, pt.size))
            p3=p2+gg.xax()+gg.yax()+gg.legend.pts(6)+ggplot.legend.text(12)+no.legend.title+theme_bw()+nogrid
            if (do.return) {
              if (do.bare) return(p)
              return(p3)
            }
            print(p3)
          }
)

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
            
            diff.genes=rownames(subset(diffExp.test(p15,marker.cells,embed.diff,genes.use=genes.use),myP<(1e-5)))
            diff.genes=subset(diffExp.test(p15,marker.cells,embed.diff,genes.use = diff.genes),myP<(1e-10))
            return(diff.genes)
          }
)
    

setGeneric("Mclust_dimension", function(object,pc.1=1,pc.2=2,cells.use=NULL,pt.size=4,reduction.use="pca",G.use=NULL,set.stat=FALSE,seed.use=1,...) standardGeneric("Mclust_dimension"))
setMethod("Mclust_dimension", "seurat", 
          function(object,pc.1=1,pc.2=2,cells.use=NULL,pt.size=4,reduction.use="pca",G.use=NULL,set.stat=FALSE,seed.use=1,...) {
            cells.use=set.ifnull(cells.use,colnames(object@data))
            dim.code="PC"
            if (reduction.use=="pca") data.plot=object@pca.rot[cells.use,]
            if (reduction.use=="tsne") {
              data.plot=object@tsne.rot[cells.use,]
              dim.code="TSNE_"
            }
            if (reduction.use=="ica") {
              data.plot=object@ica.rot[cells.use,]
              dim.code="IC"
            }
            x1=paste(dim.code,pc.1,sep=""); x2=paste(dim.code,pc.2,sep="")
            data.plot$x=data.plot[,x1]; data.plot$y=data.plot[,x2]
            set.seed(seed.use); data.mclust=ds <- dbscan(data.plot[,c("x","y")],eps = G.use,...)

            to.set=as.numeric(data.mclust$cluster+1)
            data.names=names(object@data.stat)
            object@data.info[data.names,"m"]=to.set
            if (set.stat) {
              object@data.stat=to.set; names(object@data.stat)=data.names              
            }
            
            return(object)
          }
)

setGeneric("pca.sig.genes", function(object,pcs.use,pval.cut=0.1,use.full=TRUE) standardGeneric("pca.sig.genes"))
setMethod("pca.sig.genes", "seurat", 
          function(object,pcs.use,pval.cut=0.1,use.full=TRUE) {
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
            genes.use=genes.use[genes.use%in%rownames(object@scale.data)]
            return(genes.use)       
          }
)



same=function(x) return(x)

setGeneric("doKMeans", function(object,pcs.use=1,pval.cut=0.1,k.num=NULL,k.seed=1,do.plot=TRUE,clust.cut=2.5,disp.cut=2.5,k.cols=pyCols,do.one=FALSE,do.k.col=FALSE,
                                k.col=NULL,pc.row.order=NULL,pc.col.order=NULL, rev.pc.order=FALSE, cluster.zoom=0, use.full=FALSE,clust.col=TRUE,do.annot=FALSE,
                                only.k.annot=FALSE,do.recalc=TRUE,use.imputed=FALSE,col.annot.show=NULL,genes.use=NULL,print.genes=FALSE) standardGeneric("doKMeans"))

setMethod("doKMeans","seurat",
          function(object,pcs.use=1,pval.cut=0.1,k.num=NULL,k.seed=1,do.plot=TRUE,clust.cut=2.5,disp.cut=2.5,k.cols=pyCols,do.one=FALSE,do.k.col=FALSE,
                   k.col=NULL,pc.row.order=NULL,pc.col.order=NULL,rev.pc.order=FALSE, cluster.zoom=0,use.full=FALSE,clust.col=TRUE,do.annot=FALSE,
                   only.k.annot=FALSE,do.recalc=TRUE,use.imputed=FALSE,col.annot.show=NULL,genes.use=NULL,print.genes=FALSE) {
            require(gplots)
            require(NMF)
            
            data.use.orig=object@scale.data
            if (use.imputed) data.use.orig=data.frame(t(scale(t(object@imputed))))
            data.use=minmax(data.use.orig,min=disp.cut*(-1),max=clust.cut)
            pc.row.order=set.ifnull(pc.row.order,pcs.use[1])
            pc.col.order=set.ifnull(pc.col.order,pcs.use[1])
            
            pvals.use=object@jackStraw.empP
            pcx.use=object@pca.x
            if (use.full)  {
              pvals.use=object@jackStraw.empP.full
              pcx.use=object@pca.x.full
            }
            revFxn=same; if (rev.pc.order) revFxn=function(x)max(x)+1-x;
            k.num=set.ifnull(k.num,object@.best)
            #if (length(pcs.use)==1) pvals.min=pvals.use[,pcs.use]
            #if (length(pcs.use)>1) pvals.min=apply(pvals.use[,pcs.use],1,min)
            #names(pvals.min)=rownames(pvals.use)
            #genes.use=set.ifnull(genes.use,names(pvals.min)[pvals.min<pval.cut])
            genes.use=set.ifnull(genes.use,pca.sig.genes(object,pcs.use,pval.cut,use.full))
            
            genes.use=genes.use[genes.use%in%rownames(data.use)]
            cells.use=colnames(data.use)
            kmeans.data=data.use[genes.use,cells.use]      
            if (print.genes) print(genes.use)
            
            k.col=set.ifnull(k.col, k.num)
            if (do.recalc) {
              set.seed(k.seed); kmeans.obj=kmeans(kmeans.data,k.num); kmeans.col=kmeans(t(kmeans.data),k.col)
              kmeans.obj$cluster=as.numeric(revFxn(rank(tapply(pcx.use[genes.use,pc.row.order],as.numeric(kmeans.obj$cluster),mean)))[as.numeric(kmeans.obj$cluster)])
              names(kmeans.obj$cluster)=genes.use
              kmeans.col$cluster=as.numeric(revFxn(rank(tapply(object@pca.rot[cells.use,pc.col.order],kmeans.col$cluster,mean)))[as.numeric(kmeans.col$cluster)])
              names(kmeans.col$cluster)=cells.use
              object@kmeans.obj=list(kmeans.obj)
              object@kmeans.col=list(kmeans.col)
            }    
            kmeans.obj=object@kmeans.obj[[1]]
            kmeans.col=object@kmeans.col[[1]]
            
            if (do.plot) {       
              disp.data=minmax(kmeans.data[order(kmeans.obj$cluster[genes.use]),],min=disp.cut*(-1),max=disp.cut)
              if (do.one)  {
                disp.data=disp.data[order(pcx.use[rownames(disp.data),pcs.use[1]]),order(object@pca.rot[colnames(disp.data),pcs.use[1]])]
              }
              if (do.k.col) {
                disp.data=disp.data[,order(kmeans.col$cluster)]
              }
              if (cluster.zoom[1] != 0) {
                genes.use=names(kmeans.obj$cluster)[which(kmeans.obj$cluster%in%cluster.zoom)]
                disp.data=disp.data[genes.use[order(kmeans.obj$cluster[genes.use])],]
              }  
              row.annot=data.frame(cbind(kmeans.obj$cluster[rownames(disp.data)],pcx.use[rownames(disp.data),pcs.use]))
              colnames(row.annot)=c("K",paste("PC",pcs.use,sep=""))
              k.names=colnames(disp.data)
              pcrot.info=object@pca.rot[k.names,pcs.use]
              pcrot.info=log(pcrot.info+min(pcrot.info)+1)
              col.annot=data.frame(cbind(object@data.ngene[k.names],pcrot.info,object@data.stat[k.names]))
              #col.annot.data=data.frame(cbind(object@data.ngene[k.names],pcrot.info,object@data.stat[k.names]),t(object@data.metrics),t(object@gene.scores))
              for(i in c(1:(ncol(col.annot)-1))) {
                col.annot[,i]=as.numeric(as.character(col.annot[,i]))
              }     
              colnames(col.annot)=c("nGene",paste("cPC",pcs.use,sep=""),"stat")     
              do.Colv="TRUE"
              colnames(disp.data)=sub.string(colnames(disp.data),"_rsem","")
              if (!do.annot) {
                col.annot=NULL
                row.annot=NULL
              }
              if (do.k.col) {
                col.annot=cbind(col.annot,kmeans.col$cluster[colnames(disp.data)])
                colnames(col.annot)[ncol(col.annot)]="K"
              }
              if (only.k.annot) {
                row.annot=minusc(row.annot,"PC")
                col.annot=minusc(col.annot,"PC")
              }
              if (do.one)  aheatmap(disp.data,annRow = row.annot,annCol=col.annot,Rowv=NA,Colv=NA,col=k.cols)
              if (!(do.one)&& !(do.k.col)) aheatmap(disp.data,annRow = row.annot,annCol=col.annot,Rowv=NA,Colv=clust.col,col=k.cols)
              if (!(do.one)&& (do.k.col)) aheatmap(disp.data,annRow = row.annot,annCol=col.annot,Rowv=NA,Colv=NA,col=k.cols)
              
            }
            return(object)
          }
)

setGeneric("genes.in.cluster", function(object, cluster.num)  standardGeneric("genes.in.cluster"))
setMethod("genes.in.cluster", signature = "seurat",
          function(object, cluster.num) {
            print(unlist(lapply(cluster.num,function(x)sort(names(which(object@kmeans.obj[[1]]$cluster==x))))))
          }    
)

setGeneric("cells.in.cluster", function(object, cluster.num)  standardGeneric("cells.in.cluster"))
setMethod("cells.in.cluster", signature = "seurat",
          function(object, cluster.num) {
#            print(sort(names(which(object@kmeans.col[[1]]$cluster==cluster.num))))
            return(sort(unlist(lapply(cluster.num,function(x) names(which(object@kmeans.col[[1]]$cluster==x))))))
          }    
)

setGeneric("cell.cor.matrix", function(object, cor.genes=NULL,cell.inds=NULL, do.k=FALSE,k.seed=1,k.num=4,vis.low=(-1),vis.high=1,vis.one=0.8,pcs.use=1:3,col.use=pyCols)  standardGeneric("cell.cor.matrix"))
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

setGeneric("gene.cor.matrix", function(object, cor.genes=NULL,cell.inds=NULL, do.k=FALSE,k.seed=1,k.num=4,vis.low=(-1),vis.high=1,vis.one=0.8,pcs.use=1:3,col.use=pyCols)  standardGeneric("gene.cor.matrix"))
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

setGeneric("calinskiPlot", function(object,pcs.use,pval.cut=0.1,gene.max=15,col.max=25,use.full=TRUE)  standardGeneric("calinskiPlot"))
setMethod("calinskiPlot","seurat",
          function(object,pcs.use,pval.cut=0.1,gene.max=15,col.max=25,use.full=TRUE) {
            require(vegan)
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

setMethod("show", "seurat",
          function(object) {
            cat("An object of class ", class(object), " in project ", object@project.name, "\n", sep = "")
            cat(" ", nrow(object@data), " features by ",
                ncol(object@data), " samples.\n", sep = "")
            invisible(NULL)
          }
)

plot.Vln=function(gene,data,cell.stat,ylab.max=12,do.ret=FALSE,do.sort=FALSE,size.x.use=16,size.y.use=16,size.title.use=20,adjust.use=1,size.use=1,cols.use=NULL) {
  data$gene=as.character(rownames(data))
  data.use=data.frame(data[gene,])
  if (length(gene)==1) {
    data.melt=data.frame(rep(gene,length(cell.stat))); colnames(data.melt)[1]="gene"
    data.melt$value=as.numeric(data[1,1:length(cell.stat)])
    data.melt$id=names(data)[1:length(cell.stat)]
  }
  if (length(gene)>1) data.melt=melt(data.use,id="gene")
  data.melt$stat=cell.stat
  noise <- rnorm(length(data.melt$value))/100000
  data.melt$value=as.numeric(as.character(data.melt$value))+noise
  if(do.sort) {
    data.melt$stat=factor(data.melt$stat,levels=names(rev(sort(tapply(data.melt$value,data.melt$stat,mean)))))
  }
  p=ggplot(data.melt,aes(factor(stat),value))
  p2=p + geom_violin(scale="width",adjust=adjust.use,trim=TRUE,aes(fill=factor(stat))) + ylab("Expression level (log TPM)")
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

setGeneric("dot.plot", function(object,genes.plot,cex.use=2,cols.use=NULL,thresh.col=2.5,dot.min=0.05)  standardGeneric("dot.plot"))
setMethod("dot.plot","seurat",
          function(object,genes.plot,cex.use=2,cols.use=NULL,thresh.col=2.5,dot.min=0.05) {
            genes.plot=ainb(genes.plot,rownames(object@data))
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
            plot(data.x,data.y,cex=data.cex,pch=16,col=exp.col,xaxt="n",xlab="",ylab="Cluster")
            axis(1,at = 1:length(genes.plot),genes.plot)
          }

) 

setGeneric("vlnPlot", function(object,genes.plot,nCol=NULL,ylab.max=12,do.ret=TRUE,do.sort=FALSE,
                               size.x.use=16,size.y.use=16,size.title.use=20, by.k=FALSE,use.imputed=FALSE,adjust.use=1,size.use=1,cols.use=NULL,...)  standardGeneric("vlnPlot"))
setMethod("vlnPlot","seurat",
          function(object,genes.plot,nCol=NULL,ylab.max=12,do.ret=FALSE,do.sort=FALSE,size.x.use=16,size.y.use=16,size.title.use=20,by.k=FALSE,use.imputed=FALSE,adjust.use=1,
                   size.use=1,cols.use=NULL,...) {
            if (is.null(nCol)) {
              nCol=2
              if (length(genes.plot)>6) nCol=3
              if (length(genes.plot)>9) nCol=4
            }
            data.use=object@data
            if (use.imputed) data.use = object@imputed
            if (length(which(!(genes.plot%in%rownames(data.use))))>0) {
              if (ncol(object@pca.rot)>=2) {
                data.use=data.frame(t(object@pca.rot))
              }
              if (length(ainb(genes.plot,"nGene")>0)) {
                data.use=rbind(data.use,object@data.ngene)
                rownames(data.use)[nrow(data.use)]="nGene"
              }
            }
            data.use=na.omit(data.frame(data.use[genes.plot,]))
            if (length(which(!(genes.plot%in%rownames(data.use))))>0) {
              print(paste("ERROR ", genes.plot[which(!(genes.plot%in%rownames(data.use)))], " not found", sep=""))
            }
            stat.use=object@data.stat
            if (by.k) stat.use=retreiveCluster(object)
            pList=lapply(genes.plot,function(x) plot.Vln(x,data.use[x,],stat.use,ylab.max,TRUE,do.sort,size.x.use,size.y.use,size.title.use,adjust.use,size.use,cols.use))
  
            if(do.ret) {
              return(pList)
            }
            else {
              multiplotList(pList,cols = nCol)
              rp()
            }
          }
)  

setGeneric("jackStrawPlot", function(object,plot.lim=0.4,num.pc=5,score.thresh=1e-5,...)  standardGeneric("jackStrawPlot"))
setMethod("jackStrawPlot","seurat",
          function(object,plot.lim=0.4,num.pc=5,score.thresh=0.01,...) {
            pAll=object@jackStraw.empP
            nCol=2
            if (num.pc>6) nCol=3
            if (num.pc>9) nCol=4
            num.row=floor(num.pc/nCol-1e-5)+1
            par(mfrow=c(nCol,num.row))
            qq.x=sapply(1:num.pc,function(x) {qqplot(pAll[,x],runif(1000),plot.it=FALSE)$x})
            qq.y=sapply(1:num.pc,function(x) {qqplot(pAll[,x],runif(1000),plot.it=FALSE)$y})
            pc.score=unlist(lapply(1:num.pc, function(x) mean(qq.y[which(qq.x[,x]<=score.thresh),x])))
            pc.score=unlist(lapply(1:num.pc, function(x) prop.test(c(length(which(pAll[,x]<=score.thresh)),floor(nrow(pAll)*score.thresh)),c(nrow(pAll),nrow(pAll)))$p.val))
            
           unlist(lapply(1:num.pc,function(x) {qqplot(pAll[,x],runif(1000),,xlim=c(0,plot.lim),ylim=c(0,plot.lim),xlab=colnames(pAll)[x],pch=16,main=round(pc.score[x],4)); 
                                                lines(seq(0,1,0.01),seq(0,1,0.01),lty=2,lwd=2)}))
#            unlist(lapply(1:num.pc,function(x) {qqplot(pAll[,x],runif(1000),,xlim=c(0,plot.lim),ylim=c(0,plot.lim),xlab=colnames(pAll)[x],pch=16); 
#                                                lines(seq(0,1,0.01),seq(0,1,0.01),lty=2,lwd=2)}))
            rp()
          }
)


setGeneric("jackStrawPlot2", function(object,plot.lim=0.4,num.pc=5,...)  standardGeneric("jackStrawPlot2"))
setMethod("jackStrawPlot2","seurat",
          function(object,plot.lim=0.4,num.pc=5,...) {
            pAll=object@emp.pval
            num.pc=ncol(pAll)
            num.row=floor(num.pc/2-1e-5)+1
            par(mfrow=c(2,num.row))
            unlist(lapply(1:num.pc,function(x) {qqplot(pAll[,x],runif(1000),,xlim=c(0,plot.lim),ylim=c(0,plot.lim),xlab=colnames(pAll)[x],pch=16); 
                                                lines(seq(0,1,0.01),seq(0,1,0.01),lty=2,lwd=2)}))
            rp()
          }
)

setGeneric("genePlot", function(object, gene1, gene2, cell.ids=NULL,col.use=NULL,
                                pch.use=16,cex.use=2,no.col=FALSE,inc.first=FALSE,use.imputed=FALSE,by.k=FALSE,do.ident=FALSE,...)  standardGeneric("genePlot"))
setMethod("genePlot","seurat",
          function(object, gene1, gene2, cell.ids=NULL,col.use=NULL,
                   pch.use=16,cex.use=2,no.col=FALSE,inc.first=FALSE,use.imputed=FALSE,by.k=FALSE,do.ident=FALSE,...) {
            data.use=object@data
            if (use.imputed) data.use= object@imputed
            if (length(which(!(c(gene1,gene2)%in%rownames(data.use))))>0) {
              if (ncol(object@pca.rot>=2)) {
                data.use=rbind(data.use,t(object@pca.rot))
              }
              if (length(ainb(genes.plot,"nGene")>0)) {
                data.use=rbind(data.use,object@data.ngene)
                rownames(data.use)[nrow(data.use)]="nGene"
              }
            }
            
            cell.ids=set.ifnull(cell.ids,colnames(data.use))
            g1=as.numeric(data.use[gene1,cell.ids])
            g2=as.numeric(data.use[gene2,cell.ids])
            stat.use=as.factor(object@data.stat[cell.ids])
            if (by.k) stat.use=as.factor(retreiveCluster(object))
            if (length(col.use)>1) {
              col.use=col.use[as.numeric(stat.use)]
            }
            else {
              col.use=set.ifnull(col.use,as.numeric(stat.use))
            }
            gene.cor=round(cor(g1,g2),2)
            plot(g1,g2,xlab=gene1,ylab=gene2,col=col.use,cex=cex.use,main=gene.cor,pch=pch.use)
            if (do.ident) {
              return(identify(g1,g2,labels = cell.ids))
            }
          }
)

setGeneric("geneScorePlot", function(object, gene1, score.name, cell.ids=NULL,col.use=NULL,nrpoints.use=Inf,pch.use=16,cex.use=2,...)  standardGeneric("geneScorePlot"))
setMethod("geneScorePlot","seurat",
          function(object, gene1, score.name, cell.ids=NULL,col.use=NULL,nrpoints.use=Inf,pch.use=16,cex.use=2,...) {
            cell.ids=set.ifnull(cell.ids,colnames(object@data))
            g1=as.numeric(object@data[gene1,cell.ids])
            my.score=retreiveScore(object,score.name)
            s1=as.numeric(my.score[cell.ids])
            col.use=set.ifnull(as.numeric(as.factor(object@data.stat[cell.ids])))
            gene.cor=round(cor(g1,s1),2)
            smoothScatter(g1,s1,xlab=gene1,ylab=score.name,col=col.use,nrpoints=nrpoints.use,cex=cex.use,main=gene.cor,pch=pch.use)
          }
)

setGeneric("cellPlot", function(object, cell1, cell2, gene.ids=NULL,col.use="black",nrpoints.use=Inf,pch.use=16,cex.use=0.5,do.ident=FALSE,...)  standardGeneric("cellPlot"))
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

setGeneric("jackStraw", function(object,num.pc=8,num.replicate=100,prop.freq=0.01,do.print=FALSE)  standardGeneric("jackStraw"))
setMethod("jackStraw","seurat",
          function(object,num.pc=5,num.replicate=100,prop.freq=0.01,do.print=FALSE) {
            pc.genes=rownames(object@pca.x)
            if (length(pc.genes)<200) prop.freq=max(prop.freq,0.015)
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

jackRandom=function(scaled.data,prop.use=0.01,r1.use=1,r2.use=5, seed.use=1) {
  set.seed(seed.use); rand.genes=sample(rownames(scaled.data),nrow(scaled.data)*prop.use)
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


setGeneric("mean.var.plot", function(object, fxn.x=humpMean, fxn.y=sd,do.plot=TRUE,set.var.genes=TRUE,do.text=TRUE,
                                     x.low.cutoff=4,x.high.cutoff=8,y.cutoff=2,y.high.cutoff=12,cex.use=0.5,cex.text.use=0.5,do.spike=FALSE, 
                                     pch.use=16, col.use="black", spike.col.use="red",use.imputed=FALSE,plot.both=FALSE,do.contour=TRUE,
                                     contour.lwd=3, contour.col="white", contour.lty=2) standardGeneric("mean.var.plot"))
setMethod("mean.var.plot", signature = "seurat",
          function(object, fxn.x=humpMean, fxn.y=sd,do.plot=TRUE,set.var.genes=TRUE,do.text=TRUE,
                   x.low.cutoff=4,x.high.cutoff=8,y.cutoff=1,y.high.cutoff=12,cex.use=0.5,cex.text.use=0.5,do.spike=FALSE, 
                   pch.use=16, col.use="black", spike.col.use="red",use.imputed=FALSE,plot.both=FALSE,do.contour=TRUE,
                   contour.lwd=3, contour.col="white", contour.lty=2) {
            data=object@data
            data.x=apply(data,1,fxn.x); data.y=apply(data,1,fxn.y)
            data.norm.y=meanNormFunction(data,fxn.x,fxn.y)
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




