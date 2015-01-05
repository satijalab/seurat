nmf.options(grid.patch=TRUE)
nogrid=theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

sort.column=function(x, col) {
  return(x[order(x[,col]),])
}

tsplot=function(object,x=1,cex.use=0.6) {
  cols.use=rainbow(length(levels(object@ident))); cols.use[x]="lightgrey"
  plot(object@tsne.rot[,1],object@tsne.rot[,2],col=cols.use[as.integer(object@ident)],pch=16,xlab="TSNE_1",ylab="TSNE_2",cex=cex.use)
  k.centers=t(sapply(levels(object@ident),function(x) apply(object@tsne.rot[which.cells(object,x),],2,mean)))
  points(k.centers[,1],k.centers[,2],cex=1.3,col="white",pch=16); text(k.centers[,1],k.centers[,2],levels(object@ident),cex=1)
}

getLeftDecendants=function(tree,node) {
  daughters<-tree$edge[which(tree$edge[,1]==node),2]
  if (daughters[1] <= (tree$Nnode+1)) return(daughters[1])
  daughter.use=getDescendants(tree,daughters[1])
  daughter.use=daughter.use[daughter.use<=(tree$Nnode+1)]
  return(daughter.use)
}

getRightDecendants=function(tree,node) {
  daughters<-tree$edge[which(tree$edge[,1]==node),2]
  if (daughters[2] <= (tree$Nnode+1)) return(daughters[2])
  daughter.use=getDescendants(tree,daughters[2])
  daughter.use=daughter.use[daughter.use<=(tree$Nnode+1)]
  return(daughter.use)
}

getDescendants<-function(tree,node,curr=NULL){
  if(is.null(curr)) curr<-vector()
  daughters<-tree$edge[which(tree$edge[,1]==node),2]
  curr<-c(curr,daughters)
  w<-which(daughters>=length(tree$tip))
  if(length(w)>0) for(i in 1:length(w)) 
    curr<-getDescendants(tree,daughters[w[i]],curr)
  return(curr)
}

situ3d=function(data, label=NULL, ...) {
  # Call Seurat function to get the in situ values out.
  exp.1=data
  exp.1=(exp.1-min(exp.1))/(max(exp.1)-min(exp.1))
  # Reformat them into an expression matrix as expected by the plotting function
  expression.matrix <- data.frame(matrix(exp.1, nrow=8, ncol=8))
  rownames(expression.matrix) <- c("24-30", "17-23", "13-16", "9-12", "7-8", "5-6", "3-4", "1-2")
  names(expression.matrix) <- c("1-4", "5-8", "9-12", "13-16", "17-20", "21-24", "25-28", "29-32")
  
  # Call the plotting function.
  zf.insitu.side(expression.matrix)
  par3d(windowRect=c(0, 0, 800, 800))
  
  # Label or not and then set the view.
  if (!(is.null(label))) {
    text3d(x=0, y=0, z=1.5, text=label, cex=3)
  }
  view3d(zoom=.75, theta=0, phi=-90, fov=0)
}

aucFxn=function(preds,truth,do.plot=FALSE,lab.main="",...) {
  pred.use=prediction(preds,truth,0:1)
  perf.use1=performance(pred.use,"tpr","fpr")
  perf.use=performance(pred.use,"auc")
  auc.use=round(perf.use@y.values[[1]],3)
  if (do.plot) plot(perf.use1@x.values[[1]],perf.use1@y.values[[1]],type="l",xlab="FPR",ylab="TPR",main=paste(lab.main, auc.use))
  return(auc.use)
}

humpMean=function(x, min=0) {
  return(mean(x[x>min]))
}

humpVar=function(x, min=0) {
  return(var(x[x>min]))
} 

debugdmvnorm=function (x, mean = rep(0, p), sigma = diag(p), log = FALSE) 
{
  if (is.vector(x)) 
    x <- matrix(x, ncol = length(x))
  p <- ncol(x)
  dec <- tryCatch(chol(sigma), error = function(e) e)
  tmp <- backsolve(dec, t(x) - mean, transpose = TRUE)
  rss <- colSums(tmp^2)
  logretval <- -(log(diag(dec))) - 0.5 * p * log(2 * pi) - 0.5 * rss
  names(logretval) <- names(x)
  logretval
}

slimdmvnorm=function (x, mean = rep(0, p), sigma = diag(p), log = FALSE) 
{
  x <- matrix(x, ncol = length(x))
  p <- ncol(x)
  dec <- tryCatch(chol(sigma), error = function(e) e)
  tmp <- backsolve(dec, t(x) - mean, transpose = TRUE)
  rss <- colSums(tmp^2)
  logretval <- -sum(log(diag(dec))) - 0.5 * p * log(2 * pi) - 0.5 * rss
  names(logretval) <- rownames(x)
  logretval
}

slimdmvnorm_nosum=function (x, mean = rep(0, p), sigma = diag(p), log = FALSE) 
{
  x <- matrix(x, ncol = length(x))
  p <- ncol(x)
  dec <- tryCatch(chol(sigma), error = function(e) e)
  tmp <- backsolve(dec, t(x) - mean, transpose = TRUE)
  rss <- colSums(tmp^2)
  logretval <- -(log(diag(dec))) - 0.5 * p * log(2 * pi) - 0.5 * rss
  names(logretval) <- rownames(x)
  logretval
}

compareBins=function(object,cell.use,bin.1,bin.2,bins.mu,bins.cov) {
  num.genes=ncol(bins.mu)
  genes.use=colnames(bins.mu)
  to.par=floor(sqrt(num.genes))+1
  par(mfrow=c(to.par,to.par))
  data.use=object@imputed[genes.use,cell.use]; names(data.use)=genes.use
  lik.diff=sort(unlist(lapply(genes.use,function(g)dnorm(data.use[g],bins.mu[bin.1,g],sqrt(bins.cov[[bin.1]][g,g]))/dnorm(data.use[g],bins.mu[bin.2,g],sqrt(bins.cov[[bin.2]][g,g])))))
  for(g in names(lik.diff)) {
    plot(0,0,type="n",xlim=c(0,8),ylim=c(0,2),main=paste(g, round(lik.diff[g],2)))
    lines(density(rnorm(10000,bins.mu[bin.1,g],sqrt(bins.cov[[bin.1]][g,g]))),lwd=2,col="black")
    lines(density(rnorm(10000,bins.mu[bin.2,g],sqrt(bins.cov[[bin.2]][g,g]))),lwd=2,col="red")
    points(data.use[g],0,cex=1.6,col="darkgreen",pch=16)
  }
  #rp()
}

subr=function(data,code) {
  return(data[grep(code,rownames(data)),])
}

cv=function(x)sd(x)/mean(x)

humpCt=function(x, min=0) {
  return(length(x[x>min]))
}

getAB=function(cn="lps_t1_S1_rsem",code="lps_t1",data=cell,data2=cs,code2="avg",status="",ncut=25,hasBin=FALSE,doPlot=FALSE,myfunc=plot,func2=lines,...) {
  if (status == "") {
    status=getStatus(cn)
  }
  if (!(hasBin)) {
    data[data>1]=1
    data[data<1]=0
    data$avg=data2[rownames(data),paste(status,"_",code2,sep="")]
    data[data>9]=9
    data$bin=cut(data$avg,20)
  }
  data$val=data[,cn]
  x1=(tapply(data[,cn],data$bin,mean))
  x2=(tapply(data[,"avg"],data$bin,mean))
  #glm.out = glm(val~avg,data=data,family=binomial)
  glm.out = glm(x1~x2,data=data,family=binomial)
  
  if (doPlot==TRUE) {
    
    pred=(predict(glm.out,data.frame(avg=as.numeric(x2)),type="response"))
    myfunc(x2,x1,pch=16,xlab="Log TPM",ylab="Detection Rate",main=cn,ylim=c(0,1))
    func2(x2,pred,...)
  }
  return(glm.out$coefficients)
}

crushRNANewNoLog=function(mynames,mycodes) {
  r=c()
  for(i in 1:length(mynames)) {
    name=mynames[i]
    file = paste("~/big/",name,"/",name,".rsem.all.6.res",sep="")
    input=read.table(file,header=TRUE)
    #input[,8:ncol(input)]=log(input[,8:ncol(input)]+1)
    code=mycodes[i]
    colnames(input)[8:ncol(input)]=paste(code,colnames(input)[8:ncol(input)],sep="")
    if (i==1) {
      r=input
    }
    if (i>1) {
      r=cbind(r,input[,8:ncol(input)])
    }
    print(name)
  }
  return(r)
}

log_add=function(x) {
  mpi=max(x)
  return(mpi+log(sum(exp(x-mpi))))
}

minusr=function(data,code) {
  matchCode=rownames(data)[grep(code,rownames(data))]
  toIgnore=which(rownames(data) %in% matchCode)
  if (length(toIgnore)==0) return(data)
  toRet=data.frame(data[-toIgnore,])
  rownames(toRet)=rownames(data)[-toIgnore]
  colnames(toRet)=colnames(data)
  return(toRet)                   
}

minusc=function(data,code) {
  matchCode=colnames(data)[grep(code,colnames(data))]
  toIgnore=which(colnames(data) %in% matchCode)
  if (length(toIgnore)==0) return(data)
  return(data[,-toIgnore])
}

getUMI=function(proj,code,cols=4,fullproj="none",nameCol=7,nameCode="rsem",subCode="umi",rmLast=FALSE, verbose=TRUE) {
  myCmd = paste("find ~/big/",proj,"/bam -name ", code, sep="")  
  if (fullproj != "none") {
    myCmd = paste(fullproj," -name ", code, sep="")  
  }
  myUMI=system(myCmd,intern=TRUE)
  if (verbose) print(myUMI)
  myNames=unlist(lapply(myUMI,function(x)gsub(nameCode,subCode,strsplit(x,"/")[[1]][nameCol])))
  data=read.table(myUMI[1],sep="\t",fill=TRUE)[,c(1,cols)]
  if (rmLast) {
    data=data[-nrow(data),]
  }
  colnames(data)[ncol(data)]=myNames[1]
  rownames(data)=as.character(data$V1)
  for(i in 2:length(myUMI)) {
    if (verbose) print(paste(myUMI[i],i))
    #print(data)
    if ((file.info(myUMI[i])$size) > 0) {
      data2=read.table(myUMI[i],sep="\t",fill=TRUE)[,c(1,cols)]
      if (rmLast) {
        data2=data2[-nrow(data2),]
      }
      data=merge(data,data2,by="V1",all=TRUE)
      colnames(data)[ncol(data)]=myNames[i]
    }
  }
  #data2=data
  data[is.na(data)]=0
  data2=aggregate(data[,-1],list(Gene = toupper(data$V1)),sum)
  rownames(data2)=toupper(data2$Gene)
  colnames(data2)=unlist(lapply(colnames(data2),function(x)gsub("-","_",x)))
  data2=subc(data2,subCode)
  return(data2)
}

ainb=function(a,b) {
  a2=a[a%in%b]
  return(a2)
}
meanNormFunction=function(data,myfuncX,myfuncY,nBin=20) {
  data_x=apply(data,1,myfuncX)
  data_y=apply(data,1,myfuncY)
  data_x_bin=cut(data_x,nBin)
  names(data_x_bin)=names(data_x)
  mean_y=tapply(data_y,data_x_bin,mean)
  sd_y=tapply(data_y,data_x_bin,sd)
  return((data_y-mean_y[as.numeric(data_x_bin)])/sd_y[as.numeric(data_x_bin)])
}

shift.cell =function(bin,x,y) {
  bin.y=(bin-1)%/%8+1
  bin.x=(bin-1)%%8+1
  new.x=minmax(bin.x+x,min = 1,max=8)
  new.y=minmax(bin.y+y,min = 1,max=8)
  new.bin=8*(new.y-1)+new.x
  return(new.bin)
}

empP=function(x,nullval) {
  return(length(which(nullval>x))/length(nullval))
}

neighbor.cells=function(bin) {
  return(unique(c(bin,shift.cell(bin,0,1),shift.cell(bin,1,0),shift.cell(bin,-1,0),shift.cell(bin,0,-1))))
}

all.neighbor.cells=function(bin,dist=1) {
  all.comb=expand.grid(rep(list(-dist:dist), 2)) 
  return(unique(unlist(lapply(1:nrow(all.comb),function(x)shift.cell(bin,all.comb[x,1],all.comb[x,2])))))
}

no.legend.title=theme(legend.title=element_blank())
ggplot.legend.text=function(x=12,y="bold") return(theme(legend.text = element_text(size = x, face = y)))
gg.legend.pts=function(x=6) guides(colour = guide_legend(override.aes = list(size=x)))
gg.xax=function(x=16,y="#990000",z="bold",x2=12) return(theme(axis.title.x = element_text(face=z, colour=y, size=x), 
                                                              axis.text.x  = element_text(angle=90, vjust=0.5, size=x2)))

gg.yax=function(x=16,y="#990000",z="bold",x2=12) return(theme(axis.title.y = element_text(face=z, colour=y, size=x), 
                                                              axis.text.y  = element_text(angle=90, vjust=0.5, size=x2)))

sub.string=function(x,s1,s2) return(unlist(lapply(x,function(y)gsub(s1,s2,y))))

jackStrawF=function(prop=0.1,myR1,myR2=3,data=smD) {
  randGenes=sample(rownames(data),nrow(data)*prop)
  smD.mod=data
  smD.mod[randGenes,]=shuffleMatRow(data[randGenes,])
  fmd.pca=prcomp(smD.mod)
  fmd.x=fmd.pca$x
  fmd.rot=fmd.pca$rotation
  fakeF=unlist(lapply(randGenes,jackF,r1=myR1,r2=myR2,x=fmd.x,rot=fmd.rot))
  return(fakeF)
}

jackF=function(gene,r1=1,r2=2,x=md.x,rot=md.rot) {
  if (r2==1) { #assuming r1, r2=1
    mod.x=x[,r1]
    mod.x[gene]=0
    return(var.test((x[,r1]%*%t(rot[,r1])),(mod.x%*%t(rot[,r1])))$statistic)
  }
  mod.x=x[,1:r2]
  mod.x[gene,r1:r2]=rep(0,r2-r1+1)
  return(var.test((x[,1:r2]%*%t(rot[,1:r2])),(mod.x[,1:r2]%*%t(rot[,1:r2])))$statistic)
}

shuffleMatRow=function(x) {
  x2=x
  x2 <- t(x)
  ind <- order(c(col(x2)), runif(length(x2)))
  x2 <- matrix(x2[ind], nrow=nrow(x), ncol=ncol(x), byrow=TRUE)
  return(x2)
}




logMeanMinus= function(x)log(mean(exp(as.numeric(x))-1)+1)
logVarMinus= function(x)(mean(var(as.numeric(x))-1)+1)

logVarMinus2= function(x)(var(exp(as.numeric(x))-1)+1)

quickRNAHuman=function(x) {
  dataFile=paste("~/big/",x,"/summary/",x,".expMatrix.txt",sep="")
  data=log(read.table(dataFile,header=TRUE,sep="\t")[,-1]+1)
  data=subc(data,"rsem")
  return(data)
}

expVar=function(x) {
  return(log(var(exp(x)-1)+1))
}

expSD=function(x) {
  return(log(sd(exp(x)-1)+1))
}

expMean=function(x) {
  return(log(mean(exp(x)-1)+1))
}

quickRNAZfish=function(x) {
  zdict=read.table("~/window/annotate/danRer7ensgene110512.txt.short.txt.tbl.txt",sep="\t",header=TRUE)
  rownames(zdict)=as.character(zdict$name)
  #zcan=subset(zdict,iscanonical==TRUE)
  #zdt=zcan[!duplicated(toupper(zcan$geneSymbol)),]
  #rownames(zdt)=toupper(zdt$geneSymbol)
  #zdt=zdt[,1:8]
  dataFile=paste("~/big/",x,"/",x,".rsem.iso.all.3.res",sep="")
  rawFin=read.table(dataFile,sep="\t",header=TRUE)
  rawFin$gene=toupper(zdict[as.character(rawFin$V1),"geneSymbol"])
  rF=rawFin
  rF=subc(rawFin,"rsem|gene")
  rF=aggregate(subc(rawFin,"rsem"),list(Gene = rawFin$gene),sum)
  rownames(rF)=rF$Gene
  fin=log(subc(rF,"rsem")*1e6+1)
  return(fin)
}

normal.sample=function(x) {
  return(rnorm(10000,mean(x),sd(x)))
  
}

fetch.closest=function(bin,all.centroids,num.cell) {
  bin.y=(bin-1)%/%8+1
  bin.x=(bin-1)%%8+1
  all.centroids=rbind(all.centroids,c(bin.x,bin.y))
  all.dist=as.matrix(dist(all.centroids))
  return(names(sort(all.dist[nrow(all.dist),]))[2:(num.cell+2)])
}

fetch.mincells=function(bin,cells.max,min.cells) {
  for(i in 1:5) {
    my.names=names(ainb(cells.max,all.neighbor.cells(bin,i)))
    if (length(my.names) > min.cells) break;
  }
  return(my.names)
}

cell.centroid=function(cell.probs) {
  centroid.x=round(sum(sapply(1:64,function(x)(x-1)%%8+1)*cell.probs))
  centroid.y=round(sum(sapply(1:64,function(x)(x-1)%/%8+1)*cell.probs))
  centroid.bin=8*(centroid.y-1)+centroid.x
  return(centroid.bin)
}

cell.centroid.x=function(cell.probs) {
  return(centroid.x=round(sum(sapply(1:64,function(x)(x-1)%%8+1)*cell.probs)))
}

cell.centroid.y=function(cell.probs) {
  return(centroid.y=round(sum(sapply(1:64,function(x)(x-1)%/%8+1)*cell.probs)))
}

exact.cell.centroid=function(cell.probs) {
  centroid.x=(sum(sapply(1:64,function(x)(x-1)%%8+1)*cell.probs))
  centroid.y=(sum(sapply(1:64,function(x)(x-1)%/%8+1)*cell.probs))
  return(c(centroid.x,centroid.y))
}


marker.auc.test=function(data1,data2,mygenes) {
  myAUC=unlist(lapply(mygenes,function(x)diffAUC(as.numeric(data1[x,]),as.numeric(data2[x,]))))
  myAUC[is.na(myAUC)]=0
  myDiff=unlist(lapply(mygenes,function(x)(expMean(as.numeric(data1[x,]))-expMean(as.numeric(data2[x,])))))
  toRet=data.frame(cbind(myAUC,myDiff),row.names=mygenes)
  toRet=toRet[rev(order(toRet$myAUC)),]
  return(toRet)
}  

#credit to Cole Trapnell for this
tobit_fitter <- function(x, modelFormulaStr, lower=1, upper=Inf){
  tryCatch({
    FM_fit <-  suppressWarnings(vgam(as.formula(modelFormulaStr), family=tobit(Lower=lower, Upper=upper),data = x))
    FM_fit
  }, 
  #warning = function(w) { FM_fit },
  error = function(e) { NULL }
  )
}

anotinb=function(x,y) {
  x2=x[!x%in%y]
  return(x2)
}
diffTobit=function(x1,x2,lower=1,upper=Inf) {
  my.df=data.frame(c(x1,x2),c(rep(0,length(x1)),rep(1,length(x2))))
  colnames(my.df)=c("Expression","Stat")
  #model.v1=vgam(Expression~1,family = tobit(Lower = lower,Upper = upper),data = my.df)
  model.v1=tobit_fitter(my.df,"Expression~1",lower,upper)
  #model.v2=vgam(Expression~Stat+1,family = tobit(Lower = lower,Upper = upper),data = my.df)
  model.v2=tobit_fitter(my.df,"Expression~Stat+1",lower,upper)
  p=1
  if (is.null(model.v1) == FALSE && is.null(model.v2) == FALSE) {
      (p <- pchisq(2 * (logLik(model.v2) - logLik(model.v1)), df = 1, lower.tail = FALSE))
  }
  return(p)
}

tobit.diffExp.test=function(data1,data2,mygenes) {
  myP=unlist(lapply(mygenes,function(x)diffTobit(as.numeric(data1[x,]),as.numeric(data2[x,]))))
  myP[is.na(myP)]=1
  myDiff=unlist(lapply(mygenes,function(x)(expMean(as.numeric(data1[x,]))-expMean(as.numeric(data2[x,])))))
  toRet=data.frame(cbind(myP,myDiff),row.names=mygenes)
  toRet=toRet[order(toRet$myP),]
  return(toRet)
}

bimod.diffExp.test=function(data1,data2,mygenes) {
  myP=unlist(lapply(mygenes,function(x)diffLRT(as.numeric(data1[x,]),as.numeric(data2[x,]))))
  myP[is.na(myP)]=1
  myDiff=unlist(lapply(mygenes,function(x)(expMean(as.numeric(data1[x,]))-expMean(as.numeric(data2[x,])))))
  toRet=data.frame(cbind(myP,myDiff),row.names=mygenes)
  toRet=toRet[order(toRet$myP),]
  return(toRet)
}

diffAUC = function(x,y) {
  prediction.use=prediction(c(x,y),c(rep(1,length(x)),rep(0,length(y))),0:1)
  perf.use=performance(prediction.use,"auc")
  auc.use=round(perf.use@y.values[[1]],3)
  return(auc.use)
}

diffLRT = function(x,y,xmin=1) {
  lrtX=bimodLikData(x)
  lrtY=bimodLikData(y)
  lrtZ=bimodLikData(c(x,y))
  lrt_diff=2*(lrtX+lrtY-lrtZ)
  return(1-pchisq(lrt_diff,3))
}

bimodLikData=function(x,xmin=0) {
  x1=x[x<=xmin]
  x2=x[x>xmin]
  xal=minmax(length(x2)/length(x),min=1e-5,max=(1-1e-5))
  likA=length(x1)*log(1-xal)
  mysd=sd(x2)
  if(length(x2)<2) {
    mysd=1
  }
  likB=length(x2)*log(xal)+sum(dnorm(x2,mean(x2),mysd,log=TRUE))
  return(likA+likB)
}

makeAlnPlot=function(proj) {
  alnFile=paste("~/big/",proj,"/summary/",proj,".all.aln.metrics.txt",sep="")
  alnData=read.table(alnFile)
  par(mar=c(10,5,4,1))
  mymax=max(100,max(apply(alnData[,4:7],1,sum)))
  x=barplot(as.matrix(t(alnData[,4:7])),names.arg=alnData$V1,las=2,col=1:4,ylim=c(0,mymax))
  text(x,mymax-10,alnData$V2)
  rp()
}

getCoefs=function(data) {
  my_stats=data.frame(data[,1])
  my_stats$code_humpAvg=apply(data,1,humpMean,min=1)
  rownames(my_stats)=rownames(data)
  my_coefs=data.frame(t(sapply(colnames(data),getAB,data=data,data2=my_stats,status="code",code2="humpAvg",doPlot=FALSE)))
  colnames(my_coefs)=c("a","b")
  return(my_coefs)
}

meanVarPlot=function(x,...) {
  myMean=apply(x,1,logMeanMinus)
  myVar=apply(x,1,logVarMinus)
  plot(myMean,myVar)
}

vsubc=function(data,code) {
  return(data[grep(code,names(data))])
}

vminusc=function(data,code) {
  matchCode=names(data)[grep(code,names(data))]
  toIgnore=which(names(data) %in% matchCode)
  if (length(toIgnore)==0) return(data)
  return(data[-toIgnore])
}

plotVln=function(gene,data=dc2,code="rsem",mmax=12,getStat=getStat1,doRet=FALSE,doSort=FALSE) {
  data$GENE=as.character(rownames(data))
  a1=data[gene,c(colnames(data)[grep(code,colnames(data))],"GENE")]
  a2=melt(a1,id="GENE")
  a2$stat=unlist(lapply(as.character(a2$variable),getStat))
  noise <- rnorm(length(a2$value))/100000
  a2$value=a2$value+noise
  if(doSort) {
    a2$stat=factor(a2$stat,levels=names(rev(sort(tapply(a2$value,a2$stat,mean)))))
  }
  p=ggplot(a2,aes(factor(stat),value))
  p2=p + geom_violin(scale="width",adjust=0.75,aes(fill=factor(stat))) + ylab("Expression level (log TPM)")
  p3=p2+theme(legend.position="top")+guides(fill=guide_legend(title=NULL))+geom_jitter(height=0)+blackbg+xlab("Cell Type")
  p4=p3+ theme(axis.title.x = element_text(face="bold", colour="#990000", size=16), axis.text.x  = element_text(angle=90, vjust=0.5, size=12))
  p5=(p4+theme(axis.title.y = element_text(face="bold", colour="#990000", size=16), axis.text.y  = element_text(angle=90, vjust=0.5, size=12))+ggtitle(gene)+theme(plot.title = element_text(size=20, face="bold")))
  if(doRet==TRUE) {
    return(p5)
  }
  else {
    print(p5)
  }
}

getStat1=function(x)return(strsplit(x,"_")[[1]][1])

getStat=function(x,y=1) return(strsplit(x,"_")[[1]][y])
getStat2=function(x,y=2) return(strsplit(x,"_")[[1]][y])
getStat3=function(x,y=3) return(strsplit(x,"_")[[1]][y])


multiplotList <- function(plots, file, cols=1, layout=NULL) {
  require(grid)
  
  # Make a list from the ... arguments and plotlist
  #plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

subc=function(data,code) {
  return(data[,grep(code,colnames(data))])
}

minmax=function(data,min,max) {
  data2=data
  data2[data2>max]=max
  data2[data2<min]=min
  return(data2)
}

arrange <- function(..., nrow=NULL, ncol=NULL, as.table=FALSE) {
  dots <- list(...)
  n <- length(dots)
  if(is.null(nrow) & is.null(ncol)) { nrow = floor(n/2) ; ncol = ceiling(n/nrow)}
  if(is.null(nrow)) { nrow = ceiling(n/ncol)}
  if(is.null(ncol)) { ncol = ceiling(n/nrow)}
  ## NOTE see n2mfrow in grDevices for possible alternative
  grid.newpage()
  pushViewport(viewport(layout=grid.layout(nrow,ncol) ) )
  ii.p <- 1
  for(ii.row in seq(1, nrow)){
    ii.table.row <- ii.row
    if(as.table) {ii.table.row <- nrow - ii.table.row + 1}
    for(ii.col in seq(1, ncol)){
      ii.table <- ii.p
      if(ii.p > n) break
      print(dots[[ii.table]], vp=vp.layout(ii.table.row, ii.col))
      ii.p <- ii.p + 1
    }
  }
}

vp.layout <- function(x, y) viewport(layout.pos.row=x, layout.pos.col=y)
arrange <- function(..., nrow=NULL, ncol=NULL, as.table=FALSE) {
  dots <- list(...)
  n <- length(dots)
  if(is.null(nrow) & is.null(ncol)) { nrow = floor(n/2) ; ncol = ceiling(n/nrow)}
  if(is.null(nrow)) { nrow = ceiling(n/ncol)}
  if(is.null(ncol)) { ncol = ceiling(n/nrow)}
  ## NOTE see n2mfrow in grDevices for possible alternative
  grid.newpage()
  pushViewport(viewport(layout=grid.layout(nrow,ncol) ) )
  ii.p <- 1
  for(ii.row in seq(1, nrow)){
    ii.table.row <- ii.row
    if(as.table) {ii.table.row <- nrow - ii.table.row + 1}
    for(ii.col in seq(1, ncol)){
      ii.table <- ii.p
      if(ii.p > n) break
      print(dots[[ii.table]], vp=vp.layout(ii.table.row, ii.col))
      ii.p <- ii.p + 1
    }
  }
}


subSort = function(vdat, my, fNum, sortBy) {
  vNames=colnames(vdat)
  mycol=which(vNames==sortBy)
  v2 = vdat[order(vdat[,mycol]),]
  vsort=v2
  return(v2)
}

calcMedians = function(p, start, end) {
  medians=c()
  for(i in start:end) {
    scores=p[which(p$factorNum==i),i+6]
    medians=c(medians,mean(scores))
  }
  return(medians)
}

myPalette=
  function (low = "white", high = c("green", "red"), mid = NULL,
            k = 50)
  {
    low <- col2rgb(low)/255
    high <- col2rgb(high)/255
    if (is.null(mid)) {
      r <- seq(low[1], high[1], len = k)
      g <- seq(low[2], high[2], len = k)
      b <- seq(low[3], high[3], len = k)
    }
    if (!is.null(mid)) {
      k2 <- round(k/2)
      mid <- col2rgb(mid)/255
      r <- c(seq(low[1], mid[1], len = k2), seq(mid[1], high[1],
                                                len = k2))
      g <- c(seq(low[2], mid[2], len = k2), seq(mid[2], high[2],
                                                len = k2))
      b <- c(seq(low[3], mid[3], len = k2), seq(mid[3], high[3],
                                                len = k2))
    }
    rgb(r, g, b)
  }

bwCols=myPalette(low = "white",high="black",k = 50)


comparePCA = function(a, b) {
  inds=c(6:26)
  pa=prcomp(a[,inds],scale=TRUE,center=TRUE)
  pb=prcomp(b[,inds],scale=TRUE,center=TRUE)
  print(summary(pa))
  print(summary(pb))
}

getSmooth = function(vsort, myBin, n=0, smooth=0, overlap=0.9, type=0) {
  if (smooth==0) {
    delta=1/(1-overlap)
    nbin=round((length(vsort[,1])-delta)/(n-1+delta))
    smooth=(nbin+1)*delta
    
  }
  return(smooth)
}

condInt = function(vsort, myCtrl, myX, myY, n=0, smooth=0, over=0.9, type=0) {
  overlap=over
  vNames=colnames(vsort)
  myCol=which(vNames==myCtrl)
  xCol=which(vNames==myX)
  yCol=which(vNames==myY)
  binTotal = c()
  nbin=round((length(vsort[,myCol])-smooth)/(n-1)-0.5)
  if (smooth==0) {
    delta=1/(1-overlap)
    nbin=round((length(vsort[,myCol])-delta)/(n-1+delta))
    smooth=(nbin+1)*delta
    
  }
  if (n==0) {
    n=round((length(vsort[,myCol])-smooth)/(smooth*(1-overlap)))
    nbin=round((length(vsort[,myCol])-smooth)/(n-1)-0.5)
  }
  vsort=vsort[order(vsort[,myCtrl]),]
  for (i in 1:n) {
    q1=nbin*(i-1)+1
    q2=q1+smooth
    if (q2 > length(vsort[,1])) {
      q2= length(vsort[,1])
    }
    if (type==0) {
      w=vsort[q1:q2,]
      w=w[order(w[,myX]),]
      x=calcBin(w,myX,n=100,overlap=over)
      y=calcBin(w,myY,n=100,overlap=over)
      val=cor(x,y)
      binTotal=c(binTotal,val)
    }
    if (type==7) {
      w=vsort[q1:q2,]
      w=w[order(w[,myX]),]
      x=calcBin(w,myX,n=100,overlap=over)
      y=calcBin(w,myY,n=100,overlap=over)
      ct=cor.test(x,y)
      width=(ct$conf.int[2]-ct$conf.int[1])/4
      binTotal=c(binTotal,width)
    }
    if (type==5) {
      w=vsort[q1:q2,]
      w=w[order(w[,myX]),]
      x=calcBin(w,myX,n=100)
      y=calcBin(w,myY,n=100)
      y2=calcBin(w,vNames[yCol+1],n=100)
      val=cor(x,y/(y+y2))
      binTotal=c(binTotal,val)
    }
    if (type==6) {
      w=vsort[q1:q2,]
      w=w[order(w[,myX]),]
      x=calcBin(w,myX,n=100)
      y=calcBin(w,myY,n=100,type=2)
      val=cor(x,y)
      binTotal=c(binTotal,val)
    }
  }
  return(binTotal)
}


calcBin = function(vsort, myBin, n=0, smooth=0, overlap=0.9, type=0,cut=0) {
  vNames=colnames(vsort)
  myCol=which(vNames==myBin)
  binTotal = c()
  nbin=round((length(vsort[,myCol])-smooth)/(n-1)-0.5)
  if (smooth==0) {
    delta=1/(1-overlap)
    nbin=round((length(vsort[,myCol])-delta)/(n-1+delta))
    smooth=(nbin+1)*delta
    
  }
  if (n==0) {
    n=round((length(vsort[,myCol])-smooth)/(smooth*(1-overlap)))
    nbin=round((length(vsort[,myCol])-smooth)/(n-1)-0.5)
  }
  
  for (i in 1:n) {
    q1=nbin*(i-1)+1
    q2=q1+smooth
    if (q2 > length(vsort[,1])) {
      q2= length(vsort[,1])
    }
    if (type==0) {
      binTotal=c(binTotal,mean(vsort[q1:q2,myCol]))
    }
    if (type==1) {
      binTotal=c(binTotal,median(vsort[q1:q2,myCol]))
    }
    if (type==2) {
      binTotal=c(binTotal,length(which(vsort[q1:q2,myCol]>cut))/(q2-q1+1))
    }
    if (type==3) {
      binTotal=c(binTotal,sd(vsort[q1:q2,myCol])/sqrt(q2-q1))
    }
    if (type==4) {
      myBin=vsort[q1:q2,]
      myBin=myBin[order(myBin[,myCol]),]
      yv=calcBin(myBin,"CagPkCons",n=100,overlap=0.9)
      xv=calcBin(myBin,vNames[myCol],n=100,overlap=0.9)
      lmod=coef(lm(yv ~ xv))
      slope=cor(xv,yv)
      binTotal=c(binTotal,slope)
    }
    if (type==5) {
      myBin=vsort[q1:q2,]
      myBin=myBin[order(myBin[,myCol]),]
      yv=calcBin(myBin,"CagPkCons",n=100,overlap=0.7)
      xv=calcBin(myBin,vNames[myCol],n=100,overlap=0.7)
      lmod=coef(lm(yv ~ xv))
      ct=cor.test(xv,yv)
      width=(ct$conf.int[2]-ct$conf.int[1])/4
      binTotal=c(binTotal,width)
    }
    
  }
  return(binTotal)
}

genCols=function(al=50) {
  cols=c("darkblue","darkred","darkgreen", "black","orange","purple","khaki","grey","gold4","seagreen3","chocolate")
  tcols=c()
  for (i in 1:6) {
    tcols=c(tcols,rgb (t (col2rgb (cols[i])), alpha = al, maxColorValue = 255))
  }
  return(tcols)
}

plos = function(xvals,yvals,dev,colNum,lwid=4,al=50,ltty=1) {
  #if (length(xvals)>100) {
  #  inds=seq(1,length(xvals),round(length(xvals)/60))+1
  #	xvals=xvals[inds]
  #	yvals=yvals[inds]
  #	dev=dev[inds]
  #}
  cols=c("blue","red","green", "black","orange","purple","khaki","grey","gold4","seagreen3","chocolate")
  tcols=genCols(al)
  lines(xvals,yvals,col=cols[colNum],lwd=lwid,lty=ltty)
  x=c(xvals,rev(xvals))
  y=c(yvals-dev,rev(yvals+dev))
  if (ltty==1) {
    polygon(x,y,col=tcols[colNum],border=NA)
  }
}



getBin = function(vsort, i, n=0, smooth=0,overlap) {
  vNames=colnames(vsort)
  binTotal = c()
  nbin=round((length(vsort[,myCol])-smooth)/(n-1)-0.5)
  if (smooth==0) {
    delta=1/(1-overlap)
    nbin=round((length(vsort[,myCol])-delta)/(n-1+delta))
    smooth=(nbin+1)*delta
    
  }
  if (n==0) {
    n=round((length(vsort[,myCol])-smooth)/(smooth*(1-overlap)))
    nbin=round((length(vsort[,myCol])-smooth)/(n-1)-0.5)
  }
  q1=nbin*(i-1)+1
  q2=q1+smooth
  #print(q1)
  #??print(q2)
  return(vsort[q1:q2,])
}


rHeatMap=function(vdat,ctrlDat,xax,yax,cond,nr=50,nc=50,xwid=0.25,xstep=0.075,xstart=-2,minNum=50,ov=0.9,maxVal=1,ret=1,minVal=0) {
  
}

plosHeatMap=function(vdat,ctrlDat,xax,yax,cond,nr=50,nc=50,xwid=0.25,xstep=0.075,xstart=-2,minNum=50,ov=0.9,maxVal=1,ret=1,minVal=0) {
  
  cNames=colnames(vdat)
  condCol=which(cNames==cond)
  xCol=which(cNames==xax)
  yCol=which(cNames==yax)
  
  ctrl=1:nc
  for(j in 1:nc) {
    start=xstart+xstep*(j-1)
    end=start+xwid
    mySpot= ctrlDat[which(ctrlDat[,xCol]>start),]
    mySpot=mySpot[which(mySpot[,xCol]<=end),]
    ctrl[j]=mean(mySpot[,yCol])
  }
  sum=0
  x=matrix(nrow=nr,ncol=nc,0)
  xDiff=x
  vsort=vdat[order(vdat[,condCol]),]
  numSpots=0
  for(i in 1:nr) {
    myBin=getBin(vsort,i,n=nr,overlap=ov)
    print(paste(i, sum, numSpots, sum/numSpots))
    for (j in 1:nc) {
      start=xstart+xstep*(j-1)
      end=start+xwid
      mySpot=myBin[which(myBin[,xCol]>start),]
      mySpot=mySpot[which(mySpot[,xCol]<=end),]
      myPct=mean(mySpot[,yCol])
      myLen=length(mySpot[,1])
      if (myLen > minNum) {
        myVar=sqrt(myPct*(1-myPct)/myLen)
        x[i,j]=myPct
        obs=myPct*myLen
        
        
        exp=myLen*ctrl[j]
        sum = sum + ((obs-exp)*(obs-exp)/exp)
        xDiff[i,j]=((obs-exp)*(obs-exp)/(exp*myLen))
        #sum = sum + abs(myLen*log(myPct/ctrl[j]))
        numSpots=numSpots+1
        
      }
    }
  }
  x[1,nc]=maxVal
  x[which(x==0)]=1000
  x[which(x<minVal)]=minVal
  x[which(x==1000)]=minVal-0.001
  print(sum)
  print(numSpots)
  if(ret==2) {
    rlist=list(1,2,sum,numSpots)
    rlist[[1]]=x
    rlist[[2]]=xDiff
    return(rlist)
  }
  else {
    return(x)
  }
}


processCutoffs = function(vsort, strongCut, weakCut, colP1, colP2) {
  len=length(vsort[,1])
  hasStrong=rep(0,len)
  hasWeak=hasStrong
  consWeak=hasStrong-1
  fNum = vsort$factorNum
  vNames=colnames(vsort)
  myCol=which(vNames==colP1)
  myCol2=which(vNames==colP2)
  pVals = vsort[,myCol]
  p2Vals=vsort[,myCol2]
  for(i in 1:len) {
    if (((fNum[i] != 5) && (fNum[i] > 0)) && !(is.na(pVals[i]))){
      if (pVals[i]<weakCut[fNum[i]]) {
        hasWeak[i]=1
      }
      if (pVals[i]<strongCut[fNum[i]]) {
        hasWeak[i]=1
        hasStrong[i]=1
      }
      if (hasWeak[i]==1) {
        if (p2Vals[i]<=weakCut[fNum[i]]) {
          consWeak[i]=1
        }
        else {
          consWeak[i]=0
        }
      }
    }
  }
  cutData=cbind(hasStrong,hasWeak, consWeak)
  return(cutData)
}

setupAll = function(strongCut=c(-11.2,-13.5,-11.2,-11,-11,-10),
                    weakCut=c(-9,-11.65,-9.4,-9,-11,-8)) {
  vplus=setup()
  
  vret=cbind(vplus,processCutoffs(vplus,strongCut,weakCut,"PWMP1","PWMP2"),processCutoffs(vplus,strongCut,weakCut,"DistP1","DistP2"))
  return(vret)
}


myCor=function(vdat, toSort, toCalc, n=30,type=1,over=0.9,off=0,ymin=-1,ymax=-1,cut=0) {
  vNames=colnames(vdat)
  mycol=which(vNames==toSort)
  col2=which(vNames==toCalc)
  
  adat=vdat[order(vdat[,mycol]),]
  mdat=calcBin(adat,toCalc,n=n,overlap=over,type=type,cut=cut)
  if (type==2) {
    m1=calcBin(adat,toCalc,n=n,overlap=over)
    m2=calcBin(adat,vNames[col2+1],n=n,overlap=over)
    mdat=(m1/(m1+m2))
  }
  if (type==3) {
    mdat=calcBin(adat,toCalc,n=n,type=2,cut=cut)
  }
  xvals=calcBin(adat,toSort,n=n,overlap=over)
  if(off==0) {
    plot(xvals,mdat,xlab=toSort,ylab=toCalc)
  }
  if(off==1) {
    if (ymax != (-1)) {
      plot(xvals,mdat,xlab="",ylab="",type="n",ylim=c(ymin,ymax))
    }
    else {
      plot(xvals,mdat,xlab="",ylab="",type="n")
    }
  }
  #print(cor.test(xvals,mdat))
  #plot(mdat)
  return(cor(xvals,mdat))
}


plosCor=function(vdat, toSort, toCalc, n=30,col=1,type=0,over=0.9,flip=0,pct=1,xov=0.9,plty=1) {
  vNames=colnames(vdat)
  mycol=which(vNames==toSort)
  col2=which(vNames==toCalc)
  
  adat=vdat[order(vdat[,mycol]),]
  mdat=calcBin(adat,toCalc,n=n,overlap=over,type=type)
  if (flip==1) {
    mdat=1-mdat
  }
  nsmth=getSmooth(adat,toCalc,n=n,overlap=over,type=type)
  if (pct==1) {
    adev=sqrt(mdat*(1-mdat)/nsmth)
  }
  if(pct==0) {
    adev=calcBin(adat,toCalc,n=n,overlap=over,type=3)
  }
  xvals=calcBin(adat,toSort,n=n,overlap=xov)
  plos(xvals,mdat,adev,col,lwid=lwid,ltty=plty)
  #print(cor.test(xvals,mdat))
  #plot(mdat)
  return(cor(xvals,mdat))
}

rd= function(dat,n=2) {
  return(round(dat, n))
}
src= function() {
  source("rahulFxns.R")
  source("RobFig1.txt")
  source("RobFig2.txt")
  #source("RobFig3.txt")
  source("RobFig4.txt")
  #source("RobFig5.txt")
  #source("RobFig6.txt")
  #source("gelShift.R.txt")
}

setupP = function() {
  p1=read.table("pocketOut/summits.1FDR.50",header=TRUE)
  p2=read.table("mappedMotifs/summits.1FDR.9.2.100",header=TRUE)
  p2b=read.table("mappedMotifs/summits.1FDR.11.5.100",header=TRUE)
  p3=read.table("mappedTTM/TAGTM01.summits.1FDR.9.8.500",header=TRUE)
  pca=prcomp(p1[,6:29],scale=TRUE,center=TRUE)
  pcs=-pca$x[,1:5]
  pcmel=-pca$x[,1]
  p=cbind(p1,p2,p2b,p3,pcmel,pcs)
}

setupQ = function() {
  q1=read.table("pocketOut/summits.25ORC.200",header=TRUE)
  q2=read.table("mappedMotifs/summits.25FDR.9.2.100",header=TRUE)
  q2b=read.table("mappedMotifs/summits.25FDR.11.5.100",header=TRUE)
  q2c=read.table("mappedMotifs/summits.25FDR.9.25.100",header=TRUE)
  q2d=read.table("mappedMotifs/summits.25FDR.11.55.100",header=TRUE)
  
  q3=read.table("mappedTTM/TAGTM01.summits.25FDR.9.8.500",header=TRUE)
  q4=read.table("mappedTTM/TAGTM27.summits.25FDR.9.500",header=TRUE)
  
  
  cn=colnames(q1)
  
  tfs=cn[6:29]
  tfs[24]="Z_Stage14"
  tfMeans=c()
  tfVars=c()
  bindingData=matrix(nrow=nrow(q1),ncol=24)
  melRatio=c()
  for(i in 1:24) {
    tmp=subset(q1,factorName==tfs[i])
    tfMeans[i]=mean(tmp[,i+5])
    tfVars[i]=sqrt(var(tmp[,i+5]))
    newm=(tmp[,i+5]-tfMeans[i])
    newsd=tfVars[i]
    melRatio=c(melRatio,newm/newsd);
    bindingData[,i]=(q[,i+5]-tfMeans[i])/tfVars[i]
  }
  pca=prcomp(q1[,6:29],scale=TRUE,center=TRUE)
  pcs=-pca$x[,1:5]
  pcmel=-pca$x[,1]
  
  ap=subset(q1,factorNum<6)
  pca=prcomp(ap[,6:11],scale=TRUE,center=TRUE)
  pcsAP=-pca$x[,1:5]
  pcAP=rep(0,nrow(q1))
  pcAP[1:nrow(ap)]=-pca$x[,1]
  pcAP[(nrow(ap)+1):nrow(q1)]=0
  
  q=cbind(q1,q2,q2b,q2c,q2d,q3,q4,pcmel,pcs,pcAP)
  p=cbind(q,melRatio)
  
}

setupV = function() {
  v1=read.table("summitData/summits.chipSeq",header=TRUE)
  v2=read.table("mappedMotifs/summits.chipSeq.9.2.100",header=TRUE)
  v3=read.table("mappedTTM/TAGTM01.summits.chipSeq.9.8.250",header=TRUE)
  v=cbind(v1,v2,v3)
}

gF1 = function() {
  pdf("1.pdf",height=6,width=8)
  genFig1()
  dev.off()
}

gF2 = function() {
  pdf("2.pdf",height=6,width=8)
  genFig2()
  dev.off()
}

gF3 = function() {
  pdf("33.pdf",height=9,width=5)
  genFig3.3()
  dev.off()
}
gF4 = function() {
  pdf("4.pdf",height=9,width=11)
  genFig4.1()
  dev.off()
}	

gF6 = function() {
  pdf("6.pdf",height=6,width=12)
  genFig6()
  dev.off()
}	


combineIntersection = function (prefix, peakList, toCombine, nDelim=3, returnBinary=0, returnColumn=11, returnReadCount=0,addSuffix=".bed.intersect.bed", promGeneName=0) {
  intersect = read.table(paste(prefix,"/",toCombine,addSuffix,sep=""))
  returnVals=rep(0,nrow(peakList))
  
  uniqueInt = intersect[!duplicated(intersect$V4),]
  spec1IDs = (uniqueInt$V6)
  spec2IDs = as.character(uniqueInt$V10)
  
  if (returnBinary==1) {
    returnVals[spec1IDs]=1
    return(returnVals)
  }
  
  returnVals[spec1IDs]=uniqueInt[,returnColumn]
  return(returnVals)
}

meanmin=function(x,min=1,val=0) {
  if (length(x)>min) {
    return(mean(x))
  } 
  return(val)
}


medmin=function(x,min=1,val=0) {
  if (length(x)>min) {
    return(median(x))
  } 
  return(val)
}

combineIntersection = function (prefix, peakList, toCombine, nDelim=3, returnBinary=0, returnColumn=10, returnReadCount=0,addSuffix=".bed.intersect.bed", promGeneName=0, uniqueCol=4,id1Col=6) {
  intersect = read.table(paste(prefix,"/",toCombine,addSuffix,sep=""))
  returnVals=rep(0,nrow(peakList))
  
  uniqueInt = intersect[!duplicated(intersect[,uniqueCol]),]
  spec1IDs = (uniqueInt[,id1Col])
  spec2IDs = as.character(uniqueInt$V10)
  
  if (returnBinary==1) {
    returnVals[spec1IDs]=1
    return(returnVals)
  }
  
  spec2IDList = unlist(strsplit(spec2IDs,"_"))
  spec2PeakNums=as.numeric(spec2IDList[seq(nDelim,length(spec2IDList),nDelim)])
  
  
  if (promGeneName == 1) {
    geneNames=substr(as.character(spec2IDs),1,10)
    returnVals[spec1PeakNums]=geneNames
    return(returnVals)
  }
  if (returnReadCount > 0) {
    summitFile = read.table(paste("../",toCombine,"_summits.bed",sep=""))
    returnVals[spec1PeakNums]=summitFile$V5[spec2PeakNums]
    return(returnVals)
  }
  
  
  returnVals[spec1IDs]=uniqueInt[,returnColumn]
  return(returnVals)
}

init = function() {
  library(ggplot2)
  
  
  opt <-  opts(legend.title = theme_blank(), # switch off the legend title
               legend.text = theme_text(size=12,face="bold"),        
               legend.key.size = unit(2.5, "lines"),
               legend.key = theme_blank(),
               axis.title.x = theme_text(size = 14, vjust = -0.5),
               axis.title.y = theme_text(size = 14, angle = 90),
               axis.text.x = theme_text(size = 12),
               axis.text.y = theme_text(size = 12),
               plot.title = theme_text(size = 18,vjust=2.5,face="bold"),
               plot.margin = unit(c(2,2,0.75,0.75), "lines"))
  
  lwid=2
}

writ.table=function(a, b) {
  write.table(a,b,quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t")
}


sp=function(data,x,y,z,n=50,min=10,minval=0,maxval=0,qx=0.975,qy=0.975,func="mean") {
  xs=seq(quantile(data[,x],1-qx),quantile(data[,x],qx),length.out=n)
  ys=seq(quantile(data[,y],1-qy),quantile(data[,y],qy),length.out=n)
  t=data
  t$x=as.numeric(factor(cut(t[,x],xs)))
  t$y=as.numeric(factor(cut(t[,y],ys)))
  t=subset(t,!is.na(x)&!is.na(y))
  t$code=paste(t$x,t$y,sep="_")
  res=matrix(0,nrow=n,ncol=n)
  mylev=levels(factor(t$code))
  myres=tapply(t[,z],t$code,meanmin,min)
  if (func=="median") {
    myres=tapply(t[,z],t$code,medmin,min)
  }
  for(i in 1:length(mylev)) {
    inds = as.numeric(unlist(strsplit(mylev[i],"_")))
    res[inds[1],inds[2]]=myres[i]
  }
  mapCols=myPalette(low="white",high="black")
  mapCols[1]="white"
  res[which(res<minval)]=minval
  
  res[which(res==0)]=min(res)-0.01
  if (maxval ==0) {
    maxval=max(res)
  }
  res[which(res>maxval)]=maxval
  n=n-1
  yqs=quantile(t[,y],seq(0,1,1/n))
  xqs=quantile(t[,x],seq(0,1,1/n))
  labRow=as.character(round(xs,2))
  labCol=as.character(round(ys,2))
  #heatmap.2((res),Rowv=NA,Colv=NA,scale="none",col=mapCols,ylab=x,xlab=y,main=z,,,density.info="none",trace="none",keysize=1)
  #,scales=list(x=list(at=1:n, labels=labRow)), x=list(at=1:n, labels=labCol)
  #print(levelplot(res[1:(n-1),1:(n-1)],col.regions=mapCols,xlab=list(label=y,cex=1.5),ylab=list(label=x,cex=1.5),main=list(label=z,cex=2)),scales=list(x=list(at=1:(n-1), labels=labRow[1:(n-1)]), y=list(at=1:(n-1), labels=labCol[1:(n-1)])))
  print(levelplot(res[1:n,1:n],col.regions=mapCols,,xlab=list(label=x,cex=1.3),ylab=list(label=y,cex=1.3),main=list(label=z,cex=2),scales=list(x=list(at=1:n, labels=labRow),y=list(at=1:n, labels=labCol))))
  #levelplot(res)
  return(list(t,res))
}

calcTP = function(cutoff,data,score,real,nTP) {
  return(length(which((data[,score]>cutoff)&(data[,real]>0)))/nTP)
}

calcFP = function(cutoff,data,score,real,nFP) {
  return(length(which((data[,score]>cutoff)&(data[,real]==0)))/nFP)
}

auc=function(data,score,real,n=20) {
  totalPos=length(which(data[,real]==1))
  totalNeg=length(which(data[,real]==0))
  scores=data[,score]
  data$myScore=(scores+min(scores))/(max(scores)+min(scores))
  tp=unlist(lapply(seq(-0.0001,0.9999,1/n),calcTP,data,"myScore",real,totalPos))
  fp=unlist(lapply(seq(-0.0001,0.9999,1/n),calcFP,data,"myScore",real,totalNeg))
  plot(c(fp,1),c(tp,1),xlim=c(0,1),ylim=c(0,1))
  x1=c(1,fp)
  x2=c(1,tp)
  print(sum(diff(rev(x2))*diff(rev(x1)))/2+sum(diff(rev(x1))*rev(x2[-1])))
  return(list(c(1,fp),c(1,tp)))
}

bplot=function(s,CR,q=0.98,save=0) {
  s$curdata=s[,CR]
  org = s[1,"Org"]
  x=qplot(curdata,data=s,geom="density",col=stat,xlim=c(0,quantile(s$curdata,q,na.rm=TRUE)),main=paste(CR,org),xlab="",size=I(1))+opt
  print(x)
  return(x)
}

bplotSave=function(CR) {
  a=bplot(s,CR)
  ggsave(file=paste("pic/",CR,".jpg",sep=""))
  print(paste("pic/",CR,".jpg",sep=""))
}

bcplot=function(s,CR,q=0.98,...) {
  s$curdata=s[,CR]
  qplot(curdata,data=subset(s,cg!="NA"),geom="density",col=stat,xlim=c(0.01,quantile(s$curdata,q),xlab="",size=I(1)),main=CR,...)+facet_grid(. ~ cg )+ opt
}

gea = function(data,file) {
  data[,1]=paste("chr",data[,1],sep="")
  writ.table(data[,1:3],file) 
}


bbplot = function(data,name,start=0.2,doSort=TRUE,...) {
  if (doSort) {
    x=boxplot(split(data[,name],data$stat)[order(tapply(data[,name],data$stat,median,na.rm=TRUE))],...)
    text(1:length(summary(data$stat)),start,summary(data$stat)[order(tapply(data[,name],data$stat,median,na.rm=TRUE))])
  }
  else {
    x=boxplot(split(data[,name],data$stat),...)
    text(1:length(summary(data$stat)),start,summary(data$stat))
  }
}

pyCols=myPalette(low = "magenta",high = "yellow",mid = "black")

rp=function() {par(mfrow=c(1,1))}

calcResidLog=function(x1,y1,mcut=30,toAdd=1) {
  touse=which((x1>mcut)&(y1>mcut))
  x=log(x1+toAdd,2)
  y=log(y1+toAdd,2)
  a=x[touse]
  b=y[touse]
  myLM=lm(b~a)
  myResid=y-predict(myLM,data.frame(a=x))
  return(myResid)
}

logVarDivMean=function(x) return(log(var(exp(x)-1)/mean(exp(x)-1)))

calcResid=function(x1,y1,mcut=30,toAdd=1) {
  touse=which((x1>mcut)&(y1>mcut))
  x=x1
  y=y1
  a=x[touse]
  b=y[touse]
  myLM=lm(b~a)
  myResid=y-predict(myLM,data.frame(a=x))
  return(myResid)
}

gtCut=function(x, cutoff=1) {
  return(length(which(x>cutoff)))
}

findNGene=function(data,is.expr=1) {
  toRet=unlist(lapply(1:ncol(data),function(x)gtCut(data[,x],cutoff=is.expr)))  
  names(toRet)=colnames(data)
  return(toRet)
}

getCoefs=function(data,nbin=20,mycut=1) {
  my_stats=data.frame(data[,1])
  code_humpAvg=apply(data,1,humpMean,min=mycut)
  code_humpAvg[code_humpAvg>9]=9
  code_humpAvg[is.na(code_humpAvg)]=0
  my_stats$code_humpAvg=code_humpAvg
  data[data>mycut]=1
  data[data<mycut]=0
  data$bin=cut(code_humpAvg,nbin)
  data$avg=code_humpAvg
  rownames(my_stats)=rownames(data)
  my_coefs=data.frame(t(sapply(colnames(data[1:(ncol(data)-2)]),getAB,data=data,data2=my_stats,status="code",code2="humpAvg",hasBin=TRUE,doPlot=FALSE)))
  colnames(my_coefs)=c("a","b")
  return(my_coefs)
}

makeScorePlot2=function(allscores,getStatFxn=getStat2,mytitle="Title") {
  alls=data.frame(allscores)
  alls$stat=unlist(lapply(names(allscores),getStatFxn))
  p=ggplot(alls,aes(factor(stat),allscores))
  p2=p + geom_violin(scale="width",adjust=0.75,aes(fill=factor(stat))) + ylab("Expression level (log TPM)")
  p3=p2+theme(legend.position="top")+guides(fill=guide_legend(title=NULL))+geom_jitter(height=0)+blackbg+xlab("Cell Type")
  p4=p3+ theme(axis.title.x = element_text(face="bold", colour="#990000", size=16), axis.text.x  = element_text(angle=90, vjust=0.5, size=12))
  print(p4+theme(axis.title.y = element_text(face="bold", colour="#990000", size=16), axis.text.y  = element_text(angle=90, vjust=0.5, size=12))+ggtitle(mytitle)+theme(plot.title = element_text(size=20, face="bold")))
}

getNewScore=function(mygenes,data=mcell,wt_matrix,myfunc=weighted.mean,code="rsem",scramble=FALSE) {
  my_data=subc(data,code)
  if(scramble) {
    my_data=my_data[,sample(ncol(my_data))]
  }
  myscores=unlist(lapply(colnames(my_data),function(x)myfunc(my_data[mygenes,x],wt_matrix[mygenes,x])))
  names(myscores)=colnames(my_data)
  return(myscores)
}

corCellWeightFast=function(cell1,cell2,wt_matrix,data=subc(cell,"lps_t2"),spear=FALSE) {
  cell1Drop=wt_matrix[rownames(data),cell1]
  cell2Drop=wt_matrix[rownames(data),cell2]
  my_weights=cell1Drop*cell2Drop
  cell1Data=an(data[,cell1])
  cell2Data=an(data[,cell2])
  # print(my_weights)
  my_weights[is.na(my_weights)]=0
  ret=corr(d=as.matrix(cbind(cell1Data,cell2Data)),w=my_weights)
  if (spear==TRUE) {
    ret=corr(d=as.matrix(cbind(rank(cell1Data),rank(cell2Data))),w=my_weights)
  }
  return(ret)
}

covCellWeightFast=function(cell1,cell2,wt_matrix,data=subc(cell,"lps_t2"),spear=FALSE) {
  cell1Drop=wt_matrix[rownames(data),cell1]
  cell2Drop=wt_matrix[rownames(data),cell2]
  my_weights=cell1Drop*cell2Drop
  cell1Data=an(data[,cell1])
  cell2Data=an(data[,cell2])
  
  # print(my_weights)
  ret=wtCov(cell1Data,cell2Data,my_weights)
  if (spear==TRUE) {
    ret=wtCov(rank(cell1Data),rank(cell2Data),my_weights)
  }
  return(ret)
}

scaleSCMatrix2=function(data,wt_matrix,code="rsem") {
  wtX=unlist(lapply(rownames(data),function(x)(sum(data[x,]*wt_matrix[x,])/sum(wt_matrix[x,]))))
  sdX=unlist(lapply(rownames(data),function(x)(wtCov(data[x,],data[x,],wt_matrix[x,]))))
  sData=(data-wtX)/sqrt(sdX)
  return(sData)
}

wtCov=function(x,y,w) {
  w=w/sum(w)
  wtX=sum(x*w)
  wtY=sum(y*w)
  wt_cov=sum(w*(x-wtX)*(y-wtY))
  return(wt_cov)
}
expAlpha=function(mu,coefs) {
  logA=coefs$a
  logB=coefs$b
  return(exp(logA+logB*mu)/(1+(exp(logA+logB*mu))))
}

setWt1=function(x,wts,min=1) {
  wts[x>min]=1
  return(wts)
}

gtCut=function(x, cutoff=1) {
  return(length(which(x>cutoff)))
}

findNGene=function(data,mycut=1) {
  toRet=unlist(lapply(1:ncol(data),function(x)gtCut(data[,x],cutoff=mycut)))  
  names(toRet)=colnames(data)
  return(toRet)
}

setWtMatrix1=function(data,wt_matrix,mycut=1) {
  wt1_matrix=sapply(1:ncol(data),function(x)setWt1(data[,x],wt_matrix[,x],min=mycut))
  colnames(wt1_matrix)=colnames(data)
  return(wt1_matrix)
}


getAB=function(cn="lps_t1_S1_rsem",code="lps_t1",data=cell,data2=cs,code2="avg",status="",ncut=25,hasBin=FALSE,doPlot=FALSE,myfunc=plot,func2=lines,...) {
  if (status == "") {
    status=getStatus(cn)
  }
  if (!(hasBin)) {
    data[data>1]=1
    data[data<1]=0
    data$avg=data2[rownames(data),paste(status,"_",code2,sep="")]
    data[data>9]=9
    data$bin=cut(data$avg,20)
  }
  data$val=data[,cn]
  x1=(tapply(data[,cn],data$bin,mean))
  x2=(tapply(data[,"avg"],data$bin,mean))
  #glm.out = glm(val~avg,data=data,family=binomial)
  glm.out = glm(x1~x2,data=data,family=binomial)
  
  if (doPlot==TRUE) {
    
    pred=(predict(glm.out,data.frame(avg=as.numeric(x2)),type="response"))
    myfunc(x2,x1,pch=16,xlab="Log TPM",ylab="Detection Rate",main=cn,ylim=c(0,1))
    func2(x2,pred,...)
  }
  return(glm.out$coefficients)
}


sensitivityCurve=function(cellName,scData,bulkData,mycut=1,mycex=1,mynew=TRUE,...) {
  cutLocs=cut2(bulkData,g=100,onlycuts=TRUE)
  bulkBin=cut2(bulkData,g=100)
  binaryData=scData[,cellName]
  binaryData[binaryData>mycut]=1
  binaryData[binaryData<mycut]=0
  yBin=tapply(binaryData,bulkBin,mean)
  xBin=tapply(bulkData,bulkBin,mean)
  #glm.out = glm(val~avg,data=data,family=binomial)
  options(warn=-1) #otherwise glm throws an unnecessary error
  glm.out = glm(binaryData~bulkData,family=binomial)
  options(warn=0)  
  x_vals=seq(0,10,0.1)
  y_vals=predict(glm.out,data.frame(bulkData=x_vals),type="response")
  if (mynew) {
    plot(xBin,yBin,pch=16,xlab="Average expression",ylab="Probability of detection",...)
  }
  lines(x_vals,y_vals,lwd=2,...)
}


