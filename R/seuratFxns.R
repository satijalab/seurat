nogrid=theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

sort.column=function(x, col) {
  return(x[order(x[,col]),])
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

#' @export
humpMean=function(x, min=0) {
  return(mean(x[x>min]))
}

#' @export
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

#' @export
subr=function(data,code) {
  return(data[grep(code,rownames(data)),])
}

#' @export
cv=function(x)sd(x)/mean(x)

#' @export
humpCt=function(x, min=0) {
  return(length(x[x>min]))
}

#' @export
log_add=function(x) {
  mpi=max(x)
  return(mpi+log(sum(exp(x-mpi))))
}

#' @export
minusr=function(data,code) {
  matchCode=rownames(data)[grep(code,rownames(data))]
  toIgnore=which(rownames(data) %in% matchCode)
  if (length(toIgnore)==0) return(data)
  toRet=data.frame(data[-toIgnore,])
  rownames(toRet)=rownames(data)[-toIgnore]
  colnames(toRet)=colnames(data)
  return(toRet)
}

#' @export
minusc=function(data,code) {
  matchCode=colnames(data)[grep(code,colnames(data))]
  toIgnore=which(colnames(data) %in% matchCode)
  if (length(toIgnore)==0) return(data)
  return(data[,-toIgnore])
}


#' @export
ainb=function(a,b) {
  a2=a[a%in%b]
  return(a2)
}

#' @export
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

#' @export
no.legend.title=theme(legend.title=element_blank())
#' @export
gg.legend.text=function(x=12,y="bold") return(theme(legend.text = element_text(size = x, face = y)))
#' @export
gg.legend.pts=function(x=6) guides(colour = guide_legend(override.aes = list(size=x)))
#' @export
gg.xax=function(x=16,y="#990000",z="bold",x2=12) return(theme(axis.title.x = element_text(face=z, colour=y, size=x),
                                                              axis.text.x  = element_text(angle=90, vjust=0.5, size=x2)))

#' @export
gg.yax=function(x=16,y="#990000",z="bold",x2=12) return(theme(axis.title.y = element_text(face=z, colour=y, size=x),
                                                              axis.text.y  = element_text(angle=90, vjust=0.5, size=x2)))

#' @export
sub.string=function(x,s1,s2) return(unlist(lapply(x,function(y)gsub(s1,s2,y))))

#' @export
jackStrawF=function(prop=0.1,myR1,myR2=3,data=smD) {
  randGenes=sample(rownames(data),nrow(data)*prop)
  smD.mod=data
  smD.mod[randGenes,]=shuffleMatRow(data[randGenes,])
  fmd.pca=prcomp(smD.mod)
  fmd.x=fmd.pca$x
  fmd.rot=fmd.pca$rotation
  fakeF=unlist(lapply(randGenes,jackF,r1=myR1,r2=myR2,x=fmd.x,rot=fmd.rot))
}

#' @export
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

#' @export
shuffleMatRow=function(x) {
  x2=x
  x2 <- t(x)
  ind <- order(c(col(x2)), runif(length(x2)))
  x2 <- matrix(x2[ind], nrow=nrow(x), ncol=ncol(x), byrow=TRUE)
  return(x2)
}




#' @export
logMeanMinus= function(x)log(mean(exp(as.numeric(x))-1)+1)
#' @export
logVarMinus= function(x)(mean(var(as.numeric(x))-1)+1)

#' @export
logVarMinus2= function(x)(var(exp(as.numeric(x))-1)+1)

quickRNAHuman=function(x) {
  dataFile=paste("~/big/",x,"/summary/",x,".expMatrix.txt",sep="")
  data=log(read.table(dataFile,header=TRUE,sep="\t")[,-1]+1)
  data=subc(data,"rsem")
  return(data)
}

#' @export
expVar=function(x) {
  return(log(var(exp(x)-1)+1))
}

#' @export
expSD=function(x) {
  return(log(sd(exp(x)-1)+1))
}

#' @export
expMean=function(x) {
  return(log(mean(exp(x)-1)+1))
}

#' @export
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

#' @export
cell.centroid=function(cell.probs) {
  centroid.x=round(sum(sapply(1:64,function(x)(x-1)%%8+1)*cell.probs))
  centroid.y=round(sum(sapply(1:64,function(x)(x-1)%/%8+1)*cell.probs))
  centroid.bin=8*(centroid.y-1)+centroid.x
  return(centroid.bin)
}

#' @export
cell.centroid.x=function(cell.probs) {
  return(centroid.x=round(sum(sapply(1:64,function(x)(x-1)%%8+1)*cell.probs)))
}

#' @export
cell.centroid.y=function(cell.probs) {
  return(centroid.y=round(sum(sapply(1:64,function(x)(x-1)%/%8+1)*cell.probs)))
}

#' @export
exact.cell.centroid=function(cell.probs) {
  centroid.x=(sum(sapply(1:64,function(x)(x-1)%%8+1)*cell.probs))
  centroid.y=(sum(sapply(1:64,function(x)(x-1)%/%8+1)*cell.probs))
  return(c(centroid.x,centroid.y))
}

#' @export
marker.auc.test=function(data1,data2,mygenes,print.bar=TRUE) {
  myAUC=unlist(lapply(mygenes,function(x)diffAUC(as.numeric(data1[x,]),as.numeric(data2[x,]))))
  myAUC[is.na(myAUC)]=0
  iterate.fxn=lapply; if (print.bar) iterate.fxn=pblapply
  avg_diff=unlist(iterate.fxn(mygenes,function(x)(expMean(as.numeric(data1[x,]))-expMean(as.numeric(data2[x,])))))
  toRet=data.frame(cbind(myAUC,avg_diff),row.names=mygenes)
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

#' @export
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

TobitDiffExpTest=function(data1,data2,mygenes,print.bar) {
  p_val=unlist(lapply(mygenes,function(x)diffTobit(as.numeric(data1[x,]),as.numeric(data2[x,]))))
  p_val[is.na(p_val)]=1
  iterate.fxn=lapply; if (print.bar) iterate.fxn=pblapply
  toRet=data.frame(p_val,row.names=mygenes)
  return(toRet)
}



NegBinomDiffExpTest=function(data1,data2,mygenes,print.bar) {
  p_val=unlist(lapply(mygenes,function(x)diffLRT(as.numeric(data1[x,]),as.numeric(data2[x,]))))
  p_val[is.na(p_val)]=1
  iterate.fxn=lapply; if (print.bar) iterate.fxn=pblapply
  avg_diff=unlist(iterate.fxn(mygenes,function(x)(expMean(as.numeric(data1[x,]))-expMean(as.numeric(data2[x,])))))
  toRet=data.frame(cbind(p_val,avg_diff),row.names=mygenes)
  toRet=toRet[order(toRet$p_val),]
  return(toRet)
}

diffNegBinom=function(x,data,vars) {
  
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
  return(pchisq(lrt_diff,3,lower.tail = F))
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

meanVarPlot=function(x,...) {
  myMean=apply(x,1,logMeanMinus)
  myVar=apply(x,1,logVarMinus)
  plot(myMean,myVar)
}

#' @export
vsubc=function(data,code) {
  return(data[grep(code,names(data))])
}

#' @export
vminusc=function(data,code) {
  matchCode=names(data)[grep(code,names(data))]
  toIgnore=which(names(data) %in% matchCode)
  if (length(toIgnore)==0) return(data)
  return(data[-toIgnore])
}

#' @export
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

#' @export
extract_field=function(string,field=1,delim="_") {
  fields=as.numeric(unlist(strsplit(as.character(field),",")))
  if (length(fields)==1)  return(strsplit(string,delim)[[1]][field])
  return(paste(strsplit(string,delim)[[1]][fields],collapse = delim))
}

#' @export
getStat1=function(x)return(strsplit(x,"_")[[1]][1])

#' @export
genes.ca.range=function(object,my.min,my.max) {
  ca=ClusterAlpha(object)
  ca.min=apply(ca,1,min)
  ca.max=apply(ca,1,max)
  genes.1=names(ca.min[ca.min<my.max])
  genes.2=names(ca.max[ca.max>my.min])
  return(ainb(genes.1,genes.2))
}


#' @export
getStat=function(x,y=1) return(strsplit(x,"_")[[1]][y])

#' @export
getStat2=function(x,y=2) return(strsplit(x,"_")[[1]][y])

#' @export
getStat3=function(x,y=3) return(strsplit(x,"_")[[1]][y])

#' @export
MultiPlotList <- function(plots, file, cols=1, layout=NULL) {

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

#' @export
subc=function(data,code) {
  return(data[,grep(code,colnames(data))])
}

#' @export
minmax=function(data,min,max) {
  data2=data
  data2[data2>max]=max
  data2[data2<min]=min
  return(data2)
}

#' @export
vp.layout <- function(x, y) viewport(layout.pos.row=x, layout.pos.col=y)

#' @export
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

#' @export
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

#' @export
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


genCols=function(al=50) {
  cols=c("darkblue","darkred","darkgreen", "black","orange","purple","khaki","grey","gold4","seagreen3","chocolate")
  tcols=c()
  for (i in 1:6) {
    tcols=c(tcols,rgb (t (col2rgb (cols[i])), alpha = al, maxColorValue = 255))
  }
  return(tcols)
}

init = function() {
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

#' @export
writ.table=function(a, b) {
  write.table(a,b,quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t")
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

#' @export
pyCols=myPalette(low = "magenta",high = "yellow",mid = "black")

#' @export
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

#' @export
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

#' @export
findNGene=function(data,is.expr=1) {
  toRet=unlist(lapply(1:ncol(data),function(x)gtCut(data[,x],cutoff=is.expr)))
  names(toRet)=colnames(data)
  return(toRet)
}

#' @export
returnTopX=function(data,group.by,n.return,col.return=NA) {
  to.ret=c()
  levels.use=unique(group.by); if (is.factor(group.by)) levels.use=levels(group.by)
  if (!is.na(col.return)) return(unlist(lapply(levels.use, function(x) head(data[group.by==x,col.return],n.return)))) else {
    return(unlist(lapply(levels.use, function(x) head(rownames(data[group.by==x,])))))
  }
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

#' @export
gtCut=function(x, cutoff=1) {
  return(length(which(x>cutoff)))
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

#' @export
heatmap2NoKey=function (x, Rowv = TRUE, Colv = if (symm) "Rowv" else TRUE,
                        distfun = dist, hclustfun = hclust, dendrogram = c("both",
                                                                           "row", "column", "none"), symm = FALSE, scale = c("none",
                                                                                                                             "row", "column"), na.rm = TRUE, revC = identical(Colv,
                                                                                                                                                                              "Rowv"), add.expr, breaks, symbreaks = min(x < 0, na.rm = TRUE) ||
                          scale != "none", col = "heat.colors", colsep, rowsep,
                        sepcolor = "white", sepwidth = c(0.05, 0.05), cellnote, notecex = 1,
                        notecol = "cyan", na.color = par("bg"), trace = c("column",
                                                                          "row", "both", "none"), tracecol = "cyan", hline = median(breaks),
                        vline = median(breaks), linecol = tracecol, margins = c(5,
                                                                                5), ColSideColors, RowSideColors, cexRow = 0.2 + 1/log10(nr),
                        cexCol = 0.2 + 1/log10(nc), labRow = NULL, labCol = NULL,
                        key = TRUE, keysize = 1.5, density.info = c("histogram",
                                                                    "density", "none"), denscol = tracecol, symkey = min(x <
                                                                                                                           0, na.rm = TRUE) || symbreaks, densadj = 0.25, main = NULL,
                        xlab = NULL, ylab = NULL, lmat = NULL, lhei = NULL, axRowCol="black",lwid = NULL, pc = NULL,
                        ...)
{
  scale01 <- function(x, low = min(x), high = max(x)) {
    x <- (x - low)/(high - low)
    x
  }
  retval <- list()
  scale <- if (symm && missing(scale))
    "none"
  else match.arg(scale)
  dendrogram <- match.arg(dendrogram)
  trace <- match.arg(trace)
  density.info <- match.arg(density.info)
  if (length(col) == 1 && is.character(col))
    col <- get(col, mode = "function")
  if (!missing(breaks) && (scale != "none"))
    warning("Using scale=\"row\" or scale=\"column\" when breaks are",
            "specified can produce unpredictable results.", "Please consider using only one or the other.")
  if (is.null(Rowv) || is.na(Rowv))
    Rowv <- FALSE
  if (is.null(Colv) || is.na(Colv))
    Colv <- FALSE
  else if (Colv == "Rowv" && !isTRUE(Rowv))
    Colv <- FALSE
  if (length(di <- dim(x)) != 2 || !is.numeric(x))
    stop("`x' must be a numeric matrix")
  nr <- di[1]
  nc <- di[2]
  if (nr <= 1 || nc <= 1)
    stop("`x' must have at least 2 rows and 2 columns")
  if (!is.numeric(margins) || length(margins) != 2)
    stop("`margins' must be a numeric vector of length 2")
  if (missing(cellnote))
    cellnote <- matrix("", ncol = ncol(x), nrow = nrow(x))
  if (!inherits(Rowv, "dendrogram")) {
    if (((!isTRUE(Rowv)) || (is.null(Rowv))) && (dendrogram %in%
                                                 c("both", "row"))) {
      if (is.logical(Colv) && (Colv))
        dendrogram <- "column"
      else dedrogram <- "none"
      warning("Discrepancy: Rowv is FALSE, while dendrogram is `",
              dendrogram, "'. Omitting row dendogram.")
    }
  }
  if (!inherits(Colv, "dendrogram")) {
    if (((!isTRUE(Colv)) || (is.null(Colv))) && (dendrogram %in%
                                                 c("both", "column"))) {
      if (is.logical(Rowv) && (Rowv))
        dendrogram <- "row"
      else dendrogram <- "none"
      warning("Discrepancy: Colv is FALSE, while dendrogram is `",
              dendrogram, "'. Omitting column dendogram.")
    }
  }
  if (inherits(Rowv, "dendrogram")) {
    ddr <- Rowv
    rowInd <- order.dendrogram(ddr)
  }
  else if (is.integer(Rowv)) {
    hcr <- hclustfun(distfun(x))
    ddr <- as.dendrogram(hcr)
    ddr <- reorder(ddr, Rowv)
    rowInd <- order.dendrogram(ddr)
    if (nr != length(rowInd))
      stop("row dendrogram ordering gave index of wrong length")
  }
  else if (isTRUE(Rowv)) {
    Rowv <- rowMeans(x, na.rm = na.rm)
    hcr <- hclustfun(distfun(x))
    ddr <- as.dendrogram(hcr)
    ddr <- reorder(ddr, Rowv)
    rowInd <- order.dendrogram(ddr)
    if (nr != length(rowInd))
      stop("row dendrogram ordering gave index of wrong length")
  }
  else {
    rowInd <- nr:1
  }
  if (inherits(Colv, "dendrogram")) {
    ddc <- Colv
    colInd <- order.dendrogram(ddc)
  }
  else if (identical(Colv, "Rowv")) {
    if (nr != nc)
      stop("Colv = \"Rowv\" but nrow(x) != ncol(x)")
    if (exists("ddr")) {
      ddc <- ddr
      colInd <- order.dendrogram(ddc)
    }
    else colInd <- rowInd
  }
  else if (is.integer(Colv)) {
    hcc <- hclustfun(distfun(if (symm)
      x
      else t(x)))
    ddc <- as.dendrogram(hcc)
    ddc <- reorder(ddc, Colv)
    colInd <- order.dendrogram(ddc)
    if (nc != length(colInd))
      stop("column dendrogram ordering gave index of wrong length")
  }
  else if (isTRUE(Colv)) {
    Colv <- colMeans(x, na.rm = na.rm)
    hcc <- hclustfun(distfun(if (symm)
      x
      else t(x)))
    ddc <- as.dendrogram(hcc)
    ddc <- reorder(ddc, Colv)
    colInd <- order.dendrogram(ddc)
    if (nc != length(colInd))
      stop("column dendrogram ordering gave index of wrong length")
  }
  else {
    colInd <- 1:nc
  }
  retval$rowInd <- rowInd
  retval$colInd <- colInd
  retval$call <- match.call()
  x <- x[rowInd, colInd]
  x.unscaled <- x
  cellnote <- cellnote[rowInd, colInd]
  if (is.null(labRow))
    labRow <- if (is.null(rownames(x)))
      (1:nr)[rowInd]
  else rownames(x)
  else labRow <- labRow[rowInd]
  if (is.null(labCol))
    labCol <- if (is.null(colnames(x)))
      (1:nc)[colInd]
  else colnames(x)
  else labCol <- labCol[colInd]
  if (scale == "row") {
    retval$rowMeans <- rm <- rowMeans(x, na.rm = na.rm)
    x <- sweep(x, 1, rm)
    retval$rowSDs <- sx <- apply(x, 1, sd, na.rm = na.rm)
    x <- sweep(x, 1, sx, "/")
  }
  else if (scale == "column") {
    retval$colMeans <- rm <- colMeans(x, na.rm = na.rm)
    x <- sweep(x, 2, rm)
    retval$colSDs <- sx <- apply(x, 2, sd, na.rm = na.rm)
    x <- sweep(x, 2, sx, "/")
  }
  if (missing(breaks) || is.null(breaks) || length(breaks) <
      1) {
    if (missing(col) || is.function(col))
      breaks <- 16
    else breaks <- length(col) + 1
  }
  if (length(breaks) == 1) {
    if (!symbreaks)
      breaks <- seq(min(x, na.rm = na.rm), max(x, na.rm = na.rm),
                    length = breaks)
    else {
      extreme <- max(abs(x), na.rm = TRUE)
      breaks <- seq(-extreme, extreme, length = breaks)
    }
  }
  nbr <- length(breaks)
  ncol <- length(breaks) - 1
  if (class(col) == "function")
    col <- col(ncol)
  min.breaks <- min(breaks)
  max.breaks <- max(breaks)
  x[x < min.breaks] <- min.breaks
  x[x > max.breaks] <- max.breaks
  #  if (missing(lhei) || is.null(lhei))
  #    lhei <- c(keysize, 4)
  #  if (missing(lwid) || is.null(lwid))
  #    lwid <- c(keysize, 4)
  #  if (missing(lmat) || is.null(lmat)) {
  #    lmat <- rbind(4:3, 2:1)
  #    if (!missing(ColSideColors)) {
  #      if (!is.character(ColSideColors) || length(ColSideColors) !=
  #        nc)
  #        stop("'ColSideColors' must be a character vector of length ncol(x)")
  #      lmat <- rbind(lmat[1, ] + 1, c(NA, 1), lmat[2, ] +
  #        1)
  #      lhei <- c(lhei[1], 0.2, lhei[2])
  #    }
  #    if (!missing(RowSideColors)) {
  #      if (!is.character(RowSideColors) || length(RowSideColors) !=
  #        nr)
  #        stop("'RowSideColors' must be a character vector of length nrow(x)")
  #      lmat <- cbind(lmat[, 1] + 1, c(rep(NA, nrow(lmat) -
  #        1), 1), lmat[, 2] + 1)
  #      lwid <- c(lwid[1], 0.2, lwid[2])
  #    }
  #    lmat[is.na(lmat)] <- 0
  #  }
  #  if (length(lhei) != nrow(lmat))
  #    stop("lhei must have length = nrow(lmat) = ", nrow(lmat))
  #  if (length(lwid) != ncol(lmat))
  #    stop("lwid must have length = ncol(lmat) =", ncol(lmat))
  #  op <- par(no.readonly = TRUE)
  #  on.exit(par(op))
  #  layout(lmat, widths = lwid, heights = lhei, respect = FALSE)

  if (!missing(RowSideColors)) {
    par(mar = c(margins[1], 0, 0, 0.5))
    image(rbind(1:nr), col = RowSideColors[rowInd], axes = FALSE)
  }
  if (!missing(ColSideColors)) {
    par(mar = c(0.5, 0, 0, margins[2]))
    image(cbind(1:nc), col = ColSideColors[colInd], axes = FALSE)
  }
  oldMar=par()$mar
  if (labCol[1] == ""){
    par(mar = c(margins[1]-3, margins[2]-2, margins[1]-3, margins[2]))
  }
  else{
    par(mar = c(margins[1], margins[2], margins[1], margins[2]))
  }


  x <- t(x)
  cellnote <- t(cellnote)
  if (revC) {
    iy <- nr:1
    if (exists("ddr"))
      ddr <- rev(ddr)
    x <- x[, iy]
    cellnote <- cellnote[, iy]
  }
  else iy <- 1:nr

  # add pc number as title if plotting pc heatmaps
  if(!is.null(pc)){
    pc_title <- paste("PC", pc)
  }
  else{
    pc_title <- ""
  }

  image(1:nc, 1:nr, x, xlim = 0.5 + c(0, nc), ylim = 0.5 +
          c(0, nr), axes = FALSE, xlab = "", ylab = "", main = pc_title , col = col,
        breaks = breaks, ...)
  retval$carpet <- x
  if (exists("ddr"))
    retval$rowDendrogram <- ddr
  if (exists("ddc"))
    retval$colDendrogram <- ddc
  retval$breaks <- breaks
  retval$col <- col
  if (any(is.na(x))) {
    mmat <- ifelse(is.na(x), 1, NA)
    image(1:nc, 1:nr, mmat, axes = FALSE, xlab = "", ylab = "", main = pc_title,
          col = na.color, add = TRUE)
  }
  axis(1, 1:nc, labels = labCol, las = 2, line = -0.5, tick = 0,
       cex.axis = cexCol)
  if (!is.null(xlab))
    mtext(xlab, side = 1, line = margins[1] - 1.25)
  axis(4, iy, labels = labRow, las = 2, line = -0.5, tick = 0,
       cex.axis = cexRow,col=axRowCol)
  if (!is.null(ylab))
    mtext(ylab, side = 4, line = margins[2] - 1.25)
  if (!missing(add.expr))
    eval(substitute(add.expr))
  if (!missing(colsep))
    for (csep in colsep) rect(xleft = csep + 0.5, ybottom = rep(0,
                                                                length(csep)), xright = csep + 0.5 + sepwidth[1],
                              ytop = rep(ncol(x) + 1, csep), lty = 1, lwd = 1,
                              col = sepcolor, border = sepcolor)
  if (!missing(rowsep))
    for (rsep in rowsep) rect(xleft = 0, ybottom = (ncol(x) +
                                                      1 - rsep) - 0.5, xright = nrow(x) + 1, ytop = (ncol(x) +
                                                                                                       1 - rsep) - 0.5 - sepwidth[2], lty = 1, lwd = 1,
                              col = sepcolor, border = sepcolor)
  min.scale <- min(breaks)
  max.scale <- max(breaks)
  x.scaled <- scale01(t(x), min.scale, max.scale)
  if (trace %in% c("both", "column")) {
    retval$vline <- vline
    vline.vals <- scale01(vline, min.scale, max.scale)
    for (i in colInd) {
      if (!is.null(vline)) {
        abline(v = i - 0.5 + vline.vals, col = linecol,
               lty = 2)
      }
      xv <- rep(i, nrow(x.scaled)) + x.scaled[, i] - 0.5
      xv <- c(xv[1], xv)
      yv <- 1:length(xv) - 0.5
      ##lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
    }
  }
  if (trace %in% c("both", "row")) {
    retval$hline <- hline
    hline.vals <- scale01(hline, min.scale, max.scale)
    for (i in rowInd) {
      if (!is.null(hline)) {
        abline(h = i + hline, col = linecol, lty = 2)
      }
      yv <- rep(i, ncol(x.scaled)) + x.scaled[i, ] - 0.5
      yv <- rev(c(yv[1], yv))
      xv <- length(yv):1 - 0.5
      ##lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
    }
  }
  if (!missing(cellnote))
    text(x = c(row(cellnote)), y = c(col(cellnote)), labels = c(cellnote),
         col = notecol, cex = notecex)
  #par(mar = c(margins[1], 0, 0, 0))
  if (dendrogram %in% c("both", "row")) {
    ##plot(ddr, horiz = TRUE, axes = FALSE, yaxs = "i", leaflab = "none")
  }
  ##else plot.new()
  #par(mar = c(0, 0, if (!is.null(main)) 5 else 0, margins[2]))
  if (dendrogram %in% c("both", "column")) {
    ##plot(ddc, axes = FALSE, xaxs = "i", leaflab = "none")
  }
  ##else plot.new()
  key=FALSE
  if (!is.null(main))
    title(main, cex.main = 1.5 * op[["cex.main"]])
  if (key) {
    par(mar = c(5, 4, 2, 1), cex = 0.75)
    tmpbreaks <- breaks
    if (symkey) {
      max.raw <- max(abs(c(x, breaks)), na.rm = TRUE)
      min.raw <- -max.raw
      tmpbreaks[1] <- -max(abs(x), na.rm = TRUE)
      tmpbreaks[length(tmpbreaks)] <- max(abs(x), na.rm = TRUE)
    }
    else {
      min.raw <- min(x, na.rm = TRUE)
      max.raw <- max(x, na.rm = TRUE)
    }
    z <- seq(min.raw, max.raw, length = length(col))
    #image(z = matrix(z, ncol = 1), col = col, breaks = tmpbreaks,
    #      xaxt = "n", yaxt = "n")
    par(usr = c(0, 1, 0, 1))
    lv <- pretty(breaks)
    xv <- scale01(as.numeric(lv), min.raw, max.raw)
    axis(1, at = xv, labels = lv)
    if (scale == "row")
      mtext(side = 1, "Row Z-Score", line = 2)
    else if (scale == "column")
      mtext(side = 1, "Column Z-Score", line = 2)
    else mtext(side = 1, "Value", line = 2)
    if (density.info == "density") {
      dens <- density(x, adjust = densadj, na.rm = TRUE)
      omit <- dens$x < min(breaks) | dens$x > max(breaks)
      dens$x <- dens$x[-omit]
      dens$y <- dens$y[-omit]
      dens$x <- scale01(dens$x, min.raw, max.raw)
      lines(dens$x, dens$y/max(dens$y) * 0.95, col = denscol,
            lwd = 1)
      axis(2, at = pretty(dens$y)/max(dens$y) * 0.95, pretty(dens$y))
      title("Color Key\nand Density Plot")
      par(cex = 0.5)
      mtext(side = 2, "Density", line = 2)
    }
    else if (density.info == "histogram") {
      h <- hist(x, plot = FALSE, breaks = breaks)
      hx <- scale01(breaks, min.raw, max.raw)
      hy <- c(h$counts, h$counts[length(h$counts)])
      lines(hx, hy/max(hy) * 0.95, lwd = 1, type = "s",
            col = denscol)
      axis(2, at = pretty(hy)/max(hy) * 0.95, pretty(hy))
      title("Color Key\nand Histogram")
      par(cex = 0.5)
      mtext(side = 2, "Count", line = 2)
    }
    else title("Color Key")
  }
  ##else plot.new()
  retval$colorTable <- data.frame(low = retval$breaks[-length(retval$breaks)],
                                  high = retval$breaks[-1], color = retval$col)
  invisible(retval)
  par(mar=oldMar)
}

mergeDescendents = function(object, tree, node, pcs, top.genes, acc.cutoff){
  # find daughter cells of given node in given tree
  daughters = tree$edge[which(tree$edge[,1]==node),2]
  # get the children of both daughters
  childNodes = 1:(tree$Nnode+1)
  if(length(ainb(c(daughters[1], daughters[2]), childNodes))==2){
    d1 = WhichCells(object,daughters[1])
    d2 = WhichCells(object,daughters[2])
    y  = as.numeric(object@ident[c(d1,d2)])-1
    x  = data.frame(t(object@data[PCTopGenes(object,pcs,num.genes = top.genes ),c(d1,d2)]));
    xv = apply(x,2,var)
    x  = x[,names(xv>0)]
    # run k-fold cross validation
    ctrl = trainControl(method = "repeatedcv", repeats = 5)
    set.seed(1500)
    model = train(as.factor(y)~., data=x, method = "svmLinear", trControl = ctrl)
    acc = model$results[,2]
    # if classifier can't classify them well enough, merge clusters
    if(acc<acc.cutoff){
      object = SetIdent(object,cells.use = WhichCells(object,daughters[1]), ident.use = daughters[2])
    }
    return(object)
  }
  # recursion to traverse tree
  if(daughters[1]%in%childNodes){
    object = mergeDescendents(object, tree, daughters[2], pcs, top.genes, acc.cutoff)
    return(object)
  }
  if(daughters[2]%in%childNodes){
    object = mergeDescendents(object, tree, daughters[1], pcs, top.genes, acc.cutoff)
    return(object)
  }
  object = mergeDescendents(object, tree, daughters[1], pcs, top.genes, acc.cutoff)
  object = mergeDescendents(object, tree, daughters[2], pcs, top.genes, acc.cutoff)
  return(object)
}

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

same=function(x) return(x)

DFT <- function(tree, node, path = NULL, include.children = F, only.children = F){
  children <- which(tree$edge[,1] == node)
  child1 <- tree$edge[children[1], 2]
  child2 <- tree$edge[children[2], 2]
  if(child1 %in% tree$edge[,1]){
    if(!only.children){
      path <- c(path, child1)
    }
    path <- DFT(tree, child1, path, include.children)
  }
  else{
    if(include.children){
      path <-c(path, child1)
    }
  }
  if(child2 %in% tree$edge[,1]){
    if(!only.children){
      path <- c(path, child2)
    }
    path <- DFT(tree, child2, path, include.children)
  }
  else{
    if(include.children){
      path <-c(path, child2)
    }
  }
  return(path)
}
