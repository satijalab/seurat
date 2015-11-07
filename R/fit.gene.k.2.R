#' Build mixture models of gene expression
#'
#' Proof-of-concept to accelarate runing speed of fit.gene.k function
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
setGeneric("fit.gene.k.fast", function(object, gene, do.k=2,num.iter=1,do.plot=FALSE,genes.use=NULL,start.pct=NULL) standardGeneric("fit.gene.k.fast"))
#' @export
setMethod("fit.gene.k.fast", "seurat",
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
                      cell.ident=iter.k.fit.fast(scale.data,cell.ident,data.use)
                      ident.table=table(cell.ident)
                  }
              }
              # ident.table=table(cell.ident)
              # raw.probs=t(sapply(data.use,function(y) unlist(lapply(1:do.k,function(x) ((ident.table[x]/sum(ident.table))*dnorm(y,mean(as.numeric(data.use[cell.ident==x])),sd(as.numeric(data.use[cell.ident==x]))))))))
              # norm.probs=raw.probs/apply(raw.probs,1,sum)
              # colnames(norm.probs)=unlist(lapply(1:do.k,function(x)paste(gene,x-1,"post",sep=".")))
              # norm.probs=cbind(norm.probs,cell.ident); colnames(norm.probs)[ncol(norm.probs)]=paste(gene,".ident",sep="")
              data.use.t <- as.data.frame(t(data.use))
              data.use.t <- cbind(data.use.t, cell_ident = cell.ident)
              data.use.t.melt <- melt(data.use.t, id.vars = "cell_ident")
              kmodal.mu <- dcast(data.use.t.melt, cell_ident ~ variable, mean)
              kmodal.sd <- dcast(data.use.t.melt, cell_ident ~ variable, sd)
              kmodal.norm_factor <- as.numeric(ident.table / sum(ident.table))
              raw.probs <- sapply(1:do.k, function(k) {
                                            factor.k <- (kmodal.norm_factor[k])
                                            mean.k <- kmodal.mu[k, gene]
                                            sd.k <- kmodal.sd[k, gene]
                                            prob.k <- (dnorm(data.use.t[, gene], mean = mean.k, sd = sd.k))
                                            return(prob.k)
                                            })
              norm.probs <- as.data.frame(raw.probs / rowSums(raw.probs))
              #colnames(norm.probs)=unlist(lapply(1:do.k,function(x)paste(gene,x-1,"post",sep=".")))
              colnames(norm.probs) <- paste(gene, (1:do.k)-1, "post", sep = ".")
              row.names(norm.probs) <- row.names(data.use.t)
              norm.probs <- cbind(norm.probs, cell.ident)
              colnames(norm.probs)[-1] <- paste0(gene, ".ident")

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

calc.dist <- function(x, v, method = c('euclidean')) {
    ## x: matrix. F x N
    ## v: vector. 1 x F
    if (method == 'euclidean') {
        return(sqrt(colSums((x - v) ^ 2)))
    }
}

iter.k.fit.fast <- function(scale.data, cell.ident, data.use) {
    cell.ident.K <- sort(unique(cell.ident))
    means.all <- sapply(cell.ident.K, function(x) 
            rowMeans(scale.data[, cell.ident == x]))
    all.dist <- lapply(cell.ident.K, function(i) 
            calc.dist(scale.data, means.all[, i]))
    all.dist <- matrix(unlist(all.dist), ncol = length(cell.ident.K))
    cell.ident <- apply(all.dist, 1, which.min)
    cell.ident <- order(tapply(as.numeric(data.use), cell.ident, mean))[cell.ident] 
    return(cell.ident)
}