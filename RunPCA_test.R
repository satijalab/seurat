devtools::install('.') ##IMPORTANT: otherwise Eigen will not be compiled correctly
library(SeuratData)
library(Seurat)

datasets <- list('bmcite','pbmc3k','cbmc',
                 'hcabm40k','pbmcsca','ifnb')
old_time <- rep(NA,length(datasets))
new_time <- rep(NA,length(datasets))
diffs <- rep(NA,length(datasets))
names(old_time) <- names(new_time) <- names(diffs) <- datasets
for (data in datasets) {
  obj <- LoadData(data)
  obj <- NormalizeData(obj)
  obj <- ScaleData_fast(obj,scale.max=Inf) #required for consistency in original RunPCA btwn approx=T vs F
  obj <- FindVariableFeatures(obj)
  
  st <- Sys.time()
  obj1 <- RunPCA(obj,approx=T)
  et <- Sys.time() - st
  old_time[data] <- et
  
  obj2 <- RunPCA(obj,approx=F)
  
  st <- Sys.time()
  obj3 <- RunPCA_fast(obj)
  et <- Sys.time() - st
  new_time[data] <- et
  
  emb2 <- Embeddings(obj2, reduction = 'pca')
  emb3 <- Embeddings(obj3, reduction = 'pca')
  npc <- min(ncol(emb2), ncol(emb3))
  emb2 <- emb2[, seq_len(npc), drop = FALSE]
  emb3 <- emb3[, seq_len(npc), drop = FALSE]
  s <- sign(colSums(emb2 * emb3))
  s[s == 0] <- 1
  emb3 <- sweep(emb3, 2, s, '*')
  diffs[data] <- max(abs(emb2 - emb3))
}

old_time
new_time
diffs
