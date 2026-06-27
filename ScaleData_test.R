devtools::install('.') ##must do this instead of load_all('.')
library(SeuratData)
library(Seurat)

datasets <- list('bmcite','pbmc3k','cbmc',
                 'hcabm40k','pbmcsca','ifnb')
old_time <- rep(NA,length(datasets))
new_time <- rep(NA,length(datasets))
new_time_thread <- rep(NA,length(datasets))
diffs <- rep(NA,length(datasets))
diffs_thread <- rep(NA,length(datasets))
names(old_time) <- names(new_time) <- names(diffs) <- 
  names(new_time_thread) <- names(diffs_thread) <- datasets
for (data in datasets) {
  obj <- LoadData(data)
  obj <- NormalizeData(obj)
  
  st <- Sys.time()
  obj1 <- ScaleData(obj)
  et <- Sys.time() - st
  old_time[data] <- et
  
  st <- Sys.time()
  obj2 <- ScaleData_fast(obj)
  et <- Sys.time() - st
  new_time[data] <- et
  
  st <- Sys.time()
  obj3 <- ScaleData_fast(obj,nthread=4)
  et <- Sys.time() - st
  new_time_thread[data] <- et
  
  diffs[data] <- max(abs(GetAssayData(obj1,layer='scale.data')-
    GetAssayData(obj2,layer='scale.data')))
  diffs_thread[data] <- max(abs(GetAssayData(obj1,layer='scale.data')-
                           GetAssayData(obj3,layer='scale.data')))
}

old_time
new_time
old_time/new_time
diffs
