devtools::load_all('.')
library(SeuratData)

datasets <- list('bmcite','pbmc3k','cbmc',
                 'hcabm40k','pbmcsca','ifnb')
old_time <- rep(NA,length(datasets))
new_time <- rep(NA,length(datasets))
diffs <- rep(NA,length(datasets))
names(old_time) <- names(new_time) <- names(diffs) <- datasets
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
  
  diffs[data] <- max(abs(GetAssayData(obj1,layer='scale.data')-
    GetAssayData(obj2,layer='scale.data')))
}

old_time
new_time
old_time[c('hcabm40k','pbmcsca')] <- 60*old_time[c('hcabm40k','pbmcsca')] # in seconds
old_time/new_time
diffs
