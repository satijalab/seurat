# Logs a command run, storing the name, timestamp, and argument list. Stores in the Seurat object
# object is the name of the object
# not sure if this should be internal (I think it should be)
LogSeuratCommand <- function(object) {
  call.string <- deparse(sys.call(which = 1))
  time.stamp <- Sys.time()
  
  #capture function name
  command.name <- as.character(deparse(sys.calls()[[sys.nframe()-1]]))
  command.name <- gsub(pattern = ".Seurat",replacement = "",x = command.name)
  command.name <- ExtractField(string = command.name,field = 1,delim = "\\(")
  #capture function arguments
  arg_list <- formals(sys.function(sys.parent(n = 1)))
  
  #check if function works on the Assay Level
  assay.flag=FALSE
  
  #set any default arguments
  for(i in names(x = arg_list)) {
    #don't want to store this (and its usually not passed)
    if(i=="object") {
      arg_list[[i]] <- NULL
    }
    if(i=="assay.use") {
      assay.flag <- TRUE
      if (is.null(arg_list[[i]])) {
        arg_list[[i]] <- DefaultAssay(object = object)
      }
    }
  }
  #rename function name to include assay info (if needed)
  if(assay.flag) command.name <- paste0(command.name,".",arg_list[["assay.use"]])
  
  #Store results
  seurat.command <- new(
    Class = 'SeuratCommand',
    name = command.name,
    params = arg_list,
    time.stamp=time.stamp,
    call.string=call.string
  )
  object[[command.name]] <- seurat.command
  return(object)
}


