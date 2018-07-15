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
  argnames <- names(formals(sys.function(sys.parent(n = 1))))
  
  argnames <- grep("object",argnames,invert = T,value = T)
  argnames <- grep("\\.\\.\\.",argnames,invert = T,value = T)
  params <- list()

  #check if function works on the Assay Level
  assay.flag="assay.use"%in%argnames
  p.env <- parent.frame(1)
  argnames <- intersect(argnames, ls(p.env))
  #fill in params list
  for(i in argnames) {
    param_value <- get(x = i,envir = p.env) 
    #TODO Institute some check of object size?
    params[[i]] <- param_value
  }
  #rename function name to include assay info (if needed)
  if(assay.flag) command.name <- paste0(command.name,".",params[["assay.use"]])
  #Store results
  seurat.command <- new(
    Class = 'SeuratCommand',
    name = command.name,
    params = params,
    time.stamp=time.stamp,
    call.string=call.string
  )
  object[[command.name]] <- seurat.command
  return(object)
}


