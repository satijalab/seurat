# Logs a command run, storing the name, timestamp, and argument list. Stores in
# the Seurat object
# @param object Name of Seurat object
#
# @return returns the Seurat object with command stored
#
LogSeuratCommand <- function(object) {
  time.stamp <- Sys.time()
  #capture function name
  command.name <- as.character(deparse(sys.calls()[[sys.nframe()-1]]))
  command.name <- gsub(pattern = ".Seurat",replacement = "",x = command.name)
  call.string <- command.name
  command.name <- ExtractField(string = command.name,field = 1,delim = "\\(")
  
  #capture function arguments
  argnames <- names(formals(sys.function(sys.parent(n = 1))))
  argnames <- grep(pattern = "object", x = argnames, invert = TRUE, value = TRUE)
  argnames <- grep(pattern = "\\.\\.\\.", x = argnames, invert = TRUE, value = TRUE)
  params <- list()
  p.env <- parent.frame(1)
  argnames <- intersect(argnames, ls(p.env))

  # fill in params list
  for(arg in argnames) {
    param_value <- get(x = arg, envir = p.env)
    #TODO Institute some check of object size?
    params[[arg]] <- param_value
  }
  # check if function works on the Assay and/or the DimReduc Level
  assay.use <- params[["assay.use"]]
  reduction.use <- params[["reduction.use"]]

  # rename function name to include Assay/DimReduc info
  command.name <- paste0(command.name, ".", assay.use, ".", reduction.use)
  command.name <- sub(pattern = "[\\.]+$", replacement = "", x = command.name, perl = TRUE)
  command.name <- sub(pattern = "\\.\\.", replacement = "\\.", x = command.name, perl = TRUE)

  # store results
  seurat.command <- new(
    Class = 'SeuratCommand',
    name = command.name,
    params = params,
    time.stamp = time.stamp,
    call.string = call.string
  )
  object[[command.name]] <- seurat.command
  return(object)
}


