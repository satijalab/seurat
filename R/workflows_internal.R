# Validates the workflow file that was provided.
#
# @param config results of reading in the ini file
#
# @return No return
#
ValidateWorkflowFile <- function(config){
 sections <- names(config)
 if (! "dependencies" %in% sections) {
   stop("Workflow file is missing the dependencies section")
 }
}


# Read Seurat workflow parameters
#
# Reads parameters for a workflow and assigns them to the parent function call.
#
# @param depth is the depth of the call in the function stack. If called
# directly from, for example, ScaleData, then depth=1. If called indirectly
# from PrepareWorkflow, then depth=2
#
ReadWorkflowParams <- function(object, workflow.name, depth = 2) {
  CheckWorkflow(object = object, workflow.name = workflow.name)
  param.list <- names(formals(fun = sys.function(sys.parent(depth))))
  workflow.params <- slot(object = object[[workflow.name]], name = "params")

  # global variables
  to.set <- intersect(x = param.list, y = names(workflow.params$global))
  p.env <- parent.frame(depth)
  for(i in to.set) {
    assign(x = names(workflow.params$global[i]),
           value = workflow.params$global[[i]],
           envir = p.env)
  }

  # parameter-specific variables
  command.name <- as.character(deparse(sys.calls()[[sys.nframe()-depth]]))
  command.name <- gsub(pattern = ".Seurat",replacement = "",x = command.name)
  command.name <- ExtractField(string = command.name, field = 1, delim = "\\(")
  to.set <- intersect(x = param.list, y = names(workflow.params[[command.name]]))
  for(i in to.set) {
    assign(x = names(workflow.params[[command.name]][i]),
           value = workflow.params[[command.name]][[i]],
           envir = p.env)
  }

  # overwrite any arguments passed in on the command line
  argnames <- sys.call(which = depth)
  argList <- as.list(argnames[-1])
  args_ignore <- c("", "object", "workflow.name")
  args_use <- setdiff(x = names(argList), y = args_ignore)
  for(i in args_use) {
    if(as.character(unlist(argList[i])[[1]]) == "F") {
      arg.val <- FALSE
    } else if(as.character(unlist(argList[i])[[1]]) == "T") {
      arg.val <- TRUE
    } else {
      arg.val <- unlist(argList[i])[[1]]
    }
    assign(x = i, value = arg.val, envir = p.env)
  }
}

# Prepare a Seurat function for a workflow run
#
# Checks dependencies (and runs them if necessary), and then reads in parameter values
#
PrepareWorkflow <- function(object, workflow.name) {
  CheckWorkflow(object = object, workflow.name = workflow.name)
  command.name <- as.character(deparse(sys.calls()[[sys.nframe()-1]]))
  command.name <- gsub(pattern = ".Seurat",replacement = "",x = command.name)
  command.name <- ExtractField(string = command.name, field = 1, delim = "\\(")
  depends <- slot(object = object[[workflow.name]], name = "depends")
  prereq.commands <- colnames(depends)[which(depends[command.name, ] == 1)]
  for(i in prereq.commands) {
    check.prereqs <- CheckWorkflowUpdate(
      object = object,
      workflow.name = workflow.name,
      command.name = i)
    if (check.prereqs) {
      # run the dependency
      workflow.name.quotes <- paste0('\"', workflow.name, '\"')
      new.cmd <- paste0("object <- ", i, "(object, workflow = ", workflow.name.quotes, ")")
      message(paste0("Updating ", i))
      message(new.cmd)
      eval(expr = parse(text = new.cmd))
    }
  }
  ReadWorkflowParams( object = object, workflow.name = workflow.name, depth = 2)
  return(object)
}

# UpdateWorkflow
#
# Updates a workflow object after a command is run, and makes sure timestamps are properly set.
#
UpdateWorkflow <- function(object, workflow.name, command.name = NULL) {
  CheckWorkflow(object = object,workflow.name = workflow.name)
  if(is.null(command.name)) {
    command.name <- as.character(deparse(sys.calls()[[sys.nframe()-1]]))
    command.name <- gsub(pattern = ".Seurat",replacement = "",x = command.name)
    command.name <- ExtractField(string = command.name,field = 1,delim = "\\(")
    command.name.seurat <- intersect(c(command.name, paste0(command.name,".",DefaultAssay(object))), names(object))
  } else {
    command.name.seurat <- command.name
    command.name <- ExtractField(string = command.name, field = 1, delim = "\\.")
  }

  #TODO - Deal with Assay better
  seurat.timestamp <- Sys.time()
  if (length(x = command.name)==1) {
    seurat.timestamp <- slot(object = object[[command.name.seurat]], name = "time.stamp")
  }
  object <- TouchWorkflow(object = object,workflow.name = workflow.name, command.name = command.name,time.stamp = seurat.timestamp)
  return(object)
}
