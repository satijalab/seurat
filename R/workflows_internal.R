
#' Read Seurat workflow parameters
#'
#' Reads parameters for a workflow and assigns them to the parent function call. No return
#' depth is the depth of the call in the function stack. If called directly from, for example, ScaleData, then depth=1.
#' If called indirectly from PrepareWorkflow, then depth=2
#'
ReadWorkflowParams <- function(object, workflow.name,depth=2) {
  CheckWorkflow(object = object,workflow.name = workflow.name)
  param.list <- names(formals(fun = sys.function(sys.parent(depth))))
  workflow.params <- (object[[workflow.name]]@params)
  exact.match <- intersect(param.list,names(workflow.params))
  generic.params <- grep(pattern = "\\.$",x = names(workflow.params),value = T)
  p.env <- parent.frame(depth)
  for(i in exact.match) {
    assign(x = i,value = unlist(workflow.params[i])[[1]],envir = p.env)
  }
  for(i in generic.params) {
    generic.match <- grep(pattern = i,x = param.list,value = T)
    for(j in generic.match) {
      assign(x = j,value = unlist(workflow.params[i])[[1]],envir = p.env)
    }
  }
}

#' Prepare a Seurat function for a workflow run
#'
#' Checks dependencies (and runs them if necessary), and then reads in parameter values
#'
PrepareWorkflow <- function(object, workflow.name) {
  CheckWorkflow(object = object,workflow.name = workflow.name)
  command.name <- as.character(deparse(sys.calls()[[sys.nframe()-1]]))
  command.name <- gsub(pattern = ".Seurat",replacement = "",x = command.name)
  command.name <- ExtractField(string = command.name,field = 1,delim = "\\(")
  depends <- slot(object = object[[workflow.name]],name = "depends")
  prereq.commands <- colnames(depends)[which(depends[command.name,]==1)]
  for(i in prereq.commands) {
    check.prereqs <- CheckWorkflowUpdate(object = object,workflow.name = workflow.name,command.name = i)
    if (check.prereqs) {
      # run the dependency
      workflow.name.quotes <- paste0('\"',workflow.name,'\"')
      new.cmd <- paste0("object <- ", i, "(object, workflow = ", workflow.name.quotes, ")")
      print(paste0("Updating ", i))
      print(new.cmd)
      eval(expr = parse(text = new.cmd))
    }
  }
  ReadWorkflowParams(object = object,workflow.name = workflow.name,depth = 2)
}

#' UpdateWorkflow
#'
#' Updates a workflow object after a command is run, and makes sure timestamps are properly set. 
#'
UpdateWorkflow <- function(object, workflow.name) {
  CheckWorkflow(object = object,workflow.name = workflow.name)
  command.name <- as.character(deparse(sys.calls()[[sys.nframe()-1]]))
  command.name <- gsub(pattern = ".Seurat",replacement = "",x = command.name)
  command.name <- ExtractField(string = command.name,field = 1,delim = "\\(")
  object <- TouchWorkflow(object = object,workflow.name = workflow.name, command.name = command.name)
  return(object)
}


