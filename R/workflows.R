#' Initialize a Seurat workflow
#'
#' Reads dependencies from a file and initializes a Seurat workflow
#'
#' @param object Seurat object
#' @param file Ini configuration file. See cluster.workflow.ini for an example.
#'
#' @return Object with modified workflows
#'
#' @importFrom ini read.ini
#' @export
#'
#' @examples
#' pbmc_small <- InitializeWorkflow(object = pbmc_small, file = 'workflows/cluster.workflow.txt')
#'
InitializeWorkflow <- function(object, file) {
  if(!file.exists(... = file)) {
    stop("Provided workflow file does not exist.")
  }
  config <- read.ini(filepath = file)
  ValidateWorkflowFile(config = config)
  workflow.name <- gsub(
    pattern = ".workflow.ini",
    replacement = "",
    x = basename(path = file))

  depend.fxns <- unlist(strsplit(x = unname(unlist(config$dependencies)), split = ","))
  fxns <- union(depend.fxns, names(config$dependencies))
  depends <- matrix(nrow = length(x = fxns), ncol = length(x = fxns))
  rownames(depends) <- colnames(depends) <- fxns
  for(cmd in 1:length(config$dependencies)) {
    cmd.name <- names(config$dependencies[cmd])
    cmd.vals <- unlist(strsplit(x = config$dependencies[[cmd]], split = ","))
    for(cv in cmd.vals){
      depends[cmd.name, cv] <- 1
    }
  }
  mostRecent <- rep(x = as.POSIXct("1900-01-01"), length(fxns))
  names(mostRecent) <- fxns

  for(mr in names(mostRecent)) {
    assay.use <- config$global$assay.use
    assay.use <- assay.use %iff% config[mr]$assay.use
    assay.use <- assay.use %||% DefaultAssay(object = object)
    reduction.use <- config$global$reduction.use
    reduction.use <- config[mr]$reduction.use
    reduction.use <- reduction.use %||% formals(fun = paste0(mr, ".Seurat"))$reduction.use
    command.name <- paste0(mr, ".", assay.use, ".", reduction.use)
    command.name <- sub(pattern = "[\\.]+$", replacement = "", x = command.name, perl = TRUE)
    command.name <- sub(pattern = "\\.\\.", replacement = "\\.", x = command.name, perl = TRUE)
    if(command.name %in% names(x = object)) {
      seurat.timestamp <- slot(object = object[[command.name]], name = "time.stamp")
      mostRecent[mr] <- seurat.timestamp
    }
  }
  params <- list()
  if (!is.null(config$global)){
    params[["global"]] <- config$global
    for(p in 1:length(params$global)){
      params$global[names(params$global[p])] <- ToNumeric(x = params$global[[p]])
    }
  }
  # set fxn specific params
  fxn.param.names <- setdiff(names(config), c("dependencies", "global"))
  if(length(x = fxn.param.names) > 0) {
    for(i in 1:length(fxn.param.names)) {
      params[fxn.param.names[i]] <- config[fxn.param.names[i]]
      for(p in 1:length(params[[fxn.param.names[i]]])){
        params[[fxn.param.names[i]]][[p]] <- ToNumeric(x = params[[fxn.param.names[i]]][[p]])
      }
    }
  }
  seurat.workflow <- new(
    Class = 'SeuratWorkflow',
    name = workflow.name,
    depends = depends,
    params = params,
    mostRecent = mostRecent
  )
  object[[workflow.name]] <- seurat.workflow
  return(object)
}

#' Checks if a workflow is defined for a Seurat object
#'
#' Checks if a workflow is defined for a Seurat object
#'
#' @param object Seurat object
#' @param workflow.name Workflow name, should already be initialized using InitializeWorkflow#'
#' @return TRUE if workflow is defined. STOP otherwise
#'
#' @export
#'
#' @examples
#' CheckWorkflow(pbmc_small, "cluster")
#'
CheckWorkflow <- function(object, workflow.name) {
  # Check if workflow is there
  workflow.present <- FALSE
  if (workflow.name%in%names(object)) {
    if (class(x = object[[workflow.name]])[[1]] == "SeuratWorkflow") {
      workflow.present <- TRUE
    }
  }
  if (!workflow.present) {
    stop("Workflow not present, initialize first.")
  }
  return(TRUE)
}

#' Set Seurat workflow parameters
#'
#' Sets parameters for a workflow
#'
#' @param object Seurat object
#' @param workflow.name Workflow name, should already be initialized using
#' InitializeWorkflow
#' @param fxn.param.names Name of the function and parameter to set (formatted as
#' FXNNAME_PARAM). Can take a vector to set multiple functions
#' @param fxn.param.values Value of the parameter to set. Should be of equal length
#' to fxn.param.names
#' @param \dots Global parameters to set, will be fed into workflow Seurat functions.
#' Parameters ending with a "." will populate all similar variable names (i.e.
#'  setting dims. will set both dims.compute, and dims.cluster)
#'
#' @return Object with modified workflows
#'
#' @export
#'
#' @examples
#' pbmc_small <- SetWorkflowParams(object = pbmc_small, seed.use = 31, dims. = 20)
#'
SetWorkflowParams <- function(
  object,
  workflow.name = NULL,
  ...,
  fxn.param.names = NULL,
  fxn.param.values = NULL
  ) {
  CheckWorkflow(object = object, workflow.name = workflow.name)
  if(length(x = fxn.param.names) != length(x = fxn.param.values)) {
    stop("length of fxn.parameter.names needs to equal length of fxn.parameter.values")
  }
  params <- slot(object = object[[workflow.name]], name = "params")
  # set global params
  global.params <- list(...)
  for(gp in 1:length(global.params)) {
    params[["global"]][names(global.params[gp])] <- global.params[gp]
  }
  # set fxn specific params
  if(!is.null(fxn.param.names)) {
    for(i in 1:length(fxn.param.names)) {
      fxn <- unlist(strsplit(x = fxn.param.names[i], split = "_"))
      params[[fxn[1]]][fxn[2]] <- fxn.param.values[i]
    }
  }
  slot(object = object[[workflow.name]], name = "params") <- params
  return(object)
}

#' Check if workflow command needs update
#'
#' Compares the stored timestamp with the most recently recorded timestamp to see if a dependency has been updated
#'
#' @param object Seurat object
#' @param workflow.name Workflow name, should already be initialized using InitializeWorkflow
#' @param command.name Name of the command to check
#'
#' @return Returns TRUE if the dependency has changed (or has not been run), and an update is needed. FALSE otherwise
#'
#' @export
#'
#' @examples
#' CheckWorkflowUpdate(object = pbmc_small,workflow.name = "cluster", command.name = "ScaleData")
#' #'
CheckWorkflowUpdate <- function(object, workflow.name, command.name) {
  CheckWorkflow(object = object, workflow.name = workflow.name)

  # According to the workflow, the most recent update
  mostRecent <- slot(object[[workflow.name]],"mostRecent")
  workflow.timestamp <- mostRecent[command.name]
  seurat.timestamp <- as.POSIXct("1900-01-01");

  #means Seurat command has never been run in the workflow
  if (workflow.timestamp==seurat.timestamp) {
    return(TRUE)
  }
  # According to SeuratCommand, the most recent update
  # go to workflow to look up assay and DR
  params <- slot(object = object[[workflow.name]], name = "params")
  assay.use <- params$global$assay.use
  assay.use <- assay.use %iff% params[[command.name]]$assay.use
  assay.use <- assay.use %||% DefaultAssay(object)
  reduction.use <- params$global$reduction.use
  reduction.use <- reduction.use %iff% params[[command.name]]$reduction.use
  reduction.use <- reduction.use %||% formals(fun = paste0(command.name, ".Seurat"))$reduction.use

  command.name <- paste0(command.name, ".", assay.use, ".", reduction.use)
  command.name <- sub(pattern = "[\\.]+$", replacement = "", x = command.name, perl = TRUE)
  command.name <- sub(pattern = "\\.\\.", replacement = "\\.", x = command.name, perl = TRUE)

  if (length(x = command.name)==1) {
    seurat.timestamp <- slot(object = object[[command.name]], name = "time.stamp")
  }
  if (seurat.timestamp == workflow.timestamp) {
    return(FALSE)
  }
  return(TRUE)
}


#' Updates workflow timestamps
#'
#' Like the touch command in linux. Updates a workflow command's timestamp, and its dependencies
#'
#' @param object Seurat object
#' @param workflow.name Workflow name, should already be initialized using InitializeWorkflow
#' @param command.name Name of the command to touch
#' @param time.stamp Timestamp to assign
#'
#' @return Seurat object with updated workflow
#'
#' @export
#'
#' @examples
#' TouchWorkflow(object = pbmc_small,workflow.name = "cluster", command.name = "ScaleData")
#' #'
TouchWorkflow <- function(object, workflow.name, command.name, time.stamp = Sys.time()) {
  CheckWorkflow(object = object, workflow.name = workflow.name)
  #Now update all dependencies, recursively
  depends <- slot(object = object[[workflow.name]],name = "depends")
  depend.commands <- colnames(depends)[which(depends[,command.name]==1)]
  mostRecent <- slot(object[[workflow.name]],"mostRecent")
  mostRecent[command.name] <- time.stamp
  slot(object[[workflow.name]],"mostRecent") <- mostRecent
  for(i in depend.commands) {
    object <- TouchWorkflow(object,workflow.name = workflow.name,command.name = i,time.stamp = time.stamp)
  }
  return(object)
}

#' Output status of each command in the workflow
#'
#' For each command in the workflow, indicate whether it is up-to-date.
#'
#' @param object Seurat object
#' @param workflow.name Workflow name, should already be initialized using InitializeWorkflow
#' @param command.name Name of the command at the end of the workflow
#'
#' @return Seurat object with updated workflow
#'
#' @export
#'
#' @examples
#' WorkflowStatus(object = pbmc_small,workflow.name = "cluster")
#' #
WorkflowStatus <- function(object, workflow.name, command.name) {
  CheckWorkflow(object = object, workflow.name = workflow.name)
  message(paste0("Status  for ", workflow.name, " workflow"))
  depends <- slot(object = object[[workflow.name]],name = "depends")
  all.cmds <- rownames(depends)
  for(i in all.cmds) {
    is.updated <- (!CheckWorkflowUpdate(object = object,workflow.name = workflow.name,command.name = i))
    if (is.updated) {
      message(paste0("\t",i, " up to date"))
    }
    else {
      message(paste0("\t\t",i, " is out of date"))
    }
  }
}

#' Output individual function calls to recreate workflow
#'
#' Output all commands to reproduce your analysis without shortcuts. Should enhance reproducibility, but can be confused by custom modifcations, usage of SubsetData, etc.
#' We hope this will be very useful, but use with care and verify that it does indeed reproduce your work.
#'
#' @param object Seurat object
#' @param workflow.name Workflow name, should already be initialized using InitializeWorkflow
#' @param command.name Name of the command at the end of the workflow
#' @param depth depth of the recursive call. Only depth 1 outputs the parameters
#'
#' @return Seurat object with updated workflow
#'
#' @export
#'
#' @examples
#' RecreateWorkflow(object = pbmc_small,workflow.name = "cluster", command.name = "FindClusters")
#' #
RecreateWorkflows <- function(object, workflow.name, command.name,depth=1) {
  CheckWorkflow(object = object, workflow.name = workflow.name)
  depends <- slot(object = object[[workflow.name]],name = "depends")
  prereq.commands <- colnames(depends)[which(depends[command.name,]==1)]
  if (depth == 1) {
    message(paste0("\tNeed to output SetParams"))
  }
  for(i in prereq.commands) {
    RecreateWorkflows(object = object,workflow.name = workflow.name,command.name = i,depth = depth + 1)
  }
  #TODO deal with Assay better
  command.name <- intersect(c(command.name, paste0(command.name,".",DefaultAssay(object))), names(object))
  if (length(x = command.name)==1) {
    call.string <- slot(object[[command.name]],"call.string")
    #browser()
    message(paste0("\t",call.string))
  }
}
