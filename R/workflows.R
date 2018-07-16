#' Initialize a Seurat workflow
#'
#' Reads dependencies from a file and initializes a Seurat workflow
#'
#' @param object Seurat object
#' @param file Tab-delimited text file with two columns : First is the command, second is the dependency. Name should be workflowName.workflow.txt
#'
#' @return Object with modified workflows
#'
#' @export
#'
#' @examples
#' pbmc_small <- InitializeWorkflow(object = pbmc_small, file = 'workflows/cluster.workflow.txt)
#'
InitializeWorkflow <- function(object, file) {
  worklow.name <- gsub(".workflow.txt","",basename(file))
  workflow.data <- read.table(file)
  cmds <- sort(unique(as.vector(as.matrix(workflow.data))))
  depends <- data.frame(matrix(nrow=length(cmds),ncol=length(cmds)))
  rownames(depends) <- cmds; colnames(depends) <- cmds
  for(i in 1:nrow(workflow.data)) {
    depends[as.character(workflow.data[i,1]),as.character(workflow.data[i,2])]=1
  }
  mostRecent <- rep(NA,length(cmds)); names(updates) <- cmds 
  seurat.workflow <- new(
    Class = 'SeuratWorkflow',
    name = worklow.name,
    depends = depends,
    updates = updates
  )
  object[[worklow.name]] <- seurat.workflow
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
  workflow.present=F
  if (workflow.name%in%names(object)) {
    if (class(object[[workflow.name]])[[1]]=="SeuratWorkflow") {
      workflow.present=T
    }
  }
  if (!(workflow.present)) stop("Workflow not present, initialize first.")
  return(TRUE)
}

#' Set Seurat workflow parameters
#'
#' Sets parameters for a workflow
#'
#' @param object Seurat object
#' @param workflow.name Workflow name, should already be initialized using InitializeWorkflow
#' @param \dots Parameters to set, will be fed into workflow Seurat functions. Parameters ending with a "." will populate all similar variable names (i.e. setting dims. will set both dims.compute, and dims.cluster)
#'
#' @return Object with modified workflows
#'
#' @export
#'
#' @examples
#' pbmc_small <- SetWorkflowParams(object = pbmc_small, seed.use = 31, dims. = 20)
#'
SetWorkflowParams <- function(object, workflow.name, ...) {
  CheckWorkflow(object = object, workflow.name = workflow.name)
  # Set Params
  params <- (list(...))
  slot(object[[workflow.name]],"params") <- params
  return(object)
}


