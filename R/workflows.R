# Logs a command run, storing the name, timestamp, and argument list. Stores in the Seurat object
# object is the name of the object
# not sure if this should be internal (I think it should be)
InitializeWorkflow <- function(object, file) {
  worklow.name <- gsub(".workflow.txt","",basename(file))
  workflow.data <- read.table(file)
  cmds <- sort(unique(as.vector(as.matrix(workflow.data))))
  depends <- data.frame(matrix(nrow=length(cmds),ncol=length(cmds)))
  rownames(depends) <- cmds; colnames(depends) <- cmds
  for(i in 1:nrow(workflow.data)) {
    depends[as.character(workflow.data[i,1]),as.character(workflow.data[i,2])]=1
  }
}


