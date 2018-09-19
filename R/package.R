.onLoad <- function(libname, pkgname) {
  delayedAssign(
    x = 'cluster.workflow',
    value = CreateWorkflowObject(file = system.file(
      'extdata/cluster.workflow.ini',
      package = 'Seurat'
    )),
    assign.env = asNamespace(ns = 'Seurat')
  )
  # autoload(name = 'cluster.workflow', package = 'Seurat')
}
