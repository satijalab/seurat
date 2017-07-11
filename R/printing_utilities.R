#' Print PCA Calculation Parameters
#'
#' Print the parameters chosen for the latest stored PCA calculation.
#' 
#' @param object Seurat object
#' @param raw Print the entire contents of the calculation metadata slot (calc.params) for the PCA
#' calculation. Default (FALSE) will print a nicely formatted summary.
#' @return No return value. Only prints to console.
#' @export
PrintPCAParams <- function(object, raw = FALSE){
  if (raw){
    print(object@calc.params$PCA)
  }
  else{
    cat(paste0("Parameters used in latest PCA calculation run on: ", 
               GetCalcParam(object, "PCA", "time"), "\n"))
    cat("=============================================================================\n")
    cat(paste0("PCs computed    Genes used in calculation    PCs Scaled by Variance Explained\n"))
    cat(paste0("    ", GetCalcParam(object, "PCA", "pcs.compute"), "\t\t\t\t", 
               length(GetCalcParam(object, "PCA", "pc.genes")), "\t\t\t",
               GetCalcParam(object, "PCA", "scale.by.varexp"), "\n"))
    cat("-----------------------------------------------------------------------------\n")
    cat("rev.pca \n")
    cat(paste0(" ", GetCalcParam(object, "PCA", "rev.pca"), "\n"))
    cat("-----------------------------------------------------------------------------\n")
    cat("Full gene list can be accessed at object@@calc.params$PCA$pc.genes")
  }
}

#' Print ICA Calculation Parameters
#'
#' Print the parameters chosen for the latest stored ICA calculation.
#' 
#' @param object Seurat object
#' @param raw Print the entire contents of the calculation metadata slot (calc.params) for the ICA
#' calculation. Default (FALSE) will print a nicely formatted summary.
#' @return No return value. Only prints to console.
#' @export
PrintICAParams <- function(object, raw = FALSE){
  if (raw){
    print(object@calc.params$ICA)
  }
  else{
    cat(paste0("Parameters used in latest ICA calculation run on: ", 
               GetCalcParam(object, "ICA", "time"), "\n"))
    cat("=============================================================================\n")
    cat(paste0("ICs computed \t Genes used in calculation \t rev.ica \n"))
    cat(paste0("    ", GetCalcParam(object, "ICA", "ics.compute"), "\t\t\t\t", 
               length(GetCalcParam(object, "ICA", "ic.genes")), "\t\t ",
               GetCalcParam(object, "ICA", "rev.ica"), "\n"))
    cat("-----------------------------------------------------------------------------\n")
    cat("Full gene list can be accessed at object@@calc.params$ICA$pc.genes")
  }
}

#' Print SNN Construction Calculation Parameters
#'
#' Print the parameters chosen for the latest stored SNN calculation (via BuildSNN or FindClusters).
#' 
#' @param object Seurat object
#' @param raw Print the entire contents of the calculation metadata slot (calc.params) for the 
#' BuildSNN calculation. Default (FALSE) will print a nicely formatted summary.
#' @return No return value. Only prints to console.
#' @export
PrintSNNParams <- function(object, raw = FALSE){
  if (raw){
    print(object@calc.params$BuildSNN)
  }
  else{
    cat(paste0("Parameters used in latest SNN calculation run on: ", 
               GetCalcParam(object, "BuildSNN", "time"), "\n"))
    cat("=============================================================================\n")
    if(is.null(object@calc.params$BuildSNN$genes.use)) {
      reduction <- GetCalcParam(object, "BuildSNN", "reduction.type")
      dim <- "Dims"
      n.dim <- GetCalcParam(object, "BuildSNN", "dims.use")
    } else if (!is.null(GetCalcParam(object, "BuildSNN", "distance.matrix"))){
      reduction <- "custom"
      dim <- "Custom distance matrix"
    } else {
      reduction <- "None"
      dim <- "Genes"
      n.dim <- length(GetCalcParam(object, "BuildSNN", "genes.use"))
    }
    cat(paste0("Reduction use          k.param          k.scale          prune.SNN\n"))
    cat(paste0("     ", reduction, "                 ", 
               GetCalcParam(object, "BuildSNN", "k.param"),  "               ", 
               GetCalcParam(object, "BuildSNN", "k.scale"), "              ",
               round(GetCalcParam(object, "BuildSNN", "prune.SNN"), 4), "\n")) 
    cat("-----------------------------------------------------------------------------\n")
    
    cat(paste0(dim, " used in calculation\n"))
    cat("=============================================================================\n")
    if(reduction == "None"){
      cat(paste0(n.dim, ": Full gene list can be accessed at object@@calc.params$BuildSNN$genes.use"))
    } else if (reduction == "custom") {
        cat("Full matrix can be acccessed at object@@calc.params$BuildSNN$distance.matrix")
    } else {
      cat(paste0(strwrap(paste(n.dim, "\n", collapse = " "), width = 80), 
                 collapse = "\n"))
      cat("\n\n")
    }
  }
}
