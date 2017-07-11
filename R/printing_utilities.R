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
               GetCalcParam(object = object, calculation = "PCA", parameter = "time"), "\n"))
    cat("=============================================================================\n")
    cat(paste0("PCs computed    Genes used in calculation    PCs Scaled by Variance Explained\n"))
    cat(paste0("    ", GetCalcParam(object = object, calculation = "PCA", parameter = "pcs.compute"), 
               "\t\t\t\t", length(GetCalcParam(object = object, calculation = "PCA", 
                                               parameter = "pc.genes")), "\t\t\t",
               GetCalcParam(object = object, calculation = "PCA", parameter = "scale.by.varexp"), "\n"))
    cat("-----------------------------------------------------------------------------\n")
    cat("rev.pca \n")
    cat(paste0(" ", GetCalcParam(object = object, calculation = "PCA", parameter = "rev.pca"), "\n"))
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
               GetCalcParam(object = object, calculation = "ICA", parameter = "time"), "\n"))
    cat("=============================================================================\n")
    cat(paste0("ICs computed \t Genes used in calculation \t ICA function \t rev.ica \n"))
    cat(paste0("    ", GetCalcParam(object = object, calculation = "ICA", parameter = "ics.compute"), 
               "\t\t\t\t", length(GetCalcParam(object = object, calculation = "ICA", 
                                               parameter = "ic.genes")), "\t\t  ",
               GetCalcParam(object = object, calculation = "ICA", parameter = "ica.function"), "\t  ",
               GetCalcParam(object = object, calculation = "ICA", parameter = "rev.ica"), "\n"))
    cat("-----------------------------------------------------------------------------\n")
    cat("Full gene list can be accessed at object@@calc.params$ICA$pc.genes")
  }
}

#' Print TSNE Calculation Parameters
#'
#' Print the parameters chosen for the latest stored TSNE calculation.
#' 
#' @param object Seurat object
#' @param raw Print the entire contents of the calculation metadata slot (calc.params) for the 
#' RunTSNE calculation. Default (FALSE) will print a nicely formatted summary.
#' @return No return value. Only prints to console.
#' @export
PrintTSNEParams <- function(object, raw = FALSE){
  if (raw){
    print(object@calc.params$RunTSNE)
  }
  else{
    cat(paste0("Parameters used in latest TSNE calculation run on: ", 
               GetCalcParam(object = object, calculation = "RunTSNE", parameter = "time"), "\n"))
    cat("=============================================================================\n")
    if(is.null(GetCalcParam(object = object, calculation = "RunTSNE", parameter = "genes.use"))) {
      reduction <- GetCalcParam(object = object, calculation = "RunTSNE", parameter = "reduction.use")
      dim <- "Dims"
      n.dim <- GetCalcParam(object = object, calculation = "RunTSNE", parameter = "dims.use")
    } else if (!is.null(GetCalcParam(object = object, calculation = "RunTSNE", 
                                     parameter = "distance.matrix"))){
      reduction <- "custom"
      dim <- "Custom distance matrix"
    } else {
      reduction <- "None"
      dim <- "Genes"
      n.dim <- length(GetCalcParam(object = object, calculation = "RunTSNE", parameter = "genes.use"))
    }
    cat(paste0("Reduction use          do.fast          dim.embed\n"))
    cat(paste0("     ", reduction, "                 ", 
               GetCalcParam(object, "RunTSNE", "do.fast"),  "              ", 
               GetCalcParam(object, "RunTSNE", "dim.embed"), "              ", "\n")) 
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
