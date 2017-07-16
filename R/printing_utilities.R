#' Print the calculation
#'
#' Print entire contents of calculation settings slot (calc.params) for given
#' calculation.
#'
#' @param object Seurat object
#' @param calculation Name of calculation (function name) to check parameters
#' for
#' @param raw Print the entire contents of the calculation settings slot (calc.params)
#' for the RunPCA calculation.
#' @param return.list Return the calculation parameters as a list
#' @return Prints the calculation settings and optionally returns them as a list
#' @export
PrintCalcParams <- function(object, calculation, raw = FALSE,
                            return.list = FALSE) {
  if(is.null(object@calc.params[[calculation]])){
    stop(paste0(calculation, " not computed yet."))
  }
  if (!raw){
    if(calculation == "RunPCA"){
      PrintPCAParams(object)
    }
    if(calculation == "ICA"){
      PrintICAParams(object)
    }
    if(calculation == "RunTSNE"){
      PrintTSNEParams(object)
    }
    if(calculation == "RunCCA"){
      PrintCCAParams(object)
    }
    if(calculation == "CalcVarExpRatioParams"){
      PrintCalcVarExpRatioParams(object)
    }
    if(calculation == "AlignSubspace"){
      PrintAlignSubspaceParams(object)
    }
    if(calculation == "RunDiffusion"){
      PrintDMParams(object)
    }
    if(calculation == "BuildSNN"){
      PrintSNNParams(object)
    }
    if(calculation == "FindClusters"){
      PrintFindClustersParams(object)
    }
  }
  print(object@calc.params[[calculation]])
  if(return.list){
    return(object@calc.params[[calculation]])
  }
}

#' Print PCA Calculation Parameters
#'
#' Print the parameters chosen for the latest stored PCA calculation.
#'
#' @param object Seurat object
#' @param raw Print the entire contents of the calculation settings slot
#' (calc.params) for the RunPCA calculation. Default (FALSE) will print a nicely
#' formatted summary.
#' @return No return value. Only prints to console.
#' @export
PrintPCAParams <- function(object, raw = FALSE){
  if(is.null(object@calc.params$RunPCA)){
    stop("PCA has not been computed yet")
  }
  if (raw){
    print(object@calc.params$RunPCA)
  }
  else{
    cat(paste0("Parameters used in latest PCA calculation run on: ",
               GetCalcParam(object = object,
                            calculation = "RunPCA",
                            parameter = "time"),
               "\n"))
    cat("=============================================================================\n")
    cat(paste0("PCs computed    Genes used in calculation    PCs Scaled by Variance Explained\n"))
    pcs.compute <- GetCalcParam(object = object,
                                calculation = "RunPCA",
                                parameter = "pcs.compute")
    n.gene <- length(GetCalcParam(object = object,
                                  calculation = "RunPCA",
                                  parameter = "pc.genes"))
    cat(paste0("    ",
               pcs.compute,
               FillWhiteSpace(n = 20 - nchar(pcs.compute)),
               n.gene, FillWhiteSpace(n = 35 - nchar(n.gene)),
               GetCalcParam(object = object,
                            calculation = "RunPCA",
                            parameter = "weight.by.var"),
               "\n"))
    cat("-----------------------------------------------------------------------------\n")
    cat("rev.pca \n")
    cat(paste0(" ",
               GetCalcParam(object = object,
                            calculation = "RunPCA",
                            parameter = "rev.pca"),
               "\n"))
    cat("-----------------------------------------------------------------------------\n")
    cat("Full gene list can be accessed using \n GetCalcParam(object = object, calculation = \"RunPCA\", parameter = \"pc.genes\")")
  }
}

#' Print ICA Calculation Parameters
#'
#' Print the parameters chosen for the latest stored ICA calculation.
#'
#' @param object Seurat object
#' @param raw Print the entire contents of the calculation settings slot (calc.params) for the ICA
#' calculation. Default (FALSE) will print a nicely formatted summary.
#' @return No return value. Only prints to console.
#' @export
PrintICAParams <- function(object, raw = FALSE){
  if(is.null(object@calc.params$ICA)){
    stop("ICA has not been computed yet")
  }
  if (raw){
    print(object@calc.params$ICA)
  }
  else{
    cat(paste0("Parameters used in latest ICA calculation run on: ",
               GetCalcParam(object = object,
                            calculation = "ICA",
                            parameter = "time"),
               "\n"))
    cat("=============================================================================\n")
    cat(paste0("ICs computed \t Genes used in calculation \t ICA function \t rev.ica \n"))
    ics.compute <- GetCalcParam(object = object,
                                calculation = "ICA",
                                parameter = "ics.compute")
    n.genes <- length(GetCalcParam(object = object,
                                   calculation = "ICA",
                                   parameter = "ic.genes"))
    ica.fxn <- GetCalcParam(object = object,
                            calculation = "ICA",
                            parameter = "ica.function")
    rev.ica <- GetCalcParam(object = object,
                            calculation = "ICA",
                            parameter = "rev.ica")
    cat(paste0("    ",
               ics.compute,
               FillWhiteSpace(n = 25 - nchar(ics.compute)),
               n.genes,
               FillWhiteSpace(n = 22 - nchar(n.genes)),
               ica.fxn,
               FillWhiteSpace(n = 15 - nchar(ica.fxn)),
               rev.ica,"\n"))
    cat("-----------------------------------------------------------------------------\n")
    cat("Full gene list can be accessed using \n GetCalcParam(object = object, calculation = \"ICA\", parameter = \"ic.genes\")")
  }
}

#' Print TSNE Calculation Parameters
#'
#' Print the parameters chosen for the latest stored TSNE calculation.
#'
#' @param object Seurat object
#' @param raw Print the entire contents of the calculation settings slot (calc.params) for the
#' RunTSNE calculation. Default (FALSE) will print a nicely formatted summary.
#' @return No return value. Only prints to console.
#' @export
PrintTSNEParams <- function(object, raw = FALSE){
  if(is.null(object@calc.params$RunTSNE)){
    stop("TSNE has not been computed yet")
  }
  if (raw){
    print(object@calc.params$RunTSNE)
  }
  else{
    cat(paste0("Parameters used in latest TSNE calculation run on: ",
               GetCalcParam(object = object,
                            calculation = "RunTSNE",
                            parameter = "time"),
               "\n"))
    cat("=============================================================================\n")

    if(is.null(GetCalcParam(object = object,
                            calculation = "RunTSNE",
                            parameter = "genes.use"))) {
      reduction <- GetCalcParam(object = object,
                                calculation = "RunTSNE",
                                parameter = "reduction.use")
      dim <- "Dims"
      n.dim <- GetCalcParam(object = object,
                            calculation = "RunTSNE",
                            parameter = "dims.use")
    } else if (!is.null(GetCalcParam(object = object,
                                     calculation = "RunTSNE",
                                     parameter = "distance.matrix"))){
      reduction <- "custom"
      dim <- "Custom distance matrix"
    } else {
      reduction <- "None"
      dim <- "Genes"
      n.dim <- length(GetCalcParam(object = object,
                                   calculation = "RunTSNE",
                                   parameter = "genes.use"))
    }
    do.fast <- GetCalcParam(object = object,
                            calculation = "RunTSNE",
                            parameter = "do.fast")
    dim.embed <- GetCalcParam(object = object,
                              calculation = "RunTSNE",
                              parameter = "dim.embed")
    cat(paste0("Reduction use          do.fast          dim.embed\n"))
    cat(paste0("     ",
               reduction,
               FillWhiteSpace(n = 19 - nchar(reduction)),
               do.fast,
               FillWhiteSpace(n = 20 - nchar(do.fast)),
               dim.embed,
               "\n"))
    cat("-----------------------------------------------------------------------------\n")
    cat(paste0(dim, " used in calculation\n"))
    cat("=============================================================================\n")
    if(reduction == "None"){
      cat(paste0(n.dim, " genes used: Full gene list can be accessed using \n GetCalcParam(object = object, calculation = \"RunTSNE\", parameter = \"genes.use\")"))
    } else if (reduction == "custom") {
      cat("Full matrix can be acccessed using \n GetCalcParam(object = object, calculation = \"RunTSNE\", parameter = \"distance.matrix\")")
    } else {
      cat(paste0(strwrap(paste(n.dim, "\n", collapse = " "), width = 80),
                 collapse = "\n"))
      cat("\n\n")
    }
  }
}

#' Print CCA Calculation Parameters
#'
#' Print the parameters chosen for the latest stored CCA calculation.
#'
#' @param object Seurat object
#' @param raw Print the entire contents of the calculation settings slot
#' (calc.params) for the RunCCA calculation. Default (FALSE) will print a nicely
#' formatted summary.
#' @return No return value. Only prints to console.
#' @export
PrintCCAParams <- function(object, raw = FALSE){
  if(is.null(object@calc.params$RunCCA)){
    stop("CCA has not been computed yet")
  }
  if (raw){
    print(object@calc.params$RunCCA)
  }
  else{
    cat(paste0("Parameters used in latest CCA calculation run on: ",
               GetCalcParam(object = object,
                            calculation = "RunCCA",
                            parameter = "time"), "\n"))
    cat("=============================================================================\n")

    cat(paste0("CCs computed        Genes used in calculation        scale.data\n"))
    num.cc <- GetCalcParam(object = object,
                           calculation = "RunCCA",
                           parameter = "num.cc")
    num.genes <- length(GetCalcParam(object = object,
                                     calculation = "RunCCA",
                                     parameter = "genes.use"))
    cat(paste0("    ",
               num.cc ,
               FillWhiteSpace(28 - nchar(num.cc)),
               num.genes,
               FillWhiteSpace(n = 24 - nchar(num.genes)),
               GetCalcParam(object = object,
                            calculation = "RunCCA",
                            parameter = "scale.data"),
               "\n"))
    cat("-----------------------------------------------------------------------------\n")
    g1 <- GetCalcParam(object = object,
                       calculation = "RunCCA",
                       parameter = "group1")
    g2 <- GetCalcParam(object = object,
                       calculation = "RunCCA",
                       parameter = "group2")
    if(nchar(g1) > 0){
      cat(paste0("group1",
                 FillWhiteSpace(n = 10),
                 "group2",
                 FillWhiteSpace(n = 10),
                 "group.by",
                 FillWhiteSpace(n = 10),
                 "rescale.groups\n"))
      gb <- GetCalcParam(object = object,
                         calculation = "RunCCA",
                         parameter = "group.by")
      if(length(g1) > 1) {
        g1 <- "custom group"
      }
      if(length(g2) > 1){
        g2 <- "custom group"
      }
      rsg <- GetCalcParam(object = object,
                          calculation = "RunCCA",
                          parameter = "rescale.groups")
      cat(paste0(g1,
                 FillWhiteSpace(n = 15 - nchar(g1)),
                 g2,
                 FillWhiteSpace(n = 18 - nchar(g2)),
                 gb,
                 FillWhiteSpace(n = 18 - nchar(rsg)),
                 rsg ,
                 "\n"))
      cat("-----------------------------------------------------------------------------\n")
    }
    if(!is.null(GetCalcParam(object = object,
                             calculation = "RunCCA",
                             parameter = "object.project"))){
      n1 <- GetCalcParam(object = object,
                         calculation = "RunCCA",
                         parameter = "object.project")
      n2 <- GetCalcParam(object = object,
                         calculation = "RunCCA",
                         parameter = "object2.project")
      cat("Object 1 Project Name        Object 2 Project Name\n")
      cat(paste0("  ",
                 n1,
                 FillWhiteSpace(n = 30 - nchar(n1)),
                 n2,
                 "\n"))
      cat("-----------------------------------------------------------------------------\n")
    }
    if(g2 == "custom group" | g2 == "custom group"){
      cat("Group membership lists can be accessed using \n GetCalcParam(object = object, calculation = \"RunCCA\", parameter = \"group1/2\")\n")
    }
    cat("Full gene list can be accessed using \n GetCalcParam(object = object, calculation = \"RunCCA\", parameter = \"genes.use\")")
  }
}

#' Print Parameters Associated with CalcVarExpRatio
#'
#' Print the parameters chosen for CalcVarExpRatio.
#'
#' @param object Seurat object
#' @param raw Print the entire contents of the calculation settings slot
#' (calc.params) for CalcVarExpRatio. Default (FALSE) will print a nicely
#' formatted summary.
#' @return No return value. Only prints to console.
#' @export
PrintCalcVarExpRatioParams <- function(object, raw = FALSE){
  if(is.null(object@calc.params$CalcVarExpRatio)){
    stop("CalcVarExpRatio has not been computed yet")
  }
  if (raw){
    print(object@calc.params$CalcVarExpRatio)
  }
  else{
    cat(paste0("Parameters used in latest CalcVarExpRatio run on: ",
               GetCalcParam(object = object,
                            calculation = "PCA",
                            parameter = "time"),
               "\n"))
    cat("=============================================================================\n")
    cat(paste0("Reduction Type    Grouping Variable \n"))
    reduction <- GetCalcParam(object = object,
                              calculation = "CalcVarExpRatio",
                              parameter = "reduction.type")
    grouping.var <- GetCalcParam(object = object,
                                 calculation = "CalcVarExpRatio",
                                 parameter = "grouping.var")
    dims.use <- GetCalcParam(object = object,
                             calculation = "CalcVarExpRatio",
                             parameter = "dims.use")
    cat(paste0(" ",
               reduction,
               FillWhiteSpace(n = 20 - nchar(reduction)),
               grouping.var,
               "\n"))
    cat("-----------------------------------------------------------------------------\n")
    cat("Dims used in calculation\n")
    cat(paste0(strwrap(paste(dims.use, "\n", collapse = " "), width = 80),
               collapse = "\n"))
  }
}

#' Print AlignSubspace Calculation Parameters
#'
#' Print the parameters chosen for the latest AlignSubspace calculation for each
#' stored aligned subspace.
#'
#' @param object Seurat object
#' @param raw Print the entire contents of the calculation settings slot
#' (calc.params) for the AlignSubspace calculation. Default (FALSE) will print a
#' nicely formatted summary.
#' @return No return value. Only prints to console.
#' @export
PrintAlignSubspaceParams <- function(object, raw = FALSE){
  to.print <- names(object@calc.params)[grepl("AlignSubspace.",
                                              names(object@calc.params))]
  if(length(to.print) == 0){
    stop("No stored aligned subspaces.")
  }
  for (i in to.print){
    if (raw){
      print(object@calc.params[[i]])
    }
    else{
      cat(paste0("Parameters used in latest AlignSubspace calculation run on: ",
                 GetCalcParam(object = object,
                              calculation = i,
                              parameter = "time"),
                 "\n"))
      cat("=============================================================================\n")
      reduction <- GetCalcParam(object = object,
                                calculation = i,
                                parameter = "reduction.type")
      grouping.var <- GetCalcParam(object = object,
                                   calculation = i,
                                   parameter = "grouping.var")
      dims <- GetCalcParam(object = object,
                           calculation = i,
                           parameter = "dims.align")
      n.genes <- GetCalcParam(object = object,
                              calculation = i,
                              parameter = "num.genes")
      cat(paste0("Reduction use          grouping.var          num.genes\n"))
      cat(paste0("     ",
                 reduction,
                 FillWhiteSpace(n = 19 - nchar(reduction)),
                 grouping.var,
                 FillWhiteSpace(n = 15 - nchar(dims)),
                 n.genes,
                 "\n"))
      cat("-----------------------------------------------------------------------------\n")
      cat("Dims aligned\n")
      cat("=============================================================================\n")
      cat(paste0(strwrap(paste(dims, "\n", collapse = " "), width = 80),
                 collapse = "\n"))
      cat("\n")
    }
  }
}

#' Print Diffusion Map Calculation Parameters
#'
#' Print the parameters chosen for the latest stored diffusion map calculation.
#'
#' @param object Seurat object
#' @param raw Print the entire contents of the calculation settings slot
#' (calc.params) for the RunDiffusion calculation. Default (FALSE) will print a
#' nicely formatted summary.
#' @return No return value. Only prints to console.
#' @export
PrintDMParams <- function(object, raw = FALSE){
  if(is.null(object@calc.params$RunDiffusion)){
    stop("Diffusion map has not been computed yet")
  }
  if (raw){
    print(object@calc.params$RunDiffusion)
  }
  else{
    cat(paste0("Parameters used in latest diffusion map calculation run on: ",
               GetCalcParam(object = object,
                            calculation = "RunDiffusion",
                            parameter = "time"), "\n"))
    cat("=============================================================================\n")
    max.dim <- GetCalcParam(object = object,
                            calculation = "RunDiffusion",
                            parameter = "max.dim")
    reduction <- GetCalcParam(object = object,
                              calculation = "RunDiffusion",
                              parameter = "reduction.use")
    n.genes <- length(GetCalcParam(object = object,
                            calculation = "RunDiffusion",
                            parameter = "genes.use"))
    scale.clip <- length(GetCalcParam(object = object,
                                      calculation = "RunDiffusion",
                                      parameter = "scale.clip"))
    q.use <- GetCalcParam(object = object,
                          calculation = "RunDiffusion",
                          parameter = "q.use")
    dims.use <- GetCalcParam(object = object,
                             calculation = "RunDiffusion",
                             parameter = "dims.use")
    if(n.genes > 0){
      reduction <- "None"
    }
    cat(paste0("Reduction used    DMs computed    Quantile    scale.clip \n"))
    cat(paste0("    ",
               reduction ,
               FillWhiteSpace(20 - nchar(reduction)),
               max.dim,
               FillWhiteSpace(n = 12 - nchar(max.dim)),
               q.use,
               FillWhiteSpace(n = 15 - nchar(q.use)),
               scale.clip,
               "\n"))
    cat("-----------------------------------------------------------------------------\n")
    if(reduction == "None"){
      dim <- "Genes"
    }
    else{
      dim <- "Dims"
    }
    cat(paste0(dim, " used in calculation\n"))
    cat("=============================================================================\n")
    if(reduction == "None"){
      cat(paste0(n.genes, " genes used: Full gene list can be accessed using \n GetCalcParam(object = object, calculation = \"RunDiffusion\", parameter = \"genes.use\")"))
    } else {
      cat(paste0(strwrap(paste(dims.use, "\n", collapse = " "), width = 80),
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
#' @param raw Print the entire contents of the calculation settings slot (calc.params) for the
#' BuildSNN calculation. Default (FALSE) will print a nicely formatted summary.
#' @return No return value. Only prints to console.
#' @export
PrintSNNParams <- function(object, raw = FALSE){
  if(is.null(object@calc.params$BuildSNN)){
    stop("SNN has not been computed yet")
  }
  if (raw){
    print(object@calc.params$BuildSNN)
  }
  else{
    cat(paste0("Parameters used in latest SNN calculation run on: ",
               GetCalcParam(object = object,
                            calculation = "BuildSNN",
                            parameter = "time"),
               "\n"))
    cat("=============================================================================\n")
    if(is.null(GetCalcParam(object = object,
                            calculation = "BuildSNN",
                            parameter = "genes.use")))
      {
      reduction <- GetCalcParam(object = object,
                                calculation = "BuildSNN",
                                parameter = "reduction.type")
      dim <- "Dims"
      n.dim <- GetCalcParam(object = object,
                            calculation = "BuildSNN",
                            parameter = "dims.use")
    } else if (!is.null(GetCalcParam(object = object,
                                     calculation = "BuildSNN",
                                     parameter = "distance.matrix")))
      {
        reduction <- "custom"
        dim <- "Custom distance matrix"
    } else {
        reduction <- "None"
        dim <- "Genes"
        n.dim <- length(GetCalcParam(object = object,
                                     calculation = "BuildSNN",
                                     parameter = "genes.use"))
    }
    cat(paste0("Reduction used          k.param          k.scale          prune.SNN\n"))
    k.param <- GetCalcParam(object = object,
                            calculation = "BuildSNN",
                            parameter = "k.param")
    k.scale <- GetCalcParam(object = object,
                            calculation = "BuildSNN",
                            parameter = "k.scale")
    prune.SNN <- GetCalcParam(object = object,
                              calculation = "BuildSNN",
                              parameter = "prune.SNN")
    cat(paste0("     ",
               reduction,
               FillWhiteSpace(n = 20 - nchar(reduction)),
               k.param,
               FillWhiteSpace(n = 18 - nchar(k.param)),
               k.scale,
               FillWhiteSpace(n = 16 - nchar(k.scale)),
               round(prune.SNN, 4),
               "\n"))
    cat("-----------------------------------------------------------------------------\n")
    cat(paste0(dim, " used in calculation\n"))
    cat("=============================================================================\n")
    if(reduction == "None"){
      cat(paste0(n.dim, " genes used: Full gene list can be accessed using \n GetCalcParam(object = object, calculation = \"BuildSNN\", parameter = \"genes.use\")"))
    } else if (reduction == "custom") {
        cat("Full matrix can be acccessed using \n GetCalcParam(object = object, calculation = \"RunTSNE\", parameter = \"distance.matrix\")")
    } else {
      cat(paste0(strwrap(paste(n.dim, "\n", collapse = " "), width = 80),
                 collapse = "\n"))
      cat("\n\n")
    }
  }
}

#' Print FindClusters Calculation Parameters
#'
#' Print the parameters chosen for the latest FindClusters calculation for each
#' stored resolution.
#'
#' @param object Seurat object
#' @param resolution Optionally specify only a subset of resolutions to print
#' parameters for.
#' @param raw Print the entire contents of the calculation settings slot
#' (calc.params) for the FindClusters calculation. Default (FALSE) will print a
#' nicely formatted summary.
#' @return No return value. Only prints to console.
#' @export
PrintFindClustersParams <- function(object, resolution, raw = FALSE){
  to.print <- names(object@calc.params)[grepl("FindClusters",
                                              names(object@calc.params))]
  if(length(to.print) == 0){
    stop("No stored clusterings.")
  }
  for (i in to.print){
    if(!missing(resolution)){
      if(!ExtractField(i, 2, "res.") %in% resolution){
        next
      }
    }
    if (raw){
      print(object@calc.params[[i]])
    }
    else{
      cat(paste0("Parameters used in latest FindClusters calculation run on: ",
                 GetCalcParam(object = object,
                              calculation = i,
                              parameter = "time"),
                 "\n"))
      resolution <- GetCalcParam(object = object,
                                calculation = i,
                                parameter = "resolution")
      cat("=============================================================================\n")
      cat(paste0("Resolution: ", resolution, "\n"))
      cat("-----------------------------------------------------------------------------\n")
      cat("Modularity Function    Algorithm         n.start         n.iter\n")
      modularity.fxn <- GetCalcParam(object = object,
                                     calculation = i,
                                     parameter = "modularity.fxn")
      algorithm <- GetCalcParam(object = object,
                                calculation = i,
                                parameter = "algorithm")
      n.start <- GetCalcParam(object = object,
                              calculation = i,
                              parameter = "n.start")
      n.iter <- GetCalcParam(object = object,
                             calculation = i,
                             parameter = "n.iter")
      cat(paste0("     ",
                 modularity.fxn,
                 FillWhiteSpace(n = 20 - nchar(modularity.fxn)),
                 algorithm,
                 FillWhiteSpace(n = 18 - nchar(algorithm)),
                 n.start,
                 FillWhiteSpace(n = 16 - nchar(n.start)),
                 n.iter,
                 "\n"))
      cat("-----------------------------------------------------------------------------\n")
      if (is.null(GetCalcParam(object = object,
                              calculation = i,
                              parameter = "genes.use")))
      {
        reduction <- GetCalcParam(object = object,
                                  calculation = i,
                                  parameter = "reduction.type")
        dim <- "Dims"
        n.dim <- GetCalcParam(object = object,
                              calculation = i,
                              parameter = "dims.use")
      } else if (!is.null(GetCalcParam(object = object,
                                       calculation = "BuildSNN",
                                       parameter = "distance.matrix")))
      {
        reduction <- "custom"
        dim <- "Custom distance matrix"
      } else {
        reduction <- "None"
        dim <- "Genes"
        n.dim <- length(GetCalcParam(object = object,
                                     calculation = "BuildSNN",
                                     parameter = "genes.use"))
      }
      cat(paste0("Reduction used          k.param          k.scale          prune.SNN\n"))
      k.param <- GetCalcParam(object = object,
                              calculation = "BuildSNN",
                              parameter = "k.param")
      k.scale <- GetCalcParam(object = object,
                              calculation = "BuildSNN",
                              parameter = "k.scale")
      prune.SNN <- GetCalcParam(object = object,
                                calculation = "BuildSNN",
                                parameter = "prune.SNN")
      cat(paste0("     ",
                 reduction,
                 FillWhiteSpace(n = 20 - nchar(reduction)),
                 k.param,
                 FillWhiteSpace(n = 18 - nchar(k.param)),
                 k.scale,
                 FillWhiteSpace(n = 16 - nchar(k.scale)),
                 round(prune.SNN, 4),
                 "\n"))
      cat("-----------------------------------------------------------------------------\n")
      cat(paste0(dim, " used in calculation\n"))
      cat("=============================================================================\n")
      if(reduction == "None"){
        cat(paste0(n.dim, " genes used: Full gene list can be accessed using \n GetCalcParam(object = object, calculation = \"BuildSNN\", parameter = \"genes.use\")"))
      } else if (reduction == "custom") {
        cat("Full matrix can be acccessed using \n GetCalcParam(object = object, calculation = \"RunTSNE\", parameter = \"distance.matrix\")")
      } else {
        cat(paste0(strwrap(paste(n.dim, "\n", collapse = " "), width = 80),
                   collapse = "\n"))
        cat("\n\n")
      }
    }
  }
}

