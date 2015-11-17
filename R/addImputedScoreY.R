## =======
## Improve imputation results
## - CV-similar approach to improve lasso
##
## =======



#' Calculate imputed expression values
#'
#' Uses L1-constrained linear models (LASSO) to impute single cell gene
#' expression values. Optional: Use cross-validation similar approach to 
#' improve LASSO model. 
#'
#' @param object Seurat object
#' @param genes.use A vector of genes (predictors) that can be used for
#' building the LASSO models.
#' @param genes.fit A vector of genes to impute values for
#' @param s.use Maximum number of steps taken by the algorithm (lower values
#' indicate a greater degree of smoothing)
#' @param cross Logic. Determine whether apply CV-similar approach to improve LASSO (default: FALSE). 
#' @param cross.step Numeric. Determine how many rounds of iteration of CV-similar approach is performed (default: 0). 
#' @param cross.threshold Numeric. If difference of imputed gene expression level is 
#' not greater than \code{cross.threshold} after each iteration, looping ternimates. (default: Inf). 
#' @param do.print Print progress (output the name of each gene after it has
#' been imputed).
#' @param gram The use.gram argument passed to lars
#' @return Returns a Seurat object where the imputed values have been added to
#' object@@data
#' @import lars
#' @export
setGeneric("addImputedScoreY", function(object, genes.use = NULL, 
  genes.fit = NULL, s.use = 20, cross = FALSE, cross.step = 0, cross.threshold = Inf, do.print = FALSE, gram = TRUE) standardGeneric("addImputedScoreY"))
#' @export
setMethod("addImputedScoreY", "seurat", function(object, genes.use = NULL,
                                                 genes.fit = NULL, s.use = 20, cross = FALSE, cross.step = 0, cross.threshold = Inf, do.print = FALSE, gram = TRUE) {
  genes.use <- set.ifnull(genes.use, object@var.genes)
  genes.fit <- set.ifnull(genes.fit, object@var.genes)
  genes.use <- genes.use[genes.use %in% rownames(object@data)] ## lasso-genes
  genes.fit <- genes.fit[genes.fit %in% rownames(object@data)] ## insitu-genes

  lasso.input <- t(object@data[genes.use, ])
  lasso.fits <- data.frame(t(sapply(genes.fit, 
              function(x) lasso.fxn(t(object@data[genes.use[genes.use != x], ]), 
                                    object@data[x, ], 
                                    s.use = s.use, 
                                    x, 
                                    do.print, gram))))
  genes.old <- genes.fit[genes.fit %in% rownames(object@imputed)]
  genes.new <- genes.fit[!(genes.fit %in% rownames(object@imputed))]
  
  if (length(genes.old) > 0) {
    object@imputed[genes.old, ] <- lasso.fits[genes.old, ]
  }
  object@imputed <- rbind(object@imputed, lasso.fits[genes.new, ])
  
  if (cross == TRUE) {
    ## Perform CV-similar method
    cat("#== Perform CV-similar approach\n")
    step <- 0
    imputed_data <- object@imputed
    imputed_new <- imputed_data
    while (step < cross.step) {
      imputed_new <- data.frame(t(sapply(genes.fit,
                                        function(x) {
                                          genes_use_data <- object@data[genes.use[genes.use != x],]
                                          genes_help <- genes.fit[genes.fit != x]
                                          
                                          lasso.fxn(t(rbind(genes_use_data, imputed_data[genes_help, ])),
                                                    imputed_data[x, ],
                                                    s.use = s.use,
                                                    x,
                                                    do.print, gram)
                                        })))
      imputed_data <- imputed_new
      step <- step + 1
    }
    object@imputed <- imputed_new
  } ## ends improvement
  return(object)
}) 

## Inner part of addImputedScoreY
## Could use wrapper to run recursively, without running from scratch
setGeneric("addImputedScoreSup", function(object, genes.use = NULL,
                                        genes.fit = NULL, s.use = 20, cross = FALSE, cross.step = 0, cross.threshold = Inf, do.print = FALSE, gram = TRUE)
  standardGeneric("addImputedScoreSup"))
setMethod("addImputedScoreSup", "seurat", function(object, genes.use = NULL,
                                                   genes.fit = NULL, s.use = 20, cross = FALSE, cross.step = 0, cross.threshold = Inf, do.print = FALSE, gram = TRUE) {
  if (cross == TRUE) {
    ## Perform CV-similar method
    cat("#== Perform CV-similar approach\n")
    step <- 0
    imputed_data <- object@imputed
    imputed_new <- imputed_data
    while (step < cross.step) {
      imputed_new <- data.frame(t(sapply(genes.fit,
                                         function(x) {
                                           genes_use_data <- object@data[genes.use[genes.use != x],]
                                           genes_help <-
                                             genes.fit[genes.fit != x]
                                           
                                           lasso.fxn(t(rbind(genes_use_data, imputed_data[genes_help,])),
                                                     imputed_data[x,],
                                                     s.use = s.use,
                                                     x,
                                                     do.print, gram)
                                         })))
      imputed_data <- imputed_new
      step <- step + 1
    }
    object@imputed <- imputed_new
  } ## ends improvement
  return(object)
})