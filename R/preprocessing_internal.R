#' Regress out technical effects and cell cycle
#'
#' Remove unwanted effects from scale.data
#'
#' @keywords internal
#' @param object Seurat object
#' @param latent.vars effects to regress out
#' @param genes.regress gene to run regression for (default is all genes)
#' @param model.use Use a linear model or generalized linear model (poisson, negative binomial) for the regression. Options are 'linear' (default), 'poisson', and 'negbinom'
#' @param use.umi Regress on UMI count data. Default is FALSE for linear modeling, but automatically set to TRUE if model.use is 'negbinom' or 'poisson'
#'
#' @return Returns the residuals from the regression model
#'
#' @import Matrix
#'
RegressOut <- function(
  object,
  latent.vars,
  genes.regress = NULL,
  model.use = 'linear',
  use.umi = FALSE
) {
  possible.models <- c("linear", "poisson", "negbinom")
  if (! model.use %in% possible.models){
    stop(
      paste0(
        model.use,
        " is not a valid model. Please use one the following: ",
        paste0(possible.models, collapse = ", "),
        "."
      )
    )
  }
  genes.regress <- SetIfNull(x = genes.regress, default = rownames(x = object@data))
  genes.regress <- ainb(genes.regress,rownames(object@data))
  latent.data <- FetchData(object = object, vars.all = latent.vars)
  bin.size <- 100
  if (model.use == 'negbinom') {
    bin.size <- 5
  }
  bin.ind <- ceiling(x = 1:length(x = genes.regress) / bin.size)
  max.bin <- max(bin.ind)
  print(paste("Regressing out", latent.vars))
  pb <- txtProgressBar(min = 0, max = max.bin, style = 3)
  data.resid <- c()
  data.use <- object@data[genes.regress, , drop = FALSE];
  if (model.use != "linear") {
    use.umi <- TRUE
  }
  if (use.umi) {
    data.use <- object@raw.data[genes.regress, object@cell.names, drop = FALSE]
  }
  for (i in 1:max.bin) {
    genes.bin.regress <- rownames(x = data.use)[bin.ind == i]
    gene.expr <- as.matrix(x = data.use[genes.bin.regress, , drop = FALSE])
    new.data <- do.call(
      rbind,
      lapply(
        X = genes.bin.regress,
        FUN = function(x) {
          regression.mat <- cbind(latent.data, gene.expr[x,])
          colnames(x = regression.mat) <- c(colnames(x = latent.data), "GENE")
          fmla <- as.formula(
            object = paste0(
              "GENE ",
              " ~ ",
              paste(latent.vars, collapse = "+")
            )
          )
          if (model.use == 'linear') {
            return(lm(formula = fmla, data = regression.mat)$residuals)
          }
          if (model.use == 'poisson') {
            return(residuals(
              object = glm(
                formula = fmla,
                data = regression.mat,
                family = "poisson"
              ),
              type='pearson'
            ))
          }
          if (model.use == 'negbinom') {
            return(nb.residuals(
              formula = fmla,
              regression.mat = regression.mat,
              gene = x
            ))
          }
        }
      )
    )
    if (i == 1) {
      data.resid=new.data
    }
    if (i > 1) {
      data.resid=rbind(data.resid,new.data)
    }
    setTxtProgressBar(pb, i)
  }
  close(pb)
  rownames(x = data.resid) <- genes.regress
  if (use.umi) {
    data.resid <- log1p(
      x = sweep(
        x = data.resid,
        MARGIN = 1,
        STATS = apply(X = data.resid, MARGIN = 1, FUN = min),
        FUN = "-"
      )
    )
  }
  return(data.resid)
}