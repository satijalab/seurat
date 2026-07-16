.GetSeuratNThreads <- function() {
  nthreads <- getOption(x = "Seurat.nthreads", default = 1L)
  nthreads <- suppressWarnings(expr = as.integer(x = nthreads[[1L]]))
  if (length(x = nthreads) != 1L || is.na(x = nthreads) || nthreads < 1L) {
    nthreads <- 1L
  }
  return(nthreads)
}

.FindVariableFeaturesVSTInfo <- function(
  object,
  loess.span = 0.3,
  clip.max = "auto",
  verbose = TRUE
) {
  if (clip.max == "auto" || is.null(x = clip.max)) {
    clip.max <- sqrt(x = ncol(x = object))
  }
  hvf.info <- as.data.frame(
    x = SparseRowMeanVar(
      x = object@x,
      i = object@i,
      p = object@p,
      rows = nrow(x = object),
      cols = ncol(x = object),
      nthreads = .GetSeuratNThreads(),
      display_progress = verbose
    )
  )
  hvf.info$variance.expected <- 0
  not.const <- hvf.info$variance > 0
  fit <- loess(
    formula = log10(x = variance) ~ log10(x = mean),
    data = hvf.info[not.const, ],
    span = loess.span
  )
  hvf.info$variance.expected[not.const] <- 10 ^ fit$fitted
  hvf.info$variance.standardized <- SparseRowVarStd(
    x = object@x,
    i = object@i,
    p = object@p,
    mu = hvf.info$mean,
    sd = sqrt(x = hvf.info$variance.expected),
    rows = nrow(x = object),
    cols = ncol(x = object),
    vmax = clip.max,
    nthreads = .GetSeuratNThreads(),
    display_progress = verbose
  )
  return(hvf.info)
}
