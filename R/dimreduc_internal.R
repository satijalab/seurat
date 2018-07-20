Top <- function(
  data.use,
  dim.use,
  num.use,
  do.balanced
) {
  top <- if (do.balanced) {
    num.use <- round(x = num.use / 2)
    data.use <- data.use[order(data.use), , drop = FALSE]
    positive <- head(x = rownames(x = data.use), n = num.use)
    negative <- rev(x = tail(x = rownames(x = data.use), n = num.use))
    list(positive = positive, negative = negative)
  } else {
    data.use <- data.use[rev(x = order(abs(x = data.use))), , drop = FALSE]
    top <- head(x = rownames(x = data.use), n = num.use)
    top[order(data.use[top, ])]
  }
  return(top)
}
