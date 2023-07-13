#' @include zzz.R
#'
NULL

#' @importFrom utils lsf.str
#' @importFrom rlang is_scalar_character
#'
#' @export
#'
.rd_methods <- function(method = 'integration') {
  methods <- sapply(
    X = grep(pattern = '^package:', x = search(), value = TRUE),
    FUN = function(x) {
      fxns <- as.character(x = lsf.str(pos = x))
      attrs <- vector(mode = 'logical', length = length(x = fxns))
      for (i in seq_along(along.with = fxns)) {
        mthd <- attr(x = get(x = fxns[i], pos = x), which = 'Seurat.method')
        attrs[i] <- is_scalar_character(x = mthd) && isTRUE(x = mthd == method)
      }
      return(fxns[attrs])
    },
    simplify = FALSE,
    USE.NAMES = TRUE
  )
  methods <- Filter(f = length, x = methods)
  names(x = methods) <- gsub(
    pattern = '^package:',
    replacement = '',
    x = names(x = methods)
  )
  if (!length(x = methods)) {
    return('')
  }
  ret <- vector(
    mode = 'character',
    length = sum(vapply(
      X = methods,
      FUN = length,
      FUN.VALUE = integer(length = 1L)
    ))
  )
  j <- 1L
  for (pkg in names(x = methods)) {
    for (fxn in methods[[pkg]]) {
      ret[j] <- ifelse(
        test = pkg == 'Seurat',
        yes = paste0('\\item \\code{\\link{', fxn, '}}'),
        no = paste0(
          '\\item \\code{\\link[',
          pkg,
          ':',
          fxn,
          ']{',
          pkg,
          '::',
          fxn, '}}'
        )
      )
      j <- j + 1L
    }
  }
  return(paste('\\itemize{', paste0(' ', ret, collapse = '\n'), '}', sep = '\n'))
}
