#%% spam coercion / interoperability %%#
#
# The spam package (with spam64) provides in-memory sparse matrices with 64-bit
# indexing, which can represent more than 2^31 non-zero entries. However, spam
# objects do NOT carry dimnames, while Seurat assays require named features and
# cells. spam therefore cannot serve as a Seurat layer backend directly; for
# out-of-memory ultra-large layers use the DelayedMatrix or BPCells backends
# instead. We support spam here only for coercion/interoperability, so a spam
# matrix can be converted into a Seurat-compatible sparse matrix when it fits in
# the 2^31 limit.
NULL

#' @method as.sparse spam
#' @export
#'
as.sparse.spam <- function(x, ...) {
  if (!requireNamespace('spam', quietly = TRUE)) {
    stop("Package 'spam' must be installed to coerce a spam matrix.",
         call. = FALSE)
  }
  # spam stores values row-wise (CSR); build a dgCMatrix from its slots without
  # densifying when possible, falling back to a dense round-trip otherwise.
  converted <- tryCatch(
    expr = as(object = x, Class = 'dgCMatrix'),
    error = function(e) as(object = as.matrix(x = x), Class = 'dgCMatrix')
  )
  return(converted)
}
