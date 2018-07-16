#' @export
#'
'.DollarNames.SeuratCommand' <- function(x, pattern = ''){
  params <- slot(object = x, name = "params")
  utils:::findMatches(pattern, names(params))
}

#' @export
#'
'$.SeuratCommand' <- function(x, i, ...) {
  params <- slot(object = x, name = "params")
  return(params[[i]])
}


#' @export
#'
'[.SeuratCommand' <- function(x, i, ...) {
  slot.use <- c("name", "timestamp", "call_string", "params")
  if (!i %in% slot.use) {
    stop("Invalid slot")
  }
  return(slot(object = x, name = i))
}

setMethod(
  f = 'show',
  signature = 'SeuratCommand',
  definition = function(object) {
    params <- slot(object = object, name = "params")
    cat(
      "Command: ", slot(object = object, name = "call.string"), '\n',
      "Time: ", as.character(slot(object = object, name = "time.stamp")), '\n',
      sep = ""
    )
    for(p in 1:length(params)){
      cat(
       names(params[p]), ":", params[[p]], "\n"
      )
    }
  }
)
