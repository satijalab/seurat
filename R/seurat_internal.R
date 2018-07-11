FilterObjects <- function(object, classes.keep = c('Assay', 'DimReduc')) {
  object.classes <- sapply(
    X = names(x = object),
    FUN = function(i) {
      return(class(x = object[[i]]))
    }
  )
  object.classes <- object.classes[object.classes %in% classes.keep]
  return(names(x = object.classes))
}
