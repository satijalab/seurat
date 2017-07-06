# Internal function for merging two matrices by rowname
#
# @param mat1 First matrix
# @param mat2 Second matrix
#
# @return A merged matrix
#
RowMergeSparseMatrices <- function(mat1, mat2){
  if (inherits(x = mat1, what = "data.frame")) {
    mat1 <- as.matrix(x = mat1)
  }
  if (inherits(x = mat2, what = "data.frame")) {
    mat2 <- as.matrix(x = mat2)
  }
  mat1 <- as(object = mat1, Class = "RsparseMatrix")
  mat2 <- as(object = mat2, Class = "RsparseMatrix")
  mat1.names <- rownames(x = mat1)
  mat2.names <- rownames(x = mat2)
  all.names <- union(x = mat1.names, y = mat2.names)
  new.mat <- RowMergeMatrices(
    mat1 = mat1,
    mat2 = mat2,
    mat1_rownames = mat1.names,
    mat2_rownames = mat2.names,
    all_rownames = all.names
  )
  rownames(x = new.mat) <- make.unique(names = all.names)
  colnames(x = new.mat) <- make.unique(names = c(
    colnames(x = mat1),
    colnames(x = mat2)
  ))
  return(new.mat)
}
