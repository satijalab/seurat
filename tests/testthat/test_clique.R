# Unit testing for functions in clique.cpp

# --------------------------------------------------------------------------------

context("SNN-clique Rcpp testing")
# Tests for removeRedundantClique(IntegerVector x, IntegerVector y)
# --------------------------------------------------------------------------------
test_that("removeRedundantClique works correctly",{
  x1 <- c(2, 4)
  x2 <- c(1, 2, 3, 4)
  
  expect_that(removeRedundantClique(x1,x2), is_a("logical"))
  expect_that(removeRedundantClique(x1,x2), is_true())
  expect_that(removeRedundantClique(c(5),x2), is_false())
  expect_that(removeRedundantClique(x1,x2), is_true())
})

# Tests for sizeCliqueIntersection(IntegerVector x, IntegerVector y)
# --------------------------------------------------------------------------------
test_that("sizeCliqueIntersection works correctly",{
  x1 <- c(1, 2, 3)
  x2 <- c(1, 4, 5)
  
  expect_that(sizeCliqueIntersection(x1,x1), is_a("integer"))
  expect_that(sizeCliqueIntersection(x1,x1), equals(length(x1)))
  expect_that(sizeCliqueIntersection(x1,x2), equals(1))
  expect_that(sizeCliqueIntersection(x1,c(0)), equals(0))
})

# Tests for IntegerVector removeNode(IntegerVector x, int y)
# --------------------------------------------------------------------------------
test_that("removeNode works correctly",{
  # remember that C++ indexes from 0
  x <- c(1, 2, 3)
  y <- 2
  
  expect_that(removeNode(x,y), is_a("integer"))
  expect_that(length(removeNode(x,y)), equals(length(x) - 1))
  expect_that(removeNode(x,y), equals(c(1, 2)))
})

# Tests for whichNotZero(NumericVector X)
# --------------------------------------------------------------------------------
test_that("whichNotZero works correctly",{
  v1 <- c(0, 1, 0, 1)
  non_zero <- c(1, 3) # c++ indexes from 0
  v2 <- c(0, 0, 0)
  
  expect_that(whichNotZero(v1), is_a("integer"))
  expect_that(whichNotZero(v1), equals(non_zero))
  expect_that(whichNotZero(v1), equals(c(1, 3)))
  expect_that(length(whichNotZero(v2)), equals(0))
})

# Tests for subsetMatrix(NumericMatrix m, NumericVector rows, Numeric vector cols)
# --------------------------------------------------------------------------------
test_that("subsetMatrix works correctly",{
  m <- matrix(1:16, 4, 4)
  
  expect_that(subsetMatrix(m, c(1, 2), c(1, 2)), is_a("matrix"))
  expect_that(subsetMatrix(m, c(0:3), c(0:3)), equals(m))
  expect_that(subsetMatrix(m, c(0, 1), c(0, 1)), equals(matrix(c(1, 2, 5, 6), 2, 2)))
  expect_that(subsetMatrix(m, c(0), c(0, 1)), equals(matrix(c(1, 5), 1, 2)))
  expect_that(subsetMatrix(m, c(0, 1), c(0)), equals(matrix(c(1, 2), 2, 1)))
  expect_that(subsetMatrix(m, c(0, 1), c(2, 3)), equals(matrix(c(9, 10, 13, 14), 2, 2)))
})

# Tests for NumericMatrix setRow(NumericMatrix m, int r, int n)
# --------------------------------------------------------------------------------
test_that("setRow works correctly",{
  m <- matrix(1:16, 4, 4)
  m2 <- matrix(1:16, 4, 4)
  m2[1, ] <- 0
  m3 <- matrix(1:16, 4, 4)
  m3[2, ] <- 1
  
  expect_that(setRow(m, 0, 0), is_a("matrix"))
  expect_that(setRow(m, 0, 0), equals(m2))
  expect_that(setRow(m, 1, 1), equals(m3))
})

# Tests for NumericMatrix setCol(NumericMatrix m, int c, int n)
# --------------------------------------------------------------------------------
test_that("setCol works correctly",{
  m <- matrix(1:16, 4, 4)
  m2 <- matrix(1:16, 4, 4)
  m2[ ,1] <- 0
  m3 <- matrix(1:16, 4, 4)
  m3[ ,2] <- 1
  
  expect_that(setCol(m, 0, 0), is_a("matrix"))
  expect_that(setCol(m, 0, 0), equals(m2))
  expect_that(setCol(m, 1, 1), equals(m3))
})
