# Tests to verify the RCpp version of ModularityOptimizer produces the same
# results as the java version.
# Equivalent java commands are given above.
context("ModularityOptimizer")

# The "karate club" network available from the ModularityOptimizer website at:
# http://www.ludowaltman.nl/slm/
node1 <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1,
           1, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 4, 4, 5, 5, 5, 6, 8, 8, 8, 9, 13,
           14, 14, 15, 15, 18, 18, 19, 20, 20, 22, 22, 23, 23, 23, 23, 23, 24,
           24, 24, 25, 26, 26, 27, 28, 28, 29, 29, 30, 30, 31, 31, 32)
node2 <- c(1, 2, 3, 4, 5, 6, 7, 8, 10, 11, 12, 13, 17, 19, 21, 31, 2, 3, 7, 13,
           17, 19, 21, 30, 3, 7, 8, 9, 13, 27, 28, 32, 7, 12, 13, 6, 10, 6, 10,
           16, 16, 30, 32, 33, 33, 33, 32, 33, 32, 33, 32, 33, 33, 32, 33, 32,
           33, 25, 27, 29, 32, 33, 25, 27, 31, 31, 29, 33, 33, 31, 33, 32, 33,
           32, 33, 32, 33, 33)
dim_s <- max(max(node1), max(node2)) + 1
# Note we want to represent network in the lower diagonal.
connections <- sparseMatrix(i = node2 + 1, j = node1 + 1, x = 1.0)

# Result from equivalent command to
# java -jar ModularityOptimizer.jar karate_club_network.txt communities.txt 1 1.0 1 1 1 564 0
test_that("Algorithm 1", {
  expected <- c(1, 1, 1, 1, 2, 2, 2, 1, 0, 1, 2, 1, 1, 1, 0, 0, 2, 1, 0, 1, 0, 1,
               0, 0, 3, 3, 0, 0, 3, 0, 0, 3, 0, 0)
  s <- Seurat:::RunModularityClusteringCpp(
    SNN = connections,
    modularityFunction = 1,
    resolution = 1.0,
    algorithm = 1,
    nRandomStarts = 1,
    nIterations = 1,
    randomSeed = 564,
    printOutput = 0,
    ""
  )
  expect_equal(expected, s)
})


#java -jar ModularityOptimizer.jar karate_club_network.txt communities.txt 1 1.0 2 1 1 2 0
test_that("Algorithm 2", {
  expected <- c(1, 1, 1, 1, 3, 3, 3, 1, 0, 0, 3, 1, 1, 1, 0, 0, 3, 1, 0, 1, 0, 1,
               0, 2, 2, 2, 0, 2, 2, 0, 0, 2, 0, 0)
  s <- Seurat:::RunModularityClusteringCpp(
    SNN = connections,
    modularityFunction = 1,
    resolution = 1.0,
    algorithm = 2,
    nRandomStarts = 1,
    nIterations = 1,
    randomSeed = 2,
    printOutput = 0,
    ""
  )
  expect_equal(expected, s)
})

#java -jar ModularityOptimizer.jar karate_club_network.txt communities.txt 1 1.0 3 1 1 56464 0
test_that("Algorithm 3", {
  expected <- c(1, 1, 1, 1, 3, 3, 3, 1, 0, 0, 3, 1, 1, 1, 0, 0, 3, 1, 0, 1, 0, 1,
               0, 2, 2, 2, 0, 2, 2, 0, 0, 2, 0, 0)
  s <- Seurat:::RunModularityClusteringCpp(
    SNN = connections,
    modularityFunction = 1,
    resolution = 1.0,
    algorithm = 3,
    nRandomStarts = 1,
    nIterations = 1,
    randomSeed = 56464,
    printOutput = 0,
    "")
  expect_equal(expected, s)
})

test_that("Low Resolution", {
  e1 <- rep(0, 34)
  # java -jar ModularityOptimizer.jar karate_club_network.txt outjava.txt  1 0.05 3 1 10 10 0
  s <- Seurat:::RunModularityClusteringCpp(
    SNN = connections,
    modularityFunction = 1,
    resolution = 0.05,
    algorithm = 3,
    nRandomStarts = 1,
    nIterations = 10,
    randomSeed = 10,
    printOutput = 0,
    ""
  )
  expect_equal(s, e1)
  # java -jar ModularityOptimizer.jar karate_club_network.txt outjava.txt 2 0.05 3 1 10 10 0
  s2 <- Seurat:::RunModularityClusteringCpp(
    SNN = connections,
    modularityFunction = 2,
    resolution=0.05,
    algorithm = 3,
    nRandomStarts = 1,
    nIterations = 10,
    randomSeed = 10,
    printOutput = 0,
    ""
  )
  e2 = c(0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 0, 1, 0, 1, 1,
         1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
  expect_equal(s2, e2)
})

test_that("EdgeWeights", {
  # Make 1, 4, 5 and 20 a community by weighting them
  c2 <- connections
  c2[5, 4] <- 3.0
  c2[5, 1] <- 5.0
  c2[4, 1] <- 8.0
  c2[20, 5] <- 8.0
  c2[20, 4] <- 5.0
  c2[20, 1] <- 5.0
  # java -jar ModularityOptimizer.jar weighted_karate_club_network.txt outjava.txt  1 1.0 3 1 10 40 1
  s2 <- Seurat:::RunModularityClusteringCpp(
    SNN = c2,
    modularityFunction = 1,
    resolution = 1.0,
    algorithm = 3,
    nRandomStarts = 1,
    nIterations = 10,
    randomSeed = 40,
    printOutput = 0,
    ""
  )
  exp <- c(2, 1, 1, 2, 2, 3, 3, 1, 0, 1, 3, 2, 2, 1, 0, 0, 3, 1, 0, 2, 0, 1, 0,
          0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
  expect_equal(s2, exp)
})

# test_that("pbmc_small network", {
#   observed <- as.numeric(FindClusters(
#     object = pbmc_small,
#     reduction.type = "pca",
#     dims.use = 1:10,
#     resolution = 1.1,
#     save.SNN = TRUE,
#     print.output = 0)@ident)
#   expected = c(1,1,1,1,1,1,1,1,1,1,6,1,6,1,2,2,1,6,2,1,2,2,2,2,2,2,2,2,2,6,3,5,3,3,3,3,3,3,3,3,5,1,1,1,1,1,3,1,3,1,2,1,2,2,6,2,3,2,1,3,5,2,5,5,2,2,2,2,5,3,4,4,4,4,4,4,4,4,4,4)
#   expect_equal(observed, expected)
# })
