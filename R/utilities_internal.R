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

# Calculate the percentage of a vector above some threshold
#
# @param x          Vector of values
# @param threshold  Threshold to use when calculating percentage
#
# @return           Returns the percentage of `x` values above the given
#                   threshold
#
PercentAbove <- function(x, threshold){
  return(length(x = x[x > threshold]) / length(x = x))
}

# Calculate position along a defined reference range for a given vector of
# numerics. Will range from 0 to 1.
#
# @param x      Vector of numeric type
# @param lower  Lower end of reference range
# @param upper  Upper end of reference range
#
# @return       Returns a vector that describes the position of each element in
#               x along the defined reference range

ReferenceRange <- function(x, lower = 0.025, upper = 0.975) {
  return((x - quantile(x = x, probs = lower)) /
           (quantile(x = x, probs = upper) - quantile(x = x, probs = lower)))
}

# Function to map values in a vector `v` as defined in `from`` to the values
# defined in `to`.
#
# @param v     vector of values to map
# @param from  vector of original values
# @param to    vector of values to map original values to (should be of equal
#              length as from)
# @return      returns vector of mapped values
#
MapVals <- function(v, from, to){
  if (length(from) != length(to)) {
    stop("from and to vectors are not the equal length.")
  }
  vals.to.match <- match(v, from)
  vals.to.match.idx  <- !is.na(vals.to.match)
  v[vals.to.match.idx] <- to[vals.to.match[vals.to.match.idx]]
  return(v)
}

# Fills slot in new object with equivalent slot in old object if it still exists
#
# @param slot.name   slot to fill
# @param old.object  object to get slot value from
# @param new.slot    object to set slot value in
#
# @return            returns new object with slot filled
#
FillSlot <- function(slot.name, old.object, new.object){
  new.slot <- tryCatch(
    {
      slot(object = old.object, name = slot.name)
    },
    error = function(err){
      return(NULL)
    }
  )
  if(!is.null(x = new.slot)) {
    slot(new.object, slot.name) <- new.slot
  }
  return(new.object)
}


# Use Fisher's method (Fisher's combined probability test) to combine p-values
# into single statistic
#
# @param pvals vector of p-values
#
# @returns integrated value
#
FisherIntegrate <- function(pvals) {
  return(1 - pchisq(q = -2 * sum(log(x = pvals)), df = 2 * length(x = pvals)))
}

####################### Tree Related Utilities #################################

# Function to get all the descendants on a tree left of a given node
#
# @param tree  Tree object (from ape package)
# @param node  Internal node in the tree
#
# @return      Returns all descendants left of the given node
#
GetLeftDescendants <- function(tree, node) {
  daughters <- tree$edge[which(tree$edge[, 1] == node), 2]
  if (daughters[1] <= (tree$Nnode+1)) {
    return(daughters[1])
  }
  daughter.use <- GetDescendants(tree, daughters[1])
  daughter.use <- daughter.use[daughter.use <= (tree$Nnode + 1)]
  return(daughter.use)
}

# Function to get all the descendants on a tree right of a given node
#
# @param tree  Tree object (from ape package)
# @param node  Internal node in the tree
#
# @return      Returns all descendants right of the given node
#
GetRightDescendants <- function(tree, node) {
  daughters <- tree$edge[which(x = tree$edge[, 1] == node), 2]
  if (daughters[2] <= (tree$Nnode + 1)) {
    return(daughters[2])
  }
  daughter.use <- GetDescendants(tree = tree, node = daughters[2])
  daughter.use <- daughter.use[daughter.use <= (tree$Nnode + 1)]
  return(daughter.use)
}

# Function to get all the descendants on a tree of a given node
#
# @param tree  Tree object (from ape package)
# @param node  Internal node in the tree
#
# @return      Returns all descendants of the given node
#
GetDescendants <- function(tree, node, curr = NULL) {
  if (is.null(x = curr)) {
    curr <- vector()
  }
  daughters <- tree$edge[which(x = tree$edge[, 1] == node), 2]
  curr <- c(curr, daughters)
  w <- which(x = daughters >= length(x = tree$tip))
  if (length(x = w) > 0) {
    for (i in 1:length(x = w)) {
      curr <- GetDescendants(tree = tree, node = daughters[w[i]], curr = curr)
    }
  }
  return(curr)
}

# Depth first traversal path of a given tree
#
# @param tree              Tree object (from ape package)
# @param node              Internal node in the tree
# @param path              Path through the tree (for recursion)
# @param include.children  Include children in the output path
# @param only.children     Only include children in the output path
# @return                  Returns a vector representing the depth first
#                          traversal path
#
DFT <- function(
  tree,
  node,
  path = NULL,
  include.children = FALSE,
  only.children = FALSE
) {
  if (only.children) {
    include.children = TRUE
  }
  children <- which(x = tree$edge[, 1] == node)
  child1 <- tree$edge[children[1], 2]
  child2 <- tree$edge[children[2], 2]
  if (child1 %in% tree$edge[, 1]) {
    if(! only.children){
      path <- c(path, child1)
    }
    path <- DFT(
      tree = tree,
      node = child1,
      path = path,
      include.children = include.children,
      only.children = only.children
    )
  } else {
    if (include.children) {
      path <-c(path, child1)
    }
  }
  if (child2 %in% tree$edge[, 1]) {
    if (! only.children) {
      path <- c(path, child2)
    }
    path <- DFT(
      tree = tree,
      node = child2,
      path = path,
      include.children = include.children,
      only.children = only.children
    )
  } else {
    if (include.children) {
      path <- c(path, child2)
    }
  }
  return(path)
}

# Function to check whether a given node in a tree has a child (leaf node)
#
# @param tree   Tree object (from ape package)
# @param node   Internal node in the tree
#
# @return       Returns a Boolean of whether the given node is connected to a
#               terminal leaf node

NodeHasChild <- function(tree, node) {
  children <- tree$edge[which(x = tree$edge[, 1] == node), ][, 2]
  return(any(children %in% tree$edge[, 2] && ! children %in% tree$edge[, 1]))
}

# Function to check whether a given node in a tree has only children(leaf nodes)
#
# @param tree   Tree object (from ape package)
# @param node   Internal node in the tree
#
# @return       Returns a Boolean of whether the given node is connected to only
#               terminal leaf nodes

NodeHasOnlyChildren <- function(tree, node) {
  children <- tree$edge[which(x = tree$edge[, 1] == node), ][, 2]
  return(! any(children %in% tree$edge[, 1]))
}

# Function to return all internal (non-terminal) nodes in a given tree
#
# @param tree   Tree object (from ape package)
#
# @return       Returns a vector of all internal nodes for the given tree
#
GetAllInternalNodes <- function(tree) {
  return(c(tree$edge[1, 1], DFT(tree = tree, node = tree$edge[1, 1])))
}

# Weighted Euclidean Distance
#
# @param x Dataset 1
# @param y Dataset 2
# @param w Weights
#
# @return The Weighted Euclidian Distance (numeric)
#
WeightedEuclideanDistance <- function(x, y, w) {
  v.dist <- sum(sqrt(x = w * (x - y) ^ 2))
  return(v.dist)
}

# Set a default value if an object is null
#
# @param x An object to set if it's null
# @param default The value to provide if x is null
#
# @return default if x is null, else x
#
SetIfNull <- function(x, default) {
  return(ifelse(
    test = is.null(x = x),
    yes = default,
    no = x
  ))
}
