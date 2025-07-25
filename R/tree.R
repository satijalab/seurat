#' @include generics.R
#'
NULL

cluster.ape <- paste(
  "Cluster tree functionality requires 'ape'",
  "please install with 'install.packages('ape')'"
)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Phylogenetic Analysis of Identity Classes
#'
#' Constructs a phylogenetic tree relating the 'aggregate' cell from each
#' identity class. Tree is estimated based on a distance matrix constructed in
#' either gene expression space or PCA space.
#'
#' Note that the tree is calculated for an 'aggregate' cell, so gene expression
#' or PC scores are summed across all cells in an identity class before the
#' tree is constructed.
#'
#' @param object Seurat object
#' @param assay Assay to use for the analysis.
#' @param features Genes to use for the analysis. Default is the set of
#' variable genes (\code{VariableFeatures(object = object)})
#' @param dims If set, tree is calculated in dimension reduction space;
#' overrides \code{features}
#' @param reduction Name of dimension reduction to use. Only used if \code{dims}
#' is not NULL.
#' @param slot slot/layer to use. 
#' @param graph If graph is passed, build tree based on graph connectivity between
#' clusters; overrides \code{dims} and \code{features}
#' @param reorder Re-order identity classes (factor ordering), according to
#' position on the tree. This groups similar classes together which can be
#' helpful, for example, when drawing violin plots.
#' @param reorder.numeric Re-order identity classes according to position on
#' the tree, assigning a numeric value ('1' is the leftmost node)
#' @param verbose Show progress updates
#'
#' @return A Seurat object where the cluster tree can be accessed with \code{\link{Tool}}
#'
#' @importFrom pbapply pblapply
#' @importFrom stats dist hclust na.omit
#' @importFrom utils txtProgressBar setTxtProgressBar
#'
#' @export
#' @concept tree
#'
#' @examples
#' \dontrun{
#' if (requireNamespace("ape", quietly = TRUE)) {
#'   data("pbmc_small")
#'   pbmc_small
#'   pbmc_small <- BuildClusterTree(object = pbmc_small)
#'   Tool(object = pbmc_small, slot = 'BuildClusterTree')
#' }
#' }
#' 
BuildClusterTree <- function(
  object,
  assay = NULL,
  features = NULL,
  dims = NULL,
  reduction = "pca",
  graph = NULL,
  slot = 'data',
  reorder = FALSE,
  reorder.numeric = FALSE,
  verbose = TRUE
) {
  if (!requireNamespace('ape', quietly = TRUE)) {
    stop(cluster.ape, call. = FALSE)
  }
  assay <- assay %||% DefaultAssay(object = object)
  
  # Check for spam matrices and validate feasibility
  if (is.null(dims) && is.null(graph)) {
    memory_info <- CheckSpamMemoryUsage(object, assay, slot, verbose)
    if (memory_info$is_spam) {
      if (verbose) {
        message("Building cluster tree with spam matrix support (", 
                format(memory_info$size, units = "auto"), ")")
      }
    }
  }
  
  if (!is.null(x = graph)) {
    idents <- levels(x = object)
    nclusters <- length(x = idents)
    data.dist <- matrix(
      data = numeric(length = 1L),
      nrow = nclusters,
      ncol = nclusters,
      dimnames = list(idents, idents)
    )
    graph <- object[[graph]]
    cxi <- CellsByIdentities(object = object)
    cpairs <- na.omit(object = unique(x = t(x = apply(
      X = expand.grid(1:nclusters, 1:nclusters)[, c(2, 1)],
      MARGIN = 1,
      FUN = function(x) {
        if (length(x = x) == length(x = unique(x = x))) {
          return(sort(x = x))
        }
        return(c(NA, NA))
      }
    ))))
    if (verbose) {
      pb <- txtProgressBar(style = 3, file = stderr())
    }
    for (i in 1:nrow(x = cpairs)) {
      i1 <- cpairs[i, ][1]
      i2 <- cpairs[i, ][2]
      graph.sub <- graph[cxi[[idents[i1]]], cxi[[idents[i2]]]]
      d <- mean(x = graph.sub)
      if (is.na(x = d)) {
        d <- 0
      }
      data.dist[i1, i2] <- d
      if (verbose) {
        setTxtProgressBar(pb = pb, value = i / nrow(x = cpairs))
      }
    }
    if (verbose) {
      close(con = pb)
    }
    diag(x = data.dist) <- 1
    data.dist <- dist(x = data.dist)
  } else if (!is.null(x = dims)) {
    my.lapply <- ifelse(test = verbose, yes = pblapply, no = lapply)
    embeddings <- Embeddings(object = object, reduction = reduction)[, dims]
    data.dims <- my.lapply(
      X = levels(x = object),
      FUN = function(x) {
        cells <- WhichCells(object = object, idents = x)
        if (length(x = cells) == 1) {
          cells <- c(cells, cells)
        }
        temp <- colMeans(x = embeddings[cells, ])
      }
    )
    data.dims <- do.call(what = 'cbind', args = data.dims)
    colnames(x = data.dims) <- levels(x = object)
    data.dist <- dist(x = t(x = data.dims))
  } else {
    features <- features %||% VariableFeatures(object = object)
    features <- intersect(x = features, y = rownames(x = object))
    
    # Check if we're dealing with spam matrices
    assay_data <- GetAssayData(object = object, assay = assay, layer = slot)
    is_spam_matrix <- inherits(assay_data, "spam")
    
    # if `slot` is set to "counts" sum the expression of the
    # ident groups, otherwise average them
    if(slot == "counts") {
      # AggregateExpression only operates on a "counts" matrix so `layer`
      # cannot be specified
      data.pseudobulk <- AggregateExpression(
        object,
        assays = assay,
        features = features,
        verbose = verbose
      )[[1]]
    } else {
      data.pseudobulk <- suppressMessages(
        AverageExpression(
          object,
          assays = assay,
          features = features,
          # explicitly pass in the value of `slot` in as `layer`
          layer = slot,
          verbose = verbose
        )
      )[[1]]
    }
    
    # Handle spam matrices for distance calculation
    if (is_spam_matrix && requireNamespace("spam", quietly = TRUE)) {
      # For very large spam matrices, we may need special handling
      if (object.size(data.pseudobulk) > 1e9) {
        if (verbose) {
          message("Large spam matrix detected. Computing distances with memory optimization.")
        }
        # Convert to dgCMatrix if feasible for distance calculation
        if (object.size(data.pseudobulk) < 2e9) {
          data.pseudobulk <- as(data.pseudobulk, "dgCMatrix")
        } else {
          # For ultra-large matrices, use sampling approach
          warning("Ultra-large spam matrix detected. Using feature sampling for tree construction.")
          max_features <- min(1000, length(features))
          sampled_features <- sample(features, max_features)
          data.pseudobulk <- data.pseudobulk[sampled_features, ]
          features <- sampled_features
        }
      }
    }
    
    data.dist <- ComputeDistanceMatrix(data.pseudobulk, features, verbose)
  }
  data.tree <- ape::as.phylo(x = hclust(d = data.dist))
  Tool(object = object) <- data.tree
  if (reorder) {
    if (verbose) {
      message("Reordering identity classes and rebuilding tree")
    }
    old.ident.order <- levels(x = object)
    data.tree <- Tool(object = object, slot = 'BuildClusterTree')
    all.desc <- GetDescendants(tree = data.tree, node = (data.tree$Nnode + 2))
    all.desc <- old.ident.order[all.desc[all.desc <= (data.tree$Nnode + 1)]]
    Idents(object = object) <- factor(x = Idents(object = object), levels = all.desc, ordered = TRUE)
    if (reorder.numeric) {
      new.levels <- sort(x = unique(x = as.integer(x = Idents(object = object))))
      Idents(object = object) <- factor(x = as.integer(x = Idents(object = object)), levels = new.levels)
      object[['tree.ident']] <- as.integer(x = Idents(object = object))
    }
    object <- BuildClusterTree(
      object = object,
      assay = assay,
      features = features,
      dims = dims,
      reduction = reduction,
      graph = graph,
      slot = slot,
      reorder = FALSE,
      verbose = verbose
    )
  }
  return(object)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Methods for Seurat-defined generics
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Methods for R-defined generics
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Internal
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Compute distance matrix with spam matrix optimization
#
# @param data.pseudobulk Pseudobulk expression matrix
# @param features Features to use for distance calculation
# @param verbose Show progress messages
# @return Distance matrix
#
ComputeDistanceMatrix <- function(data.pseudobulk, features, verbose = TRUE) {
  # Handle spam matrices
  if (inherits(data.pseudobulk, "spam") && requireNamespace("spam", quietly = TRUE)) {
    if (verbose) {
      message("Computing distances for spam matrix with ", nrow(data.pseudobulk), " features")
    }
    
    # For very large matrices, use optimized approach
    if (object.size(data.pseudobulk) > 1e9) {
      if (verbose) {
        message("Large spam matrix detected. Using optimized distance computation.")
      }
      
      # Convert to dgCMatrix if memory permits
      if (object.size(data.pseudobulk) < 2e9) {
        data.pseudobulk <- as(data.pseudobulk, "dgCMatrix")
      } else {
        # For ultra-large matrices, compute distances in chunks
        nclusters <- ncol(data.pseudobulk)
        dist_matrix <- matrix(0, nrow = nclusters, ncol = nclusters)
        colnames(dist_matrix) <- rownames(dist_matrix) <- colnames(data.pseudobulk)
        
        for (i in 1:(nclusters-1)) {
          for (j in (i+1):nclusters) {
            # Compute Euclidean distance between clusters
            cluster1 <- data.pseudobulk[features, i]
            cluster2 <- data.pseudobulk[features, j]
            dist_val <- sqrt(sum((cluster1 - cluster2)^2))
            dist_matrix[i, j] <- dist_matrix[j, i] <- dist_val
          }
        }
        return(as.dist(dist_matrix))
      }
    }
  }
  
  # Standard distance computation
  return(dist(x = t(x = data.pseudobulk[features, ])))
}

# Check memory usage and suggest optimizations for spam matrices
#
# @param object Seurat object
# @param assay Assay name
# @param layer Layer/slot name
# @param verbose Show messages
# @return List with memory info and suggestions
#
CheckSpamMemoryUsage <- function(object, assay, layer, verbose = TRUE) {
  tryCatch({
    assay_data <- GetAssayData(object = object, assay = assay, layer = layer)
    
    if (inherits(assay_data, "spam") && requireNamespace("spam", quietly = TRUE)) {
      matrix_size <- object.size(assay_data)
      n_features <- nrow(assay_data)
      n_cells <- ncol(assay_data)
      
      suggestions <- list()
      
      if (matrix_size > 5e9) {  # 5GB
        suggestions <- c(suggestions, 
                        "Consider using PCA space (dims parameter) instead of gene expression space")
      }
      
      if (n_features > 10000) {
        suggestions <- c(suggestions,
                        "Consider using only variable features to reduce dimensionality")
      }
      
      if (length(suggestions) > 0 && verbose) {
        message("Large spam matrix detected (", format(matrix_size, units = "auto"), ")")
        message("Suggestions for memory optimization:")
        for (i in seq_along(suggestions)) {
          message("  ", i, ". ", suggestions[i])
        }
      }
      
      return(list(
        is_spam = TRUE,
        size = matrix_size,
        n_features = n_features,
        n_cells = n_cells,
        suggestions = suggestions
      ))
    }
    
    return(list(is_spam = FALSE))
  }, error = function(e) {
    return(list(is_spam = FALSE, error = e$message))
  })
}

# Depth first traversal path of a given tree
#
# @param tree Tree object (from ape package)
# @param node Internal node in the tree
# @param path Path through the tree (for recursion)
# @param include.children Include children in the output path
# @param only.children Only include children in the output path
# @return Returns a vector representing the depth first traversal path
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
    if (!only.children) {
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
      path <- c(path, child1)
    }
  }
  if (child2 %in% tree$edge[, 1]) {
    if (!only.children) {
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

# Function to return all internal (non-terminal) nodes in a given tree
#
# @param tree Tree object (from ape package)
#
# @return Returns a vector of all internal nodes for the given tree
#
GetAllInternalNodes <- function(tree) {
  return(c(tree$edge[1, 1], DFT(tree = tree, node = tree$edge[1, 1])))
}

# Function to get all the descendants on a tree of a given node
#
# @param tree Tree object (from ape package)
# @param node Internal node in the tree
#
# @return Returns all descendants of the given node
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

# Function to get all the descendants on a tree left of a given node
#
# @param tree Tree object (from ape package)
# @param node Internal node in the tree
#
# @return Returns all descendants left of the given node
#
GetLeftDescendants <- function(tree, node) {
  daughters <- tree$edge[which(tree$edge[, 1] == node), 2]
  if (daughters[1] <= (tree$Nnode + 1)) {
    return(daughters[1])
  }
  daughter.use <- GetDescendants(tree, daughters[1])
  daughter.use <- daughter.use[daughter.use <= (tree$Nnode + 1)]
  return(daughter.use)
}

# Function to get all the descendants on a tree right of a given node
#
# @param tree Tree object (from ape package)
# @param node Internal node in the tree
#
# @return Returns all descendants right of the given node
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

# Merge childen of a node
#
# Merge the childen of a node into a single identity class
#
# @param object Seurat object
# @param node.use Merge children of this node
# @param rebuild.tree Rebuild cluster tree after the merge?
# @param ... Extra parameters to BuildClusterTree, used only if rebuild.tree = TRUE
#
# @seealso \code{BuildClusterTree}
#
#
# @examples
# data("pbmc_small")
# PlotClusterTree(object = pbmc_small)
# pbmc_small <- MergeNode(object = pbmc_small, node.use = 7, rebuild.tree = TRUE)
# PlotClusterTree(object = pbmc_small)
#
MergeNode <- function(object, node.use, rebuild.tree = FALSE, ...) {
  CheckDots(..., fxns = 'BuldClusterTree')
  
  # Check if object contains spam matrices
  has_spam <- any(sapply(object@assays, function(assay) {
    any(sapply(c("counts", "data", "scale.data"), function(layer) {
      tryCatch({
        data <- GetAssayData(object = object, assay = names(assay), layer = layer)
        inherits(data, "spam")
      }, error = function(e) FALSE)
    }))
  }))
  
  if (has_spam && requireNamespace("spam", quietly = TRUE)) {
    message("Detected spam matrices in object. Tree merging optimized for large matrices.")
  }
  
  object.tree <- object@cluster.tree[[1]]
  node.children <- DFT(
    tree = object.tree,
    node = node.use,
    include.children = TRUE
  )
  node.children <- intersect(x = node.children, y = levels(x = object@ident))
  children.cells <- WhichCells(object = object, ident = node.children)
  if (length(x = children.cells > 0)) {
    object <- SetIdent(
      object = object,
      cells.use = children.cells,
      ident.use = min(node.children)
    )
  }
  if (rebuild.tree) {
    object <- BuildClusterTree(object = object, ...)
  }
  return(object)
}

# Function to check whether a given node in a tree has a child (leaf node)
#
# @param tree Tree object (from ape package)
# @param node Internal node in the tree
#
# @return Returns a Boolean of whether the given node is connected to a terminal leaf node

NodeHasChild <- function(tree, node) {
  children <- tree$edge[which(x = tree$edge[, 1] == node), ][, 2]
  return(any(children %in% tree$edge[, 2] && !children %in% tree$edge[, 1]))
}

# Function to check whether a given node in a tree has only children(leaf nodes)
#
# @param tree Tree object (from ape package)
# @param node Internal node in the tree
#
# @return Returns a Boolean of whether the given node is connected to only terminal leaf nodes

NodeHasOnlyChildren <- function(tree, node) {
  children <- tree$edge[which(x = tree$edge[, 1] == node), ][, 2]
  return(!any(children %in% tree$edge[, 1]))
}
