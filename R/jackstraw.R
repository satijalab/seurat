
#' Determine statistical significance of PCA scores.
#'
#' Randomly permutes a subset of data, and calculates projected PCA scores for
#' these 'random' genes. Then compares the PCA scores for the 'random' genes
#' with the observed PCA scores to determine statistical signifance. End result
#' is a p-value for each gene's association with each principal component.
#'
#' @param object Seurat object
#' @param num.pc Number of PCs to compute significance for
#' @param num.replicate Number of replicate samplings to perform
#' @param prop.freq Proportion of the data to randomly permute for each
#' replicate
#' @param do.print Print the number of replicates that have been processed.
#' @param rev.pca By default computes the PCA on the cell x gene matrix. Setting to
#' true will compute it on gene x cell matrix. This should match what was set when the intial PCA was run.
#' @param do.fast Compute the PCA using the fast approximate calculation from the IRLBA package. Values stored with object
#' must also have been computed using the PCAFast() function.
#'
#' @return Returns a Seurat object where object@@jackStraw.empP represents
#' p-values for each gene in the PCA analysis. If ProjectPCA is subsequently
#' run, object@@jackStraw.empP.full then represents p-values for all genes.
#'
#' @importFrom pbapply pbsapply
#'
#' @references Inspired by Chung et al, Bioinformatics (2014)
#'
#' @export
#'
JackStraw <- function(
  object,
  num.pc = 30,
  num.replicate = 100,
  prop.freq = 0.01,
  do.print = FALSE,
  rev.pca = FALSE,
  do.fast = FALSE
) {
  # check that PCA calculation method matches
  if (do.fast) {
    if (is.null(x = object@pca.obj[[1]]$d)) {
      stop("For fast JackStraw, store PCA values computed with PCAFast()")
    }
  } else {
    if (is.null(x = object@pca.obj[[1]]$sdev)) {
      stop("For regular JackStraw, store PCA values computed with PCA()")
    }
  }
  # error checking for number of PCs
  if (num.pc > ncol(x = object@pca.rot)) {
    num.pc <- ncol(x = object@pca.rot)
    warning("Number of PCs specified is greater than PCs available. Setting num.pc to ", num.pc, " and continuing.")
  }
  if (num.pc > length(x = object@cell.names)) {
    num.pc <- length(x = object@cell.names)
    warning("Number of PCs specified is greater than number of cells. Setting num.pc to ", num.pc, " and continuing.")
  }
  pc.genes <- rownames(x = object@pca.x)
  if (length(x = pc.genes) < 3) {
    stop("Too few variable genes")
  }
  if (length(x = pc.genes) * prop.freq < 3) {
    warning(
      "Number of variable genes given ",
      prop.freq,
      " as the prop.freq is low. Consider including more variable genes and/or increasing prop.freq. ",
      "Continuing with 3 genes in every random sampling."
    )
  }
  md.x <- as.matrix(x = object@pca.x)
  md.rot <- as.matrix(x = object@pca.rot)
  if (do.print) {
    applyFunction <- pbsapply
  } else {
    applyFunction <- sapply
  }
  fake.pcVals.raw <- applyFunction(
    X = 1:num.replicate,
    FUN = function(x)
      return(JackRandom(
        scaled.data = object@scale.data[pc.genes, ],
        prop = prop.freq,
        r1.use = 1,
        r2.use = num.pc,
        seed.use = x,
        rev.pca = rev.pca,
        do.fast = do.fast
      )),
    simplify = FALSE
  )
  fake.pcVals <- sapply(
    X = 1:num.pc,
    FUN = function(x) {
      return(as.numeric(x = unlist(x = lapply(
        X = 1:num.replicate,
        FUN = function(y) {
          return(fake.pcVals.raw[[y]][, x])
        }
      ))))
    }
  )
  object@jackStraw.fakePC = data.frame(fake.pcVals)
  object@jackStraw.empP <- data.frame(
    sapply(
      X = 1:num.pc,
      FUN = function(x) {
        return(unlist(x = lapply(
          X = abs(md.x[, x]),
          FUN = empP,
          nullval = abs(fake.pcVals[,x])
        )))
      }
    )
  )
  colnames(x = object@jackStraw.empP) <- paste0("PC", 1:ncol(x = object@jackStraw.empP))
  return(object)
}

# Documentation
#' @export
#'
JackStrawPermutationTest <- function(
  object,
  genes.use = NULL,
  num.iter = 100,
  thresh.use = 0.05,
  do.print = TRUE,
  k.seed = 1
) {
  genes.use <- set.ifnull(x = genes.use, y = rownames(x = object@pca.x))
  genes.use <- ainb(a = genes.use, b = rownames(x = object@scale.data))
  data.use <- t(x = as.matrix(x = object@scale.data[genes.use, ]))
  if (do.print) {
    print(paste("Running", num.iter, "iterations"))
  }
  pa.object <- permutationPA(
    data.use,
    B = num.iter,
    threshold = thresh.use,
    verbose = do.print,
    seed = k.seed
  )
  if (do.print) {
    cat("\n\n")
  }
  if (do.print) {
    print(paste("JackStraw returns", pa.object$r, "significant components"))
  }
  return(pa.object)
}

# Documentation
#multicore version of jackstraw
#DOES NOT WORK WITH WINDOWS
#' @export
#'
JackStrawMC <- function(
  object,
  num.pc = 30,
  num.replicate = 100,
  prop.freq = 0.01,
  do.print = FALSE,
  num.cores = 8
) {
  pc.genes <- rownames(x = object@pca.x)
  if (length(x = pc.genes) < 200) {
    prop.freq <- max(prop.freq, 0.015)
  }
  md.x <- as.matrix(x = object@pca.x)
  md.rot <- as.matrix(x = object@pca.rot)
  if (do.print) {
    fake.pcVals.raw <- mclapply(
      X = 1:num.replicate,
      FUN = function(x) {
        print(x)
        return(JackRandom(
          scaled.data = object@scale.data[pc.genes, ],
          prop = prop.freq,
          r1.use = 1,
          r2.use = num.pc,
          seed.use = x
        ))
      },
      mc.cores = num.cores
    )
  } else {
    fake.pcVals.raw <- mclapply(
      X = 1:num.replicate,
      FUN = function(x) {
        return(JackRandom(
          scaled.data = object@scale.data[pc.genes, ],
          prop = prop.freq,
          r1.use = 1,
          r2.use = num.pc,
          seed.use=x
        ))
      }, mc.cores = num.cores
    )
  }
  fake.pcVals <- simplify2array(
    x = mclapply(
      X = 1:num.pc,
      FUN = function(x) {
        return(as.numeric(x = unlist(x = lapply(
          X = 1:num.replicate,
          FUN = function(y) {
            return(fake.pcVals.raw[[y]][, x])
          }
        ))))
      },
      mc.cores = num.cores
    )
  )
  object@jackStraw.fakePC <- data.frame(fake.pcVals)
  object@jackStraw.empP <- data.frame(
    simplify2array(
      x = mclapply(
        X = 1:num.pc,
        FUN = function(x) {
          return(unlist(x = lapply(
            X = abs(md.x[, x]),
            FUN = empP,
            nullval = abs(x = fake.pcVals[, x])
          )))
        },
        mc.cores = num.cores
      )
    )
  )
  colnames(x = object@jackStraw.empP) <- paste0("PC", 1:ncol(x = object@jackStraw.empP))
  return(object)
}

# Documentation
###############
JackStrawFull <- function(
  object,
  num.pc = 5,
  num.replicate = 100,
  prop.freq = 0.01
) {
  pc.genes <- rownames(x = object@pca.x)
  if (length(x = pc.genes) < 200) {
    prop.freq <- max(prop.freq, 0.015)
  }
  md.x <- as.matrix(x = object@pca.x)
  md.rot <- as.matrix(x = object@pca.rot)
  real.fval <- sapply(
    X = 1:num.pc,
    FUN = function(x) {
      return(unlist(x = lapply(
        X = pc.genes,
        FUN = jackF,
        r1 = x,
        r2 = x,
        x = md.x,
        rot = md.rot
      )))
    }
  )
  rownames(x = real.fval) <- pc.genes
  object@real.fval <- data.frame(real.fval)
  fake.fval <- sapply(
    X = 1:num.pc,
    FUN = function(x) {
      return(unlist(x = replicate(
        n = num.replicate,
        expr = jackStrawF(
          prop = prop.freq,
          data = object@scale.data[pc.genes, ],
          myR1 = x,
          myR2 = x
        ),
        simplify = FALSE
      )))
    }
  )
  rownames(x = fake.fval) <- 1:nrow(x = fake.fval)
  object@fake.fval <- data.frame(fake.fval)
  object@emp.pval <- data.frame(
    sapply(
      X = 1:num.pc,
      FUN = function(x) {
        return(unlist(x = lapply(
          X = object@real.fval[, x],
          FUN = empP,
          nullval = object@fake.fval[, x]
        )))
      }
    )
  )
  rownames(x = object@emp.pval) <- pc.genes
  colnames(x = object@emp.pval) <- paste0("PC", 1:ncol(x = object@emp.pval),)
  return(object)
}

# Documentatin
##############
#' @export
#'
JackRandom <- function(
  scaled.data,
  prop.use = 0.01,
  r1.use = 1,
  r2.use = 5,
  seed.use = 1,
  rev.pca = FALSE,
  do.fast = FALSE
) {
  set.seed(seed = seed.use)
  rand.genes <- sample(
    x = rownames(x = scaled.data),
    size = nrow(x = scaled.data) * prop.use
  )
  # make sure that rand.genes is at least 3
  if (length(x = rand.genes) < 3){
    rand.genes <- sample(x = rownames(x = scaled.data), size = 3)
  }
  data.mod <- scaled.data
  data.mod[rand.genes, ] <- shuffleMatRow(x = scaled.data[rand.genes, ])
  if (rev.pca){
    if (do.fast) {
      fake.pca <- irlba(A = data.mod, nv = r2.use)
      fake.rot <- fake.pca$v[, 1:r2.use]
      rownames(x = fake.rot) <- colnames(x = data.mod)
      colnames(x = fake.rot) <- paste0("PC", 1:r2.use)
      fake.x <- fake.pca$u[, 1:r2.use]
      rownames(x = fake.x) <- rownames(x = data.mod)
      colnames(x = fake.x) <- colnames(x = fake.rot)
    } else {
      fake.pca <- prcomp(x = data.mod)
      fake.x <- fake.pca$x
      fake.rot <- fake.pca$rotation
    }
  } else {
    data.mod <- t(x = data.mod)
    if (do.fast) {
      fake.pca <- irlba(A = data.mod, nv = r2.use)
      fake.rot <- fake.pca$u[, 1:r2.use]
      rownames(x = fake.rot) <- rownames(x = data.mod)
      colnames(x = fake.rot) <- paste0("PC", 1:r2.use)
      fake.x <- fake.pca$v[, 1:r2.use]
      rownames(x = fake.x) <- colnames(x = data.mod)
      colnames(x = fake.x) <- colnames(x = fake.rot)
    } else {
      fake.pca <- prcomp(x = data.mod)
      fake.x <- fake.pca$rotation
      fake.rot <- fake.pca$x
    }
  }
  return(fake.x[rand.genes, r1.use:r2.use])
}
