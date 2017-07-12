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
  num.pc = 20,
  num.replicate = 100,
  prop.freq = 0.01,
  do.print = FALSE
) {
  if (is.null(object@dr$pca)) {
    stop("PCA has not been computed yet. Please run PCA().")
  }
  # error checking for number of PCs
  if (num.pc > ncol(x = object@dr$pca@cell.embeddings)) {
    num.pc <- ncol(x = object@dr$pca@cell.embeddings)
    warning("Number of PCs specified is greater than PCs available. Setting num.pc to ", num.pc, " and continuing.")
  }
  if (num.pc > length(x = object@cell.names)) {
    num.pc <- length(x = object@cell.names)
    warning("Number of PCs specified is greater than number of cells. Setting num.pc to ", num.pc, " and continuing.")
  }
  pc.genes <- rownames(x = object@dr$pca@gene.loadings)
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
  md.x <- as.matrix(x = object@dr$pca@gene.loadings)
  md.rot <- as.matrix(x = object@dr$pca@cell.embeddings)
  if (do.print) {
    applyFunction <- pbsapply
  } else {
    applyFunction <- sapply
  }
  rev.pca <- GetCalcParam(object = object,
                          calculation = "PCA",
                          parameter = "rev.pca")
  scale.by.varexp <- GetCalcParam(object = object,
                                  calculation = "PCA",
                                  parameter = "scale.by.varexp")
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
        scale.by.varexp
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
  scale.by.varexp = scale.by.varexp
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
  temp.object <- new("seurat")
  temp.object@scale.data <- data.mod
  temp.object <- PCA(temp.object, pcs.compute = r2.use, pc.genes = rownames(data.mod),
                     rev.pca = rev.pca, scale.by.varexp = scale.by.varexp,
                     do.print = F)
  fake.x <- PCALoad(temp.object)
  fake.rot <- PCAEmbed(temp.object)
  return(fake.x[rand.genes, r1.use:r2.use])
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
  genes.use <- SetIfNull(x = genes.use, default = rownames(x = object@pca.x))
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


#internal
jackStrawF <- function(prop = 0.1, myR1, myR2 = 3, data = smD) {
  randGenes <- sample(x = rownames(x = data), size = nrow(x = data) * prop)
  smD.mod <- data
  smD.mod[randGenes, ] <- shuffleMatRow(x = data[randGenes, ])
  fmd.pca <- prcomp(x = smD.mod)
  fmd.x <- fmd.pca$x
  fmd.rot <- fmd.pca$rotation
  fakeF <- unlist(x = lapply(
    X = randGenes,
    FUN = jackF,
    r1 = myR1,
    r2 = myR2,
    x = fmd.x,
    rot = fmd.rot
  ))
}

#internal
jackF <- function(gene, r1 = 1,r2 = 2, x = md.x, rot = md.rot) {
  if (r2 == 1) { #assuming r1, r2=1
    mod.x <- x[, r1]
    mod.x[gene] <- 0
    return(var.test(
      x = (x[, r1] %*% t(x = rot[, r1])),
      y = (mod.x %*% t(x = rot[, r1]))
    )$statistic)
  }
  mod.x <- x[, 1:r2]
  mod.x[gene, r1:r2] <- rep(x = 0, r2 - r1 + 1)
  return(var.test(
    x = (x[, 1:r2] %*% t(x = rot[, 1:r2])),
    y = (mod.x[, 1:r2] %*% t(x = rot[, 1:r2]))
  )$statistic)
}

#internal
empP <- function(x, nullval) {
  return(sum(nullval > x) / length(x = nullval))
}
