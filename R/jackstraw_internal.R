#define class to store jackstraw data
jackstraw.data <- setClass(
  Class = "jackstraw.data",
  slots = list(
    emperical.p.value = "matrix",
    fake.pc.scores = "matrix",
    emperical.p.value.full = "matrix"
  )
)

#internal
#
#' @importFrom stats prcomp
#
JackstrawF <- function(prop = 0.1, myR1, myR2 = 3, data = smD) {
  randGenes <- sample(x = rownames(x = data), size = nrow(x = data) * prop)
  smD.mod <- data
  smD.mod[randGenes, ] <- MatrixRowShuffle(x = data[randGenes, ])
  fmd.pca <- prcomp(x = smD.mod)
  fmd.x <- fmd.pca$x
  fmd.rot <- fmd.pca$rotation
  fakeF <- unlist(x = lapply(
    X = randGenes,
    FUN = JackF,
    r1 = myR1,
    r2 = myR2,
    x = fmd.x,
    rot = fmd.rot
  ))
}

#internal
#
#' @importFrom stats var.test
#
JackF <- function(gene, r1 = 1,r2 = 2, x = md.x, rot = md.rot) {
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
EmpiricalP <- function(x, nullval) {
  return(sum(nullval > x) / length(x = nullval))
}

#internal
JackRandom <- function(
  scaled.data,
  prop.use = 0.01,
  r1.use = 1,
  r2.use = 5,
  seed.use = 1,
  rev.pca = FALSE,
  weight.by.var = weight.by.var
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
  data.mod[rand.genes, ] <- MatrixRowShuffle(x = scaled.data[rand.genes, ])
  temp.object <- new("seurat")
  temp.object@cell.names <- colnames(x = data.mod)
  temp.object@scale.data <- data.mod
  temp.object <- RunPCA(
    object = temp.object,
    pcs.compute = r2.use,
    pc.genes = rownames(x = data.mod),
    rev.pca = rev.pca,
    weight.by.var = weight.by.var,
    do.print = FALSE
  )
  fake.x <- PCALoad(object = temp.object)
  fake.rot <- PCAEmbed(object = temp.object)
  return(fake.x[rand.genes, r1.use:r2.use])
}
