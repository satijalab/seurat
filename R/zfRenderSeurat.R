#' Zebrafish analysis functions
#'
#' These functions are for the zebrafish analysis from Satija et al. 2015
#' @rdname zfRenderSeurat
#' @name zfRenderSeurat
#' @aliases zf.cells.render zf.anchor.render zf.insitu.vec.lateral zf.insitu.lateral zf.insitu.dorsal zf.insitu.ventral zf.insitu.side
#'
#' @references Satija R*, Farrell JA*, Gennert D, Schier AF, Regev A.
#' Spatial reconstruction of single-cell gene expression data.
#' Nature Biotechnology. 2015
#'
#' @importFrom grDevices colorRampPalette
#'


# Draw 3D in situ predictions from Zebrafish dataset
#
# From Jeff Farrell
#
# @param data Predicted expression levels across Zebrafish bins
# @param label Plot label
# @param ... Extra parameters
#
# @export
situ3d <- function(data, label = NULL, ...) {
  # Call Seurat function to get the in situ values out.
  exp.1 <- data
  exp.1 <- (exp.1 - min(exp.1)) / (max(exp.1) - min(exp.1))
  # Reformat them into an expression matrix as expected by the plotting function
  expression.matrix <- data.frame(matrix(data = exp.1, nrow = 8, ncol = 8))
  rownames(x = expression.matrix) <- c(
    "24-30",
    "17-23",
    "13-16",
    "9-12",
    "7-8",
    "5-6",
    "3-4",
    "1-2"
  )
  names(x = expression.matrix) <- c(
    "1-4",
    "5-8",
    "9-12",
    "13-16",
    "17-20",
    "21-24",
    "25-28",
    "29-32"
  )
  # Call the plotting function.
  zf.insitu.side(expression.matrix)
  par3d(windowRect = c(0, 0, 800, 800))
  # Label or not and then set the view.
  if (! is.null(x = label)) {
    text3d(x = 0, y = 0, z = 1.5, text = label, cex = 3)
  }
  view3d(zoom = .75, theta = 0, phi = -90, fov = 0)
}


#' @export
zf.cells.render <- function(
  seuratObject,
  cells.use,
  do.rotate = TRUE,
  label = TRUE,
  calc.new = FALSE,
  col.use = "red",
  radius.use = 0.05,
  col.prob = FALSE,
  do.new = TRUE,
  ...
) {
  tierBins <- 30 # 1 bin per cell tier.
  DVBins <- 64 # 1 bin every 5.625 degrees; compatible with our current 8-bin system.
  phiPerTier <- pi / (-2 * tierBins)
  thetaPerDV <- (2 * pi) / DVBins
  if (length(x = col.use) == 1) {
    col.use <- rep(x = col.use, length(x = cells.use))
  }
  # Reformat that probability into an expression matrix as expected by the plotting function
  if (col.prob) {
    prob.matrix <- data.frame(matrix(
      data = apply(
        X = seuratObject@final.prob[, cells.use],
        MARGIN = 1,
        FUN = sum
      ),
      nrow = 8,
      ncol = 8
    ))
  } else {
    prob.matrix <- data.frame(matrix(data = 0, nrow = 8, ncol = 8))
  }
  rownames(x = prob.matrix) <- c(
    "24-30",
    "17-23",
    "13-16",
    "9-12",
    "7-8",
    "5-6",
    "3-4",
    "1-2"
  )
  names(x = prob.matrix) <- c(
    "1-4",
    "5-8",
    "9-12",
    "13-16",
    "17-20",
    "21-24",
    "25-28",
    "29-32"
  )
  # Call the plotting function.
  if (do.new)
  {
    zf.insitu.side(
      expressionMatrix = prob.matrix,
      mirror = TRUE,
      nonmirror = FALSE
    )
  }
  i <- 1
  for(cell.use in cells.use) {
    #add the centroid
    anchor.centroid <- ExactCellCentroid(
      cell.probs = seuratObject@final.prob[, cell.use]
    )
    tiers.min <- c(30, 24, 16, 12, 8, 6, 4, 2, 0)
    tiers.size <- diff(x = tiers.min)
    anchor.dorsality <- DVBins - ((anchor.centroid[2] - 1) / 7) * DVBins / 2
    anchor.tier.bin <- anchor.centroid[1]
    anchor.tier.bin <- anchor.centroid[1]
    anchor.tier.floor <- floor(x = anchor.tier.bin)
    anchor.tier.left <- anchor.tier.bin - anchor.tier.floor
    anchor.tier <- tiers.min[anchor.tier.bin] +
      (tiers.size[anchor.tier.floor] * anchor.tier.left)
    x1 <- cos(x = pi - thetaPerDV * anchor.dorsality) *
      sin(x = 0.5 * pi + phiPerTier * anchor.tier)
    y1 <- sin(x = pi - thetaPerDV * anchor.dorsality) *
      sin(0.5 * pi + phiPerTier * anchor.tier)
    z1 <- cos(x = 0.5 * pi + phiPerTier * anchor.tier)
    spheres3d(
      x = x1,
      y = y1,
      z = z1,
      radius = radius.use,
      color = col.use[i],
      alpha = .65,
      lit = FALSE
    )
    i <- i + 1;
  }
  view3d(zoom = .75, theta = 0, phi = -90, fov = 0)
  # Format the plot
  if (do.new) {
    if (label) {
      # text3d(x=0, y=0, z=1.5, text=paste(this.anchor, anchor.distance),cex=3)
      text3d(x = 1.4, y = 0, z = -0.3, cex = 2.25, text = "Dor")
      text3d(x = -1.4, y = 0, z = -0.3, cex = 2.25, text = "Ven")
      text3d(x = 0, y = 0, z = 1.2, cex = 2.25, text = "An")
      text3d(x = 0, y = 0, z = -1.2, cex = 2.25, text = "Veg")
    }
    #view3d(zoom=.75, theta=0, phi=-90, fov=0) # This makes you look at dorsality 48.
    to.rotate <- (anchor.dorsality - 47.5) / 64
    if (to.rotate < 0) {
      to.rotate <- 1 + to.rotate
    }
    if (do.rotate) {
      play3d(spin3d(axis = c(0, 0, 1), rpm = 60), duration = to.rotate)
    }
    par3d(windowRect = c(0, 0, 800, 800))
  }
}

#' @export
zf.anchor.render <- function(
  seuratObject,
  this.anchor,
  anchors,
  label = TRUE,
  do.rotate = TRUE,
  calc.new = FALSE,
  ...
) {
  # Determine geometry
  tierBins <- 30 # 1 bin per cell tier.
  DVBins <- 64 # 1 bin every 5.625 degrees; compatible with our current 8-bin system.
  phiPerTier <- pi / (-2 * tierBins)
  thetaPerDV <- 2 * pi / DVBins
  cellColor <- "#EB008B"
  centroidColor <- "green"
  # Get probability of anchor being in each bin.
  if (calc.new) {
    anchor.prob <- data.frame(
      prob = as.numeric(x = project.cell(
        seuratObject,
        this.anchor,
        do.plot = FALSE,
        safe = FALSE
      ))
    )
  }
  if (!(calc.new)) {
    anchor.prob.raw <- data.frame(prob = seuratObject@final.prob[, this.anchor])
  }
  # Normalize so that strongest probability is 1.
  anchor.prob <- round(x = anchor.prob.raw / max(anchor.prob.raw), digits = 5)
  # Reformat that probability into an expression matrix as expected by the plotting function
  prob.matrix <- data.frame(matrix(data = anchor.prob[, 1], nrow = 8, ncol = 8))
  rownames(x = prob.matrix) <- c(
    "24-30",
    "17-23",
    "13-16",
    "9-12",
    "7-8",
    "5-6",
    "3-4",
    "1-2"
  )
  names(x = prob.matrix) <- c(
    "1-4",
    "5-8",
    "9-12",
    "13-16",
    "17-20",
    "21-24",
    "25-28",
    "29-32"
  )
  # Call the plotting function.
  zf.insitu.side(expressionMatrix = prob.matrix, mirror = TRUE, nonmirror = FALSE)
  # Add the anchor cell
  anchor.dorsality <- DVBins - anchors[this.anchor, "dorsality"] * DVBins / 2
  anchor.tier <- anchors[this.anchor, "tier"]
  x <- cos(x = pi - thetaPerDV * anchor.dorsality) *
    sin(x = 0.5 * pi + phiPerTier * anchor.tier)
  y <- sin(x = pi - thetaPerDV * anchor.dorsality) *
    sin(x = 0.5 * pi + phiPerTier * anchor.tier)
  z <- cos(x = 0.5 * pi + phiPerTier * anchor.tier)
  spheres3d(
    x = x,
    y = y,
    z = z,
    radius = 0.16,
    color = cellColor,
    alpha = .65,
    lit = FALSE
  )
  #add the centroid
  anchor.centroid <- ExactCellCentroid(cell.probs = anchor.prob.raw)
  tiers.min <- c(30, 24, 16, 12, 8, 6, 2, 1)
  tiers.size <- diff(x = tiers.min)
  anchor.dorsality <- DVBins - ((anchor.centroid[2] - 1) / 7) * DVBins / 2
  anchor.tier.bin <- anchor.centroid[1]
  anchor.tier.bin <- anchor.centroid[1]
  anchor.tier.floor <- floor(x = anchor.tier.bin)
  anchor.tier.left <- anchor.tier.bin - anchor.tier.floor
  anchor.tier <- tiers.min[anchor.tier.bin] +
    (tiers.size[anchor.tier.floor] * anchor.tier.left)
  x1 <- cos(x = pi - thetaPerDV * anchor.dorsality) *
    sin(x = 0.5 * pi + phiPerTier * anchor.tier)
  y1 <- sin(x = pi - thetaPerDV * anchor.dorsality) *
    sin(x = 0.5 * pi + phiPerTier * anchor.tier)
  z1 <- cos(x = 0.5 * pi + phiPerTier * anchor.tier)
  spheres3d(
    x = x1,
    y = y1,
    z = z1,
    radius = 0.08,
    color = centroidColor,
    alpha  =.65,
    lit = FALSE
  )
  #anchor.tier.true.bin=anchors[this.anchor,"tier.bin"]
  #anchor.tier.true <- anchors[this.anchor, "tier"]
  #tier.min.distance=anchor.tier.true-tiers.min[anchor.tier.true.bin]
  #if (tier.min.distance > 0) {
  #  anchor.tier.true.bin=anchors[this.anchor,"tier.bin"]-tier.min.distance/(tiers.size[anchor.tier.true.bin])
  #}
  anchor.distance <- round(x = dist(x = rbind(c(x, y), c(x1, y1))), digits = 2)
  # Format the plot
  if (label) {
    text3d(
      x = 0,
      y = 0,
      z = 1.5,
      text = paste(this.anchor, anchor.distance),
      cex=3
    )
    text3d(x = 1.4, y = 0, z = -0.3, cex = 2.25, text = "Dor")
    text3d(x = -1.4, y = 0, z = -0.3, cex = 2.25, text = "Ven")
    text3d(x = 0, y = 0, z = 1.2, cex = 2.25, text = "An")
    text3d(x = 0, y = 0, z = -1.2, cex = 2.25, text = "Veg")
  }
  view3d(zoom = .75, theta = 0, phi = -90, fov = 0) # This makes you look at dorsality 48.
  to.rotate <- (anchor.dorsality - 47.5) / 64
  if (to.rotate < 0) {
    to.rotate <- 1 + to.rotate
  }
  if (do.rotate) {
    play3d(spin3d(axis = c(0, 0, 1), rpm = 60), duration = to.rotate)
  }
  par3d(windowRect = c(0, 0, 800, 800))
}

zf.anchor.map <- function(
  seuratObject,
  this.anchor,
  anchors,
  calc.new = FALSE,
  ...
) {
  # Determine geometry
  if (calc.new) {
    anchor.prob <- (prob = as.numeric(x = project.cell(
      seuratObject,
      this.anchor,
      do.plot = FALSE,
      safe = FALSE
    )))
  }
  if (!(calc.new)) {
    anchor.prob.raw <- (prob = seuratObject@final.prob[, this.anchor])
  }
  # Normalize so that strongest probability is 1.
  anchor.prob <- round(x = anchor.prob.raw / max(anchor.prob.raw), digits = 5)
  hm4(
    matrix(data = as.numeric(x = anchor.prob), nrow = 8),
    trace = "none",
    Rowv = NA,
    Colv = NA
  )
  text(
    x = anchors[this.anchor, "dv.bin"],
    y = 9-anchors[this.anchor, "tier.bin"],
    labels = "X",
    cex = 1.5
  )
  text(x = 5, y = 1, labels = this.anchor)
  # Reformat that probability into an expression matrix as expected by the plotting function
}

#' @export
zf.insitu.vec.lateral <- function(
  expression.vector,
  label = TRUE,
  title = NULL,
  ...
) {
  # Reformat them into an expression matrix as expected by the plotting function
  expression.matrix <- data.frame(matrix(
    data = expression.vector,
    nrow = 8,
    ncol = 8
  ))
  rownames(x = expression.matrix) <- c(
    "24-30",
    "17-23",
    "13-16",
    "9-12",
    "7-8",
    "5-6",
    "3-4",
    "1-2"
  )
  names(x = expression.matrix) <- c(
    "1-4",
    "5-8",
    "9-12",
    "13-16",
    "17-20",
    "21-24",
    "25-28",
    "29-32"
  )
  # Call the plotting function.
  zf.insitu.side(expressionMatrix = expression.matrix, ...)
  par3d(windowRect = c(0, 0, 800, 800))
  # Label or not and then set the view.
  if (! is.null(x = title) & ! label) {
    text3d(x = 0, y = 0, z = 1.2, text = title, cex = 4.5)
  }
  if (! is.null(x = title) & label) {
    text3d(x = 0, y = 0, z = 1.5, text = title, cex = 3)
  }
  if (label) {
    text3d(x = 1.4, y = 0, z = -0.3, cex = 2.25, text = "Dor")
    text3d(x = -1.4, y = 0, z = -0.3, cex = 2.25, text = "Ven")
    text3d(x = 0, y = 0, z = 1.2, cex = 2.25, text = "An")
    text3d(x = 0, y = 0, z = -1.2, cex = 2.25, text = "Veg")
  }
  view3d(zoom = .75, theta = 0, phi = -90, fov = 0)
}

#' @export
zf.insitu.lateral <- function(seuratObject, gene, label = TRUE, ...) {
  # Call Seurat function to get the in situ values out.
  expression <- CalcInsitu(
    seuratObject,
    gene,
    do.plot = FALSE,
    do.return = TRUE,
    do.norm = TRUE,
    ...
  )
  # Reformat them into an expression matrix as expected by the plotting function
  expression.matrix <- data.frame(matrix(data = expression, nrow = 8, ncol = 8))
  rownames(x = expression.matrix) <- c(
    "24-30",
    "17-23",
    "13-16",
    "9-12",
    "7-8",
    "5-6",
    "3-4",
    "1-2"
  )
  names(x = expression.matrix) <- c(
    "1-4",
    "5-8",
    "9-12",
    "13-16",
    "17-20",
    "21-24",
    "25-28",
    "29-32"
  )
  # Call the plotting function.
  zf.insitu.side(expressionMatrix = expression.matrix)
  par3d(windowRect = c(0, 0, 800, 800))
  # Label or not and then set the view.
  text3d(x = 0, y = 0, z = 1.5, text = gene, cex = 3)
  if (label) {
    text3d(x = 1.4, y = 0, z = -0.3, cex = 2.25, text = "Dor")
    text3d(x = -1.4, y = 0, z = -0.3, cex = 2.25, text = "Ven")
    text3d(x = 0, y = 0, z = 1.2, cex = 2.25, text = "An")
    text3d(x = 0, y = 0, z = -1.2, cex = 2.25, text = "Veg")
  }
  view3d(zoom = .75, theta = 0, phi = -90, fov = 0)
}

#' @export
zf.insitu.dorsal <- function(seuratObject, gene, label=TRUE, ...) {
  # Call Seurat function to get the in situ values out.
  expression <- CalcInsitu(
    seuratObject,
    gene,
    do.plot = FALSE,
    do.return = TRUE,
    do.norm = TRUE,
    ...
  )
  # Reformat them into an expression matrix as expected by the plotting function
  expression.matrix <- data.frame(matrix(data = expression, nrow = 8, ncol = 8))
  rownames(x = expression.matrix) <- c(
    "24-30",
    "17-23",
    "13-16",
    "9-12",
    "7-8",
    "5-6",
    "3-4",
    "1-2"
  )
  names(x = expression.matrix) <- c(
    "1-4",
    "5-8",
    "9-12",
    "13-16",
    "17-20",
    "21-24",
    "25-28",
    "29-32"
  )
  # Call the plotting function.
  zf.insitu.side(expressionMatrix = expression.matrix)
  par3d(windowRect = c(0, 0, 800, 800))
  # Label or not and then set the view.
  if (label) {
    text3d(x = 0, y = 0, z = 1.5, text = gene, cex = 3)
    text3d(x = 1.4, y = 0, z = -0.3, cex = 2.25, text = "Dor")
    text3d(x = -1.4, y = 0, z = -0.3, cex = 2.25, text = "Ven")
    text3d(x = 0, y = 0, z = 1.2, cex = 2.25, text = "An")
    text3d(x = 0, y = 0, z = -1.2, cex = 2.25, text = "Veg")
  }
  rotMat <- rotationMatrix(-pi / 2, 0, 0, 1)  %*% rotationMatrix(-pi / 2, 0, 1, 0)
  view3d(zoom = .75, userMatrix = rotMat, fov = 0)
}

#' @export
zf.insitu.ventral <- function(seuratObject, gene, label=TRUE, ...) {
  # Call Seurat function to get the in situ values out.
  expression <- CalcInsitu(
    seuratObject,
    gene,
    do.plot = FALSE,
    do.return = TRUE,
    do.norm = TRUE,
    ...
  )
  # Reformat them into an expression matrix as expected by the plotting function
  expression.matrix <- data.frame(matrix(data = expression, nrow = 8, ncol = 8))
  rownames(x = expression.matrix) <- c(
    "24-30",
    "17-23",
    "13-16",
    "9-12",
    "7-8",
    "5-6",
    "3-4",
    "1-2"
  )
  names(x = expression.matrix) <- c(
    "1-4",
    "5-8",
    "9-12",
    "13-16",
    "17-20",
    "21-24",
    "25-28",
    "29-32"
  )
  # Call the plotting function.
  zf.insitu.side(expressionMatrix = expression.matrix)
  par3d(windowRect = c(0, 0, 800, 800))
  # Label or not and then set the view.
  if (label) {
    text3d(x = 0, y = 0, z = 1.5, text = gene, cex = 3)
    text3d(x = 1.4, y = 0, z = -0.3, cex = 2.25, text = "Dor")
    text3d(x = -1.4, y = 0, z = -0.3, cex = 2.25, text = "Ven")
    text3d(x = 0, y = 0, z = 1.2, cex = 2.25, text = "An")
    text3d(x = 0, y = 0, z = -1.2, cex = 2.25, text = "Veg")
  }
  rotMat <- rotationMatrix(pi / 2, 0, 0, 1) %*% rotationMatrix(pi / 2, 0, 1, 0)
  view3d(zoom = .75, userMatrix = rotMat, fov = 0)
}

# @importFrom gdata interleave
#' @importFrom utils installed.packages
#' @export
zf.insitu.side <- function(expressionMatrix, nonmirror = TRUE, mirror = TRUE) {
  if (!'gdata' %in% rownames(x = installed.packages())) {
    stop("Please install gdata")
  }
  # Determine geometry
  tierBins <- 30 # 1 bin per cell tier.
  DVBins <- 64 # 1 bin every 5.625 degrees; compatible with our current 8-bin system.
  phiPerTier <- pi / (-2 * tierBins)
  thetaPerDV <- 2 * pi / DVBins
  # Determine colors
  yolkColor <- "#FDF5E6"
  marginColor <- "#CDC8B1"
  insituPalette <- colorRampPalette(
    colors = c("#FDF5E6", "#483D8B"),
    space = "Lab"
  )
  insituColors <- insituPalette(51)
  # Make a dataframe that will hold the position of each quadrilateral for the drawing, default to yolk-colored.
  # Top of the embryo
  drawEmbryo <- data.frame(
    tier = rep(x = 1:tierBins, DVBins),
    DV = rep(x = 1:DVBins, each = tierBins),
    color = rep(x = yolkColor, tierBins * DVBins),
    stringsAsFactors = FALSE
  )
  # The yolk part
  drawEmbryo <- rbind(
    drawEmbryo,
    data.frame(
      tier = rep(x = -tierBins:-1, DVBins),
      DV = rep(x = 1:DVBins, each = tierBins),
      color = rep(x = yolkColor, tierBins * DVBins),
      stringsAsFactors = FALSE
    )
  )
  # Add the margin
  drawEmbryo <- rbind(
    drawEmbryo,
    data.frame(
      tier = rep(x = 0, DVBins),
      DV = 1:DVBins,
      color = rep(x = marginColor, DVBins),
      stringsAsFactors = FALSE
    )
  )
  # Determine the 4 coordinates for each quadrilateral defined by a bin
  drawEmbryo$x1 <- cos(x = pi - thetaPerDV * drawEmbryo$DV) *
    sin(x = 0.5 * pi + phiPerTier * (drawEmbryo$tier - 1))
  drawEmbryo$x2 <- cos(x = pi - thetaPerDV * (drawEmbryo$DV - 1)) *
    sin(x = 0.5 * pi + phiPerTier * (drawEmbryo$tier - 1))
  drawEmbryo$x3 <- cos(x = pi - thetaPerDV * (drawEmbryo$DV - 1)) *
    sin(x = 0.5 * pi + phiPerTier * drawEmbryo$tier)
  drawEmbryo$x4 <- cos(x = pi - thetaPerDV * drawEmbryo$DV) *
    sin(x = 0.5 * pi + phiPerTier * drawEmbryo$tier)
  drawEmbryo$y1 <- sin(x = pi - thetaPerDV * drawEmbryo$DV) *
    sin(x = 0.5 * pi + phiPerTier * (drawEmbryo$tier - 1))
  drawEmbryo$y2 <- sin(x = pi - thetaPerDV * (drawEmbryo$DV - 1)) *
    sin(x = 0.5 * pi + phiPerTier * (drawEmbryo$tier - 1))
  drawEmbryo$y3 <- sin(x = pi - thetaPerDV * (drawEmbryo$DV - 1)) *
    sin(0.5 * pi + phiPerTier * drawEmbryo$tier)
  drawEmbryo$y4 <- sin(x = pi - thetaPerDV * drawEmbryo$DV) *
    sin(x = 0.5 * pi + phiPerTier * drawEmbryo$tier)
  drawEmbryo$z1 <- cos(x = 0.5 * pi + phiPerTier * (drawEmbryo$tier - 1))
  drawEmbryo$z2 <- cos(x = 0.5 * pi + phiPerTier * (drawEmbryo$tier - 1))
  drawEmbryo$z3 <- cos(x = 0.5 * pi + phiPerTier * drawEmbryo$tier)
  drawEmbryo$z4 <- cos(x = 0.5 * pi + phiPerTier * drawEmbryo$tier)
  # Now, reassign the color for each of the bins that has expression >0.
  for (tier in 1:dim(x = expressionMatrix)[1]) {
    for (DV in 1:dim(x = expressionMatrix)[2]) {
      if (! expressionMatrix[tier, DV] == 0 ) {
        # Figure out limits of the bins desired from the names of the row & col of this table cell
        tierLimits <- as.numeric(x = unlist(x = strsplit(
          x = row.names(x = expressionMatrix)[tier],
          split = "-"
        )))
        DVLimits <- as.numeric(x = unlist(x = strsplit(
          x = names(x = expressionMatrix)[DV],
          split = "-"
        )))
        # Figure out the value for this color.
        thisColor <- insituColors[(floor(
          x = as.numeric(x = expressionMatrix[tier, DV]) * 50
        )) + 1]
        # Loop through and assign the color to every bin in the limits
        for (thisTier in min(tierLimits):max(tierLimits)) {
          if (nonmirror) {
            for (thisDV in min(DVLimits):max(DVLimits)) {
              thisRow <- (thisDV - 1) * tierBins + thisTier
              drawEmbryo[thisRow, ]$color <- thisColor
            }
          }
          # If mirror is on, also assign the other side of the embryo.
          if (mirror) {
            for (thisDV in (DVBins - max(DVLimits) + 1):(DVBins - min(DVLimits) + 1)) {
              thisRow <- (thisDV - 1) * tierBins + thisTier
              drawEmbryo[thisRow, ]$color <- thisColor
            }
          }
        }

      }
    }
  }
  # Take the coordinates and reformat the lists to pass to RGL
  quadX <- gdata::interleave(
    drawEmbryo$x1,
    drawEmbryo$x2,
    drawEmbryo$x3,
    drawEmbryo$x4,
    drop = TRUE
  )
  dim(x = quadX) <- c(dim(x = quadX)[1] * dim(x = quadX)[2], 1)
  quadY <- gdata::interleave(
    drawEmbryo$y1,
    drawEmbryo$y2,
    drawEmbryo$y3,
    drawEmbryo$y4,
    drop = TRUE
  )
  dim(x = quadY) <- c(dim(x = quadY)[1] * dim(x = quadY)[2], 1)
  quadZ <- gdata::interleave(
    drawEmbryo$z1,
    drawEmbryo$z2,
    drawEmbryo$z3,
    drawEmbryo$z4,
    drop = TRUE
  )
  dim(x = quadZ) <- c(dim(x = quadZ)[1] * dim(x = quadZ)[2], 1)
  quadColor <- rep(x = drawEmbryo$color, each = 4)
  # Initialize an RGL view
  open3d()
  # Call quads to plot the embryo.
  quads3d(
    x = quadX,
    y = quadY,
    z = quadZ,
    color = quadColor,
    alpha = 1,
    lit = FALSE
  )
}

#used for zebrafish plotting
vp.layout <- function(x, y) {
  viewport(layout.pos.row = x, layout.pos.col = y)
}
