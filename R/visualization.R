#' @importFrom utils globalVariables
#' @importFrom ggplot2 ggproto GeomViolin
#'
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Heatmaps
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Dimensional reduction heatmap
#'
#' Draws a heatmap focusing on a principal component. Both cells and genes are sorted by their
#' principal component scores. Allows for nice visualization of sources of heterogeneity in the dataset.
#'
#' @inheritParams DoHeatmap
#' @param dims Dimensions to plot
#' @param nfeatures Number of genes to plot
#' @param cells A list of cells to plot. If numeric, just plots the top cells.
#' @param reduction Which dimmensional reduction to use
#' @param balanced Plot an equal number of genes with both + and - scores.
#' @param projected Use the full projected dimensional reduction
#' @param ncol Number of columns to plot
#' @param fast If true, use \code{image} to generate plots; faster than using ggplot2, but not customizable
#' @param assays A vector of assays to pull data from
#' @param combine Combine plots into a single \code{\link[patchwork]{patchwork}ed}
#' ggplot object. If \code{FALSE}, return a list of ggplot objects
#'
#' @return No return value by default. If using fast = FALSE, will return a
#' \code{\link[patchwork]{patchwork}ed} ggplot object if combine = TRUE, otherwise
#' returns a list of ggplot objects
#'
#' @importFrom patchwork wrap_plots
#' @export
#'
#' @seealso \code{\link[graphics]{image}} \code{\link[ggplot2]{geom_raster}}
#'
#' @examples
#' DimHeatmap(object = pbmc_small)
#'
DimHeatmap <- function(
  object,
  dims = 1,
  nfeatures = 30,
  cells = NULL,
  reduction = 'pca',
  disp.min = -2.5,
  disp.max = NULL,
  balanced = TRUE,
  projected = FALSE,
  ncol = NULL,
  fast = TRUE,
  raster = TRUE,
  slot = 'scale.data',
  assays = NULL,
  combine = TRUE
) {
  ncol <- ncol %||% ifelse(test = length(x = dims) > 2, yes = 3, no = length(x = dims))
  plots <- vector(mode = 'list', length = length(x = dims))
  assays <- assays %||% DefaultAssay(object = object)
  disp.max <- disp.max %||% ifelse(
    test = slot == 'scale.data',
    yes = 2.5,
    no = 6
  )
  if (!DefaultAssay(object = object[[reduction]]) %in% assays) {
    warning("The original assay that the reduction was computed on is different than the assay specified")
  }
  cells <- cells %||% ncol(x = object)
  if (is.numeric(x = cells)) {
    cells <- lapply(
      X = dims,
      FUN = function(x) {
        cells <- TopCells(
          object = object[[reduction]],
          dim = x,
          ncells = cells,
          balanced = balanced
        )
        if (balanced) {
          cells$negative <- rev(x = cells$negative)
        }
        cells <- unlist(x = unname(obj = cells))
        return(cells)
      }
    )
  }
  if (!is.list(x = cells)) {
    cells <- lapply(X = 1:length(x = dims), FUN = function(x) {return(cells)})
  }
  features <- lapply(
    X = dims,
    FUN = TopFeatures,
    object = object[[reduction]],
    nfeatures = nfeatures,
    balanced = balanced,
    projected = projected
  )
  features.all <- unique(x = unlist(x = features))
  if (length(x = assays) > 1) {
    features.keyed <- lapply(
      X = assays,
      FUN = function(assay) {
        features <- features.all[features.all %in% rownames(x = object[[assay]])]
        if (length(x = features) > 0) {
          return(paste0(Key(object = object[[assay]]), features))
        }
      }
    )
    features.keyed <- Filter(f = Negate(f = is.null), x = features.keyed)
    features.keyed <- unlist(x = features.keyed)
  } else {
    features.keyed <- features.all
    DefaultAssay(object = object) <- assays
  }
  data.all <- FetchData(
    object = object,
    vars = features.keyed,
    cells = unique(x = unlist(x = cells)),
    slot = slot
  )
  data.all <- MinMax(data = data.all, min = disp.min, max = disp.max)
  data.limits <- c(min(data.all), max(data.all))
  # if (check.plot && any(c(length(x = features.keyed), length(x = cells[[1]])) > 700)) {
  #   choice <- menu(c("Continue with plotting", "Quit"), title = "Plot(s) requested will likely take a while to plot.")
  #   if (choice != 1) {
  #     return(invisible(x = NULL))
  #   }
  # }
  if (fast) {
    nrow <- floor(x = length(x = dims) / 3.01) + 1
    orig.par <- par()$mfrow
    par(mfrow = c(nrow, ncol))
  }
  for (i in 1:length(x = dims)) {
    dim.features <- c(features[[i]][[2]], rev(x = features[[i]][[1]]))
    dim.features <- rev(x = unlist(x = lapply(
      X = dim.features,
      FUN = function(feat) {
        return(grep(pattern = paste0(feat, '$'), x = features.keyed, value = TRUE))
      }
    )))
    dim.cells <- cells[[i]]
    data.plot <- data.all[dim.cells, dim.features]
    if (fast) {
      SingleImageMap(
        data = data.plot,
        title = paste0(Key(object = object[[reduction]]), dims[i]),
        order = dim.cells
      )
    } else {
      plots[[i]] <- SingleRasterMap(
        data = data.plot,
        raster = raster,
        limits = data.limits,
        cell.order = dim.cells,
        feature.order = dim.features
      )
    }
  }
  if (fast) {
    par(mfrow = orig.par)
    return(invisible(x = NULL))
  }
  if (combine) {
    plots <- wrap_plots(plots, ncol = ncol, guides = "collect")
  }
  return(plots)
}

#' Feature expression heatmap
#'
#' Draws a heatmap of single cell feature expression.
#'
#' @param object Seurat object
#' @param features A vector of features to plot, defaults to \code{VariableFeatures(object = object)}
#' @param cells A vector of cells to plot
#' @param disp.min Minimum display value (all values below are clipped)
#' @param disp.max Maximum display value (all values above are clipped); defaults to 2.5
#' if \code{slot} is 'scale.data', 6 otherwise
#' @param group.by A vector of variables to group cells by; pass 'ident' to group by cell identity classes
#' @param group.bar Add a color bar showing group status for cells
#' @param group.colors Colors to use for the color bar
#' @param slot Data slot to use, choose from 'raw.data', 'data', or 'scale.data'
#' @param assay Assay to pull from
# @param check.plot Check that plotting will finish in a reasonable amount of time
#' @param label Label the cell identies above the color bar
#' @param size Size of text above color bar
#' @param hjust Horizontal justification of text above color bar
#' @param angle Angle of text above color bar
#' @param raster If true, plot with geom_raster, else use geom_tile. geom_raster may look blurry on
#' some viewing applications such as Preview due to how the raster is interpolated. Set this to FALSE
#' if you are encountering that issue (note that plots may take longer to produce/render).
#' @param draw.lines Include white lines to separate the groups
#' @param lines.width Integer number to adjust the width of the separating white lines.
#' Corresponds to the number of "cells" between each group.
#' @param group.bar.height Scale the height of the color bar
#' @param combine Combine plots into a single \code{\link[patchwork]{patchwork}ed}
#' ggplot object. If \code{FALSE}, return a list of ggplot objects
#'
#' @return A \code{\link[patchwork]{patchwork}ed} ggplot object if
#' \code{combine = TRUE}; otherwise, a list of ggplot objects
#'
#' @importFrom stats median
#' @importFrom scales hue_pal
#' @importFrom ggplot2 annotation_raster coord_cartesian scale_color_manual
#' ggplot_build aes_string
#' @importFrom patchwork wrap_plots
#' @export
#'
#' @examples
#' DoHeatmap(object = pbmc_small)
#'
DoHeatmap <- function(
  object,
  features = NULL,
  cells = NULL,
  group.by = 'ident',
  group.bar = TRUE,
  group.colors = NULL,
  disp.min = -2.5,
  disp.max = NULL,
  slot = 'scale.data',
  assay = NULL,
  label = TRUE,
  size = 5.5,
  hjust = 0,
  angle = 45,
  raster = TRUE,
  draw.lines = TRUE,
  lines.width = NULL,
  group.bar.height = 0.02,
  combine = TRUE
) {
  cells <- cells %||% colnames(x = object)
  if (is.numeric(x = cells)) {
    cells <- colnames(x = object)[cells]
  }
  assay <- assay %||% DefaultAssay(object = object)
  DefaultAssay(object = object) <- assay
  features <- features %||% VariableFeatures(object = object)
  features <- rev(x = unique(x = features))
  disp.max <- disp.max %||% ifelse(
    test = slot == 'scale.data',
    yes = 2.5,
    no = 6
  )
  # make sure features are present
  possible.features <- rownames(x = GetAssayData(object = object, slot = slot))
  if (any(!features %in% possible.features)) {
    bad.features <- features[!features %in% possible.features]
    features <- features[features %in% possible.features]
    if(length(x = features) == 0) {
      stop("No requested features found in the ", slot, " slot for the ", assay, " assay.")
    }
    warning("The following features were omitted as they were not found in the ", slot,
            " slot for the ", assay, " assay: ", paste(bad.features, collapse = ", "))
  }
  data <- as.data.frame(x = as.matrix(x = t(x = GetAssayData(
    object = object,
    slot = slot)[features, cells, drop = FALSE])))

  object <- suppressMessages(expr = StashIdent(object = object, save.name = 'ident'))
  group.by <- group.by %||% 'ident'
  groups.use <- object[[group.by]][cells, , drop = FALSE]
  # group.use <- switch(
  #   EXPR = group.by,
  #   'ident' = Idents(object = object),
  #   object[[group.by, drop = TRUE]]
  # )
  # group.use <- factor(x = group.use[cells])
  plots <- vector(mode = 'list', length = ncol(x = groups.use))
  for (i in 1:ncol(x = groups.use)) {
    data.group <- data
    group.use <- groups.use[, i, drop = TRUE]
    if (!is.factor(x = group.use)) {
      group.use <- factor(x = group.use)
    }
    names(x = group.use) <- cells
    if (draw.lines) {
      # create fake cells to serve as the white lines, fill with NAs
      lines.width <- lines.width %||% ceiling(x = nrow(x = data.group) * 0.0025)
      placeholder.cells <- sapply(
        X = 1:(length(x = levels(x = group.use)) * lines.width),
        FUN = function(x) {
          return(RandomName(length = 20))
        }
      )
      placeholder.groups <- rep(x = levels(x = group.use), times = lines.width)
      group.levels <- levels(x = group.use)
      names(x = placeholder.groups) <- placeholder.cells
      group.use <- as.vector(x = group.use)
      names(x = group.use) <- cells
      group.use <- factor(x = c(group.use, placeholder.groups), levels = group.levels)
      na.data.group <- matrix(
        data = NA,
        nrow = length(x = placeholder.cells),
        ncol = ncol(x = data.group),
        dimnames = list(placeholder.cells, colnames(x = data.group))
      )
      data.group <- rbind(data.group, na.data.group)
    }
    lgroup <- length(levels(group.use))
    plot <- SingleRasterMap(
      data = data.group,
      raster = raster,
      disp.min = disp.min,
      disp.max = disp.max,
      feature.order = features,
      cell.order = names(x = sort(x = group.use)),
      group.by = group.use
    )
    if (group.bar) {
      # TODO: Change group.bar to annotation.bar
      default.colors <- c(hue_pal()(length(x = levels(x = group.use))))
      cols <- group.colors[1:length(x = levels(x = group.use))] %||% default.colors
      if (any(is.na(x = cols))) {
        cols[is.na(x = cols)] <- default.colors[is.na(x = cols)]
        cols <- Col2Hex(cols)
        col.dups <- sort(x = unique(x = which(x = duplicated(x = substr(
          x = cols,
          start = 1,
          stop = 7
        )))))
        through <- length(x = default.colors)
        while (length(x = col.dups) > 0) {
          pal.max <- length(x = col.dups) + through
          cols.extra <- hue_pal()(pal.max)[(through + 1):pal.max]
          cols[col.dups] <- cols.extra
          col.dups <- sort(x = unique(x = which(x = duplicated(x = substr(
            x = cols,
            start = 1,
            stop = 7
          )))))
        }
      }
      group.use2 <- sort(x = group.use)
      if (draw.lines) {
        na.group <- RandomName(length = 20)
        levels(x = group.use2) <- c(levels(x = group.use2), na.group)
        group.use2[placeholder.cells] <- na.group
        cols <- c(cols, "#FFFFFF")
      }
      pbuild <- ggplot_build(plot = plot)
      names(x = cols) <- levels(x = group.use2)
      # scale the height of the bar
      y.range <- diff(x = pbuild$layout$panel_params[[1]]$y.range)
      y.pos <- max(pbuild$layout$panel_params[[1]]$y.range) + y.range * 0.015
      y.max <- y.pos + group.bar.height * y.range
      plot <- plot +
        annotation_raster(
          raster = t(x = cols[group.use2]),
          xmin = -Inf,
          xmax = Inf,
          ymin = y.pos,
          ymax = y.max
        ) +
        coord_cartesian(ylim = c(0, y.max), clip = 'off') +
        scale_color_manual(values = cols)
      if (label) {
        x.max <- max(pbuild$layout$panel_params[[1]]$x.range)
        # Attempt to pull xdivs from x.major in ggplot2 < 3.3.0; if NULL, pull from the >= 3.3.0 slot
        x.divs <- pbuild$layout$panel_params[[1]]$x.major %||% pbuild$layout$panel_params[[1]]$x$break_positions()
        x <- data.frame(group = sort(x = group.use), x = x.divs)
        label.x.pos <- tapply(X = x$x, INDEX = x$group, FUN = median) * x.max
        label.x.pos <- data.frame(group = names(x = label.x.pos), label.x.pos)
        plot <- plot + geom_text(
          stat = "identity",
          data = label.x.pos,
          aes_string(label = 'group', x = 'label.x.pos'),
          y = y.max + y.max * 0.03 * 0.5,
          angle = angle,
          hjust = hjust,
          size = size
        )
        plot <- suppressMessages(plot + coord_cartesian(
          ylim = c(0, y.max + y.max * 0.002 * max(nchar(x = levels(x = group.use))) * size),
          clip = 'off')
        )
      }
    }
    plot <- plot + theme(line = element_blank())
    plots[[i]] <- plot
  }
  if (combine) {
    plots <- wrap_plots(plots)
  }
  return(plots)
}

#' Hashtag oligo heatmap
#'
#' Draws a heatmap of hashtag oligo signals across singlets/doublets/negative cells. Allows for the visualization of HTO demultiplexing results.
#'
#' @param object Seurat object. Assumes that the hash tag oligo (HTO) data has been added and normalized, and demultiplexing has been run with HTODemux().
#' @param classification The naming for metadata column with classification result from HTODemux().
#' @param global.classification The slot for metadata column specifying a cell as singlet/doublet/negative.
#' @param assay Hashtag assay name.
#' @param ncells Number of cells to plot. Default is to choose 5000 cells by random subsampling, to avoid having to draw exceptionally large heatmaps.
#' @param singlet.names Namings for the singlets. Default is to use the same names as HTOs.
#' @param raster If true, plot with geom_raster, else use geom_tile. geom_raster may look blurry on
#' some viewing applications such as Preview due to how the raster is interpolated. Set this to FALSE
#' if you are encountering that issue (note that plots may take longer to produce/render).
#' @return Returns a ggplot2 plot object.
#'
#' @importFrom ggplot2 guides
#' @export
#'
#' @seealso \code{\link{HTODemux}}
#'
#' @examples
#' \dontrun{
#' object <- HTODemux(object)
#' HTOHeatmap(object)
#' }
#'
HTOHeatmap <- function(
  object,
  assay = 'HTO',
  classification = paste0(assay, '_classification'),
  global.classification = paste0(assay, '_classification.global'),
  ncells = 5000,
  singlet.names = NULL,
  raster = TRUE
) {
  DefaultAssay(object = object) <- assay
  Idents(object = object) <- object[[classification, drop = TRUE]]
  if (ncells > ncol(x = object)) {
    warning("ncells (", ncells, ") is larger than the number of cells present in the provided object (", ncol(x = object), "). Plotting heatmap for all cells.")
  } else {
    object <- subset(
      x = object,
      cells = sample(x = colnames(x = object), size = ncells)
    )
  }
  classification <- object[[classification]]
  singlets <- which(x = object[[global.classification]] == 'Singlet')
  singlet.ids <- sort(x = unique(x = as.character(x = classification[singlets, ])))
  doublets <- which(object[[global.classification]] == 'Doublet')
  doublet.ids <- sort(x = unique(x = as.character(x = classification[doublets, ])))
  heatmap.levels <- c(singlet.ids, doublet.ids, 'Negative')
  object <- ScaleData(object = object, assay = assay, verbose = FALSE)
  data <- FetchData(object = object, vars = singlet.ids)
  Idents(object = object) <- factor(x = classification[, 1], levels = heatmap.levels)
  plot <- SingleRasterMap(
    data = data,
    raster = raster,
    feature.order = rev(x = singlet.ids),
    cell.order = names(x = sort(x = Idents(object = object))),
    group.by = Idents(object = object)
  ) + guides(color = FALSE)
  return(plot)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Expression by identity plots
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Single cell ridge plot
#'
#' Draws a ridge plot of single cell data (gene expression, metrics, PC
#' scores, etc.)
#'
#' @param object Seurat object
#' @param features Features to plot (gene expression, metrics, PC scores,
#' anything that can be retreived by FetchData)
#' @param cols Colors to use for plotting
#' @param idents Which classes to include in the plot (default is all)
#' @param sort Sort identity classes (on the x-axis) by the average
#' expression of the attribute being potted, can also pass 'increasing' or 'decreasing' to change sort direction
#' @param assay Name of assay to use, defaults to the active assay
#' @param group.by Group (color) cells in different ways (for example, orig.ident)
#' @param y.max Maximum y axis value
#' @param same.y.lims Set all the y-axis limits to the same values
#' @param log plot the feature axis on log scale
#' @param ncol Number of columns if multiple plots are displayed
#' @param slot Use non-normalized counts data for plotting
#' @param combine Combine plots into a single \code{\link[patchwork]{patchwork}ed}
#' ggplot object. If \code{FALSE}, return a list of ggplot objects
#'
#' @return A \code{\link[patchwork]{patchwork}ed} ggplot object if
#' \code{combine = TRUE}; otherwise, a list of ggplot objects
#'
#' @export
#'
#' @examples
#' RidgePlot(object = pbmc_small, features = 'PC_1')
#'
RidgePlot <- function(
  object,
  features,
  cols = NULL,
  idents = NULL,
  sort = FALSE,
  assay = NULL,
  group.by = NULL,
  y.max = NULL,
  same.y.lims = FALSE,
  log = FALSE,
  ncol = NULL,
  slot = 'data',
  combine = TRUE
) {
  return(ExIPlot(
    object = object,
    type = 'ridge',
    features = features,
    idents = idents,
    ncol = ncol,
    sort = sort,
    assay = assay,
    y.max = y.max,
    same.y.lims = same.y.lims,
    cols = cols,
    group.by = group.by,
    log = log,
    slot = slot,
    combine = combine
  ))
}

#' Single cell violin plot
#'
#' Draws a violin plot of single cell data (gene expression, metrics, PC
#' scores, etc.)
#'
#' @inheritParams RidgePlot
#' @param pt.size Point size for geom_violin
#' @param split.by A variable to split the violin plots by,
#' @param multi.group  plot each group of the split violin plots by multiple or single violin shapes
#' see \code{\link{FetchData}} for more details
#' @param adjust Adjust parameter for geom_violin
#'
#' @return A \code{\link[patchwork]{patchwork}ed} ggplot object if
#' \code{combine = TRUE}; otherwise, a list of ggplot objects
#'
#' @export
#'
#' @seealso \code{\link{FetchData}}
#'
#' @examples
#' VlnPlot(object = pbmc_small, features = 'PC_1')
#' VlnPlot(object = pbmc_small, features = 'LYZ', split.by = 'groups')
#'
VlnPlot <- function(
  object,
  features,
  cols = NULL,
  pt.size = 1,
  idents = NULL,
  sort = FALSE,
  assay = NULL,
  group.by = NULL,
  split.by = NULL,
  adjust = 1,
  y.max = NULL,
  same.y.lims = FALSE,
  log = FALSE,
  ncol = NULL,
  slot = 'data',
  multi.group = FALSE,
  combine = TRUE
) {
  return(ExIPlot(
    object = object,
    type = ifelse(test = multi.group, yes = 'multiViolin', no = 'violin'),
    features = features,
    idents = idents,
    ncol = ncol,
    sort = sort,
    assay = assay,
    y.max = y.max,
    same.y.lims = same.y.lims,
    adjust = adjust,
    pt.size = pt.size,
    cols = cols,
    group.by = group.by,
    split.by = split.by,
    log = log,
    slot = slot,
    combine = combine
  ))
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Dimensional reduction plots
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Color dimensional reduction plot by tree split
#'
#' Returns a DimPlot colored based on whether the cells fall in clusters
#' to the left or to the right of a node split in the cluster tree.
#'
#' @param object Seurat object
#' @param node Node in cluster tree on which to base the split
#' @param left.color Color for the left side of the split
#' @param right.color Color for the right side of the split
#' @param other.color Color for all other cells
#' @inheritDotParams DimPlot -object
#'
#' @return Returns a DimPlot
#'
#' @export
#'
#' @seealso \code{\link{DimPlot}}
#'
#' @examples
#' pbmc_small
#' pbmc_small <- BuildClusterTree(object = pbmc_small, verbose = FALSE)
#' PlotClusterTree(pbmc_small)
#' ColorDimSplit(pbmc_small, node = 5)
#'
ColorDimSplit <- function(
  object,
  node,
  left.color = 'red',
  right.color = 'blue',
  other.color = 'grey50',
  ...
) {
  CheckDots(..., fxns = 'DimPlot')
  tree <- Tool(object = object, slot = "BuildClusterTree")
  split <- tree$edge[which(x = tree$edge[, 1] == node), ][, 2]
  all.children <- sort(x = tree$edge[, 2][! tree$edge[, 2] %in% tree$edge[, 1]])
  left.group <- DFT(tree = tree, node = split[1], only.children = TRUE)
  right.group <- DFT(tree = tree, node = split[2], only.children = TRUE)
  if (any(is.na(x = left.group))) {
    left.group <- split[1]
  }
  if (any(is.na(x = right.group))) {
    right.group <- split[2]
  }
  left.group <- MapVals(v = left.group, from = all.children, to = tree$tip.label)
  right.group <- MapVals(v = right.group, from = all.children, to = tree$tip.label)
  remaining.group <- setdiff(x = tree$tip.label, y = c(left.group, right.group))
  left.cells <- WhichCells(object = object, ident = left.group)
  right.cells <- WhichCells(object = object, ident = right.group)
  remaining.cells <- WhichCells(object = object, ident = remaining.group)
  object <- SetIdent(
    object = object,
    cells = left.cells,
    value = "Left Split"
  )
  object <- SetIdent(
    object = object,
    cells = right.cells,
    value = "Right Split"
  )
  object <- SetIdent(
    object = object,
    cells = remaining.cells,
    value = "Not in Split"
  )
  levels(x = object) <- c("Left Split", "Right Split", "Not in Split")
  colors.use = c(left.color, right.color, other.color)
  return(DimPlot(object = object, cols = colors.use, ...))
}

#' Dimensional reduction plot
#'
#' Graphs the output of a dimensional reduction technique on a 2D scatter plot where each point is a
#' cell and it's positioned based on the cell embeddings determined by the reduction technique. By
#' default, cells are colored by their identity class (can be changed with the group.by parameter).
#'
#' @param object Seurat object
#' @param dims Dimensions to plot, must be a two-length numeric vector specifying x- and y-dimensions
#' @param cells Vector of cells to plot (default is all cells)
#' @param cols Vector of colors, each color corresponds to an identity class. This may also be a single character
#' or numeric value corresponding to a palette as specified by \code{\link[RColorBrewer]{brewer.pal.info}}.
#' By default, ggplot2 assigns colors. We also include a number of palettes from the pals package.
#' See \code{\link{DiscretePalette}} for details.
#' @param pt.size Adjust point size for plotting
#' @param reduction Which dimensionality reduction to use. If not specified, first searches for umap, then tsne, then pca
#' @param group.by Name of one or more metadata columns to group (color) cells by
#' (for example, orig.ident); pass 'ident' to group by identity class
#' @param split.by Name of a metadata column to split plot by;
#' see \code{\link{FetchData}} for more details
#' @param shape.by If NULL, all points are circles (default). You can specify any
#' cell attribute (that can be pulled with FetchData) allowing for both
#' different colors and different shapes on cells
#' @param order Specify the order of plotting for the idents. This can be
#' useful for crowded plots if points of interest are being buried. Provide
#' either a full list of valid idents or a subset to be plotted last (on top)
#' @param label Whether to label the clusters
#' @param label.size Sets size of labels
#' @param repel Repel labels
#' @param cells.highlight A list of character or numeric vectors of cells to
#' highlight. If only one group of cells desired, can simply
#' pass a vector instead of a list. If set, colors selected cells to the color(s)
#' in \code{cols.highlight} and other cells black (white if dark.theme = TRUE);
#' will also resize to the size(s) passed to \code{sizes.highlight}
#' @param cols.highlight A vector of colors to highlight the cells as; will
#' repeat to the length groups in cells.highlight
#' @param sizes.highlight Size of highlighted cells; will repeat to the length
#' groups in cells.highlight
#' @param na.value Color value for NA points when using custom scale
#' @param ncol Number of columns for display when combining plots
#' @param combine Combine plots into a single \code{\link[patchwork]{patchwork}ed}
#' ggplot object. If \code{FALSE}, return a list of ggplot objects
#'
#' @return A \code{\link[patchwork]{patchwork}ed} ggplot object if
#' \code{combine = TRUE}; otherwise, a list of ggplot objects
#'
#' @importFrom rlang !!
#' @importFrom ggplot2 facet_wrap vars sym
#' @importFrom patchwork wrap_plots
#'
#' @export
#'
#' @note For the old \code{do.hover} and \code{do.identify} functionality, please see
#' \code{HoverLocator} and \code{CellSelector}, respectively.
#'
#' @aliases TSNEPlot PCAPlot ICAPlot
#' @seealso \code{\link{FeaturePlot}} \code{\link{HoverLocator}}
#' \code{\link{CellSelector}} \code{\link{FetchData}}
#'
#' @examples
#' DimPlot(object = pbmc_small)
#' DimPlot(object = pbmc_small, split.by = 'ident')
#'
DimPlot <- function(
  object,
  dims = c(1, 2),
  cells = NULL,
  cols = NULL,
  pt.size = NULL,
  reduction = NULL,
  group.by = NULL,
  split.by = NULL,
  shape.by = NULL,
  order = NULL,
  label = FALSE,
  label.size = 4,
  repel = FALSE,
  cells.highlight = NULL,
  cols.highlight = '#DE2D26',
  sizes.highlight = 1,
  na.value = 'grey50',
  ncol = NULL,
  combine = TRUE
) {
  if (length(x = dims) != 2) {
    stop("'dims' must be a two-length vector")
  }
  reduction <- reduction %||% DefaultDimReduc(object = object)
  cells <- cells %||% colnames(x = object)
  data <- Embeddings(object = object[[reduction]])[cells, dims]
  data <- as.data.frame(x = data)
  dims <- paste0(Key(object = object[[reduction]]), dims)
  object[['ident']] <- Idents(object = object)
  orig.groups <- group.by
  group.by <- group.by %||% 'ident'
  data[, group.by] <- object[[group.by]][cells, , drop = FALSE]
  for (group in group.by) {
    if (!is.factor(x = data[, group])) {
      data[, group] <- factor(x = data[, group])
    }
  }
  if (!is.null(x = shape.by)) {
    data[, shape.by] <- object[[shape.by, drop = TRUE]]
  }
  if (!is.null(x = split.by)) {
    data[, split.by] <- object[[split.by, drop = TRUE]]
  }
  plots <- lapply(
    X = group.by,
    FUN = function(x) {
      plot <- SingleDimPlot(
        data = data[, c(dims, x, split.by, shape.by)],
        dims = dims,
        col.by = x,
        cols = cols,
        pt.size = pt.size,
        shape.by = shape.by,
        order = order,
        label = FALSE,
        cells.highlight = cells.highlight,
        cols.highlight = cols.highlight,
        sizes.highlight = sizes.highlight,
        na.value = na.value
      )
      if (label) {
        plot <- LabelClusters(
          plot = plot,
          id = x,
          repel = repel,
          size = label.size,
          split.by = split.by
        )
      }
      if (!is.null(x = split.by)) {
        plot <- plot + FacetTheme() +
          facet_wrap(
            facets = vars(!!sym(x = split.by)),
            ncol = if (length(x = group.by) > 1 || is.null(x = ncol)) {
              length(x = unique(x = data[, split.by]))
            } else {
              ncol
            }
          )
      }
      return(plot)
    }
  )
  if (!is.null(x = split.by)) {
    ncol <- 1
  }
  if (combine) {
    plots <- wrap_plots(plots, ncol = orig.groups %iff% ncol)
  }
  return(plots)
}

#' Visualize 'features' on a dimensional reduction plot
#'
#' Colors single cells on a dimensional reduction plot according to a 'feature'
#' (i.e. gene expression, PC scores, number of genes detected, etc.)
#'
#' @inheritParams DimPlot
#' @param order Boolean determining whether to plot cells in order of expression. Can be useful if
#' cells expressing given feature are getting buried.
#' @param features Vector of features to plot. Features can come from:
#' \itemize{
#'     \item An \code{Assay} feature (e.g. a gene name - "MS4A1")
#'     \item A column name from meta.data (e.g. mitochondrial percentage - "percent.mito")
#'     \item A column name from a \code{DimReduc} object corresponding to the cell embedding values
#'     (e.g. the PC 1 scores - "PC_1")
#' }
#' @param cols The two colors to form the gradient over. Provide as string vector with
#' the first color corresponding to low values, the second to high. Also accepts a Brewer
#' color scale or vector of colors. Note: this will bin the data into number of colors provided.
#' When blend is \code{TRUE}, takes anywhere from 1-3 colors:
#' \describe{
#'   \item{1 color:}{Treated as color for double-negatives, will use default colors 2 and 3 for per-feature expression}
#'   \item{2 colors:}{Treated as colors for per-feature expression, will use default color 1 for double-negatives}
#'   \item{3+ colors:}{First color used for double-negatives, colors 2 and 3 used for per-feature expression, all others ignored}
#' }
#' @param min.cutoff,max.cutoff Vector of minimum and maximum cutoff values for each feature,
#'  may specify quantile in the form of 'q##' where '##' is the quantile (eg, 'q1', 'q10')
#' @param split.by A factor in object metadata to split the feature plot by, pass 'ident'
#'  to split by cell identity'; similar to the old \code{FeatureHeatmap}
#' @param slot Which slot to pull expression data from?
#' @param blend Scale and blend expression values to visualize coexpression of two features
#' @param blend.threshold The color cutoff from weak signal to strong signal; ranges from 0 to 1.
#' @param ncol Number of columns to combine multiple feature plots to, ignored if \code{split.by} is not \code{NULL}
#' @param coord.fixed Plot cartesian coordinates with fixed aspect ratio
#' @param by.col If splitting by a factor, plot the splits per column with the features as rows; ignored if \code{blend = TRUE}
#' @param sort.cell If \code{TRUE}, the positive cells will overlap the negative cells
#' @param combine Combine plots into a single \code{\link[patchwork]{patchwork}ed}
#' ggplot object. If \code{FALSE}, return a list of ggplot objects
#'
#' @return A \code{\link[patchwork]{patchwork}ed} ggplot object if
#' \code{combine = TRUE}; otherwise, a list of ggplot objects
#'
#' @importFrom grDevices rgb
#' @importFrom patchwork wrap_plots
#' @importFrom cowplot theme_cowplot
#' @importFrom RColorBrewer brewer.pal.info
#' @importFrom ggplot2 labs scale_x_continuous scale_y_continuous theme element_rect
#' dup_axis guides element_blank element_text margin scale_color_brewer scale_color_gradientn
#' scale_color_manual coord_fixed ggtitle
#'
#' @export
#'
#' @note For the old \code{do.hover} and \code{do.identify} functionality, please see
#' \code{HoverLocator} and \code{CellSelector}, respectively.
#'
#' @aliases FeatureHeatmap
#' @seealso \code{\link{DimPlot}} \code{\link{HoverLocator}}
#' \code{\link{CellSelector}}
#'
#' @examples
#' FeaturePlot(object = pbmc_small, features = 'PC_1')
#'
FeaturePlot <- function(
  object,
  features,
  dims = c(1, 2),
  cells = NULL,
  cols = if (blend) {
    c('lightgrey', '#ff0000', '#00ff00')
  } else {
    c('lightgrey', 'blue')
  },
  pt.size = NULL,
  order = FALSE,
  min.cutoff = NA,
  max.cutoff = NA,
  reduction = NULL,
  split.by = NULL,
  shape.by = NULL,
  slot = 'data',
  blend = FALSE,
  blend.threshold = 0.5,
  label = FALSE,
  label.size = 4,
  repel = FALSE,
  ncol = NULL,
  coord.fixed = FALSE,
  by.col = TRUE,
  sort.cell = FALSE,
  combine = TRUE
) {
  # Set a theme to remove right-hand Y axis lines
  # Also sets right-hand Y axis text label formatting
  no.right <- theme(
    axis.line.y.right = element_blank(),
    axis.ticks.y.right = element_blank(),
    axis.text.y.right = element_blank(),
    axis.title.y.right = element_text(
      face = "bold",
      size = 14,
      margin = margin(r = 7)
    )
  )
  # Get the DimReduc to use
  reduction <- reduction %||% DefaultDimReduc(object = object)
  if (length(x = dims) != 2 || !is.numeric(x = dims)) {
    stop("'dims' must be a two-length integer vector")
  }
  # Figure out blending stuff
  if (blend && length(x = features) != 2) {
    stop("Blending feature plots only works with two features")
  }
  # Set color scheme for blended FeaturePlots
  if (blend) {
    default.colors <- eval(expr = formals(fun = FeaturePlot)$cols)
    cols <- switch(
      EXPR = as.character(x = length(x = cols)),
      '0' = {
        warning(
          "No colors provided, using default colors",
          call. = FALSE,
          immediate. = TRUE
        )
        default.colors
      },
      '1' = {
        warning(
          "Only one color provided, assuming specified is double-negative and augmenting with default colors",
          call. = FALSE,
          immediate. = TRUE
        )
        c(cols, default.colors[2:3])
      },
      '2' = {
        warning(
          "Only two colors provided, assuming specified are for features and agumenting with '",
          default.colors[1],
          "' for double-negatives",
          call. = FALSE,
          immediate. = TRUE
        )
        c(default.colors[1], cols)
      },
      '3' = cols,
      {
        warning(
          "More than three colors provided, using only first three",
          call. = FALSE,
          immediate. = TRUE
        )
        cols[1:3]
      }
    )
  }
  if (blend && length(x = cols) != 3) {
    stop("Blending feature plots only works with three colors; first one for negative cells")
  }
  # Name the reductions
  dims <- paste0(Key(object = object[[reduction]]), dims)
  cells <- cells %||% colnames(x = object)
  # Get plotting data
  data <- FetchData(
    object = object,
    vars = c(dims, 'ident', features),
    cells = cells,
    slot = slot
  )
  # Check presence of features/dimensions
  if (ncol(x = data) < 4) {
    stop(
      "None of the requested features were found: ",
      paste(features, collapse = ', '),
      " in slot ",
      slot,
      call. = FALSE
    )
  } else if (!all(dims %in% colnames(x = data))) {
    stop("The dimensions requested were not found", call. = FALSE)
  }
  features <- colnames(x = data)[4:ncol(x = data)]
  # Determine cutoffs
  min.cutoff <- mapply(
    FUN = function(cutoff, feature) {
      return(ifelse(
        test = is.na(x = cutoff),
        yes = min(data[, feature]),
        no = cutoff
      ))
    },
    cutoff = min.cutoff,
    feature = features
  )
  max.cutoff <- mapply(
    FUN = function(cutoff, feature) {
      return(ifelse(
        test = is.na(x = cutoff),
        yes = max(data[, feature]),
        no = cutoff
      ))
    },
    cutoff = max.cutoff,
    feature = features
  )
  check.lengths <- unique(x = vapply(
    X = list(features, min.cutoff, max.cutoff),
    FUN = length,
    FUN.VALUE = numeric(length = 1)
  ))
  if (length(x = check.lengths) != 1) {
    stop("There must be the same number of minimum and maximum cuttoffs as there are features")
  }
  brewer.gran <- ifelse(
    test = length(x = cols) == 1,
    yes = brewer.pal.info[cols, ]$maxcolors,
    no = length(x = cols)
  )
  # Apply cutoffs
  data[, 4:ncol(x = data)] <- sapply(
    X = 4:ncol(x = data),
    FUN = function(index) {
      data.feature <- as.vector(x = data[, index])
      min.use <- SetQuantile(cutoff = min.cutoff[index - 3], data.feature)
      max.use <- SetQuantile(cutoff = max.cutoff[index - 3], data.feature)
      data.feature[data.feature < min.use] <- min.use
      data.feature[data.feature > max.use] <- max.use
      if (brewer.gran == 2) {
        return(data.feature)
      }
      data.cut <- if (all(data.feature == 0)) {
        0
      }
      else {
        as.numeric(x = as.factor(x = cut(
          x = as.numeric(x = data.feature),
          breaks = brewer.gran
        )))
      }
      return(data.cut)
    }
  )
  colnames(x = data)[4:ncol(x = data)] <- features
  rownames(x = data) <- cells
  # Figure out splits (FeatureHeatmap)
  data$split <- if (is.null(x = split.by)) {
    RandomName()
  } else {
    switch(
      EXPR = split.by,
      ident = Idents(object = object)[cells],
      object[[split.by, drop = TRUE]][cells]
    )
  }
  if (!is.factor(x = data$split)) {
    data$split <- factor(x = data$split)
  }
  # Set shaping variable
  if (!is.null(x = shape.by)) {
    data[, shape.by] <- object[[shape.by, drop = TRUE]]
  }
  # Make list of plots
  plots <- vector(
    mode = "list",
    length = ifelse(
      test = blend,
      yes = 4,
      no = length(x = features) * length(x = levels(x = data$split))
    )
  )
  # Apply common limits
  xlims <- c(floor(x = min(data[, dims[1]])), ceiling(x = max(data[, dims[1]])))
  ylims <- c(floor(min(data[, dims[2]])), ceiling(x = max(data[, dims[2]])))
  # Set blended colors
  if (blend) {
    ncol <- 4
    color.matrix <- BlendMatrix(
      two.colors = cols[2:3],
      col.threshold = blend.threshold,
      negative.color = cols[1]
    )
    cols <- cols[2:3]
    colors <- list(
      color.matrix[, 1],
      color.matrix[1, ],
      as.vector(x = color.matrix)
    )
  }
  # Make the plots
  for (i in 1:length(x = levels(x = data$split))) {
    # Figre out which split we're working with
    ident <- levels(x = data$split)[i]
    data.plot <- data[as.character(x = data$split) == ident, , drop = FALSE]
    # Blend expression values
    if (blend) {
      features <- features[1:2]
      no.expression <- features[colMeans(x = data.plot[, features]) == 0]
      if (length(x = no.expression) != 0) {
        stop(
          "The following features have no value: ",
          paste(no.expression, collapse = ', '),
          call. = FALSE
        )
      }
      data.plot <- cbind(data.plot[, c(dims, 'ident')], BlendExpression(data = data.plot[, features[1:2]]))
      features <- colnames(x = data.plot)[4:ncol(x = data.plot)]
    }
    # Make per-feature plots
    for (j in 1:length(x = features)) {
      feature <- features[j]
      # Get blended colors
      if (blend) {
        cols.use <- as.numeric(x = as.character(x = data.plot[, feature])) + 1
        cols.use <- colors[[j]][sort(x = unique(x = cols.use))]
      } else {
        cols.use <- NULL
      }
      data.single <- data.plot[, c(dims, 'ident', feature, shape.by)]
      if (sort.cell) {
        data.single <- data.single[order(data.single[, feature]),]
      }
      # Make the plot
      plot <- SingleDimPlot(
        data = data.single,
        dims = dims,
        col.by = feature,
        order = order,
        pt.size = pt.size,
        cols = cols.use,
        shape.by = shape.by,
        label = FALSE
      ) +
        scale_x_continuous(limits = xlims) +
        scale_y_continuous(limits = ylims) +
        theme_cowplot() +
        theme(plot.title = element_text(hjust = 0.5))
      # Add labels
      if (label) {
        plot <- LabelClusters(
          plot = plot,
          id = 'ident',
          repel = repel,
          size = label.size
        )
      }
      # Make FeatureHeatmaps look nice(ish)
      if (length(x = levels(x = data$split)) > 1) {
        plot <- plot + theme(panel.border = element_rect(fill = NA, colour = 'black'))
        # Add title
        plot <- plot + if (i == 1) {
          labs(title = feature)
        } else {
          labs(title = NULL)
        }
        # Add second axis
        if (j == length(x = features) && !blend) {
          suppressMessages(
            expr = plot <- plot +
              scale_y_continuous(
                sec.axis = dup_axis(name = ident),
                limits = ylims
              ) +
              no.right
          )
        }
        # Remove left Y axis
        if (j != 1) {
          plot <- plot + theme(
            axis.line.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.text.y = element_blank(),
            axis.title.y.left = element_blank()
          )
        }
        # Remove bottom X axis
        if (i != length(x = levels(x = data$split))) {
          plot <- plot + theme(
            axis.line.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.text.x = element_blank(),
            axis.title.x = element_blank()
          )
        }
      } else {
        plot <- plot + labs(title = feature)
      }
      # Add colors scale for normal FeaturePlots
      if (!blend) {
        plot <- plot + guides(color = NULL)
        cols.grad <- cols
        if (length(x = cols) == 1) {
          plot <- plot + scale_color_brewer(palette = cols)
        } else if (length(x = cols) > 1) {
          unique.feature.exp <- unique(data.plot[, feature])
          if (length(unique.feature.exp) == 1) {
            warning("All cells have the same value (", unique.feature.exp, ") of ", feature, ".")
            if (unique.feature.exp == 0) {
              cols.grad <- cols[1]
            } else{
            cols.grad <- cols
            }
          }
          plot <- suppressMessages(
            expr = plot + scale_color_gradientn(
              colors = cols.grad,
              guide = "colorbar"
            )
          )
        }
      }
      # Add coord_fixed
      if (coord.fixed) {
        plot <- plot + coord_fixed()
      }
      # I'm not sure why, but sometimes the damn thing fails without this
      # Thanks ggplot2
      plot <- plot
      # Place the plot
      plots[[(length(x = features) * (i - 1)) + j]] <- plot
    }
  }
  # Add blended color key
  if (blend) {
    blend.legend <- BlendMap(color.matrix = color.matrix)
    for (ii in 1:length(x = levels(x = data$split))) {
      suppressMessages(expr = plots <- append(
        x = plots,
        values = list(
          blend.legend +
            scale_y_continuous(
              sec.axis = dup_axis(name = ifelse(
                test = length(x = levels(x = data$split)) > 1,
                yes = levels(x = data$split)[ii],
                no = ''
              )),
              expand = c(0, 0)
            ) +
            labs(
              x = features[1],
              y = features[2],
              title = if (ii == 1) {
                paste('Color threshold:', blend.threshold)
              } else {
                NULL
              }
            ) +
            no.right
        ),
        after = 4 * ii - 1
      ))
    }
  }
  # Remove NULL plots
  plots <- Filter(f = Negate(f = is.null), x = plots)
  # Combine the plots
  if (is.null(x = ncol)) {
    ncol <- 2
    if (length(x = features) == 1) {
      ncol <- 1
    }
    if (length(x = features) > 6) {
      ncol <- 3
    }
    if (length(x = features) > 9) {
      ncol <- 4
    }
  }
  ncol <- ifelse(
    test = is.null(x = split.by) || blend,
    yes = ncol,
    no = length(x = features)
  )
  legend <- if (blend) {
    'none'
  } else {
    split.by %iff% 'none'
  }
  # Transpose the FeatureHeatmap matrix (not applicable for blended FeaturePlots)
  if (combine) {
    if (by.col && !is.null(x = split.by) && !blend) {
      plots <- lapply(
        X = plots,
        FUN = function(x) {
          return(suppressMessages(
            expr = x +
              theme_cowplot() +
              ggtitle("") +
              scale_y_continuous(sec.axis = dup_axis(name = ""), limits = ylims) +
              no.right
          ))
        }
      )
      nsplits <- length(x = levels(x = data$split))
      idx <- 1
      for (i in (length(x = features) * (nsplits - 1) + 1):(length(x = features) * nsplits)) {
        plots[[i]] <- suppressMessages(
          expr = plots[[i]] +
            scale_y_continuous(
              sec.axis = dup_axis(name = features[[idx]]),
              limits = ylims
            ) +
            no.right
        )
        idx <- idx + 1
      }
      idx <- 1
      for (i in which(x = 1:length(x = plots) %% length(x = features) == 1)) {
        plots[[i]] <- plots[[i]] +
          ggtitle(levels(x = data$split)[[idx]]) +
          theme(plot.title = element_text(hjust = 0.5))
        idx <- idx + 1
      }
      idx <- 1
      if (length(x = features) == 1) {
        for (i in 1:length(x = plots)) {
          plots[[i]] <- plots[[i]] +
            ggtitle(levels(x = data$split)[[idx]]) +
            theme(plot.title = element_text(hjust = 0.5))
          idx <- idx + 1
        }
        ncol <- nsplits
        nrow <- 1
      } else {
        nrow <- split.by %iff% length(x = levels(x = data$split))
      }
      plots <- plots[c(do.call(
        what = rbind,
        args = split(x = 1:length(x = plots), f = ceiling(x = seq_along(along.with = 1:length(x = plots)) / length(x = features)))
      ))]
      plots <- wrap_plots(plots, ncol = ncol, nrow = nrow)
      if (!is.null(x = legend) && legend == 'none') {
        plots <- plots & NoLegend()
      }
    } else {
      plots <- wrap_plots(plots, ncol = ncol, nrow = split.by %iff% length(x = levels(x = data$split)))
    }
    if (!is.null(x = legend) && legend == 'none') {
      plots <- plots & NoLegend()
    }
  }
  return(plots)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Scatter plots
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Cell-cell scatter plot
#'
#' Creates a plot of scatter plot of features across two single cells. Pearson
#' correlation between the two cells is displayed above the plot.
#'
#' @inheritParams FeatureScatter
#' @param cell1 Cell 1 name
#' @param cell2 Cell 2 name
#' @param features Features to plot (default, all features)
#' @param highlight Features to highlight
#' @return A ggplot object
#'
#' @export
#'
#' @aliases CellPlot
#'
#' @examples
#' CellScatter(object = pbmc_small, cell1 = 'ATAGGAGAAACAGA', cell2 = 'CATCAGGATGCACA')
#'
CellScatter <- function(
  object,
  cell1,
  cell2,
  features = NULL,
  highlight = NULL,
  cols = NULL,
  pt.size = 1,
  smooth = FALSE
) {
  features <- features %||% rownames(x = object)
  data <- FetchData(
    object = object,
    vars = features,
    cells = c(cell1, cell2)
  )
  data <- as.data.frame(x = t(x = data))
  plot <- SingleCorPlot(
    data = data,
    cols = cols,
    pt.size = pt.size,
    rows.highlight = highlight,
    smooth = smooth
  )
  return(plot)
}

#' Scatter plot of single cell data
#'
#' Creates a scatter plot of two features (typically feature expression), across a
#' set of single cells. Cells are colored by their identity class. Pearson
#' correlation between the two features is displayed above the plot.
#'
#' @param object Seurat object
#' @param feature1 First feature to plot. Typically feature expression but can also
#' be metrics, PC scores, etc. - anything that can be retreived with FetchData
#' @param feature2 Second feature to plot.
#' @param cells Cells to include on the scatter plot.
#' @param group.by Name of one or more metadata columns to group (color) cells by
#' (for example, orig.ident); pass 'ident' to group by identity class
#' @param cols Colors to use for identity class plotting.
#' @param pt.size Size of the points on the plot
#' @param shape.by Ignored for now
#' @param span Spline span in loess function call, if \code{NULL}, no spline added
#' @param smooth Smooth the graph (similar to smoothScatter)
#' @param slot Slot to pull data from, should be one of 'counts', 'data', or 'scale.data'
#'
#' @return A ggplot object
#'
#' @importFrom ggplot2 geom_smooth aes_string
#' @export
#'
#' @aliases GenePlot
#'
#' @examples
#' FeatureScatter(object = pbmc_small, feature1 = 'CD9', feature2 = 'CD3E')
#'
FeatureScatter <- function(
  object,
  feature1,
  feature2,
  cells = NULL,
  group.by = NULL,
  cols = NULL,
  pt.size = 1,
  shape.by = NULL,
  span = NULL,
  smooth = FALSE,
  slot = 'data'
) {
  cells <- cells %||% colnames(x = object)
  group.by <- group.by %||% Idents(object = object)[cells]
  if (length(x = group.by) == 1) {
    group.by <- object[[]][, group.by]
  }
  plot <- SingleCorPlot(
    data = FetchData(
      object = object,
      vars = c(feature1, feature2),
      cells = cells,
      slot = slot
    ),
    col.by = group.by,
    cols = cols,
    pt.size = pt.size,
    smooth = smooth,
    legend.title = 'Identity',
    span = span
  )
  return(plot)
}

#' View variable features
#'
#' @inheritParams FeatureScatter
#' @inheritParams HVFInfo
#' @param cols Colors to specify non-variable/variable status
#' @param assay Assay to pull variable features from
#' @param log Plot the x-axis in log scale
#'
#' @return A ggplot object
#'
#' @importFrom ggplot2 labs scale_color_manual scale_x_log10
#' @export
#'
#' @aliases VariableGenePlot MeanVarPlot
#'
#' @seealso \code{\link{FindVariableFeatures}}
#'
#' @examples
#' VariableFeaturePlot(object = pbmc_small)
#'
VariableFeaturePlot <- function(
  object,
  cols = c('black', 'red'),
  pt.size = 1,
  log = NULL,
  selection.method = NULL,
  assay = NULL
) {
  if (length(x = cols) != 2) {
    stop("'cols' must be of length 2")
  }
  hvf.info <- HVFInfo(
    object = object,
    assay = assay,
    selection.method = selection.method,
    status = TRUE
    )
  var.status <- c('no', 'yes')[unlist(x = hvf.info[, ncol(x = hvf.info)]) + 1]
  hvf.info <- hvf.info[, c(1, 3)]
  axis.labels <- switch(
    EXPR = colnames(x = hvf.info)[2],
    'variance.standardized' = c('Average Expression', 'Standardized Variance'),
    'dispersion.scaled' = c('Average Expression', 'Dispersion'),
    'residual_variance' = c('Geometric Mean of Expression', 'Residual Variance')
  )
  log <- log %||% (any(c('variance.standardized', 'residual_variance') %in% colnames(x = hvf.info)))
  # var.features <- VariableFeatures(object = object, assay = assay)
  # var.status <- ifelse(
  #   test = rownames(x = hvf.info) %in% var.features,
  #   yes = 'yes',
  #   no = 'no'
  # )
  plot <- SingleCorPlot(
    data = hvf.info,
    col.by = var.status,
    pt.size = pt.size
  )
  plot <- plot +
    labs(title = NULL, x = axis.labels[1], y = axis.labels[2]) +
    scale_color_manual(
      labels = paste(c('Non-variable', 'Variable'), 'count:', table(var.status)),
      values = cols
    )
  if (log) {
    plot <- plot + scale_x_log10()
  }
  return(plot)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Polygon Plots
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Polygon DimPlot
#'
#' Plot cells as polygons, rather than single points. Color cells by identity, or a categorical variable
#' in metadata
#'
#' @inheritParams PolyFeaturePlot
#' @param group.by A grouping variable present in the metadata. Default is to use the groupings present
#' in the current cell identities (\code{Idents(object = object)})
#'
#' @return Returns a ggplot object
#'
#' @export
#'
PolyDimPlot <- function(
  object,
  group.by = NULL,
  cells = NULL,
  poly.data = 'spatial',
  flip.coords = FALSE
) {
  polygons <- Misc(object = object, slot = poly.data)
  if (is.null(x = polygons)) {
    stop("Could not find polygon data in misc slot")
  }
  group.by <- group.by %||% 'ident'
  group.data <- FetchData(
    object = object,
    vars = group.by,
    cells = cells
  )
  group.data$cell <- rownames(x = group.data)
  data <- merge(x = polygons, y = group.data, by = 'cell')
  if (flip.coords) {
    coord.x <- data$x
    data$x <- data$y
    data$y <- coord.x
  }
  plot <- SinglePolyPlot(data = data, group.by = group.by)
  return(plot)
}

#' Polygon FeaturePlot
#'
#' Plot cells as polygons, rather than single points. Color cells by any value accessible by \code{\link{FetchData}}.
#'
#' @inheritParams FeaturePlot
#' @param poly.data Name of the polygon dataframe in the misc slot
#' @param ncol Number of columns to split the plot into
#' @param common.scale ...
#' @param flip.coords Flip x and y coordinates
#'
#' @return Returns a ggplot object
#'
#' @importFrom ggplot2 scale_fill_viridis_c facet_wrap
#'
#' @export
#'
PolyFeaturePlot <- function(
  object,
  features,
  cells = NULL,
  poly.data = 'spatial',
  ncol = ceiling(x = length(x = features) / 2),
  min.cutoff = 0,
  max.cutoff = NA,
  common.scale = TRUE,
  flip.coords = FALSE
) {
  polygons <- Misc(object = object, slot = poly.data)
  if (is.null(x = polygons)) {
    stop("Could not find polygon data in misc slot")
  }
  assay.data <- FetchData(
    object = object,
    vars = features,
    cells = cells
  )
  features <- colnames(x = assay.data)
  cells <- rownames(x = assay.data)
  min.cutoff <- mapply(
    FUN = function(cutoff, feature) {
      return(ifelse(
        test = is.na(x = cutoff),
        yes = min(assay.data[, feature]),
        no = cutoff
      ))
    },
    cutoff = min.cutoff,
    feature = features
  )
  max.cutoff <- mapply(
    FUN = function(cutoff, feature) {
      return(ifelse(
        test = is.na(x = cutoff),
        yes = max(assay.data[, feature]),
        no = cutoff
      ))
    },
    cutoff = max.cutoff,
    feature = features
  )
  check.lengths <- unique(x = vapply(
    X = list(features, min.cutoff, max.cutoff),
    FUN = length,
    FUN.VALUE = numeric(length = 1)
  ))
  if (length(x = check.lengths) != 1) {
    stop("There must be the same number of minimum and maximum cuttoffs as there are features")
  }
  assay.data <- mapply(
    FUN = function(feature, min, max) {
      return(ScaleColumn(vec = assay.data[, feature], cutoffs = c(min, max)))
    },
    feature = features,
    min = min.cutoff,
    max = max.cutoff
  )
  if (common.scale) {
    assay.data <- apply(
      X = assay.data,
      MARGIN = 2,
      FUN = function(x) {
        return(x - min(x))
      }
    )
    assay.data <- t(
      x = t(x = assay.data) / apply(X = assay.data, MARGIN = 2, FUN = max)
    )
  }
  assay.data <- as.data.frame(x = assay.data)
  assay.data <- data.frame(
    cell = as.vector(x = replicate(n = length(x = features), expr = cells)),
    feature = as.vector(x = t(x = replicate(n = length(x = cells), expr = features))),
    expression = unlist(x = assay.data, use.names = FALSE)
  )
  data <- merge(x = polygons, y = assay.data, by = 'cell')
  data$feature <- factor(x = data$feature, levels = features)
  if (flip.coords) {
    coord.x <- data$x
    data$x <- data$y
    data$y <- coord.x
  }
  plot <- SinglePolyPlot(data = data, group.by = 'expression', font_size = 8) +
    scale_fill_viridis_c() +
    facet_wrap(facets = 'feature', ncol = ncol)
  return(plot)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Other plotting functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' ALRA Approximate Rank Selection Plot
#'
#' Plots the results of the approximate rank selection process for ALRA.
#'
#' @note ALRAChooseKPlot and associated functions are being moved to SeuratWrappers;
#' for more information on SeuratWrappers, please see \url{https://github.com/satijalab/seurat-wrappers}
#'
#' @param object Seurat object
#' @param start Index to start plotting singular value spacings from.
#' The transition from "signal" to "noise" in the is hard to see because the
#' first singular value spacings are so large. Nicer visualizations result from
#' skipping the first few. If set to 0 (default) starts from k/2.
#' @param combine Combine plots into a single \code{\link[patchwork]{patchwork}ed}
#' ggplot object. If \code{FALSE}, return a list of ggplot objects
#'
#' @return A list of 3 \code{\link[patchwork]{patchwork}ed} ggplot objects
#' splotting the singular values, the spacings of the singular values, and the
#' p-values of the singular values.
#'
#' @author Jun Zhao, George Linderman
#' @seealso \code{\link{RunALRA}}
#'
#' @importFrom cowplot theme_cowplot
#' @importFrom ggplot2 ggplot aes_string geom_point geom_line
#' geom_vline scale_x_continuous labs
#' @importFrom patchwork wrap_plots
#' @export
#'
ALRAChooseKPlot <- function(object, start = 0, combine = TRUE) {
  .Deprecated(
    new = 'SeruatWrappers::ALRAChooseKPlot',
    msg = paste(
      'ALRAChooseKPlot and associated functions are being moved to SeuratWrappers;',
      'for more information on SeuratWrappers, please see https://github.com/satijalab/seurat-wrappers'
    )
  )
  alra.data <- Tool(object = object, slot = 'RunALRA')
  if (is.null(x = alra.data)) {
    stop('RunALRA should be run prior to using this function.')
  }
  d <- alra.data[["d"]]
  diffs <- alra.data[["diffs"]]
  k <- alra.data[["k"]]
  if (start == 0) {
    start <- floor(x = k / 2)
  }
  if (start > k) {
    stop("Plots should include k (i.e. starting.from should be less than k)")
  }
  breaks <- seq(from = 10, to = length(x = d), by = 10)
  ggdata <- data.frame(x = 1:length(x = d), y = d)
  gg1 <- ggplot(data = ggdata, mapping = aes_string(x = 'x', y = 'y')) +
    geom_point(size = 1) +
    geom_line(size = 0.5) +
    geom_vline(xintercept = k) +
    theme_cowplot() +
    scale_x_continuous(breaks = breaks) +
    labs(x = NULL, y = 's_i', title = 'Singular values')
  ggdata <- data.frame(x = 1:(length(x = d) - 1), y = diffs)[-(1:(start - 1)), ]
  gg2 <- ggplot(data = ggdata, mapping = aes_string(x = 'x', y = 'y')) +
    geom_point(size = 1) +
    geom_line(size = 0.5) +
    geom_vline(xintercept = k + 1) +
    theme_cowplot() +
    scale_x_continuous(breaks = breaks) +
    labs(x = NULL, y = 's_{i} - s_{i-1}', title = 'Singular value spacings')
  plots <- list(spectrum = gg1, spacings = gg2)
  if (combine) {
    plots <- wrap_plots(plots)
  }
  return(plots)
}

#' Plot the Barcode Distribution and Calculated Inflection Points
#'
#' This function plots the calculated inflection points derived from the barcode-rank
#' distribution.
#'
#' See [CalculateBarcodeInflections()] to calculate inflection points and
#' [SubsetByBarcodeInflections()] to subsequently subset the Seurat object.
#'
#' @param object Seurat object
#'
#' @return Returns a `ggplot2` object showing the by-group inflection points and provided
#'   (or default) rank threshold values in grey.
#'
#' @importFrom methods slot
#' @importFrom cowplot theme_cowplot
#' @importFrom ggplot2 ggplot geom_line geom_vline aes_string
#'
#' @export
#'
#' @author Robert A. Amezquita, \email{robert.amezquita@fredhutch.org}
#' @seealso \code{\link{CalculateBarcodeInflections}} \code{\link{SubsetByBarcodeInflections}}
#'
#' @examples
#' pbmc_small <- CalculateBarcodeInflections(pbmc_small, group.column = 'groups')
#' BarcodeInflectionsPlot(pbmc_small)
#'
BarcodeInflectionsPlot <- function(object) {
  cbi.data <- Tool(object = object, slot = 'CalculateBarcodeInflections')
  if (is.null(x = cbi.data)) {
    stop("Barcode inflections not calculated, please run CalculateBarcodeInflections")
  }
  ## Extract necessary data frames
  inflection_points <- cbi.data$inflection_points
  barcode_distribution <- cbi.data$barcode_distribution
  threshold_values <- cbi.data$threshold_values
  # Set a cap to max rank to avoid plot being overextended
  if (threshold_values$rank[[2]] > max(barcode_distribution$rank, na.rm = TRUE)) {
    threshold_values$rank[[2]] <- max(barcode_distribution$rank, na.rm = TRUE)
  }
  ## Infer the grouping/barcode variables
  group_var <- colnames(x = barcode_distribution)[1]
  barcode_var <- colnames(x = barcode_distribution)[2]
  barcode_distribution[, barcode_var] <- log10(x = barcode_distribution[, barcode_var] + 1)
  ## Make the plot
  plot <- ggplot(
    data = barcode_distribution,
    mapping = aes_string(
      x = 'rank',
      y = barcode_var,
      group = group_var,
      colour = group_var
    )
  ) +
    geom_line() +
    geom_vline(
      data = threshold_values,
      aes_string(xintercept = 'rank'),
      linetype = "dashed",
      colour = 'grey60',
      size = 0.5
    ) +
    geom_vline(
      data = inflection_points,
      mapping = aes_string(
        xintercept = 'rank',
        group = group_var,
        colour = group_var
      ),
      linetype = "dashed"
    ) +
    theme_cowplot()
  return(plot)
}

#' Dot plot visualization
#'
#' Intuitive way of visualizing how feature expression changes across different
#' identity classes (clusters). The size of the dot encodes the percentage of
#' cells within a class, while the color encodes the AverageExpression level
#' across all cells within a class (blue is high).
#'
#' @param object Seurat object
#' @param assay Name of assay to use, defaults to the active assay
#' @param features Input vector of features
#' @param cols Colors to plot, can pass a single character giving the name of
#' a palette from \code{RColorBrewer::brewer.pal.info}
#' @param col.min Minimum scaled average expression threshold (everything smaller
#'  will be set to this)
#' @param col.max Maximum scaled average expression threshold (everything larger
#' will be set to this)
#' @param dot.min The fraction of cells at which to draw the smallest dot
#' (default is 0). All cell groups with less than this expressing the given
#' gene will have no dot drawn.
#' @param dot.scale Scale the size of the points, similar to cex
#' @param group.by Factor to group the cells by
#' @param split.by Factor to split the groups by (replicates the functionality of the old SplitDotPlotGG);
#' see \code{\link{FetchData}} for more details
#' @param scale.by Scale the size of the points by 'size' or by 'radius'
#' @param scale.min Set lower limit for scaling, use NA for default
#' @param scale.max Set upper limit for scaling, use NA for default
#'
#' @return A ggplot object
#'
#' @importFrom grDevices colorRampPalette
#' @importFrom cowplot theme_cowplot
#' @importFrom ggplot2 ggplot aes_string scale_size scale_radius geom_point theme element_blank labs
#' scale_color_identity scale_color_distiller scale_color_gradient guides guide_legend guide_colorbar
#' @export
#'
#' @aliases SplitDotPlotGG
#' @seealso \code{RColorBrewer::brewer.pal.info}
#'
#' @examples
#' cd_genes <- c("CD247", "CD3E", "CD9")
#' DotPlot(object = pbmc_small, features = cd_genes)
#' pbmc_small[['groups']] <- sample(x = c('g1', 'g2'), size = ncol(x = pbmc_small), replace = TRUE)
#' DotPlot(object = pbmc_small, features = cd_genes, split.by = 'groups')
#'
DotPlot <- function(
  object,
  assay = NULL,
  features,
  cols = c("lightgrey", "blue"),
  col.min = -2.5,
  col.max = 2.5,
  dot.min = 0,
  dot.scale = 6,
  group.by = NULL,
  split.by = NULL,
  scale.by = 'radius',
  scale.min = NA,
  scale.max = NA
) {
  assay <- assay %||% DefaultAssay(object = object)
  DefaultAssay(object = object) <- assay
  scale.func <- switch(
    EXPR = scale.by,
    'size' = scale_size,
    'radius' = scale_radius,
    stop("'scale.by' must be either 'size' or 'radius'")
  )
  data.features <- FetchData(object = object, vars = features)
  data.features$id <- if (is.null(x = group.by)) {
    Idents(object = object)
  } else {
    object[[group.by, drop = TRUE]]
  }
  if (!is.factor(x = data.features$id)) {
    data.features$id <- factor(x = data.features$id)
  }
  id.levels <- levels(x = data.features$id)
  data.features$id <- as.vector(x = data.features$id)
  if (!is.null(x = split.by)) {
    splits <- object[[split.by, drop = TRUE]]
    if (length(x = unique(x = splits)) > length(x = cols)) {
      stop("Not enought colors for the number of groups")
    }
    cols <- cols[1:length(x = unique(x = splits))]
    names(x = cols) <- unique(x = splits)
    data.features$id <- paste(data.features$id, splits, sep = '_')
    unique.splits <- unique(x = splits)
    id.levels <- paste0(rep(x = id.levels, each = length(x = unique.splits)), "_", rep(x = unique(x = splits), times = length(x = id.levels)))
  }
  data.plot <- lapply(
    X = unique(x = data.features$id),
    FUN = function(ident) {
      data.use <- data.features[data.features$id == ident, 1:(ncol(x = data.features) - 1), drop = FALSE]
      avg.exp <- apply(
        X = data.use,
        MARGIN = 2,
        FUN = function(x) {
          return(mean(x = expm1(x = x)))
        }
      )
      pct.exp <- apply(X = data.use, MARGIN = 2, FUN = PercentAbove, threshold = 0)
      return(list(avg.exp = avg.exp, pct.exp = pct.exp))
    }
  )
  names(x = data.plot) <- unique(x = data.features$id)
  data.plot <- lapply(
    X = names(x = data.plot),
    FUN = function(x) {
      data.use <- as.data.frame(x = data.plot[[x]])
      data.use$features.plot <- rownames(x = data.use)
      data.use$id <- x
      return(data.use)
    }
  )
  data.plot <- do.call(what = 'rbind', args = data.plot)
  if (!is.null(x = id.levels)) {
    data.plot$id <- factor(x = data.plot$id, levels = id.levels)
  }
  avg.exp.scaled <- sapply(
    X = unique(x = data.plot$features.plot),
    FUN = function(x) {
      data.use <- data.plot[data.plot$features.plot == x, 'avg.exp']
      data.use <- scale(x = data.use)
      data.use <- MinMax(data = data.use, min = col.min, max = col.max)
      return(data.use)
    }
  )
  avg.exp.scaled <- as.vector(x = t(x = avg.exp.scaled))
  if (!is.null(x = split.by)) {
    avg.exp.scaled <- as.numeric(x = cut(x = avg.exp.scaled, breaks = 20))
  }
  data.plot$avg.exp.scaled <- avg.exp.scaled
  data.plot$features.plot <- factor(
    x = data.plot$features.plot,
    levels = rev(x = features)
  )
  data.plot$pct.exp[data.plot$pct.exp < dot.min] <- NA
  data.plot$pct.exp <- data.plot$pct.exp * 100
  if (!is.null(x = split.by)) {
    splits.use <- vapply(
      X = strsplit(x = as.character(x = data.plot$id), split = '_'),
      FUN = '[[',
      FUN.VALUE = character(length = 1L),
      2
    )
    data.plot$colors <- mapply(
      FUN = function(color, value) {
        return(colorRampPalette(colors = c('grey', color))(20)[value])
      },
      color = cols[splits.use],
      value = avg.exp.scaled
    )
  }
  color.by <- ifelse(test = is.null(x = split.by), yes = 'avg.exp.scaled', no = 'colors')
  if (!is.na(x = scale.min)) {
    data.plot[data.plot$pct.exp < scale.min, 'pct.exp'] <- scale.min
  }
  if (!is.na(x = scale.max)) {
    data.plot[data.plot$pct.exp > scale.max, 'pct.exp'] <- scale.max
  }
  plot <- ggplot(data = data.plot, mapping = aes_string(x = 'features.plot', y = 'id')) +
    geom_point(mapping = aes_string(size = 'pct.exp', color = color.by)) +
    scale.func(range = c(0, dot.scale), limits = c(scale.min, scale.max)) +
    theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
    guides(size = guide_legend(title = 'Percent Expressed')) +
    labs(
      x = 'Features',
      y = ifelse(test = is.null(x = split.by), yes = 'Identity', no = 'Split Identity')
    ) +
    theme_cowplot()
  if (!is.null(x = split.by)) {
    plot <- plot + scale_color_identity()
  } else if (length(x = cols) == 1) {
    plot <- plot + scale_color_distiller(palette = cols)
  } else {
    plot <- plot + scale_color_gradient(low = cols[1], high = cols[2])
  }
  if (is.null(x = split.by)) {
    plot <- plot + guides(color = guide_colorbar(title = 'Average Expression'))
  }
  return(plot)
}

#' Quickly Pick Relevant Dimensions
#'
#' Plots the standard deviations (or approximate singular values if running PCAFast)
#' of the principle components for easy identification of an elbow in the graph.
#' This elbow often corresponds well with the significant dims and is much faster to run than
#' Jackstraw
#'
#' @param object Seurat object
#' @param ndims Number of dimensions to plot standard deviation for
#' @param reduction Reduction technique to plot standard deviation for
#'
#' @return A ggplot object
#'
#' @importFrom cowplot theme_cowplot
#' @importFrom ggplot2 ggplot aes_string geom_point labs element_line
#' @export
#'
#' @examples
#' ElbowPlot(object = pbmc_small)
#'
ElbowPlot <- function(object, ndims = 20, reduction = 'pca') {
  data.use <- Stdev(object = object, reduction = reduction)
  if (length(x = data.use) == 0) {
    stop(paste("No standard deviation info stored for", reduction))
  }
  if (ndims > length(x = data.use)) {
    warning("The object only has information for ", length(x = data.use), " reductions")
    ndims <- length(x = data.use)
  }
  stdev <- 'Standard Deviation'
  plot <- ggplot(data = data.frame(dims = 1:ndims, stdev = data.use[1:ndims])) +
    geom_point(mapping = aes_string(x = 'dims', y = 'stdev')) +
    labs(
      x = gsub(
        pattern = '_$',
        replacement = '',
        x = Key(object = object[[reduction]])
      ),
      y = stdev
    ) +
    theme_cowplot()
  return(plot)
}

#' JackStraw Plot
#'
#' Plots the results of the JackStraw analysis for PCA significance. For each
#' PC, plots a QQ-plot comparing the distribution of p-values for all genes
#' across each PC, compared with a uniform distribution. Also determines a
#' p-value for the overall significance of each PC (see Details).
#'
#' Significant PCs should show a p-value distribution (black curve) that is
#' strongly skewed to the left compared to the null distribution (dashed line)
#' The p-value for each PC is based on a proportion test comparing the number
#' of genes with a p-value below a particular threshold (score.thresh), compared with the
#' proportion of genes expected under a uniform distribution of p-values.
#'
#' @param object Seurat object
#' @param dims Dims to plot
#' @param reduction reduction to pull jackstraw info from
#' @param xmax X-axis maximum on each QQ plot.
#' @param ymax Y-axis maximum on each QQ plot.
#'
#' @return A ggplot object
#'
#' @author Omri Wurtzel
#' @seealso \code{\link{ScoreJackStraw}}
#'
#' @importFrom stats qunif
#' @importFrom ggplot2 ggplot aes_string stat_qq labs xlim ylim
#' coord_flip geom_abline guides guide_legend
#' @importFrom cowplot theme_cowplot
#'
#' @export
#'
#' @examples
#' JackStrawPlot(object = pbmc_small)
#'
JackStrawPlot <- function(
  object,
  dims = 1:5,
  reduction = 'pca',
  xmax = 0.1,
  ymax = 0.3
) {
  pAll <- JS(object = object[[reduction]], slot = 'empirical')
  if (max(dims) > ncol(x = pAll)) {
    stop("Max dimension is ", ncol(x = pAll))
  }
  pAll <- pAll[, dims, drop = FALSE]
  pAll <- as.data.frame(x = pAll)
  data.plot <- Melt(x = pAll)
  colnames(x = data.plot) <- c("Contig", "PC", "Value")
  score.df <- JS(object = object[[reduction]], slot = 'overall')
  if (nrow(x = score.df) < max(dims)) {
    stop("Jackstraw procedure not scored for all the provided dims. Please run ScoreJackStraw.")
  }
  score.df <- score.df[dims, , drop = FALSE]
  if (nrow(x = score.df) == 0) {
    stop(paste0("JackStraw hasn't been scored. Please run ScoreJackStraw before plotting."))
  }
  data.plot$PC.Score <- rep(
    x = paste0("PC ", score.df[ ,"PC"], ": ", sprintf("%1.3g", score.df[ ,"Score"])),
    each = length(x = unique(x = data.plot$Contig))
  )
  data.plot$PC.Score <- factor(
    x = data.plot$PC.Score,
    levels = paste0("PC ", score.df[, "PC"], ": ", sprintf("%1.3g", score.df[, "Score"]))
  )
  gp <- ggplot(data = data.plot, mapping = aes_string(sample = 'Value', color = 'PC.Score')) +
    stat_qq(distribution = qunif) +
    labs(x = "Theoretical [runif(1000)]", y = "Empirical") +
    xlim(0, ymax) +
    ylim(0, xmax) +
    coord_flip() +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", na.rm = TRUE) +
    guides(color = guide_legend(title = "PC: p-value")) +
    theme_cowplot()
  return(gp)
}

#' Plot clusters as a tree
#'
#' Plots previously computed tree (from BuildClusterTree)
#'
#' @param object Seurat object
#' @param \dots Additional arguments to ape::plot.phylo
#'
#' @return Plots dendogram (must be precomputed using BuildClusterTree), returns no value
#'
#' @importFrom ape plot.phylo
#' @importFrom ape nodelabels
#'
#' @export
#'
#' @examples
#' pbmc_small <- BuildClusterTree(object = pbmc_small)
#' PlotClusterTree(object = pbmc_small)
#'
PlotClusterTree <- function(object, ...) {
  if (is.null(x = Tool(object = object, slot = "BuildClusterTree"))) {
    stop("Phylogenetic tree does not exist, build using BuildClusterTree")
  }
  data.tree <- Tool(object = object, slot = "BuildClusterTree")
  plot.phylo(x = data.tree, direction = "downwards", ...)
  nodelabels()
}

#' Visualize Dimensional Reduction genes
#'
#' Visualize top genes associated with reduction components
#'
#' @param object Seurat object
#' @param reduction Reduction technique to visualize results for
#' @param dims Number of dimensions to display
#' @param nfeatures Number of genes to display
#' @param col Color of points to use
#' @param projected Use reduction values for full dataset (i.e. projected
#' dimensional reduction values)
#' @param balanced Return an equal number of genes with + and - scores. If
#' FALSE (default), returns the top genes ranked by the scores absolute values
#' @param ncol Number of columns to display
#' @param combine Combine plots into a single \code{\link[patchwork]{patchwork}ed}
#' ggplot object. If \code{FALSE}, return a list of ggplot objects
#'
#' @return A \code{\link[patchwork]{patchwork}ed} ggplot object if
#' \code{combine = TRUE}; otherwise, a list of ggplot objects
#'
#' @importFrom patchwork wrap_plots
#' @importFrom cowplot theme_cowplot
#' @importFrom ggplot2 ggplot aes_string geom_point labs
#' @export
#'
#' @examples
#' VizDimLoadings(object = pbmc_small)
#'
VizDimLoadings <- function(
  object,
  dims = 1:5,
  nfeatures = 30,
  col = 'blue',
  reduction = 'pca',
  projected = FALSE,
  balanced = FALSE,
  ncol = NULL,
  combine = TRUE
) {
  ncol <- ncol %||% 2
  if (length(x = dims) == 1) {
    ncol <- 1
  }
  if (length(x = dims) > 6) {
    ncol <- 3
  }
  if (length(x = dims) > 9) {
    ncol <- 4
  }
  loadings <- Loadings(object = object[[reduction]], projected = projected)
  features <- lapply(
    X = dims,
    FUN = TopFeatures,
    object = object[[reduction]],
    nfeatures = nfeatures,
    projected = projected,
    balanced = balanced
  )
  features <- lapply(
    X = features,
    FUN = unlist,
    use.names = FALSE
  )
  loadings <- loadings[unlist(x = features), dims, drop = FALSE]
  names(x = features) <- colnames(x = loadings) <- as.character(x = dims)
  plots <- lapply(
    X = as.character(x = dims),
    FUN = function(i) {
      data.plot <- as.data.frame(x = loadings[features[[i]], i, drop = FALSE])
      colnames(x = data.plot) <- paste0(Key(object = object[[reduction]]), i)
      data.plot$feature <- factor(x = rownames(x = data.plot), levels = rownames(x = data.plot))
      plot <- ggplot(
        data = data.plot,
        mapping = aes_string(x = colnames(x = data.plot)[1], y = 'feature')
      ) +
        geom_point(col = col) +
        labs(y = NULL) + theme_cowplot()
      return(plot)
    }
  )
  if (combine) {
    plots <- wrap_plots(plots, ncol = ncol)
  }
  return(plots)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Exported utility functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Augments ggplot2-based plot with a PNG image.
#'
#' Creates "vector-friendly" plots. Does this by saving a copy of the plot as a PNG file,
#' then adding the PNG image with \code{\link[ggplot2]{annotation_raster}} to a blank plot
#' of the same dimensions as \code{plot}. Please note: original legends and axes will be lost
#' during augmentation.
#'
#' @param plot A ggplot object
#' @param width,height Width and height of PNG version of plot
#' @param dpi Plot resolution
#'
#' @return A ggplot object
#'
#' @importFrom png readPNG
#' @importFrom ggplot2 ggplot_build ggsave ggplot aes_string geom_blank annotation_raster ggtitle
#'
#' @export
#'
#' @examples
#' \dontrun{
#' plot <- DimPlot(object = pbmc_small)
#' AugmentPlot(plot = plot)
#' }
#'
AugmentPlot <- function(plot, width = 10, height = 10, dpi = 100) {
  pbuild.params <- ggplot_build(plot = plot)$layout$panel_params[[1]]
  range.values <- c(
    pbuild.params$x.range,
    pbuild.params$y.range
  )
  xyparams <- GetXYAesthetics(
    plot = plot,
    geom = class(x = plot$layers[[1]]$geom)[1]
  )
  title <- plot$labels$title
  tmpfile <- tempfile(fileext = '.png')
  ggsave(
    filename = tmpfile,
    plot = plot + NoLegend() + NoAxes() + theme(plot.title = element_blank()),
    width = width,
    height = height,
    dpi = dpi
  )
  img <- readPNG(source = tmpfile)
  file.remove(tmpfile)
  blank <- ggplot(
    data = plot$data,
    mapping = aes_string(x = xyparams$x, y = xyparams$y)
  ) + geom_blank()
  blank <- blank + plot$theme + ggtitle(label = title)
  blank <- blank + annotation_raster(
    raster = img,
    xmin = range.values[1],
    xmax = range.values[2],
    ymin = range.values[3],
    ymax = range.values[4]
  )
  return(blank)
}

#' @inheritParams CustomPalette
#'
#' @export
#'
#' @rdname CustomPalette
#' @aliases BlackAndWhite
#'
#' @examples
#' df <- data.frame(x = rnorm(n = 100, mean = 20, sd = 2), y = rbinom(n = 100, size = 100, prob = 0.2))
#' plot(df, col = BlackAndWhite())
#'
BlackAndWhite <- function(mid = NULL, k = 50) {
  return(CustomPalette(low = "white", high = "black", mid = mid, k = k))
}

#' @inheritParams CustomPalette
#'
#' @export
#'
#' @rdname CustomPalette
#' @aliases BlueAndRed
#'
#' @examples
#' df <- data.frame(x = rnorm(n = 100, mean = 20, sd = 2), y = rbinom(n = 100, size = 100, prob = 0.2))
#' plot(df, col = BlueAndRed())
#'
BlueAndRed <- function(k = 50) {
  return(CustomPalette(low = "#313695" , high = "#A50026", mid = "#FFFFBF", k = k))
}

#' Cell selector
#'
#' Select points on a scatterplot and get information about them
#'
#' @param plot A ggplot2 plot
#' @param object An optional Seurat object; if passes, will return an object with
#' the identities of selected cells set to \code{ident}
#' @param ident An optional new identity class to assign the selected cells
#' @param ... Extra parameters, such as dark.theme, recolor, or smooth for using a dark theme,
#' recoloring based on selected cells, or using a smooth scatterplot, respectively
#'
#' @return If \code{object} is \code{NULL}, the names of the points selected; otherwise,
#' a Seurat object with the selected cells identity classes set to \code{ident}
#'
#' @importFrom ggplot2 ggplot_build
#' @export
#'
# @aliases FeatureLocator
#' @seealso \code{\link[graphics]{locator}} \code{\link[ggplot2]{ggplot_build}}
#' \code{\link[SDMTools]{pnt.in.poly}} \code{\link{DimPlot}} \code{\link{FeaturePlot}}
#'
#' @examples
#' \dontrun{
#' plot <- DimPlot(object = pbmc_small)
#' # Follow instructions in the terminal to select points
#' cells.located <- CellSelector(plot = plot)
#' cells.located
#' # Automatically set the identity class of selected cells and return a new Seurat object
#' pbmc_small <- CellSelector(plot = plot, object = pbmc_small, ident = 'SelectedCells')
#' }
#'
CellSelector <- function(plot, object = NULL, ident = 'SelectedCells', ...) {
  located <- PointLocator(plot = plot, ...)
  data <- ggplot_build(plot = plot)$plot$data
  selected <- rownames(x = data[as.numeric(x = rownames(x = located)), ])
  if (inherits(x = object, what = 'Seurat')) {
    if (!all(selected %in% Cells(x = object))) {
      stop("Cannot find selected cells in the Seurat object, please be sure you pass the same object used to generate the plot", call. = FALSE)
    }
    Idents(object = object, cells = selected) <- ident
    return(object)
  }
  return(selected)
}

#' Move outliers towards center on dimension reduction plot
#'
#' @param object Seurat object
#' @param reduction Name of DimReduc to adjust
#' @param dims Dimensions to visualize
#' @param group.by Group (color) cells in different ways (for example, orig.ident)
#' @param outlier.sd Controls the outlier distance
#' @param reduction.key Key for DimReduc that is returned
#'
#' @return Returns a DimReduc object with the modified embeddings
#'
#' @export
#'
#' @examples
#' \dontrun{
#' pbmc_small <- FindClusters(pbmc_small, resolution = 1.1)
#' pbmc_small <- RunUMAP(pbmc_small, dims = 1:5)
#' DimPlot(pbmc_small, reduction = "umap")
#' pbmc_small[["umap_new"]] <- CollapseEmbeddingOutliers(pbmc_small,
#'     reduction = "umap", reduction.key = 'umap_', outlier.sd = 0.5)
#' DimPlot(pbmc_small, reduction = "umap_new")
#' }
#'
CollapseEmbeddingOutliers <- function(
  object,
  reduction = 'umap',
  dims = 1:2,
  group.by = 'ident',
  outlier.sd = 2,
  reduction.key = 'UMAP_'
) {
  embeddings <- Embeddings(object = object[[reduction]])[, dims]
  idents <- FetchData(object = object, vars = group.by)
  data.medians <- sapply(X = dims, FUN = function(x) {
    tapply(X = embeddings[, x], INDEX = idents, FUN = median)
  })
  data.sd <- apply(X = data.medians, MARGIN = 2, FUN = sd)
  data.medians.scale <- as.matrix(x = scale(x = data.medians, center = TRUE, scale = TRUE))
  data.medians.scale[abs(x = data.medians.scale) < outlier.sd] <- 0
  data.medians.scale <- sign(x = data.medians.scale) * (abs(x = data.medians.scale) - outlier.sd)
  data.correct <- Sweep(
    x = data.medians.scale,
    MARGIN = 2,
    STATS = data.sd,
    FUN = "*"
  )
  data.correct <- data.correct[abs(x = apply(X = data.correct, MARGIN = 1, FUN = min)) > 0, ]
  new.embeddings <- embeddings
  for (i in rownames(x = data.correct)) {
    cells.correct <- rownames(x = idents)[idents[, "ident"] == i]
    new.embeddings[cells.correct, ] <- Sweep(
      x = new.embeddings[cells.correct,],
      MARGIN = 2,
      STATS = data.correct[i, ],
      FUN = "-"
    )
  }
  reduc <- CreateDimReducObject(
    embeddings = new.embeddings,
    loadings = Loadings(object = object[[reduction]]),
    assay = slot(object = object[[reduction]], name = "assay.used"),
    key = reduction.key
  )
  return(reduc)
}

#' Combine ggplot2-based plots into a single plot
#'
#' @param plots A list of gg objects
#' @param ncol Number of columns
#' @param legend Combine legends into a single legend
#' choose from 'right' or 'bottom'; pass 'none' to remove legends, or \code{NULL}
#' to leave legends as they are
#' @param ... Extra parameters passed to plot_grid
#'
#' @return A combined plot
#'
#' @importFrom cowplot plot_grid get_legend
#' @export
#'
#' @examples
#' pbmc_small[['group']] <- sample(
#'   x = c('g1', 'g2'),
#'   size = ncol(x = pbmc_small),
#'   replace = TRUE
#' )
#' plot1 <- FeaturePlot(
#'   object = pbmc_small,
#'   features = 'MS4A1',
#'   split.by = 'group'
#' )
#' plot2 <- FeaturePlot(
#'   object = pbmc_small,
#'   features = 'FCN1',
#'   split.by = 'group'
#' )
#' CombinePlots(
#'   plots = list(plot1, plot2),
#'   legend = 'none',
#'   nrow = length(x = unique(x = pbmc_small[['group', drop = TRUE]]))
#' )
#'
CombinePlots <- function(plots, ncol = NULL, legend = NULL, ...) {
  .Deprecated(msg = "CombinePlots is being deprecated. Plots should now be combined using the patchwork system.")
  plots.combined <- if (length(x = plots) > 1) {
    if (!is.null(x = legend)) {
      if (legend != 'none') {
        plot.legend <- get_legend(plot = plots[[1]] + theme(legend.position = legend))
      }
      plots <- lapply(
        X = plots,
        FUN = function(x) {
          return(x + NoLegend())
        }
      )
    }
    plots.combined <- plot_grid(
      plotlist = plots,
      ncol = ncol,
      align = 'hv',
      ...
    )
    if (!is.null(x = legend)) {
      plots.combined <- switch(
        EXPR = legend,
        'bottom' = plot_grid(
          plots.combined,
          plot.legend,
          ncol = 1,
          rel_heights = c(1, 0.2)
        ),
        'right' = plot_grid(
          plots.combined,
          plot.legend,
          rel_widths = c(3, 0.3)
        ),
        plots.combined
      )
    }
    plots.combined
  } else {
    plots[[1]]
  }
  return(plots.combined)
}

#' Create a custom color palette
#'
#' Creates a custom color palette based on low, middle, and high color values
#'
#' @param low low color
#' @param high high color
#' @param mid middle color. Optional.
#' @param k number of steps (colors levels) to include between low and high values
#'
#' @return A color palette for plotting
#'
#' @importFrom grDevices col2rgb rgb
#' @export
#'
#' @rdname CustomPalette
#' @examples
#' myPalette <- CustomPalette()
#' myPalette
#'
CustomPalette <- function(
  low = "white",
  high = "red",
  mid = NULL,
  k = 50
) {
  low <- col2rgb(col = low) / 255
  high <- col2rgb(col = high) / 255
  if (is.null(x = mid)) {
    r <- seq(from = low[1], to = high[1], len = k)
    g <- seq(from = low[2], to = high[2], len = k)
    b <- seq(from = low[3], to = high[3], len = k)
  } else {
    k2 <- round(x = k / 2)
    mid <- col2rgb(col = mid) / 255
    r <- c(
      seq(from = low[1], to = mid[1], len = k2),
      seq(from = mid[1], to = high[1], len = k2)
    )
    g <- c(
      seq(from = low[2], to = mid[2], len = k2),
      seq(from = mid[2], to = high[2],len = k2)
    )
    b <- c(
      seq(from = low[3], to = mid[3], len = k2),
      seq(from = mid[3], to = high[3], len = k2)
    )
  }
  return(rgb(red = r, green = g, blue = b))
}

#' Discrete colour palettes from the pals package
#'
#' These are included here because pals depends on a number of compiled
#' packages, and this can lead to increases in run time for Travis,
#' and generally should be avoided when possible.
#'
#' These palettes are a much better default for data with many classes
#' than the default ggplot2 palette.
#'
#' Many thanks to Kevin Wright for writing the pals package.
#'
#' @param n Number of colours to be generated.
#' @param palette Options are
#' "alphabet", "alphabet2", "glasbey", "polychrome", and "stepped".
#' Can be omitted and the function will use the one based on the requested n.
#'
#' @return A vector of colors
#'
#' @details
#' Taken from the pals package (Licence: GPL-3).
#' \url{https://cran.r-project.org/package=pals}
#' Credit: Kevin Wright
#'
#' @export
#'
DiscretePalette <- function(n, palette = NULL) {
  palettes <- list(
    alphabet = c(
      "#F0A0FF", "#0075DC", "#993F00", "#4C005C", "#191919", "#005C31",
      "#2BCE48", "#FFCC99", "#808080", "#94FFB5", "#8F7C00", "#9DCC00",
      "#C20088", "#003380", "#FFA405", "#FFA8BB", "#426600", "#FF0010",
      "#5EF1F2", "#00998F", "#E0FF66", "#740AFF", "#990000", "#FFFF80",
      "#FFE100", "#FF5005"
    ),
    alphabet2 = c(
      "#AA0DFE", "#3283FE", "#85660D", "#782AB6", "#565656", "#1C8356",
      "#16FF32", "#F7E1A0", "#E2E2E2", "#1CBE4F", "#C4451C", "#DEA0FD",
      "#FE00FA", "#325A9B", "#FEAF16", "#F8A19F", "#90AD1C", "#F6222E",
      "#1CFFCE", "#2ED9FF", "#B10DA1", "#C075A6", "#FC1CBF", "#B00068",
      "#FBE426", "#FA0087"
    ),
    glasbey = c(
      "#0000FF", "#FF0000", "#00FF00", "#000033", "#FF00B6", "#005300",
      "#FFD300", "#009FFF", "#9A4D42", "#00FFBE", "#783FC1", "#1F9698",
      "#FFACFD", "#B1CC71", "#F1085C", "#FE8F42", "#DD00FF", "#201A01",
      "#720055", "#766C95", "#02AD24", "#C8FF00", "#886C00", "#FFB79F",
      "#858567", "#A10300", "#14F9FF", "#00479E", "#DC5E93", "#93D4FF",
      "#004CFF", "#F2F318"
    ),
    polychrome = c(
      "#5A5156", "#E4E1E3", "#F6222E", "#FE00FA", "#16FF32", "#3283FE",
      "#FEAF16", "#B00068", "#1CFFCE", "#90AD1C", "#2ED9FF", "#DEA0FD",
      "#AA0DFE", "#F8A19F", "#325A9B", "#C4451C", "#1C8356", "#85660D",
      "#B10DA1", "#FBE426", "#1CBE4F", "#FA0087", "#FC1CBF", "#F7E1A0",
      "#C075A6", "#782AB6", "#AAF400", "#BDCDFF", "#822E1C", "#B5EFB5",
      "#7ED7D1", "#1C7F93", "#D85FF7", "#683B79", "#66B0FF", "#3B00FB"
    ),
    stepped = c(
      "#990F26", "#B33E52", "#CC7A88", "#E6B8BF", "#99600F", "#B3823E",
      "#CCAA7A", "#E6D2B8", "#54990F", "#78B33E", "#A3CC7A", "#CFE6B8",
      "#0F8299", "#3E9FB3", "#7ABECC", "#B8DEE6", "#3D0F99", "#653EB3",
      "#967ACC", "#C7B8E6", "#333333", "#666666", "#999999", "#CCCCCC"
    )
  )
  if (is.null(x = palette)) {
    if (n <= 26) {
      palette <- "alphabet"
    } else if (n <= 32) {
      palette <- "glasbey"
    } else {
      palette <- "polychrome"
    }
  }
  palette.vec <- palettes[[palette]]
  if (n > length(x = palette.vec)) {
    warning("Not enough colours in specified palette")
  }
  palette.vec[seq_len(length.out = n)]
}

#' @rdname CellSelector
#' @export
#'
FeatureLocator <- function(plot, ...) {
  .Defunct(
    new = 'CellSelector',
    package = 'Seurat',
    msg = "'FeatureLocator' has been replaced by 'CellSelector'"
  )
}

#' Hover Locator
#'
#' Get quick information from a scatterplot by hovering over points
#'
#' @param plot A ggplot2 plot
#' @param information An optional dataframe or matrix of extra information to be displayed on hover
#' @param dark.theme Plot using a dark theme?
#' @param ... Extra parameters to be passed to \code{plotly::layout}
#'
#' @importFrom ggplot2 ggplot_build
#' @importFrom plotly plot_ly layout
#' @export
#'
#' @seealso \code{\link[plotly]{layout}} \code{\link[ggplot2]{ggplot_build}}
#' \code{\link{DimPlot}} \code{\link{FeaturePlot}}
#'
#' @examples
#' \dontrun{
#' plot <- DimPlot(object = pbmc_small)
#' HoverLocator(plot = plot, information = FetchData(object = pbmc_small, vars = 'percent.mito'))
#' }
#'
HoverLocator <- function(
  plot,
  information = NULL,
  dark.theme = FALSE,
  ...
) {
  #   Use GGpointToBase because we already have ggplot objects
  #   with colors (which are annoying in plotly)
  plot.build <- GGpointToBase(plot = plot, do.plot = FALSE)
  data <- ggplot_build(plot = plot)$plot$data
  rownames(x = plot.build) <- rownames(x = data)
  #   Reset the names to 'x' and 'y'
  names(x = plot.build) <- c(
    'x',
    'y',
    names(x = plot.build)[3:length(x = plot.build)]
  )
  #   Add the names we're looking for (eg. cell name, gene name)
  if (is.null(x = information)) {
    plot.build$feature <- rownames(x = data)
  } else {
    info <- apply(
      X = information,
      MARGIN = 1,
      FUN = function(x, names) {
        return(paste0(names, ': ', x, collapse = '<br>'))
      },
      names = colnames(x = information)
    )
    data.info <- data.frame(
      feature = paste(rownames(x = information), info, sep = '<br>'),
      row.names = rownames(x = information)
    )
    plot.build <- merge(x = plot.build, y = data.info, by = 0)
  }
  #   Set up axis labels here
  #   Also, a bunch of stuff to get axis lines done properly
  xaxis <- list(
    title = names(x = data)[1],
    showgrid = FALSE,
    zeroline = FALSE,
    showline = TRUE
  )
  yaxis <- list(
    title = names(x = data)[2],
    showgrid = FALSE,
    zeroline = FALSE,
    showline = TRUE
  )
  #   Check for dark theme
  if (dark.theme) {
    title <- list(color = 'white')
    xaxis <- c(xaxis, color = 'white')
    yaxis <- c(yaxis, color = 'white')
    plotbg <- 'black'
  } else {
    title = list(color = 'black')
    plotbg = 'white'
  }
  #   The `~' means pull from the data passed (this is why we reset the names)
  #   Use I() to get plotly to accept the colors from the data as is
  #   Set hoverinfo to 'text' to override the default hover information
  #   rather than append to it
  p <- plotly::layout(
    p = plot_ly(
      data = plot.build,
      x = ~x,
      y = ~y,
      type = 'scatter',
      mode = 'markers',
      color = ~I(color),
      hoverinfo = 'text',
      text = ~feature
    ),
    xaxis = xaxis,
    yaxis = yaxis,
    title = plot$labels$title,
    titlefont = title,
    paper_bgcolor = plotbg,
    plot_bgcolor = plotbg,
    ...
  )
  # add labels
  label.layer <- which(x = sapply(
    X = plot$layers,
    FUN = function(x) class(x$geom)[1] == "GeomText")
  )
  if (length(x = label.layer) == 1) {
    p <- plotly::add_annotations(
      p = p,
      x = plot$layers[[label.layer]]$data[, 1],
      y = plot$layers[[label.layer]]$data[, 2],
      xref = "x",
      yref = "y",
      text = plot$layers[[label.layer]]$data[, 3],
      xanchor = 'right',
      showarrow = FALSE,
      font = list(size = plot$layers[[label.layer]]$aes_params$size * 4)
    )
  }
  return(p)
}

#' Label clusters on a ggplot2-based scatter plot
#'
#' @param plot A ggplot2-based scatter plot
#' @param id Name of variable used for coloring scatter plot
#' @param clusters Vector of cluster ids to label
#' @param labels Custom labels for the clusters
#' @param split.by Split labels by some grouping label, useful when using
#' \code{\link[ggplot2]{facet_wrap}} or \code{\link[ggplot2]{facet_grid}}
#' @param repel Use \code{geom_text_repel} to create nicely-repelled labels
#' @param ... Extra parameters to \code{\link[ggrepel]{geom_text_repel}}, such as \code{size}
#'
#' @return A ggplot2-based scatter plot with cluster labels
#'
#' @importFrom stats median
#' @importFrom ggrepel geom_text_repel
#' @importFrom ggplot2 aes_string geom_text
#' @export
#'
#' @seealso \code{\link[ggrepel]{geom_text_repel}} \code{\link[ggplot2]{geom_text}}
#'
#' @examples
#' plot <- DimPlot(object = pbmc_small)
#' LabelClusters(plot = plot, id = 'ident')
#'
LabelClusters <- function(
  plot,
  id,
  clusters = NULL,
  labels = NULL,
  split.by = NULL,
  repel = TRUE,
  ...
) {
  xynames <- unlist(x = GetXYAesthetics(plot = plot), use.names = TRUE)
  if (!id %in% colnames(x = plot$data)) {
    stop("Cannot find variable ", id, " in plotting data")
  }
  if (!is.null(x = split.by) && !split.by %in% colnames(x = plot$data)) {
    warning("Cannot find splitting variable ", id, " in plotting data")
    split.by <- NULL
  }
  data <- plot$data[, c(xynames, id, split.by)]
  possible.clusters <- as.character(x = na.omit(object = unique(x = data[, id])))
  groups <- clusters %||% as.character(x = na.omit(object = unique(x = data[, id])))
  if (any(!groups %in% possible.clusters)) {
    stop("The following clusters were not found: ", paste(groups[!groups %in% possible.clusters], collapse = ","))
  }
  labels.loc <- lapply(
    X = groups,
    FUN = function(group) {
      data.use <- data[data[, id] == group, , drop = FALSE]
      data.medians <- if (!is.null(x = split.by)) {
        do.call(
          what = 'rbind',
          args = lapply(
            X = unique(x = data.use[, split.by]),
            FUN = function(split) {
              medians <- apply(
                X = data.use[data.use[, split.by] == split, xynames, drop = FALSE],
                MARGIN = 2,
                FUN = median,
                na.rm = TRUE
              )
              medians <- as.data.frame(x = t(x = medians))
              medians[, split.by] <- split
              return(medians)
            }
          )
        )
      } else {
        as.data.frame(x = t(x = apply(
          X = data.use[, xynames, drop = FALSE],
          MARGIN = 2,
          FUN = median,
          na.rm = TRUE
        )))
      }
      data.medians[, id] <- group
      return(data.medians)
    }
  )
  labels.loc <- do.call(what = 'rbind', args = labels.loc)
  labels <- labels %||% groups
  if (length(x = unique(x = labels.loc[, id])) != length(x = labels)) {
    stop("Length of labels (", length(x = labels),  ") must be equal to the number of clusters being labeled (", length(x = labels.loc), ").")
  }
  names(x = labels) <- groups
  for (group in groups) {
    labels.loc[labels.loc[, id] == group, id] <- labels[group]
  }
  geom.use <- ifelse(test = repel, yes = geom_text_repel, no = geom_text)
  plot <- plot + geom.use(
    data = labels.loc,
    mapping = aes_string(x = xynames['x'], y = xynames['y'], label = id),
    ...
  )
  return(plot)
}

#' Add text labels to a ggplot2 plot
#'
#' @param plot A ggplot2 plot with a GeomPoint layer
#' @param points A vector of points to label; if \code{NULL}, will use all points in the plot
#' @param labels A vector of labels for the points; if \code{NULL}, will use
#' rownames of the data provided to the plot at the points selected
#' @param repel Use \code{geom_text_repel} to create a nicely-repelled labels; this
#' is slow when a lot of points are being plotted. If using \code{repel}, set \code{xnudge}
#' and \code{ynudge} to 0
#' @param xnudge,ynudge Amount to nudge X and Y coordinates of labels by
#' @param ... Extra parameters passed to \code{geom_text}
#'
#' @return A ggplot object
#'
#' @importFrom ggrepel geom_text_repel
#' @importFrom ggplot2 geom_text aes_string
#' @export
#'
#' @aliases Labeler
#' @seealso \code{\link[ggplot2]{geom_text}}
#'
#' @examples
#' ff <- TopFeatures(object = pbmc_small[['pca']])
#' cc <- TopCells(object = pbmc_small[['pca']])
#' plot <- FeatureScatter(object = pbmc_small, feature1 = ff[1], feature2 = ff[2])
#' LabelPoints(plot = plot, points = cc)
#'
LabelPoints <- function(
  plot,
  points,
  labels = NULL,
  repel = FALSE,
  xnudge = 0.3,
  ynudge = 0.05,
  ...
) {
  xynames <- GetXYAesthetics(plot = plot)
  points <- points %||% rownames(x = plot$data)
  if (is.numeric(x = points)) {
    points <- rownames(x = plot$data)
  }
  points <- intersect(x = points, y = rownames(x = plot$data))
  if (length(x = points) == 0) {
    stop("Cannot find points provided")
  }
  labels <- labels %||% points
  labels <- as.character(x = labels)
  label.data <- plot$data[points, ]
  label.data$labels <- labels
  geom.use <- ifelse(test = repel, yes = geom_text_repel, no = geom_text)
  if (repel) {
    if (!all(c(xnudge, ynudge) == 0)) {
      message("When using repel, set xnudge and ynudge to 0 for optimal results")
    }
  }
  plot <- plot + geom.use(
    mapping = aes_string(x = xynames$x, y = xynames$y, label = 'labels'),
    data = label.data,
    nudge_x = xnudge,
    nudge_y = ynudge,
    ...
  )
  return(plot)
}

#' @inheritParams CustomPalette
#'
#' @export
#'
#' @rdname CustomPalette
#' @aliases PurpleAndYellow
#'
#' @examples
#' df <- data.frame(x = rnorm(n = 100, mean = 20, sd = 2), y = rbinom(n = 100, size = 100, prob = 0.2))
#' plot(df, col = PurpleAndYellow())
#'
PurpleAndYellow <- function(k = 50) {
  return(CustomPalette(low = "magenta", high = "yellow", mid = "black", k = k))
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Seurat themes
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Seurat Themes
#'
#' Various themes to be applied to ggplot2-based plots
#' \describe{
#'   \item{\code{SeuratTheme}}{The curated Seurat theme, consists of ...}
#'   \item{\code{DarkTheme}}{A dark theme, axes and text turn to white, the background becomes black}
#'   \item{\code{NoAxes}}{Removes axis lines, text, and ticks}
#'   \item{\code{NoLegend}}{Removes the legend}
#'   \item{\code{FontSize}}{Sets axis and title font sizes}
#'   \item{\code{NoGrid}}{Removes grid lines}
#'   \item{\code{SeuratAxes}}{Set Seurat-style axes}
#'   \item{\code{SpatialTheme}}{A theme designed for spatial visualizations (eg \code{\link{PolyFeaturePlot}}, \code{\link{PolyDimPlot}})}
#'   \item{\code{RestoreLegend}}{Restore a legend after removal}
#'   \item{\code{RotatedAxis}}{Rotate X axis text 45 degrees}
#'   \item{\code{BoldTitle}}{Enlarges and emphasizes the title}
#' }
#'
#' @param ... Extra parameters to be passed to \code{theme}
#'
#' @return A ggplot2 theme object
#'
#' @export
#'
#' @rdname SeuratTheme
#' @seealso \code{\link[ggplot2]{theme}}
#' @aliases SeuratTheme
#'
SeuratTheme <- function() {
  return(DarkTheme() + NoLegend() + NoGrid() + SeuratAxes())
}

#' @inheritParams SeuratTheme
#'
#' @importFrom ggplot2 theme element_rect element_text element_line margin
#' @export
#'
#' @rdname SeuratTheme
#' @aliases DarkTheme
#'
#' @examples
#' # Generate a plot with a dark theme
#' library(ggplot2)
#' df <- data.frame(x = rnorm(n = 100, mean = 20, sd = 2), y = rbinom(n = 100, size = 100, prob = 0.2))
#' p <- ggplot(data = df, mapping = aes(x = x, y = y)) + geom_point(mapping = aes(color = 'red'))
#' p + DarkTheme(legend.position = 'none')
#'
DarkTheme <- function(...) {
  #   Some constants for easier changing in the future
  black.background <- element_rect(fill = 'black')
  black.background.no.border <- element_rect(fill = 'black', size = 0)
  font.margin <- 4
  white.text <- element_text(
    colour = 'white',
    margin = margin(
      t = font.margin,
      r = font.margin,
      b = font.margin,
      l = font.margin
    )
  )
  white.line <- element_line(colour = 'white', size = 1)
  no.line <- element_line(size = 0)
  #   Create the dark theme
  dark.theme <- theme(
    #   Set background colors
    plot.background = black.background,
    panel.background = black.background,
    legend.background = black.background,
    legend.box.background = black.background.no.border,
    legend.key = black.background.no.border,
    strip.background = element_rect(fill = 'grey50', colour = NA),
    #   Set text colors
    plot.title = white.text,
    plot.subtitle = white.text,
    axis.title = white.text,
    axis.text = white.text,
    legend.title = white.text,
    legend.text = white.text,
    strip.text = white.text,
    #   Set line colors
    axis.line.x = white.line,
    axis.line.y = white.line,
    panel.grid = no.line,
    panel.grid.minor = no.line,
    #   Validate the theme
    validate = TRUE,
    #   Extra parameters
    ...
  )
  return(dark.theme)
}

#' @inheritParams SeuratTheme
#' @param x.text,y.text X and Y axis text sizes
#' @param x.title,y.title X and Y axis title sizes
#' @param main Plot title size
#'
#' @importFrom ggplot2 theme element_text
#' @export
#'
#' @rdname SeuratTheme
#' @aliases FontSize
#'
FontSize <- function(
  x.text = NULL,
  y.text = NULL,
  x.title = NULL,
  y.title = NULL,
  main = NULL,
  ...
) {
  font.size <- theme(
    # Set font sizes
    axis.text.x = element_text(size = x.text),
    axis.text.y = element_text(size = y.text),
    axis.title.x = element_text(size = x.title),
    axis.title.y = element_text(size = y.title),
    plot.title = element_text(size = main),
    # Validate the theme
    validate = TRUE,
    # Extra parameters
    ...
  )
}

#' @inheritParams SeuratTheme
#' @param keep.text Keep axis text
#' @param keep.ticks Keep axis ticks
#'
#' @importFrom ggplot2 theme element_blank
#' @export
#'
#' @rdname SeuratTheme
#' @aliases NoAxes
#'
#' @examples
#' # Generate a plot with no axes
#' library(ggplot2)
#' df <- data.frame(x = rnorm(n = 100, mean = 20, sd = 2), y = rbinom(n = 100, size = 100, prob = 0.2))
#' p <- ggplot(data = df, mapping = aes(x = x, y = y)) + geom_point(mapping = aes(color = 'red'))
#' p + NoAxes()
#'
NoAxes <- function(..., keep.text = FALSE, keep.ticks = FALSE) {
  blank <- element_blank()
  no.axes.theme <- theme(
    # Remove the axis elements
    axis.line.x = blank,
    axis.line.y = blank,
    # Validate the theme
    validate = TRUE,
    ...
  )
  if (!keep.text) {
    no.axes.theme <- no.axes.theme + theme(
      axis.text.x = blank,
      axis.text.y = blank,
      axis.title.x = blank,
      axis.title.y = blank,
      validate = TRUE,
      ...
    )
  }
  if (!keep.ticks){
    no.axes.theme <- no.axes.theme + theme(
      axis.ticks.x = blank,
      axis.ticks.y = blank,
      validate = TRUE,
      ...
    )
  }
  return(no.axes.theme)
}

#' @inheritParams SeuratTheme
#'
#' @importFrom ggplot2 theme
#' @export
#'
#' @rdname SeuratTheme
#' @aliases NoLegend
#'
#' @examples
#' # Generate a plot with no legend
#' library(ggplot2)
#' df <- data.frame(x = rnorm(n = 100, mean = 20, sd = 2), y = rbinom(n = 100, size = 100, prob = 0.2))
#' p <- ggplot(data = df, mapping = aes(x = x, y = y)) + geom_point(mapping = aes(color = 'red'))
#' p + NoLegend()
#'
NoLegend <- function(...) {
  no.legend.theme <- theme(
    # Remove the legend
    legend.position = 'none',
    # Validate the theme
    validate = TRUE,
    ...
  )
  return(no.legend.theme)
}

#' @inheritParams SeuratTheme
#'
#' @importFrom ggplot2 theme element_blank
#' @export
#'
#' @rdname SeuratTheme
#' @aliases NoGrid
#'
#' @examples
#' # Generate a plot with no grid lines
#' library(ggplot2)
#' df <- data.frame(x = rnorm(n = 100, mean = 20, sd = 2), y = rbinom(n = 100, size = 100, prob = 0.2))
#' p <- ggplot(data = df, mapping = aes(x = x, y = y)) + geom_point(mapping = aes(color = 'red'))
#' p + NoGrid()
#'
NoGrid <- function(...) {
  no.grid.theme <- theme(
    # Set grid lines to blank
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    # Validate the theme
    validate = TRUE,
    ...
  )
  return(no.grid.theme)
}

#' @inheritParams SeuratTheme
#'
#' @importFrom ggplot2 theme element_text
#' @export
#'
#' @rdname SeuratTheme
#' @aliases SeuratAxes
#'
SeuratAxes <- function(...) {
  axes.theme <- theme(
    # Set axis things
    axis.title = element_text(face = 'bold', color = '#990000', size = 16),
    axis.text = element_text(vjust = 0.5, size = 12),
    # Validate the theme
    validate = TRUE,
    ...
  )
  return(axes.theme)
}

#' @inheritParams SeuratTheme
#'
#' @export
#'
#' @rdname SeuratTheme
#' @aliases SpatialTheme
#'
SpatialTheme <- function(...) {
  return(DarkTheme() + NoAxes() + NoGrid() + NoLegend(...))
}

#' @inheritParams SeuratTheme
#' @param position A position to restore the legend to
#'
#' @importFrom ggplot2 theme
#' @export
#'
#' @rdname SeuratTheme
#' @aliases RestoreLegend
#'
RestoreLegend <- function(..., position = 'right') {
  restored.theme <- theme(
    # Restore legend position
    legend.position = 'right',
    # Validate the theme
    validate = TRUE,
    ...
  )
  return(restored.theme)
}

#' @inheritParams SeuratTheme
#'
#' @importFrom ggplot2 theme element_text
#' @export
#'
#' @rdname SeuratTheme
#' @aliases RotatedAxis
#'
RotatedAxis <- function(...) {
  rotated.theme <- theme(
    # Rotate X axis text
    axis.text.x = element_text(angle = 45, hjust = 1),
    # Validate the theme
    validate = TRUE,
    ...
  )
  return(rotated.theme)
}

#' @inheritParams SeuratTheme
#'
#' @importFrom ggplot2 theme element_text
#' @export
#'
#' @rdname SeuratTheme
#' @aliases BoldTitle
#'
BoldTitle <- function(...) {
  bold.theme <- theme(
    # Make the title bold
    plot.title = element_text(size = 20, face = 'bold'),
    # Validate the theme
    validate = TRUE,
    ...
  )
  return(bold.theme)
}

#' @inheritParams SeuratTheme
#'
#' @importFrom ggplot2 theme element_rect
#' @export
#'
#' @rdname SeuratTheme
#' @aliases WhiteBackground
#'
WhiteBackground <- function(...) {
  white.rect = element_rect(fill = 'white')
  white.theme <- theme(
    # Make the plot, panel, and legend key backgrounds white
    plot.background = white.rect,
    panel.background = white.rect,
    legend.key = white.rect,
    # Validate the theme
    validate = TRUE,
    ...
  )
  return(white.theme)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Internal
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Automagically calculate a point size for ggplot2-based scatter plots
#
# It happens to look good
#
# @param data A data frame being passed to ggplot2
#
# @return The "optimal" point size for visualizing these data
#
# @examples
# df <- data.frame(x = rnorm(n = 10000), y = runif(n = 10000))
# AutoPointSize(data = df)
#
AutoPointSize <- function(data) {
     return(min(1583 / nrow(x = data), 1))
}

# Calculate bandwidth for use in ggplot2-based smooth scatter plots
#
# Inspired by MASS::bandwidth.nrd and graphics:::.smoothScatterCalcDensity
#
# @param data A two-column data frame with X and Y coordinates for a plot
#
# @return The calculated bandwidth
#
#' @importFrom stats quantile var
#
Bandwidth <- function(data) {
  r <- diff(x = apply(
    X = data,
    MARGIN = 2,
    FUN = quantile,
    probs = c(0.05, 0.95),
    na.rm = TRUE,
    names = FALSE
  ))
  h <- abs(x = r[2L] - r[1L]) / 1.34
  h <- ifelse(test = h == 0, yes = 1, no = h)
  bandwidth <- 4 * 1.06 *
    min(sqrt(x = apply(X = data, MARGIN = 2, FUN = var)), h) *
    nrow(x = data) ^ (-0.2)
  return(bandwidth)
}

# Blend expression values together
#
# @param data A two-column data frame with expression values for two features
#
# @return A three-column data frame with transformed and blended expression values
#
BlendExpression <- function(data) {
  if (ncol(x = data) != 2) {
    stop("'BlendExpression' only blends two features")
  }
  features <- colnames(x = data)
  data <- as.data.frame(x = apply(
    X = data,
    MARGIN = 2,
    FUN = function(x) {
      return(round(x = 9 * (x - min(x)) / (max(x) - min(x))))
    }
  ))
  data[, 3] <- data[, 1] + data[, 2] * 10
  colnames(x = data) <- c(features, paste(features, collapse = '_'))
  for (i in 1:ncol(x = data)) {
    data[, i] <- factor(x = data[, i])
  }
  return(data)
}

# Create a heatmap of blended colors
#
# @param color.matrix A color matrix of blended colors
#
# @return A ggplot object
#
#' @importFrom grid unit
#' @importFrom cowplot theme_cowplot
#' @importFrom ggplot2 ggplot aes_string scale_fill_manual geom_raster
#' theme scale_y_continuous scale_x_continuous scale_fill_manual
#
# @seealso \code{\link{BlendMatrix}}
#
BlendMap <- function(color.matrix) {
  color.heat <- matrix(
    data = 1:prod(dim(x = color.matrix)) - 1,
    nrow = nrow(x = color.matrix),
    ncol = ncol(x = color.matrix),
    dimnames = list(
      1:nrow(x = color.matrix),
      1:ncol(x = color.matrix)
    )
  )
  xbreaks <- seq.int(from = 0, to = nrow(x = color.matrix), by = 2)
  ybreaks <- seq.int(from = 0, to = ncol(x = color.matrix), by = 2)
  color.heat <- Melt(x = color.heat)
  color.heat$rows <- as.numeric(x = as.character(x = color.heat$rows))
  color.heat$cols <- as.numeric(x = as.character(x = color.heat$cols))
  color.heat$vals <- factor(x = color.heat$vals)
  plot <- ggplot(
    data = color.heat,
    mapping = aes_string(x = 'rows', y = 'cols', fill = 'vals')
  ) +
    geom_raster(show.legend = FALSE) +
    theme(plot.margin = unit(x = rep.int(x = 0, times = 4), units = 'cm')) +
    scale_x_continuous(breaks = xbreaks, expand = c(0, 0), labels = xbreaks) +
    scale_y_continuous(breaks = ybreaks, expand = c(0, 0), labels = ybreaks) +
    scale_fill_manual(values = as.vector(x = color.matrix)) +
    theme_cowplot()
  return(plot)
}

# Create a color matrix of blended colors
#
# @param n Dimensions of blended matrix (n x n)
# @param col.threshold The color cutoff from weak signal to strong signal; ranges from 0 to 1.
# @param two.colors Two colors used for the blend expression.
#
# @return An n x n matrix of blended colors
#
#' @importFrom grDevices rgb colorRamp
#
BlendMatrix <- function(
  n = 10,
  col.threshold = 0.5,
  two.colors = c("#ff0000", "#00ff00"),
  negative.color = "black"
) {
  if (0 > col.threshold || col.threshold > 1) {
    stop("col.threshold must be between 0 and 1")
  }
  C0 <- colorRamp(colors = negative.color)(1)
  ramp <- colorRamp(colors = two.colors)
  C1 <- ramp(x = 0)
  C2 <- ramp(x = 1)
  merge.weight <- min(255 / (C1 + C2 +  C0 + 0.01))
  sigmoid <- function(x) {
    return(1 / (1 + exp(-x)))
  }
  blend_color <- function(
    i,
    j,
    col.threshold,
    n,
    C0,
    C1,
    C2,
    merge.weight
  ) {
    c.min <- sigmoid(5 * (1 / n - col.threshold))
    c.max <- sigmoid(5 * (1 - col.threshold))
    c1_weight <- sigmoid(5 * (i / n - col.threshold))
    c2_weight <- sigmoid(5 * (j / n - col.threshold))
    c0_weight <-  sigmoid(5 * ((i + j) / (2 * n) - col.threshold))
    c1_weight <- (c1_weight - c.min) / (c.max - c.min)
    c2_weight <- (c2_weight - c.min) / (c.max - c.min)
    c0_weight <- (c0_weight - c.min) / (c.max - c.min)
    C1_length <- sqrt(sum((C1 - C0) ** 2))
    C2_length <- sqrt(sum((C2 - C0) ** 2))
    C1_unit <- (C1 - C0) / C1_length
    C2_unit <- (C2 - C0) / C2_length
    C1_weight <- C1_unit * c1_weight
    C2_weight <- C2_unit * c2_weight
    C_blend <- C1_weight * (i - 1) * C1_length / (n - 1) + C2_weight * (j - 1) * C2_length / (n - 1) + (i - 1) * (j - 1) * c0_weight * C0 / (n - 1) ** 2 + C0
    C_blend[C_blend > 255] <- 255
    C_blend[C_blend < 0] <- 0
    return(rgb(
      red = C_blend[, 1],
      green = C_blend[, 2],
      blue = C_blend[, 3],
      alpha = 255,
      maxColorValue = 255
    ))
  }
  blend_matrix <- matrix(nrow = n, ncol = n)
  for (i in 1:n) {
    for (j in 1:n) {
      blend_matrix[i, j] <- blend_color(
        i = i,
        j = j,
        col.threshold = col.threshold,
        n = n,
        C0 = C0,
        C1 = C1,
        C2 = C2,
        merge.weight = merge.weight
      )
    }
  }
  return(blend_matrix)
}

# Convert R colors to hexadecimal
#
# @param ... R colors
#
# @return The hexadecimal representations of input colors
#
#' @importFrom grDevices rgb col2rgb
#
Col2Hex <- function(...) {
  colors <- as.character(x = c(...))
  alpha <- rep.int(x = 255, times = length(x = colors))
  if (sum(sapply(X = colors, FUN = grepl, pattern = '^#')) != 0) {
    hex <- colors[which(x = grepl(pattern = '^#', x = colors))]
    hex.length <- sapply(X = hex, FUN = nchar)
    if (9 %in% hex.length) {
      hex.alpha <- hex[which(x = hex.length == 9)]
      hex.vals <- sapply(X = hex.alpha, FUN = substr, start = 8, stop = 9)
      dec.vals <- sapply(X = hex.vals, FUN = strtoi, base = 16)
      alpha[match(x = hex[which(x = hex.length == 9)], table = colors)] <- dec.vals
    }
  }
  colors <- t(x = col2rgb(col = colors))
  colors <- mapply(
    FUN = function(i, alpha) {
      return(rgb(colors[i, , drop = FALSE], alpha = alpha, maxColorValue = 255))
    },
    i = 1:nrow(x = colors),
    alpha = alpha
  )
  return(colors)
}

# Find the default DimReduc
#
# Searches for DimReducs matching 'umap', 'tsne', or 'pca', case-insensitive, and
# in that order. Priority given to DimReducs matching the DefaultAssay or assay specified
# (eg. 'pca' for the default assay weights higher than 'umap' for a non-default assay)
#
# @param object A Seurat object
# @param assay Name of assay to use; defaults to the default assay of the object
#
# @return The default DimReduc, if possible
#
DefaultDimReduc <- function(object, assay = NULL) {
  assay <- assay %||% DefaultAssay(object = object)
  drs.use <- c('umap', 'tsne', 'pca')
  dim.reducs <- FilterObjects(object = object, classes.keep = 'DimReduc')
  drs.assay <- Filter(
    f = function(x) {
      return(DefaultAssay(object = object[[x]]) == assay)
    },
    x = dim.reducs
  )
  if (length(x = drs.assay) > 0) {
    index <- lapply(
      X = drs.use,
      FUN = grep,
      x = drs.assay,
      ignore.case = TRUE
    )
    index <- Filter(f = length, x = index)
    if (length(x = index) > 0) {
      return(drs.assay[min(index[[1]])])
    }
  }
  index <- lapply(
    X = drs.use,
    FUN = grep,
    x = dim.reducs,
    ignore.case = TRUE
  )
  index <- Filter(f = length, x = index)
  if (length(x = index) < 1) {
    stop(
      "Unable to find a DimReduc matching one of '",
      paste(drs.use[1:(length(x = drs.use) - 1)], collapse = "', '"),
      "', or '",
      drs.use[length(x = drs.use)],
      "', please specify a dimensional reduction to use",
      call. = FALSE
    )
  }
  return(dim.reducs[min(index[[1]])])
}

# Plot feature expression by identity
#
# Basically combines the codebase for VlnPlot and RidgePlot
#
# @param object Seurat object
# @param type Plot type, choose from 'ridge', 'violin', or 'multiViolin'
# @param features Features to plot (gene expression, metrics, PC scores,
# anything that can be retreived by FetchData)
# @param idents Which classes to include in the plot (default is all)
# @param ncol Number of columns if multiple plots are displayed
# @param sort Sort identity classes (on the x-axis) by the average expression of the attribute being potted
# @param y.max Maximum y axis value
# @param same.y.lims Set all the y-axis limits to the same values
# @param adjust Adjust parameter for geom_violin
# @param pt.size Point size for geom_violin
# @param cols Colors to use for plotting
# @param group.by Group (color) cells in different ways (for example, orig.ident)
# @param split.by A variable to split the plot by
# @param log plot Y axis on log scale
# @param slot Use non-normalized counts data for plotting
# @param combine Combine plots into a single \code{\link[patchwork]{patchwork}ed}
# ggplot object. If \code{FALSE}, return a list of ggplot objects
#
# @return A \code{\link[patchwork]{patchwork}ed} ggplot object if
# \code{combine = TRUE}; otherwise, a list of ggplot objects
#
#' @importFrom scales hue_pal
#' @importFrom ggplot2 xlab ylab
#' @importFrom patchwork wrap_plots
#
ExIPlot <- function(
  object,
  features,
  type = 'violin',
  idents = NULL,
  ncol = NULL,
  sort = FALSE,
  assay = NULL,
  y.max = NULL,
  same.y.lims = FALSE,
  adjust = 1,
  cols = NULL,
  pt.size = 0,
  group.by = NULL,
  split.by = NULL,
  log = FALSE,
  slot = 'data',
  combine = TRUE
) {
  assay <- assay %||% DefaultAssay(object = object)
  DefaultAssay(object = object) <- assay
  ncol <- ncol %||% ifelse(
    test = length(x = features) > 9,
    yes = 4,
    no = min(length(x = features), 3)
  )
  data <- FetchData(object = object, vars = features, slot = slot)
  features <- colnames(x = data)
  if (is.null(x = idents)) {
    cells <- colnames(x = object)
  } else {
    cells <- names(x = Idents(object = object)[Idents(object = object) %in% idents])
  }
  data <- data[cells, , drop = FALSE]
  idents <- if (is.null(x = group.by)) {
    Idents(object = object)[cells]
  } else {
    object[[group.by, drop = TRUE]][cells]
  }
  if (!is.factor(x = idents)) {
    idents <- factor(x = idents)
  }
  if (is.null(x = split.by)) {
    split <- NULL
  } else {
    split <- object[[split.by, drop = TRUE]][cells]
    if (!is.factor(x = split)) {
      split <- factor(x = split)
    }
    if (is.null(x = cols)) {
      cols <- hue_pal()(length(x = levels(x = idents)))
      cols <- Interleave(cols, InvertHex(hexadecimal = cols))
    } else if (length(x = cols) == 1 && cols == 'interaction') {
      split <- interaction(idents, split)
      cols <- hue_pal()(length(x = levels(x = idents)))
    } else {
      cols <- Col2Hex(cols)
    }
    if (length(x = cols) < length(x = levels(x = split))) {
      cols <- Interleave(cols, InvertHex(hexadecimal = cols))
    }
    cols <- rep_len(x = cols, length.out = length(x = levels(x = split)))
    names(x = cols) <- sort(x = levels(x = split))
  }
  if (same.y.lims && is.null(x = y.max)) {
    y.max <- max(data)
  }
  plots <- lapply(
    X = features,
    FUN = function(x) {
      return(SingleExIPlot(
        type = type,
        data = data[, x, drop = FALSE],
        idents = idents,
        split = split,
        sort = sort,
        y.max = y.max,
        adjust = adjust,
        cols = cols,
        pt.size = pt.size,
        log = log
      ))
    }
  )
  label.fxn <- switch(
    EXPR = type,
    'violin' = ylab,
    "multiViolin" = ylab,
    'ridge' = xlab,
    stop("Unknown ExIPlot type ", type, call. = FALSE)
  )
  for (i in 1:length(x = plots)) {
    key <- paste0(unlist(x = strsplit(x = features[i], split = '_'))[1], '_')
    obj <- names(x = which(x = Key(object = object) == key))
    if (length(x = obj) == 1) {
      if (inherits(x = object[[obj]], what = 'DimReduc')) {
        plots[[i]] <- plots[[i]] + label.fxn(label = 'Embeddings Value')
      } else if (inherits(x = object[[obj]], what = 'Assay')) {
        next
      } else {
        warning("Unknown object type ", class(x = object), immediate. = TRUE, call. = FALSE)
        plots[[i]] <- plots[[i]] + label.fxn(label = NULL)
      }
    } else if (!features[i] %in% rownames(x = object)) {
      plots[[i]] <- plots[[i]] + label.fxn(label = NULL)
    }
  }
  if (combine) {
    plots <- wrap_plots(plots, ncol = ncol)
    if (length(x = features) > 1) {
      plots <- plots & NoLegend()
    }
  }
  return(plots)
}

# Make a theme for facet plots
#
# @inheritParams SeuratTheme
# @export
#
# @rdname SeuratTheme
# @aliases FacetTheme
#
FacetTheme <- function(...) {
  return(theme(
    strip.background = element_blank(),
    strip.text = element_text(face = 'bold'),
    # Validate the theme
    validate = TRUE,
    ...
  ))
}

# Convert a ggplot2 scatterplot to base R graphics
#
# @param plot A ggplot2 scatterplot
# @param do.plot Create the plot with base R graphics
# @param ... Extra parameters passed to PlotBuild
#
# @return A dataframe with the data that created the ggplot2 scatterplot
#
#' @importFrom ggplot2 ggplot_build
#
GGpointToBase <- function(plot, do.plot = TRUE, ...) {
  plot.build <- ggplot_build(plot = plot)
  cols <- c('x', 'y', 'colour', 'shape', 'size')
  build.use <- which(x = vapply(
    X = plot.build$data,
    FUN = function(dat) {
      return(all(cols %in% colnames(x = dat)))
    },
    FUN.VALUE = logical(length = 1L)
  ))
  if (length(x = build.use) == 0) {
    stop("GGpointToBase only works on geom_point ggplot objects")
  }
  build.data <- plot.build$data[[min(build.use)]]
  plot.data <- build.data[, cols]
  names(x = plot.data) <- c(
    plot.build$plot$labels$x,
    plot.build$plot$labels$y,
    'color',
    'pch',
    'cex'
  )
  if (do.plot) {
    PlotBuild(data = plot.data, ...)
  }
  return(plot.data)
}

# Get X and Y aesthetics from a plot for a certain geom
#
# @param plot A ggplot2 object
# @param geom Geom class to filter to
# @param plot.first Use plot-wide X/Y aesthetics before geom-specific aesthetics
#
# @return A named list with values 'x' for the name of the x aesthetic and 'y' for the y aesthetic
#
GetXYAesthetics <- function(plot, geom = 'GeomPoint', plot.first = TRUE) {
  geoms <- sapply(
    X = plot$layers,
    FUN = function(layer) {
      return(class(x = layer$geom)[1])
    }
  )
  geoms <- which(x = geoms == geom)
  if (length(x = geoms) == 0) {
    stop("Cannot find a geom of class ", geom)
  }
  geoms <- min(geoms)
  if (plot.first) {
    x <- as.character(x = plot$mapping$x %||% plot$layers[[geoms]]$mapping$x)[2]
    y <- as.character(x = plot$mapping$y %||% plot$layers[[geoms]]$mapping$y)[2]
  } else {
    x <- as.character(x = plot$layers[[geoms]]$mapping$x %||% plot$mapping$x)[2]
    y <- as.character(x = plot$layers[[geoms]]$mapping$y %||% plot$mapping$y)[2]
  }
  return(list('x' = x, 'y' = y))
}

# A split violin plot geom
#
#' @importFrom scales zero_range
#' @importFrom ggplot2 GeomPolygon
#' @importFrom grid grobTree grobName
#
# @author jan-glx on StackOverflow
# @references \url{https://stackoverflow.com/questions/35717353/split-violin-plot-with-ggplot2}
# @seealso \code{\link[ggplot2]{geom_violin}}
#
GeomSplitViolin <- ggproto(
  "GeomSplitViolin",
  GeomViolin,
  # setup_data = function(data, params) {
  #   data$width <- data$width %||% params$width %||% (resolution(data$x, FALSE) * 0.9)
  #   data <- plyr::ddply(data, "group", transform, xmin = x - width/2, xmax = x + width/2)
  #   e <- globalenv()
  #   name <- paste(sample(x = letters, size = 5), collapse = '')
  #   message("Saving initial data to ", name)
  #   e[[name]] <- data
  #   return(data)
  # },
  draw_group = function(self, data, ..., draw_quantiles = NULL) {
    data$xminv <- data$x - data$violinwidth * (data$x - data$xmin)
    data$xmaxv <- data$x + data$violinwidth * (data$xmax - data$x)
    grp <- data[1, 'group']
    if (grp %% 2 == 1) {
      data$x <- data$xminv
      data.order <- data$y
    } else {
      data$x <- data$xmaxv
      data.order <- -data$y
    }
    newdata <- data[order(data.order), , drop = FALSE]
    newdata <- rbind(
      newdata[1, ],
      newdata,
      newdata[nrow(x = newdata), ],
      newdata[1, ]
    )
    newdata[c(1, nrow(x = newdata) - 1, nrow(x = newdata)), 'x'] <- round(x = newdata[1, 'x'])
    grob <- if (length(x = draw_quantiles) > 0 & !zero_range(x = range(data$y))) {
      stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <= 1))
      quantiles <- QuantileSegments(data = data, draw.quantiles = draw_quantiles)
      aesthetics <- data[rep.int(x = 1, times = nrow(x = quantiles)), setdiff(x = names(x = data), y = c("x", "y")), drop = FALSE]
      aesthetics$alpha <- rep.int(x = 1, nrow(x = quantiles))
      both <- cbind(quantiles, aesthetics)
      quantile.grob <- GeomPath$draw_panel(both, ...)
      grobTree(GeomPolygon$draw_panel(newdata, ...), name = quantile.grob)
    }
    else {
      GeomPolygon$draw_panel(newdata, ...)
    }
    grob$name <- grobName(grob = grob, prefix = 'geom_split_violin')
    return(grob)
  }
)

# Create a split violin plot geom
#
# @inheritParams ggplot2::geom_violin
#
#' @importFrom ggplot2 layer
#
# @author jan-glx on StackOverflow
# @references \url{https://stackoverflow.com/questions/35717353/split-violin-plot-with-ggplot2}
# @seealso \code{\link[ggplot2]{geom_violin}}
#
geom_split_violin <- function(
  mapping = NULL,
  data = NULL,
  stat = 'ydensity',
  position = 'identity',
  ...,
  draw_quantiles = NULL,
  trim = TRUE,
  scale = 'area',
  na.rm = FALSE,
  show.legend = NA,
  inherit.aes = TRUE
) {
  return(layer(
    data = data,
    mapping = mapping,
    stat = stat,
    geom = GeomSplitViolin,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      trim = trim,
      scale = scale,
      draw_quantiles = draw_quantiles,
      na.rm = na.rm,
      ...
    )
  ))
}

# Invert a Hexadecimal color
#
# @param hexadecimal A character vector of hexadecimal colors
#
# @return Hexadecimal representations of the inverted color
#
# @author Matt Lagrandeur
# @references \url{http://www.mattlag.com/scripting/hexcolorinverter.php}
#
InvertHex <- function(hexadecimal) {
  return(vapply(
    X = toupper(x = hexadecimal),
    FUN = function(hex) {
      hex <- unlist(x = strsplit(
        x = gsub(pattern = '#', replacement = '', x = hex),
        split = ''
      ))
      key <- toupper(x = as.hexmode(x = 15:0))
      if (!all(hex %in% key)) {
        stop('All hexadecimal colors must be valid hexidecimal numbers from 0-9 and A-F')
      }
      if (length(x = hex) == 8) {
        alpha <- hex[7:8]
        hex <- hex[1:6]
      } else if (length(x = hex) == 6) {
        alpha <- NULL
      } else {
        stop("All hexidecimal colors must be either 6 or 8 characters in length, excluding the '#'")
      }
      value <- rev(x = key)
      inv.hex <- vapply(
        X = hex,
        FUN = function(x) {
          return(value[grep(pattern = x, x = key)])
        },
        FUN.VALUE = character(length = 1L)
      )
      inv.hex <- paste(inv.hex, collapse = '')
      return(paste0('#', inv.hex, paste(alpha, collapse = '')))
    },
    FUN.VALUE = character(length = 1L),
    USE.NAMES = FALSE
  ))
}

# Make label information for ggplot2-based scatter plots
#
# @param data A three-column data frame (accessed with \code{plot$data})
# The first column should be the X axis, the second the Y, and the third should be grouping information
#
# @return A dataframe with three columns: centers along the X axis, centers along the Y axis, and group information
#
#' @importFrom stats median
#
MakeLabels <- function(data) {
  groups <- as.character(x = na.omit(object = unique(x = data[, 3])))
  labels <- lapply(
    X = groups,
    FUN = function(group) {
      data.use <- data[data[, 3] == group, 1:2]
      return(apply(X = data.use, MARGIN = 2, FUN = median, na.rm = TRUE))
    }
  )
  names(x = labels) <- groups
  labels <- as.data.frame(x = t(x = as.data.frame(x = labels)))
  labels[, colnames(x = data)[3]] <- groups
  return(labels)
}

# Create a scatterplot with data from a ggplot2 scatterplot
#
# @param plot.data The original ggplot2 scatterplot data
# This is taken from ggplot2::ggplot_build
# @param dark.theme Plot using a dark theme
# @param smooth Use a smooth scatterplot instead of a standard scatterplot
# @param ... Extra parameters passed to graphics::plot or graphics::smoothScatter
#
#' @importFrom graphics axis plot smoothScatter
#
PlotBuild <- function(data, dark.theme = FALSE, smooth = FALSE, ...) {
  #   Do we use a smooth scatterplot?
  #   Take advantage of functions as first class objects
  #   to dynamically choose normal vs smooth scatterplot
  myplot <- ifelse(test = smooth, yes = smoothScatter, no = plot)
  CheckDots(..., fxns = myplot)
  if (dark.theme) {
    par(bg = 'black')
    axes = FALSE
    col.lab = 'white'
  } else {
    axes = 'TRUE'
    col.lab = 'black'
  }
  myplot(
    data[, c(1, 2)],
    col = data$color,
    pch = data$pch,
    cex = vapply(
      X = data$cex,
      FUN = function(x) {
        return(max(x / 2, 0.5))
      },
      FUN.VALUE = numeric(1)
    ),
    axes = axes,
    col.lab = col.lab,
    col.main = col.lab,
    ...
  )
  if (dark.theme) {
    axis(
      side = 1,
      at = NULL,
      labels = TRUE,
      col.axis = col.lab,
      col = col.lab
    )
    axis(
      side = 2,
      at = NULL,
      labels = TRUE,
      col.axis = col.lab,
      col = col.lab
    )
  }
}

# Locate points on a plot and return them
#
# @param plot A ggplot2 plot
# @param recolor Do we recolor the plot to highlight selected points?
# @param dark.theme Plot using a dark theme
# @param ... Exptra parameters to PlotBuild
#
# @return A dataframe of x and y coordinates for points selected
#
#' @importFrom graphics locator
# @importFrom SDMTools pnt.in.poly
#
PointLocator <- function(plot, recolor = TRUE, dark.theme = FALSE, ...) {
  #   Convert the ggplot object to a data.frame
  PackageCheck('SDMTools')
  plot.data <- GGpointToBase(plot = plot, dark.theme = dark.theme, ...)
  npoints <- nrow(x = plot.data)
  cat("Click around the cluster of points you wish to select\n")
  cat("ie. select the vertecies of a shape around the cluster you\n")
  cat("are interested in. Press <Esc> when finished (right click for R-terminal users)\n\n")
  polygon <- locator(n = npoints, type = 'l')
  polygon <- data.frame(polygon)
  #   pnt.in.poly returns a data.frame of points
  points.all <- SDMTools::pnt.in.poly(
    pnts = plot.data[, c(1, 2)],
    poly.pnts = polygon
  )
  #   Find the located points
  points.located <- points.all[which(x = points.all$pip == 1), ]
  #   If we're recoloring, do the recolor
  if (recolor) {
    no <- ifelse(test = dark.theme, yes = 'white', no = '#C3C3C3')
    points.all$color <- ifelse(test = points.all$pip == 1, yes = '#DE2D26', no = no)
    plot.data$color <- points.all$color
    PlotBuild(data = plot.data, dark.theme = dark.theme, ...)
  }
  return(points.located[, c(1, 2)])
}

# Create quantile segments for quantiles on violin plots in ggplot2
#
# @param data Data being plotted
# @param draw.quantiles Quantiles to draw
#
#' @importFrom stats approxfun
#
# @author Hadley Wickham (I presume)
# @seealso \code{\link[ggplot2]{geom_violin}}
#
QuantileSegments <- function(data, draw.quantiles) {
  densities <- cumsum(x = data$density) / sum(data$density)
  ecdf <- approxfun(x = densities, y = data$y)
  ys <- ecdf(v = draw.quantiles)
  violin.xminvs <- approxfun(x = data$y, y = data$xminv)(v = ys)
  violin.xmaxvs <- approxfun(x = data$y, y = data$xmaxv)(v = ys)
  return(data.frame(
    x = as.vector(x = t(x = data.frame(violin.xminvs, violin.xmaxvs))),
    y = rep(x = ys, each = 2),
    group = rep(x = ys, each = 2)
  ))
}

# Scale vector to min and max cutoff values
#
# @param vec a vector
# @param cutoffs A two-length vector of cutoffs to be passed to \code{\link{SetQuantile}}
#
# @return Returns a vector
#
ScaleColumn <- function(vec, cutoffs) {
  if (!length(x = cutoffs) == 2) {
    stop("Two cutoffs (a low and high) are needed")
  }
  cutoffs <- sapply(
    X = cutoffs,
    FUN = SetQuantile,
    data = vec
  )
  vec[vec < min(cutoffs)] <- min(cutoffs)
  vec[vec > max(cutoffs)] <- max(cutoffs)
  return(vec)
}

# Set highlight information
#
# @param cells.highlight Cells to highlight
# @param cells.all A character vector of all cell names
# @param sizes.highlight Sizes of cells to highlight
# @param cols.highlight Colors to highlight cells as
# @param col.base Base color to use for unselected cells
# @param pt.size Size of unselected cells
#
# @return A list will cell highlight information
# \describe{
#   \item{plot.order}{An order to plot cells in}
#   \item{highlight}{A vector giving group information for each cell}
#   \item{size}{A vector giving size information for each cell}
#   \item{color}{Colors for highlighting in the order of plot.order}
# }
#
SetHighlight <- function(
  cells.highlight,
  cells.all,
  sizes.highlight,
  cols.highlight,
  col.base = 'black',
  pt.size = 1
) {
  if (is.character(x = cells.highlight)) {
    cells.highlight <- list(cells.highlight)
  } else if (is.data.frame(x = cells.highlight) || !is.list(x = cells.highlight)) {
    cells.highlight <- as.list(x = cells.highlight)
  }
  cells.highlight <- lapply(
    X = cells.highlight,
    FUN = function(cells) {
      cells.return <- if (is.character(x = cells)) {
        cells[cells %in% cells.all]
      } else {
        cells <- as.numeric(x = cells)
        cells <- cells[cells <= length(x = cells.all)]
        cells.all[cells]
      }
      return(cells.return)
    }
  )
  cells.highlight <- Filter(f = length, x = cells.highlight)
  names.highlight <- if (is.null(x = names(x = cells.highlight))) {
    paste0('Group_', 1L:length(x = cells.highlight))
  } else {
    names(x = cells.highlight)
  }
  sizes.highlight <- rep_len(
    x = sizes.highlight,
    length.out = length(x = cells.highlight)
  )
  cols.highlight <- c(
    col.base,
    rep_len(x = cols.highlight, length.out = length(x = cells.highlight))
  )
  size <- rep_len(x = pt.size, length.out = length(x = cells.all))
  highlight <- rep_len(x = NA_character_, length.out = length(x = cells.all))
  if (length(x = cells.highlight) > 0) {
    for (i in 1:length(x = cells.highlight)) {
      cells.check <- cells.highlight[[i]]
      index.check <- match(x = cells.check, cells.all)
      highlight[index.check] <- names.highlight[i]
      size[index.check] <- sizes.highlight[i]
    }
  }
  plot.order <- sort(x = unique(x = highlight), na.last = TRUE)
  plot.order[is.na(x = plot.order)] <- 'Unselected'
  highlight[is.na(x = highlight)] <- 'Unselected'
  highlight <- as.factor(x = highlight)
  return(list(
    plot.order = plot.order,
    highlight = highlight,
    size = size,
    color = cols.highlight
  ))
}

# Find the quantile of a data
#
# Converts a quantile in character form to a number regarding some data
# String form for a quantile is represented as a number prefixed with 'q'
# For example, 10th quantile is 'q10' while 2nd quantile is 'q2'
#
# Will only take a quantile of non-zero data values
#
# @param cutoff The cutoff to turn into a quantile
# @param data The data to turn find the quantile of
#
# @return The numerical representation of the quantile
#
#' @importFrom stats quantile
#
SetQuantile <- function(cutoff, data) {
  if (grepl(pattern = '^q[0-9]{1,2}$', x = as.character(x = cutoff), perl = TRUE)) {
    this.quantile <- as.numeric(x = sub(
      pattern = 'q',
      replacement = '',
      x = as.character(x = cutoff)
    )) / 100
    data <- unlist(x = data)
    data <- data[data > 0]
    cutoff <- quantile(x = data, probs = this.quantile)
  }
  return(as.numeric(x = cutoff))
}

globalVariables(names = '..density..', package = 'Seurat')
# A single correlation plot
#
# @param data.plot A data frame with two columns to be plotted
# @param col.by A vector or factor of values to color the plot by
# @param cols An optional vector of colors to use
# @param pt.size Point size for the plot
# @param smooth Make a smoothed scatter plot
# @param rows.highight A vector of rows to highlight (like cells.highlight in SingleDimPlot)
# @param legend.title Optional legend title
# @param ... Extra parameters to MASS::kde2d
#
#' @importFrom stats cor
# #' @importFrom MASS kde2d
#' @importFrom cowplot theme_cowplot
#' @importFrom RColorBrewer brewer.pal.info
#' @importFrom ggplot2 ggplot geom_point aes_string labs scale_color_brewer
#' scale_color_manual guides stat_density2d aes scale_fill_continuous
#'
SingleCorPlot <- function(
  data,
  col.by = NULL,
  cols = NULL,
  pt.size = NULL,
  smooth = FALSE,
  rows.highlight = NULL,
  legend.title = NULL,
  na.value = 'grey50',
  span = NULL
) {
  pt.size <- pt.size <- pt.size %||% AutoPointSize(data = data)
  orig.names <- colnames(x = data)
  names.plot <- colnames(x = data) <- gsub(
    pattern = '-',
    replacement = '.',
    x = colnames(x = data),
    fixed = TRUE
  )
  names.plot <- colnames(x = data) <- gsub(
    pattern = ':',
    replacement = '.',
    x = colnames(x = data),
    fixed = TRUE
  )
  if (ncol(x = data) < 2) {
    msg <- "Too few variables passed"
    if (ncol(x = data) == 1) {
      msg <- paste0(msg, ', only have ', colnames(x = data)[1])
    }
    stop(msg, call. = FALSE)
  }
  plot.cor <- round(x = cor(x = data[, 1], y = data[, 2]), digits = 2)
  if (!is.null(x = rows.highlight)) {
    highlight.info <- SetHighlight(
      cells.highlight = rows.highlight,
      cells.all = rownames(x = data),
      sizes.highlight = pt.size,
      cols.highlight = 'red',
      col.base = 'black',
      pt.size = pt.size
    )
    cols <- highlight.info$color
    col.by <- factor(
      x = highlight.info$highlight,
      levels = rev(x = highlight.info$plot.order)
    )
    plot.order <- order(col.by)
    data <- data[plot.order, ]
    col.by <- col.by[plot.order]
  }
  if (!is.null(x = col.by)) {
    data$colors <- col.by
  }
  plot <- ggplot(
    data = data,
    mapping = aes_string(x = names.plot[1], y = names.plot[2])
  ) +
    labs(
      x = orig.names[1],
      y = orig.names[2],
      title = plot.cor,
      color = legend.title
    )
  if (smooth) {
    # density <- kde2d(x = data[, names.plot[1]], y = data[, names.plot[2]], h = Bandwidth(data = data[, names.plot]), n = 200)
    # density <- data.frame(
    #   expand.grid(
    #     x = density$x,
    #     y = density$y
    #   ),
    #   density = as.vector(x = density$z)
    # )
    plot <- plot + stat_density2d(
      mapping = aes(fill = ..density.. ^ 0.25),
        geom = 'tile',
        contour = FALSE,
        n = 200,
        h = Bandwidth(data = data[, names.plot])
      ) +
      # geom_tile(
      #   mapping = aes_string(
      #     x = 'x',
      #     y = 'y',
      #     fill = 'density'
      #   ),
      #   data = density
      # ) +
      scale_fill_continuous(low = 'white', high = 'dodgerblue4') +
      guides(fill = FALSE)
  }
  if (!is.null(x = col.by)) {
    plot <- plot + geom_point(
      mapping = aes_string(color = 'colors'),
      position = 'jitter',
      size = pt.size
    )
  } else {
    plot <- plot + geom_point(position = 'jitter', size = pt.size)
  }
  if (!is.null(x = cols)) {
    cols.scale <- if (length(x = cols) == 1 && cols %in% rownames(x = brewer.pal.info)) {
      scale_color_brewer(palette = cols)
    } else {
      scale_color_manual(values = cols, na.value = na.value)
    }
    plot <- plot + cols.scale
    if (!is.null(x = rows.highlight)) {
      plot <- plot + guides(color = FALSE)
    }
  }
  plot <- plot + theme_cowplot() + theme(plot.title = element_text(hjust = 0.5))
  if (!is.null(x = span)) {
    plot <- plot + geom_smooth(
      mapping = aes_string(x = names.plot[1], y = names.plot[2]),
      method = 'loess',
      span = span
    )
  }
  return(plot)
}

# Plot a single dimension
#
# @param data Data to plot
# @param dims A two-length numeric vector with dimensions to use
# @param pt.size Adjust point size for plotting
# @param col.by ...
# @param cols Vector of colors, each color corresponds to an identity class. This may also be a single character
# or numeric value corresponding to a palette as specified by \code{\link[RColorBrewer]{brewer.pal.info}}.
# By default, ggplot2 assigns colors
# @param shape.by If NULL, all points are circles (default). You can specify any cell attribute
# (that can be pulled with FetchData) allowing for both different colors and different shapes on
# cells.
# @param order Specify the order of plotting for the idents. This can be useful for crowded plots if
# points of interest are being buried. Provide either a full list of valid idents or a subset to be
# plotted last (on top).
# @param label Whether to label the clusters
# @param repel Repel labels
# @param label.size Sets size of labels
# @param cells.highlight A list of character or numeric vectors of cells to
# highlight. If only one group of cells desired, can simply
# pass a vector instead of a list. If set, colors selected cells to the color(s)
# in \code{cols.highlight} and other cells black (white if dark.theme = TRUE);
#  will also resize to the size(s) passed to \code{sizes.highlight}
# @param cols.highlight A vector of colors to highlight the cells as; will
# repeat to the length groups in cells.highlight
# @param sizes.highlight Size of highlighted cells; will repeat to the length
# groups in cells.highlight
# @param na.value Color value for NA points when using custom scale.
#
#' @importFrom cowplot theme_cowplot
#' @importFrom RColorBrewer brewer.pal.info
#' @importFrom ggplot2 ggplot aes_string labs geom_text guides
#'  scale_color_brewer scale_color_manual element_rect guide_legend discrete_scale
#'
SingleDimPlot <- function(
  data,
  dims,
  col.by = NULL,
  cols = NULL,
  pt.size = NULL,
  shape.by = NULL,
  order = NULL,
  label = FALSE,
  repel = FALSE,
  label.size = 4,
  cells.highlight = NULL,
  cols.highlight = '#DE2D26',
  sizes.highlight = 1,
  na.value = 'grey50'
) {
  pt.size <- pt.size %||% AutoPointSize(data = data)
  if (length(x = dims) != 2) {
    stop("'dims' must be a two-length vector")
  }
  if (!is.data.frame(x = data)) {
    data <- as.data.frame(x = data)
  }
  if (is.character(x = dims) && !all(dims %in% colnames(x = data))) {
    stop("Cannot find dimensions to plot in data")
  } else if (is.numeric(x = dims)) {
    dims <- colnames(x = data)[dims]
  }
  if (!is.null(x = cells.highlight)) {
    highlight.info <- SetHighlight(
      cells.highlight = cells.highlight,
      cells.all = rownames(x = data),
      sizes.highlight = sizes.highlight %||% pt.size,
      cols.highlight = cols.highlight,
      col.base = cols[1] %||% '#C3C3C3',
      pt.size = pt.size
    )
    order <- highlight.info$plot.order
    data$highlight <- highlight.info$highlight
    col.by <- 'highlight'
    pt.size <- highlight.info$size
    cols <- highlight.info$color
  }
  if (!is.null(x = order) && !is.null(x = col.by)) {
    if (typeof(x = order) == "logical") {
      if (order) {
        data <- data[order(data[, col.by]), ]
      }
    } else {
      order <- rev(x = c(
        order,
        setdiff(x = unique(x = data[, col.by]), y = order)
      ))
      data[, col.by] <- factor(x = data[, col.by], levels = order)
      new.order <- order(x = data[, col.by])
      data <- data[new.order, ]
      if (length(x = pt.size) == length(x = new.order)) {
        pt.size <- pt.size[new.order]
      }
    }
  }
  if (!is.null(x = col.by) && !col.by %in% colnames(x = data)) {
    warning("Cannot find ", col.by, " in plotting data, not coloring plot")
    col.by <- NULL
  } else {
    # col.index <- grep(pattern = col.by, x = colnames(x = data), fixed = TRUE)
    col.index <- match(x = col.by, table = colnames(x = data))
    if (grepl(pattern = '^\\d', x = col.by)) {
      # Do something for numbers
      col.by <- paste0('x', col.by)
    } else if (grepl(pattern = '-', x = col.by)) {
      # Do something for dashes
      col.by <- gsub(pattern = '-', replacement = '.', x = col.by)
    }
    colnames(x = data)[col.index] <- col.by
  }
  if (!is.null(x = shape.by) && !shape.by %in% colnames(x = data)) {
    warning("Cannot find ", shape.by, " in plotting data, not shaping plot")
  }
  plot <- ggplot(data = data) +
    geom_point(
      mapping = aes_string(
        x = dims[1],
        y = dims[2],
        color = paste0("`", col.by, "`"),
        shape = shape.by
      ),
      size = pt.size
    ) +
    guides(color = guide_legend(override.aes = list(size = 3))) +
    labs(color = NULL)
  if (label && !is.null(x = col.by)) {
    plot <- LabelClusters(
      plot = plot,
      id = col.by,
      repel = repel,
      size = label.size
    )
  }
  if (!is.null(x = cols)) {
    if (length(x = cols) == 1 && (is.numeric(x = cols) || cols %in% rownames(x = brewer.pal.info))) {
      scale <- scale_color_brewer(palette = cols, na.value = na.value)
    } else if (length(x = cols) == 1 && (cols %in% c('alphabet', 'alphabet2', 'glasbey', 'polychrome', 'stepped'))) {
      colors <- DiscretePalette(length(unique(data[[col.by]])), palette = cols)
      scale <- scale_color_manual(values = colors, na.value = na.value)
    } else {
      scale <- scale_color_manual(values = cols, na.value = na.value)
    }
    plot <- plot + scale
  }
  plot <- plot + theme_cowplot()
  return(plot)
}

# Plot a single expression by identity on a plot
#
# @param type Make either a 'ridge' or 'violin' plot
# @param data Data to plot
# @param idents Idents to use
# @param sort Sort identity classes (on the x-axis) by the average
# expression of the attribute being potted
# @param y.max Maximum Y value to plot
# @param adjust Adjust parameter for geom_violin
# @param cols Colors to use for plotting
# @param log plot Y axis on log scale
# @param seed.use Random seed to use. If NULL, don't set a seed
#
# @return A ggplot-based Expression-by-Identity plot
#
# @import ggplot2
#' @importFrom stats rnorm
#' @importFrom utils globalVariables
#' @importFrom ggridges geom_density_ridges theme_ridges
#' @importFrom ggplot2 ggplot aes_string theme labs geom_violin geom_jitter ylim
#' scale_fill_manual scale_y_log10 scale_x_log10 scale_y_discrete scale_x_continuous waiver
#' @importFrom cowplot theme_cowplot
#'
SingleExIPlot <- function(
  data,
  idents,
  split = NULL,
  type = 'violin',
  sort = FALSE,
  y.max = NULL,
  adjust = 1,
  pt.size = 0,
  cols = NULL,
  seed.use = 42,
  log = FALSE
) {
  if (!is.null(x = seed.use)) {
    set.seed(seed = seed.use)
  }
  if (!is.data.frame(x = data) || ncol(x = data) != 1) {
    stop("'SingleExIPlot requires a data frame with 1 column")
  }
  feature <- colnames(x = data)
  data$ident <- idents
  if ((is.character(x = sort) && nchar(x = sort) > 0) || sort) {
    data$ident <- factor(
      x = data$ident,
      levels = names(x = rev(x = sort(
        x = tapply(
          X = data[, feature],
          INDEX = data$ident,
          FUN = mean
        ),
        decreasing = grepl(pattern = paste0('^', tolower(x = sort)), x = 'decreasing')
      )))
    )
  }
  if (log) {
    noise <- rnorm(n = length(x = data[, feature])) / 200
    data[, feature] <- data[, feature] + 1
  } else {
    noise <- rnorm(n = length(x = data[, feature])) / 100000
  }
  if (all(data[, feature] == data[, feature][1])) {
    warning(paste0("All cells have the same value of ", feature, "."))
  } else{
    data[, feature] <- data[, feature] + noise
  }
  axis.label <- 'Expression Level'
  y.max <- y.max %||% max(data[, feature])
  if (type == 'violin' && !is.null(x = split)) {
    data$split <- split
    vln.geom <- geom_split_violin
    fill <- 'split'
  } else if (type == 'multiViolin' && !is.null(x = split )) {
    data$split <- split
    vln.geom <- geom_violin
    fill <- 'split'
    type <- 'violin'
  } else {
    vln.geom <- geom_violin
    fill <- 'ident'
  }
  switch(
    EXPR = type,
    'violin' = {
      x <- 'ident'
      y <- paste0("`", feature, "`")
      xlab <- 'Identity'
      ylab <- axis.label
      geom <- list(
        vln.geom(scale = 'width', adjust = adjust, trim = TRUE),
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
      )
      jitter <- geom_jitter(height = 0, size = pt.size)
      log.scale <- scale_y_log10()
      axis.scale <- ylim
    },
    'ridge' = {
      x <- paste0("`", feature, "`")
      y <- 'ident'
      xlab <- axis.label
      ylab <- 'Identity'
      geom <- list(
        geom_density_ridges(scale = 4),
        theme_ridges(),
        scale_y_discrete(expand = c(0.01, 0)),
        scale_x_continuous(expand = c(0, 0))
      )
      jitter <- geom_jitter(width = 0, size = pt.size)
      log.scale <- scale_x_log10()
      axis.scale <- function(...) {
        invisible(x = NULL)
      }
    },
    stop("Unknown plot type: ", type)
  )
  plot <- ggplot(
    data = data,
    mapping = aes_string(x = x, y = y, fill = fill)[c(2, 3, 1)]
  ) +
    labs(x = xlab, y = ylab, title = feature, fill = NULL) +
    theme_cowplot() +
    theme(plot.title = element_text(hjust = 0.5))
  plot <- do.call(what = '+', args = list(plot, geom))
  plot <- plot + if (log) {
    log.scale
  } else {
    axis.scale(min(data[, feature]), y.max)
  }
  if (pt.size > 0) {
    plot <- plot + jitter
  }
  if (!is.null(x = cols)) {
    if (!is.null(x = split)) {
      idents <- unique(x = as.vector(x = data$ident))
      splits <- unique(x = as.vector(x = data$split))
      labels <- if (length(x = splits) == 2) {
        splits
      } else {
        unlist(x = lapply(
          X = idents,
          FUN = function(pattern, x) {
            x.mod <- gsub(
              pattern = paste0(pattern, '.'),
              replacement = paste0(pattern, ': '),
              x = x,
              fixed = TRUE
            )
            x.keep <- grep(pattern = ': ', x = x.mod, fixed = TRUE)
            x.return <- x.mod[x.keep]
            names(x = x.return) <- x[x.keep]
            return(x.return)
          },
          x = unique(x = as.vector(x = data$split))
        ))
      }
      if (is.null(x = names(x = labels))) {
        names(x = labels) <- labels
      }
    } else {
      labels <- levels(x = droplevels(data$ident))
    }
    plot <- plot + scale_fill_manual(values = cols, labels = labels)
  }
  return(plot)
}

# A single heatmap from base R using image
#
# @param data matrix of data to plot
# @param order optional vector of cell names to specify order in plot
# @param title Title for plot
#
#' @importFrom graphics par plot.new
#
SingleImageMap <- function(data, order = NULL, title = NULL) {
  if (!is.null(x = order)) {
    data <- data[order, ]
  }
  par(mar = c(1, 1, 3, 3))
  plot.new()
  image(
    x = as.matrix(x = data),
    axes = FALSE,
    add = TRUE,
    col = PurpleAndYellow()
  )
  axis(
    side = 4,
    at = seq(from = 0, to = 1, length = ncol(x = data)),
    labels = colnames(x = data),
    las = 1,
    tick = FALSE,
    mgp = c(0, -0.5, 0),
    cex.axis = 0.75
  )
  title(main = title)
}

# A single polygon plot
#
# @param data Data to plot
# @param group.by Grouping variable
# @param ... Extra parameters passed to \code{\link[cowplot]{theme_cowplot}}
#
# @return A ggplot-based plot
#
#' @importFrom cowplot theme_cowplot
#' @importFrom ggplot2 ggplot aes_string geom_polygon
#
# @seealso \code{\link[cowplot]{theme_cowplot}}
#
SinglePolyPlot <- function(data, group.by, ...) {
  plot <- ggplot(data = data, mapping = aes_string(x = 'x', y = 'y')) +
    geom_polygon(mapping = aes_string(fill = group.by, group = 'cell')) +
    coord_fixed() +
    theme_cowplot(...)
  return(plot)
}

# A single heatmap from ggplot2 using geom_raster
#
# @param data A matrix or data frame with data to plot
# @param raster switch between geom_raster and geom_tile
# @param cell.order ...
# @param feature.order ...
# @param cols A vector of colors to use
# @param disp.min Minimum display value (all values below are clipped)
# @param disp.max Maximum display value (all values above are clipped)
# @param limits A two-length numeric vector with the limits for colors on the plot
# @param group.by A vector to group cells by, should be one grouping identity per cell
#
#' @importFrom ggplot2 ggplot aes_string geom_raster scale_fill_gradient
#' scale_fill_gradientn theme element_blank labs geom_point guides guide_legend geom_tile
#
SingleRasterMap <- function(
  data,
  raster = TRUE,
  cell.order = NULL,
  feature.order = NULL,
  colors = PurpleAndYellow(),
  disp.min = -2.5,
  disp.max = 2.5,
  limits = NULL,
  group.by = NULL
) {
  data <- MinMax(data = data, min = disp.min, max = disp.max)
  data <- Melt(x = t(x = data))
  colnames(x = data) <- c('Feature', 'Cell', 'Expression')
  if (!is.null(x = feature.order)) {
    data$Feature <- factor(x = data$Feature, levels = unique(x = feature.order))
  }
  if (!is.null(x = cell.order)) {
    data$Cell <- factor(x = data$Cell, levels = unique(x = cell.order))
  }
  if (!is.null(x = group.by)) {
    data$Identity <- group.by[data$Cell]
  }
  limits <- limits %||% c(min(data$Expression), max(data$Expression))
  if (length(x = limits) != 2 || !is.numeric(x = limits)) {
    stop("limits' must be a two-length numeric vector")
  }
  my_geom <- ifelse(test = raster, yes = geom_raster, no = geom_tile)
  plot <- ggplot(data = data) +
    my_geom(mapping = aes_string(x = 'Cell', y = 'Feature', fill = 'Expression')) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
    scale_fill_gradientn(limits = limits, colors = colors, na.value = "white") +
    labs(x = NULL, y = NULL, fill = group.by %iff% 'Expression') +
    WhiteBackground() + NoAxes(keep.text = TRUE)
  if (!is.null(x = group.by)) {
    plot <- plot + geom_point(
      mapping = aes_string(x = 'Cell', y = 'Feature', color = 'Identity'),
      alpha = 0
    ) +
      guides(color = guide_legend(override.aes = list(alpha = 1)))
  }
  return(plot)
}
