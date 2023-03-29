#' @importFrom utils globalVariables
#' @importFrom ggplot2 fortify GeomViolin ggproto
#' @importFrom SeuratObject DefaultDimReduc
#'
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Generics
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @importFrom methods setGeneric
#'
setGeneric(
  name = '.PrepImageData',
  def = function(data, cells, ...) {
    standardGeneric(f = '.PrepImageData')
  }
)

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
#' @param reduction Which dimensional reduction to use
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
#' @concept visualization
#'
#' @seealso \code{\link[graphics]{image}} \code{\link[ggplot2]{geom_raster}}
#'
#' @examples
#' data("pbmc_small")
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
#' @param vjust Vertical justification of text above color bar
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
#' ggplot_build aes_string geom_text
#' @importFrom patchwork wrap_plots
#' @export
#' @concept visualization
#'
#' @examples
#' data("pbmc_small")
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
  vjust = 0,
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
      if (!is.null(x = names(x = group.colors))) {
        cols <- unname(obj = group.colors[levels(x = group.use)])
      } else {
        cols <- group.colors[1:length(x = levels(x = group.use))] %||% default.colors
      }
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
      x.min <- min(pbuild$layout$panel_params[[1]]$x.range) + 0.1
      x.max <- max(pbuild$layout$panel_params[[1]]$x.range) - 0.1
      plot <- plot +
        annotation_raster(
          raster = t(x = cols[group.use2]),
          xmin = x.min,
          xmax = x.max,
          ymin = y.pos,
          ymax = y.max
        ) +
        coord_cartesian(ylim = c(0, y.max), clip = 'off') +
        scale_color_manual(
          values = cols[-length(x = cols)],
          name = "Identity",
          na.translate = FALSE
        )
      if (label) {
        x.max <- max(pbuild$layout$panel_params[[1]]$x.range)
        # Attempt to pull xdivs from x.major in ggplot2 < 3.3.0; if NULL, pull from the >= 3.3.0 slot
        x.divs <- pbuild$layout$panel_params[[1]]$x.major %||% attr(x = pbuild$layout$panel_params[[1]]$x$get_breaks(), which = "pos")
        x <- data.frame(group = sort(x = group.use), x = x.divs)
        label.x.pos <- tapply(X = x$x, INDEX = x$group, FUN = function(y) {
          if (isTRUE(x = draw.lines)) {
            mean(x = y[-length(x = y)])
          } else {
            mean(x = y)
          }
        })
        label.x.pos <- data.frame(group = names(x = label.x.pos), label.x.pos)
        plot <- plot + geom_text(
          stat = "identity",
          data = label.x.pos,
          aes_string(label = 'group', x = 'label.x.pos'),
          y = y.max + y.max * 0.03 * 0.5 + vjust,
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
#' @concept visualization
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
#' @param slot Slot to pull expression data from (e.g. "counts" or "data")
#' @param stack Horizontally stack plots for each feature
#' @param combine Combine plots into a single \code{\link[patchwork]{patchwork}ed}
#' ggplot object. If \code{FALSE}, return a list of ggplot
#' @param fill.by Color violins/ridges based on either 'feature' or 'ident'
#'
#' @return A \code{\link[patchwork]{patchwork}ed} ggplot object if
#' \code{combine = TRUE}; otherwise, a list of ggplot objects
#'
#' @export
#' @concept visualization
#'
#' @examples
#' data("pbmc_small")
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
  slot = deprecated(),
  layer = 'data',
  stack = FALSE,
  combine = TRUE,
  fill.by = 'feature'
) {
  if (is_present(arg = slot)) {
    deprecate_soft(
      when = '5.0.0',
      what = 'RidgePlot(slot = )',
      with = 'RidgePlot(layer = )'
    )
    layer <- slot %||% layer
  }
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
    layer = layer,
    stack = stack,
    combine = combine,
    fill.by = fill.by
  ))
}

#' Single cell violin plot
#'
#' Draws a violin plot of single cell data (gene expression, metrics, PC
#' scores, etc.)
#'
#' @inheritParams RidgePlot
#' @param pt.size Point size for points
#' @param alpha Alpha value for points
#' @param split.by A variable to split the violin plots by,
#' @param split.plot  plot each group of the split violin plots by multiple or
#' single violin shapes.
#' @param adjust Adjust parameter for geom_violin
#' @param flip flip plot orientation (identities on x-axis)
#' @param add.noise determine if adding a small noise for plotting
#' @param raster Convert points to raster format. Requires 'ggrastr' to be installed.
# default is \code{NULL} which automatically rasterizes if ggrastr is installed and
# number of points exceed 100,000.
#'
#' @return A \code{\link[patchwork]{patchwork}ed} ggplot object if
#' \code{combine = TRUE}; otherwise, a list of ggplot objects
#'
#' @export
#' @concept visualization
#'
#' @seealso \code{\link{FetchData}}
#'
#' @examples
#' data("pbmc_small")
#' VlnPlot(object = pbmc_small, features = 'PC_1')
#' VlnPlot(object = pbmc_small, features = 'LYZ', split.by = 'groups')
#'
VlnPlot <- function(
  object,
  features,
  cols = NULL,
  pt.size = NULL,
  alpha = 1,
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
  slot = deprecated(),
  layer = 'data',
  split.plot = FALSE,
  stack = FALSE,
  combine = TRUE,
  fill.by = 'feature',
  flip = FALSE,
  add.noise = TRUE,
  raster = NULL
) {
  if (is_present(arg = slot)) {
    deprecate_soft(
      when = '5.0.0',
      what = 'VlnPlot(slot = )',
      with = 'VlnPlot(layer = )'
    )
    layer <- slot %||% layer
  }
  if (
    !is.null(x = split.by) &
    getOption(x = 'Seurat.warn.vlnplot.split', default = TRUE)
  ) {
    message(
      "The default behaviour of split.by has changed.\n",
      "Separate violin plots are now plotted side-by-side.\n",
      "To restore the old behaviour of a single split violin,\n",
      "set split.plot = TRUE.
      \nThis message will be shown once per session."
    )
    options(Seurat.warn.vlnplot.split = FALSE)
  }
  return(ExIPlot(
    object = object,
    type = ifelse(test = split.plot, yes = 'splitViolin', no = 'violin'),
    features = features,
    idents = idents,
    ncol = ncol,
    sort = sort,
    assay = assay,
    y.max = y.max,
    same.y.lims = same.y.lims,
    adjust = adjust,
    pt.size = pt.size,
    alpha = alpha,
    cols = cols,
    group.by = group.by,
    split.by = split.by,
    log = log,
    layer = layer,
    stack = stack,
    combine = combine,
    fill.by = fill.by,
    flip = flip,
    add.noise = add.noise,
    raster = raster
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
#' @concept visualization
#'
#' @seealso \code{\link{DimPlot}}
#'
#' @examples
#' if (requireNamespace("ape", quietly = TRUE)) {
#'   data("pbmc_small")
#'   pbmc_small <- BuildClusterTree(object = pbmc_small, verbose = FALSE)
#'   PlotClusterTree(pbmc_small)
#'   ColorDimSplit(pbmc_small, node = 5)
#' }
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
  all.children <- sort(x = tree$edge[, 2][!tree$edge[, 2] %in% tree$edge[, 1]])
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
#' different colors and different shapes on cells.  Only applicable if \code{raster = FALSE}.
#' @param order Specify the order of plotting for the idents. This can be
#' useful for crowded plots if points of interest are being buried. Provide
#' either a full list of valid idents or a subset to be plotted last (on top)
#' @param shuffle Whether to randomly shuffle the order of points. This can be
#' useful for crowded plots if points of interest are being buried. (default is FALSE)
#' @param seed Sets the seed if randomly shuffling the order of points.
#' @param label Whether to label the clusters
#' @param label.size Sets size of labels
#' @param label.color Sets the color of the label text
#' @param label.box Whether to put a box around the label text (geom_text vs
#' geom_label)
#' @param alpha Alpha value for plotting (default is 1)
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
#' @param raster Convert points to raster format, default is \code{NULL} which
#' automatically rasterizes if plotting more than 100,000 cells
#' @param raster.dpi Pixel resolution for rasterized plots, passed to geom_scattermore().
#' Default is c(512, 512).
#'
#' @return A \code{\link[patchwork]{patchwork}ed} ggplot object if
#' \code{combine = TRUE}; otherwise, a list of ggplot objects
#'
#' @importFrom rlang !!
#' @importFrom ggplot2 facet_wrap vars sym labs
#' @importFrom patchwork wrap_plots
#'
#' @export
#' @concept visualization
#'
#' @note For the old \code{do.hover} and \code{do.identify} functionality, please see
#' \code{HoverLocator} and \code{CellSelector}, respectively.
#'
#' @aliases TSNEPlot PCAPlot ICAPlot
#' @seealso \code{\link{FeaturePlot}} \code{\link{HoverLocator}}
#' \code{\link{CellSelector}} \code{\link{FetchData}}
#'
#' @examples
#' data("pbmc_small")
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
  shuffle = FALSE,
  seed = 1,
  label = FALSE,
  label.size = 4,
  label.color = 'black',
  label.box = FALSE,
  repel = FALSE,
  alpha = 1,
  cells.highlight = NULL,
  cols.highlight = '#DE2D26',
  sizes.highlight = 1,
  na.value = 'grey50',
  ncol = NULL,
  combine = TRUE,
  raster = NULL,
  raster.dpi = c(512, 512)
) {
  if (!is_integerish(x = dims, n = 2L, finite = TRUE) || !all(dims > 0L)) {
    abort(message = "'dims' must be a two-length integer vector")
  }
  reduction <- reduction %||% DefaultDimReduc(object = object)
  # cells <- cells %||% colnames(x = object)
  ##### Cells for all cells in the assay.
  #### Cells function should not only get default layer
  cells <- cells %||% Cells(
    x = object,
    assay = DefaultAssay(object = object[[reduction]])
  )
  # data <- Embeddings(object = object[[reduction]])[cells, dims]
  # data <- as.data.frame(x = data)
  dims <- paste0(Key(object = object[[reduction]]), dims)
  orig.groups <- group.by
  group.by <- group.by %||% 'ident'
  data <- FetchData(
    object = object,
    vars = c(dims, group.by),
    cells = cells,
    clean = 'project'
  )
  # cells <- rownames(x = object)
  # object[['ident']] <- Idents(object = object)
  # orig.groups <- group.by
  # group.by <- group.by %||% 'ident'
  # data <- cbind(data, object[[group.by]][cells, , drop = FALSE])
  group.by <- colnames(x = data)[3:ncol(x = data)]
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
  if (isTRUE(x = shuffle)) {
    set.seed(seed = seed)
    data <- data[sample(x = 1:nrow(x = data)), ]
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
        alpha = alpha,
        label = FALSE,
        cells.highlight = cells.highlight,
        cols.highlight = cols.highlight,
        sizes.highlight = sizes.highlight,
        na.value = na.value,
        raster = raster,
        raster.dpi = raster.dpi
      )
      if (label) {
        plot <- LabelClusters(
          plot = plot,
          id = x,
          repel = repel,
          size = label.size,
          split.by = split.by,
          box = label.box,
          color = label.color
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
      plot <- if (is.null(x = orig.groups)) {
        plot + labs(title = NULL)
      } else {
        plot + CenterTitle()
      }
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
#' @param keep.scale How to handle the color scale across multiple plots. Options are:
#' \itemize{
#'   \item{"feature" (default; by row/feature scaling):}{ The plots for each individual feature are scaled to the maximum expression of the feature across the conditions provided to 'split.by'.}
#'   \item{"all" (universal scaling):}{ The plots for all features and conditions are scaled to the maximum expression value for the feature with the highest overall expression.}
#'   \item{NULL (no scaling):}{ Each individual plot is scaled to the maximum expression value of the feature in the condition provided to 'split.by'. Be aware setting NULL will result in color scales that are not comparable between plots.}
#' }
#' @param slot Which slot to pull expression data from?
#' @param blend Scale and blend expression values to visualize coexpression of two features
#' @param blend.threshold The color cutoff from weak signal to strong signal; ranges from 0 to 1.
#' @param ncol Number of columns to combine multiple feature plots to, ignored if \code{split.by} is not \code{NULL}
#' @param coord.fixed Plot cartesian coordinates with fixed aspect ratio
#' @param by.col If splitting by a factor, plot the splits per column with the features as rows; ignored if \code{blend = TRUE}
#' @param sort.cell Redundant with \code{order}. This argument is being
#' deprecated. Please use \code{order} instead.
#' @param interactive Launch an interactive \code{\link[Seurat:IFeaturePlot]{FeaturePlot}}
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
#' @concept visualization
#'
#' @note For the old \code{do.hover} and \code{do.identify} functionality, please see
#' \code{HoverLocator} and \code{CellSelector}, respectively.
#'
#' @aliases FeatureHeatmap
#' @seealso \code{\link{DimPlot}} \code{\link{HoverLocator}}
#' \code{\link{CellSelector}}
#'
#' @examples
#' data("pbmc_small")
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
  alpha = 1,
  order = FALSE,
  min.cutoff = NA,
  max.cutoff = NA,
  reduction = NULL,
  split.by = NULL,
  keep.scale = "feature",
  shape.by = NULL,
  slot = 'data',
  blend = FALSE,
  blend.threshold = 0.5,
  label = FALSE,
  label.size = 4,
  label.color = "black",
  repel = FALSE,
  ncol = NULL,
  coord.fixed = FALSE,
  by.col = TRUE,
  sort.cell = deprecated(),
  interactive = FALSE,
  combine = TRUE,
  raster = NULL,
  raster.dpi = c(512, 512)
) {
  # TODO: deprecate fully on 3.2.0
  if (is_present(arg = sort.cell)) {
    deprecate_stop(
      when = '4.9.0',
      what = 'FeaturePlot(sort.cell = )',
      with = 'FeaturePlot(order = )'
    )
  }
  if (isTRUE(x = interactive)) {
    return(IFeaturePlot(
      object = object,
      feature = features[1],
      dims = dims,
      reduction = reduction,
      slot = slot
    ))
  }
  # Check keep.scale param for valid entries
  if (!is.null(x = keep.scale)) {
    keep.scale <- arg_match0(arg = keep.scale, values = c('feature', 'all'))
  }
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
  if (!is_integerish(x = dims, n = 2L, finite = TRUE) && !all(dims > 0L)) {
    abort(message = "'dims' must be a two-length integer vector")
  }
  # Figure out blending stuff
  if (isTRUE(x = blend) && length(x = features) != 2) {
    abort(message = "Blending feature plots only works with two features")
  }
  # Set color scheme for blended FeaturePlots
  if (isTRUE(x = blend)) {
    default.colors <- eval(expr = formals(fun = FeaturePlot)$cols)
    cols <- switch(
      EXPR = as.character(x = length(x = cols)),
      '0' = {
        warn(message = "No colors provided, using default colors")
        default.colors
      },
      '1' = {
        warn(message = paste(
          "Only one color provided, assuming",
          sQuote(x = cols),
          "is double-negative and augmenting with default colors"
        ))
        c(cols, default.colors[2:3])
      },
      '2' = {
        warn(message = paste(
          "Only two colors provided, assuming specified are for features and agumenting with",
          sQuote(default.colors[1]),
          "for double-negatives",
        ))
        c(default.colors[1], cols)
      },
      '3' = cols,
      {
        warn(message = "More than three colors provided, using only first three")
        cols[1:3]
      }
    )
  }
  if (isTRUE(x = blend) && length(x = cols) != 3) {
    abort("Blending feature plots only works with three colors; first one for negative cells")
  }
  # Name the reductions
  dims <- paste0(Key(object = object[[reduction]]), dims)
  cells <- cells %||% Cells(x = object[[reduction]])
  # Get plotting data
  data <- FetchData(
    object = object,
    vars = c(dims, 'ident', features),
    cells = cells,
    slot = slot
  )
  # Check presence of features/dimensions
  if (ncol(x = data) < 4) {
    abort(message = paste(
      "None of the requested features were found:",
      paste(features, collapse = ', '),
      "in slot ",
      slot
    ))
  } else if (!all(dims %in% colnames(x = data))) {
    abort(message = "The dimensions requested were not found")
  }
  features <- setdiff(x = names(x = data), y = c(dims, 'ident'))
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
    abort(
      message = "There must be the same number of minimum and maximum cuttoffs as there are features"
    )
  }
  names(x = min.cutoff) <- names(x = max.cutoff) <- features
  brewer.gran <- ifelse(
    test = length(x = cols) == 1,
    yes = brewer.pal.info[cols, ]$maxcolors,
    no = length(x = cols)
  )
  # Apply cutoffs
  for (i in seq_along(along.with = features)) {
    f <- features[i]
    data.feature <- data[[f]]
    min.use <- SetQuantile(cutoff = min.cutoff[f], data = data.feature)
    max.use <- SetQuantile(cutoff = max.cutoff[f], data = data.feature)
    data.feature[data.feature < min.use] <- min.use
    data.feature[data.feature > max.use] <- max.use
    if (brewer.gran != 2) {
      data.feature <- if (all(data.feature == 0)) {
        rep_len(x = 0, length.out = length(x = data.feature))
      } else {
        as.numeric(x = as.factor(x = cut(
          x = as.numeric(x = data.feature),
          breaks = 2
        )))
      }
    }
    data[[f]] <- data.feature
  }
  # Figure out splits (FeatureHeatmap)
  data$split <- if (is.null(x = split.by)) {
    RandomName()
  } else {
    switch(
      EXPR = split.by,
      ident = Idents(object = object)[cells, drop = TRUE],
      object[[split.by, drop = TRUE]][cells, drop = TRUE]
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
    # Figure out which split we're working with
    ident <- levels(x = data$split)[i]
    data.plot <- data[as.character(x = data$split) == ident, , drop = FALSE]
    # Blend expression values
    if (isTRUE(x = blend)) {
      features <- features[1:2]
      no.expression <- features[colMeans(x = data.plot[, features]) == 0]
      if (length(x = no.expression) != 0) {
        abort(message = paste(
          "The following features have no value:",
          paste(no.expression, collapse = ', ')
        ))
      }
      data.plot <- cbind(data.plot[, c(dims, 'ident')], BlendExpression(data = data.plot[, features[1:2]]))
      features <- colnames(x = data.plot)[4:ncol(x = data.plot)]
    }
    # Make per-feature plots
    for (j in 1:length(x = features)) {
      feature <- features[j]
      # Get blended colors
      if (isTRUE(x = blend)) {
        cols.use <- as.numeric(x = as.character(x = data.plot[, feature])) + 1
        cols.use <- colors[[j]][sort(x = unique(x = cols.use))]
      } else {
        cols.use <- NULL
      }
      data.single <- data.plot[, c(dims, 'ident', feature, shape.by)]
      # Make the plot
      plot <- SingleDimPlot(
        data = data.single,
        dims = dims,
        col.by = feature,
        order = order,
        pt.size = pt.size,
        alpha = alpha,
        cols = cols.use,
        shape.by = shape.by,
        label = FALSE,
        raster = raster,
        raster.dpi = raster.dpi
      ) +
        scale_x_continuous(limits = xlims) +
        scale_y_continuous(limits = ylims) +
        theme_cowplot() +
        CenterTitle()
        # theme(plot.title = element_text(hjust = 0.5))
      # Add labels
      if (isTRUE(x = label)) {
        plot <- LabelClusters(
          plot = plot,
          id = 'ident',
          repel = repel,
          size = label.size,
          color = label.color
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
            warn(message = paste0(
              "All cells have the same value (",
              unique.feature.exp,
              ") of ",
              dQuote(x = feature)
            ))
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
      if (!(is.null(x = keep.scale)) && keep.scale == "feature" && !blend) {
        max.feature.value <- max(data[, feature])
        min.feature.value <- min(data[, feature])
        plot <- suppressMessages(plot & scale_color_gradientn(colors = cols, limits = c(min.feature.value, max.feature.value)))
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
  if (isTRUE(x = blend)) {
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
    test = is.null(x = split.by) || isTRUE(x = blend),
    yes = ncol,
    no = length(x = features)
  )
  legend <- if (isTRUE(x = blend)) {
    'none'
  } else {
    split.by %iff% 'none'
  }
  # Transpose the FeatureHeatmap matrix (not applicable for blended FeaturePlots)
  if (isTRUE(x = combine)) {
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
        ncol <- 1
        nrow <- nsplits
      } else {
        nrow <- split.by %iff% length(x = levels(x = data$split))
      }
      plots <- plots[c(do.call(
        what = rbind,
        args = split(
          x = 1:length(x = plots),
          f = ceiling(x = seq_along(along.with = 1:length(x = plots)) / length(x = features))
        )
      ))]
      # Set ncol to number of splits (nrow) and nrow to number of features (ncol)
      plots <- wrap_plots(plots, ncol = nrow, nrow = ncol)
      if (!is.null(x = legend) && legend == 'none') {
        plots <- plots & NoLegend()
      }
    } else {
      plots <- wrap_plots(plots, ncol = ncol, nrow = split.by %iff% length(x = levels(x = data$split)))
    }
    if (!is.null(x = legend) && legend == 'none') {
      plots <- plots & NoLegend()
    }
    if (!(is.null(x = keep.scale)) && keep.scale == "all" && !blend) {
      max.feature.value <- max(data[, features])
      min.feature.value <- min(data[, features])
      plots <- suppressMessages(plots & scale_color_gradientn(colors = cols, limits = c(min.feature.value, max.feature.value)))
    }
  }
  return(plots)
}

#' Visualize features in dimensional reduction space interactively
#'
#' @inheritParams FeaturePlot
#' @param feature Feature to plot
#'
#' @return Returns the final plot as a ggplot object
#'
#' @importFrom cowplot theme_cowplot
#' @importFrom ggplot2 theme element_text guides scale_color_gradientn
#' @importFrom miniUI miniPage miniButtonBlock miniTitleBarButton miniContentPanel
#' @importFrom shiny fillRow sidebarPanel selectInput plotOutput reactiveValues
#' observeEvent stopApp observe updateSelectInput renderPlot runGadget
#'
#' @export
#' @concept visualization
#'
IFeaturePlot <- function(object, feature, dims = c(1, 2), reduction = NULL, slot = 'data') {
  # Set initial data values
  feature.label <- 'Feature to visualize'
  assay.keys <- Key(object = object)[Assays(object = object)]
  keyed <- sapply(X = assay.keys, FUN = grepl, x = feature)
  assay <- if (any(keyed)) {
    names(x = which(x = keyed))[1]
  } else {
    DefaultAssay(object = object)
  }
  features <- sort(x = rownames(x = GetAssayData(
    object = object,
    slot = slot,
    assay = assay
  )))
  assays.use <- vapply(
    X = Assays(object = object),
    FUN = function(x) {
      return(!IsMatrixEmpty(x = GetAssayData(
        object = object,
        slot = slot,
        assay = x
      )))
    },
    FUN.VALUE = logical(length = 1L)
  )
  assays.use <- sort(x = Assays(object = object)[assays.use])
  reduction <- reduction %||% DefaultDimReduc(object = object)
  dims.reduc <- gsub(
    pattern = Key(object = object[[reduction]]),
    replacement = '',
    x = colnames(x = object[[reduction]])
  )
  # Set up the gadget UI
  ui <- miniPage(
    miniButtonBlock(miniTitleBarButton(
      inputId = 'done',
      label = 'Done',
      primary = TRUE
    )),
    miniContentPanel(
      fillRow(
        sidebarPanel(
          selectInput(
            inputId = 'assay',
            label = 'Assay',
            choices = assays.use,
            selected = assay,
            selectize = FALSE,
            width = '100%'
          ),
          selectInput(
            inputId = 'feature',
            label = feature.label,
            choices = features,
            selected = feature,
            selectize = FALSE,
            width = '100%'
          ),
          selectInput(
            inputId = 'reduction',
            label = 'Dimensional reduction',
            choices = Reductions(object = object),
            selected = reduction,
            selectize = FALSE,
            width = '100%'
          ),
          selectInput(
            inputId = 'xdim',
            label = 'X dimension',
            choices = dims.reduc,
            selected = as.character(x = dims[1]),
            selectize = FALSE,
            width = '100%'
          ),
          selectInput(
            inputId = 'ydim',
            label = 'Y dimension',
            choices = dims.reduc,
            selected = as.character(x = dims[2]),
            selectize = FALSE,
            width = '100%'
          ),
          selectInput(
            inputId = 'palette',
            label = 'Color scheme',
            choices = names(x = FeaturePalettes),
            selected = 'Seurat',
            selectize = FALSE,
            width = '100%'
          ),
          width = '100%'
        ),
        plotOutput(outputId = 'plot', height = '100%'),
        flex = c(1, 4)
      )
    )
  )
  # Prepare plotting data
  dims <- paste0(Key(object = object[[reduction]]), dims)
  plot.data <- FetchData(object = object, vars = c(dims, feature), slot = slot)
  # Shiny server
  server <- function(input, output, session) {
    plot.env <- reactiveValues(
      data = plot.data,
      dims = paste0(Key(object = object[[reduction]]), dims),
      feature = feature,
      palette = 'Seurat'
    )
    # Observe events
    observeEvent(
      eventExpr = input$done,
      handlerExpr = stopApp(returnValue = plot.env$plot)
    )
    observe(x = {
      assay <- input$assay
      feature.use <- input$feature
      features.assay <- sort(x = rownames(x = GetAssayData(
        object = object,
        slot = slot,
        assay = assay
      )))
      feature.use <- ifelse(
        test = feature.use %in% features.assay,
        yes = feature.use,
        no = features.assay[1]
      )
      reduc <- input$reduction
      dims.reduc <- gsub(
        pattern = Key(object = object[[reduc]]),
        replacement = '',
        x = colnames(x = object[[reduc]])
      )
      dims <- c(input$xdim, input$ydim)
      for (i in seq_along(along.with = dims)) {
        if (!dims[i] %in% dims.reduc) {
          dims[i] <- dims.reduc[i]
        }
      }
      updateSelectInput(
        session = session,
        inputId = 'xdim',
        label = 'X dimension',
        choices = dims.reduc,
        selected = as.character(x = dims[1])
      )
      updateSelectInput(
        session = session,
        inputId = 'ydim',
        label = 'Y dimension',
        choices = dims.reduc,
        selected = as.character(x = dims[2])
      )
      updateSelectInput(
        session = session,
        inputId = 'feature',
        label = feature.label,
        choices = features.assay,
        selected = feature.use
      )
    })
    observe(x = {
      feature.use <- input$feature
      feature.keyed <- paste0(Key(object = object[[input$assay]]), feature.use)
      reduc <- input$reduction
      dims <- c(input$xdim, input$ydim)
      dims <- paste0(Key(object = object[[reduc]]), dims)
      plot.data <- tryCatch(
        expr = FetchData(
          object = object,
          vars = c(dims, feature.keyed),
          slot = slot
        ),
        warning = function(...) {
          return(plot.env$data)
        },
        error = function(...) {
          return(plot.env$data)
        }
      )
      dims <- colnames(x = plot.data)[1:2]
      colnames(x = plot.data) <- c(dims, feature.use)
      plot.env$data <- plot.data
      plot.env$feature <- feature.use
      plot.env$dims <- dims
    })
    observe(x = {
      plot.env$palette <- input$palette
    })
    # Create the plot
    output$plot <- renderPlot(expr = {
      plot.env$plot <- SingleDimPlot(
        data = plot.env$data,
        dims = plot.env$dims,
        col.by = plot.env$feature,
        label = FALSE
      ) +
        theme_cowplot() +
        theme(plot.title = element_text(hjust = 0.5)) +
        guides(color = NULL) +
        scale_color_gradientn(
          colors = FeaturePalettes[[plot.env$palette]],
          guide = 'colorbar'
        )
      plot.env$plot
    })
  }
  runGadget(app = ui, server = server)
}

#' Highlight Neighbors in DimPlot
#'
#' It will color the query cells and the neighbors of the query cells in the
#' DimPlot
#'
#' @inheritParams DimPlot
#' @param nn.idx the neighbor index of all cells
#' @param query.cells cells used to find their neighbors
#' @param show.all.cells Show all cells or only query and neighbor cells
#'
#' @inherit DimPlot return
#'
#' @export
#' @concept visualization
#'
NNPlot <- function(
  object,
  reduction,
  nn.idx,
  query.cells,
  dims = 1:2,
  label = FALSE,
  label.size = 4,
  repel = FALSE,
  sizes.highlight = 2,
  pt.size = 1,
  cols.highlight = c("#377eb8", "#e41a1c"),
  na.value =  "#bdbdbd",
  order = c("self", "neighbors", "other"),
  show.all.cells = TRUE,
  ...
) {
  if (inherits(x = nn.idx, what = 'Neighbor')) {
    rownames(x = slot(object = nn.idx, name = 'nn.idx')) <- Cells(x = nn.idx)
    nn.idx <- Indices(object = nn.idx)
  }
  if (length(x = query.cells) > 1) {
    neighbor.cells <- apply(
      X = nn.idx[query.cells, -1],
      MARGIN = 2,
      FUN = function(x) {
        return(Cells(x = object)[x])
      }
    )
  } else {
    neighbor.cells <- Cells(x = object)[nn.idx[query.cells , -1]]
  }
  neighbor.cells <- as.vector(x = neighbor.cells)
  neighbor.cells <- neighbor.cells[!is.na(x = neighbor.cells)]
  object[["nn.col"]] <- "other"
  object[["nn.col"]][neighbor.cells, ] <- "neighbors"
  object[["nn.col"]][query.cells, ] <- "self"
  object$nn.col <- factor(
    x = object$nn.col,
    levels = c("self", "neighbors", "other")
  )
  if (!show.all.cells) {
    object <- subset(
      x = object,
      cells = Cells(x = object)[which(x = object[["nn.col"]] != "other")]
    )
    nn.cols  <- c(rev(x = cols.highlight))
    nn.pt.size <- sizes.highlight
  } else {
    highlight.info <- SetHighlight(
      cells.highlight = c(query.cells, neighbor.cells),
      cells.all = Cells(x = object),
      sizes.highlight = sizes.highlight,
      pt.size = pt.size,
      cols.highlight = "red"
    )
    nn.cols  <- c(na.value, rev(x = cols.highlight))
    nn.pt.size <- highlight.info$size
  }
  NN.plot <- DimPlot(
    object = object,
    reduction = reduction,
    dims = dims,
    group.by = "nn.col",
    cols = nn.cols,
    label = label,
    order =  order,
    pt.size = nn.pt.size ,
    label.size = label.size,
    repel = repel
  )
  return(NN.plot)
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
#' @inheritParams DimPlot
#' @param cell1 Cell 1 name
#' @param cell2 Cell 2 name
#' @param features Features to plot (default, all features)
#' @param highlight Features to highlight
#' @return A ggplot object
#'
#' @export
#' @concept visualization
#'
#' @aliases CellPlot
#'
#' @examples
#' data("pbmc_small")
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
  smooth = FALSE,
  raster = NULL,
  raster.dpi = c(512, 512)
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
    smooth = smooth,
    raster = raster,
    raster.dpi = raster.dpi
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
#' @param shuffle Whether to randomly shuffle the order of points. This can be
#' useful for crowded plots if points of interest are being buried. (default is FALSE)
#' @param seed Sets the seed if randomly shuffling the order of points.
#' @param group.by Name of one or more metadata columns to group (color) cells by
#' (for example, orig.ident); pass 'ident' to group by identity class
#' @param cols Colors to use for identity class plotting.
#' @param pt.size Size of the points on the plot
#' @param shape.by Ignored for now
#' @param span Spline span in loess function call, if \code{NULL}, no spline added
#' @param smooth Smooth the graph (similar to smoothScatter)
#' @param slot Slot to pull data from, should be one of 'counts', 'data', or 'scale.data'
#' @param combine Combine plots into a single \code{\link[patchwork]{patchwork}ed}
#' @param plot.cor Display correlation in plot title
#' @param raster Convert points to raster format, default is \code{NULL}
#' which will automatically use raster if the number of points plotted is greater than
#' 100,000
#' @param raster.dpi Pixel resolution for rasterized plots, passed to geom_scattermore().
#' Default is c(512, 512).
#' @param jitter Jitter for easier visualization of crowded points (default is FALSE)
#'
#' @return A ggplot object
#'
#' @importFrom ggplot2 geom_smooth aes_string
#' @importFrom patchwork wrap_plots
#'
#' @export
#' @concept visualization
#'
#' @aliases GenePlot
#'
#' @examples
#' data("pbmc_small")
#' FeatureScatter(object = pbmc_small, feature1 = 'CD9', feature2 = 'CD3E')
#'
FeatureScatter <- function(
  object,
  feature1,
  feature2,
  cells = NULL,
  shuffle = FALSE,
  seed = 1,
  group.by = NULL,
  cols = NULL,
  pt.size = 1,
  shape.by = NULL,
  span = NULL,
  smooth = FALSE,
  combine = TRUE,
  slot = 'data',
  plot.cor = TRUE,
  raster = NULL,
  raster.dpi = c(512, 512),
  jitter = FALSE
) {
  cells <- cells %||% colnames(x = object)
  if (isTRUE(x = shuffle)) {
    set.seed(seed = seed)
    cells <- sample(x = cells)
  }
  group.by <- group.by %||% 'ident'
  data <-  FetchData(
    object = object,
    vars = c(feature1, feature2, group.by),
    cells = cells,
    slot = slot
  )
  if (!grepl(pattern = feature1, x = names(x = data)[1])) {
    abort(message = paste("Feature 1", sQuote(x = feature1), "not found"))
  }
  if (!grepl(pattern = feature2, x = names(x = data)[2])) {
    abort(message = paste("Feature 2", sQuote(x = feature2), "not found"))
  }
  feature1 <-  names(x = data)[1]
  feature2 <-  names(x = data)[2]
  group.by <- intersect(x = group.by, y = names(x = data)[3:ncol(x = data)])
  for (group in group.by) {
    if (!is.factor(x = data[, group])) {
      data[, group] <- factor(x = data[, group])
    }
  }
  plots <- lapply(
    X = group.by,
    FUN = function(x) {
      SingleCorPlot(
        data = data[,c(feature1, feature2)],
        col.by = data[, x],
        cols = cols,
        pt.size = pt.size,
        smooth = smooth,
        legend.title = 'Identity',
        span = span,
        plot.cor = plot.cor,
        raster = raster,
        raster.dpi = raster.dpi,
        jitter = jitter
      )
    }
  )
  if (isTRUE(x = length(x = plots) == 1)) {
    return(plots[[1]])
  }
  if (isTRUE(x = combine)) {
    plots <- wrap_plots(plots, ncol = length(x = group.by))
  }
  return(plots)
}

#' View variable features
#'
#' @inheritParams FeatureScatter
#' @inheritParams SeuratObject::HVFInfo
#' @param cols Colors to specify non-variable/variable status
#' @param assay Assay to pull variable features from
#' @param log Plot the x-axis in log scale
#' @param raster Convert points to raster format, default is \code{NULL}
#' which will automatically use raster if the number of points plotted is greater than
#' 100,000
#'
#' @return A ggplot object
#'
#' @importFrom ggplot2 labs scale_color_manual scale_x_log10
#' @export
#' @concept visualization
#'
#' @aliases VariableGenePlot MeanVarPlot
#'
#' @seealso \code{\link{FindVariableFeatures}}
#'
#' @examples
#' data("pbmc_small")
#' VariableFeaturePlot(object = pbmc_small)
#'
VariableFeaturePlot <- function(
  object,
  cols = c('black', 'red'),
  pt.size = 1,
  log = NULL,
  selection.method = NULL,
  assay = NULL,
  raster = NULL,
  raster.dpi = c(512, 512)
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
  status.col <- colnames(hvf.info)[grepl("variable", colnames(hvf.info))][[1]]
  var.status <- c('no', 'yes')[unlist(hvf.info[[status.col]]) + 1]
  if (colnames(x = hvf.info)[3] == 'dispersion.scaled') {
    hvf.info <- hvf.info[, c(1, 2)]
  } else if (colnames(x = hvf.info)[3] == 'variance.expected') {
    hvf.info <- hvf.info[, c(1, 4)]
  } else {
    hvf.info <- hvf.info[, c(1, 3)]
  }
  axis.labels <- switch(
    EXPR = colnames(x = hvf.info)[2],
    'variance.standardized' = c('Average Expression', 'Standardized Variance'),
    'dispersion' = c('Average Expression', 'Dispersion'),
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
    pt.size = pt.size,
    raster = raster,
    raster.dpi = raster.dpi
  )
  if (length(x = unique(x = var.status)) == 1) {
    switch(
      EXPR = var.status[1],
      'yes' = {
        cols <- cols[2]
        labels.legend <- 'Variable'
      },
      'no' = {
        cols <- cols[1]
        labels.legend <- 'Non-variable'
      }
    )
  } else {
    labels.legend <- c('Non-variable', 'Variable')
  }
  plot <- plot +
    labs(title = NULL, x = axis.labels[1], y = axis.labels[2]) +
    scale_color_manual(
      labels = paste(labels.legend, 'count:', table(var.status)),
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
#' @concept visualization
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
#' Plot cells as polygons, rather than single points. Color cells by any value
#' accessible by \code{\link{FetchData}}.
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
#' @concept visualization
#' @concept spatial
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
# Spatial Plots
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Spatial Cluster Plots
#'
#' Visualize clusters or other categorical groupings in a spatial context
#'
#' @inheritParams DimPlot
#' @inheritParams SingleImagePlot
#' @param object A \code{\link[SeuratObject]{Seurat}} object
#' @param fov Name of FOV to plot
#' @param boundaries A vector of segmentation boundaries per image to plot;
#' can be a character vector, a named character vector, or a named list.
#' Names should be the names of FOVs and values should be the names of
#' segmentation boundaries
#' @param molecules A vector of molecules to plot
#' @param nmols Max number of each molecule specified in `molecules` to plot
#' @param dark.background Set plot background to black
#' @param crop Crop the plots to area with cells only
#' @param overlap Overlay boundaries from a single image to create a single
#' plot; if \code{TRUE}, then boundaries are stacked in the order they're
#' given (first is lowest)
#' @param axes Keep axes and panel background
#' @param combine Combine plots into a single
#' \code{patchwork} ggplot object.If \code{FALSE},
#' return a list of ggplot objects
#' @param coord.fixed Plot cartesian coordinates with fixed aspect ratio
#'
#' @return If \code{combine = TRUE}, a \code{patchwork}
#' ggplot object; otherwise, a list of ggplot objects
#'
#' @importFrom rlang !! is_na sym
#' @importFrom patchwork wrap_plots
#' @importFrom ggplot2 element_blank facet_wrap vars
#' @importFrom SeuratObject DefaultFOV Cells
#' DefaultBoundary FetchData Images Overlay
#'
#' @export
#'
ImageDimPlot <- function(
  object,
  fov = NULL,
  boundaries = NULL,
  group.by = NULL,
  split.by = NULL,
  cols = NULL,
  shuffle.cols = FALSE,
  size = 0.5,
  molecules = NULL,
  mols.size = 0.1,
  mols.cols = NULL,
  mols.alpha = 1.0,
  nmols = 1000,
  alpha = 1.0,
  border.color = 'white',
  border.size = NULL,
  na.value = 'grey50',
  dark.background = TRUE,
  crop = FALSE,
  cells = NULL,
  overlap = FALSE,
  axes = FALSE,
  combine = TRUE,
  coord.fixed = TRUE
) {
  cells <- cells %||% Cells(x = object)
  # Determine FOV to use
  fov <- fov %||% DefaultFOV(object = object)
  fov <- Filter(
    f = function(x) {
      return(
        x %in% Images(object = object) &&
          inherits(x = object[[x]], what = 'FOV')
      )
    },
    x = fov
  )
  if (!length(x = fov)) {
    stop("No compatible spatial coordinates present")
  }
  # Identify boundaries to use
  boundaries <- boundaries %||% sapply(
    X = fov,
    FUN = function(x) {
      return(DefaultBoundary(object = object[[x]]))
    },
    simplify = FALSE,
    USE.NAMES = TRUE
  )
  boundaries <- .BoundariesByImage(
    object = object,
    fov = fov,
    boundaries = boundaries
  )
  fov <- names(x = boundaries)
  overlap <- rep_len(x = overlap, length.out = length(x = fov))
  crop <- rep_len(x = crop, length.out = length(x = fov))
  names(x = crop) <- fov
  # Prepare plotting data
  group.by <- boundaries %!NA% group.by %||% 'ident'
  vars <- c(group.by, split.by)
  md <- if (!is_na(x = vars)) {
    FetchData(
      object = object,
      vars = vars[!is.na(x = vars)],
      cells = cells
    )
  } else {
    NULL
  }
  pnames <- unlist(x = lapply(
    X = seq_along(along.with = fov),
    FUN = function(i) {
      return(if (isTRUE(x = overlap[i])) {
        fov[i]
      } else {
        paste(fov[i], boundaries[[i]], sep = '_')
      })
    }
  ))
  pdata <- vector(mode = 'list', length = length(x = pnames))
  names(x = pdata) <- pnames
  for (i in names(x = pdata)) {
    ul <- unlist(x = strsplit(x = i, split = '_'))
    img <- paste(ul[1:length(ul)-1], collapse = '_')
    # Apply overlap
    lyr <- ul[length(ul)]
    if (is.na(x = lyr)) {
      lyr <- boundaries[[img]]
    }
    # TODO: Apply crop
    pdata[[i]] <- lapply(
      X = lyr,
      FUN = function(l) {
        if (l == 'NA') {
          return(NA)
        }
        df <- fortify(model = object[[img]][[l]])
        df <- df[df$cell %in% cells, , drop = FALSE]
        if (!is.null(x = md)) {
          df <- merge(x = df, y = md, by.x = 'cell', by.y = 0, all.x = TRUE)
        }
        df$cell <- paste(l, df$cell, sep = '_')
        df$boundary <- l
        return(df)
      }
    )
    pdata[[i]] <- if (!is_na(x = pdata[[i]])) {
      do.call(what = 'rbind', args = pdata[[i]])
    } else {
      unlist(x = pdata[[i]])
    }
  }
  # Fetch molecule information
  if (!is.null(x = molecules)) {
    molecules <- .MolsByFOV(
      object = object,
      fov = fov,
      molecules = molecules
    )
    mdata <- vector(mode = 'list', length = length(x = fov))
    names(x = mdata) <- fov
    for (img in names(x = mdata)) {
      idata <- object[[img]]
      if (!img %in% names(x = molecules)) {
        mdata[[img]] <- NULL
        next
      }
      if (isTRUE(x = crop[img])) {
        idata <- Overlay(x = idata, y = idata)
      }
      imols <- gsub(
        pattern = paste0('^', Key(object = idata)),
        replacement = '',
        x = molecules[[img]]
      )
      mdata[[img]] <- FetchData(
        object = idata,
        vars = imols,
        nmols = nmols
      )
    }
  } else {
    mdata <- NULL
  }
  # Build the plots
  plots <- vector(
    mode = 'list',
    length = length(x = pdata) * ifelse(
      test = length(x = group.by),
      yes = length(x = group.by),
      no = 1L
    )
  )
  idx <- 1L
  for (group in group.by) {
    for (i in seq_along(along.with = pdata)) {
      img <- unlist(x = strsplit(x = names(x = pdata)[i], split = '_'))[1L]
      p <- SingleImagePlot(
        data = pdata[[i]],
        col.by = pdata[[i]] %!NA% group,
        molecules = mdata[[img]],
        cols = cols,
        shuffle.cols = shuffle.cols,
        size = size,
        alpha = alpha,
        mols.size = mols.size,
        mols.cols = mols.cols,
        mols.alpha = mols.alpha,
        border.color = border.color,
        border.size = border.size,
        na.value = na.value,
        dark.background = dark.background
      )
      if (!is.null(x = split.by)) {
        p <- p + facet_wrap(
          facets = vars(!!sym(x = split.by))
        )
      }
      if (!isTRUE(x = axes)) {
        p <- p + NoAxes(panel.background = element_blank())
      }
      if (!anyDuplicated(x = pdata[[i]]$cell)) {
        p <- p + guides(fill = guide_legend(override.aes = list(size=4L, alpha=1)))
      }
      if (isTRUE(coord.fixed)) {
        p <- p + coord_fixed()
      }
      plots[[idx]] <- p
      idx <- idx + 1L
    }
  }
  if (isTRUE(x = combine)) {
    plots <- wrap_plots(plots)
  }
  return(plots)
}

#' Spatial Feature Plots
#'
#' Visualize expression in a spatial context
#'
#' @inheritParams FeaturePlot
#' @inheritParams ImageDimPlot
#' @param scale Set color scaling across multiple plots; choose from:
#' \itemize{
#'  \item \dQuote{\code{feature}}: Plots per-feature are scaled across splits
#'  \item \dQuote{\code{all}}: Plots per-feature are scaled across all features
#'  \item \dQuote{\code{none}}: Plots are not scaled; \strong{note}: setting
#'   \code{scale} to \dQuote{\code{none}} will result in color scales that are
#'   \emph{not} comparable between plots
#' }
#' Ignored if \code{blend = TRUE}
#'
#' @inherit ImageDimPlot return
#'
#' @importFrom patchwork wrap_plots
#' @importFrom cowplot theme_cowplot
#' @importFrom ggplot2 dup_axis element_blank element_text facet_wrap guides
#' labs margin vars scale_y_continuous theme
#' @importFrom SeuratObject DefaultFOV Cells DefaultBoundary
#' FetchData Images Overlay
#'
#' @export
#'
ImageFeaturePlot <- function(
  object,
  features,
  fov = NULL,
  boundaries = NULL,
  cols = if (isTRUE(x = blend)) {
    c("lightgrey", "#ff0000", "#00ff00")
  } else {
    c("lightgrey", "firebrick1")
  },
  size = 0.5,
  min.cutoff = NA,
  max.cutoff = NA,
  split.by = NULL,
  molecules = NULL,
  mols.size = 0.1,
  mols.cols = NULL,
  nmols = 1000,
  alpha = 1.0,
  border.color = 'white',
  border.size = NULL,
  dark.background = TRUE,
  blend = FALSE,
  blend.threshold = 0.5,
  crop = FALSE,
  cells = NULL,
  scale = c('feature', 'all', 'none'),
  overlap = FALSE,
  axes = FALSE,
  combine = TRUE,
  coord.fixed = TRUE
) {
  cells <- cells %||% Cells(x = object)
  scale <- scale[[1L]]
  scale <- match.arg(arg = scale)
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
  # Determine fov to use
  fov <- fov %||% DefaultFOV(object = object)
  fov <- Filter(
    f = function(x) {
      return(
        x %in% Images(object = object) &&
          inherits(x = object[[x]], what = 'FOV')
      )
    },
    x = fov
  )
  if (!length(x = fov)) {
    stop("No compatible spatial coordinates present")
  }
  # Identify boundaries to use
  boundaries <- boundaries %||% sapply(
    X = fov,
    FUN = function(x) {
      return(DefaultBoundary(object = object[[x]]))
    },
    simplify = FALSE,
    USE.NAMES = TRUE
  )
  boundaries <- .BoundariesByImage(
    object = object,
    fov = fov,
    boundaries = boundaries
  )
  fov <- names(x = boundaries)
  # Check overlaps/crops
  if (isTRUE(x = blend) || !is.null(x = split.by)) {
    type <- ifelse(test = isTRUE(x = 'blend'), yes = 'Blended', no = 'Split')
    if (length(x = fov) != 1L) {
      fov <- fov[1L]
      warning(
        type,
        ' image feature plots can only be done on a single image, using "',
        fov,
        '"',
        call. = FALSE,
        immediate. = TRUE
      )
    }
    if (any(!overlap) && length(x = boundaries[[fov]]) > 1L) {
      warning(
        type,
        " image feature plots require overlapped segmentations",
        call. = FALSE,
        immediate. = TRUE
      )
    }
    overlap <- TRUE
  }
  overlap <- rep_len(x = overlap, length.out = length(x = fov))
  crop <- rep_len(x = crop, length.out = length(x = fov))
  names(x = crop) <- names(x = overlap) <- fov
  # Checks for blending
  if (isTRUE(x = blend)) {
    if (length(x = features) != 2L) {
      stop("Blended feature plots only works with two features")
    }
    default.colors <- eval(expr = formals(fun = ImageFeaturePlot)$cols)
    cols <- switch(
      EXPR = as.character(x = length(x = cols)),
      '0' = {
        warning("No colors provided, using default colors", immediate. = TRUE)
        default.colors
      },
      '1' = {
        warning(
          "Only one color provided, assuming specified is double-negative and augmenting with default colors",
          immediate. = TRUE
        )
        c(cols, default.colors[2:3])
      },
      '2' = {
        warning(
          "Only two colors provided, assuming specified are for features and augmenting with '",
          default.colors[1],
          "' for double-negatives",
          immediate. = TRUE
        )
        c(default.colors[1], cols)
      },
      '3' = cols,
      {
        warning(
          "More than three colors provided, using only first three",
          immediate. = TRUE
        )
        cols[1:3]
      }
    )
  }
  # Get feature, splitting data
  md <- FetchData(
    object = object,
    vars = c(features, split.by[1L]),
    cells = cells
  )
  split.by <- intersect(x = split.by, y = colnames(x = md))
  if (!length(x = split.by)) {
    split.by <- NULL
  }
  imax <- ifelse(
    test = is.null(x = split.by),
    yes = ncol(x = md),
    no = ncol(x = md) - length(x = split.by)
  )
  features <- colnames(x = md)[1:imax]
  # Determine cutoffs
  min.cutoff <- mapply(
    FUN = function(cutoff, feature) {
      return(ifelse(
        test = is.na(x = cutoff),
        yes = min(md[[feature]]),
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
        yes = max(md[[feature]]),
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
  for (i in seq_along(along.with = features)) {
    f <- features[[i]]
    data.feature <- md[[f]]
    min.use <- SetQuantile(cutoff = min.cutoff[i], data = data.feature)
    max.use <- SetQuantile(cutoff = max.cutoff[i], data = data.feature)
    data.feature[data.feature < min.use] <- min.use
    data.feature[data.feature > max.use] <- max.use
    if (brewer.gran != 2) {
      data.feature <- if (all(data.feature == 0)) {
        rep_len(x = 0, length.out = length(x = data.feature))
      } else {
        as.numeric(x = as.factor(x = cut(
          x = as.numeric(x = data.feature),
          breaks = brewer.gran
        )))
      }
    }
    md[[f]] <- data.feature
  }
  # Figure out splits
  if (is.null(x = split.by)) {
    split.by <- RandomName()
    md[[split.by]] <- factor(x = split.by)
  }
  if (!is.factor(x = md[[split.by]])) {
    md[[split.by]] <- factor(x = md[[split.by]])
  }
  # Apply blends
  if (isTRUE(x = blend)) {
    md <- lapply(
      X = levels(x = md[[split.by]]),
      FUN = function(x) {
        df <- md[as.character(x = md[[split.by]]) == x, , drop = FALSE]
        no.expression <- features[colMeans(x = df[, features]) == 0]
        if (length(x = no.expression)) {
          stop(
            "The following features have no value: ",
            paste(no.expression, collapse = ', ')
          )
        }
        return(cbind(
          df[, split.by, drop = FALSE],
          BlendExpression(data = df[, features])
        ))
      }
    )
    md <- do.call(what = 'rbind', args = md)
    features <- setdiff(x = colnames(x = md), y = split.by)
  }
  # Prepare plotting data
  pnames <- unlist(x = lapply(
    X = seq_along(along.with = fov),
    FUN = function(i) {
      return(if (isTRUE(x = overlap[i])) {
        fov[i]
      } else {
        paste(fov[i], boundaries[[i]], sep = '_')
      })
    }
  ))
  pdata <- vector(mode = 'list', length = length(x = pnames))
  names(x = pdata) <- pnames
  for (i in names(x = pdata)) {
    ul <- unlist(x = strsplit(x = i, split = '_'))
    img <- paste(ul[1:length(ul)-1], collapse = '_')
    # Apply overlap
    lyr <- ul[length(ul)]
    if (is.na(x = lyr)) {
      lyr <- boundaries[[img]]
    }
    pdata[[i]] <- lapply(
      X = lyr,
      FUN = function(l) {
        df <- fortify(model = object[[img]][[l]])
        df <- df[df$cell %in% cells, , drop = FALSE]
        if (!is.null(x = md)) {
          df <- merge(x = df, y = md, by.x = 'cell', by.y = 0, all.x = TRUE)
        }
        df$cell <- paste(l, df$cell, sep = '_')
        df$boundary <- l
        return(df)
      }
    )
    pdata[[i]] <- if (!is_na(x = pdata[[i]])) {
      do.call(what = 'rbind', args = pdata[[i]])
    } else {
      unlist(x = pdata[[i]])
    }
  }
  # Fetch molecule information
  if (!is.null(x = molecules)) {
    molecules <- .MolsByFOV(
      object = object,
      fov = fov,
      molecules = molecules
    )
    mdata <- vector(mode = 'list', length = length(x = fov))
    names(x = mdata) <- fov
    for (img in names(x = mdata)) {
      idata <- object[[img]]
      if (!img %in% names(x = molecules)) {
        mdata[[img]] <- NULL
        next
      }
      if (isTRUE(x = crop[img])) {
        idata <- Overlay(x = idata, y = idata)
      }
      imols <- gsub(
        pattern = paste0('^', Key(object = idata)),
        replacement = '',
        x = molecules[[img]]
      )
      mdata[[img]] <- FetchData(
        object = idata,
        vars = imols,
        nmols = nmols
      )
    }
  } else {
    mdata <- NULL
  }
  # Set blended colors
  if (isTRUE(x = blend)) {
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
    blend.legend <- BlendMap(color.matrix = color.matrix)
  }
  limits <- switch(
    EXPR = scale,
    'all' = range(unlist(x = md[, features])),
    NULL
  )
  # Build the plots
  plots <- vector(
    mode = 'list',
    length = length(x = levels(x = md[[split.by]]))
  )
  names(x = plots) <- levels(x = md[[split.by]])
  for (i in seq_along(along.with = levels(x = md[[split.by]]))) {
    ident <- levels(x = md[[split.by]])[i]
    plots[[ident]] <- vector(mode = 'list', length = length(x = pdata))
    names(x = plots[[ident]]) <- names(x = pdata)
    if (isTRUE(x = blend)) {
      blend.key <- suppressMessages(
        expr = blend.legend +
          scale_y_continuous(
            sec.axis = dup_axis(name = ifelse(
              test = length(x = levels(x = md[[split.by]])) > 1,
              yes = ident,
              no = ''
            )),
            expand = c(0, 0)
          ) +
          labs(
            x = features[1L],
            y = features[2L],
            title = if (i == 1L) {
              paste('Color threshold:', blend.threshold)
            } else {
              NULL
            }
          ) +
          no.right
      )
    }
    for (j in seq_along(along.with = pdata)) {
      key <- names(x = pdata)[j]
      img <- unlist(x = strsplit(x = key, split = '_'))[1L]
      plots[[ident]][[key]] <- vector(
        mode = 'list',
        length = length(x = features) + ifelse(
          test = isTRUE(x = blend),
          yes = 1L,
          no = 0L
        )
      )
      data.plot <- pdata[[j]][as.character(x = pdata[[j]][[split.by]]) == ident, , drop = FALSE]
      for (y in seq_along(along.with = features)) {
        feature <- features[y]
        # Get blended colors
        cols.use <- if (isTRUE(x = blend)) {
          cc <- as.numeric(x = as.character(x = data.plot[, feature])) + 1
          colors[[y]][sort(unique(x = cc))]
        } else {
          NULL
        }
        colnames(data.plot) <- gsub("-", "_", colnames(data.plot))
        p <- SingleImagePlot(
          data = data.plot,
          col.by = gsub("-", "_", feature),
          size = size,
          col.factor = blend,
          cols = cols.use,
          molecules = mdata[[img]],
          mols.size = mols.size,
          mols.cols = mols.cols,
          alpha = alpha,
          border.color = border.color,
          border.size = border.size,
          dark.background = dark.background
        ) +
          CenterTitle() + labs(fill=feature)
        # Remove fill guides for blended plots
        if (isTRUE(x = blend)) {
          p <- p + guides(fill = 'none')
        }
        if (isTRUE(coord.fixed)) {
          p <- p + coord_fixed()
        }
        # Remove axes
        if (!isTRUE(x = axes)) {
          p <- p + NoAxes(panel.background = element_blank())
        } else if (isTRUE(x = blend) || length(x = levels(x = md[[split.by]])) > 1L) {
          if (y != 1L) {
            p <- p + theme(
              axis.line.y = element_blank(),
              axis.ticks.y = element_blank(),
              axis.text.y = element_blank(),
              axis.title.y.left = element_blank()
            )
          }
          if (i != length(x = levels(x = md[[split.by]]))) {
            p <- p + theme(
              axis.line.x = element_blank(),
              axis.ticks.x = element_blank(),
              axis.text.x = element_blank(),
              axis.title.x = element_blank()
            )
          }
        }
        # Add colors for unblended plots
        if (!isTRUE(x = blend)) {
          if (length(x = cols) == 1L) {
            p <- p + scale_fill_brewer(palette = cols)
          } else {
            cols.grad <- cols
            fexp <- data.plot[data.plot[[split.by]] == ident, feature, drop = TRUE]
            fexp <- unique(x = fexp)
            if (length(x = fexp) == 1L) {
              warning(
                "All cells have the same value (",
                fexp,
                ") of ",
                feature,
                call. = FALSE,
                immediate. = TRUE
              )
              if (fexp == 0) {
                cols.grad <- cols.grad[1L]
              }
            }
            # Check if we're scaling the colorbar across splits
            if (scale == 'feature') {
              limits <- range(pdata[[j]][[feature]])
            }
            p <- p + ggplot2::scale_fill_gradientn(
              colors = cols.grad,
              guide = 'colorbar',
              limits = limits
            )
          }
        }
        # Add some labels
        p <- p + if (i == 1L) {
          ggplot2::labs(title = feature)
        } else {
          ggplot2::labs(title = NULL)
        }
        plots[[ident]][[key]][[y]] <- p
      }
      if (isTRUE(x = blend)) {
        plots[[ident]][[key]][[length(x = plots[[ident]][[key]])]] <- blend.key
      } else if (length(x = levels(x = md[[split.by]])) > 1L) {
        plots[[ident]][[key]][[y]] <- suppressMessages(
          expr = plots[[ident]][[key]][[y]] +
            scale_y_continuous(sec.axis = dup_axis(name = ident)) +
            no.right
        )
      }
    }
    plots[[ident]] <- unlist(
      x = plots[[ident]],
      recursive = FALSE,
      use.names = FALSE
    )
  }
  plots <- unlist(x = plots, recursive = FALSE, use.names = FALSE)
  if (isTRUE(x = combine)) {
    if (isTRUE(x = blend) || length(x = levels(x = md[[split.by]])) > 1L) {
      plots <- wrap_plots(
        plots,
        ncol = ifelse(
          test = isTRUE(x = blend),
          yes = 4L,
          no = length(x = features)
        ),
        nrow = length(x = levels(x = md[[split.by]])),
        guides = 'collect'
      )
    } else {
      plots <- wrap_plots(plots)
    }
  }
  return(plots)
}

#' Visualize spatial and clustering (dimensional reduction) data in a linked,
#' interactive framework
#'
#' @inheritParams SpatialPlot
#' @inheritParams FeaturePlot
#' @inheritParams DimPlot
#' @param feature Feature to visualize
#' @param image Name of the image to use in the plot
#'
#' @return Returns final plots. If \code{combine}, plots are stiched together
#' using \code{\link{CombinePlots}}; otherwise, returns a list of ggplot objects
#'
#' @rdname LinkedPlots
#' @name LinkedPlots
#'
#' @importFrom scales hue_pal
#' @importFrom patchwork wrap_plots
#' @importFrom ggplot2 scale_alpha_ordinal guides
#' @importFrom miniUI miniPage gadgetTitleBar miniTitleBarButton miniContentPanel
#' @importFrom shiny fillRow plotOutput brushOpts clickOpts hoverOpts
#' verbatimTextOutput reactiveValues observeEvent stopApp nearPoints
#' brushedPoints renderPlot renderPrint runGadget
#'
#' @aliases LinkedPlot LinkedDimPlot
#'
#' @export
#' @concept visualization
#' @concept spatial
#'
#' @examples
#' \dontrun{
#' LinkedDimPlot(seurat.object)
#' LinkedFeaturePlot(seurat.object, feature = 'Hpca')
#' }
#'
LinkedDimPlot <- function(
  object,
  dims = 1:2,
  reduction = NULL,
  image = NULL,
  group.by = NULL,
  alpha = c(0.1, 1),
  combine = TRUE
) {
  # Setup gadget UI
  ui <- miniPage(
    gadgetTitleBar(
      title = 'LinkedDimPlot',
      left = miniTitleBarButton(inputId = 'reset', label = 'Reset')
    ),
    miniContentPanel(
      fillRow(
        plotOutput(
          outputId = 'spatialplot',
          height = '100%',
          # brush = brushOpts(id = 'brush', delay = 10, clip = TRUE, resetOnNew = FALSE),
          click = clickOpts(id = 'spclick', clip = TRUE),
          hover = hoverOpts(id = 'sphover', delay = 10, nullOutside = TRUE)
        ),
        plotOutput(
          outputId = 'dimplot',
          height = '100%',
          brush = brushOpts(id = 'brush', delay = 10, clip = TRUE, resetOnNew = FALSE),
          click = clickOpts(id = 'dimclick', clip = TRUE),
          hover = hoverOpts(id = 'dimhover', delay = 10, nullOutside = TRUE)
        ),
        height = '97%'
      ),
      verbatimTextOutput(outputId = 'info')
    )
  )
  # Prepare plotting data
  image <- image %||% DefaultImage(object = object)
  cells.use <- Cells(x = object[[image]])
  reduction <- reduction %||% DefaultDimReduc(object = object)
  dims <- dims[1:2]
  dims <- paste0(Key(object = object[[reduction]]), dims)
  group.by <- group.by %||% 'ident'
  group.data <- FetchData(
    object = object,
    vars = group.by,
    cells = cells.use
  )
  coords <- GetTissueCoordinates(object = object[[image]])
  embeddings <- Embeddings(object = object[[reduction]])[cells.use, dims]
  plot.data <- cbind(coords, group.data, embeddings)
  plot.data$selected_ <- FALSE
  Idents(object = object) <- group.by
  # Setup the server
  server <- function(input, output, session) {
    click <- reactiveValues(pt = NULL, invert = FALSE)
    plot.env <- reactiveValues(data = plot.data, alpha.by = NULL)
    # Handle events
    observeEvent(
      eventExpr = input$done,
      handlerExpr = {
        plots <- list(plot.env$spatialplot, plot.env$dimplot)
        if (combine) {
          plots <- wrap_plots(plots, ncol = 2)
        }
        stopApp(returnValue = plots)
      }
    )
    observeEvent(
      eventExpr = input$reset,
      handlerExpr = {
        click$pt <- NULL
        click$invert <- FALSE
        session$resetBrush(brushId = 'brush')
      }
    )
    observeEvent(eventExpr = input$brush, handlerExpr = click$pt <- NULL)
    observeEvent(
      eventExpr = input$spclick,
      handlerExpr = {
        click$pt <- input$spclick
        click$invert <- TRUE
      }
    )
    observeEvent(
      eventExpr = input$dimclick,
      handlerExpr = {
        click$pt <- input$dimclick
        click$invert <- FALSE
      }
    )
    observeEvent(
      eventExpr = c(input$brush, input$spclick, input$dimclick),
      handlerExpr = {
        plot.env$data <- if (is.null(x = input$brush)) {
          clicked <- nearPoints(
            df = plot.data,
            coordinfo = if (click$invert) {
              InvertCoordinate(x = click$pt)
            } else {
              click$pt
            },
            threshold = 10,
            maxpoints = 1
          )
          if (nrow(x = clicked) == 1) {
            cell.clicked <- rownames(x = clicked)
            group.clicked <- plot.data[cell.clicked, group.by, drop = TRUE]
            idx.group <- which(x = plot.data[[group.by]] == group.clicked)
            plot.data[idx.group, 'selected_'] <- TRUE
            plot.data
          } else {
            plot.data
          }
        } else if (input$brush$outputId == 'dimplot') {
          brushedPoints(df = plot.data, brush = input$brush, allRows = TRUE)
        } else if (input$brush$outputId == 'spatialplot') {
          brushedPoints(df = plot.data, brush = InvertCoordinate(x = input$brush), allRows = TRUE)
        }
        plot.env$alpha.by <- if (any(plot.env$data$selected_)) {
          'selected_'
        } else {
          NULL
        }
      }
    )
    # Set plots
    output$spatialplot <- renderPlot(
      expr = {
        plot.env$spatialplot <- SingleSpatialPlot(
          data = plot.env$data,
          image = object[[image]],
          col.by = group.by,
          pt.size.factor = 1.6,
          crop = TRUE,
          alpha.by = plot.env$alpha.by
        ) + scale_alpha_ordinal(range = alpha) + NoLegend()
        plot.env$spatialplot
      }
    )
    output$dimplot <- renderPlot(
      expr = {
        plot.env$dimplot <- SingleDimPlot(
          data = plot.env$data,
          dims = dims,
          col.by = group.by,
          alpha.by = plot.env$alpha.by
        ) + scale_alpha_ordinal(range = alpha) + guides(alpha = "none")
        plot.env$dimplot
      }
    )
    # Add hover text
    output$info <- renderPrint(
      expr = {
        cell.hover <- rownames(x = nearPoints(
          df = plot.data,
          coordinfo = if (is.null(x = input[['sphover']])) {
            input$dimhover
          } else {
            InvertCoordinate(x = input$sphover)
          },
          threshold = 10,
          maxpoints = 1
        ))
        # if (length(x = cell.hover) == 1) {
        #   palette <- hue_pal()(n = length(x = levels(x = object)))
        #   group <- plot.data[cell.hover, group.by, drop = TRUE]
        #   background <- palette[which(x = levels(x = object) == group)]
        #   text <- unname(obj = BGTextColor(background = background))
        #   style <- paste0(
        #     paste(
        #       paste('background-color:', background),
        #       paste('color:', text),
        #       sep = '; '
        #     ),
        #     ';'
        #   )
        #   info <- paste(cell.hover, paste('Group:', group), sep = '<br />')
        # } else {
        #   style <- 'background-color: white; color: black'
        #   info <- NULL
        # }
        # HTML(text = paste0("<div style='", style, "'>", info, "</div>"))
        # p(HTML(info), style = style)
        # paste0('<div style="', style, '">', info, '</div>')
        # TODO: Get newlines, extra information, and background color working
        if (length(x = cell.hover) == 1) {
          paste(cell.hover, paste('Group:', plot.data[cell.hover, group.by, drop = TRUE]), collapse = '<br />')
        } else {
          NULL
        }
      }
    )
  }
  # Run the thang
  runGadget(app = ui, server = server)
}

#' @rdname LinkedPlots
#'
#' @aliases LinkedFeaturePlot
#'
#' @importFrom ggplot2 scale_fill_gradientn theme scale_alpha guides
#' scale_color_gradientn guide_colorbar
#'
#' @export
#' @concept visualization
#' @concept spatial
LinkedFeaturePlot <- function(
  object,
  feature,
  dims = 1:2,
  reduction = NULL,
  image = NULL,
  slot = 'data',
  alpha = c(0.1, 1),
  combine = TRUE
) {
  # Setup gadget UI
  ui <- miniPage(
    gadgetTitleBar(
      title = 'LinkedFeaturePlot',
      left = NULL
    ),
    miniContentPanel(
      fillRow(
        plotOutput(
          outputId = 'spatialplot',
          height = '100%',
          hover = hoverOpts(id = 'sphover', delay = 10, nullOutside = TRUE)
        ),
        plotOutput(
          outputId = 'dimplot',
          height = '100%',
          hover = hoverOpts(id = 'dimhover', delay = 10, nullOutside = TRUE)
        ),
        height = '97%'
      ),
      verbatimTextOutput(outputId = 'info')
    )
  )
  # Prepare plotting data
  cols <- SpatialColors(n = 100)
  image <- image %||% DefaultImage(object = object)
  cells.use <- Cells(x = object[[image]])
  reduction <- reduction %||% DefaultDimReduc(object = object)
  dims <- dims[1:2]
  dims <- paste0(Key(object = object[[reduction]]), dims)
  group.data <- FetchData(
    object = object,
    vars = feature,
    cells = cells.use
  )
  coords <- GetTissueCoordinates(object = object[[image]])
  embeddings <- Embeddings(object = object[[reduction]])[cells.use, dims]
  plot.data <- cbind(coords, group.data, embeddings)
  # Setup the server
  server <- function(input, output, session) {
    plot.env <- reactiveValues()
    # Handle events
    observeEvent(
      eventExpr = input$done,
      handlerExpr = {
        plots <- list(plot.env$spatialplot, plot.env$dimplot)
        if (combine) {
          plots <- wrap_plots(plots, ncol = 2)
        }
        stopApp(returnValue = plots)
      }
    )
    # Set plots
    output$spatialplot <- renderPlot(
      expr = {
        plot.env$spatialplot <- SingleSpatialPlot(
          data = plot.data,
          image = object[[image]],
          col.by = feature,
          pt.size.factor = 1.6,
          crop = TRUE,
          alpha.by = feature
        ) +
          scale_fill_gradientn(name = feature, colours = cols) +
          theme(legend.position = 'top') +
          scale_alpha(range = alpha) +
          guides(alpha = "none")
        plot.env$spatialplot
      }
    )
    output$dimplot <- renderPlot(
      expr = {
        plot.env$dimplot <- SingleDimPlot(
          data = plot.data,
          dims = dims,
          col.by = feature
        ) +
          scale_color_gradientn(name = feature, colours = cols, guide = 'colorbar') +
          guides(color = guide_colorbar())
        plot.env$dimplot
      }
    )
    # Add hover text
    output$info <- renderPrint(
      expr = {
        cell.hover <- rownames(x = nearPoints(
          df = plot.data,
          coordinfo = if (is.null(x = input[['sphover']])) {
            input$dimhover
          } else {
            InvertCoordinate(x = input$sphover)
          },
          threshold = 10,
          maxpoints = 1
        ))
        # TODO: Get newlines, extra information, and background color working
        if (length(x = cell.hover) == 1) {
          paste(cell.hover, paste('Expression:', plot.data[cell.hover, feature, drop = TRUE]), collapse = '<br />')
        } else {
          NULL
        }
      }
    )
  }
  runGadget(app = ui, server = server)
}

#' Visualize clusters spatially and interactively
#'
#' @inheritParams SpatialPlot
#' @inheritParams DimPlot
#' @inheritParams LinkedPlots
#'
#' @return Returns final plot as a ggplot object
#'
#' @importFrom ggplot2 scale_alpha_ordinal
#' @importFrom miniUI miniPage miniButtonBlock miniTitleBarButton miniContentPanel
#' @importFrom shiny fillRow plotOutput verbatimTextOutput reactiveValues
#' observeEvent stopApp nearPoints renderPlot runGadget
#'
#' @export
#' @concept visualization
#' @concept spatial
#'
ISpatialDimPlot <- function(
  object,
  image = NULL,
  group.by = NULL,
  alpha = c(0.3, 1)
) {
  # Setup gadget UI
  ui <- miniPage(
    miniButtonBlock(miniTitleBarButton(
      inputId = 'done',
      label = 'Done',
      primary = TRUE
    )),
    miniContentPanel(
      fillRow(
        plotOutput(
          outputId = 'plot',
          height = '100%',
          click = clickOpts(id = 'click', clip = TRUE),
          hover = hoverOpts(id = 'hover', delay = 10, nullOutside = TRUE)
        ),
        height = '97%'
      ),
      verbatimTextOutput(outputId = 'info')
    )
  )
  # Get plotting data
  # Prepare plotting data
  image <- image %||% DefaultImage(object = object)
  cells.use <- Cells(x = object[[image]])
  group.by <- group.by %||% 'ident'
  group.data <- FetchData(
    object = object,
    vars = group.by,
    cells = cells.use
  )
  coords <- GetTissueCoordinates(object = object[[image]])
  plot.data <- cbind(coords, group.data)
  plot.data$selected_ <- FALSE
  Idents(object = object) <- group.by
  # Set up the server
  server <- function(input, output, session) {
    click <- reactiveValues(pt = NULL)
    plot.env <- reactiveValues(data = plot.data, alpha.by = NULL)
    # Handle events
    observeEvent(
      eventExpr = input$done,
      handlerExpr = stopApp(returnValue = plot.env$plot)
    )
    observeEvent(
      eventExpr = input$click,
      handlerExpr = {
        clicked <- nearPoints(
          df = plot.data,
          coordinfo = InvertCoordinate(x = input$click),
          threshold = 10,
          maxpoints = 1
        )
        plot.env$data <- if (nrow(x = clicked) == 1) {
          cell.clicked <- rownames(x = clicked)
          cell.clicked <- rownames(x = clicked)
          group.clicked <- plot.data[cell.clicked, group.by, drop = TRUE]
          idx.group <- which(x = plot.data[[group.by]] == group.clicked)
          plot.data[idx.group, 'selected_'] <- TRUE
          plot.data
        } else {
          plot.data
        }
        plot.env$alpha.by <- if (any(plot.env$data$selected_)) {
          'selected_'
        } else {
          NULL
        }
      }
    )
    # Set plot
    output$plot <- renderPlot(
      expr = {
        plot.env$plot <- SingleSpatialPlot(
          data = plot.env$data,
          image = object[[image]],
          col.by = group.by,
          crop = TRUE,
          alpha.by = plot.env$alpha.by,
          pt.size.factor = 1.6
        ) + scale_alpha_ordinal(range = alpha) + NoLegend()
        plot.env$plot
      }
    )
    # Add hover text
    output$info <- renderPrint(
      expr = {
        cell.hover <- rownames(x = nearPoints(
          df = plot.data,
          coordinfo = InvertCoordinate(x = input$hover),
          threshold = 10,
          maxpoints = 1
        ))
        if (length(x = cell.hover) == 1) {
          paste(cell.hover, paste('Group:', plot.data[cell.hover, group.by, drop = TRUE]), collapse = '<br />')
        } else {
          NULL
        }
      }
    )
  }
  runGadget(app = ui, server = server)
}

#' Visualize features spatially and interactively
#'
#' @inheritParams SpatialPlot
#' @inheritParams FeaturePlot
#' @inheritParams LinkedPlots
#'
#' @return Returns final plot as a ggplot object
#'
#' @importFrom ggplot2 scale_fill_gradientn theme scale_alpha guides
#' @importFrom miniUI miniPage miniButtonBlock miniTitleBarButton miniContentPanel
#' @importFrom shiny fillRow sidebarPanel sliderInput selectInput reactiveValues
#' observeEvent stopApp observe updateSelectInput plotOutput renderPlot runGadget
#'
#' @export
#' @concept visualization
#' @concept spatial
ISpatialFeaturePlot <- function(
  object,
  feature,
  image = NULL,
  slot = 'data',
  alpha = c(0.1, 1)
) {
  # Set inital data values
  assay.keys <- Key(object = object)[Assays(object = object)]
  keyed <- sapply(X = assay.keys, FUN = grepl, x = feature)
  assay <- if (any(keyed)) {
    names(x = which(x = keyed))[1]
  } else {
    DefaultAssay(object = object)
  }
  features <- sort(x = rownames(x = GetAssayData(
    object = object,
    slot = slot,
    assay = assay
  )))
  feature.label <- 'Feature to visualize'
  assays.use <- vapply(
    X = Assays(object = object),
    FUN = function(x) {
      return(!IsMatrixEmpty(x = GetAssayData(
        object = object,
        slot = slot,
        assay = x
      )))
    },
    FUN.VALUE = logical(length = 1L)
  )
  assays.use <- sort(x = Assays(object = object)[assays.use])
  # Setup gadget UI
  ui <- miniPage(
    miniButtonBlock(miniTitleBarButton(
      inputId = 'done',
      label = 'Done',
      primary = TRUE
    )),
    miniContentPanel(
      fillRow(
        sidebarPanel(
          sliderInput(
            inputId = 'alpha',
            label = 'Alpha intensity',
            min = 0,
            max = max(alpha),
            value = min(alpha),
            step = 0.01,
            width = '100%'
          ),
          sliderInput(
            inputId = 'pt.size',
            label = 'Point size',
            min = 0,
            max = 5,
            value = 1.6,
            step = 0.1,
            width = '100%'
          ),
          selectInput(
            inputId = 'assay',
            label = 'Assay',
            choices = assays.use,
            selected = assay,
            selectize = FALSE,
            width = '100%'
          ),
          selectInput(
            inputId = 'feature',
            label = feature.label,
            choices = features,
            selected = feature,
            selectize = FALSE,
            width = '100%'
          ),
          selectInput(
            inputId = 'palette',
            label = 'Color scheme',
            choices = names(x = FeaturePalettes),
            selected = 'Spatial',
            selectize = FALSE,
            width = '100%'
          ),
          width = '100%'
        ),
        plotOutput(outputId = 'plot', height = '100%'),
        flex = c(1, 4)
      )
    )
  )
  # Prepare plotting data
  image <- image %||% DefaultImage(object = object)
  cells.use <- Cells(x = object[[image]])
  coords <- GetTissueCoordinates(object = object[[image]])
  feature.data <- FetchData(
    object = object,
    vars = feature,
    cells = cells.use,
    slot = slot
  )
  plot.data <- cbind(coords, feature.data)
  server <- function(input, output, session) {
    plot.env <- reactiveValues(
      data = plot.data,
      feature = feature,
      palette = 'Spatial'
    )
    # Observe events
    observeEvent(
      eventExpr = input$done,
      handlerExpr = stopApp(returnValue = plot.env$plot)
    )
    observe(x = {
      assay <- input$assay
      feature.use <- input$feature
      features.assay <- sort(x = rownames(x = GetAssayData(
        object = object,
        slot = slot,
        assay = assay
      )))
      feature.use <- ifelse(
        test = feature.use %in% features.assay,
        yes = feature.use,
        no = features.assay[1]
      )
      updateSelectInput(
        session = session,
        inputId = 'assay',
        label = 'Assay',
        choices = assays.use,
        selected = assay
      )
      updateSelectInput(
        session = session,
        inputId = 'feature',
        label = feature.label,
        choices = features.assay,
        selected = feature.use
      )
    })
    observe(x = {
      feature.use <- input$feature
      try(
        expr = {
          feature.data <- FetchData(
            object = object,
            vars = paste0(Key(object = object[[input$assay]]), feature.use),
            cells = cells.use,
            slot = slot
          )
          colnames(x = feature.data) <- feature.use
          plot.env$data <- cbind(coords, feature.data)
          plot.env$feature <- feature.use
        },
        silent = TRUE
      )
    })
    observe(x = {
      plot.env$palette <- input$palette
    })
    # Create plot
    output$plot <- renderPlot(expr = {
      plot.env$plot <- SingleSpatialPlot(
        data = plot.env$data,
        image = object[[image]],
        col.by = plot.env$feature,
        pt.size.factor = input$pt.size,
        crop = TRUE,
        alpha.by = plot.env$feature
      ) +
        # scale_fill_gradientn(name = plot.env$feature, colours = cols) +
        scale_fill_gradientn(name = plot.env$feature, colours = FeaturePalettes[[plot.env$palette]]) +
        theme(legend.position = 'top') +
        scale_alpha(range = c(input$alpha, 1)) +
        guides(alpha = "none")
      plot.env$plot
    })
  }
  runGadget(app = ui, server = server)
}

#' Visualize spatial clustering and expression data.
#'
#' SpatialPlot plots a feature or discrete grouping (e.g. cluster assignments) as
#' spots over the image that was collected. We also provide SpatialFeaturePlot
#' and SpatialDimPlot as wrapper functions around SpatialPlot for a consistent
#' naming framework.
#'
#' @inheritParams HoverLocator
#' @param object A Seurat object
#' @param group.by Name of meta.data column to group the data by
#' @param features Name of the feature to visualize. Provide either group.by OR
#' features, not both.
#' @param images Name of the images to use in the plot(s)
#' @param cols Vector of colors, each color corresponds to an identity class.
#' This may also be a single character or numeric value corresponding to a
#' palette as specified by \code{\link[RColorBrewer]{brewer.pal.info}}. By
#' default, ggplot2 assigns colors
#' @param image.alpha Adjust the opacity of the background images. Set to 0 to
#' remove.
#' @param crop Crop the plot in to focus on points plotted. Set to \code{FALSE} to show
#' entire background image.
#' @param slot If plotting a feature, which data slot to pull from (counts,
#' data, or scale.data)
#' @param min.cutoff,max.cutoff Vector of minimum and maximum cutoff
#' values for each feature, may specify quantile in the form of 'q##' where '##'
#' is the quantile (eg, 'q1', 'q10')
#' @param cells.highlight A list of character or numeric vectors of cells to
#' highlight. If only one group of cells desired, can simply pass a vector
#' instead of a list. If set, colors selected cells to the color(s) in
#' cols.highlight
#' @param cols.highlight A vector of colors to highlight the cells as; ordered
#' the same as the groups in cells.highlight; last color corresponds to
#' unselected cells.
#' @param facet.highlight When highlighting certain groups of cells, split each
#' group into its own plot
#' @param label Whether to label the clusters
#' @param label.size Sets the size of the labels
#' @param label.color Sets the color of the label text
#' @param label.box Whether to put a box around the label text (geom_text vs
#' geom_label)
#' @param repel Repels the labels to prevent overlap
#' @param ncol Number of columns if plotting multiple plots
#' @param combine Combine plots into a single gg object; note that if TRUE;
#' themeing will not work when plotting multiple features/groupings
#' @param pt.size.factor Scale the size of the spots.
#' @param alpha Controls opacity of spots. Provide as a vector specifying the
#' min and max for SpatialFeaturePlot. For SpatialDimPlot, provide a single
#' alpha value for each plot.
#' @param stroke Control the width of the border around the spots
#' @param interactive Launch an interactive SpatialDimPlot or SpatialFeaturePlot
#' session, see \code{\link{ISpatialDimPlot}} or
#' \code{\link{ISpatialFeaturePlot}} for more details
#' @param do.identify,do.hover DEPRECATED in favor of \code{interactive}
#' @param identify.ident DEPRECATED
#'
#' @return If \code{do.identify}, either a vector of cells selected or the object
#' with selected cells set to the value of \code{identify.ident} (if set). Else,
#' if \code{do.hover}, a plotly object with interactive graphics. Else, a ggplot
#' object
#'
#' @importFrom ggplot2 scale_fill_gradientn ggtitle theme element_text scale_alpha
#' @importFrom patchwork wrap_plots
#' @export
#' @concept visualization
#' @concept spatial
#'
#' @examples
#' \dontrun{
#' # For functionality analagous to FeaturePlot
#' SpatialPlot(seurat.object, features = "MS4A1")
#' SpatialFeaturePlot(seurat.object, features = "MS4A1")
#'
#' # For functionality analagous to DimPlot
#' SpatialPlot(seurat.object, group.by = "clusters")
#' SpatialDimPlot(seurat.object, group.by = "clusters")
#' }
#'
SpatialPlot <- function(
  object,
  group.by = NULL,
  features = NULL,
  images = NULL,
  cols = NULL,
  image.alpha = 1,
  crop = TRUE,
  slot = 'data',
  min.cutoff = NA,
  max.cutoff = NA,
  cells.highlight = NULL,
  cols.highlight = c('#DE2D26', 'grey50'),
  facet.highlight = FALSE,
  label = FALSE,
  label.size = 5,
  label.color = 'white',
  label.box = TRUE,
  repel = FALSE,
  ncol = NULL,
  combine = TRUE,
  pt.size.factor = 1.6,
  alpha = c(1, 1),
  stroke = 0.25,
  interactive = FALSE,
  do.identify = FALSE,
  identify.ident = NULL,
  do.hover = FALSE,
  information = NULL
) {
  if (isTRUE(x = do.hover) || isTRUE(x = do.identify)) {
    warning(
      "'do.hover' and 'do.identify' are deprecated as we are removing plotly-based interactive graphics, use 'interactive' instead for Shiny-based interactivity",
      call. = FALSE,
      immediate. = TRUE
    )
    interactive <- TRUE
  }
  if (!is.null(x = group.by) & !is.null(x = features)) {
    stop("Please specific either group.by or features, not both.")
  }
  images <- images %||% Images(object = object, assay = DefaultAssay(object = object))
  if (length(x = images) == 0) {
    images <- Images(object = object)
  }
  if (length(x = images) < 1) {
    stop("Could not find any spatial image information")
  }
  if (is.null(x = features)) {
    if (interactive) {
      return(ISpatialDimPlot(
        object = object,
        image = images[1],
        group.by = group.by,
        alpha = alpha
      ))
    }
    group.by <- group.by %||% 'ident'
    object[['ident']] <- Idents(object = object)
    data <- object[[group.by]]
    for (group in group.by) {
      if (!is.factor(x = data[, group])) {
        data[, group] <- factor(x = data[, group])
      }
    }
  } else {
    if (interactive) {
      return(ISpatialFeaturePlot(
        object = object,
        feature = features[1],
        image = images[1],
        slot = slot,
        alpha = alpha
      ))
    }
    data <- FetchData(
      object = object,
      vars = features,
      slot = slot
    )
    features <- colnames(x = data)
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
    # Apply cutoffs
    data <- sapply(
      X = 1:ncol(x = data),
      FUN = function(index) {
        data.feature <- as.vector(x = data[, index])
        min.use <- SetQuantile(cutoff = min.cutoff[index], data.feature)
        max.use <- SetQuantile(cutoff = max.cutoff[index], data.feature)
        data.feature[data.feature < min.use] <- min.use
        data.feature[data.feature > max.use] <- max.use
        return(data.feature)
      }
    )
    colnames(x = data) <- features
    rownames(x = data) <- Cells(x = object)
  }
  features <- colnames(x = data)
  colnames(x = data) <- features
  rownames(x = data) <- colnames(x = object)
  facet.highlight <- facet.highlight && (!is.null(x = cells.highlight) && is.list(x = cells.highlight))
  if (do.hover) {
    if (length(x = images) > 1) {
      images <- images[1]
      warning(
        "'do.hover' requires only one image, using image ",
        images,
        call. = FALSE,
        immediate. = TRUE
      )
    }
    if (length(x = features) > 1) {
      features <- features[1]
      type <- ifelse(test = is.null(x = group.by), yes = 'feature', no = 'grouping')
      warning(
        "'do.hover' requires only one ",
        type,
        ", using ",
        features,
        call. = FALSE,
        immediate. = TRUE
      )
    }
    if (facet.highlight) {
      warning(
        "'do.hover' requires no faceting highlighted cells",
        call. = FALSE,
        immediate. = TRUE
      )
      facet.highlight <- FALSE
    }
  }
  if (facet.highlight) {
    if (length(x = images) > 1) {
      images <- images[1]
      warning(
        "Faceting the highlight only works with a single image, using image ",
        images,
        call. = FALSE,
        immediate. = TRUE
      )
    }
    ncols <- length(x = cells.highlight)
  } else {
    ncols <- length(x = images)
  }
  plots <- vector(
    mode = "list",
    length = length(x = features) * ncols
  )
  for (i in 1:ncols) {
    plot.idx <- i
    image.idx <- ifelse(test = facet.highlight, yes = 1, no = i)
    image.use <- object[[images[[image.idx]]]]
    coordinates <- GetTissueCoordinates(object = image.use)
    highlight.use <- if (facet.highlight) {
      cells.highlight[i]
    } else {
      cells.highlight
    }
    for (j in 1:length(x = features)) {
      cols.unset <- is.factor(x = data[, features[j]]) && is.null(x = cols)
      if (cols.unset) {
        cols <- hue_pal()(n = length(x = levels(x = data[, features[j]])))
        names(x = cols) <- levels(x = data[, features[j]])
      }
      plot <- SingleSpatialPlot(
        data = cbind(
          coordinates,
          data[rownames(x = coordinates), features[j], drop = FALSE]
        ),
        image = image.use,
        image.alpha = image.alpha,
        col.by = features[j],
        cols = cols,
        alpha.by = if (is.null(x = group.by)) {
          features[j]
        } else {
          NULL
        },
        pt.alpha = if (!is.null(x = group.by)) {
          alpha[j]
        } else {
          NULL
        },
        geom = if (inherits(x = image.use, what = "STARmap")) {
          'poly'
        } else {
          'spatial'
        },
        cells.highlight = highlight.use,
        cols.highlight = cols.highlight,
        pt.size.factor = pt.size.factor,
        stroke = stroke,
        crop = crop
      )
      if (is.null(x = group.by)) {
        plot <- plot +
          scale_fill_gradientn(
            name = features[j],
            colours = SpatialColors(n = 100)
          ) +
          theme(legend.position = 'top') +
          scale_alpha(range = alpha) +
          guides(alpha = "none")
      } else if (label) {
        plot <- LabelClusters(
          plot = plot,
          id = ifelse(
            test = is.null(x = cells.highlight),
            yes = features[j],
            no = 'highlight'
          ),
          geom = if (inherits(x = image.use, what = "STARmap")) {
            'GeomPolygon'
          } else {
            'GeomSpatial'
          },
          repel = repel,
          size = label.size,
          color = label.color,
          box = label.box,
          position = "nearest"
        )
      }
      if (j == 1 && length(x = images) > 1 && !facet.highlight) {
        plot <- plot +
          ggtitle(label = images[[image.idx]]) +
          theme(plot.title = element_text(hjust = 0.5))
      }
      if (facet.highlight) {
        plot <- plot +
          ggtitle(label = names(x = cells.highlight)[i]) +
          theme(plot.title = element_text(hjust = 0.5)) +
          NoLegend()
      }
      plots[[plot.idx]] <- plot
      plot.idx <- plot.idx + ncols
      if (cols.unset) {
        cols <- NULL
      }
    }
  }
  # if (do.identify) {
  #   return(CellSelector(
  #     plot = plot,
  #     object = identify.ident %iff% object,
  #     ident = identify.ident
  #   ))
  # } else if (do.hover) {
  #   return(HoverLocator(
  #     plot = plots[[1]],
  #     information = information %||% data[, features, drop = FALSE],
  #     axes = FALSE,
  #     # cols = c('size' = 'point.size.factor', 'colour' = 'fill'),
  #     images = GetImage(object = object, mode = 'plotly', image = images)
  #   ))
  # }
  if (combine) {
    if (!is.null(x = ncol)) {
      return(wrap_plots(plots = plots, ncol = ncol))
    }
    if (length(x = images) > 1) {
      return(wrap_plots(plots = plots, ncol = length(x = images)))
    }
    return(wrap_plots(plots = plots))
  }
  return(plots)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Other plotting functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
#' @concept visualization
#'
#' @author Robert A. Amezquita, \email{robert.amezquita@fredhutch.org}
#' @seealso \code{\link{CalculateBarcodeInflections}} \code{\link{SubsetByBarcodeInflections}}
#'
#' @examples
#' data("pbmc_small")
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
#' @param features Input vector of features, or named list of feature vectors
#' if feature-grouped panels are desired (replicates the functionality of the
#' old SplitDotPlotGG)
#' @param cols Colors to plot: the name of a palette from
#' \code{RColorBrewer::brewer.pal.info}, a pair of colors defining a gradient,
#' or 3+ colors defining multiple gradients (if split.by is set)
#' @param col.min Minimum scaled average expression threshold (everything
#' smaller will be set to this)
#' @param col.max Maximum scaled average expression threshold (everything larger
#' will be set to this)
#' @param dot.min The fraction of cells at which to draw the smallest dot
#' (default is 0). All cell groups with less than this expressing the given
#' gene will have no dot drawn.
#' @param dot.scale Scale the size of the points, similar to cex
#' @param idents Identity classes to include in plot (default is all)
#' @param group.by Factor to group the cells by
#' @param split.by Factor to split the groups by (replicates the functionality
#' of the old SplitDotPlotGG);
#' see \code{\link{FetchData}} for more details
#' @param cluster.idents Whether to order identities by hierarchical clusters
#' based on given features, default is FALSE
#' @param scale Determine whether the data is scaled, TRUE for default
#' @param scale.by Scale the size of the points by 'size' or by 'radius'
#' @param scale.min Set lower limit for scaling, use NA for default
#' @param scale.max Set upper limit for scaling, use NA for default
#'
#' @return A ggplot object
#'
#' @importFrom grDevices colorRampPalette
#' @importFrom cowplot theme_cowplot
#' @importFrom ggplot2 ggplot aes_string geom_point scale_size scale_radius
#' theme element_blank labs scale_color_identity scale_color_distiller
#' scale_color_gradient guides guide_legend guide_colorbar
#' facet_grid unit
#' @importFrom scattermore geom_scattermore
#' @importFrom stats dist hclust
#' @importFrom RColorBrewer brewer.pal.info
#'
#' @export
#' @concept visualization
#'
#' @aliases SplitDotPlotGG
#' @seealso \code{RColorBrewer::brewer.pal.info}
#'
#' @examples
#' data("pbmc_small")
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
  idents = NULL,
  group.by = NULL,
  split.by = NULL,
  cluster.idents = FALSE,
  scale = TRUE,
  scale.by = 'radius',
  scale.min = NA,
  scale.max = NA
) {
  assay <- assay %||% DefaultAssay(object = object)
  DefaultAssay(object = object) <- assay
  split.colors <- !is.null(x = split.by) && !any(cols %in% rownames(x = brewer.pal.info))
  scale.func <- switch(
    EXPR = scale.by,
    'size' = scale_size,
    'radius' = scale_radius,
    stop("'scale.by' must be either 'size' or 'radius'")
  )
  feature.groups <- NULL
  if (is.list(features) | any(!is.na(names(features)))) {
    feature.groups <- unlist(x = sapply(
      X = 1:length(features),
      FUN = function(x) {
        return(rep(x = names(x = features)[x], each = length(features[[x]])))
      }
    ))
    if (any(is.na(x = feature.groups))) {
      warning(
        "Some feature groups are unnamed.",
        call. = FALSE,
        immediate. = TRUE
      )
    }
    features <- unlist(x = features)
    names(x = feature.groups) <- features
  }
  cells <- unlist(x = CellsByIdentities(object = object, idents = idents))

  data.features <- FetchData(object = object, vars = features, cells = cells)
  data.features$id <- if (is.null(x = group.by)) {
    Idents(object = object)[cells, drop = TRUE]
  } else {
    object[[group.by, drop = TRUE]][cells, drop = TRUE]
  }
  if (!is.factor(x = data.features$id)) {
    data.features$id <- factor(x = data.features$id)
  }
  id.levels <- levels(x = data.features$id)
  data.features$id <- as.vector(x = data.features$id)
  if (!is.null(x = split.by)) {
    splits <- object[[split.by, drop = TRUE]][cells, drop = TRUE]
    if (split.colors) {
      if (length(x = unique(x = splits)) > length(x = cols)) {
        stop("Not enough colors for the number of groups")
      }
      cols <- cols[1:length(x = unique(x = splits))]
      names(x = cols) <- unique(x = splits)
    }
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
  if (cluster.idents) {
    mat <- do.call(
      what = rbind,
      args = lapply(X = data.plot, FUN = unlist)
    )
    mat <- scale(x = mat)
    id.levels <- id.levels[hclust(d = dist(x = mat))$order]
  }
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
  ngroup <- length(x = levels(x = data.plot$id))
  if (ngroup == 1) {
    scale <- FALSE
    warning(
      "Only one identity present, the expression values will be not scaled",
      call. = FALSE,
      immediate. = TRUE
    )
  } else if (ngroup < 5 & scale) {
    warning(
      "Scaling data with a low number of groups may produce misleading results",
      call. = FALSE,
      immediate. = TRUE
    )
  }
  avg.exp.scaled <- sapply(
    X = unique(x = data.plot$features.plot),
    FUN = function(x) {
      data.use <- data.plot[data.plot$features.plot == x, 'avg.exp']
      if (scale) {
        data.use <- scale(x = data.use)
        data.use <- MinMax(data = data.use, min = col.min, max = col.max)
      } else {
        data.use <- log1p(x = data.use)
      }
      return(data.use)
    }
  )
  avg.exp.scaled <- as.vector(x = t(x = avg.exp.scaled))
  if (split.colors) {
    avg.exp.scaled <- as.numeric(x = cut(x = avg.exp.scaled, breaks = 20))
  }
  data.plot$avg.exp.scaled <- avg.exp.scaled
  data.plot$features.plot <- factor(
    x = data.plot$features.plot,
    levels = features
  )
  data.plot$pct.exp[data.plot$pct.exp < dot.min] <- NA
  data.plot$pct.exp <- data.plot$pct.exp * 100
  if (split.colors) {
    splits.use <- vapply(
      X = as.character(x = data.plot$id),
      FUN = gsub,
      FUN.VALUE = character(length = 1L),
      pattern =  paste0(
        '^((',
        paste(sort(x = levels(x = object), decreasing = TRUE), collapse = '|'),
        ')_)'
      ),
      replacement = '',
      USE.NAMES = FALSE
    )
    data.plot$colors <- mapply(
      FUN = function(color, value) {
        return(colorRampPalette(colors = c('grey', color))(20)[value])
      },
      color = cols[splits.use],
      value = avg.exp.scaled
    )
  }
  color.by <- ifelse(test = split.colors, yes = 'colors', no = 'avg.exp.scaled')
  if (!is.na(x = scale.min)) {
    data.plot[data.plot$pct.exp < scale.min, 'pct.exp'] <- scale.min
  }
  if (!is.na(x = scale.max)) {
    data.plot[data.plot$pct.exp > scale.max, 'pct.exp'] <- scale.max
  }
  if (!is.null(x = feature.groups)) {
    data.plot$feature.groups <- factor(
      x = feature.groups[data.plot$features.plot],
      levels = unique(x = feature.groups)
    )
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
  if (!is.null(x = feature.groups)) {
    plot <- plot + facet_grid(
      facets = ~feature.groups,
      scales = "free_x",
      space = "free_x",
      switch = "y"
    ) + theme(
      panel.spacing = unit(x = 1, units = "lines"),
      strip.background = element_blank()
    )
  }
  if (split.colors) {
    plot <- plot + scale_color_identity()
  } else if (length(x = cols) == 1) {
    plot <- plot + scale_color_distiller(palette = cols)
  } else {
    plot <- plot + scale_color_gradient(low = cols[1], high = cols[2])
  }
  if (!split.colors) {
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
#' @concept visualization
#'
#' @examples
#' data("pbmc_small")
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

#' Boxplot of correlation of a variable (e.g. number of UMIs) with expression
#' data
#'
#' @param object Seurat object
#' @param assay Assay where the feature grouping info and correlations are
#' stored
#' @param feature.group Name of the column in meta.features where the feature
#' grouping info is stored
#' @param cor Name of the column in meta.features where correlation info is
#' stored
#'
#' @return Returns a ggplot boxplot of correlations split by group
#'
#' @importFrom ggplot2 geom_boxplot scale_fill_manual geom_hline
#' @importFrom cowplot theme_cowplot
#' @importFrom scales brewer_pal
#' @importFrom stats complete.cases
#'
#' @export
#' @concept visualization
#'
GroupCorrelationPlot <- function(
  object,
  assay = NULL,
  feature.group = "feature.grp",
  cor = "nCount_RNA_cor"
) {
  assay <- assay %||% DefaultAssay(object = object)
  data <- object[[assay]][c(feature.group, cor)]
  data <- data[complete.cases(data), ]
  colnames(x = data) <- c('grp', 'cor')
  data$grp <- as.character(data$grp)
  plot <- ggplot(data = data, aes_string(x = "grp", y = "cor", fill = "grp")) +
    geom_boxplot() +
    theme_cowplot() +
    scale_fill_manual(values = rev(x = brewer_pal(palette = 'YlOrRd')(n = 7))) +
    ylab(paste(
      "Correlation with",
      gsub(x = cor, pattern = "_cor", replacement = "")
    )) +
    geom_hline(yintercept = 0) +
    NoLegend() +
    theme(
      axis.line.x = element_blank(),
      axis.title.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.x = element_blank()
    )
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
#' @param cols Vector of colors, each color corresponds to an individual PC. This may also be a single character
#' or numeric value corresponding to a palette as specified by \code{\link[RColorBrewer]{brewer.pal.info}}.
#' By default, ggplot2 assigns colors. We also include a number of palettes from the pals package.
#' See \code{\link{DiscretePalette}} for details.
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
#' @importFrom scales hue_pal
#' @importFrom ggplot2 ggplot aes_string stat_qq labs xlim ylim
#' coord_flip geom_abline guides guide_legend
#' @importFrom cowplot theme_cowplot
#'
#' @export
#' @concept visualization
#'
#' @examples
#' data("pbmc_small")
#' JackStrawPlot(object = pbmc_small)
#'
JackStrawPlot <- function(
  object,
  dims = 1:5,
  cols = NULL,
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
  if (is.null(x = cols)) {
    cols <- hue_pal()(length(x = dims))
  }
  if (length(x = cols) < length(x = dims)) {
    stop("Not enough colors for the number of dims selected")
  }
  gp <- ggplot(data = data.plot, mapping = aes_string(sample = 'Value', color = 'PC.Score')) +
    stat_qq(distribution = qunif) +
    labs(x = "Theoretical [runif(1000)]", y = "Empirical") +
    scale_color_manual(values = cols) +
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
#' @param direction A character string specifying the direction of the tree (default is downwards)
#' Possible options: "rightwards", "leftwards", "upwards", and "downwards".
#' @param \dots Additional arguments to
#' \code{\link[ape:plot.phylo]{ape::plot.phylo}}
#'
#' @return Plots dendogram (must be precomputed using BuildClusterTree), returns no value
#'
#' @export
#' @concept visualization
#'
#' @examples
#' if (requireNamespace("ape", quietly = TRUE)) {
#'   data("pbmc_small")
#'   pbmc_small <- BuildClusterTree(object = pbmc_small)
#'   PlotClusterTree(object = pbmc_small)
#' }
PlotClusterTree <- function(object, direction = "downwards", ...) {
  if (!PackageCheck('ape', error = FALSE)) {
    stop(cluster.ape, call. = FALSE)
  }
  if (is.null(x = Tool(object = object, slot = "BuildClusterTree"))) {
    stop("Phylogenetic tree does not exist, build using BuildClusterTree")
  }
  data.tree <- Tool(object = object, slot = "BuildClusterTree")
  ape::plot.phylo(x = data.tree, direction = direction, ...)
  ape::nodelabels()
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
#' @param combine Combine plots into a single \code{patchwork}
#' ggplot object. If \code{FALSE}, return a list of ggplot objects
#'
#' @return A \code{patchwork} ggplot object if
#' \code{combine = TRUE}; otherwise, a list of ggplot objects
#'
#' @importFrom patchwork wrap_plots
#' @importFrom cowplot theme_cowplot
#' @importFrom ggplot2 ggplot aes_string geom_point labs
#' @export
#' @concept visualization
#'
#' @examples
#' data("pbmc_small")
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
  if (is.null(x = ncol)) {
    ncol <- 2
    if (length(x = dims) == 1) {
      ncol <- 1
    }
    if (length(x = dims) > 6) {
      ncol <- 3
    }
    if (length(x = dims) > 9) {
      ncol <- 4
    }
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
#' @concept visualization
#'
#' @examples
#' \dontrun{
#' data("pbmc_small")
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

#' Automagically calculate a point size for ggplot2-based scatter plots
#'
#' It happens to look good
#'
#' @param data A data frame being passed to ggplot2
#' @param raster If TRUE, point size is set to 1
#'
#' @return The "optimal" point size for visualizing these data
#'
#' @export
#' @concept visualization
#'
#' @examples
#' df <- data.frame(x = rnorm(n = 10000), y = runif(n = 10000))
#' AutoPointSize(data = df)
#'
AutoPointSize <- function(data, raster = NULL) {
  return(ifelse(
    test = isTRUE(x = raster),
    yes = 1,
    no = min(1583 / nrow(x = data), 1)
  ))
}

#' Determine text color based on background color
#'
#' @param background A vector of background colors; supports R color names and
#' hexadecimal codes
#' @param threshold Intensity threshold for light/dark cutoff; intensities
#' greater than \code{theshold} yield \code{dark}, others yield \code{light}
#' @param w3c Use \href{http://www.w3.org/TR/WCAG20/}{W3C} formula for calculating
#' background text color; ignores \code{threshold}
#' @param dark Color for dark text
#' @param light Color for light text
#'
#' @return A named vector of either \code{dark} or \code{light}, depending on
#' \code{background}; names of vector are \code{background}
#'
#' @export
#' @concept visualization
#'
#' @source \url{https://stackoverflow.com/questions/3942878/how-to-decide-font-color-in-white-or-black-depending-on-background-color}
#'
#' @examples
#' BGTextColor(background = c('black', 'white', '#E76BF3'))
#'
BGTextColor <- function(
  background,
  threshold = 186,
  w3c = FALSE,
  dark = 'black',
  light = 'white'
) {
  if (w3c) {
    luminance <- Luminance(color = background)
    threshold <- 179
    return(ifelse(
      test = luminance > sqrt(x = 1.05 * 0.05) - 0.05,
      yes = dark,
      no = light
    ))
  }
  return(ifelse(
    test = Intensity(color = background) > threshold,
    yes = dark,
    no = light
  ))
}

#' @export
#' @concept visualization
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

#' @export
#' @concept visualization
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

#' Cell Selector
#'
#' Select points on a scatterplot and get information about them
#'
#' @param plot A ggplot2 plot
#' @param object An optional Seurat object; if passes, will return an object
#' with the identities of selected cells set to \code{ident}
#' @param ident An optional new identity class to assign the selected cells
#' @param ... Ignored
#'
#' @return If \code{object} is \code{NULL}, the names of the points selected;
#' otherwise, a Seurat object with the selected cells identity classes set to
#' \code{ident}
#'
#' @importFrom miniUI miniPage gadgetTitleBar miniTitleBarButton
#' miniContentPanel
#' @importFrom shiny fillRow plotOutput brushOpts reactiveValues observeEvent
#' stopApp brushedPoints renderPlot runGadget
#'
#' @export
#' @concept visualization
#'
#' @seealso \code{\link{DimPlot}} \code{\link{FeaturePlot}}
#'
#' @examples
#' \dontrun{
#' data("pbmc_small")
#' plot <- DimPlot(object = pbmc_small)
#' # Follow instructions in the terminal to select points
#' cells.located <- CellSelector(plot = plot)
#' cells.located
#' # Automatically set the identity class of selected cells and return a new Seurat object
#' pbmc_small <- CellSelector(plot = plot, object = pbmc_small, ident = 'SelectedCells')
#' }
#'
CellSelector <- function(plot, object = NULL, ident = 'SelectedCells', ...) {
  # Set up the gadget UI
  ui <- miniPage(
    gadgetTitleBar(
      title = "Cell Selector",
      left = miniTitleBarButton(inputId = "reset", label = "Reset")
    ),
    miniContentPanel(
      fillRow(
        plotOutput(
          outputId = "plot",
          height = '100%',
          brush = brushOpts(
            id = 'brush',
            delay = 100,
            delayType = 'debounce',
            clip = TRUE,
            resetOnNew = FALSE
          )
        )
      ),
    )
  )
  # Get some plot information
  if (inherits(x = plot, what = 'patchwork')) {
    if (length(x = plot$patches$plots)) {
      warning(
        "Multiple plots passed, using last plot",
        call. = FALSE,
        immediate. = TRUE
      )
    }
    class(x = plot) <- grep(
      pattern = 'patchwork',
      x = class(x = plot),
      value = TRUE,
      invert = TRUE
    )
  }
  xy.aes <- GetXYAesthetics(plot = plot)
  dark.theme <- !is.null(x = plot$theme$plot.background$fill) &&
    plot$theme$plot.background$fill == 'black'
  plot.data <- GGpointToBase(plot = plot, do.plot = FALSE)
  plot.data$selected_ <- FALSE
  rownames(x = plot.data) <- rownames(x = plot$data)
  colnames(x = plot.data) <- gsub(
    pattern = '-',
    replacement = '.',
    x = colnames(x = plot.data)
  )
  # Server function
  server <- function(input, output, session) {
    plot.env <- reactiveValues(data = plot.data)
    # Event handlers
    observeEvent(
      eventExpr = input$done,
      handlerExpr = {
        PlotBuild(data = plot.env$data, dark.theme = dark.theme)
        selected <- rownames(x = plot.data)[plot.env$data$selected_]
        if (inherits(x = object, what = 'Seurat')) {
          if (!all(selected %in% Cells(x = object))) {
            stop("Cannot find the selected cells in the Seurat object, please be sure you pass the same object used to generate the plot")
          }
          Idents(object = object, cells = selected) <- ident
          selected <- object
        }
        stopApp(returnValue = selected)
      }
    )
    observeEvent(
      eventExpr = input$reset,
      handlerExpr = {
        plot.env$data <- plot.data
        session$resetBrush(brushId = 'brush')
      }
    )
    observeEvent(
      eventExpr = input$brush,
      handlerExpr = {
        plot.env$data <- brushedPoints(
          df = plot.data,
          brush = input$brush,
          xvar = xy.aes$x,
          yvar = xy.aes$y,
          allRows = TRUE
        )
        plot.env$data$color <- ifelse(
          test = plot.env$data$selected_,
          yes = '#DE2D26',
          no = '#C3C3C3'
        )
      }
    )
    # Render the plot
    output$plot <- renderPlot(expr = PlotBuild(
      data = plot.env$data,
      dark.theme = dark.theme
    ))
  }
  return(runGadget(app = ui, server = server))
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
#' @concept visualization
#'
#' @examples
#' \dontrun{
#' data("pbmc_small")
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
#' @concept visualization
#'
#' @examples
#' data("pbmc_small")
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
#' @concept visualization
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

#' Discrete colour palettes from pals
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
#' "alphabet", "alphabet2", "glasbey", "polychrome", "stepped", and "parade".
#' Can be omitted and the function will use the one based on the requested n.
#' @param shuffle Shuffle the colors in the selected palette.
#'
#' @return A vector of colors
#'
#' @details
#' Taken from the pals package (Licence: GPL-3).
#' \url{https://cran.r-project.org/package=pals}
#' Credit: Kevin Wright
#'
#' @export
#' @concept visualization
#'
DiscretePalette <- function(n, palette = NULL, shuffle = FALSE) {
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
    ),
    parade = c(
      '#ff6969', '#9b37ff', '#cd3737', '#69cdff', '#ffff69', '#69cdcd',
      '#9b379b', '#3737cd', '#ffff9b', '#cdff69', '#ff9b37', '#37ffff',
      '#9b69ff', '#37cd69', '#ff3769', '#ff3737', '#37ff9b', '#cdcd37',
      '#3769cd', '#37cdff', '#9b3737', '#ff699b', '#9b9bff', '#cd9b37',
      '#69ff37', '#cd3769', '#cd69cd', '#cd6937', '#3737ff', '#cdcd69',
      '#ff9b69', '#cd37cd', '#9bff37', '#cd379b', '#cd6969', '#69ff9b',
      '#ff379b', '#9bff9b', '#6937ff', '#69cd37', '#cdff37', '#9bff69',
      '#9b37cd', '#ff37ff', '#ff37cd', '#ffff37', '#37cd9b', '#379bff',
      '#ffcd37', '#379b37', '#ff9bff', '#379b9b', '#69ffcd', '#379bcd',
      '#ff69ff', '#ff9b9b', '#37ff69', '#ff6937', '#6969ff', '#699bff',
      '#ffcd69', '#69ffff', '#37ff37', '#6937cd', '#37cd37', '#3769ff',
      '#cd69ff', '#6969cd', '#9bcd37', '#69ff69', '#37cdcd', '#cd37ff',
      '#37379b', '#37ffcd', '#69cd69', '#ff69cd', '#9bffff', '#9b9b37'
    )
  )
  if (is.null(x = n)) {
    return(names(x = palettes))
  }
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
  if (isTRUE(shuffle)) {
    palette.vec <- sample(palette.vec)
  }
  palette <- palette.vec[seq_len(length.out = n)]
  return(palette)
}

#' @rdname CellSelector
#' @export
#' @concept visualization
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
#' @param axes Display or hide x- and y-axes
#' @param ... Extra parameters to be passed to \code{\link[plotly]{layout}}
#'
#' @importFrom ggplot2 ggplot_build
#' @importFrom plotly plot_ly layout add_annotations
#' @export
#' @concept visualization
#'
#' @seealso \code{\link[plotly]{layout}} \code{\link[ggplot2]{ggplot_build}}
#' \code{\link{DimPlot}} \code{\link{FeaturePlot}}
#'
#' @examples
#' \dontrun{
#' data("pbmc_small")
#' plot <- DimPlot(object = pbmc_small)
#' HoverLocator(plot = plot, information = FetchData(object = pbmc_small, vars = 'percent.mito'))
#' }
#'
HoverLocator <- function(
  plot,
  information = NULL,
  axes = TRUE,
  dark.theme = FALSE,
  ...
) {
  # Use GGpointToBase because we already have ggplot objects
  # with colors (which are annoying in plotly)
  plot.build <- suppressWarnings(expr = GGpointToPlotlyBuild(
    plot = plot,
    information = information,
    ...
  ))
  data <- ggplot_build(plot = plot)$plot$data
  # Set up axis labels here
  # Also, a bunch of stuff to get axis lines done properly
  if (axes) {
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
  } else {
    xaxis <- yaxis <- list(visible = FALSE)
  }
  # Check for dark theme
  if (dark.theme) {
    title <- list(color = 'white')
    xaxis <- c(xaxis, color = 'white')
    yaxis <- c(yaxis, color = 'white')
    plotbg <- 'black'
  } else {
    title = list(color = 'black')
    plotbg = 'white'
  }
  # The `~' means pull from the data passed (this is why we reset the names)
  # Use I() to get plotly to accept the colors from the data as is
  # Set hoverinfo to 'text' to override the default hover information
  # rather than append to it
  p <- layout(
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
  # Add labels
  label.layer <- which(x = sapply(
    X = plot$layers,
    FUN = function(x) {
      return(inherits(x = x$geom, what = c('GeomText', 'GeomTextRepel')))
    }
  ))
  if (length(x = label.layer) == 1) {
    p <- add_annotations(
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

#' Get the intensity and/or luminance of a color
#'
#' @param color A vector of colors
#'
#' @return A vector of intensities/luminances for each color
#'
#' @name contrast-theory
#' @rdname contrast-theory
#'
#' @importFrom grDevices col2rgb
#'
#' @export
#' @concept visualization
#'
#' @source \url{https://stackoverflow.com/questions/3942878/how-to-decide-font-color-in-white-or-black-depending-on-background-color}
#'
#' @examples
#' Intensity(color = c('black', 'white', '#E76BF3'))
#'
Intensity <- function(color) {
  intensities <- apply(
    X = col2rgb(col = color),
    MARGIN = 2,
    FUN = function(col) {
      col <- rbind(as.vector(x = col), c(0.299, 0.587, 0.114))
      return(sum(apply(X = col, MARGIN = 2, FUN = prod)))
    }
  )
  names(x = intensities) <- color
  return(intensities)
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
#' @param geom Name of geom to get X/Y aesthetic names for
#' @param box Use geom_label/geom_label_repel (includes a box around the text
#' labels)
#' @param position How to place the label if repel = FALSE. If "median", place
#' the label at the median position. If "nearest" place the label at the
#' position of the nearest data point to the median.
#' @param ... Extra parameters to \code{\link[ggrepel]{geom_text_repel}}, such as \code{size}
#'
#' @return A ggplot2-based scatter plot with cluster labels
#'
#' @importFrom stats median na.omit
#' @importFrom ggrepel geom_text_repel geom_label_repel
#' @importFrom ggplot2 aes_string geom_text geom_label layer_scales
#' @importFrom RANN nn2
#'
#' @export
#' @concept visualization
#'
#' @seealso \code{\link[ggrepel]{geom_text_repel}} \code{\link[ggplot2]{geom_text}}
#'
#' @examples
#' data("pbmc_small")
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
  box = FALSE,
  geom = 'GeomPoint',
  position = "median",
  ...
) {
  xynames <- unlist(x = GetXYAesthetics(plot = plot, geom = geom), use.names = TRUE)
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
  pb <- ggplot_build(plot = plot)
  if (geom == 'GeomSpatial') {
    xrange.save <- layer_scales(plot = plot)$x$range$range
    yrange.save <- layer_scales(plot = plot)$y$range$range
    data[, xynames["y"]] = max(data[, xynames["y"]]) - data[, xynames["y"]] + min(data[, xynames["y"]])
    if (!pb$plot$plot_env$crop) {
      y.transform <- c(0, nrow(x = pb$plot$plot_env$image)) - pb$layout$panel_params[[1]]$y.range
      data[, xynames["y"]] <- data[, xynames["y"]] + sum(y.transform)
    }
  }
  data <- cbind(data, color = pb$data[[1]][[1]])
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
      data.medians$color <- data.use$color[1]
      return(data.medians)
    }
  )
  if (position == "nearest") {
    labels.loc <- lapply(X = labels.loc, FUN = function(x) {
      group.data <- data[as.character(x = data[, id]) == as.character(x[3]), ]
      nearest.point <- nn2(data = group.data[, 1:2], query = as.matrix(x = x[c(1,2)]), k = 1)$nn.idx
      x[1:2] <- group.data[nearest.point, 1:2]
      return(x)
    })
  }
  labels.loc <- do.call(what = 'rbind', args = labels.loc)
  labels.loc[, id] <- factor(x = labels.loc[, id], levels = levels(data[, id]))
  labels <- labels %||% groups
  if (length(x = unique(x = labels.loc[, id])) != length(x = labels)) {
    stop("Length of labels (", length(x = labels),  ") must be equal to the number of clusters being labeled (", length(x = labels.loc), ").")
  }
  names(x = labels) <- groups
  for (group in groups) {
    labels.loc[labels.loc[, id] == group, id] <- labels[group]
  }
  if (box) {
    geom.use <- ifelse(test = repel, yes = geom_label_repel, no = geom_label)
    plot <- plot + geom.use(
      data = labels.loc,
      mapping = aes_string(x = xynames['x'], y = xynames['y'], label = id, fill = id),
      show.legend = FALSE,
      ...
    ) + scale_fill_manual(values = labels.loc$color[order(labels.loc[, id])])
  } else {
    geom.use <- ifelse(test = repel, yes = geom_text_repel, no = geom_text)
    plot <- plot + geom.use(
      data = labels.loc,
      mapping = aes_string(x = xynames['x'], y = xynames['y'], label = id),
      show.legend = FALSE,
      ...
    )
  }
  # restore old axis ranges
  if (geom == 'GeomSpatial') {
    plot <- suppressMessages(expr = plot + coord_fixed(xlim = xrange.save, ylim = yrange.save))
  }
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
#' @concept visualization
#'
#' @aliases Labeler
#' @seealso \code{\link[ggplot2]{geom_text}}
#'
#' @examples
#' data("pbmc_small")
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

#' @name contrast-theory
#' @rdname contrast-theory
#'
#' @importFrom grDevices col2rgb
#'
#' @export
#' @concept visualization
#'
#' @examples
#' Luminance(color = c('black', 'white', '#E76BF3'))
#'
Luminance <- function(color) {
  luminance <- apply(
    X = col2rgb(col = color),
    MARGIN = 2,
    function(col) {
      col <- as.vector(x = col) / 255
      col <- sapply(
        X = col,
        FUN = function(x) {
          return(ifelse(
            test = x <= 0.03928,
            yes = x / 12.92,
            no = ((x + 0.055) / 1.055) ^ 2.4
          ))
        }
      )
      col <- rbind(col, c(0.2126, 0.7152, 0.0722))
      return(sum(apply(X = col, MARGIN = 2, FUN = prod)))
    }
  )
  names(x = luminance) <- color
  return(luminance)
}

#' @export
#' @concept visualization
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
#' @concept visualization
#'
#' @rdname SeuratTheme
#' @seealso \code{\link[ggplot2]{theme}}
#' @aliases SeuratTheme
#'
SeuratTheme <- function() {
  return(DarkTheme() + NoLegend() + NoGrid() + SeuratAxes())
}

#' @importFrom ggplot2 theme element_text
#'
#' @rdname SeuratTheme
#' @export
#' @concept visualization
#'
#' @aliases CenterTitle
#'
CenterTitle <- function(...) {
  return(theme(plot.title = element_text(hjust = 0.5), validate = TRUE, ...))
}

#' @importFrom ggplot2 theme element_rect element_text element_line margin
#' @export
#' @concept visualization
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

#' @param x.text,y.text X and Y axis text sizes
#' @param x.title,y.title X and Y axis title sizes
#' @param main Plot title size
#'
#' @importFrom ggplot2 theme element_text
#' @export
#' @concept visualization
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

#' @param keep.text Keep axis text
#' @param keep.ticks Keep axis ticks
#'
#' @importFrom ggplot2 theme element_blank
#' @export
#' @concept visualization
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

#' @importFrom ggplot2 theme
#' @export
#' @concept visualization
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

#' @importFrom ggplot2 theme element_blank
#' @export
#' @concept visualization
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

#' @importFrom ggplot2 theme element_text
#' @export
#' @concept visualization
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

#' @export
#' @concept visualization
#'
#' @rdname SeuratTheme
#' @aliases SpatialTheme
#'
SpatialTheme <- function(...) {
  return(DarkTheme() + NoAxes() + NoGrid() + NoLegend(...))
}

#' @param position A position to restore the legend to
#'
#' @importFrom ggplot2 theme
#' @export
#' @concept visualization
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

#' @importFrom ggplot2 theme element_text
#' @export
#' @concept visualization
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

#' @importFrom ggplot2 theme element_text
#' @export
#' @concept visualization
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

#' @importFrom ggplot2 theme element_rect
#' @export
#' @concept visualization
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
# Fortify Methods
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Prepare Coordinates for Spatial Plots
#'
#' @inheritParams SeuratObject::GetTissueCoordinates
#' @param model A \code{\linkS4class{Segmentation}},
#' \code{\linkS4class{Centroids}},
#' or \code{\linkS4class{Molecules}} object
#' @param data Extra data to be used for annotating the cell segmentations; the
#' easiest way to pass data is a one-column
#' \code{\link[base:data.frame]{data frame}} with the values to color by and
#' the cell names are rownames
#' @param ... Arguments passed to other methods
#'
#' @name fortify-Spatial
#' @rdname fortify-Spatial
#'
#' @importFrom SeuratObject GetTissueCoordinates
#'
#' @keywords internal
#'
#' @method fortify Centroids
#' @export
#'
#' @aliases fortify
#'
fortify.Centroids <- function(model, data, ...) {
  df <- GetTissueCoordinates(object = model, full = FALSE)
  if (missing(x = data)) {
    data <- NULL
  }
  data <- .PrepImageData(data = data, cells = lengths(x = model), ...)
  df <- cbind(df, data)
  return(df)
}

#' @rdname fortify-Spatial
#' @method fortify Molecules
#'
#' @importFrom SeuratObject FetchData
#'
#' @export
#'
fortify.Molecules <- function(
  model,
  data,
  nmols = NULL,
  seed = NA_integer_,
  ...
) {
  return(FetchData(object = model, vars = data, nmols = nmols, seed = seed, ...))
}

#' @rdname fortify-Spatial
#' @method fortify Segmentation
#' @export
#'
fortify.Segmentation <- function(model, data, ...) {
  df <- GetTissueCoordinates(object = model, full = TRUE)
  if (missing(x = data)) {
    data <- NULL
  }
  data <- .PrepImageData(data = data, cells = lengths(x = model), ...)
  df <- cbind(df, data)
  return(df)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Internal
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @importFrom SeuratObject Features Key Keys Molecules
#'
.MolsByFOV <- function(object, fov, molecules) {
  keys <- Key(object = object)[fov]
  keyed.mols <- sapply(
    X = names(x = keys),
    FUN = function(img) {
      if (is.null(x = Molecules(object = object[[img]]))) {
        return(NULL)
      }
      key <- keys[img]
      mols <- grep(pattern = paste0('^', key), x = molecules, value = TRUE)
      names(x = mols) <- mols
      mols <- gsub(pattern = paste0('^', key), replacement = '', x = mols)
      keyed <- sapply(
        X = SeuratObject::Keys(object = object[[img]]),
        FUN = function(x) {
          return(grep(pattern = paste0('^', x), x = mols, value = TRUE))
        }
      )
      keyed <- unlist(x = keyed)
      names(x = keyed) <- gsub(
        pattern = '^.*\\.',
        replacement = '',
        x = names(x = keyed)
      )
      missing <- mols[!mols %in% keyed]
      missing <- missing[missing %in% Features(x = object[[img]])]
      if (length(x = missing)) {
        # TODO: replace with default molecules
        default <- Molecules(object = object[[img]])[1L]
        mn <- names(x = missing)
        missing <- paste0(
          SeuratObject::Key(object = object[[img]][[default]]),
          missing
        )
        names(x = missing) <- mn
      }
      return(c(missing, keyed))
    },
    simplify = FALSE,
    USE.NAMES = TRUE
  )
  found <- names(x = unlist(x = keyed.mols))
  found <- gsub(pattern = '^.*\\.', replacement = '', x = found)
  missing <- setdiff(x = molecules, y = found)
  names(x = missing) <- missing
  for (img in fov) {
    imissing <- missing
    for (i in seq_along(along.with = imissing)) {
      for (lkey in Keys(object = object[[img]])) {
        imissing[[i]] <- gsub(
          pattern = paste0('^', lkey),
          replacement = '',
          x = imissing[[i]]
        )
      }
    }
    imissing <- names(x = imissing[imissing %in% Features(x = object[[img]])])
    keyed.mols[[img]] <- c(keyed.mols[[img]], imissing)
  }
  keyed.mols <- Filter(f = length, x = keyed.mols)
  keyed.mols <- sapply(X = keyed.mols, FUN = unname, simplify = FALSE)
  return(keyed.mols)
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
#' @importFrom grDevices col2rgb
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
  C0 <- as.vector(col2rgb(negative.color, alpha = TRUE))
  C1 <- as.vector(col2rgb(two.colors[1], alpha = TRUE))
  C2 <- as.vector(col2rgb(two.colors[2], alpha = TRUE))
  blend_alpha <- (C1[4] + C2[4])/2
  C0 <- C0[-4]
  C1 <- C1[-4]
  C2 <- C2[-4]
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
    alpha,
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
      red = C_blend[1],
      green = C_blend[2],
      blue = C_blend[3],
      alpha = alpha,
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
        alpha = blend_alpha,
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

# Plot feature expression by identity
#
# Basically combines the codebase for VlnPlot and RidgePlot
#
# @param object Seurat object
# @param type Plot type, choose from 'ridge', 'violin', or 'splitViolin'
# @param features Features to plot (gene expression, metrics, PC scores,
# anything that can be retreived by FetchData)
# @param idents Which classes to include in the plot (default is all)
# @param ncol Number of columns if multiple plots are displayed
# @param sort Sort identity classes (on the x-axis) by the average expression of the attribute being potted,
# or, if stack is True, sort both identity classes and features by hierarchical clustering
# @param y.max Maximum y axis value
# @param same.y.lims Set all the y-axis limits to the same values
# @param adjust Adjust parameter for geom_violin
# @param pt.size Point size for points
# @param alpha Alpha value for points
# @param cols Colors to use for plotting
# @param group.by Group (color) cells in different ways (for example, orig.ident)
# @param split.by A variable to split the plot by
# @param log plot Y axis on log scale
# @param slot Slot to pull expression data from (e.g. "counts" or "data")
# @param stack Horizontally stack plots for multiple feature
# @param combine Combine plots into a single \code{\link[patchwork]{patchwork}ed}
# ggplot object. If \code{FALSE}, return a list of ggplot objects
# @param fill.by Color violins/ridges based on either 'feature' or 'ident'
# @param flip flip plot orientation (identities on x-axis)
# @param add.noise determine if adding a small noise for plotting
# @param raster Convert points to raster format, default is \code{NULL} which
# automatically rasterizes if plotting more than 100,000 cells
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
  alpha = 1,
  group.by = NULL,
  split.by = NULL,
  log = FALSE,
  slot = deprecated(),
  layer = 'data',
  stack = FALSE,
  combine = TRUE,
  fill.by = NULL,
  flip = FALSE,
  add.noise = TRUE,
  raster = NULL
) {
  if (is_present(arg = slot)) {
    layer <- layer %||% slot
  }
  assay <- assay %||% DefaultAssay(object = object)
  DefaultAssay(object = object) <- assay
  cells <- Cells(x = object, assay = NULL)
  if (isTRUE(x = stack)) {
    if (!is.null(x = ncol)) {
      warning(
        "'ncol' is ignored with 'stack' is TRUE",
        call. = FALSE,
        immediate. = TRUE
      )
    }
    if (!is.null(x = y.max)) {
      warning(
        "'y.max' is ignored when 'stack' is TRUE",
        call. = FALSE,
        immediate. = TRUE
      )
    }
  } else {
    ncol <- ncol %||% ifelse(
      test = length(x = features) > 9,
      yes = 4,
      no = min(length(x = features), 3)
    )
  }
  if (!is.null(x = idents)) {
    cells <- intersect(
      x = names(x = Idents(object = object)[Idents(object = object) %in% idents]),
      y = cells
    )
  }
  data <- FetchData(object = object, vars = features, slot = layer, cells = cells)
  pt.size <- pt.size %||% AutoPointSize(data = object)
  features <- colnames(x = data)
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
    names(x = cols) <- levels(x = split)
    if ((length(x = cols) > 2) & (type == "splitViolin")) {
      warning("Split violin is only supported for <3 groups, using multi-violin.")
      type <- "violin"
    }
  }
  if (same.y.lims && is.null(x = y.max)) {
    y.max <- max(data)
  }
  if (isTRUE(x = stack)) {
    return(MultiExIPlot(
      type = type,
      data = data,
      idents = idents,
      split = split,
      sort = sort,
      same.y.lims = same.y.lims,
      adjust = adjust,
      cols = cols,
      pt.size = pt.size,
      log = log,
      fill.by = fill.by,
      add.noise = add.noise,
      flip = flip
    ))
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
        alpha = alpha,
        log = log,
        add.noise = add.noise,
        raster = raster
      ))
    }
  )
  label.fxn <- switch(
    EXPR = type,
    'violin' = if (stack) {
      xlab
    } else {
      ylab
    },
    "splitViolin" = if (stack) {
      xlab
    } else {
      ylab
    },
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

#' @importFrom RColorBrewer brewer.pal
#' @importFrom grDevices colorRampPalette
#'
#'
SpatialColors <- colorRampPalette(colors = rev(x = brewer.pal(n = 11, name = "Spectral")))

# Feature plot palettes
#
FeaturePalettes <- list(
  'Spatial' = SpatialColors(n = 100),
  'Seurat' = c('lightgrey', 'blue')
)

# Splits features into groups based on log expression levels
#
# @param object Seurat object
# @param assay Assay for expression data
# @param min.cells Only compute for features in at least this many cells
# @param ngroups Number of groups to split into
#
# @return A Seurat object with the feature group stored as a factor in
# metafeatures
#
#' @importFrom Matrix rowMeans rowSums
#
GetFeatureGroups <- function(object, assay, min.cells = 5, ngroups = 6) {
  cm <- GetAssayData(object = object[[assay]], slot = "counts")
  # subset to keep only genes detected in at least min.cells cells
  cm <- cm[rowSums(cm > 0) >= min.cells, ]
  # use the geometric mean of the features to group them
  # (using the arithmetic mean would usually not change things much)
  # could use sctransform:::row_gmean here but not exported
  feature.gmean <- exp(x = rowMeans(log1p(x = cm))) - 1
  feature.grp.breaks <- seq(
    from = min(log10(x = feature.gmean)) - 10*.Machine$double.eps,
    to = max(log10(x = feature.gmean)),
    length.out = ngroups + 1
  )
  feature.grp <- cut(
    x = log10(x = feature.gmean),
    breaks = feature.grp.breaks,
    ordered_result = TRUE
  )
  feature.grp <- factor(
    x = feature.grp,
    levels = rev(x = levels(x = feature.grp)),
    ordered = TRUE
  )
  names(x = feature.grp) <- names(x = feature.gmean)
  return(feature.grp)
}

# Get X and Y aesthetics from a plot for a certain geom
#
# @param plot A ggplot2 object
# @param geom Geom class to filter to
# @param plot.first Use plot-wide X/Y aesthetics before geom-specific aesthetics
#
# @return A named list with values 'x' for the name of the x aesthetic and 'y' for the y aesthetic
#
#' @importFrom rlang as_label
#
GetXYAesthetics <- function(plot, geom = 'GeomPoint', plot.first = TRUE) {
  geoms <- sapply(
    X = plot$layers,
    FUN = function(layer) {
      return(class(x = layer$geom)[1])
    }
  )
  # handle case where raster is set to True
  if (geom == "GeomPoint" && "GeomScattermore" %in% geoms){
    geom <- "GeomScattermore"
  }
  geoms <- which(x = geoms == geom)
  if (length(x = geoms) == 0) {
    stop("Cannot find a geom of class ", geom)
  }
  geoms <- min(geoms)
  if (plot.first) {
    # x <- as.character(x = plot$mapping$x %||% plot$layers[[geoms]]$mapping$x)[2]
    x <- as_label(x = plot$mapping$x %||% plot$layers[[geoms]]$mapping$x)
    # y <- as.character(x = plot$mapping$y %||% plot$layers[[geoms]]$mapping$y)[2]
    y <- as_label(x = plot$mapping$y %||% plot$layers[[geoms]]$mapping$y)
  } else {
    x <- as_label(x = plot$layers[[geoms]]$mapping$x %||% plot$mapping$x)
    y <- as_label(x = plot$layers[[geoms]]$mapping$y %||% plot$mapping$y)
  }
  return(list('x' = x, 'y' = y))
}

# For plotting the tissue image
#' @importFrom ggplot2 ggproto Geom aes ggproto_parent alpha draw_key_point
#' @importFrom grid unit gpar editGrob pointsGrob viewport gTree addGrob grobName
#'
GeomSpatial <- ggproto(
  "GeomSpatial",
  Geom,
  required_aes = c("x", "y"),
  extra_params = c("na.rm", "image", "image.alpha", "crop"),
  default_aes = aes(
    shape = 21,
    colour = "black",
    point.size.factor = 1.0,
    fill = NA,
    alpha = NA,
    stroke = 0.25
  ),
  setup_data = function(self, data, params) {
    data <- ggproto_parent(Geom, self)$setup_data(data, params)
    # We need to flip the image as the Y coordinates are reversed
    data$y = max(data$y) - data$y + min(data$y)
    data
  },
  draw_key = draw_key_point,
  draw_panel = function(data, panel_scales, coord, image, image.alpha, crop) {
    # This should be in native units, where
    # Locations and sizes are relative to the x- and yscales for the current viewport.
    if (!crop) {
      y.transform <- c(0, nrow(x = image)) - panel_scales$y.range
      data$y <- data$y + sum(y.transform)
      panel_scales$x$continuous_range <- c(0, ncol(x = image))
      panel_scales$y$continuous_range <- c(0, nrow(x = image))
      panel_scales$y.range <- c(0, nrow(x = image))
      panel_scales$x.range <- c(0, ncol(x = image))
    }
    z <- coord$transform(
      data.frame(x = c(0, ncol(x = image)), y = c(0, nrow(x = image))),
      panel_scales
    )
    # Flip Y axis for image
    z$y <- -rev(z$y) + 1
    wdth <- z$x[2] - z$x[1]
    hgth <- z$y[2] - z$y[1]
    vp <- viewport(
      x = unit(x = z$x[1], units = "npc"),
      y = unit(x = z$y[1], units = "npc"),
      width = unit(x = wdth, units = "npc"),
      height = unit(x = hgth, units = "npc"),
      just = c("left", "bottom")
    )
    img.grob <- GetImage(object = image)

    img <- editGrob(grob = img.grob, vp = vp)
    # spot.size <- slot(object = image, name = "spot.radius")
    spot.size <- Radius(object = image)
    coords <- coord$transform(data, panel_scales)
    pts <- pointsGrob(
      x = coords$x,
      y = coords$y,
      pch = data$shape,
      size = unit(spot.size, "npc") * data$point.size.factor,
      gp = gpar(
        col = alpha(colour = coords$colour, alpha = coords$alpha),
        fill = alpha(colour = coords$fill, alpha = coords$alpha),
        lwd = coords$stroke)
    )
    vp <- viewport()
    gt <- gTree(vp = vp)
    if (image.alpha > 0) {
      if (image.alpha != 1) {
        img$raster = as.raster(
          x = matrix(
            data = alpha(colour = img$raster, alpha = image.alpha),
            nrow = nrow(x = img$raster),
            ncol = ncol(x = img$raster),
            byrow = TRUE)
        )
      }
      gt <- addGrob(gTree = gt, child = img)
    }
    gt <- addGrob(gTree = gt, child = pts)
    # Replacement for ggname
    gt$name <- grobName(grob = gt, prefix = 'geom_spatial')
    return(gt)
    # ggplot2:::ggname("geom_spatial", gt)
  }
)

# influenced by: https://stackoverflow.com/questions/49475201/adding-tables-to-ggplot2-with-facet-wrap-in-r
# https://ggplot2.tidyverse.org/articles/extending-ggplot2.html
#' @importFrom ggplot2 layer
#'
#'
geom_spatial <-  function(
  mapping = NULL,
  data = NULL,
  image = image,
  image.alpha = image.alpha,
  crop = crop,
  stat = "identity",
  position = "identity",
  na.rm = FALSE,
  show.legend = NA,
  inherit.aes = TRUE,
  ...
) {
  layer(
    geom = GeomSpatial,
    mapping = mapping,
    data = data,
    stat = stat,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(na.rm = na.rm, image = image, image.alpha = image.alpha, crop = crop, ...)
  )
}

#' @importFrom grid viewport editGrob grobName
#' @importFrom ggplot2 ggproto Geom ggproto_parent
#
GeomSpatialInteractive <- ggproto(
  "GeomSpatialInteractive",
  Geom,
  setup_data = function(self, data, params) {
    data <- ggproto_parent(parent = Geom, self = self)$setup_data(data, params)
    data
  },
  draw_group = function(data, panel_scales, coord) {
    vp <- viewport(x = data$x, y = data$y)
    g <- editGrob(grob = data$grob[[1]], vp = vp)
    # Replacement for ggname
    g$name <- grobName(grob = g, prefix = 'geom_spatial_interactive')
    return(g)
    # return(ggname(prefix = "geom_spatial", grob = g))
  },
  required_aes = c("grob","x","y")
)

#' @importFrom ggplot2 layer
#
geom_spatial_interactive <-  function(
  mapping = NULL,
  data = NULL,
  stat = "identity",
  position = "identity",
  na.rm = FALSE,
  show.legend = NA,
  inherit.aes = FALSE,
  ...
) {
  layer(
    geom = GeomSpatialInteractive,
    mapping = mapping,
    data = data,
    stat = stat,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(na.rm = na.rm, ...)
  )
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

# Convert a ggplot2 scatterplot to base R graphics
#
# @param plot A ggplot2 scatterplot
# @param do.plot Create the plot with base R graphics
# @param cols A named vector of column names to pull. Vector names must be 'x',
# 'y', 'colour', 'shape', and/or 'size'; vector values must be the names of
# columns in plot data that correspond to these values. May pass only values that
# differ from the default (eg. \code{cols = c('size' = 'point.size.factor')})
# @param ... Extra parameters passed to PlotBuild
#
# @return A dataframe with the data that created the ggplot2 scatterplot
#
#' @importFrom ggplot2 ggplot_build
#
GGpointToBase <- function(
  plot,
  do.plot = TRUE,
  cols = c(
    'x' = 'x',
    'y' = 'y',
    'colour' = 'colour',
    'shape' = 'shape',
    'size' = 'size'
  ),
  ...
) {
  plot.build <- ggplot_build(plot = plot)
  default.cols <- c(
    'x' = 'x',
    'y' = 'y',
    'colour' = 'colour',
    'shape' = 'shape',
    'size' = 'size'
  )
  cols <- cols %||% default.cols
  if (is.null(x = names(x = cols))) {
    if (length(x = cols) > length(x = default.cols)) {
      warning(
        "Too many columns provided, selecting only first ",
        length(x = default.cols),
        call. = FALSE,
        immediate. = TRUE
      )
      cols <- cols[1:length(x = default.cols)]
    }
    names(x = cols) <- names(x = default.cols)[1:length(x = cols)]
  }
  cols <- c(
    cols[intersect(x = names(x = default.cols), y = names(x = cols))],
    default.cols[setdiff(x = names(x = default.cols), y = names(x = cols))]
  )
  cols <- cols[names(x = default.cols)]
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

# Convert a ggplot2 scatterplot to plotly graphics
#
# @inheritParams GGpointToBase
# @param information Extra information for hovering
# @param ... Ignored
#
# @return A dataframe with the data that greated the ggplot2 scatterplot
#' @importFrom ggplot2 ggplot_build
#
GGpointToPlotlyBuild <- function(
  plot,
  information = NULL,
  cols = eval(expr = formals(fun = GGpointToBase)$cols),
  ...
) {
  CheckDots(...)
  plot.build <- GGpointToBase(plot = plot, do.plot = FALSE, cols = cols)
  data <- ggplot_build(plot = plot)$plot$data
  rownames(x = plot.build) <- rownames(data)
  # Reset the names to 'x' and 'y'
  names(x = plot.build) <- c(
    'x',
    'y',
    names(x = plot.build)[3:length(x = plot.build)]
  )
  # Add the hover information we're looking for
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
    rownames(x = plot.build) <- plot.build$Row.names
    plot.build <- plot.build[, which(x = colnames(x = plot.build) != 'Row.names'), drop = FALSE]
  }
  return(plot.build)
}

#' @importFrom stats quantile
#'
InvertCoordinate <- function(x, MARGIN = 2) {
  if (!is.null(x = x)) {
    switch(
      EXPR = MARGIN,
      '1' = {
        rmin <- 'left'
        rmax <- 'right'
        cmin <- 'xmin'
        cmax <- 'xmax'
      },
      '2' = {
        rmin <- 'bottom'
        rmax <- 'top'
        cmin <- 'ymin'
        cmax <- 'ymax'
      },
      stop("'MARGIN' must be either 1 or 2", call. = FALSE)
    )
    # Fix the range so that rmin becomes rmax and vice versa
    # Needed for both points and brushes
    range <- x$range
    x$range[[rmin]] <- range[[rmax]]
    x$range[[rmax]] <- range[[rmin]]
    # Fix the cmin and cmax values, if provided
    # These are used for brush boundaries
    coords <- c(x[[cmin]], x[[cmax]])
    if (all(!is.null(x = coords))) {
      names(x = coords) <- c(cmin, cmax)
      x[[cmin]] <- quantile(
        x = x$range[[rmin]]:x$range[[rmax]],
        probs = 1 - (coords[cmax] / x$range[[rmax]]),
        names = FALSE
      )
      x[[cmax]] <- quantile(
        x = x$range[[rmin]]:x$range[[rmax]],
        probs = 1 - (coords[cmin] / x$range[[rmax]]),
        names = FALSE
      )
    }
  }
  return(x)
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
#' @importFrom stats median na.omit
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

# Plot expression of multiple features by identity on a plot
#
# @param data Data to plot
# @param idents Idents to use
# @param type Make either a 'ridge' or 'violin' plot
# @param sort Sort identity classes and features based on hierarchical clustering
# @param same.y.lims Indicates whether to use the same ylim for each feature
# @param adjust Adjust parameter for geom_violin
# @param cols Colors to use for plotting
# @param log plot Y axis on log scale
# @param fill.by Color violins/ridges based on either 'feature' or 'ident'
# @param seed.use Random seed to use. If NULL, don't set a seed
# @param flip flip plot orientation (identities on x-axis)
#
# @return A ggplot-based Expression-by-Identity plot
#
#' @importFrom cowplot theme_cowplot
#' @importFrom utils globalVariables
#' @importFrom stats rnorm dist hclust
#' @importFrom ggridges geom_density_ridges theme_ridges
#' @importFrom ggplot2 ggplot aes_string facet_grid theme labs geom_rect
#' geom_violin geom_jitter ylim position_jitterdodge scale_fill_manual
#' scale_y_log10 scale_x_log10 scale_y_discrete scale_x_continuous
#' scale_y_continuous waiver
#'
MultiExIPlot <- function(
  data,
  idents,
  split = NULL,
  type = 'violin',
  sort = FALSE,
  same.y.lims = same.y.lims,
  adjust = 1,
  pt.size = 0,
  cols = NULL,
  seed.use = 42,
  log = FALSE,
  fill.by = NULL,
  add.noise = TRUE,
  flip = NULL
) {
  if (!(fill.by %in% c("feature", "ident"))) {
    stop("`fill.by` must be either `feature` or `ident`")
  }
  if (!is.null(x = seed.use)) {
    set.seed(seed = seed.use)
  }
  if (!is.data.frame(x = data) || ncol(x = data) < 2) {
    stop("MultiExIPlot requires a data frame with >1 column")
  }
  data <- Melt(x = data)
  data <- data.frame(
    feature = data$cols,
    expression = data$vals,
    ident = rep_len(x = idents, length.out = nrow(x = data))
  )
  if ((is.character(x = sort) && nchar(x = sort) > 0) || sort) {
    data$feature <- as.vector(x = data$feature)
    data$ident <- as.vector(x = data$ident)
    # build matrix of average expression (#-features by #-idents), lexical ordering
    avgs.matrix <- sapply(
      X = split(x = data, f = data$ident),
      FUN = function(df) {
        return(tapply(
          X = df$expression,
          INDEX = df$feature,
          FUN = mean
        ))
      }
    )
    idents.order <- hclust(d = dist(x = t(x = L2Norm(mat = avgs.matrix, MARGIN = 2))))$order
    avgs.matrix <- avgs.matrix[,idents.order]
    avgs.matrix <- L2Norm(mat = avgs.matrix, MARGIN = 1)
    # order feature clusters by position of their "rank-1 idents"
    position <- apply(X = avgs.matrix, MARGIN = 1, FUN = which.max)
    mat <- hclust(d = dist(x = avgs.matrix))$merge
    orderings <- list()
    for (i in 1:nrow(mat)) {
      x <- if (mat[i,1] < 0) -mat[i,1] else orderings[[mat[i,1]]]
      y <- if (mat[i,2] < 0) -mat[i,2] else orderings[[mat[i,2]]]
      x.pos <- min(x = position[x])
      y.pos <- min(x = position[y])
      orderings[[i]] <- if (x.pos < y.pos) {
        c(x, y)
      } else {
        c(y, x)
      }
    }
    features.order <- orderings[[length(x = orderings)]]
    data$feature <- factor(
      x = data$feature,
      levels = unique(x = sort(x = data$feature))[features.order]
    )
    data$ident <- factor(
      x = data$ident,
      levels = unique(x = sort(x = data$ident))[rev(x = idents.order)]
    )
  } else {
    data$feature <- factor(x = data$feature, levels = unique(x = data$feature))
  }
  if (log) {
    noise <- rnorm(n = nrow(x = data)) / 200
    data$expression <- data$expression + 1
  } else {
    noise <- rnorm(n = nrow(x = data)) / 100000
  }
  if (!add.noise) {
    noise <- noise*0
  }
  for (f in unique(x = data$feature)) {
    if (all(data$expression[(data$feature == f)] == data$expression[(data$feature == f)][1])) {
      warning(
        "All cells have the same value of ",
        f,
        call. = FALSE,
        immediate. = TRUE
      )
    } else {
      data$expression[(data$feature == f)] <- data$expression[(data$feature == f)] + noise[(data$feature == f)]
    }
  }
  if (type == 'violin' && !is.null(x = split)) {
    data$split <- rep_len(x = split, length.out = nrow(data))
    vln.geom <- geom_violin
    fill.by <- 'split'
  } else if (type == 'splitViolin' && !is.null(x = split)) {
    data$split <- rep_len(x = split, length.out = nrow(data))
    vln.geom <- geom_split_violin
    fill.by <- 'split'
    type <- 'violin'
  } else {
    vln.geom <- geom_violin
  }
  switch(
    EXPR = type,
    'violin' = {
      geom <- list(vln.geom(scale = 'width', adjust = adjust, trim = TRUE))
    },
    'ridge' = {
      geom <- list(
        geom_density_ridges(scale = 4),
        theme_ridges(),
        scale_y_discrete(expand = c(0.01, 0))
      )
    },
    stop("Unknown plot type: ", type)
  )
  if (flip) {
    x <- 'ident'
    x.label <- 'Identity'
    y <- 'expression'
    y.label <- 'Expression Level'
  } else {
    y <- 'ident'
    y.label <- 'Identity'
    x <- 'expression'
    x.label <- 'Expression Level'
  }
  plot <- ggplot(
    data = data,
    mapping = aes_string(x = x, y = y, fill = fill.by)[c(2, 3, 1)]
  ) +
    labs(x = x.label, y = y.label, fill = NULL) +
    theme_cowplot()
  plot <- do.call(what = '+', args = list(plot, geom))
  if (flip) {
    plot <- plot +
      scale_y_continuous(
        expand = c(0, 0),
        labels = function(x) c(rep(x = '', times = length(x)-2), x[length(x) - 1], '')) +
      facet_grid(feature ~ ., scales = (if (same.y.lims) 'fixed' else 'free')) +
      FacetTheme(
        panel.spacing = unit(0, 'lines'),
        panel.background = element_rect(fill = NA, color = "black"),
        axis.text.y = element_text(size = 7),
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text.y.right = element_text(angle = 0))
  } else {
    plot <- plot +
      scale_x_continuous(
        expand = c(0, 0),
        labels = function(x) c(rep(x = '', times = length(x)-2), x[length(x) - 1], '')) +
      facet_grid(. ~ feature, scales = (if (same.y.lims) 'fixed' else 'free')) +
      FacetTheme(
        panel.spacing = unit(0, 'lines'),
        panel.background = element_rect(fill = NA, color = "black"),
        axis.text.x = element_text(size = 7),
        strip.text.x = element_text(angle = -90))
  }
  if (log) {
    plot <- plot + scale_x_log10()
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
  .Defunct(new = "CellSelector")
  # #   Convert the ggplot object to a data.frame
  # PackageCheck('SDMTools')
  # plot.data <- GGpointToBase(plot = plot, dark.theme = dark.theme, ...)
  # npoints <- nrow(x = plot.data)
  # cat("Click around the cluster of points you wish to select\n")
  # cat("ie. select the vertecies of a shape around the cluster you\n")
  # cat("are interested in. Press <Esc> when finished (right click for R-terminal users)\n\n")
  # polygon <- locator(n = npoints, type = 'l')
  # polygon <- data.frame(polygon)
  # #   pnt.in.poly returns a data.frame of points
  # points.all <- SDMTools::pnt.in.poly(
  #   pnts = plot.data[, c(1, 2)],
  #   poly.pnts = polygon
  # )
  # #   Find the located points
  # points.located <- points.all[which(x = points.all$pip == 1), ]
  # #   If we're recoloring, do the recolor
  # if (recolor) {
  #   no <- ifelse(test = dark.theme, yes = 'white', no = '#C3C3C3')
  #   points.all$color <- ifelse(test = points.all$pip == 1, yes = '#DE2D26', no = no)
  #   plot.data$color <- points.all$color
  #   PlotBuild(data = plot.data, dark.theme = dark.theme, ...)
  # }
  # return(points.located[, c(1, 2)])
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
  highlight <- factor(x = highlight, levels = plot.order)
  return(list(
    plot.order = plot.order,
    highlight = highlight,
    size = size,
    color = cols.highlight
  ))
}

#' @importFrom shiny brushedPoints
#
ShinyBrush <- function(plot.data, brush, outputs, inverts = character(length = 0L)) {#}, selected = NULL) {
  selected <- NULL
  if (!is.null(x = brush)) {
    if (brush$outputId %in% outputs) {
      selected <- rownames(x = brushedPoints(df = plot.data, brush = brush))
    } else if (brush$outputId %in% inverts) {
      selected <- rownames(x = brushedPoints(
        df = plot.data,
        brush = InvertCoordinate(x = brush)
      ))
    }
  }
  return(selected)
}

globalVariables(names = '..density..', package = 'Seurat')
#' A single correlation plot
#'
#' @param data A data frame with two columns to be plotted
#' @param col.by A vector or factor of values to color the plot by
#' @param cols An optional vector of colors to use
#' @param pt.size Point size for the plot
#' @param smooth Make a smoothed scatter plot
#' @param rows.highight A vector of rows to highlight (like cells.highlight in
#' \code{\link{SingleDimPlot}})
#' @param legend.title Optional legend title
#' @param raster Convert points to raster format, default is \code{NULL}
#' which will automatically use raster if the number of points plotted is
#' greater than 100,000
#' @param raster.dpi the pixel resolution for rastered plots, passed to geom_scattermore().
#' Default is c(512, 512)
#' @param plot.cor ...
#' @param jitter Jitter for easier visualization of crowded points
#'
#' @return A ggplot2 object
#'
#' @importFrom stats cor
#' @importFrom cowplot theme_cowplot
#' @importFrom RColorBrewer brewer.pal.info
#' @importFrom ggplot2 ggplot aes_string geom_point labs scale_color_brewer
#' scale_color_manual guides stat_density2d aes scale_fill_continuous
#' @importFrom scattermore geom_scattermore
#'
#' @keywords internal
#'
#' @export
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
  span = NULL,
  raster = NULL,
  raster.dpi = NULL,
  plot.cor = TRUE,
  jitter = TRUE
) {
  pt.size <- pt.size %||% AutoPointSize(data = data, raster = raster)
  if ((nrow(x = data) > 1e5) & !isFALSE(raster)){
    message("Rasterizing points since number of points exceeds 100,000.",
            "\nTo disable this behavior set `raster=FALSE`")
  }
  raster <- raster %||% (nrow(x = data) > 1e5)
  if (!is.null(x = raster.dpi)) {
    if (!is.numeric(x = raster.dpi) || length(x = raster.dpi) != 2)
      stop("'raster.dpi' must be a two-length numeric vector")
  }
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
  plot.cor <- if (isTRUE(x = plot.cor)) {
    round(x = cor(x = data[, 1], y = data[, 2]), digits = 2)
  }
  else(
    ""
  )
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
  position <- NULL
  if (jitter) {
    position <- 'jitter'
  } else {
    position <- 'identity'
  }
  if (!is.null(x = col.by)) {
    if (raster) {
      plot <- plot + geom_scattermore(
        mapping = aes_string(color = 'colors'),
        position = position,
        pointsize = pt.size,
        pixels = raster.dpi
      )
    } else {
      plot <- plot + geom_point(
        mapping = aes_string(color = 'colors'),
        position = position,
        size = pt.size
      )
    }
  } else {
    if (raster) {
      plot <- plot + geom_scattermore(position = position, pointsize = pt.size, pixels = raster.dpi)
    } else {
      plot <- plot + geom_point(position = position, size = pt.size)
    }
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

#' Plot a single dimension
#'
#' @param data Data to plot
#' @param dims A two-length numeric vector with dimensions to use
#' @param col.by ...
#' @param cols Vector of colors, each color corresponds to an identity class.
#' This may also be a single character or numeric value corresponding to a
#' palette as specified by \code{\link[RColorBrewer]{brewer.pal.info}}.By
#' default, ggplot2 assigns colors
#' @param pt.size Adjust point size for plotting
#' @param shape.by If NULL, all points are circles (default). You can specify
#' any cell attribute (that can be pulled with \code{\link{FetchData}})
#' allowing for both different colors and different shapes on cells.
#' @param alpha Alpha value for plotting (default is 1)
#' @param alpha.by Mapping variable for the point alpha value
#' @param order Specify the order of plotting for the idents. This can be
#' useful for crowded plots if points of interest are being buried. Provide
#' either a full list of valid idents or a subset to be plotted last (on top).
#' @param label Whether to label the clusters
#' @param repel Repel labels
#' @param label.size Sets size of labels
#' @param cells.highlight A list of character or numeric vectors of cells to
#' highlight. If only one group of cells desired, can simply
#' pass a vector instead of a list. If set, colors selected cells to the color(s)
#' in \code{cols.highlight} and other cells black (white if dark.theme = TRUE);
#' will also resize to the size(s) passed to \code{sizes.highlight}
#' @param cols.highlight A vector of colors to highlight the cells as; will
#' repeat to the length groups in cells.highlight
#' @param sizes.highlight Size of highlighted cells; will repeat to the length
#' groups in cells.highlight
#' @param na.value Color value for NA points when using custom scale.
#' @param raster Convert points to raster format, default is \code{NULL}
#' which will automatically use raster if the number of points plotted is
#' greater than 100,000
#' @param raster.dpi the pixel resolution for rastered plots, passed to geom_scattermore().
#' Default is c(512, 512)
#'
#' @return A ggplot2 object
#'
#' @importFrom cowplot theme_cowplot
#' @importFrom RColorBrewer brewer.pal.info
#' @importFrom ggplot2 ggplot aes_string geom_point labs guides scale_color_brewer
#' scale_color_manual element_rect guide_legend discrete_scale
#'
#' @keywords internal
#'
#' @export
#'
SingleDimPlot <- function(
  data,
  dims,
  col.by = NULL,
  cols = NULL,
  pt.size = NULL,
  shape.by = NULL,
  alpha = 1,
  alpha.by = NULL,
  order = NULL,
  label = FALSE,
  repel = FALSE,
  label.size = 4,
  cells.highlight = NULL,
  cols.highlight = '#DE2D26',
  sizes.highlight = 1,
  na.value = 'grey50',
  raster = NULL,
  raster.dpi = NULL
) {
  pt.size <- pt.size %||% AutoPointSize(data = data, raster = raster)
  if ((nrow(x = data) > 1e5) & !isFALSE(raster)){
    message("Rasterizing points since number of points exceeds 100,000.",
            "\nTo disable this behavior set `raster=FALSE`")
  }
  raster <- raster %||% (nrow(x = data) > 1e5)
  if (!is.null(x = raster.dpi)) {
    if (!is.numeric(x = raster.dpi) || length(x = raster.dpi) != 2)
      stop("'raster.dpi' must be a two-length numeric vector")
  }
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
    if (inherits(x = cells.highlight, what = "data.frame")) {
      stop("cells.highlight cannot be a dataframe. ",
           "Please supply a vector or list")
    }
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
        data <- data[order(!is.na(x = data[, col.by]), data[, col.by]), ]
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
  if (!is.null(x = alpha.by) && !alpha.by %in% colnames(x = data)) {
    warning(
      "Cannot find alpha variable ",
      alpha.by,
      " in data, setting to NULL",
      call. = FALSE,
      immediate. = TRUE
    )
    alpha.by <- NULL
  }

  plot <- ggplot(data = data)
  plot <- if (isTRUE(x = raster)) {
    plot + geom_scattermore(
      mapping = aes_string(
        x = dims[1],
        y = dims[2],
        color = paste0("`", col.by, "`"),
        shape = shape.by,
        alpha = alpha.by
      ),
      pointsize = pt.size,
      alpha = alpha,
      pixels = raster.dpi
    )
  } else {
    plot + geom_point(
      mapping = aes_string(
        x = dims[1],
        y = dims[2],
        color = paste0("`", col.by, "`"),
        shape = shape.by,
        alpha = alpha.by
      ),
      size = pt.size,
      alpha = alpha
    )
  }
  plot <- plot +
    guides(color = guide_legend(override.aes = list(size = 3, alpha = 1))) +
    labs(color = NULL, title = col.by) +
    CenterTitle()
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

#' Plot a single expression by identity on a plot
#'
#' @param data Data to plot
#' @param idents Idents to use
#' @param split Use a split violin plot
#' @param type Make either a \dQuote{ridge} or \dQuote{violin} plot
#' @param sort Sort identity classes (on the x-axis) by the average
#' expression of the attribute being potted
#' @param y.max Maximum Y value to plot
#' @param adjust Adjust parameter for geom_violin
#' @param pt.size Size of points for violin plots
#' @param alpha Alpha vlaue for violin plots
#' @param cols Colors to use for plotting
#' @param seed.use Random seed to use. If NULL, don't set a seed
#' @param log plot Y axis on log10 scale
#' @param add.noise determine if adding small noise for plotting
#' @param raster Convert points to raster format. Requires 'ggrastr' to be installed.
#' default is \code{NULL} which automatically rasterizes if ggrastr is installed and
#' number of points exceed 100,000.
#'
#' @return A ggplot-based Expression-by-Identity plot
#'
#' @importFrom stats rnorm
#' @importFrom utils globalVariables
#' @importFrom ggridges geom_density_ridges theme_ridges
#' @importFrom ggplot2 ggplot aes_string theme labs geom_violin geom_jitter
#' ylim position_jitterdodge scale_fill_manual scale_y_log10 scale_x_log10
#' scale_y_discrete scale_x_continuous waiver
#' @importFrom cowplot theme_cowplot
#'
#' @keywords internal
#' @export
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
  alpha = 1,
  cols = NULL,
  seed.use = 42,
  log = FALSE,
  add.noise = TRUE,
  raster = NULL
) {
   if (!is.null(x = raster) && isTRUE(x = raster)){
    if (!PackageCheck('ggrastr', error = FALSE)) {
      stop("Please install ggrastr from CRAN to enable rasterization.")
    }
  }
  if (PackageCheck('ggrastr', error = FALSE)) {
    # Set rasterization to true if ggrastr is installed and
    # number of points exceeds 100,000
    if ((nrow(x = data) > 1e5) & !isFALSE(raster)){
      message("Rasterizing points since number of points exceeds 100,000.",
              "\nTo disable this behavior set `raster=FALSE`")
      # change raster to TRUE
      raster <- TRUE
    }
  }
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
  if (!add.noise) {
    noise <-  noise * 0
  }
  if (all(data[, feature] == data[, feature][1])) {
    warning(paste0("All cells have the same value of ", feature, "."))
  } else{
    data[, feature] <- data[, feature] + noise
  }
  axis.label <- 'Expression Level'
  y.max <- y.max %||% max(data[, feature][is.finite(x = data[, feature])])
  if (type == 'violin' && !is.null(x = split)) {
    data$split <- split
    vln.geom <- geom_violin
    fill <- 'split'
  } else if (type == 'splitViolin' && !is.null(x = split )) {
    data$split <- split
    vln.geom <- geom_split_violin
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
      if (is.null(x = split)) {
        if (isTRUE(x = raster)) {
          jitter <- ggrastr::rasterize(geom_jitter(height = 0, size = pt.size, alpha = alpha, show.legend = FALSE))
        } else {
          jitter <- geom_jitter(height = 0, size = pt.size, alpha = alpha, show.legend = FALSE)
        }
      } else {
        if (isTRUE(x = raster)) {
          jitter <- ggrastr::rasterize(geom_jitter(
            position = position_jitterdodge(jitter.width = 0.4, dodge.width = 0.9),
            size = pt.size,
            alpha = alpha,
            show.legend = FALSE
          ))
        } else {
          jitter <- geom_jitter(
            position = position_jitterdodge(jitter.width = 0.4, dodge.width = 0.9),
            size = pt.size,
            alpha = alpha,
            show.legend = FALSE
          )
        }
      }
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
      jitter <- geom_jitter(width = 0, size = pt.size, alpha = alpha, show.legend = FALSE)
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

#' A single heatmap from base R using \code{\link[graphics]{image}}
#'
#' @param data matrix of data to plot
#' @param order optional vector of cell names to specify order in plot
#' @param title Title for plot
#'
#' @return No return, generates a base-R heatmap using \code{\link[graphics]{image}}
#'
#' @importFrom graphics axis image par plot.new title
#'
#' @keywords internal
#'
#' @export
#'
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

#' Single Spatial Plot
#'
#' @param data A data frame with at least the following columns:
#' \itemize{
#'  \item \dQuote{\code{x}}: Spatial-resolved \emph{x} coordinates, will be
#'   plotted on the \emph{y}-axis
#'  \item \dQuote{\code{y}}: Spatially-resolved \emph{y} coordinates, will be
#'   plotted on the \emph{x}-axis
#'  \item \dQuote{\code{cell}}: Cell name
#'  \item \dQuote{\code{boundary}}: Segmentation boundary label; when plotting
#'   multiple segmentation layers, the order of boundary transparency is set by
#'   factor levels for this column
#' }
#' Can pass \code{NA} to \code{data} suppress segmentation visualization
#' @param col.by Name of column in \code{data} to color cell segmentations by;
#' pass \code{NA} to suppress coloring
#' @param col.factor Are the colors a factor or discrete?
#' @param cols Colors for cell segmentations; can be one of the
#' following:
#' \itemize{
#'  \item \code{NULL} for default ggplot2 colors
#'  \item A numeric value or name of a
#'   \link[RColorBrewer:RColorBrewer]{color brewer palette}
#'  \item Name of a palette for \code{\link{DiscretePalette}}
#'  \item A vector of colors equal to the length of unique levels of
#'   \code{data$col.by}
#' }
#' @param shuffle.cols Randomly shuffle colors when a palette or
#' vector of colors is provided to \code{cols}
#' @param size Point size for cells when plotting centroids
#' @param molecules A data frame with spatially-resolved molecule coordinates;
#' should have the following columns:
#' \itemize{
#'  \item \dQuote{\code{x}}: Spatial-resolved \emph{x} coordinates, will be
#'   plotted on the \emph{y}-axis
#'  \item \dQuote{\code{y}}: Spatially-resolved \emph{y} coordinates, will be
#'   plotted on the \emph{x}-axis
#'  \item \dQuote{\code{molecule}}: Molecule name
#' }
#' @param mols.size Point size for molecules
#' @param mols.cols A vector of color for molecules. The "Set1" palette from
#' RColorBrewer is used by default.
#' @param mols.alpha Alpha value for molecules, should be between 0 and 1
#' @param alpha Alpha value, should be between 0 and 1; when plotting multiple
#' boundaries, \code{alpha} is equivalent to max alpha
#' @param border.color Color of cell segmentation border; pass \code{NA}
#' to suppress borders for segmentation-based plots
#' @param border.size Thickness of cell segmentation borders; pass \code{NA}
#' to suppress borders for centroid-based plots
#' @param na.value Color value for \code{NA} segmentations when
#' using custom scale
#' @param ... Ignored
#'
#' @return A ggplot object
#'
#' @importFrom rlang is_na
#' @importFrom SeuratObject %NA% %!NA%
#' @importFrom RColorBrewer brewer.pal.info
#' @importFrom ggplot2 aes_string geom_point geom_polygon ggplot guides
#' guide_legend scale_alpha_manual scale_color_manual scale_fill_brewer
#' scale_fill_manual
#'
#' @keywords internal
#'
#' @export
#'
SingleImagePlot <- function(
  data,
  col.by = NA,
  col.factor = TRUE,
  cols = NULL,
  shuffle.cols = FALSE,
  size = 0.1,
  molecules = NULL,
  mols.size = 0.1,
  mols.cols = NULL,
  mols.alpha = 1.0,
  alpha = molecules %iff% 0.3 %||% 0.6,
  border.color = 'white',
  border.size = NULL,
  na.value = 'grey50',
  dark.background = TRUE,
  ...
) {
  # Check input data
  if (!is_na(x = data)) {
    if (!all(c('x', 'y', 'cell', 'boundary') %in% colnames(x = data))) {
      stop("Invalid data coordinates")
    }
    if (!is_na(x = col.by)) {
      if (!col.by  %in% colnames(x = data)) {
        warning(
          "Cannot find 'col.by' ('",
          col.by,
          "') in data coordinates",
          immediate. = TRUE
        )
        col.by <- NA
      } else if (isTRUE(x = col.factor) && !is.factor(x = data[[col.by]])) {
        data[[col.by]] <- factor(
          x = data[[col.by]],
          levels = unique(x = data[[col.by]])
        )
      } else if (isFALSE(x = col.factor) && is.factor(x = data[[col.by]])) {
        data[[col.by]] <- as.vector(x = data[[col.by]])
      }
    }
    if (is_na(x = col.by) && !is.null(x = cols)) {
      col.by <- RandomName(length = 7L)
      data[[col.by]] <- TRUE
    }
    if (!is.factor(x = data$boundary)) {
      data$boundary <- factor(
        x = data$boundary,
        levels = unique(x = data$boundary)
      )
    }
    # Determine alphas
    if (is.na(x = alpha)) {
      alpha <- 1L
    }
    alpha.min <- ifelse(
      test = alpha < 1L,
      yes = 1 * (10 ^ .FindE(x = alpha)),
      no = 0.1
    )
    if (alpha.min == alpha) {
      alpha.min <- 1 * (10 ^ (.FindE(x = alpha) - 1))
    }
    alphas <- .Cut(
      min = alpha.min,
      max = alpha,
      n = length(x = levels(x = data$boundary))
    )
  }
  # Assemble plot
  plot <- ggplot(
    data = data %NA% NULL,
    mapping = aes_string(
      x = 'y',
      y = 'x',
      alpha = 'boundary',
      fill = col.by %NA% NULL
    )
  )
  if (!is_na(x = data)) {
    plot <- plot +
      scale_alpha_manual(values = alphas) +
      if (anyDuplicated(x = data$cell)) {
        if (is.null(border.size)) {
          border.size <- 0.3
        }
        geom_polygon(
          mapping = aes_string(group = 'cell'),
          color = border.color,
          size = border.size
        )
      } else {
        # Default to no borders when plotting centroids
        if (is.null(border.size)) {
          border.size <- 0.0
        }
        geom_point(
          shape = 21,
          color = border.color,
          stroke = border.size,
          size = size
        )
      }
    if (!is.null(x = cols)) {
      plot <- plot + if (is.numeric(x = cols) || cols[1L] %in% rownames(x = brewer.pal.info)) {
        palette <- brewer.pal(n = length(x = levels(x = data[[col.by]])), cols)
        if (length(palette) < length(x = levels(x = data[[col.by]]))) {
          num.blank <- length(x = levels(x = data[[col.by]])) - length(palette)
          palette <- c(palette, rep(na.value, num.blank))
        }
        if (isTRUE(shuffle.cols)) {
          palette <- sample(palette)
        }
        scale_fill_manual(values = palette, na.value = na.value)
      } else if (cols[1] %in% DiscretePalette(n = NULL)) {
        scale_fill_manual(
          values = DiscretePalette(
            n = length(x = levels(x = data[[col.by]])),
            palette = cols,
            shuffle = shuffle.cols
          ),
          na.value = na.value
        )
      } else {
        if (isTRUE(shuffle.cols)) {
          cols <- sample(cols)
        }
        scale_fill_manual(values = cols, na.value = na.value)
      }
    }
    if (length(x = levels(x = data$boundary)) == 1L) {
      plot <- plot + guides(alpha = 'none')
    }
    # Adjust guides
    if (isTRUE(x = col.factor) && length(x = levels(x = data[[col.by]])) <= 1L) {
      plot <- plot + guides(fill = 'none')
    }
  }
  # Add molecule sets
  if (is.data.frame(x = molecules)) {
    if (all(c('x', 'y', 'molecule') %in% colnames(x = molecules))) {
      if (!is.factor(x = molecules$molecule)) {
        molecules$molecule <- factor(
          x = molecules$molecule,
          levels = unique(x = molecules$molecule)
        )
      }
      plot <- plot + geom_point(
        mapping = aes_string(fill = NULL, alpha = NULL, color = "molecule"),
        data = molecules,
        size = mols.size,
        alpha = mols.alpha,
        show.legend = c(color = TRUE, fill = FALSE, alpha = FALSE)
      ) +
        guides(color = guide_legend(override.aes = list(size = 3L))) +
        scale_color_manual(
          name = 'Molecules',
          values = mols.cols %||% brewer.pal(
            n = length(x = levels(x = molecules$molecule)),
            name = "Set1"
          ),
          guide = guide_legend()
        )
    } else {
      warning("Invalid molecule coordinates", immediate. = TRUE)
    }
  }

  if (isTRUE(dark.background)) {
    plot <- plot + DarkTheme()
  }
  return(plot)
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

#' A single heatmap from ggplot2 using geom_raster
#'
#' @param data A matrix or data frame with data to plot
#' @param raster switch between geom_raster and geom_tile
#' @param cell.order ...
#' @param feature.order ...
#' @param colors A vector of colors to use
#' @param disp.min Minimum display value (all values below are clipped)
#' @param disp.max Maximum display value (all values above are clipped)
#' @param limits A two-length numeric vector with the limits for colors on the plot
#' @param group.by A vector to group cells by, should be one grouping identity per cell
#'
#' @return A ggplot2 object
#
#' @importFrom ggplot2 ggplot aes_string geom_raster scale_fill_gradient
#' scale_fill_gradientn theme element_blank labs geom_point guides
#' guide_legend geom_tile
#'
#' @keywords internal
#'
#' @export
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

#' Base plotting function for all Spatial plots
#'
#' @param data Data.frame with info to be plotted
#' @param image \code{SpatialImage} object to be plotted
#' @param cols Vector of colors, each color corresponds to an identity class.
#' This may also be a single character
#' or numeric value corresponding to a palette as specified by
#' \code{\link[RColorBrewer]{brewer.pal.info}}. By default, ggplot2 assigns
#' colors
#' @param image.alpha Adjust the opacity of the background images. Set to 0 to
#' remove.
#' @param pt.alpha Adjust the opacity of the points if plotting a
#' \code{SpatialDimPlot}
#' @param crop Crop the plot in to focus on points plotted. Set to \code{FALSE}
#' to show entire background image.
#' @param pt.size.factor Sets the size of the points relative to spot.radius
#' @param stroke Control the width of the border around the spots
#' @param col.by Mapping variable for the point color
#' @param alpha.by Mapping variable for the point alpha value
#' @param cells.highlight A list of character or numeric vectors of cells to
#' highlight. If only one group of cells desired, can simply pass a vector
#' instead of a list. If set, colors selected cells to the color(s) in
#' cols.highlight
#' @param cols.highlight A vector of colors to highlight the cells as; ordered
#' the same as the groups in cells.highlight; last color corresponds to
#' unselected cells.
#' @param geom Switch between normal spatial geom and geom to enable hover
#' functionality
#' @param na.value Color for spots with NA values
#'
#' @return A ggplot2 object
#'
#' @importFrom tibble tibble
#' @importFrom ggplot2 ggplot aes_string coord_fixed geom_point xlim ylim
#' coord_cartesian labs theme_void theme scale_fill_brewer
#'
#' @keywords internal
#'
#' @export
#'
SingleSpatialPlot <- function(
  data,
  image,
  cols = NULL,
  image.alpha = 1,
  pt.alpha = NULL,
  crop = TRUE,
  pt.size.factor = NULL,
  stroke = 0.25,
  col.by = NULL,
  alpha.by = NULL,
  cells.highlight = NULL,
  cols.highlight = c('#DE2D26', 'grey50'),
  geom = c('spatial', 'interactive', 'poly'),
  na.value = 'grey50'
) {
  geom <- match.arg(arg = geom)
  if (!is.null(x = col.by) && !col.by %in% colnames(x = data)) {
    warning("Cannot find '", col.by, "' in data, not coloring", call. = FALSE, immediate. = TRUE)
    col.by <- NULL
  }
  col.by <- col.by %iff% paste0("`", col.by, "`")
  alpha.by <- alpha.by %iff% paste0("`", alpha.by, "`")
  if (!is.null(x = cells.highlight)) {
    highlight.info <- SetHighlight(
      cells.highlight = cells.highlight,
      cells.all = rownames(x = data),
      sizes.highlight = pt.size.factor,
      cols.highlight = cols.highlight[1],
      col.base = cols.highlight[2]
    )
    order <- highlight.info$plot.order
    data$highlight <- highlight.info$highlight
    col.by <- 'highlight'
    levels(x = data$ident) <- c(order, setdiff(x = levels(x = data$ident), y = order))
    data <- data[order(data$ident), ]
  }
  plot <- ggplot(data = data, aes_string(
    x = colnames(x = data)[2],
    y = colnames(x = data)[1],
    fill = col.by,
    alpha = alpha.by
  ))
  plot <- switch(
    EXPR = geom,
    'spatial' = {
      if (is.null(x = pt.alpha)) {
        plot <- plot + geom_spatial(
          point.size.factor = pt.size.factor,
          data = data,
          image = image,
          image.alpha = image.alpha,
          crop = crop,
          stroke = stroke,
        )
      } else {
        plot <- plot + geom_spatial(
          point.size.factor = pt.size.factor,
          data = data,
          image = image,
          image.alpha = image.alpha,
          crop = crop,
          stroke = stroke,
          alpha = pt.alpha
        )
      }
      plot + coord_fixed() + theme(aspect.ratio = 1)
    },
    'interactive' = {
      plot + geom_spatial_interactive(
        data = tibble(grob = list(GetImage(object = image, mode = 'grob'))),
        mapping = aes_string(grob = 'grob'),
        x = 0.5,
        y = 0.5
      ) +
        geom_point(mapping = aes_string(color = col.by)) +
        xlim(0, ncol(x = image)) +
        ylim(nrow(x = image), 0) +
        coord_cartesian(expand = FALSE)
    },
    'poly' = {
      data$cell <- rownames(x = data)
      data[, c('x', 'y')] <- NULL
      data <- merge(
        x = data,
        y = GetTissueCoordinates(object = image, qhulls = TRUE),
        by = "cell"
      )
      plot + geom_polygon(
        data = data,
        mapping = aes_string(fill = col.by, group = 'cell')
      ) + coord_fixed() + theme_cowplot()

    },
    stop("Unknown geom, choose from 'spatial' or 'interactive'", call. = FALSE)
  )
  if (!is.null(x = cells.highlight)) {
    plot <- plot + scale_fill_manual(values = cols.highlight)
  }
  if (!is.null(x = cols) && is.null(x = cells.highlight)) {
    if (length(x = cols) == 1 && (is.numeric(x = cols) || cols %in% rownames(x = brewer.pal.info))) {
      scale <- scale_fill_brewer(palette = cols, na.value = na.value)
    } else if (length(x = cols) == 1 && (cols %in% c('alphabet', 'alphabet2', 'glasbey', 'polychrome', 'stepped'))) {
      colors <- DiscretePalette(length(unique(data[[col.by]])), palette = cols)
      scale <- scale_fill_manual(values = colors, na.value = na.value)
    } else {
      cols <- cols[names(x = cols) %in% data[[gsub(pattern = '`', replacement = "", x = col.by)]]]
      scale <- scale_fill_manual(values = cols, na.value = na.value)
    }
    plot <- plot + scale
  }
  plot <- plot + NoAxes() + theme(panel.background = element_blank())
  return(plot)
}

# Reimplementation of ggplot2 coord$transform
#
# @param data A data frame with x-coordinates in the first column and y-coordinates
# in the second
# @param xlim,ylim X- and Y-limits for the transformation, must be two-length
# numeric vectors
#
# @return \code{data} with transformed coordinates
#
#' @importFrom ggplot2 transform_position
#' @importFrom scales rescale squish_infinite
#
Transform <- function(data, xlim = c(-Inf, Inf), ylim = c(-Inf, Inf)) {
  # Quick input argument checking
  if (!all(sapply(X = list(xlim, ylim), FUN = length) == 2)) {
    stop("'xlim' and 'ylim' must be two-length numeric vectors", call. = FALSE)
  }
  # Save original names
  df.names <- colnames(x = data)
  colnames(x = data)[1:2] <- c('x', 'y')
  # Rescale the X and Y values
  data <- transform_position(
    df = data,
    trans_x = function(df) {
      return(rescale(x = df, from = xlim))
    },
    trans_y = function(df) {
      return(rescale(x = df, from = ylim))
    }
  )
  # Something that ggplot2 does
  data <- transform_position(
    df = data,
    trans_x = squish_infinite,
    trans_y = squish_infinite
  )
  # Restore original names
  colnames(x = data) <- df.names
  return(data)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# S4 Methods
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

setMethod(
  f = '.PrepImageData',
  signature = c(data = 'data.frame', cells = 'rle'),
  definition = function(data, cells, ...) {
    data <- sapply(
      X = colnames(x = data),
      FUN = function(x) {
        j <- data[[x]]
        names(x = j) <- rownames(x = data)
        return(.PrepImageData(data = j, cells = cells, name = x))
      },
      simplify = FALSE,
      USE.NAMES = TRUE
    )
    return(do.call(what = 'cbind', args = data))
  }
)

#' @importFrom methods getMethod
setMethod(
  f = '.PrepImageData',
  signature = c(data = 'factor', cells = 'rle'),
  definition = function(data, cells, name, ...) {
    f <- getMethod(
      f = '.PrepImageData',
      signature = c(data = 'vector', cells = 'rle')
    )
    return(f(data = data, cells = cells, name = name, ...))
  }
)

setMethod(
  f = '.PrepImageData',
  signature = c(data = 'list', cells = 'rle'),
  definition = function(data, cells, ...) {
    .NotYetImplemented()
  }
)

setMethod(
  f = '.PrepImageData',
  signature = c(data = 'NULL', cells = 'rle'),
  definition = function(data, cells, ...) {
    return(SeuratObject::EmptyDF(n = sum(cells$lengths)))
  }
)

setMethod(
  f = '.PrepImageData',
  signature = c(data = 'vector', cells = 'rle'),
  definition = function(data, cells, name, ...) {
    name <- as.character(x = name)
    if (name %in% c('x', 'y', 'cell')) {
      stop("'name' cannot be 'x', 'y', or 'cell'", call. = FALSE)
    }
    cnames <- cells$values
    if (is.null(x = names(x = data))) {
      mlen <- min(sapply(X = list(data, cnames), FUN = length))
      names(x = data)[1:mlen] <- cnames[1:mlen]
    }
    if (anyDuplicated(x = names(x = data))) {
      dup <- duplicated(x = names(x = data))
      warning(
        sum(dup),
        ifelse(test = sum(dup) == 1, yes = ' cell', no = ' cells'),
        ' duplicated, using only the first occurance',
        call. = FALSE,
        immediate. = TRUE
      )
    }
    if (length(x = data) < length(x = cnames)) {
      mcells <- setdiff(x = cnames, y = names(x = data))
      warning(
        "Missing data for some cells, filling with NA",
        call. = FALSE,
        immediate. = TRUE
      )
      data[mcells] <- NA
    } else if (length(x = data) > length(x = cnames)) {
      warning(
        "More cells provided than present",
        call. = FALSE,
        immediate. = TRUE
      )
    }
    data <- data.frame(rep.int(x = data[cnames], times = cells$lengths))
    colnames(x = data) <- name
    return(data)
  }
)
