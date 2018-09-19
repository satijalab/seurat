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
#' @param ncol Number of columns to plot
#' @param fast If true, use \code{image} to generate plots; faster than using ggplot2, but not customizable
#' @param assays A vector of assays to pull data from
#'
#' @return A ggplot object
#'
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
  disp.max = 2.5,
  balanced = FALSE,
  ncol = NULL,
  combine = TRUE,
  fast = TRUE,
  slot = 'scale.data',
  assays = NULL,
  ...
) {
  ncol <- ncol %||% ifelse(test = length(x = dims) > 2, yes = 3, no = length(x = dims))
  plots <- vector(mode = 'list', length = length(x = dims))
  assays <- assays %||% DefaultAssay(object = object)
  if (!DefaultAssay(object = object[[reduction]]) %in% assays) {
    warning("assay")
  }
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
  cells <- cells %||% colnames(x = object)
  if (!is.list(x = cells)) {
    cells <- lapply(X = 1:length(x = dims), FUN = function(x) {return(cells)})
  }
  features <- lapply(
    X = dims,
    FUN = TopFeatures,
    object = object[[reduction]],
    nfeatures = nfeatures,
    balanced = balanced
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
    dim.features <- unname(obj = unlist(x = rev(x = features[[i]])))
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
    plots <- CombinePlots(
      plots = plots,
      ncol = ncol,
      legend = 'right'
    )
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
#' @param disp.max Maximum display value (all values above are clipped)
#' @param group.by A vector of variables to group cells by; pass 'ident' to group by cell identity classes
#' @param group.bar Add a color bar showing group status for cells
#' @param slot Data slot to use, choose from 'raw.data', 'data', or 'scale.data'
#' @param assay Assay to pull from
# @param check.plot Check that plotting will finish in a reasonable amount of time
#' @param label Label the cell identies above the color bar
#' @param size Size of text above color bar
#' @param hjust Horizontal justification of text above color bar
#' @param angle Angle of text above color bar
#' @param combine Combine plots into a single gg object; note that if TRUE; themeing will not work when plotting multiple dimensions
#' @param ... Ignored for now
#'
#' @return A ggplot object
#'
#' @importFrom stats median
#' @importFrom scales hue_pal
#' @importFrom ggplot2 annotation_raster coord_cartesian ggplot_build aes_string
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
  disp.min = -2.5,
  disp.max = 2.5,
  slot = 'scale.data',
  assay = NULL,
  label = TRUE,
  size = 5.5,
  hjust = 0,
  angle = 45,
  combine = TRUE,
  ...
) {
  cells <- cells %||% colnames(x = object)
  if (is.numeric(x = cells)) {
    cells <- colnames(x = object)[cells]
  }
  assay <- assay %||% DefaultAssay(object = object)
  DefaultAssay(object = object) <- assay
  features <- features %||% VariableFeatures(object = object)
  features <- rev(x = unique(x = features))
  disp.max <- ifelse(
    test = slot != 'scale.data',
    yes = max(disp.max, 10),
    no = disp.max
  )
  data <- FetchData(
    object = object,
    vars = features,
    cells = cells,
    slot = slot
  )
  object <- StashIdent(object = object, save.name = 'ident')
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
    group.use <- groups.use[, i, drop = TRUE]
    group.use <- factor(x = group.use)
    names(x = group.use) <- cells
    plot <- SingleRasterMap(
      data = data,
      disp.min = disp.min,
      disp.max = disp.max,
      feature.order = features,
      cell.order = names(x = sort(x = group.use)),
      group.by = group.use
    )
    if (group.bar) {
      # TODO: Change group.bar to annotation.bar
      pbuild <- ggplot_build(plot = plot)
      cols <- hue_pal()(length(x = levels(x = group.use)))
      names(x = cols) <- levels(x = group.use)
      y.pos <- max(pbuild$layout$panel_params[[1]]$y.range) + 0.25
      y.max <- y.pos + 0.5
      plot <- plot + annotation_raster(
        raster = t(x = cols[sort(x = group.use)]),
        xmin = -Inf,
        xmax = Inf,
        ymin = y.pos,
        ymax = y.max
      ) +
        coord_cartesian(ylim = c(0, y.max), clip = 'off')
      if (label) {
        x.max <- max(pbuild$layout$panel_params[[1]]$x.range)
        x.divs <- pbuild$layout$panel_params[[1]]$x.major
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
    plots <- CombinePlots(plots = plots)
  }
  return(plots)
}

#' Hashtag oligo heatmap
#'
#' Draws a heatmap of hashtag oligo signals across singlets/doublets/negative cells. Allows for the visualization of HTO demultiplexing results.
#'
#' @param object Seurat object. Assumes that the hash tag oligo (HTO) data has been added and normalized, and demultiplexing has been run with HTODemux().
#' @param hto.classification The naming for object@meta.data slot with classification result from HTODemux().
#' @param global.classification The slot for object@meta.data slot specifying a cell as singlet/doublet/negative.
#' @param assay Hashtag assay name.
#' @param ncells Number of cells to plot. Default is to choose 5000 cells by random subsampling, to avoid having to draw exceptionally large heatmaps.
#' @param singlet.names Namings for the singlets. Default is to use the same names as HTOs.
#' @param ... Additional arguments for DoHeatmap().
#'
#' @return A ggplot object
#'
#' @importFrom ggplot2 guides
#' @export
#'
#' @examples
#' \dontrun{
#' object <- HTODemux(object)
#' HTOHeatmap(object)
#' }
#'
HTOHeatmap <- function(
  object,
  hto.classification = "hto_classification",
  global.classification = "hto_classification_global",
  assay = "HTO",
  ncells = 5000,
  singlet.names = NULL,
  ...
) {
  DefaultAssay(object = object) <- assay
  Idents(object = object) <- object[[hto.classification, drop = TRUE]]
  object <- subset(
    x = object,
    select = sample(x = colnames(x = object), size = ncells)
  )
  classification <- object[[hto.classification]]
  singlets <- which(x = object[[global.classification]] == 'Singlet')
  singlet.ids <- sort(x = unique(x = as.character(x = classification[singlets, ])))
  doublets <- which(object[[global.classification]] == 'Doublet')
  doublet.ids <- sort(x = unique(x = as.character(x = classification[doublets, ])))
  heatmap.levels <- c(singlet.ids, doublet.ids, 'Negative')
  if (length(x = doublets) > 0) {
    Idents(object = object, cells = doublets) <- 'Multiplet'
  }
  Idents(object = object) <- factor(
    x = Idents(object = object),
    levels = c(singlet.ids, 'Multiplet', 'Negative')
  )
  object <- ScaleData(object = object, assay = assay, verbose = FALSE)
  if (!is.null(x = singlet.names)) {
    levels(x = object@ident) <- c(singlet.names, "Multiplet", "Negative")
  }
  data <- FetchData(object = object, vars = singlet.ids)
  plot <- SingleRasterMap(
    data = data,
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
#' expression of the attribute being potted
#' @param assay Name of assay to use, defaults to the active assay
#' @param group.by Group (color) cells in different ways (for example, orig.ident)
#' @param y.max Maximum y axis value
#' @param same.y.lims Set all the y-axis limits to the same values
#' @param log plot Y axis on log scale
#' @param ncol Number of columns if multiple plots are displayed
#' @param combine Combine plots into a single gg object; note that if TRUE; themeing will not work when plotting multiple features
#' @param slot Use non-normalized counts data for plotting
#' @param ... Ignored
#'
#' @return A ggplot object
#'
#' @export
#'
#' @examples
#' RidgePlot(object = pbmc_small, features = 'PC1')
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
  combine = TRUE,
  slot = 'data',
  ...
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
    combine = combine,
    slot = slot,
    ...
  ))
}

#' Single cell violin plot
#'
#' Draws a violin plot of single cell data (gene expression, metrics, PC
#' scores, etc.)
#'
#' @inheritParams RidgePlot
#' @param pt.size Point size for geom_violin
#' @param split.by A variable to split the violin plots by
#' @param adjust Adjust parameter for geom_violin
#'
#' @return A ggplot object
#'
#' @export
#'
#' @examples
#' VlnPlot(object = pbmc_small, features = 'PC1')
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
  combine = TRUE,
  slot = 'data',
  ...
) {
  return(ExIPlot(
    object = object,
    type = 'violin',
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
    ...
  ))
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Dimensional reduction plots
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# @param do.hover Enable hovering over points to view information
# @param data.hover Data to add to the hover, pass a character vector of
# features to add. Defaults to cell name and ident. Pass 'NULL' to clear extra
# information.
# @param do.identify Opens a locator session to identify clusters of cells.
# @param vector.friendly FALSE by default. If TRUE, points are flattened into
# a PNG, while axes/labels retain full vector resolution. Useful for producing
# AI-friendly plots with large numbers of cells.
# @param png.file Used only if vector.friendly is TRUE. Location for temporary
# PNG file.
# @param png.arguments Used only if vector.friendly is TRUE. Vector of three
# elements (PNG width, PNG height, PNG DPI) to be used for temporary PNG.
# Default is c(10,10,100)
# @param ... Extra parameters to FeatureLocator for do.identify = TRUE

#' Dimensional reduction plot
#'
#' Graphs the output of a dimensional reduction technique (PCA by default).
#' Cells are colored by their identity class.
#'
#' @param object Seurat object
#' @param dims Dimensions to plot, must be a two-length numeric vector specifying x- and y-dimensions
#' @param cells Vector of cells to plot (default is all cells)
#' @param cols Vector of colors, each color corresponds to an identity class. By default, ggplot2 assigns colors
#' @param pt.size Adjust point size for plotting
#' @param reduction Which dimensionality reduction to use
#' @param group.by A vector of variables to group (color) cells by (for example, orig.ident);
#' pass 'ident' to group by identity class
#' @param shape.by If NULL, all points are circles (default). You can specify any
#' cell attribute (that can be pulled with FetchData) allowing for both
#' different colors and different shapes on cells
#' @param order Specify the order of plotting for the idents. This can be
#' useful for crowded plots if points of interest are being buried. Provide
#' either a full list of valid idents or a subset to be plotted last (on top)
#' @param label Whether to label the clusters
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
#' @param na.value Color value for NA points when using custom scale
#' @param combine Combine plots into a single gg object; note that if TRUE; themeing will not work when plotting multiple features
#' @param ... Ignored for now
#'
#' @return A ggplot object
#'
#' @export
#'
#' @note For the old \code{do.hover} and \code{do.identify} functionality, please see
#' \code{HoverLocator} and \code{FeatureLocator}, respectively.
#'
#' @seealso \code{\link{FeaturePlot}} \code{\link{HoverLocator}}
#' \code{\link{FeatureLocator}}
#'
#' @examples
#' DimPlot(object = pbmc_small)
#'
DimPlot <- function(
  object,
  dims = c(1, 2),
  cells = NULL,
  cols = NULL,
  pt.size = NULL,
  reduction = 'pca',
  group.by = NULL,
  shape.by = NULL,
  order = NULL,
  label = FALSE,
  label.size = 4,
  cells.highlight = NULL,
  cols.highlight = 'red',
  sizes.highlight = 1,
  na.value = 'grey50',
  combine = TRUE,
  ...
) {
  ReadPlotParams(object)
  if (length(x = dims) != 2) {
    stop("'dims' must be a two-length vector")
  }
  cells <- cells %||% colnames(x = object)
  data <- Embeddings(object = object[[reduction]])[cells, dims]
  data <- as.data.frame(x = data)
  dims <- paste0(Key(object = object[[reduction]]), dims)
  object <- StashIdent(object = object, save.name = 'ident')
  group.by <- group.by %||% 'ident'
  data[, group.by] <- object[[group.by]][cells, , drop = FALSE]
  if (!is.null(x = shape.by)) {
    data[, shape.by] <- object[[shape.by, drop = TRUE]]
  }
  plots <- lapply(
    X = group.by,
    FUN = function(x) {
      return(SingleDimPlot(
        data = data[, c(dims, x)],
        dims = dims,
        col.by = x,
        cols = cols,
        pt.size = pt.size,
        shape.by = shape.by,
        plot.order = order,
        label = label,
        label.size = label.size,
        cells.highlight = cells.highlight,
        cols.highlight = cols.highlight,
        sizes.highlight = sizes.highlight,
        na.value = na.value
      ))
    }
  )
  if (combine) {
    plots <- CombinePlots(plots = plots)
  }
  return(plots)
}

#' Visualize 'features' on a dimensional reduction plot
#'
#' Colors single cells on a dimensional reduction plot according to a 'feature'
#' (i.e. gene expression, PC scores, number of genes detected, etc.)
#'
#' @inheritParams DimPlot
#' @param features Vector of features to plot
#' @param cols The two colors to form the gradient over. Provide as string vector with
#' the first color corresponding to low values, the second to high. Also accepts a Brewer
#' color scale or vector of colors. Note: this will bin the data into number of colors provided.
#' @param min.cutoff Vector of minimum cutoff values for each feature, may specify quantile in the form of 'q##' where '##' is the quantile (eg, 1, 10)
#' @param max.cutoff Vector of maximum cutoff values for each feature, may specify quantile in the form of 'q##' where '##' is the quantile (eg, 1, 10)
#' @param split.by A factor in object metadata to split the feature plot by, pass 'ident' to split by cell identity'; similar to the old \code{FeatureHeatmap}
#' @param blend Scale and blend expression values to visualize coexpression of two features
#' @param blend.threshold The color cutoff from weak signal to strong signal; ranges from 0 to 1.
#' @param ncol Number of columns to combine multiple feature plots to, ignored if \code{split.by} is not \code{NULL}
#' @param combine Combine plots into a single gg object; note that if TRUE; themeing will not work when plotting multiple features
#' @param coord.fixed Plot cartesian coordinates with fixed aspect ratio
#'
#' @return A ggplot object
#'
#' @importFrom grDevices rgb
#' @importFrom cowplot theme_cowplot
#' @importFrom RColorBrewer brewer.pal.info
#' @importFrom ggplot2 labs scale_x_continuous scale_y_continuous theme element_rect dup_axis
#' element_blank element_text margin scale_color_brewer scale_color_gradientn scale_color_manual
#' @export
#'
#' @note For the old \code{do.hover} and \code{do.identify} functionality, please see
#' \code{HoverLocator} and \code{FeatureLocator}, respectively.
#'
#' @aliases FeatureHeatmap
#' @seealso \code{\link{DimPlot}} \code{\link{HoverLocator}}
#' \code{\link{FeatureLocator}}
#'
#' @examples
#' FeaturePlot(object = pbmc_small, features = 'PC1')
#'
FeaturePlot <- function(
  object,
  features,
  dims = c(1, 2),
  cells = NULL,
  cols = c("lightgrey",  "blue"),
  pt.size = NULL,
  min.cutoff = NA,
  max.cutoff = NA,
  reduction = 'tsne',
  split.by = NULL,
  shape.by = NULL,
  blend = FALSE,
  blend.threshold = 0.5,
  order = NULL,
  label = FALSE,
  label.size = 4,
  ncol = NULL,
  combine = TRUE,
  coord.fixed = FALSE,
  ...
) {
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
  if (length(x = dims) != 2 || !is.numeric(x = dims)) {
    stop("'dims' must be a two-length integer vector")
  }
  if (blend && length(x = features) != 2) {
    stop("Blending feature plots only works with two features")
  }
  dims <- paste0(Key(object = object[[reduction]]), dims)
  cells <- cells %||% colnames(x = object)
  data <- FetchData(object = object, vars = c(dims, features), cells = cells)
  features <- colnames(x = data)[3:ncol(x = data)]
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
  data[, 3:ncol(x = data)] <- sapply(
    X = 3:ncol(x = data),
    FUN = function(index) {
      data.feature <- as.vector(x = data[, index])
      min.use <- SetQuantile(cutoff = min.cutoff[index - 2], data.feature)
      max.use <- SetQuantile(cutoff = max.cutoff[index - 2], data.feature)
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
  colnames(x = data)[3:ncol(x = data)] <- features
  rownames(x = data) <- cells
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
  plots <- vector(
    mode = "list",
    length = ifelse(
      test = blend,
      yes = 4,
      no = length(x = features) * length(x = levels(x = data$split))
    )
  )
  xlims <- c(floor(x = min(data[, dims[1]])), ceiling(x = max(data[, dims[1]])))
  ylims <- c(floor(min(data[, dims[2]])), ceiling(x = max(data[, dims[2]])))
  if (blend) {
    ncol <- 4
    color.matrix <- BlendMatrix(col.threshold = blend.threshold)
    colors <- list(
      color.matrix[, 1],
      color.matrix[1, ],
      as.vector(x = color.matrix)
    )
  }
  for (i in 1:length(x = levels(x = data$split))) {
    ident <- levels(x = data$split)[i]
    data.plot <- data[as.character(x = data$split) == ident, , drop = FALSE]
    if (blend) {
      data.plot <- cbind(data.plot[, dims], BlendExpression(data = data.plot[, features[1:2]]))
      features <- colnames(x = data.plot)[3:ncol(x = data.plot)]
    }
    for (j in 1:length(x = features)) {
      feature <- features[j]
      if (blend) {
        cols.use <- as.numeric(x = as.character(x = data.plot[, feature])) + 1
        cols.use <- colors[[j]][sort(x = unique(x = cols.use))]
      } else {
        cols.use <- NULL
      }
      plot <- SingleDimPlot(
        data = data.plot[, c(dims, feature)],
        dims = dims,
        col.by = feature,
        pt.size = pt.size,
        cols = cols.use,
        do.label = label,
        label.size = label.size
      ) +
        scale_x_continuous(limits = xlims) +
        scale_y_continuous(limits = ylims) +
        theme_cowplot()
      if (length(x = levels(x = data$split)) > 1) {
        plot <- plot + theme(panel.border = element_rect(fill = NA, colour = 'black'))
        plot <- plot + if (i == 1) {
          labs(title = feature)
        } else {
          labs(title = NULL)
        }
        if (j == length(x = features) && !blend) {
          suppressMessages(
            expr = plot <- plot +
              scale_y_continuous(sec.axis = dup_axis(name = ident)) +
              no.right
          )
        }
        if (j != 1) {
          plot <- plot + theme(
            axis.line.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.text.y = element_blank(),
            axis.title.y.left = element_blank()
          )
        }
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
      if (!blend) {
        if (length(x = cols) == 1) {
          plot <- plot + scale_color_brewer(palette = cols)
        } else if (length(x = cols) > 1) {
          plot <- suppressMessages(
            expr = plot + scale_color_gradientn(
              colors = cols,
              guide = "colorbar"
            )
          )
        }
      }
      if (coord.fixed) {
        plot <- plot + coord_fixed()
      }
      plot <- plot
      plots[[(length(x = features) * (i - 1)) + j]] <- plot
    }
  }
  if (blend) {
    blend.legend <- BlendMap(color.matrix = color.matrix)
    for (i in 1:length(x = levels(x = data$split))) {
      suppressMessages(expr = plots <- append(
        x = plots,
        values = list(
          blend.legend +
            scale_y_continuous(
              sec.axis = dup_axis(name = levels(x = data$split)[i]),
              expand = c(0, 0)
            ) +
            labs(
              x = features[1],
              y = features[2],
              title = if (i == 1) {
                paste('Color threshold:', blend.threshold)
              } else {
                NULL
              }
            ) +
            no.right
        ),
        after = 4 * i - 1
      ))
    }
  }
  plots <- Filter(f = Negate(f = is.null), x = plots)
  if (combine) {
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
    plots <- CombinePlots(
      plots = plots,
      ncol = ncol,
      legend = legend,
      nrow = split.by %iff% length(x = levels(x = data$split))
    )
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
#'
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
  smooth = FALSE,
  ...
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
    ...
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
#' @param cols Colors to use for identity class plotting.
#' @param pt.size Size of the points on the plot
#' @param shape.by Ignored for now
#' @param span Spline span in loess function call, if \code{NULL}, no spline added
#' @param smooth Smooth the graph (similar to smoothScatter)
#' @param slot Slot to pull data from, should be one of 'counts', 'data', or 'scale.data'
#' @param ... Ignored for now
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
  cols = NULL,
  pt.size = 1,
  shape.by = NULL,
  span = NULL,
  smooth = FALSE,
  slot = 'data',
  ...
) {
  cells <- cells %||% colnames(x = object)
  plot <- SingleCorPlot(
    data = FetchData(
      object = object,
      vars = c(feature1, feature2),
      cells = cells,
      slot = slot
    ),
    col.by = Idents(object = object)[cells],
    cols = cols,
    pt.size = pt.size,
    smooth = smooth,
    legend.title = 'Identity'
  )
  if (!is.null(x = span)) {
    plot <- plot + geom_smooth(
      mapping = aes_string(x = feature1, y = feature2),
      method = 'loess',
      span = span
    )
  }
  return(plot)
}

#' View variable features
#'
#' @inheritParams FeatureScatter
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
  assay = NULL
) {
  if (length(x = cols) != 2) {
    stop("'cols' must be of length 2")
  }
  hvf.info <- HVFInfo(object = object, assay = assay)
  vars <- c(
    'mean',
    ifelse(
      test = 'variance.standardized' %in% colnames(x = hvf.info),
      yes = 'variance.standardized',
      no = 'dispersion'
    )
  )
  hvf.info <- hvf.info[, vars]
  log <- log %||% 'variance.standardized' %in% colnames(x = hvf.info)
  var.features <- VariableFeatures(object = object, assay = assay)
  var.status <- ifelse(
    test = rownames(x = hvf.info) %in% var.features,
    yes = 'yes',
    no = 'no'
  )
  plot <- SingleCorPlot(
    data = hvf.info,
    col.by = var.status,
    pt.size = pt.size
  )
  plot <- plot +
    labs(title = NULL, x = 'Average Expression', y = 'Dispersion') +
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
# Other plotting functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' ALRA Approximate Rank Selection Plot
#'
#' Plots the results of the approximate rank selection process for ALRA.
#'
#'
#' @param object Seurat object
#' @param start Index to start plotting singular value spacings from.
#' The transition from "signal" to "noise" in the is hard to see because the
#' first singular value spacings are so large. Nicer visualizations result from
#' skipping the first few. If set to 0 (default) starts from k/2.
#' @param combine Combine plots into a single gg object; note that if TRUE,
#' themeing will not work when plotting multiple features
#'
#' @return A list of 3 ggplot objects splotting the singular values, the
#' spacings of the singular values, and the p-values of the singular values.
#'
#' @author Jun Zhao, George Linderman
#' @seealso \code{\link{RunALRA}}
#'
#' @importFrom cowplot theme_cowplot
#' @importFrom ggplot2 ggplot aes_string geom_point geom_line
#' geom_vline scale_x_continuous labs
#' @export
#'
ALRAChooseKPlot <- function(object, start = 0, combine = TRUE) {
  if ( is.null(object@tools[["alra"]])) {
    stop('RunALRA() should be run prior to using this function.')
  }
  d <- object@tools[["alra"]][["d"]]
  diffs <- object@tools[["alra"]][["diffs"]]
  pvals <- object@tools[["alra"]][["pvals"]]
  k <- object@tools[["alra"]][["k"]]
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
  ggdata <- data.frame(x = 2:length(x = d), y = diffs)[-(1:(start - 1)), ]
  gg2 <- ggplot(data = ggdata, mapping = aes_string(x = 'x', y = 'y')) +
    geom_point(size = 1) +
    geom_line(size = 0.5) +
    geom_vline(xintercept = k + 1) +
    theme_cowplot() +
    scale_x_continuous(breaks = breaks) +
    labs(x = NULL, y = 's_{i} - s_{i-1}', title = 'Singular value spacings')
  ggdata <- data.frame(x = 2:length(x = d), y = pvals)
  gg3 <- ggplot(data = ggdata, mapping = aes_string(x = 'x', y = 'y')) +
    geom_point(size = 1) +
    geom_vline(xintercept = k + 1) +
    theme_cowplot() +
    scale_x_continuous(breaks = breaks) +
    labs(x = NULL, y = 'p.val', title = 'Singular value spacing p-values')
  plots <- list(spectrum = gg1, spacings = gg2, pvals = gg3)
  if (combine) {
    plots <- CombinePlots(plots = plots)
  }
  return(plots)
}

#' Dot plot visualization
#'
#' Intuitive way of visualizing how feature expression changes across different
#' identity classes (clusters). The size of the dot encodes the percentage of
#' cells within a class, while the color encodes the AverageExpression level of
#' cells within a class (blue is high).
#'
#' @param object Seurat object
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
#' @param split.by Factor to split the groups by (replicates the functionality of the old SplitDotPlotGG)
#' @param scale.by Scale the size of the points by 'size' or by 'radius'
#' @param scale.min Set lower limit for scaling, use NA for default
#' @param scale.max Set upper limit for scaling, use NA for default
#' @param ... Ignored
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
  scale.max = NA,
  ...
) {
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
  data.features$id <- as.vector(x = data.features$id)
  if (!is.null(x = split.by)) {
    splits <- object[[split.by, drop = TRUE]]
    if (length(x = unique(x = splits)) > length(x = cols)) {
      stop("Not enought colors for the number of groups")
    }
    cols <- cols[1:length(x = unique(x = splits))]
    names(x = cols) <- unique(x = splits)
    data.features$id <- paste(data.features$id, splits, sep = '_')
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
      X = strsplit(x = data.plot$id, split = '_'),
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
ElbowPlot <- function(
  object,
  ndims = 20,
  reduction = 'pca'
) {
  data.use <- Stdev(object = object, reduction = reduction)
  if (length(data.use) == 0) {
    stop(paste("No standard deviation info stored for", reduction))
  }
  if (ndims > length(x = data.use)) {
    warning("The object only has information for ", length(x = data.use), " reductions")
    ndims <- length(x = data.use)
  }
  stdev <- 'Standard Deviation'
  p <- ggplot(data = data.frame(dims = 1:ndims, stdev = data.use[1:ndims])) +
    geom_point(mapping = aes_string(x = 'dims', y = 'stdev')) +
    labs(x = Key(object = object[[reduction]]), y = stdev) +
    theme_cowplot()
  return(p)
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
  score.df <- score.df[dims, ]
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

#' Visualize Dimensional Reduction genes
#'
#' Visualize top genes associated with reduction components
#'
#' @param object Seurat object
#' @param reduction Reduction technique to visualize results for
#' @param dims Number of dimensions to display
#' @param nfeatures Number of genes to display
#' @param col Color of points to use
#' @param projected Use reduction values for full dataset (i.e. projected dimensional reduction values)
#' @param balanced Return an equal number of genes with + and - scores. If FALSE (default), returns the top genes ranked by the scores absolute values
#' @param ncol Number of columns to display
#' @param combine Combine plots into a single gg object; note that if TRUE; themeing will not work when plotting multiple features
#' @param ... Ignored
#'
#' @return A ggplot object
#'
#' @importFrom cowplot theme_cowplot
#' @importFrom ggplot2 ggplot aes_string geom_point labs
#' @export
#'
#' @examples
#' VizDimReduction(object = pbmc_small)
#'
VizDimReduction <- function(
  object,
  dims = 1:5,
  nfeatures = 30,
  col = 'blue',
  reduction = 'pca',
  projected = FALSE,
  balanced = FALSE,
  ncol = NULL,
  combine = TRUE,
  ...
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
    plots <- CombinePlots(plots = plots, ncol = ncol, legend = NULL)
  }
  return(plots)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Exported utility functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# NEED TO TREAT PNG PACKAGE PROPERLY
#' Augments ggplot2 scatterplot with a PNG image.
#'
#' Used in to creating vector friendly plots. Exported as it may be useful to others more broadly
#'
#' @param plot ggplot2 scatterplot. Typically will have only labeled axes and no points
#' @param img location of a PNG file that contains the points to overlay onto the scatterplot.
#'
#' @return ggplot2 scatterplot that includes the original axes but also the PNG file
#'
#' @importFrom png readPNG
#' @importFrom ggplot2 annotation_raster ggplot_build
#'
#' @export
#'
#' @examples
#' \dontrun{
#' data("pbmc_small")
#' p <- PCAPlot(pbmc_small, do.return = TRUE)
#' ggsave(filename = 'pcaplot.png', plot = p, device = png)
#' pmod <- AugmentPlot(plot = p, img = 'pcaplot.png')
#' pmod
#' }
#'
AugmentPlot <- function(plot, img) {
  range.values <- c(
    ggplot_build(plot = plot)$layout$panel_params[[1]]$x.range,
    ggplot_build(plot = plot)$layout$panel_params[[1]]$y.range
  )
  img <- readPNG(source = img)
  p1mod <- plot + annotation_raster(
    img,
    xmin = range.values[1],
    xmax = range.values[2],
    ymin = range.values[3],
    ymax = range.values[4]
  )
  return(p1mod)
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
#' plots <- FeaturePlot(
#'   object = pbmc_small,
#'   features = c('MS4A1', 'FCN1'),
#'   split.by = 'group',
#'   combine = FALSE
#' )
#' CombinePlots(
#'   plots = plots,
#'   legend = 'none',
#'   nrow = length(x = unique(x = pbmc_small[['group', drop = TRUE]]))
#' )
#'
CombinePlots <- function(plots, ncol = NULL, legend = NULL, ...) {
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

#' Feature Locator
#'
#' Select points on a scatterplot and get information about them
#'
#' @param plot A ggplot2 plot
#' @param ... Extra parameters, such as dark.theme, recolor, or smooth for using a dark theme,
#' recoloring based on selected cells, or using a smooth scatterplot, respectively
#'
#' @return The names of the points selected
#'
#' @importFrom ggplot2 ggplot_build
#' @export
#'
#' @seealso \code{\link[graphics]{locator}} \code{\link[ggplot2]{ggplot_build}}
#' \code{\link[SDMTools]{pnt.in.poly}} \code{\link{DimPlot}} \code{\link{FeaturePlot}}
#'
#' @examples
#' \dontrun{
#' plot <- DimPlot(object = pbmc_small)
#' # Follow instructions in the terminal to select points
#' cells.located <- FeatureLocator(plot = plot)
#' cells.located
#' }
#'
FeatureLocator <- function(plot, ...) {
  located <- PointLocator(plot = plot, ...)
  data <- ggplot_build(plot = plot)$plot$data
  selected <- data[as.numeric(x = rownames(x = located)), ]
  return(rownames(x = selected))
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
    title = names(x = data.frame())[1],
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
  plotly::layout(
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
    titlefont = title,
    paper_bgcolor = plotbg,
    plot_bgcolor = plotbg,
    ...
  )
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
#' @seealso \code{\link[ggplot2]{geom_text}}
#'
#' @examples
#' ff <- TopFeatures(object = pbmc_small[['pca']])
#' cc <- TopCells(object = pbmc_small[['pca']])
#' plot <- FeatureScatter(object = pbmc_small, feature1 = ff[1], feature2 = ff[2])
#' Labeler(plot = plot, points = cc)
#'
Labeler <- function(
  plot,
  points,
  labels = NULL,
  repel = FALSE,
  xnudge = 0.3,
  ynudge = 0.05,
  ...
) {
  geoms <- sapply(
    X = plot$layers,
    FUN = function(layer) {
      return(class(layer$geom)[1])
    }
  )
  geoms <- which(x = geoms == 'GeomPoint')
  if (length(x = geoms) == 0) {
    stop("Labelling only work on ggplot-based plots with a GeomPoint layer")
  }
  geoms <- min(geoms)
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
  plot$data$labels <- ''
  plot$data[points, 'labels'] <- labels
  x <- as.character(x = plot$mapping$x %||% plot$layers[[geoms]]$mapping$x)[2]
  y <- as.character(x = plot$mapping$y %||% plot$layers[[geoms]]$mapping$y)[2]
  geom.use <- ifelse(test = repel, yes = geom_text_repel, no = geom_text)
  if (repel) {
    if (!all(c(xnudge, ynudge) == 0)) {
      message("When using repel, set xnudge and ynudge to 0 for optimal results")
    }
    if (nrow(x = plot$data) > 1500) {
      warning(
        "repel is slow with a large number of points",
        call. = FALSE,
        immediate. = TRUE
      )
    }
  }
  plot <- plot + geom.use(
    mapping = aes_string(x = x, y = y, label = 'labels'),
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
#'   \item{\code{NoGrid}}{Removes grid lines}
#'   \item{\code{SeuratAxes}}{Set Seurat-style axes}
#'   \item{\code{BarePlot}}{Remove all extraneous features}
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
    #   Set text colors
    plot.title = white.text,
    plot.subtitle = white.text,
    axis.title = white.text,
    axis.text = white.text,
    legend.title = white.text,
    legend.text = white.text,
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

#' @export
#'
#' @rdname SeuratTheme
#' @aliases BarePlot
#'
BarePlot <- function() {
  return(NoLegend() + NoAxes() + NoGrid())
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
  color.heat$rows <- as.numeric(x = color.heat$rows)
  color.heat$cols <- as.numeric(x = color.heat$cols)
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
#
# @return An n x n matrix of blended colors
#
#' @importFrom grDevices rgb
#
BlendMatrix <- function(n = 10, col.threshold = 0.5) {
  if (0 > col.threshold || col.threshold > 1) {
    stop("col.threshold must be between 0 and 1")
  }
  return(outer(
    X = 1:n,
    Y = 1:n,
    FUN = function(i, j) {
      red <- 1 / (1 + exp(x = -(i - col.threshold * n) / 0.9))
      green <- 1 / (1 + exp(x = -(j - col.threshold * n) / 0.9))
      blue <- 0.2
      alpha <- lapply(X = list(red, green, blue), FUN = '^', 40)
      alpha <- Reduce(f = '+', x = alpha)
      alpha <- 0.99 - 0.1 * exp(x = -alpha / 1)
      return(rgb(red = red, green = green, blue = blue, alpha = alpha))
    }
  ))
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

globalVariables(
  names = 'cc',
  package = 'Seurat',
  add = TRUE
)
# Evaluate CCs
#
# Looks at the biweight midcorrelation of the Xth gene across the specified CCs
# for each group in the grouping.var.
#
# @param object A Seurat object
# @param grouping.var Grouping variable specified in alignment procedure
# @param dims.eval dimensions to evalutate the bicor for
# @param gene.num Xth gene to look at bicor for
# @param num.possible.genes Number of possible genes to search when choosing
# genes for the metagene. Set to 2000 by default. Lowering will decrease runtime
# but may result in metagenes constructed on fewer than num.genes genes.
# @param display.progress Show progress bar
#
#' @importFrom tidyr gather
#' @importFrom utils txtProgressBar setTxtProgressBar
#
EvaluateCCs <- function(
  object,
  grouping.var,
  dims.eval,
  gene.num,
  num.possible.genes,
  display.progress
) {
  reduction.type <-  "cca"
  ident.orig <- Idents(object = object)
  # object <- SetAllIdent(object = object, id = grouping.var)
  Idents(object = object) <- grouping.var
  levels.split <- names(x = sort(
    x = table(Idents(object = object)),
    decreasing = TRUE
  ))
  num.groups <- length(x = levels.split)
  objects <- list()
  for (i in 1:num.groups) {
    objects[[i]] <- subset(x = object, select = levels.split[i])
  }
  # object@ident <- ident.orig
  Idents(object = object) <- ident.orig
  cc.loadings <- list()
  scaled.data <- list()
  cc.embeds <- list()
  for (i in 1:num.groups) {
    cat(paste0("Rescaling group ", i, "\n"), file = stderr())
    objects[[i]] <- ScaleData(
      object = objects[[i]],
      block.size = 5000,
      display.progress = display.progress
    )
    objects[[i]] <- ProjectDim(
      object = objects[[i]],
      reduction = reduction.type,
      verbose = FALSE
    )
    cc.loadings[[i]] <- Loadings(
      object = objects[[i]][[reduction.type]],
      projected = TRUE
    )
    cc.embeds[[i]] <- Embeddings(object = objects[[i]][[reduction.type]])
    scaled.data[[i]] <- GetAssayData(object = objects[[i]], slot = 'scale.data')
  }
  bc.gene <- matrix(ncol = num.groups, nrow = length(dims.eval))
  if (display.progress) {
    cat(paste0("Evaluating dims: ", paste(dims.eval, collapse = " "),  "\n"), file = stderr())
    pb <- txtProgressBar(
      min = 0,
      max = length(x = dims.eval) * (num.groups - 1),
      style = 3
    )
    pb.idx <- 0
  }
  for (cc.use in dims.eval) {
    bc.gene.g1 <- c()
    for (g in 2:num.groups) {
      if (display.progress) {
        pb.idx <- pb.idx + 1
        setTxtProgressBar(pb = pb, value = pb.idx)
      }
      genes.rank <- data.frame(
        rank(x = abs(x = cc.loadings[[1]][, cc.use])),
        rank(x = abs(x = cc.loadings[[g]][, cc.use])),
        cc.loadings[[1]][, cc.use],
        cc.loadings[[g]][, cc.use]
      )
      genes.rank$min <- apply(X = genes.rank[,1:2], MARGIN = 1, FUN = min)
      genes.rank <- genes.rank[order(genes.rank$min, decreasing = TRUE), ]
      genes.top <- rownames(x = genes.rank)[1:min(num.possible.genes, nrow(x = genes.rank))]
      bicors <- list()
      for (i in c(1, g)) {
        cc.vals <- cc.embeds[[i]][, cc.use]
        bicors[[i]] <- sapply(
          X = genes.top,
          FUN = function(x) {
            return(BiweightMidcor(x = cc.vals, y = scaled.data[[i]][x, ]))
          }
        )
      }
      genes.rank <- data.frame(
        rank(x = abs(x = bicors[[1]])),
        rank(x = abs(x = bicors[[g]])),
        bicors[[1]],
        bicors[[g]]
      )
      genes.rank$min <- apply(X = abs(x = genes.rank[, 1:2]), MARGIN = 1, FUN = min)
      # genes must be correlated in same direction in both datasets
      genes.rank <- genes.rank[sign(x = genes.rank[,3]) == sign(x = genes.rank[,4]), ]
      genes.rank <- genes.rank[order(genes.rank$min, decreasing = TRUE), ]
      genes.rank <- genes.rank[order(genes.rank$min, decreasing = TRUE), ]
      bc.gene[cc.use, g] <- mean(x = abs(x = genes.rank[1:gene.num, 4]))
      bc.gene.g1 <- c(bc.gene.g1, mean(x = abs(x = genes.rank[1:gene.num, 3])))
    }
    bc.gene[cc.use, 1] <- abs(x = mean(x = bc.gene.g1))
  }
  if (display.progress) {
    close(con = pb)
  }
  colnames(x = bc.gene) <- levels.split
  bc.gene <- as.data.frame(x = bc.gene)
  bc.gene$cc <- 1:nrow(x = bc.gene)
  bc.gene <- gather(data = bc.gene, key = "Group",  value = "bicor", -cc)
  return(bc.gene)
}

# Plot feature expression by identity
#
# Basically combines the codebase for VlnPlot and RidgePlot
#
# @param object Seurat object
# @param plot.type Plot type, choose from 'ridge' or 'violin'
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
# @param combine Combine plots using cowplot::plot_grid
# @param slot Use non-normalized counts data for plotting
# @param ... Ignored
#
#' @importFrom scales hue_pal
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
  combine = TRUE,
  slot = 'data',
  ...
) {
  assay <- assay %||% DefaultAssay(object = object)
  DefaultAssay(object = object) <- assay
  ncol <- ncol %||% ifelse(
    test = length(x = features) > 9,
    yes = 4,
    no = min(length(x = features), 3)
  )
  data <- FetchData(object = object, vars = features, slot = slot)
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
    } else if (cols == 'interaction') {
      split <- interaction(idents, split)
      cols <- hue_pal()(length(x = levels(x = idents)))
    } else {
      cols <- Col2Hex(cols)
    }
    cols <- Interleave(cols, InvertHex(hexadecimal = cols))
    names(x = cols) <- sort(x = levels(x = split))
  }
  if (same.y.lims && is.null(x = y.max)) {
    y.max <- max(data)
  }
  plots <- lapply(
    X = features,
    FUN = function(x) {
      return(SingleExIPlot(
        feature = x,
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
  if (combine) {
    plots <- CombinePlots(plots = plots, ncol = ncol, legend = 'none')
  }
  return(plots)
}

# Documentation needed
#Sets the parameter value for a plotting function as default for an object
# @param object ...
# @param function.name ...
# @param param.name ...
# @param param.value ...
#
# Quick example :   pbmc <- FixPlotParam(pbmc,"FeaturePlot","pt.size","0.1")
#
FixPlotParam <- function(object, function.name, param.name, param.value) {
  if ("PlotParams" %in% names(object@misc)) {
    object@misc[["PlotParams"]][paste(function.name, param.name, sep = ":")] <- param.value
  }
  else {
    plotParamsList <- list()
    plotParamsList[paste(function.name, param.name, sep = ":")] <- param.value
    object@misc[["PlotParams"]] <- plotParamsList
  }
  return(object)
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
    # browser()
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
    X = hexadecimal,
    FUN = function(hex) {
      hex <- unlist(x = strsplit(
        x = gsub(pattern = '#', replacement = '', x = hex),
        split = ''
      ))
      key <- unlist(x = strsplit(x = 'FEDCBA9876543210', split = ''))
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
#' @importFrom SDMTools pnt.in.poly
#
PointLocator <- function(plot, recolor = TRUE, dark.theme = FALSE, ...) {
  #   Convert the ggplot object to a data.frame
  plot.data <- GGpointToBase(plot = plot, dark.theme = dark.theme, ...)
  npoints <- nrow(x = plot.data)
  cat("Click around the cluster of points you wish to select\n")
  cat("ie. select the vertecies of a shape around the cluster you\n")
  cat("are interested in. Press <Esc> when finished (right click for R-terminal users)\n\n")
  polygon <- locator(n = npoints, type = 'l')
  polygon <- data.frame(polygon)
  #   pnt.in.poly returns a data.frame of points
  points.all <- pnt.in.poly(
    pnts = plot.data[, c(1, 2)],
    poly.pnts = polygon
  )
  #   Find the located points
  points.located <- points.all[which(x = points.all$pip == 1), ]
  #   If we're recoloring, do the recolor
  if (recolor) {
    no <- ifelse(test = dark.theme, yes = 'white', no = 'black')
    points.all$color <- ifelse(test = points.all$pip == 1, yes = 'red', no = no)
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

# Documentation needed
#
# @param object ...
# @param workflow.name ...
# @param depth ...
#
ReadPlotParams <- function(object, workflow.name, depth = 1) {
  if (!("PlotParams" %in% names(object@misc))) {
    return()
  }
  command.name <- as.character(deparse(sys.calls()[[sys.nframe()-depth]]))
  command.name <- gsub(pattern = ".Seurat",replacement = "",x = command.name)
  command.name <- ExtractField(string = command.name, field = 1, delim = "\\(")
  param.list <- names(formals(fun = sys.function(sys.parent(depth))))
  paramsList <- slot(object = object, name = "misc")[["PlotParams"]]
  p.env <- parent.frame(depth)
  for(i in names(paramsList)) {
    split_arr <- unlist(strsplit(x = i,split = ":"))
    function.name <- split_arr[1]
    param.name <- split_arr[2]
    param.value <- paramsList[[i]]
    if ((function.name == command.name) && (param.name %in% param.list)) {
      assign(
        x = param.name,
        value = param.value,
        envir = p.env
      )
    }
  }
  # overwrite any arguments passed in on the command line
  argnames <- sys.call(which = depth)
  argList <- as.list(argnames[-1])
  args_ignore <- c("", "object", "workflow.name")
  args_use <- setdiff(x = names(argList), y = args_ignore)
  for(i in args_use) {
    if(as.character(unlist(argList[i])[[1]]) == "F") {
      arg.val <- FALSE
    } else if(as.character(unlist(argList[i])[[1]]) == "T") {
      arg.val <- TRUE
    } else {
      arg.val <- unlist(argList[i])[[1]]
    }
    assign(x = i, value = arg.val, envir = p.env)
  }
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
  ...
) {
  pt.size <- pt.size <- pt.size %||% min(1583 / nrow(x = data), 1)
  orig.names <- colnames(x = data)
  names.plot <- colnames(x = data) <- gsub(
    pattern = '-',
    replacement = '.',
    x = colnames(x = data),
    fixed = TRUE
  )
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
    plot <- plot + stat_density2d(
      mapping = aes(fill = ..density.. ^ 0.25),
      geom = 'tile',
      contour = FALSE,
      n = 200,
      h = Bandwidth(data = data[, names.plot])
    ) +
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
  plot <- plot + theme_cowplot()
  return(plot)
}

# Plot a single dimension
#
# @param data Data to plot
# @param dims A two-length numeric vector with dimensions to use
# @param pt.size Adjust point size for plotting
# @param col.by ...
# @param cols Vector of colors, each color corresponds to an identity class. By default, ggplot assigns colors.
# @param shape.by If NULL, all points are circles (default). You can specify any cell attribute (that can be pulled with FetchData) allowing for both different colors and different shapes on cells.
# @param order Specify the order of plotting for the idents. This can be useful for crowded plots if points of interest are being buried. Provide either a full list of valid idents or a subset to be plotted last (on top).
# @param label Whether to label the clusters
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
# @param ... Ignored for now
#
#' @importFrom ggplot2 ggplot aes_string labs geom_text guides
#' scale_color_brewer scale_color_manual element_rect guide_legend
#' @importFrom cowplot theme_cowplot
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
  label.size = 4,
  cells.highlight = NULL,
  cols.highlight = 'red',
  sizes.highlight = 1,
  na.value = 'grey50',
  ...
) {
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
      sizes.highlight = sizes.highlight,
      cols.highlight = cols.highlight,
      col.base = cols[1] %||% 'black',
      pt.size = pt.size
    )
    order <- highlight.info$plot.order
    data$highlight <- highlight.info$highlight
    col.by <- 'highlight'
    pt.size <- highlight.info$size
    cols <- highlight.info$color
  }
  if (!is.null(x = order) && !is.null(x = col.by)) {
    order <- rev(x = c(
      order,
      setdiff(x = unique(x = data[, col.by]), y = order)
    ))
    data[, col.by] <- factor(x = data[, col.by], levels = order)
  }
  if (!is.null(x = col.by) && !col.by %in% colnames(x = data)) {
    warning("Cannot find ", col.by, " in plotting data, not coloring plot")
    col.by <- NULL
  } else {
    col.index <- grep(pattern = col.by, x = colnames(x = data), fixed = TRUE)
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
  pt.size <- pt.size %||% min(1583 / nrow(x = data), 1)
  plot <- ggplot(data = data) +
    geom_point(
      mapping = aes_string(
        x = dims[1],
        y = dims[2],
        color = col.by,
        shape = shape.by
      ),
      size = pt.size
    ) +
    guides(color = guide_legend(override.aes = list(size = 3))) +
    labs(color = NULL)
  if (label && !is.null(x = col.by)) {
    labels <- MakeLabels(data = plot$data[, c(dims, col.by)])
    plot <- plot +
      geom_point(
        data = labels,
        mapping = aes_string(
          x = colnames(x = labels)[1],
          y = colnames(x = labels)[2]
        ),
        size = 0,
        alpha = 0
      ) +
      geom_text(
        data = labels,
        mapping = aes_string(
          x = colnames(x = labels)[1],
          y = colnames(x = labels)[2],
          label = colnames(x = labels)[3]
        ),
        size = label.size
      )
  }
  if (!is.null(x = cols)) {
    plot <- plot + if (length(x = cols) == 1) {
      scale_color_brewer(palette = cols, na.value = na.value)
    } else {
      scale_color_manual(values = cols, na.value = na.value)
    }
  }
  plot <- plot + theme_cowplot()
  return(plot)
}

# Plot a single expression by identity on a plot
#
# @param feature Feature to plot
# @param type Make either a 'ridge' or 'violin' plot
# @param data Data to plot
# @param idents Idents to use
# @param sort Sort identity classes (on the x-axis) by the average
# expression of the attribute being potted
# @param y.max Maximum Y value to plot
# @param adjust Adjust parameter for geom_violin
# @param cols Colors to use for plotting
# @param feature.names
# @param log plot Y axis on log scale
# @param ... Ignored
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
  log = FALSE,
  ...
) {
  set.seed(seed = 42)
  if (!is.data.frame(x = data) || ncol(x = data) != 1) {
    stop("'SingleExIPlot requires a data frame with 1 column")
  }
  feature <- colnames(x = data)
  data$ident <- idents
  if (sort) {
    data$ident <- factor(
      x = data$ident,
      levels = names(x = rev(x = sort(x = tapply(
        X = data[, feature],
        INDEX = data$ident,
        FUN = mean
      ))))
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
  axis.label <- ifelse(test = log, yes = 'Log Expression Level', no = 'Expression Level')
  y.max <- y.max %||% max(data[, feature])
  if (is.null(x = split) || type != 'violin') {
    vln.geom <- geom_violin
    fill <- 'ident'
  } else {
    data$split <- split
    vln.geom <- geom_split_violin
    fill <- 'split'
  }
  switch(
    EXPR = type,
    'violin' = {
      x <- 'ident'
      y <- feature
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
      x <- feature
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
    theme_cowplot()
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
      labels <- unique(x = as.vector(x = data$ident))
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

# A single heatmap from ggplot2 using geom_raster
#
# @param data A matrix or data frame with data to plot
# @param cell.order ...
# @param feature.order ...
# @param cols A vector of colors to use
# @param disp.min Minimum display value (all values below are clipped)
# @param disp.max Maximum display value (all values above are clipped)
# @param limits A two-length numeric vector with the limits for colors on the plot
# @param group.by A vector to group cells by, should be one grouping identity per cell
# @param ... Ignored
#
#' @importFrom ggplot2 ggplot aes_string geom_raster scale_fill_gradient
#' scale_fill_gradientn theme element_blank labs geom_point guides guide_legend
#
SingleRasterMap <- function(
  data,
  cell.order = NULL,
  feature.order = NULL,
  colors = PurpleAndYellow(),
  disp.min = -2.5,
  disp.max = 2.5,
  limits = NULL,
  group.by = NULL,
  ...
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
  plot <- ggplot(data = data) +
    geom_raster(mapping = aes_string(x = 'Cell', y = 'Feature', fill = 'Expression')) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
    scale_fill_gradientn(limits = limits, colors = colors) +
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
