#' Gene expression heatmap
#'
#' Draws a heatmap of single cell gene expression.
#'
#' @param object Seurat object
#' @param data.use Option to pass in data to use in the heatmap. Default will pick from either 
#' object@@data or object@@scale.data depending on use.scaled parameter. Should have cells as columns
#' and genes as rows. 
#' @param use.scaled Whether to use the data or scaled data if data.use is NULL 
#' @param cells.use Cells to include in the heatmap (default is all cells)
#' @param genes.use Genes to include in the heatmap (ordered)
#' @param disp.min Minimum display value (all values below are clipped)
#' @param disp.max Maximum display value (all values above are clipped)
#' @param group.by Groups cells by this variable. Default is object@@ident
#' @param draw.line Draw vertical lines delineating different groups
#' @param col.low Color for lowest expression value
#' @param col.mid Color for mid expression value
#' @param col.high Color for highest expression value
#' @param slim.col.label display only the identity class name once for each group
#' @param remove.key Removes the color key from the plot.
#' @param cex.col Controls size of column labels (cells)
#' @param cex.row Controls size of row labels (genes)
#' @param group.label.loc Place group labels on bottom or top of plot.
#' @param group.label.rot Whether to rotate the group label.
#' @param group.cex Size of group label text
#' @param group.spacing Controls amount of space between columns.
#' @param do.plot Whether to display the plot. 
#' @return Returns a ggplot2 plot object
#' @importFrom reshape2 melt
#' @importFrom dplyr %>%
#' @export
DoHeatmapGG <- function(object, data.use = NULL, use.scaled = TRUE, cells.use = NULL, genes.use = NULL, 
                        disp.min = -2.5,  disp.max = 2.5, group.by = "ident", draw.line = TRUE, 
                        col.low = "#FF00FF", col.mid = "#000000", col.high = "#FFFF00", 
                        slim.col.label = FALSE, remove.key = FALSE, title = NULL, cex.col = 10, cex.row = 10,
                        group.label.loc = "bottom", group.label.rot = FALSE, group.cex = 15, 
                        group.spacing = 0.15, do.plot = TRUE, ...) {
  
  if (is.null(data.use)){
    if (use.scaled){
      data.use <- object@scale.data 
    }
    else{
      data.use <- object@data
    }
  }
  
  # note: data.use should have cells as column names, genes as row names
  cells.use <- set.ifnull(cells.use, object@cell.names)
  cells.use <- ainb(cells.use, colnames(data.use))
  if (length(cells.use) == 0) {
    stop("No cells given to cells.use present in object")
  }
  genes.use <- set.ifnull(genes.use, rownames(data.use))
  genes.use <- ainb(genes.use, rownames(data.use))
  if (length(genes.use) == 0) {
    stop("No genes given to genes.use present in object")
  }
  if (is.null(group.by) || group.by == "ident") {
    cells.ident <- object@ident[cells.use]
  }
  else {
    cells.ident <- factor(FetchData(object, cells.use = cells.use, vars.all = group.by)[, 1])
    names(cells.ident) <- cells.use
  }
  cells.ident <- factor(cells.ident, labels = ainb(levels(cells.ident), cells.ident))
  
  data.use <- data.use[genes.use, cells.use]
  if(use.scaled){
    data.use <- minmax(data.use, min = disp.min, max = disp.max)
  }
  
  data.use <- as.data.frame(t(data.use))
  data.use$cell <- rownames(data.use)
  
  colnames(data.use) <- make.unique(colnames(data.use))
  
  data.use %>% melt(id.vars = "cell") -> data.use
  names(data.use)[names(data.use) == 'variable'] <- 'gene'
  names(data.use)[names(data.use) == 'value'] <- 'expression'
  
  data.use$ident <- cells.ident[data.use$cell]
  breaks <- seq(min(data.use$expression), max(data.use$expression), length = length(pyCols)+1)
  data.use$gene <- with(data.use, factor(gene, levels = rev(unique(data.use$gene))))
  data.use$cell <- with(data.use, factor(cell, levels = cells.use))
  
  
  # might be a solution if we want discrete interval units, makes the legend clunky though
  #data.use$expression <- cut(data.use$expression, breaks = breaks, include.lowest = T)
  #heatmap <- ggplot(data.use, aes(x = cell, y = gene, fill = expression)) + geom_tile() +
  #                  scale_fill_manual(values = pyCols, name= "Expression") + 
  #                  scale_y_discrete(position = "right", labels = rev(genes.use)) + 
  #                  theme(axis.line=element_blank(), axis.title.y=element_blank(), 
  #                        axis.ticks.y = element_blank())
  
  heatmap <- ggplot(data.use, aes(x = cell, y = gene, fill = expression)) + geom_tile() +
    scale_fill_gradient2(low = col.low, mid = col.mid, high = col.high, name= "Expression") + 
    scale_y_discrete(position = "right", labels = rev(genes.use)) + 
    theme(axis.line=element_blank(), axis.title.y=element_blank(), 
          axis.ticks.y = element_blank(), strip.text.x = element_text(size = group.cex), 
          axis.text.y = element_text(size = cex.row), axis.text.x = element_text(size = cex.col),
          axis.title.x=element_blank())
  
  if (slim.col.label){
    heatmap <- heatmap + theme(axis.title.x=element_blank(), axis.text.x=element_blank(),
                               axis.ticks.x=element_blank(), axis.line=element_blank(), 
                               axis.title.y=element_blank(), axis.ticks.y = element_blank())
  }
  else{
    heatmap <- heatmap + theme(axis.text.x = element_text(angle = 90))
  }
  if(!is.null(group.by)){
    if(group.label.loc == "top"){
      heatmap <- heatmap + facet_grid(~ident,drop=T,space="free",scales="free") + scale_x_discrete(expand=c(0,0),drop=T)
    }
    else{
      heatmap <- heatmap + facet_grid(~ident,drop=T,space="free",scales="free", switch="x") + scale_x_discrete(expand=c(0,0),drop=T)
    }
    if(draw.line){
      heatmap <- heatmap + theme(strip.background = element_blank(), panel.spacing = unit(group.spacing, "lines"))
    }
    else{
      heatmap <- heatmap + theme(strip.background = element_blank(), panel.spacing = unit(0, "lines"))
    }
    if(group.label.rot){
      heatmap <- heatmap + theme(strip.text.x = element_text(angle = 90))
    }
  }
  if (remove.key){
    heatmap <- heatmap + theme(legend.position="none")
  }
  if (!is.null(title)){
    heatmap <- heatmap + labs(title = title)
  }
  if(do.plot) {
    heatmap
  }
  return(heatmap)
}


#' Single cell violin plot
#'
#' Draws a violin plot of single cell data (gene expression, metrics, PC
#' scores, etc.)
#'
#' @param object Seurat object
#' @param features.plot Features to plot (gene expression, metrics, PC scores,
#' anything that can be retreived by FetchData)
#' @param nCol Number of columns if multiple plots are displayed
#' @param do.ret FALSE by default. If TRUE, returns a list of ggplot objects.
#' @param do.sort Sort identity classes (on the x-axis) by the average
#' expression of the attribute being potted
#' @param y.max Maximum y axis value
#' @param same.y.lims Set all the y-axis limits to the same values
#' @param size.x.use X axis title font size
#' @param size.y.use Y axis title font size
#' @param size.title.use Main title font size
#' @param adjust.use Adjust parameter for geom_violin
#' @param point.size.use Point size for geom_violin
#' @param cols.use Colors to use for plotting
#' @param group.by Group (color) cells in different ways (for example, orig.ident)
#' @param y.log plot Y axis on log scale
#' @param x.lab.rot Rotate x-axis labels
#' @param y.lab.rot Rotate y-axis labels
#' @param legend.position Position the legend for the plot
#' @param single.legend Consolidate legend the legend for all plots
#' @param remove.legend Remove the legend from the plot
#' @param return.plotlist Return the list of individual plots instead of compiled plot.
#' @param \dots additional parameters to pass to FetchData (for example, use.imputed, use.scaled, use.raw)
#' @import ggplot2
#' @importFrom cowplot plot_grid
#' @return By default, no return, only graphical output. If do.return=TRUE,
#' returns a list of ggplot objects.
#' @export
VlnPlot <- function(object, features.plot, nCol = NULL, do.ret = FALSE, do.sort = FALSE, y.max = NULL,
                    same.y.lims = F, size.x.use = 16, size.y.use = 16, size.title.use = 20, 
                    adjust.use = 1, point.size.use = 1, cols.use = NULL, group.by = NULL, y.log = F, 
                    x.lab.rot = FALSE, y.lab.rot = FALSE, legend.position = "right", 
                    single.legend = TRUE, remove.legend = FALSE, return.plotlist = FALSE, ...){
            
            if (is.null(nCol)) {
              nCol <- min(length(features.plot), 3)
              if (length(features.plot) > 9) nCol <- 4
            }

            data.use <- data.frame(FetchData(object, features.plot, ...))
            ident.use <- object@ident
            if (!is.null(group.by)) ident.use <- as.factor(FetchData(object, group.by)[, 1])
            gene.names <- colnames(data.use)[colnames(data.use) %in% rownames(object@data)]
            if(single.legend) remove.legend <- TRUE
            if(same.y.lims && is.null(y.max)) y.max <- max(data.use)
            plots <- lapply(features.plot, function(x) PlotVln(feature = x, 
                                                               data = data.use[, x, drop = FALSE], 
                                                               cell.ident = ident.use, 
                                                               do.sort = do.sort, y.max = y.max,
                                                               size.x.use = size.x.use, 
                                                               size.y.use = size.y.use, 
                                                               size.title.use = size.title.use,
                                                               adjust.use = adjust.use, 
                                                               point.size.use = point.size.use,
                                                               cols.use = cols.use, 
                                                               gene.names = gene.names, y.log = y.log,
                                                               x.lab.rot = x.lab.rot, 
                                                               y.lab.rot = y.lab.rot, 
                                                               legend.position = legend.position,
                                                               remove.legend = remove.legend))
            plots.combined <- plot_grid(plotlist = plots,  ncol = nCol)
            if(single.legend){
              legend <- legend <- get_legend(plots[[1]] + theme(legend.position = legend.position))
              if(legend.position == "bottom") plots.combined <- plot_grid(plots.combined, legend, ncol = 1, rel_heights = c(1, .2))
              if(legend.position == "right") plots.combined <- plot_grid(plots.combined, legend, rel_widths = c(3, .3))
              else warning("Shared legends must be at the bottom or right of the plot")
            }
            plots.combined
            if(return.plotlist){
              return(plots)
            }
            else{
              return(plots.combined)
            }
}


PlotVln <- function(feature, data, cell.ident, do.sort, y.max, size.x.use, size.y.use, size.title.use, 
                    adjust.use, point.size.use, cols.use, gene.names, y.log, x.lab.rot, y.lab.rot,
                    legend.position, remove.legend) {
  set.seed(42)
  data$ident <- cell.ident
  if(do.sort) {
    data$ident <- factor(data$ident, levels = names(rev(sort(tapply(data[, feature], data$ident, mean)))))
  }
  if (y.log){
    noise <- rnorm(length(data[, feature])) / 200
    data.melt[, feature] <- data[, feature] + 1
  }
  else{
    noise <- rnorm(length(data[, feature])) / 100000
  }
  data[, feature] <- data[, feature] + noise
  y.max <- set.ifnull(y.max, max(data[, feature]))
  plot <- ggplot(data, aes(factor(ident), eval(parse(text = feature)))) + 
    geom_violin(scale = "width", adjust = adjust.use, trim = TRUE, aes(fill = factor(ident))) + 
    theme(legend.position = legend.position, axis.title.x = element_text(face="bold", colour="#990000", size = size.x.use),
          axis.title.y = element_text(face = "bold", colour = "#990000", size = size.y.use)) +
    guides(fill = guide_legend(title = NULL)) + 
    geom_jitter(height = 0, size = point.size.use) + xlab("Cell Type") + nogrid +
    ggtitle(feature) + theme(plot.title = element_text(size = size.title.use, face = "bold")) +
    ylim(min(data[, feature]), y.max)
  if (y.log) plot <- plot + scale_y_log10()
  if(feature %in% gene.names){
    if(y.log){
      plot <- plot + ylab("Log Expression level")
    }
    else {
      plot <- plot + ylab("Expression level")
    }
  }
  else{
    plot <- plot + ylab("")
  }
  if (!is.null(cols.use)) {
    plot <- plot + scale_fill_manual(values = cols.use)
  }
  if(x.lab.rot) plot <- plot + theme(axis.text.x = element_text(angle = 90, vjust=0.5))
  if(y.lab.rot) plot <- plot + theme(axis.text.x = element_text(angle = 90))
  if(remove.legend) plot <- plot + theme(legend.position="none")
  return(plot)
}
