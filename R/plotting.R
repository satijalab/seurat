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

