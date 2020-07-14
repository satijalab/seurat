#' @include objects.R
#' @include generics.R
#' @include visualization.R
#' @importFrom methods setClass setOldClass slot<- setAs setMethod new
#'
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Class definitions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @note \code{scalefactors} objects can be created with \code{scalefactors()}
#'
#' @param spot Spot full resolution scale factor
#' @param fiducial Fiducial full resolution scale factor
#' @param hires High resolutoin scale factor
#' @param lowres Low resolution scale factor
#'
#' @rdname ScaleFactors
#' @export
#'
scalefactors <- function(spot, fiducial, hires, lowres) {
  object <- list(
    spot = spot,
    fiducial = fiducial,
    hires = hires,
    lowres = lowres
  )
  object <- sapply(X = object, FUN = as.numeric, simplify = FALSE, USE.NAMES = TRUE)
  return(structure(.Data = object, class = 'scalefactors'))
}

setOldClass(Classes = c('scalefactors'))

#' The SpatialImage class
#'
#' The SpatialImage class is a virtual class representing spatial information for
#' Seurat. All spatial image information must inherit from this class for use with
#' \code{Seurat} objects
#'
#' @slot assay Name of assay to associate image data with; will give this image
#' priority for visualization when the assay is set as the active/default assay
#' in a \code{Seurat} object
#' @slot key Key for the image
#'
#' @section Provided methods:
#' These methods are defined on the \code{SpatialImage} object and should not be
#' overwritten without careful thought
#' \itemize{
#'   \item \code{\link{DefaultAssay}} and \code{\link{DefaultAssay<-}}
#'   \item \code{\link{Key}} and \code{\link{Key<-}}
#'   \item \code{\link{IsGlobal}}
#'   \item \code{\link{Radius}}; this method \emph{can} be overridden to provide
#'   a spot radius for image objects
#' }
#'
#' @section Required methods:
#' All subclasses of the \code{SpatialImage} class must define the following methods;
#' simply relying on the \code{SpatialImage} method will result in errors. For required
#' parameters and their values, see the \code{Usage} and \code{Arguments} sections
#' \describe{
#'   \item{\code{\link{Cells}}}{Return the cell/spot barcodes associated with each position}
#'   \item{\code{\link{dim}}}{Return the dimensions of the image for plotting in \code{(Y, X)} format}
#'   \item{\code{\link{GetImage}}}{Return image data; by default, must return a grob object}
#'   \item{\code{\link{GetTissueCoordinates}}}{Return tissue coordinates; by default,
#'   must return a two-column data.frame with x-coordinates in the first column and y-coordiantes
#'   in the second}
#'   \item{\code{\link{Radius}}}{Return the spot radius; returns \code{NULL} by
#'   default for use with non-spot image technologies}
#'   \item{\code{\link{RenameCells}}}{Rename the cell/spot barcodes for this image}
#'   \item{\code{\link{subset}} and \code{[}}{Subset the image data by cells/spots;
#'   \code{[} should only take \code{i} for subsetting by cells/spots}
#' }
#' These methods are used throughout Seurat, so defining them and setting the proper
#' defaults will allow subclasses of \code{SpatialImage} to work seamlessly
#'
#' @name SpatialImage-class
#' @rdname SpatialImage-class
#' @exportClass SpatialImage
#'
SpatialImage <- setClass(
  Class = 'SpatialImage',
  contains = 'VIRTUAL',
  slots = list(
    'assay' = 'character',
    'key' = 'character'
  )
)

#' The SlideSeq class
#'
#' The SlideSeq class represents spatial information from the Slide-seq platform
#'
#' @inheritSection SpatialImage Slots
#' @slot coordinates ...
#' @slot ...
#'
SlideSeq <- setClass(
  Class = 'SlideSeq',
  contains = 'SpatialImage',
  slots = list(
    'coordinates' = 'data.frame'
  )
)

#' The STARmap class
#'
#' The STARmap class represents spatial information from the STARmap platform
#'
#' @inheritSection SpatialImage Slots
#' @slot ...
#'
STARmap <- setClass(
  Class = 'STARmap',
  contains = 'SpatialImage',
  slots = list(
    'coordinates' = 'data.frame',
    'qhulls' = 'data.frame'
  )
)

#' The VisiumV1 class
#'
#' The VisiumV1 class represents spatial information from the 10X Genomics Visium
#' platform
#'
#' @slot image A three-dimensional array with PNG image data, see
#' \code{\link[png]{readPNG}} for more details
#' @slot scale.factors An object of class \code{\link{scalefactors}}; see
#' \code{\link{scalefactors}} for more information
#' @slot coordinates A data frame with tissue coordinate information
#' @slot spot.radius Single numeric value giving the radius of the spots
#'
#' @name VisiumV1-class
#' @rdname VisiumV1-class
#' @exportClass VisiumV1
#'
VisiumV1 <- setClass(
  Class = 'VisiumV1',
  contains = 'SpatialImage',
  slots = list(
    'image' = 'array',
    'scale.factors' = 'scalefactors',
    'coordinates' = 'data.frame',
    'spot.radius' = 'numeric'
  )
)

setClass(Class = 'SliceImage', contains = 'VisiumV1')

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Get a vector of cell names associated with an image (or set of images)
#'
#' @param object Seurat object
#' @param images Vector of image names
#' @param unlist Return as a single vector of cell names as opposed to a list,
#' named by image name.
#'
#' @return A vector of cell names
#'
#' @examples
#' \dontrun{
#' CellsByImage(object = object, images = "slice1")
#' }
#'
CellsByImage <- function(object, images = NULL, unlist = FALSE) {
  images <- images %||% Images(object = object)
  cells <- sapply(
    X = images,
    FUN = function(x) {
      Cells(x = object[[x]])
    },
    simplify = FALSE,
    USE.NAMES = TRUE
  )
  if (unlist) {
    cells <- unname(obj = unlist(x = cells))
  }
  return(cells)
}

#' Pull spatial image names
#'
#' List the names of \code{SpatialImage} objects present in a \code{Seurat} object.
#' If \code{assay} is provided, limits search to images associated with that assay
#'
#' @param object A \code{Seurat} object
#' @param assay Name of assay to limit search to
#'
#' @return A list of image names
#'
#' @export
#'
Images <- function(object, assay = NULL) {
  object <- UpdateSlots(object = object)
  images <- names(x = slot(object = object, name = 'images'))
  if (!is.null(x = assay)) {
    assays <- c(assay, DefaultAssay(object = object[[assay]]))
    images <- Filter(
      f = function(x) {
        return(DefaultAssay(object = object[[x]]) %in% assays)
      },
      x = images
    )
  }
  return(images)
}

#' Load a 10X Genomics Visium Image
#'
#' @param image.dir Path to directory with 10X Genomics visium image data;
#' should include files \code{tissue_lowres_iamge.png},
#' \code{scalefactors_json.json} and \code{tissue_positions_list.csv}
#' @param filter.matrix Filter spot/feature matrix to only include spots that
#' have been determined to be over tissue.
#' @param ... Ignored for now
#'
#' @return A \code{\link{VisiumV1}} object
#'
#' @importFrom png readPNG
#' @importFrom jsonlite fromJSON
#'
#' @seealso \code{\link{VisiumV1}} \code{\link{Load10X_Spatial}}
#'
#' @export
#'
Read10X_Image <- function(image.dir, filter.matrix = TRUE, ...) {
  image <- readPNG(source = file.path(image.dir, 'tissue_lowres_image.png'))
  scale.factors <- fromJSON(txt = file.path(image.dir, 'scalefactors_json.json'))
  tissue.positions <- read.csv(
    file = file.path(image.dir, 'tissue_positions_list.csv'),
    col.names = c('barcodes', 'tissue', 'row', 'col', 'imagerow', 'imagecol'),
    header = FALSE,
    as.is = TRUE,
    row.names = 1
  )
  if (filter.matrix) {
    tissue.positions <- tissue.positions[which(x = tissue.positions$tissue == 1), , drop = FALSE]
  }
  unnormalized.radius <- scale.factors$fiducial_diameter_fullres * scale.factors$tissue_lowres_scalef
  spot.radius <-  unnormalized.radius / max(dim(x = image))
  return(new(
    Class = 'VisiumV1',
    image = image,
    scale.factors = scalefactors(
      spot = scale.factors$tissue_hires_scalef,
      fiducial = scale.factors$fiducial_diameter_fullres,
      hires = scale.factors$tissue_hires_scalef,
      scale.factors$tissue_lowres_scalef
    ),
    coordinates = tissue.positions,
    spot.radius = spot.radius
  ))
}

#' Load Slide-seq spatial data
#'
#' @param coord.file Path to csv file containing bead coordinate positions
#' @param assay Name of assay to associate image to
#'
#' @return A \code{\link{SlideSeq}} object
#'
#' @importFrom utils read.csv
#'
#' @seealso \code{\link{SlideSeq}}
#'
#' @export
#'
ReadSlideSeq <- function(coord.file, assay = 'Spatial') {
  if (!file.exists(paths = coord.file)) {
    stop("Cannot find coord file ", coord.file, call. = FALSE)
  }
  slide.seq <- new(
    Class = 'SlideSeq',
    assay = assay,
    coordinates = read.csv(
      file = coord.file,
      header = TRUE,
      as.is = TRUE,
      row.names = 1
    )
  )
  return(slide.seq)
}

#' Filter stray beads from Slide-seq puck
#' 
#' This function is useful for removing stray beads that fall outside the main 
#' Slide-seq puck area. Essentially, it's a circular filter where you set a 
#' center and radius defining a circle of beads to keep. If the center is not 
#' set, it will be estimated from the bead coordinates (removing the 1st and 
#' 99th quantile to avoid skewing the center by the stray beads). By default,
#' this function will display a \code{\link{SpatialDimPlot}} showing which cells
#' were removed for easy adjustment of the center and/or radius.
#' 
#' @param object Seurat object with slide-seq data
#' @param image Name of the image where the coordinates are stored
#' @param center Vector specifying the x and y coordinates for the center of the 
#' inclusion circle
#' @param radius Radius of the circle of inclusion
#' @param do.plot Display a \code{\link{SpatialDimPlot}} with the cells being 
#' removed labeled.
#' 
#' @return Returns a Seurat object with only the subset of cells that pass the 
#' circular filter
#' 
#' @examples 
#' \dontrun{
#' # This example uses the ssHippo dataset which you can download 
#' # using the SeuratData package.
#' library(SeuratData)
#' data('ssHippo')
#' # perform filtering of beads
#' ssHippo.filtered <- FilterSlideSeq(ssHippo, radius = 2300)
#' # This radius looks to small so increase and repeat until satisfied
#' }
#' @export
#' 
FilterSlideSeq <- function(
  object, 
  image = "image", 
  center = NULL, 
  radius = NULL, 
  do.plot = TRUE
) {
  if (!inherits(x = object[[image]], what = "SlideSeq")) {
    warning("This fxn is intended for filtering SlideSeq data and is untested ",
            "outside of that context.")
  }
  dat <- GetTissueCoordinates(object[[image]])
  if (is.null(x = center)) {
    # heuristic for determining center of puck
    center <- c()
    x.vals <- dat[, 1]
    center[1] <- mean(x.vals[x.vals < quantile(x = x.vals, probs = 0.99) &
                            x.vals > quantile(x = x.vals, probs = 0.01)])
    y.vals <- dat[, 2]
    center[2] <- mean(y.vals[y.vals < quantile(x = y.vals, probs = 0.99) &
                             y.vals > quantile(x = y.vals, probs = 0.01)])
  }
  if (is.null(x = radius)) {
    stop("Please provide a radius.")
  }
  dists <- apply(X = dat, MARGIN = 1, FUN = function(x) {
    as.numeric(dist(rbind(x[c(1, 2)], center)))
  })
  cells.to.remove <- names(x = which(x = (dists > radius)))
  if (do.plot){
    Idents(object) <- "keep"
    object <- SetIdent(object = object, cells = cells.to.remove, value = "remove")
    print(SpatialDimPlot(object = object))
  }
  return(subset(x = object, cells = cells.to.remove, invert = TRUE))
}

#' Load STARmap data
#'
#' @param data.dir location of data directory that contains the counts matrix,
#' gene name, qhull, and centroid files.
#' @param counts.file name of file containing the counts matrix (csv)
#' @param gene.file name of file containing the gene names (csv)
#' @param qhull.file name of file containing the hull coordinates (tsv)
#' @param centroid.file name of file containing the centroid positions (tsv)
#' @param assay Name of assay to associate spatial data to
#' @param image Name of "image" object storing spatial coordinates
#'
#' @return A \code{\link{Seurat}} object
#'
#' @importFrom utils read.csv read.table
#'
#' @seealso \code{\link{STARmap}}
#'
#' @export
#'
LoadSTARmap <- function(
  data.dir,
  counts.file = "cell_barcode_count.csv",
  gene.file = "genes.csv",
  qhull.file = "qhulls.tsv",
  centroid.file = "centroids.tsv",
  assay = "Spatial",
  image = "image"
) {
  if (!dir.exists(paths = data.dir)) {
    stop("Cannot find directory ", data.dir, call. = FALSE)
  }
  counts <- read.csv(
    file = file.path(data.dir, counts.file),
    as.is = TRUE,
    header = FALSE
  )
  gene.names <- read.csv(
    file = file.path(data.dir, gene.file),
    as.is = TRUE,
    header = FALSE
  )
  qhulls <- read.table(
    file = file.path(data.dir, qhull.file),
    sep = '\t',
    col.names = c('cell', 'y', 'x'),
    as.is = TRUE
  )
  centroids <- read.table(
    file = file.path(data.dir, centroid.file),
    sep = '\t',
    as.is = TRUE,
    col.names = c('y', 'x')
  )
  colnames(x = counts) <- gene.names[, 1]
  rownames(x = counts) <- paste0('starmap', seq(1:nrow(x = counts)))
  counts <- as.matrix(x = counts)
  rownames(x = centroids) <- rownames(x = counts)
  qhulls$cell <- paste0('starmap', qhulls$cell)
  centroids <- as.matrix(x = centroids)
  starmap <- CreateSeuratObject(counts = t(x = counts), assay = assay)
  starmap[[image]] <- new(
    Class = 'STARmap',
    assay = assay,
    coordinates = as.data.frame(x = centroids),
    qhulls = qhulls
  )
  return(starmap)
}

#' Visualize spatial and clustering (dimensional reduction) data in a linked,
#' interactive framework
#'
#' @inheritParams DimPlot
#' @inheritParams FeaturePlot
#' @inheritParams SpatialPlot
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
        ) + scale_alpha_ordinal(range = alpha) + guides(alpha = FALSE)
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
#'
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
          guides(alpha = FALSE)
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

#' Load a 10x Genomics Visium Spatial Experiment into a \code{Seurat} object
#'
#' @inheritParams Read10X
#' @inheritParams CreateSeuratObject
#' @param filename Name of H5 file containing the feature barcode matrix
#' @param slice Name for the stored image of the tissue slice
#' @param filter.matrix Only keep spots that have been determined to be over
#' tissue
#' @param to.upper Converts all feature names to upper case. Can be useful when
#' analyses require comparisons between human and mouse gene names for example.
#' @param ... Arguments passed to \code{\link{Read10X_h5}}
#'
#' @return A \code{Seurat} object
#'
#' @importFrom png readPNG
#' @importFrom grid rasterGrob
#' @importFrom jsonlite fromJSON
#'
#' @export
#'
#' @examples
#' \dontrun{
#' data_dir <- 'path/to/data/directory'
#' list.files(data_dir) # Should show filtered_feature_bc_matrix.h5
#' Load10X_Spatial(data.dir = data_dir)
#' }
#'
Load10X_Spatial <- function(
  data.dir,
  filename = 'filtered_feature_bc_matrix.h5',
  assay = 'Spatial',
  slice = 'slice1',
  filter.matrix = TRUE,
  to.upper = FALSE,
  ...
) {
  data <- Read10X_h5(filename = file.path(data.dir, filename), ...)
  if (to.upper) {
    rownames(x = data) <- toupper(x = rownames(x = data))
  }
  object <- CreateSeuratObject(counts = data, assay = assay)
  image <- Read10X_Image(
    image.dir = file.path(data.dir, 'spatial'),
    filter.matrix = filter.matrix
  )
  image <- image[Cells(x = object)]
  DefaultAssay(object = image) <- assay
  object[[slice]] <- image
  return(object)
}

#' Visualize clusters spatially and interactively
#'
#' @inheritParams DimPlot
#' @inheritParams SpatialPlot
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
#' @inheritParams FeaturePlot
#' @inheritParams SpatialPlot
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
#'
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
        guides(alpha = FALSE)
      plot.env$plot
    })
  }
  runGadget(app = ui, server = server)
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

# For plotting the tissue image
#' @importFrom ggplot2 ggproto Geom aes ggproto_parent alpha draw_key_point
#' @importFrom grid unit gpar editGrob pointsGrob viewport gTree addGrob grobName
#' @export
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
      panel_scales$x$continuous_range <- c(0, nrow(x = image))
      panel_scales$y$continuous_range <- c(0, ncol(x = image))
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
#' @export
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

#' Run the mark variogram computation on a given position matrix and expression
#' matrix.
#'
#' Wraps the functionality of markvario from the spatstat package.
#'
#' @param spatial.location A 2 column matrix giving the spatial locations of
#' each of the data points also in data
#' @param data Matrix containing the data used as "marks" (e.g. gene expression)
#' @param ... Arguments passed to markvario
#'
#' @importFrom spatstat markvario ppp
#'
#' @export
#'
RunMarkVario <- function(
  spatial.location,
  data,
  ...
) {
  pp <- ppp(
    x = spatial.location[, 1],
    y = spatial.location[, 2],
    xrange = range(spatial.location[, 1]),
    yrange = range(spatial.location[, 2])
  )
  if (nbrOfWorkers() > 1) {
    chunks <- nbrOfWorkers()
    features <- rownames(x = data)
    features <- split(
      x = features,
      f = ceiling(x = seq_along(along.with = features) / (length(x = features) / chunks))
    )
    mv <- future_lapply(X = features, FUN = function(x) {
      pp[["marks"]] <- as.data.frame(x = t(x = data[x, ]))
      markvario(X = pp, normalise = TRUE, ...)
    })
    mv <- unlist(x = mv, recursive = FALSE)
    names(x = mv) <- rownames(x = data)
  } else {
    pp[["marks"]] <- as.data.frame(x = t(x = data))
    mv <- markvario(X = pp, normalise = TRUE, ...)
  }
  return(mv)
}

#' Compute Moran's I value.
#'
#' Wraps the functionality of the Moran.I function from the ape package.
#' Weights are computed as 1/distance.
#'
#' @param data Expression matrix
#' @param pos Position matrix
#' @param verbose Display messages/progress
#'
#' @importFrom stats dist
#' @importFrom ape Moran.I
#'
#' @export
#'
RunMoransI <- function(data, pos, verbose = TRUE) {
  mysapply <- sapply
  if (verbose) {
    message("Computing Moran's I")
    mysapply <- pbsapply
  }
  Rfast2.installed <- PackageCheck("Rfast2", error = FALSE)
  if (Rfast2.installed) {
    MyMoran <- Rfast2::moranI
  } else {
    MyMoran <- ape::Moran.I
    if (getOption('Seurat.Rfast2.msg', TRUE)) {
      message(
        "For a more efficient implementation of the Morans I calculation,", 
        "\n(selection.method = 'moransi') please install the Rfast2 package",
        "\n--------------------------------------------",
        "\ninstall.packages('Rfast2')",
        "\n--------------------------------------------",
        "\nAfter installation of Rfast2, Seurat will automatically use the more ",
        "\nefficient implementation (no further action necessary).",
        "\nThis message will be shown once per session"
      )
      options(Seurat.Rfast2.msg = FALSE)
    }
  }
  pos.dist <- dist(x = pos)
  pos.dist.mat <- as.matrix(x = pos.dist)
  # weights as 1/dist^2
  weights <- 1/pos.dist.mat^2
  diag(x = weights) <- 0
  results <- mysapply(X = 1:nrow(x = data), FUN = function(x) {
    tryCatch(
      expr = MyMoran(data[x, ], weights),
      error = function(x) c(1,1,1,1)
    )
  })
  results <- data.frame(
    observed = unlist(x = results[1, ]),
    p.value = unlist(x = results[2, ])
  )
  rownames(x = results) <- rownames(x = data)
  return(results)
}

# Base plotting function for all Spatial plots
#
# @param data Data.frame with info to be plotted
# @param image SpatialImage object to be plotted
# @param cols Vector of colors, each color corresponds to an identity class. This may also be a single character
# or numeric value corresponding to a palette as specified by \code{\link[RColorBrewer]{brewer.pal.info}}.
# By default, ggplot2 assigns colors
# @param image.alpha Adjust the opacity of the background images. Set to 0 to
# remove.
# @param crop Crop the plot in to focus on points plotted. Set to FALSE to show
# entire background image.
# @param pt.size.factor Sets the size of the points relative to spot.radius
# @param stroke Control the width of the border around the spots
# @param col.by Mapping variable for the point color
# @param alpha.by Mapping variable for the point alpha value
# @param cells.highlight A list of character or numeric vectors of cells to
# highlight. If only one group of cells desired, can simply pass a vector
# instead of a list. If set, colors selected cells to the color(s) in
# cols.highlight
# @param cols.highlight A vector of colors to highlight the cells as; ordered
# the same as the groups in cells.highlight; last color corresponds to
# unselected cells.
# @param geom Switch between normal spatial geom and geom to enable hover
# functionality
# @param na.value Color for spots with NA values

#' @importFrom tibble tibble
#' @importFrom ggplot2 ggplot aes_string coord_fixed geom_point xlim ylim
#' coord_cartesian labs theme_void theme scale_fill_brewer
#'
SingleSpatialPlot <- function(
  data,
  image,
  cols = NULL,
  image.alpha = 1,
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
      plot + geom_spatial(
        point.size.factor = pt.size.factor,
        data = data,
        image = image,
        image.alpha = image.alpha,
        crop = crop,
        stroke = stroke
      ) + coord_fixed()
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
      scale <- scale_fill_manual(values = cols, na.value = na.value)
    }
    plot <- plot + scale
  }
  plot <- plot + theme_void()
  return(plot)
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
#' @param crop Crop the plot in to focus on points plotted. Set to FALSE to show
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
#' min and max
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
  if (is.null(x = features)) {
    if (interactive) {
      return(ISpatialDimPlot(
        object = object,
        image = image,
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
  if (length(x = images) == 0) {
    images <- Images(object = object)
  }
  if (length(x = images) < 1) {
    stop("Could not find any spatial image information")
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
          guides(alpha = FALSE)
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
  if (length(x = images) > 1 && combine) {
    plots <- wrap_plots(plots = plots, ncol = length(x = images))
  } else if (length(x = images == 1) && combine) {
    plots <- wrap_plots(plots = plots, ncol = ncol)
  }
  return(plots)
}

# TODO: move to convienence.R
#' @export
#' @rdname SpatialPlot
#'
SpatialDimPlot <- function(
  object,
  group.by = NULL,
  images = NULL,
  cols = NULL,
  crop = TRUE,
  cells.highlight = NULL,
  cols.highlight = c('#DE2D26', 'grey50'),
  facet.highlight = FALSE,
  label = FALSE,
  label.size = 7,
  label.color = 'white',
  repel = FALSE,
  ncol = NULL,
  combine = TRUE,
  pt.size.factor = 1.6,
  alpha = c(1, 1),
  stroke = 0.25,
  label.box = TRUE,
  interactive = FALSE,
  information = NULL
) {
  return(SpatialPlot(
    object = object,
    group.by = group.by,
    images = images,
    cols = cols,
    crop = crop,
    cells.highlight = cells.highlight,
    cols.highlight = cols.highlight,
    facet.highlight = facet.highlight,
    label = label,
    label.size = label.size,
    label.color = label.color,
    repel = repel,
    ncol = ncol,
    combine = combine,
    pt.size.factor = pt.size.factor,
    alpha = alpha,
    stroke = stroke,
    label.box = label.box,
    interactive = interactive,
    information = information
  ))
}

# TODO: move to convienence.R
#' @export
#' @rdname SpatialPlot
#'
SpatialFeaturePlot <- function(
  object,
  features,
  images = NULL,
  crop = TRUE,
  slot = 'data',
  min.cutoff = NA,
  max.cutoff = NA,
  ncol = NULL,
  combine = TRUE,
  pt.size.factor = 1.6,
  alpha = c(1, 1),
  stroke = 0.25,
  interactive = FALSE,
  information = NULL
) {
  return(SpatialPlot(
    object = object,
    features = features,
    images = images,
    crop = crop,
    slot = slot,
    min.cutoff = min.cutoff,
    max.cutoff = max.cutoff,
    ncol = ncol,
    combine = combine,
    pt.size.factor = pt.size.factor,
    alpha = alpha,
    stroke = stroke,
    interactive = interactive,
    information = information
  ))
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Methods for Seurat-defined generics
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @method Cells SlideSeq
#' @export
#'
Cells.SlideSeq <- function(x) {
  return(rownames(x = GetTissueCoordinates(object = x)))
}

#' @param x,object An object inheriting from \code{SpatialImage}
#'
#' @rdname SpatialImage-class
#' @name SpatialImage-class
#'
#' @method Cells SpatialImage
#' @export
#'
Cells.SpatialImage <- function(x) {
  stop(
    "'Cells' must be implemented for all subclasses of 'SpatialImage'",
    call. = FALSE
  )
}

#' @method Cells STARmap
#' @export
#'
Cells.STARmap <- function(x) {
  return(rownames(x = GetTissueCoordinates(object = x)))
}

#' @rdname Cells
#' @method Cells VisiumV1
#' @export
#'
Cells.VisiumV1 <- function(x) {
  return(rownames(x = GetTissueCoordinates(object = x, scale = NULL)))
}

#' @rdname DefaultAssay
#' @method DefaultAssay SpatialImage
#' @export
#'
DefaultAssay.SpatialImage <- function(object, ...) {
  CheckDots(...)
  return(slot(object = object, name = 'assay'))
}

#' @rdname DefaultAssay
#' @method DefaultAssay<- SpatialImage
#' @export
#'
"DefaultAssay<-.SpatialImage" <- function(object, ..., value) {
  CheckDots(...)
  slot(object = object, name = 'assay') <- value
  return(object)
}


#' @param object A Seurat object, assay, or expression matrix
#' @param spatial.location Coordinates for each cell/spot/bead
#' @param selection.method Method for selecting spatially variable features.
#'  \itemize{
#'   \item \code{markvariogram}: See \code{\link{RunMarkVario}} for details
#'   \item \code{moransi}: See \code{\link{RunMoransI}} for details.
#' }
#'
#' @param r.metric r value at which to report the "trans" value of the mark
#' variogram
#' @param x.cuts Number of divisions to make in the x direction, helps define 
#' the grid over which binning is performed 
#' @param y.cuts Number of divisions to make in the y direction, helps define 
#' the grid over which binning is performed 
#' @param verbose Print messages and progress
#' 
#' @method FindSpatiallyVariableFeatures default
#' @rdname FindSpatiallyVariableFeatures
#' @export
#'
#'
FindSpatiallyVariableFeatures.default <- function(
  object,
  spatial.location,
  selection.method = c('markvariogram', 'moransi'),
  r.metric = 5,
  x.cuts = NULL,
  y.cuts = NULL,
  verbose = TRUE,
  ...
) {
  # error check dimensions
  if (ncol(x = object) != nrow(x = spatial.location)) {
    stop("Please provide the same number of observations as spatial locations.")
  }
  if (!is.null(x = x.cuts) & !is.null(x = y.cuts)) {
    binned.data <- BinData(
      data = object, 
      pos = spatial.location, 
      x.cuts = x.cuts, 
      y.cuts = y.cuts,
      verbose = verbose
    )
    object <- binned.data$data
    spatial.location <- binned.data$pos
  }
  svf.info <- switch(
    EXPR = selection.method,
    'markvariogram' = RunMarkVario(
      spatial.location = spatial.location,
      data = object
    ),
    'moransi' = RunMoransI(
      data = object,
      pos = spatial.location,
      verbose = verbose
    ),
    stop("Invalid selection method. Please choose one of: markvariogram, moransi.")
  )
  return(svf.info)
}

#' @param slot Slot in the Assay to pull data from
#' @param features If provided, only compute on given features. Otherwise,
#' compute for all features.
#' @param nfeatures Number of features to mark as the top spatially variable.
#' 
#' @method FindSpatiallyVariableFeatures Assay
#' @rdname FindSpatiallyVariableFeatures
#' @export
#'
FindSpatiallyVariableFeatures.Assay <- function(
  object,
  slot = "scale.data",
  spatial.location,
  selection.method = c('markvariogram', 'moransi'),
  features = NULL,
  r.metric = 5,
  x.cuts = NULL,
  y.cuts = NULL,
  nfeatures = nfeatures,
  verbose = TRUE,
  ...
) {
  features <- features %||% rownames(x = object)
  if (selection.method == "markvariogram" && "markvariogram" %in% names(x = Misc(object = object))) {
    features.computed <- names(x = Misc(object = object, slot = "markvariogram"))
    features <- features[! features %in% features.computed]
  }
  data <- GetAssayData(object = object, slot = slot)
  data <- as.matrix(x = data[features, ])
  data <- data[RowVar(x = data) > 0, ]
  if (nrow(x = data) != 0) {
    svf.info <- FindSpatiallyVariableFeatures(
      object = data,
      spatial.location = spatial.location,
      selection.method = selection.method,
      r.metric = r.metric,
      x.cuts = x.cuts,
      y.cuts = y.cuts,
      verbose = verbose,
      ...
    )
  } else {
    svf.info <- c()
  }
  if (selection.method == "markvariogram") {
    if ("markvariogram" %in% names(x = Misc(object = object))) {
      svf.info <- c(svf.info, Misc(object = object, slot = "markvariogram"))
    }
    suppressWarnings(expr = Misc(object = object, slot = "markvariogram") <- svf.info)
    svf.info <- ComputeRMetric(mv = svf.info, r.metric)
    svf.info <- svf.info[order(svf.info[, 1]), , drop = FALSE]
  }
  if (selection.method == "moransi") {
    colnames(x = svf.info) <- paste0("MoransI_", colnames(x = svf.info))
    svf.info <- svf.info[order(svf.info[, 2], -abs(svf.info[, 1])), , drop = FALSE]
  }
  var.name <- paste0(selection.method, ".spatially.variable")
  var.name.rank <- paste0(var.name, ".rank")
  svf.info[[var.name]] <- FALSE
  svf.info[[var.name]][1:(min(nrow(x = svf.info), nfeatures))] <- TRUE
  svf.info[[var.name.rank]] <- 1:nrow(x = svf.info)
  object[[names(x = svf.info)]] <- svf.info
  return(object)
}

#' @param assay Assay to pull the features (marks) from
#' @param image Name of image to pull the coordinates from
#'
#' @method FindSpatiallyVariableFeatures Seurat
#' @rdname FindSpatiallyVariableFeatures
#' @export
#'
FindSpatiallyVariableFeatures.Seurat <- function(
  object,
  assay = NULL,
  slot = "scale.data",
  features = NULL,
  image = NULL,
  selection.method = c('markvariogram', 'moransi'),
  r.metric = 5,
  x.cuts = NULL,
  y.cuts = NULL,
  nfeatures = 2000,
  verbose = TRUE,
  ...
) {
  assay <- assay %||% DefaultAssay(object = object)
  features <- features %||% rownames(x = object[[assay]])
  image <- image %||% DefaultImage(object = object)
  tc <- GetTissueCoordinates(object = object[[image]])
  # check if markvariogram has been run on necessary features
  # only run for new ones
  object[[assay]] <- FindSpatiallyVariableFeatures(
    object = object[[assay]],
    slot = slot,
    features = features,
    spatial.location = tc,
    selection.method = selection.method,
    r.metric = r.metric,
    x.cuts = x.cuts,
    y.cuts = y.cuts,
    nfeatures = nfeatures,
    verbose = verbose,
    ...
  )
  object <- LogSeuratCommand(object = object)
}

#' @param image Name of \code{SpatialImage} object to pull image data for; if
#' \code{NULL}, will attempt to select an image automatically
#'
#' @rdname GetImage
#' @method GetImage Seurat
#' @export
#'
GetImage.Seurat <- function(
  object,
  mode = c('grob', 'raster', 'plotly', 'raw'),
  image = NULL,
  ...
) {
  mode <- match.arg(arg = mode)
  image <- image %||% DefaultImage(object = object)
  if (is.null(x = image)) {
    stop("No images present in this Seurat object", call. = FALSE)
  }
  return(GetImage(object = object[[image]], mode = mode, ...))
}

#' @method GetImage SlideSeq
#' @export
#'
GetImage.SlideSeq <- function(
  object,
  mode = c('grob', 'raster', 'plotly', 'raw'),
  ...
) {
  mode <- match.arg(arg = mode)
  return(NullImage(mode = mode))
}

#' @inheritParams GetImage
#'
#' @rdname SpatialImage-class
#' @name SpatialImage-class
#'
#' @method GetImage SpatialImage
#' @export
#'
GetImage.SpatialImage <- function(
  object,
  mode = c('grob', 'raster', 'plotly', 'raw'),
  ...
) {
  mode <- match.arg(arg = mode)
  stop(
    "'GetImage' must be implemented for all subclasses of 'SpatialImage'",
    call. = FALSE
  )
}

#' @method GetImage STARmap
#' @export
#'
GetImage.STARmap <- function(
  object,
  mode = c('grob', 'raster', 'plotly', 'raw'),
  ...
) {
  mode <- match.arg(arg = mode)
  return(NullImage(mode = mode))
}

#' @importFrom plotly raster2uri
#' @importFrom grDevices as.raster
#' @importFrom grid rasterGrob unit
#'
#' @rdname GetImage
#' @method GetImage VisiumV1
#' @export
#'
GetImage.VisiumV1 <- function(
  object,
  mode = c('grob', 'raster', 'plotly', 'raw'),
  ...
) {
  mode <- match.arg(arg = mode)
  image <- slot(object = object, name = 'image')
  image <- switch(
    EXPR = mode,
    'grob' = rasterGrob(
      image = image,
      width = unit(x = 1, units = 'npc'),
      height = unit(x = 1, units = 'npc')
    ),
    'raster' = as.raster(x = image),
    'plotly' = list(
      source = raster2uri(r = GetImage(object = object, mode = 'raster')),
      xref = 'x',
      yref = 'y',
      # x = -7,
      # y = -7,
      sizex = ncol(x = object),
      sizey = nrow(x = object),
      sizing = 'stretch',
      opacity = 1,
      layer = 'below'
    ),
    'raw' = image,
    stop("Unknown image mode: ", mode, call. = FALSE)
  )
  return(image)
}

#' @param image Name of \code{SpatialImage} object to get coordinates for; if
#' \code{NULL}, will attempt to select an image automatically
#'
#' @rdname GetTissueCoordinates
#' @method GetTissueCoordinates Seurat
#' @export
#'
GetTissueCoordinates.Seurat <- function(object, image = NULL, ...) {
  image <- image %||% DefaultImage(object = object)
  if (is.null(x = image)) {
    stop("No images present in this Seurat object", call. = FALSE)
  }
  return(GetTissueCoordinates(object = object[[image]], ...))
}

#' @method GetTissueCoordinates SlideSeq
#' @export
#'
GetTissueCoordinates.SlideSeq <- function(object, ...) {
  coords <- slot(object = object, name = 'coordinates')
  colnames(x = coords) <- c('x', 'y')
  # coords$y <- -rev(x = coords$y) + 1
  # coords$y <- FlipCoords(x = coords$y)
  coords$cells <- rownames(x = coords)
  return(coords)
}

#' @inheritParams GetTissueCoordinates
#'
#' @rdname SpatialImage-class
#' @name SpatialImage-class
#'
#' @method GetTissueCoordinates SpatialImage
#' @export
#'
GetTissueCoordinates.SpatialImage <- function(object, ...) {
  stop(
    "'GetTissueCoordinates' must be implemented for all subclasses of 'SpatialImage'",
    call. = FALSE
  )
}

#' @param qhulls return qhulls instead of centroids
#' @method GetTissueCoordinates STARmap
#' @export
#'
GetTissueCoordinates.STARmap <- function(object, qhulls = FALSE, ...) {
  if (qhulls) {
    return(slot(object = object, name = 'qhulls'))
  }
  return(slot(object = object, name = 'coordinates'))
}

#' @param scale A factor to scale the coordinates by; choose from: 'tissue',
#' 'fiducial', 'hires', 'lowres', or \code{NULL} for no scaling
#' @param cols Columns of tissue coordinates data.frame to pull
#'
#' @rdname GetTissueCoordinates
#' @method GetTissueCoordinates VisiumV1
#' @export
#'
GetTissueCoordinates.VisiumV1 <- function(
  object,
  scale = 'lowres',
  cols = c('imagerow', 'imagecol'),
  ...
) {
  cols <- cols %||% colnames(x = slot(object = object, name = 'coordinates'))
  if (!is.null(x = scale)) {
    coordinates <- slot(object = object, name = 'coordinates')[, c('imagerow', 'imagecol')]
    scale <- match.arg(arg = scale, choices = c('spot', 'fiducial', 'hires', 'lowres'))
    scale.use <- ScaleFactors(object = object)[[scale]]
    coordinates <- coordinates * scale.use
  } else {
    coordinates <- slot(object = object, name = 'coordinates')[, cols]
  }
  return(coordinates)
}

#' @rdname IsGlobal
#' @method IsGlobal SpatialImage
#' @export
#'
IsGlobal.SpatialImage <- function(object) {
  return(TRUE)
}

#' @rdname Key
#' @method Key SpatialImage
#' @export
#'
Key.SpatialImage <- function(object, ...) {
  CheckDots(...)
  object <- UpdateSlots(object = object)
  return(slot(object = object, name = 'key'))
}

#' @rdname Key
#' @method Key<- SpatialImage
#' @export
#'
"Key<-.SpatialImage" <- function(object, ..., value) {
  CheckDots(...)
  object <- UpdateSlots(object = object)
  value <- UpdateKey(key = value)
  slot(object = object, name = 'key') <- value
  return(object)
}

#' @rdname Radius
#' @method Radius SlideSeq
#' @export
#'
Radius.SlideSeq <- function(object) {
  return(0.005)
}

#' @rdname SpatialImage-class
#' @name SpatialImage-class
#'
#' @method Radius SpatialImage
#' @export
#'
Radius.SpatialImage <- function(object) {
  return(NULL)
}

#' @rdname Radius
#' @method Radius STARmap
#' @export
#'
Radius.STARmap <- function(object) {
  return(NULL)
}

#' @rdname Radius
#' @method Radius VisiumV1
#' @export
#'
Radius.VisiumV1 <- function(object) {
  return(slot(object = object, name = 'spot.radius'))
}

#' @method RenameCells SlideSeq
#' @export
#'
RenameCells.SlideSeq <- function(object, new.names = NULL, ...) {
  return(RenameCells.VisiumV1(object = object, new.names = new.names))
}

#' @inheritParams  RenameCells
#'
#' @rdname SpatialImage-class
#' @name SpatialImage-class
#'
#' @method RenameCells SpatialImage
#' @export
#'
RenameCells.SpatialImage <- function(object, new.names = NULL, ...) {
  stop(
    "'RenameCells' must be implemented for all subclasses of 'SpatialImage'",
    call. = FALSE
  )
}

#' @method RenameCells STARmap
#' @export
#'
RenameCells.STARmap <- function(object, new.names = NULL, ...) {
  names(x = new.names) <- Cells(x = object)
  object <- RenameCells.VisiumV1(object = object, new.names = new.names)
  qhulls <- GetTissueCoordinates(object = object, qhull = TRUE)
  qhulls$cell <- new.names[qhulls$cell]
  slot(object = object, name = "qhulls") <- qhulls
  return(object)
}

#' @rdname RenameCells
#' @method RenameCells VisiumV1
#' @export
#'
RenameCells.VisiumV1 <- function(object, new.names = NULL, ...) {
  if (is.null(x = new.names)) {
    return(object)
  } else if (length(x = new.names) != length(x = Cells(x = object))) {
    stop("Wrong number of cell/spot names", call. = FALSE)
  }
  names(x = new.names) <- Cells(x = object)
  coordinates <- GetTissueCoordinates(object = object, scale = NULL, cols = NULL)
  rownames(x = coordinates) <- new.names[rownames(x = coordinates)]
  slot(object = object, name = 'coordinates') <- coordinates
  return(object)
}

#' @rdname ScaleFactors
#' @method ScaleFactors VisiumV1
#' @export
#'
ScaleFactors.VisiumV1 <- function(object, ...) {
  return(slot(object = object, name = 'scale.factors'))
}

#' @inheritParams FindSpatiallyVariableFeatures
#' @param decreasing Return features in decreasing order (most spatially 
#' variable first).
#' @rdname SpatiallyVariableFeatures
#' @export
#' @method SpatiallyVariableFeatures Assay
#'
SpatiallyVariableFeatures.Assay <- function(
  object,
  selection.method = "markvariogram",
  decreasing = TRUE,
  ...
) {
  CheckDots(...)
  vf <- SVFInfo(object = object, selection.method = selection.method, status = TRUE)
  vf <- vf[rownames(x = vf)[which(x = vf[, "variable"][, 1])], ]
  if (!is.null(x = decreasing)) {
    vf <- vf[order(x = vf[, "rank"], decreasing = !decreasing), ]
  }
  return(rownames(x = vf)[which(x = vf[, "variable"][, 1])])
}

#' @param Seurat object
#' @param assay Name of assay to pull spatially variable features for
#'
#' @rdname SpatiallyVariableFeatures
#' @export
#' @method SpatiallyVariableFeatures Seurat
#'
SpatiallyVariableFeatures.Seurat <- function(
  object,
  assay = NULL,
  selection.method = "markvariogram",
  decreasing = TRUE,
  ...
) {
  CheckDots(...)
  assay <- assay %||% DefaultAssay(object = object)
  return(SpatiallyVariableFeatures(object = object[[assay]], selection.method = selection.method, decreasing = decreasing))
}

#' @param selection.method Which method to pull. Options: markvariogram, moransi
#' @param status Add variable status to the resulting data.frame
#'
#' @rdname SVFInfo
#' @export
#' @method SVFInfo Assay
#'
SVFInfo.Assay <- function(object, selection.method = c("markvariogram", "moransi"), status = FALSE, ...) {
  CheckDots(...)
  vars <- switch(
    EXPR = selection.method,
    'markvariogram' = grep(pattern = "r.metric", x = colnames(x = object[[]]), value = TRUE),
    'moransi' = grep(pattern = 'moransi', x = colnames(x = object[[]]), value = TRUE),
    stop("Unknown method: '", selection.method, "'", call. = FALSE)
  )
  tryCatch(
    expr = svf.info <- object[[vars]],
    error = function(e) {
      stop(
        "Unable to find highly variable feature information for method '",
        selection.method,
        "'",
        call. = FALSE
      )
    }
  )
  colnames(x = svf.info) <- vars
  if (status) {
    svf.info$variable <- object[[paste0(selection.method, '.spatially.variable')]]
    svf.info$rank <- object[[paste0(selection.method, '.spatially.variable.rank')]]
  }
  return(svf.info)
}

#' @param assay Name of assay to pull highly variable feature information for
#'
#' @importFrom tools file_path_sans_ext
#'
#' @rdname SVFInfo
#' @export
#' @method SVFInfo Seurat
#'
SVFInfo.Seurat <- function(
  object,
  selection.method = c("markvariogram", "moransi"),
  assay = NULL,
  status = FALSE,
  ...
) {
  CheckDots(...)
  assay <- assay %||% DefaultAssay(object = object)
  return(SVFInfo(
    object = GetAssay(object = object, assay = assay),
    selection.method = selection.method,
    status = status
  ))
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Methods for R-defined generics
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @method [ SlideSeq
#' @export
#'
"[.SlideSeq" <- function(x, i, ...) {
  return(subset(x = x, cells = i, ...))
}

#' @param i,cells A vector of cells to keep
#'
#' @rdname SpatialImage-class
#' @name SpatialImage-class
#'
#' @method [ SpatialImage
#' @export
#'
"[.SpatialImage" <- function(x, i, ...) {
  stop(
    "'[' must be implemented for all subclasses of 'SpatialImage'",
    call. = FALSE
  )
}

#' @method [ VisiumV1
#' @export
#'
"[.VisiumV1" <- function(x, i, ...) {
  return(subset(x = x, cells = i))
}

#' @method dim SlideSeq
#' @export
#'
dim.SlideSeq <- function(x) {
  # return(dim(x = GetImage(object = x, mode = 'raw')))
  return(c(599, 600))
}

#' @rdname SpatialImage-class
#' @name SpatialImage-class
#'
#' @method dim SpatialImage
#' @export
#'
dim.SpatialImage <- function(x) {
  stop(
    "'dim' must be implemented for all subclasses of 'SpatialImage'",
    call. = FALSE
  )
}

#' @method dim STARmap
#' @export
#'
dim.STARmap <- function(x) {
  coords <- GetTissueCoordinates(object = x)
  return(c(
    max(coords[, 1]) - min(coords[, 1]),
    max(coords[, 2]) - min(coords[, 2])
  ))
}

#' @method dim VisiumV1
#' @export
#'
dim.VisiumV1 <- function(x) {
  return(dim(x = GetImage(object = x)$raster))
}

#' @method subset SlideSeq
#' @export
#'
subset.SlideSeq <- function(x, cells, ...) {
  x <- subset.VisiumV1(x = x, cells = cells, ...)
  return(x)
}

#' @method subset STARmap
#' @export
#'
subset.STARmap <- function(x, cells, ...) {
  x <- subset.VisiumV1(x = x, cells = cells, ...)
  qhulls <- GetTissueCoordinates(object = x, qhulls = TRUE)
  qhulls <- qhulls[qhulls$cell %in% cells, ]
  slot(object = x, name = 'qhulls') <- qhulls
  return(x)
}

#' @rdname SpatialImage-class
#' @name SpatialImage-class
#'
#' @method subset SpatialImage
#' @export
#'
subset.SpatialImage <- function(x, cells, ...) {
  stop("'subset' must be implemented for all subclasses of 'SpatialImage'")
}

#' @method subset VisiumV1
#' @export
#'
subset.VisiumV1 <- function(x, cells, ...) {
  coordinates <- GetTissueCoordinates(object = x, scale = NULL, cols = NULL)
  cells <- cells[cells %in% rownames(x = coordinates)]
  coordinates <- coordinates[cells, ]
  slot(object = x, name = 'coordinates') <- coordinates
  return(x)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# S4 methods
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

setMethod(
  f = 'show',
  signature = 'SpatialImage',
  definition = function(object) {
    object <- UpdateSlots(object = object)
    cat(
      "Spatial data from the",
      class(x = object),
      "technology for",
      length(x = Cells(x = object)),
      "samples\n"
    )
    cat("Associated assay:", DefaultAssay(object = object), "\n")
    cat("Image key:", Key(object = object), "\n")
  }
)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Internal
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Bin spatial regions into grid and average expression values
#
# @param dat Expression data
# @param pos Position information/coordinates for each sample
# @param x.cuts Number of cuts to make in the x direction (defines grid along
# with y.cuts)
# @param y.cuts Number of cuts to make in the y direction
#
# @return returns a list with positions as centers of the bins and average
# expression within the bins
#
#' @importFrom Matrix rowMeans
#
BinData <- function(data, pos, x.cuts = 10, y.cuts = x.cuts, verbose = TRUE) {
  if (verbose) {
    message("Binning spatial data")
  }
  pos$x.cuts <- cut(x = pos[, 1], breaks = x.cuts)
  pos$y.cuts <- cut(x = pos[, 2], breaks = y.cuts)
  pos$bin <- paste0(pos$x.cuts, "_", pos$y.cuts)
  all.bins <- unique(x = pos$bin)
  new.pos <- matrix(data = numeric(), nrow = length(x = all.bins), ncol = 2)
  new.dat <- matrix(data = numeric(), nrow = nrow(x = data), ncol = length(x = all.bins))
  for(i in 1:length(x = all.bins)) {
    samples <- rownames(x = pos)[which(x = pos$bin == all.bins[i])]
    dat <- data[, samples]
    if (is.null(x = dim(x = dat))) {
      new.dat[, i] <- dat
    } else {
      new.dat[, i] <- rowMeans(data[, samples])
    }
    new.pos[i, 1] <- mean(pos[samples, "x"])
    new.pos[i, 2] <- mean(pos[samples, "y"])
  }
  rownames(x = new.dat) <- rownames(x = data)
  colnames(x = new.dat) <- all.bins
  rownames(x = new.pos) <- all.bins
  colnames(x = new.pos) <- colnames(x = pos)[1:2]
  return(list(data = new.dat, pos = new.pos))
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
#'
#' @keywords color
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

# Computes the metric at a given r (radius) value and stores in meta.features
#
# @param mv Results of running markvario
# @param r.metric r value at which to report the "trans" value of the mark
# variogram
#
# @return Returns a data.frame with r.metric values
#
#
ComputeRMetric <- function(mv, r.metric = 5) {
  r.metric.results <- unlist(x = lapply(
    X = mv,
    FUN = function(x) {
      x$trans[which.min(x = abs(x = x$r - r.metric))]
    }
  ))
  r.metric.results <- as.data.frame(x = r.metric.results)
  colnames(r.metric.results) <- paste0("r.metric.", r.metric)
  return(r.metric.results)
}

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

#' Compute the correlation of features broken down by groups with another
#' covariate
#'
#' @param object Seurat object
#' @param assay Assay to pull the data from
#' @param slot Slot in the assay to pull feature expression data from (counts,
#' data, or scale.data)
#' @param var Variable with which to correlate the features
#' @param group.assay Compute the gene groups based off the data in this assay.
#' @param min.cells Only compute for genes in at least this many cells
#' @param ngroups Number of groups to split into
#' @param do.plot Display the group correlation boxplot (via 
#' \code{GroupCorrelationPlot})
#'
#' @return A Seurat object with the correlation stored in metafeatures
#'
#' @export
#'
GroupCorrelation <- function(
  object,
  assay = NULL,
  slot = "scale.data",
  var = NULL,
  group.assay = NULL,
  min.cells = 5,
  ngroups = 6,
  do.plot = TRUE
) {
  assay <- assay %||% DefaultAssay(object = object)
  group.assay <- group.assay %||% assay
  var <- var %||% paste0("nCount_", group.assay)
  gene.grp <- GetFeatureGroups(
    object = object,
    assay = group.assay,
    min.cells = min.cells,
    ngroups = ngroups
  )
  data <- as.matrix(x = GetAssayData(object = object[[assay]], slot = slot))
  data <- data[rowMeans(x = data) != 0, ]
  grp.cors <- apply(
    X = data,
    MARGIN = 1,
    FUN = function(x) {
      cor(x = x, y = object[[var]])
    }
  )
  grp.cors <- grp.cors[names(x = gene.grp)]
  grp.cors <- as.data.frame(x = grp.cors[which(x = !is.na(x = grp.cors))])
  grp.cors$gene_grp <- gene.grp[rownames(x = grp.cors)]
  colnames(x = grp.cors) <- c("cor", "feature_grp")
  object[[assay]][["feature.grp"]] <- grp.cors[, "feature_grp", drop = FALSE]
  object[[assay]][[paste0(var, "_cor")]] <- grp.cors[, "cor", drop = FALSE]
  if (do.plot) {
    print(GroupCorrelationPlot(
      object = object,
      assay = assay,
      feature.group = "feature.grp",
      cor = paste0(var, "_cor")
    ))
  }
  return(object)
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
#'
GroupCorrelationPlot <- function(
  object,
  assay = NULL,
  feature.group = "feature.grp",
  cor = "nCount_RNA_cor"
) {
  assay <- assay %||% DefaultAssay(object = object)
  data <- object[[assay]][[c(feature.group, cor)]]
  data <- data[complete.cases(data), ]
  colnames(x = data) <- c('grp', 'cor')
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

# Get the default image of an object
#
# Attempts to find all images associated with the default assay of the object.
# If none present, finds all images present in the object. Returns the name of
# the first image
#
# @param object A Seurat object
#
# @return The name of the default image
#
DefaultImage <- function(object) {
  object <- UpdateSlots(object = object)
  images <- Images(object = object, assay = DefaultAssay(object = object))
  if (length(x = images) < 1) {
    images <- Images(object = object)
  }
  return(images[[1]])
}

#' @importFrom ggplot2 ggplot_build
#'
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
#'
#' @keywords color
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

#' @name contrast-theory
#' @rdname contrast-theory
#'
#' @importFrom grDevices col2rgb
#'
#' @export
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

# Given a range from cut, compute the mean
#
# @x range from cut as a string (e.g. (10, 20] )
# @return returns a numeric with the mean of the range
#
MeanRange <- function(x) {
  left <- gsub(pattern = "\\]", replacement = "", x = sub(pattern = "\\([[:digit:]\\.e+]*,", x = x, replacement = ""))
  right <- gsub(pattern = "\\(", replacement = "", x = sub(pattern = ",[[:digit:]\\.e+]*]", x = x, replacement = ""))
  return(mean(c(as.numeric(x = left), as.numeric(x = right))))
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

# Return a null image
#
# @param mode Image representation to return
# see \code{\link{GetImage}} for more details
#
#' @importFrom grid nullGrob
#' @importFrom grDevices as.raster
#
NullImage <- function(mode) {
  image <- switch(
    EXPR = mode,
    'grob' = nullGrob(),
    'raster' = as.raster(x = new(Class = 'matrix')),
    'plotly' = list('visible' = FALSE),
    'raw' = NULL,
    stop("Unknown image mode: ", mode, call. = FALSE)
  )
  return(image)
}
