#' Build a landmark object with correlation analysis and UMAP embedding
#'
#' @param object Seurat object
#' @param sample.obj Sample-level Seurat object
#' @param rename.sample Optional renaming for samples
#' @param nn.name Name of the neighbor object to use
#' @param weighted.nn.name Deprecated. Use \code{nn.name} instead.
#' @param landmark.assay.name Name of landmark assay
#' @param features.to.test Features to test for correlation
#' @param k.nn Number of nearest neighbors. Default is 20.
#' @param min.dist Minimum distance for UMAP. Default is 0.3.
#' @param ... Additional arguments
#'
#' @return Landmark object with UMAP and correlation data
#' @export
#' @concept visualization
#'
#' @importFrom SeuratObject CreateSeuratObject CreateAssayObject AddMetaData FetchData LayerData
#' @importFrom methods new
#' @importFrom stats cor.test
#'
BuildLandmarkObject <- function(
    object = NULL,
    sample.obj = NULL,
    rename.sample = NULL,
    nn.name = "weighted.nn",
    weighted.nn.name = NULL,
    landmark.assay.name = "LANDMARK",
    features.to.test = NULL,
    k.nn = 20,
    min.dist = 0.3,
    ...){
  # deprecation shim: 'weighted.nn.name' -> 'nn.name'
  if (!is.null(weighted.nn.name)) {
    warning("'weighted.nn.name' is deprecated; please use 'nn.name' instead.")
    nn.name <- weighted.nn.name
  }
  # first normalize the library size for the sample-level density matrix
  sample.obj <- NormalizeData(sample.obj, ...)
  # and retrieve the normalized matrix
  landmark_mt <- LayerData(sample.obj, layer = "data")
  landmark_mt <- t(landmark_mt)

  # rename the sample if rename.sample is provided
  if(!is.null(x = rename.sample)){
    rownames(landmark_mt) <- rename.sample
  }
  # and create a new cell-level object by transposing the norm.matrix (sample is now treated as feature)
  landmark_obj <- CreateSeuratObject(counts = CreateAssayObject(data = as.matrix(landmark_mt)))

  # next we will retrive the WNN object we previously computed (stored in the single-cell object)
  landmark_nn <- object[[nn.name]]
  nn.idx <- as.matrix(landmark_nn@nn.idx)
  nn.dist <- as.matrix(landmark_nn@nn.dist)
  # keep only the landmark cells in the NN object
  slct_row <- match(colnames(object[[landmark.assay.name]]), landmark_nn@cell.names)
  #
  landmark_nn2 <- new(
    Class = 'Neighbor',
    nn.idx = landmark_nn@nn.idx[slct_row, 1:k.nn],
    nn.dist = landmark_nn@nn.dist[slct_row, 1:k.nn],
    alg.info = landmark_nn@alg.info,
    cell.names = colnames(landmark_obj)
  )
  landmark_obj[[nn.name]] <- landmark_nn2

  # and now we will generate UMAP based on the landmark-only WNN graph
  landmark_obj <- RunUMAP(landmark_obj, nn.name = nn.name, reduction.name = "wnn.umap", reduction.key = "wnnUMAP_",
                          min.dist = min.dist)
  # paste the meta-data of the landmark cells
  landmark_obj <- AddMetaData(landmark_obj, metadata = object@meta.data[colnames(object[[landmark.assay.name]]), ])

  # now we will perform cor.test between a given "trajectory" column(s) and all the landmarks
  # we aim to identify the top-associated landmarks for the given trajectory
  exp_data <- LayerData(sample.obj, layer = "data")
  feature_data <- FetchData(sample.obj, vars = features.to.test)
  #
  merge_dat <- matrix(nrow = nrow(exp_data), ncol = ncol(feature_data))
  j <- 1
  for(TRAJ_NAME in colnames(feature_data)){
    cor_res = vector(length = nrow(exp_data))
    p_res = vector(length = nrow(exp_data))
    #
    for(i in 1:nrow(exp_data)){
      tmp_res <- cor.test(exp_data[i, ], feature_data[[TRAJ_NAME]], na.action = "na.omit")
      cor_res[i] <- ifelse(is.na(tmp_res$estimate), 0, tmp_res$estimate)
      p_res[i] <- ifelse(is.na(tmp_res$p.value), 1, tmp_res$p.value)
    }
    #
    names(cor_res) = rownames(exp_data)
    names(p_res) = rownames(exp_data)
    #
    merge_dat[, j] = cor_res
    j = j + 1
  }
  rownames(merge_dat) <- rownames(exp_data)
  colnames(merge_dat) <- paste0("cor.coef_", colnames(feature_data))
  #
  landmark_obj <- AddMetaData(landmark_obj, metadata = merge_dat)

  return(landmark_obj)
}

#' Visualize landmark-trajectory correlations with UMAP plots
#'
#' @param landmark_obj Landmark object from BuildLandmarkObject
#' @param order Whether to order points by expression
#' @param pt.size Point size for plots
#' @param alpha Point transparency
#' @param plot.featureplot Whether to create feature plots
#' @param ncol Number of columns for faceting
#' 
#' @return List of ggplot objects
#' @export
#' @concept visualization
#'
#' @importFrom dplyr %>% arrange slice_sample
#' @importFrom ggplot2 ggplot aes geom_point facet_wrap scale_color_gradientn theme_minimal theme element_blank element_text labs coord_fixed
#' @importFrom rlang .data
#'
PlotLandmarkObject <- function(
    landmark_obj,
    order = TRUE,
    pt.size = 0.5,
    alpha = 0.8,
    plot.featureplot = FALSE,
    ncol = NULL){

  if (!requireNamespace("tidyr", quietly = TRUE)) {
    stop("The 'tidyr' package is required for PlotLandmarkObject() but is not installed. Please install it with install.packages('tidyr').")
  }

  # Extract correlation coefficient columns
  cor_cols <- grep("^cor.coef_", colnames(landmark_obj@meta.data), value = TRUE)
  merge_dat <- landmark_obj@meta.data[, cor_cols, drop = FALSE]

  #
  umap_coords <- landmark_obj[["wnn.umap"]]@cell.embeddings
  umap_df <- data.frame(
    cell = rownames(umap_coords),
    UMAP_1 = umap_coords[, 1],
    UMAP_2 = umap_coords[, 2]
  )
  #
  plot_data <- cbind(umap_df, merge_dat)

  # Reshape data for faceting (long format)
  plot_data_long <- plot_data %>%
    tidyr::pivot_longer(cols = colnames(merge_dat), names_to = "feature", values_to = "expression")

  # Order cells by expression if requested
  if (order) {
    # Order by absolute value: small to large absolute values (large absolute values on top)
    plot_data_long <- plot_data_long %>%
      arrange(abs(.data$expression))
  } else {
    # Random ordering to avoid any bias
    plot_data_long <- plot_data_long %>%
      slice_sample(n = nrow(plot_data_long))
  }

  if(is.null(x = ncol)) ncol = ncol(merge_dat)

  # Create the plot
  p <- ggplot(plot_data_long, aes(x = .data$UMAP_1, y = .data$UMAP_2, color = .data$expression)) +
    geom_point(size = pt.size, alpha = alpha) +
    facet_wrap(~ .data$feature, ncol = ncol) +

    # Custom diverging color scale (blue-white-red)
    scale_color_gradientn(
      colors = c("#053061", "#4393c3", "#f7f7f7", "#d6604d", "#67001f"),
      name = "Expression"
    ) +

    theme_minimal() +
    theme(
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      aspect.ratio = 1,
      strip.text = element_text(size = 12, face = "bold")
    ) +
    labs(x = "wnnUMAP 1", y = "wnnUMAP 2") +
    coord_fixed()

  # donor-level FeaturePlot
  if(plot.featureplot) {
    p2 = FeaturePlot(landmark_obj, slot = "data", reduction = "wnn.umap", features = rownames(landmark_obj),
                     order = T, ncol = ceiling(sqrt(nrow(landmark_obj))))
  } else {
    p2 = NA
  }

  # return only the plot objects
  return(list(p1 = p, p2 = p2))
}

#' Visualize the landmark-trajectory relevance and the sample-level cell density
#'
#' @param object Seurat object
#' @param sample.obj Sample-level Seurat object
#' @param rename.sample Optional sample renaming
#' @param nn.name Name of the neighbor object to use
#' @param weighted.nn.name Deprecated. Use \code{nn.name} instead.
#' @param landmark.assay.name Name of landmark assay
#' @param features.to.test Features to test for correlation
#' @param order Whether to order points
#' @param pt.size Point size
#' @param alpha Point transparency
#' @param ncol Number of columns
#' @param return.seurat.obj Whether to return Seurat object
#' @param k.nn Number of nearest neighbors
#' @param min.dist Minimum distance for UMAP
#' @param ... Additional arguments
#' 
#' @return List containing plots and optionally landmark object
#' @export
#' @concept visualization
#'
SampleLevelDimPlot <- function(
    object = NULL,
    sample.obj = NULL,
    rename.sample = NULL,
    nn.name = "weighted.nn",
    weighted.nn.name = NULL,
    landmark.assay.name = "LANDMARK",
    features.to.test = NULL,
    order = TRUE,
    pt.size = 0.5,
    alpha = 0.8,
    ncol = NULL,
    return.seurat.obj = TRUE,
    k.nn = 20,
    min.dist = 0.3,
    ...){

  # deprecation shim: 'weighted.nn.name' -> 'nn.name'
  if (!is.null(weighted.nn.name)) {
    warning("'weighted.nn.name' is deprecated; please use 'nn.name' instead.")
    nn.name <- weighted.nn.name
  }

  # Build the landmark object with correlation analysis
  landmark_obj <- BuildLandmarkObject(
    object = object,
    sample.obj = sample.obj,
    rename.sample = rename.sample,
    nn.name = nn.name,
    landmark.assay.name = landmark.assay.name,
    features.to.test = features.to.test,
    k.nn = k.nn,
    min.dist = min.dist,
    ...
  )

  # Create the correlation plots
  plots <- PlotLandmarkObject(
    landmark_obj = landmark_obj,
    order = order,
    pt.size = pt.size,
    alpha = alpha,
    ncol = ncol
  )

  # Return plots and landmark object based on return.seurat.obj parameter
  if(return.seurat.obj) {
    return(list(landmark_obj = landmark_obj, p1 = plots$p1, p2 = plots$p2))
  } else {
    return(plots)
  }
}


#' Generate a UMAP for a sketched assay and then project it to the full data.
#'
#' This is a wrapper function that builds on the ProjectUMAP() function, which allows users to
#' generate UMAP for a subset of the data (e.g., cells from a sketched assay), and then project this
#' UMAP to any other data.
#'
#' @param object a Seurat object
#' @param nn.name the name of the Neighbor object used to generate the UMAP.
#' @param k.nn the number of neighbors used to generate the UMAP. If it exceeds the number of k in the nn.name object, will be ignored.
#' @param landmark.assay.name the name of the sketched assay
#' @param reduction.name the name of the UMAP to be stored; Note that the projected UMAP will be named as "proj.'reduction.name'".
#' @param verbose Print progress status
#' @param ... other parameters that are accepted by RunUMAP() or ProjectUMAP()
#'
#' @return return a Seurat object that contains a UMAP of the sketched assay and a UMAP of the full data
#' @export
#' @concept visualization
#'
#' @importFrom SeuratObject Neighbors
#' @importFrom methods new
#'
RunAndProjectUMAP <- function(
    object,
    nn.name = NULL,
    k.nn = 5,
    landmark.assay.name = "LANDMARK",
    reduction.name = "sketched.umap",
    verbose = TRUE,
    ...
) {
  if(is.null(nn.name) | !nn.name %in% Neighbors(object)) {
    stop("Please specify the correct name of the Neighbor object that stored in the Seurat object.")
  }
  #
  landmark_nn <- object[[nn.name]]
  nn.idx <- as.matrix(landmark_nn@nn.idx)
  nn.dist <- as.matrix(landmark_nn@nn.dist)
  # keep only the landmark cells in the NN object
  slct_row <- match(colnames(object[[landmark.assay.name]]), landmark_nn@cell.names)
  #
  landmark_nn2 <- new(
    Class = 'Neighbor',
    nn.idx = landmark_nn@nn.idx[slct_row, 1:k.nn],
    nn.dist = landmark_nn@nn.dist[slct_row, 1:k.nn],
    alg.info = landmark_nn@alg.info,
    cell.names = colnames(object[[landmark.assay.name]])
  )
  #
  landmark_nn <- new(
    Class = 'Neighbor',
    nn.idx = landmark_nn@nn.idx[, 1:k.nn],
    nn.dist = landmark_nn@nn.dist[, 1:k.nn],
    alg.info = landmark_nn@alg.info,
    cell.names = landmark_nn@cell.names
  )
  object[['sketch_nn']] <- landmark_nn2

  # and now we will generate UMAP based on the landmark-only WNN graph
  object <- RunUMAP(object = object, nn.name = "sketch_nn", reduction.name = reduction.name, return.model = T, ...)

  # project this UMAP to the full data
  full_UMAP <- RunUMAP(object = landmark_nn, reduction.model = object[[reduction.name]], ...)
  object[[paste0("proj.", reduction.name)]] <- full_UMAP

  return(object)
}


# Build a ComplexHeatmap column annotation from a data frame of metadata values.
# Each column becomes one annotation bar; numeric columns get a continuous color
# ramp (min/median/max), categorical columns get one color per level. 'anno.colors'
# is a list keyed by variable name; a NULL entry falls back to a sensible default.
# Assumes ComplexHeatmap and circlize are available (callers guard for this).
.scslide_build_heatmap_annotation <- function(anno.df, anno.colors) {
  col.list <- list()
  for (v in colnames(x = anno.df)) {
    values <- anno.df[[v]]
    cols <- anno.colors[[v]]
    if (is.numeric(x = values)) {
      cols <- cols %||% c("#4575b4", "white", "#d73027")
      if (length(x = cols) != 3) {
        stop("anno.colors for the numeric variable '", v, "' must be three colors ",
             "(for the min, median, and max values).")
      }
      col.list[[v]] <- circlize::colorRamp2(
        breaks = c(min(values, na.rm = TRUE),
                   stats::median(x = values, na.rm = TRUE),
                   max(values, na.rm = TRUE)),
        colors = cols
      )
    } else {
      levs <- levels(x = factor(x = values))
      cols <- cols %||% scales::hue_pal()(length(x = levs))
      if (length(x = cols) < length(x = levs)) {
        stop("anno.colors for the categorical variable '", v, "' must supply at least ",
             length(x = levs), " colors (one per level).")
      }
      cols <- cols[seq_along(along.with = levs)]
      names(x = cols) <- levs
      col.list[[v]] <- cols
    }
  }
  ComplexHeatmap::HeatmapAnnotation(df = anno.df, col = col.list, border = TRUE)
}

#' Plot an annotated expression heatmap of features across cells
#'
#' Draws a features-by-cells expression heatmap with one or more metadata
#' annotation bars on top, built on \code{\link[ComplexHeatmap]{Heatmap}}. Columns
#' are ordered by a chosen metadata variable (e.g. a trajectory value) independently
#' of the annotation bars, which makes it suited to trajectory / differential-expression
#' heatmaps where cells are ordered along a continuous axis while still displaying a
#' categorical annotation (something the grouped \code{\link{DoHeatmap}} cannot do).
#'
#' This is an independent implementation whose interface is inspired by
#' \code{Scillus::plot_heatmap}; it does not reuse that package's code.
#'
#' @param object A Seurat object.
#' @param features Character vector of features (rows of the heatmap). Features
#'   absent from the requested layer are dropped with a warning.
#' @param sort.by One or more metadata column names used to order the cells
#'   (columns of the heatmap). Multiple names are applied in order as tie-breakers.
#' @param anno.var Optional character vector of metadata column names, each drawn
#'   as a top annotation bar. \code{NULL} (default) draws no annotation.
#' @param anno.colors Optional list of colors for \code{anno.var}, one element per
#'   variable (named by variable, or positionally matched to \code{anno.var}). For a
#'   categorical variable supply at least one color per level; for a numeric variable
#'   supply exactly three colors (min, median, max). \code{NULL} entries use defaults.
#' @param assay Assay to pull expression from. Defaults to the object's default assay.
#' @param layer Layer to pull expression from. Default is \code{"scale.data"}.
#' @param hm.limit Numeric vector of three breakpoints for the heatmap color
#'   gradient (low, mid, high). Default is \code{c(-2, 0, 2)}.
#' @param hm.colors Three colors mapped to \code{hm.limit}. Values beyond the range
#'   saturate to the endpoint colors. Default is a blue-white-red gradient.
#' @param cluster.rows Whether to cluster rows (features). Default is \code{FALSE}.
#' @param cluster.columns Whether to cluster columns (cells). Default is \code{FALSE};
#'   clustering overrides the \code{sort.by} ordering.
#' @param show.column.names Whether to show cell names. Default is \code{FALSE}.
#' @param row.font.size Font size for row (feature) names. Default is 12.
#' @param ... Additional arguments passed to \code{\link[ComplexHeatmap]{Heatmap}}.
#'
#' @return A \code{\link[ComplexHeatmap]{Heatmap-class}} object, which draws when
#'   printed. For the horizontal, bottom-placed legend layout use
#'   \code{ComplexHeatmap::draw(ht, heatmap_legend_side = "bottom")}.
#' @export
#' @concept visualization
#'
#' @importFrom SeuratObject LayerData FetchData DefaultAssay
#' @importFrom stats median
#'
PlotAnnotatedHeatmap <- function(
    object,
    features,
    sort.by,
    anno.var = NULL,
    anno.colors = NULL,
    assay = NULL,
    layer = "scale.data",
    hm.limit = c(-2, 0, 2),
    hm.colors = c("#4575b4", "white", "#d73027"),
    cluster.rows = FALSE,
    cluster.columns = FALSE,
    show.column.names = FALSE,
    row.font.size = 12,
    ...) {

  if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) {
    stop("The 'ComplexHeatmap' package is required for PlotAnnotatedHeatmap() but is not ",
         "installed. Please install it from Bioconductor with ",
         "BiocManager::install('ComplexHeatmap').")
  }
  if (!requireNamespace("circlize", quietly = TRUE)) {
    stop("The 'circlize' package is required for PlotAnnotatedHeatmap() but is not ",
         "installed. Please install it with install.packages('circlize').")
  }

  assay <- assay %||% DefaultAssay(object = object)

  # expression matrix (features x cells); LayerData avoids the GetAssayData slot= deprecation
  mat <- LayerData(object = object, assay = assay, layer = layer)
  features <- unique(x = features)
  missing.features <- setdiff(x = features, y = rownames(x = mat))
  if (length(x = missing.features) > 0) {
    warning("The following features are absent from the '", layer, "' layer of the '",
            assay, "' assay and were dropped: ", paste(missing.features, collapse = ", "))
  }
  features <- features[features %in% rownames(x = mat)]
  if (length(x = features) == 0) {
    stop("None of the requested features are present in the '", layer, "' layer of the '",
         assay, "' assay.")
  }
  mat <- as.matrix(x = mat[features, , drop = FALSE])
  cells <- colnames(x = mat)

  # order the columns by the requested metadata variable(s)
  sort.df <- FetchData(object = object, vars = sort.by)[cells, , drop = FALSE]
  ord <- do.call(what = order, args = unname(obj = as.list(x = sort.df)))
  cells <- cells[ord]
  mat <- mat[, cells, drop = FALSE]

  # build the top annotation bars, if requested
  top.anno <- NULL
  if (!is.null(x = anno.var)) {
    if (!is.null(x = anno.colors) && is.null(x = names(x = anno.colors))) {
      names(x = anno.colors) <- anno.var
    }
    anno.df <- FetchData(object = object, vars = anno.var)[cells, , drop = FALSE]
    top.anno <- .scslide_build_heatmap_annotation(anno.df = anno.df, anno.colors = anno.colors)
  }

  col.fun <- circlize::colorRamp2(breaks = hm.limit, colors = hm.colors)

  ComplexHeatmap::Heatmap(
    matrix = mat,
    col = col.fun,
    cluster_rows = cluster.rows,
    cluster_columns = cluster.columns,
    show_column_names = show.column.names,
    row_names_side = "left",
    row_names_gp = grid::gpar(fontsize = row.font.size),
    top_annotation = top.anno,
    heatmap_legend_param = list(title = "Expression"),
    ...
  )
}