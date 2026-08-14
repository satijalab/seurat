#' Sketch Data by Group
#'
#' This function uses sketching methods to downsample high-dimensional single-cell RNA expression data
#' within specified groups/categories, which can help with scalability for large datasets while
#' maintaining representation across different cell types or conditions.
#'
#' @param object A Seurat object.
#' @param group.by A metadata column name to group cells by before sketching. If NULL,
#'   falls back to standard sketching without grouping. Default is NULL.
#' @param assay Assay name. Default is NULL, in which case the default assay of the object is used.
#' @param ncells A positive integer, named vector, or list specifying the number of cells to sample.
#'   \itemize{
#'     \item If a single integer: same number of cells sampled from each group and layer
#'     \item If a named vector with group names: specific number per group
#'     \item If a named list with group names containing layer vectors: specific number per group per layer
#'   }
#'   Default is 1000.
#' @param sketched.assay Sketched assay name. A sketch assay is created or overwritten with the sketch data.
#'   Default is 'sketch'.
#' @param method Sketching method to use. Can be 'LeverageScore' or 'Uniform'.
#'   Default is 'Uniform'.
#' @param var.name A metadata column name to store the leverage scores. Default is 'leverage.score'.
#' @param cells A vector that contains the IDs of the cells that the user wants to keep. If this is set,
#'   these user-defined cells will be directly used to generate the sketch, ignoring group.by parameter.
#' @param over.write Whether to overwrite existing column in the metadata. Default is FALSE.
#' @param seed A positive integer for the seed of the random number generator. Default is 123.
#' @param cast The type to cast the resulting assay to. Default is NULL (no
#'   casting), which preserves the original matrix format (e.g. on-disk
#'   BPCells \code{IterableMatrix} objects stay on disk). Set to
#'   \code{"dgCMatrix"} to force the sketched assay into an in-memory sparse
#'   matrix.
#' @param verbose Print progress and diagnostic messages. Default is TRUE.
#' @param features A character vector of feature names to include in the sketched assay.
#' @param min.cells.per.group Minimum number of cells required per group to perform sketching.
#'   Groups with fewer cells will be included entirely. Default is 10.
#' @param ... Arguments passed to other methods
#'
#' @return A Seurat object with the sketched data added as a new assay. The metadata will contain
#'   information about which groups each sketched cell belongs to.
#'
#' @details
#' When \code{group.by} is specified, the function performs the following steps: \cr
#' 1. Splits cells into groups based on the metadata column \cr
#' 2. Calculates leverage scores (if method = 'LeverageScore') within each group \cr
#' 3. Samples the specified number of cells from each group \cr
#' 4. Combines all sampled cells into the sketched assay
#'
#' This approach ensures that rare cell types or conditions are not underrepresented
#' in the final sketched dataset.
#'
#' @examples
#' \dontrun{
#' # Basic usage with grouping by cell type
#' sketched_obj <- SketchDataByGroup(
#'   object = seurat_obj,
#'   group.by = "cell_type",
#'   ncells = 500
#' )
#'
#' # Different number of cells per group
#' sketched_obj <- SketchDataByGroup(
#'   object = seurat_obj,
#'   group.by = "condition",
#'   ncells = c("control" = 1000, "treatment" = 800)
#' )
#'
#' # Without grouping (falls back to original behavior)
#' sketched_obj <- SketchDataByGroup(
#'   object = seurat_obj,
#'   ncells = 5000
#' )
#' }
#'
#' @importFrom SeuratObject CastAssay Key Key<- Layers DefaultAssay Assays Cells
#' @importFrom rlang abort
#' @importFrom stats setNames
#'
#' @export
#' @concept sketching
#'
SketchDataByGroup <- function(
    object,
    group.by = NULL,
    assay = NULL,
    ncells = 1000L,
    sketched.assay = 'sketch',
    method = c('Uniform', 'LeverageScore'),
    var.name = "leverage.score",
    cells = NULL,
    over.write = FALSE,
    seed = 123L,
    cast = NULL,
    verbose = TRUE,
    features = NULL,
    min.cells.per.group = 10L,
    ...
) {
  # Get internal functions from Seurat
  
  # Input validation and setup
  assay <- assay[1L] %||% DefaultAssay(object = object)
  assay <- match.arg(arg = assay, choices = Assays(object = object))
  method <- match.arg(arg = method)

  if (sketched.assay == assay) {
    abort(message = "Cannot overwrite existing assays")
  }
  if (sketched.assay %in% Assays(object = object)) {
    if (sketched.assay == DefaultAssay(object = object)) {
      DefaultAssay(object = object) <- assay
    }
    object[[sketched.assay]] <- NULL
  }
  if (!over.write) {
    var.name <- CheckMetaVarName(object = object, var.name = var.name)
  }

  # If no grouping specified or cells are pre-defined, use original function logic
  if (!is.null(cells)) {
    if (verbose && !is.null(cells)) {
      message("Using pre-defined cells, ignoring group.by parameter")
    }
    return(SketchData(
      object = object,
      assay = assay,
      ncells = ncells,
      sketched.assay = sketched.assay,
      method = method,
      var.name = var.name,
      cells = cells,
      over.write = over.write,
      seed = seed,
      cast = cast,
      verbose = verbose,
      features = features,
      ...
    ))
  }

  # Validate group.by parameter
  if (!group.by %in% colnames(object[[]])) {
    abort(message = paste0("Column '", group.by, "' not found in object metadata"))
  }

  # Get grouping information
  group_data <- object[[group.by]]
  groups <- as.character(unique(group_data[, 1]))
  groups <- groups[!is.na(groups)]  # Remove NA values

  if (verbose) {
    message(paste0("Found ", length(groups), " groups in '", group.by, "': ",
                   paste(groups, collapse = ", ")))
  }

  # Process ncells parameter for groups
  if (length(ncells) == 1) {
    ncells_per_group <- setNames(rep(ncells, length(groups)), groups)
  } else if (is.null(names(ncells))) {
    if (length(ncells) != length(groups)) {
      abort(message = "When ncells is a vector without names, length must match number of groups")
    }
    ncells_per_group <- setNames(ncells, groups)
  } else {
    # Named vector - check that all groups are covered
    missing_groups <- setdiff(groups, names(ncells))
    if (length(missing_groups) > 0) {
      if (verbose) {
        message(paste0("Using default ncells (", ncells[1], ") for groups: ",
                       paste(missing_groups, collapse = ", ")))
      }
      ncells_per_group <- setNames(rep(ncells[1], length(groups)), groups)
      ncells_per_group[names(ncells)] <- ncells
    } else {
      ncells_per_group <- ncells[groups]
    }
  }

  # Calculate leverage scores if needed
  if (method == 'LeverageScore') {
    if (verbose) {
      message("Calculating Leverage Scores")
    }
    object <- LeverageScore(
      object = object,
      assay = assay,
      var.name = var.name,
      over.write = over.write,
      seed = seed,
      verbose = FALSE,
      features = features,
      ...
    )
  } else if (method == 'Uniform') {
    if (verbose) {
      message("Using uniform sampling")
    }
    object[[var.name]] <- 1
  }

  leverage.score <- object[[var.name]]
  layer.names <- Layers(object = object[[assay]], search = 'data')

  # Sample cells within each group
  all_sampled_cells <- c()

  for (group in groups) {
    if (verbose) {
      message(paste0("Processing group: ", group))
    }

    # Get cells in this group
    group_cells <- rownames(object[[]])[group_data[, 1] == group & !is.na(group_data[, 1])]

    if (length(group_cells) == 0) {
      if (verbose) {
        message(paste0("No cells found in group '", group, "', skipping"))
      }
      next
    }

    if (length(group_cells) < min.cells.per.group) {
      if (verbose) {
        message(paste0("Group '", group, "' has ", length(group_cells),
                       " cells (< min.cells.per.group), including all cells"))
      }
      all_sampled_cells <- c(all_sampled_cells, group_cells)
      next
    }

    # Sample cells from each layer within this group
    group_sampled_cells <- c()
    ncells_group <- ncells_per_group[group]

    for (layer_name in layer.names) {
      layer_cells <- Cells(object[[assay]], layer = layer_name)
      group_layer_cells <- intersect(group_cells, layer_cells)
      if (length(group_layer_cells) == 0) {
        next
      }

      # Calculate how many cells to sample from this layer for this group
      ncells_layer <- min(ncells_group, length(group_layer_cells))

      if (length(group_layer_cells) <= ncells_layer) {
        cells_to_keep <- group_layer_cells
      } else {
        if (!is.null(seed)) {
          set.seed(seed + which(groups == group) + which(layer.names == layer_name))
        }
        cells_to_keep <- sample(
          x = group_layer_cells,
          size = ncells_layer,
          prob = leverage.score[group_layer_cells, ]
        )
      }

      group_sampled_cells <- c(group_sampled_cells, cells_to_keep)
    }

    if (verbose) {
      message(paste0("Sampled ", length(group_sampled_cells), " cells from group '", group, "'"))
    }

    all_sampled_cells <- c(all_sampled_cells, group_sampled_cells)
  }

  if (length(all_sampled_cells) == 0) {
    abort(message = "No cells were sampled. Check your grouping variable and ncells parameters.")
  }

  if (verbose) {
    message(paste0("Total sampled cells: ", length(all_sampled_cells)))
  }

  # Create sketched assay
  sketched <- suppressWarnings(expr = subset(
    x = object[[assay]],
    cells = all_sampled_cells,
    layers = Layers(object = object[[assay]], search = c('counts', 'data'))
  ))

  if (!is.null(x = cast) && inherits(x = sketched, what = 'Assay5')) {
    sketched <- CastAssay(object = sketched, to = cast, ...)
  }

  Key(object = sketched) <- Key(object = sketched.assay, quiet = TRUE)
  object[[sketched.assay]] <- sketched
  DefaultAssay(object = object) <- sketched.assay

  if (verbose) {
    message(paste0("Created sketched assay '", sketched.assay, "' with ",
                   ncol(sketched), " cells"))
  }

  return(object)
}