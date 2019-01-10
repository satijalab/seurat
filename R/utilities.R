#' @include generics.R
#'
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Calculate module scores for featre expression programs in single cells
#'
#' Calculate the average expression levels of each program (cluster) on single cell level,
#' subtracted by the aggregated expression of control featre sets.
#' All analyzed featres are binned based on averaged expression, and the control featres are
#' randomly selected from each bin.
#'
#' @param object Seurat object
#' @param features Featre expression programs in list
#' @param pool List of features to check expression levels agains, defaults to \code{rownames(x = object)}
#' @param nbin Number of bins of aggregate expression levels for all analyzed features
#' @param ctrl Number of control features selected from the same bin per analyzed feature
#' @param k Use feature clusters returned from DoKMeans
#' @param assay Name of assay to use
#' @param name Name for the expression programs
#' @param seed Set a random seed
#'
#' @return Returns a Seurat object with module scores added to object meta data
#'
# @importFrom Hmisc cut2
#' @importFrom ggplot2 cut_number
#' @importFrom Matrix rowMeans colMeans
#'
#' @references Tirosh et al, Science (2016)
#'
#' @export
#'
#' @examples
#' \dontrun{
#' cd_features <- list(c(
#'   'CD79B',
#'   'CD79A',
#'   'CD19',
#'   'CD180',
#'   'CD200',
#'   'CD3D',
#'   'CD2',
#'   'CD3E',
#'   'CD7',
#'   'CD8A',
#'   'CD14',
#'   'CD1C',
#'   'CD68',
#'   'CD9',
#'   'CD247'
#' ))
#' pbmc_small <- AddModuleScore(
#'   object = pbmc_small,
#'   features = cd_features,
#'   ctrl = 5,
#'   name = 'CD_Features'
#' )
#' head(x = pbmc_small[])
#' }
#'
AddModuleScore <- function(
  object,
  features,
  pool = NULL,
  nbin = 24,
  ctrl = 100,
  k = FALSE,
  assay = NULL,
  name = 'Cluster',
  seed = 1
) {
  set.seed(seed = seed)
  assay.old <- DefaultAssay(object = object)
  assay <- assay %||% assay.old
  DefaultAssay(object = object) <- assay
  assay.data <- GetAssayData(object = object)
  features.old <- features
  if (k) {
    .NotYetUsed(arg = 'k')
    features <- list()
    for (i in as.numeric(x = names(x = table(object@kmeans.obj[[1]]$cluster)))) {
      features[[i]] <- names(x = which(x = object@kmeans.obj[[1]]$cluster == i))
    }
    cluster.length <- length(x = features)
  } else {
    if (is.null(x = features)) {
      stop("Missing input feature list")
    }
    features <- lapply(
      X = features,
      FUN = function(x) {
        return(intersect(x = x, y = rownames(x = object)))
      }
    )
    cluster.length <- length(x = features)
  }
  if (!all(LengthCheck(values = features))) {
    warning(paste(
      'Could not find enough features in the object from the following feature lists:',
      paste(names(x = which(x = !LengthCheck(values = features)))),
      'Attempting to match case...'
    ))
    features <- lapply(
      X = features.old,
      FUN = CaseMatch,
      match = rownames(x = object)
    )
  }
  if (!all(LengthCheck(values = features))) {
    stop(paste(
      'The following feature lists do not have enough features present in the object:',
      paste(names(x = which(x = !LengthCheck(values = features)))),
      'exiting...'
    ))
  }
  pool <- pool %||% rownames(x = object)
  data.avg <- Matrix::rowMeans(x = assay.data[pool, ])
  data.avg <- data.avg[order(data.avg)]
  data.cut <- cut_number(x = data.avg, n = nbin, labels = FALSE, right = FALSE)
  # data.cut <- as.numeric(x = Hmisc::cut2(x = data.avg, m = round(x = length(x = data.avg) / (nbin + 1))))
  names(x = data.cut) <- names(x = data.avg)
  ctrl.use <- vector(mode = "list", length = cluster.length)
  for (i in 1:cluster.length) {
    features.use <- features[[i]]
    for (j in 1:length(x = features.use)) {
      ctrl.use[[i]] <- c(
        ctrl.use[[i]],
        names(x = sample(
          x = data.cut[which(x = data.cut == data.cut[features.use[j]])],
          size = ctrl,
          replace = FALSE
        ))
      )
    }
  }
  ctrl.use <- lapply(X = ctrl.use, FUN = unique)
  ctrl.scores <- matrix(
    data = numeric(length = 1L),
    nrow = length(x = ctrl.use),
    ncol = ncol(x = object)
  )
  for (i in 1:length(ctrl.use)) {
    features.use <- ctrl.use[[i]]
    ctrl.scores[i, ] <- Matrix::colMeans(x = assay.data[features.use, ])
  }
  features.scores <- matrix(
    data = numeric(length = 1L),
    nrow = cluster.length,
    ncol = ncol(x = object)
  )
  for (i in 1:cluster.length) {
    features.use <- features[[i]]
    data.use <- assay.data[features.use, , drop = FALSE]
    features.scores[i, ] <- Matrix::colMeans(x = data.use)
  }
  features.scores.use <- features.scores - ctrl.scores
  rownames(x = features.scores.use) <- paste0(name, 1:cluster.length)
  features.scores.use <- as.data.frame(x = t(x = features.scores.use))
  rownames(x = features.scores.use) <- colnames(x = object)
  object[[colnames(x = features.scores.use)]] <- features.scores.use
  gc(verbose = FALSE)
  DefaultAssay(object = object) <- assay.old
  return(object)
}

#' Averaged feature expression by identity class
#'
#' Returns expression for an 'average' single cell in each identity class
#'
#' Output is in log-space when \code{return.seurat = TRUE}, otherwise it's in non-log space.
#' Averaging is done in non-log space.
#'
#' @param object Seurat object
#' @param assay Which assay to use. Default is all assays.
#' @param features Features to analyze. Default is all features in the assay.
#' @param return.seurat Whether to return the data as a Seurat object. Default is FALSE.
#' @param add.ident Place an additional label on each cell prior to averaging (very useful if you want to observe cluster averages, separated by replicate, for example).
#' @param use.scale Use scaled values for gene expression
#' @param use.counts Use count values for gene expression
#' @param verbose Print messages and show progress bar
#' @param ... Arguments to be passed to methods such as \code{\link{CreateSeuratObject}}
#'
#' @return Returns a matrix with genes as rows, identity classes as columns.
#' If return.seurat is TRUE, returns an object of class \code{\link{Seurat}}.
#'
#' @export
#'
#' @examples
#' head(AverageExpression(object = pbmc_small))
#'
AverageExpression <- function(
  object,
  assay = NULL,
  features = NULL,
  return.seurat = FALSE,
  add.ident = NULL,
  use.scale = FALSE,
  use.counts = FALSE,
  verbose = TRUE,
  ...
) {
  assay <- assay %||% DefaultAssay(object)
  ident.orig <- Idents(object)
  orig.levels <- levels(x = Idents(object))
  ident.new <- c()
  if (!all(assay %in% names(object@assays))) {
    assay <- assay[assay %in% names(object@assays)]
    if (length(assay) == 0) {
      stop("None of the requested assays are present in the object")
    } else {
      warning("Requested assays that do not exist in object. Proceeding with existing assays only.")
    }
  }
  if (!is.null(x = add.ident)) {
    new.data <- FetchData(object = object, vars = add.ident)
    new.ident <- paste(
      Idents(object)[rownames(x = new.data)],
      new.data[, 1],
      sep = '_'
    )
    Idents(object, cells = rownames(new.data)) <- new.ident
  }
  slot.use <- "data"
  fxn.average <- function(x) mean(expm1(x))
  if (use.scale) {
    slot.use <- "scale.data"
    fxn.average <- mean
  }
  if (use.counts) {
    slot.use <- "counts"
    fxn.average <- mean
  }
  data.return <- list()
  for (i in 1:length(x = assay)) {
    data.use <- GetAssayData(
      object = object,
      assay = assay[i],
      slot = slot.use
    )
    features.assay <- features
    if (length(x = intersect(x = features, y = rownames(x = data.use))) < 1 ) {
      features.assay <- rownames(x = data.use)
    }
    data.all <- data.frame(row.names = features.assay)
    for (j in levels(x = Idents(object))) {
      temp.cells <- WhichCells(object = object, idents = j)
      features.assay <- unique(x = intersect(x = features.assay, y = rownames(x = data.use)))
      if (length(x = temp.cells) == 1) {
        data.temp <- (data.use[features.assay, temp.cells])
        # transform data if needed (alternative: apply fxn.average to single value above)
        if (!(use.scale | use.counts)) { # equivalent: slot.use == "data"
          data.temp <- expm1(data.temp)
        }
      }
      if (length(x = temp.cells) > 1 ) {
        data.temp <- apply(
          X = data.use[features.assay, temp.cells, drop = FALSE],
          MARGIN = 1,
          FUN = fxn.average
        )
      }
      data.all <- cbind(data.all, data.temp)
      colnames(x = data.all)[ncol(x = data.all)] <- j
      if (verbose) {
        message(paste("Finished averaging", assay[i], "for cluster", j))
      }
      if (i == 1) {
        ident.new <- c(ident.new, as.character(x = ident.orig[temp.cells[1]]))
      }
    }
    names(x = ident.new) <- levels(x = Idents(object))
    data.return[[i]] <- data.all
    names(x = data.return)[i] <- assay[[i]]
  }

  if (return.seurat) {
    toRet <- CreateSeuratObject(
      counts = data.return[[1]],
      project = "Average",
      assay = names(x = data.return)[1],
      ...
    )

    #for multimodal data
    if (length(x = data.return) > 1) {
      for (i in 2:length(x = data.return)) {
        toRet[[names(x = data.return)[i]]] <- CreateAssayObject(counts = data.return[[i]])
      }
    }
    Idents(toRet, cells = colnames(toRet)) <- ident.new[colnames(toRet)]
    Idents(toRet) <- factor(
      x = Idents(toRet),
      levels = as.character(x = orig.levels),
      ordered = TRUE
    )

    # finish setting up object if it is to be returned
    toRet <- NormalizeData(toRet, verbose = verbose)
    toRet <- ScaleData(toRet, verbose = verbose)

    return(toRet)
  } else {
    return(data.return)
  }
}

#' Match the case of character vectors
#'
#' @param search A vector of search terms
#' @param match A vector of characters whose case should be matched
#'
#' @return Values from search present in match with the case of match
#'
#' @export
#'
#' @examples
#' cd_genes <- c('Cd79b', 'Cd19', 'Cd200')
#' CaseMatch(search = cd_genes, match = rownames(x = pbmc_small))
#'
CaseMatch <- function(search, match) {
  search.match <- sapply(
    X = search,
    FUN = function(s) {
      return(grep(
        pattern = paste0('^', s, '$'),
        x = match,
        ignore.case = TRUE,
        perl = TRUE,
        value = TRUE
      ))
    }
  )
  return(unlist(x = search.match))
}

#' Score cell cycle phases
#'
#' @param object A Seurat object
#' @param s.features A vector of features associated with S phase
#' @param g2m.features A vector of features associated with G2M phase
#' @param set.ident If true, sets identity to phase assignments
#' Stashes old identities in 'old.ident'
#'
#' @return A Seurat object with the following columns added to object meta data: S.Score, G2M.Score, and Phase
#'
#' @seealso \code{AddModuleScore}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # pbmc_small doesn't have any cell-cycle genes
#' # To run CellCycleScoring, please use a dataset with cell-cycle genes
#' # An example is available at http://satijalab.org/seurat/cell_cycle_vignette.html
#' pbmc_small <- CellCycleScoring(
#'   object = pbmc_small,
#'   g2m.features = cc.genes$g2m.genes,
#'   s.features = cc.genes$s.genes
#' )
#' head(x = pbmc_small@meta.data)
#' }
#'
CellCycleScoring <- function(
  object,
  s.features,
  g2m.features,
  set.ident = FALSE
) {
  name <- 'Cell Cycle'
  features <- list('S.Score' = s.features, 'G2M.Score' = g2m.features)
  object.cc <- AddModuleScore(
    object = object,
    features = features,
    name = name,
    ctrl = min(vapply(X = features, FUN = length, FUN.VALUE = numeric(length = 1)))
  )
  cc.columns <- grep(pattern = name, x = colnames(x = object.cc[[]]), value = TRUE)
  cc.scores <- object.cc[[cc.columns]]
  rm(object.cc)
  gc(verbose = FALSE)
  assignments <- apply(
    X = cc.scores,
    MARGIN = 1,
    FUN = function(scores, first = 'S', second = 'G2M', null = 'G1') {
      if (all(scores < 0)) {
        return(null)
      } else {
        if (length(which(x = scores == max(scores))) > 1) {
          return('Undecided')
        } else {
          return(c(first, second)[which(x = scores == max(scores))])
        }
      }
    }
  )
  cc.scores <- merge(x = cc.scores, y = data.frame(assignments), by = 0)
  colnames(x = cc.scores) <- c('rownames', 'S.Score', 'G2M.Score', 'Phase')
  rownames(x = cc.scores) <- cc.scores$rownames
  cc.scores <- cc.scores[, c('S.Score', 'G2M.Score', 'Phase')]
  object[[colnames(x = cc.scores)]] <- cc.scores
  if (set.ident) {
    object[['old.ident']] <- Idents(object = object)
    Idents(object = object) <- 'Phase'
  }
  return(object)
}

#' Slim down a multi-species expression matrix, when only one species is primarily of interenst.
#'
#' Valuable for CITE-seq analyses, where we typically spike in rare populations of 'negative control' cells from a different species.
#'
#' @param object A UMI count matrix. Should contain rownames that start with
#' the ensuing arguments prefix.1 or prefix.2
#' @param prefix The prefix denoting rownames for the species of interest.
#' Default is "HUMAN_". These rownames will have this prefix removed in the returned matrix.
#' @param controls The prefix denoting rownames for the species of 'negative
#' control' cells. Default is "MOUSE_".
#' @param ncontrols How many of the most highly expressed (average) negative
#' control features (by default, 100 mouse genes), should be kept? All other
#' rownames starting with prefix.2 are discarded.
#'
#' @return A UMI count matrix. Rownames that started with \code{prefix} have this
#' prefix discarded. For rownames starting with \code{controls}, only the
#' \code{ncontrols} most highly expressed features are kept, and the
#' prefix is kept. All other rows are retained.
#'
#' @importFrom Matrix rowSums
#'
#' @export
#'
#' @examples
#' \dontrun{
#' cbmc.rna.collapsed <- CollapseSpeciesExpressionMatrix(cbmc.rna)
#' }
#'
CollapseSpeciesExpressionMatrix <- function(
  object,
  prefix = "HUMAN_",
  controls = "MOUSE_",
  ncontrols = 100
) {
  features <- grep(pattern = prefix, x = rownames(x = object), value = TRUE)
  controls <- grep(pattern = controls, x = rownames(x = object), value = TRUE)
  others <- setdiff(x = rownames(x = object), y = c(features, controls))
  controls <- rowSums(x = object[controls, ])
  controls <- names(x = head(
    x = sort(x = controls, decreasing = TRUE),
    n = ncontrols
  ))
  object <- object[c(features, controls, others), ]
  rownames(x = object) <- gsub(
    pattern = prefix,
    replacement = '',
    x = rownames(x = object)
  )
  return(object)
}

#' Run a custom distance function on an input data matrix
#'
#' @author Jean Fan
#'
#' @param my.mat A matrix to calculate distance on
#' @param my.function A function to calculate distance
#' @param ... Extra parameters to my.function
#'
#' @return A distance matrix
#'
#' @importFrom stats as.dist
#'
#' @export
#'
#' @examples
#' # Define custom distance matrix
#' manhattan.distance <- function(x, y) return(sum(abs(x-y)))
#'
#' input.data <- GetAssayData(pbmc_small, assay.type = "RNA", slot = "scale.data")
#' cell.manhattan.dist <- CustomDistance(input.data, manhattan.distance)
#'
CustomDistance <- function(my.mat, my.function, ...) {
  n <- ncol(x = my.mat)
  mat <- matrix(data = 0, ncol = n, nrow = n)
  colnames(x = mat) <- rownames(x = mat) <- colnames(x = my.mat)
  for (i in 1:nrow(x = mat)) {
    for(j in 1:ncol(x = mat)) {
      mat[i,j] <- my.function(my.mat[, i], my.mat[, j], ...)
    }
  }
  return(as.dist(m = mat))
}

#' Calculate the mean of logged values
#'
#' Calculate mean of logged values in non-log space (return answer in log-space)
#'
#' @param x A vector of values
#'
#' @return Returns the mean in log-space
#'
#' @export
#'
#' @examples
#' ExpMean(x = c(1, 2, 3))
#'
ExpMean <- function(x) {
  return(log(x = mean(x = exp(x = x) - 1) + 1))
}

#' Export Seurat object for UCSC cell browser
#'
#' @param object Seurat object
#' @param dir path to directory where to save exported files. These are:
#' exprMatrix.tsv, tsne.coords.tsv, meta.tsv, markers.tsv and a default cellbrowser.conf
#' @param dataset.name name of the dataset. Defaults to Seurat project name
#' @param reductions vector of reduction names to export
#' @param markers.file path to file with marker genes
#' @param cluster.field name of the metadata field containing cell cluster
#' @param cb.dir path to directory where to create UCSC cellbrowser static
#' website content root, e.g. an index.html, .json files, etc. These files
#' can be copied to any webserver. If this is specified, the cellbrowser
#' package has to be accessible from R via reticulate.
#' @param port on which port to run UCSC cellbrowser webserver after export
#' @param skip.expr.matrix whether to skip exporting expression matrix
#' @param skip.metadata whether to skip exporting metadata
#' @param skip.reductions whether to skip exporting reductions
#' @param ... specifies the metadata fields to export. To supply field with
#' human readable name, pass name as \code{field="name"} parameter.
#'
#' @return This function exports Seurat object as a set of tsv files
#' to \code{dir} directory, copying the \code{markers.file} if it is
#' passed. It also creates the default \code{cellbrowser.conf} in the
#' directory. This directory could be read by \code{cbBuild} to
#' create a static website viewer for the dataset. If \code{cb.dir}
#' parameter is passed, the function runs \code{cbBuild} (if it is
#' installed) to create this static website in \code{cb.dir} directory.
#' If \code{port} parameter is passed, it also runs the webserver for
#' that directory and opens a browser.
#'
#' @author Maximilian Haeussler, Nikolay Markov
#'
#' @importFrom utils browseURL
#' @importFrom reticulate py_module_available import
#' @importFrom tools file_ext
#'
#' @export
#'
#' @examples
#' \dontrun{
#' ExportToCellbrowser(object = pbmc_small, dataset.name = "PBMC", dir = "out")
#' }
#'
ExportToCellbrowser <- function(
  object,
  dir,
  dataset.name = slot(object = object, name = 'project.name'),
  reductions = "tsne",
  markers.file = NULL,
  cluster.field = "Cluster",
  cb.dir = NULL,
  port = NULL,
  skip.expr.matrix = FALSE,
  skip.metadata = FALSE,
  skip.reductions = FALSE,
  ...
) {
  vars <- c(...)
  if (is.null(x = vars)) {
    vars <- c("nCount_RNA", "nFeature_RNA")
    if (length(x = levels(x = Idents(object = object))) > 1) {
      vars <- c(vars, cluster.field)
      names(x = vars) <- c("", "", "ident")
    }
  }
  names(x = vars) <- names(x = vars) %||% vars
  names(x = vars) <- sapply(
    X = 1:length(x = vars),
    FUN = function(i) {
      return(ifelse(
        test = nchar(x = names(x = vars)[i]) > 0,
        yes = names(x = vars[i]),
        no = vars[i]
      ))
    }
  )
  if (!is.null(x = port) && is.null(x = cb.dir)) {
    stop("cb.dir parameter is needed when port is set")
  }
  if (!dir.exists(paths = dir)) {
    dir.create(path = dir)
  }
  if (!dir.exists(paths = dir)) {
    stop("Output directory ", dir, " cannot be created or is a file")
  }
  if (dataset.name == "SeuratProject") {
    warning("Using default project name means that you may overwrite project with the same name in the cellbrowser html output folder")
  }
  order <- colnames(x = object)
  enum.fields <- c()
  # Export expression matrix:
  if (!skip.expr.matrix) {
    # Relatively memory inefficient - maybe better to convert to sparse-row and write in a loop, row-by-row?
    df <- as.data.frame(x = as.matrix(x = GetAssayData(object = object)))
    df <- data.frame(gene = rownames(x = object), df, check.names = FALSE)
    gzPath <- file.path(dir, "exprMatrix.tsv.gz")
    z <- gzfile(gzPath, "w")
    message("Writing expression matrix to ", gzPath)
    write.table(x = df, sep = "\t", file = z, quote = FALSE, row.names = FALSE)
    close(con = z)
  }
  # Export cell embeddings
  embeddings.conf <- c()
  for (reduction in reductions) {
    if (!skip.reductions) {
      df <- Embeddings(object = object, reduction = reduction)
      if (ncol(x = df) > 2) {
        warning(
          'Embedding ',
          reduction,
          ' has more than 2 coordinates, taking only the first 2'
        )
        df <- df[, 1:2]
      }
      colnames(x = df) <- c("x", "y")
      df <- data.frame(cellId = rownames(x = df), df)
      fname <- file.path(dir, paste0(reduction, '.coords.tsv'))
      message("Writing embeddings to ", fname)
      write.table(
        x = df[order, ],
        sep = "\t",
        file = fname,
        quote = FALSE,
        row.names = FALSE
      )
    }
    conf <- sprintf(
      '{"file": "%s.coords.tsv", "shortLabel": "Seurat %1$s"}',
      reduction
    )
    embeddings.conf <- c(embeddings.conf, conf)
  }
  # Export metadata
  df <- data.frame(row.names = rownames(x = object[[]]))
  df <- FetchData(object = object, vars = names(x = vars))
  colnames(x = df) <- vars
  enum.fields <- Filter(
    f = function(name) {!is.numeric(x = df[[name]])},
    x = vars
  )
  if (!skip.metadata) {
    fname <- file.path(dir, "meta.tsv")
    message("Writing meta data to ", fname)
    df <- data.frame(Cell = rownames(x = df), df, check.names = FALSE)
    write.table(
      x = df[order, ],
      sep = "\t",
      file = fname,
      quote = FALSE,
      row.names = FALSE
    )
  }
  # Export markers
  markers.string <- ''
  if (!is.null(x = markers.file)) {
    ext <- file_ext(x = markers.file)
    fname <- paste0("markers.", ext)
    file.copy(from = markers.file, to = file.path(dir, fname))
    markers.string <- sprintf(
      'markers = [{"file": "%s", "shortLabel": "Seurat Cluster Markers"}]',
      fname
    )
  }
  config <- c(
    'name="%s"',
    'shortLabel="%1$s"',
    'exprMatrix="exprMatrix.tsv.gz"',
    '#tags = ["10x", "smartseq2"]',
    'meta="meta.tsv"',
    '# possible values: "gencode-human", "gencode-mouse", "symbol" or "auto"',
    'geneIdType="auto"',
    'clusterField="%s"',
    'labelField="%2$s"',
    'enumFields=%s',
    '%s',
    'coords=%s'
  )
  config <- paste(config, collapse = '\n')
  enum.string <- paste0(
    "[",
    paste(paste0('"', enum.fields, '"'), collapse = ", "),
    "]"
  )
  coords.string <- paste0(
    "[",
    paste(embeddings.conf, collapse = ",\n"),
    "]"
  )
  config <- sprintf(
    config,
    dataset.name,
    cluster.field,
    enum.string,
    markers.string,
    coords.string
  )
  fname <- file.path(dir, "cellbrowser.conf")
  if (file.exists(fname)) {
    message(
      "`cellbrowser.conf` already exists in target directory, refusing to ",
      "overwrite it"
    )
  } else {
    cat(config, file = fname)
  }
  message("Prepared cellbrowser directory ", dir)
  if (!is.null(x = cb.dir)) {
    if (!py_module_available(module = "cellbrowser")) {
      stop(
        "The Python package `cellbrowser` is required to prepare and run ",
        "Cellbrowser. Please install it ",
        "on the Unix command line with `sudo pip install cellbrowser` (if root) ",
        "or `pip install cellbrowser --user` (as a non-root user). ",
        "To adapt the Python that is used, you can either set the env. variable RETICULATE_PYTHON ",
        "or do `require(reticulate) and use one of these functions: use_python(), use_virtualenv(), use_condaenv(). ",
        "See https://rstudio.github.io/reticulate/articles/versions.html; ",
        "at the moment, R's reticulate is using this Python: ",
        import(module = 'sys')$executable,
        ". "
      )
    }
    if (!is.null(x = port)) {
      port <- as.integer(x = port)
    }
    message("Converting cellbrowser directory to html/json files")
    cb <- import(module = "cellbrowser")
    cb$cellbrowser$build(dir, cb.dir)
    if (!is.null(port)) {
      message("Starting http server")
      cb$cellbrowser$stop()
      cb$cellbrowser$serve(cb.dir, port)
      Sys.sleep(time = 0.4)
      browseURL(url = paste0("http://localhost:", port))
    }
  }
}

#' Calculate the standard deviation of logged values
#'
#' Calculate SD of logged values in non-log space (return answer in log-space)
#'
#' @param x A vector of values
#'
#' @return Returns the standard deviation in log-space
#'
#' @importFrom stats sd
#'
#' @export
#'
#' @examples
#' ExpSD(x = c(1, 2, 3))
#'
ExpSD <- function(x) {
  return(log1p(x = sd(x = expm1(x = x))))
}

#' Calculate the variance of logged values
#'
#' Calculate variance of logged values in non-log space (return answer in
#' log-space)
#'
#' @param x A vector of values
#'
#' @return Returns the variance in log-space
#'
#' @importFrom stats var
#'
#' @export
#'
#' @examples
#' ExpVar(x = c(1, 2, 3))
#'
ExpVar <- function(x) {
  return(log1p(x = var(x = expm1(x = x))))
}

#' Calculate the variance to mean ratio of logged values
#'
#' Calculate the variance to mean ratio (VMR) in non-logspace (return answer in
#' log-space)
#'
#' @param x A vector of values
#'
#' @return Returns the VMR in log-space
#'
#' @importFrom stats var
#'
#' @export
#'
#' @examples
#' LogVMR(x = c(1, 2, 3))
#'
LogVMR <- function(x) {
  return(log(x = var(x = exp(x = x) - 1) / mean(x = exp(x = x) - 1)))
}

#' Aggregate expression of multiple features into a single feature
#'
#' Calculates relative contribution of each feature to each cell
#' for given set of features.
#'
#' @param object A Seurat object
#' @param features List of features to aggregate
#' @param meta.name Name of column in metadata to store metafeature
#' @param cells List of cells to use (default all cells)
#' @param assay Which assay to use
#' @param slot Which slot to take data from (default data)
#'
#' @return Returns a \code{Seurat} object with metafeature stored in objct metadata
#'
#' @importFrom Matrix rowSums colMeans
#'
#' @export
#'
#' @examples
#' pbmc_small <- MetaFeature(
#'   object = pbmc_small,
#'   features = c("LTB", "EAF2"),
#'   meta.name = 'var.aggregate'
#' )
#' head(pbmc_small[[]])
#'
MetaFeature <- function(
  object,
  features,
  meta.name = 'metafeature',
  cells = NULL,
  assay = NULL,
  slot = 'data'
) {
  cells <- cells %||% colnames(x = object)
  assay <- assay %||% DefaultAssay(object = object)
  newmat <- GetAssayData(object = object, assay = assay, slot = slot)
  newmat <- newmat[features, cells]
  if (slot == 'scale.data') {
    newdata <- Matrix::colMeans(newmat)
  } else {
    rowtotals <- Matrix::rowSums(newmat)
    newmat <- newmat / rowtotals
    newdata <- Matrix::colMeans(newmat)
  }
  object[[meta.name]] <- newdata
  return(object)
}

#' Apply a ceiling and floor to all values in a matrix
#'
#' @param data Matrix or data frame
#' @param min all values below this min value will be replaced with min
#' @param max all values above this max value will be replaced with max
#' @return Returns matrix after performing these floor and ceil operations
#' @export
#'
#' @examples
#' mat <- matrix(data = rbinom(n = 25, size = 20, prob = 0.2 ), nrow = 5)
#' mat
#' MinMax(data = mat, min = 4, max = 5)
#'
MinMax <- function(data, min, max) {
  data2 <- data
  data2[data2 > max] <- max
  data2[data2 < min] <- min
  return(data2)
}

#' Stop Cellbrowser web server
#'
#' @importFrom reticulate py_module_available
#' @importFrom reticulate import
#'
#' @export
#'
#' @examples
#' \dontrun{
#' StopCellbrowser()
#' }
#'
StopCellbrowser <- function() {
  if (py_module_available(module = "cellbrowser")) {
    cb <- import(module = "cellbrowser")
    cb$cellbrowser$stop()
  } else {
    stop("The `cellbrowser` package is not available in the Python used by R's reticulate")
  }
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Methods for Seurat-defined generics
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @importFrom methods as
#' @importClassesFrom Matrix dgCMatrix
#'
#' @rdname as.sparse
#' @export
#' @method as.sparse data.frame
#'
as.sparse.data.frame <- function(x, ...) {
  return(as(object = as.matrix(x = x), Class = 'dgCMatrix'))
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Methods for R-defined generics
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @inheritParams base::as.data.frame
#'
#' @return \code{as.data.frame.Matrix}: A data frame representation of the S4 Matrix
#'
#' @importFrom Matrix as.matrix
#'
#' @rdname as.sparse
#' @export
#' @method as.data.frame Matrix
#'
as.data.frame.Matrix <- function(
  x,
  row.names = NULL,
  optional = FALSE,
  ...,
  stringsAsFactors = default.stringsAsFactors()
) {
  return(as.data.frame(
    x = as.matrix(x = x),
    row.names = row.names,
    optional = optional,
    stringsAsFactors = stringsAsFactors,
    ...
  ))
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Internal
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Set a default value if an object is null
#
# @param lhs An object to set if it's null
# @param rhs The value to provide if x is null
#
# @return rhs if lhs is null, else lhs
#
# @author Hadley Wickham
# @references https://adv-r.hadley.nz/functions.html#missing-arguments
#
`%||%` <- function(lhs, rhs) {
  if (!is.null(x = lhs)) {
    return(lhs)
  } else {
    return(rhs)
  }
}

# Set a default value if an object is NOT null
#
# @param lhs An object to set if it's NOT null
# @param rhs The value to provide if x is NOT null
#
# @return lhs if lhs is null, else rhs
#
# @author Hadley Wickham
# @references https://adv-r.hadley.nz/functions.html#missing-arguments
#
`%iff%` <- function(lhs, rhs) {
  if (!is.null(x = lhs)) {
    return(rhs)
  } else {
    return(lhs)
  }
}

# Generate chunk points
#
# @param dsize How big is the data being chunked
# @param csize How big should each chunk be
#
# @return A matrix where each column is a chunk, row 1 is start points, row 2 is end points
#
ChunkPoints <- function(dsize, csize) {
  return(vapply(
    X = 1L:ceiling(x = dsize / csize),
    FUN = function(i) {
      return(c(
        start = (csize * (i - 1L)) + 1L,
        end = min(csize * i, dsize)
      ))
    },
    FUN.VALUE = numeric(length = 2L)
  ))
}

# L2 normalize the columns (or rows) of a given matrix
# @param mat Matrix to cosine normalize
# @param MARGIN Perform normalization over rows (1) or columns (2)
#
#
# @return returns l2-normalized matrix
#
#
L2Norm <- function(mat, MARGIN = 1){
  normalized <- sweep(
    x = mat,
    MARGIN = MARGIN,
    STATS = apply(
      X = mat,
      MARGIN = MARGIN,
      FUN = function(x){
        sqrt(x = sum(x ^ 2))
      }
    ),
    FUN = "/"
  )
  normalized[!is.finite(x = normalized)] <- 0
  return(normalized)
}

# Check the use of ...
#
# @param ... Arguments passed to a function that fall under ...
# @param fxns A list/vector of functions or function names
#
# @return ...
#
#' @importFrom utils argsAnywhere getAnywhere
#
# @examples
#
CheckDots <- function(..., fxns = NULL) {
  args.names <- names(x = c(...))
  if (length(x = c(...)) == 0) {
    return(invisible(x = NULL))
  }
  if (is.null(x = args.names)) {
    stop("No named arguments passed")
  }
  fxn.args <- lapply(
    X = fxns,
    FUN = function(x) {
      if (is.character(x = x)) {
        if (any(grepl(pattern = 'UseMethod', x = as.character(x = getAnywhere(x = x))))) {
          x <- paste0(x, '.default')
        }
        x <- argsAnywhere(x = x)
      }
      if (is.function(x = x)) {
        return(names(x = formals(fun = x)))
      } else {
        stop("CheckDots only works on characters or functions, not ", class(x = x))
      }
    }
  )
  if (any(grepl(pattern = '...', x = fxn.args, fixed = TRUE))) {
    dotted <- grepl(pattern = '...', x = fxn.args, fixed = TRUE)
    dfxn <- fxns[dotted]
    if (any(sapply(X = dfxn, FUN = is.character))) {
      message(
        "The following functions accept the dots: ",
        paste(Filter(f = is.character, x = dfxn), collapse = ', ')
      )
      dfxn <- Filter(f = Negate(f = is.character), x = dfxn)
      if (length(x = dfxn) > 0) {
        message(
          "In addition, there are ",
          length(x = dfxn),
          " unnamed function(s) that accept the dots"
        )
      }
    } else {
      message("There are ", length(x = dfxn), ' functions provided that accept the dots')
    }
  } else {
    unused <- Filter(
      f = function(x) {
        return(!x %in% unlist(x = fxn.args))
      },
      x = args.names
    )
    if (length(x = unused) > 0) {
      warning(
        "The following arguments are not used: ",
        paste(unused, collapse = ', '),
        call. = FALSE,
        immediate. = TRUE
      )
    }
  }
}

# Check a list of objects for duplicate cell names
#
# @param object.list List of Seurat objects
# @param verbose Print message about renaming
# @param stop Error out if any duplicate names exist
#
# @return Returns list of objects with duplicate cells renamed to be unique
#
CheckDuplicateCellNames <- function(object.list, verbose = TRUE, stop = FALSE) {
  cell.names <- unlist(
    x = sapply(
      X = 1:length(x = object.list),
      FUN = function(x) Cells(object = object.list[[x]])
      )
    )
  if (any(duplicated(x = cell.names))) {
    if (stop) {
      stop("Duplicate cell names present across objects provided.")
    }
    if (verbose) {
      warning("Some cell names are duplicated across objects provided. Renaming to enforce unique cell names.")
    }
    object.list <- sapply(
      X = 1:length(x = object.list),
      FUN = function(x) RenameCells(object = object.list[[x]], new.names = paste0(Cells(object = object.list[[x]]), "_", x))
    )
  }
  return(object.list)
}

# Extract delimiter information from a string.
#
# Parses a string (usually a cell name) and extracts fields based on a delimiter
#
# @param string String to parse.
# @param field Integer(s) indicating which field(s) to extract. Can be a vector multiple numbers.
# @param delim Delimiter to use, set to underscore by default.
#
# @return A new string, that parses out the requested fields, and (if multiple), rejoins them with the same delimiter
#
# @export
#
# @examples
# ExtractField(string = 'Hello World', field = 1, delim = '_')
#
ExtractField <- function(string, field = 1, delim = "_") {
  fields <- as.numeric(x = unlist(x = strsplit(x = as.character(x = field), split = ",")))
  if (length(x = fields) == 1) {
    return(strsplit(x = string, split = delim)[[1]][field])
  }
  return(paste(strsplit(x = string, split = delim)[[1]][fields], collapse = delim))
}

# Interleave vectors together
#
# @param ... Vectors to be interleaved
#
# @return A vector with the values from each vector in ... interleaved
#
Interleave <- function(...) {
  return(as.vector(x = t(x = as.data.frame(x = list(...)))))
}

# Check if a matrix is empty
#
# Takes a matrix and asks if it's empty (either 0x0 or 1x1 with a value of NA)
#
# @param x A matrix
#
# @return Whether or not \code{x} is empty
#
IsMatrixEmpty <- function(x) {
  matrix.dims <- dim(x = x)
  matrix.na <- all(matrix.dims == 1) && all(is.na(x = x))
  return(all(matrix.dims == 0) || matrix.na)
}

# Check the length of components of a list
#
# @param values A list whose components should be checked
# @param cutoff A minimum value to check for
#
# @return a vector of logicals
#
LengthCheck <- function(values, cutoff = 0) {
  return(vapply(
    X = values,
    FUN = function(x) {
      return(length(x = x) > cutoff)
    },
    FUN.VALUE = logical(1)
  ))
}

# Logs a command run, storing the name, timestamp, and argument list. Stores in
# the Seurat object
# @param object Name of Seurat object
#
# @return returns the Seurat object with command stored
#
LogSeuratCommand <- function(object, return.command = FALSE) {
  time.stamp <- Sys.time()
  #capture function name
  command.name <- as.character(deparse(sys.calls()[[sys.nframe() - 1]]))
  command.name <- gsub(pattern = ".Seurat", replacement = "", x = command.name)
  call.string <- command.name
  command.name <- ExtractField(string = command.name, field = 1, delim = "\\(")
  #capture function arguments
  argnames <- names(x = formals(fun = sys.function(which = sys.parent(n = 1))))
  argnames <- grep(pattern = "object", x = argnames, invert = TRUE, value = TRUE)
  argnames <- grep(pattern = "\\.\\.\\.", x = argnames, invert = TRUE, value = TRUE)
  params <- list()
  p.env <- parent.frame(n = 1)
  argnames <- intersect(x = argnames, y = ls(name = p.env))
  # fill in params list
  for (arg in argnames) {
    param_value <- get(x = arg, envir = p.env)
    #TODO Institute some check of object size?
    params[[arg]] <- param_value
  }
  # check if function works on the Assay and/or the DimReduc Level
  assay <- params[["assay"]]
  reduction <- params[["reduction"]]
  # rename function name to include Assay/DimReduc info
  command.name <- paste(command.name, assay, reduction, sep = '.')
  command.name <- sub(pattern = "[\\.]+$", replacement = "", x = command.name, perl = TRUE)
  command.name <- sub(pattern = "\\.\\.", replacement = "\\.", x = command.name, perl = TRUE)
  if (return.command) {
    return(command.name)
  }
  # store results
  seurat.command <- new(
    Class = 'SeuratCommand',
    name = command.name,
    params = params,
    time.stamp = time.stamp,
    call.string = call.string
  )
  object[[command.name]] <- seurat.command
  return(object)
}

# Independently shuffle values within each row of a matrix
#
# Creates a matrix where correlation structure has been removed, but overall values are the same
#
# @param x Matrix to shuffle
#
# @return Returns a scrambled matrix, where each row is shuffled independently
#
#' @importFrom stats runif
#
# @export
#
# @examples
# mat <- matrix(data = rbinom(n = 25, size = 20, prob = 0.2 ), nrow = 5)
# mat
# MatrixRowShuffle(x = mat)
#
MatrixRowShuffle <- function(x) {
  x2 <- x
  x2 <- t(x = x)
  ind <- order(c(col(x = x2)), runif(n = length(x = x2)))
  x2 <- matrix(
    data = x2[ind],
    nrow = nrow(x = x),
    ncol = ncol(x = x),
    byrow = TRUE
  )
  return(x2)
}

# Reverse the vector x and return the value at the Nth index. If N is larger
# than the length of the vector, return the last value in the reversed vector.
#
# @param x vector of interest
# @param N index in reversed vector
#
# @return returns element at given index
#
MaxN <- function(x, N = 2){
  len <- length(x)
  if (N > len) {
    warning('N greater than length(x).  Setting N=length(x)')
    N <- length(x)
  }
  sort(x, partial = len - N + 1)[len - N + 1]
}

# Melt a data frame
#
# @param x A data frame
#
# @return A molten data frame
#
Melt <- function(x) {
  if (!is.data.frame(x = x)) {
    x <- as.data.frame(x = x)
  }
  return(data.frame(
    rows = rep.int(x = rownames(x = x), times = ncol(x = x)),
    cols = unlist(x = lapply(X = colnames(x = x), FUN = rep.int, times = nrow(x = x))),
    vals = unlist(x = x, use.names = FALSE)
  ))
}

# Check the existence of a package
#
# @param ... Package names
# @param error If true, throw an error if the package doesn't exist
#
# @return Invisibly returns boolean denoting if the package is installed
#
PackageCheck <- function(..., error = TRUE) {
  pkgs <- unlist(x = c(...), use.names = FALSE)
  package.installed <- vapply(
    X = pkgs,
    FUN = requireNamespace,
    FUN.VALUE = logical(length = 1L),
    quietly = TRUE
  )
  if (error && any(!package.installed)) {
    stop(
      "Cannot find ",
      paste(pkgs[!package.installed], collapse = ', '),
      "; please install"
    )
  }
  invisible(x = package.installed)
}

# Parenting parameters from one environment to the next
#
# This function allows one to modifiy a parameter in a parent environement
# The primary use of this is to ensure logging functions store correct parameters
# if they've been modified by a child function or method
#
# @param parent.find Regex pattern of name of parent function call to modify.
# For example, this can be the class name for a method that was dispatched previously
# @param ... Parameter names and values to parent; both name and value must be supplied
# in the standard \code{name = value} format; any number of name/value pairs can be specified
#
# @return No return, modifies parent environment directly
#
# @examples
# Parenting(parent.find = 'Seurat', features = features[features > 7])
#
Parenting <- function(parent.find = 'Seurat', ...) {
  calls <- as.character(x = sys.calls())
  calls <- lapply(
    X = strsplit(x = calls, split = '(', fixed = TRUE),
    FUN = '[',
    1
  )
  parent.index <- grep(pattern = parent.find, x = calls)
  if (length(x = parent.index) != 1) {
    stop("Cannot find a parent environment called ", parent.find)
  }
  to.parent <- list(...)
  if (length(x = to.parent) == 0) {
    stop("Nothing to parent")
  } else if (is.null(x = names(x = to.parent))) {
    stop("All input must be in a key = value pair")
  } else if (length(x = Filter(f = nchar, x = names(x = to.parent))) != length(x = to.parent)) {
    stop("All inputs must be named")
  }
  parent.environ <- sys.frame(which = parent.index)
  for (i in 1:length(x = to.parent)) {
    parent.environ[[names(x = to.parent)[i]]] <- to.parent[[i]]
  }
}

# Calculate the percentage of a vector above some threshold
#
# @param x Vector of values
# @param threshold Threshold to use when calculating percentage
#
# @return Returns the percentage of `x` values above the given threshold
#
PercentAbove <- function(x, threshold){
  return(length(x = x[x > threshold]) / length(x = x))
}

# Get the number of threads provided by the current plan
#
# @return The number of threads (workers) for the current future plan, or 1 if no workers detected
#
#' @importFrom future plan
#
PlanThreads <- function() {
  nthreads <- eval(expr = formals(fun = plan())$workers)
  return(nthreads %||% 1)
}

# Generate a random name
#
# Make a name from randomly sampled lowercase letters,
# pasted together with no spaces or other characters
#
# @param length How long should the name be
# @param ... Extra parameters passed to sample
#
# @return A character with nchar == length of randomly sampled letters
#
# @seealso \code{\link{sample}}
#
RandomName <- function(length = 5L, ...) {
  return(paste(sample(x = letters, size = length, ...), collapse = ''))
}

# Internal function for merging two matrices by rowname
#
# @param mat1 First matrix
# @param mat2 Second matrix
#
# @return A merged matrix
#
#' @importFrom methods as
#
RowMergeSparseMatrices <- function(mat1, mat2){
  if (inherits(x = mat1, what = "data.frame")) {
    mat1 <- as.matrix(x = mat1)
  }
  if (inherits(x = mat2, what = "data.frame")) {
    mat2 <- as.matrix(x = mat2)
  }
  mat1 <- as(object = mat1, Class = "RsparseMatrix")
  mat2 <- as(object = mat2, Class = "RsparseMatrix")
  mat1.names <- rownames(x = mat1)
  mat2.names <- rownames(x = mat2)
  all.names <- union(x = mat1.names, y = mat2.names)
  new.mat <- RowMergeMatrices(
    mat1 = mat1,
    mat2 = mat2,
    mat1_rownames = mat1.names,
    mat2_rownames = mat2.names,
    all_rownames = all.names
  )
  rownames(x = new.mat) <- make.unique(names = all.names)
  colnames(x = new.mat) <- make.unique(names = c(
    colnames(x = mat1),
    colnames(x = mat2)
  ))
  return(new.mat)
}

# Return what was passed
#
# @param x anything
#
# @return Returns x
#
Same <- function(x) {
  return(x)
}

# Try to convert x to numeric, if NA's introduced return x as is
#
ToNumeric <- function(x){
  # check for x:y range
  if (is.numeric(x = x)) {
    return(x)
  }
  if (length(x = unlist(x = strsplit(x = x, split = ":"))) == 2) {
    num <- unlist(x = strsplit(x = x, split = ":"))
    return(num[1]:num[2])
  }
  num <- suppressWarnings(expr = as.numeric(x = x))
  if (!is.na(x = num)) {
    return(num)
  }
  return(x)
}
