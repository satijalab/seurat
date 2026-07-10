#' Perform sctransform-based normalization (rewrite scaffold)
#'
#' Exact-copy scaffold of \code{\link{SCTransform}} for stepwise rewrites.
#'
#' @param object An object
#' @param ... Arguments passed to methods
#'
#' @export
SCTransform_rewrite <- function(object, ...) {
  UseMethod(generic = 'SCTransform_rewrite', object = object)
}

#' @rdname SCTransform_rewrite
#' @concept preprocessing
#' @export
#' @method SCTransform_rewrite default
#'
SCTransform_rewrite.default <- function(
  object,
  cell.attr,
  reference.SCT.model = NULL,
  do.correct.umi = TRUE,
  ncells = 5000,
  residual.features = NULL,
  variable.features.n = 3000,
  variable.features.rv.th = 1.3,
  vars.to.regress = NULL,
  latent.data = NULL,
  do.scale = FALSE,
  do.center = TRUE,
  clip.range = c(-sqrt(x = ncol(x = umi) / 30), sqrt(x = ncol(x = umi) / 30)),
  vst.flavor = 'v2',
  conserve.memory = FALSE,
  return.only.var.genes = TRUE,
  defer.residual.matrix = FALSE,
  seed.use = 1448145,
  verbose = TRUE,
  ...
) {
  if (!is.null(x = seed.use)) {
    set.seed(seed = seed.use)
  }
  vst.args <- list(...)
  object <- as.sparse(x = object)
  umi <- object
  # check for batch_var in meta data
  if ('batch_var' %in% names(x = vst.args)) {
    if (!(vst.args[['batch_var']] %in% colnames(x = cell.attr))) {
      stop('batch_var not found in seurat object meta data')
    }
  }
  # parameter checking when reference.SCT.model is set
  if (!is.null(x = reference.SCT.model) ) {
    if (inherits(x = reference.SCT.model, what = "SCTModel")) {
      reference.SCT.model <- SCTModel_to_vst(SCTModel = reference.SCT.model)
    }
    if (is.list(x = reference.SCT.model) & inherits(x = reference.SCT.model[[1]], what = "SCTModel")) {
      stop("reference.SCT.model must be one SCTModel rather than a list of SCTModel")
    }
    if ('latent_var' %in% names(x = vst.args)) {
      stop('custom latent variables are not supported when reference.SCT.model is given')
    }
    if (reference.SCT.model$model_str != 'y ~ log_umi') {
      stop('reference.SCT.model must be derived using default SCT regression formula, `y ~ log_umi`')
    }

  }
  # check for latent_var in meta data
  if ('latent_var' %in% names(x = vst.args)) {
    known.attr <- c('umi', 'gene', 'log_umi', 'log_gene', 'umi_per_gene', 'log_umi_per_gene')
    if (!all(vst.args[['latent_var']] %in% c(colnames(x = cell.attr), known.attr))) {
      stop('latent_var values are not from the set of cell attributes sctransform calculates by default and cannot be found in seurat object meta data')
    }
  }
  # check for vars.to.regress in meta data
  if (any(!vars.to.regress %in% colnames(x = cell.attr))) {
    stop('problem with second non-regularized linear regression; not all variables found in seurat object meta data; check vars.to.regress parameter')
  }
  if (any(c('cell_attr', 'verbosity', 'return_cell_attr', 'return_gene_attr', 'return_corrected_umi') %in% names(x = vst.args))) {
    warning(
      'the following arguments will be ignored because they are set within this function:',
      paste(
        c(
          'cell_attr',
          'verbosity',
          'return_cell_attr',
          'return_gene_attr',
          'return_corrected_umi'
        ),
        collapse = ', '
      ),
      call. = FALSE,
      immediate. = TRUE
    )
  }

  if (!is.null(x = vst.flavor) && !vst.flavor %in% c("v1", "v2")){
    stop("vst.flavor can be 'v1' or 'v2'. Default is 'v2'")
  }
  if (!is.null(x = vst.flavor) && vst.flavor == "v1"){
    vst.flavor <- NULL
  }

  vst.args[['vst.flavor']] <- vst.flavor
  vst.args[['umi']] <- umi
  vst.args[['cell_attr']] <- cell.attr
  vst.args[['verbosity']] <- as.numeric(x = verbose) * 1
  vst.args[['return_cell_attr']] <- TRUE
  vst.args[['return_gene_attr']] <- TRUE
  vst.args[['return_corrected_umi']] <- do.correct.umi
  vst.args[['n_cells']] <- min(ncells, ncol(x = umi))
  residual.type <- vst.args[['residual_type']] %||% 'pearson'
  res.clip.range <- vst.args[['res_clip_range']] %||% c(-sqrt(x = ncol(x = umi)), sqrt(x = ncol(x = umi)))
 # set sct normalization method
  if (!is.null( reference.SCT.model)) {
    sct.method <- "reference.model"
  } else if (!is.null(x = residual.features)) {
    sct.method <- "residual.features"
  } else if (conserve.memory) {
    sct.method <- "conserve.memory"
  } else {
    sct.method <- "default"
  }
  # set vst model
  vst.out <- switch(
    EXPR = sct.method,
    'default' = {
      vst.args[['return_corrected_umi']] <- FALSE
      vst.args[['residual_type']] <- 'none'
      vst.out <- do.call(what = 'vst', args = vst.args)
      vst.out
    },
    'reference.model' = {
      if (verbose) {
        message("Using reference SCTModel to calculate pearson residuals")
      }
      do.center <- FALSE
      do.correct.umi <- FALSE
      vst.out <- reference.SCT.model
      clip.range <- vst.out$arguments$sct.clip.range
      cell_attr <-  data.frame(log_umi = log10(x = colSums(umi)))
      rownames(cell_attr) <- colnames(x = umi)
      vst.out$cell_attr <- cell_attr

      all.features  <- intersect(
        x =  rownames(x = vst.out$gene_attr),
        y = rownames(x = umi)
      )
      vst.out$gene_attr <- vst.out$gene_attr[all.features ,]
      vst.out$model_pars_fit <- vst.out$model_pars_fit[all.features,]
      vst.out
    },
    'residual.features' = {
      if (verbose) {
        message("Computing residuals for the ", length(x = residual.features), " specified features")
      }
      return.only.var.genes <- TRUE
      do.correct.umi <- FALSE
      vst.args[['return_corrected_umi']] <- FALSE
      vst.args[['residual_type']] <- 'none'
      vst.out <- do.call(what = 'vst', args = vst.args)
      vst.out$gene_attr$residual_variance <- NA_real_
      vst.out
    })

  # get residuals
  vst.out <- switch(
    EXPR = sct.method,
     # Default SCTransform behavior - compute Pearson residuals for all genes
    # now performed using optimized C++ workflow
    'default' = {
      # setup everything for the optimized C++ workflow
      model.pars <- vst.out$model_pars_fit
      genes <- rownames(x = model.pars)
      if (!identical(x = genes, y = rownames(x = umi))) {
        umi <- umi[genes, , drop = FALSE]
      }
      min.variance <- vst.out$arguments$min_variance
      min.var <- if (identical(x = min.variance, y = "umi_median")) {
        (median(umi@x) / 5) ^ 2
      } else {
        min.variance
      }
      # Persist the resolved numeric min_var in the model. Downstream residual
      # recomputation (FetchResiduals / GetResidual, old or _rewrite worker) reads
      # arguments$min_variance and only recomputes (median(nonzeros)/5)^2 when it is
      # the string "umi_median". That recompute is order/subset dependent (e.g. the
      # old worker uses only the first chunk_size cells), so it can diverge from the
      # value used here for scale.data. Storing the resolved value makes later
      # residuals deterministic and consistent with scale.data.
      vst.out$arguments$min_variance <- min.var
      # should be set already by the vst call but just fixing in case its null
      res.clip.range <- vst.out$arguments$res_clip_range %||%
        c(-sqrt(x = ncol(x = umi)), sqrt(x = ncol(x = umi)))

      # Compute residual statistics and corrected UMI counts
      # Note: does not compute residual matrix yet (saves a lot of memory)
      # Just computes the residual variance for each gene and (if asked for) corrected UMI counts
      # One of key optims was to compute corrected counts if needed at this stage, not later
      stats <- SCTResidualStatsAndCorrected_optimized(
        x = umi@x,
        i = umi@i,
        p = umi@p,
        rows = nrow(x = umi),
        cols = ncol(x = umi),
        theta = model.pars[, "theta"],
        intercept = model.pars[, "(Intercept)"],
        slope = model.pars[, "log_umi"],
        log_umi = vst.out$cell_attr[colnames(x = umi), "log_umi"],
        target_log_umi = median(vst.out$cell_attr[, "log_umi"]),
        min_var = min.var,
        residual_clip_min = min(res.clip.range),
        residual_clip_max = max(res.clip.range),
        n_threads = getOption(x = "Seurat.nthreads", default = 1L),
        compute_corrected = do.correct.umi
      )
      vst.out$gene_attr[genes, "residual_mean"] <- stats$residual_mean
      vst.out$gene_attr[genes, "residual_variance"] <- stats$residual_variance
    
      # Determine variable features
      feature.variance <- vst.out$gene_attr[, "residual_variance"]
      names(x = feature.variance) <- rownames(x = vst.out$gene_attr)
      feature.variance <- sort(x = feature.variance, decreasing = TRUE)
      feature.idx <- if (is.null(x = variable.features.n)) {
        feature.variance >= variable.features.rv.th
      } else {
        seq_len(length.out = min(variable.features.n, length(x = feature.variance)))
      }
      top.features <- names(x = feature.variance)[feature.idx]

      # Store corrected UMI counts if requested, if not just restore original counts matrix
      if (do.correct.umi) {
        vst.out$umi_corrected <- stats$corrected
        dimnames(x = vst.out$umi_corrected) <- dimnames(x = umi)
      } else {
        vst.out$umi_corrected <- umi
      }

      # Compute matrix of Pearson residuals for features to be included in scale.data
      scale.data.features <- if (return.only.var.genes) {
        top.features
      } else {
        genes
      }
      if (isTRUE(x = defer.residual.matrix)) {
        vst.out$y <- matrix(
          data = numeric(length = 0L),
          nrow = 0L,
          ncol = ncol(x = umi),
          dimnames = list(character(length = 0L), colnames(x = umi))
        )
      } else {
        vst.out$y <- SCTPearsonResidualMatrix_optimized(
          x = umi@x,
          i = umi@i,
          p = umi@p,
          rows = nrow(x = umi),
          cols = ncol(x = umi),
          theta = model.pars[, "theta"],
          intercept = model.pars[, "(Intercept)"],
          slope = model.pars[, "log_umi"],
          log_umi = vst.out$cell_attr[colnames(x = umi), "log_umi"],
          feature_index = as.integer(x = match(x = scale.data.features, table = genes) - 1L),
          min_var = min.var,
          clip_min = min(clip.range),
          clip_max = max(clip.range),
          do_center = do.center,
          n_threads = getOption(x = "Seurat.nthreads", default = 1L)
        )
        dimnames(x = vst.out$y) <- list(scale.data.features, colnames(x = umi))
      }
      
      vst.out
    },
    'reference.model' = {
      feature.variance <- vst.out$gene_attr[, "residual_variance"]
      names(x = feature.variance) <- rownames(x = vst.out$gene_attr)

      feature.variance <- sort(x = feature.variance, decreasing = TRUE)

      feature.idx <- if (is.null(x = variable.features.n)) {
        feature.variance >= variable.features.rv.th
      } else {
        seq_len(length.out = min(variable.features.n, length(x = feature.variance)))
      }
      top.features <- names(x = feature.variance)[feature.idx]
      if (is.null(x = residual.features)) {
        residual.features <- top.features
      }

      residual.features <- Reduce(
        f = intersect,
        x = list(residual.features, rownames(x = umi), rownames(x = vst.out$model_pars_fit))
      )
      residual.feature.mat <- get_residuals(
        vst_out = vst.out,
        umi = umi[residual.features, , drop = FALSE],
        verbosity = as.numeric(x = verbose)*2
      )
      vst.out$gene_attr <- vst.out$gene_attr[residual.features ,]
      ref.residuals.mean <- vst.out$gene_attr[,"residual_mean"]
      vst.out$y <- sweep(
        x = residual.feature.mat,
        MARGIN = 1,
        STATS = ref.residuals.mean,
        FUN = "-"
      )
      vst.out
    },
    'residual.features' = {
      residual.features <- intersect(
        x = residual.features,
        y = rownames(x = vst.out$gene_attr)
      )
      residual.feature.mat <- get_residuals(
        vst_out = vst.out,
        umi = umi[residual.features, , drop = FALSE],
        verbosity = as.numeric(x = verbose)*2
      )
      vst.out$y <- residual.feature.mat
      vst.out$gene_attr$residual_mean <- NA_real_
      vst.out$gene_attr$residual_variance <- NA_real_
      vst.out$gene_attr[residual.features, "residual_mean"] <- rowMeans2(x = vst.out$y)
      vst.out$gene_attr[residual.features, "residual_variance"] <- RowVar(x = vst.out$y)
      vst.out
    }
   )
  # default method already clips residuals
  if (!identical(x = sct.method, y = "default")) {
    scale.data <- vst.out$y
    scale.data[scale.data < clip.range[1]] <- clip.range[1]
    scale.data[scale.data > clip.range[2]] <- clip.range[2]
    vst.out$y <- scale.data
  }
  
  # User may (not common) want to regress out additional variables after SCTransform
  # Note that centering is already handled by the optimized residual matrix C++
  if (!is.null(x = vars.to.regress) || isTRUE(x = do.scale)) {
    vst.out$y <- ScaleData(
      vst.out$y,
      features = NULL,
      vars.to.regress = vars.to.regress,
      latent.data = latent.data,
      model.use = 'linear',
      use.umi = FALSE,
      do.scale = do.scale,
      do.center = do.center,
      scale.max = Inf,
      block.size = 750,
      min.cells.to.block = 3000,
      verbose = verbose
    )
  }

  vst.out$variable_features <- residual.features %||% top.features
  if (!do.correct.umi) {
    vst.out$umi_corrected <- umi
  }
  min_var <- vst.out$arguments$min_variance
  return(vst.out)
}

#' @rdname SCTransform_rewrite
#' @concept preprocessing
#' @export
#' @method SCTransform_rewrite Assay
#'
SCTransform_rewrite.Assay <- function(
    object,
    cell.attr,
    reference.SCT.model = NULL,
    do.correct.umi = TRUE,
    ncells = 5000,
    residual.features = NULL,
    variable.features.n = 3000,
    variable.features.rv.th = 1.3,
    vars.to.regress = NULL,
    latent.data = NULL,
    do.scale = FALSE,
    do.center = TRUE,
    clip.range = c(-sqrt(x = ncol(x = object) / 30), sqrt(x = ncol(x = object) / 30)),
    vst.flavor = 'v2',
    conserve.memory = FALSE,
    return.only.var.genes = TRUE,
    seed.use = 1448145,
    verbose = TRUE,
    ...
) {
  if (!is.null(x = seed.use)) {
    set.seed(seed = seed.use)
  }
  if (!is.null(reference.SCT.model)){
    do.correct.umi <- FALSE
    do.center <- FALSE
    rlang::warn(
      "A reference SCT model was provided, therefore counts are not corrected (regardless of do.correct.umi)",
      .frequency = "once",
      .frequency_id = "SCTransform-reference-SCTmodel-correct-counts"
    )
  }

  umi <- GetAssayData(object = object, layer = 'counts')
  vst.out <- SCTransform_rewrite(object = umi,
                         cell.attr = cell.attr,
                         reference.SCT.model = reference.SCT.model,
                         do.correct.umi = do.correct.umi,
                         ncells = ncells,
                         residual.features = residual.features,
                         variable.features.n = variable.features.n,
                         variable.features.rv.th = variable.features.rv.th,
                         vars.to.regress = vars.to.regress,
                         latent.data = latent.data,
                         do.scale = do.scale,
                         do.center = do.center,
                         clip.range = clip.range,
                         vst.flavor = vst.flavor,
                         conserve.memory = conserve.memory,
                         return.only.var.genes = return.only.var.genes,
                         seed.use = seed.use,
                         verbose = verbose,
                         ...)
  residual.type <- vst.out[['residual_type']] %||% 'pearson'
  sct.method <- vst.out[["sct.method"]]
  # create output assay and put (corrected) umi counts in count slot
  if (do.correct.umi & residual.type == 'pearson') {
    if (verbose) {
      message('Place corrected count matrix in counts slot')
    }
    assay.out <- CreateAssayObject(counts = vst.out$umi_corrected)
    vst.out$umi_corrected <- NULL
  } else {
    # TODO: restore once check.matrix is in SeuratObject
    # assay.out <- CreateAssayObject(counts = umi, check.matrix = FALSE)
    assay.out <- CreateAssayObject(counts = umi)
  }
  # set the variable genes
  VariableFeatures(object = assay.out) <- vst.out$variable_features
  # put log1p transformed counts in data
  assay.out <- SetAssayData(
    object = assay.out,
    layer = 'data',
    new.data = log1p(x = GetAssayData(object = assay.out, layer = 'counts'))
  )
  scale.data <- vst.out$y
  assay.out <- SetAssayData(
    object = assay.out,
    layer = 'scale.data',
    new.data = scale.data
  )
  vst.out$y <- NULL
  # save clip.range into vst model
  vst.out$arguments$sct.clip.range <- clip.range
  vst.out$arguments$sct.method <- sct.method
  Misc(object = assay.out, slot = 'vst.out') <- vst.out
  assay.out <- as(object = assay.out, Class = "SCTAssay")
  return(assay.out)
}

#' @param assay Name of assay to pull the count data from; default is 'RNA'
#' @param new.assay.name Name for the new assay containing the normalized data; default is 'SCT'
#'
#' @rdname SCTransform_rewrite
#' @concept preprocessing
#' @export
#' @method SCTransform_rewrite Seurat
#'
SCTransform_rewrite.Seurat <- function(
    object,
    assay = "RNA",
    new.assay.name = 'SCT',
    reference.SCT.model = NULL,
    do.correct.umi = TRUE,
    ncells = 5000,
    residual.features = NULL,
    variable.features.n = 3000,
    variable.features.rv.th = 1.3,
    vars.to.regress = NULL,
    do.scale = FALSE,
    do.center = TRUE,
    clip.range = c(-sqrt(x = ncol(x = object[[assay]]) / 30), sqrt(x = ncol(x = object[[assay]]) / 30)),
    vst.flavor = "v2",
    conserve.memory = FALSE,
    return.only.var.genes = TRUE,
    seed.use = 1448145,
    verbose = TRUE,
    ...
) {
  if (!is.null(x = seed.use)) {
    set.seed(seed = seed.use)
  }
  if (any(vars.to.regress %in% colnames(x = object[[]]))) {
    vars.to.regress.subset <- vars.to.regress[vars.to.regress %in% colnames(x = object[[]])]
    latent.data <- object[[vars.to.regress.subset]]
  } else {
    latent.data <- NULL
  }
  assay <- assay %||% DefaultAssay(object = object)
  if (assay == "SCT") {
    # if re-running SCTransform, use the RNA assay
    assay <- "RNA"
    warning("Running SCTransform on the RNA assay while default assay is SCT.")
  }

  if (verbose){
    message("Running SCTransform on assay: ", assay)
  }
  cell.attr <- slot(object = object, name = 'meta.data')[colnames(object[[assay]]),]
  assay.data <- SCTransform_rewrite(object = object[[assay]],
                            cell.attr = cell.attr,
                            reference.SCT.model = reference.SCT.model,
                            do.correct.umi = do.correct.umi,
                            ncells = ncells,
                            residual.features = residual.features,
                            variable.features.n = variable.features.n,
                            variable.features.rv.th = variable.features.rv.th,
                            vars.to.regress = vars.to.regress,
                            latent.data = latent.data,
                            do.scale = do.scale,
                            do.center = do.center,
                            clip.range = clip.range,
                            vst.flavor = vst.flavor,
                            conserve.memory = conserve.memory,
                            return.only.var.genes = return.only.var.genes,
                            seed.use = seed.use,
                            verbose = verbose,
                            ...)
  assay.data <- SCTAssay(assay.data, assay.orig = assay)
  
  # Extract all SCT models stored in assay
  sct_models <- slot(object = assay.data, name = "SCTModel.list")
  
  # Update umi.assay field for every SCT model 
  slot(object = assay.data, name = "SCTModel.list") <- lapply(sct_models, function(model) {
    slot(model, name = "umi.assay") <- assay
    model
  })

  object[[new.assay.name]] <- assay.data

  if (verbose) {
    message(paste("Set default assay to", new.assay.name))
  }
  DefaultAssay(object = object) <- new.assay.name
  object <- LogSeuratCommand(object = object)
  return(object)
}

#' @concept preprocessing
#' @export
#'
FetchResiduals_rewrite <- function(object, ...) {
  UseMethod(generic = "FetchResiduals_rewrite", object = object)
}

#' @concept preprocessing
#' @export
#' @method FetchResiduals_rewrite SCTAssay
#'
FetchResiduals_rewrite.SCTAssay <- function(
  object,
  umi.object,
  features,
  layer = "counts",
  clip.range = NULL,
  reference.SCT.model = NULL,
  replace.value = FALSE,
  na.rm = TRUE,
  verbose = TRUE,
  ...
) {
  sct.models <- levels(x = object)
  if (length(x = sct.models) == 0) {
    warning("SCT model not present in assay", call. = FALSE, immediate. = TRUE)
    return(LayerData(object, layer = "scale.data"))
  }

  model.features <- lapply(
    X = sct.models,
    FUN = function(model) {
      rownames(x = SCTResults(object = object, slot = "feature.attributes", model = model))
    }
  )
  names(x = model.features) <- sct.models

  possible.features <- Reduce(f = union, x = model.features)
  bad.features <- setdiff(x = features, y = possible.features)
  if (length(x = bad.features) > 0) {
    warning(
      "The following requested features are not present in any models: ",
      paste(bad.features, collapse = ", "),
      call. = FALSE
    )
    features <- intersect(x = features, y = possible.features)
  }

  features.orig <- features
  if (isTRUE(x = na.rm)) {
    common.features <- Reduce(f = intersect, x = model.features)
    features <- intersect(x = features.orig, y = common.features)
  }
  if (length(x = features) < 1) {
    warning(
      "The following requested features are not present in all the models: ",
      paste(features.orig, collapse = ", "),
      call. = FALSE
    )
    return(LayerData(object, layer = "scale.data"))
  }

  if (!is.null(x = reference.SCT.model)) {
    if (inherits(x = reference.SCT.model, what = "SCTModel")) {
      reference.SCT.model <- SCTModel_to_vst(SCTModel = reference.SCT.model)
    }
    if (is.list(x = reference.SCT.model) && inherits(x = reference.SCT.model[[1]], what = "SCTModel")) {
      stop("reference.SCT.model must be one SCTModel rather than a list of SCTModel")
    }
    if (reference.SCT.model$model_str != "y ~ log_umi") {
      stop("reference.SCT.model must be derived using default SCT regression formula, `y ~ log_umi`")
    }
  }

  layers <- Layers(object = umi.object, search = layer)
  if (length(x = layers) != length(x = sct.models)) {
    stop("The number of UMI layers must match the number of SCT models")
  }

  residuals.list <- vector(mode = "list", length = length(x = layers))
  for (i in seq_along(along.with = layers)) {
    residuals.list[[i]] <- FetchResidualSCTModel_rewrite(
      object = object,
      umi.object = umi.object,
      layer = layers[[i]],
      layer.cells = Cells(x = umi.object, layer = layers[[i]]),
      SCTModel = sct.models[[i]],
      reference.SCT.model = reference.SCT.model,
      new_features = features,
      replace.value = replace.value,
      clip.range = clip.range,
      verbose = verbose
    )
  }

  residuals <- if (length(x = residuals.list) == 1L) {
    residuals.list[[1L]]
  } else {
    do.call(what = cbind, args = residuals.list)
  }

  if (isTRUE(x = na.rm)) {
    has.missing.residuals <- any(vapply(
      X = residuals.list,
      FUN = function(x) {
        isTRUE(x = attr(x = x, which = "has_missing_residuals"))
      },
      FUN.VALUE = logical(length = 1L)
    ))
    if (has.missing.residuals) {
      keep.features <- !apply(X = residuals, MARGIN = 1, FUN = anyNA)
      residuals <- residuals[keep.features, , drop = FALSE]
      features <- intersect(x = features, y = rownames(x = residuals))
    }
  }

  if (identical(x = rownames(x = residuals), y = features)) {
    return(residuals)
  }
  return(residuals[features, , drop = FALSE])
}

#' @concept preprocessing
#' @export
#'
FetchResidualSCTModel_rewrite <- function(
  object,
  umi.object,
  layer = "counts",
  chunk_size = 2000,
  layer.cells = NULL,
  SCTModel = NULL,
  reference.SCT.model = NULL,
  new_features = NULL,
  clip.range = NULL,
  replace.value = FALSE,
  verbose = FALSE
) {
  layer.cells <- layer.cells %||% Cells(x = umi.object, layer = layer)
  model.features <- Features(x = object)
  model.cells <- character()
  sct.method <- "reference"

  if (is.null(x = reference.SCT.model)) {
    clip.range <- clip.range %||% SCTResults(object = object, slot = "clips", model = SCTModel)$sct
    model.features <- rownames(x = SCTResults(object = object, slot = "feature.attributes", model = SCTModel))
    model.cells <- Cells(x = slot(object = object, name = "SCTModel.list")[[SCTModel]])
    sct.method <- SCTResults(object = object, slot = "arguments", model = SCTModel)$sct.method %||% "default"
  }

  existing.scale.data <- NULL
  reusable.features <- character()
  if (is.null(x = reference.SCT.model)) {
    existing.scale.data <- suppressWarnings(GetAssayData(object = object, layer = "scale.data"))
    existing.cells <- intersect(x = colnames(x = existing.scale.data), y = layer.cells)
    has.all.layer.cells <- length(x = setdiff(x = layer.cells, y = existing.cells)) == 0
    if (has.all.layer.cells && nrow(x = existing.scale.data) > 0) {
      existing.layer.data <- existing.scale.data[, layer.cells, drop = FALSE]
      reusable.features <- rownames(x = existing.layer.data)[
        !apply(X = existing.layer.data, MARGIN = 1, FUN = anyNA)
      ]
    }
  }

  features.to.compute <- if (isTRUE(x = replace.value)) {
    new_features
  } else {
    setdiff(x = new_features, y = reusable.features)
  }

  if (identical(x = sct.method, y = "reference.model")) {
    if (isTRUE(x = verbose)) {
      message("sct.model ", SCTModel, " is from reference, so no residuals will be recalculated")
    }
    features.to.compute <- character()
  }

  old.features <- intersect(x = new_features, y = reusable.features)

  if (length(x = features.to.compute) == 0) {
    result <- matrix(
      data = NA_real_,
      nrow = length(x = new_features),
      ncol = length(x = layer.cells),
      dimnames = list(new_features, layer.cells)
    )
    if (length(x = old.features) > 0) {
      result[old.features, layer.cells] <- existing.scale.data[old.features, layer.cells, drop = FALSE]
    }
    attr(x = result, which = "has_missing_residuals") <- anyNA(x = result)
    return(result)
  }

  missing.features <- setdiff(x = features.to.compute, y = model.features)
  compute.features <- intersect(x = features.to.compute, y = model.features)
  if (length(x = missing.features) > 0) {
    warning(
      "In the SCTModel ", SCTModel, ", the following ", length(x = missing.features),
      " features do not exist in the counts slot: ", paste(missing.features, collapse = ", ")
    )
  }
  if (length(x = compute.features) == 0) {
    result <- matrix(
      data = NA_real_,
      nrow = length(x = new_features),
      ncol = length(x = layer.cells),
      dimnames = list(new_features, layer.cells)
    )
    if (length(x = old.features) > 0) {
      result[old.features, layer.cells] <- existing.scale.data[old.features, layer.cells, drop = FALSE]
    }
    attr(x = result, which = "has_missing_residuals") <- anyNA(x = result)
    return(result)
  }

  if (is.null(x = reference.SCT.model)) {
    vst.out <- SCTModel_to_vst(SCTModel = slot(object = object, name = "SCTModel.list")[[SCTModel]])
    clip.range <- clip.range %||%
      vst.out$arguments$sct.clip.range %||%
      vst.out$arguments$clip.range %||%
      SCTResults(object = object, slot = "clips", model = SCTModel)$sct
  } else {
    vst.out <- SCTModel_to_vst(SCTModel = reference.SCT.model)
    clip.range <- clip.range %||%
      vst.out$arguments$sct.clip.range %||%
      vst.out$arguments$clip.range
    vst.out$cell_attr <- NULL
    vst.features <- intersect(x = rownames(x = vst.out$gene_attr), y = compute.features)
    vst.out$gene_attr <- vst.out$gene_attr[vst.features, , drop = FALSE]
    vst.out$model_pars_fit <- vst.out$model_pars_fit[vst.features, , drop = FALSE]
  }
  clip.range <- clip.range %||% c(
    -sqrt(x = ncol(x = umi.object) / 30),
    sqrt(x = ncol(x = umi.object) / 30)
  )

  clip.max <- max(clip.range)
  clip.min <- min(clip.range)
  counts <- LayerData(umi.object, layer = layer, cells = layer.cells)

  if (is.null(x = reference.SCT.model) && identical(x = sct.method, y = "default")) {
    counts <- as.sparse(x = counts)
    model.pars <- vst.out$model_pars_fit
    genes <- rownames(x = model.pars)
    if (!identical(x = genes, y = rownames(x = counts))) {
      counts <- counts[genes, , drop = FALSE]
    }
    min.variance <- vst.out$arguments$min_variance
    min.var <- if (identical(x = min.variance, y = "umi_median")) {
      (median(counts@x) / 5) ^ 2
    } else {
      min.variance
    }
    new.residuals <- SCTPearsonResidualMatrix_optimized(
      x = counts@x,
      i = counts@i,
      p = counts@p,
      rows = nrow(x = counts),
      cols = ncol(x = counts),
      theta = model.pars[, "theta"],
      intercept = model.pars[, "(Intercept)"],
      slope = model.pars[, "log_umi"],
      log_umi = vst.out$cell_attr[colnames(x = counts), "log_umi"],
      feature_index = as.integer(x = match(x = compute.features, table = genes) - 1L),
      min_var = min.var,
      clip_min = clip.min,
      clip_max = clip.max,
      do_center = TRUE,
      n_threads = getOption(x = "Seurat.nthreads", default = 1L)
    )
    dimnames(x = new.residuals) <- list(compute.features, colnames(x = counts))
    if (
      length(x = old.features) == 0 &&
      length(x = missing.features) == 0 &&
      identical(x = compute.features, y = new_features)
    ) {
      attr(x = new.residuals, which = "has_missing_residuals") <- FALSE
      return(new.residuals)
    }
    result <- matrix(
      data = NA_real_,
      nrow = length(x = new_features),
      ncol = length(x = layer.cells),
      dimnames = list(new_features, layer.cells)
    )
    if (length(x = old.features) > 0) {
      result[old.features, layer.cells] <- existing.scale.data[old.features, layer.cells, drop = FALSE]
    }
    result[rownames(x = new.residuals), colnames(x = new.residuals)] <- new.residuals
    attr(x = result, which = "has_missing_residuals") <- anyNA(x = result)
    return(result)
  }

  result <- matrix(
    data = NA_real_,
    nrow = length(x = new_features),
    ncol = length(x = layer.cells),
    dimnames = list(new_features, layer.cells)
  )
  if (length(x = old.features) > 0) {
    result[old.features, layer.cells] <- existing.scale.data[old.features, layer.cells, drop = FALSE]
  }

  cells.vector <- seq_along(along.with = layer.cells)
  cells.grid <- split(x = cells.vector, f = ceiling(x = cells.vector / chunk_size))
  residuals.list <- vector(mode = "list", length = length(x = cells.grid))

  for (i in seq_along(along.with = cells.grid)) {
    vp <- cells.grid[[i]]
    umi.all <- as.sparse(x = counts[, vp, drop = FALSE])
    umi <- umi.all[compute.features, , drop = FALSE]

    if (i == 1L) {
      nz.median <- median(umi.all@x)
      min.var.custom <- (nz.median / 5)^2
    }

    cell.attr <- data.frame(
      umi = colSums(x = umi.all),
      log_umi = log10(x = colSums(x = umi.all))
    )
    rownames(x = cell.attr) <- colnames(x = umi.all)

    vst.out.tmp <- vst.out
    if (sct.method %in% c("reference.model", "reference")) {
      vst.out.tmp$cell_attr <- cell.attr[colnames(x = umi.all), , drop = FALSE]
    } else {
      cell.attr.existing <- vst.out.tmp$cell_attr
      cells.missing <- setdiff(x = rownames(x = cell.attr), y = rownames(x = cell.attr.existing))
      if (length(x = cells.missing) > 0) {
        cell.attr.missing <- cell.attr[cells.missing, , drop = FALSE]
        missing.cols <- setdiff(x = colnames(x = cell.attr.existing), y = colnames(x = cell.attr.missing))
        if (length(x = missing.cols) > 0) {
          cell.attr.missing[, missing.cols] <- NA
        }
        cell.attr.existing <- rbind(cell.attr.existing, cell.attr.missing)
      }
      vst.out.tmp$cell_attr <- cell.attr.existing[colnames(x = umi), , drop = FALSE]
    }

    min.var <- if (identical(x = vst.out.tmp$arguments$min_variance, y = "umi_median")) {
      min.var.custom
    } else {
      vst.out.tmp$arguments$min_variance
    }

    residuals.list[[i]] <- as.matrix(x = get_residuals(
      vst_out = vst.out.tmp,
      umi = umi,
      residual_type = "pearson",
      min_variance = min.var,
      res_clip_range = c(clip.min, clip.max),
      verbosity = as.numeric(x = verbose) * 2
    ))
  }

  new.residuals <- do.call(what = cbind, args = residuals.list)
  if (is.null(x = reference.SCT.model)) {
    new.residuals <- new.residuals - rowMeans(x = new.residuals)
  } else {
    if (isTRUE(x = verbose)) {
      message("Using residual mean from reference for centering")
    }
    ref.vst.out <- SCTModel_to_vst(SCTModel = reference.SCT.model)
    ref.residuals.mean <- ref.vst.out$gene_attr[rownames(x = new.residuals), "residual_mean"]
    new.residuals <- sweep(
      x = new.residuals,
      MARGIN = 1,
      STATS = ref.residuals.mean,
      FUN = "-"
    )
  }

  result[rownames(x = new.residuals), colnames(x = new.residuals)] <- new.residuals
  attr(x = result, which = "has_missing_residuals") <- anyNA(x = result)
  return(result)
}


#' @rdname SCTransform_rewrite
#' @concept preprocessing
#' @export
#' @method SCTransform_rewrite StdAssay
#'
SCTransform_rewrite.StdAssay <- function(
  object,
  layer = 'counts',
  cell.attr = NULL,
  reference.SCT.model = NULL,
  do.correct.umi = TRUE,
  ncells = 5000,
  residual.features = NULL,
  variable.features.n = 3000,
  variable.features.rv.th = 1.3,
  vars.to.regress = NULL,
  latent.data = NULL,
  do.scale = FALSE,
  do.center = TRUE,
  clip.range = c(-sqrt(x = ncol(x = object) / 30), sqrt(x = ncol(x = object) / 30)),
  vst.flavor = 'v2',
  conserve.memory = FALSE,
  return.only.var.genes = TRUE,
  seed.use = 1448145,
  verbose = TRUE,
  ...
) {
  

  # Extract counts layers
  layer_names <- Layers(object, search = layer)
  input_list <- lapply(
    layer_names,
    function(layer_name) {
      layer_counts <- LayerData(object, layer = layer_name)
      return(layer_counts)
    }
  )
  names(x = input_list) <- layer_names

  # Apply SCTransform to each set of counts in `input_list`.
  output_list <- lapply(
    names(x = input_list),
    function(layer_name) {
      input <- input_list[[layer_name]]
      layer_metadata <- cell.attr[colnames(x = input), , drop = FALSE]
      result <- SCTransform_rewrite(
        input,
        cell.attr = layer_metadata,
        reference.SCT.model = reference.SCT.model,
        do.correct.umi = do.correct.umi,
        ncells = ncells,
        residual.features = residual.features,
        variable.features.n = variable.features.n,
        variable.features.rv.th = variable.features.rv.th,
        vars.to.regress = vars.to.regress,
        latent.data = latent.data,
        do.scale = do.scale,
        do.center = do.center,
        clip.range = clip.range,
        vst.flavor = vst.flavor,
        conserve.memory = conserve.memory,
        return.only.var.genes = return.only.var.genes,
        defer.residual.matrix = TRUE,
        seed.use = seed.use,
        verbose = verbose,
        ...
      )
    }
  )
  names(x = output_list) <- names(x = input_list)


  # Merge counts assays into one, or take the single result.
  counts_list <- if (do.correct.umi) {
    lapply(output_list, function(vst.out) {
      vst.out$umi_corrected
    })
  } else {
    input_list
  }
  if (length(x = counts_list) == 1L) {
    counts <- counts_list[[1]]
  } else {
    same_count_features <- all(vapply(
      X = counts_list[-1],
      FUN = function(mat) {
        identical(x = rownames(x = mat), y = rownames(x = counts_list[[1]]))
      },
      FUN.VALUE = logical(length = 1L)
    ))

    counts <- if (same_count_features) {
      do.call(what = cbind, args = counts_list)
    } else {
      RowMergeSparseMatrices(
        mat1 = counts_list[[1]],
        mat2 = counts_list[-1]
      )
    }
  }
  
  # Determine which features to include in the output's scale.data slot.
  if (return.only.var.genes) {
    var_features_union <- Reduce(
      f = union,
      x = lapply(output_list, function(vst.out) {
        vst.out$variable_features
      })
    )

    all_features_intersect <- Reduce(
      f = intersect,
      x = lapply(counts_list, function(counts) {
        rownames(x = counts)
      })
    )

    scale_data_features <- intersect(
      x = all_features_intersect,
      y = var_features_union
    )
  } else {
    scale_data_features <- Reduce(
      f = union,
      x = lapply(counts_list, function(counts) {
        rownames(x = counts)
      })
      )
  }

  # Create output assay and put log1p transformed counts in data slot
  assay_out <- CreateAssayObject(counts = counts)
  LayerData(object = assay_out, layer = "data") <- log1p(x = counts)
  model.list <- lapply(
    X = output_list,
    FUN = function(vst.out) {
      PrepVSTResults(
        vst.res = vst.out,
        cell.names = rownames(x = vst.out$cell_attr)
      )
    }
  )
  names(x = model.list) <- paste0("model", seq_along(along.with = model.list))
  assay_out <- as(object = assay_out, Class = "SCTAssay")
  slot(object = assay_out, name = "SCTModel.list") <- model.list

  # pre-fill scale.data if pearson residuals are already computed for all layers.
  # In the rewrite2 branch, residual matrices are deferred per layer so the final
  # multi-layer residual matrix is computed exactly once below.
  prefill.matrices <- lapply(output_list, function(vst.out) {
    vst.out$y
  })
  prefill.features <- character()
  if (all(vapply(X = prefill.matrices, FUN = nrow, FUN.VALUE = integer(length = 1L)) > 0L)) {
    prefill.features <- Reduce(
      f = intersect,
      x = lapply(prefill.matrices, rownames)
    )
    scale.data.prefill <- do.call(
      what = cbind,
      args = lapply(prefill.matrices, function(y) {
        y[prefill.features, , drop = FALSE]
      })
    )
    LayerData(assay_out, layer = "scale.data") <- scale.data.prefill
  }

  # In reference mode the final FetchResiduals_rewrite() below is intentionally
  # called WITHOUT reference.SCT.model and instead reuses the per-layer reference
  # residuals prefilled above. This is correct only while every scale.data feature
  # is covered by the prefill: any feature not prefilled would be recomputed
  # without reference centering (query-centered) and be silently wrong. That
  # invariant provably holds today (scale_data_features is a subset of the shared
  # reference model's features present in all layers), so guard it here to fail
  # loudly if a future change ever breaks it.
  if (!is.null(x = reference.SCT.model)) {
    missing.prefill <- setdiff(x = scale_data_features, y = prefill.features)
    if (length(x = missing.prefill) > 0) {
      stop(
        "SCTransform_rewrite (reference model, multi-layer): ",
        length(x = missing.prefill),
        " scale.data feature(s) are not covered by the per-layer reference ",
        "residuals and would be recomputed without reference centering: ",
        paste(utils::head(x = missing.prefill, n = 10L), collapse = ", "),
        if (length(x = missing.prefill) > 10L) ", ..." else "",
        call. = FALSE
      )
    }
  }

  residuals <- suppressWarnings(
    FetchResiduals_rewrite(
      object = assay_out,
      umi.object = object,
      features = scale_data_features,
      verbose = FALSE
    )
  )
  LayerData(assay_out, layer = "scale.data") <- residuals

  # Set the output's variable features based on consensus of all layers
  VariableFeatures(assay_out) <- VariableFeatures(
    assay_out, 
    use.var.features = FALSE,
    nfeatures = variable.features.n
  )

  return (assay_out)
}
