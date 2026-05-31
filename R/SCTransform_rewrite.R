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
    },
    'default' = {
      vst.args[['return_corrected_umi']] <- FALSE
      vst.args[['residual_type']] <- 'none'
      vst.out <- do.call(what = 'vst', args = vst.args)
      vst.out
    })

  # get residuals
  vst.out <- switch(
    EXPR = sct.method,
    'reference.model' = {
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
    },
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
        n_threads = 1L,
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
        n_threads = 1L
      )
      dimnames(x = vst.out$y) <- list(scale.data.features, colnames(x = umi))
      
      vst.out
    })

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
      return (layer_counts)
    }
  )

  # Apply SCTransform to each set of counts in `input_list`.
  output_list <- lapply(
    input_list,
    function(input) {
      layer_counts <- LayerData(object, layer = layer_name)
      layer_metadata <- cell.attr[coln(input), , drop = FALSE]
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
        seed.use = seed.use,
        verbose = verbose,
        ...
      )
    }
  )
  
  # Merge output assays into one, or take the single result.
  if (length(output_list) > 1) {
    assay_out <- merge(
      output_list[[1]], 
      output_list[-1]
    )
  } else {
    assay_out <- output_list[[1]]
  }
  
  # Determine which features to include in the output's scale.data slot.
  if (return.only.var.genes) {
    # Take the union of variable features across all output assays/layers.
    var_features_union <- Reduce(
      union,
      lapply(
        output_list,
        function(output) {
          return(VariableFeatures(output))
        }
      )
    )
    # Take the intersection of all features across all output assays/layers.
    all_features_intersect <- Reduce(
      intersect,
      lapply(
        output_list,
        function(output) {
          return(rownames(output))
        }
      )
    )
    # Keep features that are variable in at least one output assay/layer but
    # present in all of them.
    scale_data_features <- intersect(all_features_intersect, var_features_union)
  } else {
    # Use every feature found in any output assay/layer,
    scale_data_features <- Reduce(
      union,
      lapply(
        output_list,
        function(output) {
          return(rownames(output))
        }
      )
    )
  }
  
  # Extract residuals for the selected features and store them in
  # the outputs scaled.data slot.
  residuals <- suppressWarnings(
    FetchResiduals(
      object = assay_out, 
      umi.object = object,
      features = scale_data_features,
      verbose = FALSE
    )
  )
  LayerData(assay_out, layer = "scale.data") <- residuals

  # Set the output's variable features.
  VariableFeatures(assay_out) <- VariableFeatures(
    assay_out, 
    use.var.features = FALSE,
    nfeatures = variable.features.n
  )

  return (assay_out)
}
