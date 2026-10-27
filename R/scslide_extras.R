#' Sample-level ComBat Batch Correction
#'
#' Performs a control-only batch-mean correction using ComBat-style empirical
#' Bayes shrinkage across genes, with regression or subtraction correction.
#'
#' This function estimates batch effects from control samples only and corrects
#' all samples via regression or subtraction. It standardises the control means
#' to z-scores, computes moment-based priors \emph{across genes} (one set of
#' hyper-parameters per batch), and shrinks via parametric EM
#' (or non-parametric MC) estimation.
#'
#' @param object A Seurat object (sample-level).
#' @param assay Name of the assay to correct.
#'   Default is \code{DefaultAssay(object)}.
#' @param layer Layer within the assay that contains the data to correct.
#'   Default is \code{"data"}.
#' @param batch_key Column name in \code{object@@meta.data} identifying the
#'   batch (e.g., plate, pool).
#' @param condition_key Column name in \code{object@@meta.data} identifying the
#'   experimental condition.
#' @param control_label Value(s) in the \code{condition_key} column that mark
#'   control samples.
#' @param par.prior Logical; if \code{TRUE} (default), uses parametric
#'   empirical Bayes estimation (EM). If \code{FALSE}, uses non-parametric
#'   Monte Carlo integration.
#' @param method Correction method. \code{"regression"} (default) estimates a
#'   per-sample scaling factor \eqn{\beta_i} via OLS regression of each
#'   sample's profile on the shrunk batch mean. \code{"subtraction"} directly
#'   subtracts the shrunk batch mean (\eqn{\beta = 1} for all samples).
#' @param correct.attenuation Logical; if \code{TRUE} the OLS regression slope
#'   is corrected for attenuation caused by measurement error in the shrunk
#'   batch means.  Default \code{FALSE}.  Only used when
#'   \code{method = "regression"}.
#' @param floor.beta Floor applied to per-sample \eqn{\beta} estimates.
#'   Default \code{NULL} (no floor).  Set to e.g. \code{0} to prevent negative
#'   betas.  Only used when \code{method = "regression"}.
#' @param verbose Print progress messages. Default \code{TRUE}.
#'
#' @return The Seurat object with the corrected data written back to the
#'   specified assay layer.
#'
#' @export
#' @concept batch_correction
#'
#'
#' @examples
#' \dontrun{
#' obj <- CorrectSampleComBat(obj, assay = "LMC",
#'                                   batch_key = "plate",
#'                                   condition_key = "Condition",
#'                                   control_label = "ctrl-inj")
#' }
#'
#' @importFrom SeuratObject LayerData LayerData<- DefaultAssay Cells
#' @importFrom matrixStats rowVars
#' @importFrom stats dgamma median setNames var
#'
CorrectSampleComBat <- function(object,
                                      assay             = NULL,
                                      layer             = "data",
                                      batch_key         = NULL,
                                      condition_key     = NULL,
                                      control_label     = NULL,
                                      par.prior         = TRUE,
                                      method            = c("regression", "subtraction"),
                                      correct.attenuation = FALSE,
                                      floor.beta        = NULL,
                                      verbose           = TRUE) {

  # ---- 1. Input validation + data extraction --------------------------------
  if (!inherits(object, "Seurat")) {
    stop("'object' must be a Seurat object")
  }

  assay <- assay %||% DefaultAssay(object = object)
  if (!assay %in% names(object@assays)) {
    stop("Assay '", assay, "' not found in object")
  }

  if (is.null(batch_key)) {
    stop("Please specify 'batch_key': the metadata column identifying batches")
  }
  if (!batch_key %in% colnames(object@meta.data)) {
    stop("Column '", batch_key, "' not found in object metadata")
  }

  if (is.null(condition_key)) {
    stop("Please specify 'condition_key': the metadata column identifying conditions")
  }
  if (!condition_key %in% colnames(object@meta.data)) {
    stop("Column '", condition_key, "' not found in object metadata")
  }

  if (is.null(control_label)) {
    stop("Please specify 'control_label': the value(s) marking control samples in '",
         condition_key, "'")
  }

  if (!is.null(floor.beta) && (!is.numeric(floor.beta) || length(floor.beta) != 1)) {
    stop("'floor.beta' must be NULL or a single numeric value")
  }

  method <- match.arg(method)
  if (method == "subtraction") {
    if (isTRUE(correct.attenuation)) {
      warning("'correct.attenuation' is ignored when method = \"subtraction\"")
    }
    if (!is.null(floor.beta) && floor.beta != 0) {
      warning("'floor.beta' is ignored when method = \"subtraction\"")
    }
  }

  # --- retrieve the data matrix ---
  data_mat <- LayerData(object = object[[assay]], layer = layer)

  all_cells  <- Cells(object[[assay]])
  conditions <- object@meta.data[all_cells, condition_key]
  batches    <- object@meta.data[all_cells, batch_key]
  is_control <- conditions %in% control_label

  # global control indices
  global_ctrl_idx <- all_cells[is_control]
  if (length(global_ctrl_idx) == 0) {
    stop("No control samples found with '", condition_key, "' in ",
         paste(control_label, collapse = ", "))
  }

  G <- nrow(data_mat)
  unique_batches <- unique(batches)
  B <- length(unique_batches)

  # ---- 2. Per-batch control means Z_b (G x B matrix) -----------------------
  ctrl_idx_list <- setNames(
    lapply(unique_batches, function(b) all_cells[batches == b & is_control]),
    unique_batches
  )
  n_ctrl <- vapply(ctrl_idx_list, length, integer(1))  # named integer vector

  # Z_batch: G x B matrix (NA for batches with 0 controls)
  Z_batch <- matrix(NA_real_, nrow = G, ncol = B,
                    dimnames = list(rownames(data_mat), unique_batches))
  for (b in unique_batches) {
    nb <- n_ctrl[b]
    if (nb >= 2) {
      Z_batch[, b] <- rowMeans(data_mat[, ctrl_idx_list[[b]], drop = FALSE])
    } else if (nb == 1) {
      Z_batch[, b] <- data_mat[, ctrl_idx_list[[b]]]
    }
    # nb == 0 stays NA
  }

  # Global control mean (across all individual control samples)
  if (length(global_ctrl_idx) >= 2) {
    Z_global <- rowMeans(data_mat[, global_ctrl_idx])
  } else {
    Z_global <- data_mat[, global_ctrl_idx]
  }

  # ---- 3. Pooled within-batch variance (sigma2) -----------------------------
  batches_for_sigma <- unique_batches[n_ctrl >= 2]
  df_within <- sum(n_ctrl[batches_for_sigma] - 1)

  if (df_within >= 1) {
    SS_within <- rep(0, G)
    for (b in batches_for_sigma) {
      ctrl_cells <- ctrl_idx_list[[b]]
      resid <- data_mat[, ctrl_cells, drop = FALSE] - Z_batch[, b]
      SS_within <- SS_within + rowSums(resid^2)
    }
    sigma2 <- SS_within / df_within
  } else {
    warning("df_within < 1: cannot estimate within-batch variance; ",
            "setting sigma^2 = median of row variances across all controls")
    sigma2 <- matrixStats::rowVars(as.matrix(data_mat[, global_ctrl_idx, drop = FALSE]))
  }

  # Guard against zero sigma2 (would cause division by zero in standardisation)
  sigma2_safe <- sigma2
  sigma2_safe[sigma2_safe == 0] <- median(sigma2_safe[sigma2_safe > 0])

  # ---- 4. Standardise control means (ComBat-style) --------------------------
  batches_obs <- unique_batches[n_ctrl >= 1]
  B_obs <- length(batches_obs)

  if (B_obs < 2) {
    stop("Need at least 2 batches with control samples for EB shrinkage, found ", B_obs)
  }

  # gamma_hat[b, g] = (Z_batch[g,b] - Z_global[g]) / sqrt(sigma2_safe[g] / n_b)
  gamma_hat <- matrix(NA_real_, nrow = B_obs, ncol = G,
                      dimnames = list(batches_obs, rownames(data_mat)))
  for (i in seq_along(batches_obs)) {
    b <- batches_obs[i]
    gamma_hat[i, ] <- (Z_batch[, b] - Z_global) / sqrt(sigma2_safe / n_ctrl[b])
  }

  # delta_hat[b, g] = within-batch variance of standardised controls
  delta_hat <- matrix(1, nrow = B_obs, ncol = G,
                      dimnames = list(batches_obs, rownames(data_mat)))
  for (i in seq_along(batches_obs)) {
    b <- batches_obs[i]
    if (n_ctrl[b] >= 2) {
      # Standardise individual control samples: (x - Z_global) / sqrt(sigma2_safe)
      s_ctrl <- sweep(data_mat[, ctrl_idx_list[[b]], drop = FALSE], 1, Z_global) /
        sqrt(sigma2_safe)
      delta_hat[i, ] <- matrixStats::rowVars(s_ctrl)
    }
    # n_ctrl == 1: keep default of 1
  }

  # ---- 5. ComBat-style EB priors across genes -------------------------------
  gamma_bar <- rowMeans(gamma_hat)                        # length B_obs
  t2        <- apply(gamma_hat, 1, var)                   # length B_obs
  a_prior   <- apply(delta_hat, 1, .combat_aprior)        # length B_obs
  b_prior   <- apply(delta_hat, 1, .combat_bprior)        # length B_obs

  if (isTRUE(verbose)) {
    message("Fitting ", if (par.prior) "parametric" else "non-parametric",
            " EB priors on control samples (", B_obs, " batches, ", G, " genes)")
  }

  # ---- 6. EB shrinkage via EM -----------------------------------------------
  gamma_star <- matrix(NA_real_, nrow = B_obs, ncol = G,
                       dimnames = list(batches_obs, rownames(data_mat)))
  delta_star <- matrix(NA_real_, nrow = B_obs, ncol = G,
                       dimnames = list(batches_obs, rownames(data_mat)))

  for (i in seq_along(batches_obs)) {
    b <- batches_obs[i]

    # Standardised control data for this batch: (G x n_b)
    s_ctrl <- sweep(data_mat[, ctrl_idx_list[[b]], drop = FALSE], 1, Z_global) /
      sqrt(sigma2_safe)

    if (par.prior) {
      tmp <- .combat_it_sol(
        sdat  = s_ctrl,
        g_hat = gamma_hat[i, ],
        d_hat = delta_hat[i, ],
        g_bar = gamma_bar[i],
        t2    = t2[i],
        a     = a_prior[i],
        b     = b_prior[i]
      )
    } else {
      tmp <- .combat_int_eprior(
        sdat  = s_ctrl,
        g_hat = gamma_hat[i, ],
        d_hat = delta_hat[i, ]
      )
    }
    gamma_star[i, ] <- tmp$gamma_star
    delta_star[i, ] <- tmp$delta_star

    if (isTRUE(verbose)) {
      message("  Batch ", b, " done")
    }
  }

  # ---- 7. Convert gamma_star back to original scale -------------------------
  # Z_star[g, b] = Z_global[g] + gamma_star[b, g] * sqrt(sigma2_safe[g] / n_b)
  Z_star <- matrix(NA_real_, nrow = G, ncol = B,
                   dimnames = list(rownames(data_mat), unique_batches))
  for (b in unique_batches) {
    if (b %in% batches_obs) {
      idx <- match(b, batches_obs)
      Z_star[, b] <- Z_global + gamma_star[idx, ] * sqrt(sigma2_safe / n_ctrl[b])
    } else {
      # No controls in this batch: fall back to global mean (no batch effect)
      Z_star[, b] <- Z_global
      if (isTRUE(verbose)) {
        warning("Batch '", b, "' has no control samples; using global control mean")
      }
    }
  }

  if (method == "regression") {
  # ---- 8. Per-sample regression of Z*_b with attenuation correction ----------
  beta_vec <- setNames(rep(NA_real_, length(all_cells)), all_cells)

  # Track per-batch attenuation diagnostics
  atten_raw_denom <- setNames(rep(NA_real_, B), unique_batches)
  atten_me_var    <- setNames(rep(NA_real_, B), unique_batches)
  atten_denom     <- setNames(rep(NA_real_, B), unique_batches)
  atten_clamped   <- setNames(rep(FALSE, B), unique_batches)

  for (b in unique_batches) {
    batch_cells <- all_cells[batches == b]
    if (length(batch_cells) == 0) next

    z_star_b <- Z_star[, b]       # G-length vector
    z_bar    <- mean(z_star_b)
    z_cent   <- z_star_b - z_bar  # centred Z*

    raw_denom <- sum(z_cent^2)
    atten_raw_denom[b] <- raw_denom

    # Attenuation correction
    if (isTRUE(correct.attenuation) && any(sigma2 > 0)) {
      nb <- n_ctrl[b]
      # For v2 the measurement error in Z_star comes from the EB estimate
      # Approximate: use sigma2 / n_b as upper bound on ME variance per gene
      me_var <- sum(sigma2_safe / max(nb, 1))
      corrected_denom <- raw_denom - me_var
      clamped <- corrected_denom < 0.5 * raw_denom
      denom <- max(corrected_denom, 0.5 * raw_denom)
      atten_me_var[b]  <- me_var
      atten_denom[b]   <- denom
      atten_clamped[b] <- clamped
    } else {
      atten_me_var[b] <- 0
      atten_denom[b]  <- raw_denom
      denom <- raw_denom
    }

    if (denom < .Machine$double.eps) {
      beta_vec[batch_cells] <- 0
      next
    }

    for (cell in batch_cells) {
      x_i    <- data_mat[, cell]
      x_bar  <- mean(x_i)
      numer  <- sum(z_cent * (x_i - x_bar))
      beta_vec[cell] <- numer / denom
    }
  }

  # Apply floor.beta if requested
  n_floored <- setNames(integer(B), unique_batches)
  if (!is.null(floor.beta)) {
    for (b in unique_batches) {
      batch_cells <- all_cells[batches == b]
      below <- beta_vec[batch_cells] < floor.beta
      n_floored[b] <- sum(below, na.rm = TRUE)
      beta_vec[batch_cells[below]] <- floor.beta
    }
    if (isTRUE(verbose)) {
      total_floored <- sum(n_floored)
      message("floor.beta = ", floor.beta, ": floored ", total_floored,
              " / ", length(all_cells), " samples")
    }
  }

  # ---- 9. Apply regression correction ---------------------------------------
  for (b in unique_batches) {
    batch_cells <- all_cells[batches == b]
    if (length(batch_cells) == 0) next
    z_star_b <- Z_star[, b]
    for (cell in batch_cells) {
      data_mat[, cell] <- data_mat[, cell] - beta_vec[cell] * z_star_b
    }
  }

  } else {
  # ---- 8-9. Simple subtraction of shrunk batch mean --------------------------
  beta_vec <- setNames(rep(1, length(all_cells)), all_cells)
  for (b in unique_batches) {
    batch_cells <- all_cells[batches == b]
    if (length(batch_cells) == 0) next
    data_mat[, batch_cells] <- data_mat[, batch_cells] - Z_star[, b]
  }
  }

  # ---- 10. Verbose output and write back ------------------------------------
  if (isTRUE(verbose)) {
    # EB shrinkage summary
    shrink_df <- data.frame(
      batch      = batches_obs,
      n_ctrl     = n_ctrl[batches_obs],
      gamma_bar  = round(gamma_bar, 4),
      t2         = round(t2, 4),
      a_prior    = round(a_prior, 4),
      b_prior    = round(b_prior, 4),
      stringsAsFactors = FALSE
    )
    message("--- EB shrinkage summary ---")
    message(paste(utils::capture.output(print(shrink_df, row.names = FALSE)),
                  collapse = "\n"))

    if (method == "regression") {
    # Regression summary
    reg_df <- data.frame(
      batch     = unique_batches,
      n_samples = vapply(unique_batches,
                         function(b) sum(batches == b), integer(1)),
      median_beta = vapply(unique_batches,
                           function(b) median(beta_vec[all_cells[batches == b]]), numeric(1)),
      mean_beta   = vapply(unique_batches,
                           function(b) mean(beta_vec[all_cells[batches == b]]),   numeric(1)),
      sd_beta     = vapply(unique_batches,
                           function(b) {
                             vals <- beta_vec[all_cells[batches == b]]
                             if (length(vals) < 2) return(NA_real_)
                             sd(vals)
                           }, numeric(1)),
      stringsAsFactors = FALSE
    )
    if (!is.null(floor.beta)) {
      reg_df$n_floored <- n_floored[unique_batches]
    }
    message("--- Regression summary ---")
    message(paste(utils::capture.output(print(reg_df, row.names = FALSE)),
                  collapse = "\n"))

    # Attenuation correction summary
    if (isTRUE(correct.attenuation) && any(sigma2 > 0)) {
      atten_df <- data.frame(
        batch        = unique_batches,
        total_var    = round(atten_raw_denom[unique_batches], 4),
        noise_var    = round(atten_me_var[unique_batches], 4),
        signal_var   = round(atten_denom[unique_batches], 4),
        noise_frac   = round(atten_me_var[unique_batches] /
                               pmax(atten_raw_denom[unique_batches],
                                    .Machine$double.eps), 4),
        clamped      = atten_clamped[unique_batches],
        stringsAsFactors = FALSE
      )
      message("--- Attenuation correction summary ---")
      message("  total_var  = Var_genes(Z*_b)  [raw OLS denominator]")
      message("  noise_var  = estimated measurement error in Z*")
      message("  signal_var = total_var - noise_var  [corrected denominator, clamped at 50%]")
      message(paste(utils::capture.output(print(atten_df, row.names = FALSE)),
                    collapse = "\n"))
    }
    }
  }

  LayerData(object[[assay]], layer = layer) <- data_mat
  if (isTRUE(verbose)) {
    message("EB shrinkage + ", method, " batch correction (v2) completed for assay '",
            assay, "'")
  }

  return(object)
}


# ==============================================================================
# Internal ComBat helper functions
# Adapted from sva::ComBat (Johnson et al. 2007)
# These are vendored here to avoid a dependency on the sva package.
# ==============================================================================

#' Inverse-gamma shape prior for ComBat
#' @noRd
.combat_aprior <- function(delta_hat) {
  m <- mean(delta_hat)
  s2 <- var(delta_hat)
  # When all delta_hat are identical (s2 = 0, e.g. single-control batches),
  # the prior is degenerate.  Return a large finite shape that concentrates
  # the inverse-gamma prior tightly around m (prior mean = b/(a-1) = m).
  if (s2 < .Machine$double.eps) return(1002)
  (2 * s2 + m^2) / s2
}

#' Inverse-gamma rate prior for ComBat
#' @noRd
.combat_bprior <- function(delta_hat) {
  m <- mean(delta_hat)
  s2 <- var(delta_hat)
  if (s2 < .Machine$double.eps) return(1001 * m)
  (m * s2 + m^3) / s2
}

#' Posterior mean for ComBat
#' @noRd
.combat_postmean <- function(g_hat, g_bar, n, d_star, t2) {
  (t2 * n * g_hat + d_star * g_bar) / (t2 * n + d_star)
}

#' Posterior variance for ComBat
#' @noRd
.combat_postvar <- function(sum2, n, a, b) {
  (0.5 * sum2 + b) / (n / 2 + a - 1)
}

#' Parametric EM solver for ComBat
#' @noRd
.combat_it_sol <- function(sdat, g_hat, d_hat, g_bar, t2, a, b,
                           conv = 1e-4, max_iter = 5000) {
  n <- ncol(sdat)
  g_old <- g_hat
  d_old <- d_hat
  change <- 1
  count <- 0

  while (change > conv && count < max_iter) {
    g_new <- .combat_postmean(g_hat, g_bar, n, d_old, t2)
    sum2 <- rowSums((sdat - g_new)^2)
    d_new <- .combat_postvar(sum2, n, a, b)

    change <- max(abs(g_new - g_old) / pmax(abs(g_old), 1e-10),
                  abs(d_new - d_old) / pmax(abs(d_old), 1e-10),
                  na.rm = TRUE)
    g_old <- g_new
    d_old <- d_new
    count <- count + 1
  }

  list(gamma_star = g_new, delta_star = d_new)
}

#' Non-parametric empirical prior for ComBat (Monte Carlo integration)
#' @noRd
.combat_int_eprior <- function(sdat, g_hat, d_hat) {
  n <- ncol(sdat)
  g_star <- d_star <- numeric(length(g_hat))

  for (i in seq_along(g_hat)) {
    # Density of g_hat[i] given each g_hat[j] as the "true" value
    # using normal kernel with variance d_hat[i]/n
    g_densities <- dnorm(g_hat[i], g_hat, sqrt(d_hat[i] / n))
    # Weighted mean: posterior mean of gamma
    if (sum(g_densities) > 0) {
      g_star[i] <- sum(g_hat * g_densities) / sum(g_densities)
    } else {
      g_star[i] <- g_hat[i]
    }

    # For delta: use inverse-gamma kernel
    # For each possible true gamma g_hat[j], compute sum_s (sdat[i,s] - g_hat[j])^2
    x_i <- as.numeric(sdat[i, ])
    dat_rep <- matrix(x_i, nrow = length(g_hat), ncol = n, byrow = TRUE)
    sum2_i <- rowSums((dat_rep - g_hat)^2)
    # inverse-gamma density for d_hat[i]
    d_densities <- dgamma(1 / d_hat[i], shape = n / 2,
                          rate = sum2_i / 2) / d_hat[i]^2
    if (sum(d_densities) > 0) {
      d_star[i] <- sum(d_hat * d_densities) / sum(d_densities)
    } else {
      d_star[i] <- d_hat[i]
    }
  }

  list(gamma_star = g_star, delta_star = d_star)
}


#' Perform CellAnova Batch correction
#'
#' CellAnova (Zhaojun Zhang et al, 2024 Nat Biotech) is a new method that can remove/mitigate batch effect in single-cell data
#' and return users with a "corrected" data matrix (instead of just a low-dim embedding). Here we implement the 2nd step of the CellAnova.
#' Specifically, this function is an R-implemention of the CellAnova python package's calc_BE() function
#' (https://github.com/Janezjz/cellanova).
#'
#' This function takes a Seurat object and its pre-computed integrated embedding from
#' methods like Harmony or Seurat-CCA, a batch-effect index, and a case-control index as input, to
#' estimate the batch effect from the control samples, and correct for it from the full original
#' expression data.
#' Most of the procedures are kept the same, with the following modifications:
#' \itemize{
#'   \item we have additionally implemented a future_lapply() and a more efficient regression framework to enhance the efficiency.\cr
#'   \item the procedure can be done to a "sketched" data and later project to the whole data (for the purposes of efficiency and data balance).\cr
#' }
#'
#' @param object A Seurat object
#' @param assay the name of the assay to perform the CellAnova correction.
#' @param layer the name of the layer to be used to correct for the batch effect. Should be scale.data.
#' @param integrate_key A string indicating the smallest batch unit in the meta-data (e.g., library, donor, etc.),
#' which will be used for integration later.
#' @param features features to compute corrected expression for. Defaults to the variable features set in the assay specified.
#' @param control_dict A named list (or character vector) indicating which batches are controls.
#' Each element is a character vector of batch names (values in the \code{integrate_key} column)
#' that belong to one control group. The names of the list are arbitrary control-group labels
#' (e.g. \code{list(ctrl = c("batch1", "batch2"))}). A plain character vector is treated as a
#' single control group.
#' @param reduction the name of the DimReduc object we use as the integrated embeddings. Should be from methods like
#' Harmony or Seurat-integration methods (e.g., CCA).
#' @param var_cutoff the fraction of explained variance to determine the optimal value of k in truncated SVD when calculating
#' the basis of the batch effect. Default is 0.9.
#' @param k_max the maximum of singular values and vectors to compute.
#' @param k_select the user-defined number of singular values and vectors to compute (override var_cutoff and k_max).
#' Default is NULL.
#' @param new.assay.name the name for the new assay to store the corrected expression matrix
#' @param max_core the number of cores to use for parallel processing. Default is 1 (sequential).
#' @param future.memory.per.core the maximum memory (in MB) allowed per core for parallel processing. Default is 2000.
#' @param store.diagnostics Logical; if \code{TRUE} store the intermediate quantities
#' (\code{M_overall}, \code{DD1}, \code{VV1T}) in the object's \code{@@tools} slot for inspection/debugging.
#' Default is \code{FALSE} to avoid holding these (\code{VV1T} is genes x k) in the returned object.
#' @param verbose display progress + messages
#' @return Returns a Seurat object with a new assay added containing the batch-corrected expression matrix.
#' When the input layer is an on-disk BPCells IterableMatrix, the corrected assay is stored as a lazy
#' BPCells transform (no data materialised to disk).
#'
#' @importFrom irlba irlba
#' @importFrom RSpectra svds
#' @importFrom SeuratObject LayerData Embeddings CreateAssay5Object
#' @importFrom future.apply future_lapply
#' @importFrom methods slot slot<- new
#' @export
#' @concept batch_correction
cellanova_calc_BE <- function(object = NULL, assay = NULL, layer = "scale.data", integrate_key = NULL,
                              features = NULL, control_dict = NULL, reduction = NULL, var_cutoff = 0.9, k_max = 1500, k_select = NULL,
                              new.assay.name = "CORRECTED", max_core = 1, future.memory.per.core = 2000,
                              store.diagnostics = FALSE, verbose = TRUE){

  # Input validation
  if (is.null(object) || !inherits(object, "Seurat")) {
    stop("'object' must be a Seurat object")
  }
  if (!is.null(assay) && !assay %in% names(object@assays)) {
    stop("Assay '", assay, "' not found in object")
  }
  if (is.null(x = reduction)) {
    stop('Please provide dimensionality reduction name')
  }
  if (!reduction %in% names(object@reductions)) {
    stop("Reduction '", reduction, "' not found in object")
  }
  if (is.null(x = integrate_key)) {
    stop("Please specify the colname of the integration unit")
  }
  if (!integrate_key %in% colnames(object@meta.data)) {
    stop("Column '", integrate_key, "' not found in object metadata")
  }
  if (is.null(x = control_dict)) {
    stop("Please provide a list or a vector indicating the control samples")
  }
  if (!is.numeric(var_cutoff) || var_cutoff <= 0 || var_cutoff > 1) {
    stop("'var_cutoff' must be a numeric value between 0 and 1")
  }
  if (!is.numeric(k_max) || k_max <= 0) {
    stop("'k_max' must be a positive integer")
  }

  assay <- assay %||% DefaultAssay(object = object )
  features <- features %||% VariableFeatures(object = object[[assay]])
  if (length(x = features) == 0) {
    features <- rownames(x = LayerData(object = object[[assay]], layer = layer))
  }
  if (length(features) == 0) {
    stop("No features found. Please specify features or ensure the assay has variable features or data.")
  }

  # Detect whether the layer is an on-disk BPCells IterableMatrix
  GEX_full <- LayerData(object = object[[assay]], layer = layer)
  ondisk <- inherits(GEX_full, "IterableMatrix")

  if (ondisk && !requireNamespace("BPCells", quietly = TRUE)) {
    stop("Package 'BPCells' is required for on-disk CellAnova correction. ",
         "Please install it with: BiocManager::install('BPCells')")
  }

  # Check that the resulting matrix will not exceed R's 2^31-1 element limit (dense path only)
  n_features <- nrow(GEX_full)
  n_cells <- ncol(GEX_full)
  n_elements <- as.double(n_features) * as.double(n_cells)
  if (!ondisk && n_elements > .Machine$integer.max) {
    stop("The expression matrix (", n_features, " features x ", n_cells, " cells = ",
         format(n_elements, big.mark = ",", scientific = FALSE), " elements) exceeds R's ",
         "maximum dense matrix size (2^31 - 1 = ", format(.Machine$integer.max, big.mark = ","),
         " elements). Many internal R and Seurat operations use 32-bit indexing and will fail ",
         "on matrices this large. Consider reducing the number of cells (e.g., sketching) or ",
         "features before running CellAnova.")
  }

  if (ondisk && isTRUE(verbose)) {
    message("On-disk mode: BPCells IterableMatrix detected")
  }

  if(!is.list(control_dict)){
    control_dict <- list(g1 = control_dict)
  }
  control_groups <- names(control_dict)

  # Validate control samples exist in the data
  all_batches <- unique(object@meta.data[[integrate_key]])
  missing_batches <- setdiff(unlist(control_dict), all_batches)
  if (length(missing_batches) > 0) {
    stop("Control batches not found in data: ", paste(missing_batches, collapse = ", "))
  }

  # get the integration embeddings and scale each one
  LL <- Embeddings(object = object, reduction = reduction)

  # Remove columns that contain all zeros (more efficient vectorized approach)
  zero_cols <- which(colSums(abs(LL) > .Machine$double.eps) == 0)
  if(length(zero_cols) > 0){
    LL <- LL[, -zero_cols, drop = FALSE]
    if(isTRUE(verbose)) {
      message("Removed ", length(zero_cols), " zero columns from embeddings")
    }
  }
  # if the Embeddings are from a different assay, the cells might be different from those in 'assay'
  # need to check for this
  if( nrow(LL) != ncol(object[[assay]]) ){
    warning("Some of the cells in the reduction do not match to those in the assay.")
    LL <- LL[colnames(object[[assay]]), ]
  }

  # get the cell idendity by donor id
  Idents(object) <- integrate_key
  cell_donor_list <- CellsByIdentities(object = object, cells = colnames(object[[assay]]))

  ### regression 1: calculate the overall mean effect
  batch_list <- unique(object[[integrate_key]][colnames(object[[assay]]),1])

  # Check that every batch has at least as many cells as embedding dimensions.
  # qr.coef() returns NA for underdetermined systems (n_cells < n_pcs), and
  # those NAs silently propagate through M_overall into the final correction.
  n_pcs <- ncol(LL)
  small_batches <- vapply(batch_list, function(d) {
    length(cell_donor_list[[d]]) < n_pcs
  }, logical(1))
  if (any(small_batches)) {
    bad <- batch_list[small_batches]
    bad_sizes <- vapply(bad, function(d) length(cell_donor_list[[d]]), integer(1))
    stop("The following batches in '", integrate_key, "' have fewer cells than ",
         "embedding dimensions (", n_pcs, "), which would produce NA regression ",
         "coefficients:\n",
         paste0("  ", bad, " (", bad_sizes, " cells)", collapse = "\n"),
         "\nPlease remove these batches or reduce the number of dimensions in '",
         reduction, "'.")
  }

  # GEX_full was already extracted above (during on-disk detection)
  # For the dense path, pre-transpose once: tGEX_full is (n_cells x n_genes)
  if (!ondisk) {
    tGEX_full <- t(GEX_full)
  }

  # Set up future parallelization
  original_plan <- future::plan()
  user_has_custom_plan <- !inherits(original_plan, "sequential")
  future_modified <- FALSE

  if (!user_has_custom_plan) {
    if (max_core > 1) {
      if (isTRUE(verbose)) {
        if (.Platform$OS.type == "windows") {
          message("Setting up future multisession with ", max_core, " workers (Windows)")
        } else {
          message("Setting up future multicore with ", max_core, " workers")
        }
        message("Setting future.globals.maxSize to ",
                max_core * future.memory.per.core, " MB")
      }
      if (.Platform$OS.type == "windows") {
        future::plan(future::multisession, workers = max_core)
      } else {
        future::plan(future::multicore, workers = max_core)
      }
      options(future.globals.maxSize = max_core * future.memory.per.core * 1024^2)
      future_modified <- TRUE
    } else {
      if (isTRUE(verbose)) {
        message("Using sequential future plan (no parallelization)")
      }
      future::plan("sequential")
      future_modified <- TRUE
    }
  } else {
    if (isTRUE(verbose)) {
      message("Using user-specified future plan")
    }
  }

  if (future_modified) {
    on.exit({
      future::plan(original_plan)
      if (isTRUE(verbose)) message("Restored original future plan")
    }, add = TRUE)
  }

  if (ondisk) {
    # On-disk path: materialise each batch subset individually (no full transpose)
    # Use sequential lapply -- future::multicore cannot serialise IterableMatrix pointers
    overall_bloc_coef <- lapply(batch_list, function(d) {
      cells <- cell_donor_list[[d]]
      if (length(cells) < 1) return(NULL)
      LL_tmp <- LL[cells, , drop = FALSE]
      qr_LL <- qr(LL_tmp)
      # Materialise only this batch's subset from on-disk directly in transposed form
      tGEX_mat <- as.matrix(t(GEX_full[, cells]))
      base::t(qr.coef(qr_LL, tGEX_mat))
    })
  } else {
    # Dense path: use pre-transposed full matrix with future_lapply
    overall_bloc_coef <- future_lapply(batch_list, function(d) {
      if (length(cell_donor_list[[d]]) < 1) return(NULL)
      LL_tmp <- LL[cell_donor_list[[d]], , drop = FALSE]
      qr_LL <- qr(LL_tmp)
      tGEX_mat <- tGEX_full[cell_donor_list[[d]], , drop = FALSE]
      t(qr.coef(qr_LL, tGEX_mat))
    })
    # Free the transposed matrix -- no longer needed
    rm(tGEX_full)
  }
  names(overall_bloc_coef) <- batch_list
  overall_bloc_coef <- Filter(Negate(is.null), overall_bloc_coef)

  if (isTRUE(verbose)) {
    message("Done with regression for ", length(overall_bloc_coef), " batches in 'integrate_key'")
  }
  M_overall <- Reduce(`+`, overall_bloc_coef) / length(overall_bloc_coef)

  ### regression 2: calculate the batch specific effect from controls
  # an empty object to store the regression outputs
  res <- list()
  #
  for(g in control_groups){
    control_batch <- control_dict[[g]]
    counter <- 0
    bloc_coef <- list()
    # loop through each batch to compute the regression effect
    for(d in control_batch){
      # LL_tmp <- LL[cell_donor_list[[d]], ]
      # qr_LL <- qr(LL_tmp)
      # # get the coef for all the m genes
      # coef <- qr.coef(qr_LL, t(LayerData(object = object[[assay]], layer = layer, cells = cell_donor_list[[d]])))
      bloc_coef[[d]] <- overall_bloc_coef[[d]]
      counter <- counter + 1
    }

    # calculate the mean of it
    M <- Reduce(`+`, bloc_coef) / length(control_batch)
    #
    res_temp <- do.call(cbind, bloc_coef) - do.call(cbind, replicate(length(control_batch), M, simplify = FALSE))
    res[[g]] <- res_temp
  }
  #
  res_combined <- t(do.call(cbind, res))

  # Clean up intermediate objects to free memory
  rm(overall_bloc_coef, bloc_coef, res, res_temp, M)

  ## Perform SVD decomposition for res_combined
  # Helper: try RSpectra::svds(), fall back to base::svd() when ARPACK fails
  # (e.g. "TridiagEigen" errors on numerically rank-deficient matrices)
  safe_svds <- function(A, k, verbose = FALSE) {
    tryCatch(
      RSpectra::svds(A, k = k),
      error = function(e) {
        if (isTRUE(verbose)) {
          message("RSpectra::svds() failed (", conditionMessage(e),
                  "), falling back to base::svd()")
        }
        res <- base::svd(A, nu = 0, nv = k)
        list(d = res$d[seq_len(k)], v = res$v[, seq_len(k), drop = FALSE])
      }
    )
  }

  # Compute true total variance (Frobenius norm squared = sum of ALL singular values squared)
  # This is needed because truncated SVD only returns the top-k singular values,
  # so sum(top_k_sv^2) / sum(top_k_sv^2) always reaches 1.0 at the last component.
  total_var <- sum(res_combined^2)

  if(is.null(x = k_select)){
    k_max <- min(c(k_max, min(dim(res_combined)) - 1))

    # Two-stage SVD: run a small pilot first to estimate k
    k_pilot <- min(k_max, 100)
    k_pilot <- min(k_pilot, min(dim(res_combined)) - 1)

    if(isTRUE(verbose)) {
      message("Performing pilot SVD with k_pilot = ", k_pilot)
    }

    pilot_svds <- safe_svds(res_combined, k = k_pilot, verbose = verbose)
    DD1 <- pilot_svds$d

    # Calculate cumulative variance explained from pilot
    positive_sv <- DD1[DD1 > 0]
    if(length(positive_sv) == 0) {
      stop("No positive singular values found in SVD decomposition")
    }

    variance <- cumsum(positive_sv^2) / total_var
    k <- which(variance >= var_cutoff)[1]

    if (is.na(k) || k > k_pilot) {
      # Pilot was not enough -- fall back to full k_max SVD
      k_max <- min(k_max, min(dim(res_combined)) - 1)
      if(isTRUE(verbose)) {
        message("Pilot SVD insufficient, performing full SVD with k_max = ", k_max)
      }
      final_svds <- safe_svds(res_combined, k = k_max, verbose = verbose)
      DD1 <- final_svds$d
      positive_sv <- DD1[DD1 > 0]
      variance <- cumsum(positive_sv^2) / total_var
      k <- which(variance >= var_cutoff)[1]
      if (is.na(k)) {
        k <- length(positive_sv)
        warning("Could not find k components explaining ", var_cutoff * 100, "% variance. Using all ", k, " components.")
      }
    } else {
      final_svds <- pilot_svds
    }

    if(isTRUE(verbose)) {
      message("Selected k = ", k, " components explaining ", round(variance[k] * 100, 2), "% variance")
    }
  } else {
    k <- k_select
    if(k > min(dim(res_combined)) - 1) {
      k <- min(dim(res_combined)) - 1
      warning("k_select too large, reduced to ", k)
    }

    final_svds <- safe_svds(res_combined, k = k, verbose = verbose)

    if(isTRUE(verbose)) {
      message("Using user-specified k = ", k, " components")
    }
  }
  #

  # Extract selected components
  DD1 <- final_svds$d[1:k]
  VV1T <- final_svds$v[, 1:k, drop = FALSE]

  # Store intermediate results in object tools slot for debugging/inspection
  if (isTRUE(store.diagnostics)) {
    slot(object = object, name = "tools")[["M_overall"]] <- M_overall
    slot(object = object, name = "tools")[["DD1"]] <- DD1
    slot(object = object, name = "tools")[["VV1T"]] <- VV1T
  }

  if(isTRUE(verbose)) {
    message("Computing batch effect correction...")
  }

  if (ondisk) {
    # =========================================================================
    # On-disk path: lazy correction via BPCells TransformLinearResidual
    # =========================================================================
    # corrected = GEX_full - VV1T %*% (t(VV1T) %*% (GEX_full - M_overall %*% t(LL)))
    #
    # Rearranged as:  corrected = GEX_full - VV1T %*% correction_coef
    # where correction_coef = t(VV1T) %*% GEX_full - t(VV1T) %*% M_overall %*% t(LL)
    #
    # This is exactly X - t(row_params) %*% col_params, which BPCells'
    # TransformLinearResidual computes lazily (no full matrix materialised).

    tVV1T <- t(VV1T)  # (k x genes)

    # proj = t(VV1T) %*% GEX_full  -->  (k x cells), computed by streaming GEX_full once
    # BPCells: matrix %*% IterableMatrix returns a dense matrix
    if (isTRUE(verbose)) message("  Computing projection (streaming on-disk data)...")
    proj <- tVV1T %*% GEX_full  # (k x cells)

    # small_term = t(VV1T) %*% M_overall %*% t(LL)  -->  (k x cells) via tiny intermediates
    small_term <- (tVV1T %*% M_overall) %*% base::t(LL)  # (k x n_pcs) %*% (n_pcs x cells)

    # correction_coef = proj - small_term  -->  (k x cells)
    correction_coef <- proj - small_term
    rm(proj, small_term)

    # Build lazy IterableMatrix: corrected = GEX_full - t(row_params) %*% col_params
    row_params <- tVV1T            # (k x genes), reuse of t(VV1T) computed above
    col_params <- correction_coef  # (k x cells)

    # BPCells' TransformLinearResidual applies row_params/col_params in the
    # internal (pre-transpose) coordinate system.  When the matrix is stored
    # transposed (storage_order == "row"), we must swap the params -- exactly
    # as BPCells' own regress_out() does (see BPCells/R/transforms.R).
    if (GEX_full@transpose) {
      tmp <- row_params
      row_params <- col_params
      col_params <- tmp
    }

    if(isTRUE(verbose)) {
      message("Creating corrected assay (lazy on-disk transform)...")
    }

    corrected <- new("TransformLinearResidual",
      matrix = BPCells::convert_matrix_type(GEX_full, "double"),
      transpose = GEX_full@transpose,
      dim = GEX_full@dim,
      dimnames = dimnames(GEX_full),
      row_params = row_params,
      col_params = col_params,
      global_params = numeric(0),
      vars_to_regress = character(0)
    )

    new.assay <- CreateAssay5Object(data = corrected)
    object[[new.assay.name]] <- new.assay

  } else {
    # =========================================================================
    # Dense (in-memory) path -- original logic
    # =========================================================================
    # Compute residuals after removing overall mean effect
    part1 <- GEX_full - tcrossprod(M_overall, LL)

    # Project onto batch effect subspace and back to compute batch effect
    # The intermediate crossprod(VV1T, part1) is (k x cells), which is small
    be <- VV1T %*% crossprod(VV1T, part1)
    rm(part1)

    # Apply correction by subtracting estimated batch effect
    corrected <- GEX_full - be
    rm(be, GEX_full)

    if(isTRUE(verbose)) {
      message("Creating corrected assay...")
    }

    # Create new assay with corrected data
    new.assay <- suppressWarnings(
      expr = CreateAssayObject(
        data = corrected,
        min.cells = -Inf,
        min.features = -Inf,
        check.matrix = FALSE
      )
    )

    # Add corrected assay to object
    object[[new.assay.name]] <- new.assay
  }

  if(isTRUE(verbose)) {
    message("Batch correction completed. New assay '", new.assay.name, "' added to object.")
  }

  return(object)
}

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
