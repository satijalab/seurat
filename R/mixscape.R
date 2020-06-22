#' @include generics.R
#'
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Calculate a perturbation score
#' 
#' Function to calculate perturbation score for pooled CRISPR screen datasets. 
#' For each target cell (expressing one target gRNA), we identified 20 cells 
#' from the control pool (non-targeting cells) with the most similar mRNA 
#' expression profiles. The perturbation score is calculated by subtracting the 
#' averaged mRNA expression profile of the non-targeting neighbors from the mRNA 
#' expression profile of the target cell.
#'  
#' @param object An object of class Seurat.
#' @param assay Name of Assay PRTB  signature is being calculated on.
#' @param slot Data slot to use for PRTB score calculation.
#' @param gd.class Metadata column containing target gene classification.
#' @param nt.cell.class Non-targeting gRNA cell classification identity.
#' @param split.by Provide metadata column if multiple biological replicates 
#' exist to calculate PRTB score for every replicate separately.
#' @param num.neighbors Number of nearest neighbors to consider.
#' @param ndims Number of dimensions to use from dimensionality reduction method.
#' @param reduction Reduction method used to calculate nearest neighbors.
#' @param new.assay.name Name for the new assay.
#' @return Returns a Seurat object with a new assay added containing the 
#' perturbation scores for all cells in the data slot.
#' 
#' @importFrom RANN nn2
#' @export
#' 
CalcPerturbScore <- function( 
  object, 
  assay = NULL,
  slot = "data", 
  gd.class = "guide_ID",
  nt.cell.class = "NT",
  split.by = NULL,
  num.neighbors = NULL,
  reduction = "pca", 
  ndims = 15,
  new.assay.name = "PRTB"
) {
  assay <- assay %||% DefaultAssay(object = object )
  if (is.null(x = reduction)) {
    stop('Please provide dimensionality reduction name.')
  } 
  if (is.null(x = num.neighbors)) {
    stop("Please specify number of nearest neighbors to consider")
  }
  if (is.null(x = ndims)) {
    stop("Please provide number of ", reduction, " dimensions to consider")
  }
  if (! is.null(x = split.by)) {
    Idents(object = object) <-  split.by
  } else {
    Idents(object = object) <- "rep1"
  }
  replicate <- levels(x = object)
  all_diff <- matrix(nrow = length(x = rownames(x = GetAssayData(object = object, assay = assay, slot = slot))), ncol = 0)
  for (r in replicate) {
    rep1 <- object[, WhichCells(object = object, idents = r)]
    # isolate nt cells
    all_cells <- Cells(x = rep1)
    nt_cells <- Cells(x = rep1[, grep(pattern = nt.cell.class, x = rep1[[]][, gd.class], value = FALSE)])
    # subset the objects based on guide ID
    all <- rep1[, all_cells]
    nt <- rep1[, nt_cells]
    # get pca cell embeddings
    all_mtx <- Embeddings(object = all, reduction = reduction)
    nt_mtx <- Embeddings(object = nt, reduction = reduction)
    # run nn2 to find the 20 nearest NT neighbors for all cells. Use the same 
    # number of PCs as the ones you used for umap
    mtx <- nn2(
      data = nt_mtx[, 1:ndims], 
      query = all_mtx[, 1:ndims], 
      k = num.neighbors
    )
    nt_data <- expm1(x = GetAssayData(object = nt, assay = assay, slot = slot))
    new_expr <- pbsapply(X = 1:length(x = all_cells), FUN = function(i) {
      index <- mtx$nn.idx[i,]
      nt_cells20 <- nt_cells[index]
      avg_nt <- rowMeans(x = nt_data[, nt_cells20])
      avg_nt <- as.matrix(x = avg_nt)
      colnames(x = avg_nt) <- all_cells[i]
      return(avg_nt)
    })
    new_expr <- log1p(x = new_expr)
    rownames(x = new_expr) <- rownames(x = GetAssayData(object = nt, assay = assay, slot = slot))
    colnames(x = new_expr) <- all_cells
    diff <- new_expr - GetAssayData(object = rep1, slot = slot, assay = assay)[, colnames(x = new_expr)]
    all_diff <- cbind(all_diff, diff)
  } 
  all_diff <- as(Class = "dgCMatrix", object = all_diff)
  prtb.assay <- CreateAssayObject(data =  all_diff[, colnames(x = object)], min.cells = -Inf, min.features = -Inf)
  object[[new.assay.name]] <- prtb.assay
  return(object)
}

#' DE and EnrichR pathway visualization barplot
#' 
#' @inheritParams FindMarkers
#' @param object Name of object class Seurat.
#' @param ident.1 Cell class identity 1.
#' @param ident.2 Cell class identity 2. 
#' @param balanced Option to display pathway enrichments for both negative and 
#' positive DE genes.If false, only positive DE gene will be displayed.
#' @param max.genes Maximum number of genes to use as input to enrichR.
#' @param p.val.cutoff Cutoff to select DE genes.
#' @param cols A list of colors to use for barplots.
#' @param enrich.database Database to use from enrichR.
#' @param num.pathway Number of pathways to display in barplot.
#' @param return.gene.list Return list of DE genes
#' 
#' @return Returns one (only enriched) or two (both enriched and depleted) 
#' barplots with the top enriched/depleted GO terms from EnrichR. 
#' 
#' @importFrom ggplot2 ggplot geom_bar coord_flip scale_fill_manual ylab ggtitle
#' theme_classic theme element_text 
#' @importFrom patchwork wrap_plots
#' 
#' @export

DEenrichRPlot <- function(
  object,
  ident.1 = NULL,
  ident.2 = NULL,
  balanced = TRUE,
  logfc.threshold = 0.25,
  assay = NULL,
  max.genes,
  test.use = 'wilcox',
  p.val.cutoff = 0.05,
  cols = NULL,
  enrich.database = NULL,
  num.pathway = 10,
  return.gene.list = FALSE,
  ...
) { 
  enrichr.installed <- PackageCheck("enrichR", error = FALSE)
  if (!enrichr.installed[1]) {
    stop(
      "Please install the enrichR package to use DEenrichRPlot",
      "\nThis can be accomplished with the following command: ",
      "\n----------------------------------------",
      "\ninstall.packages('enrichR')",
      "\n----------------------------------------",
      call. = FALSE
    )
  }
  if (is.null(x = enrich.database)) {
    stop("Please specify the name of enrichR database to use")
  }
  if (!is.numeric(x = max.genes)) { 
    stop("please set max.genes")
  }
  assay <- assay %||% DefaultAssay(object = object)
  DefaultAssay(object = object) <- assay
  all.markers <- FindMarkers(
    object = object, 
    ident.1 = ident.1,
    ident.2 = ident.2, 
    only.pos = FALSE, 
    logfc.threshold = logfc.threshold,
    test.use = test.use, 
    assay = assay
  )
  pos.markers <- all.markers[all.markers[, 2] > logfc.threshold & all.markers[, 1] < p.val.cutoff, , drop = FALSE]
  pos.markers.list <- rownames(x = pos.markers)[1:min(max.genes, nrow(x = pos.markers))]
  pos.er <- enrichR::enrichr(genes = pos.markers.list, databases = enrich.database)
  pos.er <- do.call(what = cbind, args = pos.er)
  pos.er$log10pval <- -log10(x = pos.er[, paste(enrich.database, sep = ".", "P.value")])
  pos.er$term <- pos.er[, paste(enrich.database, sep = ".", "Term")]
  pos.er <- pos.er[1:num.pathway, ]
  pos.er$term <- factor(x = pos.er$term, levels = pos.er$term[order(pos.er$log10pval)])
  gene.list <- list(pos = pos.er)
  
  if (isTRUE(x = balanced)) {
    neg.markers <- all.markers[all.markers[, 2] < logfc.threshold & all.markers[, 1] < p.val.cutoff, , drop = FALSE]
    neg.markers.list <- rownames(x = neg.markers)[1:min(max.genes, nrow(x = neg.markers))]
    neg.er <- enrichR::enrichr(genes = neg.markers.list, databases = enrich.database)
    neg.er <- do.call(what = cbind, args = neg.er)
    neg.er$log10pval <- -log10(x = neg.er[, paste(enrich.database, sep = ".", "P.value")])
    neg.er$term <- neg.er[, paste(enrich.database, sep = ".", "Term")]
    neg.er <- neg.er[1:num.pathway, ]
    neg.er$term <- factor(x = neg.er$term, levels = neg.er$term[order(neg.er$log10pval)])
    gene.list <- list(pos = pos.er, neg = neg.er)
  }
  if (return.gene.list) {
    return(gene.list)
  }
  
  p <- ggplot(data = pos.er, aes_string(x = "term", y = "log10pval")) +
    geom_bar(stat = "identity") +
    coord_flip() + xlab("Pathway") +
    scale_fill_manual(values = cols, drop = FALSE) +
    ylab("-log10(pval)") +
    ggtitle(paste(enrich.database, ident.1, sep = "_", "positive markers")) +
    theme_classic() +
    theme(axis.text.y = element_text(size = 12, face = "bold"))
  if (isTRUE(x = balanced)) {
    p2 <- ggplot(data = neg.er, aes_string(x = "term", y = "log10pval")) +
    geom_bar(stat = "identity") +
    coord_flip() + xlab("Pathway") +
    scale_fill_manual(values = cols, drop = FALSE) +
    ylab("-log10(pval)") +
    ggtitle(paste(enrich.database, ident.1, sep = "_", "negative markers")) +
    theme_classic() +
    theme(axis.text.y = element_text(size = 12, face = "bold"))
    p <- wrap_plots(p, p2)
  }
  return(p)
}

#' Linear discriminant analysis on pooled CRISPR screen data.
#' 
#' This function performs unsupervised PCA on each mixscape class separately and projects each subspace onto all 
#' cells in the data. Finally, it uses the first 10 principle components from each projection as input to lda in MASS package together with mixscape class labels.
#' 
#' @inheritParams PrepLDA
#' @inheritParams RunLDA
#' 
#' @return Returns a Seurat object with LDA added in the reduction slot.
#' 
#' @export
#' 
MixscapeLDA <- function(
  object,
  assay = NULL,
  ndims.print = 1:5,
  nfeatures.print = 30,
  reduction.key = "LDA_",
  seed = 42,
  pc.assay = "PRTB",
  labels = "gene",
  nt.label = "NT",
  npcs = 10,
  verbose = TRUE,
  logfc.threshold = 0.25
) {
  projected_pcs <- PrepLDA(
    object = object, 
    de.assay = assay, 
    pc.assay = pc.assay, 
    labels = labels, 
    nt.label = nt.label, 
    npcs = npcs , 
    verbose = verbose
  )
  lda.lables <- object[[labels]][,]
  object_lda <- RunLDA(
    object = projected_pcs, 
    labels = lda.lables, 
    assay = assay, 
    verbose = verbose
  )
  object[["lda"]] <- object_lda
  return(object)
}

#' Function to prepare data for Linear Discriminant Analysis.
#' 
#' This function performs unsupervised PCA on each mixscape class separately and projects each subspace onto all 
#' cells in the data.
#' 
#' @param object An object of class Seurat.
#' @param de.assay Assay to use for selection of DE genes.
#' @param pc.assay Assay to use for running Principle components analysis.
#' @param labels Meta data column with target gene class labels.
#' @param nt.label Name of non-targeting cell class.
#' @param npcs Number of principle components to use.
#' @param verbose Print progress bar.
#' @inheritParams FindMarkers
#' @return Returns a list of the first 10 PCs from each projection.
#' 
#' @export
#' 
PrepLDA <- function( 
  object,
  de.assay = "RNA",
  pc.assay = "PRTB",
  labels = "gene",
  nt.label = "NT",
  npcs = 10,
  verbose = TRUE,
  logfc.threshold = 0.25
) {
  projected_pcs <- list()
  gene_list <- setdiff(x = unique(x = object[[labels]][, 1]), y = nt.label)
  Idents(object = object) <- labels
  for (g in gene_list) {
    if (verbose) {
      message(g)
    }
    gene_subset <- subset(x = object, idents = c(g, nt.label))
    DefaultAssay(object = gene_subset) <- pc.assay
    gene_set <- TopDEGenesMixscape(
      object = gene_subset,
      ident.1 = g,
      de.assay = de.assay,
      labels = labels,
      logfc.threshold = logfc.threshold
    )
    if (length(x = gene_set) < (npcs + 1)) {
      next
    }
    gene_subset <- ScaleData(
      object = gene_subset, 
      features = gene_set, 
      verbose = FALSE
    )
    gene_subset <- RunPCA(
      object = gene_subset,
      features = gene_set,
      npcs = npcs,
      verbose = FALSE
    )
    project_pca <- ProjectCellEmbeddings(
      reference = gene_subset,
      query = object,
      dims = 1:npcs
    )
    colnames(x = project_pca) <- paste(g, colnames(x = project_pca), sep = "_")
    projected_pcs[[g]] <- project_pca
  }
  return(projected_pcs)
}

#' @param object Input values for LDA (numeric), with observations as rows
#' @param labels Observation labels for LDA
#' @param assay Name of Assay LDA is being run on
#' @param ndims.print PCs to print genes for
#' @param nfeatures.print Number of genes to print for each PC
#' @param reduction.key dimensional reduction key, specifies the string before
#' the number for the dimension names. LDA by default
#' @param seed Set a random seed. By default, sets the seed to 42. Setting
#' NULL will not set a seed.
#'
#' @importFrom MASS lda
#' @importFrom stats predict
#'
#' @rdname RunLDA
#' @export
#' @method RunLDA default
#'
RunLDA.default <- function(
  object,
  labels,
  assay = NULL,
  verbose = TRUE,
  ndims.print = 1:5,
  nfeatures.print = 30,
  reduction.key = "LDA_",
  seed = 42,
  ...
) {
  if (!is.null(x = seed)) {
    set.seed(seed = seed)
  }
  object <- data.frame(object)
  var_names <- colnames(x = object)
  object$lda_cluster_label <- labels
  lda_results <- lda(formula = lda_cluster_label ~ ., data = object)
  lda_predictions <- predict(object = lda_results, newdata = object)
  lda_cv <-lda(
    formula = lda_cluster_label ~ ., 
    data = object, 
    CV = TRUE
  )$posterior
  feature.loadings <- lda_results$scaling
  cell.embeddings <- lda_predictions$x
  lda.assignments <- lda_predictions$class
  lda.posterior <- lda_predictions$posterior
  colnames(x = lda.posterior) <- paste0("LDAP_", colnames(x = lda.posterior))
  rownames(x = feature.loadings) <- var_names
  colnames(x = feature.loadings) <- paste0(reduction.key, 1:ncol(x = cell.embeddings))
  rownames(x = cell.embeddings) <- rownames(x = object)
  colnames(x = cell.embeddings) <- colnames(x = feature.loadings)
  reduction.data <- CreateDimReducObject(
    embeddings = cell.embeddings,
    loadings = feature.loadings,
    assay = assay,
    key = reduction.key,
    misc = list(
      assignments = lda.assignments, 
      posterior = lda.posterior, 
      model = lda_results, 
      cv = lda_cv
    )
  )
  if (verbose) {
    print(x = reduction.data, dims = ndims.print, nfeatures = nfeatures.print)
  }
  return(reduction.data)
}


#' Function to perform Linear Discriminant Analysis. 
#' 
#' @param ndims.print Number of LDA dimensions to print.
#' @param nfeatures.print Number of features to print for each LDA component.
#' @param reduction.key Reduction key name.
#' 
#' @rdname RunLDA
#' @export
#' @method RunLDA Assay
#'
RunLDA.Assay <- function(
  object,
  assay,
  labels,
  features = NULL,
  verbose = TRUE,
  ndims.print = 1:5,
  nfeatures.print = 30,
  reduction.key = "LDA_",
  seed = 42,
  ...
) {
  data.use <- PrepDR(
    object = object,
    features = features,
    verbose = verbose
  )
  reduction.data <- RunLDA(
    object = t(x = data.use),
    assay = assay,
    labels = labels,
    verbose = verbose,
    ndims.print = ndims.print,
    nfeatures.print = nfeatures.print,
    reduction.key = reduction.key,
    seed = seed,
    ...
  )
  return(reduction.data)
}

#' @param object An object of class Seurat.
#' @param assay Assay to use for performing Linear Discriminant Analysis (LDA).
#' @param labels Meta data column with target gene class labels.
#' @param features Features to compute LDA on
#' @param reduction.name dimensional reduction name, lda by default
#' @param reduction.key Reduction key name.
#' @param seed Value for random seed
#' @param verbose Print the top genes associated with high/low loadings for
#' the PCs
#' @param ndims.print Number of LDA dimensions to print.
#' @param nfeatures.print Number of features to print for each LDA component.
#' 
#' @rdname RunLDA
#' @export
#' @method RunLDA Seurat
#'
RunLDA.Seurat <- function(
  object,
  assay = NULL,
  labels,
  features = NULL,
  reduction.name = "lda",
  reduction.key = "LDA_",
  seed = 42,
  verbose = TRUE,
  ndims.print = 1:5,
  nfeatures.print = 30,
  ...
) {
  assay <- assay %||% DefaultAssay(object = object)
  assay.data <- GetAssay(object = object, assay = assay)
  reduction.data <- RunLDA(
    object = assay.data,
    labels = labels,
    features = features,
    verbose = verbose,
    ndims.print = ndims.print,
    nfeatures.print = nfeatures.print,
    reduction.key = reduction.key,
    seed = seed,
    ...
  )
  object[[reduction.name]] <- reduction.data
  object$lda.assignments <- slot(object = object[[reduction.name]], name = "misc")[["assignments"]]
  object <- AddMetaData(
    object = object,
    metadata = as.data.frame(
      x = slot(object = object[[reduction.name]], name = "misc")[["posterior"]]
    )
  )
  object <- LogSeuratCommand(object = object)
  object <- ProjectDim(
    object = object,
    reduction = reduction.name,
    assay = assay,
    verbose = verbose,
    dims.print = ndims.print,
    nfeatures.print = nfeatures.print
  )
  Loadings(object = object[[reduction.name]]) <- Loadings(
    object = object[[reduction.name]], 
    projected = TRUE
  )
  return(object)
}

#' Function to identify perturbed and non-perturbed gRNA expressing cells.
#' 
#' @inheritParams FindMarkers
#' @param object An object of class Seurat.
#' @param assay Assay to use for mixscape classification.
#' @param slot Assay data slot to use.
#' @param labels metadata column with target gene classifications.
#' @param nt.class.name Classification name of non-targeting gRNA cells.
#' @param new.class.name Name of mixscape classification to be stored in 
#' metadata.
#' @param min.de.genes Required number of genes that are differentially 
#' expressed for method to separate perturbed and non-perturbed cells. 
#' @param de.assay Assay to use when performing differential expression analysis. 
#' Usually RNA.
#' @param iter.num Number of normalmixEM iterations to run if convergence does 
#' not occur.
#' @param verbose Display messages
#' @return Returns Seurat object with with the following information in the 
#' meta data:
#' \describe{
#'   \item{mixscape_class}{Classification result with cells being either 
#'   classified as perturbed (KO) or non-perturbed (NP) based on their target 
#'   gene class.}
#'   \item{mixscape_class.global}{Global classification result (KO, NP or NT)}
#'   \item{p_ko}{Posterior probabilities used to determine if a cell is KO 
#'   (>0.5) or NP}
#' }
#' 
#' @export
#' 
RunMixscape <- function( 
  object = NULL,
  assay = "PRTB",
  slot = "scale.data",
  labels = "gene",
  nt.class.name = "NT",
  new.class.name = "mixscape_class",
  min.de.genes = 5,
  de.assay = "RNA",
  logfc.threshold = 0.25,
  iter.num = 10,
  verbose = TRUE
) {
  mixtools.installed <- PackageCheck("mixtools", error = FALSE)
  if (!mixtools.installed[1]) {
    stop(
      "Please install the mixtools package to use RunMixscape",
      "\nThis can be accomplished with the following command: ",
      "\n----------------------------------------",
      "\ninstall.packages('mixtools')",
      "\n----------------------------------------",
      call. = FALSE
    )
  }
  assay <- assay %||% DefaultAssay(object = object)
  if (is.null(x = labels)) {
    stop("Please specify target gene class metadata name")
  }
  # de marker genes
  prtb_markers <- c()
  # new metadata column for mixscape classification
  object[[new.class.name]] <- object[[labels]]
  object[[new.class.name]][, 1] <- as.character(x = object[[new.class.name]][, 1])
  genes <- setdiff(
    x = unique(object[[ labels]][, 1]), 
    y = nt.class.name
  )
  # perturbation vectors storage, make list to store probabilities.
  p_ko <- list()
  for (gene in genes) {
    if (verbose) {
      message("Processing ", gene)
    }
    Idents(object = object) <- labels
    # Get object containing only guide of interest + non-targeting
    object.gene <- subset(x = object, idents = c(gene, nt.class.name))
    orig.guide.cells <- WhichCells(object = object.gene, idents = gene)
    DefaultAssay(object = object.gene) <- assay
    # find de genes between guide positive and non-targeting
    de.genes <- TopDEGenesMixscape(
      object.gene, 
      ident.1 = gene, 
      de.assay = de.assay, 
      logfc.threshold = logfc.threshold, 
      labels = labels
    )
    prtb_markers[[gene]] <- de.genes
    # if fewer than 5 DE genes, call all guide cells NP 
    if (length(x = de.genes) < min.de.genes) {
      if (verbose) {
        message("Fewer than ",min.de.genes, " DE genes for ", gene, ". Assigning cells as NP.")
      }
      object.gene[[new.class.name]][orig.guide.cells, 1] <- paste(gene, " NP", sep = "")
    } else {
      object.gene <- ScaleData(object = object.gene, features = de.genes, verbose = FALSE) 
      dat <- GetAssayData(object = object.gene[[assay]], slot = slot)[de.genes, ]
      converged <- FALSE
      n.iter <- 0
      old.classes <- object.gene[[new.class.name]]
      while (! converged & n.iter < iter.num) {
        # Define perturbation vector using only the de genes
        Idents(object = object.gene) <- new.class.name
        nt.cells <- WhichCells(object = object.gene, idents = nt.class.name)
        guide.cells <- WhichCells(object = object.gene, idents = gene)
        vec <- rowMeans(x = dat[, guide.cells]) - rowMeans(x = dat[, nt.cells])
        # project cells onto new perturbation vector
        pvec <- apply(X = dat, MARGIN = 2, FUN = ProjectVec, v2 = vec)
        # define normal distributions mu & sd for guide and nt groups
        guide.norm <- DefineNormalMixscape(pvec[guide.cells])
        nt.norm <- DefineNormalMixscape(pvec[nt.cells])
        # construct mixture model using normalmixEM from mixtools
        mm <- mixtools::normalmixEM(
          x = pvec, 
          mu = c(nt.norm$mu, guide.norm$mu), 
          sigma = c(nt.norm$sd, guide.norm$sd),
          k = 2, 
          mean.constr = c(nt.norm$mu, NA),
          sd.constr = c(nt.norm$sd, NA), 
          verb = FALSE,
          maxit = 5000,
          maxrestarts = 100
        )
        # compute posterior prob 
        lik.ratio <- dnorm(x = pvec[orig.guide.cells], mean = mm$mu[1], sd = mm$sigma[1]) /
          dnorm(x = pvec[orig.guide.cells], mean = mm$mu[2], sd = mm$sigma[2])
        post.prob <- 1 / (1 + lik.ratio)
        # update classifications
        object.gene[[new.class.name]][names(x = which(post.prob > 0.5)), 1] <- gene
        object.gene[[new.class.name]][names(x = which(post.prob < 0.5)), 1] <- paste(gene, " NP", sep = "")
        if (length(x = which(x = object.gene[[new.class.name]] == gene)) < min.de.genes ){
          if (verbose) {
            message("Fewer than ", min.de.genes, " cells assigned as ", gene, "Assigning all to NP.")
          }
          object.gene[[new.class.name]][, 1] <- "NP"
          converged <- TRUE
        }
        if (all(object.gene[[new.class.name]] == old.classes)) {
          converged <- TRUE
        }
        old.classes <- object.gene[[new.class.name]]
        n.iter <- n.iter + 1
      }
      object.gene[[new.class.name]][object.gene[[new.class.name]] == gene, 1] <- paste(gene, " KO", sep = "")
    }
    # assign classifications back to original object
    object[[new.class.name]][Cells(object.gene), 1] <- object.gene[[new.class.name]]
    # add global classifications of KO, NP and NT class
    object[[paste(new.class.name, ".global", sep = "")]] <- as.character(x = sapply(
      X = as.character(x = object[[new.class.name]][, 1]), 
      FUN = function(x) {
        strsplit(x = x, split = ' (?=[^ ]+$)', perl = TRUE)[[1]][2]
      }
    ))
    object[[paste(new.class.name, ".global", sep = "")]][which(x = is.na(x = object[[paste(new.class.name, ".global", sep = "")]])), 1] <- nt.class.name
    p_ko[[gene]] <- post.prob 
  }
  # add posterior probabilities to Seurat object as a meta data column.
  names(x = p_ko) <- NULL
  prob <- unlist(x = p_ko)
  object <- AddMetaData(object = object, metadata = prob, col.name = "p_ko")
  object$p_ko[names(x = which(x = is.na(x = object$p_ko)))] <- 0
  return(object)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Internal
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Function to define Normal distribution - 
# returns list with mu (mean) and sd (standard deviation)
DefineNormalMixscape <- function(x) {
  mu <- mean(x)
  sd <- sd(x)
  return(list(mu = mu, sd = sd))
}

# Helper function to project cells onto the perturbation vector
# @param v1 vector 1
# @param v2 vector 2
#
ProjectVec <- function(v1, v2) {
  return(as.vector(x = (v1 %*% v2) / (v2 %*% v2)))
}

# Function to find top DE genes that pass some p value cutoff between cells 
# with targeting and non-targeting gRNAs.
# 
# @param object An object of class Seurat.
# @param ident.1 Target gene class to find DE genes for.
# @param labels metadata column with target gene classification.
# @param de.assay Name of Assay DE is performed on.
# @param test.use 	Denotes which test to use. See all available tests on 
# FindMarkers documentation.
# @param pval.cut.off P-value cut-off for selection of significantly DE genes.
# @param logfc.threshold Limit testing to genes which show, on average, at 
# least X-fold difference (log-scale) between the two groups of cells. Default 
# is 0.25 Increasing logfc.threshold speeds up the function, but can miss
# weaker signals.
# @param verbose Display messages
# @return 
#
TopDEGenesMixscape <- function(
  object, 
  ident.1, 
  labels = 'gene', 
  de.assay = "RNA", 
  test.use = "LR", 
  pval.cutoff = 5e-2,
  logfc.threshold = 0.25,
  verbose = TRUE
) {
  if (verbose) {
    message("Finding new perturbation gene set")
  }
  de.genes <- data.frame()  
  tryCatch(
    expr = {
      de.genes <- FindMarkers(
        object = object, 
        ident.1 = ident.1, 
        group.by = labels,
        assay = de.assay,
        test.use = test.use,
        logfc.threshold = logfc.threshold, 
        verbose = verbose
      )
      de.genes <- de.genes[de.genes$p_val_adj < pval.cutoff, ]
    }, 
    error = function(e) {}
  )
  return(rownames(x = de.genes))
}
