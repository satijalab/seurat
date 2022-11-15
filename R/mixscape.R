#' @include generics.R
#'
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Calculate a perturbation Signature
#'
#' Function to calculate perturbation signature for pooled CRISPR screen datasets.
#' For each target cell (expressing one target gRNA), we identified 20 cells
#' from the control pool (non-targeting cells) with the most similar mRNA
#' expression profiles. The perturbation signature is calculated by subtracting the
#' averaged mRNA expression profile of the non-targeting neighbors from the mRNA
#' expression profile of the target cell.
#'
#' @param object An object of class Seurat.
#' @param assay Name of Assay PRTB  signature is being calculated on.
#' @param features Features to compute PRTB signature for. Defaults to the
#' variable features set in the assay specified.
#' @param slot Data slot to use for PRTB signature calculation.
#' @param gd.class Metadata column containing target gene classification.
#' @param nt.cell.class Non-targeting gRNA cell classification identity.
#' @param split.by Provide metadata column if multiple biological replicates
#' exist to calculate PRTB signature for every replicate separately.
#' @param num.neighbors Number of nearest neighbors to consider.
#' @param ndims Number of dimensions to use from dimensionality reduction method.
#' @param reduction Reduction method used to calculate nearest neighbors.
#' @param new.assay.name Name for the new assay.
#' @param verbose Display progress + messages
#' @return Returns a Seurat object with a new assay added containing the
#' perturbation signature for all cells in the data slot.
#'
#' @importFrom RANN nn2
#' @export
#' @concept mixscape
#'
CalcPerturbSig <- function(
  object,
  assay = NULL,
  features = NULL,
  slot = "data",
  gd.class = "guide_ID",
  nt.cell.class = "NT",
  split.by = NULL,
  num.neighbors = NULL,
  reduction = "pca",
  ndims = 15,
  new.assay.name = "PRTB",
  verbose = TRUE
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
  features <- features %||% VariableFeatures(object = object[[assay]])
  if (length(x = features) == 0) {
    features <- rownames(x = GetAssayData(object = object[[assay]], slot = slot))
  }
  if (! is.null(x = split.by)) {
    Idents(object = object) <-  split.by
  } else {
    Idents(object = object) <- "rep1"
  }
  replicate <- unique(x = Idents(object = object))
  all_diff <- list()
  all_nt_cells <- Cells(x = object)[which(x = object[[]][gd.class] == nt.cell.class)]
  all_neighbors <- list()
  for (r in replicate) {
    if (verbose) {
      message("Processing ", r)
    }
    all_cells <- WhichCells(object = object, idents = r)
    nt_cells <- intersect(x = all_nt_cells, all_cells)
    # get pca cell embeddings
    all_mtx <- Embeddings(object = object, reduction = reduction)[all_cells, ]
    nt_mtx <- Embeddings(object = object, reduction = reduction)[nt_cells, ]
    # run nn2 to find the 20 nearest NT neighbors for all cells. Use the same
    # number of PCs as the ones you used for umap
    neighbors <- NNHelper(
      data = nt_mtx[, 1:ndims],
      query = all_mtx[, 1:ndims],
      k = num.neighbors,
      method = "rann"
    )
    diff <- PerturbDiff(
      object = object,
      assay = assay,
      slot = slot,
      all_cells = all_cells,
      nt_cells = nt_cells,
      features = features,
      neighbors = neighbors,
      verbose = verbose
    )
    all_diff[[r]] <- diff
    all_neighbors[[make.names(names = paste0(new.assay.name, "_", r))]] <- neighbors
  }
  slot(object = object, name = "tools")[[paste("CalcPerturbSig", assay, reduction, sep = ".")]] <- all_neighbors
  all_diff <- do.call(what = cbind, args = all_diff)
  prtb.assay <- suppressWarnings(
    expr = CreateAssayObject(
      data =  all_diff[, colnames(x = object)],
      min.cells = -Inf,
      min.features = -Inf,
      check.matrix = FALSE
    )
  )
  object[[new.assay.name]] <- prtb.assay
  object <- LogSeuratCommand(object = object)
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
#' @importFrom ggplot2 ggplot geom_bar geom_density coord_flip scale_fill_manual
#' ylab ggtitle theme_classic theme element_text
#' @importFrom patchwork wrap_plots
#'
#' @export
#' @concept mixscape

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

  if(nrow(pos.markers) == 0){
    message("No positive markers pass the logfc.thershold")
    pos.er <- c()
  }

  else{
  pos.markers.list <- rownames(x = pos.markers)[1:min(max.genes, nrow(x = pos.markers))]
  pos.er <- enrichR::enrichr(genes = pos.markers.list, databases = enrich.database)
  pos.er <- do.call(what = cbind, args = pos.er)
  pos.er$log10pval <- -log10(x = pos.er[, paste(enrich.database, sep = ".", "P.value")])
  pos.er$term <- pos.er[, paste(enrich.database, sep = ".", "Term")]
  pos.er <- pos.er[1:num.pathway, ]
  pos.er$term <- factor(x = pos.er$term, levels = pos.er$term[order(pos.er$log10pval)])
  gene.list <- list(pos = pos.er)
  }

  if (isTRUE(x = balanced)) {
    neg.markers <- all.markers[all.markers[, 2] < logfc.threshold & all.markers[, 1] < p.val.cutoff, , drop = FALSE]
    neg.markers.list <- rownames(x = neg.markers)[1:min(max.genes, nrow(x = neg.markers))]
    neg.er <- enrichR::enrichr(genes = neg.markers.list, databases = enrich.database)
    neg.er <- do.call(what = cbind, args = neg.er)
    neg.er$log10pval <- -log10(x = neg.er[, paste(enrich.database, sep = ".", "P.value")])
    neg.er$term <- neg.er[, paste(enrich.database, sep = ".", "Term")]
    neg.er <- neg.er[1:num.pathway, ]
    neg.er$term <- factor(x = neg.er$term, levels = neg.er$term[order(neg.er$log10pval)])

      if(isTRUE(length(neg.er$term) == 0) & isTRUE(length(pos.er == 0))){
        stop("No positive or negative marker genes identified")
      }

      else{
        if(isTRUE(length(neg.er$term) == 0)){

        gene.list <- list(pos = pos.er)

        }
        else{
          gene.list <- list(pos = pos.er, neg = neg.er)
        }
      }

  }
  if (return.gene.list) {
    return(gene.list)
  }

  if(nrow(pos.markers) == 0){
    message("No positive markers to plot")

    if (isTRUE(x = balanced)) {

      p2 <- ggplot(data = neg.er, aes_string(x = "term", y = "log10pval")) +
        geom_bar(stat = "identity", fill = "indianred2") +
        coord_flip() + xlab("Pathway") +
        scale_fill_manual(values = cols, drop = FALSE) +
        ylab("-log10(pval)") +
        ggtitle(paste(enrich.database, ident.1, sep = "_", "negative markers")) +
        theme_classic() +
        geom_text(aes_string(label = "term", y = 0),
                  size = 5,
                  color = "black",
                  position = position_dodge(1),
                  hjust = 0)+
        theme(axis.title.y= element_blank(),
              axis.text.y = element_blank(),
              axis.ticks.y = element_blank())
      p <- p2

    }
    else{
      stop("Nothing to plot")
    }
  }

  else {
  p <- ggplot(data = pos.er, aes_string(x = "term", y = "log10pval")) +
    geom_bar(stat = "identity", fill = "dodgerblue") +
    coord_flip() + xlab("Pathway") +
    scale_fill_manual(values = cols, drop = FALSE) +
    ylab("-log10(pval)") +
    ggtitle(paste(enrich.database, ident.1, sep = "_", "positive markers")) +
    theme_classic() +
    geom_text(aes_string(label = "term", y = 0),
              size = 5,
              color = "black",
              position = position_dodge(1),
              hjust = 0)+
    theme(axis.title.y= element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank())
  if (isTRUE(x = balanced)) {

    p2 <- ggplot(data = neg.er, aes_string(x = "term", y = "log10pval")) +
      geom_bar(stat = "identity", fill = "indianred2") +
      coord_flip() + xlab("Pathway") +
      scale_fill_manual(values = cols, drop = FALSE) +
      ylab("-log10(pval)") +
      ggtitle(paste(enrich.database, ident.1, sep = "_", "negative markers")) +
      theme_classic() +
      geom_text(aes_string(label = "term", y = 0),
                size = 5,
                color = "black",
                position = position_dodge(1),
                hjust = 0)+
      theme(axis.title.y= element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank())
    p <- p+p2

  }
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
#' @concept mixscape
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
#' @concept mixscape
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
  DefaultAssay(object = object) <- pc.assay
  all_genes <- list()
  nt.cells <- WhichCells(object = object, idents = nt.label)
  for (g in gene_list) {
    if (verbose) {
      message(g)
    }
    gd.cells <- WhichCells(object = object, idents = g)
    gene_set <- TopDEGenesMixscape(
      object = object,
      ident.1 = gd.cells,
      ident.2 = nt.cells,
      de.assay = de.assay,
      logfc.threshold = logfc.threshold,
      labels = labels,
      verbose = verbose
    )
    if (length(x = gene_set) < (npcs + 1)) {
      all_genes[[g]] <- character()
      next
    }
    all_genes[[g]] <- gene_set
  }
  all_markers <- unique(x = unlist(x = all_genes))
  missing_genes <- all_markers[!all_markers %in% rownames(x = object[[pc.assay]])]
  object <- GetMissingPerturb(object = object, assay = pc.assay, features = missing_genes, verbose = verbose)
  for (g in gene_list) {
    if (verbose) {
      message(g)
    }
    gene_subset <- subset(x = object, idents = c(g, nt.label))
    gene_set <- all_genes[[g]]
    if (length(x = gene_set) == 0) {
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
      dims = 1:npcs,
      verbose = FALSE
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
#' @concept mixscape
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
#' @concept mixscape
#' @export
#' @method RunLDA Assay
#'
RunLDA.Assay <- function(
  object,
  assay = NULL,
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
#' @concept mixscape
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
  reduction.data <- RunLDA(
    object = object[[assay]],
    assay = assay,
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

#' Run Mixscape
#'
#' Function to identify perturbed and non-perturbed gRNA expressing cells that
#' accounts for multiple treatments/conditions/chemical perturbations.
#'
#' @inheritParams FindMarkers
#' @importFrom ggplot2 geom_density position_dodge
#' @param object An object of class Seurat.
#' @param assay Assay to use for mixscape classification.
#' @param slot Assay data slot to use.
#' @param labels metadata column with target gene labels.
#' @param nt.class.name Classification name of non-targeting gRNA cells.
#' @param new.class.name Name of mixscape classification to be stored in
#' metadata.
#' @param min.de.genes Required number of genes that are differentially
#' expressed for method to separate perturbed and non-perturbed cells.
#' @param min.cells Minimum number of cells in target gene class. If fewer than
#' this many cells are assigned to a target gene class during classification,
#' all are assigned NP.
#' @param de.assay Assay to use when performing differential expression analysis.
#' Usually RNA.
#' @param iter.num Number of normalmixEM iterations to run if convergence does
#' not occur.
#' @param verbose Display messages
#' @param split.by metadata column with experimental condition/cell type
#' classification information. This is meant to be used to account for cases a
#' perturbation is condition/cell type -specific.
#' @param fine.mode When this is equal to TRUE, DE genes for each target gene
#' class will be calculated for each gRNA separately and pooled into one DE list
#' for calculating the perturbation score of every cell and their subsequent
#' classification.
#' @param fine.mode.labels metadata column with gRNA ID labels.
#' @param prtb.type specify type of CRISPR perturbation expected for labeling mixscape classifications. Default is KO.
#' @return Returns Seurat object with with the following information in the
#' meta data and tools slots:
#' \describe{
#'   \item{mixscape_class}{Classification result with cells being either
#'   classified as perturbed (KO, by default) or non-perturbed (NP) based on their target
#'   gene class.}
#'   \item{mixscape_class.global}{Global classification result (perturbed, NP or NT)}
#'   \item{p_ko}{Posterior probabilities used to determine if a cell is KO (default). Name of this item will change to match prtb.type parameter setting.
#'   (>0.5) or NP}
#'   \item{perturbation score}{Perturbation scores for every cell calculated in
#'   the first iteration of the function.}
#' }
#'
#' @export
#' @concept mixscape
#'
RunMixscape <- function(
  object,
  assay = "PRTB",
  slot = "scale.data",
  labels = "gene",
  nt.class.name = "NT",
  new.class.name = "mixscape_class",
  min.de.genes = 5,
  min.cells = 5,
  de.assay = "RNA",
  logfc.threshold = 0.25,
  iter.num = 10,
  verbose = FALSE,
  split.by = NULL,
  fine.mode = FALSE,
  fine.mode.labels = "guide_ID",
  prtb.type = "KO"
) {
  mixtools.installed <- PackageCheck("mixtools", error = FALSE)
  if (!mixtools.installed[1]) {
    stop("Please install the mixtools package to use RunMixscape",
         "\nThis can be accomplished with the following command: ",
         "\n----------------------------------------",
         "\ninstall.packages('mixtools')",
         "\n----------------------------------------", call. = FALSE)
  }
  assay <- assay %||% DefaultAssay(object = object)
  if (is.null(x = labels)) {
    stop("Please specify target gene class metadata name")
  }
  prtb_markers <- list()
  object[[new.class.name]] <- object[[labels]]
  object[[new.class.name]][, 1] <- as.character(x = object[[new.class.name]][, 1])
  object[[paste0(new.class.name, "_p_", tolower(x = prtb.type))]] <- 0
  #create list to store perturbation scores.
  gv.list <- list()

  if (is.null(x = split.by)) {
    split.by <- splits <- "con1"
  } else {
    splits <- as.character(x = unique(x = object[[split.by]][, 1]))
  }

  # determine gene sets across all splits/groups
  cells.s.list <- list()
  for (s in splits) {
    Idents(object = object) <- split.by
    cells.s <- WhichCells(object = object, idents = s)
    cells.s.list[[s]] <- cells.s
    genes <- setdiff(x = unique(x = object[[labels]][cells.s, 1]), y = nt.class.name)
    Idents(object = object) <- labels
    for (gene in genes) {
      if (isTRUE(x = verbose)) {
        message("Processing ", gene)
      }
      orig.guide.cells <- intersect(x = WhichCells(object = object, idents = gene), y = cells.s)
      nt.cells <- intersect(x = WhichCells(object = object, idents = nt.class.name), y = cells.s)
      if (isTRUE(x = fine.mode)) {
        guides <- setdiff(x = unique(x = object[[fine.mode.labels]][orig.guide.cells, 1]), y = nt.class.name)
        all.de.genes <- c()
        for (gd in guides) {
          gd.cells <- rownames(x = object[[]][orig.guide.cells, ])[which(x = object[[]][orig.guide.cells, fine.mode.labels] == gd)]
          de.genes <- TopDEGenesMixscape(
            object = object,
            ident.1 = gd.cells,
            ident.2 = nt.cells,
            de.assay = de.assay,
            logfc.threshold = logfc.threshold,
            labels = fine.mode.labels,
            verbose = verbose
          )
          all.de.genes <- c(all.de.genes, de.genes)

        }
        all.de.genes <- unique(all.de.genes)
      } else {
        all.de.genes <- TopDEGenesMixscape(
          object = object,
          ident.1 = orig.guide.cells,
          ident.2 = nt.cells,
          de.assay = de.assay,
          logfc.threshold = logfc.threshold,
          labels = labels,
          verbose = verbose
        )
      }
      prtb_markers[[s]][[gene]] <- all.de.genes
      if (length(x = all.de.genes) < min.de.genes) {
        prtb_markers[[s]][[gene]] <- character()
      }
    }
  }
  all_markers <- unique(x = unlist(x = prtb_markers))
  missing_genes <- all_markers[!all_markers %in% rownames(x = object[[assay]])]
  object <- GetMissingPerturb(object = object, assay = assay, features = missing_genes, verbose = verbose)
  for (s in splits) {
    cells.s <- cells.s.list[[s]]
    genes <- setdiff(x = unique(x = object[[labels]][cells.s, 1]), y = nt.class.name)
    if (verbose) {
      message("Classifying cells for: ")
    }
    for (gene in genes) {
      Idents(object = object) <- labels
      post.prob <- 0
      orig.guide.cells <- intersect(x = WhichCells(object = object, idents = gene), y = cells.s)
      nt.cells <- intersect(x = WhichCells(object = object, idents = nt.class.name), y = cells.s)
      all.cells <- c(orig.guide.cells, nt.cells)
      if (length(x = prtb_markers[[s]][[gene]]) == 0) {
        if (verbose) {
          message("  Fewer than ", min.de.genes, " DE genes for ", gene,
                  ". Assigning cells as NP.")
        }
        object[[new.class.name]][orig.guide.cells, 1] <- paste0(gene, " NP")
      } else {
        if (verbose) {
          message("  ", gene)
        }
        de.genes <- prtb_markers[[s]][[gene]]
        dat <- GetAssayData(object = object[[assay]], slot = "data")[de.genes, all.cells, drop = FALSE]
        if (slot == "scale.data") {
          dat <- ScaleData(object = dat, features = de.genes, verbose = FALSE)
        }
        converged <- FALSE
        n.iter <- 0
        old.classes <- object[[new.class.name]][all.cells, ]
        while (!converged && n.iter < iter.num) {
          Idents(object = object) <- new.class.name
          guide.cells <- intersect(x = WhichCells(object = object, idents = gene), y = cells.s)
          vec <- rowMeans2(x = dat[, guide.cells, drop = FALSE]) - rowMeans2(x = dat[, nt.cells, drop = FALSE])
          pvec <- apply(X = dat, MARGIN = 2, FUN = ProjectVec, v2 = vec)
          if (n.iter == 0){
            #store pvec
            gv <- as.data.frame(x = pvec)
            gv[, labels] <- nt.class.name
            gv[intersect(x = rownames(x = gv), y = guide.cells), labels] <- gene
            gv.list[[gene]][[s]] <- gv
          }
          guide.norm <- DefineNormalMixscape(pvec[guide.cells])
          nt.norm <- DefineNormalMixscape(pvec[nt.cells])
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
          lik.ratio <- dnorm(x = pvec[orig.guide.cells], mean = mm$mu[1], sd = mm$sigma[1]) /
            dnorm(x = pvec[orig.guide.cells], mean = mm$mu[2], sd = mm$sigma[2])
          post.prob <- 1/(1 + lik.ratio)
          object[[new.class.name]][names(x = which(post.prob > 0.5)), 1] <- gene
          object[[new.class.name]][names(x = which(post.prob < 0.5)), 1] <- paste(gene, " NP", sep = "")
          if (length(x = which(x = object[[new.class.name]] == gene & Cells(x = object) %in% cells.s)) < min.de.genes) {
            if (verbose) {
              message("Fewer than ", min.cells, " cells assigned as ",
                      gene, "Assigning all to NP.")
            }
            object[[new.class.name]][guide.cells, 1] <- "NP"
            converged <- TRUE
          }
          if (all(object[[new.class.name]][all.cells, ] == old.classes)) {
            converged <- TRUE
          }
          old.classes <- object[[new.class.name]][all.cells, ]
          n.iter <- n.iter + 1
        }
        object[[new.class.name]][which(x = object[[new.class.name]] == gene & Cells(x = object) %in% cells.s), 1] <- paste(gene, prtb.type, sep = " ")
      }
      object[[paste0(new.class.name, ".global")]] <- as.character(x = sapply(X = as.character(x = object[[new.class.name]][, 1]), FUN = function(x) {strsplit(x = x, split = " (?=[^ ]+$)", perl = TRUE)[[1]][2]}))
      object[[paste0(new.class.name, ".global")]][which(x = is.na(x = object[[paste0(new.class.name, ".global")]])), 1] <- nt.class.name
      object[[paste0(new.class.name,"_p_", tolower(prtb.type))]][names(x = post.prob), 1] <- post.prob
    }
  }
  Tool(object = object) <- gv.list
  Idents(object = object) <- new.class.name
  return(object)
}

#' Differential expression heatmap for mixscape
#'
#' Draws a heatmap of single cell feature expression with cells ordered by their
#' mixscape ko probabilities.
#'
#' @inheritParams FindMarkers
#' @inheritParams DoHeatmap
#' @param max.cells.group Number of cells per identity to plot.
#' @param max.genes Total number of DE genes to plot.
#' @param balanced Plot an equal number of genes with both groups of cells.
#' @param order.by.prob Order cells on heatmap based on their mixscape knockout
#' probability from highest to lowest score.
#' @param group.by (Deprecated) Option to split densities based on mixscape
#' classification. Please use mixscape.class instead
#' @param mixscape.class metadata column with mixscape classifications.
#' @param prtb.type specify type of CRISPR perturbation expected for labeling
#' mixscape classifications. Default is KO.
#' @param fc.name Name of the fold change, average difference, or custom
#' function column in the output data.frame. Default is avg_log2FC
#' @param pval.cutoff P-value cut-off for selection of significantly DE genes.
#' @return A ggplot object.
#'
#' @importFrom stats median
#' @importFrom scales hue_pal
#' @importFrom ggplot2 annotation_raster coord_cartesian ggplot_build aes_string
#' @export
#' @concept mixscape
#'
MixscapeHeatmap <- function(
  object,
  ident.1 = NULL,
  ident.2 = NULL,
  balanced = TRUE,
  logfc.threshold = 0.25,
  assay = "RNA",
  max.genes = 100,
  test.use ='wilcox',
  max.cells.group = NULL,
  order.by.prob = TRUE,
  group.by = NULL,
  mixscape.class = "mixscape_class",
  prtb.type = "KO",
  fc.name = "avg_log2FC",
  pval.cutoff = 5e-2,
  ...
)
{
  if (!is.null(x = group.by)) {
    message("The group.by parameter is being deprecated. Please use ",
            "mixscape.class instead. Setting mixscape.class = ", group.by,
            " and continuing.")
    mixscape.class <- group.by
  }
  DefaultAssay(object = object) <- assay
  if (is.numeric(x = max.genes)) {
    all.markers <- FindMarkers(
      object = object,
      ident.1 = ident.1,
      ident.2 = ident.2,
      only.pos = FALSE,
      logfc.threshold = logfc.threshold,
      test.use = test.use
    )
    if (balanced) {
      pos.markers <- all.markers[which(x = all.markers[,fc.name] > (logfc.threshold)), ]
      neg.markers <- all.markers[which(x = all.markers[,fc.name] < (-logfc.threshold)), ]
      if (length(x = rownames(x = subset(x = pos.markers, p_val < pval.cutoff))) < max.genes ) {
        marker.list <- c(rownames(x = subset(x = pos.markers, p_val < pval.cutoff)))
        if (length(x = rownames(x = subset(x = neg.markers, p_val < pval.cutoff))) < max.genes){
          marker.list <- c(marker.list, rownames(x = subset(x = neg.markers, p_val < pval.cutoff)))
        } else {
          marker.list <- c(marker.list, rownames(x = subset(x = neg.markers, p_val < pval.cutoff))[1:max.genes])
        }
      } else {
        marker.list <- c(rownames(x = subset(x = pos.markers, p_val < pval.cutoff))[1:max.genes])
        if (length(x = rownames(x = subset(x = neg.markers, p_val < pval.cutoff))) < max.genes) {
          marker.list <- c(marker.list, rownames(x = subset(x = neg.markers, p_val < pval.cutoff)))
        } else {
          marker.list <- c(marker.list, rownames(x = subset(x = neg.markers, p_val < pval.cutoff))[1:max.genes])
        }
      }
    }
    else {
      pos.markers <- all.markers[which(x = all.markers[, fc.name] > (logfc.threshold)),]
      if (length(x = rownames(x = subset(x = pos.markers, p_val < pval.cutoff))) < max.genes ){
        marker.list <- c(rownames(x = subset(x = pos.markers, p_val < pval.cutoff)))
      } else {
        marker.list <- c(rownames(x = subset(x = pos.markers, p_val < pval.cutoff))[1:max.genes])
      }
    }
    if (is.null(x = max.cells.group)) {
      if (is.null(x = group.by)) {
        sub2 <- subset(x = object, idents = c(ident.1, ident.2))
      } else{
        sub2 <- subset(x = object, idents = c(ident.1, ident.2))
        Idents(object = sub2) <- group.by
      }
    }
    else {
      if (is.null(x = group.by)) {
        sub2 <- subset(x = object, idents = c(ident.1, ident.2), downsample = max.cells.group)
      } else {
        sub <- subset(x = object, idents = c(ident.1, ident.2))
        Idents(object = sub) <- group.by
        sub2 <- subset(x = sub, downsample = max.cells.group)
      }
    }
    sub2 <- ScaleData(object = sub2, features = marker.list, assay = assay)
    if (isTRUE(x = order.by.prob)) {
      p_ko <- sub2[[paste0(mixscape.class, "_p_", tolower(x = prtb.type) )]][, 1, drop = FALSE]
      ordered.cells <- rownames(x = p_ko)[order(p_ko[,1], decreasing = TRUE)]
      p <- DoHeatmap(object = sub2, features = marker.list, label = TRUE, cells = ordered.cells, assay = assay, ...)
    } else{
      p <- DoHeatmap(object = sub2, features = marker.list, label = TRUE, cells = sample(x = Cells(x = sub2)), assay = assay, ...)
    }
    return(p)
  }
}


#' Function to plot perturbation score distributions.
#'
#' Density plots to visualize perturbation scores calculated from RunMixscape
#' function.
#'
#' @param object An object of class Seurat.
#' @param target.gene.ident Target gene name to visualize perturbation scores for.
#' @param target.gene.class meta data column specifying all target gene names in the experiment.
#' @param before.mixscape Option to split densities based on mixscape classification (default) or original target gene classification.
#' Default is set to NULL and plots cells by original class ID.
#' @param col Specify color of target gene class or knockout cell class. For
#' control non-targeting and non-perturbed cells, colors are set to different
#' shades of grey.
#' @param mixscape.class meta data column specifying mixscape classifications.
#' @param prtb.type specify type of CRISPR perturbation expected for labeling mixscape classifications. Default is KO.
#' @param split.by For datasets with more than one cell type. Set equal TRUE to visualize perturbation scores for each cell type separately.
#' @return A ggplot object.
#'
#' @importFrom stats median
#' @importFrom scales hue_pal
#' @importFrom ggplot2 annotation_raster coord_cartesian ggplot_build aes_string
#' geom_density theme_classic
#' @export
#' @concept mixscape
#'
PlotPerturbScore <- function(
  object,
  target.gene.class = "gene",
  target.gene.ident = NULL,
  mixscape.class = "mixscape_class",
  col = "orange2",
  split.by = NULL,
  before.mixscape = FALSE,
  prtb.type = "KO"
){

  if(is.null(target.gene.ident) == TRUE){
    message("Please provide name of target gene class to plot")
  }
  prtb_score_list <- Tool(object = object, slot = "RunMixscape")[[target.gene.ident]]

  for (nm in names(prtb_score_list)){
    prtb_score_list[[nm]]['name'] <- nm
  }
  prtb_score <- do.call(rbind, prtb_score_list)
  prtb_score[, 2] <- as.factor(x = prtb_score[, 2])
  gd <- setdiff(x = unique(x = prtb_score[, target.gene.class]), y = target.gene.ident)
  colnames(x = prtb_score)[2] <- "gene"
  prtb_score$cell.bc <- sapply(rownames(prtb_score), FUN = function(x) substring(x, regexpr("[.]", x) + 1))

  if (isTRUE(x = before.mixscape)) {
    cols <- setNames(
      object = c("grey49", col),
      nm = c(gd, target.gene.ident)
    )

    p <- ggplot(data = prtb_score, mapping = aes_string(x = "pvec", color = "gene")) +
      geom_density() + theme_classic()
    top_r <- ggplot_build(p)$layout$panel_params[[1]]$y.range[2]
    prtb_score$y.jitter <- prtb_score$pvec
    prtb_score$y.jitter[prtb_score[, "gene"] == gd] <- runif(
      n = prtb_score$y.jitter[prtb_score[, "gene"] == gd],
      min = 0.001,
      max = top_r / 10
    )
    prtb_score$y.jitter[prtb_score[,"gene"] == target.gene.ident] <- runif(
      n = prtb_score$y.jitter[prtb_score[, "gene"] == target.gene.ident],
      min = -top_r / 10,
      max = 0
    )

    if(is.null(split.by)==FALSE) {
      prtb_score$split <- as.character(object[[split.by]][prtb_score$cell.bc,1])
      p2 <- p + scale_color_manual(values = cols, drop = FALSE) +
        geom_density(size = 1.5) +
        geom_point(data = prtb_score, aes_string(x = "pvec", y = "y.jitter"), size = 0.1) +
        theme(axis.text = element_text(size = 18), axis.title = element_text(size = 20)) +
        ylab("Cell density") + xlab("perturbation score") +
        theme(legend.key.size = unit(1, "cm"),
              legend.text = element_text(colour = "black", size = 14),
              legend.title = element_blank(), plot.title = element_text(size = 16, face = "bold"))+
        facet_wrap(vars(split))
    }

    else{
      p2 <- p + scale_color_manual(values = cols, drop = FALSE) +
        geom_density(size = 1.5) +
        geom_point(data = prtb_score, aes_string(x = "pvec", y = "y.jitter"), size = 0.1) +
        theme(axis.text = element_text(size = 18), axis.title = element_text(size = 20)) +
        ylab("Cell density") + xlab("perturbation score") +
        theme(legend.key.size = unit(1, "cm"),
              legend.text = element_text(colour = "black", size = 14),
              legend.title = element_blank(), plot.title = element_text(size = 16, face = "bold"))
    }
  }


  else {
    cols <- setNames(
      object = c("grey49", "grey79", col),
      nm = c(gd, paste0(target.gene.ident, " NP"), paste(target.gene.ident, prtb.type, sep = " "))
    )
    #add mixscape identities
    prtb_score$mix <- object[[mixscape.class]][prtb_score$cell.bc,]

    p <- ggplot(data = prtb_score, aes_string(x = "pvec", color = "mix")) +
      geom_density() + theme_classic()

    top_r <- ggplot_build(p)$layout$panel_params[[1]]$y.range[2]
    prtb_score$y.jitter <- prtb_score$pvec
    gd2 <- setdiff(
      x = unique(x = prtb_score[, "mix"]),
      y = c(paste0(target.gene.ident, " NP"), paste(target.gene.ident, prtb.type, sep = " "))
    )
    prtb_score$y.jitter[prtb_score[, "mix"] == gd2] <- runif(
      n = prtb_score$y.jitter[prtb_score[, "mix"] == gd2],
      min = 0.001,
      max = top_r / 10
    )
    prtb_score$y.jitter[prtb_score$mix == paste(target.gene.ident, prtb.type, sep = " ")] <- runif(
      n = prtb_score$y.jitter[prtb_score[, "mix"] == paste(target.gene.ident, prtb.type, sep = " ")],
      min = -top_r / 10,
      max = 0
    )
    prtb_score$y.jitter[prtb_score$mix == paste0(target.gene.ident, " NP")] <- runif(
      n = prtb_score$y.jitter[prtb_score[, "mix"] == paste0(target.gene.ident, " NP")],
      min = -top_r / 10,
      max = 0
    )
    prtb_score[, "mix"] <- as.factor(x = prtb_score[,"mix"])

    if(is.null(split.by) == FALSE){
      prtb_score$split <- as.character(object[[split.by]][prtb_score$cell.bc,1])
      p2 <- ggplot(data = prtb_score, aes_string(x = "pvec", color = "mix")) +
        scale_color_manual(values = cols, drop = FALSE) +
        geom_density(size = 1.5) +
        geom_point(aes_string(x = "pvec", y = "y.jitter"), size = 0.1) +
        theme_classic() +
        theme(axis.text = element_text(size = 18), axis.title = element_text(size = 20)) +
        ylab("Cell density") + xlab("perturbation score") +
        theme(legend.key.size = unit(1, "cm"),
              legend.text = element_text(colour ="black", size = 14),
              legend.title = element_blank(),
              plot.title = element_text(size = 16, face = "bold"))+
        facet_wrap(vars(split))
    }
    else{
      p2 <- p + scale_color_manual(values = cols, drop = FALSE) +
        geom_density(size = 1.5) +
        geom_point(data = prtb_score, aes_string(x = "pvec", y = "y.jitter"), size = 0.1) +
        theme(axis.text = element_text(size = 18), axis.title = element_text(size = 20)) +
        ylab("Cell density") + xlab("perturbation score") +
        theme(legend.key.size = unit(1, "cm"),
              legend.text = element_text(colour ="black", size = 14),
              legend.title = element_blank(),
              plot.title = element_text(size = 16, face = "bold"))
    }

  }
  return(p2)
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

# Get missing perturbation signature for missing features
#
# @param object Seurat object
# @param assay Perturbation signature assay name
# @param features vector of features to compute for
# @param verbose display progress
# @return Returns Seurat object with assay updated with new features
#
GetMissingPerturb <- function(object, assay, features, verbose = TRUE) {
  if (length(x = features) == 0) {
    return(object)
  }
  if (verbose) {
    message("Computing perturbation signature for missing features.")
  }
  command <- grep(pattern = "CalcPerturbSig", x = Command(object = object), value = TRUE)
  command.match <- sapply(X = command, FUN = function(x) {
    Command(object = object, command = x, value = "new.assay.name") == assay
  })
  if (length(x = which(x = command.match)) > 1) {
    stop("Ambiguous command log.")
  }
  if(length(x = which(x = command.match)) == 0) {
    stop("Cannot find previously run CalcPertubSig command. Please make sure you've run CalcPerturbSig to create the provided assay.")
  }
  command <- names(x = command.match)
  if ("split.by" %in% names(x = slot(object = Command(object = object, command = command), name ="params"))) {
    split.by <- Command(object = object, command = command, value = "split.by")
  } else {
    split.by <- NULL
  }
  gd.class <- Command(object = object, command = command, value = "gd.class")
  nt.cell.class <- Command(object = object, command = command, value = "nt.cell.class")
  slot <- Command(object = object, command = command, value = "slot")
  assay.orig <- Command(object = object, command = command, value = "assay")
  old.idents <- Idents(object = object)
  if (! is.null(x = split.by)) {
    Idents(object = object) <-  split.by
  } else {
    Idents(object = object) <- "rep1"
  }
  replicate <- unique(x = Idents(object = object))
  all_diff <- list()
  all_nt_cells <- Cells(x = object)[which(x = object[[]][gd.class] == nt.cell.class)]
  features <- setdiff(x = features, y = rownames(x = object[[assay]]))
  for (r in replicate) {
    # isolate nt cells
    all_cells <- WhichCells(object = object, idents = r)
    nt_cells <- intersect(x = all_nt_cells, all_cells)
    # pull previously computed neighbors
    neighbors <- Tool(object = object, slot = command)[[make.names(names = paste0(assay, "_", r))]]
    diff <- PerturbDiff(
      object = object,
      assay = assay.orig,
      slot = slot,
      all_cells = all_cells,
      nt_cells = nt_cells,
      features = features,
      neighbors = neighbors,
      verbose = verbose
    )
    all_diff[[r]] <- diff
  }
  all_diff <- do.call(what = cbind, args = all_diff)
  all_diff <- all_diff[, colnames(x = object[[assay]]), drop = FALSE]
  new.assay <- CreateAssayObject(
    data = rbind(
      GetAssayData(object = object[[assay]], slot = "data"),
      all_diff
    ),
    min.cells = 0,
    min.features = 0,
    check.matrix = FALSE
  )
  new.assay <- SetAssayData(
    object = new.assay,
    slot = "scale.data",
    new.data = GetAssayData(object = object[[assay]], slot = "scale.data")
  )
  object[[assay]] <- new.assay
  Idents(object = object) <- old.idents
  return(object)
}

# Helper function to compute the perturbation differences - enables reuse in
# GetMissingPerturb
#
# @param object Seurat object
# @param assay assay to use
# @param slot slot to use
# @param all_cells vector of cell names to compute difference for
# @param nt_cells vector of nt cell names
# @param features vector of features to compute for
# @param neighbors Neighbor object containing indices of nearest NT cells
# @param verbose display progress bar
# @return returns matrix of perturbation differences
#
#' @importFrom matrixStats rowMeans2
#' @importFrom Matrix sparseMatrix colSums
#'
PerturbDiff <- function(object, assay, slot, all_cells, nt_cells, features, neighbors, verbose) {
  nt_data <- as.matrix(x = expm1(x = GetAssayData(object = object, assay = assay, slot = slot)[features, nt_cells, drop = FALSE]))
  mysapply <- ifelse(test = verbose, yes = pbsapply, no = sapply)
  # new_expr <- mysapply(X = all_cells, FUN = function(i) {
  #   index <- Indices(object = neighbors)[i, ]
  #   nt_cells20 <- nt_cells[index]
  #   avg_nt <- rowMeans2(x = nt_data[, nt_cells20, drop = FALSE])
  #   avg_nt <- as.matrix(x = avg_nt)
  #   colnames(x = avg_nt) <- i
  #   return(avg_nt)
  # })
  idx <- Indices(object = neighbors)[all_cells,]
  model.matrix <- sparseMatrix(i = as.vector(idx), j = rep(1:nrow(x = idx), times = ncol(x = idx)), x = 1, dims = c(length(x = nt_cells), nrow(x = idx)))
  model.matrix <- model.matrix/rep(colSums(model.matrix), each = nrow(x = model.matrix))
  new_expr <- nt_data %*% model.matrix

  new_expr <- matrix(data = new_expr, nrow = length(x = features))
  new_expr <- log1p(x = new_expr)
  rownames(x = new_expr) <- rownames(x = nt_data)
  colnames(x = new_expr) <- all_cells
  diff <- new_expr - as.matrix(GetAssayData(object = object, slot = slot, assay = assay)[features, colnames(x = new_expr), drop = FALSE])
  return(diff)
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
# @param ident.1 Target gene class or cells to find DE genes for.
# @param ident.2 Non-targetting class or cells
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
  ident.2 = NULL,
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
        ident.2 = ident.2,
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
