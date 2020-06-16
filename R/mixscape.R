#'  Function to calculate perturbation score for pooled CRISPR screen datasets.
#'  @param object An object of class Seurat.
#'  @param assay Name of Assay PRTB  signature is being calculated on.
#'  @param slot Data slot to use for PRTB score calculation.
#'  @param gd.class Metadata column containing target gene classification.
#'  @param nt.cell.class Non-targeting gRNA cell classification identity.
#'  @param split.by Provide metadata column if multiple biological replicates exist to calculate PRTB score for every replicate separately.
#'  @param num.neighbors Number of nearest neighbors to consider.
#'  @param ndims Number of dimensions to use from dimensionality reduction method.
#'  @param reduction Reduction method used to calculate nearest neighbors.
#'  @param new.assay.name Name for the new assay.
#'  @export

CalcPerturbScore <- function ( object, 
                               assay = NULL,
                               slot = "data", 
                               gd.class = "guide_ID",
                               nt.cell.class = "NT",
                               split.by = NULL,
                               num.neighbors = NULL,
                               reduction = "pca", 
                               ndims= 15,
                               new.assay.name = "PRTB",
                               ...
) {
  if (is.null(assay) == TRUE){
    assay <- DefaultAssay(object)
  }
  
  else{
    
    if (is.null(reduction) == TRUE){
      stop('Please provide dimensionality reduction name.')
    } 
    else{
      
      if (is.null(num.neighbors) == TRUE){
        stop("Please specify number of nearest neighbors to consider")
      }
      else{
        if(is.null(ndims) == TRUE){
          stop("Please provide number of ",reduction, " dimensions to consider")
        }
        else{
          
          if(is.null(split.by) == FALSE){
            
            Idents(object) = split.by
            replicate = levels(Idents(object))
          }
          
          else{
            Idents(object) <- "rep1"
            replicate <- levels(Idents(object))
          }
          all_diff = matrix(nrow = length(rownames(GetAssayData(object = object, assay = assay, slot = slot))), ncol = 0)
          
          for (r in replicate) {
            rep1 <- object[, WhichCells(object = object, idents = r)]
            #isolate nt cells
            all_cells <- Cells(rep1)
            nt_cells <- Cells(rep1[,grep(nt.cell.class, rep1@meta.data[,gd.class], value = FALSE)])
            #subset the objects based on guide ID
            all <- rep1[,all_cells]
            nt <- rep1[,nt_cells]
            
            #get pca cell embeddings
            all_mtx <- Embeddings(all, reduction = reduction)
            nt_mtx <- Embeddings(nt, reduction = reduction)
            
            #run nn2 to find the 20 nearest NT neighbors for all cells. Use the same number of PCs as the ones you used for umap
            mtx <- nn2(data = nt_mtx[,1:ndims], query = all_mtx[,1:ndims], k = num.neighbors)
            #browser()
            nt_data <- expm1(GetAssayData(object = nt, assay = assay, slot = slot))
            new_expr <- pbsapply(X = 1:length(all_cells), FUN = function(i) {
              #browser()
              index <- mtx$nn.idx[i,]
              nt_cells20 <- nt_cells[index]
              #Idents(nt_sub) <- "nt_sub"
              avg_nt <- rowMeans(nt_data[,nt_cells20])
              avg_nt <- as.matrix(avg_nt)
              colnames(avg_nt) <- all_cells[i]
              return(avg_nt)
            })
            new_expr <- log1p(new_expr)
            rownames(new_expr) <- rownames(GetAssayData(object = nt, assay = assay, slot = slot))
            colnames(new_expr) <- all_cells
            diff <- new_expr - GetAssayData(object = rep1, slot = slot, assay = assay)[, colnames(new_expr)]
            all_diff <- cbind(all_diff, diff )
          } 
        }
      }
    }
  }
  all_diff <- as(Class = "dgCMatrix", object = all_diff)
  prtb.assay <- CreateAssayObject(data =  all_diff[, colnames(object)], min.cells = -Inf, min.features = -Inf)
  object[[new.assay.name]] <- prtb.assay
  return(object)
}

#' Helper function to project cells onto the perturbation vector
#' @param v1 vector 1
#' @param v2 vector 2
#' @export
#' 
Project <- function(v1, v2) {
  return(as.vector((v1 %*% v2) / (v2 %*% v2)))
}

#' Function to find top DE genes that pass some p value cutoff between cells with targeting and non-targeting gRNAs.
#' @param object An object of class Seurat.
#' @param ident.1 Target gene class to find DE genes for.
#' @param labels metadata column with target gene classification.
#' @param de.assay Name of Assay DE is performed on.
#' @param test.use 	Denotes which test to use. See all available tests on FindMarkers documentation.
#' @param pval.cut.off P-value cut-off for selection of significantly DE genes.
#' @param logfc.threshold Limit testing to genes which show, on average, at least X-fold difference (log-scale) between the two groups of cells. Default is 0.25 Increasing logfc.threshold speeds up the function, but can miss weaker signals.
#' @export
#' 
TopDEGenesMixscape <- function(
  object, 
  ident.1, 
  labels = 'gene', 
  de.assay = "RNA", 
  test.use = "LR", 
  pval.cutoff = 5e-2,
  logfc.threshold = 0.25
) {
  message("Finding new perturbation gene set")
  de.genes <- data.frame()  
  tryCatch(
    expr = {
      de.genes <- FindMarkers(
        object = object, 
        ident.1 = ident.1, 
        group.by = labels,
        assay = de.assay,
        test.use = test.use,
        logfc.threshold = logfc.threshold
      )
      de.genes <- subset(de.genes, p_val_adj < pval.cutoff)
    }, 
    error = function(e) {}
  )
  return(rownames(de.genes))
}


#' Function to define Normal distribution - returns list with mu (mean) and sd (standard deviation)
#' @export
DefineNormalMixscape <- function(x) {
  mu <- mean(x)
  sd <- sd(x)
  return(list(mu = mu, sd = sd))
}

#' Function to identify perturbed and non-perturbed gRNA expressing cells.
#' @param object An object of class Seurat.
#' @param assay Assay to use for mixscape classification.
#' @param slot Assay data slot to use.
#' @param labels metadata column with target gene classifications.
#' @param nt.class.name Classification name of non-targeting gRNA cells.
#' @param new.class.name Name of mixscape classification to be stored in metadata.
#' @param min.de.genes Required number of genes that are differentially expressed for method to separate perturbed and non-perturbed cells. 
#' @param de.assay Assay to use when performing differential expression analysis. Usually RNA.
#' @param iter.num Number of normalmixEM iterations to run if convergence does not occur.
#' @import mixtools
#' @import reshape2
#' @export
#' 
RunMixscape <- function( object = NULL,
                         assay = "PRTB",
                         slot = "scale.data",
                         labels = "gene",
                         nt.class.name = "NT",
                         new.class.name = "mixscape_class",
                         min.de.genes = 5,
                         de.assay = "RNA",
                         logfc.threshold = 0.25,
                         iter.num = 10,
                         ...
                         
){
  if (is.null(assay) == TRUE){
    assay <- DefaultAssay(object)
  }
  
  if (is.null(labels) == TRUE){
    stop("Please specify target gene class metadata name")
  }
  #de marker genes
  prtb_markers <- c()
  #new metadata column for mixscape classification
  object[[new.class.name]] <- object[[labels]]
  object[[new.class.name]][,1] <- as.character(object[[new.class.name]][,1])
  
  genes <- setdiff(unique(object[[ labels]][,1]), y = nt.class.name)
  
  #pertubration vectors storage
  #gv.list <- list(list())
  #make list to store probabilities.
  p_ko <- list()
  
  for (gene in genes){
    message("Processing ", gene)
    Idents(object) <- labels
    
    # Get object containing only guide of interest + non-targeting
    object.gene <- subset(object, idents = c(gene, nt.class.name))
    orig.guide.cells <- WhichCells(object.gene, idents = gene)
    DefaultAssay(object.gene) <- assay
    
    # find de genes between guide positive and non-targeting
    de.genes <- TopDEGenesMixscape(object.gene, ident.1 = gene, de.assay = de.assay, logfc.threshold = logfc.threshold, labels = labels)
    prtb_markers[[gene]] <- de.genes
    
    # if fewer than 5 DE genes, call all guide cells NP 
    if (length(de.genes) < min.de.genes) {
      message("Fewer than ",min.de.genes, " DE genes for ", gene, ". Assigning cells as NP.")
      object.gene[[new.class.name]][orig.guide.cells,1] <- paste(gene, " NP", sep = "")
    } else {
      object.gene <- ScaleData(object.gene, features = de.genes, verbose = FALSE) 
      dat <- GetAssayData(object = object.gene[[assay]], slot = slot)[de.genes, ]
      converged <- FALSE
      n.iter <- 0
      old.classes <- object.gene[[new.class.name]]
      
      while (! converged & n.iter < iter.num) {
        #message("Iteration ", n.iter + 1)
        
        # Define pertubation vector using only the de genes
        Idents(object.gene) <- new.class.name
        nt.cells <- WhichCells(object.gene, idents = nt.class.name)
        guide.cells <- WhichCells(object.gene, idents = gene)
        vec <- rowMeans(dat[, guide.cells]) - rowMeans(dat[, nt.cells])
        
        # project cells onto new perturbation vector
        pvec <- apply(X = dat, MARGIN = 2, FUN = Project, v2 = vec)
        
        #store pvec
        #gv<- melt(pvec)
        #gv$name <- nt.class.name
        #gv[intersect(rownames(gv), guide.cells),"name"] <- gene
        #gv.list[[gene]][[n.iter+1]] <- gv
        
        # define normal distributions mu & sd for guide and nt groups
        guide.norm <- DefineNormalMixscape(pvec[guide.cells])
        nt.norm <- DefineNormalMixscape(pvec[nt.cells])
        
        # construct mixture model using normalmixEM from mixtools
        mm <- normalmixEM(
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
        object.gene[[new.class.name]][names(which(post.prob > 0.5)),1] <- gene
        object.gene[[new.class.name]][names(which(post.prob < 0.5)),1] <- paste(gene, " NP", sep = "")
        
        if (length(which(object.gene[[new.class.name]] == gene)) < min.de.genes ){
          message("Fewer than ", min.de.genes, " cells assigned as ", gene, "Assigning all to NP.")
          object.gene[[new.class.name]][,1] <- "NP"
          converged <- TRUE
        }
        
        if (all(object.gene[[new.class.name]] == old.classes)) {
          converged <- TRUE
        }
        
        old.classes <- object.gene[[new.class.name]]
        n.iter <- n.iter + 1
        
      }
      object.gene[[new.class.name]][object.gene[[new.class.name]] == gene,1] <- paste(gene, " KO", sep = "")
    }
    
    # assign classifications back to original object
    object[[new.class.name]][Cells(object.gene),1] <- object.gene[[new.class.name]]
    
    #add global classifications of KO, NP and NT class
    object[[paste(new.class.name, ".global", sep = "")]] <- as.character(sapply(as.character(object[[new.class.name]][,1]), function(x) strsplit(x, ' (?=[^ ]+$)', perl=TRUE)[[1]][2]))
    object[[paste(new.class.name, ".global", sep = "")]][which(is.na(object[[paste(new.class.name, ".global", sep = "")]])),1] <- nt.class.name
    p_ko[[gene]] <- post.prob 
  }
  #add posterior probabilities to seurat object as a meta data column.
  names(p_ko) <- NULL
  prob <- unlist(p_ko)
  object <- AddMetaData(object, metadata = prob, col.name = "p_ko")
  object$p_ko[names(which(is.na(object$p_ko)))] <- 0
  return(object)
}

#' Function to prepare data for Linear Discriminant Analysis.
#' @param object An object of class Seurat.
#' @param de.assay Assay to use for selection of DE genes.
#' @param pc.assay Assay to use for running Principle components analysis.
#' @param labels Meta data column with target gene class labels.
#' @param nt.label Name of non-targeting cell class.
#' @param npcs Number of principle components to use.
#' @param verbose Print progress bar.
#' @export
PrepLDA <- function( object,
                     de.assay = "RNA",
                     pc.assay = "PRTB",
                     labels = "gene",
                     nt.label = "NT",
                     npcs = 10,
                     verbose = TRUE
){
  projected_pcs <- list()
  gene_list <- setdiff(unique(object[[labels]][,1]),nt.label)
  features <- rownames(object[[de.assay]])
  VariableFeatures(object[[de.assay]]) <- features
  object <- ScaleData(object ,assay = de.assay, do.scale = T, do.center = T)
  Idents(object) <- labels
  for(g in gene_list) {
    if (verbose){
      print(g)
    }
    gene_subset <- subset(object, idents =c(g,nt.label))
    gene_set <- TopDEGenesMixscape(object = gene_subset, ident.1 = g, de.assay = de.assay, pval.cutoff = 0.05, labels = labels)
    if (length(gene_set) < npcs) {
      next;
    }
    DefaultAssay(gene_subset) <- pc.assay
    gene_subset <- RunPCA(gene_subset, npcs = npcs ,verbose = F,features = gene_set, assay = pc.assay)
    project_pca <- t(object[[pc.assay]]@scale.data[rownames(gene_subset[["pca"]]@feature.loadings),])%*%gene_subset[["pca"]]@feature.loadings
    colnames(project_pca) <- paste(g,colnames(project_pca),sep="_")
    projected_pcs[[g]] <- project_pca
  }
  return(projected_pcs)
}

#' Function to perform Linear Discriminant Analysis. 
#' @param object An object of class Seurat.
#' @param labels Meta data column with target gene class labels.
#' @param assay Assay to use for performing Linear Discriminant Analysis (LDA).
#' @param features Features to compute LDA on
#' @param ndims.print Number of LDA dimensions to print.
#' @param nfeatures.print Number of features to print for each LDA component.
#' @param reduction.key Reduction key name.
#' @param seed.use Set to 42 for reproducibility.
#' @rdname RunLDA
#' @export
#' @method RunLDA Assay
#'
RunLDA.Assay <- function(
  object,
  labels,
  assay = NULL,
  features = NULL,
  verbose = TRUE,
  ndims.print = 1:5,
  nfeatures.print = 30,
  reduction.key = "LDA_",
  seed.use = 42,
  ...
) {
  data.use <- PrepDR(
    object = object,
    features = features,
    verbose = verbose
  )
  reduction.data <- RunLDA.default(
    object = t(data.use),
    labels = labels,
    assay = assay,
    verbose = verbose,
    ndims.print = ndims.print,
    nfeatures.print = nfeatures.print,
    reduction.key = reduction.key,
    seed.use = seed.use,
    ...
  )
  return(reduction.data)
}

#' @param reduction.name dimensional reduction name,  lda by default
#' @param object An object of class Seurat.
#' @param labels Meta data column with target gene class labels.
#' @param assay Assay to use for performing Linear Discriminant Analysis (LDA).
#' @param features Features to compute LDA on
#' @param ndims.print Number of LDA dimensions to print.
#' @param nfeatures.print Number of features to print for each LDA component.
#' @param reduction.key Reduction key name.
#' @param seed.use Set to 42 for reproducibility.
#' @param project
#' @rdname RunLDA
#' @export
#' @method RunLDA Seurat
#'
RunLDA.Seurat <- function(
  object,
  labels,
  assay = NULL,
  features = NULL,
  verbose = TRUE,
  ndims.print = 1:5,
  nfeatures.print = 30,
  reduction.name = "lda",
  reduction.key = "LDA_",
  seed.use = 42,
  project = TRUE,
  ...
) {
  assay <- assay %||% DefaultAssay(object = object)
  assay.data <- GetAssay(object = object, assay = assay)
  reduction.data <- RunLDA.Assay(
    object = assay.data,
    assay = assay,
    labels = labels,
    features = features,
    verbose = FALSE,
    ndims.print = ndims.print,
    nfeatures.print = nfeatures.print,
    reduction.key = reduction.key,
    seed.use = seed.use,
    ...
  )
  object[[reduction.name]] <- reduction.data
  object$lda.assignments <- slot(object = object[[reduction.name]], name = "misc")[["assignments"]]
  object <- AddMetaData(object = object, metadata = as.data.frame(x = slot(object = object[[reduction.name]], name = "misc")[["posterior"]]))
  object <- LogSeuratCommand(object = object)
  object <- ProjectDim(object = object,reduction = reduction.name,assay = assay,verbose = verbose,dims.print = ndims.print,nfeatures.print = nfeatures.print)
  object[[reduction.name]]@feature.loadings <- object[[reduction.name]]@feature.loadings.projected
  return(object)
}

#' @param object Input values for LDA (numeric), with observations as rows
#' @param labels Observation labels for LDA
#' @param assay Name of Assay LDA is being run on
#' @param verbose Print the top genes associated with high/low loadings for
#' the PCs
#' @param ndims.print PCs to print genes for
#' @param nfeatures.print Number of genes to print for each PC
#' @param reduction.key dimensional reduction key, specifies the string before
#' the number for the dimension names. LDA by default
#' @param seed.use Set a random seed. By default, sets the seed to 42. Setting
#' NULL will not set a seed.
#'
#' @importFrom MASS lda
#'
#' @rdname RunLDA
#' @export
#'
RunLDA.default <- function(
  object,
  labels,
  assay = NULL,
  verbose = TRUE,
  ndims.print = 1:5,
  nfeatures.print = 30,
  reduction.key = "LDA_",
  seed.use = 42,
  ...
) {
  if (!is.null(x = seed.use)) {
    set.seed(seed = seed.use)
  }
  object <- data.frame(object)
  var_names <- colnames(object)
  object$lda_cluster_label <- labels
  lda_results <- lda(lda_cluster_label ~ ., object,)
  lda_predictions <- predict(object = lda_results, newdata = object)
  lda_cv <-lda(CV=TRUE ,lda_cluster_label ~ ., object,)$posterior
  
  feature.loadings <- lda_results$scaling
  cell.embeddings <- lda_predictions$x
  lda.assignments <- lda_predictions$class
  lda.posterior <- lda_predictions$posterior
  colnames(x = lda.posterior) <- paste0("LDAP_", colnames(x = lda.posterior))
  rownames(x = feature.loadings) <- var_names
  colnames(x = feature.loadings) <- paste0(reduction.key, 1:ncol(cell.embeddings))
  rownames(x = cell.embeddings) <- rownames(x = object)
  colnames(x = cell.embeddings) <- colnames(x = feature.loadings)
  reduction.data <- CreateDimReducObject(
    embeddings = cell.embeddings,
    loadings = feature.loadings,
    assay = assay,
    key = reduction.key,
    misc = list(assignments=lda.assignments, posterior = lda.posterior, model = lda_results, cv = lda_cv)
  )
  if (verbose) {
    print(x = reduction.data, dims = ndims.print, nfeatures = nfeatures.print)
  }
  return(reduction.data)
}

#' DE and EnrichR pathway visualization barplot
#' @param object Name of object class Seurat.
#' @param ident.1 Cell class identity 1.
#' @param ident.2 Cell class identity 2. 
#' @param balanced Option to display pathway enrichments for both negative and positive DE genes.If false, only positive DE gene will be displayed.
#' @param logfc.threshold Limit testing to genes which show, on average, at least X-fold difference (log-scale) between the two groups of cells. Default is 0.25 Increasing logfc.threshold speeds up the function, but can miss weaker signals.
#' @param assay Assay to use for differential expression analysis.
#' @param max.genes Maximun number of genes to use as input to enrichR.
#' @param test.use Denotes which test to use for differential expression. options include: same options as in FindMarkers function.
#' @param p.val.cutoff Cutoff to select DE genes.
#' @param cols A list of colors to use for barplots.
#' @param enirch.database Database to use from enrichR.
#' @param num.pathway Number of pathways to display in barplot.
#' @import enrichR
#' @importFrom ggplot2
#' @importFrom enrichR
#' @inheritParams FindMarkers
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
  assay <- assay %||% DefaultAssay(object = object)
  DefaultAssay(object = object) <- assay
  logfc.threshold <- logfc.threshold
  p.val.cutoff <- p.val.cutoff
  if (is.numeric(x = max.genes)){
    all.markers <- FindMarkers(
      object = object, 
      ident.1 = ident.1,
      ident.2 = ident.2, 
      only.pos = FALSE, 
      logfc.threshold = logfc.threshold,
      test.use = test.use, 
      assay = assay
    )
    if (balanced == TRUE) {
      pos.markers <- all.markers[all.markers[, 2] > logfc.threshold, , drop = FALSE]
      neg.markers <- all.markers[all.markers[, 2] < logfc.threshold, , drop = FALSE]
      if (nrow(x = pos.markers[pos.markers[, 1] < p.val.cutoff, , drop = FALSE]) < max.genes) {
        pos.markers.list <- c(rownames(x = pos.markers[pos.markers[, 1] < p.val.cutoff, , drop = FALSE]))
        
        if (nrow(neg.markers[neg.markers[,1] < p.val.cutoff, , drop = FALSE]) < max.genes) {
          neg.markers.list <- c(rownames(x = neg.markers[neg.markers[, 1] < p.val.cutoff, , drop = FALSE]))
        } else {
          neg.markers.list <- c(rownames(x = neg.markers[neg.markers[, 1] < p.val.cutoff, , drop = FALSE])[1:max.genes])
        }
      } else {
        pos.markers.list <- c(rownames(x = pos.markers[pos.markers[, 1] < p.val.cutoff, , drop = FALSE]))[1:max.genes]
        if (nrow(x = neg.markers[neg.markers[, 1] < p.val.cutoff, , drop = FALSE]) < max.genes) {
          neg.markers.list <- c(rownames(x = neg.markers[neg.markers[, 1] < p.val.cutoff, , drop = FALSE]))
        } else {
          neg.markers.list <- c(rownames(x = neg.markers[neg.markers[, 1] < p.val.cutoff, , drop = FALSE])[1:max.genes])
        }
      }
    } 
    
    else {
      pos.markers <- all.markers[all.markers[, 2] > logfc.threshold, , drop = FALSE]
      if (nrow(pos.markers[pos.markers[, 1] < p.val.cutoff, , drop = FALSE]) < max.genes) {
        pos.markers.list <- c(rownames(x = pos.markers[pos.markers[, 1] < p.val.cutoff, , drop = FALSE]))
      } else{
        pos.markers.list<- c(rownames(x = pos.markers[pos.markers[, 1] < p.val.cutoff, , drop = FALSE]))[1:max.genes]
      }
    }
  }
  else {
    stop("please set max.genes")
  }
  
  if(is.null(enrich.database) == TRUE) {
    stop("please specify the name of enrichR database to use")
  }
  else {
    if (balanced == TRUE){
      
      pos.er <- enrichr(genes = pos.markers.list, databases = enrich.database)
      pos.er <- do.call(what = cbind, pos.er)
      pos.er$log10pval <- -log10(pos.er[,paste(enrich.database,sep = ".", "P.value")])
      pos.er$term <- pos.er[, paste(enrich.database,sep = ".", "Term")]
      pos.er <- pos.er[1:num.pathway, ]
      pos.er$term <- factor(x = pos.er$term, levels = pos.er$term[order(pos.er$log10pval)])
      
      
      neg.er <- enrichr(genes = neg.markers.list, databases = enrich.database)
      neg.er <- do.call(what = cbind, neg.er)
      neg.er$log10pval <- -log10(neg.er[,paste(enrich.database,sep = ".", "P.value")])
      neg.er$term <- neg.er[, paste(enrich.database,sep = ".", "Term")]
      neg.er <- neg.er[1:num.pathway, ]
      neg.er$term <- factor(x = neg.er$term, levels = neg.er$term[order(neg.er$log10pval)])
      
      p1 <- ggplot(data=pos.er, aes(x=term, y=log10pval)) +
        geom_bar(stat="identity") +
        coord_flip() + xlab("Pathway") +
        scale_fill_manual(values = cols, drop = FALSE)+
        ylab("-log10(pval)") +
        ggtitle(paste(enrich.database, ident.1, sep = "_", "positive markers")) + 
        theme_classic() +
        theme(axis.text.y = element_text(size = 12, face = "bold"))
      
      p2 <- ggplot(data=neg.er, aes(x=term, y=log10pval)) +
        geom_bar(stat="identity") +
        coord_flip() + xlab("Pathway") +
        scale_fill_manual(values = cols, drop = FALSE)+
        ylab("-log10(pval)") +
        ggtitle(paste(enrich.database, ident.1, sep = "_", "negative markers")) + 
        theme_classic() +
        theme(axis.text.y = element_text(size = 12, face = "bold"))
      p <-plot_grid(p1,p2)
    }
    
    else{
      pos.er <- enrichr(genes = pos.markers.list, databases = enrich.database)
      pos.er <- do.call(what = cbind, pos.er)
      pos.er$log10pval <- -log10(pos.er[,paste(enrich.database,sep = ".", "P.value")])
      pos.er$term <- pos.er[, paste(enrich.database,sep = ".", "Term")]
      pos.er <- pos.er[1:num.pathway, ]
      pos.er$term <- factor(x = pos.er$term, levels = pos.er$term[order(pos.er$log10pval)])
      
      p <- ggplot(data=pos.er, aes(x=term, y=log10pval)) +
        geom_bar(stat="identity") +
        scale_fill_manual(values = cols, drop = FALSE)+
        coord_flip() + xlab("Pathway") +
        ylab("-log10(pval)") +
        ggtitle(paste(enrich.database, ident.1, sep = "_", "positive markers")) + 
        theme_classic() +
        theme(axis.text.y = element_text(size = 12, face = "bold"))
    }
  }
  if (!return.gene.list){
    return(p)
  }
  else{
    if(balanced){
      return(list(pos = pos.er, neg = neg.er))
    } else {
      return(list(pos = pos.er))
    }
  }
}