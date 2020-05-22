#'  Calculate Perturbation score in cells with gRNAs.
#'  @param object An object of class Seurat.
#'  @param assay Name of Assay PRTB  signature is being calculated on.
#'  @param slot Data slot to use for PRTB score calculation.
#'  @param gd.class Metadata column containing gRNA classification.
#'  @param nt.cell.class Non-targeting gRNA cell classification identity.
#'  @param split.by Provide metadata column if multiple biological replicates exist to calculate PRTB score for every replicate separately.
#'  @param num.neighbors Number of nearest neigbors to consider.
#'  @param ndims Number of dimensions to use from dimensionality reduction method.
#'  @param reduction Reduction method used to calculate nearest neighbors.

CalcPerturbScore <- function ( object, 
                             assay = NULL,
                             gd.class = "guide_ID",
                             nt.cell.class = "NT",
                             num.neighbors = NULL,
                             reduction = "pca", 
                             slot = "data", 
                             split.by = NULL,
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
#' 
Project <- function(v1, v2) {
  return(as.vector((v1 %*% v2) / (v2 %*% v2)))
}

#' Find top DE genes that pass some p value cutoff between cells with targeting and non-targeting gRNAs.
#' @param object An object of class Seurat.
#' @param ident.1 Target gene class to find DE genes for.
#' @param group.by metadata column with target gene classification.
#' @param assay Name of Assay DE is performed on.
#' @param test.use 	Denotes which test to use. See all available tests on FindMarkers documentation.
#' @param pval.cut.off P-value cut-off for selection of significantly DE genes.
#' 
TopDEGenesMixscape <- function(
  object, 
  ident.1, 
  group.by = 'gene', 
  de.assay = "RNA", 
  test.use = "LR", 
  pval.cutoff = 5e-2
) {
  message("Finding new perturbation gene set")
  de.genes <- data.frame()  
  tryCatch(
    expr = {
      de.genes <- FindMarkers(
        object = object, 
        ident.1 = ident.1, 
        group.by = group.by,
        assay = de.assay,
        test.use = test.use
      )
      de.genes <- subset(de.genes, p_val_adj < pval.cutoff)
    }, 
    error = function(e) {}
  )
  return(rownames(de.genes))
}


#' Define Normal distribution - returns list with mu and sd
DefineNormalMixscape <- function(x) {
  mu <- mean(x)
  sd <- sd(x)
  return(list(mu = mu, sd = sd))
}

#' Identify perturbed and non-perturbed gRNA expressing cells.
#' @param An object of class Seurat.
#' @param assay Assay to use.
#' @param slot data slot to use.
#' @param gene.class metadata column with target gene classifications.
#' @param nt.class.name Classification name of non-targeting gRNA cells.
#' @param new.class.name Name of mixscape classification to be stored in metadata.
#' @param min.de.genes Required number of genes that are DE for method to separate perturbed and non-perturbed cells. 
#' @param de.assay Assay to use when performing differential expression analysis. Usually RNA.
#' @param iter.num Number of normalmixEM iterations to run if convergance does not occur.
#' @import mixtools
#' @import reshape2
#' 
RunMixscape <- function( object = NULL,
                         assay = "PRTB",
                         slot = "scale.data",
                         gene.class = "gene",
                         nt.class.name = "NT",
                         new.class.name = "mixscape_class",
                         min.de.genes = 5,
                         de.assay = "RNA",
                         iter.num = 10,
                         ...
                         
){
  if (is.null(assay) == TRUE){
    assay <- DefaultAssay(object)
  }
  
  if (is.null(gene.class) == TRUE){
    stop("Please specify target gene class metadata name")
  }
  #de marker genes
  prtb_markers <- c()
  
  #new metadata column for mixscape classification
  object[[new.class.name]] <- object[[gene.class]]
  
  genes <- setdiff(unique(object[["gene"]][,1]), y = nt.class.name)
  
  #pertubration vectors storage
  gv.list <- list(list())
  
  for (gene in genes){
    message("Processing ", gene)
    Idents(object) <- gene.class
    
    # Get object containing only guide of interest + non-targeting
    object.gene <- subset(object, idents = c(gene, nt.class.name))
    orig.guide.cells <- WhichCells(object.gene, idents = gene)
    DefaultAssay(object.gene) <- assay
    
    # find de genes between guide positive and non-targeting
    de.genes <- TopDEGenesMixscape(object.gene, ident.1 = gene, de.assay = de.assay)
    prtb_markers[[gene]] <- de.genes
    
    # if fewer than 5 DE genes, call all guide cells NP 
    if (length(de.genes) < min.de.genes) {
      message("Fewer than ",min.de.genes, " DE genes for ", gene, ". Assigning cells as NP.")
      object.gene[[new.class.name]][orig.guide.cells,1] <- paste(gene, "NP", sep = "_")
    } else {
      object.gene <- ScaleData(object.gene, features = de.genes, verbose = FALSE) 
      dat <- GetAssayData(object = object.gene[["PRTB"]], slot = slot)[de.genes, ]
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
        gv<- melt(pvec)
        gv$name <- nt.class.name
        gv[intersect(rownames(gv), guide.cells),"name"] <- gene
        gv.list[[gene]][[n.iter+1]] <- gv
        
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
        object.gene[[new.class.name]][names(which(post.prob < 0.5)),1] <- paste(gene, "NP", sep = "_")
        
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
      object.gene[[new.class.name]][object.gene[[new.class.name]] == gene,1] <- paste(gene, "KO", sep = "_")
    }
    
    # assign classifications back to original object
    object[[new.class.name]][Cells(object.gene),1] <- object.gene[[new.class.name]]
    
    #add global classifications of KO, NP and NT class
    object[[paste(new.class.name, ".global", sep = "")]] <- as.character(sapply(as.character(object[[new.class.name]][,1]), function(x) strsplit(x,"_")[[1]][2]))
    object[[paste(new.class.name, ".global", sep = "")]][which(is.na(object[[paste(new.class.name, ".global", sep = "")]])),1] <- nt.class.name
    
  }
  return(object)
}

PrepLDA <- function( object,
                     assay = "PRTB",
                     slot = "scale.data",
                     cut.off = 0.25,
                     labels = "gene",
                     nt.label = "NT",
                     npcs = 10,
                     reduction = "pca",
                     reduction.name = "pc",
                     reduction.key = "pc",
                     verbose = TRUE
  
){
  avg_perturb=AverageExpression(object = object,assays = assay,use.scale = T,slot = slot)
  
  #0.25 cutoff is somewhat arbitrary
  prtb_features <- rownames(which(abs(avg_perturb[[assay]])>cut.off,arr.ind = T))
  VariableFeatures(object[[assay]]) <- prtb_features
  
  # Third goal, classify non-perturbed cells and visualize gene differences
  projected_pcs <- list()
  gene_list <- setdiff(unique(object[[labels]][,1]),nt.label)
  
  for(g in gene_list) {
    
    if (verbose){
      print(g)
    }
    
    Idents(object) <- labels
    gene_subset <- subset(object, idents =c(g,nt.label))
    DefaultAssay(gene_subset) <- assay
    gene_subset <- ScaleData(gene_subset,verbose = F,do.scale = T) %>% RunPCA(npcs = npcs ,verbose = F)
    pca_embed <- gene_subset[[reduction]]@cell.embeddings
    
    project_pca=ProjectCellEmbeddings(reference = gene_subset,query = object ,reference.assay = assay ,query.assay = assay,dims = 1:npcs,verbose = F)
   
    colnames(project_pca) <- paste(g,colnames(project_pca),sep="_")
    projected_pcs[[g]] <- project_pca
  }
  return(projected_pcs)
}

#' @param features Features to compute LDA on
#'
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
#'
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
#' @rdname RunPCA
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