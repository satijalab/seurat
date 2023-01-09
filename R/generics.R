#' @include reexports.R
#'
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Generics
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' Add info to anchor matrix
#'
#' @param anchors An \code{\link{AnchorSet}} object
#' @param vars Variables to pull for each object via FetchData
#' @param slot Slot to pull feature data for
#' @param assay Specify the Assay per object if annotating with expression data
#' @param ... Arguments passed to other methods
#
#' @return Returns the anchor dataframe with additional columns for annotation
#' metadata
#'
#' @export
#'
AnnotateAnchors <- function(anchors, vars, slot, ...) {
  UseMethod(generic = 'AnnotateAnchors', object = anchors)
}

#' Convert objects to CellDataSet objects
#'
#' @param x An object to convert to class \code{CellDataSet}
#' @param ... Arguments passed to other methods
#'
#' @rdname as.CellDataSet
#' @export as.CellDataSet
#'
as.CellDataSet <- function(x, ...) {
  UseMethod(generic = 'as.CellDataSet', object = x)
}

#' Convert objects to SingleCellExperiment objects
#'
#' @param x An object to convert to class \code{SingleCellExperiment}
#' @param ... Arguments passed to other methods
#'
#' @rdname as.SingleCellExperiment
#' @export as.SingleCellExperiment
#'
as.SingleCellExperiment <- function(x, ...) {
  UseMethod(generic = 'as.SingleCellExperiment', object = x)
}

#' Cluster Determination
#'
#' Identify clusters of cells by a shared nearest neighbor (SNN) modularity
#' optimization based clustering algorithm. First calculate k-nearest neighbors
#' and construct the SNN graph. Then optimize the modularity function to
#' determine clusters. For a full description of the algorithms, see Waltman and
#' van Eck (2013) \emph{The European Physical Journal B}. Thanks to Nigel
#' Delaney (evolvedmicrobe@github) for the rewrite of the Java modularity
#' optimizer code in Rcpp!
#'
#' To run Leiden algorithm, you must first install the leidenalg python
#' package (e.g. via pip install leidenalg), see Traag et al (2018).
#'
#' @param object An object
#' @param ... Arguments passed to other methods
#'
#' @return Returns a Seurat object where the idents have been updated with new cluster info;
#' latest clustering results will be stored in object metadata under 'seurat_clusters'.
#' Note that 'seurat_clusters' will be overwritten everytime FindClusters is run
#'
#' @export
#'
#' @rdname FindClusters
#' @export FindClusters
#'
FindClusters <- function(object, ...) {
  UseMethod(generic = 'FindClusters', object = object)
}

#' Gene expression markers of identity classes
#'
#' Finds markers (differentially expressed genes) for identity classes
#'
#' @param object An object
#' @param ... Arguments passed to other methods and to specific DE methods

#' @return data.frame with a ranked list of putative markers as rows, and associated
#' statistics as columns (p-values, ROC score, etc., depending on the test used (\code{test.use})). The following columns are always present:
#' \itemize{
#'   \item \code{avg_logFC}: log fold-chage of the average expression between the two groups. Positive values indicate that the gene is more highly expressed in the first group
#'   \item \code{pct.1}: The percentage of cells where the gene is detected in the first group
#'   \item \code{pct.2}: The percentage of cells where the gene is detected in the second group
#'   \item \code{p_val_adj}: Adjusted p-value, based on bonferroni correction using all genes in the dataset
#' }
#'
#' @details p-value adjustment is performed using bonferroni correction based on
#' the total number of genes in the dataset. Other correction methods are not
#' recommended, as Seurat pre-filters genes using the arguments above, reducing
#' the number of tests performed. Lastly, as Aaron Lun has pointed out, p-values
#' should be interpreted cautiously, as the genes used for clustering are the
#' same genes tested for differential expression.
#'
#' @references McDavid A, Finak G, Chattopadyay PK, et al. Data exploration,
#' quality control and testing in single-cell qPCR-based gene expression experiments.
#' Bioinformatics. 2013;29(4):461-467. doi:10.1093/bioinformatics/bts714
#' @references Trapnell C, et al. The dynamics and regulators of cell fate
#' decisions are revealed by pseudotemporal ordering of single cells. Nature
#' Biotechnology volume 32, pages 381-386 (2014)
#' @references Andrew McDavid, Greg Finak and Masanao Yajima (2017). MAST: Model-based
#' Analysis of Single Cell Transcriptomics. R package version 1.2.1.
#' https://github.com/RGLab/MAST/
#' @references Love MI, Huber W and Anders S (2014). "Moderated estimation of
#' fold change and dispersion for RNA-seq data with DESeq2." Genome Biology.
#' https://bioconductor.org/packages/release/bioc/html/DESeq2.html
#'
#' @export
#'
#' @examples
#' data("pbmc_small")
#' # Find markers for cluster 2
#' markers <- FindMarkers(object = pbmc_small, ident.1 = 2)
#' head(x = markers)
#'
#' # Take all cells in cluster 2, and find markers that separate cells in the 'g1' group (metadata
#' # variable 'group')
#' markers <- FindMarkers(pbmc_small, ident.1 = "g1", group.by = 'groups', subset.ident = "2")
#' head(x = markers)
#'
#' # Pass 'clustertree' or an object of class phylo to ident.1 and
#' # a node to ident.2 as a replacement for FindMarkersNode
#' if (requireNamespace("ape", quietly = TRUE)) {
#'   pbmc_small <- BuildClusterTree(object = pbmc_small)
#'   markers <- FindMarkers(object = pbmc_small, ident.1 = 'clustertree', ident.2 = 5)
#'   head(x = markers)
#' }
#'
#' @rdname FindMarkers
#' @export FindMarkers
#'
#' @aliases FindMarkersNode
#' @seealso \code{FoldChange}
#'
FindMarkers <- function(object, ...) {
  UseMethod(generic = 'FindMarkers', object = object)
}

#' (Shared) Nearest-neighbor graph construction
#'
#' Computes the \code{k.param} nearest neighbors for a given dataset. Can also
#' optionally (via \code{compute.SNN}), construct a shared nearest neighbor
#' graph by calculating the neighborhood overlap (Jaccard index) between every
#' cell and its \code{k.param} nearest neighbors.
#'
#' @param object An object
#' @param ... Arguments passed to other methods
#'
#' @return This function can either return a \code{\link{Neighbor}} object
#' with the KNN information or a list of \code{\link{Graph}} objects with
#' the KNN and SNN depending on the settings of \code{return.neighbor} and
#' \code{compute.SNN}. When running on a \code{\link{Seurat}} object, this
#' returns the \code{\link{Seurat}} object with the Graphs or Neighbor objects
#' stored in their respective slots. Names of the Graph or Neighbor object can
#' be found with \code{\link{Graphs}} or \code{\link{Neighbors}}.
#'
#' @examples
#' data("pbmc_small")
#' pbmc_small
#' # Compute an SNN on the gene expression level
#' pbmc_small <- FindNeighbors(pbmc_small, features = VariableFeatures(object = pbmc_small))
#'
#' # More commonly, we build the SNN on a dimensionally reduced form of the data
#' # such as the first 10 principle components.
#'
#' pbmc_small <- FindNeighbors(pbmc_small, reduction = "pca", dims = 1:10)
#'
#' @rdname FindNeighbors
#' @export FindNeighbors
#'
FindNeighbors <- function(object, ...) {
  UseMethod(generic = 'FindNeighbors', object = object)
}

#' Find variable features
#'
#' Identifies features that are outliers on a 'mean variability plot'.
#'
#' For the mean.var.plot method:
#' Exact parameter settings may vary empirically from dataset to dataset, and
#' based on visual inspection of the plot. Setting the y.cutoff parameter to 2
#' identifies features that are more than two standard deviations away from the
#' average dispersion within a bin. The default X-axis function is the mean
#' expression level, and for Y-axis it is the log(Variance/mean). All mean/variance
#' calculations are not performed in log-space, but the results are reported in
#' log-space - see relevant functions for exact details.
#'
#' @param object An object
#' @param ... Arguments passed to other methods
#'
#' @rdname FindVariableFeatures
#' @export FindVariableFeatures
#'
#' @aliases FindVariableGenes
#'
FindVariableFeatures <- function(object, ...) {
  UseMethod(generic = 'FindVariableFeatures', object = object)
}

#' Find spatially variable features
#'
#' Identify features whose variability in expression can be explained to some
#' degree by spatial location.
#'
#' @param object An object
#' @param ... Arguments passed to other methods
#'
#' @rdname FindSpatiallyVariableFeatures
#' @export FindSpatiallyVariableFeatures
#'
FindSpatiallyVariableFeatures <- function(object, ...) {
  UseMethod(generic = 'FindSpatiallyVariableFeatures', object = object)
}

#' Fold Change
#'
#' Calculate log fold change and percentage of cells expressing each feature
#' for different identity classes.
#'
#' If the slot is \code{scale.data} or a reduction is specified, average difference
#' is returned instead of log fold change and the column is named "avg_diff".
#' Otherwise, log2 fold change is returned with column named "avg_log2_FC".
#'
#' @examples
#' data("pbmc_small")
#' FoldChange(pbmc_small, ident.1 = 1)
#'
#' @param object A Seurat object
#' @param ... Arguments passed to other methods
#' @rdname FoldChange
#' @export FoldChange
#' @return Returns a data.frame
#' @seealso \code{FindMarkers}
FoldChange <- function(object, ...) {
  UseMethod(generic = 'FoldChange', object = object)
}

#' Get an Assay object from a given Seurat object.
#'
#' @param object An object
#' @param ... Arguments passed to other methods
#'
#' @return Returns an Assay object
#'
#' @rdname GetAssay
#' @export GetAssay
#'
GetAssay <- function(object, ...) {
  UseMethod(generic = 'GetAssay', object = object)
}

#' Integrate low dimensional embeddings
#'
#' Perform dataset integration using a pre-computed Anchorset of specified low
#' dimensional representations.
#'
#' The main steps of this procedure are identical to \code{\link{IntegrateData}}
#' with one key distinction. When computing the weights matrix, the distance
#' calculations are performed in the full space of integrated embeddings when
#' integrating more than two datasets, as opposed to a reduced PCA space which
#' is the default behavior in \code{\link{IntegrateData}}.
#'
#' @param anchorset An AnchorSet object
#' @param new.reduction.name Name for new integrated dimensional reduction.
#' @param reductions Name of reductions to be integrated. For a
#' TransferAnchorSet, this should be the name of a reduction present in the
#' anchorset object (for example, "pcaproject"). For an IntegrationAnchorSet,
#' this should be a \code{\link{DimReduc}} object containing all cells present
#' in the anchorset object.
#' @param dims.to.integrate Number of dimensions to return integrated values for
#' @param weight.reduction Dimension reduction to use when calculating anchor
#' weights. This can be one of:
#' \itemize{
#'    \item{A string, specifying the name of a dimension reduction present in
#'    all objects to be integrated}
#'    \item{A vector of strings, specifying the name of a dimension reduction to
#'    use for each object to be integrated}
#'    \item{A vector of \code{\link{DimReduc}} objects, specifying the object to
#'    use for each object in the integration}
#'    \item{NULL, in which case the full corrected space is used for computing
#'    anchor weights.}
#' }
#' @param ... Reserved for internal use
#'
#' @return When called on a TransferAnchorSet (from FindTransferAnchors), this
#' will return the query object with the integrated embeddings stored in a new
#' reduction. When called on an IntegrationAnchorSet (from IntegrateData), this
#' will return a merged object with the integrated reduction stored.
#'
#' @rdname IntegrateEmbeddings
#' @export IntegrateEmbeddings
#'
IntegrateEmbeddings <- function(anchorset, ...) {
  UseMethod(generic = "IntegrateEmbeddings", object = anchorset)
}

#' @export
#'
LeverageScore <- function(object, ...) {
  UseMethod(generic = 'LeverageScore', object = object)
}

#' Normalize Raw Data
#'
#' @param data Matrix with the raw count data
#' @param scale.factor Scale the data; default is \code{1e4}
#' @param verbose Print progress
#'
#' @return A matrix with the normalized and log-transformed data
#'
#' @template param-dotsm
#'
#' @export
#' @concept preprocessing
#'
#' @examples
#' mat <- matrix(data = rbinom(n = 25, size = 5, prob = 0.2), nrow = 5)
#' mat
#' mat_norm <- LogNormalize(data = mat)
#' mat_norm
#'
LogNormalize <- function(
  data,
  scale.factor = 1e4,
  # margin = 2L,
  verbose = TRUE,
  ...
) {
  UseMethod(generic = 'LogNormalize', object = data)
}


#' Metric for evaluating mapping success
#'
#' This metric was designed to help identify query cells that aren't well
#' represented in the reference dataset. The intuition for the score is that we
#' are going to project the query cells into a reference-defined space and then
#' project them back onto the query. By comparing the neighborhoods before and
#' after projection, we identify cells who's local neighborhoods are the most
#' affected by this transformation. This could be because there is a population
#' of query cells that aren't present in the reference or the state of the cells
#' in the query is significantly different from the equivalent cell type in the
#' reference.
#'
#' @param anchors Set of anchors
#' @param ... Arguments passed to other methods
#'
#' @rdname MappingScore
#' @export MappingScore
#'
MappingScore <- function(anchors, ...) {
  UseMethod(generic = "MappingScore", object = anchors)
}

#' Normalize Data
#'
#' Normalize the count data present in a given assay.
#'
#' @param object An object
#' @param ... Arguments passed to other methods
#'
#' @return Returns object after normalization
#'
#' @rdname NormalizeData
#' @export NormalizeData
#'
NormalizeData <- function(object, ...) {
  UseMethod(generic = 'NormalizeData', object = object)
}

#' Project query data to the reference dimensional reduction
#'
#'
#' @param query An object for query cells
#' @param reference An object for reference cells
#' @param query.assay Assay name for query object
#' @param reference.assay Assay name for reference object
#' @param reduction Name of dimensional reduction from reference object
#' @param dims Dimensions used for reference dimensional reduction
#' @param scale Determine if scale query data based on reference data variance
#' @param verbose Print progress
#' @param feature.mean Mean of features in reference
#' @param feature.sd Standard variance of features in reference
#'
#' @return A matrix with projected cell embeddings
#'
#' @rdname ProjectCellEmbeddings
#' @export ProjectCellEmbeddings
#'
#' @keywords internal
#'
ProjectCellEmbeddings <- function(
  query,
  ...
) {
  UseMethod(generic = 'ProjectCellEmbeddings', object = query)
}

#' Project query into UMAP coordinates of a reference
#'
#' This function will take a query dataset and project it into the coordinates
#' of a provided reference UMAP. This is essentially a wrapper around two steps:
#' \itemize{
#'   \item{FindNeighbors - Find the nearest reference cell neighbors and their
#'   distances for each query cell.}
#'   \item{RunUMAP - Perform umap projection by providing the neighbor set
#'   calculated above and the umap model previously computed in the reference.}
#' }
#'
#' @param query Query dataset
#'
#' @rdname ProjectUMAP
#' @export ProjectUMAP
#'
ProjectUMAP <- function(query, ...) {
  UseMethod(generic = "ProjectUMAP", object = query)
}

#' Pseudobulk Expression
#'
#' Normalize the count data present in a given assay.
#'
#' @param object An assay
#' @param ... Arguments passed to other methods
#'
#' @return Returns object after normalization
#'
#' @rdname PseudobulkExpression
#' @export PseudobulkExpression
#'
PseudobulkExpression <- function(object, ...) {
  UseMethod(generic = "PseudobulkExpression", object = object)
}

#' Perform Canonical Correlation Analysis
#'
#' Runs a canonical correlation analysis using a diagonal implementation of CCA.
#' For details about stored CCA calculation parameters, see
#' \code{PrintCCAParams}.
#' @param object1 First Seurat object
#' @param object2 Second Seurat object.
# @param ... Arguments passed to other methods
#'
#' @return Returns a combined Seurat object with the CCA results stored.
#'
#' @seealso \code{\link{merge.Seurat}}
#'
#' @examples
#' data("pbmc_small")
#' pbmc_small
#' # As CCA requires two datasets, we will split our test object into two just for this example
#' pbmc1 <- subset(pbmc_small, cells = colnames(pbmc_small)[1:40])
#' pbmc2 <- subset(pbmc_small, cells = colnames(x = pbmc_small)[41:80])
#' pbmc1[["group"]] <- "group1"
#' pbmc2[["group"]] <- "group2"
#' pbmc_cca <- RunCCA(object1 = pbmc1, object2 = pbmc2)
#' # Print results
#' print(x = pbmc_cca[["cca"]])
#'
#' @rdname RunCCA
#' @export RunCCA
#'
RunCCA <- function(object1, object2, ...) {
  UseMethod(generic = 'RunCCA', object = object1)
}


#' Run Graph Laplacian Eigendecomposition
#'
#' Run a graph laplacian dimensionality reduction. It is used as a low
#' dimensional representation for a cell-cell graph. The input graph
#' should be symmetric
#'
#' @param object A Seurat object
#' @param ... Arguments passed to
#' \code{\link[RSpectra:eigs_sym]{RSpectra::eigs_sym}}
#'
#' @return Returns Seurat object with the Graph laplacian eigenvector
#' calculation stored in the reductions slot
#'
#' @rdname RunGraphLaplacian
#' @export RunGraphLaplacian
#'

RunGraphLaplacian <- function(object, ...) {
  UseMethod(generic = 'RunGraphLaplacian', object = object)
}


#' Run Independent Component Analysis on gene expression
#'
#' Run fastica algorithm from the ica package for ICA dimensionality reduction.
#' For details about stored ICA calculation parameters, see
#' \code{PrintICAParams}.
#'
#' @param object Seurat object
#'
#' @rdname RunICA
#' @export RunICA
#'
RunICA <- function(object, ...) {
  UseMethod(generic = "RunICA", object = object)
}

#' Run Linear Discriminant Analysis
#'
#'
#' @param object An object
#' @param ... Arguments passed to other methods
#'
#' @rdname RunLDA
#' @export RunLDA
#'
#' @aliases RunLDA
#'
RunLDA <- function(object, ...) {
  UseMethod(generic = 'RunLDA', object = object)
}

#' Run Principal Component Analysis
#'
#' Run a PCA dimensionality reduction. For details about stored PCA calculation
#' parameters, see \code{PrintPCAParams}.
#'
#' @param object An object
#' @param ... Arguments passed to other methods and IRLBA
#'
#' @return Returns Seurat object with the PCA calculation stored in the reductions slot
#'
#' @export
#'
#' @rdname RunPCA
#' @export RunPCA
#'
RunPCA <- function(object, ...) {
  UseMethod(generic = 'RunPCA', object = object)
}

#' Run Supervised Latent Semantic Indexing
#'
#' Run a supervised LSI (SLSI) dimensionality reduction supervised by a
#' cell-cell kernel. SLSI is used to capture a linear transformation of peaks
#' that maximizes its dependency to the given cell-cell kernel.
#'
#' @param object An object
#' @param ... Arguments passed to IRLBA irlba
#'
#' @return Returns Seurat object with the SLSI calculation stored in the
#' reductions slot
#'
#' @export
#'
#' @rdname RunSLSI
#' @export RunSLSI
#'
RunSLSI <- function(object, ...) {
  UseMethod(generic = 'RunSLSI', object = object)
}

#' Run Supervised Principal Component Analysis
#'
#' Run a supervised PCA (SPCA) dimensionality reduction supervised by a cell-cell kernel.
#' SPCA is used to capture a linear transformation which maximizes its dependency to
#' the given cell-cell kernel. We use SNN graph as the kernel to supervise the linear
#' matrix factorization.
#'
#' @param object An object
#' @param ... Arguments passed to other methods and IRLBA
#'
#' @return Returns Seurat object with the SPCA calculation stored in the reductions slot
#' @references Barshan E, Ghodsi A, Azimifar Z, Jahromi MZ.
#' Supervised principal component analysis: Visualization, classification and
#' regression on subspaces and submanifolds.
#' Pattern Recognition. 2011 Jul 1;44(7):1357-71. \url{https://www.sciencedirect.com/science/article/pii/S0031320310005819?casa_token=AZMFg5OtPnAAAAAA:_Udu7GJ7G2ed1-XSmr-3IGSISUwcHfMpNtCj-qacXH5SBC4nwzVid36GXI3r8XG8dK5WOQui};
#' @export
#'
#' @rdname RunSPCA
#' @export RunSPCA
#'
RunSPCA <- function(object, ...) {
  UseMethod(generic = 'RunSPCA', object = object)
}

#' Run t-distributed Stochastic Neighbor Embedding
#'
#' Run t-SNE dimensionality reduction on selected features. Has the option of
#' running in a reduced dimensional space (i.e. spectral tSNE, recommended),
#' or running based on a set of genes. For details about stored TSNE calculation
#' parameters, see \code{PrintTSNEParams}.
#'
#' @param object Seurat object
#' @param ... Arguments passed to other methods and to t-SNE call (most commonly used is perplexity)
#'
#' @rdname RunTSNE
#' @export RunTSNE
#'
RunTSNE <- function(object, ...) {
  UseMethod(generic = 'RunTSNE', object = object)
}

#' Run UMAP
#'
#' Runs the Uniform Manifold Approximation and Projection (UMAP) dimensional
#' reduction technique. To run using \code{umap.method="umap-learn"}, you must
#' first install the umap-learn python package (e.g. via
#' \code{pip install umap-learn}). Details on this package can be
#' found here: \url{https://github.com/lmcinnes/umap}. For a more in depth
#' discussion of the mathematics underlying UMAP, see the ArXiv paper here:
#' \url{https://arxiv.org/abs/1802.03426}.
#'
#' @param object An object
#' @param ... Arguments passed to other methods and UMAP
#'
#' @return Returns a Seurat object containing a UMAP representation
#'
#' @references McInnes, L, Healy, J, UMAP: Uniform Manifold Approximation and
#' Projection for Dimension Reduction, ArXiv e-prints 1802.03426, 2018
#'
#' @export
#'
#' @examples
#' \dontrun{
#' data("pbmc_small")
#' pbmc_small
#' # Run UMAP map on first 5 PCs
#' pbmc_small <- RunUMAP(object = pbmc_small, dims = 1:5)
#' # Plot results
#' DimPlot(object = pbmc_small, reduction = 'umap')
#' }
#'
#' @rdname RunUMAP
#' @export RunUMAP
#'
RunUMAP <- function(object, ...) {
  UseMethod(generic = 'RunUMAP', object = object)
}

#' Scale and center the data.
#'
#' Scales and centers features in the dataset. If variables are provided in vars.to.regress,
#' they are individually regressed against each feature, and the resulting residuals are
#' then scaled and centered.
#'
#' ScaleData now incorporates the functionality of the function formerly known
#' as RegressOut (which regressed out given the effects of provided variables
#' and then scaled the residuals). To make use of the regression functionality,
#' simply pass the variables you want to remove to the vars.to.regress parameter.
#'
#' Setting center to TRUE will center the expression for each feature by subtracting
#' the average expression for that feature. Setting scale to TRUE will scale the
#' expression level for each feature by dividing the centered feature expression
#' levels by their standard deviations if center is TRUE and by their root mean
#' square otherwise.
#'
#' @param object An object
#' @param ... Arguments passed to other methods
#'
#' @rdname ScaleData
#' @export ScaleData
#'
ScaleData <- function(object, ...) {
  UseMethod(generic = 'ScaleData', object = object)
}

#' Get image scale factors
#'
#' @param object An object to get scale factors from
#' @param ... Arguments passed to other methods
#'
#' @return An object of class \code{scalefactors}
#'
#' @rdname ScaleFactors
#' @export ScaleFactors
#'
ScaleFactors <- function(object, ...) {
  UseMethod(generic = 'ScaleFactors', object = object)
}

#' Compute Jackstraw scores significance.
#'
#' Significant PCs should show a p-value distribution that is
#' strongly skewed to the left compared to the null distribution.
#' The p-value for each PC is based on a proportion test comparing the number
#' of features with a p-value below a particular threshold (score.thresh), compared with the
#' proportion of features expected under a uniform distribution of p-values.
#'
#' @param object An object
#' @param ... Arguments passed to other methods
#'
#' @return Returns a Seurat object
#'
#' @author Omri Wurtzel
#' @seealso \code{\link{JackStrawPlot}}
#'
#' @rdname ScoreJackStraw
#' @export ScoreJackStraw
#'
ScoreJackStraw <- function(object, ...) {
  UseMethod(generic = 'ScoreJackStraw', object = object)
}

#' Perform sctransform-based normalization
#' @param object An object
#' @param ... Arguments passed to other methods (not used)
#'
#' @rdname SCTransform
#' @export SCTransform
#'
SCTransform <- function(object, ...) {
  UseMethod(generic = 'SCTransform', object = object)
}

#' Get SCT results from an Assay
#'
#' Pull the \code{\link{SCTResults}} information from an \code{\link{SCTAssay}}
#' object.
#'
#' @param object An object
#' @param ... Arguments passed to other methods (not used)
#'
#' @rdname SCTResults
#' @export SCTResults
#'
SCTResults <- function(object, ...) {
  UseMethod(generic = 'SCTResults', object = object)
}


#' @param value new data to set
#'
#' @rdname SCTResults
#' @export SCTResults<-
#'
"SCTResults<-" <- function(object, ..., value) {
  UseMethod(generic = 'SCTResults<-', object = object)
}

#' Variance Stabilizing Transformation
#'
#' Apply variance stabilizing transformation for selection of variable features
#'
#' @inheritParams stats::loess
#' @param data A matrix-like object
#' @param margin Unused
#' @param nselect Number of of features to select
#' @param clip Upper bound for values post-standardization; defaults to the
#' square root of the number of cells
#' @param verbose ...
#'
#' @template param-dotsm
#'
#' @return A data frame with the following columns:
#' \itemize{
#'  \item \dQuote{\code{mean}}: ...
#'  \item \dQuote{\code{variance}}: ...
#'  \item \dQuote{\code{variance.expected}}: ...
#'  \item \dQuote{\code{variance.standardized}}: ...
#'  \item \dQuote{\code{variable}}: \code{TRUE} if the feature selected as
#'   variable, otherwise \code{FALSE}
#'  \item \dQuote{\code{rank}}: If the feature is selected as variable, then how
#'   it compares to other variable features with lower ranks as more variable;
#'   otherwise, \code{NA}
#' }
#'
#' @rdname VST
#' @export VST
#'
#' @keywords internal
#'
VST <- function(
  data,
  margin = 1L,
  nselect = 2000L,
  span = 0.3,
  clip = NULL,
  ...
) {
  UseMethod(generic = 'VST', object = data)
}
