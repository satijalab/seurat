#' @rdname Convert
#' @export as.SingleCellExperiment
#' @aliases as.SingleCellExperiment
#'
as.SingleCellExperiment <- function(from) {
  UseMethod(generic = 'as.SingleCellExperiment', object = from)
}

#' @rdname Convert
#' @export as.seurat
#' @aliases as.seurat
#'
as.seurat <- function(from) {
  UseMethod(generic = 'as.seurat', object = from)
}

#' SNN Graph Construction
#'
#' Constructs a Shared Nearest Neighbor (SNN) Graph for a given dataset. We
#' first determine the k-nearest neighbors of each cell. We use this knn graph
#' to construct the SNN graph by calculating the neighborhood overlap
#' (Jaccard index) between every cell and its k.param nearest neighbors.
#'
#' @param object Seurat object
#' @param k.param Defines k for the k-nearest neighbor algorithm
#' @param prune.SNN Sets the cutoff for acceptable Jaccard index when
#' computing the neighborhood overlap for the SNN construction. Any edges with
#' values less than or equal to this will be set to 0 and removed from the SNN
#' graph. Essentially sets the strigency of pruning (0 --- no pruning, 1 ---
#' prune everything).
#' @param nn.eps Error bound when performing nearest neighbor seach using RANN;
#' default of 0.0 implies exact nearest neighbor search
#' @param verbose Whether or not to print output to the console
#' @param force.recalc Force recalculation of SNN.
#'
#' @importFrom RANN nn2
#' @importFrom igraph plot.igraph graph.adjlist graph.adjacency E
#' @importFrom Matrix sparseMatrix
#' @return Returns the object with object@@snn filled
#'
#' @examples
#' pbmc_small
#' # Compute an SNN on the gene expression level
#' pbmc_small <- BuildSNN(pbmc_small, features = VariableFeatures(object = pbmc_small))
#'
#' # More commonly, we build the SNN on a dimensionally reduced form of the data
#' # such as the first 10 principle components.
#'
#' pbmc_small <- BuildSNN(pbmc_small, reduction = "pca", dims = 1:10)
#'
#' @rdname BuildSNN
#' @export BuildSNN
#'
BuildSNN <- function(
  object,
  k.param,
  prune.SNN,
  nn.eps,
  verbose,
  force.recalc,
  ...
) {
  UseMethod(generic = 'BuildSNN', object = object)
}

#' Get SeuratCommands
#'
#' Pull information on previously run commands in the Seurat object.
#'
#' @param object An object
#' @param command Name of the command to pull
#' @param value Name of the parameter to pull the value for
#'
#' @return Either the SeuratCommand object or the paramter value
#'
#' @rdname Command
#' @export Command
#'
Command <- function(object, command, ..., value) {
  UseMethod(generic = 'Command', object = object)
}

#' Convert Seurat objects to other classes and vice versa
#'
#' Currently, we support direct conversion to/from loom (\url{http://loompy.org/}),
#' SingleCellExperiment (\url{https://bioconductor.org/packages/release/bioc/html/SingleCellExperiment.html}),
#' and Anndata(\url{https://anndata.readthedocs.io/en/latest/}) objects.
#'
#' @param from Object to convert from
#' @param to Class of object to convert to
#' @param ... Arguments passed to and from other methods
#'
#' @return An object of class \code{to}
#'
#' @rdname Convert
#' @export Convert
#'
Convert <- function(from, ...) {
  UseMethod(generic = 'Convert', object = from)
}

#' Get the default assay
#'
#' @param object An object
#' @param ... Arguments passed to other methods
#'
#' @return The name of the default assay
#'
#' @rdname DefaultAssay
#' @export DefaultAssay
#'
DefaultAssay <- function(object, ...) {
  UseMethod(generic = 'DefaultAssay', object = object)
}

#' @inheritParams DefaultAssay
#' @param value Name of assay to set as default
#'
#' @return An object with the new default assay
#'
#' @rdname DefaultAssay
#' @export DefaultAssay<-
#'
"DefaultAssay<-" <- function(object, ..., value) {
  UseMethod(generic = 'DefaultAssay<-', object = object)
}

#' Get cell embeddings
#'
#' @param object An object
#' @param ... Arguments passed to other methods
#'
#' @rdname Embeddings
#' @export Embeddings
#'
Embeddings <- function(object, ...) {
  UseMethod(generic = 'Embeddings', object = object)
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
#' @param object Seurat object
#' @param modularity.fxn Modularity function (1 = standard; 2 = alternative).
#' @param resolution Value of the resolution parameter, use a value above
#' (below) 1.0 if you want to obtain a larger (smaller) number of communities.
#' @param algorithm Algorithm for modularity optimization (1 = original Louvain
#' algorithm; 2 = Louvain algorithm with multilevel refinement; 3 = SLM
#' algorithm).
#' @param n.start Number of random starts.
#' @param n.iter Maximal number of iterations per random start.
#' @param random.seed Seed of the random number generator.
#' @param temp.file.location Directory where intermediate files will be written.
#' Specify the ABSOLUTE path.
#' @param edge.file.name Edge file to use as input for modularity optimizer jar.
#' @param verbose Print output
#'
#' @importFrom Matrix sparseMatrix
#' @importFrom methods .hasSlot
#' @importFrom igraph plot.igraph graph.adjlist
#'
#' @return Returns a Seurat object and optionally the SNN matrix,
#'         object idents have been updated with new cluster info
#'
#' @export
#'
#' @rdname FindClusters
#' @export FindClusters
#'
FindClusters <- function(
  object,
  modularity.fxn,
  resolution,
  algorithm,
  n.start,
  n.iter,
  random.seed,
  temp.file.location,
  edge.file.name,
  verbose,
  ...
) {
  UseMethod(generic = 'FindClusters', object = object)
}

#' Gene expression markers of identity classes
#'
#' Finds markers (differentially expressed genes) for identity classes
#'
#' @param object Seurat object
#' @param features Genes to test. Default is to use all genes
#' @param logfc.threshold Limit testing to genes which show, on average, at least
#' X-fold difference (log-scale) between the two groups of cells. Default is 0.25
#' Increasing logfc.threshold speeds up the function, but can miss weaker signals.
#' @param test.use Denotes which test to use. Available options are:
#' \itemize{
#'  \item{"wilcox"} : Identifies differentially expressed genes between two
#'  groups of cells using a Wilcoxon Rank Sum test (default)
#'  \item{"bimod"} : Likelihood-ratio test for single cell gene expression,
#'  (McDavid et al., Bioinformatics, 2013)
#'  \item{"roc"} : Identifies 'markers' of gene expression using ROC analysis.
#'  For each gene, evaluates (using AUC) a classifier built on that gene alone,
#'  to classify between two groups of cells. An AUC value of 1 means that
#'  expression values for this gene alone can perfectly classify the two
#'  groupings (i.e. Each of the cells in cells.1 exhibit a higher level than
#'  each of the cells in cells.2). An AUC value of 0 also means there is perfect
#'  classification, but in the other direction. A value of 0.5 implies that
#'  the gene has no predictive power to classify the two groups. Returns a
#'  'predictive power' (abs(AUC-0.5)) ranked matrix of putative differentially
#'  expressed genes.
#'  \item{"t"} : Identify differentially expressed genes between two groups of
#'  cells using the Student's t-test.
#'  \item{"negbinom"} : Identifies differentially expressed genes between two
#'   groups of cells using a negative binomial generalized linear model.
#'   Use only for UMI-based datasets
#'  \item{"poisson"} : Identifies differentially expressed genes between two
#'   groups of cells using a poisson generalized linear model.
#'   Use only for UMI-based datasets
#'  \item{"LR"} : Uses a logistic regression framework to determine differentially
#'  expressed genes. Constructs a logistic regression model predicting group
#'  membership based on each feature individually and compares this to a null
#'  model with a likelihood ratio test.
#'  \item{"MAST} : Identifies differentially expressed genes between two groups
#'  of cells using a hurdle model tailored to scRNA-seq data. Utilizes the MAST
#'  package to run the DE testing.
#'  \item{"DESeq2} : Identifies differentially expressed genes between two groups
#'  of cells based on a model using DESeq2 which uses a negative binomial
#'  distribution (Love et al, Genome Biology, 2014).This test does not support
#'  pre-filtering of genes based on average difference (or percent detection rate)
#'  between cell groups. However, genes may be pre-filtered based on their
#'  minimum detection rate (min.pct) across both cell groups. To use this method,
#'  please install DESeq2, using the instructions at
#'  https://bioconductor.org/packages/release/bioc/html/DESeq2.html
#' }
#' @param min.pct  only test genes that are detected in a minimum fraction of
#' min.pct cells in either of the two populations. Meant to speed up the function
#' by not testing genes that are very infrequently expressed. Default is 0.1
#' @param min.diff.pct  only test genes that show a minimum difference in the
#' fraction of detection between the two groups. Set to -Inf by default
#' @param only.pos Only return positive markers (FALSE by default)
#' @param print.bar Print a progress bar once expression testing begins (uses
#' pbapply to do this)
#' @param max.cells.per.ident Down sample each identity class to a max number.
#' Default is no downsampling. Not activated by default (set to Inf)
#' @param random.seed Random seed for downsampling
#' @param latent.vars Variables to test, used only when \code{test.use} is one of
#' 'negbinom', 'poisson', or 'MAST'
#' @param min.cells.features Minimum number of cells expressing the feature in at least one
#' of the two groups, currently only used for poisson and negative binomial tests
#' @param min.cells.group Minimum number of cells in one of the groups
#' @param pseudocount.use Pseudocount to add to averaged expression values when
#' calculating logFC. 1 by default.
#' @param \dots Additional parameters to pass to specific DE functions

#' @return Matrix containing a ranked list of putative markers, and associated
#' statistics (p-values, ROC score, etc.)
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
#' Biotechnology volume 32, pages 381–386 (2014)
#' @references Andrew McDavid, Greg Finak and Masanao Yajima (2017). MAST: Model-based
#' Analysis of Single Cell Transcriptomics. R package version 1.2.1.
#' https://github.com/RGLab/MAST/
#' @references Love MI, Huber W and Anders S (2014). "Moderated estimation of
#' fold change and dispersion for RNA-seq data with DESeq2." Genome Biology.
#' https://bioconductor.org/packages/release/bioc/html/DESeq2.html
#'
#' @import pbapply
#' @importFrom lmtest lrtest
#'
#' @export
#'
#' @examples
#' markers <- FindMarkers(object = pbmc_small, ident.1 = 2)
#' head(markers)
#'
#' @rdname FindMarkers
#' @export FindMarkers
#'
FindMarkers <- function(
  object,
  features,
  logfc.threshold,
  test.use,
  min.pct,
  min.diff.pct,
  verbose,
  only.pos,
  max.cells.per.ident,
  random.seed,
  latent.vars,
  min.cells.feature,
  min.cells.group,
  pseudocount.use,
  ...
) {
  UseMethod(generic = 'FindMarkers', object = object)
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
#' @param selection.method How to choose top variable features. Choose one of :
#' \itemize{
#'   \item{vst:}{ First, fits a line to the relationship of log(variance) and
#'   log(mean) using local polynomial regression (loess). Then standardizes the
#'   feature values using the observed mean and expected variance (given by the
#'   fitted line). Feature variance is the calculated on the standardized values
#'   after clipping to a maximum (default is 50).}
#'   \item{mean.var.plot:}{ First, uses a function to calculate average
#'   expression (mean.function) and dispersion (dispersion.function) for each
#'   feature. Next, divides features into num.bin (deafult 20) bins based on
#'   their average expression, and calculates z-scores for dispersion within
#'   each bin. The purpose of this is to identify variable features while
#'   controlling for the strong relationship between variability and average
#'   expression.}
#'   \item{dispersion:}{ selects the genes with the highest dispersion values}
#' }
#' @param loess.span (vst method) Loess span parameter used when fitting the
#' variance-mean relationship
#' @param clip.max (vst method) After standardization values larger than
#' clip.max will be set to clip.max; default is 'auto' which sets this value to
#' the square root of the number of cells
#' @param mean.function Function to compute x-axis value (average expression).
#'  Default is to take the mean of the detected (i.e. non-zero) values
#' @param dispersion.function Function to compute y-axis value (dispersion).
#' Default is to take the standard deviation of all values
#' @param num.bin Total number of bins to use in the scaled analysis (default
#' is 20)
#' @param binning.method Specifies how the bins should be computed. Available
#' methods are:
#' \itemize{
#'   \item{equal_width:}{ each bin is of equal width along the x-axis [default]}
#'   \item{equal_frequency:}{ each bin contains an equal number of features (can
#'   increase statistical power to detect overdispersed features at high
#'   expression values, at the cost of reduced resolution along the x-axis)}
#' }
#' @param verbose show progress bar for calculations
#'
#' @rdname FindVariableFeatures
#' @export FindVariableFeatures
#'
#' @aliases FindVariableGenes
#'
FindVariableFeatures <- function(
  object,
  selection.method,
  loess.span,
  clip.max,
  mean.function,
  dispersion.function,
  num.bin,
  binning.method,
  verbose,
  ...
) {
  UseMethod(generic = 'FindVariableFeatures', object = object)
}

#' Get an assay from an object
#'
#' @param object An object
#' @param assay Assay to get
#' @param ... Arguments passed to other methods
#'
#' @return Returns an Assay object
#'
#' @rdname GetAssay
#' @export GetAssay
#'
GetAssay <- function(object, assay, ...) {
  UseMethod(generic = 'GetAssay', object = object)
}

#' Accessor function for multimodal data
#'
#' Pull information for specified stored dimensional reduction analysis
#'
#' @param object An object
#' @param slot Specific information to pull (i.e. raw.data, data, scale.data,...). Default is data
#' @param ... Arguments passed to other methods
#'
#' @return Returns assay data
#'
#' @rdname GetAssayData
#' @export GetAssayData
#'
GetAssayData <- function(object, slot, ...) {
  UseMethod(generic = 'GetAssayData', object = object)
}

#' Get highly variable feature information
#'
#' @param object An object
#' @param ... Arguments passed to other methods
#'
#' @return A dataframe with feature means, dispersion, and scaled dispersion
#'
#' @rdname HVFInfo
#' @export HVFInfo
#'
HVFInfo <- function(object, ...) {
  UseMethod(generic = 'HVFInfo', object = object)
}

#' Get an object's cell identities
#'
#' @param object An object
#' @param ... Arguments passed to other methods
#'
#' @return \code{Idents}: The cell identies
#'
#' @rdname Idents
#' @export Idents
#'
Idents <- function(object, ... ) {
  UseMethod(generic = 'Idents', object = object)
}

#' @inheritParams Idents
#' @param value,idents The name of the identites to pull from object metadata or the identities themselves
#'
#' @return \code{Idents<-}: An object with the cell identites changed
#'
#' @rdname Idents
#' @export Idents<-
#'
"Idents<-" <- function(object, ..., value) {
  UseMethod(generic = 'Idents<-', object = object)
}

#' Get JackStraw information
#'
#' @param object An object
#' @param slot Name of slot to store JackStraw scores to
#' Can shorten to 'empirical', 'fake', 'full', or 'overall'
#' @param ... Arguments passed to other methods
#'
#' @rdname JS
#' @export JS
#'
JS <- function(object, slot, ...) {
  UseMethod(generic = 'JS', object = object)
}

#' Set JackStraw information
#'
#' @inherit JS
#' @param value JackStraw information
#'
#' @rdname JS
#' @export JS<-
#'
"JS<-" <- function(object, ..., value) {
  UseMethod(generic = 'JS<-', object = object)
}

#' Get a key
#'
#' @param object An object
#' @param ... Arguments passed to other methods
#'
#' @rdname Key
#' @export Key
#'
Key <- function(object, ...) {
  UseMethod(generic = 'Key', object = object)
}

#' Set a key
#'
#' @inheritParams Key
#' @param value Key value
#'
#' @rdname Key
#' @export Key<-
#'
"Key<-" <- function(object, ..., value) {
  UseMethod(generic = 'Key<-', object = object)
}

#' Get feature loadings
#'
#' @param object An object
#' @param projected Pull the projected feature loadings?
#' @param ... Arguments passed to other methods
#'
#' @rdname Loadings
#' @export Loadings
#'
Loadings <- function(object, projected, ...) {
  UseMethod(generic = 'Loadings', object = object)
}

#' Add feature loadings
#'
#' @inheritParams Loadings
#' @param value Feature loadings to add
#'
#' @rdname Loadings
#' @export Loadings<-
#'
"Loadings<-" <- function(object, ..., value) {
  UseMethod(generic = 'Loadings<-', object = object)
}

#' Access miscellaneous data
#'
#' @param object An object
#' @param slot Name of specific bit of meta data to pull
#' @param ... Arguments passed to other methods
#'
#' @return Miscellaneous data
#'
#' @rdname Misc
#' @export Misc
#'
Misc <- function(object, slot, ...) {
  UseMethod(generic = 'Misc', object = object)
}

#' Set miscellaneous data
#'
#' @inheritParams Misc
#' @param value Data to add
#'
#' @return An object with miscellaneous data added
#'
#' @rdname Misc
#' @export Misc<-
#'
"Misc<-" <- function(object, slot, ..., value) {
  UseMethod(generic = 'Misc<-', object = object)
}

#' Normalize Data
#'
#' Normalize the count data present in a given assay.
#'
#' @param object An object
#' @param normalization.method Method for normalization.
#'  \itemize{
#'   \item{LogNormalize: }{Feature counts for each cell are divided by the total
#'   counts for that cell and multiplied by the scale.factor. This is then
#'   natural-log transformed using log1p.}
#'   \item{CLR: }{Applies a centered log ratio transformation}
#' }
#' More methods to be added.
#' @param scale.factor Sets the scale factor for cell-level normalization
#' @param verbose display progress bar for normalization procedure.
#'
#' @return Returns object after normalization
#'
#' @rdname NormalizeData
#' @export NormalizeData
#'
NormalizeData <- function(
  object,
  normalization.method,
  scale.factor,
  verbose
) {
  UseMethod(generic = 'NormalizeData', object = object)
}

#' Identify cells matching certain criteria
#'
#' Returns a list of cells that match a particular set of criteria such as
#' identity class, high/low values for particular PCs, ect..
#'
#' @param object Seurat object
#' @param cells Subset of cell names
#' @param subset.name Parameter to subset on. Eg, the name of a gene, PC1, a
#' column name in object@@meta.data, etc. Any argument that can be retreived
#' using FetchData
#' @param low.threshold Low cutoff for the parameter (default is -Inf)
#' @param high.threshold High cutoff for the parameter (default is Inf)
#' @param accept.value Returns all cells with the subset name equal to this value
#' @param ... Arguments passed to other methods
# @param \dots Additional arguments to be passed to FetchData
#'
#' @return A vector of cell names
#'
#' @rdname WhichCells
#' @export OldWhichCells
#'
#' @examples
#' OldWhichCells(object = pbmc_small, ident.keep = 2)
#'
OldWhichCells <- function(
  object,
  cells,
  subset.name,
  low.threshold,
  high.threshold,
  accept.value,
  ...
) {
  UseMethod(generic = 'OldWhichCells', object = object)
}
#' Print the results of a dimensional reduction analysis
#'
#' Prints a set of genes that most strongly define a set of components
#'
#' @param object DimReduc object
#'
#' @return Set of features defining the components
#'
#' @rdname Print
#' @export Print
#'
Print <- function(object, ...) {
  UseMethod(generic = "Print", object = object)
}

#' Rename cells
#'
#' Change the cell names in all the different parts of an object. Can
#' be useful before combining multiple objects.
#'
#' @param object An object
#' @param new.names vector of new cell names
#'
#' @details
#' If \code{add.cell.id} is set a prefix is added to existing cell names. If
#' \code{new.names} is set these will be used to replace existing names.
#'
#' @return An object with new cell names
#'
#' @rdname RenameCells
#' @export RenameCells
#'
#' @examples
#' head(x = colnames(x = pbmc_small))
#' pbmc_small <- RenameCells(pbmc_small, add.cell.id = "Test")
#' head(x = colnames(x = pbmc_small))
#'
RenameCells <- function(object, new.names, ...) {
  UseMethod(generic = 'RenameCells', object = object)
}

#' Run Adaptively-thresholded Low Rank Approximation (ALRA)
#'
#' Runs ALRA, a method for imputation of dropped out values in scRNA-seq data.
#' Computes the k-rank approximation to A_norm and adjusts it according to the
#' error distribution learned from the negative values. Described in
#' Linderman, G. C., Zhao, J., Kluger, Y. (2018). "Zero-preserving imputation
#' of scRNA-seq data using low rank approximation." (bioRxiv:138677)
#'
#' @param object a Seurat object
#' @param k  The rank of the rank-k approximation. Set to NULL for automated choice of k.
#' @param q  The number of additional power iterations in randomized SVD when
#' computing rank k approximation. By default, q=10.
#' @rdname RunALRA
#' @export RunALRA
# @import rsvd
#'
#' @author Jun Zhao, George Linderman
#' @references Linderman, G. C., Zhao, J., Kluger, Y. (2018). "Zero-preserving imputation
#' of scRNA-seq data using low rank approximation." (bioRxiv:138677)
#' @seealso \code{\link{ALRAChooseKPlot}}
#'
#' @examples
#' pbmc_small
#' # Example 1: Simple usage, with automatic choice of k.
#' pbmc_small_alra <- RunALRA(object = pbmc_small)
#' \dontrun{
#' # Example 2: Visualize choice of k, then run ALRA
#' # First, choose K
#' pbmc_small_alra <- RunALRA(pbmc_small, k.only=TRUE)
#' # Plot the spectrum, spacings, and p-values which are used to choose k
#' ggouts <- ALRAChooseKPlot(pbmc_small_alra)
#' do.call(gridExtra::grid.arrange, c(ggouts, nrow=1))
#' # Run ALRA with the chosen k
#' pbmc_small_alra <- RunALRA(pbmc_small_alra)
#' }

RunALRA <- function(object, k, q, ...) {
  UseMethod(generic = 'RunALRA', object = object)
}

#' Perform Canonical Correlation Analysis
#'
#' Runs a canonical correlation analysis using a diagonal implementation of CCA.
#' For details about stored CCA calculation parameters, see
#' \code{PrintCCAParams}.
#' @param object1 First Seurat object
#' @param object2 Second Seurat object.
#' @param num.cc Number of canonical vectors to calculate
#'
#' @return Returns a combined Seurat object with the CCA results stored.
#'
#' @seealso \code{MergeSeurat}
#'
#' @export
#' @importFrom irlba irlba
#'
#' @examples
#' pbmc_small
#' # As CCA requires two datasets, we will split our test object into two just for this example
#' pbmc1 <- SubsetData(pbmc_small, cells = colnames(x = pbmc_small)[1:40])
#' pbmc2 <- SubsetData(pbmc_small, cells = colnames(x = pbmc_small)[41:80])
#' pbmc1["group"] <- "group1"
#' pbmc2["group"] <- "group2"
#' pbmc_cca <- RunCCA(object1 = pbmc1, object2 = pbmc2)
#' # Print results
#' Print(object = pbmc_cca[["cca"]])
#'
#' @rdname RunCCA
#' @export RunCCA
#'
RunCCA <- function(object1, object2, num.cc, ...) {
  UseMethod(generic = 'RunCCA', object = object1)
}

#' Perform Canonical Correlation Analysis with more than two groups
#'
#' Runs a canonical correlation analysis
#'
#' @param object.list List of Seurat objects
#' @param niter Number of iterations to perform. Set by default to 25.
#' @param num.ccs Number of canonical vectors to calculate
#' @param standardize standardize scale.data matrices to be centered (mean zero)
#' and scaled to have a standard deviation of 1.
#'
#' @return Returns a combined Seurat object with the CCA stored as a DimReduc
#'
#' @importFrom methods slot
#'
#' @export
#'
#' @examples
#' \dontrun{
#' pbmc_small
#' # As multi-set CCA requires more than two datasets, we will split our test object into
#' # three just for this example
#' pbmc1 <- SubsetData(pbmc_small,cells = pbmc_small@cell.names[1:30])
#' pbmc2 <- SubsetData(pbmc_small,cells = pbmc_small@cell.names[31:60])
#' pbmc3 <- SubsetData(pbmc_small,cells = pbmc_small@cell.names[61:80])
#' pbmc1@meta.data$group <- "group1"
#' pbmc2@meta.data$group <- "group2"
#' pbmc3@meta.data$group <- "group3"
#' pbmc.list <- list(pbmc1, pbmc2, pbmc3)
#' pbmc_cca <- RunMultiCCA(object.list = pbmc.list, genes.use = pbmc_small@var.genes, num.ccs = 3)
#' # Print results
#' PrintDim(pbmc_cca,reduction.type = 'cca')
#' }
#'
RunMultiCCA <- function(
  object.list,
  features,
  add.cell.ids,
  niter,
  num.ccs,
  standardize,
  verbose
) {
  if (length(x = object.list) < 3) {
    stop("Must give at least 3 objects for MultiCCA")
  }
  UseMethod(generic = 'RunMultiCCA', object = object.list[[1]])
}

#' Run Principal Component Analysis on gene expression using IRLBA
#'
#' Run a PCA dimensionality reduction. For details about stored PCA calculation
#' parameters, see \code{PrintPCAParams}.
#'
#' @param object Seurat object
#' @param assay Name of Assay PCA is being run on
#' @param compute.dims Total Number of PCs to compute and store (20 by default)
#' @param rev.pca By default computes the PCA on the cell x gene matrix. Setting
#' to true will compute it on gene x cell matrix.
#' @param weight.by.var Weight the cell embeddings by the variance of each PC
#' (weights the gene loadings if rev.pca is TRUE)
#' @param verbose Print the top genes associated with high/low loadings for
#' the PCs
#' @param print.dims PCs to print genes for
#' @param features.print Number of genes to print for each PC
#' @param reduction.name dimensional reduction name,  pca by default
#' @param reduction.key dimensional reduction key, specifies the string before
#' the number for the dimension names. PC by default
#' @param seed.use Set a random seed. By default, sets the seed to 42. Setting
#' NULL will not set a seed.
#' @param \dots Additional arguments to be passed to IRLBA
#'
#' @importFrom irlba irlba
#' @importFrom methods new
#'
#' @return Returns Seurat object with the PCA calculation stored in the reductions slot
#'
#' @importFrom irlba irlba
#'
#' @export
#'
#' @rdname RunPCA
#' @export RunPCA
#'
RunPCA <- function(
  object,
  pcs.compute,
  rev.pca,
  weight.by.var,
  verbose,
  pcs.print,
  features.print,
  reduction.name,
  reduction.key,
  seed.use,
  ...
) {
  UseMethod(generic = 'RunPCA', object = object)
}

#' Run t-distributed Stochastic Neighbor Embedding
#'
#' Run t-SNE dimensionality reduction on selected features. Has the option of
#' running in a reduced dimensional space (i.e. spectral tSNE, recommended),
#' or running based on a set of genes. For details about stored TSNE calculation
#' parameters, see \code{PrintTSNEParams}.
#'
#' @param object Seurat object
#' @param reduction Which dimensional reduction (e.g. PCA, ICA) to use for
#' the tSNE. Default is PCA
#' @param cells Which cells to analyze (default, all cells)
#' @param seed.use Random seed for the t-SNE
#' @param tsne.method Select the method to use to compute the tSNE. Available
#' methods are:
#' \itemize{
#' \item{Rtsne: }{Use the Rtsne package Barnes-Hut implementation of tSNE (default)}
# \item{tsne: }{standard tsne - not recommended for large datasets}
#' \item{FIt-SNE: }{Use the FFT-accelerated Interpolation-based t-SNE. Based on
#' Kluger Lab code found here: https://github.com/KlugerLab/FIt-SNE}
#' }
#' @param add.iter If an existing tSNE has already been computed, uses the
#' current tSNE to seed the algorithm and then adds additional iterations on top
#' of this
#' @param dim.embed The dimensional space of the resulting tSNE embedding
#' (default is 2). For example, set to 3 for a 3d tSNE
#' @param \dots Additional arguments to the tSNE call. Most commonly used is
#' perplexity (expected number of neighbors default is 30)
#' @param reduction.key dimensional reduction key, specifies the string before the number for the dimension names. tSNE_ by default
#'
#' @rdname RunTSNE
#' @export RunTSNE
#'
RunTSNE <- function(
  object,
  reduction,
  cells,
  seed.use,
  tsne.method,
  add.iter,
  dim.embed,
  reduction.key,
  ...
) {
  UseMethod(generic = 'RunTSNE', object = object)
}

#' Run UMAP
#'
#' Runs the Uniform Manifold Approximation and Projection (UMAP) dimensional
#' reduction technique. To run, you must first install the umap-learn python
#' package (e.g. via pip install umap-learn). Details on this package can be
#' found here: \url{https://github.com/lmcinnes/umap}. For a more in depth
#' discussion of the mathematics underlying UMAP, see the ArXiv paper here:
#' \url{https://arxiv.org/abs/1802.03426}.
#'
#' @param object Seurat object
#' @param cells Which cells to analyze (default, all cells)
#' @param dims Which dimensions to use as input features, used only if
#' \code{genes.use} is NULL
#' @param reduction Which dimensional reduction (PCA or ICA) to use for the
#' UMAP input. Default is PCA
#' @param features If set, run UMAP on this subset of features (instead of running on a
#' set of reduced dimensions). Not set (NULL) by default
#' @param assay Assay to pull data for when using \code{genes.use}
#' @param max.dim Max dimension to keep from UMAP procedure.
#' @param reduction.name dimensional reduction name, specifies the position in
#' the object$dr list. umap by default
#' @param reduction.key dimensional reduction key, specifies the string before
#' the number for the dimension names. UMAP by default
#' @param n_neighbors This determines the number of neighboring points used in
#' local approximations of manifold structure. Larger values will result in more
#' global structure being preserved at the loss of detailed local structure. In
#' general this parameter should often be in the range 5 to 50.
#' @param min_dist min_dist: This controls how tightly the embedding is allowed
#' compress points together. Larger values ensure embedded points are more
#' evenly distributed, while smaller values allow the algorithm to optimise more
#' accurately with regard to local structure. Sensible values are in the range
#' 0.001 to 0.5.
#' @param metric metric: This determines the choice of metric used to measure
#' distance in the input space. A wide variety of metrics are already coded, and
#' a user defined function can be passed as long as it has been JITd by numba.
#' @param seed.use Set a random seed. By default, sets the seed to 42. Setting
#' NULL will not set a seed.
#' @param ... Additional arguments to the umap
#'
#' @return Returns a Seurat object containing a UMAP representation
#'
#' @references McInnes, L, Healy, J, UMAP: Uniform Manifold Approximation and Projection for Dimension Reduction, ArXiv e-prints 1802.03426, 2018
#'
#' @importFrom reticulate import py_module_available py_set_seed
#'
#' @export
#'
#' @examples
#' \dontrun{
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
RunUMAP <- function(
  object,
  cells,
  dims,
  reduction,
  features,
  assay,
  max.dim ,
  reduction.name,
  reduction.key,
  n_neighbors ,
  min_dist,
  metric,
  seed.use,
  ...
) {
  UseMethod(generic = 'RunUMAP', object = object)
}

#' Scale and center the data.
#'
#' Scales and centers features in the dataset. If variables are provided in vars.to.regress,
#' they are individually regressed against each feautre, and the resulting residuals are
#' then scaled and centered.
#'
#' ScaleData now incorporates the functionality of the function formerly known
#' as RegressOut (which regressed out given the effects of provided variables
#' and then scaled the residuals). To make use of the regression functionality,
#' simply pass the variables you want to remove to the vars.to.regress parameter.
#'
#' Setting center to TRUE will center the expression for each feautre by subtracting
#' the average expression for that feautre. Setting scale to TRUE will scale the
#' expression level for each feautre by dividing the centered feautre expression
#' levels by their standard deviations if center is TRUE and by their root mean
#' square otherwise.
#'
#' @param object An object
#' @param features Vector of features names to scale/center. Default is all features
#' @param vars.to.regress Variables to regress out (previously latent.vars in
#' RegressOut). For example, nUMI, or percent.mito.
#' @param model.use Use a linear model or generalized linear model
#' (poisson, negative binomial) for the regression. Options are 'linear'
#' (default), 'poisson', and 'negbinom'
#' @param use.umi Regress on UMI count data. Default is FALSE for linear
#' modeling, but automatically set to TRUE if model.use is 'negbinom' or 'poisson'
#' @param do.scale Whether to scale the data.
#' @param do.center Whether to center the data.
#' @param scale.max Max value to return for scaled data. The default is 10.
#' Setting this can help reduce the effects of feautres that are only expressed in
#' a very small number of cells. If regressing out latent variables and using a
#' non-linear model, the default is 50.
#' @param block.size Default size for number of feautres to scale at in a single
#' computation. Increasing block.size may speed up calculations but at an
#' additional memory cost.
#' @param min.cells.to.block If object contains fewer than this number of cells,
#' don't block for scaling calculations.
#' @param verbose Displays a progress bar for scaling procedure
#'
#' @rdname ScaleData
#' @export ScaleData
#'
ScaleData <- function(
  object,
  features,
  vars.to.regress,
  model.use,
  use.umi,
  do.scale,
  do.center,
  scale.max,
  block.size,
  min.cells.to.block,
  verbose,
  ...
) {
  UseMethod(generic = 'ScaleData', object = object)
}

#' Compute Jackstraw scores significance.
#'
#' Significant PCs should show a p-value distribution that is
#' strongly skewed to the left compared to the null distribution.
#' The p-value for each PC is based on a proportion test comparing the number
#' of genes with a p-value below a particular threshold (score.thresh), compared with the
#' proportion of genes expected under a uniform distribution of p-values.
#'
#' @param object An object
#' @param dims Which dimensions to examine
#' @param score.thresh Threshold to use for the proportion test of PC
#' significance (see Details)
#'
#' @return Returns a Seurat object
#'
#' @author Thanks to Omri Wurtzel for integrating with ggplot
#'
#' @rdname ScoreJackStraw
#' @export ScoreJackStraw
#'
ScoreJackStraw <- function(object, dims, score.thresh, ...) {
  UseMethod(generic = 'ScoreJackStraw', object = object)
}

#' Setter for multimodal data
#'
#' @param object An object
#' @param slot Where to store the new data
#' @param new.data New data to insert
#' @param ... Arguments passed to other methods
#'
#' @return object with the assay data set
#'
#' @rdname SetAssayData
#' @export SetAssayData
#'
SetAssayData <- function(object, slot, new.data, ...) {
  UseMethod(generic = 'SetAssayData', object = object)
}

#' Set cell identity information
#'
#' @inheritParams Idents
#'
#' @return \code{SetIdent}: An object with new identity classes set
#'
#' @rdname Idents
#' @export SetIdent
#'
SetIdent <- function(object, idents, ...) {
  UseMethod(generic = 'SetIdent', object = object)
}

#' Stash an object's identity information
#'
#' @inheritParams Idents
#' @param save.name Store current identity information under this name
#'
#' @return \code{StashIdent}: An object with the identities stashed
#'
#' @rdname Idents
#' @export StashIdent
#'
#' @examples
#' head(x = pbmc_small[])
#' pbmc_small <- StashIdent(object = pbmc_small, save.name = 'cluster.ident')
#' head(x = pbmc_small[])
#'
StashIdent <- function(object, save.name, ...) {
  UseMethod(generic = 'StashIdent', object = object)
}

#' Get the standard deviations for an object
#'
#' @param object An object
#' @param ... Arguments passed to other methods
#'
#' @rdname Stdev
#' @export Stdev
#'
Stdev <- function(object, ...) {
  UseMethod(generic = 'Stdev', object = object)
}

#' Return a subset of the Seurat object
#'
#' Creates a Seurat object containing only a subset of the cells in the
#' original object. Takes either a list of cells to use as a subset, or a
#' parameter (for example, a gene), to subset on.
#'
#' @param object Seurat object
#' @param cells A vector of cell names to use as a subset. If NULL
#' (default), then this list will be computed based on the next three
#' arguments. Otherwise, will return an object consissting only of these cells
#' @param subset.name Parameter to subset on. Eg, the name of a gene, PC1, a
#' column name in object@@meta.data, etc. Any argument that can be retreived
#' using FetchData
#' @param low.threshold Low cutoff for the parameter (default is -Inf)
#' @param high.threshold High cutoff for the parameter (default is Inf)
#' @param accept.value Returns cells with the subset name equal to this value
#' @param do.center Recenter the new object@@scale.data
#' @param do.scale Rescale the new object@@scale.data. FALSE by default
#' @param do.clean Only keep object@@raw.data and object@@data. Cleans out most
#' other slots. Can be useful if you want to start a fresh analysis on just a
#' subset of the data. Also clears out stored clustering results in
#' object@@meta.data (any columns containing "res"). Will by default subset the
#' raw.data slot.
#' @param subset.raw Also subset object@@raw.data
#' @param ... Arguments passed to other methods
# @param \dots Additional arguments to be passed to FetchData (for example,
# use.imputed=TRUE)
#'
#' @return Returns a Seurat object containing only the relevant subset of cells
#'
#' @rdname SubsetData
#' @export SubsetData
#'
#' @examples
#' pbmc1 <- SubsetData(object = pbmc_small, cells = colnames(x = pbmc_small)[1:40])
#' pbmc1
#'
SubsetData <- function(object, slot, ...) {
  UseMethod(generic = 'SubsetData', object = object)
}

#' Get and set variable feature information
#'
#' @param object An object
#' @param ... Arguments passed to other methods
#'
#' @rdname VariableFeatures
#' @export VariableFeatures
#'
VariableFeatures <- function(object, ...) {
  UseMethod(generic = 'VariableFeatures', object = object)
}

#' @inheritParams VariableFeatures
#' @param value A character vector of variable features
#'
#' @rdname VariableFeatures
#' @export VariableFeatures<-
#'
"VariableFeatures<-" <- function(object, ..., value) {
  UseMethod(generic = 'VariableFeatures<-', object = object)
}

#' Identify cells matching certain criteria
#'
#' Returns a list of cells that match a particular set of criteria such as
#' identity class, high/low values for particular PCs, ect..
#'
#' @param object Seurat object
#' @param cells Subset of cell names
#' @param expression A predicate expression for feature/variable expression, can
#' evalue anything that can be pulled by \code{FetchData}
#' @param invert Invert the selection of cells
#' @param ... Arguments passed to other methods
#'
#' @return A vector of cell names
#'
#' @seealso \code{\link{FetchData}}
#' @rdname WhichCells
#' @export WhichCells
#'
#' @examples
#' WhichCells(object = pbmc_small, idents = 2)
#' WhichCells(object = pbmc_small, expression = MS4A1 > 3)
#' WhichCells(object = pbmc_small, idents = c(1, 3), invert = TRUE)
#'
WhichCells <- function(
  object,
  cells,
  expression,
  invert,
  ...
) {
  UseMethod(generic = 'WhichCells', object = object)
}
