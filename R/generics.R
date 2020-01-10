#' Add in metadata associated with either cells or features.
#'
#' Adds additional data to the object. Can be any piece of information
#' associated with a cell (examples include read depth, alignment rate,
#' experimental batch, or subpopulation identity) or feature (ENSG name,
#' variance). To add cell level information, add to the Seurat object. If adding
#' feature-level metadata, add to the Assay object (e.g. object[["RNA"]]))
#'
#' @param x,object An object
#' @param i,col.name Name to store metadata or object as
#' @param value,metadata Metadata or object to add
#' @param j Ignored
#' @param ... Arguments passed to other methods
#'
#' @return An object with metadata or and object added
#'
#' @rdname AddMetaData
#' @export AddMetaData
#'
#' @aliases SeuratAccess
#'
#' @examples
#' cluster_letters <- LETTERS[Idents(object = pbmc_small)]
#' names(cluster_letters) <- colnames(x = pbmc_small)
#' pbmc_small <- AddMetaData(
#'   object = pbmc_small,
#'   metadata = cluster_letters,
#'   col.name = 'letter.idents'
#' )
#' head(x = pbmc_small[[]])
#'
AddMetaData <- function(object, metadata, col.name = NULL) {
  UseMethod(generic = 'AddMetaData', object = object)
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

#' Convert a matrix (or Matrix) to the Graph class.
#'
#' @param x The matrix to convert
#' @param ... Arguments passed to other methods (ignored for now)
#'
#' @rdname as.Graph
#' @export as.Graph
#'
as.Graph <- function(x, ...) {
  UseMethod(generic = "as.Graph", object = x)
}

#' Convert objects to loom objects
#'
#' @param x An object to convert to class \code{loom}
#' @inheritParams loomR::create
#'
#' @seealso \code{\link[loomR]{create}}
#'
#' @rdname as.loom
#' @export as.loom
#'
as.loom <- function(x, ...) {
  UseMethod(generic = 'as.loom', object = x)
}

#' Convert objects to Seurat objects
#'
#' @param x An object to convert to class \code{Seurat}
#' @param ... Arguments passed to other methods
#'
#' @rdname as.Seurat
#' @export as.Seurat
#'
as.Seurat <- function(x, ...) {
  UseMethod(generic = 'as.Seurat', object = x)
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

#' Convert between data frames and sparse matrices
#'
#' @param x An object
#' @param ... Arguments passed to other methods
#'
#' @return \code{as.sparse}: A sparse representation of the input data
#'
#' @rdname as.sparse
#' @export as.sparse
#'
as.sparse <- function(x, ...) {
  UseMethod(generic = 'as.sparse', object = x)
}

#' Get cells present in an object
#'
#' @param x An object
#'
#' @return A vector of cell names
#'
#' @rdname Cells
#' @export Cells
#'
#' @examples
#' Cells(x = pbmc_small)
#'
Cells <- function(x) {
  UseMethod(generic = 'Cells', object = x)
}

#' Get SeuratCommands
#'
#' Pull information on previously run commands in the Seurat object.
#'
#' @param object An object
#' @param ... Arguments passed to other methods
#'
#' @return Either a SeuratCommand object or the requested paramter value
#'
#' @rdname Command
#' @export Command
#'
Command <- function(object, ...) {
  UseMethod(generic = 'Command', object = object)
}

#' Get and set the default assay
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
#' pbmc_small <- BuildClusterTree(object = pbmc_small)
#' markers <- FindMarkers(object = pbmc_small, ident.1 = 'clustertree', ident.2 = 5)
#' head(x = markers)
#'
#' @rdname FindMarkers
#' @export FindMarkers
#'
#' @aliases FindMarkersNode
#'
FindMarkers <- function(object, ...) {
  UseMethod(generic = 'FindMarkers', object = object)
}

#' SNN Graph Construction
#'
#' Constructs a Shared Nearest Neighbor (SNN) Graph for a given dataset. We
#' first determine the k-nearest neighbors of each cell. We use this knn graph
#' to construct the SNN graph by calculating the neighborhood overlap
#' (Jaccard index) between every cell and its k.param nearest neighbors.
#'
#' @param object An object
#' @param ... Arguments passed to other methods
#'
#' @return Returns the object with object@@snn filled
#'
#' @examples
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

#' General accessor function for the Assay class
#'
#' This function can be used to pull information from any of the slots in the Assay class. For
#' example, pull one of the data matrices("counts", "data", or "scale.data").
#'
#' @param object An object
#' @param ... Arguments passed to other methods
#'
#' @return Returns info from requested slot
#'
#' @rdname GetAssayData
#' @export GetAssayData
#'
GetAssayData <- function(object, ...) {
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

#' Get, set, and manipulate an object's identity classes
#'
#' @param x,object An object
#' @param ... Arguments passed to other methods; for \code{RenameIdents}: named
#' arguments as \code{old.ident = new.ident}; for \code{ReorderIdent}: arguments
#' passed on to \code{\link{FetchData}}
#'
#' @return \code{Idents}: The cell identies
#'
#' @rdname Idents
#' @export Idents
#'
#' @examples
#' # Get cell identity classes
#' Idents(object = pbmc_small)
#'
Idents <- function(object, ... ) {
  UseMethod(generic = 'Idents', object = object)
}

#' @inheritParams Idents
#' @param value The name of the identites to pull from object metadata or the identities themselves
#'
#' @return \code{Idents<-}: An object with the cell identites changed
#'
#' @rdname Idents
#' @export Idents<-
#'
#' @examples
#' # Set cell identity classes
#' # Can be used to set identities for specific cells to a new level
#' Idents(object = pbmc_small, cells = 1:4) <- 'a'
#' head(x = Idents(object = pbmc_small))
#'
#' # Can also set idents from a value in object metadata
#' colnames(x = pbmc_small[[]])
#' Idents(object = pbmc_small) <- 'RNA_snn_res.1'
#' levels(x = pbmc_small)
#'
"Idents<-" <- function(object, ..., value) {
  UseMethod(generic = 'Idents<-', object = object)
}

#' Is an object global/persistent?
#'
#' Typically, when removing \code{Assay} objects from an \code{Seurat} object,
#' all associated objects (eg. \code{DimReduc}, \code{Graph}, and \code{SeuratCommand} objects)
#' are removed as well. If an associated object is marked as global/persistent,
#' the associated object will remain even if its original assay was deleted
#'
#' @param object An object
#' @param ... Arguments passed to other methods
#'
#' @return \code{TRUE} if the object is global/persistent otherwise \code{FALSE}
#'
#' @rdname IsGlobal
#' @export IsGlobal
#'
#' @examples
#' IsGlobal(pbmc_small[['pca']])
#'
IsGlobal <- function(object, ...) {
  UseMethod(generic = 'IsGlobal', object = object)
}

#' Get JackStraw information
#'
#' @param object An object
#' @param ... Arguments passed to other methods
#'
#' @rdname JS
#' @export JS
#'
JS <- function(object, ...) {
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
#' @param ... Arguments passed to other methods
#'
#' @rdname Loadings
#' @export Loadings
#'
Loadings <- function(object, ...) {
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
#' @param ... Arguments passed to other methods
#'
#' @return Miscellaneous data
#'
#' @rdname Misc
#' @export Misc
#'
Misc <- function(object, ...) {
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
"Misc<-" <- function(object, ..., value) {
  UseMethod(generic = 'Misc<-', object = object)
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

#' Identify cells matching certain criteria
#'
#' Returns a list of cells that match a particular set of criteria such as
#' identity class, high/low values for particular PCs, ect..
#'
#' @param object An object
#' @param ... Arguments passed to other methods and \code{FetchData}
#'
#' @return A vector of cell names
#'
#' @rdname OldWhichCells
#' @export OldWhichCells
#'
#' @examples
#' \dontrun{
#' OldWhichCells(object = pbmc_small, ident.keep = 2)
#' }
#'
OldWhichCells <- function(object, ...) {
  UseMethod(generic = 'OldWhichCells', object = object)
}

#' Get and set project information
#'
#' @param object An object
#' @param ... Arguments passed to other methods
#'
#' @return Project information
#'
#' @rdname Project
#' @export Project
#'
Project <- function(object, ...) {
  UseMethod(generic = 'Project', object = object)
}

#' @param value Project information to set
#'
#' @return An object with project information added
#'
#' @rdname Project
#' @export Project<-
#'
"Project<-" <- function(object, ..., value) {
  UseMethod(generic = 'Project<-', object = object)
}

#' Read from and write to h5ad files
#'
#' Utilize the Anndata h5ad file format for storing and sharing single-cell expression
#' data. Provided are tools for writing objects to h5ad files, as well as reading
#' h5ad files into a Seurat object
#'
#' @details
#' \code{ReadH5AD} and \code{WriteH5AD} will try to automatically fill slots based
#' on data type and presence. For example, objects will be filled with scaled and
#' normalized data if \code{adata.X} is a dense matrix and \code{raw} is present
#' (when reading), or if the \code{scale.data} slot is filled (when writing). The
#' following is a list of how objects will be filled
#' \describe{
#'   \item{\code{adata.X} is dense and \code{adata.raw} is filled; \code{ScaleData} is filled}{Objects will be filled with scaled and normalized data}
#'   \item{
#'     \code{adata.X} is sparse and \code{adata.raw} is filled; \code{NormalizeData} has been run, \code{ScaleData} has not been run
#'   }{
#'     Objects will be filled with normalized and raw data
#'   }
#'   \item{\code{adata.X} is sparse and \code{adata.raw} is not filled; \code{NormalizeData} has not been run}{Objects will be filled with raw data only}
#' }
#' In addition, dimensional reduction information and nearest-neighbor graphs will
#' be searched for and added if and only if scaled data is being added.
#'
#' When reading: project name is \code{basename(file)}; identity classes will be
#' set as the project name; all cell-level metadata from \code{adata.obs} will be
#' taken; feature level metadata from \code{data.var} and \code{adata.raw.var}
#' (if present) will be merged and stored in assay \code{meta.features}; highly
#' variable features will be set if \code{highly_variable} is present in feature-level
#' metadata; dimensional reduction objects will be given the assay name provided
#' to the function call; graphs will be named \code{assay_method} if method is
#' present, otherwise \code{assay_adata}
#'
#' When writing: only one assay will be written; all dimensional reductions and
#' graphs associated with that assay will be stored, no other reductions or graphs
#' will be written; active identity classes will be stored in \code{adata.obs} as
#' \code{active_ident}
#'
#' @param file Name of h5ad file, or an H5File object for reading in
#'
#' @return \code{ReadH5AD}: A Seurat object with data from the h5ad file
#'
#' @aliases ReadH5AD
#'
#' @rdname h5ad
#' @export ReadH5AD
#'
ReadH5AD <- function(file, ...) {
  UseMethod(generic = 'ReadH5AD', object = file)
}

#' Rename cells
#'
#' Change the cell names in all the different parts of an object. Can
#' be useful before combining multiple objects.
#'
#' @param object An object
#' @param ... Arguments passed to other methods
#'
#' @return An object with new cell names
#'
#' @rdname RenameCells
#' @export RenameCells
#'
RenameCells <- function(object, ...) {
  UseMethod(generic = 'RenameCells', object = object)
}

#' @inheritParams Idents
#'
#' @return \code{RenameIdents}: An object with selected identity classes renamed
#'
#' @rdname Idents
#' @export RenameIdents
#' @aliases RenameIdent
#'
#' @examples
#' # Rename cell identity classes
#' # Can provide an arbitrary amount of idents to rename
#' levels(x = pbmc_small)
#' pbmc_small <- RenameIdents(object = pbmc_small, '0' = 'A', '2' = 'C')
#' levels(x = pbmc_small)
#'
RenameIdents <- function(object, ...) {
  UseMethod(generic = 'RenameIdents', object = object)
}

#' @inheritParams Idents
#' @param var Feature or variable to order on
#'
#' @return \code{ReorderIdent}: An object with
#'
#' @rdname Idents
#' @export ReorderIdent
#' @aliases ReorderIdent
#'
#' @examples
#' \dontrun{
#' head(x = Idents(object = pbmc_small))
#' pbmc_small <- ReorderIdent(object = pbmc_small, var = 'PC_1')
#' head(x = Idents(object = pbmc_small))
#' }
#'
ReorderIdent <- function(object, var, ...) {
  UseMethod(generic = 'ReorderIdent', object = object)
}

#' Run Adaptively-thresholded Low Rank Approximation (ALRA)
#'
#' Runs ALRA, a method for imputation of dropped out values in scRNA-seq data.
#' Computes the k-rank approximation to A_norm and adjusts it according to the
#' error distribution learned from the negative values. Described in
#' Linderman, G. C., Zhao, J., Kluger, Y. (2018). "Zero-preserving imputation
#' of scRNA-seq data using low rank approximation." (bioRxiv:138677)
#'
#' @note RunALRA and associated functions are being moved to SeuratWrappers;
#' for more information on SeuratWrappers, please see \url{https://github.com/satijalab/seurat-wrappers}
#'
#' @param object An object
#' @param ... Arguments passed to other methods
#'
#' @rdname RunALRA
#' @export RunALRA
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
#'
RunALRA <- function(object, ...) {
  .Deprecated(
    new = 'SeruatWrappers::RunALRA',
    msg = paste(
      'RunALRA and associated functions are being moved to SeuratWrappers;',
      'for more information on SeuratWrappers, please see https://github.com/satijalab/seurat-wrappers'
    )
  )
  UseMethod(generic = 'RunALRA', object = object)
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

#' Run Latent Semantic Indexing on binary count matrix
#'
#' For details about stored LSI calculation parameters, see
#' \code{PrintLSIParams}.
#'
#' @note RunLSI is being moved to Signac. Equivalent functionality can be
#' achieved via the Signac::RunTFIDF and Signac::RunSVD functions;
#' for more information on Signac, please see
#' \url{https://github.com/timoast/Signac}
#'
#' @param object Seurat object
#' @param ... Arguments passed to other methods
#'
#' @rdname RunLSI
#' @export RunLSI
#'
RunLSI <- function(object, ...) {
  .Deprecated(
    new = 'Signac::RunTFIDF',
    msg = paste(
      "RunLSI is being moved to Signac. Equivalent functionality can be",
      "achieved via the Signac::RunTFIDF and Signac::RunSVD functions; for",
      "more information on Signac, please see https://github.com/timoast/Signac"
    )
  )
  UseMethod(generic = "RunLSI", object = object)
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
#' reduction technique. To run, you must first install the umap-learn python
#' package (e.g. via \code{pip install umap-learn}). Details on this package can be
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
#' @param ... Arguments passed to other methods
#'
#' @rdname ScaleData
#' @export ScaleData
#'
ScaleData <- function(object, ...) {
  UseMethod(generic = 'ScaleData', object = object)
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

#' Setter for multimodal data
#'
#' @param object An object
#' @param ... Arguments passed to other methods
#'
#' @return object with the assay data set
#'
#' @rdname SetAssayData
#' @export SetAssayData
#'
SetAssayData <- function(object, ...) {
  UseMethod(generic = 'SetAssayData', object = object)
}

#' @return \code{SetIdent}: An object with new identity classes set
#'
#' @rdname Idents
#' @export SetIdent
#'
#' @examples
#' # Set cell identity classes using SetIdent
#' cells.use <- WhichCells(object = pbmc_small, idents = '1')
#' pbmc_small <- SetIdent(object = pbmc_small, cells = cells.use, value = 'B')
#'
SetIdent <- function(object, ...) {
  UseMethod(generic = 'SetIdent', object = object)
}

#' @return \code{StashIdent}: An object with the identities stashed
#'
#' @rdname Idents
#' @export StashIdent
#'
#' @examples
#' head(x = pbmc_small[[]])
#' pbmc_small <- StashIdent(object = pbmc_small, save.name = 'idents')
#' head(x = pbmc_small[[]])
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
#' @param object An object
#' @param ... Arguments passed to other methods
#'
#' @return Returns a Seurat object containing only the relevant subset of cells
#'
#' @rdname SubsetData
#' @export SubsetData
#'
#' @examples
#' \dontrun{
#' pbmc1 <- SubsetData(object = pbmc_small, cells = colnames(x = pbmc_small)[1:40])
#' pbmc1
#' }
#'
SubsetData <- function(object, ...) {
  UseMethod(generic = 'SubsetData', object = object)
}

#' Get and set additional tool data
#'
#' Use \code{Tool} to get tool data. If no additional arguments are provided,
#' will return a vector with the names of tools in the object.
#'
#' @param object An object
#' @param ... Arguments passed to other methods
#'
#' @return If no additional arguments, returns the names of the tools in the
#' object; otherwise returns the data placed by the tool requested
#'
#'@note For developers: set tool data using \code{Tool<-}. \code{Tool<-} will
#'automatically set the name of the tool to the function that called \code{Tool<-},
#'so each function gets one entry in the tools list and cannot overwrite another
#'function's entry. The automatic naming will also remove any method identifiers
#'(eg. RunPCA.Seurat will become RunPCA); please plan accordingly.
#'
#' @rdname Tool
#' @export Tool
#' @aliases Tools
#'
Tool <- function(object, ...) {
  UseMethod(generic = 'Tool', object = object)
}

#' @inheritParams Tool
#' @param value Information to be added to tool list
#'
#' @rdname Tool
#' @export Tool<-
#'
"Tool<-" <- function(object, ..., value) {
  UseMethod(generic = 'Tool<-', object = object)
}

#' Get and set variable feature information
#'
#' @param object An object
#' @param selection.method Method used to set variable features
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
#' @param object An object
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
#' levels(x = pbmc_small)
#' WhichCells(object = pbmc_small, idents = c(1, 2), invert = TRUE)
#'
WhichCells <- function(object, ...) {
  UseMethod(generic = 'WhichCells', object = object)
}

#' @param object An object
#' @param ... arguments passed to other methods
#'
#' @return \code{WriteH5AD}: None, writes to disk
#'
#' @rdname h5ad
#'
WriteH5AD <- function(object, ...) {
  UseMethod(generic = 'WriteH5AD', object = object)
}
