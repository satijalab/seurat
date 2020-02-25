# News
All notable changes to Seurat will be documented in this file.
The format is based on [Keep a Changelog](http://keepachangelog.com/en/1.0.0/)

## [3.1.4] - 2020-02-20
### Changes 
- Fixes to `DoHeatmap` to remain compatible with ggplot2 v3.3
- Adoption of `patchwork` framework to replace `CombinePlots`

## [3.1.3] - 2020-02-07
### Added
- New system agnostic `Which` function to address problems with FItSNE on Windows

### Changes
- Export `CellsByIdentities` and `RowMergeSparseMatrices` functions
- nCount and nFeature metadata variables retained after subset and updated properly with `UpdateSeuratObject`
- Fix uwot support for running directly on feature matrices
- Fixes for keys with underscores
- Fix issue with leiden option for `FindClusters`
- Fix for data transfer when using sctransform
- SDMTools moved to Suggests as package is orphaned 

## [3.1.2] - 2019-12-11
### Added
- New silent slot updater
- New random seed options to `RunCCA`, `RunTSNE`, `WhichCells`, `HTODemux`, `AddModuleScore`, `VlnPlot`, and `RidgePlot`
- Enhancements for dealing with `Assay`-derived objects

### Changed
- Only run `CalcN` (generates nFeatures and nCounts) when `counts` changes
- Fix issue regarding colons in feature names
- Change object class testing to use `inherits` or `is.*` for R 4.0 compatability

## [3.1.1] - 2019-09-20
### Added
- New `RegroupIdents` function to reassign idents based on metadata column majority
- `UpdateSymbolList` function to pull new gene names from HGNC
- Added support for H5AD layers as additional assays in a `Seurat` object

### Changed
- Fix rownames issue when running UMAP on dist object
- Add support for new H5AD `obsm` and `varm` stucture
- Fix issue when trying to read non-existent feature-level metadata from an H5AD file
- Fix in integration workflow when using SCTransform
- Improved error checking for `AddModuleScore`
- cbind fix in reference-based integration (`MapQuery`)
- Fix for convenience plots error hanging
- Ensure Seurat objects aren't stored in the command logs

## [3.1.0] - 2019-08-20
### Added
- New `PrepSCTIntegration` function to facilitate integration after `SCTransform`
- Reference-based integration with the `reference` parameter in `FindIntegrationAnchors`
- Reciprocal PCA as a `reduction` option in `FindIntegrationAnchors`
- New `CollapseEmbeddingOutliers` function
- Enable `FindTransferAnchors` after `SCTransform`
- Added back `ColorDimSplit` functionality
- Include a code of conduct
- Added uwot support as new default UMAP method
- Added `CheckDots` to catch unused parameters and suggest updated names
- `Reductions` and `Assays` assays functions to list stored DimReducs and Assays

### Changed
- Fix regex in `LogSeuratCommand`
- Check for NAs in feature names in `Read10X`
- Prevent dimnames for counts/data/scale.data matrices from being arrays
- Updates `ReadH5AD` to distinguish FVF methods
- Fixes to UpdateSeuratObject for v2 objects
- Sink all output from stdout to stderr
- Fix to scale.data cell ordering after subsetting 
- Enable `Assay` specification in `BuildClusterTree`
- Fix `FeaturePlot` when using both `blend` and `split.by`
- Fix to `WhichCells` when passing `cells` and `invert`
- Fix to `HoverLocator` labels and title
- Ensure features names don't contain pipes (`|`)
- Deprecation of `RunLSI` and `RunALRA`
- Fix legend bug when sorting in `ExIPlot`

## [3.0.2] - 2019-06-07
### Added
- Flag to skip singleton grouping in `FindClusters`
- New custom colors for blended `FeaturePlot`s
- New `GetResidual` function
- New Seurat/Monocle converters

### Changed
- Fix issue where certain assays weren't being shown in the `Seurat` object
- Fix issue where we weren't updating `DimReduc` object column names
- Fix line spacers in `DoHeatmap`
- Fix uninformative labels in `FeaturePlot`
- Fix unset identities when converting from SCE to Seurat
- Fix single colors being interpreted as palettes in `SingleDimPlot`
- Ensure factor levels are always numerically increasing after `FindClusters`
- Better cell highlighting colors for `DimPlot`
- Fix to `levels<-.Seurat`
- Add ability to use counts/scaled data in `BuildClusterTree`
- Minor fix to split `ScaleData`

## [3.0.1] - 2019-05-16
### Added
- Add global option (Seurat.memsafe) to skip gc() calls
- Restore draw.lines to DoHeatmap, maintain size of color bar with different number of features (#1429)
- Enable split.by parameter for ScaleData
- Add slot parameter to FeaturePlot (#1483)
- Add assay parameter to DotPlot (#1404)

### Changed
- Fix to color options for VlnPlot with split.by option (#1425)
- Improvements to conversion functions (loom, SCE)
- Fix for cluster tree reordering (#1434)
- Fix PercentageFeatureSet for single feature case
- Fix to fold change calculation and filtering for other slots in FindMarkers (#1454)
- Keep title vectorized in AugmentPlot (#1515)
- Export LogSeuratCommand function
- Fix for FindConservedMarkers when one ident is missing from a group (#1517)

## [3.0.0] - 2019-04-16
### Added
- New method for identifying anchors across single-cell datasets
- Parallelization support via future
- Additional method for demultiplexing with MULTIseqDemux
- Support normalization via sctransform
- New option for clustering with the Leiden algorithm
- Support for reading 10X v3 files
- New function to export Seurat objects for the UCSC cell browser
- Support for data import from Alevin outputs
- Imputation of dropped out values via ALRA

### Changed
- Significant code restructuring
- Most occurances of "gene(s)" in function names/arguments renamed to "feature(s)"
- Changes to the Seurat object class to facilitate multimodal data
- New BlendPlot implementation

## [2.3.4] - 2018-07-13
### Added
- GetIdent function added to pull identity info

### Changed
- DiffusionMap dependency replaced with destiny to avoid archival
- Java dependency removed and functionality rewritten in Rcpp 
- Speed and efficiency improvements for Rcpp code
- More robust duplicate handling in CellCycleScoring

## [2.3.3] - 2018-07-02
### Added
- New HTOHeatmap function
- Support for custom PNG arguments for vector-friendly plotting
- Fix for 'NA'-labeled cells disappearing with custom color scale

### Changed
- Replaced FNN with RANN
- Removed unused compiler flags
- Moved several lightly-used packages from 'imports' to 'suggests'

## [2.3.2] - 2018-06-11
### Added
- RenameCells added for easy renaming of all cells
- Read10X_h5 added to read in 10X formatted h5 files
- SetAssayData ensures cell order is the same between assay objects and the Seurat object
- Compatability updates for ggplot2 v2.3.0

## [2.3.1] - 2018-05-03
### Added
- Support for [UMAP](https://github.com/lmcinnes/umap) dimensional reduction technique
- New conversion functions for SingleCellExperiment and anndata

### Changed
- FetchData preserves cell order
- Require Matrix 1.2-14 or higher
- AddModuleScore no longer densifies sparse-matrices
- Various visualization fixes and improvements
- Default value for latent.vars in FindMarkers/FindAllMarkers changed to NULL.

## [2.3.0] - 2018-03-22
### Added
- Support for HTO demultiplexing
- Utility functions: TransferIdent, CombineIdent, SplitObject, vector.friendly
- C++ implementation for parts of BuildSNN
- Preliminary parallelization support (regression and JackStraw)
- Support for FItSNE

### Changed
- MetaDE replaced with metap for combining p-values (MetaDE was removed from CRAN)
- NMF heatmaps replaced (NMF to be archived by CRAN)

## [2.2.1] - 2018-02-14
### Changed
 - MetaDE replaced with metap for combining p-values (MetaDE was removed from CRAN)
 - NMF heatmaps replaced (NMF to be archived by CRAN)

## [2.2.0] - 2018-01-10
### Added
 - Multiple alignment functionality with RunMultiCCA and AlignSubspace extended to multiple datasets
 - CalcAlignmentScore added to evaluate alignment quality
 - MetageneBicorPlot added to guide CC selection
 - Change cluster order in DoHeatmap with group.order parameter
 - Ability to change plotting order and add a title to DimPlot
 - do.clean and subset.raw options for SubsetData

### Changed
 - JoyPlot has been replaced with RidgePlot
 - FindClusters is now more robust in making temp files
 - MetaDE support for combining p-values in DE testing

## [2.1.0] - 2017-10-12
### Added
- Support for using MAST and DESeq2 packages for differential expression testing in FindMarkers
- Support for multi-modal single-cell data via @assay slot

### Changed
- Default DE test changed to Wilcoxon rank sum test

## [2.0.1] - 2017-08-18
### Added
 - Now available on CRAN
 - Updated documentation complete with examples
 - Example datasets: `pbmc_small` and `cc.genes`
 - C++ implementation for parts of FindVariableGenes
 - Minor bug fixes

## [2.0.0] - 2017-07-26
### Added
- New method for aligning scRNA-seq datasets
- Significant code restructuring
- New methods for scoring gene expression and cell-cycle phases
- New visualization features (do.hover, do.identify)
