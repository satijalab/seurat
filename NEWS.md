# Seurat 4.4.0 (2023-09-27)

## Added
- Added parallelization support with speed improvements for `PrepSCTFindMarkers` 
- Fix bug in `LoadNanostring`([#7566](https://github.com/satijalab/seurat/pull/7566))

## Changes
- Fix bug in `as.Seurat.SingleCellExperiment()` ([#6692](https://github.com/satijalab/seurat/issues/6692))
- Support for Visium probe information introduced in Spaceranger 2.1 ([#7141](https://github.com/satijalab/seurat/pull/7141))
- Add `LoadCurioSeeker` to load sequencing-based spatial datasets generated using the Curio Seeker 
- Fix fold change calculation for assays ([#7095](https://github.com/satijalab/seurat/issues/7095))
- Fix `pt.size` bug when rasterization is set to true ([#7379](https://github.com/satijalab/seurat/issues/7379)) 
- Fix `FoldChange` and `FindMarkers` to support all normalization approaches ([#7115](https://github.com/satijalab/seurat/pull/7115),[#7110](https://github.com/satijalab/seurat/issues/7110),[#7095](https://github.com/satijalab/seurat/issues/7095),[#6976](https://github.com/satijalab/seurat/issues/6976),[#6654](https://github.com/satijalab/seurat/issues/6654),[#6701](https://github.com/satijalab/seurat/issues/6701),[#6773](https://github.com/satijalab/seurat/issues/6773), [#7107](https://github.com/satijalab/seurat/issues/7107))
- Fix for handling newer ParseBio formats in `ReadParseBio` ([#7565](https://github.com/satijalab/seurat/pull/7565))
- Fix for handling rasterization by default ([#7842](https://github.com/satijalab/seurat/pull/7842))

# Seurat 4.3.0 (2022-11-18)

## Added
- Add support for imaging-based spatial datasets

## Changes
- Fix bug in `FindMarkers()` when run post Integration/Transfer ([#6856](https://github.com/satijalab/seurat/issues/6586))

# Seurat 4.2.1 (2022-11-08)

## Changes
- Replaced import from `spatstat.core` with `spatstat.explore`
- Fix bug in `FindMarkers()` when using `SCT` assay ([#6856](https://github.com/satijalab/seurat/issues/6586))

# Seurat 4.2.0 (2022-09-21)

## Changes
- Fix legend color in `DoHeatmap()` ([#5783](https://github.com/satijalab/seurat/issues/5783))
- Fix bug in `ScaleData()` when regressing out one gene ([#5970](https://github.com/satijalab/seurat/pull/5970))
- Fix name pulling in `PlotPerturbScore()` ([#6081](https://github.com/satijalab/seurat/pull/6081))
- Support spaceranger 2.0 ([#6208](https://github.com/satijalab/seurat/pull/6208))
- Fix bug in `SpatialDimPlot()` when using `group.by` ([#6179](https://github.com/satijalab/seurat/issues/6179))
- Add `add.noise` parameter in `VlnPlot()`
([#5756](https://github.com/satijalab/seurat/issues/5756))
- Fix uwot model backwards compatibility ([#6345](https://github.com/satijalab/seurat/issues/6345))
- Allow `pseudocount.use` in differential expression functions to be set at the `Assay` level

# Seurat 4.1.1 (2022-05-01)

## Changes
- Fix `giveCsparse` related warnings in `Read10X_h5`
- Fix ident labeling for `SpatialPlot` ([#5774](https://github.com/satijalab/seurat/issues/5774))
- Fix `ReadMtx` on Windows ([#5687](https://github.com/satijalab/seurat/issues/5687))
- Fix `VlnPlot` to switch on rasterization only when required ([#5846](https://github.com/satijalab/seurat/pull/5846))
- Fix `ncol` behavior in `SpatialPlot` ([#5774](https://github.com/satijalab/seurat/issues/5774))
- Set `jitter` to FALSE in `FeatureScatter` ([#5876](https://github.com/satijalab/seurat/pull/5876))
- Update `Cells` methods to new signature (`x, ...`)
- Replace use of `default.stringsAsFactors()` with `getOption("stringsAsFactors")`

# Seurat 4.1.0 (2022-01-14)
## Added
- Add `raster.dpi` parameter to `DimPlot/FeaturePlot` to optionally rasterize individual points ([#5392](https://github.com/satijalab/seurat/pull/5392))
- Add support for sctransform v2, differential expression with SCT assay

## Changes
- Update `ReadParseBio` to support split-pipe 0.9.6p ([#5446](https://github.com/satijalab/seurat/pull/5446))
- Fixes for MAST differential expression ([#5441](https://github.com/satijalab/seurat/issues/5441))
- Fix scaling options when using `split.by` in `FeaturePlot()` ([#5243](https://github.com/satijalab/seurat/issues/5243)) 

# Seurat 4.0.6 (2021-12-16)
## Added

- Implement supervised LSI

## Changes
- Add `raster` parameter to `VlnPlot` to optionally rasterize individual points ([#5076](https://github.com/satijalab/seurat/pull/5076))
- Add `min.cells.group` parameter to `FindConservedMarkers` ([#5079](https://github.com/satijalab/seurat/pull/5079))
- Set `do.center` to FALSE for `lsiproject` in `FindTransferAnchors`
- Fix error message in `ReadMtx()` ([#5158](https://github.com/satijalab/seurat/issues/5158))
- Add `label.color` parameter to `FeaturePlot` ([#5314](https://github.com/satijalab/seurat/pull/5314))
- Fix issues in `ProjectUMAP` ([#5257](https://github.com/satijalab/seurat/issues/5257), [#5104](https://github.com/satijalab/seurat/issues/5104), [#5373](https://github.com/satijalab/seurat/issues/5373))

# Seurat 4.0.5 (2021-10-04)
## Changes
- Update documentation for `to.upper` parameter in `Load10X_Spatial()` ([#4576](https://github.com/satijalab/seurat/issues/4576))
- Update concept tags for `RunSPCA()` ([#4978](https://github.com/satijalab/seurat/discussions/4987))
- Conditionally run tests/packages that use suggested packages ([#5160](https://github.com/satijalab/seurat/pull/5160))
- Set random state in `RunUMAP()` when using the `umap-learn` method ([#5194](https://github.com/satijalab/seurat/issues/5194))

# Seurat 4.0.4 (2021-08-19)
## Added
- Add `reduction` parameter to `BuildClusterTree()` ([#4598](https://github.com/satijalab/seurat/issues/4598))
- Add DensMAP option to `RunUMAP()` ([#4630](https://github.com/satijalab/seurat/pull/4630))
- Add `image` parameter to `Load10X_Spatial()` and `image.name` parameter to `Read10X_Image()` ([#4641](https://github.com/satijalab/seurat/pull/4641))
- Add `ReadSTARsolo()` function to read output from STARsolo
- Add `densify` parameter to `FindMarkers()`
- Add `ReadParsebio()` function to read output from Parse Biosciences
- Add the `image.alpha` parameter to `SpatialDimPlot()` and `SpatialFeaturePlot()`
- Add support for the correlation metric in `RunUMAP` ([#4972](https://github.com/satijalab/seurat/issues/4972))


## Changes
- Warn and continue rather than erroring if not all features are available in `FindSpatiallyVariableFeatures()` ([#4611](https://github.com/satijalab/seurat/issues/4611))
- Bug fix for SCT-based integration in selecting proper reference model ([#4355](https://github.com/satijalab/seurat/issues/4355))
- Bug fix for reading from absolute paths in ReadMtx ([#4723](https://github.com/satijalab/seurat/issues/4723))
- Bug fix in SingleCellExperiment conversion ([#4633](https://github.com/satijalab/seurat/issues/4633))
- Bug fix in `FindVariableFeatures()` when using `selection.method = "mvp"` and `binning.method = "equal_frequency"` ([#4712](https://github.com/satijalab/seurat/issues/4712))
- Bug fix in `DoHeatmap()` to remove random characters from plot legend([#4660](https://github.com/satijalab/seurat/issues/4660))
- Fix cell renaming in `RunCCA()`
- Fix issue in SingleCellExperiment conversion where the mainExp would not be set properly
- Fix for default dispersion info displayed in `VariableFeaturePlot()`

# Seurat 4.0.3 (2021-06-10)
## Added
- Add `jitter` parameter to `FeatureScatter()`

## Changes
- Fix issues with `as.SingleCellExperiment.Seurat()` for the latest verion of SingleCellExperiment ([#4532](https://github.com/satijalab/seurat/pull/4532))
- Ensure proper reference.reduction is used in `MapQuery()`
- Fix to `UpdateSymbolList()`, no longer searches aliases and exposes the `search.types` parameter in `GeneSymbolThesarus()` ([#4545](https://github.com/satijalab/seurat/issues/4545))
- Transfer `scale.data` slot as well when converting with `as.SingleCellExperiment.Seurat()`
- Enable `alpha` parameter for `SpatialDimPlot()`
- Fix `as.SingleCellExperiment.Seurat()` conversion for atypical `reducedDim` components

# Seurat 4.0.2 (2020-05-20)
## Added
- New `AddAzimuthScores()` and `AddAzimuthResults()` functions
- Add `shuffle` parameter to `FeatureScatter()` ([#4280](https://github.com/satijalab/seurat/pull/4280))
- Add `lsiproject` and `rpca` options for `FindTransferAnchors()`
- Add `rlsi` option for `FindIntegrationAnchors()`

## Changes
- Preserve feature metadata when converting from `SingleCellExperiment` to `SeuratObject` class
([#4205](https://github.com/satijalab/seurat/issues/4205))
- Preserve multiple assays when converting from `SingleCellExperiment` to `SeuratObject` class
([#3764](https://github.com/satijalab/seurat/issues/3764))
- Fix passing of `score.thresh` parameter in `ScoreJackStraw()` ([#4268](https://github.com/satijalab/seurat/pull/4268))
- Fix FC calculation in `FindMarkers()` non-log transformed data.
- Add umap-learn version >= 0.5.0 compatibility for `RunUMAP()`
- Fix `DotPlot` to use `log1p` when `scale=False`
([#4298](https://github.com/satijalab/seurat/issues/4298))
- Fix split and shuffled `DimPlot`
- Disallow NULL or another length 0 vector for `ident.1` in `FindMarkers()`
- Fix range shift when labeling clusters on a GeomSpatial plot
- Fix SpatialPlot distortion for non-square images.
- Fix future-related warnings in `FindIntegrationAnchors()`
- Fix `fc.name` parameter in `FindMarkers()` ([#4474](https://github.com/satijalab/seurat/issues/4474))
- Deprecate `group.by` parameter in `PlotPerturbScore()` in favor of `mixscape.class`.

# Seurat 4.0.1 (2020-03-17)
## Added
- Add direction option to `PlotClusterTree()`
- Add `cols` parameter to `JackStrawPlot()`
- Add `ReadMtx()` to read local and remote mtx files with associated cell and feature name files

## Changes
- Equality added to differential expression thresholds in `FindMarkers` (e.g, >= logfc.threshold rather than >)
- `Read10X()` now prepends dataset number for first dataset when reading multiple datasets
- Bug fix for `subset.AnchorSet()`
- Bug fix for fold change values in `FindMarkers()` when setting a different pseudocount ([#4111](https://github.com/satijalab/seurat/pull/4111))
- Bug fix for `RunLDA()` related to proper passing of assay parameter.
- When using `order=TRUE` in `SingleDimPlot()`, print NA points under all others.
- Remove default parameter value for `data.dir` in `Read10X()`
- Import spatstat fxns from subpackages (spatstat.core, spatstat.geom)
- `RunUMAP` now checks for graph/neighbor consistency

# Seurat 4.0.0 (2020-01-27)
## Added
- Expose `FoldChange()` component in `FindMarkers()`.
- Add the `merge.DimReduc` method
- Add `IntegrateEmbeddings()` to correct embeddings of `DimReduc`s
- Add `ProjectUMAP()` to project query cells into a reference UMAP space
- Add `MapQuery()` as a wrapper around `IntegrateData()`, `IntegrateEmbeddings()`, and `ProjectUMAP()`
- Add `MappingScore` to compute a per-cell mapping score used in Azimuth
- Add `AggregateExpression()` for summation based pseudobulk calculations
- Add mixscape functionality via `CalcPerturbSig()`, `PrepLDA()`, `RunLDA()`, `DEenrichRPlot()`, `MixscapeHeatmap()`, `MixscapeLDA()`, `PlotPerturbScore()`, `RunMixscape()`
- Add `FindSubCluster()` to further cluster existing clusters
- Add supervised PCA functionality via `RunSPCA()`
- Add functionality to enable weighted nearest neighbor analyses via `FindMultiModalNeighbors()`
- Add neighbor visualization plot via `NNPlot()`.
- Add `PredictAssay()` to impute expression or embeddings from nearest neighbors
- Add `Graphs()` function to access the names of the stored Graph objects or pull a specific one
- Add checks for NA, NaN, logical, non-integer, and infinite values during CreateAssayObject and NormalizeData.default
- Add `AnnotateAnchors()` to aid in AnchorSet interpretation as well as `subset.AnchorSet()`
- Add flexibility of choice for cell column in `Read10X()`
- Add rasterization option to `FeatureScatter()` and `VariableFeaturePlot()`
- Add step1 feature parameters in the SCTModel via `PrepVSTResults()`

## Changes
- Default neighbor finding algorithm changed from "rann" to "annoy"
- Default `ncells` parameter in `SCTransform()` changed to 5000
- Default fold change in `FindMarkers()` changed from ln to log2
- Implementation improvements to `AverageExpression()`
- `AnchorSet` class re-implemented as a virtual class from which `IntegrationAnchorSet` and `TransferAnchorSet` now inherit.
- Point size in `VlnPlot()` now set automatically if not specified
- Return the sample.tree properly when integrating with a single reference dataset
- Replace `as.character.quosure` usage with `as_label` due to deprecation
- Minor changes to the exact calculation of the anchor weight matrix
- Default rasterization limit in `DimPlot()` and `FeaturePlot()` changed from 50,000 to 100,000
- `SCTransform()` now returns a formalized `Assay` subclass `SCTAssay()`
- When using `normalization.method='SCT'` in `FindTransferAnchors()`, normalize query using reference SCT model when possible.
- Change default Neighbor name in `FindNeighbors` to `Assay.nn`

## Removed
- `CreateGeneActivityMatrix` replaced by `Signac::GeneActivity()`
- `RunLSI` replaced by by `Signac::RunTFIDF()` and `Signac::RunSVD()`
- `ReadAlevin` and `ReadAlevinCsv` moved to SeuratWrappers
- `ExportToCellbrowser` and `StopCellbrowser` moved to SeuratWrappers

# Seurat 3.2.3 - 2020-12-14
## Added
- Titles added to `DimPlot` when specifying `group.by` parameter
- `keep.scale` parameter added to `FeaturePlot` to control scaling across multiple features and/or splits.

## Changes
- `Same` deprecated in favor of `base::identity`
- Fix in `DietSeurat` to work with specialized `Assay` objects
- Fix p-value return when using the `ape` implementation of Moran's I
- Fix bug in FindMarkers when using MAST with a latent variable
- Updates to `Key<-.DimReduc` that allow handling of empty reduction column names
- Allow setting `ctrl` in `CellCycleScoring`
- Modify subset.Seurat to allow specialized Assay subsetting methods
- Fix image selection in interactive spatial plots
- Update Rcpp functions with `export(rng=FALSE)` to avoid potential future warnings
- Fix RenameCells bug for integrated SCT assays
- Fix highlight order with proper factor levels when using `SetHighlight` in plots
- Small change in CellRanger version detection logic of h5 file to improve robustness to outside tools.
- `do.cpp` deprecated and will default to true

# Seurat 3.2.2 (2020-09-25)
## Changes
- Set the seed in `WhichCells` regardless of whether or not `idents` is passed
- Retain Graph and Neighbor objects when subsetting only on features
- Fix data.frame input to `CreateAssayObject()` when data.frame has no rownames.
- Default annoy search to sequential if not using multicore future plans.
- Require sctransform >= 0.3.0

# Seurat 3.2.1 (2020-09-04)
## Added
- Added support for nearest neighbor input and `return.model` parameter in `RunUMAP()`
- Enable named color vectors in `DoHeatmap()`
- Add `label.color` and `label.box` parameters to `DimPlot`
- Added `shuffle` and `seed` parameters to `DimPlot()` to help with overplotting
- Added new stacked violin plot functionality

## Changes
- Allow setting `slot` parameter in `RunUMAP`
- Added support for FIt-SNE v1.2+
- Fix for `Spatial*Plot` when running with interactive=TRUE
- Set max for number of items returned by `Top` and remove duplicate items when balanced=TRUE
- Fix logging bug when functions were run via `do.call()`
- Fix handling of weight.by.var parameter when approx=FALSE in `RunPCA()`
- Fix issue where feature names with dashes crashed `CellSelector`
- Fix issue where errors in subsetting were being swallowed
- Fix issue where labeling uncropped spatial plots was broken

## Deprecated
- `CreateActivityMatrix` deprecated in favor of `Signac::GeneActivity`
- `ReadAlevin` and `ReadAlevinCsv` deprecated in favor of `SeuratWrappers::ReadAlevin`
- `ExportToCellbrowser` and `StopCellbrowser` deprecated in favor of `SeuratWrappers::ExportToCellbrowser` and `SeuratWrappers::StopCellbrowser`
- `ReadH5AD` and `WriteH5AD` deprecated in favor of h5Seurat/H5AD functionality found in SeuratDisk
- `as.loom` and `as.Seurat.loom` deprecated in favor of functionality found in SeuratDisk

# Seurat 3.2.0 (2020-07-15)
## Added
- Added ability to create a Seurat object from an existing Assay object, or any
object inheriting from the Assay class
- Added ability to cluster idents and group features in `DotPlot`
- Added ability to use RColorBrewer plaettes for split `DotPlots`
- Added visualization and analysis functionality for spatially resolved datasets (Visium, Slide-seq).

## Changes
- Removed `add.iter` parameter from `RunTSNE` function
- Fixed integer overflow error in the WilcoxDETest function
- Minor visual fixes in `DoHeatmap` group bar + labels
- Efficiency improvements in anchor scoring (`ScoreAnchors`)
- Fix bug in `FindClusters()` when the last node has no edges
- Default to weighted = TRUE when constructing igraph objects in `RunLeiden`. Remove corresponding weights parameter from `FindClusters()`.
- Fix handling of keys in `FeatureScatter()`
- Change `CellSelector` to use Shiny gadgets instead of SDMTools
- Mark `PointLocator` as defunct
- Remove `SDMTools`
- Fixed data slot return in `AverageExpression` when subsetting features and returning a Seurat object

# Seurat 3.1.5 (2020-04-14)
## Added
- New `scale` parameter in `DotPlot`
- New `keep.sparse parameter in `CreateGeneActivityMatrix` for a more memory efficient option
- Added ability to store model learned by UMAP and project new data
- New `strip.suffix` option in `Read10X`. **This changes the default behavior of `Read10X`**.
  A trailing `-1` present in all cell names will not be removed by default.
- Added `group.by` parameter to `FeatureScatter`

## Changes
- Replace wilcox.test with limma implementation for a faster FindMarkers default method
- Better point separation for `VlnPlot`s when using the `split.by` option
- Efficiency improvements for anchor pairing
- Deprecate redundant `sort.cell` parameter in `FeaturePlot`
- Fixes to ensure correct class of Matrix passed to c++ functions
- Fixes for underscores in ident labels for `DotPlot`
- Ensure preservation of matrix dimnames in `SampleUMI`
- Fix non-standard evaluation problems in `subset` and `WhichCells`
- Default split violin option is now a multi group option
- Preserve alpha in `FeaturePlot` when using `blend`
- Update `assay.used` slot for `DimReduc`s when Assay is renamed

# Seurat 3.1.4 (2020-02-20)
## Changes
- Fixes to `DoHeatmap` to remain compatible with ggplot2 v3.3
- Adoption of `patchwork` framework to replace `CombinePlots`

# Seurat 3.1.3 (2020-02-07)
## Added
- New system agnostic `Which` function to address problems with FItSNE on Windows

## Changes
- Export `CellsByIdentities` and `RowMergeSparseMatrices` functions
- nCount and nFeature metadata variables retained after subset and updated properly with `UpdateSeuratObject`
- Fix uwot support for running directly on feature matrices
- Fixes for keys with underscores
- Fix issue with leiden option for `FindClusters`
- Fix for data transfer when using sctransform
- SDMTools moved to Suggests as package is orphaned

# Seurat 3.1.2 (2019-12-11)
## Added
- New silent slot updater
- New random seed options to `RunCCA`, `RunTSNE`, `WhichCells`, `HTODemux`, `AddModuleScore`, `VlnPlot`, and `RidgePlot`
- Enhancements for dealing with `Assay`-derived objects

## Changed
- Only run `CalcN` (generates nFeatures and nCounts) when `counts` changes
- Fix issue regarding colons in feature names
- Change object class testing to use `inherits` or `is.*` for R 4.0 compatability

# Seurat 3.1.1 (2019-09-20)
## Added
- New `RegroupIdents` function to reassign idents based on metadata column majority
- `UpdateSymbolList` function to pull new gene names from HGNC
- Added support for H5AD layers as additional assays in a `Seurat` object

## Changed
- Fix rownames issue when running UMAP on dist object
- Add support for new H5AD `obsm` and `varm` stucture
- Fix issue when trying to read non-existent feature-level metadata from an H5AD file
- Fix in integration workflow when using SCTransform
- Improved error checking for `AddModuleScore`
- cbind fix in reference-based integration (`MapQuery`)
- Fix for convenience plots error hanging
- Ensure Seurat objects aren't stored in the command logs

# Seurat 3.1.0 (2019-08-20)
## Added
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

## Changed
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

# Seurat 3.0.2 (2019-06-07)
## Added
- Flag to skip singleton grouping in `FindClusters`
- New custom colors for blended `FeaturePlot`s
- New `GetResidual` function
- New Seurat/Monocle converters

## Changed
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

# Seurat 3.0.1 (2019-05-16)
## Added
- Add global option (Seurat.memsafe) to skip gc() calls
- Restore draw.lines to DoHeatmap, maintain size of color bar with different number of features (#1429)
- Enable split.by parameter for ScaleData
- Add slot parameter to FeaturePlot (#1483)
- Add assay parameter to DotPlot (#1404)

## Changed
- Fix to color options for VlnPlot with split.by option (#1425)
- Improvements to conversion functions (loom, SCE)
- Fix for cluster tree reordering (#1434)
- Fix PercentageFeatureSet for single feature case
- Fix to fold change calculation and filtering for other slots in FindMarkers (#1454)
- Keep title vectorized in AugmentPlot (#1515)
- Export LogSeuratCommand function
- Fix for FindConservedMarkers when one ident is missing from a group (#1517)

# Seurat 3.0.0 (2019-04-16)
## Added
- New method for identifying anchors across single-cell datasets
- Parallelization support via future
- Additional method for demultiplexing with MULTIseqDemux
- Support normalization via sctransform
- New option for clustering with the Leiden algorithm
- Support for reading 10X v3 files
- New function to export Seurat objects for the UCSC cell browser
- Support for data import from Alevin outputs
- Imputation of dropped out values via ALRA

## Changed
- Significant code restructuring
- Most occurances of "gene(s)" in function names/arguments renamed to "feature(s)"
- Changes to the Seurat object class to facilitate multimodal data
- New BlendPlot implementation

# Seurat 2.3.4 (2018-07-13)
## Added
- GetIdent function added to pull identity info

## Changed
- DiffusionMap dependency replaced with destiny to avoid archival
- Java dependency removed and functionality rewritten in Rcpp
- Speed and efficiency improvements for Rcpp code
- More robust duplicate handling in CellCycleScoring

# Seurat 2.3.3 (2018-07-02)
## Added
- New HTOHeatmap function
- Support for custom PNG arguments for vector-friendly plotting
- Fix for 'NA'-labeled cells disappearing with custom color scale

## Changed
- Replaced FNN with RANN
- Removed unused compiler flags
- Moved several lightly-used packages from 'imports' to 'suggests'

# Seurat 2.3.2 (2018-06-11)
## Added
- RenameCells added for easy renaming of all cells
- Read10X_h5 added to read in 10X formatted h5 files
- SetAssayData ensures cell order is the same between assay objects and the Seurat object
- Compatability updates for ggplot2 v2.3.0

# Seurat 2.3.1 (2018-05-03)
## Added
- Support for [UMAP](https://github.com/lmcinnes/umap) dimensional reduction technique
- New conversion functions for SingleCellExperiment and anndata

## Changed
- FetchData preserves cell order
- Require Matrix 1.2-14 or higher
- AddModuleScore no longer densifies sparse-matrices
- Various visualization fixes and improvements
- Default value for latent.vars in FindMarkers/FindAllMarkers changed to NULL.

# Seurat 2.3.0 (2018-03-22)
## Added
- Support for HTO demultiplexing
- Utility functions: TransferIdent, CombineIdent, SplitObject, vector.friendly
- C++ implementation for parts of BuildSNN
- Preliminary parallelization support (regression and JackStraw)
- Support for FItSNE

## Changed
- MetaDE replaced with metap for combining p-values (MetaDE was removed from CRAN)
- NMF heatmaps replaced (NMF to be archived by CRAN)

# Seurat 2.2.1 (2018-02-14)
## Changed
 - MetaDE replaced with metap for combining p-values (MetaDE was removed from CRAN)
 - NMF heatmaps replaced (NMF to be archived by CRAN)

# Seurat 2.2.0 (2018-01-10)
## Added
 - Multiple alignment functionality with RunMultiCCA and AlignSubspace extended to multiple datasets
 - CalcAlignmentScore added to evaluate alignment quality
 - MetageneBicorPlot added to guide CC selection
 - Change cluster order in DoHeatmap with group.order parameter
 - Ability to change plotting order and add a title to DimPlot
 - do.clean and subset.raw options for SubsetData

## Changed
 - JoyPlot has been replaced with RidgePlot
 - FindClusters is now more robust in making temp files
 - MetaDE support for combining p-values in DE testing

# Seurat 2.1.0 (2017-10-12)
## Added
- Support for using MAST and DESeq2 packages for differential expression testing in FindMarkers
- Support for multi-modal single-cell data via @assay slot

## Changed
- Default DE test changed to Wilcoxon rank sum test

# Seurat 2.0.1 (2017-08-18)
## Added
 - Now available on CRAN
 - Updated documentation complete with examples
 - Example datasets: `pbmc_small` and `cc.genes`
 - C++ implementation for parts of FindVariableGenes
 - Minor bug fixes

# Seurat 2.0.0 (2017-07-26)
## Added
- New method for aligning scRNA-seq datasets
- Significant code restructuring
- New methods for scoring gene expression and cell-cycle phases
- New visualization features (do.hover, do.identify)
