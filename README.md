[![Build Status](https://travis-ci.com/satijalab/seurat.svg)](https://travis-ci.com/satijalab/seurat)
[![AppVeyor build status](https://ci.appveyor.com/api/projects/status/github/satijalab/seurat?svg=true)](https://ci.appveyor.com/project/satijalab/seurat)
[![CRAN Version](https://www.r-pkg.org/badges/version/Seurat)](https://cran.r-project.org/package=Seurat)
[![CRAN Downloads](https://cranlogs.r-pkg.org/badges/Seurat)](https://cran.r-project.org/package=Seurat)

# Seurat v3.1.4


Seurat is an R toolkit for single cell genomics, developed and maintained by the Satija Lab at NYGC.

Instructions, documentation, and tutorials can be found at:
* https://satijalab.org/seurat

Seurat is also hosted on GitHub, you can view and clone the repository at
* https://github.com/satijalab/seurat

Seurat has been successfully installed on Mac OS X, Linux, and Windows, using the devtools package to install directly from GitHub

Improvements and new features will be added on a regular basis, please contact seuratpackage@gmail.com with any questions or if you would like to contribute

Version History

August 20, 2019
* Version 3.1
* Changes:
  * Support for SCTransform integration workflows
  * Integration speed ups: reference-based integration + reciprocal PCA

April 12, 2019
* Version 3.0
* Changes:
  * Preprint published describing new methods for identifying anchors across single-cell datasets
  * Restructured Seurat object with native support for multimodal data
  * Parallelization support via future

July 20, 2018
* Version 2.4
* Changes:
  * Java dependency removed and functionality rewritten in Rcpp 

March 22, 2018
* Version 2.3
* Changes:
  * New utility functions
  * Speed and efficiency improvments

January 10, 2018
* Version 2.2
* Changes:
   * Support for multiple-dataset alignment with RunMultiCCA and AlignSubspace
   * New methods for evaluating alignment performance

October 12, 2017
* Version 2.1
* Changes:
   * Support for using MAST and DESeq2 packages for differential expression testing in FindMarkers
   * Support for multi-modal single-cell data via \@assay slot

July 26, 2017
* Version 2.0
* Changes:
   * Preprint released for integrated analysis of scRNA-seq across conditions, technologies and species
   * Significant restructuring of code to support clarity and dataset exploration
   * Methods for scoring gene expression and cell-cycle phase

October 4, 2016
* Version 1.4 released
* Changes:
   * Improved tools for cluster evaluation/visualizations
   * Methods for combining and adding to datasets

August 22, 2016:
* Version 1.3 released
* Changes :
    * Improved clustering approach - see FAQ for details
    * All functions support sparse matrices
    * Methods for removing unwanted sources of variation
    * Consistent function names
    * Updated visualizations

May 21, 2015:
* Drop-Seq manuscript published. Version 1.2 released
* Changes :
  * Added support for spectral t-SNE and density clustering
  * New visualizations - including pcHeatmap, dot.plot, and feature.plot
  * Expanded package documentation, reduced import package burden
  *  Seurat code is now hosted on GitHub, enables easy install through devtools
  * Small bug fixes

April 13, 2015:
* Spatial mapping manuscript published. Version 1.1 released (initial release)
