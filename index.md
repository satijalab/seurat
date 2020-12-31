# Beta release of Seurat 4.0

We are excited to release a beta version of Seurat v4.0! This update brings the following new features and functionality:

* **Integrative multimodal analysis.** The ability to make simultaneous measurements of multiple data types from the same cell, known as multimodal analysis, represents a new and exciting frontier for single-cell genomics. In Seurat v4, we introduce weighted nearest neighbor (WNN) analysis, an unsupervised strategy to learn the information content of each modality in each cell, and to define cellular state based on a weighted combination of both modalities.
    In our new preprint, we generate a CITE-seq dataset featuring paired measurements of the transcriptome and 228 surface proteins, and leverage WNN to define a multimodal reference of human PBMC. You can use WNN to analyze multimodal data from a variety of technologies, including CITE-seq, ASAP-seq, 10X Genomics ATAC + RNA, and SHARE-seq.

    - Preprint: [Integrated analysis of multimodal single-cell data](https://satijalab.org/v4preprint)
    - Vignette: [Multimodal clustering of a human bone marrow CITE-seq dataset](articles/weighted_nearest_neighbor_analysis.html)
    - Portal: [Click here](https://atlas.fredhutch.org/nygc/multimodal-pbmc/)
    - Dataset: [Download here](https://atlas.fredhutch.org/data/nygc/multimodal/pbmc_multimodal.h5seurat)

* **Rapid mapping of query datasets to references.** We introduce Azimuth, a workflow to leverage high-quality reference datasets to rapidly map new scRNA-seq datasets (queries). For example, you can map any scRNA-seq dataset of human PBMC onto our reference, automating the process of visualization, clustering annotation, and differential expression. Azimuth can be run within Seurat, or using a standalone web application that requires no installation or programming experience.

    - Vignette: [Mapping scRNA-seq queries onto reference datasets](articles/reference_mapping.html)
    - Web app: [Automated mapping, visualization, and annotation of scRNA-seq datasets from human PBMC]("../azimuth/")

Additional speed and usability updates: We have made minor changes in v4, primarily to improve the performance of Seurat v4 on large datasets. These changes substantially improve the speed and memory requirements, but do not adversely impact downstream results. We provide a detailed description of key changes [here](articles/v4_changes.html). Users who wish to fully reproduce existing results can continue to do so by continuing to install Seurat v3.  

We believe that users who are familiar with Seurat v3 should experience a smooth transition to Seurat v4. While we have introduced extensive new functionality, existing workflows, functions, and syntax are largely unchanged in this update. In addition, Seurat objects that have been previously generated in Seurat v3 can be seamlessly loaded into Seurat v4 for further analysis.

# Official release of Seurat 3.0

On April 16, 2019 - we officially updated the Seurat CRAN repository to release 3.0! 

We have been working on this update for the past year, and are excited to introduce new features and functionality, in particular:

* **Improved and expanded methods for single-cell integration.**  As described in [Stuart\*, Butler\*, et al., Cell 2019](https://www.cell.com/cell/fulltext/S0092-8674(19)30559-8) (bioRxiv preprint link [here](https://www.biorxiv.org/content/10.1101/460147v1)), Seurat v3 implements new methods to identify 'anchors' across diverse single-cell data types, in order to construct harmonized references, or to transfer information across experiments.
    - Vignette: [Stimulated vs. Control PBMCs](articles/immune_alignment.html)
    - Vignette: [Multiple Dataset Integration and Label Transfer](articles/integration.html)

* **Improved methods for normalization.** Seurat v3 includes support for sctransform, a new modeling approach for the normalization of single-cell data, described in a second [preprint](https://www.biorxiv.org/content/10.1101/576827v2). Compared to standard log-normalization, sctransform effectively removes technically-driven variation while preserving biological heterogeneity.
    - Vignette: [SCTransform vignette](articles/sctransform_vignette.html)

* **An efficiently restructured Seurat object, with an emphasis on multi-modal data.** We have carefully re-designed the structure of the Seurat object, with clearer documentation, and a flexible framework to easily switch between RNA, protein, cell hashing, batch-corrected / integrated, or imputed data. 
    - Vignette: [Multimodal vignette](articles/multimodal_vignette.html)
    - For a technical discussion of the object, please see the [developer's guide](https://github.com/satijalab/seurat/wiki)

While we are excited for users to upgrade, we are **committed to making this transition as smooth as possible, and to ensure that users can complete existing projects in Seurat v2** prior to upgrading:
* Users on all platforms can easily re-install Seurat v2, with detailed instructions [here](articles/install.html#previous)
* All website vignettes have been updated to v3, but v2 versions remain as well (look for the red button on the bottom-right of the screen).
* Seurat v3 includes an `UpdateSeuratObject()` function, so old objects can be analyzed with the upgraded version.
* We include a [command ‘cheat sheet’](articles/essential_commands.html), a brief introduction to new commands, data accessors, visualization, and multiple assays in Seurat v3.0 
* The [command ‘cheat sheet’](articles/essential_commands.html) also contains a translation guide between Seurat v2 and v3
<br><br>


# About Seurat

Seurat is an R package designed for QC, analysis, and exploration of single-cell RNA-seq data. Seurat aims to enable users to identify and interpret sources of heterogeneity from single-cell transcriptomic measurements, and to integrate diverse types of single-cell data.

If you use Seurat in your research, please considering citing:

* [Butler et al., Nature Biotechnology 2018](https://www.nature.com/articles/nbt.4096). 
* [Stuart\*, Butler\*, et al., Cell 2019](https://www.cell.com/cell/fulltext/S0092-8674(19)30559-8)

All methods emphasize clear, attractive, and interpretable visualizations, and were designed to be [easily used](articles/vignettes_overview.html) by both dry-lab and wet-lab researchers.

Seurat is developed and maintained by the Satija lab, in particular by [Andrew Butler](mailto:abutler@nygenome.org), [Paul Hoffman](mailto:phoffman@nygenome.org), [Tim Stuart](mailto:tstuart@nygenome.org), [Christoph Hafemeister](mailto:chafemeister@nygenome.org), and [Shiwei Zheng](mailto:szheng@nygenome.org), and is released under the GNU Public License (GPL 3.0). We are also grateful for significant ideas and code from [Jeff Farrell](mailto:jfarrell@g.harvard.edu), [Karthik Shekhar](mailto:karthik@broadinstitute.org), and [other generous contributors](../join_contact/).


## News

**October 13, 2020** Version 4.0 beta released

* Preprint published describing new methods for analysis of multimodal single-cell datasets

**July 16, 2020** Version 3.2 released

** Support for visualization and analysis of spatially resolved datasets

**August 20, 2019** Version 3.1 released

* Support for SCTransform integration workflows
* Integration speed ups: reference-based integration + reciprocal PCA

**April 16, 2019** Version 3.0 released

**November 2, 2018** Version 3.0 alpha released

* Preprint published describing new methods for identifying 'anchors' across single-cell datasets

**March 23, 2018** Version 2.3 released

* Improvements for speed and memory efficiency
* New vignette for analyzing ~250,000 cells from the Microwell-seq Mouse Cell Atlas dataset

**January 10, 2018:** Version 2.2 released

* Support for multiple-dataset alignment
* New methods for evaluating alignment performance

**October 16, 2017:** Version 2.1 released

* Support for multimodal single-cell data
* Support for MAST and DESeq2 packages for differential expression testing

**July 26, 2017:** Version 2.0 released

* Preprint published for integrated analysis of scRNA-seq datasets
* New methods for dataset integration, visualization, and exploration
* Significant restructuring of codebase to emphasize clarity and clear documentation

**October 4, 2016:** Version 1.4 released

* Added methods for negative binomial regression and differential expression testing for UMI count data
* New ways to merge and downsample Seurat objects

**August 22, 2016:** Version 1.3 released

* Improved clustering approach - see FAQ for details
* All functions support sparse matrices
* Methods for removing unwanted sources of variation
* Consistent function names
* Updated visualizations

**May 21, 2015:**
[Drop-Seq manuscript](http://www.cell.com/cell/abstract/S0092-8674(15)00549-8) published. Version 1.2 released

* Added support for spectral t-SNE (non-linear dimensional reduction), and density clustering
* New visualizations - including pcHeatmap, dot.plot, and feature.plot
* Expanded package documentation, reduced import package burden
* Seurat code is now hosted on GitHub, enables easy install through devtools package
* Small bug fixes

**April 13, 2015:**
[Spatial mapping manuscript](http://www.nature.com/nbt/journal/vaop/ncurrent/full/nbt.3192.html) published. Version 1.1 released

<script>
    // taken from https://github.com/twbs/bootstrap/issues/2415#issuecomment-4589184
    $(function () {
       var activeTab = $('[href=' + location.hash + ']');
       activeTab && activeTab.tab('show');
    });
</script>
