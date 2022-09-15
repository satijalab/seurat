![](articles/assets/seurat_banner.jpg)

# Official release of Seurat 4.0

We are excited to release Seurat v4.0! This update brings the following new features and functionality:

* **Integrative multimodal analysis.** The ability to make simultaneous measurements of multiple data types from the same cell, known as multimodal analysis, represents a new and exciting frontier for single-cell genomics. In Seurat v4, we introduce weighted nearest neighbor (WNN) analysis, an unsupervised strategy to learn the information content of each modality in each cell, and to define cellular state based on a weighted combination of both modalities.
    In our new paper, we generate a CITE-seq dataset featuring paired measurements of the transcriptome and 228 surface proteins, and leverage WNN to define a multimodal reference of human PBMC. You can use WNN to analyze multimodal data from a variety of technologies, including CITE-seq, ASAP-seq, 10X Genomics ATAC + RNA, and SHARE-seq.

    - Paper: [Integrated analysis of multimodal single-cell data](https://doi.org/10.1016/j.cell.2021.04.048)
    - Vignette: [Multimodal clustering of a human bone marrow CITE-seq dataset](articles/weighted_nearest_neighbor_analysis.html)
    - Portal: [Click here](https://atlas.fredhutch.org/nygc/multimodal-pbmc/)
    - Dataset: [Download here](https://atlas.fredhutch.org/data/nygc/multimodal/pbmc_multimodal.h5seurat)

* **Rapid mapping of query datasets to references.** We introduce Azimuth, a workflow to leverage high-quality reference datasets to rapidly map new scRNA-seq datasets (queries). For example, you can map any scRNA-seq dataset of human PBMC onto our reference, automating the process of visualization, clustering annotation, and differential expression. Azimuth can be run within Seurat, or using a standalone web application that requires no installation or programming experience.

    - Vignette: [Mapping scRNA-seq queries onto reference datasets](articles/multimodal_reference_mapping.html)
    - Web app: [Automated mapping, visualization, and annotation of scRNA-seq datasets from human PBMC](https://azimuth.hubmapconsortium.org/)

Additional speed and usability updates: We have made minor changes in v4, primarily to improve the performance of Seurat v4 on large datasets. These changes substantially improve the speed and memory requirements, but do not adversely impact downstream results. We provide a detailed description of key changes [here](articles/v4_changes.html). Users who wish to fully reproduce existing results can continue to do so by continuing to install Seurat v3.  

We believe that users who are familiar with Seurat v3 should experience a smooth transition to Seurat v4. While we have introduced extensive new functionality, existing workflows, functions, and syntax are largely unchanged in this update. In addition, Seurat objects that have been previously generated in Seurat v3 can be seamlessly loaded into Seurat v4 for further analysis.

# About Seurat

Seurat is an R package designed for QC, analysis, and exploration of single-cell RNA-seq data. Seurat aims to enable users to identify and interpret sources of heterogeneity from single-cell transcriptomic measurements, and to integrate diverse types of single-cell data.

If you use Seurat in your research, please considering citing:

* [Hao\*, Hao\*, et al., Cell 2021](https://doi.org/10.1016/j.cell.2021.04.048) [Seurat V4]
* [Stuart\*, Butler\*, et al., Cell 2019](https://www.cell.com/cell/fulltext/S0092-8674(19)30559-8) [Seurat V3]
* [Butler\* et al., Nat Biotechnol 2018](https://doi.org/10.1038/nbt.4096) [Seurat V2]
* [Satija\*, Farrell\*, et al., Nat Biotechnol 2015](https://doi.org/10.1038/nbt.3192) [Seurat V1]


All methods emphasize clear, attractive, and interpretable visualizations, and were designed to be [easily used](articles/get_started.html) by both dry-lab and wet-lab researchers.

Seurat is developed and maintained by the [Satija lab](authors.html) and is released under the MIT license.
