![](articles/assets/seurat_banner.jpg)

## **Beta release of Seurat v5**

We are excited to release an initial beta version of Seurat v5! This update brings the following new features and functionality:

* **Analysis of sequencing and imaging-based spatial datasets:** Spatially resolved datasets are redefining our understanding of cellular interactions and the organization of human tissues. Both sequencing-based(i.e. Visium, SLIDE-seq, etc.), and imaging-based (MERFISH/Vizgen, Xenium, CosMX, etc.) technologies have unique advantages, and require tailored analytical methods and software infrastructure. In Seurat v5, we introduce flexible and diverse support for a wide variety of spatially resolved data types, and support for analytical techniqiues for scRNA-seq integration, deconvolution, and niche identification.

    - Vignette: [Analysis of spatial datasets (Sequencing-based)](articles/seurat5_spatial_vignette.html)
    - Vignette: [Analysis of spatial datasets (Imaging-based)](articles/seurat5_spatial_vignette_2.html)\
\
* **Integrative multimodal analysis:** The cellular transcriptome is just one aspect of cellular identity, and recent technologies enable routine profiling of chromatin accessibility, histone modifications, and protein levels from single cells. In Seurat v5, we introduce 'bridge integration', a statistical method to integrate experiments measuring different modalities (i.e. separate scRNA-seq and scATAC-seq datasets), using a separate multiomic dataset as a molecular 'bridge'. For example, we demonstrate how to map scATAC-seq datasets onto scRNA-seq datasets, to assist users in interpreting and annotating data from new modalities.\
\
We recognize that while the goal of matching shared cell types across datasets may be important for many problems, users may also be concerned about which method to use, or that integration could result in a loss of biological resolution. In Seurat v5, we also introduce flexible and streamlined workflows for the integration of multiple scRNA-seq datasets. This makes it easier to explore the results of different integration methods, and to compare these results to a workflow that excludes integration steps.

    - Paper: [Dictionary learning for integrative, multimodal, and scalable single-cell analysis](https://doi.org/10.1101/2022.02.24.481684)
    - Vignette: [Streamlined integration of scRNA-seq data](articles/seurat5_integration.html)
    - Vignette: [Cross-modality bridge integration](articles/seurat5_integration_bridge.html)
    - Website: [Azimuth-ATAC, reference-mapping for scATAC-seq datasets](https://azimuth.hubmapconsortium.org/references/)\
\
* **Flexible, interactive, and highly scalable analsyis:** The size and scale of single-cell sequencing datasets is rapidly increasing, outpacing even Moore's law. In Seurat v5, we introduce new infrastructure and methods to analyze, interpret, and explore exciting datasets spanning millions of cells, even if they cannot be fully loaded into memory. We introduce support for 'sketch'-based analysis, where representative subsamples of a large dataset are stored in-memory to enable rapid and iterative analysis - while the full dataset remains accessible via on-disk storage.\
\
We enable high-performance via the BPCells package, developed by Ben Parks in the Greenleaf Lab. The BPCells package enables high-performance analysis via innovative bit-packing compression techniques, optimized C++ code, and use of streamlined and lazy operations.

    - Vignette: [Sketch-based clustering of 1.3M brain cells (10x Genomics)](articles/seurat5_sketch_analysis.html)
    - Vignette: [Sketch-based integration of 1M healthy and diabetic PBMC (Parse Biosciences)](articles/ParseBio_sketch_integration.html)
    - Vignette: [Mapping 1.5M cells from multiple studies to an Azimuth reference](articles/COVID_SCTMapping.html)
    - Vignette: [Interacting with BPCell matrices in Seurat v5](articles/seurat5_bpcells_interaction_vignette.html)
    - BPCells R Package: [Scaling Single Cell Analysis to Millions of Cells](https://bnprks.github.io/BPCells/)\
\
* **Backwards compatibility:** While Seurat v5 introduces new functionality, we have ensured that the software is backwards-compatible with previous versions, so that users will continue to be able to re-run existing workflows. As v5 is still in beta, the CRAN installation (`install.packages("Seurat")`) will continue to install Seurat v4, but users can opt-in to test Seurat v5 by following the instructions in our [install page](install.html).\

## **About Seurat**

Seurat is an R package designed for QC, analysis, and exploration of single-cell RNA-seq data. Seurat aims to enable users to identify and interpret sources of heterogeneity from single-cell transcriptomic measurements, and to integrate diverse types of single-cell data.

If you use Seurat in your research, please considering citing:

* [Hao, et al., bioRxiv 2022](https://doi.org/10.1101/2022.02.24.481684) [Seurat v5]
* [Hao\*, Hao\*, et al., Cell 2021](https://doi.org/10.1016/j.cell.2021.04.048) [Seurat v4]
* [Stuart\*, Butler\*, et al., Cell 2019](https://www.cell.com/cell/fulltext/S0092-8674(19)30559-8) [Seurat v3]
* [Butler, et al., Nat Biotechnol 2018](https://doi.org/10.1038/nbt.4096) [Seurat v2]
* [Satija\*, Farrell\*, et al., Nat Biotechnol 2015](https://doi.org/10.1038/nbt.3192) [Seurat v1]


All methods emphasize clear, attractive, and interpretable visualizations, and were designed to be [easily used](articles/get_started.html) by both dry-lab and wet-lab researchers.

Seurat is developed and maintained by the [Satija lab](authors.html) and is released under the MIT license.
