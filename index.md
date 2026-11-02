<img src="images/seurat_banner.jpg" alt="Seurat" class="seurat-banner" width="1100">

## **Seurat v6**

We are excited to release Seurat v6! To install, please follow the instructions on our [install](articles/install.html) page. This update brings the following new features and functionality:

* **Sample-level analysis with scSLIDE:** Single-cell studies routinely profile large numbers of samples across donors, conditions, and time points, highlighting the need---especially from a translational perspective---to investigate heterogeneity at the sample level rather than focusing solely on cell-to-cell variability. In Seurat v6, we introduce scSLIDE (Single-cell Sample-Level Inference using Density Estimation) (Jiang, et al., _Cell_ 2026), a framework for defining sample-level representations that captures both cell type and cell state variation to enable informative and integrative analyses across samples.

    - Paper: [Reconstructing developmental and disease progression with sample-level embeddings](https://www.cell.com/)
    - Vignette: **coming soon**\
\
* **Fast, scalable analysis:** Seurat v6 builds on previously introduced sketch-based analysis and BPCells-backed on-disk storage to make large datasets easier to analyze and explore. Performance improvements throughout the core pipeline accelerate normalization, variable feature selection, scaling, SCTransform, PCA, neighbor finding, clustering, UMAP, differential expression, and RPCA/CCA integration. These improvements work with existing workflows, with optional multithreading through `setThreads()` providing additional speedups for supported functions.

    - Release notes: [Performance improvements and multithreading](articles/announcements.html#october-2026)\
\
* **Expanded spatial support for 10x Genomics' Visium HD and Atera platforms:** Seurat supports the analysis of [Visium HD](https://www.10xgenomics.com/platforms/visium/product-family) data at both binned and segmented (~single-cell) resolutions. We also add support for loading [Atera](https://www.10xgenomics.com/platforms/atera) data, including expression matrices, cell boundaries, transcript coordinates, and morphology images. Together, these tools enable the exploration of tissue architecture across sequencing- and imaging-based platforms.

    - Vignette: [Analysis of Visium HD binned data](articles/visiumhd_analysis_vignette.html)
    - Vignette: [Analysis of Visium HD with cell segmentations](articles/visiumhd_analysis_cell_segmentations.html)\
\
* **Cell type annotation with Pan-human Azimuth:** Pan-human Azimuth (Sarkar<sup>\*</sup>, Li<sup>\*</sup>, et al., _bioRxiv_ 2026) is a neural network classifier that annotates human single-cell and single-nucleus RNA-seq data, across tissues and technologies, into a consistent and interpretable hierarchical cell typology. The companion `AzimuthAPI` R package allows users to bring these annotations into Seurat workflows by running the model on a Seurat object locally or via the cloud, returning cell type labels, confidence scores, and latent embeddings for downstream analysis.

    - Preprint: [Organism-scale annotation with Pan-human Azimuth](https://www.biorxiv.org/content/10.64898/2026.07.16.738997v1)
    - Website: [Pan-human Azimuth](https://satijalab.org/pan_human_azimuth/)
    - Vignette: [Annotating Visium HD cell segmentations](articles/visiumhd_analysis_cell_segmentations.html#cell-type-annotation-with-pan-human-azimuth)\
\
* **Interpreting cellular responses with RNA Fingerprinting:** RNA Fingerprinting (Grabski, et al., _Cell_ 2026) is a statistical framework that maps transcriptional responses in scRNA-seq experiments to reference perturbation dictionaries. The companion R package probabilistically assigns query cells or groups of cells within a Seurat object to candidate perturbations while accounting for uncertainty, allowing users to answer important questions that link observed cellular states to experimentally characterized responses.

    - Paper: [Mapping transcriptional responses to cellular perturbations with RNA fingerprinting](https://doi.org/10.1016/j.cell.2026.09.017)
    - GitHub: [RNA Fingerprinting](https://github.com/satijalab/rna-fingerprinting)\
    - Vignette: **coming soon**\

## **Older announcements**

<details>

<summary>Seurat v5</summary>

<h2 id="seurat-v5"> <img src="images/SeuratV5.png" alt="Seurat v5" width=50 height=50> <strong> Seurat v5</strong></h2>

We are excited to release Seurat v5! To install, please follow the instructions on our [install](articles/install.html) page. This update brings the following new features and functionality:

* **Integrative multimodal analysis:** The cellular transcriptome is just one aspect of cellular identity, and recent technologies enable routine profiling of chromatin accessibility, histone modifications, and protein levels from single cells. In Seurat v5, we introduce 'bridge integration', a statistical method to integrate experiments measuring different modalities (i.e. separate scRNA-seq and scATAC-seq datasets), using a separate multiomic dataset as a molecular 'bridge'. For example, we demonstrate how to map scATAC-seq datasets onto scRNA-seq datasets, to assist users in interpreting and annotating data from new modalities.\
\
We recognize that while the goal of matching shared cell types across datasets may be important for many problems, users may also be concerned about which method to use, or that integration could result in a loss of biological resolution. In Seurat v5, we also introduce flexible and streamlined workflows for the integration of multiple scRNA-seq datasets. This makes it easier to explore the results of different integration methods, and to compare these results to a workflow that excludes integration steps.

    - Paper: [Dictionary learning for integrative, multimodal, and scalable single-cell analysis](https://doi.org/10.1038/s41587-023-01767-y)
    - Vignette: [Streamlined integration of scRNA-seq data](articles/integration.html)
    - Vignette: [Cross-modality bridge integration](articles/integration_bridge.html)
    - Website: [Azimuth-ATAC, reference-mapping for scATAC-seq datasets](https://azimuth.hubmapconsortium.org/references/)\
\
* **Flexible, interactive, and highly scalable analysis:** The size and scale of single-cell sequencing datasets is rapidly increasing, outpacing even Moore's law. In Seurat v5, we introduce new infrastructure and methods to analyze, interpret, and explore exciting datasets spanning millions of cells, even if they cannot be fully loaded into memory. We introduce support for 'sketch'-based analysis, where representative subsamples of a large dataset are stored in-memory to enable rapid and iterative analysis - while the full dataset remains accessible via on-disk storage.\
\
We enable high-performance via the BPCells package, developed by Ben Parks in the Greenleaf Lab. The BPCells package enables high-performance analysis via innovative bit-packing compression techniques, optimized C++ code, and use of streamlined and lazy operations.

    - Vignette: [Sketch-based clustering of 1.3M brain cells (10x Genomics)](articles/sketch_analysis.html)
    - Vignette: [Sketch-based integration of 1M healthy and diabetic PBMC (Parse Biosciences)](articles/ParseBio_sketch_integration.html)
    - Vignette: [Mapping 1.5M cells from multiple studies to an Azimuth reference](articles/COVID_SCTMapping.html)
    - Vignette: [Interacting with BPCell matrices in Seurat v5](articles/bpcells_interaction_vignette.html)
    - BPCells R Package: [Scaling Single Cell Analysis to Millions of Cells](https://bnprks.github.io/BPCells/)\
\
* **Analysis of sequencing and imaging-based spatial datasets:** Spatially resolved datasets are redefining our understanding of cellular interactions and the organization of human tissues. Both sequencing-based(i.e. Visium, SLIDE-seq, etc.), and imaging-based (MERFISH/Vizgen, Xenium, CosMX, etc.) technologies have unique advantages, and require tailored analytical methods and software infrastructure. In Seurat v5, we introduce flexible and diverse support for a wide variety of spatially resolved data types, and support for analytical techniqiues for scRNA-seq integration, deconvolution, and niche identification.

    - Vignette: [Analysis of spatial datasets (Sequencing-based)](articles/spatial_vignette.html)
    - Vignette: [Analysis of spatial datasets (Imaging-based)](articles/spatial_vignette_2.html)\
\
</details>

## **About Seurat**

Seurat is an R package designed for QC, analysis, and exploration of single-cell RNA-seq data. Seurat aims to enable users to identify and interpret sources of heterogeneity from single-cell transcriptomic measurements, and to integrate diverse types of single-cell data.

All methods emphasize clear, attractive, and interpretable visualizations, and were designed to be [easily used](articles/get_started.html) by both dry-lab and wet-lab researchers.

If you use Seurat in your research, please consider citing:

* [Seurat v6] [Jiang, et al., Cell 2026](https://doi.org/10.1016/j.cell.2026.10.018)
* [Seurat v5] [Hao, et al., Nature Biotechnology 2023](https://www.nature.com/articles/s41587-023-01767-y)
* [Seurat v4] [Hao<sup>\*</sup>, Hao<sup>\*</sup>, et al., Cell 2021](https://doi.org/10.1016/j.cell.2021.04.048)
* [Seurat v3] [Stuart<sup>\*</sup>, Butler<sup>\*</sup>, et al., Cell 2019](https://www.cell.com/cell/fulltext/S0092-8674(19)30559-8)
* [Seurat v2] [Butler, et al., Nature Biotechnology 2018](https://doi.org/10.1038/nbt.4096)
* [Seurat v1] [Satija<sup>\*</sup>, Farrell<sup>\*</sup>, et al., Nature Biotechnology 2015](https://doi.org/10.1038/nbt.3192)

Seurat is developed and maintained by the [Satija Lab](authors.html) and is released under the MIT license.
