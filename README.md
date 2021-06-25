[![Build Status](https://travis-ci.com/satijalab/seurat.svg?branch=master)](https://travis-ci.com/github/satijalab/seurat)
[![AppVeyor build status](https://ci.appveyor.com/api/projects/status/github/satijalab/seurat?branch=master&svg=true)](https://ci.appveyor.com/project/satijalab/seurat)
[![CRAN Version](https://www.r-pkg.org/badges/version/Seurat)](https://cran.r-project.org/package=Seurat)
[![CRAN Downloads](https://cranlogs.r-pkg.org/badges/Seurat)](https://cran.r-project.org/package=Seurat)

# Seurat v4.0 + experimental optimization for WilcoxDETest

Seurat is an R toolkit for single cell genomics, developed and maintained by the Satija Lab at NYGC.

Instructions, documentation, and tutorials can be found at:

* https://satijalab.org/seurat

Seurat is also hosted on GitHub, you can view and clone the repository at

* https://github.com/satijalab/seurat

Seurat has been successfully installed on Mac OS X, Linux, and Windows, using the devtools package to install directly from GitHub

Improvements and new features will be added on a regular basis, please post on the [github page](https://github.com/satijalab/seurat) with any questions or if you would like to contribute

For a version history/changelog, please see the [NEWS file](https://github.com/satijalab/seurat/blob/master/NEWS.md).

## Original version
```
> install.packages("Seurat")
> library(Seurat)
> library(SeuratData)
> library(future)
> pbmc <- LoadData("pbmc3k", type = "pbmc3k.final")

> plan(multiprocess(workers = 1), .cleanup=TRUE)
> system.time(a0 <- FindAllMarkers(pbmc))
   user  system elapsed
 66.747   0.329  67.102
> save(a0, file="pbmc3k_all_org.RData")

> plan(multiprocess(workers = 2), .cleanup=TRUE)
> system.time(a2 <- FindAllMarkers(pbmc))
   user  system elapsed
 71.175   5.599  39.678
> all.equal(a0, a2)
[1] TRUE

> plan(multiprocess(workers = 4), .cleanup=TRUE)
> system.time(a4 <- FindAllMarkers(pbmc))
   user  system elapsed
 78.607   7.777  24.716
> all.equal(a0, a4)
[1] TRUE
>
> plan(multiprocess(workers = 8), .cleanup=TRUE)
> system.time(a8 <- FindAllMarkers(pbmc))
   user  system elapsed
 96.006  13.866  18.974
> all.equal(a0, a8)
[1] TRUE
```

## Modified version
```
> devtools::install_github("mogushi/Seurat")
> library(Seurat)
> library(SeuratData)
> library(future)
> pbmc <- LoadData("pbmc3k", type = "pbmc3k.final")

load(file="pbmc3k_all_org.RData")

> plan(multiprocess(workers = 1), .cleanup=TRUE)
> system.time(a1 <- FindAllMarkers(pbmc))
   user  system elapsed
 35.304   0.271  35.586   # original: 67.102
 > all.equal(a0, a1)
[1] TRUE

> plan(multiprocess(workers = 2), .cleanup=TRUE)
> system.time(a2 <- FindAllMarkers(pbmc))
   user  system elapsed
 39.802   5.029  25.205   # original: 39.678
> all.equal(a0, a2)
[1] TRUE

> plan(multiprocess(workers = 4), .cleanup=TRUE)
> system.time(a4 <- FindAllMarkers(pbmc))
   user  system elapsed
 42.361   6.146  16.203   # original: 24.716
> all.equal(a0, a4)
[1] TRUE

> plan(multiprocess(workers = 8), .cleanup=TRUE)
> system.time(a8 <- FindAllMarkers(pbmc))
   user  system elapsed
 50.166   9.109  12.891   # original: 18.974
> all.equal(a0, a8)
[1] TRUE
```

## sessionInfo()
```
> sessionInfo()
R version 4.1.0 (2021-05-18)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Ubuntu 20.04.2 LTS

Matrix products: default
BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.9.0
LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.9.0

locale:
 [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8
 [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8
 [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C
[10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base

other attached packages:
[1] future_1.21.0             pbmc3k.SeuratData_3.1.4
[3] panc8.SeuratData_3.0.2    ifnb.SeuratData_3.1.0
[5] hcabm40k.SeuratData_3.0.0 cbmc.SeuratData_3.1.4
[7] SeuratData_0.2.1          SeuratObject_4.0.2
[9] Seurat_4.0.3

loaded via a namespace (and not attached):
  [1] nlme_3.1-152          matrixStats_0.59.0    spatstat.sparse_2.0-0
  [4] RcppAnnoy_0.0.18      RColorBrewer_1.1-2    httr_1.4.2
  [7] sctransform_0.3.2     tools_4.1.0           utf8_1.2.1
 [10] R6_2.5.0              irlba_2.3.3           rpart_4.1-15
 [13] KernSmooth_2.23-20    uwot_0.1.10           mgcv_1.8-36
 [16] lazyeval_0.2.2        colorspace_2.0-1      tidyselect_1.1.1
 [19] gridExtra_2.3         compiler_4.1.0        cli_2.5.0
 [22] plotly_4.9.4          scales_1.1.1          lmtest_0.9-38
 [25] spatstat.data_2.1-0   ggridges_0.5.3        pbapply_1.4-3
 [28] rappdirs_0.3.3        goftest_1.2-2         stringr_1.4.0
 [31] digest_0.6.27         spatstat.utils_2.2-0  pkgconfig_2.0.3
 [34] htmltools_0.5.1.1     parallelly_1.26.0     fastmap_1.1.0
 [37] htmlwidgets_1.5.3     rlang_0.4.11          rstudioapi_0.13
 [40] shiny_1.6.0           generics_0.1.0        zoo_1.8-9
 [43] jsonlite_1.7.2        ica_1.0-2             dplyr_1.0.6
 [46] magrittr_2.0.1        patchwork_1.1.1       Matrix_1.3-4
 [49] Rcpp_1.0.6            munsell_0.5.0         fansi_0.5.0
 [52] abind_1.4-5           reticulate_1.20       lifecycle_1.0.0
 [55] stringi_1.6.2         MASS_7.3-54           Rtsne_0.15
 [58] plyr_1.8.6            grid_4.1.0            parallel_4.1.0
 [61] listenv_0.8.0         promises_1.2.0.1      ggrepel_0.9.1
 [64] crayon_1.4.1          miniUI_0.1.1.1        deldir_0.2-10
 [67] lattice_0.20-44       cowplot_1.1.1         splines_4.1.0
 [70] tensor_1.5            ps_1.6.0              pillar_1.6.1
 [73] igraph_1.2.6          spatstat.geom_2.2-0   future.apply_1.7.0
 [76] reshape2_1.4.4        codetools_0.2-18      leiden_0.3.8
 [79] glue_1.4.2            data.table_1.14.0     png_0.1-7
 [82] vctrs_0.3.8           httpuv_1.6.1          gtable_0.3.0
 [85] RANN_2.6.1            purrr_0.3.4           spatstat.core_2.1-2
 [88] polyclip_1.10-0       tidyr_1.1.3           scattermore_0.7
 [91] ggplot2_3.3.4         mime_0.10             xtable_1.8-4
 [94] later_1.2.0           survival_3.2-11       viridisLite_0.4.0
 [97] tibble_3.1.2          cluster_2.1.2         globals_0.14.0
[100] fitdistrplus_1.1-5    ellipsis_0.3.2        ROCR_1.0-11
```
