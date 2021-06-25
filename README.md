[![Build Status](https://travis-ci.com/satijalab/seurat.svg?branch=master)](https://travis-ci.com/github/satijalab/seurat)
[![AppVeyor build status](https://ci.appveyor.com/api/projects/status/github/satijalab/seurat?branch=master&svg=true)](https://ci.appveyor.com/project/satijalab/seurat)
[![CRAN Version](https://www.r-pkg.org/badges/version/Seurat)](https://cran.r-project.org/package=Seurat)
[![CRAN Downloads](https://cranlogs.r-pkg.org/badges/Seurat)](https://cran.r-project.org/package=Seurat)

# Seurat v4.0 + experimental optimization for WilcoxDETTest

Seurat is an R toolkit for single cell genomics, developed and maintained by the Satija Lab at NYGC.

Instructions, documentation, and tutorials can be found at:

* https://satijalab.org/seurat

Seurat is also hosted on GitHub, you can view and clone the repository at

* https://github.com/satijalab/seurat

Seurat has been successfully installed on Mac OS X, Linux, and Windows, using the devtools package to install directly from GitHub

Improvements and new features will be added on a regular basis, please post on the [github page](https://github.com/satijalab/seurat) with any questions or if you would like to contribute

For a version history/changelog, please see the [NEWS file](https://github.com/satijalab/seurat/blob/master/NEWS.md).

```
> install.packages("Seurat")
> library(Seurat)
> library(SeuratData)
> library(future)
> pbmc <- LoadData("pbmc3k", type = "pbmc3k.final")

> load(file = "pbmc3k_all_org.RData")

> plan(multiprocess(workers = 1), .cleanup=TRUE)
> system.time(a1 <- FindAllMarkers(pbmc))
   user  system elapsed
 66.747   0.329  67.102
> a0 <- a1
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

```
> devtools::install_github("mogushi/Seurat")
> library(Seurat)
> library(SeuratData)
> library(future)
> pbmc <- LoadData("pbmc3k", type = "pbmc3k.final")

> load(file = "pbmc3k_all_org.RData")

> plan(multiprocess(workers = 1), .cleanup=TRUE)
> system.time(a1 <- FindAllMarkers(pbmc))
   user  system elapsed
 35.304   0.271  35.586
> all.equal(a0, a1)
[1] TRUE

> plan(multiprocess(workers = 2), .cleanup=TRUE)
> system.time(a2 <- FindAllMarkers(pbmc))
   user  system elapsed
 39.802   5.029  25.205
> all.equal(a0, a2)
[1] TRUE

> plan(multiprocess(workers = 4), .cleanup=TRUE)
> system.time(a4 <- FindAllMarkers(pbmc))
   user  system elapsed
 42.361   6.146  16.203
> all.equal(a0, a4)
[1] TRUE

> plan(multiprocess(workers = 8), .cleanup=TRUE)
> system.time(a8 <- FindAllMarkers(pbmc))
   user  system elapsed
 50.166   9.109  12.891
> all.equal(a0, a8)
[1] TRUE
```
