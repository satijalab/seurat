# Seurat v5.4.0

## Test environments
* Ubuntu 20.04 (local) (R 4.3.2)
* macOS 15.6.1 (local) (R 4.5.1)
* Ubuntu 24.04 (GitHub Actions Runner): R-oldrelease, R-release
* [win-builder](https://win-builder.r-project.org/): R-oldrelease, R-release, R-devel

## R CMD check results
**Status: OK**

```
* checking CRAN incoming feasibility ... [42s] NOTE
Maintainer: 'Rahul Satija <seurat@nygenome.org>'

Suggests or Enhances not in mainstream repositories:
  BPCells, presto
Availability using Additional_repositories specification:
  BPCells   yes   https://bnprks.r-universe.dev   
  presto    yes   https://satijalab.r-universe.dev
```

The maintainer remains Rahul Satija and the email is correct. BPCells and presto are hosted on R-universe and used conditionally in Seurat.

```
* checking package dependencies ... NOTE
Packages suggested but not available for checking: 'BPCells', 'presto'
```

BPCells and presto are hosted on R-universe and used conditionally in Seurat.

```
* checking DESCRIPTION meta-information ... NOTE
Author field differs from that derived from Authors@R
```

There seems to be a slight difference in ORCID formatting between Author and Authors@R; the information in both is the same.

```
* checking Rd cross-references ... NOTE
Package unavailable to check Rd xrefs: 'BPCells'
```

BPCells is listed under "Suggests"; it is hosted on R-universe and used conditionally in Seurat.

## Reverse dependency check results
We checked 80 reverse dependencies (40 from CRAN + 40 from Bioconductor), comparing R CMD check results across CRAN and dev versions of this package, and saw no new problems.
