# Seurat v5.5.0

## Test environments
* Ubuntu 20.04 (local) (R 4.3.2)
* Ubuntu 24.04 (GitHub Actions Runner): R-oldrelease, R-release, R-devel
* macOS 15.6.1 (local) (R 4.5.1)
* [macos-builder](https://mac.r-project.org/macbuilder/submit.html): R-devel
* [win-builder](https://win-builder.r-project.org/): R-oldrelease, R-release, R-devel

## R CMD check results

**possible false positive(s) on R-devel - explained below**

* ERRORs on installation in R-devel are false positives and entirely unrelated to changes in this Seurat version; we have seen these stem from downstream dependencies (such as Rcpp) using header files that do not meet the new requirements posed in R-4.6. 
* As of 21 April, this seems to no longer be an issue, but adding here as a note.

**Status: OK**

**4 NOTEs**

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

We checked 81 reverse dependencies (47 from CRAN + 34 from Bioconductor), comparing R CMD check results across CRAN and dev versions of this package, and saw no new problems.
