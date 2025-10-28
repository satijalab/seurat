# Seurat v5.3.1

## Test environments
* ubuntu 20.04 (local) (R 4.3.2)
* macOS 15.6.1 (local) (R 4.5.1)
* ubuntu 24.04 (GitHub Actions Runner): R-oldrelease, R-release, R-devel
* [win-builder](https://win-builder.r-project.org/): R-oldrelease, R-release, R-devel
* [macOS builder](https://mac.r-project.org/macbuilder/submit.html): R-release, R-devel

## R CMD check results
**Status: OK**

There were no ERRORs or WARNINGs.

There were four NOTEs given by win-builder (oldrel):

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
We checked 145 reverse dependencies (69 from CRAN + 76 from Bioconductor), comparing R CMD check results across CRAN and dev versions of this package.

 * We saw 0 new problems
 * We failed to check 15 packages (conos, CytoSimplex, DIscBIO, DR.SC, harmony, Platypus, PopComm, PoweREST, PRECAST, ProFAST, RepeatedHighDim, RESET, rPanglaoDB, sccore, SpaCCI)
