# Seurat v5.2.1
This update is a small patch release to address failing M1mac checks ahead of the deadline on 2025-01-28.

## Test environments
* macOS 14.6.1 (Apple M3, arm64, local install): R-release
* Ubuntu 22.04 (Docker container): R-release
* [win-builder](https://win-builder.r-project.org/): R-oldrelease, R-release, R-devel
* [macOS builder](https://mac.r-project.org/macbuilder/submit.html): R-release, R-devel

## R CMD check results
There were no ERRORs or WARNINGs but there was one NOTE.

```
❯ checking CRAN incoming feasibility ... [15s/92s] NOTE
  Maintainer: ‘Rahul Satija <seurat@nygenome.org>’
  
  Suggests or Enhances not in mainstream repositories:
    BPCells, enrichR, presto
  Availability using Additional_repositories specification:
    BPCells   yes   https://bnprks.r-universe.dev   
    enrichR   yes   https://satijalab.r-universe.dev
    presto    yes   https://satijalab.r-universe.dev
```

`BPCells`, `enrichR`, and `presto` are hosted on R-universe and used conditionally in `Seurat`.

## Reverse dependency check results
We checked 65 reverse dependencies, comparing `R CMD check` results across the CRAN and dev versions of this package, and saw no new problems.
