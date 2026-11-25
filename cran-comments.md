# Seurat v6.0.1

## Resubmission

This resubmission fixes the `-Wdeprecated-declarations` warning in `src/pls.cpp` reported on `r-devel-windows-x86_64` during incoming checks. The BDCSVD computation options are now passed as a template parameter when building against Eigen >= 5.0.0, as recommended.

## Test environments
* Ubuntu 24.04 (local) (R 4.6.0)
* macOS 15.6.1 (local) (R 4.6.1)
* Ubuntu 26.04 (GitHub Actions Runner): R-oldrelease, R-release, R-devel
* [win-builder](https://win-builder.r-project.org/): R-oldrelease, R-release, R-devel

## R CMD check results

**Status: 1 NOTE**

* checking CRAN incoming feasibility ... [33s] NOTE
Maintainer: 'Rahul Satija <seurat@nygenome.org>'

Suggests or Enhances not in mainstream repositories:
  BPCells
Availability using Additional_repositories specification:
  BPCells   yes   https://bnprks.r-universe.dev   
  ?           ?   https://satijalab.r-universe.dev

Found the following (possibly) invalid URLs:
  URL: https://stackoverflow.com/questions/3942878/how-to-decide-font-color-in-white-or-black-depending-on-background-color
    From: man/BGTextColor.Rd
          man/contrast-theory.Rd
    Status: 403
    Message: Forbidden

- The maintainer remains Rahul Satija and the email is correct.
- BPCells and presto are hosted on R-universe and used conditionally in Seurat.
- The URL can be accessed normally manually so does not seem to be invalid.

## Reverse dependency check results

No changes to reverse dependencies are expected. This release involves a single bug patch to new functionality introduced in the previous version.