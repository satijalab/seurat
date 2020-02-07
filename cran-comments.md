# Seurat v3.1.3

## Test environments
* local Ubuntu 16.04.6 and 18.04.2 installs, R 3.6.1
* Ubuntu 16.04.6 (on travis-ci), R 3.6.1
* macOS 10.13.3 (on travis-ci), R 3.6.1
* Windows Server 2012 R2 (on AppVeyor), R 3.6.1 Patched
* win-builder (oldrelease, release, devel, devel_gcc8)

## R CMD check results
There were no ERRORs or WARNINGs

There were 3 NOTEs:

* checking CRAN incoming feasibility ... NOTE
    Maintainer: ‘Paul Hoffman <nygcSatijalab@nygenome.org>’

    New submission

    Package was archived on CRAN

    CRAN repository db overrides:
        X-CRAN-Comment: Archived on 2020-02-07 as check issues were not
            corrected in time.

    Suggests or Enhances not in mainstream repositories:
        loomR
    Availability using Additional_repositories specification:
        loomR   yes   https://mojaveazure.github.io/loomR

  This is a patch for an archived package Seurat. We were slow in fixing the errors and would like to get Seurat back on CRAN. In addition. the package we suggest, loomR, is currently underdevelopment and not yet available on CRAN. This package is not required for core functionality of Seurat.

* checking package dependencies ... NOTE
  Package suggested but not available for checking: 'loomR'

  This is a suggested package hosted on a custom repository and maintained by us (both the package and repository).

* checking Rd cross-references ... NOTE
  Package unavailable to check Rd xrefs: 'loomR'

  This is a suggested package that we maintain. We have checked the cross references to ensure they link to the correct help pages.

## Downstream dependencies

There is one pacakge that imports Seurat: multicross; this update does not impact its functionality

There are five packages that suggest Seurat: BisqueRNA, clustree, diem, iCellR, and Rmagic; this update does not impact their functionality.
