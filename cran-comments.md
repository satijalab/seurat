# Seurat v3.2.0

## Test environments
* local Ubuntu 16.04.6 install, R 3.6.1
* local Ubuntu 18.04.4 install, R 4.0.1
* local Windows 10 install, R 4.0.0
* Ubuntu 16.04.6 (on travis-ci), R 4.0.0, R devel
* macOS 10.13.6 (on travis-ci), R 4.0.2
* Windows Server 2012 R2 (on AppVeyor), R 4.0.2 Patched
* win-builder (oldrelease, release, devel)

## R CMD check results
There were no ERRORs or WARNINGs

There were 3 NOTEs:

* checking CRAN incoming feasibility ... NOTE
  Maintainer: ‘Paul Hoffman <nygcSatijalab@nygenome.org>’

  Suggests or Enhances not in mainstream repositories:
    loomR
  Availability using Additional_repositories specification:
    loomR      yes   https://mojaveazure.github.io/loomR

  The package we suggest, loomR, is currently under development and not yet available on CRAN. This package is not required for core functionality of Seurat.

* checking package dependencies ... NOTE
  Package suggested but not available for checking: 'loomR'

  loomR is a suggested package hosted on a custom repository and maintained by us (both the package and repository).

* checking Rd cross-references ... NOTE
  Package unavailable to check Rd xrefs: 'loomR'

  This is a suggested package that we maintain. We have checked the cross references to ensure they link to the correct help pages.

## Downstream dependencies

There are three pacakges that imports Seurat: scMappR, Signac, and SoupX; this update does not impact their functionality

There are four packages that suggest Seurat: BisqueRNA, clustree, nanny, Rmagic, singleCellHaystack, treefit; this update does not impact their functionality.
