# Seurat v3.1.5

## Test environments
* local Ubuntu 16.04.6 and 18.04.2 installs, R 3.6.1
* local Windows 10 install, R 3.5.3, R-devel (4.1.0)
* Ubuntu 16.04.6 (on travis-ci), R 3.6.1
* macOS 10.13.3 (on travis-ci), R 3.6.1
* Windows Server 2012 R2 (on AppVeyor), R 3.6.1 Patched
* win-builder (oldrelease, release, devel, devel_gcc8)

## R CMD check results
There were no ERRORs or WARNINGs

There were 3 NOTEs:

* checking CRAN incoming feasibility ... NOTE
  Maintainer: ‘Paul Hoffman <nygcSatijalab@nygenome.org>’

  Suggests or Enhances not in mainstream repositories:
    loomR, SDMTools
  Availability using Additional_repositories specification:
    loomR      yes   https://mojaveazure.github.io/loomR
    SDMTools    no   ?

  One of the packages we suggest, loomR, is currently under development and not yet available on CRAN. This package is not required for core functionality of Seurat. The other package, SDMTools is also not required for core functionality. We are working to replace the function that calls SDMTools with a new function using non-orphaned packages

* checking package dependencies ... NOTE
  Package suggested but not available for checking: 'loomR'

  Suggests orphaned package: ‘SDMTools’

  loomR is a suggested package hosted on a custom repository and maintained by us (both the package and repository).

  SDMTools is not required for any essential functionality. We are working to replace the function that calls SDMTools with a new function using non-orphaned packages.

* checking Rd cross-references ... NOTE
  Package unavailable to check Rd xrefs: 'loomR'

  This is a suggested package that we maintain. We have checked the cross references to ensure they link to the correct help pages.

## Downstream dependencies

There are three pacakges that imports Seurat: multicross, scMappR, and Signac; this update does not impact their functionality

There are four packages that suggest Seurat: BisqueRNA, clustree, Rmagic, and treefit; this update does not impact their functionality.
