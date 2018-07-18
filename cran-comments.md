# Seurat v2.3.4

## Test environments
* local Ubuntu 16.04 and 18.04 installs, R 3.4.4
* Ubuntu 12.04 (on travis-ci), R 3.5.0
* macOS 10.12.6 (on travis-ci), R 3.5.0
* Windows Server 2012 (on AppVeyor), R 3.5.0
* win-builder (devel, release)

## R CMD check results
There were no ERRORs or WARNINGs

There were 2 NOTEs:

* checking CRAN incoming feasibility ... NOTE
  Maintainer: 'Paul Hoffman <seuratpackage@gmail.com>'
  Suggests or Enhances not in mainstream repositories:
    loomR
  Availability using Additional_repositories specification:
    loomR   yes   https://mojaveazure.github.io/loomR

  The package we suggest, loomR, is currently underdevelopment and not yet available on CRAN. This package is not required for core functionality of Seurat.

* checking package dependencies ... NOTE
  Package suggested but not available for checking: 'loomR'

  This is a suggested package hosted on a custom repository and maintained by us (both the package and repository).

## Downstream dependencies

There are currently no downstream dependencies for this package.
