# Seurat v3.0.1

## Test environments
* local Ubuntu 16.04.6 and 18.04.2 installs, R 3.5.3
* Ubuntu 14.04.5 (on travis-ci), R 3.6.0
* macOS 10.13.3 (on travis-ci), R 3.6.0
* Windows Server 2012 (on AppVeyor), R 3.6.0 Patched
* win-builder (oldrelease, release, devel)

## R CMD check results
There were no ERRORs or WARNINGs

There were 2 NOTEs:

* checking CRAN incoming feasibility ... NOTE
    Maintainer: ‘Paul Hoffman <nygcSatijalab@nygenome.org>’
    Suggests or Enhances not in mainstream repositories:
        loomR
    Availability using Additional_repositories specification:
        loomR   yes   https://mojaveazure.github.io/loomR

  The package we suggest, loomR, is currently underdevelopment and not yet available on CRAN. This package is not required for core functionality of Seurat.

* checking package dependencies ... NOTE
  Package suggested but not available for checking: 'loomR'

  This is a suggested package hosted on a custom repository and maintained by us (both the package and repository).

## Downstream dependencies

There is a package that suggests Seurat (clustree), but this update does not impact their functionality.
