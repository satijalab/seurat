# Seurat v3.1.0

## Test environments
* local Ubuntu 16.04.6 and 18.04.2 installs, R 3.5.3
* Ubuntu 16.04.6 (on travis-ci), R 3.6.1
* macOS 10.13.3 (on travis-ci), R 3.6.1
* Windows Server 2012 R2 (on AppVeyor), R 3.6.1 Patched
* win-builder (oldrelease, release, devel)

## R CMD check results
There were no ERRORs or WARNINGs

There were 3 NOTEs:

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
  
* checking Rd cross-references ... NOTE
  Package unavailable to check Rd xrefs: 'loomR'
  
  This is a suggested package that we maintain. We have checked the cross references to ensure they link to the correct help pages.

## Downstream dependencies

There is a package that suggests Seurat (clustree), but this update does not impact their functionality.
