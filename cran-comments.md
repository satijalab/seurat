# Seurat v3.0.0

## Test environments
* local Ubuntu 16.04.6 and 18.04.2 installs, R 3.5.3
* Ubuntu 14.04.5 (on travis-ci), R 3.5.3
* macOS 10.13.3 (on travis-ci), R 3.5.3
* Windows Server 2012 (on AppVeyor), R 3.5.3 Patched
* win-builder (devel)

## R CMD check results
There were no ERRORs or WARNINGs

There were 2 NOTEs:

* checking CRAN incoming feasibility ... NOTE
  Maintainer: 'Paul Hoffman <nygcSatijalab@nygenome.org>'

  New maintainer:
    Paul Hoffman <nygcSatijalab@nygenome.org>
  Old maintainer(s):
    Paul Hoffman <seuratpackage@gmail.com>

  Suggests or Enhances not in mainstream repositories:
    loomR
  Availability using Additional_repositories specification:
    loomR   yes   https://mojaveazure.github.io/loomR

  We are changing our maintainer email to separate our support email from our R package email. The package we suggest, loomR, is currently underdevelopment and not yet available on CRAN. This package is not required for core functionality of Seurat.

* checking package dependencies ... NOTE
  Package suggested but not available for checking: 'loomR'

  This is a suggested package hosted on a custom repository and maintained by us (both the package and repository).
  
## Downstream dependencies

There are currently no downstream dependencies for this package. There is a package that suggests Seurat (clustree), and we have reached out to them to prepare for our update.
