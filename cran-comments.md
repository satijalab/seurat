## Test environments
* local Ubuntu 16.04 install, R 3.4.4
* Ubuntu 12.04 (on travis-ci), R 3.4.4
* Windows Server 2012 (on appveyor), R 3.4.4
* win-builder (devel)

## R CMD check results
There were no ERRORs or WARNINGs

There were 2 NOTEs:

* checking CRAN incoming feasibility ... NOTE
  Maintainer: ‘Paul Hoffman <seuratpackage@gmail.com>’
  Uses the superseded package: ‘doSNOW’

  We use doSNOW instead of doParallel because doSNOW has support for a progress bar while doParallel does not. doSNOW is still on CRAN and works for our purposes

* checking installed package size ... NOTE
  installed size is  6.7Mb
  sub-directories of 1Mb or more:
    libs   5.7Mb

  The large size is due to compiled code. All other directories in our package are under 500K each and the rest of our package is approximately 1.5M in size.

## Downstream dependencies

There are currently no downstream dependencies for this package
