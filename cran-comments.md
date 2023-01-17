# Seurat v4.2.1

## Test environments
* local Ubuntu 20.04 install, R 4.1.3
* win-builder (release, devel)

## R CMD check results
There were no ERRORs or WARNINGs

There was one NOTE:

* checking CRAN incoming feasibility ... NOTE
Maintainer: 'Paul Hoffman <seurat@nygenome.org>'

Found the following (possibly) invalid URLs:
  URL: https://www.science.org/doi/abs/10.1126/science.aad0501
    From: man/cc.genes.Rd
          man/cc.genes.updated.2019.Rd
    Status: 503
    Message: Service Unavailable

This URL is valid and the service still exists. When navigating to the URL either via the documentation or directly, you are taken to the correct article

## Downstream dependencies

There no packages that depend on Seurat

There are sixteen packages that import Seurat: CAMML, CIDER, DR.SC, DUBStepR, maple, Platypus, rPanglaoDB, scDiffCom, scMappR, SCRIP, scRNAstat, Signac, SignacX, SoupX, spruce, and tidyseurat; this update does not impact their functionality

There are twelve packages that suggest Seurat: BisqueRNA, CIARA, ClustAssess, clustree, conos, DIscBIO, dyngen, harmony, rliger, Rmagic, treefit, and VAM; this update does not impact their functionality.
