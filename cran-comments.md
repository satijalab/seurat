# Seurat v5.0.3

## Test environments
* local ubuntu 20.04 install, R 4.3.2
* win-builder (oldrelease, release, devel)
* mac-builder (release, devel)

## R CMD check results

There were no ERRORs or WARNINGs

There was one NOTE

> Suggests or Enhances not in mainstream repositories:
>  BPCells, presto
> Availability using Additional_repositories specification:
>  presto    yes   https://satijalab.r-universe.dev
>  BPCells   yes   https://satijalab.r-universe.dev

BPCells and presto are hosted on R-universe and used conditionally in Seurat.


## Downstream dependencies

There are two packages that depend on Seurat: CACIMAR and scCustomize; this update does not impact their functionality

There are 31 packages that import Seurat: AnanseSeurat, APackOfTheClones, bbknnR, CAMML, DR.SC, DWLS, GeneNMF, ggsector, mixhvg, nebula, Platypus, PRECAST, ProFAST, rPanglaoDB, scAnnotate, scaper, sccca, scDiffCom, scfetch, scGate, scGOclust, scMappR, scperturbR, scpoisson, SCRIP, scRNAstat, SignacX, SoupX, SPECK, STREAK, and tidyseurat; this update does not impact their functionality

There are 22 packages that suggest Seurat: BisqueRNA, Canek, cellpypes, CIARA, ClustAssess, clustree, combiroc, conos, countland, CRMetrics, CytoSimplex, DIscBIO, dyngen, grandR, harmony, RESET, rliger, SCORPIUS, SCpubr, Signac, treefit, and VAM; this update does not impact their functionality
