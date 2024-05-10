# Seurat v5.1.0

## Test environments
* local ubuntu 20.04 install, R 4.3.2
* local macOS 14.1, R 4.4.0
* win-builder (oldrelease, release, devel)

## R CMD check results

There were no ERRORs or WARNINGs

There were two NOTEs

> Suggests or Enhances not in mainstream repositories:
>   BPCells, presto
> Availability using Additional_repositories specification:
>   BPCells   yes   https://bnprks.r-universe.dev   
>   presto    yes   https://satijalab.r-universe.dev

> * checking package dependencies ... NOTE
> Package suggested but not available for checking: 'BPCells'

BPCells is hosted on R-universe and used conditionally in SeuratObject.

## Downstream dependencies

There are two packages that depend on Seurat: CACIMAR and scCustomize; this update does not impact their functionality

There are 30 packages that import Seurat: AnanseSeurat, APackOfTheClones, bbknnR, CAMML, DR.SC, DWLS, GeneNMF, ggsector, mixhvg, nebula, Platypus, PRECAST, ProFAST, rPanglaoDB, scAnnotate, scaper, sccca, scDiffCom, scGate, scGOclust, scMappR, scperturbR, scpoisson, SCRIP, scRNAstat, SignacX, SoupX, SPECK, STREAK, and tidyseurat; this update does not impact their functionality

There are 22 packages that suggest Seurat: BisqueRNA, Canek, cellpypes, CIARA, ClustAssess, clustree, combiroc, conos, countland, CRMetrics, CytoSimplex, DIscBIO, dyngen, grandR, harmony, RESET, rliger, SCORPIUS, SCpubr, Signac, treefit, and VAM; this update does not impact their functionality
