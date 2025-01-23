# Seurat v5.2.0

## Test environments
* local ubuntu 20.04 install, R 4.4.2
* win-builder (oldrelease, release, devel)
* mac-builder (devel)

We were unable to test on r-release on mac-builder because the portal seemed to point to the wrong version.

## R CMD check results

There were no ERRORs or WARNINGs

There was one NOTE

❯ checking CRAN incoming feasibility ... [12s/61s] NOTE
  Maintainer: ‘Rahul Satija <seurat@nygenome.org>’
  
  Suggests or Enhances not in mainstream repositories:
    BPCells, enrichR, presto
  Availability using Additional_repositories specification:
    BPCells   yes   https://bnprks.r-universe.dev   
    enrichR   yes   https://cran.r-universe.dev     
    presto    yes   https://satijalab.r-universe.dev

BPCells, enrichR, and presto are hosted on R-universe and used conditionally in Seurat.

## Downstream dependencies

There are 3 packages that depend on Seurat: CACIMAR, scCustomize, and SCdeconR; this update does not impact their functionality.

There are 34 packages that import Seurat: AnanseSeurat, APackOfTheClones, bbknnR, CAESAR.Suite, CAMML, DR.SC, DWLS, GeneNMF, ggsector, mixhvg, nebula, Platypus, PoweREST, PRECAST, ProFAST, rPanglaoDB, scAnnotate, scaper, sccca, scDiffCom, scGate, scGOclust, SCIntRuler, scMappR, scperturbR, scpoisson, SCRIP, scRNAstat, SignacX, SoupX, SpaCCI, SPECK, STREAK, and tidyseurat; this update does not impact their functionality.

There are 27 packages that suggest Seurat: BisqueRNA, Canek, cellpypes, CIARA, ClustAssess, clustree, combiroc, conos, countland, CRMetrics, CytoSimplex, DIscBIO, dyngen, easybio, grandR, harmony, laminr, mxfda, RESET, rliger, SCORPIUS, SCpubr, scregclust, Signac, SuperCell, treefit, and VAM; this update does not impact their functionality.
