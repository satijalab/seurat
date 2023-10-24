# Seurat v5.0.0

## Test environments
* local ubuntu 20.04 install, R 4.1.3
* win-builder (oldrelease, release, devel)
* mac-builder (release)

## R CMD check results

There were no ERRORs or WARNINGs

There were two NOTEs

> * checking CRAN incoming feasibility ... NOTE
> Maintainer: 'Rahul Satija <seurat@nygenome.org>'

> New maintainer:
>  Rahul Satija <seurat@nygenome.org>
> Old maintainer(s):
>  Paul Hoffman <seurat@nygenome.org>

> Suggests or Enhances not in mainstream repositories:
>  BPCells, presto
> Availability using Additional_repositories specification:
>  presto    yes   https://satijalab.r-universe.dev
>  BPCells    no   ?
>  ?           ?   https://bnprks.r-universe.dev
> Additional repositories with no packages:
>  https://bnprks.r-universe.dev

> Packages suggested but not available for checking: 'BPCells', 'presto'

BPCells and presto are hosted on R-universe and used conditionally in Seurat.


## Downstream dependencies

The following reverse dependencies are impacted by this release of Seurat:

- AnanseSeurat
    - Failure in examples and test because of changes in how default objects are created in Seurat. 
    - Functionality impacted. The author was made aware of these changes: https://github.com/JGASmits/AnanseSeurat/issues/34

- CAMML
    - Failure in example because of changes in how default objects are created in Seurat.
    - Functionality impacted. The author was made aware of these changes over email

- Canek
    - Failure in tests because of changes in how default objects are created in Seurat.
    - Functionality impacted. The author was made aware of these changes: https://github.com/MartinLoza/Canek/issues/20

- clustree:
    - Failure in tests because of changes in data accessor methods in Seurat.
    - Functionality impacted. The author was made aware of these changes: https://github.com/lazappi/clustree/issues/93
    - Author has incorporated changes

- CSCDRNA
    - Failure in tests because of changes in data accessor methods in Seurat.
    - Functionality impacted. The author was made aware of these changes: https://github.com/empiricalbayes/CSCDRNA/issues/1 

- scCustomize
    - Failure in example because of changes in how default objects are created in Seurat.
    - Functionality impacted. The author was made aware of these changes: https://github.com/samuel-marsh/scCustomize/issues/131

- SCpubr:
    - Failure in example because of changes in how default objects are created in Seurat.
    - Functionality impacted. The author was made aware of these changes: https://github.com/enblacar/SCpubr/issues/42

- Signac
    - Faulure in new tests because of SeuratObject changing the order of the results, but not the actual values. 
    - Functionality not impacted. The author was made aware of these changes over email and has made changes.

- tidyseurat
      - Faulure in new tests because of SeuratObject changing the order of the results, but not the actual values. 
      - Functionality not impacted. The author was made aware of these changes: https://github.com/stemangiola/tidyseurat/issues/74
- VAM
    - Failure in tests because of changes in data accessor methods in Seurat.
    - Functionality impacted. The author was made aware of these changes over email
