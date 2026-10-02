# Seurat v5.6.0

## Test environments
* Ubuntu 20.04 (local) (R 4.5.2)
* Ubuntu 24.04 (GitHub Actions Runner): R-oldrelease, R-release, R-devel
* macOS Sequoia 15.6.1 (local) (R 4.6.1)
* Windows 11 x64 (local) (R 4.6.1)
* [macos-builder](https://mac.r-project.org/macbuilder/submit.html): R-devel
* [win-builder](https://win-builder.r-project.org/): R-oldrelease, R-release, R-devel

## R CMD check results

**Status: OK**

**1 NOTE** (from all Windows/win-builder only; other systems show no NOTEs)

```
Maintainer: 'Rahul Satija <seurat@nygenome.org>'

Possibly misspelled words in DESCRIPTION:
  Basu (4:370)
  Gennert (4:311)
  Hao (4:506, 4:511)
  Macosko (4:359)
  Satija (4:290, 4:378)
  al (4:325, 4:391, 4:458, 4:519)
  et (4:322, 4:388, 4:455, 4:516)
  transcriptomic (4:205)

Suggests or Enhances not in mainstream repositories:
  BPCells
Availability using Additional_repositories specification:
  BPCells   yes   https://bnprks.r-universe.dev   
  ?           ?   https://satijalab.r-universe.dev

Found the following (possibly) invalid URLs:
  URL: https://stackoverflow.com/questions/3942878/how-to-decide-font-color-in-white-or-black-depending-on-background-color
    From: man/BGTextColor.Rd
          man/contrast-theory.Rd
    Status: 403
    Message: Forbidden
```

- The maintainer remains Rahul Satija and the email is correct.
- BPCells and presto are hosted on R-universe and used conditionally in Seurat.
- All words/names in DESCRIPTION are correctly spelled.
- The URL can be accessed normally manually so does not seem to be invalid.

## Reverse dependency check results

We checked 89 reverse dependencies and found no regressions with the new version of Seurat.