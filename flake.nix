{
  description = "Seurat development environment";

  inputs = {
    nixpkgs.url = "github:nixos/nixpkgs?ref=nixos-unstable";
  };

  outputs = { self, nixpkgs }:
    let
      system = "x86_64-linux";
      pkgs = nixpkgs.legacyPackages.${system};
    in
      {
        devShells.${system}.default = pkgs.mkShell {
          buildInputs = with pkgs; [
            (rWrapper.override {
              packages = with rPackages; [
                # Development tools
                devtools
                testthat
                roxygen2
                
                # For testing - SeuratData will be installed via devtools::install_github
                # during test execution since it's not in nixpkgs
                
                # Core dependencies (from Imports in DESCRIPTION)
                cluster
                cowplot
                fastDummies
                fitdistrplus
                future
                future_apply
                generics
                ggplot2
                ggrepel
                ggridges
                httr
                ica
                igraph
                irlba
                jsonlite
                KernSmooth
                leidenbase
                lifecycle
                lmtest
                MASS
                Matrix
                matrixStats
                miniUI
                patchwork
                pbapply
                plotly
                png
                progressr
                RANN
                RColorBrewer
                Rcpp
                RcppAnnoy
                RcppEigen
                RcppHNSW
                RcppProgress
                reticulate
                rlang
                ROCR
                RSpectra
                Rtsne
                scales
                scattermore
                sctransform
                shiny
                spatstat_explore
                spatstat_geom
                tibble
                uwot
                
                # SeuratObject dependency
                SeuratObject
                
                # Additional packages for testing/development (from Suggests)
                # Only including packages available in nixpkgs
                ape
                arrow
                base64enc
                Biobase
                BiocGenerics
                data_table
                DESeq2
                DelayedArray
                GenomicRanges
                GenomeInfoDb
                glmGamPoi
                hdf5r
                IRanges
                limma
                magrittr
                MAST
                metap
                mixtools
                rsvd
                rtracklayer
                S4Vectors
                sf
                sp
                SingleCellExperiment
                SummarizedExperiment
                VGAM
                
                # tidyverse for the test script
                tidyverse
                readr
              ];
            })
          ];
          
          shellHook = ''
            echo "Seurat development environment loaded"
            echo "R version: $(R --version | head -n1)"
            echo ""
            echo "To build and test Seurat:"
            echo "  R CMD build ."
            echo "  R CMD check seurat_*.tar.gz"
            echo ""
            echo "Or in R:"
            echo "  devtools::load_all()"
            echo "  devtools::test()"
            echo ""
            echo "To test the SCTransform do.call bug fix:"
            echo "  Rscript test_sctransform_issue.R"
            echo "  (Auto-installs SeuratData and celegans.embryo via dev_mode)"
          '';
        };
      };
}
