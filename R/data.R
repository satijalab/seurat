#' Cell cycle genes
#'
#' A list of genes used in cell-cycle regression
#'
#' @format A list of two vectors
#' \describe{
#'   \item{s.genes}{Genes associated with S-phase}
#'   \item{g2m.genes}{Genes associated with G2M-phase}
#' }
#' @source http://science.sciencemag.org/content/352/6282/189
"cc.genes"

#' A small example version of the PBMC dataset
#'
#' A subsetted version of 10X Genomics' 3k PBMC dataset
#'
#' @format A Seurat object with the following slots filled
#' \describe{
#'   \item{assays}{
#'   \itemize{Currently only contains one assay ("RNA" - scRNA-seq expression data)
#'   \item{counts - Raw expression data}
#'   \item{data - Normalized expression data}
#'   \item{scale.data - Scaled expression data}
#'   \item{var.features - names of the current features selected as variable}
#'   \item{meta.features - Assay level metadata such as mean and variance}
#'    }}
#'   \item{meta.data}{Cell level metadata}
#'   \item{active.assay}{Current default assay}
#'   \item{active.ident}{Current default idents}
#'   \item{graphs}{Neighbor graphs computed, currently stores the SNN}
#'   \item{reductions}{Dimensional reductions: currently PCA and tSNE}
#'   \item{version}{Seurat version used to create the object}
#'   \item{commands}{Command history}
#' }
#' @source https://support.10xgenomics.com/single-cell-gene-expression/datasets/1.1.0/pbmc3k
"pbmc_small"
