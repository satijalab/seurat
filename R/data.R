#' Cell cycle genes
#'
#' A list of genes used in cell-cycle regression
#'
#' @format A list of two vectors
#' \describe{
#'   \item{s.genes}{Genes associated with S-phase}
#'   \item{g2m.genes}{Genes associated with G2M-phase}
#' }
#' @source \url{http://science.sciencemag.org/content/352/6282/189}
#'
"cc.genes"

#' Cell cycle genes: 2019 update
#'
#' A list of genes used in cell-cycle regression, updated with 2019 symbols
#'
#' @section Updated symbols:
#' The following symbols were updated from \code{\link{cc.genes}}
#' \describe{
#'   \item{s.genes}{
#'     \itemize{
#'       \item \emph{MCM2}: \emph{MCM7}
#'       \item \emph{MLF1IP}: \emph{CENPU}
#'       \item \emph{RPA2}: \emph{POLR1B}
#'       \item \emph{BRIP1}: \emph{MRPL36}
#'     }
#'   }
#'   \item{g2m.genes}{
#'     \itemize{
#'       \item \emph{FAM64A}: \emph{PIMREG}
#'       \item \emph{HN1}: \emph{JPT1}
#'     }
#'   }
#' }
#'
#' @format A list of two vectors
#' \describe{
#'   \item{s.genes}{Genes associated with S-phase}
#'   \item{g2m.genes}{Genes associated with G2M-phase}
#' }
#' @source \url{http://science.sciencemag.org/content/352/6282/189}
#'
#' @seealso \code{\link{cc.genes}}
#'
#' @examples
#' \dontrun{
#' cc.genes.updated.2019 <- cc.genes
#' cc.genes.updated.2019$s.genes <- UpdateSymbolList(symbols = cc.genes.updated.2019$s.genes)
#' cc.genes.updated.2019$g2m.genes <- UpdateSymbolList(symbols = cc.genes.updated.2019$g2m.genes)
#' }
#'
"cc.genes.updated.2019"

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
#' @source \url{https://support.10xgenomics.com/single-cell-gene-expression/datasets/1.1.0/pbmc3k}
#'
"pbmc_small"
