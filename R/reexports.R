
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Classes
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' The Assay Class
#'
#' The \code{Assay} object is the basic unit of Seurat; for more details, please
#' see the documentation in \code{\link[SeuratObject:Assay]{SeuratObject}}
#'
#' @importClassesFrom SeuratObject Assay
#'
#' @exportClass Assay
#'
#' @docType class
#' @name Assay-class
#' @rdname Assay-class
#'
#' @seealso \code{\link[SeuratObject:Assay]{SeuratObject::Assay-class}}
#'
NULL

#' The DimReduc Class
#'
#' The \code{DimReduc} object stores a dimensionality reduction taken out in
#' Seurat; for more details, please see the documentation in
#' \code{\link[SeuratObject:DimReduc]{SeuratObject}}
#'
#' @importClassesFrom SeuratObject DimReduc
#'
#' @exportClass DimReduc
#'
#' @docType class
#' @name DimReduc-class
#' @rdname DimReduc-class
#'
#' @seealso \code{\link[SeuratObject:DimReduc]{SeuratObject::DimReduc-class}}
#'
NULL

#' The Graph Class
#'
#' For more details, please see the documentation in
#' \code{\link[SeuratObject:Graph]{SeuratObject}}
#'
#' @importClassesFrom SeuratObject Graph
#'
#' @exportClass Graph
#'
#' @docType class
#' @name Graph-class
#' @rdname Graph-class
#'
#' @seealso \code{\link[SeuratObject:Graph]{SeuratObject::Graph-class}}
#'
NULL

#' The JackStrawData Class
#'
#' For more details, please see the documentation in
#' \code{\link[SeuratObject:JackStrawData]{SeuratObject}}
#'
#' @importClassesFrom SeuratObject JackStrawData
#'
#' @exportClass JackStrawData
#'
#' @docType class
#' @name JackStrawData-class
#' @rdname JackStrawData-class
#'
#' @seealso \code{\link[SeuratObject:JackStrawData]{SeuratObject::JackStrawData-class}}
#'
NULL

#' The Neighbor Class
#'
#' For more details, please see the documentation in
#' \code{\link[SeuratObject:Neighbor]{SeuratObject}}
#'
#' @importClassesFrom SeuratObject Neighbor
#'
#' @exportClass Neighbor
#'
#' @docType class
#' @name Neighbor-class
#' @rdname Neighbor-class
#'
#' @seealso \code{\link[SeuratObject:Neighbor]{SeuratObject::Neighbor-class}}
#'
NULL

#' The Seurat Class
#'
#' The Seurat object is a representation of single-cell expression data for R;
#' for more details, please see the documentation in
#' \code{\link[SeuratObject:Seurat]{SeuratObject}}
#'
#' @importClassesFrom SeuratObject Seurat
#'
#' @exportClass Seurat
#'
#' @docType class
#' @name Seurat-class
#' @rdname Seurat-class
#'
#' @seealso \code{\link[SeuratObject:Seurat]{SeuratObject::Seurat-class}}
#'
NULL

#' The SeuratCommand Class
#'
#' For more details, please see the documentation in
#' \code{\link[SeuratObject:SeuratCommand]{SeuratObject}}
#'
#' @importClassesFrom SeuratObject SeuratCommand
#'
#' @exportClass SeuratCommand
#'
#' @docType class
#' @name SeuratCommand-class
#' @rdname SeuratCommand-class
#'
#' @seealso \code{\link[SeuratObject:SeuratCommand]{SeuratObject::SeuratCommand-class}}
#'
NULL

#' The SpatialImage Class
#'
#' For more details, please see the documentation in
#' \code{\link[SeuratObject:SpatialImage]{SeuratObject}}
#'
#' @importClassesFrom SeuratObject SpatialImage
#'
#' @exportClass SpatialImage
#'
#' @docType class
#' @name SpatialImage-class
#' @rdname SpatialImage-class
#'
#' @seealso \code{\link[SeuratObject:SpatialImage]{SeuratObject::SpatialImage-class}}
#'
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Functions and Generics
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @importFrom generics components
#' @rdname reexports
#' @export
#'
generics::components

#' @importFrom SeuratObject %||%
#' @rdname reexports
#' @export
#'
SeuratObject::`%||%`

#' @importFrom SeuratObject %iff%
#' @rdname reexports
#' @export
#'
SeuratObject::`%iff%`

#' @importFrom SeuratObject AddMetaData
#' @export
#'
SeuratObject::AddMetaData

#' @importFrom SeuratObject as.Graph
#' @export
#'
SeuratObject::as.Graph

#' @importFrom SeuratObject as.Neighbor
#' @export
#'
SeuratObject::as.Neighbor

#' @importFrom SeuratObject as.Seurat
#' @export
#'
SeuratObject::as.Seurat

#' @importFrom SeuratObject as.sparse
#' @export
#'
SeuratObject::as.sparse

#' @importFrom SeuratObject Assays
#' @export
#'
SeuratObject::Assays

#' @importFrom SeuratObject Cells
#' @export
#'
SeuratObject::Cells

#' @importFrom SeuratObject CellsByIdentities
#' @export
#'
SeuratObject::CellsByIdentities

#' @importFrom SeuratObject Command

#' @export
#'
SeuratObject::Command

#' @importFrom SeuratObject CreateAssayObject

#' @export
#'
SeuratObject::CreateAssayObject

#' @importFrom SeuratObject CreateDimReducObject
#' @export
#'
SeuratObject::CreateDimReducObject

#' @importFrom SeuratObject CreateSeuratObject

#' @export
#'
SeuratObject::CreateSeuratObject

#' @importFrom SeuratObject DefaultAssay
#' @export
#'
SeuratObject::DefaultAssay

#' @importFrom SeuratObject DefaultAssay<-
#' @export
#'
SeuratObject::`DefaultAssay<-`

#' @importFrom SeuratObject Distances
#' @export
#'
SeuratObject::Distances

#' @importFrom SeuratObject Embeddings
#' @export
#'
SeuratObject::Embeddings

#' @importFrom SeuratObject FetchData
#' @export
#'
SeuratObject::FetchData

#' @importFrom SeuratObject GetAssayData
#' @export
#'
SeuratObject::GetAssayData

#' @importFrom SeuratObject GetImage
#' @export
#'
SeuratObject::GetImage

#' @importFrom SeuratObject GetTissueCoordinates
#' @export
#'
SeuratObject::GetTissueCoordinates

#' @importFrom SeuratObject HVFInfo
#' @export
#'
SeuratObject::HVFInfo

#' @importFrom SeuratObject Idents
#' @export
#'
SeuratObject::Idents

#' @importFrom SeuratObject Idents<-
#' @export
#'
SeuratObject::`Idents<-`

#' @importFrom SeuratObject Images
#' @export
#'
SeuratObject::Images

#' @importFrom SeuratObject Index
#' @export
#'
SeuratObject::Index

#' @importFrom SeuratObject Index<-
#' @export
#'
SeuratObject::`Index<-`

#' @importFrom SeuratObject Indices
#' @export
#'
SeuratObject::Indices

#' @importFrom SeuratObject IsGlobal
#' @export
#'
SeuratObject::IsGlobal

#' @importFrom SeuratObject JS
#' @export
#'
SeuratObject::JS

#' @importFrom SeuratObject JS<-

#' @export
#'
SeuratObject::`JS<-`

#' @importFrom SeuratObject Key
#' @export
#'
SeuratObject::Key

#' @importFrom SeuratObject Key<-

#' @export
#'
SeuratObject::`Key<-`

#' @importFrom SeuratObject Loadings

#' @export
#'
SeuratObject::Loadings

#' @importFrom SeuratObject Loadings<-
#' @export
#'
SeuratObject::`Loadings<-`

#' @importFrom SeuratObject LogSeuratCommand

#' @export
#'
SeuratObject::LogSeuratCommand

#' @importFrom SeuratObject Misc
#' @export
#'
SeuratObject::Misc

#' @importFrom SeuratObject Misc<-
#' @export
#'
SeuratObject::`Misc<-`

#' @importFrom SeuratObject Neighbors
#' @export
#'
SeuratObject::Neighbors

#' @importFrom SeuratObject Project
#' @export
#'
SeuratObject::Project

#' @importFrom SeuratObject Project<-
#' @export
#'
SeuratObject::`Project<-`

#' @importFrom SeuratObject Radius
#' @export
#'
SeuratObject::Radius

#' @importFrom SeuratObject Reductions
#' @export
#'
SeuratObject::Reductions

#' @importFrom SeuratObject RenameCells
#' @export
#'
SeuratObject::RenameCells

#' @importFrom SeuratObject RenameIdents
#' @export
#'
SeuratObject::RenameIdents

#' @importFrom SeuratObject ReorderIdent
#' @export
#'
SeuratObject::ReorderIdent

#' @importFrom SeuratObject RowMergeSparseMatrices
#' @export
#'
SeuratObject::RowMergeSparseMatrices

#' @importFrom SeuratObject SetAssayData
#' @export
#'
SeuratObject::SetAssayData

#' @importFrom SeuratObject SetIdent
#' @export
#'
SeuratObject::SetIdent

#' @importFrom SeuratObject SpatiallyVariableFeatures
#' @export
#'
SeuratObject::SpatiallyVariableFeatures

#' @importFrom SeuratObject StashIdent
#' @export
#'
SeuratObject::StashIdent

#' @importFrom SeuratObject Stdev
#' @export
#'
SeuratObject::Stdev

#' @importFrom SeuratObject SVFInfo
#' @export
#'
SeuratObject::SVFInfo

#' @importFrom SeuratObject Tool
#' @export
#'
SeuratObject::Tool

#' @importFrom SeuratObject Tool<-
#' @export
#'
SeuratObject::`Tool<-`

#' @importFrom SeuratObject UpdateSeuratObject
#' @export
#'
SeuratObject::UpdateSeuratObject

#' @importFrom SeuratObject VariableFeatures

#' @export
#'
SeuratObject::VariableFeatures

#' @importFrom SeuratObject VariableFeatures<-
#' @export
#'
SeuratObject::`VariableFeatures<-`

#' @importFrom SeuratObject WhichCells
#' @export
#'
SeuratObject::WhichCells
