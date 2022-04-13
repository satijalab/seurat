
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
#' @rdname reexports
#' @export
#'
SeuratObject::AddMetaData

#' @importFrom SeuratObject as.Graph
#' @rdname reexports
#' @export
#'
SeuratObject::as.Graph

#' @importFrom SeuratObject as.Neighbor
#' @rdname reexports
#' @export
#'
SeuratObject::as.Neighbor

#' @importFrom SeuratObject as.Seurat
#' @rdname reexports
#' @export
#'
SeuratObject::as.Seurat

#' @importFrom SeuratObject as.sparse
#' @rdname reexports
#' @export
#'
SeuratObject::as.sparse

#' @importFrom SeuratObject Assays
#' @rdname reexports
#' @export
#'
SeuratObject::Assays

#' @importFrom SeuratObject Cells
#' @rdname reexports
#' @export
#'
SeuratObject::Cells

#' @importFrom SeuratObject CellsByIdentities
#' @rdname reexports
#' @export
#'
SeuratObject::CellsByIdentities

#' @importFrom SeuratObject Command
#' @rdname reexports
#' @export
#'
SeuratObject::Command

#' @importFrom SeuratObject CreateAssayObject
#' @rdname reexports
#' @export
#'
SeuratObject::CreateAssayObject

#' @importFrom SeuratObject CreateDimReducObject
#' @rdname reexports
#' @export
#'
SeuratObject::CreateDimReducObject

#' @importFrom SeuratObject CreateSeuratObject
#' @rdname reexports
#' @export
#'
SeuratObject::CreateSeuratObject

#' @importFrom SeuratObject DefaultAssay
#' @rdname reexports
#' @export
#'
SeuratObject::DefaultAssay

#' @importFrom SeuratObject DefaultAssay<-
#' @rdname reexports
#' @export
#'
SeuratObject::`DefaultAssay<-`

#' @importFrom SeuratObject Distances
#' @rdname reexports
#' @export
#'
SeuratObject::Distances

#' @importFrom SeuratObject Embeddings
#' @rdname reexports
#' @export
#'
SeuratObject::Embeddings

#' @importFrom SeuratObject FetchData
#' @rdname reexports
#' @export
#'
SeuratObject::FetchData

#' @importFrom SeuratObject GetAssayData
#' @rdname reexports
#' @export
#'
SeuratObject::GetAssayData

#' @importFrom SeuratObject GetImage
#' @rdname reexports
#' @export
#'
SeuratObject::GetImage

#' @importFrom SeuratObject GetTissueCoordinates
#' @rdname reexports
#' @export
#'
SeuratObject::GetTissueCoordinates

#' @importFrom SeuratObject HVFInfo
#' @rdname reexports
#' @export
#'
SeuratObject::HVFInfo

#' @importFrom SeuratObject Idents
#' @rdname reexports
#' @export
#'
SeuratObject::Idents

#' @importFrom SeuratObject Idents<-
#' @rdname reexports
#' @export
#'
SeuratObject::`Idents<-`

#' @importFrom SeuratObject Images
#' @rdname reexports
#' @export
#'
SeuratObject::Images

#' @importFrom SeuratObject Index
#' @rdname reexports
#' @export
#'
SeuratObject::Index

#' @importFrom SeuratObject Index<-
#' @rdname reexports
#' @export
#'
SeuratObject::`Index<-`

#' @importFrom SeuratObject Indices
#' @rdname reexports
#' @export
#'
SeuratObject::Indices

#' @importFrom SeuratObject IsGlobal
#' @rdname reexports
#' @export
#'
SeuratObject::IsGlobal

#' @importFrom SeuratObject JS
#' @rdname reexports
#' @export
#'
SeuratObject::JS

#' @importFrom SeuratObject JS<-
#' @rdname reexports
#' @export
#'
SeuratObject::`JS<-`

#' @importFrom SeuratObject Key
#' @rdname reexports
#' @export
#'
SeuratObject::Key

#' @importFrom SeuratObject Key<-
#' @rdname reexports
#' @export
#'
SeuratObject::`Key<-`

#' @importFrom SeuratObject Loadings
#' @rdname reexports
#' @export
#'
SeuratObject::Loadings

#' @importFrom SeuratObject Loadings<-
#' @rdname reexports
#' @export
#'
SeuratObject::`Loadings<-`

#' @importFrom SeuratObject LogSeuratCommand
#' @rdname reexports
#' @export
#'
SeuratObject::LogSeuratCommand

#' @importFrom SeuratObject Misc
#' @rdname reexports
#' @export
#'
SeuratObject::Misc

#' @importFrom SeuratObject Misc<-
#' @rdname reexports
#' @export
#'
SeuratObject::`Misc<-`

#' @importFrom SeuratObject Neighbors
#' @rdname reexports
#' @export
#'
SeuratObject::Neighbors

#' @importFrom SeuratObject Project
#' @rdname reexports
#' @export
#'
SeuratObject::Project

#' @importFrom SeuratObject Project<-
#' @rdname reexports
#' @export
#'
SeuratObject::`Project<-`

#' @importFrom SeuratObject Radius
#' @rdname reexports
#' @export
#'
SeuratObject::Radius

#' @importFrom SeuratObject Reductions
#' @rdname reexports
#' @export
#'
SeuratObject::Reductions

#' @importFrom SeuratObject RenameCells
#' @rdname reexports
#' @export
#'
SeuratObject::RenameCells

#' @importFrom SeuratObject RenameIdents
#' @rdname reexports
#' @export
#'
SeuratObject::RenameIdents

#' @importFrom SeuratObject ReorderIdent
#' @rdname reexports
#' @export
#'
SeuratObject::ReorderIdent

#' @importFrom SeuratObject RowMergeSparseMatrices
#' @rdname reexports
#' @export
#'
SeuratObject::RowMergeSparseMatrices

#' @importFrom SeuratObject SetAssayData
#' @rdname reexports
#' @export
#'
SeuratObject::SetAssayData

#' @importFrom SeuratObject SetIdent
#' @rdname reexports
#' @export
#'
SeuratObject::SetIdent

#' @importFrom SeuratObject SpatiallyVariableFeatures
#' @rdname reexports
#' @export
#'
SeuratObject::SpatiallyVariableFeatures

#' @importFrom SeuratObject StashIdent
#' @rdname reexports
#' @export
#'
SeuratObject::StashIdent

#' @importFrom SeuratObject Stdev
#' @rdname reexports
#' @export
#'
SeuratObject::Stdev

#' @importFrom SeuratObject SVFInfo
#' @rdname reexports
#' @export
#'
SeuratObject::SVFInfo

#' @importFrom SeuratObject Tool
#' @rdname reexports
#' @export
#'
SeuratObject::Tool

#' @importFrom SeuratObject Tool<-
#' @rdname reexports
#' @export
#'
SeuratObject::`Tool<-`

#' @importFrom SeuratObject UpdateSeuratObject
#' @rdname reexports
#' @export
#'
SeuratObject::UpdateSeuratObject

#' @importFrom SeuratObject VariableFeatures
#' @rdname reexports
#' @export
#'
SeuratObject::VariableFeatures

#' @importFrom SeuratObject VariableFeatures<-
#' @rdname reexports
#' @export
#'
SeuratObject::`VariableFeatures<-`

#' @importFrom SeuratObject WhichCells
#' @rdname reexports
#' @export
#'
SeuratObject::WhichCells
