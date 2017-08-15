#' Deprecated function(s) in the Seurat package
#'
#' These functions are provided for compatibility with older version of the Seurat package.  They may eventually be completely removed.
#' @rdname Seurat-deprecated
#' @name Seurat-deprecated
#' @param ... Parameters to be passed to the modern version of the function
#' @export vlnPlot subsetData pca PCA project.pca viz.pca set.ident pca.plot pcHeatmap jackStraw jackStrawPlot run_tsne tsne.plot find.markers find_all_markers genePlot feature.plot tsne.plot buildClusterTree plotClusterTree plotNoiseModel add_samples subsetCells project.samples run_diffusion ica ICA cluster.alpha average.pca average.expression icTopGenes pcTopGenes pcTopCells fetch.data viz.ica regulatorScore find.markers.node diffExp.test tobit.test batch.gene marker.test which.cells set.all.ident rename.ident posterior.plot map.cell get.centroids refined.mapping initial.mapping calc.insitu fit.gene.k fit.gene.mix addSmoothedScore addImputedScore getNewScore calcNoiseModels feature.plot.keynote feature.heatmap ica.plot spatial.de DBclust_dimension Kclust_dimension pca.sig.genes doHeatMap icHeatmap doKMeans genes.in.cluster kMeansHeatmap cell.cor.matrix gene.cor.matrix calinskiPlot dot.plot addMetaData removePC geneScorePlot cellPlot jackStraw.permutation.test jackStrawMC jackStrawFull writ.table jackRandom MeanVarPlot HeatmapNode minusr minusc RegressOut VizClassification
#' @aliases vlnPlot subsetData pca PCA project.pca viz.pca set.ident pca.plot pcHeatmap jackStraw jackStrawPlot run_tsne tsne.plot find.markers find_all_markers genePlot feature.plot tnse.plot buildClusterTree plotClusterTree plotNoiseModel add_samples subsetCells project.samples run_diffusion ica ICA cluster.alpha average.pca average.expression icTopGenes pcTopGenes pcTopCells fetch.data viz.ica regulatorScore find.markers.node diffExp.test tobit.test batch.gene marker.test which.cells set.all.ident rename.ident posterior.plot map.cell get.centroids refined.mapping initial.mapping calc.insitu fit.gene.k fit.gene.mix addSmoothedScore addImputedScore getNewScore calcNoiseModels feature.plot.keynote feature.heatmap ica.plot spatial.de DBclust_dimension Kclust_dimension pca.sig.genes doHeatMap icHeatmap doKMeans genes.in.cluster kMeansHeatmap cell.cor.matrix gene.cor.matrix calinskiPlot dot.plot addMetaData removePC geneScorePlot cellPlot jackStraw.permutation.test jackStrawMC jackStrawFull writ.table jackRandom MeanVarPlot HeatmapNode minusr minusc RegressOut VizClassification
#' @section Details:
#' \tabular{rl}{
#'   \code{vlnPlot} \tab now a synonym for \code{VlnPlot}\cr
#'   \code{subsetData} \tab now a synonym for \code{SubsetData}\cr
#'   \code{pca} \tab now a synonym for \code{RunPCA}\cr
#'   \code{PCA} \tab now a synonym for \code{PCA}\cr
#'   \code{project.pca} \tab now a synonym for \code{ProjectPCA}\cr
#'   \code{viz.pca} \tab now a synonym for \code{VizPCA}\cr
#'   \code{set.ident} \tab now a synonym for \code{SetIdent}\cr
#'   \code{pca.plot} \tab now a synonym for \code{PCAPlot}\cr
#'   \code{pcHeatmap} \tab now a synonym for \code{PCHeatmap}\cr
#'   \code{jackStraw} \tab now a synonym for \code{JackStraw}\cr
#'   \code{jackStrawPlot} \tab now a synonym for \code{JackStrawPlot}\cr
#'   \code{run_tsne} \tab now a synonym for \code{RunTSNE}\cr
#'   \code{tsne.plot} \tab now a synonym for \code{TSNEPlot}\cr
#'   \code{find.markers} \tab now a synonym for \code{FindMarkers}\cr
#'   \code{find_all_markers} \tab now a synonym for \code{FindAllMarkers}\cr
#'   \code{genePlot} \tab now a synonym for \code{GenePlot}\cr
#'   \code{feature.plot} \tab now a synonym for \code{FeaturePlot}\cr
#'   \code{buildClusterTree} \tab now a synonym for \code{BuildClusterTree}\cr
#'   \code{plotClusterTree} \tab now a synonym for \code{PlotClusterTree}\cr
#'   \code{plotNoiseModel} \tab has been removed and may be replaced at a later date\cr
#'   \code{add_samples} \tab now a synonym for \code{AddSamples}\cr
#'   \code{subsetCells} \tab now deleted\cr
#'   \code{project.samples} \tab has been removed and may be replaced at a later date\cr
#'   \code{run_diffusion} \tab now a synonym for \code{RunDiffusion}\cr
#'   \code{ica} \tab now a synonym for \code{RunICA}\cr
#'   \code{ICA} \tab now a synonym for \code{RunICA}\cr
#'   \code{cluster.alpha} \tab now a synonym for \code{AverageDetectionRate}\cr
#'   \code{average.pca} \tab now a synonym for \code{AveragePCA}\cr
#'   \code{average.expression} \tab now a synonym for \code{AverageExpression}\cr
#'   \code{icTopGenes} \tab now a synonym for \code{ICTopGenes}\cr
#'   \code{pcTopGenes} \tab now a synonym for \code{PCTopGenes}\cr
#'   \code{pcTopCells} \tab now a synonym for \code{PCTopCells}\cr
#'   \code{fetch.data} \tab now a synonym for \code{FetchData}\cr
#'   \code{viz.ica} \tab now a synonym for \code{VizIca}\cr
#'   \code{regulatorScore} \tab now deleted\cr
#'   \code{find.markers.node} \tab now a synonym for \code{FindMarkersNode}\cr
#'   \code{diffExp.test} \tab now a synonym for \code{DiffExpTest}\cr
#'   \code{tobit.test} \tab now a synonym for \code{TobitTest}\cr
#'   \code{batch.gene} \tab has been removed and may be restored at a later date\cr
#'   \code{marker.test} \tab now a synonym for \code{MarkerTest}\cr
#'   \code{which.cells} \tab now a synonym for \code{WhichCells}\cr
#'   \code{set.all.ident} \tab now a synonym for \code{SetAllIdent}\cr
#'   \code{rename.ident} \tab now a synonym for \code{RenameIdent}\cr
#'   \code{posterior.plot} \tab now a synonym for \code{PosteriorPlot}\cr
#'   \code{map.cell} \tab has been deprecated\cr
#'   \code{get.centroids} \tab now a synonym for \code{GetCentroids}\cr
#'   \code{refined.mapping} \tab now a synonym for \code{RefinedMapping}\cr
#'   \code{initial.mapping} \tab now a synonym for \code{InitialMapping}\cr
#'   \code{calc.insitu} \tab now a synonym for \code{CalcInsitu}\cr
#'   \code{fit.gene.k} \tab now a synonym for \code{FitGeneK}\cr
#'   \code{fit.gene.mix} \tab now a synonym for \code{FitGeneMix}\cr
#'   \code{addSmoothedScore} \tab now a synonym for \code{AddSmoothedScore}\cr
#'   \code{addImputedScore} \tab now a synonym for \code{AddImputedScore}\cr
#'   \code{getNewScore} \tab has been removed without replacement\cr
#'   \code{calcNoiseModels} \tab has been removed and may be replaced at a later date\cr
#'   \code{feature.plot.keynote} \tab has been removed without replacement\cr
#'   \code{feature.heatmap} \tab now a synonym for \code{FeatureHeatmap}\cr
#'   \code{ica.plot} \tab now a synonym for \code{ICAPlot}\cr
#'   \code{spatial.de} \tab has been removed and may be replaced at a later date\cr
#'   \code{DBclust_dimension} \tab now a synonym for \code{DBClustDimension}\cr
#'   \code{Kclust_dimension} \tab now a synonym for \code{KClustDimension}\cr
#'   \code{pca.sig.genes} \tab now a synonym for \code{PCASigGenes}\cr
#'   \code{doHeatMap} \tab now a synonym for \code{DoHeatMap}\cr
#'   \code{icHeatmap} \tab now a synonym for \code{ICHeatmap}\cr
#'   \code{doKMeans} \tab now a synonym for \code{DoKMeans}\cr
#'   \code{genes.in.cluster} \tab now a synonym for \code{GenesInCluster}\cr
#'   \code{kMeansHeatmap} \tab now a synonym for \code{KMeansHeatmap}\cr
#'   \code{cell.cor.matrix} \tab has been removed and may be replaced at a later date\cr
#'   \code{gene.cor.matrix} \tab has been removed and may be replaced at a later date\cr
#'   \code{calinskiPlot} \tab has been removed and may be replaced at a later date\cr
#'   \code{dot.plot} \tab now a synonym for \code{DotPlot}\cr
#'   \code{addMetaData} \tab now a synonym for \code{AddMetaData}\cr
#'   \code{removePC} \tab has been removed and may be replaced at a later date\cr
#'   \code{geneScorePlot} \tab now deleted\cr
#'   \code{cellPlot} \tab now a synonym for \code{CellPlot}\cr
#'   \code{jackStraw.permutation.test} \tab has been deleted\cr
#'   \code{jackStrawMC} \tab has been deleted\cr
#'   \code{jackStrawFull} \tab has been deleted\cr
#'   \code{PCAFast} \tab now a synonym for \code{PCA}\cr
#'   \code{writ.table} \tab has been removed without replacement\cr
#'   \code{jackRandom} \tab has been removed without replacement\cr
#'   \code{MeanVarPlot} \tab now a synonym for \code{FindVariableGenes}\cr
#'   \code{myPalette} \tab now a synonym for \code{CustomPalette}\cr
#'   \code{minusr} \tab now a synonym for \code{SubsetRow}\cr
#'   \code{minusc} \tab now a synonym for \code{SubsetColumn}\cr
#'   \code{RegressOut} \tab now part of \code{ScaleData}\cr
#'   \code{VizClassification} \tab has been removed without replacement\cr
#' }
#'
vlnPlot <- function(...) {
    .Deprecated("VlnPlot", package="Seurat")
    VlnPlot(...)
}

subsetData <- function(...) {
    .Deprecated("SubsetData", package="Seurat")
    SubsetData(...)
}

pca <- function(...) {
    .Deprecated("RunPCA", package="Seurat")
    RunPCA(...)
}

PCA <- function(...) {
  .Deprecated("RunPCA", package="Seurat")
  RunPCA(...)
}

project.pca <- function(...) {
    .Deprecated("ProjectPCA", package="Seurat")
    ProjectPCA(...)
}

viz.pca <- function(...) {
    .Deprecated("VizPCA", package="Seurat")
    VizPCA(...)
}

set.ident <- function(...) {
    .Deprecated("SetIdent", package="Seurat")
    SetIdent(...)
}

pca.plot <- function(...) {
    .Deprecated("PCAPlot", package="Seurat")
    PCAPlot(...)
}

pcHeatmap <- function(...) {
    .Deprecated("PCHeatmap", package="Seurat")
    PCHeatmap(...)
}

jackStraw <- function(...) {
    .Deprecated("JackStraw", package="Seurat")
    JackStraw(...)
}

jackStrawPlot <- function(...) {
    .Deprecated("JackStrawPlot", package="Seurat")
    JackStrawPlot(...)
}

run_tsne <- function(...) {
    .Deprecated("RunTSNE", package="Seurat")
    RunTSNE(...)
}

tsne.plot <- function(...) {
    .Deprecated("TSNEPlot", package="Seurat")
    TSNEPlot(...)
}

find.markers <- function(...) {
    .Deprecated("FindMarkers", package="Seurat")
    FindMarkers(...)
}

find_all_markers <- function(...) {
    .Deprecated("FindAllMarkers", package="Seurat")
    FindAllMarkers(...)
}

genePlot <- function(...) {
    .Deprecated("GenePlot", package="Seurat")
    GenePlot(...)
}

feature.plot <- function(...) {
    .Deprecated("FeaturePlot", package="Seurat")
    FeaturePlot(...)
}

buildClusterTree <- function(...) {
    .Deprecated("BuildClusterTree", package="Seurat")
    BuildClusterTree(...)
}

plotClusterTree <- function(...) {
    .Deprecated("PlotClusterTree", package="Seurat")
    PlotClusterTree(...)
}

plotNoiseModel <- function(...) {
  .Deprecated(
    package = 'Seurat',
    msg = 'plotNoiseModel has been removed and may be replaced at a later date'
  )
    # .Deprecated("PlotNoiseModel", package="Seurat")
    # PlotNoiseModel(...)
}

add_samples <- function(...) {
    .Deprecated("AddSamples", package="Seurat")
    AddSamples(...)
}

subsetCells <- function(...) {
    .Deprecated(
      package = "Seurat",
      msg = 'subsetCells is now deleted, please use SubsetData'
    )
}

project.samples <- function(...) {
  .Deprecated(
    package = 'Seurat',
    msg = 'project.samples has been removed and may be replaced at a later date'
  )
    # .Deprecated("ProjectSamples", package="Seurat")
    # ProjectSamples(...)
}

run_diffusion <- function(...) {
    .Deprecated("RunDiffusion", package="Seurat")
    RunDiffusion(...)
}

ica <- function(...) {
    .Deprecated("RunICA", package="Seurat")
    RunICA(...)
}

ICA <- function(...) {
  .Deprecated("RunICA", package="Seurat")
  RunICA(...)
}

cluster.alpha <- function(...) {
    .Deprecated("AverageDetectionRate", package="Seurat")
    AverageDetectionRate(...)
}

average.pca <- function(...) {
    .Deprecated("AveragePCA", package="Seurat")
    AveragePCA(...)
}

average.expression <- function(...) {
    .Deprecated("AverageExpression", package="Seurat")
    AverageExpression(...)
}

icTopGenes <- function(...) {
    .Deprecated("ICTopGenes", package="Seurat")
    ICTopGenes(...)
}

pcTopGenes <- function(...) {
    .Deprecated("PCTopGenes", package="Seurat")
    PCTopGenes(...)
}
pcTopCells <- function(...) {
    .Deprecated("PCTopCells", package="Seurat")
    PCTopCells(...)
}

fetch.data <- function(...) {
    .Deprecated("FetchData", package="Seurat")
    FetchData(...)
}

viz.ica <- function(...) {
    .Deprecated("VizICA", package="Seurat")
    VizICA(...)
}

regulatorScore <- function(...) {
    .Deprecated(
      package = "Seurat",
      msg = 'regulatorScore has been deleted without replacement'
    )
}

find.markers.node <- function(...) {
    .Deprecated("FindMarkersNode", package="Seurat")
    FindMarkersNode(...)
}

diffExp.test <- function(...) {
    .Deprecated("DiffExpTest", package="Seurat")
    DiffExpTest(...)
}

tobit.test <- function(...) {
    .Deprecated("TobitTest", package="Seurat")
    TobitTest(...)
}

batch.gene <- function(...) {
  .Deprecated(
    package = 'Seurat',
    msg = 'batch.gene has been removed and may be restored at a later date'
  )
    # .Deprecated("BatchGene", package="Seurat")
    # BatchGene(...)
}

marker.test <- function(...) {
    .Deprecated("MarkerTest ", package="Seurat")
    MarkerTest(...)
}

which.cells <- function(...) {
    .Deprecated("WhichCells", package="Seurat")
    WhichCells(...)
}

set.all.ident <- function(...) {
    .Deprecated("SetAllIdent", package="Seurat")
    SetAllIdent(...)
}

rename.ident <- function(...) {
    .Deprecated("RenameIdent", package="Seurat")
    RenameIdent(...)
}
posterior.plot <- function(...) {
    .Deprecated("PosteriorPlot", package="Seurat")
    PosteriorPlot(...)
}

map.cell <- function(...) {
  .Deprecated(
    package = 'Seurat',
    msg = 'map.cell has been removed and may be restored at a later date'
  )
    # .Deprecated("MapCell", package="Seurat")
    # MapCell(...)
}

get.centroids <- function(...) {
    .Deprecated("GetCentroids", package="Seurat")
    GetCentroids(...)
}

refined.mapping <- function(...) {
    .Deprecated("RefinedMapping", package="Seurat")
    RefinedMapping(...)
}

initial.mapping <- function(...) {
    .Deprecated("InitialMapping", package="Seurat")
    InitialMapping(...)
}

calc.insitu <- function(...) {
    .Deprecated("CalcInsitu ", package="Seurat")
    CalcInsitu(...)
}

fit.gene.k <- function(...) {
    .Deprecated("FitGeneK ", package="Seurat")
    FitGeneK(...)
}

fit.gene.mix <- function(...) {
    .Deprecated("FitGeneMix ", package="Seurat")
    FitGeneMix(...)
}

addSmoothedScore <- function(...) {
    .Deprecated("AddSmoothedScore ", package="Seurat")
    AddSmoothedScore(...)
}

addImputedScore <- function(...) {
    .Deprecated("AddImputedScore ", package="Seurat")
    AddImputedScore(...)
}

getNewScore <- function(...) {
  .Deprecated(
    package = 'Seurat',
    msg = 'getNewScore has been removed without replacement'
  )
}

calcNoiseModels <- function(...) {
  .Deprecated(
    package = 'Seurat',
    msg = 'calcNoiseModels has been removed and may be replaced at a later date'
  )
    # .Deprecated("CalcNoiseModels ", package="Seurat")
    # CalcNoiseModels(...)
}

feature.plot.keynote <- function(...) {
  .Deprecated(
    package = 'Seurat',
    msg = 'feature.plot.keynote has been removed without replacement'
  )
}

feature.heatmap <- function(...) {
    .Deprecated("FeatureHeatmap ", package="Seurat")
    FeatureHeatmap(...)
}

ica.plot <- function(...) {
    .Deprecated("ICAPlot ", package="Seurat")
    ICAPlot(...)
}

spatial.de <- function(...) {
  .Deprecated(
    package = 'Seurat',
    msg = 'spatial.de has been removed and may be replaced at a later date'
  )
    # .Deprecated("SpatialDe ", package="Seurat")
    # SpatialDe(...)
}

DBclust_dimension <- function(...) {
    .Deprecated("DBClustDimension ", package="Seurat")
    DBClustDimension(...)
}

Kclust_dimension <- function(...) {
    .Deprecated("KClustDimension ", package="Seurat")
    KClustDimension(...)
}

pca.sig.genes <- function(...) {
    .Deprecated("PCASigGenes ", package="Seurat")
    PCASigGenes(...)
}

doHeatMap <- function(...) {
    .Deprecated("DoHeatmap ", package="Seurat")
    DoHeatmap(...)
}

icHeatmap <- function(...) {
    .Deprecated("ICHeatmap ", package="Seurat")
    ICHeatmap(...)
}

doKMeans <- function(...) {
    .Deprecated("DoKMeans ", package="Seurat")
    DoKMeans(...)
}

genes.in.cluster <- function(...) {
    .Deprecated("GenesInCluster ", package="Seurat")
    GenesInCluster(...)
}

kMeansHeatmap <- function(...) {
    .Deprecated("KMeansHeatmap ", package="Seurat")
    KMeansHeatmap(...)
}

cell.cor.matrix <- function(...) {
  .Deprecated(
    package = 'Seurat',
    msg = 'cell.cor.matrix has been removed and may be replaced at a later date'
  )
    # .Deprecated("CellCorMatrix ", package="Seurat")
    # CellCorMatrix(...)
}

gene.cor.matrix <- function(...) {
  .Deprecated(
    package = 'Seurat',
    msg = 'gene.cor.matrix has been removed and may be replaced at a later date'
  )
    # .Deprecated("GeneCorMatrix ", package="Seurat")
    # GeneCorMatrix(...)
}

calinskiPlot <- function(...) {
  .Deprecated(
    package = 'Seurat',
    msg = 'calinksiPlot has been removed and may be replaced at a later date'
  )
    # .Deprecated("CalinskiPlot ", package="Seurat")
    # CalinskiPlot(...)
}

dot.plot <- function(...) {
    .Deprecated("DotPlot ", package="Seurat")
    DotPlot(...)
}

addMetaData <- function(...) {
    .Deprecated("AddMetaData ", package="Seurat")
    AddMetaData(...)
}

removePC <- function(...) {
  .Deprecated(
    package = 'Seurat',
    msg = 'removePC has been removed and may be replaced at a later date'
  )
    # .Deprecated("RemovePC ", package="Seurat")
    # RemovePC(...)
}

geneScorePlot <- function(...) {
    .Deprecated(
      package = "Seurat",
      msg = 'geneScorePlot has been removed without replacement'
    )
}

cellPlot <- function(...) {
    .Deprecated("CellPlot ", package="Seurat")
    CellPlot(...)
}

jackStraw.permutation.test <- function(...) {
  .Deprecated(
    package="Seurat",
    msg = 'jackStraw.permutation.test has been removed without replacement'
  )
}

jackStrawMC <- function(...) {
  .Deprecated(
    package = 'Seurat',
    msg = 'jackStrawMC has been removed without replacement'
  )
}

jackStrawFull <- function(...) {
  .Deprecated(
    package = 'Seurat',
    msg = 'jackStrawFull has been removed without replacement'
  )
}

PCAFast <- function(...) {
  .Deprecated("PCA", package= "Seurat")
  PCA(...)
}

writ.table <- function(...) {
  .Deprecated(
    new = 'write.table',
    package = 'Seurat',
    msg = "'writ.table' no longer exists, use 'write.table' instead"
  )
}

jackRandom <- function(...) {
  .Deprecated(
    package = 'Seurat',
    msg = "jackRandom has bee removed; it may be added back in the future"
  )
  # .Deprecated(new = 'JackRandom', package = 'Seurat')
  # JackRandom(...)
}

MeanVarPlot <- function(...) {
  .Deprecated(new = 'FindVariableGenes', package = 'Seurat')
  FindVariableGenes(...)
}

HeatmapNode <- function(...) {
  .Deprecated(new = 'NodeHeatmap', package = 'Seurat')
  NodeHeatmap(...)
}

myPalette <- function(...) {
  .Deprecated(new = 'CustomPalette', package = 'Seurat')
  CustomPalette(...)
}

minusr <- function(...) {
  .Deprecated(
    new = 'SubsetRow',
    package = 'Seurat',
    msg = "Use SubsetRow with 'invert = TRUE' instead"
  )
  SubsetRow(..., invert = TRUE)
}

minusc <- function(...) {
  .Deprecated(
    new = 'SubsetColumn',
    package = 'Seurat',
    msg = "Use SubsetColumn with 'invert = TRUE' instead"
  )
  SubsetColumn(..., invert = TRUE)
}

RegressOut <- function(...) {
  .Deprecated(
    new = "ScaleData",
    package = "Seurat",
    msg = "RegressOut functionality has been incorporated into ScaleData. See ?ScaleData for usage details."
  )
  stop()
}

PCAFast <- function(...) {
  .Deprecated("PCA", package= "Seurat")
  PCA(...)
}

writ.table <- function(...) {
  .Deprecated(
    new = 'write.table',
    package = 'Seurat',
    msg = "'writ.table' no longer exists, use 'write.table' instead"
  )
}

jackRandom <- function(...) {
  .Deprecated(new = 'JackRandom', package = 'Seurat')
  JackRandom(...)
}

MeanVarPlot <- function(...) {
  .Deprecated(new = 'FindVariableGenes', package = 'Seurat')
  FindVariableGenes(...)
}

HeatmapNode <- function(...) {
  .Deprecated(new = 'NodeHeatmap', package = 'Seurat')
  NodeHeatmap(...)
}

myPalette <- function(...) {
  .Deprecated(new = 'CustomPalette', package = 'Seurat')
  CustomPalette(...)
}

minusr <- function(...) {
  .Deprecated(
    new = 'SubsetRow',
    package = 'Seurat',
    msg = "Use SubsetRow with 'invert = TRUE' instead"
  )
  SubsetRow(..., invert = TRUE)
}

minusc <- function(...) {
  .Deprecated(
    new = 'SubsetColumn',
    package = 'Seurat',
    msg = "Use SubsetColumn with 'invert = TRUE' instead"
  )
  SubsetColumn(..., invert = TRUE)
}

VizClassification <- function(...) {
  .Deprecated(
    package = 'Seurat',
    msg = 'VizClassification has been removed without replacement'
  )
}
