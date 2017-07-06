#' Deprecated function(s) in the Seurat package
#'
#' These functions are provided for compatibility with older version of the Seurat package.  They may eventually be completely removed.
#' @rdname Seurat-deprecated
#' @name Seurat-deprecated
#' @param ... Parameters to be passed to the modern version of the function
#' @export setup vlnPlot subsetData mean.var.plot pca project.pca print.pca viz.pca set.ident pca.plot pcHeatmap jackStraw jackStrawPlot run_tsne tsne.plot find.markers find_all_markers genePlot feature.plot tsne.plot buildClusterTree plotClusterTree plotNoiseModel add_samples subsetCells project.samples run_diffusion ica cluster.alpha reorder.ident average.pca average.expression icTopGenes pcTopGenes pcTopCells fetch.data viz.ica getWeightMatrix regulatorScore find.markers.node diffExp.test tobit.test batch.gene marker.test diff.t.test which.cells set.all.ident rename.ident posterior.plot map.cell get.centroids refined.mapping initial.mapping calc.insitu fit.gene.k fit.gene.mix addSmoothedScore addImputedScore getNewScore calcNoiseModels feature.plot.keynote feature.heatmap ica.plot dim.plot spatial.de DBclust_dimension Kclust_dimension pca.sig.genes doHeatMap icHeatmap doKMeans genes.in.cluster kMeansHeatmap cell.cor.matrix gene.cor.matrix calinskiPlot dot.plot addMetaData removePC geneScorePlot cellPlot jackStraw.permutation.test jackStrawMC jackStrawFull writ.table jackRandom
#' @aliases setup vlnPlot subsetData mean.var.plot pca project.pca print.pca viz.pca set.ident pca.plot pcHeatmap jackStraw jackStrawPlot run_tsne tsne.plot find.markers find_all_markers genePlot feature.plot tnse.plot buildClusterTree plotClusterTree plotNoiseModel add_samples subsetCells project.samples run_diffusion ica cluster.alpha reorder.ident average.pca average.expression icTopGenes pcTopGenes pcTopCells fetch.data viz.ica getWeightMatrix regulatorScore find.markers.node diffExp.test tobit.test batch.gene marker.test diff.t.test which.cells set.all.ident rename.ident posterior.plot map.cell get.centroids refined.mapping initial.mapping calc.insitu fit.gene.k fit.gene.mix addSmoothedScore addImputedScore getNewScore calcNoiseModels feature.plot.keynote feature.heatmap ica.plot dim.plot spatial.de DBclust_dimension Kclust_dimension pca.sig.genes doHeatMap icHeatmap doKMeans genes.in.cluster kMeansHeatmap cell.cor.matrix gene.cor.matrix calinskiPlot dot.plot addMetaData removePC geneScorePlot cellPlot jackStraw.permutation.test jackStrawMC jackStrawFull writ.table jackRandom
#' @section Details:
#' \tabular{rl}{
#'   \code{setup} \tab now a synonym for \code{\link{Setup}}\cr
#'   \code{vlnPlot} \tab now a synonym for \code{\link{VlnPlot}}\cr
#'   \code{subsetData} \tab now a synonym for \code{\link{SubsetData}}\cr
#'   \code{mean.var.plot} \tab now a synonym for \code{\link{MeanVarPlot}}\cr
#'   \code{pca} \tab now a synonym for \code{\link{PCA}}\cr
#'   \code{project.pca} \tab now a synonym for \code{\link{ProjectPCA}}\cr
#'   \code{print.pca} \tab now a synonym for \code{\link{PrintPCA}}\cr
#'   \code{viz.pca} \tab now a synonym for \code{\link{VizPCA}}\cr
#'   \code{set.ident} \tab now a synonym for \code{\link{SetIdent}}\cr
#'   \code{pca.plot} \tab now a synonym for \code{\link{PCAPlot}}\cr
#'   \code{pcHeatmap} \tab now a synonym for \code{\link{PCHeatmap}}\cr
#'   \code{jackStraw} \tab now a synonym for \code{\link{JackStraw}}\cr
#'   \code{jackStrawPlot} \tab now a synonym for \code{\link{JackStrawPlot}}\cr
#'   \code{run_tsne} \tab now a synonym for \code{\link{RunTSNE}}\cr
#'   \code{tsne.plot} \tab now a synonym for \code{\link{TSNEPlot}}\cr
#'   \code{find.markers} \tab now a synonym for \code{\link{FindMarkers}}\cr
#'   \code{find_all_markers} \tab now a synonym for \code{\link{FindAllMarkers}}\cr
#'   \code{genePlot} \tab now a synonym for \code{\link{GenePlot}}\cr
#'   \code{feature.plot} \tab now a synonym for \code{\link{FeaturePlot}}\cr
#'   \code{buildClusterTree} \tab now a synonym for \code{\link{BuildClusterTree}}\cr
#'   \code{plotClusterTree} \tab now a synonym for \code{\link{PlotClusterTree}}\cr
#'   \code{plotNoiseModel} \tab now a synonym for \code{\link{PlotNoiseModel}}\cr
#'   \code{add_samples} \tab now a synonym for \code{\link{AddSamples}}\cr
#'   \code{subsetCells} \tab now a synonym for \code{\link{SubsetCells}}\cr
#'   \code{project.samples} \tab now a synonym for \code{\link{ProjectSamples}}\cr
#'   \code{run_diffusion} \tab now a synonym for \code{\link{RunDiffusion}}\cr
#'   \code{ica} \tab now a synonym for \code{\link{ICA}}\cr
#'   \code{cluster.alpha} \tab now a synonym for \code{\link{ClusterAlpha}}\cr
#'   \code{reorder.ident} \tab now a synonym for \code{\link{ReorderIdent}}\cr
#'   \code{average.pca} \tab now a synonym for \code{\link{AveragePCA}}\cr
#'   \code{average.expression} \tab now a synonym for \code{\link{AverageExpression}}\cr
#'   \code{icTopGenes} \tab now a synonym for \code{\link{ICTopGenes}}\cr
#'   \code{pcTopGenes} \tab now a synonym for \code{\link{PCTopGenes}}\cr
#'   \code{pcTopCells} \tab now a synonym for \code{\link{PCTopCells}}\cr
#'   \code{fetch.data} \tab now a synonym for \code{\link{FetchData}}\cr
#'   \code{viz.ica} \tab now a synonym for \code{\link{VizIca}}\cr
#'   \code{getWeightMatrix} \tab now a synonym for \code{\link{GetWeightMatrix}}\cr
#'   \code{regulatorScore} \tab now a synonym for \code{\link{RegulatorScore}}\cr
#'   \code{find.markers.node} \tab now a synonym for \code{\link{FindMarkersNode}}\cr
#'   \code{diffExp.test} \tab now a synonym for \code{\link{DiffExpTest}}\cr
#'   \code{tobit.test} \tab now a synonym for \code{\link{TobitTest}}\cr
#'   \code{batch.gene} \tab now a synonym for \code{\link{BatchGene}}\cr
#'   \code{marker.test} \tab now a synonym for \code{\link{MarkerTest}}\cr
#'   \code{diff.t.test} \tab now a synonym for \code{\link{DiffTTest}}\cr
#'   \code{which.cells} \tab now a synonym for \code{\link{WhichCells}}\cr
#'   \code{set.all.ident} \tab now a synonym for \code{\link{SetAllIdent}}\cr
#'   \code{rename.ident} \tab now a synonym for \code{\link{RenameIdent}}\cr
#'   \code{posterior.plot} \tab now a synonym for \code{\link{PosteriorPlot}}\cr
#'   \code{map.cell} \tab now a synonym for \code{\link{MapCell}}\cr
#'   \code{get.centroids} \tab now a synonym for \code{\link{GetCentroids}}\cr
#'   \code{refined.mapping} \tab now a synonym for \code{\link{RefinedMapping}}\cr
#'   \code{initial.mapping} \tab now a synonym for \code{\link{InitialMapping}}\cr
#'   \code{calc.insitu} \tab now a synonym for \code{\link{CalcInsitu}}\cr
#'   \code{fit.gene.k} \tab now a synonym for \code{\link{FitGeneK}}\cr
#'   \code{fit.gene.mix} \tab now a synonym for \code{\link{FitGeneMix}}\cr
#'   \code{addSmoothedScore} \tab now a synonym for \code{\link{AddSmoothedScore}}\cr
#'   \code{addImputedScore} \tab now a synonym for \code{\link{AddImputedScore}}\cr
#'   \code{getNewScore} \tab now a synonym for \code{\link{GetNewScore}}\cr
#'   \code{calcNoiseModels} \tab now a synonym for \code{\link{CalcNoiseModels}}\cr
#'   \code{feature.plot.keynote} \tab now a synonym for \code{\link{FeaturePlotKeynote}}\cr
#'   \code{feature.heatmap} \tab now a synonym for \code{\link{FeatureHeatmap}}\cr
#'   \code{ica.plot} \tab now a synonym for \code{\link{ICAPlot}}\cr
#'   \code{dim.plot} \tab now a synonym for \code{\link{DimPlot}}\cr
#'   \code{spatial.de} \tab now a synonym for \code{\link{SpatialDe}}\cr
#'   \code{DBclust_dimension} \tab now a synonym for \code{\link{DBClustDimension}}\cr
#'   \code{Kclust_dimension} \tab now a synonym for \code{\link{KClustDimension}}\cr
#'   \code{pca.sig.genes} \tab now a synonym for \code{\link{PCASigGenes}}\cr
#'   \code{doHeatMap} \tab now a synonym for \code{\link{DoHeatMap}}\cr
#'   \code{icHeatmap} \tab now a synonym for \code{\link{ICHeatmap}}\cr
#'   \code{doKMeans} \tab now a synonym for \code{\link{DoKMeans}}\cr
#'   \code{genes.in.cluster} \tab now a synonym for \code{\link{GenesInCluster}}\cr
#'   \code{kMeansHeatmap} \tab now a synonym for \code{\link{KMeansHeatmap}}\cr
#'   \code{cell.cor.matrix} \tab now a synonym for \code{\link{CellCorMatrix}}\cr
#'   \code{gene.cor.matrix} \tab now a synonym for \code{\link{GeneCorMatrix}}\cr
#'   \code{calinskiPlot} \tab now a synonym for \code{\link{CalinskiPlot}}\cr
#'   \code{dot.plot} \tab now a synonym for \code{\link{DotPlot}}\cr
#'   \code{addMetaData} \tab now a synonym for \code{\link{AddMetaData}}\cr
#'   \code{removePC} \tab now a synonym for \code{\link{RemovePC}}\cr
#'   \code{geneScorePlot} \tab now a synonym for \code{\link{GeneScorePlot}}\cr
#'   \code{cellPlot} \tab now a synonym for \code{\link{CellPlot}}\cr
#'   \code{jackStraw.permutation.test} \tab now a synonym for \code{\link{JackStrawPermutationTest}}\cr
#'   \code{jackStrawMC} \tab now a synonym for \code{\link{JackStrawMC}}\cr
#'   \code{jackStrawFull} \tab now a synonym for \code{\link{JackStrawFull}}\cr
#'   \code{writ.table} \tab is delteded without replacement\cr
#'   \code{jackRandom} \tab now a synonym for \code{\link{JackRandom}}\cr
#' }
#'

setup <- function(...) {
    .Deprecated("Setup", package="Seurat")
    Setup(...)
}

vlnPlot <- function(...) {
    .Deprecated("VlnPlot", package="Seurat")
    VlnPlot(...)
}

subsetData <- function(...) {
    .Deprecated("SubsetData", package="Seurat")
    SubsetData(...)
}

mean.var.plot <- function(...) {
    .Deprecated("MeanVarPlot", package="Seurat")
    MeanVarPlot(...)
}

pca <- function(...) {
    .Deprecated("PCA", package="Seurat")
    PCA(...)
}

project.pca <- function(...) {
    .Deprecated("ProjectPCA", package="Seurat")
    ProjectPCA(...)
}

print.pca <- function(...) {
    .Deprecated("PrintPCA", package="Seurat")
    PrintPCA(...)
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
    .Deprecated("PlotNoiseModel", package="Seurat")
    PlotNoiseModel(...)
}

add_samples <- function(...) {
    .Deprecated("AddSamples", package="Seurat")
    AddSamples(...)
}

subsetCells <- function(...) {
    .Deprecated("SubsetCells", package="Seurat")
    SubsetCells(...)
}

project.samples <- function(...) {
    .Deprecated("ProjectSamples", package="Seurat")
    ProjectSamples(...)
}

run_diffusion <- function(...) {
    .Deprecated("RunDiffusion", package="Seurat")
    RunDiffusion(...)
}

ica <- function(...) {
    .Deprecated("ICA", package="Seurat")
    ICA(...)
}

cluster.alpha <- function(...) {
    .Deprecated("ClusterAlpha", package="Seurat")
    ClusterAlpha(...)
}

reorder.ident <- function(...) {
    .Deprecated("ReorderIdent", package="Seurat")
    ReorderIdent(...)
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

getWeightMatrix <- function(...) {
    .Deprecated("GetWeightMatrix", package="Seurat")
    GetWeightMatrix(...)
}

regulatorScore <- function(...) {
    .Deprecated("RegulatorScore", package="Seurat")
    RegulatorScore(...)
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
    .Deprecated("BatchGene", package="Seurat")
    BatchGene(...)
}

marker.test <- function(...) {
    .Deprecated("MarkerTest ", package="Seurat")
    MarkerTest(...)
}

diff.t.test <- function(...) {
    .Deprecated("DiffTTest", package="Seurat")
    DiffTTest(...)
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
    .Deprecated("MapCell", package="Seurat")
    MapCell(...)
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
    .Deprecated("GetNewScore ", package="Seurat")
    GetNewScore(...)
}

calcNoiseModels <- function(...) {
    .Deprecated("CalcNoiseModels ", package="Seurat")
    CalcNoiseModels(...)
}

feature.plot.keynote <- function(...) {
    .Deprecated("FeaturePlotKeynote ", package="Seurat")
    FeaturePlotKeynote(...)
}

feature.heatmap <- function(...) {
    .Deprecated("FeatureHeatmap ", package="Seurat")
    FeatureHeatmap(...)
}

ica.plot <- function(...) {
    .Deprecated("ICAPlot ", package="Seurat")
    ICAPlot(...)
}

dim.plot <- function(...) {
    .Deprecated("DimPlot ", package="Seurat")
    DimPlot(...)
}

spatial.de <- function(...) {
    .Deprecated("SpatialDe ", package="Seurat")
    SpatialDe(...)
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
    .Deprecated("CellCorMatrix ", package="Seurat")
    CellCorMatrix(...)
}

gene.cor.matrix <- function(...) {
    .Deprecated("GeneCorMatrix ", package="Seurat")
    GeneCorMatrix(...)
}

calinskiPlot <- function(...) {
    .Deprecated("CalinskiPlot ", package="Seurat")
    CalinskiPlot(...)
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
    .Deprecated("RemovePC ", package="Seurat")
    RemovePC(...)
}

geneScorePlot <- function(...) {
    .Deprecated("GeneScorePlot ", package="Seurat")
    GeneScorePlot(...)
}

cellPlot <- function(...) {
    .Deprecated("CellPlot ", package="Seurat")
    CellPlot(...)
}

jackStraw.permutation.test <- function(...) {
    .Deprecated("JackStrawPermutationTest ", package="Seurat")
    JackStrawPermutationTest(...)
}

jackStrawMC <- function(...) {
    .Deprecated("JackStrawMC ", package="Seurat")
    JackStrawMC(...)
}

jackStrawFull <- function(...) {
    .Deprecated("JackStrawFull ", package="Seurat")
    JackStrawFull(...)
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
