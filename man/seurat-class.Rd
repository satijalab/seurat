\name{seurat-class}
\Rdversion{1.1}
\docType{class}
\alias{seurat-class}
\alias{seurat}
\alias{addImputedScore,seurat-method}
\alias{addMetaData,seurat-method}
\alias{average.expression,seurat-method}
\alias{average.pca,seurat-method}
\alias{batch.gene,seurat-method}
\alias{buildClusterTree,seurat-method}
\alias{calc.insitu,seurat-method}
\alias{calcNoiseModels,seurat-method}
\alias{calinskiPlot,seurat-method}
\alias{cell.cor.matrix,seurat-method}
\alias{cellPlot,seurat-method}
\alias{cells.in.cluster,seurat-method}
\alias{cluster.alpha,seurat-method}
\alias{diff.t.test,seurat-method}
\alias{diffExp.test,seurat-method}
\alias{doHeatMap,seurat-method}
\alias{doKMeans,seurat-method}
\alias{dot.plot,seurat-method}
\alias{feature.heatmap,seurat-method}
\alias{feature.plot,seurat-method}
\alias{fetch.data,seurat-method}
\alias{find_all_markers_node,seurat-method}
\alias{find_all_markers,seurat-method}
\alias{find.markers,seurat-method}
\alias{find.markers.node,seurat-method}
\alias{fit.gene.k,seurat-method}
\alias{fit.gene.mix,seurat-method}
\alias{gene.cor.matrix,seurat-method}
\alias{genePlot,seurat-method}
\alias{genes.in.cluster,seurat-method}
\alias{geneScorePlot,seurat-method}
\alias{getNewScore,seurat-method}
\alias{getWeightMatrix,seurat-method}
\alias{ica,seurat-method}
\alias{initial.mapping,seurat-method}
\alias{jackStraw,seurat-method}
\alias{jackStrawFull,seurat-method}
\alias{jackStrawPlot,seurat-method}
\alias{jackStrawPlot2,seurat-method}
\alias{loadMetrics,seurat-method}
\alias{map.cell,seurat-method}
\alias{marker.test,seurat-method}
\alias{Mclust_dimension,seurat-method}
\alias{mean.var.plot,seurat-method}
\alias{pca,seurat-method}
\alias{pca.plot,seurat-method}
\alias{pca.sig.genes,seurat-method}
\alias{plotClusterTree,seurat-method}
\alias{plotNoiseModel,seurat-method}
\alias{posterior.plot,seurat-method}
\alias{print.pca,seurat-method}
\alias{project.pca,seurat-method}
\alias{refined.mapping,seurat-method}
\alias{regulatorScore,seurat-method}
\alias{reorder.by.tree,seurat-method}
\alias{retreiveCellInfo,seurat-method}
\alias{retreiveCluster,seurat-method}
\alias{retreiveScore,seurat-method}
\alias{run_tsne,seurat-method}
\alias{set.all.stat,seurat-method}
\alias{set.stat,seurat-method}
\alias{setup,seurat-method}
\alias{show,seurat-method}
\alias{spatial.de,seurat-method}
\alias{subsetData,seurat-method}
\alias{tobit.test,seurat-method}
\alias{viz.pca,seurat-method}
\alias{vlnPlot,seurat-method}
\alias{which.cells,seurat-method}

\title{Class \code{"seurat"}}
\description{
%%  ~~ A concise (1-5 lines) description of what the class is. ~~
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{seurat(...)}.
%%  ~~ describe objects here ~~ 
}
\section{Slots}{
  \describe{
    \item{\code{raw.data}:}{Object of class \code{"data.frame"} ~~ }
    \item{\code{data}:}{Object of class \code{"data.frame"} ~~ }
    \item{\code{scale.data}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{var.genes}:}{Object of class \code{"vector"} ~~ }
    \item{\code{is.expr}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{stat.fxn}:}{Object of class \code{"function"} ~~ }
    \item{\code{ident}:}{Object of class \code{"vector"} ~~ }
    \item{\code{data.ngene}:}{Object of class \code{"vector"} ~~ }
    \item{\code{pca.x}:}{Object of class \code{"data.frame"} ~~ }
    \item{\code{pca.rot}:}{Object of class \code{"data.frame"} ~~ }
    \item{\code{real.fval}:}{Object of class \code{"data.frame"} ~~ }
    \item{\code{fake.fval}:}{Object of class \code{"data.frame"} ~~ }
    \item{\code{emp.pval}:}{Object of class \code{"data.frame"} ~~ }
    \item{\code{kmeans.obj}:}{Object of class \code{"list"} ~~ }
    \item{\code{pca.obj}:}{Object of class \code{"list"} ~~ }
    \item{\code{calinski.best}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{gene.scores}:}{Object of class \code{"data.frame"} ~~ }
    \item{\code{k.num}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{drop.coefs}:}{Object of class \code{"data.frame"} ~~ }
    \item{\code{wt.matrix}:}{Object of class \code{"data.frame"} ~~ }
    \item{\code{drop.wt.matrix}:}{Object of class \code{"data.frame"} ~~ }
    \item{\code{trusted.genes}:}{Object of class \code{"vector"} ~~ }
    \item{\code{drop.expr}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{data.info}:}{Object of class \code{"data.frame"} ~~ }
    \item{\code{project.name}:}{Object of class \code{"character"} ~~ }
    \item{\code{project.dir}:}{Object of class \code{"character"} ~~ }
    \item{\code{kmeans.gene}:}{Object of class \code{"list"} ~~ }
    \item{\code{kmeans.cell}:}{Object of class \code{"list"} ~~ }
    \item{\code{jackStraw.empP}:}{Object of class \code{"data.frame"} ~~ }
    \item{\code{pc.x.full}:}{Object of class \code{"data.frame"} ~~ }
    \item{\code{jackStraw.fakePC}:}{Object of class \code{"data.frame"} ~~ }
    \item{\code{jackStraw.empP.full}:}{Object of class \code{"data.frame"} ~~ }
    \item{\code{pca.x.full}:}{Object of class \code{"data.frame"} ~~ }
    \item{\code{kmeans.col}:}{Object of class \code{"list"} ~~ }
    \item{\code{mean.var}:}{Object of class \code{"data.frame"} ~~ }
    \item{\code{imputed}:}{Object of class \code{"data.frame"} ~~ }
    \item{\code{mix.probs}:}{Object of class \code{"data.frame"} ~~ }
    \item{\code{mix.mu}:}{Object of class \code{"data.frame"} ~~ }
    \item{\code{mix.sigma}:}{Object of class \code{"data.frame"} ~~ }
    \item{\code{mu.alpha}:}{Object of class \code{"data.frame"} ~~ }
    \item{\code{mix.param}:}{Object of class \code{"data.frame"} ~~ }
    \item{\code{final.prob}:}{Object of class \code{"data.frame"} ~~ }
    \item{\code{insitu.matrix}:}{Object of class \code{"data.frame"} ~~ }
    \item{\code{tsne.rot}:}{Object of class \code{"data.frame"} ~~ }
    \item{\code{ica.rot}:}{Object of class \code{"data.frame"} ~~ }
    \item{\code{ica.x}:}{Object of class \code{"data.frame"} ~~ }
    \item{\code{ica.obj}:}{Object of class \code{"list"} ~~ }
    \item{\code{cell.names}:}{Object of class \code{"vector"} ~~ }
    \item{\code{cluster.tree}:}{Object of class \code{"list"} ~~ }
  }
}
\section{Methods}{
  \describe{
    \item{addImputedScore}{\code{signature(object = "seurat")}: ... }
    \item{addMetaData}{\code{signature(object = "seurat")}: ... }
    \item{average.expression}{\code{signature(object = "seurat")}: ... }
    \item{average.pca}{\code{signature(object = "seurat")}: ... }
    \item{batch.gene}{\code{signature(object = "seurat")}: ... }
    \item{buildClusterTree}{\code{signature(object = "seurat")}: ... }
    \item{calc.insitu}{\code{signature(object = "seurat")}: ... }
    \item{calcNoiseModels}{\code{signature(object = "seurat")}: ... }
    \item{calinskiPlot}{\code{signature(object = "seurat")}: ... }
    \item{cell.cor.matrix}{\code{signature(object = "seurat")}: ... }
    \item{cellPlot}{\code{signature(object = "seurat")}: ... }
    \item{cells.in.cluster}{\code{signature(object = "seurat")}: ... }
    \item{cluster.alpha}{\code{signature(object = "seurat")}: ... }
    \item{diff.t.test}{\code{signature(object = "seurat")}: ... }
    \item{diffExp.test}{\code{signature(object = "seurat")}: ... }
    \item{doHeatMap}{\code{signature(object = "seurat")}: ... }
    \item{doKMeans}{\code{signature(object = "seurat")}: ... }
    \item{dot.plot}{\code{signature(object = "seurat")}: ... }
    \item{feature.heatmap}{\code{signature(object = "seurat")}: ... }
    \item{feature.plot}{\code{signature(object = "seurat")}: ... }
    \item{fetch.data}{\code{signature(object = "seurat")}: ... }
    \item{find_all_markers_node}{\code{signature(object = "seurat")}: ... }
    \item{find_all_markers}{\code{signature(object = "seurat")}: ... }
    \item{find.markers}{\code{signature(object = "seurat")}: ... }
    \item{find.markers.node}{\code{signature(object = "seurat")}: ... }
    \item{fit.gene.k}{\code{signature(object = "seurat")}: ... }
    \item{fit.gene.mix}{\code{signature(object = "seurat")}: ... }
    \item{gene.cor.matrix}{\code{signature(object = "seurat")}: ... }
    \item{genePlot}{\code{signature(object = "seurat")}: ... }
    \item{genes.in.cluster}{\code{signature(object = "seurat")}: ... }
    \item{geneScorePlot}{\code{signature(object = "seurat")}: ... }
    \item{getNewScore}{\code{signature(object = "seurat")}: ... }
    \item{getWeightMatrix}{\code{signature(object = "seurat")}: ... }
    \item{ica}{\code{signature(object = "seurat")}: ... }
    \item{initial.mapping}{\code{signature(object = "seurat")}: ... }
    \item{jackStraw}{\code{signature(object = "seurat")}: ... }
    \item{jackStrawFull}{\code{signature(object = "seurat")}: ... }
    \item{jackStrawPlot}{\code{signature(object = "seurat")}: ... }
    \item{jackStrawPlot2}{\code{signature(object = "seurat")}: ... }
    \item{loadMetrics}{\code{signature(object = "seurat")}: ... }
    \item{map.cell}{\code{signature(object = "seurat")}: ... }
    \item{marker.test}{\code{signature(object = "seurat")}: ... }
    \item{Mclust_dimension}{\code{signature(object = "seurat")}: ... }
    \item{mean.var.plot}{\code{signature(object = "seurat")}: ... }
    \item{pca}{\code{signature(object = "seurat")}: ... }
    \item{pca.plot}{\code{signature(object = "seurat")}: ... }
    \item{pca.sig.genes}{\code{signature(object = "seurat")}: ... }
    \item{plotClusterTree}{\code{signature(object = "seurat")}: ... }
    \item{plotNoiseModel}{\code{signature(object = "seurat")}: ... }
    \item{posterior.plot}{\code{signature(object = "seurat")}: ... }
    \item{print.pca}{\code{signature(object = "seurat")}: ... }
    \item{project.pca}{\code{signature(object = "seurat")}: ... }
    \item{refined.mapping}{\code{signature(object = "seurat")}: ... }
    \item{regulatorScore}{\code{signature(object = "seurat")}: ... }
    \item{reorder.by.tree}{\code{signature(object = "seurat")}: ... }
    \item{retreiveCellInfo}{\code{signature(object = "seurat")}: ... }
    \item{retreiveCluster}{\code{signature(object = "seurat")}: ... }
    \item{retreiveScore}{\code{signature(object = "seurat")}: ... }
    \item{run_tsne}{\code{signature(object = "seurat")}: ... }
    \item{set.all.stat}{\code{signature(object = "seurat")}: ... }
    \item{set.stat}{\code{signature(object = "seurat")}: ... }
    \item{setup}{\code{signature(object = "seurat")}: ... }
    \item{show}{\code{signature(object = "seurat")}: ... }
    \item{spatial.de}{\code{signature(object = "seurat")}: ... }
    \item{subsetData}{\code{signature(object = "seurat")}: ... }
    \item{tobit.test}{\code{signature(object = "seurat")}: ... }
    \item{viz.pca}{\code{signature(object = "seurat")}: ... }
    \item{vlnPlot}{\code{signature(object = "seurat")}: ... }
    \item{which.cells}{\code{signature(object = "seurat")}: ... }
	 }
}
\references{
%%  ~~put references to the literature/web site here~~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%%  ~~objects to See Also as \code{\link{~~fun~~}}, ~~~
%%  ~~or \code{\linkS4class{CLASSNAME}} for links to other classes ~~~
}
\examples{
showClass("seurat")
}
\keyword{classes}
\usage{
seurat(...)
}
\arguments{
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
}
}
