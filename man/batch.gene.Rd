\name{batch.gene}
\alias{batch.gene}
\title{
Identify potential genes associated with batch effects
}
\description{
Test for genes whose expression value is strongly predictive of batch (based on ROC classification). Important note: Assumes that the 'batch' of each cell is assigned to the cell's identity class (will be improved in a future release)
}
\usage{
batch.gene(object, idents.use, genes.use = NULL, use.imputed = FALSE, auc.cutoff = 0.6, thresh.use = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
Seurat object
}
  \item{idents.use}{
Batch names to test
}
  \item{genes.use}{
Gene list to test
}
  \item{auc.cutoff}{
Minimum AUC needed to qualify as a 'batch gene'
}
  \item{thresh.use}{
  Limit testing to genes which show, on average, at least X-fold difference (log-scale) in any one batch}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Returns a list of genes that are strongly correlated with batch.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}
