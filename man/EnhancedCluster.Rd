% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/enhanced_cluster.R
\name{EnhancedCluster}
\alias{EnhancedCluster}
\title{Enhance spot resolution}
\usage{
EnhancedCluster(
  sce,
  q,
  d,
  platform = c("Visium", "ST"),
  init = NULL,
  init.method = c("spatialcluster", "mclust", "kmeans"),
  mu0 = NULL,
  Lambda0 = NULL,
  gamma = 3,
  alpha = 1,
  beta = 0.01,
  jitter_scale = 5,
  jitter_prior = 0.3,
  n.rep = 1000,
  ncores = 2,
  save_para = FALSE
)
}
\arguments{
\item{sce}{A SingleCellExperiment object containing the spatial data.}

\item{q}{The number of clusters.}

\item{d}{Number of top principal components to use when clustering.}

\item{platform}{Number of top principal components to use when clustering.}

\item{init}{Initial cluster assignments for spots.}

\item{init.method}{If \code{init} is not provided, cluster the top \code{d}
PCs with this method to obtain initial cluster assignments.}

\item{mu0}{Prior mean hyperparameter for mu. If not provided, mu0 is set to
the mean of PCs over all spots.}

\item{Lambda0}{Prior precision hyperparam for mu. If not provided, lambda0
is set to a diagonal matrix \eqn{0.01 I}.}

\item{gamma}{Smoothing parameter.}

\item{alpha}{Hyperparameter for Wishart distributed precision lambda.}

\item{beta}{Hyperparameter for Wishart distributed precision lambda.}

\item{jitter_scale}{Controls the amount of jittering. Small amounts of
jittering are more likely to be accepted but result in exploring the space
more slowly. We suggest tuning \code{jitter_scale} so that Ychange is on
average around 25\%-40\%.}

\item{jitter_prior}{Scale factor for the prior variance, parameterized as the
proportion of the mean variance of the PCs.}

\item{n.rep}{The number of MCMC iterations.}

\item{ncores}{The number of cores used for parallel calculation.}

\item{save_para}{If true, save the list of clustering parameters.}
}
\value{
Returns a new SingleCellExperiment object. By default, the
  \code{assays} of this object are empty, and the enhanced resolution PCs
  are stored as a reduced dimensionality result accessible with
  \code{reducedDim(sce, 'PCA')}.
}
\description{
Enhanced clustering of a spatial expression dataset to subspot resolution.
}
\details{
The enhanced \code{SingleCellExperiment} has most of the properties of the
  input SCE - \code{rowData}, \code{colData}, \code{reducedDims} - but does
  not include expression data in \code{counts} or \code{logcounts}. To impute
  enhanced expression vectors, please use [EnhanceFeatures()] after
  running \code{spatialEnhance}.

The \code{colData} of the enhanced \code{SingleCellExperiment} includes the
  following columns to permit referencing the subspots in spatial context and
  linking back to the original spots:
  \itemize{
  \item \code{spot.idx}: Index of the spot this subspot belongs to (with
    respect to the input SCE).
  \item \code{subspot.idx}: Index of the subspot within its parent spot.
  \item \code{spot.row}: Array row of the subspot's parent spot.
  \item \code{spot.col}: Array col of the subspot's parent spot.
  \item \code{row}: Array row of the subspot. This is the parent spot's row
    plus an offset based on the subspot's position within the spot.
  \item \code{col}: Array col of the subspot. This is the parent spot's col
    plus an offset based on the subspot's position within the spot.
    }

(optional)
The BayesSpace of the enhanced \code{SingleCellExperiment} includes
the list of clustering parameters.
}
