% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/enhanced_cluster.R
\name{enhanced.results}
\alias{enhanced.results}
\title{Get the subspots cluster assignments of the BayesSpace}
\usage{
enhanced.results(
  Y,
  q,
  sub.neighbors,
  platform,
  init,
  mu0 = colMeans(Y),
  Lambda0 = diag(ncol(Y)) * 0.01,
  gamma = 3,
  alpha = 1,
  beta = 0.01,
  jitter_scale,
  jitter_prior,
  ncores = 2,
  n.rep = 1000,
  save_para = FALSE
)
}
\arguments{
\item{Y}{PCs of original sce}

\item{q}{The number of clusters}

\item{sub.neighbors}{A list of neighbor indices for each subspot}

\item{platform}{Spatial transcriptomic platform.}

\item{init}{Initial cluster assignments for spots.}

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
proportion  of the mean variance of the PCs.}

\item{ncores}{The number of cores used for parallel calculation.}

\item{n.rep}{The number of MCMC iterations.}

\item{save_para}{If true, save the list of clustering parameters.}
}
\value{
If \code{save_para = FALSE}, return the vector of subspots cluster assignments.
If \code{save_para = TRUE}, return the list of clustering parameter values
 apart from the vector of subspots cluster assignments.
}
\description{
Get the subspots cluster assignments of the BayesSpace
}
\keyword{internal}
