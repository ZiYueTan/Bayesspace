% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/prepare.R
\name{getRDS}
\alias{getRDS}
\title{Download a processed sample}
\usage{
getRDS(dataset, sample, cache = TRUE)
}
\arguments{
\item{dataset}{identify Dataset}

\item{sample}{identify Sample}

\item{cache}{e If true, cache the dataset locally with \code{BiocFileCache}.}
}
\value{
sce A specify SingleCellExpreriment
}
\description{
Datasets are cached locally using \code{BiocFileCache}. The first time using
this function, you may need to consent to creating a BiocFileCache directory
if one does not already exist.
}
\details{
#' | Dataset  | Sample(s) |
| ------------- | ------------- |
| 2018_thrane_melanoma | ST_mel1_rep2 |
| 2020_maynard_prefrontal-cortex  | 151507, 151508, 151509, 151510, 151669, 151670, 151671, 151672, 151673, 151674, 151675, 151676  |
| 2020_ji_squamous-cell-carcinoma | P4_rep1 |
| 2020_10X-IDC | IDC1 |
| 2020_10X-demo_ovarian-cancer | whole_transcriptome |
}
