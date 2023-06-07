
#' Preprocess SingleCellExperiemnt for BayesSpace
#'
#' Add metadata for further analysis, and optionally performs PCA on log-normalized expression on top HVGs.
#'
#' @param sce  SingleCellExperiment to preprocess
#' @param platform  Spatial sequencing platform.
#' @param n.PCs  Number of principal components to compute.
#' @param n.HVGs Number of highly variable genes to run PCA upon.
#' @param skip.PCA Whether to skip PCA if dimentionality reduction was previously computed.
#' @param log.normalize Whether to log-normalize the input data with scater.
#' @param assay.log  Name of assay in \code{sce} containing normalized counts.
#'
#' @return SingleCellExpriment with PCA and BayesSpace metadata
#'
#' @export
#'
#' @importFrom scater logNormCounts runPCA
#' @importFrom scran modelGeneVar getTopHVGs
#' 
preprocess <- function(sce, platform = c("Visium", "ST"), n.PCs = 15, n.HVGs = 2000, skip.PCA = FALSE, log.normalize = TRUE, assay.log = "logcounts"){

  ## Set BayesSpace.data metadata
  metadata(sce)$BayesSpace.data <- list()
  metadata(sce)$BayesSpace.data$platform <- match.arg(platform)
  metadata(sce)$BayesSpace.data$is.enhanced <- FALSE

  ## Run PCA on top HVGs
  if (!skip.PCA){
    if(log.normalize)
      sce <- scater::logNormCounts(sce)
    stats <- scran::modelGeneVar(sce, assay.type = assay.log)
    top <- scran::getTopHVGs(sce, n = n.HVGs)
    sce <- scater::runPCA(sce, subset_row = top, ncomponent = n.PCs)
  }

  return(sce)
}




#' Download a processed sample
#'
#' Datasets are cached locally using \code{BiocFileCache}. The first time using
#' this function, you may need to consent to creating a BiocFileCache directory
#' if one does not already exist.
#'
#' @param dataset identify Dataset
#' @param sample identify Sample
#' @param cache e If true, cache the dataset locally with \code{BiocFileCache}.
#'
#' @return sce A specify SingleCellExpreriment
#' @details
#'
#' #' | Dataset  | Sample(s) |
#' | ------------- | ------------- |
#' | 2018_thrane_melanoma | ST_mel1_rep2 |
#' | 2020_maynard_prefrontal-cortex  | 151507, 151508, 151509, 151510, 151669, 151670, 151671, 151672, 151673, 151674, 151675, 151676  |
#' | 2020_ji_squamous-cell-carcinoma | P4_rep1 |
#' | 2020_10X-IDC | IDC1 |
#' | 2020_10X-demo_ovarian-cancer | whole_transcriptome |
#' @export
#' @importFrom RCurl url.exists
#' @importFrom utils download.file
#' @importFrom assertthat assert_that
#' @importFrom BiocFileCache BiocFileCache bfcrpath



getRDS <- function(dataset, sample, cache = TRUE) {


  url <- "https://fh-pi-gottardo-r-eco-public.s3.amazonaws.com/SpatialTranscriptomes/%s/%s.rds"
  url <- sprintf(url, dataset, sample)
  assert_that(url.exists(url), msg="Dataset/sample not available")

  if (cache) {
    bfc <- BiocFileCache()
    local.path <- bfcrpath(bfc, url)
  } else {
    local.path <- tempfile(fileext=".rds")
    download.file(url, local.path, quiet=TRUE, mode="wb")
  }

  readRDS(local.path)
}

