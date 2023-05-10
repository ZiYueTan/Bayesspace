

#preprocess sce 
preprocess <- function(sce, platform = c("Visium", "ST"), n.PCs = 15, n.HVGs = 2000, skip.PCA = FALSE, log.normalize = TRUE, assay.log = "logcounts"){
  
  # Set BayesSpace.data metadata
  metadata(sce)$BayesSpace.data <- list()
  metadata(sce)$BayesSpace.data$platform <- match.arg(platform)
  metadata(sce)$BayesSpace.data$is.enhanced <- FALSE
  
  # Run PCA on top HVGs
  if (!skip.PCA){
    if(log.normalize)
      sce <- logNormCounts(sce)
    stats <- scran::modelGeneVar(sce, assay.type = assay.log)
    top <- scran::getTopHVGs(sce, n = n.HVGs)
    sce <- scater::runPCA(sce, subset_row = top, ncomponent = n.PCs)
  }
  
  return(sce)
}

library(RCurl)
library(BiocFileCache)
library(assertthat)

#cache If true, cache the dataset locally with BiocFileCache. Otherwise, download directly from our S3 bucket. 

#' | Dataset  | Sample(s) |
#' | ------------- | ------------- |
#' | 2018_thrane_melanoma | ST_mel1_rep2 |
#' | 2020_maynard_prefrontal-cortex  | 151507, 151508, 151509, 151510, 151669, 151670, 151671, 151672, 151673, 151674, 151675, 151676  |
#' | 2020_ji_squamous-cell-carcinoma | P4_rep1 |
#' | 2020_10X-IDC | IDC1 |
#' | 2020_10X-demo_ovarian-cancer | whole_transcriptome |

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

