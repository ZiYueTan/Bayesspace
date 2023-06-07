



#' Define offsets for each subspot layout.
#'
#' Hex spots are divided into 6 triangular subspots, square spots are divided
#' into 9 squares. Offsets are relative to the spot center.
#'
#' @param platform If "Visium", devide hex spots into 6 triangular subspots ; if
#'   "ST", devide square spots into 9 squares.
#'
#' @return Matrix of x and y offsets, one row per subspot
#'
#' @keywords internal
subplot.offsets <- function(platform = c("Visium", "ST")){
  if(platform == "Visium"){
    offsets <- rbind(expand.grid(c(1/3, -1/3), c(1/3, -1/3)), expand.grid(c(2/3, -2/3), 0))
  }
  if(platform == "ST"){
    offsets <- rbind(expand.grid(c(1/3, -1/3, 0), c(1/3, -1/3, 0)))
  }
  colnames(offsets) <- c("offsets.row", "offsets.col")
  return(offsets)
}


#' Add subspot labels and offset row/col locations.
#'
#' Subspots are stored as (1.1, 2.1, 3.1, ..., 1.2, 2.2, 3.2, ...)
#'
#' @param sce  Original sce
#' @param platform If "Visium", make 6 suplots; if "ST", make  subplots.
#'
#' @return Data frame with added subplot names, parent spot indices, subplot indice,
#' and offset row/column coordinates.
#'
#' @keywords  internal

make_subplots <- function(sce,  platform = c("Visium", "ST")){

  n <- ncol(sce)
  ## Number of subplots
  subplots <- ifelse(platform == "Visium", 6, 9)

  offsets <- subplot.offsets(platform)

  positions <- as.data.frame(colData(sce)[, c("row", "col")])
  positions$spot.index <- c(1 : nrow(positions))

  ## Segment spots into subspots
  sub.positions <- positions[rep(seq_len(n), subplots), ]

  ## Index of the spots
  colnames(sub.positions)[colnames(sub.positions) %in% c("row", "col")] <- c("spot.row", "spot.col")

  ## Array row and col of the subplots
  sub.positions$sub.row <- sub.positions$spot.row + rep(offsets[, 1], each = n)
  sub.positions$sub.col <- sub.positions$spot.col + rep(offsets[, 2], each = n)

  ## Index of the subspot within its parent spot
  sub.positions$sub.index <- rep(seq_len(subplots), each = n)

  ## Rename the row
  spot_idx <- sub.positions$spot.index
  sub_idx <- sub.positions$sub.index
  rownames(sub.positions) <- paste("subspot_", spot_idx, ".", sub_idx, sep = "")

  return(sub.positions)
}

#' Compute pairwise distances between all subspots and return list of neighbors
#' for each subspot.
#'
#' @param sub.positions  Data frame with added subplot names, parent spot indices, subplot indice,
#' and offset row/column coordinates.
#' @param platform if "Visium",  subplots with distances less than 0.7 are neighbors.
#'
#' @return \code{sub.neighbors} a list of neighbor indices for each subspot
#'
#' @keywords internal

find.neighbors.sub <- function(sub.positions, platform){
  if(platform == "ST"){
    radius <- 0.35
    pdist <- as.matrix(stats::dist(sub.positions[, c("sub.row", "sub.col")], method = "manhattan"))
    neighbors <- (pdist <= radius & pdist > 0)
    sub.neighbors <- sapply(seq_len(nrow(sub.positions)), function(x) as.vector(which(neighbors[x, ])))
  }
  if(platform == "Visium"){
    radius <- 0.7
    pdist <- as.matrix(stats::dist(sub.positions[, c("sub.row", "sub.col")], method = "manhattan"))
    neighbors <- (pdist <= radius & pdist > 0)
    sub.neighbors <- sapply(seq_len(nrow(sub.positions)),  function(x) as.vector(which(neighbors[x, ])))

    for(i in 1 : nrow(sub.positions)){
      get.neighbors <-  c(sub.neighbors[[i]])
      get.n <- length(get.neighbors)

      neighbors.dist <- as.matrix(stats::dist(sub.positions[get.neighbors, c("sub.row", "sub.col")], method = "manhattan")) > 1
      neighbors.dist <- colSums(neighbors.dist)

      if(length(which(neighbors.dist == get.n - 1)) > 0){
        surplus <- which(neighbors.dist == get.n - 1)
        sub.neighbors[[i]] <- sub.neighbors[[i]][-surplus]
      }
    }
  }

  return(sub.neighbors)
}

#' Get the subspots cluster assignments of the BayesSpace
#'
#
#' @param Y  PCs of original sce
#' @param q The number of clusters
#' @param sub.neighbors A list of neighbor indices for each subspot
#' @param platform Spatial transcriptomic platform.
#' @param init Initial cluster assignments for spots.
#' @param mu0 Prior mean hyperparameter for mu. If not provided, mu0 is set to
#'   the mean of PCs over all spots.
#' @param Lambda0  Prior precision hyperparam for mu. If not provided, lambda0
#'   is set to a diagonal matrix \eqn{0.01 I}.
#' @param gamma Smoothing parameter.
#' @param alpha Hyperparameter for Wishart distributed precision lambda.
#' @param beta Hyperparameter for Wishart distributed precision lambda.
#' @param jitter_scale Controls the amount of jittering. Small amounts of
#'   jittering are more likely to be accepted but result in exploring the space
#'   more slowly. We suggest tuning \code{jitter_scale} so that Ychange is on
#'   average around 25\%-40\%.
#' @param jitter_prior  Scale factor for the prior variance, parameterized as the
#'   proportion  of the mean variance of the PCs.
#' @param ncores The number of cores used for parallel calculation.
#' @param n.rep The number of MCMC iterations.
#' @param save_para If true, save the list of clustering parameters.
#'
#' @return  If \code{save_para = FALSE}, return the vector of subspots cluster assignments.
#' If \code{save_para = TRUE}, return the list of clustering parameter values
#'  apart from the vector of subspots cluster assignments.
#'
#' @keywords internal
#'
#' @importFrom stats cov
#' @importFrom mvtnorm dmvnorm
#' @importFrom stats rWishart runif
#' @importFrom parallel makeCluster parLapply stopCluster
#' @importFrom mvtnorm rmvnorm
enhanced.results <- function(Y, q, sub.neighbors, platform, init,
                             mu0 = colMeans(Y), Lambda0 = diag(ncol(Y)) * 0.01,
                             gamma = 3, alpha = 1, beta = 0.01,
                             jitter_scale, jitter_prior,
                             ncores = 2, n.rep = 1000, save_para = FALSE){

  subplots <- ifelse(platform == "Visium", 6, 9)


  d <- ncol(Y)
  n <- nrow(Y)

  Y.sub <- Y[rep(seq_len(n), subplots), ]
  n.sub <- nrow(Y.sub)

  ## The number of subplots
  n.sub <- nrow(Y.sub)
  ## The clusters of subplots
  z.sub <- rep(init, subplots)

  ## Initialize parameters
  ## Initialize mu_k
  ## The ith row of mu is mu_k
  mu <- matrix(0, nrow = q, ncol =  d)
  ## Initialize Lambda
  Lambda <- diag(d) * beta
  ## Initialize wi
  w <- rep(1, n.sub)

  ## The fixed degrees-of-freedom parameter v
  v <- 4

  ## Jumping scaling factor
  ## error_var= scaling / d
  error_var <- diag(d) * jitter_scale / d

  ## Scale factor for the prior variance
  ## Be parameterized as the proportion (default = 0.3) of the mean variance of the PCs.
  c <- jitter_prior /  mean(diag(cov(Y)))

  ## Enhanced Cluster
  for(i in 1 : n.rep){
    ## Update y_ij
    for(j in  1 : n){

      vector.sub <- 0 + c(0 : (subplots - 1)) * n

      Y0 <- Y[j, ]
      Y.prev <- Y.sub[c(j + vector.sub), ]

      ## Error value
      error <- rmvnorm(subplots, rep(0, d), error_var)
      error_mean <- colSums(error) / subplots
      error <- t(t(error) - error_mean)
      Y.new <- Y.prev + error
      Sigma <- solve(Lambda)

      for(r in 1 : subplots){

        label <- z.sub[j + vector.sub[r]]

        likeli_prev <- mvtnorm::dmvnorm(Y.prev[r, ], mu[label, ], Sigma / w[j + vector.sub[r]])
        likeli_new <- mvtnorm:: dmvnorm(Y.new[r, ], mu[label, ], Sigma / w[j + vector.sub[r]])

        prior_prev <- exp(- c / 2 * sum((Y.prev[r, ] - Y0) ^ 2))
        prior_new <- exp(- c / 2 * sum((Y.new[r, ] - Y0) ^ 2))

        prob <- (likeli_new / likeli_prev) * (prior_new / prior_prev)

        ## Accept or reject
        prob <- min(prob, 1)
        u <- runif(1, 0, 1)
        if (u <=  prob)
          Y.sub[(j + vector.sub[r]), ] <- Y.new[r, ]
      }
    }

    ## Update mu
    mu <-t((sapply(seq_len(q),
                   function(x) as.vector(update.mu(j = x, z = z.sub, Y =Y.sub, mu0 = mu0, Lambda0 = Lambda0, Lambda = Lambda, w = w, mu = mu)))))

    ## Update Lambda
    df_1 <- alpha + n.sub
    W_1 <- diag(d) * beta
    for (j in 1 : q){
      cluster.label <- which(z.sub == j)
      S_weight <- crossprod(((Y.sub[cluster.label, ] - mu[j, ]) * w[cluster.label]), (Y.sub[cluster.label, ] - mu[j, ]))

      W_1 <- W_1 + S_weight
    }
    Lambda <- rWishart(1, df_1, solve(W_1))[, , 1]

    if(ncores > 1){
      ## Set core numbers
      cl <- makeCluster(ncores)

      ## Assigning tasks to each core for parallel computing
      ## Update w
      w <- parLapply(cl, 1 : n.sub, update.w, z = z.sub, Y = Y.sub, d = d, v = v, Lambda = Lambda, mu = mu, w = w)
      w <- unlist(w)

      ## Use M-H algorithm to update z
      z.sub <- parLapply(cl, 1 : n.sub, update.z, z = z.sub, Y = Y.sub, q = q, spots.neighbors = sub.neighbors, gamma = gamma, Lambda = Lambda, mu = mu, w = w)

      ## End parallel computing
      stopCluster(cl)

      z.sub <- unlist(z.sub)
    }

    if (ncores == 1){
      ## Update w_i
      w <- sapply(seq_len(n.sub), function(x) as.vector(update.w(j = x, z = z.sub, Y = Y.sub, d, v, Lambda, mu, w)))

      ## Use M-H algorithm to update z_j
      ## New cluster label for z_j
      z.sub <- sapply(seq_len(n.sub), function(x) as.vector(update.z(j = x, z = z.sub, Y = Y.sub, q, sub.neighbors, gamma, Lambda, mu, w)))
    }

  }

  if(save_para){
    return(list(z.sub = z.sub, Y.sub = Y.sub, w = w, Lambda = Lambda, mu = mu))

  }
  else
    return(list(z.sub = z.sub, Y.sub = Y.sub))

}

#' Enhance spot resolution
#'
#' Enhanced clustering of a spatial expression dataset to subspot resolution.
#'
#' @param sce A SingleCellExperiment object containing the spatial data.
#' @param q The number of clusters.
#' @param d Number of top principal components to use when clustering.
#' @param platform Number of top principal components to use when clustering.
#' @param init Initial cluster assignments for spots.
#' @param init.method If \code{init} is not provided, cluster the top \code{d}
#'   PCs with this method to obtain initial cluster assignments.
#' @param mu0 Prior mean hyperparameter for mu. If not provided, mu0 is set to
#'   the mean of PCs over all spots.
#' @param Lambda0  Prior precision hyperparam for mu. If not provided, lambda0
#'   is set to a diagonal matrix \eqn{0.01 I}.
#' @param gamma Smoothing parameter.
#' @param alpha Hyperparameter for Wishart distributed precision lambda.
#' @param beta  Hyperparameter for Wishart distributed precision lambda.
#' @param jitter_scale Controls the amount of jittering. Small amounts of
#'   jittering are more likely to be accepted but result in exploring the space
#'   more slowly. We suggest tuning \code{jitter_scale} so that Ychange is on
#'   average around 25\%-40\%.
#' @param jitter_prior Scale factor for the prior variance, parameterized as the
#'   proportion of the mean variance of the PCs.
#' @param n.rep The number of MCMC iterations.
#' @param ncores The number of cores used for parallel calculation.
#' @param save_para If true, save the list of clustering parameters.
#'
#' @return  Returns a new SingleCellExperiment object. By default, the
#'   \code{assays} of this object are empty, and the enhanced resolution PCs
#'   are stored as a reduced dimensionality result accessible with
#'   \code{reducedDim(sce, 'PCA')}.
#'
#' @details
#' The enhanced \code{SingleCellExperiment} has most of the properties of the
#'   input SCE - \code{rowData}, \code{colData}, \code{reducedDims} - but does
#'   not include expression data in \code{counts} or \code{logcounts}. To impute
#'   enhanced expression vectors, please use [EnhanceFeatures()] after
#'   running \code{spatialEnhance}.
#'
#' The \code{colData} of the enhanced \code{SingleCellExperiment} includes the
#'   following columns to permit referencing the subspots in spatial context and
#'   linking back to the original spots:
#'   \itemize{
#'   \item \code{spot.idx}: Index of the spot this subspot belongs to (with
#'     respect to the input SCE).
#'   \item \code{subspot.idx}: Index of the subspot within its parent spot.
#'   \item \code{spot.row}: Array row of the subspot's parent spot.
#'   \item \code{spot.col}: Array col of the subspot's parent spot.
#'   \item \code{row}: Array row of the subspot. This is the parent spot's row
#'     plus an offset based on the subspot's position within the spot.
#'   \item \code{col}: Array col of the subspot. This is the parent spot's col
#'     plus an offset based on the subspot's position within the spot.
#'     }
#'
#' (optional)
#' The BayesSpace of the enhanced \code{SingleCellExperiment} includes
#' the list of clustering parameters.
#'
#' @export
#'
#' @importFrom SingleCellExperiment reducedDim  reducedDim<-  SingleCellExperiment
#' @importFrom SummarizedExperiment colData colData rowData
#' @importFrom S4Vectors metadata metadata<-
#' @importFrom assertthat assert_that

EnhancedCluster <- function(sce, q, d, platform = c("Visium", "ST"),
                            init = NULL, init.method = c("spatialcluster", "mclust", "kmeans"),
                            mu0 = NULL, Lambda0 = NULL,
                            gamma = 3, alpha = 1, beta = 0.01,
                            jitter_scale = 5, jitter_prior = 0.3,
                            n.rep = 1000, ncores = 2, save_para = FALSE){

  platform <- match.arg(platform)

  ## Get PCs
  Y <- reducedDim(sce, "PCA")
  Y <- as.matrix(Y[,seq_len(d)])

  ## Initialize cluster
  if (is.null(init)) {
    init.method <- match.arg(init.method)
    if (init.method == "spatialcluster") {
      init <- sce$spatial.cluster
    } else
      init <- init.cluster(Y, q, init, init.method)
  }


  ## Make subplots and find neighbors
  sub.positions <- make_subplots(sce, platform = platform)
  sub.neighbors <-  find.neighbors.sub(sub.positions, platform)

  ## Initialize hyperparameters
  ## Set mu0 to be the empirical mean vector of the data, which is generally the zero vector for PCA
  if (is.null(mu0))
    mu0 <- rep(0, d)
  ## Lambda0 is set to 0.01 times the identity matrix
  if (is.null(Lambda0))
    Lambda0 <- diag(d) * 0.01

  ## gamma = 3 for Visium and gamma = 2 for ST
  if (is.null(gamma))
    gamma <- ifelse(platform == "Visium", 3, 2)

  enhanced <- enhanced.results(Y, q, sub.neighbors, platform, init,
                               mu0, Lambda0,
                               gamma, alpha, beta,
                               jitter_scale, jitter_prior,
                               ncores, n.rep, save_para)

  enhanced.cluster <- enhanced$z.sub
  counts.pca <- enhanced$Y.sub

  rownames(counts.pca) <- rownames(sub.positions)

  cdata <- data.frame(sub.positions, enhanced.cluster)
  colnames(cdata) <- c(colnames(sub.positions), "spatial.cluster")

  enhanced.sce <- SingleCellExperiment(assays = list(), rowData = rowData(sce), colData = cdata)

  reducedDim(enhanced.sce, "PCA") <- counts.pca

  ## Add metadata to new SingleCellExperiment object
  metadata(enhanced.sce)$BayesSpace.data <- list()
  metadata(enhanced.sce)$BayesSpace.data$platform <- platform
  metadata(enhanced.sce)$BayesSpace.data$is.enhanced <- TRUE

  ## Add parameters
  if(save_para){
    para <- list(mu = enhanced$mu, w = enhanced$w, Lambda = enhanced$Lambda)
    metadata(enhanced.sce)$para <- para
  }


  return(enhanced.sce)

}


