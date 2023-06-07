

#' Initialize cluster
#'
#' @param Y PCs of SingleCellExperiment
#' @param q The number of clusters
#' @param init Vector of initial cluster assignments
#' @param init.method Initialization clustering algorithm
#'
#' @return Vector of initial cluster assignments
#'
#' @keywords internal
#' @importFrom stats kmeans
#' @importFrom mclust Mclust mclustBIC

init.cluster <- function(Y, q, init = NULL, init.method = c("mclust", "kmeans")){

  if(is.null(init)){
    init.method <- match.arg(init.method)
    if (init.method == "mclust"){
      init <- Mclust(Y, G = q, "EEE", verbose = FALSE)$classification
    }
    else
      if (init.method == "kmeans"){
        init <- kmeans(Y, q)$cluster
      }
  }

  return(init)
}


#' Find neighboring spots based on array coordinates
#'
#'
#' @param sce SingleCellExpreriment
#' @param platform If "Visium", select six neighboring spots around center; if
#'   "ST", select four adjacent spots.
#'
#' @return \code{spots.neighbors} a list of neighbor indices for each spot
#' @keywords internal
#'
find.neighbors <- function(sce, platform){
  if (platform == "ST"){
    ## ST platform: squares
    ## Neighbors: below, above, right, left
    origin <- data.frame(origin.row = c(-1, 1, 0, 0),
                         origin.col = c(0, 0, 1, -1))
  }
  if (platform == "Visium")
    ## Visium platform: hexagon
    ## Neighbors: left, right, upper left, upper right, lower left, lower right
    origin <- data.frame(origin.row = c(0, 0, 1, 1, -1, -1),
                         origin.col =  c(-2, 2, -1, 1, -1, 1))


  ## Get coordinates
  positions <- as.data.frame(colData(sce)[, c("row", "col")])
  positions$spot.index <- c(1 : nrow(positions))
  n <- nrow(positions)

  ## Get each possible spot neighbor
  neighbor.cor <- merge(positions, origin)
  neighbor.cor$neighbor.row <- neighbor.cor$row + neighbor.cor$origin.row
  neighbor.cor$neighbor.col <- neighbor.cor$col + neighbor.cor$origin.col

  ## Get each possible spot neighbor
  neighbor.cor <- merge(positions, origin)
  neighbor.cor$neighbor.row <- neighbor.cor$row + neighbor.cor$origin.row
  neighbor.cor$neighbor.col <- neighbor.cor$col + neighbor.cor$origin.col

  ## Match each possible spot neighbor with existing spots
  spots.neighbors <- merge(neighbor.cor, positions,
                           by.x = c("neighbor.row", "neighbor.col"), by.y = c("row", "col"),
                           suffixes = c(".spots", ".neighbors"))
  spots.neighbors <- lapply(seq_len(n),
                            function(x) as.vector(spots.neighbors[which(spots.neighbors$spot.index.spots == x), ]$spot.index.neighbors))

  return(spots.neighbors)
}


#' Update parameter mu
#'
#' @param j Update mean of cluster j
#' @param z Vector of cluster assignments
#' @param Y PCs of SingleCellExperiment
#' @param mu0  Prior mean hyperparameter for mu. If not provided, mu0 is set to
#'   the mean of PCs over all spots.
#' @param Lambda0 Prior precision hyperparam for mu. If not provided, lambda0
#'   is set to a diagonal matrix \eqn{0.01 I}.
#' @param Lambda Precision matrix
#' @param w Proportion factor
#' @param mu Mean of cluster
#'
#' @return New mean of cluster j
#'
#' @keywords internal
#'
#' @importFrom MASS mvrnorm

update.mu <- function(j, z, Y, mu0, Lambda0, Lambda, w, mu){
  cluster.label <- which(z == j)
  Lambda_1 <- Lambda0 + Lambda * sum(w[cluster.label])
  Sigma_1 <- solve(Lambda_1)
  if(length(cluster.label) > 1){
    mu_1 <- Sigma_1  %*% (Lambda0  %*% mu0 +  Lambda %*% colSums(w[cluster.label] * Y[cluster.label, ]))
  }else
      if (length(cluster.label) == 1){
        mu_1 <- Sigma_1  %*% (Lambda0  %*% mu0 +  Lambda %*% (w[cluster.label] * Y[cluster.label, ]))
      }

  if(length(cluster.label) >  0){
    mu[j, ] <- mvrnorm(n = 1, mu = mu_1, Sigma = solve(Lambda_1))
  }else{
    mu[j, ] <- mu[j, ]
  }
}



#' Update parameter w
#'
#' @param j Update jth proportion factor
#' @param z Vector of cluster assignments
#' @param Y  Vector of cluster assignments
#' @param d Dimention of PCs
#' @param v A fixed degree of freedom
#' @param Lambda Precision matrix
#' @param mu Mean of cluster
#' @param w Proportion factor
#'
#' @return New jth proportion factor
#'
#' @keywords internal
#'
#' @importFrom stats rgamma
update.w <- function(j, z, Y, d, v, Lambda, mu, w){
  shape_1 <- (d + v) / 2
  S <- (Y[j, ] - mu[z[j], ]) %*% Lambda %*% (Y[j, ] - mu[z[j], ])
  rate_1 <- (v + S) / 2
  w[j] <- rgamma(n = 1, shape = shape_1, rate = rate_1)
}


#' Update parameter z
#'
#' @param j Update jth cluster assignments
#' @param z Vector of cluster assignments
#' @param Y Vector of cluster assignments
#' @param q Number of clusters
#' @param spots.neighbors List of neighbor indices for each spot
#' @param gamma Smoothing parameter. Defaults to 2 for \code{platform="ST"} and
#'   3 for \code{platform="Visium"}.
#' @param Lambda  Precision matrix
#' @param mu Mean of cluster
#' @param w Proportion factor
#'
#' @return New jth cluster assignments
#'
#' @keywords internal
#'
#' @importFrom mvtnorm dmvnorm
#' @importFrom stats runif


update.z <- function(j, z, Y, q, spots.neighbors, gamma, Lambda, mu, w){

  z_new <- sample(c(1 : q)[-z[j]], size = 1)
  neighbors <- spots.neighbors[[j]]
  z.neighbors <- z[neighbors]
  n.neighbors <- length(neighbors)

  Sigma <- solve(Lambda)

  if(n.neighbors > 0){


    ## Potts model prior
    pi_j <-  gamma / n.neighbors * 2 * length(which(z.neighbors == z[j]))
    pi_new <- gamma / n.neighbors * 2 * length(which(z.neighbors == z_new))

    ## Likelihoood

    density_j <- mvtnorm::dmvnorm(Y[j, ], mu[z[j], ], Sigma / w[j], log = TRUE)
    density_new <- mvtnorm::dmvnorm(Y[j, ], mu[z_new, ], Sigma / w[j], log = TRUE)

    ## Calculate the ratio of posterior distributions
    prob <- exp((density_new + pi_new) - (density_j + pi_j))

  }
  else {
    ## Likelihoood
    density_j <- mvtnorm::dmvnorm(Y[j, ], mu[z[j], ], Sigma / w[j])
    density_new <- mvtnorm::dmvnorm(Y[j, ], mu[z_new, ], Sigma / w[j])

    ## Calculate the ratio of  distributions
    prob <- (density_new / density_j)
  }

  ## Accept or reject
  prob <- min(prob, 1)
  u <- runif(1, 0, 1)
  if (u <=  prob)
    z[j] <- z_new

  return(z[j])
}


#' Get the cluster assignments of the BayesSpace
#'
#' @param Y The PCs of SingleCellExperiment
#' @param q The numbers of clusters
#' @param spots.neighbors List of neighbor indices for each spot
#' @param init Vector of initial cluster assignments
#' @param mu0  Prior mean hyperparameter for mu. If not provided, mu0 is set to
#'   the mean of PCs over all spots.
#' @param Lambda0  Prior precision hyperparam for mu. If not provided, lambda0
#'   is set to a diagonal matrix \eqn{0.01 I}.
#' @param gamma  Smoothing parameter. Defaults to 2 for \code{platform="ST"} and
#'   3 for \code{platform="Visium"}.
#' @param alpha Hyperparameter for Wishart distributed precision lambda.
#' @param beta Hyperparameter for Wishart distributed precision lambda.
#' @param ncores The number of cores used for parallel calculation.
#' @param n.rep The number of MCMC iterations
#' @param save_para If true, return the parameter mu, w and Lambda apart from z.
#'
#' @return If \code{save_para = FALSE}, return the vector of cluster assignments.
#' If \code{save_para = TRUE}, return the list of clustering parameter values apart from the vector of cluster assignments.
#'
#' @keywords internal
#'
#' @importFrom stats rWishart
#' @importFrom parallel makeCluster parLapply stopCluster
cluster.results <- function(Y, q, spots.neighbors, init, mu0 = colMeans(Y),
                            Lambda0 = diag(ncol(Y)) * 0.01, gamma = 3,
                            alpha = 1, beta = 0.01, ncores = 2, n.rep =1000, save_para = FALSE){

  d <- ncol(Y)
  n <- nrow(Y)

  # Initialize cluster
  z <- init

  # Initialize parameters
  # Initialize mu_k
  # The ith row of mu is mu_k
  mu <- matrix(0, nrow = q, ncol =  d)
  # Initialize Lambda
  Lambda <- diag(d) * beta
  # Initialize wi
  w <- rep(1, n)

  # The fixed degrees-of-freedom parameter v
  v <- 4

  for(i in 1 : n.rep){

    # Use Gibbs sampling update most of parameters
    # Update mu
    mu <-t((sapply(seq_len(q),
                   function(x) as.vector(update.mu(j = x, z = z, Y = Y,
                                                   mu0 = mu0, Lambda0 = Lambda0,
                                                   Lambda = Lambda, w = w, mu = mu)))))

    # Update Lambda
    df_1 <- alpha + n
    W_1 <- diag(d) * beta

    for (j in 1 : q){

      cluster.label <- which(z == j)
      if(length(cluster.label) > 1){
        S_weight <- crossprod(((Y[cluster.label, ] - mu[j, ]) * w[cluster.label]), (Y[cluster.label, ] - mu[j, ]))
        W_1 <- W_1 + S_weight
      }
      else if(length(cluster.label) == 1){
        S_weight <- ((Y[cluster.label, ] - mu[j, ]) * w[cluster.label])  %o% (Y[cluster.label, ] - mu[j, ])
        W_1 <- W_1 + S_weight
      }
     }


    Lambda <- rWishart(1, df_1, solve(W_1))[, , 1]


    if(ncores > 1){
      # Set core numbers
      cl <- makeCluster(ncores)

      # Assigning tasks to each core for parallel computing.
      # Update w
      w <- parLapply(cl, 1 : n, update.w, z = z, Y = Y, d = d, v = v,
                     Lambda = Lambda, mu = mu, w = w)
      w <- unlist(w)

      # Use M-H algorithm to update z
      z <- parLapply(cl, 1 : n, update.z, z = z, Y = Y, q = q, spots.neighbors = spots.neighbors,
                     gamma = gamma, Lambda = Lambda, mu = mu, w = w)

      # End parallel computing
      stopCluster(cl)

      z <- unlist(z)
    }
    else if(ncores == 1){
      ## Update w
      w <- sapply(seq_len(n), function(x) as.vector(update.w(j = x, z = z, Y = Y,
                                                             d, v, Lambda, mu, w)))

      ## Use M-H algorithm to update z
      z <- sapply(seq_len(n), function(x) as.vector(update.z(j = x, z = z, Y = Y, q,
                                                             spots.neighbors, gamma, Lambda, mu, w)))
    }

  }
  if(save_para){
    return(list(z = z, w = w, Lambda = Lambda, mu = mu))

  }
  else
    return(z)

}


#' Spatial clustering
#'
#' Cluster a spatial expression dataset.
#'
#' @param sce A SingleCellExperiment object containing the spatial data.
#' @param q The numbers of clusters.
#' @param d The Number of top principal components to use when clustering.
#' @param platform Spatial transcriptomic platform.
#' @param init Vector of initial cluster assignments
#' @param init.method If \code{init} is not provided, cluster the top \code{d}
#'   PCs with this method to obtain initial cluster assignments.
#' @param mu0  Prior mean hyperparameter for mu. If not provided, mu0 is set to
#'   the mean of PCs over all spots.
#' @param Lambda0 Prior precision hyperparam for mu. If not provided, lambda0
#'   is set to a diagonal matrix \eqn{0.01 I}.
#' @param alpha Hyperparameter for Wishart distributed precision lambda.
#' @param beta Hyperparameter for Wishart distributed precision lambda.
#' @param gamma  Smoothing parameter. Defaults to 2 for \code{platform="ST"} and
#'   3 for \code{platform="Visium"}.
#' @param ncores The number of cores used for parallel calculation.
#' @param n.rep The number of MCMC iterations.
#' @param save_para  If true, save the parameter mu, w and Lambda apart from z.
#'
#' @return  Returns a modified \code{sce} with cluster assignments
#' stored in \code{colData} under the name \code{spatial.cluster},
#' and (optional) the list of clustering parameters,
#' stored BayesSpace metadata under the name \code{para}.
#'
#' @export
#'
#' @importFrom SingleCellExperiment reducedDim
#' @importFrom SummarizedExperiment colData colData<-
#' @importFrom S4Vectors metadata metadata<-
#' @importFrom assertthat assert_that
SpatialCluster <- function(sce, q, d, platform = c("Visium", "ST"), init = NULL, init.method = c("mclust", "kmeans"),
                           mu0 = NULL, Lambda0 = NULL,
                           alpha = 1, beta = 0.01, gamma = NULL, ncores = 2,
                           n.rep = 1000, save_para = FALSE) {

  platform <- match.arg(platform)

  ## Get PCs
  Y <- reducedDim(sce, "PCA")
  Y <- as.matrix(Y[,seq_len(d)])

  ## Get the index of neighbors
  spots.neighbors <- find.neighbors(sce, platform)

  ## Initialize cluster
  init <- init.cluster(Y, q, init, init.method)
  colData(sce)$init.cluster <- init

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


  results <- cluster.results(Y, q = q, spots.neighbors, init, mu0 = mu0, Lambda0 = Lambda0,
                             gamma = gamma, alpha = alpha, beta = beta, ncores = ncores,
                             n.rep = n.rep, save_para = save_para)

  if(save_para){
    para <- list(mu = results$mu, w = results$w, Lambda = results$Lambda)
    metadata(sce)$para <- para

    colData(sce)$spatial.cluster <- results$z
  }
  else{
    colData(sce)$spatial.cluster <- results
  }


  return(sce)
}


