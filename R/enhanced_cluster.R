library(mvtnorm)
library(MASS)
library(parallel)

library(scuttle)
library(scater)



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


make_subplots <- function(sce,  platform = c("Visium", "ST")){
  
  n <- ncol(sce)
  #Number of subplots
  subplots <- ifelse(platform == "Visium", 6, 9)
  
  offsets <- subplot.offsets(platform)
  
  positions <- as.data.frame(colData(sce)[, c("row", "col")])
  positions$spot.index <- c(1 : nrow(positions))
  
  #Segment spots into subspots
  sub.positions <- positions[rep(seq_len(n), subplots), ]
  
  #index of the spots 
  colnames(sub.positions)[colnames(sub.positions) %in% c("row", "col")] <- c("spot.row", "spot.col")
  
  #array row and col of the subplots
  sub.positions$sub.row <- sub.positions$spot.row + rep(offsets[, 1], each = n)
  sub.positions$sub.col <- sub.positions$spot.col + rep(offsets[, 2], each = n)
  
  #index of the subspot within its parent spot
  sub.positions$sub.index <- rep(seq_len(subplots), each = n)
  
  #rename the row
  spot_idx <- sub.positions$spot.index
  sub_idx <- sub.positions$sub.index
  rownames(sub.positions) <- paste("subspot_", spot_idx, ".", sub_idx, sep = "")
  
  return(sub.positions)
}

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

enhanced.results <- function(Y, q, sub.neighbors, platform, init, 
                             mu0 = colMeans(Y), Lambda0 = diag(ncol(Y)) * 0.01, 
                             gamma = 3, alpha = 1, beta = 0.01,
                             jitter_scale, jitter_prior,
                             ncores = 2, n.rep = 1000, save_para = FALSE){
  
  d <- ncol(Y)
  n <- nrow(Y)
  
  subplots <- ifelse(platform == "Visium", 6, 9)
  
  Y.sub <- Y[rep(seq_len(n), subplots), ]
  
  #The number of subplots
  n.sub <- nrow(Y.sub)
  #The clusters of subplots
  z.sub <- rep(init, subplots)
  
  #Initialize parameters
  #Initialize mu_k
  #The ith row of mu is mu_k
  mu <- matrix(0, nrow = q, ncol =  d)
  #Initialize Lambda
  Lambda <- diag(d) * beta
  #Initialize wi
  w <- rep(1, n.sub)
  
  #The fixed degrees-of-freedom parameter v
  v <- 4
  
  #Jumping scaling factor
  #error_var= scaling / d
  error_var <- diag(d) * jitter_scale / d
  
  #Scale factor for the prior variance
  #Be parameterized as the proportion (default = 0.3) of the mean variance of the PCs.
  c <- jitter_prior /  mean(diag(cov(Y)))
  
  #Enhanced Cluster
  for(i in 1 : n.rep){
    #update y_ij
    for(j in  1 : n){
      
      vector.sub <- 0 + c(0 : (subplots - 1)) * n
      
      Y0 <- Y[j, ]
      Y.prev <- Y.sub[c(j + vector.sub), ]
      
      #Error value
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
        
        #accept or reject
        prob <- min(prob, 1)
        u <- runif(1, 0, 1)
        if (u <=  prob)
          Y.sub[(j + vector.sub[r]), ] <- Y.new[r, ]
      }
    }
    
    #update mu
    mu <-t((sapply(seq_len(q), 
                   function(x) as.vector(update.mu(j = x, z = z.sub, Y =Y.sub, mu0 = mu0, Lambda0 = Lambda0, Lambda = Lambda, w = w, mu = mu)))))
    
    #update Lambda
    df_1 <- alpha + n.sub
    W_1 <- diag(d) * beta
    for (j in 1 : q){
      cluster.label <- which(z.sub == j)
      S_weight <- crossprod(((Y.sub[cluster.label, ] - mu[j, ]) * w[cluster.label]), (Y.sub[cluster.label, ] - mu[j, ]))
      
      W_1 <- W_1 + S_weight
    }
    Lambda <- rWishart(1, df_1, solve(W_1))[, , 1]
    
    if(ncores > 1){
      # Set core numbers
      cl <- makeCluster(ncores)
      
      # Assigning tasks to each core for parallel computing
      # Update w
      w <- parLapply(cl, 1 : n.sub, update.w, z = z.sub, Y = Y.sub, d = d, v = v, Lambda = Lambda, mu = mu, w = w)
      w <- unlist(w)
      
      # Use M-H algorithm to update z
      z.sub <- parLapply(cl, 1 : n.sub, update.z, z = z.sub, Y = Y.sub, q = q, spots.neighbors = sub.neighbors, gamma = gamma, Lambda = Lambda, mu = mu, w = w)
      
      # End parallel computing
      stopCluster(cl)
      
      z.sub <- unlist(z.sub)
    }
    
    if (ncores == 1){
      #Update w_i
      w <- sapply(seq_len(n.sub), function(x) as.vector(update.w(j = x, z = z.sub, Y = Y.sub, d, v, Lambda, mu, w)))
      
      #Use M-H algorithm to update z_j
      #new cluster label for z_j
      z.sub <- sapply(seq_len(n.sub), function(x) as.vector(update.z(j = x, z = z.sub, Y = Y.sub, q, sub.neighbors, gamma, Lambda, mu, w)))
    }
    
  }
  
  if(save_para){
    return(list(z.sub = z.sub, Y.sub = Y.sub, w = w, Lambda = Lambda, mu = mu))
    
  }
  else
    return(list(z.sub = z.sub, Y.sub = Y.sub))

}

EnhancedCluster <- function(sce, q, d, platform = c("Visium", "ST"), 
                            init = NULL, init.method = c("spatialcluster", "mclust", "kmeans"),
                            mu0 = NULL, Lambda0 = NULL,
                            gamma = 3, alpha = 1, beta = 0.01, 
                            jitter_scale = 5, jitter_prior = 0.3,   
                            n.rep = 1000, ncores = 2, save_para = FALSE){
  
  platform <- match.arg(platform)
  
  #Initialize cluster
  if (is.null(init)) {
    init.method <- match.arg(init.method)
    if (init.method == "spatialcluster") {
      init <- sce$spatial.cluster
    } else 
      init <- init.cluster(inputs$PCs, q, init, init.method)
  }
  
  #Get PCs
  Y <- reducedDim(sce, "PCA")
  Y <- as.matrix(Y[,seq_len(d)])
  
  #Make subplots and find neighbors
  sub.positions <- make_subplots(sce, platform = platform)
  sub.neighbors <-  find.neighbors.sub(sub.positions, platform)
  
  #Initialize hyperparameters
  #Set mu0 to be the empirical mean vector of the data, which is generally the zero vector for PCA
  if (is.null(mu0))
    mu0 <- rep(0, d)
  #Lambda0 is set to 0.01 times the identity matrix
  if (is.null(Lambda0))
    Lambda0 <- diag(d) * 0.01
  
  # gamma = 3 for Visium and gamma = 2 for ST
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
  
  #Add metadata to new SingleCellExperiment object
  metadata(enhanced.sce)$BayesSpace.data <- list()
  metadata(enhanced.sce)$BayesSpace.data$platform <- platform
  metadata(enhanced.sce)$BayesSpace.data$is.enhanced <- TRUE
  
  #Add parameters
  if(save_para){
    para <- list(mu = enhanced$mu, w = enhanced$w, Lambda = enhanced$Lambda)
    metadata(enhanced.sce)$para <- para
  }
  
    
  return(enhanced.sce)
  
}


