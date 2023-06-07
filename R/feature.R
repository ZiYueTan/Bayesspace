

#' Predict feature vectors from enhanced PCs by linear model
#'
#' @param X  Original PCs
#' @param X.enhanced Enhanced PCs
#' @param Y Original log expression.
#' @param feature_names  List of genes/features to predict expression/values for.
#'
#' @return Vector of enhanced features.
#'
#' @keywords internal
#'
#' @importFrom stats lm predict
lm.feature <- function(X, X.enhanced, Y, feature_names) {
  r.squared <- numeric(length(feature_names))
  names(r.squared) <- feature_names

  Y.enhanced <- matrix(nrow=length(feature_names), ncol=nrow(X.enhanced))
  rownames(Y.enhanced) <- feature_names
  colnames(Y.enhanced) <- rownames(X.enhanced)

  for (feature in feature_names) {
    fit <- lm(Y[feature, ] ~ ., data=X)
    r.squared[feature] <- summary(fit)$r.squared
    Y.enhanced[feature, ] <- predict(fit, newdata=X.enhanced)
  }

  diagnostic <- list("r.squared" = r.squared)
  attr(Y.enhanced, "diagnostic") <- diagnostic
  Y.enhanced
}


#' Predict feature vectors from enhanced PCs by xgboost
#'
#' @param X  Original PCs
#' @param X.enhanced Enhanced PCs
#' @param Y Original log expression.
#' @param feature_names List of genes/features to predict expression/values for.
#' @param nrounds  Nonnegative integer to set the \code{nrounds} parameter for xgboost.
#'  If \code{nrounds} is set to 0, the parameter will be tuned using a train-test split.
#'  We recommend tuning \code{nrounds} for improved feature prediction, but note this will increase
#'  runtime.
#' @param train.n Number of spots to use in the training dataset for tuning
#'   nrounds. By default, 2/3 the total number of spots are used.
#'
#'
#' @return Vector of enhanced features.
#'
#' @keywords internal
#'
#' @importFrom xgboost xgboost xgb.DMatrix xgb.train
#'
xgboost.feature <- function(X, X.enhanced, Y, feature_names, nrounds, train.n){

  Y.enhanced <- matrix(nrow = length(feature_names), ncol = nrow(X.enhanced))
  rownames(Y.enhanced) <- feature_names
  colnames(Y.enhanced) <- rownames(X.enhanced)

  rmse <- numeric(length(feature_names))
  names(rmse) <- feature_names

  train.index <- sample(seq_len(ncol(Y)), train.n)

  default.nround <- nrounds

  for(feature in feature_names){
    nrounds <- default.nround

    ## If nrounds is set to 0, the parameter will be tuned using a train-test split.
    if (nrounds == 0){
      data.train <- xgb.DMatrix(data = X[train.index, ],label = Y[feature, train.index])
      data.test <- xgb.DMatrix(data = X[-train.index,], label = Y[feature, -train.index])

      watchlist <- list(train = data.train, test = data.test)

      fit.train <- xgb.train(data = data.train, max_depth=2,
                             watchlist=watchlist, eta=0.03, nrounds=500,
                             objective="reg:squarederror",
                             verbose=FALSE)
      nrounds <- which.min(fit.train$evaluation_log$test_rmse)
    }

    fit.xgboost <- xgboost(data = X, label = Y[feature, ], max_depth=2, eta = 0.03, nrounds = nrounds, nthread = 2, objective="reg:squarederror", verbose = FALSE)

    Y.enhanced[feature, ] <- predict(fit.xgboost, X.enhanced)
    rmse[feature] <- fit.xgboost$evaluation_log$train_rmse[nrounds]

  }


  diagnostic <- list("rmse"=rmse)
  attr(Y.enhanced, "diagnostic") <- diagnostic

  return(Y.enhanced)
}


#' Get predicted features
#'
#' @param X.enhanced Enhanced PCs
#' @param X  Original PCs
#' @param Y Original log expression.
#' @param feature_names List of genes/features to predict expression/values for.
#' @param model  Model used to predict enhanced values.
#' @param nrounds  Nonnegative integer to set the \code{nrounds} parameter for xgboost.
#'  If \code{nrounds} is set to 0, the parameter will be tuned using a train-test split.
#'  We recommend tuning \code{nrounds} for improved feature prediction, but note this will increase
#'  runtime.
#' @param train.n Number of spots to use in the training dataset for tuning
#'   nrounds. By default, 2/3 the total number of spots are used.
#'
#' @return Vector of enhanced features.
#' 
#
features <- function(X.enhanced, X, Y, feature_names = rownames(Y), model = c("xgboost", "lm"), nrounds, train.n){

  if (model == "xgboost")
    Y.enhanced <- xgboost.feature(X, X.enhanced, Y, feature_names, nrounds, train.n)
  if (model == "lm")
    Y.enhanced <- lm.feature(X, X.enhanced, Y, feature_names)

  return(Y.enhanced)

}



#' Predict feature vectors from enhanced PCs.
#'
#' @param sce.enhanced  SingleCellExperiment object with enhanced PCs.
#' @param sce SingleCellExperiment object with original PCs and log expression.
#' @param feature_names List of genes/features to predict expression/values for.
#' @param model  Model used to predict enhanced values.
#' @param nrounds Nonnegative integer to set the \code{nrounds} parameter for xgboost.
#'  If \code{nrounds} is set to 0, the parameter will be tuned using a train-test split.
#'  We recommend tuning \code{nrounds} for improved feature prediction, but note this will increase
#'  runtime.
#' @param train.n Number of spots to use in the training dataset for tuning
#'   nrounds. By default, 2/3 the total number of spots are used.
#'
#' @return Th enhanced features are stored in the \code{logcounts} of
#'   \code{sce.enhanced} and the modified SingleCellExperiment object is
#'   returned.
#'
#' @export
#'
#' @importFrom SingleCellExperiment reducedDim
#' @importFrom SummarizedExperiment assay assay<- SummarizedExperiment rowData<- rowData

EnhanceFeatures <- function(sce.enhanced, sce, feature_names, model = c("xgboost", "lm"), nrounds = 0, train.n = round(ncol(sce) * 2 / 3)){

  model <- match.arg(model)

  X.enhanced  <- reducedDim(sce.enhanced, "PCA")
  X <- reducedDim(sce, "PCA")

  Y <- assay(sce, "logcounts")

  if (is.null(feature_names)){
    feature_names <- rownames(Y)
  }else{
    feature_names <- intersect(feature_names, rownames(Y))
  }

  Y.enhanced <- features(X.enhanced, X, Y, feature_names, model, nrounds, train.n)

  Y.enhanced <- pmax(Y.enhanced, 0)

  diagnostic <- attr(Y.enhanced, "diagnostic")

  if(length(feature_names) != nrow(Y)){
    Y.enhanced.full <- matrix(data = NA, nrow = nrow(Y), ncol = ncol(sce.enhanced))
    rownames(Y.enhanced.full) <- rownames(Y)
    colnames(Y.enhanced.full) <- colnames(Y.enhanced)
    Y.enhanced.full[feature_names, ] <- Y.enhanced
    assay(sce.enhanced, "logcounts") <- Y.enhanced.full

    for(name in names(diagnostic)){
      diagnostic.full <- rep(NA, nrow(sce))
      names(diagnostic.full) <- rownames(sce)
      diagnostic.full[feature_names] <- diagnostic[[name]]
      col.name <- sprintf("enhanceFeatures.%s", name)
      rowData(sce.enhanced)[[col.name]] <- diagnostic.full
    }
  }
  else
  {
    assay(sce.enhanced, "logcounts") <- Y.enhanced
    for (name in names(diagnostic)) {
      col.name <- sprintf("enhanceFeatures.%s", name)
      rowData(sce.enhanced)[[col.name]] <- diagnostic[[name]]
    }
  }

  return(sce.enhanced)
}
