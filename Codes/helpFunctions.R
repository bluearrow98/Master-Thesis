library("Rlab")
library(MASS)
library(Matrix)
library(parallel)
library(doParallel)
library(gurobi)

stabSignal <- function(p){
  k = sample(c(1,2,3,4),size = 1)
  M = matrix(rnorm(p**2)*rbern(p**2,prob = 1 / p), nrow = p)
  nonDiagAbsSum = sum(abs(M)) - sum(diag(abs(M)))
  diag(M) = - nonDiagAbsSum - abs(rnorm(p))
  return (M)
}


getPopCovar <- function(M, C, lambda = 0){
  p = dim(M)[1]
  Sigma <- solve(kronecker(diag(p), M) + kronecker(M, diag(p))
                  + lambda * diag(p**2), -c(C))
  Sigma <- matrix(Sigma, nrow = p)
  
  return(Sigma)
}

invCov <- function(Mvec,Cvec){
  # Get number of variables
  p <- sqrt(length(Mvec))

  # Get the drift matrix and volatility matrix
  M <- matrix(Mvec, nrow = p)
  C <- matrix(Cvec, nrow = p)

  # Get the recovered covariance matrix
  Sigma <- getPopCovar(M, C, lambda = 10^-5)

  # Find the inverse using decomposition
  E <- eigen(Sigma)
  Lambda <- diag(E$values^(-1))
  Gamma <- E$vectors
  invCov <- Gamma %*% Lambda %*% t(Gamma)

  return(invCov)
}

sampleX <- function(n,Sigma){
  X = mvrnorm(n, mu = rep(0, dim(Sigma)[1]), Sigma)
  return (X)
}

getDataCovar <- function(X){
  n = dim(X)[1]
  Sigma = t(X) %*% X / n
  
  return (Sigma)
}

standardize <- function(X){
  n = dim(X)[1]
  p = dim(X)[2]
  mean = matrix(rep(t(apply(X, 2, function(x) {mean(x)})),n) , ncol = p, byrow = T)
  X_shift = X - mean
  X_norm = X_shift%*%diag(1/apply(X_shift, 2, function(x) {norm(x,type = "2")}))
  return (X_shift)
}

applyBS <- function(j, ASigma, vecC) {

  n <- dim(ASigma)[1]
  p <- length(vecC)
  result <- customBS(ASigma, -vecC, j)

  # result <- bs(ASigma, -vecC, j, intercept = F,
  #             params = list(NonConvex = 2), form = 2)

  # Store objVal and constructed signal
  if (p >= n) {
    result$beta <- result$x[(n + 1):(n + p)]
    }
  else {
    result$beta <- result$x[1:p]
    }

  return(result)

}

customBS <- function(X,y,k,polish = T){

  n <- dim(X)[1]
  p <- dim(X)[2]

  # Discrete First order method
  algo2time <- system.time(algo2Sol <- proj_grad(X, y, k, polish = polish))
  M_p <- algo2Sol$bm

  ## MIO Formulation

  Mu <- 2 * norm(as.matrix(M_p), type = "i")


  ## Initial guess (Not used for cold start)
  z0 <- rep(0,p)
  z0[M_p != 0] <- 1

  model <- list()
  params <- list()
  if (p < n) {
    # Warm-start model creation (Eq 2.5)
    model$A <- spMatrix(nrow = 2, ncol = 2 * p + 1, i = c(rep(1, p), 2),
                       j = c(seq(p + 1, 2 * p + 1)),
                       x = c(rep(1, p + 1)))
    # model$A = spMatrix(nrow = 2*p+1, ncol = 2*p,
    #                    i= c(rep(seq(1,p),2),rep(seq(p+1,2*p),2),rep(2*p+1,p)),
    #                    j=c(rep(seq(1,2*p),2),seq(p+1,2*p)),
    #                    x=c(rep(-1,p),rep(-Mu,p),rep(1,p),rep(-Mu,p),rep(1,p)))
    model$Q <- spMatrix(nrow = 2 * p + 1,ncol = 2 * p + 1,
                        i = c(rep(seq(1, p), each = p)),
                        j = c(rep(seq(1, p), p)), x = c(0.5 * t(X) %*% X))
    model$obj <- c(-t(X) %*% y, rep(0, p + 1))
    model$ub <- c(rep(Mu, p), rep(1, p), Inf)
    model$lb <- c(rep(-Mu, p), rep(0, p), 0)
    model$rhs <- c(p - k, Mu)
    model$sense <- c(">=", "<=")
    model$vtype <- c(rep("C", p), rep("B", p), "C")
    model$start <- c(M_p, 1-z0, norm(as.matrix(M_p), type = "i"))
    model$genconnorm <- list(list(resvar = 2 * p + 1,
                                 vars = c(seq(1, p)),
                                 which = Inf))

    ## Creating SOS constraints
    sos <- list()
    for (j in 1:p){
      sos[[j]] <- list()
      sos[[j]]$type <- 1
      sos[[j]]$index <- c(j, (p + j))
      sos[[j]]$weights <- c(1, 2)
    }
    model$sos <- sos

  }else{
    # Warm-start model creation (Eq 2.6)
    Mu_zeta <- max(colSums(apply(abs(X), 1, sort,
                    decreasing=TRUE)[, 1:k, drop = F])) * Mu
    model$A <- spMatrix(nrow = n + 1, ncol = n + 2 * p,
                       i = c(seq(1, n), rep(seq(1, n), each = p), rep(n + 1, p)),
                       j = c(seq(1, n), rep(seq(n+1, n+p), n), seq(n + p + 1, n + 2 * p)),
                       x = c(rep(1, n), c(-t(X)), rep(1, p)))
    model$Q <- spMatrix(nrow = n + 2 * p, ncol = n + 2 * p,
                       i = seq(1, n), j = seq(1, n), x = rep(0.5, n))
    model$obj <- c(rep(0, n), -t(X) %*% y, rep(0, p))
    model$ub <- c(rep(Mu_zeta, n), rep(Mu, p), rep(1, p))
    model$lb <- c(rep(-Mu_zeta, n), rep(-Mu, p), rep(0, p))
    model$rhs <- c(rep(0, n), p - k)
    model$sense <- c(rep("=", n),">=")
    model$vtype <- c(rep("C", n), rep("C", p), rep("B", p))
    model$start <- c(X %*% M_p, M_p, 1 - z0)


    ## Creating SOS constraints
    sos <- list()
    for (j in 1:p){
      sos[[j]] <- list()
      sos[[j]]$type <- 1
      sos[[j]]$index <- c((n + j), (n + p + j))
      sos[[j]]$weights <- c(1, 2)
    }
    model$sos <- sos

  }

  ## Parameters for Gurobi optimization
  params$OutputFlag <- 0
  params$timeLimit <- 100
  params$NonConvex <- 2

  result <- gurobi(model, params)

  result$algo2Time <- algo2time
  result$algo2Sol <- algo2Sol

  return(result)
}

recoverDriftMatrix <- function(i, p, n, method) {

set.seed(i)

# True signal and Volatility matrix
M_star <- stabSignal(p)
trueSparsity <- (M_star != 0)
C <- diag(runif(p))
# Population covariance
Sigma_star <- getPopCovar(M_star, C)
# Sampled data and covariance
if (n == Inf) {
  Sigma <- Sigma_star
}else {
  X <- sampleX(n, Sigma_star)
  Sigma <- getDataCovar(X)
}
# Commutaion matrix that transforms column vetorization of
# a matrix to its row vectorization
KPermute <- spMatrix(nrow = p**2, ncol = p**2, i = c(seq(1, p**2)),
                    j = c(sapply(seq(1, p),
                          function(x) {seq(x, p**2, by = p)})),
                    x = c(rep(1, p**2)))
# Data matrix for regression
ASigma <- as.matrix(kronecker(Sigma, diag(p)) +
                    kronecker(diag(p), Sigma) %*% KPermute)
vecC <- c(C)

# True number of best subset
k0 <- sum(M_star != 0)

# Storing the best result after hyperparameter search
bestResult <- list()

if (method != "tLasso") {

  if (method == "bs") {
    # Number of best subset to select (min. val set to p -
    # Atleast diagonal should be present)
    # k = seq(k0-k0%/%2,k0 + k0%/%2)
    k <- floor(seq(p, p**2, length.out = 10))


    # resultBS <- lapply(k, applyBS, ASigma, vecC)
    foreach(i = k) %dopar% {applyBS(i, ASigma, vecC)} -> resultBS

    betas <- sapply(resultBS, function(x) {return(x$beta)})

  }else if (method == "lasso") {
    # Lasso
    penalty <- rep(1, p**2)
    penalty[seq(1, p**2, by = p + 1)] <- 0

    initGrid <- seq(10, 10^-5, length.out = 100)
    coarseLasso <- glmnet(ASigma, -vecC, intercept = FALSE, alpha = 1,
                            standardize = TRUE, penalty.factor = penalty,
                            lambda = initGrid)
    # Find min lambda such that M is diagonal
    lambda <- min(initGrid[(colSums(penalty * (coarseLasso$beta!=0)) == 0 &
                              colSums((1 - penalty) *
                              (coarseLasso$beta!=0)) == p)])
    fineGrid <- seq(lambda, lambda / 10^4, length.out = 100)
    fineLasso <- glmnet(ASigma, -vecC, intercept = FALSE, alpha = 1,
                        standardize = TRUE, penalty.factor = penalty,
                        lambda = fineGrid)

    betas <- fineLasso$beta
  }

  # Choose the best drift matrix based on posterior model probabilities
  Siginv_comp <- apply(betas, 2, invCov, vecC)
  bic_scores <- apply(Siginv_comp, 2, bic_score, Sigma, n)
  posterior_prob <- postprb(bic_scores)
  bestM <- betas[, which.max(posterior_prob)]

  # Save the best result
  bestk <- sum(bestM!=0)
  bestResult$M_star <- M_star
  bestResult$M_comp <- bestM
  bestResult$bestk <- bestk
  bestResult$k <- k0

  # Evaluate Metrics

  bestResult$acc <- acc(bestM, trueSparsity)
  bestResult$tpr <- tpr(bestM, trueSparsity)
  bestResult$fpr <- fpr(bestM, trueSparsity)
  bestResult$f1score <- f1score(bestM, trueSparsity)
  bestResult$aucroc <- AUCROC(bestM, trueSparsity)
  bestResult$aucpr <- AUCPR(bestM, trueSparsity)

  return(bestResult)

}else {
  # Save Dataset to be later loaded in Matlab
writeMat(con = paste0("D:/TUM/Master-Thesis/Dataset/", n, "_", p, "/X_",
                      n, "_", p, "_", i, ".mat"), X = ASigma)
writeMat(con = paste0("D:/TUM/Master-Thesis/Dataset/", n, "_", p, "/y_",
                      n, "_", p, "_", i, ".mat"), y = -vecC)
writeMat(con = paste0("D:/TUM/Master-Thesis/Dataset/", n, "_", p, "/M_",
                      n, "_", p, "_", i, ".mat"), M = M_star)
}
}

tp <- function(comp,gt){return (length(intersect(which(comp!=0),which(gt!=0))))}
fp <- function(comp,gt){return (length(setdiff(which(comp!=0),which(gt!=0))))}
tn <- function(comp,gt){return (length(intersect(which(comp==0),which(gt==0))))}
fn <- function(comp,gt){return (length(setdiff(which(comp==0),which(gt==0))))}
acc <- function(comp,gt){
  accuracy <- (tp(comp,gt) + tn(comp,gt)) /(tp(comp,gt) + tn(comp,gt) 
                                            + fp(comp,gt) + fn(comp,gt))
  return(accuracy)
  }
tpr <- function(comp,gt){
  return(tp(comp,gt)/(tp(comp,gt) + fn(comp,gt)))
}
fpr <- function(comp,gt){
  return(fp(comp,gt)/(fp(comp,gt) + tn(comp,gt)))
}
f1score <- function(comp,gt){
  return(tp(comp,gt)/(tp(comp,gt) + 0.5*(fp(comp,gt) + fn(comp,gt))))
}

sigmoid <- function(x){
  return(exp(abs(x)) / (1 + exp(abs(x))))
}
precision <- function(comp,gt){
  return(tp(comp,gt) / (tp(comp,gt) + fp(comp,gt)))
}


AUCPR <- function(comp,gt){

  comp <- sigmoid(abs(comp))

  x = list()
  comp_sort <- c(0,sort(comp))
  for (i in 1:length(comp_sort)) {
    comp_i <- (comp > comp_sort[i])
    x$precision[i] <- precision(comp_i, gt)
    if (is.nan(x$precision[i])) x$precision[i] <- 1
    x$recall[i] <- tpr(comp_i,gt)
    if (is.nan(x$recall[i])) x$recall[i] <- 1
  }

  temp <- 0
  n <- length(x$precision)
  for (i in 1:(n-1)){
    ## trapezoidal quadrature rule
    temp <- temp + (x$precision[i + 1] + x$precision[i]) *
      (x$recall[i] - x$recall[i + 1]) / 2
  }


  return(temp)
}

AUCROC <- function(comp,gt) {

  comp <- sigmoid(abs(comp))
  x = list()
  comp_sort <- c(0,sort(comp))
  for (i in 1:length(comp_sort)) {
    comp_i <- (comp > comp_sort[i])
    x$tpr[i] <- tpr(comp_i,gt)
    if (is.nan(x$tpr[i])) x$tpr[i] <- 1
    x$fpr[i] <- fpr(comp_i,gt)
    if (is.nan(x$fpr[i])) x$fpr[i] <- 1
  }

  temp <- 0
  n <- length(x$tpr)
  for (i in 1:(n-1)){
    ## trapezoidal quadrature rule
    temp <- temp + (x$tpr[i + 1] + x$tpr[i]) *
      (x$fpr[i] - x$fpr[i + 1]) / 2
  }
  return(temp)
}

# Compute BIC Scores 
bic_score <- function(K, S, n){
  if (is.matrix(K)){
    # do nothing
  }else{
    K <- matrix(K, nrow = sqrt(length(K)))
  }
  if (n == Inf) n <- 10^10
  bic <- n*(-determinant(K)$modulus + sum(diag(S%*%K))) + sum(unique(K) != 0) * log(n)

  return(bic)
}

# Compute extended BIC scores
ebic_score <- function(K, S, n, gamma){
  if (is.matrix(K)){
    # do nothing
  }else{
    K <- matrix(K,nrow = sqrt(length(K)))
  }
  if (n==Inf) n <- 10^10
  p <- dim(K)[1]
  ebic = n*(-determinant(K)$modulus + sum(diag(S%*%K))) + sum(unique(K)!=0)*log(n) + 4*gamma*sum(unique(K)!=0)*log(p)

  return(ebic)
}

aic_score <- function(K, S, n) {
  if (is.matrix(K)){
    # do nothing
  }else{
    K <- matrix(K,nrow = sqrt(length(K)))
  }
  if (n==Inf) n <- 10^10
  aic <- n * (-determinant(K)$modulus + sum(diag(S %*% K))) + 2 * sum(unique(K)!=0)

  return(aic)
}

# Posterior model probability
postprb <- function(scores) {
  weights <- -scores / 2
  post_prb <- sapply(weights, function(x) {x / sum(weights)})
}