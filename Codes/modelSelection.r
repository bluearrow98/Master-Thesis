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

lyap_bic <- function(M, C_true, n) {
  if (is.matrix(M) && is.matrix(C_true)){
    # do nothing
  }else{
    M <- matrix(M, nrow = sqrt(length(M)))
    C_true <- matrix(C_true, nrow = sqrt(length(C_true)))
  }
  if (n == Inf) n <- 10^10
  Sigma <- getPopCovar(M, C_true)
  bic <- 2 * (log(norm(M %*% Sigma + Sigma %*% t(M) + C_true, type = "F") + 1e-10)) + sum(M != 0) * log(n)

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
  ebic = n*(-determinant(K)$modulus + sum(diag(S%*%K))) + sum(unique(K) != 0) * log(n)
          + 4 * gamma * sum(unique(K) != 0) * log(p)

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
  weights <- exp(-(scores - min(scores) + 1e-10) / 2)
  post_prb <- sapply(weights, function(x) {x / sum(weights)})
}
