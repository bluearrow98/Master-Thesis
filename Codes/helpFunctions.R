library("Rlab")
library(MASS)
library(Matrix)
library(parallel)
library(doParallel)
library(gurobi)

stabSignal <- function(p){
  k = 0.05 * p
  M = matrix(rnorm(p**2) * rbern(p ** 2, prob = k / p), nrow = p)
  nonDiagAbsSum = sum(abs(M)) - sum(diag(abs(M)))
  diag(M) = - nonDiagAbsSum - abs(rnorm(p))
  return(M)
}


getPopCovar <- function(M, C){

  p <- dim(M)[1]
  AM <- kronecker(diag(p), M) + kronecker(M, diag(p))
  E <- eigen(AM)
  Sigma <- solve(AM, -c(C))
  Sigma <- matrix(Sigma, nrow = p)

  return(Sigma)
}

getSignals <- function(i, p) {

  set.seed(i)

  # True signal and Volatility matrix
  M_star <- stabSignal(p)
  trueSparsity <- (M_star != 0)
  C <- diag(runif(p, min = 0.2))

  # True number of subsets
  k0 <- sum(trueSparsity)

  # Population covariance
  Sigma_star <- getPopCovar(M_star, C)

  return(list(M = M_star, C = C, popSigma = Sigma_star, gt = trueSparsity, k0 = k0))
}
 
invCov <- function(Mvec,Cvec){
  # Get number of variables
  p <- sqrt(length(Mvec))

  # Get the drift matrix and volatility matrix
  M <- matrix(Mvec, nrow = p)
  C <- matrix(Cvec, nrow = p)

  # Get the recovered covariance matrix
  Sigma <- getPopCovar(M, C)
  
  # Find the inverse using decomposition
  E <- eigen(Sigma)
  Lambda <- diag(E$values^(-1))
  Gamma <- E$vectors
  invCov <- Gamma %*% Lambda %*% t(Gamma)

  return(invCov)
}

sampleX <- function(n,Sigma){
  X = mvrnorm(n, mu = rep(0, dim(Sigma)[1]), Sigma)
  return(X)
}

getDataCovar <- function(X){
  n <- dim(X)[1]
  Sigma <- t(X) %*% X / n
  
  return(Sigma)
}

# standardize <- function(X){
#   n = dim(X)[1]
#   p = dim(X)[2]
#   mean = matrix(rep(t(apply(X, 2, function(x) {mean(x)})),n) , ncol = p, byrow = T)
#   X_shift = X - mean
#   X_norm = X_shift%*%diag(1/apply(X_shift, 2, function(x) {norm(x,type = "2")}))
#   return (X_shift)
# }

applyBS <- function(i, j, signals) {

  ASigma <- signals[[i]]$ASigma
  vecC <- signals[[i]]$vecC

  n <- dim(ASigma)[1]
  p <- length(vecC)
  result <- customBS(ASigma, -vecC, j, polish = FALSE)

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

customBS <- function(X, y, k, polish = TRUE) {

  n <- dim(X)[1]
  p <- dim(X)[2]

  # Discrete First order method
  algo2time <- system.time(algo2Sol <- proj_grad_oneEdge(X, y, k, polish = polish))
  # algo2time <- system.time(algo2Sol <- applyLasso(X, y))
  M_p <- algo2Sol$bm

  ## MIO Formulation
  eps <- 1e-3
  d_ind <- seq(1, p, by = sqrt(p) + 1)
  Mu_d <- 2 * norm(as.matrix(M_p[d_ind]), type = "i")
  if (min(abs(M_p[d_ind])) < 1e-8) Mu_d_upper <- eps else Mu_d_upper <- min(abs(M_p[d_ind])) - eps
  Mu_nd <- 2 * norm(as.matrix(M_p[-d_ind]), type = "i")
  if (Mu_nd == 0) Mu_nd <- Mu_d

  ## Initial guess (Not used for cold start)
  z0 <- rep(0,p)
  z0[M_p != 0] <- 1

  diagA <- rep(-Mu_nd, p)
  diagA[d_ind] <- -Mu_d

  model <- list()
  params <- list()
  if (p < n) {
    # Warm-start model creation (Eq 2.5)
    model$A <- spMatrix(nrow = 2 + 2 * p, ncol = 2 * p, i = c(rep(1, p), rep(2, sqrt(p)), rep(seq(3, 2 + p), 2), rep(seq(2 + p + 1, 2 + 2 * p), 2)),
                       j = c(seq(p + 1, 2 * p), seq(p + 1, 2 * p, by = sqrt(p) + 1), rep(seq(1, 2 * p), 2)),
                       x = c(rep(1, p + sqrt(p)), rep(1, p), diagA, rep(-1, p), diagA))
    model$Q <- spMatrix(nrow = 2 * p,ncol = 2 * p,
                        i = c(rep(seq(1, p), each = p)),
                        j = c(rep(seq(1, p), p)), x = c(0.5 * t(X) %*% X))
    model$obj <- c(-t(X) %*% y, rep(0, p))
    model$ub <- c(rep(-Mu_d_upper, p), rep(1, p))
    model$ub[-c(d_ind, seq(p + 1, 2 * p))] <- Mu_nd
    model$lb <- c(rep(-Mu_d + 1e-6, p), rep(0, p))
    model$lb[-c(d_ind, seq(p + 1, 2 * p))] <- -Mu_nd
    model$rhs <- c(k, sqrt(p), rep(0, 2 * p))
    model$sense <- c("<=", "=", rep("<=", 2 * p))
    model$vtype <- c(rep("C", p), rep("B", p))
    model$start <- c(M_p, z0)

    ## Creating SOS constraints
    # sos <- list()
    # for (j in 1:p){
    #   sos[[j]] <- list()
    #   sos[[j]]$type <- 1
    #   sos[[j]]$index <- c(j, (p + j))
    #   sos[[j]]$weights <- c(1, 2)
    # }
    # model$sos <- sos

  }else{
    # Warm-start model creation (Eq 2.6)
    M_zeta <- max(colSums(apply(abs(X), 1, sort,
                    decreasing=TRUE)[, 1:k, drop = FALSE])) * Mu_d
    model$A <- spMatrix(nrow = n + 2 + 2 * p, ncol = n + 2 * p,
                       i = c(seq(1, n), rep(seq(1, n), each = p), rep(n + 1, p), rep(n + 2, sqrt(p)), rep(seq(n + 3, n + 2 + p), 2), rep(seq(n + 3 + p, n + 2 + 2 * p),2)),
                       j = c(seq(1, n), rep(seq(n+1, n+p), n), seq(n + p + 1, n + 2 * p), seq(n + p + 1, n + 2 * p, by = sqrt(p) + 1), rep(seq(n + 1, n + 2 * p), 2)),
                       x = c(rep(1, n), c(-t(X)), rep(1, p), rep(1, sqrt(p)), rep(1, p), diagA, rep(-1, p), diagA))
    model$Q <- spMatrix(nrow = n + 2 * p, ncol = n + 2 * p,
                       i = seq(1, n), j = seq(1, n), x = rep(0.5, n))
    model$obj <- c(rep(0, n), -t(X) %*% y, rep(0, p))
    model$ub <- c(rep(M_zeta, n), rep(-Mu_d_upper, p), rep(1, p))
    model$ub[-c(seq(1, n), n + d_ind, seq(n + p + 1, n + 2 * p))] <-
      c(rep(Mu_nd, p - sqrt(p)))
    model$lb <- c(rep(-M_zeta, n), rep(-Mu_d, p), rep(0, p))
    model$lb[-c(seq(1, n), n + d_ind, seq(n + p + 1, n + 2 * p))] <-
      c(rep(-Mu_nd, p - sqrt(p)))
    model$rhs <- c(rep(0, n), k, sqrt(p), rep(0, 2 * p))
    model$sense <- c(rep("=", n), "<=", "=", rep("<=", 2 * p))
    model$vtype <- c(rep("C", n), rep("C", p), rep("B", p))
    model$start <- c(X %*% M_p, M_p, z0)

    ## Creating SOS constraints
    # sos <- list()
    # for (j in 1:p){
    #   sos[[j]] <- list()
    #   sos[[j]]$type <- 1
    #   sos[[j]]$index <- c((n + j), (n + p + j))
    #   sos[[j]]$weights <- c(1, 2)
    # }
    # model$sos <- sos

  }

  ## Parameters for Gurobi optimization
  params$OutputFlag <- 1
  params$timeLimit <- 100
  params$NonConvex <- 2
  params$MIPGap <- 0.01

  result <- gurobi(model, params)
  
  result$algo2Time <- algo2time
  result$algo2Sol <- algo2Sol

  return(result)
}

recoverDriftMatrix <- function(signals, p, method,
                                k = lapply(1: length(signals), function(x){
                                  floor(seq(p, p ** 2, length.out = 10))})) {

  if (method == "bs") {
    # resultBS <- lapply(k, applyBS, ASigma, vecC)
    resultBS <- foreach(i = 1:length(signals)) %:%
                  foreach(j = k[[i]]) %dopar% {
                    applyBS(i, j, signals)
                  }

    # Check if all solutions obtained are feasible
    feasibleResults <- lapply(resultBS, checkMIPGap, k)

    resultBS <- lapply(feasibleResults, function(z) {
          return(z$results)
    })
    # Get solutions
    betas <- lapply(resultBS, sapply, function(x) {
          return(x$beta)
          })
    # betas_fix <- lapply(betas, function(x){
    #   apply(x, 2, fixInvertibility)
    # })

    objvals <- lapply(resultBS, sapply, function(x) {
          return(x$objval)
          })


  }else if (method == "lasso") {
      # Lasso
      penalty <- rep(1, p ** 2)
      penalty[seq(1, p ** 2, by = p + 1)] <- 0

      initGrid <- seq(10, 10^-5, length.out = 100)
      coarseLasso <- foreach(i = 1:length(signals)) %dopar% {
                        glmnet(signals[[i]]$ASigma, -signals[[i]]$vecC, 
                          intercept = FALSE, alpha = 1,
                          standardize = TRUE, penalty.factor = penalty,
                          lambda = initGrid)
                      }

      # Find min lambda such that M is diagonal
      lambdas <- lapply(coarseLasso, function(z) {
            return(min(initGrid[(colSums(penalty * (z$beta != 0)) == 0 &
          colSums((1 - penalty) * (z$beta != 0)) == p)]))
            })
      fineGrid <- lapply(lambdas, function(z) {
            return(seq(z, z / 10^4, length.out = 100))
            })
      fineLasso <- foreach(i = 1:length(signals)) %dopar% {
                        glmnet(signals[[i]]$ASigma, -signals[[i]]$vecC, 
                          intercept = FALSE, alpha = 1,
                          standardize = TRUE, penalty.factor = penalty,
                          lambda = fineGrid[[i]])
                      }
      # Get solution                    }
      betas <- lapply(fineLasso, function(z) {
            return(z$beta)
            })

      # Get objective values, if necessary
      objvals <- lapply(1:length(betas), function(i) {
            return(apply(betas[[i]], 2, function(z){
              norm(signals[[i]]$ASigma %*% z + signals[[i]]$vecC, type = "2")}))
            })

    }

    return(betas)
}

problem <- function(signal, p, n){

  # Get the population covariance matrix and the Volatility matrix
  C <- signal$C
  Sigma_star <- signal$popSigma


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

  # Vectorize volatility matrix for regression
  vecC <- c(C)

  signal$ASigma <- ASigma
  signal$Sigma <- Sigma
  signal$p <- p
  signal$n <- n
  signal$X <- X
  signal$vecC <- vecC

  return(signal)

}

setupProblems <- function(signals, p, n) {

  signals <- lapply(signals, problem, p, n)

  return(signals)
}

# Check for small diagonal entries matrices
is.diagonal <- function(Mvec, p) {

  diag_ind <- seq(1, p ** 2, by = p + 1)
  # M <- matrix(Mvec, nrow = p)
  # AM <- kronecker(diag(p), M) + kronecker(M, diag(p))
  # lambda <- eigen(AM)$values
  return(sum(abs(Mvec[diag_ind]) >= 1e-3) == p)
}

# Accept results only with an MIP gap of less than 1 % 
checkMIPGap <- function(results, k){
  #  Get the indices at which we get a feasible solution
  fInd <- which(sapply(results,
                function(x){if(!is.null(x$mipgap)) {
                  return(x$mipgap <= 0.05)
                }else{
                  return(FALSE)
                }}))

  # Keep only the feasible results
  resultsFeasible <- results[fInd]

  feasibleK <- k[fInd]

  return(list(fInd = fInd, results = resultsFeasible, k = feasibleK))
}

# Fix singular matrices using nearPD (not used)
fixInvertibility <- function(beta){

  p <- sqrt(length(beta))
  M <- matrix(beta, nrow = p)
  E <- eigen(M)
  if(all(abs(E$values)) < 1e-3){
    print(TRUE)
    if(all(E$values) >= 0){
      M <- nearPD(M, base.matrix = TRUE, conv.tol = 1e-4)$mat
    }else{
      M <- -nearPD(-M, base.matrix = TRUE, conv.tol = 1e-4)$mat
    }
  }

  return(c(M))
}

getPreSelectionMetrics <- function(signals, betas, p, n, methodDir) {

  results <- list()

  results$betas <- betas

  # Compute metrics for only off-diagonal elements
  diag <- seq(1, p ** 2, by = p + 1)
  # Evaluate Metrics
  results$max_acc <- mean(sapply(1:length(betas), function(x){
    max(apply(betas[[x]], 2, function(y){
      acc(y[-diag], signals[[x]]$gt[-diag])
    }))
  }))
  results$max_f1 <- mean(sapply(1:length(betas), function(x){
    max(apply(betas[[x]], 2, function(y){
      score <- f1score(y[-diag], signals[[x]]$gt[-diag])
      if (is.nan(score)) score <- 1
      return(score)
    }))
  }))
  results$tpr <- sapply(1:length(betas), function(x){
    apply(betas[[x]], 2, function(y){
      score <- tpr(y[-diag], signals[[x]]$gt[-diag])
      if (is.nan(score)) score <- 1
      return(score)
    })
  })
  results$fpr <- sapply(1:length(betas), function(x){
    apply(betas[[x]], 2, function(y){
      score <- fpr(y[-diag], signals[[x]]$gt[-diag])
      if (is.nan(score)) score <- 1
      return(score)
    })
  })
  results$precision <- sapply(1:length(betas), function(x){
    apply(betas[[x]], 2, function(y){
      score <- precision(y[-diag], signals[[x]]$gt[-diag])
      if (is.nan(score)) score <- 1
      return(score)
    })
  })

  # Area under the Curve along the regularization path
  results$aucroc <- mean(sapply(1:length(betas), function(x){
    AUCROC_path(results$tpr[, x], results$fpr[, x])
  }))
  results$aucpr <- mean(sapply(1:length(betas), function(x){
    AUCPR_path(results$precision[, x], results$tpr[, x])
  }))

  # Save the results
    save(results,
      file = paste0(getwd(), "/../Results/", methodDir, "/allModels_edgeProb5/", n, "_", p, "_res.RData"))
}

saveResults <- function(signals, results, p, n, elapsedTime, methodDir) {

      betas <- results

      # Evaluate Pre model selection metrics
      getPreSelectionMetrics(signals, betas, p, n, methodDir)

      ## Model Selection

      # Choose the best drift matrix based on posterior model probabilities
      Siginv_comp <- lapply(1:length(betas), function(x) {
        return(apply(betas[[x]], 2, invCov, signals[[x]]$vecC))
        })
      # bic_scores <- lapply(1:length(betas), function(x) {
      # return(apply(betas[[x]], 2, lyap_bic, signals[[x]]$vecC, n))
      # })
      bic_scores <- lapply(1:length(betas), function(x) {
        return(apply(Siginv_comp[[x]], 2, bic_score, signals[[x]]$Sigma, signals[[x]]$n))
        }) 
      posterior_prob <- lapply(bic_scores, postprb)
      bestM <- lapply(1:length(betas), function(x) {
        return(b <- betas[[x]][, which.max(posterior_prob[[x]])])
      })

      # Choose the best drift matrix based on objective values
      # bestM <- betas[, which.min(objvals)]

      bestk <- lapply(bestM, function(x) {return(sum(x != 0))})
      bestResult <- lapply(1:length(signals), function(x) {
        return(list(M_star = signals[[x]]$M, M_comp = bestM[[x]],
        bestk = bestk[[x]], k = signals[[x]]$k,
        t = elapsedTime))
      })

      # Save the results
      save(bestResult,
      file = paste0(getwd(), "/../Results/test/", n, "_", p, "_res.RData"))

      # Return, if required for further calculation of metrics
      return(bestResult)
}

getMetrics <- function(signals, bestResult, p, n, method, methodDir){

    # results <- recoverDriftMatrix(signals, p, method, k = lapply(bestResult, function(x){return(x$bestk)}))

    # Evaluate Metrics
    diag <- seq(1, p ** 2, by = p + 1)
    bestResult <- apply(as.matrix(1:length(signals)), 1,  function(x) {
      # Remove noise
      bestResult[[x]]$M_noiseFree <- bestResult[[x]]$M_comp
      # bestResult[[x]]$M_noiseFree[abs(bestResult[[x]]$M_noiseFree) < 1e-6] <- 0

      bestResult[[x]]$acc <- acc(bestResult[[x]]$M_noiseFree[-diag],
        signals[[x]]$gt[-diag])
      bestResult[[x]]$tpr <- tpr(bestResult[[x]]$M_noiseFree[-diag],
        signals[[x]]$gt[-diag])
      if (is.nan(bestResult[[x]]$tpr)) bestResult[[x]]$tpr <- 1
      bestResult[[x]]$fpr <- fpr(bestResult[[x]]$M_noiseFree[-diag],
        signals[[x]]$gt[-diag])
      if (is.nan(bestResult[[x]]$fpr)) bestResult[[x]]$fpr <- 1
      bestResult[[x]]$f1score <- f1score(bestResult[[x]]$M_noiseFree[-diag],
        signals[[x]]$gt[-diag])
      if (is.nan(bestResult[[x]]$f1score)) bestResult[[x]]$f1score <- 1
      bestResult[[x]]$aucroc <- AUCROC(bestResult[[x]]$M_noiseFree[-diag],
        signals[[x]]$gt[-diag])
      bestResult[[x]]$aucpr <- AUCPR(bestResult[[x]]$M_noiseFree[-diag],
        signals[[x]]$gt[-diag])
      return(bestResult[[x]])
    })

    res <- bestResult

    # Get average over all signals
    bestResult$accMax <- max(unlist(lapply(res, function(x) {x$acc})))
    bestResult$tprAvg <- mean(unlist(lapply(res, function(x) {x$tpr})))
    bestResult$fprAvg <- mean(unlist(lapply(res, function(x) {x$fpr})))
    bestResult$f1scoreAvg <- mean(unlist(lapply(res, function(x) {x$f1score})))
    bestResult$aucrocAvg <- mean(unlist(lapply(res, function(x) {x$aucroc})))
    bestResult$aucprAvg <- mean(unlist(lapply(res, function(x) {x$aucpr})))
    bestResult$timeAvg <- res[[1]]$t / length(signals)

    # Save the results
    save(bestResult,
      file = paste0(getwd(), "/../Results/", methodDir, "/completeResultsWithDiagMetric_edgeProb5/", n, "_", p, "_res.RData"))

}